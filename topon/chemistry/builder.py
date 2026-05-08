"""
Chemistry Builder for Topon.

Builds molecular structures based on assigned graph attributes.
Supports atomistic (DREIDING) and coarse-grained (Kremer-Grest) models.

Key features:
- Smart connection chemistry (auto-bridge atoms)
- Node type → molecule mapping
- Edge type → monomer mapping
- Copolymer sequences
- Grafted chains
"""

import warnings
from typing import Optional
import networkx as nx
import numpy as np

from topon.config.schema import ChemistryConfig


# Feature-detect marker for downstream projects.
# True since the fix that removed the trailing "[O]" placeholder from
# `_create_chain_from_smiles`, which produced a peroxide -O-O- bond at the
# chain tail for O-terminal monomers (e.g. PDMS, PTFPMS) under the atomistic
# auto-bridge path. See docs/topon_issues_from_solubility_v2.md §P0-1.
_PEROXIDE_FIX_APPLIED = True


class ChemistryBuilder:
    """
    Builds molecular structure from attributed graph.
    
    Uses configuration to map:
    - Node types → crosslinker molecules (Si, POSS, custom)
    - Edge types → chain monomers (PDMS, FPDMS, etc.)
    - Auto-detects when bridge atoms are needed
    """
    
    def __init__(
        self, 
        G: nx.MultiGraph, 
        dims: Optional[np.ndarray],
        config: ChemistryConfig
    ):
        """
        Initialize the chemistry builder.
        
        Args:
            G: Attributed graph from assignment module.
            dims: Box dimensions.
            config: Chemistry configuration.
        """
        self.G = G
        self.dims = dims
        self.config = config
        
        # RDKit molecule (lazy import to avoid hard dependency)
        self.chemical_space = None
        
        # Atom mappings
        self.node_map = {}  # node_id -> atom_idx or list of atom_idxs (for POSS)
        self.edge_atom_map = {}  # edge_id -> list of atom_idxs
        self.edge_backbone_map = {}  # edge_id -> (head_idx, tail_idx)
        self.graft_atom_map = {}  # edge_id -> {backbone_pos: [graft_atom_idxs]}

        # For POSS corner tracking
        self.poss_usage = {}  # node_id -> set of used corner indices

        # Entangled pairs
        self.entangled_pairs = []
    
    def build(self):
        """
        Build the molecular structure.
        
        Returns:
            RDKit RWMol object with complete molecular structure.
        """
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError("RDKit is required for chemistry building. Install with: pip install rdkit")
        
        print("Building chemistry...")
        print(f"  Model type: {self.config.model_type}")
        
        # Initialize RWMol
        self.chemical_space = Chem.RWMol()
        
        # Get entangled pairs from graph
        self._extract_entangled_pairs()
        
        # Build nodes (crosslinkers)
        self._build_nodes()
        
        # Build chains (edges)
        self._build_chains()
        
        # Report
        print(f"  Total atoms: {self.chemical_space.GetNumAtoms()}")
        print(f"  Total bonds: {self.chemical_space.GetNumBonds()}")
        
        return self.chemical_space
    
    def _extract_entangled_pairs(self):
        """Extract entangled edge pairs from graph attributes."""
        processed = set()
        for u, v, key, data in self.G.edges(keys=True, data=True):
            partner = data.get("entangled_with")
            if partner:
                edge = (u, v, key)
                pair = tuple(sorted([edge, partner]))
                if pair not in processed:
                    self.entangled_pairs.append(pair)
                    processed.add(pair)
        
        if self.entangled_pairs:
            print(f"  Entangled pairs: {len(self.entangled_pairs)}")
    
    def _build_nodes(self):
        """Build crosslinker structures for each node."""
        from rdkit import Chem
        
        print("  Building nodes...")
        
        for node in self.G.nodes():
            degree = self.G.degree(node)
            if degree == 0:
                continue

            # Get node type from graph.  For backward compatibility with
            # legacy callers (notably generate_matrix.py in solvent_effects),
            # accept the attribute name ``type`` as a fallback.  Emit a
            # DeprecationWarning so new code migrates to ``node_type``.
            node_attrs = self.G.nodes[node]
            if "node_type" in node_attrs:
                node_type = node_attrs["node_type"]
            elif "type" in node_attrs:
                import warnings
                warnings.warn(
                    "Graph node attribute 'type' is read for backward-compat; "
                    "please use 'node_type' in new code.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                node_type = node_attrs["type"]
            else:
                node_type = "A"

            # Get molecule config for this type
            node_config = self.config.node_type_map.get(node_type)
            if not node_config:
                # Default to Si
                molecule = "Si"
                is_end_cap = degree == 1
            else:
                molecule = node_config.molecule
                is_end_cap = node_config.is_end_cap
            
            # Build the appropriate structure
            if molecule.upper() == "POSS" or molecule.upper() == "SI8O12":
                self._place_poss_cage(node)
            elif molecule.upper() == "POSS_AM0270":
                self._place_poss_am0270(node)
            elif is_end_cap or degree == 1:
                self._place_end_cap(node, molecule)
            else:
                self._place_simple_atom(node, molecule)
        
        print(f"    Placed {len(self.node_map)} node structures")
    
    def _place_simple_atom(self, node: int, atom_symbol: str = "Si"):
        """Place a simple atom crosslinker."""
        from rdkit import Chem
        
        # Handle SMILES vs atom symbol
        if len(atom_symbol) <= 2 and atom_symbol.isalpha():
            atom = Chem.Atom(atom_symbol)
            idx = self.chemical_space.AddAtom(atom)
            self.node_map[node] = idx
        else:
            # It's a SMILES string
            mol = Chem.MolFromSmiles(atom_symbol)
            if mol:
                idxs = [self.chemical_space.AddAtom(a) for a in mol.GetAtoms()]
                for b in mol.GetBonds():
                    self.chemical_space.AddBond(
                        idxs[b.GetBeginAtomIdx()],
                        idxs[b.GetEndAtomIdx()],
                        b.GetBondType()
                    )
                self.node_map[node] = idxs[0]  # Use first atom as attachment point
            else:
                # Fallback to Si
                idx = self.chemical_space.AddAtom(Chem.Atom("Si"))
                self.node_map[node] = idx
    
    def _place_end_cap(self, node: int, molecule: str):
        """Place an end-cap molecule."""
        from rdkit import Chem
        
        mol = Chem.MolFromSmiles(molecule)
        if mol:
            mol = Chem.RemoveHs(mol)
            idxs = [self.chemical_space.AddAtom(a) for a in mol.GetAtoms()]
            for b in mol.GetBonds():
                self.chemical_space.AddBond(
                    idxs[b.GetBeginAtomIdx()],
                    idxs[b.GetEndAtomIdx()],
                    b.GetBondType()
                )
            # Find attachment point (usually Si)
            for i, a in enumerate(mol.GetAtoms()):
                if a.GetSymbol() == "Si":
                    self.node_map[node] = idxs[i]
                    return
            self.node_map[node] = idxs[0]
        else:
            self._place_simple_atom(node, "Si")
    
    def _place_poss_cage(self, node: int):
        """Place a POSS cage structure (Si8O12)."""
        from rdkit import Chem
        
        # Create 8 Si atoms at corners
        ids = [self.chemical_space.AddAtom(Chem.Atom("Si")) for _ in range(8)]
        
        # Connect with O bridges (cube edges)
        # Cube connectivity: edges of a cube
        conns = [
            (0, 1), (0, 2), (0, 4),
            (1, 3), (1, 5),
            (2, 3), (2, 6),
            (3, 7),
            (4, 5), (4, 6),
            (5, 7),
            (6, 7)
        ]
        
        for a, b in conns:
            o = self.chemical_space.AddAtom(Chem.Atom("O"))
            self.chemical_space.AddBond(ids[a], o, Chem.BondType.SINGLE)
            self.chemical_space.AddBond(o, ids[b], Chem.BondType.SINGLE)
        
        self.node_map[node] = ids  # List of 8 corner atoms
        self.poss_usage[node] = set()

    def _place_poss_am0270(self, node: int):
        """
        Place AM0270 POSS (AminopropylIsooctyl POSS).
        
        Structure:
        - Si8O12 core (cubic cage)
        - 7 corners (1-7) functionalized with IsoOctyl (2,4,4-trimethylpentyl)
        - 1 corner (0) functionalized with Propyl linker -> connects to network
        
        Corner layout (cube vertices):
            Corner 0: (-1, -1, -1)  <- Propyl linker (network connection)
            Corner 1: (-1, -1, +1)
            Corner 2: (-1, +1, -1)
            Corner 3: (-1, +1, +1)
            Corner 4: (+1, -1, -1)
            Corner 5: (+1, -1, +1)
            Corner 6: (+1, +1, -1)
            Corner 7: (+1, +1, +1)
        
        The IsoOctyl arms extend along the space diagonal (outward from cage center).
        """
        from rdkit import Chem
        
        # Initialize structure tracking if not exists
        if not hasattr(self, 'poss_structure'):
            self.poss_structure = {}
        
        # 1. Build Base Cage (Si8O12)
        corner_si_ids = [self.chemical_space.AddAtom(Chem.Atom("Si")) for _ in range(8)]
        
        # Oxygen bridges (12 edges of cube)
        cage_oxygen_ids = []
        conns = [
            (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6),
            (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)
        ]
        
        for a, b in conns:
            o = self.chemical_space.AddAtom(Chem.Atom("O"))
            cage_oxygen_ids.append(o)
            self.chemical_space.AddBond(corner_si_ids[a], o, Chem.BondType.SINGLE)
            self.chemical_space.AddBond(o, corner_si_ids[b], Chem.BondType.SINGLE)
        
        # 2. Prepare Functional Groups
        isooctyl_smiles = "CC(C)CC(C)(C)C"  # 2,4,4-trimethylpentyl (8 carbons)
        propyl_smiles = "CCC"  # Propyl linker
        
        # 3. Attach Propyl to Corner 0 (Network Connection)
        mol_prop = Chem.MolFromSmiles(propyl_smiles)
        mol_prop = Chem.RemoveHs(mol_prop)
        
        propyl_idxs = [self.chemical_space.AddAtom(a) for a in mol_prop.GetAtoms()]
        for b in mol_prop.GetBonds():
            self.chemical_space.AddBond(
                propyl_idxs[b.GetBeginAtomIdx()],
                propyl_idxs[b.GetEndAtomIdx()],
                b.GetBondType()
            )
        
        # Connect first C to Si (Corner 0)
        self.chemical_space.AddBond(corner_si_ids[0], propyl_idxs[0], Chem.BondType.SINGLE)
        
        # Register the LAST carbon as the network attachment point
        self.node_map[node] = propyl_idxs[-1]
        
        # 4. Attach IsoOctyl to Corners 1-7 (Dangling Ends)
        isooctyl_arms = {}  # corner_idx -> list of atom indices
        
        for corner_idx in range(1, 8):
            mol_iso = Chem.MolFromSmiles(isooctyl_smiles)
            mol_iso = Chem.RemoveHs(mol_iso)
            
            iso_idxs = [self.chemical_space.AddAtom(a) for a in mol_iso.GetAtoms()]
            
            for b in mol_iso.GetBonds():
                self.chemical_space.AddBond(
                    iso_idxs[b.GetBeginAtomIdx()],
                    iso_idxs[b.GetEndAtomIdx()],
                    b.GetBondType()
                )
            
            # Connect first C to Si (Corner)
            self.chemical_space.AddBond(corner_si_ids[corner_idx], iso_idxs[0], Chem.BondType.SINGLE)
            
            isooctyl_arms[corner_idx] = iso_idxs
        
        # 5. Store structure metadata for coordinate generation
        self.poss_structure[node] = {
            'corner_si_ids': corner_si_ids,          # [8 Si atom indices]
            'cage_oxygen_ids': cage_oxygen_ids,      # [12 O atom indices]
            'propyl_arm': {
                'corner_idx': 0,
                'atom_ids': propyl_idxs              # [3 C atom indices]
            },
            'isooctyl_arms': isooctyl_arms           # {corner_idx: [8 C atom indices]}
        }

    
    def _get_attachment_atom(self, node: int, vec: np.ndarray) -> int:
        """
        Get the attachment atom for a node.
        
        For POSS, selects corner based on direction vector.
        For simple atoms, returns the atom index directly.
        """
        target = self.node_map.get(node)
        
        if target is None:
            return None
        
        if isinstance(target, int):
            return target
        
        # POSS cage - select corner based on direction
        # Corner vectors (normalized cube corners)
        corner_vecs = np.array([
            [-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1],
            [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1]
        ], dtype=float)
        corner_vecs = corner_vecs / np.linalg.norm(corner_vecs, axis=1, keepdims=True)
        
        # Normalize direction
        v_dir = vec / (np.linalg.norm(vec) + 1e-9)
        
        # Find best matching corner not yet used
        used = self.poss_usage.get(node, set())
        dots = corner_vecs @ v_dir
        ranked = np.argsort(dots)[::-1]  # Best match first
        
        for idx in ranked:
            if idx not in used:
                used.add(idx)
                self.poss_usage[node] = used
                return target[idx]
        
        # All used, return first one again
        return target[0]
    
    def _build_chains(self):
        """Build chain structures for each edge."""
        print("  Building chains...")
        
        chain_count = 0
        for u, v, key, data in self.G.edges(keys=True, data=True):
            dp = data.get("dp", 25)
            edge_type = data.get("edge_type", "A")

            # Get monomer for this edge type
            edge_config = self.config.edge_type_map.get(edge_type)
            if edge_config:
                monomer_name = edge_config.monomer
            else:
                monomer_name = "PDMS"  # Default

            # Get monomer definition
            monomer_config = self.config.monomers.get(monomer_name)
            if not monomer_config:
                print(f"    Warning: Monomer {monomer_name} not found, using PDMS")
                monomer_config = self.config.monomers.get("PDMS")

            # Build the chain (pass full edge data for graft/copolymer attributes)
            self._build_chain(u, v, key, dp, monomer_config, data)
            chain_count += 1
        
        print(f"    Built {chain_count} chains")
    
    def _build_chain(self, u: int, v: int, key: int, dp: int, monomer_config, edge_data: dict = None):
        """Build a single chain between two nodes."""
        from rdkit import Chem
        
        # Get positions for direction calculation
        pos_u = self.G.nodes[u].get("pos")
        pos_v = self.G.nodes[v].get("pos")
        
        if pos_u is not None and pos_v is not None:
            pos_u = np.array(pos_u)
            pos_v = np.array(pos_v)
            
            # MIC vector
            if self.dims is not None:
                raw = pos_v - pos_u
                vec = raw - self.dims * np.round(raw / self.dims)
            else:
                vec = pos_v - pos_u
        else:
            vec = np.array([1, 0, 0])  # Default direction
        
        # Get attachment atoms
        att_u = self._get_attachment_atom(u, vec)
        att_v = self._get_attachment_atom(v, -vec)
        
        if att_u is None or att_v is None:
            return  # Skip if nodes not built
        
        # For coarse-grained, just add beads
        if self.config.model_type == "coarse_grained":
            self._build_chain_cg(u, v, key, dp, att_u, att_v, edge_data or {})
        else:
            self._build_chain_atomistic(u, v, key, dp, monomer_config, att_u, att_v)
    
    def _build_chain_cg(self, u, v, key, dp, att_u, att_v, edge_data: dict = None):
        """Build a coarse-grained chain (simple bead chain).

        Reads edge attributes (set by AssignmentManager) to:
        - Set per-bead ``bead_type`` from ``monomer_sequence`` if present.
        - Attach graft side chains at positions listed in ``graft_positions``.
        """
        from rdkit import Chem

        if edge_data is None:
            edge_data = {}

        edge_id = (u, v, key)
        chain_idxs = []

        # Per-bead monomer sequence (length == dp when set by copolymer assignment)
        monomer_seq = edge_data.get("monomer_sequence")  # list[str] | None

        # Create backbone beads
        for i in range(dp):
            bead = Chem.Atom("C")  # Use C as generic CG bead element
            bead_type = (monomer_seq[i] if monomer_seq and i < len(monomer_seq) else "B")
            bead.SetProp("bead_type", bead_type)
            idx = self.chemical_space.AddAtom(bead)
            chain_idxs.append(idx)

        # Connect backbone beads
        for i in range(len(chain_idxs) - 1):
            self.chemical_space.AddBond(chain_idxs[i], chain_idxs[i + 1], Chem.BondType.SINGLE)

        # Connect backbone to nodes
        if chain_idxs:
            self.chemical_space.AddBond(att_u, chain_idxs[0], Chem.BondType.SINGLE)
            self.chemical_space.AddBond(chain_idxs[-1], att_v, Chem.BondType.SINGLE)

        # Graft side chains
        graft_positions = edge_data.get("graft_positions", [])
        graft_dp = edge_data.get("graft_dp", 5)
        graft_monomer = edge_data.get("graft_monomer", "G")

        graft_map = {}  # backbone_pos -> [side-chain atom indices]
        for pos in graft_positions:
            if pos < 0 or pos >= len(chain_idxs):
                continue
            backbone_idx = chain_idxs[pos]
            side_idxs = []
            prev = backbone_idx
            for _ in range(graft_dp):
                g_bead = Chem.Atom("C")
                g_bead.SetProp("bead_type", graft_monomer)
                g_idx = self.chemical_space.AddAtom(g_bead)
                self.chemical_space.AddBond(prev, g_idx, Chem.BondType.SINGLE)
                side_idxs.append(g_idx)
                prev = g_idx
            graft_map[pos] = side_idxs

        self.edge_atom_map[edge_id] = chain_idxs
        self.edge_backbone_map[edge_id] = (chain_idxs[0], chain_idxs[-1]) if chain_idxs else (None, None)
        if graft_map:
            self.graft_atom_map[edge_id] = graft_map
    
    def _build_chain_atomistic(self, u, v, key, dp, monomer_config, att_u, att_v):
        """Build an atomistic chain with smart bridge detection."""
        from rdkit import Chem
        
        edge_id = (u, v, key)
        smiles = monomer_config.smiles
        chain_head_atom = monomer_config.chain_head
        chain_tail_atom = monomer_config.chain_tail
        
        # Create chain from SMILES
        try:
            chain_mol = self._create_chain_from_smiles(smiles, dp)
        except ValueError as exc:
            warnings.warn(
                f"Skipping edge ({u}, {v}): {exc}",
                RuntimeWarning,
                stacklevel=2,
            )
            return
        if chain_mol is None:
            return
        
        # Add chain atoms to chemical space
        chain_idxs = [self.chemical_space.AddAtom(a) for a in chain_mol.GetAtoms()]
        for b in chain_mol.GetBonds():
            self.chemical_space.AddBond(
                chain_idxs[b.GetBeginAtomIdx()],
                chain_idxs[b.GetEndAtomIdx()],
                b.GetBondType()
            )
        
        if not chain_idxs:
            return
        
        # Find head and tail (first and last atoms matching expected types)
        chain_head = chain_idxs[0]
        chain_tail = chain_idxs[-1]
        
        # Check if need bridge on left side
        node_u_symbol = self.chemical_space.GetAtomWithIdx(att_u).GetSymbol()
        head_symbol = chain_head_atom
        
        if self.config.connection.auto_bridge and node_u_symbol == head_symbol:
            # Same atom type - need bridge
            bridge = self.chemical_space.AddAtom(Chem.Atom(self.config.connection.default_bridge_atom))
            self.chemical_space.AddBond(att_u, bridge, Chem.BondType.SINGLE)
            self.chemical_space.AddBond(bridge, chain_head, Chem.BondType.SINGLE)
        else:
            # Direct bond OK
            self.chemical_space.AddBond(att_u, chain_head, Chem.BondType.SINGLE)
        
        # Check if need bridge on right side (usually not for PDMS)
        node_v_symbol = self.chemical_space.GetAtomWithIdx(att_v).GetSymbol()
        tail_symbol = chain_tail_atom
        
        if self.config.connection.auto_bridge and node_v_symbol == tail_symbol:
            # Same atom type - need bridge
            bridge = self.chemical_space.AddAtom(Chem.Atom(self.config.connection.default_bridge_atom))
            self.chemical_space.AddBond(chain_tail, bridge, Chem.BondType.SINGLE)
            self.chemical_space.AddBond(bridge, att_v, Chem.BondType.SINGLE)
        else:
            # Direct bond OK
            self.chemical_space.AddBond(chain_tail, att_v, Chem.BondType.SINGLE)
        
        self.edge_atom_map[edge_id] = chain_idxs
        self.edge_backbone_map[edge_id] = (chain_head, chain_tail)
    
    def _create_chain_from_smiles(self, smiles: str, dp: int):
        """Create a polymer chain from repeating SMILES unit.

        Returns the RDKit molecule for a chain of *dp* repeat units, or
        ``None`` if the chain cannot be constructed.  Raises ``ValueError``
        if the SMILES concatenation fails and falls back to a single monomer
        (which would silently produce an under-length chain).
        """
        from rdkit import Chem

        # Remove trailing O for linking (if present)
        if smiles.endswith("O"):
            unit = smiles[:-1]
            linker = "O"
        else:
            unit = smiles
            linker = ""

        # Build full chain SMILES.
        # The last repeat unit's linker O serves as the chain tail; the
        # auto-bridge in `_build_chain_atomistic` direct-bonds it to the
        # network/end-cap Si node. Do NOT append a trailing "[O]" — that
        # produces a spurious -O-O- peroxide bond at the chain tail.
        if linker:
            chain_smiles = (unit + linker) * dp
        else:
            chain_smiles = unit * dp

        mol = Chem.MolFromSmiles(chain_smiles)
        if mol is not None:
            return Chem.RemoveHs(mol)

        # Chain SMILES failed — do NOT silently return a single monomer.
        # That would embed a 1-unit chain instead of a dp-unit chain with
        # no warning, corrupting the molecular structure.
        raise ValueError(
            f"_create_chain_from_smiles: failed to parse chain SMILES for "
            f"dp={dp}, smiles={smiles!r}.\n"
            f"  Attempted chain SMILES: {chain_smiles!r}\n"
            f"  The monomer SMILES may not support simple string concatenation "
            f"for chain building. Use a monomer with a terminal 'O' linker, "
            f"or check that repeated concatenation produces valid SMILES."
        )
