"""
Command-line interface for Topon.
"""

import json
import sys
from pathlib import Path

import click

from topon import __version__


@click.group()
@click.version_option(version=__version__, prog_name="topon")
def main():
    """Topon: Polymer network generation for molecular simulations."""
    pass


@main.command()
@click.argument("config_path", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Override output directory")
@click.option("--dry-run", is_flag=True, help="Validate config without running pipeline")
def generate(config_path: str, output: str, dry_run: bool):
    """
    Run the full pipeline from a configuration file.
    
    CONFIG_PATH: Path to the JSON configuration file.
    """
    from topon.config import load_config, validate_config
    from topon.pipeline import Pipeline
    
    click.echo(f"Loading configuration from: {config_path}")
    
    try:
        config = load_config(config_path)
    except Exception as e:
        click.echo(f"Error loading configuration: {e}", err=True)
        sys.exit(1)
    
    # Override output directory if specified
    if output:
        config.study.output_dir = output
    
    # Validate configuration
    errors = validate_config(config)
    if errors:
        click.echo("Configuration validation errors:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)
    
    click.echo("Configuration is valid.")
    
    if dry_run:
        click.echo("Dry run - not executing pipeline.")
        return
    
    # Run pipeline
    click.echo("Running pipeline...")
    pipeline = Pipeline(config)
    pipeline.run()
    
    click.echo(f"Pipeline complete. Output written to: {config.study.output_dir}")


@main.command()
@click.argument("config_path", type=click.Path(exists=True))
def validate(config_path: str):
    """
    Validate a configuration file without running the pipeline.
    
    CONFIG_PATH: Path to the JSON configuration file.
    """
    from topon.config import load_config, validate_config
    
    try:
        config = load_config(config_path)
    except Exception as e:
        click.echo(f"Error loading configuration: {e}", err=True)
        sys.exit(1)
    
    errors = validate_config(config)
    
    if errors:
        click.echo("Configuration validation errors:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)
    else:
        click.echo("Configuration is valid!")


@main.command()
@click.argument("graph_path", type=click.Path(exists=True))
@click.option("--format", "-f", type=click.Choice(["text", "json"]), default="text")
@click.option("--nodes", type=click.Path(exists=True), default=None,
              help="Companion .nodes file (when GRAPH_PATH is a .edges file)")
def analyze(graph_path: str, format: str, nodes: str):
    """
    Analyze a topology graph and report statistics.

    GRAPH_PATH: Path to a .gpickle file, or a .nodes file (pair it with --nodes
    when passing a .edges file).

    Examples:

        topon analyze network.gpickle

        topon analyze network.nodes

        topon analyze network.edges --nodes network.nodes
    """
    from topon.topology.loader import load_graph
    from topon.analysis.report import analyze_graph

    p = Path(graph_path)
    try:
        if p.suffix == ".gpickle":
            G, dims = load_graph(gpickle_path=graph_path)
        elif p.suffix == ".nodes":
            edges_path = str(p.with_suffix(".edges"))
            if not Path(edges_path).exists():
                click.echo(f"Error: companion .edges file not found: {edges_path}", err=True)
                sys.exit(1)
            G, dims = load_graph(nodes_path=graph_path, edges_path=edges_path)
        elif p.suffix == ".edges":
            if not nodes:
                nodes_path = str(p.with_suffix(".nodes"))
                if not Path(nodes_path).exists():
                    click.echo("Error: provide --nodes <path> for a .edges file.", err=True)
                    sys.exit(1)
                nodes = nodes_path
            G, dims = load_graph(nodes_path=nodes, edges_path=graph_path)
        else:
            click.echo(f"Error: unsupported file type '{p.suffix}'. Use .gpickle, .nodes, or .edges.", err=True)
            sys.exit(1)
    except Exception as e:
        click.echo(f"Error loading graph: {e}", err=True)
        sys.exit(1)

    report = analyze_graph(G, dims, verbose=(format == "text"))

    if format == "json":
        click.echo(json.dumps(report, indent=2))


@main.command()
@click.option("--output", "-o", type=click.Path(), default="config.json")
@click.option("--full", is_flag=True, help="Include all options with defaults")
def init(output: str, full: bool):
    """
    Create a new configuration file with defaults.
    """
    from topon.config.loader import create_default_config, save_config
    
    config = create_default_config()
    save_config(config, output)
    
    click.echo(f"Created configuration file: {output}")


@main.command()
@click.option("--output", "-o", type=click.Path(), default="simbox_output",
              show_default=True, help="Output directory for LAMMPS files")
@click.option("--n-epoxy", type=int, default=50, show_default=True,
              help="Number of Epoxy-PDMS molecules")
@click.option("--n-amino", type=int, default=25, show_default=True,
              help="Number of Amino-PDMS molecules")
@click.option("--n-poss", type=int, default=10, show_default=True,
              help="Number of AM0270-POSS molecules")
@click.option("--density", type=float, default=0.85, show_default=True,
              help="Target packing density (g/cm³)")
@click.option("--seed", type=int, default=42, show_default=True,
              help="Random seed for reproducible packing")
def simbox(output: str, n_epoxy: int, n_amino: int, n_poss: int,
           density: float, seed: int):
    """
    Pack a crosslink simulation box and write LAMMPS input files.

    Builds Epoxy-PDMS, Amino-PDMS, and AM0270-POSS molecules, packs them
    into a periodic box at the target density, and writes DREIDING-
    parameterised LAMMPS data + input scripts ready to run.

    Example:

        topon simbox --output my_system --n-epoxy 600 --n-amino 300

    Then run LAMMPS:

        cd my_system && lmp -in 1_minimize.in
    """
    from topon.simbox.workflow import run_workflow

    try:
        run_workflow(
            output_dir=output,
            n_epoxy=n_epoxy,
            n_amino=n_amino,
            n_poss=n_poss,
            density=density,
            seed=seed,
            verbose=True,
        )
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.option("--output", "-o", type=click.Path(), default="chain_output",
              show_default=True, help="Output directory for LAMMPS files")
@click.option("--chain-smiles", required=True,
              help="SMILES for the polymer repeat unit (e.g. \"[Si](C)(C)O\" for PDMS)")
@click.option("--dp", type=int, required=True,
              help="Degree of polymerization (number of repeat units)")
@click.option("--solvent-smiles", default=None,
              help="SMILES for single solvent (default: toluene). Ignored if --solvent-mixture is set.")
@click.option("--n-solvent", type=int, default=None,
              help="Number of solvent molecules (auto if omitted)")
@click.option("--solvent-mixture", default=None,
              help='Multi-solvent JSON: \'{"smiles":"...","weight_fraction":0.5}\'')
@click.option("--graft-density", type=float, default=0.0, show_default=True,
              help="Graft density: probability of side-chain attachment per backbone unit (0–1)")
@click.option("--graft-smiles", default=None,
              help="SMILES for graft repeat unit (required if --graft-density > 0)")
@click.option("--graft-dp", type=int, default=5, show_default=True,
              help="Number of repeat units per side chain")
@click.option("--density", type=float, default=0.85, show_default=True,
              help="Target packing density (g/cm³)")
@click.option("--seed", type=int, default=42, show_default=True,
              help="Random seed")
def chain(
    output, chain_smiles, dp, solvent_smiles, n_solvent,
    solvent_mixture, graft_density, graft_smiles, graft_dp, density, seed,
):
    """
    Build a single polymer chain in solvent and write DREIDING LAMMPS files.

    The chain is built as a linear atomistic polymer from the given repeat
    unit SMILES and packed with the specified solvent in a periodic box.
    Optional side-chain grafts are supported.

    Examples:

        # PDMS chain in toluene
        topon chain --chain-smiles "[Si](C)(C)O" --dp 20 \\
                    --solvent-smiles "Cc1ccccc1" --n-solvent 200

        # Grafted chain
        topon chain --chain-smiles "[Si](C)(C)O" --dp 30 \\
                    --graft-density 0.1 --graft-smiles "[Si](C)(C)O" --graft-dp 5 \\
                    --solvent-smiles "Cc1ccccc1" --n-solvent 150

    Then run LAMMPS:

        cd chain_output && lmp -in 1_minimize.in
    """
    from topon.singlechain.workflow import run_workflow

    if graft_density > 0 and not graft_smiles:
        click.echo("Error: --graft-smiles is required when --graft-density > 0", err=True)
        sys.exit(1)

    # Parse --solvent-mixture JSON if provided
    parsed_mixture = None
    if solvent_mixture:
        try:
            parsed_mixture = json.loads(solvent_mixture)
            if isinstance(parsed_mixture, dict):
                parsed_mixture = [parsed_mixture]  # wrap single entry
        except json.JSONDecodeError as e:
            click.echo(f"Error parsing --solvent-mixture JSON: {e}", err=True)
            sys.exit(1)

    try:
        result = run_workflow(
            output_dir=output,
            chain_smiles=chain_smiles,
            dp=dp,
            solvent_smiles=solvent_smiles,
            n_solvent=n_solvent,
            solvent_mixture=parsed_mixture,
            graft_density=graft_density,
            graft_smiles=graft_smiles,
            graft_dp=graft_dp,
            density=density,
            seed=seed,
            verbose=True,
        )
        click.echo(f"Chain atoms  : {result['chain_atoms']}")
        click.echo(f"Box length   : {result['box_length_ang']:.2f} Å")
        click.echo(f"Data file    : {result.get('data', result.get('data_file', ''))}")
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.option("--port", "-p", type=int, default=8501)
def gui(port: int):
    """
    Launch the Streamlit GUI for interactive configuration.
    """
    try:
        import streamlit.web.cli as stcli
    except ImportError:
        click.echo("Streamlit not installed. Install with: pip install topon[gui]", err=True)
        sys.exit(1)
    
    gui_path = Path(__file__).parent / "gui" / "app.py"
    
    if not gui_path.exists():
        click.echo("GUI module not yet implemented.", err=True)
        sys.exit(1)
    
    sys.argv = ["streamlit", "run", str(gui_path), "--server.port", str(port)]
    stcli.main()


if __name__ == "__main__":
    main()
