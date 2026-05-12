"""
Command-line interface for Topon.
"""

import json
import sys
from pathlib import Path

import click

from topon import __version__


_BANNER = r"""
   +==================================================================+
   |     ########   ########   ########    ########   ###     ##      |
   |        ##      ##    ##   ##     ##   ##    ##   ####    ##      |
   |        ##      ##    ##   ########    ##    ##   ## ##   ##      |
   |        ##      ##    ##   ##          ##    ##   ##  ##  ##      |
   |        ##      ########   ##          ########   ##   #####      |
   |                                                                  |
   |   Topological polymer & protein network generator for LAMMPS     |
   |                                                       v{version:<7s}    |
   +==================================================================+

   Get started in 3 commands:
     init                  Make a starter config that works out of the box
     doctor <config>       Lint it for common mistakes
     generate <config>     Run the 6-stage pipeline -> LAMMPS data files

   All commands (type `help <cmd>` for details):
     init       validate    doctor      generate    inspect
     analyze    simbox      chain       topro       recipes
     gui

   Pipeline (`generate`):
      Topology -> Analysis -> Assignment -> Chemistry -> Conformation -> Output

   Quick start:
     topon> init                       # writes config.json (atomistic PDMS)
     topon> doctor config.json         # checks for footguns
     topon> generate config.json       # builds LAMMPS files
     topon> inspect <output_dir>       # summarises what landed
     topon> recipes                    # "I want X -> run Y" cheatsheet

   Shell built-ins:
     help [cmd]            Top-level command list, or detail on one command
     exit | quit | Ctrl-D  Leave the shell

   Docs: AGENTS.md  |  docs/USAGE.md  |  examples/demos/  |  examples/workflows/
"""


_BANNER_ONE_SHOT_FOOTER = r"""
   Tips:
     `python -m topon`            drop into the interactive shell
     `python -m topon <command>`  run a single command and exit
     `python -m topon <cmd> --help`  flags for one command
"""


@click.group(invoke_without_command=True)
@click.version_option(version=__version__, prog_name="topon")
@click.option(
    "--no-shell", is_flag=True,
    help="Print the banner and exit (don't drop into the interactive shell).",
)
@click.pass_context
def main(ctx, no_shell: bool):
    """Topon: topological polymer + protein network generator for LAMMPS.

    With no subcommand on a TTY, drops into an interactive shell where
    you can type `help`, `init`, `doctor`, etc. directly. Pipe input or
    pass --no-shell to get the old banner-and-exit behaviour.
    """
    if ctx.invoked_subcommand is not None:
        return

    banner = _BANNER.format(version=__version__)
    interactive = sys.stdin.isatty() and not no_shell

    if not interactive:
        click.echo(banner)
        click.echo(_BANNER_ONE_SHOT_FOOTER)
        return

    # Interactive shell. Banner first, then drop into the REPL.
    from topon.shell import run_shell
    intro = banner + "\n   Type `help` for the command list, or `exit` to leave.\n"
    sys.exit(run_shell(main, intro=intro))


@main.command()
@click.argument("config_path", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Override output directory")
@click.option("--dry-run", is_flag=True, help="Validate config without running pipeline")
@click.option(
    "--export-graphml",
    is_flag=True,
    help="Also export the dual-graph (chains + entanglement edges) as <name>.graphml",
)
@click.option(
    "--export-npz",
    is_flag=True,
    help="Also export the graph as <name>.npz for downstream GNN pipelines (planned; falls back gracefully if writer not yet implemented)",
)
def generate(
    config_path: str,
    output: str,
    dry_run: bool,
    export_graphml: bool,
    export_npz: bool,
):
    """
    Run the full pipeline from a configuration file.

    CONFIG_PATH: Path to the JSON configuration file.
    """
    from topon.config import load_config_full, validate_config
    from topon.pipeline import Pipeline

    click.echo(f"Loading configuration from: {config_path}")

    try:
        config, raw_cfg = load_config_full(config_path)
    except Exception as e:
        click.echo(f"Error loading configuration: {e}", err=True)
        sys.exit(1)

    # Override output directory if specified
    if output:
        config.study.output_dir = output

    # CLI flags override config.output.*
    if export_graphml:
        config.output.export_graphml = True
    if export_npz:
        config.output.export_npz = True

    # Validate configuration
    errors = validate_config(config)
    if errors:
        click.echo("Configuration validation errors:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)

    click.echo("Configuration is valid.")
    if raw_cfg:
        click.echo(
            f"  Raw extras forwarded to Pipeline: {sorted(raw_cfg.keys())}"
        )

    if dry_run:
        click.echo("Dry run - not executing pipeline.")
        return

    # Run pipeline (raw_config carries conformation / simulation /
    # execution / experimental sections that aren't in ToponConfig).
    click.echo("Running pipeline...")
    pipeline = Pipeline(config, raw_config=raw_cfg)
    pipeline.run()

    click.echo(f"Pipeline complete. Output written to: {config.study.output_dir}")


@main.command()
@click.option("--port", default=8501, type=int, help="Streamlit server port (default 8501)")
def gui(port: int):
    """Launch the Streamlit GUI (scaffold; install with `pip install topon[gui]`)."""
    import shutil
    import subprocess

    try:
        import streamlit  # noqa: F401
    except ImportError:
        click.echo(
            "streamlit is not installed. Install with `pip install topon[gui]` "
            "or `pip install streamlit`, then re-run `topon gui`.",
            err=True,
        )
        sys.exit(1)

    app_path = Path(__file__).parent / "gui" / "app.py"
    if not app_path.exists():
        click.echo(f"GUI app not found at: {app_path}", err=True)
        sys.exit(1)

    streamlit_bin = shutil.which("streamlit")
    if streamlit_bin is None:
        click.echo("streamlit binary not on PATH; falling back to `python -m streamlit`.")
        cmd = [sys.executable, "-m", "streamlit", "run", str(app_path),
               "--server.port", str(port)]
    else:
        cmd = [streamlit_bin, "run", str(app_path), "--server.port", str(port)]

    click.echo(f"Launching topon GUI on http://localhost:{port}/ ...")
    sys.exit(subprocess.call(cmd))


@main.command()
@click.argument("config_path", type=click.Path(exists=True))
def validate(config_path: str):
    """
    Validate a configuration file without running the pipeline.

    CONFIG_PATH: Path to the JSON configuration file.
    """
    from topon.config import validate_config
    from topon.utils.errors import load_config_or_die

    config, _raw = load_config_or_die(config_path)

    errors = validate_config(config)

    if errors:
        click.echo("Configuration validation errors:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)
    else:
        click.echo("Configuration is valid!")


@main.command()
@click.argument("config_path", type=click.Path(exists=True))
@click.option(
    "--strict", is_flag=True,
    help="Treat WARNs as errors (exit 1 if any warn-or-worse fires).",
)
def doctor(config_path: str, strict: bool):
    """Lint a config for known footguns beyond what `validate` checks.

    Where `topon validate` is a Pydantic schema check, `topon doctor` runs
    a small registry of semantic rules sourced from
    `internal/DEVELOPMENT_INTERNAL.md` known issues + things we've watched
    new users trip on. Each rule prints one of:

    \b
        [ok]    ... informational
        [warn]  ... likely surprise; recommended fix included
        [error] ... will crash or silently produce wrong output

    Exit code: 0 if no errors (and no warns when --strict), 1 otherwise.
    """
    from topon.diagnostics import run_all_rules
    from topon.utils.errors import load_config_or_die

    cfg, raw = load_config_or_die(config_path)

    issues = run_all_rules(cfg, raw)
    if not issues:
        click.echo("No issues found.")
        return

    n_err = sum(1 for i in issues if i.level == "error")
    n_warn = sum(1 for i in issues if i.level == "warn")
    n_ok = sum(1 for i in issues if i.level == "ok")

    icon = {"ok": "[ok]   ", "warn": "[warn] ", "error": "[error]"}
    for issue in issues:
        click.echo(f"{icon[issue.level]} {issue.rule}: {issue.message}")
        if issue.fix:
            click.echo(f"        fix: {issue.fix}")

    click.echo()
    click.echo(f"Summary: {n_err} error / {n_warn} warn / {n_ok} ok")
    if n_err or (strict and n_warn):
        sys.exit(1)


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


_PRESET_MAP = {
    "atomistic_pdms": "examples/demos/polymer/atomistic/basic/config.json",
    "cg_kg":         "examples/demos/polymer/coarse_grained/basic/config.json",
    "poss":          "examples/demos/poss/config.json",
    "martini_resilin": None,  # MARTINI uses a separate CLI; print instructions instead
    "charmm_resilin":  None,  # ditto
}


def _resolve_preset_path(preset: str) -> Path:
    """Locate a preset config.json by walking up from the topon package."""
    rel = _PRESET_MAP[preset]
    if rel is None:
        return None
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        candidate = parent / rel
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find bundled preset '{preset}' (looked for {rel}).")


def _interactive_config() -> dict:
    """Prompt the user for the 5-6 knobs that actually vary; return a raw dict.

    The dict is in the on-disk JSON shape (not a Pydantic object); demos already
    use this shape so the result is a drop-in replacement for any
    examples/demos/.../config.json.
    """
    click.echo()
    click.echo("topon init — interactive config builder")
    click.echo("Press Enter to accept the default in [brackets].")
    click.echo()

    name = click.prompt("Study name", default="my_run").strip()
    out_dir = click.prompt("Output directory", default=f"output_{name}").strip()

    model = click.prompt(
        "Chemistry model",
        type=click.Choice(["atomistic", "coarse_grained"]),
        default="atomistic",
    )
    lattice_type = click.prompt(
        "Lattice type",
        type=click.Choice(["SC", "BCC", "FCC"]),
        default="SC",
    )
    lattice_size = click.prompt("Lattice size (NxNxN)", default="5x5x5").strip()
    max_func = click.prompt("Max functionality", type=int, default=4)
    dp_mean = click.prompt("Degree of polymerization (mean)", type=int, default=10)
    density = click.prompt(
        "Target density",
        type=float,
        default=(1.1 if model == "atomistic" else 0.85),
    )

    return {
        "study": {"name": name, "output_dir": out_dir},
        "topology": {
            "source": "generate",
            "generator": {
                "lattice_size": lattice_size,
                "lattice_type": lattice_type,
                "max_functionality": max_func,
                "degree_distribution": "0:0,1:25",
            },
        },
        "chemistry": {
            "model_type": model,
            "target_density": density,
        },
        "assignment": {
            "dp_distribution": {"default": {"mean": float(dp_mean), "pdi": 1.0}},
        },
        "conformation": {"overlap_cutoff": 0.2, "overlap_max_iters": 20},
        "simulation": {"run_steps": 2000},
        "execution": {"auto_run": False, "executable": "lmp", "n_procs": 1},
    }


@main.command()
@click.option("--output", "-o", type=click.Path(), default="config.json",
              show_default=True, help="Path for the new config file")
@click.option(
    "--preset",
    type=click.Choice(list(_PRESET_MAP.keys())),
    default=None,
    help="Start from a bundled demo (atomistic_pdms / cg_kg / poss / "
         "martini_resilin / charmm_resilin). Without --preset and without "
         "--interactive, writes a working atomistic_pdms-style starter.",
)
@click.option(
    "--interactive", "-i",
    is_flag=True,
    help="Prompt for the 5-6 knobs that actually vary and write the result.",
)
def init(output: str, preset: str, interactive: bool):
    """Create a starter config.json that runs as-is.

    Three modes:

      \b
      topon init                    fast path - copy the atomistic_pdms preset
      topon init --preset cg_kg     start from a different bundled demo
      topon init --interactive      prompt-driven, walks through 5-6 knobs

    The MARTINI and CHARMM protein-network paths use a separate CLI
    (`python -m topon.protein_network ...`); `--preset martini_resilin` /
    `charmm_resilin` print the right invocation rather than write a JSON.
    """
    import json as _json
    import shutil

    if preset in ("martini_resilin", "charmm_resilin"):
        click.echo()
        if preset == "martini_resilin":
            click.echo("MARTINI 3 protein networks use a separate CLI. Try:")
            click.echo()
            click.echo("  python -m topon.protein_network generate \\")
            click.echo("      --block-seq GGRPSDSYGAPGGGN \\")
            click.echo("      --n-chains 4 --n-repeats 6 \\")
            click.echo("      --water-density 0 --output runs/resilin_dry --seed 42")
        else:
            click.echo("CHARMM36m protein networks have a worked example at:")
            click.echo()
            click.echo("  python examples/demos/protein/charmm/run.py "
                       "--output runs/charmm_demo")
        click.echo()
        click.echo("See examples/demos/protein/{martini,charmm}/README.md for details.")
        return

    out_path = Path(output)
    if out_path.exists():
        if not click.confirm(f"{output} already exists. Overwrite?", default=False):
            click.echo("Aborted; no file written.")
            sys.exit(1)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if interactive:
        data = _interactive_config()
        out_path.write_text(_json.dumps(data, indent=2), encoding="utf-8")
        click.echo()
        click.echo(f"Wrote {out_path}")
    else:
        chosen = preset or "atomistic_pdms"
        src = _resolve_preset_path(chosen)
        shutil.copy2(src, out_path)
        click.echo(f"Wrote {out_path} (preset: {chosen}, copied from {src.name})")

    click.echo()
    click.echo("Next steps:")
    click.echo(f"  topon doctor {out_path}        # lint for known footguns")
    click.echo(f"  topon validate {out_path}      # schema check")
    click.echo(f"  topon generate {out_path}      # run the 6-stage pipeline")


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
@click.argument("run_dir", type=click.Path(exists=True, file_okay=False))
def inspect(run_dir: str):
    """Summarise what's inside a `topon generate` output directory.

    RUN_DIR: path either to the study folder (containing 02_Chemistry/,
    03_Conformation/, 04_Simulation/) or to its parent. Prints atom counts,
    box dimensions, which displacement files landed, and the next LAMMPS
    commands to run. Useful for confirming a long pipeline finished cleanly
    without grepping `system.data` headers by hand.
    """
    from topon.analysis.run_summary import summarise, format_summary

    summary = summarise(Path(run_dir))
    click.echo(format_summary(summary))


@main.command()
def shell():
    """Drop into the interactive `topon>` shell explicitly.

    Same as running `python -m topon` on a real terminal; use this
    when you want to force the REPL even if stdin isn't a TTY (e.g.
    inside a script or from an IDE terminal that confuses isatty).
    """
    from topon.shell import run_shell
    banner = _BANNER.format(version=__version__)
    intro = banner + "\n   Type `help` for the command list, or `exit` to leave.\n"
    sys.exit(run_shell(main, intro=intro))


@main.command(
    name="topro",
    context_settings={
        "ignore_unknown_options": True,
        "allow_extra_args": True,
        "help_option_names": [],  # let inner argparse handle --help
    },
)
@click.argument("topro_args", nargs=-1, type=click.UNPROCESSED)
def topro(topro_args):
    """Protein-network sub-system (MARTINI 3 + CHARMM36m via topon.protein_network).

    Passes all arguments through to the existing argparse CLI:

    \b
        topon topro generate --block-seq ... --n-chains 4 ...
        topon topro sweep    --water-densities 0,4,8 ...
        topon topro topology --n-chains 16 ...

    Run `topon topro --help` for the full subcommand list and flags.
    For CHARMM atomistic, the recipe is `python -m topon.protein_network.charmm.build_systems ...`
    (see `topon recipes`).
    """
    from topon.protein_network.cli import main as topro_main
    sys.exit(topro_main(list(topro_args)))


@main.command()
def recipes():
    """Print a 'I want X -> run Y' table for the most common use cases.

    Quick orientation for new users: which subcommand or script handles
    which kind of network. Edit `topon/cli.py:_RECIPES` to add rows.
    """
    rows = [
        ("Polymer network from a JSON config",
         "topon init && topon doctor config.json && topon generate config.json"),
        ("Atomistic PDMS network",
         "topon init --preset atomistic_pdms"),
        ("Coarse-grained Kremer-Grest network",
         "topon init --preset cg_kg"),
        ("POSS chain-cap demo",
         "topon init --preset poss"),
        ("MARTINI 3 protein network (resilin)",
         "python -m topon.protein_network generate "
         "--block-seq GGRPSDSYGAPGGGN --n-chains 4 --n-repeats 6 \\\n"
         "      --water-density 0 --output runs/resilin_dry"),
        ("CHARMM36m atomistic protein network",
         "python examples/demos/protein/charmm/run.py --output runs/charmm"),
        ("Crosslink simulation box (no graph)",
         "topon simbox --n-epoxy 50 --n-amino 25 --n-poss 10"),
        ("Single polymer chain in solvent",
         "topon chain --chain-smiles \"[Si](C)(C)O\" --dp 50"),
        ("Batch of 25 lattice graphs + CSV summary",
         "python examples/workflows/batch_polymer_topology/run.py"),
        ("BFM gel-point parameter sweep",
         "python examples/workflows/bfm_gel_point_sweep/run.py"),
        ("Inspect a finished run directory",
         "topon inspect <run_dir>"),
        ("Graph statistics (no chemistry build)",
         "topon analyze graph.gpickle"),
    ]
    click.echo()
    click.echo("topon recipes — common use cases")
    click.echo("=" * 70)
    for goal, cmd in rows:
        click.echo()
        click.echo(f"  {goal}")
        for line in cmd.split("\n"):
            click.echo(f"    {line}")
    click.echo()
    click.echo("More: docs/USAGE.md, examples/demos/, examples/workflows/")


if __name__ == "__main__":
    main()
