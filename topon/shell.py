"""Interactive REPL behind `python -m topon` (with no subcommand on a TTY).

The shell wraps the existing click subcommands — no command logic is
duplicated. Each line typed at the `topon>` prompt is `shlex.split` and
dispatched through `topon.cli.main` with `standalone_mode=False`, so
click handles option parsing, help generation, error messages and exit
codes identically to the one-shot CLI.

Built on `cmd.Cmd` (stdlib) so it works on Windows cmd.exe without
extra dependencies. If `readline` (or `pyreadline3` on Windows) is
installed, you get arrow-key history and basic line editing for free.
"""
from __future__ import annotations

import cmd
import shlex
import sys


_GOODBYE = "bye."


def _build_completions(main) -> list[str]:
    """Read the command names from the click group for tab completion."""
    try:
        return sorted(main.commands.keys())
    except AttributeError:
        return []


class TopOnShell(cmd.Cmd):
    """Single-screen REPL wrapping the topon click CLI."""

    prompt = "topon> "
    ruler = "-"
    doc_header = "Commands (type `help <command>` for details):"
    misc_header = "Shell built-ins:"
    undoc_header = "Other commands:"

    def __init__(self, click_main, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._click_main = click_main
        self._completions = _build_completions(click_main)

    # ---- shell built-ins ----

    def do_exit(self, _arg):
        """Leave the topon shell."""
        print(_GOODBYE)
        return True

    do_quit = do_exit
    do_q = do_exit
    do_EOF = do_exit  # Ctrl-D

    def emptyline(self):
        """Pressing Enter on an empty line should do nothing (default behavior
        in cmd.Cmd is to re-run the last command, which is surprising)."""
        return False

    def do_help(self, arg):
        """help [command] — show the top-level command list, or detail on one command."""
        if arg.strip():
            self._dispatch([arg.strip(), "--help"])
        else:
            self._dispatch(["--help"])

    # ---- fallback dispatch: anything else goes to click ----

    def default(self, line: str):
        text = line.strip()
        if not text:
            return False
        try:
            args = shlex.split(text, posix=False)
        except ValueError as exc:
            print(f"parse error: {exc}", file=sys.stderr)
            return False
        self._dispatch(args)
        return False

    # ---- tab completion: subcommand names on the first token ----

    def completedefault(self, text, line, begidx, endidx):
        # Empty input or extending the first token -> command name completion.
        prefix = line[:begidx].strip()
        if not prefix:
            return [c for c in self._completions if c.startswith(text)]
        return []

    def completenames(self, text, line, begidx, endidx):
        return [c for c in self._completions if c.startswith(text)]

    # ---- private ----

    def _dispatch(self, args: list[str]) -> None:
        """Invoke the click main with `standalone_mode=False` so exits and
        errors land back in the shell instead of killing the process."""
        try:
            self._click_main.main(args, standalone_mode=False)
        except SystemExit:
            # click maps --help and validation errors to SystemExit; swallow
            # so the shell keeps running.
            pass
        except KeyboardInterrupt:
            print("(interrupted)")
        except Exception as exc:  # never let a sub-command crash the shell
            print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)


def run_shell(click_main, intro: str | None = None) -> int:
    """Entry point — drop into the interactive shell.

    Returns the exit code for the process (currently always 0).
    """
    shell = TopOnShell(click_main)
    try:
        shell.cmdloop(intro=intro)
    except KeyboardInterrupt:
        print()
        print(_GOODBYE)
    return 0
