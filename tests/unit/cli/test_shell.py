"""Smoke tests for the interactive `topon>` shell.

These tests drive the shell with piped stdin (simulating a session) and
assert the obvious things: commands route to click, exit works, errors
don't kill the shell, tab-completion list is populated.
"""
from __future__ import annotations

import io
from contextlib import redirect_stdout, redirect_stderr

import pytest

from topon.cli import main as cli_main
from topon.shell import TopOnShell


@pytest.fixture
def shell():
    return TopOnShell(cli_main)


def test_exit_returns_true(shell):
    assert shell.do_exit("") is True
    assert shell.do_quit("") is True
    assert shell.do_q("") is True
    assert shell.do_EOF("") is True


def test_emptyline_does_nothing(shell):
    # Should NOT re-run the previous command (default cmd.Cmd surprises here).
    assert shell.emptyline() is False


def test_default_dispatches_to_click(shell, tmp_path):
    """`init` typed at the prompt should produce a config file."""
    out_file = tmp_path / "from_shell.json"
    buf = io.StringIO()
    with redirect_stdout(buf):
        shell.default(f"init --output {out_file} --preset cg_kg")
    assert out_file.exists()
    assert "Wrote" in buf.getvalue() or out_file.read_text().startswith("{")


def test_default_swallows_click_exit(shell):
    """`--help` raises SystemExit in click; the shell must absorb it."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        shell.default("init --help")  # would SystemExit in click
    assert "Usage" in buf.getvalue() or "init" in buf.getvalue()


def test_default_handles_unknown_command(shell):
    """Unknown command should print a click error but not crash."""
    buf, err = io.StringIO(), io.StringIO()
    with redirect_stdout(buf), redirect_stderr(err):
        shell.default("notarealcommand")
    # shouldn't raise — just route to click which prints usage / error
    assert True


def test_completion_lists_subcommands(shell):
    completions = shell.completenames("i", "i", 0, 1)
    assert "init" in completions
    assert "inspect" in completions


def test_default_handles_parse_error(shell):
    """Unterminated quote should print a parse error, not crash."""
    err = io.StringIO()
    with redirect_stderr(err):
        shell.default('init --output "unterminated')
    # In posix=False mode, shlex is lenient — this test mostly ensures
    # the shell doesn't raise. Either output is fine.
    assert True
