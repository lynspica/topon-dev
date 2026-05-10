"""Topon — Streamlit GUI (scaffold).

Run via ``topon gui`` (or ``streamlit run topon/gui/app.py``). The full UI
is a planned future deliverable; this scaffold ships the navigation
skeleton + a working "load and validate a config" page so the CLI flow
isn't the only entry point.

Roadmap (track in internal/DEVELOPMENT_INTERNAL.md):
  - Page: Generate (config builder + run-button + log streaming)
  - Page: Inspect (load .gpickle / .nodes/.edges, render dual graph)
  - Page: simbox (epoxy/amine/POSS sliders + LAMMPS log tail)
  - Page: topro (sequence + n_chains/n_repeats/water + run)
"""
from __future__ import annotations

from pathlib import Path

try:
    import streamlit as st  # type: ignore[import-untyped]
except ImportError as exc:  # pragma: no cover — handled by `topon gui` CLI
    raise SystemExit(
        "streamlit is not installed. Install with `pip install topon[gui]` "
        "or `pip install streamlit`."
    ) from exc


# --------------------------------------------------------------------------
# Page setup
# --------------------------------------------------------------------------
st.set_page_config(
    page_title="topon — polymer & protein network generator",
    page_icon="🧬",
    layout="wide",
)

st.title("topon")
st.caption(
    "Topological polymer & protein network generator for LAMMPS. "
    "**This GUI is a scaffold — most pages are placeholders.**"
)

with st.sidebar:
    st.header("Navigation")
    page = st.radio(
        "Page",
        ["Validate config", "Generate (planned)", "Inspect graph (planned)",
         "simbox (planned)", "topro (planned)"],
    )
    st.divider()
    st.markdown(
        "**Docs**\n\n"
        "- `AGENTS.md` — for AI agents starting fresh\n"
        "- `docs/USAGE.md` — CLI + APIs + recipes\n"
        "- `docs/ARCHITECTURE.md` — package layout"
    )


# --------------------------------------------------------------------------
# Pages
# --------------------------------------------------------------------------
def page_validate_config() -> None:
    """Working page: drag a JSON config in and view the schema-vs-raw split."""
    st.header("Validate a JSON config")
    st.write(
        "Pick or upload a topon JSON config; the page reports which "
        "top-level sections fall under the validated schema vs. the "
        "raw extras (forwarded to `Pipeline(..., raw_config=...)`)."
    )

    uploaded = st.file_uploader("Upload config.json", type=["json"])
    if not uploaded:
        st.info("Drop a config.json above to validate it.")
        return

    import json
    import tempfile

    from topon.config import load_config_full, validate_config

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, encoding="utf-8"
    ) as tmp:
        tmp.write(uploaded.getvalue().decode("utf-8"))
        tmp_path = Path(tmp.name)

    try:
        config, raw = load_config_full(tmp_path)
    except Exception as exc:
        st.error(f"Config failed Pydantic validation: {exc}")
        return

    errors = validate_config(config)
    if errors:
        st.warning("Business-rule warnings:")
        for e in errors:
            st.markdown(f"- {e}")
    else:
        st.success("No business-rule errors.")

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Schema-validated sections")
        st.json(json.loads(config.model_dump_json()))
    with col2:
        st.subheader("Raw extras (forwarded to Pipeline)")
        if raw:
            st.json(raw)
        else:
            st.caption("(none)")


def page_placeholder(name: str) -> None:
    """Stub page until the corresponding feature lands."""
    st.header(f"{name} — planned")
    st.info(
        f"The `{name}` page is on the roadmap. Use the CLI for now; "
        "see `docs/USAGE.md`."
    )


# --------------------------------------------------------------------------
# Dispatch
# --------------------------------------------------------------------------
if page == "Validate config":
    page_validate_config()
else:
    page_placeholder(page)
