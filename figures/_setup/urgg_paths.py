"""Path helpers for Python figure scripts (same rules as urgg_paths.R).

Order: environment variables URGG_PROJECT / URGG_FIGURES, then config/local_paths.py
(git-ignored), then repository defaults data/ and figures_out/.

    from urgg_paths import urgg_in, urgg_out
    df = pd.read_csv(urgg_in("website_release/web_data/web_votus.tsv"), sep="\\t")
    out = urgg_out("fig6")
"""
import importlib.util, os
from pathlib import Path

def _repo_root(start=Path(__file__).resolve()):
    for d in [start, *start.parents]:
        if (d / ".urgg-root").exists():
            return d
    raise FileNotFoundError("Repository root not found (no .urgg-root file)")

ROOT = _repo_root()
_cfg = {}
_local = ROOT / "config" / "local_paths.py"
if _local.exists():
    spec = importlib.util.spec_from_file_location("local_paths", _local)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    _cfg = {k: getattr(mod, k) for k in ("URGG_PROJECT", "URGG_FIGURES") if hasattr(mod, k)}

URGG_PROJECT = Path(os.environ.get("URGG_PROJECT") or _cfg.get("URGG_PROJECT") or ROOT / "data")
URGG_FIGURES = Path(os.environ.get("URGG_FIGURES") or _cfg.get("URGG_FIGURES") or ROOT / "figures_out")

def urgg_in(rel, must_exist=True):
    p = URGG_PROJECT / rel
    if must_exist and not p.exists():
        raise FileNotFoundError(f"Input not found: {p} (set URGG_PROJECT or config/local_paths.py)")
    return p

def urgg_out(fig):
    d = URGG_FIGURES / fig
    d.mkdir(parents=True, exist_ok=True)
    return d
