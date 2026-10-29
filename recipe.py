"""recipe.py - the record of the options a run used.

Every run writes (or updates) PlotExtractor_recipe.json in its output folder, one
section per step:

    plot_extraction        tab 4 (Metashape project or orthomosaic) and tab 5 (IIQ)
    plot_tiles             tab 6a - fixed-size ground tiles of the plot images
    rows_and_row_tiles     tab 6b - whole rows and fixed-size tiles along each row

A tiles run also copies the plot_extraction section of the run its plot images came
from (the recipe next to RawCrops), so one file in the tiles folder records the whole
chain: how the plots were cut, then how the tiles were cut from them.

A recipe can be loaded back into the tab that wrote it ("Load recipe..."): the
section's "gui" values are the tab's own fields.
"""
from __future__ import annotations

import datetime
import json
import os

APP_VERSION = "2026.10.30"
RECIPE_NAME = "PlotExtractor_recipe.json"


def _plain(v):
    """JSON-safe copy (paths, numbers, bools, lists, dicts; anything else as str)."""
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    return str(v)


def read_recipe(path):
    """The recipe dict from a file or a folder holding PlotExtractor_recipe.json;
    {} when there is none."""
    if path and os.path.isdir(path):
        path = os.path.join(path, RECIPE_NAME)
    try:
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        return d if isinstance(d, dict) else {}
    except Exception:
        return {}


def upstream_extraction(src_dir):
    """plot_extraction section of the run whose plot images are in src_dir
    (recipe in the run folder above RawCrops, or in src_dir itself)."""
    if not src_dir:
        return None
    for p in (os.path.dirname(os.path.abspath(src_dir)), os.path.abspath(src_dir)):
        sec = read_recipe(p).get("plot_extraction")
        if sec:
            return sec
    return None


def write_recipe(out_dir, section, gui=None, params=None, results=None, src_dir=None,
                 log=None):
    """Add / replace `section` in out_dir/PlotExtractor_recipe.json and return the
    file's path. gui = the tab's field values (what "Load recipe" puts back),
    params = the values the run was given, results = counts after the run."""
    try:
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, RECIPE_NAME)
        d = read_recipe(path)
        d["tool"] = "PlotExtractor"
        d["version"] = APP_VERSION
        d["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
        sec = dict(run_at=d["updated"], version=APP_VERSION)
        if gui is not None:
            sec["gui"] = _plain(gui)
        if params is not None:
            sec["params"] = _plain(params)
        if results is not None:
            sec["results"] = _plain(results)
        d[section] = sec
        if src_dir and section != "plot_extraction":
            up = upstream_extraction(src_dir)
            if up:
                d["plot_extraction"] = up
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(d, fh, indent=1)
        if log:
            log(f"Recipe ({section}): {path}")
        return path
    except Exception as e:
        if log:
            log(f"[!] could not write the recipe: {e}")
        return ""


def gui_section(path, section):
    """The 'gui' values of `section` in a recipe file (or folder); {} if absent."""
    return (read_recipe(path).get(section) or {}).get("gui") or {}


def apply_gui(V, values):
    """Put recipe values into a dict of Tk variables; returns how many were set."""
    n = 0
    for k, v in (values or {}).items():
        var = V.get(k)
        if var is None:
            continue
        try:
            var.set(v)
            n += 1
        except Exception:
            pass
    return n
