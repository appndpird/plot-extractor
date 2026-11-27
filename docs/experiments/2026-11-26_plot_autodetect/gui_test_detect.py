"""GUI test of Auto-detect (1a) -> Next to check / Checked / Save (1b) in the real app.

usage: python gui_test_detect.py <case>      (northam | york2026)
"""
import csv
import os
import shutil
import sys
import time
import tkinter as tk

APP = r"E:\Software\PlotExtractor-2026-11-26"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, APP)
sys.path.insert(0, HERE)
os.chdir(APP)
from eval_detect import CASES, rings_of, ref_quad  # noqa: E402

case = sys.argv[1] if len(sys.argv) > 1 else "northam"
C = CASES[case]
OUT = os.path.join(HERE, "gui_out", case)
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT)
cfg = os.path.join(APP, "config.json")
shutil.copyfile(cfg, os.path.join(HERE, "config_backup.json"))
results = []


def check(name, ok, info=""):
    results.append((name, bool(ok), info))
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info else ""), flush=True)


def pump(root, s=0.3):
    t = time.time()
    while time.time() - t < s:
        root.update()
        time.sleep(0.02)


def test(root):
    pe = root._pe
    gd = pe["grid_designer"].designer
    gd.state["suppress_dialogs"] = True
    pe["vars"]["edit_basemap"].set("")       # as a fresh user: 1a sets it on save
    pe["nb"].select(pe["tabs"]["grid_top"])
    pe["nb_grid"].select(pe["tabs"]["grid"])
    pump(root)
    check("1a load ortho", gd.load_ortho(C["ortho"]))
    ref, _ = rings_of(C["ref"])
    q = dict(zip(("BL", "BR", "TR", "TL"), ref_quad(ref)))
    for lab in ("BL", "BR", "TR", "TL"):              # default click order
        gd.add_corner(*q[lab])
    t = time.time()
    r = gd.auto_detect(wait=True)
    pump(root)
    check("1a auto-detect ran", r is not None, "%.1f s" % (time.time() - t))
    if r is None:
        return
    v = gd.vars
    check("1a plots drawn", len(gd.grid) == len(r["plots"]), len(gd.grid))
    check("1a fields filled", v["ranges"].get() == str(r["info"]["n_ranges"])
          and v["rows"].get() == str(r["info"]["n_cols"]),
          "%s x %s, %s x %s m, rows %s" % (v["ranges"].get(), v["rows"].get(), v["plot_w"].get(),
                                           v["plot_l"].get(), v["crop_rows"].get()))
    check("1a corners set", gd.ordered_corners() is not None)
    n_flag = sum(1 for p in gd.grid if p["conf"] < 0.5)
    path = gd.save(True, path=os.path.join(OUT, "grid_auto.shp"))
    pump(root, 1.0)
    check("1a saved grid", path and os.path.exists(path))
    rep = os.path.join(OUT, "grid_auto_detect.csv")
    check("1a saved confidence csv", os.path.exists(rep))
    with open(rep) as f:
        rows = list(csv.DictReader(f))
    check("1a csv flags", sum(1 for x in rows if x["check"] == "yes") == n_flag, n_flag)
    check("1a crop rows file", os.path.exists(os.path.join(OUT, "grid_auto_rows.shp")))
    check("current grid set", pe["vars"]["current_grid"].get() == path)
    check("1b basemap = the 1a ortho", pe["vars"]["edit_basemap"].get() == C["ortho"])
    # ---- 1b
    pump(root, 1.5)
    ed = pe["editor"]["w"]
    check("1b editor loaded the grid", ed is not None)
    api = ed.api
    chk = api["checks"]
    check("1b list of plots to check", len(chk["order"]) == n_flag, len(chk["order"]))
    if n_flag:
        api["next_check"](1)
        pump(root)
        i0 = next(iter(api["sel"]))
        check("1b Next selects a flagged plot", i0 in chk["order"], api["ids"][i0])
        api["next_check"](1)
        i1 = next(iter(api["sel"]))
        check("1b Next moves on", i1 != i0 or n_flag == 1)
        api["next_check"](-1)
        check("1b Shift+N goes back", next(iter(api["sel"])) == i0)
        api["shift"]("width", 0.05)                    # fix it a little
        api["mark_checked"]()
        check("1b Checked removes it", i0 in chk["done"] and len(chk["order"]) - len(chk["done"]) == n_flag - 1)
    outp = api["save"](quiet=True)
    pump(root, 0.5)
    check("1b saved edited grid", outp and os.path.exists(outp))
    rep2 = os.path.splitext(outp)[0] + "_detect.csv"
    check("1b list carried over", os.path.exists(rep2))
    if os.path.exists(rep2):
        with open(rep2) as f:
            rr = list(csv.DictReader(f))
        check("1b checked plot leaves the list", sum(1 for x in rr if x["check"] == "yes") == max(n_flag - 1, 0))
    # the edited grid opened again keeps the list
    pe["load_editor"]()
    pump(root, 1.5)
    ed2 = pe["editor"]["w"]
    if ed2 is not None and n_flag:
        check("1b reopened edited grid keeps the list",
              len([i for i in ed2.api["checks"]["order"]]) == n_flag - 1)


def main():
    import tkinter.messagebox as mb
    for nm in ("showerror", "showinfo", "showwarning"):
        setattr(mb, nm, lambda *a, _n=nm, **k: print("DIALOG %s: %s" % (_n, a), flush=True))
    for nm in ("askyesno", "askokcancel"):
        setattr(mb, nm, lambda *a, _n=nm, **k: (print("DIALOG %s: %s" % (_n, a), flush=True), True)[1])
    import plot_extractor
    real = tk.Tk.mainloop

    def fake(self, n=0):
        try:
            test(self)
        except Exception:
            import traceback
            traceback.print_exc()
            check("no exception", False)
        finally:
            shutil.copyfile(os.path.join(HERE, "config_backup.json"), cfg)
            self.destroy()
    tk.Tk.mainloop = fake
    plot_extractor.main()
    tk.Tk.mainloop = real
    bad = [r for r in results if not r[1]]
    print("\n%d / %d passed" % (len(results) - len(bad), len(results)))


if __name__ == "__main__":
    main()
