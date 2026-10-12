"""Coming back to the trap: where recaptures happen and whether voles learn to avoid traps.

Uses grids whose station layout is exact (excludes BG 7 and the Trelease
Prairie grids 12-14) and captures with a recorded check (AM/PM/Eve), ordered
AM < PM < Eve within a day.

  return distance  distance (m) from each capture to the animal's previous one,
                   split into "later in the same session" and "a later session"
  experience       for animals caught in >= COHORT_MIN sessions, mean captures per
                   session in their 1st..COHORT_MIN-th session, and how often the
                   first capture of a session is in the same trap as the last
                   capture of the previous one

    python scripts/traps.py <processed_dir> <web_data_dir>
"""
import json
import sys

import numpy as np
import pandas as pd

EXACT_GRIDS = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 15, 16, 17, 18, 19]
HABS = ["all", "alfalfa", "bluegrass", "tallgrass"]
DIST_EDGES = [0.1, 10.1, 20.1, 30.1]          # same trap | 10 m | 14-20 m | 22-30 m | > 30 m
DIST_LABELS = ["same trap", "next trap (10 m)", "14–20 m", "22–30 m", "> 30 m"]
COHORT_MIN = 5
MIN_N = 20


def load(proc):
    c = pd.read_csv(f"{proc}/captures.csv", parse_dates=["date"], low_memory=False)
    ind = pd.read_csv(f"{proc}/individuals.csv")
    c["isex"] = c["ind"].map(ind.set_index("ind")["sex"])
    g = c[c["site"].isin(EXACT_GRIDS) & c["grid_col"].notna() & c["ampm"].isin([1, 2, 3])
          & c["isex"].isin([1, 2])].copy()
    g["isex"] = g["isex"].astype(int)
    g["ord"] = g["ampm"].map({1: 0, 2: 1, 3: 2})
    return g.sort_values(["ind", "date", "ord"]).reset_index(drop=True)


def returns(g):
    p = g.groupby("ind")
    g["pc"], g["pr"], g["ps"] = p["grid_col"].shift(), p["grid_row"].shift(), p["session"].shift()
    r = g.dropna(subset=["pc"]).copy()
    r["dist"] = 10 * np.hypot(r["grid_col"] - r["pc"], r["grid_row"] - r["pr"])
    r["within"] = r["session"] == r["ps"]
    r["dbin"] = np.searchsorted(DIST_EDGES, r["dist"].to_numpy(), side="left")
    return r


def experience(g, r):
    s = (g.groupby(["ind", "session"])
          .agg(n=("date", "size"), d=("date", "min"), species=("species", "first"),
               habitat=("habitat", "first"), sex=("isex", "first"))
          .reset_index().sort_values(["ind", "d"]))
    s["k"] = s.groupby("ind").cumcount() + 1
    s["total"] = s.groupby("ind")["k"].transform("max")
    # first capture of each later session vs. last capture of the previous one
    first = r[~r["within"]][["ind", "session", "dist"]]
    s = s.merge(first, on=["ind", "session"], how="left")
    return s


def summarize(r, s):
    out = {"labels": DIST_LABELS, "cohort_min": COHORT_MIN, "dist": {}, "exp": {}}
    for hab in HABS:
        for sp in ["MO", "MP"]:
            for sx in [1, 2]:
                key = f"{sp}|{hab}|{sx}"
                q = r[(r["species"] == sp) & (r["isex"] == sx) & ((r["habitat"] == hab) | (hab == "all"))]
                d = {}
                for name, w in [("within", True), ("next", False)]:
                    v = q[q["within"] == w]
                    if len(v) >= MIN_N:
                        d[name] = {"n": int(len(v)),
                                   "p": [round(float((v["dbin"] == i).mean()), 3) for i in range(len(DIST_LABELS))],
                                   "median": float(v["dist"].median())}
                out["dist"][key] = d
                e = s[(s["species"] == sp) & (s["sex"] == sx) & ((s["habitat"] == hab) | (hab == "all"))
                      & (s["total"] >= COHORT_MIN) & (s["k"] <= COHORT_MIN)]
                rows = []
                for k in range(1, COHORT_MIN + 1):
                    v = e[e["k"] == k]
                    vd = v.dropna(subset=["dist"])
                    rows.append({"n": int(len(v)),
                                 "caps": round(float(v["n"].mean()), 2) if len(v) >= MIN_N else None,
                                 "same": round(float((vd["dist"] == 0).mean()), 3) if len(vd) >= MIN_N else None})
                out["exp"][key] = {"animals": int(e["ind"].nunique()), "rows": rows}
    return out


def main(proc, out):
    g = load(proc)
    r = returns(g)
    s = experience(g, r)
    res = summarize(r, s)
    with open(f"{out}/traps.json", "w") as f:
        json.dump(res, f, separators=(",", ":"))
    for k in ["MO|all|1", "MO|all|2", "MP|all|1", "MP|all|2"]:
        print(k, res["dist"][k]["within"]["p"][0], res["dist"][k]["next"]["p"][0],
              [x["caps"] for x in res["exp"][k]["rows"]], res["exp"][k]["animals"])
    return res


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
