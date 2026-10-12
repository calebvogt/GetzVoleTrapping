"""Dataset-structure summaries for the explorer.

  coverage   animals caught per grid per calendar month (None = not trapped)
  turnover   per grid x species x session: animals caught for the first time
             vs. caught before (aligned with population_by_site, by date)
  residence  share of animals still known alive t months after first capture,
             by species x habitat x sex (animals first caught >= HORIZON months
             before their grid's last session, so every animal could be followed)
  recruit    by calendar month: share of animals caught that were new
  breeding   year x month calendar on the unmanipulated grids: share of females
             pregnant and of males in breeding condition (one record per animal
             per session)

    python scripts/structure.py <processed_dir> <web_data_dir>
"""
import json
import sys

import numpy as np
import pandas as pd

from sites import SITES

HABS = ["all", "alfalfa", "bluegrass", "tallgrass"]
HORIZON = 18          # months followed for residence curves
MONTH_DAYS = 30.44
JUV_G = 20            # first-capture mass below this is reported as "small" in turnover tooltips
MIN_N = 5
CONTROL = [s for s, v in SITES.items() if v["treatment"] == "control"]


def load(proc):
    c = pd.read_csv(f"{proc}/captures.csv", parse_dates=["date"], low_memory=False)
    ind = pd.read_csv(f"{proc}/individuals.csv", parse_dates=["first", "last"])
    pop = pd.read_csv(f"{proc}/population_by_site.csv", parse_dates=["date"])
    sess = pd.read_csv(f"{proc}/sessions.csv", parse_dates=["date"])
    c["isex"] = c["ind"].map(ind.set_index("ind")["sex"]).fillna(0).astype(int)
    return c, ind, pop, sess


def coverage(c, sess):
    sess["ym"] = sess["year"] * 12 + sess["month"] - 1
    y0 = int(sess["ym"].min())
    y1 = int(sess["ym"].max())
    caught = (c.merge(sess[["site", "session", "ym"]], on=["site", "session"])
                .groupby(["site", "ym", "species"])["ind"].nunique().unstack(fill_value=0))
    out = {"start": y0, "n": y1 - y0 + 1, "sites": {}}
    for site in sorted(sess["site"].unique()):
        trapped = set(sess.loc[sess["site"] == site, "ym"])
        row = []
        for ym in range(y0, y1 + 1):
            if ym not in trapped:
                row.append(None)
                continue
            v = caught.loc[(site, ym)] if (site, ym) in caught.index else None
            row.append([int(v.get("MO", 0)), int(v.get("MP", 0))] if v is not None else [0, 0])
        out["sites"][str(site)] = row
    return out


def turnover(c, pop):
    s = (c.groupby(["species", "site", "session", "ind"])
          .agg(date=("date", "min"), mass=("mass", "max")).reset_index())
    first = s.groupby("ind")["date"].transform("min")
    s["new"] = s["date"] == first
    s["small"] = s["new"] & (s["mass"] > 0) & (s["mass"] < JUV_G)
    agg = (s.groupby(["species", "site", "session"])
            .agg(caught=("ind", "size"), new=("new", "sum"), small=("small", "sum")).reset_index())
    out = {}
    for (sp, site), p in pop.sort_values("date").groupby(["species", "site"]):
        a = p[["session"]].merge(agg[(agg["species"] == sp) & (agg["site"] == site)], on="session", how="left").fillna(0)
        out.setdefault(str(site), {})[sp] = {"new": [int(v) for v in a["new"]],
                                             "ret": [int(v) for v in a["caught"] - a["new"]],
                                             "small": [int(v) for v in a["small"]]}
    return out


def residence(ind, sess):
    end = sess.groupby("site")["date"].max()
    g = ind[(ind["n_captures"] > 0) & ind["sex"].isin([1, 2])].copy()
    g["end"] = g["site"].map(end)
    g = g[(g["end"] - g["first"]).dt.days >= HORIZON * MONTH_DAYS]
    span = (g["last"] - g["first"]).dt.days / MONTH_DAYS
    g["span"] = span
    out = {}
    for hab in HABS:
        for sp in ["MO", "MP"]:
            for sx in [1, 2]:
                q = g[(g["species"] == sp) & (g["sex"] == sx) & ((g["habitat"] == hab) | (hab == "all"))]
                if len(q) < 20:
                    continue
                out[f"{sp}|{hab}|{sx}"] = {"n": int(len(q)),
                                           "s": [round(float((q["span"] >= t).mean()), 4) for t in range(HORIZON + 1)],
                                           "median_days": int((q["last"] - q["first"]).dt.days.median())}
    return out


def recruit(c):
    s = c.groupby(["species", "habitat", "isex", "session", "ind"]).agg(date=("date", "min")).reset_index()
    s["new"] = s["date"] == s.groupby("ind")["date"].transform("min")
    s["month"] = s["date"].dt.month
    out = {}
    for hab in HABS:
        for sp in ["MO", "MP"]:
            for sx in [1, 2]:
                q = s[(s["species"] == sp) & (s["isex"] == sx) & ((s["habitat"] == hab) | (hab == "all"))]
                rows = []
                for m in range(1, 13):
                    v = q[q["month"] == m]
                    rows.append({"n": int(len(v)), "p": round(float(v["new"].mean()), 3) if len(v) >= 20 else None})
                out[f"{sp}|{hab}|{sx}"] = rows
    return out


def breeding(c, sess):
    q = c[c["site"].isin(CONTROL) & c["isex"].isin([1, 2])]
    s = (q.groupby(["species", "habitat", "site", "session", "ind"])
          .agg(sex=("isex", "first"), repro=("repro", "max")).reset_index()
          .merge(sess[["site", "session", "year", "month"]], on=["site", "session"]))
    s = s[s["repro"].isin([1, 2, 3])]
    years = list(range(int(s["year"].min()), int(s["year"].max()) + 1))
    out = {"years": years, "cells": {}}
    for hab in HABS:
        for sp in ["MO", "MP"]:
            v = s[(s["species"] == sp) & ((s["habitat"] == hab) | (hab == "all"))]
            for sx, hit in [(2, [3]), (1, [2])]:
                w = v[v["sex"] == sx]
                g = w.groupby(["year", "month"])["repro"].agg(n="size", k=lambda r: int(r.isin(hit).sum()))
                grid = []
                for y in years:
                    row = []
                    for m in range(1, 13):
                        if (y, m) in g.index and g.loc[(y, m), "n"] >= MIN_N:
                            row.append([int(g.loc[(y, m), "k"]), int(g.loc[(y, m), "n"])])
                        else:
                            row.append(None)
                    grid.append(row)
                out["cells"][f"{sp}|{hab}|{sx}"] = grid
    return out


def main(proc, out):
    c, ind, pop, sess = load(proc)
    res = {"coverage": coverage(c, sess.copy()), "turnover": turnover(c, pop),
           "residence": residence(ind, sess), "recruit": recruit(c), "breeding": breeding(c, sess),
           "horizon": HORIZON, "juv_g": JUV_G}
    with open(f"{out}/structure.json", "w") as f:
        json.dump(res, f, separators=(",", ":"))
    for k in ["MO|all|1", "MO|all|2", "MP|all|1", "MP|all|2"]:
        r = res["residence"][k]
        print(k, "n", r["n"], "median days", r["median_days"], "known alive at 3/6/12 mo",
              r["s"][3], r["s"][6], r["s"][12])
    return res


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
