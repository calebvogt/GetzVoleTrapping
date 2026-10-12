"""Export the processed tables as compact JSON for the web explorer.

    python scripts/export_web.py <processed_dir> <web_data_dir>

<web_data_dir> is calebvogt.github.io/projects/getz-vole-trapping/data.

Writes meta.json (sites, population series), mass.json (body-mass summaries)
and one site_NN.json per study site (individuals + captures as columnar arrays; dates are days since
1970-01-01).
"""
import json
import os
import sys
from datetime import date

import numpy as np
import pandas as pd

EPOCH = pd.Timestamp("1970-01-01")


def days(s):
    return ((pd.to_datetime(s) - EPOCH).dt.days).astype("Int64")


def col(s, fill=-1):
    """Series -> plain list of ints with NaN replaced by `fill`."""
    return [int(v) for v in pd.to_numeric(s, errors="coerce").fillna(fill)]


MASS_BIN = 2      # g, histogram bin width
MASS_MAX = 80
GROWTH_DAYS = 400 # trajectories are clipped to this many days after first weighing


def mass_summaries(caps):
    """Body-mass summaries from the monthly grid captures (25-year study).

    One value per animal per session (mean mass, highest reproductive code),
    so animals caught several times in a session are not over-counted.
    """
    c = caps[(caps["mass"] > 0) & (caps["mass"] <= 100)]
    # "Pregnant" is a female-only code; the ~20 male records carrying it are entry errors.
    c = c[~((c["repro"] == 3) & (c["sex"] == 1))]
    s = (c.groupby(["ind", "session"])
          .agg(species=("species", "first"), habitat=("habitat", "first"), sex=("sex", "max"),
               date=("date", "min"), mass=("mass", "mean"), repro=("repro", "max"))
          .reset_index())
    s["month"] = s["date"].dt.month
    habs = ["all", "alfalfa", "bluegrass", "tallgrass"]
    out = {"bin": MASS_BIN, "max": MASS_MAX, "n_records": int(len(s)), "season": {}, "repro": {},
           "growth": {}}

    for sp in ["MO", "MP"]:
        for hab in habs:
            for sex in [0, 1, 2]:
                q = s[(s["species"] == sp) & ((s["habitat"] == hab) | (hab == "all"))
                      & ((s["sex"] == sex) | (sex == 0))]
                key = f"{sp}|{hab}|{sex}"
                rows = []
                for mth in range(1, 13):
                    m = q[q["month"] == mth]
                    known = m[m["repro"].isin([1, 2, 3])]
                    qs = m["mass"].quantile([.1, .25, .5, .75, .9]).round(1).tolist() if len(m) >= 5 else None
                    rows.append({"n": int(len(m)), "q": qs,
                                 "pr": round(float(known["repro"].isin([2, 3]).mean()), 3) if len(known) >= 5 else None})
                out["season"][key] = rows
                bins = np.arange(0, MASS_MAX + MASS_BIN, MASS_BIN)
                for rc in [1, 2, 3]:
                    v = q[q["repro"] == rc]["mass"].clip(upper=MASS_MAX - 0.01)
                    if len(v) < 10:
                        continue
                    h, _ = np.histogram(v, bins=bins)
                    out["repro"][f"{key}|{rc}"] = {"n": int(len(v)), "h": [int(x) for x in h],
                                                   "med": round(float(v.median()), 1)}

    # Growth trajectories: animals weighed in 3+ sessions.
    g = s.sort_values("date")
    g["t"] = (g["date"] - g.groupby("ind")["date"].transform("min")).dt.days
    g = g[g["t"] <= GROWTH_DAYS]
    keep = g.groupby("ind")["session"].transform("size") >= 3
    g = g[keep]
    hab_i = {h: i for i, h in enumerate(habs[1:])}
    for sp in ["MO", "MP"]:
        q = g[g["species"] == sp]
        traj = {"hab": [], "sex": [], "t": [], "m": []}
        for _, a in q.groupby("ind"):
            traj["hab"].append(hab_i[a["habitat"].iloc[0]])
            traj["sex"].append(int(a["sex"].iloc[0]))
            traj["t"].append([int(x) for x in a["t"]])
            traj["m"].append([int(round(x)) for x in a["mass"]])
        out["growth"][sp] = traj
    return out


def main(proc, out):
    os.makedirs(out, exist_ok=True)
    sites = pd.read_csv(f"{proc}/sites.csv")
    caps = pd.read_csv(f"{proc}/captures.csv", parse_dates=["date"], low_memory=False)
    caps = caps.dropna(subset=["date"])
    inds = pd.read_csv(f"{proc}/individuals.csv", parse_dates=["first", "last"])
    pop = pd.read_csv(f"{proc}/population_by_site.csv", parse_dates=["date"])
    main_s = pd.read_csv(f"{proc}/population_main_series.csv")
    rem = pd.read_csv(f"{proc}/removals.csv")
    nest = pd.read_csv(f"{proc}/nest_captures.csv", parse_dates=["date"])
    sess = pd.read_csv(f"{proc}/sessions.csv", parse_dates=["date"])

    meta = {"built": date.today().isoformat(), "epoch": "1970-01-01", "sites": [], "main": {}}

    for _, s in sites.iterrows():
        site = int(s["site"])
        c = caps[caps["site"] == site]
        g = c.dropna(subset=["grid_col"])
        # Grid extent from stations used 3+ times, so one-off typos in the
        # station column don't stretch the map.
        g = g[g.groupby("station")["station"].transform("size") >= 3]
        p = pop[pop["site"] == site]
        meta["sites"].append({
            "site": site, "name": s["name"], "habitat": s["habitat"], "area_ha": s["area_ha"],
            "start": s["start"], "end": s["end"], "treatment": s["treatment"],
            "location": s["location"],
            "n_individuals": int(c["ind"].nunique()), "n_captures": int(len(c)),
            "n_nest": int((nest["site"] == site).sum()),
            "n_nest_only": int(((inds["site"] == site) & (inds["n_captures"] == 0)).sum()),
            "grid": {"cols": int(g["grid_col"].max()) + 1 if len(g) else 0,
                     "rows": int(g["grid_row"].max()) + 1 if len(g) else 0,
                     "col0": int(g["grid_col"].min()) if len(g) else 0,
                     "row0": int(g["grid_row"].min()) if len(g) else 0},
            "pop": {sp: {"d": col(days(q["date"])), "mnka": col(q["mnka"]),
                         "m": col(q["mnka_m"]), "f": col(q["mnka_f"]), "caught": col(q["n_caught"])}
                    for sp, q in p.sort_values("date").groupby("species")},
            "removed": {sp: {"session": col(q["session"]), "n": col(q["n_removed"])}
                        for sp, q in rem[rem["site"] == site].groupby("species")},
        })

        # Per-site detail file
        ii = inds[inds["site"] == site].sort_values("first").reset_index(drop=True)
        index = {k: i for i, k in enumerate(ii["ind"])}
        c = c.sort_values(["date", "ampm"])
        ss = sess[sess["site"] == site].sort_values("date")
        detail = {
            "site": site,
            "sessions": {"session": col(ss["session"]), "d": col(days(ss["date"])),
                         "year": col(ss["year"]), "month": col(ss["month"])},
            "inds": {"id": ii["ind"].tolist(), "sp": ii["species"].tolist(),
                     "sex": col(ii["sex"], 0), "first": col(days(ii["first"])),
                     "last": col(days(ii["last"])), "n": col(ii["n_captures"]),
                     "nsess": col(ii["n_sessions"]), "nn": col(ii["n_nest"].fillna(0))},
            "caps": {"i": [index[k] for k in c["ind"]], "d": col(days(c["date"])),
                     "t": col(c["ampm"]), "s": col(c["session"]), "st": col(c["station"]),
                     "gc": col(c["grid_col"]), "gr": col(c["grid_row"]),
                     "m": col(c["mass"]), "rp": col(c["repro"])},
        }
        nn = nest[nest["site"] == site].sort_values("date")
        if len(nn):
            detail["nest"] = {"i": [index[k] for k in nn["ind"]], "d": col(days(nn["date"])),
                              "m": col(nn["mass"]), "rp": col(nn["repro"])}
        with open(f"{out}/site_{site:02d}.json", "w") as f:
            json.dump(detail, f, separators=(",", ":"))

    for (hab, sp), q in main_s.groupby(["habitat", "species"]):
        q = q.sort_values(["year", "month"])
        d = days(pd.to_datetime(dict(year=q["year"], month=q["month"], day=15)))
        meta["main"].setdefault(hab, {})[sp] = {
            "d": col(d), "density": [round(float(v), 1) for v in q["density_ha"]],
            "m": [round(float(v), 1) for v in q["density_m"]],
            "f": [round(float(v), 1) for v in q["density_f"]], "site": col(q["site"])}

    with open(f"{out}/meta.json", "w") as f:
        json.dump(meta, f, separators=(",", ":"))
    with open(f"{out}/mass.json", "w") as f:
        json.dump(mass_summaries(caps), f, separators=(",", ":"))

    total = sum(os.path.getsize(f"{out}/{n}") for n in os.listdir(out))
    print(f"wrote {len(os.listdir(out))} files, {total / 1e6:.1f} MB to {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
