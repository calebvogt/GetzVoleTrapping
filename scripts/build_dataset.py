"""Build tidy tables from the mirrored Getz website.

    python scripts/build_dataset.py <mirror_dir> <out_dir>

Writes captures.csv, nest_captures.csv, individuals.csv, sessions.csv,
sites.csv, removals.csv and the population tables to <out_dir>.
"""
import os
import sys

import numpy as np
import pandas as pd

from load_raw import load_all
from nest_data import link_nest, load_nest
from sites import HABITAT, MAIN_SERIES, REMOVAL_FILES, SITES, grid_xy

# Gaps between consecutive captures of the same number at the same site longer
# than this are treated as a re-used number (a different animal). The gap
# distribution has a clean break: ~60 gaps exceed 180 d, only ~16 exceed 300 d.
ID_REUSE_GAP_DAYS = 300
# A capture date more than this far from its site-session median is treated as
# a data-entry error and replaced with the median.
DATE_TOLERANCE_DAYS = 30

SPECIES_NAME = {"MO": "M. ochrogaster", "MP": "M. pennsylvanicus"}


def to_num(s):
    return pd.to_numeric(s.str.strip(), errors="coerce")


def clean_captures(df):
    for c in ["animal", "sex", "repro", "mass", "ampm", "station", "new_recap", "session"]:
        df[c] = to_num(df[c])
    df["species"] = df.pop("file_species")
    df["site"] = df.pop("file_site")
    df = df.drop(columns=["new_recap"])

    d = pd.to_datetime(df["date"].str.strip().str.rstrip("`"), format="%y/%m/%d", errors="coerce")
    d = d.where(d.dt.year < 2000, d - pd.DateOffset(years=100))
    df["date"] = d

    # Session-level date: median date of all rows for that site-session (both
    # species), used to repair bad dates and to date zero-capture sessions.
    med = df.groupby(["site", "session"])["date"].median().dt.floor("D").rename("session_date")
    df = df.join(med, on=["site", "session"])
    off = (df["date"] - df["session_date"]).dt.days.abs()
    df["date_fixed"] = df["date"].isna() | (off > DATE_TOLERANCE_DAYS)
    df.loc[df["date_fixed"], "date"] = df.loc[df["date_fixed"], "session_date"]
    return df


def sessions_table(df, months):
    """One row per site x session that was trapped (from either species file)."""
    s = (df.groupby(["site", "session"])
           .agg(date=("session_date", "first"))
           .reset_index()
           .dropna(subset=["date"]))
    s = s.merge(months, on="session", how="left")
    # Fall back to the calendar month of the session date when the
    # session->month table has no entry.
    s["year"] = s["year"].fillna(s["date"].dt.year).astype(int)
    s["month"] = s["month"].fillna(s["date"].dt.month).astype(int)
    return s.sort_values(["site", "date"]).reset_index(drop=True)


def assign_individuals(cap):
    """Give each marked animal a stable id, splitting re-used numbers."""
    cap = cap.sort_values(["species", "site", "animal", "date"]).copy()
    gap = cap.groupby(["species", "site", "animal"])["date"].diff().dt.days
    epoch = (gap > ID_REUSE_GAP_DAYS).groupby([cap["species"], cap["site"], cap["animal"]]).cumsum()
    base = cap["species"] + cap["site"].map("{:02d}".format) + "-" + cap["animal"].astype(int).astype(str)
    cap["ind"] = base + np.where(epoch > 0, "." + (epoch + 1).astype(str), "")
    return cap


def individuals_table(cap):
    def mode_pos(s):
        s = s[s > 0]
        return int(s.mode().iloc[0]) if len(s) else 0

    g = cap.groupby("ind")
    ind = g.agg(species=("species", "first"), site=("site", "first"), animal=("animal", "first"),
                first=("date", "min"), last=("date", "max"), n_captures=("date", "size"),
                n_sessions=("session", "nunique"), n_stations=("station", "nunique"),
                max_mass=("mass", "max"), ever_pregnant=("repro", lambda r: bool((r == 3).any())))
    ind["sex"] = g["sex"].agg(mode_pos)
    ind["sex_conflict"] = g["sex"].agg(lambda s: s[s > 0].nunique() > 1)
    ind["first_mass"] = g["mass"].first()
    ind["span_days"] = (ind["last"] - ind["first"]).dt.days
    ind["habitat"] = ind["site"].map(HABITAT)
    return ind.reset_index()


def add_nest(ind, nest):
    """Fold twice-weekly nest records into the individuals table.

    Adds n_nest / n_nest_sessions, extends first/last to cover nest records,
    and appends nest-only individuals (never caught on the grid).
    """
    g = nest.groupby("ind")
    agg = g.agg(nfirst=("date", "min"), nlast=("date", "max"), n_nest=("date", "size"),
                nsex=("sex", lambda s: int(s[s > 0].mode().iloc[0]) if (s > 0).any() else 0),
                nmass=("mass", "max"), site=("site", "first"), animal=("animal", "first"))
    ind = ind.set_index("ind")
    ind["n_nest"] = agg["n_nest"].reindex(ind.index).fillna(0).astype(int)
    both = ind.index.intersection(agg.index)
    ind.loc[both, "first"] = np.minimum(ind.loc[both, "first"], agg.loc[both, "nfirst"])
    ind.loc[both, "last"] = np.maximum(ind.loc[both, "last"], agg.loc[both, "nlast"])
    new = agg.loc[agg.index.difference(ind.index)]
    add = pd.DataFrame({
        "species": "MO", "site": new["site"], "animal": new["animal"], "first": new["nfirst"],
        "last": new["nlast"], "n_captures": 0, "n_sessions": 0, "n_stations": 0,
        "max_mass": new["nmass"], "ever_pregnant": False, "sex": new["nsex"], "sex_conflict": False,
        "first_mass": np.nan, "habitat": "alfalfa", "n_nest": new["n_nest"]}, index=new.index)
    ind = pd.concat([ind, add])
    ind["n_captures"] = ind["n_captures"].astype(int)
    ind["n_sessions"] = ind["n_sessions"].astype(int)
    ind["span_days"] = (ind["last"] - ind["first"]).dt.days
    ind.index.name = "ind"
    return ind.reset_index()


def mnka(cap, sess, sex_of):
    """Minimum number known alive per site x species x session, overall and by sex.

    An animal is known alive in every trapped session from its first through
    its last capture (Krebs 1966), the estimator Getz used. Sex is each
    individual's modal recorded sex (`sex_of`: ind -> 1 male, 2 female, 0 unknown).
    """
    rows = []
    for (sp, site), c in cap.groupby(["species", "site"]):
        ss = sess[sess["site"] == site].reset_index(drop=True)
        order = {s: i for i, s in enumerate(ss["session"])}
        idx = c["session"].map(order)
        span = idx.groupby(c["ind"]).agg(["min", "max"]).dropna().astype(int)
        sex = span.index.map(sex_of)
        known = {k: np.zeros(len(ss), int) for k in ("all", 1, 2)}
        for (lo, hi), sx in zip(span.itertuples(index=False), sex):
            known["all"][lo:hi + 1] += 1
            if sx in (1, 2):
                known[sx][lo:hi + 1] += 1
        caught = c.groupby("session")["ind"].nunique().reindex(ss["session"]).fillna(0).astype(int)
        rows.append(pd.DataFrame({"species": sp, "site": site, "session": ss["session"],
                                  "date": ss["date"], "year": ss["year"], "month": ss["month"],
                                  "n_caught": caught.values, "mnka": known["all"],
                                  "mnka_m": known[1], "mnka_f": known[2]}))
    out = pd.concat(rows, ignore_index=True)
    area = out["site"].map(lambda s: SITES[s]["area_ha"])
    out["area_ha"] = area
    out["density_ha"] = out["mnka"] / area
    out["density_m"] = out["mnka_m"] / area
    out["density_f"] = out["mnka_f"] / area
    return out


def main_series(pop):
    """Stitch Getz's unmanipulated site sequence into one monthly series per habitat."""
    parts = []
    for hab, seq in MAIN_SERIES.items():
        for site, lo, hi in seq:
            p = pop[pop["site"] == site].copy()
            ym = p["year"].astype(str) + "-" + p["month"].map("{:02d}".format)
            if lo:
                p = p[ym >= lo]
            if hi:
                p = p[ym.loc[p.index] <= hi]
            p["habitat"] = hab
            parts.append(p)
    m = pd.concat(parts)
    return (m.groupby(["habitat", "species", "year", "month"])
             .agg(density_ha=("density_ha", "mean"), density_m=("density_m", "mean"),
                  density_f=("density_f", "mean"), mnka=("mnka", "mean"), site=("site", "first"))
             .reset_index())


def main(mirror, out):
    os.makedirs(out, exist_ok=True)
    raw = load_all(mirror)
    df = clean_captures(raw)

    months = pd.read_excel(f"{mirror}/25-year/trapping-montha/trapping-months.xls", engine="xlrd")
    months.columns = ["year", "month", "session"]
    # 1975-76 was trapped at 3-week intervals, so some months list two
    # sessions ("45-46"); both map to that month.
    months["session"] = months["session"].astype(str).str.split("-")
    months = months.explode("session")
    months["session"] = pd.to_numeric(months["session"], errors="coerce")
    months = months.dropna().astype(int)
    months["year"] += 1900
    months = months.drop_duplicates("session")

    sess = sessions_table(df, months)

    removal = [(sp, st) in REMOVAL_FILES for sp, st in zip(df["species"], df["site"])]
    df["removal"] = removal
    marked = df[(df["animal"] > 0) & ~df["removal"] & df["date"].notna()].copy()
    marked = assign_individuals(marked)
    xy = [grid_xy(s, st) for s, st in zip(marked["site"], marked["station"])]
    marked["grid_col"] = [p[0] if p else np.nan for p in xy]
    marked["grid_row"] = [p[1] if p else np.nan for p in xy]
    marked["habitat"] = marked["site"].map(HABITAT)

    ind = individuals_table(marked)
    nest = link_nest(load_nest(mirror), ind)
    ind = add_nest(ind, nest)
    pop = mnka(marked, sess, ind.set_index("ind")["sex"].to_dict())
    main = main_series(pop)

    # Removal files: count animals removed per session (marked or not).
    rem = df[df["removal"] & (df["animal"].fillna(0) != -1)]
    rem = (rem.groupby(["species", "site", "session"]).size().rename("n_removed").reset_index())

    sites = pd.DataFrame([{"site": k, "habitat": HABITAT[k], **v} for k, v in SITES.items()])

    cols = ["ind", "species", "site", "habitat", "animal", "session", "date", "date_fixed", "ampm",
            "station", "grid_col", "grid_row", "sex", "repro", "mass", "file"]
    marked[cols].sort_values(["date", "site", "ind"]).to_csv(f"{out}/captures.csv", index=False)
    ind.to_csv(f"{out}/individuals.csv", index=False)
    nest.to_csv(f"{out}/nest_captures.csv", index=False)
    pop.to_csv(f"{out}/population_by_site.csv", index=False)
    main.to_csv(f"{out}/population_main_series.csv", index=False)
    sess.to_csv(f"{out}/sessions.csv", index=False)
    rem.to_csv(f"{out}/removals.csv", index=False)
    sites.to_csv(f"{out}/sites.csv", index=False)

    print(f"raw rows            {len(raw):>7}")
    print(f"marked captures     {len(marked):>7}  (dates repaired: {int(marked['date_fixed'].sum())})")
    print(f"individuals         {len(ind):>7}  (id re-use splits: {int(ind['ind'].str.contains('[.]').sum())})")
    print(f"nest records        {len(nest):>7}  (linked to grid animals: {int(nest['linked'].sum())}, "
          f"nest-only animals: {int((ind['n_captures'] == 0).sum())})")
    print(f"site-sessions       {len(sess):>7}")
    print(f"removal rows        {int(rem['n_removed'].sum()):>7}")
    print(f"stations w/o grid   {int(marked['grid_col'].isna().sum()):>7}")
    return marked, ind, pop, main, sess, nest


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
