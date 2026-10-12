"""Load the twice-weekly (3.5-day) prairie-vole nest-trapping data, 1982-87.

Source: https://www.life.illinois.edu/getz/3-5-data/
Alfalfa2_datafile.xls (Alf II, Mar 1982-Jul 1984) and Alfalfa3_datafile.xls
(Alf III, Oct 1983-May 1987). Traps were set around every known nest each
Mon-Tue and Thu-Fri; one row = one animal caught in one 2-day session, dated
by the session's first day. No station is recorded (nest trapping, not grid).

The same toe-clip numbers were used as in the monthly grid trapping of the
same fields (sites 2 and 3), so records are linked to grid individuals by
animal number when the dates overlap.
"""
import numpy as np
import pandas as pd

FILES = {2: "Alfalfa2_datafile.xls", 3: "Alfalfa3_datafile.xls"}
MONTHS = {m: i + 1 for i, m in enumerate(
    ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"])}
# A nest record links to a grid individual if it falls within this many days
# of that individual's first-last grid capture window.
LINK_WINDOW_DAYS = 60
ID_REUSE_GAP_DAYS = 300


def load_nest(mirror_dir):
    parts = []
    for site, name in FILES.items():
        df = pd.read_excel(f"{mirror_dir}/3-5-data/{name}", engine="xlrd", dtype=str)
        df.columns = ["disk", "month_index", "day", "mon", "yr", "animal", "sex", "repro", "mass", "new"]
        df["site"] = site
        parts.append(df)
    df = pd.concat(parts, ignore_index=True)
    df["animal"] = pd.to_numeric(df["animal"], errors="coerce")
    df["date"] = pd.to_datetime(dict(
        year=pd.to_numeric(df["yr"], errors="coerce") + 1900,
        month=df["mon"].str.strip().str.upper().map(MONTHS),
        day=pd.to_numeric(df["day"], errors="coerce")), errors="coerce")
    df["sex"] = df["sex"].str.strip().map({"M": 1, "F": 2}).fillna(0).astype(int)
    # 3.5-day coding: 0 non-reproductive, 1 reproductive, 2 "also reproductive"
    # (per Getz). Map onto the grid coding: 1 non-repro, 2 repro.
    r = pd.to_numeric(df["repro"], errors="coerce")
    df["repro"] = np.select([r == 0, r.isin([1, 2])], [1, 2], -1)
    m = pd.to_numeric(df["mass"], errors="coerce")
    df["mass"] = m.where(m > 0)
    df = df.dropna(subset=["animal", "date"])
    df = df[df["animal"] > 0]
    df["animal"] = df["animal"].astype(int)
    # Same animal entered twice for one session (different data disks): keep one.
    df = (df.sort_values("mass", na_position="last")
            .drop_duplicates(["site", "animal", "date"])
            .sort_values(["site", "date", "animal"]))
    return df[["site", "date", "animal", "sex", "repro", "mass"]].reset_index(drop=True)


def link_nest(nest, individuals):
    """Assign each nest record an individual id.

    Linked to a grid individual (prairie vole, same site, same number) whose
    grid window +/- LINK_WINDOW_DAYS contains the record; otherwise grouped
    into a nest-only individual (split on gaps > ID_REUSE_GAP_DAYS).
    """
    g = individuals[individuals["species"] == "MO"]
    win = pd.Timedelta(days=LINK_WINDOW_DAYS)
    out = []
    for (site, animal), rec in nest.groupby(["site", "animal"]):
        cand = g[(g["site"] == site) & (g["animal"] == animal)]
        ids = pd.Series(pd.NA, index=rec.index, dtype="object")
        for _, c in cand.iterrows():
            hit = (rec["date"] >= c["first"] - win) & (rec["date"] <= c["last"] + win) & ids.isna()
            ids[hit] = c["ind"]
        rest = rec[ids.isna()]
        if len(rest):
            gap = rest["date"].diff().dt.days > ID_REUSE_GAP_DAYS
            epoch = gap.cumsum()
            base = f"MO{site:02d}-{animal}"
            taken = set(cand["ind"])
            for e, grp in rest.groupby(epoch):
                nid = base + "n" + ("" if e == 0 else f".{e + 1}")
                while nid in taken:
                    nid += "'"
                ids[grp.index] = nid
        out.append(rec.assign(ind=ids, linked=ids.isin(set(cand["ind"]))))
    return pd.concat(out).sort_values(["site", "date", "ind"]).reset_index(drop=True)
