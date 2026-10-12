"""Social-structure summaries from same-trap co-captures and trap sharing.

Getz used multiple-capture live-traps, so two animals found in the same trap
at the same check were physically together. From the monthly grid captures:

  co-capture    same grid, date, check (AM/PM/Eve) and station
  trap sharing  two animals caught at >=1 common station in the same session
                (any check): a coarse within-session range overlap

Only same-species associations among animals of known sex are counted.
Density is the species' MNKA per hectare on that grid in that session.

    python scripts/social.py <processed_dir> <web_data_dir>
"""
import itertools
import json
import sys

import numpy as np
import pandas as pd

HABS = ["all", "alfalfa", "bluegrass", "tallgrass"]
DENSITY_BINS = [0, 10, 20, 40, 80, 160, 1000]
MASS_BINS = [25, 30, 35, 40, 45, 50]          # lower edges; < 25 g pooled in, >= 50 g pooled in
MIN_N = 20


def load(proc):
    c = pd.read_csv(f"{proc}/captures.csv", parse_dates=["date"], low_memory=False)
    ind = pd.read_csv(f"{proc}/individuals.csv")
    pop = pd.read_csv(f"{proc}/population_by_site.csv")
    sess = pd.read_csv(f"{proc}/sessions.csv", parse_dates=["date"])
    c["isex"] = c["ind"].map(ind.set_index("ind")["sex"])
    c = c[(c["station"] > 0) & (c["isex"].isin([1, 2]))].copy()
    c["isex"] = c["isex"].astype(int)
    c["key"] = (c["site"].astype(str) + "|" + c["date"].dt.strftime("%Y%m%d") + "|"
                + c["ampm"].fillna(-1).astype(int).astype(str) + "|" + c["station"].astype(int).astype(str))
    c["has_check"] = c["ampm"].isin([1, 2, 3])
    return c, pop, sess


def dyads(c):
    """All same-species male-female pairs caught in the same trap at the same check."""
    cc = c[c["has_check"]].drop_duplicates(["key", "ind"])
    rows = []
    for (sp, key), g in cc.groupby(["species", "key"]):
        if len(g) < 2:
            continue
        recs = g[["ind", "isex", "site", "session", "habitat"]].to_numpy()
        for a, b in itertools.combinations(recs, 2):
            rows.append((sp, a[4], a[2], a[3], key, a[0], a[1], b[0], b[1]))
    return pd.DataFrame(rows, columns=["species", "habitat", "site", "session", "key", "a", "sa", "b", "sb"])


def animal_sessions(c, d, pop):
    """One row per animal per session with its co-capture and trap-sharing partner counts."""
    a = (c.groupby(["species", "habitat", "site", "session", "ind"])
          .agg(sex=("isex", "first"), repro=("repro", "max"), mass=("mass", "mean"),
               ncap=("key", "size"), st=("station", lambda s: frozenset(s.astype(int))))
          .reset_index())
    # co-capture partners
    both = pd.concat([d.rename(columns={"a": "ind", "sa": "s1", "b": "o", "sb": "s2"}),
                      d.rename(columns={"b": "ind", "sb": "s1", "a": "o", "sa": "s2"})])
    opp = both[both["s1"] != both["s2"]].groupby(["site", "session", "ind"])["o"].nunique().rename("co_opp")
    anyp = both.groupby(["site", "session", "ind"])["o"].nunique().rename("co_any")
    a = a.join(opp, on=["site", "session", "ind"]).join(anyp, on=["site", "session", "ind"])
    a[["co_opp", "co_any"]] = a[["co_opp", "co_any"]].fillna(0).astype(int)
    # trap sharing (common station within the session)
    a["ov_opp"] = 0
    a["ov_same"] = 0
    for _, g in a.groupby(["species", "site", "session"]):
        sts, sx = list(g["st"]), list(g["sex"])
        opp, same = [], []
        for i in range(len(g)):
            o = s = 0
            for j in range(len(g)):
                if i != j and sts[i] & sts[j]:
                    if sx[i] != sx[j]:
                        o += 1
                    else:
                        s += 1
            opp.append(o)
            same.append(s)
        a.loc[g.index, "ov_opp"] = opp
        a.loc[g.index, "ov_same"] = same
    a = a.merge(pop[["species", "site", "session", "density_ha"]], on=["species", "site", "session"], how="left")
    a["breeding"] = a["repro"].isin([2, 3])
    return a.drop(columns="st")


def fidelity(d, c, sess):
    """Do male-female pairs caught together stay together the next session?

    For each MF pair co-captured in session t where both animals are caught
    again in the grid's next trapped session, the share co-captured again.
    Two comparisons, computed the same way:
      neighbors   MF pairs that shared a trap in session t (same station, any
                  check) but were never caught in it together
      any pair    every MF pair caught on the grid in session t
    """
    nxt = {}
    for site, s in sess.sort_values("date").groupby("site"):
        ss = list(s["session"])
        nxt.update({(site, a): b for a, b in zip(ss, ss[1:])})
    caught = c.groupby(["species", "site", "session"])["ind"].apply(set).to_dict()
    mf = d[d["sa"] != d["sb"]].copy()
    mf["m"] = np.where(mf["sa"] == 1, mf["a"], mf["b"])
    mf["f"] = np.where(mf["sa"] == 1, mf["b"], mf["a"])
    together = set(zip(mf["species"], mf["site"], mf["session"], mf["m"], mf["f"]))
    st = c.groupby(["species", "habitat", "site", "session", "ind"]).agg(
        sex=("isex", "first"), st=("station", lambda s: frozenset(s.astype(int)))).reset_index()

    def rate(pairs):
        followed = again = 0
        for sp, site, ses, m, f in pairs:
            n = nxt.get((site, ses))
            if n is None:
                continue
            here = caught.get((sp, site, n), set())
            if m in here and f in here:
                followed += 1
                again += (sp, site, n, m, f) in together
        return followed, again

    out = {}
    for hab in HABS:
        for sp in ["MO", "MP"]:
            q = st[(st["species"] == sp) & ((st["habitat"] == hab) | (hab == "all"))]
            pairs_t, neigh, anyp = [], [], []
            for (site, ses), g in q.groupby(["site", "session"]):
                ms, fs = g[g["sex"] == 1], g[g["sex"] == 2]
                for m, ms_st in zip(ms["ind"], ms["st"]):
                    for f, fs_st in zip(fs["ind"], fs["st"]):
                        key = (sp, site, ses, m, f)
                        anyp.append(key)
                        if key in together:
                            pairs_t.append(key)
                        elif ms_st & fs_st:
                            neigh.append(key)
            res = {}
            for name, pairs in [("together", pairs_t), ("neighbors", neigh), ("any", anyp)]:
                fo, ag = rate(pairs)
                res[name] = {"followed": fo, "again": ag, "p": round(ag / fo, 4) if fo else None}
            out[f"{sp}|{hab}"] = res
    return out


def summarize(a, d, c, sess):
    out = {"density_bins": DENSITY_BINS, "mass_bins": MASS_BINS, "comp": {}, "together": {},
           "partners": {}, "mass": {}, "fidelity": fidelity(d, c, sess)}
    a["dbin"] = pd.cut(a["density_ha"], DENSITY_BINS, labels=False, right=True)
    a["mbin"] = np.clip(np.floor(a["mass"] / 5) * 5, MASS_BINS[0], MASS_BINS[-1])
    for hab in HABS:
        for sp in ["MO", "MP"]:
            # composition of two-animal co-captures
            dd = d[(d["species"] == sp) & ((d["habitat"] == hab) | (hab == "all"))]
            size = dd.groupby("key").size()
            two = dd[dd["key"].isin(size[size == 1].index)]   # exactly one dyad == exactly two animals
            nm = (two["sa"] == 1).astype(int) + (two["sb"] == 1).astype(int)
            n = len(two)
            if n:
                pm = nm.sum() / (2 * n)
                out["comp"][f"{sp}|{hab}"] = {"n": int(n), "MF": int((nm == 1).sum()), "MM": int((nm == 2).sum()),
                                              "FF": int((nm == 0).sum()), "expMF": round(2 * pm * (1 - pm), 3)}
            q = a[(a["species"] == sp) & ((a["habitat"] == hab) | (hab == "all"))]
            for sx in [1, 2]:
                k = f"{sp}|{hab}|{sx}"
                s = q[q["sex"] == sx]
                rows = []
                for b in range(len(DENSITY_BINS) - 1):
                    v = s[s["dbin"] == b]
                    rows.append({"n": int(len(v)),
                                 "opp": round(float((v["co_opp"] > 0).mean()), 3) if len(v) >= MIN_N else None,
                                 "any": round(float((v["co_any"] > 0).mean()), 3) if len(v) >= MIN_N else None})
                out["together"][k] = rows
                br = s[s["breeding"]]
                if len(br) >= MIN_N:
                    out["partners"][k] = {
                        "n": int(len(br)),
                        "co": [int((br["co_opp"].clip(upper=3) == i).sum()) for i in range(4)],
                        "ov": [int((br["ov_opp"].clip(upper=3) == i).sum()) for i in range(4)],
                        "ov_same": round(float((br["ov_same"] > 0).mean()), 3)}
            males = q[(q["sex"] == 1) & q["breeding"] & (q["mass"] > 0)]
            rows = []
            for mb in MASS_BINS:
                v = males[males["mbin"] == mb]
                ok = len(v) >= MIN_N
                rows.append({"n": int(len(v)),
                             "co1": round(float((v["co_opp"] >= 1).mean()), 3) if ok else None,
                             "co2": round(float((v["co_opp"] >= 2).mean()), 3) if ok else None,
                             "ov1": round(float((v["ov_opp"] >= 1).mean()), 3) if ok else None,
                             "ov2": round(float((v["ov_opp"] >= 2).mean()), 3) if ok else None,
                             "ncap": round(float(v["ncap"].mean()), 2) if ok else None})
            out["mass"][f"{sp}|{hab}"] = rows
    return out


def main(proc, out):
    c, pop, sess = load(proc)
    d = dyads(c)
    a = animal_sessions(c, d, pop)
    a.to_csv(f"{proc}/social_animal_sessions.csv", index=False)
    s = summarize(a, d, c, sess)
    with open(f"{out}/social.json", "w") as f:
        json.dump(s, f, separators=(",", ":"))
    for k in ["MO|all", "MP|all"]:
        print(k, "composition", s["comp"][k], "fidelity", s["fidelity"][k])
    return s


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
