"""Study-site metadata for the Getz 25-year vole study.

Transcribed from the Getz site pages (study_sites/list.html, size.html,
sequence.html, grid-maps/index.html, microtus/manipulative/index.html).
Where those pages disagree, the choice made is noted inline.
"""

HABITAT = {**{s: "alfalfa" for s in range(1, 6)},
           **{s: "bluegrass" for s in range(6, 12)},
           **{s: "tallgrass" for s in range(12, 20)}}

# treatment codes:
#   control      both species present, unmanipulated
#   supp         supplementally fed (Purina rabbit chow at every station)
#   mp_removed   all M. pennsylvanicus removed monthly (M. ochrogaster alone)
#   mo_removed   all M. ochrogaster removed monthly (M. pennsylvanicus alone)
SITES = {
    1:  dict(name="Alf I",   area_ha=1.4,  start="1972-05", end="1977-02", treatment="control",
             location="Phillips Tract"),
    2:  dict(name="Alf II",  area_ha=1.0,  start="1977-03", end="1984-07", treatment="control",
             location="Phillips Tract"),
    3:  dict(name="Alf III", area_ha=1.0,  start="1983-09", end="1989-12", treatment="control",
             location="Phillips Tract"),
    4:  dict(name="Alf IV",  area_ha=1.0,  start="1990-01", end="1993-05", treatment="control",
             location="Phillips Tract"),
    5:  dict(name="Alf V",   area_ha=1.0,  start="1993-06", end="1997-05", treatment="control",
             location="Phillips Tract"),
    6:  dict(name="BG Phillips N+S", area_ha=2.0, start="1971-08", end="1977-03", treatment="control",
             location="Phillips Tract"),
    7:  dict(name="BG Control", area_ha=0.8, start="1977-03", end="1987-07", treatment="control",
             location="Phillips Tract"),
    8:  dict(name="BG Supplemental", area_ha=0.5, start="1977-06", end="1983-11", treatment="supp",
             location="Phillips Tract"),
    # list.html + manipulative page agree site 9 is the M. ochrogaster-only
    # (Mp removed) site; grid-maps page says the reverse. Majority wins.
    9:  dict(name="BG North (Mp removed)", area_ha=1.0, start="1977-05", end="1997-05",
             treatment="mp_removed", location="Phillips Tract"),
    10: dict(name="BG South (Mo removed)", area_ha=1.0, start="1977-05", end="1987-05",
             treatment="mo_removed", location="Phillips Tract"),
    11: dict(name="BG South", area_ha=1.0, start="1987-06", end="1997-05", treatment="control",
             location="Phillips Tract"),
    12: dict(name="TG Trelease (72-77)", area_ha=2.0, start="1972-03", end="1977-03",
             treatment="control", location="Trelease Prairie"),
    13: dict(name="TG Trelease (84-87)", area_ha=0.7, start="1984-09", end="1987-04",
             treatment="control", location="Trelease Prairie"),
    14: dict(name="TG Trelease (Mp removed)", area_ha=0.67, start="1981-01", end="1983-11",
             treatment="mp_removed", location="Trelease Prairie"),
    15: dict(name="TG Control", area_ha=0.5, start="1977-09", end="1984-05", treatment="control",
             location="Phillips Tract"),
    16: dict(name="TG Phillips Tract", area_ha=0.5, start="1987-09", end="1997-05", treatment="control",
             location="Phillips Tract"),
    17: dict(name="TG Supplemental", area_ha=0.5, start="1977-09", end="1984-05", treatment="supp",
             location="Phillips Tract"),
    18: dict(name="TG Supplemental (Mp removed)", area_ha=0.5, start="1984-09", end="1987-03",
             treatment="supp_mp_removed", location="Phillips Tract"),
    19: dict(name="TG West (Mp removed)", area_ha=0.5, start="1984-09", end="1997-05",
             treatment="mp_removed", location="Phillips Tract"),
}

# Species-site files that are removal records (animals captured and taken off
# the site), not mark-recapture data.
REMOVAL_FILES = {("MP", 9), ("MO", 10), ("MP", 14), ("MP", 18), ("MP", 19)}

# Getz's "sequence of study sites for data analysis of main data files"
# (sequence.html): the unmanipulated series stitched into one 25-yr record per
# habitat. Alf 2 is used through Sep 1983, Alf 3 from Oct 1983.
MAIN_SERIES = {
    "alfalfa":   [(1, None, None), (2, None, "1983-09"), (3, "1983-10", None), (4, None, None),
                  (5, None, None)],
    "bluegrass": [(6, None, None), (7, None, None), (11, None, None)],
    "tallgrass": [(12, None, None), (15, None, None), (13, None, None), (16, None, None)],
}


def grid_xy(site, station, year=None):
    """Map a trap-station code to (col, row) grid coordinates, or None.

    Encodings from grid-maps/index.html and the scanned grid PDFs:
      Alf 1-5       col*100 + row (rows 01-14)
      BG 6,9,10,11  serpentine 1-225, 25 stations per column
      BG 7          1977-Oct 1980: serpentine 1-72 (6 per column assumed);
                    after: col*100 + row on a 5-m (later 10-m) grid
      BG 8          serpentine 1-60, 10 per column (column length assumed)
      TG 15-19      serpentine 1-49, 7 per column
      TG 12-14      irregular L-shaped serpentine; approximated as 25 per column
    """
    if station is None or not station > 0:
        return None
    s = int(station)
    if site <= 5:
        c, r = divmod(s, 100)
        return (c, r) if 1 <= r <= 14 else None
    if site in (6, 9, 10, 11):
        return _serp(s, 25) if s <= 225 else None
    if site == 7:
        if s <= 72:
            return _serp(s, 6)
        c, r = divmod(s, 100)
        return (c, r) if 1 <= c <= 17 and 1 <= r <= 17 else None
    if site == 8:
        return _serp(s, 10) if s <= 60 else None
    if 15 <= site <= 19:
        return _serp(s, 7) if s <= 49 else None
    if site in (12, 13, 14):
        return _serp(s, 25) if s <= 300 else None
    return None


def _serp(s, n):
    c, i = divmod(s - 1, n)
    r = i if c % 2 == 0 else n - 1 - i
    return (c, r)
