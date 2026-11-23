# Getz vole trapping data: cleaning and analysis

Code that cleans and summarises Lowell L. Getz's long-term live-trapping study of prairie voles
(*Microtus ochrogaster*) and meadow voles (*M. pennsylvanicus*) in central Illinois, 1971–1997, and exports
the results for the interactive explorer on the project page:
**<https://calebvogt.github.io/projects/getz-vole-trapping/>**

The data are Getz's, published on his site at <https://www.life.illinois.edu/getz/>. The project page
lists the data sources, credits and his terms of use.

## Run

```bash
pip install -r requirements.txt

# 1. Mirror Getz's site (about 20 MB)
wget --mirror --no-parent -e robots=off --no-host-directories --cut-dirs=1 -P getz_mirror https://www.life.illinois.edu/getz/

# 2. Clean the capture records into tidy tables
python scripts/build_dataset.py getz_mirror processed

# 3. Summaries for the explorer (JSON), and the project tile
python scripts/export_web.py processed web_data
python scripts/social.py processed web_data
python scripts/traps.py processed web_data
python scripts/structure.py processed web_data
python scripts/make_tile.py processed tile.png
```

| Script | Writes |
|---|---|
| `build_dataset.py` | `captures.csv`, `nest_captures.csv`, `individuals.csv`, `sessions.csv`, `sites.csv`, `removals.csv`, `population_by_site.csv`, `population_main_series.csv` |
| `export_web.py` | `meta.json` (sites, population series), `mass.json` (body mass), `site_NN.json` (individuals and captures per site) |
| `social.py` | `social.json`, plus `social_animal_sessions.csv` in the processed folder |
| `traps.py` | `traps.json` |
| `structure.py` | `structure.json` |
| `make_tile.py` | the project tile: 25 years of vole density, one strip per habitat |

`load_raw.py`, `nest_data.py` and `sites.py` are helpers: the raw-file reader, the twice-weekly nest data,
and study-site metadata transcribed from Getz's pages.

## Cleaning decisions

- **Source**: the 38 "original field trap records" files (`MOnn`/`MPnn Final Cln.txt`): one row per capture
  with animal, sex, reproductive condition, mass, date, AM/PM check, station, species, site and session.
  Some files use old-Mac CR line endings.
- **Individuals**: species + site + animal number. Getz re-used number sequences, so a number unseen at
  a site for more than 300 days becomes a new individual (id suffix `.2`). The gap distribution has a clean
  break, and only one split occurs after date repair.
- **Dates**: a date more than 30 days from its site-session median (typos like `75/05/30` in a 1972 session)
  is replaced with the median and flagged `date_fixed`.
- **Twice-weekly nest data** (`nest_data.py`): Alf II/III, 1982–87, one row per animal per Mon–Tue or
  Thu–Fri session. Records link to a grid individual with the same number if within 60 days of its grid
  captures (2,077 of 3,114 Alf III numbers link; recorded sex agrees for 95%). Unlinked numbers become
  nest-only individuals (id suffix `n`). Reproductive codes are mapped 0 → non-reproductive,
  1/2 → reproductive.
- **Density**: minimum number known alive (MNKA) per site × species × session, divided by site area,
  overall and by sex (`mnka_m`, `mnka_f`; sex = each animal's modal recorded sex). The 25-year means
  reproduce Getz et al. (2001): bluegrass *M. ochrogaster* 19.6/ha (published 18), bluegrass
  *M. pennsylvanicus* 13.5 (14), tallgrass *M. ochrogaster* 8.0 (8).
- **Main series**: Getz's stitched site sequence per habitat (`sites.py: MAIN_SERIES`).
- **Removal sites** (MP09, MO10, MP14, MP18, MP19): counted per session in `removals.csv`, not treated as
  mark-recapture data.
- **Grid coordinates** (`sites.py: grid_xy`): exact for alfalfa (col × 100 + row), bluegrass 6/9/10/11
  (25-row serpentine) and tallgrass 15–19 (7 × 7 serpentine). Bluegrass 8's column length, bluegrass 7's
  three layouts and the L-shaped Trelease grid (tallgrass 12–14) are approximate.
- **AM/PM code 3** appears only in Alf II/III from 1982 to 1987, and is probably the 20:00 check of the
  social-behaviour study.
- **Body mass** (`export_web.py: mass_summaries`): monthly grid captures only, one value per animal per
  session (mean mass, highest reproductive code), masses 1–100 g; the ~20 male records coded pregnant are
  dropped.
- **Social structure** (`social.py`): a *co-capture* is two same-species animals of known sex in the same
  trap at the same check (Getz used multiple-capture traps); *trap sharing* is any common station within a
  session. Outputs the opposite-sex co-capture rate by density, two-animal trap composition against random
  expectation, partner counts for breeding animals, partners against male mass, and pair fidelity (pairs
  caught together again the next session, against neighbour pairs and any pair).
- **Trap return** (`traps.py`): distance to the previous capture (same session or a later one) and, for
  animals caught in 5 or more sessions, captures per session and same-trap return over their first five
  sessions. Exact-layout grids only (excludes 7 and 12–14).
- **Dataset structure** (`structure.py`): animals caught per grid per month (coverage), first-capture
  against caught-before counts per session, residence curves (share still known alive t months after first
  capture, for animals first caught at least 18 months before their grid's last session), share of new
  animals by month, and a year × month breeding calendar on the unmanipulated grids (females pregnant,
  males in breeding condition).

Not used: the home-range index files, the *Blarina* monthly densities, the raw-data compilations, and the
unpublished manuscripts and appendices.

## License

The code is MIT-licensed (`LICENSE`). The data are not part of this repository.
