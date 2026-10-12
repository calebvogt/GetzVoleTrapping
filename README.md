# GetzVoleTrapping

Cleaning pipeline and web export for Lowell L. Getz's long-term vole live-trapping data
(University of Illinois, 1971–1997), which he released at <https://www.life.illinois.edu/getz/>.
The output powers the interactive explorer at
`calebvogt.github.io/projects/getz-vole-trapping/`.

All data credit belongs to Getz and colleagues. Archived versions are on EDI:
[25-year study](https://doi.org/10.6073/pasta/62617e944d053ba553d37cc4568c2367),
[twice-weekly alfalfa study](https://doi.org/10.6073/pasta/e62ec4197e8e2a687f2ef9d2c13cc3c0),
[*Blarina* densities](https://doi.org/10.6073/pasta/920c200ed6b3dd3dea83e8a62e5dff44).
Read Getz's usage request on his site before publishing analyses.

## Layout

| Where | What |
|---|---|
| OneDrive `2026_GetzVoleTrapping/Data/getz_website_mirror/` | Full `wget` mirror of the Getz site (167 files, untouched) |
| OneDrive `2026_GetzVoleTrapping/Data/processed/` | Tidy CSVs produced by `build_dataset.py` |
| `calebvogt.github.io/projects/getz-vole-trapping/data/` | Compact JSON for the explorer (`export_web.py`) |
| `calebvogt.github.io/_projects/2026-getz-vole-trapping.html` | The explorer page (JS/CSS live next to the data) |
| `scripts/` | The pipeline (this repo) |

## Run

```bash
python3 -m venv .venv && .venv/bin/pip install pandas xlrd matplotlib
D="$HOME/Library/CloudStorage/OneDrive-UCB-O365/Work/Projects/2026_GetzVoleTrapping/Data"
.venv/bin/python scripts/build_dataset.py "$D/getz_website_mirror" "$D/processed"
.venv/bin/python scripts/export_web.py "$D/processed" ../calebvogt.github.io/projects/getz-vole-trapping/data
.venv/bin/python scripts/social.py "$D/processed" ../calebvogt.github.io/projects/getz-vole-trapping/data
.venv/bin/python scripts/traps.py "$D/processed" ../calebvogt.github.io/projects/getz-vole-trapping/data
.venv/bin/python scripts/make_tile.py "$D/processed" ../calebvogt.github.io/projects/getz-vole-trapping/tile.png
```

To refresh the mirror: `wget --mirror --no-parent -e robots=off --no-host-directories --cut-dirs=1 https://www.life.illinois.edu/getz/`

## Cleaning decisions

- **Source**: the 38 "original field trap records" files (`MOnn`/`MPnn Final Cln.txt`): one row per capture
  with animal, sex, reproductive condition, mass, date, AM/PM check, station, species, site, session.
  Some files use old-Mac CR line endings.
- **Individuals**: species + site + animal number. Getz re-used number sequences, so a number unseen at
  a site for > 300 days becomes a new individual (id suffix `.2`). The gap distribution has a clean break, and only
  one split occurs after date repair.
- **Dates**: a date > 30 days from its site-session median (typos like `75/05/30` in a 1972 session) is
  replaced with the median and flagged `date_fixed`.
- **Twice-weekly nest data** (`nest_data.py`): Alf II/III, 1982–87, one row per animal per Mon–Tue or
  Thu–Fri session. Records link to a grid individual with the same number if within 60 days of its grid
  captures (2,077 of 3,114 Alf III numbers link; recorded sex agrees for 95%). Unlinked numbers become
  nest-only individuals (id suffix `n`). Repro codes are mapped 0 → non-repro, 1/2 → reproductive.
- **Density**: minimum number known alive (MNKA) per site × species × session, divided by site area,
  overall and by sex (`mnka_m`, `mnka_f`; sex = each animal's modal recorded sex).
  The 25-year means reproduce Getz et al. (2001): bluegrass *M. ochrogaster* 19.6/ha (published 18),
  bluegrass *M. pennsylvanicus* 13.5 (14), tallgrass *M. ochrogaster* 8.0 (8).
- **Main series**: Getz's stitched site sequence per habitat (`sites.py: MAIN_SERIES`).
- **Removal sites** (MP09, MO10, MP14, MP18, MP19): counted per session in `removals.csv`, not treated as
  mark-recapture data.
- **Grid coordinates** (`sites.py: grid_xy`): exact for alfalfa (col×100+row), BG 6/9/10/11 (25-row
  serpentine), and TG 15–19 (7×7 serpentine). BG 8 column length, BG 7's three layouts, and the L-shaped
  Trelease grid (TG 12–14) are approximate.
- **AM/PM code 3** appears only in Alf II/III from 1982 to 1987, and is probably the 20:00 check from the
  social-behavior study.

- **Body mass** (`export_web.py: mass_summaries`): monthly grid captures only, one value per animal per
  session (mean mass, highest repro code), masses 1–100 g; the ~20 male records coded pregnant are dropped.

- **Social structure** (`social.py`): a *co-capture* is two same-species animals of known sex in the same
  trap at the same check (multiple-capture traps); *trap sharing* is any common station within a session.
  Outputs opposite-sex co-capture rate by density, two-animal trap composition vs. random expectation,
  partner counts for breeding animals, partners vs. male mass, and pair fidelity (pairs caught together
  again next session vs. neighbor pairs vs. any pair).
- **Trap return** (`traps.py`): distance to previous capture (within session / later session) and, for
  animals caught in 5+ sessions, captures per session and same-trap return over their first five sessions.
  Exact-layout grids only (excludes 7, 12–14).

## Not yet used

- Home-range index files, *Blarina* monthly densities, raw-data compilations
- Unpublished manuscripts and appendices (PDF)
