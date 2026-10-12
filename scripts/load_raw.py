"""Load Getz's 25-year original field trap records into one tidy DataFrame.

Source: https://www.life.illinois.edu/getz/25-year/microtus/original_trap_records/
Each file (MOnn / MPnn) is tab-delimited with 11 columns; some use old-Mac CR
line endings, so we normalise before parsing.
"""
import glob
import io
import os
import re

import pandas as pd

COLS = ["animal", "sex", "repro", "mass", "date", "ampm", "station",
        "new_recap", "species", "site", "session"]


def read_trap_file(path):
    raw = open(path, "rb").read().decode("latin-1")
    raw = raw.replace("\r\n", "\n").replace("\r", "\n")
    df = pd.read_csv(io.StringIO(raw), sep="\t", header=0, names=COLS,
                     dtype=str, keep_default_na=False, on_bad_lines="warn",
                     index_col=False)
    name = os.path.basename(path)
    m = re.match(r"(M[OP])\s*(\d+)", name, re.I)
    df["file"] = name
    df["file_species"] = m.group(1).upper()
    df["file_site"] = int(m.group(2))
    return df


def load_all(mirror_dir):
    base = os.path.join(mirror_dir, "25-year/microtus/original_trap_records")
    files = sorted(glob.glob(f"{base}/M?_files/*.txt"))
    return pd.concat([read_trap_file(f) for f in files], ignore_index=True)
