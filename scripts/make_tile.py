"""Render the Projects-page tile: 25 years of vole density, one strip per habitat.

    python scripts/make_tile.py <processed_dir> <out_png>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

COL = {"MO": "#199e70", "MP": "#9085e9"}
BG = "#0d0d0d"


def main(proc, out):
    m = pd.read_csv(f"{proc}/population_main_series.csv")
    m["t"] = m["year"] + (m["month"] - 0.5) / 12
    fig, axes = plt.subplots(3, 1, figsize=(8, 4.5), dpi=150, sharex=True)
    fig.patch.set_facecolor(BG)
    for ax, hab in zip(axes, ["alfalfa", "bluegrass", "tallgrass"]):
        ax.set_facecolor(BG)
        for sp in ["MP", "MO"]:
            q = m[(m.habitat == hab) & (m.species == sp)].sort_values("t")
            # break the line across missing months
            q = q.set_index("t")["density_ha"]
            gaps = q.index.to_series().diff() > 0.2
            q[gaps] = float("nan")
            ax.plot(q.index, q.values, color=COL[sp], lw=1.4, solid_joinstyle="round")
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(colors="#898781", length=0, labelsize=8)
        ax.text(1971.8, ax.get_ylim()[1] * 0.85, hab.capitalize(), color="#c3c2b7", fontsize=9)
    axes[-1].set_xlim(1971.6, 1997.5)
    fig.subplots_adjust(left=0.03, right=0.98, top=0.97, bottom=0.07, hspace=0.08)
    fig.savefig(out, facecolor=BG)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
