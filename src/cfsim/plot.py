import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sb
from dataclasses import dataclass
from cfsim.colours import colours


@dataclass
class AxesSettings:
    df: pd.DataFrame
    yvar: str
    ylabel: str | None = None
    title: str | None = None
    legend: bool = False

    def __post_init__(self):
        if self.yvar == 'baf':
            self.df['baf'] = self.df['b'] / (self.df['b'] + self.df['a'])


def plot_cfdna(
    cfdna: pd.DataFrame,
    out_file: str | None = None,
    show: bool = False
) -> None:
    
    axes = [
        AxesSettings(cfdna, yvar='reads'),
        AxesSettings(cfdna, yvar='rdr'),
        AxesSettings(cfdna, yvar='baf'),

    ]
    
    fig = plt.figure(figsize=(16, 2 * len(axes)))

    grid = fig.add_gridspec(len(axes), 1, hspace=0.1)

    chroms = sort_chroms(cfdna["chrom"].unique())

    chroms_size = cfdna["chrom"].value_counts()

    width_ratios = [chroms_size[x] for x in chroms]

    for i, a in enumerate(axes):
        sub_grid = grid[i].subgridspec(
            nrows=1, ncols=len(chroms), width_ratios=width_ratios, wspace=0.05
        )

        plot_data(
            df=a.df,
            yvar=a.yvar,
            chroms=chroms,
            fig=fig,
            grid=sub_grid,
            title=a.title,
            ylabel=a.ylabel,
            legend=a.legend,
        )

    fig.align_labels()

    fig.supxlabel("Chromosome", x=0.525)

    grid.tight_layout(fig)

    fig.savefig(out_file, dpi=150, bbox_inches="tight")

    if show:

        plt.close()

    plt.close()


def plot_data(
    df: pd.DataFrame,
    yvar: str, 
    chroms: list[str],
    fig: plt.Figure,
    grid: plt.GridSpec,
    ylabel: str | None = None,
    title: str | None = None,
    legend: bool = False,
):
    y_max = df[yvar].max()

    y_min = df[yvar].min()
    
    grouped = df.groupby("chrom")

    for i, chrom in enumerate(chroms):

        chrom_df = grouped.get_group(chrom)

        chrom_df = chrom_df.sort_values(by=["start"])

        num_bins = chrom_df.shape[0]

        chrom_df["idx"] = np.arange(num_bins)

        ax = fig.add_subplot(grid[0, i])
        
        ax.scatter(
            chrom_df["idx"],
            chrom_df[yvar],
            c=colours["orange"],
            s=1,
        )

        # SETUP SPLINES

        ax.spines["left"].set_position(("outward", 10))

        ax.spines["bottom"].set_position(("outward", 10))

        ax.spines["left"].set_color("black")

        ax.spines["bottom"].set_color("black")

        ax.spines["top"].set_visible(False)

        ax.spines["right"].set_visible(False)

        ax.xaxis.tick_bottom()

        ax.yaxis.tick_left()

        ax.xaxis.grid(True, which="major", linestyle=":")

        ax.yaxis.grid(True, which="major", linestyle=":")

        sb.despine(ax=ax, offset=10)

        ax.spines["top"].set_visible(False)

        ax.spines["right"].set_visible(False)

        ax.xaxis.grid(False)

        if i != 0:
            ax.spines["left"].set_visible(False)

            ax.tick_params(axis="y", labelleft=False, left=False)

        else:
            ax.tick_params(axis="x", which="major", labelsize=12)

            ax.set_ylabel(yvar if ylabel is None else ylabel)

        ax.set_xticks([num_bins / 2])

        ax.set_xticklabels([chrom.replace("chr", "")], fontsize=12)

        ax.set_ylim(y_min, y_max)

    if title is not None:
        ax = fig.add_subplot(grid[:])

        ax.axis("off")

        ax.set_title(title)



def sort_chroms(chroms: list[str]) -> list[str]:
    """sort_chroms adapted from: https://github.com/Roth-Lab/hapclone-smk/blob/main/scripts/plot_clone_pseudobulk.py
    """
    numeric = []
    string = []

    if chroms[0].startswith("chr"):
        chr_prefix = True
    else:
        chr_prefix = False

    for c in chroms:
        if chr_prefix:
            c = c.replace("chr", "")
        try:
            numeric.append(int(c))
        except ValueError:
            string.append(c)

    chroms = [str(x) for x in sorted(numeric)] + list(sorted(string))

    if chr_prefix:
        chroms = ["chr{}".format(x) for x in chroms]

    return chroms