"""
This module provides functionalities for plotting data in the taxmyphage package.
It includes functions to plot heatmaps of similarity matrices.
"""

####################################################################################################
# Imports
####################################################################################################

from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.cluster.hierarchy as sch
import pandas as pd
import seaborn as sns

sns.set_theme(style="white")


####################################################################################################
# Functions
####################################################################################################


def heatmap(
    dfM: pd.DataFrame,
    outfile: str,
    matrix_out: str,
    accession_genus_dict: Dict[str, str],
):
    """
    Plot a heatmap of the similarity matrix

    Args:
        dfM (pd.DataFrame): Similarity matrix
        outfile (str): Path to the output file
        matrix_out (str): Path to the output matrix
        accession_genus_dict (Dict[str, str]): Dictionary of accessions and genera

    Returns:
        None
    """

    # Set up matplotlib
    plt.rcParams["text.color"] = "#131516"
    plt.rcParams["svg.fonttype"] = "none"  # Editable SVG text
    # plt.rcParams["font.family"] = "sans-serif"
    # plt.rcParams["font.sans-serif"] = "Arial"
    # plt.rcParams["font.weight"] = "light"

    # define output files
    svg_out = outfile + ".svg"
    pdf_out = outfile + ".pdf"
    png_out = outfile + ".png"

    # Add the genus to the accession
    dfM["A"] = dfM["A"].map(lambda x: x + ":" + accession_genus_dict.get(x, ""))
    dfM["B"] = dfM["B"].map(lambda x: x + ":" + accession_genus_dict.get(x, ""))

    # Remove the ':' at the end of the index and columns names
    dfM["A"] = dfM["A"].str.rstrip(":")
    dfM["B"] = dfM["B"].str.rstrip(":")

    # Make the matrix symmetric
    dfM.update(dfM.loc[dfM.A > dfM.B].rename({"A": "B", "B": "A"}, axis=1))

    # Round the similarity to 2 decimals
    dfM = dfM.round(2)

    # Pivot the dataframe
    df = dfM.pivot(index="A", columns="B", values="sim").fillna(0)

    # Make the matrix symmetric
    df = df + df.T - np.diag(df.values.diagonal())

    # Perform hierarchical clustering
    Z = sch.linkage(df, method="ward")

    # Plot the dendrogram
    dendrogram = sch.dendrogram(Z, labels=df.index, no_plot=True)

    # Looking for the query leave to put at the end
    leaves_order = []

    # Get the query leave else empty string
    query_leave = ""

    for leave in dendrogram["ivl"]:
        if "query" in leave:
            query_leave = leave
        else:
            leaves_order.append(leave)

    if query_leave != "":
        leaves_order.append(query_leave)

    # Reorder the matrix
    df = df.loc[leaves_order, leaves_order]
    df.iloc[:, :] = np.triu(df.values, k=0)
    # Maybe the following method is faster
    # df = df.where(np.triu(np.ones(df.shape)).astype(np.bool))

    # Save the matrix
    df.to_csv(matrix_out, sep="\t", index=True)

    # Define the colors
    colors = ["white", "lightgray", "skyblue", "steelblue", "darkgreen"]
    boundaries = [0, 1, 50, 70, 95, 100]

    # Create a colormap
    norm = mcolors.BoundaryNorm(boundaries, len(colors))
    # Create the colormap
    custom_cmap = mcolors.ListedColormap(colors)

    mask = np.zeros_like(df)
    mask[np.tril_indices_from(mask, k=-1)] = True

    # adjust figure size based on number of rows and columns
    fig_width = max(3, df.shape[1] * 0.75)
    fig_height = max(3, df.shape[0] * 0.75)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # plot heatmap using seaborn
    sns.heatmap(
        df,
        cmap=custom_cmap,
        norm=norm,
        annot=True,
        fmt=".1f",
        cbar=False,
        square=True,
        linewidth=2,
        linecolor="white",
        mask=mask,
        ax=ax,
        # If you want to change the color of the text to white
        # annot_kws={"color": "white"},
    )

    ax.set_ylabel("")
    ax.set_xlabel("")

    ax.xaxis.set_ticks_position("top")

    # Rotate the tick labels and set their alignment.
    ax.set_xticklabels(ax.get_xticklabels(), rotation=-30, horizontalalignment="right")

    sns.despine(left=True, top=True, right=True, bottom=True)

    plt.savefig(svg_out, bbox_inches="tight")
    plt.savefig(pdf_out, bbox_inches="tight")
    plt.savefig(png_out, dpi=300, bbox_inches="tight")
    plt.close()

    return


####################################################################################################
