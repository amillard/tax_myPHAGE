# Description: Plotting functions for tax_myPHAGE

############################################################################################################
# Imports
############################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import matplotlib.colors as mcolors
import pandas as pd

from typing import Dict

############################################################################################################
# Functions
############################################################################################################

def heatmap(dfM: pd.DataFrame, outfile: str, matrix_out:str, accession_genus_dict: Dict[str, str]):
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

    #define output files
    svg_out = outfile+".svg"
    pdf_out = outfile+".pdf"
    png_out = outfile+".png"

    # Plot the heatmap
    ax = plt.gca()

    # Add the genus to the accession
    dfM["A"] = dfM["A"].map(lambda x: x + ":" + accession_genus_dict.get(x, ""))
    dfM["B"] = dfM["B"].map(lambda x: x + ":" + accession_genus_dict.get(x, ""))

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

    for leave in dendrogram["ivl"]:
        if "query" in leave:
            query_leave = leave
        else:
            leaves_order.append(leave)

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

    # image
    # im = plt.imshow(df.values, cmap=cmap)
    im = plt.imshow(df.values, cmap=custom_cmap, norm=norm)

    ax.set_xticks(np.arange(df.shape[1]), labels=df.columns.tolist())
    ax.set_yticks(np.arange(df.shape[0]), labels=df.index.tolist())

    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(df.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(df.shape[0] + 1) - 0.5, minor=True)
    
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Set the size of the figure
    fig_width = max(4, df.shape[1] * 0.75)
    fig_height = max(4, df.shape[0] * 0.75)
    plt.gcf().set_size_inches(fig_width, fig_height)

    # Add the text
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            font_size = (
                min(fig_width, fig_height) / max(df.shape[0], df.shape[1])
            ) * 10

            ax.text(
                j,
                i,
                df.iloc[i, j],
                ha="center",
                va="center",
                color="w",
                fontsize=font_size,
            )

    # adjust figure size based on number of rows and columns

    # plot with padding
    plt.tight_layout(pad=2.0)

    plt.savefig(svg_out, bbox_inches="tight")
    plt.savefig(pdf_out, bbox_inches="tight")
    plt.savefig(png_out, dpi=300, bbox_inches="tight")
    plt.close()

    return

############################################################################################################
