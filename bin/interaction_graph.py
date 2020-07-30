#!/usr/bin/env python

import os
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout

translation = {"NTR": "NTR+",
               "LYS": "LYS+",
               "ARG": "ARG+",
               "GLU": "GLU-",
               "HEM": "HEM+",
               "HEC": "HEC+",
               "ASP": "ASP-",
               "CTR": "CTR-",
               "HIS": "HIS+",
               "TYR": "TYR-"}

def network(df, cutoff):
    column_names = df.keys()
    residues = list(column_names[2:])

    row_names = df["Terms"]
    i_start = np.where(row_names == "TOTAL")[0][0] + 1
    influencers = list(row_names[i_start:])

    df.set_index('Terms', inplace=True)

    G = nx.DiGraph()

    for influencer in influencers:
        if influencer[:3] in translation:
            influencer_name = translation[influencer[:3]] + influencer[3:]
        else:
            influencer_name = influencer
        for residue in residues:
            weight = abs(df.at[influencer, residue])
            if weight > cutoff:
                G.add_edge(influencer_name, residue, weight=weight)

    plt.figure(figsize=(12,8))
    pos = graphviz_layout(G)
    nx.draw(G, pos=pos, with_labels=True)
    plt.show()
    return


if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Plot the residue pKa influence graph. A residue's pKa is influenced by other residues, especially the charged ones. This kind of influence is plotted by an weighted arrow. This program requires mfe_all.py output."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-c", metavar="pw_cutoff", default=0.1, help="Cutoff pairwise interaction, default is 0.1.", type=float)
    parser.add_argument("mfe_all", metavar="mfe_all.xlsx", default="mfe_all.xlsx", help="mfe output file", nargs="?")
    args = parser.parse_args()

    pathname = args.mfe_all
    filename, file_extension = os.path.splitext(pathname)

    mfe_data = ""
    if file_extension == ".csv":
        mfe_data = pd.read_csv(pathname)
    elif file_extension == ".xlsx":
        mfe_data = pd.read_excel(pathname)
    elif file_extension == ".htm" or file_extension == ".html":
        mfe_data = pd.read_html(pathname)
    else:
        print("Can not interpret data file %s. Only xlsx, csv, and html files are acceppted." % pathname)

    cutoff = args.c

    network(mfe_data, cutoff)
