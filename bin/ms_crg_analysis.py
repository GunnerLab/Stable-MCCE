#!/usr/bin/env python

__doc__ = """
Module `ms_crg_analysis.py` holds functions/classes to process MCCE
microstates into charge microstates for residues of interest (default:
ionizable residues), or for all residues if list variable `res_of_interest`
is empty.

There are two tasks defined as functions that the cli can call:
  - Produce a charged microstate analysis using weighted correlation,
    along with various plots.
  - Return the <top_n> microstates into a csv file.

Note:
The -top_n command line argument acts as a switch: if a value is given, the
top n charge microstates are calculated and saved in a csv file. To obtain the
charge-microstates analysis with correlation, do not include -top_n .

Output files:
Plots:
  Histogram figure: enthalpy_dis.pdf
  Plot of unique charge distribution: all_en_cr_vs_log(count).pdf
  Correlation heat map: corr.pdf
CSV files:
  all_res_crg.csv
  all_res_crg_count.csv
  top_<top_n>_crg_ms.csv
"""

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict
import logging
import math
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import pandas as pd
from pathlib import Path
import seaborn as sns
from scipy.stats import skewnorm, rankdata
import sys
from typing import Tuple, Union


logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s]: %(name)s, %(funcName)s:\n\t%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    filename="crgms.log",
    encoding="utf-8",
)
logger = logging.getLogger(__name__)


try:
    import ms_analysis as msa
except Exception as e1:
    try:
        if Path(__file__).parent.name == "crgms":
            from crgms import ms_analysis as msa
    except Exception as e2:
        logger.exception(f"Error importing ms_analysis as msa:\n{e1}\n")
        logger.exception(f"Error importing ms_analysis from crgms:\n{e2}\n")
        sys.exit(2)



IONIZABLES = ["ASP", "GLU", "ARG", "HIS", "LYS", "CYS", "TYR", "NTR", "CTR"]
ACIDS = ["ASP","GLU"]
BASES = ["ARG","HIS","LYS"]
POLARS = ["CYS","TYR"]
res3_to_res1 = {"ASP":"D",
                "GLU":"E",
                "ARG":"R",
                "HIS":"H",
                "LYS":"K",
                "CYS":"C",
                "TYR":"Y",
               }


class WeightedCorr:
    def __init__(self,
                 xyw_df: pd.DataFrame=None,
                 x: pd.Series=None, y: pd.Series=None, w: pd.Series=None,
                 df: pd.DataFrame=None, wcol: str=None):
        """Class for Weighted Correlation.
        Either supply xyw_df, or (x, y, w), or (df, wcol).
        Call example:
            WeightedCorr(xyw_df=mydata[[x, y, w]])(method='pearson')

        :param xyw_df: pd.DataFrame with shape(n, 3) containing x, y, and w columns;
                       (column names irrelevant)
        :param x: pd.Series (n, ) containing values for x
        :param y: pd.Series (n, ) containing values for y
        :param w: pd.Series (n, ) containing weights
        :param df: pd.Dataframe (n, m+1) containing m phenotypes and a weight column
        :param wcol: column name of the weight column in the df dataframe.
        """

        self.x = None
        self.y = None
        self.w = None
        self.df = None

        no_args = all(arg is None for arg in [xyw_df, x, y, w, df, wcol])
        if no_args:
            logger.error("ValueError: No data supplied")
            sys.exit(2)

        no_series = all(arg is None for arg in [x, y, w])
        no_full_df = (df is None) and (wcol is None)

        if no_full_df:
            if no_series:
                xyw_df.dropna()
            else:
                xyw_df = pd.concat([x, y, w], axis=1).dropna()

            self.x, self.y, self.w = (pd.to_numeric(xyw_df[c], errors="coerce").values
                                      for c in xyw_df.columns)
            self.df = None
        else:
            if wcol not in df.columns:
                logger.error(f"KeyError: {wcol!r} not found in columns of df")
                sys.exit(2)
            self.df = df.loc[:, [c for c in df.columns if c != wcol]]
            self.w = pd.to_numeric(df.loc[:, wcol], errors="coerce")

    def _wcov(self, x, y, ms):
        return np.sum(self.w * (x - ms[0]) * (y - ms[1]))

    def _pearson(self, x=None, y=None):
        x, y = (self.x, self.y) if (x is None) and (y is None) else (x, y)
        mx, my = (np.sum(i*self.w)/np.sum(self.w) for i in [x, y])
        return self._wcov(x,y, [mx, my])/np.sqrt(self._wcov(x, x,
                                                             [mx, mx])*self._wcov(y,y, [my, my]))

    def _wrank(self, x):
        (unique, arr_inv, counts) = np.unique(rankdata(x),
                                              return_counts=True,
                                              return_inverse=True)
        a = np.bincount(arr_inv, self.w)
        return (np.cumsum(a) - a)[arr_inv]+((counts + 1)/2 * (a/counts))[arr_inv]

    def _spearman(self, x=None, y=None):
        x, y = (self.x, self.y) if (x is None) and (y is None) else (x, y)
        return self._pearson(self._wrank(x), self._wrank(y))

    def __call__(self, method: str="pearson") -> Union[float, pd.DataFrame]:
        """
        Args:
          method (str="pearson"): Correlation method to be used:
          "pearson" for pearson r, "spearman" for spearman rank-order correlation.
        Return:
          The correlation value (float) if class initialized with xyw_df, or (x, y, w).
          The correlation matrix as pd.DataFrame (m, m), if class initialized with (df, wcol).
        """

        if method not in ["pearson", "spearman"]:
            logger.error(f"ValueError: correlation method not in ['pearson','spearman']")
            sys.exit(2)

        cor = {"pearson": self._pearson, "spearman": self._spearman}[method]
        if self.df is None:
            return cor()
        else:
            out = pd.DataFrame(np.nan, index=self.df.columns, columns=self.df.columns)
            for i, x in enumerate(self.df.columns):
                for j, y in enumerate(self.df.columns):
                    if i >= j:
                        out.loc[x, y] = cor(x=pd.to_numeric(self.df[x], errors='coerce'),
                                            y=pd.to_numeric(self.df[y], errors='coerce'))
                        out.loc[y, x] = out.loc[x, y]
            return out


# >>> plotting functions ....................................
def plot_hist_by_ms_energy(ms_by_enrg: list, save_dir:Path):

    energy_lst_count = np.asarray([a for a,f in zip([x[0] for x in ms_by_enrg],
                                                    [x[1] for x in ms_by_enrg]) for _ in range(f)])

    #mu, sigma = norm.fit(energy_lst_count)  #unused
    skewness, mean, std = skewnorm.fit(energy_lst_count)

    fig = plt.figure(figsize = (10,8))
    fs = 15 # fontsize
    graph_hist = plt.hist(energy_lst_count, bins=100, alpha=0.6)
    Y = graph_hist[0]
    #y = skewnorm.pdf(np.array(energy_lst_count), skewness, mean, std)
    y = skewnorm.pdf(energy_lst_count, skewness, mean, std)

    pdf_data = Y.max()/max(y)*y
    plt.plot(energy_lst_count, pdf_data, label="approximated skewnorm", color="k")
    plt.title(f"{skewness = :.2f} {mean = :.2f} {std = :.2f}", fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.xlabel("Microstate Energy (Kcal/mol)", fontsize=fs)
    plt.ylabel("Count", fontsize=fs)
    plt.tick_params(axis="x", direction="out", length=8, width=2)
    plt.tick_params(axis="y", direction="out", length=8, width=2)

    fig_fp = save_dir.joinpath("enthalpy_dis.pdf")
    fig.savefig(fig_fp,
                dpi=300, bbox_inches="tight")
    logging.info(f"Histogram figure saved as {fig_fp}")

    return


def plots_unique_crg_histogram(charge_ms_info: tuple,
                               background_charge: float,
                               save_dir: Path):
    """
    Visualize which tautomer charge state is most populated.
    This includes the background charge.
    """

    x_av = [sum(x) + background_charge for x in charge_ms_info[0]]
    y_av = [math.log10(x) for x in charge_ms_info[1]]
    energy_diff_all_fl = [float(x) for x in charge_ms_info[3]]

    g1 = sns.JointGrid(marginal_ticks=True, height=6)
    ax = sns.scatterplot(x=x_av, y=y_av,
                         hue=energy_diff_all_fl,
                         palette='viridis',
                         size=energy_diff_all_fl,
                         sizes=(10, 200),
                         ax=g1.ax_joint)

    ax.set_xticks(range(int(min(x_av)), int(max(x_av)) + 1))
    ax.set_xlabel("Charge",fontsize=15)
    ax.set_ylabel("log$_{10}$(Count)",fontsize=16)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    ax2 = sns.histplot(x=x_av, linewidth=2, discrete=True, ax=g1.ax_marg_x)
    ax2.set_ylabel(None,fontsize=16)
    g1.ax_marg_y.set_axis_off()
    g1.fig.subplots_adjust(top=0.9)
    g1.fig.suptitle("All microstate energy", fontsize = 16)

    fig_fp = save_dir.joinpath("all_en_cr_vs_log(count).pdf")
    g1.savefig(fig_fp, dpi=300, bbox_inches="tight")
    logging.info("Plot of unique charge distribution saved as {fig_fp}")

    return


def corr_heat_map(df_corr: pd.DataFrame, save_dir: Path):

    plt.figure(figsize=(25, 8))
    cmap = ListedColormap(["darkred","red","orange","lightgray","skyblue","blue","darkblue"])
    norm = BoundaryNorm([-1.0, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 1.0], cmap.N)
    heatmap = sns.heatmap(df_corr, cmap=cmap, norm=norm, square=True,
                          linecolor="gray", linewidths=.01, fmt=".2f",
                          annot=True, annot_kws={"fontsize":12}
                         )
    plt.ylabel(None)
    plt.xlabel(None)
    plt.yticks(fontsize=15, rotation=0)
    plt.xticks(fontsize=15, rotation=90)
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)

    fig_fp = save_dir.joinpath("corr.pdf")
    plt.savefig(fig_fp, dpi=300, bbox_inches="tight")
    logging.info("Correlation heat map saved as {fig_fp}")

    return
# <<< plotting functions - end ....................................


def fixed_residues(fixed_iconfs: list,
                   conformers: list,
                   res_of_interest: list=IONIZABLES) -> Tuple[float, pd.DataFrame]:
    """
     Return a 2-tuple: fixed res background charge,
                       fixed res crg df for res listed in res_of_interest.
    """

    dd = defaultdict(float)
    for conf in conformers:
        if conf.iconf in fixed_iconfs:
            dd[conf.resid] = conf.crg
    fixed_backgrd_charge = sum(dd.values())
    if res_of_interest:
        dd = {k:dd[k] for k in dd if k[:3] in res_of_interest}
    fixed_res_crg_df = pd.DataFrame(dd.items(), columns=["Residue","crg"])

    return fixed_backgrd_charge, fixed_res_crg_df


def id_crg_relation(conformers: list) -> dict:

    return dict((conf.iconf, conf.crg) for conf in conformers)


def convert_ms_crg(l, d):

    crg_lst =[[y[0], y[1], [convert_ms_crg(x, d)
                            if isinstance(x, list) else d.get(x, x)
                            for x in y[2]
                            ]] for y in l]
    return crg_lst


def find_uniq_crgms_count_order(crg_list_ms: list,
                                E_range: Tuple[float,float]=None) -> Tuple[list,list,list,list]:
    """
    Args:
      crg_list_ms: the charge microstates list; assumed to be sorted in increasing order.
      E_range: To filter charge id based on the energy.

    Return:
      A 4-tuple: all_crg_ms_unique, all_count, unique_crg_state_order, energy_diff_all.

      `unique_crg_state_order` gives the order/rank of unique charge state based on energy.
      The lowest energy charge state will give the order 1 and then second unique charge
      state will give the order 2.
    """

    if E_range is None:
        logging.info("All energy microstates are selected.")
        begin_energy = crg_list_ms[0][0]
        end_energy = crg_list_ms[-1][0]
    else:
        begin_energy, end_energy = E_range
        if begin_energy and end_energy:
            crg_list_ms = [[x[0], x[1], x[2]]
                           for x in crg_list_ms
                           if x[0] >= begin_energy and x[0] <= end_energy
                           ]
        else:
            logging.critical("No energy bound is defined")
            sys.exit('Give the lower or upper energy bound.')


    # unique charge as key and energy, count and order
    crg_all_count = {}
    unique_crg_state_order = 1

    for x, array in enumerate(crg_list_ms):
        if tuple(array[2]) not in crg_all_count:
            crg_all_count[(tuple(array[2]))] = [array[1], [array[0]], [unique_crg_state_order]]
            unique_crg_state_order += 1
        else:
            crg_all_count[(tuple(array[2]))][0] += array[1]

            # add the maximum and minimum energy
            min_energy = min(min(crg_all_count[(tuple(array[2]))][1]), array[0])
            max_energy = max(max(crg_all_count[(tuple(array[2]))][1]), array[0])

            # clear energy list and append minimum and maximum energy
            crg_all_count[(tuple(array[2]))][1].clear()
            crg_all_count[(tuple(array[2]))][1].append(min_energy)
            crg_all_count[(tuple(array[2]))][1].append(max_energy)

    # make a list of count, unique charge microstate, energy difference and order.
    all_crg_ms_unique = []
    all_count = []
    energy_diff_all = []
    unique_crg_state_order = []

    for k in crg_all_count:
        v = crg_all_count[k]
        all_crg_ms_unique.append(list(k))
        all_count.append(v[0])
        unique_crg_state_order.append(v[2][0])
        if len(v[1]) == 2:
            energy_diff_all.append(round(v[1][1]-v[1][0], 6))
        elif len(v[1]) == 1:
            energy_diff_all.append(0)
        else:
            logging.critical("Error in unique charge state creation.")
            sys.exit("Error in unique charge state creation.")

    logging.info(f"Count of charge-ms: {len(crg_list_ms)}")
    logging.info(f"Count of unique charge-ms: {len(all_crg_ms_unique)}")

    return all_crg_ms_unique, all_count, unique_crg_state_order, energy_diff_all


def concat_crgms_dfs(unique_crgms: list,
                     ms_count: list,
                     ms_order: list,
                     free_res_df: pd.DataFrame,
                     background_charge: float,
                     res_of_interest: list=IONIZABLES) -> pd.DataFrame:

    uniq_crg_df = pd.DataFrame(unique_crgms).T
    ms_count_df = pd.DataFrame(ms_count, columns=["Count"]).T
    ms_order_df = pd.DataFrame(ms_order, columns=["Order"]).T

    crg_ms_count_df = pd.concat([uniq_crg_df, ms_count_df, ms_order_df])
    crg_count_res_1 = pd.concat([free_res_df, crg_ms_count_df], axis=1)
    crg_count_res_1.loc["Count", "Residue"] = "Count"
    crg_count_res_1.loc["Order", "Residue"] = "Order"

    all_crg_count_res = crg_count_res_1.set_index("Residue")
    all_crg_count_res = all_crg_count_res.sort_values(by="Count", axis=1,
                                                      ascending=False)
    # rename columns
    all_crg_count_res.columns = range(all_crg_count_res.shape[1])
    all_crg_count_res = all_crg_count_res.T.set_index("Order")
    all_crg_count_res.index = all_crg_count_res.index.astype(int)
    all_crg_count_res["Occupancy"] = ((all_crg_count_res["Count"]/sum(all_crg_count_res["Count"]))
                                      .round(3)
    )
    all_crg_count_res['Sum_crg_protein'] = (all_crg_count_res.iloc[:,:-2].sum(axis=1)
                                            + background_charge
    )

    if res_of_interest:
        cols_to_drop = [c for c in all_crg_count_res.columns
                        if (c[:3] not in res_of_interest)
                        and (c not in ["Occupancy","Count","Sum_crg_protein"])
                       ]
        out_df = all_crg_count_res.drop(cols_to_drop, axis=1).sort_index()
        return out_df

    return all_crg_count_res.sort_index()


def combine_free_fixed_residues(fixed_residue_df: pd.DataFrame,
                                free_res_crg_count_df: pd.DataFrame) -> pd.DataFrame:

    df_fixed_res_crg = fixed_residue_df.T
    df_fixed_res_crg.columns = df_fixed_res_crg.iloc[0]
    df_fixed_res_crg = df_fixed_res_crg.iloc[1:,:].reset_index(drop=True)
    df_fixed_res_dup = pd.concat([df_fixed_res_crg]*len(free_res_crg_count_df),
                                 ignore_index=True)
    df_all = free_res_crg_count_df.join(df_fixed_res_dup)

    return df_all


def correlation_data_parsing(df: pd.DataFrame) -> pd.DataFrame:

    all_crg_count = df.iloc[:,:-2]
    all_crg_count = all_crg_count.T
    all_crg_count["std"] = all_crg_count.std(axis=1).round(3)
    all_crg_count_std = (all_crg_count.loc[all_crg_count["std"] != 0].T[:-1]
                         .reset_index(drop=True)
    )
    logging.info(("Number of residues that change the protonation state: ",
                  f"{len(all_crg_count_std.columns)-1}"))

    return all_crg_count_std


def rename_order_residues(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return the amended dataframe:
        Rename the residues with shorter name;
        Reorder columns by acid first, then polar,
        then base residues and then MQ in columns.
    """

    rename_dict = {}
    acid_list = []
    base_list = []
    polar_rest_list = []
    ub_q_list = []
    non_res_list = []

    for c in df.columns[:-1]:
        rename_dict[c] = f"{c[3]}_{c[:3]}{c[4:8]}{c[8:]}"

    rename_dict["Count"] = "Count"

    for k in rename_dict:
        v = rename_dict[k]
        RES = v[2:5]

        if RES in IONIZABLES[:-2]:  # w/o terminii
            rename_dict[k] = v[:1] + res3_to_res1[RES] + v[5:]
            if RES in ACIDS:
                acid_list.append(rename_dict[k])
            elif RES in BASES:
                base_list.append(rename_dict[k])
            elif RES in POLARS:
                polar_rest_list.append(rename_dict[k])
        else:
            if RES == 'MQ8':
                rename_dict[k] = "MQ" + v[5:]
                ub_q_list.append(rename_dict[k])
            else:
                non_res_list.append(v)

    df = df.rename(rename_dict, axis=1)
    col_order_list = [*acid_list, *polar_rest_list, *base_list, *ub_q_list, *non_res_list]

    return df[col_order_list]


def filter_weightedcorr(df: pd.DataFrame, cutoff: float=0.0) -> pd.DataFrame:
    """Filter the weighted correlation of df by cutoff.
    Uses WeightedCorr class.
    """

    wc_df = WeightedCorr(df=df, wcol="Count")(method="pearson")
    for i in wc_df.columns:
        if list(abs(wc_df[i]) >= cutoff).count(True) == 1:
            wc_df.drop(i, inplace=True)
            wc_df.drop(i, axis=1, inplace=True)

    return wc_df


def check_mcce_input_files(mcce_dir: Path, ph_pt: str) -> tuple:

    h3_fp = mcce_dir.joinpath("head3.lst")
    if not h3_fp.exists():
        sys.exit(f"FileNotFoundError: head3.lst not found in {h3_fp.parent}.")

    msout_fp = mcce_dir.joinpath("ms_out", f"pH{ph_pt}eH0ms.txt")
    if not msout_fp.exists():
        sys.exit(f"FileNotFoundError: ms_out file {msout_fp} not found in {msout_fp.parent}.")

    return h3_fp, msout_fp


def get_topN_crg_microstates(mcce_dir: Path,
                             ph_pt: int=7,
                             top_n=3,
                             res_of_interest: list=IONIZABLES,
                             to_csv: bool=True,
                             return_df: bool=False):
    """Process MCCE microstates and return the top N charge microstates
    for the residues in res_of_interest.
    """

    h3_fp, msout_fp = check_mcce_input_files(mcce_dir, ph_pt)

    conformers = msa.read_conformers(h3_fp)
    mc = msa.MSout(msout_fp)
    msg = (f"{len(conformers) = }\n",
           f"Total microstates (total count): {mc.N_ms:,}\n",
           f"Total unique conformers ms: {mc.N_uniq:,}\n"
           )
    logger.info(msg)

    ms_by_E = mc.sort_microstates()
    ms_orig_lst = [[ms.E, ms.count, ms.state] for ms in ms_by_E]
    ms_free_res_df = msa.free_residues_df(mc.free_residues,
                                          conformers,
                                          colname="Residue")
    background_charge, fixed_res_crg_df = fixed_residues(mc.fixed_iconfs,
                                                         conformers,
                                                         res_of_interest)
    id_vs_charge = id_crg_relation(conformers)
    crg_orig_lst = convert_ms_crg(ms_orig_lst, id_vs_charge)
    charge_ms_info = find_uniq_crgms_count_order(crg_orig_lst)  #4-tuple
    all_res_crg_count_df = concat_crgms_dfs(charge_ms_info[0],
                                            charge_ms_info[1],
                                            charge_ms_info[2],
                                            ms_free_res_df,
                                            background_charge,
                                            res_of_interest)
    if top_n <= 0:
        top_n = 1
        logger.info("Negative top_n reset to 1.")

    top_df = all_res_crg_count_df[:top_n].T
    if to_csv:
        top_fp = mcce_dir.joinpath(f"top_{top_n}_crg_ms.csv")
        top_df.to_csv(top_fp)
        logger.info(f"Top {top_n} crg ms file: {top_fp}")

    if return_df:
        return top_df

    return


def crg_msa_with_correlation(mcce_dir: Path,
                             ph_pt: int=7,
                             res_of_interest: list=IONIZABLES,
                             corr_cutoff: float=0.):

    h3_fp, msout_fp = check_mcce_input_files(mcce_dir, ph_pt)

    conformers = msa.read_conformers(h3_fp)
    mc = msa.MSout(msout_fp)
    msg = (f"{len(conformers) = }\n",
           f"Total microstates (total count): {mc.N_ms:,}\n",
           f"Total unique conformers ms: {mc.N_uniq:,}\n"
           )
    ms_by_E = mc.sort_microstates()
    ms_orig_lst = [[ms.E, ms.count, ms.state] for ms in ms_by_E]

    plot_hist_by_ms_energy(ms_orig_lst, mcce_dir)

    ms_free_res_df = msa.free_residues_df(mc.free_residues,
                                          conformers,
                                          colname="Residue")

    background_charge, fixed_res_crg_df = fixed_residues(mc.fixed_iconfs,
                                                         conformers,
                                                         res_of_interest)

    id_vs_charge = id_crg_relation(conformers)
    crg_orig_lst = convert_ms_crg(ms_orig_lst, id_vs_charge)

    charge_ms_info = find_uniq_crgms_count_order(crg_orig_lst)
    plots_unique_crg_histogram(charge_ms_info,
                               background_charge,
                               mcce_dir)

    free_res_crg_count_df = concat_crgms_dfs(charge_ms_info[0],
                                             charge_ms_info[1],
                                             charge_ms_info[2],
                                             ms_free_res_df,
                                             background_charge,
                                             res_of_interest)

    res_crg_csv = mcce_dir.joinpath("all_res_crg.csv")
    combine_free_fixed_residues(fixed_res_crg_df,
                                free_res_crg_count_df).to_csv(res_crg_csv)

    res_crg_count_csv = mcce_dir.joinpath("all_res_crg_count.csv")
    free_res_crg_count_df.to_csv(res_crg_count_csv, header=True)

    all_crg_count_std = correlation_data_parsing(free_res_crg_count_df)
    df = rename_order_residues(all_crg_count_std)
    df_correlation = filter_weightedcorr(df, cutoff=corr_cutoff)
    corr_heat_map(df_correlation, mcce_dir)

    logger.info("Charge microstate analysis with correlation over.") 

    return


def crgmsa_parser() -> ArgumentParser:

    p = ArgumentParser(
        prog="ms_crg_analysis.py",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-run_dir",
        type=str,
        default=".",
        help="titration point; default: %(default)s."
    )
    p.add_argument(
        "-ph_pt",
        type=int,
        default=7,
        help="titration point; default: %(default)s."
    )
    p.add_argument(
        "-res_of_interest",
        nargs='+',  # 1 or more
        type=str,
        default=IONIZABLES,
        help="List of residues of interest; default: %(default)s."
    )
    p.add_argument(
        "-corr_cutoff",
        type=float,
        default=0.,
        help="Correlation cutoff value; default: %(default)s."
    )
    p.add_argument(
        "-top_n",
        type=int,
        # if not None: call get_topN_crgma
        default=None,
        help="""Save the top n charge ms to a file named top_<top_n>_crg_ms.csv.
        This option is a switch: if present, get_topN_crg_microstates is called.
        if not, crg_msa_with_correlation is called; %(default)s.
        """
    )
    return p


def crgmsa_cli(argv=None):

    parser = crgmsa_parser()
    args = parser.parse_args(argv)

    args.run_dir = Path(args.run_dir).resolve()
    if args.top_n is not None:
        get_topN_crg_microstates(args.run_dir,
                                 ph_pt=args.ph_pt,
                                 top_n=args.top_n,
                                 res_of_interest=args.res_of_interest)
    else:
        crg_msa_with_correlation(args.run_dir,
                                 ph_pt=args.ph_pt,
                                 res_of_interest=args.res_of_interest,
                                 corr_cutoff=args.corr_cutoff)


if __name__ == "__main__":

    crgmsa_cli(sys.argv[1:])
