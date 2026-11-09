import os
import json
import shutil
from collections import Counter
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
import math
from itertools import combinations

from formulaic.utils.iterators import peekable_iter
from scipy.stats import chi2_contingency

from scipy.stats import ranksums
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.spatial.distance import pdist, squareform
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import pandas as pd
plt.rcParams["svg.fonttype"] = "none"


def reconstruct(nums, seps, start, end):
    """
    根据数字位置和分隔符位置重新组合字符串
    nums[start:end+1]
    seps[start:end]
    """
    s = str(nums[start])
    for i in range(start, end):
        s += seps[i] + str(nums[i+1])
    return s


def parse_hor(hor_name, nums, seps):

    START = {1, 2}

    results = []
    n = len(nums)

    i = 0

    while i < n:

        if nums[i] in START:

            end_pos = None

            # 只在下一个START之前搜索
            j = i + 1
            while j < n and nums[j] not in START:

                if nums[j] == 10:
                    end_pos = j   # 10优先

                elif nums[j] == 9 and end_pos is None:
                    end_pos = j   # 暂存9

                elif nums[j] == 8 and end_pos is None:
                    end_pos = j   # 暂存8

                j += 1


            if end_pos is not None:
                results.append(reconstruct(nums, seps, i, end_pos))
                i = end_pos + 1

            else:
                results.append(reconstruct(nums, seps, i, i))
                i += 1

        else:

            start = i
            i += 1

            while i < n and nums[i] not in START:
                i += 1

            results.append(reconstruct(nums, seps, start, i-1))

    for i in range(len(results)):
        results[i] = "{}.{}".format(hor_name, results[i])
    return results


def deal_with_chr16_hor(hor_type, hor_cnt):

    split_hor_dict = {}

    hor_name = hor_type.split(".")[0]
    hor_combine = hor_type.split(".")[1]
    if "2-9" in hor_type or "2-10" in hor_type or "2-6_9" in hor_type or "1-6_9" in hor_type:

        hor_num_list = list(map(int, re.findall(r'\d+', hor_combine)))
        hor_sep_list = re.findall(r'[-_]', hor_combine)

        # print(hor_type)

        split_hor_list = parse_hor(hor_name, hor_num_list, hor_sep_list)

        for split_hor in split_hor_list:
            if split_hor not in split_hor_dict:
                split_hor_dict[split_hor] = 0

            split_hor_dict[split_hor] += hor_cnt

    else:
        split_hor_dict[hor_type] = hor_cnt

    return split_hor_dict


def load_hor_stats_file(input_path):
    stats_dict = {}

    current_chr = None

    with open(input_path) as fin:
        for line in fin:
            line = line.strip()

            if line.startswith("chr"):
                current_chr = line.split("_")[0]
                stats_dict[current_chr] = {}
            else:
                hor_type, hor_cnt = line.split("\t")
                if current_chr == "chr16":
                    if hor_type == "total":
                        # stats_dict[current_chr][hor_type] = int(hor_cnt)
                        stats_dict[current_chr][hor_type] = 0
                    else:

                        split_hor_dict = deal_with_chr16_hor(hor_type, int(hor_cnt))

                        for split_hor in split_hor_dict:
                            if split_hor not in stats_dict[current_chr]:
                                stats_dict[current_chr][split_hor] = 0
                            stats_dict[current_chr][split_hor] += split_hor_dict[split_hor]
                            stats_dict[current_chr]["total"] += split_hor_dict[split_hor]

                elif current_chr == "chr1":

                    if "S1C1/5/19H1L.1-5" in hor_type:
                        hor_type = "S1C1/5/19H1L.1-5(_6/4_5){n}-6"
                    if "S1C1/5/19H1L.1-3" in hor_type:
                        hor_type = "S1C1/5/19H1L.1-3(_6/4_5){n}-6"

                    if hor_type not in stats_dict[current_chr]:
                        stats_dict[current_chr][hor_type] = 0
                    stats_dict[current_chr][hor_type] += int(hor_cnt)


                else:
                    stats_dict[current_chr][hor_type] = int(hor_cnt)

    # print(dict(sorted(stats_dict["chr16"].items(), key=lambda x: x[1], reverse=True)))
    return stats_dict



def plot_principal_hor_analysis(feature_file, spec_chrom, cluster_col):

    cluster_colors = {-1: "#000000", 0: "#4C72B0", 1: "#DD8452", 2: "#55A868", 3: "#C44E52", 4: "#4a3d5c", 5: "#7f7f7f"}

    # print(spec_chrom)

    df = pd.read_csv(feature_file, sep="\t")

    # df = df[df[cluster_col] != -1].reset_index(drop=True)
    # df[cluster_col] = df[cluster_col].astype(str)

    features = c1_principle_hor_dict[spec_chrom] + ["aSat_length"]
    # features = chrom_selected_hor_dict[spec_chrom]

    sns.pairplot(
        df,
        vars=features,
        hue=cluster_col,
        palette=cluster_colors,
        diag_kind="hist",
        plot_kws={
            "s": 5, "edgecolor": "none",
        }
    )

    plt.tight_layout()
    plt.savefig(feature_file.replace(".txt", ".pha.png"), dpi=600)
    plt.close()


def plot_feature_compare(feature_file, cluster_col):
    def calc_frequency(df, popu_col, cluster_col):
        count = pd.crosstab(df[cluster_col], df[popu_col])
        freq = count.div(count.sum(axis=0), axis=1)
        return freq.T

    cluster_colors = {-1: "#dddddd", 0: "#4C72B0", 1: "#DD8452", 2: "#55A868", 3: "#C44E52", 4: "#4a3d5c", 5: "#7f7f7f"}

    df = pd.read_csv(feature_file, sep="\t")
    df = df[df[cluster_col] != -1].reset_index(drop=True)

    features = df.columns.tolist()[3: -1]
    # # optional: merge EAS_CPC and EAS_Other
    df_merge_EAS = df.copy()
    df_merge_EAS["hyper_popu"] = df_merge_EAS["hyper_popu"].replace({"EAS_CPC": "EAS", "EAS_Other": "EAS"})
    df_merge_EAS["detail_popu"] = df_merge_EAS["detail_popu"].replace({"EAS_CPC": "EAS", "EAS_Other": "EAS"})

    hyper_freq = calc_frequency(df_merge_EAS, "hyper_popu", cluster_col)
    detail_freq = calc_frequency(df_merge_EAS, "detail_popu", cluster_col)

    if len(features) + 2 > 9:
        fig, axes = plt.subplots(4, 4, figsize=(20, 20))
    else:
        fig, axes = plt.subplots(3, 3, figsize=(16, 16))

    axes = axes.flatten()

    # # hyper_popu grouped bar
    hyper_freq.plot(kind="bar", stacked=True, ax=axes[0], color=[cluster_colors[c] for c in hyper_freq.columns], width=0.8)
    axes[0].set_title("Hyper population", fontsize=18)
    axes[0].legend(title="", fontsize=12)
    axes[0].tick_params(axis="both", labelsize=18)

    # detail_popu grouped bar
    detail_freq.plot(kind="bar", stacked=True, ax=axes[1], color=[cluster_colors[c] for c in detail_freq.columns], width=0.8)
    axes[1].set_title("Detail population", fontsize=18)
    axes[1].legend(title="", fontsize=12)
    axes[1].tick_params(axis="both", labelsize=18)

    # # box plot
    clusters = sorted(df[cluster_col].unique())

    for i, feature in enumerate(features):
        ax = axes[i + 2]
        data = [df.loc[df[cluster_col] == c, feature].values for c in clusters]

        bp = ax.boxplot(data, tick_labels=clusters, widths=0.6, showfliers=False, patch_artist=True)
        for box, c in zip(bp["boxes"], clusters):
            box.set_facecolor(cluster_colors[c])
            box.set_alpha(0.6)

        ax.set_title(feature, fontsize=18)
        ax.tick_params(axis="both", labelsize=18)

        # ax.set_xlabel("Cluster")
        # ax.set_ylabel(feature)

    # # remove unused subplot
    # for ax in axes[7:]:
    #     ax.axis("off")

    plt.tight_layout()
    plt.savefig(feature_file.replace(".txt", ".compare.png"), dpi=600)
    plt.close()


def principal_hor_analysis_classify_c1(chrom_order, output_folder):


    for spec_chrom in chrom_order:


        # # initiate
        feature_dict = {}
        for hor in c1_principle_hor_dict[spec_chrom]:
            feature_dict[hor] = []
        feature_dict["aSat_length"] = []
        for oSat in oSat_list:
            feature_dict["{}_length".format(oSat)] = []

        feature_dict["cluster_1"] = []

        # # begin to traverse
        sample_cnt = 0
        sample_list = []
        sample_hyper_popu_list = []
        sample_detail_popu_list = []

        for sample_stats_file in os.listdir(stats_res_folder):
            if not sample_stats_file.endswith("#stv_stats.tsv"):
                continue

            sample = sample_stats_file.replace("#stv_stats.tsv", "")
            if sample not in all_included_samples:
                continue
            if "{}_RagTag".format(spec_chrom) not in whole_cent_dict[sample]:
                continue

            # if sample != "HG02622_mat_hprc_r2_v1.0.1.fa.gz":
            #     continue

            sample_stats_file = os.path.join(stats_res_folder, sample_stats_file)
            try:
                sample_hor_stats_dict = load_hor_stats_file(sample_stats_file)[spec_chrom]
            except KeyError:
                continue

            # print(sample_hor_stats_dict)
            # # calculate ratio
            ratio_dict = {}
            for hor in c1_principle_hor_dict[spec_chrom]:
                # ratio_dict[hor] = sample_hor_stats_dict.get(hor, 0)
                ratio_dict[hor] = round(sample_hor_stats_dict.get(hor, 0) / sample_hor_stats_dict["total"], 3)

            # # add feature of this sample to all
            for hor in c1_principle_hor_dict[spec_chrom]:
                feature_dict[hor].append(ratio_dict[hor])

            feature_dict["aSat_length"].append(sample_sSat_length_dict[sample][spec_chrom])

            for oSat in oSat_list:
                feature_dict["{}_length".format(oSat)].append(sample_oSat_length_dict[oSat][sample][spec_chrom])

            sample_cnt += 1
            sample_list.append(sample)

            sample_hyper_popu_list.append(hyper_popu_dict[sample])
            sample_detail_popu_list.append(detail_popu_dict[sample])


            # # classify rules
            if spec_chrom == "chr1":

                if ratio_dict["S1C1/5/19H1L.1-3(_6/4_5){n}-6"] > -1.01 * ratio_dict["S1C1/5/19H1L.1-5(_6/4_5){n}-6"] + 0.83:
                    feature_dict["cluster_1"].append(0)
                else:
                    feature_dict["cluster_1"].append(1)


            elif spec_chrom == "chr5":

                if ratio_dict["S1C1/5/19H1L.1_2/4_5_2/6"] < 0.3:
                    feature_dict["cluster_1"].append(0)
                else:
                    if ratio_dict["S1C1/5/19H1L.1_2/4_5_2/6"] > 1.3 * ratio_dict["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6"] + 0.11:
                        feature_dict["cluster_1"].append(2)
                    else:
                        feature_dict["cluster_1"].append(1)


            elif spec_chrom == "chr19":

                if ratio_dict["S1C1/5/19H1L.5-6"] > 0.18:
                    feature_dict["cluster_1"].append(0)
                elif ratio_dict["S1C1/5/19H1L.5-6"] > 0.04:
                    if ratio_dict["S1C1/5/19H1L.1-3(_6/4_5){n}_2/6"] < 0.1:
                        feature_dict["cluster_1"].append(-1)
                    else:
                        feature_dict["cluster_1"].append(1)
                else:
                    if ratio_dict["S1C1/5/19H1L.1_2/6"] > 0.09:
                        feature_dict["cluster_1"].append(2)
                    else:
                        if ratio_dict["S1C1/5/19H1L.1-3(_6/4_5){n}_2/6"] > 0.1:
                            feature_dict["cluster_1"].append(3)
                        else:
                            if ratio_dict["S1C1/5/19H1L.5(_6/4_5){n}-6"] > 0.2:
                                feature_dict["cluster_1"].append(4)
                            else:
                                feature_dict["cluster_1"].append(5)

            elif spec_chrom == "chr2":

                if ratio_dict["S2C2H1L.1-2"] > -0.41 * ratio_dict["S2C2H1L.1-4"] + 0.373:
                    feature_dict["cluster_1"].append(0)
                else:
                    if ratio_dict["S2C2H1L.1/3_4"] >= 0.06:
                        feature_dict["cluster_1"].append(1)
                    else:
                        feature_dict["cluster_1"].append(2)

            elif spec_chrom == "chr3":

                feature_dict["cluster_1"].append(0)


            elif spec_chrom == "chr4":

                if ratio_dict["S2C4H1L.1-2_7-19"] < 0.02:
                        feature_dict["cluster_1"].append(0)

                else:
                    if ratio_dict["S2C4H1L.1-2_7-19"] < -0.6 * ratio_dict["S2C4H1L.1-19"] + 0.5:
                        feature_dict["cluster_1"].append(1)
                    else:
                        feature_dict["cluster_1"].append(2)

            elif spec_chrom == "chr6":
                if ratio_dict["S01C6H1L.1-12_16-18"] >= 0.05:
                        feature_dict["cluster_1"].append(1)
                else:


                    if ratio_dict["S01C6H1L.1-18"] < 0.7:
                        feature_dict["cluster_1"].append(-1)
                    else:
                        feature_dict["cluster_1"].append(0)


            elif spec_chrom == "chr7":
                if ratio_dict["S1C7H1L.1-5_2-4_3-6"] > 0.1:
                    feature_dict["cluster_1"].append(1)
                else:
                    feature_dict["cluster_1"].append(0)


            elif spec_chrom == "chr8":

                if ratio_dict["S2C8H1L.1-7"] >= -0.6 * ratio_dict["S2C8H1L.1-11"] + 0.45:
                    if ratio_dict["S2C8H1L.1-11"] < 0.23:
                        feature_dict["cluster_1"].append(2)
                    else:
                        feature_dict["cluster_1"].append(1)
                else:
                    feature_dict["cluster_1"].append(0)

            elif spec_chrom == "chr9":
                if ratio_dict["S2C9H1L.1-7"] < 0.1:
                    feature_dict["cluster_1"].append(-1)
                else:
                    if ratio_dict["S2C9H1L.1-2_6-7"] >= ratio_dict["S2C9H1L.1-7"]:
                        feature_dict["cluster_1"].append(1)
                    else:
                        feature_dict["cluster_1"].append(0)

            elif spec_chrom == "chr10":

                if ratio_dict["S1C10H1L.1-2_5-8"] < 0.1:
                    feature_dict["cluster_1"].append(0)
                else:
                    if ratio_dict["S1C10H1L.1-2_5_6/8"] < -0.48 * ratio_dict["S1C10H1L.1-2_5-8"] + 0.225:
                        feature_dict["cluster_1"].append(1)
                    else:
                        if ratio_dict["S1C10H1L.1-2_5_6/8"] < -0.2 * ratio_dict["S1C10H1L.1-2_5-8"] + 0.145:
                            feature_dict["cluster_1"].append(2)
                        else:
                            feature_dict["cluster_1"].append(3)

            elif spec_chrom == "chr11":
                # # done
                if ratio_dict["S3C11H1L.1"] > 0.05:
                    feature_dict["cluster_1"].append(1)
                else:
                    feature_dict["cluster_1"].append(0)

            elif spec_chrom == "chr12":

                if ratio_dict["S1C12H1L.1_6-8"] < 0.06:
                    feature_dict["cluster_1"].append(0)
                else:
                    feature_dict["cluster_1"].append(1)


            elif spec_chrom == "chr13":

                if ratio_dict["S2C13/21H1L.1-11"] > 0.795:
                    feature_dict["cluster_1"].append(0)
                elif ratio_dict["S2C13/21H1L.1-11"] > 0.577:
                    feature_dict["cluster_1"].append(1)

                elif ratio_dict["S2C13/21H1L.1-11"] > 0.48:
                    feature_dict["cluster_1"].append(2)
                else:
                    if ratio_dict["S2C13/21H1L.1-7"] > 0.3:
                        feature_dict["cluster_1"].append(3)
                    else:
                        feature_dict["cluster_1"].append(4)

            elif spec_chrom == "chr14":
                if ratio_dict["S2C14/22H1L.1-8"] < 0.8:
                    feature_dict["cluster_1"].append(-1)
                else:
                    feature_dict["cluster_1"].append(0)

            elif spec_chrom == "chr15":
                feature_dict["cluster_1"].append(0)


            elif spec_chrom == "chr16":
                # # done
                if ratio_dict["S1C16H1L.1-10"] == 0:
                    feature_dict["cluster_1"].append(-1)
                else:

                    if ratio_dict["S1C16H1L.2-9"] < 0.2:
                        feature_dict["cluster_1"].append(0)
                    else:
                        feature_dict["cluster_1"].append(1)

            elif spec_chrom == "chr17":
                # # done
                if ratio_dict["S3C17H1L.1-16"] > 0.4:
                    if ratio_dict["S3C17H1L.1-13_15_15-16"] > 0.05:
                        feature_dict["cluster_1"].append(-1)
                    else:
                        feature_dict["cluster_1"].append(0)
                elif ratio_dict["S3C17H1L.1-16"] > 0.03:
                    feature_dict["cluster_1"].append(1)
                else:
                    feature_dict["cluster_1"].append(2)


            elif spec_chrom == "chr18":

                feature_dict["cluster_1"].append(0)

            elif spec_chrom == "chr20":
                if ratio_dict["S2C20H1L.7-16"] > 0.02:
                    if ratio_dict["S2C20H1L.1-16"] < 0.77:

                        if ratio_dict["S2C20H1L.7-16"] > -0.57 * ratio_dict["S2C20H1L.1-16"] +0.5:
                            feature_dict["cluster_1"].append(3)
                        else:
                            feature_dict["cluster_1"].append(4)
                    else:
                        feature_dict["cluster_1"].append(2)
                else:
                    if ratio_dict["S2C20H1L.1-16"] < 0.46:
                        feature_dict["cluster_1"].append(1)
                    else:
                        feature_dict["cluster_1"].append(0)


            elif spec_chrom == "chr21":

                if ratio_dict["S2C13/21H1L.1-11"] > 0.4:

                    if ratio_dict["S2C13/21H1L.8-11"] > 0.06:
                        feature_dict["cluster_1"].append(1)
                    else:
                        feature_dict["cluster_1"].append(0)

                else:
                    if ratio_dict["S2C13/21H1L.1-5"] > 0.3:
                        feature_dict["cluster_1"].append(3)
                    else:
                        feature_dict["cluster_1"].append(2)

            elif spec_chrom == "chr22":

                feature_dict["cluster_1"].append(0)

            elif spec_chrom == "chrX":
                feature_dict["cluster_1"].append(0)

            elif spec_chrom == "chrY":
                feature_dict["cluster_1"].append(0)

            else:
                feature_dict["cluster_1"].append(-1)

        # # calculate the ratio of -1 cluster
        total_cnt, outlier_cnt = 0, 0
        for cluster in feature_dict["cluster_1"]:
            total_cnt += 1

            if cluster == -1:
                outlier_cnt += 1

        print(spec_chrom, "c1", round((total_cnt - outlier_cnt) / total_cnt, 3), outlier_cnt, total_cnt, sep="\t")


        # # plot the res for this chromosome
        feature_names = list(feature_dict.keys())

        # order = np.argsort(feature_dict["aSat_length"])[::-1]
        order = np.argsort(feature_dict[feature_names[0]])[::-1]
        feature_sorted = {k: np.asarray(v)[order] for k, v in feature_dict.items()}
        sample_list_sorted = [sample_list[i] for i in order]
        sample_hyper_popu_list_sorted = [sample_hyper_popu_list[i] for i in order]
        sample_detail_popu_list_sorted = [sample_detail_popu_list[i] for i in order]

        # # output to file
        feature_file = os.path.join(output_folder, "{}.features.c1.txt".format(spec_chrom))
        with open(feature_file, "w") as fout:
            fout.write("Sample\thyper_popu\tdetail_popu\t" + "\t".join(feature_sorted.keys()) + "\n")

            for sample, hyper_popu, detail_popu, values in zip(sample_list_sorted, sample_hyper_popu_list_sorted, sample_detail_popu_list_sorted, zip(*feature_sorted.values())):
                fout.write(sample + "\t" + hyper_popu + "\t" + detail_popu + "\t" + "\t".join(map(str, values)) + "\n")

        # # # plot img
        # for oSat in oSat_list:
        #     feature_names.remove("{}_length".format(oSat))
        # fig, axes = plt.subplots(len(feature_names), 1, figsize=(12, 2.5 * len(feature_names)), sharex=True)
        #
        # x = np.arange(sample_cnt)
        # for ax, feature in zip(axes, feature_names):
        #     # # plot by all
        #     # ax.scatter(x, feature_sorted[feature], s=5, alpha=0.8)
        #     # ax.set_ylabel(feature)
        #
        #     # # plot by popu
        #     for popu, color in detail_popu_color.items():
        #
        #         idx = [i for i, p in enumerate(sample_detail_popu_list_sorted) if p == popu]
        #         if len(idx) != 0:
        #             ax.scatter(x[idx], np.asarray(feature_sorted[feature])[idx], color=color, s=5, label=popu)
        #             ax.set_title(feature)
        #
        # handles, labels = axes[0].get_legend_handles_labels()
        # fig.legend(handles, labels, loc="upper right", frameon=False)
        #
        # plt.tight_layout()
        # plt.savefig(os.path.join(output_folder, "{}.features.png".format(spec_chrom)), dpi=600)
        # plt.close()

        # #
        plot_principal_hor_analysis(feature_file, spec_chrom, "cluster_1")
        plot_feature_compare(feature_file, "cluster_1")


def principal_hor_analysis_classify_c2(chrom_order, output_folder):

    for spec_chrom in chrom_order:

        chrom_c1_file = os.path.join(output_folder, "{}.features.c1.txt".format(spec_chrom))
        chrom_c2_file = os.path.join(output_folder, "{}.features.c2.txt".format(spec_chrom))

        df = pd.read_csv(chrom_c1_file, sep="\t")


        # # initiate
        df["cluster_2"] = df["cluster_1"]
        df.loc[df["cluster_1"] == -1, "cluster_2"] = -1     # # maintain -1 if cluster_1 is -1

        # # deal with each cluster_1
        for cluster, features in c2_feature_dict[spec_chrom].items():

            cluster_mask = df["cluster_1"] == cluster
            if not cluster_mask.any():
                continue

            for feature in features:

                # # calculate boxplox
                feature_values = df.loc[cluster_mask, feature]

                Q1 = feature_values.quantile(0.25)
                Q3 = feature_values.quantile(0.75)
                IQR = Q3 - Q1

                lower = Q1 - 1.5 * IQR
                upper = Q3 + 1.5 * IQR

                # # find outlier and set cluster_2 to -1
                outlier_mask = (cluster_mask & ((df[feature] < lower) | (df[feature] > upper)))

                df.loc[outlier_mask, "cluster_2"] = -1


        # # stats and output
        outlier_cnt = (df["cluster_2"] == -1).sum()
        total_cnt = len(df)

        print(spec_chrom, "c2", round((total_cnt - outlier_cnt) / total_cnt, 3), outlier_cnt, total_cnt, sep="\t")
        df.to_csv(chrom_c2_file, sep="\t", index=False)

        plot_principal_hor_analysis(chrom_c2_file, spec_chrom, "cluster_2")
        # plot_feature_compare(chrom_c2_file, "cluster_2")


def principal_hor_analysis_classify_c3_c4_c5_c6(chrom_order, output_folder, level=5):

    for spec_chrom in chrom_order:

        chrom_c2_file = os.path.join(output_folder, "{}.features.c2.txt".format(spec_chrom))
        chrom_c3_c4_c5_c6_file = os.path.join(output_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom))

        df = pd.read_csv(chrom_c2_file, sep="\t")

        # # initiate
        df["cluster_3"] = df["cluster_2"]
        df.loc[df["cluster_2"] == -1, "cluster_3"] = -1     # # maintain -1 if cluster_2 is -1

        df["cluster_4"] = df["cluster_2"]
        df.loc[df["cluster_2"] == -1, "cluster_4"] = -1     # # maintain -1 if cluster_2 is -1

        df["cluster_5"] = df["cluster_2"]
        df.loc[df["cluster_2"] == -1, "cluster_5"] = -1     # # maintain -1 if cluster_2 is -1

        df["cluster_6"] = df["cluster_2"]
        df.loc[df["cluster_2"] == -1, "cluster_6"] = -1     # # maintain -1 if cluster_2 is -1

        # # deal with each cluster_2
        for c2_cluster, features in c3_c4_c5_c6_feature_dict[spec_chrom].items():

            c2_cluster_mask = df["cluster_2"] == c2_cluster
            if not c2_cluster_mask.any():
                continue

            feature_cnt = 0
            for feature in features:
                if feature_cnt == 0:
                    c3_c4_cluster_criteria = "HOR"
                    c3_c4_cluster = "cluster_3"

                    c5_c6_cluster_criteria = "HOR"
                    c5_c6_cluster = "cluster_5"

                else:
                    c3_c4_cluster_criteria = "Length"
                    c3_c4_cluster = "cluster_4"

                    c5_c6_cluster_criteria = "Length"
                    c5_c6_cluster = "cluster_6"

                feature_cnt += 1

                # # STEP: For C3_C4, calculate boxplox
                feature_values = df.loc[c2_cluster_mask, feature]
                Q1 = feature_values.quantile(0.25)
                Q3 = feature_values.quantile(0.75)
                IQR = Q3 - Q1

                lower = Q1 - 1.5 * IQR
                upper = Q3 + 1.5 * IQR
                median = feature_values.quantile(0.50)

                # # 1. > Upper
                c3_c4_cluster_mask = (c2_cluster_mask & (df[feature] > upper))
                df.loc[c3_c4_cluster_mask, c3_c4_cluster] =  f"{c2_cluster}_{c3_c4_cluster_criteria}_Larger_Upper"

                # # 2. Upper  ~ Q3
                c3_c4_cluster_mask = (c2_cluster_mask & (df[feature] <= upper) & (df[feature] > Q3))
                df.loc[c3_c4_cluster_mask, c3_c4_cluster] = f"{c2_cluster}_{c3_c4_cluster_criteria}_Upper_Q3"

                if level == 5:
                    # # 3. Upper Q3 ~ Q1
                    c3_c4_cluster_mask = (c2_cluster_mask & (df[feature] <= Q3) & (df[feature] >= Q1))
                    df.loc[c3_c4_cluster_mask, c3_c4_cluster] = f"{c2_cluster}_{c3_c4_cluster_criteria}_Q3_Q1"

                elif level == 6:
                    # # 3_1. Upper Q3 ~ Median
                    c3_c4_cluster_mask = (c2_cluster_mask & (df[feature] <= Q3) & (df[feature] > median))
                    df.loc[c3_c4_cluster_mask, c3_c4_cluster] = f"{c2_cluster}_{c3_c4_cluster_criteria}_Q3_Median"

                    # # 3_2. Upper Median ~ Q1
                    c3_c4_cluster_mask = (c2_cluster_mask & (df[feature] <= median) & (df[feature] >= Q1))
                    df.loc[c3_c4_cluster_mask, c3_c4_cluster] = f"{c2_cluster}_{c3_c4_cluster_criteria}_Median_Q1"

                else:
                    print("No this level", level)
                    exit()

                # # 4. Upper Q1 ~ Lower
                c3_c4_cluster_mask = (c2_cluster_mask & (df[feature] < Q1) & (df[feature] >= lower))
                df.loc[c3_c4_cluster_mask, c3_c4_cluster] = f"{c2_cluster}_{c3_c4_cluster_criteria}_Q1_Lower"

                # # 5. < Lower
                c3_c4_cluster_mask = (c2_cluster_mask & (df[feature] < lower))
                df.loc[c3_c4_cluster_mask, c3_c4_cluster] =  f"{c2_cluster}_{c3_c4_cluster_criteria}_Less_Lower"
                ######################################################################################################


                # # STEP: For C3_C4, calculate boxplox
                feature_values = df.loc[c2_cluster_mask, feature]

                Q_20 = feature_values.quantile(0.20)
                Q_40 = feature_values.quantile(0.40)
                Q_60 = feature_values.quantile(0.60)
                Q_80 = feature_values.quantile(0.80)

                c5_c6_cluster_mask = (c2_cluster_mask & (df[feature] >= Q_80))
                df.loc[c5_c6_cluster_mask, c5_c6_cluster] = f"{c2_cluster}_{c5_c6_cluster_criteria}_Q100_Q80"

                c5_c6_cluster_mask = (c2_cluster_mask & (df[feature] < Q_80) & (df[feature] > Q_60))
                df.loc[c5_c6_cluster_mask, c5_c6_cluster] = f"{c2_cluster}_{c5_c6_cluster_criteria}_Q80_Q60"

                c5_c6_cluster_mask = (c2_cluster_mask & (df[feature] <= Q_60) & (df[feature] > Q_40))
                df.loc[c5_c6_cluster_mask, c5_c6_cluster] = f"{c2_cluster}_{c5_c6_cluster_criteria}_Q60_Q40"

                c5_c6_cluster_mask = (c2_cluster_mask & (df[feature] <= Q_40) & (df[feature] > Q_20))
                df.loc[c5_c6_cluster_mask, c5_c6_cluster] = f"{c2_cluster}_{c5_c6_cluster_criteria}_Q40_Q20"

                c5_c6_cluster_mask = (c2_cluster_mask & (df[feature] <= Q_20))
                df.loc[c5_c6_cluster_mask, c5_c6_cluster] = f"{c2_cluster}_{c5_c6_cluster_criteria}_Q20_Q0"
                ######################################################################################################


        df.to_csv(chrom_c3_c4_c5_c6_file, sep="\t", index=False)


def significance_label(p):
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def plot_bin_level_hor_frequence_among_cluster(bin_count=20, spec_cluster="cluster_2"):

    """
    Similar with plot_bin_level_hor_frequence_among_popu() from step0_principle_hor_frequency.py
    """

    color_dict = {

        "HOR0": "#6688A6",  # 雾霾蓝
        "HOR1": "#D18B72",  # 陶土橙
        "HOR2": "#7FA18A",  # 鼠尾草绿
        "HOR3": "#A18BB5",  # 灰紫色
        "HOR4": "#C3A45D",  # 暗金色
        "HOR5": "#C47F91",  # 柔和玫瑰红
        "HOROther": "#858B91",  # 石板灰

    }

    bin_level_output_folder = os.path.join(work_folder, "bin_{}_principal_hor_frequency_{}".format(bin_count, spec_cluster))
    if not os.path.exists(bin_level_output_folder):
        os.mkdir(bin_level_output_folder)

    for spec_chrom in chrom_order:
        # if spec_chrom in ["chr1", "chr3", "chr14", "chr15", "chr18", "chr22", "chrX"]:
        #     continue

        chrom_input_folder = os.path.join(work_folder, "bin_{}_principal_hor_frequency/{}".format(bin_count, spec_chrom))
        if not os.path.exists(chrom_input_folder):
            continue

        chrom_color_dict = {}
        for i in range(len(c1_principle_hor_dict[spec_chrom])):
            chrom_color_dict[c1_principle_hor_dict[spec_chrom][i]] =  color_dict["HOR{}".format(i)]

        chrom_principle_hor_list = c1_principle_hor_dict[spec_chrom]

        chrom_features = chrom_principle_hor_list + ["Other"]
        bin_columns = ["bin", "bin_start", "bin_end", "bin_length"]

        # # load sample cluster info
        classify_file = os.path.join(pha_output_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom))

        cluster_dict = {}
        all_possible_clusters = []
        with open(classify_file) as fin:
            for line in fin:
                if "Sample" in line:
                    continue
                line_split = line.strip().split()

                sample = line_split[0]
                if spec_cluster == "cluster_2":
                    cluster = line_split[-5]
                elif spec_cluster == "cluster_3":
                    cluster = line_split[-4]
                elif spec_cluster == "cluster_4":
                    cluster = line_split[-3]
                elif spec_cluster == "cluster_5":
                    cluster = line_split[-2]
                elif spec_cluster == "cluster_6":
                    cluster = line_split[-1]
                else:
                    print("no this spec cluster ", spec_cluster)
                    exit()

                if cluster == "-1":
                    continue

                if cluster not in all_possible_clusters:
                    all_possible_clusters.append(cluster)
                cluster_dict[sample] = cluster


        # # collect from all samples
        all_df = []
        for sample, cluster in cluster_dict.items():

            sample_bin_file = os.path.join(chrom_input_folder, "{}.{}.bin_{}.txt".format(cluster, sample, bin_count))

            df = pd.read_csv(sample_bin_file, sep="\t", header=None)
            df.columns = bin_columns + chrom_features

            records = []
            for _, row in df.iterrows():
                for feature in chrom_features:
                    value = str(row[feature]).strip().split(":")[1]

                    count, freq = value.split(",")

                    records.append({"Sample": sample, "Cluster": cluster, "Bin": row["bin"], "Feature": feature, "Count": int(count), "Frequency": float(freq)})

            all_df.append(pd.DataFrame(records))

        all_df = pd.concat(all_df, ignore_index=True)
        print(all_df)
        all_df.to_csv(os.path.join(bin_level_output_folder, "{}.data.txt".format(spec_chrom)), sep="\t")

        # # STEP: plot a boxplot
        for cluster, cluster_df in all_df.groupby("Cluster"):
            cluster_df = cluster_df[cluster_df["Feature"] != "Other"].copy()

            plt.figure(figsize=(16, 7))

            # sns.boxplot(data=cluster_df, x="Bin", y="Frequency", hue="Feature", width=0.8, flierprops={"marker": "o", "markersize": 0.5, "linestyle": "none"}, )
            sns.boxplot(data=cluster_df, x="Bin", y="Frequency", hue="Feature", palette=chrom_color_dict, width=0.8, showfliers=False,    saturation=1,)

            plt.legend(loc="upper left", frameon=False)

            plt.tight_layout()
            plt.savefig(os.path.join(bin_level_output_folder, "{}.{}.boxplot.png".format(spec_chrom, cluster)), dpi=600)
            plt.savefig(os.path.join(bin_level_output_folder, "{}.{}.boxplot.svg".format(spec_chrom, cluster)), dpi=600)

            plt.close()

        # # STEP: calculate a p value for each feature in each bin
        if len(all_possible_clusters) > 1:
            stats_results = []
            for (bin_id, feature), group in all_df.groupby(["Bin", "Feature"]):
                clusters = sorted(group["Cluster"].dropna().unique())

                # # output for each cluster
                for cluster in clusters:
                    values = group.loc[group["Cluster"] == cluster, "Frequency"].dropna().values

                    q1 = np.percentile(values, 25)
                    median = np.percentile(values, 50)
                    q3 = np.percentile(values, 75)
                    iqr = round(q3 - q1, 3)
                    lower_bound = round(q1 - 1.5 * iqr, 3)
                    upper_bound = round(q3 + 1.5 * iqr, 3)

                    print(spec_chrom, bin_id, feature, cluster, median, q1, q3, iqr, lower_bound, upper_bound, sep="\t")

                # # compare between clusters
                for cluster1, cluster2 in combinations(clusters, 2):
                    values1 = group.loc[group["Cluster"] == cluster1, "Frequency"].dropna().values
                    values2 = group.loc[group["Cluster"] == cluster2, "Frequency"].dropna().values
                    if len(values1) == 0 or len(values2) == 0:
                        continue

                    values1_median = np.percentile(values1, 50)
                    values2_median = np.percentile(values2, 50)

                    stat, pvalue = mannwhitneyu(values1, values2, alternative="two-sided")

                    stats_results.append({"Bin": bin_id, "Feature": feature, "Cluster1": cluster1, "Cluster2": cluster2,
                                          "N1": len(values1), "N2": len(values2), "Median1": values1_median, "Median2": values2_median,
                                          "U": stat, "P_value": pvalue})

            stats_results_df = pd.DataFrame(stats_results)

            # FDR correction
            # stats_results_df["FDR"] = multipletests(stats_results_df["P_value"], method="fdr_bh" )[1]
            stats_results_df["FDR"] = (stats_results_df.groupby(["Bin", "Feature"])["P_value"].transform(
                lambda x: multipletests(x, method="fdr_bh")[1]))
            stats_results_df["Significance"] = stats_results_df["FDR"].apply(significance_label)
            print(stats_results_df)

            stats_results_df.to_csv(os.path.join(bin_level_output_folder, "{}.stats.txt".format(spec_chrom)), sep="\t", float_format="%.2e", index=False)




def stats_popu_freq_among_clusters(popu_mode="hyper"):

    output_freq_file = open(os.path.join(pha_output_folder, "chrall.popu_freq.{}.txt".format(popu_mode)), "w")
    output_cnt_file = open(os.path.join(pha_output_folder, "chrall.popu_cnt.{}.txt".format(popu_mode)), "w")

    if popu_mode == "hyper":
        popu_order = hyper_popu_order
    else:
        popu_order = detail_popu_order

    all_chrom_cluster_popu_cnt_dict = {}
    for spec_chrom in chrom_order:
        chrom_cluster_file = os.path.join(pha_output_folder, "{}.features.c2.txt".format(spec_chrom))

        cluster_popu_cnt_dict = {}
        with open(chrom_cluster_file) as fin:
            for line in fin:
                if "Sample" in line:
                    continue

                line_split = line.strip().split()
                if popu_mode == "hyper":
                    sample, sample_popu, sample_cluster = line_split[0], line_split[1], line_split[-1]
                else:
                    sample, sample_popu, sample_cluster = line_split[0], line_split[2], line_split[-1]

                if sample_cluster == "-1":
                    continue

                if sample_cluster not in cluster_popu_cnt_dict:
                    cluster_popu_cnt_dict[sample_cluster] = {popu: 0 for popu in popu_order}

                # cluster_popu_cnt_dict[sample_cluster]["All"] += 1
                cluster_popu_cnt_dict[sample_cluster][sample_popu] += 1

        print(spec_chrom, cluster_popu_cnt_dict)

        cluster_popu_proportion_dict = {
            cluster: {
                pop: round(count / sum(cluster_popu_cnt_dict[c][pop] for c in cluster_popu_cnt_dict), 3)
                for pop, count in counts.items()
            }
            for cluster, counts in cluster_popu_cnt_dict.items()
        }
        print(cluster_popu_proportion_dict)

        for cluster in sorted(cluster_popu_cnt_dict.keys()):
            output_cnt_file.write("{}\t{}\t{}\n".format(spec_chrom, cluster, "\t".join([str(cluster_popu_cnt_dict[cluster][popu]) for popu in popu_order])))

        for cluster in sorted(cluster_popu_proportion_dict.keys()):
            output_freq_file.write("{}\t{}\t{}\n".format(spec_chrom, cluster, "\t".join([str(cluster_popu_proportion_dict[cluster][popu]) for popu in popu_order])))

        all_chrom_cluster_popu_cnt_dict[spec_chrom] = cluster_popu_cnt_dict


    # # do fisher test
    if popu_mode == "hyper":
        populations = hyper_popu_order
    else:
        populations = detail_popu_order

    results = []

    for chrom, clusters in all_chrom_cluster_popu_cnt_dict.items():
        cluster_ids = sorted(list(clusters.keys()))

        for cluster in cluster_ids:
            other_clusters = [
                c for c in cluster_ids if c != cluster
            ]

            for pop in populations:
                # 当前 cluster 在该民族中的数量
                a = sum(clusters[cluster][pop] for _ in [0])

                # 其他 cluster 在该民族中的数量
                b = sum(
                    clusters[c][pop] for c in other_clusters
                )

                # 当前 cluster 在其他民族中的数量
                c = sum(
                    clusters[cluster][p]
                    for p in populations if p != pop
                )

                # 其他 cluster 在其他民族中的数量
                d = sum(
                    clusters[cl][p]
                    for cl in other_clusters
                    for p in populations if p != pop
                )

                # Fisher's exact test
                odds_ratio, p_value = fisher_exact(
                    [[a, b], [c, d]],
                    alternative='two-sided'
                )

                results.append({
                    'Chromosome': chrom,
                    'Cluster': cluster,
                    'Population': pop,
                    'Cluster_count': a,
                    'Other_cluster_count': b,
                    'Odds_ratio': odds_ratio,
                    'P_value': p_value
                })

    df = pd.DataFrame(results)

    # 对所有染色体、cluster、population 的检验统一进行 BH-FDR
    df['FDR'] = multipletests(
        df['P_value'], method='fdr_bh'
    )[1]

    # 显著富集 / 显著缺失
    df['Enrichment'] = 'Not significant'
    df.loc[
        (df['FDR'] < 0.05) & (df['Odds_ratio'] > 1),
        'Enrichment'
    ] = 'Enriched'
    df.loc[
        (df['FDR'] < 0.05) & (df['Odds_ratio'] < 1),
        'Enrichment'
    ] = 'Depleted'

    print(df)

    df.to_csv(os.path.join(pha_output_folder, "chrall.popu_cnt.{}.fisher.txt".format(popu_mode)), sep="\t")



def stats_popu_freq_among_clusters_rmcpc():
    popu_mode = "detail"
    output_freq_file = open(os.path.join(pha_output_folder, "chrall.popu_freq.{}.rmcpc.txt".format(popu_mode)), "w")
    output_cnt_file = open(os.path.join(pha_output_folder, "chrall.popu_cnt.{}.rmcpc.txt".format(popu_mode)), "w")


    popu_order = ["AFR", "EUR", "AMR", "SAS", "EAS"]

    all_chrom_cluster_popu_cnt_dict = {}
    for spec_chrom in chrom_order:
        chrom_cluster_file = os.path.join(pha_output_folder, "{}.features.c2.txt".format(spec_chrom))

        cluster_popu_cnt_dict = {}
        with open(chrom_cluster_file) as fin:
            for line in fin:
                if "Sample" in line:
                    continue

                line_split = line.strip().split()

                sample, sample_popu, sample_cluster = line_split[0], line_split[2], line_split[-1]


                if sample_cluster == "-1":
                    continue

                if "CPC" in sample:
                    continue
                if sample_popu == "CHN":
                    sample_popu = "EAS"

                if sample_cluster not in cluster_popu_cnt_dict:
                    cluster_popu_cnt_dict[sample_cluster] = {popu: 0 for popu in popu_order}

                cluster_popu_cnt_dict[sample_cluster][sample_popu] += 1

        print(spec_chrom, cluster_popu_cnt_dict)

        cluster_popu_proportion_dict = {
            cluster: {
                pop: round(count / sum(cluster_popu_cnt_dict[c][pop] for c in cluster_popu_cnt_dict), 3)
                for pop, count in counts.items()
            }
            for cluster, counts in cluster_popu_cnt_dict.items()
        }
        print(cluster_popu_proportion_dict)

        for cluster in sorted(cluster_popu_cnt_dict.keys()):
            output_cnt_file.write("{}\t{}\t{}\n".format(spec_chrom, cluster, "\t".join([str(cluster_popu_cnt_dict[cluster][popu]) for popu in popu_order])))

        for cluster in sorted(cluster_popu_proportion_dict.keys()):
            output_freq_file.write("{}\t{}\t{}\n".format(spec_chrom, cluster, "\t".join([str(cluster_popu_proportion_dict[cluster][popu]) for popu in popu_order])))

        all_chrom_cluster_popu_cnt_dict[spec_chrom] = cluster_popu_cnt_dict


    # # do fisher test

    populations = popu_order

    results = []

    for chrom, clusters in all_chrom_cluster_popu_cnt_dict.items():
        cluster_ids = sorted(list(clusters.keys()))

        for cluster in cluster_ids:
            other_clusters = [
                c for c in cluster_ids if c != cluster
            ]

            for pop in populations:
                # 当前 cluster 在该民族中的数量
                a = sum(clusters[cluster][pop] for _ in [0])

                # 其他 cluster 在该民族中的数量
                b = sum(
                    clusters[c][pop] for c in other_clusters
                )

                # 当前 cluster 在其他民族中的数量
                c = sum(
                    clusters[cluster][p]
                    for p in populations if p != pop
                )

                # 其他 cluster 在其他民族中的数量
                d = sum(
                    clusters[cl][p]
                    for cl in other_clusters
                    for p in populations if p != pop
                )

                # Fisher's exact test
                odds_ratio, p_value = fisher_exact(
                    [[a, b], [c, d]],
                    alternative='two-sided'
                )

                results.append({
                    'Chromosome': chrom,
                    'Cluster': cluster,
                    'Population': pop,
                    'Cluster_count': a,
                    'Other_cluster_count': b,
                    'Odds_ratio': odds_ratio,
                    'P_value': p_value
                })

    df = pd.DataFrame(results)

    # 对所有染色体、cluster、population 的检验统一进行 BH-FDR
    df['FDR'] = multipletests(
        df['P_value'], method='fdr_bh'
    )[1]

    # 显著富集 / 显著缺失
    df['Enrichment'] = 'Not significant'
    df.loc[
        (df['FDR'] < 0.05) & (df['Odds_ratio'] > 1),
        'Enrichment'
    ] = 'Enriched'
    df.loc[
        (df['FDR'] < 0.05) & (df['Odds_ratio'] < 1),
        'Enrichment'
    ] = 'Depleted'

    print(df)

    df.to_csv(os.path.join(pha_output_folder, "chrall.popu_cnt.{}.fisher.rmcpc.txt".format(popu_mode)), sep="\t")


def generate_sample_HOR_freq_matrix():
    with open(os.path.join(work_folder, "stats.aSat_length_monomer.merge_0k.json")) as file:
            sample_sat_length_dict = dict(json.load(file))

    all_pricipal_HORs = c1_principle_hor_dict.values()

    sample_HOR_dict = {}
    for sample in sample_sat_length_dict:
        sample_HOR_dict[sample] = ["." for hor in all_pricipal_HORs]

    for i in range(len(chrom_order)):
        chrom = chrom_order[i]

        chrom_cluster_file = os.path.join(pha_output_folder, "{}.features.c2.txt".format(chrom))


def generate_sample_cluster_matrix():

    with open(os.path.join(work_folder, "stats.aSat_length_monomer.merge_0k.json")) as file:
            sample_sat_length_dict = dict(json.load(file))

    sample_cluster_dict = {}
    for sample in sample_sat_length_dict:
        sample_cluster_dict[sample] = ["." for chrom in chrom_order]

    for i in range(len(chrom_order)):
        chrom = chrom_order[i]

        chrom_cluster_file = os.path.join(pha_output_folder, "{}.features.c2.txt".format(chrom))
        with open(chrom_cluster_file) as fin:
            for line in fin:
                if "Sample" in line:
                    continue

                line_split = line.strip().split()

                sample, cluster = line_split[0], line_split[-1]
                sample_cluster_dict[sample][i] = cluster

    print("{}\t{}".format("Sample", "\t".join(chrom_order)))
    for sample in sample_cluster_dict:
        print("{}\t{}".format(sample, "\t".join(sample_cluster_dict[sample])))


if __name__ == '__main__':

    work_folder = "/mnt/e/workspace/cpc2_centro/res_all_post"
    stats_res_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot/chrall"

    pha_output_folder = os.path.join(work_folder, "principal_hor_analysis_outlier")
    if not os.path.exists(pha_output_folder):
        os.mkdir(pha_output_folder)

    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    # chrom_order = ["chrY"]


    with open(os.path.join(work_folder, "stats.whole_core_cent.json")) as file:
        whole_cent_dict = dict(json.load(file))

    with open(os.path.join(work_folder, "stats.aSat_length_monomer.merge_0k.json")) as file:
        sample_sSat_length_dict = dict(json.load(file))

    # oSat_list = ["HSat1", "HSat2", "HSat3", "BSat", "GSat"]
    oSat_list = []
    sample_oSat_length_dict = {}
    for oSat in oSat_list:
        with open(os.path.join(work_folder, "stats.{}_length_rpmask.merge_0k.json".format(oSat))) as file:
            sample_oSat_length_dict[oSat] = dict(json.load(file))

    # # load popu info
    detail_popu_dict = {}  # # what is detail: AFR, AMR, EUR, SAS, EAS-other, EAS-CPC
    detail_popu_order = ["AFR", "EUR", "AMR", "SAS", "EAS", "CHN"]
    detail_popu_color = {"AFR": "#156082", "EUR": "#4ea72e", "AMR": "#0f9ed5", "SAS": "#e97132", "EAS": "#70ad47", "CHN": "#a02b93"}

    hyper_popu_dict = {}  # # what is hyper: non-EAS (other), EAS-other, and EAS-CPC
    hyper_popu_order = ["Other", "CHN"]
    hyper_popu_color = {"Other": "#156082", "CHN": "#a02b93"}

    all_included_samples = []

    # with open("data_population_hgsvc3_hprc2_cpc2.txt") as fin:
    with open("../data/data_population_hgsvc3_hprc2_cpc2.txt") as fin:
        for line in fin:
            line_split = line.strip().split("\t")
            if line_split[6] == "NOUSE" or line_split[5] == "NOUSE":
                continue
            all_included_samples.append(line_split[1])
            detail_popu_dict[line_split[1]] = line_split[6]
            hyper_popu_dict[line_split[1]] = line_split[5]

    # # set hors
    c1_principle_hor_dict = {"chr1": ["S1C1/5/19H1L.1-5(_6/4_5){n}-6", "S1C1/5/19H1L.1-3(_6/4_5){n}-6"], "chr2": ["S2C2H1L.1-4", "S2C2H1L.1-2", "S2C2H1L.1_4", "S2C2H1L.1/3_4"],
                               "chr3": ["S01/1C3H1L.1-17", "S01/1C3H1L.1-4_12/5_6-17", "S01/1C3H1L.1-6_5/7_8-17", "S01/1C3H1L.1-4_5/7_8-17"], "chr4": ["S2C4H1L.1-19", "S2C4H1L.1-2_7-19", "S2C4H1L.1-9_3-4_12-19", "S2C4H1L.1-18", "S2C4H1L.1-3_16-19"],
                               "chr5": ["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1_2/4_5_2/6", "S1C1/5/19H1L.1-2_5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1-5_6/4"], "chr6": ["S01C6H1L.1-18", "S01C6H1L.1-12_16-18"],
                               "chr7": ["S1C7H1L.1-6", "S1C7H1L.1-5_2-4_3-6"], "chr8": ["S2C8H1L.1-11", "S2C8H1L.1-7", "S2C8H1L.1-3_4/7_8-11"],
                               "chr9": ["S2C9H1L.1-7", "S2C9H1L.1-2_6-7", "S2C9H1L.1-3_7"], "chr10": ["S1C10H1L.1-8", "S1C10H1L.1-2_5-8", "S1C10H1L.1-2_5_6/8"],
                               "chr11": ["S3C11H1L.1-5", "S3C11H1L.1"], "chr12": ["S1C12H1L.1-8", "S1C12H1L.1_6-8"],
                               "chr13": ["S2C13/21H1L.1-11", "S2C13/21H1L.8-11", "S2C13/21H1L.1-7", "S2C13/21H1L.1-2_9-11"], "chr14": ["S2C14/22H1L.1-8"],
                               "chr15": ["S2C15H1L.1-11", "S2C15H1L.1-4"], "chr16": ["S1C16H1L.1-10", "S1C16H1L.2-9"],
                               "chr17": ["S3C17H1L.1-16", "S3C17H1L.1-10_15_15-16", "S3C17H1L.1-13_15_15-16", "S3C17H1L.1-13_15-16"], "chr18": ["S2C18H1L.1-12", "S2C18H1L.1-2_3/7_8-12", "S2C18H1L.1-2_7_4-12"],
                               "chr19": ["S1C1/5/19H1L.5(_6/4_5){n}-6", "S1C1/5/19H1L.5-6", "S1C1/5/19H1L.1-3(_6/4_5){n}_2/6", "S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1_2/6", "S1C1/5/19H1L.1(_6/4_5){n}_2/6"], "chr20": ["S2C20H1L.1-16", "S2C20H1L.7-16"],
                               "chr21": ["S2C13/21H1L.1-11", "S2C13/21H1L.1-2_9-11", "S2C13/21H1L.1-5", "S2C13/21H1L.8-11"], "chr22": ["S2C14/22H1L.1-8"],
                               "chrX": ["S3CXH1L.1-12"], "chrY": ["S4CYH1L.1-15_26-31_32/34_35-46", "S4CYH1L.1-15_26-46"]
                               }



    c2_feature_dict = {
        # # 长度规则：如果长度有散点，则加入

        "chr1": {0: [], 1: ["aSat_length"]},
        "chr2": {0: ["S2C2H1L.1_4", "aSat_length"], 1: ["S2C2H1L.1_4", "S2C2H1L.1/3_4", "aSat_length"], 2: ["S2C2H1L.1_4", "S2C2H1L.1/3_4", "aSat_length"]},
        "chr3": {0: []},
        "chr4": {0: [], 1: [], 2: []},
        "chr5": {0: ["S1C1/5/19H1L.1-5_6/4"], 1: ["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1_2/4_5_2/6"], 2:["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "S1C1/5/19H1L.1-2_5(_6/4_5){n}_2/6"]},
        "chr6": {0: ["aSat_length"], 1: ["aSat_length"], 2: ["aSat_length"]},
        "chr7": {0: [], 1: []},
        "chr8": {0: ["S2C8H1L.1-3_4/7_8-11"], 1: ["S2C8H1L.1-3_4/7_8-11", "S2C8H1L.1-7"], 2: ["S2C8H1L.1-3_4/7_8-11"]},
        "chr9": {0: [], 1: []},
        "chr10": {0: ["S1C10H1L.1-2_5_6/8"], 1: ["S1C10H1L.1-8"], 2: ["S1C10H1L.1-2_5_6/8"], 3: ["S1C10H1L.1-2_5_6/8"]},
        "chr11": {0: ["aSat_length"], 1: ["aSat_length"]},
        "chr12": {0: ["aSat_length"], 1: ["aSat_length"]},
        "chr13": {0: ["S2C13/21H1L.8-11"], 1: ["S2C13/21H1L.8-11", "S2C13/21H1L.1-11"], 2: ["S2C13/21H1L.1-7"], 3: ["S2C13/21H1L.1-7"], 4: ["S2C13/21H1L.1-11"]},
        "chr14": {0: []},
        "chr15": {0: []},
        "chr16": {0: [], 1: []},
        "chr17": {0: ["S3C17H1L.1-13_15-16"], 1: [], 2:[]},
        "chr18": {0: []},
        "chr19": {0: [], 1: ["S1C1/5/19H1L.5-6"], 2: ["S1C1/5/19H1L.5(_6/4_5){n}-6"], 3: ["S1C1/5/19H1L.5(_6/4_5){n}-6"], 4: [], 5: ["S1C1/5/19H1L.5(_6/4_5){n}-6"]},
        "chr20": {0: ["S2C20H1L.7-16",  "aSat_length"], 1: ["aSat_length"], 2: ["S2C20H1L.1-16", "S2C20H1L.7-16", "aSat_length"], 3:["aSat_length"], 4:["S2C20H1L.1-16", "S2C20H1L.7-16", "aSat_length"]},
        "chr21": {0: ["S2C13/21H1L.1-11"], 1: ["S2C13/21H1L.1-2_9-11"], 2: ["S2C13/21H1L.1-5"], 3: ["S2C13/21H1L.1-5"]},
        "chr22": {0: []},
        "chrX": {0: []},
        "chrY": {0: []}

    }


    c3_c4_c5_c6_feature_dict = {
        "chr1": {0: ["S1C1/5/19H1L.1-5(_6/4_5){n}-6", "aSat_length"], 1: ["S1C1/5/19H1L.1-5(_6/4_5){n}-6", "aSat_length"]},
        "chr2": {0: ["S2C2H1L.1-4", "aSat_length"], 1: ["S2C2H1L.1-4", "aSat_length"], 2: ["S2C2H1L.1-4", "aSat_length"]},
        "chr3": {0: ["S01/1C3H1L.1-17", "aSat_length"]},
        "chr4": {0: ["S2C4H1L.1-19", "aSat_length"], 1: ["S2C4H1L.1-19", "aSat_length"], 2: ["S2C4H1L.1-19", "aSat_length"]},
        "chr5": {0: ["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "aSat_length"], 1: ["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "aSat_length"], 2: ["S1C1/5/19H1L.1-5(_6/4_5){n}_2/6", "aSat_length"]},
        "chr6": {0: ["S01C6H1L.1-18", "aSat_length"], 1: ["S01C6H1L.1-18", "aSat_length"], 2: ["S01C6H1L.1-18", "aSat_length"]},
        "chr7": {0: ["S1C7H1L.1-6", "aSat_length"], 1: ["S1C7H1L.1-6", "aSat_length"]},
        "chr8": {0: ["S2C8H1L.1-11", "aSat_length"], 1: ["S2C8H1L.1-11", "aSat_length"], 2: ["S2C8H1L.1-11", "aSat_length"]},
        "chr9": {0: ["S2C9H1L.1-7", "aSat_length"], 1: ["S2C9H1L.1-7", "aSat_length"]},
        "chr10": {0: ["S1C10H1L.1-8", "aSat_length"], 1: ["S1C10H1L.1-8", "aSat_length"], 2: ["S1C10H1L.1-8", "aSat_length"], 3: ["S1C10H1L.1-8", "aSat_length"]},
        "chr11": {0: ["S3C11H1L.1-5", "aSat_length"], 1: ["S3C11H1L.1-5", "aSat_length"]},
        "chr12": {0: ["S1C12H1L.1-8", "aSat_length"], 1: ["S1C12H1L.1-8", "aSat_length"]},
        "chr13": {0: ["S2C13/21H1L.1-11", "aSat_length"], 1: ["S2C13/21H1L.1-11", "aSat_length"], 2: ["S2C13/21H1L.1-11", "aSat_length"], 3: ["S2C13/21H1L.1-11", "aSat_length"], 4: ["S2C13/21H1L.1-11", "aSat_length"]},
        "chr14": {0: ["S2C14/22H1L.1-8", "aSat_length"]},
        "chr15": {0: ["S2C15H1L.1-11", "aSat_length"], 1: ["S2C15H1L.1-11", "aSat_length"]},
        "chr16": {0: ["S1C16H1L.1-10", "aSat_length"], 1: ["S1C16H1L.1-10", "aSat_length"]},
        "chr17": {0: ["S3C17H1L.1-16", "aSat_length"], 1: ["S3C17H1L.1-16", "aSat_length"], 2: ["S3C17H1L.1-16", "aSat_length"]},
        "chr18": {0: ["S2C18H1L.1-12", "aSat_length"]},
        "chr19": {0: ["S1C1/5/19H1L.5(_6/4_5){n}-6", "aSat_length"], 1: ["S1C1/5/19H1L.5(_6/4_5){n}-6", "aSat_length"], 2: ["S1C1/5/19H1L.5(_6/4_5){n}-6", "aSat_length"], 3: ["S1C1/5/19H1L.5(_6/4_5){n}-6", "aSat_length"], 4: ["S1C1/5/19H1L.5(_6/4_5){n}-6", "aSat_length"], 5: ["S1C1/5/19H1L.5(_6/4_5){n}-6", "aSat_length"]},
        "chr20": {0: ["S2C20H1L.1-16", "aSat_length"], 1: ["S2C20H1L.1-16", "aSat_length"], 2: ["S2C20H1L.1-16", "aSat_length"], 3: ["S2C20H1L.1-16", "aSat_length"], 4: ["S2C20H1L.1-16", "aSat_length"]},
        "chr21": {0: ["S2C13/21H1L.1-11", "aSat_length"], 1: ["S2C13/21H1L.1-11", "aSat_length"], 2: ["S2C13/21H1L.1-11", "aSat_length"], 3: ["S2C13/21H1L.1-11", "aSat_length"]},
        "chr22": {0: ["S2C14/22H1L.1-8", "aSat_length"]},
        "chrX": {0: ["S3CXH1L.1-12", "aSat_length"]},
        "chrY": {0: ["S4CYH1L.1-15_26-31_32/34_35-46", "aSat_length"]}

    }


    # # STEP: classify
    principal_hor_analysis_classify_c1(chrom_order, pha_output_folder)
    principal_hor_analysis_classify_c2(chrom_order, pha_output_folder)
    principal_hor_analysis_classify_c3_c4_c5_c6(chrom_order, pha_output_folder)

    # # STEP: stats
    # generate_sample_HOR_freq_matrix()
    # generate_sample_cluster_matrix()
    stats_popu_freq_among_clusters(popu_mode="detail")
    stats_popu_freq_among_clusters_rmcpc()

    # # STEP: plot difference among cluster in bin_level
    plot_bin_level_hor_frequence_among_cluster(bin_count=10, spec_cluster="cluster_2")



