import multiprocessing
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
from scipy.stats import ranksums
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.stats import mannwhitneyu
from matplotlib.ticker import PercentFormatter
from statsmodels.stats.multitest import multipletests
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import kruskal
import diptest
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
                        stats_dict[current_chr][hor_type] = int(hor_cnt)
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


def load_hor_file(hor_file):

    hor_list = []
    with open(hor_file) as fin:

        for line in fin:
            line_split = line.strip().split("\t")
            chrom, start, end, hor, orient = line_split[0], int(line_split[1]), int(line_split[2]), line_split[3], line_split[4]

            hor_list.append([chrom, start, end, hor, orient])

    return hor_list


def reverse_hor(in_hor):
    name = in_hor.split(".")[0]
    combine = in_hor.split(".")[1]

    # # reverse combine
    ids = re.split(r'[-_]', combine)
    seps = re.findall(r'[-_]', combine)

    ids.reverse()
    seps.reverse()

    return name + "." + ids[0] + ''.join(sep + id_ for sep, id_ in zip(seps, ids[1:]))



def find_principle_hor():
    total_hor_cnt = 0
    for spec_chrom in chrom_order:
        candidate_hors = {}

        total_sample_num = 0
        for sample_stats_file in os.listdir(stats_res_folder):
            if not sample_stats_file.endswith("#stv_stats.tsv"):
                continue

            sample = sample_stats_file.replace("#stv_stats.tsv", "")

            if sample not in all_included_samples:
                continue

            if "{}_RagTag".format(spec_chrom) not in whole_core_cent[sample]:
                continue

            sample_stats_file = os.path.join(stats_res_folder, sample_stats_file)

            try:
                sample_hor_stats_dict = load_hor_stats_file(sample_stats_file)[spec_chrom]
            except KeyError:
                continue

            for hor in sample_hor_stats_dict:
                if hor == "total":
                    continue

                hor_ratio = round(sample_hor_stats_dict[hor] / sample_hor_stats_dict["total"], 3)

                if hor_ratio > 0.1:
                    if hor not in candidate_hors:
                        candidate_hors[hor] = 0
                    candidate_hors[hor] += 1

            total_sample_num += 1

        final_hors = {}
        for hor in candidate_hors:
            if candidate_hors[hor] > 10:
                final_hors[hor] = candidate_hors[hor]
        print(spec_chrom, len(final_hors), final_hors)
        # print(spec_chrom, len(final_hors))

        total_hor_cnt += len(final_hors)

    print(total_hor_cnt)


def get_cnt_and_freq_for_one_HOR(spec_chrom, spec_HOR):

    hor_cnt = 0

    total_sample_num = 0
    for sample_stats_file in os.listdir(stats_res_folder):
        if not sample_stats_file.endswith("#stv_stats.tsv"):
            continue

        sample = sample_stats_file.replace("#stv_stats.tsv", "")

        if sample not in all_included_samples:
            continue

        if "{}_RagTag".format(spec_chrom) not in whole_core_cent[sample]:
            continue

        sample_stats_file = os.path.join(stats_res_folder, sample_stats_file)

        try:
            sample_hor_stats_dict = load_hor_stats_file(sample_stats_file)[spec_chrom]
        except KeyError:
            continue

        for hor in sample_hor_stats_dict:
            if hor != spec_HOR:
                continue

            hor_ratio = round(sample_hor_stats_dict[hor] / sample_hor_stats_dict["total"], 3)

            hor_cnt += sample_hor_stats_dict[hor]

        total_sample_num += 1


    print(spec_chrom, spec_HOR, hor_cnt, hor_cnt / total_sample_num)

def plot_hor_frequence_among_popu(popu_mode):

    if popu_mode == "hyper":
        popu_dict = hyper_popu_dict
        popu_order = hyper_popu_order

    else:
        popu_dict = detail_popu_dict
        popu_order = detail_popu_order

    output_folder = os.path.join(work_folder, "principal_hor_frequency_{}_popu".format(popu_mode))
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_stats_file = open(os.path.join(output_folder, "chrall.stats.txt"), "w")
    output_all_file = open(os.path.join(output_folder, "chrall.all.txt"), "w")

    # # traverse each sample
    for spec_chrom in chrom_order:
        collect_dict = {}
        collect_bins_dict = {}
        chrom_selected_hors = c1_principle_hor_dict[spec_chrom]

        # # initiate
        for hor in chrom_selected_hors:
            collect_dict[hor] = {popu: [] for popu in popu_order}
            # collect_bins_dict[hor] = {popu: [0 for i in range(10)] for popu in popu_order}
            collect_bins_dict[hor] = {popu: [0 for i in range(4)] for popu in popu_order}

        for sample_stats_file in os.listdir(stats_res_folder):
            if not sample_stats_file.endswith("#stv_stats.tsv"):
                continue

            sample = sample_stats_file.replace("#stv_stats.tsv", "")

            if sample not in popu_dict:
                continue

            sample_popu = popu_dict[sample]

            # # load hor stats res
            sample_stats_file = os.path.join(stats_res_folder, sample_stats_file)
            try:
                sample_hor_stats_dict = load_hor_stats_file(sample_stats_file)[spec_chrom]
            except KeyError:
                continue

            for hor in chrom_selected_hors:
                ratio = round(sample_hor_stats_dict.get(hor, 0) / sample_hor_stats_dict["total"], 3)
                collect_dict[hor][sample_popu].append(ratio)

                if ratio < 0.1:
                    bin_index = 0
                elif 0.1 <= ratio < 0.5:
                    bin_index = 1
                elif 0.5 <= ratio < 0.9:
                    bin_index = 2
                else:
                    bin_index = 3

                # bin_index = int(ratio * 10)

                collect_bins_dict[hor][sample_popu][bin_index] += 1

        for hor in collect_dict:

            hor_dict = collect_dict[hor]
            hor_bin_dict = collect_bins_dict[hor]

            # # for each hor, do stat test
            groups = {popu: np.asarray(lengths, dtype=float) for popu, lengths in hor_dict.items()}
            _, kruskal_p = kruskal(*groups.values())

            output_stats_file.write("{}\t{}, Kruskal-Wallis test, P={}\n".format(spec_chrom, hor, round(kruskal_p, 3) if kruskal_p > 0.01 else format(kruskal_p, '.1e')))

            for base_popu in popu_order:

                q1 = round(np.percentile(hor_dict[base_popu], 25), 4)
                median = round(np.percentile(hor_dict[base_popu], 50), 4)
                q3 = round(np.percentile(hor_dict[base_popu], 75), 4)
                iqr = q3 - q1
                lower = int(q1 - 1.5 * iqr)
                upper = int(q3 + 1.5 * iqr)

                dip, dip_p = diptest.diptest(np.asarray(hor_dict[base_popu]))
                dip_p_label = significance_label(dip_p)
                dip_p = round(dip_p, 3) if dip_p > 0.01 else format(dip_p, '.1e')
                #

                p_values = []
                for target_popu in popu_order:

                    stat, muw_p = mannwhitneyu(hor_dict[base_popu], hor_dict[target_popu], alternative="two-sided", method="auto")
                    muw_p_label = significance_label(muw_p)
                    muw_p = round(muw_p, 3) if muw_p > 0.01 else format(muw_p, '.1e')

                    p_values.append(str(muw_p))

                # output_stats_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(spec_chrom, hor, base_popu, len(hor_dict[base_popu]), median, dip_p, "\t".join(p_values)))
                output_stats_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(spec_chrom, hor, base_popu, len(hor_dict[base_popu]), "\t".join([str(round(var / len(hor_dict[base_popu]), 3)) for var in hor_bin_dict[base_popu]]), "\t".join(p_values)))

            # # output raw data
            for popu in hor_dict:
                output_all_file.write("{}\t{}\t{}\t{}\n".format(spec_chrom, hor, popu, "\t".join([str(var) for var in hor_dict[popu]])))

        # # plot
        records = []

        for hor, popu_info in collect_dict.items():
            for popu, values in popu_info.items():
                for value in values:
                    records.append({'hor': hor, 'popu': popu, 'value': value})

        df = pd.DataFrame(records)

        n = len(collect_dict)
        ncols = 5
        nrows = math.ceil(n / ncols)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4 * ncols, 4 * nrows))
        axes = axes.flatten()

        for i, hor in enumerate(collect_dict.keys()):
            ax = axes[i]

            sns.violinplot(data=df[df['hor'] == hor], x='popu', y='value', ax=ax)

            ax.set_title(hor)
            ax.set_xlabel('')

        # remove non-used subplot
        for i in range(n, len(axes)):
            fig.delaxes(axes[i])

        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, "{}.hor_freq.png".format(spec_chrom)), dpi=600)
        plt.close()

    output_stats_file.close()
    output_all_file.close()

def plot_bin_level_hor_frequence_each_sample(spec_chrom, bin_count=50):


    output_folder = os.path.join(work_folder, "bin_{}_principal_hor_frequency".format(bin_count))
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # # load cluste
    cluster_dict = {}
    classify_file = os.path.join(pha_input_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom))
    with open(classify_file) as fin:
        for line in fin:
            line_split = line.strip().split()
            sample, cluster = line_split[0], line_split[-5]

            cluster_dict[sample] = cluster

    chrom_output_folder = os.path.join(output_folder, spec_chrom)
    if not os.path.exists(chrom_output_folder):
        os.mkdir(chrom_output_folder)

    chrom_principle_hor_list = c1_principle_hor_dict[spec_chrom]
    chrom_annot_res_folder = os.path.join(annot_res_folder, spec_chrom)

    for sample_hor_file in os.listdir(chrom_annot_res_folder):
        if not sample_hor_file.endswith("#stv_hor_osat_te.pha.bed"):
            continue

        sample_name = sample_hor_file.replace("#stv_hor_osat_te.pha.bed", "")

        if sample_name == "T2T-CHM13.fa.gz":
            sample_cluster = "ref"

        else:
            if sample_name not in all_included_samples:
                continue
            if "{}_RagTag".format(spec_chrom) not in whole_core_cent[sample_name]:
                continue
            if sample_name not in cluster_dict:
                continue

            sample_cluster = cluster_dict[sample_name]

        sample_hor_list = load_hor_file(os.path.join(chrom_annot_res_folder, sample_hor_file))
        if len(sample_hor_list) == 0:
            print("Empty HOR for sample", sample_name, spec_chrom)
            continue


        # # STEP: begin to stats in each bin
        # # init
        sample_bin_list = []
        for bin_index in range(bin_count):
            sample_bin_list.append({hor: 0 for hor in chrom_principle_hor_list + ["Other"]})

        # # for each hor, check and add into the bin
        # ref_start, ref_end = int(sample_hor_list[0][0].split(":")[1].split("-")[0]), int(sample_hor_list[0][0].split(":")[1].split("-")[1])
        ref_start, ref_end = int(sample_hor_list[0][1]), int(sample_hor_list[-1][2])

        bin_size = int((ref_end - ref_start + 1) / bin_count)
        # print(spec_chrom, sample_name, ref_start, ref_end, bin_size)

        for _, hor_start, hor_end, hor_type, hor_orient in sample_hor_list:

            if hor_type in TE_types or hor_type in oSat_types:
                continue

            if hor_orient == "-":
                hor_type = reverse_hor(hor_type)
            hor_type = re.sub(r'\{\d+\}', '{n}', hor_type)

            hor_bin_index = int((hor_start - ref_start) / bin_size)

            if hor_type in chrom_principle_hor_list:
                sample_bin_list[hor_bin_index][hor_type] += 1
            else:
                sample_bin_list[hor_bin_index]["Other"] += 1


        # # output into file
        sample_bin_freq_list = []
        sample_bin_output_file = os.path.join(chrom_output_folder, "{}.{}.bin_{}.txt".format(sample_cluster, sample_name, bin_count))
        with open(sample_bin_output_file, "w") as fout:
            for bin_index in range(len(sample_bin_list)):
                bin_start, bin_end = int(bin_index * bin_size), int((bin_index + 1) * bin_size)

                fout.write("{}\t{}\t{}\t{}".format(bin_index, bin_start, bin_end, bin_size))

                bin_total_hor_cnt = sum(sample_bin_list[bin_index].values())
                bin_freq_dict = {}

                for hor_type in sample_bin_list[bin_index].keys():
                    hor_cnt = sample_bin_list[bin_index][hor_type]

                    # # calcudate freq
                    hor_freq = hor_cnt / bin_total_hor_cnt if bin_total_hor_cnt != 0 else 0
                    bin_freq_dict[hor_type] = hor_freq

                    fout.write("\t{}:{},{}".format(hor_type, hor_cnt, round(hor_freq, 3)))

                fout.write("\n")

                sample_bin_freq_list.append(bin_freq_dict)

        # # plot
        plot_features = chrom_principle_hor_list + ["Other"]

        x = np.arange(len(sample_bin_freq_list))
        bottom = np.zeros(len(sample_bin_freq_list))

        fig, ax = plt.subplots(figsize=(12, 3))

        for feature in plot_features:
            values = np.array([d.get(feature, 0) for d in sample_bin_freq_list])
            ax.bar(x, values, bottom=bottom, label=feature, width=0.8)

            bottom += values

        ax.set_ylim(0, 1)
        ax.yaxis.set_major_formatter(PercentFormatter(1))

        ax.set_xticks(x)
        ax.set_xticklabels(np.arange(1, len(sample_bin_freq_list) + 1))
        ax.legend(loc="upper left", frameon=False)

        plt.tight_layout()
        plt.savefig(os.path.join(chrom_output_folder, "{}.{}.bin_{}.png".format(sample_cluster, sample_name, bin_count)), dpi=600)
        plt.close()



def plot_bin_size_level_hor_frequence_each_sample(bin_size=50):
    color_dict = {

        "HOR0": "#6688A6",  # 雾霾蓝
        "HOR1": "#D18B72",  # 陶土橙
        "HOR2": "#7FA18A",  # 鼠尾草绿
        "HOR3": "#A18BB5",  # 灰紫色
        "HOR4": "#C3A45D",  # 暗金色
        "HOR5": "#C47F91",  # 柔和玫瑰红
        "HOROther": "#858B91",  # 石板灰

    }

    output_folder = os.path.join(work_folder, "bin_size_{}_principal_hor_frequency".format(bin_size))
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for spec_chrom in chrom_order:

        # # load cluste
        cluster_dict = {}
        classify_file = os.path.join(pha_input_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom))
        with open(classify_file) as fin:
            for line in fin:
                line_split = line.strip().split()
                sample, cluster = line_split[0], line_split[-5]

                cluster_dict[sample] = cluster

        chrom_output_folder = os.path.join(output_folder, spec_chrom)
        if not os.path.exists(chrom_output_folder):
            os.mkdir(chrom_output_folder)

        chrom_principle_hor_list = c1_principle_hor_dict[spec_chrom]
        chrom_annot_res_folder = os.path.join(annot_res_folder, spec_chrom)

        chrom_color_dict = {"Other": color_dict["HOROther"]}
        for i in range(len(chrom_principle_hor_list)):
            chrom_color_dict[chrom_principle_hor_list[i]] = color_dict["HOR{}".format(i)]


        for sample_hor_file in os.listdir(chrom_annot_res_folder):
            if not sample_hor_file.endswith("#stv_hor_osat_te.pha.bed"):
                continue

            sample_name = sample_hor_file.replace("#stv_hor_osat_te.pha.bed", "")


            if sample_name == "T2T-CHM13.fa.gz":
                sample_cluster = "ref"

            else:
                if sample_name not in all_included_samples:
                    continue
                if "{}_RagTag".format(spec_chrom) not in whole_core_cent[sample_name]:
                    continue
                if sample_name not in cluster_dict:
                    continue

                sample_cluster = cluster_dict[sample_name]

            sample_hor_list = load_hor_file(os.path.join(chrom_annot_res_folder, sample_hor_file))
            if len(sample_hor_list) == 0:
                print("Empty HOR for sample", sample_name)
                continue


            # # STEP: begin to stats in each bin
            # # init
            bin_count = int(len(sample_hor_list) / bin_size) + 1
            sample_bin_list = []
            for bin_index in range(bin_count):
                sample_bin_list.append({hor: 0 for hor in chrom_principle_hor_list + ["Other"]})

            # # for each hor, check and add into the bin
            print(spec_chrom, sample_name, bin_size)

            hor_index = 0
            for _, hor_start, hor_end, hor_type, hor_orient in sample_hor_list:
                if hor_type in TE_types or hor_type in oSat_types:
                    continue

                if hor_orient == "-":
                    hor_type = reverse_hor(hor_type)
                hor_type = re.sub(r'\{\d+\}', '{n}', hor_type)

                hor_bin_index = int(hor_index / bin_size)

                if hor_type in chrom_principle_hor_list:
                    sample_bin_list[hor_bin_index][hor_type] += 1
                else:
                    sample_bin_list[hor_bin_index]["Other"] += 1

                hor_index += 1

            # # output into file
            sample_bin_freq_list = []
            sample_bin_output_file = os.path.join(chrom_output_folder, "{}.{}.bin_{}.txt".format(sample_cluster, sample_name, bin_size))
            with open(sample_bin_output_file, "w") as fout:
                for bin_index in range(len(sample_bin_list)):
                    bin_start, bin_end = int(bin_index * bin_size), int((bin_index + 1) * bin_size)

                    fout.write("{}\t{}\t{}\t{}".format(bin_index, bin_start, bin_end, bin_size))

                    bin_total_hor_cnt = sum(sample_bin_list[bin_index].values())
                    bin_freq_dict = {}

                    for hor_type in sample_bin_list[bin_index].keys():
                        hor_cnt = sample_bin_list[bin_index][hor_type]

                        # # calcudate freq
                        hor_freq = hor_cnt / bin_total_hor_cnt if bin_total_hor_cnt != 0 else 0
                        bin_freq_dict[hor_type] = hor_freq

                        fout.write("\t{}:{},{}".format(hor_type, hor_cnt, round(hor_freq, 3)))

                    fout.write("\n")

                    sample_bin_freq_list.append(bin_freq_dict)

            # # plot
            plot_features = chrom_principle_hor_list + ["Other"]

            x = np.arange(len(sample_bin_freq_list))
            bottom = np.zeros(len(sample_bin_freq_list))

            fig, ax = plt.subplots(figsize=(12, 3))

            for feature in plot_features:
                values = np.array([d.get(feature, 0) for d in sample_bin_freq_list])
                ax.bar(x, values, bottom=bottom, label=feature, width=0.8, color=chrom_color_dict[feature],     rasterized=True)

                bottom += values

            ax.set_ylim(0, 1)
            ax.yaxis.set_major_formatter(PercentFormatter(1))

            ax.set_xticks(x)
            ax.set_xticklabels(np.arange(1, len(sample_bin_freq_list) + 1))
            ax.legend(loc="upper left", frameon=False)

            plt.tight_layout()
            plt.savefig(os.path.join(chrom_output_folder, "{}.{}.bin_{}.png".format(sample_cluster, sample_name, bin_size)))
            plt.savefig(os.path.join(chrom_output_folder, "{}.{}.bin_{}.svg".format(sample_cluster, sample_name, bin_size)))

            plt.close()

            if spec_chrom == "chr9":
                fig, ax = plt.subplots(figsize=(12, 3))
                feature = "S2C9H1L.1-3_7"
                x = np.arange(len(sample_bin_freq_list))
                values = np.array([d.get(feature, 0) for d in sample_bin_freq_list])
                ax.plot( x, values, marker="o", markersize=5, linewidth=1.5, color=chrom_color_dict[feature], rasterized=True)

                ax.set_ylim(0, 1)
                ax.yaxis.set_major_formatter(PercentFormatter(1))

                ax.set_xticks(x)
                ax.set_xticklabels(np.arange(1, len(sample_bin_freq_list) + 1))

                plt.tight_layout()
                plt.savefig(os.path.join(chrom_output_folder, "{}.{}.bin_{}.var3.svg".format(sample_cluster, sample_name, bin_size)))
                plt.close()



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


def plot_bin_level_hor_frequence_among_popu(bin_count=20, popu_mode="hyper"):

    if popu_mode == "hyper":
        popu_dict = hyper_popu_dict
    else:
        popu_dict = detail_popu_dict


    output_folder = os.path.join(work_folder, "bin_{}_principal_hor_frequency_{}_popu".format(bin_count, popu_mode))
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for spec_chrom in chrom_order:

        chrom_input_folder = os.path.join(work_folder, "bin_{}_principal_hor_frequency/{}".format(bin_count, spec_chrom))
        if not os.path.exists(chrom_input_folder):
            continue

        chrom_principle_hor_list = c1_principle_hor_dict[spec_chrom]

        chrom_features = chrom_principle_hor_list + ["Other"]
        bin_columns = ["bin", "bin_start", "bin_end", "bin_length"]

        # # collect from all samples
        all_df = []
        for sample, population in popu_dict.items():
            if "{}_RagTag".format(spec_chrom) not in whole_core_cent[sample]:
                continue

            sample_bin_file = os.path.join(chrom_input_folder, "{}.bin_{}.txt".format(sample, bin_count))
            if not os.path.exists(sample_bin_file):
                continue

            df = pd.read_csv(sample_bin_file, sep="\t", header=None)
            df.columns = bin_columns + chrom_features

            records = []
            for _, row in df.iterrows():
                for feature in chrom_features:
                    value = str(row[feature]).strip().split(":")[1]

                    count, freq = value.split(",")

                    records.append({ "Sample": sample, "Population": population, "Bin": row["bin"], "Feature": feature, "Count": int(count), "Frequency": float(freq) })

            all_df.append(pd.DataFrame(records))

        all_df = pd.concat(all_df, ignore_index=True)
        # print(all_df)

        # # STEP: plot a boxplot
        for population, pop_df in all_df.groupby("Population"):
            plt.figure(figsize=(16, 7))

            sns.boxplot( data=pop_df, x="Bin", y="Frequency", hue="Feature", width=0.8, flierprops={"marker": "o", "markersize": 0.5, "linestyle": "none"}, )

            plt.legend(loc="upper left", frameon=False)

            plt.tight_layout()
            plt.savefig(os.path.join(output_folder, "{}.{}.boxplot.png".format(spec_chrom, population)), dpi=600)
            plt.close()


        # # STEP: calculate a p value for each feature in each bin
        stats_results = []
        for (bin_id, feature), group in all_df.groupby(["Bin", "Feature"]):
            populations = sorted(group["Population"].dropna().unique())

            for pop1, pop2 in combinations(populations, 2):
                values1 = group.loc[group["Population"] == pop1, "Frequency"].dropna().values
                values2 = group.loc[group["Population"] == pop2, "Frequency"].dropna().values
                if len(values1) == 0 or len(values2) == 0:
                    continue

                stat, pvalue = mannwhitneyu(values1, values2, alternative="two-sided" )

                # # # #
                n1 = len(values1)
                n2 = len(values2)

                combined = np.concatenate([values1, values2])
                ranks = pd.Series(combined).rank(method="average").values

                R1 = ranks[:n1].sum()

                # U for x
                U1 = R1 - n1 * (n1 + 1) / 2

                # P(X > Y) - P(X < Y)
                cliff_delta = (
                        2 * U1 / (n1 * n2) - 1
                )

                print(spec_chrom, bin_id, feature, len(values1), len(values2), pvalue, cliff_delta, sep="\t")
                stats_results.append({"Bin": bin_id, "Feature": feature, "Population1": pop1, "Population2": pop2,
                    "N1": len(values1), "N2": len(values2), "U": stat, "P_value": pvalue})

        stats_results_df = pd.DataFrame(stats_results)

        # FDR correction
        stats_results_df["FDR"] = (stats_results_df.groupby(["Bin", "Feature"])["P_value"].transform(lambda x: multipletests(x, method="fdr_bh")[1]))
        stats_results_df["Significance"] = stats_results_df["FDR"].apply(significance_label)

        stats_results_df.to_csv(os.path.join(output_folder, "{}.stats.txt".format(spec_chrom)), sep="\t", float_format="%.2e", index=False)


        # # STEP: plot a cluster heatmap
        sample_matrix = (all_df.pivot_table(index="Sample", columns=["Bin", "Feature"], values="Frequency", aggfunc="mean"))

        populations = sorted(set(popu_dict.values()))
        palette = dict(zip(populations, sns.color_palette("tab10", len(populations))))
        row_colors = sample_matrix.index.map( lambda x: palette[popu_dict[x]])

        distance_df = pd.DataFrame(squareform(pdist(sample_matrix.values, metric="euclidean")), index=sample_matrix.index, columns=sample_matrix.index)

        g = sns.clustermap(distance_df, cmap="Reds", figsize=(12, 12), row_colors=row_colors, col_colors=row_colors, row_cluster=True, col_cluster=True, linewidths=0 )

        # # add legend for row colors
        handles = [Patch(facecolor=palette[pop], label=pop) for pop in populations]
        g.ax_row_dendrogram.legend(handles=handles, title="Population", loc="best", frameon=False )

        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, "{}.heatmap.png".format(spec_chrom)), dpi=600)
        plt.close()


if __name__ == '__main__':

    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    # chrom_order = ["chr9"]

    work_folder = "/mnt/e/workspace/cpc2_centro/res_all_post"
    annot_res_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot"
    stats_res_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot/chrall"
    pha_input_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier"

    # # STEP: load whole centro info
    with open(os.path.join(work_folder, "stats.whole_core_cent.json")) as file:
        whole_core_cent = dict(json.load(file))

    # # STEP: load population info
    all_included_samples = []
    detail_popu_dict = {}
    detail_popu_order = ["AFR", "EUR", "AMR", "SAS", "EAS", "CHN"]

    hyper_popu_dict = {}
    hyper_popu_order = ["Other", "CHN"]

    with open("../data/data_population_hgsvc3_hprc2_cpc2.txt") as fin:
        for line in fin:
            line_split = line.strip().split("\t")

            all_included_samples.append(line_split[1])
            detail_popu_dict[line_split[1]] = line_split[6]
            hyper_popu_dict[line_split[1]] = line_split[5]

    TE_types = ["LINE", "SINE"]
    oSat_types = ["HSat1A", "HSat1B", "HSat2", "HSat3", "BSat", "GSat"]

    # # this is calculated by the find_principle_hor(), by removing single-monomer HOR, but S3C11H1L.1 is saved, because it actually is S3C11H1L.1_3-5
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


    # # STEP: find the principle hor for each chrom
    find_principle_hor()
    # get_cnt_and_freq_for_one_HOR(spec_chrom="chr7", spec_HOR="S1C7H1L.1-5_2-4_3-6")

    # # STEP: plot the principle HOR freq among popu
    plot_hor_frequence_among_popu(popu_mode="detail")

    # # STEP: plot the bin level principle HOR frequency for each sample
    process_pool = multiprocessing.Pool(10)
    for spec_chrom in chrom_order:
        process_pool.apply_async(plot_bin_level_hor_frequence_each_sample, (spec_chrom, 10))
    process_pool.close()
    process_pool.join()

    # # STEP: plot the bin level principle HOR frequency among popu
    plot_bin_level_hor_frequence_among_popu(bin_count=10)



    # # STEP: plot bin size levle
    # plot_bin_size_level_hor_frequence_each_sample(100)

