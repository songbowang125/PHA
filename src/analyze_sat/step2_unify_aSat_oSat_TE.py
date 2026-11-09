import multiprocessing
import os
import json
import re
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

def reverse_hor(in_hor):
    name = in_hor.split(".")[0]
    combine = in_hor.split(".")[1]

    # # reverse combine
    ids = re.split(r'[-_]', combine)
    seps = re.findall(r'[-_]', combine)

    ids.reverse()
    seps.reverse()

    return name + "." + ids[0] + ''.join(sep + id_ for sep, id_ in zip(seps, ids[1:]))



def overlap_oSat_and_TE_with_aSat(spec_chrom):

    with open(os.path.join(work_folder, "stats.whole_core_cent.json")) as file:
        whole_core_cent = dict(json.load(file))

    chrom_annot_folder = os.path.join(annot_folder, spec_chrom)


    for sample_aSat_file in os.listdir(chrom_annot_folder):

        if not sample_aSat_file.endswith("#stv_hor.bed"):
            continue

        sample = sample_aSat_file.replace("#stv_hor.bed", "")


        if sample != "T2T-CHM13.fa.gz":
            if sample not in whole_core_cent:
                continue
            if "{}_RagTag".format(spec_chrom) not in whole_core_cent[sample]:
                continue

        print(spec_chrom, sample)
        sample_aSat_file = os.path.join(chrom_annot_folder, "{}#stv_hor.bed".format(sample))
        sample_oSat_file = os.path.join(chrom_annot_folder, "{}#oSAT.bed".format(sample))
        sample_TE_file = os.path.join(chrom_annot_folder, "{}#rpmask.bed".format(sample))


        if within_PHA_flag:
            if limit_types:
                sample_output_tmp_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.pha.tmpbed")
                sample_output_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.pha.bed")
            else:
                sample_output_tmp_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.pha.full.tmpbed")
                sample_output_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.pha.full.bed")
        else:
            if limit_types:
                sample_output_tmp_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.tmpbed")
                sample_output_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.bed")
            else:
                sample_output_tmp_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.full.tmpbed")
                sample_output_file = sample_aSat_file.replace("#stv_hor.bed", "#stv_hor_osat_te.full.bed")

        with open(sample_output_tmp_file, "w") as fout:

            # # deal with aSat and find the relative
            aSat_contig, aSat_contig_start, aSat_start, aSat_end = None, None, np.inf, -np.inf
            with open(sample_aSat_file) as fin:
                for line in fin:

                    line_split = line.strip().split()

                    aSat_type, aSat_orient = line_split[3], line_split[4]
                    if aSat_orient == "-":
                        aSat_type = reverse_hor(aSat_type)
                    aSat_type = re.sub(r'\{\d+\}', '{n}', aSat_type)

                    aSat_contig, aSat_contig_start, aSat_contig_end, relative_start, relative_end = line_split[0], int(line_split[0].split(":")[1].split("-")[0]), int(line_split[0].split(":")[1].split("-")[1]), int(line_split[1]), int(line_split[2])

                    if not within_PHA_flag:
                        aSat_start = aSat_contig_start
                        aSat_end = aSat_contig_end

                        fout.write(line)

                    else:
                        if aSat_type not in c1_principle_hor_dict[spec_chrom]:
                            continue

                        if aSat_type in c1_principle_hor_single_list:
                            continue

                        fout.write(line)

                        relative_start, relative_end = aSat_contig_start + relative_start, aSat_contig_start + relative_end

                        if relative_start < aSat_start:
                            aSat_start = relative_start
                        if relative_end > aSat_end:
                            aSat_end = relative_end


            # # deal with oSat
            with open(sample_oSat_file) as fin:
                for line in fin:
                    line_split = line.strip().split()

                    contig, relative_start, relative_end, oSat_type = line_split[0], int(line_split[1]), int(line_split[2]), line_split[3]

                    oSat_chrom = contig.split(":")[0].replace("_RagTag", "")
                    if oSat_chrom != spec_chrom:
                        continue

                    relative_start = int(line_split[0].split(":")[1].split("-")[0]) + relative_start
                    relative_end = int(line_split[0].split(":")[1].split("-")[0]) + relative_end

                    # # overlap
                    if relative_start < aSat_end and aSat_start < relative_end:
                        fout.write("{}\t{}\t{}\t{}\tNA\n".format(aSat_contig, relative_start - aSat_contig_start, relative_end - aSat_contig_start, oSat_type))


            # # deal with TE
            with open(sample_TE_file) as fin:
                for line in fin:
                    line_split = line.strip().split()

                    contig, relative_start, relative_end, TE_type = line_split[0], int(line_split[1]), int(line_split[2]), line_split[3].split("/")[0]

                    if limit_types:
                        if TE_type not in TE_types:
                            continue
                    else:
                        if TE_type in ["Satellite"]:
                            continue

                    oSat_chrom = contig.split(":")[0].replace("_RagTag", "")
                    if oSat_chrom != spec_chrom:
                        continue

                    relative_start = int(line_split[0].split(":")[1].split("-")[0]) + relative_start
                    relative_end = int(line_split[0].split(":")[1].split("-")[0]) + relative_end

                    # # overlap
                    if relative_start < aSat_end and aSat_start < relative_end:
                        fout.write("{}\t{}\t{}\t{}\tNA\n".format(aSat_contig, max(relative_start - aSat_contig_start, 0), relative_end - aSat_contig_start, TE_type))

        os.system("bedtools sort -i {} > {}".format(sample_output_tmp_file, sample_output_file))
        os.remove(sample_output_tmp_file)




def plot_hor_osat_te_distribution_each_sample(spec_chrom, ):

    output_folder = os.path.join(work_folder, "distribution_hor_osat_te_new")
    if not os.path.exists(output_folder):
        os.system("mkdir {}".format(output_folder))


    # # # init color dict
    # color_dict = {
    #
    #     # "HOR0": "#501549", "HOR1": "#7f340d", "HOR2": "#78206e", "HOR3": "#bf4f14", "HOR4": "#d76dcc", "HOR5": "#f1a984", "HOR6": "#e49edd", "HOR7": "#f6c6ac", "HOR8": "#f1ceee", "HOR9": "#fae2d6",
    #     "HOR0": "#156082", "HOR1": "#e97132", "HOR2": "#196b24", "HOR3": "#0f9ed5", "HOR4": "#a02b93", "HOR5": "#4ea72e", "HOR6": "#000000", "HOR7": "#f6c6ac", "HOR8": "#f1ceee", "HOR9": "#fae2d6",
    #
    #     "HSat1A": "#18988b", "HSat1B": "#18988b",
    #     "HSat2": "#323366", "HSat3": "#77a0ba",
    #     "BSat": "#0f9ed5", "GSat": "#758660",
    #
    #     "SINE": "#ff0000", "LINE": "#0080ff",
    #
    #     "+": "#ffffff",
    #     "NA": "#ffffff",
    #     "-": "#000000"
    #
    # }
    color_dict = {

        "HOR0": "#6688A6",  # 雾霾蓝
        "HOR1": "#D18B72",  # 陶土橙
        "HOR2": "#7FA18A",  # 鼠尾草绿
        "HOR3": "#A18BB5",  # 灰紫色
        "HOR4": "#C3A45D",  # 暗金色
        "HOR5": "#C47F91",  # 柔和玫瑰红
        "HOROther": "#858B91",  # 石板灰

        "HSat1A": "#18988b", "HSat1B": "#18988b",
        "HSat2": "#323366", "HSat3": "#77a0ba",
        "BSat": "#dcaed0", "GSat": "#758660",

        "SINE": "#ff0000", "LINE": "#0080ff", "TE_Other": "#959595",

        "+": "#ffffff",
        "NA": "#ffffff",
        "-": "#000000"

    }

    # # load cluster
    cluster_dict = {}
    classify_file = os.path.join(pha_input_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom))
    with open(classify_file) as fin:
        for line in fin:
            line_split = line.strip().split()
            sample, cluster = line_split[0], line_split[-5]

            cluster_dict[sample] = cluster
    # #
    chrom_output_folder = os.path.join(output_folder, spec_chrom)
    if not os.path.exists(chrom_output_folder):
        os.mkdir(chrom_output_folder)

    chrom_principle_hor_list = c1_principle_hor_dict[spec_chrom]
    chrom_annot_res_folder = os.path.join(annot_folder, spec_chrom)

    # # for each sample
    for sample_annot_file in os.listdir(chrom_annot_res_folder):
        if not sample_annot_file.endswith("#stv_hor_osat_te.pha.full.bed"):
            continue

        sample_name = sample_annot_file.replace("#stv_hor_osat_te.pha.full.bed", "")
        if sample_name == "T2T-CHM13.fa.gz":

            sample_cluster = "ref"
        else:
            if sample_name not in cluster_dict:
                continue
            sample_cluster = cluster_dict[sample_name]


        print(spec_chrom, sample_name)
        # # STEP: load info
        sample_info_list = []
        with open(os.path.join(chrom_annot_res_folder, sample_annot_file)) as fin:
            for line in fin:
                line_split = line.strip().split("\t")
                chrom, start, end, type, orient = line_split[0], int(line_split[1]), int(line_split[2]), line_split[3], line_split[4]

                orient_color = color_dict[orient]

                if "." in type:
                    # # belong to hor
                    if orient == "-":
                        type = reverse_hor(type)
                    type = re.sub(r'\{\d+\}', '{n}', type)

                    if type in chrom_principle_hor_list:
                        type = "HOR{}".format(chrom_principle_hor_list.index(type))
                    else:
                        type = "HOROther"

                    type_color = color_dict[type]
                elif type in oSat_types:
                    type_color = color_dict[type]

                else:
                    if type in color_dict:
                        type_color = color_dict[type]
                    else:
                        type_color = color_dict["TE_Other"]


                sample_info_list.append([chrom, start, end, type, type_color, orient_color])

        if len(sample_info_list) == 0:
            continue

        # # STEP: plot info
        starts, ends = [x[1] for x in sample_info_list], [x[2] for x in sample_info_list]

        xmin, xmax = min(starts), max(ends)

        genome_length = xmax - xmin
        margin = genome_length * 0.02

        plot_min, plot_max = xmin - margin, xmax + margin

        fig, axes = plt.subplots( 3, 1, figsize=(14, 2), sharex=True, gridspec_kw={"height_ratios": [1, 1, 1]})
        ax_hor, ax_orient, ax_te = axes[0], axes[1], axes[2]

        # # plot sattlelites

        # # track 1
        for chrom, start, end, feature, type_color, orient_color in sample_info_list:
            if feature in TE_types + oSat_types:
                continue

            rect = Rectangle((start, 0), end - start, 1, facecolor=type_color, rasterized=True )
            ax_hor.add_patch(rect)

        ax_hor.set_ylim(0, 1)
        ax_hor.set_yticks([])

        # # track 2
        for chrom, start, end, feature, type_color, orient_color in sample_info_list:
            if feature not in oSat_types:
                continue

            rect = Rectangle((start, 0), end - start, 1, facecolor=type_color,  rasterized=True  )
            ax_orient.add_patch(rect)

        ax_orient.set_ylim(0, 1)
        ax_orient.set_yticks([])

        # for chrom, start, end, feature, type_color, orient_color in sample_info_list:
        #
        #     # 中间画一条线
        #     y = 0.5
        #
        #     ax_orient.plot( [start, end], [y, y], color=orient_color, linewidth=2 )
        #
        # ax_orient.set_ylim(0, 1)
        # ax_orient.set_yticks([])

        # # track 3
        for chrom, start, end, feature, type_color, orient_color in sample_info_list:

            if feature not in TE_types:
                continue

            x = (start + end) / 2

            ax_te.scatter(x, 0.5, marker="^", s=45, facecolor=type_color, edgecolor="black", linewidth=0.4, zorder=5)

        ax_te.set_ylim(0, 1)
        ax_te.set_yticks([])

        # # set xaxis
        ax_te.set_xlim(plot_min, plot_max)
        ax_te.set_xticks([])

        # # remove boarders
        for ax in axes:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)

        # # add legend
        # legend_handles = [ Patch( facecolor=color, edgecolor="none", label=feature ) for feature, color in color_dict.items() ]
        # ax_hor.legend( handles=legend_handles, title="HOR", loc="upper center", bbox_to_anchor=(0.5, 1.8), ncol=4, frameon=False)

        plt.tight_layout()
        plt.savefig(os.path.join(chrom_output_folder, "{}_{}.distribution.png".format(sample_cluster, sample_name)))
        plt.savefig(os.path.join(chrom_output_folder, "{}_{}.distribution.svg".format(sample_cluster, sample_name)))

        plt.close()



def plot_type_length_among_cluster():


    plot_output_folder = os.path.join(work_folder, "principal_hor_type_length_cluster_2")
    if not os.path.exists(plot_output_folder):
        os.mkdir(plot_output_folder)


    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

    type_order = ["HOR"] + TE_types + oSat_types


    for spec_type in type_order:
        type_output_fout = open(os.path.join(plot_output_folder, "res.{}_length.txt".format(spec_type)), "w")

        for spec_chrom in chrom_order:
            chrom_cluster_file = os.path.join(pha_input_folder, "{}.features.c2.txt".format(spec_chrom))
            chrom_input_folder = os.path.join(annot_folder, spec_chrom)

            print(spec_type, spec_chrom)
            # # load cluster info
            collect_dict = {}
            sample_cluster_dict = {}
            with open(chrom_cluster_file) as fin:
                for line in fin:
                    if "Sample" in line:
                        continue

                    line_split = line.strip().split()
                    sample, cluster = line_split[0], line_split[-1]

                    sample_cluster_dict[sample] = cluster

                    if cluster not in collect_dict:
                        collect_dict[cluster] = []
            cluster_order = sorted(list(collect_dict.keys()))

            # # traverse each file
            for sample_bed in os.listdir(chrom_input_folder):
                if not sample_bed.endswith("#stv_hor_osat_te.pha.full.bed"):
                    continue

                sample = sample_bed.replace("#stv_hor_osat_te.pha.full.bed", "")
                if sample not in sample_cluster_dict:
                    continue

                sample_feature_length = 0
                with open(os.path.join(chrom_input_folder, sample_bed)) as fin:
                    for line in fin:
                        line_split = line.strip().split()

                        start, end, type = int(line_split[1]), int(line_split[2]), line_split[3]

                        if type not in TE_types + oSat_types:
                            type = "HOR"

                        if type == spec_type:
                            sample_feature_length += (end - start + 1)

                collect_dict[sample_cluster_dict[sample]].append(sample_feature_length)

            # # begin to output
            for cluster in cluster_order:
                length_avg, length_std = int(np.average(collect_dict[cluster])), int(np.std(collect_dict[cluster]))

                length_min, length_max = np.min(collect_dict[cluster]), np.max(collect_dict[cluster])

                type_output_fout.write("{}\t{}\t{}\tAvg:{}\tMin:{}\tMax:{}\tStd:{}\t{}\n".format(spec_chrom, cluster, len(collect_dict[cluster]), length_avg, length_min, length_max, length_std, "\t".join([str(cnt) for cnt in collect_dict[cluster]])))


            # # begin to plot
            all_values, all_groups = [], []

            for cluster in cluster_order:
                all_values.extend(collect_dict[cluster])
                all_groups.extend([cluster] * len(collect_dict[cluster]))

            df = pd.DataFrame({"value": all_values, "group": all_groups})

            fig, ax = plt.subplots(figsize=(3, 5))

            # sns.violinplot(data=df, x="group", y="value", inner=None, cut=0, linewidth=1, color="#5975a5", ax=ax)
            sns.boxplot(data=df, x="group", y="value", color="#5975a5", showfliers=False, ax=ax)
            # sns.swarmplot(data=df, x="group", y="value", size=1, color="#fe8d24", ax=ax)

            sns.despine(ax=ax)
            ax.set_xlabel("")
            ax.set_ylabel("")
            # ax.set_ylim(0, 10000000)
            ax.yaxis.get_offset_text().set_visible(False)

            plt.tight_layout()
            plt.savefig(os.path.join(plot_output_folder, "res.{}_length.{}.png".format(spec_type, spec_chrom)))
            plt.close()

        type_output_fout.close()

def plot_type_length_among_popu():


    plot_output_folder = os.path.join(work_folder, "principal_hor_type_length_hyper_popu")
    if not os.path.exists(plot_output_folder):
        os.mkdir(plot_output_folder)


    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

    type_order = ["HOR"] + TE_types + oSat_types


    for spec_type in type_order:
        type_output_fout = open(os.path.join(plot_output_folder, "res.{}_length.txt".format(spec_type)), "w")

        for spec_chrom in chrom_order:
            chrom_input_folder = os.path.join(annot_folder, spec_chrom)

            print(spec_type, spec_chrom)
            # # load popu info
            collect_dict = {}
            sample_popu_dict = {}
            with open("../data/data_population_hgsvc3_hprc2_cpc2.txt") as fin:
                for line in fin:
                    if "Sample" in line:
                        continue

                    line_split = line.strip().split("\t")
                    sample, popu = line_split[1], line_split[5]

                    sample_popu_dict[sample] = popu

                    if popu not in collect_dict:
                        collect_dict[popu] = []
            popu_order = sorted(list(collect_dict.keys()))

            # # traverse each file
            for sample_bed in os.listdir(chrom_input_folder):
                if not sample_bed.endswith("#stv_hor_osat_te.pha.full.bed"):
                    continue

                sample = sample_bed.replace("#stv_hor_osat_te.pha.full.bed", "")
                if sample not in sample_popu_dict:
                    continue

                sample_feature_length = 0
                with open(os.path.join(chrom_input_folder, sample_bed)) as fin:
                    for line in fin:
                        line_split = line.strip().split()

                        start, end, type = int(line_split[1]), int(line_split[2]), line_split[3]

                        if type not in TE_types + oSat_types:
                            type = "HOR"

                        if type == spec_type:
                            sample_feature_length += (end - start + 1)

                collect_dict[sample_popu_dict[sample]].append(sample_feature_length)

            # # begin to output
            for popu in popu_order:

                length_avg, length_std = int(np.average(collect_dict[popu])), int(np.std(collect_dict[popu]))

                length_min, length_max = np.min(collect_dict[popu]), np.max(collect_dict[popu])

                type_output_fout.write("{}\t{}\t{}\tAvg:{}\tMin:{}\tMax:{}\tStd:{}\t{}\n".format(spec_chrom, popu, len(collect_dict[popu]), length_avg, length_min, length_max, length_std, "\t".join([str(cnt) for cnt in collect_dict[popu]])))


            # # begin to plot
            all_values, all_groups = [], []

            for popu in popu_order:
                all_values.extend(collect_dict[popu])
                all_groups.extend([popu] * len(collect_dict[popu]))

            df = pd.DataFrame({"value": all_values, "group": all_groups})

            fig, ax = plt.subplots(figsize=(3, 5))

            # sns.violinplot(data=df, x="group", y="value", inner=None, cut=0, linewidth=1, color="#5975a5", ax=ax)
            sns.boxplot(data=df, x="group", y="value", color="#5975a5", showfliers=False, ax=ax)
            # sns.swarmplot(data=df, x="group", y="value", size=1, color="#fe8d24", ax=ax)

            sns.despine(ax=ax)
            ax.set_xlabel("")
            ax.set_ylabel("")
            # ax.set_ylim(0, 10000000)
            ax.yaxis.get_offset_text().set_visible(False)

            plt.tight_layout()
            plt.savefig(os.path.join(plot_output_folder, "res.{}_length.{}.png".format(spec_type, spec_chrom)))
            plt.close()

        type_output_fout.close()
if __name__ == '__main__':

    work_folder = "/mnt/e/workspace/cpc2_centro/res_all_post"
    annot_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot"
    pha_input_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier"

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

    # # STEP: load cent info


    within_PHA_flag = True
    limit_types = False

    chrom_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
    # chrom_order = ["chr4"]

    TE_types = ["LINE", "SINE"]
    oSat_types = ["HSat1A", "HSat1B", "HSat2", "HSat3", "BSat", "GSat"]

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
                               "chrX": ["S3CXH1L.1-12"]
                               }

    c1_principle_hor_single_list = ["S2C9H1L.7", "S3C11H1L.1", "S2C15H1L.1", "S2C15H1L.7", "S1C1/5/19H1L.1", "S2C14/22H1L.7"]


    # # STEP: find osat and te around HOR
    process_pool = multiprocessing.Pool(24)
    for spec_chrom in chrom_order:
        process_pool.apply_async(overlap_oSat_and_TE_with_aSat, (spec_chrom, ))
    process_pool.close()
    process_pool.join()

    # # STEP: plot the distribution
    process_pool = multiprocessing.Pool(24)
    for spec_chrom in chrom_order:
        process_pool.apply_async(plot_hor_osat_te_distribution_each_sample, (spec_chrom,))
        # plot_hor_osat_te_distribution_each_sample(spec_chrom,)
    process_pool.close()
    process_pool.join()


    # # STEP:
    # plot_type_length_among_cluster()
    # plot_type_length_among_popu()