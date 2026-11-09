import multiprocessing
import pandas as pd
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO
import numpy as np
import os
import pandas as pd
from ete3 import Tree, TreeStyle, faces, TextFace
from collections import defaultdict
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from scipy import stats
from scipy.stats import mannwhitneyu
from collections import Counter
from ete3 import Tree, TreeStyle, NodeStyle, faces, CircleFace
import re
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from scipy.stats import kruskal


def generate_mafft_distance_matrxi(spec_chrom):
    if "flank" in work_folder:
        mafft_file = os.path.join(work_folder, "{}.flank.mafft.fa".format(spec_chrom))
    else:
        mafft_file = os.path.join(work_folder, "{}.hor.mafft.fa".format(spec_chrom))

    missing_chars = {"-"}

    alignment = AlignIO.read(mafft_file, "fasta")

    names = [record.id for record in alignment]
    sequences = [str(record.seq).upper() for record in alignment]
    n = len(names)

    order = sorted(range(n), key=lambda i: names[i])
    names = [names[i] for i in order]
    sequences = [sequences[i] for i in order]

    # =========================
    # 3. Calculate p-distance
    # =========================
    distance_matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        seq1 = np.array(list(sequences[i]))

        for j in range(i + 1, n):
            seq2 = np.array(list(sequences[j]))

            # 排除缺失字符，仅比较双方都有有效字符的位置
            valid = (
                    ~np.isin(seq1, list(missing_chars)) &
                    ~np.isin(seq2, list(missing_chars))
            )

            n_valid = valid.sum()

            if n_valid == 0:
                distance = np.nan
            else:
                distance = np.mean(seq1[valid] != seq2[valid])

            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    # =========================
    # 4. Save matrix
    # =========================
    df = pd.DataFrame(
        distance_matrix,
        index=names,
        columns=names
    )

    df.to_csv(os.path.join(work_folder, "{}.mafft_dist.matrix.txt".format(spec_chrom)), float_format="%.8f")


def generate_tree_distance_matrix(spec_chrom):


    if "flank" in work_folder:
        tree = Tree(os.path.join(work_folder, "{}.flank.mafft.fa.treefile".format(spec_chrom)), format=1)
    else:
        tree = Tree(os.path.join(work_folder, "{}.hor.mafft.fa.treefile".format(spec_chrom)), format=1)
    # midpoint = tree.get_midpoint_outgroup()
    # tree.set_outgroup(midpoint)

    names = sorted(tree.get_leaf_names())
    if "8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz" in names:
        names.remove("8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz")

    clusters = [name.split("_")[0] for name in names]
    n = len(names)

    distance_matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            distance = tree.get_distance(names[i], names[j])

            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    # np.save(os.path.join(work_folder, "{}.tree_dist.npy".format(spec_chrom)), distance_matrix)

    # # plot a heatmap
    print(spec_chrom, "plot heatmap")

    distance_df = pd.DataFrame(distance_matrix, index=names, columns=names)
    distance_df.to_csv(os.path.join(work_folder, "{}.tree_dist.matrix.txt".format(spec_chrom)))

    plt.figure(figsize=(12, 10))
    sns.heatmap(distance_df, cmap="viridis", square=True, xticklabels=True, yticklabels=True, linewidths=0, cbar_kws={"label": "Pairwise p-distance"}, annot=False, fmt=".2f")
    plt.xlabel("")
    plt.ylabel("")
    plt.xticks(rotation=90, fontsize=6)
    plt.yticks(rotation=0, fontsize=6)
    plt.tight_layout()
    plt.savefig( os.path.join(work_folder, "{}.tree_dist.png".format(spec_chrom)), dpi=600)
    plt.close()


    # # plot the boxplot
    print(spec_chrom, "plot boxplot")

    output_file = open(os.path.join(work_folder, "{}.tree_dist.txt".format(spec_chrom)), "w")

    distance_matrix = np.asarray(distance_matrix)
    clusters = np.asarray(clusters)

    unique_clusters = np.unique(clusters)

    for from_cluster in unique_clusters:

        from_idx = np.where(clusters == from_cluster)[0]

        plot_data = []


        to_cluster_values_dict = {}

        for to_cluster in unique_clusters:

            to_idx = np.where(clusters == to_cluster)[0]

            values = distance_matrix[np.ix_(from_idx, to_idx)]

            if from_cluster == to_cluster:
                # Only upper triangle
                values = values[np.triu_indices(len(from_idx), k=1)]

            else:
                # All pairwise distances between two clusters
                values = values.flatten()

            values = values[~np.isnan(values)]
            if len(values) == 0:
                continue


            q1 = np.percentile(values, 25)
            median = np.percentile(values, 50)
            q3 = np.percentile(values, 75)
            mean = np.average(values)

            to_cluster_values_dict[to_cluster] = {"Values": values, "Mean": mean, "Median": median, "Q3": q3, "Q1": q1, "N_dist": len(values), "N_sample": len(to_idx)}


            plot_data.extend([{"Target": str(to_cluster), "Distance": x } for x in values])


        for to_cluster in to_cluster_values_dict:
            from_cluster_info = to_cluster_values_dict[from_cluster]
            to_cluster_info = to_cluster_values_dict[to_cluster]
            # # self cluster, directly output
            if to_cluster == from_cluster:
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(from_cluster, to_cluster, to_cluster_info["N_sample"], to_cluster_info["N_dist"], to_cluster_info["Q1"], to_cluster_info["Q3"], to_cluster_info["Median"], to_cluster_info["Mean"], ))

            else:

                average_fold = to_cluster_info["Mean"] / from_cluster_info["Mean"]
                stat, p = mannwhitneyu(from_cluster_info["Values"], to_cluster_info["Values"], alternative="two-sided")

                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(from_cluster, to_cluster, to_cluster_info["N_sample"], to_cluster_info["N_dist"], to_cluster_info["Q1"], to_cluster_info["Q3"], to_cluster_info["Median"], to_cluster_info["Mean"], average_fold, p))


        # # plot a boxplot
        plot_df = pd.DataFrame(plot_data)
        plt.figure(figsize=(7, 5))
        sns.boxplot(data=plot_df, x="Target", y="Distance", width=0.6, showfliers=False)
        plt.tight_layout()
        plt.savefig(os.path.join(work_folder, "{}.tree_dist.C{}.dist.png".format(spec_chrom, from_cluster)), dpi=600)



def assess_chr1_inversion_cluster():

    chr1_cluster_file = os.path.join(pha_folder, "{}.features.c2.txt".format("chr1"))
    chr1_annot_folder = os.path.join(annot_folder, "chr1")

    # # load cluster
    sample_cluster_dict = {}
    with open(chr1_cluster_file) as fin:
        for line in fin:
            if "Sample" in line:
                continue

            line_split = line.strip().split()
            sample, cluster = line_split[0], line_split[-1]
            if cluster == "-1":
                continue

            sample_cluster_dict[sample] = cluster

    # # load inv
    sample_inv_dict = {}
    sample_total_dict = {}
    for sample in sample_cluster_dict:
        sample_inv_dict[sample] = 0
        sample_total_dict[sample] = 0

        with open(os.path.join(chr1_annot_folder, "{}#monomer.bed".format(sample))) as fin:

            for line in fin:
                monomer, orient = line.strip().split()[3], line.strip().split()[4]

                if "L." not in monomer:
                    continue

                sample_total_dict[sample] += 1
                if orient == "-":
                    sample_inv_dict[sample] += 1

    for sample in sample_cluster_dict:
        print("{}\t{}\t{}\t{}\n".format(sample, sample_cluster_dict[sample], sample_inv_dict[sample], sample_total_dict[sample]), end="")


def assign_branch_with_time_assess(tree_path):


    chimp_name = "8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"
    HC_TIME = 6.2

    # # load tree
    t = Tree(tree_path, format=1)

    chimp_nodes = t.search_nodes(name=chimp_name)
    if len(chimp_nodes) != 1:
        raise ValueError(f"Found {len(chimp_nodes)} chimp nodes, expected exactly 1.")

    chimp = chimp_nodes[0]
    t.set_outgroup(chimp)


    human_leaves = [ leaf for leaf in t.get_leaves() if leaf.name != chimp_name]

    # # calucate Boostsuport
    support_cnt = 0
    all_cnt = 0
    for node in t.traverse():
        if node.is_leaf() or node.is_root():
            continue

        all_cnt += 1


        if node.name == "":
            continue
        value = float(node.name)

        if value > 95:
            support_cnt += 1

    print(support_cnt, all_cnt, round(support_cnt / all_cnt, 3))

    # ============================================================
    # 4. Calculate root-to-tip distances
    # ============================================================
    records = []

    for leaf in human_leaves:
        root_to_tip = t.get_distance(t, leaf, topology_only=False )
        records.append({"Sample": leaf.name, "root_to_tip_distance": root_to_tip})

    df = pd.DataFrame(records)

    # ============================================================
    # 5. Calculate clock statistics
    # ============================================================
    distances = df["root_to_tip_distance"].values

    mean_distance = np.mean(distances)
    median_distance = np.median(distances)
    std_distance = np.std(distances, ddof=1)

    min_distance = np.min(distances)
    max_distance = np.max(distances)

    range_distance = max_distance - min_distance

    cv = std_distance / mean_distance

    # relative range
    relative_range = range_distance / mean_distance

    # coefficient of variation in %
    cv_percent = cv * 100

    relative_range_percent = relative_range * 100

    # ============================================================
    # 6. Root-to-tip regression
    #
    # If the tree is perfectly clock-like, ordering leaves should
    # not strongly affect distance. Therefore the regression itself
    # is not a biological clock test.
    #
    # We use it mainly as a diagnostic plot.
    # ============================================================
    df = df.sort_values(
        "root_to_tip_distance"
    ).reset_index(drop=True)

    df["sample_index"] = np.arange(1, len(df) + 1)

    slope, intercept, r_value, p_value, stderr = stats.linregress(
        df["sample_index"],
        df["root_to_tip_distance"]
    )

    r_squared = r_value ** 2

    # ============================================================
    # 7. Print statistics
    # ============================================================

    print("\n" + "=" * 70)
    print("Root-to-tip distance statistics")
    print("=" * 70)

    print(f"Number of human tips : {len(distances)}")
    print(f"Mean                  : {mean_distance:.8f}")
    print(f"Median                : {median_distance:.8f}")
    print(f"SD                    : {std_distance:.8f}")
    print(f"Min                   : {min_distance:.8f}")
    print(f"Max                   : {max_distance:.8f}")
    print(f"Range                 : {range_distance:.8f}")
    print(f"CV                    : {cv:.4f}")
    print(f"CV (%)                : {cv_percent:.2f}%")
    print(f"Relative range (%)    : {relative_range_percent:.2f}%")

    # ============================================================
    # 8. Clock-likeness assessment
    # ============================================================
    print("\n" + "=" * 70)
    print("Clock-likeness assessment")
    print("=" * 70)

    if cv < 0.05:
        clock_status = "Strongly clock-like"

    elif cv < 0.10:

        clock_status = "Approximately clock-like"

    elif cv < 0.20:
        clock_status = "Moderate rate variation"

    else:
        clock_status = "Strong rate variation / non-clock-like"

    print("Status:", clock_status)

    # # ============================================================
    # # 9. Recommendation
    # # ============================================================
    #
    # print("\n" + "=" * 70)
    # print("Dating recommendation")
    # print("=" * 70)
    #
    # if cv < 0.10:
    #
    #     print(
    #         "The root-to-tip distances are relatively homogeneous."
    #     )
    #
    #     print(
    #         "Simple proportional scaling may be reasonable "
    #         "as a rough estimate."
    #     )
    #
    #     print(
    #         "However, a formal clock test / dating analysis "
    #         "is still preferable for publication."
    #     )
    #
    # elif cv < 0.20:
    #
    #     print(
    #         "There is noticeable rate variation among lineages."
    #     )
    #
    #     print(
    #         "Simple scaling should be treated cautiously."
    #     )
    #
    #     print(
    #         "LSD2 or r8s is recommended."
    #     )
    #
    # else:
    #
    #     print(
    #         "The tree shows substantial root-to-tip variation."
    #     )
    #
    #     print(
    #         "Simple scaling using the chimpanzee branch "
    #         "is NOT recommended."
    #     )
    #
    #     print(
    #         "Use a relaxed-rate dating approach such as "
    #         "LSD2 or penalized likelihood with r8s."
    #     )


def get_cluster_label(node, sample_cluster):
    """统计子树中各个 cluster 的样本数量。"""
    counts = Counter(
        sample_cluster[leaf.name]
        for leaf in node.get_leaves()
        if leaf.name in sample_cluster
    )

    if not counts:
        return "Unknown"

    return ", ".join(
        f"{cluster}: {count}"
        for cluster, count in sorted(
            counts.items(),
            key=lambda x: (-x[1], x[0])
        )
    )

# def cut_tree_by_dist(work_folder, spec_chrom, threshold=0.01):
#
#     tree = Tree(os.path.join(work_folder, "{}.flank.mafft.fa.treefile.time".format(spec_chrom)), format=0)
#
#     # Save the original number of terminal samples
#     n_tips_before = len(tree)
#
#     # Postorder traversal: process descendants before ancestors
#     for node in list(tree.traverse("postorder")):
#         # Skip terminal nodes and root
#         if node.is_leaf() or node.is_root():
#             continue
#
#         # Branch length from this node to its parent
#         if node.dist < threshold:
#             node.delete(prevent_nondicotomic=False)
#
#     tree.write(outfile=os.path.join(work_folder, "{}.flank.mafft.fa.treefile.time.cut.test".format(spec_chrom)), format=1)

def cut_tree_by_dist(work_folder, spec_chrom, threshold=0.01):

    # # load tree and sample cluster dict
    tree = Tree(os.path.join(work_folder, "{}.flank.mafft.fa.treefile.time".format(spec_chrom)), format=1)

    sample_cluster = {}
    for sample in tree.get_leaf_names():
        sample_cluster[sample] = sample.split("_")[0]

        # 使用显式栈实现前序遍历
    stack = [tree]

    while stack:
        node = stack.pop()

        if node.is_root():
            # 根节点不作为待折叠分支
            stack.extend(reversed(node.children))
            continue

        if node.dist < threshold:
            parent = node.up

            # 汇总折叠前的所有后代样本
            label = get_cluster_label(node, sample_cluster)

            # 创建替代叶节点，分支长度设为 0
            replacement = parent.add_child(
                name=label,
                dist=0.0
            )

            # 删除原节点及其整个子树
            parent.remove_child(node)

            # 不再访问该子树的任何后代
            continue

        # 当前分支足够长，继续检查其子节点
        stack.extend(reversed(node.children))

    tree.write(outfile=os.path.join(work_folder, "{}.flank.mafft.fa.treefile.time.cut".format(spec_chrom)), format=1)



# def assess_time(spec_chrom):
#
#     fa


def collapse_tree_by_cluster(input_tree, max_other_ratio=0.1):
    print(input_tree)
    if input_tree.endswith(".nex"):
        tree = Phylo.read(input_tree, "nexus")
    else:
        tree = Phylo.read(input_tree, "newick")

    def get_cluster(name):
        """提取样本名称中第一个下划线之前的 cluster。"""
        return name.split("_", 1)[0]

    def collapse(clade):
        """
        递归缩减分支。

        如果非主要 cluster 的样本比例不超过 max_other_ratio，
        则将当前分支缩减为一个节点。
        """
        leaves = clade.get_terminals()

        if not leaves:
            return

        # 统计每个 cluster 的样本数
        cluster_counts = Counter(
            get_cluster(leaf.name) for leaf in leaves
        )

        total_count = len(leaves)

        # 找到样本数最多的 cluster
        major_cluster, major_count = cluster_counts.most_common(1)[0]

        # 计算其他 cluster 的样本比例
        other_ratio = (total_count - major_count) / total_count

        # 满足比例阈值，缩减当前分支
        if other_ratio <= max_other_ratio:
            clade.name = f"{major_cluster}_n{total_count}_o{(total_count - major_count)}"
            clade.clades = []
            return

        # 不满足阈值，继续向下遍历子分支
        for child in list(clade.clades):
            collapse(child)

    collapse(tree.root)

    # 输出文件
    if input_tree.endswith(".nex"):
        output_tree = input_tree + ".collapsed.nex"
        Phylo.write(tree, output_tree, "nexus")
    else:
        output_tree = input_tree + ".collapsed.treefile"
        Phylo.write(tree, output_tree, "newick")

    return output_tree


def collapse_tree_by_cluster_pure(input_tree):

    if input_tree.endswith(".nex"):
        tree = Phylo.read(input_tree, "nexus")
    else:
        tree = Phylo.read(input_tree, "newick")

    def get_cluster(name):
        """提取样本名称中第一个下划线之前的 cluster。"""
        return name.split("_", 1)[0]

    def collapse(clade):
        """
        递归缩减分支。
        返回该分支下的叶节点所属 cluster 集合。
        """
        # 当前分支下的所有叶节点
        leaves = clade.get_terminals()
        clusters = {get_cluster(leaf.name) for leaf in leaves}

        # 所有叶节点均来自同一个 cluster
        if len(clusters) == 1:
            cluster = next(iter(clusters))
            sample_count = len(leaves)

            # 将当前分支缩减为一个叶节点
            clade.name = f"{cluster}_n{sample_count}"
            clade.clades = []

            return clusters

        # 不属于同一个 cluster，继续遍历子分支
        for child in list(clade.clades):
            collapse(child)

        return clusters

    collapse(tree.root)

    # 输出 NEXUS 格式
    if input_tree.endswith(".nex"):

        Phylo.write(tree, input_tree + ".collapsed.nex", "nexus")

        return input_tree + ".collapsed.nex"

    else:
        Phylo.write(tree, input_tree + ".collapsed.treefile", "newick")

        return input_tree + ".collapsed.treefile"


def cut_tree_by_depth(tree_path, max_depth=3):

    # # load tree and sample cluster dict
    tree = Tree(tree_path, format=1)

    sample_cluster = {}
    for sample in tree.get_leaf_names():
        sample_cluster[sample] = sample.split("_")[0]


    stack = [(tree, 0)]

    while stack:
        node, depth = stack.pop()

        # 最深保留层的节点，其子节点将被折叠
        if depth >= max_depth:
            for child in list(node.children):
                label = get_cluster_label(child, sample_cluster)

                # 添加代表子树的叶节点
                node.add_child(name=label, dist=0.0)

                # 删除原子树
                node.remove_child(child)

            continue

        # 继续向下遍历
        for child in reversed(node.children):
            stack.append((child, depth + 1))

    tree.write(outfile=tree_path + ".cut_depth{}".format(max_depth), format=1)

    return tree_path + ".cut_depth{}".format(max_depth)


def visualize_tree_circular(tree_path, rm_ape=True):
    color_dict = {
        "0": "#547A95",  # 灰蓝
        "1": "#E5A06A",  # 柔和橙
        "2": "#8D80B8",  # 灰紫
        "3": "#6FA69A",  # 青绿
        "4": "#D4776A",  # 珊瑚红
        "5": "#C5A65A",  # 金色
    }
    outgroup_name = "8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"

    circle_radius = 3

    # =========================
    # 2. 读取树
    # =========================
    tree = Tree(tree_path, format=1)
    if rm_ape:
        tree.prune(
            [leaf for leaf in tree.get_leaves()
             if leaf.name != outgroup_name],
            preserve_branch_length=False
        )

    # =========================
    # 3. 设置树的样式
    # =========================
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_scale = False
    ts.rotation = 0

    # 分支线条样式
    for node in tree.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        nstyle["vt_line_color"] = "#555555"
        nstyle["hz_line_color"] = "#555555"
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        node.set_style(nstyle)

    # =========================
    # 4. 最外圈标记不同 cluster
    # =========================
    circle_radius = 3

    def layout(node):
        # 只给叶节点添加颜色标记
        if not node.is_leaf():
            return

        # 节点名称格式：cluster_sample
        cluster = node.name.split("_", 1)[0]

        if cluster in color_dict:
            face = CircleFace(
                radius=circle_radius,
                color=color_dict[cluster]
            )

            faces.add_face_to_node(
                face,
                node,
                column=0,
                position="branch-right"
            )

    ts.layout_fn = layout

    # =========================
    # 5. 导出图片
    # =========================
    tree.render(
        tree_path + ".circle.png",
        w=2000,
        units="px",
        tree_style=ts
    )



def visualize_tree_rectangular_ete3(tree_path, rm_ape=True):


    # # STEP: visualize
    # midpoint = tree.get_midpoint_outgroup()
    # tree.set_outgroup(midpoint)
    # node = tree & "8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"
    # node.dist /= 100
    # tree.set_outgroup("8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz")

    outgroup_name = "8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_scale = False

    tree = Tree(tree_path, format=1)

    if rm_ape:


        tree.prune(
            [leaf for leaf in tree.get_leaves()
             if leaf.name != outgroup_name],
            preserve_branch_length=False
        )

    #     tree.write(outfile=tree_path + ".rm_ape", format=1)
    #
    #     tree.render( tree_path + ".rm_ape.png", w=2000, units="px", tree_style=ts)
    #
    # else:
    tree.render( tree_path + ".png", w=2000, units="px", tree_style=ts)


def visualize_tree_rectangular(tree_path, rm_ape=True):


    if tree_path.endswith(".nex"):
        tree = Phylo.read(tree_path, "nexus")
    else:
        tree = Phylo.read(tree_path, "newick")


    # if rm_ape:
    if rm_ape:
        outgroup_name = "8_n1_o0"
        target = next( (clade for clade in tree.find_clades() if clade.name == outgroup_name), None )
        if target is not None:
            tree.prune(target)

        outgroup_name = "8_GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"
        target = next((clade for clade in tree.find_clades() if clade.name == outgroup_name), None)
        if target is not None:
            tree.prune(target)

    # 绘制树
    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(1, 1, 1)

    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True,     branch_labels=lambda clade: (f"{clade.branch_length:.3f}" if clade.branch_length is not None else ""))
    # Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)


    plt.tight_layout()

    if rm_ape:
        plt.savefig(tree_path + ".rmape.png", bbox_inches="tight")
    else:
        plt.savefig(tree_path + ".png", bbox_inches="tight")
        plt.savefig(tree_path + ".svg", bbox_inches="tight")

    plt.close()
    plt.show()

    # ts = TreeStyle()
    # ts.show_leaf_name = True
    # ts.show_scale = False
    #
    # tree = Tree(tree_path, format=1)
    #
    # if rm_ape:
    #
    #
    #     tree.prune(
    #         [leaf for leaf in tree.get_leaves()
    #          if leaf.name != outgroup_name],
    #         preserve_branch_length=False
    #     )
    #
    #     tree.write(outfile=tree_path + ".rm_ape", format=1)
    #
    #
    #     tree.render( tree_path + ".rm_ape.png", w=2000, units="px", tree_style=ts)
    #
    # else:
    #     tree.render( tree_path + ".png", w=2000, units="px", tree_style=ts)

def generate_iTOL_cluster_color_file():

    color_dict = {
        "0": "#547A95",  # 灰蓝
        "1": "#E5A06A",  # 柔和橙
        "2": "#8D80B8",  # 灰紫
        "3": "#6FA69A",  # 青绿
        "4": "#D4776A",  # 珊瑚红
        "5": "#C5A65A",  # 金色
    }

    for spec_chrom in chrom_order:


        chrom_fa_file = os.path.join(work_folder, "{}.hor.fa".format(spec_chrom))
        chrom_annot_file = os.path.join(work_folder, "{}.hor.mafft.fa.treefile.annot.txt".format(spec_chrom))

        with open(chrom_annot_file, "w") as fout, pysam.FastaFile(chrom_fa_file) as fin:
            fout.write("DATASET_COLORSTRIP\n")
            fout.write("SEPARATOR\tTAB\n")
            fout.write("DATASET_LABEL	Cluster\n")
            fout.write("COLOR\t#ff0000\n\n")
            fout.write("STRIP_WIDTH\t40\n")

            fout.write("DATA\n\n")

            for sequence_name in fin.references:

                sample_cluster = sequence_name.split("_")[0]
                if sample_cluster == "8":
                    continue

                sample_color = color_dict[sample_cluster]

                fout.write("{}\t{}\n".format(sequence_name, sample_color))



def assess_flanking_quality_among_cluster(mode="match"):
    CLUSTER_COLORS = {
        "0": "#547A95", "1": "#E5A06A", "2": "#8D80B8",
        "3": "#6FA69A", "4": "#D4776A", "5": "#C5A65A",
    }

    arm = "p"

    flank = "50k"

    spec_cluster = "cluster_2"

    if flank == "all":
        input_folder = "/mnt/g/res_all_flanking/{}-arm".format(arm)
    else:
        input_folder = "/mnt/e/workspace/cpc2_centro/res_all_flanking_{}/{}-arm".format(flank, arm)


    if mode == "length":
        output_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/flanking_{}_{}_{}_length".format(arm, flank, spec_cluster)

    else:
        output_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/flanking_{}_{}_{}".format(arm, flank, spec_cluster)

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    data_matrix = {}
    with open(os.path.join(output_folder, "{}.max_align.stats.txt".format("chraa")), "w") as stats_fout, open(os.path.join(output_folder, "{}.max_align.data.txt".format("chrall")), "w") as data_fout, open(os.path.join(output_folder, "{}.max_align.matrix.txt".format("chrall")), "w") as matrix_fout:
        stats_fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Chromosome", "Cluster", "Count", "Median", "Q1", "Q3", "IQR", "Lower bound", "Upper bound",  "Min", "Max", "Std"))

        for spec_chrom in chrom_order:
            print(spec_chrom)
            chrom_input_folder = os.path.join(input_folder, spec_chrom)
            if spec_cluster == "cluster_2":
                chrom_cluster_file = os.path.join(pha_folder, "{}.features.c2.txt".format(spec_chrom))
            elif spec_cluster in ["cluster_3", "cluster_4"]:
                chrom_cluster_file = os.path.join(pha_folder, "{}.features.c3_c4_c5_c6.txt".format(spec_chrom))
            else:
                print("No this spec cluster", spec_cluster)
                exit()


            # # load sample dict
            sample_cluster_dict = {}
            cluster_order = []
            with open(chrom_cluster_file) as fin:
                for line in fin:
                    if "Sample" in line:
                        continue

                    line_split = line.strip().split()
                    if spec_cluster == "cluster_2":
                        sample, cluster = line_split[0], line_split[-1]
                        if cluster == "-1":
                            continue
                        sample_cluster_dict[sample] = cluster

                        if cluster not in cluster_order:
                            cluster_order.append(cluster)

                    else:
                        sample, cluster_2, cluster_3, cluster_4 = line_split[0], line_split[-5], line_split[-4], line_split[-3]
                        if cluster_2 == "-1":
                            continue

                        if spec_cluster == "cluster_3":
                            if "Q3_Q1" not in cluster_3:
                                continue

                            sample_cluster_dict[sample] = cluster_3

                            if cluster_3 not in cluster_order:
                                cluster_order.append(cluster_3)

                        if spec_cluster == "cluster_4":
                            if "Q3_Q1" not in cluster_4:
                                continue

                            sample_cluster_dict[sample] = cluster_4

                            if cluster_4 not in cluster_order:
                                cluster_order.append(cluster_4)

                    if sample not in data_matrix:
                        data_matrix[sample] = {tmp_chrom: "." for tmp_chrom in chrom_order}



            cluster_order = sorted(cluster_order)

            # # traverse each sample

            collect_dict = {cluster: [] for cluster in cluster_order}


            for sample in sample_cluster_dict:
                cluster = sample_cluster_dict[sample]

                sample_paf = os.path.join(chrom_input_folder, "{}.paf".format(sample))
                max_align_ratio = 0
                max_alt_length = 0
                with open(sample_paf) as fin:
                    for line in fin:
                        line_split = line.strip().split()


                        target_length, align_length = int(line_split[6]), int(line_split[9])

                        align_ratio = align_length / target_length

                        if align_ratio > max_align_ratio:
                            max_align_ratio = align_ratio
                            max_alt_length = int(line_split[3]) - int(line_split[2])

                if mode == "length":
                    collect_dict[cluster].append(max_alt_length)

                    data_matrix[sample][spec_chrom] = max_alt_length

                else:
                    collect_dict[cluster].append(max_align_ratio)

                    data_matrix[sample][spec_chrom] = max_align_ratio

            # # output
            groups = {cluster: np.asarray(collect_dict[cluster], dtype=float) for cluster in collect_dict}
            _, kruskal_p = kruskal(*groups.values())
            kruskal_p = round(kruskal_p, 3) if kruskal_p > 0.01 else format(kruskal_p, '.1e')

            stats_fout.write("{}, Kruskal-Wallis test, P={}\n".format(spec_chrom, kruskal_p))

            with open(os.path.join(output_folder, "{}.max_align.stats.txt".format(spec_chrom)), "w") as chrom_stats_fout:
                chrom_stats_fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("chrom", "cluster", "min", "max", "median", "q1", "q3", "iqr", "lower_bound", "upper_bound", "std", "var"))

                for cluster in collect_dict:

                    values = collect_dict[cluster]

                    data_fout.write("{}\t{}\t{}\n".format(spec_chrom, cluster, "\t".join([str(val) for val in values])))
                    q1 = np.percentile(values, 25)
                    median = round(np.percentile(values, 50), 4)
                    q3 = np.percentile(values, 75)
                    iqr = q3 - q1

                    lower_bound = q1 - 1.5 * iqr
                    upper_bound = q3 + 1.5 * iqr

                    std = round(np.std(values), 4)
                    var = round(np.var(values), 4)

                    min = np.min(values)
                    max = np.max(values)


                    p_values = []
                    for target_cluster in collect_dict:
                        try:
                            p_value = mannwhitneyu(collect_dict[cluster], collect_dict[target_cluster], alternative="two-sided").pvalue
                            p_value = round(p_value, 3) if p_value > 0.01 else format(p_value, '.1e')
                        except ValueError:
                            p_value = "NA"
                        p_values.append(p_value)

                    chrom_stats_fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(spec_chrom, cluster, min, max, median, q1, q3, iqr, lower_bound, upper_bound, std, var, "\t".join([str(p) for p in p_values])))

                    stats_fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(spec_chrom, cluster, len(values), round(median, 4), round(q1, 4), round(q3, 4), round(iqr, 4), round(lower_bound, 4), round(upper_bound, 4), min, max,  std, "\t".join([str(p) for p in p_values])))


            # # plot
            df = pd.DataFrame([ {"cluster": k, "value": x} for k, values in collect_dict.items() for x in values ])

            plt.figure(figsize=(10, 6))

            sns.boxplot( data=df, x="cluster", y="value", showfliers=False, palette=CLUSTER_COLORS)

            plt.tight_layout()
            plt.savefig(os.path.join(output_folder, "{}.max_align.png".format(spec_chrom)))
            plt.close()

        # # output data matrix
        matrix_fout.write("Sample\t{}\n".format("\t".join(chrom_order)))

        for sample in data_matrix:
            matrix_fout.write("{}\t{}\n".format(sample, "\t".join([str(data_matrix[sample][chrom]) for chrom in chrom_order])))


def cigar_to_list(cigar, rm_clip=True):
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    if rm_clip == True:
        if ops[0] == "S" or ops[0] == "H":
            ops = ops[1:]
            lengths = lengths[1:]

        if len(ops) > 1:  # # fix on v1.9 for VACmap aligner, there are aligns with CIGAR like '13520S' and there will be no index -1
            if ops[-1] == "S" or ops[-1] == "H":
                ops = ops[: -1]
                lengths = lengths[: -1]

    return ops, lengths


def merge_flanking_fa(spec_chrom, flank_arm, flank_length, output_folder):

    chr_cluster_file = os.path.join(pha_folder, "{}.features.c2.txt".format(spec_chrom))

    chr_flanking_folder = "/mnt/e/workspace/cpc2_centro/res_all_flanking_{}/{}-arm/{}".format(flank_length, flank_arm, spec_chrom)

    output_fa = os.path.join(output_folder, "all.merged.fa")

    with open(output_fa, "w") as fout:

        # # output ref and ape
        for sample in ["T2T-CHM13.fa.gz", "GCA_028858775.3_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"]:
            sample_fa_file = os.path.join(chr_flanking_folder, "{}.flanking.fa".format(sample))
            fout.write(">{}_{}".format("ref", sample) + "\n")
            with pysam.FastaFile(sample_fa_file) as fin:
                fout.write(fin.fetch(fin.references[0]) + "\n")

            if sample == "T2T-CHM13.fa.gz":
                os.system("cp {} {}".format(sample_fa_file, output_folder))

        # # load cluster
        with open(chr_cluster_file) as fin:
            for line in fin:
                if "Sample" in line:
                    continue

                line_split = line.strip().split()
                sample, cluster = line_split[0], line_split[-1]
                if cluster == "-1":
                    continue

                sample_fa_file = os.path.join(chr_flanking_folder, "{}.flanking.fa".format(sample))
                fout.write(">{}_{}".format(cluster, sample) + "\n")

                with pysam.FastaFile(sample_fa_file) as fin:
                    fout.write(fin.fetch(fin.references[0]) + "\n")



def detect_flank_sv(spec_chrom, flank_arm, flank_length):

    output_path = os.path.join("/mnt/e/workspace/cpc2_centro/res_all_post/flanking_{}_{}_cluster_2/{}".format(flank_arm, flank_length, spec_chrom))
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # # generate merged fa
    merge_flanking_fa(spec_chrom, flank_arm, flank_length, output_path)

    # # align
    merged_fa = os.path.join(output_path, "all.merged.fa")
    ref_fa = os.path.join(output_path, "T2T-CHM13.fa.gz.flanking.fa")
    paf_file = merged_fa + ".realign.paf"

    os.system("minimap2 -x asm20 -c --eqx {} {} > {}".format(ref_fa, merged_fa, paf_file))

    # # parse paf for sv
    sv_file = merged_fa + ".realign.paf.sv.txt"

    with open(sv_file, "w") as fout, open(paf_file) as fin:

        for line in fin:
            line_split = line.strip().split()

            sample, ref_start, read_start, cigar = line_split[0], int(line_split[7]), int(line_split[2]), line_split[-1].replace("cg:Z:", "")
            # print(sample, cigar)
            cigar_ops, cigar_lengths = cigar_to_list(cigar)

            for i in range(len(cigar_ops)):
                op, length = cigar_ops[i], cigar_lengths[i]

                if op in ["=", "X"]:
                    ref_start += length
                    read_start += length
                elif op == "I":
                    if length >= 50:
                        fout.write("{}\t{}\t{}\t{}\t{}\n".format(sample, ref_start, read_start, op, length))
                    read_start += length

                elif op == "D":
                    if length >= 50:
                        fout.write("{}\t{}\t{}\t{}\t{}\n".format(sample, ref_start, read_start, op, length))
                    ref_start += length
                else:
                    print("ERROR, no this cigar", op, length)
                    exit()







def assess_flanking_tree_monophyletic(tree_file):

    tree = Tree(tree_file, format=0)

    MAX_FOREIGN_RATIO = 0.0

    def get_cluster(name):
        return name.split("_", 1)[0].strip()

    # 1. 统计各 cluster 的总样本数
    cluster_total = defaultdict(int)

    for leaf in tree.iter_leaves():
        cluster_total[get_cluster(leaf.name)] += 1

    # 2. 统计每个节点的后代 cluster 数量
    node_counts = {}

    for node in tree.traverse("postorder"):
        counts = defaultdict(int)

        if node.is_leaf():
            counts[get_cluster(node.name)] = 1
        else:
            for child in node.children:
                for cluster, n in node_counts[child].items():
                    counts[cluster] += n

        node_counts[node] = dict(counts)

    # 3. 寻找每个 cluster 的最大容错分支
    best_clade = {
        cluster: {
            "target_count": 0,
            "foreign_count": 0,
            "purity": 0.0,
        }
        for cluster in cluster_total
    }

    for node, counts in node_counts.items():
        n_descendants = sum(counts.values())

        for cluster, target_count in counts.items():
            foreign_count = n_descendants - target_count
            foreign_ratio = foreign_count / n_descendants

            # 要求异类样本比例不超过 5%
            if foreign_ratio > MAX_FOREIGN_RATIO:
                continue

            purity = target_count / n_descendants
            old = best_clade[cluster]

            # 优先最大化目标 cluster 的覆盖样本数
            # 覆盖数相同时，优先选择纯度更高的分支
            if (
                    target_count > old["target_count"]
                    or (
                    target_count == old["target_count"]
                    and purity > old["purity"]
            )
            ):
                best_clade[cluster] = {
                    "target_count": target_count,
                    "foreign_count": foreign_count,
                    "purity": purity,
                }

    # 4. 输出统计结果
    results = []

    for cluster, total in sorted(
            cluster_total.items(),
            key=lambda x: str(x[0])
    ):
        if cluster == "8":
            continue

        best = best_clade[cluster]
        coverage = best["target_count"] / total

        results.append({
            "cluster": cluster,
            "n_samples": total,
            "max_clade_samples": best["target_count"],
            "foreign_samples": best["foreign_count"],
            "coverage": coverage,
            "purity": best["purity"],
        })

        print(os.path.basename(tree_file), cluster, total, best["target_count"], round(coverage, 3), coverage == 1.0, sep="\t")


def assess_chr4_HSat_LINE():
    def merge_intervals(intervals, merge_threshold=10000):
        if not intervals:
            return []

        intervals = sorted(intervals, key=lambda x: x[0])
        merged = [list(intervals[0])]

        for start, end in intervals[1:]:
            if start - merged[-1][1] <= merge_threshold:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])

        return merged

    chr4_cluster_file = os.path.join(pha_folder, "{}.features.c2.txt".format("chr4"))
    chr4_annot_folder = os.path.join(annot_folder, "chr4")

    type_list = ["LINE", "SINE", "HSat1A"]

    # # load cluster
    sample_cluster_dict = {}
    sample_length_dict = {}
    with open(chr4_cluster_file) as fin:
        for line in fin:
            if "Sample" in line:
                continue

            line_split = line.strip().split()
            sample, cluster, aSat_length = line_split[0], line_split[-1], line_split[-3]
            if cluster == "-1":
                continue

            sample_cluster_dict[sample] = cluster
            sample_length_dict[sample] = aSat_length

    print("Sample\tCluster\taSat\t{}".format("\t".join(type_list)))

    for sample in sample_cluster_dict:
        sample_region_dict = {type: [] for type in type_list}

        with open(os.path.join(chr4_annot_folder, "{}#stv_hor_osat_te.pha.full.bed".format(sample))) as fin:

            for line in fin:
                line_split = line.strip().split()
                start, end, type = int(line_split[1]), int(line_split[2]), line_split[3]

                if type not in type_list:
                    continue

                sample_region_dict[type].append([start, end])


        sample_merged_region_dict = {type: [] for type in type_list}
        for type in sample_region_dict:
            sample_merged_region_dict[type] = merge_intervals(sample_region_dict[type])

        print("{}\t{}\t{}".format(sample, sample_cluster_dict[sample], sample_length_dict[sample]), end="")

        for type in sample_merged_region_dict:


            print("\t", end="")
            if len(sample_merged_region_dict[type]) == 0:
                print(".", end="")
            else:
                print(";".join(["{},{}".format(start, end) for start, end in sample_merged_region_dict[type]]), end="")
        print()


def assess_chr12_SINE_LINE():
    def merge_intervals(intervals, merge_threshold=10000):
        if not intervals:
            return []

        intervals = sorted(intervals, key=lambda x: x[0])
        merged = [list(intervals[0])]

        for start, end in intervals[1:]:
            if start - merged[-1][1] <= merge_threshold:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])

        return merged

    chr_cluster_file = os.path.join(pha_folder, "{}.features.c2.txt".format("chr12"))
    chr_annot_folder = os.path.join(annot_folder, "chr12")

    type_list = ["LINE"]

    # # load cluster
    sample_cluster_dict = {}
    sample_length_dict = {}
    sample_canonical_dict = {}
    sample_variant1_dict = {}
    sample_variant2_dict = {}
    sample_variant3_dict = {}

    with open(chr_cluster_file) as fin:
        for line in fin:
            if "Sample" in line:
                continue

            line_split = line.strip().split()
            sample, cluster, aSat_length, canonical, variant1, variant2, variant3 = line_split[0], line_split[-1], line_split[-3], float(line_split[3]), float(line_split[4]), float(line_split[5]), float(line_split[6])
            if cluster == "-1":
                continue

            sample_cluster_dict[sample] = cluster
            sample_length_dict[sample] = aSat_length
            sample_canonical_dict[sample] = canonical
            sample_variant1_dict[sample] = variant1
            sample_variant2_dict[sample] = variant2
            sample_variant3_dict[sample] = variant3


    print("Sample\tCluster\taSat\tCanonical\tVariant1\tVariant2\tVariant3\t{}".format("\t".join(type_list)))

    for sample in sample_cluster_dict:
        sample_region_dict = {type: [] for type in type_list}

        with open(os.path.join(chr_annot_folder, "{}#stv_hor_osat_te.pha.full.bed".format(sample))) as fin:

            for line in fin:
                line_split = line.strip().split()
                start, end, type = int(line_split[1]), int(line_split[2]), line_split[3]

                if type not in type_list:
                    continue

                sample_region_dict[type].append([start, end])


        sample_merged_region_dict = {type: [] for type in type_list}
        for type in sample_region_dict:
            sample_merged_region_dict[type] = merge_intervals(sample_region_dict[type], merge_threshold=1000)

        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sample, sample_cluster_dict[sample], sample_length_dict[sample], sample_canonical_dict[sample], sample_variant1_dict[sample], sample_variant2_dict[sample], sample_variant3_dict[sample]), end="")

        for type in sample_merged_region_dict:


            # print("\t{}\t{}".format(sum([ end - start for start, end in sample_merged_region_dict[type]]), len(sample_merged_region_dict[type])), end="")
            print("\t{}".format(len(sample_merged_region_dict[type])), end="")

        print()


def assess_chr9_dup():
    target = "S2C9H1L.1-3_7"

    chr9_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/bin_size_100_principal_hor_frequency/chr9"

    for file in os.listdir(chr9_folder):

        if not file.endswith(".bin_100.txt"):
            continue

        file_name = file.replace(".bin_100.txt", "")

        sample_name = ".".join(file_name.split(".")[1: ])

        sample_cluster = file_name.split(".")[0]
        if sample_cluster == "-1":
            continue

        sample_bin_file = os.path.join(chr9_folder, file)

        bin_indices = []

        with open(sample_bin_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                fields = line.split("\t")

                # bin_index
                bin_index = int(fields[0])

                # 每个特征占两列：count 和 proportion
                # 前4列是基础信息
                # 第5-6列：S2C9H1L.1-7
                # 第7-8列：S2C9H1L.1-2_6-7
                # 第9-10列：S2C9H1L.1-3_7
                count = int(fields[6].split(":")[-1].split(",")[0])

                if count == 0:
                    bin_indices.append(bin_index)

        # 合并连续的 bin_index
        regions = []
        if bin_indices:
            start = prev = bin_indices[0]

            for idx in bin_indices[1:]:
                if idx == prev + 1:
                    prev = idx
                else:
                    regions.append((start, prev))
                    start = prev = idx

            regions.append((start, prev))

        # 输出格式：
        # 区间数量 \t 总 bin 数量 \t 各区间信息
        total_bins = sum(end - start + 1 for start, end in regions)

        region_info = [
            f"{start}-{end}" if start != end
            else f"{start}"
            for start, end in regions
        ]
        print(sample_name, sample_cluster, len(regions), total_bins, ";".join(region_info), sep="\t")
        # print(
        #     f"{len(regions)}\t{total_bins}\t"
        #     + "\t".join(region_info)
        # )

def assess_chr21_insdel():


    cigar_string = "105=1X78=1X917=1X358=1X438=1X24=1X216=1X250=1X280=1X1062=1X797=1X216=1X38=1X112=1X1067=2D1040=1X497=1X164=1X204=1X235=1X505=1X295=1X101=1X426=1X295=1X161=508I186=1X112=1X67=3X19=1X2=1X2=3X4=1X2=2X3=1X1=1X1=1X3=1X1=1X2=1X3=1X2=1X4=1X7=1X2=1X1=1X621=1X248=1D203=1D144=1X243=1X1092=1X292=1X123=1X177=2X3=1X551=1X13=1X479=1X402=1X477=1X139=1X331=1X36=1X74=1X25=1X229=1X358=1X145=1X1239=1X720=1X175=1X80=1X181=1X626=1X145=1X35="


    ops, lengths = cigar_to_list(cigar_string)

    for i in range(len(ops)):

        op, length = ops[i], lengths[i]

        if op in ["I", "D"] and length >= 50:
            print(op, length)



def assess_chr1_chr5_chr19():

    # for chrom in ["chr1", "chr5", "chr19"]:

    chrom = "chr1"
    chrom_pha_file = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier/{}.features.c2.txt".format(chrom)

    cluster_dict = {}

    with open(chrom_pha_file) as chrom_fin:
        for line in chrom_fin:
            if "Sample" in line:
                continue

            line_split = line.strip().split()

            sample, cluster = line_split[0], line_split[-1]
            if cluster == "-1":
                continue

            if cluster not in cluster_dict:
                cluster_dict[cluster] = []

            sample_annot_file = "/mnt/e/workspace/cpc2_centro/res_all_annot/{}/{}#stv_hor.bed".format(chrom, sample)


            count = 0
            with open(sample_annot_file) as sample_fin:
                for hor_line in sample_fin:
                    hor_type = hor_line.strip().split()[3]

                    if hor_type == "S1C1/5/19H1L.6/4_5":
                        count += 1
                    else:
                        match = re.search(r'\(_6/4_5\)\{(\d+)\}', hor_type)

                        if match:
                            count += int(match.group(1))
            cluster_dict[cluster].append(count)
            # print(sample, cluster, count, sep="\t")

    all_values = []
    for cluster in cluster_dict:

        values = cluster_dict[cluster]
        all_values.extend(values)

        q1 = np.percentile(values, 25)
        median = round(np.percentile(values, 50), 4)
        q3 = np.percentile(values, 75)
        iqr = q3 - q1

        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr


        print(chrom, cluster, median, q1, q3, iqr, lower_bound, upper_bound, sep="\t")

    q1 = np.percentile(all_values, 25)
    median = round(np.percentile(all_values, 50), 4)
    q3 = np.percentile(all_values, 75)
    iqr = q3 - q1

    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    print(chrom, "all", median, q1, q3, iqr, lower_bound, upper_bound, sep="\t")


def compare_p_q_tree(p_tree, q_tree):

    tree_p = Tree(p_tree, format=1)
    tree_q = Tree(q_tree, format=1)

    # 必须使用相同的样本集合
    names_p = set(tree_p.get_leaf_names())
    names_q = set(tree_q.get_leaf_names())

    if names_p != names_q:
        raise ValueError(
            f"Sample sets differ: "
            f"p-only={names_p - names_q}, "
            f"q-only={names_q - names_p}"
        )

    # unrooted RF：忽略根的位置
    rf_result = tree_p.robinson_foulds(
        tree_q,
        unrooted_trees=True,
        expand_polytomies=False
    )


    rf, max_rf = rf_result[0], rf_result[1]


    normalized_rf = rf / max_rf if max_rf else 0.0

    print(f"RF distance: {rf}")
    print(f"Maximum RF distance: {max_rf}")
    print(f"Normalized RF distance: {normalized_rf:.4f}")

if __name__ == '__main__':


    chrom_order = ["chr1", "chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr16", "chr17", "chr19", "chr20", "chr21"]
    # chrom_order = ["chr12", "chr13", "chr16", "chr17", "chr19", "chr20", "chr21"]
    #
    # chrom_order = ["chr10"]

    work_folder = \
        "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier/cluster_2-flank_q_20k-random-nall-ape/"

    annot_folder = "/mnt/e/workspace/cpc2_centro/res_all_annot"
    pha_folder = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier"


    # # STEP: plot distance_matrix
    process_pool = multiprocessing.Pool(23)
    for spec_chrom in chrom_order:
        process_pool.apply_async(generate_tree_distance_matrix, (spec_chrom,))
        # process_pool.apply_async(generate_mafft_distance_matrxi, (spec_chrom,))

        # plot_tree_distance_matrix(spec_chrom)
    process_pool.close()
    process_pool.join()


    # # STEP: assess chr1 and chr4
    assess_chr1_inversion_cluster()
    assess_chr4_HSat_LINE()
    assess_chr12_SINE_LINE()
    assess_chr9_dup()
    assess_chr21_insdel()
    assess_chr1_chr5_chr19()


    # # STEP:
    # generate_iTOL_cluster_color_file()


    # # # compare two_trees
    # old_res = "/mnt/e/workspace/cpc2_centro/res_all_post/principal_hor_analysis_outlier/cluster_2-flank_q_20k-random-nall-ape-old/"
    # for spec_chrom in chrom_order:
    #     print(spec_chrom)
    #     nwk_tree1 = os.path.join(work_folder, "{}.flank.mafft.fa.treefile".format(spec_chrom))
    #
    #     nwk_tree2 = os.path.join(old_res, "{}.flank.mafft.fa.treefile".format(spec_chrom))
    #     compare_p_q_tree(nwk_tree1, nwk_tree2)

    # # STEP: cut tree
    for spec_chrom in chrom_order:
        print(spec_chrom)

        nex_tree = os.path.join(work_folder, "{}.flank.mafft.fa.treefile.lsd.timetree.nex".format(spec_chrom))
        collapsed_nex_tree = collapse_tree_by_cluster(nex_tree)
        visualize_tree_rectangular(collapsed_nex_tree, rm_ape=True)
        visualize_tree_rectangular(collapsed_nex_tree, rm_ape=False)


        nwk_tree = os.path.join(work_folder, "{}.flank.mafft.fa.treefile".format(spec_chrom))
        collapsed_nwk_tree = collapse_tree_by_cluster(nwk_tree)
        visualize_tree_rectangular(collapsed_nwk_tree, rm_ape=True)
        visualize_tree_rectangular(collapsed_nwk_tree, rm_ape=False)


        cut_nwk_tree = cut_tree_by_depth(nwk_tree, max_depth=4)
        # visualize_tree_rectangular(cut_nwk_tree, rm_ape=True)
        visualize_tree_rectangular_ete3(cut_nwk_tree, rm_ape=True)


    # # STEP: for flanking
    assess_flanking_quality_among_cluster()

    for chrom in chrom_order:
        detect_flank_sv(chrom, "p", "20k")


    # # STEP: assess monophyletic
    for spec_chrom in chrom_order:
        assess_flanking_tree_monophyletic(os.path.join(work_folder, "{}.flank.mafft.fa.treefile".format(spec_chrom)))



