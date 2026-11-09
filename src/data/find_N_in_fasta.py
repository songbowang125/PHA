import multiprocessing
import os
import pysam
import re


def process_one_fasta(project, input_folder, sample_fasta):

    print(sample_fasta)
    with open("tmp_N_{}".format(sample_fasta), "w") as fout:
        with pysam.FastaFile(os.path.join(input_folder, sample_fasta)) as fin:

            for contig in sorted(fin.references):

                contig_seq = fin.fetch(contig)

                # regions = [(m.start(), m.end()) for m in re.finditer(r'N+', contig_seq)]
                regions = [(m.start(), m.end()) for m in  re.finditer(r'[^Nn]+', contig_seq)]

                if len(regions) == 0:
                    continue

                for start, end in regions:
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(project, sample_fasta, contig, start, end, end - start))
                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(project, sample_fasta, contig, start, end, end - start))



if __name__ == '__main__':
    candidate_chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    projects =  ["cpc2", "hgsvc3", "hprc2", "ref", "ape", "apg1", "apr", "sea3k1", "tbn"]

    process_pool = multiprocessing.Pool(40)
    for project in projects:
        input_folder = "/home/data/user/songbo/workspace/pangenome-compare/data_{}".format(project)

        for sample_fasta in os.listdir(input_folder):
            if not sample_fasta.endswith(".gz"):
                continue

            # process_one_fasta(project, input_folder, sample_fasta)
            process_pool.apply_async(process_one_fasta, (project, input_folder, sample_fasta))

    process_pool.close()
    process_pool.join()

    with open("data_non-N_regions.txt", "w") as fout:
        for project in projects:
            input_folder = "/home/data/user/songbo/workspace/pangenome-compare/data_{}".format(project)

            for sample_fasta in os.listdir(input_folder):
                if not sample_fasta.endswith(".gz"):
                    continue

                sample_N_file = "tmp_N_{}".format(sample_fasta)

                with open(sample_N_file) as fin:
                    for line in fin:
                        fout.write(line)


