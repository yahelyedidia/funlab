import sys
import os
import pandas as pd
from Bio import SeqIO

MOTIV_SIZE = 200

MATRIX = "{0}/{1}.tsv"
FASTA = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/hg38.analysisSet.chroms/chr{0}.fa"
OUTPUT = "{0}/fastas/{1}.fasta"


def singal_action(matrix, output_file, cnum):
    """
    create fasta format from one chr data only.
    :param matrix: the relevant CTCF binding sites matrix
    :param output_file: fata file to write into
    :param cnum: the chromosome number

    """
    fasta = SeqIO.read(FASTA.format(cnum), "fasta")
    temp_matrix = matrix[matrix["chr"] == "chr{0}".format(cnum)]
    for site in temp_matrix.iterrows():
        title = "chr{0}|id{1}|start{2}|end{3}".format(cnum, site[1]["Unnamed: 0"], site[1]["start"], site[1]["end"])
        seq = fasta.seq[site[1]["start"]: site[1]["end"]]
        if len(seq) > MOTIV_SIZE:
            y = (len(seq) - MOTIV_SIZE) // 2
            seq = seq[y: y + 200]
        output_file.write(">{0}\n{1}\n".format(title, seq))


def fasta_creator(dir, file_name):
    """
    creates a new fasta file with sequence of CTCF binding sites by matrix
    :param matrix_name: the matrix name to open
    :param output_name:  the new output file name.
    """
    matrix = pd.read_csv(MATRIX.format(dir, file_name), sep="\t")
    output_file = open(OUTPUT.format(dir, file_name), "w+")
    for cnum in range(1, 23):
        singal_action(matrix, output_file, cnum)
    singal_action(matrix, output_file, "X")
    singal_action(matrix, output_file, "Y")
    output_file.close()
    print("done")


def convert_dir(dir):
    for file in os.listdir(dir):
        if file.endswith(".tsv"):
            fasta_creator(dir, file[:-4])


if __name__ == '__main__':
    """
    the arguments files: [1] the matrix name, [2] the output name.
    """
    dir_name = sys.argv[1]
    convert_dir(dir_name)
    # fasta_creator(matrix_name, output_name)
