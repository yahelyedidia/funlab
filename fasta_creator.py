import pandas as pd
from Bio import SeqIO


MATRIX = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/sites_groups/stable_met.tsv"
FASTA = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/hg38.analysisSet.chroms/chr{0}.fa"
OUTPUT = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/seq_stable_met.fasta"


def singal_action(matrix, output_file, cnum):
    fasta = SeqIO.read(FASTA.format(cnum), "fasta")
    temp_matrix = matrix[matrix["chr"] == "chr{0}".format(cnum)]
    for site in temp_matrix.iterrows():
        title = "chr{0}|id {1}|start {2}|end {3}".format(cnum, site[1]["Unnamed: 0"], site[1]["start"], site[1]["end"])
        seq = fasta.seq[site[1]["start"]: site[1]["end"]]
        output_file.write(">{0}\n{1}\n".format(title, seq))


def main():
    matrix = pd.read_csv(MATRIX, sep="\t")
    output_file = open(OUTPUT, "w")
    for cnum in range(1, 23):
        singal_action(matrix, output_file, cnum)
    singal_action(matrix, output_file, "X")
    singal_action(matrix, output_file, "Y")
    output_file.close()
    print("done")


main()
