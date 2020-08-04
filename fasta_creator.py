import pandas as pd
from Bio import SeqIO


MATRIX = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/the_big_matrix.tsv"
FASTA = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/hg38.analysisSet.chroms/chr1.fa"

def main():
    matrix = pd.read_csv(MATRIX, sep="\t")
    fasta = SeqIO.read(FASTA, "fasta")
    matrix = matrix[matrix["chr"] == "chr1"]
    with open("example.fasta", "w") as output_handle:
        for site in matrix.iterrows():
            # temp_seq = Seq(fasta.seq[site[1]["start"]: site[1]["end"]])
            # SeqIO.write(temp_seq, output_handle, "fasta")
            title = "chr{0}|id {1}|start {2}|end {3}".format(1, site[1]["Unnamed: 0"], site[1]["start"], site[1]["end"])
            seq = fasta.seq[site[1]["start"]: site[1]["end"]]
            output_handle.write(">{0}\n{1}\n".format(title, seq))

    print(matrix.head())

main()
