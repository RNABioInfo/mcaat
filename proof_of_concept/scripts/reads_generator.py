from Bio import SeqIO
import random
import sys

#write a read generator that takes a fasta file and generates reads from it of length 125  with 10x coverage
def generate_reads(input_file, output_file, read_length, coverage):
    with open(input_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    count = 0
    id = records[0].id
    sequence = records[0].seq
    print(len(sequence))
    # cut the sequence into reads of length 125
    # generate 10x coverage
    with open(output_file, "w") as handle:
        for i in range(coverage):
            for j in range(len(sequence)-read_length):
                count += 1
                handle.write("@read_"+str(count)+"\n")
                handle.write(str(sequence[j:j+read_length])+"\n")
                handle.write("+\n")
                handle.write("I"*read_length+"\n")
    print(count)

if __name__ == '__main__':
    input_file = "../data/"+sys.argv[1]+"/genome/"+sys.argv[1]+".fasta"
    output_file = "../data/"+sys.argv[1]+"/reads/"+sys.argv[2]+".fastq"
    read_length = 125
    coverage = 10
    generate_reads(input_file, output_file, read_length, coverage)
    