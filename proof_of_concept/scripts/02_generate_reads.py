from Bio import SeqIO
import sys
import random
file = "../data/" + sys.argv[1]

def split_genome_125(genome, k=125):
    """Split genome into k-mers."""
    k_mers = []
    for i in range(0, len(genome) - k + 1):
        k_mers.append(genome[i:i+k])
    #random.shuffle(k_mers)
    return k_mers

def read_genome(genome_file):
    """Read in genome file."""
    genome = []
    for record in SeqIO.parse(genome_file, "fasta"):
        genome.append(str(record.seq))
    return genome

def write_reads(k_mers):
    """Write k-mers to file."""
    with open(file+"/reads/R1.fastq", "w") as f:
        for i in range(len(k_mers)):
            f.write("@"+str(i) +"\n"+
                    k_mers[i] + "\n" +
                    "+" + "\n" +
                    "I"*len(k_mers[i]) + "\n")

def create_reads():
    genome = read_genome(file+"/genome/"+sys.argv[1]+".fasta")
    k_mers = split_genome_125(genome[0])
    #print(genome[0][0:125])
    #print(k_mers[0])
    write_reads(k_mers)


#read from input 

create_reads()