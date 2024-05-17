import os
import random
import sys

t = os.listdir("../data/")
count = 0
list_of_dirs = []
for folder in t:
    if "_" in folder or "Genomes" in folder or "metagenome" in folder:
        continue
    list_of_dirs.append(folder)

print(list_of_dirs)
folder_list = random.sample(list_of_dirs,int(sys.argv[2]))

for folder in folder_list:
    #copy files from each of the folders to a sys.argv[1] folder
    os.system("cp -r ../data/"+folder+"/genome/"+folder+".fasta"+" ../data/"+sys.argv[1]+"/genome/")

