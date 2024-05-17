import os, sys

#read all the files in ..data/sys.argv[1]/genome folder
files = os.listdir("../data/"+sys.argv[1]+"/genome/")
files_names = []

for file in files:
    if "metagenome" in file:
        continue
    files_names.append(file[:-6])

for file in files_names:
    folder_for_reference = "../data/Genomes/CRISPR_Seq/" + file[:-2] +"/"
    if not os.path.isdir(folder_for_reference):
        folder_for_reference = "../data/" + file +"/"+"reference/"
    os.system("cp -r "+folder_for_reference+" ../data/"+sys.argv[1]+"/reference/")

print(files_names)
