# read folder names from ../proof_of_concept/data/ 
# if it starts with _ then skip
# else add to list


import os

folder = "../data/"
names = []
for filename in os.listdir(folder):
    if filename[0] != "_" and filename != "Genomes" and filename!="graphviz.gv":
        names.append(filename)

#for each folder/filename/genome create a .txt file with the genome name
def create_libs():
    for filename in names:
        #if (folder + "/" + filename+"/genome/"+filename+".lib") exists remove it
        if os.path.exists(folder + "/" + filename+"/genome/"+filename+".lib"):
            os.remove(folder + "/" + filename+"/genome/"+filename+".lib")
        with open(folder + "/" + filename+"/genome/"+filename+".lib", "w") as f:
            f.write("# "+filename+"\nse "+"/vol/d/development/git/mCAAT/proof_of_concept/data/"+ filename+"/genome/"+filename+".fasta")
            f.close()
    print("Created all libraries")

def build_libs():
    for filename in names:
        os.system("./megahit_core buildlib "+folder + filename+"/genome/"+filename+".lib "+folder + "/" + filename+"/genome/"+filename)
    print("Built all libraries")

def rename_libs():
    
    for filename in names:
        if os.path.exists(folder + "/" + filename+"/genome/"+filename+".info"):
            os.remove(folder + "/" + filename+"/genome/"+filename+".info")
        os.rename(folder + filename+"/genome/"+filename+".lib_info", folder + filename+"/genome/"+filename+".fasta.info")
    print("Renamed all .lib_info files to .info")

def build_graph():
    for filename in names:
        os.system("./megahit_core seq2sdbg --host_mem 150000000000 -k 23 --contig "+folder + filename+"/genome/"+filename+".fasta -o"+folder + filename+"/genome/graph/")
    print("Built all graphs")
build_graph()
print(names)
    