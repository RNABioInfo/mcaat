import os

folder = "../data/"
names = []

#add all folders to the list
for filename in os.listdir(folder):
    if filename[0] != "_" and filename != "Genomes" and filename!="graphviz.gv":
        names.append(filename)

#rename all .lib_info files to .info if you want to utilise a true genome->sdbg graph
def rename_libs():
    for filename in names:
        os.rename(folder + filename+"/genome/"+filename+".fasta.info", folder + filename+"/genome/"+filename+".lib.info")
    print("Renamed all .lib_info files to .info")

#for each folder/filename/genome create a .txt file with the genome name
def create_libs():
    for filename in names:
        if os.path.exists(folder + "/" + filename+"/genome/"+filename+".lib"):
            os.remove(folder + "/" + filename+"/genome/"+filename+".lib")
        with open(folder + "/" + filename+"/genome/"+filename+".lib", "w") as f:
            f.write("# "+filename+"\nse "+"/vol/d/development/git/mCAAT/proof_of_concept/data/"+ filename+"/genome/"+filename+".fasta")
            f.close()
    print("Created all libraries")

#select the command to execute based on the command parameter
def select_megahit_command(command:bool,name:str) -> str:
    if command:
        return "./megahit_core read2sdbg --host_mem 150000000000 -m 1 --num_cpu_threads 24 -k 23 --read_lib_file "+folder + name+"/genome/"+name+" --o "+folder + name+"/genome/graph/graph")
    return "./megahit_core buildlib "+folder + name+"/genome/"+name+".lib "+folder + "/" + name+"/genome/"+name

#execute the build_lib function for each folder/filename/genome
def build_libs()-> None:
    for filename in names:
        os.system(select_megahit_command(False,filename))
    print("Built all libraries")

def build_graph(name:str="") -> None:
    if name != "":
        os.system(select_megahit_command(True,name))
        return
    
    for filename in names:
        os.system(select_megahit_command(True,filename))
    

build_graph()
