import sys

folder = "/vol/d/development/git/mCAAT/proof_of_concept/data/"+sys.argv[1]+"/cycles_genome/"
file = "fp_before_filter.txt"
def find_second_occurence(txt, str1):
  return txt.find(str1, txt.find(str1)+1)

def read_file(folder, file):
    output_list = []
    with open(folder+file, "r") as f:
        output_list = f.readlines()
    for i in range(len(output_list)):
        output_list[i] = output_list[i][1:-1]
        a = find_second_occurence(output_list[i],output_list[i][0:23])
        output_list[i] = output_list[i][0:a]
    output_list = output_list[2:]
    print(output_list[0:10])
    #write resulting list to file
    output_list.sort()
    with open(folder + "FP_"+sys.argv[1]+".fasta", "w") as f:
        for line_index in range(len(output_list)):
            f.write(">" + str(line_index) + "\n")
            f.write(output_list[line_index] + "\n")
    return output_list

read_file(folder, file)