import os
import sys

def merge(folder):
    folder = "../data/"+folder+"/cycles/"
    #folder = "../data/"+sys.argv[1]+"/cycles/"
    file = folder + "str_paths.txt"
    #if file exists, delete it
    try:
        os.remove(folder+"results.txt")
    except OSError:
        pass
    #create true_positives.txt
    open(folder+"results.txt", "w").close()

    with open(file, "r") as f:
        lines = f.readlines()
        true_positives = []
        for line in lines:
            if line.startswith("F"):
                continue
            words = line.split()
            #print(len(words))
            first_word = words[0]
            last_chars = ""
            for i in range(1,len(words)):
                last_chars += words[i][-1]
            true_positives.append(first_word+last_chars)

    true_positives = list(set(true_positives))
    true_positives.sort()
    #true_positives.sort()
    with open(folder + "results.txt", "a") as f:
        for element in true_positives:
            f.write(element+"\n")

    print(len(true_positives))

merge(sys.argv[1])