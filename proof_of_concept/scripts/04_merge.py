import os
import sys
folder = "../data/"+sys.argv[1]+"/cycles/"
file = folder + "str_true_positives.txt"
#if file exists, delete it
try:
    os.remove(folder+"true_positives_after_filter.txt")
except OSError:
    pass
#create true_positives.txt
open(folder+"true_positives_after_filter.txt", "w").close()

with open(file, "r") as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith("F"):
            continue
        words = line.split()
        first_word = words[0]
        last_chars = ""
        for i in range(1,len(words)):
            last_chars += words[i][-1]
        with open(folder + "true_positives_after_filter.txt", "a") as f2:
            f2.write(first_word+last_chars+last_chars+"\n")
