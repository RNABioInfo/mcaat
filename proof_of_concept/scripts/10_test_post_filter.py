import sys
import os
folder = "../data/"
for filename in os.listdir(folder):
    print(filename)
    os.system("python3 09_post_filtering_i.py " + filename)
    print("Done with", filename)
    print("")
