import sys
import os

folder = "../data/"
for filename in os.listdir(folder):
    print(filename)
    os.system("python3 05_compare_to_ref.py " + filename)
    print("Done with", filename)
    print("")