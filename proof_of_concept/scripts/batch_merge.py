# pick all the file names from the folder and run .04_merge.py $file_name

# Path: development/git/mCAAT/proof_of_concept/scripts/compare_to_ref.py

import sys
import os

#list all the files in the folder]
folder = "../data/"
for filename in os.listdir(folder):
    print(filename)
    os.system("python3 04_merge.py " + filename)
    print("Done with", filename)
    print("")
