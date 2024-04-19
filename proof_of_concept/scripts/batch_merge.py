# pick all the file names from the folder and run .04_merge.py $file_name

# Path: development/git/mCAAT/proof_of_concept/scripts/compare_to_ref.py

import sys
import os

#list all the files in the folder]
lines = []

with open("/vol/d/development/git/mCAAT/proof_of_concept/data/data.txt", "r") as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines]
print(lines[5])

report = ""
for line in lines:
    if line.startswith("#"): continue
    os.system("python3 04_merge.py " + line)
    print("Done with", line)
    print("")