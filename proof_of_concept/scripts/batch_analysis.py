import sys
import os

lines = []

with open("/vol/d/development/git/mCAAT/proof_of_concept/scripts/folders.txt", "r") as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines]
print(lines[5])

report = ""

for line in lines:
    if line.startswith("#"): continue
    os.system("python3 05_compare_to_ref.py " + line + " >> " "no_filter_report.txt")
    print("Done with", line)
    print("")
