import os
#read text line by line
lines = []

with open("/vol/d/development/git/mCAAT/proof_of_concept/data/data.txt", "r") as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines]
print(lines[5])

report = ""
for line in lines:
    if line.startswith("#"): continue
    os.system("./mcaat " + line)
    print("Done with", line)
    print("")

