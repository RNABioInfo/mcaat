import os

folder = "proof_of_concept/data/"
for filename in os.listdir(folder):
    print(filename)
    os.system("./mcaat " + filename)
    print("Done with", filename)
    print("")
