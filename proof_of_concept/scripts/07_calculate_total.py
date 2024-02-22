'''
AP019695.1
- Spacers in genome: 36
- Elements in cycles: 44
- False positives: 8
- Match: 36
- Match percentage: 100.0 

AP019711.1
- Spacers in genome: 26
- Elements in cycles: 30
- False positives: 4
- Match: 26
- Match percentage: 100.0
'''

# above is the format of input file for this script
# pick lines that start with spacers in genome and read only integer after the colon
# pick lines that start with match and read only integer after the colon
# divide match by spacers in genome and multiply by 100
import sys

def calculate():
    before_after_filter = sys.argv[1]
    with open("../new_method_"+before_after_filter, "r") as f:
        lines = f.readlines()
        
    spacers = 0
    matched = 0
    false_positives = 0
    counter=0
    matched_percentage = 0
    
    for line in lines:
        if line.startswith("- Spacers in genome:"):
            spacers += int(line.split(":")[1].strip())
            counter+=1
        elif line.startswith("- Match:"):
            matched += int(line.split(":")[1].strip())
        elif line.startswith("- False positives:"):
            false_positives += int(line.split(":")[1].strip())
        elif line.startswith("- Match percentage:"):
            matched_percentage += float(line.split(":")[1].strip())
    
    print("Results for", counter, "genomes", ":")
    print("-- Spacers in genome:", spacers)
    print("-- Matched in cycles:", matched)
    print("-- False positives:", false_positives)
    print("-- Match percentage:", matched/spacers*100)
    print("-- CRISPR System match average:", matched_percentage/counter)

if __name__ == "__main__":
    calculate()