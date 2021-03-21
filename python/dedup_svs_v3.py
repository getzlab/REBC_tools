import sys
import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
            required=True, help="The input SV file.", dest="input")
    parser.add_argument("-s", "--slop",
            type=int, required=False,
            default=200, help="The amount of slop to add to breakpoints for overlap checking", dest="slop")
    parser.add_argument("-V", "--invert",
            action="store_true",
            help="Return events that would normally be filtered, rather than passing events.", required=False)
    parser.add_argument("-N", "--max-n-alt", type=int, required=False,
            default=0,dest="max_alt",
            help="The maximum number of alt reads in the normal for a variant to pass.")
    parser.add_argument("-A", "--min-n-alg", type=int, required=False,
            default=3,dest="min_alg",
            help="The minimum number of algorithms to override max-n-alt criteria")
    parser.add_argument("-v", "--min-vaf", type=float, required=False,
            default=0.09,dest="min_vaf",
            help="minimum VAF")
    parser.add_argument("-T", "--min-t-alt", type=int, required=False,
            default=3,dest="min_talt",
            help="The minimum number of alt counts in tumor")

    return parser.parse_args()

def pos_diff(pos_first, pos_second):
    return (abs(pos_first - pos_second))

def lines_duplicated(first, seconds, header_d, bp_tolerance=200):

    for second in seconds:
        tokens_first = first.strip().split("\t")
        tokens_second = second.strip().split("\t")
    
        ## Check the sample ID, the chroms, the class, and finally whether
        ## pos1 / pos2 are within bp_tolerance basepairs of each other
        chroms_first = sorted([tokens_first[header_d["chr1"]], tokens_first[header_d["chr2"]]])
        chroms_second = sorted([tokens_second[header_d["chr1"]], tokens_second[header_d["chr2"]]])

        pos_first = sorted([int(tokens_first[header_d["pos1"]]), int(tokens_first[header_d["pos2"]])])
        pos_second = sorted([int(tokens_second[header_d["pos1"]]), int(tokens_second[header_d["pos2"]])])
        

        strands_first = "".join([tokens_first[header_d["str1"]], tokens_first[header_d["str2"]]])
        strands_second = "".join([tokens_second[header_d["str1"]], tokens_second[header_d["str2"]]])

        class_first = tokens_first[header_d["class"]]
        class_second = tokens_second[header_d["class"]]

        strands_same = strands_first == strands_second

        chroms_same = chroms_first[0] == chroms_second[0] and chroms_first[1] == chroms_second[1]

        pos_same = pos_diff(pos_first[0], pos_second[0]) < bp_tolerance and \
                pos_diff(pos_first[1], pos_second[1]) < bp_tolerance

        class_same = class_first == class_second

        if chroms_same and pos_same and strands_same and class_same:
            return True
    
    return False


def gen_line_id(line, header_d):
    
    tokens = line.strip().split("\t")
    g_id = "_".join([tokens[header_d["individual"]], tokens[header_d["chr1"]], tokens[header_d["chr2"]], tokens[header_d["class"]]])
    return g_id, line, tokens

if __name__ == "__main__":
    
    header_d = defaultdict(int)

    ## A dictionary that holds a first-pass filter for each line
    ## This only works because SVs are generally sparse in our dataset.
    sv_d = defaultdict(list)

    lines = []

    args = parse_args()
    with open(args.input, "r") as ifi:
        for line in ifi:
            #line = line.strip()
            if not line.startswith("individual"):
                d_id, l, tokens = gen_line_id(line, header_d)
                pass_nalt = int(tokens[header_d["VCF_NALT"]]) <= args.max_alt 
                pass_nalg = int(tokens[header_d["NALG"]]) >= args.min_alg
                TALT=float(tokens[header_d["VCF_TALT"]])
                TREF=float(tokens[header_d["VCF_TREF"]])
                VAF=TALT/(TALT+TREF/2)
                pass_vaf = (VAF >= args.min_vaf)
                pass_talt = (TALT >= args.min_talt)
                pass1=pass_nalt
                if pass_nalg:
                    pass1=True
                if not pass_vaf:
                    pass1=False
                if not pass_talt:
                    pass1=False

                is_dup =  lines_duplicated(line, sv_d[d_id], header_d, args.slop)
                #print(is_dup, args.invert, pass_nalt)
                ## If a sample has no entries in the dictionary,
                ## the line being processed cannot be a duplicate
                if not d_id in sv_d:
                    lines.append(line)
                elif (not is_dup and not args.invert) and pass1:
                    lines.append(line)
                elif args.invert and (is_dup or not pass1):
                    lines.append(line)
                else:
                    #print('pass',int(tokens[header_d["VCF_NALT"]]),int(tokens[header_d["NALG"]]) )
                    pass
                        #print("Removed:", line, file=sys.stderr)
                sv_d[d_id].append(line)
            else:
                tokens = line.strip().split("\t")
                for i in range(0, len(tokens)):
                    header_d[tokens[i]] = i
    
    print("\t".join(header_d.keys()))
    for i in lines:
        print(i)
