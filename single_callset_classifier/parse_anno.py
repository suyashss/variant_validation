#!/usr/bin/python

#Example input to parse
#chr22   16052238        16052239        chr22   16052239        .       A       G  78.77    .       AC=1;AF=0.500;AN=2;BaseQRankSum=0.940;DP=45;Dels=0.00;FS=8.104;HaplotypeScore=0.9967;MLEAC=1;MLEAF=0.500;MQ=166.13;MQ0=0;MQRankSum=0.496;QD=1.75;ReadPosRankSum=-0.120  GT:AB:AD:DP:GQ:PL       0/1:0.870:39,6:45:99:107,0,1295 1


import sys

def main():
    with open(sys.argv[1]) as f:
        for line in f.readlines():
            if line[0]=="#":
                continue
            fields=line.strip().split()
            if len(fields)==4:
                outarr=[fields[0],str(int(fields[1])-1),fields[1],"1"]
            else:
                if fields[4]==".":
                    outarr=[fields[0],str(int(fields[1])-1),fields[1],"1"]
                else:
                    sample=fields[9].split(":")
                    format=fields[8].split(":")
                    try:
                        idx=format.index("AB")
                        ab=sample[idx]
                    except ValueError:
                        ab="1"
                    outarr=[fields[0],str(int(fields[1])-1),fields[1],str(ab)]
            print "\t".join(outarr)

if __name__ == "__main__":
    main()
