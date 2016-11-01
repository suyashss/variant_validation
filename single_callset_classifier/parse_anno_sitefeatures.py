#!/usr/bin/python

import sys

def main():

    features=["AC","Dels","FS","HaplotypeScore","MQ","QD","BaseQRankSum","MQRankSum","ReadPosRankSum","GC","DP"]
    header=["CHR","START","END","QUAL"]
    header.extend(features)
    header[-1]="Coverage"
    print "\t".join(header)

    with open(sys.argv[1]) as f:
        for line in f.readlines():
            if line[0]=="#":
                continue
            fields=line.strip().split()
            outarr=[fields[0],str(int(fields[1])-1),fields[1]]
            if len(fields)==2:
                outarr.append("0")
                for f in features:
                    outarr.append("0")
            else:
                if fields[5]==".":
                    outarr.append("0") #Quality score has .
                    for f in features:
                        outarr.append("0")
                else:
                    outarr.append(fields[5])
                    info=fields[7]
                    iterinfo=iter(info.replace("=",";").split(";"))
                    d=dict(zip(iterinfo,iterinfo))
                    for f in features:
                        if f in d:
                            outarr.append(d[f])
                        else:
                            outarr.append("0")
            print "\t".join(outarr)


if __name__ == "__main__":
    main()
