import sys,os

afasta=sys.argv[2]
gff=sys.argv[1]
gffcomp=gff+".complement"

os.system("/usr/global/blp/bin/samtools faidx "+afasta);
os.system("/usr/global/blp/bin/complementBed -i "+gff+" -g "+afasta+".fai > "+gffcomp);

gffCompFile=open(gffcomp)
gffclines=gffCompFile.readlines()
gffCompFile.close()

seqname=""
igDistances=[]
for l in gffclines:
    cols=l.split("\t")
    if cols[0]!=seqname:
        seqname=cols[0]
        continue
    else:
        igDistances.append(int(cols[2])-int(cols[1]))

avgigdist=sum(igDistances)*1.0/len(igDistances)
print avgigdist