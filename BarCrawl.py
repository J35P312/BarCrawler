import BarLib
import argparse
import os
import gzip

parser = argparse.ArgumentParser("""BarCrawler - analysis package for 10X genomics linked reads data""",add_help=False)
parser.add_argument('--coverage',action="store_true",help="Analyse the coverage of a 10X genomics bam file")
parser.add_argument('--pssearch',action="store_true",help="test if the ps of a sample is found in another sample")
args, unknown = parser.parse_known_args()

if args.coverage:

    parser = argparse.ArgumentParser("""BarCrawler - analysis package for 10X genomics linked reads data""")
    parser.add_argument('--coverage',action="store_true",help="Analyse the coverage of a 10X genomics bam file")
    parser.add_argument('--bin_size',default=100,type=int,help="size of the reported bins")
    parser.add_argument('--medfilt',default=51,type=int,help="size of the median filter for the coverage plots")
    parser.add_argument('--bam',required=True,type=str,help="The bam file")
    parser.add_argument('--working_dir',default="",required=True,type=str,help="the files are put int his folder")
    parser.add_argument('--sample',default=0.05,type=float,help="compute barcode statistics based on a subset of barcodes (default=0.05)")
    args = parser.parse_args()

    os.system("mkdir {}".format(args.working_dir))
    args.prefix=os.path.join(args.working_dir,args.bam.split("/")[-1].replace(".bam",""))
    BarLib.compute_coverage(args)
    BarLib.generate_summary(args)

elif args.pssearch:
    parser = argparse.ArgumentParser("""BarCrawler - analysis package for 10X genomics linked reads data""")
    parser.add_argument('--pssearch',action="store_true",help="test if the ps of a sample(A) is found in another sample (B)")
    parser.add_argument('--vcfa',required=True,type=str,help="the (gzipped) vcf of sample A. This sample should contain the PS of interest")
    parser.add_argument('--vcfb',required=True,type=str,help="The (gzipped) vcf of sample B. This sample will be searched for variants similar to the ps of sample A")
    parser.add_argument('--ps',required=True,type=str,help="the phase set (PS) of interest")
    args = parser.parse_args()

    snvsa=0
    snv_structure={}


    chromps=0
    startps=0
    endps=0

    for line in gzip.open(args.vcfa):
        if line[0] == "#":
            continue
    
        format=line.strip().split()[8].split(":")
        ps=-1
        i=0
        for col in format:
            if col == "PS":
                ps=i
            i+=1

        if ps > -1 and not "\tPHASING_INC" in line and "HAPLOCALLED=1" in line:
            if not line.strip().split()[9].split(":")[ps] == args.ps:
                continue

            content=line.strip().split()
            chrom=content[0]
            pos=int(content[1])
            change=content[4]


            if not chrom in snv_structure:
                snv_structure[chrom]={}
                chromps=chrom
                startps=pos

            if not pos in snv_structure[chrom]:
                snv_structure[chrom][pos]=[]
            snv_structure[chrom][pos].append(change)
            endps=pos
            snvs+=1


    psb={}
    snvsb=0
    for line in gzip.open(args.vcfb):

        #skip header
        if line[0] == "#" or "HAPLOCALLED=0" in line:
            continue
    
        #skip inconsistent phasing
        if "\tPHASING_INC" in line:
            continue

        content=line.strip().split()
        #variant is on wrong chromosome
        if not content[0] in snv_structure:
            continue

        #this position is not in the phase set
        if not int(content[1]) in snv_structure[content[0]]:
            continue

        #check if this DNA change is found on the current position
        if content[4] in snv_structure[content[0]][int(content[1])]:
            ps=-1
            i=0
            for col in format:
                if col == "PS":
                    ps=i
                i+=1
            if ps  == -1:
                continue
                        
            ps=content[9].split(":")[ps]
            if not ps in psb:
                psb[ps]={}
                psb[ps]["chrom"]=content[0]
                psb[ps]["start"]=int(content[1])
                psb[ps]["end"]=int(content[1])
                psb[ps]["snvs"]=1
            else:
                psb[ps]["end"]=int(content[1])
                psb[ps]["snvs"]+=1

            snvsb+=1



    print "span of phase set {}: {}:{}-{}".format(chromps,startps,endps)
    print "SNVs on phase set {}:     {}".format(args.ps,snvsa)
    print "SNVs matched to sample b: {}".format(snvsb)
    print "percentage matched to sample b: {}".format(100*snvsa/float(snvsb))
    for ps in psb:
        print "span of ps {}: span: {}:{}-{} snvs : {} percentage: {}".format(ps,psb[ps]["chrom"],psb[ps]["start"],psb[ps]["end"],psb[ps]["snvs"],psb[ps]["snvs"]/float(snvsb)*100)

else:
    print "hi there, and welcome!"
    print "To produce a help message, type:"
    print "python BarCrawl.py --help"

