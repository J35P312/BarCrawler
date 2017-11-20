import BarLib
import argparse
import os

parser = argparse.ArgumentParser("""BarCrawl - analyse the coverage of haplotypes""")
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

