# BarCrawler
An easy to use QC package for 10X genomics barcoded reads.

#Install

Run the install script, it will download the dependencies via pip and compile the python code

#operation

Bacrawler analyses 10X genomics bam files, and collect metrics such as molecule length, number of barcodes, coverage of the haplotypes etc. Barcrawler generates a histograms and plots of the coverage, as well as text files.

#Running

python BarCrawl.py --bam <input.bam> --working_dir <output_directory> 

The output files are stored in the working_dir directory. The filenames of the files will be set to the filename  of the bam file.
