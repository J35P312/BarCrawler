# BarCrawler
An easy to use QC package for 10X genomics barcoded reads.

#Install

Run the install script, it will download the dependencies via pip and compile the python code

#operation

Barcrawler contains two separate modules, the coverage module, and the pssearch module.

The coverage module analyses 10X genomics bam files, and collect metrics such as molecule length, number of barcodes, coverage of the haplotypes etc. Barcrawler generates a histograms and plots of the coverage, as well as text files.

The pssearch module accepts two 10x barrcoded SNV vcf files, as well as a ps tag of a phase set of snvs in one of those samples. pssearch will then check if the same ps exists in the other sample.
 
#Running

The coverage module is run through the following command:

    python BarCrawl.py --coverage --bam <input.bam> --working_dir <output_directory> 

The output files are stored in the working_dir directory. The filenames of the files will be set to the filename  of the bam file.

The pssearch module is run through the following command:

        python BarCrawl.py --pssearch --vcfa sampleA.vcf --vcfb sampleB.vcf --ps PSID > output.tab


The PS is expected to be in sampleA, and sampleB will be searched for a similar haplotype. The results are printed to stdout.
