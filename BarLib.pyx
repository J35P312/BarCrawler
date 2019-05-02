import sys
import os
import numpy
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.signal

def read_cigar(cigar):
    current_str=[]
    length=0
    for i in range(0,len(cigar)):
        if cigar[i].isdigit():
            current_str.append(cigar[i])
        elif cigar[i] == "M":
            length += int("".join(current_str))
            current_str=[]
        else:
            current_str=[]        
    return length

def hist(data,main_title,x_title,y_title,bins,filename,jumps,logX):

    weights= numpy.ones_like(numpy.array(data))/float(len(data))
    plt.hist(data, bins=bins, facecolor='grey', alpha=1,histtype="bar",weights=weights,edgecolor="black")
    plt.xticks(jumps)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.title(main_title)
    plt.grid(False)
    if logX:
        plt.xscale('log')
    figure = plt.gcf()
    figure.set_size_inches(16, 10)
    plt.savefig(filename)
    plt.close()


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """

    arr = numpy.ma.array(arr).compressed()
    med = numpy.median(arr)
    return numpy.median(numpy.abs(arr - med))


def coverage_summary(args,f):

    coverage_tot=[]
    coverage_h1=[]
    coverage_h2=[]
    coverage_null=[]
    no_hap=0
    hom=0
    for line in open(args.prefix+"_Coverage.tab"):
        if line[0] == "#":
            continue
        content=line.strip().split()
        tot=float(content[-4])
        h1=float(content[-3])
        h2=float(content[-2])
        null=tot-h1-h2
        
        coverage_tot.append(tot)
        coverage_h1.append(h1)
        coverage_h2.append(h2)
        coverage_null.append(null)

        if tot and not h1 and not h2:
            no_hap +=1
        if tot and h1 > 0.8*tot or h2 > 0.8*tot:
            hom+=1

    coverage_tot=numpy.array(coverage_tot)
    coverage_h1=numpy.array(coverage_h1)
    coverage_h2=numpy.array(coverage_h2)
    coverage_null=numpy.array(coverage_null)
    f.write("#Coverage summary\n" )
    f.write("average coverage:{}\n".format(numpy.average(coverage_tot)) )
    f.write("coverage over 10X:{}\n".format(len(coverage_tot[coverage_tot > 10])/float(len(coverage_tot))) )
    f.write("coverage over 20X:{}\n".format(len(coverage_tot[coverage_tot > 20])/float(len(coverage_tot))) )
    f.write("coverage over 30X:{}\n".format(len(coverage_tot[coverage_tot > 30])/float(len(coverage_tot))) )
    f.write("average depth h1:{}\n".format(numpy.average(coverage_h1)))
    f.write("average depth h2:{}\n".format(numpy.average(coverage_h2)))
    f.write("average depth non-haplotyped:{}\n".format(numpy.average(coverage_null)))
    f.write("fraction non-haplotyped:{}\n".format(numpy.average(no_hap/float(len(coverage_tot)))))
    f.write("fraction homozygous:{}\n".format(numpy.average(hom/float(len(coverage_tot)))))

    f.write("Median coverage:{}\n".format(numpy.median(coverage_tot)))
    f.write("MAD coverage:{}\n".format(mad(coverage_tot)))
    f.write("Median coverage h1:{}\n".format(numpy.median(coverage_h1)))
    f.write("MAD coverage h1:{}\n".format(mad(coverage_h1)))
    f.write("Median coverage h2:{}\n".format(numpy.median(coverage_h2)))
    f.write("MAD coverage h2:{}\n".format(mad(coverage_h2)))
    f.write("median non-haplotyped coverage:{}\n".format(numpy.median(coverage_null)))
    f.write("MAD no-haplotyped coverage:{}\n".format(mad(coverage_null)))

def generate_coverage_plot(args):
    coverage_structure={ "total_coverage":{},"1":{},"2":{}}

    with os.popen("samtools view -H {}".format(args.bam)) as pipe:
        for line in pipe:
            if line[0] == "@":
                if "SN:" in line:
                    content=line.strip().split()
                    chromosome=content[1].split("SN:")[-1]
                    length=int(content[2].split("LN:")[-1])
                    bins=int( math.ceil(length/float(args.bin_size)) )
                    coverage_structure["total_coverage"][chromosome]=numpy.zeros(bins)
                    coverage_structure["1"][chromosome]=numpy.zeros(bins)
                    coverage_structure["2"][chromosome]=numpy.zeros(bins)

    i=0
    chromosome=False
    os.system("mkdir {}/coverage_plots/".format(args.working_dir))
    for line in open(args.prefix+"_Coverage.tab"):
        if line[0] == "#":
            continue

        content=line.strip().split()
        if not chromosome or content[0] != chromosome:
            print chromosome
            chromosome=content[0]
            i=0
        coverage_structure["total_coverage"][content[0]][i]=float(content[-4])
        coverage_structure["1"][content[0]][i]=float(content[-3])
        coverage_structure["2"][content[0]][i]=float(content[-2])
        i+=1

    coverage_values=[]
    for chromosome in coverage_structure["total_coverage"]:
        for i in range(0,len( coverage_structure["total_coverage"][chromosome] )):
            if coverage_structure["total_coverage"][chromosome][i]:
               coverage_values.append(coverage_structure["total_coverage"][chromosome][i])


    median_coverage=numpy.median(coverage_values)

    for chromosome in coverage_structure["total_coverage"]:
        if "chrUn" in chromosome or "gl0" in chromosome or "KI" in chromosome:
            continue
        total_cov = scipy.signal.medfilt(coverage_structure["total_coverage"][chromosome],args.medfilt)
        h1 = scipy.signal.medfilt(coverage_structure["1"][chromosome],args.medfilt)
        h2 = scipy.signal.medfilt(coverage_structure["2"][chromosome],args.medfilt)

        median_coverage=numpy.median(coverage_structure["total_coverage"][chromosome])
        posvector=numpy.array( range(0,len(coverage_structure["1"][chromosome])))*args.bin_size/1000

        no_hap=total_cov-h1-h2

        total_cov[total_cov > 3*median_coverage] = median_coverage*3
        h1[h1 > 4*median_coverage] = 3*median_coverage
        h2[h2 > 4*median_coverage] =3*median_coverage
        no_hap[no_hap > 4*median_coverage] = 3*median_coverage


        median_no_hap=numpy.median(no_hap)

        plt.figure(1)
        plt.subplot(211)

        total=plt.scatter(posvector, total_cov, c="green", alpha=0.5,marker = 'o',label='Total coverage')
        non_haplotyped=plt.scatter(posvector, no_hap, c="grey", alpha=0.25,marker = 'o',label='non-haplotyped coverage')

        total_line, =plt.plot([0,max(posvector)], [median_coverage,median_coverage],c="black")
        non_haplotyped_line, =plt.plot([0,max(posvector)], [median_no_hap,median_no_hap],c="black",linestyle="--")

        plt.legend([total,non_haplotyped,total_line,non_haplotyped_line], ["Total", "Non-haplotyped","Median Total", "Median Non-haplotyped"])

        plt.ylim(ymax = 4*median_coverage, ymin = 0)
        plt.title('Coverage on chromosome {}'.format(chromosome))

       
        plt.xlabel('Positions(Kb)')
        plt.ylabel('Coverage')

        plt.subplot(212)
        hap1 =plt.scatter(posvector, h1, c="blue", alpha=0.5,marker = 'o')
        hap2 =plt.scatter(posvector, h2, c="red", alpha=0.25,marker = 'o')

        total_line, =plt.plot([0,max(posvector)], [median_coverage,median_coverage],c="black")
        non_haplotyped_line, =plt.plot([0,max(posvector)], [median_no_hap,median_no_hap],c="black",linestyle="--")

        plt.legend([hap1,hap2,total_line,non_haplotyped_line], ["Haplotype 1", "Haplotype 2","Median Total", "Median Non-haplotyped"])
        plt.title('Haplotyped coverage on chromosome {}'.format(chromosome))
        plt.ylim(ymax = 3*median_coverage, ymin = 0)
        plt.xlabel('Positions(Kb)')
        plt.ylabel('Coverage')
        figure = plt.gcf()
        figure.set_size_inches(16, 10)
        plt.savefig("{}/coverage_plots/{}_{}.png".format(args.working_dir,args.prefix.split("/")[-1],chromosome),dpi=100)
        plt.close()




def molecule_summary(args,f):
    reads=[]
    length=[]
    number_of_molecules=0

    for line in open(args.prefix+"_MoleculeReadHist.tab"):
        if line[0] == "#":
            continue
        content=line.strip().split()
        number_of_molecules+=int(content[0])*int(content[1])
        numbers=numpy.random.uniform(0, 1, size=int(content[0]))
        number_to_sample=len(numbers[numbers <= args.sample])
        reads+=number_to_sample*[int(content[1])]

    for line in open(args.prefix+"_MoleculeLenHist.tab"):  
        if line[0] == "#":
            continue
        content=line.strip().split()
        numbers=numpy.random.uniform(0, 1, size=int(content[1]))
        number_to_sample=len(numbers[numbers <= args.sample])
        length+=number_to_sample*[int(content[0])]            

    reads=numpy.array(reads)
    length=numpy.array(length)
    f.write("#Molecule summary\n" )
    f.write("Median reads per molecule:{}\n".format(numpy.median(reads)))
    f.write("MAD reads per molecule:{}\n".format(mad(reads)))
    f.write("Median molecule length:{}\n".format(numpy.median(length)))
    f.write("MAD molecule length:{}\n".format(mad(length)))
    f.write("Median molecule density(reads/kb):{}\n".format(1000*numpy.median(reads)/numpy.median(length)))
    f.write("Number of molecules:{}\n".format(int(number_of_molecules*args.sample)))
    f.write("Fraction of molecules longer than 10Kb:{}\n".format(len(length[length > 10000])/(number_of_molecules*args.sample) ))
    f.write("Fraction of molecules longer than 20KB:{}\n".format(len(length[length > 20000])/(number_of_molecules*args.sample) ))
    f.write("Fraction of molecules longer than 50KB:{}\n".format(len(length[length > 50000])/(number_of_molecules*args.sample) ))
    f.write("Fraction of molecules longer than 100Kb:{}\n".format(len(length[length > 100000])/(number_of_molecules*args.sample) ))

    hist(length/1000,"Molecule length distribution",'Molecule Length(Kb)','Count',numpy.arange(min(length/1000), 201, 5),args.prefix+"_MoleculeLenHist.png",numpy.arange(min(length/1000), 201, 5),False)
    hist(reads,"Molecule linked read distribution",'Linked reads','Count',numpy.arange(2, 101, 2),args.prefix+"_MoleculeReadHist.png",[2,10,20,40,100],False)

def generate_summary(args):
    f = open( args.prefix+"_Summary.txt", 'w')
    generate_coverage_plot(args)
    coverage_summary(args,f)
    molecule_summary(args,f)

    linked_reads=[]
    barcodes=[]
    number_of_reads_per_barcode=[]
    number_of_barcodes_per_read=[]
    not_barcoded=0
    number_of_reads=0
    number_of_barcodes=0

    for line in open(args.prefix+"_BarcodeHist.tab"):  
        if line[0] == "#":
            continue
        content=line.strip().split()

        if content[0] == "0":
            not_barcoded=int(content[1])
            number_of_reads+=int(content[1])    
        else:
            number_of_reads+=int(content[0])*int(content[1])
            number_of_barcodes+=int(content[1])

            numbers=numpy.random.uniform(0, 1, size=int(content[1]))

            number_to_sample=len(numbers[numbers <= args.sample])
            number_of_reads_per_barcode+= number_to_sample*[int(content[0])]

            numbers=numpy.random.uniform(0, 1, size=int(content[0]))

            number_to_sample=len(numbers[numbers <= args.sample])
            number_of_barcodes_per_read+= number_to_sample*[int(content[1])]

    number_of_reads_per_barcode=numpy.array(number_of_reads_per_barcode)
    number_of_barcodes_per_read=numpy.array(number_of_barcodes_per_read)

    f.write("#barcode summary\n" )
    f.write("Median reads per barcode:{}\n".format(numpy.median(number_of_reads_per_barcode)))
    f.write("MAD reads per barcode:{}\n".format(mad(number_of_reads_per_barcode)))
    f.write("Number of reads:{}\n".format(number_of_reads))
    f.write("Number of barcodes:{}\n".format(number_of_barcodes))
    f.write("Fraction of reads with no barcode:{}\n".format(not_barcoded/float(number_of_reads)))
    f.write("Fraction of barcodes with more than 4 linked reads:{}\n".format(len( number_of_barcodes_per_read[number_of_barcodes_per_read > 4])/float(len(number_of_barcodes_per_read) )))
    f.write("Fraction of barcodes with more than 10 linked reads:{}\n".format(len( number_of_barcodes_per_read[number_of_barcodes_per_read > 10])/float(len(number_of_barcodes_per_read))))    
    f.write("Fraction of barcodes with more than 20 linked reads:{}\n".format(len( number_of_barcodes_per_read[number_of_barcodes_per_read > 20])/float(len(number_of_barcodes_per_read))))
    f.write("Fraction of barcodes with more than 40 linked reads:{}\n".format(len( number_of_barcodes_per_read[number_of_barcodes_per_read > 40])/float(len(number_of_barcodes_per_read))))

    del number_of_reads_per_barcode
    del number_of_barcodes_per_read

    barcodes={}
    barcode_contigs={}

    for line in open(args.prefix+"_MoleculeSpan.tab"):  
        if line[0] == "#":
            continue
        content=line.strip().split()
        if "chrUn" in content[0] or "gl0" in content[0]:
            continue

        if not content[4] in barcodes:
            barcodes[content[4]] =0
            barcode_contigs[content[4]]=set([])
        barcodes[content[4]]+=1                
        barcode_contigs[content[4]].add(content[0])    

    bx_hist=[]
    for bx in barcodes:
        bx_hist.append(barcodes[bx])

    bx_hist=numpy.array(bx_hist)
    f.write("fraction of non-fragmented barcoddes:{}\n".format(len(bx_hist[bx_hist == 1])/float(len(bx_hist))))
    f.write("Median molecules per barcode:{}\n".format(numpy.median(bx_hist)))
    f.write("MAD molecules per barcode:{}\n".format(mad(bx_hist)))
    f.write("Fraction of barcodes split in more than 2 molecules:{}\n".format(len(bx_hist[bx_hist > 1])/float(len(bx_hist)) ))
    f.write("Fraction of barcodes split in more than 5 molecules:{}\n".format(len(bx_hist[bx_hist > 5])/float(len(bx_hist)) ))
    f.write("Fraction of barcodes split in more than 10 molecules:{}\n".format(len(bx_hist[bx_hist > 10])/float(len(bx_hist)) ))
    f.write("Fraction of barcodes split in more than 20 molecules:{}\n".format(len(bx_hist[bx_hist > 20])/float(len(bx_hist)) ))
    f.write("Fraction of barcodes split in more than 40 molecules:{}\n".format(len(bx_hist[bx_hist > 40])/float(len(bx_hist)) ))



    hist(bx_hist,"Molecules per barcode distribution",'Molecules per barcode','Count',numpy.arange(1, 51, 1),args.prefix+"_MoleculesPerBarcode.png",numpy.arange(1, 51, 1),False)

    chr_list=[]
    for barcode in barcode_contigs:
        chr_list.append( len(barcode_contigs[barcode]))

    chr_list=numpy.array(chr_list)
    f.write("Fraction of barcodes mapped to only 1 chromosome:{}\n".format(len(chr_list[chr_list == 1])/float(len(chr_list)) ))
    f.write("Fraction of barcodes mapped to more than 1 chromosome:{}\n".format(len(chr_list[chr_list > 1])/float(len(chr_list)) ))
    f.write("Fraction of barcodes mapped to more than 5 chromosomes:{}\n".format(len(chr_list[chr_list > 5])/float(len(chr_list)) ))

    hist(chr_list,"Chromosomes per barcode distibution",'Chromosome per barcode','Count',numpy.arange(1, 25, 1), args.prefix+"_ChromosomesPerBarcode.png",numpy.arange(1, 25, 5),False)

def compute_coverage(args):
    chromosome=-1
    barcode_dictionary={}
    barcode_dictionary["none"]=0
    coverage_structure={ "total_coverage":{},"1":{},"2":{},"Q":{} }
    global_molecule_dictionary={}
    length_dictionary={}
    read_dictionary={}
    f = open(args.prefix+"_MoleculeSpan.tab", 'w')
    f.write("#chromosome\tstart\tend\tMI\tBX\tHP\tn_reads\n")

    chromosome_order=[]
    with os.popen("samtools view -h -F 2048 -F 256 {}".format(args.bam)) as pipe:
        for line in pipe:
            if line[0] == "@":
                if "SN:" in line:
                    content=line.strip().split()
                    chromosome=content[1].split("SN:")[-1]
                    length=int(content[2].split("LN:")[-1])
                    bins=int( math.ceil(length/float(args.bin_size)) )
                    coverage_structure["total_coverage"][chromosome]=numpy.zeros(bins)
                    coverage_structure["1"][chromosome]=numpy.zeros(bins)
                    coverage_structure["2"][chromosome]=numpy.zeros(bins)
                    coverage_structure["Q"][chromosome]=numpy.zeros((bins,2))
                    chromosome_order.append(chromosome)
                continue
    
            content=line.strip().split()
            position=int(content[3])
            quality=int(content[4])
            read_length=read_cigar(content[5])
    
            if chromosome != content[2]:
                for MI in global_molecule_dictionary:
                    start=global_molecule_dictionary[MI][0][0]
                    end=global_molecule_dictionary[MI][0][1]
                    BX=global_molecule_dictionary[MI][1]
                    HP=global_molecule_dictionary[MI][0][2]
                    reads=global_molecule_dictionary[MI][0][3]
                    if not reads in read_dictionary:
                        read_dictionary[reads]=1
                    else:
                        read_dictionary[reads]+=1
    
                    length=(end-start)/500
                    length=int(round(length)*500)
                    if not length in length_dictionary:
                        length_dictionary[length]=1
                    else:
                        length_dictionary[length]+=1
    
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chromosome,start,end,MI,BX,HP,reads))
                global_molecule_dictionary={}
            chromosome=content[2]
            #read_length=150

            HP=0
            if "HP:i:" in line:
               HP=int(line.strip().split("HP:i:")[-1].split()[0])
    
            if "MI:" in line:
               MI=int(line.strip().split("MI:i:")[-1].split("\t")[0])
               if not MI in global_molecule_dictionary:
                   global_molecule_dictionary[MI]=[numpy.array([position,position,HP,1]),line.strip().split("BX:Z:")[-1].split("\t")[0]]
               else:
                   global_molecule_dictionary[MI][0][1]=position
                   global_molecule_dictionary[MI][0][-1]+=1
    

            if chromosome in chromosome_order:
               pos=int(math.floor(position/float(args.bin_size)))
               to_add=int(math.floor(read_length/float(args.bin_size)))
               for i in range(0,to_add):
                   coverage_structure["total_coverage"][chromosome][pos]+= 1
                   coverage_structure["Q"][chromosome][pos][0] += quality
                   coverage_structure["Q"][chromosome][pos][1] += 1
    
                   if "HP:i:" in line:
                       haplotype=line.strip().split("HP:i:")[-1].split()[0]
                       coverage_structure[haplotype][chromosome][pos] += 1
            
                   pos+=1

               rest=read_length-to_add
               if rest > 0:
                   coverage_structure["total_coverage"][chromosome][pos]+= rest/float(args.bin_size)
                   coverage_structure["Q"][chromosome][pos][0] += quality
                   coverage_structure["Q"][chromosome][pos][1] += 1

                   if HP:
                       coverage_structure[str(HP)][chromosome][pos] += rest/float(args.bin_size)

            if "BX:Z" in line:
               bx=line.strip().split("BX:Z:")[-1].split("\t")[0]
               if not bx  in barcode_dictionary:
                   barcode_dictionary[bx]=0
               barcode_dictionary[bx]+=1
    
            else:
                barcode_dictionary["none"]+=1
         

    for MI in global_molecule_dictionary:
        start=global_molecule_dictionary[MI][0][0]
        end=global_molecule_dictionary[MI][0][1]
        BX=global_molecule_dictionary[MI][1]
        HP=global_molecule_dictionary[MI][0][2]
        reads=global_molecule_dictionary[MI][0][3]
        if not reads in read_dictionary:
            read_dictionary[reads]=1
        else:
            read_dictionary[reads]+=1

        length=(end-start)/500
        length=int(round(length)*500)
        if not length in length_dictionary:
            length_dictionary[length]=1
        else:
            length_dictionary[length]+=1

        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chromosome,start,end,MI,BX,HP,reads))

    global_molecule_dictionary={}
    f.close()        
    f = open( args.prefix+"_Coverage.tab", 'w')
    f.write("#chromosome\tstart\tend\ttotal_coverage\tHL1_coverage\tHL2_coverage\tquality\n")
    for chromosome in chromosome_order:
        print chromosome
        for i in range(0,len(coverage_structure["total_coverage"][chromosome])):
            total=coverage_structure["total_coverage"][chromosome][i]
            hl1=coverage_structure["1"][chromosome][i]
            hl2=coverage_structure["2"][chromosome][i]
            if coverage_structure["Q"][chromosome][i][1]:
                qual=coverage_structure["Q"][chromosome][i][0]/coverage_structure["Q"][chromosome][i][1]
            else:
                qual=0
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chromosome,i*args.bin_size,(i+1)*args.bin_size,total,hl1,hl2,qual))

    
    f.close()
    barcode_hist={}
    barcode_hist[0]=barcode_dictionary["none"]
    for bx in barcode_dictionary:
        if not barcode_dictionary[bx] in barcode_hist:
            barcode_hist[barcode_dictionary[bx]]=0
        barcode_hist[barcode_dictionary[bx]]+=1
    f = open(args.prefix+"_BarcodeHist.tab", 'w')
    f.write("#linked_reads\tn_barcodes\n")
    for bx in sorted(barcode_hist):
        f.write("{}\t{}\n".format(bx,int(math.ceil(barcode_hist[bx]/2.0))))
    f.close()


    f = open(args.prefix+"_MoleculeReadHist.tab", 'w')  
    f.write("#linked_reads\tn_molecules\n")
    for mi in sorted(read_dictionary):
        f.write("{}\t{}\n".format(mi,read_dictionary[mi]))
    f.close()

    f = open( args.prefix+"_MoleculeLenHist.tab", 'w')
    f.write("#length\tmolecules\n")
    for bx in sorted(length_dictionary):
        f.write("{}\t{}\n".format(bx,length_dictionary[bx]))
    f.close()

    
