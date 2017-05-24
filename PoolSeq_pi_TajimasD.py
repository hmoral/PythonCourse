#!/usr/bin/python
# -*- coding: utf-8 -*-
##### TODO
'''
1) Separate per contig
2) Automatic process to identify sample name and do count per sample ([pops])
3) write output table
4) Add SNP filter to ignore MSP and indels
5) Figure out how monomorphic positions are terated
'''
import sys

'''
This scripts take a single allelic frequency table and calculates pi and Tajima's D for sliding widows

'''

print '''
This scripts take a single allelic frequency table and calculates pi and Tajima's D for sliding widows
Usage: PoolSeq_pi_TajimasD.py infile [pops] nPOOL. Writes to StdOut
infile = allelic frequency table (see example below)
[pops] = name of populations included in the infile. The name has to be an exact match
nPOOL = number of individuals in the pool in the same order as [pops]

---------   Allelic frequency table example -----
contig	pos	ref	alt	type	qual	INFOestAllfreq	P2202_101_1to4_GQ	P2202_101_1to4_RO	P2202_101_1to4_AO	P2202_101_1to4_TotCount	P2202_101_1to4_ReFreq	P2202_101_1to4_AltFreq	P2202_103_1to4_GQ	P2202_103_1to4_RO	P2202_103_1to4_AO	P2202_103_1to4_TotCount	P2202_103_1to4_ReFreq	P2202_103_1to4_AltFreq
Contig0	1	T		.	0	NA	NA	NA	NA	0	NA	NA	160.002	1	NA	1	1	NA
Contig0	2	T		.	0	NA	NA	NA	NA	0	NA	NA	160.002	1	NA	1	1	NA
Contig0	3	A		.	0	NA	NA	NA	NA	0	NA	NA	160.002	1	NA	1	1	NA
....
Contig0	126	T		.	2.86661e-15	NA	151.804	1	NA	1	1	NA	151.804	6	NA	6	1	NA
Contig0	127	CTG	CG	del	79.264	0.75	0.0464733	0	1	1	0	1	0.0460044	3	3	6	0.5	0.5
Contig0	130	G		.	2.86661e-15	NA	151.804	1	NA	1	1	NA	151.804	6	NA	6	1	NA
Contig0	131	T		.	2.86661e-15	NA	151.804	1	NA	1	1	NA	151.804	6	NA	6	1	NA
Contig0	132	T		.	2.86661e-15	NA	151.804	1	NA	1	1	NA	151.804	6	NA	6	1	NA
Contig0	133	TT	TGT	ins	79.2638	0.75	0.0462747	0	1	1	0	1	0.0460065	3	3	6	0.5	0.5
Contig0	135	T		.	0	NA	160.002	1	NA	1	1	NA	160.002	7	NA	7	1	NA
...
Contig0	3964	A		.	1.26568e-14	NA	145.355	65	NA	65	1	NA	145.355	97	NA	97	1	NA
Contig0	3965	G	A	snp	2331.78	0.635	2.81866	20	35	55	0.363636363636364	0.636363636363636	3.64236	35	59	94	0.372340425531915	0.627659574468085
Contig0	3966	G		.	4.6712e-15	NA	149.684	55	NA	55	1	NA	149.684	90	NA	90	1	NA
Contig0	3967	A		.	0	NA	160.002	56	NA	56	1	NA	160.002	93	NA	93	1	NA
Contig0	3968	A		.	0	NA	160.002	55	NA	55	1	NA	160.002	91	NA	91	1	NA
Contig0	3969	C		.	0	NA	160.002	55	NA	55	1	NA	160.002	90	NA	90	1	NA
Contig0	3970	T	G	snp	411.033	0.19	2.82409	44	11	55	0.8	0.2	0.248166	74	16	90	0.822222222222222	0.177777777777778
Contig0	3971	C		.	0	NA	160.002	54	NA	54	1	NA	160.002	90	NA	90	1	NA
'''

### Arguments

#infile = sys.argv[0]
#pops = sys.argv[1]
#nPOOL = sys.argv[2]
#nPOOL = sys.argv[1]

#freq_file = open("DUMMYFreqTable.txt", "r")
#freq_file = open("OUT_multiThreat.ALL.GQ.FreqTable.txt", "r")
#freq_file = pandas.read_table("OUT_multiThreat.ALL.GQ.FreqTable.txt",nrows=10)
#freq_file = pandas.read_table(infile)


freq_file = open(sys.argv[1], "r")
freq_file.readline()
print freq_file
#Pi_output = open("Pi_windows.txt", "w")
#Pi_output.write()
#sys.stdout.write("contig" + "\t" + "win.start" + "\t" + "win.end" + "\t" + "win.len" + "\t" + "freq" + "\t"+ "pi"+ "\t" + "TajD"+ "\t" + "missing_data" + "\t" + "win_depth" +"\n" )

pi_list = []
tajD_list = []
freq_list = []
depth_list=[]
positions = []
pos = 0

for line in freq_file:
    print(line)
#firstLine = freq_file.readline()
    position = line.split()[1]
    #print position
    #print("Processing position: " + str(position))

    try:
        contig = line.split()[0]
        last_field = line.split()[17]
        qual = float(line.split()[5]) # QUAL: prob of error = 1−10^−(Qual/10)
        Vartype =  line.split()[4]
        depth =  float(line.split()[10])
    except ValueError: # ignore exception from having "NA" QUAL.
        continue
    if qual < 20 or 100 > depth < 4 : # if quality is below 20 or depth is below 4 or above 100 skip position
        continue
    else:
        try:
        	freq_Alt = float(line.split()[12])  #expected allele frequencies, will be 0 for invariant reference sites
    	except ValueError: # ignore exception from having "NA" QUAL.
        	continue

        site_pi = 2 * freq_Alt * (1-freq_Alt) * (float(sys.argv[2])/(float(sys.argv[2]) - 1))  #Tajima 1989
        print("PASS QUAL: "+ str(qual) + " Contig= " + str(contig) + " Pos= " + str(position) + " Pi = " + str(site_pi))
        site_depth = float(depth)  #depth reference plus depth alt in IA and CP respectively.
        freq_list.append(freq_Alt) #enter the value into the list
        pi_list.append(site_pi)     #enter the value into the list, including zeros
        depth_list.append(depth)

        positions.append(position)
        if qual > 20 or 100 < depth > 4:
        	pos += 1
        	#print("pos count: " + str(pos))
            
    if pos == 100:     #Every 100 snps, calculate end position of window and avg. pi for window
        contig = line.split()[0]
        window_start = float(positions[0])
        window_end = float(positions[-1])
        window_lenght = float(window_end-window_start)
        window_freq = float(sum(freq_list))/len(freq_list)  #window avg alt allele freq equals avg of all incl. fixed sites
        window_pi = float(sum(pi_list))/len(pi_list)  #window pi is equal to average of site pi across all sites in window, including invariant ones
        print("#####    Window: contig = " + str(contig) + " start=" + str(window_start) + " end=" + str(window_end) + " len=" + str(window_lenght) + " Pi = " + str(window_pi) + "    #####")
        missing_data = (window_end - window_start) - len(pi_list)
        window_tajD = 0
        window_depth = float(sum(depth_list))/len(depth_list)
        
        #sys.stdout.write(str(contig) + "\t" + str(window_start) + "\t" + str(window_end) + "\t" + str(window_lenght) + "\t" + str(window_freq) + "\t" + str(window_freq) + "\t" + str(window_pi) + "\t" + str(window_tajD) + "\t" + str(missing_data) + "\t" + str(window_depth))

        pi_list = []
        tajD_list = []
        freq_list = []
        depth_list=[]
        positions = []
        pos = 0
#get the last window
if pos > 1:
	contig = line.split()[0]
	window_start = float(positions[0])
	window_end = float(positions[-1])
	window_lenght = float(window_end-window_start)
	window_freq = float(sum(freq_list))/len(freq_list)  #window avg alt allele freq equals avg of all incl. fixed sites
	window_pi = float(sum(pi_list))/len(pi_list)  #window pi is equal to average of site pi across all sites in window, including invariant ones
	print "Last window"
	print("#####    Window: contig = " + str(contig) + " start=" + str(window_start) + " end=" + str(window_end) + " len=" + str(window_lenght) + " Pi = " + str(window_pi) + "    #####")
	window_tajD = 0
	missing_data = (window_end - window_start) - len(pi_list)
	window_depth = float(sum(depth_list))/len(depth_list)

    
    #sys.stdout.write(str(contig) + "\t" + str(window_start) + "\t" + str(window_end) + "\t" + str(window_lenght) + "\t" + str(window_freq) + "\t" + str(window_freq) + "\t" + str(window_pi) + "\t" + str(window_tajD) + "\t" + str(missing_data) + "\t" + str(window_depth))

freq_file.close()
#print( "Completed pi windows for scaffold " + str(scaff) +"\t" )