#!/usr/bin/python
# -*- coding: utf-8 -*-
##### TODO
'''
1) Separate per contig, do it by listing samples
2) Automatic process to identify sample name and do count per sample ([pops])
3) write output table
4) Add SNP filter to ignore MSP and indels
5) Figure out how monomorphic positions are terated
6) expected allele frequencies need to be 0 for invariant reference sites
'''
import sys

'''
This scripts take a single allelic frequency table and calculates pi and Tajima's D for sliding widows

'''

'''
This scripts take a single allelic frequency table and calculates pi and Tajima's D for sliding widows
Usage: PoolSeq_pi_TajimasD.py infile [pops] nPOOL. Writes to StdOut
infile = allelic frequency table (see example below)
[pops] = name of populations included in the infile. The name has to be an exact match
nPOOL = number of individuals in the pool in the same order as [pops]

---------   Allelic frequency table example -----
contig	pos	ref	alt	type	qual	INFOestAllfreq	GQ	RO	AO	TotCount	ReFreq	AltFreq	sample
Contig0	1	T	NA	.	0	NA	NA	NA	NA	0	NA	NA	P2202_101_1to4
Contig0	2	T	NA	.	0	NA	NA	NA	NA	0	NA	NA	P2202_101_1to4
Contig0	3	A	NA	.	0	NA	NA	NA	NA	0	NA	NA	P2202_101_1to4
Contig0	4	A	NA	.	0	NA	NA	NA	NA	0	NA	NA	P2202_101_1to4
....
Contig0	761	T	NA	.	160.002	NA	160.002	9	NA	9	1	NA	P2202_101_1to4
Contig0	762	C	T	snp	251.04	0.535	0.0434226	1	8	9	0.111111111111111	0.888888888888889	P2202_101_1to4
Contig0	763	A	NA	.	160.002	NA	160.002	9	NA	9	1	NA	P2202_101_1to4
Contig0	765	C	G	snp	356.169	0.625	0.103332	0	9	9	0	1	P2202_101_1to4
Contig0	766	T	NA	.	160.002	NA	160.002	8	NA	8	1	NA	P2202_101_1to4
Contig0	769	A	G	snp	82.0797	0.24	0.059195	5	4	9	0.555555555555556	0.444444444444444	P2202_101_1to4
Contig0	770	A	NA	.	160.002	NA	160.002	9	NA	9	1	NA	P2202_101_1to4
...
Contig10	39092	A	NA	.	153.351	NA	153.351	1	NA	1	1	NA	P2202_103_1to4
Contig10	39093	C	NA	.	160.002	NA	160.002	1	NA	1	1	NA	P2202_103_1to4
Contig10	39094	A	C	snp	24.4162	0.25	3.20017	1	0	1	1	0	P2202_103_1to4
Contig10	39095	C	NA	.	160.002	NA	160.002	1	NA	1	1	NA	P2202_103_1to4
Contig10	39096	C	A	snp	24.4145	0.25	3.20198	1	0	1	1	0	P2202_103_1to4
Contig10	39097	C	NA	.	160.002	NA	160.002	1	NA	1	1	NA	P2202_103_1to4
Contig10	39098	C	NA	.	160.002	NA	160.002	1	NA	1	1	NA	P2202_103_1to4
Contig10	39099	C	NA	.	160.002	NA	160.002	1	NA	1	1	NA	P2202_103_1to4
Contig10	39100	C	G	snp	24.4145	0.25	3.20198	1	0	1	1	0	P2202_103_1to4
Contig10	39101	C	NA	.	160.002	NA	160.002	1	NA	1	1	NA	P2202_103_1to4

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
#print freq_file
#Pi_output = open("Pi_windows.txt", "w")
#Pi_output.write()
sys.stdout.write("sample" + "\t" + "contig" + "\t" + "win.start" + "\t" + "win.end" + "\t" + "win.len" + "\t" + "freq" + "\t"+ "pi"+ "\t" + "TajD"+ "\t" + "missing_data" + "\t" + "win_depth" +"\n" )

pi_list = []
freq_list = []
depth_list=[]
positions = []
last_contig = []
pos = 0

for line in freq_file:
	#sample=line.split()[13] 
	#if  sample == sys.argv[2]:
	#	continue
	#else:
	#print(line)
#firstLine = freq_file.readline()
    position = line.split()[1]
    #print position
    #print("Processing position: " + str(position))

    try:
        contig = line.split()[0]
        last_field = line.split()[13]
        qual = float(line.split()[5]) # QUAL: prob of error = 1−10^−(Qual/10)
        Vartype =  line.split()[4]
        depth =  float(line.split()[10])
        sample_name = line.split()[13]
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
        #print("PASS QUAL: "+ str(qual) + " sample=" + str(sample_name) + " Contig= " + str(contig) + " Pos= " + str(position) + " Pi = " + str(site_pi))
        site_depth = float(depth)  #depth reference plus depth alt in IA and CP respectively.
        freq_list.append(freq_Alt) #enter the value into the list
        pi_list.append(site_pi)     #enter the value into the list, including zeros
        depth_list.append(depth)
        last_contig.append(str(contig))

        positions.append(position)
        
        if qual > 20 or 100 < depth > 4:
        	pos += 1
        	#print("Contigs: current = " + contig + " last = " + last_contig[0] + last_contig[-1] + " pos count: " + str(pos))
            
    if pos == 100 and str(last_contig[0]) == str(last_contig[-1]):     #Every 100 snps, calculate end position of window and avg. pi for window
        contig = last_contig[0]
        window_start = float(positions[0])
        window_end = float(positions[-1])
        window_lenght = float(window_end-window_start)
        window_freq = float(sum(freq_list))/len(freq_list)  #window avg alt allele freq equals avg of all incl. fixed sites
        window_pi = float(sum(pi_list))/len(pi_list)  #window pi is equal to average of site pi across all sites in window, including invariant ones
        #print("Lenght pi list " + str(len(pi_list)))
        #print last_contig
        #print("##############    Window: contig = " + str(contig) + " sample=" + str(sample_name) + " start=" + str(window_start) + " end=" + str(window_end) + " len=" + str(window_lenght) + " Pi = " + str(window_pi) + "    #####################")
        missing_data = window_lenght - len(pi_list) + 1
        window_tajD = 0
        window_depth = float(sum(depth_list))/len(depth_list)
        last_contig = contig
        
        sys.stdout.write(str(sample_name) + "\t" + str(contig) + "\t" + str(window_start) + "\t" + str(window_end) + "\t" + str(window_lenght) + "\t" + str(window_freq) + "\t" + str(window_pi) + "\t" + str(window_tajD) + "\t" + str(missing_data) + "\t" + str(window_depth) + "\n")

        pi_list = []
        freq_list = []
        depth_list=[]
        positions = []
        last_contig = []
        pos = 0
     
    elif str(last_contig[0]) != str(last_contig[-1]):     #Every 100 snps, calculate end position of window and avg. pi for window
        contig = last_contig[0]
        window_start = float(positions[0])
        window_end = float(positions[-2])
        window_lenght = float(window_end-window_start)
        window_freq = float(sum(freq_list[0:-1]))/len(freq_list[0:-1])  #window avg alt allele freq equals avg of all incl. fixed sites
        window_pi = float(sum(pi_list[0:-1]))/len(pi_list[0:-1])  #window pi is equal to average of site pi across all sites in window, including invariant ones
        #print("Lenght pi list " + str(len(pi_list)))
        #print last_contig
        #print("####################  LAST  Window: contig = " + str(contig) + " sample=" + str(sample_name) + " start=" + str(window_start) + " end=" + str(window_end) + " len=" + str(window_lenght) + " Pi = " + str(window_pi) + "    #####################")
        missing_data = window_lenght - len(pi_list[0:-1]) + 1
        window_tajD = 0
        window_depth = float(sum(depth_list[0:-1]))/len(depth_list[0:-1])
        last_contig = contig
        
        sys.stdout.write(str(sample_name) + "\t" + str(contig) + "\t" + str(window_start) + "\t" + str(window_end) + "\t" + str(window_lenght) + "\t" + str(window_freq) + "\t" + str(window_pi) + "\t" + str(window_tajD) + "\t" + str(missing_data) + "\t" + str(window_depth) + "\n")

        pi_list = [pi_list[-1]]
        freq_list = [freq_list[-1]]
        depth_list=[depth_list[-1]]
        positions = [positions[-1]]
        last_contig = []
        pos = 0
        
#get the last window
if pos > 1:
	contig = last_contig[0]
	window_start = float(positions[0])
	window_end = float(positions[-1])
	window_lenght = float(window_end-window_start)
	window_freq = float(sum(freq_list))/len(freq_list)  #window avg alt allele freq equals avg of all incl. fixed sites
	window_pi = float(sum(pi_list))/len(pi_list)  #window pi is equal to average of site pi across all sites in window, including invariant ones
	#print "Last window"
	#print("Lenght pi list " + str(len(pi_list)))
	#print("##############     Window: contig = " + str(contig) +" sample=" + str(sample_name) + " start=" + str(window_start) + " end=" + str(window_end) + " len=" + str(window_lenght) + " Pi = " + str(window_pi) + "   ############## ")
	window_tajD = 0
	missing_data = (window_end - window_start) - len(pi_list)
	window_depth = float(sum(depth_list))/len(depth_list)

    
	sys.stdout.write(str(sample_name) + "\t" +  str(contig) + "\t" + str(window_start) + "\t" + str(window_end) + "\t" + str(window_lenght) + "\t" + str(window_freq) + "\t" + str(window_pi) + "\t" + str(window_tajD) + "\t" + str(missing_data) + "\t" + str(window_depth) + "\n")

freq_file.close()
#print( "Completed pi windows for scaffold " + str(scaff) +"\t" )