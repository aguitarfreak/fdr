'''
Created on May 29, 2011

@author: Ahwan
'''
from __future__ import division
import sys, csv
import numpy as np
import nipy.algorithms.statistics.empirical_pvalue as fdr
import optparse
import matplotlib.pylab as mp

#=====================================================================================
def switch_snp_key_order(snps):
    return (snps[0],snps[1]), (snps[1], snps[0])
    
def create_I3_dictionary(gain_file, threshold, verbose):
    snps =None
    I3_dict = {} #threshold filtered
    all_I3_dict = {} #all
    reader = csv.reader(gain_file,delimiter = '\t')
    snps = reader.next()
    for row in reader:
        for col in enumerate(row):
            i = col[0]
            j = reader.line_num -2
            val = float(col[1])
            if(i > j):
                all_I3_dict[snps[j],snps[i]] = val
                if(val >= threshold):
                    I3_dict[snps[j],snps[i]] = val
    
    return I3_dict, all_I3_dict

def create_epi_dictionary(epistasis_file, verbose):
    header =None
    epi_dict = {}
    reader = csv.reader(epistasis_file, delimiter = ' ')
    header = reader.next()
    #remove '' elements from list
    header[:] = (value for value in header if value != '')
    #get indices
    snp1_idx = header.index('SNP1')
    snp2_idx = header.index('SNP2')
    p_val_idx = header.index('P')
    temp_row=[]
    for row in reader:
        temp_row[:] = (value for value in row if value != '')
        epi_dict[temp_row[snp1_idx],temp_row[snp2_idx]] = float(temp_row[p_val_idx])
    return epi_dict

def create_fdr_dictionary(epi_dict,verbose):
    f = fdr.FDR()
    fdrs =  f.all_fdr(epi_dict.values(), 0)
    fdr_dict = dict(zip(epi_dict.keys(),fdrs))
    return fdr_dict

def filter_epi_dict(epi, max_p):
    q_val_snp = ''
    filtered_dict = {}
    #sort the dictionary low to high
    sorted_dictionary_list = sorted(epi.iteritems(), key=lambda (k,v): (v,k))
    for key,value in sorted_dictionary_list:
        if(value<max_p):
            filtered_dict[key] = value
        else:
            q_val_snp = key
            break
    return filtered_dict, q_val_snp

def create_filtered_I3_dictionary(p_val_filtered_epi,threshold_filtered_I3):
    #to do
    #max of two dictionaries
    fdr_filtered_I3_dict = {}
    for snps in threshold_filtered_I3:
        order = switch_snp_key_order(snps)
        if p_val_filtered_epi.has_key(order[0]) or p_val_filtered_epi.has_key(order[1]):
            fdr_filtered_I3_dict[snps] = threshold_filtered_I3[snps]
    return fdr_filtered_I3_dict

def fdr_plot(epi,fdrs):
    epi_x = []
    fdrs_y = []
    for snps in epi.keys():
        epi_x.append(epi[snps])
        fdrs_y.append(fdrs[snps])
    mp.figure(1)
    mp.xlabel("P-VALUE")
    mp.ylabel("Q-VALUE")
    mp.title("FDR CURVE")
    mp.plot(epi_x, fdrs_y, '.')
    mp.show()

def epi_vs_gain_volcano_plot(filtered_gain_snps, filtered_epi_snps, gain_vals, epi_vals, max_p, min_I3):
    gain_I3 = []
    gain_log_p = []
    for snps in filtered_gain_snps:
        gain_I3.append(gain_vals[snps])
        
        order = switch_snp_key_order(snps)
        if epi_vals.has_key(order[0]):
            gain_log_p.append(epi_vals[order[0]])
        elif epi_vals.has_key(order[1]):
            gain_log_p.append(epi_vals[order[1]])
    gain_log_p = -1*np.log10(gain_log_p)
    
    epi_I3 = []
    epi_log_p = []
    for snps in filtered_epi_snps:
        order = switch_snp_key_order(snps)
        if gain_vals.has_key(order[0]):
            epi_I3.append(gain_vals[order[0]])
        elif gain_vals.has_key(order[1]):
            epi_I3.append(gain_vals[order[1]])
            
        epi_log_p.append(epi_vals[snps])
    epi_log_p = -1*np.log10(epi_log_p)
    
    mp.figure(1)
    mp.xlabel("I3")
    mp.ylabel("-log10(P)")
    mp.title("Volcano plot - EPISTASIS and GAIN")
    mp.plot(epi_I3,epi_log_p,'bo')
    mp.plot(gain_I3,gain_log_p, 'ro')
    mp.axhline(y=(-1*np.log10(max_p)),linewidth=2, color='g')
    mp.axvline(x=min_I3, linewidth=2, color='g')
    #label max point
    max_x = np.max(gain_I3)
    max_y = np.max(gain_log_p)
    best_connection = ''
    #label best edge
    for snps in epi_vals:
        if(-1*np.log10(epi_vals[snps]) == max_y):
            best_connection = str(snps)
    mp.text(np.max(gain_I3),np.max(gain_log_p), best_connection, fontsize=10, horizontalalignment='center',verticalalignment='center',)
    mp.show()
    
    print 

def write_sif(dictionary,max_p,filename):
    outputfile =  filename+'_'+str(max_p)+'.sif'
    sifwriter = csv.writer(open(outputfile, 'wb'), delimiter='\t')
    for snps in dictionary:
        sifwriter.writerow([snps[0],str(dictionary[snps]),snps[1]])
    print 'wrote to file: ',outputfile 
#====================================================================================

# Create option parser
parser = optparse.OptionParser(usage="%prog [OPTIONS]", version="%prog 0.1")

# Add options to parser; use defaults if none specified
parser.add_option("--p-val-file", dest="p_val_filename", help="read data from PLINK epistasis file with p-vals", default="")
parser.add_option("--I3-file", dest="gain_filename", help="read data from gain file with I3 values", default="")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="verbose mode", default=False)
parser.add_option("--max-p", dest="max_p", help="maximum p-value", default=0.001)
parser.add_option("--min-I3", dest="min_I3", help="I3 threshold", default=0.0019)
parser.add_option("-o", "--output-sif", dest="outputfile", help="output sif filename", default="")
(options, args) = parser.parse_args()

#p_val_filename file check
if(options.p_val_filename!=""):
    try:
        epistasis_file = open(options.p_val_filename, 'r')
    except IOError as err:
        if err.errno != 2:
            raise err
        print "Error: PLINK epistasis file could not be opened."
else:
    print "Must specify PLINK epistasis filename!"
    parser.print_help()

#sif file check
if(options.gain_filename!=""):
    try:
        gain_file = open(options.gain_filename, 'r')
    except IOError as err:
        if err.errno != 2:
            raise err
        print "Error: gain file could not be opened."

verbose = options.verbose
max_p = float(options.max_p)
min_I3 = float(options.min_I3)

#create dictionaries
I3_dict, all_I3_dict = create_I3_dictionary(gain_file,min_I3,verbose)
epi_dict = create_epi_dictionary(epistasis_file,verbose)
fdr_dict = create_fdr_dictionary(epi_dict,verbose)
#draw fdr plot
#fdr_plot(epi_dict,fdr_dict)

#create filtered dictionaries
filtered_epi_dict, q_snp = filter_epi_dict(epi_dict, max_p)
q_value = fdr_dict[q_snp]
filtered_I3_dict = create_filtered_I3_dictionary(filtered_epi_dict,I3_dict)

#draw volcano plot
#epi_vs_gain_volcano_plot(filtered_I3_dict.keys(),filtered_epi_dict.keys(), all_I3_dict, epi_dict, max_p, min_I3)
epi_vs_gain_volcano_plot(filtered_I3_dict.keys(),epi_dict.keys(), all_I3_dict, epi_dict, max_p, min_I3)

if(options.outputfile=="epi"):
    write_sif(filtered_epi_dict,max_p,options.outputfile)
elif(options.outputfile=="gain"):
    write_sif(filtered_I3_dict,max_p,options.outputfile)
    
print 'GAIN - minimum I3 threshold = ', min_I3
print 'GAIN - initial # of connections = ', len(I3_dict)
print 'EPI - maximum p-value threshold = ', max_p
print 'EPI - initial # of connections = ', len(epi_dict)
print 'EPI - filtered by p-val connections left = ', len(filtered_epi_dict)
print 'GAIN - filtered by epi connections, final connections left = ', len(filtered_I3_dict)
print '% false = ', q_value*100
