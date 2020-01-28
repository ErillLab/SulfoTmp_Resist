# -*- coding: utf-8 -*-

# Load the necessary modules
import warnings
from scipy.stats import mannwhitneyu
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SeqIO, AlignIO
from Bio.Emboss.Applications import NeedleCommandline
from StringIO import StringIO

# Load the configuration file

import json

with open("test_arg_origin_mwtest.json") as json_conf : 
    conf = json.load(json_conf)


def arg_origin_mwtest(infile, outfile, arg_source):

    ''' Computes systematically pair-wise amino acid identity among the products of chromosomal folA orthologs
        and compares this distribution with the pair-wise amino acid identity of the putative origin versus the chromosomal folA orthologs
        through one-sided Mann-Whitney U test
    '''
    
    out_handle = open(outfile,"w")
    
    source = arg_source
    mw_source = []
    mw_genus = []
    
    validated_ARG_sequence, validated_ARG_id = [], []
    for record in SeqIO.parse(infile, "fasta"):
        validated_ARG_sequence.append(record.seq)
        validated_ARG_id.append(record.id.split(" ")[0])
    
    for i in range(0,len(validated_ARG_sequence)):
        for j in range (i+1, len(validated_ARG_sequence)):
            if validated_ARG_id[i] != validated_ARG_id[j]:
                if len(validated_ARG_sequence[j]) > 10:                  
                    needle_cline = NeedleCommandline(asequence="asis:"+validated_ARG_sequence[i], bsequence="asis:"+validated_ARG_sequence[j], gapopen=10, gapextend=0.5, outfile = 'stdout')
                    stdout, stderr = needle_cline()
                    alignment = AlignIO.parse(StringIO(stdout), "emboss")
                    for needle_records in alignment:
                        query = list(needle_records[0].seq)
                        subject = list(needle_records[1].seq)
                        gaps = subject.count("-")+query.count("-")
                        matches = [h for h, k in zip(query, subject) if h == k]
                        while '-' in matches: matches.remove('-')   
                        similarity = (float(len(matches))/(len(query)-gaps)*100)
                    
                    out_handle.write(validated_ARG_id[i]+"\t"+validated_ARG_id[j]+"\t"+str(similarity)+"\n")
                    
                    if validated_ARG_id[i] == source or validated_ARG_id[j] == source:
                        mw_source.append(float(similarity))
                    else:
                        mw_genus.append(float(similarity))
    
    u_statistic, pVal = mannwhitneyu(mw_source, mw_genus)
    print pVal
    out_handle.write("Mann-Whitney U Test, p-value = "+str(pVal)+"\n")
    
    out_handle.close()

# Run the function
    
arg_origin_mwtest(conf["infile"], conf["outfile"], conf["arg_source"])
