# -*- coding: utf-8 -*-



# Load the necessary modules

import warnings, json
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML

# Load the configuration file

with open("test_arg_origin_detection.json") as json_conf : 
    conf = json.load(json_conf)

# Tell to NCBI who I am

Entrez.email = conf["email"]

def arg_origin_detection(idsfile,outfile,evalue):
    
    '''
    idsfile: TBLASTN Query Protein Id;
    outfile: output file;
    evalue: expected number of chance matches for tblastn;
    '''
    
    with open(idsfile) as f:
        prot_ids = f.readlines()
        prot_ids = [x.strip() for x in prot_ids]
    
    out_handle = open(outfile, "w")
    
    for query in prot_ids:
        result_handle = NCBIWWW.qblast('tblastn', 'nt', query, format_type="XML", hitlist_size= 500, alignments = 500, descriptions = 500, entrez_query="complete[Title] NOT plasmid[Title] NOT vector[Title] NOT integron[Title] NOT transposon[Title] NOT phage[Title] NOT virus[Title] NOT cds[Title] NOT island[Title] NOT conjugati[Title]")
        blast_records = NCBIXML.parse(result_handle)
        
        E_VALUE_THRESH = evalue
    
        hit_id, hit_description, evalue = [], [], []
        
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hit_id.append(alignment.accession)
                    hit_description.append(alignment.title)
                    evalue.append(hsp.expect)
                    
        for i in range(0,len(hit_id)):
            if i == 0:
                if hit_id[i] != hit_id[i+1]:
                    if evalue[i] < E_VALUE_THRESH:
                        out_handle.write(query+"\t"+hit_id[i]+"\t"+hit_description[i]+"\n")
            elif i == len(hit_id)-1:
                if hit_id[i] != hit_id[i-1]:
                    if evalue[i] < E_VALUE_THRESH:
                        out_handle.write(query+"\t"+hit_id[i]+"\t"+hit_description[i]+"\n")
            else:
                if hit_id[i] != hit_id[i-1] and hit_id[i] != hit_id[i+1]:
                    if evalue[i] < E_VALUE_THRESH:
                        out_handle.write(query+"\t"+hit_id[i]+"\t"+hit_description[i]+"\n")
        
    out_handle.close()

# Run the function

arg_origin_detection(conf["idsfile"],conf["outfile"],conf["evalue"])
