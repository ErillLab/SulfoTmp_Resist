# -*- coding: utf-8 -*-



# Load the necessary modules

import warnings, subprocess, random, urllib2, httplib, json
from Bio.Blast.Applications import NcbiblastpCommandline
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO, SeqIO, Entrez
from Bio.SeqUtils import GC
from StringIO import StringIO

# Load the configuration file

import json

with open("test_arg_gc_amelioration.json") as json_conf : 
    conf = json.load(json_conf)

### Tell to NCBI who I am

Entrez.email = conf["email"]

def gc_amelioration(infile, outfile, blastpdatabase):
    
    ''' Computes systematically the GC content of ARG genes, MGE genome & host genome, 
        ARG genes are detected through BLASTP in local blastpdatabase using queryfile
    '''
    
    out_handle_F = open(outfile, "w")


    out_handle_F.write("Query\tOrganism\tOrganism %GC\tMGE Nucleotide Id\tMGE %GC\tARG Protein Id\tARG %GC\tIntegrase/Transposase %GC\n")

    plasmid_gc, organism_gc = {}, {}
    
    # Input multi-seq FASTA format
    
    sequences = {}
    for record in SeqIO.parse(infile, "fasta"):
        sequences[record.id.split(" [")[0]] = str(record.seq)
    
    hits = []
    
    for key1 in sequences.keys():
        
        out_handle = open("blastp_input.fasta", "w")
        out_handle.write(">"+key1+"\n")
        out_handle.write(sequences[key1])
        out_handle.close()
        
        # Each individual sequence is blasted against MGE complete sequences
        
        COVERAGE = 75
           
        blastp_cline = NcbiblastpCommandline(query="blastp_input.fasta", db=blastpdatabase, evalue=1e-20, outfmt='"6 std qcovs"', num_alignments = 1000)
        stdout, sterr = blastp_cline()
             
        blast_records = SearchIO.parse(StringIO(stdout),"blast-tab",fields=['std','qcovs'])
        prot_nucleotide_dict = {}
        for blast_record in blast_records:
            for alignment in blast_record:
                if alignment.query_coverage >= COVERAGE:
                    if "P_" in alignment.id:
                        acc_prot = alignment.id.split("_prot_")[1].split(".")[0]
                        acc_prot = acc_prot+".1"
                    else:
                        acc_prot = alignment.id.split("_prot_")[1].split("_")[0]
                        
                    acc_nuc = alignment.id.split("|")[1].split("_prot")[0]
                    
                    if acc_prot not in hits:
                        prot_nucleotide_dict[acc_prot] = acc_nuc
                        hits.append(acc_prot)
        
        
        # Get hit GC, MGE GC & organism GC
        
        for key2 in prot_nucleotide_dict.keys():
            
            c = None
            
            # Compute MGE GC
                
            if prot_nucleotide_dict[key2] not in plasmid_gc.keys():
                
                try:
                    
                    seq_record_mge = SeqIO.read(Entrez.efetch(db="nucleotide", id = prot_nucleotide_dict[key2], rettype="fasta"),"fasta")
                    plasmid_gc[prot_nucleotide_dict[key2]] = [float(GC(seq_record_mge.seq))]
    
                
                 # Compute Integron/Transposon GC
                    
                    handle_control = Entrez.efetch(db="nucleotide", id=prot_nucleotide_dict[key2], rettype="gbwithparts", retmode="txt")
                    records_control = SeqIO.parse(handle_control, "gb")
                    for record_control in records_control:
                        for f in record_control.features:
                            if f.type == "CDS":
                                if f.strand == 1 and "<" in str(f.location) or f.strand == -1 and ">" in str(f.location):
                                    continue
                                else:
                                    
                                    try:                                
                                        if "ntegr" in f.qualifiers["product"][0]:
                                            control_gc = float(GC(f.location.extract(record_control).seq))
                                            plasmid_gc[prot_nucleotide_dict[key2]].append(control_gc)
                                            break
                                        elif "ranspos" in f.qualifiers["product"][0]:
                                            control_gc = float(GC(f.location.extract(record_control).seq))
                                            plasmid_gc[prot_nucleotide_dict[key2]].append(control_gc)
                                            break
                                        
                                    except KeyError:
                                        continue
                                    
                    handle_control.close()
                    
                except (urllib2.HTTPError,httplib.IncompleteRead):
                    break
                
            # Compute the organism GC    
            
            try:
                
                handle_organism = Entrez.efetch(db="protein", id=key2, rettype="gb", retmode="text")
                records = SeqIO.parse(handle_organism, "gb")   
                for record in records:
                    organism = record.annotations["organism"]
                    break
                handle_organism.close()
                
                if organism not in organism_gc.keys():
                        
                        handle_organismgc = Entrez.esearch(db="nucleotide", term=organism+'[Primary Organism] AND complete genome[title]' , retmax = 500)
                        record = Entrez.read(handle_organismgc)
                        ids = record['IdList']
                        handle_organismgc.close()
                        
                        if len(ids) == 0:
                            handle_organismgc = Entrez.esearch(db="nucleotide", term=organism.split(" ")[0]+'[Primary Organism] AND complete genome[title]' , retmax = 500)
                            record = Entrez.read(handle_organismgc)
                            ids = record['IdList']
                            handle_organismgc.close()
                            
                        if len(ids) != 0:                   
                            seq_record_mge_organismgc = SeqIO.read(Entrez.efetch(db="nucleotide", id = ids[random.randint(0,len(ids)-1)], rettype="fasta"),"fasta")
                            organism_gc[organism] = float(GC(seq_record_mge_organismgc.seq))
                
            except (urllib2.HTTPError, httplib.IncompleteRead):
                break
                
            # Compute ARG GC
            
            protein_info = {}
            
            if "WP_" in key2:
                
                try:
                    
                    c = 0
                    records_ARG_ipg_gc = Entrez.read(Entrez.efetch(db="protein", id=key2, rettype='ipg', retmode='xml',retmax=1))
                    for record_ARG_ipg_gc in records_ARG_ipg_gc[0]['ProteinList']:
                        for nucleotide_id in record_ARG_ipg_gc:
                            r = nucleotide_id["CDSList"][0][0]
                            if "_" in r.attributes["accver"]:
                                protein_info[r.attributes["accver"]] = [key2,int(r.attributes["start"]),int(r.attributes["stop"])]
                                c += 1
                                break
                        if c == 1:
                            break
                        
                except (urllib2.HTTPError, IndexError, httplib.IncompleteRead):
                  break
                
            if not "WP_" in key2:
                
                try:
                    
                    records_ARG_gc = SeqIO.parse(Entrez.efetch(db="protein",id = key2, rettype="gb",retmode = "text"),"gb")
                    for record_ARG_gc in records_ARG_gc:
                        if "db_source" in record_ARG_gc.annotations:
                            try:
                                protein_info[record_ARG_gc.annotations["db_source"].split()[-1]] = [key2, int(record_ARG_gc.features[-1].qualifiers["coded_by"][0].split(":")[1].split("..")[0]),int(record_ARG_gc.features[-1].qualifiers["coded_by"][0].split(":")[1].split("..")[1].replace(")",""))]
                            except ValueError:
                                c = "Error"
                    
                    if c == "Error":
                        break
                
                except urllib2.HTTPError:
                    break
            
            try:  
                seq_record_ARG_gc = SeqIO.read(Entrez.efetch(db="nucleotide", id = protein_info.keys()[0], seq_start=protein_info[protein_info.keys()[0]][1], seq_stop=protein_info[protein_info.keys()[0]][2], rettype = "fasta"),"fasta")
                ARG_gc = float(GC(seq_record_ARG_gc.seq))
            
            except (urllib2.HTTPError,IndexError,httplib.IncompleteRead):
                break
                
            
            # Print All GC values
            
            try:
                if len(plasmid_gc[prot_nucleotide_dict[key2]]) == 1:
                    out_handle_F.write(key1+"\t"+organism+"\t"+str(organism_gc[organism])+"\t"+prot_nucleotide_dict[key2]+"\t"+str(plasmid_gc[prot_nucleotide_dict[key2]][0])+"\t"+key2+"\t"+str(ARG_gc)+"\n")
                else:
                    out_handle_F.write(key1+"\t"+organism+"\t"+str(organism_gc[organism])+"\t"+prot_nucleotide_dict[key2]+"\t"+str(plasmid_gc[prot_nucleotide_dict[key2]][0])+"\t"+key2+"\t"+str(ARG_gc)+"\t"+str(plasmid_gc[prot_nucleotide_dict[key2]][1])+"\n")
            except KeyError:
                continue
            
                    
        bashCommand = "rm blastp_input.fasta"
        subprocess.call(bashCommand.split())
    
    out_handle_F.close()


## Run the function

gc_amelioration(conf["infile"], conf["outfile"], conf["blastpdatabase"])
