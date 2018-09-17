# -*- coding: utf-8 -*-


## Load the necessary modules

import os, subprocess, warnings, shutil
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SeqIO, Entrez
from Bio.Data import CodonTable
from datetime import datetime
startTime=datetime.now()

## Tell to NCBI who I am

Entrez.email = ""

### Codon alignment & Syn NoSyn computation

file1 = open('/home/micromol/Output.txt', 'w')
file1.write("S\tN\tdN/dS\tdN\tdS\n")
        
species, sequences = [], []

# Input CDS in FASTA format. Multi-seq allowed.

for record in SeqIO.parse("/home/micromol/input.fasta", "fasta"):
    species.append(record.id)
    sequences.append(record.seq)
    
orthologs_cds, species_cds = {}, {}

for a1 in range(0,len(sequences)-1):
    for a2 in range(a1+1,len(sequences)):
        orthologs_cds, species_cds = {}, {}
        species_cds[species[a1]] = species[a2]
        orthologs_cds[sequences[a1]] = sequences[a2]
        
        
        for (key,element), (k2,v2) in zip(orthologs_cds.items(), species_cds.items()):
            
            cmd = "echo %s %s | tr '\n' '\t' >> /home/micromol/Output.txt"% (k2,v2)
            os.system(cmd)
            
            bashCommand = "mkdir ../TEMP"
            subprocess.call(bashCommand.split())
                
            try:
        
                # PAL2NAL codon alignment 
                # Mikita Suyama, David Torrents, and Peer Bork (2006)
                # PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments.
                # Nucleic Acids Res. 34, W609-W612.
            
                ofile1 = open("../TEMP/cds.fasta", "w")
                ofile1.write(">species2\n" + str(key) + "\n>species1\n"+str(element)+"\n")
                ofile1.close()
                
                ofile2 = open("../TEMP/prot2.fasta", "w")
                ofile2.write(">species2\n" + str(key.translate(table=11,stop_symbol=""))+"\n")
                ofile2.close()
                
                ofile3 = open("../TEMP/prot1.fasta", "w")
                ofile3.write(">species1\n"+ str(element.translate(table=11,stop_symbol="")) +"\n")
                ofile3.close()
                
                bashCommand = "needle -asequence ../TEMP/prot2.fasta -bsequence ../TEMP/prot1.fasta -outfile ../TEMP/align.fasta -aformat fasta -auto"
                subprocess.call(bashCommand.split())
                
                cmd = "perl ../PAL2NAL/pal2nal.pl ../TEMP/align.fasta ../TEMP/cds.fasta -nogap -codontable 11 -output paml > ../TEMP/input.paml"
                os.system(cmd)               
                
                # codeml Syn NoSyn computation
                
                os.chdir("/home/micromol/PAML")
                
                cmd = "bin/codeml"
                os.system(cmd) 
                
                cmd = "grep 't=' mlc | awk -F'=' '{print $3,$4,$5,$6,$7}' OFS='\t'| tr -d '[:alpha:]' | sed -e 's/\///gi'| sed -e 's/  //gi' | sed -e 's/ //gi' >> /home/micromol/Output.txt"
                os.system(cmd) 
                
                os.chdir("/home/micromol/.spyder2")
                                                   
            except CodonTable.TranslationError:
                continue
        
            shutil.rmtree('../TEMP')
            
file1.close()
