# -*- coding: utf-8 -*-



# Load the necessary modules

from Bio import Entrez, SeqIO
from Bio.SeqUtils import GC
import numpy


# Tell to NCBI who I am

Entrez.email = ""

# Input Nucleotide IDs
    
with open("/home/micromol/input.txt") as f:
    content = f.readlines()
    content = [x.strip() for x in content]
    
for nuc_id in content:
    records = SeqIO.parse(Entrez.efetch(db="nucleotide", id=nuc_id, rettype="gb", retmode="txt"),"gb")
    organism, GC123, GC1, GC2, GC3 = "", [], [], [], []  
    for record in records:
        organism = record.annotations["source"]
        if int(record.annotations["contig"].split("..")[-1].split(")")[0]) >= 1000000:
            print "Species\tIdentifier\tProduct\tGC\tGC1\tGC2\tGC3"            
            handle = Entrez.efetch(db="nucleotide", id=nuc_id, rettype="gbwithparts", retmode="txt")
            records_2 = SeqIO.parse(handle, "gb")
            for record_2 in records_2:
                for f in record_2.features:
                    if f.type == "CDS":
                        GC1_cds, GC2_cds, GC3_cds, prot_id, product, sequence = [], [], [], "", "", ""
                        try:
                            protein_id, product = f.qualifiers["protein_id"][0], f.qualifiers["product"][0]
                            sequence = f.location.extract(record_2).seq
                            if len(sequence) % 3 == 0:
                                for j in range(0,len(sequence),3):
                                    codon = sequence[j:j+3]
                                    if codon[0] == "G" or codon[0] == "C":
                                        GC1_cds.append(1)
                                    if codon[1] == "G" or codon[1] == "C":
                                        GC2_cds.append(1)
                                    if codon[2] == "G" or codon[2] == "C":
                                        GC3_cds.append(1)
                                GC123.append(GC(sequence))
                                GC1.append(sum(GC1_cds)/(len(sequence)/3.0)*100)
                                GC2.append(sum(GC2_cds)/(len(sequence)/3.0)*100)
                                GC3.append(sum(GC3_cds)/(len(sequence)/3.0)*100)
                                print organism+"\t"+protein_id+"\t"+product+"\t"+str(GC(sequence))+"\t"+str(sum(GC1_cds)/(len(sequence)/3.0)*100)+"\t"+str(sum(GC2_cds)/(len(sequence)/3.0)*100)+"\t"+str(sum(GC3_cds)/(len(sequence)/3.0)*100)
                        except KeyError:
                            continue
    print organism+"\t"+nuc_id+"\tMean\t"+str(numpy.mean(GC123))+"\t"+str(numpy.mean(GC1))+"\t"+str(numpy.mean(GC2))+"\t"+str(numpy.mean(GC3))
    print organism+"\t"+nuc_id+"\tStandard Deviation\t"+str(numpy.std(GC123))+"\t"+str(numpy.std(GC1))+"\t"+str(numpy.std(GC2))+"\t"+str(numpy.std(GC3))    
    print organism+"\t"+nuc_id+"\tVariance\t"+str(numpy.var(GC123))+"\t"+str(numpy.var(GC1))+"\t"+str(numpy.var(GC2))+"\t"+str(numpy.var(GC3))+"\n"

                        
### 
# CSV format: species, gene, %GC1, %GC2, %GC3
###
