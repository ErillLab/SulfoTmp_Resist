# -*- coding: utf-8 -*-


## Load the necessary modules

import os, subprocess, warnings, shutil
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SeqIO, Entrez
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from itertools import permutations

### Tell to NCBI who I am

Entrez.email = ""

### Codon alignment & Syn NoSyn computation

file1 = open('../Output.txt', 'w')
file1.write("Clade 1\tClade 2\tS\t↑%GC\tNS\t↑%GC\tS\t↑%GC\tNS\t↑%GC\tS\t↑%GC\tNS\t↑%GC\tS\t↑%GC\tNS\t↑%GC\n")

count = 1

SPECIES, SEQS = [], []

# Input CDS in fasta format. Multi-seq allowed.

for record in SeqIO.parse("/home/micromol/input.fasta", "fasta"):
    SPECIES.append(record.id)
    SEQS.append(record.seq)

for A in range(0,len(SEQS)-1):
    for B in range(A+1,len(SEQS)):
        file1.write(SPECIES[A]+"\t"+SPECIES[B]+"\t")
            
        key = SEQS[A]
        element= SEQS[B]
        
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
            ofile2.write(">species2\n" + str(key.translate(table=11,stop_symbol="",cds=True))+"\n")
            ofile2.close()
            
            ofile3 = open("../TEMP/prot1.fasta", "w")
            ofile3.write(">species1\n"+ str(element.translate(table=11,stop_symbol="",cds=True)) +"\n")
            ofile3.close()
            
            bashCommand = "needle -asequence ../TEMP/prot2.fasta -bsequence ../TEMP/prot1.fasta -outfile ../TEMP/align.fasta -aformat fasta -auto"
            subprocess.call(bashCommand.split())
            
            cmd = "perl ../PAL2NAL/pal2nal.pl ../TEMP/align.fasta ../TEMP/cds.fasta -nogap -codontable 11 -output fasta > ../TEMP/input.fasta"
            os.system(cmd) 
            
            # Syn NoSyn computation
    
            species, sequences = [], []        
            
            for record in SeqIO.parse("../TEMP/input.fasta", "fasta"):
                species.append(record.id)
                sequences.append(record.seq)
            
            for i in range(0,len(sequences)-1):
                for j in range(i+1,len(sequences)):
                    two_sequences = [sequences[i],sequences[j]]
                    syn_subst, nonsyn_subst, syn_subst_to_gc, nonsyn_subst_to_gc = [],[],[],[]
                    score = {"1":1,"2":0.5,"3":(3/18.0)}
                    codon_gc_dict, gc_dict, syn_position, nonsyn_position, syn_position_to_gc, nonsyn_position_to_gc = {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}
                    c = 1
                    for k in range(0,len(two_sequences[0]),3):
                        CODON, codon_1, codon_2 = two_sequences[0][k:k+3], two_sequences[0][k:k+3], two_sequences[1][k:k+3]
                        codon_1_ms, codon_2_ms = codon_1.tomutable(), codon_2.tomutable()
                        syn_subst_int, nonsyn_subst_int, syn_subst_to_gc_int, nonsyn_subst_to_gc_int = [], [], [], []
                        codon_gc_dict_int,gc_dict_int, syn_position_int, nonsyn_position_int, syn_position_to_gc_int, nonsyn_position_to_gc_int = {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}, {"1":0,"2":0,"3":0}
                        
                        #print "Codon %s"%str(c)+"\t"+codon_1+"->"+codon_1.translate(table="Bacterial")+"\t"+codon_2+"->"+codon_2.translate(table="Bacterial")
                        
                        c+=1
                        if codon_1 != codon_2:
                            perm = list(permutations([a for a in xrange(len(codon_1)) if codon_1[a] != codon_2[a]]))
                            for l in range(0,len(perm)):
                                CODON, codon_1, codon_2 = two_sequences[0][k:k+3], two_sequences[0][k:k+3], two_sequences[1][k:k+3]
                                codon_1_ms, codon_2_ms = codon_1.tomutable(), codon_2.tomutable()
                                for l2 in range(0,len(perm[l])):
                                    codon_1_ms[perm[l][l2]] = codon_2_ms[perm[l][l2]]
                                    if perm >= 2 and (Seq(str(codon_1_ms),generic_dna)).translate(table="Bacterial") == "*":
                                        break
                                    else:
                                        codon_gc_dict_int[str(perm[l][l2]+1)] = codon_gc_dict_int[str(perm[l][l2]+1)]+1
                                        if codon_1.translate(table="Bacterial") != (Seq(str(codon_1_ms),generic_dna)).translate(table="Bacterial"):
                                            nonsyn_subst_int.append(1)
                                            nonsyn_position_int[str(perm[l][l2]+1)] = nonsyn_position_int[str(perm[l][l2]+1)]+1
                                            if (CODON[perm[l][l2]] == "A" or CODON[perm[l][l2]] == "T") and (codon_1_ms[perm[l][l2]] == "G" or codon_1_ms[perm[l][l2]] == "C"):
                                                nonsyn_subst_to_gc_int.append(1)
                                                nonsyn_position_to_gc_int[str(perm[l][l2]+1)] = nonsyn_position_to_gc_int[str(perm[l][l2]+1)]+1
                                                gc_dict_int[str(perm[l][l2]+1)] = gc_dict_int[str(perm[l][l2]+1)]+1
                                        else:
                                            syn_subst_int.append(1)
                                            syn_position_int[str(perm[l][l2]+1)] = syn_position_int[str(perm[l][l2]+1)]+1
                                            if (CODON[perm[l][l2]] == "A" or CODON[perm[l][l2]] == "T") and (codon_1_ms[perm[l][l2]] == "G" or codon_1_ms[perm[l][l2]] == "C"):
                                                syn_subst_to_gc_int.append(1)
                                                syn_position_to_gc_int[str(perm[l][l2]+1)] = syn_position_to_gc_int[str(perm[l][l2]+1)]+1
                                                gc_dict_int[str(perm[l][l2]+1)] = gc_dict_int[str(perm[l][l2]+1)]+1
                                    codon_1 = Seq(str(codon_1_ms))
                                codon_1_ms, codon_2_ms = codon_1.tomutable(), codon_2.tomutable()
                            
                            # Score correction
                            
                            syn_subst.append(sum(syn_subst_int)*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int))))
                            nonsyn_subst.append((sum(nonsyn_subst_int)*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int)))))
                            syn_subst_to_gc.append((sum(syn_subst_to_gc_int)*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int)))))
                            nonsyn_subst_to_gc.append((sum(nonsyn_subst_to_gc_int)*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int)))))
                            for x in range(0,3):
                                gc_dict[str(x+1)] = gc_dict[str(x+1)]+(gc_dict_int[str(x+1)]*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int))))
                                codon_gc_dict[str(x+1)] = codon_gc_dict[str(x+1)]+(codon_gc_dict_int[str(x+1)]*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int))))
                                syn_position[str(x+1)] = syn_position[str(x+1)]+(syn_position_int[str(x+1)]*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int))))
                                nonsyn_position[str(x+1)] = nonsyn_position[str(x+1)]+(nonsyn_position_int[str(x+1)]*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int))))
                                syn_position_to_gc[str(x+1)] = syn_position_to_gc[str(x+1)]+(syn_position_to_gc_int[str(x+1)]*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int))))
                                nonsyn_position_to_gc[str(x+1)] = nonsyn_position_to_gc[str(x+1)]+(nonsyn_position_to_gc_int[str(x+1)]*(float(len(perm))/(sum(syn_subst_int)+sum(nonsyn_subst_int))))
                                
                    if int(sum(syn_subst)) != 0:
                        if int(sum(syn_subst_to_gc)) != 0:
                            file1.write(str(sum(syn_subst))+"\t"+(str(((sum(syn_subst_to_gc))/float(sum(syn_subst)))*100)+"\t"))
                        else:
                            file1.write(str(sum(syn_subst))+"\t0.0\t")
                            
                    if  int(sum(nonsyn_subst)) != 0:
                        if int(sum(nonsyn_subst_to_gc)) != 0:
                            file1.write(str(sum(nonsyn_subst))+"\t"+(str(((sum(nonsyn_subst_to_gc))/float(sum(nonsyn_subst)))*100)+"\t"))
                        else:
                            file1.write(str(sum(nonsyn_subst))+"\t0.0\t")
                    for y in range(0,3):
                        
                        if int(syn_position[str(y+1)]) != 0:
                            if int(syn_position_to_gc[str(y+1)]) != 0:
                                file1.write(str(syn_position[str(y+1)])+"\t"+str((syn_position_to_gc[str(y+1)]/float(syn_position[str(y+1)]))*100)+"\t")
                            else:
                                file1.write(str(syn_position[str(y+1)])+"\t"+str("0.0")+"\t")
                        if int(syn_position[str(y+1)]) == 0:
                            file1.write(str("0.0")+"\t"+str("0.0")+"\t")
                            
                        if int(nonsyn_position[str(y+1)]) != 0:
                            if int(nonsyn_position_to_gc[str(y+1)]) != 0:
                                file1.write(str(nonsyn_position[str(y+1)])+"\t"+str((nonsyn_position_to_gc[str(y+1)]/float(nonsyn_position[str(y+1)]))*100)+"\t")
                            else:
                                 file1.write(str(nonsyn_position[str(y+1)])+"\t"+str("0.0")+"\t")
                        if int(nonsyn_position[str(y+1)]) == 0:
                            file1.write(str("0.0")+"\t"+str("0.0")+"\t")
                            
                        if y == 2:
                            file1.write("\n")
                            
                    
        except CodonTable.TranslationError:
            continue
    
        shutil.rmtree('../TEMP')
        
file1.close()
