# -*- coding: utf-8 -*-



# Load the necessary modules

import warnings
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Emboss.Applications import NeedleCommandline
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO, SeqIO, Entrez, AlignIO
from StringIO import StringIO

### Tell to NCBI who I am

Entrez.email = ""

 
def first_blast(database):

    """FolP, GlmM and Sul1-3 homologs are identified in complete GenBank sequences
       through BLASTP passing stringent e-value (1<e-20) and query coverage (75%) thresholds"""        

    COVERAGE = 75
       
    blastp_cline = NcbiblastpCommandline(query="query.fasta", db=database, evalue=1e-20, outfmt='"6 std qcovs"', num_alignments = 1000)
    stdout, sterr = blastp_cline()
         
    blast_records = SearchIO.parse(StringIO(stdout),"blast-tab",fields=['std','qcovs'])
    array = []
    for blast_record in blast_records:
        for alignment in blast_record:
            if alignment.query_coverage >= COVERAGE:
                if "P_" in alignment.id:
                    acc = alignment.id.split("_prot_")[1].split(".")[0]
                    acc = acc+".1"
                else:
                    acc = alignment.id.split("_prot_")[1].split("_")[0]                
                handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text")
                records = SeqIO.parse(handle, "gb")   
                for record in records:
                    print acc + "\t" + record.annotations["organism"]
                    break
                array.append(acc)
                array = list(set(array))
    
    return array


def remove_similars(ident): 

    """Remove sequences with desired % identity"""

    # Input multi-seq FASTA format

    validated_sequence, validated_id = [], []
    for record in SeqIO.parse("../input.fasta", "fasta"):
        validated_sequence.append(record.seq)
        validated_id.append(record.id.split(" ")[0])
    
    
    remove_id = []

    for i in range(0,len(validated_sequence)-1):
        
        if len(validated_sequence[i]) > 10:
        
    
            if validated_id[i] not in remove_id:
                for j in range (i+1, len(validated_sequence)):
                    if len(validated_sequence[j]) > 10:                  
                        needle_cline = NeedleCommandline(asequence="asis:"+validated_sequence[i], bsequence="asis:"+validated_sequence[j], gapopen=10, gapextend=0.5, outfile = 'stdout')
                        stdout, stderr = needle_cline()
                        alignment = AlignIO.parse(StringIO(stdout), "emboss")
                        for needle_records in alignment:
                            query = list(needle_records[0].seq)
                            subject = list(needle_records[1].seq)
                            matches = [h for h, k in zip(query, subject) if h == k]
                            while '-' in matches: matches.remove('-')   
                            similarity = (float(len(matches))/len(query))*100
        
                            
                            if similarity >= ident:    
                                remove_id.append(validated_id[j])
                                remove_id = list(set(remove_id))
        else:
            remove_id.append(validated_id[i])
            remove_id = list(set(remove_id))
                
                        
    records = (r for r in SeqIO.parse("../input.fasta", "fasta") if r.id.split(" ")[0] not in remove_id)
    SeqIO.write(records, "input2.fasta", "fasta")
    os.rename("input2.fasta","../output.fasta")
