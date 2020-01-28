# -*- coding: utf-8 -*-



# Load the necessary modules


import warnings, json
from Bio.SeqUtils import GC
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SeqIO, Entrez

# Load the configuration file

with open("test_arg_summary-env.json") as json_conf : 
    conf = json.load(json_conf)


# Tell to NCBI who I am

Entrez.email = conf["email"]


## Protein info, genetic environment & GC% computation

# Input: Txt file containing Protein Ids

file_name = conf["idsfile"]

with open(file_name) as f:
    ids_array = f.readlines()
    ids_array = [x.strip() for x in ids_array]

# Construction of a dictionary with GenBank information of each Protein Id
# protein_info[NucId] = [species, description, ProteinId, prot_seq_len, seq_start, seq_end, strand]

protein_info = {}

for id_protein in ids_array:
    if "WP" not in id_protein:
        recordsA = SeqIO.parse(Entrez.efetch(db="protein",id = id_protein, rettype="gb",retmode = "text"),"gb")
        for recordA in recordsA:
            if "db_source" in recordA.annotations:
                protein_info[recordA.annotations["db_source"].split()[-1]] = [recordA.annotations["source"], recordA.description.split(" [")[0], id_protein, int(len(recordA.seq)), int(recordA.features[-1].qualifiers["coded_by"][0].split(":")[1].split("..")[0]),int(recordA.features[-1].qualifiers["coded_by"][0].split(":")[1].split("..")[1].replace(")",""))]
                if not "complement" in recordA.features[-1].qualifiers["coded_by"][0]:
                    protein_info[recordA.annotations["db_source"].split()[-1]].append("1")
                else:
                    protein_info[recordA.annotations["db_source"].split()[-1]].append("2")
    else:
        c = 0
        recordsA = Entrez.read(Entrez.efetch(db="protein", id=id_protein, rettype='ipg', retmode='xml',retmax=1))
        for recordA in recordsA[0]['ProteinList']:
            for nucleotide_id in recordA:
                r = nucleotide_id["CDSList"][0][0]
                if "_" in r.attributes["accver"]:
                    protein_info[r.attributes["accver"]] = [r.attributes["org"],"",id_protein,"",str(r.attributes["start"]),str(r.attributes["stop"])]
                    c += 1
                    break
            if c == 1:
                break
            
        recordsB = SeqIO.parse(Entrez.efetch(db="protein",id = id_protein, rettype="gb",retmode = "text"),"gb")
        for recordB in recordsB:
            protein_info[r.attributes["accver"]][1] = recordB.description.split(" [")[0]
            protein_info[r.attributes["accver"]][3] = int(len(recordB.seq))
            if "+" in r.attributes["strand"]:
                protein_info[r.attributes["accver"]].append("1")
            else:
                protein_info[r.attributes["accver"]].append("2")


def gene_summary():
    
    ''' This function prints protein_description, species, the Nucleotide Id, protein length and GC% of each protein contained in the input file '''
    
    print "Protein_description,Species,Nucc_id,Protein_id,Protein_length,GC%"
        
    for ids_CDS in protein_info.keys():
        seq_recordB = SeqIO.read(Entrez.efetch(db="nucleotide", id = ids_CDS, rettype="fasta", seq_start=protein_info[ids_CDS][-3], seq_stop=protein_info[ids_CDS][-2]),"fasta")
        GC_cds = GC(seq_recordB.seq)
        
        print protein_info[ids_CDS][1].replace(",","")+","+protein_info[ids_CDS][0]+","+ids_CDS+","+protein_info[ids_CDS][2]+","+str(protein_info[ids_CDS][3])+","+str(round(GC_cds,2))


def gene_env(flanking_region):
    
    ''' This function prints the genetic environment of each protein contained in the input file
        flanking_region option indicates the 3' & 5' bp retained for each gene'''
    
    for ids_CDS in protein_info.keys():
        product = []
        seq_records = SeqIO.parse(Entrez.efetch(db="nucleotide", id = ids_CDS, rettype="gb", seq_start=int(protein_info[ids_CDS][-3])-int(flanking_region), seq_stop=int(protein_info[ids_CDS][-2])+int(flanking_region),strand=int(protein_info[ids_CDS][-1])),"gb")
        for seq_record in seq_records:
            for f in seq_record.features:
                if f.type == "CDS":
                    if f.strand == 1 and "<" in str(f.location) or f.strand == -1 and ">" in str(f.location):
                        continue
                    else:
                        try:                                
                            product.append(f.qualifiers["product"][0])
                        except KeyError:
                            continue
        
        print protein_info[ids_CDS][0]+","+ids_CDS+","+protein_info[ids_CDS][2]+","+",".join(product)

# Run the functions
        
gene_summary()
gene_env(conf["flanking_region"])
