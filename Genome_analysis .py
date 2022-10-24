#!/usr/bin/env python
# coding: utf-8

# In[197]:


# Retrieving sequences from NCBI
# We choose Escherichia coli JE86-ST05 DNA, complete genome
# The sequence could be found at https://www.ncbi.nlm.nih.gov/nuccore/AP022815
from Bio import Entrez


# In[198]:


Entrez.email = "wanzix@andrew.cmu.edu"


# In[199]:


handle = Entrez.efetch(db="nucleotide", id="AP022815 ", rettype="fasta", retmode="text")


# In[68]:


#read sequence data and get basic information
from Bio import SeqIO


# In[69]:


record = SeqIO.read(handle, "fasta")


# In[70]:


print(record.description)


# In[71]:


print(record.id)


# In[72]:


print(record.name)


# In[73]:


record.seq


# In[74]:


len(record.seq)


# In[75]:


handle.close()


# In[76]:


from Bio.Seq import Seq


# In[110]:


#since the genome sequence is very large, we do a modification for current sequence
seq = record.seq[:1000]


# In[111]:


print(len(seq))


# In[112]:


import re


# In[113]:


pattern = re.compile(r'(?=(ATG(?:...)*?)(?=TAG|TGA|TAA))')


# In[115]:


#reverse complement
revcompseq = seq.reverse_complement()


# In[116]:


#find sequences with number of codons larger than 50
open_read_frame = []
l1 = len(pattern.findall(str(seq)))
l2 = len(pattern.findall(str(revcompseq)))


# In[119]:


#forward search
for i in range(l1):
    if len(pattern.findall(str(seq))[i]) > 50:
        open_read_frame.append(pattern.findall(str(seq))[i])


# In[121]:


#backward search
for i in range(l2):
    if len(pattern.findall(str(revcompseq))[i]) > 50:
        open_read_frame.append(pattern.findall(str(revcompseq))[i])


# In[122]:


print (open_read_frame[:3])


# In[126]:


import csv


# In[127]:


#Translate the potential open reading frames into protein sequences.
#obtain the complement or reverse complement of a Seq object using its built-in methods
#write into a new csv file 
table = 1
max_pro_len = 1000

f = open('protein_translation.csv','w')
header = ['sequence','length','strand','frame']
writer = csv.writer(f)
writer.writerow(header)

protein_list = []
count = 0
for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
    for frame in range(3):
        for pro in nuc[frame:].translate(table).split("*"):
             if len(pro) < max_pro_len:
                count += 1
                protein_list.append(pro)
                data = [str(pro[:30])+'...'+str(pro[-3:]),len(pro),strand,frame]               
                writer.writerow(data)

f.close()


# In[129]:


with open('protein_translation.csv','r') as f:
    #create interator
    for i , line in enumerate(f):
        #print each line
        print(line)
        if i > 5:
            break


# In[56]:


import pandas as pd


# In[171]:


df = pd.read_csv('protein_translation.csv')
# We have to clean the data, removing those whose length is equal to 0 for further calculation.
for x in df.index:
    if df.loc[x,"length"] == 0:
        df.drop(x,inplace = True)
print(df)


# In[172]:


#Now we want to run some protein analysis and store it in the file.
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# In[191]:


#first we creat lists as new columns in the csv file
#We calculate the molecular weight of each protein, the GRAVY(grand average of hydropathy) value of each protein 
#(this represents the property of hydrophobicity),number of each type of amino acid and their composition in each protein.
#We can also easily compute the basic information about secondary structure of each protein including beta sheets,
#alpha helixes, and turns.
molecular_weight = []
gravy = []
AA_count = []
AA_percent = []
Secondary_structure = []


# In[192]:


for pro in protein_list:
    if len(str(pro)) == 0:
        continue
    #get the sequence
    analyzed_seq = ProteinAnalysis(str(pro))
    #compute molecular weight of each protein
    molecular_weight.append(round(analyzed_seq.molecular_weight(),2))
    #compute the gravy value of each protein
    gravy.append(round(analyzed_seq.gravy(),2))
    #compute number of each type of amino acid of each protein
    AA_count.append(analyzed_seq.count_amino_acids())
    #compute composition of each type of amino acid in each protein
    AA_percent.append(analyzed_seq.get_amino_acids_percent())
    #compute the fraction of amino acids in secondary structure
    Secondary_structure.append(analyzed_seq.secondary_structure_fraction())


# In[193]:


#add new column to existing dataframe in pandas
df['molecular_weight'] = molecular_weight
df['gravy'] = gravy
df['AA_count'] = AA_count
df['AA_percent'] = AA_percent
df['Secondary_structure'] = Secondary_structure


# In[195]:


print(df)


# In[196]:


#save the modified content into file
df.to_csv('protein_translation.csv',index = False)


# In[200]:


#blast 5 protein coding sequences at NCBI and return most similar hits.
import Bio


# In[205]:


from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# In[206]:


Entrez.email = "wanzix@andrew.cmu.edu"


# In[207]:


#select first five protein coding sequences
selected_protein_list = protein_list[0:5]


# In[217]:


Blast_results = {}
for i in selected_protein_list:
    #store each blast result into a dictionary
    Blast_results[i] = {}
    results = NCBIWWW.qblast(program = "blastp", database = "nr", sequence = i, format_type= "XML")
    records = NCBIXML.parse(results)
    for r in records:
        for alignment in r.alignments:
            for hsp in alignment.hsps:
                matched_seq = hsp.sbjct
                Blast_results[i]["score"] = hsp.score
                #compare the scores in blast results, select the largest one(most similar hit)
                if matched_seq not in Blast_results[i].keys() and hsp.score >= Blast_results[i]["score"]:
                    Blast_results[i]["matched_sequence"] = matched_seq
                    Blast_results[i]["score"] = hsp.score
                    Blast_results[i]["description"] = alignment.title
                    Blast_results[i]["length:"] = alignment.length
    #print most similar result for each protein sequence
    print("****Alignment****")
    print("sequence:", alignment.title)
    print("length:", alignment.length)
    print("score:", hsp.score)


# In[ ]:




