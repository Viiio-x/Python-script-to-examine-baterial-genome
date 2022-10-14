#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Retrieving sequences from NCBI
# We choose Escherichia coli JE86-ST05 DNA, complete genome
# The sequence could be found at https://www.ncbi.nlm.nih.gov/nuccore/AP022815
from Bio import Entrez


# In[2]:


Entrez.email = "wanzix@andrew.cmu.edu"


# In[3]:


handle = Entrez.efetch(db="nucleotide", id="AP022815 ", rettype="fasta", retmode="text")


# In[4]:


#read sequence data and get basic information
from Bio import SeqIO


# In[5]:


record = SeqIO.read(handle, "fasta")


# In[6]:


print(record.description)


# In[7]:


print(record.id)


# In[8]:


print(record.name)


# In[10]:


record.seq


# In[11]:


#count the number of condons
len(record.seq)


# In[12]:


handle.close()


# In[13]:


from Bio.Seq import Seq


# In[19]:


#since the genome sequence is very large, we do a modification for current sequence
seq = str(record.seq[:1000])


# In[25]:


print(seq)


# In[20]:


#scan both genomic DNA strands for potential open reading frames longer than 50 codons
#use regular expression to find start codon and stop codon
import re


# In[21]:


pattern = re.compile(r'(?=(ATG(?:...)*?)(?=TAG|TGA|TAA))')


# In[22]:


#reverse complement
revcompseq = seq[::-1].maketrans("ATGC", "TACG")


# In[29]:


#find sequences with number of codons larger than 50
open_read_frame = []
l1 = len(pattern.findall(seq))
l2 = len(pattern.findall(seq[::-1].translate(revcompseq)))

#forward search
for i in range(l1):
    if len(pattern.findall(seq)[i]) > 50:
        open_read_frame.append(pattern.findall(seq)[i])

#backward search
for i in range(l2):
    if len(pattern.findall(seq[::-1].translate(revcompseq))[i]) > 50:
        open_read_frame.append(pattern.findall(seq[::-1].translate(revcompseq))[i])


# In[31]:


#print first three results
print (open_read_frame[:3])


# In[111]:


#Translate the potential open reading frames into protein sequences.
#obtain the complement or reverse complement of a Seq object using its built-in methods
#write into a new file 
table = 1
max_pro_len = 1000

protein_list = []
count = 0
for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
    for frame in range(3):
        for pro in nuc[frame:].translate(table).split("*"):
             if len(pro) < max_pro_len:
                count += 1
                protein_list.append(pro)
                print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame))


# In[113]:


print(protein_list[:10])


# In[114]:


from Bio import SeqUtils


# In[115]:


#write a function to evaluate the weight of protein
def fragment_weights(fragments):
    nums = []  
    for i in fragments:                                             
        fragment_weights = 0.0                                                     
        for pro in i:
            #change the unit to kD
            fragment_weights += SeqUtils.molecular_weight(pro, "protein")/1000
        #convert to two decimal numbers
        nums.append(round(fragment_weights, 2))
    print(nums)


# In[116]:


fragment_weights(protein_list)


# In[117]:


print(len(protein_list))


# In[118]:


#blast 5 protein coding sequences at NCBI and return most similar hits.
import Bio


# In[119]:


from Bio.Blast import NCBIWWW


# In[120]:


Entrez.email = "wanzix@andrew.cmu.edu"


# In[121]:


#select first five protein coding sequences
selected_protein_list = protein_list[0:5]


# In[122]:


print(selected_protein_list)


# In[123]:


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




