# Python-script-to-examine-baterial-genome

This project is based on the course 03701 Practical Computing for Biological Sciences at Carnegie Mellon University.

➡The script Final_Porject.py run a basic analysis of a bacterial genome analysis.
In this project, I select the species Escherichia coli JE86-ST05 as my analysis object. This script successfully did:

  a. Scan both genomic DNA strands for potential open reading frames longer than 50 codons.
  b. Translate the potential open reading frames into protein sequences.
  c. Calculate the predicted molecular mass of each protein (in kD).
  d. Blast 5 protein coding sequences at NCBI and return most similar hits.
    i.	**** Choose sequences < 1000 amino acids long.

➡️The report Project_Report.pdf is a further explanation of Final_Porject.py.
In this report, I presented my methods to achieve tasks above and code structure as well.
Results of this script and summary can also be found in this report.

➡️The script Genome_analysis.py explores a bit more basked on the Final_Porject.py
This Python script accomplishes the following tasks:
 a. Scan both genomic DNA strands for potential open reading frames longer than 50 codons querying PubMed
 b. Translate the potential open reading frames into protein sequences and write it to a csv file using Biopython
 c. Analyzing Residues Sequences including molecular weight, GRAVY, Amino acids and secondary structure using Protein Analysis 
    and add them to the file using Pandas
 d. Blast protein coding sequences at NCBI and return most similar hits using BLAST
 
➡️The script Genome_visualization.py is a visualization of whole genomes.
 a. import GenomeDiagram module visualizing genomes in a circular diagram from GenBank
 b. pick EcoRI recognition sites as a sensible caption for the features
 c. draw the diagram showing the position of restriction digest sites

Thank you for seeing! If you have any questions, feel free to contact with me! 😄
