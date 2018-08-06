# My-prot-pipeline
Takes in a FASTA file containing proteomic sequences and outputs accession number, Prosite id, Pfam id, as well as KEGG and GO information in a tab-separated file.  

Run by driver.py, which calls the script in the appropriate order. blaster.py -> init_table.py -> table_extender.py
- blaster.py queries the first set of databases and outputs raw xml files for BLAST and PROSITE and isolated prosite ids in a line separated file.
- init_table filters the xml files for relevant results, and queries a few additional databases. Outputs a table that can be read by seq.py's controller object. 
- table_extender.py uses seq.py's structures to read this table, add the final touches and outputs final table. 

-seq.py contains a structure to hold each query sequence and the relevant data pertaining to it. A controller class is also defined that  makes it easy to manage and perform operations on these seq objects (as well as file io). This approach is used so that I can easily modify the capabilities of this program after basic database querying, corresponding to the current dataset and objective I have at hand.

- A local go slim database is used to replace GO terms with GO-slim terms if possible. (GO-slim terms tend to lean less towards inaccuracy as they are more generalized.) 
