''' 
USAGE: init.py [inputfile] [outputfile]
Parses blaster.py output file. Collects top 3 BLAST hits (and finds kegg ids for these) for each sequence. 
If no kegg ids are found, continues to collect BLAST results and searching for kegg ids until one is found or we reach insignificant hits (E >= 1)
-Stores each hits' data as a 'seq' object defined in seqs.py 
-Outputs a tab deliminated table with some columns containing temp values to be filled later. 
-This table can easily be read and edited via the controller object defined in seqs.py 
'''
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys
from seqs import *
from Bio import SeqIO


HEADER = 'Query_Number' + '\t' + 'Uniprot_ID'+ '\t' + 'Description' + '\t' + 'e-Value' + '\t' + 'KEGG_ID'+ '\t' + 'PFAM'+ '\t' + 'Prosite' + '\t'

if len(sys.argv) >= 3:
	input_file = sys.argv[1]
	output_file = sys.argv[2]
else:
	sys.exit()

f = open(output_file,'w')
print(HEADER,  file=f)

#list where each index contains data about the corresponding query
#The data is in the form of a list of 3 or so alignments, in seq objects
queries = list() 

result=open(input_file,"r")
records= NCBIXML.parse(result)


querycount = 1
alignmentcount = 0

for item in records:
	good_alignments = list() #holds 3 best alignments (as seq items) for each query to be stored
	keggFound = False #whether we have a kegg id for a match, so long as matches are above a certain eval keep searching till one gives a keggid
	for alignment in item.alignments:

		eval = alignment.hsps[0].expect
		
		#process components to build seq obj
		strlist = alignment.title.split('|')
		sp = strlist[3]
		description = strlist[4].split('Full=')[1].split(';')[0]
		description = description.split('>gi')[0].strip() #Fixes weird issue where some descs have '>gi' at end
		newseq = seq(querycount, sp, description,eval)
		
		#doing it like this so that once keggfound is set true it wont change for the rest of the curr query processing. 1 kegg id is enough
		currfound= newseq.getKEGG()
		if keggFound == False:
			keggFound = currfound
		
		
		good_alignments.append(newseq)
		f.write(newseq.printLine() + '\n')
		
		
		#max 3 alignments stored per query before moving on, dont need more than that i think
		alignmentcount+=1
		
		#if 3 alignments and found a kegg id we can be done with this query
		if alignmentcount >= 1:
			if keggFound == True:
					keggFound = False
					alignmentcount = 0 
					break
			# if no kegg id found, 3 alignments saved and getting bad evals its time stop
			elif eval > 1:
					alignmentcount = 0
					break

				
	queries.append(good_alignments)
	querycount+=1
	
		
f.close()