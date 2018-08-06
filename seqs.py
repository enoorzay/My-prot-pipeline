import urllib.request
from bioservices.kegg import KEGG


HEADER = 'Query_Number' + '\t' + 'Uniprot_ID'+ '\t' + 'Description' + '\t' + 'e-Value' + '\t' + 'KEGG_ID'+ '\t' + 'PFAM'+ '\t' + 'Prosite' + '\t'
FINALHEADER = "Query\tDescription\tE-val\tPfam_Acc\tPfam\tProsite\tKEGGPathways\tGO\tComments"

#easy file io
class controller:
	#input file is a tab deliminated file matching the header above. output is generally same thing, just updates
	def __init__(self,inputfile, outputfile=""):
		self.inputfile = inputfile
		self.outputfile = outputfile
		self.seqs = list()
		if outputfile != 'a':
			self.readInput()
		else:
			self.readResults()
			
	#initialize swq objects for each 
	def readInput(self):
		nums_seen = set()
		with open(self.inputfile) as f:
			rawdata = f.read().split('\n')
		for i in rawdata:
			data= i.split('\t')
			if len(data) == 7:
				newseq = seq(data[0],data[1],data[2],data[3],data[4],data[5],data[6].strip())
				#alignments per query are ordered by score
				if data[0] not in nums_seen:
					newseq.isbest = True
					nums_seen.add(data[0])
				self.seqs.append(newseq)
	def readResults(self):
		with open(self.inputfile) as f:
			rawdata = f.read().split('\n')
		for i in rawdata:
			data= i.split('\t')
			if len(data) >= 5:
				newseq = seq(data[0],',',data[1],data[2],'.',data[4],data[5].strip())
				newseq.pfam_acc = data[3]
				newseq.pathways = data[6]
				newseq.go = data[7].split(',')
				newseq.isbest = True
				self.seqs.append(newseq)
		
	def writeOut(self,clean=False):
		if len(self.outputfile) == 0:
			print("Set output file with controllername.setOutputFile(f)")
			return
		with open(self.outputfile,'w') as f:
			print(HEADER,file=f)
			for i in self.seqs:
				if clean == True:
					if i.isbest:
						print(i.printLine(), file = f)
				else:
					print(i.printLine(), file = f)
	def setOutputFile(self,f):
		self.outputfile =f
	def getBestSeqs(self):
		
		outpt = list()
		for i in self.seqs:
			if i.isbest == True:
				outpt.append(i)
		return outpt
	
	def prepareData(self):
		outseqs = self.getBestSeqs()
		for i in outseqs:
			if i.kegg_id == '.':
				curr = self.getSeqs(i.id)
				min_e = 999
				compatible_kegg = '.'
				for j in curr:
					if j.kegg_id != '.': #and float(j.eval)< min_e:
						compatible_kegg = j.kegg_id
						#min_e = j.eval
				i.kegg_id = compatible_kegg
		self.seqs = outseqs		
		
	def getSeqs(self, id):
		results = list()
		for i in self.seqs:
			if i.id == id:
				results.append(i)
		return results
			
def createSeqs(self, inputfile):
	with open(filename) as f:
		rawdata = f.read().split('\n')
	data = list()
	for i in rawdata:
		data.append(i.split('\t'))
		
	return data
class seq:
	def __init__(self, id, accession,description,eval, kegg_id=".",pfam=".",prosite="."):
		self.id = id
		self.accession = accession
		self.description = description
		self.eval = eval
		self.kegg_id = kegg_id
		self.pfam = pfam
		self.prosite = prosite
		self.isbest = False
		self.go = '.'
		self.pfam_acc = '.'
		self.pathways = '.'
	def printLine(self):
		c = str(self.id)
		e = str(self.eval)
		return c + '\t' + self.accession + '\t' + self.description + '\t' + e+'\t' + self.kegg_id + '\t' + self.pfam + '\t' + self.prosite
	def printOut(self):
		print(self.id,self.accession,self.description, self.eval,self.kegg_id,self.pfam,self.prosite,sep='\t')
	'''
	def printKEGGline(self):
		return self.kegg_id+'\t'+self.orth + '\t' + self.path + '\t' + self.module + '\t' + self.fxn + '\t' + self.enzymes
	def setKEGGData(self,orth,path,module,fxn,enzymes):
		self.orth = orth 
		self.pathways = path
		self.module = module
		self.fxn = fxn
		self.enzymes = enzymes
'''
	def getKEGG(self):
		
		url = 'http://www.uniprot.org/uploadlists/'
		params = {
			'from':'ACC+ID',
			'to':'KEGG_ID',
			'format':'tab',
			'query':self.accession #'P13368 P20806 Q9UM73 P97793 Q17192'
		}
		data = urllib.parse.urlencode(params).encode("utf-8")
		request = urllib.request.Request(url, data)
		contact = "enoorzay@gmail.com" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
		request.add_header('User-Agent', 'Python %s' % contact)
		with urllib.request.urlopen(url,data=data) as url:
			s = url.read(200000)
		result = s.decode('utf-8').split('\t')
		if len(result) > 2:
			self.kegg_id= result[-1].strip()
			return True
		else:
			return False
	'''
	def requestKEGG(self):
		if self.kegg_id == '.':
			return False
		request = REST.kegg_get(self.kegg_id)
		open("kegg_outputs/"+self.kegg_id,'w').write(request.read())
		records = Enzyme.parse(open("kegg_outputs/"+self.kegg_id))
		record = list(records)[0]
		return record
		#url = "http://rest.kegg.jp/get/" + self.kegg_id
		#with urllib.request.urlopen(url) as f:
		#	s = f.read()
		#return s
	'''
	def __repr__(self):
		return self.printLine()
