#! /usr/bin/env python

# -*- coding: iso-8859-1 -*-

# Importation de modules

from DialignFragments import *
from DialignAnchors import *
import re, sys, os, altgraph,itertools
#import boost.graph as bgl
# import networkx as NX
import string,math
import mySyntax

WFACTOR=10
## Zu aendern
BLOSUM='/gobics/home/eduardo/svn/unmask/branches/python_codes3/blosum62.txt'
##
class dialignWeights:
	""" Urgent need to simplify and rationalise this class """
	def __init__(self,scoreFile,fastaFile):
		self.fasta = FastaObject(fastaFile)
		#try:
		#	self.weights = scores(self.fasta,scoreFile)
		#except:
		#	cmd = "%s -wgtpr %s | sed 's/  */ /g' > %s" % ("$DIALIGN",fastaFile,"scores.txt")
		#	print cmd
		#	mySyntax.runAndReadAll(cmd)
		#	self.weights = scores(self.fasta,scoreFile)
		self.weights = pairscore("scores.txt")
		self.blosum = mySyntax.parseBlosum(BLOSUM)

	def weights(self,a,b):
		f = open('fasta.tfa','w')
		print >> f,'>s1'
		print >> f,'a'*int(a)
		print >> f,'>s2'
		print >> f,'a'*int(b)
		f.close()
		cmd = "%s -wgtpr %s | sed 's/  */ /g' > %s" % ("$DIALIGN","fasta.tfa","scores.txt")
		#print cmd
		mySyntax.runAndReadAll(cmd)
		f = FastaObject("fasta.tfa")
		return scores(f,"scores.txt")[(0,1)]

	def scorePairs(self,s1,s2,n,seqtype='prot'):
		""" Computes the Dialign score of the length n fragment starting at s1,s2 -- to be continued... """
		sq1 = self.fasta.printSeq(s1,n)
		sq2 = self.fasta.printSeq(s2,n)
		if seqtype == 'prot':
			bl = mySyntax.parseBlosum(BLOSUM)
			b = mySyntax.blosumScore(sq1,sq2,bl)  #blosumScore
		elif seqtype == 'dna':
			b = mySyntax.matchNb(sq1,sq2)
		else:
			raise TypeError
		scorePair = (n,b)
		a,b = int(s1[0]),int(s2[0])
		if a > b:
			a,b = b,a
		#idScore = self.weights(len(self.fasta.sites[a]),len(self.fasta.sites[b]))
		idScore = self.weights
		try:
			weight = idScore[scorePair]
			return float(weight)
		except KeyError:
			print "have to interpolate weight for fragment of length %d and \
			blosum score %d and sequences lengths %d and %d" % (n,b,len(self.fasta.sites[s1[0]]),len(self.fasta.sites[s2[0]]))
			return interpolate(idScore,n,b)

def scorePairs(f,s1,s2,n):
	""" Computes the Dialign score of the length n fragment starting at s1,s2 -- to be continued... """
	bl = mySyntax.parseBlosum(BLOSUM)
	sq1 = f.printSeq(s1,n)
	sq2 = f.printSeq(s2,n)
	b = mySyntax.blosumScore(sq1,sq2,bl)  #blosumScore
	scorePair = (n,b)
	#print scorePair
	#scorePair = "%3d %3d" % (n,b)
	#seq = f.sequences[:]            ### For computing simplicity, only use the reference scores.txt on the first 2 seq
	#seq.remove(f.sequences[s1[0]])  # remove the sequences where s1 and s2 are
	#seq.remove(f.sequences[s2[0]])
	#f.removeSeq(seq,separator=' ')             # construct the Fasta file with only the two sequences
	#fname = f.filename+".rm"
	#cmd = "dialign2-2 -wgtpr %s > %s" % (fname,"scores.txt")
	#mySyntax.runAndReadAll(cmd)
	a,b = int(s1[0]),int(s2[0])
	if a > b:
		a,b = b,a
	idScore = scores(f,"scores.txt")[(a,b)]
	#cmd1 = """grep "%s" scores.txt""" % str(scorePair)
	#try:
		#score = float(mySyntax.runAndReadAll(cmd1)[0].split()[2])
		#return score
	#except ValueError:
		#print "Offending command %s" % cmd1
	try:
		weight = idScore[scorePair]
		return float(weight)
	except KeyError:
		print "have to interpolate weight for fragment of length %d and \
blosum score %d and sequences lengths %d and %d" % (n,b,len(f.sites[s1[0]]),len(f.sites[s2[0]]))
		return interpolate(idScore,n,b)

def interpolate(idScore,n,b):
	weightList = []
	for i,j in idScore.keys():
		if j == b:
			weightList.append((i,idScore[(i,j)]))
	weightList.sort()
	w = weightList[-5:]
	x1,y1 = w[0]
	x2,y2 = w[-1]
	return (n-x1)*float(y2-y1)/float(x2-x1)+y1

def pairscore(scoreFile):
	try:
		sF = open(scoreFile,'r')
	except IOError:
		raise IOError	
	for line in sF.xreadlines():
		line = line.strip()
		if line and line.startswith("sequence"):
			valueDict = {}
		elif not line or line.startswith("weight"):
			pass
		else:
			l,bs,w = line.split()
			valueDict[(int(l),int(bs))] = float(w)
	return valueDict

def scores(f,scoreFile):
	a = map(len,f.sites.values())     # lengths of sequences
	indexPairs = []
	for i in range(len(a)):
		for j in range(len(a))[i+1:]:
			p = (i,j)
			indexPairs.append(p)
	indexScores = {}
	try:
		sF = open(scoreFile,'r')
	except IOError:
		raise IOError	
	for line in sF.xreadlines():
		line = line.strip()
		if line and line.startswith("sequence"):
			try:
				v = valueDict
				p = indexPairs.pop(0)
				indexScores[p] = v
			except NameError:
				pass
			valueDict = {}
		elif not line or line.startswith("weight"):
			pass
		else:
			l,bs,w = line.split()
			valueDict[(int(l),int(bs))] = float(w)
	p = indexPairs.pop(0)
	indexScores[p] = valueDict
	return indexScores

def main(fname=None,format=None,fasta=None,seqtype=None):
	if len(sys.argv) < 2:
		print "Usage: pairwiseAnchors.py AnchorFile(.gdf) alignerType(clustal/dialign)"
		sys.exit(1)
	fname = sys.argv[1]
	format = sys.argv[2]
	if len(sys.argv) >= 4:
		fasta = sys.argv[3]
	if len(sys.argv) >= 5:
		seqtype = sys.argv[4]
	else:
		seqtype = 'prot'
	f = FastaObject(fname,sep=' ')
	#print f.sites
	a = Anchors(filename=fname,anchorDict='suffix',sep=' ')
	#print a.anchors
	#for anc in a.anchors:
	#	print anc,a.anchors[anc]
	pairList = a.outputAnchorPairs()
	#for i in pairList:
	#	print i
	if format == 'clustal':
		ancName = mySyntax.subExt(fname,'dbc')
		g = open(ancName,'w')
		print >> g, "Ballast 1000"
		for anc in pairList:
			s1,s2,l = anc
			seq1,seq2 = f.sequences[s1[0]],f.sequences[s2[0]]
			pos1,pos2 = s1[1]+1,s2[1]+1
			length = l
			weight = WFACTOR*l
			line = "seq: %s %s pos: 0 beg: %d %d len: %d weight: %d" % (seq1,seq2,pos1,pos2,length,weight)
			print >> g, line
		g.close()
	elif format == 'dialign':
		if not fasta:
			ancName = mySyntax.subExt(fname,'dig')
			g = open(ancName,'w')
			for anc in pairList:
				s1,s2,l = anc
				seq1,seq2 = s1[0]+1,s2[0]+1
				pos1,pos2 = s1[1]+1,s2[1]+1
				length = l
				weight = l
				line = "%d %d %d %d %d %d" % (seq1,seq2,pos1,pos2,length,weight)
				print >> g, line
		else:        # if we have given a 3rd argument like "score" (or even something else... for the moment, give the undecoded Fasta file
			origFasta = FastaObject(fasta)
			ancName = mySyntax.fileBase(fname)+'w'+'.dig'  # the nicknames should then be rw and flw
			g = open(ancName,'w')
			#cmd = "dialign2-2 -wgtpr %s | sed 's/  */ /g' > %s" % (fname,"scores.txt") # replaces the commented part up there in scorePairs
			#mySyntax.runAndReadAll(cmd)
			lenList = map(len,origFasta.sites.values())
			meanLength = math.floor(reduce(lambda x,y:x+y,lenList)/len(lenList))
			if seqtype:
				cmd = "DialignScoreGen.py %d %d %s" % (meanLength,meanLength,seqtype)
			else:
				cmd = "DialignScoreGen.py %d %d" % (meanLength,meanLength)
			os.system(cmd)
			scoreWeights = dialignWeights("scores.txt",sys.argv[3])
			r1,r2 = len(origFasta.sites[0]),len(origFasta.sites[1])
			currentWeights = scoreWeights.weights
			for anc in pairList:
				s1,s2,l = anc
				seq1,seq2 = s1[0]+1,s2[0]+1
				R1,R2 = len(origFasta.sites[s1[0]]),len(origFasta.sites[s2[0]])
				pos1,pos2 = s1[1]+1,s2[1]+1
				length = l
				weight = scoreWeights.scorePairs(s1,s2,l,seqtype=seqtype)
				line = "%d %d %d %d %d %f" % (seq1,seq2,pos1,pos2,length,weight)
				print >> g, line
		g.close()
	else:
		msg = "Wrong type of format (clustal/dialign)"
		sys.exit(msg)
	print "Created file %s" % ancName

if __name__ == '__main__':
    main()
