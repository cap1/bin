#! /usr/bin/env python

# -*- coding: iso-8859-1 -*-

"""
This module performs the computation of Dialign-formatted anchor files that
correspond to the concepts of partial columns and consistent partial columns
defined in "A min-cut Algorithm for the Consistency Problem in Multiple
Sequence Alignment" by E. Corel, F. Pitschi and B. Morgenstern. For
installation and use, see the INSTALL.txt file attached.
"""

# Python Modules required

# Local -- 

import DialignAnchors,IncidenceGraph
import IncidenceGraphBoost
from fragments2anchors import viewGraph

# to install -- networkX 0.99

import networkx as nx

# standard -- with Python v2.4 or higher

import re, sys, os, itertools, time
from optparse import OptionParser
import inspect,math

cmd = "export DIALIGN=/usr/bin/dialign2-2"
os.system(cmd)

class AnchorPairs:

	def __init__(self,frgFile=None,thr=0):
		self.frgFile = frgFile
		self.frgPairs = frgPairs(self.frgFile,thr=thr)
		self.frgGraph = self.frgGraph()

	def frgGraph(self,anchor=None):
		g = nx.Graph()
		if anchor:
			g.add_edges_from(anchor)
		else:
			g.add_edges_from(self.frgPairs)
		return g

def frgPairs(fname,thr=0):
	""" For a dialign anchor formatted file, compute the pairs of sites """
	f = open(fname,'U')
	pairwiseFrg = []
	for line in f.xreadlines():
		pairwiseFrg.extend(frgList(line,thr=thr))
	f.close()
	pairwiseFrg.sort()
	return pairwiseFrg

def frgList(line,thr=0):
	" Reads a line in Dialign's .anc format"
	outList = []
	anclist = line.strip().split()
	s1,s2 = map(lambda x:x-1,map(lambda x:int(x),anclist[0:2]))
	b1,b2 = map(lambda x:x-1,map(lambda x:int(x),anclist[2:4]))
	lg = int(anclist[4])
	w = float(anclist[5])     # dialign weight
	if w >= thr:
		for i in range(lg):
			frg = [(s1,b1+i),(s2,b2+i),w]
			outList.append(frg)
	return outList


class IncGraph(nx.Graph):
	""" Computes the incidence graph corresponding to an existing .gdf file """

	def __init__(self,gdfFile=None,env=False):
		nx.Graph.__init__(self)
		if gdfFile:
			a = DialignAnchors.Anchors(filename=gdfFile,anchorDict='suffix',sep=' ')
			n = 0
			for key,value in a.anchors.iteritems():
				if env:
					n = int(key.split("_")[1])-1  # order of decoding if the file comes from a decoding procedure !!!
				for i in range(2*n+1):
					c = map(lambda (x,y):(x,y-n+i),value)
					clique = []
					for site in c:
						seq,pos = site
						if 0 <= pos < len(a.sites[seq]):
							clique.append(site)
					self.add_clique(clique)
	def add_clique(self,clique):
		while clique:
			vertex = clique.pop()
			for tail in clique:
				self.add_edge(vertex,tail)

	def copy(self):
		g = IncGraph()
		for edge in self.edges(data=True):
			g.add_edge(*edge)
		return g
				  

##== Syntaxic functions

def sortList(itemList,function=None):
	""" Goes around a funny behaviour of Python's sort -- when list elements are lists or tuples """
	itemDict = dict([(i,itemList[i]) for i in range(len(itemList))])
	functionDict = {}
	for i in range(len(itemList)):
		functionDict[i] = function(itemList[i])
	vDict = dictInvert(functionDict)
	valueDict = dict([(k,[itemDict[j] for j in vDict[k]]) for k in vDict.keys()])
	ks = valueDict.keys()
	ks.sort()
	outList = []
	for i in ks:
		for L in valueDict[i]:
			outList.append(L)
	return outList

def dictInvert(d):
	"""Dictionary inversion routine"""
	inv = {}
	for k, v in d.iteritems():
		keys = inv.setdefault(v, [])
		keys.append(k)
	return inv

def subExt(fname,newExt=None):
	"""Change the extension of a file"""
	lst = fname.split(".")[0:-1]
	string = ''
	for l in lst:
		string += l+"."
	string += newExt
	return string

def filterGraph(graph,t=0,mode='degree'):
	if mode == 'degree':
		sg = filterNodes(graph,deg=t)
	elif mode == 'edge':
		sg = filterEdges(graph,weight=t)
	return sg

def filterNodes(graph,deg=0):
	""" Returns the subgraph filtered by degree value"""
	nodes = (n for n in graph.nodes() if graph.degree(n) >= deg)
	subgraph = graph.subgraph(nodes)
	return subgraph

def filterEdges(graph,weight=0):
	""" Returns the subgraph filtered by edge weight"""	
	subgraph = nx.Graph()
	edges = (edge for edge in graph.edges(data=True) if graph.get_edge(*edge) >= weight)
	subgraph.add_edges_from(edges)
	return subgraph

def hasRep(graph):
	"""Checks for colour repetitions in graph -- Boolean"""
	seqList = []
	for n in graph.nodes():
		seqList.append(n[0])   # add the sequence where n is to seqList
	return len(seqList) != len(set(seqList))  # if there is a repetition, then set and list have different cardinalities

def getRep(graph):
	""" Dictionary seq : list of repeated positions in seq in connected graph """
	colour = {}
	for n in graph.nodes():
		try:
			colour[n[0]].append(n)
		except KeyError:
			colour[n[0]] = [n]
	delKeys = []
	for i,v in colour.iteritems():
		if len(v) <= 1:
			delKeys.append(i)
		else:
			colour[i].sort()
	for i in delKeys:
		del colour[i]
	return colour

def edgeRemoval(g,m=0,filt='degree',weight=None,mode=None,view=None):
	"""Main procedure: iteratively filters and calls the max-flow-min-cut algorithm for spurious edge detection
	-- if weight, uses Dialign weights as flow capacities
	-- if mode, uses Boost Edmods-Karp read-write DIMACS procedure from IncidenceGraphBoost.py"""
	if not getRep(g):
		return []
	print "Processing ambiguous connected component of size %d" % len(g)
	graph = g.copy()
	removedEdges = []
	if filt == 'degree':
		maxDegree = max([graph.degree(n) for n in graph.nodes()])
		k = maxDegree
	elif filt == 'edge':
		maxWeight = max([graph.get_edge(*e) for e in graph.edges()])
		k = math.floor(maxWeight)
	else:
		raise TypeError
	sgList = []
	while k >= m:
		subgraph = filterGraph(graph,t=k,mode=filt)
		compList = nx.connected_component_subgraphs(subgraph)
		goodComps = [s for s in compList if not hasRep(s)]
		badComps = [s for s in compList if hasRep(s)]
		if badComps or k<=0:
			if filt == 'degree':
				print "Filtering nodes of degree >= %d" % k			
			elif filt == 'edge':
				print "Filtering edges of weight >= %d" % k
			if (k < 0 and filt == 'edge') or (k == 0 and filt == 'degree'):
				break
			else:
				print "Processing %d ambiguous subgraphs of sizes" % len(badComps)
				print map(len,badComps)
				while badComps:
					s = badComps.pop()
					#edges,good,bad = EK_edge_removal(s,weight=weight) #V0
					edges,good,bad = EK_edge_removal(graph,s,weight=weight,mode=mode,view=view) #V1
					removedEdges.extend(edges)
					graph.remove_edges_from(edges)
					if bad:
						badComps.extend(bad)
		graph.remove_edges_from(removedEdges)
		k-= 1
	return removedEdges

def EK_edge_removal(graph,subgraph,weight=None,mode=None,view=None):
	"""Computes the set of edges to remove according to the Edmonds-Karp max-flow algorithm -- variant with graph is experimental"""
	removedEdges = []
	while True:
		reps = getRep(subgraph)
		if not reps:
			break
		repetitions = sortList(reps.values(),function=len) # sorts the repetitions by cardinality
		kRep = repetitions.pop()   # choose which sequence having most number of repetitions
		rep = sortList(kRep,function=subgraph.degree) # sorts the repetitions by cardinality
		s,t = rep[0],rep[-1]
		if mode == 'subgraph':
			del_edges = IncidenceGraphBoost.EdmondsKarp(subgraph,source=s,sink=t,weight=weight)  # variant with C++ Boost
		elif mode == 'graph':
			del_edges = IncidenceGraphBoost.EdmondsKarp(graph,source=s,sink=t,weight=weight)  # variant with C++ Boost
		else:
			del_edges = IncidenceGraph.EdmondsKarp(subgraph,source=s,sink=t,weight=weight)
		if view:       # temporary viewing options
			sN = subgraph.nodes()
			if view == 'boundary':
				eN = nx.node_boundary(graph,sN)
				sN.extend(eN)
			extSg = nx.subgraph(graph,sN)
			viewGraph(extSg,nodes=[s,t],subset=del_edges)
		removedEdges.extend(del_edges)
		subgraph.remove_edges_from(del_edges)
		graph.remove_edges_from(del_edges)
		if not nx.is_connected(subgraph):
			compList = nx.connected_component_subgraphs(subgraph)			
			goodComp = [s for s in compList if not hasRep(s)]
			badComp = [s for s in compList if hasRep(s)]
			break
	return removedEdges,goodComp,badComp

def Main(*args,**kwargs):
	"""Takes a fasta file and an optional threshold value (default 0), and returns the dictionary of
	positions that form an equivalence class"""
	fname = args[0]
	thr = kwargs['thr']
	t0 = time.clock()
	rmFiles = []
	print "Computing incidence graph..."
	if kwargs['gdf']:
		g = IncGraph(gdfFile=kwargs['gdf'],env=kwargs['env'])
	else:
		if kwargs['frg']:   # if we use an already defined anchor file
			frgName = kwargs['frg']
		else:               # compute the pairwise fragments as an anchor file
			cmd = "$DIALIGN -ff %s" % fname
			print "Running Dialign..."
			os.system(cmd)
			print "Finished"
			frgFile = fname+".frg"
			frgName = fname+".ff"
			cmd1 = """grep ")" %s | sed 's/  */ /g' | cut -d" " -f4,5,7,8,10,12 | sort -n > %s""" % (frgFile,frgName)
			os.system(cmd1)
			aliName = fname+".ali"
			rmFiles = [frgFile,aliName]
		ap = AnchorPairs(frgName,thr=thr)
		g = ap.frgGraph
	print "Finished"
	N,E = (len(g.nodes()),len(g.edges())) 
	del_edges = []
	sList = nx.connected_component_subgraphs(g)
	kk = 0
	print """The initial incidence graph has %d connected components of sizes""" % len(sList)
	print map(len,sList)
	for sg in sList:
		del_edges.extend(edgeRemoval(sg,filt=kwargs['filt'],weight=kwargs['wg'],mode=kwargs['bs'],view=kwargs['vw']))
	g.remove_edges_from(del_edges)
	print """The initial incidence graph has %d nodes and %d edges""" % (N,E)
	print """The final incidence graph has %d nodes and %d edges""" % (len(g.nodes()),len(g.edges())) 
	print "Removed %d edges" % len(del_edges)
	comp = nx.connected_component_subgraphs(g)
	k = 0
	anchorDict = {}
	for s in comp:
		if len(s) > 1:
			key = str(k)+"_"+str(int(thr))
			anchorDict[key] = s.nodes()
			k += 1
	a = DialignAnchors.Anchors(filename=fname,anchorDict=anchorDict)
	outfile = fname+"_t"+str(int(thr))+"_pc.gdf"
	a.recoverAnchors(outbase=fname+"_t"+str(int(thr))+"_pc")
	anchorFiles = [outfile]
	if kwargs['cs']:
		cmd = "my_pairs_miscreateanchors %s %s" % (fname,outfile)
		os.system(cmd)
		c_outfile = fname+"_t"+str(int(thr))+"_cpc.gdf"
		mvCmd = "mv %s %s" % (outfile+".C.pairs.anc",c_outfile)
		os.system(mvCmd)
		anchorFiles.append(c_outfile)
		rmFiles.append(outfile+".heu")
	for f in anchorFiles: ## writing the output files
		f_out = subExt(f,newExt='anc')
		DialignAnchors.makeAnchorFile(f,f_out)
	rmFiles.extend(anchorFiles)
	for ff in rmFiles:   ## cleaning up
		rmCmd = "rm %s" % ff
		os.system(rmCmd)

def processArgs():
	""" Parser function of main """
	parser = OptionParser()
	parser.add_option("-c", "--consistent", action="store_true", dest="cs",help="Computes consistent columns",default=False)
	parser.add_option("-w", "--weight", action="store_true", dest="wg",help="Uses weighted incidence graph",default=False)
	parser.add_option("-b", "--boost", dest="bs",help="Uses Boost on subgraph/graph")
	parser.add_option("-f", "--filter", dest="filt",help="type of filtering (degree/edge weights)",default='degree')
	parser.add_option("-v", "--view", dest="vw",help="Views progress of edge removal (subgraph/boundary)")	
	parser.add_option("-t", "--threshold", dest="thr",help="Threshold value -- default 0",default=0)
	parser.add_option("-i", "--input", dest="frg",help="Optional alternative fragment input file")
	parser.add_option("-I", "--gdf-input", dest="gdf",help="Optional alternative fragment input file -- in gdf format")
	parser.add_option("-e", "--environment", action="store_true", dest="env",help="When used with gdf, adds the environment \
	to the anchor point",default=False)
	return parser

#========= Main program

if __name__ == '__main__':
	prog = sys.argv[0]
	parser = processArgs()
	options, args = parser.parse_args()
	usage = "usage: %s [options] fastaFile; Try -h for details" % prog
	if len(args) < 1:
		parser.error(usage)
	OPTIONS = {}
	for opt, value in options.__dict__.items():
		print opt,value
		OPTIONS[opt] = value
	fname = args[0]
	Main(fname,**OPTIONS)	
