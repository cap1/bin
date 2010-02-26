#! /usr/bin/env python

# -*- coding: iso-8859-1 -*-


import sys,getopt,os,re,tG
from subprocess import *
from altgraph import Graph
import networkx as nx
# la commande magique: for f in $(ls *.llr); do gawk '{if ($2 > 2.0); print $0}' $f > $f.filter2.0; done

def myCompare(x,y):
    if x < y:
        return -1
    elif x > y:
        return 1
    else:
        return 0

def lexicoTuple(x,y):
    """Defines the lexicographic ordering of couples (seq,pos)"""
    try:
        x0,x1 = x
        y0,y1 = y
    except (ValueError,TypeError):
        return myCompare(x,y)
    a = lexicoTuple(x0,y0)
    if a:
        return a
    else:
        return lexicoTuple(x1,y1)

def tupleToDict(siteList):
    """For a lexicoTuple ordered list of sites (tuples), convert it into the dictionary seq -> ordered list of positions in seq"""
    seqSite = {}
    for site in siteList:
        s,p = site
        try:
            seqSite[s].append(p)
        except KeyError:
            seqSite[s] = [p]
    return seqSite

def extractTuple(a,sep='()'):
    return map(lambda x:int(x),a.strip(sep).split(', '))

def extractTuple0(a,sep='()'):
    return a.strip(sep).split(', ')

def tuplify(x):
    """A syntactic routine to tackle a format problem"""
    y = extractTuple(x)
    return (y[0],y[1])

def parseTupleList(line):
	""" Convert a string formatted as a Python list into the Python list itself -- a pickle would be better... """
	l = line.strip('[]').split('), (')
	outList = map(tuplify,l)
	return outList
    
def matchTuple(a):
    p = re.compile('\(\d+, \d+\)')
    return map(extractTuple,p.findall(a))

def matchNb(s1,s2):
    """ Counts the nb of matches between strings of same length """
    if len(s1) != len(s2):
        return None
    else:
        Z = zip(list(s1),list(s2))
        M = map(lambda x:x[0]==x[1],Z)
        return reduce(lambda x,y:x+y,M)

def myDiff(itemList):
    resList = []
    for i in range(len(itemList))[1:]:
        for j in range(i):
            resList.append(itemList[i]-itemList[j])
    return resList

def normList(inputString):
    """Routine to convert eg 1-3 5 9-11 6
    into the list [1,2,3,5,6,9,10,11]"""
    endSet = set([])
    for t in inputString.split():
        j = t.split("-")
        j = map(lambda x:int(x),j)
        j.sort()
        jList = range(j[-1]+1)[j[0]:]
        for i in jList:
            endSet.add(i)
    endList = list(endSet)
    endList.sort()
    return endList

def sortList(itemList,function=None):
    """ Goes around a funny behaviour of Python's sort -- when list elements are lists or tuples """
    itemDict = dict([(i,itemList[i]) for i in range(len(itemList))])
    functionDict = {}
    for i in range(len(itemList)):
        #print itemList[i]
        functionDict[i] = function(itemList[i])
    #print functionDict
    vDict = tG.invertDico1(functionDict)
    #print vDict
    valueDict = dict([(k,[itemDict[j] for j in vDict[k]]) for k in vDict.keys()])
    #print valueDict
    ks = valueDict.keys()
    #print ks
    ks.sort()
    outList = []
    for i in ks:
        for L in valueDict[i]:
            outList.append(L)
    return outList

def myZip(l1,l2):
    l = []
    for i in range(len(l1)):
        try:
            ll = l1[i]
            ll.append(l2[i])
        except AttributeError:
            ll = [l1[i]]
            ll.append(l2[i])
        l.append(ll)
    return l

def executeCmd(*args):
	cmd = "%s"
	for i in range(len(args)-1):
		cmd += " %s"
	return cmd % args
	
def subExt(filename,newExt):
	""" Substitute newExt to current filename extension """
	path,name = os.path.split(os.path.abspath(filename))
	lst = name.split(".")[0:-1]
	string = ''
	for l in lst:
		string += l+"."
	string += newExt
	return os.path.join(path,string)

def fileBase(filename):
	""" Get name without extension, """
	path,name = os.path.split(os.path.abspath(filename))
	lst = name.split(".")[0:-1]
	string = ''
	for l in lst:
		string += l+"."
	string = string.strip(".")
	return os.path.join(path,string)

def fileSplit(filename):
	""" Get name without extension, and the extension """
	path,name = os.path.split(filename)
        lst = name.split(".")
        ext = lst.pop()
	string = ''
	for l in lst:
		string += l+"."
	string = string.strip(".")
	return os.path.join(path,string),ext

def typeDict(filename):
	""" Exclusively for anchorFile format Fasta_dectype_sminX_anctype/csType.ext files """
	dectype = {'ex':'exact','mm':'mismatch','fg':'fragments'}
        anctype = {'b':'b','bw':'bw','r':'raw','rw':'raw','flw':'consistent','csDAG':'consistent','csDFS':'consistent','fl':'consistent'}
        cstype = {'b':None,'bw':None,'r':None,'rw':None,'csDAG':'DAG','csDFS':'DFS','fl':'fl','flw':'fl'}
	path,fname = os.path.split(filename)   # extract file name
	f = fname.split(".")[0]	               # extract file base name
	chars = f.split("_")[1:]                   # extract anchor characteristics
	outDict = {'dectype':dectype[chars[0]],'anctype':anctype[chars[2]],'cstype':cstype[chars[2]],'smin':int(chars[1].strip('smin'))}
	return outDict

def removeAll(item='', item_list=[]):
    """ This function removes all occurrences of 'item' in 'item_list' """
    while True:
        try:
            item_list.remove(item)
        except ValueError:
            break
    return item_list

def subDivide(itemList):
	""" Replaces an integer list by a list of lists of successive integers """
	cList,sList = [],[]
	while itemList:
		cv = itemList.pop(0)
		try: 
			v = cv == start+1
		except NameError:
			v = True
		if v:
			cList.append(cv)
			start = cv
		else:
			sList.append(cList)
			cList = [cv]
			start = cv
	sList.append(cList)
	return sList

def suffixFile(anchorFile,suffix="0"):
    """ another joke for pairwiseAnchors """
    f = open(anchorFile,'U')
    g = open("tmp.anc",'w')
    for line in f.xreadlines():
        l = line.strip().split(' ')
        newline = ''
        for s in l:
            try:
                int(s)
                s = str(s)+"_0"
            except ValueError:
                pass
            newline += s+' '
        print >> g,newline
    f.close()
    g.close()
    cmd = "mv tmp.anc %s" % anchorFile
    os.system(cmd)

def stringify(itemList):
    s = ''
    while True:
        i = itemList.pop(0)
        if itemList:
            s += str(i)+','
        else:
            s += str(i)
            break
    return s

def multiNomial(itemList,itemProbas):
    """ Given an itemList and the elementary probabilities, return the probability of the input vector -- uses an Rscript file """
    s1 = stringify(itemList)
    s2 = stringify(itemProbas)
    f = open('tmp.r','w')
    print >>f,"#! /usr/bin/Rscript"
    str1 = "x <- c(%s)" % s1
    str2 = "p <- c(%s)" % s2
    print >>f, str1
    print >>f, str2
    print >>f,"dmultinom(x, prob = p)"
    f.close()
    cmd1 = "chmod +x tmp.r"
    cmd2 = "./tmp.r"    
    os.system(cmd1)
    return float(runAndRead(cmd2)[1])

def chiNorm(refDist,testDist):
    """ Computes the chi-sqared test of goodness-of-fit for the testDistribution wrt the refDistribution """
    s = sum(refDist.values())
    refList,testList = [],[]
    for i in refDist.keys():
        refList.append(refDist[i])
        try:
            testList.append(testDist[i])
        except KeyError:
            testList.append(0)
    s1 = stringify(refList)
    s2 = stringify(testList)
    f = open('tmp.r','w')
    print >>f,"#! /usr/bin/Rscript"
    str1 = "pr <- c(%s)/%s" % (s1,str(s))
    str2 = "x <- c(%s)" % s2
    print >>f, str1
    print >>f, str2
    print >>f,"chisq.test(x, p = pr)"
    f.close()
    cmd1 = "chmod +x tmp.r"
    cmd2 = "./tmp.r"    
    os.system(cmd1)
    try:
        lines = runAndReadAll(cmd2)
        return float(lines[4].split(",")[2].strip().strip('\n').split()[2].strip())
    except ValueError:
        return 1.5

def countList(itemList):
    """ Replaces repeated item list by list of tuples (item,counts) """
    newList = []
    tmpList = itemList[:]
    tmpList.sort()
    start = tmpList.pop(0)
    k = 1
    while tmpList:
        try:
            nx = tmpList.pop(0)
            if nx == start:
                k += 1
            else:
                newItem = (start,k)
                newList.append(newItem)
                start = nx
                k = 1
        except IndexError:
            newItem = (start,k)
            newList.append(newItem)
    newItem = (start,k)
    newList.append(newItem)
    return newList

def invertDico(dico):
    """
    Builds the reverse dictionary value -> (list of) keys -- not tested for velocity"""
    ocid = {}
    for i,v in dico.iteritems():
        try:
            ocid[v].append(i)
        except IndexError:
            ocid[v] = [i]
    return ocid

def invertDico1(dico):
    """ builds the reverse dictionary value -> (list of) keys """
    lst1 = [(v,i) for i,v in dico.iteritems()]
    print lst1
    gr = Graph.Graph()
    for e in lst1:
        gr.add_edge(*e)
    lst2 = [(e,gr.out_nbrs(e)) for e in gr.nodes if gr.out_nbrs(e)]
    ocid = dict(lst2)
    return ocid

def subDict(dico,condition):
    subDico = {}
    for key,value in dico.iteritems():
        if condition(key):
            subDico[key] = value
    return subDico

def copyDict(dico):
    """Performs a deep copy of dico"""
    g = {}
    for i,v in dico.iteritems():
        try:
            g[i] = v[:]
        except TypeError:
            g[i] = v
    return g

def partitionUnion(setList,newSet):
    """Computes the partition obtained by adding newSet (or list of newSets) to the setList"""
    #initial test
    if not setList and newSet:
        return [newSet]
    # from now on, we have a non-empty setList
    att = {}                   # this is the element -> index dictionary
    k = 0                       # initialise the first index to 0
    # start by attributing an index to each set in setList (ie the same to all of its elements)
    for s in setList:    
        for v in s:
            att[v] = k
        k += 1
    # the dictionary att is isomorphic to setList
    tta = tG.invertDico1(att)        # tta is the reverse dictionary index -> list of elements in set labelled by 'index'
    # here we compute the list of labels such that set indexed by label has a non-empty intersection with newSet
    labels = set([])
    for i in newSet:
        try:
            labels.add(att[i])
        except KeyError:
            pass
    # if no label is found, then the newSet is disjoint from sets in setList, and will be simply added
    newClass = newSet
    # if an intersection exists, we construct the reunion of sets with indices in labels
    for l in labels:
        nx = set(tta[l])
        newClass.update(nx)
        del tta[l]                 # and destroy the old labels
    # we relabel all elements in newClass with the new label k
    for j in newClass:
        att[j] = k
    # and construct the new description of sets tta
    tta = tG.invertDico1(att)
    newSets = []
    # we now use tta to generate the new setList
    for m in tta.keys():
        s = set(tta[m])
        newSets.append(s)
    return newSets
               
def partitionUnionList(setList,newSetList):
    """List of newSets version of partitionUnion"""
    if not isinstance(newSetList,list):
        newSetList = [newSetList]
    for newSet in newSetList:
        setList1 = partitionUnion(setList,newSet)
        setList = setList1
    return setList

def printDict(dico):
    for i,v in dico.iteritems():
        v.sort(lexicoTuple)
        print i,v

def concatList(sList,sep=''):
    """Convert list into string of concatenated elements"""
    line = ''
    for s in sList:
        line += sep+str(s)
    return line

#======== Matrix functions

def initMatrix(n):
        mm = [[0]*n for i in range(n)]
        return mm

def printMatrix(a,format='(%3d,%3d)'):
    """Prints a matrix, implemented as a list of lines"""
    l,c = len(a),len(a[0])
    # Formatage de la matrice
    form = '"'+format
    for i in range(c-1):
        form += ' '+format
    form += '"'
    for j in range(l):
        line = ()
        for i in range(c):
            line += (a[j][i],)
        print form % line    

def writeNexus(filename,distances,names,format='%.3f ',comment='',ext='.nexus'):
    """Fonction qui ecrit une matrice au format nexus (cf doc SplitsTree) -- For any complaints, cf. Claudine"""
    # ouverture fichier 
    fname = filename+ext               # rajout extension
    f = open(fname,'w')
    nm = {}
    # ecriture resultat
    
    # header nexus
    nbseq = len(distances)
    f.write('#NEXUS\n\n')
    titre = '[Distance matrix calculated by '+comment+']\n\n'
    f.write(titre)                                  # titre
    f.write('BEGIN taxa;\n')                # labels
    dimension = '\tDIMENSIONS ntax = %d;\n' % (nbseq) # - elt[0]
    f.write(dimension)
    f.write('TAXLABELS\n')
    for i in xrange(0,nbseq):
        nm[i] = names[i].replace(' ','_').replace('/','-')
        f.write('\t'+nm[i]+'\n')
    f.write(';\n')
    f.write('END;\n\n')
    f.write('BEGIN distances;\n')   # header distances
    f.write(dimension)
    f.write('\tFORMAT\n')
    f.write('\t\ttriangle=BOTH\n')
    f.write('\t\tdiagonal\n')
    f.write('\t\tlabels\n')
    f.write(';\n')
    f.write('MATRIX\n')
    
    # matrice de distances
    for i in xrange(0,nbseq):
        ligne = nm[i]+'\t'
        for d in distances[i]:  # parcours une ligne de la matrice (sauf elt0)
            w = format % (d)               # formatage du reel
            ligne = ligne + w       # rajout a la ligne en cours
            ligne = ligne+'\n'
        f.write(ligne)
    f.write(';\n')
    f.write('END;\n')
    f.close()

def msf2fasta(filein,fileout):
    f = open(filein,'U')
    names = []
    fasta = {}
    for line in f.xreadlines():
        l = line.split(' ')
        if l[0] in names:
            try:
                fasta[l[0]].extend(l[1:])
            except:
                fasta[l[0]] = l[1:]
        if line.startswith(' Name:'):
            l = line.split(' ')
            names.append(l[2])
    #print names
    for i in names:
        string = ''
        for el in fasta[i]:
            for letter in el:
                if re.match('\w',letter):
                    string += letter.upper()
        if string:
            fasta[i] = string
        else:
            fasta[i] = string
            msg = 'Attention, sequence %s in %s was removed because it is empty\n' %(i,filein)
            sys.stderr.write(msg)
    #print fasta
    f.close()
    g = open(fileout,'w')
    for i in names:
        if fasta[i]:
            print >> g,'>'+i
            print >> g,fasta[i]
    g.close()
    sys.stderr.write('transformed '+filein+'...\n'+'...into '+fileout+'\n')	

def msf2anchors(filein,fileout,mode='exec'):
    """ Convert an msf file into a collection of anchor columns """
    f = open(filein,'U')
    names = []
    fasta = {}
    for line in f.xreadlines():
        l = line.split(' ')         # split read line into letter blocks
        if l[0] in names:
            try:
                fasta[l[0]].extend(l[1:])
            except:
                fasta[l[0]] = l[1:]   # here fasta[i] is a list of strings (letter blocks)
        if line.startswith(' Name:'):  
            l = line.split(' ')
            names.append(l[2])
    #print names
    newNames = []
    for i in names:
        k = 0                         # index of character in sequence
        gappedSequence = []           # here we stack the line of msf as [(i,0),(i,1),'.',(i,2)....]
        for el in fasta[i]:           # for each letter block
            for letter in el:         # for each element in the block
                if re.match('\w',letter):  # if element is a letter
                    #print "Character %s is a letter " % letter
                    char = k
                    k += 1            # increment letter index
                elif letter == '.':        # if it is a '.'
                    char = letter
                try:
                    gappedSequence.append(char)
                except:
                    pass
        if k:
            fasta[i] = gappedSequence
            newNames.append(i)
        else:
            del fasta[i]
            msg = 'Attention, sequence %s in %s was removed because it is empty\n' %(i,filein)
            sys.stderr.write(msg)
    names = newNames[:]
    for ii in range(len(names)):
        fasta[ii] = fasta[names[ii]]
        del fasta[names[ii]]                  # here we replace keys of fasta by integers
    m = len(fasta[0])                 # normally all values of fasta dictionary have same length (the length of the multiple alignment)
    anchors = {}
    for j in range(m):                # for each column in the alignment
        anchor = []
        for iii in range(len(names)):   # list all positions of the sequences that have been aligned in this column
            char = fasta[iii][j]
            if re.match('\d+',str(char)):
                site = (iii,char)
                anchor.append(site)
                #print "Added character %3d as site (%3d,%3d)" % (char,iii,char)
            #else:
                #print "Not added character %s " % char
        anchors[j] = anchor        # so far we are building the list of the equivalence classes corresponding to the msa
    if mode == 'exec':
        g = open(fileout,'w')
        for j in range(m):
            anchorString = ''
            for site in anchors[j]:
                anchorString += str(site)+' '
            print >> g,anchorString+'\n'
        g.close()
    sys.stderr.write('transformed '+filein+'...\n'+'...into '+fileout+'\n')

def runAndReadAll(cmd,lineIndex=-1):
	""" Will be deprecated as soon as I understand the new subprocess Module """
	p = Popen([cmd], shell=True,stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
	(child_std_in, child_stdout_and_err) = (p.stdin, p.stdout)
	ret = p.wait()
	#print "Returned %s code" % ret
	#sout,sin = popen2.popen2(cmd)
	if ret:
		raise ValueError
		return
	return child_stdout_and_err.readlines()

def parseBlosum(blosumFile):
    pairValue = {}
    f = open(blosumFile)
    letters = f.readline().strip().split()
    for line in f.xreadlines():
        line = line.strip('\n').strip()
        if line:
            values = line.strip().split()
            letter = values.pop()
            for i in range(len(values)):
                key = frozenset([letter,letters[i]])
                pairValue[key] = int(values[i])
            letters.pop(0)  # discard current letter (tr sup)
    return pairValue

def blosumScore(s1,s2,pairValue):
    if len(s1) != len(s2):
        return
    score = 0
    ss1,ss2 = list(s1),list(s2)
    spair = []
    while ss1:
        lpair = frozenset([ss1.pop(),ss2.pop()])
        spair.append(lpair)
    for i in spair:
        try:
            score += pairValue[i]
        except KeyError:
            print (s1,s2)
    return score

def main():
    return

#
# Programme principal
#
#

if __name__  ==  '__main__':
        
    #========
    #
    #                Lecture des arguments
    #
    #========
    
    m = len(sys.argv)                            # lecture des arguments passes en ligne de commande (nombre d'iceux)
    if m >= 1:                                      # mode script : les arguments ont ete passes en ligne de commande
        main(*sys.argv[1:])
    else:                                                  # nombre incorrect d'arguments : fin du script
        sys.exit('Usage : mySyntax.py')


