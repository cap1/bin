# -*- coding: iso-8859-1 -*-

import sys,getopt,os,time,re

##
# Fasta
#

class FastaObject:

    def __init__(self,fname=None,sep=''):
        if fname:
            self.sequences, self.sites = readFastaSites(fname,sep=sep)
            self.code = dict([(seq,self.sequences.index(seq)) for seq in self.sequences])
            self.site = self.positionsToSites()
            self.position = self.sitesToPositions()
        else:
            self.sequences = []
            self.sites = {}
            self.code = {}
            self.site = {}
            self.position = {}
        self.filename = fname

     #=========  Constructor methods

    def siteIterator(self):
        """Returns an iterable on (seq,pos) with seq repeated as many times as sites in the sequence"""
        sit = {}
        ss = []
        dec = 0
        s = len(self.sequences)
        for i in range(s):
            sit[i] = len(self.sites[i])+dec
            a = [i for t in range(sit[i])]
            ss.extend(zip(a,range(sit[i])))
        return iter(ss)

    def positionsToSites(self):
        """Constructs the dictionary pos -> (seq,position) for self  """
        sites = [s for s in self.siteIterator()]
        site = dict(zip(range(len(sites)),sites))
        position = dict(zip(sites,range(len(sites))))
        return site

    def sitesToPositions(self):
        """Constructs the dictionary (seq,position) -> pos for self """
        sites = [s for s in self.siteIterator()]
        position = dict(zip(sites,range(len(sites))))
        return position    

    #====== Attribute readers

    def getSequences(self):
        return self.sequences

    def getSites(self):
        return self.sites

    #======= Syntaxic methods

    def copyFastaObject(self):
        """Copy method"""
        g = FastaObject()
        g.sequences = self.sequences[:]
        g.sites = copyDict(self.sites)
        g.code = copyDict(self.code)
        g.site = copyDict(self.site)
        g.position = copyDict(self.position)
        g.filename = self.filename[:]
        return g

    def recoverFasta(self,name=None,ext=None,sep=' ',lineLength=80):
        """Reconstructs a fasta file from the FastaObject
        -- dir = directory, ext = extension, to recognize it from the input, sep=to separate characters """
        ext = "."+ext
        if name == None: 
            name = self.filename
        fname = str(name)+str(ext)
        f = open(fname,'w')
        for seq in self.sequences:
            line1 = '>%s\n'% (seq)   
            f.write(line1)                    
            lineList = self.sites[self.code[seq]]  
            i = 0                                                         
            newline = ''                                                         
            for l in lineList:                                                
                if l:                                                                  
                    newline += str(l)+sep                               
                    i += 1
                    try:
                        if i % int(lineLength) == 0:                                           
                            newline += '\n'                                    
                            f.write(newline)                                   
                            newline = ''
                    except ZeroDivisionError:
                        pass
            f.write(newline+'\n')                                     
        f.close()

def readFastaSites(name,sep=''):
    """This the main function of FastaObject constructor!!
    Reads a FASTA or enriched FASTA file (separator option) 
    Returns
    1. the list of sequence names
    2. the dictionary index of name -> list of sites of sequence"""	
    sites = {}
    f = open(name,'U')
    s = []
    seqs = []
    for line in f.xreadlines():
        line = line.strip().strip('\n')
        if line.startswith('>'):
            if s:
                sites[sq] = s[:]
            sq = line.strip('>')
            s = []
            seqs.append(sq)
        elif line:
            if sep == '':
                s.extend(list(line))
            elif sep == ' ':
                s.extend(line.split(sep))
            else:
                print "Not a valid Fasta file"
                sys.exit(1)
    sites[sq] = s[:]
    f.close()
    newsites = {}
    for a in sites.keys():
        newsites[seqs.index(a)] = sites[a]
    return seqs,newsites

##
# Anchors
#

class Anchors(FastaObject):

    def __init__(self,filename=None,anchorDict=None,sep=''):
        """ In this new way, we can cast a FastaObject into an Anchors object --
        filename = FastaFile
        anchorDict = dictionary of anchors """
        if filename:
            FastaObject.__init__(self,filename,sep=sep)
        else:
            FastaObject.__init__(self)
        if type(anchorDict) == dict:
            anchorDict.values().sort(lexicoTuple)
            self.anchors = anchorDict
            self.classes = self.getClasses()
        elif anchorDict == 'suffix':
            self.anchors = self.letterLocation()
            self.classes = self.getClasses()
        self.filename = filename
        
    #======= Constructor methods

    def getClasses(self):
        """ Computes the inverse dictionary of anchors ie site -> name of anchor """
        classes = {}
        for node,values in self.anchors.iteritems():
            if len(values)>1:
                for i in values:
                    classes[i] = node
        return classes

    def letterLocation(self):
        """Forms the dictionary letter -> list of sites (seq,position) bearing letter"""
        classes = {}
        for s,l in self.sites.iteritems():
            for i in range(len(l)):
                try:
                    classes[l[i]].append((s,i))
                except KeyError:
                    classes[l[i]] = [(s,i)]        
        for k in classes.keys()[:]:
            if not re.search('\d+_\d',k):
                del classes[k]
        return classes
        
    #======== Syntaxic methods

    def recoverAnchors(self,outbase=None,outext='gdf',lineLength=80):
        """ Rewrites a fasta GD-style file with the identifiers of anchors as letters """
        print self.filename,outext
        ff = FastaObject(self.filename)
        f = ff.copyFastaObject()
        for n in self.anchors.keys():
            for site in self.anchors[n]:
                seq,pos = site[0],site[1]
                f.sites[seq][pos] = n
        f.recoverFasta(name=outbase,ext=outext,lineLength=lineLength)  # Outputs the file of partition
        print "Wrote file",outbase+"."+outext,"by recoverAnchors"

    def outputAnchorPairs(self):
        """ Creates pairwise segments of anchors """
        pairList = []
        initialAnchorList = self.anchors.values()
        seqSet = set([])
        for a in initialAnchorList:
            for s in a:
                seqSet.add(s[0])
        seqList = list(seqSet)
        seqList.sort()  # up to here, we only got the sequence indices
        anchorDict = {}
        for i in seqList: # create anchorDict seq: anchors starting with seq (when sorted)
            anchorDict[i] = [l for l in initialAnchorList if l[0][0] == i]
            anchorDict[i].sort()
        for i in seqList:  # for each sequence do create anchors and update anchorDict
            anchorList = anchorDict[i]
            if anchorList:
                newAnchorList,segmentList = self.createSegments(anchorList)
                del anchorDict[i]
                anchorDict = self.updateAnchorDict(anchorDict,newAnchorList)
                pairList.extend(segmentList)
        pairList.sort()
        return pairList
            
    def updateAnchorDict(self,anchorDict,newAnchorList):
        """ Auxiliary function of outputAnchorPairs -- updates anchorDict by redistributing
        truncated anchors generated by createSegments from newAnchorList
        into anchorDict according to their (new) first element """
        for l in newAnchorList:
            anchorDict[l[0][0]].append(l)
        for i in anchorDict.keys():
            anchorDict[i].sort()
        return anchorDict

    def createSegments(self,anchorList):
        """ For a list of homogeneous anchors (starting with a site in the same sequence i), sorted according to the 
        sequence order of sequence i, outputs the set of pairwise segments involving sequence i only, and returns
        the unordered list of remaining anchors, ie the same anchors without the site in sequence i, and the set of
        segments """
        segmentList = []
        newAnchorList = [] # storage for anchor list without first sequence i
        startList = anchorList.pop(0)
        i0 = startList[0][0]  # index of the current start sequence
        if len(startList) > 2: # if the anchor minus the first site is still an anchor (ie at least 2 elements + the first)
            newAnchorList.append(startList[1:])  # keep anchor minus first site for future reference
        start = dict(startList)
        length = dict([(i,1) for i in start.keys()])
        current = copyDict(start)
        while anchorList:
            nextList = anchorList.pop(0)
            if len(nextList) > 2:
                newAnchorList.append(nextList[1:]) # keep anchor minus first site for future reference
            next = dict(nextList)
            nextSet = set(next.keys())
            currentSet = set(current.keys())
            stack = list(nextSet.intersection(currentSet))  
            stack.sort()
            j0 = stack.pop(0) #this one is necessarily equal to i
            i0 = j0
            d0 = next[j0]-current[j0]
            if d0 > 1: # start all afresh because segment is discontinued in reference sequence i
                #output all pairwise current segments (i,x) present in current
                segmentList.extend([((j0,start[j0]+length[j0]-length[j]),(j,start[j]),length[j]) for j in current.keys() if j != j0]) 
                start = copyDict(next)
                length = dict([(i,1) for i in start.keys()])
                current = copyDict(start)
            else: # if segment on sequence i goes on
                length[j0] += 1
                currentSet.remove(j0)
                nextSet.remove(j0)
                while stack:  # for all subsequent sequence indices
                    j = stack.pop(0)
                    d = next[j]-current[j]
                    if d == 1:   # segment for j continues
                        length[j] += 1
                    else:   # reinitialise data for sequence j only
                        segmentList.append(((j0,start[j0]+length[j0]-length[j]-1),(j,start[j]),length[j])) # output segment (i,j)
                        start[j] = next[j]
                        length[j] = 1
                    currentSet.remove(j)
                    nextSet.remove(j)
                currentList = list(currentSet) # this contains seq indices of current that are not in next
                nextList = list(nextSet)   # this contains seq indices of next that are not in current 
                #print currentList,nextList
                # for keys only in current (fragments that stop there on sequence j != i)
                while currentList: # these are only present in current, therefore segment stops
                    j = currentList.pop(0)
                    segmentList.append(((j0,start[j0]+length[j0]-length[j]-1),(j,start[j]),length[j])) # output segment (i,j)
                    del start[j] # no current segment on j remains
                    del length[j]
                while nextList: # these are only present in next, therefore new segment starts
                    j = nextList.pop(0)
                    start[j] = next[j]  # create key j in start 
                    length[j] = 1       # and start segment length index for j
            current = copyDict(next)
        segmentList.extend([((i0,start[i0]+length[i0]-length[j]),(j,start[j]),length[j]) for j in current.keys() if j != i0])
        newAnchorList.sort()
        return (newAnchorList,segmentList)

## == Other syntaxic functions

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

def copyDict(dico):
    """Performs a deep copy of dico"""
    g = {}
    for i,v in dico.iteritems():
        try:
            g[i] = v[:]
        except TypeError:
            g[i] = v
    return g

def makeAnchorFile(fname,outName):
	a = Anchors(filename=fname,anchorDict='suffix',sep=' ')
	pairList = a.outputAnchorPairs()
        g = open(outName,'w')
        for anc in pairList:
            s1,s2,l = anc
            seq1,seq2 = s1[0]+1,s2[0]+1
            pos1,pos2 = s1[1]+1,s2[1]+1
            length = l
            weight = l
            line = "%d %d %d %d %d %d" % (seq1,seq2,pos1,pos2,length,weight)
            print >> g, line
        g.close()
        print "Created file %s" % outName
