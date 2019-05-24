#!/usr/bin/python

############################################################
#
#            Python code for PeakMatch algorithm
#
############################################################

import os,sys,re,math

EPSILON     = 0.000001
PSEUDO_RAD  = 2
PSEUDO_COEF = 2.0
INTV_RAD    = 50
INTV_COEF   = 1.1

#################### functions ####################
# Malucelli et al.'s algorithm for computing MWNCM
def MWNCM(W): # W[i][j]: weight of edge {i,j}
    n = len(W)
    m = len(W[0])
    # labelling phase
    LN = [0] * m
    L = [[0]*m for i in range(n)]
    for i in range(n):
        # Step 1
        for j in range(m):
            if j==0:
                L[i][j] = W[i][j]
            else:
                L[i][j] = W[i][j] + max(LN[:j])
        # Step 2
        for j in range(m):
            LN[j] = max(LN[j], L[i][j])
    A = []
    for i in range(n):
        for j in range(m):
            A.append((i,j,L[i][j]))
    A.sort(key=lambda a: -a[2])
    # selection phase
    E = []
    i_min = n
    j_min = m
    for a in A:
        if a[0]<i_min and a[1]<j_min:
            i_min = a[0]
            j_min = a[1]
            E.append((a[0],a[1],W[a[0]][a[1]]))
    return E

# compute average of numbers in array
def getAvg(array):
    return float(sum(array)) / float(len(array))

# compute variance of numbers in array
def getVar(array):
    avg = getAvg(array)
    A = [float(a-avg)*float(a-avg) for a in array]
    return sum(A)/float(len(A))

# compute standard deviation of numbers in array
def getStdev(array):
    return math.sqrt(getVar(array))

# compute exponential moving average of numbers in array; rad=(radius), coef=(smoothing factor)
def getMovAvg(array, rad, coef):
    movavg = []
    n = len(array)
    for i in range(n):
        num = 0.0
        den = 0.0
        for k in range(-rad,rad+1):
            if i+k<0 or i+k>=n:
                continue
            w = float(pow(coef, abs(k)))
            num += array[i+k]/w
            den += 1.0/w
        movavg.append(num/den)
    return movavg

# compute array [a[0]/b[0], a[1]/b[1], ...]
def getRatioArray(a, b):
    if len(a)!=len(b):
        sys.stderr.write('error: the sizes of a and b should be equivalent.\n')
        sys.exit(1)
    array = []
    for (x,y) in zip(a,b):
        if y < EPSILON:
            array.append(-1)
        else:
            array.append(float(x)/float(y))
    return array

# compute linear interpolation of x; inter=(# of inserted points)
def interpolate(x, inter):
    if inter==0:
        return x
    if len(x)==0:
        sys.stderr.write('error: interpolate is invoked although x=[]\n')
        sys.exit(1)
    y = []
    for i in range(len(x)-1):
        y.append(x[i])
        for k in range(1, inter+1):
            a = ((inter+1-k)*x[i] + k*x[i+1]) / float(inter+1)
            y.append(a)
    y.append(x[-1])
    return y

#################### classes ####################
# representation of interval
class Interval:
    def __init__(self, left, right, max, r=0):
        self.left = left
        self.right = right
        self.max = max
        self.r = r

# single/bulk data storage
class Data:
    def __init__(self, filename, fr=0, to=-1):
        if not os.path.isfile(filename):
            sys.stderr.write(filename+' is not found.\n')
            sys.exit(1)
        self.rowlabel = []
        self.collabel = None
        self.src = [] # src[k][t]=(gene k's expression level at pseudo time t)
        f = open(filename)
        for (t,row) in enumerate(f):
            a = re.split('\t', row)
            if t==0: # in head row, read column labels
                if to == -1:
                    self.collabel = a[fr+1:]
                else:
                    self.collabel = a[fr+1:to+2]
                self.cols = len(self.collabel)
                self.src = [[] for j in range(self.cols)]
                continue
            self.rowlabel.append(a[0])
            if to == -1:
                v = a[fr+1:]
            else:
                v = a[fr+1:to+2]
            if len(v) != self.cols:
                sys.stderr.write('error: len(v)='+str(len(v))+', but self.cols='+str(self.cols)+'\n')
                sys.stderr.write('       in '+filename+', the line '+str(t)+' is illegal.\n')
                print(v)
                sys.exit(1)
            for k in range(len(v)):
                if v[k] == '':
                    sys.stderr.write('warning: line '+str(t)+' contains an empty string. It is interpreted as zero.\n')
                    self.src[k].append(0.0)
                else:
                    self.src[k].append(float(v[k]))
        f.close()

    # compute exp moving average; the result is substituted to self.val
    def getMovAvg(self, rad=PSEUDO_RAD, coef=PSEUDO_COEF):
        self.val = []
        for j in range(self.cols):
            self.val.append(getMovAvg(self.src[j], rad, coef))

    # compute peak intervals
    def getPeaks(self, val, rad=INTV_RAD, coef=INTV_COEF, theta=1.0):
        self.I = []
        for j in range(len(val)):
            mval = getMovAvg(val[j], rad, coef) # exp moving average of gene j
            diff = [a-b for (a,b) in zip(val[j], mval)] # diff between actual value and moving average for gene j
            stdev = getStdev(diff) # diff's standard deviation
            I = []
            in_flag = False
            for k in range(len(val[j])):
                if val[j][k] >= mval[k]+theta*stdev:
                    if in_flag == False:
                        interval = Interval(left=k, right=k, max=k)
                        I.append(interval)
                        in_flag = True
                    I[-1].right = k
                    if val[j][k] > val[j][I[-1].max]:
                        I[-1].max = k
                else:
                    in_flag = False
            self.I.append(I)

    # compute linear interpolation; the result is substituted to self.val
    def interpolate(self, inter):
        self.val = []
        for j in range(self.cols):
            self.val.append(interpolate(self.src[j], inter))

    # row is a property
    def _get_rows(self):
        return len(self.val[0])

    def _set_rows(self):
        sys.stderr.write('warning: rows is read-only property.\n')

    def _del_rows(self):
        pass

    rows = property(_get_rows, _set_rows, _del_rows)


#################### main ####################

# arguments
if len(sys.argv)<7:
    sys.stdout.write('''

=== PeakMatch algorithm ===

The program outputs:
(pseudo time ID) (pseudo time name) (estimated real time) (delay from last time)


usage: main.py (single)(bulk)(T)(last)(intv)(inter)[(from)(to)]

  [mandatory]
  - single ... tabu-separated file for single cell
  - bulk   ... tabu-separated file for bulk
  - T      ... threshold for deciding peak; it is ratio. >=1.0 is recommended.
  - last   ... whether last times are forcibly matched (1) or not (0)
  - intv   ... whether all points in an interval are regarded as peak points (1) or not (0; i.e., only max is considered)
  - inter  ... number of points that are added between real time points for interpolation (e.g., 7 for 30min) 

  [optional]
  - from   ... ID of first gene in the target range; 0 is the leftmost gene.
  - to     ... ID of last gene in the target range; -1 is the rightmost gene.

''')    
    sys.exit(1)

if len(sys.argv)>=9:
    Single = Data(sys.argv[1], fr=int(sys.argv[7]), to=int(sys.argv[8]))
    Bulk = Data(sys.argv[2], fr=int(sys.argv[7]), to=int(sys.argv[8]))
else:
    Single = Data(sys.argv[1])
    Bulk = Data(sys.argv[2])

T = float(sys.argv[3])
(Last, Intv, Inter) = map(int, sys.argv[4:7])

# check the format
if Single.cols != Bulk.cols:
    sys.stderr.write('error: the numbers of columns are not equivalent.\n')
    sys.exit(1)

# compute exp moving average for single data
Single.getMovAvg()      # by this, Single.val becomes the moving average

# compute linear interpolation for bulk data
Bulk.interpolate(Inter) # by this, new points are inserted to Bulk.val

# compute peak intervals
Single.getPeaks(Single.val, theta=T) 
Bulk.getPeaks(Bulk.val, theta=T)

# assign weights to complete bipartite graph K_{n,m}
n = Single.rows 
m = Bulk.rows
W = [[0]*m for i in range(n)]
W[0][0] = Bulk.rows+1

if Last:
    W[Single.rows-1][Bulk.rows-1] = Bulk.rows+1

if Intv==0:
    for c in range(Single.cols):
        comb = len(Single.I[c]) * len(Bulk.I[c])
        for s_intv in Single.I[c]:
            i = s_intv.max
            for b_intv in Bulk.I[c]:
                j = b_intv.max
                W[i][j] += 1.0 / float(comb)
else:
    for c in range(Single.cols):
        comb = 0
        for s_intv in Single.I[c]:
            for b_intv in Bulk.I[c]:
                comb += (s_intv.right-s_intv.left+1) * (b_intv.right-b_intv.left+1)
        for s_intv in Single.I[c]:
            for b_intv in Bulk.I[c]:
                for i in range(s_intv.left, s_intv.right+1):
                    for j in range(b_intv.left, b_intv.right+1):
                        W[i][j] += 1.0 / float(comb)


# compute maximum weighted non-crossing matching
M = MWNCM(W)
M.sort(key=lambda e: e[0]) # M is ordered w.r.t. pseudo times
St = [e[0] for e in M]     # St is the set of matched pseudo times
I = [St[k+1]-St[k] for k in range(len(M)-1)] 

# output the result
k = -1
A = []
for i in range(Single.rows):
    if i in St:
        k += 1
        A.append(M[k][1])
    else:
        x = 0
        if k < len(St)-1:
            x = float( (i-St[k])*M[k+1][1] + (St[k+1]-i)*M[k][1] ) / float(St[k+1]-St[k])
        else:
            x = A[i-1] + EPSILON
        A.append(x)
    if i>0:
        print(i,Single.rowlabel[i],A[i],A[i]-A[i-1])
    else:
        print(i,Single.rowlabel[i],A[i])

