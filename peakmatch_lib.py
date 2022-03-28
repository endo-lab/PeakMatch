import os,sys,re,math

EPSILON     = 0.000001
PSEUDO_RAD  = 2
PSEUDO_COEF = 2.0
INTV_RAD    = 50
INTV_COEF   = 1.1

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
    def contains(self, intv, i):
        if intv==0:
            if i==self.max:
                return True
        else:
            if i>=self.left and i<=self.right:
                return True
        return False
                
            

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
                #print(v)
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
        self.inter_rowlabel = []
        for label in self.rowlabel:
            self.inter_rowlabel.append(label)
            if label != self.rowlabel[-1]:
                for i in range(1,inter+1):
                    self.inter_rowlabel.append(f"{label}_{i}")
        #print(self.inter_rowlabel)
        #print(len(self.inter_rowlabel))
        #print(len(self.val[0]))
        #print(len(self.val[1]))
              

    # row is a property
    def _get_rows(self):
        return len(self.val[0])

    def _set_rows(self):
        sys.stderr.write('warning: rows is read-only property.\n')

    def _del_rows(self):
        pass

    rows = property(_get_rows, _set_rows, _del_rows)

