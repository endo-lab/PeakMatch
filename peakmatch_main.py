#!/usr/bin/python
############################################################
#
#   Python code for PeakMatch algorithm (v04: March 2022)
#
############################################################

import os,sys,time,random
from peakmatch_lib import *
        
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
    opt = 0
    opt_all = 0
    for a in A:
        if a[0]<i_min and a[1]<j_min:
            i_min = a[0]
            j_min = a[1]
            E.append((a[0],a[1],W[a[0]][a[1]]))
            if W[a[0]][a[1]] < 1.00001:
                opt += W[a[0]][a[1]]
            opt_all += W[a[0]][a[1]]
    return E,opt,opt_all,0

# conventional maximum weighted matching using networkx
def MWM(W):
    try:
        import networkx as nx
    except:
        sys.stderr.write("error: failed to import networkx; confirm whether networkx is installed on your system.\n")
        exit(1)
    n,m = len(W),len(W[0])    
    G = nx.Graph()
    G.add_nodes_from([f"single{i}" for i in range(n)])
    G.add_nodes_from([f"bulk{j}" for j in range(m)])
    for i in range(n):
        for j in range(m):
            G.add_edge(f"single{i}", f"bulk{j}", weight=W[i][j])
    nxM_raw = nx.max_weight_matching(G)
    Att = nx.get_edge_attributes(G,"weight")
    opt = 0
    opt_all = 0
    nxM = []
    for e in nxM_raw:
        if not "single" in e[0]:
            e = (e[1],e[0])
        s = e[0].replace('single','')
        b = e[1].replace('bulk','')
        nxM.append((int(s), int(b), Att[e]))
        if Att[e] < 1.00001:
            opt += Att[e]
        opt_all += Att[e]
    cross = 0
    for e1 in nxM:
        for e2 in nxM:
            if e1[0] < e2[0] and e1[1] > e2[1]:
                cross += 1
    return nxM, opt, opt_all, cross

# generate a non-crossing matching randomly
def RAND(W):
    n,m = len(W),len(W[0])
    if n>=m:
        sys.stderr.write(f"error: the bulk side should contain nodes more than the single side; but {n}>={m}.\n")
        sys.exit(1)
    R = list(range(1,m-1))
    random.shuffle(R)
    B = [0] + R[:n-2] + [m-1]
    B.sort()
    M = []
    opt = 0
    opt_all = 0
    for k in range(n):
        w = W[k][B[k]]
        M.append((k,B[k],w))
        if w < 1.00001:
            opt += w
        opt_all += w
    return M,opt,opt_all,0
#################### main ####################
def main(argv):
    start_time = time.time()
    
    # read mandatory parameters
    T = float(argv[3])
    (Last, Intv, Inter) = map(int, argv[4:7])

    # read options
    pred_file = ''
    sum_file = ''
    from_col = 0
    to_col = -1
    relax = False
    rand = False
    for arg in argv[7:]:
        try:
            if "-pred=" in arg:
                arr = arg.split('=')
                pred_file = arr[1]
            elif "-sum=" in arg:
                arr = arg.split('=')
                sum_file = arr[1]
            elif "-from=" in arg:
                arr = arg.split('=')
                from_col = int(arr[1])
            elif "-to=" in arg:
                arr = arg.split('=')
                to_col = int(arr[1])
            elif "-relax" in arg:
                relax = not relax
            elif "-rand" in arg:
                arr = arg.split('=')
                random.seed(int(arr[1]))
                rand = not rand
        except:
            sys.stderr.write("error: {} cannot be read.\n".format(arg))
            sys.exit(1)
    
    if pred_file == '':
        pred_fp = None
    else:
        pred_fp = open(pred_file, 'w')
    if sum_file == '':
        sum_fp = None
    else:
        sum_fp = open(sum_file, 'w')

    activated = [1 if a==True else 0 for a in (relax,rand)]
    if sum(activated)>=2:
        sys.stderr.write("error: illegal option.\n")
        sys.exit(1)
        
    #read data
    Single = Data(argv[1], fr=from_col, to=to_col)
    Bulk = Data(argv[2], fr=from_col, to=to_col)

    # check the format
    if Single.cols != Bulk.cols:
        sys.stderr.write('error: the numbers of columns are not equivalent.\n')
        sys.exit(1)

    # compute exp moving average for single data
    Single.getMovAvg()  # by this, Single.val becomes the moving average

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

    # count the number of edges with non-zero weights
    nonzero = 0
    for i in range(n):
        for j in range(m):
            if W[i][j]>0:
                nonzero += 1
    
    # compute maximum weighted non-crossing matching
    if relax == True:
        M,opt,opt_all,cross = MWM(W)
    elif rand == True:
        M,opt,opt_all,cross = RAND(W)
    else:
        M,opt,opt_all,cross = MWNCM(W)
    M.sort(key=lambda e: e[0]) # M is ordered w.r.t. pseudo times
    St = [e[0] for e in M]     # St is the set of matched pseudo times
    I = [St[k+1]-St[k] for k in range(len(M)-1)]
 
    # output summary
    for dt in (Single, Bulk):
        if dt == Single:
            if sum_fp == None:
                print("=== Single ===")
            else:
                sum_fp.write("=== Single ===\n")
        else:
            if sum_fp == None:
                print("=== Bulk ===")
            else:
                sum_fp.write("=== Bulk ===\n")
        arr_peak = []
        arr_intv = []
        for I in dt.I:
            num_peak = 0
            for intv in I:
                num_peak += intv.right - intv.left + 1
            arr_peak.append(num_peak)
            arr_intv.append(len(I))
        if sum_fp == None:
            print(f"1. Number of peaks:\t{getAvg(arr_peak)}")
            print(f"2. Number of intervals:\t{getAvg(arr_intv)}")
            print(f"3. Number of all data points:\t{dt.rows}")
            print(f"4. Number of genes:\t{dt.cols}")
            print(f"5. 1/3:\t{getAvg(arr_peak)/dt.rows}")
            print(f"6. Nonzero-weight edges:\t{nonzero}")
            print(f"7. all possible edges:\t{n*m}")
            print(f"8. 6/7:\t{nonzero/(n*m)}")
            print(f"9. processing time:\t{time.time()-start_time}")
        else:
            sum_fp.write(f"1. Number of peaks:\t{getAvg(arr_peak)}\n")
            sum_fp.write(f"2. Number of intervals:\t{getAvg(arr_intv)}\n")
            sum_fp.write(f"3. Number of all data points:\t{dt.rows}\n")
            sum_fp.write(f"4. Number of genes:\t{dt.cols}\n")
            sum_fp.write(f"5. 1/3:\t{getAvg(arr_peak)/dt.rows}\n")
            sum_fp.write(f"6. Nonzero-weight edges:\t{nonzero}\n")
            sum_fp.write(f"7. all possible edges:\t{n*m}\n")
            sum_fp.write(f"8. 6/7:\t{nonzero/(n*m)}\n")
            sum_fp.write(f"9. processing time:\t{time.time()-start_time}\n")
            
    if sum_fp == None:
        print("=== Matching ===")
        print(f"1. Matching weight:\t{opt_all}")
        print(f"2. Matching weight without borders:\t{opt}")
        print(f"3. Matching edges:\t{len(M)}")
        print(f"4. Matching edges without borders:\t{len(M)-2}")
        print(f"5. Edge crossings in matching:\t{cross}")
    else:
        sum_fp.write(f"=== Matching ===\n")
        sum_fp.write(f"1. Matching weight:\t{opt_all}\n")
        sum_fp.write(f"2. Matching weight without borders:\t{opt}\n")
        sum_fp.write(f"3. Matching edges:\t{len(M)}\n")
        sum_fp.write(f"4. Matching edges without borders:\t{len(M)-2}\n")
        sum_fp.write(f"5. Edge crossings in matching:\t{cross}\n")

        
    # output the result
    k = -1
    A = []
    B = []
    for i in range(Single.rows):
        if i in St:
            k += 1
            A.append(M[k][1])
            B.append(Bulk.inter_rowlabel[M[k][1]])
        else:
            x = 0
            if k < len(St)-1:
                x = float( (i-St[k])*M[k+1][1] + (St[k+1]-i)*M[k][1] ) / float(St[k+1]-St[k])
            else:
                x = A[i-1] + EPSILON
            A.append(x)
            B.append('__NA__')
        if i>0:
            if pred_fp == None:
                print('{}\t{}\t{}\t{}\t{}'.format(i,Single.rowlabel[i],A[i],A[i]-A[i-1],B[i]))
            else:
                pred_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(i,Single.rowlabel[i],A[i],A[i]-A[i-1],B[i]))
        else:
            if pred_fp == None:
                print('{}\t{}\t{}\t__NA__\t{}'.format(i,Single.rowlabel[i],A[i],B[i]))
            else:
                pred_fp.write('{}\t{}\t{}\t__NA__\t{}\n'.format(i,Single.rowlabel[i],A[i],B[i]))


    if pred_fp != None:
        pred_fp.close()
    if sum_fp != None:
        sum_fp.close()
    
############################################################
if __name__=='__main__':
    if len(sys.argv)>=7:
        main(sys.argv)
        exit(0)
    sys.stdout.write('''
=== PeakMatch algorithm (v04) ===

The program outputs prediction in the following form:
(pseudo time ID) (pseudo time name) (estimated real time) (delay from last time) (label of matched bulk time) 


usage: {} (single)(bulk)(T)(last)(intv)(inter)[(option_1)(option_2)...]

  [mandatory]
  single ... tabu-separated file for single cell
  bulk   ... tabu-separated file for bulk
  T      ... threshold for deciding peak; it is ratio. >=1.0 is recommended.
  last   ... whether last times are forcibly matched (1) or not (0)
  intv   ... whether all points in an interval are regarded as peak points (1) or not (0; i.e., only max is considered)
  inter  ... number of points that are added between real time points for interpolation (e.g., 7 for 30min) 

  [options]
  -pred=<TEXT>    ... filename for outputting predicted times;
                      if not specified, output to stdout.
  -sum=<TEXT>     ... filename for outputting summary;
                      if not specified, output to stdout.
  -from=<INT>     ... ID of first gene in the target range;
                      0 is the leftmost gene.
  -to=<INT>       ... ID of last gene in the target range;
                      -1 is the rightmost gene.
  -relax          ... relax the non-crossing constraint;
                      i.e., typical maximum matching is computed
  -rand=<INT>     ... generate a non-crossing matching randomly;
                      the integer is used as a random seed. 

  --example--
  $ {} sample_single.txt sample_bulk.txt 1 1 1 7 -pred=pred.txt
  will output predicted times to pred.txt and data summary to stdout. 
'''.format(sys.argv[0], sys.argv[0]))
