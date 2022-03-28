#!/usr/bin/python
############################################################
#
#   Python code for computing evaluators
#
############################################################

import os,sys
from peakmatch_lib import *

Stime_stp = dict()
Stime_stp["stp0.1"] = -24
Stime_stp["stp0.2"] = -24
Stime_stp["stp0.3"] = -24
Stime_stp["stp0.4"] = -24
Stime_stp["stp0.5"] = -24
Stime_stp["stp0.6"] = -24
Stime_stp["stp0.7"] = -24
Stime_stp["stp1.1"] = -20
Stime_stp["stp1.2"] = -20
Stime_stp["stp1.3"] = -20
Stime_stp["stp1.4"] = -20
Stime_stp["stp1.5"] = -20
Stime_stp["stp1.6"] = -20
Stime_stp["stp2.1"] = -16
Stime_stp["stp2.2"] = -16
Stime_stp["stp2.3"] = -16
Stime_stp["stp2.4"] = -16
Stime_stp["stp2.5"] = -16
Stime_stp["stp2.6"] = -16
Stime_stp["stp3.1"] = -12
Stime_stp["stp3.2"] = -12
Stime_stp["stp3.3"] = -12
Stime_stp["stp3.4"] = -12
Stime_stp["stp3.5"] = -12
Stime_stp["stp3.6"] = -12
Stime_stp["stp3.7"] = -12
Stime_stp["stp3.8"] = -12
Stime_stp["stp4.1"] = -8
Stime_stp["stp4.10"] = -8
Stime_stp["stp4.11"] = -8
Stime_stp["stp4.2"] = -8
Stime_stp["stp4.3"] = -8
Stime_stp["stp4.4"] = -8
Stime_stp["stp4.5"] = -8
Stime_stp["stp4.6"] = -8
Stime_stp["stp4.7"] = -8
Stime_stp["stp4.8"] = -8
Stime_stp["stp4.9"] = -8
Stime_stp["stp5.1"] = -4
Stime_stp["stp5.2"] = -4
Stime_stp["stp5.3"] = -4
Stime_stp["stp5.4"] = -4
Stime_stp["stp5.5"] = -4
Stime_stp["stp5.6"] = -4
Stime_stp["stp6.7"] = 0
Stime_stp["stp6.9"] = 0
Stime_stp["stp7.6"] = 0
Stime_stp["stp6.2"] = 0
Stime_stp["stp8.3"] = 0
Stime_stp["stp6.4"] = 0
Stime_stp["stp6.5"] = 0
Stime_stp["stp27.1"] = 4
Stime_stp["stp7.1"] = 4
Stime_stp["stp24.3"] = 4
Stime_stp["stp6.10"] = 4
Stime_stp["stp10.6"] = 4
Stime_stp["stp24.4"] = 4
Stime_stp["stp6.1"] = 8
Stime_stp["stp8.9"] = 8
Stime_stp["stp24.7"] = 8
Stime_stp["stp9.8"] = 8
Stime_stp["stp7.4"] = 8
Stime_stp["stp7.5"] = 8
Stime_stp["stp22.3"] = 8
Stime_stp["stp9.5"] = 8
Stime_stp["stp7.2"] = 12
Stime_stp["stp8.5"] = 12
Stime_stp["stp8.4"] = 12
Stime_stp["stp8.8"] = 12
Stime_stp["stp26.3"] = 12
Stime_stp["stp9.4"] = 12
Stime_stp["stp10.3"] = 12
Stime_stp["stp9.7"] = 12
Stime_stp["stp9.2"] = 12
Stime_stp["stp9.10"] = 16
Stime_stp["stp12.5"] = 16
Stime_stp["stp8.1"] = 16
Stime_stp["stp7.3"] = 16
Stime_stp["stp10.5"] = 16
Stime_stp["stp9.3"] = 20
Stime_stp["stp9.6"] = 20
Stime_stp["stp11.6"] = 20
Stime_stp["stp14.3"] = 20
Stime_stp["stp9.1"] = 20
Stime_stp["stp10.4"] = 20
Stime_stp["stp10.1"] = 20
Stime_stp["stp11.2"] = 20
Stime_stp["stp12.2"] = 24
Stime_stp["stp12.1"] = 24
Stime_stp["stp11.4"] = 24
Stime_stp["stp13.9"] = 24
Stime_stp["stp8.2"] = 24
Stime_stp["stp11.3"] = 24
Stime_stp["stp12.8"] = 24
Stime_stp["stp12.9"] = 24
Stime_stp["stp12.10"] = 24
Stime_stp["stp8.6"] = 28
Stime_stp["stp11.1"] = 28
Stime_stp["stp12.7"] = 28
Stime_stp["stp14.2"] = 28
Stime_stp["stp11.5"] = 28
Stime_stp["stp14.8"] = 28
Stime_stp["stp11.7"] = 28
Stime_stp["stp13.3"] = 28
Stime_stp["stp14.6"] = 32
Stime_stp["stp12.6"] = 32
Stime_stp["stp12.4"] = 32
Stime_stp["stp11.11"] = 32
Stime_stp["stp13.8"] = 32
Stime_stp["stp14.9"] = 36
Stime_stp["stp15.6"] = 36
Stime_stp["stp13.4"] = 36
Stime_stp["stp15.8"] = 36
Stime_stp["stp26.2"] = 36
Stime_stp["stp13.5"] = 36
Stime_stp["stp19.4"] = 36
Stime_stp["stp19.7"] = 40
Stime_stp["stp17.2"] = 40
Stime_stp["stp15.1"] = 40
Stime_stp["stp13.1"] = 40
Stime_stp["stp13.2"] = 44
Stime_stp["stp16.2"] = 44
Stime_stp["stp17.4"] = 44
Stime_stp["stp15.4"] = 44
Stime_stp["stp17.3"] = 44
Stime_stp["stp15.3"] = 44
Stime_stp["stp17.5"] = 48
Stime_stp["stp17.1"] = 48
Stime_stp["stp13.7"] = 48
Stime_stp["stp18.5"] = 48
Stime_stp["stp16.5"] = 48
Stime_stp["stp23.7"] = 48
Stime_stp["stp21.1"] = 48
Stime_stp["stp27.2"] = 52
Stime_stp["stp18.6"] = 52
Stime_stp["stp23.3"] = 52
Stime_stp["stp21.4"] = 52
Stime_stp["stp24.1"] = 52
Stime_stp["stp24.5"] = 52
Stime_stp["stp27.4"] = 56
Stime_stp["stp25.2"] = 56
Stime_stp["stp15.2"] = 56
Stime_stp["stp23.6"] = 60
Stime_stp["stp16.1"] = 60
Stime_stp["stp15.5"] = 60
Stime_stp["stp16.9"] = 60
Stime_stp["stp24.2"] = 60
Stime_stp["stp20.2"] = 64
Stime_stp["stp19.8"] = 64
Stime_stp["stp21.2"] = 64
Stime_stp["stp21.3"] = 64
Stime_stp["stp19.5"] = 64
Stime_stp["stp19.2"] = 64
Stime_stp["stp25.4"] = 68
Stime_stp["stp17.6"] = 68
Stime_stp["stp19.6"] = 68
Stime_stp["stp23.1"] = 68
Stime_stp["stp27.3"] = 72
Stime_stp["stp20.4"] = 72
Stime_stp["stp18.2"] = 72
Stime_stp["stp26.7"] = 72
Stime_stp["stp25.1"] = 72
Stime_stp["stp21.5"] = 72
Stime_stp["stp18.9"] = 72
Stime_stp["stp22.6"] = 76
Stime_stp["stp22.4"] = 76
Stime_stp["stp18.10"] = 76
Stime_stp["stp24.8"] = 76
Stime_stp["stp22.5"] = 80
Stime_stp["stp26.4"] = 80
Stime_stp["stp26.6"] = 80
Stime_stp["stp22.2"] = 80
Stime_stp["stp18.7"] = 80
Stime_stp["stp25.3"] = 80
Stime_stp["stp27.5"] = 80
Stime_stp["stp18.8"] = 84
Stime_stp["stp22.7"] = 84
Stime_stp["stp26.1"] = 84
Stime_stp["stp26.8"] = 84
Stime_stp["stp20.3"] = 84

Stime_icpsc = dict()
Stime_icpsc["O01"] = 0
Stime_icpsc["O02"] = 0
Stime_icpsc["O03"] = 0
Stime_icpsc["O04"] = 0
Stime_icpsc["O05"] = 0
Stime_icpsc["O06"] = 0
Stime_icpsc["O07"] = 0
Stime_icpsc["O12"] = 0
Stime_icpsc["B02"] = 9
Stime_icpsc["B03"] = 9
Stime_icpsc["B04"] = 9
Stime_icpsc["B05"] = 9
Stime_icpsc["B06"] = 9
Stime_icpsc["B08"] = 9
Stime_icpsc["B09"] = 9
Stime_icpsc["B13"] = 9
Stime_icpsc["D01"] = 38
Stime_icpsc["D02"] = 38
Stime_icpsc["D04"] = 38
Stime_icpsc["D06"] = 38
Stime_icpsc["D07"] = 38
Stime_icpsc["D08"] = 38
Stime_icpsc["D09"] = 38
Stime_icpsc["D13"] = 38
Stime_icpsc["F01"] = 60
Stime_icpsc["F02"] = 60
Stime_icpsc["F04"] = 60
Stime_icpsc["F05"] = 60
Stime_icpsc["F06"] = 60
Stime_icpsc["F08"] = 60
Stime_icpsc["F09"] = 60
Stime_icpsc["F10"] = 60
Stime_icpsc["H01"] = 84
Stime_icpsc["H02"] = 84
Stime_icpsc["H03"] = 84
Stime_icpsc["H04"] = 84
Stime_icpsc["H05"] = 84
Stime_icpsc["H06"] = 84
Stime_icpsc["H08"] = 84
Stime_icpsc["H09"] = 84
Stime_icpsc["J01"] = 106
Stime_icpsc["J03"] = 106
Stime_icpsc["J04"] = 106
Stime_icpsc["J05"] = 106
Stime_icpsc["J08"] = 106
Stime_icpsc["J09"] = 106
Stime_icpsc["J11"] = 106
Stime_icpsc["J12"] = 106
Stime_icpsc["L02"] = 144
Stime_icpsc["L04"] = 144
Stime_icpsc["L08"] = 144
Stime_icpsc["L10"] = 144
Stime_icpsc["L11"] = 144
Stime_icpsc["L14"] = 144
Stime_icpsc["L16"] = 144
Stime_icpsc["L17"] = 144
Stime_icpsc["M01"] = 216
Stime_icpsc["M02"] = 216
Stime_icpsc["M04"] = 216
Stime_icpsc["M06"] = 216
Stime_icpsc["M07"] = 216
Stime_icpsc["M12"] = 216
Stime_icpsc["M13"] = 216
Stime_icpsc["M16"] = 216

Stime_sample = dict()
for t in range(1,209):
    Stime_sample["PT{}".format(t)] = t

#################### functions ####################
def get_ratio(Single, Bulk, Intv, M):
    n = len(M)
    comp = 0
    cov = 0.0
    zero = 0
    for c in range(Single.cols):
        peaks = 0
        peak_match = 0
        for edge in M:
            i = edge['index']
            is_single_peak = False
            for s_intv in Single.I[c]:
                if s_intv.contains(Intv,i):
                    is_single_peak = True
                    break
            if is_single_peak:
                peaks += 1
                if "__NA__" not in edge['BulkLabel']:
                    is_bulk_peak = False
                    for b_intv in Bulk.I[c]:
                        if b_intv.contains(Intv,edge['BulkIndex']):
                            is_bulk_peak = True
                            break
                    if is_bulk_peak == True:
                        peak_match += 1            
        #if peaks == 0:
        #    print(f"{c} {peak_match} {peaks} N/A")
        #else:
        #    print(f"{c} {peak_match} {peaks} {float(peak_match)/peaks}")
        if peaks == peak_match:
            comp += 1
        if peaks == 0:
            zero += 1
        else:
            cov += float(peak_match)/float(peaks)
            
    return float(comp) / float(Single.cols), cov/(Single.cols-zero)

def get_kendall_tau(M):
    n = len(M)
    match = 0
    non_match = 0
    for i in range(n):
        pt_i = M[i]
        for j in range(i+1,n):
            pt_j = M[j]
            if pt_i['BulkIndex'] < pt_j['BulkIndex']:
                match += 1
            elif pt_i['BulkIndex'] > pt_j['BulkIndex']:
                non_match += 1
            #else:
            #    sys.stderr.write("warning: there are two PTs for which same real time is predicted.\n")
                
    return float(match-non_match) / float(n*(n-1)/2)

def get_R2(M):
    BulkIndexMax = None
    Stime = None
    for data in M:
        if Stime == None:
            if "stp" in data['PT']:
                Stime = Stime_stp
                PTMin,PTMax = -24,84
                fam_name = "PM_data"
            elif "PT" in data['PT']:
                Stime = Stime_sample
                PTMin,PTMax = 1,208
                fam_name = "Sample"
            else:
                Stime = Stime_icpsc
                PTMin,PTMax = 0,216
                fam_name =" iCpSc_data"
        '''
        else:
            if Stime == Stime_stp and ("stp" not in data["PT"]):
                sys.stderr.write("error: unknown family of pseudo time family.\n")
                exit(1)
            elif Stime == Stime_icpsc and ("stp" in data["PT"]):
                sys.stderr.write("error: unknown family of pseudo time family.\n")
                exit(1)
        '''
        if BulkIndexMax == None or data['BulkIndex'] > BulkIndexMax:
            BulkIndexMax = data['BulkIndex']
            
    if Stime == None:
        sys.stderr.write("warning: unknown pseudo time family.\n")
        sys.stderr.write("         R2 is not computed\n")
        return None, "None"
    x = []
    y = []    
    for data in M:
        ptlabel = data['PT']
        if ptlabel not in Stime:
            sys.stderr.write("warning: unknown pseudo time family.\n")
            sys.stderr.write("         R2 is not computed\n")
            return None, "None"
        x.append(Stime[ptlabel])
        y.append(float(PTMin) + (data['BulkIndex']/BulkIndexMax) * float(PTMax-PTMin))

    y_mean = getAvg(y)
    num = sum([(x[k]-y[k])**2 for k in range(len(M))])
    den = sum([(a-y_mean)**2 for a in y])
    R2 = num/den
    return 1.0-R2, fam_name

#################### main ####################
def main(argv):
    # read mandatory parameters
    T = float(argv[4])
    (Last, Intv, Inter) = map(int, argv[5:8])

    # read options
    from_col = 0
    to_col = -1
    for arg in argv[8:]:
        if "-from=" in arg:
            arr = arg.split('=')
            from_col = int(arr[1])
        elif "-to=" in arg:
            arr = arg.split('=')
            to_col = int(arr[1])
        else:
            sys.stderr.write("error: {} is not a valid option.\n".format(arg))
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

    # read matching from argv[3]
    fp = open(argv[3])
    M = []
    while 1:
        line = fp.readline()
        if line=="":
            break
        line = line.replace('\n','')
        s = line.split('\t')
        i = int(s[0])
        x = {}
        x['index'] = i
        x['PT'] = s[1]
        x['BulkIndex'] = float(s[2])
        x['diff'] = s[3]
        x['BulkLabel'] = s[4]
        M.append(x)
        
    comp, cov = get_ratio(Single, Bulk, Intv, M)
    print("Complete ratio:\t{}".format(comp))
    print("Coverage ratio:\t{}".format(cov))
    print("Kendall-tau distance:\t{}".format(get_kendall_tau(M)))
    R2,fam_name = get_R2(M)
    print("R^2:\t{}\t{}".format(R2,fam_name))
'''
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
'''

                
############################################################
if __name__=='__main__':
    if len(sys.argv)>=8:
        main(sys.argv)
        exit(0)
    sys.stdout.write('''
=== matching evaluator for PeakMatch ===

usage: {} (single)(bulk)(pred)(T)(last)(intv)(inter)[(option_1)(option_2)...]

  [mandatory]
  single   ... tabu-separated file for single cell
  bulk     ... tabu-separated file for bulk
  pred ... tabu-separated file for matching 
           (can be obtained by using -pred option of peakmatch_main.py)
  T      ... threshold for deciding peak; it is ratio. >=1.0 is recommended.
  last   ... whether last times are forcibly matched (1) or not (0); NOT USED in this program
  intv   ... whether all points in an interval are regarded as peak points (1) or not (0; i.e., only max is considered)
  inter  ... number of points that are added between real time points for interpolation (e.g., 7 for 30min) 

  [options]
  -from=<INT>     ... ID of first gene in the target range;
                      0 is the leftmost gene.
  -to=<INT>       ... ID of last gene in the target range;
                      -1 is the rightmost gene.

'''.format(sys.argv[0], sys.argv[0]))
