# PeakMatch
PeakMatch is a Python program that estimates the actual time when gene expressions in pseudo time based scRNA-seq data.

# Package
The package contains the following five files.  
• peakmatch_main.py: The main Python script for PeakMatch.  
• peakmatch_lib.py: A Python script that contains subroutines.  
• peakmatch_eval.py: A Python script that computes evaluation values for a matching computed by peakmatch_main.py.  
• sample_single.txt: A sample of single cell data (i.e., scRNA-seq data).  
• sample_bulk.txt: A sample of bulk data (i.e., cpRNA-seq data).  

# DEMO
To run the program, launch the terminal, and run **peakmatch_main.py** with some additional arguments as follows  
`$ python3 peakmatch_main . py sample_single . txt sample_bulk . txt \1 1 1 7 - pred = sample_pred . txt`
    
• By the 1st argument, you should specify the name of the file that contains scRNA-seq
data. In the above example, it is **sample single.txt**.  
• Similarly, by the 2nd argument, you should specify the name of the file that contains
cpRNA-seq data. In the above example, it is **sample bulk.txt**.  
• The subsequent four integer arguments are a bit technical. See Section 5 for detail.
If you are not interested in the detail, you may consider **1 1 1 7** as default values
for the time being.  
• The last argument **-pred=sample_pred.txt** is optional. The option **-pred** specifies
the file name to which prediction result is written. The generated file is required to
compute evaluation values of the matching by running **peakmatch_eval.py**. If the
option is not given, prediction result is presented to the standard output.  

The program will output the following to the standard output.  

    === Single ===
    1. Number of peaks : 21.8
    2. Number of intervals : 9.0
    3. Number of all data points : 208
    4. Number of genes : 5
    5. 1/3: 0.10480769230769231
    6. Nonzero - weight edges : 2693
    7. all possible edges : 45136
    8. 6/7: 0.05966412619638426
    9. processing time : 0.1721482276916504
    === Bulk ===
    1. Number of peaks : 25.8
    2. Number of intervals : 4.2
    3. Number of all data points : 217
    4. Number of genes : 5
    5. 1/3: 0.11889400921658987
    6. Nonzero - weight edges : 2693
    7. all possible edges : 45136
    8. 6/7: 0.05966412619638426
    9. processing time : 0.17216205596923828
    === Matching ===
    1. Matching weight : 436.11427663811986
    2. Matching weight without borders : 0.11133244969149594
    3. Matching edges : 58
    4. Matching edges without borders : 56
    5. Edge crossings in matching : 0  
    
The above shows some statistics on the input and the output. It can be output to a file,
where the file name can be specified by option **-sum**.  
Besides, the program generates a file named **sample_pred.txt**, which was specified by
**-pred**, whose contents are as follows.

    0 PT1 0 __NA__ AT1
    1 PT2 1 1 AT1_1
    2 PT3 2 1 AT1_2
    3 PT4 3 1 AT1_3
    4 PT5 4 1 AT1_4
    5 PT6 7.2 3.2 __NA__
    ...
    
• The 1st column indicates the index of a pseudo time. The index begins from zero.  
• The 2nd column indicates the label of the pseudo time.  
• The 3rd column indicates the estimated actual time, which is normalized to a value
between 0 and the number of bulk nodes minus one.  
• The 4th column indicates the delay; in the above example, because PT5 is estimated
as 4 and PT6 is estimated as 7.2, the delay shown in the PT6’s row is 7.2 − 4 = 3.2.
3  

See PDF for further details.

# Author
* HARAGUCHI, Kazuya
* As of writing, the author is with Department of Applied Mathematics and Physics, Graduate School
of Informatics, Kyoto University, Japan. 
* E-mail: dr.kazuya.haraguchi@gmail.com
