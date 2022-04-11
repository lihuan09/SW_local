#!/usr/bin/python
__author__ = "Huan Li"
__email__ = "huan.li@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import pandas as pd
import numpy as np

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()


### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    ### calculation
    
    #read in the score file, which is BLOSUM62.txt in this case
    file_scoreFile =open(scoreFile)
    scoreFile_r = file_scoreFile.read()
    scoreFile_list = scoreFile_r.split("\n")
    for i in range(len(scoreFile_list)):
        scoreFile_list[i] = scoreFile_list[i].split(" ")
    for i in range(len(scoreFile_list)):
        while True:
            try:
                scoreFile_list[i].remove("")
            except ValueError:
                break
    scoreFile_mat = []
    for line in range(len(scoreFile_list)):
        if line == 0:
            item = scoreFile_list[line][0:]
        else:
            item = scoreFile_list[line][1:]
        scoreFile_mat.append(item)
    # get the scores for corresponding characters
    score_list = scoreFile_mat[0]
    scoreFile_mat = scoreFile_mat[1:(len(scoreFile_list)-2)]
    scoreFile_mat = np.array(scoreFile_mat)
    
    #read in input file
    seq_file =open(inputFile)
    seq = seq_file.read()
    seq_list = seq.split("\n")
    
    a = seq_list[1] #rows 
    b = seq_list[0] #columns
    
    # create score matrix H and backtrace (position) matrix P
    H = np.zeros((len(a) + 1, len(b) + 1))
    P = np.zeros((len(a) + 1, len(b) + 1))
    
    # set up the position score
    DELETION, INSERTION, MATCH = range(3)
    
    # calculating the score matrix
    max_score = 0
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            x = score_list.index(a[i - 1])
            y = score_list.index(b[j - 1])

            alignment_score = int(scoreFile_mat[x][y])
            # or missmatch
            match = (H[i - 1, j - 1] + alignment_score, MATCH)

            gap_cost1 = [openGap + (k -1)*extGap for k in range(1,i+1)][::-1]
            delete = (max([H[0:i, j][n]+gap_cost1[n] for n in range(len(H[0:i, j]))]), DELETION)

            gap_cost2 = [openGap + (l -1)*extGap for l in range(1,j+1)][::-1]
            insert = (max([H[i, 0:j][n]+gap_cost2[n] for n in range(len(H[i, 0:j]))]), INSERTION)

            H[i, j], P[i, j] = max(match, delete, insert, (0,0))
            if H[i][j] >= max_score:
                max_i = i
                max_j = j
                #max score
                max_score = H[i][j]
                
    # backtrace and find the matching alignment
    m, n = max_i, max_j
    a_str = ''
    b_str = ''
    while m > 0 or n > 0 and H[m][n] != 0:
        if P[m][n] == MATCH:
            m -= 1
            n -= 1
            a_str += a[m]
            b_str += b[n]

        elif P[m][n] == INSERTION:
            n -= 1
            a_str += '-'
            b_str += b[n]

        elif P[m][n] == DELETION:
            m -=1
            a_str += a[m]
            b_str += '-'
    
    # consider the starting edge    
    if n == 0 or m == 0:
        if m >= n:
            start_a=a[0:(m-n)]
            start_b=' '*(n-m)
        elif m < n:
            start_b=b[0:(n-m)]
            start_a=' '*(m-n)
    else:
        if m >= n:
            start_a=a[0:(m)]
            start_b=' '*(m-n) + b[0:n]
        elif m < n:
            start_b=b[0:(n)]
            start_a=' '*(n-m) + a[0:m]
    
    # consider the ending edge          
    ending_a = len(a) - max_i
    ending_b = len(b) - max_j
    if ending_a==0 and ending_b==0:
        end_a = ""
        end_b = ""
    elif ending_b==0:
        end_a = a[-ending_a:]
        end_b = ' '*ending_a
    elif ending_a==0:
        end_a = ' '*ending_b
        end_b = b[-ending_b:]
    else:
        end_a = a[-ending_a:]
        end_b = b[-ending_b:]
    
    # putting all together
    a_string = start_a + "(" + a_str[::-1] + ")" + end_a
    b_string = start_b + "(" + b_str[::-1] + ")" + end_b
    
    H = H.astype(int)
    
    # insert the "|" accordingly
    def insert_pipe(align_a, align_b):
        pipe = ""
        if len(align_a) >= len(align_b):
            for i in range(0, len(align_b)):
                if align_b[i] == align_a[i] and (align_b[i] != "(" and align_b[i] != ")"):
                    pipe += "|"
                else:
                    pipe += " "
        else:   
            for i in range(0, len(align_a)):
                if align_a[i] == align_b[i] and (align_a[i] != "(" and align_a[i] != ")"):
                    pipe += "|"
                else:
                    pipe += " "

        return pipe
    ### write output
    print("""-----------
|Sequences|
-----------""")
    print("sequence1")
    print(b)
    print("sequence2")
    print(a)
    print("""--------------
|Score Matrix|
--------------""")
    x="\t".join([char for char in b])
    print("\t\t" + x)
    print("\t" + "\t".join([str(char) for char in H[0]]))
    for i in range(0,len(a)):
        print(a[i] + "\t" + "\t".join([str(char) for char in H[i+1]]))
    print("""----------------------
|Best Local Alignment|
----------------------""")
    print("Alignment Score:%s" %(int(max_score)))
    print("Alignment Results:")
    print(b_string)
    print(insert_pipe(a_string, b_string))
    print(a_string)

### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)