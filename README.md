# SW_local

This program is to implement the Smith-Watermen local alighmnet algorithm for protein sequences.

The dependecies are including argparse, pandas, and numpy

blosum62.txt is provided as the score matrix for match and mismatch

The default gap penalities are as follows: openning gap = -2, extension gap = -1

There are two sample input and output files, additionally, one input text file and one output text file by running the algo

Please follow the instructions to run the program

Usage: python hw1.py -i <input file> -s <score file>

Example: python hw1.py -i input.txt -s blosum62.txt
  
