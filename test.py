#!/usr/local/bin/python3
import bio
import sys

print("HERE WE GO!!!!!!\n")

args = sys.argv

lines = []
params = []
sequence = ""
kmer = ""

if len(args) > 1:
  file_name = args[1]
  lines = [ line.rstrip() for line in open(file_name, "r").readlines() ]

  sequence = lines[1]
  kmer = lines[0]
  params = lines[2].split(" ")




###############


#print("seq len: ", len(sequence), ", kmer: ", kmer, ", fuzziness: ", int(params[0]), "\n\n")

#answer = bio.fuzzy_occurrences(sequence, kmer, int(params[0]))

#print(" ".join([str(i) for i in answer]))

#print("\n\nRESULTS: ", len(answer))

print(bio.mutations("AAAAA", 1))
