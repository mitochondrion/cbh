#!/usr/local/bin/python3
import bio
import sys

print("HERE WE GO!!!!!!\n")

args = sys.argv

lines = []
params = []
sequence = ""

if len(args) > 1:
  file_name = args[1]
  lines = [ line.rstrip() for line in open(file_name, "r").readlines() ]
  sequence = lines[0]

if len(lines) > 1:
  params = lines[1].split(" ")



#####


answer = bio.min_gc_skews(sequence)

print(" ".join([str(i) for i in answer]))

