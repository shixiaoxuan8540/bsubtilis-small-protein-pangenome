#!/usr/bin/env python3
"""Descriptive exact-AA-sequence overlap of annotated short proteins. NOT orthology."""
import csv,sys,itertools,collections
input_csv,out_csv=sys.argv[1:3]
sets=collections.defaultdict(set)
with open(input_csv,newline='') as f:
 for r in csv.DictReader(f):
  sets[r['accession']].add(r['sequence'])
names=sorted(sets)
if len(names)<2:sys.exit('Need at least two genomes')
with open(out_csv,'w',newline='') as f:
 w=csv.writer(f);w.writerow(['presence_pattern','accessions','distinct_exact_sequences'])
 for bits in itertools.product('01',repeat=len(names)):
  if not any(x=='1' for x in bits):continue
  on=[sets[names[i]] for i,b in enumerate(bits) if b=='1']
  off=[sets[names[i]] for i,b in enumerate(bits) if b=='0']
  seqs=set.intersection(*on)
  if off:seqs-=set.union(*off)
  w.writerow([''.join(bits),';'.join(names[i] for i,b in enumerate(bits) if b=='1'),len(seqs)])
print('Input genome order:',names)
