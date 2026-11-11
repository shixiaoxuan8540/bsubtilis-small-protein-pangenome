#!/usr/bin/env python3
"""Preliminary annotated <=100 aa proteins, not homolog-rescued pangenome."""
import csv,sys,pathlib
from Bio import SeqIO
root=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as handle:
 w=csv.writer(handle);w.writerow(['assembly_accession','protein_id','description','length_aa','sequence'])
 for faa in sorted(root.glob('*/protein.faa')):
  for p in SeqIO.parse(faa,'fasta'):
   if 0<len(p.seq)<=100 and '*' not in str(p.seq):
    w.writerow([faa.parent.name,p.id,p.description,len(p.seq),str(p.seq)])
