#!/usr/bin/env python3
import csv, sys, re, collections
rows=list(csv.DictReader(open(sys.argv[1],newline='')))
expected={'reference':2,'aquatic':3,'soil':3,'human_gut':3,'chicken_gut':3,'fish_gut':3,'food':3}
counts=collections.Counter(r['niche'] for r in rows)
errors=[]
if counts!=expected: errors.append(f'Unexpected niche counts: {dict(counts)}')
for r in rows:
    if r['metadata_verified'].lower()=='yes':
        for field in ['strain','assembly_accession','biosample_accession','original_isolation_source','source_evidence_url','taxonomic_verification']:
            if not r[field].strip(): errors.append(f"{r['slot']}: missing {field}")
        if not re.fullmatch(r'GC[AF]_\d+\.\d+',r['assembly_accession']): errors.append(f"{r['slot']}: malformed accession")
seen=[r['assembly_accession'] for r in rows if r['assembly_accession']]
if len(seen)!=len(set(seen)): errors.append('Duplicate assembly accession')
print(f'{len(rows)} rows; {sum(r["metadata_verified"].lower()=="yes" for r in rows)} verified; {len(errors)} errors')
for e in errors: print('ERROR:',e)
sys.exit(bool(errors))
