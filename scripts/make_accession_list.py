#!/usr/bin/env python3
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1],newline='')))
if len(rows)!=20 or any(r['metadata_verified'].lower()!='yes' for r in rows):
    sys.exit('All 20 rows must be reviewed and marked metadata_verified=yes before downloading.')
acc=[r['assembly_accession'] for r in rows]
if len(set(acc))!=20: sys.exit('Accessions must be unique')
print('\n'.join(acc))
