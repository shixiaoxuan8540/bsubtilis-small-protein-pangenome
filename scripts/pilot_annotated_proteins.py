#!/usr/bin/env python3
"""Preliminary annotated short protein inventory; NOT an orthology/pangenome analysis."""
import csv, pathlib, collections, sys, json
root=pathlib.Path(sys.argv[1]); out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
rows=[]; seqgroups=collections.defaultdict(set); qc=[]
for d in sorted((root/'ncbi_dataset'/'data').glob('GC*')):
    acc=d.name
    f=d/'protein.faa'
    g=next(d.glob('*_genomic.fna'),None)
    if not f.exists():
        qc.append({'accession':acc,'status':'missing protein.faa','genome_bp':'','contigs':''});continue
    length=0; nseq=0
    if g:
        for line in g.open():
            if line.startswith('>'):nseq+=1
            else:length+=len(line.strip())
    name=None; sequence=[]
    def process():
        if name is None:return
        s=''.join(sequence).upper().rstrip('*')
        if 0<len(s)<=100 and set(s)<=set('ACDEFGHIKLMNPQRSTVWYUXO'):
            rows.append({'accession':acc,'protein_id':name,'length_aa':len(s),'sequence':s})
            seqgroups[s].add(acc)
    for line in f.open():
        if line.startswith('>'):
            process();name=line[1:].strip().split()[0];sequence=[]
        else:sequence.append(line.strip())
    process()
    qc.append({'accession':acc,'status':'ok','genome_bp':length,'contigs':nseq})
with (out/'annotated_small_proteins.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['accession','protein_id','length_aa','sequence']);w.writeheader();w.writerows(rows)
with (out/'assembly_qc.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['accession','status','genome_bp','contigs'],extrasaction='ignore');w.writeheader();w.writerows(qc)
counts=collections.Counter(r['accession'] for r in rows)
summary={'genomes_with_protein_faa':sum(r['status']=='ok' for r in qc),'annotated_small_protein_counts':dict(counts),'distinct_exact_amino_acid_sequences':len(seqgroups),'note':'Exact sequence matches are NOT ortholog families; annotation-based counts are a preliminary baseline.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
