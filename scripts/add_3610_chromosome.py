#!/usr/bin/env python3
"""Add NCIB 3610 chromosome CP020102 only to an existing annotated-protein pilot."""
import csv, json, pathlib, sys
from Bio import SeqIO
from Bio.Seq import Seq
gb=pathlib.Path(sys.argv[1]); out=pathlib.Path(sys.argv[2])
records=list(SeqIO.parse(gb,"genbank"))
if len(records)!=1 or not records[0].id.startswith("CP020102"):
    raise SystemExit("Expected exactly one CP020102 chromosome record; refusing plasmid or mixed records")
rec=records[0]
if not 4_200_000 < len(rec.seq) < 4_230_000:
    raise SystemExit("Unexpected NCIB 3610 chromosome length")
rows=[]
for feat in rec.features:
    if feat.type!="CDS": continue
    translations=feat.qualifiers.get("translation",[])
    if not translations:continue
    seq=translations[0].upper().rstrip("*")
    if not (0<len(seq)<=100 and set(seq)<=set("ACDEFGHIKLMNPQRSTVWYUXO")):continue
    protein_id=feat.qualifiers.get("protein_id",feat.qualifiers.get("locus_tag",["unknown"]))[0]
    rows.append({"accession":"CP020102","protein_id":protein_id,"length_aa":len(seq),"sequence":seq})
p=out/"annotated_small_proteins.csv"
with p.open("a",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["accession","protein_id","length_aa","sequence"]);w.writerows(rows)
q=out/"assembly_qc.csv"
with q.open("a",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["accession","status","genome_bp","contigs"])
    w.writerow({"accession":"CP020102","status":"ok (chromosome only)","genome_bp":len(rec.seq),"contigs":1})
s=out/"summary.json"
data=json.loads(s.read_text())
data["genomes_with_protein_faa"]+=1
data["annotated_small_protein_counts"]["CP020102"]=len(rows)
with p.open(newline="") as f: data["distinct_exact_amino_acid_sequences"]=len({r["sequence"] for r in csv.DictReader(f)})
data["note"]="Chromosome CP020102 only for NCIB 3610; pBS32 CP020103 excluded. Annotation-based counts, not orthology."
s.write_text(json.dumps(data,indent=2)+"\n")
print(f"NCIB 3610 chromosome: {len(rec.seq)} bp, {len(rows)} annotated short proteins; plasmid excluded")
