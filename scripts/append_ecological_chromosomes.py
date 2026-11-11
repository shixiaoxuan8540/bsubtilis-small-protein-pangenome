#!/usr/bin/env python3
"""Append additional chromosome-only GenBank records to the expanded exact-sequence dataset.
No plasmids. Failed records remain in QC and are excluded from comparisons."""
import csv, json, pathlib, sys, time, urllib.request
from Bio import SeqIO

out=pathlib.Path(sys.argv[1])
new=[
 ("Soil 73","CP045826",4166516),
 ("Soil 75","CP045825",4156459),
 ("Soil MB8_B1","CP045823",4221278),
 ("Plant root UD1022","CP011534",4025326),
 ("Plant endophyte BSn5","CP002468",4093599),
]
# CGMCC 2108 excluded: available metadata labels isolation source as soil, not verified natto isolation.
p=out/"annotated_small_proteins.csv"
qcpath=out/"assembly_qc.csv"
with p.open(newline="") as f:
    reader=csv.DictReader(f); fields=reader.fieldnames; rows=list(reader)
with qcpath.open(newline="") as f:
    reader=csv.DictReader(f); qfields=reader.fieldnames; qc=list(reader)
assert fields==["accession","protein_id","length_aa","sequence"]
for label,acc,expected_bp in new:
    try:
        url=f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text"
        dest=out/(acc+".gb")
        for attempt in range(4):
            try:
                req=urllib.request.Request(url,headers={"User-Agent":"academic-comparative-genomics/1.0"})
                with urllib.request.urlopen(req,timeout=120) as r: dest.write_bytes(r.read())
                break
            except Exception:
                if attempt==3: raise
                time.sleep(3*(attempt+1))
        recs=list(SeqIO.parse(dest,"genbank"))
        if len(recs)!=1 or not recs[0].id.startswith(acc) or abs(len(recs[0].seq)-expected_bp)>2000:
            raise ValueError(f"Unexpected accession/length/records: {[(r.id,len(r.seq)) for r in recs]}")
        rec=recs[0]
        count=0
        for feat in rec.features:
            if feat.type!="CDS":continue
            aa=feat.qualifiers.get("translation",[""])[0].upper().rstrip("*")
            if not (0<len(aa)<=100 and set(aa)<=set("ACDEFGHIKLMNPQRSTVWYUXO")):continue
            pid=feat.qualifiers.get("protein_id",feat.qualifiers.get("locus_tag",["unknown"]))[0]
            rows.append({"accession":acc,"protein_id":pid,"length_aa":len(aa),"sequence":aa});count+=1
        qc.append({"accession":acc,"status":"ok (chromosome only)","genome_bp":len(rec.seq),"contigs":1})
        print(label,acc,"proteins:",count)
    except Exception as e:
        qc.append({"accession":acc,"status":"FAILED: "+str(e)[:150],"genome_bp":"","contigs":""})
        print(label,acc,"FAILED:",e)
with p.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
with qcpath.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=qfields);w.writeheader();w.writerows(qc)
counts={a:len({r["sequence"] for r in rows if r["accession"]==a}) for a in sorted({r["accession"] for r in rows})}
summary={"distinct_sequences_per_genome":counts,"distinct_exact_amino_acid_sequences":len({r["sequence"] for r in rows}),"genomes_compared":len(counts),"warning":"Exact protein sequence identity only, not homology. Native annotations are not standardized; failed records excluded. CP chromosomes exclude plasmids."}
(out/"expanded_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
(out/"genome_labels.csv").write_text("accession,label,niche\n"+ "\n".join(f"{a},{b},{c}" for a,b,c in [
("GCF_000009045.1","Reference 168","reference"),
("CP020102","Ancestral NCIB 3610","reference"),
("GCA_000321395.1","Chicken-associated BSP1","chicken"),
("GCA_000209795.2","Natto-associated BEST195","food"),
("GCA_025837035.1","Aquatic sediment MARUCo01","aquatic"),
("GCA_030297075.1","Natto-associated NARUSE","food"),
("CP045826","Soil 73","soil"),
("CP045825","Soil 75","soil"),
("CP045823","Soil MB8_B1","soil"),
("CP011534","Plant root UD1022","plant"),
("CP002468","Plant endophyte BSn5","plant"),
])+"\n")
