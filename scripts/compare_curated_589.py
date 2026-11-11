#!/usr/bin/env python3
"""Curated 589 reference peptide exact-match presence in annotated proteins.
One curated entry counted once per genome; no homology inference.
"""
import csv,collections,pathlib,sys,json,re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
out=pathlib.Path(sys.argv[1]);ref=pathlib.Path(sys.argv[2])
with ref.open(newline="") as f: curated=list(csv.DictReader(f))
assert len(curated)==589, f"Expected 589 curated rows; got {len(curated)}"
for i,r in enumerate(curated):
    seq=re.sub(r"\s+","",r["reference_sequence"]).upper().rstrip("*")
    if not seq or not set(seq)<=set("ACDEFGHIKLMNPQRSTVWYUXO"):raise ValueError(f"Invalid curated sequence row {i+1}")
    r["reference_sequence"]=seq
with (out/"annotated_small_proteins.csv").open(newline="") as f: annotations=list(csv.DictReader(f))
with (out/"genome_labels.csv").open(newline="") as f: labels={r["accession"]:r["label"] for r in csv.DictReader(f)}
genomes=collections.defaultdict(set)
for r in annotations:
    if r["accession"] in labels:genomes[r["accession"]].add(r["sequence"])
accessions=[a for a in labels if a in genomes]
# Exact equality to annotated translated protein, no <=100 filter on curated reference.
# Current source annotation inventory contains <=100aa proteins; flag reference >100 if present.
too_long=[r["name_synonyms"] for r in curated if len(r["reference_sequence"])>100]
if too_long:raise ValueError(f"Curated reference contains proteins >100 aa; baseline inventory must be regenerated: {too_long[:8]}")
matrix=[]
for i,r in enumerate(curated,1):
    present={a:int(r["reference_sequence"] in genomes[a]) for a in accessions}
    matrix.append({"curated_id":i,"name_synonyms":r["name_synonyms"],"reference_sequence":r["reference_sequence"],**present})
with (out/"curated_589_presence_matrix.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["curated_id","name_synonyms","reference_sequence"]+accessions);w.writeheader();w.writerows(matrix)
summary=[]
for a in accessions:
    n=sum(r[a] for r in matrix)
    summary.append({"accession":a,"label":labels[a],"matched_curated_entries":n,"total_curated_entries":589,"percent":round(n/589*100,2)})
with (out/"curated_589_genome_summary.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
# The same reference peptide may occur in more than one curated row. Count entries, not distinct peptides.
patterns=collections.Counter(tuple(a for a in accessions if r[a]) for r in matrix)
with (out/"curated_589_intersections.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["accessions","labels","curated_entries","number_of_genomes"])
    for combo,count in sorted(patterns.items(),key=lambda x:-x[1]):
        w.writerow([";".join(combo),";".join(labels[a] for a in combo),count,len(combo)])
# Top 30 observed nonempty exact presence patterns. Unmatched in all is reported separately.
top=[(k,v) for k,v in sorted(patterns.items(),key=lambda x:-x[1]) if k][:30]
fig=plt.figure(figsize=(19,max(9,5+len(accessions)*.48)))
gs=fig.add_gridspec(2,1,height_ratios=[2.2,1.8],hspace=.05)
ax=fig.add_subplot(gs[0]);bx=fig.add_subplot(gs[1],sharex=ax)
ax.bar(range(len(top)),[v for _,v in top],color="#55738a")
ax.set_title("589 curated B. subtilis proteins: exact annotated-protein matches by genome")
ax.set_ylabel("Curated reference entries")
for i,(_,v) in enumerate(top):ax.text(i,v+1,str(v),ha="center",fontsize=8)
ax.set_ylim(0,max(v for _,v in top)*1.18);ax.tick_params(axis="x",bottom=False,labelbottom=False)
for j,a in enumerate(accessions):
    y=len(accessions)-1-j
    bx.text(-.8,y,labels[a],ha="right",va="center",fontsize=9)
    for i,(combo,_) in enumerate(top):bx.scatter(i,y,s=35,color="#233f5b" if a in combo else "#d9dfe5",zorder=3)
for i,(combo,_) in enumerate(top):
    ys=[len(accessions)-1-j for j,a in enumerate(accessions) if a in combo]
    if len(ys)>1:bx.plot([i,i],[min(ys),max(ys)],color="#233f5b",lw=1.4)
bx.set_xlim(-.5,len(top)-.5);bx.set_ylim(-.6,len(accessions)-.4);bx.set_xticks([]);bx.set_yticks([])
for sp in bx.spines.values():sp.set_visible(False)
fig.subplots_adjust(left=.25,right=.98,bottom=.06,top=.93)
fig.savefig(out/"curated_589_upset.png",dpi=220);plt.close(fig)
selected=["GCF_000009045.1","CP020102","CP045826"]
if all(a in accessions for a in selected):
    groups=[{r["curated_id"] for r in matrix if r[a]} for a in selected]
    fig,ax=plt.subplots(figsize=(10,8))
    venn3(groups,set_labels=[labels[a] for a in selected],ax=ax)
    ax.set_title("589 curated proteins: exact-match presence (area-proportional approximation)")
    fig.tight_layout();fig.savefig(out/"curated_589_proportional_venn.png",dpi=220);plt.close(fig)
report={"reference_entries":len(curated),"unique_reference_sequences":len({r["reference_sequence"] for r in curated}),"genomes_compared":len(accessions),"unmatched_in_all":patterns.get((),0),"shared_in_all":patterns.get(tuple(accessions),0),"method":"Exact full-length peptide equality to native annotated proteins, not genomic absence or homolog detection; curated entries counted individually.","genomes":summary}
(out/"curated_589_summary.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2))
