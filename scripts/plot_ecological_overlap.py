#!/usr/bin/env python3
"""Generate exact-AA overlap figures. The full-cohort UpSet is proportional; 3-set Venn uses area fitting."""
import csv, collections, pathlib, sys, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out=pathlib.Path(sys.argv[1])
with (out/"annotated_small_proteins.csv").open(newline="") as f:
    rows=list(csv.DictReader(f))
with (out/"genome_labels.csv").open(newline="") as f:
    labels={r["accession"]:r["label"] for r in csv.DictReader(f)}
sets=collections.defaultdict(set)
for r in rows:
    if r["accession"] in labels:sets[r["accession"]].add(r["sequence"])
names=list(sets)
presence=collections.Counter()
for s in set.union(*sets.values()):
    presence[tuple(n for n in names if s in sets[n])]+=1
with (out/"all_intersections.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["genomes","labels","exact_sequence_count","intersection_size"])
    for combo,count in sorted(presence.items(),key=lambda kv:-kv[1]):
        w.writerow([";".join(combo),";".join(labels[n] for n in combo),count,len(combo)])
# UpSet-style plot: top 30 intersection patterns with exact proportional bars and membership matrix.
top=sorted(presence.items(),key=lambda kv:-kv[1])[:30]
fig=plt.figure(figsize=(19,10))
gs=fig.add_gridspec(2,1,height_ratios=[2.4,1.8],hspace=.04)
ax=fig.add_subplot(gs[0]);bx=fig.add_subplot(gs[1],sharex=ax)
ax.bar(range(len(top)),[v for _,v in top],color="#526d8a")
ax.set_ylabel("Exact amino-acid sequences")
ax.set_title(f"Annotated small proteins (≤100 aa): top {len(top)} exact intersections across {len(names)} genomes")
for i,(_,v) in enumerate(top):ax.text(i,v+max(2,max(v for _,v in top)*.01),str(v),ha="center",fontsize=7,rotation=90)
ax.set_ylim(0,max(v for _,v in top)*1.22);ax.tick_params(axis="x",bottom=False,labelbottom=False)
for j,n in enumerate(names):
    y=len(names)-1-j
    bx.text(-.7,y,labels[n],ha="right",va="center",fontsize=9)
    for i,(combo,_) in enumerate(top):
        active=n in combo
        bx.scatter(i,y,s=38,color="#263c56" if active else "#d6dce3",zorder=3)
for i,(combo,_) in enumerate(top):
    ys=[len(names)-1-j for j,n in enumerate(names) if n in combo]
    if len(ys)>1:bx.plot([i,i],[min(ys),max(ys)],color="#263c56",lw=1.5)
bx.set_xlim(-.5,len(top)-.5);bx.set_ylim(-.6,len(names)-.4);bx.set_yticks([]);bx.set_xticks([])
for spine in bx.spines.values():spine.set_visible(False)
fig.subplots_adjust(left=.23,right=.98,bottom=.07,top=.92)
fig.savefig(out/"all_genomes_upset.png",dpi=220);plt.close(fig)
# True 3-set area-fitting Venn with explicit warning if optional package unavailable.
try:
    from matplotlib_venn import venn3
    selected=["GCF_000009045.1","CP020102","CP045826"]
    if all(n in sets for n in selected):
        fig,ax=plt.subplots(figsize=(9,7))
        venn3([sets[n] for n in selected],set_labels=[labels[n] for n in selected],ax=ax)
        ax.set_title("Exact amino-acid sequence overlap (area-proportional approximation)")
        fig.tight_layout();fig.savefig(out/"three_genome_proportional_venn.png",dpi=220);plt.close(fig)
except ImportError:
    print("matplotlib-venn unavailable; skipping three-genome Venn")
print(json.dumps({"genomes":len(names),"distinct_sequences":len(set.union(*sets.values())),"shared_by_all":presence.get(tuple(names),0),"patterns":len(presence)},indent=2))
