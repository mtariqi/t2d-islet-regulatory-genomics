import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

plt.rcParams.update({"font.family":"Liberation Sans","font.size":8,"axes.spines.top":False,
    "axes.spines.right":False,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "pdf.fonttype":42,"svg.fonttype":"none"})
# Okabe-Ito (colour-blind safe)
ND="#0072B2"; PD="#009E73"; T2D="#D55E00"; GREY="#8C8C8C"; INK="#222222"
DES="#56B4E9"; EDG="#E69F00"; BOTH="#4D4D4D"
MM=1/25.4

def save(fig, name):
    for ext in ("png","pdf","tiff"):
        kw = {"dpi":300}
        if ext=="tiff": kw["pil_kwargs"]={"compression":"tiff_lzw"}
        fig.savefig(f"{name}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)

def label(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

# ---------------- Figure 1: study design + cohort ----------------
fig = plt.figure(figsize=(180*MM, 95*MM))
axA = fig.add_axes([0.0, 0.0, 0.66, 0.95]); axA.axis("off"); axA.set_xlim(0,100); axA.set_ylim(0,100)
label(axA, "A", x=0.0, y=0.98)
def box(x,y,w,h,title,body,fc="#F2F5FA",ec="#9AA7BD",ls="-"):
    axA.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.4,rounding_size=1.8",fc=fc,ec=ec,lw=0.8,ls=ls))
    axA.text(x+w/2,y+h-3.2,title,ha="center",va="top",fontsize=7.6,fontweight="bold",color=INK)
    axA.text(x+w/2,y+h-8.6,body,ha="center",va="top",fontsize=6.6,color="#333",linespacing=1.3)
def arrow(x1,y1,x2,y2,ls="-"):
    axA.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=8,lw=0.8,color="#5A6478",ls=ls))
box(2,74,30,22,"Data and audit","GSE221156 (CELLxGENE)\n48 donors, 54 libraries\n99,029 beta cells\nraw/X integer counts")
box(36,74,28,22,"Donor pseudobulk","26,936 genes \u00d7 48 donors\nall 2,081,765,364\ncounts preserved")
box(68,74,30,22,"QC","PCA and 18 UMAPs\nchemistry structure\nno donor excluded")
box(68,40,30,26,"Three models","Primary (34 donors)\nBMI-adjusted (34)\nV3 only (23)\n17 T2D vs 17 ND")
box(36,40,28,26,"Two methods","DESeq2 (Wald)\nedgeR robust QL\nsame genes per model\nFDR < 0.05")
box(2,40,30,26,"Robustness","913 jointly significant\n285 core genes\n(6/6 combinations)\nleave-one-donor-out")
box(2,6,46,26,"Pathways","Ranked fgsea (Reactome, GO BP)\nredundancy clustering\nexcitability program (up)\nproteasome\u2013ubiquitin module (down)")
box(52,6,46,26,"Pre-registered replication","frozen 29-gene module + full pathway\ndecision rules, tag prereg-validation-v1\nindependent donor cohorts\n(not yet analysed)",fc="#FFFFFF",ec="#9AA7BD",ls="--")
arrow(32.8,85,35.2,85); arrow(64.8,85,67.2,85)
arrow(83,73.2,83,67.2)
arrow(67.2,53,64.8,53); arrow(35.2,53,32.8,53)
arrow(17,39.2,17,33.2); arrow(48.8,19,51.2,19,ls="--")

axB = fig.add_axes([0.74, 0.14, 0.25, 0.72])
label(axB, "B", x=-0.32, y=1.02)
groups=["ND","PD","T2D"]; v2=np.array([4,1,7]); v3=np.array([13,13,10])
cols=[ND,PD,T2D]
axB.bar(groups, v3, color=cols, alpha=0.45, width=0.62, edgecolor="none", label="V3")
axB.bar(groups, v2, bottom=v3, color=cols, width=0.62, edgecolor="none", label="V2")
for i,(a,b) in enumerate(zip(v3,v2)):
    axB.text(i, a/2, str(a), ha="center", va="center", fontsize=7, color=INK)
    axB.text(i, a+b/2, str(b), ha="center", va="center", fontsize=7, color="white", fontweight="bold")
axB.set_ylabel("Donors"); axB.set_ylim(0,19.5); axB.set_yticks([0,5,10,15])
axB.text(0.5,1.02,"Donors by chemistry",transform=axB.transAxes,ha="center",fontsize=7.6,fontweight="bold")
axB.text(1.02,0.15,"lighter: V3\ndarker: V2",transform=axB.transAxes,fontsize=6.4,color="#444",va="bottom")
save(fig,"Figure1_study_design")

# ---------------- Figure 3: method agreement + robustness funnel ----------------
fig, (a, b) = plt.subplots(1,2, figsize=(180*MM, 70*MM), gridspec_kw={"width_ratios":[1.05,1],"wspace":0.45})
models=["Primary","BMI-adjusted","V3 only"]
both=np.array([913,656,364]); donly=np.array([1144,892,662])-both; eonly=np.array([983,691,378])-both
y=np.arange(3)[::-1]
a.barh(y, both, color=BOTH, height=0.6, label="Both methods")
a.barh(y, donly, left=both, color=DES, height=0.6, label="DESeq2 only")
a.barh(y, eonly, left=both+donly, color=EDG, height=0.6, label="edgeR only")
for yi,bv,dv,ev in zip(y,both,donly,eonly):
    a.text(bv/2,yi,f"{bv:,}",ha="center",va="center",color="white",fontsize=7,fontweight="bold")
    a.text(bv+dv+ev+15,yi,f"{dv} / {ev}",va="center",fontsize=6.6,color="#444")
a.set_yticks(y); a.set_yticklabels(models); a.set_xlabel("Genes with FDR < 0.05"); a.set_xlim(0,1600)
a.legend(frameon=False, fontsize=6.6, loc="lower right", bbox_to_anchor=(1.0,1.0), ncol=3, handlelength=1, columnspacing=0.8)
label(a,"A",x=-0.28,y=1.1)
a.text(1.0,-0.3,"Right of bars: DESeq2-only / edgeR-only",transform=a.transAxes,ha="right",fontsize=6.2,color="#555")

steps=["Primary,\nboth methods","+ BMI-adjusted","+ V3 only\n(core genes)","+ strict LODO"]
vals=[913,611,285,None]
yy=np.arange(4)[::-1]
for yi,s,v in zip(yy,steps,vals):
    if v is None:
        b.barh(yi, 285, color="none", edgecolor="#888", ls="--", lw=0.8, height=0.6)
        b.text(142,yi,"pending",ha="center",va="center",fontsize=6.8,color="#666",style="italic")
    else:
        b.barh(yi, v, color=T2D if v==285 else "#B9C3D3", height=0.6)
        b.text(v+12,yi,f"{v:,}",va="center",fontsize=7.4,fontweight="bold" if v==285 else "normal")
b.set_yticks(yy); b.set_yticklabels(steps, fontsize=7); b.set_xlabel("Jointly significant genes, same direction"); b.set_xlim(0,1050)
b.text(1.0,-0.3,"Primary and V3-only overlap: 322 genes",transform=b.transAxes,ha="right",fontsize=6.2,color="#555")
label(b,"B",x=-0.36,y=1.1)
save(fig,"Figure3_method_agreement_and_robustness")

# ---------------- Figure 4A: core genes in excitability pathways ----------------
genes={"KCNA1":1.719,"KCNAB1":1.394,"KCNH8":1.416,"KCNAB2":1.095,"KCNJ2":0.986,"SCN3A":0.875,"NTRK2":1.423,
       "SYNDIG1":2.354,"ADGRB3":1.515,"GAP43":1.327,"CNTN5":1.093,"DSCAM":1.032,"LSAMP":0.846,"EFNA5":0.379,
       "GOT1":-0.620,"COA3":-0.440}
paths={"Reactome\npotassium channels":["KCNA1","KCNAB1","KCNH8","KCNAB2","KCNJ2"],
       "Reactome voltage-gated\npotassium channels":["KCNA1","KCNAB1","KCNAB2","KCNH8"],
       "GO action\npotential":["KCNA1","NTRK2","KCNJ2","SCN3A"],
       "GO synapse\nassembly":["SYNDIG1","ADGRB3","NTRK2","GAP43","CNTN5","DSCAM","LSAMP","EFNA5"],
       "Reactome aerobic\nrespiration":["GOT1","COA3"]}
order=list(genes.keys())
fig, (ax, axl) = plt.subplots(1,2, figsize=(180*MM, 72*MM), gridspec_kw={"width_ratios":[4.2,1],"wspace":0.05})
pn=list(paths.keys())
for j,g in enumerate(order):
    for i,p in enumerate(pn):
        if g in paths[p]:
            v=genes[g]; c=T2D if v>0 else ND
            ax.scatter(j, i, s=12+110*abs(v), color=c, alpha=0.9, edgecolor="white", lw=0.5)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=55, ha="right", fontstyle="italic")
ax.set_yticks(range(len(pn))); ax.set_yticklabels(pn, fontsize=7); ax.invert_yaxis()
ax.set_xlim(-0.7,len(order)-0.3); ax.set_ylim(len(pn)-0.4,-0.6)
ax.grid(color="#EEE", lw=0.5); ax.set_axisbelow(True)
for s in ("left","bottom"): ax.spines[s].set_visible(False)
ax.tick_params(length=0)
axl.axis("off")
for k,(v,lab) in enumerate([(2.0,"2.0"),(1.0,"1.0"),(0.5,"0.5")]):
    axl.scatter(0.25, 0.80-0.13*k, s=12+110*v, color="#999", transform=axl.transAxes, clip_on=False)
    axl.text(0.45, 0.80-0.13*k, lab, transform=axl.transAxes, va="center", fontsize=6.8)
axl.text(0.05,0.95,"|log2 FC|",transform=axl.transAxes,fontsize=7,fontweight="bold")
axl.scatter(0.25,0.32,s=60,color=T2D,transform=axl.transAxes,clip_on=False); axl.text(0.45,0.32,"higher in T2D",transform=axl.transAxes,va="center",fontsize=6.8)
axl.scatter(0.25,0.20,s=60,color=ND,transform=axl.transAxes,clip_on=False); axl.text(0.45,0.20,"lower in T2D",transform=axl.transAxes,va="center",fontsize=6.8)
label(ax,"A",x=-0.2,y=1.02)
save(fig,"Figure4A_excitability_core_genes")

# ---------------- Figure 5A: proteasome module across models ----------------
fig, (a, b) = plt.subplots(1,2, figsize=(150*MM, 62*MM), gridspec_kw={"wspace":0.55})
models=["Primary","BMI-\nadjusted","V3 only"]
dmed=[-0.182,-0.181,-0.082]; emed=[-0.181,-0.180,-0.088]
x=np.arange(3); w=0.34
a.bar(x-w/2, dmed, w, color=DES, label="DESeq2"); a.bar(x+w/2, emed, w, color=EDG, label="edgeR")
a.axhline(0,color="#444",lw=0.6)
for xi,dv,ev in zip(x,dmed,emed):
    a.text(xi,min(dv,ev)-0.008,f"{dv:.2f} / {ev:.2f}".replace("-","\u2212"),ha="center",va="top",fontsize=6.2)
a.set_xticks(x); a.set_xticklabels(models); a.set_ylabel("Median log2 FC of 29 genes\n(T2D vs ND)"); a.set_ylim(-0.23,0.02)
a.legend(frameon=False,fontsize=6.6,loc="lower left", bbox_to_anchor=(0,1.0), ncol=2)
label(a,"A",x=-0.42,y=1.08)
neg=[29,29,28]; sig=[2,0,0]
b.bar(x, neg, 0.55, color="#B9C3D3", label="Negative in both methods")
b.bar(x, sig, 0.55, color=T2D, label="Significant in both methods")
for xi,n,s in zip(x,neg,sig):
    b.text(xi,n+0.6,f"{n}/29",ha="center",fontsize=6.8)
    b.text(xi,max(s,0)+0.6,str(s),ha="center",fontsize=6.8,color=T2D,fontweight="bold")
b.set_xticks(x); b.set_xticklabels(models); b.set_ylabel("Genes"); b.set_ylim(0,33)
b.legend(frameon=False,fontsize=6.4,loc="lower left", bbox_to_anchor=(0,1.0))
label(b,"B",x=-0.32,y=1.08)
save(fig,"Figure5A_proteasome_module_summary")
print("done")
