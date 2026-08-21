#!/usr/bin/env python3
"""Summarize and plot VCFtools windowed nucleotide diversity."""
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def key(x): return [int(v) if v.isdigit() else v.lower() for v in re.split(r"(\d+)", str(x))]

def main():
    p=argparse.ArgumentParser(); p.add_argument("--input-dir",required=True,type=Path); p.add_argument("--pattern",default="*.windowed.pi"); p.add_argument("--minimum-variants",type=int,default=1); p.add_argument("--contigs-file",type=Path); p.add_argument("--output-prefix",required=True,type=Path); a=p.parse_args()
    files=sorted(a.input_dir.glob(a.pattern))
    if not files or a.minimum_variants < 1: raise ValueError("No inputs or invalid minimum variants")
    frames=[]
    for f in files:
        d=pd.read_csv(f,sep="\t"); needed={"CHROM","BIN_START","BIN_END","N_VARIANTS","PI"}
        if not needed.issubset(d.columns): raise ValueError(f"Missing columns in {f}: {sorted(needed-set(d.columns))}")
        for c in ["BIN_START","BIN_END","N_VARIANTS","PI"]: d[c]=pd.to_numeric(d[c],errors="coerce")
        d=d.dropna(subset=list(needed)); d=d[d.N_VARIANTS>=a.minimum_variants].copy(); d["population"]=f.name.removesuffix(".windowed.pi"); frames.append(d)
    data=pd.concat(frames,ignore_index=True)
    if data.empty: raise ValueError("No windows remain")
    a.output_prefix.parent.mkdir(parents=True,exist_ok=True)
    summary=data.groupby("population").agg(n_windows=("PI","size"),mean_pi=("PI","mean"),median_pi=("PI","median"),sd_pi=("PI","std"),min_pi=("PI","min"),max_pi=("PI","max")).reset_index()
    summary.to_csv(f"{a.output_prefix}.summary.tsv",sep="\t",index=False)
    order=summary.sort_values("median_pi").population
    fig,ax=plt.subplots(figsize=(max(7,len(order)*1.2),6)); sns.violinplot(data=data,x="population",y="PI",order=order,inner=None,cut=0,color="#B8C4D9",ax=ax); sns.boxplot(data=data,x="population",y="PI",order=order,width=.18,showfliers=False,color="white",ax=ax); ax.tick_params(axis="x",rotation=35); ax.set_xlabel(""); ax.set_ylabel("Nucleotide diversity (π)"); fig.tight_layout()
    for ext in ("png","pdf"): fig.savefig(f"{a.output_prefix}.distribution.{ext}",dpi=300,bbox_inches="tight")
    plt.close(fig)
    contigs=[x.strip() for x in a.contigs_file.read_text().splitlines() if x.strip()] if a.contigs_file else sorted(data.CHROM.astype(str).unique(),key=key)
    pops=sorted(data.population.unique()); ncols=min(2,len(pops)); nrows=int(np.ceil(len(pops)/ncols)); fig,axes=plt.subplots(nrows,ncols,figsize=(8*ncols,4*nrows),squeeze=False)
    for ax,pop in zip(axes.flat,pops):
        sub=data[data.population==pop]; offset=0; ticks=[]; labels=[]
        for i,chrom in enumerate(contigs):
            c=sub[sub.CHROM.astype(str)==chrom].copy()
            if c.empty: continue
            c["POS"]=(c.BIN_START+c.BIN_END)/2+offset; ax.scatter(c.POS,c.PI,s=6,color=["#5C66A8","#AAB5D5"][i%2],linewidths=0); ticks.append((c.POS.min()+c.POS.max())/2); labels.append(chrom); offset=c.POS.max()+1_000_000
        ax.set_title(pop); ax.set_ylabel("π"); ax.set_xticks(ticks,labels,rotation=45,ha="right",fontsize=7)
    for ax in axes.flat[len(pops):]: ax.set_visible(False)
    fig.tight_layout()
    for ext in ("png","pdf"): fig.savefig(f"{a.output_prefix}.genome_scan.{ext}",dpi=300,bbox_inches="tight")
    plt.close(fig)

if __name__=="__main__": main()
