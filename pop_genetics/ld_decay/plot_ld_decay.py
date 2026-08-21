#!/usr/bin/env python3
"""Bin PLINK pairwise r-squared reports and plot LD decay."""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

COLORS=["#5C66A8","#C97C7C","#6F9969","#EFC86E","#454A74","#8FB8B8"]

def main():
    p=argparse.ArgumentParser(); p.add_argument("--input-dir",required=True,type=Path); p.add_argument("--pattern",default="*.ld.gz"); p.add_argument("--bin-size-bp",type=int,default=10000); p.add_argument("--maximum-distance-bp",type=int); p.add_argument("--chunksize",type=int,default=1000000); p.add_argument("--output-prefix",required=True,type=Path); a=p.parse_args()
    if a.bin_size_bp<1 or a.chunksize<1: raise ValueError("Bin size and chunksize must be positive")
    files=sorted(a.input_dir.glob(a.pattern));
    if not files: raise ValueError("No PLINK LD reports found")
    rows=[]
    for f in files:
        sums={}; counts={}
        for chunk in pd.read_csv(f,sep=r"\s+",compression="infer",chunksize=a.chunksize,usecols=["BP_A","BP_B","R2"]):
            for c in ["BP_A","BP_B","R2"]: chunk[c]=pd.to_numeric(chunk[c],errors="coerce")
            chunk=chunk.dropna(); chunk["distance_bp"]=(chunk.BP_B-chunk.BP_A).abs(); chunk=chunk[chunk.distance_bp>0]
            if a.maximum_distance_bp is not None: chunk=chunk[chunk.distance_bp<=a.maximum_distance_bp]
            chunk["bin_start_bp"]=(chunk.distance_bp//a.bin_size_bp)*a.bin_size_bp
            g=chunk.groupby("bin_start_bp").R2.agg(["sum","count"])
            for b,r in g.iterrows(): sums[b]=sums.get(b,0.0)+r["sum"]; counts[b]=counts.get(b,0)+int(r["count"])
        label=f.name.removesuffix(".ld.gz")
        rows.extend({"population":label,"bin_start_bp":int(b),"bin_center_bp":int(b)+a.bin_size_bp/2,"mean_r2":sums[b]/counts[b],"n_pairs":counts[b]} for b in sorted(sums))
    out=pd.DataFrame(rows)
    if out.empty: raise ValueError("No SNP pairs remain")
    a.output_prefix.parent.mkdir(parents=True,exist_ok=True); out.to_csv(f"{a.output_prefix}.binned.tsv",sep="\t",index=False)
    fig,ax=plt.subplots(figsize=(8,6))
    for i,(label,d) in enumerate(out.groupby("population")): ax.plot(d.bin_center_bp,d.mean_r2,lw=1.8,color=COLORS[i%len(COLORS)],label=label)
    ax.set_xlabel("Physical distance (bp)"); ax.set_ylabel("Mean $r^2$"); ax.set_ylim(bottom=0); ax.legend(frameon=False); ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    for ext in ("png","pdf"): fig.savefig(f"{a.output_prefix}.{ext}",dpi=300,bbox_inches="tight")
    plt.close(fig)

if __name__=="__main__": main()
