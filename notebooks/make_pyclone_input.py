#!/usr/bin/env python3
"""Build PyClone-VI input from a joint multi-sample Mutect2 VCF + ASCAT output.

Each --sample is NAME=CNV_DIR, where CNV_DIR holds ASCAT's segments.tsv
(columns: sample, chr, startpos, endpos, nMajor, nMinor) and purity_ploidy.tsv.
Use --purity NAME=VALUE to override purity (e.g., a low-purity border sample
where ASCAT fails and purity is estimated from clonal SNV VAFs).

Only PASS, biallelic, autosomal SNV/indels with copy number available in
every sample are kept (PyClone-VI needs each mutation in all samples).
"""
import argparse
from pathlib import Path

import pandas as pd
from cyvcf2 import VCF

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def kv(items):
    return dict(x.split("=", 1) for x in items or [])


def load_segments(cnv_dir: Path) -> pd.DataFrame:
    seg = pd.read_csv(cnv_dir / "segments.tsv", sep="\t")
    seg["chr"] = seg["chr"].astype(str)
    seg.loc[~seg["chr"].str.startswith("chr"), "chr"] = "chr" + seg["chr"]
    return seg


def load_purity(cnv_dir: Path) -> float:
    return float(pd.read_csv(cnv_dir / "purity_ploidy.tsv", sep="\t")["purity"].iloc[0])


def cn_at(seg: pd.DataFrame, chrom: str, pos: int):
    hit = seg[(seg["chr"] == chrom) & (seg["startpos"] <= pos) & (seg["endpos"] >= pos)]
    if hit.empty:
        return None
    r = hit.iloc[0]
    return int(r["nMajor"]), int(r["nMinor"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample", action="append", required=True, help="NAME=CNV_DIR")
    ap.add_argument("--purity", action="append", help="NAME=VALUE override")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    samples = kv(a.sample)
    overrides = {k: float(v) for k, v in kv(a.purity).items()}
    segs = {s: load_segments(Path(d)) for s, d in samples.items()}
    purity = {s: overrides.get(s) or load_purity(Path(d)) for s, d in samples.items()}

    vcf = VCF(a.vcf)
    idx = {s: vcf.samples.index(s) for s in samples}

    rows = []
    for v in vcf:
        if v.FILTER is not None or len(v.ALT) != 1 or v.CHROM not in AUTOSOMES:
            continue
        ad = v.format("AD")
        mut_id = f"{v.CHROM}:{v.POS}:{v.REF}>{v.ALT[0]}"
        recs = []
        for s, i in idx.items():
            cn = cn_at(segs[s], v.CHROM, v.POS)
            if cn is None or cn[0] == 0:
                break
            recs.append(dict(mutation_id=mut_id, sample_id=s,
                             ref_counts=int(ad[i][0]), alt_counts=int(ad[i][1]),
                             normal_cn=2, major_cn=cn[0], minor_cn=cn[1],
                             tumour_content=purity[s]))
        if len(recs) == len(idx):
            rows.extend(recs)

    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print(f"{len(rows) // len(idx)} mutations x {len(idx)} samples -> {a.out}")


if __name__ == "__main__":
    main()
