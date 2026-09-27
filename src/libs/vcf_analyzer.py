#!/usr/bin/env python3
"""
vcf_analyzer.py - analyse VCF/BCF files with bcftools, driven from Python.

bcftools does the heavy lifting (filtering, splitting multi-allelics, optional
REF check against a FASTA, stats, field extraction); Python classifies variants
and summarises genotypes, allele fractions, ploidy, copy number and aneuploidy.

What it reports, per input file (outdir/<vcf_name>/):
  bcftools_stats.txt   raw `bcftools stats -s -` output
  records.tsv          one row per (split) record: class, size, scale, flags
  genotypes.tsv        one row per record x sample: GT, ploidy, dosage, zygosity, VAF, CN
  ploidy.tsv           per-sample genotype ploidy profile (haploid ... polyploid)
  cn_events.tsv        per-sample copy-number events (gain/loss/amp/homdel/LOH)
  aneuploidy.tsv       per-sample x chromosome fraction gained/lost -> whole-chrom / arm calls
  somatic.tsv          tumor-vs-normal calls, when a tumor/normal pair is detected
  summary.json         headline numbers
And across all inputs (outdir/):
  summary.tsv, variant_classes.png, copy_number.png

Usage:
  python vcf_analyzer.py sample.vcf.gz
  python vcf_analyzer.py *.vcf -o results --pass-only --ref hg38.fa
  python vcf_analyzer.py tumor.vcf.gz -i 'QUAL>30 && INFO/DP>10' --baseline-ploidy 4

Requires: bcftools on PATH (or --bcftools), pandas; matplotlib for plots.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

# ----------------------------------------------------------------------------- bcftools

class BcftoolsError(RuntimeError):
    pass


class Bcftools:
    def __init__(self, exe: str = "bcftools", threads: int = 2):
        path = shutil.which(exe)
        if path is None:
            sys.exit(f"bcftools not found ('{exe}'). Install it or pass --bcftools /path/to/bcftools")
        self.exe, self.threads = path, threads
        self.version = self.run(["--version"]).splitlines()[0]

    def run(self, args: list[str]) -> str:
        p = subprocess.run([self.exe, *args], capture_output=True, text=True)
        if p.returncode != 0:
            raise BcftoolsError(f"bcftools {' '.join(args)}\n{p.stderr.strip()}")
        return p.stdout

    def pipe(self, first: list[str], second: list[str]) -> None:
        """Run `bcftools first | bcftools second` without a temp file."""
        p1 = subprocess.Popen([self.exe, *first], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.run([self.exe, *second], stdin=p1.stdout, capture_output=True, text=True)
        p1.stdout.close()
        err1 = p1.stderr.read().decode()
        if p1.wait() != 0 or p2.returncode != 0:
            raise BcftoolsError(f"{' '.join(first)} | {' '.join(second)}\n{err1}\n{p2.stderr}")


# ----------------------------------------------------------------------------- header

@dataclass
class Header:
    fileformat: str = ""
    samples: list[str] = field(default_factory=list)
    contigs: dict[str, int] = field(default_factory=dict)
    info: set[str] = field(default_factory=set)
    fmt: set[str] = field(default_factory=set)
    tumor: str | None = None
    normal: str | None = None


def read_header(bt: Bcftools, vcf: str) -> Header:
    h = Header()
    for line in bt.run(["view", "-h", vcf]).splitlines():
        if line.startswith("##fileformat="):
            h.fileformat = line.split("=", 1)[1]
        elif line.startswith("##contig="):
            cid = re.search(r"ID=([^,>]+)", line)
            ln = re.search(r"length=(\d+)", line)
            if cid:
                h.contigs[cid.group(1)] = int(ln.group(1)) if ln else 0
        elif line.startswith("##INFO=<ID="):
            h.info.add(re.search(r"ID=([^,>]+)", line).group(1))
        elif line.startswith("##FORMAT=<ID="):
            h.fmt.add(re.search(r"ID=([^,>]+)", line).group(1))
        elif line.startswith("##tumor_sample="):
            h.tumor = line.split("=", 1)[1]
        elif line.startswith("##normal_sample="):
            h.normal = line.split("=", 1)[1]
    h.samples = [s for s in bt.run(["query", "-l", vcf]).split("\n") if s]
    if h.tumor is None or h.normal is None:  # fall back to sample names
        t = [s for s in h.samples if re.search(r"tum", s, re.I)]
        n = [s for s in h.samples if re.search(r"norm|germ|blood", s, re.I)]
        if len(t) == 1 and len(n) == 1:
            h.tumor, h.normal = t[0], n[0]
    return h


# ----------------------------------------------------------------------------- bcftools stats

def parse_stats(text: str) -> dict:
    sn, tstv, idd, psc = {}, {}, [], []
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "SN":
            sn[f[2].rstrip(":").replace("number of ", "")] = int(f[3])
        elif f[0] == "TSTV":
            tstv = {"ts": int(f[2]), "tv": int(f[3]), "ts/tv": float(f[4])}
        elif f[0] == "IDD":
            idd.append({"indel_length": int(f[2]), "count": int(f[3])})
        elif f[0] == "PSC":  # per-sample counts
            psc.append({"sample": f[2], "nRefHom": int(f[3]), "nNonRefHom": int(f[4]),
                        "nHets": int(f[5]), "nTransitions": int(f[6]), "nTransversions": int(f[7]),
                        "nIndels": int(f[8]), "nSingletons": int(f[10]), "nMissing": int(f[13])})
    return {"SN": sn, "TSTV": tstv, "IDD": idd, "PSC": psc}


# ----------------------------------------------------------------------------- classification

SYMBOLIC = re.compile(r"^<([^>]+)>$")
BREAKEND = re.compile(r"[\[\]]")


def classify(ref: str, alt: str, svtype: str | None) -> str:
    m = SYMBOLIC.match(alt)
    if m:
        tag = m.group(1).split(":")[0].upper()
        if tag == "*" or tag == "NON_REF":
            return "REF_BLOCK"
        if tag.startswith("STR"):
            return "STR"
        return tag if tag in {"DEL", "DUP", "INV", "INS", "CNV", "BND"} else (svtype or tag)
    if BREAKEND.search(alt):
        return "BND"
    if alt in {".", "*"}:
        return "NO_ALT"
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    if len(ref) == 1 and alt.startswith(ref):
        return "INS"
    if len(alt) == 1 and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def event_size(cls: str, pos: int, ref: str, alt: str, end, svlen) -> int | None:
    if cls in {"SNV", "MNV"}:
        return len(ref)
    if cls in {"INS", "DEL", "COMPLEX"} and not SYMBOLIC.match(alt):
        return abs(len(alt) - len(ref)) or len(ref)
    if svlen not in (None, "."):
        return abs(int(str(svlen).split(",")[0]))
    if end not in (None, "."):
        return int(end) - pos + 1
    return None


def event_scale(span: int | None, contig_len: int) -> str:
    """Size class for copy-number / SV events, relative to the chromosome."""
    if not span:
        return "."
    if contig_len and span >= 0.9 * contig_len:
        return "WHOLE_CHROMOSOME"
    if (contig_len and span >= 0.25 * contig_len) or span >= 10_000_000:
        return "ARM_LEVEL"
    if span >= 1_000_000:
        return "LARGE"
    if span >= 50:
        return "FOCAL_SV"
    return "SMALL"


def repeat_flag(cls: str, ref: str, alt: str, ru) -> str:
    """Label repeat-associated indels.
    MSI_SLIPPAGE      1-2 repeat units gained/lost (typical of mismatch-repair deficiency)
    REPEAT_EXPANSION  >=3 units gained/lost (e.g. HTT CAG expansion), or a symbolic <STR..> allele
    Without a reference FASTA the repeat context is unknown, so INFO/RU is used when present and a
    mononucleotide 1-2 bp indel is otherwise treated as a homopolymer slippage candidate."""
    if cls == "STR":
        return "REPEAT_EXPANSION"
    if cls not in {"INS", "DEL"}:
        return ""
    delta = alt[1:] if cls == "INS" else ref[1:]
    unit = ru if ru not in (None, ".", "") else (delta[0] if len(set(delta)) == 1 else None)
    if not unit or len(delta) % len(unit) or delta != unit * (len(delta) // len(unit)):
        return ""
    n = len(delta) // len(unit)
    return "MSI_SLIPPAGE" if n <= 2 else "REPEAT_EXPANSION"


# ----------------------------------------------------------------------------- genotypes

PLOIDY_NAME = {0: "missing", 1: "haploid", 2: "diploid", 3: "triploid", 4: "tetraploid",
               5: "pentaploid", 6: "hexaploid"}


def parse_gt(gt: str) -> tuple[int, int, str]:
    """Return (ploidy, alt dosage, zygosity) for any ploidy, e.g. '0/0/1/1' -> (4, 2, 'HET')."""
    if gt in (".", "", "./.", ".|."):
        return 0, 0, "MISSING"
    alleles = re.split(r"[/|]", gt)
    called = [a for a in alleles if a != "."]
    if not called:
        return len(alleles), 0, "MISSING"
    dosage = sum(a != "0" for a in called)
    if dosage == 0:
        zyg = "HOM_REF"
    elif dosage == len(called):
        zyg = "HEMI_ALT" if len(alleles) == 1 else "HOM_ALT"
    else:
        zyg = "HET"
    return len(alleles), dosage, zyg


def vaf_from(ad: str, af: str) -> float | None:
    if af not in (None, ".", ""):
        try:
            return float(af.split(",")[0])
        except ValueError:
            pass
    if ad not in (None, ".", ""):
        try:
            v = [int(x) for x in ad.split(",") if x != "."]
            if len(v) >= 2 and sum(v) > 0:
                return round(v[1] / sum(v), 4)
        except ValueError:
            pass
    return None


def cn_state(cn: int, base: int, mcn: int | None) -> str:
    if cn == 0:
        return "HOMDEL"
    if cn >= base * 2 + 1:          # high-level amplification wins over LOH
        return "AMP"
    if mcn == 0 and cn >= base:
        return "CN_LOH" if cn == base else "GAIN_LOH"
    if cn < base:
        return "LOSS"
    if cn > base:
        return "GAIN"
    return "NEUTRAL"


# ----------------------------------------------------------------------------- analysis

SITE_INFO = ["END", "SVTYPE", "SVLEN", "GENE", "RU", "SOMATIC"]
SAMPLE_FMT = ["GT", "AD", "DP", "AF", "CN", "MCN"]


def prepare(bt: Bcftools, vcf: str, work: Path, args) -> Path:
    """Filter (PASS / -i expression), split multi-allelics, optionally check REF, write indexed BCF."""
    out = work / "normalized.bcf"
    view = ["view", "-Ou", "--threads", str(bt.threads)]
    if args.pass_only:
        view += ["-f", "PASS,."]
    if args.include:
        view += ["-i", args.include]
    if args.regions:
        view += ["-r", args.regions]
        vcf = _indexed(bt, vcf, work)
    norm = ["norm", "-m", "-any", "-Ob", "-o", str(out), "--threads", str(bt.threads)]
    if args.ref:
        norm += ["-f", args.ref, "--check-ref", "w"]  # warn about REF mismatches, keep going
    bt.pipe(view + [vcf], norm + ["-"])
    bt.run(["index", "-f", str(out)])
    return out


def _indexed(bt: Bcftools, vcf: str, work: Path) -> str:
    """Region queries need an index; if the input has none, make an indexed BCF copy in the work dir."""
    if vcf.endswith((".gz", ".bcf")) and any(Path(vcf + s).exists() for s in (".tbi", ".csi")):
        return vcf
    copy = work / "input.bcf"
    bt.run(["view", "-Ob", "-o", str(copy), vcf])
    bt.run(["index", "-f", str(copy)])
    return str(copy)


def query_table(bt: Bcftools, bcf: Path, h: Header) -> tuple[pd.DataFrame, pd.DataFrame]:
    info = [t for t in SITE_INFO if t in h.info]
    site_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + info
    site_fmt = "\t".join(["%CHROM", "%POS", "%ID", "%REF", "%ALT", "%QUAL", "%FILTER"]
                         + [f"%INFO/{t}" for t in info])
    sites = _query(bt, bcf, site_fmt + "\n", site_cols)
    sites.insert(0, "VID", range(len(sites)))

    fmt = [t for t in SAMPLE_FMT if t in h.fmt]
    if not h.samples or not fmt:
        return sites, pd.DataFrame()
    # one line per record x sample; the running record index is re-attached below
    gfmt = "[" + "\t".join(["%CHROM", "%POS", "%REF", "%ALT", "%SAMPLE"] + [f"%{t}" for t in fmt]) + "\n]"
    geno = _query(bt, bcf, gfmt, ["CHROM", "POS", "REF", "ALT", "SAMPLE"] + fmt)
    geno.insert(0, "VID", [i // len(h.samples) for i in range(len(geno))])
    return sites, geno


def _query(bt: Bcftools, bcf: Path, fmt: str, cols: list[str]) -> pd.DataFrame:
    out = bt.run(["query", "-f", fmt, str(bcf)])
    rows = [r.split("\t") for r in out.splitlines() if r]
    df = pd.DataFrame(rows, columns=cols)
    if not df.empty:
        df["POS"] = df["POS"].astype(int)
    return df


def annotate_sites(sites: pd.DataFrame, h: Header) -> pd.DataFrame:
    if sites.empty:
        return sites
    g = lambda r, k: r.get(k, ".")
    cls, size, scale, rep = [], [], [], []
    for _, r in sites.iterrows():
        c = classify(r.REF, r.ALT, None if g(r, "SVTYPE") == "." else g(r, "SVTYPE"))
        s = event_size(c, r.POS, r.REF, r.ALT, g(r, "END"), g(r, "SVLEN"))
        cls.append(c)
        size.append(s)
        scale.append("SMALL" if c in {"SNV", "MNV"} else "." if c in {"STR", "REF_BLOCK", "NO_ALT"}
                     else event_scale(s, h.contigs.get(r.CHROM, 0)))
        rep.append(repeat_flag(c, r.REF, r.ALT, g(r, "RU")))
    sites["CLASS"], sites["SIZE"], sites["SCALE"], sites["REPEAT"] = cls, size, scale, rep
    return sites


def annotate_genotypes(geno: pd.DataFrame, sites: pd.DataFrame, base_ploidy: int) -> pd.DataFrame:
    if geno.empty:
        return geno
    geno = geno.merge(sites[["VID", "CLASS", "SIZE", "SCALE"] + [c for c in ("END", "GENE") if c in sites]],
                      on="VID", how="left")
    gts = geno["GT"] if "GT" in geno else pd.Series(["."] * len(geno))
    parsed = [parse_gt(x) for x in gts]
    geno["PLOIDY"] = [p for p, _, _ in parsed]
    geno["PLOIDY_NAME"] = [PLOIDY_NAME.get(p, f"{p}-ploid") for p, _, _ in parsed]
    geno["ALT_DOSAGE"] = [d for _, d, _ in parsed]
    geno["ZYGOSITY"] = [z for _, _, z in parsed]
    geno["VAF"] = [vaf_from(geno.at[i, "AD"] if "AD" in geno else None,
                            geno.at[i, "AF"] if "AF" in geno else None) for i in geno.index]
    if "CN" in geno:
        mcn = geno["MCN"] if "MCN" in geno else pd.Series(["."] * len(geno), index=geno.index)
        geno["CN_STATE"] = [
            cn_state(int(c), base_ploidy, None if m in (".", "") else int(m)) if c not in (".", "") else "."
            for c, m in zip(geno["CN"], mcn)]
    return geno


def ploidy_profile(geno: pd.DataFrame) -> pd.DataFrame:
    if geno.empty or "GT" not in geno:
        return pd.DataFrame()
    called = geno[geno["PLOIDY"] > 0]
    tab = called.groupby(["SAMPLE", "PLOIDY_NAME"]).size().unstack(fill_value=0)
    tab["dominant_ploidy"] = tab.idxmax(axis=1)
    return tab.reset_index()


def cn_events(geno: pd.DataFrame) -> pd.DataFrame:
    if geno.empty or "CN_STATE" not in geno:
        return pd.DataFrame()
    ev = geno[geno["CN_STATE"].isin(["HOMDEL", "LOSS", "GAIN", "AMP", "CN_LOH", "GAIN_LOH"])].copy()
    ev["END"] = pd.to_numeric(ev.get("END", ev["POS"]), errors="coerce").fillna(ev["POS"]).astype(int)
    cols = ["SAMPLE", "CHROM", "POS", "END", "SIZE", "SCALE", "CN", "CN_STATE"] + \
           [c for c in ("MCN", "GENE") if c in ev]
    return ev[cols].sort_values(["SAMPLE", "CHROM", "POS"])


def aneuploidy_calls(ev: pd.DataFrame, h: Header) -> pd.DataFrame:
    """Fraction of each chromosome gained / lost per sample; >=90% -> whole-chromosome call."""
    if ev.empty:
        return pd.DataFrame()
    rows = []
    for (s, chrom), d in ev.groupby(["SAMPLE", "CHROM"]):
        L = h.contigs.get(chrom, 0)
        if not L:
            continue
        frac = lambda states: d[d.CN_STATE.isin(states)].eval("END - POS + 1").sum() / L
        gain, loss = frac(["GAIN", "AMP", "GAIN_LOH"]), frac(["LOSS", "HOMDEL"])
        call = ("WHOLE_CHROMOSOME_GAIN" if gain >= 0.9 else "WHOLE_CHROMOSOME_LOSS" if loss >= 0.9
                else "ARM_LEVEL_GAIN" if gain >= 0.25 else "ARM_LEVEL_LOSS" if loss >= 0.25
                else "FOCAL_ONLY")
        rows.append({"SAMPLE": s, "CHROM": chrom, "frac_gained": round(gain, 3),
                     "frac_lost": round(loss, 3), "call": call})
    return pd.DataFrame(rows)


def somatic_calls(geno: pd.DataFrame, h: Header, min_vaf: float, max_normal_vaf: float) -> pd.DataFrame:
    if geno.empty or not (h.tumor and h.normal) or geno["VAF"].isna().all():
        return pd.DataFrame()
    keep = ["VID", "CHROM", "POS", "REF", "ALT", "CLASS"] + (["GENE"] if "GENE" in geno else [])
    t = geno[geno.SAMPLE == h.tumor][keep + ["VAF"]].rename(columns={"VAF": "TUMOR_VAF"})
    n = geno[geno.SAMPLE == h.normal][["VID", "VAF"]].rename(columns={"VAF": "NORMAL_VAF"})
    m = t.merge(n, on="VID")
    m["SOMATIC_CALL"] = (m.TUMOR_VAF >= min_vaf) & (m.NORMAL_VAF.fillna(0) <= max_normal_vaf)
    return m.drop(columns="VID")


# ----------------------------------------------------------------------------- driver

def analyse(bt: Bcftools, vcf: str, outdir: Path, args) -> dict:
    name = re.sub(r"\.(vcf|bcf)(\.gz)?$", "", Path(vcf).name)
    work = outdir / name
    work.mkdir(parents=True, exist_ok=True)

    h = read_header(bt, vcf)
    bcf = prepare(bt, vcf, work, args)
    stats_txt = bt.run(["stats", "-s", "-", str(bcf)])
    (work / "bcftools_stats.txt").write_text(stats_txt)
    stats = parse_stats(stats_txt)

    sites, geno = query_table(bt, bcf, h)
    sites = annotate_sites(sites, h)
    geno = annotate_genotypes(geno, sites, args.baseline_ploidy)
    ploidy = ploidy_profile(geno)
    ev = cn_events(geno)
    aneu = aneuploidy_calls(ev, h)
    som = somatic_calls(geno, h, args.min_vaf, args.max_normal_vaf)

    for df, fn in [(sites, "records"), (geno, "genotypes"), (ploidy, "ploidy"),
                   (ev, "cn_events"), (aneu, "aneuploidy"), (som, "somatic")]:
        if not df.empty:
            df.to_csv(work / f"{fn}.tsv", sep="\t", index=False)

    summary = {
        "file": vcf, "fileformat": h.fileformat, "samples": h.samples,
        "tumor_normal": [h.tumor, h.normal] if h.tumor and h.normal else None,
        "records_after_filters": int(len(sites)),
        "variant_classes": sites["CLASS"].value_counts().to_dict() if not sites.empty else {},
        "event_scales": sites.loc[~sites.CLASS.isin(["SNV", "MNV"]), "SCALE"].value_counts().to_dict()
        if not sites.empty else {},
        "msi_slippage_indels": int((sites["REPEAT"] == "MSI_SLIPPAGE").sum()) if not sites.empty else 0,
        "repeat_expansions": int((sites["REPEAT"] == "REPEAT_EXPANSION").sum()) if not sites.empty else 0,
        "ts_tv": stats["TSTV"].get("ts/tv") if stats["TSTV"].get("ts", 0) + stats["TSTV"].get("tv", 0) else None,
        "bcftools_SN": stats["SN"],
        "dominant_ploidy": dict(zip(ploidy.SAMPLE, ploidy.dominant_ploidy)) if not ploidy.empty else {},
        "aneuploid_chromosomes": aneu[aneu.call.str.startswith("WHOLE")].apply(
            lambda r: f"{r.SAMPLE}:{r.CHROM}:{r.call}", axis=1).tolist() if not aneu.empty else [],
        "arm_level_events": aneu[aneu.call.str.startswith("ARM")].apply(
            lambda r: f"{r.SAMPLE}:{r.CHROM}:{r.call}", axis=1).tolist() if not aneu.empty else [],
        "somatic_calls": int(som.SOMATIC_CALL.sum()) if not som.empty else None,
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return {"summary": summary, "sites": sites, "cn": ev, "contigs": h.contigs}


# ----------------------------------------------------------------------------- plots

INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e0", "#fcfcfb"
BLUE, ORANGE, NEUTRAL = "#2a78d6", "#eb6834", "#9a9990"


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)


def plot_classes(results: list[dict], out: Path):
    import matplotlib.pyplot as plt
    counts = pd.concat([r["sites"]["CLASS"] for r in results if not r["sites"].empty]).value_counts()
    if counts.empty:
        return
    counts = counts.sort_values()
    fig, ax = plt.subplots(figsize=(7, 0.38 * len(counts) + 1.2), facecolor=SURFACE)
    _style(ax)
    ax.barh(counts.index, counts.values, color=BLUE, height=0.6)
    for y, v in enumerate(counts.values):
        ax.text(v, y, f" {v}", va="center", color=INK, fontsize=9)
    ax.set_xlabel("Records (after filtering and splitting)", color=INK2)
    ax.set_title("Variant classes across all input files", loc="left", color=INK, fontsize=12)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_copy_number(results: list[dict], out: Path, base: int):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    ev = [(r["cn"], r["contigs"]) for r in results if not r["cn"].empty]
    if not ev:
        return
    contigs = {}
    for _, c in ev:
        contigs.update(c)
    order = sorted(contigs, key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c))
    offset, acc = {}, 0
    for c in order:
        offset[c], acc = acc, acc + contigs[c]
    samples = [(s, d) for df, _ in ev for s, d in df.groupby("SAMPLE")]
    fig, axes = plt.subplots(len(samples), 1, figsize=(10, 1.6 * len(samples) + 0.6),
                             sharex=True, squeeze=False, facecolor=SURFACE)
    ymax = max(pd.to_numeric(d.CN).max() for _, d in samples) + 1
    for ax, (s, d) in zip(axes[:, 0], samples):
        _style(ax)
        for c in order:
            ax.axvline(offset[c], color=GRID, lw=0.6)
        ax.axhline(base, color=NEUTRAL, lw=1, ls="--")
        for _, r in d.iterrows():
            cn = int(r.CN)
            col = BLUE if cn < base else ORANGE if cn > base else NEUTRAL
            x0 = offset[r.CHROM] + r.POS
            x1 = offset[r.CHROM] + max(r.END, r.POS + contigs[r.CHROM] * 0.004)  # keep focal events visible
            ax.plot([x0, x1], [cn, cn], color=col, lw=4, solid_capstyle="round")
        ax.set_ylim(-0.5, ymax)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylabel("CN", color=INK2, fontsize=9)
        ax.set_title(s, loc="left", color=INK, fontsize=10)
    axes[-1, 0].set_xticks([offset[c] + contigs[c] / 2 for c in order])
    axes[-1, 0].set_xticklabels([c.replace("chr", "") for c in order])
    axes[-1, 0].set_xlabel("Chromosome", color=INK2)
    handles = [Line2D([], [], color=BLUE, lw=4, label=f"Loss (CN < {base})"),
               Line2D([], [], color=ORANGE, lw=4, label=f"Gain (CN > {base})"),
               Line2D([], [], color=NEUTRAL, lw=4, label="Copy-neutral (e.g. LOH)"),
               Line2D([], [], color=NEUTRAL, lw=1, ls="--", label=f"Baseline CN {base}")]
    fig.legend(handles=handles, loc="upper right", ncol=4, frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle("Copy-number events by sample", x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def run(vcfs:list, outdir:Path, bcftools:str, ref:Path|str|None=None, include:str|None=None,
        regions:Path|str|None=None, pass_only:bool=False, baseline_ploidy:int=2,
        min_vaf:float=0.05, max_normal_vaf:float=0.02, threads:int=2, no_plots:bool=False):
    
    args = SimpleNamespace(
        vcf=[str(v) for v in vcfs], outdir=outdir, bcftools=bcftools, ref=ref,
        include=include, regions=regions, pass_only=pass_only,
        baseline_ploidy=baseline_ploidy, min_vaf=min_vaf,
        max_normal_vaf=max_normal_vaf, threads=threads, no_plots=no_plots,
    )

    bt = Bcftools(args.bcftools, args.threads)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    results = []
    for vcf in args.vcf:
        try:
            results.append(analyse(bt, vcf, out, args))
        except BcftoolsError as e:
            print(f"[skip] {vcf}: {e}")

    if results and not args.no_plots:
        plot_classes(results, out / "variant_classes.png")
        plot_copy_number(results, out / "copy_number.png", args.baseline_ploidy)
    return results

