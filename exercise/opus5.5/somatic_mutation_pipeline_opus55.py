#!/usr/bin/env python3
"""Multi-region somatic WGS pipeline (Python port of pipeline.sh).

Per patient: tumour (<PT>_T), border (<PT>_B), non-tumour tissue (<PT>_N),
plus GTEx donors (GTEX-*) as panel-of-normals / population controls.
Reference: GRCh38, Broad resource bundle file names.

Default is a DRY RUN: every command is printed, nothing is executed and no
directory is created. Add --run to execute.

    python pipeline.py                         # print all commands
    python pipeline.py --steps align,markdup   # print a subset
    python pipeline.py --patients P01 --run    # execute for one patient

Production alternative for steps align..annotate: nf-core/sarek
(--tools haplotypecaller,mutect2,strelka,manta,ascat,msisensorpro,vep).
"""
from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
THREADS = 16
MEM = "32g"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

R = Path("ref")
REF = R / "Homo_sapiens_assembly38.fasta"   # + .fai, .dict, bwa-mem2 index, .alt
DBSNP = R / "Homo_sapiens_assembly38.dbsnp138.vcf"
MILLS = R / "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KINDELS = R / "Homo_sapiens_assembly38.known_indels.vcf.gz"
GNOMAD = R / "af-only-gnomad.hg38.vcf.gz"
COMMON = R / "small_exac_common_3.hg38.vcf.gz"
INTERVALS = R / "wgs_calling_regions.hg38.interval_list"
VBID_SVD = R / "verifybamid/1000g.phase3.100k.b38.vcf.gz.dat"
SOMALIER_SITES = R / "somalier/sites.hg38.vcf.gz"
DELLY_EXCL = R / "delly/human.hg38.excl.tsv"
MSI_LIST = R / "hg38.msisensor.list"   # once: msisensor-pro scan -d REF -o MSI_LIST
ASCAT_REF = R / "ascat"
VEP_CACHE = R / "vep_cache"
MELT_JAR = Path("MELT.jar")

# mitochondria (Broad hg38/v0/chrM bundle; GRCh38 chrM = rCRS)
MT = R / "chrM"
MT_REF = MT / "Homo_sapiens_assembly38.chrM.fasta"                        # + bwa index, .fai, .dict
MT_SHIFT_REF = MT / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
MT_CHAIN = MT / "ShiftBack.chain"
MT_BLACKLIST = MT / "blacklist_sites.hg38.chrM.bed"
MT_NONCONTROL = "chrM:576-16024"     # called on the standard chrM
MT_CONTROL_SHIFTED = "chrM:8025-9144"  # control region (D-loop) called on shifted chrM
MT_MIN_AF = 0.01                     # heteroplasmy reporting threshold
MT_HOMOPLASMIC = 0.95

OUT_DIRS = ["qc/somalier", "bam", "vcf/pon", "sv", "cnv", "te", "clonality", "sig/vcfs",
            "tmp", "mito"]


# ---------------------------------------------------------------------------
# sample sheet
# ---------------------------------------------------------------------------
@dataclass
class Lane:
    lane: str
    fq1: str
    fq2: str


@dataclass
class Sample:
    patient: str
    sex: str
    status: int
    name: str
    lanes: list[Lane] = field(default_factory=list)

    @property
    def is_gtex(self) -> bool:
        return self.patient.startswith("GTEX")


def read_sheet(path: Path) -> dict[str, Sample]:
    samples: dict[str, Sample] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            s = samples.setdefault(
                row["sample"],
                Sample(row["patient"], row["sex"], int(row["status"]), row["sample"]),
            )
            s.lanes.append(Lane(row["lane"], row["fastq_1"], row["fastq_2"]))
    return samples


# ---------------------------------------------------------------------------
# command runner
# ---------------------------------------------------------------------------
class Runner:
    def __init__(self, execute: bool):
        self.execute = execute

    @staticmethod
    def _fmt(cmd) -> str:
        return " ".join(shlex.quote(str(c)) for c in cmd)

    def run(self, *cmds, stdin: str | None = None, stdout: Path | None = None):
        """Run one command, or several joined by pipes (cmd1 | cmd2 | ...)."""
        cmds = [[str(c) for c in cmd] for cmd in cmds]
        text = " | ".join(self._fmt(c) for c in cmds)
        if stdin is not None:
            text += f" <<'EOF'\n{stdin.rstrip()}\nEOF"
        if stdout is not None:
            text += f" > {stdout}"
        print(text + "\n", flush=True)
        if not self.execute:
            return

        out_fh = open(stdout, "w") if stdout else None
        procs, prev = [], None
        for i, cmd in enumerate(cmds):
            last = i == len(cmds) - 1
            p = subprocess.Popen(
                cmd,
                stdin=prev.stdout if prev else (subprocess.PIPE if stdin else None),
                stdout=(out_fh if last else subprocess.PIPE),
                text=True,
            )
            if prev:
                prev.stdout.close()          # let upstream get SIGPIPE
            procs.append(p)
            prev = p
        if stdin is not None:
            procs[0].stdin.write(stdin)
            procs[0].stdin.close()
        codes = [p.wait() for p in procs]
        if out_fh:
            out_fh.close()
        if any(codes):
            sys.exit(f"FAILED ({codes}): {text}")

    def write(self, path: Path, content: str):
        print(f"# write {path}:\n{content}", flush=True)
        if self.execute:
            path.write_text(content)

    def mkdir(self, *paths):
        if self.execute:
            for p in paths:
                Path(p).mkdir(parents=True, exist_ok=True)


def gatk(tool: str, *args) -> list[str]:
    return ["gatk", "--java-options", f"-Xmx{MEM}", tool, *map(str, args)]


def bqsr_bam(sm: str) -> Path:
    return Path("bam") / f"{sm}.bqsr.bam"


# ---------------------------------------------------------------------------
# 0. GTEx CRAM -> FASTQ (realign everything with the same aligner/version)
# ---------------------------------------------------------------------------
def gtex_to_fastq(r: Runner, donor: str, cram: Path):
    r.run(
        ["samtools", "collate", "-@", THREADS, "-Ou", "--reference", REF, cram, f"tmp/{donor}"],
        ["samtools", "fastq", "-@", THREADS, "-n",
         "-1", f"fastq/{donor}_L001_R1.fastq.gz", "-2", f"fastq/{donor}_L001_R2.fastq.gz",
         "-0", "/dev/null", "-s", "/dev/null", "-"],
    )


# ---------------------------------------------------------------------------
# 1-2. FASTQ QC + alignment per lane (no trimming: BWA soft-clips adapters)
# ---------------------------------------------------------------------------
def align(r: Runner, s: Sample):
    for ln in s.lanes:
        rg = (f"@RG\\tID:{s.name}.{ln.lane}\\tSM:{s.name}\\tLB:{s.name}_lib1"
              f"\\tPL:ILLUMINA\\tPU:{s.name}.{ln.lane}")
        r.run(["fastqc", "-t", 4, "-o", "qc", ln.fq1, ln.fq2])
        r.run(
            ["bwa-mem2", "mem", "-t", THREADS, "-Y", "-K", 100000000, "-R", rg, REF, ln.fq1, ln.fq2],
            ["samtools", "sort", "-@", 4, "-m", "2G", "-o", f"bam/{s.name}.{ln.lane}.bam", "-"],
        )


# ---------------------------------------------------------------------------
# 3. merge lanes + MarkDuplicates + BQSR
# ---------------------------------------------------------------------------
def markdup_bqsr(r: Runner, s: Sample):
    sm = s.name
    ins = [x for ln in s.lanes for x in ("-I", f"bam/{sm}.{ln.lane}.bam")]
    r.run(gatk("MarkDuplicates", *ins, "-O", f"bam/{sm}.md.bam",
               "-M", f"qc/{sm}.md_metrics.txt",
               "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", 2500, "--CREATE_INDEX", "true"))
    r.run(gatk("BaseRecalibrator", "-R", REF, "-I", f"bam/{sm}.md.bam",
               "--known-sites", DBSNP, "--known-sites", MILLS, "--known-sites", KINDELS,
               "-O", f"bam/{sm}.recal.table"))
    r.run(gatk("ApplyBQSR", "-R", REF, "-I", f"bam/{sm}.md.bam",
               "--bqsr-recal-file", f"bam/{sm}.recal.table", "-O", bqsr_bam(sm)))


# ---------------------------------------------------------------------------
# 4. BAM QC: depth, WGS metrics, insert size, contamination, identity
# ---------------------------------------------------------------------------
def bam_qc(r: Runner, s: Sample):
    sm, b = s.name, bqsr_bam(s.name)
    r.run(["mosdepth", "-t", 4, "-n", "--fast-mode", "--by", 1000000, f"qc/{sm}", b])
    r.run(gatk("CollectWgsMetrics", "-R", REF, "-I", b, "-O", f"qc/{sm}.wgs_metrics.txt"))
    r.run(gatk("CollectInsertSizeMetrics", "-I", b, "-O", f"qc/{sm}.insert.txt",
               "-H", f"qc/{sm}.insert.pdf"))
    r.run(["verifybamid2", "--SVDPrefix", VBID_SVD, "--Reference", REF,
           "--BamFile", b, "--Output", f"qc/{sm}.vbid2"])
    r.run(["somalier", "extract", "-d", "qc/somalier/", "--sites", SOMALIER_SITES, "-f", REF, b])


def identity_check(r: Runner, samples: list[Sample]):
    """T/B/N of a patient must show relatedness ~1; GTEx donors unrelated."""
    r.run(["somalier", "relate", "--infer", "-o", "qc/somalier_relate",
           *[f"qc/somalier/{s.name}.somalier" for s in samples]])


# ---------------------------------------------------------------------------
# 5. germline (N + GTEx)
# ---------------------------------------------------------------------------
def germline_gvcf(r: Runner, sm: str):
    r.run(gatk("HaplotypeCaller", "-R", REF, "-I", bqsr_bam(sm), "-L", INTERVALS,
               "-ERC", "GVCF", "-O", f"vcf/{sm}.g.vcf.gz"))


def germline_joint(r: Runner, names: list[str]):
    vs = [x for n in names for x in ("-V", f"vcf/{n}.g.vcf.gz")]
    for c in CHROMS:
        r.run(gatk("GenomicsDBImport", *vs, "-L", c, "--genomicsdb-workspace-path", f"vcf/gdb_{c}"))
        r.run(gatk("GenotypeGVCFs", "-R", REF, "-V", f"gendb://vcf/gdb_{c}",
                   "-O", f"vcf/germline.{c}.vcf.gz"))
    # then: MergeVcfs -> VariantRecalibrator/ApplyVQSR (WGS cohort) or hard filters


# ---------------------------------------------------------------------------
# 6. panel of normals from GTEx only (adjacent N may carry tumour hotspots)
# ---------------------------------------------------------------------------
def pon_normal(r: Runner, sm: str):
    r.run(gatk("Mutect2", "-R", REF, "-I", bqsr_bam(sm), "-L", INTERVALS,
               "--max-mnp-distance", 0, "-O", f"vcf/pon/{sm}.vcf.gz"))


def pon_build(r: Runner, gtex: list[str]):
    vs = [x for g in gtex for x in ("-V", f"vcf/pon/{g}.vcf.gz")]
    r.run(gatk("GenomicsDBImport", "-R", REF, "-L", INTERVALS,
               "--genomicsdb-workspace-path", "vcf/pon_db", *vs))
    r.run(gatk("CreateSomaticPanelOfNormals", "-R", REF, "-V", "gendb://vcf/pon_db",
               "--germline-resource", GNOMAD, "-O", "vcf/pon.vcf.gz"))


# ---------------------------------------------------------------------------
# 7. somatic SNV/indel: joint Mutect2 (T + B vs N), scattered by chromosome
# ---------------------------------------------------------------------------
def mutect2_multi(r: Runner, pt: str):
    d = Path("vcf") / pt
    r.mkdir(d)
    T, B, N = (bqsr_bam(f"{pt}_{x}") for x in "TBN")
    for c in CHROMS:
        r.run(gatk("Mutect2", "-R", REF, "-I", T, "-I", B, "-I", N, "-normal", f"{pt}_N",
                   "--germline-resource", GNOMAD, "--panel-of-normals", "vcf/pon.vcf.gz",
                   "-L", c, "--f1r2-tar-gz", d / f"f1r2.{c}.tar.gz",
                   "-O", d / f"unf.{c}.vcf.gz"))

    vs = [x for c in CHROMS for x in ("-I", d / f"unf.{c}.vcf.gz")]
    ss = [x for c in CHROMS for x in ("-stats", d / f"unf.{c}.vcf.gz.stats")]
    fs = [x for c in CHROMS for x in ("-I", d / f"f1r2.{c}.tar.gz")]
    r.run(gatk("MergeVcfs", *vs, "-O", d / "unfiltered.vcf.gz"))
    r.run(gatk("MergeMutectStats", *ss, "-O", d / "unfiltered.vcf.gz.stats"))
    r.run(gatk("LearnReadOrientationModel", *fs, "-O", d / "rom.tar.gz"))

    for x in "TBN":
        r.run(gatk("GetPileupSummaries", "-I", bqsr_bam(f"{pt}_{x}"), "-V", COMMON,
                   "-L", COMMON, "-O", d / f"{x}.pileups.table"))
    for x in "TB":
        r.run(gatk("CalculateContamination", "-I", d / f"{x}.pileups.table",
                   "-matched", d / "N.pileups.table",
                   "--tumor-segmentation", d / f"{x}.segments.table",
                   "-O", d / f"{x}.contamination.table"))

    r.run(gatk("FilterMutectCalls", "-R", REF, "-V", d / "unfiltered.vcf.gz",
               "--contamination-table", d / "T.contamination.table",
               "--contamination-table", d / "B.contamination.table",
               "--tumor-segmentation", d / "T.segments.table",
               "--tumor-segmentation", d / "B.segments.table",
               "--ob-priors", d / "rom.tar.gz", "-O", d / "mutect2.filtered.vcf.gz"))
    r.run(["bcftools", "view", "-f", "PASS", d / "mutect2.filtered.vcf.gz",
           "-Oz", "-o", d / "mutect2.pass.vcf.gz"])
    r.run(["bcftools", "index", "-t", d / "mutect2.pass.vcf.gz"])


# ---------------------------------------------------------------------------
# 7b. mitochondrial genome (GATK mitochondria pipeline, adapted)
#   1. take reads that mapped to chrM in the whole-genome BAM (NuMT-derived
#      reads mostly stay on their nuclear loci)
#   2. realign them to chrM and to chrM shifted by 8 kb, so the control region
#      (D-loop, which spans the chrM start/end junction) is called away from
#      the linear-reference edge
#   3. Mutect2 --mitochondria-mode on each; lift the shifted calls back; merge
#   4. filter twice: pass 1 without contamination -> haplocheck -> pass 2 with it
#   5. blacklist mask, split multiallelics, haplogroup (Haplogrep3)
#   Simplification vs the WDL: FASTQ realignment instead of RevertSam +
#   MergeBamAlignment (original read-group and tag metadata not carried over).
# ---------------------------------------------------------------------------
def mt_dir(sm: str) -> Path:
    return Path("mito") / sm


def mosdepth_means(sm: str) -> tuple[float, float]:
    """(autosomal mean depth, chrM mean depth) from mosdepth summary."""
    auto_bases = auto_len = 0.0
    mt = None
    with open(f"qc/{sm}.mosdepth.summary.txt") as fh:
        next(fh)
        for line in fh:
            chrom, length, bases, mean, *_ = line.split("\t")
            if chrom in {f"chr{i}" for i in range(1, 23)}:
                auto_len += float(length)
                auto_bases += float(bases)
            elif chrom == "chrM":
                mt = float(mean)
    if mt is None or auto_len == 0:
        sys.exit(f"chrM/autosomes missing in qc/{sm}.mosdepth.summary.txt")
    return auto_bases / auto_len, mt


def parse_haplocheck(path: Path) -> float:
    """Contamination level from haplocheck output ('ND' / absent -> 0)."""
    with open(path, newline="") as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    key = next((k for k in row if "contamination level" in k.lower().strip('"')), None)
    try:
        return float(row[key].strip('"')) if key else 0.0
    except ValueError:
        return 0.0


def mito_extract_align(r: Runner, sm: str):
    d = mt_dir(sm)
    r.mkdir(d)
    r.run(gatk("PrintReads", "-R", REF, "-I", f"bam/{sm}.md.bam", "-L", "chrM",
               "--read-filter", "MateOnSameContigOrNoMappedMateReadFilter",
               "--read-filter", "MateUnmappedAndUnmappedReadFilter",
               "-O", d / "chrM.subset.bam"))
    r.run(["samtools", "collate", "-Ou", d / "chrM.subset.bam", f"tmp/{sm}.mt"],
          ["samtools", "fastq", "-n", "-0", "/dev/null", "-s", "/dev/null",
           "-o", d / "chrM.interleaved.fq.gz", "-"])
    rg = f"@RG\\tID:{sm}.chrM\\tSM:{sm}\\tLB:{sm}_lib1\\tPL:ILLUMINA"
    for tag, ref in (("chrM", MT_REF), ("shifted", MT_SHIFT_REF)):
        r.run(["bwa", "mem", "-K", 100000000, "-p", "-v", 3, "-t", 4, "-Y", "-R", rg,
               ref, d / "chrM.interleaved.fq.gz"],
              ["samtools", "sort", "-o", d / f"{tag}.bam", "-"])
        r.run(gatk("MarkDuplicates", "-I", d / f"{tag}.bam", "-O", d / f"{tag}.md.bam",
                   "-M", d / f"{tag}.md_metrics.txt", "--CREATE_INDEX", "true"))


def mito_call(r: Runner, sm: str):
    d = mt_dir(sm)
    common = ["--mitochondria-mode", "--annotation", "StrandBiasBySample",
              "--max-reads-per-alignment-start", 75, "--max-mnp-distance", 0,
              "--read-filter", "MateOnSameContigOrNoMappedMateReadFilter",
              "--read-filter", "MateUnmappedAndUnmappedReadFilter"]
    r.run(gatk("Mutect2", "-R", MT_REF, "-I", d / "chrM.md.bam", "-L", MT_NONCONTROL,
               *common, "-O", d / "chrM.vcf.gz"))
    r.run(gatk("Mutect2", "-R", MT_SHIFT_REF, "-I", d / "shifted.md.bam",
               "-L", MT_CONTROL_SHIFTED, *common, "-O", d / "shifted.vcf.gz"))
    r.run(gatk("LiftoverVcf", "-I", d / "shifted.vcf.gz", "-O", d / "shifted.lifted.vcf.gz",
               "-R", MT_REF, "--CHAIN", MT_CHAIN, "--REJECT", d / "shifted.rejected.vcf.gz"))
    r.run(gatk("MergeVcfs", "-I", d / "chrM.vcf.gz", "-I", d / "shifted.lifted.vcf.gz",
               "-O", d / "raw.vcf.gz"))
    r.run(gatk("MergeMutectStats", "-stats", d / "chrM.vcf.gz.stats",
               "-stats", d / "shifted.vcf.gz.stats", "-O", d / "raw.vcf.gz.stats"))


def mito_filter(r: Runner, sm: str, tag: str, contamination, autosomal_cov):
    d = mt_dir(sm)
    r.run(gatk("FilterMutectCalls", "-R", MT_REF, "-V", d / "raw.vcf.gz",
               "--stats", d / "raw.vcf.gz.stats", "--mitochondria-mode",
               "--max-alt-allele-count", 4, "--min-allele-fraction", 0,
               "--autosomal-coverage", autosomal_cov,
               "--contamination-estimate", contamination,
               "-O", d / f"{tag}.filtered.vcf.gz"))
    r.run(gatk("VariantFiltration", "-R", MT_REF, "-V", d / f"{tag}.filtered.vcf.gz",
               "--apply-allele-specific-filters",
               "--mask", MT_BLACKLIST, "--mask-name", "blacklisted_site",
               "-O", d / f"{tag}.masked.vcf.gz"))
    r.run(gatk("LeftAlignAndTrimVariants", "-R", MT_REF, "-V", d / f"{tag}.masked.vcf.gz",
               "--split-multi-allelics", "--dont-trim-alleles", "--keep-original-ac",
               "-O", d / f"{tag}.split.vcf.gz"))


def mito(r: Runner, sm: str):
    """Full per-sample mtDNA workflow. Needs markdup + qc (mosdepth) done."""
    d = mt_dir(sm)
    mito_extract_align(r, sm)
    mito_call(r, sm)

    if r.execute:
        auto_cov, mt_cov = mosdepth_means(sm)
        r.write(d / "mtcn.tsv", "sample\tautosomal_mean\tchrM_mean\tmtDNA_copies_per_cell\n"
                f"{sm}\t{auto_cov:.2f}\t{mt_cov:.2f}\t{2 * mt_cov / auto_cov:.1f}\n")
    else:
        auto_cov = f"<autosomal mean from qc/{sm}.mosdepth.summary.txt>"
        print(f"# mtDNA copies/cell = 2 * chrM_mean / autosomal_mean -> {d}/mtcn.tsv\n")

    mito_filter(r, sm, "pass1", 0, auto_cov)
    r.run(["haplocheck", "--out", d / "haplocheck.txt", d / "pass1.split.vcf.gz"])
    contamination = (parse_haplocheck(d / "haplocheck.txt") if r.execute
                     else f"<contamination level from {d}/haplocheck.txt>")
    mito_filter(r, sm, "final", contamination, auto_cov)
    r.run(["haplogrep3", "classify", "--in", d / "final.split.vcf.gz",
           "--out", d / "haplogroup.txt", "--tree", "phylotree-rcrs@17.2"])


def read_mt_calls(sm: str) -> dict[tuple, dict]:
    from cyvcf2 import VCF
    calls = {}
    for v in VCF(str(mt_dir(sm) / "final.split.vcf.gz")):
        calls[(v.POS, v.REF, v.ALT[0])] = dict(
            AF=float(v.format("AF")[0][0]), DP=int(v.format("DP")[0][0]),
            FILTER=v.FILTER or "PASS")
    return calls


def mito_gtex_recurrent(r: Runner, gtex: list[str]):
    """Heteroplasmic PASS calls recurring in >=2 unrelated GTEx donors:
    likely NuMT/alignment artefacts or known heteroplasmy hotspots. Flagged, not removed.
    Note: GTEx WGS is blood, so this is a technical background, not tissue biology."""
    print(f"# python: recurrent heteroplasmic sites across GTEx {gtex} -> mito/gtex_recurrent.tsv\n")
    if not r.execute:
        return
    import pandas as pd
    rows = []
    for g in gtex:
        for (pos, ref, alt), c in read_mt_calls(g).items():
            if c["FILTER"] == "PASS" and MT_MIN_AF <= c["AF"] < MT_HOMOPLASMIC:
                rows.append(dict(POS=pos, REF=ref, ALT=alt, donor=g))
    df = pd.DataFrame(rows, columns=["POS", "REF", "ALT", "donor"])
    rec = df.groupby(["POS", "REF", "ALT"]).donor.nunique().rename("n_donors").reset_index()
    rec[rec.n_donors >= 2].to_csv("mito/gtex_recurrent.tsv", sep="\t", index=False)


def mito_compare(r: Runner, pt: str):
    """Per-patient T/B/N heteroplasmy table + mtDNA copy number."""
    out = Path("mito") / f"{pt}.heteroplasmy.tsv"
    print(f"# python: merge T/B/N chrM calls, classify, add mtCN -> {out}\n")
    if not r.execute:
        return
    import pandas as pd
    names = {x: f"{pt}_{x}" for x in "TBN"}
    calls = {x: read_mt_calls(sm) for x, sm in names.items()}
    keys = sorted(set().union(*calls.values()))
    rec_path = Path("mito/gtex_recurrent.tsv")
    recurrent = set()
    if rec_path.exists():
        g = pd.read_csv(rec_path, sep="\t")
        recurrent = set(zip(g.POS, g.REF, g.ALT))

    rows = []
    for k in keys:
        row = dict(POS=k[0], REF=k[1], ALT=k[2])
        for x in "TBN":
            c = calls[x].get(k, dict(AF=0.0, DP=None, FILTER="absent"))
            row.update({f"AF_{x}": c["AF"], f"DP_{x}": c["DP"], f"FILTER_{x}": c["FILTER"]})
        n, t, b = row["AF_N"], row["AF_T"], row["AF_B"]
        pass_tb = "PASS" in (row["FILTER_T"], row["FILTER_B"])
        if min(n, t, b) >= MT_HOMOPLASMIC:
            cls = "germline_homoplasmic"
        elif n < MT_MIN_AF and pass_tb and max(t, b) >= MT_MIN_AF:
            cls = "somatic_candidate"      # absent from N, present in T and/or B
        elif n >= MT_MIN_AF and abs(t - n) >= 0.10:
            cls = "heteroplasmy_shift"     # inherited heteroplasmy drifting in tumour
        else:
            cls = "other"
        row["class"] = cls
        row["gtex_recurrent"] = k in recurrent
        rows.append(row)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)

    cn = pd.concat(pd.read_csv(mt_dir(sm) / "mtcn.tsv", sep="\t") for sm in names.values())
    cn.to_csv(Path("mito") / f"{pt}.mtcn.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# 8. Manta -> Strelka2, Delly, SV consensus (pairwise vs N)
# ---------------------------------------------------------------------------
def manta_strelka(r: Runner, pt: str, s: str):
    d = Path("sv") / f"{pt}_{s}"
    T, N = bqsr_bam(f"{pt}_{s}"), bqsr_bam(f"{pt}_N")
    r.run(["configManta.py", "--normalBam", N, "--tumorBam", T,
           "--referenceFasta", REF, "--runDir", d / "manta"])
    r.run([d / "manta/runWorkflow.py", "-m", "local", "-j", THREADS])
    r.run(["configureStrelkaSomaticWorkflow.py", "--normalBam", N, "--tumorBam", T,
           "--referenceFasta", REF,
           "--indelCandidates", d / "manta/results/variants/candidateSmallIndels.vcf.gz",
           "--runDir", d / "strelka"])
    r.run([d / "strelka/runWorkflow.py", "-m", "local", "-j", THREADS])


def delly_somatic(r: Runner, pt: str, s: str):
    d = Path("sv") / f"{pt}_{s}"
    r.mkdir(d)
    r.run(["delly", "call", "-x", DELLY_EXCL, "-g", REF, "-o", d / "delly.bcf",
           bqsr_bam(f"{pt}_{s}"), bqsr_bam(f"{pt}_N")])
    r.write(d / "delly_samples.tsv", f"{pt}_{s}\ttumor\n{pt}_N\tcontrol\n")
    r.run(["delly", "filter", "-f", "somatic", "-s", d / "delly_samples.tsv",
           "-o", d / "delly.somatic.bcf", d / "delly.bcf"])


def sv_consensus(r: Runner, pt: str, s: str):
    """SVs supported by >= 2 callers, breakpoints within 1 kb, same type, min 50 bp."""
    d = Path("sv") / f"{pt}_{s}"
    r.run(["bcftools", "view", "-f", "PASS", d / "manta/results/variants/somaticSV.vcf.gz",
           "-Ov", "-o", d / "manta.vcf"])
    r.run(["bcftools", "view", "-f", "PASS", d / "delly.somatic.bcf", "-Ov", "-o", d / "delly.vcf"])
    r.write(d / "sv_list.txt", f"{d / 'manta.vcf'}\n{d / 'delly.vcf'}\n")
    r.run(["SURVIVOR", "merge", d / "sv_list.txt", 1000, 2, 1, 1, 0, 50, d / "sv.consensus.vcf"])


# ---------------------------------------------------------------------------
# 9. copy number, purity, ploidy (ASCAT via Rscript)
#    Low-purity border samples may not fit: use --purity in make_pyclone_input.py
# ---------------------------------------------------------------------------
ASCAT_R = r"""
args <- commandArgs(trailingOnly = TRUE)
pt <- args[1]; s <- args[2]; sex <- args[3]; r <- args[4]
library(ASCAT)
tum <- paste0(pt, "_", s); nor <- paste0(pt, "_N")
out <- file.path("cnv", tum); dir.create(out, recursive = TRUE, showWarnings = FALSE)
f <- function(x) file.path(out, x)
ascat.prepareHTS(
  tumourseqfile = file.path("bam", paste0(tum, ".bqsr.bam")),
  normalseqfile = file.path("bam", paste0(nor, ".bqsr.bam")),
  tumourname = tum, normalname = nor, allelecounter_exe = "alleleCounter",
  alleles.prefix = file.path(r, "G1000_alleles_hg38_chr"),
  loci.prefix    = file.path(r, "G1000_loci_hg38_chr"),
  gender = sex, genomeVersion = "hg38", nthreads = 8,
  tumourLogR_file = f("T_LogR.txt"), tumourBAF_file = f("T_BAF.txt"),
  normalLogR_file = f("N_LogR.txt"), normalBAF_file = f("N_BAF.txt"))
bc <- ascat.loadData(Tumor_LogR_file = f("T_LogR.txt"), Tumor_BAF_file = f("T_BAF.txt"),
                     Germline_LogR_file = f("N_LogR.txt"), Germline_BAF_file = f("N_BAF.txt"),
                     gender = sex, genomeVersion = "hg38")
bc <- ascat.correctLogR(bc, GCcontentfile = file.path(r, "GC_G1000_hg38.txt"),
                        replictimingfile = file.path(r, "RT_G1000_hg38.txt"))
bc <- ascat.aspcf(bc)
res <- ascat.runAscat(bc, gamma = 1, write_segments = TRUE)
write.table(res$segments, f("segments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = tum, purity = res$aberrantcellfraction, ploidy = res$ploidy),
            f("purity_ploidy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
"""


def ascat_run(r: Runner, pt: str, s: str, sex: str):
    r.run(["Rscript", "-", pt, s, sex, ASCAT_REF], stdin=ASCAT_R)


# ---------------------------------------------------------------------------
# 10. MSI, annotation
# ---------------------------------------------------------------------------
def msi(r: Runner, pt: str, s: str):
    r.run(["msisensor-pro", "msi", "-d", MSI_LIST, "-n", bqsr_bam(f"{pt}_N"),
           "-t", bqsr_bam(f"{pt}_{s}"), "-o", f"qc/{pt}_{s}.msi", "-b", THREADS])


def annotate(r: Runner, pt: str):
    d = Path("vcf") / pt
    r.run(["vep", "--offline", "--cache", "--dir_cache", VEP_CACHE, "--assembly", "GRCh38",
           "--fasta", REF, "--everything", "--vcf", "--compress_output", "bgzip",
           "--fork", 8, "-i", d / "mutect2.pass.vcf.gz", "-o", d / "mutect2.pass.vep.vcf.gz"])


# ---------------------------------------------------------------------------
# 11. transposable elements
#     Germline MEIs: MELT. Somatic L1/Alu/SVA (T, B vs N): xTea -- its flags
#     change between releases; generate its run script following its README.
# ---------------------------------------------------------------------------
def melt_germline(r: Runner, sm: str):
    r.run(["java", "-Xmx8g", "-jar", MELT_JAR, "Single", "-bamfile", bqsr_bam(sm), "-h", REF,
           "-t", R / "melt/mei_list.txt", "-n", R / "melt/hg38.genes.bed", "-w", f"te/{sm}"])


# ---------------------------------------------------------------------------
# 12. clonality: CCF clusters across T and B (PyClone-VI)
# ---------------------------------------------------------------------------
def clonality(r: Runner, pt: str):
    r.run([sys.executable, "make_pyclone_input.py",
           "--vcf", f"vcf/{pt}/mutect2.pass.vcf.gz",
           "--sample", f"{pt}_T=cnv/{pt}_T", "--sample", f"{pt}_B=cnv/{pt}_B",
           "-o", f"clonality/{pt}.pyclone_in.tsv"])
    r.run(["pyclone-vi", "fit", "-i", f"clonality/{pt}.pyclone_in.tsv",
           "-o", f"clonality/{pt}.h5", "-c", 40, "-d", "beta-binomial", "-r", 10])
    r.run(["pyclone-vi", "write-results-file", "-i", f"clonality/{pt}.h5",
           "-o", f"clonality/{pt}.pyclone_out.tsv"])


# ---------------------------------------------------------------------------
# 13. mutational signatures (needs SigProfilerMatrixGenerator GRCh38 installed)
# ---------------------------------------------------------------------------
def signatures(r: Runner, patients: list[str]):
    for pt in patients:
        for s in "TB":
            r.run(["bcftools", "view", "-s", f"{pt}_{s}", f"vcf/{pt}/mutect2.pass.vcf.gz"],
                  ["bcftools", "view", "-i", "FMT/AD[0:1]>=3", "-Ov",
                   "-o", f"sig/vcfs/{pt}_{s}.vcf"])
    print("# SigProfilerAssignment.cosmic_fit(samples='sig/vcfs', output='sig/out', "
          "input_type='vcf', genome_build='GRCh38')\n")
    if r.execute:
        from SigProfilerAssignment import Analyzer as Analyze
        Analyze.cosmic_fit(samples="sig/vcfs", output="sig/out",
                           input_type="vcf", genome_build="GRCh38")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
STEPS = ["align", "markdup", "qc", "germline", "pon", "somatic", "mito",
         "sv", "cnv", "msi", "annotate", "clonality", "signatures"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", type=Path, default=Path("samplesheet.csv"))
    ap.add_argument("--run", action="store_true", help="execute (default: dry run)")
    ap.add_argument("--steps", default=",".join(STEPS),
                    help=f"comma-separated subset of: {','.join(STEPS)}")
    ap.add_argument("--patients", help="comma-separated subset of patients (GTEx always kept)")
    a = ap.parse_args()

    steps = set(a.steps.split(","))
    if unknown := steps - set(STEPS):
        ap.error(f"unknown steps: {sorted(unknown)}")

    r = Runner(a.run)
    if not a.run:
        print("# DRY RUN: commands are printed, nothing is executed. Add --run to execute.\n")
    r.mkdir(*OUT_DIRS)

    samples = read_sheet(a.sheet)
    gtex = sorted(n for n, s in samples.items() if s.is_gtex)
    patients = sorted({s.patient for s in samples.values() if not s.is_gtex})
    if a.patients:
        keep = set(a.patients.split(","))
        patients = [p for p in patients if p in keep]
        samples = {n: s for n, s in samples.items() if s.is_gtex or s.patient in keep}
    sex_of = {s.patient: s.sex for s in samples.values()}
    normals = [f"{pt}_N" for pt in patients]

    if "align" in steps:
        for s in samples.values():
            align(r, s)
    if "markdup" in steps:
        for s in samples.values():
            markdup_bqsr(r, s)
    if "qc" in steps:
        for s in samples.values():
            bam_qc(r, s)
        identity_check(r, list(samples.values()))
    if "germline" in steps:
        for sm in normals + gtex:
            germline_gvcf(r, sm)
            melt_germline(r, sm)
        germline_joint(r, normals + gtex)
    if "pon" in steps:
        for g in gtex:
            pon_normal(r, g)
        pon_build(r, gtex)
    if "mito" in steps:
        for s in samples.values():
            mito(r, s.name)
        mito_gtex_recurrent(r, gtex)
        for pt in patients:
            mito_compare(r, pt)

    for pt in patients:
        if "somatic" in steps:
            mutect2_multi(r, pt)
        for s in "TB":
            if "sv" in steps:
                manta_strelka(r, pt, s)
                delly_somatic(r, pt, s)
                sv_consensus(r, pt, s)
            if "cnv" in steps:
                ascat_run(r, pt, s, sex_of[pt])
            if "msi" in steps:
                msi(r, pt, s)
        if "annotate" in steps:
            annotate(r, pt)
        if "clonality" in steps:
            clonality(r, pt)
    if "signatures" in steps:
        signatures(r, patients)


if __name__ == "__main__":
    main()
