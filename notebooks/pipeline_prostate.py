#!/usr/bin/env python3
"""Prostate cancer field-effect WGS pipeline (CRUK-ICGC prostate, EGA).

Data: EGAD00001004125 (30 prostatectomy patients: tumour regions, morphologically
normal prostate, blood; plus 7 cancer-free men: prostate + blood) and
EGAD00001000689 (3 patients, 3-5 tumour regions each). Controlled access via
DAC EGAC00001000010. Files are BAMs; they are reverted and realigned to GRCh38.

Groups (samplesheet_prostate.csv):
  T  tumour region (one or more per patient: <PT>_T1, <PT>_T2, ...)
  B  morphologically normal prostate from a cancer patient ("border")
  H  prostate from a cancer-free man (healthy-tissue baseline)
  N  blood (germline control for everyone)

Clinical filter (clinical_prostate.csv): cancer patients are kept only if
ISUP grade group is 2 or 3 and M == M0 (use --pN0 to also require pN0).
Fill the sheet from the EGA/ICGC clinical metadata; NA excludes the patient
unless --allow-missing-clinical is given (useful for dry runs).

Default is a DRY RUN: commands are printed, nothing is executed, nothing is
downloaded. Add --run to execute.

    python pipeline_prostate.py --allow-missing-clinical
    python pipeline_prostate.py --steps revert,align,markdup --run
    python pipeline_prostate.py --steps fetch           # print EGA download commands
"""
from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
THREADS = 16
MEM = "32g"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
EGA_DATASETS = ["EGAD00001004125", "EGAD00001000689"]
KEEP_GRADE_GROUPS = {"2", "3"}

R = Path("ref")
REF = R / "Homo_sapiens_assembly38.fasta"            # + .fai, .dict, bwa-mem2 index, .alt
DBSNP = R / "Homo_sapiens_assembly38.dbsnp138.vcf"
MILLS = R / "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KINDELS = R / "Homo_sapiens_assembly38.known_indels.vcf.gz"
GNOMAD = R / "af-only-gnomad.hg38.vcf.gz"
COMMON = R / "small_exac_common_3.hg38.vcf.gz"
INTERVALS = R / "wgs_calling_regions.hg38.interval_list"
VBID_SVD = R / "verifybamid/1000g.phase3.100k.b38.vcf.gz.dat"
SOMALIER_SITES = R / "somalier/sites.hg38.vcf.gz"
DELLY_EXCL = R / "delly/human.hg38.excl.tsv"
MSI_LIST = R / "hg38.msisensor.list"
ASCAT_REF = R / "ascat"
VEP_CACHE = R / "vep_cache"

MT = R / "chrM"
MT_REF = MT / "Homo_sapiens_assembly38.chrM.fasta"
MT_SHIFT_REF = MT / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
MT_CHAIN = MT / "ShiftBack.chain"
MT_BLACKLIST = MT / "blacklist_sites.hg38.chrM.bed"
MT_NONCONTROL = "chrM:576-16024"
MT_CONTROL_SHIFTED = "chrM:8025-9144"
MT_MIN_AF = 0.01
MT_HOMOPLASMIC = 0.95

# prostate cancer genes: SNV/indel drivers, CN targets and HR/MMR genes
PROSTATE_GENES = [
    "SPOP", "FOXA1", "TP53", "PTEN", "RB1", "CHD1", "NKX3-1", "MYC", "AR",
    "BRCA1", "BRCA2", "ATM", "CDK12", "PALB2", "CHEK2", "HOXB13",
    "KMT2C", "KMT2D", "IDH1", "CTNNB1", "PIK3CA", "APC", "SPEN", "ZMYM3", "MED12",
    "MLH1", "MSH2", "MSH6", "PMS2",
    "ERG", "ETV1", "ETV4", "TMPRSS2",
]
# hg38: ERG ~chr21:38.38 Mb, TMPRSS2 ~chr21:41.46-41.53 Mb (interstitial deletion or BND)
TMPRSS2_ERG_REGION = "chr21:38300000-41600000"

OUT_DIRS = ["ubam", "qc/somalier", "bam", "vcf/pon", "sv", "cnv", "clonality",
            "sig/vcfs", "tmp", "mito", "drivers"]


# ---------------------------------------------------------------------------
# sample + clinical sheets
# ---------------------------------------------------------------------------
@dataclass
class Sample:
    patient: str
    cohort: str      # cancer | healthy
    group: str       # T | B | H | N
    name: str
    bam: str


def read_sheet(path: Path) -> dict[str, Sample]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    samples = {r["sample"]: Sample(r["patient"], r["cohort"], r["group"], r["sample"], r["bam"])
               for r in rows}
    bad = [s.name for s in samples.values() if s.group not in "TBHN"]
    if bad:
        sys.exit(f"unknown group for {bad}")
    return samples


def select_patients(samples: dict[str, Sample], clinical: Path,
                    allow_missing: bool, require_pn0: bool) -> set[str]:
    with open(clinical, newline="") as fh:
        clin = {r["patient"]: r for r in csv.DictReader(fh)}
    keep = set()
    for pt in sorted({s.patient for s in samples.values()}):
        cohort = next(s.cohort for s in samples.values() if s.patient == pt)
        if cohort == "healthy":
            keep.add(pt)
            continue
        c = clin.get(pt, {})
        gg, m, pn = c.get("isup_grade_group", "NA"), c.get("M", "NA"), c.get("pN", "NA")
        missing = "NA" in (gg, m) or (require_pn0 and pn == "NA")
        ok = (gg in KEEP_GRADE_GROUPS and m == "M0" and (not require_pn0 or pn == "pN0"))
        if ok or (missing and allow_missing):
            keep.add(pt)
            if missing:
                print(f"# WARNING {pt}: clinical data missing, kept (--allow-missing-clinical)")
        else:
            print(f"# excluded {pt}: grade group={gg}, M={m}, pN={pn}")
    return keep


# ---------------------------------------------------------------------------
# command runner (dry run by default)
# ---------------------------------------------------------------------------
class Runner:
    def __init__(self, execute: bool):
        self.execute = execute

    @staticmethod
    def _fmt(cmd) -> str:
        return " ".join(shlex.quote(str(c)) for c in cmd)

    def run(self, *cmds, stdin: str | None = None):
        """Run one command, or several joined by pipes."""
        cmds = [[str(c) for c in cmd] for cmd in cmds]
        text = " | ".join(self._fmt(c) for c in cmds)
        if stdin is not None:
            text += f" <<'EOF'\n{stdin.strip()}\nEOF"
        print(text + "\n", flush=True)
        if not self.execute:
            return
        procs, prev = [], None
        for cmd in cmds:
            p = subprocess.Popen(cmd, text=True,
                                 stdin=prev.stdout if prev else (subprocess.PIPE if stdin else None),
                                 stdout=None if cmd is cmds[-1] else subprocess.PIPE)
            if prev:
                prev.stdout.close()
            procs.append(p)
            prev = p
        if stdin is not None:
            procs[0].stdin.write(stdin)
            procs[0].stdin.close()
        codes = [p.wait() for p in procs]
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
# fetch (explicit only): EGA download with pyega3
# ---------------------------------------------------------------------------
def fetch(r: Runner):
    for d in EGA_DATASETS:
        r.run(["pyega3", "-cf", "ega_credentials.json", "files", d])
        r.run(["pyega3", "-cf", "ega_credentials.json", "-c", 8, "fetch", d, "--output-dir", "ega"])
    print("# then map EGA sample aliases -> samplesheet_prostate.csv (T/B/H/N)\n")


# ---------------------------------------------------------------------------
# 1. revert EGA BAMs to unmapped BAMs, one per read group
#    (original BAMs are on an older reference; everything is realigned to GRCh38)
# ---------------------------------------------------------------------------
def revert(r: Runner, s: Sample):
    d = Path("ubam") / s.name
    r.mkdir(d)
    r.run(gatk("RevertSam", "-I", s.bam, "-O", d, "--OUTPUT_BY_READGROUP", "true",
               "--SANITIZE", "true", "--SORT_ORDER", "queryname",
               "--RESTORE_ORIGINAL_QUALITIES", "true", "--REMOVE_DUPLICATE_INFORMATION", "true",
               "--REMOVE_ALIGNMENT_INFORMATION", "true", "--ATTRIBUTE_TO_CLEAR", "XS",
               "--ATTRIBUTE_TO_CLEAR", "XA", "--TMP_DIR", "tmp"))


def ubams(r: Runner, sm: str) -> list[Path]:
    if not r.execute:
        return [Path("ubam") / sm / "<RG>.bam"]
    found = sorted((Path("ubam") / sm).glob("*.bam"))
    if not found:
        sys.exit(f"no unmapped BAMs in ubam/{sm}; run the revert step first")
    return found


def rg_line(r: Runner, ubam: Path, sm: str) -> tuple[str, str]:
    """Read group from the uBAM header, with SM set to the pipeline sample name."""
    if not r.execute:
        return "<RG>", f"@RG\\tID:<RG>\\tSM:{sm}\\tLB:<LB>\\tPL:ILLUMINA\\tPU:<PU>"
    hdr = subprocess.run(["samtools", "view", "-H", str(ubam)],
                         capture_output=True, text=True, check=True).stdout
    rg = next(l for l in hdr.splitlines() if l.startswith("@RG"))
    fields = [f for f in rg.split("\t")[1:] if not f.startswith("SM:")] + [f"SM:{sm}"]
    rgid = next(f[3:] for f in fields if f.startswith("ID:"))
    return rgid, "@RG\\t" + "\\t".join(fields)


# ---------------------------------------------------------------------------
# 2. align each read group; 3. merge + MarkDuplicates + BQSR
# ---------------------------------------------------------------------------
def align(r: Runner, s: Sample):
    for ub in ubams(r, s.name):
        rgid, rg = rg_line(r, ub, s.name)
        r.run(["samtools", "fastq", "-n", "-0", "/dev/null", "-s", "/dev/null", ub],
              ["bwa-mem2", "mem", "-t", THREADS, "-p", "-Y", "-K", 100000000, "-R", rg, REF, "-"],
              ["samtools", "sort", "-@", 4, "-m", "2G", "-o", f"bam/{s.name}.{rgid}.bam", "-"])


def markdup_bqsr(r: Runner, s: Sample):
    sm = s.name
    if r.execute:
        parts = sorted(p for p in Path("bam").glob(f"{sm}.*.bam")
                       if not p.name.endswith((".md.bam", ".bqsr.bam")))
    else:
        parts = [Path("bam") / f"{sm}.<RG>.bam"]
    ins = [x for p in parts for x in ("-I", p)]
    r.run(gatk("MarkDuplicates", *ins, "-O", f"bam/{sm}.md.bam", "-M", f"qc/{sm}.md_metrics.txt",
               "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", 100, "--CREATE_INDEX", "true"))
    #   pixel distance 100: HiSeq 2000 (unpatterned); use 2500 for patterned flowcells
    r.run(gatk("BaseRecalibrator", "-R", REF, "-I", f"bam/{sm}.md.bam",
               "--known-sites", DBSNP, "--known-sites", MILLS, "--known-sites", KINDELS,
               "-O", f"bam/{sm}.recal.table"))
    r.run(gatk("ApplyBQSR", "-R", REF, "-I", f"bam/{sm}.md.bam",
               "--bqsr-recal-file", f"bam/{sm}.recal.table", "-O", bqsr_bam(sm)))


# ---------------------------------------------------------------------------
# 4. QC
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
    """All samples of one man must be relatedness ~1 to each other, ~0 to others."""
    r.run(["somalier", "relate", "--infer", "-o", "qc/somalier_relate",
           *[f"qc/somalier/{s.name}.somalier" for s in samples]])


# ---------------------------------------------------------------------------
# 5. germline (bloods): predisposition (BRCA2, ATM, CHEK2, HOXB13 G84E, ...)
# ---------------------------------------------------------------------------
def germline(r: Runner, bloods: list[str]):
    for sm in bloods:
        r.run(gatk("HaplotypeCaller", "-R", REF, "-I", bqsr_bam(sm), "-L", INTERVALS,
                   "-ERC", "GVCF", "-O", f"vcf/{sm}.g.vcf.gz"))
    vs = [x for sm in bloods for x in ("-V", f"vcf/{sm}.g.vcf.gz")]
    for c in CHROMS:
        r.run(gatk("GenomicsDBImport", *vs, "-L", c, "--genomicsdb-workspace-path", f"vcf/gdb_{c}"))
        r.run(gatk("GenotypeGVCFs", "-R", REF, "-V", f"gendb://vcf/gdb_{c}",
                   "-O", f"vcf/germline.{c}.vcf.gz"))
    r.run(gatk("GatherVcfs", *[x for c in CHROMS for x in ("-I", f"vcf/germline.{c}.vcf.gz")],
               "-O", "vcf/germline.vcf.gz"))
    r.run(["bcftools", "index", "-t", "vcf/germline.vcf.gz"])
    # then VQSR (37 WGS is enough) or hard filters, and VEP on the result


# ---------------------------------------------------------------------------
# 6. panel of normals: BLOOD ONLY (same platform/batch, all 37 men).
#    Prostate tissue (B, H) is deliberately excluded: clonal expansions shared
#    across prostates are the field-effect signal and must not be filtered out.
# ---------------------------------------------------------------------------
def pon(r: Runner, bloods: list[str]):
    for sm in bloods:
        r.run(gatk("Mutect2", "-R", REF, "-I", bqsr_bam(sm), "-L", INTERVALS,
                   "--max-mnp-distance", 0, "-O", f"vcf/pon/{sm}.vcf.gz"))
    vs = [x for sm in bloods for x in ("-V", f"vcf/pon/{sm}.vcf.gz")]
    r.run(gatk("GenomicsDBImport", "-R", REF, "-L", INTERVALS,
               "--genomicsdb-workspace-path", "vcf/pon_db", *vs))
    r.run(gatk("CreateSomaticPanelOfNormals", "-R", REF, "-V", "gendb://vcf/pon_db",
               "--germline-resource", GNOMAD, "-O", "vcf/pon.vcf.gz"))


# ---------------------------------------------------------------------------
# 7. somatic SNV/indel: one joint Mutect2 per man
#    cancer: all tumour regions + B vs blood; healthy: H vs blood
# ---------------------------------------------------------------------------
def mutect2_joint(r: Runner, pt: str, cases: list[str], normal: str):
    d = Path("vcf") / pt
    r.mkdir(d)
    ins = [x for sm in cases + [normal] for x in ("-I", bqsr_bam(sm))]
    for c in CHROMS:
        r.run(gatk("Mutect2", "-R", REF, *ins, "-normal", normal,
                   "--germline-resource", GNOMAD, "--panel-of-normals", "vcf/pon.vcf.gz",
                   "-L", c, "--f1r2-tar-gz", d / f"f1r2.{c}.tar.gz", "-O", d / f"unf.{c}.vcf.gz"))
    r.run(gatk("MergeVcfs", *[x for c in CHROMS for x in ("-I", d / f"unf.{c}.vcf.gz")],
               "-O", d / "unfiltered.vcf.gz"))
    r.run(gatk("MergeMutectStats",
               *[x for c in CHROMS for x in ("-stats", d / f"unf.{c}.vcf.gz.stats")],
               "-O", d / "unfiltered.vcf.gz.stats"))
    r.run(gatk("LearnReadOrientationModel",
               *[x for c in CHROMS for x in ("-I", d / f"f1r2.{c}.tar.gz")], "-O", d / "rom.tar.gz"))

    for sm in cases + [normal]:
        r.run(gatk("GetPileupSummaries", "-I", bqsr_bam(sm), "-V", COMMON, "-L", COMMON,
                   "-O", d / f"{sm}.pileups.table"))
    contam = []
    for sm in cases:
        r.run(gatk("CalculateContamination", "-I", d / f"{sm}.pileups.table",
                   "-matched", d / f"{normal}.pileups.table",
                   "--tumor-segmentation", d / f"{sm}.segments.table",
                   "-O", d / f"{sm}.contamination.table"))
        contam += ["--contamination-table", d / f"{sm}.contamination.table",
                   "--tumor-segmentation", d / f"{sm}.segments.table"]
    r.run(gatk("FilterMutectCalls", "-R", REF, "-V", d / "unfiltered.vcf.gz", *contam,
               "--ob-priors", d / "rom.tar.gz", "-O", d / "mutect2.filtered.vcf.gz"))
    r.run(["bcftools", "view", "-f", "PASS", d / "mutect2.filtered.vcf.gz",
           "-Oz", "-o", d / "mutect2.pass.vcf.gz"])
    r.run(["bcftools", "index", "-t", d / "mutect2.pass.vcf.gz"])


# ---------------------------------------------------------------------------
# 8. SVs (Manta + Strelka2, Delly, consensus); 9. CNV/purity (ASCAT); 10. MSI
# ---------------------------------------------------------------------------
def sv(r: Runner, case: str, normal: str):
    d = Path("sv") / case
    T, N = bqsr_bam(case), bqsr_bam(normal)
    r.mkdir(d)
    r.run(["configManta.py", "--normalBam", N, "--tumorBam", T, "--referenceFasta", REF,
           "--runDir", d / "manta"])
    r.run([d / "manta/runWorkflow.py", "-m", "local", "-j", THREADS])
    r.run(["configureStrelkaSomaticWorkflow.py", "--normalBam", N, "--tumorBam", T,
           "--referenceFasta", REF,
           "--indelCandidates", d / "manta/results/variants/candidateSmallIndels.vcf.gz",
           "--runDir", d / "strelka"])
    r.run([d / "strelka/runWorkflow.py", "-m", "local", "-j", THREADS])
    r.run(["delly", "call", "-x", DELLY_EXCL, "-g", REF, "-o", d / "delly.bcf", T, N])
    r.write(d / "delly_samples.tsv", f"{case}\ttumor\n{normal}\tcontrol\n")
    r.run(["delly", "filter", "-f", "somatic", "-s", d / "delly_samples.tsv",
           "-o", d / "delly.somatic.bcf", d / "delly.bcf"])
    r.run(["bcftools", "view", "-f", "PASS", d / "manta/results/variants/somaticSV.vcf.gz",
           "-Ov", "-o", d / "manta.vcf"])
    r.run(["bcftools", "view", "-f", "PASS", d / "delly.somatic.bcf", "-Ov", "-o", d / "delly.vcf"])
    r.write(d / "sv_list.txt", f"{d / 'manta.vcf'}\n{d / 'delly.vcf'}\n")
    r.run(["SURVIVOR", "merge", d / "sv_list.txt", 1000, 2, 1, 1, 0, 50, d / "sv.consensus.vcf"])


ASCAT_R = r"""
args <- commandArgs(trailingOnly = TRUE)
tum <- args[1]; nor <- args[2]; r <- args[3]
library(ASCAT)
out <- file.path("cnv", tum); dir.create(out, recursive = TRUE, showWarnings = FALSE)
f <- function(x) file.path(out, x)
ascat.prepareHTS(
  tumourseqfile = file.path("bam", paste0(tum, ".bqsr.bam")),
  normalseqfile = file.path("bam", paste0(nor, ".bqsr.bam")),
  tumourname = tum, normalname = nor, allelecounter_exe = "alleleCounter",
  alleles.prefix = file.path(r, "G1000_alleles_hg38_chr"),
  loci.prefix    = file.path(r, "G1000_loci_hg38_chr"),
  gender = "XY", genomeVersion = "hg38", nthreads = 8,
  tumourLogR_file = f("T_LogR.txt"), tumourBAF_file = f("T_BAF.txt"),
  normalLogR_file = f("N_LogR.txt"), normalBAF_file = f("N_BAF.txt"))
bc <- ascat.loadData(Tumor_LogR_file = f("T_LogR.txt"), Tumor_BAF_file = f("T_BAF.txt"),
                     Germline_LogR_file = f("N_LogR.txt"), Germline_BAF_file = f("N_BAF.txt"),
                     gender = "XY", genomeVersion = "hg38")
bc <- ascat.correctLogR(bc, GCcontentfile = file.path(r, "GC_G1000_hg38.txt"),
                        replictimingfile = file.path(r, "RT_G1000_hg38.txt"))
bc <- ascat.aspcf(bc)
res <- ascat.runAscat(bc, gamma = 1, write_segments = TRUE)
write.table(res$segments, f("segments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = tum, purity = res$aberrantcellfraction, ploidy = res$ploidy),
            f("purity_ploidy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
"""


def ascat(r: Runner, case: str, normal: str):
    r.run(["Rscript", "-", case, normal, ASCAT_REF], stdin=ASCAT_R)


def msi(r: Runner, case: str, normal: str):
    r.run(["msisensor-pro", "msi", "-d", MSI_LIST, "-n", bqsr_bam(normal), "-t", bqsr_bam(case),
           "-o", f"qc/{case}.msi", "-b", THREADS])


# ---------------------------------------------------------------------------
# 11. annotation + prostate driver report
# ---------------------------------------------------------------------------
def annotate(r: Runner, pt: str):
    d = Path("vcf") / pt
    r.run(["vep", "--offline", "--cache", "--dir_cache", VEP_CACHE, "--assembly", "GRCh38",
           "--fasta", REF, "--everything", "--pick", "--vcf", "--compress_output", "bgzip",
           "--fork", 8, "-i", d / "mutect2.pass.vcf.gz", "-o", d / "mutect2.pass.vep.vcf.gz"])


def drivers(r: Runner, pt: str, cases: list[str]):
    d = Path("vcf") / pt
    genes = "|".join(PROSTATE_GENES)
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%HGVSp\\t%IMPACT[\\t%SAMPLE=%AD]\\n"
    r.run(["bash", "-c",
           f"bcftools +split-vep {d}/mutect2.pass.vep.vcf.gz -d -f '{fmt}' "
           f"-i 'SYMBOL~\"^({genes})$\" && IMPACT~\"HIGH|MODERATE\"' "
           f"> drivers/{pt}.snv_indel.tsv"])
    for case in cases:   # TMPRSS2-ERG: deletion or breakends in the chr21 window
        r.run(["bash", "-c",
               f"bcftools view -t {TMPRSS2_ERG_REGION} sv/{case}/sv.consensus.vcf "
               f"| bcftools query -f '%CHROM\\t%POS\\t%INFO/SVTYPE\\t%INFO/END\\t%ALT\\n' "
               f"> drivers/{case}.tmprss2_erg.tsv"])


# ---------------------------------------------------------------------------
# 12. mitochondria (same workflow as pipeline.py, generalised to N samples)
# ---------------------------------------------------------------------------
def mt_dir(sm: str) -> Path:
    return Path("mito") / sm


def mosdepth_means(sm: str) -> tuple[float, float]:
    auto_bases = auto_len = 0.0
    mt = None
    autosomes = {f"chr{i}" for i in range(1, 23)}
    with open(f"qc/{sm}.mosdepth.summary.txt") as fh:
        next(fh)
        for line in fh:
            chrom, length, bases, mean, *_ = line.split("\t")
            if chrom in autosomes:
                auto_len += float(length)
                auto_bases += float(bases)
            elif chrom == "chrM":
                mt = float(mean)
    if mt is None or auto_len == 0:
        sys.exit(f"chrM/autosomes missing in qc/{sm}.mosdepth.summary.txt")
    return auto_bases / auto_len, mt


def parse_haplocheck(path: Path) -> float:
    with open(path, newline="") as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    key = next((k for k in row if "contamination level" in k.lower().strip('"')), None)
    try:
        return float(row[key].strip('"')) if key else 0.0
    except ValueError:
        return 0.0


def mito_filter(r: Runner, sm: str, tag: str, contamination, autosomal_cov):
    d = mt_dir(sm)
    r.run(gatk("FilterMutectCalls", "-R", MT_REF, "-V", d / "raw.vcf.gz",
               "--stats", d / "raw.vcf.gz.stats", "--mitochondria-mode",
               "--max-alt-allele-count", 4, "--min-allele-fraction", 0,
               "--autosomal-coverage", autosomal_cov, "--contamination-estimate", contamination,
               "-O", d / f"{tag}.filtered.vcf.gz"))
    r.run(gatk("VariantFiltration", "-R", MT_REF, "-V", d / f"{tag}.filtered.vcf.gz",
               "--apply-allele-specific-filters", "--mask", MT_BLACKLIST,
               "--mask-name", "blacklisted_site", "-O", d / f"{tag}.masked.vcf.gz"))
    r.run(gatk("LeftAlignAndTrimVariants", "-R", MT_REF, "-V", d / f"{tag}.masked.vcf.gz",
               "--split-multi-allelics", "--dont-trim-alleles", "--keep-original-ac",
               "-O", d / f"{tag}.split.vcf.gz"))


def mito(r: Runner, sm: str):
    d = mt_dir(sm)
    r.mkdir(d)
    r.run(gatk("PrintReads", "-R", REF, "-I", f"bam/{sm}.md.bam", "-L", "chrM",
               "--read-filter", "MateOnSameContigOrNoMappedMateReadFilter",
               "--read-filter", "MateUnmappedAndUnmappedReadFilter", "-O", d / "chrM.subset.bam"))
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
    return {(v.POS, v.REF, v.ALT[0]): dict(AF=float(v.format("AF")[0][0]),
                                          DP=int(v.format("DP")[0][0]),
                                          FILTER=v.FILTER or "PASS")
            for v in VCF(str(mt_dir(sm) / "final.split.vcf.gz"))}


def mito_recurrent(r: Runner, bloods: list[str]):
    """Heteroplasmic PASS sites in >=2 unrelated men's blood: artefact/hotspot flag."""
    print(f"# python: recurrent heteroplasmy across {len(bloods)} bloods -> mito/recurrent.tsv\n")
    if not r.execute:
        return
    import pandas as pd
    rows = [dict(POS=p, REF=a, ALT=b, sample=sm)
            for sm in bloods for (p, a, b), c in read_mt_calls(sm).items()
            if c["FILTER"] == "PASS" and MT_MIN_AF <= c["AF"] < MT_HOMOPLASMIC]
    df = pd.DataFrame(rows, columns=["POS", "REF", "ALT", "sample"])
    rec = df.groupby(["POS", "REF", "ALT"])["sample"].nunique().rename("n_men").reset_index()
    rec[rec.n_men >= 2].to_csv("mito/recurrent.tsv", sep="\t", index=False)


def mito_compare(r: Runner, pt: str, cases: list[str], normal: str):
    out = Path("mito") / f"{pt}.heteroplasmy.tsv"
    print(f"# python: heteroplasmy {cases} vs {normal}, classes + mtCN -> {out}\n")
    if not r.execute:
        return
    import pandas as pd
    allsm = cases + [normal]
    calls = {sm: read_mt_calls(sm) for sm in allsm}
    rec = set()
    if Path("mito/recurrent.tsv").exists():
        g = pd.read_csv("mito/recurrent.tsv", sep="\t")
        rec = set(zip(g.POS, g.REF, g.ALT))
    rows = []
    for k in sorted(set().union(*calls.values())):
        row = dict(POS=k[0], REF=k[1], ALT=k[2])
        for sm in allsm:
            c = calls[sm].get(k, dict(AF=0.0, DP=None, FILTER="absent"))
            row.update({f"AF_{sm}": c["AF"], f"DP_{sm}": c["DP"], f"FILTER_{sm}": c["FILTER"]})
        n = row[f"AF_{normal}"]
        case_af = [row[f"AF_{sm}"] for sm in cases]
        case_pass = [row[f"FILTER_{sm}"] == "PASS" for sm in cases]
        if min(case_af + [n]) >= MT_HOMOPLASMIC:
            cls = "germline_homoplasmic"
        elif n < MT_MIN_AF and any(p and af >= MT_MIN_AF for af, p in zip(case_af, case_pass)):
            present = [sm for sm, af, p in zip(cases, case_af, case_pass) if p and af >= MT_MIN_AF]
            cls = "somatic:" + ",".join(present)   # which tissues share it (lineage marker)
        elif n >= MT_MIN_AF and any(abs(af - n) >= 0.10 for af in case_af):
            cls = "heteroplasmy_shift"
        else:
            cls = "other"
        row.update({"class": cls, "recurrent_in_bloods": k in rec})
        rows.append(row)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    pd.concat(pd.read_csv(mt_dir(sm) / "mtcn.tsv", sep="\t") for sm in allsm) \
      .to_csv(Path("mito") / f"{pt}.mtcn.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# 13. clonality (all tumour regions + border); 14. signatures (incl. healthy)
# ---------------------------------------------------------------------------
def clonality(r: Runner, pt: str, cases: list[str]):
    # If ASCAT fails on the low-purity border sample, point it at a tumour region's
    # segments and give its purity explicitly, e.g.
    #   --sample PR01_B=cnv/PR01_T1 --purity PR01_B=0.08
    samples = [x for sm in cases for x in ("--sample", f"{sm}=cnv/{sm}")]
    r.run([sys.executable, "make_pyclone_input.py", "--vcf", f"vcf/{pt}/mutect2.pass.vcf.gz",
           *samples, "-o", f"clonality/{pt}.pyclone_in.tsv"])
    r.run(["pyclone-vi", "fit", "-i", f"clonality/{pt}.pyclone_in.tsv",
           "-o", f"clonality/{pt}.h5", "-c", 40, "-d", "beta-binomial", "-r", 10])
    r.run(["pyclone-vi", "write-results-file", "-i", f"clonality/{pt}.h5",
           "-o", f"clonality/{pt}.pyclone_out.tsv"])


def signatures(r: Runner, per_man: dict[str, list[str]]):
    for pt, cases in per_man.items():
        for sm in cases:
            r.run(["bcftools", "view", "-s", sm, f"vcf/{pt}/mutect2.pass.vcf.gz"],
                  ["bcftools", "view", "-i", "FMT/AD[0:1]>=3", "-Ov", "-o", f"sig/vcfs/{sm}.vcf"])
    print("# SigProfilerAssignment.cosmic_fit(samples='sig/vcfs', output='sig/out', "
          "input_type='vcf', genome_build='GRCh38')\n")
    if r.execute:
        from SigProfilerAssignment import Analyzer as Analyze
        Analyze.cosmic_fit(samples="sig/vcfs", output="sig/out",
                           input_type="vcf", genome_build="GRCh38")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
STEPS = ["revert", "align", "markdup", "qc", "germline", "pon", "somatic", "mito",
         "sv", "cnv", "msi", "annotate", "drivers", "clonality", "signatures"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", type=Path, default=Path("samplesheet_prostate.csv"))
    ap.add_argument("--clinical", type=Path, default=Path("clinical_prostate.csv"))
    ap.add_argument("--run", action="store_true", help="execute (default: dry run)")
    ap.add_argument("--steps", default=",".join(STEPS),
                    help=f"comma-separated subset of: fetch,{','.join(STEPS)}")
    ap.add_argument("--allow-missing-clinical", action="store_true")
    ap.add_argument("--pN0", action="store_true", help="also require pN0")
    a = ap.parse_args()

    steps = set(a.steps.split(","))
    if unknown := steps - set(STEPS) - {"fetch"}:
        ap.error(f"unknown steps: {sorted(unknown)}")
    r = Runner(a.run)
    if not a.run:
        print("# DRY RUN: commands are printed, nothing is executed. Add --run to execute.\n")
    if "fetch" in steps:
        fetch(r)
        if steps == {"fetch"}:
            return
    r.mkdir(*OUT_DIRS)

    samples = read_sheet(a.sheet)
    keep = select_patients(samples, a.clinical, a.allow_missing_clinical, a.pN0)
    samples = {n: s for n, s in samples.items() if s.patient in keep}
    by_man: dict[str, list[Sample]] = {}
    for s in samples.values():
        by_man.setdefault(s.patient, []).append(s)

    blood = {}
    cases = {}          # man -> non-blood samples (T regions + B, or H)
    tumours = {}        # man -> tumour regions only
    for pt, ss in sorted(by_man.items()):
        ns = [s.name for s in ss if s.group == "N"]
        if len(ns) != 1:
            sys.exit(f"{pt}: expected exactly one blood (N) sample, found {ns}")
        blood[pt] = ns[0]
        cases[pt] = sorted(s.name for s in ss if s.group in "TBH")
        tumours[pt] = sorted(s.name for s in ss if s.group == "T")
    bloods = sorted(blood.values())
    print(f"# men kept: {len(by_man)} "
          f"(cancer {sum(1 for p in tumours if tumours[p])}, healthy "
          f"{sum(1 for p in tumours if not tumours[p])}); samples: {len(samples)}\n")

    for s in samples.values():
        if "revert" in steps:
            revert(r, s)
        if "align" in steps:
            align(r, s)
        if "markdup" in steps:
            markdup_bqsr(r, s)
        if "qc" in steps:
            bam_qc(r, s)
    if "qc" in steps:
        identity_check(r, list(samples.values()))
    if "germline" in steps:
        germline(r, bloods)
    if "pon" in steps:
        pon(r, bloods)
    if "mito" in steps:
        for s in samples.values():
            mito(r, s.name)
        mito_recurrent(r, bloods)
        for pt in by_man:
            mito_compare(r, pt, cases[pt], blood[pt])

    for pt in sorted(by_man):
        N = blood[pt]
        if "somatic" in steps:
            mutect2_joint(r, pt, cases[pt], N)
        for c in cases[pt]:
            if "sv" in steps:
                sv(r, c, N)
            if "cnv" in steps and samples[c].group != "H":
                ascat(r, c, N)            # tumour regions and border; not healthy prostate
        for t in tumours[pt]:
            if "msi" in steps:
                msi(r, t, N)
        if "annotate" in steps:
            annotate(r, pt)
        if "drivers" in steps:
            drivers(r, pt, cases[pt])
        if "clonality" in steps and tumours[pt]:
            clonality(r, pt, cases[pt])
    if "signatures" in steps:
        signatures(r, cases)


if __name__ == "__main__":
    main()
