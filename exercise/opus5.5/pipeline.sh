#!/usr/bin/env bash
# =============================================================================
# Multi-region somatic WGS
#   per patient: tumour (<PT>_T), border (<PT>_B), non-tumour tissue (<PT>_N)
#   plus GTEx donors (GTEX-*) as population / panel-of-normals controls
# Reference: GRCh38, Broad resource bundle file names.
# NOT EXECUTED. Review paths, resources and tool versions first.
# Guard: the script only runs with ./pipeline.sh --run
#
# Production alternative covering steps 1-10 (same sample sheet format):
#   nextflow run nf-core/sarek -r <version> -profile docker \
#     --input samplesheet.csv --genome GATK.GRCh38 --outdir results \
#     --tools haplotypecaller,mutect2,strelka,manta,ascat,msisensorpro,vep \
#     --joint_mutect2
# =============================================================================
set -euo pipefail
[[ "${1:-}" == "--run" ]] || { echo "Not running. Review, then: ./pipeline.sh --run"; exit 0; }

THREADS=16
MEM=32g
SHEET=samplesheet.csv
GATK="gatk --java-options -Xmx${MEM}"

# ---- reference & resources --------------------------------------------------
REF=ref/Homo_sapiens_assembly38.fasta              # + .fai, .dict, bwa-mem2 index, .alt
DBSNP=ref/Homo_sapiens_assembly38.dbsnp138.vcf
MILLS=ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
KINDELS=ref/Homo_sapiens_assembly38.known_indels.vcf.gz
GNOMAD=ref/af-only-gnomad.hg38.vcf.gz               # Mutect2 germline resource
COMMON=ref/small_exac_common_3.hg38.vcf.gz          # contamination sites
INTERVALS=ref/wgs_calling_regions.hg38.interval_list
VBID_SVD=ref/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat
SOMALIER_SITES=ref/somalier/sites.hg38.vcf.gz
DELLY_EXCL=ref/delly/human.hg38.excl.tsv
MSI_LIST=ref/hg38.msisensor.list                    # once: msisensor-pro scan -d $REF -o $MSI_LIST
ASCAT_REF=ref/ascat                                 # G1000 alleles/loci, GC, RT files (hg38)
VEP_CACHE=ref/vep_cache

CHROMS=( $(printf 'chr%s ' {1..22} X Y) )
mkdir -p qc/somalier bam vcf/pon sv cnv te clonality sig tmp

# =============================================================================
# 0. GTEx: WGS comes as aligned CRAMs (dbGaP phs000424, via AnVIL/Terra).
#    Revert to FASTQ so every sample goes through the SAME aligner and version.
#    (To keep original read groups, use gatk RevertSam + SamToFastq per RG instead.)
# =============================================================================
gtex_to_fastq() {   # $1 donor id  $2 cram
  local id=$1 cram=$2
  samtools collate -@ "$THREADS" -Ou --reference "$REF" "$cram" tmp/"$id" \
    | samtools fastq -@ "$THREADS" -n \
        -1 fastq/"$id"_L001_R1.fastq.gz -2 fastq/"$id"_L001_R2.fastq.gz \
        -0 /dev/null -s /dev/null -
}

# =============================================================================
# 1-2. FASTQ QC + alignment per lane (no adapter trimming: BWA soft-clips)
# =============================================================================
align_lane() {      # patient sample lane fq1 fq2
  local sm=$2 ln=$3 fq1=$4 fq2=$5
  local rg="@RG\tID:${sm}.${ln}\tSM:${sm}\tLB:${sm}_lib1\tPL:ILLUMINA\tPU:${sm}.${ln}"
  fastqc -t 4 -o qc "$fq1" "$fq2"
  bwa-mem2 mem -t "$THREADS" -Y -K 100000000 -R "$rg" "$REF" "$fq1" "$fq2" \
    | samtools sort -@ 4 -m 2G -o bam/"${sm}.${ln}".bam -
}

# =============================================================================
# 3. Merge lanes + MarkDuplicates + BQSR
# =============================================================================
markdup_bqsr() {    # sample
  local sm=$1 ins=()
  for b in bam/"${sm}".L*.bam; do ins+=( -I "$b" ); done
  $GATK MarkDuplicates "${ins[@]}" -O bam/"$sm".md.bam -M qc/"$sm".md_metrics.txt \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --CREATE_INDEX true
  $GATK BaseRecalibrator -R "$REF" -I bam/"$sm".md.bam \
      --known-sites "$DBSNP" --known-sites "$MILLS" --known-sites "$KINDELS" \
      -O bam/"$sm".recal.table
  $GATK ApplyBQSR -R "$REF" -I bam/"$sm".md.bam \
      --bqsr-recal-file bam/"$sm".recal.table -O bam/"$sm".bqsr.bam
}

# =============================================================================
# 4. BAM QC: depth, WGS metrics, insert size, contamination, sample identity
# =============================================================================
bam_qc() {          # sample
  local sm=$1 b=bam/$1.bqsr.bam
  mosdepth -t 4 -n --fast-mode --by 1000000 qc/"$sm" "$b"
  $GATK CollectWgsMetrics -R "$REF" -I "$b" -O qc/"$sm".wgs_metrics.txt
  $GATK CollectInsertSizeMetrics -I "$b" -O qc/"$sm".insert.txt -H qc/"$sm".insert.pdf
  verifybamid2 --SVDPrefix "$VBID_SVD" --Reference "$REF" --BamFile "$b" --Output qc/"$sm".vbid2
  somalier extract -d qc/somalier/ --sites "$SOMALIER_SITES" -f "$REF" "$b"
}
identity_check() {  # T/B/N of one patient must be relatedness ~1; GTEx unrelated
  somalier relate --infer -o qc/somalier_relate qc/somalier/*.somalier
}

# =============================================================================
# 5. Germline (N + GTEx): predisposition variants, population comparison
# =============================================================================
germline_gvcf() {   # sample
  $GATK HaplotypeCaller -R "$REF" -I bam/"$1".bqsr.bam -L "$INTERVALS" \
      -ERC GVCF -O vcf/"$1".g.vcf.gz
}
germline_joint() {
  local vs=()
  for g in vcf/*_N.g.vcf.gz vcf/GTEX-*.g.vcf.gz; do vs+=( -V "$g" ); done
  for c in "${CHROMS[@]}"; do
    $GATK GenomicsDBImport "${vs[@]}" -L "$c" --genomicsdb-workspace-path vcf/gdb_"$c"
    $GATK GenotypeGVCFs -R "$REF" -V gendb://vcf/gdb_"$c" -O vcf/germline."$c".vcf.gz
  done
  # then: MergeVcfs -> VariantRecalibrator/ApplyVQSR (WGS cohort) or hard filters
}

# =============================================================================
# 6. Panel of Normals
#    Built from GTEx (processed identically). Adjacent non-tumour tissue (N) is
#    left out on purpose: tumour-in-normal contamination could put recurrent
#    hotspot mutations (e.g., KRAS G12) into the PoN and filter them everywhere.
#    Scatter GenomicsDBImport by chromosome for WGS.
# =============================================================================
pon_normal() {      # sample (GTEx donor)
  $GATK Mutect2 -R "$REF" -I bam/"$1".bqsr.bam -L "$INTERVALS" \
      --max-mnp-distance 0 -O vcf/pon/"$1".vcf.gz
}
pon_build() {
  local vs=()
  for v in vcf/pon/GTEX-*.vcf.gz; do vs+=( -V "$v" ); done
  $GATK GenomicsDBImport -R "$REF" -L "$INTERVALS" \
      --genomicsdb-workspace-path vcf/pon_db "${vs[@]}"
  $GATK CreateSomaticPanelOfNormals -R "$REF" -V gendb://vcf/pon_db \
      --germline-resource "$GNOMAD" -O vcf/pon.vcf.gz
}

# =============================================================================
# 7. Somatic SNV/indel: joint multi-sample Mutect2 (T + B vs N)
#    Joint calling gives ref/alt counts in BOTH T and B at every site, which is
#    what the clonality step needs.
# =============================================================================
mutect2_multi() {   # patient
  local pt=$1 d=vcf/$1
  local T=bam/${pt}_T.bqsr.bam B=bam/${pt}_B.bqsr.bam N=bam/${pt}_N.bqsr.bam
  mkdir -p "$d"
  for c in "${CHROMS[@]}"; do
    $GATK Mutect2 -R "$REF" -I "$T" -I "$B" -I "$N" -normal "${pt}_N" \
        --germline-resource "$GNOMAD" --panel-of-normals vcf/pon.vcf.gz -L "$c" \
        --f1r2-tar-gz "$d"/f1r2."$c".tar.gz -O "$d"/unf."$c".vcf.gz
  done
  local vs=() ss=() fs=()
  for c in "${CHROMS[@]}"; do
    vs+=( -I "$d"/unf."$c".vcf.gz )
    ss+=( -stats "$d"/unf."$c".vcf.gz.stats )
    fs+=( -I "$d"/f1r2."$c".tar.gz )
  done
  $GATK MergeVcfs "${vs[@]}" -O "$d"/unfiltered.vcf.gz
  $GATK MergeMutectStats "${ss[@]}" -O "$d"/unfiltered.vcf.gz.stats
  $GATK LearnReadOrientationModel "${fs[@]}" -O "$d"/rom.tar.gz

  for s in T B N; do
    $GATK GetPileupSummaries -I bam/"${pt}_$s".bqsr.bam -V "$COMMON" -L "$COMMON" \
        -O "$d"/"$s".pileups.table
  done
  for s in T B; do
    $GATK CalculateContamination -I "$d"/"$s".pileups.table -matched "$d"/N.pileups.table \
        --tumor-segmentation "$d"/"$s".segments.table -O "$d"/"$s".contamination.table
  done
  $GATK FilterMutectCalls -R "$REF" -V "$d"/unfiltered.vcf.gz \
      --contamination-table "$d"/T.contamination.table \
      --contamination-table "$d"/B.contamination.table \
      --tumor-segmentation "$d"/T.segments.table \
      --tumor-segmentation "$d"/B.segments.table \
      --ob-priors "$d"/rom.tar.gz -O "$d"/mutect2.filtered.vcf.gz
  bcftools view -f PASS "$d"/mutect2.filtered.vcf.gz -Oz -o "$d"/mutect2.pass.vcf.gz
  bcftools index -t "$d"/mutect2.pass.vcf.gz
}

# Mitochondrial heteroplasmy: lineage markers shared between T and B.
# (The full GATK mito workflow also calls on a shifted chrM for the control region.)
mito() {            # sample
  $GATK Mutect2 -R "$REF" -I bam/"$1".bqsr.bam -L chrM --mitochondria-mode \
      -O vcf/"$1".chrM.vcf.gz
  $GATK FilterMutectCalls -R "$REF" --mitochondria-mode \
      -V vcf/"$1".chrM.vcf.gz -O vcf/"$1".chrM.filt.vcf.gz
}

# =============================================================================
# 8. Second SNV/indel caller + SVs: Manta -> Strelka2, Delly (pairwise vs N)
# =============================================================================
manta_strelka() {   # patient T|B
  local pt=$1 s=$2 d=sv/$1_$2
  local T=bam/${pt}_${s}.bqsr.bam N=bam/${pt}_N.bqsr.bam
  configManta.py --normalBam "$N" --tumorBam "$T" --referenceFasta "$REF" --runDir "$d"/manta
  "$d"/manta/runWorkflow.py -m local -j "$THREADS"
  configureStrelkaSomaticWorkflow.py --normalBam "$N" --tumorBam "$T" --referenceFasta "$REF" \
      --indelCandidates "$d"/manta/results/variants/candidateSmallIndels.vcf.gz \
      --runDir "$d"/strelka
  "$d"/strelka/runWorkflow.py -m local -j "$THREADS"
}
delly_somatic() {   # patient T|B
  local pt=$1 s=$2 d=sv/$1_$2
  mkdir -p "$d"
  delly call -x "$DELLY_EXCL" -g "$REF" -o "$d"/delly.bcf \
      bam/"${pt}_${s}".bqsr.bam bam/"${pt}_N".bqsr.bam
  printf '%s\ttumor\n%s\tcontrol\n' "${pt}_${s}" "${pt}_N" > "$d"/delly_samples.tsv
  delly filter -f somatic -s "$d"/delly_samples.tsv -o "$d"/delly.somatic.bcf "$d"/delly.bcf
}
sv_consensus() {    # patient T|B : SVs supported by >=2 callers, within 1 kb
  local d=sv/$1_$2
  bcftools view -f PASS "$d"/manta/results/variants/somaticSV.vcf.gz -Ov -o "$d"/manta.vcf
  bcftools view -f PASS "$d"/delly.somatic.bcf -Ov -o "$d"/delly.vcf
  printf '%s\n' "$d"/manta.vcf "$d"/delly.vcf > "$d"/sv_list.txt
  SURVIVOR merge "$d"/sv_list.txt 1000 2 1 1 0 50 "$d"/sv.consensus.vcf
}

# =============================================================================
# 9. Copy number, purity, ploidy (ASCAT). Border samples with very low tumour
#    fraction may not fit; estimate purity from clonal SNV VAFs instead.
# =============================================================================
ascat_run() {       # patient T|B sex(XX|XY)
Rscript - "$1" "$2" "$3" "$ASCAT_REF" <<'EOF'
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
EOF
}

# =============================================================================
# 10. MSI, annotation
# =============================================================================
msi() {             # patient T|B
  msisensor-pro msi -d "$MSI_LIST" -n bam/"$1_N".bqsr.bam -t bam/"$1_$2".bqsr.bam \
      -o qc/"$1_$2".msi -b "$THREADS"
}
annotate() {        # patient
  local d=vcf/$1
  vep --offline --cache --dir_cache "$VEP_CACHE" --assembly GRCh38 --fasta "$REF" \
      --everything --vcf --compress_output bgzip --fork 8 \
      -i "$d"/mutect2.pass.vcf.gz -o "$d"/mutect2.pass.vep.vcf.gz
}

# =============================================================================
# 11. Transposable elements
#     Germline MEIs (N + GTEx): MELT.  Somatic L1/Alu/SVA (T, B vs N): xTea.
#     xTea flags differ between releases; generate its run script per its README, e.g.:
#       xtea --case_ctrl --tumor -i ids.txt -b bams.txt -p te/work -o te/submit.sh \
#            -l ref/xtea/rep_lib_annotation -r "$REF" -g ref/gencode.gtf \
#            --xtea <xtea_dir> -y 7 -f 5907
# =============================================================================
melt_germline() {   # sample
  java -Xmx8g -jar MELT.jar Single -bamfile bam/"$1".bqsr.bam -h "$REF" \
      -t ref/melt/mei_list.txt -n ref/melt/hg38.genes.bed -w te/"$1"
}

# =============================================================================
# 12. Clonality: CCF clusters across T and B (PyClone-VI)
# =============================================================================
clonality() {       # patient
  local pt=$1
  python make_pyclone_input.py --vcf vcf/"$pt"/mutect2.pass.vcf.gz \
      --sample "${pt}_T=cnv/${pt}_T" --sample "${pt}_B=cnv/${pt}_B" \
      -o clonality/"$pt".pyclone_in.tsv
  pyclone-vi fit -i clonality/"$pt".pyclone_in.tsv -o clonality/"$pt".h5 \
      -c 40 -d beta-binomial -r 10
  pyclone-vi write-results-file -i clonality/"$pt".h5 -o clonality/"$pt".pyclone_out.tsv
}

# =============================================================================
# 13. Mutational signatures (SBS/ID) per sample
#     Needs the GRCh38 reference installed once for SigProfilerMatrixGenerator.
# =============================================================================
signatures() {      # patients...
  mkdir -p sig/vcfs
  for pt in "$@"; do
    for s in T B; do
      bcftools view -s "${pt}_${s}" vcf/"$pt"/mutect2.pass.vcf.gz \
        | bcftools view -i 'FMT/AD[0:1]>=3' -Ov -o sig/vcfs/"${pt}_${s}".vcf
    done
  done
  python - <<'EOF'
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples="sig/vcfs", output="sig/out",
                   input_type="vcf", genome_build="GRCh38")
EOF
}

# =============================================================================
# main
# =============================================================================
mapfile -t PATIENTS < <(awk -F, 'NR>1 && $1!~/^GTEX/ {print $1}' "$SHEET" | sort -u)
mapfile -t SAMPLES  < <(awk -F, 'NR>1 {print $4}' "$SHEET" | sort -u)
mapfile -t GTEX     < <(printf '%s\n' "${SAMPLES[@]}" | grep '^GTEX' || true)
sex_of() { awk -F, -v p="$1" 'NR>1 && $1==p {print $2; exit}' "$SHEET"; }

while IFS=, read -r pt sex status sm ln fq1 fq2; do
  align_lane "$pt" "$sm" "$ln" "$fq1" "$fq2"
done < <(tail -n +2 "$SHEET")

for sm in "${SAMPLES[@]}"; do markdup_bqsr "$sm"; bam_qc "$sm"; done
identity_check

for pt in "${PATIENTS[@]}"; do germline_gvcf "${pt}_N"; melt_germline "${pt}_N"; done
for g in "${GTEX[@]}"; do germline_gvcf "$g"; melt_germline "$g"; pon_normal "$g"; done
germline_joint
pon_build

for pt in "${PATIENTS[@]}"; do
  mutect2_multi "$pt"
  for sm in "${pt}_T" "${pt}_B" "${pt}_N"; do mito "$sm"; done
  for s in T B; do
    manta_strelka "$pt" "$s"
    delly_somatic "$pt" "$s"
    sv_consensus "$pt" "$s"
    ascat_run "$pt" "$s" "$(sex_of "$pt")"
    msi "$pt" "$s"
  done
  annotate "$pt"
  clonality "$pt"
done
signatures "${PATIENTS[@]}"
