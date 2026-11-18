{%- set key = "sarek-cineca" if cookiecutter.project_type == "nf-core/sarek" and not cookiecutter.use_aws_batch else cookiecutter.__pipeline -%}
---
name: estimate-cpu-hours
description: Estimate the CPU-hours (budget) needed to run {{ cookiecutter.project_type }} on the samples of this project, from the number of reads, using a model calibrated on real runs. Use when the user asks how many CPU hours / how much budget the analysis needs, before launching it.
---

# Estimate CPU-hours ({{ cookiecutter.project_type }})

Goal: CPU-hours per process and in total, for the budget request. Cost is the
**requested** cpus x walltime (what SLURM bills), not the CPU actually used.
The model `conf/cpu_model.json` (one key per pipeline, here `{{ key }}`) was fitted with
`scripts/estimate_cpu_hours.py fit` on real runs: CPU-h scale linearly with the number of reads.
{%- if cookiecutter.project_type == "nf-core/methylseq" %}
Two RRBS runs (98 samples, bismark, rat genome): ~1 CPU-h per million read pairs,
80% of it in `BISMARK_ALIGN`; per-sample error ~5%; a model fitted on one run predicted the other within 2%.
{%- elif cookiecutter.project_type == "cnr-ibba/nf-resequencing-mem" %}
One run on CINECA (10 sheep samples, 2x100, 3078 M read pairs, v0.5.2): 1626 CPU-h, ~0.53 CPU-h per million read pairs,
64% in `BWA_MEM`, then `FREEBAYES_CHUNK` (15%), `TRIMGALORE`, `PICARD_MARKDUPLICATES`. Single run: the per-sample error
(6%) only says how the 10 samples scatter, not how another genome or number of samples behaves.
{%- elif cookiecutter.project_type == "nf-core/sarek" %}
{%- if not cookiecutter.use_aws_batch %}
Key `sarek-cineca`: one CINECA run of sarek 3.4.0 (10 sheep samples, 2x100, 3078 M read pairs, BQSR + joint calling):
9447 CPU-h, ~3.1 CPU-h per million read pairs, 76% in `GATK4_HAPLOTYPECALLER` (requested 6-12 cpus, but only ~1.4 cores used on average: `%cpu` of the trace),
11% in `BWAMEM1_MEM`, 8% in joint calling (`GENOTYPEGVCFS`, `GENOMICSDBIMPORT`, fixed per run).
Key `sarek` (AWS, bean 521 Mb, 2x150, HaplotypeCaller with 2 cpus) predicts the sheep run 4x too low (2294 vs 9447):
alignment and MarkDuplicates agree per sequenced base (BWA 1.6-1.7 CPU-h/Gbp at both sites), the variant calling does not
(genome size, cpus requested, CPU speed: AWS wall time per Mb of genome was ~2x shorter). Don't use the AWS key on another
genome. The CPU-h billed on CINECA are the requested ones, so capping `GATK4_HAPLOTYPECALLER` at 2 cpus in `conf/custom-*.config`
(as in the AWS runs) would cut the largest item (not tested here).
{%- endif %}
Nine AWS batches (541 samples, bwa-mem + HaplotypeCaller + joint germline calling, *Phaseolus vulgaris* 521 Mb, 2x150,
pool-seq, AWS Batch): ~0.75 CPU-h per million read pairs (~2.2-2.9 CPU-h per Gbp),
70% in `BWAMEM1_MEM`, 22% in `GATK4_HAPLOTYPECALLER`. Held-out batches are predicted within ~10-18%.
{%- else %}
Two runs (STAR + Salmon, paired-end 150, 129 samples, bos taurus and rattus norvegicus): ~0.35 CPU-h per million
read pairs (STAR_ALIGN, DUPRADAR, TRIMGALORE, PICARD_MARKDUPLICATES dominate). The total of a run is predicted
within ~10% (a model fitted on one species predicted the other within 12-24%, the pooled one within 6%),
but **per-sample cost varies a lot** (0.06-1.2 CPU-h per million reads), so trust the total, not single samples.
{%- endif %}

## 1. Reads

The predictor is the number of **read pairs** (not genome size x coverage).
{%- if cookiecutter.project_type == "nf-core/sarek" %} Alignment and variant calling scale with
the sequenced bases (genome x coverage), but the calibration has a single genome and read length (150):
for another genome or read length, scale the estimate by bases (read pairs x read length), not by pairs.
{%- endif %}
Use R1 `num_seqs` from `conf/seq_stats.tsv` (see the `split-samplesheet` skill; run it
first if missing) and sum it in millions:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $h["file"] ~ /_(R)?1[._]/ {s+=$h["num_seqs"]} END{print s/1e6}' conf/seq_stats.tsv
```

(check the R1 pattern on the file names; for single-end use all files).

## 2. Estimate

```bash
python scripts/estimate_cpu_hours.py predict --pipeline {{ key }} --million-reads <M>   # FASTQ read pairs
```

`--no-genome-prep` if the genome index already exists.
{%- if cookiecutter.project_type == "nf-core/sarek" %} On AWS add `--usd-per-cpu-h 0.0126 0.0188`
(EC2 only, all services and taxes; from the consolidated bill of the calibration runs: 41k CPU-h, $772, $1.33 per sample)
to print the cost range; the model counts only successful tasks, so for the billed total use the upper price
or add the retries of Spot losses.
{%- endif %}
{%- if cookiecutter.project_type == "nf-core/methylseq" %} `--aligned` if the
number is already after trimming.
{%- endif %} Report the table and the total.

## 3. Caveats to tell the user

{% if cookiecutter.project_type == "nf-core/methylseq" -%}
- In RRBS the trimming length filter drops reads: the model applies the average
  retention of the calibration runs (0.88), but it was 0.17-0.98 per sample, so
  libraries with poor quality/adapter content cost less than predicted.
- Calibrated on **RRBS** with the `cpus` requested in those runs (bismark_align 24
  cpus); WGBS (e.g. 30X) or different resources need a new calibration, so say the
  estimate is only indicative. To recalibrate: `fit --pipeline methylseq --run
  pipeline_info/execution_report_*.html <multiqc_data>/bismark_alignment.txt
  <multiqc_data>/multiqc_fastqc.txt`
{%- elif cookiecutter.project_type == "nf-core/sarek" -%}
- Key `sarek` = AWS: calibrated on pool-seq of a 521 Mb genome with the `cpus` requested in those runs (BWA 24, MarkDuplicates 6, HaplotypeCaller 2);
  per-batch error is systematic: 2.2 CPU-h/Gbp in the BMK batches, 2.8 in the older ones (CPU type, retries or data
  quality cannot be separated), so give a +/-15% range. Joint calling (`GENOMICSDBIMPORT`, `GENOTYPEGVCFS`) is a
  fixed cost per run (~35 CPU-h) that grows with the number of samples, not with the reads.
  To recalibrate: `fit --pipeline <key> --genome-mb <sum of the .fai> --read-length <bp> --calibration "<text>" --run
  pipeline_info/execution_report_*.html <multiqc_data>/multiqc_fastqc.txt` (use a new key per site/version, e.g. `sarek-cineca`; lanes are summed per sample; a sample sequenced again in another run is kept apart)
{%- elif cookiecutter.project_type == "cnr-ibba/nf-resequencing-mem" -%}
- Calibrated on a single run (10 sheep samples, 2x100, CINECA, 2.6 Gb genome): `BWA_MEM` scales with the bases, but
  `FREEBAYES_CHUNK` (one task per genome chunk, joint over all samples) is counted as fixed: it grows with the genome size and the
  number of samples, which a single run cannot separate, so it is only indicative for another genome.
  To recalibrate: `fit --pipeline nf-resequencing-mem --genome-mb <sum of the .fai> --read-length <bp> --run
  pipeline_info/execution_report_*.html <multiqc_data>/multiqc_fastqc.txt`
{%- else -%}
- Calibrated on mammalian PE150 samples with the `cpus` requested in those runs (STAR 12 cpus
  or similar); different read length, genome size or resources make the estimate only indicative.
  To recalibrate: `fit --pipeline rnaseq --run pipeline_info/execution_report_*.html
  <multiqc_data>/multiqc_fastqc_fastqc_raw.txt`
{%- endif %}
  (repeat `--run` to pool runs; sample names must be unique; only the `--pipeline` key is
  rewritten). The report HTML has `cpus`, the trace .txt does not. After a `-resume`, use the
  last report (cached tasks keep their realtime); failed retries are not counted.
- Check the result against the budget: round up (+20%) for reruns of failed tasks.
