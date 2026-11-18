#!/usr/bin/env python
"""Estimate the CPU-hours (requested cpus x realtime) of a nf-core/methylseq, rnaseq, sarek or nf-resequencing-mem run.

fit:     calibrate CPU-h per million reads, per process, on a finished run
         methylseq: execution_report_*.html + MultiQC bismark_alignment.txt + multiqc_fastqc.txt
         rnaseq:    execution_report_*.html + MultiQC multiqc_fastqc_fastqc_raw.txt (raw reads)
         sarek, nf-resequencing-mem: execution_report_*.html + MultiQC multiqc_fastqc.txt (raw reads; sarek has one row per lane)
         --pipeline is the key in cpu_model.json: use one key per site and version (e.g. sarek = AWS, sarek-cineca),
         because the requested cpus of each process differ; sarek* keys read sarek reports
predict: apply conf/cpu_model.json to the number of reads of a new run
"""

import re
import csv
import json
import argparse
import collections
import pathlib

CPU_RE = re.compile(r'"name":"([^"]+)","status":"[^"]*","exit":"[^"]*","submit":"\d+",'
                    r'"start":"\d+","complete":"\d+","duration":"\d+","realtime":"(\d+)".*?"cpus":"(\d+)"')
GENOME_PREP = ("BISMARK_GENOMEPREPARATION_BOWTIE", "STAR_GENOMEGENERATE", "SALMON_INDEX", "BWAMEM1_INDEX")
RESEQ_RE = re.compile(r"_T\d+$")  # nf-resequencing-mem: MultiQC names are <sample>_T1_1, tasks <sample>
LANE_RE = re.compile(r"[-_]lane.*$")  # sarek: FASTQC/FASTP tasks and FastQC rows are per lane, the other tasks per sample
TASK_RE = re.compile(r"^(?:\S+:)*(\w+) \((.+)\)$")


def read_tasks(report, suffix=None):
    """yield (process, sample or None, requested cpus x realtime in hours)"""
    for name, realtime, cpus in CPU_RE.findall(pathlib.Path(report).read_text()):
        m = TASK_RE.match(name)
        if m:
            yield m.group(1), (suffix.sub("", m.group(2)) if suffix else m.group(2).rstrip("_")), int(cpus) * int(realtime) / 3.6e6
        else:
            yield name.split(":")[-1], None, int(cpus) * int(realtime) / 3.6e6


def read_millions(bismark_alignment):
    """reads given to bismark per sample, in millions (R1/R2 split names merged)"""
    reads = collections.defaultdict(float)
    with open(bismark_alignment) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sample = re.sub(r"_\d+$", "", re.sub(r"_val_\d+$", "", row["Sample"])).rstrip("_")
            reads[sample] += sum(float(v) for k, v in row.items() if k != "Sample" and v) / 1e6
    return reads


def read_raw_millions(multiqc_fastqc, suffix=None):
    """raw read pairs per sample, in millions (R1 rows of multiqc_fastqc.txt; the suffix regex, e.g. the lane, is removed and summed)"""
    reads = collections.defaultdict(float)
    with open(multiqc_fastqc) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Sample"].endswith("_1"):
                name = r["Sample"][:-2]
                reads[suffix.sub("", name) if suffix else name.rstrip("_")] += float(r["Total Sequences"]) / 1e6
    return reads


def fit(runs, out, pipeline, calibration=None, genome_mb=None, read_length=None):
    """runs: list of [report, files...], pooled; sample names must be unique across runs.
    methylseq files: bismark_alignment.txt, multiqc_fastqc.txt; rnaseq files: multiqc_fastqc_fastqc_raw.txt"""
    sarek = pipeline.startswith("sarek")
    suffix = LANE_RE if sarek else RESEQ_RE if pipeline == "nf-resequencing-mem" else None
    reads, raw = {}, {}
    per_sample = collections.defaultdict(dict)  # process -> sample -> cpu-h
    fixed = collections.defaultdict(float)
    for i, (report, *files) in enumerate(runs):
        if sarek or pipeline in ("rnaseq", "nf-resequencing-mem"):  # predictor = raw read pairs, known before the run (no retention)
            assert len(files) == 1, f"{pipeline}: --run REPORT multiqc raw-reads file"
            run_raw = read_raw_millions(files[0], suffix)
            run_reads = dict(run_raw)
        else:
            assert len(files) == 2, "methylseq: --run REPORT bismark_alignment.txt multiqc_fastqc.txt"
            run_reads, run_raw = read_millions(files[0]), read_raw_millions(files[1])
        if sarek:  # a sample can be resequenced in a later batch: keep its runs apart
            run_reads = {f"{i}:{k}": v for k, v in run_reads.items()}
            run_raw = dict(run_reads)
        assert not set(run_reads) & set(reads), "sample names repeated across runs"
        reads.update(run_reads)
        raw.update(run_raw)
        for proc, sample, cpuh in read_tasks(report, suffix):
            if sarek and sample:
                sample = f"{i}:{sample}"
            if sarek and proc == "BWAMEM1_MEM":
                assert sample in run_reads, f"{sample}: task without FastQC reads in {files[0]}"
            if run_reads.get(sample, 0) >= 0.1:  # empty FASTQ (<0.1 M reads) would break the per-read fit: counted as fixed
                per_sample[proc][sample] = per_sample[proc].get(sample, 0) + cpuh
            else:
                fixed[proc] += cpuh / len(runs)  # once-per-run processes: mean over runs

    def coef(procs, skip=None):
        num = sum(per_sample[p][s] * reads[s] for p in procs for s in per_sample[p] if s != skip)
        den = sum(reads[s] ** 2 for p in procs for s in per_sample[p] if s != skip)
        return num / den

    model = {"unit": "CPU-h (requested cpus x realtime) per million reads", "n_samples": len(reads),
             "per_million_reads": {p: coef([p]) for p in per_sample}, "fixed_cpu_h": dict(fixed)}

    # leave-one-out relative error of the per-sample total
    samples = sorted(set(s for p in per_sample.values() for s in p))
    errs = []
    for s in samples:
        pred = sum(coef([p], skip=s) * reads[s] for p in per_sample if s in per_sample[p])
        true = sum(per_sample[p][s] for p in per_sample if s in per_sample[p])
        errs.append((pred - true) / true)
    model["loo_sample_rel_err_rmse"] = (sum(e * e for e in errs) / len(errs)) ** 0.5
    total_true = sum(sum(p.values()) for p in per_sample.values())
    model["loo_total_rel_err"] = abs(sum(
        coef([p], skip=s) * reads[s] for s in samples for p in per_sample if s in per_sample[p]
    ) - total_true) / total_true
    # trimming (RRBS length filter) drops reads: fraction of FASTQ reads reaching bismark
    if pipeline == "methylseq":
        ratios = [reads[s] / raw[s] for s in samples if s in raw]
        model["retention_total"] = sum(reads[s] for s in samples) / sum(raw[s] for s in samples if s in reads)
        model["retention_range"] = [min(ratios), max(ratios)]
    if calibration:
        model["calibration"] = calibration
    if genome_mb:  # informational: the per-process coefficients hold for this genome and read length only
        model["genome_mb"], model["read_length"] = genome_mb, read_length
    path = pathlib.Path(out)
    models = json.loads(path.read_text()) if path.exists() else {}
    models[pipeline] = model  # other pipelines in the file are kept
    path.write_text(json.dumps(models, indent=2) + "\n")
    print(f"{len(runs)} run(s), {len(samples)} samples, {total_true:.0f} per-sample CPU-h; "
          f"LOO error per sample {model['loo_sample_rel_err_rmse']:.1%}, total {model['loo_total_rel_err']:.1%}")


def predict(model_path, pipeline, million_reads, genome_prep, raw, usd=None):
    model = json.loads(pathlib.Path(model_path).read_text())[pipeline]
    if raw and "retention_total" in model:  # FASTQ read pairs (seqkit num_seqs of R1) -> reads reaching bismark
        million_reads *= model["retention_total"]
        lo, hi = model["retention_range"]
        print(f"reads reaching bismark = FASTQ x {model['retention_total']:.2f} (per sample {lo:.2f}-{hi:.2f} in the calibration run)")
    rows = {p: c * million_reads for p, c in model["per_million_reads"].items()}
    for p, c in model["fixed_cpu_h"].items():
        # genome preparation runs once per genome, and can be skipped if the index exists
        if p not in GENOME_PREP or genome_prep:
            rows[p] = c
    total = sum(rows.values())
    for p, v in sorted(rows.items(), key=lambda i: -i[1]):
        print(f"{p:32s}{v:10.1f}")
    if "genome_mb" in model:
        print(f"calibrated on a {model['genome_mb']:.0f} Mb genome, {model['read_length']} bp reads: genome-bound processes (variant calling) do not "
              f"scale with the reads alone, a different genome or read length is only indicative")
    err = model["loo_sample_rel_err_rmse"] / max(model["n_samples"], 1) ** 0.5
    print(f"{'TOTAL CPU-h':32s}{total:10.1f}  (+/- {total * max(err, 0.05):.0f} from the fit alone; calibrated on {model.get('calibration', 'the fitted runs')} only{'; add the read-retention uncertainty' if 'retention_total' in model else ''})")
    if usd:  # price per requested CPU-h, e.g. from Cost Explorer / total CPU-h of past runs
        print(f"{'COST':32s}{total * usd[0]:10.0f} - {total * usd[-1]:.0f} USD  ({usd[0]} - {usd[-1]} USD per CPU-h)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--run", nargs="+", action="append", required=True, metavar="FILE",
                   help="execution_report_*.html followed by the MultiQC files of the pipeline (see top); repeat to pool runs")
    f.add_argument("--calibration", help="text shown by predict, e.g. 'RRBS' or 'rnaseq STAR+Salmon, bos taurus + rattus norvegicus'")
    f.add_argument("--out", default="conf/cpu_model.json")
    f.add_argument("--genome-mb", type=float, help="genome size of the fitted runs, in Mb (sum of the .fai), stored for the warning of predict")
    f.add_argument("--read-length", type=int, help="read length (bp) of the fitted runs")
    for sp in (f, p := sub.add_parser("predict")):
        sp.add_argument("--pipeline", default="methylseq", help="key in cpu_model.json")
    p.add_argument("--million-reads", type=float, required=True, help="total read pairs of all samples (R1 num_seqs), in millions")
    p.add_argument("--aligned", action="store_true", help="the number is already the reads given to bismark (after trimming)")
    p.add_argument("--model", default="conf/cpu_model.json")
    p.add_argument("--usd-per-cpu-h", type=float, nargs="+", metavar="USD", help="price per CPU-h (one value, or low high) to print the cost")
    p.add_argument("--no-genome-prep", action="store_true", help="bismark index already available")
    a = ap.parse_args()
    if a.cmd == "fit":
        fit(a.run, a.out, a.pipeline, a.calibration, a.genome_mb, a.read_length)
    else:
        predict(a.model, a.pipeline, a.million_reads, not a.no_genome_prep, not a.aligned, a.usd_per_cpu_h)
