#!/usr/bin/env bash
set -euo pipefail
ROOT=/public3/labmember/qianqh/soil_amp_analysis/10_scAMP_8samples
REF_JOB=${1:?Usage: submit_pipeline.sh REF_JOB_ID}
cd "$ROOT"
printf 'stage\tsample\tjob_id\tdepends_on\n' > logs/job_manifest.tsv
for sample in 1T P S12A SAV FLC SD XM1 XM2; do
  aid=$(qsub -N "ampA_${sample}" -o "$ROOT/logs/02_align_${sample}.log" -v "SAMPLE=$sample" -W "depend=afterok:$REF_JOB" scripts/02_align_sample.pbs)
  printf 'alignment\t%s\t%s\t%s\n' "$sample" "$aid" "$REF_JOB" >> logs/job_manifest.tsv
  cid=$(qsub -N "ampC_${sample}" -o "$ROOT/logs/03_count_${sample}.log" -v "SAMPLE=$sample" -W "depend=afterok:$aid" scripts/03_count_sample.pbs)
  printf 'count\t%s\t%s\t%s\n' "$sample" "$cid" "$aid" >> logs/job_manifest.tsv
done
cat logs/job_manifest.tsv
