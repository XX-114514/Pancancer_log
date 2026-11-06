#!/bin/bash
# E2 脚本的合成数据验证。用法: bash run_tests.sh <可写的临时目录>
set -u; S=$1; here=$(cd "$(dirname "$0")" && pwd); cd "$here"
rm -rf "$S"; mkdir -p "$S"
Rscript --vanilla make_synthetic.R "$S/data" >/dev/null 2>&1
run() { E2_SYN="$S/data" E2_OUT="$1" E2_AXES="$2" E2_B="$3" Rscript --vanilla ../e2_axis_increment.R config_synthetic.R "$4" 2>&1; }
echo "### T1 未过复现门时拒绝计算新效应"
run "$S/t1" axes.tsv 50 null | grep -E 'repro has not been run' >/dev/null && echo PASS || echo FAIL
echo "### T2 已发表值被改动时复现门拦截"
mkdir -p "$S/bad"; cp "$S"/data/* "$S/bad/"
awk -F'\t' 'BEGIN{OFS="\t"} NR==2{$4=$4*1.05} 1' "$S/data/published.tsv" > "$S/bad/published.tsv"
E2_SYN="$S/bad" E2_OUT="$S/t2" Rscript --vanilla ../e2_axis_increment.R config_synthetic.R repro 2>&1 | grep -E 'repro FAILED' >/dev/null && echo PASS || echo FAIL
echo "### T3 完整流程(埋入答案)"
run "$S/out" axes.tsv 2000 all | sed -n '/axis /,$p' | cut -c1-170
echo "### T4 输出不覆盖"
run "$S/out" axes.tsv 50 repro | grep -E 'refuse to overwrite' >/dev/null && echo PASS || echo FAIL
echo "### T5 假阳性校准: 30 条随机轴"
mkdir -p "$S/cal"; cp "$S"/out/repro_gate.tsv "$S"/out/repro_mode.txt "$S"/out/gene_marginals.tsv "$S/cal/"
run "$S/cal" axes_calibration.tsv 500 null >/dev/null
Rscript --vanilla -e 'x<-read.delim(commandArgs(TRUE)[1]);print(aggregate(p_emp~family,x,function(p) c(frac05=mean(p<.05),median=median(p))))' "$S/cal/null_summary.tsv"
