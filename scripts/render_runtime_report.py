#!/usr/bin/env python3
"""Render a readable, timestamped report from an aggregate-only snapshot."""
import argparse
import json
from pathlib import Path
from validate_runtime_snapshot import validate

def render(s, snapshot_name="runtime_snapshot.json"):
    errors=validate(s)
    if errors:raise ValueError("; ".join(errors))
    e=s["expansion"];v=s["project_v3"]
    text=f"# 下载与Project_v3运行快照\n\n采集时间（UTC）：**{s['collected_utc']}**。这是一次实测快照，不是GitHub实时监控。\n\n"
    text+=f"## 新公共数据收集\n\n- 已验证下载/源矩阵库存：**{e['source_cells_verified_download_or_inventory']:,}源细胞**。\n- QC分片：**{e['qc_shards']}份，{e['qc_cells_pending_final_source_acceptance']:,}细胞**，仍待最终来源验收。\n- 正式验收新增：**{e['final_accepted_new_cells']:,}**；目标为排重、QC后450万—550万。\n- 用户已授权下载与逐来源前处理；跨研究整合等待单独决定。\n\n"
    text+="| 队列 | 源细胞 | 已完成的处理 |\n|---|---:|---|\n"
    for d in e["datasets"]:
        parts=[]
        if d.get("QC"):parts.append(f"QC保留{d['QC']['n_qc_nonduplicate']:,}；待最终来源验收")
        if d.get("native_preprocessing"):parts.append("原始计数Zarr读回验证完成；待QC/来源验收")
        if d.get("baseline_overlap"):parts.append(f"完整计数对基线精确匹配{d['baseline_overlap']['baseline_exact_matches']}个；不能排除不同基因覆盖下的重叠")
        if not parts:parts.append("已下载；来源/文库/QC策略审计中")
        text+=f"| {d['name']} | {d['source_cells']:,} | {'；'.join(parts)} |\n"
    text+="\n### 第二批文件传输\n\n| 队列 | 已校验文件 | 未完成文件进度 |\n|---|---:|---|\n"
    for d in e["batch2"]:
        done=sum(f["download_verified"] for f in d["files"]);transfers=[]
        for f in d["files"]:
            if f["download_verified"]:continue
            t=f["transfer"]
            label=f"{t.get('percent','未知')}%，{t.get('rate_per_second','未知')}/s" if "percent" in t else "等待/未观测到传输"
            if t.get("stale_over_10min"):label+="（日志超过10分钟未更新，需检查作业）"
            transfers.append(f["name"]+"："+label)
        text+=f"| {d['dataset']} | {done}/{len(d['files'])} | {'<br>'.join(transfers) or '文件校验完成，计数与临床验收另计'} |\n"
    z=e["zenodo"];t=z["transfer"]
    text+=f"\nZenodo 19586163：{t.get('percent','未知')}%，{t.get('rate_per_second','未知')}/s；官方校验完成标记：{bool(z['download_proof'])}。传输进度不作为完成证明。\n"
    text+="\n### 当前审计限制\n\n- TNBC的7个orig.ident是分组，不能作为原始文库；472种barcode前缀也尚未获得作者文库映射。\n- CTCL的13个来源需追溯原研究别名、文库及基因覆盖；570,282个非primary标记细胞处于来源隔离审计。\n- CAR-T为治疗产品，原位肿瘤生态资格单列。临床缺失不从表达矩阵推断。\n- CRC元数据有975,275条细胞记录，但计数文件未完整校验时不计入上面的完成下载规模。\n"
    text+="\n## Project_v3当前动态\n\n后续独立审计和最终交接登记优先于生产者早期manifest的pending检查点；保留早期状态用于追溯。\n\n"
    text+=f"源状态表生成时间：`{v['source_state_generated_utc']}`；本次重新读取了当前登记和轻量审计文件，并保存SHA-256。报告中的PASS是服务器已有审计的记录，本次没有重算数百万细胞。\n\n"
    state=v["identity_release_readiness"]
    text+=f"- 核心范围：{state.get('full_cells',0):,}细胞；登记细身份非空{state.get('current_nonempty_fine_cells',0):,}。\n- 最终科学身份发布：`{state.get('status')}`；GSE299623仍有{state.get('remaining_GSE299623_groups',0)}群/{state.get('remaining_GSE299623_cells',0):,}细胞未获最终身份决定。\n- 表达注释与新下游交接：`{v['expression_annotation'].get('artifact_task_status','未登记')}`。描述性单基因命名保留既有身份，不自动成为新细胞亚型。\n\n"
    text+="| 运行/证据 | 已记录状态或规模 |\n|---|---|\n"
    for r in v["run_evidence"]:
        if not r["file"]["exists"]:continue
        summary=r["summary"];parts=[str(summary.get("status",summary.get("canonical_write_status","已登记")))]
        for key in ["included_rows","samples_min100","samples","nodes","K","all_LIANA_rows","within_CM_rows","reviewed_samples","reviewed_patient_units","primary_model_BH_significant_strata"]:
            if key in summary:parts.append(key+"="+str(summary[key]))
        text+=f"| {r['run_id']}/{r['artifact']} | {'；'.join(parts)} |\n"
    scientific=[p for p in s["observed_processes"] if p["scope"]=="Project_v3" and p.get("role")=="pipeline"]
    text+=f"\n本次观察到Project_v3科学/流水线进程：**{len(scientific)}**。代理、日志刷新等基础进程不计为科学计算。进程存在不等于结果验收；观察范围为当前用户可读进程与所登记路径。\n"
    text+="\n### 下一步\n\n1. 补齐CM候选来源的患者、癌种、组织、分选和终点映射；六个可估计主模型分层没有BH显著结果，不升级为已验证机制。\n2. 复核K与来源/分母敏感性；9个CM是当前模型，不声称已证明全局最优。\n3. 新临床队列完成来源、文库、QC和临床端点验收后，冻结独立验证集；不在验证集中重新选择marker或重拟合发现模型。\n4. GSE299623身份审阅单独闭环；表达状态、细胞身份、CNV恶性和生态/临床验证分别登记。\n"
    text+=f"\n完整机器快照与文件哈希：[{snapshot_name}]({snapshot_name})。\n"
    return text

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("snapshot",type=Path);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error("Preserve prior report; choose a new output filename")
    a.output.write_text(render(json.loads(a.snapshot.read_text(encoding="utf-8")), a.snapshot.name),encoding="utf-8")
if __name__=="__main__":main()
