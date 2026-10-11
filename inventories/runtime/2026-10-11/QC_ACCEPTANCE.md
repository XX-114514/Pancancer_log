# QC驗收結果

**結論：4份分片、993,312個細胞通過全量技術QC驗收；可以進行來源內前處理。正式新增來源驗收仍為0，原研究別名與基線文庫重疊審查尚未閉環。**

完整計數審計完成：2026-10-11T13:15:27.440923+00:00（2026-10-11 21:15:27 UTC+8）；登記：2026-10-11T13:18:52.909137+00:00（21:18:52 UTC+8）。

## 逐來源結果

| 來源 | 原始細胞 | QC保留 | 保留率 | 文庫 | 已連患者 | 判定 |
|---|---:|---:|---:|---:|---:|---|
| CAR-T（GSE241783 / CXG） | 417,167 | 416,800 | 99.91% | 59 | 59 | 技術PASS，保留用途限制 |
| GSE286827 | 416,211 | 384,628 | 92.41% | 57 | 29 | 技術PASS，保留用途限制 |
| GSE292046 | 55,570 | 54,156 | 97.46% | 8 | 2 | 技術PASS，保留用途限制 |
| GSE294109 | 145,316 | 137,728 | 94.78% | 25 | 15 | 技術PASS，保留用途限制 |

患者數為各來源內公開映射的非空代碼數，不能直接跨來源相加當作獨立患者數。頭頸癌須使用patientfix_v2覆蓋表：29患者、57文庫已連，原v1缺失的7,463細胞已修正。

## 本次實際驗收內容

- 全部993,312個保留細胞、全部稀疏非零計數逐塊對原始H5AD，完全一致；不是抽樣。
- X與counts層全部數值/結構一致；有限、非負、整數；源H5AD和官方metadata SHA與既有憑證一致。
- 穩定ID來源位置/原始barcode、原始ID順序、跨分片唯一性、基因順序/唯一性核對通過。
- 原始與保留账本行數、149文庫registry/阈值集合、原始數合計一致；每文庫floor和median+3MAD線粒體阈值重新計算通過。
- counts重新計算的UMI、基因數、線粒體比例與ledger/obs一致；保留細胞均過floor/各自mt門檻，doublet_score有限、標記singlet。
- 原v1的sample/patient/treatment/timepoint/tissue/clinical_tier等字段與registry一致；修復覆蓋表另獨立驗收。

全部4×36項技術條件PASS。用5000行塊、單線程執行；沒有重算Scrublet或追加刪除細胞。總過濾40,952：低基因/UMI5,919、高mt33,103、doublet1,930（原因為流程最終記錄，不是互相獨立的原始超閾值計數）。

## Doublet與分布限制

| 來源 | 刪除doublet | mt之後評估細胞 | 實際刪除率 | 零檢出文庫 | 保留基因數中位 | 保留UMI中位 |
|---|---:|---:|---:|---:|---:|---:|
| CAR-T | 103 | 416,903 | 0.025% | 20/59 | 2605 | 8186 |
| GSE286827 | 289 | 384,917 | 0.075% | 28/57 | 1622 | 4630 |
| GSE292046 | 1,352 | 55,508 | 2.436% | 0/8 | 3855.5 | 10229 |
| GSE294109 | 186 | 137,914 | 0.135% | 7/25 | 1699 | 3996 |

CAR-T、頭頸癌、腎癌的doublet檢出很低。這是觀測結果；作者既有QC、細胞組成或方法敏感性均可能影響，尚未定位唯一原因。不能把computed或singlet標記當作雙細胞完全清除的證明，也不應為湊預期4%而硬刪。適合先來源內探索，正式組成/狀態結果仍需原作者QC與marker/深度檢查。

## 已處理的患者映射問題

頭頸癌一個官方GEX標題含post2，舊解析規則未識別。依官方標題的患者token與同來源已有前後樣本的唯一patient_uid建立v2 registry與stable-ID元數據覆蓋表，修復7,463細胞；患者數仍29，不製造新患者。timepoint保留官方字面code，不推斷時間間隔。response、生存、治療及tier沒有改動。

獨立核對用split/token查找而非producer regex，全部overlay ID精確對應原QC該文庫，無額外/重複ID，覆蓋後未映射文庫為0。原矩陣、v1 registry與舊audit保留。新版registry和細胞行只在伺服器，GitHub僅有SHA及匯總。

## 用途與技術分層

- 腫瘤組織QC共576,512；CAR-T產品QC416,800單列。CAR-T原研究為[GSE241783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241783)的輸注產品，不能冒充原位腫瘤生態。
- [GSE292046](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292046)是T細胞富集snRNA，2患者前後+重複共8文庫；不得當8獨立患者，也不能用富集比例直接代表全組織細胞組成。
- [GSE286827原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC12629806/)與[GSE294109原研究](https://doi.org/10.1016/j.ccell.2026.05.007)追溯明確；本次核對不宣稱已關閉54基線研究全部別名/文庫排重。
- clinical_sidecar_sha256=null表示本流程沒有使用外部TSV，preflight明確用官方H5AD/metadata；不能單憑null認定無臨床來源。仍缺患者結局者保持C層，不由標題或治療推斷A/B。

## 下一步

1. 以patientfix_v2覆蓋表加入頭頸癌元數據；四來源分別做normalize/log、來源內圖/marker審查，counts保持原始。
2. 三個低doublet檢出來源核查作者歷史QC與高深度/跨譜系混合marker，按文庫審查，不重複硬過濾。
3. 完成原研究別名和原始文庫對基線重疊審查，簽發正式source acceptance後才增加450萬—550萬目標數。
4. 跨研究整合仍待用户決定；本輪沒有啟動新注釋、CNV或模型。

## 證據

- [完整數值審計](independent_qc_acceptance_v1.json)、[技術验收登記](technical_QC_acceptance_v1.json)。
- [QC分布/舊患者映射背景](qc_context_review_v1.json)；其中旧缺失數由v2覆蓋修正，不重寫歷史。
- [v2元數據修復](qc_metadata_overlay_acceptance_v2.json)、[獨立覆核](qc_independent_overlay_audit_v1.json)。
- [正式run](../../../runs/20261011_expansion_QC_acceptance_v1.md)、[任務日志](../../../logs/codex/2026/2026-10-11_QC-acceptance.md)。
