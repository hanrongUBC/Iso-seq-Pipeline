# LAPA 参数整理与重新过滤

不重跑 LAPA，直接用已完成的运行 `lapa_end_flnc_job63782941_20261009_023104_n3qqA0`：

1. 从 LAPA 的计数 BigWig 恢复 9 个样本在每个原始簇（`raw_polyA_clusters.bed`，158,179 个）上的整数计数；
2. 用基因组序列重新判定内部引发：不区分大小写，也不再因为有人类 PAS 六联体就放过；
3. 用明确的逐簇重复规则代替 LAPA 的 1,000 簇滚动重复率；
4. 输出统一的 9 样本计数矩阵，用 DEXSeq 做 3 对 3 检验（WL、DR 各与 Control 比较）。

## 运行

```bash
cd /home/hanrong/scratch/hanrong_LR_stress/APA/LAPA_refilter   # 放这三个文件的目录
sbatch run.sbatch
```

需要 `movapa_westar` 环境（已含 pyBigWig、DEXSeq、data.table）。`FASTA` 必须是 LAPA 运行时用的那份（见运行目录的 `command.sh`）。

## LAPA 原参数与本脚本的替代规则

| 参数 | 本次运行值 | LAPA 默认 | 问题 | 本脚本的处理 |
|---|---|---|---|---|
| `--counting_method` | end | end | poly(A) 尾已剪，只能用 read 末端 | 不变 |
| `--mapq` | 30 | 10 | 去掉 A/C homeolog 上的 reads | 需要重跑 LAPA 才能改；建议另做一次 `--mapq 10` 敏感性分析 |
| `--cluster_window` | 25 | 25 | — | 不变 |
| `--cluster_extent_cutoff` | 3 | 3 | 用九个样本合并后的计数 | 不变（由逐样本重复规则把关） |
| `--cluster_ratio_cutoff` | 0.05 | 0.05 | 同时用于样本内“计数 > 5% 覆盖度”过滤 | 不再使用；直接用原始簇 |
| 内部引发过滤 | 内置，下游 10 nt 中 ≥ 8 个 A 且无 PAS 六联体 | 同左 | 只统计大写 A；有六联体就豁免（人类列表） | 下游 10 nt ≥ 7 个 A，或 20 nt ≥ 12 个 A，或连续 A ≥ 6；不区分大小写，无豁免 |
| `--min_replication_rate` | 0.75 | 0.95 | 滚动窗口比例，不是逐簇规则；前 999 个簇自动通过 | 逐簇：至少一组中 ≥ 2 个重复各有 ≥ 3 条 reads |
| `--replication_rolling_size` | 1000 | 1000 | 同上 | 不再使用 |
| `--replication_num_sample` | 2 | 2 | — | 等同于 `--min-reps 2` |
| `--replication_min_count` | 1 | 1 | 1 条 read 就算“观察到” | `--min-count 3` |
| `--non_replicates_read_threhold` | 10 | 10 | 只对单样本组生效，本实验无影响 | — |
| 差异检验 | 无 | 无 | LAPA 不做组间检验 | DEXSeq 3 对 3；簇层面 padj < 0.05 且 \|ΔPAU\| ≥ 0.1；基因层面用 `perGeneQValue` |

阈值都可以在 `run.sbatch` 里改。

## 输出

- `refilter/clusters_refiltered.tsv.gz`：全部原始簇，含每个样本计数、内部引发得分、各组观察到的重复数、是否在 LAPA 最终结果中、是否保留。
  - `last_base_1based`：最后一个转录碱基（+ 链 = polyA_site，− 链 = polyA_site + 1）。
  - `site_bed_start/site_bed_end`：该碱基的 0-based BED 区间，修正了 LAPA 负链偏 1 bp 的问题。
- `refilter/count_matrix.tsv`：保留簇 × 9 样本计数。
- `refilter/refilter_summary.json`：各步数量，以及与 LAPA 最终 112,249 个簇的重合情况。
- `DEXSeq/cluster_results.tsv`、`gene_results.tsv`、`summary.tsv`：差异 PAS 使用结果。`padj_all_contrasts` 是把两个对比合在一起做 BH 的参考值。

## 注意

- 只检验被唯一分配到一个注释基因的簇；intergenic 簇和跨两个基因的簇不进入 DEXSeq。若要纳入 1,000 bp 延伸规则分配的位点，需要先把它们的 `gene_id` 合并进 `count_matrix.tsv`，并在结果中标明。
- `lapa_refilter.py` 已用模拟数据测试（内部引发、小写碱基、负链坐标、重复规则）；`run_dexseq_apa.R` 未在本地运行（无 R 环境），首次在 Fir 上运行时请检查日志。
