# v1.1.0 阶段 Q 自动验收报告

状态：**2026-09-11 本地自动卡口通过；候选 `v1.1.0-q-candidate.1` 等待 Q→R 人工验收。**

阶段 P 自动通过后已按项目负责人授权直接合入 `release/v1.1.0` 并连续进入 Q，中间没有设置
人工停点。本报告只放行 Q 候选，不授权合并 Q，也不授权开始阶段 R。

## 1. 范围与规格证据

- Q 实现 SSPRK3 候选态容许性检查、troubled-cell 识别、离散支持传播、单调面通量降阶、
  owner 级别一致化、确定性重算和整步缩步重试。
- Case07 的四块 48x32 Mach 5 圆柱纳入生产验收；配置从 `WENO-Z + characteristic` 高阶方案
  起步，未把全域方案改为 `zero_order`。
- `cases/manual/case06_2d_naca0012/` 明确不属于 P/Q。该目录保持用户未跟踪状态，本阶段没有
  读取、修改或加入 Git。
- 《算法补充》11.2.3 是 troubled cell 与传播的权威定义，已补充其功能边界、容许集合、实际
  残差支持、一次转置支持闭包、PH/SCMM6/二点精确 stencil、MPI owner 合并、伪代码与终止性。
  实现中的 `inviscid_residual_stencil(...)` 同时服务残差装配和传播，避免两套支持规则漂移。
- 阶段设计冻结见 [`stage-q-design.md`](stage-q-design.md)。`verify_algorithm_spec.py` 的 6 项
  规格核验全部通过。

## 2. 实现与 Git 审计

Q 阶段生产提交为：

| 提交 | 内容 |
|---|---|
| `b55c882` | 冻结 Q 的稳健化设计 |
| `8288970` | 候选检查、离散支持传播、面级别和 MPI 请求合并 |
| `2895f6a` | SSPRK 局部重算、整步恢复与缩步重试接入运行时 |
| `c4513fe` | 配置、history/manifest 诊断和用户文档 |
| `f0249f6` | 高阶起步的 Case07 Mach 5 稳健配置 |

当前工作分支为 `stage/v1.1.0-q`；P 已以非快进合并 `651ceae` 进入 `release/v1.1.0`。Q 候选
标签创建在本报告对应提交上，但在人工批准前不合并 `release/v1.1.0`。远程推送因当前环境未获
外部网络写入授权而没有执行；本地提交和标签完整，远程状态明确为未同步。

## 3. 自动测试卡口

| 卡口 | 结果 | 关键覆盖 |
|---|---:|---|
| 算法规格核验 | 6/6 通过 | 插值/导数矩、度量、无量纲恒等式、文档契约 |
| Release 串行 CTest | 49/49 通过，59.62 s | 含 robustness、Sod、双马赫、输出/重启和失败路径 |
| Release MPI CTest | 87/87 通过，144.93 s | 含 robustness 的 1/2/4 rank、共享面 owner 合并和无死锁 |
| Case07 生成/拓扑检查 | 通过 | 四块 O 网格与 4-rank 划分 |
| Case07 正式运行 | 通过 | 4 rank，12114 步，精确到达 `t=8` |

`robustness.enabled=false` 是默认值，既有 v1.0 配置仍走原 SSPRK 路径；全量回归没有出现退化。
独立 retry 测试验证 `dt: 1 -> 0.5 -> 0.25`、步首状态精确恢复、只由接受时间步推进时间以及
诊断计数。跨 rank 测试验证 receiver 请求反向传给 donor owner、按最大级别合并，再由正常
通量 halo 正向发布唯一权威通量。

## 4. Case07 迭代与正式结果

开发期间保留了三次失败结论，用于说明保护带设计不是任意调参：

1. profile 通量差分、`CFL=0.08`、最多 4 次 retry，在 `t=0.0798667` 返回数值失败；
2. 改为守恒两点通量差分后，在 `t=0.647820` 返回数值失败；
3. 使用 `CFL=0.02`、最多 8 次 retry，但仅传播直接支持时，在 `t=0.649718` 出现时间步趋零。

因此冻结一次转置离散支持闭包

$$
\mathcal M_0=\bigcup_{i\in\mathcal T}\mathcal S(i),\qquad
\mathcal A=\{j:\mathcal S(j)\cap\mathcal M_0\ne\varnothing\},\qquad
\mathcal M=\mathcal M_0\cup\bigcup_{j\in\mathcal A}\mathcal S(j).
$$

最终配置使用 `scmm6_wcns + conservative_two_point + WENO-Z/characteristic + Rusanov`，
`CFL=0.02`，局部最多重算 3 轮，整步最多 retry 8 次，缩步因子 0.5。4-rank 正式运行信息为：

| 指标 | 结果 |
|---|---:|
| 停止原因 | `physical_time_reached` |
| 步数 / 终止时间 | 12114 / 8.0 |
| 墙钟时间 | 221.379306 s |
| 末步合成 L2 | 4.6916441310e-2 |
| 末态 `rho` 范围 | 0.0188791 -- 4.87646 |
| 末态 `p` 范围 | 5.07856e-4 -- 0.891297 |
| 末态 `T` 范围 | 0.710781 -- 7.13739 |
| 末态 Mach 范围 | 0.0949840 -- 6.20624 |

正式结果的 manifest 记录 `git_commit=f0249f69a680`。归档只保留末态场、完整 history、manifest
和 statistics；检查点、半时刻场及失败试跑目录已在数值核对后清理。

## 5. 局部性、并行一致性与守恒诊断

- history 共 1213 行，其中 469 行抽样出现局部降阶；单步最大低阶 owner 面为 168/9360，
  即 1.7949%，没有全域零阶。
- 抽样最大 level-1/2/3 owner 面数分别为 99/36/156，最大 troubled-cell 数为 51。
- 单步最多 1 次整步 retry；最小接受时间步抽样值为 `8.579138238e-5`。局部重算累计最大 16，
  该数跨三个 RK 阶段及一次 retry 累加，每个 RK 阶段仍受 3 轮上限约束。
- 终态没有触发降阶，候选诊断最小值为 `rho=0.0188791035`、`p=0.000507823885`。
- 串行和 4 rank 在 `t=1` 的 1536 个单元、10 列 Tecplot 值完全相同；第 1760 步诊断、时间和
  `dt` 相同，L2 归约量的最大差为 `2.220446049e-16`。
- 开放远场允许守恒量通过边界。从 `t=0` 到 `t=8`，总质量变化 +0.1353%，总能量变化
  +0.1080%；这些数据用于人工趋势判断，不误报为周期域守恒误差。
- 已归档 [`末态密度`](../../cases/manual/case07_2d_cylinder/figures/mach5-robust-density.png)、
  [`末态 Mach`](../../cases/manual/case07_2d_cylinder/figures/mach5-robust-mach.png) 和
  [`低阶面比例/troubled cell/接受时间步时序`](../../cases/manual/case07_2d_cylinder/figures/mach5-robust-history.png)。
  当前 history 不含面坐标，不能事后重建降阶面的空间图；因此该项明确留在 Q→R 人工卡口，
  不计作自动卡口通过项。

## 6. Q→R 人工卡口

自动卡口结论为通过。人工验收仍需检查 Mach 5 密度/Mach 场、激波位置与厚度、降阶面的时空
分布、最小 `rho/p` 历史、开放边界守恒历史以及与 v1.0 全域零阶基线的差异。人工批准前：

- 不合并 `stage/v1.1.0-q` 到 `release/v1.1.0`；
- 不创建 Q-approved 标签；
- 不建立或实现阶段 R 分支。
