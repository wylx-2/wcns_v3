# v1.1.0 阶段 S 自动验收报告

状态：**2026-09-11 自动卡口通过；按项目负责人 R→S→T 连续推进授权，S 不设置单独人工
停点，候选 `v1.1.0-s-candidate.1` 将直接非快进合入 `release/v1.1.0` 并进入 T。**

## 1. 范围、接口与数学证据

- [`stage-s-design.md`](stage-s-design.md) 冻结 `constant`/`sutherland`、Prandtl 数、Kelvin 与
  无量纲 Sutherland 温度二选一、默认迁移、合法域和重启签名语义。
- 根目录《算法补充》11.5 冻结
  `mu=mu_ratio`、`mu=mu_ref_ratio*T^(3/2)*(1+S)/(T+S)` 和
  `chi=mu/[(gamma-1) Ma^2 Pr]`；`1/Re` 没有移入输运模型，仍只在粘性散度和时间步限制中
  各使用一次。
- `CaseConfig` 是输运参数的唯一所有者；求解器、场输出和边界输出由同一个
  `TransportConfig` 创建模型。启动摘要打印模型、Pr、`mu(1)`、`S/T_ref` 和初始化/恢复后
  全局真实单元温度、黏度范围。
- 未给 `transport.*` 的 v1.0 输入映射为 `constant, viscosity_ratio=1, Pr=0.72`；只有该精确
  默认配置可读取旧签名检查点。任一有效输运参数改变均在推进前拒绝恢复。
- `cases/manual/case06_2d_naca0012/` 不属于本阶段，保持用户未跟踪状态；本阶段没有读取、修改
  或加入 Git。

## 2. 实现与 Git 审计

| 提交 | 内容 |
|---|---|
| `9c45ba2` | 冻结 S 输运配置、数学、迁移和验收设计 |
| `8c9f51e` | 将输运模型接入配置、求解器、输出、摘要和重启签名 |
| `b0444ce` | 验证解析公式、合法单位等价及严格配置失败路径 |
| `c963d3a` | 在启动输运扫描前同步初始化后的 primitive 温度 |
| `c2c2f89` | 增加生产 Sutherland 串行及 2-rank 运行测试 |
| `22d1474` | 增加改变输运参数时的检查点拒绝测试 |
| `5002947` | 更新完整模板、运行指南、用户手册和 README |
| `11b126a` | 据前验冻结内部、壁带和全局空间阶卡口 |
| `28fc91b` | 增加两 profile、二维/三维、三层网格的粘性截断误差矩阵 |
| `fa417ef` | 增加四档时间步 SSPRK3 自收敛卡口 |
| `e1e84ed` | 允许粘性解析矩阵选择生产输运模型 |
| `c04e2f1` | 增加经验 CFL 分类扫描与稳健性诊断 |
| `ffe3c5c` | 生成 v1.0 默认输运签名夹具并验证迁移/拒绝策略 |
| `02d982a` | 保存空间、时间和 CFL 机器可读原始证据 |

远程推送没有执行：此前外部网络写入授权未获批准。本地分支、提交和候选标签完整，远程状态
明确为未同步。

## 3. 公式、配置与重启卡口

| 卡口 | 结果 |
|---|---:|
| 常黏度、Sutherland 点值与 `T=1` 归一化 | 相对误差不超过 `1e-14` |
| Kelvin/无量纲 Sutherland 常数等价 | 规范化签名相同 |
| `chi`、Pr 和粘性时间步尺度 | 相对误差不超过 `1e-13`，`1/Re` 未重复使用 |
| 未知模型、错配模型键、重复/双温度表示、非正/非有限参数 | 全部在推进及正常输出前拒绝 |
| 生产 Sutherland 串行/2-rank | 通过，输出的温度与黏度范围有限、为正且一致 |
| v1.0 默认检查点迁移 | 串行与 MPI 恢复通过 |
| 非默认黏度比或错误热力参数恢复 | 一致拒绝，无伪正常结果 |

## 4. 空间与时间收敛

空间卡口在 24/48/96 三层均匀网格对正弦温度扩散算子计算截断误差，同时保存内部、等温壁
带和全局的 L1/L2/Linf。完整数值见
[`stage-s-spatial-convergence.csv`](stage-s-spatial-convergence.csv)，摘要见
[`stage-s-validation-summary.json`](stage-s-validation-summary.json)。

| profile | 维数 | 内部 L2 阶 | 壁带 L2 阶 | 全局 L2 阶 | 冻结下限（内部/壁带/全局） |
|---|---:|---|---|---|---|
| `phenglei_wcns` | 2D | 3.9589, 4.0527 | 0.6077, 0.9043 | 1.1194, 1.4046 | 3.5 / 0.55 / 1.0 |
| `phenglei_wcns` | 3D | 3.9589, 4.0527 | 0.6077, 0.9043 | 1.1194, 1.4046 | 3.5 / 0.55 / 1.0 |
| `scmm6_wcns` | 2D | 5.6918, 6.0290 | 3.3778, 3.1057 | 3.8780, 3.6057 | 5.0 / 3.0 / 3.0 |
| `scmm6_wcns` | 3D | 5.6918, 6.0397 | 3.3778, 3.1057 | 3.8780, 3.6057 | 5.0 / 3.0 / 3.0 |

`phenglei_wcns` 的内部算子达到约四阶，但现有等温壁闭合使壁带和全局阶分别只有约
0.61--0.90 与 1.12--1.40。这是明确保留的算法限制，不被内部阶掩盖，也没有在 S 中改写
profile 或事后放宽阈值；是否在后续版本提高边界闭合阶留给 T→U 人工判断。

SSPRK3 在 48x48 周期等熵涡上使用 CFL `0.4/0.2/0.1/0.05`。相邻末态密度 L2 差依次为
`3.2962925571e-7`、`4.0852262328e-8`、`4.6211727582e-9`，观测阶为 `3.01236`、
`3.14409`，均高于冻结下限 2.8。

## 5. 经验 CFL 分类

以下仅是固定网格和终止时刻上的经验结果，不是稳定性定理：

| 问题 | CFL | 结果与诊断 |
|---|---|---|
| 周期无粘光滑流 | 0.1, 0.4, 0.8, 1.6, 3.2 | 全部成功 |
| Sutherland 黏性 smoke | 0.05, 0.2, 0.8, 2.0 | 全部成功 |
| Sod，关闭事后稳健化 | 2 | 成功 |
| Sod，关闭事后稳健化 | 4, 8 | `numerical_failure` |
| Sod，启用事后稳健化 | 2, 4, 8, 16 | 全部成功；最大 troubled cell 为 0/32/120/72，最大整步重试为 0/1/2/2 |
| Case07 曲线四块圆柱 | 0.1, 0.4, 1.6, 4.0 | 全部成功 |

扫描明确覆盖直接失败、局部面重算和整步缩步边界；运行工具不会把一次成功外推为普适保证。

## 6. 正式自动验收

| 卡口 | 结果 | 证据 |
|---|---:|---|
| 算法规格核验 | 6/6 通过 | `python tools/verify_algorithm_spec.py` |
| 输运生产运行 | 串行、2-rank 通过 | Sutherland 温度/黏度范围及最终场检查 |
| 空间收敛 | 8/8 profile/维数/区带组合满足阈值 | 三层网格完整 L1/L2/Linf 原始表 |
| SSPRK3 时间收敛 | 通过 | 两个观测阶均 `>=2.8` |
| CFL 分类 | 20/20 结果符合预期分类 | 成功、数值失败、troubled cell、局部重算及缩步诊断 |
| Release 串行 CTest | 59/59 通过，58.96 s | 含 v1.0 默认检查点迁移和非默认输运拒绝 |
| Release MPI CTest | 107/107 通过，104.62 s | 含 1/2/4/8-rank、输运及检查点路径 |
| 差异审计 | 通过 | `git diff --check` 无错误；仅保留明确排除的用户未跟踪目录 |

正式命令为：

```text
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DWCNS_ENABLE_MPI=OFF
cmake --build build -j 4
ctest --test-dir build --output-on-failure -j 4
cmake -S . -B build-mpi -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DWCNS_ENABLE_MPI=ON
cmake --build build-mpi -j 4
ctest --test-dir build-mpi --output-on-failure -j 4
python tools/verify_algorithm_spec.py
python tools/run_stage_s_time_convergence.py
python tools/run_stage_s_cfl_scan.py
git diff --check
```

## 7. 连续推进结论

S 的全部冻结自动卡口通过，没有把未执行的人工审查写成批准。空间/时间 log-log、PH 壁闭合
限制、CFL 推荐值和旧检查点策略等人工判断项保留到连续 R→S→T 完成后的 T→U 人工卡口。
依照现有授权，本候选直接合入 `release/v1.1.0`，随后从该合并提交建立
`stage/v1.1.0-t`。
