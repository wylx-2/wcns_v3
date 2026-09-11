# 阶段 M 正式自动验收报告（候选 v1）

## 1. 验收结论

阶段 M 的冻结范围 M0--M5 已完成，正式自动验收通过。按照发布开发流程，本里程碑不等待
逐阶段人工检查、不合并 `main`，验收提交创建不可移动标签后直接进入阶段 N。

- 被测实现提交：`7a36ee5`；
- 发布开发分支：`stage/release`；
- 候选标签：`stage-m-candidate-v1`（指向包含本报告的提交）；
- 验收日期：2026-09-02（Asia/Shanghai）。

## 2. 完成内容

| 小任务 | 完成内容 | 主要提交 |
|---|---|---|
| M0 | 冻结严格配置、叶块分区、CGNS hyperslab、初场和短步入口 | `ef19e2c` |
| M1 | 严格 UTF-8 键值 schema、配置摘要、模型/profile/源项转换 | `ecc414a` |
| M2a | 确定性受约束结构二分、可行性与负载诊断、叶块分配 | `e734c81` |
| M2b | owner-only CGNS 坐标子区间读取、边界/原连接切片、兄弟连接 | `b6332d7` |
| M3 | 六种生产初场 registry，按单元中心只初始化真实单元 | `15855c6` |
| M4 | `wcns_run` 配置广播、分区、度量、初场、边界和短步推进闭环 | `f76f59a` |
| M5 | 单 zone 1/2/4/8 rank 运行及原生多块/轴变换/MPI 回归 | `7a36ee5` |

正式入口在 rank 0 读取配置并广播原始文本，各 rank 独立严格解析并比较 64 位摘要；CGNS
元数据在各 rank 一致读取，但坐标数组只按 owner 的叶块 hyperslab 分配和读取。单 zone 少于
rank 时，32 x 16 测试网格分别形成足量合法叶块，并完成 1、2、4、8 rank 的相同短步路径。

## 3. 正式验收环境与命令

- Windows；CMake 3.28.0；
- MinGW-w64 GCC 8.1.0；
- 串行及 Intel MPI 2021.10 Release 构建；
- 仓库内固定 CGNS 4.4.0 静态库。

```powershell
git status --short
cmake --build build --config Release --parallel 4
cmake --build build-mpi --config Release --parallel 4
python tools/verify_algorithm_spec.py
ctest --test-dir build -C Release --output-on-failure
ctest --test-dir build-mpi -C Release --output-on-failure
git diff --check
```

结果：两套构建通过；算法规格 6/6 通过；串行 CTest 10/10 通过；MPI CTest 19/19
通过；正式验收开始时工作树干净，`git diff --check` 通过。

## 4. 阶段专项证据

| 专项 | 验证内容 | 结果 |
|---|---|---|
| 配置一致性 | 字符串广播、摘要全 rank 相等、未知/重复/缺失/非法值拒绝 | 通过 |
| 分区覆盖 | 原 zone 单元无重叠且恰好覆盖、受 profile 最小模板约束 | 通过 |
| rank 可用性 | 单 zone 32 x 16 单元在 1/2/4/8 rank 均成功推进且无空闲 rank | 通过 |
| CGNS 子区间 | owner-only 坐标读取、物理边界切片、非 owner 坐标保持 NaN | 通过 |
| 多块连接 | 原生双 zone、轴置换/反向连接切片、兄弟连接及消息计划唯一 | 通过 |
| 周期/变换 | 阶段 J/K 的旋转周期状态、通量和梯度变换回归继续通过 | 通过 |
| 生产闭环 | 配置到 CGNS、度量、初场、halo、边界、CFL 与 SSPRK3 短步 | 通过 |
| 串并行基线 | 自由来流在 1/2/4/8 rank 全部完成两步且回归无失败 | 通过 |

## 5. 卡口状态

| 卡口 | 状态 | 证据 |
|---|---|---|
| G0 范围冻结 | 通过 | `stage-m-design.md` 与发布计划 |
| G1 接口/合法域 | 通过 | 强类型叶块范围、严格配置和坐标可选拓扑验证 |
| G2 构建/静态 | 通过 | 两套 Release 构建、干净树和 diff 检查 |
| G3 单元/公式 | 通过 | 规格 6/6、配置/分区/初始化/CGNS 单元测试 |
| G4 串行回归 | 通过 | CTest 10/10 |
| G5 MPI 回归 | 通过 | CTest 19/19，含生产入口 2/4/8 rank |
| G6 数值专项 | 通过 | 冻结的阶段 M 自由流短步与串并行基线 |
| G7 文档/追溯 | 通过 | 设计、报告、分目标提交和候选标签 |
| G8 人工放行 | 本阶段不适用 | 仅阶段 O 集中执行 |

## 6. 已知边界与阶段 N 输入

1. 阶段 M 只提供 `max_steps` 控制的短步循环；steady/unsteady 状态机、残差判停、准确
   `t_end`、wall-time 和统一停止原因属于阶段 N。
2. 最终/周期流场、历史、统计量、manifest、CGNS/Tecplot 输出和检查点重启尚未实现，按
   阶段 N 冻结任务执行。
3. 当前生产边界数据足以覆盖自由来流、外推和固定温度静止壁 smoke run。Couette 移动壁、
   分边界非均匀目标和发布算例参数将在 N/O 配置扩展中显式提供，不得依赖测试夹具默认值。
4. 冻结的 M3 registry 覆盖 uniform、Sod、四象限、等熵涡、Couette 和制造周期初场；从
   CGNS FlowSolution 恢复由阶段 N 的检查点/重启唯一实现，不建立第二条初始化旁路。

## 7. Git/GitHub 状态

阶段 M 各小目标均已形成独立本地提交，并已创建 `stage-m-candidate-v1` 标签。2026-09-02
向 `https://github.com/wylx-2/wcns_v3.git` 推送 `stage/release` 及标签时连接被重置，故当前
状态为“远程未同步”。按流程这不阻塞本地进入阶段 N；阶段 O 集中人工验收前必须补齐远程
同步。

后续状态更新：2026-09-02 网络恢复后，`stage/release`、`stage-m-candidate-v1` 以及后续
阶段 N 候选均已成功推送到上述 GitHub 仓库，远程同步缺口已经关闭。
