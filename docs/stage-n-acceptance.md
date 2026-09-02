# 阶段 N 正式自动验收报告（候选 v1）

## 1. 验收结论

阶段 N 的冻结范围 N0--N8 已完成，正式自动验收通过。生产入口已经由阶段 M 的受控短步
扩展为唯一的定常/非定常运行状态机，并具备残差判停、精确时间事件、可选择输出、原 zone
重组和可改变 rank 数的检查点重启。按负责人批准的连续发布流程，本里程碑不等待逐阶段人工
检查、不合并 `main`，创建不可移动候选标签后直接进入阶段 O。

- 被测实现提交：`17f47cf284dc`；
- 发布开发分支：`stage/release`；
- 候选标签：`stage-n-candidate-v1`（指向包含本报告的提交）；
- 验收日期：2026-09-02（Asia/Shanghai）。

## 2. 完成内容

| 小任务 | 完成内容 | 主要提交 |
|---|---|---|
| N0 | 冻结运行状态机、停止优先级、输出/重启格式和小任务边界 | `fb77a85` |
| N1 | MPI 全局加权五分量残差、冻结参考值和统一停止原因 | `9df0ece` |
| N2 | 无粘/黏性适配器、唯一 `SimulationDriver`、精确事件时间步裁剪 | `f3ebba7` |
| N3 | 通用调度、历史、manifest、目录与原子提交规则 | `fffd9a6` |
| N4 | 可选择全场量/统计量 registry、量纲缩放和依赖校验 | `fefb37e` |
| N5 | 原 zone MPI 重组、CGNS/Tecplot 流场及 TXT/Tecplot 序列 | `87295cb` |
| N6 | CGNS 检查点、网格/数值签名和同/异 rank 重启 | `24803b3` |
| N7 | 停止码、I/O 冲突、非法重启和非平凡重启连续性测试 | `a58f745` |
| N8 | 生产 registry 注入、运行/配置/输出/重启和诊断文档 | `5c9f731`, `17f47cf` |

`SimulationDriver` 在所有 rank 使用同一个完整步状态，只在完整 SSPRK3 步后判定停止和触发
输出。非定常步会裁剪到 `t_end` 和最近输出时刻；定常残差只覆盖真实单元，以 profile 守恒
权重和正 Jacobian 作全局归约。数值失败不写貌似正常的最终流场或检查点。

## 3. 正式验收环境与命令

- Windows；CMake 3.28.0；
- MinGW-w64 GCC 8.1.0；
- 串行及 Intel MPI 2021.10 Release 构建；
- 仓库内固定 CGNS 4.4.0 静态 ADF 后端。

```powershell
cmake -S . -B build -DWCNS_ENABLE_MPI=OFF -DWCNS_ENABLE_CGNS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
cmake -S . -B build-mpi -DWCNS_ENABLE_MPI=ON -DWCNS_ENABLE_CGNS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-mpi --parallel 4
python tools/verify_algorithm_spec.py
ctest --test-dir build -C Release --output-on-failure
ctest --test-dir build-mpi -C Release --output-on-failure
git diff --check
git status --short
```

结果：两套 Release 构建通过；算法规格 6/6 通过；串行 CTest 22/22 通过；MPI CTest 36/36
通过；正式验收开始和结束时工作树干净，`git diff --check` 通过。重新配置后 manifest 正确
记录 `program_version=0.1.0`、被测 Git 提交、编译器、构建类型、配置/网格/重启签名和输出
文件清单。

## 4. 阶段专项证据

| 专项 | 验证内容 | 结果 |
|---|---|---|
| 定常停止 | 五分量 `L2/Linf`、参考冻结、连续通过、正常完成优先于最大步数 | 通过 |
| 非定常停止 | 一步被精确裁剪到 `t_end=0.001`，不由残差提前结束 | 通过 |
| 安全/失败 | 最大步数退出码 2、墙钟/信号单元路径、数值失败原因、非法签名退出码 1 | 通过 |
| 调度 | 步/时间/显式/初末事件取并集并对相同 step/time 去重 | 通过 |
| I/O 冲突 | `allow_existing=false` 首次成功、第二次在既有目录一致失败 | 通过 |
| 流场输出 | CGNS 与 Tecplot 按原 zone 重组，独立重读 12 个字段及文本结构 | 通过 |
| 历史/统计 | TXT 字符串停止原因、Tecplot 数值停止码、全局守恒统计独立重读 | 通过 |
| 检查点 | 独立重读网格、五个守恒场、格式/网格/数值签名和残差状态 | 通过 |
| 重启连续性 | Sod 非平凡场串行连续算与 1→1、1→2 rank 续算逐点比较 | 通过 |
| 串并行 | 生产入口、输出和重启覆盖 1/2 rank；自由流继续覆盖 1/2/4/8 rank | 通过 |

Sod 连续计算从 `t=0` 到 `0.0006`；中途在 `t=0.0003` 写检查点，然后以新的 1 rank 和
2 rank 分区续算。五个守恒场独立 CGNS 重读后的最大绝对差为 `4.44089e-16`，小于冻结的
异 rank 容差 `1e-12`。

CGNS 4.4.0 `cgnscheck` 对生产流场和检查点均完成检查且无错误。它仍报告输入 zone 无 family、
坐标/字段未写 dimension exponents、Mach 数据类提示以及自定义 `Jacobian` 不是标准数据名等
8--15 条元数据警告；字段位置、数值和检查点重读不受影响，阶段 O 发布整理时补充或明确保留。

## 5. 卡口状态

| 卡口 | 状态 | 证据 |
|---|---|---|
| G0 范围冻结 | 通过 | `stage-n-design.md` 与发布计划 |
| G1 接口/合法域 | 通过 | 唯一状态机、严格停止原因、registry 和检查点签名 |
| G2 构建/静态 | 通过 | 两套 Release 构建、干净树和 diff 检查 |
| G3 单元/公式 | 通过 | 规格 6/6、停止/调度/量纲/重组/重启单元与系统测试 |
| G4 串行回归 | 通过 | CTest 22/22 |
| G5 MPI 回归 | 通过 | CTest 36/36，含 1→2 rank 非平凡续算 |
| G6 数值专项 | 通过 | 精确终止时间、自由流、Sod 连续性、输出独立重读 |
| G7 文档/追溯 | 通过 | 运行指南、分目标提交、报告和候选标签 |
| G8 人工放行 | 本阶段不适用 | 仅阶段 O 发布候选集中执行 |

## 6. 已知边界与阶段 O 输入

1. 当前内建全场量覆盖原始量、守恒量、声速、Mach、总焓、熵代理、黏度和 Jacobian；涡量、
   散度、Q、壁面力/热流等需要梯度或边界几何的量尚未注册，请求时明确失败。扩展后的字段和
   统计 registry 已可注入生产输出对象，不需要修改状态机。
2. 当前内建统计量为全局质量、三分量动量和能量。点探针、区域 min/max/mean/RMS/integral、
   壁面系数和守恒通量账目需要随阶段 O 的实际算例定义坐标/边界选择 schema 后实现，不能用
   未冻结的隐式选择规则。
3. 首版 I/O 由 rank 0 重组单文件，不提供并行 HDF5；检查点仍与 rank 数和叶块切分解耦。
4. 标准入口的 farfield/inflow 和等温壁数据取自通用初场参数，尚不覆盖移动壁、分边界非均匀
   入口和复杂壁温；Couette/导热等阶段 O 算例必须先补充明确的边界数据 schema。
5. 阶段 N 验收的是运行、停止、输出和重启机制，不替代阶段 O 的等熵涡收敛阶、Sod/四象限
   结构、黏性解析解、源项制造解、三维几何、周期连接和发布性能矩阵。

## 7. Git/GitHub 状态

阶段 N 各小目标均形成独立本地提交；包含本报告的提交将创建
`stage-n-candidate-v1`。候选建立后立即尝试向
`https://github.com/wylx-2/wcns_v3.git` 同步 `stage/release` 和候选标签；同步结果若失败，
将像阶段 M 一样使用后续独立提交记录，不移动候选标签，也不阻塞本地进入阶段 O。阶段 O
集中人工验收前必须补齐远程同步。

后续状态更新：候选建立后本次同步成功；远程已包含 `stage/release`、
`stage-n-candidate-v1`，并补齐先前缺失的阶段 F、L、M 里程碑标签。
