# 阶段 I 正式验收报告（候选 v2）

## 1. 验收结论

阶段 I 候选 v2 的正式自动验收通过；项目负责人于 2026-09-01 批准完成修订后进入阶段 J，人工卡口 G8 已通过。候选仅包含静态几何 profile、几何连接计划、共享面同步和守恒权重；没有接入阶段 J 的高阶流场推进、面通量或实际源项。

- 被测实现提交：`3e133ec97490b090cfcdc4671fadbcc6e2a0b508`
- 阶段分支：`stage/i-geometry-profiles`
- 候选标签：`stage-i-candidate-v2`（指向包含本报告的提交）
- 验收日期：2026-09-01（Asia/Shanghai）

## 2. 环境

- 操作系统：Windows
- CMake：3.28.0
- C/C++ 编译器：MinGW-w64 GCC 8.1.0
- 构建类型：串行和 MPI 均为 `Release`
- MPI：Intel MPI Library 2021.10，Build 20230619
- CGNS：仓库内固定的 CGNS 4.4.0 源码包，静态链接

## 3. 正式验收命令与结果

```powershell
cmake --build build --config Release --parallel
cmake --build build-mpi --config Release --parallel
python tools/verify_algorithm_spec.py
ctest --test-dir build -C Release --output-on-failure
ctest --test-dir build-mpi -C Release --output-on-failure
git diff --check
git status --short
```

结果：

- 串行 Release 构建：通过；
- MPI Release 构建：通过；
- 阶段 G 精确规格检查：6/6 通过；
- 串行 CTest：7/7 通过；
- MPI CTest：11/11 通过，其中新增的双 rank 几何连接计划验收通过；
- `git diff --check`：通过；
- 开始正式验收时工作树：干净。

## 4. 阶段 I 专项结果

| 专项 | 验证内容 | 结果 |
|---|---|---|
| I1 profile 绑定 | 两套 profile 工厂、名称/重启签名、部件同源校验和交叉组合拒绝 | 通过 |
| I2 PH 度量 | `GridDelta`、加密交错对称度量、二维/三维恒等及一般仿射映射 | 通过 |
| I3 SCMM6 度量 | `I6/D6-D4`、共同中心 `delta=D6*I6`、二维/三维解析映射 | 通过 |
| I4 连接几何 | 分阶段消息、层宽、唯一标签、donor 路径、轴置换/反向、旋转周期和共享面所有者 | 通过 |
| I5 守恒权重 | PH `N=4..12`、SCMM6 `N=5..12` 正权重及共享面加权抵消 | 通过 |
| I6 冻结输出 | `MetricField` 仅允许受控构造，发布后只读；参考有限体积诊断独立保存 | 通过 |
| 二维 SCMM z 分量复核 | `Array3D<Real>` 原本已值初始化为零；现进一步显式清零并逐面断言严格为零 | 通过 |
| 测试说明 | 所有 `tests/*.cpp` 验收函数定义前均有简短中文目的注释 | 通过 |

数值专项阈值及实测判定：

- 光滑扭曲网格从 `N=16` 加密到 `N=32`：PH Jacobian L2 误差比大于 3，SCMM6 大于 20；
- 离散度量恒等式最大残差：PH 小于 `2e-12`，SCMM6 小于 `2e-11`；
- 一维守恒权重最大残差和多块共享面不匹配量均小于 `1e-11`；
- 仿射网格 Jacobian 与面积矢量按各专项舍入误差容差通过。

## 5. 公共卡口状态

| 卡口 | 状态 | 证据 |
|---|---|---|
| G0 范围冻结 | 通过 | `docs/stage-i-design.md` |
| G1 接口与合法域 | 通过 | profile 同源校验、真实几何域和受控构造接口 |
| G2 构建与静态检查 | 通过 | 两套 Release 构建和 `git diff --check` |
| G3 单元与公式 | 通过 | 规格检查 6/6，PH/SCMM6 解析与非法路径测试 |
| G4 串行回归 | 通过 | CTest 7/7 |
| G5 MPI 回归 | 通过 | CTest 11/11，双 rank 几何计划、halo、Euler 和 runtime 通过 |
| G6 数值验收 | 通过 | 仿射精确性、扭曲网格收敛、度量恒等式和守恒权重 |
| G7 文档与可追溯性 | 通过 | 本报告、阶段提交、远程分支和候选标签共同保存 |
| G8 人工放行 | 通过 | 项目负责人于 2026-09-01 明确批准完成本次修订后进入阶段 J |

## 6. 已知边界与后续事项

1. `GeometryHaloPlan` 固定几何初始化所需的分阶段描述符、donor 路径和 MPI 消息身份；阶段 I 的 MPI 验收验证各 rank 生成一致计划，不包含阶段 J 的流场面通量传输。
2. `SharedMetricSynchronizer` 当前执行同一进程内持有的多块共享面发布；分布式运行使用同一所有者、变换和消息描述语义，后续随高阶初始化驱动接入实际缓冲区交换。
3. 阶段 E 的低阶 Euler residual 保持兼容；阶段 I 尚未让求解时间推进使用新 `MetricField`。
4. 物理状态 ghost、强边界面状态开关、Riemann 面通量和显式源项属于阶段 J。

## 7. 人工验收记录

项目负责人已批准完成二维 SCMM z 分量复核和修订后进入阶段 J。本报告记录的修订已完成且全量回归通过，可将候选 v2 合并到 `main`、创建 `stage-i-accepted` 标签并建立阶段 J 分支。
