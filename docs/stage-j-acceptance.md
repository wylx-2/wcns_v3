# 阶段 J 正式验收报告（候选 v1）

## 1. 验收结论

阶段 J 候选 v1 的正式自动验收通过，当前等待项目负责人人工检测。候选实现了两套几何 profile 对应的无粘 WCNS 空间离散、公共物理边界 ghost、Rusanov 界面通量、块间 `FaceFluxHalo`、显式源项和 SSPRK3 子步装配。阶段 K 尚未开始。

- 被测实现提交：`edf18e667c4554c19e0b6a48b92794c31038faef`
- 阶段分支：`stage/j-inviscid-wcns`
- 候选标签：`stage-j-candidate-v1`（指向包含本报告的提交）
- 验收日期：2026-09-02（Asia/Shanghai）

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

- 串行与 MPI Release 构建：通过；
- 阶段 G 算法规格检查：6/6 通过；
- 串行 CTest：8/8 通过；
- MPI CTest：13/13 通过，包含双 rank 阶段 J 多块求解测试；
- `git diff --check`：通过；
- 正式验收开始时工作树：干净。

## 4. 阶段 J 专项结果

| 专项 | 验证内容 | 结果 |
|---|---|---|
| 公共状态 ghost | 三个法向 ghost 层只填面状 slab；温度型 primitive 为权威值并同步派生压力和守恒量；边角保持 NaN | 通过 |
| 边界面状态 | 强/弱开关、固壁零法向质量通量、入口/远场线性化特征混合及超声速分支 | 通过 |
| 重构 | `linear5` 四次多项式精确性、scaled WCNS-JS 尺度不变性、有限性与正性逐级回退 | 通过 |
| Riemann | 仅提供严格单位法向的 Rusanov 接口，真实面积由调用方乘一次 | 通过 |
| 高阶散度 | PH `D4-PH` 与 SCMM6 `D6/D4` 各自独立，恒定流与光滑通量收敛测试通过 | 通过 |
| 面通量通信 | PH 层 0--1、SCMM6 层 0--2；所有者、方向、索引变换、消息版本及旋转周期动量通量变换 | 通过 |
| 多块求解 | 单进程双块和双 rank 双块的残差及一个 SSPRK3 步与单块自由流结果一致 | 通过 |
| 显式源项 | 常量守恒源、体力和制造源；子步时间；二维约束；MPI 全局体积加权平衡 | 通过 |
| 二维 SCMM 面度量 | 底层数组原本值初始化为零；现对 `s_i.z`、`s_j.z` 显式清零并逐面断言严格为零 | 通过 |

光滑正弦面通量从 `N=16` 加密到 `N=32` 时，PH 误差比要求大于 3，SCMM6 误差比要求大于 10，实测均满足门槛。两个 profile 的自由流残差、串行多块与双 rank 一致性均按测试规定容差通过。

## 5. 公共卡口状态

| 卡口 | 状态 | 证据 |
|---|---|---|
| G0 范围冻结 | 通过 | `docs/stage-j-design.md` |
| G1 接口与合法域 | 通过 | profile 同源校验、NaN 非法 ghost、单位法向及版本检查 |
| G2 构建与静态检查 | 通过 | 两套 Release 构建及 `git diff --check` |
| G3 单元与公式 | 通过 | 规格 6/6、重构/Riemann/边界/源项专项测试 |
| G4 串行回归 | 通过 | CTest 8/8 |
| G5 MPI 回归 | 通过 | CTest 13/13，双 rank 通量交换、求解及源项平衡通过 |
| G6 数值验收 | 通过 | 自由流、多块等价、光滑通量收敛和一个时间步状态 |
| G7 文档与可追溯性 | 通过 | 本报告、阶段分支和候选标签已同步至 GitHub |
| G8 人工放行 | 待项目负责人检测 | 未进入阶段 K |

## 6. 已知边界与后续事项

1. 亚声速入口/远场首版采用以内部状态冻结的五波线性化特征混合，不是完整非线性特征边界；后续基准需要继续验证反射水平。
2. 阶段 J 已验证旋转周期面通量的方向及 Cartesian 动量分量变换；旋转周期的完整状态 halo、残差和时间推进系统等价性按路线图留到阶段 M。
3. 当前唯一 Riemann 求解器为 Rusanov，HLLC 与特征变量重构属于阶段 L。
4. 黏性梯度、输运模型、应力和热流均未实现，属于阶段 K。
5. 阶段 E 的有限体积 Euler 路径继续作为独立基线；阶段 J 的实际源项通过 `InviscidWcnsSolver` 路径执行。

## 7. GitHub 同步状态

2026-09-02 已将 `main`、`stage/i-geometry-profiles`、`stage-i-candidate-v2`、`stage-i-accepted`、`stage/j-inviscid-wcns` 和 `stage-j-candidate-v1` 成功同步至 `https://github.com/wylx-2/wcns_v3.git`。阶段 J 未合并到 `main`。

## 8. 人工验收记录

候选 v1 等待项目负责人检测。只有人工批准后，才能将阶段分支合并到 `main`、创建阶段 J accepted 标签并进入阶段 K。
