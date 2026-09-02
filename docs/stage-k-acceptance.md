# 阶段 K 正式验收报告（候选 v1）

## 1. 验收结论

阶段 K 候选 v1 的正式自动验收通过，等待项目负责人人工检查。候选实现了两套算法 profile 对应的层流 Navier--Stokes 粘性 WCNS 路径，包括公共物理状态 ghost、守恒型 primitive 梯度、两级梯度通信、常黏度/Sutherland 输运、Stokes 应力、Fourier 热流、粘性壁面强约束、粘性面通量通信、`1/Re` 单次缩放、显式粘性时间步和 SSPRK3 子步集成。阶段 L 未开始。

- 被测实现提交：`159aa3b7a43e034fd1b14cae5c428db3aaf17f28`
- 阶段分支：`stage/k-viscous-wcns`
- 候选标签：`stage-k-candidate-v1`（指向包含本报告的提交）
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
git status --porcelain
cmake --build build --config Release --clean-first --parallel
cmake --build build-mpi --config Release --clean-first --parallel
python tools/verify_algorithm_spec.py
ctest --test-dir build -C Release --output-on-failure
ctest --test-dir build-mpi -C Release --output-on-failure
git diff --check
git status --porcelain
```

结果：

- 正式验收开始时工作树：干净；
- 串行与 MPI clean-first Release 构建：通过；
- 阶段 G 算法规格检查：6/6 通过；
- 串行 CTest：9/9 通过；
- MPI CTest：15/15 通过，其中 `wcns.stage_k_solver.2` 实际使用双 rank；
- `git diff --check`：通过；正式测试结束后工作树仍干净。

## 4. 阶段 K 专项结果

| 专项 | 验证内容 | 结果 |
|---|---|---|
| 输运与本构 | 常黏度、Sutherland `T=1` 恢复参考比值、固定 Pr、Stokes 应力、Fourier 能量通量以及二维 z 严格约束 | 通过 |
| 二维/三维梯度 | 两套 profile 在仿射网格上恢复二维和三维线性 `u,v,w,T`；物理边/角梯度存储保持 NaN | 通过 |
| 梯度 operand 通信 | `q*S` 独立于状态 halo 交换；普通连接与旋转周期的面方向、速度张量和温度矢量变换 | 通过 |
| 梯度通信 | PH 两层、SCMM6 三层连接 cell halo；同 rank 与双 rank 多块路径使用同一连接计划和版本 | 通过 |
| 壁面迹 | PH 0--4 次、SCMM6 0--6 次 Dirichlet 法向导数矩；无滑移、等温和绝热真实面强约束 | 通过 |
| 粘性面通量 | 真实面 Cartesian 本构只乘一次面矢量；PH 层 0--1、SCMM6 层 0--2 的连接通量及旋转/方向变换 | 通过 |
| Reynolds 缩放 | 相同通量下 `Re` 从 10 改为 20，粘性残差严格减半；质量粘性残差严格为零 | 通过 |
| 时间步 | 四项 `(profile,dimension,SSPRK3)` `C_v` 独立配置、严格正值检查、局部无粘/粘性分母及 MPI 全局最小值 | 通过 |
| 求解闭环 | 每个 RK 子步执行状态/物理 ghost、无粘通量、operand、梯度、粘性通量、两类散度和源项；阶段 J 路径保持独立回归 | 通过 |
| 多块物理夹具 | 两套 profile 的双块自由流、Couette 和线性导热在单 rank 与双 rank 下通过冻结容差 | 通过 |
| 粘性制造解 | 正弦温度解由 `24^2` 加密到 `48^2`，PH 误差下降大于 12，SCMM6 大于 40 | 通过 |

Couette 夹具验证质量/动量平衡及解析粘性耗散 `(du/dy)^2/Re`；线性导热夹具验证能量方程中的常热流零散度。导热状态是常压、变密度静止接触面，Rusanov 的无粘接触耗散不计入导热误差，判据已在正式验收前冻结于 `docs/stage-k-design.md`。

## 5. 公共卡口状态

| 卡口 | 状态 | 证据 |
|---|---|---|
| G0 范围冻结 | 通过 | `docs/stage-k-design.md` |
| G1 接口与合法域 | 通过 | profile/version 检查、连接专用 halo、NaN 非法物理 ghost 二级量 |
| G2 构建与静态检查 | 通过 | 两套 clean-first Release 构建及 `git diff --check` |
| G3 单元与公式 | 通过 | 规格 6/6、二维/三维梯度、本构、壁面、缩放及稳定系数测试 |
| G4 串行回归 | 通过 | CTest 9/9 |
| G5 MPI 回归 | 通过 | CTest 15/15，双 rank 阶段 K 求解测试通过 |
| G6 数值验收 | 通过 | 自由流、Couette、线性导热及粘性制造解加密测试 |
| G7 文档与可追溯性 | 通过 | 设计、冻结门槛、本报告、阶段分支和候选标签同步 GitHub |
| G8 人工放行 | 待项目负责人 | 未合并 `main`，未进入阶段 L |

## 6. 已知边界与后续事项

1. 四项 `C_v` 的阶段 K 默认值均为 4.0，是保守起步配置；生产算例仍须针对具体网格、物性和 CFL 验证稳定裕量，不能把默认值解释为所有问题的最优上限。
2. 当前壁面数据在一个 patch 上是常值，因此壁面速度和等温壁温的切向导数为零；空间变化壁面数据需要扩展 `BoundaryData` 后另行设计和验收。
3. 旋转周期的梯度 operand、梯度和粘性通量变换已做纯函数测试；包含完整状态 halo、残差和推进的旋转周期系统等价性仍按路线图留在阶段 M。
4. 生产算例配置、初场、运行控制、监测、残差历史、CGNS 流场和重启输出均属于阶段 M；本阶段使用确定性测试夹具，不把测试代码伪装成用户算例驱动。
5. HLLC、特征变量重构及其他获批高阶格式属于阶段 L；阶段 K 保留阶段 J 的 Rusanov 无粘路径。

## 7. GitHub 同步与人工验收请求

阶段分支和 `stage-k-candidate-v1` 标签同步至 `https://github.com/wylx-2/wcns_v3.git` 后，本候选提交人工验收。自动验收通过不授权合并到 `main`，也不授权开始阶段 L；必须等待项目负责人明确同意。
