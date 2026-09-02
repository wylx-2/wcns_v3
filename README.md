# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

项目的数学与算法约定见 [`算法补充.md`](算法补充.md)，阶段 H--O 的开发流程与 Git 规则见 [`docs/development-roadmap.md`](docs/development-roadmap.md)，阶段 M--O 的结构分区、停止/输出契约和发布算例矩阵见 [`docs/release-development-plan.md`](docs/release-development-plan.md)。

阶段 A--J 已按各自卡口完成，其中 H--J 已通过项目负责人人工验收。阶段 K、L 候选版本的正式自动验收已通过，K/L 尚未进行的人工检查已按负责人批示并入最终发布算例验收，不视为追认通过，也不再阻塞 M--O 连续开发。阶段 M 已通过正式自动验收；阶段 N 已完成实现，正在形成正式候选。

当前程序具备 CGNS 结构多块网格读取、两套独立高阶几何 profile、单 zone 受约束二次剖分、同 rank/MPI 非阻塞 halo 交换、四种界面重构、Rusanov/HLLC/Roe、WCNS-Euler 空间离散、层流 Navier--Stokes 粘性通量、显式源项和 SSPRK3 推进。正式入口支持严格配置、定常/非定常停止、MPI 全局残差、精确时间事件、CGNS/Tecplot 流场、TXT/Tecplot 历史与统计、manifest，以及可改变 rank 数和叶块划分的 CGNS 检查点重启。具体用法见 [`docs/runtime-guide.md`](docs/runtime-guide.md)；阶段 O 将用发布算例矩阵集中验收数值精度、守恒和全部生产组合。

## 构建与测试

```powershell
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
ctest --test-dir build --output-on-failure
```

启用 MPI（Windows/MinGW 下使用 Intel MPI）：

```powershell
cmake -S . -B build-mpi -G "MinGW Makefiles" -DWCNS_ENABLE_MPI=ON
cmake --build build-mpi
ctest --test-dir build-mpi --output-on-failure
```

CGNS 4.4.0 源码归档随仓库提供，CMake 会以静态 ADF 后端构建，不需要联网下载或单独安装 HDF5。

## 设计约定

- 内部索引从零开始。
- 物理单元范围为 `[0, n)`，ghost 索引允许为负数。
- 二维网格仍使用三维索引和五分量 Euler 状态。
- 全局块编号、MPI 所属进程和进程内数组下标相互独立。
- Euler 状态采用五分量守恒量；重构空间和 Riemann 求解器由严格配置选择。
- WCNS 求解需要至少三层 cell-centered ghost，块连接通信守恒量，接收后转换原始量。
