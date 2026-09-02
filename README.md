# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

项目的数学与算法约定见 [`算法补充.md`](算法补充.md)，阶段 H--M 的开发流程、Git 规则和人工检查卡口见 [`docs/development-roadmap.md`](docs/development-roadmap.md)。

阶段 A--J 已按各自卡口完成，其中 H--J 已通过项目负责人人工验收。阶段 K 的粘性 WCNS 实现和专项数值测试已完成，正在执行正式自动验收；阶段 L 尚未获准开始。
当前程序具备 CGNS 结构多块网格读取、两套高阶几何 profile、确定性块分配、同 rank/MPI 非阻塞 halo 交换、WCNS-Euler 空间离散、层流 Navier--Stokes 粘性通量、显式源项和 SSPRK3 时间推进闭环。生产算例配置、初始化、监测以及 CGNS 流场/重启输出仍属于阶段 M。

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
- Euler 状态采用五分量守恒量，WCNS 首版逐原始变量重构并使用 Rusanov 通量。
- WCNS 求解需要至少三层 cell-centered ghost，块连接通信守恒量，接收后转换原始量。
