# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

阶段 A、B 已通过项目负责人验收。阶段 C 候选版本已经支持 CGNS 共形结构多块
连接、轴交换/反向索引映射、全局互易拓扑校验和 halo 单元配对计划，目前等待检测。
MPI 块分配和跨进程通信将在阶段 C 获批后的阶段进入。

## 构建与测试

```powershell
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
ctest --test-dir build --output-on-failure
```

CGNS 4.4.0 源码归档随仓库提供，CMake 会以静态 ADF 后端构建，不需要联网下载或单独安装 HDF5。

## 设计约定

- 内部索引从零开始。
- 物理单元范围为 `[0, n)`，ghost 索引允许为负数。
- 二维网格仍使用三维索引和五分量 Euler 状态。
- 全局块编号、MPI 所属进程和进程内数组下标相互独立。
