# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

阶段 A、B 已通过项目负责人验收，阶段 C 已通过自动验收并获准跳过人工门禁。
阶段 D 已通过自动验收，支持确定性块分配以及同 rank/MPI 非阻塞 halo 交换。
当前连续执行阶段 E，建立可运行的 WCNS Euler 求解骨架。

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
