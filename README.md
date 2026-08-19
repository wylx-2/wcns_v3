# WCNS

一个面向结构多块网格、CGNS 和 MPI 并行设计的小型高阶 CFD 程序。

阶段 A 已通过项目负责人验收。阶段 B 候选版本已经支持 CGNS 4.4.0 单块结构网格的元数据、节点坐标和物理边界读取，以及二维/三维基础几何量计算。多块连接和 MPI 将在阶段 B 获批后进入后续阶段。

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
