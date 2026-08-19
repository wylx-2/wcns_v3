# 阶段 A 候选验收报告

日期：2026-08-19

状态：已由项目负责人于 2026-08-19 检测通过，允许进入阶段 B。

## 已完成范围

- CMake/C++20 工程骨架与 CTest 测试入口。
- `Index3`、`Extent3` 和支持反向范围的 `IndexRange3`。
- 带可配置 ghost 层、负索引和连续内存的 `Array3D<T>`。
- 以分量为最快变化维度的多分量 `Field<T>`。
- 二维/三维 `StructuredBlock`，包含节点坐标、单元几何量、方向面几何量、Euler 流场、边界 patch 和块连接元数据。
- 块全局编号、owner rank 和本地存储概念相互独立。

## 综合验收命令

```powershell
cmake --build build --clean-first
ctest --test-dir build --output-on-failure --verbose
```

结果：构建成功；`wcns.stage_a` 通过，1/1 tests passed。

## 已记录的失败验收

提交 `ffd1388` 记录了 MinGW GCC 8.1 不支持 C++20 默认比较运算符造成的构建失败。后续提交 `18040cd` 改为显式比较运算符并通过验收。保留该失败检查点是为了遵守“每次验收无论成功或失败均进行 Git 管理”的项目规则。

## 阶段边界

本阶段尚未实现：

- CGNS 文件读取；
- 全局 Mesh/本地块管理；
- MPI 初始化和块分配；
- 接口索引映射与 halo 通信；
- 几何量计算；
- WCNS 数值格式。

项目负责人已经批准阶段 A。后续内容按阶段闸门继续实施。

## 远程仓库状态

本地仓库位于 `wcns/.git`，分支为 `main`。截至本报告生成时尚未配置 GitHub remote，因此当前提交和候选标签仅保存在本地。
