# 阶段 B 候选验收报告

日期：2026-08-19（2026-08-24 补充范围校验复验）

状态：等待项目负责人检测。

## 已完成范围

- 从 PHengLEI 随附归档接入 CGNS 4.4.0，使用固定 SHA-256 校验。
- 采用 MinGW 可构建的静态 ADF 后端，不依赖网络、HDF5 或 Fortran。
- 提供二维、三维单块结构 CGNS 测试网格生成器。
- 读取所有 Base 和 Structured Zone 元数据，并分配稳定的全局 `BlockId`。
- 读取 `CoordinateX/Y/Z` 节点坐标；二维缺省 Z 坐标置零。
- 读取 `Vertex + PointRange` 物理边界，支持切向范围正向或反向排列。
- 将 CGNS 1-based 节点范围转换为内部 0-based 节点范围和相邻单元范围。
- 转换后逐轴验证节点范围位于 `vertex_extent` 内，并对推导后的单元面范围按
  `cell_extent` 再做一次防御性校验。
- 计算二维单元面积、三维单元体积、单元中心、方向面面积和单位法向。
- 拒绝退化网格面和非正面积/体积单元。

## 综合验收

```powershell
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build --clean-first -j 4
ctest --test-dir build --output-on-failure --verbose
```

结果：全量构建成功，3/3 tests passed。

测试包括：

- `wcns.generate_test_cgns`：生成并重读二维、三维及边界越界 CGNS 文件；
- `wcns.cgns_reader`：验证元数据、坐标、边界索引转换、越界拒绝和解析几何量；
- `wcns.unit`：验证基础数组、字段、结构块、CGNS 链接和非法几何拒绝。

生成文件还通过 CGNS 官方 `cgnscheck` 检查，二维和三维文件返回值均为 0。测试文件有未设置 family/dataclass 的推荐性警告；用于有向范围覆盖的反向 PointRange 也有方向提示，均不是合法性错误。

补充负例 `generated_invalid_2d.cgns` 的 `I-min` 边界将 J 终点写为 5，而网格
J 顶点范围为 1..4。CGNS 官方 `cgnscheck` 报告 `1 elements are out of range`；
读取器稳定拒绝该文件，并报告内部索引 4 超出 `[0, 3]`。2026-08-24 已执行
clean-first 综合构建，回归测试结果仍为 3/3 tests passed。

## 已记录的失败验收

- `c241c82`：CGNS 解包入口多指定一层，CMake 目标未创建。
- `9723375`：CGNS Windows 64 位类型与 MinGW C++ ABI 配置不匹配。
- `7075148`：生成器未继承 WCNS C++ 编译特性。
- `8f933c2`：CGNS 4.4.0 Windows 头使用 MSVC `__int64`，MinGW 需要兼容定义。

每个问题均保留独立 Git 检查点，后续提交已完成修复并通过验收。

## 当前约定

- CGNS 文件索引只存在于 I/O 层，内部一律使用 0-based 索引。
- CGNS Zone 映射为一个 `StructuredBlock`。
- 二维节点维度 `(Ni, Nj)` 映射为内部 `(Ni, Nj, 1)`。
- `cell_metrics.jacobian` 当前表示单位计算单元对应的物理面积或体积，与 `volume` 相同。
- 面法向沿正计算坐标方向，而不是统一采用块外法向。

## 阶段边界

本阶段尚未实现：

- `GridConnectivity1to1_t` 多块连接；
- CGNS Transform 轴交换和反向映射；
- 全局 Mesh 与本地块集合；
- MPI 块分配和 halo 通信；
- 非共形连接；
- HDF5 后端；
- WCNS 空间离散与时间推进。

这些内容只有在阶段 B 经项目负责人检测同意后才进入下一阶段。

## GitHub 同步状态

远程仓库为 `https://github.com/wylx-2/wcns_v3.git`。阶段 B 候选提交及标签已同步；
本次范围校验修复在验收提交后继续同步 `main`，并保留新的候选标签，不移动旧标签。
