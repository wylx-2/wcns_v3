# 阶段 C 候选验收报告

日期：2026-08-24

状态：等待项目负责人检测。

## 已完成范围

- 参考 PHengLEI 的 `StructBC`、`InterfaceInfo` 和 CGNS 结构网格转换路径，明确
  全局块编号、邻块编号、owner rank、两侧范围和方向变换相互独立。
- `IndexTransform` 支持二维/三维有符号轴置换校验、点映射和求逆。
- `CgnsReader` 读取 `GridConnectivity1to1_t` 的连接名称、donor Zone、两侧
  `PointRange` 和 Transform，并统一转换为内部 0-based 索引。
- receiver/donor 节点与单元范围均验证 extent、面位置、端点映射和维数一致性。
- `StructuredMesh` 建立稳定的 `BlockId -> 本地容器下标` 映射，阻止重复块编号。
- 全局拓扑验证连接互易性、donor rank、ghost 宽度、两侧物理坐标共点以及
  方向相反但等价的有向范围。
- `make_halo_exchange_plan` 生成 receiver ghost 单元与 donor 内部单元的确定性
  配对表；后续本地复制与 MPI 通信不再解释 CGNS 数据。
- 物理边界保存在 `boundaries`，块连接保存在 `connectivities`，两者没有混用。

## 测试网格与覆盖

- 原有二维、三维单块 CGNS 及边界越界负例；
- 二维两块网格：receiver I 面映射 donor J 面，Transform 为 `{2, -1, 3}`；
- 三维两块网格：I 面对接且 J 切向反转，Transform 为 `{1, -2, 3}`；
- 单边连接、未知 donor、缩小 extent 后的连接范围越界；
- 重复 Transform 轴、donor rank/ghost 宽度不一致、接口坐标不共点；
- donor 法向厚度小于 ghost 宽度时拒绝生成 halo 计划。

二维和三维有效多块文件均通过 CGNS 官方 `cgnscheck`，无错误。输出中的
family/dataclass 缺省及反向有向范围提示属于测试文件的非强制警告。

## 综合验收命令

```powershell
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build --clean-first -j 4
ctest --test-dir build --output-on-failure --verbose
```

综合验收结果：clean-first 全量构建成功，3/3 tests passed。

## 已记录的失败验收

- `623eb26`：把合法二维轴交换 `{2, 1, 3}` 误写成非法 Transform 测试用例；
- `76e129a`：生成器验证函数仍把 Zone 数量写死为 1，多块 fixture 被误判；
- `4a2b63e`：尝试用公开 CGNS API 写越界 1-to-1 范围，库在写入时即拒绝。

对应修复均保留后续独立 Git 检查点并通过回归。

## 阶段边界

本阶段没有实现：

- MPI 初始化、块到 rank 的自动分配；
- 非阻塞发送/接收、消息标签和通信缓冲区；
- halo 数值打包、解包及本地块直接复制；
- 周期连接的旋转/平移物理量变换；
- 非共形、滑移、混合面和 Overset 连接；
- WCNS 离散格式、边界数值处理和时间推进。

上述内容必须在阶段 C 经项目负责人检测同意后再进入下一阶段。
