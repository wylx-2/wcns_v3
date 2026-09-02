# 阶段 M 设计：生产配置、结构分区与初始化

状态：M0--M5 已实现，正式自动验收通过；详见
[`stage-m-acceptance.md`](stage-m-acceptance.md)。

本阶段只建立生产入口到“可受控推进短步”的闭环；停止、监测、输出和重启在阶段 N
实现。通用约束见 [`release-development-plan.md`](release-development-plan.md)。

## 1. 配置格式

首版使用 UTF-8 严格键值格式，每行是 `key = value`，支持空行和以 `#` 开头的整行
注释。键区分大小写，字符串无需引号且会裁剪首尾空白；布尔值仅接受 `true|false`，
整数和浮点数必须完整消费。重复键、未知键、缺失必填键和非有限数确定性失败。

发布 schema 版本为 1。M 阶段冻结以下键：

- `schema_version`、`case.name`、`mesh.path`；
- `algorithm.profile`、`algorithm.reconstruction`、
  `algorithm.reconstruction_variables`、`algorithm.riemann`；
- `gas.gamma` 和 `gas.molar_mass` 或 `gas.specific_gas_constant` 二选一；
- 五个 `reference.*` 基准量；禁止输入 `Re` 或 `Ma`；
- `partition.mode`、`partition.allow_idle_ranks`、
  `partition.max_load_ratio`、`partition.min_cells_per_active_direction`；
- `initial.type` 与对应初场参数；
- `run.viscous`、`run.cfl`、`run.max_steps` 和 M 阶段短步参数；
- `boundary.default` 及 `boundary.<CGNS-name>.type`；
- `source.enabled` 及已冻结源项字段。

阶段 N 在同一 parser/registry 上增加 `run.mode`、停止、输出和监测键，不建立第二套配置
读取器。根 rank 读取原始文本，广播后各 rank 独立严格解析并比较 FNV-1a 64 位摘要。

## 2. 分区数据结构

`StructuredPartitionPlan` 保存原 zone 元数据和 `PartitionLeaf`。叶块使用原 zone
零基 cell 半开区间，局部坐标从零开始；由区间唯一导出 CGNS 一基顶点 hyperslab。
全局叶块 ID 按原 zone ID、cell begin 的字典序重新编号，不依赖切分发现顺序或 rank。

分区器不修改 `BlockDistribution`；它先产生叶块，后者再按叶块单元数分配 rank。兄弟
连接与原连接切片必须在完整计划层生成，之后才构造 `StructuredMesh`。

首版切分评分按以下顺序决定：

1. 切分后两侧均满足活动方向最小单元数；
2. 最小化两子块与目标单元数的最大偏差；
3. 最小化新增切面单元数；
4. 优先较长方向，最后按 I/J/K 和较小 cut 打破平局。

`auto_split` 至少产生 `min(rank_count, maximum_feasible_leaves)` 个叶块，并在负载比超过
阈值时继续细化。若请求 rank 数不可行且禁止空闲，返回包含最大可行叶块数的错误。

## 3. CGNS 与拓扑切片

`CgnsReader::read_block_range` 使用 `cg_coord_read` 的范围参数直接读取叶块顶点。
物理边界与叶块外表面取交并转换为局部范围。原一对一连接先在原 zone 坐标中取交，再
通过既有 `SignedAxisTransform` 映射到 donor 叶块。兄弟切面建立单位变换的互逆连接。
所有生成结果继续通过 `StructuredMesh::validate_connectivities`。

M2 分两步实施：

1. 先完成与 CGNS 无关的矩形叶块规划、覆盖验证、兄弟连接和负载分配；
2. 再接入 CGNS hyperslab、物理边界和原连接切片。

## 4. 初场与生产启动链

`FlowInitializer` 只写真实单元守恒量，初场坐标使用 `MetricField::cell_coordinates`。
首版内建：

- `uniform`；
- `quadrant_riemann`；
- `sod_x`；
- `isentropic_vortex`；
- `couette`；
- `manufactured_periodic`。

所有初场先得到 `rho,u,v,w,T` 或 `rho,u,v,w,p`，再通过唯一 `GasModel` 路径生成守恒
量。初始化后由现有 solver 统一交换连接 halo、填充物理 ghost，不直接填几何或二级量。

正式可执行程序启动顺序固定为：

```text
MPI 初始化 -> 广播并解析配置 -> CGNS 元数据 -> 分区计划 -> 叶块分配/读取
-> 拓扑与度量 -> 初场/边界数据 -> solver -> M 阶段短步运行
```

## 5. 小任务与提交边界

1. M0：本文和配置示例；
2. M1a：严格键值 parser、schema 和摘要；
3. M1b：`CaseConfig` 到现有 gas/profile/solver 配置的转换；
4. M2a：矩形结构分区计划与可行性；
5. M2b：叶块拓扑切片、CGNS 子区间读取和 MPI 分配；
6. M3：初场 registry；
7. M4：`wcns_run` 正式入口；
8. M5：单 zone 1/2/4/8 rank 和多 zone 自动验收。

每项通过局部测试后独立提交。M5 正式验收通过后形成阶段 M 里程碑并直接进入 N。

## 6. 非目标

- 阶段 M 不实现最终 CGNS/Tecplot 输出或重启；
- 不实现低 Mach 预处理 L6；
- 不实现一般表达式解释器，解析/制造初场由显式注册模型提供；
- 不为填满 rank 而违反 profile 模板尺寸；
- 不以所有 rank 读取完整原 zone 作为 hyperslab 的临时替代。
