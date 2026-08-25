# 阶段 D 自动验收报告

日期：2026-08-25

状态：自动验收通过，按项目负责人授权连续进入阶段 E。

## 已完成范围

- 可选 `WCNS_ENABLE_MPI` 构建；关闭 MPI 时运行时退化为单 rank。
- Windows/MinGW 自动接入 Intel MPI，其他平台使用 CMake `FindMPI`。
- MPI 生命周期、线程级别、rank/size、屏障和全局 sum/min/max/all-true 规约。
- 按块单元数执行确定性最大优先贪心分配，并保存每个 rank 的计算负载。
- `LocalBlockSet` 只注册 owner rank 对应的块，自动更新连接 donor rank。
- 将互易连接配成稳定 `ConnectionId`，两个方向使用不同且确定的消息标签。
- `DistributedTopology` 分离本地复制、远端接收和远端发送计划。
- `BlockFieldRegistry` 对任意多分量 cell-centered `Field<Real>` 建立块字段视图。
- `HaloExchanger` 使用同一 `HaloCellPair` 顺序完成本地复制或
  `MPI_Irecv/MPI_Isend/MPI_Waitall` 打包交换。
- MPI count、消息标签、字段组件、owner rank、ghost 索引和 donor 内部索引均检查。

## 验收结果

串行 clean-first：5/5 tests passed。

MPI clean-first：7/7 tests passed，其中包括：

- 1 rank 与 2 rank MPI 生命周期及全局规约；
- 二维轴交换、切向反转、三层 ghost 的同 rank 直接复制；
- 相同算例在 2 rank 上的跨进程非阻塞交换；
- 每个 receiver ghost 值逐项等于预期 donor 内部值。

## 当前边界

- 只通信 cell-centered `Real` 字段；节点场和面场尚未接入。
- 当前每轮交换创建临时缓冲区，长期运行时可进一步缓存计划和缓冲区。
- 负载权重暂用单元数，尚未考虑不同物理模型或加速设备权重。
- CGNS 拓扑解析阶段仍会在各 rank 读取全局小网格以完成严格验证；大规模元数据
  广播和真正的并行 CGNS I/O 留待后续优化。

