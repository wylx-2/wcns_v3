# 阶段 D：MPI 块分配与 halo 通信

日期：2026-08-25

状态：自动验收通过，连续进入阶段 E。

## 目标

在阶段 C 的稳定 `BlockId`、连接 Transform 和 halo 单元配对表之上建立 MPI
执行层。同一套交换计划必须同时支持同 rank 块间直接复制和跨 rank 非阻塞通信，
数值模块不能感知 CGNS 文件索引。

## 子任务

1. 增加可选 `WCNS_ENABLE_MPI` 构建，封装 MPI 初始化、结束、rank/size 和错误。
2. 根据块单元数执行确定性贪心负载分配，保持 `BlockId` 与 rank 独立。
3. 由全局连接生成本地接收计划、远端发送计划、稳定通信 ID 和消息标签。
4. 为多分量 cell-centered 字段实现本地复制及 `MPI_Irecv/MPI_Isend/Waitall`。
5. 在 1 rank、2 rank 下验证轴交换和反向连接的三层 halo 数据完全一致。

## 约束

- MPI 关闭时核心、网格、I/O 和数值库仍可构建测试；
- 不在 MPI 消息中发送裸 C++ 对象，只发送连续 `Real` 缓冲区；
- 所有 count/tag/rank 在调用 MPI 前验证；
- 每轮交换按连接和字段组件确定大小，不依赖进程到达顺序；
- MPI 生命周期不重复初始化，也不结束外部已初始化的 MPI。

## 验收

- clean-first 串行构建与测试通过；
- MPI clean-first 构建通过；
- `mpiexec -n 1` 与 `mpiexec -n 2` 回归通过；
- 本地与远端交换使用同一配对顺序，结果逐值相同；
- 异常分配、过大消息或不一致缓冲区得到明确错误。
