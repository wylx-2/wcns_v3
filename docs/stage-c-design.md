# 阶段 C：共形多块连接设计

日期：2026-08-24

状态：执行中。

## 目标

读取结构 CGNS 的 `GridConnectivity1to1_t`，将文件中的 donor 名称、1-based
有向节点范围和 Transform 转换为可供后续本地复制及 MPI halo 通信共同使用的
0-based 多块拓扑。本阶段不发送 MPI 消息。

## PHengLEI 对照结论

参考 `Pre_CGNSConversion_Struct.cpp`、`Geo_StructBC.h` 和 `Geo_Interface.h`：

- 连接需要同时保存当前块与邻块编号，不能用 MPI rank 替代块编号；
- 当前面和邻块面都需要轴、侧别、起止范围和方向关系；
- CGNS Transform 是有符号轴置换，应直接保存并验证，而不是只保存“同向/反向”；
- 邻块按全局块编号组织，owner rank 是后续分配结果；
- 通信层应消费已经验证的拓扑，不应再次解释 CGNS 名称和 1-based 索引。

WCNS 保留上述分层思想，但采用值类型 `ConnectivityPatch`、`IndexRange3`、
`FaceLocation` 和 `IndexTransform`，避免 PHengLEI 中分散的裸数组及隐式下标语义。

## 子任务

1. 为 `IndexTransform` 增加维度感知的合法性、求逆和点映射。
2. 读取每个 Zone 的 `cg_n1to1/cg_1to1_read` 数据。
3. 按同一 Base 内唯一 Zone 名称解析 donor，转换并验证两侧节点/单元范围。
4. 建立 `StructuredMesh` 全局块容器，检查块编号、名称及互易连接。
5. 生成含轴交换、反向映射的二维/三维多块 CGNS 文件，并覆盖未知 donor、
   越界范围、非法 Transform 和缺失互易连接等异常。

## 数据流

```text
CGNS Zone/GridConnectivity1to1
        │ name, donor, ranges, transform
        ▼
CgnsReader（1-based → 0-based，范围与映射验证）
        ▼
StructuredBlock::connectivities
        │ receiver/donor BlockId + face/range/transform
        ▼
StructuredMesh（全局注册、互易性校验）
        ▼
后续阶段：块分配 → 本地复制或 MPI halo 通信
```

## 阶段 C 验收条件

- 二维和三维共形多块 CGNS 可完整读取；
- 轴交换与反向索引映射的端点和单元数量一致；
- 每条连接都能解析为稳定 `BlockId`，双向记录互相匹配；
- 物理边界和连接面保持分离；
- 异常 donor、范围、Transform 或单边连接会抛出明确错误；
- clean-first 构建成功，全部测试通过。

