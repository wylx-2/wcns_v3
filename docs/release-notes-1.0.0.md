# WCNS v1.0.0 发布说明

WCNS v1.0.0 是项目首个整理发布版本。它面向结构多块 CGNS 网格和 MPI 并行计算，提供统一的
`wcns_run` 生产入口、严格配置、可重启时间推进以及可审计输出。

## 主要能力

- 二维/三维结构多块 ADF-CGNS 读取，共形普通、平移周期和旋转周期连接；zone 数少于 MPI
  rank 数时可执行满足高阶模板下限的确定性二次剖分。
- 相互独立的 `phenglei_wcns` 与 `scmm6_wcns` 度量/差分 profile。
- `zero_order`、`linear5`、WENO-JS、WENO-Z、MDCD-LINEAR 和 MDCD-HYBRID 重构；Rusanov、
  HLLC 和 Roe Riemann 求解器。
- 热完全理想气体 Euler 与层流 Navier--Stokes 方程、SSPRK3、均匀守恒源、体力、压力梯度和
  制造源。
- 定常残差、非定常物理时间、最大步数、墙钟及终止信号停止；完整步 CGNS 检查点支持改变
  合法 MPI rank 数和叶块划分后续算。
- CGNS/Tecplot 流场、TXT/Tecplot history/statistics、manifest，以及槽道截面流量和壁面摩擦
  等内建统计。
- 经典双马赫反射专用初场和时变边界。

## 验证状态

自动发布矩阵覆盖自由流、等熵涡、间断、黏性流、源项、三维、多块/MPI、输出、停止、故障和
跨 rank 重启。随后完成二维 Riemann、三维泊肃叶流、扭曲网格等熵涡、双马赫反射和三维槽道
流迁移/长算等人工检查。槽道长算用于并行稳定性与统计链路验证，不构成高分辨率 DNS 数据集。

## 发行内容

开发仓库保留完整测试、生成器、人工算例结果和阶段记录。独立 `wcns_v3_release` 仓库仅保留：

- `include/`、`src/` 和顶层构建文件；
- 构建所需的 CGNS 4.4.0 源码归档；
- 项目介绍、用户手册、运行速查、算法说明、开发指南和已知限制；
- 两个带小型 CGNS 网格的可运行示例。

## 兼容性与限制

配置 schema 仍为 1。检查点只承诺相同网格和数值签名下的 WCNS 内部兼容性，不是第三方交换
格式。当前没有稳定库 ABI、湍流模型、低 Mach 预处理、隐式推进或并行 CGNS I/O。全部边界见
[`known-limitations.md`](known-limitations.md)。

WCNS 自有代码尚未选择对外开源许可证，精简源码仓库默认保持私有；CGNS 继续遵循其独立许可。
