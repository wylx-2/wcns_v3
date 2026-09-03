# WCNS 0.1.0 发布候选说明

状态：阶段 O5 资料候选。只有 O6 在全新构建目录完成冻结发布矩阵、创建候选标签并由项目
负责人集中验收后，本文才代表一个可接受的发布候选。

## 主要能力

- 二维/三维结构多块 ADF-CGNS 读取，共形普通/平移周期/旋转周期连接；当 zone 少于 MPI
  rank 时，可按高阶模板下限确定性二次剖分。
- 独立的 `phenglei_wcns` 与 `scmm6_wcns` 几何/线性算子 profile；WENO-JS、WENO-Z、
  MDCD-LINEAR、MDCD-HYBRID 重构及 Rusanov、HLLC、Roe 求解器。
- 热完全 Euler 和层流 Navier--Stokes 显式 SSPRK3 推进，以及均匀守恒源、体力与制造源。
- 定常残差连续判据、非定常物理时间、最大步数、墙钟和信号安全停止。
- CGNS/Tecplot 场、TXT/Tecplot 历史/统计、原子 manifest/检查点，以及跨 rank/合法叶块划分
  的无损重启。
- 正式生产入口驱动的解析、守恒、串并行等价、输出/重启和确定性失败矩阵。

## 构建与安装

已验证工具链与完整命令见 [`runtime-guide.md`](runtime-guide.md)。`cmake --install` 安装
`wcns_run`、网格生成器、独立验证器、配置模板和文档；上游 CGNS 自身的安装规则还会安装
其静态库、头文件和工具。编译器与 MPI 运行库须由目标环境提供。

## 基线和兼容性

提交级数值与代表性性能记录冻结在
[`cases/baselines/stage-o-quick-v1.json`](../cases/baselines/stage-o-quick-v1.json)。性能值只用于
同环境回归参考，不是跨机器速度门槛。检查点格式版本和数值签名受严格校验，但当前不提供
稳定 WCNS 库 ABI，也不保证检查点与第三方软件互换。

所有已知模型、网格、I/O、扩展、平台与许可限制见
[`known-limitations.md`](known-limitations.md)。尤其需要注意：当前仅构建 ADF-CGNS、没有低
Mach 预处理/湍流/隐式推进/并行 I/O，且 WCNS 自有代码尚未选择对外许可证。

## O6 前仍需完成

- 从空目录重建串行和 MPI Release，并运行全部 CTest；
- 执行冻结的中等规模自由流、等熵涡、间断、黏性、源项、三维、输出/重启和失败矩阵；
- 记录最终提交、工具链、命令、误差、守恒、负载、墙钟和峰值 RSS；
- 验证安装树、工作树/远程同步，创建 `stage-o-release-candidate-v1` 标签后停止等待人工验收。
