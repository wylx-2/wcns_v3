# v1.1.0 阶段 S 设计冻结：输运配置与数值收敛验证

状态：2026-09-11 在 R→S→T 连续推进授权下自动冻结。该授权替代 S 前后的人工等待，
不取消本文件定义的自动卡口。

## 1. 范围与非目标

S 只把现有常黏度和 Sutherland 模型接入生产配置、求解器、输出、检查点签名和验证矩阵，
不增加新的本构关系，不改变无粘离散、`1/Re` 的放置或 SSPRK3 公式。NACA0012 用户目录不在
本阶段范围内。

## 2. 配置接口、单位和合法域

公共键为：

```text
transport.model = constant | sutherland
transport.prandtl = <finite positive>
```

`constant` 唯一允许的模型专属键为：

```text
transport.constant.viscosity_ratio = mu_const*/mu_ref
```

`sutherland` 要求参考温度黏度比，并要求以下两种 Sutherland 常数表示恰好出现一种：

```text
transport.sutherland.reference_viscosity_ratio = mu(T_ref)*/mu_ref
transport.sutherland.temperature = S*              # K
# 或
transport.sutherland.temperature_ratio = S*/T_ref # 无量纲
```

两种表示同时出现、模型专属键错配、未知模型、重复键、非有限值或非正值均在配置解析/验证期
失败。Kelvin 表示在解析后除以 `reference.temperature`，因此后续只保存无量纲
`constant_temperature_ratio`。未出现任何 `transport.*` 键时固定迁移为
`constant, viscosity_ratio=1, Pr=0.72`；只出现部分公共键时其余公共默认值仍适用，但选中模型
必须满足其专属键规则。常黏度专属键可缺省为 1；Sutherland 的两个物性量必须显式给出，避免
把空气经验常数静默用于其他气体。

## 3. 数学与尺度

模型公式唯一引用根目录《算法补充》11.5。常黏度为
`mu(T)=mu_ratio`；Sutherland 为
`mu(T)=mu_Tref_ratio*T^(3/2)*(1+S)/(T+S)`，其中温度均相对 `T_ref` 无量纲化。
热传导系数仍为 `chi=mu/[(gamma-1) Ma^2 Pr]`。粘性面通量先使用上述无量纲系数，
`1/Re` 只在粘性散度和粘性时间步限制中各出现一次，配置层不得预乘。

## 4. 数据流、启动诊断与重启

`CaseConfig` 拥有唯一 `TransportConfig`，同一个配置副本构造求解器 `TransportModel` 和
`QuantityContext`，保证场输出、边界输出、统计与求解器使用同一有效黏度。root 启动摘要打印
模型、Pr、`mu(1)`、`S/T_ref`（常黏度为 `not_applicable`），并在初始化/恢复后的全局真实单元
最小、最大温度处打印黏度范围。扫描会调用真实模型，因此非法运行温度在推进和正常输出前由
所有 rank 一致失败。

输运摘要进入配置摘要和 manifest；规范化后的模型、Pr、黏度比和 `S/T_ref` 无条件进入
restart signature。K 与无量纲比值若表示同一物理参数，签名完全相同；任何有效参数改变都必须
在恢复推进前被检查点服务拒绝。v1.0 未配置输运参数的检查点与缺省迁移配置保持旧默认签名
语义下的兼容路径，新增显式非默认输运不得伪装成兼容。

## 5. 验证矩阵和冻结阈值

| 卡口 | 自动证据 | 阈值 |
|---|---|---|
| S-V1 配置/公式 | 点值、T=1、尺度等价、互斥/非法键 | 点值相对误差 `<=1e-14`；所有非法输入拒绝 |
| S-V2 热系数/CFL | `chi` 与粘性 dt 手算对照 | 相对误差 `<=1e-13`；Re 只出现一次 |
| S-V3 旧配置 | 未增加键的 v1.0 配置与 S 前基线 | 配置签名与数值场逐位一致 |
| S-V4 空间收敛 | 2D/3D、两 profile、至少三层网格 | 每组保存内部/边界/全局 L1/L2/Linf；渐近阶 `>=1.8` |
| S-V5 时间收敛 | 周期光滑问题至少四档 dt | 最后三档拟合/逐级观测阶 `>=2.8` |
| S-V6 解析流 | Couette、导热，constant/Sutherland | 既有解析阈值不放宽；串行/MPI场差 `<=5e-10` |
| S-V7 重启/并行 | 单/多块，1/2/4 rank，同/异 rank 重启 | 场差 `<=5e-10`；参数改变恢复失败且无正常末态 |
| S-V8 CFL 扫描 | 无粘、黏性、强激波、曲线网格 | 报告每档成功/失败、降阶、缩步边界，不外推为定理 |
| S-V9 全回归 | Release 串行/MPI CTest、规格检查、diff | 全部通过；工作树仅保留已登记用户目录 |

现有工程还没有闭式粘性制造源来保证高阶精确解，S 的空间矩阵必须输出三类误差原始表；若
边界闭合把全局阶限制到二阶，仍以 `1.8` 为卡口且不得只展示内部高阶。时间阶矩阵采用空间
离散误差至少低一个数量级的周期光滑解，并保存所有 dt/误差数据。

## 6. 提交与停止条件

设计、配置接线、公式/失败测试、收敛与 CFL 工具、文档、正式验收分别独立提交。每个提交先
运行最小相关测试和 `git diff --check`。S 全卡口通过后创建不可移动的
`v1.1.0-s-candidate.N`，以 `--no-ff` 合入 `release/v1.1.0` 并直接进入 T；失败则留在 S 修复，
不得通过放宽阈值继续。
