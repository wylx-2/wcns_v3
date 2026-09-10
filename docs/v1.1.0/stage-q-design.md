# v1.1.0 阶段 Q 设计冻结：SSPRK 事后稳健化

状态：2026-09-11 自动卡口通过，候选 `v1.1.0-q-candidate.1` 等待 Q→R 人工验收；未进入 R。

## 1. 目标、权威公式与非目标

本阶段实现《算法补充》11.1--11.3，尤其以 11.2.3 的 troubled cell 与离散支持传播为权威
算法。目标是在每个 SSPRK3 候选提交前验证物理容许性，升级能影响失败单元的无粘面通量及
一层转置离散支持保护带；
共享面仍只有 owner 计算的一份权威通量。禁止状态截断、全域静默降阶、改变 profile 差分、
关闭黏性/源项或把 rejected step 写入输出/checkpoint。

Q 不包含 NACA0012，也不实现新的激波传感器、保正缩放器、人工黏性、隐式推进、AMR 或并行
CGNS。Case07 四块 48x32 Mach 5 圆柱是生产验证算例。

## 2. 配置冻结

新增可选键及默认值：

```text
robustness.enabled = false
robustness.max_local_recomputations = 3
robustness.max_step_retries = 4
robustness.time_step_reduction = 0.5
robustness.minimum_time_step = 1.0e-12
```

全部键影响数值轨迹，进入配置摘要和 restart signature。缺省 `enabled=false`，因此 v1.0 配置、
求解路径和逐位结果保持不变。整数上限非负；缩步因子严格位于 `(0,1)`；最小时间步为有限正数。

每个 RK 阶段和每次整步 retry 均把面级别重新初始化为用户配置级别 0，不跨阶段缓存提示。
这样相同 `(U^n,dt,配置,划分)` 不依赖上一次失败的本地历史。达到上限时拒绝整步；第 `r+1`
次尝试使用 `dt_{r+1}=0.5 dt_r`。最终采用的 `dt` 返回运行驱动，只有 accepted step 增加时间。

## 3. 数据与接口

- `RobustnessConfig`：验证、摘要和数值重启签名。
- `CandidateState`/`TroubledCell`：按真实单元存放尚未提交的候选及失败量；检测不写权威状态。
- `FaceRobustnessField`：每块、每方向单分量整数面级别；真实面加 PH 1 层/SCMM6 2 层连接 halo。
- `inviscid_residual_stencil(...)`：残差与传播共同调用的唯一面导数行查询。
- `FaceRobustnessExchanger`：沿 `FaceFluxHaloPlan` 反向把 receiver 请求以 `max` 合并到 donor
  真实面；随后正常通量交换正向发布 owner 结果。
- 求解器 `advance` 返回 accepted `dt`；稳健模式关闭时继续调用原 SSPRK3 路径。

MPI 级别消息包含 plan version 和每个 descriptor 的固定数量级别，使用独立 tag base。级别为
无方向整数，不做周期旋转或符号变换。同 rank 与跨 rank 使用相同 pair 列表和 `max` 语义。

## 4. 有效降阶梯子

概念梯子为：

1. 用户配置 scheme/variables/riemann；
2. 同 scheme 的 primitive 变量；
3. `linear5 + primitive`；
4. `zero_order + conservative + rusanov`。

运行时比较 `(scheme,variables,riemann)` 删除重复组合；例如原配置已经是
`linear5 + primitive` 时，下一级直接是最终零阶 Rusanov。每轮同一物理面最多前进一个有效
级别，同一阶段只升不降。物理边界仍先产生对应级别的内外迹，再解 Riemann 问题。

## 5. 候选、传播与重算顺序

每阶段保存 `U^n` 和本阶段输入，不直接覆盖后者。形成候选缓冲后遍历真实单元：五个守恒分量
有限，且 `rho/p/T/e` 达到统一 `NumericalFloors` 才可一次性提交。失败时：

1. 由权威残差 stencil 得到每个 troubled cell 的全部直接面支持；
2. 查找残差支持与直接面相交的真实单元，再把这些单元的完整支持并入一层保护带；
3. 本地面请求前进一个有效级别；
4. 所有 rank 反向交换请求并在 owner 真实面取最大值；
5. 若全局没有 owner 面变化，则本阶段不能再局部修复，拒绝整步；
6. 以未改变的阶段输入重新计算状态/边界/完整面通量/halo、残差和候选；
7. 候选合法则提交，否则继续至局部上限。

首版为保证清晰和确定性，升级后重新计算该阶段完整无粘面场及完整残差；未升级面必须得到与
上一轮相同的数值。阶段 T 才允许缓存未变面，且必须证明等价。

## 6. 诊断

每个 accepted step 保存：每级 owner 真实面计数、troubled cell 总数、局部重算轮数、整步
retry 数、初始/接受 dt、最小 `rho/p/T/e`。失败异常包含 rank、block、`i/j/k`、RK stage、
stage time、失败量和最终 dt。共享面只在 owner 计数。运行 history/manifest 追加稳定字段，
但诊断输出频率不进入重启签名。

## 7. 自动验收映射

- 单元：配置合法/非法/default/signature；候选有限性与 floors；PH/SCMM6/二点 stencil 精确面集；
  有效梯子去重和单调终止；SSPRK 候选提交/恢复/缩步。
- 多块/MPI：连接 halo 请求反向传播、owner 最大级别、2/4 rank 计数和状态一致，无 deadlock。
- 回归：`enabled=false` 时 P 的 48/84 CTest 全通过；自由流不降阶；既有光滑涡误差不退化。
- 强激波：Sod、双马赫 smoke 和 Case07 48x32；Case07 不允许全域配置 `zero_order`，必须从
  高阶配置启动并记录局部最终级别、retry、最小状态与守恒历史。
- 失败：达到局部/整步上限或低于 `minimum_time_step` 时返回 numerical failure；rejected step
  不增加 step/time，不产生正常终场/checkpoint。

正式验收记录构建、测试数、Case07 配置/时间、各级比例及远程状态。自动通过后创建
`v1.1.0-q-candidate.N`，但不合并 `release/v1.1.0`，等待项目负责人对 Q→R 的人工判断。
本次自动验收证据见 [`stage-q-acceptance.md`](stage-q-acceptance.md)。
