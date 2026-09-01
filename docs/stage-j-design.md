# 阶段 J 设计：无粘 WCNS、边界状态与显式源项

状态：候选 v1 已于 2026-09-02 通过正式自动验收并获项目负责人批准，允许进入阶段 K。

## 1. 目标与范围

阶段 J 在阶段 I 两套静态几何 profile 上实现无粘高阶空间离散和单元局部显式源项。保留阶段 E 有限体积 Euler 路径作为独立基线，不允许将其静默用作 WCNS profile 的替代。

本阶段不实现黏性梯度、应力、热流、HLLC、Roe、特征重构、动态网格或非共形连接；这些分别属于阶段 K/L/M。

## 2. J0：配置与公共状态语义

- 固定 `SpatialScheme=finite_volume_baseline | wcns`、`Reconstruction=linear5 | wcns_js`、`ReconstructionVariables=conservative | primitive`、`RiemannSolver=rusanov`。
- WCNS 配置必须携带阶段 I 的完整 `AlgorithmProfile`；profile、度量、差分和面通量 halo 层宽不可交叉组合。
- 新路径使用温度型 primitive `(rho,u,v,w,T)` 作为物理边界 ghost 真值，并在同一操作中派生压力型 primitive 与守恒量。
- `strong_boundary_face_state` 默认开启，配置摘要和重启签名必须包含该开关及源项顺序。

## 3. J1：重构与正性回退

- 实现文档固定的五阶线性左右重构和带 `ReconstructionScaling` 的 `wcns_js`。
- 第一版支持 conservative 或 pressure-primitive 分量重构；characteristic 明确拒绝并留待阶段 L。
- 每个面按“所选非线性方案、linear5、相邻单元一阶”检查有限性、正密度和正压力；每级回退具有确定性计数。
- 重构只沿面法向访问三层状态，切向索引始终位于真实单元范围，不读取边/角 ghost。

## 4. J2：Riemann 与真实面通量

- `RiemannSolver` 只接收单位法向并返回每单位物理面积通量；首个实现仅为 Rusanov。
- 调用方用真实面面积矢量生成单位法向，并且只乘一次面积。
- 无通量分裂接口；非单位法向、非物理状态和非有限通量确定性失败。
- `InviscidFaceFluxField` 只拥有真实面和块连接所需法向层，不拥有物理边界外面通量。

## 5. J3：公共物理 ghost 与边界面状态

- `PhysicalGhostStateOperator` 对每个 patch 仅填法向三层、切向真实范围的面状 slab；边和角保持 NaN/非法。
- 每层先产生 `(rho,u,v,w,T)`，再在同一事务派生 `p,U` 并共享版本。
- 首版覆盖已有 CGNS 边界枚举：入口/远场目标状态、出口外推、滑移/对称反射、无滑移绝热/等温壁；缺少必需 patch 数据失败。
- `InviscidBoundaryFaceState` 在重构后按统一开关处理外侧迹。固壁强修正仅反射法向速度并保持内部迹的密度、压力和切向速度；弱模式保留原始重构迹。

## 6. J4：高阶散度与 `FaceFluxHalo`

- PH 使用 `D4-PH`，连接面层为 0--1；SCMM6 使用 `D6/D4`，连接面层为 0--2。
- 共享连接面由阶段 I 相同所有者计算；内部层由其所在块计算，接收块不得重复重构。
- 同进程和 MPI 复用同一面索引/变换计划。普通连接只重排逻辑索引和通量方向；旋转周期额外旋转三个 Cartesian 动量通量分量。
- 面通量版本必须匹配当前 RK 子步和 profile，散度前检查 halo 完整性。

## 7. J5：显式局部源项

- 实现 `uniform_conservative`、`body_force` 和可注入回调的 `manufactured`。
- 模型只读取当前真实单元的 `U`、冻结中心坐标和当前 RK 子步时间；源项直接累加到 `Residual=dU/dt`，不再乘 Jacobian。
- 二维源项严格保持 z 动量分量为零。体力能量源固定为 `rho*u dot a`。
- 关闭源项时不创建模型、不分配临时源字段，并与基线逐位一致；模型顺序和参数进入摘要/重启签名。
- SSPRK3 残差回调显式接收子步时间，三个阶段分别使用 `t_n`、`t_n+dt`、`t_n+dt/2`。

## 8. 分段 Git 节点

1. J0：本设计和配置冻结；
2. J1：温度型公共 ghost、边界面强/弱约束及测试；
3. J2：linear5、scaled WCNS-JS、Rusanov 接口和正性回退；
4. J3：真实面通量、profile 散度和 `FaceFluxHalo`；
5. J4：三个显式源项、RK 子步时间和关闭路径兼容性；
6. J5：单块/多块/MPI 数值回归及正式验收报告。

每个小目标测试通过后提交并尝试同步 GitHub。失败的正式验收按项目 Git 规则记录；局部开发测试失败不单独提交。

## 9. 正式验收卡口

- linear5 精确重构四次多项式；scaled WCNS-JS 对常量/光滑数据正确，并能确定性降级。
- Rusanov 自由流通量等于解析 Euler 通量，单位法向与面积只处理一次。
- 强壁面模式逐面零法向质量通量且无无粘切向耗散；弱模式路径可辨识。
- 物理边界边/角 ghost 保持 NaN，几何 ghost 从未创建或读取。
- 两套 profile 的单块恒定流残差满足舍入误差；相同网格的串行多块和双 rank 结果一致。
- PH 层 0--1、SCMM6 层 0--2 面通量 halo 的索引、方向、周期旋转、消息身份和版本通过测试。
- 关闭源项与阶段 E 无源基线逐位一致；常量源、体力能量功、制造源的 RK 时间和全局体积加权平衡通过。
- 阶段 G 规格检查、串行 Release、MPI Release、`git diff --check` 全部通过，之后等待项目负责人人工验收。
