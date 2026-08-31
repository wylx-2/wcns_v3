# 阶段 G 自动验收报告

日期：2026-08-31

状态：算法规格、精确公式核验及既有程序回归通过，等待人工检查。

## 本阶段范围

阶段 G 只审查、修订并冻结 `算法补充.md`，补充可重复执行的规格核验工具；没有修改 `include/`、`src/`、`tests/` 或 `CMakeLists.txt`，也没有把 WCNS 算法提前接入现有求解器。

## 已完成内容

- 将 `algorithm_profile = phenglei_wcns | scmm6_wcns` 定义为两套独立且不可拆分的完整 profile；度量、中心到面插值、通量差分和物理边界闭合禁止交叉组合。
- 补齐无量纲守恒量、Euler 通量、量热完全理想气体 `GasModel`、参考量一致性、CGNS 长度单位、数值阈值和重构变量尺度。
- 无粘通量只保留 Riemann 求解器路线；明确 Riemann 求解器输出为单位法向的每单位面积通量，计算空间通量必须再乘真实面面积。
- 将边界处理拆成 `PhysicalGhostStateOperator`、`InviscidBoundaryFaceState` 和 `ViscousBoundaryTrace`。`strong_boundary_face_state` 只控制重构后的无粘外侧迹，不能关闭无滑移、等温或绝热粘性壁面条件。
- 审计确认当前逐方向重构和梯度流程不读取物理边界边/角状态 ghost；只允许边界条件生成的三层面状状态 ghost，物理 ghost 坐标、度量、梯度和通量均为非法数据。
- 区分 `AdjacentCellRange` 与 `BoundaryFaceRange`，并为块连接定义相邻单元范围、唯一共享面范围、轴变换、周期旋转/平移和确定性所有者。
- 明确高阶块连接在状态 halo 之外还必须同步 `FaceFluxHalo`：PH profile 使用层 0--1，SCMM6 profile 使用层 0--2；粘性梯度连接 halo 分别为 2 层和 3 层。
- 明确两套度量的嵌套离散不能只交换一次坐标：每个 `GeometryOperand` 必须按算法阶段补齐连接 halo；需要跨两个或三个连接的模板必须解析唯一 donor 路径，不能穿过物理边界。
- 补齐 PH profile 的 $I_4^{PH}$、两套壁面 Dirichlet 法向导数生成规则、温度面值与 Sutherland 输运系数的计算顺序。
- 给出单块和多块全局守恒权重的离散条件，并把唯一性、严格正性和残差检查定义为运行前置条件。
- 冻结字段的拓扑位置、合法索引域、生产者、消费者、生命周期、交换消息种类和版本规则，并将后续实现重新划分为阶段 H--M。

## 事实、公式与可行性结论

- 文档中的 $J$ 始终表示无量纲物理坐标映射的 $\partial(x,y,z)/\partial(\xi,\eta,\zeta)$；它无量纲，但数值上具有单位计算单元面积/体积意义。有量纲体积因子为 $L_{ref}^dJ$。
- 动压压力尺度下的 $p=\rho T/(\gamma Ma^2)$、$e=T/[\gamma(\gamma-1)Ma^2]$ 与 $a=\sqrt{\gamma p/\rho}$ 相互一致；粘性散度中的 $1/Re$ 只出现一次。
- 文档列出的 $I_4^{PH}$、$I_6$、顶点到中心 $I_6^{v\to c}$、$D_4^{PH}/D_2^{PH}$、$D_6$/单边闭合系数均满足其声明的精确矩条件。
- 在当前逐方向算法假设下，不填物理边界边/角状态 ghost 是合理且可实现的；若未来加入会读取多法向 ghost 的多维重构或滤波，必须重新设计 corner operator，不能扩展现有填充器的隐含职责。
- 连接界面作为模板完整的内点仍需交换半节点通量值：状态 halo 足够形成接口 Riemann 状态，但不足以保证高阶散度所需的邻块内部面通量具有唯一权威值。
- 多块正守恒权重不是对任意拓扑无条件成立的数学假设。实现必须对实际组装的 $\mathcal D^TW_c=B^TW_b$ 检查一致性、唯一性、正性和残差；不满足时拒绝当前 topology/profile，而不能静默改用等权或局部权重。

## 自动验收结果

执行：

```text
python tools/verify_algorithm_spec.py
```

结果为 6 组全部通过：

1. 插值与一阶差分的精确有理数矩条件；
2. PH 四中心和 SCMM6 六中心壁面法向导数系数生成；
3. PH 的 $N=4\ldots12$、SCMM6 的 $N=5\ldots12$ 单块正守恒权重及 $D^Tw=e_R-e_L$；
4. 一般三维仿射映射的 Jacobian/余因子恒等式；
5. 无量纲压力、内能、Reynolds/Mach 尺度关系；
6. 两 profile、三边界职责、`FaceFluxHalo`、`GeometryOperand` 和禁用通量分裂等冻结文档契约。

既有程序回归：

- `cmake --build build --config Release`：通过；
- `ctest --test-dir build -C Release --output-on-failure`：6/6 通过；
- `cmake --build build-mpi --config Release`：通过；
- `ctest --test-dir build-mpi -C Release --output-on-failure`：9/9 通过，其中包括双 rank halo、Euler 和 MPI runtime 测试；
- `git diff --check`：通过。

## 后续实现前仍需实证的项目

以下内容已被定义为接口或验收条件，但尚未由当前 C++ 程序实现，不能把阶段 G 的文档通过误认为求解能力已经具备：

- 两套高阶 `MetricMethod` 及其分阶段几何通信；
- 实际 CGNS 多块拓扑上的全局正守恒权重；
- WCNS 非线性重构、Riemann 面通量和高阶差分；
- 特征入口/出口、物理壁面梯度和输运模型；
- 旋转周期连接的状态、梯度、面矢量和通量变换；
- 粘性时间步系数 $C_v$ 的 profile/维数/时间推进组合标定。

人工检查通过前不进入阶段 H。

## 2026-08-31 源项补充

阶段 G 候选完成后按批示补充了方程源项规格，仍未修改求解器 C++：

- 定义 $\mathcal S$ 为单位物理体积内五个守恒量的无量纲变化率；映射守恒方程使用 $J\mathcal S$，当前 `Residual=dU/dt` 只直接累加 $\mathcal S$，禁止重复乘 Jacobian。
- 增加 `enable_source_terms=false|true` 和有序 `source_models`；默认关闭，关闭路径不分配或调用源项对象。
- 给出质量、动量和总能源的有量纲到无量纲尺度，以及体力的动量/能量功一致形式。
- 规定非刚性局部源项在每个 Runge--Kutta 子步按当前 $U^{(s)},\boldsymbol x_c,t^{(s)}$ 重新计算，只访问真实单元且不隐含 MPI 通信。
- 把全局守恒验收扩展为物理边界净通量与 $\sum W_cJ_c\mathcal S_c$ 的共同平衡，并增加关闭等价性、常量源、体力、制造解、多块/MPI 和重启测试要求。
- 阶段 H 将建立 `SourceTermConfig/Registry` 空框架；阶段 J 将实现 `uniform_conservative`、`body_force`、`manufactured` 及逐 RK 子步装配。

`tools/verify_algorithm_spec.py` 已增加体力源动量/能量无量纲关系和源项文档契约检查，全部规格核验继续通过。

## 2026-08-31 Re/Ma 输入与公共 ghost 修订

- `Re`、`Ma` 改为只能由五个参考量和唯一 `GasModel` 输入派生；输入文件出现任一字段即失败，不再执行 Re/Ma 输入一致性断言。派生结果写入 log 和检查点供核查。
- `molar_mass` 与 `specific_gas_constant` 改为严格二选一，从源头避免两套气体常数输入的一致性问题。
- 物理边界采用唯一权威 primitive ghost：边界条件逐层给出 $\rho,u,v,w,T$，随后在同一事务中计算 $p$ 和 $U$，两者共享同一版本。
- 无粘重构和粘性原始量导数读取同一份 ghost；仍不允许 ghost 坐标、度量、梯度或通量，也不需要边/角 ghost。
- 采用“公共 ghost + 分开的真实面修正”方案：无粘面状态按 `strong_boundary_face_state` 决定是否覆写；粘性壁面迹在生产算法中始终强修正。
- 不采用“同一 ghost 直接作为两种最终边界结果”的方案。无滑移 ghost 的切向速度跳跃若直接进入 Rusanov/HLLC，会产生非物理无粘切向耗散；默认无粘强修正会把固壁外侧迹改为只反射法向速度，从而移除该影响。
