# 阶段 F 人工复核修订报告（第 1 次）

日期：2026-08-28

状态：修订与自动验算通过，等待人工复查。

## PHengLEI WCNS 核对结论

PHengLEI 在 `isFVMOrFDM == FD_METHOD` 时进入结构高阶度量；`str_highorder_solver` 为 WCNS 时使用非 HDCS 的 `ComputeMetricsStructHighOrder2D/3D`。其核心流程与原文以下部分一致：

- CGNS/结构网格顶点被放入 `2N_v-1` 加密交错工作区；
- 偶数位置由相邻位置算术平均；
- 坐标导数、中间乘积、面矢量和 Jacobian 使用同一 `GridDelta`；
- 面矢量采用 SCMM 等价的对称守恒形式；
- Jacobian 用 $\frac13\nabla_\xi\cdot(\boldsymbol r\cdot\boldsymbol S)$ 计算，并存为物理体积尺度。

原文与 PHengLEI 实际实现存在以下差异，现已修正：

- 原文把输入坐标描述为单元中心坐标；PHengLEI 实际从顶点坐标开始，偶-偶-偶位置才是单元中心。
- PHengLEI WCNS 的内部通量差分为 $9/8-1/24$ 四阶形式，物理边界相邻单元降为相邻两面差，并非原文的 D6。
- PHengLEI 在块连接通信后另行重算连接邻域度量和 Jacobian；Jacobian 异常时可回退到有限体积体积。
- PHengLEI 会镜像生成部分 ghost 单元中心/体积。本项目根据人工约束不继承该策略。

## 本次设计修订

- 增加受约束的 `phenglei_wcns | scmm6_wcns` profile，避免度量与不相容通量差分任意组合。
- 物理边界 ghost 只允许边界条件生成的物理状态；坐标、度量、梯度、输运系数、通量和残差均无 ghost 合法区。
- 无粘通量删除通量分裂路线，只保留可切换的 Riemann 求解器。
- `scmm6_wcns` 使用共同中心工作网格，并强制 SCMM 各嵌套层次满足 $\delta^1=\delta^2=\delta^3$。
- CGNS 顶点到 SCMM6 单元中心改用六点五次 Lagrange 插值；不再使用会限制几何阶数的四/八顶点算术平均。
- 增加真实边界面状态顺序、Riemann 面通量操作、边角状态算子和无 ghost 粘性梯度流程。
- 增加 $D^Tw=e_R-e_L$ 的加权全局守恒条件，并明确多块高阶接口需要全局耦合守恒验收。
- PHengLEI 有限体积 Jacobian 回退改为显式配置，默认 `strict`，禁止静默替换。

## 自动验算

- 原有 $I_6$、两端插值闭合、$D_6$ 和两端导数闭合的精确多项式矩条件继续通过。
- 新增顶点到中心的三组六点插值对 0 至 5 次多项式全部精确。
- PHengLEI 四阶内部差分对 0 至 4 次多项式精确，边界相邻两面差满足其声明阶数。
- 对 `phenglei_wcns` 的 $N=2\ldots128$ 和 `scmm6_wcns` 的 $N=5\ldots128$，均得到正的守恒积分权重，且 $\|D^Tw-(e_R-e_L)\|_\infty<2\times10^{-11}$。
- 使用包含两端单边行的三维张量积 SCMM 算子，在随机坐标数据上检查离散面守恒律，相对抵消误差为 $1.97\times10^{-16}$。
- 九个 SCMM 分量和对称 Jacobian 的随机仿射映射精确验算继续通过。
- `git diff --check` 通过；本次未修改 C++，未执行求解器回归测试。

## 后续人工复查重点

- 是否认可两个 profile 同时约束度量方法和通量差分。
- 是否认可 `scmm6_wcns` 使用六阶顶点到中心插值，而 `phenglei_wcns` 为对照保留算术平均。
- 是否认可用正积分权重定义边界闭合下的全局守恒量。
- 多块高阶接口采用全局耦合算子还是保守接口闭合，仍需在实现该子任务前选定并单独验收。

人工复查通过前不进入下一阶段或修改求解器实现。
