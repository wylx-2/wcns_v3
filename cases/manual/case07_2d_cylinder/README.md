# Case07：二维圆柱绕流与 Mach 5 钝体绕流

## 1. 算例目的与结论边界

Case07 用同一套四块结构 O 网格和 WCNS 正式生产入口 `wcns_run` 检查两类外流问题：

1. `Re=20`、`Re=40`：亚临界二维层流的稳定、近似对称回流区；
2. `Re=100`、`Re=200`：超过圆柱尾迹首次全局失稳阈值后的非定常尾迹和交替涡脱落；
3. `Ma=5`：无粘高超声速圆柱钝体前方的脱体弓形激波、驻点压缩和背风膨胀区。

本目录中的实际结果是**当前程序的粗网格功能/定性验收结果**，不是网格收敛的工程数据。
尤其需要注意：

- 低速算例采用可压缩层流 Navier--Stokes，来流 `Ma=0.2`，而引用的经典数据通常采用
  不可压缩方程；
- 实际低速网格只有 `32 x 20` 个周向/径向单元，远低于参考计算；
- v1.1 配置启用真实边界面 `Cp`、压力/黏性牵引、热流及力/力矩积分；逐面量来自与求解器
  一致的高阶壁面迹，面积使用全局守恒边界权重，并支持切分后的确定性 MPI 汇总；
- `Re=200` 仍按二维方程计算，只用于清楚显示数值涡街，不能代替真实三维转捩尾迹；
- Mach 5 算例是 Euler 滑移壁计算，因此没有黏性阻力、边界层或气动热。

圆柱尾迹的首个全局失稳临界值约为 `Re=47`，这解释了为什么 Re=20/40 应趋向稳定对称态，
而 Re=100/200 应出现周期性涡脱落。参见
[Cantwell 与 Barkley (2010)](https://doi.org/10.1103/PhysRevE.82.026315)。

## 2. 控制参数与无量纲化

圆柱直径为参考长度：

\[
D=L_{ref}=1,
\qquad
Re_D=\frac{\rho_\infty U_\infty D}{\mu_\infty},
\qquad
Ma_\infty=\frac{U_\infty}{\sqrt{\gamma R T_\infty}}.
\]

所有算例采用热完全理想气体：

\[
\gamma=1.4,
\qquad
\widetilde M=0.029\;\mathrm{kg/mol},
\qquad
R=\frac{8.314}{0.029}=286.689655\;\mathrm{J/(kg\,K)}.
\]

取 `rho_ref=1`、`T_ref=300 K`、`D=1 m`。低速组采用

\[
U_{ref}=U_\infty=69.4001888\;\mathrm{m/s},
\qquad Ma_\infty\simeq0.2,
\]

并按目标雷诺数计算 `mu_ref=rho_ref U_ref D/Re`。Mach 5 组采用

\[
U_{ref}=1735.004719\;\mathrm{m/s},
\]

配置中的 `Re=10^6` 仅为满足五参考量输入而选取；Euler 方程不使用该黏性尺度。

由于 `U_ref=U_inf`，无量纲来流速度为 `(u,v)=(1,0)`、`rho=T=1`，来流压力为

\[
p_\infty=\frac{1}{\gamma Ma_\infty^2}.
\]

因此 Ma=0.2 时 `p_inf=17.857142857`，Ma=5 时 `p_inf=0.0285714286`。

## 3. 算例矩阵

| 配置 | 方程 | `Re_D` | `Ma` | 网格 | 重构/变量 | Riemann | CFL | 终止时间 |
|---|---|---:|---:|---|---|---|---:|---:|
| `cylinder_re20.wcns` | 层流 NS | 20 | 0.2 | 32x20 | linear5/primitive | HLLC | 0.25 | 15 |
| `cylinder_re40.wcns` | 层流 NS | 40 | 0.2 | 32x20 | linear5/primitive | HLLC | 0.25 | 20 |
| `cylinder_re100.wcns` | 层流 NS | 100 | 0.2 | 32x20 | linear5/primitive | HLLC | 0.25 | 30 |
| `cylinder_re200.wcns` | 层流 NS | 200 | 0.2 | 32x20 | linear5/primitive | HLLC | 0.25 | 20 |
| `cylinder_mach5_euler.wcns` | Euler | -- | 5 | 48x32 | zero_order/conservative | Rusanov | 0.02 | 8 |
| `cylinder_mach5_robust.wcns` | Euler | -- | 5 | 48x32 | WENO-Z/characteristic | Rusanov | 0.02 | 8 |

全部计算使用独立的 `scmm6_wcns` profile，没有与 `phenglei_wcns` 交叉使用度量或差分算子。
低速流动平滑，选择五阶线性重构以减小本次实际计算成本；v1.0 强激波基线采用六点模板兼容
的全域零阶重构和 Rusanov 通量。v1.1 稳健算例采用高阶 WENO-Z 起步、守恒两点通量差分和
事后局部离散支持降阶，不把全域方案预先改成零阶。

Re=100 初场设置 `v=0.01`，Re=200 设置 `v=0.02`，外远场仍严格保持 `v=0`。这一小扰动只用于
缩短全局不稳定模态从舍入误差增长到可见振幅所需的时间，不能解释为持续横向来流。

## 4. 四块结构 O 网格

本任务向 `wcns_generate_release_cgns` 增加了通用模式：

```text
wcns_generate_release_cgns cylinder-o output.cgns \
  cells_theta cells_radial zones_theta diameter outer_radius radial_cluster_strength
```

网格采用

\[
\theta=-2\pi\xi,
\qquad
r(\eta)=\frac{D}{2}+
\left(R_o-\frac{D}{2}\right)
\frac{\exp(\beta\eta)-1}{\exp(\beta)-1}.
\]

负号使 `(theta,r)` 的计算坐标次序产生正 Jacobian。周向拆为四个原生 CGNS zone，相邻块和
首尾块均使用成对 `GridConnectivity1to1`；`j-min` 是名为 `cylinder` 的壁面，`j-max` 是
`farfield`。物理 ghost 只由边界条件生成，块接口 ghost 由正式 halo 通信生成。

目录保存三套网格：

| 文件 | 总单元数 | 外半径 | beta | 用途 |
|---|---:|---:|---:|---|
| `cylinder_o_32x20_r8.cgns` | 640 | 8D | 2.5 | 本次低速定性计算 |
| `cylinder_o_48x32_r8.cgns` | 1536 | 8D | 2.5 | 本次 Mach 5 计算/低速复核候选 |
| `cylinder_o_96x48_r10.cgns` | 4608 | 10D | 3.0 | 后续精化复核 |

中等网格的 4-rank `--dry-run` 已确认每个原 zone 分配给一个 rank；粗网格在当前 SCMM6 配置下
最多可形成 8 个合法叶块。

## 5. 边界和初始条件

低速组：

- 圆柱：`no_slip_adiabatic_wall`，壁速为零；
- 外圆：带显式目标状态的 `farfield`；
- 初场：全场均匀来流，Re=100/200 叠加上述极小横向速度；
- 黏度：当前正式配置入口的常黏度模型；
- 计算模式：非定常 SSPRK3，准确裁剪到指定 `t_end`。

Mach 5 组：

- 圆柱：`slip_wall`；
- 外圆：特征远场；
- 初场：全场均匀 Mach 5 来流；
- 方程：Euler；
- 计算模式：非定常推进至弓形激波近似稳定。

所有配置均关闭源项。边界 patch 名和类型可在各 `.wcns` 文件中直接核对。

## 6. 从零复现

以下命令从仓库 `wcns` 根目录执行。

### 6.1 构建

```powershell
cmake -S . -B build-case07 `
  -DWCNS_ENABLE_CGNS=ON `
  -DWCNS_ENABLE_MPI=ON `
  -DWCNS_BUILD_TESTS=ON `
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-case07 --parallel 4
ctest --test-dir build-case07 -R "cylinder_o" --output-on-failure
```

新增回归包含：

- `wcns.generate_cylinder_o_mesh`：实际写出四块 ADF-CGNS O 网格；
- `wcns.run.cylinder_o.dry_run.serial`：正式 reader 读取网格、建立连接、计算 SCMM6 度量；
- `wcns.run.cylinder_o.boundary.serial`：输出闭合圆柱的真实面量和载荷历史；
- `wcns.check.cylinder_o.boundary.mpi_equivalence`：核对 1/2/4-rank 面键、逐面量和载荷等价。

### 6.2 重新生成网格

```powershell
powershell -ExecutionPolicy Bypass -File `
  cases\manual\case07_2d_cylinder\scripts\generate_grids.ps1 `
  -BuildDirectory build-case07
```

### 6.3 启动前检查

```powershell
powershell -ExecutionPolicy Bypass -File `
  cases\manual\case07_2d_cylinder\scripts\run_case07.ps1 `
  -Case re20 -BuildDirectory build-case07 -Ranks 1 -DryRun
```

将 `re20` 替换为 `re40`、`re100`、`re200`、`mach5-euler` 或 `mach5-robust`。MPI 拓扑检查
可把 `Ranks` 改为 2 或 4。`--dry-run` 不推进、不创建正式结果目录。

### 6.4 正式运行

```powershell
powershell -ExecutionPolicy Bypass -File `
  cases\manual\case07_2d_cylinder\scripts\run_case07.ps1 `
  -Case re100 -BuildDirectory build-case07 -Ranks 1
```

当前实际结果使用一个串行 rank。不要在已有结果目录上重复运行：配置采用
`output.allow_existing=false`，防止不同行程的文件混合。需要续算时应复制配置、设置新的
`output.directory`，并增加 `restart.path=<旧 checkpoint.latest.cgns>`。

### 6.5 后处理

```powershell
python cases\manual\case07_2d_cylinder\scripts\analyze_case07.py `
  --case-root cases\manual\case07_2d_cylinder
```

需要 Python、NumPy 和 Matplotlib。脚本读取正式 Tecplot 单元中心流场、边界面文件和
`*.loads.rN.txt`，输出：

- `results/analysis-summary.json`；
- 各算例 `surface-history.csv` 和 `surface-final.csv`；
- 残差、阻力、探针信号、涡量、密度和 Mach 数图。

## 7. 后处理定义

压力系数采用

\[
C_p=\frac{p-p_\infty}{\tfrac12\rho_\infty U_\infty^2}.
\]

程序直接从当前算法 profile 的真实边界迹输出所请求的壁面量。Mach 5 Euler 算例只请求
$p_w$、$C_p$、压力牵引和总牵引，不计算与无粘载荷无关的热学量；黏性算例还输出
$T_w,\mu_w$，并使用强修正后的 $\nabla\boldsymbol u_w$ 和 $\nabla T_w$。每个离散面使用

\[
\Delta A_f=w_{b,f}|\boldsymbol S_f|,
\qquad
\Delta\boldsymbol F_f=
\left[p_w\boldsymbol I-\frac{1}{Re}\boldsymbol\tau_w\right]
\boldsymbol n_{\Omega,f}\Delta A_f
\]

积分得到 `Cd_pressure`、`Cd_viscous`、`Cd_total` 及对应升力、力矩。分析脚本不再从第一层
单元重建壁压或剪切；单元中心流场只用于体场图、尾迹长度和探针。

尾迹长度指标 `reverse_flow_x_max` 是中心线附近 `u<0` 单元的最大 `x/D`，从圆柱中心起算。
Re=100/200 的频率分别从真实边界载荷升力与 `(x/D,y/D)≈(2,0.5)` 的横向速度探针计算：

\[
St=\frac{fD}{U_\infty}.
\]

经典关系给出 Re=100 时 `St≈0.1666`；参考高阶计算给出平均 `Cd≈1.32`。参见
[Stålberg 等 (2006)](https://doi.org/10.1007/s10915-005-9043-y)。

## 8. 实际运行过程中的数值调整

第一次 Mach 5 尝试采用 `WENO-Z + characteristic + Rusanov`、`CFL=0.08`。均匀来流在滑移壁上
瞬时建立强激波，第 29 步圆柱迎风面附近出现负压候选状态，程序按数值失败机制返回 3，且没有
把非法状态写成正常终场。

最终配置改为：

```text
algorithm.reconstruction = zero_order
algorithm.reconstruction_variables = conservative
algorithm.riemann = rusanov
run.cfl = 0.02
```

修改后运行到 `t=8`，停止原因是 `physical_time_reached`，没有 reconstruction/Riemann fallback。
这说明当前高阶重构在“全场均匀流直接撞击钝体”的强启动瞬态下还不够正性保持；零阶结果适合
验证弓形激波拓扑，但不适合精确激波位置或表面压力积分。

v1.1 阶段 Q 重新使用 `WENO-Z + characteristic + Rusanov` 验证事后稳健化。实现迭代保留了
三项失败证据：profile 差分、`CFL=0.08` 在 `t=0.0798667` 失败；改用守恒两点差分后在
`t=0.647820` 失败；进一步使用 `CFL=0.02` 和 8 次缩步，但只传播 troubled cell 的直接支持
时仍在 `t=0.649718` 发生时间步趋零。三者都返回 `numerical_failure`，没有写正常终场。

最后冻结的一层转置离散支持保护带按“直接支持面 → 读取这些面的真实单元 → 这些单元的完整
支持面”扩展，不使用任意几何半径。正式配置为：

```text
algorithm.flux_difference = conservative_two_point
algorithm.reconstruction = weno_z
algorithm.reconstruction_variables = characteristic
algorithm.riemann = rusanov
robustness.enabled = true
robustness.max_local_recomputations = 3
robustness.max_step_retries = 8
robustness.time_step_reduction = 0.5
run.cfl = 0.02
```

该配置在 4 rank 下以 12114 步到达 `t=8`。串行与 4 rank 的 `t=1` 共 1536 个单元、10 列
Tecplot 数值逐值相同；同一步历史除墙钟外只有 MPI 归约造成的 `2.22e-16` 以内 L2 差异。

## 9. 结果概览

六个配置均正常推进到配置的物理时间，停止原因为 `physical_time_reached`，没有触发重构或
Riemann fallback。实际串行运行汇总如下；`Residual` 是最后一步五个守恒方程增量的合成 L2，
非定常算例不以该值收敛为目标。

| 算例 | rank | 步数 | 终止时间 | 墙钟时间/s | 末步 Residual | 重构/Riemann fallback |
|---|---:|---:|---:|---:|---:|---:|
| Re=20 | 1 | 10114 | 15 | 1356.87 | 3.8322e-4 | 0 |
| Re=40 | 1 | 11106 | 20 | 1429.30 | 1.1276e-4 | 0 |
| Re=100 | 1 | 14767 | 30 | 1772.55 | 7.0846e-3 | 0 |
| Re=200 | 1 | 9614 | 20 | 912.41 | 2.7252e-2 | 0 |
| Mach 5 Euler（v1.0 全域零阶） | 1 | 10355 | 8 | 330.57 | 1.5625e-2 | 0 |
| Mach 5 robust（v1.1 高阶起步） | 4 | 12114 | 8 | 221.38 | 4.6916e-2 | 0 |

### 9.1 低速黏性圆柱绕流

以下数值是 v1.0 已归档运行的历史表，仍保留旧的第一层估算标签；用 v1.1 配置重新运行后，
后处理会以真实边界面输出重建同名报告，不能把下表作为 R 阶段新算法的量化证据。

| Re | `u_min`（尾迹） | 回流终点 `x/D` | `Cd_p` 历史估算 | `Cd_v` 历史估算 | `Cd` 历史估算 | 晚期 `Cl_rms` | `St(Cl)` | `St(probe)` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | -0.0263 | 1.158 | 1.348 | 0.807 | 2.156 | -- | -- | -- |
| 40 | -0.0877 | 1.761 | 1.134 | 0.513 | 1.647 | -- | -- | -- |
| 100 | -0.2063 | 2.310 | 1.000 | 0.269 | 1.268 | 0.0173 | 0.1403 | 0.1392 |
| 200 | -0.3311 | 2.640 | 0.989 | 0.160 | 1.149 | 0.0602 | 0.1479 | 0.1371 |

Re=20 和 Re=40 的上下剪切层与尾迹保持近似镜像对称，Re 增大时回流更强、回流区更长；
最终估计升力绝对值仅为 0.0029 和 0.0045。相对于第 10 节的不可压缩高阶参考值，Re=20 的
估计总阻力高 5.1%、回流终点短 17.3%，Re=40 的估计总阻力高 7.7%、回流终点短 33.0%。
方向和数量级正确，但径向只有 20 个单元且外边界仅为 8D，尾迹长度尚不能作为定量结果。

Re=100 和 Re=200 的升力与尾迹横向速度都出现一致的周期信号，流场上下对称性被打破；
`Cl_rms` 从 0.0173 增到 0.0602，说明更高 Re 下不稳定模态增长更明显。Re=100 的估计阻力
1.268 比参考均值 1.32 低 3.9%；两条信号给出的 `St=0.139--0.140`，比参考关系的 0.1666
低约 16%。目前 Re=100 信号仍在增长、Re=200 采样周期偏少，因此这里只能认定已经捕捉到
涡脱落，而不能认定极限环振幅和频率已经时间收敛。

对应图像：

- [`Re=20 涡量`](figures/re20-vorticity.png) 与 [`Re=40 涡量`](figures/re40-vorticity.png)；
- [`Re=100 涡量`](figures/re100-vorticity.png) 与 [`Re=200 涡量`](figures/re200-vorticity.png)；
- [`阻力、升力与尾迹探针历史`](figures/force-and-probe-history.png)；
- [`表面 Cp`](figures/surface-cp.png) 和 [`残差历史`](figures/residual-history.png)。

### 9.2 Mach 5 无粘钝体绕流

终场出现关于来流轴近似对称的脱体弓形激波、迎风驻点压缩区和背风低密度区：

| 指标 | 计算值 | 解释 |
|---|---:|---|
| `rho_min / rho_max` | 0.1250 / 4.2785 | 全场单元中心范围 |
| `p_min / p_max` | 0.01257 / 0.82076 | `p_max/p_inf=28.73` |
| `Ma_min / Ma_max` | 0.1081 / 5.0000 | 驻点附近显著减速 |
| 弓形激波迎风轴位置 | `x/D=-0.7094` | 由轴线上最大密度梯度估计 |
| 激波脱体距离 | `0.2094D` | 相对圆柱迎风点 `x/D=-0.5` |
| 历史压力阻力估算 | 1.2172 | v1.0 第一层单元归档值，仅供定性参考 |
| 历史升力估算 | 1.4e-5 | v1.0 归档值，用于观察上下对称性 |

计算的 `p_max/p_inf=28.73` 是正激波理论值 29 的 99.1%，而全场最大密度 4.279 是理论
密度比 5 的 85.6%。压力跃迁量级很好，密度和脱体距离则明显受 48x32 粗网格、零阶耗散和
粗网格与激波数值厚度影响。图像见 [`Mach 5 密度`](figures/mach5-density.png) 与
[`Mach 5 马赫数`](figures/mach5-mach.png)。

机器可读的全精度汇总见 [`results/analysis-summary.json`](results/analysis-summary.json)，每个
算例还保留末态 Tecplot 流场、完整残差/统计历史、运行 manifest、表面终值及抽样历史 CSV。
网格局部形状见 [`圆柱附近网格`](figures/grid-near-cylinder.png)。

### 9.3 v1.1 稳健化诊断

4-rank 历史共记录 1213 行抽样及终态。469 行出现局部降阶；单步最大低阶 owner 面数为
168/9360，即 `1.7949%`，没有出现全域零阶。抽样最大 level-1/2/3 面数分别为 99/36/156，
最大 troubled-cell 数 51，单步最多一次整步 retry；接受时间步最小抽样值为 `8.5791e-5`。
单步诊断的局部重算轮数最大为 16，这是三个 RK 阶段及一次 retry 的累计值；每个 RK 阶段仍
满足配置的最多 3 轮限制。
终态不再触发降阶，`rho_min=0.0188791`、`p_min=5.07856e-4`，全部真实单元保持物理容许。
人工复核图见 [`v1.1 密度`](figures/mach5-robust-density.png)、
[`v1.1 Mach`](figures/mach5-robust-mach.png) 和
[`降阶比例/troubled cell/接受时间步时序`](figures/mach5-robust-history.png)。history 没有保存
面坐标，因而不能从归档结果重建降阶面的空间位置；该空间图仍是 Q→R 人工卡口的待补证据。

开放远场会携带质量、动量和能量，不能把全域积分要求为常数；从 `t=0` 到 `t=8`，总质量和
总能量变化分别为 `+0.1353%` 和 `+0.1080%`，与远场通量方向一致。末态最大密度 4.8765、
最大压力 0.89130，对应 `p_max/p_inf=31.20`；相较 v1.0 全域零阶基线，更接近正激波密度比 5，
但该粗网格结果仍只作定性验收。

## 10. 参考尺度

高阶不可压缩参考计算给出：

| Re | `Cd_pressure` | `Cd_viscous` | `Cd_total` | 回流终点 `L/D` |
|---:|---:|---:|---:|---:|
| 20 | 1.229 | 0.823 | 2.052 | 1.40 |
| 40 | 0.994 | 0.536 | 1.530 | 2.63 |

这些数据使用约 `90 x 45` 网格、40D 外边界和不可压缩方程，不能要求本算例的 640 单元、8D 外
边界和 Ma=0.2 结果逐点一致，只用于判断数量级和随 Re 变化的方向。

Mach 5 迎风驻点附近可用正激波关系给出上限尺度：

\[
\frac{\rho_2}{\rho_1}=5,
\qquad
\frac{p_2}{p_1}=29,
\qquad
M_2\approx0.4152.
\]

公式来自 [NASA Glenn 正激波关系](https://www.grc.nasa.gov/www/k-12/airplane/normal.html)。圆柱弓形
激波在驻点流线上接近正激波，但第一层单元中心值和激波数值厚度不会精确达到理论跳跃端点。

## 11. 已知限制和下一步

若把 Case07 升级为定量验证，必须至少完成：

1. 用 48x32、96x48 以及更细网格做系统网格收敛；
2. 把外边界从 8--10D 推到 20--40D，检查低速尾迹长度和阻力的域敏感性；
3. 用 v1.1 真实边界面输出重新完成全部长时算例，并对升阻力做独立守恒复核；
4. Re=100 至少计算多个饱和周期，并同时由升力和尾迹探针验证 Strouhal 数；
5. 在更细网格上评估 v1.1 局部受控降阶的空间分布，并继续研究严格保正通量或缩放器；
6. Mach 5 使用更细的迎风区域网格，复核弓形激波脱体距离；
7. 若研究高雷诺数真实圆柱尾迹，需要三维网格、足够展向长度以及 DNS/LES/RANS 能力。

v1.0 的 Case07 建模没有修改求解器核心方程或边界实现；v1.1 阶段 Q 增加候选态验收、局部
面策略和时间步事务控制，阶段 R 增加只读边界工程量输出，均没有改变 Euler/NS 方程、几何
profile 或边界公式。所有失败、局部降阶和结果精度边界均在本文明确记录。
