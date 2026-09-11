# 人工验收 case02：三维 Poiseuille 通道流

本目录是 case02 的完整、可复现人工验收包。它在

$$
(x,y,z)\in[0,2\pi]\times[0,1]\times[0,\pi]
$$

上使用 $36\times48\times36$ 个三维结构单元，分别计算均匀网格和向上下壁面对称加密的
网格。$x,z$ 为平移周期方向，$y=0,1$ 为等温无滑移壁面。计算求解可压缩层流
Navier--Stokes 方程，用常量压力梯度驱动，以连续 Poiseuille 解析解作为初场，并且必须由
五个守恒分量的 MPI 全局残差判据返回 `steady_converged`；达到最大步数不算通过。

本例覆盖三维度量、i/k 双向周期连接、物理壁面 ghost、黏性梯度/热流、能量一致源项、
CGNS zone 少于 MPI rank 时的运行时剖分、稳态停止、三维流场重组和独立解析剖面检查。

## 1. 目录内容

```text
case02_3d_poiseuille/
|-- README.md                         本文
|-- uniform_36x48x36.wcns             均匀网格正式配置
|-- wall_clustered_36x48x36.wcns      壁面加密网格正式配置
|-- run_case02.py                     生成、运行、验证和汇总驱动
|-- case02-summary.json               机器可读的实测汇总（本次为强制停止失败状态）
|-- files.sha256                      生成文件 SHA-256 清单
|-- grids/
|   |-- uniform_36x48x36.cgns
|   `-- wall_clustered_36x48x36.cgns
|-- results/
|   |-- uniform/                      均匀网格初/终场、历史、统计和 manifest
|   `-- wall_clustered/               加密网格初场和强停时保留的历史/统计
|-- logs/                             网格生成及两次求解的完整输出
`-- validation/                       有限性与解析剖面检查报告
```

`grids/`、`results/`、`logs/`、`validation/`、`case02-summary.json` 和
`files.sha256` 都是验收产物；其中本次加密网格在第 3600 步收到人工强制停止，故其 history
和 statistics 以 `forced-step3600` 命名，且没有终场 CGNS/manifest。两份配置、驱动和本文
定义算例，不属于清理目标。

## 2. 三维多块网格

### 2.1 共同拓扑

两套网格均有 36×48×36 个单元和 37×49×37 个全局逻辑顶点。生成器在 $x,z$ 各分成
两个原生 CGNS zone，即 $2\times1\times2=4$ 个 zone；每个原生 zone 含
18×48×18 个单元。每个 zone 都有：

- `imin/imax` 两条 i 向连接，域首尾连接带 $2\pi$ 平移周期变换；
- `kmin/kmax` 两条 k 向连接，域首尾连接带 $\pi$ 平移周期变换；
- `bottom/top` 两个 y 向物理边界 patch。

正式计算使用 8 个 MPI rank。因为原生 zone 数 4 小于 rank 数 8，`auto_split` 会把每个
zone 沿 y 确定性二分，形成八个 18×24×18 叶块，rank 0--7 各持有一个。这既保留原生 i/k
周期通信，又实际验收运行时新增的兄弟连接和物理壁面切片。

### 2.2 均匀网格

顶点坐标为

$$
x_i=2\pi\frac{i}{36},\qquad
y_j=\frac{j}{48},\qquad
z_k=\pi\frac{k}{36}.
$$

对应间距为

| 方向 | 间距 |
|---|---:|
| $\Delta x$ | $\pi/18=0.17453292519943295$ |
| $\Delta y$ | $1/48=0.02083333333333333$ |
| $\Delta z$ | $\pi/36=0.08726646259971647$ |

### 2.3 壁面加密网格

$x,z$ 仍均匀，只把逻辑坐标 $s=j/48$ 通过对称 tanh 映射到 y：

$$
y(s)=\frac12\left[1+
\frac{\tanh\!\left(\beta(2s-1)\right)}{\tanh\beta}\right],
\qquad \beta=1.
$$

该映射严格保持 $y(0)=0$、$y(1/2)=1/2$、$y(1)=1$，上下半区镜像对称，坐标和一阶导数
连续。当前网格指标为

| 指标 | 数值 |
|---|---:|
| 壁面首层/最小 $\Delta y$ | 0.0118577795007742 |
| 中心附近/最大 $\Delta y$ | 0.0273390823838626 |
| 最大/最小间距比 | 2.30558194998293 |

`wall_cluster_strength=0` 就退化为均匀网格；增大 $\beta$ 会进一步减小壁面首层，但同时显著
收紧显式黏性时间步。本例选择 1.0 是“明确壁面加密”和“可承受稳态迭代代价”之间的基线，
不是湍流壁面 $y^+$ 设计。若改变强度，必须重新报告首层厚度、Jacobian 正性、稳定步长和
收敛步数。

## 3. 物理模型与解析初场

### 3.1 无量纲参数

本例采用

$$
\gamma=1.4,\qquad Re=10,\qquad Ma=0.2,\qquad Pr=0.72.
$$

五个参考量为 $U_{ref}=\rho_{ref}=T_{ref}=L_{ref}=1$、$\mu_{ref}=0.1$，所以
$Re=\rho_{ref}U_{ref}L_{ref}/\mu_{ref}=10$。令比气体常数

$$R=\frac{1}{\gamma Ma^2}=17.857142857142858,$$

程序由参考量得到 $Ma=1/\sqrt{\gamma R}=0.2$。压力按动压标度无量纲化，基准常压取

$$p_0=\frac{1}{\gamma Ma^2}=17.857142857142858,$$

从而当 $T=1$ 时由 $\rho=\gamma Ma^2p/T$ 得到 $\rho=1$。当前生产配置尚未开放
`transport.prandtl`，黏性路径使用程序内已验收默认值 $Pr=0.72$ 和常黏度
$\mu/\mu_{ref}=1$；修改 Pr 需要先扩展配置协议，不能只改本文公式。

### 3.2 压力梯度驱动及能量闭合

周期方向不能直接在首尾指定不同压力。本例把常压降写成单位体积驱动力

$$
\boldsymbol G=-\nabla p=(0.8,0,0).
$$

`pressure_gradient` 模型在每个 SSPRK3 子步加入

$$
\mathcal S=(0,G_x,G_y,G_z,\boldsymbol u\cdot\boldsymbol G)^T.
$$

最后一项是压力功，不能省略；只向动量方程加 0.8 会破坏总能量守恒形式。它与“单位质量
体力”不同：`body_force` 的动量源为 $\rho\boldsymbol a$，密度随黏性升温变化后不能严格代表
常压力梯度。这里使用常体积力使解析抛物线在可压缩变密度状态下仍满足 x 动量平衡。

### 3.3 速度、温度、密度解析式

令 $\eta=(y-y_0)/(y_1-y_0)=y$、中心线速度 $U_c=1$。速度取

$$
u(\eta)=4U_c\eta(1-\eta),\qquad v=w=0.
$$

因此 $u''=-8U_c$，常黏度 x 动量方程满足

$$
\frac{1}{Re}u''+G_x=-\frac8{10}+0.8=0.
$$

等温壁温 $T_w=1$。动量方程与压力功相消后，总能量方程要求

$$
T''=-(\gamma-1)Ma^2Pr\,(u')^2.
$$

定义

$$
P(\eta)=\frac{\eta}{6}-\frac{\eta^2}{2}
 +\frac{2\eta^3}{3}-\frac{\eta^4}{3},
\qquad P''(\eta)=-(1-2\eta)^2,
$$

则热平衡初场为

$$
T(\eta)=T_w+A_TP(\eta),
$$

其中

$$
A_T=(\gamma-1)Ma^2Pr(4U_c)^2=0.18432.
$$

中心温度为 $T(1/2)=1+A_T/48=1.00384$。最后取

$$p=p_0,\qquad \rho=\frac{\gamma Ma^2p_0}{T}=\frac1T.$$

初始化器在真实单元的物理中心逐点计算以上原始量，再统一转换成五个守恒量；ghost 仍只由
周期连接或壁面边界算子填充。连续解析解与有限宽度物理边界闭合的离散定常解会有小差别，
所以必须实际迭代到残差收敛，并在终场独立测量剖面误差，而不能因为使用解析初场就跳过推进。

## 4. 配置文件逐项说明

两份 `.wcns` 仅在 `case.name`、`mesh.path` 和 `output.directory` 不同，数值和物理参数一致。
配置是严格 UTF-8 `key = value`：未知键、重复键、非法枚举、NaN/Inf 会在分配流场前失败。

### 4.1 算例、网格和算法

```text
schema_version = 1
case.name = case02-poiseuille-uniform-36x48x36
mesh.path = grids/uniform_36x48x36.cgns
```

- `schema_version` 固定当前协议版本 1。
- `case.name` 成为输出文件名前缀；加密网格配置使用 `case02-poiseuille-wall-clustered-36x48x36`。
- `mesh.path` 相对配置文件所在目录解析；加密配置指向对应 CGNS。

```text
algorithm.profile = phenglei_wcns
algorithm.reconstruction = mdcd_linear
algorithm.reconstruction_variables = primitive
algorithm.riemann = hllc
```

- `profile` 选择 PHengLEI-WCNS 成套度量与差分闭合。可整体换成 `scmm6_wcns`，但不得交叉
  混用两套 profile 的度量和差分部件。
- Poiseuille 解全域光滑，所以选择低色散、可控耗散的 `mdcd_linear`。强间断问题不得照搬；
  可替换为 `weno_js`、`weno_z` 或 `mdcd_hybrid`。
- `primitive` 表示分别重构原始变量。可换 `conservative` 或 `characteristic`；后者对激波更
  合适，但本光滑三维稳态例会增加明显开销。
- `hllc` 计算无粘面通量。可换 `rusanov` 或 `roe`，任何替换均须重新收敛和校验。
- 生产时间推进固定 SSPRK3；当前没有可在非定常/稳态间任意混用的低 Mach 预处理键。

### 4.2 气体和参考量

```text
gas.gamma = 1.4
gas.specific_gas_constant = 17.857142857142858
reference.velocity = 1.0
reference.density = 1.0
reference.temperature = 1.0
reference.length = 1.0
reference.viscosity = 0.1
```

- `specific_gas_constant` 与 `gas.molar_mass` 严格二选一；本例直接给 R 以精确导出 Ma=0.2。
- 五个参考量分别控制速度、密度、温度、长度和黏度标度。
- Re 和 Ma 禁止直接输入，只在日志/manifest 中报告。
- 若要改 Re，至少同时修改 `reference.viscosity` 和压力梯度
  $G_x=8U_c/Re$；只改一项就不再对应同一解析解。
- 若要改 Ma，须通过 R 或其他参考量导出，并同步更新 $p_0=1/(\gamma Ma^2)$ 和
  $A_T=(\gamma-1)Ma^2Pr(4U_c)^2$。

### 4.3 分区

```text
partition.mode = auto_split
partition.allow_idle_ranks = false
partition.max_load_ratio = 1.2
partition.min_cells_per_active_direction = 8
```

- `auto_split` 允许 4 个原生 zone 继续切成 8 个叶块。
- `allow_idle_ranks=false` 要求八个 rank 都参与；不可行时启动失败而不静默空闲。
- `max_load_ratio=1.2` 是最大负载/平均负载目标。
- `min_cells_per_active_direction=8` 防止切出无法容纳高阶模板的窄块；最终叶块三方向为
  18×24×18，满足约束。

`--ranks` 可改成 4，使每个 rank 对应一个原生 zone；也可在可行范围内改成其他数量。更改
rank 会改变分区与性能，应先做 MPI `--dry-run`，并把结果视为新的并行验收记录。

### 4.4 解析初场

```text
initial.type = poiseuille
initial.y0 = 0.0
initial.y1 = 1.0
initial.centerline_velocity = 1.0
initial.temperature = 1.0
initial.temperature_curvature = 0.18432
initial.pressure = 17.857142857142858
```

- `poiseuille` 选择第 3.3 节的解析初始化器。
- `y0/y1` 定义两壁坐标和归一化 $\eta$。
- `centerline_velocity` 是 $U_c$，不是抛物线前的系数；代码自动使用 $4U_c$。
- `temperature` 是两侧相同壁温，`temperature_curvature` 是 $P(\eta)$ 的系数 $A_T$。
- `pressure` 是全场常压 $p_0$。

改变 $U_c$ 时必须同步改 $G_x=8U_c/Re$ 和
$A_T=(\gamma-1)Ma^2Pr(4U_c)^2$。所有温度、压力和密度必须保持正值。

### 4.5 壁面和周期边界

```text
boundary.default = no_slip_isothermal_wall
boundary.bottom.type = no_slip_isothermal_wall
boundary.bottom.wall_velocity_x = 0.0
boundary.bottom.wall_velocity_y = 0.0
boundary.bottom.wall_velocity_z = 0.0
boundary.bottom.wall_temperature = 1.0
boundary.top.type = no_slip_isothermal_wall
...
```

- `bottom/top` 与 CGNS patch 名严格对应，显式设静止、等温、无滑移壁。
- `default` 为遗漏 patch 提供同类保护，但不会把 CGNS 连接误当物理壁面。
- x/z 周期性由 CGNS `1to1` connectivity 和周期平移定义，不能再写
  `boundary.left.type=periodic`；周期边界必须是连接而不是物理 patch。
- 若改成绝热壁，热解析式和终场温度都会改变，不能继续使用本例温度设置与基线。

物理壁面 ghost 只保证由边界条件得到的原始/守恒物理量有效；算法不读取壁面 ghost 坐标、
度量、梯度或二级量。边和角 ghost 也不作为当前面模板的输入。

### 4.6 压力梯度源

```text
source.enabled = true
source.models = pressure_gradient
source.pressure_gradient.x = 0.8
source.pressure_gradient.y = 0.0
source.pressure_gradient.z = 0.0
```

- 开关打开源项装配；模型列表不能为空。
- 三个分量保存 $\boldsymbol G=-\nabla p$，不是 $\nabla p$ 本身，因此正 x 流动使用正 0.8。
- 每个 RK 子步都用当前局部速度计算能量功 $\boldsymbol u\cdot\boldsymbol G$。
- 若改用 `body_force`，键应改成 `source.body.ax/ay/az`，其动量源是
  $\rho\boldsymbol a$，物理问题随之改变。

### 4.7 稳态迭代和停止判据

```text
run.mode = steady
run.viscous = true
run.cfl = 0.5
run.max_steps = 20000
run.max_wall_time = 0
```

- `steady` 中累计的 `time` 是伪时间，不是物理终止时刻。
- `viscous=true` 启用 Navier--Stokes 黏性应力与 Fourier 导热；关闭后本例不存在稳态平衡。
- `cfl=0.5` 控制含声学和黏性限制的全局显式时间步。该值经均匀网格预运行确认残差稳定
  衰减；若继续增大 CFL，必须重新检查两套网格的稳定性，不能只凭均匀网格判断。
- `max_steps` 是硬保护上限，不是成功条件。
- `max_wall_time=0` 关闭墙钟停止；若设正值，必须同时启用检查点。

```text
steady.min_steps = 20
steady.check_interval_steps = 20
steady.consecutive_checks = 3
steady.reference_floor = 1e-30
steady.l2_absolute = 5e-5
steady.l2_relative = 1e-6
steady.linf_enabled = true
steady.linf_absolute = 2e-4
steady.linf_relative = 1e-5
```

- 每 20 步进行一次判定；第一次检查冻结五个守恒分量各自的参考 L2/Linf。
- 对每个分量，L2 满足“绝对值不大于 $5\times10^{-5}$”或“相对参考不大于
  $10^{-6}$”之一；Linf 同理采用 $2\times10^{-4}$/$10^{-5}$。
- 五个分量必须同时通过 L2 和 Linf，并连续通过 3 次，才返回 `steady_converged`。
- `min_steps=20` 防止初始残差异常小导致零步成功。
- 绝对阈值按本例的无量纲方程尺度确定：相对驱动源 $G_x=0.8$，上述 L2/Linf 上限分别为
  $6.25\times10^{-5}$ 和 $2.5\times10^{-4}$；在正式时间步约 $O(10^{-4})$ 时，Linf 上限
  对应单步守恒量改变量约 $O(10^{-8})$。解析初场使冻结参考残差本来就很小，因此不采用
  “相对该小量再下降很多数量级”作为唯一尺度。若返回 `maximum_steps`、
  `numerical_failure`、墙钟或信号停止，算例仍都不通过。

### 4.8 流场输出

```text
output.directory = results/uniform
output.allow_existing = false
output.dimensional = false
output.field.enabled = true
output.field.format = cgns
output.field.every_steps = 0
output.field.every_time = 0
output.field.write_initial = true
output.field.write_final = true
output.field.quantities = rho,u,v,w,p,T,rho_u,rho_v,rho_w,rho_E,sound_speed,mach,total_enthalpy,entropy_proxy,viscosity,jacobian
```

- 两套配置分别写 `results/uniform` 和 `results/wall_clustered`。输出目录相对进程当前目录解析，
  驱动会先切换到 case 目录。
- `allow_existing=false` 防止覆盖；重复正式运行应使用驱动 `--clean` 或新的输出目录。
- `dimensional=false` 保存内部无量纲值。
- `format=cgns` 按原四个 zone 重组三维 `CellCenter` 解。可改为 `tecplot` 或 `both`；Tecplot
  ASCII 便于直接作图，但不是检查点。
- 中间场调度均为 0，所以只保存初场和收敛终场。长算例若需观察过程，可设
  `every_steps=200`；这只改变输出，不应改变数值轨迹。
- 16 个字段包括原始量、守恒量、声速/Mach/总焓/熵代理、常黏度和 Jacobian。可删减；未知
  字段会明确失败。

### 4.9 残差、统计和检查点

```text
output.history.enabled = true
output.history.format = txt
output.history.every_steps = 20
output.history.write_initial = true
output.history.write_final = true
```

历史与收敛检查同频。固定列含 step、伪时间、dt、CFL、墙钟、总 L2、五分量 L2/Linf、冻结
参考值、归一化值、连续通过次数、重构/Riemann 回退计数、判定标志和停止原因。

```text
output.statistics.enabled = true
output.statistics.format = txt
output.statistics.every_steps = 20
output.statistics.write_initial = true
output.statistics.write_final = true
output.statistics.quantities = total_mass,total_momentum_x,total_momentum_y,total_momentum_z,total_energy
output.statistics.xz_planes.enabled = true
output.statistics.xz_planes.cell_j_indices = 0,11,23,35,47
```

统计量使用正 Jacobian 和原 zone 守恒权重作 MPI 全局积分。压力梯度持续向流体输入动量和功，
因此不能要求总 x 动量/总能量在伪时间中保持常数；y/z 总动量应保持近零。

新增的 x-z 平面监测以原始 CGNS zone 的零基单元 J 索引选层；以上五层覆盖近下壁、四分之一、
中心附近、四分之三和近上壁位置。每个索引自动向 statistics 追加两列：
`xz_mean_u_jN` 为面积加权流向平均速度，`xz_mass_flow_x_jN` 为
$\int_{xz}\rho u\,\mathrm dA$。后者是在壁平行截面上对 x 向质量通量密度作面积积分，用于监测
不同 y 层的一致性，不应解释为穿过该 x-z 面的法向流量。开关关闭时不注册这些列；可修改
索引列表来监测更多层，但各原 zone 必须具有相同合法 J 索引，且选中层必须几何共面。

本目录已归档的长时运行结果生成于该监测功能加入之前，现有 `*.statistics*.txt` 因而只含
原来的五个全场总量列。两份配置现已升级以便下一次重算直接产生十个新增截面列；不能用旧
结果文件声称已经完成截面统计验收。

```text
output.checkpoint.enabled = false
```

本次连续验收不写检查点。若实际运行环境需要墙钟分段，启用检查点后可由不同合法 rank 数
重启，但必须保留相同 profile、源项、物性、壁面和网格签名。

## 5. 从构建到验收的逐步命令

除特别说明外，命令从仓库 `wcns/` 根目录运行。

### 步骤 1：配置并构建 Release MPI 程序

已有构建目录时：

```powershell
cmake -S . -B build-rc-mpi
cmake --build build-rc-mpi --parallel 4
```

第一条重新生成构建系统并把当前 Git 提交写入版本 manifest；第二条编译求解器、网格生成器
和校验器。从空目录开始时应明确启用 MPI/CGNS：

```powershell
cmake -S . -B build-rc-mpi -G "MinGW Makefiles" `
  -DCMAKE_BUILD_TYPE=Release -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_MPI=ON
cmake --build build-rc-mpi --parallel 4
```

### 步骤 2：运行基础自动测试

```powershell
ctest --test-dir build-rc-mpi -R "wcns\.(unit|cgns_reader)" --output-on-failure
```

它检查解析初始化/压力源单元契约以及 CGNS 读取和周期连接既有回归。通过只说明基础组件没有
回归，不能替代本 36×48×36 正式计算。

### 步骤 3：一键执行完整 case02

```powershell
python -B cases/manual/case02_3d_poiseuille/run_case02.py --clean
```

`-B` 禁止 Python 字节码缓存；`--clean` 只删除本 case 内的生成目录、汇总和散列。脚本依次：

1. 生成均匀和 $\beta=1$ 壁面加密 CGNS；
2. 使用 Intel MPI 8 rank 计算均匀网格，直到程序报告残差收敛；
3. 校验其初/终场有限性，严格核对初始解析剖面，再按冻结容差核对终场；
4. 对壁面加密网格重复计算和校验；
5. 仅在全部成功后写 `case02-summary.json` 和 `files.sha256`。

替换构建目录或 MPI 启动器：

```powershell
python -B cases/manual/case02_3d_poiseuille/run_case02.py --clean `
  --run build-other/wcns_run.exe `
  --generator build-other/wcns_generate_release_cgns.exe `
  --validator build-other/wcns_validate_release_case.exe `
  --mpi-exec "C:/path/to/mpiexec.exe" --ranks 8
```

脚本还接受终场 `--velocity-l2-tolerance`、`--crossflow-tolerance`、
`--pressure-span-tolerance`、`--homogeneity-tolerance`，用于更严格的专项试验。不得仅为让失败
结果通过而放宽；任何更改必须在新结果文档中给出物理理由。

若求解已经以 `steady_converged` 完成、但后处理或校验阶段中断，可续跑：

```powershell
python -B cases/manual/case02_3d_poiseuille/run_case02.py --resume
```

`--resume` 与 `--clean` 互斥。它不会仅凭目录存在就跳过计算：脚本要求恰好一个与当前 rank
数相符的 manifest、要求其中 `stop_reason=steady_converged`，并重新执行 MPI `--dry-run`，把
当前配置 digest 与 manifest 中的 digest 比较。只有完全一致才复用该终场；随后所有有限性和
剖面检查仍会重跑，缺失的另一套网格仍正常求解。续跑无法恢复此前驱动进程测得的峰值 RSS，
汇总中该项记为 `null`，但求解器 manifest 的实际墙钟时间会保留。

### 步骤 4（可选）：单独生成网格

从 case 目录运行：

```powershell
cd cases/manual/case02_3d_poiseuille
../../../build-rc-mpi/wcns_generate_release_cgns.exe periodic-channel `
  grids/uniform_36x48x36.cgns 36 48 36 2 2 `
  6.283185307179586 1.0 3.141592653589793 0.0
../../../build-rc-mpi/wcns_generate_release_cgns.exe periodic-channel `
  grids/wall_clustered_36x48x36.cgns 36 48 36 2 2 `
  6.283185307179586 1.0 3.141592653589793 1.0
```

参数依次为输出文件、i/j/k 总单元数、i/k 原生 zone 数、x/y/z 长度和壁面加密强度。

### 步骤 5（可选）：MPI 启动检查

```powershell
& "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe" -n 8 `
  ../../../build-rc-mpi/wcns_run.exe --config uniform_36x48x36.wcns --dry-run
& "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe" -n 8 `
  ../../../build-rc-mpi/wcns_run.exe --config wall_clustered_36x48x36.wcns --dry-run
```

`--dry-run` 读取/广播配置和 CGNS，完成八叶块剖分、周期拓扑、度量、初场、边界和输出注册
校验，但不推进或生成正式结果。摘要必须显示 8 个 18×24×18 叶块、Re=10、Ma=0.2。

### 步骤 6（可选）：拆开正式求解

```powershell
& "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe" -n 8 `
  ../../../build-rc-mpi/wcns_run.exe --config uniform_36x48x36.wcns
& "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe" -n 8 `
  ../../../build-rc-mpi/wcns_run.exe --config wall_clustered_36x48x36.wcns
```

两个命令应顺序执行，避免资源争用。运行中每步在日志打印残差；每 20 步的完整分量数据写
history。最终退出码必须为 0，末行必须含 `reason=steady_converged`。

### 步骤 7：独立检查结果

将 `<field.cgns>` 换成实际初/终场名：

```powershell
../../../build-rc-mpi/wcns_validate_release_case.exe finite <field.cgns>
../../../build-rc-mpi/wcns_validate_release_case.exe poiseuille-profile `
  <field.cgns> 0.0 1.0 1.0 5e-3 1e-6 1e-2 1e-7
```

`finite` 遍历全部 16 个场并拒绝 NaN/Inf 或非正 $\rho,p,T$。`poiseuille-profile` 使用
CGNS 单元中心和 Jacobian 独立计算：

- $u=4y(1-y)$ 的体积加权 L2/Linf 与体积平均速度；
- $v,w$ 最大绝对值；
- 压力最大值减最小值；
- 相同 y 层上 $\rho,u,v,w,p,T$ 的最大 x/z 非均匀度；
- $\rho,T$ 实际范围。

初场由脚本使用 $10^{-12}$ 速度/横流和 $10^{-10}$ 压力/均匀度容差严格检查；终场采用命令中
较适合离散稳态解的冻结容差。

### 步骤 8：核对散列并返回仓库根

```powershell
Get-Content files.sha256 | ForEach-Object {
  $expected, $relative = $_ -split '  ', 2
  $actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $relative).Hash.ToLower()
  if ($actual -ne $expected) { throw "checksum mismatch: $relative" }
}
cd ../../..
```

散列用于发现生成后文件被截断或修改，不替代收敛和物理剖面检查。

## 6. 正式运行结果

本节由同一次 8-rank Release 正式验收填写。机器原始记录见 `case02-summary.json`，最终步数
和停止原因以 manifest 为准。

<!-- CASE02_RESULTS_BEGIN -->

### 6.1 验收结论

本次整体结论为 **未完全通过**，不是程序崩溃：

- 均匀网格在第 680 步由五分量残差连续三次通过，`stop_reason=steady_converged`，有限性和
  解析剖面检查均通过；
- 壁面加密网格按人工指令在观察到第 3600 步 history 开始写入后强制停止。最后一条完整
  history 是第 3580 步，其中全部 L∞ 通过，但 x 动量和总能量 L2 尚未通过，所以该网格不能
  标记为稳态；
- 强制信号恰在第 3600 步 history/statistics 行写入期间到达，两个末行被截断，求解器也没有
  提交终场 CGNS 或 manifest。因此第 3600 步只能证明“推进到该步后被外部终止”，场级分析
  必须止于已有初场，残差分析以第 3580 步最后完整行为准。

机器可读结论见 `case02-summary.json`，明确保存
`status=forced_stop_not_fully_accepted` 和 `overall_acceptance=false`。本次可执行程序内嵌 Git
提交为 `51b2820fb39d`，两套网格均使用 8 个 MPI rank 和八个 18×24×18 运行时叶块。

### 6.2 均匀网格：已通过

| 项目 | 实测值 |
|---|---:|
| 网格签名 | `11936517813990972397` |
| 配置 digest | `0xe6b998156e40704e` |
| 停止原因 | `steady_converged` |
| 最终步数 | 680 |
| 最终伪时间 | 0.19897579806505086 |
| 最终 dt | 0.00029261147527331785 |
| 求解器墙钟 | 2233.353254 s（约 37 min 13.35 s） |
| 连续通过次数 | 3 |
| 重构/Riemann 回退 | 0 / 0 |

最终完整残差为：

| 守恒分量 | L2 | L∞ | L2 是否通过 | L∞ 是否通过 |
|---|---:|---:|:---:|:---:|
| $\rho$ | 1.06520960213937e-6 | 1.98752687016462e-6 | 是 | 是 |
| $\rho u$ | 3.14585972952092e-5 | 5.16239758818760e-5 | 是 | 是 |
| $\rho v$ | 4.71621522681942e-6 | 6.82638159434487e-6 | 是 | 是 |
| $\rho w$ | 6.45691839961769e-14 | 2.91432161311462e-13 | 是 | 是 |
| $\rho E$ | 3.13394504996680e-5 | 5.95590508346011e-5 | 是 | 是 |

总 L2 为 1.99758947682724e-5。最后三次完整检查为第 640、660、680 步，`consecutive` 依次
为 1、2、3，最终停止原因由 `running` 变为 `steady_converged`。

终场 CGNS 的独立检查得到：

| 量 | 实测值 | 容差/参考 | 结论 |
|---|---:|---:|:---:|
| $u=4y(1-y)$ 体积加权 L2 | 2.64892987320790e-5 | 5e-3 | 通过 |
| 速度 Linf | 5.65770431804336e-5 | 仅报告 | — |
| 体积平均速度 | 0.666807386476045 | 精确体积平均 2/3 | 偏差 1.40720e-4 |
| $\max(|v|,|w|)$ | 2.10657267401372e-7 | 1e-6 | 通过 |
| 压力跨度 | 1.95128111002418e-6 | 1e-2 | 通过 |
| 同 y 层 x/z 最大非均匀度 | 4.61852778244065e-14 | 1e-7 | 通过 |
| 密度范围 | [0.996174625796517, 0.999689308296897] | 正值 | 通过 |
| 温度范围 | [1.00031056773286, 1.00383995215693] | 正值 | 通过 |

有限性校验实际遍历 995328 个“单元×字段”样本，最小密度、压力、温度分别为
0.996174625796517、17.8571389203382、1.00031056773286。最终场文件为
`results/uniform/case02-poiseuille-uniform-36x48x36.field.step00000680.time1p989757981eM01.cgns`。

积分量从初场到终场的变化为：总质量 -1.64135e-12、x 动量 -7.81063e-5、y 动量
-2.89476e-17、z 动量 -1.19043e-15、总能量 -1.96861e-4。质量相对漂移约 8.34e-14，横向
动量保持机器精度；x 动量和能量的小调整来自连续解析初场向离散定常解的松弛，不能按封闭
无源系统要求严格守恒，因为存在压力梯度功和壁面黏性/热通量。

### 6.3 壁面加密网格：第 3600 步强制停止，未通过

| 项目 | 实测值 |
|---|---:|
| 网格签名 | `15472241633829854185` |
| 配置 digest | `0x63d999e3986e2636` |
| 人工要求停止步 | 3600 |
| history 已观察但不完整的末步 | 3600 |
| 最后完整 history 步 | 3580 |
| 最后完整 statistics 步 | 3560 |
| 第 3580 步伪时间 | 0.39333542008419314 |
| 第 3580 步 dt | 0.00010987023617565823 |
| 第 3580 步求解器墙钟 | 11244.79784 s（约 3 h 7 min 24.80 s） |
| 第 3600 步不完整行记录的墙钟 | 11327.235926 s（约 3 h 8 min 47.24 s） |
| 第 3580 步停止状态/连续计数 | `running` / 0 |
| 重构/Riemann 回退 | 0 / 0 |
| 终场 CGNS / manifest | 无 / 无 |

第 3580 步最后完整残差为：

| 守恒分量 | L2 | L2/绝对上限 | L∞ | 结论 |
|---|---:|---:|---:|:---:|
| $\rho$ | 1.55699369816747e-6 | 0.0311 | 3.03824910282819e-6 | 通过 |
| $\rho u$ | 7.56068810288939e-5 | 1.5121 | 1.31578677458499e-4 | **L2 未通过** |
| $\rho v$ | 1.15034691046300e-6 | 0.0230 | 1.64365653527838e-6 | 通过 |
| $\rho w$ | 8.25557054072345e-14 | 1.65e-9 | 3.85198687561228e-13 | 通过 |
| $\rho E$ | 1.16625822440475e-4 | 2.3325 | 1.81426238803462e-4 | **L2 未通过** |

五个 L∞ 都已小于 2e-4，但 x 动量和能量 L2 分别高于 5e-5 上限 51.2% 和 133.3%；总 L2
为 6.21638648190813e-5。相对于第 20 步冻结参考，x 动量 L2 降到 0.0084768 倍，能量降到
0.0277222 倍；它们确实在收敛，但尚未满足本例定义的稳态。

若只看若干代表步，衰减过程为：

| step | 伪时间 | 总 L2 | $\rho u$ L2 | $\rho E$ L2 | $\rho u$ L∞ | $\rho E$ L∞ |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 2.167679e-2 | 3.841066e-2 | 2.956413e-2 | 2.133166e-1 | 1.678785e-1 |
| 20 | 2.197403e-3 | 4.410458e-3 | 8.919236e-3 | 4.206946e-3 | 4.326242e-2 | 1.391973e-2 |
| 500 | 5.493510e-2 | 5.945263e-4 | 1.225717e-3 | 5.135896e-4 | 2.752387e-3 | 1.119419e-3 |
| 1000 | 1.098702e-1 | 3.335728e-4 | 6.710727e-4 | 3.238880e-4 | 1.341741e-3 | 6.793935e-4 |
| 1500 | 1.648053e-1 | 2.836107e-4 | 4.078350e-4 | 4.854845e-4 | 7.871173e-4 | 7.676101e-4 |
| 2000 | 2.197405e-1 | 1.449293e-4 | 2.524329e-4 | 2.023048e-4 | 4.813349e-4 | 3.397260e-4 |
| 2500 | 2.746756e-1 | 1.019449e-4 | 1.597139e-4 | 1.625605e-4 | 3.048382e-4 | 2.755329e-4 |
| 3000 | 3.296107e-1 | 5.549975e-5 | 1.082529e-4 | 5.934880e-5 | 2.017658e-4 | 1.193406e-4 |
| 3500 | 3.845458e-1 | 6.532815e-5 | 7.905883e-5 | 1.228170e-4 | 1.393617e-4 | 1.917417e-4 |
| 3580 | 3.933354e-1 | 6.216386e-5 | 7.560688e-5 | 1.166258e-4 | 1.315787e-4 | 1.814262e-4 |

表中 3000--3580 步总 L2 的非单调性来自能量阻尼振荡；x 动量残差仍单调下降。该事实再次
说明必须逐分量并连续判断，不能因为第 3000 步总 L2 接近阈值就宣称收敛。

加密初场本身通过严格解析检查：速度误差和横流为 0，压力跨度 3.55e-15，层内非均匀度 0，
体积平均速度为 0.666857162081033。初场有限性检查的最小密度、压力、温度为
0.996174691322312、17.8571428571429、1.00017892144295。但由于没有终场 CGNS，无法计算
第 3600 步的解析速度误差、横流、压力跨度、x/z 均匀度或热力学范围；这些量必须标为“无
数据”，不能用初场值或第 3580 步残差替代。

第 3560 步最后完整积分记录表明：质量相对初场改变 -4.97735e-12（约 -2.53e-13 相对量），
x 动量改变 -1.80053e-3，总能量改变 -8.18792e-4，y/z 总动量仍为 1.62e-18/-1.62e-15。
质量与横向对称性保持良好；主流和能量的变化与残差仍未完全收敛的结论一致。

### 6.4 两套网格的对比和本次强停暴露的问题

- 加密网格的最小 $\Delta y$ 是均匀网格的 0.5692 倍，其实际 dt 是均匀网格的 0.3755 倍；
  因此达到相同伪时间至少需要约 2.66 倍步数。
- 加密网格连续解析初场在离散黏性算子上的初始总 L2 为 2.16768e-2，是均匀网格
  3.09786e-3 的约 7.00 倍；x 动量初始 L2 约为均匀网格的 13.17 倍。这进一步增加了显式
  松弛代价。
- 加密网格运行到 3580 步的墙钟已约为均匀网格完整收敛的 5.04 倍，却仍未稳态。这不是
  MPI 分区失败：两者同为八个等单元数叶块，主要限制来自最小网格尺度、非均匀网格离散
  失配和全局显式时间步。
- 当前程序未启用检查点。外部 Ctrl+C 在 MPI/输出写入期间属于硬中断，本次没有获得终场，
  同时留下被截断的第 3600 步 history 和第 3580 步 statistics 末行。若预先知道要在 3600
  步取场，应在启动前令 `run.max_steps=3600` 并使用新的输出目录，让求解器以
  `maximum_steps` 正常走完 `on_final`；它仍属于“未收敛失败”，但会有终场 CGNS/manifest。
- 若需从人工停止点续算，必须在停止前启用 checkpoint。当前文件不能恢复 3600 步流场；
  `--resume` 只会复用带 `steady_converged` manifest 且配置 digest 一致的完整结果，故会拒绝
  本次加密目录。要重新继续，需先归档/移走该部分结果，再用 `--resume` 复用均匀网格并从
  加密初场重新计算。

<!-- CASE02_RESULTS_END -->

## 7. 通过条件

以下是 case02 的目标通过条件；本次只有均匀网格满足，加密网格因第 3600 步人工强停未满足
第 1、2、3、4、6 项，所以整体状态必须保持失败：

1. 两套求解均退出码 0，manifest 停止原因为 `steady_converged`，不是最大步数；
2. history 最后三次检查中五分量 L2/Linf 均按第 4.7 节规则连续通过；
3. 两个初场和两个终场共四个 CGNS 中所有请求字段有限，$\rho,p,T>0$；
4. 初始解析剖面通过严格检查；终场速度 L2 不超过 $5\times10^{-3}$，横流不超过
   $10^{-6}$，压力跨度不超过 $10^{-2}$，同 y 层 x/z 非均匀度不超过 $10^{-7}$；横流上限
   相对 $U_c=1$ 为百万分之一，用于容纳残差卡口内的微小阻尼声学扰动；
5. 重构/Riemann 回退无异常增长，日志无网格、周期连接、壁面、MPI、黏性或 I/O 错误；
6. 输出完整且 SHA-256 清单逐项复算一致。

解析剖面通过不等于网格收敛证明。均匀与壁面加密网格的 y 单元中心不同，不能逐数组索引
相减；若要比较两套解，应按公共物理 y 位置插值，或另做系统网格加密研究。

## 8. 每类输出文件包含什么

本次目录共 28 个文件、33434881 字节。`files.sha256` 覆盖其中 23 个生成/验收产物；不包含
本文、两份配置、驱动脚本和散列文件自身。

- `*.field.step........time....cgns`：四个原生 zone 上的三维 `CellCenter` 场，含配置中 16 个
  无量纲字段及原网格坐标。本次均匀网格有初/终场，加密网格只有初场。
- `*.history.r8.txt`：正常提交的步数、伪时间、dt、CFL、墙钟、五分量 L2/Linf、参考/归一化
  残差、连续通过次数、回退计数、判定标志与停止原因。本次加密网格的原始临时文件重命名为
  `*.history.r8.forced-step3600.txt`，第 3600 步末行不完整，最后完整行为 3580。
- `*.statistics.r8.txt`：旧归档文件记录总质量、三个方向总动量和总能量随伪时间积分；它们
  早于 x-z 截面监测功能，不含 `xz_mean_u_jN`/`xz_mass_flow_x_jN`。按当前配置重算时会追加
  五个 J 层各两列。本次加密文件名为 `*.statistics.r8.forced-step3600.txt`，最后完整行为 3560。
- `*.manifest.r8.txt`：程序版本、Git 提交、编译器/构建类型、MPI 数、配置/分区摘要、网格与
  重启签名、最终步/伪时间/dt/墙钟、停止原因和成功提交文件清单。
- `logs/generate-*.log`：两套网格生成命令结果。
- `logs/run-*.log`：完整 rank-0 配置/分区摘要及每步残差，定位长稳态运行问题的首要证据。
- `validation/*-finite.txt`：场样本总数和 $\rho,p,T$ 下界。
- `validation/*-profile.txt`：解析速度误差、体积平均速度、横流、压力跨度、层内均匀度及
  $\rho,T$ 范围。
- `validation/wall-clustered-forced-stop.txt`：人工停止点、最后完整步、残差阈值比值、回退计数
  以及没有终场 CGNS/manifest 的原因。
- `case02-summary.json`：域、网格间距、物理/算法参数、容差、均匀网格完整结果，以及加密网格
  的强停状态、最后完整残差/统计和缺失终场标志。
- `files.sha256`：所有可重新生成文件的 SHA-256；不包含配置、驱动、本文和散列文件自身。

## 9. 合理替换与常见错误

- **网格**：可改分辨率和加密强度，但 2×2 周期 zone 要求 i/k 单元数分别可被 2 整除；每个
  运行时叶块活动方向还必须满足模板最小宽度。
- **物理参数**：Re、Ma、$U_c$、Pr 任一改变都要按第 4.2/4.4 节联动更新黏度、R、压力梯度、
  基准压力和热曲率，不能孤立修改。
- **算法**：WENO-Z/特征重构可作为更昂贵交叉检查；`MDCD_LINEAR` 只适合本类光滑场，不能
  从本例成功推导其适合激波。
- **边界**：去掉 x/z 周期连接、把上下壁改成 outflow，或只设 no-slip 不设热边界，都会定义
  不同问题。
- **停止**：降低 `max_steps` 只会更早失败；放宽残差阈值必须作为新的验收标准审查，不能把
  `maximum_steps` 文件改名为“收敛结果”。
- **输出**：`allow_existing=true` 可能混入旧文件；正式复现优先使用 `--clean` 或新目录。

任何修改后的派生算例都应更换 `case.name` 和 `output.directory`，保留独立 manifest、日志、
配置副本和容差说明，避免与本次冻结基线混淆。

## 10. 本算例引入的程序修改审计

case02 不只是输入文件补充。提交 `51b2820` 明确修改了以下程序能力：

- `src/runtime/flow_initializer.cpp` 与配置解析增加 `poiseuille` 解析初场，包括抛物线速度和
  可选黏性耗散—导热平衡温度多项式；
- 源项注册、配置与计算增加 `pressure_gradient`，把配置量定义为
  $\boldsymbol G=-\nabla p$，同步加入动量源 $\boldsymbol G$ 和能量功
  $\boldsymbol u\cdot\boldsymbol G$；
- `src/mesh/structured_mesh.cpp` 对来自 CGNS 单精度周期变换参数的坐标校验采用与其存储精度
  相称的容差，非周期连接仍保留原双精度严格容差；
- 发布网格生成器增加三维 `periodic-channel`，验证器增加体积加权 Poiseuille 剖面检查；
- 对上述配置、初场、源项和生成路径增加单元测试，并同步运行时/算法文档。

同一提交还加入本目录的两份配置、驱动和初始说明。提交 `5367d9a` 没有继续改变求解器离散，
但扩充 `run_case02.py` 的完整结果复用检查，调整正式 CFL/残差阈值，并写入实际运行产物和
第 3600 步人工强停分析。这里列出的修改是 case02 结果可复现性的一部分；它们不能在比较旧
可执行文件时忽略，也不能描述成纯粹的算例后处理变化。
