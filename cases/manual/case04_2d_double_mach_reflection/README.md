# case04：二维经典双马赫反射

本目录给出 Woodward--Colella 经典双马赫反射问题的可执行 WCNS 配置、统一 CGNS 网格生成
方法和运行脚本。它既检验强激波捕捉，也检验专用的分段/时变物理边界、SSPRK 子步时间传递、
多块连接和 MPI。仓库现已附带 480×120、16-rank、`t=0.2` 的 WENO5 与 MDCD_HYBRID
实际结果、验证文件和图像；完整结论及“SCMM6 度量已通过、纯 D6 尚未通过”的适用范围见
[`COMPARISON_REPORT.md`](COMPARISON_REPORT.md)。960×240 文件保留为高分辨率复算模板。

实际比较的统一入口是：

```powershell
python run_case04_comparison.py --method both --ranks 16 `
  --mpi-exec "C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin\mpiexec.exe"
```

两组使用 `algorithm.profile=scmm6_wcns`、特征重构和 Roe；为克服纯 D6 在初始激波足的
负内能问题，正式结果显式使用 `algorithm.flux_difference=conservative_two_point`。
这不是两套 profile 的交叉混用，但通量散度已经降阶，不得把结果标为纯 SCMM6-D6。

## 1. 物理问题

计算域为

\[
(x,y)\in[0,4]\times[0,1],\qquad \gamma=1.4.
\]

一条与 x 轴成 60 度的 Mach 10 激波在初始时刻与下壁交于
\(x_0=1/6\)。以激波线

\[
x_s(y,t)=x_0+\frac{y}{\sqrt 3}+\frac{20t}{\sqrt 3}
\]

划分状态；满足 \(x\le x_s\) 的一侧为激波后态，另一侧为激波前态。这里不把这个初始分片
称作问题的全时刻解析解：激波与下壁作用后会产生双马赫反射等复杂结构；上式只用于初场和
外边界随时间给定。

激波前态为

\[
(\rho,u,v,w,p)=(1.4,0,0,0,1),
\]

激波后态为

\[
(\rho,u,v,w,p)=
\left(8,\;8.25\frac{\sqrt3}{2},\;-4.125,\;0,\;116.5\right).
\]

全部数值是无量纲量。当前专用模型有意冻结为经典参数，不把 Mach 数、角度和两侧状态暴露
为可任意组合的配置项；这样可以避免产生不满足 Rankine--Hugoniot 关系的伪算例。唯一输入
参数 `initial.x0` 是激波足的 x 坐标，默认和本配置均为 `1/6`。

## 2. 边界条件及其实现含义

CGNS 网格外边界必须命名为 `left`、`right`、`bottom`、`top`。配置采用：

| 边界 | 配置类型 | 处理 |
|---|---|---|
| `left` | `double_mach_reflection` | 三层 ghost 和强约束面状态均取固定激波后态 |
| `top` | `double_mach_reflection` | 按真实边界面中心坐标及当前 SSPRK 子步时间，用移动激波线选择前/后态 |
| `bottom`, \(x<x_0\) | `double_mach_reflection` | 固定激波后态 |
| `bottom`, \(x\ge x_0\) | `double_mach_reflection` | 静止滑移反射壁，法向速度镜像、切向速度保留 |
| `right` | `outflow`（默认） | 无目标态时从内部外推 |

实现只平均真实边界面的顶点坐标来获得面中心，不生成或读取 ghost 坐标、ghost 度量或其他
ghost 派生量。先从边界条件生成三层 ghost 上的 \(\rho,u,v,w,T\)，再由气体模型一致地换算
\(p\) 和守恒量；无粘重构和粘性导数若存在都会看到同一份物理 ghost 内容。本算例无粘，
强面约束仍默认开启，以在重构后再次严格施加上述真实边界面状态。

专用边界只允许放在 i-lower、j-lower、j-upper；把它赋给右边界或 k 面会在推进前明确失败。
配置还会拒绝非二维、`gamma != 1.4`、粘性开启、源项开启、专用初场与边界不成套等组合。

## 3. 网格

目标网格是 960×240 个二维单元，物理长度 4×1。发布生成器把 x 方向分成 4 个原生
CGNS zone，每块 240×240 单元；相邻块用 `GridConnectivity1to1` 连接。y 方向不分块，保证
`bottom/top` 边界命名和分段边界操作直观。`periodic_x=false`，所以左右边界均为物理边界。

先在仓库根目录完成 MPI Release 构建，例如：

```powershell
cmake -S . -B build-rc-mpi -G "MinGW Makefiles" -DWCNS_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-rc-mpi --parallel 4
```

再进入本目录并生成网格：

```powershell
Set-Location cases\manual\case04_2d_double_mach_reflection
New-Item -ItemType Directory -Force grids,logs | Out-Null
..\..\..\build-rc-mpi\wcns_generate_release_cgns.exe rectangle grids\double_mach_960x240.cgns 960 240 4 4.0 1.0 false
```

参数从左到右依次是：子命令、输出文件、总 i 单元数、总 j 单元数、x 向 zone 数、x 长度、
y 长度、x 是否周期。若改变分辨率，i 单元数必须能被 zone 数整除；若改变物理域、边界名或
拓扑，就不再是本配置已经检查过的经典布置，必须同步复核边界和 `x0`。

等价的自动命令是：

```powershell
python run_case04.py --generate-only
```

脚本默认在 `build-rc-mpi` 查找程序；构建目录不同可使用 `--generator`、`--run` 指定绝对路径。

## 4. 配置文件逐项解释

配置文件是 `double_mach_960x240.wcns`。`mesh.path` 和 `restart.path` 相对配置文件所在目录
解析，`output.directory` 相对启动进程的工作目录解析；本脚本固定以本目录为工作目录，因此
三者不会混淆。手工运行时也建议先进入本目录。

### 4.1 schema、算例名和网格

```text
schema_version = 1
case.name = case04-double-mach-reflection-960x240
mesh.path = grids/double_mach_960x240.cgns
```

`schema_version` 必须是程序支持的 1；`case.name` 进入输出文件名、日志和 manifest；
`mesh.path` 指向刚生成的 CGNS。未知键、重复键和缺少必需键都会终止启动。

### 4.2 空间算法

```text
algorithm.profile = phenglei_wcns
algorithm.reconstruction = weno_z
algorithm.reconstruction_variables = characteristic
algorithm.riemann = hllc
```

- `phenglei_wcns` 选择成套 PHengLEI 风格几何度量路径；可整体改成 `scmm6_wcns` 对比，但
  不能在同一次运行中交叉混用两套 profile 的面度量、体 Jacobian 或模板。
- `weno_z` 用六点调用模板重构左右面状态；强激波算例不建议把正式结果改成 `linear5` 或
  `zero_order`。可替换为 `weno_js`、`mdcd_linear`、`mdcd_hybrid` 做算法比较。
- `characteristic` 在共同面特征基内重构；可替换 `primitive` 或 `conservative`，但应当把
  结果差异作为算法差异记录，而不是覆盖原结果。
- `hllc` 分辨接触波优于 Rusanov；可换 `rusanov` 或 `roe`。Roe 异常时有已记录的回退链。

若选择 MDCD，可在这组配置后增加：

```text
algorithm.mdcd.disp = 0.0463783
algorithm.mdcd.diss = 0.01
```

两项分别控制 MDCD 色散与耗散系数，必须满足有限、`disp>0`、`0<=diss<disp` 和
`3*disp+9*diss<1`；即使当前重构不是 MDCD，出现非法显式值仍会被拒绝。

### 4.3 气体和参考量

```text
gas.gamma = 1.4
gas.specific_gas_constant = 1.0
reference.velocity = 1.0
reference.density = 1.0
reference.temperature = 1.0
reference.length = 1.0
reference.viscosity = 1.0e-5
```

专用模型要求 `gamma=1.4`。这里所有参考尺度取便于无量纲解释的值；黏度参考值仍是 schema
必需项，但 `run.viscous=false`，不会产生粘性通量。不要通过开启黏性或源项把这个经典无粘
验证改成另一个问题。

### 4.4 分区

```text
partition.mode = auto_split
partition.allow_idle_ranks = false
partition.max_load_ratio = 1.2
partition.min_cells_per_active_direction = 8
```

网格本身已有四块。4 rank 时通常每 rank 一块；rank 数更多时，`auto_split` 可继续沿结构
方向二分，只要子块仍满足 WCNS 模板最小宽度。`allow_idle_ranks=false` 要求每个 rank 都有
合法叶块，否则在分配流场前失败。改变 rank 数不应改变原 zone 物理定义或共享面唯一通量。

### 4.5 初场和物理边界

```text
initial.type = double_mach_reflection
initial.x0 = 0.16666666666666667

boundary.default = outflow
boundary.left.type = double_mach_reflection
boundary.bottom.type = double_mach_reflection
boundary.top.type = double_mach_reflection
```

初场按每个真实单元中心和 `t=0` 激波线分片。三个显式覆盖名必须与 CGNS 边界名完全一致；
`right` 没有覆盖，继承 `outflow`。不能把普通 `inflow` 的 `rho/u/v/temperature` 数据混入专用
边界。`x0` 同时进入初场、三类专用边界数据和重启签名；修改后旧 checkpoint 会被拒绝。

### 4.6 方程和时间推进

```text
source.enabled = false
run.mode = unsteady
run.viscous = false
run.cfl = 0.2
run.max_steps = 1000000
run.t_end = 0.2
```

物理停止目标是无量纲时间 0.2。求解器会裁剪最后一步使时间事件准确落在目标时刻；
`max_steps` 只是异常保护，不是本算例的预期停止理由。若最终 manifest 是
`maximum_steps`，算例没有完成。CFL 0.2 是保守起点，不保证任意算法替换都稳定；出现正性
失败时先检查配置、边界、网格和日志中的降阶/回退，再决定是否减小 CFL。

### 4.7 流场输出

```text
output.directory = results/hllc_weno_z
output.allow_existing = false
output.dimensional = false
output.field.enabled = true
output.field.format = cgns
output.field.every_time = 0.02
output.field.write_initial = true
output.field.write_final = true
output.field.quantities = rho,u,v,p,T,mach,entropy_proxy
```

每 0.02 时间单位输出一次原 zone 重组后的 CGNS，并另写初场和终场。`allow_existing=false`
防止覆盖旧证据；重跑前应归档结果或显式安全清理本算例的 `results`。`dimensional=false` 表示
字段保持求解器无量纲值。可以增加已注册字段，但未知字段会立即失败。

### 4.8 历史、统计和检查点

```text
output.history.enabled = true
output.history.format = txt
output.history.every_steps = 20
output.history.write_initial = true
output.history.write_final = true

output.statistics.enabled = true
output.statistics.format = txt
output.statistics.every_steps = 20
output.statistics.write_initial = true
output.statistics.write_final = true
output.statistics.quantities = total_mass,total_momentum_x,total_momentum_y,total_energy

output.checkpoint.enabled = true
output.checkpoint.every_time = 0.05
output.checkpoint.write_initial = false
output.checkpoint.write_final = true
```

history 记录步数、时间、dt、五分量残差、回退计数和停止状态；statistics 记录全场积分，便于
识别异常质量/能量漂移。由于边界持续有通量交换，这些总量不应被误判为周期域守恒常数。
checkpoint 每 0.05 时间单位保存一次，并在终止时保存最终状态；重启仍需匹配网格和数值签名。

## 5. 建议运行流程

第一步只验证网格、CGNS 边界名、分区、度量、初场、配置和输出装配：

```powershell
python run_case04.py --dry-run-only
```

`--dry-run` 不推进时间，不能证明激波反射结果正确，但应当以 0 退出，并在
`logs/dry-run.log` 给出解析后的算法、分区和网格摘要。

第二步进行 4 rank 正式运行：

```powershell
python run_case04.py --ranks 4
```

若 MPI 启动器不在 PATH：

```powershell
python run_case04.py --ranks 4 --mpi-exec "C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin\mpiexec.exe"
```

脚本的 `--clean` 只删除本目录下精确命名的 `grids`、`logs`、`results`，属于显式破坏性操作；
使用前先归档需要保留的结果。没有 `--clean` 时会复用已有网格，且求解器会因
`output.allow_existing=false` 拒绝覆盖已有结果。

## 6. 完成后应得到什么

文件名含 case name、step、time 和 rank 数，具体格式以程序 manifest 为准。至少应有：

- `grids/double_mach_960x240.cgns`：四个原生结构 zone、真实顶点坐标、物理边界和块连接；
- `logs/generate-grid.log`、`logs/dry-run.log`、`logs/run-r4.log`：三阶段完整标准输出；
- `results/hllc_weno_z/*.field.*.cgns`：初场、每 0.02 时刻和 `t=0.2` 终场；
- `*.history.r4.txt`：时间步、残差、CFL、回退和停止原因；
- `*.statistics.r4.txt`：配置选择的全场积分随时间记录；
- `*.checkpoint.*.cgns`：0.05、0.10、0.15、0.20 附近由精确事件调度提交的检查点；
- `*.manifest.r4.txt`：版本、配置/网格/分区/重启签名、最终 step/time、停止原因和成功文件清单。

人工接受前至少检查：所有输出有限且 \(\rho,p,T>0\)；最终时间为 0.2 且停止理由是
`physical_time_reached`；初始激波足位于 1/6；下壁右段无穿透；顶面移动激波位置使用实际 RK 时间；块连接
两侧无接缝；三重点、滑移线和射流结构位置合理；网格加密时主要结构收敛；改变合法 rank 数
后积分、残差和终场在规定容差内一致。可视化图像是必要的物理检查，但不能替代上述机器
可读证据。

## 7. 可替换项与不可静默替换项

可以建立新的输出目录分别比较重构、Riemann、几何 profile、CFL、分辨率和 rank 数。每次
比较必须保存配置、Git 提交、日志、manifest 和结果，不能在同一目录覆盖。下列变化会改变
问题定义，必须在新报告中明确说明：域尺寸、`x0`、激波角/Mach/状态、边界名或边界类型、
开启粘性/源项、把顶部改为普通外推、把下壁全部设为滑移壁。当前程序不接受后几类“看似
相近”的设置作为经典双马赫反射专用模型的一部分。
