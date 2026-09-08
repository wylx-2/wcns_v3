# 人工验收 case01：二维四象限 Riemann 问题

本目录是 case01 的完整、可复现验收包。它在同一物理区域
$[0,1]\times[0,1]$ 上比较两套 $256\times256$ 结构多块网格：

- `uniform_256x256.cgns`：两个方向均匀；
- `clustered_256x256_x0p5_y0p5.cgns`：分别向 $x=0.5$、$y=0.5$ 光滑加密。

两套计算采用相同的四象限 Riemann 初值，间断交点是 $(0.8,0.8)$，用 4 个 MPI rank
推进到无量纲物理时间 $t=0.3$。本算例主要人工检查网格读取、多块连接、MPI 分区与 halo、
物理边界、间断初值、WCNS 无粘通量、非定常停止以及 CGNS/文本输出是否能在较大二维网格上
连成一条完整生产路径。

## 1. 目录内容

```text
case01_2d_riemann/
|-- README.md                         本文：设置、命令、结果和修改方法
|-- uniform_256x256.wcns              均匀网格计算配置
|-- clustered_256x256.wcns            中心加密网格计算配置
|-- run_case01.py                     生成、计算、校验和汇总驱动
|-- case01-summary.json               机器可读的参数、耗时和校验汇总
|-- files.sha256                      本次生成文件的 SHA-256 清单
|-- grids/
|   |-- uniform_256x256.cgns
|   `-- clustered_256x256_x0p5_y0p5.cgns
|-- results/
|   |-- uniform/                      均匀网格的场、历史、统计和 manifest
|   `-- clustered/                    加密网格的场、历史、统计和 manifest
|-- logs/                             网格生成和两次求解的完整标准输出
`-- validation/                       初/终场有限性和终场对角对称性检查
```

`grids/`、`results/`、`logs/`、`validation/`、`case01-summary.json` 和
`files.sha256` 都可由驱动重新生成。配置、驱动和本文是算例定义，不应由清理操作删除。

## 2. 网格设置

### 2.1 拓扑和 CGNS 内容

两套网格都有 256×256 个二维单元、257×257 个顶点。生成器沿 $i$ 方向建立四个原生
CGNS zone，每个 zone 是 64×256 个单元。相邻 zone 的 `1to1` 连接覆盖整条公共面；外边界
写入 `left`、`right`、`bottom`、`top` 四类物理边界。这样 4-rank 运行时每个 rank 恰好持有
一个 64×256 叶块，同时仍然实际经过跨块状态和面通量通信，而不是把问题当作单块串行计算。

网格坐标和连接由 `wcns_generate_release_cgns` 写成 CGNS ADF 文件。求解器读取后会重新检查
点范围、连接变换、覆盖关系、正 Jacobian 以及全局网格签名。

### 2.2 均匀网格

逻辑顶点 $s_i=i/256$，两个方向都采用

$$x_i=s_i,\qquad y_j=s_j.$$

因此

$$\Delta x=\Delta y=\frac1{256}=0.00390625.$$

$0.5$ 是第 128 个顶点；$(0.8,0.8)$ 不要求落在顶点上，初值根据**单元中心坐标**判断。

### 2.3 向 `(0.5,0.5)` 加密的网格

令加密中心 $c=0.5$、强度 $a=2.5$，对 $x$ 和 $y$ 各自独立使用同一个分段映射

$$
g(s)=
\begin{cases}
c\left[1-\dfrac{\sinh\!\left(a(c-s)/c\right)}{\sinh(a)}\right],&s\le c,\\[6pt]
c+(1-c)\dfrac{\sinh\!\left(a(s-c)/(1-c)\right)}{\sinh(a)},&s>c.
\end{cases}
$$

坐标为 $x_i=g(i/256)$、$y_j=g(j/256)$。该映射在 $s=c$ 两侧连续且一阶导数连续，端点仍
严格为 0 和 1，并且第 128 个顶点严格落在 0.5。当前参数给出

| 指标 | 数值 |
|---|---:|
| 最小网格间距（靠近 0.5） | 0.0016142009619946607 |
| 最大网格间距（靠近外边界） | 0.009803377113974776 |
| 最大/最小间距比 | 6.073204753… |

这不是向初始间断 $(0.8,0.8)$ 加密；它有意检查高阶算法在光滑非均匀坐标和跨 zone
连接上的一致性。如果要把加密中心移到间断处，可把生成命令的两个 `0.5` 改成 `0.8`，但应
同步改网格文件名，并重新运行 `--dry-run` 和全部有限性检查。改变强度参数 `2.5` 可以调整
疏密比；值越大中心附近越密、边界附近越疏，过大会造成极小 CFL 时间步和较差网格质量。

## 3. 四象限 Riemann 初值

程序以单元中心 $(x_c,y_c)$ 判定区域：`x >= x0` 为东侧，`y >= y0` 为北侧。当前
`x0 = y0 = 0.8`，无量纲原始变量为

| 区域 | 判定 | $\rho$ | $u$ | $v$ | $p$ |
|---|---|---:|---:|---:|---:|
| NE | $x_c\ge0.8,\ y_c\ge0.8$ | 1.5 | 0 | 0 | 1.5 |
| NW | $x_c<0.8,\ y_c\ge0.8$ | 0.5323 | 1.206 | 0 | 0.3 |
| SW | $x_c<0.8,\ y_c<0.8$ | 0.138 | 1.206 | 1.206 | 0.029 |
| SE | $x_c\ge0.8,\ y_c<0.8$ | 0.5323 | 0 | 1.206 | 0.3 |

这是常用的 Lax--Liu 二维 Riemann 配置 3。$w=0$；程序由理想气体关系计算温度，再计算
五个守恒量。初值和网格关于 $x=y$ 对称，所以终场也应在浮点与并行误差范围内保持对角
对称，这是本例除有限性以外的一项强检查。

如果修改任一象限的状态，必须保持 $\rho>0$、$p>0$。若将 `initial.x0` 与 `initial.y0`
改成不同数值，或给 $x/y$ 使用不同网格映射，当前 `diagonal-symmetry` 验收就不再适用，需从
`run_case01.py` 删除或替换该检查，不能简单放宽容差来掩盖不对称设置。

## 4. 配置文件逐项说明

两份 `.wcns` 除 `case.name`、`mesh.path` 和 `output.directory` 外完全相同。配置采用严格的
`key = value` 语法：未知键、重复键、非法枚举和非有限数字会在分配流场前失败。

### 4.1 标识、网格和无粘算法

```text
schema_version = 1
case.name = case01-uniform-256x256
mesh.path = grids/uniform_256x256.cgns
```

- `schema_version` 固定当前配置协议；不能随意改成别的整数。
- `case.name` 是所有输出文件前缀；加密算例对应 `case01-clustered-256x256`。
- `mesh.path` 相对于配置文件所在目录解析；第二份配置指向加密 CGNS。

```text
algorithm.profile = phenglei_wcns
algorithm.reconstruction = weno_z
algorithm.reconstruction_variables = characteristic
algorithm.riemann = hllc
```

- `profile` 选择 PHengLEI-WCNS 成套的度量、中心到面插值、通量差分和边界闭合。
  可整体替换为 `scmm6_wcns`，但两套 profile 的部件不可交叉拼装；切换后应重新验收。
- `reconstruction` 是左右面状态重构。可替换为 `weno_js`、`mdcd_linear` 或
  `mdcd_hybrid`。`mdcd_linear` 对强间断不具备非线性抑振能力，本例不推荐单独使用。
- `reconstruction_variables` 可选 `conservative`、`primitive`、`characteristic`。当前特征
  重构开销更高，但更适合含多波系间断的验收。
- `riemann` 可替换为 `rusanov`、`roe`。Rusanov 更耗散但稳健；Roe 分辨率较高并带熵修正；
  当前 HLLC 在接触间断分辨率和稳健性间折中。
- 时间推进器当前由生产求解路径固定为 SSPRK3，不是本 schema 的可选键。
- 低 Mach 预处理尚未实现，所以配置中没有预处理开关，本非定常激波算例也不应使用稳态
  低 Mach 预处理。

### 4.2 气体模型和无量纲参考量

```text
gas.gamma = 1.4
gas.specific_gas_constant = 1.0
reference.velocity = 1.0
reference.density = 1.0
reference.temperature = 1.0
reference.length = 1.0
reference.viscosity = 1.0
```

- `gamma=1.4` 表示量热完全理想气体。`gas.specific_gas_constant` 与 `gas.molar_mass` 必须
  二选一；若改成摩尔质量，应删除前者再设置后者，不能同时保留。
- 五个 `reference.*` 分别是 $U_{ref},\rho_{ref},T_{ref},L_{ref},\mu_{ref}$，均必须为正。
  本例全部取 1，所以配置表中的数值就是内部无量纲值。
- 压力标度是 $\rho_{ref}U_{ref}^2$。程序只由参考输入计算并在日志报告
  $Re=\rho_{ref}U_{ref}L_{ref}/\mu_{ref}$ 和
  $Ma=U_{ref}/\sqrt{\gamma R T_{ref}}$；配置中不能直接输入 Re 或 Ma。
- 本例 `run.viscous=false`，所以参考黏度和导出的 Re 不进入无粘右端，但参考黏度仍是完整
  无量纲配置的必填项。

### 4.3 MPI 分区

```text
partition.mode = auto_split
partition.allow_idle_ranks = false
partition.max_load_ratio = 1.2
partition.min_cells_per_active_direction = 8
```

- `auto_split` 在原 CGNS zone 少于 rank 数时继续确定性剖分；本例 4 zone/4 rank，无需再切。
  可选 `zones_only`（禁止切分）或 `force_split`（即使 zone 已足够也尝试切分）。
- `allow_idle_ranks=false` 要求每个 rank 都获得叶块；rank 过多或块过窄会明确失败。
- `max_load_ratio=1.2` 把最大叶块负载限制为平均负载的 1.2 倍目标。
- `min_cells_per_active_direction=8` 防止切出小于高阶模板需要的窄块。当前 profile 的硬下限
  是 4；保留 8 提供安全余量。

用不同 MPI 数运行是允许的，但输出文件名、manifest 和性能结果会改变。若 rank 数不是 4，
分区器可能合并多个 zone 到同一 rank，或进一步切分现有 zone；先执行 `--dry-run` 检查摘要。

### 4.4 初场、边界和源项

```text
initial.type = quadrant_riemann
initial.x0 = 0.8
initial.y0 = 0.8
initial.ne_rho = ...
...
initial.se_p = ...
```

- `quadrant_riemann` 启用第 3 节的四区初始化器。
- `x0/y0` 是间断判定坐标。
- 每个 `ne/nw/sw/se` 前缀依次显式给出 `rho,u,v,p`。显式列出全部值可以避免依赖内建
  默认值，使算例定义可审计。

```text
boundary.default = outflow
source.enabled = false
```

- 所有未单独覆盖的 CGNS 物理边界使用外流 ghost 状态。该设置适合本例在有限时间内让波系
  接近或穿过边界；因为存在边界通量，全域质量、动量和能量并不要求恒定。
- 边界可按 CGNS patch 名改成 `boundary.<patch>.type = ...`，支持 `farfield`、`inflow`、
  `outflow`、`slip_wall`、两类无滑移壁、`symmetry` 和 `periodic`。改变边界类型会改变物理
  问题，不能继续与本次基线数值直接比较。
- `source.enabled=false` 关闭方程源项。若开启，必须同时按运行指南设置
  `source.models` 和相应的守恒源、体力或制造解参数；二维 z 动量源必须为零。

物理边界 ghost 只由边界条件填充可用物理量及守恒量；算法不会依赖物理边界 ghost 的坐标、
度量或二级派生量。

### 4.5 时间推进和停止条件

```text
run.mode = unsteady
run.viscous = false
run.cfl = 0.4
run.max_steps = 1000000
run.t_end = 0.3
run.max_wall_time = 0
```

- `unsteady` 表示真实物理时间推进，残差只监测、不作为收敛停止条件。
- `viscous=false` 求解二维无粘 Euler 方程。
- `cfl=0.4` 用全局稳定步长推进；可减小以增强稳健性，但计算步数增加。增大 CFL 必须重新
  做正性、对称性和网格收敛检查，不能只看程序是否退出。
- `max_steps` 是防失控硬上限。正常情况下不会接近一百万步。
- `t_end=0.3` 是本例的首要停止条件；最后一步会裁剪，保证最终时间严格命中 0.3。
- `max_wall_time=0` 关闭墙钟安全停止。若设为正值，运行配置还必须启用检查点，才能在完整
  时间步边界安全退出。

本次正常停止原因应为 `physical_time_reached`。如果 manifest 显示 `maximum_steps`、
`numerical_failure` 或墙钟/信号停止，就不能判定 case01 通过。

### 4.6 流场输出

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
output.field.quantities = rho,u,v,p,T,rho_u,rho_v,rho_w,rho_E,sound_speed,mach,total_enthalpy,entropy_proxy,jacobian
```

- 两个配置写入不同子目录；该路径相对于启动进程的工作目录解析，所以驱动先切换到 case
  目录再运行。
- `allow_existing=false` 防止误覆盖已有结果；驱动的 `--clean` 会先安全删除本 case 的生成
  子目录。手工重复运行前也必须选择新输出目录或明确清理旧结果。
- `dimensional=false` 输出无量纲数据。改成 `true` 时，程序按参考量恢复量纲。
- 格式可替换为 `tecplot` 或 `both`。CGNS 保留原 zone 和 `CellCenter` 场；Tecplot 是可直接
  可视化的 ASCII ordered zone，但不是重启文件。
- `every_steps=0/every_time=0` 关闭中间场快照，只写初场和终场，避免 256² 算例产生大量
  文件。可设 `every_steps=100`，或 `every_time=0.05`；非定常求解会裁剪步长以精确命中时间
  输出事件。
- `quantities` 依次输出密度、速度、压力、温度、四个动量/能量守恒量、声速、Mach 数、总焓、
  熵代理和 Jacobian。可以删减列表；加入未注册名称会在首次输出时失败而不是写伪数据。

### 4.7 残差、全域统计和检查点

```text
output.history.enabled = true
output.history.format = txt
output.history.every_steps = 100
output.history.write_initial = true
output.history.write_final = true
```

残差历史每 100 步以及初/终状态写一行。固定列包含 step、time、dt、CFL、墙钟时间，总残差、
五守恒分量的 L2/Linf、参考/归一化残差、连续满足次数、重构和 Riemann 回退计数、检查标志与
停止原因。格式可改为 `tecplot`，但该固定 schema 不接受 quantities 列表。

```text
output.statistics.enabled = true
output.statistics.format = txt
output.statistics.every_steps = 100
output.statistics.write_initial = true
output.statistics.write_final = true
output.statistics.quantities = total_mass,total_momentum_x,total_momentum_y,total_energy
```

全域统计与残差同频输出，使用正 Jacobian 和原 zone 守恒积分权重进行 MPI 归约。它们用于
发现异常漂移并观察开放边界的净通量，不应在本例的 `outflow` 边界下强行要求守恒为初值。
可加入已注册的 `total_momentum_z`，或删去不关心的量。

```text
output.checkpoint.enabled = false
```

本次连续短算例不写重启检查点。长时间或有墙钟上限的计算应启用检查点并设置步/时间调度；
检查点保存无量纲守恒场和重启签名，与普通流场 CGNS 的用途不同。

## 5. 从构建到验收的运行步骤

以下命令从仓库 `wcns/` 根目录运行，PowerShell 路径以本项目已验证的 Windows 环境为例。

### 步骤 1：配置并编译 Release MPI 程序

已有 `build-rc-mpi` 时只需增量构建：

```powershell
cmake --build build-rc-mpi --parallel 4
```

该命令编译求解器、CGNS 网格生成器和发布算例校验器。若从空构建目录开始，先执行：

```powershell
cmake -S . -B build-rc-mpi -G "MinGW Makefiles" `
  -DCMAKE_BUILD_TYPE=Release -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_MPI=ON
cmake --build build-rc-mpi --parallel 4
```

第一条生成启用 CGNS 与 MPI 的 Release 构建系统，第二条真正编译。编译器、MPI 和 CGNS
依赖必须已按项目运行指南配置；也可以使用另一个构建目录，并通过后述脚本参数传入可执行文件。

### 步骤 2：一键生成、计算和校验

```powershell
python -B cases/manual/case01_2d_riemann/run_case01.py --clean
```

`-B` 禁止生成 Python 字节码缓存；`--clean` 只删除本 case 目录内的六类生成物，然后按以下
顺序执行：

1. 生成均匀和中心加密 CGNS 网格；
2. 用 Intel MPI 各以 4 rank 运行均匀、加密算例；两次求解串行执行，避免互相争抢资源；
3. 对两个初场和两个终场逐变量检查有限值、正密度和正压力；
4. 对两个终场执行 $x/y$ 交换后的对角对称性检查，容差为 `5e-4`；
5. 写入 `case01-summary.json` 和所有生成文件的 `files.sha256`。

驱动默认使用 `build-rc-mpi` 和 Intel MPI。替换构建目录、MPI 启动器或 rank 数的示例：

```powershell
python -B cases/manual/case01_2d_riemann/run_case01.py --clean `
  --run build-other/wcns_run.exe `
  --generator build-other/wcns_generate_release_cgns.exe `
  --validator build-other/wcns_validate_release_case.exe `
  --mpi-exec "C:/path/to/mpiexec.exe" --ranks 8
```

改变 rank 数后生成的性能数据不再与本次 4-rank 基线等价，但数值有限性与合理误差范围内的
解应保持一致。

### 步骤 3（可选）：只检查启动装配

从 case 目录执行：

```powershell
cd cases/manual/case01_2d_riemann
../../../build-rc-mpi/wcns_run.exe --config uniform_256x256.wcns --dry-run
../../../build-rc-mpi/wcns_run.exe --config clustered_256x256.wcns --dry-run
cd ../../..
```

`--dry-run` 会真正读取配置和 CGNS、完成 MPI/分区装配、度量、初场与边界启动检查，但不做
时间推进，也不写正式结果。它适合在更换网格、rank 数或算法后先快速发现配置和几何错误。
若检查多 rank 分区，应在命令前加 `mpiexec -n <N>`。

### 步骤 4（可选）：拆开执行以定位问题

在 case 目录下单独生成网格：

```powershell
../../../build-rc-mpi/wcns_generate_release_cgns.exe rectangle `
  grids/uniform_256x256.cgns 256 256 4 1.0 1.0 false
../../../build-rc-mpi/wcns_generate_release_cgns.exe clustered-rectangle `
  grids/clustered_256x256_x0p5_y0p5.cgns `
  256 256 4 1.0 1.0 0.5 0.5 2.5 false
```

第一条参数依次是输出、$i/j$ 单元数、$i$ 向 zone 数、$x/y$ 长度和是否令 x 周期；第二条
额外给出 x/y 加密中心和映射强度。

单独求解：

```powershell
& "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe" -n 4 `
  ../../../build-rc-mpi/wcns_run.exe --config uniform_256x256.wcns
& "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin/mpiexec.exe" -n 4 `
  ../../../build-rc-mpi/wcns_run.exe --config clustered_256x256.wcns
```

`mpiexec -n 4` 启动四个进程；每个进程读取 rank 0 广播的统一配置，获得一个叶块，迭代中
交换多块 halo 和共享面通量。程序退出码 0 且 manifest 的停止原因正确，才表示推进完成。

单独校验终场时，把 `<final.cgns>` 替换为 manifest 记录的实际文件名：

```powershell
../../../build-rc-mpi/wcns_validate_release_case.exe finite <final.cgns>
../../../build-rc-mpi/wcns_validate_release_case.exe diagonal-symmetry <final.cgns> 5e-4
```

第一条遍历 CGNS CellCenter 解，检查全部请求字段有限，并报告 $\rho,p,T$ 的最小值，拒绝
NaN/Inf 和非正热力学状态；
第二条匹配 $(x,y)$ 与 $(y,x)$ 单元并比较应互换的速度/动量分量和应相同的标量。

### 步骤 5：核对文件完整性

在 case 目录用 PowerShell 可逐项复算：

```powershell
Get-Content files.sha256 | ForEach-Object {
  $expected, $relative = $_ -split '  ', 2
  $actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $relative).Hash.ToLower()
  if ($actual -ne $expected) { throw "checksum mismatch: $relative" }
}
```

这验证文件在验收生成后未被截断或修改；它不替代物理解检查。

## 6. 本次正式运行结果

本节的数值由 `run_case01.py` 在 4-rank Release 构建上实测，机器可读原始记录见
`case01-summary.json`，详细过程见 `logs/`，最终文件名以各自 manifest 为准。

<!-- CASE01_RESULTS_BEGIN -->

运行环境为 Windows 10、MinGW GCC 8.1.0 Release、Intel MPI 2021.10，求解器 manifest
记录的源码提交为 `6ca6e6ec062f`。两个网格都形成 4 个 64×256 叶块，rank 0--3 各持有一个，
分区摘要完全一致；网格签名分别为 `17874315645299610072` 和 `1192861213358331032`。

| 项目 | 均匀网格 | 中心加密网格 |
|---|---:|---:|
| MPI rank | 4 | 4 |
| 最终步数 | 841 | 1625 |
| 最终时间 | 0.3 | 0.3 |
| 停止原因 | `physical_time_reached` | `physical_time_reached` |
| 求解器内部墙钟时间 | 465.332 s | 911.065 s |
| 驱动测得进程总时间 | 466.190 s | 911.730 s |
| 进程树峰值 RSS | 676.234 MiB | 666.012 MiB |
| 重构回退累计数 | 0 | 0 |
| Riemann 回退累计数 | 0 | 0 |

加密网格用了约 1.93 倍步数，符合其更小物理网格间距对全局 CFL 时间步的限制。两次运行都
没有重构降阶或 HLLC 回退，说明本次稳定结果来自所配置的 WENO-Z 特征重构和 HLLC 主路径。

有限性检查遍历了每个终场的 $65536\times14=917504$ 个场样本。热力学下界为

| 状态 | $\min\rho$ | $\min p$ | $\min T$ |
|---|---:|---:|---:|
| 两套初场 | 0.1380000000 | 0.0290000000 | 0.2101449275 |
| 均匀终场 | 0.1293159781 | 0.01846925934 | 0.1428227170 |
| 加密终场 | 0.1225963117 | 0.01674918940 | 0.1328807455 |

所有请求字段均有限，且两套终场的密度、压力和温度保持严格为正。终场对角对称性为

| 网格 | L1 误差 | Linf 误差 | 验收容差 |
|---|---:|---:|---:|
| 均匀 | $5.8335\times10^{-15}$ | $2.0384\times10^{-13}$ | $5\times10^{-4}$ |
| 加密 | $2.1825\times10^{-7}$ | $1.0620\times10^{-4}$ | $5\times10^{-4}$ |

两者均通过。加密网格误差明显大于均匀网格的舍入量级，但其 L1 和 Linf 仍同时低于本例阈值；
后续若改变度量算法、分区或编译器，应把该数值作为需要持续观察的基线，而不是把容差继续
放宽。

全域积分的初值和终值如下：

| 网格/时间 | 总质量 | x 总动量 | y 总动量 | 总能量 |
|---|---:|---:|---:|---:|
| 均匀，$t=0$ | 0.3178610687 | 0.2091333569 | 0.2091333569 | 0.6868318710 |
| 均匀，$t=0.3$ | 0.4943018650 | 0.1923808798 | 0.1923808798 | 1.1268857052 |
| 加密，$t=0$ | 0.3214028297 | 0.2095434504 | 0.2095434185 | 0.6952799037 |
| 加密，$t=0.3$ | 0.4988988970 | 0.1922675933 | 0.1922675414 | 1.1383838110 |

两套离散初始积分略有不同，是因为 $(0.8,0.8)$ 不落在同一物理单元边界，按单元中心赋值后
两套网格对四个常状态区的离散覆盖面积不同。积分随时间改变则主要反映 `outflow` ghost
条件下边界上的实际净通量；本例不是封闭或周期守恒箱。

本次两个最终场文件为：

```text
results/uniform/case01-uniform-256x256.field.step00000841.time3p000000000eM01.cgns
results/clustered/case01-clustered-256x256.field.step00001625.time3p000000000eM01.cgns
```

连同初场、网格、文本历史、统计、manifest、日志和校验报告，目录共有 28 个文件，合计
36,383,963 bytes（约 34.70 MiB）。散列清单覆盖其中 23 个可重新生成的结果文件并已逐项
复算通过；算例定义文件和散列清单自身不在清单内。

<!-- CASE01_RESULTS_END -->

### 6.1 如何判断通过

本例同时满足以下条件才通过：

1. 两次求解进程退出码均为 0，manifest 均正常提交且停止原因均为
   `physical_time_reached`，最终时间为 0.3；
2. 初场、终场所有请求字段有限，$\rho,p,T$ 严格为正；
3. 终场对角对称误差不超过 `5e-4`；
4. 日志没有 MPI、CGNS、网格 Jacobian、halo/共享面、正性或 I/O 错误；
5. 输出清单完整，SHA-256 全部复算一致。

对称性通过说明数值离散没有明显破坏该特定对称关系，但不等价于证明解已网格收敛。均匀与
加密网格的单元中心不同，不能用逐索引相减作为误差；如需定量比较，须先在共同物理采样点上
保守投影或插值，并另设网格收敛算例。

## 7. 每类结果文件包含什么

- `*.field.step........time....cgns`：按原四个 zone 重组的 CGNS `CellCenter` 初场或终场。
  每个 FlowSolution 含配置列出的 14 个无量纲字段；网格坐标也随原 CGNS zone 可供后处理。
- `*.history.r4.txt`：残差和诊断时间序列。首行是固定列名，最后一行记录最终 step/time 与
  `physical_time_reached`。
- `*.statistics.r4.txt`：全域质量、x/y 动量和总能随时间的 MPI 归约结果。
- `*.manifest.r4.txt`：本次可执行版本/Git 提交、编译器、构建类型、rank 数、配置和分区摘要、
  网格签名、最终步/时间/步长/墙钟、停止原因及成功提交的输出文件列表。
- `logs/generate-*.log`：网格生成器确认的尺寸、zone 和映射参数。
- `logs/run-*.log`：启动摘要以及每一步的 time、dt、残差、检查标记与停止状态，是故障定位的
  首要证据。
- `validation/*-finite.txt`：初/终场的全局 min/max 和有限性判定。
- `validation/*-final-symmetry.txt`：对角对称性误差和容差判定。
- `case01-summary.json`：网格间距元数据、初始间断、算法、最终场路径，以及每个外部命令的
  命令行、返回码、墙钟时间、进程树峰值 RSS、日志和末行摘要。
- `files.sha256`：上述所有生成物（包括网格和 summary）的相对路径与 SHA-256；它本身不把
  自己列入清单，以避免递归摘要。

## 8. 后处理建议与允许的派生试验

可用支持 CGNS 的 ParaView、Tecplot 或自编 CGNS 读取器打开终场，至少绘制 `rho`、`p`、
`mach` 的二维等值图，并检查四区相互作用产生的激波、接触间断和稀疏波是否关于 $x=y$ 对称。
如果直接希望 Tecplot 文件，把 `output.field.format` 改为 `both` 后重算；不要把旧 CGNS 和
不同配置的新 Tecplot 混作同一次 manifest 结果。

在保留一个变量的基线时，可做以下派生试验：

- 固定网格，依次更换 WENO-JS/WENO-Z/MDCD-HYBRID，比较间断厚度、极值和回退计数；
- 固定重构，更换 Rusanov/HLLC/Roe，比较耗散和稳健性；
- 把输出调度改为 `every_time=0.05`，观察波系演化；
- 改用 1、2、8 rank，比较最终统计和共同物理点上的场，检查并行/分区不变性；
- 增加 128²、512² 网格形成真正的网格收敛序列。

每个派生试验都应使用新的 `case.name` 和 `output.directory`，保留独立 manifest，并在记录中
明确它已不再是本 README 的 case01 基线。
