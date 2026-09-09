# case05：\(Re_\tau=180\) 三维周期湍流槽道（Linux 迁移前可行性卡口）

本目录配置一个由定常体积力驱动的三维可压缩周期槽道：

\[
(x,y,z)\in[0,2\pi]\times[-1,1]\times[0,\pi],\qquad
(N_x,N_y,N_z)=(36,48,36),\qquad Re_\tau=180.
\]

当前阶段只回答“该算例能否被读入、分区、并行推进五步并生成完整输出”。
它**不是已达统计定常的湍流计算，也不构成 DNS 精度验收**。在 Linux
服务器长算之前，必须先通过本次人工审阅。

## 1. 目录内容

| 路径 | 用途 |
|---|---|
| `channel_retau180_feasibility.wcns` | 冻结的 4-rank、5 步可行性配置 |
| `run_case05.py` | 生成网格、MPI dry-run、5 步推进和自动校验 |
| `validate_case05.py` | 检查统计列、五步终止、壁面统计正性与初始 \(Re_\tau\) |
| `grids/*.cgns` | 2×2 原生 zone 的 x/z 双周期、y 向壁面加密 CGNS 网格 |
| `logs/` | 网格生成、dry-run 和真实短算的标准输出 |
| `results/feasibility-r4/` | 初/终场、终止检查点、history、statistics 和 manifest |
| `validation/` | 终场有限性及槽道统计检查结果 |

## 2. 物理量、无量纲化与 \(Re_\tau\)

按本次订正，取槽道半高 \(h\)、摩擦速度 \(u_\tau\) 和参考密度作基本尺度：

\[
L_{ref}=h=1,\qquad U_{ref}=u_\tau=1,\qquad
\rho_{ref}=1,\qquad T_{ref}=71.42857142857143.
\]

摩擦速度、摩擦雷诺数和速度壁面单位定义为

\[
u_\tau=\sqrt{\frac{|\bar\tau_w|}{\rho_w}},\qquad
Re_\tau=\frac{\rho_w u_\tau h}{\mu_w},\qquad
U_b^+=\frac{U_b}{u_\tau}.
\]

因为速度尺度已经是 \(u_\tau\)，程序从参考量导出的 Reynolds 数就是目标摩擦
Reynolds 数：

\[
Re_{ref}=\frac{\rho_{ref}u_\tau h}{\mu_{ref}}=Re_\tau=180,
\qquad \mu_{ref}=\frac1{180}=0.005555555555555556.
\]

初始速度型积分仍给出 \(U_b^+=U_b/u_\tau=15.481978793165828\)，所以配置中的
无量纲体积平均速度是 \(U_b^+\)，不再是 1。以半高和体积平均速度定义的另一个 Reynolds
数可以从它导出：

\[
Re_b^{(h)}=\frac{\rho_bU_bh}{\mu_b}
\simeq U_b^+Re_\tau=2786.756182769849.
\]

请注意不同文献可能以全高 \(2h\) 定义体系雷诺数；那个数是本文的两倍，
即 5573.512365539698。配置和程序全程使用半高定义，不应在中途换约定。

平均壁面距离和速度的黏性无量纲量为

\[
y_w=h-|y|,\qquad y^+=\frac{\rho_wu_\tau y_w}{\mu_w}
\simeq Re_\tau\frac{y_w}{h},\qquad u^+=\frac{\bar u}{u_\tau}.
\]

这里最后一个等号只在初始常密度、常黏度设定下成立。后续长算分析应使用
程序实测的 \(\rho_w,\mu_w,u_\tau\)，不应永远把 180 代入横坐标。

本次订正要求程序日志中的参考马赫数为 0.1。配置取

\[
R=1,\qquad \gamma=1.4,\qquad
T_{ref}=71.42857142857143,\qquad U_{ref}=u_\tau=1,
\]

使

\[
c_{ref}=\sqrt{\gamma RT_{ref}}=10,\qquad
Ma_{ref}=\frac{U_{ref}}{c_{ref}}=0.1.
\]

启动日志因此应打印 `Re=180 Ma=0.1`。压力以 \(\rho_{ref}u_\tau^2\) 无量纲化，
无量纲初始 \(\rho=T=1\) 对应

\[
p=\frac{\rho T}{\gamma Ma_{ref}^2}=71.42857142857143.
\]

**必须区分参考马赫数与体积平均马赫数。**因 \(U_b/U_{ref}=U_b^+\)，当前设定导致

\[
Ma_b=\frac{U_b}{c_{ref}}=U_b^+Ma_{ref}=1.548197879316583,
\]

所以它不是通常意义上“体积平均 Ma=0.1”的低速经典槽道。实测五步终场最大局部
Mach 约为 1.914。如果物理目标其实是 \(Ma_b=0.1\)，在保持 \(U_{ref}=u_\tau,R=1\) 时应改为
\(Ma_{ref}=0.1/U_b^+=0.0064591226571\)、\(T_{ref}\simeq17120.83338\)，或重新改回以 \(U_b\)
为速度尺度。本次按批示保留 \(U_{ref}=u_\tau,Ma_{ref}=0.1\) 的组合，不自行替换物理目标。

## 3. 定常体积力

无量纲质量力加速度设为

\[
a_x^*=\frac{a_xh}{u_\tau^2}=1.
\]

因为 \(h=1\)，它与理想充分发展槽道的平均压力梯度平衡
\(|\mathrm d\bar p/\mathrm dx|=\bar\tau_w/h\) 一致。`body_force` 在方程中加入

\[
S=(0,\rho a_x,0,0,\rho u a_x)^T,
\]

即同时处理 x 动量和外力做功，不是只对速度每步人工加一个常数。这是
定体积力模式；它不会反馈调节体积流量。若长算必须恒流量，需新增闭环迫力
模型，不能把当前常量迫力误称为恒流量算法。

## 4. CGNS 网格、周期连接与分辨率

生成命令的完整参数是

```text
periodic-channel grids/channel_36x48x36_beta1p75.cgns \
  36 48 36 2 2 6.283185307179586 2.0 3.141592653589793 1.75 -1.0
```

从左到右依次为：总单元数 \(N_x,N_y,N_z\)，x/z 原生 zone 数，长度
\(L_x,L_y,L_z\)，双壁加密强度 \(\beta\)，y 起点。末参数 `-1.0` 不能遗漏；
如果省略，生成器为保持旧用法会默认产生 \([0,2]\) 而不是 \([-1,1]\) 的通道。

x/z 均匀；y 顶点映射是

\[
\eta_j=\frac{j}{N_y},\qquad
y_j=-1+2\,\frac12\left[1+
\frac{\tanh\{\beta(2\eta_j-1)\}}{\tanh\beta}\right],
\quad j=0,\ldots,N_y,\quad\beta=1.75.
\]

网格包含 2×2=4 个原生 CGNS zone，每块 18×48×18 单元。x 和 z 首尾及块间均
用 `GridConnectivity1to1` 表达，不是物理 `periodic` BC 字符串；y 下/上壁物理 patch
名为 `bottom`/`top`。4 rank dry-run 确认每 rank 分配一个原生 zone。

按目标 180 估算的网格间距为

\[
\Delta x^+=Re_\tau\frac{2\pi}{36}=31.416,qquad
\Delta z^+=Re_\tau\frac{\pi}{36}=15.708.
\]

第一个单元中心位于 \(y_w^+\simeq0.821\)，第一个顶点间隔
\(\Delta y_w^+\simeq1.701\)，中心附近最大 \(\Delta y^+\simeq13.918\)。这些数值表明网格
适合并行与功能可行性测试，但 36×48×36 明显是稀疏网格；不应用它宣称获得定量
DNS 统计。

本例不再采样 x-z 壁平行层。两个目标 y-z 截面是 \(x=0\) 和 \(x=\pi\)；
因为流场存在单元中心，程序统一选择目标正 x 侧最近的共面单元中心层。对本均匀 x 网格，
实际位置为

\[
x_{sample,0}=\frac{\pi}{36}=0.08726646259971647,\qquad
x_{sample,1}=\pi+\frac{\pi}{36}=3.2288591161895095.
\]

每个真实 y-z 截面面积为 \(A_{yz}=2\pi\)，而不是 x-z 壁平行面。

## 5. 初始流场：参考思路后的重新设计

用户提供的片段采用壁律型平均速度加三角扰动，这一总体思路被保留；
但初场并非复制该代码。实际实现做了下列适配：

1. 原式是 Reichardt 类复合壁律，不将其误称为 Spalding 隐式壁律。
2. 两面壁用 \(y_w=h-|y|\) 对称构造，而不在代码中复制两个分支。
3. 增加很小的中心线修正，使镜像平均速度在中心线一阶导数为零。
4. 扰动只用 1--3 阶整数 x/z 模态，确保在 36×36 网格上充分分辨且严格周期。
   参考片段直接使用了 20--50 的频率系数；按其分母和本例 x/z 周期长度化简后，
   其中有较高波数且部分是半整数周期模态，不满足本 CGNS 首尾周期连接的函数相容性，
   因而没有照搬。
5. 三个速度扰动均乘壁面包络，在两面壁上为零；x-z 平面平均为零，
   不会人为改变指定平均速度。

令

\[
\hat y=y/h,\quad s=1-|\hat y|,\quad y^+=Re_\tau s,
\]

平均壁律为

\[
f(y^+)=\frac{\ln(1+0.41y^+)}{0.41}
+7.8\left[1-e^{-y^+/11}-\frac{y^+}{11}e^{-y^+/3}\right].
\]

令 \(C=Re_\tau f'(Re_\tau)\)，则中心修正后

\[
u_m^+(s)=f(Re_\tau s)-\frac{C}{8}s^8,qquad
\left.\frac{\mathrm d u_m^+}{\mathrm ds}\right|_{s=1}=0,qquad
\bar u_m=U_b\frac{u_m^+}{U_b^+}.
\]

`bulk_velocity_plus=15.481978793165828` 是对这个修正型的连续半槽积分值。离散网格初始
体平均速度为 15.4871222867，与连续目标 15.4819787932 的差别是稀疏网格离散积分误差。

记 \(\theta_x=2\pi x/L_x\)、\(\theta_z=2\pi z/L_z\)、
\(A=0.05U_b(1-\hat y^2)^2\)，则扰动为

\[
\begin{aligned}
u'&=A[\sin(2\theta_x)\cos\theta_z+\tfrac12\sin(3\theta_x+2\theta_z)],\\
v'&=A[\sin\theta_x\sin\theta_z+\tfrac12\cos(2\theta_x-\theta_z)],\\
w'&=A[\cos\theta_x\sin(2\theta_z)-\tfrac12\sin(3\theta_x-\theta_z)].
\end{aligned}
\]

初始 \(\rho=1,T=1\)，压力和守恒量均由程序的同一气体模型转换得到。扰动用于
触发三维非定常性，不保证离散无散、不是湍流库采样，也不代替足够长的过渡过程。

## 6. 边界、粘性和空间算法

- x/z 方向周期性由 CGNS 1-to-1 连接及 MPI halo 通信实现。
- `bottom`/`top` 均是静止、无滑移、等温 \(T_w=1\) 壁。物理 ghost 先获得
  \(\rho,u,v,w,T\)，再统一闭合压力与守恒量；无粘通量和粘性导数使用同一份 ghost
  物理量。
- `run.viscous=true` 启用常粘性 Navier--Stokes 通量。当前没有湍流/RANS/LES 模型；
  这是一个稀疏的直接非定常离散可行性算例。
- 几何路径用 `phenglei_wcns`，重构用 `mdcd_hybrid` 原始变量，无粘面通量用 HLLC。
  若替换为 `scmm6_wcns`，必须整体切换 profile，不能交叉组合两套度量。
- `algorithm.mdcd.disp=0.0463783`、`diss=0.01` 是本例冻结的色散/耗散参数。

## 7. 配置文件逐组说明与可替换项

### 7.1 算例、气体和参考量

`schema_version=1` 是当前输入 schema；`case.name` 进入每个输出文件名；`mesh.path`
相对配置文件解析。参考量不能为达到某个 Re/Ma 而直接添加 `Re=` 或 `Ma=` 键：
程序只从 `reference.*`、`gas.*` 导出并在启动日志中打印 Re/Ma。修改任一参考量都
必须重新做第 2--3 节的整体推导，不可只改 `initial.re_tau`。

### 7.2 初场

`initial.type=turbulent_channel` 必须用三维网格；`y0/y1`、`x0/z0`、`period_x/period_z`
必须与网格一致。`re_tau`、`bulk_velocity`、`bulk_velocity_plus` 联合定义平均型缩放。
`perturbation_amplitude` 目前为 0.05，允许范围是 `[0,0.5]`；设为 0 可做层流化诊断，但不是
湍流起始条件。改动初场公式后必须重新计算 `bulk_velocity_plus`。

### 7.3 时间推进和停止设置

当前是专用短测：

```text
run.mode = unsteady
run.viscous = true
run.cfl = 0.15
run.max_steps = 5
run.t_end = 0.01
run.max_wall_time = 0
```

`t_end=0.01` 被故意设得比五步预期时间大，因此本地证据应以 `maximum_steps`
停止，`wcns_run` 对该安全保护停止返回码 2，运行脚本仅在这一个预期命令上接受 2。

服务器长算前至少要复制配置到新文件、使用新输出目录，并同时改大：

```text
run.max_steps = 5000000
run.t_end = 500
run.max_wall_time = <略小于作业墙钟上限的秒数>
```

这三个数字只是配置示例，不是已证明足以统计收敛的终值。非定常计算的物理目标是
`t_end`；`max_steps` 是异常保护；`max_wall_time>0` 会在墙钟到限前写 checkpoint 并安全返回 2。
长算不能改成 `steady`：槽道湍流瞬时场本身不会收敛到定常解。

### 7.4 流场、history、statistics 和 checkpoint 间隔

当前五步测试为了最强可见性，history/statistics 每步写，流场只写初场与终场，
checkpoint 只会因 `write_final=true` 在终止时写。长算应按存储容量与物理时标调整：

| 配置 | 当前 | 长算起点建议 | 含义 |
|---|---:|---:|---|
| `output.field.every_steps` | 0 | 0 或较大步数 | 0 关闭步数场输出 |
| `output.field.every_time` | 0 | 10 或更大 | 按物理时间写全场，会裁剪时间步对齐事件 |
| `output.history.every_steps` | 1 | 100--1000 | 残差、CFL、回退和停止状态 |
| `output.statistics.every_steps` | 1 | 0 或稀疏保险值 | 按步数采样统计 |
| `output.statistics.every_time` | 0 | 0.1--1 | 按物理时间采样统计 |
| `output.checkpoint.every_time` | 5 | 5--20 | 可重启 CGNS 状态 |

表中只是起点。不要同时设置过密的 `every_steps` 和 `every_time`；两类事件取并集。
每次新运行修改 `output.directory`，并保持 `output.allow_existing=false` 来防止覆盖证据。

## 8. 统计量的准确含义

### 8.1 \(x=0\) 和 \(x=\pi\) 附近的真实 y-z 截面

配置为

```text
output.statistics.yz_planes.enabled = true
output.statistics.yz_planes.target_x_coordinates = 0.0,3.141592653589793
```

对每个目标 \(x_t\)，程序在所有 MPI 叶块中寻找满足 \(x_c\ge x_t\) 的最近常 x
单元中心面。若存在与目标重合的单元中心面，则直接选中它。选择“正 x 侧”而非绝对距离最小，
可避免在周期端点 \(x=0\) 同时累加首尾两个等距平面。当前实现要求截面为几何平面；
非常 x 曲面会明确失败。

两个目标按配置顺序自动追加：

| 目标 | 实际单元中心 x | 平均速度列 | 质量流量列 |
|---|---:|---|---|
| \(x=0\) | \(\pi/36=0.0872664626\) | `yz_mean_u_plane0` | `yz_mass_flow_x_plane0` |
| \(x=\pi\) | \(\pi+\pi/36=3.2288591162\) | `yz_mean_u_plane1` | `yz_mass_flow_x_plane1` |

数学定义是

\[
\langle u\rangle_{yz,n}
=\frac{\displaystyle\int_{A_{yz,n}}u\,\mathrm dA}
       {\displaystyle\int_{A_{yz,n}}\mathrm dA},\qquad
\dot m_{x,n}=\int_{A_{yz,n}}\rho u\,\mathrm dA.
\]

面积在离散上使用该 cell-i 层左右两个 I 面面积的平均。`yz_mean_u_plane*`
是真实 y-z 截面上的面积加权流向平均速度；`yz_mass_flow_x_plane*` 是穿过该面的
x 向质量流量，不再是旧 x-z 壁平行积分指标。无量纲输出时，后者的量纲尺度是
\(\rho_{ref}U_{ref}L_{ref}^2\)。

当前短测初始两截面质量流量分别为 97.31403751990630 和
97.31403751990624，周期初场一致到浮点精度。五步后分别为 97.31386368944678 和
97.31386329908302；小差异来自非定常瞬时场在不同 x 截面的局部演化。

### 8.2 两面壁摩擦

`output.statistics.channel_walls.enabled=true` 会自动追加：

- `channel_wall_shear_lower/upper`：下/上壁面积加权的 \(|\tau_{wx}|\)；
- `channel_wall_shear_mean`：两壁再按面积合并的平均值；
- `channel_friction_velocity`：\(\sqrt{\bar\tau_w/\bar\rho_w}\)；
- `channel_re_tau`：按实测壁面密度、壁温黏性和 `half_height` 得到的 \(Re_\tau\)。

速度法向导数使用与当前 algorithm profile 相匹配的壁面 Dirichlet 单边模板；壁压由内部
高阶迹量插值，壁密度由等温壁温与状态方程得到。该统计只接受三维、平面 x-z、
J-lower/J-upper 的等温无滑移壁，不符合时明确失败，不会静默给出错误数据。

以上都是**瞬时空间统计**。当前程序不保存从某个起始时间开始的累积时间平均、
RMS 或 Reynolds 应力。长算应保留 statistics 时序列后离线做时间平均，且只在初始过渡被排除
后开始统计。

## 9. Windows 本地可行性操作

在仓库根目录中配置并构建 MPI Release：

```powershell
cmake -S . -B build-rc-mpi -G "MinGW Makefiles" -DWCNS_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-rc-mpi --parallel 4
ctest --test-dir build-rc-mpi -C Release --output-on-failure
```

然后一键运行冻结卡口：

```powershell
python cases\manual\case05_3d_turbulent_channel\run_case05.py --clean --ranks 4
```

`--clean` 只删除该 case 目录下的 `grids/logs/results/validation`；它会清除上次证据，
所以只在已归档或确认可重生时使用。脚本顺序是：

1. 生成 36×48×36 CGNS 网格；
2. 在 4 rank 上执行 `wcns_run --dry-run`，检查输入、连接、分区、度量和初场；
3. 真实推进 5 步，预期由 `maximum_steps` 安全停止；
4. 用发布校验器检查终场每个样本有限且 \(\rho,p,T>0\)；
5. 检查槽道统计文件的列集、步数、时间、壁摩擦正性与初始 \(Re_\tau\) 范围。

若还要复核报告中的终场最大局部 Mach，在根目录执行：

```powershell
build-rc-mpi\wcns_validate_release_case.exe nonzero cases\manual\case05_3d_turbulent_channel\results\feasibility-r4\case05-channel-retau180-feasibility.field.step00000005.time4p727996516eM04.cgns Mach 0
```

这里 `nonzero ... Mach 0` 会遍历 CGNS 中的 `Mach` 场并打印 `max_abs`；阈值 0 仅要求该场
不是全零，本命令的主要用途是独立读取并报告最大值。

可以只生成网格：

```powershell
python cases\manual\case05_3d_turbulent_channel\run_case05.py --generate-only
```

或只做 MPI dry-run：

```powershell
python cases\manual\case05_3d_turbulent_channel\run_case05.py --dry-run-only --ranks 4
```

构建目录不是 `build-rc-mpi` 时，用 `--run`、`--generator`、`--validator` 传入实际绝对路径。
Windows 找不到 `mpiexec` 时，用 `--mpi-exec` 指定启动器。已归档配置的输出目录固定含 `r4`，
因而录制验收只允许 `--ranks 4`；其他 rank 数需复制配置并使用新目录。

## 10. Linux 迁移预检流程（本阶段不执行长算）

审阅通过后，在 Linux 服务器上首先只做构建、单元测试和同样的 5 步卡口：

```bash
git clone https://github.com/wylx-2/wcns_v3.git
cd wcns_v3/wcns
git checkout stage/release
cmake -S . -B build-linux-mpi -DCMAKE_BUILD_TYPE=Release -DWCNS_ENABLE_MPI=ON
cmake --build build-linux-mpi --parallel 8
ctest --test-dir build-linux-mpi --output-on-failure
python3 cases/manual/case05_3d_turbulent_channel/run_case05.py \
  --clean --ranks 4 \
  --run build-linux-mpi/wcns_run \
  --generator build-linux-mpi/wcns_generate_release_cgns \
  --validator build-linux-mpi/wcns_validate_release_case \
  --mpi-exec mpiexec
```

注意 `run_case05.py` 会先把传入路径解析成绝对路径，上述命令应在 `wcns` 根目录
执行。若集群使用 Slurm，请将 MPI 启动器换成管理员支持的 `srun`/`mpiexec`，并先用
一个短作业确认 MPI 与编译器 ABI、CGNS 读写、节点间共享文件系统、栈大小和返回码行为。
不要把 Windows 下的 `.exe`、CMake cache 或构建目录复制到 Linux；只迁移 Git 追踪的源码与算例资产。

短卡口通过后再从本配置复制出专用长算配置，改变输出目录、时间上限、步数上限、
墙钟上限和输出间隔，但先不改物理、网格、算法和扰动，以便进行跨平台对照。

## 11. 当前实际可行性结果

已在 Windows、Release、4 MPI rank 上运行完整脚本。实际证据为：

| 项目 | 结果 |
|---|---:|
| 终止原因 | `maximum_steps`（预期的安全短测停止） |
| 终止步数/时间 | 5 / 0.0004727996515916612 |
| 每步 dt | 约 0.00009456 |
| 重构/Riemann 回退 | 0 / 0（所有五步） |
| 终场最小 \(\rho\) | 0.9969360990527756 |
| 终场最小 \(p\) | 71.12238987088477 |
| 终场最小 \(T\) | 0.9987736015763192 |
| 终场最大局部 Mach | 1.913982395650608 |
| 初始/终止实测 \(Re_\tau\) | 177.9758103 / 178.7424517 |
| 初始/终止平均 \(|\tau_w|\) | 0.9776354646 / 0.9856895927 |
| 初始 y-z 截面流量（0/1） | 97.31403751990630 / 97.31403751990624 |
| 终止 y-z 截面流量（0/1） | 97.31386368944678 / 97.31386329908302 |
| 初始/终止总质量 | 39.4784176043553 / 39.4784176043563 |

初始离散实测 \(Re_\tau\) 比设计值低约 1.12%，来自稀疏壁法向网格上的单边导数离散。
下/上壁初值分别为 0.9776354646456509 和 0.9776354646456257，对称性达到浮点精度。

这些数据支持如下结论：CGNS 周期多块网格、三维槽道初场、两面等温壁、常体积力、
粘性 WCNS、4-rank 通信、两个 y-z 截面流量/平均速度统计、壁摩擦统计、CGNS 场/检查点输出在短时推进中可用。
这些数据**不支持**湍流已发展、统计已收敛、平均速度符合参考 DNS，或稀疏网格达到 DNS 分辨率。

## 12. 长算的物理验收原则

服务器作业不应仅以“跑到 `t_end`”为成功标准。建议分成过渡段和采样段：

1. 过渡段持续监测体积平均流量指标、两壁摩擦、\(u_\tau\)、\(Re_\tau\) 与总量；
2. 只在这些量的长时滑动平均无系统漂移、且两壁对称后，记录采样起始时间；
3. 采样段至少要覆盖多个大涡周转时间，并用分块平均或延长采样检查统计不确定度；
4. 用实测 \(u_\tau\) 构造 \(y^+\)、\(u^+\)，将时间平均速度型、壁摩擦和对称性与
   \(Re_\tau\approx180\) 的公开基准数据比较；
5. 查看 manifest/history 确认没有数值失败，并报告所有降阶/黎曼回退计数。

在当前 36×48×36 网格上，即使时间平均看似平稳，也只应视为稀疏网格数值试验。
要进行经典 \(Re_\tau=180\) DNS 定量对照，必须先做网格收敛/加密研究并补充时间累积湍流统计量。

## 13. 本 case 为程序增加的功能

为使配置真正可执行，本阶段对程序做了四类显式修改，不存在隐藏修改：

1. 新增 `turbulent_channel` 初场类型、严格参数校验和初场单元测试；
2. `periodic-channel` 网格生成器增加可选 `origin_y`，用于直接生成 \([-1,1]\) 法向区间；
3. 新增两面等温 J 壁的平均摩擦应力、摩擦速度和实测 \(Re_\tau\) 统计，包含配置检查和注册测试；
4. 新增按目标 x 坐标选择正 x 侧最近平面的 y-z 截面平均流向速度和质量流量统计，
   并补充配置与注册单元测试。

体积力、周期连接、粘性通量、MPI 分区、输出与重启均是已有通用功能，
本 case 只对它们进行组合和验证。

## 14. 人工审阅卡口

请在批准 Linux 迁移/长算前至少检查：

- 物理域、网格数、4 个 zone、x/z 周期和 y 双等温壁是否与要求一致；
- \(Re_\tau\)--\(Re_b^{(h)}\)--\(U_b^+\)--常体积力推导是否接受；
- 重新设计的低模态、壁面衰减扰动是否符合预期；
- 是否接受当前为定体积力而非恒流量闭环驱动；
- 是否接受 y-z 采样面按目标正 x 侧最近单元中心选取；
- 是否接受“5 步工程可行性通过，湍流/DNS 物理验收未开始”的结论；
- 长算前是否需要先新增累积时间平均/RMS/Reynolds 应力或恒流量控制器。

人工审阅通过前，不启动 Linux 长时间计算。
