# 人工验收 case03：二维等熵涡在周期扭曲网格上的一周期输运

本目录用二维等熵涡检验扭曲结构网格、四块 CGNS 周期连接、4-rank MPI 通信和两套完整
WCNS 几何路径。同一个 $100\times100$ 网格分别采用 `phenglei_wcns` 和 `scmm6_wcns`，从
$t=0$ 推进到 $t=10$。背景速度为 $(1,1)$，因此涡心在 x、y 两个方向都恰好跨过一个长度为
10 的周期域，终场解析解与初场相同。

本次两套计算均在第 1352 步由 `physical_time_reached` 正常停止，有限性、全局周期守恒、
解析密度 L1 以及度量卡口全部通过。SCMM6 在本分辨率、本扭曲映射和本算法组合下的初终场
密度 L1 为 $1.43474\times10^{-5}$，PHengLEI 路径为 $8.78036\times10^{-5}$；这说明本次
SCMM6 结果更准确，但单一网格不能证明一般精度阶或所有曲线网格上的优劣。

## 1. 目录和文件职责

```text
case03_2d_vortex/
|-- README.md                         本文：定义、命令、结果和修改审计
|-- Isentropic_curl.cpp               用户提供的参考 CGNS 网格生成代码
|-- phenglei_100x100.wcns             PHengLEI-WCNS 正式配置
|-- scmm6_100x100.wcns                SCMM6-WCNS 正式配置
|-- run_case03.py                     安全生成、运行、验证和汇总驱动
|-- case03-summary.json               机器可读的实测汇总
|-- files.sha256                      28 个可重建产物的 SHA-256
|-- grids/
|   `-- warped_100x100.cgns           四 zone、双周期二维扭曲网格
|-- results/
|   |-- phenglei/                     PHengLEI 初/终场、历史、统计、manifest
|   `-- scmm6/                        SCMM6 初/终场、历史、统计、manifest
|-- logs/
|   |-- generate-grid.log
|   |-- run-phenglei.log
|   `-- run-scmm6.log
`-- validation/
    |-- metric-profile-comparison.txt 两套完整度量场的直接比较
    |-- *-initial-finite.txt           初场有限性和正性
    |-- *-final-finite.txt             终场有限性和正性
    |-- *-initial-analytic-error.txt   初场相对解析式的基线误差
    |-- *-final-analytic-error.txt     t=10 相对解析式的误差
    |-- *-initial-final-error.txt      终场相对同路径初场的逐变量误差
    `-- *-conservation.txt             周期域全局守恒漂移
```

`run_case03.py --clean` 只删除 `grids/`、`results/`、`logs/`、`validation/`、summary 和散列
清单；不会删除本文、两个配置、驱动或参考代码。没有 `--clean` 时，如果已有 `results/`，驱动
会拒绝覆盖，避免把两次计算混入同一 manifest。

## 2. 如何从参考代码得到正式二维网格

### 2.1 参考代码中的坐标思路

`Isentropic_curl.cpp` 使用 101×101 个 x/y 顶点，并额外沿 z 挤出 8 层。case03 只采用其
x、y 映射思想，把逻辑坐标记为 $\xi,\eta\in[0,10]$：

$$
x(\xi,\eta)=\xi+0.5\sin\left(\frac{2\pi\eta}{10}\right),
$$

$$
y(\xi,\eta)=\eta+0.5\sin\left(\frac{4\pi\xi}{10}\right).
$$

第一式沿 y 方向有一个正弦周期，第二式沿 x 方向有两个正弦周期。这不是把 x、y 分别做
一维拉伸，而是含交叉导数的真正二维扭曲映射。

参考程序生成三维薄层单 zone，并用自连接表示周期；正式算例改成二维 2×2 原生 zone：

- 全局 100×100 单元、101×101 顶点；
- 每个 zone 为 50×50 单元；
- 内部块面写互逆 `1to1` 连接；
- 域首尾块面写带 $(10,0)$ 或 $(0,10)$ 平移的 CGNS 周期连接；
- 4 个 MPI rank 各持有一个原生 zone，`partition.mode=zones_only`，实际经过跨 rank halo 和
  共享面通量同步。

### 2.2 网格可逆性

连续映射 Jacobian 为

$$
\frac{\partial(x,y)}{\partial(\xi,\eta)}
=1-\frac{\pi^2}{50}
\cos\left(\frac{2\pi\eta}{10}\right)
\cos\left(\frac{4\pi\xi}{10}\right).
$$

因为 $\pi^2/50=0.1973920880$，连续 Jacobian 范围是

$$[0.80260791198,\ 1.19739208802],$$

严格为正。计算代码约定计算空间步长为 1，所以 100×100 离散网格的 J 具有单元面积意义，
还要乘 $\Delta\xi\Delta\eta=0.1^2$；因此离散 J 约在 0.008--0.012，而不是 0.8--1.2。

生成器接受一般振幅 $A_x,A_y$，并在写文件前要求

$$\frac{8\pi^2A_xA_y}{L^2}<1,$$

否则直接拒绝可能翻转的映射。生成正式网格的命令是：

```powershell
./build-rc-mpi/wcns_generate_release_cgns.exe `
  warped-periodic-square `
  cases/manual/case03_2d_vortex/grids/warped_100x100.cgns `
  100 100 10.0 0.5 0.5
```

参数依次是模式、输出文件、x/y 单元数、周期长度、x 扰动振幅和 y 扰动振幅。单元数必须为
偶数才能均分成 2×2 zone；还需保证每个活动方向至少满足 SCMM6 和分区配置要求，本例每块
50 个单元，远高于下限 8。

## 3. 等熵涡及完整周期

令周期最短有符号距离为

$$
\bar x=\operatorname{remainder}(x-x_c,10),\qquad
\bar y=\operatorname{remainder}(y-y_c,10),\qquad
r^2=\bar x^2+\bar y^2.
$$

背景状态为 $\rho_\infty=T_\infty=1$，背景速度 $(u_\infty,v_\infty)=(1,1)$，涡强
$\beta=5$。初始原始变量为

$$
u=1-\frac{\beta}{2\pi}e^{(1-r^2)/2}\bar y,
\qquad
v=1+\frac{\beta}{2\pi}e^{(1-r^2)/2}\bar x,
$$

$$
T=1-\frac{(\gamma-1)Ma^2\beta^2}{8\pi^2}e^{1-r^2},
\qquad
\rho=T^{1/(\gamma-1)}.
$$

由状态方程可得 $p=\rho T/(\gamma Ma^2)$。这里 $\gamma=1.4$，参考量不是直接输入 Ma：

$$
Ma=\frac{340}{\sqrt{1.4(8.314/0.029)(288.15)}}
=0.9997687921778579.
$$

涡心轨迹为

$$
(x_c(t),y_c(t))=(5+t,5+t)\pmod {10}.
$$

所以 $t=10$ 时涡心和完整解析场都回到初始周期位置。本例求解无粘可压缩 Euler 方程，不
启用源项；配置中由参考黏度导出的 Re 只写入日志，不参与无粘右端。

## 4. 为什么初值必须显式设置周期

扭曲网格在逻辑周期边附近可能出现 $x<0$、$x>10$、$y<0$ 或 $y>10$ 的物理坐标，但对应
网格点仍通过平移周期连接。原有 `isentropic_vortex` 初值只使用 `x-x0,y-y0`，在这种网格上
会把周期映像误认为距离更远的点，使接口两侧涡尾不完全一致。

case03 因此为主程序新增：

```text
initial.period_x = 10.0
initial.period_y = 10.0
```

正周期值让初值使用 `std::remainder` 得到最短周期距离；缺省或零值保留非周期旧行为；负值
在配置验证阶段拒绝。修正后 PHengLEI 初场相对解析式的密度误差为逐位 0，primitive L∞
仅 $7.77\times10^{-16}$。该修改、测试和提交号在第 10 节公开列出。

## 5. 两套配置逐段说明

两份配置只允许 `case.name`、`algorithm.profile` 和 `output.directory` 不同；其余网格、物理、
重构、Riemann、CFL、停止和输出设置完全一致。

### 5.1 网格和完整算法 profile

```text
schema_version = 1
case.name = case03-vortex-phenglei-100x100
mesh.path = grids/warped_100x100.cgns
algorithm.profile = phenglei_wcns
algorithm.reconstruction = weno_z
algorithm.reconstruction_variables = characteristic
algorithm.riemann = hllc
```

- schema 1 开启严格键检查，未知键、重复键会报错；
- `profile` 不是只切换 J：它整体绑定度量、公共线性插值、通量差分和边界闭合。另一配置用
  `scmm6_wcns`，两套部件不能交叉拼接；
- WENO-Z 在特征变量上重构，HLLC 计算无粘面通量。若更换 WENO-JS、MDCD 或 Roe/Rusanov，
  必须两份配置同步替换并另存为派生算例，不能把算法差异误记为纯度量差异。

### 5.2 气体和参考量

```text
gas.gamma = 1.4
gas.molar_mass = 0.029
reference.velocity = 340.0
reference.density = 1.225
reference.temperature = 288.15
reference.length = 1.0
reference.viscosity = 1.7894e-5
```

程序由摩尔质量计算比气体常数，再导出 Ma 和 Re；不能额外输入 Re/Ma。改变任一参考量会
改变无量纲涡温度和压力，因此必须重新生成两套结果，不能只在后处理中换算。

### 5.3 MPI 分区、初场和边界

```text
partition.mode = zones_only
partition.allow_idle_ranks = false
partition.max_load_ratio = 1.2
partition.min_cells_per_active_direction = 8
initial.type = isentropic_vortex
initial.x0 = 5.0
initial.y0 = 5.0
initial.beta = 5.0
initial.background_u = 1.0
initial.background_v = 1.0
initial.period_x = 10.0
initial.period_y = 10.0
boundary.default = farfield
source.enabled = false
```

4 zone/4 rank 下 `zones_only` 不进行二次切块。网格所有外表面均已成为周期 connectivity，
所以 `boundary.default=farfield` 实际不产生物理边界 ghost；若日志出现物理 farfield patch，
说明 CGNS 周期连接缺失，应判失败，不能靠边界默认值继续算。

### 5.4 时间推进和停止

```text
run.mode = unsteady
run.viscous = false
run.cfl = 0.5
run.max_steps = 100000
run.t_end = 10.0
run.max_wall_time = 0
```

SSPRK3 每步根据当前场和扭曲网格尺度更新全局 dt；末步自动截短，保证时间恰为 10。成功
停止原因必须是 `physical_time_reached`。`max_steps=100000` 只是失控保护，若先触发则本例
失败；`max_wall_time=0` 表示不设墙钟上限。CFL 可降低做时间误差检查，增大则必须重新验证
稳定性和误差。

### 5.5 输出

```text
output.directory = results/phenglei
output.allow_existing = false
output.dimensional = false
output.field.enabled = true
output.field.format = cgns
output.field.write_initial = true
output.field.write_final = true
output.field.quantities = rho,u,v,p,T,rho_u,rho_v,rho_w,rho_E,sound_speed,mach,total_enthalpy,entropy_proxy,jacobian
output.history.enabled = true
output.history.format = txt
output.history.every_steps = 100
output.history.write_initial = true
output.history.write_final = true
output.statistics.enabled = true
output.statistics.format = txt
output.statistics.every_steps = 100
output.statistics.write_initial = true
output.statistics.write_final = true
output.statistics.quantities = total_mass,total_momentum_x,total_momentum_y,total_energy
output.checkpoint.enabled = false
```

初场与终场都必须写出，才能直接按同一 profile 的离散自由度比较。场保持无量纲 CGNS；历史
每 100 步记录残差，统计每 100 步记录全局守恒量，并强制写末步。本例约两分钟一条路径，
未启用 checkpoint；若扩展到更细网格或长时间多周期，应另设 checkpoint 调度。

## 6. 构建和运行命令

从仓库 `wcns/` 根目录执行：

```powershell
cmake -S . -B build-rc-mpi `
  -DCMAKE_BUILD_TYPE=Release `
  -DWCNS_ENABLE_MPI=ON `
  -DWCNS_ENABLE_CGNS=ON `
  -DWCNS_BUILD_TESTS=ON

cmake --build build-rc-mpi --target `
  wcns_run wcns_generate_release_cgns `
  wcns_validate_release_case wcns_compare_metric_profiles `
  wcns_unit_tests -j 4

ctest --test-dir build-rc-mpi `
  -R "wcns\.(unit|cgns_reader)" --output-on-failure
```

第一条固定 Release、MPI、CGNS 和测试开关；第二条构建求解器、网格生成器、通用验证器、
度量比较器和单元测试；第三条运行本次相关回归测试。

只检查网格、配置、MPI 分配和两套度量构造而不推进：

```powershell
python -B cases/manual/case03_2d_vortex/run_case03.py --clean --dry-run
```

正式重建全部结果：

```powershell
python -B cases/manual/case03_2d_vortex/run_case03.py --clean
```

驱动依次执行：生成网格；独立构造/比较两套度量；PHengLEI 4-rank 一周期计算；初终场有限性、
解析式、初终场和守恒验证；再对 SCMM6 重复完全相同流程；最后写 JSON 和 SHA-256。默认密度
解析 L1 上限为 `1e-3`、守恒相对漂移上限为 `2e-10`。可以用相应命令行参数收紧阈值，但
正式基线不应通过放宽阈值覆盖失败。

## 7. 度量系数比较方法和实测结果

`wcns_compare_metric_profiles` 对每个**完整原生 CGNS zone**读取两份独立 block，一份只交给
`phenglei_wcns`，另一份只交给 `scmm6_wcns`。没有先算一套再转换为另一套，也没有共享线性
算子或用有限体积 J 替代高阶结果。

- PHengLEI 路径在顶点间构造二倍加密网格，用 `grid_delta` 对称形式计算 J 和面面积向量，
  再采样到真实单元/面；
- SCMM6 路径先以六阶顶点到中心插值生成公共 cell-center 坐标，用 SCMM D6/D4 算子计算
  对称守恒度量，再用 SCMM I6 得到面向量；
- 两者均独立与低阶有限体积单元面积比较，并检查 J 正性；本次 fallback 必须为 0；
- 每个单元还检查
  $\boldsymbol S_i(i+1)-\boldsymbol S_i(i)+\boldsymbol S_j(j+1)-\boldsymbol S_j(j)$ 的 L∞，
  这是离散几何守恒闭合诊断。

| 指标 | PHengLEI | SCMM6 |
|---|---:|---:|
| J 最小值 | 8.0312724700e-3 | 8.0309525382e-3 |
| J 最大值 | 1.1968727530e-2 | 1.1969047462e-2 |
| $\sum J$ | 100.00000000000061 | 99.99999999999969 |
| 相对有限体积 J 最大差 | 2.0026217007e-4 | 2.0077637255e-4 |
| fallback 单元 | 0 | 0 |
| GCL 闭合 L∞ | 2.2204460493e-16 | 4.2779027700e-14 |

两套路径彼此的差异为：

| 量 | L2 差 | L∞ 差 | 相对 L∞ |
|---|---:|---:|---:|
| J | 1.7934982543e-7 | 1.5180083071e-6 | 1.2682782919e-4 |
| i 面 x 分量 | 2.3130341599e-15 | 6.5392136150e-14 | 6.5392136150e-13 |
| i 面 y 分量 | 1.0133285876e-6 | 5.0665666097e-6 | 1.6135394761e-4 |
| j 面 x 分量 | 7.6483808360e-6 | 3.8232201702e-5 | 6.0968749262e-4 |
| j 面 y 分量 | 2.6361485517e-15 | 7.0901617910e-14 | 7.0901617910e-13 |

主要对角分量接近机器精度，交叉分量和 J 有预期的离散算子差异；两者总体积均正确、GCL
闭合远低于 $10^{-10}$，且跨 profile J 相对 L∞ 低于本例 $10^{-3}$ 卡口。

## 8. 一周期计算结果

### 8.1 运行和停止

| 路径 | 步数 | 最终时间 | 末步 dt | 求解器墙钟 | 停止原因 |
|---|---:|---:|---:|---:|---|
| PHengLEI | 1352 | 10 | 3.9036065586e-3 | 112.252314 s | `physical_time_reached` |
| SCMM6 | 1352 | 10 | 4.2223103629e-3 | 119.481111 s | `physical_time_reached` |

两套 history 的重构回退和 Riemann 回退始终为 0。末步 dt 不同是到达同一终止时间时的自动
截短量，不代表 SCMM6 使用更大的常规 CFL。

### 8.2 相对解析式的密度误差

| 路径/时刻 | 密度 L1 | 密度 L2 | 密度 L∞ | primitive L∞ |
|---|---:|---:|---:|---:|
| PHengLEI，初场 | 0 | 0 | 0 | 7.77156e-16 |
| PHengLEI，$t=10$ | 8.5989193e-5 | 1.1313617e-4 | 5.0680561e-4 | 2.0579789e-3 |
| SCMM6，初场 | 1.1076422e-5 | 3.9687158e-5 | 3.0471886e-4 | 7.1866360e-4 |
| SCMM6，$t=10$ | 2.0188364e-5 | 5.1383973e-5 | 4.5476495e-4 | 1.7911724e-3 |

解析验证器从输出 CGNS 顶点做算术平均确定显示用 CellCenter。PHengLEI 的内部中心与该约定
一致，因此周期修正后初场解析误差为零；SCMM6 内部中心由 I6 插值得到，与算术平均中心有
最多约 $9.87\times10^{-4}$ 的坐标差，所以 SCMM6 的“初场解析误差”含采样位置基线，不能
全部解释成流场初始化错误。为避免这一位置约定混入输运误差，本例的主比较使用下一节同一
profile 初场—终场误差。

### 8.3 终场相对同路径初场的体积加权误差

L1/L2 用初场 J 加权，L∞ 为逐单元最大值：

| 变量 | PHengLEI L1 | PHengLEI L2 | PHengLEI L∞ | SCMM6 L1 | SCMM6 L2 | SCMM6 L∞ |
|---|---:|---:|---:|---:|---:|---:|
| $\rho$ | 8.78036e-5 | 1.15312e-4 | 5.06806e-4 | 1.43474e-5 | 2.79675e-5 | 2.80381e-4 |
| $u$ | 7.68729e-5 | 1.53640e-4 | 2.05798e-3 | 2.93951e-5 | 1.06315e-4 | 1.53938e-3 |
| $v$ | 1.25747e-4 | 1.95679e-4 | 1.41733e-3 | 2.56336e-5 | 8.60379e-5 | 1.08949e-3 |
| $p$ | 8.54182e-5 | 1.11011e-4 | 4.60840e-4 | 1.34560e-5 | 3.05636e-5 | 3.79487e-4 |
| $T$ | 4.21699e-5 | 8.96503e-5 | 1.60469e-3 | 1.50731e-5 | 7.54603e-5 | 1.53227e-3 |

SCMM6 的密度 L1 是 PHengLEI 的约 16.34%，即降低约 83.66%；压力 L2 降低约 72.47%，
y 速度 L2 降低约 56.03%。温度 L2 的改善较小，约 15.83%。这些结果表明两条实现都稳定并
保持涡结构，SCMM6 在这一个周期、这一个网格上总体耗散/相位误差更低。要声明正式收敛阶，
仍需至少增加 50²、200² 等相似扭曲网格并控制时间误差。

### 8.4 守恒和正性

| 路径 | 统计行数 | 最大相对守恒漂移 | 容差 |
|---|---:|---:|---:|
| PHengLEI | 15 | 7.9473445533e-14 | 2e-10 |
| SCMM6 | 15 | 8.0200735726e-14 | 2e-10 |

PHengLEI 质量从 97.631688451856675 变为 97.631688451849087；SCMM6 从
97.632283260154253 变为 97.632283260146806。两套初始积分略有不同，是两套 cell-center 坐标
和 J 离散不同，不是运行中失守恒。所有 14 个输出字段均有限，初终场最小值为：

| 路径/时刻 | min $\rho$ | min $p$ | min $T$ |
|---|---:|---:|---:|
| PHengLEI 初场 | 0.3505230607 | 0.1646933431 | 0.6574861943 |
| PHengLEI 终场 | 0.3505193691 | 0.1643184797 | 0.6559144644 |
| SCMM6 初场 | 0.3505224716 | 0.1646929556 | 0.6574857523 |
| SCMM6 终场 | 0.3505986232 | 0.1643749350 | 0.6559534856 |

## 9. 通过条件和结论边界

本次通过条件为：

1. 两套 manifest 均为 4 rank、同一 mesh signature、终止步 1352、时间 10，停止原因为
   `physical_time_reached`；
2. 两套 J 严格为正，fallback 为 0，$|\sum J-100|\le10^{-10}$，相对有限体积 J 差不超过
   $10^{-3}$，GCL 闭合不超过 $10^{-10}$；两套 J 彼此相对 L∞ 不超过 $10^{-3}$；
3. 初终场字段有限且 $\rho,p,T>0$；
4. $t=10$ 解析密度 L1 不超过 $10^{-3}$；
5. 周期域总质量、两分量动量、总能量的最大相对漂移不超过 $2\times10^{-10}$；
6. 重构/Riemann 回退为 0，日志无 MPI、CGNS、度量、正性或输出错误；
7. 28 个生成产物的 SHA-256 全部复算一致。

全部条件通过。该结论验证“在此 100² 扭曲网格上，两套完整 profile 均能完成一个周期并保持
精度/守恒”，不等同于证明任一算法达到理论阶，也不覆盖黏性曲线网格、非匹配接口、三维
扭曲或更强网格畸变。

## 10. 程序修改审计：本例没有隐藏改动

### 10.1 case03 新增或修改的程序

基础设施提交 `edbce01` 做了以下公开改动：

- `tools/generate_release_cgns.cpp` 增加 `warped-periodic-square`，严格采用参考代码的 x/y 映射，
  输出 2×2 周期 CGNS，并增加映射正 Jacobian 参数检查；
- 新增 `wcns_compare_metric_profiles` 及 CMake 目标，直接输出两套 J、面向量和 GCL 比较；
- `wcns_validate_release_case` 增加 `field-error`，按参考 J 输出逐字段 L1/L2/L∞；
- 增加两份 case03 配置、驱动、参考代码入库，更新发布验证说明；
- 追补 case01、case02 报告中的历史程序修改审计。

周期初场修正提交 `71786ec` 修改了：

- `src/runtime/flow_initializer.cpp`：等熵涡可按 `period_x/period_y` 使用周期最短距离；
- `src/runtime/case_config.cpp`：注册并验证两个新输入键；
- `tests/test_flow_initializer.cpp`、`tests/test_case_config.cpp`：覆盖周期等价点和负周期拒绝；
- `docs/runtime-guide.md` 和两份 case03 配置：公开输入语义并实际启用。

没有修改 `src/mesh/high_order_metrics.cpp`、WENO、Riemann、通量差分、SSPRK3、MPI halo、共享
面同步或停止控制来迎合本算例结果。正式 manifest 中嵌入 `git_commit=71786ec0544d`。

### 10.2 追补 case01 的程序改动

case01 的报告现已增加第 9 节审计。提交 `6ca6e6e` 只给发布 CGNS 生成器增加
`clustered-rectangle`，再加入配置和驱动；没有修改 `src/`/`include/` 的求解器算法。提交
`9472350` 记录实际结果和报告，不改主程序。

### 10.3 追补 case02 的程序改动

case02 的报告现已增加第 10 节审计。提交 `51b2820` 增加 Poiseuille 初场、压力梯度动量/能量
一致源项、三维周期通道生成器、Poiseuille 验证器，并针对 CGNS 单精度周期变换调整连接坐标
校验容差；同时增加测试和文档。提交 `5367d9a` 主要记录实际结果、强停分析和安全复用逻辑，
没有继续修改空间离散。

## 11. 结果文件如何使用

- `*.field.step00000000...cgns` 是初场，`*.field.step00001352...cgns` 是 $t=10$ 终场；每个
  文件按四个原生 zone 重组，含密度、速度、压力、温度、守恒量、Mach、声速、总焓、熵代理
  和 J；
- `*.history.r4.txt` 含 step/time/dt/CFL、五分量残差、回退数和停止原因；
- `*.statistics.r4.txt` 含全域质量、x/y 动量和总能；
- `*.manifest.r4.txt` 固定程序版本/Git、编译器、rank、配置/分区 digest、网格签名、终止状态
  和输出清单；
- `case03-summary.json` 汇总网格、物理、度量差、两套 manifest、解析误差、初终场误差和全部
  外部命令资源记录；
- `files.sha256` 用相对路径保护全部可重建产物。配置、驱动、参考代码、本文和散列清单自身
  不放入清单，避免把算例定义误当运行产物或产生递归散列。

后处理可在 ParaView/Tecplot 打开初终场，绘制 `Density`、`Pressure`、`Mach` 和 `Jacobian`。
由于两套内部 cell-center 定义不同，不应把 PHengLEI 与 SCMM6 CGNS 按物理坐标直接插值后
宣称严格逐点相等；优先比较各自相对解析式和各自初终场，跨 profile 几何差使用本目录专门的
metric 报告。
