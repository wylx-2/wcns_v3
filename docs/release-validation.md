# 发布算例基础设施

阶段 O 的数值验收只调用正式 `wcns_run`，不建立简化求解器。共同入口如下：

- `wcns_generate_release_cgns` 按总单元数、维数、x 向原生 zone 数、扭曲幅值和 x 周期
  开关生成确定性结构 CGNS；x 周期至少使用两个 zone，避免生成当前拓扑契约不支持的
  self-connectivity；
- `cases/config/*.wcns.in` 保存严格 schema v1 模板；驱动器只替换 `@...@` 占位符；
- `wcns_validate_release_case` 独立重读最终 CGNS，可检查有限性/正性、均匀解析值和两个
  串并行结果的逐场最大差；
- `tools/run_release_matrix.py` 建立隔离工作目录，生成网格和配置，调用正式程序及验证器，
  最终写出 `matrix-summary.json` 和每条命令的独立日志；每条记录同时包含端到端墙钟和采样到
  的进程树聚合峰值 RSS。Windows MPI 工作进程可能由服务进程启动，采样器同时追踪本次命令
  启动后出现的同名求解进程；Linux 使用 `/proc`。不支持的平台把 RSS 记为 `null`，不得伪造。

## 1. 网格生成器

```text
wcns_generate_release_cgns output.cgns dimension cells_i cells_j cells_k zones_i warp periodic_x
wcns_generate_release_cgns periodic-square output.cgns cells_i cells_j length
wcns_generate_release_cgns rectangle output.cgns cells_i cells_j zones_i length_x length_y periodic_x
wcns_generate_release_cgns invalid-one-sided output.cgns cells_i cells_j
```

`dimension` 为 2 或 3；二维要求 `cells_k=1`；`cells_i` 必须能被 `zones_i` 整除；
`abs(warp)<0.2`。非周期网格的 x 外边界名为 `left/right`，其余物理边界为
`bottom/top[/front/back]`；原生多区连接及周期连接均成对写入。扭曲使用跨 zone 连续的
解析坐标映射，因而同一参数总是产生逐位相同的网格。
`periodic-square` 生成 2×2 原生多 zone 且 x/y 双向平移周期的方形网格；`rectangle` 用于
Sod 薄域和四象限问题，可指定二维物理长度及 x 向原生 zone 数。
`invalid-one-sided` 仅用于失败路径验收：它生成两个 zone，但故意省略右区指回左区的互逆
连接；不得作为物理解算网格使用。

## 2. 独立验证器

```text
wcns_validate_release_case finite field.cgns
wcns_validate_release_case uniform field.cgns rho u v w T tolerance
wcns_validate_release_case compare lhs.cgns rhs.cgns tolerance
wcns_validate_release_case vortex field.cgns time length x0 y0 beta u0 v0 gamma Mach L1-tolerance
wcns_validate_release_case sod field.cgns time x0 gamma rho-L1-tolerance position-cell-tolerance
wcns_validate_release_case diagonal-symmetry field.cgns L1-tolerance
wcns_validate_release_case viscous-profile field.cgns couette|conduction Reynolds L2-tolerance pressure-tolerance
wcns_validate_release_case uniform-source field.cgns time U0[5] source[5] tolerance
wcns_validate_release_case tecplot-consistency field.cgns field.dat tolerance
wcns_validate_release_case derived field.cgns gamma viscosity Jacobian tolerance
wcns_validate_release_case nonzero field.cgns field-name minimum-maximum-absolute-value
```

验证器只使用 CGNS API 重读输出，不链接求解器或其内存对象。`compare` 要求 zone 名、尺寸、
字段集合完全一致；所有比较同时拒绝非有限值。后续 O2--O4 在这个可执行程序中增加光滑误差、
守恒、截面和事件表等子命令，不改变已有命令语义。

## 3. 矩阵驱动

最小串行复现示例（路径按实际构建目录替换）：

```text
python tools/run_release_matrix.py \
  --run build/wcns_run \
  --generator build/wcns_generate_release_cgns \
  --validator build/wcns_validate_release_case \
  --template cases/config/freestream.wcns.in \
  --work-dir build/release-smoke --ranks 1
```

MPI 使用 `--mpi-exec <mpiexec> --ranks 1,2,4`。原生多区/周期/三维扭曲分别通过
`--zones-i`、`--periodic-x`、`--dimension 3 --cells-k ... --warp ...` 选择。驱动器只会
清理含自身 `.wcns-release-matrix` 标记的既有目录；对未标记目录立即失败，防止误删用户数据。
每个 rank 使用不同输出目录，随后以 rank 列表首项为参考逐场比较。

`--profile`、`--reconstruction`、`--riemann` 和 `--steps` 选择冻结算法组合与固定推进步数；
`--reference-directory` 指向另一矩阵的输出目录时，验证器以 cell-center 物理坐标排序后比较，
因此可严格比较同一逻辑网格的“单 zone 运行时切分”和“原生多 zone”结果，而不要求输出
zone 布局相同。统计文件同时通过 `series-constant` 检查全程守恒漂移。

CTest 中固定三条串行 smoke（二维、二维原生多区周期、三维扭曲），MPI 构建另执行
1/2 rank 等价性 smoke。它们验证基础设施可用，不替代 O2--O6 冻结的中等网格发布矩阵。

`run_vortex_matrix.py` 生成 2×2 原生多 zone 双向周期方形网格，执行等熵涡解析误差、全局
守恒和多 rank 最终场比较。解析式显式使用由模板参考量导出的 Mach 数；`--resolutions`
给出逗号分隔网格序列，`--minimum-order` 约束每一对相邻网格的观测阶，
`--minimum-finest-order` 只约束最细两级网格的观测阶，`--finest-l1` 约束最细网格误差。
阶段 O6 按设计卡口使用后两项检查 64²→128² 和 128²；保留 `--minimum-order` 供要求
所有加密区间均进入渐近区的专项测试使用。

`run_shock_matrix.py --case sod|quadrant` 生成矩形 CGNS 网格并运行间断算例。Sod 检查精确
Riemann 解的密度/速度/压力 L1、接触与激波位置以及正性；四象限算例检查 x-y 对角交换下
的密度、压力和速度分量对称性。两种模式均可用 `--ranks` 和 `--mpi-exec` 比较串并行最终场，
所有工作目录继续遵守标记文件清理协议。

`run_viscous_matrix.py --case couette|conduction` 生成 x 周期、上下壁面的二维通道。Couette
模式可用 `--velocity-curvature` 在不改变壁面值的前提下加入非平衡速度扰动，要求程序通过
残差判据停止，再独立检查线性速度、压力和壁面剪切；线性导热检查 `T=1+y`、常压和由状态
方程确定的密度。两种模式均比较多 rank 最终场。

`run_source_matrix.py` 在双向周期网格上检验常守恒源项的逐点解析更新和全局积分随时间的
线性变化，默认绝对误差阈值为 `2e-12`。`run_manufactured_matrix.py` 驱动二维或三维非均匀
制造初场及制造源项，至少推进指定步数，检查最终场有限性、正性和多 rank 等价性；三维模式
同时启用非零 z 动量源、扭曲网格以及 CGNS/Tecplot 双格式流场输出。

`run_output_restart_matrix.py` 固定执行非平凡周期场的 100 步连续计算、40+60 步检查点续算、
1/2/4 rank 异分区恢复及仅最终场/全监测对照。事件组同时启用 `every_steps=10`、
`every_time=0.01`、重合显式时刻和初末输出，要求 CGNS/Tecplot 事件成对且无重复；检查点每
20 步更新滚动 `latest`。验证器逐值比较两种场格式，并从基本场独立重算动量、总能、温度、
声速、Mach、总焓、熵代理、黏度和 Jacobian。残差历史与全局统计只比较可复现数值列，明确
排除 wall time 及重启初始行没有物理意义的前一步 `dt`。
驱动器的 `--dimension 3 --cells-k ...` 模式改用三维双 zone、x 周期网格，启用非零 z 动量
源，并额外确认最终 `VelocityZ` 确实非零；用于证明三维状态、格式和异 rank 重启路径被实际
执行，而不是二维零分量的形式覆盖。

## 4. 失败矩阵

```text
python tools/run_failure_matrix.py \
  --run build/wcns_run \
  --generator build/wcns_generate_release_cgns \
  --template cases/config/failure_base.wcns.in \
  --work-dir build/release-failures --ranks 1
```

MPI 模式再提供 `--mpi-exec <mpiexec> --ranks 1,2,4`。驱动逐项覆盖未知/重复配置键、直接
输入派生量 Re/Ma、未知/重复字段、未知统计量、不可行分区、单向连接、截断 CGNS、输出
路径类型冲突、既有输出目录、重启签名不匹配和非有限数值状态。每个进程数使用独立的失败
输出生命周期，并用命令超时把潜在 MPI 死锁判为失败。

配置/网格/重启等启动失败必须返回 `1`；正常到达最大步数返回 `2`；非有限数值状态返回
专用退出码 `3` 并写明 `numerical_failure`。失败前产生的 `.tmp` 只允许作为诊断残留，矩阵
明确拒绝任何伪装成场、检查点或成功 manifest 的最终文件。输出目录冲突还会比较失败前后
目录的 SHA-256 汇总，保证拒绝运行不会改写既有结果。observer 或 I/O 失败由所有 rank
共同退出，并报告首个失败 rank 的原始异常文字。

自定义字段/统计量的成功注册和求值、重复注册、未知依赖及循环依赖属于进程内扩展接口，
由 `wcns.unit` 覆盖；生产配置只能选择已经由程序注册的量，因此未知选择同时进入上述正式
可执行程序矩阵。
