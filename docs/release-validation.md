# 发布算例基础设施

阶段 O 的数值验收只调用正式 `wcns_run`，不建立简化求解器。共同入口如下：

- `wcns_generate_release_cgns` 按总单元数、维数、x 向原生 zone 数、扭曲幅值和 x 周期
  开关生成确定性结构 CGNS；x 周期至少使用两个 zone，避免生成当前拓扑契约不支持的
  self-connectivity；
- `cases/config/*.wcns.in` 保存严格 schema v1 模板；驱动器只替换 `@...@` 占位符；
- `wcns_validate_release_case` 独立重读最终 CGNS，可检查有限性/正性、均匀解析值和两个
  串并行结果的逐场最大差；
- `tools/run_release_matrix.py` 建立隔离工作目录，生成网格和配置，调用正式程序及验证器，
  最终写出 `matrix-summary.json` 和每条命令的独立日志。

## 1. 网格生成器

```text
wcns_generate_release_cgns output.cgns dimension cells_i cells_j cells_k zones_i warp periodic_x
wcns_generate_release_cgns periodic-square output.cgns cells_i cells_j length
wcns_generate_release_cgns rectangle output.cgns cells_i cells_j zones_i length_x length_y periodic_x
```

`dimension` 为 2 或 3；二维要求 `cells_k=1`；`cells_i` 必须能被 `zones_i` 整除；
`abs(warp)<0.2`。非周期网格的 x 外边界名为 `left/right`，其余物理边界为
`bottom/top[/front/back]`；原生多区连接及周期连接均成对写入。扭曲使用跨 zone 连续的
解析坐标映射，因而同一参数总是产生逐位相同的网格。
`periodic-square` 生成 2×2 原生多 zone 且 x/y 双向平移周期的方形网格；`rectangle` 用于
Sod 薄域和四象限问题，可指定二维物理长度及 x 向原生 zone 数。

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
给出逗号分隔网格序列，`--minimum-order` 和 `--finest-l1` 分别约束观测阶与最细网格误差。

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
