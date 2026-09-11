# WCNS v1.1.0 版本说明

WCNS v1.1.0 是在 v1.0.0 数值路径上完成的稳健性、壁面工程量、输运验证和性能版本。
`v1.1.0-rc.1` 已于 2026-09-11 通过人工核验；发布前又完成一次无功能变更的源码清理和
本地全量回归。当前版本用于本机独立开发并进入 `main`，不创建对外发行页面，不上传公开
发行产物，代码许可证决定暂缓。

## 主要变化

### 强激波事后稳健化

- 新增默认关闭的 SSPRK 候选态物理容许性检查；候选状态在独立缓冲区验证，失败步不会污染
  已接受状态。
- 新增 troubled-cell 识别、与实际 profile stencil 一致的直接/转置离散支持传播、共享面
  owner 最大级别归并，以及“原方案→同重构 primitive→linear5→zero-order+Rusanov”的
  单调降阶梯子。
- 局部重算仍失败时从完整 `U^n` 回滚并缩小时间步。history 和 manifest 记录每级面数、
  troubled cell、局部重算、整步 retry、候选最小物理状态与 proposed/accepted dt。
- Case07 四块圆柱 Mach 5 完整运行到 `t=8`，作为强激波候选的人工和自动验收对象。

原理、算法、传播方向、伪代码与终止性见《算法补充》11.1--11.3，尤其是 11.2.3；实现卡口
见 [`v1.1.0/stage-q-acceptance.md`](v1.1.0/stage-q-acceptance.md)。

### 边界面与载荷

- 新增权威真实边界面数据流，可输出 `p_w,T_w,mu_w,Cp,Cf,q_wall` 以及压力、黏性和总牵引。
- 使用全局守恒边界权重积分力/力矩及 `Cd/Cl/Cm`，MPI 结果按稳定面键排序、查重并写固定
  schema；二维明确为单位展向长度。
- 无粘压力/载荷输出只计算请求量；高阶压力（或显式请求的温度）迹发生非正输出侧超调时，
  仅对诊断值退回最近内部真实单元正值，不改变求解状态、通量、时间步或重启签名。
- Case07 的 Re=20/40/100 层流圆柱结果已使用生产 face-based 载荷链路，不再依赖后处理近似
  作为权威值。

### 输运和验证

- schema 1 正式暴露 `constant|sutherland`、Prandtl 数、参考温度黏度比，以及 Kelvin 或
  无量纲 Sutherland 常数二选一表示；规范化参数进入摘要、manifest 和重启签名。
- 增加黏性制造解空间阶、SSPRK3 时间阶、解析 Couette/导热/Poiseuille、Sutherland 生产
  smoke 和 CFL 分类扫描。新增公式见《算法补充》11.5。

### 性能与并行

- 缓存 line operator 和固定 stencil 行，使用显式 Roe 左特征行，复用 SSPRK/solver/halo
  工作区及 MPI 消息缓冲，并将 root 场输出和检查点改为分块汇聚/写出。
- 保留原同步交换和归约控制语义；没有宣称尚未实现的计算通信重叠或并行 CGNS。
- Windows/MinGW/Intel MPI 固定协议中，热路径分配从 946,214 次降至 6 次；二维等熵涡和
  三维黏性工作负载中位数加速分别为 1.589x 和 7.942x，4-rank 强扩展效率为 84.17%。
  这些数字只适用于记录的机器和构建，不是跨平台性能承诺。

### 发布前源码清理

- 使用仓库级 `.clang-format`、`.editorconfig`、`.gitattributes` 和 `pyproject.toml` 统一
  C/C++、Python、行尾及基础编辑格式；格式规则随内部源码包保存。
- 内部源码包明确包含 Case07 的配置、网格、证据和脚本，且继续排除不属于本版本范围的
  NACA0012。
- 机械整理覆盖 157 个 C/C++ 文件和 32 个 Python 文件。C/C++ 清理前后词法记号等价，
  Python 清理前后抽象语法树等价。
- 删除未使用的 Python 导入/中间变量、冗余默认成员初始化和恒等条件，限制仅供 MPI 使用的
  常量只在 MPI 构建出现，并消除少量无必要参数复制；未改变配置、数值公式、离散路径、输出
  schema、重启签名或外部可观察结果。
- 清理验收、保留项及命令见
  [`v1.1.0/release-cleanup-acceptance.md`](v1.1.0/release-cleanup-acceptance.md)。

## 配置兼容性

配置 `schema_version` 仍为 1。旧 v1.0 配置不含以下新键时保持原默认数值路径：稳健化关闭、
`Pr=0.72`、常黏度 `mu/mu_ref=1`、边界面输出关闭。新增配置组为：

- `robustness.*`；
- `transport.*`；
- `output.boundary.*`；
- 可选 `algorithm.mdcd.disp` 和 `algorithm.mdcd.diss`。

解析器仍严格拒绝未知/重复键、非法枚举、非有限值和模型专属键冲突。完整键、默认值与示例见
[`user-manual.md`](user-manual.md) 和 [`../examples/full_case_template.wcns`](../examples/full_case_template.wcns)。

## 检查点兼容性

v1.1 检查点继续保存五个无量纲守恒场、网格/数值签名和推进状态，并允许在合法条件下改变
rank 数和叶块划分。缺省输运的 v1.0 检查点按固定的常黏度 `Pr=0.72,mu/mu_ref=1` 迁移；
显式改变输运、稳健化、MDCD、气体、参考量、边界或源项等会影响轨迹的配置后，签名必须
不匹配并拒绝恢复。检查点仍是 WCNS 内部格式，不承诺与第三方 CFD 软件互换。

## 验证与支持平台

P--T 阶段候选均通过各自冻结自动卡口并有可追溯验收报告。U 阶段从空目录执行 Windows
Release 串行、Intel MPI、多类发布算例、错误路径、安装/解包复现和 T 性能复验。清理后又从
两个新目录完成无警告的串行/MPI Release 构建，CTest 为 60/60 和 108/108，算法规格为 6/6。
Linux GCC/Clang、OpenMPI 及 ASan/UBSan 作业仍保留在 CI 配置中，但依据本机独立开发决定，
其远端结果不作为 v1.1.0 内部版本完成条件；未实际运行的平台不会标为通过或受支持。
内部源码包经双次生成逐字节一致性检查，并从全新解压目录通过独立构建和 60/60 串行 CTest。

## 已知限制和许可

本版本仍不包含低 Mach 预处理、湍流/转捩、化学反应、隐式推进、AMR、非共形接口、并行
CGNS、动态插件、稳定库 ABI 或涡量/Q 体场输出。完整边界见
[`known-limitations.md`](known-limitations.md)。NACA0012 不属于本版本验收范围。

WCNS 自有代码尚未获得对外再分发许可证，v1.1.0 仅作为本机/内部版本使用并暂缓对外发布；
随附 CGNS 4.4.0 遵循独立许可和
[`../THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md)。
