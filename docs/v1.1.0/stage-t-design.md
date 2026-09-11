# v1.1.0 阶段 T 设计冻结：热路径性能与并行扩展

状态：2026-09-11 在 R→S→T 连续推进授权下自动冻结；候选
`v1.1.0-t-candidate.1` 已于同日获项目负责人人工核验通过，允许进入 U。

## 1. 基线、剖析和范围

固定机器、Release 编译器和三项 P 基线保持不变：100x100 二维等熵涡 4.5033814 s、Case07
四块 48x32 圆柱 0.5177693 s、36x48x36 三维黏性算例 67.9583291 s；均为一次预热、五次
正式重复且 CV 小于 5%。T 只接受相同输入、步数、输出开关和工具链的中位数对照。

额外的 `-pg` 采样构建只用于热点排序，不用于绝对加速比。三维黏性样本中：

- `LineOperators::build` 被调用 28,345,034 次，占 30.28% 自身采样时间；
- `compute_primitive_gradients` 占 11.83%；
- 源码审计确认 `LineOperators::build` 位于真实单元/变量/方向循环内；
- 每次残差还重新构造无粘面通量、梯度操作数、梯度、黏性面通量、四类 halo plan、registry
  和 MPI 临时消息缓冲。

因此优化顺序冻结为：先消除重复离散算子装配，再持久化 plan/字段/registry 和 SSPRK 状态，
最后持久化消息缓冲及分块 root I/O。不得把重构、通量、壁面闭合、归约顺序或输出精度变化
混入优化提交。

采样还显示每个 characteristic face 都对固定结构的 5x5 Roe 右特征矩阵执行通用 Gauss--Jordan
求逆。允许按《算法补充》4.1 新增的显式左特征行替换该求逆，但必须保留右特征列顺序、局部
法切基、`L*R=I` 运行时检查及投影/恢复对照；这是独立可回退的代数等价优化。

本阶段不增加通信压缩、GPU、OpenMP、并行 CGNS/HDF5、异步 I/O 或新数值算法。NACA0012
用户目录不在本阶段范围内。

## 2. 不变量与生命周期

热路径缓存键为《算法补充》11.6 的
`K_workspace=(mesh_signature,partition_digest,profile,dimension,local_layout)`。初始化时为每个本地
块和有效轴构造一次 `LineOperators`，为当前全局拓扑构造一次 face/operand/gradient/viscous
halo plan，并创建固定 block-id registry。每个字段持有的 profile、维数、extent 和 halo 宽度
必须与该键一致。

运行状态 `version` 不属于结构缓存键。每个残差调用先单调递增版本，再把同一版本写入所有
工作字段和 plan 消息头；消费者仍逐项验证 profile 与版本。工作字段在改写前用 NaN 重置存储，
使漏写、漏收或提前读取不能被上一步的有限值掩盖。版本溢出、网格/布局不匹配或错误重启
必须立即失败，不能隐式重新解释旧存储。

SSPRK 初始真实单元状态按本地块大小只扩容一次，之后覆盖复用。稳健模式的整步快照仍服从
Q 阶段回滚契约；优化不能使被拒绝步污染接受状态。

## 3. 接口和任务拆分

### T1：离散算子与拓扑计划

- 给残差装配和梯度算子提供已构造的轴向 `LineOperators`；公开单点 stencil API 保持兼容。
- 离散支持传播保留公开 stencil 查询；生产残差按同一 profile 行表引用，避免每单元构造
  `StencilRow`。WENO-Z 的默认平方幂使用乘法求值，块级不变量只验证一次，逐面有限性、权重
  和正性检查保留。
- halo plan 增加受校验的 `set_version`/等价接口，只更新描述符版本，不改变 pairs、owner、tag、
  消息长度或周期变换。
- 增加测试证明同一缓存跨版本的 pairs/tag/owner 不变、旧版本字段被拒绝、错误 profile/extent
  被拒绝。

### T2：持久字段、SSPRK 状态和 registry

- 无粘 solver 持久化面通量、registry、face plan 和本地 block 指针。
- 黏性 solver 另持久化梯度操作数、真实梯度、黏性面通量及其 registry/plan。
- 计算核增加写入既有字段的 `*_into` 入口；原返回值 API 作为薄包装保留给现有测试和独立调用。
- `advance_ssprk3` 增加显式 workspace 入口，现有 API 保持兼容并只在独立调用时临时创建。

### T3：通信缓冲和重叠决策

- state/face/operand/gradient/viscous exchanger 在初始化时按 descriptor 精确长度建立 send/receive
  缓冲和 request 槽，调用时仅 pack、发布、等待、校验版本并 unpack。
- 本阶段保留同步 `exchange` 语义。现有计算核以整个 block 为粒度，边界附近 stencil 与 owner
  通量会影响随后全块残差；在没有内部/连接邻域的可验证分区前，不发布伪 begin/finish API，
  也不声称通信计算重叠。
- NaN/旧版本/错误长度及 1/2/4/8-rank 注入测试证明等待前不消费 receive buffer。

### T4：集合通信决策

时间步最小值、推进成功、残差范数、停止信号和 observer 成功分别位于不同的控制依赖点；
强行合并会改变失败与停止语义。本阶段不合并这些归约。性能报告记录每步归约次数和载荷，
把该项作为“审计后保留”的显式结论，而不是未完成实现。

### T5：ADF/root 分块输出

字段与检查点在 root 上按 source zone，再按 quantity/component 分块汇聚和立即写入；任一时刻
只允许保留当前块、CGNS 必需坐标和固定元数据。最终仍由 root 使用 ADF 串行原子提交，不称为
并行 CGNS。独立重读必须与内存值及旧 schema 一致。

rank 拼接次序到原 zone 线性次序使用 partition 定义的置换在 gathered payload 内原地完成；
CGNS 和 checkpoint 每写完一个 field 即释放该 payload。Tecplot 使用标准 BLOCK 布局逐变量写入，
校验器同时保留旧 POINT 布局的读取兼容。root 异常在下一次集合调用前广播给所有 rank。

### T6/T7：计时和性能证据

`tools/run_stage_t_performance.py` 生成可选 JSON manifest，记录提交、机器、编译器、命令、预热、
五次样本、median/min/max/CV、RSS、cell-stage/s、分配探针、1/2/4/8-rank 强扩展、等负载弱
扩展以及输出/checkpoint 字节率。`--detailed` 必须在控制台和 JSON 写
`detailed_timing=true`；默认模式不启用采样器或逐面计时。

正式 Windows/Intel MPI 扩展性测量使用 5 个时间步、一次预热和五次正式样本。对于本机
12 逻辑处理器的混合核拓扑，处理器列表固定为 `0,2,4,6,8,9,10,11`，rank 从左至右取用，
并把本轮新建的求解进程设为 High priority；两项设置必须写入 manifest。该策略只约束测量
环境，不改变离散、网格、分区或通信路径。2 步样本只有约 3 s，单次系统抢占足以使 CV 超过
5%；5 步样本仍使用同一问题且所有 rank 工作量一致，用于使系统噪声占比可验证地降低。

默认计时由主线程直接等待子进程，不枚举进程树；`--detailed` 的 RSS 采样在独立线程完成，
进程 `wait` 返回的时刻立即截取墙钟终点，采样线程停止和日志关闭不计入求解墙钟。默认与
detailed 的正式开销对照在五次预热后各保留五个样本，防止混合核电源状态的冷启动漂移。
复用已有 manifest 时，计时模式、预热/重复次数、扩展步数、绑核和优先级协议必须完全一致。

## 4. 分配探针和性能公式

单独的串行 allocation probe 覆盖一次黏性 RK 残差。探针只在被测区间打开全局
`new/new[]` 计数，并在关闭后输出次数和字节；初始化、日志和测试框架分配不计入。S 候选上的
计数作为 `A_hot,baseline`，T 候选使用同一探针得到 `A_hot,new`，按《算法补充》11.6 计算
减少比例。候选必须至少减少 90%，且连续第二次 residual 不得因缓存扩容增加分配。

## 5. 冻结自动卡口

| 卡口 | 自动证据 | 阈值 |
|---|---|---|
| T-V1 数值等价 | P--S 全量、基线签名及输出重读 | 串行同配置逐位一致；MPI 既有冻结容差不放宽 |
| T-V2 分配 | 同一 allocation probe 的 S/T 对照 | 每 RK 残差成功分配次数减少 `>=90%`，第二次不回升 |
| T-V3 串行吞吐 | 三项 P 固定工作负载，一次预热、五次样本 | 每项 CV `<=5%`；代表性纯计算中位数加速 `>=1.5x` |
| T-V4 计时开销 | 默认与 detailed 各五次 | detailed 关闭时相对无计时探针开销 `<=2%` |
| T-V5 强扩展 | 足够大三维问题，1/2/4/8 rank | 每 rank 至少 15,552 真实单元；4-rank 效率 `>=70%` |
| T-V6 弱扩展 | 每 rank 固定局部 cell 数，1/2/4/8 rank | 保存效率和通信分解；不设掩盖单机内存带宽的虚假硬阈值 |
| T-V7 I/O | 场输出/checkpoint、分块峰值及独立重读 | 数值/schema 一致；root 暂存不超过最大单 zone/quantity payload 加固定元数据 |
| T-V8 安全 | 版本、错误长度、NaN、重启/网格签名、错误注入 | 1/2/4/8 rank 无死锁且一致失败 |
| T-V9 回归 | Release 串行/MPI CTest、规格、diff | 全部通过；仅保留已登记用户目录 |

代表性串行目标解释为三项 P 工作负载均报告，其中二维等熵涡和三维黏性两项均须达到 1.5x；
Case07 20 步短算固定启动成本较高，仍必须无性能回退，但不单独用其总墙钟否决计算核目标。
4-rank 强扩展使用三维工作负载且每 rank 不少于 15,552 单元。性能失败不能通过验收后改阈值；
必须修复或保留失败候选并停止。

## 6. 提交和停止条件

设计/探针、离散算子缓存、持久 solver 工作区、持久通信缓冲、分块 I/O、性能工具和正式验收
分别独立提交。每个优化提交运行最小等价测试和 `git diff --check`；阶段末执行完整串行/MPI
回归和五重复性能矩阵。全部通过后创建不可移动的 `v1.1.0-t-candidate.N`，但不得合入或进入
U；在 T→U 人工卡口等待项目负责人判断。
