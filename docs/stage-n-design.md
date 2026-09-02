# 阶段 N 设计：运行控制、输出监测与检查点重启

状态：N0 冻结规格。

本阶段把阶段 M 的短步生产入口扩展为可停止、可观测、可重读、可重启的唯一正式运行链。
数学定义和停止优先级以 [`release-development-plan.md`](release-development-plan.md) 第 4--5
节为准；本文冻结首版接口、配置键和小任务提交边界。

## 1. 唯一运行状态机

`SimulationDriver` 是生产推进的唯一循环。求解器通过小型 `ISimulationSolver` 适配器向它提供
`global_time_step`、`advance`、当前残差和诊断计数；无粘/黏性分支不得各自复制停止逻辑。

每个完整步严格执行：

```text
计算全局 dt -> 裁剪到 t_end/下一精确时间事件 -> 完整 SSPRK3 步
-> step/time 更新 -> 全局残差与有限性检查 -> StopController
-> 历史/统计/流场/检查点调度 -> 若停止则唯一最终输出和 manifest
```

`StopReason` 为 `running | steady_converged | physical_time_reached | maximum_steps |
wall_time_checkpoint | user_signal_checkpoint | numerical_failure`。同一步优先级为数值失败、模式正常
完成、安全停止、最大步数。所有判定输入先做 MPI 归约，因此所有 rank 得到相同原因。

定常模式按五个守恒分量分别计算真实单元、正 `J` 加权的全局 `L2/Linf`；连续检查只在
`check_interval_steps` 发生，参考值首次检查冻结并可检查点恢复。非定常模式的残差不参与
停止，时间步精确裁剪至 `t_end`，`max_steps` 对两种模式均为硬上限。

## 2. schema 版本 1 的阶段 N 扩展

仍使用阶段 M 的严格 parser；未知、重复、缺失和非有限值启动即失败。增加以下键：

- 运行：`run.mode`、`run.t_end`、`run.max_wall_time`；
- 定常判停：`steady.min_steps`、`steady.check_interval_steps`、
  `steady.consecutive_checks`、`steady.reference_floor`、`steady.l2_absolute`、
  `steady.l2_relative`、`steady.linf_enabled`、`steady.linf_absolute`、
  `steady.linf_relative`；
- 输出根：`output.directory`、`output.allow_existing`、`output.dimensional`；
- 流场：`output.field.enabled`、`output.field.format = cgns | tecplot | both`、
  五个通用调度键及 `output.field.quantities`；
- 历史：`output.history.enabled`、`output.history.format = txt | tecplot`、
  `output.history.every_steps`；
- 统计：`output.statistics.enabled`、`output.statistics.format`、
  `output.statistics.quantities` 及通用调度键；
- 检查点：`output.checkpoint.enabled`、通用调度键、`restart.path`。

通用调度键为 `every_steps`、`every_time`、逗号分隔 `explicit_times`、`write_initial`、
`write_final`。零间隔表示关闭该类周期事件；显式时刻必须严格递增且非负。`restart.path` 为空
表示从配置初场启动。初场和重启互斥，重启数据是新的初始状态。

## 3. 输出数据模型

`OutputSchedule` 只生成事件，不执行 I/O；步事件、时间事件、显式事件、初场和最终场取并集，
以 `(kind,step,time)` 容差去重。非定常驱动向所有启用 schedule 查询下一个精确时间事件，统一
裁剪时间步。

`OriginalZoneField` 使用原 CGNS zone 的 cell extent 和固定 I-fastest 存储。每个 owner 将叶块
真实单元按 `source_zone + source cell range` 打包；I/O rank 重组时使用覆盖位图，遗漏、重复或
越界均失败。不得把当前叶块编号或 rank 分片写成重启语义。

首版流场内建量为 `rho,u,v,w,p,T,rho_u,rho_v,rho_w,rho_E`，并提供声速、Mach、总焓、
熵代理、黏度和 Jacobian；需要梯度的涡量、散度、Q、热流在依赖合法时注册，否则请求即明确
失败。`IFieldQuantity`/`IStatisticQuantity` 声明唯一名称、cell/face 位置、依赖、单位和归约，
registry 拒绝未知、重名及循环依赖。

CGNS 输出按原 zone 写 `CellCenter` FlowSolution；Tecplot 为每个原 zone 的 cell-center ordered
zone。历史与统计只由 rank 0 写固定列 TXT/Tecplot。文件先写同目录临时名，关闭成功后原子
改名；默认拒绝既有输出根，续算只有显式 `allow_existing=true` 才可追加新 step 文件。

## 4. 检查点与重启

检查点使用项目内部 CGNS 布局，按原 zone 保存五个守恒量，并保存格式版本、完整步 step/time/
dt、配置摘要、网格和分区无关签名、算法/参考量签名、定常参考残差与连续通过次数。文件仅在
完整步边界产生。

重启先按原 zone 读取，再根据本次新建的 `StructuredPartitionPlan` 以全局 cell range 恢复到
owner 叶块，因而支持不同合法 rank 数和不同切分。网格维数/zone extent、算法 profile、气体、
参考量、边界、源项或格式版本不兼容时，在第一次 halo 通信前全 rank 一致失败。输出目录和
调度可以改变，不属于数值重启签名。

## 5. 小任务与提交边界

1. N0：冻结本文和完整配置示例；
2. N1：扩展 `CaseConfig`，实现 `ResidualNorms`、`StopController` 和停止原因；
3. N2：实现唯一 `SimulationDriver`、时间事件裁剪和无粘/黏性适配器；
4. N3：实现 `OutputSchedule`、目录策略、历史和 manifest；
5. N4：实现全场量/统计量 registry 与内建基础量；
6. N5：实现原 zone MPI 重组、CGNS/Tecplot 流场和 TXT/Tecplot 序列；
7. N6：实现 CGNS 检查点、严格签名及同/异 rank 重启；
8. N7：完成停止、调度去重、输出无扰动、I/O 失败和重启连续性自动测试；
9. N8：完成配置、运行、输出、扩展和故障诊断文档，执行阶段 N 正式自动验收。

每个小任务达到局部测试后独立提交并尝试推送。N8 自动卡口通过后创建
`stage-n-candidate-v1`，不等待人工回复，直接进入阶段 O。

## 6. 非目标与首版限制

- 不实现并行 HDF5；首版由 rank 0 重组单文件，字段语义不因此改变；
- 不把 Tecplot 当作检查点，不从可视化舍入文本重启；
- 不在半个 RK 步或部分 rank 上写可重启文件；
- 不用残差提前终止非定常计算，不用 `t_end` 终止定常伪时间计算；
- 不允许输出回调修改求解场，关闭/开启输出后的内存状态必须一致；
- 信号处理仅设置无锁安全标志，真正 MPI 通信和 I/O 留到完整步边界。
