# 阶段 H 设计：核心热力学、合法域与源项框架

状态：已获项目负责人批准，允许在 `stage/h-core-data` 分支实施。

## 1. 阶段目标

阶段 H 建立后续高阶度量、边界和源项依赖的基础设施，同时保证默认配置下现有阶段 E Euler 求解结果逐位不变。本阶段不实现两套高阶度量、不替换现有 WCNS 空间离散、不实现粘性通量，也不启用实际源项。

## 2. 当前代码兼容边界

现有 `PrimitiveState` 的第五分量是压力，且阶段 E 的重构、边界和 Riemann 代码依赖此布局。`算法补充.md` 为后续物理边界定义的权威 primitive ghost 是 $(\rho,u,v,w,T)$。阶段 H 采用以下过渡方案：

- 保留现有压力型 `PrimitiveState` 和 `IdealGas` 接口，确保旧求解路径不变；
- 新增明确命名的 `TemperaturePrimitiveState` 和 `GasModel`，负责温度型状态、压力型状态与守恒量之间的统一转换；
- 阶段 H 只验证新旧接口在相同无量纲状态下给出一致守恒量；
- 物理边界存储正式迁移到公共温度型 ghost、同步派生 $p/U$ 的工作属于阶段 J，不在 H 中提前实施。

仓库目前没有通用算例配置文件解析器、日志器或检查点写入器。本阶段以强类型输入结构表示配置，并提供确定性的 log/restart 摘要字符串；未来接入具体输入格式时必须调用相同验证入口，不能复制验证逻辑。

## 3. 小任务与文件规划

### H1：参考量与气体模型

- 新增 `include/wcns/physics/thermodynamics.hpp` 和对应实现。
- `ReferenceInput` 只接受五个参考量，并保留专门的禁用 `Re/Ma` 字段检测入口。
- `ReferenceScales::derive` 计算 $Re$、$Ma$、参考动压和时间尺度，提供固定顺序摘要。
- `GasModelInput` 严格执行 `molar_mass`/`specific_gas_constant` 二选一。
- `GasModel` 实现温度型 primitive、压力型 primitive、守恒量、$p,T,a,h$ 的换算。

### H2：数值阈值与重构尺度

- 实现集中式 `NumericalFloors`，包含状态、Jacobian、面积和重构阈值。
- 实现 `ReconstructionScaling` 基础容器和按分量归一化接口。
- 本阶段不改变现有 WCNS 权重公式；阶段 J 接入新尺度前保持旧 `WcnsParameters` 行为。

### H3：强类型范围和字段合法域

- 新增 `TopologyLocation`、`AccessRegion`、强类型范围包装和 `TopologyField`。
- 物理边界 slab 只允许一个活动方向越出真实区，其他活动方向必须是真实索引；二维非活动 K 方向必须保持 0。
- 连接 halo 与物理边界 slab 必须由调用者显式区分。
- 现有低层 `Field` 作为存储实现保留；新算法只能通过带策略的 `TopologyField` 访问需要严格合法域的数据。
- 阶段 H 将拓扑描述符中的边界/连接范围迁移为强类型，编译期禁止不同范围标签直接替换。

### H4：源项空框架

- 新增 `SourceTermConfig`、`SourceModelKind` 和 `SourceTermRegistry`。
- 默认关闭且模型列表为空；关闭路径不分配模型、不执行源项。
- 开启但列表为空、关闭但列表非空、未知模型或当前阶段尝试创建实际模型都必须失败。
- `SpatialParameters` 持有默认关闭的源项配置并参与验证，但 `compute_euler_residual` 在 H 阶段不增加任何源项运算。

### H5：兼容接入与回归

- 新增 H 专项单元测试并加入 `wcns_unit_tests`。
- 记录阶段 E 代表状态转换、自由流残差和一个 SSPRK3 步结果，验证默认新配置前后逐位一致。
- 运行精确规格检查、串行全部测试和 MPI 全部测试。

## 4. 接口约束

1. 所有配置值在构造/派生边界一次性验证；对象创建后只读。
2. `Re`、`Ma` 没有公共 setter，也不能从输入覆盖。
3. `GasModel` 与 `ReferenceScales` 必须显式配对；温度换算不得从全局变量读取 Mach 数。
4. 二维温度型状态要求 $w=0$；转换失败抛出明确的物理配置异常。
5. `TopologyField` 不允许通过裸 `Field&` 绕过合法域检查。
6. 源项关闭路径不进入虚函数、注册表循环或 residual 写入。

## 5. 阶段 H 专项验收

- 五个参考量为正且有限；输入携带 `Re` 或 `Ma` 时失败。
- 气体常数输入严格二选一；派生 $\widetilde R,Re,Ma$ 与解析值一致。
- 温度型 primitive、压力型 primitive 和守恒量往返达到舍入误差，二维非零 $w$ 失败。
- 所有 `NumericalFloors` 和 `ReconstructionScaling` 非法值路径有测试。
- 强类型范围不可隐式互换；物理边界角 ghost、非活动轴 ghost 和错误拓扑位置访问失败。
- `enable_source_terms=false` 的 registry 为空；所有不允许的开启组合失败。
- 阶段 G 规格检查、串行测试、MPI 测试全部通过，默认 Euler 代表结果逐位不变。

## 6. 明确不在阶段 H 范围内

- CGNS 高阶度量 profile 和几何分阶段通信；
- 公共温度型物理 ghost 对现有边界求解器的替换；
- `strong_boundary_face_state` 新实现；
- WCNS profile 专属插值/差分、`FaceFluxHalo`；
- 任何实际源项求值；
- 粘性梯度、输运模型和壁面迹。
