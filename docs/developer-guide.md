# WCNS 用户自定义开发指南

本文面向需要修改或扩展 WCNS 的开发者，对应当前 `0.1.0`/schema 1 源码。目标不是只告诉
读者“改哪个文件”，而是说明一次扩展必须穿过哪些数据、验证、并行、重启、输出和测试路径，
避免新增代码在串行小算例中可运行、到多块/MPI/重启时失效。

先按 [`user-manual.md`](user-manual.md) 完成串行与 MPI 构建，并阅读
[`算法补充.md`](../算法补充.md) 和 [`known-limitations.md`](known-limitations.md)。本项目目前
没有稳定的对外 ABI；“接口”是源码扩展点，不是无需重编译的动态插件接口。

## 1. 扩展前的不可破坏契约

任何算法或模型扩展都必须保持以下约束，除非该任务明确包含一次经过设计和验收的契约升级。

1. 内部单元索引是 0-based 半开区间；CGNS 是 1-based，转换后必须检查上下界。
2. 二维仍使用三维索引容器和五分量守恒状态，但 k extent 为 1、w/z 动量为零。
3. WCNS 重构当前冻结为六点标量模板：首偏移 -2、点数 6、halo 宽度 3。
4. 块 ID、rank ID 和本地容器下标不是同一概念，不得互相替代。
5. 原生连接和运行时切分连接都要支持同 rank 与跨 rank；周期连接还要正确旋转矢量/通量。
6. 一个共享连接面的数值通量只由确定性 owner 计算，再通信给另一侧；接收侧不可独立重算。
7. `phenglei_wcns` 和 `scmm6_wcns` 是不可交叉的完整 profile bundle。
8. 物理边界 ghost 上只有边界条件构造的物理量、压力和守恒量有效；边/角 ghost、ghost
   坐标、度量、梯度及任意派生量不可读。
9. 密度、压力、温度、Jacobian、面面积和重构尺度必须经过统一 `NumericalFloors`/有限性检查；
   不允许用无记录的截断或零值掩盖错误。
10. 所有 rank 必须以相同顺序进入集合通信；某个 rank 的异常必须转成一致失败，不能让其他
    rank 停在 MPI 调用中。
11. 会改变数值轨迹或状态解释的配置必须进入摘要和重启签名；所有输入键进入配置 digest。
12. 新功能必须同时覆盖 2D/3D 适用性、物理边界、原生多块、运行时切分、串行/MPI 和重启影响。

## 2. 代码架构和一次时间步的数据流

生产入口在 `src/app/wcns_run.cpp`。启动和推进关系如下：

```text
.wcns 配置 --严格解析/广播--> CaseConfig
CGNS 元数据 -----------------> StructuredPartitionPlan
CGNS zone + 分区叶块 --------> StructuredMesh + LocalBlockSet
连接/所属 rank --------------> DistributedTopology / halo plan
原 zone 坐标 + profile ------> MetricField，再按叶块切片分发
初场或 checkpoint -----------> 每个真实单元的守恒场
边界配置 --------------------> BlockBoundaryDataMap
                                      |
                                      v
每个 SSPRK3 子步：更新原始量 -> 守恒 halo -> 物理 ghost
                -> 左右面重构 -> 物理面强约束 -> Riemann 通量
                -> 共享面通量通信 -> profile 通量散度
                -> 粘性项（可选）+ 源项 -> residual -> RK 更新
                                      |
                                      v
SimulationDriver -> 停止判据 -> OutputSchedule -> 场/历史/统计/checkpoint/manifest
```

主要数据对象：

| 对象 | 位置 | 职责 |
|---|---|---|
| `StructuredBlock` | `include/wcns/mesh/structured_block.hpp` | 坐标、度量、流场、边界和连接的本地块容器 |
| `StructuredMesh` | `include/wcns/mesh/structured_mesh.hpp` | 全局块描述及拓扑 |
| `LocalBlockSet` | `include/wcns/parallel/block_distribution.hpp` | 当前 rank 拥有的块 |
| `MetricField` | `include/wcns/mesh/high_order_metrics.hpp` | profile 绑定的单元坐标、J、三向面面积矢量 |
| `CaseConfig` | `include/wcns/runtime/case_config.hpp` | 严格配置的唯一结构化表示 |
| `InviscidWcnsConfig` | `include/wcns/solver/inviscid_wcns_solver.hpp` | 重构、Riemann、无粘边界、源项配置 |
| `SimulationDriver` | `include/wcns/runtime/simulation_driver.hpp` | SSPRK 外层步进、时间裁剪、停止和 observer |
| `RuntimeOutputManager` | `include/wcns/runtime/output_manager.hpp` | 输出事件、历史/统计、manifest 和原子提交 |

## 3. 推荐开发流程

### 3.1 建立可重复基线

```powershell
git status --short --branch
cmake -S . -B build-dev-serial -G "MinGW Makefiles" `
  -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_MPI=OFF `
  -DWCNS_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-dev-serial --parallel 4
ctest --test-dir build-dev-serial --output-on-failure
```

再建立独立 MPI 目录并运行 CTest。保存基线的 Git commit、测试结果和一个代表性算例 manifest。
已有未提交修改属于用户工作，新增功能不能覆盖或顺带格式化无关文件。

### 3.2 把扩展拆成小卡口

建议顺序：

1. 写数学/接口约定和适用范围；
2. 增加数据结构与严格验证，暂不接生产路径；
3. 写最小单元测试和失败测试；
4. 接入串行生产路径；
5. 接入多块/MPI、周期变换和共享面 owner 路径；
6. 更新配置摘要、restart signature、模板和文档；
7. 运行串行/MPI CTest；
8. 用正式 `wcns_run` 做具体算例，保存配置、日志、manifest 和独立验证结果；
9. 每个被认可的小卡口单独提交，提交中不要混入算例生成物或用户未审核修改。

### 3.3 完成定义

“代码能编译”不是完成。一个扩展至少应有：合法输入测试、非法输入测试、解析或制造结果、
自由流/守恒检查、2D/3D 适用性、物理边界、原生连接、运行时切分、1/多 rank 等价、输出选择、
重启兼容/拒绝、README/模板更新和 Release 构建实测。

## 4. 增加一个配置键

配置 parser 不是自动反射，新增键必须完整经过以下步骤。

1. 在合适的 config 结构中增加带安全默认值的成员，例如
   `include/wcns/runtime/case_config.hpp`。
2. 把完整键名加入 `src/runtime/case_config.cpp::fixed_keys()`。按 patch 动态命名的边界键需
   更新 `dynamic_boundary_key()` 属性集合。
3. 在 `CaseConfig::from_text()` 中用 `require` 或 `optional_*` 解析。枚举写专用解析函数；布尔
   只接受 `true|false`；实数必须有限。
4. 在对应 `validate()` 中检查数值范围、互斥/依赖条件、维数约束和功能开关状态。
5. 在 `summary()` 中输出用户可审查的最终值。
6. 如果该键影响离散算子、物理状态、边界、源项或恢复后轨迹，把它加入
   `CaseConfig::restart_signature()` 或相应子配置的 `restart_signature()`。
7. 把它传到真正消费该值的生产对象；只解析但不使用属于错误。
8. 在 `tests/test_case_config.cpp` 增加默认、合法、非法、未知/重复键和摘要/签名测试。
9. 更新 `examples/full_case_template.wcns`、用户手册和相关算例。

配置 digest 对规范化后的全部输入键值做 FNV-1a，用于各 rank 一致性，不等同于重启兼容
签名。输出调度、目录等可不进入重启签名；数值算法参数必须进入。

## 5. 增加或修改界面重构算法

### 5.1 公共接口

接口位于 `include/wcns/solver/wcns_reconstruction.hpp`：

```cpp
class IReconstructionScheme {
public:
    virtual ~IReconstructionScheme() = default;
    virtual std::string_view name() const noexcept = 0;
    virtual StencilRequirement stencil_requirement() const noexcept = 0;
    virtual Real reconstruct_scalar(
        ScalarStencilView stencil,
        TraceSide side,
        const ReconstructionContext& context) const = 0;
};
```

当前 registry 强制 `first_offset=-2,point_count=6,halo_width=3`。如果新方法需要其他模板宽度，
不能只改这个检查；还必须同时改块 ghost 分配、CGNS 叶块读取、物理 ghost、状态 halo、分区
最小宽度、面模板提取、测试和所有 profile 边界闭合。这应作为独立架构阶段处理。

### 5.2 六点模板算法接入步骤

1. 在 `src/solver/wcns_reconstruction.cpp` 或新的 solver 源文件实现
   `IReconstructionScheme`。算法名只能用小写字母、数字和下划线。
2. `TraceSide::Right` 必须与左侧镜像方向一致；可复用现有 `orient_stencil` 思路，不能简单用
   同一权重直接作用于未反转模板。
3. 所有输入先检查有限性；`context.scale` 和 `context.parameters` 用于无量纲 smoothness，
   不能绕过统一 epsilon/scale floor。
4. 在 `ReconstructionRegistry::with_builtins()` 注册工厂：

```cpp
result.register_scheme("my_scheme", [] {
    return std::make_unique<MyScheme>();
});
```

5. 当前 `reconstruct_thermodynamic_face()` 在内部创建 `with_builtins()`，所以只在应用层另建
   registry 不会生效。要让标准求解器识别新算法，必须把工厂加入该 built-in 注册路径，或先
   重构求解器使 registry 可注入。
6. `algorithm.reconstruction` 本身保存字符串，不需增加枚举即可选择注册名；若还维护
   `ReconstructionKind/reconstruction_name()`，则同步增加枚举分支。
7. 算法常数若要从配置输入，按第 4 章扩展 `WcnsParameters`、parser、summary 和签名。

### 5.3 正性与回退契约

重构结果必须通过守恒/原始态转换和 floor。特征重构失败时先回退同算法 primitive，再回退
`linear5`，最后一阶；不要在新算法内部悄悄夹断状态而不记录。若要改变回退链，应同步修改
`ReconstructionDiagnostics`、历史输出、重启签名和专项测试。

### 5.4 必需测试

- 常数保持、左右镜像、线性/高阶多项式精度；
- 光滑正弦网格收敛阶和间断非振荡性；
- conservative/primitive/characteristic 三种变量空间；
- 非法特征基和非物理解的回退事件、位置、计数；
- 扭曲网格自由流、物理边界面、原生多块和人工切分面；
- 1/2/4 rank 最终场与回退统计一致。

## 6. 增加 Riemann 求解器

### 6.1 实现接口

`include/wcns/solver/riemann_solver.hpp` 的 `IRiemannSolver::solve()` 接收左右压力原始态、单位
法向、气体和 floors，返回**单位面积**守恒通量与谱半径。面积由上层乘入。

```cpp
class MyRiemann final : public wcns::IRiemannSolver {
public:
    std::string_view name() const noexcept override { return "my_riemann"; }

    wcns::RiemannResult solve(
        const wcns::PressurePrimitiveState& left,
        const wcns::PressurePrimitiveState& right,
        wcns::Normal3 unit_normal,
        const wcns::GasModel& gas,
        const wcns::NumericalFloors& floors) const override;
};
```

返回结果必须满足：通量全部有限、谱半径有限且非负、`requested_solver` 为注册名、无回退时
`used_solver` 相同且 path 为空。有回退时必须给非 `None` 原因和连续的
`from_solver -> to_solver` 路径，最终到达 `used_solver`。

### 6.2 注册与配置

在 `RiemannSolverRegistry::with_builtins()` 注册。当前无粘和粘性 WCNS solver 构造函数都在
内部创建 built-in registry，所以标准入口必须修改这一注册路径，或重构构造函数接受注入的
registry。`algorithm.riemann` 是字符串，注册后即可由 schema 选择。

若有新参数：扩展 `RiemannSolverParameters`、验证、summary/restart signature 和配置 parser。
低 Mach 预处理不是简单开一个布尔值：需要定义预处理矩阵、适用 Mach、稳态/非定常一致性、
时间尺度和回退行为，并建立低 Mach 解析/收敛验收后才可暴露配置。

### 6.3 测试

至少覆盖相同状态物理通量、接触/激波、法向旋转不变性、二维 z 分量、熵修正、非物理中间态
回退、强壁面无穿透、周期旋转连接通量变换、自由流和串并行共享面唯一计算。

## 7. 增加算法 profile 或修改度量

仅增加重构/Riemann 不需要新 profile。profile 用于绑定一整套几何和线性算子；新 profile
必须作为不可拆 bundle 设计。

需要审查的入口至少包括：

- `include/wcns/mesh/algorithm_profile.hpp`：profile 与各组件 enum；
- `src/mesh/algorithm_profile.cpp`：期望组件、字符串工厂和禁止混用验证；
- `src/mesh/linear_operators.cpp`：插值、通量导数和单边闭合；
- `src/mesh/high_order_metrics.cpp`：单元坐标、Jacobian、面面积矢量和几何诊断；
- `src/mesh/conservation_weights.cpp`：残差/统计积分权重；
- `src/runtime/structured_partition.cpp`：该模板所需的叶块最小活动方向宽度；
- `src/app/wcns_run.cpp::initialize_partitioned_metrics()`：在原 zone 上计算后切片/打包/分发；
- 面通量 halo 层数、共享面 owner、metric pack/unpack 和 restart signature。

度量约定固定为计算步长 1，`J=partial(x,y,z)/partial(xi,eta,zeta)`，二维有面积、三维有体积
意义，必须为正。人工 MPI 切面不能变成几何单边界；高阶度量应先在完整原 CGNS zone 上计算，
再切片给运行时叶块。物理边界外没有可用坐标，不得通过伪造坐标把中心模板伸入 ghost。

验收包括笛卡尔精确值、解析扭曲网格收敛、离散几何守恒、正 J、2D 所有 z 分量初始化、3D
三方向、原生连接两侧一致、切分前后逐值一致、pack/unpack 无损和两套旧 profile 回归不变。

## 8. 增加源项模型

当前源项不是多态插件，而是 `SourceModelKind + SourceTermConfig + evaluate switch`。接入步骤：

1. 在 `include/wcns/physics/source_terms.hpp` 增加 enum 和参数存储；复杂参数可建独立结构。
2. 在 `source_model_name()`、`validate_kind()`、`SourceTermConfig::validate/summary/
   restart_signature()` 增加分支。
3. 在 `src/runtime/case_config.cpp` 增加允许键、字符串解析和数值读取。
4. 在 `SourceTermRegistry::evaluate()` 实现每单位无量纲物理体积的五分量守恒源。
5. 明确源是状态依赖、坐标依赖还是时间依赖；当前函数已提供 `U,x,t,dimension`。
6. 二维必须保证 z 动量源严格为零；能量源要与动量功一致。
7. 新参数和模型顺序进入摘要/重启签名，关闭开关时非零参数必须拒绝。
8. 若模型需要气体、参考量或其他场，扩展 registry 上下文，而不是读取全局变量。

测试源项自身的解析值、多个模型叠加、每个 SSPRK stage time、均匀场解析演化、守恒积分、
2D/3D、MPI 等价和 restart 连续性。若是制造源，必须把“初场/解析场”和“源”分别推导，不能
因为二者名字相同就假定自动匹配。

## 9. 增加或修改边界条件

边界扩展至少跨越以下层次：

1. `include/wcns/mesh/topology.hpp::BoundaryType` 增加类型；
2. `src/runtime/case_config.cpp` 增加字符串解析、名称输出、动态物理数据键和组合验证；
3. 如需读取 CGNS 原类型，更新 `src/io/cgns_reader.cpp` 的类型映射；注意生产入口最终仍以
   配置的 default/patch override 为准；
4. 扩展 `BoundaryPhysicalDataConfig` 和运行时 `BoundaryData`；
5. 在 `src/app/wcns_run.cpp::make_boundary_data()` 完成无量纲数据闭合和默认策略；
6. 在 `PhysicalGhostStateOperator` 的 `make_ghost()` 构造三层面 ghost 的
   `rho,u,v,w,T`，再由同一气体模型计算 p 和守恒量；
7. 在 `apply_inviscid_boundary_face_state()` 定义重构后的无粘真实面强/弱约束；
8. 在 `src/solver/viscous_boundary.cpp` 定义粘性真实面状态与梯度约束；
9. 把新增物理数据写入 summary 和 restart signature；
10. 更新配置模板、CGNS patch 名示例和边界测试。

不要给物理 boundary ghost 生成虚假坐标或度量。需要边界法向时使用真实面 `FaceMetric`；需要
高阶单边导数时只使用真实域可用数据和 profile 的单边闭合。边/角 ghost 当前不填，新算法若
读取它们必须先提出并实现清晰的多边界优先级，而不是依赖未初始化内存。

周期不是普通物理边界类型扩展。周期必须通过 `GridConnectivity1to1` 进入拓扑和 halo 路径，
处理坐标平移、旋转下的速度/动量/梯度/通量变换及反向连接一致性。

边界测试至少包含静止/移动滑移壁、无滑移绝热/等温壁、入流/远场亚声速与超声速特征选择、
出口有/无目标态、两侧法向符号、扭曲面法向、三层 ghost、无粘壁通量、粘性剪切/热流、
切分 patch 继承和 MPI。

## 10. 增加气体、输运或其他物理模型

### 10.1 当前气体模型

`GasModel` 只表示常 gamma 的热完全理想气体。所有温度/压力/守恒转换集中在
`thermodynamics.hpp/.cpp`。增加真实气体、多组分或变比热会影响状态变量、声速、焓、
Riemann、特征矩阵、时间步、边界、输出和 checkpoint，不适合只在一个函数中加 switch。

推荐先把 EOS/热力学改成明确的模型接口或封闭 variant，并保持转换函数为唯一入口。新增模型
必须有状态域验证、正性 floor、无量纲关系、声速/焓一致性和序列化签名。

### 10.2 当前输运模型

底层 `TransportConfig` 已支持常黏度和 Sutherland、Prandtl，但 schema 1 的生产入口只实例化
默认 `Pr=0.72,mu/mu_ref=1` 常黏度。若暴露配置：

1. 在 `CaseConfig` 增加 transport 配置与键；
2. 解析 `constant|sutherland` 及参数并验证正值；
3. 把**同一份** config 同时传给 `ViscousWcnsConfig::transport` 和输出
   `QuantityContext::transport`，否则求解黏度与输出黏度不一致；
4. 把 transport summary 纳入真正被 checkpoint 使用的 `CaseConfig::restart_signature()`；
5. 更新 Re/Pr/Sutherland 解析解、时间步稳定性和输出测试。

若增加湍流或多方程模型，还需扩展流场分量、halo payload、checkpoint 和输出，不应假装为
一个等效黏度标量完成。

## 11. 初始化新的算例流场

### 11.1 接入步骤

1. 在 `InitialConditionConfig::validate()` 的 `valid_types` 加入新名字。
2. 如有新参数，把 `initial.<name>` 加入 `fixed_keys()`；否则严格 parser 会在启动时拒绝。
3. 在验证中检查必需参数、正性、维数和参数组合，而不是全部依赖内建默认值。
4. 在 `src/runtime/flow_initializer.cpp` 写返回 `TemperaturePrimitiveState` 的纯函数。
5. 在 `FlowInitializer::evaluate()` 分派新类型。
6. 使用 `pressure_primitive()`、`temperature_primitive()` 和
   `thermodynamic_conservative()` 做闭合，不复制另一套状态方程。
7. `initialize_block()` 只遍历真实单元；不要初始化 ghost。halo 与物理边界在空间算子前填充。
8. 对周期解析场使用明确的周期距离，不能根据扭曲网格坐标包围盒猜周期。
9. 更新模板、用户手册和一个正式算例。

### 11.2 函数骨架

```cpp
TemperaturePrimitiveState my_initial_state(
    const InitialConditionConfig& config,
    Real x, Real y, Real z,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    int dimension)
{
    const Real rho = config.parameter("rho", 1.0);
    const Real u = /* analytic u(x,y,z) */;
    const Real v = /* analytic v(x,y,z) */;
    const Real w = dimension == 3 ? /* analytic w */ : 0.0;
    const Real temperature = /* positive nondimensional T */;
    TemperaturePrimitiveState state {rho, u, v, w, temperature};
    static_cast<void>(pressure_primitive(
        state, gas, reference, floors, dimension));
    return state;
}
```

测试坐标分区边界两侧、恰好位于间断/中心的点、默认参数、非法温度/压力、2D z 分量、多个
叶块初始化一致和解析输出。重启会跳过初场写入；因此修改初场不会改变已有 checkpoint 的状态。

## 12. 增加新的流场输出量

### 12.1 接口和描述符

实现 `IFieldQuantity`：

```cpp
class MyFieldQuantity final : public wcns::IFieldQuantity {
public:
    MyFieldQuantity()
    {
        descriptor_.name = "my_quantity";
        descriptor_.location = wcns::TopologyLocation::Cell;
        descriptor_.dimensional_unit = "1";
        descriptor_.scale = wcns::QuantityScale::Dimensionless;
    }

    const wcns::QuantityDescriptor& descriptor() const override
    {
        return descriptor_;
    }

    wcns::Real evaluate_cell(
        const wcns::StructuredBlock& block,
        const wcns::MetricField& metric,
        wcns::Index3 index,
        const wcns::QuantityContext& context) const override
    {
        // 只返回该真实单元上的有限无量纲值。
        return 0.0;
    }

private:
    wcns::QuantityDescriptor descriptor_;
};
```

当前只接受 `Cell` 位置。可选尺度为 Dimensionless、Density、Velocity、Pressure、Temperature、
Momentum、Energy、SpecificEnergy、Viscosity、LengthPower。量纲输出时 registry 自动乘尺度。

`descriptor.dependencies` 当前用于检查未知依赖、重复依赖和环，并不会把已计算依赖值传给
`evaluate_cell()`；实现仍需从 block/metric/context 自行计算，或以后扩展缓存接口。

### 12.2 注入标准入口

`ProductionFieldWriter` 支持传入 registry，但当前 `wcns_run` 使用默认参数。修改启动装配：

```cpp
auto field_registry = wcns::FieldQuantityRegistry::create_builtin();
field_registry.register_quantity(
    std::make_shared<MyFieldQuantity>());

wcns::ProductionFieldWriter field_writer(
    mpi, config, plan, local_blocks, metrics,
    quantity_context, mesh_name, std::move(field_registry));
```

然后配置 `output.field.quantities` 加入 `my_quantity`。还要把新 `.cpp` 加入 CMake 对应库。

### 12.3 梯度型二级量的特别要求

涡量、散度、Q、壁面剪切等不能在 `evaluate_cell()` 中无条件读取三方向 ghost：物理边界的
ghost 坐标/度量/梯度无效，输出时也没有独立的“刷新所有派生 halo”契约。可选设计是：

- 只用真实单元和 profile 单边闭合计算边界附近值；
- 在输出前增加明确的同步/precompute 阶段，把二级量写入有版本的 `TopologyField`；
- 对壁面量建立 face/boundary quantity 接口，而不是冒充 cell-centered 量。

任何方案都要测试原生接口和人工 MPI 切面的值一致，且不可把分区切面误当物理边界。

### 12.4 输出测试

测试唯一名称、重复注册、未知/循环依赖、有限性、无量纲/量纲尺度、CGNS 字段名、Tecplot 列、
多 zone 重组、串并行逐值一致和未知配置选择明确失败。

## 13. 增加新的全场统计量

实现 `IStatisticQuantity`：

```cpp
class MyStatistic final : public wcns::IStatisticQuantity {
public:
    MyStatistic()
    {
        descriptor_.name = "my_statistic";
        descriptor_.location = wcns::TopologyLocation::Cell;
        descriptor_.scale = wcns::QuantityScale::Dimensionless;
    }

    const wcns::QuantityDescriptor& descriptor() const override
    {
        return descriptor_;
    }

    wcns::Real evaluate(const wcns::StatisticContext& context) const override
    {
        wcns::Real local = 0.0;
        for (const auto& block : context.local_blocks.blocks()) {
            // 只累计本 rank 拥有的真实单元；需要积分时使用全局守恒权重和 J。
        }
        return context.mpi.sum(local);
    }

private:
    wcns::QuantityDescriptor descriptor_;
};
```

在应用层注入：

```cpp
auto statistics = wcns::StatisticRegistry::create_builtin();
statistics.register_quantity(std::make_shared<MyStatistic>());

wcns::RuntimeOutputManager output(
    mpi, config, plan, checkpoint.mesh_signature(),
    &statistic_context, event_writer, std::move(statistics));
```

每个统计事件所有 rank 都会调用 `evaluate()`，所以实现中集合通信的次数和顺序必须与名称
选择无关并在所有 rank 一致。只遍历 owned real cells，不能累计 ghost。体积积分使用
`GlobalConservationWeights` 与对应 `MetricField::jacobian()`，避免原 zone 物理边界闭合权重或
人工切分重复计数。

当前 registry 在 `output.dimensional=true` 时还会自动乘 descriptor 的场尺度和
`L_ref^dimension`。若统计量已经在实现中自行恢复量纲，descriptor 应保持 Dimensionless，
否则会二次缩放；更推荐实现始终返回无量纲值，由 descriptor 统一缩放。

测试解析积分、网格加密、原生多 zone/单 zone 切分等价、1/多 rank 归约、量纲尺度、重复名、
未知选择、输出列顺序和 restart 前后时间序列拼接。

## 14. 增加新的输出格式或观察量类别

字段格式位于 `ProductionFieldWriter::write()`，历史/统计和 manifest 位于
`RuntimeOutputManager`。新增格式需：

1. 扩展配置枚举、parser、名称、summary 和模板；
2. 明确 cell/face 位置、zone 顺序、坐标、单位和时间元数据；
3. 先写同目录 `.tmp`，完整关闭后再按 `allow_existing` 原子提交；
4. 只由确定 rank 写文件，但所有 rank 仍以一致顺序参加 gather/reduction；
5. 返回最终文件路径，让 manifest 记录；
6. 写独立重读器测试，不只检查“文件存在”；
7. 保持原 zone 重组，不把 MPI 叶块布局泄漏成用户物理网格，除非格式明确声明为分区输出。

若新增随时间统计类别，优先复用 `OutputSchedule` 和 observer；不要把 I/O 直接塞进求解器子步。

## 15. 修改停止机制或时间推进

停止判据在 `StopController`，时间步裁剪和 observer 通知在 `SimulationDriver`，SSPRK3 状态更新
在 solver/time integrator。新增停止条件必须定义：

- steady/unsteady 哪些模式生效；
- 与数值失败、正常完成、信号、墙钟、最大步数的优先级；
- 是否要求安全 checkpoint；
- manifest 名称和进程退出码；
- MPI 如何归约为全局一致决定；
- restart 需要保存哪些内部计数。

修改时间推进需重新核对源项 stage time、输出精确时间事件、t_end 不越界、checkpoint dt、
粘性稳定限制和残差定义。隐式、本地时间步或双时间步不能复用“非定常物理时间”语义而不
区分内外迭代。

## 16. 修改 CGNS 读取、连接或自动分区

### 16.1 CGNS 读取

新增 CGNS 特性时，元数据和实际块读取必须使用同一验证。所有 `cgsize_t` 转换到 `int` 前检查
范围；所有 Vertex PointRange 转换到 0-based 后检查 `vertex_extent`，再推导 cell/face range
并检查 `cell_extent`。不要信任供体名、transform 或边界定位。

### 16.2 连接

连接需验证接收/供体计数、轴置换、方向、互逆性和唯一覆盖。周期刚体变换必须是有限的正确
旋转/平移，并在状态 halo、梯度和面通量上用相符变换。消息 tag 必须确定且无冲突。

### 16.3 分区

自动分区当前是按单元数的确定性递归二分。改进负载模型时必须保持：同一元数据在所有 rank
产生相同 digest、叶块不重叠且完整覆盖原 zone、物理 patch 正确切片、跨原连接正确切片、
兄弟连接成对、每个活动方向满足 profile 下限、rank owner 稳定。

分区变化不应改变原 zone 上的度量、守恒积分和最终重组字段。用同一逻辑网格验证：单 zone
运行时切分、原生多 zone、不同 rank 数最终结果一致。

## 17. 检查点格式与兼容性开发

检查点版本在 `src/runtime/checkpoint.cpp`。新增持久状态时：

1. 明确它是否可由恢复后的守恒量重建；不可重建则写入文件；
2. 增加版本或向后兼容读取规则，绝不把缺失字段默认为看似合理的零；
3. 数值模型变化加入 restart signature；
4. 原 zone 存储，不绑定当前叶块/rank；
5. 保持无量纲守恒状态为恢复权威；
6. 测试 1→2、2→1、2→4 rank 和不同合法分区；
7. 测试网格、profile、气体、边界、源项、模型参数不匹配时在推进前失败；
8. 连续 N 步与 K 步 checkpoint + N-K 步恢复逐场比较。

当前 mesh signature 主要覆盖 base/zone 名、维数、extent 和坐标，不应把它当成所有 CGNS
拓扑语义的密码学证明。若扩展连接/边界语义，考虑把规范化拓扑也纳入签名。

## 18. 测试如何加入工程

### 18.1 普通单元测试

1. 在 `tests/test_<feature>.cpp` 写一个清晰命名的入口函数；
2. 在测试函数附近用简短注释说明它验证的目的、输入和主要失败条件；
3. 把源文件加入 `tests/CMakeLists.txt` 的 `wcns_unit_tests`；
4. 在 `tests/test_main.cpp` 声明并调用入口；
5. 保证测试没有外部顺序依赖、不会写源码目录、失败时抛出包含上下文的异常。

### 18.2 独立/MPI 测试

需要独立进程、CGNS fixture 或 MPI 的测试单独建 executable/add_test。MPI 测试至少运行 2 rank，
对可能死锁的失败路径设置合理 TIMEOUT。临时目录必须位于构建树；清理脚本只允许删除带自身
标记的目录。

### 18.3 生产入口算例测试

最终验收必须调用 `wcns_run --config`，不能用简化测试 solver 代替。驱动器应保存每条命令、
退出码、stdout/stderr、配置、网格、manifest 和独立验证结果。比较时排除 wall time 等天然
非确定列，但不能排除物理场、step/time、残差和回退计数。

### 18.4 常用命令

```powershell
cmake --build build-dev-serial --parallel 4
ctest --test-dir build-dev-serial -R wcns.unit --output-on-failure
ctest --test-dir build-dev-serial --output-on-failure

cmake --build build-dev-mpi --parallel 4
ctest --test-dir build-dev-mpi --output-on-failure
```

提交前再做全新 Release 构建，避免旧对象或 CMake cache 隐藏依赖遗漏。

## 19. 自定义扩展审查清单

### 配置和文档

- [ ] 名称、默认值、单位、范围和互斥条件已定义；未知/重复键仍失败。
- [ ] 生产代码实际消费该值；summary 和日志可看见最终选择。
- [ ] 影响数值轨迹的值进入 restart signature。
- [ ] 完整模板、用户手册、算法文档和已知限制同步。

### 数值和物理

- [ ] 推导与无量纲方程一致；2D/3D 分量和单位一致。
- [ ] 常数保持、解析值、精度/收敛和正性已测试。
- [ ] 回退、floor 或限制器全部可诊断，没有静默夹断。
- [ ] 物理边界只读取可用数据，边/角 ghost 未被意外访问。

### 多块和 MPI

- [ ] 原生连接、周期连接、人工切分面和同 rank/跨 rank 均覆盖。
- [ ] 共享面只计算一次，旋转/方向变换正确。
- [ ] 所有集合通信顺序一致，失败路径不会死锁。
- [ ] 1/多 rank 结果及诊断在规定容差内等价。

### 输出和重启

- [ ] 新字段/统计量可选择、量纲缩放正确、未知名失败。
- [ ] 文件经独立读回验证，原子提交和 manifest 记录完整。
- [ ] restart 连续性通过，不兼容配置在推进前拒绝。
- [ ] 不把 MPI 叶块布局错误地当成物理原 zone。

### 工程

- [ ] 串行/MPI、Debug/Release 中相关目标编译无警告升级。
- [ ] 单元、失败、端到端和人工算例均有记录。
- [ ] Git 提交只包含本小目标；没有覆盖用户未提交文件。

## 20. 常见开发错误

- 只把名字加入 parser，忘记 built-in registry，导致配置能读但 solver 报 unknown。
- 只在一个 rank 上注册统计量，其他 rank 进入不同集合通信顺序。
- 新重构声称需要更宽模板，却仍用三层 ghost 和旧分区下限。
- 在人工 MPI 切面两侧各算一次 Riemann 通量，破坏逐位守恒和串并行等价。
- 为物理 ghost 推造坐标/度量，让高阶中心差分越过真实边界。
- 给粘性 solver 改了 TransportConfig，却让 `viscosity` 输出仍使用默认模型。
- 新数值参数没有进入 restart signature，旧 checkpoint 被错误接受。
- 在输出 quantity 中读取“上次残留”的 halo 或派生场，没有版本和预输出同步。
- 用 `maximum_steps` 退出当作定常收敛，或只看 total L2 不看五分量判据。
- 测试只调用底层函数，没有通过正式配置、生产装配、输出和重启路径。

按照本指南完成扩展后，仍应以具体算例报告作为最终证据。报告必须公开说明对程序源码的所有
修改，不得把为算例添加的算法、边界、初场或验证逻辑隐藏在结果目录中。
