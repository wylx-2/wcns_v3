# 阶段 H 正式验收报告（候选 v1）

## 1. 验收结论

阶段 H 的正式自动验收通过；项目负责人于 2026-08-31 批准候选 v1，人工卡口 G8 已通过。阶段 H 没有实现阶段 I 的高阶度量，也没有实现阶段 J 的无粘 WCNS/实际源项求值。

- 被测实现提交：`80f4fffc790269ec71a8a8af6b69e4370cffa7e4`
- 阶段分支：`stage/h-core-data`
- 候选标签：`stage-h-candidate-v1`（指向包含本报告的提交）
- 验收日期：2026-08-31（Asia/Shanghai）

## 2. 环境

- 操作系统：Windows
- CMake：3.28.0
- C/C++ 编译器：MinGW-w64 GCC 8.1.0
- 构建类型：串行和 MPI 均为 `Release`
- MPI：Intel MPI Library 2021.10，Build 20230619
- CGNS：仓库内固定的 CGNS 4.4.0 源码包，静态链接

## 3. 正式验收命令与结果

```powershell
cmake --build build --config Release --parallel
cmake --build build-mpi --config Release --parallel
python tools/verify_algorithm_spec.py
ctest --test-dir build -C Release --output-on-failure
ctest --test-dir build-mpi -C Release --output-on-failure
git diff --check
git status --short
```

结果：

- 串行 Release 构建：通过；
- MPI Release 构建：通过；
- 阶段 G 精确规格检查：6/6 通过；
- 串行 CTest：6/6 通过；
- MPI CTest：9/9 通过，其中 2-rank halo、Euler solver 和 MPI runtime 均通过；
- `git diff --check`：通过；
- 开始验收时工作树：干净。

## 4. 阶段 H 专项结果

| 专项 | 验证内容 | 结果 |
|---|---|---|
| H1 参考量 | 五个参考量正且有限；拒绝输入 Re/Ma；派生 Re、Ma、动压与时间；确定性摘要/重启签名 | 通过 |
| H2 气体与状态 | 气体常数输入二选一；温度/压力 primitive 与守恒量往返；声速、焓；二维 `w=0` | 通过 |
| H3 数值阈值 | 状态、Jacobian、面积、重构阈值集中校验及非法路径 | 通过 |
| H4 拓扑合法域 | range 标签不可互换；cell/vertex/face 分离；物理边界边角 ghost 与二维 K ghost 被拒绝 | 通过 |
| H5 源项框架 | 默认关闭 registry 为空；矛盾、未知、重复及阶段 H 实际启用路径确定性失败 | 通过 |
| H6 兼容回归 | 旧 residual 接口与新配置接口的自由流 residual、一次 SSPRK3 后守恒量逐位相等 | 通过 |

CGNS 边界拓扑在本阶段同时消除了旧字段的双重含义：`adjacent_cell_range` 的上边界法向索引为 `N-1`，`boundary_face_range/shared_face_range` 的上边界法向索引为 `N`。串行和 MPI 的 CGNS、多块连接、halo 与 Euler 回归全部通过。

## 5. 公共卡口状态

| 卡口 | 状态 | 证据 |
|---|---|---|
| G0 范围冻结 | 通过 | `docs/stage-h-design.md` |
| G1 接口与合法域 | 通过 | 强类型 range、`FieldDomain/TopologyField` 及异常测试 |
| G2 构建与静态检查 | 通过 | 两套 Release 构建和 `git diff --check` |
| G3 单元与公式 | 通过 | 精确规格 6/6，新增专项均进入单元测试 |
| G4 串行回归 | 通过 | CTest 6/6 |
| G5 MPI 回归 | 通过 | CTest 9/9 |
| G6 数值验收 | 通过 | 解析派生/往返容差测试及新旧路径逐位等价 |
| G7 文档与可追溯性 | 通过条件 | 本报告、阶段提交、远程分支及候选标签需共同存在；以远程 Git 记录为准 |
| G8 人工放行 | 通过 | 项目负责人于 2026-08-31 明确批准进入阶段 I |

## 6. 已知限制与后续边界

1. 现有阶段 E 求解路径仍使用压力型 `PrimitiveState`；物理边界统一迁移到温度型 ghost 属于阶段 J。
2. `TopologyField` 已冻结合法域接口，但现有低阶存量字段不会在阶段 H 全量迁移；新增高阶算法必须使用该合法域接口。
3. 源项模型只定义配置与注册契约，不执行任何源项；阶段 J 实现前所有开启配置均失败。
4. 两套度量 profile、几何通信、共享面所有权和高阶几何验收属于阶段 I，当前候选没有提前实现。

## 7. 人工验收请求

项目负责人已批准 `stage-h-candidate-v1`。按流程合并 `main`、创建 `stage-h-accepted` 标签后建立阶段 I 分支。
