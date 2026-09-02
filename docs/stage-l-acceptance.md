# 阶段 L 正式验收报告（候选 v1）

## 1. 验收结论

阶段 L 候选 v1 的必做范围 L1--L5、L7 已完成，正式自动验收通过，现停在人工检查卡口。
可选 L6 `weiss_smith_roe` 没有获得冻结文档要求的单独批准，因此没有实现，也没有提供可开启
但数学不完整的占位模式。阶段 K 的人工检查仍为暂缓状态；本次结果不追认 K，也不授权进入
阶段 M 或合并 `main`。

- 被测实现提交：`952dff6`；
- 阶段分支：`stage/l-algorithm-extensions`；
- 候选标签：`stage-l-candidate-v1`（指向包含本报告的提交）；
- 验收日期：2026-09-02（Asia/Shanghai）。

## 2. 完成内容

| 小任务 | 完成内容 | 主要提交 |
|---|---|---|
| L0 | 冻结四种重构、三种 Riemann、回退和扩展接口 | `fc6d00b` |
| L1 | `IReconstructionScheme`/`IRiemannSolver` 注册表、自定义策略、严格名称与配置 | `3141326` |
| L2 | `weno_js`、`weno_z`、`mdcd_linear`、`mdcd_hybrid`，修正 MDCD 指标首项 `271779` | `c1f12c0` |
| L3 | conservative、primitive、characteristic 六点同面 Roe 特征基重构及确定性降阶 | `340658d` |
| L4 | Rusanov、Batten 波速 HLLC、带 Harten 熵修正 Roe，以及 `Roe -> HLLC -> Rusanov` 回退 | `3e6d04d` |
| L5 | 可定位降阶/回退事件、逐级 Riemann 回退路径、MPI 全局计数、共享面唯一 owner 计算 | `9861b07` |
| L6 | 可选且需单独批准；本候选未实施 | 不适用 |
| L7 | 光滑、镜像、接触、激波、低压、高 Mach、单块/多块/MPI 及 H--K 回归 | `952dff6` |

重构与 Riemann 都由注册表创建。新增算法只需实现相应接口、注册唯一键、定义参数校验与
重启签名并补齐冻结文档列出的测试，不需要修改通用面通量调用链。

## 3. 正式验收环境与命令

- Windows，CMake 3.28.0，MinGW-w64 GCC 8.1.0；
- 串行和 Intel MPI 2021.10 Release 构建；
- 仓库内固定 CGNS 4.4.0 静态库。

```powershell
cmake --build build --config Release --parallel 4
cmake --build build-mpi --config Release --parallel 4
python tools/verify_algorithm_spec.py
ctest --test-dir build -C Release --output-on-failure
ctest --test-dir build-mpi -C Release --output-on-failure
git diff --check
git status --short
```

结果：两套构建通过；阶段 G 规格检查 6/6 通过；串行 CTest 9/9 通过；MPI CTest 15/15
通过，其中包含双 rank 多块阶段 J/K 回归以及阶段 L 三种 Riemann 与 characteristic WENO-Z
路径；`git diff --check` 通过，正式验收开始时工作树干净。

## 4. 数值与健壮性专项

| 专项 | 验证内容 | 结果 |
|---|---|---|
| 四种重构 | 常数保持、尺度不变、左右镜像、光滑特征变量加密误差比大于 12 | 通过 |
| MDCD | LINEAR 系数、HYBRID 传感器/权、修正后的六点指标、50 组确定性样本非负 | 通过 |
| 非法输入 | 错误模板、NaN/Inf、未知名称、重复注册和不兼容参数确定性失败 | 通过 |
| 特征重构 | 同面 Roe 基、确定切向、`L_n R_n=I`、往返变换、均匀流和完整降阶链 | 通过 |
| Riemann | 三种求解器一致性和法向反转；HLLC 静止接触与超声速迎风；Roe 单特征波和熵修正 | 通过 |
| 健壮状态 | Sod 激波、高 Mach、近 floor HLLC 回退以及两级回退事件逐级计数 | 通过 |
| 共享面/MPI | 128 个单元的双块网格只计算 280 个唯一二维面；串行与双 rank 计数、残差一致 | 通过 |
| 全量回归 | CGNS、halo、Euler、两套度量 profile、阶段 J 无粘、阶段 K 黏性均通过 | 通过 |

降阶/回退事件记录请求算法、前后策略、原因、轴向、面索引、块、rank、残差求值版本和
RK 子步。非 owner 的块间连接面不再先计算后覆盖，而是直接等待 `FaceFluxHaloExchanger`
发布 owner 的唯一面通量；这同时避免重复回退计数。

## 5. 卡口状态

| 卡口 | 状态 | 证据 |
|---|---|---|
| L0 设计冻结 | 通过 | `docs/stage-l-design.md` |
| L1--L4 算法与接口 | 通过 | 四个实现提交及单元测试 |
| L5 健壮性与诊断 | 通过 | 可定位事件、逐级路径、唯一 owner 和 MPI 归约测试 |
| L6 可选预处理 | 未纳入 | 未获单独批准；经典三种求解器保持无预处理基线 |
| L7 自动验收 | 通过 | 规格 6/6、串行 9/9、MPI 15/15 |
| 人工放行 | 待项目负责人检查 | 未进入阶段 M，未合并 `main` |

## 6. 已知边界与后续事项

1. 本阶段是面算法和多块/MPI 路径专项验收，不包含正式算例配置、初场、监测、残差文件、
   CGNS 流场输出或重启；这些仍按路线图属于阶段 M。
2. 因此 Sod 状态在本阶段用于 Riemann 局部激波健壮性测试，不等同于通过生产驱动运行完整
   二维 Riemann 算例。
3. Weiss--Smith 预处理必须同时实现稳态伪时间项、预处理特征系统、谱半径/CFL 和远场边界；
   在单独批准前保持不存在，物理时间 SSPRK3 不会误用不完整预处理。
4. 诊断事件当前保存在每次残差求值的本地对象中并提供 MPI 全局计数；阶段 M 的
   `MonitorManager` 负责持久化日志和历史输出。

## 7. Git/GitHub 状态

阶段实现按 L1--L5、L7 小目标分别提交，本地分支和候选标签完整。向
`https://github.com/wylx-2/wcns_v3.git` 推送时，GitHub 443 连接先被重置、随后连接超时；
因此远程分支和 `stage-l-candidate-v1` 标签尚不能确认同步，网络恢复后必须重试。阶段 L
未合并到 `main`。

后续状态更新：2026-09-02 网络恢复后，`stage-l-candidate-v1` 已成功补推到上述 GitHub
仓库；本节前述失败记录保留为候选建立时的历史事实。

## 8. 后续流程变更（2026-09-02）

项目负责人在本报告形成后明确取消 K/L 到发布候选之间的逐阶段人工等待，要求先连续完成
可发布程序，再通过具体生产算例进行集中验收。因此第 5 节“未进入阶段 M”的表述保留为本次
验收当时的历史状态，不再是当前前置卡口；它也不追认阶段 K 或 L 已人工通过。后续按
[`development-roadmap.md`](development-roadmap.md) 和
[`release-development-plan.md`](release-development-plan.md) 直接进入阶段 M，最终只在阶段 O
发布候选执行一次集中人工放行。
