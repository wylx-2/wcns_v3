# v1.1.0 阶段 T 自动验收报告

状态：**2026-09-11 自动卡口通过；候选 `v1.1.0-t-candidate.1` 完成后停止，不合入
`release/v1.1.0`、不建立 U 分支，等待 T→U 人工判断。**

## 1. 范围与设计结论

- [`stage-t-design.md`](stage-t-design.md) 和根目录《算法补充》11.6 冻结工作区键、版本语义、
  分配减少率、CV、串行加速比和强/弱扩展公式。
- 热路径缓存了 profile 轴向离散算子，以解析 Roe 左特征行替换固定结构的通用 5x5 求逆，
  并持久化 SSPRK 状态、无粘/黏性中间场、registry、halo plan、消息缓冲和 MPI request 槽。
- 生产残差直接引用缓存的 stencil 行；公开 stencil API 仍保留给 Q 阶段 troubled-cell 离散
  支持传播，没有把性能优化变成另一套数值算法。
- 通信保持同步 `exchange`。现有核没有经过验证的块内部/连接邻域拆分，因此没有发布伪
  `begin/finish` 重叠接口。
- 停止、错误、残差和 observer 归约位于不同控制依赖点，本阶段审计后不合并集合通信，避免
  改变失败和停止语义。
- 场输出与 checkpoint 按 source zone、quantity/component 流式写出；root 复用本地 gather
  缓冲并用 `MPI_IN_PLACE` 接收，浮点暂存只有一份全局 quantity payload。Tecplot 使用标准
  `DATAPACKING=BLOCK`，读取器兼容既有 POINT 文件。
- `cases/manual/case06_2d_naca0012/` 不属于本阶段；没有读取、修改或加入 Git。

## 2. 实现与 Git 审计

| 提交 | 内容 |
|---|---|
| `996d5e6` | 冻结 T 性能设计、热点和卡口 |
| `2b78344` | 冻结 S 候选热路径分配基线 |
| `91c3d46` | 缓存不变的轴向离散算子 |
| `4bb21e2` | 使用解析 Roe characteristic inverse |
| `942a6eb` | 复用 solver 场和 SSPRK 工作区 |
| `2893aa9` | 消除生产重构中的不变 stencil/参数开销 |
| `31380a7` | 持久化 halo 消息和 request 缓冲 |
| `4996207` | 流式 root 场输出和 checkpoint |
| `8d2a0e9`--`e97d718` | 建立性能矩阵并修复工作目录、复用、弱扩展网格和定向复测 |
| `cabb928` | 固化无采样默认计时、绑核、高优先级、5 步扩展及精确终点协议 |
| `995d6ff` | root 原地复用 gather payload，并验证非零 root |

每个主题均独立提交。最后的性能原始矩阵记录被测提交 `e97d718`；随后 `cabb928` 只修改测量
工具和文档，`995d6ff` 只修改输出启用时的 root gather。正式串行/扩展性能配置关闭输出，故
从 `e97d718` 到最终实现提交 `995d6ff` 的被测 solver 热路径没有变化。最终实现另行完成完整
Release 回归。

远程推送没有执行：此前外部网络写入授权未获批准。本地阶段分支、提交和候选标签完整，远程
状态明确为未同步。

## 3. 测量协议与失败复测记录

固定机器为 Windows 11、Intel64 Family 6 Model 151、12 logical CPU、MinGW Release 和
Intel MPI。每项正式结果保留一次预热和五次样本；扩展工作负载统一执行 5 步。混合核映射
固定为 `0,2,4,6,8,9,10,11`，MPI worker 使用 High priority，策略写入 JSON 和每条命令。

初次 2 步 strong 矩阵没有通过：r1 CV 为 5.03%，r8 CV 为 19.4%；单独增加预热、仅绑核或
去掉 RSS 采样后，r8 仍会因一个 3--6 s 的系统抢占样本超过 5%。这些失败样本没有删除或当作
通过。最终修复包括：

1. 默认模式完全关闭进程树 RSS 采样，符合冻结设计；
2. detailed 采样移到独立线程，墙钟在子进程 `wait` 返回时截取，不计采样线程停止和日志关闭；
3. 固定混合核映射和 worker 优先级，找不齐 rank 或设置失败时硬失败；
4. strong/weak 统一由 2 步增至 5 步，使系统抢占不再主导短样本；阈值、网格和统计公式未改；
5. manifest 复用必须匹配 detailed、预热/重复、步数、绑核和优先级协议。

最终 r8 五次为 16.2395、15.0857、15.1584、15.0284、14.6808 s，CV 3.46%。全部 40 个
strong/weak 正式日志均到达 step=5，没有提前终止样本。

## 4. 分配、串行吞吐和详细计时

分配减少率、加速比和计时开销分别使用
`1-A_new/A_baseline`、`T_baseline/T_new` 和 `T_detailed/T_default-1`。完整五样本数据见
[`stage-t-validation-summary.json`](stage-t-validation-summary.json)。

| 卡口 | 结果 | 阈值 |
|---|---:|---:|
| 基线每残差分配 | 946,214 次 / 108,367,664 B | 冻结值 |
| T 首次残差 | 6 次 / 208 B | -- |
| T 第二次残差 | 6 次 / 208 B | 不回升 |
| 分配减少率 | 99.99937% | ≥90% |
| 100x100 等熵涡 | 2.7862 s，CV 0.38%，1.616× | CV≤5%，加速≥1.5× |
| Case07 四块圆柱 | 0.4496 s，CV 1.85%，1.152× | CV≤5%，不回退 |
| 36x48x36 三维黏性 | 17.3659 s，CV 4.17%，3.913× | CV≤5%，加速≥1.5× |

详细计时开销对照使用默认和 detailed 各五次预热、五次正式样本：

| 算例 | default 中位数 | detailed 中位数 | 开销 | 两组最大 CV |
|---|---:|---:|---:|---:|
| 等熵涡 | 2.8306 s | 2.8297 s | -0.03% | 0.47% |
| Case07 | 0.4458 s | 0.4453 s | -0.10% | 1.09% |
| 三维黏性 | 8.5797 s | 8.6991 s | +1.39% | 0.41% |

三项均满足 2% 上限；负值视为测量噪声，不声称 detailed 能加速求解。

## 5. 强扩展与弱扩展

强扩展效率使用 `E_p=T_1/(p*T_p)`；弱扩展使用 `E_p^weak=T_1^weak/T_p^weak`。

| rank | strong 中位数 | CV | cells/rank | strong 效率 |
|---:|---:|---:|---:|---:|
| 1 | 65.3508 s | 3.59% | 165,888 | 100.0% |
| 2 | 30.1058 s | 0.30% | 82,944 | 108.5% |
| 4 | 19.2048 s | 1.10% | 41,472 | 85.1% |
| 8 | 15.0857 s | 3.46% | 20,736 | 54.1% |

4-rank 的 41,472 cells/rank 高于 15,552 下限，85.1% 高于 70% 下限。超线性的 2-rank
结果只按实记录，不外推为普适结论。

| rank | weak 中位数 | CV | cells/rank | weak 效率 |
|---:|---:|---:|---:|---:|
| 1 | 4.1781 s | 1.31% | 16,128 | 100.0% |
| 2 | 4.1917 s | 0.77% | 16,128 | 99.7% |
| 4 | 4.8437 s | 0.74% | 16,128 | 86.3% |
| 8 | 9.9730 s | 0.24% | 16,128 | 41.9% |

弱扩展按设计只保存效率，不设置掩盖本机内存带宽和混合核限制的事后硬阈值。8-rank 效率
下降是提交 T→U 人工判断的明确限制。

## 6. I/O、数值安全与完整回归

- root 浮点数据暂存由两份降为一份 source-zone/quantity payload；非零 root 和可变 rank
  payload 的串行、2、4、8-rank 测试通过。
- 场输出、CGNS/Tecplot 独立重读、checkpoint、同/异 rank restart 连续性专项测试：串行
  13/13、MPI 20/20；标准 BLOCK schema 和旧 POINT 兼容均通过。
- 版本、profile/extent、消息长度、NaN、错误注入、重启签名以及 1/2/4/8-rank 无死锁路径
  包含在单元、release smoke 和完整 MPI 矩阵中。
- 最终提交上的 Release 串行 CTest 60/60（25.65 s），MPI CTest 108/108（73.89 s）；
  `python tools/verify_algorithm_spec.py` 6/6；`git diff --check` 通过。

## 7. 正式自动卡口结论

| 卡口 | 结果 |
|---|---|
| T-V1 数值等价 | 通过；P--S 回归、输出/restart 及 MPI 冻结容差未放宽 |
| T-V2 分配 | 通过；减少 99.99937%，第二次不回升 |
| T-V3 串行吞吐 | 通过；两项代表计算核 ≥1.5×，Case07 无回退，全部 CV≤5% |
| T-V4 计时开销 | 通过；最坏 +1.39% |
| T-V5 强扩展 | 通过；4-rank 85.1%，41,472 cells/rank |
| T-V6 弱扩展 | 完整报告；8-rank 限制明确保留 |
| T-V7 I/O | 通过；单 payload、独立重读和 restart 均通过 |
| T-V8 安全 | 通过；版本/长度/NaN/错误注入无死锁且一致失败 |
| T-V9 回归 | 通过；60/60、108/108、规格 6/6、diff clean |

正式命令包括：

```text
python tools/run_stage_t_performance.py ... --repetitions 5 --warmups 1 \
  --scaling-repetitions 5 --scaling-warmups 1 --scaling-steps 5 \
  --mpi-pin-processor-list 0,2,4,6,8,9,10,11 --mpi-high-priority
python tools/run_stage_t_performance.py ... --repetitions 5 --warmups 5 --skip-scaling
python tools/run_stage_t_performance.py ... --repetitions 5 --warmups 5 --detailed --skip-scaling
cmake --build build --config Release -j 2
cmake --build build-mpi --config Release -j 2
ctest --test-dir build --output-on-failure -j 1
ctest --test-dir build-mpi --output-on-failure -j 1
python tools/verify_algorithm_spec.py
git diff --check
```

## 8. T→U 人工判断与停止点

T 的全部自动卡口通过，但尚未获得人工批准。请项目负责人重点判断：解析 Roe 逆矩阵和缓存
复杂度是否可维护；不实施通信重叠/归约合并的取舍；8-rank strong 54.1%、weak 41.9% 是否
可接受；root 分块 I/O 的接口复杂度；以及优化收益是否足以进入发布级 U。

候选创建后，本轮严格停在 T→U：不合入 `release/v1.1.0`、不建立 `stage/v1.1.0-u`、不创建
任何 U/RC/正式版本标签。人工批准后才能继续阶段 U。
