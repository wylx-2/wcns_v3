# v1.1.0 阶段 P 设计冻结：规格、CI 与可重复基线

状态：2026-09-11 已按项目负责人“P/Q 连续推进”指令冻结。P 自动验收通过后直接进入 Q，
不另设 P→Q 人工停点。

## 1. 基线与范围

- 正式发布基线：`v1.0.0`；v1.1.0 计划基线：`0032176`。
- v1.0 后纳入基线的 Case07 提交：`9f0d9fa`（四块圆柱 O 网格生成）和 `3a50dc0`
  （低速及 Mach 5 圆柱验证）。P 必须运行 `wcns.generate_cylinder_o_mesh` 与
  `wcns.run.cylinder_o.dry_run.serial`，后续 Q 使用 48x32 Mach 5 配置。
- `cases/manual/case06_2d_naca0012/` 是既有未跟踪用户目录，P/Q 暂不包含，不读取、不修改、
  不暂存。验收前后 `git status --short` 均允许它作为唯一已登记例外。
- P 不改变方程、数值离散、配置 schema、输出 schema 或重启签名。

## 2. CI 矩阵

权威入口为 `.github/workflows/ci.yml`：

| 作业 | 编译器/模式 | 必须覆盖 |
|---|---|---|
| Linux serial | GCC、Clang，Release | 规格核验和完整串行 CTest |
| Linux MPI | GCC + OpenMPI，Release | 完整 MPI 构建中的串行、2-rank、4-rank CTest |
| Windows serial | MSVC x64，Release | 规格核验和完整串行 CTest |
| Sanitizers | Clang Debug + ASan/UBSan | 完整串行 CTest，首个错误即失败 |

CI 使用仓库内带哈希校验的 CGNS 4.4.0 包，不访问运行期外部网格；只保留失败的 JUnit、CTest
临时日志，不上传成功场文件。普通 CI 不设置墙钟性能阈值。

## 3. 功能与数值基线

本机正式基线使用两个既有干净 Release 构建：串行 48 项、MPI 84 项。矩阵覆盖自由流、等熵
涡、Sod、四象限、双马赫反射 smoke、Couette、导热、制造源、周期多块、输出/检查点、异 rank
重启和失败路径。测试数量或命令若变化，验收报告必须列出差异，不能只记录“全部通过”。

Case07 属于 P 的明确卡口：测试必须实际生成 32x20、四原生 zone 的 O 网格，再由正式
`wcns_run --dry-run` 读取、建立连接并生成 SCMM6 度量。NACA0012 不作为替代或补充。

## 4. 性能基线协议

P 冻结三种工作负载，正式阈值仅在固定机器使用：

1. 100x100 周期二维等熵涡，`scmm6_wcns + weno_z + hllc`，无中间场输出；
2. Case07 四块 48x32 圆柱 O 网格，短推进且关闭中间场/统计/checkpoint 输出；
3. 36x48x36 三维黏性短推进，关闭中间输出。

每项先预热一次，再执行五次；记录每次 wall time、进程树峰值 RSS、中位数、标准差和变异系数
`CV=sigma/mean`。`CV>5%` 时报告环境不稳定，不冻结退化阈值。P 只冻结可由现有外部采样得到的
端到端时间/RSS；初始化、度量、无粘、粘性、halo、归约和 I/O 分段计时在阶段 T 加入生产
instrumentation 后再建立，不用估算值冒充测量。

## 5. 自动卡口与 Git

1. `python tools/verify_algorithm_spec.py` 全通过；
2. 串行与 MPI Release 配置/构建成功，完整 CTest 分别 48/48、84/84；
3. Case07 两项测试通过，测试清单中不存在 NACA0012；
4. CI YAML 可解析，四类作业均有失败证据保留；远程连续两次运行结果在网络同步后补录；
5. `git diff --check` 通过，状态只含已登记的 NACA0012 用户目录；
6. 自动验收报告提交后创建不可移动的 `v1.1.0-p-candidate.1`，再以 `--no-ff` 合入
   `release/v1.1.0`，建立 `stage/v1.1.0-q`。

若 GitHub 不可访问，本地候选和合并可按路线图继续，但验收报告必须写明“远程未同步”。
