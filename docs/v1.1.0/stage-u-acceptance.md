# v1.1.0 阶段 U 自动验收报告

状态：**2026-09-11 本地全部冻结自动卡口通过并创建候选 `v1.1.0-rc.1`。项目负责人随后
完成人工核验并接受当前计算结果，批准在无功能变更清理通过后合入 release/main、创建正式
标签并同步远程。Linux/Clang/OpenMPI/Sanitizer 外部 CI 不作为本机独立开发的完成卡口；
对外发布和许可证决定暂缓。**

机器可读结果见 [`stage-u-validation-summary.json`](stage-u-validation-summary.json)，发布说明见
[`../release-notes-1.1.0.md`](../release-notes-1.1.0.md)，冻结范围和阈值见
[`stage-u-design.md`](stage-u-design.md)。

## 1. 来源、环境和范围

- T 已由 `v1.1.0-t-approved.1` 记录人工批准，并由非快进提交 `9d276b2` 合入
  `release/v1.1.0`；U 分支为 `stage/v1.1.0-u`。
- 最终完整构建、CTest、规格和性能绑定实现/文档提交 `843b88777b66`；最终源码包来源
  `c7cb353de21d` 只比它多一个发布说明自引用修正，不改变代码、配置、测试或性能路径。
- 本机：Windows 11 10.0.26200、MinGW-w64 GCC 8.1.0、CMake 3.28.0、Intel MPI 2021.10，
  Release，12 logical CPU。
- `cases/manual/case06_2d_naca0012/` 不属于 v1.1.0。U 未读取、修改、暂存或打包该目录；
  Case07 明确进入验收。

## 2. 自动卡口结论

| 卡口 | 结果 | 主要证据 |
|---|---:|---|
| U-V0 范围/规格 | 通过 | U 设计冻结；Case07 输出保护公式已补入《算法补充》11.4.4 |
| U-V1 洁净构建 | 通过 | 新目录串行/MPI Release；编译输出无 warning |
| U-V2 全量回归 | 通过 | 串行 60/60（28.47 s），MPI 108/108（78.09 s），规格 6/6 |
| U-V3 数值/MPI | 通过 | 自由流、涡、Sod、四象限、黏性、源项、制造解、1/2/4/8 rank |
| U-V4 强激波 | 通过 | Case07 Mach 5，4 rank，`t=8`，12084 步；最终场独立重读 |
| U-V5 输出/重启 | 通过 | CGNS/Tecplot/history/statistics/boundary/load/manifest；100=40+60 |
| U-V6 失败路径 | 通过 | 16 类错误、1/2/4 rank、单命令 30 s 上限，无死锁/伪成功 |
| U-V7 性能 | 通过 | 热路径两项加速均大于 1.5x；Case07 不回退；4-rank 效率 84.17% |
| U-V8 打包/安装 | 通过 | 两包字节一致、包内/sidecar 哈希通过、无 Git/manual case、无 Git 构建安装运行 |
| U-V9 Git/文档 | 通过 | P--T 标签与非快进合并完整；相对链接 0 断链；`git diff --check` 通过 |

## 3. Case07 完整发布卡口与验收修复

最终 `cylinder_mach5_robust` 保持冻结的四块 48×32 O 网格、WENO-Z characteristic、Rusanov、
`CFL=0.02`、4 ranks 和 `t=8`。通过结果：

- 12084 步命中 `physical_time_reached`，运行约 134.29 s；17 行载荷事件、最终 checkpoint 和
  最终 CGNS 均存在，无 `.tmp`；manifest 为 `program_version=1.1.0`。
- 最大 level-1/2/3 面数为 111/65/156；最大 troubled cell 60，局部重算 18，整步 retry 3。
- history 中 436 个采样行包含被拒绝的非正“尝试候选态”，每行均有 troubled-cell 与重算或
  重试证据；这正是阶段 Q 的恢复语义，不是接受态正性失败。
- 最终接受候选最小值：`rho=1.8946665e-2`、`p=5.2859167e-4`、`T=7.3348284e-1`、
  `e=5.2391615e-2`；独立 CGNS 重读同样有限且为正。

首轮验收在 `t=0.5` 发现输出侧缺陷，未降低时长或放宽阈值：

1. Mach 5 Euler 配置不再请求载荷不需要的 `T_w,mu_w`；写出器只计算实际请求的热学量。
2. 激波附近无粘高阶压力/温度迹若非有限或不大于 floor，诊断值按《算法补充》11.4.4 退回
   最近内部真实单元正值；不回写状态，不进入通量、残差、时间步或 restart signature。
3. 边界测试改按列名而非固定列号验收；Case07 卡口正确区分“最差尝试候选”与最终接受态，
   并允许 manifest 合法的重复 `file=` 清单。

修复提交为 `a571ee4`、`eaa6857`、`2096a17`、`fe48922`。修复后定向边界/Sutherland 测试、
规格检查、完整 Case07、最终串行/MPI CTest 全部重跑通过。

## 4. 数值、并行、输运和失败路径

- 自由流：128×64 二维和 64×32×24 扭曲三维、两 profile、1/2/4/8 rank 以及原生四 zone
  全部通过；最大原始量误差 `3.39e-14`，rank 间最终场差为 0。
- 等熵涡：六个 profile/重构组合的 32²/64²/128²、`t=10` 全部通过；最小 64→128 阶
  3.799，最大 128² rho-L1 `5.368e-6`；32² 的 1/2/4/8-rank 最终场差为 0。
- 强间断：五个 Sod 组合最大 rho-L1 `1.484e-3`、位置误差 0.197 单元；200² 四象限对称
  L1 为 0；双马赫反射 smoke 在最终洁净 CTest 中通过。
- 黏性/输运：Couette 1/2/4 rank 壁切相对误差 `2.77e-11`；线性导热 Linf
  `4.44e-16`；constant/Sutherland 串行与 2-rank 生产路径均通过。Re=20/40/100 圆柱沿用
  R 阶段有来源提交和统计窗的 face-based 载荷证据，U 未用重复长算覆盖该权威结果。
- 源项/制造解：均匀源项最大误差 `1.78e-15`；二维制造解 1/2/4 rank、三维扭曲制造解
  1/2/4/8 rank 场差为 0；SSPRK3 时间阶为 3.012/3.144（卡口 2.8）。
- 16 类失败矩阵在 1/2/4 rank 全部按预期一致退出，每命令冻结超时 30 s。

## 5. 性能复验

协议为同机器、Release、固定逻辑核 `0,2,4,6,8,9,10,11`、High priority、一次预热、五次
样本；扩展组固定 5 步。结果：

| 工作负载 | 相对 P 加速 | CV |
|---|---:|---:|
| 二维等熵涡 100×100 | 1.589x | 0.202% |
| Case07 四 zone 48×32、20 步 | 1.164x | 0.799% |
| 三维黏性 36×48×36 | 7.942x | 0.194% |

allocation probe 为 `946214→6`，减少 99.9994%，第二次 residual 不回升。强扩展效率为
2-rank 96.29%、4-rank 84.17%、8-rank 57.01%；4-rank 高于 70% 硬卡口。弱扩展效率为
99.30%/87.21%/60.85%，8-rank 限制按计划保留。5 预热 + 5 样本的 detailed 对照最大正开销
仅 0.374%（卡口 2%）。

## 6. 最终源码包、安装与重启

最终包：`wcns-1.1.0-source-c7cb353de21d.tar.gz`

```text
SHA-256 562e9b2ef31aee88267f57ae94440f8ae45cc8b279f1f6ec171816b41457ba8a
size    2644573 bytes
files   271 tracked payload files
```

两个独立输出目录生成的 archive 字节完全相同，sidecar 和 `PACKAGE_CONTENTS.sha256` 均通过。
包不含 `.git`、构建目录、运行输出或任何 `cases/manual` 内容，顶层
`WCNS_SOURCE_REVISION=c7cb353de21d`。

在此前不存在且无 `.git` 的解包目录完成 Release+Intel MPI 构建、安装和运行。安装树文件数为
bin/include/lib/share = 11/8/3/64，含 README、LICENSE、第三方通知、算法补充、发布说明、
完整配置和可运行模板。安装二进制的 1/2-rank smoke 通过；100 步连续计算与 40+60 重启在
1/2 rank 的 history 最大差 `3.89e-16`、statistics 最大差 `1.11e-14`。

## 7. CI、许可证和已知限制

`.github/workflows/ci.yml` 静态审计确认以下作业仍使用当前 CMake/CTest/规格入口：Linux GCC、
Linux Clang、Linux OpenMPI、Windows MSVC 和 Clang ASan+UBSan。当前机器没有 Clang，本机
WSL 枚举被系统拒绝，且阶段分支未同步到 origin，因此这些外部结果是**待执行**，不是通过。

WCNS 自有代码尚未选择对外许可证；RC 只能私有/内部使用。包已包含 `LICENSE.md` 和
`THIRD_PARTY_NOTICES.md`，没有声称可以自由再分发。并行 CGNS、弱扩展/8-rank 强扩展、低 Mach、
湍流/转捩等限制继续按 [`../known-limitations.md`](../known-limitations.md) 公开。

## 8. 正式验收中被卡口拒绝的尝试

以下尝试均未记为通过，且修复后从冻结入口重跑：

- 打包初版对 Git Unicode 路径和 Windows archive 名称处理错误，修复后两份包重新生成并验证；
- 输出/重启工具拒绝了错误的 12/5 参数，因为协议冻结为 100/40+60；
- 自由流参考目录、涡 rank-only 最细网格阈值、二维/三维制造解最少步数等调用参数错误均由
  工具拒绝，改正后重跑；尝试把非线性 Sutherland 导热与线性温度解析解比较被识别为模型
  假设错误并安全终止，Sutherland 证据改用 S 阶段解析参数测试、生产 transient/CFL 和 CTest；
- Case07 的运行期壁面输出、manifest 重复 `file=`、缩写提交号、候选态语义四类问题均被卡口
  拒绝；最后一次完整运行同时满足终时、恢复、最终正性、独立重读和精确溯源；
- 性能入口拒绝了错误基线文件；allocation 探针旧路径导致的一组不完整样本没有复用，修正后
  整个正式性能矩阵从头重跑。

## 9. Git 审计和人工停止点

`release/v1.1.0` 第一父链包含 P/Q/R/S/T 的五个非快进合并：`651ceae`、`47ce10b`、
`652a102`、`775a48e`、`9d276b2`。P--T 候选标签、Q/T 批准标签均存在；U 来源包含 T 合并。
最终 `git diff --check` 和 62 个 tracked Markdown 文档相对链接审计通过。

远程写入授权此前未获批准，因此分支和标签保持本地未同步。工作树唯一例外应为已登记且未触碰
的 `cases/manual/case06_2d_naca0012/`。

人工核验建议逐项检查：Case07 完整结果及输出迹保护、Re=20/40/100 face-based 载荷、S 阶段
制造解/时间阶、T/U 性能与扩展、源码包哈希、外部 CI 待执行状态和许可证限制。

本报告提交时的停止条件已由 2026-09-11 人工批准解除。其后的源码清理和复验不改写本报告
记录的 RC 数值证据，独立记录在
[`release-cleanup-acceptance.md`](release-cleanup-acceptance.md)。正式标签仅表示内部版本
基线，不表示已经取得对外发布或再分发授权。
