# v1.1.0 阶段 P 自动验收报告

状态：**本地自动卡口通过；依据 2026-09-11 连续推进授权进入 Q。远程尚未同步。**

验收日期：2026-09-11（Asia/Shanghai）

阶段实现提交：`6a0076eab02c8508fab18a86f5a8385620990090`

基线/范围提交：`0032176`
正式版本基线：`v1.0.0`

## 1. 环境

- Windows 11 `10.0.26200`，Intel Family 6 Model 151，12 logical CPUs；
- CMake 3.28.0，MinGW-w64 `C:/mingw64/bin/c++.exe`，Release；
- MPI 构建使用本机 Intel MPI，串行构建 `WCNS_ENABLE_MPI=OFF`，MPI 构建为 `ON`；
- Python 3.14.5；仓库内 CGNS 4.4.0 固定压缩包。

## 2. 正式命令与结果

```powershell
cmake -S . -B build-rc-serial -DWCNS_ENABLE_CGNS=ON `
  -DWCNS_ENABLE_MPI=OFF -DWCNS_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-rc-serial --config Release --parallel 4
python tools\verify_algorithm_spec.py
ctest --test-dir build-rc-serial -C Release --output-on-failure
```

结果：配置/构建成功；规格检查 6/6；CTest **48/48**，总计 36.10 s。

```powershell
cmake -S . -B build-rc-mpi -DWCNS_ENABLE_CGNS=ON `
  -DWCNS_ENABLE_MPI=ON -DWCNS_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-rc-mpi --config Release --parallel 4
ctest --test-dir build-rc-mpi -C Release --output-on-failure
```

结果：配置/构建成功；CTest **84/84**，总计 95.93 s。矩阵实际包含 2-rank、4-rank 和既有
8-rank 拓扑测试。

```powershell
python tools\run_v110_performance_baseline.py `
  --run build-rc-serial\wcns_run.exe `
  --generator build-rc-serial\wcns_generate_release_cgns.exe `
  --work-dir build-p-performance --warmups 1 --repetitions 5
```

结果：三项均完成 1 次预热、5 次正式采样，且 CV 均小于 5%。原始数值见
`stage-p-performance-baseline.json`。

| 工作负载 | 中位 wall time | CV | 峰值进程树 RSS |
|---|---:|---:|---:|
| 100x100 二维等熵涡，20 步 | 4.5033814 s | 0.3513% | 50,573,312 B |
| Case07 四块 48x32 Mach 5 圆柱，20 步 | 0.5177693 s | 1.5873% | 17,244,160 B |
| 36x48x36 三维黏性，3 步 | 67.9583291 s | 0.3885% | 151,310,336 B |

以上性能值只适用于本机本构建；阶段 T 在加入分段计时后再冻结热路径指标。

## 3. Case07 与 NACA0012 范围证据

- `wcns.generate_cylinder_o_mesh` 和 `wcns.run.cylinder_o.dry_run.serial` 在串行、MPI 构建的
  CTest 清单中均通过；性能基线还实际生成并推进了 Case07 四块 48x32 网格 6 次。
- 基线包含提交 `9f0d9fa`、`3a50dc0`。
- `cases/manual/case06_2d_naca0012/` 始终是唯一未跟踪项，未被读取、修改或暂存；CI、测试与
  性能基线没有 NACA0012 条目。

## 4. CI 与远程状态

`.github/workflows/ci.yml` 已定义 Linux GCC/Clang 串行、GCC/OpenMPI 2/4-rank、Windows MSVC
串行和 Clang ASan/UBSan 作业，失败保留 JUnit/CTest 日志。本地没有 Linux/MSVC runner，
因此不能把本地结果冒充跨平台结果。

阶段分支推送请求被当前执行环境的外部写入审批拒绝，故本报告明确标记 **远程未同步，远程
CI 未执行**。按路线图 3.3 的离线条款，这不阻塞 P→Q 的本地连续开发；在获得明确推送授权后
必须补推阶段分支/候选标签，并连续核对两次远程矩阵。若任一远程作业失败，P 候选转为失败
候选并在 Q 分支继续集成前修复；不得移动现有标签。

## 5. 卡口结论

- V0/V1：范围、非目标、Case07/NACA 边界和 Q 数学输入已冻结；
- V2--V6：本机构建、规格、48/84 回归、Case07 和性能协议全部通过；
- V7：本地设计、脚本、JSON、报告和提交证据完整，远程状态已如实记录；
- V8：P→Q 的人工停点已由 2026-09-11 项目负责人连续推进指令预先豁免。

结论：允许创建本地 `v1.1.0-p-candidate.1`，非快进合入 `release/v1.1.0`，并立即开始 Q。
