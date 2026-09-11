# case05 Linux 服务器编译与长算操作手册

## 1. 文档目的与计算边界

本手册对应最小源码包 **wcns-case05-linux-<revision>.tar.gz**。包中只保留：

- 编译 WCNS MPI/CGNS 可执行程序所需的 CMake、C++ 源码和头文件；
- CGNS 4.4.0 原始源码归档及必须保留的许可证通知；
- case05 的五步迁移卡口配置、长算配置、网格生成/结果校验脚本；
- 本操作手册和源文件 SHA-256 清单。

包中没有其他算例、完整回归测试、历史计算结果、Windows 构建目录或通用开发文档。因此，
服务器上的验收方式是“从源码编译 + 4-rank 五步卡口 + 目标 rank dry-run + 正式长算”，
不是运行完整 CTest 矩阵。完整的 87 项 MPI Release 回归已经在打包前通过。

本例网格为 36×48×36，计算域为

\[
[0,2\pi]\times[-1,1]\times[0,\pi],
\]

x/z 周期，y 两面为等温无滑移壁，流向由定常体积力驱动。参考速度是初始体积平均速度
\(U_{b,0}\)，初始体积平均 Mach 数为 0.1。它用于验证服务器迁移、并行运行、重启和统计
链路；该网格不足以声明得到定量 DNS 结果。

### 本次低马赫参考量变更

当前配置按已审核约定取

\[
U_{ref}=U_{b,0}=1,\quad L_{ref}=h=1,\quad
\rho_{ref}=1,\quad T_{ref}=71.42857142857143,\quad
\mu_{ref}=1/(180U_b^+)=0.0003588401476178181.
\]

其中 \(U_b^+=15.4819787932\)，所以

\[
Re_b^{(h)}=U_b^+Re_\tau=2786.756182769849,
\qquad Ma_{b,0}=0.1.
\]

程序启动日志必须显示 `Re=2786.7561827698491 Ma=0.1`。旧包采用
\(U_{ref}=u_\tau\)，其日志虽同样显示 `Ma=0.1`，实际体积平均 Mach 约为 1.548；
旧包、旧配置和旧 checkpoint 都不能作为本次低马赫算例的续算起点。

## 2. 包内目录及各文件用途

解压后的目录结构如下：

~~~text
wcns-case05-linux-<revision>/
├── CMakeLists.txt
├── WCNS_SOURCE_REVISION
├── PACKAGE_CONTENTS.sha256
├── LICENSE.md
├── THIRD_PARTY_NOTICES.md
├── include/wcns/...
├── src/...
├── tools/
│   ├── generate_release_cgns.cpp
│   ├── validate_release_case.cpp
│   └── compare_metric_profiles.cpp
├── third_party/cgns/
│   ├── CGNS-4.4.0.zip
│   └── README.md
└── cases/manual/case05_3d_turbulent_channel/
    ├── LINUX_SERVER_GUIDE.md
    ├── channel_retau180_feasibility.wcns
    ├── channel_retau180_longrun.wcns
    ├── run_case05.py
    └── validate_case05.py
~~~

重要文件：

- **WCNS_SOURCE_REVISION**：打包时的 12 位 Git 提交号。源码包没有 .git 目录，CMake 会读取
  此文件并把版本写进每次运行的 manifest。
- **PACKAGE_CONTENTS.sha256**：包内所有有效载荷文件的校验值，不包含该清单自身。
- **channel_retau180_feasibility.wcns**：固定 4-rank、5 步迁移卡口，不能当作长算配置。
- **channel_retau180_longrun.wcns**：低马赫长算第 01 段模板，默认 \(t_{end}=500\)，最多
  6,000,000 步，程序墙钟 23 h。按本地五步卡口的初始步长粗估约需 418 万步；实际步长
  会随流场变化，600 万只是硬保护余量，正式提交前必须用服务器实测复核。
- **run_case05.py**：生成网格、执行 dry-run、推进五步并校验短算结果。
- **validate_case05.py**：检查统计列、壁摩擦、实测 \(Re_\tau\) 和两个 y-z 截面流量。

## 3. 第一步：在本地确认并发送压缩包

发布目录同时提供 tar.gz 和 tar.gz.sha256。Windows PowerShell 中先检查：

~~~powershell
Get-ChildItem .\dist\wcns-case05-linux-*.tar.gz*
Get-FileHash .\dist\wcns-case05-linux-<revision>.tar.gz -Algorithm SHA256
Get-Content .\dist\wcns-case05-linux-<revision>.tar.gz.sha256
~~~

两处 SHA-256 必须相同。然后用 scp 发送；将用户名、主机和目录替换成实际值：

~~~powershell
scp .\dist\wcns-case05-linux-<revision>.tar.gz user@login.example:/home/user/packages/
scp .\dist\wcns-case05-linux-<revision>.tar.gz.sha256 user@login.example:/home/user/packages/
~~~

如果单位使用 WinSCP、SFTP 网关或对象存储，也必须同时传输校验文件。不要传输
build-rc-mpi、.exe、CMakeCache.txt 或 Windows 生成的结果目录。

### 3.1 已有旧 Case05 服务器目录的升级方法

旧目录可能已经包含高体积平均 Mach 配置、构建缓存和计算结果。按以下顺序处理：

1. 如果旧作业仍在运行，先用调度器正常取消，使程序在完整时间步边界写 checkpoint；
   记录作业号、停止原因和旧 revision。不要使用 `kill -9`。
2. 将旧目录改成只读归档用途，保留其日志、manifest、statistics、field 和 checkpoint；
   不要在其中直接覆盖配置，也不要删除旧结果来腾出同名目录。
3. 把新的 `wcns-case05-linux-<new-revision>.tar.gz` 和校验文件上传到原 packages 目录，
   校验后解压为新的并列目录。新包名中的 revision 必须不同于旧包。
4. 在新目录中使用全新的 `build-linux-mpi`。CMake cache、目标文件和可执行程序都不能从
   旧目录复制，因为参考量已进入 restart signature，源代码版本和 CGNS/MPI ABI 也可能变化。
5. 在新目录重新执行 4-rank 五步卡口。新日志必须同时满足：

~~~text
case=case05-channel-retau180-mab0p1-feasibility
derived Re=2786.7561827698491 Ma=0.1
WCNS run stopped: reason=maximum_steps step=5
~~~

6. 新结果必须写入 `results/lowmach-*`。任何 `results/feasibility-r4`、
   `results/longrun-segment*` 或旧 case 名均属于旧基准，不能与新结果拼接。
7. **不得在新配置中设置旧 checkpoint 的 `restart.path`。**改变 `reference.viscosity`、
   `initial.bulk_velocity` 和 `source.body.ax` 既改变无量纲状态定义，也改变物理问题；
   本次必须从新的 `turbulent_channel` 初场冷启动。

如果服务器上保存的是完整 Git 克隆而不是最小包，推荐保留原 checkout，并建立新 worktree：

~~~bash
cd /path/to/existing/wcns_v3
git status --short
git fetch origin
git worktree add ../wcns-case05-lowmach <new-commit-or-tag>
cd ../wcns-case05-lowmach
~~~

`git status --short` 非空时不得强制 reset；先归档或人工处理服务器本地改动。若不熟悉
worktree，直接重新 clone 到新目录更安全。无论哪种路径，都必须从零配置构建目录并执行
相同的卡口，不能只复制两份 `.wcns` 文件后继续使用旧二进制。

## 4. 第二步：登录服务器并检查硬件、空间和软件

登录后不要立即在登录节点进行长算。先查看系统：

~~~bash
ssh user@login.example
uname -a
lscpu
df -h .
quota -s 2>/dev/null || true
ulimit -s
~~~

最低软件要求：

- CMake 3.20 或更高；
- 支持 C++20 的 GCC/G++，建议 GCC 10 或更高；
- 同一套 MPI 的编译器封装和启动器，例如 OpenMPI 的 mpicxx/mpiexec，或 MPICH 对应工具；
- Python 3.8 或更高；
- make 或 Ninja；
- tar 和 sha256sum。

先检查现有环境：

~~~bash
cmake --version
g++ --version
mpicxx --version
mpiexec --version
python3 --version
which cmake g++ mpicxx mpiexec python3
~~~

集群通常使用 environment modules。下面只是示例，实际模块名由管理员决定：

~~~bash
module purge
module avail gcc
module avail openmpi
module load gcc/12
module load openmpi/4
~~~

如果是独立 Ubuntu 主机且具有管理员权限，可安装：

~~~bash
sudo apt update
sudo apt install -y build-essential cmake python3 openmpi-bin libopenmpi-dev
~~~

Rocky/Alma/RHEL 主机可由管理员安装：

~~~bash
sudo dnf install -y gcc gcc-c++ cmake python3 openmpi openmpi-devel
~~~

不要混用不同 MPI：例如不能用 MPICH 的 mpicxx 编译后再用 OpenMPI 的 mpiexec 启动。
再次执行 **which mpicxx** 和 **which mpiexec**，确认二者属于同一模块前缀。

## 5. 第三步：校验压缩包并解压

进入接收目录：

~~~bash
cd /home/user/packages
sha256sum -c wcns-case05-linux-<revision>.tar.gz.sha256
~~~

必须看到 **OK**。失败时不要继续编译，应重新传输文件。

解压到独立目录：

~~~bash
tar -xzf wcns-case05-linux-<revision>.tar.gz
cd wcns-case05-linux-<revision>
pwd
~~~

校验包内每个文件：

~~~bash
sha256sum -c PACKAGE_CONTENTS.sha256
cat WCNS_SOURCE_REVISION
~~~

所有行都必须显示 **OK**。不要直接在压缩包所在目录混放不同 revision 的源码。

## 6. 第四步：从零配置并编译 MPI Release

在包根目录执行：

~~~bash
cmake -S . -B build-linux-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DWCNS_ENABLE_MPI=ON \
  -DWCNS_ENABLE_CGNS=ON \
  -DWCNS_BUILD_TESTS=OFF
cmake --build build-linux-mpi --parallel 8
~~~

参数含义：

- **-S .**：源码根目录是当前解压目录；
- **-B build-linux-mpi**：所有构建产物写入单独目录，不污染源码；
- **CMAKE_BUILD_TYPE=Release**：开启编译器优化，长算不得使用 Debug；
- **WCNS_ENABLE_MPI=ON**：要求 CMake 找到 MPI C++；
- **WCNS_ENABLE_CGNS=ON**：从包内 CGNS-4.4.0.zip 构建静态 ADF 后端，不访问网络；
- **WCNS_BUILD_TESTS=OFF**：最小包没有完整 tests 目录；
- **--parallel 8**：只控制编译并发，可按登录节点政策调小，不是 CFD 的 MPI rank 数。

成功后检查三个运行相关程序：

~~~bash
ls -lh build-linux-mpi/wcns_run \
       build-linux-mpi/wcns_generate_release_cgns \
       build-linux-mpi/wcns_validate_release_case
ldd build-linux-mpi/wcns_run | grep -Ei 'mpi|not found'
~~~

**not found** 表示动态库环境不完整，不能提交作业。若管理员要求在计算节点加载模块，
作业脚本里必须重复相同的 module purge/load。

更换编译器或 MPI 后不要复用旧 CMake cache，应改用新的构建目录，例如：

~~~bash
cmake -S . -B build-linux-mpi-gcc12 ...
~~~

## 7. 第五步：运行固定 4-rank 五步卡口

仍在包根目录执行：

~~~bash
python3 cases/manual/case05_3d_turbulent_channel/run_case05.py \
  --clean \
  --ranks 4 \
  --run build-linux-mpi/wcns_run \
  --generator build-linux-mpi/wcns_generate_release_cgns \
  --validator build-linux-mpi/wcns_validate_release_case \
  --mpi-exec mpiexec
~~~

脚本依次完成：

1. 生成 36×48×36、四个原生 zone 的周期 CGNS 网格；
2. 用 4 rank 执行配置、连接、分区、度量、初场 dry-run；
3. 真正推进 5 步；
4. 检查终场所有 \(\rho,p,T\) 有限且为正；
5. 检查统计文件包含两个 y-z 截面和两壁摩擦数据。

脚本中的 **--clean** 只清除本 case 的 grids/logs/results/validation。它会删除上一次短测证据，
所以第二次执行前应先归档需要保留的结果。

成功末尾应出现：

~~~text
derived Re=2786.7561827698491 Ma=0.1
WCNS run stopped: reason=maximum_steps step=5
case05 short feasibility validation completed; no turbulence acceptance claimed
~~~

程序把 **maximum_steps** 定义为安全上限，直接运行 wcns_run 时退出码是 2；
run_case05.py 已把短卡口的 2 视为预期。退出码 3 才是数值失败，退出码 1 是配置/CGNS/I/O
错误。脚本自身成功时返回 0。

查看关键证据：

~~~bash
CASE=cases/manual/case05_3d_turbulent_channel
tail -n 20 $CASE/logs/run-feasibility-r4.log
cat $CASE/validation/final-field-finite.txt
cat $CASE/validation/statistics-check.json
grep -E '^(git_commit|mpi_ranks|step|time|stop_reason)=' \
  $CASE/results/lowmach-feasibility-r4/*.manifest.r4.txt
~~~

manifest 中的 **git_commit** 必须等于 WCNS_SOURCE_REVISION。短测仅验证工程可运行性，
不验证湍流达到统计稳定。

## 8. 第六步：逐项理解长算配置

长算文件是：

~~~text
cases/manual/case05_3d_turbulent_channel/channel_retau180_longrun.wcns
~~~

所有相对路径都以该配置文件所在目录为基准。

### 8.1 身份和网格

| 键 | 当前值 | 含义及允许修改 |
|---|---:|---|
| schema_version | 1 | 配置语法版本，不可改 |
| case.name | ...segment01 | 输出文件前缀；每个续算段必须使用新名称 |
| mesh.path | grids/...cgns | 输入网格；短卡口会生成。改变网格后旧 checkpoint 不能重启 |

### 8.2 数值算法

| 键 | 当前值 | 含义 |
|---|---|---|
| algorithm.profile | phenglei_wcns | 度量算法路径；不可和 scmm6_wcns 混用 |
| algorithm.reconstruction | mdcd_hybrid | 六点界面重构 |
| algorithm.reconstruction_variables | primitive | 对原始变量重构 |
| algorithm.riemann | hllc | 面通量 Riemann 求解器 |
| algorithm.mdcd.disp | 0.0463783 | MDCD 色散参数 |
| algorithm.mdcd.diss | 0.01 | MDCD 耗散参数 |

算法、MDCD 参数或 Riemann 求解器属于数值重启签名。续算时改变它们会被拒绝；若要比较算法，
应从初场建立新 case、新输出目录并重新验收。

### 8.3 气体与参考尺度

| 键 | 当前值 | 作用 |
|---|---:|---|
| gas.gamma | 1.4 | 比热比 |
| gas.specific_gas_constant | 1.0 | 当前无量纲气体常数 |
| reference.velocity | 1 | \(U_{ref}=U_{b,0}\) |
| reference.density | 1 | \(\rho_{ref}\) |
| reference.temperature | 71.428571... | 使参考声速为 \(10U_{b,0}\)，即 \(Ma_{b,0}=0.1\) |
| reference.length | 1 | \(L_{ref}=h\) |
| reference.viscosity | 0.0003588401476 | 使 \(Re_b^{(h)}=2786.75618\) 与 \(Re_\tau=180\) 同时成立 |

这些量决定 Re、Ma、压力尺度、时间尺度和黏性，长算中途绝对不能修改。

### 8.4 MPI 自动分区

| 键 | 当前值 | 含义 |
|---|---|---|
| partition.mode | auto_split | zone 少于 rank 时递归切分结构块 |
| allow_idle_ranks | false | 不允许某些 rank 无计算块 |
| max_load_ratio | 1.2 | 最大/平均单元负载限制 |
| min_cells_per_active_direction | 8 | 每个活动方向的最小叶块宽度 |

当前四个原生 zone 的理论最大可行叶块数为 96。不要直接使用 96；先测试 4、8、12、16、
24 或 36 rank 的 dry-run，并结合服务器核数和通信效率选择。36×48×36 只有 62,208 个单元，
rank 过多通常会因通信成本变慢。

### 8.5 初始场

**initial.type=turbulent_channel** 使用复合壁律速度型和低波数三维扰动。

- y0/y1、x0/z0 和 period_x/period_z 必须和网格范围一致；
- re_tau=180 是目标摩擦 Reynolds 数；
- bulk_velocity=1，因为速度尺度是初始体积平均速度；
- bulk_velocity_plus=15.4819787932，定义 \(U_{b,0}/u_{\tau,0}\)；
- perturbation_amplitude=0.05 控制初始扰动；
- rho=1、temperature=1 是无量纲初始热力学状态。

使用 restart.path 后，守恒场来自 checkpoint，而不是重新采用 initial 初场。

### 8.6 边界和体积力

bottom/top 都是无滑移等温壁，壁速为零、壁温为 1；x/z 周期关系存储在 CGNS 网格中。
source.models=body_force 且 source.body.ax=0.004172026549973031，代表

\[
a_x^*=a_xh/U_{b,0}^2=1/(U_b^+)^2.
\]

体积力同时进入 x 动量和总能量。当前没有恒流量闭环控制器；流量会随瞬时流场变化。

### 8.7 推进和停止

| 键 | 长算值 | 说明 |
|---|---:|---|
| run.mode | unsteady | 非定常，残差只监测、不触发收敛 |
| run.viscous | true | 开启层流黏性通量 |
| run.cfl | 0.15 | 已过短测的保守值；提高前必须做独立稳定性试验 |
| run.max_steps | 6,000,000 | 硬步数上限；本地初始步长粗估约需 418 万步，服务器必须复核 |
| run.t_end | 500 | 目标物理时间，单位 \(h/U_{b,0}\)，约 79.6 个流向域穿越时间 |
| run.max_wall_time | 82,800 s | 23 h 后安全写 checkpoint，适配示例 24 h 作业 |

若调度墙钟不是 24 h，应令程序 max_wall_time 比调度器上限至少提前 10--60 min。例如：

- 调度 12 h：可设 39,600 s（11 h）；
- 调度 48 h：可设 169,200 s（47 h）。

正常达到 t_end 返回 0；达到 max_steps、程序墙钟或收到 SIGINT/SIGTERM 时在完整时间步边界
安全写 checkpoint 并返回 2；数值失败返回 3。

### 8.8 输出

- **field**：每 10 个物理时间单位写一次 CGNS 全场，并写初场/终场；
- **history**：每 100 步写残差、dt、回退计数和停止状态；
- **statistics**：每 0.1 时间单位写总量、两个 y-z 截面和壁摩擦；
- **checkpoint**：每 5 时间单位写可重启 CGNS，并维护 latest 别名；
- **output.allow_existing=false**：拒绝覆盖已有结果；
- **output.dimensional=false**：保存无量纲量。

两个 y-z 截面目标是 x=0 和 x=pi，自动列为：

~~~text
yz_mean_u_plane0       yz_mass_flow_x_plane0
yz_mean_u_plane1       yz_mass_flow_x_plane1
~~~

壁面统计列为：

~~~text
channel_wall_shear_lower  channel_wall_shear_upper
channel_wall_shear_mean   channel_friction_velocity
channel_re_tau
~~~

这些都是瞬时空间统计，不是累积时间平均、RMS 或 Reynolds 应力。

## 9. 第七步：用目标 rank 做长算前 dry-run

进入 case 目录：

~~~bash
cd cases/manual/case05_3d_turbulent_channel
ROOT=$(cd ../../.. && pwd)
RANKS=36
mpiexec -n $RANKS $ROOT/build-linux-mpi/wcns_run \
  --config channel_retau180_longrun.wcns --dry-run \
  > longrun-dry-run-r$RANKS.log 2>&1
echo $?
~~~

退出码必须为 0。检查：

~~~bash
grep -E 'derived Re=|partition_plan|WCNS dry-run' longrun-dry-run-r$RANKS.log
~~~

应看到 Re=2786.7561827698491、Ma=0.1、完整 partition_plan 和 dry-run completed。若提示块过窄、空闲 rank、
负 Jacobian、周期连接或 MPI digest 不一致，不得开始长算。

## 10. 第八步：在已分配资源中启动正式长算

不要在共享登录节点直接运行。先取得交互作业或编写调度脚本。通用 MPI 命令为：

~~~bash
cd cases/manual/case05_3d_turbulent_channel
ROOT=$(cd ../../.. && pwd)
RANKS=36
mkdir -p job-logs
set +e
mpiexec -n $RANKS $ROOT/build-linux-mpi/wcns_run \
  --config channel_retau180_longrun.wcns \
  > job-logs/segment01-r$RANKS.log 2>&1
RC=$?
set -e
echo "wcns exit code=$RC"
~~~

RC=0 表示达到 t_end；RC=2 表示安全检查点停止，需要查看 manifest 后续算；RC=1 或 3
必须先诊断，不能自动循环重启。

### Slurm 示例

将下面内容保存为包外的 **case05-segment01.sbatch**，并按集群修改 account、partition、
module 和 MPI 类型：

~~~bash
#!/bin/bash
#SBATCH --job-name=wcns-c05-s01
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail
module purge
module load gcc/12
module load openmpi/4

ROOT=/home/user/packages/wcns-case05-linux-<revision>
CASE=$ROOT/cases/manual/case05_3d_turbulent_channel
cd $CASE
mkdir -p job-logs

set +e
srun -n $SLURM_NTASKS $ROOT/build-linux-mpi/wcns_run \
  --config channel_retau180_longrun.wcns \
  > job-logs/segment01-slurm-$SLURM_JOB_ID.log 2>&1
RC=$?
set -e

if [ $RC -eq 0 ]; then
  echo "WCNS reached t_end"
elif [ $RC -eq 2 ]; then
  echo "WCNS stopped safely; inspect manifest and restart checkpoint"
else
  echo "WCNS failed with exit code $RC"
  exit $RC
fi
~~~

提交和查看：

~~~bash
sbatch case05-segment01.sbatch
squeue -u $USER
sacct -j <jobid> --format=JobID,State,Elapsed,ExitCode,MaxRSS
~~~

某些集群要求 srun --mpi=pmix、mpiexec 或 mpirun。只能按管理员提供的 MPI 示例修改，
不能猜测 MPI 插件。

### PBS 示例

~~~bash
#!/bin/bash
#PBS -N wcns-c05-s01
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -l walltime=24:00:00
#PBS -j oe

set -euo pipefail
cd $PBS_O_WORKDIR
module purge
module load gcc/12
module load openmpi/4

ROOT=/home/user/packages/wcns-case05-linux-<revision>
CASE=$ROOT/cases/manual/case05_3d_turbulent_channel
cd $CASE
mkdir -p job-logs

set +e
mpiexec -n 36 $ROOT/build-linux-mpi/wcns_run \
  --config channel_retau180_longrun.wcns \
  > job-logs/segment01-pbs-$PBS_JOBID.log 2>&1
RC=$?
set -e
test $RC -eq 0 -o $RC -eq 2
~~~

## 11. 第九步：运行中监控

查看 stdout、history 和 statistics：

~~~bash
tail -f job-logs/segment01-*.log
tail -n 20 results/lowmach-longrun-segment01/*.history.r*.txt
tail -n 20 results/lowmach-longrun-segment01/*.statistics.r*.txt
du -sh results/lowmach-longrun-segment01
df -h .
~~~

重点检查：

1. time 和 step 持续增加，dt 为正且没有异常快速缩小；
2. reconstruction_fallbacks 和 riemann_fallbacks 不持续爆发；
3. 总质量无非物理漂移；
4. 两个 y-z 截面流量数量级一致；
5. 上下壁剪切长期平均应逐渐接近；
6. \(\rho,p,T\) 不出现非有限或非正；
7. 磁盘空间足够容纳下一次 field/checkpoint。

需要提前停止时，优先让调度器发送 SIGTERM，或对 mpiexec 作业执行正常 cancel。程序在下一完整
时间步边界写安全 checkpoint。不要直接 kill -9；SIGKILL 无法保证 checkpoint 完整。

## 12. 第十步：检查停止原因和结果文件

运行结束后：

~~~bash
cd cases/manual/case05_3d_turbulent_channel
MANIFEST=$(ls -1t results/lowmach-longrun-segment01/*.manifest.r*.txt | head -n 1)
grep -E '^(git_commit|mpi_ranks|step|time|time_step|wall_time|stop_reason)=' $MANIFEST
grep -Ei 'nan|inf|numerical_failure|error' job-logs/segment01-*.log || true
ls -lh results/lowmach-longrun-segment01
~~~

输出文件：

- **field.step*.cgns**：可视化流场，cell-center 数据；
- **history.rN.txt**：步长、残差、五分量范数、回退计数和停止状态；
- **statistics.rN.txt**：总量、截面平均/流量、壁面摩擦和实测 \(Re_\tau\)；
- **checkpoint.step*.cgns**：不可变的分段检查点；
- **checkpoint.latest.cgns**：最近安全检查点别名；
- **manifest.rN.txt**：版本、配置、网格/重启签名、rank 数、终止原因和输出清单。

可独立检查最终场：

~~~bash
FINAL=$(ls -1t results/lowmach-longrun-segment01/*.field.*.cgns | head -n 1)
$ROOT/build-linux-mpi/wcns_validate_release_case finite $FINAL
$ROOT/build-linux-mpi/wcns_validate_release_case nonzero $FINAL Mach 0
~~~

## 13. 第十一步：从安全检查点续算

假设第 01 段因 wall_time_checkpoint 返回 2：

1. 确认 latest 文件存在且 manifest 的 stop_reason 是 wall_time_checkpoint 或
   user_signal_checkpoint；
2. 复制配置：

~~~bash
cp channel_retau180_longrun.wcns channel_retau180_longrun_segment02.wcns
~~~

3. 编辑 segment02 文件，只修改下列运行/输出项：

~~~text
case.name = case05-channel-retau180-mab0p1-longrun-segment02
output.directory = results/lowmach-longrun-segment02
restart.path = results/lowmach-longrun-segment01/case05-channel-retau180-mab0p1-longrun-segment01.checkpoint.latest.cgns
~~~

4. t_end 是从初始时刻计数的最终绝对时间，必须大于 checkpoint 内的 time；
5. 可以改变 rank 数、输出频率和 max_wall_time；
6. 不得改变网格、profile、重构、Riemann、气体、参考量、边界、源项或 viscous 开关；
7. 先 dry-run：

~~~bash
mpiexec -n 36 $ROOT/build-linux-mpi/wcns_run \
  --config channel_retau180_longrun_segment02.wcns --dry-run
~~~

8. dry-run 通过后再提交 segment02 作业。

每段使用新 case.name 和新 output.directory，可以防止覆盖前段证据。不同合法 rank 数之间允许
重启，因为 checkpoint 按原始 CGNS zone 保存，启动时再分发到新的叶块。

## 14. 统计稳定与物理验收

程序达到 t_end 只表示时间推进完成，不表示槽道湍流统计收敛。至少还需：

1. 从 statistics 中剔除初始过渡段；
2. 检查两个截面的流量时间序列无系统漂移；
3. 检查上下壁平均剪切和实测 \(Re_\tau\) 的长时平均及对称性；
4. 用实测 \(u_\tau\) 构造 \(y^+\) 和 \(u^+\)；
5. 对不同时间窗做分块平均，估计统计不确定度；
6. 使用更细网格做网格收敛研究；
7. 与公开 \(Re_\tau\approx180\) 基准比较平均速度、RMS 和 Reynolds 应力。

当前程序只直接输出瞬时空间统计，没有内建累积时间平均、RMS 或 Reynolds 应力。若需要这些
量，应保存足够密集的 statistics/field 数据后离线处理，或在下一开发阶段新增统计模块。

## 15. 常见故障

| 现象 | 原因 | 处理 |
|---|---|---|
| CMake 找不到 MPI_CXX | MPI 开发包或 module 未加载 | 检查 mpicxx，清空构建目录后重新配置 |
| libmpi.so not found | 运行节点没有相同 MPI 环境 | 在作业脚本重新 module load，检查 ldd |
| CGNS zip hash失败 | 传输损坏或文件被替换 | 用包级和内容级 SHA-256 重新校验 |
| output directory already exists | allow_existing=false 防覆盖 | 使用新 segment 名和新输出目录 |
| maximum_steps，退出 2 | 未达到 t_end 就触发硬上限 | 检查稳定性后增加步数并从 checkpoint 续算 |
| wall_time_checkpoint，退出 2 | 程序墙钟安全停止 | 建立新 segment，从 latest checkpoint 续算 |
| numerical_failure，退出 3 | 非有限量、正性或算子失败 | 保留全部日志，不自动重启；检查 CFL、初边值和网格 |
| rank 数过多/块过窄 | 违反最小模板宽度 | 减少 rank，先做 dry-run |
| 作业看似挂起 | MPI/文件系统/某 rank 异常 | 查看调度器状态和每节点日志，不要反复盲目重启 |
| 磁盘快速增长 | field/checkpoint 间隔过密 | 停止后建立新段，增大 every_time |

## 16. 每次长算前检查清单

- [ ] 压缩包和包内 SHA-256 全部通过；
- [ ] 编译器与 mpiexec 属于同一 MPI 环境；
- [ ] Release MPI 构建完成，ldd 无 not found；
- [ ] 4-rank 五步卡口通过；
- [ ] 目标 rank dry-run 通过；
- [ ] manifest 版本等于 WCNS_SOURCE_REVISION；
- [ ] case.name 和 output.directory 对本段唯一；
- [ ] 程序 max_wall_time 早于调度器墙钟；
- [ ] checkpoint 间隔和磁盘配额合理；
- [ ] 重启路径指向完整的 latest checkpoint；
- [ ] 没有修改重启签名中的物理或数值设置；
- [ ] 已记录编译模块、节点、rank、作业号和配置副本；
- [ ] 明确本稀疏网格不构成 DNS 定量验收。
