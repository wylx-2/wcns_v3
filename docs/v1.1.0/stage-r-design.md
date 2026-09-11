# v1.1.0 阶段 R 设计冻结：真实边界面工程量

状态：2026-09-11 根据 R→S→T 连续推进授权冻结；自动卡口通过后直接合入并进入 S。

## 1. 目标与非目标

本阶段实现《算法补充》11.4：从求解器真实边界迹计算逐面壁压、温度、黏度、压力/黏性/总
牵引、进入壁面的热流、`Cp/Cf`，再按命名 patch 积分力和力矩。它替换 Case07 第一层单元估算，
不改变 Euler/NS 方程、边界条件、输运模型或数值 restart signature，也不实现并行 CGNS、任意
曲面插值或非物理连接面积分。

## 2. 配置冻结

新增可选输出块，缺省完全关闭：

```text
output.boundary.enabled = false
output.boundary.format = txt | tecplot
output.boundary.every_steps = 0
output.boundary.every_time = 0
output.boundary.explicit_times =
output.boundary.write_initial = false
output.boundary.write_final = true
output.boundary.patches = cylinder
output.boundary.quantities = p_w,T_w,mu_w,Cp,Cf,q_wall,...
output.boundary.reference_pressure = <required>
output.boundary.reference_density = <required>
output.boundary.reference_velocity_x/y/z = <required>
output.boundary.reference_area = <required>
output.boundary.reference_length = <required>
output.boundary.moment_center_x/y/z = 0
output.boundary.drag_direction_x/y/z = <required unit vector>
output.boundary.lift_direction_x/y/z = <required unit vector>
output.boundary.tangent_direction_x/y/z = <required unit vector>
```

`patches` 和 `quantities` 使用逗号分隔且禁止重复。启用时参考密度、面积、长度为有限正数，参考
速度有限且产生大于压力 floor 的动压；三个方向归一化，drag/lift 正交。调度和上述参考值只
影响输出，进入完整 config digest、summary 和 manifest，但不进入 restart signature。

## 3. 权威数据流

1. 每次 accepted state 的 residual refresh 使用现有 profile、公共 ghost、梯度 halo 和
   `ViscousBoundaryTrace` 形成边界迹；求解器保存只读 `BoundaryFaceSample`。
2. 样本包含 leaf block、patch、方向/侧、局部面索引、面心、外法向、壁面 primitive、
   $\tau\boldsymbol n$、$\chi\nabla T\cdot\boldsymbol n$ 和是否具有黏性迹。
3. 输出器从 `GlobalConservationWeights` 取 $w_b$，计算 $\Delta A=w_b|S|$；从 partition leaf
   恢复 source zone 及原面索引。连接和周期面没有物理 patch，天然不进入样本。
4. 所有 rank 把固定宽度数值 payload 收集到 root；root 按
   `(patch,source_zone,k,j,i,axis,side)` 排序、查重并计算逐面量与 patch/load 总量。
5. 每个事件写一个逐面文件；同一事务追加一行固定 schema 的 load history。最终文件原子提交，
   与既有覆盖和 manifest 规则一致。

## 4. 逐面量与载荷 schema

内建量固定为 `x,y,z,area,nx,ny,nz,p_w,T_w,mu_w,Cp,pressure_traction_[xyz],
viscous_traction_[xyz],traction_[xyz],Cf,q_wall`。无粘壁允许几何、`p_w/T_w/Cp` 和压力牵引；
请求黏性量必须在初始化时拒绝。载荷历史始终保存 pressure/viscous/total 三组 $F_x,F_y,F_z$、
$M_x,M_y,M_z$，以及 `Cd_pressure/Cd_viscous/Cd_total`、`Cl_*`、`Cm_[xyz]`。二维文件头明确
单位展向长度。

## 5. 自动验收阈值

- 单元公式：平面单位面上给定 $p,\tau,\nabla T$ 的逐面量相对误差不超过 `1e-13`；法向翻转时
  压力力、黏性力、热流符号按公式翻转。
- 闭合曲面常压净力与总面积尺度之比不超过 `1e-12`；压力/黏性分量逐面相加与总量差不超过
  `1e-13`。
- Couette 剪切和线性导热壁热流相对解析误差不超过 `5e-11`；绝热热流绝对值不超过 `1e-13`。
- 同一多块算例 1/2/4 rank 的逐面排序键、载荷和力矩差不超过 `1e-12`；运行时切分不得改变
  source zone/patch/全局索引或重复物理面。
- 无量纲/有量纲缩放往返相对误差不超过 `1e-12`；非法引用、未知 patch/quantity、无粘请求
  黏性量均在正常输出前失败。
- Case07 生成真实 `Cp` 和压力载荷历史，配置及文档删除对第一层单元载荷作为生产依据的依赖。
- 公式核验、Release 串行/MPI 全量 CTest、`git diff --check` 和工作树审计全部通过。

## 6. Git 与阶段边界

R 的设计、核心 quantity、MPI/输出接入、测试/Case07 证据和自动验收分别提交。自动验收成功
创建 `v1.1.0-r-candidate.N`，按连续授权以 `--no-ff` 合入 `release/v1.1.0` 并立即建立 S
分支；不创建虚构的人工批准标签。远程不可写时在验收报告明确记录本地状态。
