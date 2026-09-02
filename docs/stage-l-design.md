# 阶段 L：界面重构、Riemann 求解器与健壮性设计

状态：设计已由项目负责人批准；必做 L1--L5、L7 候选 v1 已完成自动验收，等待人工检查。
可选 L6 未获单独批准，未实施。详见 [`stage-l-acceptance.md`](stage-l-acceptance.md)。

本文是阶段 L 的实施前冻结规格。阶段 L 只扩展无粘界面算法和健壮性机制，不改变
`phenglei_wcns`、`scmm6_wcns` 各自绑定的度量、线性插值、通量差分、边界闭合和
`FaceFluxHalo` 语义。[`interpolation.cpp`](interpolation.cpp) 仅作为 MDCD 公式来源，不参与构建；其中
未检查模板长度、用浮点符号选择左右方向和未知枚举静默返回等行为不得复制到生产代码。

## 1. 配置与统一调用约定

阶段 L 必须提供以下规范名称：

```text
reconstruction = weno_js | weno_z | mdcd_linear | mdcd_hybrid
reconstruction_variables = conservative | primitive | characteristic
riemann_solver = rusanov | hllc | roe
low_mach_preconditioning = off | weiss_smith_roe
```

`linear5` 保留为测试、正性回退和光滑基线方案，但不替代上述四种必做重构。
早期名称 `wcns_js` 只能作为带弃用警告的输入迁移别名，内部配置摘要和重启签名一律写
`weno_js`。算法名称不区分大小写的便利转换只能发生在配置读取层；注册表键使用上述小写
规范名称。未知名称、重复注册、非法参数或不兼容组合必须在推进前报错。

每个方向的面法向指向逻辑坐标增加方向，左状态来自低索引侧，右状态来自高索引侧。
重构器只返回状态，不计算通量；Riemann 求解器只接收单位法向并返回每单位物理面积的
通量，调用方唯一地乘面面积。所有内建方案均复用 `ReconstructionScaling`、
`NumericalFloors` 和同一套逐级回退/诊断设施。

## 2. 重构策略接口与扩展方法

生产接口采用显式策略和注册表，语义等价于：

```cpp
enum class TraceSide { Left, Right };

struct ReconstructionContext {
    ReconstructionScaling scaling;
    ReconstructionParameters parameters;
};

class IReconstructionScheme {
public:
    virtual std::string_view name() const noexcept = 0;
    virtual StencilRequirement stencil_requirement() const noexcept = 0;
    virtual double reconstruct_scalar(ScalarStencilView stencil,
                                      TraceSide side,
                                      const ReconstructionContext& context) const = 0;
    virtual ~IReconstructionScheme() = default;
};
```

`ScalarStencilView` 携带值的数量和逻辑偏移；本阶段四个方案统一要求六点
`[q_{j-2},q_{j-1},q_j,q_{j+1},q_{j+2},q_{j+3}]`。入口必须检查长度、偏移和全部值有限，
不得越界读取或把错误模板降级成任意单元值。`TraceSide::Right` 通过反转六点模板后调用同一
左状态核得到，不使用 `double flag` 的正负号表达方向。

新增自定义重构算法的固定步骤为：

1. 实现 `IReconstructionScheme`，声明唯一规范名称、左右所需半径和精确模板布局。
2. 定义参数、默认值、合法域、配置序列化和重启签名；禁止在算法内部隐藏常数。
3. 在 `ReconstructionRegistry` 注册工厂。调用代码只依赖接口和注册表，不增加中央
   `switch`，也不修改已有算法。
4. 若所需 halo 超过当前三层，必须先扩展并验收状态交换和边界闭合；不能读取不存在的
   物理边界边/角 ghost 或任何 ghost 几何量。
5. 增加常数保持、声明阶数、左右镜像、尺度不变性、非有限输入、正性回退、物理边界、
   多块和 MPI 一致性测试，并补充算法公式与限制。

本阶段注册表是稳定的源码级扩展点，不承诺跨编译器的动态库 ABI；若后续需要运行期加载
外部插件，必须另行设计 ABI 版本、所有权、异常边界和 MPI 各 rank 插件一致性协议。

## 3. 四种必做界面重构

以下均给出标量左状态。令

$$
(f_0,f_1,f_2,f_3,f_4,f_5)
=(q_{j-2},q_{j-1},q_j,q_{j+1},q_{j+2},q_{j+3}),
$$

三个三点候选值为

$$
q_0=\frac{3f_0-10f_1+15f_2}{8},\quad
q_1=\frac{-f_1+6f_2+3f_3}{8},\quad
q_2=\frac{3f_2+6f_3-f_4}{8},
$$

线性权为 $(d_0,d_1,d_2)=(1,10,5)/16$。对应线性组合正好是

$$
q^{L}_{j+1/2}=\frac{3f_0-20f_1+90f_2+60f_3-5f_4}{128}.
$$

令 $s=\max(s_q,s_{floor})$，所有光滑度指标用无量纲值
$\bar f_m=f_m/s$ 计算。这样同一个严格正的无量纲 `reconstruction_epsilon`
$\epsilon_w$ 可以用于不同物理量和特征分量。

### 3.1 WENO-JS

$$
\begin{aligned}
\bar\beta_0&=\frac{13}{12}(\bar f_0-2\bar f_1+\bar f_2)^2
+\frac14(\bar f_0-4\bar f_1+3\bar f_2)^2,\\
\bar\beta_1&=\frac{13}{12}(\bar f_1-2\bar f_2+\bar f_3)^2
+\frac14(\bar f_1-\bar f_3)^2,\\
\bar\beta_2&=\frac{13}{12}(\bar f_2-2\bar f_3+\bar f_4)^2
+\frac14(3\bar f_2-4\bar f_3+\bar f_4)^2,
\end{aligned}
$$

$$
\alpha_k=\frac{d_k}{(\epsilon_w+\bar\beta_k)^{p_{JS}}},\qquad
\omega_k=\frac{\alpha_k}{\sum_l\alpha_l},\qquad
q_L=\sum_k\omega_kq_k,
$$

其中项目默认 $p_{JS}=2$。

### 3.2 WENO-Z

WENO-Z 复用上述候选值和光滑度指标，定义

$$
\tau_5=|\bar\beta_0-\bar\beta_2|,\qquad
\alpha_k=d_k\left[1+
\left(\frac{\tau_5}{\epsilon_w+\bar\beta_k}\right)^{p_Z}\right],
\qquad q_L=\frac{\sum_k\alpha_kq_k}{\sum_k\alpha_k},
$$

项目默认 $p_Z=2$。WENO-Z 不能复用 WENO-JS 的 $d_k/(\epsilon+\beta_k)^p$ 权重后
只更换名称。

### 3.3 MDCD_LINEAR

`docs/interpolation.cpp` 中的六点线性 MDCD 公式冻结为

$$
q_L=\sum_{m=0}^{5}c_mf_m,
$$

$$
\begin{aligned}
c_0&=\frac{3(\gamma_{disp}+\gamma_{diss})}{8},&
c_1&=\frac{-18\gamma_{disp}-30\gamma_{diss}-1}{16},\\
c_2&=\frac{12\gamma_{disp}+60\gamma_{diss}+9}{16},&
c_3&=\frac{12\gamma_{disp}-60\gamma_{diss}+9}{16},\\
c_4&=\frac{-18\gamma_{disp}+30\gamma_{diss}-1}{16},&
c_5&=\frac{3(\gamma_{disp}-\gamma_{diss})}{8}.
\end{aligned}
$$

项目默认 `mdcd_dispersion=0.0463783`、`mdcd_dissipation=0.01`。二者必须有限，且满足

$$
\gamma_{disp}>0,\qquad 0\le\gamma_{diss}<\gamma_{disp},\qquad
3\gamma_{disp}+9\gamma_{diss}<1,
$$

从而下一节四个最优权均严格为正。参数进入日志、配置摘要和重启签名。

### 3.4 MDCD_HYBRID

混合方案先在无量纲 $\bar f_m$ 上计算光滑区传感器：

$$
\begin{aligned}
a_1&=|\bar f_2-\bar f_1|+|\bar f_2-2\bar f_1+\bar f_0|,&
b_1&=|\bar f_2-\bar f_3|+|\bar f_2-2\bar f_3+\bar f_4|,\\
a_2&=|\bar f_3-\bar f_2|+|\bar f_3-2\bar f_2+\bar f_1|,&
b_2&=|\bar f_3-\bar f_4|+|\bar f_3-2\bar f_4+\bar f_5|,\\
\psi_r&=\frac{2a_rb_r+\epsilon_s}{a_r^2+b_r^2+\epsilon_s},&
\psi&=\min(\psi_1,\psi_2).
\end{aligned}
$$

冻结 `docs/interpolation.cpp` 的默认值

$$
\epsilon_s=\frac{0.9\times0.4}{1-0.9\times0.4}\,10^{-4}=5.625\times10^{-5},
\qquad \psi_{threshold}=0.4.
$$

当 $\psi>\psi_{threshold}$ 时判为光滑区并调用 `mdcd_linear`；否则进入下述非线性
MDCD-WENO 分支。这个布尔量是“光滑区开关”，不得反向解释为间断开关。

非光滑分支增加候选

$$
q_3=\frac{15f_3-10f_4+3f_5}{8},
$$

参考文件非光滑分支的四个权为

$$
\begin{aligned}
d_0&=1.5(\gamma_{disp}+\gamma_{diss}),&
d_1&=0.5-1.5(\gamma_{disp}-3\gamma_{diss}),\\
d_2&=0.5-1.5(\gamma_{disp}+3\gamma_{diss}),&
d_3&=1.5(\gamma_{disp}-\gamma_{diss}).
\end{aligned}
$$

$\bar\beta_{0,1,2}$ 与 WENO-JS 相同；六点指标严格采用参考实现的二次型

$$
\begin{aligned}
120960\bar\beta_3={}&271779\bar f_0^2
+\bar f_0(-2380800\bar f_1+4086352\bar f_2-3462252\bar f_3
+1458762\bar f_4-245620\bar f_5)\\
&+\bar f_1(5653317\bar f_1-20427884\bar f_2+17905032\bar f_3
-7727988\bar f_4+1325006\bar f_5)\\
&+\bar f_2(19510972\bar f_2-35817664\bar f_3+15929912\bar f_4
-2792660\bar f_5)\\
&+\bar f_3(17195652\bar f_3-15880404\bar f_4+2863984\bar f_5)\\
&+\bar f_4(3824847\bar f_4-1429976\bar f_5)+139633\bar f_5^2.
\end{aligned}
$$

最后计算

$$
\tau_6=\left|\bar\beta_3-\frac{\bar\beta_0+4\bar\beta_1+\bar\beta_2}{6}\right|,
\quad
\alpha_k=d_k\left(C_{MDCD}+\frac{\tau_6}{\bar\beta_k+\epsilon_{MDCD}}\right)^2,
\quad
q_L=\frac{\sum_{k=0}^3\alpha_kq_k}{\sum_{k=0}^3\alpha_k}.
$$

随项目收到的参考文件曾把首项误写为 271799，使常数模板得到
$\bar\beta_3=20/120960$；这里和参考文件现均修正为 271779，从而常数模板严格得到零且
二次型为半正定。参考值冻结为 $C_{MDCD}=20$、$\epsilon_{MDCD}=10^{-40}$。它们作用于已经无量纲化的
指标，不与 WENO 的 $\epsilon_w$ 混用。实现必须检查每个 $\bar\beta_k$、$\alpha_k$ 和权重和
有限且分母为正；失败时走统一的正性/非有限回退链，不能在 MDCD 内部静默返回单点值。

还需注意，参考文件上述四权与候选值组合后，一般不等于第 3.3 节独立
`MDCD_LINEAR` 的六点系数（例如 $f_0$ 系数分别为
$9(\gamma_{disp}+\gamma_{diss})/16$ 和 $3(\gamma_{disp}+\gamma_{diss})/8$）。因此本项目按
参考文件把它们定义为两条独立路径：光滑区显式调用第 3.3 节线性核，非光滑区逐项采用本节
四权和 $\bar\beta_3$；不得把非光滑分支重构成一个声称与 `MDCD_LINEAR` 具有相同线性极限的
统一权重公式。若以后要改写这组权，须作为算法规格变更单独论证、测试和审批。

## 4. 特征重构

`reconstruction_variables=characteristic` 时，每个面先用该面单位法向和 Roe 平均状态建立
同一套法向 Euler 左右特征矩阵 $L_n,R_n$，并验证 $L_nR_n=I$。选择确定性的正交切向
$\boldsymbol t_1,\boldsymbol t_2$ 后，把整个六点模板用同一 $L_n$ 投影，分别调用所选标量
重构器，再用 $R_n$ 还原守恒状态。不得为模板内每个单元建立不同特征基，也不得把物理通量
Jacobian $A_n$ 与网格 Jacobian $J$ 混称。

Roe 平均无效、特征矩阵病态、还原结果非有限或不满足正密度/正压力时，按
“所选特征重构 -> 同方案 primitive 重构 -> `linear5` -> 相邻单元一阶状态”确定性回退，
记录面、方向、块、rank、RK 子步和原因。

## 5. Riemann 策略接口与扩展方法

生产接口语义等价于：

```cpp
struct RiemannResult {
    Conservative flux_per_unit_area;
    double spectral_radius;
    RiemannDiagnostics diagnostics;
};

class IRiemannSolver {
public:
    virtual std::string_view name() const noexcept = 0;
    virtual RiemannResult solve(const Conservative& left,
                                const Conservative& right,
                                const UnitNormal& normal,
                                const RiemannContext& context) const = 0;
    virtual ~IRiemannSolver() = default;
};
```

入口检查法向模长在容差内为 1、左右状态有限且满足 floors。结果必须给出通量、用于 CFL 的
谱半径、实际使用的求解器以及回退原因。新增自定义 Riemann 求解器时必须：实现此接口；
定义并校验参数及重启签名；在 `RiemannSolverRegistry` 注册唯一键；声明可兼容的预处理器；
加入一致性 $\mathcal R(U,U,n)=F(U)\cdot n$、反向法向/交换左右、接触与激波、非物理状态
回退、多块共享面唯一所有者和 MPI 一致性测试。调用方不增加中央 `switch`。

## 6. 三种必做 Riemann 求解器

以下 $F_n(U)=F(U)\cdot\boldsymbol n$，$u_n=\boldsymbol u\cdot\boldsymbol n$，
$a=\sqrt{\gamma p/\rho}$。

### 6.1 Rusanov

$$
\mathcal R_{Rus}=\frac12\left[F_n(U_L)+F_n(U_R)
-\alpha(U_R-U_L)\right],\qquad
\alpha=\max(|u_{n,L}|+a_L,|u_{n,R}|+a_R).
$$

Rusanov 是所有重构和 Riemann 异常的最终确定性回退；如果连相邻单元一阶状态也不合法，
立即终止并报告源单元，不得继续推进。

### 6.2 HLLC

以平方根密度权得到 Roe 平均 $\widetilde u_n,\widetilde a$，波速采用

$$
S_L=\min(u_{n,L}-a_L,\widetilde u_n-\widetilde a),\qquad
S_R=\max(u_{n,R}+a_R,\widetilde u_n+\widetilde a),
$$

$$
S_M=\frac{p_R-p_L+\rho_Lu_{n,L}(S_L-u_{n,L})
-\rho_Ru_{n,R}(S_R-u_{n,R})}
{\rho_L(S_L-u_{n,L})-\rho_R(S_R-u_{n,R})},
$$

$$
p_*=p_K+\rho_K(S_K-u_{n,K})(S_M-u_{n,K}),\quad K\in\{L,R\}.
$$

星区状态为

$$
\rho_{*K}=\rho_K\frac{S_K-u_{n,K}}{S_K-S_M},\quad
\boldsymbol u_{*K}=\boldsymbol u_K+(S_M-u_{n,K})\boldsymbol n,
$$

$$
(\rho E)_{*K}=
\frac{(S_K-u_{n,K})(\rho E)_K-p_Ku_{n,K}+p_*S_M}{S_K-S_M}.
$$

按 $0\le S_L$、$S_L\le0\le S_M$、$S_M\le0\le S_R$、$S_R\le0$ 四段分别返回
$F_L$、$F_L+S_L(U_{*L}-U_L)$、$F_R+S_R(U_{*R}-U_R)$、$F_R$。
波速、任一分母、星区状态或最终通量无效时回退到 Rusanov并计数。

### 6.3 Roe

令 $r_K=\sqrt{\rho_K}$，采用

$$
\widetilde{\boldsymbol u}=\frac{r_L\boldsymbol u_L+r_R\boldsymbol u_R}{r_L+r_R},\qquad
\widetilde H=\frac{r_LH_L+r_RH_R}{r_L+r_R},\qquad
\widetilde a^2=(\gamma-1)\left(\widetilde H-\frac12|\widetilde{\boldsymbol u}|^2\right).
$$

在与特征重构相同的法向/切向正交基中，使用标准五波 Euler 特征分解

$$
\mathcal R_{Roe}=\frac12(F_L+F_R)
-\frac12\sum_{m=1}^{5}|\lambda_m|_{fix}\,\alpha_m\boldsymbol r_m,
$$

其中特征值为
$(\widetilde u_n-\widetilde a,\widetilde u_n,\widetilde u_n,
\widetilde u_n,\widetilde u_n+\widetilde a)$，波强度由
$U_R-U_L=\sum_m\alpha_m\boldsymbol r_m$ 唯一确定。默认启用 Harten 型熵修正：

$$
|\lambda|_{fix}=
\begin{cases}
|\lambda|,&|\lambda|\ge\delta,\\
\frac12(\lambda^2/\delta+\delta),&|\lambda|<\delta,
\end{cases}\qquad
\delta=C_{entropy}\max(a_L,a_R,\widetilde a),
$$

项目默认 $C_{entropy}=0.1$。Roe 平均、声速平方、特征分解、预测更新或通量无效时，
按 `Roe -> HLLC -> Rusanov` 回退并分别计数。

## 7. 低 Mach 预处理开关

默认 `low_mach_preconditioning=off`，此时三种求解器必须与经典公式逐位一致。阶段 L 预留
`IPreconditioningModel`，可选实现 `weiss_smith_roe`；它只与 `riemann_solver=roe` 和稳态
伪时间模式组合。开启后必须一致地修改伪时间导数预处理矩阵、Roe 特征系统、谱半径/CFL 和
远场特征边界，不能只在某个公式里替换声速。

局部截止量采用与 PHengLEI 思路一致的有界形式

$$
\theta_c=\min\left(1,\max(M_c^2,K_{prec}Ma_{ref}^2)\right),\qquad
0<K_{prec}\le1,
$$

面上使用左右 $\theta_c$ 的算术平均。具体 Weiss--Smith 矩阵、特征向量和双时间项必须作为
该可选子功能的同一提交给出并单独验收。当前物理时间 SSPRK3 若请求非 `off`，或 Rusanov/
HLLC 与 `weiss_smith_roe` 组合，必须在启动时拒绝；不得把预处理后的伪时间解冒充正确的
非定常物理解。预处理是否纳入阶段 L 候选由项目负责人单独批准，不影响三种经典求解器的
必做验收。

## 8. 回退、诊断与实施顺序

一次面的确定性处理链为：校验输入 -> 选择变量空间并重构 -> 检查重构状态 -> 必要时按第
4 节降阶 -> 调用所选 Riemann 求解器 -> 必要时按第 6 节回退 -> 面所有者发布唯一通量。
诊断至少按算法、回退前后方案、原因、方向、块和 rank 计数，并可归约为全局日志；计数本身
不改变浮点计算顺序。物理边界的重构后强约束仍由 `strong_boundary_face_state` 唯一控制。

阶段 L 实施按 [`development-roadmap.md`](development-roadmap.md) 的小任务和卡口进行。
项目负责人已批准本设计并启动实现；候选 v1 自动验收完成后仍须等待人工放行。阶段 K 的
人工检查被暂缓，不等同于验收通过，也不改变其候选状态。

## 9. 主要算法依据

- Jiang--Shu WENO：G.-S. Jiang, C.-W. Shu, *Efficient Implementation of Weighted ENO
  Schemes*, JCP 126 (1996), DOI [10.1006/jcph.1996.0130](https://doi.org/10.1006/jcph.1996.0130)。
- WENO-Z：R. Borges et al., *An improved weighted essentially non-oscillatory scheme for
  hyperbolic conservation laws*, JCP 227 (2008), DOI
  [10.1016/j.jcp.2007.11.038](https://doi.org/10.1016/j.jcp.2007.11.038)。
- MDCD：Z.-S. Sun et al., *A class of finite difference schemes with low dispersion and
  controllable dissipation for DNS of compressible turbulence*, JCP 230 (2011),
  [期刊页面](https://www.sciencedirect.com/science/article/pii/S0021999111001276)；本项目具体混合
  公式以随仓库提供的 `interpolation.cpp` 为直接冻结参考。
- HLLC 波速：P. Batten et al., *On the Choice of Wavespeeds for the HLLC Riemann Solver*,
  SIAM J. Sci. Comput. 18 (1997), DOI
  [10.1137/S1064827593260140](https://doi.org/10.1137/S1064827593260140)。
- Roe：P. L. Roe, *Approximate Riemann Solvers, Parameter Vectors, and Difference Schemes*,
  JCP 43 (1981), DOI
  [10.1016/0021-9991(81)90128-5](https://doi.org/10.1016/0021-9991(81)90128-5)。
- 低 Mach 可选预处理：J. M. Weiss, W. A. Smith, *Preconditioning Applied to Variable and
  Constant Density Flows*, AIAA Journal 33 (1995), DOI
  [10.2514/3.12946](https://doi.org/10.2514/3.12946)。
