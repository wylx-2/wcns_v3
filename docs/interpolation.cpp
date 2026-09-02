#include "ns3d_func.h"
#include <array>
#include <vector>
#include <algorithm>

// Stage-L algorithm reference only. This file is not compiled by wcns.
// The production implementation must use the checked strategy/registry contract
// documented in stage-l-design.md and must not copy the silent fallback below.

// wrapper that selects recon method

// Public wrapper that accepts a runtime-sized stencil (std::vector)
// and forwards to the fixed-size implementations. This historical reference
// assumes (but does not validate) a six-value stencil. The Stage-L production
// wrapper must reject a wrong stencil size before indexing it.
double interpolate_select(const std::vector<double> &vstencil, double flag, const SolverParams P)
{
    SolverParams::Interpolation r = P.interpolation;

    // helper to optionally reverse according to sign(flag)
    if (r == SolverParams::Interpolation::WENO5)
    {
        std::array<double, 6> a5;
        for (int i = 0; i < 6; ++i)
            a5[i] = vstencil[i];
        if (flag < 0.0)
            std::reverse(a5.begin(), a5.end());
        return weno5_interpolate(a5);
    }

    if (r == SolverParams::Interpolation::ZERO)
    {
        std::array<double, 6> a5;
        for (int i = 0; i < 6; ++i)
            a5[i] = vstencil[i];
        if (flag < 0.0)
            std::reverse(a5.begin(), a5.end());
        return zero_interpolate(a5);
    }

    if (r == SolverParams::Interpolation::MDCD_LINEAR)
    {
        std::array<double, 6> a5;
        for (int i = 0; i < 6; ++i)
            a5[i] = vstencil[i];
        if (flag < 0.0)
            std::reverse(a5.begin(), a5.end());
        return mdcd_linear_interpolate(a5, P);
    }

    if (r == SolverParams::Interpolation::MDCD_HYBRID)
    {
        std::array<double, 6> a5;
        for (int i = 0; i < 6; ++i)
            a5[i] = vstencil[i];
        if (flag < 0.0)
            std::reverse(a5.begin(), a5.end());
        return mdcd_hybrid_interpolate(a5, P);
    }
    // last resort: return the single value or 0
    return vstencil.empty() ? 0.0 : vstencil[0];
}

// 简单的零阶插值
double zero_interpolate(const std::array<double, 6> &stencil)
{
    return stencil[2]; // upwind point for left state, downwind point for right state
}

// WENO5 插值（标量，6点模板）
double weno5_interpolate(const std::array<double, 6> &stencil)
{
    double f0 = stencil[0];
    double f1 = stencil[1];
    double f2 = stencil[2];
    double f3 = stencil[3];
    double f4 = stencil[4];

    double eps = 1e-6;
    double alpha0 = 0.0, alpha1 = 0.0, alpha2 = 0.0;

    double beta0 = (13.0 / 12.0) * (f0 - 2 * f1 + f2) * (f0 - 2 * f1 + f2) + 0.25 * (f0 - 4 * f1 + 3 * f2) * (f0 - 4 * f1 + 3 * f2);
    double beta1 = (13.0 / 12.0) * (f1 - 2 * f2 + f3) * (f1 - 2 * f2 + f3) + 0.25 * (f1 - f3) * (f1 - f3);
    double beta2 = (13.0 / 12.0) * (f2 - 2 * f3 + f4) * (f2 - 2 * f3 + f4) + 0.25 * (3 * f2 - 4 * f3 + f4) * (3 * f2 - 4 * f3 + f4);

    alpha0 = (1.0 / 16.0) / ((eps + beta0) * (eps + beta0));
    alpha1 = (5.0 / 8.0) / ((eps + beta1) * (eps + beta1));
    alpha2 = (5.0 / 16.0) / ((eps + beta2) * (eps + beta2));

    double asum = alpha0 + alpha1 + alpha2;
    double w0 = alpha0 / asum;
    double w1 = alpha1 / asum;
    double w2 = alpha2 / asum;

    double p0 = (3 * f0 - 10 * f1 + 15 * f2) / 8.0;
    double p1 = (-f1 + 6 * f2 + 3 * f3) / 8.0;
    double p2 = (3 * f2 + 6 * f3 - f4) / 8.0;

    return w0 * p0 + w1 * p1 + w2 * p2;
}

// MDCD 插值
double mdcd_linear_interpolate(const std::array<double, 6> &stencil, const SolverParams P)
{
    double f0 = stencil[0];
    double f1 = stencil[1];
    double f2 = stencil[2];
    double f3 = stencil[3];
    double f4 = stencil[4];
    double f5 = stencil[5];
    double diss = P.mdcd_diss;
    double disp = P.mdcd_disp;

    // 线性格式
    return (3.0 * disp + 3.0 * diss) / 8.0 * f0 + (-18.0 * disp - 30 * diss - 1.0) / 16.0 * f1 +
           (12.0 * disp + 60.0 * diss + 9.0) / 16.0 * f2 + (12.0 * disp - 60.0 * diss + 9.0) / 16.0 * f3 +
           (-18.0 * disp + 30.0 * diss - 1.0) / 16.0 * f4 + (3.0 * disp - 3.0 * diss) / 8.0 * f5;
}

double mdcd_hybrid_interpolate(const std::array<double, 6> &stencil, const SolverParams P)
{
    double diss = P.mdcd_diss;
    double disp = P.mdcd_disp;
    double eps_small = 1e-40;

    double f0 = stencil[0];
    double f1 = stencil[1];
    double f2 = stencil[2];
    double f3 = stencil[3];
    double f4 = stencil[4];
    double f5 = stencil[5];

    /*---------------------------------------------
     * 1. Smooth-region detector (sigma)
    *--------------------------------------------*/
    double eps = 0.9 * 0.4 / (1.0 - 0.9 * 0.4) * 1e-4;

    double a1 = std::abs(f2 - f1) + std::abs(f2 - 2*f1 + f0);
    double b1 = std::abs(f2 - f3) + std::abs(f2 - 2*f3 + f4);
    double a2 = std::abs(f3 - f2) + std::abs(f3 - 2*f2 + f1);
    double b2 = std::abs(f3 - f4) + std::abs(f3 - 2*f4 + f5);
    double sai1 = (2*a1*b1 + eps) / (a1*a1 + b1*b1 + eps);
    double sai2 = (2*a2*b2 + eps) / (a2*a2 + b2*b2 + eps);
    double sai  = std::min(sai1, sai2);
    bool sigma = (sai > 0.4);

    /*---------------------------------------------
     * 2. MDCD interpolation
    *--------------------------------------------*/
    if (sigma)
        // sai 接近 1 表示两侧变化对称/光滑；光滑区使用线性 MDCD。
        return mdcd_linear_interpolate(stencil, P);
    else
    {
    /*---------------------------------------------
     * 3. Nonsmooth-region MDCD-WENO reconstruction
     *--------------------------------------------*/
        double q0 = (3 * f0 - 10 * f1 + 15 * f2) / 8.0;
        double q1 = (-f1 + 6 * f2 + 3 * f3) / 8.0;
        double q2 = (3 * f2 + 6 * f3 - f4) / 8.0;
        double q3 = (15 * f3 - 10 * f4 + 3 * f5) / 8.0;

        double d0 = 1.5 * (disp + diss);
        double d1 = 0.5 - 1.5 * (disp - 3*diss);
        double d2 = 0.5 - 1.5 * (disp + 3*diss);
        double d3 = 1.5 * (disp - diss);

        double beta0 =
            13.0/12.0*std::pow(f0-2*f1+f2,2) +
            0.25*std::pow(f0-4*f1+3*f2,2);
        double beta1 =
            13.0/12.0*std::pow(f1-2*f2+f3,2) +
            0.25*std::pow(f1-f3,2);
        double beta2 =
            13.0/12.0*std::pow(f2-2*f3+f4,2) +
            0.25*std::pow(3*f2-4*f3+f4,2);
        // 271779 is required for beta3 to vanish for a constant stencil.
        // The original supplied draft used 271799, a 20/120960 constant-state typo.
        double beta3 =
            (271779*std::pow(f0,2)
            + f0*(-2380800*f1+4086352*f2-3462252*f3+1458762*f4-245620*f5)
            + f1*(5653317*f1-20427884*f2+17905032*f3-7727988*f4+1325006*f5)
            + f2*(19510972*f2-35817664*f3+15929912*f4-2792660*f5)
            + f3*(17195652*f3-15880404*f4+2863984*f5)
            + f4*(3824847*f4-1429976*f5)
            + 139633*std::pow(f5,2)) / 120960.0;
        double tau6 = std::abs(beta3 - (beta0 + 4*beta1 + beta2)/6.0);

        double a0 = d0 * std::pow(20 + tau6/(beta0+eps_small),2);
        double a1 = d1 * std::pow(20 + tau6/(beta1+eps_small),2);
        double a2 = d2 * std::pow(20 + tau6/(beta2+eps_small),2);
        double a3 = d3 * std::pow(20 + tau6/(beta3+eps_small),2);

        return (a0*q0 + a1*q1 + a2*q2 + a3*q3) / (a0 + a1 + a2 + a3);
    }
}
