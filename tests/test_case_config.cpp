#include "test_support.hpp"

#include <wcns/runtime/case_config.hpp>

#include <string>

namespace {

std::string valid_config()
{
    return R"(
schema_version = 1
case.name = parser-test
mesh.path = mesh.cgns
algorithm.profile = phenglei_wcns
algorithm.reconstruction = weno_z
algorithm.reconstruction_variables = characteristic
algorithm.riemann = hllc
algorithm.mdcd.disp = 0.04
algorithm.mdcd.diss = 0.005
gas.gamma = 1.4
gas.molar_mass = 0.029
reference.velocity = 340
reference.density = 1.225
reference.temperature = 288.15
reference.length = 1
reference.viscosity = 1.7894e-5
partition.mode = auto_split
partition.allow_idle_ranks = false
partition.max_load_ratio = 1.2
partition.min_cells_per_active_direction = 8
initial.type = uniform
initial.rho = 1
initial.u = 0.2
initial.v = 0
initial.w = 0
initial.temperature = 1
boundary.default = farfield
boundary.wall.type = no_slip_isothermal_wall
boundary.wall.wall_velocity_x = 0.75
boundary.wall.wall_temperature = 1.1
boundary.inlet.type = inflow
boundary.inlet.rho = 1.05
boundary.inlet.u = 0.3
boundary.inlet.temperature = 0.95
source.enabled = false
run.mode = steady
run.viscous = false
run.cfl = 0.2
run.max_steps = 10
output.directory = output
output.allow_existing = false
output.dimensional = false
output.field.enabled = false
output.history.enabled = false
output.statistics.enabled = false
output.checkpoint.enabled = false
)";
}

} // namespace

void test_case_config()
{
    {
        const auto config = wcns::CaseConfig::from_text(valid_config());
        WCNS_REQUIRE(config.schema_version == 1);
        WCNS_REQUIRE(config.case_name == "parser-test");
        WCNS_REQUIRE(config.mesh_path == "mesh.cgns");
        WCNS_REQUIRE(config.profile == wcns::AlgorithmProfileKind::PhengleiWcns);
        WCNS_REQUIRE(config.reconstruction.scheme == "weno_z");
        WCNS_REQUIRE(
            config.reconstruction.variables
            == wcns::ReconstructionVariables::Characteristic);
        WCNS_REQUIRE(config.riemann.scheme == "hllc");
        WCNS_REQUIRE_NEAR(
            config.reconstruction.nonlinear.mdcd_dispersion, 0.04, 1.0e-15);
        WCNS_REQUIRE_NEAR(
            config.reconstruction.nonlinear.mdcd_dissipation, 0.005, 1.0e-15);
        WCNS_REQUIRE(config.restart_signature().find("mdcd_dispersion=0.04")
            != std::string::npos);
        WCNS_REQUIRE(config.partition.mode == wcns::PartitionMode::AutoSplit);
        WCNS_REQUIRE(
            config.boundary_overrides.at("wall")
            == wcns::BoundaryType::NoSlipIsothermalWall);
        WCNS_REQUIRE_NEAR(
            *config.boundary_data.at("wall").wall_velocity[0], 0.75, 1.0e-15);
        WCNS_REQUIRE_NEAR(
            *config.boundary_data.at("wall").wall_temperature, 1.1, 1.0e-15);
        WCNS_REQUIRE_NEAR(
            *config.boundary_data.at("inlet").rho, 1.05, 1.0e-15);
        WCNS_REQUIRE(config.restart_signature().find("wall_velocity_x=0.75")
            != std::string::npos);
        WCNS_REQUIRE(config.run.max_steps == 10);
        WCNS_REQUIRE(config.run.mode == wcns::RunMode::Steady);
        WCNS_REQUIRE(config.output.directory == "output");
        WCNS_REQUIRE(config.digest() != 0);
        WCNS_REQUIRE(config.summary().find("Re=") != std::string::npos);
        WCNS_REQUIRE(config.summary().find("Ma=") != std::string::npos);
    }
    {
        auto monitored = valid_config();
        const auto enabled = monitored.find("output.statistics.enabled = false");
        monitored.replace(
            enabled,
            std::string("output.statistics.enabled = false").size(),
            "output.statistics.enabled = true\n"
            "output.statistics.every_steps = 2\n"
            "output.statistics.quantities = total_mass\n"
            "output.statistics.xz_planes.enabled = true\n"
            "output.statistics.xz_planes.cell_j_indices = 0, 23, 47");
        const auto config = wcns::CaseConfig::from_text(monitored);
        WCNS_REQUIRE(config.output.xz_planes.enabled);
        WCNS_REQUIRE(config.output.xz_planes.cell_j_indices
            == std::vector<int>({0, 23, 47}));
        WCNS_REQUIRE(config.output.statistics.quantities
            == std::vector<std::string>({
                "total_mass", "xz_mean_u_j0", "xz_mass_flow_x_j0",
                "xz_mean_u_j23", "xz_mass_flow_x_j23",
                "xz_mean_u_j47", "xz_mass_flow_x_j47"}));

        auto duplicate = monitored;
        const auto indices = duplicate.find(
            "output.statistics.xz_planes.cell_j_indices = 0, 23, 47");
        duplicate.replace(
            indices,
            std::string(
                "output.statistics.xz_planes.cell_j_indices = 0, 23, 47").size(),
            "output.statistics.xz_planes.cell_j_indices = 0, 23, 23");
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(duplicate));
    }
    {
        auto unsteady = valid_config();
        const auto mode = unsteady.find("run.mode = steady");
        unsteady.replace(mode, std::string("run.mode = steady").size(),
            "run.mode = unsteady\nrun.t_end = 0.25");
        const auto config = wcns::CaseConfig::from_text(unsteady);
        WCNS_REQUIRE(config.run.mode == wcns::RunMode::Unsteady);
        WCNS_REQUIRE_NEAR(config.run.end_time, 0.25, 1.0e-15);

        auto missing_end = valid_config();
        const auto missing_mode = missing_end.find("run.mode = steady");
        missing_end.replace(
            missing_mode,
            std::string("run.mode = steady").size(),
            "run.mode = unsteady");
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(missing_end));
    }
    {
        auto reordered = valid_config();
        reordered = "# ignored\n" + reordered;
        const auto first = wcns::CaseConfig::from_text(valid_config());
        const auto second = wcns::CaseConfig::from_text(reordered);
        WCNS_REQUIRE(first.digest() == second.digest());
    }
    {
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(valid_config() + "run.cfl = 0.1\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(valid_config() + "unknown.field = 1\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(valid_config() + "Re = 1000\n"));
        auto invalid_mdcd = valid_config();
        const auto diss_position = invalid_mdcd.find("algorithm.mdcd.diss = 0.005");
        invalid_mdcd.replace(
            diss_position,
            std::string("algorithm.mdcd.diss = 0.005").size(),
            "algorithm.mdcd.diss = 0.04");
        WCNS_REQUIRE_THROWS(
            std::invalid_argument,
            wcns::CaseConfig::from_text(invalid_mdcd));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(
                valid_config() + "output.history.quantities = rho\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(
                valid_config() + "boundary.bad.rho = 1\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(
                valid_config()
                + "boundary.bad.type = inflow\n"
                  "boundary.bad.rho = 1\n"
                  "boundary.bad.temperature = 1\n"
                  "boundary.bad.pressure = 1\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(
                valid_config()
                + "boundary.bad.type = farfield\n"
                  "boundary.bad.wall_temperature = 1\n"));

        auto invalid_vortex_period = valid_config();
        const auto initial_type = invalid_vortex_period.find("initial.type = uniform");
        invalid_vortex_period.replace(
            initial_type,
            std::string("initial.type = uniform").size(),
            "initial.type = isentropic_vortex\ninitial.period_x = -10");
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(invalid_vortex_period));
    }
    {
        auto two_gas_inputs = valid_config()
            + "gas.specific_gas_constant = 287\n";
        WCNS_REQUIRE_THROWS(
            wcns::PhysicsConfigurationError,
            wcns::CaseConfig::from_text(two_gas_inputs));

        auto scmm = valid_config();
        const auto profile_position = scmm.find("phenglei_wcns");
        scmm.replace(
            profile_position,
            std::string("phenglei_wcns").size(),
            "scmm6_wcns");
        const auto minimum_position = scmm.find(
            "partition.min_cells_per_active_direction = 8");
        scmm.replace(
            minimum_position,
            std::string("partition.min_cells_per_active_direction = 8").size(),
            "partition.min_cells_per_active_direction = 4");
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(scmm));
    }
    {
        auto double_mach = valid_config();
        const auto initial = double_mach.find("initial.type = uniform");
        double_mach.replace(
            initial,
            std::string("initial.type = uniform").size(),
            "initial.type = double_mach_reflection\ninitial.x0 = 0.16666666666666667");
        double_mach += "boundary.double_mach.type = double_mach_reflection\n";
        const auto config = wcns::CaseConfig::from_text(double_mach);
        WCNS_REQUIRE(
            config.boundary_overrides.at("double_mach")
            == wcns::BoundaryType::DoubleMachReflection);
        WCNS_REQUIRE(config.restart_signature().find(
            "double_mach_reflection_v1") != std::string::npos);

        auto wrong_gamma = double_mach;
        const auto gamma = wrong_gamma.find("gas.gamma = 1.4");
        wrong_gamma.replace(
            gamma, std::string("gas.gamma = 1.4").size(), "gas.gamma = 1.3");
        WCNS_REQUIRE_THROWS(
            wcns::PhysicsConfigurationError,
            wcns::CaseConfig::from_text(wrong_gamma));
    }
}
