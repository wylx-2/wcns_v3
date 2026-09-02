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
boundary.wall.type = slip_wall
source.enabled = false
run.viscous = false
run.cfl = 0.2
run.max_steps = 10
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
        WCNS_REQUIRE(config.partition.mode == wcns::PartitionMode::AutoSplit);
        WCNS_REQUIRE(
            config.boundary_overrides.at("wall") == wcns::BoundaryType::SlipWall);
        WCNS_REQUIRE(config.run.max_steps == 10);
        WCNS_REQUIRE(config.digest() != 0);
        WCNS_REQUIRE(config.summary().find("Re=") != std::string::npos);
        WCNS_REQUIRE(config.summary().find("Ma=") != std::string::npos);
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
}
