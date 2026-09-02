#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/runtime/case_config.hpp>
#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/runtime/structured_partition.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/viscous_wcns_solver.hpp>

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

struct CommandLine {
    std::string config_path;
    bool dry_run = false;
};

CommandLine parse_command_line(int argc, char** argv)
{
    CommandLine result;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--dry-run") {
            result.dry_run = true;
        } else if (argument == "--config" && index + 1 < argc) {
            result.config_path = argv[++index];
        } else {
            throw std::invalid_argument(
                "usage: wcns_run --config <case.wcns> [--dry-run]");
        }
    }
    if (result.config_path.empty()) {
        throw std::invalid_argument(
            "usage: wcns_run --config <case.wcns> [--dry-run]");
    }
    return result;
}

std::string read_text(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open configuration: " + path);
    std::ostringstream result;
    result << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed to read configuration: " + path);
    }
    return result.str();
}

bool is_absolute_path(const std::string& path)
{
    if (path.empty()) return false;
    if (path.front() == '/' || path.front() == '\\') return true;
    return path.size() >= 3
        && std::isalpha(static_cast<unsigned char>(path[0])) != 0
        && path[1] == ':'
        && (path[2] == '/' || path[2] == '\\');
}

std::string resolve_mesh_path(
    const std::string& config_path,
    const std::string& mesh_path)
{
    if (is_absolute_path(mesh_path)) return mesh_path;
    const auto separator = config_path.find_last_of("/\\");
    if (separator == std::string::npos) return mesh_path;
    return config_path.substr(0, separator + 1) + mesh_path;
}

std::string broadcast_configuration(
    const wcns::MpiRuntime& mpi,
    const std::string& path)
{
    std::string payload;
    if (mpi.rank() == 0) {
        try {
            payload = "OK\n" + read_text(path);
        } catch (const std::exception& error) {
            payload = "ERROR\n" + std::string(error.what());
        }
    }
    payload = mpi.broadcast_string(std::move(payload));
    if (payload.rfind("ERROR\n", 0) == 0) {
        throw std::runtime_error(payload.substr(6));
    }
    if (payload.rfind("OK\n", 0) != 0) {
        throw std::runtime_error("invalid broadcast configuration envelope");
    }
    return payload.substr(3);
}

std::vector<wcns::PartitionZone> partition_zones(
    const wcns::CgnsMeshMetadata& metadata)
{
    std::vector<wcns::PartitionZone> result;
    result.reserve(metadata.zones.size());
    for (const auto& zone : metadata.zones) {
        result.push_back({
            zone.block_id,
            zone.name,
            zone.cell_dimension,
            zone.cell_extent,
        });
    }
    return result;
}

std::vector<wcns::CgnsPartitionLeaf> cgns_leaves(
    const wcns::StructuredPartitionPlan& plan)
{
    std::vector<wcns::CgnsPartitionLeaf> result;
    result.reserve(plan.leaves().size());
    for (const auto& leaf : plan.leaves()) {
        result.push_back({
            leaf.block,
            leaf.source_zone,
            leaf.cells.begin,
            leaf.cells.end,
            leaf.owner,
        });
    }
    return result;
}

wcns::BoundaryType configured_boundary_type(
    const wcns::CaseConfig& config,
    const std::string& name)
{
    const auto iterator = config.boundary_overrides.find(name);
    return iterator == config.boundary_overrides.end()
        ? config.default_boundary : iterator->second;
}

void configure_boundaries(
    wcns::StructuredMesh& global_mesh,
    std::vector<wcns::StructuredBlock>& local_blocks,
    const wcns::CaseConfig& config)
{
    for (const auto& descriptor : global_mesh.blocks()) {
        auto& block = global_mesh.block(descriptor.id());
        for (auto& patch : block.boundaries) {
            patch.type = configured_boundary_type(config, patch.name);
        }
    }
    for (auto& block : local_blocks) {
        for (auto& patch : block.boundaries) {
            patch.type = configured_boundary_type(config, patch.name);
        }
    }
}

wcns::BlockBoundaryDataMap make_boundary_data(
    const wcns::LocalBlockSet& local_blocks,
    const wcns::CaseConfig& config,
    const wcns::GasModel& gas,
    const wcns::ReferenceScales& reference,
    const wcns::NumericalFloors& floors)
{
    wcns::BlockBoundaryDataMap result;
    for (const auto& block : local_blocks.blocks()) {
        const auto target = wcns::FlowInitializer::evaluate(
            config.initial,
            {0.0, 0.0, 0.0},
            gas,
            reference,
            floors,
            block.cell_dimension());
        wcns::BoundaryDataMap data;
        for (const auto& patch : block.boundaries) {
            wcns::BoundaryData patch_data;
            if (patch.type == wcns::BoundaryType::Farfield
                || patch.type == wcns::BoundaryType::Inflow) {
                patch_data.target_state = target;
            }
            if (patch.type == wcns::BoundaryType::NoSlipIsothermalWall) {
                patch_data.wall_temperature = config.initial.parameter(
                    "temperature", 1.0);
            }
            patch_data.validate(patch.type, block.cell_dimension());
            data.emplace(patch.name, patch_data);
        }
        result.emplace(block.id(), std::move(data));
    }
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const auto command = parse_command_line(argc, argv);
        wcns::MpiRuntime mpi(argc, argv);
        const auto text = broadcast_configuration(mpi, command.config_path);
        auto config = wcns::CaseConfig::from_text(text);
        if (!mpi.all_equal(config.digest())) {
            throw std::runtime_error(
                "case configuration digest differs across MPI ranks");
        }

        const auto mesh_name = resolve_mesh_path(
            command.config_path,
            config.mesh_path);
        wcns::CgnsReader reader;
        const auto metadata = reader.read_metadata(mesh_name);
        const auto plan = wcns::StructuredPartitionPlan::build(
            partition_zones(metadata),
            mpi.size(),
            config.partition);
        if (!mpi.all_equal(plan.digest())) {
            throw std::runtime_error(
                "structured partition digest differs across MPI ranks");
        }
        auto partitioned = reader.read_partitioned_mesh(
            mesh_name,
            cgns_leaves(plan),
            mpi.rank(),
            3);
        configure_boundaries(
            partitioned.global_mesh,
            partitioned.local_blocks,
            config);
        const auto topology = wcns::DistributedTopology::build(
            partitioned.global_mesh,
            plan.distribution(),
            false);
        wcns::LocalBlockSet local_blocks(
            mpi.rank(),
            std::move(partitioned.local_blocks),
            plan.distribution());

        const auto gas = config.make_gas_model();
        const auto reference = config.make_reference_scales(gas);
        const auto profile = config.make_profile();
        const wcns::NumericalFloors floors;
        wcns::BlockMetricMap metrics;
        for (auto& block : local_blocks.blocks()) {
            metrics.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block.id()),
                std::forward_as_tuple(
                    wcns::initialize_metric_field(block, profile).metric));
        }
        wcns::FlowInitializer::initialize_local_blocks(
            local_blocks,
            metrics,
            config.initial,
            gas,
            reference,
            floors);
        const auto boundary_data = make_boundary_data(
            local_blocks,
            config,
            gas,
            reference,
            floors);

        if (mpi.rank() == 0) {
            std::cout << config.summary() << '\n'
                      << plan.summary() << '\n'
                      << "derived Re=" << std::setprecision(17)
                      << reference.reynolds()
                      << " Ma=" << reference.mach() << '\n';
        }
        if (command.dry_run) {
            if (mpi.rank() == 0) {
                std::cout << "WCNS dry-run completed\n";
            }
            return EXIT_SUCCESS;
        }

        wcns::Real time = 0.0;
        if (config.run.viscous) {
            wcns::ViscousWcnsConfig solver_config;
            solver_config.inviscid = config.make_inviscid_config();
            wcns::ViscousWcnsSolver solver(
                mpi,
                local_blocks,
                partitioned.global_mesh,
                topology,
                plan.distribution().rank_count(),
                metrics,
                boundary_data,
                profile,
                gas,
                reference,
                floors,
                solver_config);
            for (std::size_t step = 1; step <= config.run.max_steps; ++step) {
                const auto time_step = solver.global_time_step(config.run.cfl);
                solver.advance(time_step, time);
                time += time_step;
                if (mpi.rank() == 0) {
                    std::cout << "step=" << step << " time="
                              << std::setprecision(17) << time
                              << " dt=" << time_step
                              << " residual=" << solver.global_residual_l2()
                              << '\n';
                } else {
                    static_cast<void>(solver.global_residual_l2());
                }
            }
        } else {
            auto solver_config = config.make_inviscid_config();
            wcns::InviscidWcnsSolver solver(
                mpi,
                local_blocks,
                partitioned.global_mesh,
                topology,
                plan.distribution().rank_count(),
                metrics,
                boundary_data,
                profile,
                gas,
                reference,
                floors,
                solver_config);
            for (std::size_t step = 1; step <= config.run.max_steps; ++step) {
                const auto time_step = solver.global_time_step(config.run.cfl);
                solver.advance(time_step, time);
                time += time_step;
                if (mpi.rank() == 0) {
                    std::cout << "step=" << step << " time="
                              << std::setprecision(17) << time
                              << " dt=" << time_step
                              << " residual=" << solver.global_residual_l2()
                              << '\n';
                } else {
                    static_cast<void>(solver.global_residual_l2());
                }
            }
        }
        if (mpi.rank() == 0) {
            std::cout << "WCNS run completed at time=" << std::setprecision(17)
                      << time << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "wcns_run: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
