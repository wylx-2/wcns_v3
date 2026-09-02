#include <wcns/runtime/field_output.hpp>

#include <wcns/io/cgns_reader.hpp>

#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace wcns {
namespace {

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(
            std::string(operation) + ": " + cg_get_error());
    }
}

class CgnsOutputFile {
public:
    explicit CgnsOutputFile(const std::string& path)
    {
        check_cgns(cg_open(path.c_str(), CG_MODE_WRITE, &file_), "cg_open output");
    }

    ~CgnsOutputFile()
    {
        if (file_ >= 0) cg_close(file_);
    }

    int id() const noexcept { return file_; }

    void close()
    {
        if (file_ < 0) return;
        check_cgns(cg_close(file_), "cg_close output");
        file_ = -1;
    }

private:
    int file_ = -1;
};

std::string join_path(const std::string& directory, const std::string& name)
{
    if (directory.empty()) return name;
    const char last = directory.back();
    return directory + ((last == '/' || last == '\\') ? "" : "/") + name;
}

std::string safe_name(std::string name)
{
    for (char& character : name) {
        const bool safe = (character >= 'a' && character <= 'z')
            || (character >= 'A' && character <= 'Z')
            || (character >= '0' && character <= '9')
            || character == '-' || character == '_';
        if (!safe) character = '_';
    }
    return name.empty() ? "case" : name;
}

std::string time_tag(Real time)
{
    std::ostringstream result;
    result << std::scientific << std::setprecision(9) << time;
    auto tag = result.str();
    for (char& character : tag) {
        if (character == '.') character = 'p';
        else if (character == '+') character = 'P';
        else if (character == '-') character = 'M';
    }
    return tag;
}

std::string field_stem(const CaseConfig& config, const SimulationState& state)
{
    std::ostringstream result;
    result << safe_name(config.case_name) << ".field.step"
           << std::setw(8) << std::setfill('0') << state.step
           << ".time" << time_tag(state.time);
    return result.str();
}

bool file_exists(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    return static_cast<bool>(input);
}

void commit_file(
    const std::string& temporary,
    const std::string& target,
    bool allow_existing)
{
    if (file_exists(target)) {
        if (!allow_existing) {
            throw std::runtime_error("field output already exists: " + target);
        }
        if (std::remove(target.c_str()) != 0) {
            throw std::runtime_error("cannot replace field output: " + target);
        }
    }
    if (std::rename(temporary.c_str(), target.c_str()) != 0) {
        throw std::runtime_error("cannot commit field output: " + target);
    }
}

std::size_t flat_index(Extent3 extent, int i, int j, int k)
{
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(extent.nj)
        + static_cast<std::size_t>(j)) * static_cast<std::size_t>(extent.ni)
        + static_cast<std::size_t>(i);
}

const StructuredBlock& local_block(
    const std::unordered_map<BlockId, const StructuredBlock*>& blocks,
    BlockId id)
{
    const auto iterator = blocks.find(id);
    if (iterator == blocks.end()) {
        throw std::runtime_error(
            "field output cannot find owned block " + std::to_string(id));
    }
    return *iterator->second;
}

std::string cgns_field_name(const std::string& name)
{
    static const std::unordered_map<std::string, std::string> standard {
        {"rho", "Density"}, {"u", "VelocityX"}, {"v", "VelocityY"},
        {"w", "VelocityZ"}, {"p", "Pressure"}, {"T", "Temperature"},
        {"rho_u", "MomentumX"}, {"rho_v", "MomentumY"},
        {"rho_w", "MomentumZ"}, {"rho_E", "EnergyStagnationDensity"},
        {"mach", "Mach"}, {"jacobian", "Jacobian"},
    };
    const auto iterator = standard.find(name);
    return iterator == standard.end() ? safe_name(name) : iterator->second;
}

} // namespace

OriginalFieldSnapshot gather_original_zone_fields(
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    const StructuredPartitionPlan& partition,
    const FieldQuantityRegistry& registry,
    const std::vector<std::string>& quantities,
    const QuantityContext& context,
    RankId root)
{
    registry.validate_selection(quantities);
    std::unordered_map<BlockId, const StructuredBlock*> blocks;
    for (const auto& block : local_blocks.blocks()) {
        blocks.emplace(block.id(), &block);
    }
    std::vector<Real> local_payload;
    const Real coordinate_scale = context.dimensional
        ? context.reference.length() : 1.0;
    for (const auto& leaf : partition.leaves()) {
        if (leaf.owner != mpi.rank()) continue;
        const auto& block = local_block(blocks, leaf.block);
        const auto metric_iterator = metrics.find(block.id());
        if (metric_iterator == metrics.end()) {
            throw std::runtime_error("field output is missing block metrics");
        }
        const auto& metric = metric_iterator->second;
        const auto extent = block.cell_extent();
        const auto append_coordinates = [&](const Array3D<Real>& values) {
            for (int k = 0; k < extent.nk; ++k) {
                for (int j = 0; j < extent.nj; ++j) {
                    for (int i = 0; i < extent.ni; ++i) {
                        local_payload.push_back(
                            coordinate_scale * values(i, j, k));
                    }
                }
            }
        };
        append_coordinates(metric.cell_coordinates().x);
        append_coordinates(metric.cell_coordinates().y);
        append_coordinates(metric.cell_coordinates().z);
        for (const auto& name : quantities) {
            const auto values = registry.evaluate(
                name, block, metric, context);
            local_payload.insert(
                local_payload.end(), values.values.begin(), values.values.end());
        }
    }
    const auto gathered = mpi.gather_reals(local_payload, root);
    if (mpi.rank() != root) return {};

    OriginalFieldSnapshot result;
    result.zones.reserve(partition.zones().size());
    std::unordered_map<BlockId, std::size_t> zone_indices;
    std::vector<std::vector<unsigned char>> coverage;
    for (const auto& zone : partition.zones()) {
        OriginalZoneField output;
        output.source_zone = zone.source_zone;
        output.name = zone.name;
        output.cell_dimension = zone.cell_dimension;
        output.cell_extent = zone.cell_extent;
        const auto count = zone.cell_extent.size();
        const Real nan = std::numeric_limits<Real>::quiet_NaN();
        output.x.assign(count, nan);
        output.y.assign(count, nan);
        output.z.assign(count, nan);
        for (const auto& name : quantities) {
            output.quantities.emplace(name, std::vector<Real>(count, nan));
        }
        zone_indices.emplace(zone.source_zone, result.zones.size());
        result.zones.push_back(std::move(output));
        coverage.emplace_back(count, 0);
    }

    std::size_t offset = 0;
    for (int rank = 0; rank < mpi.size(); ++rank) {
        for (const auto& leaf : partition.leaves()) {
            if (leaf.owner != rank) continue;
            auto& zone = result.zones.at(zone_indices.at(leaf.source_zone));
            auto& zone_coverage = coverage.at(zone_indices.at(leaf.source_zone));
            const auto extent = leaf.cell_extent();
            const auto count = extent.size();
            const std::size_t leaf_payload = (3 + quantities.size()) * count;
            if (offset + leaf_payload > gathered.size()) {
                throw std::runtime_error("field gather payload is truncated");
            }
            for (int k = 0; k < extent.nk; ++k) {
                for (int j = 0; j < extent.nj; ++j) {
                    for (int i = 0; i < extent.ni; ++i) {
                        const auto local = flat_index(extent, i, j, k);
                        const auto global = flat_index(
                            zone.cell_extent,
                            leaf.cells.begin.i + i,
                            leaf.cells.begin.j + j,
                            leaf.cells.begin.k + k);
                        if (zone_coverage[global] != 0) {
                            throw std::runtime_error(
                                "field gather overlaps an original-zone cell");
                        }
                        zone_coverage[global] = 1;
                        zone.x[global] = gathered[offset + local];
                        zone.y[global] = gathered[offset + count + local];
                        zone.z[global] = gathered[offset + 2 * count + local];
                        for (std::size_t quantity = 0;
                             quantity < quantities.size(); ++quantity) {
                            zone.quantities.at(quantities[quantity])[global]
                                = gathered[offset + (3 + quantity) * count + local];
                        }
                    }
                }
            }
            offset += leaf_payload;
        }
    }
    if (offset != gathered.size()) {
        throw std::runtime_error("field gather payload has trailing values");
    }
    for (std::size_t zone = 0; zone < result.zones.size(); ++zone) {
        if (std::find(coverage[zone].begin(), coverage[zone].end(), 0)
            != coverage[zone].end()) {
            throw std::runtime_error("field gather leaves original-zone cells missing");
        }
    }
    return result;
}

ProductionFieldWriter::ProductionFieldWriter(
    const MpiRuntime& mpi,
    const CaseConfig& config,
    const StructuredPartitionPlan& partition,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    QuantityContext quantity_context,
    std::string mesh_path,
    FieldQuantityRegistry registry)
    : mpi_(mpi)
    , config_(config)
    , partition_(partition)
    , local_blocks_(local_blocks)
    , metrics_(metrics)
    , quantity_context_(std::move(quantity_context))
    , mesh_path_(std::move(mesh_path))
    , registry_(std::move(registry))
{
    registry_.validate_selection(config_.output.field.quantities);
}

std::vector<std::string> ProductionFieldWriter::write(
    const SimulationState& state) const
{
    const auto snapshot = gather_original_zone_fields(
        mpi_, local_blocks_, metrics_, partition_, registry_,
        config_.output.field.quantities, quantity_context_);
    if (mpi_.rank() != 0) return {};
    std::vector<std::string> result;
    const auto stem = field_stem(config_, state);
    if (config_.output.field.format == FieldOutputFormat::Cgns
        || config_.output.field.format == FieldOutputFormat::Both) {
        const auto path = join_path(config_.output.directory, stem + ".cgns");
        write_cgns(snapshot, state, path);
        result.push_back(path);
    }
    if (config_.output.field.format == FieldOutputFormat::Tecplot
        || config_.output.field.format == FieldOutputFormat::Both) {
        const auto path = join_path(config_.output.directory, stem + ".dat");
        write_tecplot(snapshot, state, path);
        result.push_back(path);
    }
    return result;
}

void ProductionFieldWriter::write_cgns(
    const OriginalFieldSnapshot& snapshot,
    const SimulationState& state,
    const std::string& path) const
{
    const auto temporary = path + ".tmp";
    CgnsReader reader;
    const auto metadata = reader.read_metadata(mesh_path_);
    std::unordered_map<BlockId, const OriginalZoneField*> snapshot_zones;
    for (const auto& zone : snapshot.zones) {
        snapshot_zones.emplace(zone.source_zone, &zone);
    }
    CgnsOutputFile file(temporary);
    std::unordered_map<int, int> output_bases;
    for (const auto& base : metadata.bases) {
        int output_base = 0;
        check_cgns(
            cg_base_write(
                file.id(), base.name.c_str(), base.cell_dimension,
                base.physical_dimension, &output_base),
            "cg_base_write output");
        output_bases.emplace(base.file_index, output_base);
        check_cgns(
            cg_goto(file.id(), output_base, "end"),
            "cg_goto output base");
        check_cgns(
            cg_dataclass_write(
                config_.output.dimensional
                ? Dimensional : NormalizedByDimensional),
            "cg_dataclass_write output");
        if (config_.output.dimensional) {
            check_cgns(
                cg_units_write(Kilogram, Meter, Second, Kelvin, Radian),
                "cg_units_write output");
        }
    }
    const Real coordinate_scale = quantity_context_.dimensional
        ? quantity_context_.reference.length() : 1.0;
    for (const auto& zone : metadata.zones) {
        const auto snapshot_iterator = snapshot_zones.find(zone.block_id);
        if (snapshot_iterator == snapshot_zones.end()) {
            throw std::runtime_error("CGNS output is missing an original zone");
        }
        const auto& fields = *snapshot_iterator->second;
        auto block = reader.read_block(mesh_path_, zone, 0, 0);
        std::array<cgsize_t, 9> size {{}};
        for (int axis = 0; axis < zone.cell_dimension; ++axis) {
            size[static_cast<std::size_t>(axis)]
                = static_cast<cgsize_t>(zone.vertex_extent[axis]);
            size[static_cast<std::size_t>(zone.cell_dimension + axis)]
                = static_cast<cgsize_t>(zone.cell_extent[axis]);
        }
        const int base = output_bases.at(zone.base_file_index);
        int output_zone = 0;
        check_cgns(
            cg_zone_write(
                file.id(), base, zone.name.c_str(), size.data(),
                Structured, &output_zone),
            "cg_zone_write output");
        const auto vertex_count = zone.vertex_extent.size();
        std::vector<Real> coordinate(vertex_count);
        const std::array<std::pair<const char*, const Array3D<Real>*>, 3> coordinates {{
            {"CoordinateX", &block.coordinates.x},
            {"CoordinateY", &block.coordinates.y},
            {"CoordinateZ", &block.coordinates.z},
        }};
        for (int axis = 0; axis < zone.physical_dimension; ++axis) {
            std::size_t offset = 0;
            for (int k = 0; k < zone.vertex_extent.nk; ++k) {
                for (int j = 0; j < zone.vertex_extent.nj; ++j) {
                    for (int i = 0; i < zone.vertex_extent.ni; ++i) {
                        coordinate[offset++] = coordinate_scale
                            * (*coordinates[static_cast<std::size_t>(axis)].second)(i, j, k);
                    }
                }
            }
            int coordinate_index = 0;
            check_cgns(
                cg_coord_write(
                    file.id(), base, output_zone, RealDouble,
                    coordinates[static_cast<std::size_t>(axis)].first,
                    coordinate.data(), &coordinate_index),
                "cg_coord_write output");
        }
        int solution = 0;
        check_cgns(
            cg_sol_write(
                file.id(), base, output_zone, "FlowSolution",
                CellCenter, &solution),
            "cg_sol_write output");
        for (const auto& name : config_.output.field.quantities) {
            int field_index = 0;
            check_cgns(
                cg_field_write(
                    file.id(), base, output_zone, solution, RealDouble,
                    cgns_field_name(name).c_str(),
                    fields.quantities.at(name).data(), &field_index),
                "cg_field_write output");
        }
        check_cgns(
            cg_goto(
                file.id(), base, "Zone_t", output_zone,
                "FlowSolution_t", solution, "end"),
            "cg_goto output solution");
        std::ostringstream time;
        time << std::setprecision(17) << state.time;
        check_cgns(
            cg_descriptor_write("WCNS_Time", time.str().c_str()),
            "cg_descriptor_write output time");
    }
    file.close();
    commit_file(temporary, path, config_.output.allow_existing);
}

void ProductionFieldWriter::write_tecplot(
    const OriginalFieldSnapshot& snapshot,
    const SimulationState& state,
    const std::string& path) const
{
    const auto temporary = path + ".tmp";
    std::ofstream output(temporary, std::ios::out | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot open Tecplot output: " + temporary);
    output << "TITLE=\"WCNS step " << state.step << " time "
           << std::setprecision(17) << state.time << "\"\n"
           << "VARIABLES=\"X\",\"Y\"";
    const bool three_dimensional = !snapshot.zones.empty()
        && snapshot.zones.front().cell_dimension == 3;
    if (three_dimensional) output << ",\"Z\"";
    for (const auto& name : config_.output.field.quantities) {
        output << ",\"" << name << "\"";
    }
    output << '\n';
    for (const auto& zone : snapshot.zones) {
        output << "ZONE T=\"" << zone.name << "\", I="
               << zone.cell_extent.ni << ", J=" << zone.cell_extent.nj;
        if (zone.cell_dimension == 3) output << ", K=" << zone.cell_extent.nk;
        output << ", DATAPACKING=POINT, SOLUTIONTIME="
               << std::setprecision(17) << state.time << '\n';
        const auto count = zone.cell_extent.size();
        for (std::size_t cell = 0; cell < count; ++cell) {
            output << std::setprecision(17) << zone.x[cell] << ' '
                   << zone.y[cell];
            if (zone.cell_dimension == 3) output << ' ' << zone.z[cell];
            for (const auto& name : config_.output.field.quantities) {
                output << ' ' << zone.quantities.at(name)[cell];
            }
            output << '\n';
        }
    }
    output.close();
    if (!output) throw std::runtime_error("failed to write Tecplot output");
    commit_file(temporary, path, config_.output.allow_existing);
}

} // namespace wcns
