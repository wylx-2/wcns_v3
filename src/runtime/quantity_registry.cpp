#include <wcns/runtime/quantity_registry.hpp>

#include <wcns/mesh/conservation_weights.hpp>
#include <wcns/solver/euler.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace wcns {
namespace {

using FieldEvaluator = std::function<Real(
    const StructuredBlock&,
    const MetricField&,
    Index3,
    const QuantityContext&)>;

class FunctionalFieldQuantity final : public IFieldQuantity {
public:
    FunctionalFieldQuantity(QuantityDescriptor descriptor, FieldEvaluator evaluator)
        : descriptor_(std::move(descriptor))
        , evaluator_(std::move(evaluator))
    {
        descriptor_.validate();
        if (!evaluator_) throw std::invalid_argument("field evaluator is empty");
    }

    const QuantityDescriptor& descriptor() const override { return descriptor_; }
    Real evaluate_cell(
        const StructuredBlock& block,
        const MetricField& metric,
        Index3 index,
        const QuantityContext& context) const override
    {
        return evaluator_(block, metric, index, context);
    }

private:
    QuantityDescriptor descriptor_;
    FieldEvaluator evaluator_;
};

using StatisticEvaluator = std::function<Real(const StatisticContext&)>;

class FunctionalStatisticQuantity final : public IStatisticQuantity {
public:
    FunctionalStatisticQuantity(
        QuantityDescriptor descriptor,
        StatisticEvaluator evaluator)
        : descriptor_(std::move(descriptor))
        , evaluator_(std::move(evaluator))
    {
        descriptor_.validate();
        if (!evaluator_) throw std::invalid_argument("statistic evaluator is empty");
    }

    const QuantityDescriptor& descriptor() const override { return descriptor_; }
    Real evaluate(const StatisticContext& context) const override
    {
        return evaluator_(context);
    }

private:
    QuantityDescriptor descriptor_;
    StatisticEvaluator evaluator_;
};

TemperaturePrimitiveState cell_temperature_state(
    const StructuredBlock& block,
    Index3 index,
    const QuantityContext& context)
{
    return temperature_primitive_from_conservative(
        load_conservative(block.flow.conservative, index),
        context.gas,
        context.reference,
        context.floors,
        block.cell_dimension());
}

PressurePrimitiveState cell_pressure_state(
    const StructuredBlock& block,
    Index3 index,
    const QuantityContext& context)
{
    return pressure_primitive(
        cell_temperature_state(block, index, context),
        context.gas,
        context.reference,
        context.floors,
        block.cell_dimension());
}

QuantityDescriptor make_descriptor(
    std::string name,
    std::string dimensional_unit,
    QuantityScale scale)
{
    QuantityDescriptor result;
    result.name = std::move(name);
    result.dimensional_unit = std::move(dimensional_unit);
    result.scale = scale;
    return result;
}

std::shared_ptr<const IFieldQuantity> field(
    QuantityDescriptor information,
    FieldEvaluator evaluator)
{
    return std::make_shared<FunctionalFieldQuantity>(
        std::move(information), std::move(evaluator));
}

std::shared_ptr<const IStatisticQuantity> statistic(
    QuantityDescriptor information,
    StatisticEvaluator evaluator)
{
    return std::make_shared<FunctionalStatisticQuantity>(
        std::move(information), std::move(evaluator));
}

Real weighted_conservative_integral(
    const StatisticContext& context,
    int component)
{
    std::unordered_map<BlockId, const PartitionLeaf*> leaves;
    for (const auto& leaf : context.partition.leaves()) {
        leaves.emplace(leaf.block, &leaf);
    }
    std::unordered_map<BlockId, std::array<LineConservationWeights, 3>> weights;
    for (const auto& zone : context.partition.zones()) {
        std::array<LineConservationWeights, 3> lines;
        lines[0] = build_line_conservation_weights(
            context.profile, zone.cell_extent.ni);
        lines[1] = build_line_conservation_weights(
            context.profile, zone.cell_extent.nj);
        lines[2] = build_line_conservation_weights(
            context.profile,
            zone.cell_dimension == 3 ? zone.cell_extent.nk : 1,
            zone.cell_dimension == 2);
        weights.emplace(zone.source_zone, std::move(lines));
    }
    Real local = 0.0;
    for (const auto& block : context.local_blocks.blocks()) {
        const auto leaf_iterator = leaves.find(block.id());
        const auto metric_iterator = context.metrics.find(block.id());
        if (leaf_iterator == leaves.end() || metric_iterator == context.metrics.end()) {
            throw std::invalid_argument("statistic is missing leaf or metric data");
        }
        const auto& leaf = *leaf_iterator->second;
        const auto& lines = weights.at(leaf.source_zone);
        const auto& jacobian = metric_iterator->second.jacobian();
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real weight
                        = lines[0].cell_weights[static_cast<std::size_t>(
                              leaf.cells.begin.i + i)]
                        * lines[1].cell_weights[static_cast<std::size_t>(
                              leaf.cells.begin.j + j)]
                        * lines[2].cell_weights[static_cast<std::size_t>(
                              leaf.cells.begin.k + k)]
                        * jacobian(i, j, k);
                    const Real value = block.flow.conservative(
                        i, j, k, component);
                    if (!std::isfinite(weight) || weight <= 0.0
                        || !std::isfinite(value)) {
                        throw PhysicsError(
                            "statistic encountered an invalid cell value");
                    }
                    local += weight * value;
                }
            }
        }
    }
    return context.mpi.sum(local);
}

} // namespace

void QuantityDescriptor::validate() const
{
    if (name.empty()) throw std::invalid_argument("quantity name is empty");
    if (location != TopologyLocation::Cell) {
        throw std::invalid_argument(
            "stage N quantity registry currently accepts cell quantities only");
    }
    if (scale == QuantityScale::LengthPower
        && length_power != -1 && length_power <= 0) {
        throw std::invalid_argument(
            "length-power quantity requires a positive power or dimension sentinel");
    }
    std::set<std::string> unique;
    for (const auto& dependency : dependencies) {
        if (dependency.empty() || !unique.insert(dependency).second) {
            throw std::invalid_argument("quantity dependency list is invalid");
        }
    }
}

Real quantity_scale_factor(
    const QuantityDescriptor& information,
    const QuantityContext& context,
    int dimension)
{
    if (!context.dimensional) return 1.0;
    switch (information.scale) {
    case QuantityScale::Dimensionless: return 1.0;
    case QuantityScale::Density: return context.reference.density();
    case QuantityScale::Velocity: return context.reference.velocity();
    case QuantityScale::Pressure:
    case QuantityScale::Energy:
        return context.reference.dynamic_pressure();
    case QuantityScale::SpecificEnergy:
        return context.reference.velocity() * context.reference.velocity();
    case QuantityScale::Temperature: return context.reference.temperature();
    case QuantityScale::Momentum:
        return context.reference.density() * context.reference.velocity();
    case QuantityScale::Viscosity: return context.reference.viscosity();
    case QuantityScale::LengthPower:
        return std::pow(
            context.reference.length(),
            information.length_power == -1 ? dimension : information.length_power);
    }
    throw std::invalid_argument("invalid quantity scale");
}

FieldQuantityRegistry FieldQuantityRegistry::create_builtin()
{
    FieldQuantityRegistry result;
    const auto conservative = [](int component) {
        return [component](const StructuredBlock& block, const MetricField&,
                           Index3 index, const QuantityContext&) {
            return block.flow.conservative(
                index.i, index.j, index.k, component);
        };
    };
    result.register_quantity(field(
        make_descriptor("rho", "kg/m^3", QuantityScale::Density),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            return cell_temperature_state(block, index, context)
                [temperature_density];
        }));
    const std::array<std::pair<const char*, int>, 3> velocities {{
        {"u", temperature_velocity_x},
        {"v", temperature_velocity_y},
        {"w", temperature_velocity_z},
    }};
    for (const auto& item : velocities) {
        result.register_quantity(field(
            make_descriptor(item.first, "m/s", QuantityScale::Velocity),
            [component = item.second](
                const StructuredBlock& block, const MetricField&, Index3 index,
                const QuantityContext& context) {
                return cell_temperature_state(block, index, context)
                    [static_cast<std::size_t>(component)];
            }));
    }
    result.register_quantity(field(
        make_descriptor("p", "Pa", QuantityScale::Pressure),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            return cell_pressure_state(block, index, context)[pressure];
        }));
    result.register_quantity(field(
        make_descriptor("T", "K", QuantityScale::Temperature),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            return cell_temperature_state(block, index, context)[temperature_value];
        }));
    result.register_quantity(field(
        make_descriptor("rho_u", "kg/(m^2*s)", QuantityScale::Momentum),
        conservative(momentum_x)));
    result.register_quantity(field(
        make_descriptor("rho_v", "kg/(m^2*s)", QuantityScale::Momentum),
        conservative(momentum_y)));
    result.register_quantity(field(
        make_descriptor("rho_w", "kg/(m^2*s)", QuantityScale::Momentum),
        conservative(momentum_z)));
    result.register_quantity(field(
        make_descriptor("rho_E", "Pa", QuantityScale::Energy),
        conservative(total_energy)));
    result.register_quantity(field(
        make_descriptor("sound_speed", "m/s", QuantityScale::Velocity),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            return thermodynamic_sound_speed(
                cell_temperature_state(block, index, context),
                context.gas, context.reference, context.floors,
                block.cell_dimension());
        }));
    result.register_quantity(field(
        make_descriptor("mach", "1", QuantityScale::Dimensionless),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            const auto state = cell_temperature_state(block, index, context);
            const Real speed = std::sqrt(
                state[temperature_velocity_x] * state[temperature_velocity_x]
                + state[temperature_velocity_y] * state[temperature_velocity_y]
                + state[temperature_velocity_z] * state[temperature_velocity_z]);
            return speed / thermodynamic_sound_speed(
                state, context.gas, context.reference, context.floors,
                block.cell_dimension());
        }));
    result.register_quantity(field(
        make_descriptor("total_enthalpy", "m^2/s^2", QuantityScale::SpecificEnergy),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            return thermodynamic_total_enthalpy(
                cell_temperature_state(block, index, context),
                context.gas, context.reference, context.floors,
                block.cell_dimension());
        }));
    result.register_quantity(field(
        make_descriptor("entropy_proxy", "1", QuantityScale::Dimensionless),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            const auto state = cell_pressure_state(block, index, context);
            return state[pressure]
                / std::pow(state[primitive_density], context.gas.gamma());
        }));
    result.register_quantity(field(
        make_descriptor("viscosity", "Pa*s", QuantityScale::Viscosity),
        [](const StructuredBlock& block, const MetricField&, Index3 index,
           const QuantityContext& context) {
            return context.transport.viscosity(
                cell_temperature_state(block, index, context)[temperature_value]);
        }));
    auto jacobian = make_descriptor(
        "jacobian", "m^d", QuantityScale::LengthPower);
    jacobian.length_power = -1;
    result.register_quantity(field(
        std::move(jacobian),
        [](const StructuredBlock&, const MetricField& metric, Index3 index,
           const QuantityContext&) {
            return metric.jacobian()(index.i, index.j, index.k);
        }));
    return result;
}

void FieldQuantityRegistry::register_quantity(
    std::shared_ptr<const IFieldQuantity> quantity)
{
    if (!quantity) throw std::invalid_argument("field quantity is null");
    quantity->descriptor().validate();
    const auto name = quantity->descriptor().name;
    if (!quantities_.emplace(name, std::move(quantity)).second) {
        throw std::invalid_argument("duplicate field quantity: " + name);
    }
}

bool FieldQuantityRegistry::contains(const std::string& name) const noexcept
{
    return quantities_.find(name) != quantities_.end();
}

const QuantityDescriptor& FieldQuantityRegistry::descriptor(
    const std::string& name) const
{
    const auto iterator = quantities_.find(name);
    if (iterator == quantities_.end()) {
        throw std::invalid_argument("unknown field quantity: " + name);
    }
    return iterator->second->descriptor();
}

void FieldQuantityRegistry::validate_dependencies(const std::string& root) const
{
    std::set<std::string> active;
    std::set<std::string> finished;
    std::function<void(const std::string&)> visit = [&](const std::string& name) {
        if (finished.find(name) != finished.end()) return;
        if (!active.insert(name).second) {
            throw std::invalid_argument("cyclic field quantity dependency at " + name);
        }
        const auto iterator = quantities_.find(name);
        if (iterator == quantities_.end()) {
            throw std::invalid_argument("unknown field quantity dependency: " + name);
        }
        for (const auto& dependency : iterator->second->descriptor().dependencies) {
            visit(dependency);
        }
        active.erase(name);
        finished.insert(name);
    };
    visit(root);
}

void FieldQuantityRegistry::validate_selection(
    const std::vector<std::string>& names) const
{
    std::set<std::string> unique;
    for (const auto& name : names) {
        if (!unique.insert(name).second) {
            throw std::invalid_argument("field quantity was selected twice: " + name);
        }
        validate_dependencies(name);
    }
}

QuantityField FieldQuantityRegistry::evaluate(
    const std::string& name,
    const StructuredBlock& block,
    const MetricField& metric,
    const QuantityContext& context) const
{
    validate_dependencies(name);
    const auto& quantity = *quantities_.at(name);
    QuantityField result;
    result.descriptor = quantity.descriptor();
    result.extent = block.cell_extent();
    result.values.resize(result.extent.size());
    const Real scale = quantity_scale_factor(
        result.descriptor, context, block.cell_dimension());
    std::size_t offset = 0;
    for (int k = 0; k < result.extent.nk; ++k) {
        for (int j = 0; j < result.extent.nj; ++j) {
            for (int i = 0; i < result.extent.ni; ++i) {
                const Real value = quantity.evaluate_cell(
                    block, metric, {i, j, k}, context) * scale;
                if (!std::isfinite(value)) {
                    throw PhysicsError(
                        "field quantity produced a non-finite value: " + name);
                }
                result.values[offset++] = value;
            }
        }
    }
    return result;
}

StatisticRegistry StatisticRegistry::create_builtin()
{
    StatisticRegistry result;
    const std::array<std::tuple<const char*, int, QuantityScale>, 5> items {{
        {"total_mass", density, QuantityScale::Density},
        {"total_momentum_x", momentum_x, QuantityScale::Momentum},
        {"total_momentum_y", momentum_y, QuantityScale::Momentum},
        {"total_momentum_z", momentum_z, QuantityScale::Momentum},
        {"total_energy", total_energy, QuantityScale::Energy},
    }};
    for (const auto& item : items) {
        auto information = make_descriptor(
            std::get<0>(item), "integral", std::get<2>(item));
        result.register_quantity(statistic(
            std::move(information),
            [component = std::get<1>(item)](const StatisticContext& context) {
                return weighted_conservative_integral(context, component);
            }));
    }
    return result;
}

void StatisticRegistry::register_quantity(
    std::shared_ptr<const IStatisticQuantity> quantity)
{
    if (!quantity) throw std::invalid_argument("statistic quantity is null");
    quantity->descriptor().validate();
    const auto name = quantity->descriptor().name;
    if (!quantities_.emplace(name, std::move(quantity)).second) {
        throw std::invalid_argument("duplicate statistic quantity: " + name);
    }
}

bool StatisticRegistry::contains(const std::string& name) const noexcept
{
    return quantities_.find(name) != quantities_.end();
}

Real StatisticRegistry::evaluate(
    const std::string& name,
    const StatisticContext& context) const
{
    const auto iterator = quantities_.find(name);
    if (iterator == quantities_.end()) {
        throw std::invalid_argument("unknown statistic quantity: " + name);
    }
    const Real value = iterator->second->evaluate(context);
    if (!std::isfinite(value)) {
        throw PhysicsError("statistic produced a non-finite value: " + name);
    }
    if (!context.quantities.dimensional) return value;
    if (context.partition.zones().empty()) {
        throw std::invalid_argument("statistic partition has no source zones");
    }
    const int dimension = context.partition.zones().front().cell_dimension;
    for (const auto& zone : context.partition.zones()) {
        if (zone.cell_dimension != dimension) {
            throw std::invalid_argument(
                "statistic source zones use mixed cell dimensions");
        }
    }
    const Real field_scale = quantity_scale_factor(
        iterator->second->descriptor(), context.quantities, dimension);
    return value * field_scale
        * std::pow(context.quantities.reference.length(), dimension);
}

void StatisticRegistry::validate_selection(
    const std::vector<std::string>& names) const
{
    std::set<std::string> unique;
    for (const auto& name : names) {
        if (!unique.insert(name).second) {
            throw std::invalid_argument("statistic was selected twice: " + name);
        }
        if (!contains(name)) {
            throw std::invalid_argument("unknown statistic quantity: " + name);
        }
    }
}

} // namespace wcns
