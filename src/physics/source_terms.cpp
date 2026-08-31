#include <wcns/physics/source_terms.hpp>

#include <stdexcept>
#include <unordered_set>

namespace wcns {
namespace {

void validate_kind(SourceModelKind kind)
{
    switch (kind) {
    case SourceModelKind::UniformConservative:
    case SourceModelKind::BodyForce:
    case SourceModelKind::ManufacturedSolution:
        return;
    }
    throw std::invalid_argument("source-term configuration contains an unknown model");
}

} // namespace

const char* source_model_name(SourceModelKind kind)
{
    validate_kind(kind);
    switch (kind) {
    case SourceModelKind::UniformConservative:
        return "uniform_conservative";
    case SourceModelKind::BodyForce:
        return "body_force";
    case SourceModelKind::ManufacturedSolution:
        return "manufactured_solution";
    }
    throw std::logic_error("unreachable source model");
}

void SourceTermConfig::validate() const
{
    if (!enable_source_terms) {
        if (!models.empty()) {
            throw std::invalid_argument(
                "disabled source terms require an empty model list");
        }
        return;
    }
    if (models.empty()) {
        throw std::invalid_argument(
            "enabled source terms require at least one model");
    }

    std::unordered_set<int> seen;
    for (const auto model : models) {
        validate_kind(model);
        if (!seen.insert(static_cast<int>(model)).second) {
            throw std::invalid_argument(
                "source-term configuration contains a duplicate model");
        }
    }
}

std::string SourceTermConfig::summary() const
{
    validate();
    std::string result = "enable_source_terms=";
    result += enable_source_terms ? "true;models=" : "false;models=";
    for (std::size_t index = 0; index < models.size(); ++index) {
        if (index != 0) {
            result += ',';
        }
        result += source_model_name(models[index]);
    }
    return result;
}

std::string SourceTermConfig::restart_signature() const
{
    return "source_terms_v1;" + summary();
}

SourceTermRegistry SourceTermRegistry::create_stage_h(
    const SourceTermConfig& config)
{
    config.validate();
    if (config.enable_source_terms) {
        throw std::logic_error(
            "source-term models cannot be instantiated before stage J");
    }
    return {};
}

} // namespace wcns
