#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace wcns {

enum class SourceModelKind {
    UniformConservative,
    BodyForce,
    ManufacturedSolution,
};

[[nodiscard]] const char* source_model_name(SourceModelKind kind);

struct SourceTermConfig {
    bool enable_source_terms = false;
    std::vector<SourceModelKind> models;

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

// Stage H deliberately owns no concrete source model. The factory validates
// the complete configuration and rejects every enabled configuration until
// actual source evaluation is implemented in stage J.
class SourceTermRegistry {
public:
    [[nodiscard]] static SourceTermRegistry create_stage_h(
        const SourceTermConfig& config);

    [[nodiscard]] constexpr bool empty() const noexcept { return model_count_ == 0; }
    [[nodiscard]] constexpr std::size_t size() const noexcept { return model_count_; }

private:
    std::size_t model_count_ = 0;
};

} // namespace wcns
