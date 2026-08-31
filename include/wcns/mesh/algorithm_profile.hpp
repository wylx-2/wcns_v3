#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

namespace wcns {

class ProfileError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class AlgorithmProfileKind {
    PhengleiWcns,
    Scmm6Wcns,
};

enum class MetricMethodKind {
    RefinedGridDelta,
    CommonCenterScmm6,
};

enum class LinearInterpolationKind {
    PhengleiI4,
    ScmmI6,
};

enum class FluxDerivativeKind {
    PhengleiD4D2,
    ScmmD6D4,
};

enum class PhysicalClosureKind {
    PhengleiOneSided,
    Scmm6OneSided,
};

struct ProfileComponents {
    AlgorithmProfileKind profile = AlgorithmProfileKind::PhengleiWcns;
    MetricMethodKind metric = MetricMethodKind::RefinedGridDelta;
    LinearInterpolationKind interpolation = LinearInterpolationKind::PhengleiI4;
    FluxDerivativeKind derivative = FluxDerivativeKind::PhengleiD4D2;
    PhysicalClosureKind physical_closure = PhysicalClosureKind::PhengleiOneSided;
};

class AlgorithmProfile {
public:
    [[nodiscard]] AlgorithmProfileKind kind() const noexcept
    {
        return components_.profile;
    }

    [[nodiscard]] const ProfileComponents& components() const noexcept
    {
        return components_;
    }

    [[nodiscard]] std::string name() const;
    [[nodiscard]] std::string restart_signature() const;

    void require_compatible(const ProfileComponents& components) const;

private:
    friend class ProfileFactory;
    explicit AlgorithmProfile(ProfileComponents components)
        : components_(components)
    {
    }

    ProfileComponents components_;
};

class ProfileFactory {
public:
    [[nodiscard]] static AlgorithmProfile create(AlgorithmProfileKind kind);
    [[nodiscard]] static AlgorithmProfile from_string(std::string_view name);
    static void validate_bundle(const ProfileComponents& components);
};

} // namespace wcns
