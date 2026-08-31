#include <wcns/mesh/algorithm_profile.hpp>

namespace wcns {
namespace {

ProfileComponents expected_components(AlgorithmProfileKind kind)
{
    switch (kind) {
    case AlgorithmProfileKind::PhengleiWcns:
        return {
            kind,
            MetricMethodKind::RefinedGridDelta,
            LinearInterpolationKind::PhengleiI4,
            FluxDerivativeKind::PhengleiD4D2,
            PhysicalClosureKind::PhengleiOneSided,
        };
    case AlgorithmProfileKind::Scmm6Wcns:
        return {
            kind,
            MetricMethodKind::CommonCenterScmm6,
            LinearInterpolationKind::ScmmI6,
            FluxDerivativeKind::ScmmD6D4,
            PhysicalClosureKind::Scmm6OneSided,
        };
    }
    throw ProfileError("unknown algorithm profile");
}

bool same_components(const ProfileComponents& lhs, const ProfileComponents& rhs)
{
    return lhs.profile == rhs.profile && lhs.metric == rhs.metric
        && lhs.interpolation == rhs.interpolation
        && lhs.derivative == rhs.derivative
        && lhs.physical_closure == rhs.physical_closure;
}

} // namespace

std::string AlgorithmProfile::name() const
{
    switch (kind()) {
    case AlgorithmProfileKind::PhengleiWcns:
        return "phenglei_wcns";
    case AlgorithmProfileKind::Scmm6Wcns:
        return "scmm6_wcns";
    }
    throw ProfileError("unknown algorithm profile");
}

std::string AlgorithmProfile::restart_signature() const
{
    return "algorithm_profile_v1;name=" + name();
}

void AlgorithmProfile::require_compatible(
    const ProfileComponents& components) const
{
    ProfileFactory::validate_bundle(components);
    if (!same_components(components_, components)) {
        throw ProfileError("algorithm component belongs to a different profile");
    }
}

AlgorithmProfile ProfileFactory::create(AlgorithmProfileKind kind)
{
    return AlgorithmProfile(expected_components(kind));
}

AlgorithmProfile ProfileFactory::from_string(std::string_view name)
{
    if (name == "phenglei_wcns") {
        return create(AlgorithmProfileKind::PhengleiWcns);
    }
    if (name == "scmm6_wcns") {
        return create(AlgorithmProfileKind::Scmm6Wcns);
    }
    throw ProfileError("algorithm_profile must be phenglei_wcns or scmm6_wcns");
}

void ProfileFactory::validate_bundle(const ProfileComponents& components)
{
    if (!same_components(components, expected_components(components.profile))) {
        throw ProfileError("metric, interpolation, derivative, and closure profiles cannot be mixed");
    }
}

} // namespace wcns
