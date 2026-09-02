#include <wcns/io/cgns_reader.hpp>

#include <cgnslib.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {
namespace {

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw CgnsError(std::string(operation) + ": " + cg_get_error());
    }
}

class FileHandle {
public:
    explicit FileHandle(const std::string& path)
    {
        check_cgns(cg_open(path.c_str(), CG_MODE_READ, &id_), "cg_open");
    }

    ~FileHandle()
    {
        if (id_ >= 0) {
            cg_close(id_);
        }
    }

    FileHandle(const FileHandle&) = delete;
    FileHandle& operator=(const FileHandle&) = delete;

    [[nodiscard]] int id() const noexcept { return id_; }

private:
    int id_ = -1;
};

int checked_dimension(cgsize_t value, const char* label)
{
    if (value <= 0 || value > static_cast<cgsize_t>(std::numeric_limits<int>::max())) {
        throw CgnsError(std::string(label) + " is outside the supported integer range");
    }
    return static_cast<int>(value);
}

Extent3 make_vertex_extent(const cgsize_t* size, int cell_dimension)
{
    return {
        checked_dimension(size[0], "zone vertex extent I"),
        checked_dimension(size[1], "zone vertex extent J"),
        cell_dimension == 3 ? checked_dimension(size[2], "zone vertex extent K") : 1,
    };
}

Extent3 make_cell_extent(Extent3 vertex_extent, int cell_dimension)
{
    return {
        vertex_extent.ni - 1,
        vertex_extent.nj - 1,
        cell_dimension == 3 ? vertex_extent.nk - 1 : 1,
    };
}

void validate_base_dimensions(int cell_dimension, int physical_dimension)
{
    if ((cell_dimension != 2 && cell_dimension != 3)
        || (physical_dimension != 2 && physical_dimension != 3)
        || cell_dimension > physical_dimension) {
        throw CgnsError("only 2D/3D structured CFD bases are supported");
    }
}

std::vector<double> read_coordinate(
    int file,
    const CgnsZoneMetadata& zone,
    const char* coordinate_name,
    Index3 vertex_begin,
    Index3 vertex_end)
{
    std::array<cgsize_t, 3> lower {{
        static_cast<cgsize_t>(vertex_begin.i + 1),
        static_cast<cgsize_t>(vertex_begin.j + 1),
        static_cast<cgsize_t>(vertex_begin.k + 1),
    }};
    std::array<cgsize_t, 3> upper {{
        static_cast<cgsize_t>(vertex_end.i + 1),
        static_cast<cgsize_t>(vertex_end.j + 1),
        static_cast<cgsize_t>(vertex_end.k + 1),
    }};
    const Extent3 extent {
        vertex_end.i - vertex_begin.i + 1,
        vertex_end.j - vertex_begin.j + 1,
        vertex_end.k - vertex_begin.k + 1,
    };
    if (vertex_begin.i < 0 || vertex_begin.j < 0 || vertex_begin.k < 0
        || vertex_end.i < vertex_begin.i || vertex_end.j < vertex_begin.j
        || vertex_end.k < vertex_begin.k
        || vertex_end.i >= zone.vertex_extent.ni
        || vertex_end.j >= zone.vertex_extent.nj
        || vertex_end.k >= zone.vertex_extent.nk) {
        throw CgnsError("coordinate read range is outside the zone vertex extent");
    }
    std::vector<double> values(extent.size());
    check_cgns(
        cg_coord_read(
            file,
            zone.base_file_index,
            zone.zone_file_index,
            coordinate_name,
            RealDouble,
            lower.data(),
            upper.data(),
            values.data()),
        "cg_coord_read");
    return values;
}

using RotationMatrix = std::array<std::array<Real, 3>, 3>;

RotationMatrix multiply_rotation(
    const RotationMatrix& lhs,
    const RotationMatrix& rhs)
{
    RotationMatrix result {{}};
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            for (int inner = 0; inner < 3; ++inner) {
                result[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)]
                    += lhs[static_cast<std::size_t>(row)]
                          [static_cast<std::size_t>(inner)]
                     * rhs[static_cast<std::size_t>(inner)]
                           [static_cast<std::size_t>(column)];
            }
        }
    }
    return result;
}

PeriodicTransform make_periodic_transform(
    const std::array<float, 3>& center,
    const std::array<float, 3>& angle,
    const std::array<float, 3>& translation,
    int physical_dimension)
{
    if (physical_dimension != 2 && physical_dimension != 3) {
        throw CgnsError("periodic transform requires physical dimension 2 or 3");
    }
    if (physical_dimension == 2
        && (angle[0] != 0.0F || angle[1] != 0.0F)) {
        throw CgnsError(
            "two-dimensional CGNS periodic rotation is ambiguous; "
            "use translation or a three-dimensional embedding");
    }
    const Real cx = std::cos(static_cast<Real>(angle[0]));
    const Real sx = std::sin(static_cast<Real>(angle[0]));
    const Real cy = std::cos(static_cast<Real>(angle[1]));
    const Real sy = std::sin(static_cast<Real>(angle[1]));
    // A two-dimensional periodic rotation is represented by its z angle.
    const Real cz = std::cos(static_cast<Real>(angle[2]));
    const Real sz = std::sin(static_cast<Real>(angle[2]));
    const RotationMatrix rx {{{{1.0, 0.0, 0.0}},
                              {{0.0, cx, -sx}},
                              {{0.0, sx, cx}}}};
    const RotationMatrix ry {{{{cy, 0.0, sy}},
                              {{0.0, 1.0, 0.0}},
                              {{-sy, 0.0, cy}}}};
    const RotationMatrix rz {{{{cz, -sz, 0.0}},
                              {{sz, cz, 0.0}},
                              {{0.0, 0.0, 1.0}}}};
    PeriodicTransform result;
    // CGNS applies rotations current->donor in x, then y, then z order.
    result.rotation = multiply_rotation(rz, multiply_rotation(ry, rx));
    const std::array<Real, 3> origin {{
        static_cast<Real>(center[0]),
        static_cast<Real>(center[1]),
        physical_dimension == 3 ? static_cast<Real>(center[2]) : 0.0,
    }};
    const auto rotated_origin = result.apply_vector(origin);
    for (int component = 0; component < 3; ++component) {
        const auto index = static_cast<std::size_t>(component);
        result.translation[index] = origin[index] - rotated_origin[index]
            + (component < physical_dimension
                ? static_cast<Real>(translation[index]) : 0.0);
    }
    if (!result.valid(physical_dimension)) {
        throw CgnsError("CGNS periodic transform is not a proper rigid transform");
    }
    return result;
}

PeriodicTransform read_periodic_transform(
    int file,
    const CgnsZoneMetadata& zone,
    int connection_index)
{
    std::array<float, 3> center {{}};
    std::array<float, 3> angle {{}};
    std::array<float, 3> translation {{}};
    const int status = cg_1to1_periodic_read(
        file,
        zone.base_file_index,
        zone.zone_file_index,
        connection_index,
        center.data(),
        angle.data(),
        translation.data());
    if (status == CG_NODE_NOT_FOUND) return {};
    check_cgns(status, "cg_1to1_periodic_read");
    return make_periodic_transform(
        center, angle, translation, zone.physical_dimension);
}

std::vector<double> read_coordinate(
    int file,
    const CgnsZoneMetadata& zone,
    const char* coordinate_name)
{
    return read_coordinate(
        file,
        zone,
        coordinate_name,
        {0, 0, 0},
        {
            zone.vertex_extent.ni - 1,
            zone.vertex_extent.nj - 1,
            zone.vertex_extent.nk - 1,
        });
}

void copy_coordinate(const std::vector<double>& source, Array3D<Real>& destination)
{
    const auto extent = destination.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const auto source_index = static_cast<std::size_t>((k * extent.nj + j) * extent.ni + i);
                destination(i, j, k) = source[source_index];
            }
        }
    }
}

bool coordinate_exists(int file, const CgnsZoneMetadata& zone, const char* expected_name)
{
    int coordinate_count = 0;
    check_cgns(
        cg_ncoords(file, zone.base_file_index, zone.zone_file_index, &coordinate_count),
        "cg_ncoords");
    for (int coordinate = 1; coordinate <= coordinate_count; ++coordinate) {
        DataType_t data_type = DataTypeNull;
        char name[33] = {};
        check_cgns(
            cg_coord_info(
                file,
                zone.base_file_index,
                zone.zone_file_index,
                coordinate,
                &data_type,
                name),
            "cg_coord_info");
        if (expected_name == std::string(name)) {
            return true;
        }
    }
    return false;
}

BoundaryType convert_boundary_type(BCType_t type)
{
    switch (type) {
    case BCFarfield:
        return BoundaryType::Farfield;
    case BCInflow:
    case BCInflowSubsonic:
    case BCInflowSupersonic:
        return BoundaryType::Inflow;
    case BCOutflow:
    case BCOutflowSubsonic:
    case BCOutflowSupersonic:
        return BoundaryType::Outflow;
    case BCSymmetryPlane:
    case BCSymmetryPolar:
        return BoundaryType::Symmetry;
    case BCWall:
    case BCWallInviscid:
        return BoundaryType::SlipWall;
    case BCWallViscous:
        return BoundaryType::NoSlipAdiabaticWall;
    case BCWallViscousIsothermal:
        return BoundaryType::NoSlipIsothermalWall;
    default:
        return BoundaryType::Undefined;
    }
}

void validate_range_within_extent(
    const IndexRange3& range,
    Extent3 extent,
    int dimension,
    const char* label)
{
    constexpr std::array<const char*, 3> axis_names {{"I", "J", "K"}};
    for (int axis = 0; axis < dimension; ++axis) {
        const auto axis_index = static_cast<std::size_t>(axis);
        const int axis_extent = extent[axis_index];
        if (axis_extent <= 0) {
            throw CgnsError(
                std::string(label) + " has a non-positive " + axis_names[axis_index]
                + " extent");
        }

        const auto validate_endpoint = [&](int coordinate, const char* endpoint) {
            if (coordinate < 0 || coordinate >= axis_extent) {
                throw CgnsError(
                    std::string(label) + " " + endpoint + " " + axis_names[axis_index]
                    + " index " + std::to_string(coordinate)
                    + " is outside [0, " + std::to_string(axis_extent - 1) + "]");
            }
        };
        validate_endpoint(range.begin[axis_index], "begin");
        validate_endpoint(range.end[axis_index], "end");
    }
}

IndexRange3 convert_vertex_range(
    const std::vector<cgsize_t>& points,
    Extent3 vertex_extent,
    int cell_dimension,
    const char* label)
{
    const auto expected_point_values = static_cast<std::size_t>(2 * cell_dimension);
    if (points.size() != expected_point_values) {
        throw CgnsError(std::string(label) + " has an unexpected number of index values");
    }

    IndexRange3 range;
    for (int axis = 0; axis < cell_dimension; ++axis) {
        range.begin[static_cast<std::size_t>(axis)] = checked_dimension(
            points[static_cast<std::size_t>(axis)], label)
            - 1;
        range.end[static_cast<std::size_t>(axis)] = checked_dimension(
            points[static_cast<std::size_t>(cell_dimension + axis)],
            label)
            - 1;
    }
    if (cell_dimension == 2) {
        range.begin.k = 0;
        range.end.k = 0;
    }
    validate_range_within_extent(range, vertex_extent, cell_dimension, label);
    return range;
}

FaceLocation identify_boundary_face(
    const IndexRange3& vertex_range,
    Extent3 vertex_extent,
    int cell_dimension)
{
    int normal_axis = -1;
    Side side = Side::Lower;
    for (int axis = 0; axis < cell_dimension; ++axis) {
        const auto coordinate = vertex_range.begin[static_cast<std::size_t>(axis)];
        if (coordinate != vertex_range.end[static_cast<std::size_t>(axis)]) {
            continue;
        }
        if (coordinate == 0) {
            if (normal_axis >= 0) {
                throw CgnsError("boundary PointRange describes an edge or corner, not a face");
            }
            normal_axis = axis;
            side = Side::Lower;
        } else if (coordinate == vertex_extent[static_cast<std::size_t>(axis)] - 1) {
            if (normal_axis >= 0) {
                throw CgnsError("boundary PointRange describes an edge or corner, not a face");
            }
            normal_axis = axis;
            side = Side::Upper;
        }
    }
    if (normal_axis < 0) {
        throw CgnsError("boundary PointRange is not located on a block face");
    }
    return {static_cast<Axis>(normal_axis), side};
}

IndexRange3 make_adjacent_cell_range(
    const IndexRange3& vertex_range,
    FaceLocation face,
    Extent3 cell_extent,
    int cell_dimension,
    const char* label)
{
    IndexRange3 cells;
    const int normal_axis = static_cast<int>(face.axis);
    for (int axis = 0; axis < cell_dimension; ++axis) {
        if (axis == normal_axis) {
            const int coordinate = face.side == Side::Lower
                ? 0
                : cell_extent[static_cast<std::size_t>(axis)] - 1;
            cells.begin[static_cast<std::size_t>(axis)] = coordinate;
            cells.end[static_cast<std::size_t>(axis)] = coordinate;
            continue;
        }

        const int first = vertex_range.begin[static_cast<std::size_t>(axis)];
        const int last = vertex_range.end[static_cast<std::size_t>(axis)];
        if (first == last) {
            throw CgnsError("boundary PointRange must span cells in every tangential direction");
        }
        if (last > first) {
            cells.begin[static_cast<std::size_t>(axis)] = first;
            cells.end[static_cast<std::size_t>(axis)] = last - 1;
        } else {
            cells.begin[static_cast<std::size_t>(axis)] = first - 1;
            cells.end[static_cast<std::size_t>(axis)] = last;
        }
    }
    if (cell_dimension == 2) {
        cells.begin.k = 0;
        cells.end.k = 0;
    }
    validate_range_within_extent(cells, cell_extent, cell_dimension, label);
    return cells;
}

IndexRange3 make_boundary_face_range(
    const IndexRange3& vertex_range,
    FaceLocation face,
    Extent3 cell_extent,
    int cell_dimension,
    const char* label)
{
    auto faces = make_adjacent_cell_range(
        vertex_range,
        face,
        cell_extent,
        cell_dimension,
        label);
    const auto normal_axis = static_cast<std::size_t>(face.axis);
    const int normal_coordinate
        = face.side == Side::Lower ? 0 : cell_extent[normal_axis];
    faces.begin[normal_axis] = normal_coordinate;
    faces.end[normal_axis] = normal_coordinate;

    auto face_extent = cell_extent;
    ++face_extent[normal_axis];
    validate_range_within_extent(faces, face_extent, cell_dimension, label);
    return faces;
}

std::vector<BoundaryPatch> read_boundaries(
    int file,
    const CgnsZoneMetadata& zone)
{
    int boundary_count = 0;
    check_cgns(
        cg_nbocos(file, zone.base_file_index, zone.zone_file_index, &boundary_count),
        "cg_nbocos");
    std::vector<BoundaryPatch> boundaries;
    boundaries.reserve(static_cast<std::size_t>(boundary_count));

    for (int boundary_index = 1; boundary_index <= boundary_count; ++boundary_index) {
        char name[33] = {};
        BCType_t boundary_type = BCTypeNull;
        PointSetType_t point_set_type = PointSetTypeNull;
        cgsize_t point_count = 0;
        int normal_index[3] = {};
        cgsize_t normal_list_size = 0;
        DataType_t normal_data_type = DataTypeNull;
        int data_set_count = 0;
        check_cgns(
            cg_boco_info(
                file,
                zone.base_file_index,
                zone.zone_file_index,
                boundary_index,
                name,
                &boundary_type,
                &point_set_type,
                &point_count,
                normal_index,
                &normal_list_size,
                &normal_data_type,
                &data_set_count),
            "cg_boco_info");
        if (point_set_type != PointRange || point_count != 2) {
            throw CgnsError("only PointRange physical boundaries are supported");
        }

        GridLocation_t location = GridLocationNull;
        check_cgns(
            cg_boco_gridlocation_read(
                file,
                zone.base_file_index,
                zone.zone_file_index,
                boundary_index,
                &location),
            "cg_boco_gridlocation_read");
        if (location != Vertex) {
            throw CgnsError("only Vertex-located boundary PointRange is supported");
        }

        std::vector<cgsize_t> points(
            static_cast<std::size_t>(point_count * zone.cell_dimension));
        check_cgns(
            cg_boco_read(
                file,
                zone.base_file_index,
                zone.zone_file_index,
                boundary_index,
                points.data(),
                nullptr),
            "cg_boco_read");

        const auto vertex_range = convert_vertex_range(
            points,
            zone.vertex_extent,
            zone.cell_dimension,
            "boundary vertex PointRange");
        const auto face = identify_boundary_face(
            vertex_range, zone.vertex_extent, zone.cell_dimension);
        boundaries.push_back({
            name,
            convert_boundary_type(boundary_type),
            face,
            VertexRange(vertex_range),
            AdjacentCellRange(make_adjacent_cell_range(
                vertex_range,
                face,
                zone.cell_extent,
                zone.cell_dimension,
                "boundary adjacent-cell range")),
            BoundaryFaceRange(make_boundary_face_range(
                vertex_range,
                face,
                zone.cell_extent,
                zone.cell_dimension,
                "boundary face range")),
            {},
        });
    }
    return boundaries;
}

const CgnsZoneMetadata& find_donor_zone(
    const CgnsMeshMetadata& metadata,
    int base_file_index,
    const char* donor_name)
{
    const CgnsZoneMetadata* result = nullptr;
    for (const auto& candidate : metadata.zones) {
        if (candidate.base_file_index == base_file_index && candidate.name == donor_name) {
            if (result != nullptr) {
                throw CgnsError(
                    std::string("CGNS donor zone name is ambiguous: ") + donor_name);
            }
            result = &candidate;
        }
    }
    if (result == nullptr) {
        throw CgnsError(std::string("CGNS connectivity references unknown donor zone: ")
            + donor_name);
    }
    return *result;
}

std::vector<ConnectivityPatch> read_connectivities(
    int file,
    const CgnsZoneMetadata& zone,
    const CgnsMeshMetadata& metadata,
    int ghost_width)
{
    int connection_count = 0;
    check_cgns(
        cg_n1to1(file, zone.base_file_index, zone.zone_file_index, &connection_count),
        "cg_n1to1");
    std::vector<ConnectivityPatch> connectivities;
    connectivities.reserve(static_cast<std::size_t>(connection_count));

    for (int connection_index = 1; connection_index <= connection_count;
         ++connection_index) {
        char name[33] = {};
        char donor_name[33] = {};
        std::array<cgsize_t, 6> receiver_points {};
        std::array<cgsize_t, 6> donor_points {};
        std::array<int, 3> raw_transform {{1, 2, 3}};
        check_cgns(
            cg_1to1_read(
                file,
                zone.base_file_index,
                zone.zone_file_index,
                connection_index,
                name,
                donor_name,
                receiver_points.data(),
                donor_points.data(),
                raw_transform.data()),
            "cg_1to1_read");

        const auto& donor_zone = find_donor_zone(
            metadata, zone.base_file_index, donor_name);
        if (zone.cell_dimension != donor_zone.cell_dimension) {
            throw CgnsError(
                std::string("connectivity ") + name
                + " joins zones with different cell dimensions");
        }

        const auto value_count = static_cast<std::size_t>(2 * zone.cell_dimension);
        const std::vector<cgsize_t> receiver_values(
            receiver_points.begin(), receiver_points.begin() + value_count);
        const std::vector<cgsize_t> donor_values(
            donor_points.begin(), donor_points.begin() + value_count);
        const auto receiver_vertex_range = convert_vertex_range(
            receiver_values,
            zone.vertex_extent,
            zone.cell_dimension,
            "connectivity receiver vertex range");
        const auto donor_vertex_range = convert_vertex_range(
            donor_values,
            donor_zone.vertex_extent,
            zone.cell_dimension,
            "connectivity donor vertex range");
        const auto receiver_face = identify_boundary_face(
            receiver_vertex_range, zone.vertex_extent, zone.cell_dimension);
        const auto donor_face = identify_boundary_face(
            donor_vertex_range, donor_zone.vertex_extent, zone.cell_dimension);

        IndexTransform transform;
        for (int axis = 0; axis < zone.cell_dimension; ++axis) {
            transform.receiver_to_donor[static_cast<std::size_t>(axis)]
                = raw_transform[static_cast<std::size_t>(axis)];
        }
        if (!transform.valid(zone.cell_dimension)) {
            throw CgnsError(std::string("connectivity ") + name
                + " has an invalid CGNS Transform");
        }
        if (transform.map(
                receiver_vertex_range.end,
                receiver_vertex_range.begin,
                donor_vertex_range.begin,
                zone.cell_dimension)
            != donor_vertex_range.end) {
            throw CgnsError(std::string("connectivity ") + name
                + " Transform does not map receiver range onto donor range");
        }
        const auto periodic = read_periodic_transform(
            file, zone, connection_index);

        connectivities.push_back({
            name,
            zone.block_id,
            donor_zone.block_id,
            invalid_rank_id,
            receiver_face,
            donor_face,
            ReceiverVertexRange(receiver_vertex_range),
            DonorVertexRange(donor_vertex_range),
            ReceiverAdjacentCellRange(make_adjacent_cell_range(
                receiver_vertex_range,
                receiver_face,
                zone.cell_extent,
                zone.cell_dimension,
                "connectivity receiver adjacent-cell range")),
            DonorAdjacentCellRange(make_adjacent_cell_range(
                donor_vertex_range,
                donor_face,
                donor_zone.cell_extent,
                zone.cell_dimension,
                "connectivity donor adjacent-cell range")),
            SharedFaceRange(make_boundary_face_range(
                receiver_vertex_range,
                receiver_face,
                zone.cell_extent,
                zone.cell_dimension,
                "connectivity shared-face range")),
            transform,
            ghost_width,
            invalid_connection_id,
            periodic,
        });
    }
    return connectivities;
}

StructuredBlock read_block_data(
    int file,
    const CgnsZoneMetadata& zone,
    const CgnsMeshMetadata& metadata,
    RankId owner_rank,
    int ghost_width)
{
    StructuredBlock block(
        zone.block_id,
        zone.name,
        owner_rank,
        zone.cell_dimension,
        zone.physical_dimension,
        zone.vertex_extent,
        ghost_width);

    if (!coordinate_exists(file, zone, "CoordinateX")
        || !coordinate_exists(file, zone, "CoordinateY")) {
        throw CgnsError("structured zone must provide CoordinateX and CoordinateY");
    }
    copy_coordinate(read_coordinate(file, zone, "CoordinateX"), block.coordinates.x);
    copy_coordinate(read_coordinate(file, zone, "CoordinateY"), block.coordinates.y);

    if (coordinate_exists(file, zone, "CoordinateZ")) {
        copy_coordinate(read_coordinate(file, zone, "CoordinateZ"), block.coordinates.z);
    } else if (zone.physical_dimension == 3) {
        throw CgnsError("three-dimensional physical space requires CoordinateZ");
    } else {
        block.coordinates.z.fill(0.0);
    }

    block.boundaries = read_boundaries(file, zone);
    block.connectivities = read_connectivities(
        file, zone, metadata, ghost_width);
    return block;
}

const CgnsZoneMetadata& find_source_zone(
    const CgnsMeshMetadata& metadata,
    BlockId source_zone)
{
    for (const auto& zone : metadata.zones) {
        if (zone.block_id == source_zone) return zone;
    }
    throw CgnsError(
        "partition leaf references unknown source zone "
        + std::to_string(source_zone));
}

Extent3 leaf_cell_extent(
    const CgnsPartitionLeaf& leaf,
    int dimension)
{
    return {
        leaf.cell_end.i - leaf.cell_begin.i,
        leaf.cell_end.j - leaf.cell_begin.j,
        dimension == 3 ? leaf.cell_end.k - leaf.cell_begin.k : 1,
    };
}

void validate_leaf(
    const CgnsPartitionLeaf& leaf,
    const CgnsZoneMetadata& zone)
{
    if (leaf.block_id < 0 || leaf.owner_rank < 0
        || leaf.cell_begin.i < 0 || leaf.cell_begin.j < 0
        || leaf.cell_end.i <= leaf.cell_begin.i
        || leaf.cell_end.j <= leaf.cell_begin.j
        || leaf.cell_end.i > zone.cell_extent.ni
        || leaf.cell_end.j > zone.cell_extent.nj) {
        throw CgnsError("partition leaf cell range is invalid");
    }
    if (zone.cell_dimension == 3) {
        if (leaf.cell_begin.k < 0 || leaf.cell_end.k <= leaf.cell_begin.k
            || leaf.cell_end.k > zone.cell_extent.nk) {
            throw CgnsError("three-dimensional partition K range is invalid");
        }
    } else if (leaf.cell_begin.k != 0 || leaf.cell_end.k != 1) {
        throw CgnsError("two-dimensional partition K range must be [0,1)");
    }
}

StructuredBlock make_leaf_block(
    int file,
    const CgnsZoneMetadata& zone,
    const CgnsPartitionLeaf& leaf,
    int ghost_width,
    bool read_coordinates)
{
    const auto cells = leaf_cell_extent(leaf, zone.cell_dimension);
    const Extent3 vertices {
        cells.ni + 1,
        cells.nj + 1,
        zone.cell_dimension == 3 ? cells.nk + 1 : 1,
    };
    StructuredBlock block(
        leaf.block_id,
        zone.name + "__part_" + std::to_string(leaf.block_id),
        leaf.owner_rank,
        zone.cell_dimension,
        zone.physical_dimension,
        vertices,
        ghost_width);
    if (!read_coordinates) {
        const auto invalid = std::numeric_limits<Real>::quiet_NaN();
        block.coordinates.x.fill(invalid);
        block.coordinates.y.fill(invalid);
        block.coordinates.z.fill(invalid);
        return block;
    }
    if (!coordinate_exists(file, zone, "CoordinateX")
        || !coordinate_exists(file, zone, "CoordinateY")) {
        throw CgnsError("structured zone must provide CoordinateX and CoordinateY");
    }
    const Index3 vertex_end {
        leaf.cell_end.i,
        leaf.cell_end.j,
        zone.cell_dimension == 3 ? leaf.cell_end.k : 0,
    };
    copy_coordinate(
        read_coordinate(
            file, zone, "CoordinateX", leaf.cell_begin, vertex_end),
        block.coordinates.x);
    copy_coordinate(
        read_coordinate(
            file, zone, "CoordinateY", leaf.cell_begin, vertex_end),
        block.coordinates.y);
    if (coordinate_exists(file, zone, "CoordinateZ")) {
        copy_coordinate(
            read_coordinate(
                file, zone, "CoordinateZ", leaf.cell_begin, vertex_end),
            block.coordinates.z);
    } else if (zone.physical_dimension == 3) {
        throw CgnsError("three-dimensional physical space requires CoordinateZ");
    } else {
        block.coordinates.z.fill(0.0);
    }
    return block;
}

std::optional<IndexRange3> intersect_leaf_vertices(
    const IndexRange3& range,
    const CgnsPartitionLeaf& leaf,
    int dimension)
{
    IndexRange3 result;
    for (int axis = 0; axis < dimension; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int range_low = std::min(range.begin[a], range.end[a]);
        const int range_high = std::max(range.begin[a], range.end[a]);
        const int low = std::max(range_low, leaf.cell_begin[a]);
        const int high = std::min(range_high, leaf.cell_end[a]);
        if (low > high) return std::nullopt;
        if (range.end[a] >= range.begin[a]) {
            result.begin[a] = low;
            result.end[a] = high;
        } else {
            result.begin[a] = high;
            result.end[a] = low;
        }
    }
    if (dimension == 2) {
        result.begin.k = 0;
        result.end.k = 0;
    }
    return result;
}

IndexRange3 localize_range(
    IndexRange3 range,
    const CgnsPartitionLeaf& leaf,
    int dimension)
{
    for (int axis = 0; axis < dimension; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        range.begin[a] -= leaf.cell_begin[a];
        range.end[a] -= leaf.cell_begin[a];
    }
    if (dimension == 2) {
        range.begin.k = 0;
        range.end.k = 0;
    }
    return range;
}

bool is_face_range(
    const IndexRange3& range,
    FaceLocation face,
    int dimension)
{
    for (int axis = 0; axis < dimension; ++axis) {
        const auto count = range.counts()[static_cast<std::size_t>(axis)];
        if (axis == static_cast<int>(face.axis)) {
            if (count != 1) return false;
        } else if (count < 2) {
            return false;
        }
    }
    return true;
}

std::vector<BoundaryPatch> slice_boundaries(
    const std::vector<BoundaryPatch>& source,
    const CgnsPartitionLeaf& leaf,
    Extent3 local_cells,
    int dimension)
{
    std::vector<BoundaryPatch> result;
    for (const auto& patch : source) {
        const auto intersection = intersect_leaf_vertices(
            patch.vertex_range.untyped(), leaf, dimension);
        if (!intersection || !is_face_range(*intersection, patch.face, dimension)) {
            continue;
        }
        const auto local_vertices = localize_range(*intersection, leaf, dimension);
        result.push_back({
            patch.name,
            patch.type,
            patch.face,
            VertexRange(local_vertices),
            AdjacentCellRange(make_adjacent_cell_range(
                local_vertices,
                patch.face,
                local_cells,
                dimension,
                "partitioned boundary adjacent-cell range")),
            BoundaryFaceRange(make_boundary_face_range(
                local_vertices,
                patch.face,
                local_cells,
                dimension,
                "partitioned boundary face range")),
            patch.parameters,
        });
    }
    return result;
}

void append_original_connectivities(
    std::vector<StructuredBlock>& blocks,
    const std::unordered_map<BlockId, std::vector<ConnectivityPatch>>& source_connections,
    const std::vector<CgnsPartitionLeaf>& leaves)
{
    std::unordered_map<BlockId, std::size_t> block_index;
    for (std::size_t index = 0; index < blocks.size(); ++index) {
        block_index.emplace(blocks[index].id(), index);
    }
    for (const auto& receiver_leaf : leaves) {
        auto& receiver_block = blocks.at(block_index.at(receiver_leaf.block_id));
        const auto connection_iterator = source_connections.find(
            receiver_leaf.source_zone);
        if (connection_iterator == source_connections.end()) continue;
        for (const auto& source : connection_iterator->second) {
            const auto receiver_intersection = intersect_leaf_vertices(
                source.receiver_vertex_range.untyped(),
                receiver_leaf,
                receiver_block.cell_dimension());
            if (!receiver_intersection
                || !is_face_range(
                    *receiver_intersection,
                    source.receiver_face,
                    receiver_block.cell_dimension())) {
                continue;
            }
            const IndexRange3 donor_candidate {
                source.transform.map(
                    receiver_intersection->begin,
                    source.receiver_vertex_range.begin,
                    source.donor_vertex_range.begin,
                    receiver_block.cell_dimension()),
                source.transform.map(
                    receiver_intersection->end,
                    source.receiver_vertex_range.begin,
                    source.donor_vertex_range.begin,
                    receiver_block.cell_dimension()),
            };
            for (const auto& donor_leaf : leaves) {
                if (donor_leaf.source_zone != source.donor_block) continue;
                const auto donor_intersection = intersect_leaf_vertices(
                    donor_candidate,
                    donor_leaf,
                    receiver_block.cell_dimension());
                if (!donor_intersection
                    || !is_face_range(
                        *donor_intersection,
                        source.donor_face,
                        receiver_block.cell_dimension())) {
                    continue;
                }
                const auto inverse = source.transform.inverse(
                    receiver_block.cell_dimension());
                const IndexRange3 receiver_piece {
                    inverse.map(
                        donor_intersection->begin,
                        source.donor_vertex_range.begin,
                        source.receiver_vertex_range.begin,
                        receiver_block.cell_dimension()),
                    inverse.map(
                        donor_intersection->end,
                        source.donor_vertex_range.begin,
                        source.receiver_vertex_range.begin,
                        receiver_block.cell_dimension()),
                };
                const auto local_receiver = localize_range(
                    receiver_piece,
                    receiver_leaf,
                    receiver_block.cell_dimension());
                const auto local_donor = localize_range(
                    *donor_intersection,
                    donor_leaf,
                    receiver_block.cell_dimension());
                const auto& donor_block = blocks.at(
                    block_index.at(donor_leaf.block_id));
                receiver_block.connectivities.push_back({
                    source.name + "__part_" + std::to_string(receiver_leaf.block_id)
                        + "_to_" + std::to_string(donor_leaf.block_id),
                    receiver_leaf.block_id,
                    donor_leaf.block_id,
                    donor_leaf.owner_rank,
                    source.receiver_face,
                    source.donor_face,
                    ReceiverVertexRange(local_receiver),
                    DonorVertexRange(local_donor),
                    ReceiverAdjacentCellRange(make_adjacent_cell_range(
                        local_receiver,
                        source.receiver_face,
                        receiver_block.cell_extent(),
                        receiver_block.cell_dimension(),
                        "partitioned receiver adjacent-cell range")),
                    DonorAdjacentCellRange(make_adjacent_cell_range(
                        local_donor,
                        source.donor_face,
                        donor_block.cell_extent(),
                        donor_block.cell_dimension(),
                        "partitioned donor adjacent-cell range")),
                    SharedFaceRange(make_boundary_face_range(
                        local_receiver,
                        source.receiver_face,
                        receiver_block.cell_extent(),
                        receiver_block.cell_dimension(),
                        "partitioned shared-face range")),
                    source.transform,
                    receiver_block.ghost_width(),
                    invalid_connection_id,
                    source.periodic,
                });
            }
        }
    }
}

void append_sibling_connectivities(
    std::vector<StructuredBlock>& blocks,
    const std::vector<CgnsPartitionLeaf>& leaves)
{
    std::unordered_map<BlockId, std::size_t> block_index;
    for (std::size_t index = 0; index < blocks.size(); ++index) {
        block_index.emplace(blocks[index].id(), index);
    }
    for (std::size_t first = 0; first < leaves.size(); ++first) {
        for (std::size_t second = first + 1; second < leaves.size(); ++second) {
            const auto& lhs = leaves[first];
            const auto& rhs = leaves[second];
            if (lhs.source_zone != rhs.source_zone) continue;
            const int dimension = blocks.at(block_index.at(lhs.block_id)).cell_dimension();
            int normal_axis = -1;
            Side lhs_side = Side::Upper;
            IndexRange3 global_vertices;
            bool valid = true;
            for (int axis = 0; axis < dimension; ++axis) {
                const auto a = static_cast<std::size_t>(axis);
                if (lhs.cell_end[a] == rhs.cell_begin[a]) {
                    if (normal_axis >= 0) {
                        valid = false;
                        break;
                    }
                    normal_axis = axis;
                    lhs_side = Side::Upper;
                    global_vertices.begin[a] = lhs.cell_end[a];
                    global_vertices.end[a] = lhs.cell_end[a];
                } else if (rhs.cell_end[a] == lhs.cell_begin[a]) {
                    if (normal_axis >= 0) {
                        valid = false;
                        break;
                    }
                    normal_axis = axis;
                    lhs_side = Side::Lower;
                    global_vertices.begin[a] = lhs.cell_begin[a];
                    global_vertices.end[a] = lhs.cell_begin[a];
                } else {
                    const int low = std::max(lhs.cell_begin[a], rhs.cell_begin[a]);
                    const int high = std::min(lhs.cell_end[a], rhs.cell_end[a]);
                    if (high <= low) {
                        valid = false;
                        break;
                    }
                    global_vertices.begin[a] = low;
                    global_vertices.end[a] = high;
                }
            }
            if (!valid || normal_axis < 0) continue;
            if (dimension == 2) {
                global_vertices.begin.k = 0;
                global_vertices.end.k = 0;
            }
            auto& lhs_block = blocks.at(block_index.at(lhs.block_id));
            auto& rhs_block = blocks.at(block_index.at(rhs.block_id));
            const FaceLocation lhs_face {
                static_cast<Axis>(normal_axis), lhs_side};
            const FaceLocation rhs_face {
                static_cast<Axis>(normal_axis),
                lhs_side == Side::Upper ? Side::Lower : Side::Upper};
            const auto lhs_vertices = localize_range(
                global_vertices, lhs, dimension);
            const auto rhs_vertices = localize_range(
                global_vertices, rhs, dimension);
            const auto add = [&](StructuredBlock& receiver,
                                 const CgnsPartitionLeaf& receiver_leaf,
                                 const IndexRange3& receiver_vertices,
                                 FaceLocation receiver_face,
                                 StructuredBlock& donor,
                                 const CgnsPartitionLeaf& donor_leaf,
                                 const IndexRange3& donor_vertices,
                                 FaceLocation donor_face) {
                receiver.connectivities.push_back({
                    "partition_internal_" + std::to_string(receiver_leaf.block_id)
                        + "_to_" + std::to_string(donor_leaf.block_id),
                    receiver_leaf.block_id,
                    donor_leaf.block_id,
                    donor_leaf.owner_rank,
                    receiver_face,
                    donor_face,
                    ReceiverVertexRange(receiver_vertices),
                    DonorVertexRange(donor_vertices),
                    ReceiverAdjacentCellRange(make_adjacent_cell_range(
                        receiver_vertices,
                        receiver_face,
                        receiver.cell_extent(),
                        dimension,
                        "partition sibling receiver cells")),
                    DonorAdjacentCellRange(make_adjacent_cell_range(
                        donor_vertices,
                        donor_face,
                        donor.cell_extent(),
                        dimension,
                        "partition sibling donor cells")),
                    SharedFaceRange(make_boundary_face_range(
                        receiver_vertices,
                        receiver_face,
                        receiver.cell_extent(),
                        dimension,
                        "partition sibling shared faces")),
                    {},
                    receiver.ghost_width(),
                });
            };
            add(
                lhs_block, lhs, lhs_vertices, lhs_face,
                rhs_block, rhs, rhs_vertices, rhs_face);
            add(
                rhs_block, rhs, rhs_vertices, rhs_face,
                lhs_block, lhs, lhs_vertices, lhs_face);
        }
    }
}

void copy_topology_to_local(
    const std::vector<StructuredBlock>& global_blocks,
    std::vector<StructuredBlock>& local_blocks)
{
    std::unordered_map<BlockId, const StructuredBlock*> sources;
    for (const auto& block : global_blocks) sources.emplace(block.id(), &block);
    for (auto& local : local_blocks) {
        const auto* source = sources.at(local.id());
        local.boundaries = source->boundaries;
        local.connectivities = source->connectivities;
    }
}

} // namespace

CgnsMeshMetadata CgnsReader::read_metadata(const std::string& path) const
{
    FileHandle file(path);
    CgnsMeshMetadata metadata;
    int base_count = 0;
    check_cgns(cg_nbases(file.id(), &base_count), "cg_nbases");
    BlockId next_block_id = 0;

    for (int base_index = 1; base_index <= base_count; ++base_index) {
        char base_name[33] = {};
        int cell_dimension = 0;
        int physical_dimension = 0;
        check_cgns(
            cg_base_read(
                file.id(), base_index, base_name, &cell_dimension, &physical_dimension),
            "cg_base_read");
        validate_base_dimensions(cell_dimension, physical_dimension);
        metadata.bases.push_back(
            {base_index, base_name, cell_dimension, physical_dimension});

        int zone_count = 0;
        check_cgns(cg_nzones(file.id(), base_index, &zone_count), "cg_nzones");
        for (int zone_index = 1; zone_index <= zone_count; ++zone_index) {
            ZoneType_t zone_type = ZoneTypeNull;
            check_cgns(
                cg_zone_type(file.id(), base_index, zone_index, &zone_type),
                "cg_zone_type");
            if (zone_type != Structured) {
                throw CgnsError("only Structured CGNS zones are supported");
            }

            char zone_name[33] = {};
            std::array<cgsize_t, 9> size {};
            check_cgns(
                cg_zone_read(file.id(), base_index, zone_index, zone_name, size.data()),
                "cg_zone_read");
            const auto vertex_extent = make_vertex_extent(size.data(), cell_dimension);
            metadata.zones.push_back({
                next_block_id++,
                base_index,
                zone_index,
                base_name,
                zone_name,
                cell_dimension,
                physical_dimension,
                vertex_extent,
                make_cell_extent(vertex_extent, cell_dimension),
            });
        }
    }
    return metadata;
}

StructuredBlock CgnsReader::read_block(
    const std::string& path,
    const CgnsZoneMetadata& zone,
    RankId owner_rank,
    int ghost_width) const
{
    const auto metadata = read_metadata(path);
    FileHandle file(path);
    return read_block_data(file.id(), zone, metadata, owner_rank, ghost_width);
}

StructuredMesh CgnsReader::read_mesh(
    const std::string& path,
    RankId owner_rank,
    int ghost_width) const
{
    const auto metadata = read_metadata(path);
    FileHandle file(path);
    std::vector<StructuredBlock> blocks;
    blocks.reserve(metadata.zones.size());
    for (const auto& zone : metadata.zones) {
        blocks.push_back(
            read_block_data(file.id(), zone, metadata, owner_rank, ghost_width));
    }
    for (auto& block : blocks) {
        for (auto& connection : block.connectivities) {
            connection.donor_rank = owner_rank;
        }
    }

    StructuredMesh mesh(std::move(blocks));
    try {
        mesh.validate_connectivities();
    } catch (const TopologyError& error) {
        throw CgnsError(std::string("invalid CGNS multiblock topology: ") + error.what());
    }
    return mesh;
}

CgnsPartitionedMesh CgnsReader::read_partitioned_mesh(
    const std::string& path,
    const std::vector<CgnsPartitionLeaf>& leaves,
    RankId local_rank,
    int ghost_width) const
{
    if (local_rank < 0 || ghost_width < 0 || leaves.empty()) {
        throw CgnsError("partitioned CGNS read arguments are invalid");
    }
    const auto metadata = read_metadata(path);
    FileHandle file(path);
    std::vector<StructuredBlock> global_blocks;
    std::vector<StructuredBlock> local_blocks;
    global_blocks.reserve(leaves.size());
    std::unordered_map<BlockId, std::vector<BoundaryPatch>> source_boundaries;
    std::unordered_map<BlockId, std::vector<ConnectivityPatch>> source_connections;
    std::unordered_map<BlockId, bool> block_ids;

    for (const auto& zone : metadata.zones) {
        source_boundaries.emplace(
            zone.block_id, read_boundaries(file.id(), zone));
        source_connections.emplace(
            zone.block_id,
            read_connectivities(file.id(), zone, metadata, ghost_width));
    }
    for (const auto& leaf : leaves) {
        const auto& zone = find_source_zone(metadata, leaf.source_zone);
        validate_leaf(leaf, zone);
        if (!block_ids.emplace(leaf.block_id, true).second) {
            throw CgnsError("partition leaf block ids must be unique");
        }
        auto global = make_leaf_block(
            file.id(), zone, leaf, ghost_width, false);
        global.boundaries = slice_boundaries(
            source_boundaries.at(zone.block_id),
            leaf,
            global.cell_extent(),
            zone.cell_dimension);
        global_blocks.push_back(std::move(global));
        if (leaf.owner_rank == local_rank) {
            local_blocks.push_back(make_leaf_block(
                file.id(), zone, leaf, ghost_width, true));
        }
    }

    append_original_connectivities(
        global_blocks, source_connections, leaves);
    append_sibling_connectivities(global_blocks, leaves);
    copy_topology_to_local(global_blocks, local_blocks);

    StructuredMesh global_mesh(std::move(global_blocks));
    try {
        global_mesh.validate_connectivities(false);
    } catch (const TopologyError& error) {
        throw CgnsError(
            std::string("invalid partitioned CGNS topology: ") + error.what());
    }
    return {std::move(global_mesh), std::move(local_blocks)};
}

} // namespace wcns
