#include <wcns/io/cgns_reader.hpp>

#include <cgnslib.h>

#include <array>
#include <cstdlib>
#include <cstddef>
#include <limits>
#include <string>
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
    const char* coordinate_name)
{
    std::array<cgsize_t, 3> lower {{1, 1, 1}};
    std::array<cgsize_t, 3> upper {{
        static_cast<cgsize_t>(zone.vertex_extent.ni),
        static_cast<cgsize_t>(zone.vertex_extent.nj),
        static_cast<cgsize_t>(zone.vertex_extent.nk),
    }};
    std::vector<double> values(zone.vertex_extent.size());
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
    int cell_dimension)
{
    const auto expected_point_values = static_cast<std::size_t>(2 * cell_dimension);
    if (points.size() != expected_point_values) {
        throw CgnsError("boundary PointRange has an unexpected number of index values");
    }

    IndexRange3 range;
    for (int axis = 0; axis < cell_dimension; ++axis) {
        range.begin[static_cast<std::size_t>(axis)] = checked_dimension(
            points[static_cast<std::size_t>(axis)], "boundary range begin")
            - 1;
        range.end[static_cast<std::size_t>(axis)] = checked_dimension(
            points[static_cast<std::size_t>(cell_dimension + axis)],
            "boundary range end")
            - 1;
    }
    if (cell_dimension == 2) {
        range.begin.k = 0;
        range.end.k = 0;
    }
    validate_range_within_extent(
        range, vertex_extent, cell_dimension, "boundary vertex PointRange");
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

IndexRange3 make_cell_face_range(
    const IndexRange3& vertex_range,
    FaceLocation face,
    Extent3 cell_extent,
    int cell_dimension)
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
    validate_range_within_extent(
        cells, cell_extent, cell_dimension, "boundary cell-face range");
    return cells;
}

void read_boundaries(int file, const CgnsZoneMetadata& zone, StructuredBlock& block)
{
    int boundary_count = 0;
    check_cgns(
        cg_nbocos(file, zone.base_file_index, zone.zone_file_index, &boundary_count),
        "cg_nbocos");
    block.boundaries.reserve(static_cast<std::size_t>(boundary_count));

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
            points, block.vertex_extent(), zone.cell_dimension);
        const auto face = identify_boundary_face(
            vertex_range, block.vertex_extent(), zone.cell_dimension);
        block.boundaries.push_back({
            name,
            convert_boundary_type(boundary_type),
            face,
            vertex_range,
            make_cell_face_range(
                vertex_range, face, block.cell_extent(), zone.cell_dimension),
            {},
        });
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
    FileHandle file(path);
    StructuredBlock block(
        zone.block_id,
        zone.name,
        owner_rank,
        zone.cell_dimension,
        zone.physical_dimension,
        zone.vertex_extent,
        ghost_width);

    if (!coordinate_exists(file.id(), zone, "CoordinateX")
        || !coordinate_exists(file.id(), zone, "CoordinateY")) {
        throw CgnsError("structured zone must provide CoordinateX and CoordinateY");
    }
    copy_coordinate(read_coordinate(file.id(), zone, "CoordinateX"), block.coordinates.x);
    copy_coordinate(read_coordinate(file.id(), zone, "CoordinateY"), block.coordinates.y);

    if (coordinate_exists(file.id(), zone, "CoordinateZ")) {
        copy_coordinate(read_coordinate(file.id(), zone, "CoordinateZ"), block.coordinates.z);
    } else if (zone.physical_dimension == 3) {
        throw CgnsError("three-dimensional physical space requires CoordinateZ");
    } else {
        block.coordinates.z.fill(0.0);
    }

    read_boundaries(file.id(), zone, block);

    return block;
}

} // namespace wcns
