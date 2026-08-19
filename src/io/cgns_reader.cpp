#include <wcns/io/cgns_reader.hpp>

#include <cgnslib.h>

#include <array>
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

    return block;
}

} // namespace wcns

