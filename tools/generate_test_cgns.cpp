#include <cgnslib.h>

#include <array>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(std::string(operation) + ": " + cg_get_error());
    }
}

class CgnsFile {
public:
    CgnsFile(const std::string& path, int mode)
    {
        check_cgns(cg_open(path.c_str(), mode, &id_), "cg_open");
    }

    ~CgnsFile()
    {
        if (id_ >= 0) {
            cg_close(id_);
        }
    }

    CgnsFile(const CgnsFile&) = delete;
    CgnsFile& operator=(const CgnsFile&) = delete;

    int id() const { return id_; }

    void close()
    {
        if (id_ >= 0) {
            check_cgns(cg_close(id_), "cg_close");
            id_ = -1;
        }
    }

private:
    int id_ = -1;
};

void write_boundary(
    int file,
    int base,
    int zone,
    const char* name,
    BCType_t type,
    const std::vector<cgsize_t>& point_range)
{
    int boundary = 0;
    check_cgns(
        cg_boco_write(
            file,
            base,
            zone,
            name,
            type,
            PointRange,
            2,
            point_range.data(),
            &boundary),
        "cg_boco_write");
    check_cgns(
        cg_boco_gridlocation_write(file, base, zone, boundary, Vertex),
        "cg_boco_gridlocation_write");
}

void write_one_to_one(
    int file,
    int base,
    int zone,
    const char* name,
    const char* donor_name,
    const std::vector<cgsize_t>& receiver_range,
    const std::vector<cgsize_t>& donor_range,
    const std::vector<int>& transform)
{
    int connection = 0;
    check_cgns(
        cg_1to1_write(
            file,
            base,
            zone,
            name,
            donor_name,
            receiver_range.data(),
            donor_range.data(),
            transform.data(),
            &connection),
        "cg_1to1_write");
}

void write_2d(const std::string& path, bool add_invalid_boundary = false)
{
    constexpr int ni = 5;
    constexpr int nj = 4;

    CgnsFile file(path, CG_MODE_WRITE);
    int base = 0;
    int zone = 0;
    int coordinate = 0;
    check_cgns(cg_base_write(file.id(), "Base2D", 2, 2, &base), "cg_base_write");

    cgsize_t size[6] = {ni, nj, ni - 1, nj - 1, 0, 0};
    check_cgns(
        cg_zone_write(file.id(), base, "Zone2D", size, Structured, &zone),
        "cg_zone_write");

    std::vector<double> x(static_cast<std::size_t>(ni * nj));
    std::vector<double> y(x.size());
    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            const auto index = static_cast<std::size_t>(j * ni + i);
            x[index] = 0.25 * static_cast<double>(i);
            y[index] = (1.0 / 3.0) * static_cast<double>(j);
        }
    }
    check_cgns(
        cg_coord_write(file.id(), base, zone, RealDouble, "CoordinateX", x.data(), &coordinate),
        "cg_coord_write CoordinateX");
    check_cgns(
        cg_coord_write(file.id(), base, zone, RealDouble, "CoordinateY", y.data(), &coordinate),
        "cg_coord_write CoordinateY");

    write_boundary(file.id(), base, zone, "imin", BCInflow, {1, 1, 1, nj});
    write_boundary(file.id(), base, zone, "imax", BCOutflow, {ni, 1, ni, nj});
    write_boundary(file.id(), base, zone, "jmin", BCWall, {1, 1, ni, 1});
    // Reverse the tangential direction to exercise directed PointRange handling.
    write_boundary(file.id(), base, zone, "jmax", BCWall, {ni, nj, 1, nj});
    if (add_invalid_boundary) {
        // I-min is a recognizable face, but its J endpoint exceeds the zone vertex extent.
        write_boundary(
            file.id(), base, zone, "invalid-j-extent", BCWall, {1, 1, 1, nj + 1});
    }
    file.close();
}

void write_3d(const std::string& path)
{
    constexpr int ni = 4;
    constexpr int nj = 3;
    constexpr int nk = 3;

    CgnsFile file(path, CG_MODE_WRITE);
    int base = 0;
    int zone = 0;
    int coordinate = 0;
    check_cgns(cg_base_write(file.id(), "Base3D", 3, 3, &base), "cg_base_write");

    cgsize_t size[9] = {ni, nj, nk, ni - 1, nj - 1, nk - 1, 0, 0, 0};
    check_cgns(
        cg_zone_write(file.id(), base, "Zone3D", size, Structured, &zone),
        "cg_zone_write");

    const auto point_count = static_cast<std::size_t>(ni * nj * nk);
    std::vector<double> x(point_count);
    std::vector<double> y(point_count);
    std::vector<double> z(point_count);
    for (int k = 0; k < nk; ++k) {
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                const auto index = static_cast<std::size_t>((k * nj + j) * ni + i);
                x[index] = static_cast<double>(i);
                y[index] = 0.5 * static_cast<double>(j);
                z[index] = 0.25 * static_cast<double>(k);
            }
        }
    }
    check_cgns(
        cg_coord_write(file.id(), base, zone, RealDouble, "CoordinateX", x.data(), &coordinate),
        "cg_coord_write CoordinateX");
    check_cgns(
        cg_coord_write(file.id(), base, zone, RealDouble, "CoordinateY", y.data(), &coordinate),
        "cg_coord_write CoordinateY");
    check_cgns(
        cg_coord_write(file.id(), base, zone, RealDouble, "CoordinateZ", z.data(), &coordinate),
        "cg_coord_write CoordinateZ");

    write_boundary(file.id(), base, zone, "imin", BCFarfield, {1, 1, 1, 1, nj, nk});
    write_boundary(file.id(), base, zone, "imax", BCFarfield, {ni, 1, 1, ni, nj, nk});
    write_boundary(file.id(), base, zone, "jmin", BCWall, {1, 1, 1, ni, 1, nk});
    write_boundary(file.id(), base, zone, "jmax", BCWall, {ni, nj, 1, 1, nj, nk});
    write_boundary(file.id(), base, zone, "kmin", BCWall, {1, 1, 1, ni, nj, 1});
    write_boundary(file.id(), base, zone, "kmax", BCWall, {1, 1, nk, ni, nj, nk});
    file.close();
}

void write_multiblock_2d(
    const std::string& path,
    bool reciprocal = true,
    const char* first_donor_name = "Right2D",
    bool invalid_receiver_range = false)
{
    constexpr int left_ni = 5;
    constexpr int left_nj = 4;
    constexpr int right_ni = 4;
    constexpr int right_nj = 5;

    CgnsFile file(path, CG_MODE_WRITE);
    int base = 0;
    int left_zone = 0;
    int right_zone = 0;
    int coordinate = 0;
    check_cgns(cg_base_write(file.id(), "BaseMulti2D", 2, 2, &base), "cg_base_write");

    cgsize_t left_size[6] = {
        left_ni, left_nj, left_ni - 1, left_nj - 1, 0, 0};
    cgsize_t right_size[6] = {
        right_ni, right_nj, right_ni - 1, right_nj - 1, 0, 0};
    check_cgns(
        cg_zone_write(file.id(), base, "Left2D", left_size, Structured, &left_zone),
        "cg_zone_write Left2D");
    check_cgns(
        cg_zone_write(file.id(), base, "Right2D", right_size, Structured, &right_zone),
        "cg_zone_write Right2D");

    std::vector<double> left_x(static_cast<std::size_t>(left_ni * left_nj));
    std::vector<double> left_y(left_x.size());
    for (int j = 0; j < left_nj; ++j) {
        for (int i = 0; i < left_ni; ++i) {
            const auto index = static_cast<std::size_t>(j * left_ni + i);
            left_x[index] = 0.25 * static_cast<double>(i);
            left_y[index] = (1.0 / 3.0) * static_cast<double>(j);
        }
    }
    check_cgns(
        cg_coord_write(
            file.id(), base, left_zone, RealDouble, "CoordinateX", left_x.data(), &coordinate),
        "cg_coord_write Left2D CoordinateX");
    check_cgns(
        cg_coord_write(
            file.id(), base, left_zone, RealDouble, "CoordinateY", left_y.data(), &coordinate),
        "cg_coord_write Left2D CoordinateY");

    std::vector<double> right_x(static_cast<std::size_t>(right_ni * right_nj));
    std::vector<double> right_y(right_x.size());
    for (int j = 0; j < right_nj; ++j) {
        for (int i = 0; i < right_ni; ++i) {
            const auto index = static_cast<std::size_t>(j * right_ni + i);
            right_x[index] = 1.0 + 0.25 * static_cast<double>(j);
            right_y[index] = 1.0 - (1.0 / 3.0) * static_cast<double>(i);
        }
    }
    check_cgns(
        cg_coord_write(
            file.id(), base, right_zone, RealDouble, "CoordinateX", right_x.data(), &coordinate),
        "cg_coord_write Right2D CoordinateX");
    check_cgns(
        cg_coord_write(
            file.id(), base, right_zone, RealDouble, "CoordinateY", right_y.data(), &coordinate),
        "cg_coord_write Right2D CoordinateY");

    write_one_to_one(
        file.id(),
        base,
        left_zone,
        "left-to-right",
        first_donor_name,
        {left_ni, 1, left_ni, left_nj + (invalid_receiver_range ? 1 : 0)},
        {right_ni, 1, 1, 1},
        {2, -1});
    if (reciprocal) {
        write_one_to_one(
            file.id(),
            base,
            right_zone,
            "right-to-left",
            "Left2D",
            {right_ni, 1, 1, 1},
            {left_ni, 1, left_ni, left_nj},
            {-2, 1});
    }
    file.close();
}

void write_multiblock_3d(const std::string& path)
{
    constexpr int ni = 4;
    constexpr int nj = 3;
    constexpr int nk = 3;

    CgnsFile file(path, CG_MODE_WRITE);
    int base = 0;
    int first_zone = 0;
    int second_zone = 0;
    int coordinate = 0;
    check_cgns(cg_base_write(file.id(), "BaseMulti3D", 3, 3, &base), "cg_base_write");
    cgsize_t size[9] = {ni, nj, nk, ni - 1, nj - 1, nk - 1, 0, 0, 0};
    check_cgns(
        cg_zone_write(file.id(), base, "First3D", size, Structured, &first_zone),
        "cg_zone_write First3D");
    check_cgns(
        cg_zone_write(file.id(), base, "Second3D", size, Structured, &second_zone),
        "cg_zone_write Second3D");

    const auto point_count = static_cast<std::size_t>(ni * nj * nk);
    std::vector<double> x(point_count);
    std::vector<double> y(point_count);
    std::vector<double> z(point_count);
    for (int k = 0; k < nk; ++k) {
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                const auto index = static_cast<std::size_t>((k * nj + j) * ni + i);
                x[index] = static_cast<double>(i) / static_cast<double>(ni - 1);
                y[index] = 0.5 * static_cast<double>(j);
                z[index] = 0.5 * static_cast<double>(k);
            }
        }
    }
    check_cgns(
        cg_coord_write(
            file.id(), base, first_zone, RealDouble, "CoordinateX", x.data(), &coordinate),
        "cg_coord_write First3D CoordinateX");
    check_cgns(
        cg_coord_write(
            file.id(), base, first_zone, RealDouble, "CoordinateY", y.data(), &coordinate),
        "cg_coord_write First3D CoordinateY");
    check_cgns(
        cg_coord_write(
            file.id(), base, first_zone, RealDouble, "CoordinateZ", z.data(), &coordinate),
        "cg_coord_write First3D CoordinateZ");

    for (int k = 0; k < nk; ++k) {
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                const auto index = static_cast<std::size_t>((k * nj + j) * ni + i);
                x[index] = 1.0 + static_cast<double>(i) / static_cast<double>(ni - 1);
                y[index] = 0.5 * static_cast<double>(nj - 1 - j);
                z[index] = 0.5 * static_cast<double>(k);
            }
        }
    }
    check_cgns(
        cg_coord_write(
            file.id(), base, second_zone, RealDouble, "CoordinateX", x.data(), &coordinate),
        "cg_coord_write Second3D CoordinateX");
    check_cgns(
        cg_coord_write(
            file.id(), base, second_zone, RealDouble, "CoordinateY", y.data(), &coordinate),
        "cg_coord_write Second3D CoordinateY");
    check_cgns(
        cg_coord_write(
            file.id(), base, second_zone, RealDouble, "CoordinateZ", z.data(), &coordinate),
        "cg_coord_write Second3D CoordinateZ");

    write_one_to_one(
        file.id(),
        base,
        first_zone,
        "first-to-second",
        "Second3D",
        {ni, 1, 1, ni, nj, nk},
        {1, nj, 1, 1, 1, nk},
        {1, -2, 3});
    write_one_to_one(
        file.id(),
        base,
        second_zone,
        "second-to-first",
        "First3D",
        {1, nj, 1, 1, 1, nk},
        {ni, 1, 1, ni, nj, nk},
        {1, -2, 3});
    file.close();
}

void verify(
    const std::string& path,
    int expected_dimension,
    int expected_zones,
    int expected_boundaries)
{
    CgnsFile file(path, CG_MODE_READ);
    int bases = 0;
    int zones = 0;
    int cell_dimension = 0;
    int physical_dimension = 0;
    int boundaries = 0;
    char base_name[33] = {};
    check_cgns(cg_nbases(file.id(), &bases), "cg_nbases");
    if (bases != 1) {
        throw std::runtime_error("generated CGNS file must contain one base");
    }
    check_cgns(
        cg_base_read(file.id(), 1, base_name, &cell_dimension, &physical_dimension),
        "cg_base_read");
    check_cgns(cg_nzones(file.id(), 1, &zones), "cg_nzones");
    check_cgns(cg_nbocos(file.id(), 1, 1, &boundaries), "cg_nbocos");
    if (cell_dimension != expected_dimension || physical_dimension != expected_dimension
        || zones != expected_zones || boundaries != expected_boundaries) {
        throw std::runtime_error("generated CGNS metadata does not match expectations");
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 9) {
        std::cerr
            << "usage: wcns_generate_test_cgns <2d-output.cgns> <3d-output.cgns> "
               "<invalid-2d-output.cgns> <multi-2d-output.cgns> "
               "<multi-3d-output.cgns> <one-sided-2d-output.cgns> "
               "<unknown-donor-2d-output.cgns> "
               "<invalid-connectivity-range-2d-output.cgns>\n";
        return EXIT_FAILURE;
    }

    try {
        check_cgns(cg_set_file_type(CG_FILE_ADF), "cg_set_file_type");
        write_2d(argv[1]);
        write_3d(argv[2]);
        write_2d(argv[3], true);
        write_multiblock_2d(argv[4]);
        write_multiblock_3d(argv[5]);
        write_multiblock_2d(argv[6], false);
        write_multiblock_2d(argv[7], false, "Missing2D");
        write_multiblock_2d(argv[8], false, "Right2D", true);
        verify(argv[1], 2, 1, 4);
        verify(argv[2], 3, 1, 6);
        verify(argv[3], 2, 1, 5);
        verify(argv[4], 2, 2, 0);
        verify(argv[5], 3, 2, 0);
        verify(argv[6], 2, 2, 0);
        verify(argv[7], 2, 2, 0);
        verify(argv[8], 2, 2, 0);
        std::cout << "generated and verified CGNS test grids\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
