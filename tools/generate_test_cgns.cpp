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

void write_2d(const std::string& path)
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
    write_boundary(file.id(), base, zone, "jmax", BCWall, {1, nj, ni, nj});
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
    write_boundary(file.id(), base, zone, "jmax", BCWall, {1, nj, 1, ni, nj, nk});
    write_boundary(file.id(), base, zone, "kmin", BCWall, {1, 1, 1, ni, nj, 1});
    write_boundary(file.id(), base, zone, "kmax", BCWall, {1, 1, nk, ni, nj, nk});
    file.close();
}

void verify(const std::string& path, int expected_dimension, int expected_boundaries)
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
        || zones != 1 || boundaries != expected_boundaries) {
        throw std::runtime_error("generated CGNS metadata does not match expectations");
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: wcns_generate_test_cgns <2d-output.cgns> <3d-output.cgns>\n";
        return EXIT_FAILURE;
    }

    try {
        check_cgns(cg_set_file_type(CG_FILE_ADF), "cg_set_file_type");
        write_2d(argv[1]);
        write_3d(argv[2]);
        verify(argv[1], 2, 4);
        verify(argv[2], 3, 6);
        std::cout << "generated and verified CGNS test grids\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

