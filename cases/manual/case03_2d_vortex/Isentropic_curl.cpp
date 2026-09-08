#include <cgnslib.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct PeriodicFace {
    const char *name;
    const char *donor_name;
    cgsize_t range[6];
    cgsize_t donor_range[6];
    int transform[3];
    float translation[3];
};

[[noreturn]] void fail(const std::string &message) {
    std::cerr << message << std::endl;
    std::exit(EXIT_FAILURE);
}

template <typename Fn>
void cg_check(Fn fn, const char *step) {
    if (fn != CG_OK) {
        std::cerr << "CGNS error while " << step << std::endl;
        cg_error_exit();
    }
}

cgsize_t node_index(cgsize_t i, cgsize_t j, cgsize_t k, cgsize_t ni, cgsize_t nj) {
    return (k * nj + j) * ni + i;
}

PeriodicFace make_periodic_face(const char *name,
                                const char *donor_name,
                                cgsize_t i0, cgsize_t j0, cgsize_t k0,
                                cgsize_t i1, cgsize_t j1, cgsize_t k1,
                                cgsize_t di0, cgsize_t dj0, cgsize_t dk0,
                                cgsize_t di1, cgsize_t dj1, cgsize_t dk1,
                                float tx, float ty, float tz) {
    PeriodicFace face{};
    face.name = name;
    face.donor_name = donor_name;
    face.range[0] = i0;
    face.range[1] = j0;
    face.range[2] = k0;
    face.range[3] = i1;
    face.range[4] = j1;
    face.range[5] = k1;
    face.donor_range[0] = di0;
    face.donor_range[1] = dj0;
    face.donor_range[2] = dk0;
    face.donor_range[3] = di1;
    face.donor_range[4] = dj1;
    face.donor_range[5] = dk1;
    face.transform[0] = 1;
    face.transform[1] = 2;
    face.transform[2] = 3;
    face.translation[0] = tx;
    face.translation[1] = ty;
    face.translation[2] = tz;
    return face;
}

}  // namespace

int main(int argc, char **argv) {
    std::string output_name = "Isentropic_curl.cgns";
    if (argc > 1) {
        output_name = argv[1];
    }

    const cgsize_t ni = 101;
    const cgsize_t nj = 101;
    const cgsize_t nk = 8;

    const cgsize_t n_nodes = ni * nj * nk;
    std::vector<double> x(n_nodes);
    std::vector<double> y(n_nodes);
    std::vector<double> z(n_nodes);

    const double pi = std::acos(-1.0);
    const double dx = 10.0 / 100.0;
    const double dy = 10.0 / 100.0;
    const double zmin = 0.0;
    const double zmax = 0.1;
    const double dz = (zmax - zmin) / static_cast<double>(nk - 1);

    for (cgsize_t k = 0; k < nk; ++k) {
        const double zk = zmin + static_cast<double>(k) * dz;
        for (cgsize_t j = 0; j < nj; ++j) {
            const double jj = static_cast<double>(j);
            const double x_shift = 0.5 * std::sin(2.0 * pi * jj / 100.0);
            for (cgsize_t i = 0; i < ni; ++i) {
                const double ii = static_cast<double>(i);
                const double y_shift = 0.5 * std::sin(4.0 * pi * ii / 100.0);
                const cgsize_t idx = node_index(i, j, k, ni, nj);
                x[idx] = ii * dx + x_shift;
                y[idx] = jj * dy + y_shift;
                z[idx] = zk;
            }
        }
    }

    cgsize_t size[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    size[0] = ni;
    size[1] = nj;
    size[2] = nk;
    size[3] = ni - 1;
    size[4] = nj - 1;
    size[5] = nk - 1;

    int file = 0;
    int base = 0;
    int zone = 0;
    int coord = 0;

    cg_check(cg_open(output_name.c_str(), CG_MODE_WRITE, &file), "opening output file");
    cg_check(cg_base_write(file, "Base", 3, 3, &base), "writing base");
    cg_check(cg_zone_write(file, base, "IsentropicCurlZone", size, CGNS_ENUMV(Structured), &zone),
             "writing structured zone");

    cg_check(cg_coord_write(file, base, zone, CGNS_ENUMV(RealDouble), "CoordinateX", x.data(), &coord),
             "writing CoordinateX");
    cg_check(cg_coord_write(file, base, zone, CGNS_ENUMV(RealDouble), "CoordinateY", y.data(), &coord),
             "writing CoordinateY");
    cg_check(cg_coord_write(file, base, zone, CGNS_ENUMV(RealDouble), "CoordinateZ", z.data(), &coord),
             "writing CoordinateZ");

    const PeriodicFace faces[] = {
        make_periodic_face("Imin_to_Imax", "IsentropicCurlZone",
                           1, 1, 1, 1, nj, nk,
                           ni, 1, 1, ni, nj, nk,
                           10.0f, 0.0f, 0.0f),
        make_periodic_face("Imax_to_Imin", "IsentropicCurlZone",
                           ni, 1, 1, ni, nj, nk,
                           1, 1, 1, 1, nj, nk,
                           -10.0f, 0.0f, 0.0f),
        make_periodic_face("Jmin_to_Jmax", "IsentropicCurlZone",
                           1, 1, 1, ni, 1, nk,
                           1, nj, 1, ni, nj, nk,
                           0.0f, 10.0f, 0.0f),
        make_periodic_face("Jmax_to_Jmin", "IsentropicCurlZone",
                           1, nj, 1, ni, nj, nk,
                           1, 1, 1, ni, 1, nk,
                           0.0f, -10.0f, 0.0f),
        make_periodic_face("Kmin_to_Kmax", "IsentropicCurlZone",
                           1, 1, 1, ni, nj, 1,
                           1, 1, nk, ni, nj, nk,
                           0.0f, 0.0f, 0.1f),
        make_periodic_face("Kmax_to_Kmin", "IsentropicCurlZone",
                           1, 1, nk, ni, nj, nk,
                           1, 1, 1, ni, nj, 1,
                           0.0f, 0.0f, -0.1f),
    };

    for (const PeriodicFace &face : faces) {
        int conn = 0;
        cg_check(cg_1to1_write(file, base, zone, face.name, face.donor_name,
                               face.range, face.donor_range, face.transform, &conn),
                 "writing 1-to-1 connectivity");
        float rotation_center[3] = {0.0f, 0.0f, 0.0f};
        float rotation_angle[3] = {0.0f, 0.0f, 0.0f};
        cg_check(cg_1to1_periodic_write(file, base, zone, conn,
                                        rotation_center, rotation_angle, face.translation),
                 "writing periodic connectivity property");
    }

    cg_check(cg_close(file), "closing file");

    std::cout << "Wrote " << output_name << " with an Isentropic curl grid "
              << ni << " x " << nj << " x " << nk << std::endl;
    return EXIT_SUCCESS;
}