#pragma once

#include <wcns/core/array3d.hpp>
#include <wcns/core/types.hpp>

namespace wcns {

struct NodeCoordinates {
    explicit NodeCoordinates(Extent3 vertex_extent)
        : x(vertex_extent)
        , y(vertex_extent)
        , z(vertex_extent)
    {
    }

    Array3D<Real> x;
    Array3D<Real> y;
    Array3D<Real> z;
};

struct CellMetrics {
    CellMetrics(Extent3 cell_extent, int ghost_width)
        : center_x(cell_extent, ghost_width)
        , center_y(cell_extent, ghost_width)
        , center_z(cell_extent, ghost_width)
        , volume(cell_extent, ghost_width)
        , jacobian(cell_extent, ghost_width)
    {
    }

    Array3D<Real> center_x;
    Array3D<Real> center_y;
    Array3D<Real> center_z;
    Array3D<Real> volume;
    Array3D<Real> jacobian;
};

struct FaceMetric {
    explicit FaceMetric(Extent3 face_extent)
        : normal_x(face_extent)
        , normal_y(face_extent)
        , normal_z(face_extent)
        , area(face_extent)
    {
    }

    Array3D<Real> normal_x;
    Array3D<Real> normal_y;
    Array3D<Real> normal_z;
    Array3D<Real> area;
};

struct FaceMetrics {
    explicit FaceMetrics(Extent3 cell_extent)
        : i_faces({cell_extent.ni + 1, cell_extent.nj, cell_extent.nk})
        , j_faces({cell_extent.ni, cell_extent.nj + 1, cell_extent.nk})
        , k_faces({cell_extent.ni, cell_extent.nj, cell_extent.nk + 1})
    {
    }

    FaceMetric i_faces;
    FaceMetric j_faces;
    FaceMetric k_faces;
};

} // namespace wcns

