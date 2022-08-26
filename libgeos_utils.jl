#
# Helper functions for LibGEOS
#

import LibGEOS

function make_multipolygon(
    multipolygon::Vector{Vector{Array{Float64, 2}}},
    )::LibGEOS.MultiPolygon
    return LibGEOS.MultiPolygon(
        LibGEOS.createCollection(
            LibGEOS.GEOS_MULTIPOLYGON,
            LibGEOS.GEOSGeom[
                LibGEOS.createPolygon(
                    create_linear_ring(coords[1]),
                    LibGEOS.GEOSGeom[create_linear_ring(c) for c in coords[2:end]],
                ) for coords in multipolygon
            ],
        ),
    )
end

function create_linear_ring(
    coords::Array{Float64, 2},
    )::Ptr{LibGEOS.GEOSGeometry}
    context = LibGEOS._context
    ncoords = size(coords, 2)
    coord_seq = LibGEOS.createCoordSeq(ncoords, context, ndim = 2)
    for j in 1:ncoords
        LibGEOS.GEOSCoordSeq_setX_r(
            context.ptr,
            coord_seq,
            j-1,
            coords[1,j],
        )
        LibGEOS.GEOSCoordSeq_setY_r(
            context.ptr,
            coord_seq,
            j-1,
            coords[2,j],
        )
    end
    return LibGEOS.GEOSGeom_createLinearRing_r(context.ptr, coord_seq)
end
