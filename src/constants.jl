module Constants

const earth_radius::Float64 = 6370.286  # km (computed at (37.411764°N, 92.394544°W))

# A polygon is composed of one or more lines, each of which is
# composed of 2D coordinates. Each line forms a complete ring. The
# polygon's first line is its external border, and other lines are
# internal borders.
const PolygonCoords = Vector{Vector{Vector{Float64}}}
const MultiPolygonCoords = Vector{PolygonCoords}

end  # module Constants
