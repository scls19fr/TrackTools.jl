using Test
using TrackTools
using TrackTools: interpolate_position
using TrackTools: matrix_line, calc_intersect_coordinates, hasintersect
using TrackTools: interpolate_time, find_intersections, TrackPoint
using TimeZones
using KMLTracks

@testset "interpolate_position" begin
    x1, y1, t1 = 0.0, 0.0, ZonedDateTime(2019, 11, 1, 10, 0, 0, tz"UTC")
    x2, y2, t2 = 2.0, 1.0, ZonedDateTime(2019, 11, 1, 10, 0, 40, tz"UTC")
    @test interpolate_position(x1, y1, t1, x2, y2, t2, ZonedDateTime(2019, 11, 1, 10, 0, 20, tz"UTC")) == (1.0, 0.5)
end

@testset "intersect sample" begin
    # positions
    x1, y1, t1 = 0.0, 0.2, ZonedDateTime(2019, 11, 1, 10, 0, 0, tz"UTC")
    x2, y2, t2 = 1.0, 0.2, ZonedDateTime(2019, 11, 1, 10, 0, 40, tz"UTC")

    # line
    x1L, y1L = 0.5, -0.2
    x2L, y2L = 0.5, 0.4

    # calculate matrix for line
    AL = matrix_line(x1L, y1L, x2L, y2L)

    # calculate coordinates of intersect
    xI, yI = calc_intersect_coordinates(x1, y1, x2, y2, AL)

    @test (xI, yI) == (0.5, 0.2)

    @test hasintersect(x1, y1, x2, y2, xI, yI)

    t = interpolate_time(x1, y1, t1, x2, y2, t2, xI, yI)
    @test t == ZonedDateTime(2019, 11, 1, 10, 0, 20, tz"UTC")

end

@testset "intersect sample" begin
    fname = "sample.kml"
    kmldoc = read_kml_file(fname)
    @test length(kmldoc.placemark.track.points) == 45
    #y1L, x1L = (46.585904, 0.309138)  # LBFI ALPHA2 SUD (lat/lon)
    #y2L, x2L = (46.585963, 0.309425)  # LFBI ALPHA2 NORD (lat/lon)

    y1L, x1L = (46.586296, 0.308861)  # LFBI ALPHA SUD
    y2L, x2L = (46.586480, 0.309045)  # LFBI ALPHA NORD

    points = ((p.time, p.longitude, p.latitude) for p in kmldoc.placemark.track.points)
    lst_I = find_intersections(x1L, y1L, x2L, y2L, points)

    @test length(lst_I) == 1
    @test lst_I[1].time == ZonedDateTime(2019, 11, 7, 8, 44, 49, 473, tz"UTC")

end
