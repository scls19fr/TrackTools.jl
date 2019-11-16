"""
A set of tools to deal with GPS track (find intersection with line, interpolate time of intersection...).
"""
module TrackTools

using Dates


const gamma = 1.0
const G = [gamma; gamma]

"""
    TrackPoint(time::T, x::Float64, y::Float64) where {T <: Dates.AbstractDateTime}

A point on a track.
"""
struct TrackPoint{T<:Dates.AbstractDateTime}
    time::T
    x::Float64
    y::Float64
end

"""
    matrix_line(x1, y1, x2, y2)

Compute A matrix parameters (ie parameters α and β) of a line between 2 points (x1, y1) and (x2, y2).

Equation of line (with `𝛾 = 1`)

``
α \\cdot x + β \\cdot y = 𝛾
``

with 2 points 

``
\\begin{cases}
α \\cdot x_1 + β \\cdot y_1 = 𝛾 \\
α \\cdot x_2 + β \\cdot y_1 = 𝛾
\\end{cases}
``

with matrix notation

``
P \\cdot A = G

A = P^-1 \\cdot G
``

where

``
G =
\\begin{bmatrix}
𝛾 \\
𝛾
\\end{bmatrix}
``

and

``
P =
\\begin{bmatrix}
x_1 & y_1 \\
x_2 & y_2 
\\end{bmatrix}
``

and

``
A =
\\begin{bmatrix}
α \\
β
\\end{bmatrix}
``
"""
function matrix_line(x1, y1, x2, y2)
    P = [x1 y1; x2 y2]
    return inv(P) * G
end

"""
    hasintersect(x1, y1, x2, y2, xI, yI) -> Bool

Return `true` if (`xI`, `yI`) is inside a rectangle defined by (`x1`, `y1`) and (`x2`, `y2`).
"""
function hasintersect(x1, y1, x2, y2, xI, yI)
    isinxrange(x1, x2, xI) = min(x1, x2) <= xI <= max(x1, x2)
    isinyrange(y1, y2, yI) = min(y1, y2) <= yI <= max(y1, y2)
    return isinxrange(x1, x2, xI) && isinyrange(y1, y2, yI)
end

"""
    interpolate_time(x1, t1::T, x2, t2::T, xI) -> T where {T<:Dates.AbstractDateTime}

Interpolate time of intersection (with a millisecond resolution) given time and position on an x-axis (may also be y)

x(t) = (x2 - x1) / (t2 - t1) * (t - t1) + x1

Δt = t2 - t1

t - t1 = (x - x1) / (x2 - x1) * Δt

t = (x - x1) / (x2 - x1) * Δt + t1
"""
function interpolate_time(x1, t1::T, x2, t2::T, xI) where {T<:Dates.AbstractDateTime}
    Δt = Dates.value(t2 - t1)
    return Millisecond(round(Int, Δt * (xI - x1) / (x2 - x1))) + t1
end


"""
    interpolate_time(x1, y1, t1, x2, y2, t2, xI, yI) -> Dates.AbstractDateTime

Interpolate time of intersection (with a millisecond resolution) given time and axis position where difference was the most significative.
"""
function interpolate_time(x1, y1, t1::T, x2, y2, t2::T, xI, yI) where {T<:Dates.AbstractDateTime}
    if abs(x2 - x1) >= abs(y2 - y1)
        return interpolate_time(x1, t1, x2, t2, xI)
    else
        return interpolate_time(y1, t1, y2, t2, yI)
    end
end

"""
    interpolate_position(x1, y1, t1, x2, y2, t2, t) -> (x, y)

Interpolate position (x, y) at time t given 2 positions (x1, y1) and (x2, y2) at time t1 and t2.

x(t) = (x2 - x1) / (t2 - t1) * (t - t1) + x1

"""
function interpolate_position(x1, y1, t1::T, x2, y2, t2::T, t::T) where {T<:Dates.AbstractDateTime}
    x = (x2 - x1) * (t - t1) / (t2 - t1) + x1
    y = (y2 - y1) * (t - t1) / (t2 - t1) + y1
    return x, y
end

"""
    calc_intersect_coordinates(x1, y1, x2, y2, AL) -> (xI, yI)

Calculate intersect coordinates given 2 positions coordinates (`x1`, `y1`) and (`x2`, `y2`)
and matrix line parameters `AL`.


``
\begin{cases}
α \\cdot x + β \\cdot y = 𝛾
α_L \\cdot x + β_L \\cdot y = 𝛾
\end{cases}
``

or with matrices

``
T \\cdot X = G
X = T^-1 \\cdot G
``

with

``
T =
\begin{bmatrix}
α & β \\
α_L & β_L 
\end{bmatrix}
``

and

``
X =
\begin{bmatrix}
xI \\
yI
\end{bmatrix}
``

"""
function calc_intersect_coordinates(x1, y1, x2, y2, AL)
    A = matrix_line(x1, y1, x2, y2)
    T = [transpose(A); transpose(AL)]
    X = inv(T) * G
    xI, yI = X
    return xI, yI
end

"""
    find_intersections(x1L, y1L, x2L, y2L, points) -> TrackPoint[]

Find all intersections between track and a line defined by 2 points (`x1L`, `y1L`) and (`x2L`, `y2L`).
Tracks points are defined using a generator of tuples `(t, x, y)`` named `points`
where t is time (`DateTime` or `ZonedDateTime`) and position (`x`, `y`)
"""
function find_intersections(x1L, y1L, x2L, y2L, points)
    AL = matrix_line(x1L, y1L, x2L, y2L)
    x1, y1, x2, y2 = zero([x1L, y1L, x2L, y2L])
    ((t0, _, _), _) = iterate(points)
    t1 = t2 = t0
    I = TrackPoint(t1, x1, y1)
    lst_I = TrackPoint[]
    for (i, (t2, x2, y2)) in enumerate(points)
        if i > 1
            if (x1 != x2) || (y1 != y2)
                (xI, yI) = calc_intersect_coordinates(x1, y1, x2, y2, AL)
                if hasintersect(x1, y1, x2, y2, xI, yI)
                    tI = interpolate_time(x1, y1, t1, x2, y2, t2, xI, yI)
                    I = TrackPoint(tI, x1, y1)
                    push!(lst_I, I)
                end
            end
        end
        (x1, y1, t1) = (x2, y2, t2)
    end
    return lst_I
end

end # module
