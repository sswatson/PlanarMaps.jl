using PlanarMaps
using Base.Test

# write your own tests here
N = NeighborCycle([3,7,9,4])

@assert length(N) == 4
@assert over(N,9,3) == 7
@assert rotate(N,9) == [9,4,3,7]
@assert NeighborCycle([1,4,2,7]) == NeighborCycle([2,7,1,4])
@assert NeighborCycle([1,4,2,7]) ≠ NeighborCycle([1,4,7,2])

@assert Face([2,1,5,4,2,3]) == Face([4, 2, 3, 2, 1, 5])

P = PlanarMap([[2,5],[3,4,1],[2],[5,2],[4,1]])
F = face(P,1,2)
@assert F == Face([1,2,4,5])

@assert cw(P,1,2) == 5
@assert cw(P,2,4) == 3
@assert ccw(P,4,5) == 2

W = UWT(100)

draw(W)
draw(W,fillfaces=true,linecolor="black") 
draw(P)


# Below is the wooded triangulation from the article
# Schnyder woods, SLE(16), and Liouville Quantum Gravity
# by Li, Sun, and Watson, with the vertices numbered
# in order counterclockwise around the blue tree,
# starting with the first interior vertex and concluding
# with the three exterior vertices (10 = blue, 11 = red,
# 12 = green)
P = PlanarMap([
    [2,12,10,3],
    [12,1,3,4,5],
    [10, 11, 8, 5, 4, 2, 1],
    [5,2,3],
    [6,12,2,4,3,8,9,7],
    [12,5,7],
    [12,6,5,9,11],
    [9,5,3,11],
    [7,5,8,11],
    [1,12,11,3],
    [12,7,9,8,3,10],
    [10,1,2,5,6,7,11]])

bluetree = RootedTree(
    Dict(1=>10,2=>1,3=>10,4=>3,5=>3,6=>5,7=>5,8=>3,9=>8),
    10)

redtree = RootedTree(
    Dict(1=>3,2=>5,3=>11,4=>5,5=>9,6=>7,7=>11,8=>11,9=>11),
    11)

greentree = RootedTree(
    Dict(1=>12,2=>12,3=>2,4=>2,5=>12,6=>12,7=>12,8=>5,9=>7),
    12)

WT = WoodedTriangulation(P,bluetree,redtree,greentree)
