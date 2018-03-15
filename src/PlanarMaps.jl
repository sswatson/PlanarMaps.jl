__precompile__(true)

module PlanarMaps

import AsyPlots,
       DataStructures

export NeighborCycle, CyclicRange, to,
       PlanarMap,
       Face,
       rotate, over,
       cw, ccw, ccwrange, insertccw, insertcw,
       isedge, add_edge!, delete_edge!,
       edges, face, faces, neighbors,
       interiorfaces, outerface,
       triangulation,
       RootedTree, root,
       isblue, isred, isgreen,
       vertices, interiorvertices,
       WoodedTriangulation, PartialWoodedTriangulation, shell!, shell,
       outerface, blueroot, redroot, greenroot,
       schnyderwood, allwoods,
       schnyder_blue, schnyder_red, schnyder_green, schnyder_coords,
       EmbeddedMap, draw, GraphPenSet,
       UWT

"""
    NeighborCycle(elements::Vector)

A circular permutation of distinct integers, representing
the (counterclockwise) neighbors of a vertex in a `PlanarMap`
"""
mutable struct NeighborCycle{T<:Integer}
    elements::Vector{T}
    lookup::Dict{T,T}
end

function remap!(N::NeighborCycle)
    N.lookup = Dict(map(reverse,enumerate(N.elements)));
end

function NeighborCycle(elements::Vector)
    NeighborCycle(elements,Dict(map(reverse,enumerate(elements))))
end

import Base.position
function position(N::NeighborCycle{T},u::T) where T<:Integer
    N.lookup[u]
end

"""
    insertccw(N::NeighborCycle{T},u::T,v::T)

Return the `NeighborCycle` obtained by inserting the vertex
`v` counterclockwise of the vertex `u` in `N`
"""
function insertccw(N::NeighborCycle{T},u::T,v::T) where T<:Integer
    nbs = copy(N.elements)
    k = position(N,u)
    insert!(nbs,k+1,v)
    NeighborCycle(nbs)
end

function insertcw(N::NeighborCycle{T},u::T,v::T) where T<:Integer
    nbs = copy(N.elements)
    k = position(N,u)
    insert!(nbs,k,v)
    NeighborCycle(nbs)
end

"""
    PlanarMap(nbs::Vector{NeighborCycle})

A planar map with vertex set `1:length(nbs)` and
for which the neighbors of vertex `i` listed in
counterclockwise order are `nbs[i]`
"""
struct PlanarMap{T<:Integer}
    nbs::Vector{NeighborCycle{T}}
end

function PlanarMap(a::Vector{<:Vector{<:Integer}})
    return PlanarMap(map(NeighborCycle,a))
end

vertices(P::PlanarMap) = 1:length(P)

function neighbors(P::PlanarMap)
    P.nbs
end

function neighbors(P::PlanarMap,i::Integer)
    P.nbs[i]
end

"""
    edges(P::PlanarMap)

Return a vector containing both `(a,b)` and `(b,a)`
for every edge `(a,b)` in `P`
"""
function edges(P::PlanarMap)
    vcat([[(u,v) for v in N.elements]
                for (u,N) in enumerate(neighbors(P))]...)
end

struct Face{T<:Integer}
    elements::Vector{T}
end

Face(args...) = Face([args...])

function Base.show(io::IO,F::Face)
    print(io,"Face(")
    show(io,F.elements)
    print(io,")")
end

function isedge(P::PlanarMap{T},
                i::T,
                j::T) where T<:Integer
    j in neighbors(P,i)
end

"""
Find a face which has (i,j) as a diagonal
"""
function containing_face(P::PlanarMap{T},
                         i::T,
                         j::T) where T<:Integer
    for v in neighbors(P,i)
        F = Face(P,i,v)
        k = findfirst(F,j)
        if k ≠ 0
            return F
        end
    end
    error("No face found which has a diagonal ($i,$j)")
end

"""
    add_edge!(P::PlanarMap,i::Integer,j::Integer,F::Face)
    add_edge!(P::PlanarMap,i::Integer,j::Integer)

Add edge from `i` to `j` through the face `F`
in the planar map `P`. If `F` is not given, the
face discovered first in the `NeighborCycle` of `i`
will be used
"""
function add_edge!(P::PlanarMap{T},
                   i::T,
                   j::T) where T<:Integer
    if isedge(P,i,j)
        error("($i,$j) already an edge of $P")
    end
    for v in neighbors(P,i)
        F = face(P,i,v)
        k = findfirst(F,j)
        if k ≠ 0
            P.nbs[i] = insertccw(P.nbs[i],F[2],j)
            P.nbs[j] = insertccw(P.nbs[j],F[k+1],i)
            return
        end
    end
end

function add_edge!(P::PlanarMap{T},
                   i::T,
                   j::T,
                   F::Face) where T<:Integer
    if isedge(P,i,j)
        error("($i,$j) already an edge of $P")
    end
    k = findfirst(F,i)
    l = findfirst(F,j)
    if k == 0 || l == 0
        error("$i and $j are not both in $F")
    end
    P.nbs[i] = insertccw(P.nbs[i],F[k+1],j)
    P.nbs[j] = insertccw(P.nbs[j],F[l+1],i);
end

import Base: filter, filter!
Base.filter!(f::Function,N::NeighborCycle) = (Base.filter!(f,N.elements); remap!(N))
Base.filter(f::Function,N::NeighborCycle) = NeighborCycle(filter(f,N.elements))
Base.filter!(f::Function,F::Face) = Base.filter!(f,F.elements)
Base.filter(f::Function,F::Face) = Base.filter(f,F.elements)

function delete_edge!(P::PlanarMap{T},
                      i::T,
                      j::T) where T<:Integer
    if ~(i in neighbors(P,j) || j in neighbors(P,i))
        error("No edge to delete")
    end
    filter!(x->x≠j,P.nbs[i])
    filter!(x->x≠i,P.nbs[j])
    foreach(remap!,(P.nbs[i],P.nbs[j]))
end
"""
    OrientedPlanarMap(nbs::Vector{NeighborCycle{T}},
                      outgoing::Vector{Set{T}}) where T<:Integer

An oriented planar map with vertex set `1:length(nbs)` and
for which the neighbors of vertex `i` listed in counterclockwise order
are `nbs[i]` and for which the outneighbors are `outgoing[i]`
"""
struct OrientedPlanarMap{T<:Integer}
    nbs::Vector{NeighborCycle{T}}
    outneighbors::Vector{Set{T}}
end

function neighbors(O::OrientedPlanarMap)
    return O.nbs
end

function outneighbors(O::OrientedPlanarMap)
    return O.outneighbors
end

function PlanarMap(O::OrientedPlanarMap)
    return PlanarMap(neighbors(O))
end

#-----DISPLAY-----------------------------------
import Base.show
function show(io::IO,N::NeighborCycle)
    print(io,"NeighborCycle(")
    print(io,string(N.elements))
    print(io,")")
end

function show(io::IO,M::Union{PlanarMap,OrientedPlanarMap})
    s = isa(M,PlanarMap) ? "" : "Oriented"
    print(io,"$(s)PlanarMap(")
    print(io,string(length(M)))
    print(io,")")
end
#-------------------------------------------------

import Base.length
length(N::NeighborCycle) = length(N.elements)
length(P::PlanarMap) = length(P.nbs)
length(O::OrientedPlanarMap) = length(O.nbs)

#--- CYCLE FUNCTIONS -----------------------------
"""
    over(N::NeighborCycle,u::T,k::T) where T<: Integer

Return the element which is `k` positions
to the right of `u` in the neighbor cycle `N`
"""
function over(N::NeighborCycle,u::T,k::T) where T<: Integer
    return N.elements[mod1(position(N,u)+k,length(N))]
end

"""
    rotate(N::NeighborCycle,k::Integer)

Return a vector containing the elements of `N`
in order beginning with position k
"""
function rotate(N::NeighborCycle{T},u::T) where T
    k = position(N,u)
    return [N.elements[k:end]; N.elements[1:k-1]]
end

prv(N::NeighborCycle,k::Integer) = over(N,k,-1)
nxt(N::NeighborCycle,k::Integer) = over(N,k,1)

import Base.==
function ==(N::NeighborCycle,M::NeighborCycle)
    if length(N) == 0 || length(M) == 0
        return length(N) == 0 && length(M) == 0
    end
    a = first(N.elements)
    N.elements == rotate(M,a)
end

function pairs(N::NeighborCycle)
    return zip(N.elements,rotate(N,2))
end

cw(P::PlanarMap{T},u::T,v::T) where T<:Integer =
                        over(neighbors(P,u),v,-1)
ccw(P::PlanarMap{T},u::T,v::T) where T<:Integer =
                        over(neighbors(P,u),v,1)

struct CyclicRange
    start::Int64
    done::Int64
end

to(a,b) = CyclicRange(a,b)

length(F::Face) = length(F.elements)
import Base: getindex, start, next, done, endof
Base.getindex(F::Face,i::Integer) = Base.getindex(F.elements,mod1(i,length(F)))
Base.getindex(F::Face,r::UnitRange) = Base.getindex(F.elements,r)
function Base.getindex(F::Face,C::CyclicRange)
    a,b = mod1.((C.start,C.done),length(F))
    if a ≤ b
        F[a:b]
    else
        vcat(F[a:end], F[1:b])
    end
end
function Base.getindex(N::NeighborCycle,C::CyclicRange)
    a,b = map(mod1,(C.start, C.done),(length(N),length(N)))
    if a ≤ b
        N[a:b]
    else
        vcat(N[a:end], N[1:b])
    end
end
Base.getindex(N::PlanarMaps.NeighborCycle{Int64},
              r::UnitRange{Int64}) = Base.getindex(N.elements,r)
Base.start(F::Face) = Base.start(F.elements)
Base.next(F::Face,i::Integer) = Base.next(F.elements,i)
Base.done(F::Face,i::Integer) = Base.done(F.elements,i)
Base.endof(F::Face) = Base.endof(F.elements)
Base.getindex(N::NeighborCycle,i::Integer) = Base.getindex(N.elements,mod1(i,length(N)))
Base.start(N::NeighborCycle) = Base.start(N.elements)
Base.next(N::NeighborCycle,i::Integer) = Base.next(N.elements,i)
Base.done(N::NeighborCycle,i::Integer) = Base.done(N.elements,i)
Base.endof(N::NeighborCycle) = Base.endof(N.elements)


"""
    face(P::PlanarMap,u::Int64,v::Int64,rotate=cw)

Return the face of the planar map P which includes the
edge (u,v) and finds each new edge by rotating in the
specified direction about the current vertex. Note
that this traces out the face in the direction *opposite*
to `rotate`
"""
function face(P::PlanarMap,u::Int64,v::Int64;rotate=cw)
    f = [u,v]
    while length(f) == 2 || (f[end-1],f[end]) ≠ (u,v)
        push!(f,rotate(P,f[end],f[end-1]))
    end
    return Face(f[1:end-2])
end

function edges(F::Face)
    ((F[i],F[mod1(i+1,length(F))]) for i=1:length(F))
end

function rotate(F::Face,k::Integer)
    return circshift(F.elements,-(k-1))
end

function ==(F::Face,G::Face)
    a = first(F)
    places = find(x->x==a,G)
    any(rotate(G,p) == F.elements for p in places)
end

struct FaceList
    interiorfaces::Vector{Face}
    outerface::Face
end

"""
    faces(P::PlanarMap)

Return a vector of faces of `P`
"""
function faces(RNG::AbstractRNG,
               P::PlanarMap;
               outeredge=(1,neighbors(P,1)[1]),
               rotate=cw)
    rem_edges = Set(edges(P))
    int_faces = Face[]
    outerface = face(P,pop!(rem_edges,outeredge)...;rotate=rotate)
    for e in edges(outerface)
        delete!(rem_edges,e)
    end
    while length(rem_edges) > 0
        e = rand(RNG,rem_edges)
        f = face(P,pop!(rem_edges,e)...;rotate=rotate)
        for e in edges(f)
	        delete!(rem_edges,e)
        end
        push!(int_faces,f)
    end
    return FaceList(int_faces,outerface)
end

faces(P::PlanarMap;kwargs...) = faces(MersenneTwister(0),P;kwargs...)

function allfaces(RNG::AbstractRNG,
                  P::PlanarMap;
                  outeredge=(1,neighbors(P,1)[1]),
                  rotate=cw)
    L = faces(RNG,P;outeredge=outeredge,rotate=rotate)
    return vcat(interiorfaces(L),[outerface(L)])
end

allfaces(P::PlanarMap;kwargs...) =
            allfaces(MersenneTwister(0),P;kwargs...)

function interiorfaces(L::FaceList)
    L.interiorfaces
end

function outerface(L::FaceList)
    L.outerface
end

function find_drawable_edge(P::PlanarMap,F::Face)
    if length(F) == 3
        error("No drawable edge; $F is a triangle")
    end
    for (i,v) in enumerate(F)
        j = findfirst(x-> x≠v && x∉neighbors(P,v), F)
        j > 0 && return (i,j)
    end
    error("No drawable edge found")
end

function poprand!(G::AbstractRNG,collection)
    splice!(collection,rand(G,1:length(collection)))
end

function facewithcorners(P::PlanarMap{T},blue::T,green::T,red::T) where T<:Integer
    for v in neighbors(P,blue)
        F = face(P,blue,v)
        k = findfirst(F.elements,green)
        # look for a red *after* the green already found:
        l = k+findfirst(F.elements[k+1:end],red)
        if 1 < k < l
            return F
        end
    end
    error("No face containing $blue, $green, $red")
end


poprand!(collection) = poprand(Base.Random.globalRNG(),collection)

"""
    triangulation(P::PlanarMap)

Find a simple plane triangulation obtained from P
by adding edges
"""
function triangulation(RNG::AbstractRNG,
                       P::PlanarMap;
                       outeredge=(1,neighbors(P,1)[1]),
                       corners=nothing)
    T = deepcopy(P)
    if corners ≠ nothing
        F = PlanarMaps.facewithcorners(T,corners...)
        for (a,b) in edges(Face(corners...))
            if !isedge(T,a,b)
                add_edge!(T,a,b,F)
            end
        end
        outeredge=corners[1:2]
    end
    facelist = [F for F in allfaces(T;outeredge=outeredge) if length(F) > 3]
    while !isempty(facelist)
        F = poprand!(RNG,facelist)
        (i,j) = find_drawable_edge(T,F)
        T.nbs[F[i]] = insertccw(T.nbs[F[i]],F[i+1],F[j])
        T.nbs[F[j]] = insertccw(T.nbs[F[j]],F[j+1],F[i])
        newfaces = [f for f in map(Face,(F[to(i,j)],F[to(j,i)])) if length(f) > 3]
        append!(facelist,newfaces)
    end
    return T
end

function triangulation(P::PlanarMap;kwargs...)
    triangulation(MersenneTwister(0),P;kwargs...)
end

import Base.in
Base.in(a,N::NeighborCycle) = Base.in(a,keys(N.lookup))

pairs(F::Face) = zip(F.elements,circshift(F.elements,-1))

include("woodedtriangulations.jl")
include("uniformwoodedtriangulations.jl")
include("graphing.jl")

end # module
