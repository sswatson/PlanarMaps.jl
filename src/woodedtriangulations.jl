
struct RootedTree{T<:Integer}
    ancestors::Dict{T,T}
    root::T
end

root(T::RootedTree) = T.root
ancestor(T::RootedTree{S},i::S) where S<:Integer = T.ancestors[i]
vertices(T::RootedTree) = keys(T.ancestors)∪[root(T)]

function ancestormap(T::RootedTree{S}) where S<:Integer
    D = Dict{S,Set{S}}(k=>Set() for k in vertices(T))
    for (child,ancestor) in T.ancestors
        push!(D[ancestor],child)
    end
    return D
end

function treesum(T::RootedTree{S},f::Dict{S,U}) where {S<:Integer, U<:Real}
    D = Dict{S,U}(k=>zero(U) for k in vertices(T))
    A = ancestormap(T)
    R = root(T)
    function updatebranch!(k)
        D[k] += f[k]
        if k ≠ R
            D[k] += D[ancestor(T,k)]
        end
        foreach(updatebranch!,A[k])
    end
    updatebranch!(R)
    return D
end

function backsum(T::RootedTree{S},f::Dict{S,U}) where {S<:Integer, U<:Real}
    D = Dict{S,U}()
    A = ancestormap(T)
    function updatevertex!(k)
        if k in keys(D)
            return D[k]
        else
            D[k] = f[k] + reduce(+,0,updatevertex!(v) for v in A[k])
        end
        D[k]
    end
    updatevertex!(root(T))
    D
end

import Base: +, -
+(D::Dict{<:Any,<:Real},t::Real) = Dict(k=>v+t for (k,v) in D)
-(D::Dict{<:Any,<:Real},t::Real) = D + (-t)

struct WoodedTriangulation{T<:Integer}
    P::PlanarMap{T}
    bluetree::RootedTree{T}
    redtree::RootedTree{T}
    greentree::RootedTree{T}
end

PlanarMap(W::WoodedTriangulation) = W.P

length(WT::WoodedTriangulation) = length(WT.P)

function show(io::IO,WT::WoodedTriangulation)
    print(io,"WoodedTriangulation(")
    print(io,string(length(WT.P)))
    print(io,")")
end

outerface(WT::WoodedTriangulation) =
    map(root,(WT.bluetree,WT.redtree,WT.greentree))

interiorvertices(WT) = setdiff(1:length(WT.P),outerface(WT))

blueroot(WT::WoodedTriangulation) = root(WT.bluetree)
redroot(WT::WoodedTriangulation) = root(WT.redtree)
greenroot(WT::WoodedTriangulation) = root(WT.greentree)

isblue(WT::WoodedTriangulation{T},i::T,j::T) where T<:Integer =
    i ∉ outerface(WT) && ancestor(WT.bluetree,i) == j
isred(WT::WoodedTriangulation{T},i::T,j::T) where T<:Integer =
    i ∉ outerface(WT) && ancestor(WT.redtree,i) == j
isgreen(WT::WoodedTriangulation{T},i::T,j::T) where T<:Integer =
    i ∉ outerface(WT) && ancestor(WT.greentree,i) == j

struct ColoredEdge
    tail::Int64
    head::Int64
    color::Symbol
end

function edges(WT::WoodedTriangulation)
    interior_edges = ColoredEdge[]
    for (d,c) in zip((WT.bluetree,WT.redtree,WT.greentree),(:blue,:red,:green))
        append!(interior_edges,vcat([[ColoredEdge(k,d[k],c),
                                 ColoredEdge(d[k],k,c)] for k in keys(d)]...))
    end
    roots = [f(WT) for f in (blueroot,redroot,greenroot)]
    vcat(interior_edges,
    [ColoredEdge(a,b,:gray) for a in roots for b in roots if a ≠ b])
end

function faces(WT::WoodedTriangulation{T}) where T<:Integer
    n = length(WT.P)-3
    all_edges = edges(WT)
    colors = Dict((E.tail,E.head)=>E.color for E in all_edges);
    interior_edges = filter(E->E.color ≠ :gray, all_edges)
    found_edges = Set{NTuple{2,T}}()
    all_faces = Face{T}[]
    u = blueroot(WT)
    v = cw(WT.P,u,greenroot(WT))
    while length(found_edges) < length(interior_edges)
        if colors[(u,v)] == :blue
            (u,v) = v,cw(WT.P,v,u)
        else
            (u,v) = u,cw(WT.P,u,v)
        end
        if (u,v) in found_edges
            continue
        end
        f = face(WT.P,u,v)
        push!(found_edges,edges(f)...)
        push!(all_faces,f)
    end
    return all_faces
end

"""
    schnyderwood(P::PlanarMap,
                 outface=face(P,1,neighbors(P,1)[1]))

Find a Schnyder wood associated with the plane
triangulation P, with given outer face
"""
function schnyderwood(RNG::AbstractRNG,
                      P::PlanarMap{T};
                      outface::Face=face(P,1,neighbors(P,1)[1])) where T<: Integer
    PWT = PartialWoodedTriangulation(P;outface=outface)
    while !allfound(PWT)
        v = rand(RNG,PWT.shellcandidates)
        shell!(PWT,v)
    end
    return WoodedTriangulation(PWT)
end

function schnyderwood(P::PlanarMap{T};
                      outface::Face=face(P,1,neighbors(P,1)[1])) where T<: Integer
      schnyderwood(Base.Random.globalRNG(),P,outface=outface)
end

struct PartialWoodedTriangulation{T<:Integer}
    P::PlanarMap{T}
    bluetree::Dict{T,T}
    redtree::Dict{T,T}
    greentree::Dict{T,T}
    blueroot::T
    redroot::T
    greenroot::T
    discovered::Set{T}
    shellcandidates::Set{T}
    fneighborcount::Vector{Int64}
end

function PartialWoodedTriangulation(P::PlanarMap{T};
                        outface::Face=face(P,1,neighbors(P,1)[1])) where T<:Integer
    blueroot, greenroot, redroot = outface # bgr order not a mistake
    bluetree = Dict{T,T}()
    redtree = Dict{T,T}()
    greentree = Dict{T,T}()
    discovered = Set([blueroot])
    shellcandidates = Set([blueroot])
    affectedvertices = Set{T}()
    # number of frontier neighbors of each vertex in P
    fneighborcount = zeros(Int64,length(P))
    for v in outface
        for w in neighbors(P,v)
            fneighborcount[w] += 1
        end
    end
    PartialWoodedTriangulation(
    P,bluetree,redtree,greentree,blueroot,redroot,greenroot,
    discovered,shellcandidates,fneighborcount
    )
end

function WoodedTriangulation(PWT::PartialWoodedTriangulation)
    if !allfound(PWT)
        error("Partial wooded triangulation not complete")
    end
    WoodedTriangulation(PWT.P,
            map(RootedTree,(PWT.bluetree, PWT.redtree, PWT.greentree),
                           (PWT.blueroot, PWT.redroot, PWT.greenroot))...)
end

function allfound(PWT::PartialWoodedTriangulation)
    length(PWT.P) == length(PWT.discovered) + 2
end

function shell!(PWT::PartialWoodedTriangulation{T},v::T) where T<:Integer
    push!(PWT.discovered,v)
    affectedvertices = Set{T}()
    if v == PWT.blueroot
        N = collect(rotate(neighbors(PWT.P,v),PWT.redroot))
    else
        N = downwardneighbors(neighbors(PWT.P,v),PWT.bluetree[v])
        filter!(x->x ∉ PWT.discovered,N)
        PWT.redtree[v] = N[1]
        PWT.greentree[v] = N[end]
    end
    for w in N[2:end-1]
        PWT.bluetree[w] = v
    end
    # update frontier neighbor count
    for w in neighbors(PWT.P,v)
        PWT.fneighborcount[w] -= 1
        push!(affectedvertices,w)
    end
    for u in N[2:end-1]
        for w in neighbors(PWT.P,u)
            PWT.fneighborcount[w] += 1
            push!(affectedvertices,w)
        end
    end
    # update active vertex list
    pop!(PWT.shellcandidates,v)
    for w in N
        if PWT.fneighborcount[w] == 2 && w ≠ PWT.redroot && w ≠ PWT.greenroot
            push!(PWT.shellcandidates,w)
        end
    end
    for w in affectedvertices
        if PWT.fneighborcount[w] ≠ 2 && w in PWT.shellcandidates
            pop!(PWT.shellcandidates,w)
        end
    end
    if 0 in PWT.shellcandidates
        error("hey!")
    end
end

function downwardneighbors(N::NeighborCycle{T},v::T) where T<:Integer
    p = position(N,v)
    N[to(p+1,p-1)]
end

"""
    schnyder_red(WT::WoodedTriangulation)

Walk the green tree to determine the red coordinate
(the number of faces between the green and blue flow lines
from each vertex)
"""
function schnyder_red(WT::WoodedTriangulation{T}) where T<:Integer
    # RD = red descendants, GA = green ancestors, BA = blue ancestors
    RD = backsum(WT.redtree,Dict(k=>1 for k in vertices(WT.redtree)))-1
    GA = treesum(WT.greentree,Dict(k=>1 for k in vertices(WT.greentree)))-1
    BA = treesum(WT.bluetree,Dict(k=>1 for k in vertices(WT.bluetree)))-1
    # BR = sum of red descendants with respect to blue tree
    # GR = sum of red descendants with respect to
    RD[blueroot(WT)] = 0; RD[greenroot(WT)] = 0
    BR = treesum(WT.bluetree,RD)
    GR = treesum(WT.greentree,RD)
    Dict(k=>2(BR[k] + GR[k] - RD[k]) + GA[k] + BA[k] - 1 for k in interiorvertices(WT))
end

function schnyder_green(WT::WoodedTriangulation{T}) where T<:Integer
    schnyder_red(WoodedTriangulation(WT.P,WT.redtree,WT.greentree,WT.bluetree))
end

function schnyder_blue(WT::WoodedTriangulation{T}) where T<:Integer
    schnyder_red(WoodedTriangulation(WT.P,WT.greentree,WT.bluetree,WT.redtree))
end

function schnyder_coords(WT::WoodedTriangulation)
    redcoords = schnyder_red(WT)
    greencoords = schnyder_green(WT)
    D = Dict(k=>(redcoords[k],greencoords[k]) for k in keys(redcoords))
    n = length(WT)
    merge!(D,Dict(blueroot(WT)=>(0,0),redroot(WT)=>(2n-5,0),greenroot(WT)=>(0,2n-5)))
    [D[k] for k=1:length(D)]
end

"""
    schnyderaverage(P::PlanarMap,
                    outface::Face,
                    n::Integer=100)
"""
function schnyderaverage(RNG::AbstractRNG,
                         P::PlanarMap{T};
                         outface::Face=face(P,1,neighbors(P,1)[1]),
                         n::Integer=100) where T<: Integer
    N = length(P)
    A = zeros(N)
    for i=1:n
        D = schnyder_coords(schnyderwood(RNG,P,outface=outface))
        A = [D[i] .+ A[i] for i=1:N]
    end
    [a./n for a in A]
end

function schnyderaverage(P::PlanarMap;kwargs...)
    schnyderaverage(Base.Random.globalRNG(),P;kwargs...)
end

"""
Non-mutating version of shell!
"""
function shell(P::PartialWoodedTriangulation{T},v::T) where T<:Integer
    Q = deepcopy(P)
    shell!(Q,v)
    Q
end

"""
    allwoods(PWT::PartialWoodedTriangulation,C::Channel)

Send to `C` a sequence of all `WoodedTriangulation`s on `P`
obtained by further shelling `PWT`
"""
function allwoods(PWT::PartialWoodedTriangulation,C::Channel)
    if allfound(PWT)
        put!(C,WoodedTriangulation(PWT))
    else
        for v in copy(PWT.shellcandidates)
            allwoods(shell(PWT,v),C)
        end
    end
end

"""
    allwoods(P::PlanarMap,C::Channel)

Send a sequence of all `WoodedTriangulation`s on `P`
to the channel `C`

# Examples:

```julia
# Count the number of Schnyder woods on a
# a randomly sampled wooded triangulation of size 10
sum(1 for W in Channel(C->allwoods(PlanarMap(UWT(10)),C)))
```
"""
function allwoods(P::PlanarMap,C::Channel;
                  outface::Face=face(P,1,neighbors(P,1)[1]))
    allwoods(PartialWoodedTriangulation(P;outface=outface),C)
end
