
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

function downwardneighbors(N::NeighborCycle{T},v::T) where T<:Integer
    p = position(N,v)
    N[to(p+1,p-1)]
end

"""
    schnyderwood(P::PlanarMap,
                 outface=face(P,1,neighbors(P,1)[1]))

Find a Schnyder wood associated with the plane
triangulation P, with given outer face

The algorithm used is the shelling method
"""
function schnyderwood(RNG::AbstractRNG,
                      P::PlanarMap{T};
                      outface::Face=face(P,1,neighbors(P,1)[1])) where T<: Integer
    blueroot, greenroot, redroot = outface # bgr order not a mistake
    bluetree = Dict{T,T}()
    redtree = Dict{T,T}()
    greentree = Dict{T,T}()
    last_discovered = blueroot
    discovered = Set([blueroot])
    affectedvertices = Set{T}()
    # number of frontier neighbors of each vertex in P
    fneighborcount = zeros(Int64,length(P))
    for v in outface
        for w in neighbors(P,v)
            fneighborcount[w] += 1
        end
    end
    activevertices = Set([blueroot])
    while length(discovered) < length(P) - 2
        v = rand(RNG,activevertices)
        push!(discovered,v)
        if v == blueroot
            N = collect(rotate(neighbors(P,v),redroot))
        else
            N = downwardneighbors(neighbors(P,v),bluetree[v])
            filter!(x->x ∉ discovered,N)
            redtree[v] = N[1]
            greentree[v] = N[end]
        end
        for w in N[2:end-1]
            bluetree[w] = v
        end
        # update frontier neighbor count
        for w in neighbors(P,v)
            fneighborcount[w] -= 1
            push!(affectedvertices,w)
        end
        for u in N[2:end-1]
            for w in neighbors(P,u)
                fneighborcount[w] += 1
                push!(affectedvertices,w)
            end
        end
        # update active vertex list
        pop!(activevertices,v)
        for w in N
            if fneighborcount[w] == 2 && w ≠ redroot && w ≠ greenroot
                push!(activevertices,w)
            end
        end
        for w in affectedvertices
            if fneighborcount[w] ≠ 2 && w in activevertices
                pop!(activevertices,w)
            end
        end
    end
    return WoodedTriangulation(P,
            map(RootedTree,(bluetree, redtree, greentree),
                            (blueroot, redroot, greenroot))...)
end

function schnyderwood(P::PlanarMap{T};
                      outface::Face=face(P,1,neighbors(P,1)[1])) where T<: Integer
      schnyderwood(Base.Random.globalRNG(),P,outface=outface)
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
    merge(D,Dict(blueroot(WT)=>(0,0),redroot(WT)=>(2n-5,0),greenroot(WT)=>(0,2n-5)))
end

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
