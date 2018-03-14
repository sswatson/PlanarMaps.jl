
#--------------------------------------------------------

function sample(RNG::AbstractRNG,w::Vector)
    w /= sum(w)
    p = rand(RNG)
    runningsum = 0
    for i=1:length(w)
        runningsum += w[i]
        if p < runningsum
            return i
        end
    end
end

sample(w::Vector) = sample(Base.Random.globalRNG(),w)

###---- FUNCTIONS FOR GENERATING CONDITIONED RANDOM WALK ----------

lbinom(n,k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
p0(a,c,n::BigInt) = mod(n+c-a,2) ≠ 0 ? 0.0 : exp(lbinom(n,(n+c-a)/2) + n * log(1/2))
p0float(a,c,n) = exp(lbinom(n,(n+c-a)/2) + n * log(1/2))
p0precise(a,c,n) = Float64(binomial(big(n),big(div(n+c-a,2))) * big(1/2)^n)
function p0(a,c,n::Int64)
    if mod(n+c-a,2) ≠ 0 || abs(a-c) > n
        return 0.0
    elseif n < 15_000
        return p0float(a,c,n)
    else
        # Use binomial approxmation
        k = div(n+c-a,2)
        return sqrt(2/(pi*n))*exp(-2*(k-n/2)^2/n)
    end
end
p(a,c,n) = p0(a,c,n) - p0(a,-c,n)
q(a,b,n) = p(b,3,n)*p(a,1,n) - p(b,1,n)*p(a,3,n)

function normalize(v)
    return v / sum(v)
end

function transweights(a,b,n)
    return normalize([q(i,j,n) for i=a-1:2:a+1, j=b-1:2:b+1])
end

function takestep(RNG::AbstractRNG,a,b,n)
    return [(i,j) for j=b-1:2:b+1 for i=a-1:2:a+1][
                            sample(RNG,transweights(a,b,n)[:])]
end

takestep(a,b,n) = takestep(Base.Random.globalRNG(),a,b,n)

function contour_to_tree(v::Array{Int64,1})
    n = div(length(v)-1,2)
    nbs = [Int64[] for _=1:n+3]
    stack = Int64[n+1]
    heads = Int64[]
    num_ones = 0
    bluetree = Dict{Int64,Int64}()
    for s in diff(v)
        if s == 1
            push!(stack,num_ones+1)
            num_ones += 1
        else
            e = pop!(stack)
            push!(heads,e)
            bluetree[e] = stack[end]
            push!(nbs[e],stack[end])
            push!(nbs[stack[end]],e)
        end
    end
    return PlanarMap(map(reverse,nbs)), bluetree, heads
end

function pair_dyck_paths(RNG::AbstractRNG,
                         n::Integer;
                         verbose=false)
    X = Int64[1]
    Y = Int64[3]
    power_of_ten = 10^floor(Integer,log10(n/2))
    for k=1:n
        (a,b) = takestep(RNG,X[end],Y[end],(n-k))
        push!(X,a)
        push!(Y,b)
    end
    X = X-1
    Y = Y-3
    return X,Y
end

pair_dyck_paths(n::Integer;kwargs...) =
    pair_dyck_paths(Base.Random.globalRNG(),n)
#----------------------------------------------------

#--- FUNCTIONS FOR GENERATING WOODS -----------------

import Base: keys, getindex
keys(R::RootedTree) = keys(R.ancestors)
getindex(R::RootedTree,k) = getindex(R.ancestors,k)

function add_red_tree(M::PlanarMap,
                      bluetree::Dict{Int64,Int64},
                      Y::Array{Int64,1},
                      heads::Array{Int64,1})
    nbs = deepcopy(M.nbs)
    headDict = Dict(map(reverse,enumerate(heads)))
    n = (length(Y)-1) ÷ 2
    tails = Int64[]
    num_ones = 0
    redtree = Dict{Int64,Int64}()
    tails = cumsum(diff(Y).==1)[find(diff(Y) .== -1)]
    tails = [t == n ? n+2 : t+1 for t in tails]
    unmatchedHeadInds = Set(1:length(heads))
    for tail in tails
        matchInd = tail == n+2 ? n : findfirst(x -> x ≥ tail, heads) - 1
        while matchInd ∉ unmatchedHeadInds
            matchInd -= 1
        end
        pop!(unmatchedHeadInds,matchInd)
        head = heads[matchInd]
        redtree[head] = tail
        nbs[head] = insertccw(nbs[head],bluetree[head],tail)
        if tail == n+2
            if length(nbs[tail]) == 0
                nbs[tail] = NeighborCycle([head])
            else
                nbs[tail] = insertcw(nbs[tail],nbs[tail][1],head)
            end
        else
            cwneighbor = bluetree[tail]
            while true
                cwneighbor = cw(M,tail,cwneighbor)
                if (cwneighbor == n+1 ||
                    bluetree[cwneighbor] == tail ||
                        bluetree[tail] == cwneighbor)
                    break
                end
            end
            nbs[tail] = insertccw(nbs[tail],cwneighbor,head)
        end
    end
    PlanarMap(nbs), redtree
end

function find_split(F::Face,
                    bluetree::Dict{Int64,Int64},
                    redtree::Dict{Int64,Int64})
    # finds which vertex in each face has two incident faces outgoing red and blue
    f = [F.elements; F.elements[1:1]]
    for (i,v) in enumerate(f)
        prv = f[i == 1 ? length(f)-1 : i - 1]
        nxt = f[i == length(f)-1 ? 1 : i + 1]
        if v in keys(redtree) && redtree[v] == prv && bluetree[v] == nxt
            return v, prv, nxt
        end
    end
    error("Split not found for face "*string(f))
    (0,0,0)
end

function rotate(a::Array{Int64,1},n::Integer)
    vcat(a[n:end],a[1:n-1])
end

import Base.getindex, Base.findfirst
getindex(Z::Base.Iterators.Zip2,k::Integer) = (Z.a[k],Z.b[k])
findfirst(C::Face,t::Tuple) = findfirst(pairs(C),t)

function add_green_tree(M::PlanarMap,
                        bluetree::Dict{Int64,Int64},
                        redtree::Dict{Int64,Int64})
    n = maximum(keys(bluetree))
    nbs = deepcopy(M.nbs)
    greentree = Dict{Int64,Int64}()
    found = false # whether outer face has been found yet
    FL = faces(M;outeredge=(n+1,1))
    # connect edges on left side of the planar map to
    # the green root:
    F = outerface(FL)
    leftedges = F[1:findfirst(F.elements,n+2)-1]
    for (i,w) in enumerate(leftedges[2:end])
        greentree[w] = n+3
        nbs[w] = insertcw(nbs[w],leftedges[i],n+3)
        if length(nbs[n+3]) == 0
            nbs[n+3] = NeighborCycle([w])
        else
            nbs[n+3] = insertccw(nbs[n+3],nbs[n+3][end],w)
        end
    end
    # triangulate each interior face with green edges
    for F in interiorfaces(FL)
        if length(F) > 3
            v, prv, nxt = find_split(F,bluetree,redtree)
            # idx stores index in the neighbor list of v
            # where incoming green edges should be inserted
            fc = rotate(F,findfirst(F,nxt))[2:end-1]
            for (i,w) in enumerate(fc[1:end-1])
                greentree[w] = v
                nbs[v] = insertcw(nbs[v],prv,w)
                nbs[w] = insertccw(nbs[w],fc[i+1],v)
            end
        end
    end
    return PlanarMap(nbs), greentree
end

function UWT(RNG::AbstractRNG,n::Integer)
    X,Y = pair_dyck_paths(RNG,2n)
    M, bluetree, heads = contour_to_tree(X)
    M, redtree = add_red_tree(M,bluetree,Y,heads)
    M, greentree = add_green_tree(M,bluetree,redtree)
    add_edge!(M,n+1,n+2)
    add_edge!(M,n+2,n+3)
    add_edge!(M,n+3,n+1)
    return WoodedTriangulation(M,map(RootedTree,
                    (bluetree,redtree,greentree),
                    (n+1,n+2,n+3))...)
end

UWT(n::Integer) = UWT(Base.Random.globalRNG(),n)

function descendants(d::Dict{Int64,Int64})
    n = maximum(keys(d))
    desc = zeros(Int64,n)
    for j = 1:n
        k = j
        while true
            k = d[k]
            if k < n+1
                desc[k] += 1
            else
                break
            end
        end
    end
    return desc
end

function flowline(v::Int64,d::Dict{Int64,Int64})
    n = maximum(keys(d))
    branch = Int64[v]
    while branch[end] < n+1
        push!(branch,d[branch[end]])
    end
    return branch
end

#------------------------------------------
