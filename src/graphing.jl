
struct EmbeddedMap{T<:Integer}
    P::PlanarMap{T}
    locs::Array{<:Tuple{<:Real,<:Real}}
end

function EmbeddedMap(RNG::Random.AbstractRNG,
                     P::PlanarMap;
                     outeredge=(1,neighbors(P,1)[1]),
                     corners=nothing)
    T = triangulation(RNG,P,outeredge=outeredge,corners=corners)
    outface = corners ≠ nothing ? Face(corners...) : face(T,outeredge...)
    WT = schnyderwood(RNG,T,outface=outface)
    n = length(P)
    map((a,b)->ϕ(a,b,n),EmbeddedMap(P,schnyder_coords(WT)))
end

EmbeddedMap(P::PlanarMap;kwargs...) =
        EmbeddedMap(Random.MersenneTwister(0),P;kwargs...)

length(E::EmbeddedMap) = length(E.P)

import Base.map
function map(f::Function,E::EmbeddedMap)
    EmbeddedMap(E.P,[f(a,b) for (a,b) in E.locs])
end

ϕ(a::Real,b::Real,n::Real) = tuple(([1 1/2; 0 sqrt(3)/2]*[a; b]/n)...)
ϕ(a::Real,b::Real) = ϕ(a,b,1)
ϕ(t::NTuple{2}) = ϕ(t...)

function bbox_size(E::EmbeddedMap)
    My = maximum(y for (x,y) in E.locs)
    my = minimum(y for (x,y) in E.locs)
    Mx = maximum(x for (x,y) in E.locs)
    mx = minimum(x for (x,y) in E.locs)
    max(Mx-mx,My-my)
end

struct GraphPenSet
    pointpen
    diskpen
    labelpen
    edgepen
    facepen
    pointpens
    diskpens
    labelpens
    edgepens
    facepens
end

function GraphPenSet(;kwargs...)
    args = DataStructures.OrderedDict(
            :pointpen=>AsyPlots.Pen(color="DarkGreen"),
            :diskpen=>AsyPlots.Pen(color="Black"),
            :labelpen=>AsyPlots.Pen(fontsize=12,color="White"),
            :edgepen=>AsyPlots.Pen(color="Navy",linewidth=1.2),
            :facepen=>AsyPlots.Pen(color="LightBlue"),
            :pointpens=>nothing,
            :diskpens=>nothing,
            :labelpens=>nothing,
            :edgepens=>nothing,
            :facepens=>nothing)

    merge!(args,Dict(kwargs))
    GraphPenSet(values(args)...)
end

for name in (:point,:disk,:label,:edge,:face)
    namepen = Symbol(string(name)*"pen")
    namepens = Symbol(string(name)*"pens")
    @eval function $namepen(G::GraphPenSet,k)
        if hasmethod(getindex,map(typeof,(G.$namepens,k)))
            G.$namepens[k]
        elseif !isempty(methods(G.$namepens))
            G.$namepens(k)
        else
            G.$namepen
        end
    end
end

function draw(E::EmbeddedMap;
              outeredge=(1,neighbors(E.P,1)[1]),
              drawpoints=true,
              drawedges=true,
              fillfaces=false,
              drawlabels=true,
              graphpenset=GraphPenSet(),
              kwargs...)

    if length(kwargs) > 0
        G = GraphPenSet(;kwargs...)
    else
        G = graphpenset
    end
    s = bbox_size(E) # for sizing the disk that the label goes in
    grlist = AsyPlots.GraphicElement2D[]
    if fillfaces
        append!(grlist,[AsyPlots.Polygon2D([E.locs[i] for i in fc],
                            pen=AsyPlots.NoPen(),
                            fillpen=facepen(G,fc)) for fc in interiorfaces(
                                faces(E.P;outeredge=outeredge))])
    end
    if drawedges
        append!(grlist,[AsyPlots.Path([E.locs[i],E.locs[j]],pen=edgepen(G,(i,j)))
                                                    for (i,j) in edges(E.P)])
    end
    if drawlabels
        append!(grlist,vcat([[AsyPlots.Circle(E.locs[i],
                                              s*labelpen(G,i).fontsize/500,
                                     pen=diskpen(G,i),
                                     fillpen=pointpen(G,i)),
                              AsyPlots.Label(string(i),E.locs[i];
                                     pen=labelpen(G,i))] for i in 1:length(E.P)]...))
    elseif drawpoints
        append!(grlist,[AsyPlots.Point(E.locs[i],
                                  pen=pointpen(G,i)) for i in 1:length(E.P)])
    end
    AsyPlots.Plot(grlist)
end

"""
    draw(P,outeredge,corners)

Draw the planar map P with outer face either:
    - the face to the left of the `outeredge`
    - the face containing the sequence of vertices `corners`

- `outeredge` is an ordered pair
- `corners` is either `nothing` (to indicate deference to
    specification via `outeredge`) or an ordered triple
"""
function draw(P::PlanarMap;
              outeredge=(1,neighbors(P,1)[1]),
              corners=nothing,
              kwargs...)
    draw(EmbeddedMap(P;outeredge=outeredge,corners=corners);kwargs...)
end

function draw(WT::WoodedTriangulation;
              rotation=0,
              linewidth=10/length(WT),
              pointsize=3,
              pointcolor=AsyPlots.NamedColor(0.65,0.65,0.65),
              linecolor="default",
              labelcolor="black",
              labelsize=max(1,floor(log10(length(WT))))*min(12,max(144/length(WT))),
              kwargs...)
    ASYCOLORS =     Dict(zip((:blue, :red, :green,     :gray),
        AsyPlots.NamedColor.(("Navy","Red","DarkGreen","Gray"))))
    edgedict = Dict((e.tail,e.head)=>e.color for e in edges(WT))
    edgepens(t) = linecolor == "default" ?
                        AsyPlots.Pen(linewidth=linewidth,
                                     color=ASYCOLORS[edgedict[t]]) :
                                AsyPlots.Pen(linewidth=linewidth,
                                             color=linecolor)
    function facepens(F::Face)
        AsyPlots.Pen(color=sum(1/3 * AsyPlots.NamedColor(
                            string(edgedict[e])) for e in edges(F)))
    end
    GPS = GraphPenSet(;edgepens=edgepens,
                       facepens=facepens,
                       diskpen=AsyPlots.Pen(linewidth=10/length(WT)),
                       pointpen=AsyPlots.Pen(linewidth=pointsize,
                                             color=pointcolor),
                       labelpen=AsyPlots.Pen(fontsize=labelsize,color=labelcolor))
    draw(EmbeddedMap(WT.P,
         map((t->reim(cis(rotation)*(t[1]+im*t[2])))∘ϕ,schnyder_coords(WT))),
         outeredge=(blueroot(WT),greenroot(WT));
         graphpenset=GPS,kwargs...)
end
