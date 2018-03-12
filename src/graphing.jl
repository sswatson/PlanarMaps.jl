
struct EmbeddedMap{T<:Integer}
    P::PlanarMap{T}
    locs::Dict{T,<:Tuple{<:Real,<:Real}}
end

function EmbeddedMap(P::PlanarMap,
                     locs::Array{<:Tuple{<:Real,<:Real}})
    EmbeddedMap(P,Dict(k=>v for (k,v) in enumerate(locs)))
end

function EmbeddedMap(RNG::AbstractRNG,
                     P::PlanarMap;
                     outeredge=(1,neighbors(P,1)[1]))
    T = triangulation(P,outeredge=outeredge)
    WT = schnyderwood(T)
    n = length(P)
    map((a,b)->ϕ(a,b,n),EmbeddedMap(P,schnyder_coords(WT)))
end

EmbeddedMap(P::PlanarMap;kwargs...) =
        EmbeddedMap(MersenneTwister(0),P;kwargs...)

import Base.map
function map(f::Function,E::EmbeddedMap)
    EmbeddedMap(E.P,Dict(k=>f(v...) for (k,v) in E.locs))
end

ϕ(a::Real,b::Real,n::Real) = tuple(([1 1/2; 0 sqrt(3)/2]*[a; b]/n)...)
ϕ(a::Real,b::Real) = ϕ(a,b,1)

function draw(E::EmbeddedMap;
              outeredge=(1,neighbors(E.P,1)[1]),
              edgepen=AsyPlots.Pen(),
              pointpen=AsyPlots.Pen(),
              labels=true,
              fontsize=12)
    Edges = [AsyPlots.Path([E.locs[i],E.locs[j]],pen=edgepen) for (i,j) in edges(E.P)]
    if labels
        Disks = [AsyPlots.Circle(E.locs[i],fontsize/250,
                fillpen=AsyPlots.Pen(color="white")) for i in 1:length(E.P)]
        Labels = [AsyPlots.Label(string(i),E.locs[i];
                pen=AsyPlots.Pen(fontsize=fontsize)) for i in 1:length(E.P)]
        AsyPlots.Plot(vcat(Edges,Disks,Labels))
    else
        Points = [AsyPlots.Point(E.locs[i],pen=pointpen) for i in 1:length(E.P)]
        AsyPlots.Plot(vcat(Edges,Points))
    end
end

function draw(P::PlanarMap;kwargs...)
    draw(EmbeddedMap(P);kwargs...)
end

function facecolor(s::Array{String,1})
    return sum([1/3 * AsyPlots.NamedColor(c) for c in s])
end

function draw(WT::WoodedTriangulation;
	      rotation=0,
              linewidth=1.0,
              pointsize=3.0,
              pointcolor="black",
              fillfaces=false,
              labels=true,
              coords=schnyder_coords(WT),
	      fontsize=12,
              kwargs...)
    ϕ(z) = cis(rotation)*(z[1] + 0.5*z[2] + im * sqrt(3)/2 * z[2])
    n = length(WT.P) - 3
    colors = Dict((a,b)=>c for (a,b,c) in edges(WT))
    for i=1:3
        colors[(n+mod1(i,3),n+mod1(i+1,3))] = :gray
        colors[(n+mod1(i+1,3),n+mod1(i,3))] = :gray
    end
    grlist = AsyPlots.GraphicElement2D[]
    # EDGES:
    for (tree,color) in zip((WT.bluetree,WT.redtree,WT.greentree),
                            ("Navy","Red","DarkGreen"))
        for k=1:n
            push!(grlist,AsyPlots.Path2D([ϕ(coords[k]),
                    ϕ(coords[tree[k]])];pen=AsyPlots.Pen(
                    color=AsyPlots.NamedColor(color),linewidth=linewidth)))
        end
    end
    # POINTS:
    if labels
        append!(grlist,vcat([[AsyPlots.Circle2D(ϕ(coords[k]),n*fontsize/250;
                pen=AsyPlots.Pen(linewidth=linewidth),
                fillpen=AsyPlots.Pen(color="white")),
                AsyPlots.Label2D("$k",ϕ(coords[k]);
                pen=AsyPlots.Pen(linewidth=pointsize,
                        color=AsyPlots.NamedColor(pointcolor),
                        fontsize=fontsize))] for k=1:length(coords)]...))
    else
        append!(grlist,[AsyPlots.Point2D(ϕ(coords[k]);
                pen=AsyPlots.Pen(linewidth=pointsize,
                color=AsyPlots.NamedColor(pointcolor))) for k=1:length(coords)])
    end
    if fillfaces
        append!(grlist,
                [AsyPlots.Polygon2D([ϕ(coords[k]) for k in fc.elements];
                    pen=AsyPlots.Pen(linewidth=linewidth),
                    fillpen=AsyPlots.Pen(color=facecolor(sort([colors[p] for p in pairs(fc)]))))
                        for fc in interiorfaces(faces(WT.P;outeredge=(n+1,n+3)))])
    end
    return AsyPlots.Plot(grlist;kwargs...)
end
