## An animation of a random walk on a graph.

## We are going to use a GNR graph.

##
using Plots
using PlotRecipes
using NearestNeighbors
#using Roots

##
function gnr(n,r)
  xy = rand(2,n)
  T = BallTree(xy)
  idxs = inrange(T, xy, r)
  # form the edges for sparse
  ei = Int[]
  ej = Int[]
  for i=1:n
    for j=idxs[i]
      if i > j
        push!(ei,i)
        push!(ej,j)
      end
    end
  end
  return xy, ei, ej
end


srand(9)
xy, ei, ej = gnr(50,0.2)

pyplot(size=(300,300))
graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
  markercolor=:black, markerstrokecolor=:white,
  markersize=4, linecolor=1, linealpha=0.8, linewidth=0.7,
  axis_buffer=0.02, background=nothing)
gui()
#savefig("gnr-figure.pdf")

## Simulate the random walk

srand(1)
function rw_on_graph(ei, ej, xy, start, nstep)
  eirw = [ei; ej]
  ejrw = [ej; ei]

  X = start
  # random walk step
  p = graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
    markercolor=:black, markerstrokecolor=:white,
    markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
    axis_buffer=0.02, background=:white, framestyle=:none)
  p1 = scatter!([xy[1,X]],[xy[2,X]],markersize=4, color=1)
  p2 = scatter!([xy[1,X]],[xy[2,X]],markersize=8, color=:orange)

  anim = @animate for i=1:nstep
    Xn = rand(ejrw[eirw .== X])
    p[3] = [xy[1,X]],[xy[2,X]]
    p[4] = [xy[1,Xn]],[xy[2,Xn]]
    X = Xn
  end
  return anim
end

anim = rw_on_graph(ei, ej, xy, 1, 100)
gif(anim, "rw-on-gnr.gif", fps=4)

##

srand(1)
function nbrw_on_graph(ei, ej, xy, start, nstep)
  eirw = [ei; ej]
  ejrw = [ej; ei]

  X = start
  Xp = X
  X = rand(ejrw[eirw .== X]) # next step

  # random walk step
  p = graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
    markercolor=:black, markerstrokecolor=:white,
    markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
    axis_buffer=0.02, background=:white, framestyle=:none)
  p1 = scatter!([xy[1,X]],[xy[2,X]],markersize=4, color=1)
  p2 = scatter!([xy[1,X]],[xy[2,X]],markersize=8, color=:orange)

  anim = @animate for i=1:nstep
    Xn = rand(ejrw[eirw .== X])
    while Xn == Xp
      Xn = rand(ejrw[eirw .== X])
    end

    p[3] = [xy[1,X]],[xy[2,X]]
    p[4] = [xy[1,Xn]],[xy[2,Xn]]
    Xp = X
    X = Xn
  end
  return anim
end

anim = nbrw_on_graph(ei, ej, xy, 1, 100)
gif(anim, "nbrw-on-gnr.gif", fps=4)


srand(1)
function nbrw_on_graph(ei, ej, xy, start, nstep)
  eirw = [ei; ej]
  ejrw = [ej; ei]

  X = start
  Xp = X
  X = rand(ejrw[eirw .== X]) # next step

  # random walk step
  p = graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
    markercolor=:black, markerstrokecolor=:white,
    markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
    axis_buffer=0.02, background=:white, framestyle=:none)
  p1 = scatter!([xy[1,X]],[xy[2,X]],markersize=4, color=1)
  p2 = scatter!([xy[1,X]],[xy[2,X]],markersize=8, color=:orange)

  anim = @animate for i=1:nstep
    Xn = rand(ejrw[eirw .== X])
    while Xn == Xp
      Xn = rand(ejrw[eirw .== X])
    end

    p[3] = [xy[1,X]],[xy[2,X]]
    p[4] = [xy[1,Xn]],[xy[2,Xn]]
    Xp = X
    X = Xn
  end
  return anim
end

anim = nbrw_on_graph(ei, ej, xy, 1, 100)
gif(anim, "nbrw-on-gnr.gif", fps=4)
