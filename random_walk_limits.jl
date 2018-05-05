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

function spectral_order(ei,ej,xy)
  n = maximum(ei)
  A = sparse(ei,ej,1,n,n)
  A = A + A'
  L = diagm(vec(sum(A,2))) - A
  lams,Vs = eigs(L,nev=2, which=:SR)
  v = Vs[:,2]
  p = sortperm(v)
  @show lams
  ip = Vector{Int}(n)
  ip[p] = 1:n
  return p,ip
end
p,ip = spectral_order(ei, ej, xy)

xy1 = copy(xy)
ei1 = copy(ei)
ei = ip[ei]
ej = ip[ej]
xy = xy[:,p]

pyplot(size=(300,300))
graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
  markercolor=:black, markerstrokecolor=:white,
  markersize=4, linecolor=1, linealpha=0.8, linewidth=0.7,
  axis_buffer=0.02, background=nothing)
for i=1:50
  annotate!(xy[1,i],xy[2,i],"$i")
end
gui()
#savefig("gnr-figure.pdf")

## Simulate the random walk and show the distribution
pyplot(size=(400,300))
l = @layout([a{0.75w} b])
x = ones(50)/50
plot(graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
  markercolor=:black, markerstrokecolor=:white,
  markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
  axis_buffer=0.02, background=:white, framestyle=:none),
  scatter(x, 1:50, label="", framestyle=:none), layout=l)
gui()
##
srand(1)
pyplot(size=(400,300))
function rw_on_graph(ei, ej, xy, start, nstep)
  eirw = [ei; ej]
  ejrw = [ej; ei]
  N = maximum(ei)
  X = start
  # random walk step
  ## Simulate the random walk and show the distribution

  l = @layout([a{0.75w} b])
  x = ones(N)
  ps = plot(graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
    markercolor=:black, markerstrokecolor=:white,
    markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
    axis_buffer=0.02, background=:white, framestyle=:none),
    scatter(x/N, 1:50, label="", framestyle=:none), layout=l)
  p1 = ps[1]
  p2 = ps[2]
  p1a = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=4, color=1)
  p1b = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=8, color=:orange)
  anim = @animate for i=1:nstep
    Xn = rand(ejrw[eirw .== X])
    x[Xn] += 1
    ps[4] = [xy[1,X]],[xy[2,X]]
    ps[5] = [xy[1,Xn]],[xy[2,Xn]]
    ps[3] = x/(N+i),collect(1:50)
    X = Xn
  end
  return anim
end

anim = rw_on_graph(ei, ej, xy, 1, 500)
gif(anim, "rw-with-dist-on-gnr.gif", fps=25)

##
symmat = x->x+x'
d = vec(sum(symmat(sparse(ei,ej,1,50,50)),1))
scatter(d/sum(d), 1:50,  label="", framestyle=:none, size=(75,300))
savefig("rw-stationary-dist.pdf")

## Show what happens with a VRRW (the classic)
using StatsBase
srand(1)
pyplot(size=(400,300))
function vrrw_on_graph(ei, ej, xy, start, nstep)
  eirw = [ei; ej]
  ejrw = [ej; ei]
  N = maximum(ei)
  X = start
  # random walk step
  ## Simulate the random walk and show the distribution

  l = @layout([a{0.75w} b])
  x = ones(N)
  ps = plot(graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
    markercolor=:black, markerstrokecolor=:white,
    markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
    axis_buffer=0.02, background=:white, framestyle=:none),
    scatter(x/N, 1:50, label="", framestyle=:none), layout=l)
  p1 = ps[1]
  p2 = ps[2]
  p1a = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=4, color=1)
  p1b = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=8, color=:orange)
  anim = @animate for i=1:nstep

    neighs = ejrw[eirw .== X]
    weights = x[neights]


    Xn = neighs[sample(Weights(weights))]
    x[Xn] += 1
    ps[4] = [xy[1,X]],[xy[2,X]]
    ps[5] = [xy[1,Xn]],[xy[2,Xn]]
    ps[3] = x/(N+i),collect(1:50)
    X = Xn
  end
  return anim
end

anim = rw_on_graph(ei, ej, xy, 1, 500)
gif(anim, "vrrw-with-dist-on-gnr.gif", fps=25)


##
srand(1)
pyplot(size=(400,300))
function nbrw_on_graph(ei, ej, xy, start, nstep)
  eirw = [ei; ej]
  ejrw = [ej; ei]
  N = maximum(ei)
  X = start
  Xp = X
  X = rand(ejrw[eirw .== X]) # next step
  # random walk step
  l = @layout([a{0.75w} b])
  x = ones(N)
  x[X] += 1
  ps = plot(graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
    markercolor=:black, markerstrokecolor=:white,
    markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
    axis_buffer=0.02, background=:white, framestyle=:none),
    scatter(x/N, 1:50, label="", framestyle=:none), layout=l)
  p1 = ps[1]
  p2 = ps[2]
  p1a = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=4, color=1)
  p1b = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=8, color=:orange)
  anim = @animate for i=1:nstep
    Xn = rand(ejrw[eirw .== X])
    while Xn == Xp
      Xn = rand(ejrw[eirw .== X])
    end
    x[Xn] += 1
    ps[4] = [xy[1,X]],[xy[2,X]]
    ps[5] = [xy[1,Xn]],[xy[2,Xn]]
    ps[3] = x/(N+i+1),collect(1:50)
    Xp = X
    X = Xn
  end
  return anim
end

anim = nbrw_on_graph(ei, ej, xy, 1, 50)
gif(anim, "nbrw-with-dist-on-gnr.gif", fps=25)

##
##
srand(1)
pyplot(size=(400,300))
function nbsrw_on_graph(ei, ej, xy, start, nstep)
  eirw = [ei; ej]
  ejrw = [ej; ei]
  N = maximum(ei)
  X = start
  Xp = X
  X = rand(ejrw[eirw .== X]) # next step
  # random walk step
  l = @layout([a{0.75w} b])
  x = ones(N)
  x[X] += 1
  ps = plot(graphplot(ei, ej, x =xy[1,:], y=xy[2,:],
    markercolor=:black, markerstrokecolor=:white,
    markersize=4, linecolor=2, linealpha=0.8, linewidth=0.7,
    axis_buffer=0.02, background=:white, framestyle=:none),
    scatter(x/N, 1:50, label="", framestyle=:none), layout=l)
  p1 = ps[1]
  p2 = ps[2]
  p1a = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=4, color=1)
  p1b = scatter!(p1, [xy[1,X]],[xy[2,X]],markersize=8, color=:orange)
  anim = @animate for i=1:nstep
    Xn = rand(ejrw[eirw .== X])
    while Xn == Xp
      Xn = rand(ejrw[eirw .== X])
    end
    x[Xn] += 1
    ps[4] = [xy[1,X]],[xy[2,X]]
    ps[5] = [xy[1,Xn]],[xy[2,Xn]]
    ps[3] = x/(N+i+1),collect(1:50)
    Xp = X
    X = Xn
  end
  return anim
end

anim = nbrw_on_graph(ei, ej, xy, 1, 50)
gif(anim, "nbrw-with-dist-on-gnr.gif", fps=25)
