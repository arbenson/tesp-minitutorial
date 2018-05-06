##
using Plots
pyplot()

## Illustrate Peterson's case where the SRW does not converge.

M1 = [0.98 0.98 0.98; 0.01 0.01 0.01; 0.01 0.01 0.01 ]
M2 = [0.01 0.01 0.01; 0.98 0.98 0.98; 0.01 0.01 0.01 ]
M3 = [0.01 0.01 0.01; 0.01 0.01 0.01; 0.98 0.98 0.98 ]
Md = [0.98 0.01 0.01; 0.01 0.98 0.01; 0.01 0.01 0.98 ]

P = zeros(3,3,3,3,3)
P[:,:,3,3,3] = M1
P[:,:,1,1,2] = M1
P[:,:,1,2,1] = M1
P[:,:,2,1,1] = M1
P[:,:,1,1,3] = M1
P[:,:,1,3,1] = M1
P[:,:,3,1,1] = M1

P[:,:,1,1,1] = M2
P[:,:,2,2,1] = M2
P[:,:,2,1,2] = M2
P[:,:,1,2,2] = M2
P[:,:,2,2,3] = M2
P[:,:,2,3,2] = M2
P[:,:,3,2,2] = M2

P[:,:,2,2,2] = M3
P[:,:,3,3,1] = M3
P[:,:,3,1,3] = M3
P[:,:,1,3,3] = M3
P[:,:,3,3,2] = M3
P[:,:,3,2,3] = M3
P[:,:,2,3,3] = M3

P[:,:,1,2,3] = Md
P[:,:,1,3,2] = Md
P[:,:,2,1,3] = Md
P[:,:,2,3,1] = Md
P[:,:,3,1,2] = Md
P[:,:,3,2,1] = Md

##
function example_tensor()
    P = zeros(3,3,3,3,3)
    P1 = [[0.98, 0.01, 0.01] [0.98, 0.01, 0.01] [0.98, 0.01, 0.01]]
    P2 = [[0.01, 0.98, 0.01] [0.01, 0.98, 0.01] [0.01, 0.98, 0.01]]
    P3 = [[0.01, 0.01, 0.98] [0.01, 0.01, 0.98] [0.01, 0.01, 0.98]]
    P4 = [[0.98, 0.01, 0.01] [0.01, 0.98, 0.01] [0.01, 0.01, 0.98]]

    for trip in [(2,2,2), (0,0,1), (0,1,0), (1,0,0), (0,0,2), (0,2,0), (2,0,0)]
        tripp1 = (trip[1] + 1, trip[2] + 1, trip[3] + 1)
        P[:,:,tripp1...] = P1
    end
    for trip in [(0,0,0), (1,1,0), (1,0,1), (0,1,1), (1,1,2), (1,2,1), (2,1,1)]
        tripp1 = (trip[1] + 1, trip[2] + 1, trip[3] + 1)
        P[:,:,tripp1...] = P2
    end
    for trip in [(1,1,1), (2,2,0), (2,0,2), (0,2,2), (2,2,1), (2,1,2), (1,2,2)]
        tripp1 = (trip[1] + 1, trip[2] + 1, trip[3] + 1)
        P[:,:,tripp1...] = P3
    end
    for trip in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]
        tripp1 = (trip[1] + 1, trip[2] + 1, trip[3] + 1)
        P[:,:,tripp1...] = P4
    end

    return P
end
P = example_tensor()

## Simulate the dynamical system


## Simulate the SRW
# take a step on a SRW model with tensor P, history x, and current state i
function srw_step(P,x,i)
    nhist = ndims(P)-2
    inds = zeros(Int,nhist+1)
    inds[1] = i
    w = Weights(x)
    for j=2:length(inds)
        inds[j] = sample(w)
    end
    #@show inds
    #@show P[:,inds...]
    return sample(Weights(P[:,inds...]))
end

x = ones(3)
N = 3
plt = plot(1,
    background=:white,  legend=false, marker=(1.0,stroke(0)), line=(0),
    alpha=0.1, size=(300,300), xlim=(0,1), ylim=(0,1))
#scatter3d!([x[1]/N],[x[2]/N],[x[3]/N], markersize=1.5)
X = 1
step = 0
#steps = [1:25; 25:5:100; 100:25:1000; 1000:100:10000; 10000:1000:100000; 100000:10000:1000000]]
steps = [1:25; 25:5:100; 100:25:1000; 1000:100:10000; 10000:1000:100000; 100000:10000:1000000; 1000000:100000:10000000]
anim = @animate for i=1:length(steps)
    while step < steps[i]
        x[X] += 1
        N += 1
        X = srw_step(P,x/N,X)
        step += 1
    end

    push!(plt, x[1]/N, x[2]/N)

    title!("Step $(step)")
    #plt[2] = [x[1]/N],[x[2]/N],[x[3]/N]
end
gif(anim, "peterson-counterexample.gif", fps=30)

## This doesn't seem to be a counterexample -- what have I done wrong?
# This was all okay, there was a typo in srw_step


function mult5(T::Array{Float64,5}, x::Vector{Float64})
    dims = size(T)
    M = zeros(dims[1],dims[2])
    for i=1:dims[3], j=1:dims[4], k=1:dims[5]
        M += T[:,:,i,j,k] * x[i] * x[j] * x[k]
    end
    return M
end

function perron_vector(M::Array{Float64,2})
    d,V = eig(M)
    d = real.(d)
    V = real.(V)
    inds = sortperm(d, rev=true)[1:2]
    ind = inds[1]
    v = V[:, ind]
    if d[ind] - d[inds[2]] <= 1.0e-8
      warn("Non unique perron vector")
    end
    if v[1] < 0; v = -v; end
    v /= norm(v, 1)
end
M = mult5(P, x/N)
perron_vector(M)

##
x = ones(3)
X = 1
for i=1:10^8
    x[X] += 1
    N += 1
    X = srw_step(P,x/N,X)
end

##

x = ones(3)
plt = plot(1,
    background=:white,  legend=false, marker=(1.0,stroke(0)), line=(0),
    alpha=0.1, size=(300,300), xlim=(0,1), ylim=(0,1))
#scatter3d!([x[1]/N],[x[2]/N],[x[3]/N], markersize=1.5)
x = [0.85;0.1;0.1]
x = x/sum(x)
h = 0.01
N = 10000
anim = @animate for i=1:N
    x += h*(perron_vector(mult5(P,x))-x)
    x = x / sum(x) # keep stochastic

    push!(plt, x[1], x[2])

    title!("Step $(i)")
    #plt[2] = [x[1]/N],[x[2]/N],[x[3]/N]
end every 100
gif(anim, "peterson-counterexample-system.gif", fps=30)


## Make the picture for another version of the stochastic process
# take a step on a SRW model with tensor P, history x, and current state i
function srw_step_2(P,x,i)
    nhist = ndims(P)-1
    inds = zeros(Int,nhist)
    w = Weights(x)
    for j=1:length(inds)
        inds[j] = sample(w)
    end
    # randomly place our state
    inds[rand(1:nhist)] = i
    return sample(Weights(P[:,inds...]))
end

x = ones(3)
N = 3
plt = plot(1,
    background=:white,  legend=false, marker=(1.0,stroke(0)), line=(0),
    alpha=0.1, size=(300,300), xlim=(0,1), ylim=(0,1))
#scatter3d!([x[1]/N],[x[2]/N],[x[3]/N], markersize=1.5)
X = 1
step = 0
#steps = [1:25; 25:5:100; 100:25:1000; 1000:100:10000; 10000:1000:100000; 100000:10000:1000000]]
steps = [1:25; 25:5:100; 100:25:1000; 1000:100:10000; 10000:1000:100000; 100000:10000:1000000; 1000000:100000:10000000]
anim = @animate for i=1:length(steps)
    while step < steps[i]
        x[X] += 1
        N += 1
        X = srw_step_2(P,x/N,X)
        step += 1
    end

    push!(plt, x[1]/N, x[2]/N)

    title!("Step $(step)")
    #plt[2] = [x[1]/N],[x[2]/N],[x[3]/N]
end
gif(anim, "peterson-example-2.gif", fps=30)

## Show the dynamical system for this second example

function mult5_average(T,x)
    dims = size(T)
    M = zeros(dims[1],dims[2])
    for i=1:dims[3], j=1:dims[4], k=1:dims[5]
        M += T[:,:,i,j,k] * x[i] * x[j] * x[k]
    end

    for i=1:dims[3], j=1:dims[4], k=1:dims[5]
        M += T[:,i,:,j,k] * x[i] * x[j] * x[k]
    end

    for i=1:dims[3], j=1:dims[4], k=1:dims[5]
        M += T[:,i,j,:,k] * x[i] * x[j] * x[k]
    end

    for i=1:dims[3], j=1:dims[4], k=1:dims[5]
        M += T[:,i,j,k,:] * x[i] * x[j] * x[k]
    end

    return M/4
end

x = ones(3)
plt = plot(1,
    background=:white,  legend=false, marker=(1.0,stroke(0)), line=(0),
    alpha=0.1, size=(300,300), xlim=(0,1), ylim=(0,1))
#scatter3d!([x[1]/N],[x[2]/N],[x[3]/N], markersize=1.5)
x = [0.85;0.1;0.1]
x = x/sum(x)
h = 0.01
N = 10000
anim = @animate for i=1:N
    x += h*(perron_vector(mult5_average(P,x))-x)
    x = x / sum(x) # keep stochastic

    push!(plt, x[1], x[2])

    title!("Step $(i)")
    #plt[2] = [x[1]/N],[x[2]/N],[x[3]/N]
end every 100
gif(anim, "peterson-example-2-system.gif", fps=30)
