using StatsBase
using MAT

""" Simulate a first-order Markov chain with transition probabilities given by
the column stochastic matrix P. """
function sim_FOMC(P::Array{Float64,2}, start::Int64, niter::Int64=10_000_000)
    n = size(P, 1)
    counts = zeros(Int64, n)
    states = [start]
    for _ in 1:niter
        col = P[:, states[end]]
        next = sample(1:n, Weights(col))
        counts[next] += 1
        push!(states, next)
    end
    return (counts, states)
end

""" Simulate a second-order Markov chain with transition probabilities given by
the transition probability tensor P. """
function sim_SOMC(P::Array{Float64,3}, start::NTuple{2,Int64}, niter::Int64=10_000_000)
    n = size(P, 1)
    counts = zeros(Int64, n)
    states = [start[1], start[2]]
    for _ in 1:niter
        col = vec(P[:, states[end], states[end - 1]])
        assert(sum(col) ≈ 1)
        next = sample(1:n, Weights(col))
        counts[next] += 1
        push!(states, next)
    end
    return (counts, states)
end

""" Simulate a spacey random walk from a second-order Markov chain with
transition probabilities given by the transition probability tensor P. The
column vector zerocol specifies what transitions to make in the case of an
undefined transition. """
function sim_SOSRW(P::Array{Float64,3}, start::Int64, zerocol::Vector{Float64},
                   super_spacey::Bool=false, niter::Int64=10_000_000)
    if !super_spacey; assert(sum(zerocol) ≈ 1); end
    n = size(P, 1)
    counts = zeros(Int64, n)
    w_counts = ones(Int64, n)
    states = [start]
    for _ in 1:niter
        second_last = sample(1:n, Weights(w_counts))
        col = vec(P[:, states[end], second_last])
        next = nothing
        if maximum(col) > 0.0
            next = sample(1:n, Weights(col))
        else
            if super_spacey; next = sample(1:n, Weights(w_counts))
            else             next = sample(1:n, Weights(zerocol))
            end
        end
        counts[next] += 1
        w_counts[next] += 1
        push!(states, next)
    end
    return (counts, states)
end

""" Compute Tx^2. """
function apply(T::Array{Float64,3}, x::Vector{Float64})
    n = size(T, 1)
    y = zeros(Float64, n)
    for j = 1:n, k = 1:n
        y += T[:, j, k] * x[j] * x[k]
    end
    return y
end
