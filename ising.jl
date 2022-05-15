using BenchmarkTools;
using Debugger;
import Base.show
#using GLMakie <- this is causing problems

# Define a lattice struct using bools to make it as compact as possible in memory.
mutable struct Lattice{N} <: Unsigned
    a::BitArray{N}
end
# Define how show prints the lattice to the output.
function Base.show(io::IO, x::Lattice)
    print(io::IO, x.a)
end
# Define a constructor for conveniently initializing the lattice of arbitrary shape.
function Lattice(dims::Tuple{Vararg{Int}}, bool::Bool)
    out = bool ? trues(dims) : falses(dims)
    Lattice(out)
end
# For completeness, extend the methods to include working for boolean vectors
function Lattice(len::Int, bool::Bool)
    out = bool ? trues(len) : falses(len)
    Lattice(out)
end
# Now do the same for randomly generating spins up and down.
function hot_lattice(dims::Tuple{Vararg{Int}})
    Lattice(convert(BitArray, rand(Bool, dims)))
end
# Entending the function above to work for vectors.
function hot_lattice(len::Int)
    Lattice(convert(BitArray, rand(Bool, len)))
end

## Functions for working with the lattice and calculating energy and so forth.
function binmap(a::Bool)
    a ? Int8(1) : Int8(-1)
end

function energy2D(lat::Lattice)
    # This computes the equation:
    # H(σ) = -Σ Jᵢⱼ σᵢσⱼ; where the interation 
    # is only felt between nearest neighbors.
    array = lat.a
    J = 1.0
    dims = size(array)
    if length(dims) != 2
        error("This function only works for 2D lattices!")
    end

    local energy = Int8(0)
    for i in 1:dims[2], j in 1:dims[1]
        energy += (binmap(array[j, i]) * binmap(array[mod1(j + 1, dims[1]), mod1(i + 1, dims[2])])
                   + binmap(array[j, i]) * binmap(array[mod1(j - 1, dims[1]), mod1(i + 1, dims[2])])
                   + binmap(array[j, i]) * binmap(array[mod1(j + 1, dims[1]), mod1(i - 1, dims[2])])
                   + binmap(array[j, i]) * binmap(array[mod1(j - 1, dims[1]), mod1(i - 1, dims[2])]))
    end

    return -J * energy / 2
end

function energy2D(array::BitMatrix)
    # This computes the equation:
    # H(σ) = -Σ Jᵢⱼ σᵢσⱼ - hσᵢ; where the interation 
    # is only felt between nearest neighbors.
    J = 1.0
    dims = size(array)
    if length(dims) != 2
        error("This function only works for 2D lattices!")
    end

    local energy = Int8(0)
    for i in 1:dims[2], j in 1:dims[1]
        energy += (binmap(array[j, i]) * binmap(array[mod1(j + 1, dims[1]), mod1(i + 1, dims[2])])
                   + binmap(array[j, i]) * binmap(array[mod1(j - 1, dims[1]), mod1(i + 1, dims[2])])
                   + binmap(array[j, i]) * binmap(array[mod1(j + 1, dims[1]), mod1(i - 1, dims[2])])
                   + binmap(array[j, i]) * binmap(array[mod1(j - 1, dims[1]), mod1(i - 1, dims[2])]))
    end

    return -J * energy / 2
end

function step!(lat::Lattice, temp::Float64)
    # This function implements the metropolis algorithm and 
    # iterates for a given temperature.
    β = 1 / temp
    array = lat.a
    dims = size(array)
    enpast = energy2D(array)
    x, y = rand(1:dims[1]), rand(1:dims[2])
    # Try flipping the bit at the x, y location
    array[x, y] = !array[x, y]
    en = energy2D(array)
    p = rand()
    # If the new energy is lower, keep it. Otherwise,
    # keep it with probability exp(-(en - enpast) * β).
    if (en - enpast) > 0
        #println(p)
        #println(en - enpast)
        if p <= exp(-(en - enpast) * β)
            #println("keep the change!")
        else
            array[x, y] = !array[x, y]
            #println("undo the change!")
        end
    else
        #println("the energy lowered!")
    end
    return array
end

function steps!(lat::Lattice, temps::AbstractVector{Float64})
    for temp in temps
        step!(lat, temp)
    end
    #display(heatmap(lat.a, legend=false))
    return lat.a
end

function steps!(lat::Lattice, temp::Float64, nsteps::Int)
    for i in 1:nsteps
        step!(lat, temp)
    end
    #display(heatmap(lat.a, legend=false))
    return lat.a
end

function plot_ising(lat::Lattice)
    array = lat.a
    heatmap(array, legend=false, ticks=false)
end