# author: Jianqi Huang
# Restuccia and Rogerson (2008)
begin
    using Plots
    using Parameters
    using SparseArrays
    using LinearAlgebra
    using DelimitedFiles
    using DataFrames
end


alpha=.85/3; #capital
gamma=.85/3*2; #labor
beta_dis=0.96; # the discount factor
delta=0.08;  #depreciation rate
nu=1.0; #energy price
S=1; #scale 

r=1/beta_dis-(1-delta)
R=r-delta

# read the data
firms_data=DataFrame(readdlm("data/establishment_dist.txt", '\t', header=false),[:firms_size,:firms_proportion])


ns = 100
ndata=size(firms_data,1)
s_max = firms_data[ndata,"firms_size"]^(1-gamma-alpha)  # z_ns represents the upper bound of firm's productivity
s_grid = 10.0.^(range(log10(1),log10(s_max),length=ns))
n_grid = s_grid.^(1.0/(1. - gamma -alpha))
theo_hs = zeros(ns)

# add additional new row: firms_size = 0, firms_proportion = 0
firms_data = vcat(firms_data, DataFrame(firms_size = 0, firms_proportion = 0))
# sort the firms_data by firms_size
firms_data = sort(firms_data, [:firms_size])

# calculate the theoretical firm size distribution
for i in 2:ndata+1  # Fixed: account for the added row
    println(i)
    I = findall((n_grid .<= firms_data[i,"firms_size"]) .& (n_grid .> firms_data[i-1,"firms_size"]) )
    println(I)
    if length(I) > 0  # Added safety check
        theo_hs[I] .= firms_data[i,"firms_proportion"] / length(I)
    end
end

# Calculate the CDF of the firms' data 
firms_data[!,:cum_firms_proportion] .= cumsum(firms_data[!,:firms_proportion])

# visualize the theoretical and empirical firm size distribution
plot(n_grid,cumsum(theo_hs),xaxis = :log10,xticks=[1, 10, 100, 1000, 10000],xlims=(1,10000),label="Theoretical",xlabel="Number of Employees (log scale)",ylabel="Cumulative Proportion of Firms")
scatter!(firms_data[!,:firms_size],firms_data[!,:cum_firms_proportion],label="Data",xlabel="Number of Employees (log scale)",ylabel="Cumulative Proportion of Firms")

# save the plot
savefig("figures/firm_size_distribution.pdf")

include("benchmark.jl")
include("plots.jl")
include("distorted.jl")  # Comment out for now to test benchmark first

