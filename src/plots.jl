using Plots

# Compute the normalized invariant distribution of plants.
muhat = 1 / lambda_exit .* xbar .* g

# Find mass of entry that clears the labor market.
N_hat = sum(n_bar .* muhat) 
# since here multiply the muhat, We still need to multiply by the E to get the total labor
# N*E represents the total number of workers in the economy
E = 1 / N_hat
mu = muhat * E
mu_s = sum(mu, dims=1) # Note sum(.,2) sums over rows #100*1

# Compute aggregate statistics.
Y = sum(smatrix .* k_bar.^alpha .* n_bar.^gamma .* mu)
K = sum(k_bar .* mu)
KY = K / Y
Kbe = K
# Average employment per plant (AEPP)
normalization_emp = n_bar[1, 3]
# let the optimal employee number of the lowest productivity plant with no tax to be 1
rnbar = n_bar / normalization_emp
AEPP = 1 / (M * normalization_emp)
# since the total labor supply is 1 and the total numbers of the operating firms is M,
# so 1/M is the average number of employee. Finally, divided by the normalization factor.

# Distribution and plant statistics
relmus = sum(mu, dims=2) / M #100*1
# Recall M = np.sum(np.sum(mu))
sK = (k_bar .* mu) ./ K
sN = (n_bar .* mu) ./ (N_hat .* E)
sY = (smatrix .* k_bar.^alpha .* n_bar.^gamma .* mu) ./ Y

ntax=5 # Number of experiments: tau=0,0.1,0.2,0.3,0.4.
ntau=3 # Define ntau which is used in benchmark.jl

# Fixed: Julia array indexing and syntax
modeln = rnbar[:, 1]  # number of workers with lowest s in BE norm to 1
vec_rnbar = vec(rnbar')
I = sortperm(vec_rnbar)  # Fixed: use vec() to flatten
modelnp = zeros(ns * ntax, ntax)
modelnp[I, 1] = vec_rnbar[I]  # Fixed: Julia uses 1-based indexing
modelsNp = zeros(ns * ntax, ntax)
modelsNp[I, 1] = vec(sN)[I]  # Fixed: use vec() instead of flatten()
modelmup = zeros(ns * ntax, ntax)
modelmup[I, 1] = vec(mu)[I] / sum(mu)  # Fixed: use vec() instead of flatten()
modelmu = mu[:, 2] / sum(mu[:, 2])  # prob dist of plants
modelsN = sN[:, 2]


# initialize matrices for printing 
Yp = zeros(1, ntax)
Kp = zeros(1, ntax)
KYp = zeros(1, ntax)
Ap = zeros(1, ntax)
Ep = zeros(1, ntax)
Mp = zeros(1, ntax)
wp = zeros(1, ntax)
sgdpp = zeros(1, ntax)
tausp = zeros(1, ntax)
SYp = zeros(1, ntax)
relmup = zeros(ns, ntax * ntau)
relmusp = zeros(ns, ntax)
sKp = zeros(ns, ntax * ntau)
sNp = zeros(ns, ntax * ntau)
sYp = zeros(ns, ntax * ntau)
kbarp = zeros(ns, ntax * ntau)
nbarp = zeros(ns, ntax * ntau)
xbarp = zeros(ns, ntax * ntau)
AEPPp = zeros(1, ntax)
KsKp = zeros(1, ntax)

Yp[1, 1] = Y  # save the benchmark result in the first index in the Yp vector.
Kp[1, 1] = K
KYp[1, 1] = KY
Mp[1, 1] = M
Ap[1, 1] = A
Ep[1, 1] = E
wp[1, 1] = w
sgdpp[1, 1] = sgdp
AEPPp[1, 1] = AEPP
relmusp[:, 1] = relmus

I5 = findall(modeln .< 5)  # plants with less than 5 workers
I50 = findall((modeln .>= 5) .& (modeln .< 50))  # plants with 5 to less than 50 workers
Irest = findall(modeln .>= 50)  # plants with 50 workers or more
dis_plants = zeros(3, ntax)
dis_plants[:, 1] = [sum(modelmu[I5]), sum(modelmu[I50]), sum(modelmu[Irest])]
dis_sN = zeros(3, ntax)
dis_sN[:, 1] = [sum(modelsN[I5]), sum(modelsN[I50]), sum(modelsN[Irest])]



fig = 1
if fig == 1
    # 第一个图形：双子图柱状图
    p1 = plot(layout = @layout([a; b]), size=(600,800))
    
    # 子图1（顶部）
    bar!(p1[1], ["<5", "5-50", "50 or more"], dis_plants[:,1], 
         label="", title="Benchmark Economy", ylabel="Share of Establishments")
    
    # 子图2（底部）
    bar!(p1[2], ["<5", "5-50", "50 or more"], dis_sN[:,1], 
         label="", xlabel="Establishment Size (by Number of Workers)", 
         ylabel="Share of Employment")
    
    # 第二个图形：累积产出分布
    p2 = plot(log.(n_grid), cumsum(sY[:,2]), 
         linewidth=2.5, label="", 
         title="Cumulative Share of Output", 
         legend=:bottomright)
    
    # 第三个图形：累积就业分布
    p3 = plot(log.(n_grid), cumsum(sN[:,2]), 
         linewidth=2.5, label="", 
         title="Cumulative Share of Employment", 
         xlabel="log(employment)")
    
    # 同时显示所有图形
    savefig(p1,"figures/plant_distribution.pdf")
    savefig(p2,"figures/cumulative_output_distribution.pdf")
    savefig(p3,"figures/cumulative_employment_distribution.pdf")
end