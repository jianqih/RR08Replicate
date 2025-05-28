# Benchmark model 
ntau=3;# number of different tax(capital,labor and the total output)
# Assume constant separation rates across plants: lambda (if separation
# rate depends on s need to create a matrix rho)
cf=0; # the fuxed cost of the production
ce=1; # the entrance cost to enter the market
lambda_exit =0.10; # the exit rate to get out of the market
rho=(1-lambda_exit)/(1+R);  # the discount factor to discount the steady state welfare

p = hcat(zeros(ns), ones(ns), zeros(ns)) #300*3
hsmatrix = theo_hs*ones(1,ntau)
g = hsmatrix.*p  #hs 100*3 ; p 100*3  elementwise multiplication  (s,tau)
smatrix = S*s_grid*ones(1,ntau)
taus = 0.0 #the subsidy
tau = 0.0 #the tax
#In the benchmark case, the subsidy and the tax are equal to zero.
tauv = [-taus, 0, tau]'
mtau = ones(ns,1)*tauv
mtauo = mtau
mtauk = zeros(ns, ntau) #tax on capital
mtaun = zeros(ns, ntau) #tax on labor
#mtaue = np.zeros((ns, ntau)) #tax on energy

function We_function(w, r, smatrix, alpha, gamma, mtauk, mtaun, mtauo, rho, g, ce, cf, ns, ntau)
    # demand for capital 
    k_bar=@. (alpha./(r.*(1+mtauk))).^((1-gamma)/(1-gamma-alpha)).*(gamma./(w.*(1+mtaun))).^(gamma/(1-gamma-alpha)).*(smatrix.*(1-mtauo)).^(1/(1-alpha-gamma))
    # demand for labor
    n_bar=@. (1+mtauk)*r*gamma./((1+mtaun)*w*alpha).*k_bar;

    # substitute kabr, nbar in pi and compute W(s,kbar,s,th)
    pi_bar=@. (1-mtauo).*smatrix.*k_bar.^alpha.*n_bar.^gamma-(1+mtaun).*w.*n_bar-(1+mtauk).*r.*k_bar-cf
    # Calculate W(s,tau)
    # initial guess of W(s,tau) value function
    W = pi_bar./(1-rho)
    xbar = zeros(ns,ntau)
    for i in 1:ns 
        for j in 1:ntau
            if W[i,j] >=0 
                xbar[i,j] = 1
            else
                xbar[i,j] = 0
            end
        end
    end
    # compute expected value function making a draw (s,theta) from G(s,theta)
    We = sum(W.*g.*xbar) - ce
    return We, k_bar, n_bar, pi_bar,W, xbar
end

# bisection on equilibrium wage
## initial guess
global w_0 = 1.0 
global Output0=We_function(w_0,r, smatrix, alpha, gamma, mtauk, mtaun, mtauo, rho, g, ce, cf, ns, ntau)
global We0=Output0[1]
print(We0)
# the output sequence of the function We 
# We,k_bar,n_bar,pi_bar,W,xbar,
global w_1=2;
global Output1=We_function(w_1,r, smatrix, alpha, gamma, mtauk, mtaun, mtauo, rho, g, ce, cf, ns, ntau)
global We1=Output1[1]

# the actual Bisection part
global converged = 0
tol2 = 0.0001
maxit2 = 100
global w_0 = 1.0
global w_1 = 2.0
global Output = zeros(6)
global w = 1.0

# w = (w_0 + w_1) / 2
for i in 1:maxit2
    global w = (w_0 + w_1) / 2
    global Output = We_function(w, r, smatrix, alpha, gamma, mtauk, mtaun, mtauo, rho, g, ce, cf, ns, ntau)
    global We=Output[1]
    if abs(We) < tol2
        global converged = 1
        println("wage has converged in ", i)
        println("wage is ", w)
        break
    else
        if We * We1 > 0
            global w_1 = w
        else
            global w_0 = w
        end
    end
    if i == maxit2
        println("Warning: Zero profit condition not satisfied")
    end
end

global We=Output[1]
global k_bar=Output[2] # optimal capital
global n_bar=Output[3] # optimal labor
global pi_bar=Output[4] # optimal profit
global W=Output[5] # optimal value function
global xbar=Output[6] # indicator function

# aggregate output
global muhat = 1 ./lambda_exit .* (xbar .* g)
global N_hat = sum(n_bar .* muhat)
global E = 1/ N_hat
global mu = muhat * E
global mu_s = sum(mu, dims=2)

global Y = sum(smatrix .* k_bar.^alpha .* n_bar.^gamma .* mu)
global K = sum(k_bar .* mu)
global KY = K / Y
global Kbe = K

global A = Y / (N_hat * E) / (K / (N_hat * E)) ^ alpha
global M = sum(mu)
global dgp = @. -mtau*smatrix*(k_bar)^alpha* n_bar^gamma*mu- mtauk*r*k_bar*mu- mtaun*w*n_bar*mu
global sgdp = sum(dgp)/Y

# Average employment per plant (AEPP)
global normalization_emp = n_bar[1, 3]
# let the optimal employee number of the lowest productivity firm be 1 
global rnbar = n_bar / normalization_emp
global AEPP = 1/ (M* normalization_emp)

global relmus = sum(mu, dims=2)
global sK = (k_bar .* mu)/K
global sN = (n_bar .* mu)/(N_hat*E)
global sY = (smatrix .* mu)/Y

# relative price of capital
# rK = alpha * (sK ./ sY) ./ (sK ./ sL)


