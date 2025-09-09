# Overview

This codebase implements three Monte Carlo components in base R:

## 1) Rejection & Importance Sampling

- **Target**: heavy-tailed density \( f(x) = \frac{\alpha}{k}\left(\frac{k-\mu+x}{k}\right)^{-\alpha-1} \) on \(x \ge \mu\)  
- **Proposal**: Uniform \(U[\mu, \mu+10k]\) with bound \(f(x) \le c\,g(x)\), using \(c = 10\alpha\)  
- **Rejection sampling**: samples accepted according to \(U \le f(Y)/(c\,g(Y))\); histogram + theoretical PDF overlay; KS tests across sample sizes  
- **Importance sampling**: estimate \(E[|X|]\) via weights \(w_i = f(x_i)/g(x_i)\); notes on weight stability and variance

## 2) Markov Chains / MCMC

- Build a **5×5** transition matrix and normalise rows  
- Compute the **stationary distribution** from the left eigenvector (eigenvalue 1)  
- **Simulate** a long chain (e.g., 100k steps) and compare **empirical state frequencies** to the stationary distribution

## 3) Random-walk PDE solver (Laplace/Poisson)

- **1D**: estimate \(u(x)\) with random walks on \([0,1]\); compare to analytic \(u(x)=x\); report max error vs number of walks  
- **2D**: on \([0,1]^2\) with **Dirichlet** left/right (0/1) and **reflecting** top/bottom boundaries; report max error vs number of walks  
- Empirical error decreases with more walks (≈ \(O(N^{-1/2})\) behaviour), consistent with standard MC convergence

## What to run

Open `src/monte_carlo_project.R`, then run sequentially:
- Part 1: rejection sampling (hist + PDF + KS), importance sampling for \(E[|X|]\)  
- Part 2: transition matrix, stationary distribution, long-run simulation  
- Part 3: `mc1d()` and `mc2d()` experiments; print error tables and plot 1D log–log convergence

## Notes

- Runs are seeded for reproducibility
- Reduced grids/walk counts are included for reasonable runtimes; increase \(N\) / resolution for tighter errors
