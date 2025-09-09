# Monte Carlo Simulation Project

Clean, runnable R code for three core Monte Carlo topics:

1) **Rejection & Importance Sampling** — simulate from a heavy-tailed target when inverse CDF is unavailable; KS goodness-of-fit; importance weights and stability  
2) **Markov Chains / MCMC** — build a valid transition matrix, compute stationary distribution, compare empirical vs theoretical long-run behaviour  
3) **Random-walk PDE solver (Laplace/Poisson)** — estimate harmonic function values via random walks; error vs number of walks in 1D and 2D

> This repo is a tidy re-implementation of my undergraduate project (module mark: 68%). Code is organised for readability and reproducibility. See [docs/overview.md](docs/overview.md) for context.

## Quickstart
- Open `src/monte_carlo_project.R` and run **Part 1 → Part 2 → Part 3** in order  
- Plots print inline; you can direct outputs to `figures/` / `results/`  
- Base R only (no extra packages needed)

## Structure
What you’ll see when you run it
- **Part 1:** Histogram with theoretical PDF overlay; a KS table across sample sizes; an IS estimate for E|X| with a short note on weight stability.
- **Part 2:** The stationary distribution from the transition matrix and the empirical state frequencies from a long chain, shown side by side.
- **Part 3:** A small error table (1D/2D) and a 1D log–log convergence plot illustrating error decreasing with the number of walks.


## Reproduce
- R ≥ 4.1  
- Seeds set via `set.seed()` where relevant  
- Optionally append `sessionInfo()` to record environment


## Challenges & lessons
- **Proposal vs acceptance.** Uniform on [μ, μ+10k] with c=10α is easy to code but not super efficient. Good reminder that simplicity ≠ speed.
- **KS at small n.** With tiny samples the KS test jumps around; scaling n up made the fit settle.
- **Importance weights.** Badly matched proposals blow up IS variance. Watching the weight distribution (not just the estimate) was key.
- **Markov chain sanity.** Added a small self-loop to keep the chain irreducible/aperiodic; long runs then matched the stationary vector. Short runs can mislead.
- **Compute limits (real talk).** Full grids / very large N **stalled or didn’t finish** on my laptop. I reduced settings and labelled them. Behaviour still showed the expected trend (~O(N⁻¹/²)), just with shorter runtimes.


## Run settings & performance notes
- **What I actually ran:** 1D `n=31`, N ∈ {200, 500, 1000}; 2D `15×15`, N ∈ {1000, 2000}.  
- **What I tried but cut:** 1D `n=51`, 2D `21×21`, and N up to 10k — **too slow on my machine**.  
- **Tip:** If you’ve got more compute, feel free to bump N/grid; you’ll get smaller errors but expect longer runs.


## License
MIT
