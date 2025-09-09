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

src/
monte_carlo_project.R # Part 1: rejection/importance; Part 2: Markov chain; Part 3: PDE random walks
docs/
overview.md # short project context + results snapshots
notebooks/ # (optional) demo notebooks if added later
figures/ # generated plots (empty placeholder)
results/ # generated tables/files (empty placeholder)
data/ # (empty placeholder — no private data committed)


## Reproduce
- R ≥ 4.1  
- Seeds set via `set.seed()` where relevant  
- Optionally append `sessionInfo()` to record environment

## License
MIT
