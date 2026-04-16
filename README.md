# JointDoublyIntervalCensored

An implementation of

Abe, K., & Sakumura, T. (2026). Empirical Joint Estimation for Incubation Period and Infection under the Doubly Interval-Censored with Right Truncated Observation. Journal of Computational and Graphical Statistics, 1–13. https://doi.org/10.1080/10618600.2026.2652940

## Structure of this repository

- Gibbs sampler: `jl/JointTools_Gibbs.jl`
- VB and EM algortithm: `jl/JointTools.jl`
- Run simulation study in Section 3: `jl/run_sim_joint.jl`
- Run case study in Section 4: `jl/dataanalysis.jl`

