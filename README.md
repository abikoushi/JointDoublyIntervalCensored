# JointDoublyIntervalCensored

An implementation of

Abe & Sakumura (2025) "Empirical Joint Estimation for Infection and Incubation Period under the Doubly Interval-Censored with Right Truncated Observation." (in preparation)

## Structure of this repository

- Gibbs sampler: `jl/JointTools_Gibbs.jl`
- VB and EM algortithm: `jl/JointTools.jl`
- Run simulation study in Section 3: `jl/run_sim_joint.jl`
- Run case study in Section 4: `jl/dataanalysis.jl`

