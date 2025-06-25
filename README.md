# Entanglement-accelerates-quantum-simulation
code in the manuscript Entanglement accelerates quantum simulation arXiv:2406.02379

simulation_entanglement.jl:  all the used functions 
entanglement.jl: entanglement entropy increase and Trotter error for each segment (curves in Fig 1)

PF1theory.jl, PF2theory.jl: theoretical entanglement (distanced-based bound ) for PF1 and PF2 (curves in Fig 2 a,b)

Adaptivecheckpoint.jl: minimum Trotter step in adaptive protocol， theoretical worst-case bound, theoretical average-case bound (curves in Fig 2d, Fig 2c)
The theoretical bound in Fig 2c can be obtained by inserting many checkpoints and estimating the converged number of Trotter steps. 

Empiricalsystemsize.m: empirical worst-case, average-case, and entangled-case bounds in Fig. 4












