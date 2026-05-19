# Dicke_phase_transition_classical_drive
This repository contains a Julia simulation of the superradiant phase transition in the driven-dissipative Dicke model — a canonical quantum optics system in which NN N two-level emitters are collectively coupled to a single bosonic cavity mode. 

The simulation captures the full crossover from the normal phase (α=0) to the superradiant phase (α≠0) as the atom-cavity coupling λ is swept across the critical point λc​, and tracks how quantum fluctuations diverge near the transition for finite N. 

This code directly supports my master's thesis research on how non-classical (squeezed) light modifies this phase transition.

# Physical Model
The system is governed by the open Tavis-Cummings / Dicke Hamiltonian:

$$H = \omega_a \, a^\dagger a + \omega_0 \, c^\dagger c + \lambda (a^\dagger + a)(c^\dagger f + f \, c)$$
