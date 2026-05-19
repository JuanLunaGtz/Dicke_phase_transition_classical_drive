# Dicke_phase_transition_classical_drive
This repository contains a Julia simulation of the superradiant phase transition in the driven-dissipative Dicke model — a canonical quantum optics system in which NN N two-level emitters are collectively coupled to a single bosonic cavity mode. 

The simulation captures the full crossover from the normal phase (α=0) to the superradiant phase (α≠0) as the atom-cavity coupling λ is swept across the critical point λc​, and tracks how quantum fluctuations diverge near the transition for finite N. 

This code directly supports my master's thesis research on how non-classical (squeezed) light modifies this phase transition.

# Physical Model
The system is governed by the open Tavis-Cummings / Dicke Hamiltonian:

$$\hat{H} = \omega_a \hat{a}^\dagger \hat{a} + \omega_0 \hat{c}^\dagger \hat{c} + \lambda (\hat{a}^\dagger + \hat{a})(\hat{c}^\dagger f(\hat{n}) + f(\hat{n}) \hat{c})$$

with cavity photon loss described by the Lindblad master equation:

$$\frac{d\hat{\rho}}{dt} = -i[\hat{H},\hat{\rho}] + \kappa \left( \hat{a} \hat{\rho} \hat{a}^\dagger - \tfrac{1}{2} \hat{a}^\dagger \hat{a} \hat{\rho} - \tfrac{1}{2} \hat{\rho} \hat{a}^\dagger \hat{a} \right)$$

The classical critical coupling in the thermodynamic limit is:
$$\lambda_c = \sqrt{\frac{(\omega_a^2 + \kappa^2)\omega_0}{4\omega_a}}$$


Results
The simulation produces a three-panel figure showing, as a function of λ\lambda
λ:

∣α∣2|\alpha|^2
∣α∣2 — cavity mean-field occupation (order parameter)
∣γ∣2|\gamma|^2
∣γ∣2 — atomic mean-field occupation
⟨δa†δa⟩/N\langle \delta a^\dagger \delta a \rangle / N
⟨δa†δa⟩/N — normalized photon number fluctuations

Each quantity is shown for N∈{10,50,500,1000}N \in \{10, 50, 500, 1000\}
N∈{10,50,500,1000}, illustrating the convergence to the thermodynamic limit and the finite-NN
N rounding of the critical point.​​
