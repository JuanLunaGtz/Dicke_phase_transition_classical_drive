# Dicke_phase_transition_classical_drive
This repository contains a Julia simulation of the superradiant phase transition in the driven-dissipative Dicke model — a canonical quantum optics system in which N two-level emitters are collectively coupled to a single bosonic cavity mode. 

The simulation captures the full crossover from the normal phase (α=0) to the superradiant phase (α≠0) as the atom-cavity coupling λ is swept across the critical point λc​, and tracks how quantum fluctuations diverge near the transition for finite N. 

This code directly supports my master's thesis research on how non-classical (squeezed) light modifies this phase transition.

# Physical Model
The system is governed by the open Tavis-Cummings / Dicke Hamiltonian:

$$\hat{H} = \omega_a \hat{a}^\dagger \hat{a} + \omega_0 \hat{c}^\dagger \hat{c} + \lambda (\hat{a}^\dagger + \hat{a})(\hat{c}^\dagger f(\hat{n}) + f(\hat{n}) \hat{c})$$

with cavity photon loss described by the Lindblad master equation:

$$\frac{d\hat{\rho}}{dt} = -i[\hat{H},\hat{\rho}] + \kappa \left( \hat{a} \hat{\rho} \hat{a}^\dagger - \tfrac{1}{2} \hat{a}^\dagger \hat{a} \hat{\rho} - \tfrac{1}{2} \hat{\rho} \hat{a}^\dagger \hat{a} \right)$$

The classical critical coupling in the thermodynamic limit is:
$$\lambda_c = \sqrt{\frac{(\omega_a^2 + \kappa^2)\omega_0}{4\omega_a}}$$

## Method
 
The simulation uses a semi-classical expansion combined with 2nd-order cumulant truncation:
 
| Step | Technique | Purpose |
|------|-----------|---------|
| 1 | Mean-field displacement: $a = \sqrt{N}\alpha + \delta a$ | Separate macroscopic order parameter from quantum noise |
| 2 | Holstein-Primakoff expansion of atomic inversion | Enforce spin constraint $\|\gamma\|^2 \leq 1$ perturbatively in $1/N$ |
| 3 | Cumulant truncation at 2nd order | Close the hierarchy of moment equations |
| 4 | Symbolic EOM generation via `QuantumCumulants.jl` | Derive all equations of motion automatically from the Hamiltonian |
| 5 | Real/imaginary splitting into `ModelingToolkit.jl` ODESystem | Compile to efficient native code |
| 6 | Adiabatic sweep + stiff ODE integration | Trace the physical branch of the transition for each $N$ |


# Results
The simulation produces a three-panel figure showing, as a function of λ\lambda
λ: $|\alpha|^2$ — cavity mean-field occupation (order parameter), $|\gamma|^2$ — atomic mean-field occupation, and
$\langle \delta a^\dagger \delta a \rangle / N$ — normalized photon number fluctuations.

Each quantity is shown for $N \in \{10, 50, 500, 1000\}$, illustrating the convergence to the thermodynamic limit and the finite-N rounding of the critical point.​

## Dependencies
 
Install all dependencies from the Julia REPL:
 
```julia
using Pkg
Pkg.add([
    "QuantumCumulants",
    "ModelingToolkit",
    "OrdinaryDiffEq",
    "Plots",
    "Latexify"
])
```

## References
 
- Dicke, R.H. (1954). *Coherence in spontaneous radiation processes*. Phys. Rev. **93**, 99.
- Emary, C. & Brandes, T. (2003). *Chaos and the quantum phase transition in the Dicke model*. Phys. Rev. E **67**, 066203.
- Kirton, P. et al. (2019). *Introduction to the Dicke Model*. Adv. Quantum Technol. **2**, 1800043.
- Plankensteiner, D. et al. (2022). *QuantumCumulants.jl*. Quantum **6**, 617.​
