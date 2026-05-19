#=
================================================================================
  Dissipative Dicke Model: Phase Transition via Cumulant Expansion
================================================================================

  Author : Juan Antonio Luna Gutierrez
  Contact: j.luna@edu.rptu.de
  Date   : 2026

  OVERVIEW
  --------
  This script simulates the superradiant phase transition of the driven-dissipative
  Dicke model, a paradigmatic quantum optics system in which N two-level atoms
are collectively coupled to a single bosonic cavity mode.

  Above a critical coupling strength λ_c, the system undergoes a phase transition
  from a normal phase (α = 0, γ = 0) to a superradiant phase (α ≠ 0, γ ≠ 0),
  where the cavity becomes macroscopically populated and the atomic ensemble
  develops a coherent collective excitation.

  PHYSICAL MODEL
  --------------
  The system is governed by the open Tavis-Cummings / Dicke Hamiltonian:

      H = ω_a a†a + ω_0 c†c + λ (a† + a)(c† f + f c)

  where:
    - a†, a       : creation/annihilation operators for the cavity photon mode
    - c†, c       : operators for the collective atomic pseudo-spin (Holstein-Primakoff bosons)
    - f           : Holstein-Primakoff square-root factor enforcing the spin constraint
    - ω_a         : cavity frequency
    - ω_0         : atomic transition frequency
    - λ           : atom-cavity coupling strength
    - N           : total number of atoms

  Dissipation (cavity photon loss at rate κ) is included via the Lindblad master equation:
  
      dρ/dt = -i[H, ρ] + κ (a ρ a† - ½ a†a ρ - ½ ρ a†a)

  APPROACH: SEMI-CLASSICAL EXPANSION + CUMULANT TRUNCATION
  ---------------------------------------------------------
  1. Mean-field displacement: Each operator is decomposed into a macroscopic
     mean field (scaling as √N) and a zero-mean quantum fluctuation:
         a = √N α + δa,    c = √N γ + δc

  2. Holstein-Primakoff expansion: The constraint |γ|² ≤ 1 (atomic inversion
     bounded by s = √(1 - |γ|²)) is expanded perturbatively in 1/N, capturing
     leading quantum back-action corrections.

  3. Cumulant truncation at 2nd order: Expectation values ⟨δa† δa⟩, ⟨δc† δc⟩,
     and off-diagonal correlators are tracked as dynamic variables, closing the
     hierarchy of moment equations at second order.

  4. Equations of motion are generated symbolically via QuantumCumulants.jl,
     then converted to a ModelingToolkit.jl ODESystem for efficient numerical
     integration.

  NUMERICAL STRATEGY: ADIABATIC SWEEP
  ------------------------------------
  Rather than finding steady states analytically, this script performs an
  adiabatic parameter sweep: for each value of λ, it integrates the ODE system
  to steady state, then uses that solution as the initial condition for the next
  λ value. This naturally traces the physical branch of the phase transition
  and captures the finite-N rounding of the critical point.

================================================================================
=#

using QuantumCumulants   # Symbolic cumulant hierarchy generation
using ModelingToolkit    # Symbolic ODE system assembly and simplification
using OrdinaryDiffEq     # Adaptive ODE solvers (Tsit5, Rosenbrock23)
using Plots
using Latexify

import ModelingToolkit: t_nounits as t, D_nounits as D


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Hilbert Space, Operators, and Hamiltonian
# ─────────────────────────────────────────────────────────────────────────────

# Define two independent bosonic Fock spaces:
#   :cavity — the photon mode
#   :atoms  — the collective Holstein-Primakoff bosons for the atomic pseudo-spin
ha = FockSpace(:cavity)
hc = FockSpace(:atoms)
h  = ha ⊗ hc   # Full composite Hilbert space

# Quantum fluctuation operators (zero-mean by construction after displacement)
da = Destroy(h, :da, 1)   # Cavity fluctuation annihilation operator
dc = Destroy(h, :dc, 2)   # Atomic pseudo-spin fluctuation annihilation operator

# Commuting numbers: mean fields (complex amplitudes) and physical parameters
@cnumbers α γ ω_a ω_0 κ λ s N

# ── Mean-field displacement ──────────────────────────────────────────────────
# Decompose the full field operators as:
#   a_full = √N · α  +  δa
#   c_full = √N · γ  +  δc
# where α, γ ~ O(1) are the macroscopic order parameters and δa, δc are
# the quantum fluctuations (with ⟨δa⟩ = ⟨δc⟩ = 0 by definition).
# The √N prefactor ensures the energy is extensive (scales as N).
a_full = sqrt(N) * α + da
c_full = sqrt(N) * γ + dc

# ── Holstein-Primakoff expansion of the atomic inversion ────────────────────
# The collective spin constraint requires |γ|² ≤ 1. The inversion operator
# s_z = s - c†c/N is expanded perturbatively in 1/√N:
#
#   f ≈ s - (1/2s) [ (γ* δc + γ δc†) / √N  +  δc† δc / N ]
#
# where s = √(1 - |γ|²). Using exact rational 1//2 prevents floating-point
# errors from propagating through the symbolic algebra engine (Symbolics.jl).
f_exp = s - (1//(2*s)) * ( (conj(γ)*dc + γ*dc') / sqrt(N) + (dc'*dc) / N )

# ── Effective Hamiltonian ────────────────────────────────────────────────────
# We insert the displaced + expanded operators into the Tavis-Cummings Hamiltonian.
# Keeping the expression un-expanded at this stage is intentional: Symbolics.jl
# evaluates lazy polynomial trees ~10x faster than fully expanded forms during
# the cumulant equation generation step.
H_shift = ω_a * (a_full' * a_full) +
          ω_0 * (c_full' * c_full) +
          λ * (a_full' + a_full) * (c_full' * f_exp + f_exp * c_full)


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Generate Cumulant Equations of Motion
# ─────────────────────────────────────────────────────────────────────────────

# Lindblad jump operator: cavity photon loss (rate κ).
# Only the cavity field a_full contributes to dissipation; the atoms are assumed
# to have negligible spontaneous emission on the timescale of interest.
J_sym = [a_full]

# Observable list: mean fields (1st order) and correlators (2nd order).
# Together these form a closed hierarchy at 2nd-order cumulant truncation.
ops = [
    da,        # ⟨δa⟩  — cavity fluctuation mean (enforced = 0)
    dc,        # ⟨δc⟩  — atomic fluctuation mean (enforced = 0)
    da'*da,    # ⟨δa†δa⟩ — cavity photon number fluctuation
    dc'*dc,    # ⟨δc†δc⟩ — atomic excitation number fluctuation
    da*da,     # ⟨δa δa⟩ — cavity squeezing correlator
    dc*dc,     # ⟨δc δc⟩ — atomic squeezing correlator
    da'*dc,    # ⟨δa†δc⟩ — cross-mode coherence (cavity–atom)
    da*dc      # ⟨δa δc⟩ — cross-mode anomalous correlator
]

# Generate Heisenberg-Langevin equations of motion via the Master Equation,
# truncated at 2nd cumulant order (factorizing all 3rd-order moments).
# `complete` adds any additional operator averages needed to close the hierarchy.
eqs_raw      = meanfield(ops, H_shift, J_sym; rates=[κ], order=2, simplify=false)
eqs_complete = complete(eqs_raw, simplify=false)


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Enforce ⟨δa⟩ = ⟨δc⟩ = 0 (Mean-Field Self-Consistency)
# ─────────────────────────────────────────────────────────────────────────────

# By the displacement construction, the fluctuations have exactly zero mean.
# We substitute this constraint explicitly before extracting the mean-field ODEs.
zero_dict = Dict(
    average(da)  => 0,
    average(da') => 0,
    average(dc)  => 0,
    average(dc') => 0,
    conj(N)      => N    # N is real-valued
)

# Extract right-hand sides for d⟨a_full⟩/dt and d⟨c_full⟩/dt, then divide
# by √N to isolate dα/dt and dγ/dt (since ⟨a_full⟩ = √N α).
rhs_α_raw = substitute(eqs_complete[1].rhs, zero_dict)
rhs_γ_raw = substitute(eqs_complete[2].rhs, zero_dict)

d_alpha_dt_qc = rhs_α_raw / sqrt(N)
d_gamma_dt_qc = rhs_γ_raw / sqrt(N)

# Apply the same zero-mean substitution to all 2nd-order fluctuation equations
fluct_eqs_clean = [substitute(eq, zero_dict) for eq in eqs_complete[3:end]]
n_fluct = length(fluct_eqs_clean)


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: MTK State Variables and Physical Safety Constraints
# ─────────────────────────────────────────────────────────────────────────────

# We split every complex ODE into two real ODEs (real and imaginary parts).
# ModelingToolkit requires explicit real-valued state variables.
@variables α_re(t) α_im(t) γ_re(t) γ_im(t)
@parameters ω_a_p ω_0_p κ_p λ_p N_p

# Reconstruct complex scalars for substitution
α_c = α_re + im * α_im   # α  (mean field)
α_s = α_re - im * α_im   # α* (complex conjugate)
γ_c = γ_re + im * γ_im   # γ  (atomic mean field)
γ_s = γ_re - im * γ_im   # γ* (complex conjugate)

# ── Physical safety clamp on the atomic inversion ───────────────────────────
# The Holstein-Primakoff inversion reads s = √(1 - |γ|²), which requires |γ|² < 1.
# Numerical noise near the phase transition can momentarily push |γ|² ≥ 1, causing
# a domain error (√ of negative). We clamp |γ|² to 0.9999 as a physical ceiling
# (fully inverted ensemble) to keep the solver in the valid domain.
γ_sq      = γ_re^2 + γ_im^2
γ_sq_safe = ifelse(γ_sq > 0.9999, 0.9999, γ_sq)
s_expr    = sqrt(1 - γ_sq_safe)

# Create one real and one imaginary MTK variable per fluctuation correlator
fluct_re = [only(@variables $(Symbol("f_re_$i"))(t)) for i in 1:n_fluct]
fluct_im = [only(@variables $(Symbol("f_im_$i"))(t)) for i in 1:n_fluct]


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Substitution Map — QuantumCumulants Symbols → MTK Variables
# ─────────────────────────────────────────────────────────────────────────────

# This dictionary maps every symbolic quantity from QuantumCumulants.jl
# (abstract operator averages and parameters) to concrete MTK variables
# that the ODE solver can integrate numerically.
mtk_map = Dict{Any, Any}(
    α => α_c,  conj(α) => α_s,
    γ => γ_c,  conj(γ) => γ_s,
    s => s_expr, conj(s) => s_expr,   # s is real
    ω_a => ω_a_p, conj(ω_a) => ω_a_p,
    ω_0 => ω_0_p, conj(ω_0) => ω_0_p,
    κ   => κ_p,   conj(κ)   => κ_p,
    λ   => λ_p,   conj(λ)   => λ_p,
    N   => N_p,   conj(N)   => N_p
)

# Map each fluctuation correlator ⟨O⟩ and its conjugate ⟨O†⟩ to real/imaginary MTK variables.
# We use !isequal (not !=) because != on symbolic expressions produces an un-compilable
# symbolic inequality; !isequal gives a plain Bool at compile time.
for i in 1:n_fluct
    lhs_expr = fluct_eqs_clean[i].lhs
    avg_term = (Symbolics.istree(lhs_expr) && Symbolics.operation(lhs_expr) isa Differential) ?
               lhs_expr.arguments[1] : lhs_expr

    O   = avg_term.arguments[1]
    f_c = fluct_re[i] + im * fluct_im[i]
    f_s = fluct_re[i] - im * fluct_im[i]

    mtk_map[avg_term]       = f_c
    mtk_map[conj(avg_term)] = f_s

    # Only map ⟨O†⟩ separately if it is genuinely distinct from ⟨O⟩
    # (i.e., O is not self-adjoint, e.g. da†da is self-adjoint, da†dc is not).
    adj_avg = average(adjoint(O))
    if !isequal(adj_avg, avg_term)
        mtk_map[adj_avg] = f_s
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Complex → Real ODE Splitting
# ─────────────────────────────────────────────────────────────────────────────

# Each complex-valued cumulant equation dz/dt = f(z) is split into two real ODEs:
#   d(Re z)/dt = Re(f)    and    d(Im z)/dt = Im(f)
# `expand` forces full polynomial expansion before real/imag extraction, which is
# necessary because Symbolics.jl cannot split unexpanded symbolic products correctly.
function split_complex_eq(var_re, var_im, rhs_qc)
    rhs_sub = substitute(rhs_qc, mtk_map)
    rhs_exp = expand(rhs_sub)
    return [D(var_re) ~ real(rhs_exp), D(var_im) ~ imag(rhs_exp)]
end

# Mean-field ODEs (α and γ, 4 real equations total)
mf_eqs = vcat(
    split_complex_eq(α_re, α_im, d_alpha_dt_qc),
    split_complex_eq(γ_re, γ_im, d_gamma_dt_qc)
)

# Fluctuation ODEs (one pair per correlator)
fluct_eqs_mtk = Equation[]
for i in 1:n_fluct
    for eq_item in split_complex_eq(fluct_re[i], fluct_im[i], fluct_eqs_clean[i].rhs)
        push!(fluct_eqs_mtk, eq_item isa AbstractArray ? eq_item[1] : eq_item)
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Assemble and Structurally Simplify the ODE System
# ─────────────────────────────────────────────────────────────────────────────

all_eqs_raw = vcat(mf_eqs, fluct_eqs_mtk)

# Flatten to a strict Vector{Equation} (split_complex_eq may return nested arrays)
strict_eqs = Equation[]
for item in all_eqs_raw
    push!(strict_eqs, item isa AbstractArray ? item[1] : item)
end

# `structural_simplify` performs alias elimination, index reduction, and
# tearing, then generates optimized Julia code for the ODE right-hand side.
# This step can take several minutes on first run (it JIT-compiles a large
# symbolic polynomial into native machine code via ModelingToolkit's codegen).
@named sys = ODESystem(strict_eqs, t)
sys_s = structural_simplify(sys)


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Adiabatic Sweep Across the Phase Transition
# ─────────────────────────────────────────────────────────────────────────────

# Physical parameters (resonance condition: ω_a = ω_0 = 1, weak decay κ = 0.1)
# Classical critical coupling in the thermodynamic limit for the dissipative
# Dicke model:  λ_c = √(ω_a² + κ²/4) · √ω_0 / 2
#             = √(1 + 0.0025) / 2  ≈ 0.5006
λ_values      = range(0.0, 2.0, length=80)
N_values      = [10.0, 50.0, 500.0, 1000.0]
λ_c_classical = sqrt(1.0^2 + 0.1^2 / 4) / 2

p1 = plot(xlabel="λ / λ_c", ylabel="|α|²",       title="Cavity Mean Field",       legend=:topleft)
p2 = plot(xlabel="λ / λ_c", ylabel="|γ|²",       title="Atomic Mean Field",       legend=:topleft)
p3 = plot(xlabel="λ / λ_c", ylabel="⟨δa†δa⟩/N", title="Photon Number Fluctuations", legend=:topleft)

for N_val in N_values
    println("\n── Adiabatic sweep: N = $N_val ──────────────────────────────────")

    # Pre-allocate observable storage
    α2_ss          = Float64[]
    γ2_ss          = Float64[]
    δa_dag_δa_ss   = Float64[]

    # ── Initial conditions (normal phase: all fields zero) ───────────────────
    # We start from zero and let the noise seed (below) break the symmetry
    # spontaneously once λ exceeds λ_c.
    u_current = Dict{Any, Float64}(
        α_re => 0.0, α_im => 0.0,
        γ_re => 0.0, γ_im => 0.0
    )
    for i in 1:n_fluct
        u_current[fluct_re[i]] = 0.0
        u_current[fluct_im[i]] = 0.0
    end

    for λ_val in λ_values
        params = Dict(
            ω_a_p => 1.0,
            ω_0_p => 1.0,
            κ_p   => 0.1,
            λ_p   => λ_val,
            N_p   => N_val
        )

        # ── Spontaneous symmetry breaking seed ───────────────────────────────
        # The normal phase (α = 0) is a mathematical saddle point above λ_c:
        # the ODE will remain stuck at zero unless perturbed. We inject a small
        # random complex displacement to select a symmetry-broken direction.
        # The random phase θ ensures no preferred direction is imposed by the code.
        noise_amp = 0.05
        if (u_current[α_re]^2 + u_current[α_im]^2) < noise_amp^2
            θ = 2π * rand()
            u_current[α_re] = noise_amp * cos(θ)
            u_current[α_im] = noise_amp * sin(θ)
        end

        # ── Time integration to steady state ─────────────────────────────────
        # T = 300 >> 1/κ = 10, so the system reliably reaches steady state.
        # AutoTsit5(Rosenbrock23): non-stiff solver (Tsit5) with automatic
        # fallback to a stiff solver (Rosenbrock23) near the phase transition,
        # where the Jacobian becomes ill-conditioned.
        prob = ODEProblem(sys_s, merge(u_current, params), (0.0, 300.0))
        sol  = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-8, abstol=1e-8)

        # ── Adiabatic update ─────────────────────────────────────────────────
        # Use the steady-state solution as the starting point for the next λ.
        # This is the core of the adiabatic sweep: we continuously follow the
        # physical branch of the solution rather than re-initializing from zero.
        for var in unknowns(sys_s)
            u_current[var] = sol[var][end]
        end

        # ── Fluctuation safety reset ─────────────────────────────────────────
        # Right at the phase transition, the fluctuations can diverge briefly
        # (critical slowing down). If any fluctuation exceeds the mean-field
        # scale √N (unphysical), we zero them out to avoid corrupting the sweep.
        if any(abs(u_current[fluct_re[i]]) > sqrt(N_val) for i in 1:n_fluct)
            for i in 1:n_fluct
                u_current[fluct_re[i]] = 0.0
                u_current[fluct_im[i]] = 0.0
            end
        end

        # Record steady-state observables
        push!(α2_ss,        sol[α_re][end]^2 + sol[α_im][end]^2)
        push!(γ2_ss,        sol[γ_re][end]^2 + sol[γ_im][end]^2)
        push!(δa_dag_δa_ss, sol[fluct_re[1]][end])
    end

    # Plot normalized by N for finite-size comparison
    plot!(p1, λ_values, α2_ss,                  lw=2, label="N = $(Int(N_val))")
    plot!(p2, λ_values, γ2_ss,                  lw=2, label="N = $(Int(N_val))")
    plot!(p3, λ_values, δa_dag_δa_ss ./ N_val,  lw=2, label="N = $(Int(N_val))")
end


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Annotate and Display Results
# ─────────────────────────────────────────────────────────────────────────────

# Mark the thermodynamic-limit critical coupling λ_c on all panels
for p in (p1, p2, p3)
    vline!(p, [λ_c_classical], ls=:dash, color=:red, lw=1.5,
           label="λ_c (N → ∞)")
end

fig = plot(p1, p2, p3,
           layout = (1, 3),
           size   = (1400, 450),
           margin = 5Plots.mm)

display(fig)
savefig(fig, "dicke_phase_transition.png")
println("\nFigure saved to dicke_phase_transition.png")
