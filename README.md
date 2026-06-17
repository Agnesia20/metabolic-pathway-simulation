import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ----------------------------
# Model definition
# ----------------------------
def network(y, t, params):
    A, B, P = y

    V1max = params['V1max']
    Km1   = params['Km1']
    Ki    = params['Ki']
    X     = params['X']
    k2    = params['k2']
    k3    = params['k3']
    k4    = params['k4']

    # v1: Michaelis-Menten with non-competitive allosteric inhibition by P
    v1 = (V1max * X) / ((Km1 + X) * (1 + P / Ki))

    # v2: A -> B  (first-order)
    v2 = k2 * A

    # v3: B -> P  (first-order)
    v3 = k3 * B

    # v4: A -> byproduct (first-order, escape pathway)
    v4 = k4 * A

    # ODEs (from stoichiometric matrix S . v)
    dA = v1 - v2 - v4
    dB = v2 - v3
    dP = v3

    return [dA, dB, dP]


# ----------------------------
# Parameters (from exam table)
# ----------------------------
params = {
    'V1max': 5.0,
    'Km1':   2.0,
    'Ki':    3.0,
    'X':     10,
    'k2':    1.0,
    'k3':    0.8,
    'k4':    0.3,
}

# ----------------------------
# Initial conditions & time span
# ----------------------------
# Assume system starts empty (no A, B, P at t=0)
init = [0.0, 0.0, 0.0]
t = np.linspace(0, 30, 300)

# ----------------------------
# Solve ODEs
# ----------------------------
sol = odeint(network, init, t, args=(params,))
A, B, P = sol[:, 0], sol[:, 1], sol[:, 2]

# ----------------------------
# Plot results
# ----------------------------
plt.figure(figsize=(8, 5.5))
plt.plot(t, A, label='[A]', linewidth=2)
plt.plot(t, B, label='[B]', linewidth=2)
plt.plot(t, P, label='[P]', linewidth=2)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.title('Kinetic Model: 4-Reaction Pathway with\nAllosteric Feedback Inhibition of v1 by P')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('simulation_q2c.png', dpi=150, bbox_inches='tight')
plt.show()

# ----------------------------
# Print steady-state values
# ----------------------------
print("Steady-state (final time point) concentrations:")
print(f"  [A] = {A[-1]:.4f}")
print(f"  [B] = {B[-1]:.4f}")
print(f"  [P] = {P[-1]:.4f}")

# Check v1 at steady state to confirm inhibition effect
P_ss = P[-1]
v1_ss = (params['V1max'] * params['X']) / \
        ((params['Km1'] + params['X']) * (1 + P_ss / params['Ki']))
v1_uninhibited = (params['V1max'] * params['X']) / \
                 (params['Km1'] + params['X'])
print(f"\nv1 at steady state (with inhibition) = {v1_ss:.4f}")
print(f"v1 without inhibition (P=0)          = {v1_uninhibited:.4f}")
print(f"Reduction due to allosteric inhibition: {(1 - v1_ss/v1_uninhibited)*100:.1f}%")
