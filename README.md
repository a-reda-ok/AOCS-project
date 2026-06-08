# Orbital Rendezvous GNC — MATLAB Simulations

> SCAO/AOCS design exercise conducted in collaboration with **Thales Alenia Space**.  
> Focus: Guidance, Navigation & Control (GNC) for a hypothetical satellite capture mission.

---

## Mission Overview

A chaser satellite must approach and capture a target in a **770 km circular LEO orbit**.  
The approach is split into two distinct phases, each handled by a dedicated script:

| Phase | Script | Range | Method |
|---|---|---|---|
| Far-range transfer | `approche_loin.m` | 100 km → 1 km | Ballistic V-bar hop (2 impulses) |
| Close-range approach | `hcw_forced.m` | 100 m → 2 m | LQR stop-and-go |

Both scripts model relative motion in the **LVLH frame** (Local Vertical Local Horizontal) using the **Hill-Clohessy-Wiltshire (HCW)** linearised equations, valid for close, nearly-circular orbits.

---

## Reference Frame

```
        Z (H-bar, cross-track)
        |
        |
        +-------> Y (V-bar, along-track / velocity direction)
       /
      /
     X (R-bar, radial outward)

Target satellite sits at the origin.
```

> **Note:** variable naming inside the MATLAB code (`x`, `y`, `z`) follows a different internal convention. All plots are labelled with the correct physical axis names.

---

## Scripts

### `approche_loin.m` — Far-range ballistic transfer

**Strategy:** Two-impulse open-loop V-bar hop.

The chaser starts 100 km behind the target and must close to 1 km.  
A first retrograde impulse `dV_1` raises the chaser's orbit slightly, inducing a natural forward drift. After a chosen number of orbits a second prograde impulse `dV_2 = -dV_1` cancels the drift exactly at the target distance.

**Key physics used:**

For a tangential impulse `dV_y` applied on the V-bar, HCW predicts a mean drift:

```
vy_drift = -3 * dV_y
```

so the required impulse for a given drift velocity is:

```
dV_1 = -vy_drift_required / 3
```

**User-tunable parameters:**

| Parameter | Variable | Default |
|---|---|---|
| Transfer duration | `t_transfert` | 2 orbital periods |
| Final hold distance | `dist_finale` | −100 m |

**Outputs:**
- 3-D LVLH trajectory plot
- V-bar distance vs. time plot
- Console summary (transfer time, drift velocity, impulse magnitudes, total ΔV)

---

### `hcw_forced.m` — Close-range LQR stop-and-go approach

**Strategy:** Optimal (LQR) feedback control through a sequence of waypoints.

The chaser starts 100 m away on the V-bar. It stops at intermediate waypoints (50 m, 10 m) before reaching the final hold point at 2 m. Thruster force is saturated at ±2 N per axis.

A 7th ODE state accumulates the integral of `‖u_applied‖` over time, giving the total ΔV without post-processing.

**Controller:**

```
u = -K·(x - x_ref) + u_ff
```

where `K` is the LQR gain computed offline from `lqr(A, B, Q, R)` and `u_ff = -(pinv(B)·A)·x_ref` is a feedforward term that places the closed-loop equilibrium exactly at the target waypoint.

**User-tunable parameters:**

| Parameter | Variable | Default |
|---|---|---|
| Waypoints (V-bar) | `waypoints_y` | `[50, 10, 2]` m |
| Waypoints (R-bar) | `waypoints_x` | `[0, 0, 0]` m |
| Max time per segment | `t_per_segment` | 3000 s |
| Thruster saturation | `F_max_per_axis` | 2 N |
| LQR position weight | `Q` (diagonal) | see code |
| LQR control weight | `R` (diagonal) | see code |

**Outputs:**
- 3-D LVLH trajectory
- Position and velocity time histories
- Applied thruster force history (showing saturation)
- Cumulative ΔV vs. time
- Cross-track error (deviation from ideal V-bar corridor)
- Console budget: total ΔV and propellant mass (Tsiolkovsky)

---

## Requirements

- **MATLAB R2018b or later** (uses `lqr`, `ode45`, `yline`)
- Control System Toolbox (for `lqr`)

No additional toolboxes are required.

---

## Quick Start

```matlab
% Far-range phase
run('approche_loin.m')

% Close-range phase
run('hcw_forced.m')
```

---

## Results (default parameters)

**Far-range transfer** (`approche_loin.m`):
- Transfer duration: ~3.3 hours (2 orbital periods)
- Total |ΔV|: ~5.5 m/s (two equal impulses)

**Close-range approach** (`hcw_forced.m`):
- Total ΔV: ~0.23 m/s
- Propellant mass (Isp = 220 s): ~0.17 kg

---

## Repository Structure

```
.
├── approche_loin.m     # Far-range ballistic V-bar hop simulation
├── hcw_forced.m        # Close-range LQR guided approach simulation
└── README.md           # This file
```

---

## Background

This work was produced during the **EI SCAO** intensive design exercise at Thales Alenia Space, covering the full AOCS design chain for an agile Earth observation satellite including attitude control, disturbance torque estimation, actuator/sensor sizing, and (here) orbital rendezvous GNC.
