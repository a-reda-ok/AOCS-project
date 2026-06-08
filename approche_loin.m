%% =========================================================================
%  approche_loin.m
%  Far-range ballistic V-bar hop: 100 km --> 1 km
%
%  Context: Orbital rendezvous GNC — SCAO / AOCS design exercise
%           in collaboration with Thales Alenia Space.
%
%  Strategy: Two-impulse, open-loop (ballistic) transfer along the V-bar
%            using Hill-Clohessy-Wiltshire (HCW) relative motion equations.
%            The chaser starts 100 km behind the target on the V-bar and
%            reaches 1 km behind it after a user-chosen number of orbits.
%
%  Reference frame: LVLH (Local Vertical Local Horizontal)
%    x  — R-bar (radial, positive outward)
%    y  — V-bar (along-track, positive in velocity direction)
%    z  — H-bar (cross-track, positive in angular momentum direction)
%
%  Key assumptions:
%    - Circular reference orbit (HCW validity)
%    - Impulsive manoeuvres (instantaneous delta-V)
%    - Simplified V-bar hop formula: vy_drift = -3 * dV_y
%      (valid for multi-orbit transfers where x_offset << orbit radius)
%
%  Author : [Your name / team]
%  Date   : 2024
%% =========================================================================

clear; clc; close all;

%% ---- 1. MISSION CONSTANTS ----------------------------------------------

mu       = 3.986e14;   % Earth gravitational parameter [m^3/s^2]
r_earth  = 6371e3;     % Earth mean radius             [m]
alt_cible = 770e3;     % Target orbit altitude         [m]

% Derived orbital parameters
a       = r_earth + alt_cible;  % Semi-major axis of reference orbit [m]
n       = sqrt(mu / a^3);       % Mean motion                        [rad/s]
T_orbit = 2 * pi / n;           % Orbital period                     [s]

fprintf('=== Far-Range V-bar Hop Simulation (100 km -> 1 km) ===\n');
fprintf('Orbital period : %.1f min\n\n', T_orbit / 60);

%% ---- 2. HCW STATE-SPACE MODEL ------------------------------------------
%
%  State vector: S = [x; y; z; x_dot; y_dot; z_dot]
%
%  HCW equations (unforced):
%    x_ddot = 3*n^2*x + 2*n*y_dot
%    y_ddot =          -2*n*x_dot
%    z_ddot = -n^2*z
%
%  Written as a first-order ODE for ode45:

ode_hcw = @(t, S) [ S(4);                       % dx/dt   = vx
                    S(5);                       % dy/dt   = vy
                    S(6);                       % dz/dt   = vz
                    3*n^2*S(1) + 2*n*S(5);     % dvx/dt  (HCW)
                   -2*n*S(4);                   % dvy/dt  (HCW)
                   -n^2*S(3) ];                 % dvz/dt  (HCW)

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

%% ---- 3. TRANSFER DESIGN (TWO-IMPULSE V-BAR HOP) -----------------------
%
%  The V-bar hop exploits a well-known HCW property:
%    A single tangential impulse dV_y applied at a point on the V-bar
%    produces a mean drift velocity:
%         vy_drift = -3 * dV_y
%    (positive dV_y -> retrograde relative to target -> chaser drifts forward)
%
%  Manoeuvre plan:
%    t = 0          : Impulse 1 (dV_1 < 0, retrograde) initiates drift.
%    t = t_transfer : Impulse 2 (dV_2 = -dV_1, prograde) stops drift at target.
%
%  The transfer duration is chosen by the user (here: 2 orbital periods).

dist_initiale = -100e3;   % Initial y-position (behind target) [m]
dist_finale   = -1e2;     % Final   y-position (still behind)  [m]
delta_y_total = dist_finale - dist_initiale;  % Net y-displacement required [m]  (+99 000 m)

t_transfert = 2 * T_orbit;   % Transfer duration: 2 orbits [s]  <-- USER PARAMETER

% Required mean drift velocity along V-bar
vy_drift_req = delta_y_total / t_transfert;   % [m/s]

% From vy_drift = -3 * dV_y  =>  dV_y = -vy_drift / 3
dV_1_y = -vy_drift_req / 3;   % First  impulse (retrograde) [m/s]
dV_2_y = -dV_1_y;             % Second impulse (prograde)   [m/s]

fprintf('Transfer duration        : %.1f h\n',   t_transfert / 3600);
fprintf('Required drift velocity  : %.2f m/s\n', vy_drift_req);
fprintf('Impulse dV_1 (retrograde): %.2f m/s\n', dV_1_y);
fprintf('Impulse dV_2 (prograde)  : %.2f m/s\n', dV_2_y);
fprintf('Total |Delta-V|          : %.2f m/s\n', 2 * abs(dV_1_y));

%% ---- 4. NUMERICAL SIMULATION (TWO PHASES) ------------------------------

% --- Phase 1: Drift phase (0 to t_transfer) ---
% Initial state: chaser at y = -100 km, all velocities zero except dV_1_y
t_span1 = [0, t_transfert];
S0_phase1 = [0; dist_initiale; 0;   % x, y, z  [m]
             0; dV_1_y;       0];   % vx, vy, vz [m/s]

[t1, S1] = ode45(ode_hcw, t_span1, S0_phase1, options);

% --- Phase 2: Station-keeping at 1 km (one extra orbit for verification) ---
% Apply dV_2 at the end of the drift phase, then coast
t_span2 = [t_transfert, t_transfert + T_orbit];
S0_phase2      = S1(end, :)';          % Carry over final state from phase 1
S0_phase2(5)   = S0_phase2(5) + dV_2_y;   % Apply second impulse (vy only)

[t2, S2] = ode45(ode_hcw, t_span2, S0_phase2, options);

% Concatenate results for plotting
t_total = [t1; t2];
S_total = [S1; S2];

%% ---- 5. VISUALISATION --------------------------------------------------

figure('Name', 'Far-Range V-bar Hop (100 km -> 1 km)', ...
       'NumberTitle', 'off', 'Color', 'w');

% --- 5a. 3-D LVLH trajectory (Y, X, Z axes for intuitive V-bar view) ---
subplot(2, 1, 1);
plot3(S_total(:,2)/1e3, S_total(:,1)/1e3, S_total(:,3)/1e3, ...
      'b-', 'LineWidth', 1.5);
hold on;

% Key event markers
plot3(S1(1,   2)/1e3, S1(1,   1)/1e3, S1(1,   3)/1e3, ...
      'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');   % Departure
plot3(S1(end, 2)/1e3, S1(end, 1)/1e3, S1(end, 3)/1e3, ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');   % dV_2 application
plot3(S2(end, 2)/1e3, S2(end, 1)/1e3, S2(end, 3)/1e3, ...
      'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');   % End of simulation

title('3-D LVLH Trajectory');
xlabel('V-bar / Y axis (km)');
ylabel('R-bar / X axis (km)');
zlabel('H-bar / Z axis (km)');
legend('Trajectory', 'Departure (-100 km)', ...
       'dV_2 applied (-1 km)', 'End of simulation', ...
       'Location', 'best');
grid on; axis equal; view(3);
hold off;

% --- 5b. V-bar position vs. time ---
subplot(2, 1, 2);
plot(t_total/3600, S_total(:,2)/1e3, 'r-', 'LineWidth', 1.5);
hold on;

% Mark the moment of second impulse
idx_dV2 = length(t1);
plot(t_total(idx_dV2)/3600, S_total(idx_dV2, 2)/1e3, ...
     'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Target line
yline(dist_finale/1e3, 'k:', 'LineWidth', 1.2, ...
      'Label', 'Target: -1 km');

title('V-bar Position vs. Time');
xlabel('Time (hours)');
ylabel('Y position (km)');
legend('V-bar distance', 'dV_2 applied', 'Location', 'best');
grid on;
hold off;
