%% =========================================================================
%  hcw_forced.m
%  Close-range LQR-guided stop-and-go approach: 100 m --> 2 m
%
%  Context: Orbital rendezvous GNC — SCAO / AOCS design exercise
%           in collaboration with Thales Alenia Space.
%
%  Strategy: Optimal (LQR) feedback control driving the chaser through a
%            sequence of waypoints along the V-bar (Stop-and-Go profile).
%            The simulation uses Hill-Clohessy-Wiltshire (HCW) equations
%            as the plant model inside an ode45 integration loop.
%
%  Reference frame: LVLH (Local Vertical Local Horizontal)
%    x  — R-bar (radial)
%    y  — V-bar (along-track)   <-- approach axis
%    z  — H-bar (cross-track)
%
%  NOTE: Variable naming inside the code follows MATLAB convention
%        (x/y/z = first/second/third state). The axes are relabelled
%        correctly in all plots. See comment in Section 3.
%
%  Features:
%    - LQR optimal controller (Q/R tunable)
%    - Thruster saturation (+/- F_max_per_axis per axis)
%    - Cumulative Delta-V tracking via a 7th ODE state
%    - Propellant mass estimate via Tsiolkovsky equation
%
%  Author : [Your name / team]
%  Date   : 2024
%% =========================================================================

clear; clc; close all;

%% ---- 1. SYSTEM CONSTANTS -----------------------------------------------

mu        = 3.986e14;   % Earth gravitational parameter [m^3/s^2]
r_earth   = 6371e3;     % Earth mean radius             [m]
alt_cible = 770e3;      % Target orbit altitude         [m]

a = r_earth + alt_cible;  % Reference orbit semi-major axis [m]
n = sqrt(mu / a^3);       % Mean motion                     [rad/s]

m_satellite = 1600;       % Initial chaser mass  [kg]
Isp         = 220;        % Specific impulse     [s]  (assumption)
g0          = 9.80665;    % Standard gravity     [m/s^2]

fprintf('=== Close-Range LQR Stop-and-Go Approach (100 m -> 2 m) ===\n');
fprintf('Mean motion (n)  : %.6e rad/s\n', n);
fprintf('Chaser mass      : %.1f kg\n',    m_satellite);
fprintf('Isp (assumed)    : %.1f s\n\n',   Isp);

%% ---- 2. SIMULATION PARAMETERS ------------------------------------------

% Initial state: [x(R-bar); y(V-bar); z(H-bar); vx; vy; vz; delta_V_cumul]
%   y = 100 m  ->  chaser is 100 m ahead along V-bar (target is at y = 0)
%   All velocities zero (starting from rest relative to target)
%   7th state accumulates delta-V for propellant budget
s0 = [0; 100; 0;   % Position  [m]
      0;   0; 0;   % Velocity  [m/s]
      0];          % Delta-V cumulator [m/s]

% Waypoints along V-bar (y-axis).  The chaser stops at each before continuing.
% Add or remove points to modify the approach corridor.
waypoints_y = [50, 10, 2];   % [m]  intermediate + final stop
waypoints_x = [ 0,  0, 0];   % [m]  keep on V-bar (R-bar offset = 0)

t_per_segment  = 3000;   % Max integration time per segment [s]
F_max_per_axis = 2.0;    % Thruster saturation per axis      [N]

%% ---- 3. HCW LINEAR STATE-SPACE MODEL -----------------------------------
%
%  Full 6-state HCW model:
%    d/dt [x; y; z; vx; vy; vz] = A * [x; y; z; vx; vy; vz] + B * u
%
%  where u = [ax; ay; az] is the applied acceleration (control input).
%
%  The 7th state (delta-V integrator) is handled inside hcw_dynamics_lqr
%  and is NOT part of A/B (those remain 6x6 and 6x3).

A = [ 0      0    0    1    0    0  ;
      0      0    0    0    1    0  ;
      0      0    0    0    0    1  ;
      3*n^2  0    0    0    2*n  0  ;
      0      0    0   -2*n  0    0  ;
      0      0   -n^2  0    0    0  ];

B = [ 0 0 0 ;
      0 0 0 ;
      0 0 0 ;
      1 0 0 ;
      0 1 0 ;
      0 0 1 ];

% Sanity check: HCW has two zero eigenvalues (drift modes) and two
% purely imaginary pairs. Print eigenvalues for verification.
fprintf('A matrix eigenvalues:\n');
disp(eig(A));

%% ---- 4. LQR CONTROLLER DESIGN ------------------------------------------
%
%  Cost function: J = integral( x'*Q*x + u'*R*u ) dt
%
%  Q penalises state error (position more than velocity here).
%  R penalises control effort (delta-V cost).
%
%  Tuning guidelines:
%    - Increase Q  -> faster convergence, more thrust usage.
%    - Increase R  -> more fuel-efficient, slower convergence.

max_pos_err = 100;    % Normalisation: max expected position error [m]
max_vel_err = 0.1;    % Normalisation: max expected velocity error [m/s]
max_accel   = F_max_per_axis / m_satellite;   % Max acceleration per axis [m/s^2]

Q = diag([ 100/max_pos_err^2,  ...    % x (R-bar)  -- weighted x10 vs y,z
             1/max_pos_err^2,  ...    % y (V-bar)
             1/max_pos_err^2,  ...    % z (H-bar)
             1/max_vel_err^2,  ...    % vx
             1/max_vel_err^2,  ...    % vy
             1/max_vel_err^2 ]);      % vz

R = diag([ 1/max_accel^2, ...        % ax
           1/max_accel^2, ...        % ay
           1/max_accel^2 ]);         % az

K = lqr(A, B, Q, R);
fprintf('LQR gain matrix K computed.\n\n');

%% ---- 5. SEQUENTIAL SEGMENT SIMULATION ----------------------------------
%
%  Each segment runs ode45 independently from the current state to the next
%  waypoint.  The final state of segment i becomes the initial state of i+1,
%  guaranteeing a continuous trajectory.
%
%  The feedforward term u_ff cancels the constant terms in A*x_ref so that
%  the equilibrium of the closed-loop system is exactly x_ref (not zero).
%  For waypoints on V-bar with zero velocity, u_ff is numerically ~0 but
%  is included for generality.

fprintf('Starting sequential segment simulation...\n');

n_seg  = length(waypoints_y);
T_list = cell(n_seg, 1);
S_list = cell(n_seg, 1);

s_current = s0;
options   = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);

for i = 1:n_seg

    % Target state for this segment (6-element mechanical state)
    x_ref = [waypoints_x(i); waypoints_y(i); 0;   % Position target [m]
             0;              0;              0];    % Velocity target [m/s]

    % Feedforward acceleration to hold the target point at equilibrium
    u_ff_accel = -pinv(B) * A * x_ref;   % (3x1) [m/s^2]

    % Build ODE function handle with captured parameters
    ode_func = @(t, s) hcw_dynamics_lqr(t, s, A, B, K, x_ref, ...
                                         u_ff_accel, m_satellite, F_max_per_axis);

    [T_seg, S_seg] = ode45(ode_func, [0, t_per_segment], s_current, options);

    T_list{i} = T_seg;
    S_list{i} = S_seg;

    % The end state (7x1) becomes the start of the next segment
    s_current = S_seg(end, :)';

    fprintf('  Segment %d/%d  (target y = %4.0f m)  done. '  , ...
            i, n_seg, waypoints_y(i));
    fprintf('Final y = %.3f m,  |dv_seg| = %.4f m/s\n', ...
            S_seg(end, 2), S_seg(end, 7) - S_seg(1, 7));
end

fprintf('All segments complete.\n\n');

%% ---- 6. POST-PROCESSING ------------------------------------------------

% --- 6a. Assemble full time history with absolute timestamps ---
T_full = [];  S_full = [];
Fx_hist = []; Fy_hist = []; Fz_hist = [];
t_offset = 0;
T_waypoints = [0];   % Timestamps of segment boundaries for plot markers

for i = 1:n_seg

    T_seg = T_list{i};
    S_seg = S_list{i};

    x_ref      = [waypoints_x(i); waypoints_y(i); 0; 0; 0; 0];
    u_ff_accel = -pinv(B) * A * x_ref;

    if i == 1
        T_full = T_seg;
        S_full = S_seg;
    else
        T_full = [T_full; T_seg(2:end) + t_offset];   %#ok<AGROW>
        S_full = [S_full; S_seg(2:end, :)];            %#ok<AGROW>
    end
    t_offset    = T_full(end);
    T_waypoints = [T_waypoints; t_offset];             %#ok<AGROW>

    % Re-compute applied force history for plotting
    Fx_seg = zeros(length(T_seg), 1);
    Fy_seg = zeros(length(T_seg), 1);
    Fz_seg = zeros(length(T_seg), 1);

    for j = 1:length(T_seg)
        s_mech   = S_seg(j, 1:6)';
        s_err    = s_mech - x_ref;
        u_cmd    = (-K * s_err + u_ff_accel) * m_satellite;   % [N]
        Fx_seg(j) = max(-F_max_per_axis, min(F_max_per_axis, u_cmd(1)));
        Fy_seg(j) = max(-F_max_per_axis, min(F_max_per_axis, u_cmd(2)));
        Fz_seg(j) = max(-F_max_per_axis, min(F_max_per_axis, u_cmd(3)));
    end

    if i == 1
        Fx_hist = Fx_seg;  Fy_hist = Fy_seg;  Fz_hist = Fz_seg;
    else
        Fx_hist = [Fx_hist; Fx_seg(2:end)];   %#ok<AGROW>
        Fy_hist = [Fy_hist; Fy_seg(2:end)];   %#ok<AGROW>
        Fz_hist = [Fz_hist; Fz_seg(2:end)];   %#ok<AGROW>
    end
end

% --- 6b. Unpack state columns ---
pos_x = S_full(:,1);  pos_y = S_full(:,2);  pos_z = S_full(:,3);
vel_x = S_full(:,4);  vel_y = S_full(:,5);  vel_z = S_full(:,6);
delta_v_cumul = S_full(:,7);

%% ---- 7. FIGURES --------------------------------------------------------

% Helper: draw vertical dashed lines at each waypoint transition
draw_segments = @(ax) arrayfun(@(t) xline(t, 'k--', 'LineWidth', 1, ...
                                'Parent', ax), T_waypoints);

% --- Figure 1: 3-D LVLH Trajectory ---
figure('Name', '3-D Trajectory (LVLH)', 'NumberTitle', 'off', 'Color', 'w');
plot3(pos_y, pos_x, pos_z, 'b-', 'LineWidth', 1.5); hold on;
plot3(s0(2), s0(1), s0(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
plot3(2,  0, 0, 'rx',  'LineWidth', 2.5,  'MarkerSize', 12);
plot3(50, 0, 0, 'ms',  'MarkerFaceColor', 'm', 'MarkerSize', 10);
plot3(10, 0, 0, 'cd',  'MarkerFaceColor', 'c', 'MarkerSize', 10);
xlabel('V-bar / Y (m)'); ylabel('R-bar / X (m)'); zlabel('H-bar / Z (m)');
title('Close-Range Approach — 3-D LVLH Trajectory');
legend('Trajectory', 'Start (100 m)', 'Final target (2 m)', ...
       'Waypoint (50 m)', 'Waypoint (10 m)', 'Location', 'best');
xlim([-10 105]); ylim([-10 10]); zlim([-1 1]);
grid on; view(45, 25); hold off;

% --- Figure 2: Position and Velocity vs. Time ---
figure('Name', 'States vs. Time', 'NumberTitle', 'off', 'Color', 'w');

ax1 = subplot(2, 1, 1);
plot(T_full, pos_x, 'r-', T_full, pos_y, 'g-', T_full, pos_z, 'b-', ...
     'LineWidth', 1.2);
draw_segments(ax1);
ylabel('Position [m]'); title('Position vs. Time');
legend('x (R-bar)', 'y (V-bar)', 'z (H-bar)', 'Segment boundary');
grid on;

ax2 = subplot(2, 1, 2);
plot(T_full, vel_x, 'r-', T_full, vel_y, 'g-', T_full, vel_z, 'b-', ...
     'LineWidth', 1.2);
draw_segments(ax2);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
legend('vx', 'vy', 'vz');  grid on;

% --- Figure 3: Thruster Force History ---
figure('Name', 'Control Effort', 'NumberTitle', 'off', 'Color', 'w');
ax3 = axes;
plot(T_full, Fx_hist, 'r-', T_full, Fy_hist, 'g-', T_full, Fz_hist, 'b-', ...
     'LineWidth', 1.2);
hold on;
yline( F_max_per_axis, 'k:', 'LineWidth', 1.5, 'Label', '+F_{max}');
yline(-F_max_per_axis, 'k:', 'LineWidth', 1.5, 'Label', '-F_{max}');
draw_segments(ax3);
xlabel('Time [s]'); ylabel('Force [N]');
title(sprintf('LQR Thrust Commands (saturated at \\pm%.1f N per axis)', F_max_per_axis));
legend('Fx (R-bar)', 'Fy (V-bar)', 'Fz (H-bar)', 'Saturation limit');
ylim([-F_max_per_axis*1.3, F_max_per_axis*1.3]);
grid on; hold off;

% --- Figure 4: Cumulative Delta-V ---
figure('Name', 'Cumulative Delta-V', 'NumberTitle', 'off', 'Color', 'w');
ax4 = axes;
plot(T_full, delta_v_cumul, 'k-', 'LineWidth', 1.5);
draw_segments(ax4);
xlabel('Time [s]'); ylabel('\DeltaV [m/s]');
title('Cumulative Delta-V vs. Time');
legend('\DeltaV total', 'Segment boundary', 'Interpreter', 'tex');
grid on;

% --- Figure 5: Cross-track error (deviation from ideal V-bar line) ---
%  Ideal V-bar approach stays at x = 0, z = 0.
%  Any excursion in x or z is an undesired cross-track deviation.
error_x          = pos_x;   % R-bar deviation [m]
error_z          = pos_z;   % H-bar deviation [m]
error_cross_track = sqrt(error_x.^2 + error_z.^2);  % Radial error norm [m]

figure('Name', 'Cross-Track Error', 'NumberTitle', 'off', 'Color', 'w');
ax5 = axes;
plot(T_full, error_x,          'r--', 'LineWidth', 1.2); hold on;
plot(T_full, error_z,          'b--', 'LineWidth', 1.2);
plot(T_full, error_cross_track,'k-',  'LineWidth', 2.0);
draw_segments(ax5);
xlabel('Time [s]'); ylabel('Position error [m]');
title('Cross-Track Error (deviation from ideal V-bar line)');
legend('x error (R-bar)', 'z error (H-bar)', ...
       '\surd(x^2 + z^2)  radial error', 'Segment boundary', ...
       'Interpreter', 'tex');
grid on; hold off;

%% ---- 8. MISSION BUDGET -------------------------------------------------

total_delta_v  = delta_v_cumul(end);
ve             = Isp * g0;                              % Effective exhaust velocity [m/s]
m_final        = m_satellite / exp(total_delta_v / ve); % Tsiolkovsky
m_propellant   = m_satellite - m_final;

fprintf('=== MISSION BUDGET ===\n');
fprintf('  Total Delta-V            : %.4f m/s\n', total_delta_v);
fprintf('  Propellant mass (Isp=%g s): %.3f kg\n', Isp, m_propellant);
fprintf('======================\n');

%% =========================================================================
%%  LOCAL FUNCTION: hcw_dynamics_lqr
%%
%%  Called at each ODE step. Computes the time derivative of the 7-state
%%  vector [x; y; z; vx; vy; vz; cumulative_delta_v].
%%
%%  Inputs:
%%    t           - current time (unused, autonomous system)
%%    s           - 7x1 state vector
%%    A, B        - HCW state-space matrices (6x6, 6x3)
%%    K           - LQR gain matrix (3x6)
%%    x_ref       - 6x1 target state
%%    u_ff_accel  - 3x1 feedforward acceleration [m/s^2]
%%    m           - chaser mass [kg]
%%    F_max       - saturation threshold per axis [N]
%%
%%  Outputs:
%%    dsdt        - 7x1 state derivative
%% =========================================================================

function dsdt = hcw_dynamics_lqr(~, s, A, B, K, x_ref, u_ff_accel, m, F_max)

    % Unpack mechanical states (indices 1-6)
    s_mech = s(1:6);

    % --- LQR control law ---
    s_err      = s_mech - x_ref;                    % State error
    u_cmd_accel = -K * s_err + u_ff_accel;           % Commanded acceleration [m/s^2]

    % --- Thruster saturation (per axis) ---
    F_cmd     = u_cmd_accel * m;                     % Convert to force [N]
    F_applied = max(-F_max, min(F_max, F_cmd));      % Clamp each element

    % --- Applied acceleration (after saturation) ---
    u_applied = F_applied / m;                       % [m/s^2]

    % --- State derivatives ---
    dsdt        = zeros(7, 1);
    dsdt(1:6)   = A * s_mech + B * u_applied;        % HCW dynamics
    dsdt(7)     = norm(u_applied);                   % Delta-V rate [m/s^2]  (integrates to m/s)

end
