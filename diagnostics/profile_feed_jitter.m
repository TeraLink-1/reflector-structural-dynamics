clc; clear; close all;

% initial parameters, all in SI
m = 2;
g = 9.81;
R = 1.0;
f = 0.91;
n_truss = 3;
c = 0.02; % damping constant per truss, initial estimate

% excitation force
F_mag = 1; % amplitude
F_hz = 1; % frequency
k_vals = linspace(1, 500, 800);
E_al = 200e9;

di = 0.05; % 0.05 m diameter rods
I_rod = pi * (di)^2/ 64;
A_rod = pi * (di)^2/ 2;

% truss length
L = sqrt(R^2 + f^2);  % assume roughly 45 degree angle trusses

omega = 2 * pi * F_hz;
period = 1/F_hz;

k_rod = E_al * A_rod/L; 

angles = 2 * pi * (0:n_truss-1)'/n_truss;
q0 = [0; 0; f]; % position @ equilibrium

% Angles are 0, 2*pi/n_truss, 4*pi/n_truss, ... around the z axis
angles = 2*pi*(0:n_truss-1)'/n_truss;

% One row per fixed point: [x_i, y_i, z_i]. z = 0 for every anchor
fixed_pts = [-R*cos(angles), R*sin(angles), zeros(n_truss,1)];


d_ref = q0 - fixed_pts.';               % matrix = n_trusses x 3.
u_ref = d_ref ./ vecnorm(d_ref,2,1);    % unit vectors
% small motion damping matrix:
C_eq = c*(u_ref * u_ref.');             % 3 x 3 matrix [N*s/m]

L0_from_k = @(kNow) L + m*g*L/(n_truss * kNow * f);
fprintf('Point-load CoG [m]: (0, 0, %.3f)\n',f);
fprintf('Current truss length L = %.6f m\n',L);
fprintf('Ideal rod EA/L = %.3e N/m (not the sweep k)\n',k_rod);
fprintf('Example free length L0 at k=100 N/m: %.6f m\n',L0_from_k(100));


relTol = 1e-8;
absTol = 1e-10;
options = odeset('RelTol',relTol,'AbsTol',absTol);
maxPeriods = 250;                  % Maximum cycles for initial settling
orbitTolerance = 1e-7;             % Scaled state repeat tolerance
nOutput = 401;                     % Points over one forcing cycle

sway_rms = nan(size(k_vals));
sway_max = nan(size(k_vals));
mean_x = nan(size(k_vals));
valid = false(size(k_vals));

k_vals;
size = size(k_vals);
sway_rms = zeros(size); sway_max = zeros(size); mean_x = zeros(size);

% Number of forcing periods to simulate before recording sway
nTransientCycles = 200; % we ignore the transient period for a bit

for j = [1 2 100 400 800]
    fprintf('START j=%d k=%.6g\n', j, k_vals(j));
    caseTimer = tic;

    % Current spring stiffness [N/m]
    kj = k_vals(j);

    % Calculate spring length at-rest.
    % This ensures that the initial position is in static equilibrium.
    % The springs must balance the gravitational force.
    L0 = L + m*g*L/(n_truss*kj*f);

    % Initial state vector:
    % z = [x; y; z; vx; vy; vz]
    % Mass starts at static equilibrium with zero velocity.
    z0 = [q0; 0; 0; 0];

    % Define the nonlinear equations of motion: z0 = rhs.
    rhs = @(t,z) trussODE(t, z, m, g, kj, c, ...
                         L0, fixed_pts, F_mag, omega);

    % Simulate initial motion to allow transients to decay.
    tTransient = nTransientCycles * period;
    [tTrans, Ztransient] = ode45(rhs, [0 tTransient], z0, options);

    fprintf('TRANSIENT seconds=%.3f rows=%d final_t=%.6g\n', toc(caseTimer), numel(tTrans), tTrans(end));
    % Use final state as initial condition for the measurement period.
    zSteady = Ztransient(end,:).';

    % Simulate one additional forcing period.
    tOne = linspace(tTransient, tTransient + period, nOutput);
    [tOne, Z] = ode45(rhs, tOne, zSteady, options);

    % Calculate horizontal sway using pythagorean
    sway = hypot(Z(:,1)-q0(1), Z(:,2)-q0(2));

    % Calculate RMS sway.
    sway_rms(j) = sqrt(trapz(tOne, sway.^2)/period);

    % Calculate maximum sway.
    sway_max(j) = max(sway);

    % Calculate mean horizontal x-position.
    mean_x(j) = trapz(tOne, Z(:,1))/period;

    % Mark simulation as completed.
    valid(j) = true;
    fprintf('DONE seconds=%.3f sway_rms=%.6g\n', toc(caseTimer), sway_rms(j));

end

function dz = trussODE(t, z, m, g, kj, c, L0, fixed_pts, F_mag, omega)

    % Extract position and velocity from state vector
    r = z(1:3);       % Position [x; y; z]
    v = z(4:6);       % Velocity [vx; vy; vz]

    % Initialize force vector with gravity and applied force
    F_vec = [F_mag*cos(omega*t); 0; -m*g];

    n_truss = size(fixed_pts,1);

    % Calculate spring and damping forces for each truss
    for i = 1:n_truss

        % Vector from fixed point to mass
        d = r - fixed_pts(i,:).';

        % Current truss length
        ell = norm(d);

        % Unit vector along truss
        u = d / ell;

        extension = ell - L0;

        v_axial = dot(v,u);

        F_vec = F_vec - (kj*extension + c*v_axial)*u;

    end

    a = F_vec / m;

    % Return state derivatives
    dz = [v; a];

end