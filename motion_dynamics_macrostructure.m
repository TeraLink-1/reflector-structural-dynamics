%[text] # Structural Dynamics of a Cantilever-Supported Reflector
%[text] This script builds a reduced-order model of a large reflector dish ground station to answer the question of how much wind and gravity deform the primary reflector surface. In other words, we want to determine and extract a preliminary *residual* surface RMS, a useful metric for tolerance stack ups. The point of creating a symbolic walk-through is to help us as a team determine at a quick glance the parametric relationships (i.e. what factors beyond the obvious may be driving the tower's resonance) from first principles.
%[text] 
%[text] ### Conventions
%[text] SI units throughout. The tower axis is $z$ (up from the clamped base), and lateral tower deflection is $y(z,t)$. The reflector axis (boresight) is tilted by the elevation angle $\\epsilon$ above the horizon, and $r$ is the radial coordinate measured from the reflector axis. A positive reflector displacement points along the boresight, toward the focus.
%[text] 
%[text] ### Prerequisites
%[text] Assuming a background in Statics and beam bending (Euler-Bernoulli), ordinary differential equations, and multivariable calculus. Plate theory is introduced as needed.
%[text] 
%[text] ### Assumptions / Constraints
%[text] The subreflector, feed, and trusses - that is, everything above the primary reflector, are structural excluded here. They are instead presented as a point-load mass. This is unlike in CAD where we can model a point-load mass's position and distribution. Thermal effects are also omitted from this study.
%%
%[text] ## Inputs & Symbolics
%[text] - Hollow tube: integrating over the annular cross-section gives $A=\\pi(d\_o^2-d\_i^2)/4$ and $I=\\int y^2\\,\\mathrm dA=\\pi(d\_o^4-d\_i^4)/64$. The mass per unit length is $\\mu=\\rho A$.
%[text] - Thin plate: bending rigidity $D=Eh^3/\[12(1-\\nu^2)\]$. Bending stiffness grows with $h^3$, so halving the thickness makes the plate eight times more flexible.
%[text] - Paraboloid $z\_0=r^2/(4f)$ only its shape parameter $\\beta=R^2/(4f^2)$ is used here, to weight the RMS surface error by true surface area.
%[text] - Wind: dynamic pressure $p\_\\infty=\\frac12\\rho\_{air}V^2$, modulated by a sinusoidal gust, $p(t)=p\_\\infty\[1+g\_f\\sin(\\Omega t)\]$. The wind is horizontal and blows in the elevation plane: `wind_direction` $s=+1$ means it strikes the concave face, $s=-1$ the back.
%[text] - Ruze's equation gives the aperture efficiency of a reflector with RMS surface error $\\sigma$: $\\eta=\\exp\\left\[-(4\\pi\\sigma/\\lambda)^2\\right\]$. Requiring $\\eta\\ge\\eta\_{min}$ gives $\\sigma\_{max}=\\frac{\\lambda}{4\\pi}\\sqrt{-\\ln\\eta\_{min}}$. \
%[text] **Modeling options.** Six option variables (set in the Simulation Model section) select between alternative modeling choices; each option is explained in the section that uses it. The defaults are the Ritz-consistent choices; the legacy choices of earlier versions of this script are kept for comparison.
clearvars; close all; clc;
assert(license('test','Symbolic_Toolbox') == 1, ...
    'This Live Script requires MATLAB Symbolic Math Toolbox.');

% =====================================================================
% RAW INPUTS. Every physical quantity is one named variable with its SI unit.
% =====================================================================
% ---- Tower: hollow aluminium tube, clamped at the ground ----
tower_length              = 2.5;       % L      [m]
tower_youngs_modulus      = 69e9;      % E_t    [Pa]
tower_density             = 2700;      % rho    [kg/m^3]
tower_outer_diameter      = 0.2;       % d_o    [m]
tower_wall_thickness      = 0.05;      % [m]
tower_tip_mass            = 45;        % M_tip  [kg], mount + reflector hub
tower_damping_ratio       = 0.02;      % zeta_t [-]
tower_drag_coefficient    = 1;         % C_d    [-]
tower_wind_width          = tower_outer_diameter;   % width the wind sees [m]

% ---- Reflector: CFRP, isotropic idealization ----
dish_radius               = 1.0;       % R      [m]
dish_focal_length         = 1.0;       % f      [m]
dish_thickness            = 0.0127;    % h      [m] (0.5 in)
dish_youngs_modulus       = 60e9;      % E_d    [Pa]
dish_poisson_ratio        = 0.28;      % nu     [-]
dish_structure_mass       = 45;        % M_dish [kg]
dish_damping_ratio        = 0.015;     % zeta_d [-] (assumed)
dish_hub_radius           = 0.10;      % a      [m], rigid central hub
dish_pressure_coefficient = 1.0;       % C_p    [-]
dish_wind_normal_factor   = 1.0;       % k_n0   [-], normal-pressure factor at face-on incidence
dish_elevation_deg        = 90;        % epsilon for single-angle runs [deg]

% ---- Environment ----
gravity_accel             = 9.81;      % g      [m/s^2]
air_density               = 1.225;     % rho_air [kg/m^3]
wind_speed                = 6;         % V      [m/s]
gust_fraction             = 0.30;      % g_f    [-]
gust_frequency            = 0.5;       % [Hz]
wind_direction            = +1;        % s: +1 onto the concave face, -1 from behind

% ---- RF requirement ---- (only for efficiency calculation)
rf_frequency_Hz           = 225e9;     % [Hz]
rf_min_efficiency         = 0.60;      % eta_min [-]
speed_of_light            = 299792458; % [m/s], exact

% ---- Simulation settings ----
sim_duration              = 5;         % [s]
sim_min_num_samples       = 1001;      % minimum number of output samples
sim_samples_per_period    = 20;        % output samples per period of the fastest mode
ode_rel_tol               = 1e-8;      % ode45 relative tolerance
ode_abs_tol               = 1e-10;     % ode45 absolute tolerance
plot_elevations_deg       = [0 45 90]; % elevations compared in the dish checks and the wind sweep [deg]
wind_speed_sweep          = 0:0.5:30;  % wind speeds for the RMS-versus-wind study [m/s]
rms_mesh_rings            = 4000;      % rings in the mesh-based RMS check [-]


% Derived variables ---------

tower_inner_diameter   = tower_outer_diameter - 2*tower_wall_thickness;
tower_area             = pi*(tower_outer_diameter^2 - tower_inner_diameter^2)/4;    % A   [m^2]
tower_second_moment    = pi*(tower_outer_diameter^4 - tower_inner_diameter^4)/64;   % I   [m^4]
tower_line_mass        = tower_density*tower_area;                                  % mu  [kg/m]
dish_flexural_rigidity = dish_youngs_modulus*dish_thickness^3/(12*(1 - dish_poisson_ratio^2)); % D [N m]
dish_depth_beta        = dish_radius^2/(4*dish_focal_length^2);                     % beta = (rim slope)^2 [-]
dish_flat_area         = pi*dish_radius^2;                                          % projected area [m^2]
dish_areal_density     = dish_structure_mass/dish_flat_area;                        % m_A [kg/m^2]
dynamic_pressure       = 0.5*air_density*wind_speed^2;                              % p_inf [Pa]
gust_omega             = 2*pi*gust_frequency;                                       % Omega [rad/s]
rf_wavelength          = speed_of_light/rf_frequency_Hz;                            % lambda [m]
max_rms_error          = rf_wavelength/(4*pi)*sqrt(-log(rf_min_efficiency));        % Ruze sigma_max [m]


% Gusting dynamic pressure p(t) = p_inf*(1 + g_f*sin(Omega*t)) [Pa].
% This is a function of time because ode45 needs the load at whatever t it chooses.
wind_pressure_fun = @(t) dynamic_pressure*(1 + gust_fraction*sin(gust_omega*t));


disp(table(tower_area,tower_second_moment,tower_line_mass,dish_flexural_rigidity, ... %[output:group:6b8c12a4] %[output:392983c0]
    dynamic_pressure,1e6*max_rms_error, ... %[output:392983c0]
    'VariableNames',["TubeArea_m2","TubeI_m4","LineMass_kg_m", ... %[output:392983c0]
    "PlateD_Nm","DynPressure_Pa","RuzeLimit_um"])); %[output:group:6b8c12a4] %[output:392983c0]
%[text] 
%[text] 
% Symbolic variables: the derivations below are done in letters, numbers come in only at the end.
syms z s r real                                % coordinates (tower height, dummy, radius)
syms L E_t I_t A_t F P_dist M_tip mu g real    % tower
syms nu R a p0 D M_dish beta real              % reflector
syms C1 C2 C3 C4 Ca real                       % integration constants
assumeAlso([L E_t I_t A_t R D a] > 0);


% Each symbol and the number that replaces it. Converting a symbolic result to a number is
% always the same pattern:  double(subs(expression, phys_symbols, phys_values)).
phys_symbols = [L E_t I_t A_t M_tip mu g nu R a D M_dish beta];
phys_values  = [tower_length tower_youngs_modulus tower_second_moment tower_area ...
                tower_tip_mass tower_line_mass gravity_accel ...
                dish_poisson_ratio dish_radius dish_hub_radius ...
                dish_flexural_rigidity dish_structure_mass dish_depth_beta];
%%
%[text] ## 
%[text] Euler–Bernoulli theory (Model 1) models the tower as a continuous beam whose exact vibration modes come from solving a partial differential equation, while the Rayleigh–Ritz method (Models 2 and 3) assumes a single deflection shape and uses energy to reduce that equation to one mass-spring-damper ODE. 
%[text] 
%[text] ## Tower Model \#1: cantilever with a lateral tip force *F*
%[text] **Goal.** Find the deflected shape of the tower when a single force acts at its tip, and use that shape to build a one-degree-of-freedom dynamic model.
%[text] **Theory.** For small slopes the beam curvature is $\\kappa\\simeq y''(z)$. Integrating Hooke's law $\\sigma=E\\epsilon$ over the cross-section gives the moment-curvature law $M\_b=EIy''$. Equilibrium of a beam element without distributed load gives $EIy^{(4)}=0$, so the deflection is a cubic with four unknown constants:
%[text] $y(z)=C\_1+C\_2z+C\_3z^2+C\_4z^3$.
%[text] Four boundary conditions fix the constants. At the clamped base there is no deflection and no slope, $y(0)=y'(0)=0$. At the free tip the moment vanishes, $EIy''(L)=0$, and the internal shear balances the applied force, $EIy'''(L)=-F$.
%[text] **Implementation.** Each boundary condition is a symbolic equation built with `subs` and `diff`; `solve`, called with one output per unknown, returns the constants, which are collected into a row vector. The normalized shape $\\phi\_t=y/y(L)$ satisfies $\\phi\_t(L)=1$, so the tip deflection itself becomes the generalized coordinate later.
% General solution of EI*y'''' = 0: a cubic with four unknown constants
deflection_general = C1 + C2*z + C3*z^2 + C4*z^3;

% Boundary conditions: clamped base (y = y' = 0), free tip with force F (M = 0, shear = -F)
bc_tip_load = [subs(deflection_general,z,0) == 0; ...
               subs(diff(deflection_general,z),z,0) == 0; ...
               subs(E_t*I_t*diff(deflection_general,z,2),z,L) == 0; ...
               subs(E_t*I_t*diff(deflection_general,z,3),z,L) == -F];

% solve returns one output per unknown, in the order the unknowns are listed
[sol_C1,sol_C2,sol_C3,sol_C4] = solve(bc_tip_load,[C1 C2 C3 C4]);
tipload_constants      = simplify([sol_C1 sol_C2 sol_C3 sol_C4]);
tipload_deflection     = simplify(subs(deflection_general,[C1 C2 C3 C4],tipload_constants));  % y(z)
tipload_tip_deflection = simplify(subs(tipload_deflection,z,L));                              % delta = y(L)
tipload_shape          = simplify(tipload_deflection/tipload_tip_deflection);                 % phi_t, phi_t(L) = 1
tipload_stiffness      = simplify(F/tipload_tip_deflection);                                  % static spring F/delta
disp('Integration constants [C1 C2 C3 C4]:'); disp(tipload_constants); %[output:39226161] %[output:8337e5bb]
disp('Deflection y(z):');                    disp(tipload_deflection); %[output:42c5d1f8] %[output:534651bc]
disp('Tip deflection delta = y(L):');        disp(tipload_tip_deflection); %[output:01d3e9f7] %[output:401ac056]
disp('Normalized shape phi_t(z) = y/delta:'); disp(tipload_shape); %[output:119e2848] %[output:7a20901a]
%[text] ### 
%[text] ### Parametric Model
%[text] **Goal.** Simulate the tower as a dynamic model when the wind switches on.
%[text] **Theory (Rayleigh-Ritz).** Assume the moving tower keeps the static shape, $y(z,t)=q(t)\\phi\_t(z)$, where $q$ is the tip deflection. Writing kinetic energy $T$ and strain energy $V$ in terms of $q$ gives
%[text] $T=\\frac12m\_1\\dot q^2,\\quad m\_1=M\_{tip}+\\mu\\int\_0^L\\phi\_t^2\\,\\mathrm dz;\\qquad V=\\frac12k\_1q^2,\\quad k\_1=EI\\int\_0^L(\\phi\_t'')^2\\,\\mathrm dz$,
%[text] and the dynamical equation gives $m\_1\\ddot q+c\_1\\dot q+k\_1q=F(t)$. The damping coefficient follows from the assumed damping ratio, $c=2\\zeta m\\omega\_n$ with $\\omega\_n=\\sqrt{k/m}$. For the tip-shape the stiffness equals the static spring constant, $k\_1=F/\\delta=3EI/L^3$.
%[text] This first model is a lumped model. All forces are exerted on that point, i.e. $F(t)=P\_\\ell(t)L$, and drag per unit height $P\_\\ell=C\_d d\_o\\,p(t)$.
%[text] **Implementation.**
%[text] - The Ritz integrals are written out symbolically for the trial shape, numbers are substituted with `subs` and `double`, and `make_sdof` turns mass, stiffness and damping ratio into the natural frequency and the damping coefficient.
%[text] - `simulate_sdof` integrates the ODE with `ode45`. The peak sway includes the overshoot caused by switching the wind on at $t=0$; the steady-state peak is taken after four decay time constants $4/(\\zeta\\omega\_n)$ for comparison. \
% ---- Rayleigh-Ritz integrals for the tip-load shape phi_t (symbolic) ----
tip_shape_mass_sym      = M_tip + mu*int(tipload_shape^2,z,0,L);          % m_1 = M_tip + mu*int(phi^2) dz
tip_shape_stiffness_sym = E_t*I_t*int(diff(tipload_shape,z,2)^2,z,0,L);   % k_1 = EI*int(phi''^2) dz
tip_shape_force_sym     = int(tipload_shape,z,0,L);                       % Q per unit line load = int(phi) dz
% For the tip-load shape, Ritz must reproduce the static spring constant 3EI/L^3
assert(isAlways(simplify(tip_shape_stiffness_sym - tipload_stiffness) == 0));

% ---- Substitute numbers and build the one-node model ----
tip_shape_mass      = double(subs(tip_shape_mass_sym,phys_symbols,phys_values));        % [kg]
tip_shape_stiffness = double(subs(tip_shape_stiffness_sym,phys_symbols,phys_values));   % [N/m]
[tip_shape_omega_n,tip_shape_freq_Hz,tip_shape_damping] = ...
    make_sdof(tip_shape_mass,tip_shape_stiffness,tower_damping_ratio);

% Wind drag per unit tower height, P_l(t) = C_d * d_o * p(t) [N/m]
tower_drag_per_length = @(t) tower_drag_coefficient*tower_wind_width*wind_pressure_fun(t);

% ---- Tower model #1, lumped-model: all drag lumped at the tip, so Q = P_l * L ----
tipeq_force_per_line_load = tower_length;                % [m]

% Output time grid shared by the phi_t tower models: >= sim_samples_per_period points per period
tower_num_samples = max(sim_min_num_samples,ceil(sim_samples_per_period*tip_shape_freq_Hz*sim_duration) + 1);
tower_time_grid   = linspace(0,sim_duration,tower_num_samples).';

% Integrate m*q'' + c*q' + k*q = Q*P_l(t), starting at rest (the wind switches on at t = 0)
[tower1_time,tower1_sway] = simulate_sdof(tip_shape_mass,tip_shape_damping,tip_shape_stiffness, ...
    tower_time_grid,@(t) tipeq_force_per_line_load*tower_drag_per_length(t),[0;0],ode_rel_tol,ode_abs_tol);



% Steady-state peak: ignore the start-up wobble.
% When the wind switches on at t = 0, the tower rings for a while before it settles.
% The ringing shrinks by a factor e (about 2.7) every 1/(zeta*omega_n) seconds,
% so after 4/(zeta*omega_n) seconds less than 2% of it is left.
% Only the motion after that time is used for the steady-state peak.

tower1_settle_time = 4/(tower_damping_ratio*tip_shape_omega_n);       % [s]
is_steady = tower1_time > tower1_settle_time;
if ~any(is_steady), is_steady(:) = true; end                         % run too short: use it all
tower1_steady_peak = max(abs(tower1_sway(is_steady)));               % [m]
fprintf('Tower #1: k1 = %.4g N/m, m1 = %.4g kg, f1 = %.3f Hz, peak tip sway = %.3f um (steady state %.3f um)\n', ... %[output:group:4e18bd20] %[output:22e6a343]
    tip_shape_stiffness,tip_shape_mass,tip_shape_freq_Hz,1e6*max(abs(tower1_sway)), ... %[output:22e6a343]
    1e6*tower1_steady_peak); %[output:group:4e18bd20] %[output:22e6a343]
%%
%[text] ## 
%[text] ## Tower Model \#2: cantilever under a distributed lateral load
%[text] **Goal.** Wind actually acts along the whole height. Solve that case exactly and find out how much error the tip-force idealization of Model \#1 introduces.
%[text] **Theory.** With a constant load per unit height $P\_\\ell$ (N/m) equilibrium becomes $EIy^{(4)}=P\_\\ell$. A particular solution is $P\_\\ell z^4/(24EI)$, added to the same cubic as before. The base is still clamped, but the tip now carries **no** applied force or moment, so $y''(L)=y'''(L)=0$.
%[text] **Implementation.** Identical to Model \#1 except for the particular solution and the tip conditions. That is the benefit of the symbolic approach: changing the physics means changing two lines.
% Same procedure as Model #1, with the particular solution P*z^4/(24EI) added
deflection_general = P_dist*z^4/(24*E_t*I_t) + C1 + C2*z + C3*z^2 + C4*z^3;
% Clamped base; free tip with no applied force or moment (y'' = y''' = 0)
bc_line_load = [subs(deflection_general,z,0) == 0; ...
                subs(diff(deflection_general,z),z,0) == 0; ...
                subs(diff(deflection_general,z,2),z,L) == 0; ...
                subs(diff(deflection_general,z,3),z,L) == 0];
[sol_C1,sol_C2,sol_C3,sol_C4] = solve(bc_line_load,[C1 C2 C3 C4]);
lineload_constants          = simplify([sol_C1 sol_C2 sol_C3 sol_C4]);
lineload_deflection         = simplify(subs(deflection_general,[C1 C2 C3 C4],lineload_constants));
lineload_tip_deflection     = simplify(subs(lineload_deflection,z,L));
lineload_shape              = simplify(lineload_deflection/lineload_tip_deflection);        % phi_line
lineload_base_moment        = simplify(E_t*I_t*subs(diff(lineload_deflection,z,2),z,0));   % exact base moment
lineload_force_on_tip_shape = simplify(int(P_dist*tipload_shape,z,0,L));                    % virtual work on phi_t
disp('Uniform-load constants [C1 C2 C3 C4]:'); disp(lineload_constants); %[output:084e8dce] %[output:9fd24dea]
disp('Uniform-load deflection y(z):');        disp(lineload_deflection); %[output:47fbc80d] %[output:4791f0f2]
disp('Exact base moment M_b(0):');            disp(lineload_base_moment); %[output:7c04463d] %[output:539c5648]
disp('Generalized force of the line load on the tip shape:'); disp(lineload_force_on_tip_shape); %[output:3b4ef50b] %[output:82196bdf]
%[text] ### 
%[text] **Theory.** The generalized force in a Ritz model comes from **virtual work**: if the tower moves by $\\delta q$, a distributed load does work $\\int\_0^LP\_\\ell\\,\\phi(z)\\,\\delta q\\,\\mathrm dz$, so $Q=P\_\\ell\\int\_0^L\\phi\\,\\mathrm dz$. Three models are compared:
%[text] - lumped-model (Model \#1): all drag at the tip, $Q=P\_\\ell L$.
%[text] - **Ritz with the tip-load shape:** $\\phi\_t$ with $Q=\\frac38P\_\\ell L$.
%[text]     - We use this for vibration and dynamics. (i.e. periodic, gusty wind)
%[text] - **Ritz with the uniform-load shape:** $\\phi\_\\ell=y\_\\ell/\\delta\_\\ell$ with its own $m$, $k$ and $Q$.
%[text]     - We use this for static analysis (i.e. quasi-static wind) \
% ---- Rayleigh-Ritz integrals for the uniform-load shape phi_line (symbolic) ----
line_shape_mass_sym      = M_tip + mu*int(lineload_shape^2,z,0,L);
line_shape_stiffness_sym = E_t*I_t*int(diff(lineload_shape,z,2)^2,z,0,L);
line_shape_force_sym     = int(lineload_shape,z,0,L);

% Exact static results per unit line load, for reference
exact_tip_deflection_per_load = double(subs(lineload_tip_deflection/P_dist,phys_symbols,phys_values)); % [m per N/m]
exact_base_moment_per_load    = double(subs(lineload_base_moment/P_dist,phys_symbols,phys_values));    % [N m per N/m]

% ---- Tower model #2a, Ritz phi_t: same mass and stiffness as Model #1, virtual-work force ----
ritz_tip_force_per_line_load = double(subs(tip_shape_force_sym,phys_symbols,phys_values));   % = 3L/8 [m]

% ---- Tower model #2b, Ritz phi_line: its own mass, stiffness and force ----
ritz_line_mass      = double(subs(line_shape_mass_sym,phys_symbols,phys_values));        % [kg]
ritz_line_stiffness = double(subs(line_shape_stiffness_sym,phys_symbols,phys_values));   % [N/m]
[ritz_line_omega_n,ritz_line_freq_Hz,ritz_line_damping] = ...
    make_sdof(ritz_line_mass,ritz_line_stiffness,tower_damping_ratio);
ritz_line_force_per_line_load = double(subs(line_shape_force_sym,phys_symbols,phys_values));  % [m]

fprintf('%-16s f = %8.4f Hz\n','lumped-model',tip_shape_freq_Hz); %[output:26fee239]
fprintf('%-16s f = %8.4f Hz\n','Ritz phi_t',tip_shape_freq_Hz); %[output:3f79f281]
fprintf('%-16s f = %8.4f Hz\n','Ritz phi_{line}',ritz_line_freq_Hz); %[output:977f159f]

% Ritz phi_t time history: needed for bending shortening (Tower Model #3) and the tower plot
[~,tower2_sway] = simulate_sdof(tip_shape_mass,tip_shape_damping,tip_shape_stiffness, ...
    tower_time_grid,@(t) ritz_tip_force_per_line_load*tower_drag_per_length(t),[0;0],ode_rel_tol,ode_abs_tol);

% Ritz phi_line time history, on its own grid because its frequency differs slightly
line_num_samples = max(sim_min_num_samples,ceil(sim_samples_per_period*ritz_line_freq_Hz*sim_duration) + 1);
[ritz_line_time,ritz_line_sway] = simulate_sdof(ritz_line_mass,ritz_line_damping,ritz_line_stiffness, ...
    linspace(0,sim_duration,line_num_samples).', ...
    @(t) ritz_line_force_per_line_load*tower_drag_per_length(t),[0;0],ode_rel_tol,ode_abs_tol);
%%
%[text] ## 
%[text] ## Tower Model \#3: gravitational compression and bending shortening
%[text] **Goal.** Find how much the tower shortens under its own weight and the tip mass, and how much additional shortening the wind sway produces.
%[text] **Theory.** Cut the tower at height $z$: the compressive force is the weight above the cut, $N(z)=M\_{tip}g+\\mu g(L-z)$. Hooke's law for axial deformation is $EA\\,u\_a'(z)=-N(z)$, a first-order equation with one constant $C\_a$, fixed by the clamped base $u\_a(0)=0$.
%[text] A bent beam also shortens geometrically. Integrating gives the second-order shortening $u\_b(z,t)=-B(z)q(t)^2$ with
%[text] $B(z)=\\frac12\\int\_0^z\\left(\\phi\_t'(s)\\right)^2\\mathrm ds$.
% Axial force at height z: weight of the tip mass plus the tube above z
axial_force_sym        = M_tip*g + mu*g*(L-z);

% Hooke's law EA*u' = -N(z), integrated once; the constant Ca comes from u(0) = 0
axial_disp_general     = simplify(int(-axial_force_sym/(E_t*A_t),z) + Ca);
axial_constant         = simplify(solve(subs(axial_disp_general,z,0) == 0,Ca));
axial_displacement_sym = simplify(subs(axial_disp_general,Ca,axial_constant));

% Bending shortening coefficient B(z) = (1/2) int_0^z phi_t'(s)^2 ds
shortening_coeff_sym   = simplify(int(subs(diff(tipload_shape,z)^2,z,s),s,0,z)/2);
disp('Axial constant Ca:'); disp(axial_constant); %[output:060b5a6f] %[output:63616f02]
disp('Gravity-induced axial displacement u_a(z):'); disp(axial_displacement_sym); %[output:14c49967] %[output:2b48a988]
disp('Bending-shortening coefficient B(z):');     disp(shortening_coeff_sym); %[output:08550b27] %[output:10d40b67]
%[text] ### 
%[text] ### Parametric Model
%[text] **Implementation**
% Evaluate the symbolic fields at the tip (z = L) and the base (z = 0)
gravity_tip_settlement = double(subs(subs(axial_displacement_sym,z,L),phys_symbols,phys_values));  % u_a(L) [m]
shortening_coeff_tip   = double(subs(subs(shortening_coeff_sym,z,L),phys_symbols,phys_values));    % B(L) [1/m]
bending_shortening_tip = -shortening_coeff_tip*tower2_sway.^2;      % u_b(L,t) = -B(L) q(t)^2 [m]
base_axial_stress      = -double(subs(subs(axial_force_sym,z,0),phys_symbols,phys_values))/tower_area;  % [Pa]

% First axial frequency with a linear axial shape (one-term estimate)
axial_freq_Hz = sqrt((tower_youngs_modulus*tower_area/tower_length)/ ...
    (tower_tip_mass + tower_line_mass*tower_length/3))/(2*pi);

fprintf(['Tower #3: gravity tip settlement = %.4f um, base axial stress = %.4f MPa\n' ... %[output:group:33760346] %[output:5f8f2385]
    '          peak bending shortening = %.4g pm, axial f_a = %.1f Hz (%.0fx twice the lateral f1)\n'], ... %[output:5f8f2385]
    1e6*gravity_tip_settlement,1e-6*base_axial_stress,1e12*max(abs(bending_shortening_tip)), ... %[output:5f8f2385]
    axial_freq_Hz,axial_freq_Hz/(2*tip_shape_freq_Hz)); %[output:group:33760346] %[output:5f8f2385]
%[text] ## 
%[text] ### Gravity-adjusted (geometric) stiffness of the tower
%[text] **Goal.** Gravity compresses the tower, and a compressed column is easier to bend sideways (in the limit it buckles). Quantify this compressing.
%[text] **Theory.** With $y=q\\phi\_t$, the generalized mass is $M\_t=M\_{tip}+\\mu\\int\_0^L\\phi\_t^2\\,\\mathrm dz$ and the elastic bending stiffness is $K\_b=EI\\int\_0^L(\\phi\_t'')^2\\,\\mathrm dz$. When the column bends, the weight above each section moves down by the bending shortening of Tower Model \#3, so gravity does **positive** work and the total potential energy drops by $\\frac12K\_gq^2$, with
%[text] $K\_g=\\int\_0^LN(z)\\left(\\phi\_t'\\right)^2\\mathrm dz$.
%[text] The effective stiffness is $K\_t=K\_b-K\_g$. Multiplying gravity by a load factor $\\lambda$ gives $\\omega(\\lambda)=\\sqrt{(K\_b-\\lambda K\_g)/M\_t}$, which reaches zero at $\\lambda\_{cr}=K\_b/K\_g$: the one-term Ritz estimate of the self-weight buckling load. Rayleigh-Ritz always over-estimates stiffness, so $\\lambda\_{cr}$ is an **upper bound** on the true buckling multiplier.
% Generalized mass, elastic stiffness and geometric (gravity) stiffness for y = q*phi_t
tower_gen_mass_sym            = simplify(M_tip + mu*int(tipload_shape^2,z,0,L));
tower_bending_stiffness_sym   = simplify(E_t*I_t*int(diff(tipload_shape,z,2)^2,z,0,L));
tower_geometric_stiffness_sym = simplify(int(axial_force_sym*diff(tipload_shape,z)^2,z,0,L));
disp('Generalized mass Mt:');            disp(tower_gen_mass_sym); %[output:413b5a00] %[output:16d986d4]
disp('Elastic bending stiffness Kb:');   disp(tower_bending_stiffness_sym); %[output:2f6cc261] %[output:6527b51d]
disp('New, reduced stiffness Kg:'); disp(tower_geometric_stiffness_sym); %[output:5cca9466] %[output:7f70666a]
%[text] ### 
%[text] ## Simulation Model
% ---- Modelling options ----
tower_wind_model   = 'ritz_distributed';    % 'ritz_distributed' | 'tip_equivalent'
base_mass_model    = 'ritz_participation';  % 'ritz_participation' | 'modal_mass'
initial_condition  = 'gravity_sag';         % 'gravity_sag' | 'mean_wind'
coupled_plate_case = 'B';                   % 'B' clamped hub (real support) | 'A' clamped rim
dish_wind_model    = 'cos_elevation';       % 'cos_elevation' | 'constant' (legacy)
dish_drag_on_tower = true;                  % apply the dish's wind drag at the tower tip
assert(any(strcmp(coupled_plate_case,["A","B"])),'coupled_plate_case must be ''A'' or ''B''.');

% ---- Normal-pressure factor k_n of a horizontal wind on the tilted dish ----
% Evaluated once: at the design elevation, and at each comparison elevation.
switch dish_wind_model
    case 'cos_elevation'    % only the wind component along the boresight pushes on the surface
        design_normal_factor  = dish_wind_normal_factor*abs(cosd(dish_elevation_deg));
        compare_normal_factor = dish_wind_normal_factor*abs(cosd(plot_elevations_deg));
    case 'constant'         % legacy: the same push at every elevation
        design_normal_factor  = dish_wind_normal_factor;
        compare_normal_factor = dish_wind_normal_factor*ones(size(plot_elevations_deg));
    otherwise
        error('Unknown dish_wind_model "%s".',dish_wind_model);
end

%[text] **Implementation.** The gravity-compressed tower ($K\_t$) is simulated with the same Ritz $\\phi\_t$ wind force as Tower Model \#2, and a figure at the end of this section compares the tip sway of every tower model. This section also assembles the `coupled_tower_*` variables, the tower model used in the coupled tower-dish simulation. Its wind participation follows `tower_wind_model`:
%[text] - `'ritz_distributed'` (default): $Q=\\frac38P\_\\ell L$, the virtual-work force of the real distributed wind.
%[text] - `'lumped-model'`: $Q=P\_\\ell L$, the legacy choice, which over-predicts the tower sway. \
%[text] The remaining options (`coupled_plate_case`, `dish_wind_model`, `dish_drag_on_tower`) are explained in the coupled tower-dish section.
%[text] 
% ---- Numbers for M_t, K_b, K_g ----
tower_gen_mass            = double(subs(tower_gen_mass_sym,phys_symbols,phys_values));             % M_t [kg]
tower_bending_stiffness   = double(subs(tower_bending_stiffness_sym,phys_symbols,phys_values));    % K_b [N/m]
tower_geometric_stiffness = double(subs(tower_geometric_stiffness_sym,phys_symbols,phys_values));  % K_g [N/m]
tower_effective_stiffness = tower_bending_stiffness - tower_geometric_stiffness;                    % K_t [N/m]
buckling_load_factor      = tower_bending_stiffness/tower_geometric_stiffness;                      % lambda_cr [-]
assert(tower_effective_stiffness > 0,'Tower effective stiffness is nonpositive (the tower buckles).');
% M_t and K_b are the same integrals as the Ritz phi_t model, so that model is the "K_b only" case
assert(abs(tower_gen_mass/tip_shape_mass - 1) < 1e-12 && abs(tower_bending_stiffness/tip_shape_stiffness - 1) < 1e-12);
[kt_omega_n,kt_freq_Hz,kt_damping] = make_sdof(tower_gen_mass,tower_effective_stiffness,tower_damping_ratio);

% gravity-compressed tower under the same Ritz phi_t wind force
[~,kt_sway] = simulate_sdof(tower_gen_mass,kt_damping,tower_effective_stiffness, ...
    tower_time_grid,@(t) ritz_tip_force_per_line_load*tower_drag_per_length(t),[0;0],ode_rel_tol,ode_abs_tol);
fprintf(['Without gravity the stiffness is %.4f%% too high (Kb/Kt - 1).\n' ... %[output:group:4a2c197a] %[output:918f967a]
    'f_n: %.4f Hz (Kb) -> %.4f Hz (Kt); Ritz self-weight buckling at %.0f g\n'], ... %[output:918f967a]
    100*(tower_bending_stiffness/tower_effective_stiffness-1),tip_shape_freq_Hz,kt_freq_Hz,buckling_load_factor); %[output:group:4a2c197a] %[output:918f967a]
%[text] 
% ---- Tower model used in the coupled tower-dish simulation: Ritz phi_t with K_t ----
switch tower_wind_model
    case 'ritz_distributed', coupled_tower_force_per_line_load = ritz_tip_force_per_line_load;
    case 'tip_equivalent',   coupled_tower_force_per_line_load = tipeq_force_per_line_load;
    otherwise, error('Unknown tower_wind_model "%s".',tower_wind_model);
end
coupled_tower_mass      = tower_gen_mass;
coupled_tower_stiffness = tower_effective_stiffness;
coupled_tower_damping   = kt_damping;
coupled_tower_omega_n   = kt_omega_n;
coupled_tower_freq_Hz   = kt_freq_Hz;
coupled_tower_tip_slope_per_q = double(subs(subs(diff(tipload_shape,z),z,L),phys_symbols,phys_values));  % = 3/(2L) [1/m]
%[text] ### 
%[text] ## Tip sway of every tower model
%[text] **Goal.** Put every tower model on one plot, all driven by the same gusting wind and all started at rest, to see how much each modelling choice matters.
%[text] - **\#1 lumped-model:** all drag at the tip, $Q=P\_\\ell L$. It sways $8/3$ times more than the Ritz $\\phi\_t$ model, because it gives the whole wind load the leverage of the tip.
%[text] - **\#2 Ritz** $\\phi\_t$ **and Ritz** $\\phi\_{line}$**:** both use the virtual-work force of the distributed wind. The two trial shapes give almost the same curve: both reproduce the exact static tip deflection $P\_\\ell L^4/(8EI)$, and their frequencies differ by about 1%.
%[text] - **Ritz** $\\phi\_t$ **with** $K\_t=K\_b-K\_g$**:** the gravity-compressed tower used in the coupled simulation. It lies on top of the $\\phi\_t$ curve; the printout above gives the difference in percent. \
%[text] The first fraction of a second shows the start-up overshoot (the wind switches on at $t=0$); after about $4/(\\zeta\\omega\_n)$ every curve settles into steady forced motion at the 0.5 Hz gust frequency.
figure('Color','w','Name','Tower models: tip sway q(t)'); %[output:66480527]
plot(tower1_time,1e6*tower1_sway,'LineWidth',1.2); hold on;               %[output:66480527]
plot(tower_time_grid,1e6*tower2_sway,'--','LineWidth',1.2); %[output:66480527]
plot(ritz_line_time,1e6*ritz_line_sway,'-.','LineWidth',1.2); %[output:66480527]
plot(tower_time_grid,1e6*kt_sway,':','LineWidth',1.8); %[output:66480527]
grid on; xlabel('Time [s]'); ylabel('Tip sway q [μm]'); %[output:66480527]
title('Tower-tip sway of each tower model under the same gusting wind'); %[output:66480527]
legend('#1 lumped-model (Q = P_\ell L)','#2 Ritz \phi_t (Q = 3P_\ell L/8)','#2 Ritz \phi_{line}', ... %[output:group:234a25cb] %[output:66480527]
       'Ritz \phi_t with K_b - K_g (coupled-model tower)','Location','best'); %[output:group:234a25cb] %[output:66480527]
%[text] ### Peak tip sway of every tower model
%[text] **Peak** is the largest sway over the whole run, so it includes the start-up overshoot. **Steady peak** is the largest sway after four decay time constants of that model. The last column compares the steady peak with the exact static tip deflection of the real distributed wind at the gust peak: a value near 1 means the model gets the physics right, and the lumped-model model shows its 8/3 over-prediction.
% ---- Table: peak tip sway of each tower model ----
tower_model_names  = ["#1 lumped-model";"#2 Ritz phi_t";"#2 Ritz phi_line";"Ritz phi_t, K_b - K_g"];
tower_freqs_Hz     = [tip_shape_freq_Hz; tip_shape_freq_Hz; ritz_line_freq_Hz; kt_freq_Hz];          % f_n [Hz]
tower_force_factor = [tipeq_force_per_line_load; ritz_tip_force_per_line_load; ...
                      ritz_line_force_per_line_load; ritz_tip_force_per_line_load];                 % Q per unit line load [m]
% Settling time of each model: four decay time constants, 4/(zeta*omega_n) [s]
tower_settle_times = 4./(tower_damping_ratio*[tip_shape_omega_n; tip_shape_omega_n; ritz_line_omega_n; kt_omega_n]);
tower_peak_um   = zeros(4,1);       % largest |q| over the whole run [um]
tower_steady_um = zeros(4,1);       % largest |q| after settling [um]
for model_idx = 1:4
    switch model_idx                % pick this model's time history
        case 1, sway_time = tower1_time;     sway = tower1_sway;
        case 2, sway_time = tower_time_grid; sway = tower2_sway;
        case 3, sway_time = ritz_line_time;  sway = ritz_line_sway;
        case 4, sway_time = tower_time_grid; sway = kt_sway;
    end
    is_steady = sway_time > tower_settle_times(model_idx);
    if ~any(is_steady), is_steady(:) = true; end                 % run too short: use it all
    tower_peak_um(model_idx)   = 1e6*max(abs(sway));
    tower_steady_um(model_idx) = 1e6*max(abs(sway(is_steady)));
end
% Reference: exact static tip deflection of the distributed wind at the gust peak,
% (tip deflection per unit line load) x (drag per unit height at the gust peak)
exact_static_tip_um = 1e6*exact_tip_deflection_per_load* ...
    tower_drag_coefficient*tower_wind_width*dynamic_pressure*(1 + gust_fraction);
fprintf('Exact static tip deflection of the distributed wind at the gust peak: %.3f um\n',exact_static_tip_um); %[output:54806142]
disp(table(tower_freqs_Hz,tower_force_factor,tower_peak_um,tower_steady_um,tower_steady_um/exact_static_tip_um, ... %[output:group:1399f49a] %[output:91478ab6]
    'RowNames',tower_model_names, ... %[output:91478ab6]
    'VariableNames',["f_n_Hz","Q_per_line_load_m","PeakSway_um","SteadyPeak_um","SteadyPeak_over_exact_static"])); %[output:group:1399f49a] %[output:91478ab6]
%%
%[text] # Surface RMS Calculations
%[text] ## Dish Model \#1: flat circular plate under uniform pressure
%[text] **Goal.** Before tackling the curved reflector, solve the simplest surface model exactly: a flat disk of the same material and thickness.
%[text] **Theory (Kirchhoff-Love plate).** A thin plate obeys $D\\nabla^4w=p\_0$, the two-dimensional analogue of $EIy^{(4)}=P$. For a field that depends only on $r$, the Laplacian is $\\nabla^2(\\cdot)=\\frac1r\\frac{\\mathrm d}{\\mathrm dr}\\left(r\\frac{\\mathrm d(\\cdot)}{\\mathrm dr}\\right)$. Integrating the resulting ordinary differential equation four times gives
%[text] $w(r)=\\frac{p\_0r^4}{64D}+C\_1+C\_2r^2+C\_3\\ln r+C\_4r^2\\ln r$.
%[text] The plate's internal forces follow from $w$: radial moment $M\_r=-D\\left(w''+\\nu w'/r\\right)$ and transverse shear $Q\_r=-D\\,\\frac{\\mathrm d}{\\mathrm dr}\\nabla^2w$.
%[text] #### 
%[text] #### Case A: solid disk, clamped rim
%[text] At $r=0$ the field must be finite and smooth. The $\\ln r$ term makes $w(0)$ infinite, so $C\_3=0$. The $r^2\\ln r$ term is finite, but its second derivative $2\\ln r+3$ makes $M\_r(0)$ infinite (it represents a point force at the centre), so $C\_4=0$. The clamped rim gives $w(R)=w'(R)=0$, and
%[text] $w\_A(r)=\\frac{p\_0(R^2-r^2)^2}{64D}$.
%[text] 
%[text] #### Case B: annulus, clamped at center with a free rim
%[text] This is how the real reflector is supported. Now $r=0$ is not part of the plate, so all four basis functions are allowed and all four constants come from the edges:
%[text] - **Clamped hub:** no deflection and no slope, $w(a)=0$ and $w'(a)=0$.
%[text] - **Free rim:** no moment and no shear, $M\_r(R)=0$ and $Q\_r(R)=0$. (In general the free-edge condition uses Kirchhoff's effective shear $V\_r=Q\_r+\\frac1r\\partial M\_{r\\theta}/\\partial\\theta$, which reduces to $Q\_r$ for axisymmetric deformation.) \
assumeAlso(r > 0); assumeAlso(R > a);
% General axisymmetric solution of D*lap(lap(w)) = p0
plate_deflection_general = p0*r^4/(64*D) + C1 + C2*r^2 + C3*log(r) + C4*r^2*log(r);
% Check it: apply the axisymmetric Laplacian lap(w) = w'' + w'/r twice
plate_laplacian   = diff(plate_deflection_general,r,2) + diff(plate_deflection_general,r)/r;
plate_bilaplacian = diff(plate_laplacian,r,2) + diff(plate_laplacian,r)/r;
assert(isAlways(simplify(D*plate_bilaplacian - p0) == 0));

% ---- Case A: solid disk, clamped rim (regularity removes C3 and C4) ----
caseA_deflection_general = subs(plate_deflection_general,[C3 C4],[0 0]);
bc_caseA = [subs(caseA_deflection_general,r,R) == 0; ...            % w(R) = 0
            subs(diff(caseA_deflection_general,r),r,R) == 0];        % w'(R) = 0
[sol_C1,sol_C2] = solve(bc_caseA,[C1 C2]);
plateA_constants  = simplify([sol_C1 sol_C2]);
plateA_deflection = simplify(subs(caseA_deflection_general,[C1 C2],plateA_constants));
disp('Case A (solid disk, clamped rim) deflection w_A(r):'); disp(plateA_deflection); %[output:3dae280a] %[output:1af315bb]

% ---- Case B: annulus, clamped hub at r = a, free rim at r = R ----
% Dimensionless form: x = r/R, rho = a/R, wbar = w*D/(p0*R^4). Fewer symbols make solve faster.
syms x rho positive
caseB_dimless_general = x^4/64 + C1 + C2*x^2 + C3*log(x) + C4*x^2*log(x);
bc_caseB = [subs(caseB_dimless_general,x,rho) == 0; ...                                  % w(a) = 0
    subs(diff(caseB_dimless_general,x),x,rho) == 0; ...                                  % w'(a) = 0
    subs(diff(caseB_dimless_general,x,2) + nu*diff(caseB_dimless_general,x)/x,x,1) == 0; ... % M_r(R) = 0
    subs(diff(diff(caseB_dimless_general,x,2) + diff(caseB_dimless_general,x)/x,x),x,1) == 0]; % Q_r(R) = 0
[sol_C1,sol_C2,sol_C3,sol_C4] = solve(bc_caseB,[C1 C2 C3 C4]);
plateB_dimless_constants = simplify([sol_C1 sol_C2 sol_C3 sol_C4]);
assert(isAlways(plateB_dimless_constants(4) == -sym(1)/8));

% Dimensional field for the rest of the script: w = (p0 R^4 / D) * wbar(r/R, a/R)
plateB_deflection = simplify(p0*R^4/D*subs(subs(caseB_dimless_general,[C1 C2 C3 C4], ...
    plateB_dimless_constants),[x rho],[r/R a/R]));
disp('Case B (annulus, clamped hub, free rim) deflection w_B(r):'); disp(plateB_deflection); %[output:7c3b4b46] %[output:10153e16]
% ---- Plot w(r) of both support cases under the reflector's self-weight ----
plot_elevation_deg = 0;   % measured from the horizon: 0 = dish faces the horizon, 90 = points straight up

% Self-weight of the reflector, spread over its projected area
self_weight_N        = dish_structure_mass*gravity_accel;           % total weight [N]
self_weight_pressure = self_weight_N/dish_flat_area;               % m_A*g [Pa]
% Split the weight: the normal part bends the plate, the in-plane part does not (in this model).
% Sign convention of this script: positive pressure and deflection point toward the focus.
normal_pressure  = -self_weight_pressure*sind(plot_elevation_deg);   % bending load [Pa]
inplane_pressure =  self_weight_pressure*cosd(plot_elevation_deg);   % in-plane load, not modelled [Pa]
fprintf(['Self-weight %.1f N = %.1f Pa over the dish; at %g deg elevation: ' ... %[output:group:524a4d9b] %[output:5c4fcaf6]
         'normal part %.1f Pa (bends the plate), in-plane part %.1f Pa (not modelled)\n'], ... %[output:5c4fcaf6]
        self_weight_N,self_weight_pressure,plot_elevation_deg,normal_pressure,inplane_pressure); %[output:group:524a4d9b] %[output:5c4fcaf6]

% Deflection per unit pressure (p0 = 1 Pa) for each case; the plate is linear, so w = p * w_unit(r)
caseA_w_unit = matlabFunction(subs(plateA_deflection,[phys_symbols p0],[phys_values 1]),'Vars',r);   % [m/Pa]
caseB_w_unit = matlabFunction(subs(plateB_deflection,[phys_symbols p0],[phys_values 1]),'Vars',r);   % [m/Pa]
% Deflection per pascal at each plate's reference point [m/Pa] (used by the Ritz checks in Dish Model #2)
caseA_compliance = caseA_w_unit(0);              % Case A: centre of the disk
caseB_compliance = caseB_w_unit(dish_radius);    % Case B: free rim

% Radius grids: Case A covers the whole disk; Case B only exists outside the hub (r >= a)
r_caseA = linspace(0,dish_radius,400);
r_caseB = linspace(dish_hub_radius,dish_radius,400);
% Deflection at the chosen elevation, and the full self-weight sag (dish pointing straight up) [um]
w_caseA_um      = 1e6*normal_pressure*caseA_w_unit(r_caseA);
w_caseB_um      = 1e6*normal_pressure*caseB_w_unit(r_caseB);
w_caseA_full_um = -1e6*self_weight_pressure*caseA_w_unit(r_caseA);
w_caseB_full_um = -1e6*self_weight_pressure*caseB_w_unit(r_caseB);

figure('Color','w','Name','Plate deflection w(r), self-weight');
plot(r_caseA,w_caseA_full_um,'LineWidth',1.2); hold on
plot(r_caseB,w_caseB_full_um,'LineWidth',1.2);
grid on; xlim([0 dish_radius]);
xlabel('Radius r [m]'); ylabel('Deflection w [μm] (positive toward the focus)');
title(sprintf('Self-weight deflection, elevation %g deg (normal pressure %.1f Pa)', ...
      plot_elevation_deg,normal_pressure));
legend(sprintf('Case A: clamped rim, %g deg',plot_elevation_deg), ... %[output:group:6e5ec5bc] %[output:28f3b431] %[output:8170fec3] %[output:62200cae] %[output:964446b7] %[output:8ad4db5f] %[output:5628c2f1] %[output:83eb0070]
       sprintf('Case B: clamped hub, %g deg',plot_elevation_deg), ... %[output:28f3b431] %[output:8170fec3] %[output:62200cae] %[output:964446b7] %[output:8ad4db5f] %[output:5628c2f1] %[output:83eb0070]
       'Case A: full self-weight (90 deg)','Case B: full self-weight (90 deg)','Location','southwest'); %[output:group:6e5ec5bc] %[output:28f3b431] %[output:8170fec3] %[output:62200cae] %[output:964446b7] %[output:8ad4db5f] %[output:5628c2f1] %[output:83eb0070]
%%
%[text] ## 
%[text] ## Dish Model \#2: Rayleigh-Ritz plate dynamics
%[text] **Goal.** Give each plate inertia, so it has a natural frequency. Case B, the real hub-mounted support, is the model used in the coupled simulation; Case A is kept for comparison.
%[text] **Theory.** Exactly as for the tower, assume the plate keeps its static shape, $w(r,t)=\\eta(t)\\psi(r)$, with $\\psi$ normalized to 1 at a reference point. Energy integrals over the plate area $\\mathrm dA=2\\pi r\\,\\mathrm dr$ give
%[text] - modal mass $M\_d=m\_A\\int\_A\\psi^2\\,\\mathrm dA$ (from kinetic energy),
%[text] - modal area $A\_m=\\int\_A\\psi\\,\\mathrm dA$, so a uniform pressure produces the generalized force $pA\_m$ (from virtual work; note $\\psi$, not $\\psi^2$),
%[text] - stiffness $K\_d=D\\int\_A\\left\[(\\psi'')^2+(\\psi'/r)^2+2\\nu\\psi''\\psi'/r\\right\]\\mathrm dA$ (from plate bending energy). \
%[text] For Case A the mixed $\\nu$ term integrates to $\\left\[(\\psi')^2/2\\right\]\_0^R=0$ because $\\psi'$ vanishes at both ends, so $\\nu$ enters only through $D$. For Case B it does not vanish, because $\\psi'(R)\\neq0$ at the free rim.
%[text] **Implementation.** Case A uses the closed-form integrals. Case B evaluates the same three integrals by adaptive quadrature (`integral`, inline) with the exact annulus field normalized at the rim, $\\psi\_B=w\_B/w\_B(R)$. Because each trial function is the exact static field of its case, Ritz must reproduce the exact static deflection, $A\_m/K\_d=w(r\_{ref})|\_{p\_0=1}$; the script asserts this for both cases.
% ---- Case A: Ritz integrals in closed form, with psi_A normalized to 1 at the centre ----
caseA_shape = simplify(plateA_deflection/subs(plateA_deflection,r,0)); %[output:3ae2b2de]
plateA_modal_mass_sym = simplify((M_dish/(pi*R^2))*int(caseA_shape^2*2*pi*r,r,0,R));   % m_A * int psi^2 dA %[output:96f52379] %[output:997e78f1]
plateA_modal_area_sym = simplify(int(caseA_shape*2*pi*r,r,0,R));                       % int psi dA %[output:87ec659a] %[output:59c68b18]
plateA_stiffness_sym  = simplify(D*int((diff(caseA_shape,r,2)^2 + (diff(caseA_shape,r)/r)^2 + ... %[output:group:1929771a] %[output:0caae702] %[output:9b5502ad] %[output:63bee956]
    2*nu*diff(caseA_shape,r,2)*diff(caseA_shape,r)/r)*2*pi*r,r,0,R));                   % bending energy %[output:group:1929771a] %[output:0caae702] %[output:9b5502ad] %[output:63bee956]
disp('Case A modal mass:'); disp(plateA_modal_mass_sym); %[output:4f8846f7] %[output:93140652] %[output:21ae7169] %[output:17bcaa4a] %[output:60c92f0a] %[output:29584962]
disp('Case A modal area:'); disp(plateA_modal_area_sym); %[output:69412218] %[output:42a17e7f] %[output:2581fbca] %[output:243f62ab] %[output:71f6998f]
disp('Case A stiffness:');  disp(plateA_stiffness_sym); %[output:3417f6f3] %[output:4581e01c] %[output:280d53ee]

% Case A numbers
dishA_mass       = double(subs(plateA_modal_mass_sym,phys_symbols,phys_values));   % M_d [kg] %[output:6e0a9860]
dishA_stiffness  = double(subs(plateA_stiffness_sym,phys_symbols,phys_values));    % K_d [N/m] %[output:222a7955] %[output:53128d67]
dishA_modal_area = double(subs(plateA_modal_area_sym,phys_symbols,phys_values));   % A_m [m^2] %[output:4b9019a7] %[output:1b76cd87]
[dishA_omega_n,dishA_freq_Hz,dishA_damping] = make_sdof(dishA_mass,dishA_stiffness,dish_damping_ratio); %[output:11014414] %[output:0b5d8d26]
% Numeric shape psi_A(r); matlabFunction turns the symbolic expression into a function of r
dishA_shape_fun  = matlabFunction(subs(caseA_shape,phys_symbols,phys_values),'Vars',r); %[output:94e1bb69] %[output:2dae094e]
dishA_label      = 'Disk A: clamped rim, one mode'; %[output:1983974e] %[output:42ebfabb]

% ---- Case B: same three integrals, by adaptive quadrature on the exact annulus shape ----
caseB_unit_deflection_sym = subs(subs(plateB_deflection,p0,1),phys_symbols,phys_values); %[output:13ead4b6] %[output:1967a5a5]
caseB_shape_sym     = caseB_unit_deflection_sym/subs(caseB_unit_deflection_sym,r,dish_radius);   % psi_B(R) = 1 %[output:72d5c4cb] %[output:6c3c8a5c]
dishB_shape_fun     = matlabFunction(caseB_shape_sym,'Vars',r);              % psi_B(r) %[output:8040fa78] %[output:8f34bbc8]
caseB_slope_fun     = matlabFunction(diff(caseB_shape_sym,r),'Vars',r);      % psi_B'(r) %[output:6e146d4e] %[output:206000bc]
caseB_curvature_fun = matlabFunction(diff(caseB_shape_sym,r,2),'Vars',r);    % psi_B''(r) %[output:71311b91]

% Each integrand carries the ring area element dA = 2*pi*r dr
dishB_mass = dish_areal_density*integral(@(rad) 2*pi*rad.*dishB_shape_fun(rad).^2, ... %[output:group:055f24ef] %[output:547e7ab3]
    dish_hub_radius,dish_radius,'RelTol',1e-12,'AbsTol',1e-16);                              % M_d [kg] %[output:group:055f24ef] %[output:547e7ab3]
dishB_stiffness = dish_flexural_rigidity*integral(@(rad) 2*pi*rad.*( ... %[output:group:152aefb3] %[output:29c1aaa4]
    caseB_curvature_fun(rad).^2 + (caseB_slope_fun(rad)./rad).^2 + ... %[output:29c1aaa4]
    2*dish_poisson_ratio*caseB_curvature_fun(rad).*caseB_slope_fun(rad)./rad), ... %[output:29c1aaa4]
    dish_hub_radius,dish_radius,'RelTol',1e-12,'AbsTol',1e-16);                              % K_d [N/m] %[output:group:152aefb3] %[output:29c1aaa4]
dishB_modal_area = integral(@(rad) 2*pi*rad.*dishB_shape_fun(rad), ... %[output:group:712781d5] %[output:4e06b12e]
    dish_hub_radius,dish_radius,'RelTol',1e-12,'AbsTol',1e-16);                              % A_m [m^2] %[output:group:712781d5] %[output:4e06b12e]
[dishB_omega_n,dishB_freq_Hz,dishB_damping] = make_sdof(dishB_mass,dishB_stiffness,dish_damping_ratio); %[output:6992b9ea]
dishB_label = 'Disk B: clamped hub, one mode'; %[output:5e0f0fed]

% Ritz with the exact static field must reproduce the exact static deflection: A_m/K_d = w(r_ref)|p0=1
assert(abs(dishA_modal_area/dishA_stiffness/caseA_compliance - 1) < 1e-8, ... %[output:group:2691baa7] %[output:8b325637]
    'Ritz static deflection for Case A does not match the exact field.'); %[output:group:2691baa7] %[output:8b325637]
assert(abs(dishB_modal_area/dishB_stiffness/caseB_compliance - 1) < 1e-8, ... %[output:group:38f74f0b] %[output:78d435e8]
    'Ritz static deflection for Case B does not match the exact field.'); %[output:group:38f74f0b] %[output:78d435e8]
fprintf('Dish #2 (Case %s): M_d = %.4g kg, K_d = %.4g N/m, A_m = %.4f m^2, f_n = %.3f Hz\n', ... %[output:group:0790a1e4] %[output:5162f37e] %[output:5871ab63]
    'A',dishA_mass,dishA_stiffness,dishA_modal_area,dishA_freq_Hz); %[output:group:0790a1e4] %[output:5162f37e] %[output:5871ab63]
fprintf('Dish #2 (Case %s): M_d = %.4g kg, K_d = %.4g N/m, A_m = %.4f m^2, f_n = %.3f Hz\n', ... %[output:group:988e55dc] %[output:6a4faa08] %[output:96c63b70]
    'B',dishB_mass,dishB_stiffness,dishB_modal_area,dishB_freq_Hz); %[output:group:988e55dc] %[output:6a4faa08] %[output:96c63b70]
%[text] ### 
%[text] ### Dynamic plate deflection versus time
%[text] **Goal.** Repeat the analysis of Dish Model \#1, but now let the plate have mass, so it can vibrate, and check whether the quasi-static answer is good enough.
%[text] **Idea.** The one-mode model behaves like a mass on a spring, $M\_d\\ddot\\eta+C\_d\\dot\\eta+K\_d\\eta=A\_m\\,p\_n(t)$, where $\\eta$ is the deflection at the reference point. The force is the same normal pressure as before, multiplied by the modal area $A\_m$.
%[text] **What to look for in the table.**
%[text] - **Start-up wobble.** The wind is switched on at $t=0$ (with `initial_condition = 'gravity_sag'`), so the plate rings at its natural frequency for a fraction of a second before it settles; the start-up peak column includes this ringing.
%[text] - **After settling, the dynamic peak matches the quasi-static peak** (ratio column close to 1). The gust frequency (0.5 Hz) is far below the plate's natural frequency (tens of Hz), so the plate has plenty of time to follow the load. That is why the quasi-static answer of Dish Model \#1 is trustworthy here. \
disp('Dish #2 check: peak plate deflection, dynamic vs quasi-static'); %[output:4aff5ece] %[output:045d0012]
% One row per (plate case, elevation); filled in by the loop, displayed as a table after it
num_check_rows      = 2*numel(plot_elevations_deg); %[output:61339828]
check_plate         = strings(num_check_rows,1); %[output:9b1f4b9c]
check_elevation_deg = zeros(num_check_rows,1); %[output:549dc5cd]
check_startup_um    = zeros(num_check_rows,1);    % largest |eta| over the whole run (includes ringing) %[output:7697b005]
check_steady_um     = zeros(num_check_rows,1);    % largest |eta| after settling, dynamic model %[output:08c1fc54]
check_quasi_um      = zeros(num_check_rows,1);    % largest |eta| after settling, quasi-static model %[output:214808a1]
check_row = 0; %[output:6fa828b1]
for case_idx = 1:2 %[output:group:309b5cea] %[output:3843f62b]
    % Pick the plate model for this pass of the loop
    if case_idx == 1
        plate_label = dishA_label;          plate_mass       = dishA_mass;
        plate_damping = dishA_damping;      plate_stiffness  = dishA_stiffness;
        plate_omega_n = dishA_omega_n;      plate_freq_Hz    = dishA_freq_Hz;
        plate_modal_area = dishA_modal_area;
        plate_short_name = "Case A (clamped rim)";
    else
        plate_label = dishB_label;          plate_mass       = dishB_mass;
        plate_damping = dishB_damping;      plate_stiffness  = dishB_stiffness;
        plate_omega_n = dishB_omega_n;      plate_freq_Hz    = dishB_freq_Hz;
        plate_modal_area = dishB_modal_area;
        plate_short_name = "Case B (clamped hub)";
    end
    % Time grid fine enough to resolve the plate vibration; settling time 4/(zeta*omega_n)
    plate_num_samples = max(sim_min_num_samples,ceil(sim_samples_per_period*plate_freq_Hz*sim_duration) + 1);
    plate_time_grid   = linspace(0,sim_duration,plate_num_samples).';
    plate_settle_time = 4/(dish_damping_ratio*plate_omega_n);
    for elev_idx = 1:numel(plot_elevations_deg)
        elev_deg      = plot_elevations_deg(elev_idx);
        normal_factor = compare_normal_factor(elev_idx);
        % Generalized force on the mode = A_m * p_n(t), where the normal pressure p_n has a
        % wind part -s*C_p*k_n*p(t) and a weight part -m_A*g*sin(epsilon)  [N]
        modal_force_fun = @(t) plate_modal_area*( ...
            -wind_direction*dish_pressure_coefficient*normal_factor*wind_pressure_fun(t) ...
            - dish_areal_density*gravity_accel*sind(elev_deg));
        % Starting state: gravity sag only, or already balanced under the mean wind
        switch initial_condition
            case 'gravity_sag'
                initial_state = [-plate_modal_area*dish_areal_density*gravity_accel*sind(elev_deg)/plate_stiffness; 0];
            case 'mean_wind'
                initial_state = [modal_force_fun(0)/plate_stiffness; 0];
            otherwise
                error('Unknown initial condition "%s".',initial_condition);
        end
        % Dynamic answer: integrate M_d*eta'' + C_d*eta' + K_d*eta = A_m*p_n(t)
        [response_time,response_deflection] = simulate_sdof(plate_mass,plate_damping,plate_stiffness, ...
            plate_time_grid,modal_force_fun,initial_state,ode_rel_tol,ode_abs_tol);
        % Quasi-static answer: the plate follows the load instantly, eta = A_m*p_n/K_d
        quasi_static_deflection = modal_force_fun(response_time)/plate_stiffness;
        % Compare peaks after the start-up transient
        is_steady = response_time > plate_settle_time;
        if ~any(is_steady), is_steady(:) = true; end
        dynamic_peak      = max(abs(response_deflection(is_steady)));
        quasi_static_peak = max(abs(quasi_static_deflection(is_steady)));
        % Store this row of the table
        check_row = check_row + 1;
        check_plate(check_row)         = plate_short_name;
        check_elevation_deg(check_row) = elev_deg;
        check_startup_um(check_row)    = 1e6*max(abs(response_deflection));
        check_steady_um(check_row)     = 1e6*dynamic_peak;
        check_quasi_um(check_row)      = 1e6*quasi_static_peak;
    end
end %[output:group:309b5cea]
disp(table(check_plate,check_elevation_deg,check_startup_um,check_steady_um,check_quasi_um, ... %[output:group:843ffefd] %[output:8d6cae1f] %[output:11d3d358] %[output:6f22ee26] %[output:3b0eddb6] %[output:71287fb0] %[output:48fa193f] %[output:146670de] %[output:6ba23b93] %[output:6cef442b] %[output:994ad0f4] %[output:6374ad02] %[output:7e9c66a3] %[output:43cac026] %[output:2d52a723] %[output:29ff9b26] %[output:3073f8ae] %[output:89017ec5] %[output:3b22a42b]
    check_steady_um./check_quasi_um, ... %[output:8d6cae1f] %[output:11d3d358] %[output:6f22ee26] %[output:3b0eddb6] %[output:71287fb0] %[output:48fa193f] %[output:146670de] %[output:6ba23b93] %[output:6cef442b] %[output:994ad0f4] %[output:6374ad02] %[output:7e9c66a3] %[output:43cac026] %[output:2d52a723] %[output:29ff9b26] %[output:3073f8ae] %[output:89017ec5] %[output:3b22a42b]
    'VariableNames',["Plate","Elevation_deg","StartupPeak_um","SteadyPeak_um","QuasiStaticPeak_um","Dynamic_over_QuasiStatic"])); %[output:group:843ffefd] %[output:8d6cae1f] %[output:11d3d358] %[output:6f22ee26] %[output:3b0eddb6] %[output:71287fb0] %[output:48fa193f] %[output:146670de] %[output:6ba23b93] %[output:6cef442b] %[output:994ad0f4] %[output:6374ad02] %[output:7e9c66a3] %[output:43cac026] %[output:2d52a723] %[output:29ff9b26] %[output:3073f8ae] %[output:89017ec5] %[output:3b22a42b]
%%
%[text] ## 
%[text] ## RMS over a true paraboloid surface
%[text] 
%[text] **Theory.** For $z\_0=r^2/(4f)$ the surface element is $\\mathrm dA\_s=2\\pi r\\sqrt{1+(r/2f)^2}\\,\\mathrm dr$. In the normalized area variable $u=(r/R)^2$ (equal steps in $u$ are equal projected areas) this is the projected element multiplied by the weight $\\gamma(u)=\\sqrt{1+\\beta u}$ with $\\beta=R^2/(4f^2)$. Every average in the fit is built from five moments $J\_k=\\int\_0^1u^k\\gamma(u)\\,\\mathrm du$, $k=0,\\ldots,4$, because the Case A shape $\\psi=(1-u)^2$ and $\\psi^2$ are polynomials of degree 2 and 4. The substitution $t=1+\\beta u$ turns each moment into a finite sum:
%[text] $J\_k=\\beta^{-(k+1)}\\sum\_{j=0}^{k}\\binom{k}{j}(-1)^{k-j}\\frac{(1+\\beta)^{j+3/2}-1}{j+3/2}$.
%[text] 
%[text] **Implementation.** The moments are summed symbolically and evaluated with 32-digit arithmetic (`vpa`), because the formula subtracts nearly equal numbers when $\\beta$ is small. The resulting residual factor is stored in `dishA_rms_factor`: the fitted RMS of Case A is $\\sigma=|\\eta|\\times$ this factor. The closed form exists only because $\\psi\_A$ is a polynomial in $u$. For Case B the same weighted fit, over the annulus $u\\in\[(a/R)^2,1\]$, is evaluated by adaptive quadrature (`integral`) in the loop below; running the same quadrature on Case A must reproduce the closed form, which checks it. The loop also gives the raw and piston-removed RMS and the flat-weighted values, which feed the surface-model comparison table later.
% ---- Case A closed form: the five area moments J_k, k = 0..4 ----
area_moments = sym(zeros(1,5)); %[output:8bbc42f2]
for k = 0:4 %[output:group:69e518c2] %[output:282d09a9]
    moment_sum = sym(0);
    for j = 0:k
        moment_sum = moment_sum + sym(nchoosek(k,j))*(-1)^(k-j)* ...
            ((1+beta)^(sym(j)+sym(3)/2) - 1)/(sym(j)+sym(3)/2);
    end
    area_moments(k+1) = moment_sum/beta^(k+1);
end %[output:group:69e518c2]
% Weighted averages for psi_A = (1-u)^2:  <f> = int f*gamma du / int gamma du
weighted_mean_u        = area_moments(2)/area_moments(1);                                   % <u> %[output:544dbbc2]
weighted_mean_u2       = area_moments(3)/area_moments(1);                                   % <u^2> %[output:4a9bae7b]
weighted_mean_psi      = (area_moments(1) - 2*area_moments(2) + area_moments(3))/area_moments(1);   % <psi> %[output:46ac682e]
weighted_mean_u_psi    = (area_moments(2) - 2*area_moments(3) + area_moments(4))/area_moments(1);   % <u psi> %[output:5f74f4ef]
weighted_mean_psi2     = (area_moments(1) - 4*area_moments(2) + 6*area_moments(3) ... %[output:group:2796b74d] %[output:05202769]
                          - 4*area_moments(4) + area_moments(5))/area_moments(1);            % <psi^2> %[output:group:2796b74d] %[output:05202769]
weighted_var_u         = weighted_mean_u2 - weighted_mean_u^2; %[output:56ab5d16]
weighted_defocus_coeff = (weighted_mean_u_psi - weighted_mean_u*weighted_mean_psi)/weighted_var_u;  % least-squares slope %[output:1acf5954] %[output:56d0c7b9]
weighted_residual_var  = weighted_mean_psi2 - weighted_mean_psi^2 - weighted_defocus_coeff^2*weighted_var_u; %[output:02d3cfb5] %[output:2d9b656a]
% Evaluate at the dish's beta with 32 significant digits (the formula cancels badly for small beta)
beta_vpa = vpa(sym(dish_depth_beta),32); %[output:5ed0aba6] %[output:579ceedb]
dishA_rms_raw     = double(vpa(subs(sqrt(weighted_mean_psi2),beta,beta_vpa),25)); %[output:63fcb2cd] %[output:3dc6cb3d]
dishA_rms_piston  = double(vpa(subs(abs(weighted_mean_psi),beta,beta_vpa),25)); %[output:3a91cc4e] %[output:5bb9428a]
dishA_rms_defocus = double(vpa(subs(sqrt(weighted_defocus_coeff^2*weighted_var_u),beta,beta_vpa),25)); %[output:965ca6fb] %[output:09815d6f]
dishA_rms_factor  = double(vpa(subs(sqrt(weighted_residual_var),beta,beta_vpa),25)); %[output:90180336] %[output:8e9ab67b]
disp('Physical-area RMS factors [raw, piston, defocus, residual]:'); %[output:85f9a4e2] %[output:64b5005d] %[output:2833c2a6] %[output:1920200f]
disp([dishA_rms_raw, dishA_rms_piston, dishA_rms_defocus, dishA_rms_factor]); %[output:09bb7a2a] %[output:35e4d1e8]

% ---- Same weighted fit by adaptive quadrature, for every reflector surface model ----
% A surface model = a support case (A or B) plus an area weighting
% (flat: beta = 0, or true paraboloid: beta = R^2/4f^2).
% Rows: A flat, A paraboloid, B flat, B paraboloid.
% Columns: raw RMS (nothing removed), RMS after piston removal, residual after piston + defocus.
surface_rms_factors = zeros(4,3); %[output:612eba1f]
for case_idx = 1:2 %[output:group:2e5ba625] %[output:009a3847]
    if case_idx == 1
        shape_fun = dishA_shape_fun;  u_inner = 0;                                 % solid disk
    else
        shape_fun = dishB_shape_fun;  u_inner = (dish_hub_radius/dish_radius)^2;   % annulus starts at the hub
    end
    shape_of_u = @(u) shape_fun(dish_radius*sqrt(u));        % psi written in terms of u = (r/R)^2
    for weight_idx = 1:2
        if weight_idx == 1, weight_beta = 0; else, weight_beta = dish_depth_beta; end
        area_weight = @(u) sqrt(1 + weight_beta*u);           % gamma(u): true area per projected area
        % Weighted averages <f> = int f*gamma du / int gamma du over the plate, u in [u_inner, 1]
        total_weight = integral(area_weight,u_inner,1,'RelTol',1e-12,'AbsTol',1e-15);
        mean_u     = integral(@(u) u.*area_weight(u),u_inner,1,'RelTol',1e-12,'AbsTol',1e-15)/total_weight;
        mean_u2    = integral(@(u) u.^2.*area_weight(u),u_inner,1,'RelTol',1e-12,'AbsTol',1e-15)/total_weight;
        mean_psi   = integral(@(u) shape_of_u(u).*area_weight(u),u_inner,1,'RelTol',1e-12,'AbsTol',1e-15)/total_weight;
        mean_u_psi = integral(@(u) u.*shape_of_u(u).*area_weight(u),u_inner,1,'RelTol',1e-12,'AbsTol',1e-15)/total_weight;
        mean_psi2  = integral(@(u) shape_of_u(u).^2.*area_weight(u),u_inner,1,'RelTol',1e-12,'AbsTol',1e-15)/total_weight;
        var_u         = mean_u2 - mean_u^2;
        defocus_coeff = (mean_u_psi - mean_u*mean_psi)/var_u;       % best-fit defocus coefficient
        row = 2*(case_idx - 1) + weight_idx;
        surface_rms_factors(row,:) = [sqrt(mean_psi2), ...                                    % raw
                                      sqrt(mean_psi2 - mean_psi^2), ...                       % piston removed
                                      sqrt(mean_psi2 - mean_psi^2 - defocus_coeff^2*var_u)];  % residual
    end
end %[output:group:2e5ba625]
% Residual factors under each assumption, by name
rmsA_flat_quad  = surface_rms_factors(1,3); %[output:2c6f6342]
rmsA_parab_quad = surface_rms_factors(2,3); %[output:7e616dc5]
rmsB_flat_quad  = surface_rms_factors(3,3); %[output:45cdd578]
rmsB_parab_quad = surface_rms_factors(4,3); %[output:00e7bf71]
dishB_rms_factor = rmsB_parab_quad;           % true-area residual factor of Case B %[output:9fa129d9]
% Check: quadrature on Case A must reproduce the closed form (raw RMS and residual)
assert(abs(surface_rms_factors(2,1)/dishA_rms_raw - 1) < 1e-8,'Raw RMS quadrature does not match the closed form.'); %[output:724c27a5]
assert(abs(rmsA_parab_quad/dishA_rms_factor - 1) < 1e-8,'Residual RMS quadrature does not match the closed form.'); %[output:7282728d]
fprintf('Residual RMS factor per unit amplitude: Case A %.5f (centre), Case B %.5f (rim)\n', ... %[output:group:10128288] %[output:40be6060] %[output:27015193]
    dishA_rms_factor,dishB_rms_factor); %[output:group:10128288] %[output:40be6060] %[output:27015193]
%%
%[text] ## Summary of the reduced-order models
%[text] **Goal.** See every model side by side before using them. This is the "derive once, evaluate many times" habit: each row was obtained by substituting numbers into a symbolic result, never by re-deriving.
%[text] **Check yourself.** Every forcing frequency in this problem is 0.5 Hz. Which rows tell you a quasi-static analysis is safe, and which model is closest to violating that?
model_names = ["Tower, lumped-model";"Tower, Ritz phi_t";"Tower, Ritz phi_line"; ... %[output:group:8b16660e] %[output:9a13c586]
               "Tower, coupled (Kb-Kg)";"Disk A, one mode";"Disk B, one mode"]; %[output:group:8b16660e] %[output:9a13c586]
natural_freqs_Hz = [tip_shape_freq_Hz; tip_shape_freq_Hz; ritz_line_freq_Hz; ... %[output:group:4d27a1f4] %[output:1f715cf8]
                    coupled_tower_freq_Hz; dishA_freq_Hz; dishB_freq_Hz]; %[output:group:4d27a1f4] %[output:1f715cf8]
disp(table(model_names,natural_freqs_Hz,natural_freqs_Hz/gust_frequency, ... %[output:group:1daba2d4] %[output:4bbd129d] %[output:36b22933]
    'VariableNames',["Model","f1_Hz","f1_over_f_gust"])); %[output:group:1daba2d4] %[output:4bbd129d] %[output:36b22933]
fprintf('Ruze screening RMS limit at %.0f GHz for eta >= %.2f: %.2f um\n', ... %[output:group:23a894d0] %[output:8f17dd5c] %[output:46a0aa6b]
    rf_frequency_Hz/1e9,rf_min_efficiency,1e6*max_rms_error); %[output:group:23a894d0] %[output:8f17dd5c] %[output:46a0aa6b]
%%
%[text] ## 
%[text] ## Coupled tower-dish dynamics
%[text] **Goal.** Connect the tower and the plate: the plate sits on top of a swaying tower, so the tower's acceleration shakes the plate.
%[text] #### Loads on the dish
%[text] - **Wind.** The wind is horizontal and the dish is tilted $\\epsilon$ above it, so only part of it pushes on the surface: $k\_n=k\_{n0}|\\cos\\epsilon|$.
%[text] - **Drag on the tower.** The dish passes its wind drag to the tower tip: $F\_{dish}=C\_pk\_n\\cos\\epsilon\\,\\pi R^2p(t)$.
%[text] - **Base motion.** The swaying tower shakes the plate, the way a passenger is thrown forward when a car brakes. This adds a load $L\_b\\cos\\epsilon\\,\\ddot q$, with $L\_b=m\_AA\_m$.
%[text] - **Gravity.** $F\_g=-m\_Ag\\sin\\epsilon\\,A\_m$. \
%[text] #### Equations of motion
%[text] Tower:
%[text] $M\_t\\ddot q+C\_t\\dot q+K\_tq=\\left(Q\_tC\_dd\_o+C\_pk\_n\\cos\\epsilon\\,\\pi R^2\\right)p(t)$
%[text] Plate:
%[text] $M\_d\\ddot\\eta+C\_d\\dot\\eta+K\_d\\eta=F\_g-s\\,A\_mC\_pk\_np(t)+s\\,L\_b\\cos\\epsilon\\,\\ddot q$
%[text] Both equations are integrated together with `ode45`, once for each plate case; the case named by `coupled_plate_case` is the one reported in detail below.
%[text] #### Signs to remember
%[text] - $s=+1$: wind hits the concave face and pushes the same way as gravity. This is the worst case.
%[text] - $\\epsilon=90^\\circ$ (zenith): no wind pressure and no base motion on the dish; only gravity acts. \
%[text] #### Options
%[text] - `base_mass_model`: `'ritz_participation'` (default, $L\_b=M/3$ for Case A) or `'modal_mass'` (legacy, $M/5$).
%[text] - `initial_condition`: `'gravity_sag'` (wind switched on at $t=0$) or `'mean_wind'` (starts at equilibrium).
%[text] - `coupled_plate_case`: `'B'` (hub-mounted, the real support) or `'A'`.
%[text] - `dish_wind_model`: `'cos_elevation'` (default) or `'constant'` (legacy).
%[text] - `dish_drag_on_tower`: `true` (default) or `false`. \
%[text] #### What is printed
%[text] - **Start state** at $t=0$: gravity only.
%[text] - **Steady-state peak:** after the start-up transient dies out.
%[text] - **Worst case:** gust peak adding to gravity. \
%[text] All three use flat-plate stiffness. A real curved dish is stiffer, so treat them as upper bounds.
% Quasi-static worst-case pressure at the design elevation: gust peak with the wind ADDING to the weight
design_peak_wind_pressure = dish_pressure_coefficient*design_normal_factor*dynamic_pressure*(1 + gust_fraction); % [Pa] %[output:014a6c03]
design_gravity_pressure   = dish_areal_density*gravity_accel*abs(sind(dish_elevation_deg));                       % [Pa] %[output:8b5c8b95]

% The coupled simulation runs once per plate case. Both cases feed the summary table;
% the case named by coupled_plate_case is also reported in detail.
summary_rms_um = zeros(2,5);   % rows A, B: [steady flat, steady parab, worst flat, worst parab] in um, parab/flat %[output:166ed07f]
for case_idx = 1:2 %[output:group:437ff73c] %[output:9cd37f86]
    % ---- Pick the plate model for this pass of the loop ----
    if case_idx == 1
        plate_case = 'A';                   plate_label      = dishA_label;
        plate_mass = dishA_mass;            plate_damping    = dishA_damping;
        plate_stiffness = dishA_stiffness;  plate_omega_n    = dishA_omega_n;
        plate_freq_Hz = dishA_freq_Hz;      plate_modal_area = dishA_modal_area;
        plate_rms_factor = dishA_rms_factor;
        flat_factor = rmsA_flat_quad;       parab_factor     = rmsA_parab_quad;
    else
        plate_case = 'B';                   plate_label      = dishB_label;
        plate_mass = dishB_mass;            plate_damping    = dishB_damping;
        plate_stiffness = dishB_stiffness;  plate_omega_n    = dishB_omega_n;
        plate_freq_Hz = dishB_freq_Hz;      plate_modal_area = dishB_modal_area;
        plate_rms_factor = dishB_rms_factor;
        flat_factor = rmsB_flat_quad;       parab_factor     = rmsB_parab_quad;
    end

    % ---- Base-participation mass L_b: how strongly the tower acceleration shakes the mode ----
    switch base_mass_model
        case 'ritz_participation', plate_base_mass = dish_areal_density*plate_modal_area;  % m_A * int psi dA
        case 'modal_mass',         plate_base_mass = plate_mass;                           % legacy choice
        otherwise, error('Unknown base_mass_model "%s".',base_mass_model);
    end

    % ---- Time grid resolving the faster of the two modes ----
    coupled_num_samples = max(sim_min_num_samples, ...
        ceil(sim_samples_per_period*max(coupled_tower_freq_Hz,plate_freq_Hz)*sim_duration) + 1);
    coupled_time_grid = linspace(0,sim_duration,coupled_num_samples).';

    % ---- Loads ----
    pressure_factor = dish_pressure_coefficient*design_normal_factor;                    % C_p * k_n
    gravity_force   = -dish_areal_density*gravity_accel*sind(dish_elevation_deg)*plate_modal_area;   % F_g [N]
    tower_drag_area = coupled_tower_force_per_line_load*tower_drag_coefficient*tower_wind_width;     % Q_t*C_d*d_o [m^2]
    dish_drag_area  = dish_drag_on_tower*pressure_factor*abs(cosd(dish_elevation_deg))*dish_flat_area;  % [m^2]
    tower_force_fun      = @(t) (tower_drag_area + dish_drag_area)*wind_pressure_fun(t);      % on the tower [N]
    plate_wind_force_fun = @(t) -wind_direction*plate_modal_area*pressure_factor*wind_pressure_fun(t);  % on the plate [N]
    base_accel_factor    = -wind_direction*cosd(dish_elevation_deg);   % boresight share of the tower acceleration

    % ---- Initial state [q; qdot; eta; etadot] ----
    switch initial_condition
        case 'gravity_sag'      % wind switched on at t = 0, plate already sagging under its weight
            initial_state = [0; 0; gravity_force/plate_stiffness; 0];
        case 'mean_wind'        % both modes start balanced under the mean wind
            initial_state = [tower_force_fun(0)/coupled_tower_stiffness; 0; ...
                             (gravity_force + plate_wind_force_fun(0))/plate_stiffness; 0];
        otherwise
            error('Unknown initial condition "%s".',initial_condition);
    end

    % ---- Equations of motion in first-order form ----
    % Tower:  M_t*q'' = F_t(t) - C_t*q' - K_t*q
    % Plate:  M_d*eta'' = F_g + F_wind(t) - L_b*a_n*q'' - C_d*eta' - K_d*eta
    % The tower does not feel the plate (one-way coupling), so q'' is computed first.
    tower_accel_fun = @(t,state) (tower_force_fun(t) - coupled_tower_damping*state(2) ...
                                  - coupled_tower_stiffness*state(1))/coupled_tower_mass;
    coupled_rhs = @(t,state) [state(2); ...
                              tower_accel_fun(t,state); ...
                              state(4); ...
                              (gravity_force + plate_wind_force_fun(t) ...
                               - plate_base_mass*base_accel_factor*tower_accel_fun(t,state) ...
                               - plate_damping*state(4) - plate_stiffness*state(3))/plate_mass];
    [response_time,response_state] = ode45(coupled_rhs,coupled_time_grid,initial_state, ...
        odeset('RelTol',ode_rel_tol,'AbsTol',ode_abs_tol));
    response_tower_sway       = response_state(:,1);     % q(t)   [m]
    response_plate_deflection = response_state(:,3);     % eta(t) [m]
    response_tower_accel = (tower_force_fun(response_time) - coupled_tower_damping*response_state(:,2) ...
                            - coupled_tower_stiffness*response_tower_sway)/coupled_tower_mass;   % q''(t)

    % ---- Peaks: steady state (after 4 decay time constants of the slower mode) and worst case ----
    settle_time_s = max(4/(dish_damping_ratio*plate_omega_n),4/(tower_damping_ratio*coupled_tower_omega_n));
    is_steady = response_time > settle_time_s;
    if ~any(is_steady), is_steady(:) = true; end
    eta_steady_peak = max(abs(response_plate_deflection(is_steady)));                                  % [m]
    eta_worst_case  = plate_modal_area/plate_stiffness*(design_gravity_pressure + design_peak_wind_pressure);  % [m]
    summary_rms_um(case_idx,:) = [1e6*eta_steady_peak*flat_factor, 1e6*eta_steady_peak*parab_factor, ...
                                  1e6*eta_worst_case*flat_factor,  1e6*eta_worst_case*parab_factor, ...
                                  parab_factor/flat_factor];

    % ---- Keep the detailed results of the selected case ----
    if strcmp(plate_case,coupled_plate_case)
        coupled_plate_label      = plate_label;
        coupled_plate_base_mass  = plate_base_mass;
        coupled_plate_modal_area = plate_modal_area;
        coupled_plate_stiffness  = plate_stiffness;
        coupled_flat_factor      = flat_factor;
        coupled_settle_time      = settle_time_s;
        % Residual surface RMS after removing piston and defocus: sigma(t) = |eta(t)| x residual factor
        coupled_rms_start        = abs(response_plate_deflection(1))*plate_rms_factor;    % [m]
        coupled_rms_steady_peak  = eta_steady_peak*plate_rms_factor;                       % [m]
        coupled_rms_worst_case   = eta_worst_case*plate_rms_factor;                        % [m]
        % Base inertia relative to the full weight, worst orientation (cos(epsilon) = 1)
        coupled_inertia_ratio    = max(abs(plate_base_mass*response_tower_accel))/ ...
                                   (dish_areal_density*gravity_accel*plate_modal_area);
        % Tower-tip slope (pointing error) = phi_t'(L) * q(t)
        coupled_peak_tip_slope   = coupled_tower_tip_slope_per_q*max(abs(response_tower_sway));   % [rad]
    end
end %[output:group:437ff73c]

fprintf(['Residual surface RMS, %s, elevation %g deg (L_b = %.3g kg):\n' ... %[output:group:0bf9644a] %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    '  start state, gravity only (t = 0)          = %7.3f um\n' ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    '  steady-state peak (t > %.2f s)              = %7.3f um\n' ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    '  worst case, gust peak adding to weight      = %7.3f um\n' ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    '  Ruze limit                                  = %7.2f um\n' ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    '  peak base-inertia force / weight (cos = 1)  = %.2e\n' ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    '  peak tower-tip slope (pointing error)       = %.3g arcmin\n'], ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    coupled_plate_label,dish_elevation_deg,coupled_plate_base_mass,1e6*coupled_rms_start,coupled_settle_time, ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    1e6*coupled_rms_steady_peak,1e6*coupled_rms_worst_case, ... %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]
    1e6*max_rms_error,coupled_inertia_ratio,coupled_peak_tip_slope*180*60/pi); %[output:group:0bf9644a] %[output:85dcc2ac] %[output:9d90d5fe] %[output:8c1ddb14] %[output:879659dd] %[output:502f0208] %[output:3a469d4f] %[output:432d1e2d] %[output:9c88c1be] %[output:982e4ef0] %[output:2554db69] %[output:0791430f] %[output:7abdf61b] %[output:5a74e9ce] %[output:0cbc117f] %[output:873516b4] %[output:7d898731]

%%
%[text] ## Surface RMS residual for the true paraboloid shape
%[text] 
%[text] The plate bends by (size `eta`) x (shape `psi`). Moving the feed horn cancels two parts of that for free: **piston**, a uniform shift along the axis, and **defocus**, a change of focal length. Only what is left, the **residual**, costs gain. Picture a straight line fitted through scattered data: the line is what the feed fixes, the scatter is the residual.
%[text] **Derived below with** `sym`**:**
%[text] - Defocus is a straight line in `u = (r/R)^2`, so removing piston and defocus is a straight-line fit in `u`.
%[text] - Equal steps in `u` are equal projected areas, so an area-weighted average is a plain average over `u`.
%[text] - The curved-surface area per projected area, `area_weight_u`, depends only on `u` and `beta`. \
%[text] **Two averaging assumptions:** *flat* (each ring weighted by its projected area) or *true paraboloid* (weighted by its curved area). The deflection shape is the same in both.
%[text] **Recipe:** minimize the weighted mean-square error of `psi - piston - defocus*u` (the normal equations, solved with `solve`). The RMS of what is left is the **residual factor**; multiplied by `eta` it is the residual RMS in metres.
%[text] **Three independent checks:** exact symbolic (Case A), adaptive quadrature (earlier section) and a 4000-ring mesh. Quadrature should match the exact value to about 1e-8, and the mesh should match quadrature to about 1e-5 (midpoint-rule error).
% =====================================================================
% 0. DERIVE THE GEOMETRY USED BY THE FIT (symbolic)
% =====================================================================
syms u real                                   % normalized-area coordinate u = (r/R)^2: 0 at centre, 1 at rim %[output:933e7ed1]
syms focal_len delta_focal real               % focal length of the ideal surface and a small change to it %[output:7cb153b7]
assumeAlso(focal_len > 0); assumeAlso(beta > 0);   % beta = R^2/(4*focal_len^2) is positive for a real dish %[output:678a484b]
paraboloid_height = r^2/(4*focal_len);        % ideal reflector surface: height above the vertex at radius r %[output:9d1e6c88]
r_of_u = R*sqrt(u);                           % radius written in terms of u %[output:53bd2ecb]

% (a) Defocus: first-order change in height when the focal length changes by delta_focal
defocus_height   = diff(paraboloid_height,focal_len)*delta_focal;   % d(height)/d(focal_len) * delta_focal %[output:94e2c576]
defocus_height_u = simplify(subs(defocus_height,r,r_of_u));         % the same change, written in u %[output:1fa97fac]
disp('(a) Height change of the paraboloid for a focal-length change, in terms of u:'); disp(defocus_height_u); %[output:3aa67aa0] %[output:3dead744] %[output:31ba9613]
% A straight line in u through the origin (zero second derivative, zero at u = 0): the basis function u
assert(isAlways(diff(defocus_height_u,u,2) == 0) && isAlways(subs(defocus_height_u,u,0) == 0)); %[output:46acf3ad]

% (b) Projected area per unit u: a ring of radius r and width dr has area 2*pi*r*dr
projected_area_per_u = simplify(2*pi*r_of_u*diff(r_of_u,u)); %[output:13d15e49]
disp('(b) Projected ring area per unit u (a constant, so equal steps in u are equal areas):'); %[output:00567bb7] %[output:781cd027] %[output:37fcbea9] %[output:7f175c70]
disp(projected_area_per_u); %[output:65688e6b] %[output:67210f7d] %[output:4db80dab] %[output:1fb6392a]
assert(isAlways(projected_area_per_u == pi*R^2)); %[output:32e1ae44] %[output:20eeb69d]

% (c) True surface area per projected area: a ring on the curved surface is wider by the
%     length of the meridian per unit radius, sqrt(1 + (d height/dr)^2).
area_stretch_r = sqrt(1 + diff(paraboloid_height,r)^2); %[output:0d9e26e5] %[output:7d894329]
% Write it in u, and replace the focal length by the depth parameter beta = R^2/(4*focal_len^2)
focal_len_of_beta = solve(beta == R^2/(4*focal_len^2),focal_len);   % focal_len in terms of beta (positive root) %[output:588609c8] %[output:9aa30371]
area_weight_u = simplify(subs(area_stretch_r,[r focal_len],[r_of_u focal_len_of_beta])); %[output:70a84856] %[output:589293a0]
disp('(c) True-area weight per projected area, in terms of u and beta:'); disp(area_weight_u); %[output:151ab937] %[output:891b8cfa] %[output:0d2a0f7f] %[output:8967f26e]
assert(isAlways(simplify(area_weight_u^2 - (1 + beta*u)) == 0)); %[output:3ead5be2]

% =====================================================================
% 1. SYMBOLIC least-squares fit, Case A (an exact answer to test the numerics against)
% =====================================================================
caseA_shape_u = (1-u)^2;                     % Case A shape psi_A, written in terms of u %[output:9dd974da]
syms piston defocus real                     % unknown fit coefficients: psi ~ piston + defocus*u %[output:8e06af0b]
dish_beta_exact = sym(dish_depth_beta);      % this dish's beta as an exact fraction (1/4) %[output:1c84ff78]
caseA_rms_exact = sym(zeros(1,2));           % residual factor: [flat weighting, true paraboloid weighting] %[output:1b4ab0fb]
for weight_idx = 1:2 %[output:group:1d06f8ad] %[output:19a14edc]
    if weight_idx == 1
        weight_u = sym(1);                                    % flat: every unit of u weighs the same
    else
        weight_u = subs(area_weight_u,beta,dish_beta_exact);  % true paraboloid area of this dish
    end
    total_weight = int(weight_u,u,0,1);
    % Weighted mean-square error of the fit, as a function of the two coefficients
    mean_square_error = int((caseA_shape_u - piston - defocus*u)^2*weight_u,u,0,1)/total_weight;
    % Best fit: both partial derivatives are zero (the normal equations); solve them
    [piston_best,defocus_best] = solve([diff(mean_square_error,piston) == 0, ...
                                        diff(mean_square_error,defocus) == 0],[piston defocus]);
    % Residual factor = square root of the smallest mean-square error
    caseA_rms_exact(weight_idx) = simplify(sqrt(subs(mean_square_error,[piston defocus],[piston_best defocus_best])));
    if weight_idx == 1                       % keep the flat-weighting fit for display
        caseA_piston_flat  = simplify(piston_best);
        caseA_defocus_flat = simplify(defocus_best);
        caseA_residual_u   = expand(caseA_shape_u - piston_best - defocus_best*u);   % e(u) = what is left over
    end
end %[output:group:1d06f8ad]
caseA_rms_flat_exact = caseA_rms_exact(1); %[output:0655788b] %[output:7d0aa295]
fprintf('Case A best fit (flat weighting): piston = %s, defocus = %s\n', ... %[output:group:7094ff82] %[output:7de62d5c] %[output:4bfc4ea5] %[output:84aeac3c] %[output:8e2a1e71] %[output:15acea20] %[output:1846de84] %[output:650a4d28] %[output:38462457] %[output:0c7e67f1] %[output:25d6776c] %[output:5fe63709] %[output:58cdf1e0]
    char(caseA_piston_flat),char(caseA_defocus_flat)); %[output:group:7094ff82] %[output:7de62d5c] %[output:4bfc4ea5] %[output:84aeac3c] %[output:8e2a1e71] %[output:15acea20] %[output:1846de84] %[output:650a4d28] %[output:38462457] %[output:0c7e67f1] %[output:25d6776c] %[output:5fe63709] %[output:58cdf1e0]
disp('Case A residual e(u):'); disp(caseA_residual_u); %[output:168e6b08] %[output:8b10c884] %[output:9ed670e7] %[output:88859b46]
disp('Case A residual factor, flat-plate assumption (exact):'); disp(caseA_rms_flat_exact); %[output:82db0c2c] %[output:98c1beb2] %[output:00212b2f]
assert(isAlways(caseA_rms_flat_exact == sqrt(sym(5))/30)); %[output:7f5211b5]

% How much of the raw error does refocusing remove? (raw RMS = nothing removed)
caseA_raw_rms_flat     = simplify(sqrt(int(caseA_shape_u^2,u,0,1))); %[output:1e19f30f]
caseA_share_removed    = simplify(1 - caseA_rms_flat_exact/caseA_raw_rms_flat); %[output:2966a5b0]
fprintf('Case A raw RMS factor = %s; refocusing removes %s = %.1f%% of it\n', ... %[output:group:8e6e0618] %[output:8e020f81] %[output:1a0e2275]
    char(caseA_raw_rms_flat),char(caseA_share_removed),100*double(caseA_share_removed)); %[output:group:8e6e0618] %[output:8e020f81] %[output:1a0e2275]

% Where is the residual largest? Compare the ends with the stationary point de/du = 0
caseA_residual_stationary_u = solve(diff(caseA_residual_u,u) == 0,u); %[output:212f101d]
fprintf('Case A residual: e(0) = %s, e(1) = %s, e(%s) = %s (stationary point)\n', ... %[output:group:2c95fb2b] %[output:54371559] %[output:5d1bc094] %[output:1cf9a091] %[output:2b7a8153] %[output:888efc43] %[output:1b7ec761] %[output:7d1ec4e3] %[output:4fade0bb] %[output:9e89eebb] %[output:4093966c] %[output:9aead909] %[output:7bbae824] %[output:45c22e81] %[output:50cc853a] %[output:73d2e459] %[output:2dcaa3cf] %[output:260dbaf7]
    char(subs(caseA_residual_u,u,0)),char(subs(caseA_residual_u,u,1)), ... %[output:54371559] %[output:5d1bc094] %[output:1cf9a091] %[output:2b7a8153] %[output:888efc43] %[output:1b7ec761] %[output:7d1ec4e3] %[output:4fade0bb] %[output:9e89eebb] %[output:4093966c] %[output:9aead909] %[output:7bbae824] %[output:45c22e81] %[output:50cc853a] %[output:73d2e459] %[output:2dcaa3cf] %[output:260dbaf7]
    char(caseA_residual_stationary_u),char(subs(caseA_residual_u,u,caseA_residual_stationary_u))); %[output:group:2c95fb2b] %[output:54371559] %[output:5d1bc094] %[output:1cf9a091] %[output:2b7a8153] %[output:888efc43] %[output:1b7ec761] %[output:7d1ec4e3] %[output:4fade0bb] %[output:9e89eebb] %[output:4093966c] %[output:9aead909] %[output:7bbae824] %[output:45c22e81] %[output:50cc853a] %[output:73d2e459] %[output:2dcaa3cf] %[output:260dbaf7]

% Consistency checks with the earlier closed form (built from the area moments, a formula in beta):
caseA_rms_paraboloid_exact = sqrt(weighted_residual_var); %[output:8670b1e3]
% (i) this section's least-squares derivation gives the same true-area answer at this dish's beta
assert(abs(double(caseA_rms_exact(2))/dishA_rms_factor - 1) < 1e-10, ... %[output:group:5a781273] %[output:67ee1a3d]
    'The least-squares derivation disagrees with the area-moment closed form.'); %[output:group:5a781273] %[output:67ee1a3d]
% (ii) a very shallow paraboloid (beta -> 0) is a flat plate, so the formula must give the flat answer
assert(abs(double(vpa(subs(caseA_rms_paraboloid_exact,beta,sym(1)/1000),50)) - double(caseA_rms_flat_exact)) < 1e-6, ... %[output:group:13bc03c2] %[output:20218990]
    'The paraboloid formula does not reduce to the flat one as beta -> 0.'); %[output:group:13bc03c2] %[output:20218990]

% =====================================================================
% 2. NUMERICAL, adaptive quadrature: already computed in the "RMS over a true paraboloid
%    surface" section, stored in surface_rms_factors and rmsA/B_flat/parab_quad.
% =====================================================================

% =====================================================================
% 3. NUMERICAL, surface mesh: thin rings, each weighted by its area
% =====================================================================
% Turn the derived true-area weight into a numeric function of u for this dish
area_weight_fun = matlabFunction(subs(area_weight_u,beta,dish_beta_exact),'Vars',u); %[output:50d9dd25]
% The four reflector surface models, in the row order used by every table and plot below
surface_model_names = ["Case A, flat";"Case A, paraboloid";"Case B, flat";"Case B, paraboloid"]; %[output:0265610e]
surface_rms_mesh    = zeros(4,1);                   % residual factor of each surface model %[output:14790efa]
ring_radius_all     = zeros(rms_mesh_rings,4);      % ring-centre radius of each model [m] (kept for the plot) %[output:9add2ad8]
residual_shape_all  = zeros(rms_mesh_rings,4);      % residual per unit eta (kept for the plot) %[output:16ad56fb]

for model_idx = 1:4 %[output:group:208438cf] %[output:50f34e01]
    % ---- Which support case and which weighting is this surface model? ----
    is_case_A     = model_idx <= 2;                 % rows 1-2: Case A, rows 3-4: Case B
    is_paraboloid = mod(model_idx,2) == 0;          % rows 2 and 4 use the true paraboloid area
    if is_case_A
        shape_fun = dishA_shape_fun;  inner_radius = 0;               % solid disk: rings start at the centre
    else
        shape_fun = dishB_shape_fun;  inner_radius = dish_hub_radius; % annulus: rings start at the hub
    end

    % ---- Build the ring mesh ----
    ring_edges   = linspace(inner_radius,dish_radius,rms_mesh_rings + 1).';    % ring boundaries [m]
    ring_centres = (ring_edges(1:end-1) + ring_edges(2:end))/2;                % ring mid-radii [m]
    ring_widths  = diff(ring_edges);                                           % dr [m]
    ring_u       = (ring_centres/dish_radius).^2;                              % u at each ring
    % Ring area = circumference x width x area weight (1 for flat, derived weight for the paraboloid) [m^2]
    if is_paraboloid
        ring_area = 2*pi*ring_centres.*ring_widths.*area_weight_fun(ring_u);
    else
        ring_area = 2*pi*ring_centres.*ring_widths;
    end
    ring_shape   = shape_fun(ring_centres);                                    % psi at each ring

    % ---- Weighted straight-line fit psi ~ piston + defocus*u ----
    % Minimizing sum(area .* (psi - piston - defocus*u).^2) gives the same 2x2 normal equations
    % as the symbolic derivation:  (X' W X) c = X' W psi,  with X = [1, u] and W = diag(area).
    fit_basis  = [ones(size(ring_u)) ring_u];                  % column 1: piston, column 2: defocus
    fit_coeffs = (fit_basis.'*(ring_area.*fit_basis)) \ (fit_basis.'*(ring_area.*ring_shape));   % [piston; defocus]

    % ---- Residual and its area-weighted RMS ----
    fit_residual = ring_shape - fit_basis*fit_coeffs;          % what refocusing cannot remove
    surface_rms_mesh(model_idx) = sqrt(sum(ring_area.*fit_residual.^2)/sum(ring_area));   % residual factor

    ring_radius_all(:,model_idx)    = ring_centres;
    residual_shape_all(:,model_idx) = fit_residual;
end %[output:group:208438cf]
rmsA_flat_mesh = surface_rms_mesh(1);   rmsA_parab_mesh = surface_rms_mesh(2); %[output:2f0325b5]
rmsB_flat_mesh = surface_rms_mesh(3);   rmsB_parab_mesh = surface_rms_mesh(4); %[output:87e0d49f]

% ---- Do the three methods agree? ----
assert(abs(rmsA_flat_quad/double(caseA_rms_flat_exact) - 1) < 1e-8,'Flat quadrature disagrees with the exact value.'); %[output:685295c2]
assert(abs(rmsA_parab_quad/double(caseA_rms_exact(2)) - 1) < 1e-8,'Paraboloid quadrature disagrees with the exact value.'); %[output:826003bb]
assert(all(abs(surface_rms_mesh./surface_rms_factors(:,3) - 1) < 1e-5),'Mesh and quadrature disagree.'); %[output:882420af]
disp(array2table([rmsA_flat_quad rmsA_parab_quad rmsA_flat_mesh rmsA_parab_mesh; ... %[output:group:8e539061] %[output:65f0dab6] %[output:32af5974]
                  rmsB_flat_quad rmsB_parab_quad rmsB_flat_mesh rmsB_parab_mesh], ... %[output:65f0dab6] %[output:32af5974]
    'RowNames',["Case A","Case B"],'VariableNames', ... %[output:65f0dab6] %[output:32af5974]
    ["Flat_quadrature","Paraboloid_quadrature","Flat_mesh","Paraboloid_mesh"])); %[output:group:8e539061] %[output:65f0dab6] %[output:32af5974]
%[text] ### Surface RMS of each reflector surface model, in micrometres
%[text] **Goal.** Turn the residual factors (pure numbers) into physical surface errors for the four reflector surface models, and see how much of each model's error the receiver can correct.
%[text] For each model the table shows three numbers, each removing more of the error than the one before:
%[text] - **Raw RMS:** the RMS of the deflection with nothing removed.
%[text] - **Piston removed:** after cancelling the uniform shift along the axis.
%[text] - **Residual (piston + defocus removed):** after refocusing the feed. This is the only one that enters Ruze's equation, so it is the one compared with the Ruze limit. \
%[text] To get metres, each factor is multiplied by the bending size. The code derives the **worst-case bending size** `eta_worst_sym` from the static balance of the one-mode plate (stiffness force equals the generalized pressure force), with the largest normal pressure: the weight plus the gust-peak wind pushing the same way. The Residual column therefore matches the *WorstCase* columns of the summary table below.
%[text] The last column is the Ruze aperture efficiency of the residual RMS. The code also re-derives the Ruze RMS limit by solving the Ruze equation for the RMS at the minimum efficiency.
%[text] **Bar chart: what the raw RMS is made of.** Because the best-fit piston and defocus are found by least squares, the three parts of the error are independent and add *in quadrature*: raw RMS squared = piston squared + defocus squared + residual squared. The left plot shows each part in micrometres next to the raw RMS; the right plot shows each part's share of the raw mean-square error, so the bars always add up to 100%. Only the residual part is compared with the Ruze limit.
% ---- Worst-case bending size, derived from the static one-mode equation (symbolic) ----
syms eta_amp p_normal real                                  % bending size and uniform normal pressure %[output:998ba623]
syms K_d A_m m_A C_p k_n p_inf g_f sin_elev_abs real        % plate stiffness, modal area, areal density, wind terms %[output:50cc1cf9]
assumeAlso([K_d A_m m_A C_p p_inf] > 0); %[output:29d5b5ef] %[output:6e793794]
assumeAlso([k_n g_f sin_elev_abs] >= 0); %[output:8c14aaeb] %[output:85a027cd]
% Static balance of one mode: stiffness force K_d*eta = generalized force A_m*p
eta_static_sym = solve(K_d*eta_amp == A_m*p_normal,eta_amp); %[output:64c5f01c] %[output:69257ae9]
% Largest normal pressure: weight m_A*g*|sin(elevation)| plus gust-peak wind C_p*k_n*p_inf*(1 + g_f)
peak_pressure_sym = m_A*g*sin_elev_abs + C_p*k_n*p_inf*(1 + g_f); %[output:96d870bc] %[output:68069839]
eta_worst_sym = subs(eta_static_sym,p_normal,peak_pressure_sym); %[output:929141a4] %[output:20650835]
disp('Worst-case bending size eta_worst:'); disp(eta_worst_sym); %[output:642ea8d6] %[output:373ca1fa] %[output:18602282] %[output:5b0c1c55] %[output:1e4b4919] %[output:9f2e5f2e]

% ---- Numbers: the design condition is the same for both plates ----
design_symbols = [m_A g sin_elev_abs C_p k_n p_inf g_f]; %[output:4c8ead0c] %[output:4ca67804]
design_values  = [dish_areal_density gravity_accel abs(sind(dish_elevation_deg)) ... %[output:group:96baed3a] %[output:071c7e02]
                  dish_pressure_coefficient design_normal_factor dynamic_pressure gust_fraction]; %[output:group:96baed3a] %[output:071c7e02]
dishA_worst_amplitude = double(subs(eta_worst_sym,[K_d A_m design_symbols],[dishA_stiffness dishA_modal_area design_values])); %[output:351156b1]
dishB_worst_amplitude = double(subs(eta_worst_sym,[K_d A_m design_symbols],[dishB_stiffness dishB_modal_area design_values])); %[output:610f613a]
% Check: same as the pressures used in the coupled-simulation loop
assert(abs(dishA_worst_amplitude*dishA_stiffness/dishA_modal_area/(design_gravity_pressure + design_peak_wind_pressure) - 1) < 1e-12); %[output:14f200be]
surface_amplitudes = [dishA_worst_amplitude; dishA_worst_amplitude; ...          % one amplitude per surface model %[output:group:87cd9a50] %[output:720065c8]
                      dishB_worst_amplitude; dishB_worst_amplitude];             % (rows as surface_model_names) %[output:group:87cd9a50] %[output:720065c8]

% ---- RMS in micrometres = factor x amplitude ----
% surface_rms_factors columns: [raw, piston removed, residual]; one row per surface model
surface_rms_um = 1e6*surface_rms_factors.*surface_amplitudes;                    % [um] %[output:092c5c33]

% ---- Ruze equation (symbolic): efficiency of a surface with RMS error sigma at wavelength lambda ----
syms sigma_rms lambda_rf eta_min real %[output:422bbed5]
assumeAlso(sigma_rms > 0); assumeAlso(lambda_rf > 0); assumeAlso(eta_min > 0 & eta_min < 1); %[output:321eb212]
ruze_efficiency_sym = exp(-(4*pi*sigma_rms/lambda_rf)^2); %[output:258db12e]
% Largest allowed RMS: solve efficiency = eta_min for sigma
sigma_limit_sym = solve(ruze_efficiency_sym == eta_min,sigma_rms); %[output:42189739]
disp('Ruze RMS limit, solved from efficiency = eta_min:'); disp(sigma_limit_sym); %[output:9a243612] %[output:67369cb3] %[output:394e828f] %[output:63454957]
assert(abs(double(subs(sigma_limit_sym,[lambda_rf eta_min],[rf_wavelength rf_min_efficiency]))/max_rms_error - 1) < 1e-12); %[output:2bfd36a8] %[output:2361add0]
% Efficiency of each surface model's residual RMS
ruze_efficiency_fun = matlabFunction(subs(ruze_efficiency_sym,lambda_rf,rf_wavelength),'Vars',sigma_rms); %[output:722846cd] %[output:9e68fee2]
surface_ruze_efficiency = ruze_efficiency_fun(1e-6*surface_rms_um(:,3)); %[output:9e8bcdfe] %[output:518ff14f]

fprintf('Surface RMS per reflector surface model, elevation %g deg, wind %g m/s, gust peak (Ruze limit %.2f um):\n', ... %[output:group:57142bb0] %[output:3d26dac0] %[output:0826dd97] %[output:147cd0a9] %[output:8ef7b83e] %[output:50f574b8] %[output:11e576da] %[output:08e109e1] %[output:77f1b85a] %[output:41d8c077] %[output:56f7cc79]
    dish_elevation_deg,wind_speed,1e6*max_rms_error); %[output:group:57142bb0] %[output:3d26dac0] %[output:0826dd97] %[output:147cd0a9] %[output:8ef7b83e] %[output:50f574b8] %[output:11e576da] %[output:08e109e1] %[output:77f1b85a] %[output:41d8c077] %[output:56f7cc79]
disp(table(surface_rms_um(:,1),surface_rms_um(:,2),surface_rms_um(:,3),surface_ruze_efficiency, ... %[output:group:92531640] %[output:1974face] %[output:8d8d80f2]
    'RowNames',surface_model_names, ... %[output:1974face] %[output:8d8d80f2]
    'VariableNames',["Raw_um","PistonRemoved_um","Residual_um","RuzeEfficiency"])); %[output:group:92531640] %[output:1974face] %[output:8d8d80f2]

% ---- Break the raw RMS of each surface model into its three independent parts [um] ----
% surface_rms_um columns: [raw, piston removed, residual]. The least-squares fit makes the parts
% orthogonal, so they add in quadrature:  raw^2 = piston^2 + defocus^2 + residual^2.
component_raw_um      = surface_rms_um(:,1); %[output:24215b1e]
component_piston_um   = sqrt(surface_rms_um(:,1).^2 - surface_rms_um(:,2).^2);   % removed by moving the feed axially %[output:6fd05a4d]
component_defocus_um  = sqrt(surface_rms_um(:,2).^2 - surface_rms_um(:,3).^2);   % removed by refocusing %[output:523a015f]
component_residual_um = surface_rms_um(:,3);                                      % what refocusing cannot remove %[output:8a73abdf]
component_rms_um = [component_piston_um component_defocus_um component_residual_um]; %[output:7ebf4fc2]
% Check: the parts rebuild the raw RMS exactly
assert(all(abs(sqrt(sum(component_rms_um.^2,2))./component_raw_um - 1) < 1e-10),'RMS components do not add up.'); %[output:3d19546e]
% Share of the raw mean-square error carried by each part [%]
component_share_pct = 100*component_rms_um.^2./component_raw_um.^2; %[output:98adf10a]

disp(table(component_raw_um,component_piston_um,component_defocus_um,component_residual_um,component_share_pct(:,3), ... %[output:group:26f1702f] %[output:610aedc7] %[output:381e3a65]
    'RowNames',surface_model_names, ... %[output:610aedc7] %[output:381e3a65]
    'VariableNames',["Raw_um","Piston_um","Defocus_um","Residual_um","Residual_share_pct"])); %[output:group:26f1702f] %[output:610aedc7] %[output:381e3a65]

% ---- Figure: bar chart of the RMS breakdown for each surface model ----
figure('Color','w','Name','Breakdown of raw surface RMS into piston, defocus and residual'); %[output:5f2ab6c1]
tiledlayout(1,2,'TileSpacing','compact'); %[output:3f87dc83]

% Left: each part in micrometres, next to the raw RMS (grouped bars, one group per surface model)
nexttile; %[output:9ad12b24]
bar([component_raw_um component_rms_um]); %[output:8045a21c]
hold on; %[output:780e45a5] %[output:6409e4a4] %[output:0e251706] %[output:8f6e4863] %[output:4e6b99f6] %[output:85702a69] %[output:3e2fa706] %[output:309bcd49]
yline(1e6*max_rms_error,'--k','Ruze limit (residual only)','LabelHorizontalAlignment','left'); %[output:0cc2cc3b] %[output:63cc2df3] %[output:1fd339d1]
set(gca,'XTickLabel',surface_model_names); %[output:4b20fe49]
grid on; ylabel('RMS [μm]'); %[output:2c6997e9]
title('Worst-case RMS and its parts'); %[output:40107e6f]
legend('Raw (total)','Piston','Defocus','Residual','Location','northwest'); %[output:21db08d4]

% Right: each part's share of the raw mean-square error (stacked bars, each adds up to 100%)
nexttile; %[output:07343fc3]
bar(component_share_pct,'stacked'); %[output:81b5ba64] %[output:965451e5] %[output:76ef2bf6] %[output:4eb96057] %[output:21dbdaed] %[output:1f9671ee] %[output:3137f95b] %[output:09e23664] %[output:07e3db6c] %[output:14695a36]
set(gca,'XTickLabel',surface_model_names); %[output:27711343]
grid on; ylim([0 100]); ylabel('Share of raw mean-square error [%]'); %[output:0913fea3]
title('What the raw RMS is made of'); %[output:12c12072]
legend('Piston','Defocus','Residual','Location','southoutside','Orientation','horizontal'); %[output:56543b81] %[output:4b2aa982] %[output:69089d7a] %[output:990a597f] %[output:0b423e72] %[output:0788f97d] %[output:89b1f39e] %[output:7db6b46f] %[output:36e63e55] %[output:8f3e015a]

%[text] 
%%
%[text] ## Summary: maximum residual RMS, and how it grows with wind
%[text] ### Flat-plate vs true paraboloid assumption
%[text] **Goal.** Put the largest residual RMS of each plate model side by side, under both averaging assumptions.
%[text] Two "maximum" values are reported for each case:
%[text] - **Steady-state peak:** the largest residual RMS in the coupled simulation after the start-up wobble has died out.
%[text] - **Worst case:** a quick hand-calculation with the gust at its peak and the wind pushing the same way as gravity. The code builds it symbolically as (residual factor) times the worst-case bending size derived above, and checks it against the numbers from the coupled-simulation loop. \
% ---- Worst-case residual RMS (symbolic): residual factor x worst-case bending size ----
syms residual_factor real %[output:1463273b]
assumeAlso(residual_factor > 0); %[output:9c78f762] %[output:91c7f05f]
sigma_worst_sym = residual_factor*eta_worst_sym; %[output:6dc10915] %[output:9b4ebb0c]
disp('Worst-case residual RMS sigma_worst:'); disp(sigma_worst_sym); %[output:080920b8] %[output:3c8fe4a4] %[output:22317ff6] %[output:0fa313d9] %[output:6d8360fe] %[output:0ebf9f35]
% Check: the symbolic formula reproduces the flat-weighting worst cases filled in during the coupled loop
worst_check_um = 1e6*[ ... %[output:group:401a5efc] %[output:98ba01c2] %[output:3e98e5aa] %[output:6985a994] %[output:768f3467]
    double(subs(sigma_worst_sym,[residual_factor K_d A_m design_symbols],[rmsA_flat_quad dishA_stiffness dishA_modal_area design_values])); ... %[output:98ba01c2] %[output:3e98e5aa] %[output:6985a994] %[output:768f3467]
    double(subs(sigma_worst_sym,[residual_factor K_d A_m design_symbols],[rmsB_flat_quad dishB_stiffness dishB_modal_area design_values]))]; %[output:group:401a5efc] %[output:98ba01c2] %[output:3e98e5aa] %[output:6985a994] %[output:768f3467]
assert(all(abs(worst_check_um./summary_rms_um(:,3) - 1) < 1e-10),'Symbolic worst case disagrees with the coupled-loop numbers.'); %[output:29fc20f7] %[output:92397eac]

% ---- Table: peak and steady-state residual RMS of every surface model ----
% summary_rms_um (filled in the coupled-simulation loop) has columns:
%   [steady flat, steady paraboloid, worst flat, worst paraboloid, paraboloid/flat]
% Steady state = largest residual RMS in the coupled run after the start-up transient has died out.
% Peak         = worst case: gust peak with the wind pushing the same way as the weight.
rms_model_names = ["Case A, flat plate";"Case A, true paraboloid"; ... %[output:group:84598411] %[output:4dc70225]
                   "Case B, flat plate";"Case B, true paraboloid"]; %[output:group:84598411] %[output:4dc70225]
rms_steady_um   = [summary_rms_um(1,1); summary_rms_um(1,2); summary_rms_um(2,1); summary_rms_um(2,2)]; %[output:4e80adac]
rms_peak_um     = [summary_rms_um(1,3); summary_rms_um(1,4); summary_rms_um(2,3); summary_rms_um(2,4)]; %[output:0e933e29]
% Change from the flat-plate value of the same support case [%]
flat_peak_um    = [summary_rms_um(1,3); summary_rms_um(1,3); summary_rms_um(2,3); summary_rms_um(2,3)]; %[output:853287b0]
rms_change_pct  = 100*(rms_peak_um./flat_peak_um - 1); %[output:964f9592]
% Ruze efficiency at the peak, and whether the peak stays within the limit
rms_peak_efficiency = ruze_efficiency_fun(1e-6*rms_peak_um); %[output:9fce0104]

fprintf('Residual RMS at elevation %g deg, wind %g m/s (Ruze limit %.2f um):\n', ... %[output:group:11c97fff] %[output:33c7f5f5] %[output:203c694d]
    dish_elevation_deg,wind_speed,1e6*max_rms_error); %[output:group:11c97fff] %[output:33c7f5f5] %[output:203c694d]
disp(table(rms_steady_um,rms_peak_um,rms_change_pct,rms_peak_efficiency,rms_peak_um <= 1e6*max_rms_error, ... %[output:group:68c0e92a] %[output:7234077c] %[output:2578dbec]
    'RowNames',rms_model_names, ... %[output:7234077c] %[output:2578dbec]
    'VariableNames',["SteadyState_um","Peak_um","Change_vs_flat_pct","RuzeEfficiency_at_peak","Within_Ruze_limit"])); %[output:group:68c0e92a] %[output:7234077c] %[output:2578dbec]
%[text] 
%[text] 
%[text] ### Residual RMS versus wind speed (flat plate and true paraboloid)
%[text] **Goal.** Find out how strong the wind can get before the surface error reaches the Ruze limit.
%[text] 
%[text] At zenith the wind does not press on the dish in this model (the normal-pressure factor is zero), so the curve is flat and the critical wind speed is infinite.
% ---- Worst-case residual RMS as a function of wind speed (symbolic) ----
syms wind_V rho_air real %[output:7690d45b]
assumeAlso(wind_V >= 0); assumeAlso(rho_air > 0); %[output:16dfc757]
sigma_worst_V_sym = subs(sigma_worst_sym,p_inf,rho_air*wind_V^2/2);    % dynamic pressure = rho*V^2/2 %[output:2c75789c]
disp('Worst-case residual RMS versus wind speed V:'); disp(sigma_worst_V_sym); %[output:669fbf34] %[output:2da62fe4] %[output:3d8db93a] %[output:869e28a2]

% Gravity part (no wind) and wind part
gravity_part_sym = simplify(subs(sigma_worst_V_sym,wind_V,0)); %[output:87960db1] %[output:518e2c61]
wind_part_sym    = simplify(sigma_worst_V_sym - gravity_part_sym); %[output:41ee0220] %[output:9268efeb]
disp('Gravity part (V = 0):'); disp(gravity_part_sym); %[output:485ef88a] %[output:7bb895bf] %[output:7d18e504] %[output:098c5a17] %[output:8e6db08f] %[output:6535cc22]
disp('Wind part:');            disp(wind_part_sym); %[output:44cc79d2] %[output:338b4cc8] %[output:7587017e] %[output:29e87d51] %[output:57f6b99a] %[output:270d2ed2]
doubling_factor = simplify(subs(wind_part_sym,wind_V,2*wind_V)/wind_part_sym); %[output:9aa94e7b]
fprintf('Doubling the wind speed multiplies the wind part by %s\n',char(doubling_factor)); %[output:1771c7d4] %[output:1882bb95]

% Gust equivalence: the steady speed whose pressure equals the gust-peak pressure of a gusting wind
syms gust_mean_V steady_V positive %[output:2f6d9085]
equivalent_speed_ratio = simplify(solve((1 + g_f)*gust_mean_V^2 == steady_V^2,steady_V)/gust_mean_V); %[output:89effa97]
fprintf('Equivalent steady wind = %s x gust mean speed = %.3f x for g_f = %.2f\n', ... %[output:group:5c7af4b2] %[output:3b132792] %[output:3d9d8beb]
    char(equivalent_speed_ratio),double(subs(equivalent_speed_ratio,g_f,gust_fraction)),gust_fraction); %[output:group:5c7af4b2] %[output:3b132792] %[output:3d9d8beb]

% Critical wind speed: solve sigma_worst(V) = sigma_max. The equation is linear in V^2,
% so solve for V_squared and take the square root.
syms sigma_max_sym V_squared real %[output:98607d6d]
assumeAlso(sigma_max_sym > 0); %[output:588540fc]
V_squared_sol = solve(subs(sigma_worst_V_sym,wind_V,sqrt(V_squared)) == sigma_max_sym,V_squared); %[output:7185a8fd]
critical_speed_sym = sqrt(V_squared_sol); %[output:17c5017d] %[output:70ae8c77]
disp('Critical wind speed V_crit:'); disp(critical_speed_sym); %[output:11b255a9] %[output:222c4708] %[output:2bada1e5] %[output:082e5e86] %[output:1e5ee7de] %[output:02d6fb1a]

% ---- Numbers: evaluate the symbolic results for the coupled plate at each elevation ----
sweep_speeds    = wind_speed_sweep(:);                                         % wind speeds [m/s] %[output:67dd7d90] %[output:1e6474ea]
sweep_symbols   = [residual_factor K_d A_m m_A g sin_elev_abs C_p k_n g_f rho_air]; %[output:8dbe2bf8] %[output:8b357688]
sweep_rms_worst = zeros(numel(sweep_speeds),numel(plot_elevations_deg));      % [m], rows: speed, columns: elevation %[output:075409db] %[output:0e6f3390]
critical_speed  = zeros(1,numel(plot_elevations_deg));                        % V_crit [m/s] per elevation %[output:794085c8] %[output:3fda3d1c]
for elev_idx = 1:numel(plot_elevations_deg) %[output:group:36e8bf06] %[output:7aa8c2f4] %[output:1efcdc49]
    elev_deg = plot_elevations_deg(elev_idx);
    sweep_values = [coupled_flat_factor coupled_plate_stiffness coupled_plate_modal_area dish_areal_density ...
                    gravity_accel abs(sind(elev_deg)) dish_pressure_coefficient compare_normal_factor(elev_idx) ...
                    gust_fraction air_density];
    % RMS curve: substitute this elevation's numbers, then evaluate at every wind speed
    sigma_of_V = subs(sigma_worst_V_sym,sweep_symbols,sweep_values);           % depends only on wind_V now
    sweep_rms_worst(:,elev_idx) = double(subs(sigma_of_V,wind_V,sweep_speeds)).*ones(size(sweep_speeds));
    % Critical speed
    if compare_normal_factor(elev_idx) == 0
        critical_speed(elev_idx) = Inf;                                        % wind cannot reach the dish
    else
        V_squared_value = double(subs(V_squared_sol,[sweep_symbols sigma_max_sym],[sweep_values max_rms_error]));
        if V_squared_value < 0
            critical_speed(elev_idx) = 0;                                      % gravity alone already too much
        else
            critical_speed(elev_idx) = sqrt(V_squared_value);
        end
    end
end %[output:group:36e8bf06]
disp(table(plot_elevations_deg(:),1e6*sweep_rms_worst(1,:).',1e6*interp1(sweep_speeds,sweep_rms_worst,wind_speed).',critical_speed(:), ... %[output:group:5e94f3c2] %[output:5d81d09c] %[output:005a3d41] %[output:536ddcb5] %[output:6e01a8ca] %[output:6b0e9ce6] %[output:6406317a] %[output:4769d9db] %[output:6e14cf74] %[output:64a526e5] %[output:85437dc7] %[output:85e31600] %[output:7ccd7db4] %[output:5bb2b523]
    'VariableNames',["Elevation_deg","RMS_at_0ms_um","RMS_at_design_wind_um","Vcrit_ms"])); %[output:group:5e94f3c2] %[output:5d81d09c] %[output:005a3d41] %[output:536ddcb5] %[output:6e01a8ca] %[output:6b0e9ce6] %[output:6406317a] %[output:4769d9db] %[output:6e14cf74] %[output:64a526e5] %[output:85437dc7] %[output:85e31600] %[output:7ccd7db4] %[output:5bb2b523]

% ---- True-paraboloid curves: the model is linear in pressure, so they are the flat-plate
%      curves scaled by the ratio of the two residual factors of the coupled plate ----
switch coupled_plate_case %[output:group:3697050a] %[output:93b82c65] %[output:105a419c]
    case 'A', coupled_parab_factor = rmsA_parab_quad; %[output:1b3d5686]
    case 'B', coupled_parab_factor = rmsB_parab_quad;
end %[output:group:3697050a]
sweep_rms_worst_parab = sweep_rms_worst*(coupled_parab_factor/coupled_flat_factor);   % [m] %[output:979f87b9]
critical_speed_parab  = critical_speed;           % V_crit with paraboloid weighting [m/s] %[output:43eb49c2]
for elev_idx = 1:numel(plot_elevations_deg) %[output:group:437aeeb8] %[output:5c45bbe8]
    if isfinite(critical_speed(elev_idx)) && critical_speed(elev_idx) > 0
        % Same equation as above with the paraboloid residual factor
        sweep_values = [coupled_parab_factor coupled_plate_stiffness coupled_plate_modal_area dish_areal_density ...
                        gravity_accel abs(sind(plot_elevations_deg(elev_idx))) dish_pressure_coefficient ...
                        compare_normal_factor(elev_idx) gust_fraction air_density];
        V_squared_value = double(subs(V_squared_sol,[sweep_symbols sigma_max_sym],[sweep_values max_rms_error]));
        critical_speed_parab(elev_idx) = sqrt(max(V_squared_value,0));
    end
end %[output:group:437aeeb8]
disp(table(plot_elevations_deg(:),critical_speed(:),critical_speed_parab(:), ... %[output:group:2d41dbb7] %[output:3ff92b1b] %[output:8b2af639]
    'VariableNames',["Elevation_deg","Vcrit_flat_ms","Vcrit_paraboloid_ms"])); %[output:group:2d41dbb7] %[output:3ff92b1b] %[output:8b2af639]

% ---- Figure: residual RMS versus wind speed ----
figure('Color','w','Name',sprintf('Residual RMS versus wind speed, %s',coupled_plate_label)); %[output:9e4fd169]
hold on; %[output:673cc0e3] %[output:62233bd6] %[output:79ef6561] %[output:44e1d1c0] %[output:40fe681c] %[output:5e8d7843] %[output:52874a2f] %[output:7f342048]
line_colors = lines(numel(plot_elevations_deg));     % one colour per elevation %[output:188c1e91] %[output:06db1b34]
for elev_idx = 1:numel(plot_elevations_deg) %[output:group:4faa3fdc] %[output:0ad942b8] %[output:9d614c91]
    plot(sweep_speeds,1e6*sweep_rms_worst(:,elev_idx),'-','Color',line_colors(elev_idx,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('\\epsilon = %g^\\circ, flat plate',plot_elevations_deg(elev_idx)));
    plot(sweep_speeds,1e6*sweep_rms_worst_parab(:,elev_idx),'--','Color',line_colors(elev_idx,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('\\epsilon = %g^\\circ, true paraboloid',plot_elevations_deg(elev_idx)));
    % Mark where the flat-plate curve reaches the Ruze limit
    if isfinite(critical_speed(elev_idx)) && critical_speed(elev_idx) > 0 && critical_speed(elev_idx) <= max(sweep_speeds)
        plot(critical_speed(elev_idx),1e6*max_rms_error,'o','Color',line_colors(elev_idx,:), ...
            'MarkerFaceColor',line_colors(elev_idx,:),'HandleVisibility','off');
    end
end %[output:group:4faa3fdc]
yline(1e6*max_rms_error,'--k','Ruze limit','HandleVisibility','off'); %[output:4b437120]
xline(wind_speed,':k','design wind','HandleVisibility','off'); %[output:60948342]
grid on; xlabel('Wind speed V [m/s]'); ylabel('Worst-case residual RMS \sigma_{worst} [μm]'); %[output:97872dfe]
title(sprintf('Residual RMS versus wind speed, %s (dots: critical wind speed)',coupled_plate_label)); %[output:72979fed]
legend('Location','northwest'); %[output:511f9c5e] %[output:6b87809c] %[output:320b45b1] %[output:0422c1de] %[output:26d1eee1] %[output:90bafddd] %[output:74f09957] %[output:5657a1a5] %[output:5b02579a] %[output:664e5bf1] %[output:302186c2] %[output:6b3a8940] %[output:9e46565a]

%[text] 
%%
%[text] ## Local functions
%[text] Only two local functions remain: `make_sdof` (natural frequency and damping coefficient of one mode) and `simulate_sdof` (time integration of one mode with `ode45`). Everything else is written in line where it is used. MATLAB requires local functions in a script to appear after all other code.
% ================== Single-degree-of-freedom building blocks =====================
function [omega_n,freq_Hz,damping] = make_sdof(mass,stiffness,damping_ratio)
% One vibration mode from its generalized mass m and stiffness k:
%   natural frequency omega_n = sqrt(k/m) [rad/s], f_n = omega_n/(2*pi) [Hz],
%   viscous damping coefficient c = 2*zeta*m*omega_n [N s/m] for the damping ratio zeta.
    omega_n = sqrt(stiffness/mass);
    freq_Hz = omega_n/(2*pi);
    damping = 2*damping_ratio*mass*omega_n;
end

function [time_out,displacement,velocity] = simulate_sdof(mass,damping,stiffness,time_vec, ...
        force_fun,initial_state,rel_tol,abs_tol)
% Solve m*x'' + c*x' + k*x = F(t) with ode45 and return x(t) and x'(t) on time_vec.
% ode45 needs first-order form, so the state is [x; x'] and its derivative is [x'; x''].
% force_fun must accept a vector t (it is also evaluated on whole time vectors elsewhere).
    rhs = @(t,state) [state(2); (force_fun(t) - damping*state(2) - stiffness*state(1))/mass];
    [time_out,state_out] = ode45(rhs,time_vec,initial_state,odeset('RelTol',rel_tol,'AbsTol',abs_tol));
    displacement = state_out(:,1);
    velocity     = state_out(:,2);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:392983c0]
%   data: {"dataType":"text","outputData":{"text":"    <strong>TubeArea_m2<\/strong>     <strong>TubeI_m4<\/strong>     <strong>LineMass_kg_m<\/strong>    <strong>PlateD_Nm<\/strong>    <strong>DynPressure_Pa<\/strong>    <strong>RuzeLimit_um<\/strong>\n    <strong>___________<\/strong>    <strong>__________<\/strong>    <strong>_____________<\/strong>    <strong>_________<\/strong>    <strong>______________<\/strong>    <strong>____________<\/strong>\n\n     0.023562      7.3631e-05       63.617          11113          22.05            75.782   \n\n","truncated":false}}
%---
%[output:39226161]
%   data: {"dataType":"text","outputData":{"text":"Integration constants [C1 C2 C3 C4]:\n","truncated":false}}
%---
%[output:8337e5bb]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{cccc}\n0 & 0 & \\frac{F\\,L}{2\\,E_t \\,I_t } & -\\frac{F}{6\\,E_t \\,I_t }\n\\end{array}\\right)"}}
%---
%[output:42c5d1f8]
%   data: {"dataType":"text","outputData":{"text":"Deflection y(z):\n","truncated":false}}
%---
%[output:534651bc]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{F\\,z^2 \\,{\\left(3\\,L-z\\right)}}{6\\,E_t \\,I_t }"}}
%---
%[output:01d3e9f7]
%   data: {"dataType":"text","outputData":{"text":"Tip deflection delta = y(L):\n","truncated":false}}
%---
%[output:401ac056]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{F\\,L^3 }{3\\,E_t \\,I_t }"}}
%---
%[output:119e2848]
%   data: {"dataType":"text","outputData":{"text":"Normalized shape phi_t(z) = y\/delta:\n","truncated":false}}
%---
%[output:7a20901a]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{z^2 \\,{\\left(3\\,L-z\\right)}}{2\\,L^3 }"}}
%---
%[output:22e6a343]
%   data: {"dataType":"text","outputData":{"text":"Tower #1: k1 = 9.755e+05 N\/m, m1 = 82.49 kg, f1 = 17.307 Hz, peak tip sway = 22.216 um (steady state 14.741 um)\n","truncated":false}}
%---
%[output:084e8dce]
%   data: {"dataType":"text","outputData":{"text":"Uniform-load constants [C1 C2 C3 C4]:\n","truncated":false}}
%---
%[output:9fd24dea]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\left(\\begin{array}{cccc}\n0 & 0 & \\frac{L^2 \\,P_{\\textrm{dist}} }{4\\,E_t \\,I_t } & -\\frac{L\\,P_{\\textrm{dist}} }{6\\,E_t \\,I_t }\n\\end{array}\\right)"}}
%---
%[output:47fbc80d]
%   data: {"dataType":"text","outputData":{"text":"Uniform-load deflection y(z):\n","truncated":false}}
%---
%[output:4791f0f2]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{P_{\\textrm{dist}} \\,z^2 \\,{\\left(6\\,L^2 -4\\,L\\,z+z^2 \\right)}}{24\\,E_t \\,I_t }"}}
%---
%[output:7c04463d]
%   data: {"dataType":"text","outputData":{"text":"Exact base moment M_b(0):\n","truncated":false}}
%---
%[output:539c5648]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{L^2 \\,P_{\\textrm{dist}} }{2}"}}
%---
%[output:3b4ef50b]
%   data: {"dataType":"text","outputData":{"text":"Generalized force of the line load on the tip shape:\n","truncated":false}}
%---
%[output:82196bdf]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{3\\,L\\,P_{\\textrm{dist}} }{8}"}}
%---
%[output:26fee239]
%   data: {"dataType":"text","outputData":{"text":"lumped-model     f =  17.3073 Hz\n","truncated":false}}
%---
%[output:3f79f281]
%   data: {"dataType":"text","outputData":{"text":"Ritz phi_t       f =  17.3073 Hz\n","truncated":false}}
%---
%[output:977f159f]
%   data: {"dataType":"text","outputData":{"text":"Ritz phi_{line}  f =  17.5224 Hz\n","truncated":false}}
%---
%[output:060b5a6f]
%   data: {"dataType":"text","outputData":{"text":"Axial constant Ca:\n","truncated":false}}
%---
%[output:63616f02]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"-\\frac{L\\,g\\,{\\left(2\\,M_{\\textrm{tip}} +L\\,\\mu \\right)}}{2\\,A_t \\,E_t }"}}
%---
%[output:14c49967]
%   data: {"dataType":"text","outputData":{"text":"Gravity-induced axial displacement u_a(z):\n","truncated":false}}
%---
%[output:2b48a988]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"-\\frac{g\\,z\\,{\\left(2\\,M_{\\textrm{tip}} +2\\,L\\,\\mu -\\mu \\,z\\right)}}{2\\,A_t \\,E_t }"}}
%---
%[output:08550b27]
%   data: {"dataType":"text","outputData":{"text":"Bending-shortening coefficient B(z):\n","truncated":false}}
%---
%[output:10d40b67]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{3\\,z^3 \\,{\\left(20\\,L^2 -15\\,L\\,z+3\\,z^2 \\right)}}{40\\,L^6 }"}}
%---
%[output:5f8f2385]
%   data: {"dataType":"text","outputData":{"text":"Tower #3: gravity tip settlement = -1.8784 um, base axial stress = -0.0850 MPa\n          peak bending shortening = 16.66 pm, axial f_a = 410.0 Hz (12x twice the lateral f1)\n","truncated":false}}
%---
%[output:413b5a00]
%   data: {"dataType":"text","outputData":{"text":"Generalized mass Mt:\n","truncated":false}}
%---
%[output:16d986d4]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"M_{\\textrm{tip}} +\\frac{33\\,L\\,\\mu }{140}"}}
%---
%[output:2f6cc261]
%   data: {"dataType":"text","outputData":{"text":"Elastic bending stiffness Kb:\n","truncated":false}}
%---
%[output:6527b51d]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{3\\,E_t \\,I_t }{L^3 }"}}
%---
%[output:5cca9466]
%   data: {"dataType":"text","outputData":{"text":"New, reduced stiffness Kg:\n","truncated":false}}
%---
%[output:7f70666a]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{3\\,g\\,\\mu }{8}+\\frac{6\\,M_{\\textrm{tip}} \\,g}{5\\,L}"}}
%---
%[output:918f967a]
%   data: {"dataType":"text","outputData":{"text":"Without gravity the stiffness is 0.0457% too high (Kb\/Kt - 1).\nf_n: 17.3073 Hz (Kb) -> 17.3033 Hz (Kt); Ritz self-weight buckling at 2187 g\n","truncated":false}}
%---
%[output:66480527]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAArwAAAGmCAYAAACeH1cjAAAQAElEQVR4AeydCZxUxbX\/z5Ah7JsigsAEFKPE3SeCQvQZ8UnceJr4FDRBjWvQYKI40QQlMTFRcYsaleQZSVDM30SDRp9bohAxKiaKSwgJ+6qAMGyCgvLvb02fnpo7t3u6e3q73YcPZ6rq1KmqU7+6y+\/Wra7baqf9MwQMAUPAEDAEDAFDwBAwBMoYgVZi\/wwBQ8AQMARExEAwBAwBQ8AQKFcEjPCW68havwwBQ8AQMAQMAUPAEMgGgTIsY4S3DAfVumQIGAKGgCFgCBgChoAh0ICAEd4GLCxmCBgC6SNgloaAIWAIGAKGQGQQMMIbmaEyRw0BQ8AQMAQMAUOg9BAwj6KAgBHeKIyS+WgIGAKGgCFgCBgChoAhkDUCRnizhs4KGgLpI2CWhoAhYAgYAoaAIVA8BIzwFg97a9kQMAQMAUPAEKg0BKy\/hkBREDDCWxTYrVFDwBAwBAwBQ8AQMAQMgUIhYIS3UEhbO+kjYJaGgCFgCBgChoAhYAjkEAEjvDkE06oyBAwBQ8AQMARyiYDVZQgYArlBwAhvljheeeWVsueeezYr9957b5YtlFYx+jF69GjZvHmzc4yQNBiAhVPan6wReP755xsdS2ALxllXmIeCjDPj3RLf6BPlOZ7y4GJFV9mS8dHj76CDDpK33367pHBsSb+y6UjYMQomYMPxD1bZ1GtliodA2JiqjjHlGCued5m3nM\/jESzAhOs0GGXuXemWMMJbumOTpmf5N+MEuOmmm\/LfUIW2sHr1arn++usb9X6PPfaQjh07NtJFPUE\/TzzxRHnllVei3hXzv0wRsGO0\/AbWxrT8xjTbHhnhzRa5Ci8HGXvooYdk4cKFMmnSpApHI3fdv+qqqwzT3MFpNVUaAtZfQyBNBKJ8DzvggANkzpw57l4xfPjwNHtsZkZ4szwGIHmQPeStt96SIUOGuJpOO+00dxCiRy6++GK3DIDXA7wmUGHW1BWI\/eG1geb7ev+1hf8KmFdq1MMrNmxiVYS2gR15KtRBuaOPPloefPDBxCt0v021JeTJGNtHH32UpJuZO\/DAAwX7MJ\/VnjZoCzviiO+rqyzFH8pSxhd0WoR+kResU9vDZ3zBPsxP9GpLPSrUSx6YUjd61aFHUuWRH5RgO4wzPqkd9XPsLFu2zKmYSQ9r12V6f4L1+vh4Zm6sqM8X2vRtiIMXuPl2yerE\/p133hHFiDL4gz6Z0GZz\/fSxpU6Eclqn76PvGzbYIsTVHht09Iuy6INt+HnkI\/SFcozVLbfckjhPqI\/8oPj2zF77uKg\/akO9YW1SJ8cFbWKjQjnygkK9akNIOmij6XT6rLbNhWBAe8E+0D56hDj1+P2hH+jJV6Eu7HwJ2pD28\/14Ov2iXdoD1+bGkraaO0a1fa2XuhlvfNE8DamPfBV8ABPNTxZig62WI0zWBu2Sh40KZalD6wdn8hgzzltC0gh52BGSRqiPetH7km1\/qIP6qJf6EerSNvFHz0\/FNVkfqIO6qBOhHurzhTrIQ8hPNqZgRDuU1TJBHeXJV8Fn6vUFf\/BLbSiDHWm\/b34Z4mE2fvv4Rhpb8KEu6qSc3ybtYYMtZbChP+iJq1AOO1+CNqT9\/HKLVxrhLfj4cWAOHTrUkUW\/cUgkBzAHMk+aRx11lMteuXKlI68kZs2aJZs2bSIqM2fOTOiffvppp9tvv\/2kf\/\/+Qh1hr4ovvPBCCTvIIVcTJkxwdfBnxIgRBDkViBt91Erpx8iRI6W5Ewp\/KavlNERHHmlId9++fR0277\/\/PiqHwezZs12c\/vEQQmLRokXy7rvvEhXtJxcD3zeXGfsDXvjH0\/Nxxx0X04go1i4R+6NjQvv4EVOF\/mdMGN9gOxAijgeOi9CCzSj1ohasF3zol1+cdNCOfO0ncQRf6C+4kVYJq5M8+sCFlTEljdCOjg\/pTIWyHB9+ndSBr\/SDeI8ePWTQoEFEZf78+S7kjz9GGgcnzhnyKUNZxjbYBn2m72CArS\/08+6773aqTp06CePmEkn+YB\/E5YorrhDObfDRYrR5+eWXJ85n9LRP\/dRBWoVyHEccT6oDK3DRNCFpbIn7kmmf\/bK5jOMbPvp1cnzRF9URD9qQpqzaaJhpv8A1k7HUdsLCoE8cs5deeqm7Bql9WF\/wgeu0P5ZqryHHLfVjqzpC2jj77LMbra\/mmEFHHjYqlKUO6lIdIcfdKaecIoSkEcbgpJNOEkLSCPXlqj\/UFzZW+Oe3iV2mQr3UEyzH8aLXjGBeJmnqCdaPz4yt1oMPwWsKZbBTm2Qh5zvXlbq6OtH7mH+\/4r5FmvLczxg37CmHLpWk4zv9wFe\/HtKU9XXlFDfCm8fR5ILzk5\/8xBEzCBIXImZ9p0+fLhy4HMB6YnAQo\/MPcv+mvmLFCvnwww\/dTRJSjNvcSCHL1EFdfhu8Gsdm2rRpjS7E6BB\/JjrZKxFIwowZMwRbygyJzWJz4jG7Tbo5mTx5spvtpt\/4hv3999\/v+kA8KOClJIU2wQohjq32Bb8gMeiU4HDB4MKBDlHs0HMBp30IKjeJ5557zuHPOFC\/75\/Wp+QYEq03KN8\/2scP2goTHRPGNNgO\/nBcUB\/Y++0zbviEPqzeqVOnuocnv15wxpZ+0T\/ihKR9O78d7Sc+4As+gRE2tK+YUwd1Uacv6ifHA8cFeYwd9REPCv2hbtogT8ujp\/577rkHtTvWaB\/BBiUXYG4sxIPjQnt6PpBPHB3nCucMOsowhrpOmr5Rv\/pO36dMmYJpE8EHbHl9yINQE4OAQu11TKh7w4YNoscAbVME3\/CROP6GjQFlGD\/ObY4nbOkH5wFxcKcP+Ee76HzBNnWfw\/vs15HLOJjgq\/aLuvWYwdd89wuMaD\/ZWHIsJjtG8VWFYxg76qJP6BkjxoK4fzxrm9hTDjuOZ+zCBILDPYA86qYNLcuxNMU7Tomj45jCDqE9ylIHdRH3RW39MfjHP\/4h2paW51rKtZOyLekPxzbXfOqh\/\/QFP2kPXUtEr2Fh54Fet9Id02R+4Cf++njpMdvSvjFZxaQVY8hECj6AOWnihKSJ6\/0Me8qha05S+Z7J+dZcO1HKN8Kbx9HigsOFhyaYUVWCxI3zkksuQS1KqDiIOZg5yDn4OSDJc0axP1wouaBqndwIIcm+3ahRo0Tb4MLGBUbLxapI\/KfsmDFjEul8RGifiw114xP9Jw4e9IF4KuGmwBMoNhBsLjozYuSbutBBYgjBCAzADOzQIXpR0ouiElSw54aHEKcNLpjgRDkVyHEQP\/zGf2y0feJBwR\/8Qs840w5xfM8UB8qpcIGlX6SZldR68R+h\/+BAPnn0ESGerJ9+n\/ANHylP\/+g\/Nw7Ko1NBz\/iS5oGLBy\/iPokjna7gM75Tr95wKcvsFf0iruMYHBfff+wYH3ScK4wpdVJG0\/6xj+\/nnXcexRLnoUvE\/1BW+xlXpQx8e9okTQF\/rMAVnU8o8Be\/0ftjAO4cP+g5njiutB\/o8J0+EPexIo2obSZ9plw+hHFEqJt+gQlxFfWVdD76xVhkMpb4kUz866w\/zkpK\/ONZ2+S8ohx1cg5zLhNPJbwdgGxSlmsf10CuhVqGuOo4NngToA9GauOH\/nGg9xvyGReEOPcU7DgflWy1pD\/+sU3\/6QvtcG9QbEi3RCDREydOdFVcfPHFbpJFr3tOmeUfMEEoHnbMtrRvnLt67dRjQq9ztImQ5lghnzT2lCOeSvAbwSbM90zON+ooF0lJeMulk8XqBxcMLhxcQHbfffdGbgwYMMCl9cbHQcwv81Fy4aQseVyo9cBVPXVCjrloqR3luNjpeh3KcMNHTzlCla5du0rQHwiRltVQZ9W0XCah9k\/L0B444Ds+q94PwYCbneq0P1zIuaCrnlDrAyPq0z5yI6UdyNeaNWuEGT\/slWgQ9\/tKG+iCwoUZkoyeiw6hXvjBFkEXJviDX+RlgwPlwoQZQfpFHg8EOk7cdLnoo1cciKfTT+zChBsSN1lwCOb37t1b2rdvH1RnnVafg\/VyPOg5wThy4add7GiMcmDNMcX5gBBHRx42jCF90DT5vIJU7HiFhx3jRTniKrRDe5puLkxmHzwGgvXQLn5x3HJc+\/laNuhf0NbHSstn02ctm+uQccTH5urNV7+SjU1z\/oTl65iE5aFT3Ln+cp3QY02vNZzDnMvYBsUnJxwTeqyGvaLnmsi1kfpph\/aC9fnpsOs++c2NTUv6Q\/0qQdyCabVLN2TShuMFe70espaWhwR0LZXmcPHrD\/YlmPZt\/bjacUz49yt9GOC6h558+soDiV8+WTxd36nTv+ZwjlI2Wb1R1xvhLaERVFLGbM7vfvc7txSCG7aSQJ7yHnnkEecxByUHp0uU0R+IFrMW\/kwfF3Iu6KyPhPTQXci+EhywAjNO3nPPPVfQU4aLIDNn6PWkhsTrjQc9r6ogizxYUK8v\/njwNA\/+5Kf7lI1tsSSTfhbLx2za5ZgHf8oyHno+sDYSQY+OPOJ6QyHeQrHihkBBEGDmlhk4rnnaINcyiC0Psei4DrIOnOscaQhS8LqJvpyFhwNmcnl1r\/3UhwQeBHggUH2phkxWcO9hHBlj7lekzznnHCEkjZ78ZA8spdq3UvTLCG8eRwWSBaniJGQWx29Kn5r9g9g\/+H\/zm984c0iX6iFmvGImAz0hbVAHcUgiF72g8JqH\/EKK9k\/bpP\/gAB74rPpkIT5rP\/SCRv8RyvjEB6y4IEB0999\/f\/cjIWz4kQptoocgo1O\/uJkwY8tFE32YYINQ909\/+lP34zf8b+4pm\/7pmGh7Wn+mOGg5QmYbmakirjc4xUhDbpbka7v431w\/sUe0DHFmSZgtgTiTzqcoKWUWw5\/54qbODAdt+w944M84cCxwPhBHhxBHRx5xdJTXNtDxkKN4aciNM9WxQB35Eo4X\/OJY5fjw29Ex4XjCTvOCtj5WalPKfVYfg2E59Etxh7BwHOoxpmGyNyc+FlzfdNtHynG+k886Z4gc5wnnCzquj3rek8615KI\/+KTHMnEkmEaXjegkCThxbnMucc2GKGZTXzZlgn0JppPVydsnJrXI1\/sVae5jhJwP6MknjT3xXAn1+9ecsOtIeFvR1BrhzeO4QbIgWzTBj0e4UBGHTOiPdPyDmIOZNDYIJy43OZ\/ooOdCCgkm7pfRiyF6iAozApAW2kOXSnyCyYUD4UKSqkyqPC42+IAN\/ab\/xMEDXIgHBTuezPHbf4UHBmCBvX8hgcyonjwlRXqBRocwI8gNhLgKT87M2pLGVy6QxH2hDGXRQaK4OKTyHzvEHxPGWfGnf+ngQB1hEvRH6yVknMFNMdfyzfWTsaBP2PvHDySZ\/vLKP1gntrkUHUfGQGffqV9\/oEdcH\/CI+z6TVjLIcUIcHUK\/sCXO+cJ5Q5\/4sQ86hOMM3Py3B+gLKfiIr7TJ8cFxQpxx5fghznWB40r7gc639bEiD1HbXPZZzy3GivOGdvAXX4hnK+or5amLOokXql+0lSvxj2fFCCLBMcaxxjGX2EZgIgAAEABJREFUrC3ONWwQ4mqnuAeXtpCvS67AjHMYXS6lJf3xj218w0d8o2+KDWkV7ScPCgh6\/zwgjfh4gitp9P41wL9XkJdrybRvydr3r23YKAYaokOCduiylUzOt2zbKMVyRnjzOCoQlKuvvtrtCMANgtk2LmSsy+ImRFoX26sb\/kHNTZCTinqUeGGnNz\/iiK5l8tuAqJDHmtaWzlzpiccFiBMl1QWbNlXwgf7ST3xDz\/IM+kM8KNzQ+WEDei6GlEUULwiLznRgAzZgRBxR7PARW3QQYi7YxBHi6MCfeqnfJ1nY+KL2qmMckvmvNoTMtuOD347igP7222\/P6ktq9J\/yfr30gzT1I7SvfqMnP1k\/6UvYMaqYUB9CnZlImK3\/4Eb9+MSNj+NTf5zljzs21AOW\/sMXPjMO5CF6PnD8EEeHYIMtcfLCji3aIz\/VcUl+PgUfw8aAcWP8wF+vE34\/OKfIA0fFyvfTt6Wf2CHEscumzxAKzh\/K0yb14QO+oMtWfF+pizqpmzaCdfq29AU7hDi22fSLckiyY5S8dIXjmesu9viPb1yTuH6CHddr8sKEfiPk6fWT8tSDjnqpHwz0WKff2FAO7LDLpdAe7VInftBWuv3h2GY8KItv+Eh5+vaFL3wBdSPRew1KbLDV8wCdil8vuOIPttRPOz7OuRhTbdcPfR9ok7bxAb\/D+uaX9eP4zjUdHX5z7SZOSJo4+dgRz4Vw\/Oj10Ped8c1F\/aVahxHePI8MFwtmyzgZ\/KYgLryy4qTx9f4NRWcsyfcvBErs0CPJ2oAo5OJVF75ywtEWwmvmLVu2EE0qtI2oAScur5t84qJ5fshMM3bY+3rwe\/LJJxO7UJAHdmBEHHuwI+5f4CDEEGP0CFixAT1xFV4Jqq+sBdZZCPIpSx3EaYOLEPHmhAsKry7BzrcN64ef31w8Wb204x9PmfQTW2aw\/THGj2Cd6FoijJfe\/LQenYVJNu6MDXlqr6F\/PoTFw8aKeoLHFnbomjsutd18hYxButcJ+gEuvi8cv4yXryOOLf2jn6QR4uiy6TN+MutKHdSFcEzfeuutRFsk+FqsfqnjqY5RtUkn5Lob7AvnF+cZGCarg\/Y5jxnPoA066lU9D0Fgr2nif\/rTn9zaTx6UOJ40r6Uh7WbTH9rlOON402OGkDT7\/5LvC7bBdjiu6btvRxxbyC64klYh7eMMpsmuO1om2xAf6At9og5C0oG+kZVU\/PsVb6j0PkZImoI83HDtJ54ryeR8y1Wbxa7HCG8ORoATiosUywC4MASr9POxQcLsKMfFkPWEQRtOLHQIcWx9CWuDA9q3IU15iFgmJw+2lKEsQl85GQlJJ+uLtocNfaJvvj\/J4thhTzkV2qKPwTK0jQ32lCMfO+zRE5JGrwJ+5KmQVl\/pJ\/1VWz+E+EKAfV1zcfVP2wrzh\/ZoFxv8aK5O8oP1kkbvC\/2iThXS1E+a9mhX7YmjI08lWCdp8oJ9SFan1u2H+EAdKpTVfMaPcdQ8Quw13w\/Rk4\/4dRBHRz3U55chjo48bBDi6MhTSdZPzQ+GYfY+nvikZdTvsHY5TsEWv1SoW8v6odajdrSBLWnqoC61p3+0Rx5CHJ3mE2p9YXnk+0JZ7KgLob3\/\/u\/\/dttBkaYu7PGBPHT4hk6FNHrysVM9ZdGrtLRfydrR9sLCMB\/8PpOv5ZKNM\/nYaT8IOb+wJ685od+U8QWdXw7cwE9tiHN9oh10ak9IGr2275cFI603WT\/Jb0l\/\/Ho5dkhTZ5gE28E\/7UOwLP2hX\/RPhTR6v+5gndQXhkGYTuvBD9oAZ+xUT1\/wizxC0pqXTkhd1El533f6QBo9bft10QZtkUffyPPrCdqTxpZ2sMMeoSx6FXBJZot91MUIrzeCrBXStZC8miCOzjNxX7pBT74v\/joi397i0UaAGQSEXvivyEnnTKwiQ8AQMATKEAHelunvMghJ003W3OpOKvmYvaSNfAt9oU\/wAELStFkOfaMf5ShGeOOjylpC1gqxllCfdoijIy9ulvgEIK8t1I4w+OSk9hZGEwHGnAsZ67HoAa\/JeLVG3MQQMAQMAUOgeQSYpYTQYumvFWU9qk4kBJfoYRsFKee+RQH\/bHw0whtHjV+6sgaKLxbFVUIcHXmqY81hNq+2tbyF0UOAdVl33XVXo\/XD0euFeWwIGAKGQOER4BV52BpcrqtMHPFavfBe5abFcu5bbhAqrVqM8MbHgwM32SwtP9LiNQWmEF5+KOWvg0FvIo4Q6poj1gLlB5PC1MpFmJl7hLVSrJkqTMvWiiFgCBgC5YUA9wOupb6Uy3W1nPtWXkehiBHeFCOqm3srwWWNDr\/i\/+c\/\/ym87lZJd5uuFE1ZliFgCBgChoAhkBkCZm0IGAJpI2CENwVU7EnHRt+6byJfJCG97777Jn6RzDokSHA6P1pbvny5YG\/yiuHwimFg54EdA3YM2DFgx0B+jwF4RwqaU1FZRniTDDffLGdDb364pq+zCXkNw\/IHLcbC9QkTJiQInOqDIQfd+PHjBWJcQLH2Ro82DAwDOwbsGLBjwI6BijwG4B3wjyAnqcS0Ed6QUYfsMrvLQnvW54SYNFKxJy0L8Fnf2yjDS3DAvfrqq3LzzTcLa4VNHippHMaNG+dGz8artMfJP49szKIzVoxb4ccrWviAUamJjVm0jiHGC94B\/3A3tAr\/Y4Q3cACwHjcTshso3myyd+\/ews4PJkNKGofBgwe7sbTxKu1x8s8jG7PojBXjZuMVrfGyMYvueLmbmf2xH635xwBkl2UMkydPlrCZXd2bldAvx9pe0kOHDiUwMQQMAUPAEDAEDAFDwBAoIQRshjc+GCxjULLLllRxdaOAJ1zk+uuvF3ZsIJMvsV1xxRVy3HHHCWt80ZkYAoaAIRBxBMx9Q8AQMATKCgEjvLHhhLxOmzYtFhPhy1q63ZiG+tlA9t5lTRVfjoH4ks+X2Phhm\/9DNleR\/Yk0Aj179hR25yCMdEcqyHnGysYsOgNu4xWdsVJPbcwUCQujiED2hDeKvU3iMzst6AcT\/I2xNU4eNloccqt5hGHLH9TWwmgiwIX961\/\/uhBGsweV5zVjZWMWnXG38YrOWKmnNmaKhIVRRMAIbxRHzXw2BAyBkkLAnDEEDAFDwBAobQSM8Jb2+Jh3hoAhYAgYAoaAIWAIRAWBkvXTCG\/JDo05ZggYAoaAIWAIGAKGgCGQCwSM8OYCRavDEDAE0kfALA0BQ8AQMAQMgQIjYIS3wIBbc4aAIWAIGAKGgCFgCICASeEQMMJbOKytJUPAEDAEDAFDwBAwBAyBIiBghLcIoFuThkD6CJilIWAIGAKGgCFgCLQUASO8LUXQyhsChoAhYAgYAoZA\/hGwFgyBFiBghLcF4FlRQ8AQMAQMAUPAEDAEDIHSR8AIb+mPkXmYPgIZWy5fvlxeeeWVJvL666\/LnDlzmujDbE3XFL9iYGJjVhrjkO7Y23hFa7wY1+bGjOtpxhdhK2AIFAgBI7wFAtqaKT0EuDiPHz9eRo8e3UT4RO13vvMdOe+885rkhdmbrimGhcbExqz4Y5DJmNt45Xu8cl9\/c2PG9ZTrauld7c0jQ0DECK8dBRWLABfmV199VW6++WZ56KGHTAwDOwbsGLBjIMtjYNy4ccL1lOtqxd5UrOMljYAR3pIenvw6Z7XXI9C7d28ZMmSIiWFgx4AdA3YMZHkMDB48uP6Can8NgRJFwAhviQ6MuWUI5BuBK6+8Up5\/\/nnXzOrVq93SDUKn8P5gd++993qaplFskKY5pa3BZ6TYXoIvywE2b97crCuM0dFHH50Yu2QFGNtgnVp2zz33FBXaTlZHS\/TUq20EQ\/KyqZty2iftC\/0Mqws9OGEXlv\/222\/LySefLPH8MBPTGQKGQBkhYIS3jAbTumIIpIsAxGrTpk2y++67uyLvv\/++dOrUSdq3b+\/S+gcy+Oijj2rSwoggAIm744475Oqrr5aOHTs6ryGLvMmYMGGCLFy40Mlbb70lM2fOlFTE0BXO8g\/t0Ya2Rzh9+nS55557miXsWTaZdrEDDjhATjzxRLnpppvSLmOGhoAhEF0EjPCmO3ZmZwiUEQIffvih+IR31qxZcsghhyTIEYQJEvTcc89Jnz59yqjnldEVHlL22WcfgdTRY2YzWacO2Rw+fDgqJ5DhyZMnC8t6CkX88Om4446Tp59+2vlQzD+nnXaazJs3T8CnmH5Y24aAIZB\/BIzw5h9ja8EQKBkEmOXj9TIzb2wzREgasoMwo4uzEKZBgwYJRDhTwgt5OOqooxqRCF+nZPqWW26Rgw46yL1aJ8RG\/cMnfXWNP+hJX3TRRc4+mI8Ngv\/kIVonehU\/n\/og\/ZoXFoa1Sx3aB9pBeH3ul8cGPRLmB\/WSh\/BgsWrVKr+4MAOPf+Qj2NBmI6MkCeymTZsmI0aMSFgwjocffniCACcyYhFILzPBED\/KxlRN\/vv9wR9fgn1vUjhDBfX59YNVhlWkbd6jRw83yztlypS0y5ihIWAIRBMBI7zRHDfz2hDICoGLL77Yvcq+6qqrBOEVM6+cmXFj9m\/SpEmuXuw07hR5+PP4448LM8j4QPsjR450rZCGjK9YsUKmTp3qdPxBx7IL8vEZ3cSJEwkSBJEE+QiE+uyzz04Qb4jT7Nmz3d7K5EPKaZ8yqYR2mf2mDBhR5vTTT5dHHnkkgeX111\/v1oIqUV25cqXgI2UuueQS8f2A0PFKn7rIZ4nBb37zG9F\/kE5eteMf+cioUaOENslTu2Qh7ZJ34IEHEjhsWLYwYMAAlw77w9KWjRs3Op\/D8jkW8CNM\/BnjsLJBHf3ngcon5GrDGIElmNMWIeQdvdrkOgQXjot0sM1121afIWAIFA6BPBHewnXAWjIEDIHMEZg\/f75wo6dkcHkDukIIJI4ZNtqC\/PTt21d4xUyatcS8ZieuQj4knbTOSr722muO0C5atMgt0dB8bCBiEGlm7yAzECfIpbYJCR0yZAimKYV21a\/+\/fvLfvvtJ77viiOV4AfbMjFjio\/oaIcyzLJCiO+\/\/36BBPNqn3z81PpJQ1jpO+VII8TRQRRJpxLGFlswTGVXiDwIK8Tbn7G98MILhWUU9Nv3AWwg5v4YMVakGTvG0LfPVRyyv3PnTmEde67qtHoMAUOg9BAwwlt6Y2IeGQJ5RQBiwat8bvQ0BMHaY489BHJBulDiE8V02gySOPyvqqpyRAWy8u6777otpXxypQSRfEgNZbQtCCn91nSyMNhuMjv0tFNVVZX4MSC6jh07CrO1EFEeLpi5DvbdT2MXJIqQRnTUFyXhgYLji9laX5TsMnPLjC99UmwgxP4YkiY\/X8IxUVVVfxzlqw2r1xAwBIqPgBHe4o+BeVBiCJx89xsSdRk7bW4oqhAMyBOv5VlCALGAUEAMibNWM7RgiSshiczEQgp9YkWc1\/HNuc\/sIetkwQAhjq65cunm418mtsmIIktN0q1H7SD2SrhVFwwh6sEHAt+G4wJcwvApj7MAABAASURBVEQJq2+fbhyirz9ew4e6ujo3+8u4+TJjxoyCP5Cl2wezMwQMgfQRKKalEd5iom9tlyQCsxbUSTlIGLgQJl4n8wpdCQVxdKTTIYdh9Tang8xs2LChObOU+cyMMguoRtSpJA3iBFlCp\/l+GDaLx0w3a22xY3YbUgUGCHF05GUitINPvh\/aDj6yzIAZ4yAB9tPYBfuaiQ9htkOHDhVd\/qH5PPzwwzj8mzJlituNgyUbmu+HHBfgEiY6W+vbpxvn4Ysx4OEC7Lp27So+FunW0xI7xooxo\/2W1GNlDQFDoLQRMMJb2uNj3hUBgce\/eYhEXe46c2BS5CAUkCoMIDv+8gZ0LRWIQ1VVlUCiqIs2WLdKvHlJbrFs2bLEnqnU+ZOf\/ER05wFmRFkni448aoFEMVMLsYO8su6WH0ShJ58fxDEjTDxXAmFkVwvfD9phuQWkk9nW8847z+1Dy64UtMsMKTPsxBEeQAjZNYMQwZbdHrAlnUpohzH1Hw5YL3zttdcKs\/paBw8\/LOmAdDLj7687TlV\/LvMYF3xg2QNxxogf9NFfbYfZZSXmqstlCOGtqmq8DCWX9VtdhoAhUBoIGOEtjXEwL0oIgWEDuko5SDJIfcLLj6ywg6gR5kIgLnfddZfbgYFX4BAwdh3o0qVLi6qH1ELkqBOSBlFi5pFKIZLMUhMnDxvsIVAQO\/SEpNGTzw+k+FEbebkS\/GC\/W3xTPyBwkF5IJ+0wI8oOEpBP\/ICEf+1rXyPLCfixAwQ7B5CPYEsZyjqjFH8YS3azgET6ZpSF4NMedSJKtMGVH8b5RNMvm884P1jkmKQNxogf9NFf\/EOYAWZswRabTISHJB1v6lJR0k9dLKlgCz5wJ21iCJQMAuZIThEwwptTOK0yQ6D0EYAkQn7wFBJ23333JT44gc4XSAYEDiLi64Nx6kRUT71z5sxx23YRnnXWWe6LXughFjNmzBD1gTLE0ZFHOlm7kD59re6355fRfMKg36TRI\/SLvgfroS4V7LHDH3SEpNGTRoK+o6NO2kDoP\/1Gr0IZ8hD6\/YMf\/ECol\/qxAQf05KtQxs\/TNDpfqIM1uxA5X088rF6tH1IeJMmUyVbAyO9Tsnp4I\/DGG2+47dOwoZz6ROjXQZ6mtS\/JcEBP+TAhj7aY7Wf\/4TFjxpA0MQQMgTJGwAhvGQ+uda3sEbAOGgKhCLAsAiKXyYwtpJwHk9AK86hkRprq9W0D8UIJM9z+F+kK1a61YwgYAoVHwAhv4TG3Fg0BQ8AQyCsCzH6OGzdO\/LXEeW2wBZUzI82HPdinuAXVZFyUh4Enn3zSfYAl48JWoAQRMJcMgdQIGOFNjY\/lGgKGQAkg4L\/KLgF3IuECr+319X+pO8z4IoX0kxntJ554wrY7KyTo1pYhUEQEjPAWEXxrurAIWGuGgCFgCBgChoAhUJkIGOGtzHG3XhsChoAhYAhULgLWc0Og4hAwwluEIZ82+z1J9iWsIrhjTRoChoAhYAgYAoaAIVDWCBjhLcLw3vjMIoH0Ll23rQitp9mkmZU9Amzor\/uRsj0Tm\/sT0nF+0MOHDnTfUvL0gw7kB4WPO6htMKQd3540orp33nlHtF3VZRvSHz42kUl92J588slCn\/128dHvS3MY+GUzjdO2jzf90DrAnbZ9XzROGcqqbSYh\/UMow\/jRBm2R9gUdedj4+rA4Ngh5+IV\/6qvqyUuGOXkmhoAhYAjkAwEjvPlANc06jfCmCZSZ5RwBSAwfG2APVCrna1N8rIBP30JU+AiBv+ctH1K48MILE3ulUiYobPDPPq7+vqd86IAPKCixogx71CLEIXZjx44lmhPhAwZ8upd+pFshXzTjwxj8iIky9B+iRtzvC3vb8jEJfCYvV0J7Pt7Tp0+XH\/7wh00I+FVXXeX2NfZ94iMNl156ac4eGFrSJ44p9tPlQyMQWvzCP\/ylT3yAQ7GL0i4SYGJiCBgC0UfACG8RxtCIbhFAtyYbIcBnZ33Cy5ZQbA3FFlHE+UwvBFYLsTH\/8uXLJdO9UiE2EyZMEEgvJEjry0cI4eLraRBT+pFOG5BN9qtl31rsqYOtvCBqSsrRI+wiAOnkS2W57At488U3dlWgHYg3n0xGTzqVqN88aKSyK0SeHhvsq6v+qH\/0CUzvv\/\/+xEOTHl88FBXCP2vDEDAEKhsBI7w5Gf\/sKlm23pY0ZIeclcoWAV4r84oZsgHRICTNLCfCTCzELpfbWTHr6vtLGwizfcwaL1u2TPCDtPqHT75g79cRjFOO2Vf6RD+oK2gTlp4yZYrw4QGIOfmQNh4ElKih84WZWGaQldD5ecRp1\/fbjyfrA3gHyTV1BXFDl4lA3lmKoD5kutQjk7aw1bcE6T5sYHfeeeeJT4Kpx8QQMAQMgXwgYIQ3H6hanYZAiSIAueIVMzOVCHHIGzOMvHYOI150hdnGqqoq0SUQ6NIRZkKnTZsmo0aNarLfKTOakydPlr59+wpElbT6h18I+Sy1YIY5WXsQSdrAtk+fPsIM7BVXXCGQz2Rl0OMbM88jRowg6YR+0l6yJRGQNJZ3hH22lwroA36HSTJsKecLfj\/33HOSqs9qD7knDtknVKFvLNNgtlt9YQxOP\/30vC1\/ABOHZcwJ9YcvmcWSrk3GCH\/AEB3C8ZTNmwPKmhgChoAhkAkCRngzQSsHtjarmwMQrYoWI8Ds4YABA1w9weUNTun9gYBBrFiaoDOhXnYiCmmF6OiMIiEztxAtiGzCMM0Iyw0grqwl5pV4WDFsXnvtNbnrrrtcdlVVlUDeEfrolEn+MCO5c+fOjEl8kuparKYvrB1m1hv\/g31mDMDUFwj7I4880uRhgocYZqKZkVbHiKNTEqr6XIQQ7JUrVwrjT30cJ\/gFycVfjgOOn+BxAOGtqqoSxoJyJoaAIWAI5AuBYhDefPUlcvUuXbc1cj6bw9FHgFfdvLaHbNAbyBGzlpAU0r5AdiFgzAYze+nnBeOQGupiRhHyy8wtSwOCJCdYLiwNgdIfPaVqlxlZZnVZN8oM46BBgxLkD8ILiWQHBuoLa6eUdBDcOXPmuB+m4RfLERgr4ghjALa+zJgxw\/UXO2a6tZ\/0nTGAgEI4EeLoqCvXAmH1Z8Y5biDtPIjgL8cFSxeCfWImHRKOv7n2yeozBAwBQ8BHwAivj4bFDYEYAgtPr5JMJVYs6f9M61L7pBXGM9QuLFz6zf5xq8aBrnXllfnIkSMFIgShZdaPOKRJS0BayINoZUpaIc+QHdrx69S6U4WQt8svv1wgr821C1GCrFMfM4zMWlNe4xCxjRs3kp2WsMMADwPMeocV8OsOywczcAyTcBzCahG3nCGTV\/26TABiSY3g4j+AQDpVmsOU8pkKDx76o0cwgtzyIzVIPHXh39VXXy2Z9IlyJoaAIWAI5AoBI7y5QjIH9bA378E\/+qu8NL8uB7VZFYZAUwQgO6x1ZeZVCRBxdKR1nSnEDbKLnjJNa2peA9mB9ECmqa\/5EvUWEydOdBENXaKZPxDUFStWCISXWcx3331XBg4c6Nbz8qM4lgM0U4XLZqaYmUp8dorYH2aJWXtKqHVDjGNZTf4zGw2OYaLYBgtBhJGgPtM062eZ5aYcOIAHuJDOp0BwdTuyTNvBP\/zE30zLmr0hYAiUGQJ57o4R3jwDHKx+aYqPTbw0f72QP2vB+mAxSxcQgT0f2SmZSir3Mq1L7VPVSZ7ahYU1P1+ESagw+6cEA7LCjKYub6AAxI61s5BdCBy6bIV1o8w0MuNHW8F6aJd1tMzEkscMNLPCzAYyK4gulUDysEeoB1t+tAbRZlcC1o2ytIJZavJ8gdyyHELbJo82b7\/9dmHtqZJQiPu1114rzIjzEEDd6LDPhdAHCLY+FIATW6OxNVkm7bBkgZltljXwEINvPtFnXFkjrO2QnwuBtPrHEBiy+wL77tImbWifwBvc0SFgz7hxHJA2MQQMAUMgXwgY4c0XslavIVCiCPiEl224cNMnIWzVBYGB3AVfzWdKlpT8MDMaNmNLu5AgyCQzoBBN2ibttx1f+4mrjQRCDgGF2PK6HJ+JpzMrjW\/M3OqsqFbMcgzWxpJWH6iXNAKJhJgTz4XQBx4uaIP2IK4s0wCPTOrH77333tv9AIw4PxrjR23UiYApPwCkvUzqVVv6TT2+MC60gb+0qba0QVu0ib32KbjdHUshGH+OAy1roSFgCBgC+UDACG8+UM2yTmZ3KaohcRNDINcIQKQgJNTLDOJ9990nkD\/SCPlhr+TRaTnsfIFgBsmM5lOGstSLjhAhTruUI58ZVYgm8aBggy1lgkLbzGgilKO9oE2yNEsT+PAEs6JBG3ykvqCwTnbLli2JDygEy2WTxme\/HdrWeug3\/aefqksW9urVSyCR5ENAg3jSDnkIbSDEqZs2aIu0L+jI8\/3TOHq2P9N6\/HK0pXaEQRtmfflQCLPBtOGXtbghYAg0h4DlZ4qAEd5METN7Q8AQKCkEIE68ytdlGpk4B+HnwxMsKUi3HOSMJR+E6ZYplB0EnvW0YFKoNrNth1l\/yrLkhdDEEDAEDIF8ImCEN5\/oNlO3zeQ2A5BltwiBSinMGtJkP3xibWhdXZ3wOj4ZHqzvffLJJ0XXmyazi4JelwboUpVS9ZkZ9TvuuEPSXatdqv0wvwwBQyA6CBjhLaGx0o9SGBEuoUExV0oeAX11zyv0oLPM4LK3bfB1um9H+SeeeEKw9fVRjDPrzBKVUu9LOWEexeOkAn22LhsCYoQ3AgcBBPjGZxbZdmURGCtz0RAwBAwBQ8AQMARKDwEjvKU3Jk08mrWgTm58ZrFMm72qSZ4pcoSAVWMIGAKGgCFgCBgCZYuAEd4IDS0zvRFy11w1BAwBQ8AQiCAC5rIhUI4IGOEt8KimIq3J8pau2+q81DW+LmF\/DAFDwBAwBAwBQ8AQMATSQsAIb1owFdYoGbFNRogL6x2tmZQDAux7qx+S4FfzfESAkL6xYwFf5eKjAQh5qba64kMM2IUJ7VCnCmlE0++8845ou6rLNqQ\/fGEtV\/Vl64eVMwQMAUPAECgtBIzwltZ4mDeGQEEQgLzyRTO27aJBPvHaqVMnad++vduei08C86UsPhiA8CUtvgRGOezDhP1U+SgD9irstcqXuHyCy44JCHVAUMeOHUs0J8JX5Hr37u36kZMKrRJDoDkELN8QMAQigYAR3iIOU7oztqnsps1+z3ZvKOIYRrVp9q71CS9f5zrkkEPcF9eI77fffgKB1f6NGTNG+HRvpvu7sv0Un\/qF9OZ71hUyzpe7+Fww23Op7xYaAoaAIWAIGAJGePN\/DOStBYjw2Glz5ZSfv5G3Nqzi8kJAlx5AZpl9JWQZAh9mQJiJvfjii+Whhx5y5DcXvWfW1a+HNhBmd5k1XrZsmSPXpNU\/fPIFe7+OYJxyBx54oNAn+kFdQRtLGwKGgCFgCFQuAkZ4y2TsIb9l0hXrRh4RgMyy3ICviyHEWYZw3HHHyfTp00WXGgRdYNa7FAl8AAAQAElEQVS3qqpKdAlEMD9ZmlndadOmyahRo4TZXt+OD0VMnjxZ+vbt64gqafUPvxDyWWrBDLNf1o9DhmkD2z59+sj1118vfPrXSK+PUqnEzQ9DwBAwBIqDgBFeD\/fgD3X40Q46z8RFuZH6s0\/MLrmMAvwxYlsAkCugCWZdBwwY4HoaXN7glN4fjndmTVmaECStnpkjrcyy+ucGM8iQXYisb5tOnHMP4spa4mRfDsPmtddek7vuustVWVVVJZB3hD46pf0xBAwBQ8AQqHgESo7wFmtEuKmPHDlSLrnkEmFmCSGOjjz1izg3YWbDsCG85557pKWk1yeyfpx2g2l0iK\/34+SZZI\/A9rWvSDritxC09\/OIB\/OTpbFVCdqoXsNgvp\/+9MPlatYkZK2rv36XGV5+lBZGZjneWXbAbDAzsE0q8xSQW+rivGBpATO3p512mmRDdpkZvvTSS935mKpdZp6Z1e3fv788\/fTTMmjQoMRMMoSX85IZYPrR3E4TXlcsaggYAoaAIVBmCBjhjQ8oN0tu2Pw6Pa4S4ujIQ8dNmNelzB7pjBMhxJhXquRjl0oKRUxpB0nli+WFI7Dp5VGSjvilg\/Z+HvFgfrI0tipBG9VrGMz30xtjfVA7P4QAMgv73HPPCQ9zzMZCaB999FEhDjlUe0gieZDdiy++WNVphZBnZl1px68zncIQ8ssvv9yR1+bahdRC1ql35cqVwqw15TWOHoE053JdMnUWQKwJQ8AQMAQMgRwhYIQ3DiRrF5PdELl5chNl66a6ujoZMWJEvFR9MHToUEHP7Fa9pvh\/p81e5X7MZqQ387Go3nWIpCN+zUF7P494MD9ZGluVoI3qNQzm++lW7fqoWaMQAslaV2ZemYlFiKMjznlAASW76CmDLlPRh0HINPWlW37ixInOVEOXaOYPyzJWrFjhCC+zy++++65wXmox2meGl3OYEBIOwUfIw47lESxjQmd7+YKIiSFgCBgCpYJAy\/0wwpsCQ72JMoPENkfcLDFP9sMdZpvIL5T4H6jw47T\/0vw6gezOWlBH0iQDBDoPnSbpiF9l0N7PIx7MT5bGViVoo3oNg\/nBtNoFQ45TZkLR8yDnL29AB\/Fj2Q5kl5lRdNmKviW5\/\/77hbaC9XAu7dy5U\/TcYgaaWeGrr746rV0iePjEHqEe6uctDG9dINykkwkEn9lrfFu8eLFcc801MnXqVLekiXXHrFtOVtb0hoAhYAgYAtFCwAhvivHihsfMbapfiFOcm3bXrl2JNisrVixP2BDnRrtt2zbZsWNHQk8EnQppFdUR7tixXdWx8tsFnYoSYGxUZ+G2RhgFMU+AWeYRn\/DqvrqsgdVuT5kyRSDBLGdgttMXnQ1V2+ZCHhTPO+8894O2sBlb2mUNLssrmF1maRBtk\/bbZVY2jDBDyCG3\/KCOfYLxmXhzs9JK+DVkOdLSpUsTyzw49wuxd3Bz+Fm+IRA1BOw+0\/g+U0g84BO84VJZEXvrFbXjJ5\/+GuFNgi4zTbyK5Wba3ExRkipC1dyUNWPpum0CuVi1apUgqickrcIJgw55a2G9LXnr1zfM3hJHp0Ld2L+7ZK2rW\/UWNuC3du1aIKo4gVhCFOk4x\/Z9993XaDaVfGY\/w0TLUdYXCGayJUGUoS7qpQwhQhxCTDnyWWYwY8YMN8NK2hdssKVMUGibZRkIZWgvaJNOuqamxhFz6kDwhbXI6ZQ1G0PAEKhHgEkiu8803GcKiQV8gsmB0aNHC+H48ePrB8X+OgSM8DoYGv+B7DLDw+tObqaNc5umeB3LSd40p6lm8OAhCSUzW8we9+rVS3bbbbeEngg6lerqalROunffTVTfrVvDrDJx1RM649if9durE\/boTXol8Eh3Vj4Go\/0vYQSY+WWdvc7WZuMqxHbjxo3CQy7luQZww6Bu0iaGgCGQHgLdu3dPXGPtftNwvykEFvAJJgdUxo0bl96gVYiVEd7AQDPLlIzssnQBcwguYVCyueH269dP2rZt68SvT3WEPuEljg5ZubFhGUR1dWtXB\/rVHzbUVB0jy+hUyNvjmr\/K68u3Jew1rxLDBqQqLFZG3dW19tmcfwoDN2l2lWCLQZZSsLTi9ttvbzTzrbYWGgKGQHIEgvecSryvFKvP8IkhQ4a4L1cSDh48OPlAVWCOEV5v0CG7zPDwY52wmV0Ib9euXYU1kF4xYS9Q9Gz35OtbEtclCdSh63GJJ5Ol67Ymy2qkZ\/cGFBoSNzEEoowAs7MsPwhbysB5zPIJ8pj14BwmRE+fVc9yCZZ3zJkzxy2poD7qxcbEEDAEyhsB611lIGCENz7OvMJUsstNMK5uFHAD1F9v6493+EU7s0LoyW9UIIOET3AzKJaxaaHaydgxK2AIGAKGgCFgCBgChkCeEDDCGwOWX2jzCjMWFX7lzStNX\/w9OZkZYgZY7fg1OT9sQ0\/5fEs6s72pfFDCq2EqW8tTBCw0BAwBQ8AQMAQMgSgjYIQ3NnrMzPIKk19mhwl52MRM3X9mgH27QpFd13iO\/rSUOOfIDavGEDAEDAFDIEoImK+GQEQRMMJb4IErBNH0Z2\/9eKquYnfjM4uED1aksrM8Q8AQMAQMAUPAEDAEooaAEd6ojZjnLyRVk35cdWFhMsJN+RufWSyXPjw3rFgmOrMtcQRYf+4v2fHjwa24+CEnol1i+c8777yjyRaF+OEvF2pRZV5h1tUfddRRQoiarcXYiJ04Qn8Q4ukIfcZP\/PXtqRe8wI\/fAPh5LY1r3WH14jttEmbajtYb7Eum9eTKnv6BIX5lWyc4INmWz1U5fEDSqS\/ZMeWX5TzDztcVM07fGK9i+mBtGwItQcAIb0vQK1BZyGium0pWZzJ9rtu3+oqLQN++fRt9ZIElOm+99ZZzivXpSkDY4QAhg5vv6aefLu+99x7JFgu7nfTu3Vvat2\/f4rr8CthtYebMmUJIP+jPm2++6Zu0OK718iUjyHShljVBOviM8vTp00XHJZPO8OlkPpXOsqxMyplttghkV44HkrFjx2ZXOE+l2Jf+ySefTDxI5qkZq9YQyBsCRnjzBm20KvZnfo30RmvscuUtW3NdffXVwtcAFy1alKtqQ+uBMEJKmYml3VCjElXiOyQasvvII4+Iv74\/ny4r2YW0QuYzbYsHFggLm9NnWtbsDQGO8xNPPNF9HdTQMASiiIAR3iKPWr6b94lsqrbS3cc3VR2WF30E+KjKhg0bEh2BZCGQJWZ3ly1b5nYyQcfrTV6tB4W8RAUhEcqxZzUzo3zkhdmsELOECoLJa2\/KqZIyBx10UKPZJtpFWMoAkX7ttdecr9oOdVAXdWzatMl9elN9pxz65oTyUSS79IttFzt16iT9+\/cn6YT+gIviEMQULNFpvo8TeeBM6CqL\/SHu67C\/6KKL0sY66A9LSTj2YlUn\/lOn+oPvjGUiMyTCcYMdfmg56qBe6lcdx5RfHBvNAwP6liyf+oN+pNMXvz6N4wfHGOcaHw8gTV6y+lRPP7FDKBP0mf4g5CPEk\/WPvOuuu04UH+qjzGmnnSbz5s1rdN6hNzEEooCAEd4ijlIxZlKzaTObMkWE1ZrOEgFunPfff78cd9xxbjmAXw2zO8xmshSCbfl4nc5rfJZCqKCHUKWaQeRGyhaA2PJp7euvv16uuOIK0Ruq36bGmQGGRDEjjI\/oWQ4BwYCgk4a8zJ49W0aMGEHSSbt27YR2hgwZIryO5YMT1EUmywLOO+8895EJCDFlfcKATVBoGyIStZld+oHv4AeOigE6+kM+y1kYR7ZYvPTSSwU8GRO2Xbzlllsa4cQYUiZdAWvapf5UWK9evVqYQVRb7NnfnAct\/KE9xoixoh7ysaV+8lIJ9occcojrB8tBKEO9HNPUw\/HBsUg74AKBXblypfi4nH322Qmi15wf1NNcX5L5y3ITjlvONfwmrfWxHAV\/kUGDBgl94EuD4MD44jv1pjo\/sKF\/2FEPwhj7\/SPvD3\/4g\/D1QfLxAR3LjzjH+dgSaRNDIEoIGOGN0mil6Wu6s7p+dUvXbfOTjeKn\/PwNOfhHf22ks0S0EdDZI53hIWTWlRsqZDbT3jH7BXHlxpnsdTs2zLpyE6X+qqoqR64h2Nyg0SWToUOHJpZacMN+4403hDWOTz\/9tCsCMSFCHwibE0gwgh1kHvLQnA\/f+c53nA91dXVy+eWXC35Inv8xAw4569KliygRzaZJSBFE3f\/8MstWWL7CMhYlwTzE6DaMYMuMnpIdcJowYYJAOCFg6foBzpAp7LUOHnqCdTCGrOlWW+yJo2N2GnvK4QP1aD71E08lkEf6gk3\/2Az3fvvtJ5Bprac5XPCDMhC9dPxori\/4kYlQH\/YQc0JE4+Rlcn4w7jwsannqYow5D6dMmULSCf0FK5eI\/+E44RrR3LkSN7fAECgpBKJFeEsKuvw6k4qAasvp2KhtWJisvK8n7ktYPeWmW\/LyZlHRvmlaQ9UTqo6QNELcF3QqYXpfR1xtCUmrkEY0rSE6FXQbln2sydAQAsDsEbM33DAhDcjEiRND7VMpIQCQMWYHuXEms4UsMKvLTRQyBclUwsENFEJ88sknu9nFYB2UoSwzupA3bti8suX1Ku1T3q8vWD6Y5qbNzTuoT5UGn6eeekpYQ\/vuu++60Ldn1i\/T2U+\/fFiccYLw\/u53v3PZEGAXycEfsKyqqpLdd9+9SW2QeWY4fSKIEQ8UnTt3FsqSTkeCWNPezp07m9TBGHJM0gYPYAhxdLRDm5SjPGmEMaR+4qkE0szsZCobzaOdqqrGuNAOs6j4SH5zfmCH3\/hPPxDi6LSdTELqC\/aBc4djnrxMzg\/85\/jleMYvFR4qfJ\/AlX77OuIcExwbHCOkTQyBqCBghLdERirdWdlirrWdtaCuRNDKrxsPnrpQVLQlTWuoekLVEZJGiPuCTiVM7+uIqy0haRXSiKY1RKeCbmqsD5puLuSmdvvttwszgLzizuRGhi2zndx4mR1M1RY3Zm6i2HDD5MZJeY1zI964cSPZTQQfKQtRhqAT54ZNCGnnda6\/nKFJBTlQsAQCP5jBhtxDPnntn2nVkGIlGYSp6tBZSMgNM5uQEog1bRZDGCN\/jXcufeD4YEwZX8bUl+aOLfzgwUfXnIIrcXTk5ULwL916sG1JX9JpR88dbDkuORfSOT\/wjQcpyLePMfFs3u7QvokhEAUEjPBGYZRy5CMztZlWlU2ZTNsw++IjoISKmyAzmOl6pDPCGqZTjhlayDWEl\/aYbRo4cKCwhpKlFhDJsHogtLxOZ50xcW7yrMv81a9+5ZYa+DN\/YeVzqeMVN4QGn31SBXmHaEG4ILZhbUIqIBcqqWbF\/fLY8Vr+nnvuSawlJZ92aI\/ZcSQVgcZeBbyYqYTEqk5DsIVAQY5UR4htVVXj2U\/0KuQHCTGYQM58m6qqpnVwPHBccHyorR\/ib1VVVaOZYeqlfuw4hmfMmOHW6YItcXTkZSK0E8RF28FH8quqkvtBW9il6gs2cApiXwAAEABJREFUmUhYfeBEG+RRF+dEOucH9izLYawol6lwTHBscIxkWtbsDYE8IZBWtUZ404Kp9I2SEdNk+mCPfLtl65Ov5w2WK8d0zZEdREX7p2kNVU+oOkLSCHFf0KmE6X0dcbUlJK1CGtG0huhU0HXp21qTaYdKqFgn6ZM4rYBXwrxW5YaHjplGXrf7a0DRJxNuyNgjEArsIIzMlkISmcFk5slfW4iNCq+EiUOQIR3EuXlTH8sdeK2Lzhduytyc1Wc\/ryVx6qXfEAdmuCFE1IdvrFGGyEM+0iWflE1HwKZr166J9bzUTzu0d8MNN8jSpUtDq4H4MQvv4wBe4Oav22RZCa\/uCRkvZpRpg0o5Jhgv6qE+xqCqqiqxTRUY8DCCrS\/4pg9RWofOXPt2kHnS\/gMPfrB0BR9ok3L4QD3YUi\/1E8+VKC4\/+clPEuu0aYexZa1sOn4015fmfAVbzhElpXrs+9hoXPM0xE\/K00bY+cGDGutz\/f6BJ+cg5zTlkgljzAMG9SazMb0hUKoIGOEt1ZEpQb\/CllP4RLkEXc7KpbMf20tUtAJNa6h6QtURkkaI+4JOJUzv64irLSFpFdKIpjVEpxKm07zmQnZYCJI4LQPJgwxxo\/3qV78q7HrAWlp+yc8Mowq\/AOfGqOU0hFBDbiG2\/FiK5RPE03ldTR0QDcgWN2tICTpu8pBk\/MI\/dEFxxO3RR90WS9zYg\/nZplnawI\/0IFw6w62+qa8+wcy2Hb8c9UKodYyeeOIJARP0YEL7vr0fBweWfujYgBe7AUBgdOyYub777rvdLh2MFzsa8GNE8iFKtMUMNfXSJr7wwEE+ZJCdCfiBHfkq\/BiKdrGhDkhr2JhTH7smQOCxRTi2wBhfqI9ylKce8qmX+snLlYALxzYPShxftMOsOqSXMaed5vxIpy\/Uk0wYSx5G6D8klPrYQ9kfK3ACL\/Koh5Dx4RigPDr8D54f9I9x13z6B57gSr\/QJxN+8Ma5y1gnszG9IVCqCBjhLdWRifuVDaEMI6bx6hoF6dSdarb3xmcWud0bXppfGWt7G4EX4QTkIdnrXm7oc+bMcWSWGyPkBtHuckPkdTE\/ooJsEA8KZIGyWsYPKc\/sF0I5fPHzm4vji18\/N3n6Qr1alj7gGyE62qAt7LCnDoQ8FdKIpv2QMpSlHl9PHB11JyuLTTYCfvTT75fWQ790jFq3Tn8mH\/LDAwqkRevSdugDQr3Ur\/nE0ZGHBPvp52N31llniY899bCNFX2hPOL3iTh5+IGtYo2dChiTp0IZzaPsfffdl\/Krc9hjp20QkkavddIGY0z7qqOv2g59o6+aR0h5zae+oB\/URZ1qQ0g7lNU8TaPzRX2kDO2Q5+vQUzf1kKeCz\/iCLTrysdM60CHkY0c9Kr4N9SDY+sIPTw8\/\/HD3QOTrLR4tBCrVWyO8BR75dEhmgV3KujkluqlIcdaVW8GyRIDZRWapyvWVKK+TIZTMJDMDl+9+MmtLO7RHu7Sf7MCB\/DAD6y9hSGZrekMgiADHGLPMvAUK5lnaEIgCAkZ4ozBKBfYx2QxxkKwr0X1p\/voCe2jN5ReB\/NUe\/KGN3xLrDnlVz5IJXx+leE1NjVxzzTXCK2JeLyebwctVn6ifdmiPdmk\/Vd0sWeCBgzWxqewszxAIIsB5yQNTcKY7aGdpQ6BUETDCW6ojU2C\/lLwGm01GfrELEmB0JoZAKgSYZeQVK0QtaMeNlFfHYa9Sg7almOaVMGtqEV4TF6oftEN7\/GiNHRJ4cEiGj77KDsM\/WZmW6PENaUkdVrY0EGAcOcZLw5sCemFNlQ0CRnjLZiib70i+CGq+6m2+R2ZhCFQ2AszU8qMjhB848YNBHhwqGxXrvSFgCBgCTREwwtsUk8hoks3KFqIDRnITKFvEECgaAszUMrurYjNwRRsKa9gQMARKHAEjvEUeoGTEMRmZzdSe7mVThnJBSVZP0I70tNnvEZgYAoaAIWAIVAwC1lFDoHQRMMJbImOTaq1sMV1Mh+QGyTm7N4ydNleQYvpubRsChoAhYAgYAoaAIQACRnhBocwkHZIa7LJfxo8H7dJJz1pQv2tDWD3plC+0zauvvip8PMDkFcPhFcPAzgM7BrI5BlasWFHoS7e1ZwhkhIAR3ozgMmNFwCezfpx8Tc9aUNofpOBLRoMHD5Y77rhD+LGPyWjDYbRhYOdBwY6Bsjrfxo8fL1xPua5yHzAxBEoNASO8pTYiJexPcOlCMleV8CbLLxU9F+abb77ZfVWMrw75cv\/998utt94qv\/71r0PzfVuLP1QSGNmYlcY4pHs+2HhFa7wY1+bGjOsp19VSucabH4aAj4ARXh+NUovH\/Fm6blvsr\/3PFwJcnNm0P0wOOuggOeyww9xHBMLyTTek5LCxMSu9MUl1nth4RWu8GMtUY8b1NF\/XaqvXEGgpAkZ4W4pghMqnO0Prd6mlhLul5X1fLG4IGAKVi4D13BAwBAyBliBghLcl6BW5bEvJZEt2hkhFnlPlFRkya94QMAQMAUPAEDAEKhCBMiK8FTh6XpdTkd9UeVQRzA+msWlO0imDzY3PLBLbo7c5NC3fEDAEDAFDwBAwBHKJgBHeXKKZRl2ftt81DavSNIGwpuOZbxeM3\/jMYtufNx0QzcYQaAkCVtYQMAQMAUOgEQJGeBvBUR4Jn2SWao+i4GOpYmd+GQKGgCFgCBgChkB6CKiVEV5FIiJhqRDFbNb\/6gcpgLpU+oEvJoaAIWAIGAKGgCFQ3ggY4S3v8W3Uu2xIaro\/QFMCq6E2nG55tbfQECg8AtaiIWAIGAKGQLkjYIS3REZYiWLNLm2dR9mQU1cwIn+MCEdkoMxNQ8AQMAQMgcpBoIx7GhnCu2bNGjnzzDPlqKOOylooX8ZjmbOu5YuMKqnPmaNWkSFgCBgChoAhYAgYAmkgEBnCu3PnTlm1apVs25bdl8coR\/k0MCmIScNMbnb9aYmTwbaViKo+Vd1qm8omnbxyn8FOB4OI2pjbhoAhYAgYAoZA5BCIDOFVZG+44QaZOXNmxkI5raNSQiWw+ZixTZf4+sQ2VRny2KO3UsbG+mkIGAKGgCEQdQTM\/yghEBnC27ZtWznhhBOkZ8+eWeFLOcpnVbiIhZSsKnktoispm1Y\/Ia4pDZNknvLzN4Q9erMtn6RaUxsChoAhYAgYAoaAISCRIbydO3eW2tpa2X\/\/\/bMaNspRPqvCJVhIiWG+iXC+2tF6gZY4QnzWgjqCshHriCFgCBgChoAhYAgUH4HIEN4gVJ988onMnj1bfvCDH8g111wTKjfeeKNs3LgxWLTs0kp6lTRGuYP+Eogo98N8NwQMAUPAEGiEgCUMgaIiEEnCyw\/Y\/vd\/\/9ft2jBlyhR5+OGHQ+Wpp57K+kduRR2VFjauywv6dqvf4qyF1bWouE\/C1a9ghb6NHw\/aWdoQMAQMAUPAEDAEDIFsEIgk4V2\/fr384Q9\/kCOPPFL+9re\/ycKFC0NlxowZ0qNHj2xwKcsySiaTEWHVq10qENRGZ5fVNhmp1fyUoWUaAoaAIWAIGAKGgCGQBwQiSXh37NghW7ZskZEjR0q3bt3yAItVGUSgb3y2WIluMD\/TtF+PkeRM0TN7Q8AQKHcErH+GgCGQWwQiSXi7dOkiAwcOlAULFgjLG3ILSWnWpoQzF+RQZ2WVdGqdqs8EAfUrnbW32l6q+tOxSVXe8gwBQ8AQMAQMAUPAEAgiEEnC26ZNG7n00kvl+eefl3vvvVdWrlwpq1evbiJ8nY0ftwU7XR7pwvUiGyKciXfpkOVM6jNbQ8AQMAQMAUPAEDAEfAQiSXjpQIcOHaRjx45y8803y7Bhw2TIkCFN5Ktf\/ap88MEHmJtkgUA+Z1t1VjnoVlCPD9Nmvxc0s7QhYAgYAg0IWMwQMAQMgWYQiCTh\/eijj+THP\/6xzJ8\/X0499VS59dZb5c4772wiEydOlE6dOjUDQWGzP23f3TWoSwFcwvujeoiep84qqjOzhZxB1baG7tXV+az90NApM\/jD19fGTpsrhBkUM1NDwBAwBAwBQ8AQMAQSCESS8G7YsEH+\/e9\/y9ixY+WWW26R\/\/7v\/5YTTzyxiRxzzDHSrl07OmtSwgj4ZNiP47KmNURnYggYAoaAIWAIGAKGQCYIRJLwVldXC0saSm32NhPgg7ZK6HRWNphf7LT6pUsONBw2oPFMbq791HYUn1zXb\/UZApWHgPXYEDAEDIHKQyCShHeXXXZxSxmeeOIJYU\/ech42XSKghDPffdV28kUwM61f\/VDim+\/+W\/2GgCFgCBgChoAhUH4IhBLeUu\/m1q1b5XOf+5xs2rRJjjrqKLnooovkmmuuaSI33nijbNy4sSS7o8QvWyKnRJA1vwidVJ2G2gZ5SHL9VrIlaO+UWf7RutLpn\/oVbMrX+\/GgnaUNAUPAEDAEDAFDwBBIhUAkCS9Elx+tzZ07132A4rnnnrNPC6ca5RLIU8Kq5DyVS2qbysbyDIECIWDNGAKGgCFgCJQBApEkvD169JAZM2bIwoULUwo22JbBOOW0CzW71P+Qr7nZVyWe6ZBUdVDL1OzSVlVZh1qXVhBMq95CQ8AQMAQMAUPAEMg3AtGuP5KEN9qQm\/eGgCFgCBgChoAhYAgYAoVEwAhvIdFuQVs6Y1oqs5zqR018tjiTrtXEZ3+1juZmmjOp22yLi4C1bggYAoaAIWAIlCICkSS8q1evlqOPPlr23HPPlIINtqUIfCF8qomTUSWWLWmzJkFS63\/gpnXpcodgGzVZtF2TaGObVp80pL2T735DXppfl9TGMgwBQ8AQMAQMgSIhYM2WGAKRJLxt27aVE044Qc4888wmMmzYMLdHb8+ePeXyyy8vuS+tldj4l7Q7wZlfSK46PG32Kpm1oM6+wKaAWGgIGAKGgCFgCBgCSRGIJOHt3Lmz1NbWyg033NBEfv3rX8szzzwju+66q7z44ovSunXrpJ1PlXHvvffK6NGjZfPmzY3M3n77bTnooIOazCyH2TYqGJFEshnbdN1XUqr1aDklrzqLq\/psQ20H0pttHUUvZw4YAoaAIWAIGAKGQEEQiCThbQ6ZPfbYw838vvnmm7Ju3brmzJvkP\/\/883LTTTc10aN4\/\/33CWT69OmNdoh46KGHpGPHji7P\/mSGgJLXoXvVf7VNyfHSdY2XT6ie2rVMME7axBAwBAwBQyBaCJi3hkC+EShLwgtofIGtQ4cOUl1dTTJtufLKK+XCCy9Maj9\/\/nzZb7\/9pH\/\/\/kltopah5LEmvu423\/7r7K+2m+\/2rH5DwBAwBAwBQ8AQqGwEyo7wfvTRR\/Loo48KSxtYZsBniNMdYsguH7Fg9va0004LLQbhZQa5VGdzlUTWxH8A5ndCZ0h1JlVtfRs\/rvbBurRcUO+XzX28cY3qG1r1h7iJIWAIGAKGgCFgCBgCQQQiSXjZeYEdGMJ2aRg4cKBAXNeuXSvXXntto7W2lAkC4KcnTZokc+bMkQMOOMBXJ+K0O3v2bAUT58AAABAASURBVPnnP\/\/ZqF7aSxi1MKIkMhMSp2V8EthCN9IunqptzaMvujxBddoAecSDenSaRxzROoj7Uox+++1b3BAwBAyBgiJgjRkChkDGCESS8KbapSFs5wbVsbNDxgh5BVi\/W1dXJ\/vuu29i\/e4rr7wikGBmk4M\/cPOKuuin7bu7sNWHa2X58hUuvmPHdtm2rWEbrh07dsT1O5zeT7uM2B\/sF63ZFIuJ7NG5WhpstnvxHbH4dmdDPmVcIvYnGCc\/pk7YE9+2bWsiTT7toF+4enNM3+AjeegJG9fbsP4WPflqR4jUtxFWVz0mWkbJMGnqQnwyrBiiz6XwtmDr1q1uHHJZr9W1LW+Y2pjlD9t8HLc2XtEaL44BG7PSHbPFixcLnERlxYp6nsH91kQkkoQ31S4NYTs3qI6dHVoy6Mz8MgPMTLDWw6eLJ0yYkDjIVN9c+OqrrziTuXP\/KatWrZJFa+p3g+hUVU9+IXfo\/7WyztmpngT69evr9VyAEPToyBMRkkKaCPnoqZP0mjVrCJygJ58E9hqvq6trVJ48bMinDHEkGNd+tPukgRjT3qZN9f0j3OWz9ST3H0vXur5TD9Ktdb3+3SX1etpCr35TlvYQ9Cr4hi6X8t577wl+85CTy3qtrlVuzPOBg41Z\/rC18YoWtvkYL+q0c6y0j4MpU6a43aWYgEPGjx+vt0kLYwhEhvBCPpipfemll2JuZ\/6fcpTPvGTzJXbffXe33y\/re5u3FunTp48MHjzEmQ4cuK\/06tUr8eM6XRvMLDZ6Qgy7davfwYA4ek2Tj6BHt9tuuxF19ZEmQT5l9Ad8akMeevKJY6823bvvJqTRk+\/HtTy2fpy60FEGvR\/v1Kl+BwtC6sOGOrEjji15xAmpa\/WHpET679ZQFv32Nt3qM+J\/N+1s6zAkL1fCmOJbz549c153rny0eno1Ghsbs8Z4lPrxUbjxihYupTxuNmalfSyNGTNG2DFKZdy4cfG7pAUgEBnCu3PnTjczpLN+OJ+JUI4n1EzK5NO2T5\/ervrq6taiBBBFu3btCJygr47vMlEds3PK2J96ff3+wtWxfCSmjpHcxnVVx8tUx2wogw0SjJOPvjpmv2ePjkRl1aYdrj4S1bHy1bE8jQfLo0dS6atjdWBDiNTHg\/5Wo461W+0wUTtd0rBy4w6n99uR2D\/s0KncMWOVvL58W8JW9ZmEjEObNm1aVEcm7Zlt2xZjbWPWcgwLeRzaeEVrvDg2bMxKe8z69esnQ4YMScjgwYNjd0j7rwhEhvCqw9dcc40cddRRGQvltI5gmG76+eefdz9WI\/TL8Nqb9NChQwlKRvp2a+t88de7OkUR\/qgPNYGtz1SvvmbjmtZBWeI3PrNYTvn5G0IcnYkhEBUEOGZPvvsNufGZRfLS\/DrnNjrSu3znBTn4R391eS7D\/hgChkDWCHBecY5Nm\/2eIJxjCHH0SNaVW8GSRCAyhLeqqsq9vuQpMxskKceromzKahl9crr++uuFHRvQ8+W1K664Qo477rikuztgF5Sa+LZhyXYeCNrnM82JT\/1B0qm+qa\/YIL695qmOfBWtLyxPbfywJk6Gm7NPle\/n+XG\/HYsbAqWIADdbCC1fD9SHNtIIaXzmmCZer1skY6fNdSQYksyNGpsWihU3BMoSAc4dzhHOFX14ZGKEcwjhvEKIo0fUjnOTsmUJTAV1KjKEl\/WUDz\/8sMycOTNroXxLxpb1tayNGTRokHtlwLZoI0eOlEsuuUT8H7K1pI1UZWsSJHlbKrOKzpu1YH2i\/7ZdWQIKi5QYAtxAIa0IcW6y3Gxxs\/b4foIQ5yZNOGpQT3n8m4fIm98\/QrgOoMeemzBxSDJ1INibGAKGgLi3fJxfkFzONc4PzhWw4TxiT3rOLYRzDiGOHsGO84tzjbIQYOrivCPPJFoIZE54o9W\/rL2FwEJuIbnBSshbuHBhYmuyiy++OGhSkmlOXBzjRCfMh6TThravM8jJ\/FDCOmxA\/Y\/UtO6gfXJ9w9ZowTKWNgSKhQA3YG6gHLcIcW6gnBfccGuP7x8jvP1l3a3HOIIL0b171EAZNqCrI7uQ3to4KUZPPmn6Qz3c2OvbaFgWQZ6JIVApCHBeQVDrz4XFAsnl\/ILMcs5wDiFPjD1ESCO18fOOOHpEz8Ha2PmmBJi6GupeVCmQlkU\/jfCWxTBm3glOfkpxYSBsTjJdntBcfS3JVyKsffDrSrc\/fhmLGwLZIpBpOdYFQnA5drmJctMlROqJa\/9GVWIH0W2kjCVq4zdnbuDkk6Z8LCs+q7VYaIfXssxIoTcxBModAa7\/SkZ5+OP84Rzh3OBcg8ySRp8uFthyfkGAqaM2Rn7R0RbnWD2pNuKbLp7FtDPCW0z0K6htLg50V4kzcSQVeSU\/U9F2KOfHSZsYAoVGgJsuhJObIjdiCCg+cNOtjZFWbpyECHHyshWILzNStbEbsgp1MSNF+8RNDIFyRIBrPW81OM455ziXOAcgupBczo1c9Lu+3v6JNy+k69te7NbS80Cbi3YiWkfJu22Et+SHKDMHOfkoESSW6AotXAxoU0kt8VSiPmsfNNRXSVqPLoUI6lPVbXmGQKER4MYLyYVwciyTxgduxLUxsks8H0LdKsxIcZ7QPmSAMB9tWp2GQLEQ4LziQZLZVu45tbEHPo772tg5RjpffkGiaQdSTTucW\/jBA26+2rR6W4aAEd6W4VcSpTnZcESJIPFsJZu6GspsS6vZmjR3Y0irshRGXIBSZFtWoRGooPY49iC7dLk2dgPmpshME2Ft7EaMvhBSs0tbt0bRJ70QX2bDCtG+tWEI5AsBzjGOY84z4hzrhT6\/6JsS39rYeU6aB1zOMZvtBY3SEiO8pTUeaXuj5JaTPO1COTLUmdZgdVx00BXDJ9pFkvlGnokhUCgEuBHTFkSzNkZwuSmyjIEQfSGF85H1h7XxGzLnKbNh3JQL6Ye1ZQjkCgGOYWZTOY45vjm2mW0lnqs2Mq2nNnae4wPnvPoHGU9Wj+kLj0AkCe\/GjRvlxhtvFD4mkYlQpvAQZ9aivtYvBeLGSYv3NfEZWeK+1MRmj\/x0NvGWEHdtX\/0Mtl8KGAZ9snR5IsAxyKtVBLJLyPHJrG6p9FhvyMyC4Rs+G+ktldExP9JFgPOL45bjl+OY45ljO93y+bTDH\/\/hkuuA+prPdq3u9BCIJOHlM8FPPfWUsK9uJkKZ9GApHStOarxRIkw8XeHkwzYXxC+XdeFTmPhtaL\/RIdirTsOaECKuedgHxb9QBvNKL20eRQUBjjlmm5jNQXTWiRndmhw8FOYSB\/xhlhmSoDNRtuYwlwhbXflCgPOMa7ieX7WxNxbMqHJM56vNbOutjc\/24ht+c32wJQ7Zopm7cpEkvD169JCf\/exnQnj++efL3\/72t8SeuG+99ZaMHz\/efZXt0UcfTejZN3fGjBm5Q66MalJCzMmZr24Voo2g71xoVEecCyUha6xUb6Eh0BIEOJ4guYScP5BchBse0pK681kWX3X2mfPBSG8+0ba6W4rA0nXbBNLINZxjV8+xltabz\/L4CSHXB0v8Z8Y3n21a3akRaJU6uzRzt27dKrfddpscc8wxctVVV0m3bvUfJsBbPhTBhyCOO+44uf322wVb9CaFR4ATnlZ12QLxfAqkg\/q1XeJh8tL8hq+xheWbzhBIFwHIIsIxx80NEolwQ063jmLZ4TMzvbRPH3j1ygwaaRNDoFQQ4LoOWSTUY7aUHyaDuPlLHHg4tnMsiFDh0pEkvJs2bZJFixbJsGHDpLq6uglaVVVVsv\/++zsbbJsYlKCCkxm3OKEJcy3B+vt2a+uaUL1LFOhPS9tWn6knjpfbbD+V+1omlY3lGQKZIMAx9dBrq1yRKN2AncPxPyxvgKjrLBQzaDbbGwfHgqIjwDnGgxgh13oe0AiL7liGDnB94EGYYpxjEF\/iJoVFIJKEt02bNrLLLrvIiy++KNu3b2+CGDqWL3Tv3l3atWvXJL9Yik\/b7+qahqi5SJn80f7osoVMusWFDPuaXdq5z6YSR4cQ17pr4mshVU9eKklln24dqeq3PEOAWVGEYy0KM7rJRgz\/mYWCTGBDn4z0gkRUpTz8Zs0rZJfecH7xYMaxSjqKQh\/0HGNpg51jhR\/FSBLeLl26yAknnCDTp0+X6667zs3kfvTRR478zp07V8aOHSvPPPOMnHXWWdKpU6fCo5qnFpX85ZOw1cSIJ+6nakOXKNTESSj2iOrVT3SFFPU5rP1ZCxqWMWRDzAvZD2urNBHg+OJ1JDcqRGdp7jpzYGk6nKFXOttLMUgvfSVuYggUGgEIIcsYaLf2+H5uL2niURfOMSW9nGNcR6Lepyj5H0nCC8DnnHOOfOMb35BHHnlEjj32WBk4cKDss88+cuKJJ8pf\/\/pX+d73vicjR47E1CSAQLpJbvDYhhFI9JmI1hUkyZnUEWarvmn9YTZBXSa2wbKWrlwEps1eJTc+s1i4USEcy9yMuYmVCyr0SW\/I9JVZtnLpm\/UjGghAdvVhkvOr9vj+0XA8TS+5XjBbjTnXESO9IFEYaVWYZnLfSuvWraW2tlZefvllueuuu2TMmDFOJk+e7HSkP\/OZz+S+4RKvUWcvuXHl0lWtLxOyWJNktrgmPjOsvmbqp5bTesLKZ0OEw+oxnSEAAhz3EECOudrYjBOkkJtWbZndjOkrN+TaWB+JM8tW5qSXbpqUCAIca+VMdhVmriNcP0gb6QWFwkhkCa\/Cs9tuu7nlDSxtQIYPHx6JZQw1Scig9iuTMJd1ZdJuLmwhEtSjBJV4mGi+2ofZpCLCwXLBdFh9pjMEFAF9vc+Pu2pjJBdSqHnlGNJH1hzSN0iv9p+0iSGQDwS4JnOsUXdt7IGrNnaeES9XqYlN\/BjpLezoRp7w5h2uCmiACw3d5AQkzFS0nNaTaXnfXuvSNbea9m38uOYr2fXz\/HgufPPrs3jlIMCxw2tWejxqUC+CihB+VV4bIx50ltlte\/UKEib5QIBzTH+gxjFXW+ZkVzHk\/mWkV9HIf2iEN\/8YF60FTiYa52JCWCzR9tWflvqRTn3aVnNEuKW+WPnyR4BXjvSS2d1yn9mln75APFi+gQ4cjPSChEkuEeB6rmR31KCewjGXy\/pLvS7uVf45pks6St3vKPpnhDeKoxbzmYtELJCa+NII4sWWBp\/q9\/hN5k9N7FUOeWpPPBeiu0Ro\/anq9NsmztqxVPaWV1kI6DFBqDegSrsR64hD8v1ZKDtXFBkLW4oA55cuY4Ds8lahpXVGsTznmJJe3ibZEqL8jGKOCW9+nLRac4+AksJsZkC5SOFRTR7JdkMb9eS5Jk6Ss\/EXX7VcWD0QGi66diMHKROOA44HxJ954qZUqehw3ugNGVz0\/KxUPKzfLUeAY0iPJY6vSiW7iiTXFz3HWELEdUjzLMwNApEkvJ988onbczc3EFRGLVxc6Kn++It4mDSXH1YmmU7rUrIZZqd5XPDI1zLqLzpfgrO4Wk63I2ysAAAQAElEQVT1aptMr\/WrHaH6wLZTpE0qFwGOu0sfnuu+3McxhNQen+U+oGUGIzdklnXQLR4SCU0MgWwQ4DzjGCLkHNM3CNnUVU5lOMeU+PMwYKQ3t6MbScL7wQcfyIgRI+TrX\/+6PPvss8JHJ3ILS\/Fq4+SndS4EhMUS2lciqD5pSF6x\/ErVrvpVEzLz3JBXP2Os9aBHSPMqidCkchFgnSrHA8SOmzBSWyE\/oEln1J8Ye4jUxN62gJOt500HMbMJQ4DJBY4hjiWd1Qyzq0QdSztqYw\/Z9L2cSC\/9KbZEkvDy9bQzzzxTFi5cKBdffLEcdNBBctlll8ns2bNt5reIR5QS5KALEAh0YbOr6IshwRlh9UF91bSFlYXAQ6+tch2uNZLrcAj7owQFwmKkNwwh06VCgPWpvLKH7PKVQsJU9pWYx\/WHh276zhsnQpOWIxBJwtuuXTu54IIL5KWXXpK\/\/OUv8q1vfUveeOMNOeOMM+TQQw+V8ePHC58YZulDyyEqbg1KIkv5oqC+KVnUdCrk1EbLBG0b+t3OZTVnH6ynJj7Lix6hkprYzBRCXEXzkqVVb2G+ESh+\/bw+hMRxjPBqsfgelaYH4GOktzTHptS94i2akl1InZ1nyUeMtymQXu5R9mCZHKdMciJJeP0O9u7dW775zW864jtz5kw30ztnzhw56aSTZMiQIXLDDTfIihUr\/CIWbyECSkZzOWPLTRS3NOQkJx0U1WvbNXFii10wD52JIZAuAsw8YcsrRUKT5AhAVFjugQUPCXZDBgmTVAhwfWbdLjacYwhxkxAE4irW89bEJmrsHIsD0sIg8oQ32P9NmzYlljVs27ZNpk2bJkcddZSMGTNG1qxZEzSPVJoDH4eVcBIPir6qV9tgfjCdqX2wfKp0uj6kqiPXeTVxgszFl7qDWAbT2JiUPwIcD9xUamI3l1EV9HGJlowsWPmklxnyltRnZcsXAc4v1qPSw1EVuNcu\/c5GOMdY9kFZrk\/6UE7aJHMEIk14d+7cKYsWLZIbb7xRjjjiCEdsf\/GLX8j+++8vDz\/8sLz55ptOfvnLX8o\/\/vEP+d73vpc5QgUowcWAZji4CU0aEAjO5GZDSBVfrauhdhEl\/L4uAnFzMQcIcFwwM8nWY3oz5hWinYfpgwtWurwBDME0\/dJmWQkIcEwws0vI8cKsZSX0O1d95G2KnmMsB7EHy+yRjSTh3bp1q0yaNEkOP\/xwOfbYYwVC+\/nPf14mT54sb7\/9ttx5550yaNAg+cxnPuPkmGOOkaOPPlrmzZuXPVIVUpILUr666tftk09fn2nbWg8XUyXD1OfrM6nTCHAmaEXXluMFgsasCXGE48ZuxpmPKTdkHhQoCbEhNDEEFIH878igLZVvyDlWazs3tHiAW7W4hiJUwLKFJ598UvbYYw9HfPnB2q9\/\/WsZPny4tG7dOtQjCPB5550Xmhd1ZSpyx02c\/ikZJI6onhs96aBgr3lqm6ydmvgyAa1D0zWx18OqSzdsrkzQp3TrVbtkfdB8CysDASW6HG\/MnkB09fV8ZSCQ217yAxtqBFd77QoSJiDg\/0iNV\/Ocb+hNMkeg9vj+og+WtnND5vhRIpKEd9ddd5WnnnpKpk+fLqeddpp07NiRvqQUdnBgHW9KI8ssOwQg7kiyjimB1guJprFfuXGHDJn0htgrJNAoL9Htx7gJM3vCusLy6mHhe8ODA63aa1dQMOFaqjP+nF+cZ4ZKyxDgwZJ7FdiyHKtltVVe6UgSXpYqsDVZquFiS7Lly5cLyx9S2RUy79P23V1zuXjK5YB3lUX0TzIMaprMFrd1PdSZWZeI\/amJzx5rCKlVTNAhMbNG\/30d9mT6ZUj7Mvm1OiGfV3K+3uLRRoAHGGYiOR7sJpy7sQTLWnvtmjtA81NTQWrlusmSIRqD7NbGZieJm7QcAd5GUQvXMHubAhLpSyQJL91jWcO1114rBxxwgOy5555NZO+995azzjpLsMO+XCVIBJP1U9emcpNPZpOpPtO6ktmn24dM\/VN77TvpZD6QFxRmeNHxWo7QpDwQ0AcYbsTl0aPS6QXERnGF8PBwUTremSeFQACyy8wuIddbJWiFaLsS2gBTe5uS3UhHkvCyO8Ntt90mU6dOlR49esjJJ58sfH2NL64NGzZMOnToIP3795cJEyZIt27dskMmz6X6dmvrWuCi4CKl8icHfmjfclCVcHJTj4bEw0RxVDsNVU8Z1REPiub59qs27UiYGelNQBHpCOOrYzl0r9K8NkQa4JjzEBwlvbbWMAZIhf3ngZLZR66pSswqDIK8d9d\/m2LnWPpwR5Lwbty40X1Z7cQTT5Snn35abrnlFvniF78o++67r\/Djtccff1yqqqoEu2Q\/YksfovK2hADQw3RIKhcwbIOSqmyyMsE6Mk3XxJc+1MSXNqRbXu2137q0oSZen6apT2d4ifuzxKRNookAN2I8H7pXV+GmQdwk9whAesGY8yzKaw1zj0x518jDJGu4a2LXZdbHE5Z3j4vXO96m2DmWGf6RJLwfffSRrFu3zn1JDUJbXV3tyC5fVNuyZYub3R0xYoT7YRu2mUFS2dZ6geJGlWskauKkknprYhdEQsSPpyLP2GYj2pcar\/1gPcF2fbKLrdZB3CRaCHAThnQhuuZt9OG9otWJCHoL6cVtHjJsaQNIlLdwjWQpA71kht8eKEEiv+KfY3pty2+L0a49koS3qqpKILl8SU3h32uvvWTt2rWiBHfgwIHy73\/\/WzZs2KAmeQjLv8oaj5gWurd+234cP4IEFR3i67UMF2LywkTz9OKs6TBb00UPAYgWN2FIF8L4MivCDTl6vYmWx5x\/ekNmPW+0vDdvM0GA84rzjDKcW8w+EjfJLwKcY7ps5MZnFtuOQs3AHUnC27lzZzeLy\/67ugtD9+7dpa6uTpjlpc9z584liKRwEOM4FxHCfEhNfLbTf4WfaTs1aZBhn4D69dfE20dX49Xjx8lLR5oro31UX9Q+Fb7BGd50\/DCb0kNAZz1qj+8nCASMrX1Kz9My8CikC5AfHjDIYoad0KT8EPDX7XKOlV8PS7dHTNZwbcNDHixT3dewqWSJJOFt06aNnHnmmfL888\/LOeecI+vXr5cBAwbILrvsIldccYVceumlMmXKFGGWt0uXLpEeXz14laQ115kG+3bNmeYsP+hbMK0NKeHUdHNh0D5Zvc3Vkyy\/Aav6HxAms0OvtsRNooEAY8asLsdN7fH9Y4S3v0DAouF9+XipBIix0AeQ8umd9YQlQ8wucp6xbtcQKTwCXN\/0wVJn2gvvRem3WEjCm1M0jj32WGGnBojuZz\/7WbcbwzXXXOOWMPBRCnZvgPhCjnPacMQqU9LIzb8YrnMRbK7dmiSzvanK1Xizwr6dr8+071rWx0p1Okvst2Xx0kaAWSc81BsBcZPCI8A5ZK9dC497IVrkWqkEi4dJZhsL0a610RQB3lxxrtmDZVNsVNNKI1ELq6qq5Mtf\/rL8+Mc\/dtuQ4f+RRx4pL7\/8srz++uvy7LPPyv7774+6bISDmc4sXbeVICGqT0XKuDBRQG2JZyotKRtsS8loUE9a29EQHaJl0CPoENUTTyVaRu1T4aV5auvXC5bMVL00v85XW7zEENDxGTWoFH+gVmJg5dkdiJDO9PLaVccmz81a9XlEgOugT3aZZcxjc1Z1GgjoDDsz7naONQUskoR3zZo1btnCb3\/7WyHud4uvsLG0gdDXWzy3CCh5pFY\/HpZGh\/h23ADRIWGkEr1vH5ZGh\/h2yeLYBUUfHMLa1zz1k4u7lifOBYUbt+osLC0EGCNmOjgedAxLy8PK82bUoJ6JJSWcO4xR5aFQPj3mDYqeY\/owUz69i2ZPuNbVHt\/POc85phM3TlHKfwrkWyQJLzs0sCMDSxgGDx4sfGzixhtvdLsyfPLJJwWCrryagRhoj\/y46rINU9XFq2byOUn9+nlK5aStPb6\/r07EwwhqIjNFhLbI1jB4w9V6g3rK+DJrwfpE0p6iE1CUVISbMQ7pmBI3KT4CECPOezwZ98i\/CUwiiMBDr60SHvq5lnK9jmAXytZl7pt6jv3kL5vKtp\/ZdCyShLdbt27y8MMPy5tvvimTJ0+WQw891H117fjjj5eDDz7Yzf6yjrfcPyuczYAHyzRH7vrGvwgXLJduuia+Prdml7aC+OWeGHuIvPn9I3yVi0OAOWldwvszKvZqmjqCe6iiU7OhGXw9S\/vul9d6GvLaJfxWnYbYKrEibpI3BJqtmDF55I01ct+rdULID2koFHYcoTcpHgKc97T+6pLN8sTczURNIoQAO9hc+YfFzuNRsVl7rtcuYX9KBgHOMe5ri7d1lG37nlIyfhXbkUgSXgWtU6dOMnz4cPnZz37myO+f\/\/xnueyyy+Ttt992pPekk06S1atXq3nJhtysca6l5JI68ik1cfIa1kZNjNCi15C4L7noGxdWfvzCRdavu8bzqybuB\/nJ4uSFido39xpIx4s6\/Dhpk+IgwKvVcY\/Ml8mv1Qkh48J4cswUxyNrNRUCnMfkT3x+rUCgiJuUPgKMFWOGp1yH7YESJAot6bWnM+\/b9h0pSz7qmF6hMreKNOHVsVm3bp288MILctddd8kDDzwgy5Ytk6qqKrdNWVVVlZqVVaikjJt6KXes\/qLYT3iVmQs\/U\/WXPETbqfGIsOoIa+KkGJLkp4mrQJiIN0fUdRywNSkOAoyV\/njm5IEdZeheXd2sfNjbg+J4aK0GEeBBhHFCb0sbQCEawtuTv63Y5s6vXF3To9Hz6HnJOVZ7fD9p9eFaWbytQ\/Q6kAePI0l4+djE7373O\/n6178uBxxwgBx22GFyySWXyD\/+8Q857bTT5NFHH3Xxxx57THbbbTexf00RaI7INS2Rvab2+P7uApl9DalLQqq5eRIms2xJf7UsxIr6fZKrOvSlIpXmh47BHp2rZeLw7vLINwaGLpWpNFxKvb+3nNrfuciDJ7ueuIT9KVkEWCbEul3OM509LFlnzTGHwNFd3peOL90shE5R4X9aRbH\/rM2988475aWXXpJ+\/fq5j0y88847wrrdK6+80q3jrfT9d\/M9rkoCtR1Na6j6QoWsWYJYp9NeTWDmV9M18Znfpeu2ipJaLu6qV50SrHTaMpv8I6Bk6fRDuue\/MWshZwhwXk0+raerDyL10vw6F7c\/pYcA1zx9i3LSvh3lsD5tS8\/JcI8qXssMb8WDEAcgkoSXtbsQ2xNOOEGWL18uY8aMcbO8zPg++OCDsnjxYqnU3Rq4MDG2Sjy5qZBWska8OdGy2Gl54r6etC9qp6GfV6y47++wAV1z4obiSz8RKlUdcZPCIgD2zBDSKrP8hCbRQeA\/ereVcUfXk162UTLSW3pjxzmmZPfUA7rKRYNzcy0tvZ6aR+WOQCQJb7t27eTkk092a3b\/\/ve\/u49N\/OAHP5DOnTvLpEmT5Etf+pLst99+cv7558vGjRtLZgw\/bb+r80WJkkuUyJ+W0LJnpQAAEABJREFU+jQqvoPCsAHdSqRHIsM8kjs02e4NMW+VGNfEZ365wCOxLGGGlzBd0XLp2ptdyxBQvGtis\/PDSujYa1mvKqs0b2aG7lVPoi59eG5ldT4CvWUnGh4qOcds3W4EBsxcTIpAJAmv35uqqirp2bOnW7t73XXXyYQJE4S9eT\/++GO3L++2bdt885KJc\/HAmUxmXrEPE61Lb\/5hNppXEyd1YTaZ6GpiBMO3HxYjl\/zyetSg+tkaP6+YcW6k+Ip\/6ocSXE03F1IeG5Y6KI7UgdTrG44xbg4H\/+ivqE0KgAB400ypHXf4ZJI+AixJ4lzl\/Dr57jfSL2iWOUUgWJmu2+UaaOt2g+hYOmoIRJrwbt68WWbMmCEsbxg0aJAjusTZpeGiiy5yM8C77lo\/qxq1gSlFf7nopfKrufxUZfOVx40UIu7XH\/QzmOami31zs7taTh9aKMdaREJuFNRhklsEFFvW7WqcFoammMEn36T0EdDZQ2YT7fwp\/nhxfulSBh4o\/UmD4ntnHhgCmSMQScLL54RPPfVUOeigg+Tcc8+VZ555xi1hYFuyOXPmuB+z1dbWuh0covyJYS44DKkSK+K5kuzqbGi9peUbasp\/LOhrsrTO2CqBVc9q4rPijEcwT20IySdEXpq\/nsAkxwhwA0Z4sGDNJ9XXxN422M0YJKItjKM+nDLG\/vkU7Z5Fz3uwZwzwHLJbe3z9jhqkTQyBqCIQScIL2F27dpUrrrjCkV2+uPbAAw8IP2LjB23kR1nSIV757B83Hq0\/WVzzox6G9Y+LPf3q1amaIFT8cmrgk2GboVJUchcyLsz+Kfakidur1txhXOyaeHAZGl\/Pq4Sr2D4lbb+MM1gqpOeazryXcXetaxWCQCQJL3vr\/upXv5JvfvObsvfee4s\/i8va3e3bt1fI8Fk3s0UAokRZfbgg3pxAsFjHix3lEeKq0xAdgj2hSW4Q4CZMTRAiPiqBQHYhSehNygMBliHREwgXS1eImxQOAR7WeYPC9Y3zq3AtW0uGQH4RiCThBZKdO3e69bu33HKLbNmyRSC53\/\/+92XgwIHyhS98Qa6++mphjS+2RZKybJaLIB0bGvE1k0Pjs0ijD+9Fd5xo31wi9kfX8CYjxTXeUoeYuRjBBYX8iW5ZNSq+I0iNLWXIH9hFrlmXNkC8dNyL7FJFNM81TGfWRw3q2Winm4oAwDpZ1ghElvDyFbXzzjtPXnzxRfnoo4\/c1mSPPPKIHHrooXL88ccL+VOnTi3rwctl52ri5K25OiGKNTGigTRnW8r5vKZDuKgn81OXNGhfl63fliC1NSF4cbPw62KGyk9bPHsEwBY8a2LHns3oZo9j8Upm1jJjXHt8P1eItdpGeh0Uef3DOeaT3Vpbt5tXvK3ywiPQqvBNtrxFZnT\/8Ic\/yNChQ+Whhx6Sbt26yR\/\/+Efp0aOH3HrrrcJX2M4880y3vnfDhg0tb7AEa+DihFsQAMJk0lx+WDm\/TE2A2EESmX2piRGPsLJR0oWR3eb6pbj37Zb8S0PN1REljErFV13OkAr3UvHV\/MgNAhAuPUchvXru5aZ2qyWIAMtH9KGS63ww39KGQNQRKBnCmwmQEN4lS5a4mVx+pLZ+\/Xp599133a4NvXr1kqqqKjn44INl3bp1bvY3k7rV9t5775XRo0eHLot4\/vnnZc8990wItlqu1EMlDC0hZS0pGxV88FOXNGh\/gzdcxVL1zABTTvXETXKDwEvz61xFkCAXsT8VgQDEi7dKdFZnH4mb5BYByC5rd7nWMaGR29qtNkOgNBCIJOGtqqqS6upq0Y9KLFiwwH1O+IgjjnB6oOWTw9hUVVWRzEggtDfddFNoGfLYHWL69OmycOFCIbznnnukVEivki4uXKEdMGVKBHzclPD6BVLhq8SX17GUWbpuK4ETbih8kEKJm1Pan6QIgCU3YEJEZ54U26QFyyPDeuEh8MTYQwTSyzFgpNcDJkdRrkk3PrPY1cYDpX8NdEr7YwiUCQKRJLzdunVzP0576qmnZO7cufLwww9Lhw4d5LDDDnPDsmjRIkdE+\/fvL507d3a6dP\/w4YoLL7ww1Hz16tVy\/fXXy3HHHef2+MXogAMOkEsuuUSmTZsm5KPLVFKRqEzryrV93xSv7nPdVinVB9n9j94Nyxb0JgD5Uj9Vx\/ipHl1NfBmI6rDnpkJaX82jMwlHAJx4OIDcEPKwgGWlHov0vdIFIgYGPARxLhE3aTkCnGv6Oefa4\/uJLiFpec1WgyGQKwRyV08kCS8zt+eff74sXbpUTjzxRGE9L+Fee+3lljGMHTtW3nvvPTnnnHOkTZs2aaMF2X3uueccWT7ttNOalHv\/\/felrq5ORowY0SiPtcTo33rrrUb6QiZq4mtquYC1tN2aOGFraT1RLD9qUC\/n9uDPdXRhsj81MbwR8pvDnHxmp7Dlhk2auEk4AkpwmdXDAswIa+1HNMBQkTJsQFepjREyOq8EjbhJ9ghwHeKhkpBzrdbOr+zBtJKRQCCShBdkWaP77LPPyuTJk4VlBtddd51bzsBM79lnn+1I67BhwzBNWyZNmiR8qY1Z27BCEF70u+++O0ETmT9\/fhNdmGL58hWyfPlyl7Vjxw7ZsaN+32DiukyDzKDeT5OPYI8Q9wUdQp2N9VvdUhDyfP1u7XYm9NoO+du2NdhTptzlsD5t5ZUrD5FbTu3v1n9v3Vrff392kdlfxQGMkIWrNxMIOzsgJNBhFxyDGfPWJLAm32RbIzz04WDc0T0FOfWArvL\/zttXGJswrHwdO7bomPl6izfGuFTwyGS8xh3dS3gQhaCdfPcbjY6ZUulPlPzgbRPnGtczrnfp+p7JmKVbp9nl5vxcvHixvPLKKwlZsWIFtyKTOAKt4mEkg27dusnw4cOFpQtVVfVrdZnR5cdm\/KisUJ2CAHft2jXt5n7\/+9+5H8RpgfXr61yUk37VqlXuQo5C9cTRaxo7dAh6JEiq1qxZI+g3baonYtgizESjR0ir+Pp2nzSUIY5tJUnrj9a7NwRgyEMOfd\/lszsUKheiQ7hZoHhr4SqC2MPLDqE8CZY6YPP3fy0jmRDGEb3JKneM+jj8bd4yt\/UbuO7VYZt8\/cC28v3\/7CrEfbuwOG91\/DELszFdU8yLhUk248WxwIkEUYP0Fsv3qLc7+c\/\/Ftbtcp5NHN7dXbPS6VM2Y5ZOvWaTm\/NyypQpjlvAgZDx48dzupjEEYg04Y33IXLBV77yVbn55vofxbE8o1u3erLctm1bYZcJdHRK9aTRaxo78hH0CDakVfgaHfr9PtddVS7s3n031wZ5+kqeDF9\/4J69UDn5j336JuwpUynCQwwY9uzZ0\/Xfx7z\/bh2dDiwU900769f77tmjo1BOYv\/Iw2brZxovj8AWvUmvBI6KxYv1Lz7k6H26N8lTm2RhcMyS2Zm+Ke7FwCSb8eJ69OiFB8TOLpG\/rdgmP31pc8bHSTH6Wkptbm\/TTSY+v9Zh+LUj+shJh\/VPG8NsxqyU+l7uvowZM8Zt1cp2rci4cePcONufegSM8Nbj0KK\/zAIyQ5puJX369JbBg4eI\/quubu2i1dXVArEilNi\/6rg+FnX6\/rt1IiorNzbMNmKPuAzvDzrEr4Ps6ngb5JFW8fWf36OegEOIsatEadeunVv\/rX3\/z313U6jEx6Vv\/Ed9ryyunxWH8Cp+vHqlvI4X5aiENPqiSOyhqpTbVRxZS52pn8Exy7S82bd115lC4ZDteP3nvt1Ft8565I218tjbdQX1u1D45KOd1R+K\/M\/9c7kMuTXRrNvNpJ1sxyyTNsw2+\/OwX79+MmTIkIQMHjzYjbX9qUfACG89Dmn95ekWQwguYVAGDBgQVCVNK\/mBFCUz0m2tlFQls8tEr+02V+bN7x+RuKk0Z1tp+ZAx7bPiyfIF1QVDHeNcjmOwjaimFRv8J86rauKKK3ETQyCIAD9iU9LLD69eml+\/LCxoZ+kGBDi\/9AMeowb1jBHe\/g2ZFRqzblcWAkZ4MxhvCG\/Xrl0l+OO0WbNmCfoDDzwwg9qSm+rNngtUcqvc5mibWitpRNOVHnKDAAMw4WZL3Bcdq5r4Dhc1u7R12egRElpO0+jYgWCX77wg3LRJV5LQ94N\/9Feh\/6zHVFz4xbjiV0l4WF8zQ4DzqTa+c4MSucxqqBxrzi2uMYScW3zQo3J6bz01BOoRMMJbj0Naf\/l08ahRo4SPUrAzBIXefvtt4cMT6MlHl0w+bd\/dZdXESZFLZPFHZxNr4qQqrIpUeWH2pgtDoLGOGSXE1wZxTjWLWxMfdx0\/6nnotfofu0H+uBmhqxTRvtNfZnZ1uymIDDoTQ6A5BGqP7+8+SoEdhI7QpCkCuiNDTeyeEbyGNbU2jSFQnghEmvDu3LlT\/vWvf8lPf\/pTueaaa5zcdtttwocn8jVcF198sdsKjY9TsBPEyJEj3Ycn0OerzXzVm4qc5avNKNcLEeOGkU4fFFtIrBJcZi4pi44wKNyUgrpyTYMBJBc8dbYJHelR8b2Qy7Xv1q\/cIuB\/iY03BbmtPfq1sa81OzJwbt115kD3G4SsemWFDIGIIxBZwrtp0ybhF4hf\/vKXHQHla2vInXfe6bYq+\/73vy\/bt9fvb5vNGLEnL79y7NixY5PibIW2cOFCUYki2W3SKVNkhUBNfNZWC9fEZlCIawjZhcipTvWqg\/SRh6iOeLmL9psHg1GDegqkF2zshlzuI5+f\/jHTS80cVzbTCxL1wpsjyC4pMOKhnbiJIVCJCESS8DKzywcn+PAEX0fjYxFKPolfdtll8sgjj8gDDzxQMWMKcWhJZyEbLSkvIhVZ3McdDBEfCP3hYVCPTZDgBtPYlKvocgZuwvQR0ssPJe2GDBommSLAcaOv6iF5RnrF7WetONQe3084xzLF1ewNgXJCIJKEd\/369fKnP\/1JzjzzTGF2tVOn+u26GBjizPyypvapp56SDRs2oK5Y8QlZxYKQx45zow2rXglukMTqeKBHKKu2zAaTrgSppL5WwniWQh85F0uH9BYXEa4t\/CAULyC7tcfbjgxgYVLZCESS8PJVsS1btsj+++8vVVVVTUawqqpKBg0aJOvWrRM+g9jEoIIVSq6AwI+TNskOAcXRv6nUxJc6zFpQ5ypVolsTX\/LgEz7N4ybljON\/2GopqItnRSpgxo2bL2sJcZw+IWABSUFnYgjkAgGOJ5\/06jGXi7qjUgfnFucb\/jKr61+X0JkYApWKQCQJb5cuXWTgwIEyY8aM0HW6rN198cUXnQ22pTi4peQTxKOU\/ImaL6w7rU3yypCbD\/0Jw3jWgvVkCTdpzVd7QrZa4sZF3BlG8A++81qVkLWExPXHeUr0I9gtc7mEEeB8UtLLMceDYwm7m1PXOM+4blApZJe18cRNDAFDQCSShLdNmzYyYcIEtxsDyxfYleGTTz4R1va+9957cv311zsy\/D\/\/8z9uScPq1asFWbNmjY25h4CSLE9l0SwQ4AZbG3hlGCRzihU1+KYAABAASURBVLWGS9dtTdmSkkKMojxLxQ2YPugOFcz2IuiCmKEzqUgEct5pzkklvRDASiC9nGv0ldDIbs4PKauwDBCIJOGFvI4ePVr+8Y9\/yNNPPy3HHnus7L333rLXXnvJkUceKVOnTpW1a9fK+eefn\/jEHp\/b++pXv1qSQ6bkR8lQoZysib92D5KzQrVfzu3UxJcuBPtYE8ecmxJCPjodA10CoXnkK0EkHjVRsg4BUdJL32pi+KCLWn\/M3+ggwPFVG3vzgscQwXImvZxT9JGQ88xmdhl1E0OgMQLRILyNfRZ+mDZx4kRhC7JMhDKBqiKVhCTgMBc1wmSidsnyVa8ka\/ThvVRlYY4QCI7B0L26NalZx5FxCNprXpNCEVXQf70J01edfYtod8ztiCBQG3vzMnSvrs5bCGE5kl6uFfSNkHPribGHuP7aH0PAEGiMQCQJb7t27eSYY46RE088MSOhTOPul3+KC6D20o+jYwaEraB4\/UXaJLcI+HhrHHJLK9yc\/B+uab7O9vt5ak+IUJaw1ITZXDb+V1KBn8xY0zeONUKON4R4qfkfFX\/Mz8wQgADqNQ5iyHGaWQ2la805Rp8IOac4t0rXW\/PMECguApEhvKzRZQ0uwo\/SCFnakEqwoVxxIS5M61zssmkp23LZtFVpZXRmCYwRv\/8QWm5S6OrJYDuibu9M9AhltA6IozOI\/WF9Lz9mU2IZUxX9P\/7yAyH85AaMQ+gIleQTp0+EJoZAIRHg7QJCmxynPJgRj7Kw1InrAOcZhN7IbpRH03zPEoGMikWG8H7wwQfCGlxk3rx5Ls663FSCLeUyQqSEjLmQ4U4ykuATCexMSgsBXqcydv7re8gtXjY3ttgglCdUoRw3bMJLH56r6qKHkHDfCW7G\/i4Ufp7FDYFiIAAp1HORB7Mok15mqdnxBBxrj+\/nvlRI3MQQMASSIxAZwqvrdlmH27t3byFsbv0uNpRL3v3SzFGSwyxgaXpoXqWDAOPIrAthc\/b68AKRRbBHp2V1qYPmkU8cIV5s0dlmSAW+8CU11Q0NWb+MTUHEGjEEPAR44Iwy6eV8h+zy0Eu3ILu1x9tHJcDCxBBoDoHIEF6WMbz22mtue7Fu3bqltYaXNbus920OhHLOV8JUzn2MWt\/8Menbra1zX3U85CAo0dV4uzqg01lT4ggzVYQq3BA1no+Q+nmNqrNLtIEOP2p2aetmmghJI8QhGdiZGAKlgADHIw+i+MIxyvGsD2foSlU4z3iTAtnlvKqNzezWGtkt1eEqSb8q3anIEN5t27bJU089JbNnzy6bMeOiRWe4kBHmW5Rc5bsdqz81AkPjvxrHipsvoX8s6Gyu6shXCR4raks+ywi4eftkFH0uhRsuPtCWtkOaNvT40lledH6ctIkhUAoIcG5BejkXOX5Zd87MaSn4FuZDg4+LBd\/rP3ZjM7thWJnOEEiGQGQIb7IORFFft+zjKLptPucBgaEhr\/u5udFUTWx2V0kkM1HoNE+JpKbJYxkBIWQUIZ6ZNFhTb1gd\/kyY5uuss5J3Zp1q3exTP6m1GagGUC1WUgjUxN5IsINDbexYxTFmTnlg5NgnXQqCLxBx9QufWZKh51op+Gg+GAJRQcAIb4FHqvOWVlL37Fa5+7B\/SufN+YFfSVKBu2bNpYnAqEH1+x5z8\/JvXKSpQskt46j53PjI85c7kFYhX8uhe2n+egIn5PEDHW6cTuH9IS9Mj44ZXMqpOba0gZ8IekivkmCfvNfGiC6CjYkhUMoIcJxCIjmmOcaZ7dVjuph+qy8QcXyrjRFzZqWJF9OvimnbOlp2COSHceURpuXLl8uTTz6ZlrzwwguydevWPHqTWdWftu8uGzt8Kp\/8fbFsiM3y9pkROfgz67BZhyIAieXGhfgGEFzS3OgIgzc29Ah6Jc0QUGzRE4YJyxCw48YJiVUbynBzRw9x9fWappzGsccGP0cN6klUmFXGhgR+EZoYAlFDgHMS0jt0r65ua0DOCx729JgvZH9okwfO4KxubewhspB+WFuGQLkhEDnGxQ\/XLrvsMklH2KVh06ZNJTVmzPC2WlI\/w9clNtubiXPZEopsy2XiWwXatqjLYWPi64gjNKIhN0LSSFCnywq4YZOvJJW4P1ul5BQ99SHEg0QYnYrOFmsbkIOaXer3Ddb6amKvhxEtY6EhEDUEOH79JQ4c25BOyKeeJ\/nsE23QFmSbh1D8qbVZ3XxCbnVXGAKRI7yHH3542p8UhvCW2rZkR77VPnGI9X2\/dSJOZMOy7fLHby0jmrFwcWyuUDo2zdVh+flDYNiAhs8PM4uqLWlcCaemdTy5UartsAH1n1HVNCE3bkLEt9X60CNKklXPzdbXK3Eeulc3YYZX28eGNKGJIRB1BGpjM6m8famNkU36AvmE+OZrxpdzEqJLG7RFmnMrWj9MAykTQ6C0EYgc4e3Tp4+ceOKJaUkpbkuWbFaXmd+O39smb\/12vay7ZENOj5ogQcpp5VZZXhDwiSs3PxpRQurnoecG6ZNRtUePYIMOIa46LdOgr1\/+o3qIreahU+Ks7XNDpj7Ibm2MJBA3MQTKAQGOe45pn\/hy\/ENKEQhqS\/rJOcg5xZuV+voWu+pol6UVtKvnmcuwP4aAIdBiBCJHeFvc4xKqYP4XZya8YW3vst23uzSkGALsEjn4w0U0B9W0qAor3DwCEMehe3V12w4Njc2iagkdP26SqiMcGrMl5IdsCHFsVT9rQZ1bj4iehx5fj458Qm7shFp\/eF3rMXG+uUjsDzdkbsx3jxoYS9l\/Q6D8EOB84vzgOK+Nz\/hynjATu8t3XpB6srpIIMDokyFAHg+sEFxmiinH0gV0tMG5b0Q3GXqmNwRyg4AR3tzgmFUt1bvXf3RAC3fe\/BmNNgpZ6sCuDn+5+f1Gei6UjRTNJIZ6JKoZU8suEgKsIeTmCplUF2ri62U1reNYs0v98cMaW26opBG1I9TlCdQ3LL5kQu3Jx56bLXFuvtSDoEe0DDd4bCDOhCrYaNxCQ6BcEeA4r429xeDcrI0R36Hxh03OFc4NBBILCQ4T8sZOmyucY7MW1DmY6uvsJxBdHhqHhSxHcob2xxAwBHKCQGQIb9u2beWEE06QQYMG5aTjxapk+9n\/Jxf\/YqSTQQe9KMzkctHcb2EbYWZX4v+6eOS3972ful0d\/jLpfQmu+42bpwy4mCJ2QU0JU8lm+iSTm6SOY02cCHP8+M4P84gtr03JU1viiJbRumvi5FlnfVWvZJgySG3spk9oYghUIgKcJ5wDPJiuu\/UY8QkweckwIQ+SXBsjyxDchrL9G701SVbe9IaAIdByBCJDeDt37iy1tbVyxhlntLzXRa7hwe\/+Qu69YLqsmXBVwpMu3p683Q98QljeoKSEuBpCkDVOyKwvQlyFi6vGNQwSF9VbWPoIDIvN\/HCzZFy5WQY9DpLUYD5pCKweA\/WzTPVLFKibfOon5NUsoeqJc5MmpH1fj87EEKhkBDgnamMPgRBgyC9ENkzIwwZbO4cq+YixvhcTgcgQ3mKClMu2d733JNn0QY9ElV28mVxVtupYrVE3A7zHjvcS6S4eMa57bqsw64vwMYuEkUXKDgFm6LlpcoPVzimB1bTeSFVfT2zrXLbmuUTsD3mxQIbGl7lovfqQpXpsuEnTNkLaxBAoJALWliFgCBgCuUDACG8uUMxBHS8fuFX6PLDQLXU47OSHG9W4y9LPJdKdtzSs8+38+s6EvmZox0TcIuWHgBLSYM9UTwgxDeaTJo8Q0bgSW3SIT3CxCRJkdNiZGAKGgCFgCBgCUUSgDAhvtGD31+n+v8u2uKUL9GBw139Kz9+\/6pY6PHX1FDnyrfqN\/bt4M8A1R3aQp4\/YjHkT+eBv\/2qiM0X5IwDJhYzqFmHa46HxH9WQ1hlf4r6eckpsCTVPQ+xNDAFDwBAwBAyBckDACG\/ERhGSoi77uzp0aLtS1aI2N+z+liCJjFhkw7KPZeqpC2TJy+HEOWZi\/yOEAGSWpQYQVt9tTUNeIcWaN2pQ\/Vf+SJNHqMIaQ5ZOIKqzMGIImLuGgCFgCBgCoQgY4Q2FJX\/Kj793v+xz5J9dA\/9zZwdhdwYS\/X53uLz+xCiiTnTpAj9YO+XK7znd0pe3yBH\/77Muzp\/\/nPhdAidb7zrYhfzhR2w+0YXkoidkezPqefDUhW7nB\/Qm5YcAJBfiCon1ewcR5odvCPl+HvFRg3oSmBgChoAhYAgYApFGIOi8Ed4gInlO9267tkkLutG\/ZnzpwG\/InEHzXTLM3mWk+KNkOYWJy1oSI9AuEvsDGYYo86W3WNL+lwECycgrpBcpgy5aFwwBQ8AQMAQMgbQQMMKbFky5Meq8pZXbimzey19KWWHH\/+opq6obZto6dV+dsNdlDNT1xvibEnoi6Aj9dcKk65Z9TCAakmA98IFndCPq5IlvLXM7PvwxFi6x5Q4OE\/tT7ghY\/wwBQ8AQMAQqBQEjvCU40o9P+rGM+Gv9rgufe62vsG9vKjdZvzvwrGtlY4dPndmhl94m+nGLCyd9Tz53ZH1de3z+Hae\/9NGfyNmP7eVs9Q\/LHDTetW\/Dsgl0zPwixE0MAUPAEDAEDAFDoMwQqIDuGOEt0iB32nW1XHjbifLunh85D+Z\/caYM3meCiyf7c9jJ02TN2X9Ilh2q\/\/TD5Ql96+5D5I1JB8pdp13tfszGbC6ZLGcgRCDPCHGEJQ7s84uoPXoTQ8AQMAQMAUPAEDAEooKAEd4SGamh+yyR3U+un6HFpeCyBHTIHp\/Uf4SC2dzP\/LqTdInNxm7ZtofMffCH7iMV2OzcsL9sm1P\/wQHSKszSvjrvek2GhtT19PUNe\/2+9fC6hN0BZzYsgUgoLVIpCFg\/DQFDwBAwBAyByCJghLeAQ1e3Wzfp8OWr5MyjDpKzfnqB+F9U23r7OcJShlTusItD+19+M2FywUdfkg7t\/55I6569C+4ZIQ\/cNcPt6fvLa6cIRBcjfyaXtEr7dn93\/mh6ifdjtk8\/XKHqJuGmFx+Qpd\/sL1vffbFJnikMAUPAEDAEDIHyRMB6FUUEWkXR6Sj7\/Pkviczv8l1HRu+9YHriAxN+n7rd01V+O3yjU31y1Btu3a1LpPkn2exwsLj\/IzaWWGi+T4z9uL+V2ZKXN8udZxwqU38\/XXasWaxFLTQEDAFDwBAwBAwBQ6DkEDDCW8Ah6ft+a2FnBWZq\/WaDn3n185LFqQvCvHLe\/gkT3aUhoYhHlLT6BPekS2rlKz9a5izm\/mkfR8Bdwvvz6YfL3cyvT4Y1e+msLRqVN2d\/udGevuvmLRT2+4UUJ4wqLGLdNQQMAUPAEDAEDIHSQcAIb5HHovOW+vWy7S5\/wM3kfnvmn+WAAzsnvOr2f4c1IaOdt7RKrNfF8MAzusk17x\/ofgBHHjoEonr+j46XERM+ISmnXnuva4MdHNixIWxmlo9iXHDtbc5e\/\/jboilpXvevhZrttjNbEl8GAdm996jNjgD7M8IJY4vMqMiiAAAQAElEQVQYAoaAIWAIVBIC1ldDoCQQMMJbpGHYY593HPn0PzCx6YMecttRX5KZQ5fKBX+o\/4GYvzyBL651u\/SHwg\/WwtweuldXYUuxbz441tXNOmHf7pPYjC1tMDM8+dtPyq+vHOKyNyyt36eXRKdd30+sLYbcYuvPIm9Yth0z2bC0PnSJ2J8N8TrayeOJL8nF1KJEmLiJIWAIGAKGgCFgCBgCxUDACG8xUI+1ueO9bVL3m8WxWPr\/Nz\/7nny8YHPSAnxGlk\/GshSBj1tAVn\/5\/WdEd134zbjrGu3puyH+QQq\/QpZb3HfVZUIe5Fjz2A2CWWRmk9EdfeZ4OebcO4g6wd5FYn8gzbHA\/V86q8HfrfNul9UvT3P6Rn8sYQgYAoaAIWAIGAKGQB4RMMKbR3CDVS\/bfbucdOc3HFFcveEweXjmHDlo9gBntmbCVY3IqFPG\/mw\/+\/+k\/x97uO3H\/vzW\/8qW\/7vJLWdg\/152e9j\/c\/cIe+XesPtbbjcGZoQhui\/8alysdP1\/ZmrrY+F\/jzhnq9s5okeX1xMGwTJd+rZO5EGoWS6x94HPCssfIMIn\/ayvy2+3z+Vy2CkPuxnmcy49OrGcgjLsFvHLU\/dz63udsf0xBAwBQ8AQaISAJQwBQyA\/CBjhzQ+uGdXqr7ulIGts9QMTfdquleNXHdFk+7HBXf8pXb\/WT9oe2JUiCQkS1URGishHy+8W6qru2baR1afelmR8iU1JNUaQapZF\/OKH3yaZEEgtecjfXzg7sWUZyyn2GfpnZ8dsMHYuYX8MAUPAEDAEDAFDwBDIMwJGePMMsF\/9fgvbyB8v+1\/xZ1\/9fOLMnBJCdAmbSr3mMzMPcT9mYwlCvab5vxBNtYJUX3TTnS752A8vdnUF1+r27Pt\/bqY22dIFVzj2ByI89dQFsZi4pRAuEvvz1rsXuK3LYlHZ8cErsmltD6JO+HKbi8T\/+L7FVRYYAoaAIWAIGAKGgCGQEwSM8OYExuwrYQmCX5p1s8yOtr3jQqfmgxSkfTLqMrw\/fHKY3Ri+OH536XtwnSOpfLYYUouc\/dhezlp\/zPb1MV9x64c\/3bzD6f0\/XzrwG3L2V0aKrtX18zQeNou8If5jtrA8yq2cd4DwQz0l9Og2LPtYPv1wuVvXyzZmyJKXG9b8YmNiCBgChkBSBCzDEDAEDIE0ETDCmyZQuTZj+y\/I6EsnLXS7Lkw66wM5alZNymbYdWG3628S1gJ3jm9nRgHdNixsmzHyVSCXxH\/\/2zvd+mGWI0A60al0\/K+esvvVN7vk64+PcjO\/qWaknaH3Z8cHr3qp+ihtbI\/N8ELcIfRoWfvLD+Fate8jT313Jyo3O\/y5Izu6uP0xBAwBQ8AQMAQMAUMgVwiUO+HNFU55q6e5pQvJGm5\/6jNuJpflBhBSfzeGsDKQTggnsmXbHmEmTvf4pB\/LLYd93v0Azinif5hFPu+7x8qxFy1w64kh68wix7MdWYVQs3UaM8mQWc0jXOl9IKP\/7o\/Ldu\/rbFs+PBQTJ4tv+Z4L7Y8hYAgYAoaAIWAIGAK5QsAIb66QTKMedlZgJvfsuZ+Xrzx2uUzY8pS8Wrev23Xhygd3dfvvpqrmwe\/+QtjNgbXAyezeeaKDm5XlB2WQW4Q1tsnsWYJw0thzHHk+cL9fJDML1X+6ZYcrB\/G95M7HGtl8fdIrbucHncltvWv9nr8YdWi7QnasbtiSDTKOHmm14CECJ2vfmCkzfvg3t6uDLXVwkNgfQ6AFCFhRQ8AQMAQqFwEjvEUY+6UvbxZIaO\/Lt8sZz3du5EGnXVc7Erlr7Q1Ov+q4hx1xZFbUKfL4h7a1ekgooml+HHf\/T\/8k94\/pKe8tG+FINbPKzAhjw3IKdmKAYP\/8rLvdj9X2fGSn7Db2V2QL64vZvmzg8HnyzpJLZOrvp7vt1MhEf\/gZfyHq9LqDw79\/O0Vm3d3azR77W6Y5Q\/tjCBgChoAhYAgYAoZAmgg0IrxpljGzLBE48q324mZyv7EyUUPnzZ+RLjFRRbsdS90Pyra9VedUX\/yP5W7LsM\/u1Xhta+upX3akk+UMzjD2h9la\/eJZLNnoP3mq4MdjLEm44vV\/CWtmmTmGqM6a2Xg5wVFnXOXI9\/mPvatFXbgh\/gM1EpDkVh2ricqHWxuWJqB45+aTHbFnn2DSCPaECL5ueWO8rHu8v1smgQ7ZEP8xm7+rw2vXNp5Bxs7EEDAEDAFDwBAwBAyBdBBolY6R2eQXAX6ExgcpIKGn395AhmmVJQyQUb6cRjpMWDdL2a\/8aFmj7MNOnuYI69hfPZ7YvYH1taff8UYjOwimKphtRfRDEuj9L6f5tuThF\/7xAzifVJOns7\/Eddsyf20vdS17s2sT4o49s8Ur\/3UAUSezXvqG\/PFb9f1jrTDiMuyPIZAfBKxWQ8AQMAQMgTJCwAhvkQZTyah+YELd+PWVQ9wOCvpFNdUTQmrZfoy1wKTD5PAz\/+JI7rinfyNDzq1fJ\/up9wEJyrzwq285ksmP0\/zZV\/J8eXzSj5zdbUd9yVcLM7ONFF7inEvGuSUYh37pNwnt+nkLhd0bbtj9LXlu0mCnP+MXa9wyB38W12XE\/ihxhsjHkon\/qv9o2e+k7vkvyur7Dkl82CJhZBFDwBAwBAwBQ8AQyCEC5VGVEd4CjuOn7Xdt0pouB5j3+Ha3lIEfgjUxakax+dn3XNmg2W\/GXecI6z1jLxH9MRvra32SGySvkFKEGVmfKEM+9cdphw6f5Eg1SyL6HFy\/9IK2t867Q9oeuth9\/S24x69PbMe+vq\/sdcqxwmwvs7yURXgI+Pq5Q2XkuL\/ECPIrqOT0C0a5kD8b4ksplry8haRU96qTD6Zc4OL2xxAwBAwBQ8AQMAQMgWQIGOFNhkwe9EsOX+aIIvvp8iMwlgLseu9J0vf91vLv535dP7PrbRnG2t5Dbr7KlTn48F858soPxY58q538dvhGYbcHyOOf3\/pfV5YflD1382Bnd8eIr8mGZR8nehEktomMWIQlDMhZj+0ZSzX81z1z0XRsdaBsfGy56I\/T0LVq10dOGf895x8fs\/jckR1QOznmpyOFOpEjv7ld\/Logukqq2\/T9qrPXP6061K8H1nSfrzUQWvrz6YfLZemsza6P4PeXp8e4XRzU3sLiIGCtGgKGgCFgCBgCpYyAEd4Cjg6znnW\/WSwbH13eqNXOWxqGQX9Qph+YUMPqHv00WpTwoZvOcKSatbp\/+OEljnDevF97mfv8550\/zFSvnHeA0z9w1wz5\/ff7Or3+YTmGxiG7xDfEZmyPOGerWwLBzC4PAZS984xDY2R9u6vLX04R3PqMOha9f0rMtp7YQ4bZyowZ7A0e2cfOxBAoJQSCxydpjlvWqW968QG3VGf72lfcftjoeOPi+89OJtgjlPXzLG4IGAINCHB+IEtebviKJ+cUwnmFcK6pkEbI11r0fPPr0LyQ0FQlikAD0ypRB8vJrYNmD3Ck8YnX\/i\/RrS4e2UXp72JAWoVZUY133vIZueAP3WTqwH+lPbvJCa\/l9StvfEjiwDO7uZ0UIKEPnrpQTWIkcruMvuq3joxCNP32E0axyAu\/GueI6QN3\/zbxxbSYWpa+vCVR7zM\/7C2fLPtPV9dlv\/072U7waXt8P95WbQ9zSxxcRuzP6g\/Oj\/2t\/89DAL7q7LL\/Y7Z6C4n5+7GwtveP45a7H7fdfdg\/xS5Oio6FpYIAxzw3T45PluNseunbsnHWKHcec4OFwK65+1xZNfEY2fTyKPnLpPfd9n2cT34f0GOPLPnfkfLebbslzjfqRs\/54pexuCFQ7ghwfnEOMfmBcC6ocH\/jtx\/sCoQNwnmFcK6pkEbIxxbR8406uFdunXe7QJBpj\/PZ7jXROLKM8JbYOOmuB+zOwAcm3hh\/kyOUf7pnaFJPWSLBDOolMzvKiAmfOGJ59ldGyrinf+OWG1w46XvCrgtnHnWQyzvpZ30SdekPwVD0Obiu3v62E4UPR+xs+4DbEo08n4gvf7MrKic13jIGf9kCJNUZxP4siZFfZoZj0Rgp\/T1BQn57yUr3EPDAPXc40prICESYQebrbJ9uXS6sJw5kC\/3YvvbVRnlcnIJ2JZE2JyoCAW6YU09d4IgoN0glttw8AeC+qy6T7ev+kFivjg55eOYcabvffwrbBZJW4eaqcT9kGZC\/bSF2tF333Bdl9S96ydJv9nc+cGNG\/LIWNwTKAQEe8CCikFviEFvk0OE3N+oeJLiRIp6o3nWIICyTi6tSBlvn3eEeSJc9\/BX3UPpgbLKI9re++6Irx\/nnIvanpBAwwpvhcLz99tty0EEHyZ577tlIRo8eLZs3N7wyaa5aPgkMSd02brKw6wJxZI+914YWhXCSzw\/H5gyaH2qDUmdMg0sgmP3s+rV+0v7I7jLv5WMdieZDEn\/8VuPlFdSBcGHg5EfuuexUVBlJl76tE\/bcgElwI5\/2k3FEnZw\/sVb8H7M5ZfzPWw+vi8dEVs7bX1i7DCl+81dvOz27VYAHIeSe\/YQhw2SCFSGibRM3MQTyhQDHGTc5zputsdmfqTGiy42XmSLa3Bq7Qe74oP6HmKQR3pqwhr1114vc2w3SX7xyd2Et\/R4TX5Cjr\/sPIX3Sz\/o6HfmUWxJ7NYv+wDO6CbLHKXdIh0NuFv\/hEzseEqt3b5t4kIRoI3XPni5bYz5CDowAg5RJMRHIpG09zzi3IJhalnNP4xq2at8ndn7skji3OD+6Dv+L7HLKIjcBpOcV51vnodMEIa5CPraInm\/UwXnYbp9xjiA\/PunH2pwLeTOz8PQq95YR\/+wcc7CUzB8jvBkOxfvvv+9KTJ8+XRYuXJiQhx56SDp2bPxxCGeY5A8nY1iWTxS7bE4+PP858btuNnb8ux+6mSCI6U\/7x0hi9ZsCsW17YNdG1b\/28BcTJDfZTa5V+97OZvK3n5S7zz0lUZ4TXH+c9t1Fu8i3\/7zOzRSzBdkXhv8rYadLJc4Ze4a7UScyvEin7qud3+f\/6Hj3AzgvS5iphtBT79mP7ZX40Rs3dbWDUNDXLiPudCpu6DvWLBb0rCcmT2eX8Rshjw9cfLTsd66M\/TEEcoUAN2DILTc2bsLr\/nBEjEzeISecO7xRE9w0Ox05zd18uWlyU+UHpxDVTkd9V4gjfJGQhzcKc9yTJlQdeuLouSEju+yzp0CcOWeog7rRc3OnzZmr\/kSxhHy67XXnI75DgLkx04eEgUUMgRJCgOMUd7jOf\/bDnzsyGUZwseF6z\/nF+cbxz\/nFOYFwfuh9l3MK4VxCKIsQVyEfHaLnG3VQV7t9LncEGSKM0CZ2vJkhVMF3zjFN4zcPrJq2sLAIJGdUhfUjMq3Nnz9f9ttvP+nfv3\/GPq\/\/8uuO7F1830WxJ89urjyzMCxdgKghc\/\/cxen5w1pdQoTlAuRDRk9esL8M3au+PHm+\/P67JzvSeu8lP5AHLviGi0++8sduHaBvp\/E99nnb+XTpoz8RTmbVB0N2ZEDHByG2LjhXB8Hd6wAAEABJREFUINWtd+snfCqYmdZvz\/yzHDu2frkCP85TYsva2xGnNszqaj28hoWI0i\/qRZiZBY+Pl7zkvr7GKykRcSSBfOSDDycSCHsIgwfCzPFT390pO3a9weXNe\/lLcuiXfiNf7j\/QpbfHZtYgu5Bep7A\/hkCWCHAD44aFcDxBdHUWlyo5frmpcqPl5svNkJ1KyGvdfYgjtpxn3FTR5Vpok7q5WeOHtslNGl\/Q6+yU3zZ94HzzZ37pq29jcUOgUAhw7EEOeRDjHFs58Ri3NGfTzBsbuUAeCs4xHvQ4zjm\/0BVCIMIIbeIDb2Y2f3WdKAEO+sBD8YOnLnRr9pNNPAXLWDp3CBjhzRBLCO8ee+yR0WxuWBPsPgBZ++yPzxN\/JnfIuYscAT3hJ2Pk6SM2y6SzPnDbjzFjo\/UM26urfLJiX2HHh81\/fkbVaYf+LPKWDw+VzU+LrL394SYzrn6FM397kyPP7MyA3wjrcjcs254weyH+QYtfT\/m9vPt4B6eHAPT5j3+6Pp19+mAZPbVGeA0Ecf\/1r2a5NbfnXHq08MO4F+I\/gHvoqb+5ttgzmNmnDfEdF7igHHBmONGnMbUj3uVzsdnuWIQLzEdLf+\/axGe\/vli2\/TcE0kKAY4tjiRssNyyEhyjWxVMBRFPJJDNL9TOuewo3Q\/KLLfiHL1xHdHYKgoDP+EY+D6Bb593hlkAwK0VfdQkRNiYlgEAFuMBxODVGCjnfeBCjyyxnq45NsFR9PMJNgPDwxvE79vV9yXbCg56LFPkP5xQEGN\/wkftW0CWuJ5xbW+fdLuvmLXS7sQRtLJ17BIzwZoDp6tWrZfbs2fLPf\/6z0frdK6+8Mq1a+v3ucEfk7r3ovrTs1Wjbtm3SoVeVJmXtog+FdbXMbBImMgKRHTt2JDScYJr4wvB\/OQJ6\/o+Olz0H\/l2m\/n66++HYH6\/b6S4mascSA2xG37dY\/Lo0n7DTXvXbj328cXHsRtlAfjv0qj+0WHuL6Ezujk49hZlcyiKql116x9puWPdLnsreI9vJd5Z8XsAhmR9gwg\/XtAzkGXywpw3WAWseN\/IXfrLC1UedyeSjjz6SrVu3NmuXrLzptxUcu1yPGcc1xLZ67feavCXZNqdOqvcdID36r5KvTNtDht+0m1R97uJEn9vs9mkiXorHwu6HVjufz3upn3wjJm2GPRfz\/yI9TVy44PE\/ua8abvnbJJl66gKHQS77kuvxyqVvVlf4+ZvrMfsgdj9jtpPrMhMSvGkYde2X3fHHH+4lPQ79jPS4da50+sY97lgddFkX4fgt5THi\/A\/6OOjShje43F+3xh4w7z1qs\/CASd8h+Qtm1mV93Vi8eLG88sorCVmxYgUQmsQRqGcl8YQFqRFg\/W5dXZ3su+++ibW7HFyQ4Ex\/tOa3VLffSkdAWb+699mfJrL4IMWVD+4qvztkqfz29KUJPSQukYhFzn\/oWlf+7Btvkg8WbY1p6v8f++N1wgwUcs70R53NOdeeJbsOr1\/fy7KCtQs\/rDeO\/W3V7SM540fnOruvf+eMBDGtW79eyIuZNPm\/aNUgR+LvPuk89+MyNdj4UQe3tljTrz8+ykVXrVrlQv0DMWXv3d\/fcoS0+vxxqm4Ubv\/7ufLhswPl1s\/9SyCqo17pJMgpj3VI2G2IzTQPvbWV0yeUsciSl7fI+u7fisUa\/3\/v3xvl\/SVvyMZZZ0rd7Mtl1T9mC76pvPfee7JmzRphzFVn4apGGJUaHrkas5d\/udSROwjefUPmC0sXIL0s3eFhref+a4Qb8EdH\/Fw+Oe9ZOXba52Vr9QcljU2qsVLfV2\/4jNR1PkM4r+gfZww\/rmXG7Y4T\/kuYbeOGDCZznlze4v7marxS9c3ycnvO5nrMuAZzPP1l0vuxCZP4fuqf7Skfdzs+NsnTSvY\/\/7PueDz25+1bfLyVwrEwIHZ\/13vXpn0ekm09xnCaJYRlHJ22nyMf\/OOBrPo7ZcoUgYuojB8\/PlG3RUSM8GZwFBxwwAEyZ84cmTRpUqJUjx49ZMKECYknqkRGMxE+nwsR7fDlq2RDh08S1n\/7+cGOQD519RRhba9msCbo7AvGOgI7YuwSVcdmRT8rbTr3c+kO3TvK1+74gSOsF9zxMxlw7OFOz58On6kjkFbbRNq02s+1wbKCd144xOn506ljJ2nd5SLZ+NhyWflKV2fD7Oysuw8X8rDxpUvfz0r33XZLqEgzE4sccOphCT0R1tVOfeRVWT+rY8zn1qjkpEtqXciftm3bSteu3UT\/9d\/9cY2KdD3M+YJi65JV0m3B16VXr14y4NDebuaX9s75xS+ly9tfkh7yqtNhi2xZ9am03tqFaCP5ZP1nZZcOH8nqv38in13\/jLSbcZp0XTfP1Uvdu+++u+wW61vPnj0TOvQmvUoWj1yN2SvXb3PHBccGbwc2LNlFPv53R2nTe6xcPLOjjH5yqFvv3mvYqSWLRUuOU84rXsle9MoAaX\/Yc8KaX\/\/k4Zzq\/NfvNTpfsmkvV+OVvO1eZTk+xexvS8es3Y5dZc2MNvKnb34ond67RTrNG+3uV3p8cQ\/Zsdcz0nXQ7cLx918T+pXlGHKO7f65Q2SXg78rLM+g34pB9ZY50uafP5adPzlKHj91i8yf2krALZ1xHzNmjPADepVx48ZptRbGEDDCGwOhpf+5CHTq1ElY35tuXZ\/pVSf88Gv9ofs2KrJr\/\/aJdJfARyk6\/2eM+32tn7Tfu3GZ6up6AklYXV3tylfHdLPv3CC81keYgWUN6\/03Pehu5s4o8Kc6VvaJm74orKF9\/uX7GuUOOe0Rd2FiiQMzXYTs1dtz784Juw3LPhaIqwo\/XEtkxiPVravlsFMednVtWtsjro31K0aeDz17t8TODCO+\/4kj96zt7X\/Y5oQdBIQEbexccq+b+SXcvGqn0L87hnxBPlrTcFgPv2iB7HNsd1cXPlOWi8uG2Ixw9ZY3G23R9vdv3yPVm95zfWjXrp20adPGxWnLpG3JY9GSMds897VE\/1h3x3GiMu1Hv5LuX39COvzHldK6+5CEXbkfE1yL2u4xIEZ4L3c\/tuPGrJi0WfL\/ZP1Pv+zeuvDm5fmr1rjzLhNMWjJembRjtrk7d1syZtyPfj9qpTxfu9bdg7ofcZfww8o2fb8amwT5rLC+nHWv3FMqacxYW0+\/EXZB4gGTSaeHpz8pPFzOvmuD3D9scVrXnX79+smQIUMSMnjwYD1lLYwh0MAMYolK+1\/o\/i7+6muO6LEjwp\/u\/oojaPxo7YAX+rk4M6682lG\/lvXY7uwhmJwU\/IgLUseOBN98MDbb9IuRcs49ryd2Y7h7zAVy12lXu7ruGXtJ4hWR1hcW8oqW+s\/\/4RihjTCbrjEyGqbnRwToXfkfHe9+eEZa5b+vvcf5r9uEoe\/St7W0ateHqMx45AoX8gcSyqtjvmrD62O2oEFf3eNzsYtib6IJAQfWOy1\/o5vT3X\/WCJlywTdcnD+\/Pn6ma\/ecS4+WPtXXuB0fOv5XT2H3CHwdOXyQ\/M+5NwofqlACTblX510v9315G1Fp9fF7AiH21wW7DPtTNghwjG168QHhWJo8oqsL6Rw\/foH0ImNf39c9hOmxTn4lCucn1wfwuPCpttLt9InS+dT681hi\/3gVyxrMWNT+GwIJBDbEJkFYEsM1nXsbac3keOl85DS3hzTHFT\/20rxKDDnHeADgR6V7PrLTPQQoDty3uDdunXe72+GB803zLEwfASO86WMlzz\/\/vPuxGqFfjHWepIcOHUqQsfgzuUPHbndk7cLbTpSXD2xYj+uTzg2x2cl0GqkZ2jHUDNKpGT7hQxc2K4v+32+OcUSaJQ6QbsLff7+v2\/+W\/K3VNe6HdJAHBJ3KKVd+T\/Y+8FnZI\/5RjcXLf+7qIn\/g8HnCUy03U90HmC3H2AGCmenpd3wxRnj7yGEnT3O4+OSZ8g9+9xdNiL32oVXHapGq99wODfjLhzb48RprlyX+j3Q86oINsQs0ke0L7pblU\/\/onqy5MKuePJNoI8BYQnT5Atmmv3y7UWf02IX0ItyEGhmUbyKtnoFH9\/\/4vHT7n+uk++hljbZf0l+jgy8\/Qlry8ua06jSj8kNAjwGunewPz2QGpI2ecgzxtuDsx\/Z013Z0Jk0R4CEA4UGAWXAs2N0JbFlLD7ZGfEElfTHCmz5WidcE119\/vbBjA0X58toVV1whx\/1\/9s4Ezori3NsvO7IoIkJUQAEj+YICGheEuKNyFTFICKAguKAgKPoBLoniQlQS8SIqgoILCKIXlQBXxQQXNCAoQYKSZFxAHFBBQFAW2eSep2aqp+fMdmaYc+Ys\/\/lNnaquqq6ufqq76u2336o+91zDxpe4otxRuas0oIXloi0sH52Bjz9wW2V7e\/yxtnbiLvNfUCONfdFyIng+M\/DEyJNgjkkDabE6TAWwIcZWlnImjZjsZmBXydW+Ug6vmBjE8NEEExfteCJF6MT5NH8OH334ohNsKb92zbV2Ya9nLUeIyKvvtu0nBJ0eZfkywj4aYcwg9qz70bCn9GkHt2xm0QIwaQj0TIZDwGVCXFioJZ44BGoEb85vyNxngwl67M8yMWgkyMs2vBnACculLgE0TesXTjeWPHqkxwnuRKq3yHso5LrlWncJ+omJQO+ZLZy5A3MM\/A4MxGjzpnVdKW2Uh5JBPvcZ9xjXgD9t+nbWpkbQRYhDwcH95tPlF04ARvRJsGO5w3AuxiWE4d0bFoWjFS6GgATeYuBEJ\/ElNYzBTzrpJCf88nnhiy++2AYOHJhvIlv0fmz\/VKuBHbi1CkHnuJBdIPITjl8wrpoTEjFv6PG3g+w\/S++1Oe+\/ZnP\/mLdvZJfgH+Gu26g5Tvs5YPyd7vUrQhzuyPa13SdJ6WDY9g6hc1\/NZwwb4rpn9g3KInDWFQ+7soYu+ST4OAbxhbmw1rmw9IOaVg+it6w+xPxrYY7vE5j5jVYNbRCCbbRgTd4qtXJenaKZPeWYO9yu8Gtw\/Oku7H9OaXmH8Znhiyeeag3b9\/LR+fzCBOSqDU7JJzh\/MS\/L7dPr9iucz4+epKGQeo5BAdOFLa\/+0aac94lN6toqeCPw2qp\/Oy0lbxi4N7hPGEBS7ywrtsbciwi+4Vr4MPyPanydaVD2RNLTp51RCvDKffO804KP\/nC2XB9fbpxtaCkRdImTKz0BHhroo3ho8Htj5vDDwl7OHMs9aD6wzifJL4SABN5CoJQUxSoNK1euDJYmGzBgQEm7FEjvHXmdU2VKXfdhieX93nNCJra9RQ24CJcMyjhsCwsUWEgEHQ0CI350MlpXXFiQQ3j2+ei8EERxhMNpPo\/3o8svLO+qdV3suT\/38Ls4P1q4ZUJQYQIpNzraVh4C3lz+pB175Hi7cvI3rgyWTeJVmddWE0n+sLDdsv2b1rv7KU6LW\/eQvA6B10LkpyPuNvNGN7GNcjCzcPGbV1rtml85xznSqROPI4wjLJecBLi2GQTW3Pxz27l2nFnlnGuG2tKenR\/OfRD21LQAABAASURBVJDKfaAiXm7\/CUQPyns2LrLwoEy77P9RVEIyEKAPZHzgPkOjuyNrrKtW3TP7Rd48Vg8mohU1rrnM+omZAP0WDw3cY8gCTHDzO7u2GL3OmMPi4+TnJ1BBAm\/+SmTqVtP6BxQ49bCgFrbt5UL3mRFiEfJw3f6YbX4C3ISBdzuzBJ8v7NMhIbziuDHC5YXzLZ03zGmY33\/htCCaOiGccjxWOrj+haVOk9z54SYuT90GeUIkWtaeT+Roo2sfsNSlR\/+g7WHRbW8\/TF04J\/IhwGK\/3G\/wGa7DbNqhNtHW+PjvnB\/9g3BL3DvrprkVKao2PIrNfI4yMYe4ekZV8\/a9PgM8cHTaPm73hsXO7hfTh20\/Hm49buxjXTsPdPUhD9poeOIIEyeXPAT2bl9ja2cPMf9Aw6ojVRvVDCrI4MsDp7\/mggQFyoUA97MflNGeRw\/KtAv3TrkcTIVUGAFMF2hHBF1fCZQoddtPt0MHPW2DlvyixLeEfj\/5pSPAPcYeTHBbkZXfpKHzoL6Glp2Pxqx6Ps98kPyZ7ipnOoBEnv\/GAf\/rNLnYjPoLluMf\/FrOOrPY9i5\/fhNRhToEM+8KzRCJZDAhDx0RwhhhHEJuJNn9h49NBMLn1XfdkmNfGzJDIM27sFnBnvV56wCTfs51LxtCJVpUbIK96QJCKpPNyIND0MBHk0sdEKDZDmuEsVOq0XiQcSPTYXqh5Oc9+pLVueyfhtinH7dx4R+tixPQvZ0un6AkgWMNX7HdmDC3ZXV9oqxaRBj+f7lfmbvy1nPymTEgWFNvJsr9bfQp+T6igWB++F1vuTL4mX93niCPnSJxchVLgAeWF3utdatr7Hi\/j32XtTKoENfarlrX2eVzTw8e1IgLMlR0IE2PD2MeSLmXF857J99ZntFjuP3493Nt5dRVtvb9H\/OlaSN5CWz7+if7aNJOtxwdQpXvw6kx7c3bR\/p3tuUSQ8A\/XGLqQBtwVB5GXuh\/qK2aXs16bH\/CNi3IURyRlsmuciaffEWe+9Sun9uqzuuNL6kdu7JGUBUuXuxQebU+v+98JyCj8SQ+yJQbYJmt3OB+e9UOLagZ9YUikDJwMTENs4Lx13d1NkP+1STCcNbCc4wPSyAw+nhMBZhsdt20QU7QCJ+DX46G1zK9Z7bwhzKOU\/f0W90s8CAyEqjfsrl7PYYg23tm80BzcEjbPE10JFughWVWMK9RZ4++11hHlXq998wBwZJoTFYKrwOc9dS9NvWlWRRhWQvPdj4\/dQ9Zb4ecvN2dLw8TxHUZ+ogzjyCM8+dLWC7xBBB20TJhD\/5W\/1VO6G18\/GZXEQYArpefdRzuri0XqZ+EE+BtEA+wDMocHNMlHsIXP9rAfUWS\/lD3EWSS13Gf8VW0jyftcpVEEUAfz4RFBF3a1ysoXAb9JIwA\/RzjK21Qv8sq27j9ruDYdfc1tFanNgu2kz0Qz\/pJ4I0n3aiyD5nQ2Wkjx3bqY1uKWVqsTiczJpS5ZbWs8D9eHSF8vvF4nrBYeM7CY7kxEDa5URBkC1u3F9tXXkkW1YkhCFM6y4mFhUcfTxqvtxBuCYcdgi2dZTiupDDCLgMndfZ50SbQ2dLpMpgi3JDGRDTOizCOfegQEM7Zhh8dNmHSTu75d4IFHHnIS8KqeVluTd89GxYZAzZxOCYH8kSNkL0z+0Wi5BJAAOYfPDDREHb94fbUbmO1zvu3cc1xfXM90L4+XX7FEaAduAdpF97khO9PHlZ4oNy5uvD7sOJqndlH3pK7TCMUaDt871iFh9UDes9s4d4O+nj5FU+gYYdexnhITX46YJl9XfljghnvJPAm8BI4aFseboQkf+j1JyxwE6b6DRySr+Oo8s7xTkBGsOXp2ueP9jvdsdftf3brq4IkhE4EVV45hR2DTZCpkEDL9m84rTKv\/I+vepNtfPqOINfBLQt\/Sqx6SDtDMCQjgxrHJYxDIEW4JRwvx\/HodOmQOT7HQZhmmzBxnXMnKNFBT7\/nNaKdw+aYCXDVIueANtdFRn669+\/lmPLZ4sim++cc0RjTHmiBWckCnpTPOSL4I\/R+OfgA27HibbePfsqPAIPvuw+sM7SBmOx8MHqitWh5X3CAn198gJ3zWC23TXsQoG3w5ZKLAO1z5d+PsmYX5tkY8hC57Z997PsFvZz9NSZZyVXrzKmNv9e4zzCJC5957cMqW48ZTd0bt3C8wslDgH6P8e9Xj26x\/6n0WPJUrIJrkieBVXBFMu3w549Y6wTLC+7vayua59mw8VoPjSKu\/Uc5g7dnM2Tus26fs\/q+4KOcj+0pgZ81XubSMSFAACTOu83PfmG4H95+xkWhTaEjo2NzEZEfbhKWBosEbcbE6W4SGOvV+jysFEEajrwImoQRatEYo4HtPbM5UUnh2l5xXLAsm6+rRf7CDxtN29exA1qdaWiEsWX2Dwe1a3xldnDjSG5zZhKcLwJxWJP9TXZEFe9y5Py8dOtF7gGlynnX2Jiz6xc5gTAnt35LQwBtLitmoM1FG8g1yWTNXZ9vdeYytM9FuZMoS1Ou8lYcAfqTdnfUNARftFE9Jq53lWFVB\/rBaVrH1\/FI9A\/3GoIu91r0sVEAdJlZ2444uWZ0krYrgkAJx0TxtbVSzn1VQtaMSJbAmwTNfOSR25wJQ632DSLCVZ7GI1w1Vkrw23u+\/cIHnb8z8hodE4hDb7jRflx6lG16\/O9Ow4hQgOCMw4YVN2VYO7dPUT\/bd5zghDa0mT4PwgThzhGBwgmEj1\/rZuAS5x15SMf3ccnghwVdX5\/jeh7sg5a1fZQLoxHG9qnK8VnGq1ZWaHh+1iu2buked64D36ljXYb\/IdBksxMD8n2Nljv7XrRRCGLET+6fo2n\/coG+NAWP\/XVc398v7GXfv3eW07qHy9v1aR23yQOXC+gn5Qgg+KKNYnJbtKkDfRgP0ghhxb3lSrmTTrIKwxklyHf\/c7d90Tdn9R1fxWTr03295ItAaQlULu0Oyl92At93\/m83YPOVs9fvOcIJlq\/eNtmaP9jFhR\/rc4d9uWBbcABe8XkBk07\/mYgghfD67qu3GNpINLkIme8\/f5rbH9tgVilAsOWzv7F2VOcMeslpyThw9D7R29UOvNO2\/vUbsqasQwhGQML2F37+RNAqodnwwj6v7uo3XWK8YmWZK\/KxRBk+XHCEcV\/NHoKXz1EeEQwmDNaUzaBCXIJcyh6GhzreRuQwW2MIPJwMbzMwI6H9aLt9vT8iWi5NCGDqcPnrpxsaX39KLK7PigBoHLke\/H3l0+WXnQB9E32S5\/rdjJzJTte\/9LTRv9FHcq+V\/QjaUwSSh0Dl5KlK+tekZpt6TpNbs81BRZ5s0w61nfBZZUpde6HjliLzkeDNDwhHO7Qm0XFs84ojx6+O51z4Nb2LKOan7pn9rOljq4rJkRpJdOYIvuHa8lDhtxF2ed3K9oo5tezZIXc6ze\/707u7h5bLH3rcftZkLsnOLZndy\/nRPwwoHKtpvRFGmMGaASY6n7bNCbUItnwV7cvrmtnGyf0ds6yFZ7sHOh72Dr\/7LYcKprRXdBu6RP2kNAHaFo0vghYPN6zjS9tzUtxD3D8IaGzLlZ0A91rVDX\/IZ3qFsuSwu96yur8e495spdf9VXZW2jM9CEjgTWA7Vr\/3Sjdwjzn97CKPGhZUj6i50eUrbPkxXruzbi+aw7C5g9shhp+Te75rrFFLVj4yMbXr5wSdY8BxgcgPE7oiXsb8o9FAkGJSRsMTqhj2yWEe67ecaGgZqx\/ZwQ5vmaNdZNk4NO6sQ8zEQ\/YHGDa\/+7LOceYlVnUZUc4h9BJgwMHPZIcAw\/nDYsOyd21Kp3ds3BVd7MCuje0fmx8gKXC0A+yDCAXSmgDtzQli6hDW+BLHOr5ofXdkPWT0gf46Ik2ucAIwou+hr+drXPDDXIh+Ctb0W\/R\/zGkovATFikBqE5DAm2Ttt3rhNmcTuvfyHwJThyeGPumewr12NrrKdFRoQoav2O4mrfE1NFZuIB9mD84sYuLFRny3P2YTXazrN\/FJV06\/QT3chK5iM6dRIqeCRoNX5eEHjxZdziHJOb68hpYbm1+WbYMty8fxULJn3Y92wC\/PdLOXaY8B79RxGkrMHZgEx8DiCon8bMpa6b6GM++Gca690VitXrg1kpIZ\/wy+CCqcN2fMq+zVC7YGH\/xg7eQTu0y3k3u8615vMxCj8SOvXOYR8BpfL\/hi7sVD0vqF092ydFxHXE9cV5lHp+Qzpm+BEdpx5hqsWVbPuOfQnnNfMdmYfo\/+r+TSlEMEUpOABN4EtdtBoSXJeJq+YFQlJ1SiHdz1h6dcGLup8AoC4aoh7HYbNcflC3+9LJynsHC02QNfGyNf1sJzjCW2CEc7lu7CNRyYt3xXdJ5M2ka4ReCi3XCVT73NnT4DBoGsyCt3lirjdWDVn9Vwa\/Vi98tsc8xFXnjoWSOdvEOXfOJMVmpsH8+mndBxtPMZqKd1XenC6f6ze8Mi+8e415ygEj7Xk4b3DzZZYaNJj5es48ODDGFHA3GAJmMD3HtcCwhoTDBFWONNlweCje9P791v2H+jyfTxmerTpyDgosU9aPt\/5cMAyy83znZftKQfYztfBjNtikDaEZDAm6Am\/alWg3xHopNhmbDvZ67JF1\/sxp62bmmxX1R\/3tmR8kU2BgA6d1YLeKBVLXvr6SGuCP\/qt9Lu8+3D0a2dKcWk21+3uSOruHQEaBfI\/YkWKBDyfBm5WTLagw8DLc4PDqzdy4DrmQNo\/ouX4RnmJpiuhB8qEOLgSoa923PanQmADQ9aQpRzaGLQXLGeL4Khi0zxHzRvaJdwnN\/nc96wRU83C86KgTnYiATQ4nFdc49ENvUvAvkI+PsPUwfuRz464zMw6Qr7bwQ9rjeuvejry+dNV5\/+Y2rXz43zZ2zIiig3iDvrirFuIhr3F9x4M5iuDHReIlAYAQm8hVGJQ9yW2j\/ZoSP\/7DS0V8\/M\/9WTHQ\/1cwLpI92usM2hL9tEV4OVF9Ai8gnc9ZtPDJLDwmtWRNuIYEtHT4ZZY0+zxVkjCTrnB4vwK3sSEDDwY3bKaGjBW4eWOIMtHBHUwuYLHhUTEnkwwSEQIwxPnbHYWF\/W51n29EfO1GHBuKq2bdlwW9m9Uo4NsM+QQj6aJWwsP53+ujPtQPCYFtFi19y0PTgLzELg5iN4rQpDvy1fBIojwLXTe2YLN8Hqmldr2sHd77IP94wx\/rje0Poi+LGdro7zRLDlS5Dcbyzhh9mCP1+Ef754eeqoJx0n3V+ejPxMIyCBt4JaHIEU4RUXNndAEB388v1OMB7w+LXOv+mdN90X2Ejz1d33q9bGur1+u7Q+A4XfJ6wh8XHyYyNwxohfOZtdBDW0JnCt0aSbYWMYLoF4tMThuK+yjnWbq9Z1cX7P09vYqSde59b6XTKnlz122Tir3qKOfTvuCpce\/cNAFx2XLNsMvGipd2SNtd\/cdWfAA+H+p0X3uTcUMDlrVM65+3pL6+QRGuHxAAAQAElEQVRJyC8NAa6lBr86xg7+3Z3WsH0vY9vvjz29D\/OwiQBM\/+vjysOviDK4\/xFmOR\/86ZftNO43tLnUBwbcT5hjMfmWODkRyGQCEngT1PoItd\/ecbPT5E7qmiPoFHZoOqnC4qPjeI3uJ7NFp7HtV24IC8nEhx2CNJPaYpnIFt5P4fwEGFRwPpZXrX4FB4RYTE+ufS3ny0SFPVwgHDP5ze1f+RvLWpC3isfcl8caGv1Z\/d9zyf7HD9yYCPi4RPgMshwT549HHIOu32bA5RpnPWO+hEY8r1M5z0tvfsFqtjrTaeJ4QKjfMnm+zEc95VKfABpMri1c+L7kOuXs8NH87sh6yDAbShXhl3qjyUW43Zn9olX6sq+bzMw5eYfJFNpchHzOnwfx6Adtn1e+CGQaAQm8CWnx4g\/ycfOdtmZAZeNJPJxz2u9HOgEZW1A6u3CaD2POQIeGwMTrYR\/vfT4qQRru6nv6uglAPg3BjC8baRkaT6T8\/FNG3G03vbnJsIM+dNDTzmeQatHyvuAg3na3Zfs3XNxhfxxjVRvVNLS7LiLywzJoEc9WzK5tq58+xdCasu0dgzXlZj\/fzS3PhCA8tevnFn29hAVUv29pfIRYVpaYc0O2Teu60jkv5C5\/\/jt3PI5NvlHNNjktNfbNTNirtKuT07p1H\/uhHXLFSDv8rrecJq40x1deESgtAZQHCHx+P65dH8bfEXn78MPCXm7yJNcu1zP3E2nJ4Pw9zD2FYE79EHYReukHmBTr68m5Ypu7YcedbilFHy9fBEQgj4AE3jwWCQ2hhUAIxc09datlN9ztzBaoxKOX3OYEXTRkbHt37vDFzsQBTZmP8z52o2jRKA\/B1ms2iA\/n8WH58SUAdx4k+EgHPkdD+9K47WZjZY5+g8+ws9tcZde\/sNRa96hvS2b3NB5s\/AL75McOONzWS2b3sqoNTiEpcNjqrV8w3WrXWuoGbhKImxoRShkov1\/Qy+bf8w8noDKoRw\/ofM2Mz4n6eIRnlnryH38gjXK2fjjcljww0SibY+AYkP969W3BcYmbedmLeIFjGbfsvfe52eDYPMMlSFQgMwlU0Fl7W18EQ960sMpDuCpcz2h+WZ92R9ZDhnDJfYEjLZw3XmGEWY47NfLQioDL\/Ut9EMzDx6Qv4X4atOQXTlGCz5iC8iOcT2EREIE8AhJ481jENXRgvmXJqpXpWNFLjBVVCFpFn\/b+86c54RlB6olh99r+avp8ufLLRuCILmOtbofpZnvbutf6rOlbvUk3q9tgfVAg2hpeSV457XXz2l8SeQB68MRj3Lq9YW3+TzvWGGsAYz6BjSx5GaARLtECHdf2t0Q5x4C+\/IXvXBhhF\/tgZrYTzyD7o3WxV2\/d5z7+sODvVxlp68f\/l333yUo74eypzvbW7Rz54eHqxC7PR0I5\/wjzPHT1GXt3RIg\/2Dk0bP7hKyeXfkWg4ghwbyEYIvzyhmvNjqWG8BuuEQ942\/7xoDMX4L7A8SCHEIzjPkEoRRDG+X2558IuOj6676UMysL5e5I4wv7BkvuXcriX8V39hzUyNLkIvWxLyIWMnAiUTCAZBd6Sa52CObIb7rHRl22090bsNTpbOjWEUNywaYfYqfdUcYJM9KkhxPQfMcZNwmApJ\/Iz2Ym1eAeMvzPQCjPBiTTnBt4dCLbelteXW5xNr88jP34EGLhqHPnrfK\/1mVDCDGo\/8Pq1mNHmEk\/797ixT\/BRBmqHjfCVt55jfEyk7iHrzdvKntp3h11+RQcnmO5Y8bZhssIxuV7YD\/fp85PtvVuvsh\/ensymMXGSAIMsAzjXF9sb9rQzPqox48kxhnnCujmVjYcpTDWGzH3W+GPQRThn4L14UjtjfdQmPV8KJvJJ2IWSXLIS4PqkP+YaRkuKX7vNs1al6tX5qlw3co\/tyBrrJoVxn9B\/Iwjj0MCiiUUjG3aYHeB83LTIW5dPb73MWDZt0+xmTqCmLByCLgfk+PjeTb\/nNXdPcR+TRh0R2CXkekLyRSB2AhJ4Y2cV95wIDRwEYQLtGY7tolyVWo2LSgri6dD9BuXj\/Lb85CLAQEaNVi\/chmcIkwRqtjrDvn2zOkHnGrfdHNH8vulWcMDml4ecKU8vcNtHn\/j\/nd98xj7bdWB7y15Wz3h1iyYWoRdhmC+YZS08x54ZP9YJu2FtMQO5O0jkB3OKOYtfc\/a4kU2b9eYcwzwDt3vjIldu7eMfICnyENfcPZS5Df2UIwEVlSgCvm\/kgZTVHpxwGdGmYgLRZWI7d71zr4Xrwz5ohL0mNpyGeRAuHDfj6VvchzHCcdFhjsebES\/gRqdrWwREoGwEJPCWjVup92qyvqp5Te7Urp8Xu\/9PW1q5D0ygtavX5yir1a5VgfxMbBrXt7+bLEQiHS++d2FNLh03HSi+T5effATQ2tBGuHDtNlX\/rb25\/EkXRTujSa16SDv3kZHwADxr3gdOC\/zuq7e4twVolqZFtEqYQTgtb5fn3VJ2fP3Na3HRXGGGQOEcl0GWMMc5d9hi82k+7j\/bRxF0mmNeCbuNyA\/5I57+RSBtCHBN8xCK4y0M1zv3CJOLEUpxbKN9ZWUE8oddjSMGGS4cx1sc7Pp5E0KfjKM8f98Bj+OFFRXEyYmACJjZfkKQwLufAOOx+5Rh7Zzm7alRbzj727EXnOcOE22e4CJzf\/pNfNJNaBs0eaLRedLJ5iY5zZs6UE8juf1wu\/maekGYgZEBljx1IprVDoN2+yzOp41bXvEHQ4PrInJ\/yP\/thqtzt8xmj743CDO4IgxfN22Q1bTZLp7rh\/gDO0y3n\/fs6xarZ2Du\/HBjI95l0o8IZCAB7iXuR+4DHAi4fxCIuTfDrvavhhkuHIf5BGZB7Mf9iqM8tuVEQATiS6ByfItX6UURoKNDyMBsgVn7zp+7uajsLv7Yi7Y5odZt5P7QAecGnbd31\/+6tSXdhn7ShgDtHB4YGWR\/1nG4E0bRNHV+uIkTRhl4m\/R4ydnQsg\/XWeeIoIrgihZqS6257gGItM6RfdpecZzT1qKl8iYUQGM\/fOLIy3b4+KQloVOVREAEREAERKBQAhJ4C8WS2MhP3jTb\/OwX9uPH\/4npwAjIN398oNPk+h1m3jPAaYN5xf3vN1r6aPlpTgBhFE0TAqk\/VYRhttEsIdR6QRVhuFnHlk5IJo08fh\/5IiACIiAC6URA5xJNQAJvNJE4bbNKQ+dHrnIa2kun1QiOwlJTb895xZkwYMoQJIQCCDVsfjynthNqn7jpFZs7sgpRgQvb7DLTPkhQQAREQAREQAREQAQynIAE3gRdANnf\/Wjfz1zjNLk\/bdvjlg1jKTGWe\/JV4KtphK+e+bETjDFzwF357DKizaezsTl7F17gmnaoE4TD+YJIBUQgioA2RUAEREAERCBTCEjgTVBLt1pZw95583WnyY3WzpalCqzdyCx8vy+vp7HlxGdihI+XLwIiIAIiIAIiUCwBJWYAAQm8SdDIrHeKJrf7Q1+72jBRqMrxWc58AZvcp\/q0dfFhswUXEfWDLSc2m1HR2hQBERABERABERCBjCYggbeCmp+JRINfvt+ZLnQZ9od8teCrO2HtrTdRwJb3NyPGB3l9fBChQPwIqGQREAEREAEREIGUJSCBN8Wajhn4KVZlVVcEREAERCCNCOhURCAVCUjgTVCrrWi+02lzMV3A7GBL9i579JLbArMFTBdm3HhYgdqQN2yT+\/KtFwV5wvFBpAIiIAIiIAIiIAIiIAL5CEjgzYejYje8je75I9bmCMePX2tMQgvXCrMGv43Q7MPJ5as2IiACIiACIiACIpA8BCTwJqgtGr9dKdDmsk5uWHAtqgp7vv2iQNKR7WsHcasXbgvCCoiACIiACCQhAVVJBEQgKQhI4E2KZshfiUVPN3PC8aTbX7f7Gi03JrH5HKzEgLCM5hfn4+WLgAiIgAiIgAiIgAgUTkACb+FcEhI78JGZ1vP0NjnmCxMvtk537C3xuAi7vWc2N2x7S8ysDCIgAiIgAiIgAiIgAiaBtwIvgpptDrJ6fY4KalD10JzwQU2rB3EEEHLxvYve9vHyRUAERCC1Caj2IiACIhAfAhJ448O1QKmrT852mtzrpg0yzBJ8hiWzezrzBVZpwLbXx8sXAREQAREQAREQAREoHwIpJ\/CWz2knRylPXdbJlszpVaAy2OayfJlzj1+bT0AukFkRIiACIiACIiACIiACxRKQwFssnvJLrPdSI6fJfeyycbb8he9KLPib7E72Vdax9tnkErMqgwiIQGYS0FmLgAiIgAjESEACb4yg4p1t0JJf5NPkZi0822aPvtfmLXzcVi\/cGu\/Dq3wREAEREAEREAERSFECJVdbAm\/JjOKWAyH31lX1nW1vpS\/75jvOv+e1DLb9BymCCAVEQAREQAREQAREQARiJiCBN2ZU+5dxc\/auoICDmlQLwnu3rwnC4YBfduy0YY1MqzKEySgsAmUjoL1EQAREQAQyl4AE3gS1feNbP3Wa3MEv329Htq\/jjjq16+f2QKtazrb30Utuy2e6wMQ1NMDhFR3cTvoRAREQAREQAREQgbITyMg9JfAmcbNLs5vEjaOqiYAIiIAIiIAIpAwBCbwJaqo1o35eQJMr29wEwddhSk9Ae4iACIiACIhAGhGQwJtGjalTEQEREAEREAERKF8CKi09CEjgrcB2vGBUJWfXe82YC63\/iDGBbW8FVkmHFgEREAEREAEREIG0IyCBN0FNuiV7d3Akb8pQuVZjF1e5TlWr2aae6S9VCajeIiACIiACIiACyUxAAm+CWqfFgLnW8\/Q2TpMbnow2of8sZ9s75cZrE1QTHUYEREAEREAE4kRAxYpAkhKQwJughqneoo7V63OU1WrfIDji8he+C8IKiIAIiIAIiIAIiIAIxIeABN74cC1Q6r+H\/9ppcsd26mNbcj9CwVq7TdvXdh+WIFxgp\/SM0FmJgAiIgAiIgAiIQEIJSOBNEG4v5HK4sElD75ktrPfM5qYPTEBGTgREQAQyiYDOVQREIFEEJPAminQxxwkLwMVkU5IIiIAIiIAIiIAIiEAZCEjgLQO00u5SlHY3lnKURwREQAREQAREQAREYP8ISOAtA7958+ZZ8+bNAzdhwoRiS0GD2\/W5J9yau1fP\/LjYvEpMDgLffPONTZkyxfCTo0aqRUkEaCu1WUmUkie9DO2VPJXP0JqozTK04dPktCXwlrIhEXaHDh1qs2bNspUrVzp\/\/PjxVpLQW8rDKHsFE6Bjnzx5sgTeCm6H0hxebVYaWhWfV+1V8W1Q2hqozUpLTPmTiUB6CbxxJrt+\/XobOXKknXvuuXbccce5o+EPHDjQpk+fbqS7yKgfTBpmXnqNW6VhUtdjo1K1KQIiIAIiIAIiIAIiEE8CEnhLQXfdunW2efNm69SpU769OnTo4OKXL1+eL95vbM5dhoztg5pUw5MTAREQgbgSUOEiIAIiIAJ5BCTw5rEoMYTAS6ZGjRrhFXCfREd3JQAADaJJREFUffZZgbjoiLVr1tqiRYvkkpzB2rVrXdMtXrxYbZXkbeXvJ7VZavUraq\/Uai\/uM7VZarWZby83mGX2jzt7CbwOw\/79IADXq1evyEK2ZO+2Sr+o4tI\/+HquXXrppXJJzmD48OGuvcaOHau2SvK28veT2iy1+hW1V2q1F\/eZ2iy12oz2OuWUU6xx48ZuPMv0Hwm8CbgC+IrabfNb2eX\/qG\/DXuxizz33nJwY6BpIpmtAddH1qGtA10AaXgMPPPCABN5cOU8Cby6I\/fEwdcC2t6QyeMpq166dyYmBrgFdA7oGdA3oGtA1EO9rALmjJNkkOj1dtyXwlqJlMV0gOwIufrQ7+uijo6O0LQIiIAIiIAIiIAIiUMEEJPCWogEQeLHVjZ6ctmDBAiO+devWpShNWUUgVQmo3iIgAiIgAiKQWgQk8JaivRo2bGi9evWyP\/\/5z8YHKNj1o48+Mj48QTzpxMmJgAiIgAiIgAhkAAGdYsoQkMBbyqYaMGCAPfHEE3bNNddY8+bN7eKLLzY+PEF8KYtSdhEQAREQAREQAREQgQQQkMBbBsgdO3Z0nxXm08I4CbtlgJg5u+hMRUAEREAEREAEKpiABN4KbgAdXgREQAREQAQyg4DOUgQqjoAE3opjryOLgAiIgAiIgAiIgAgkgIAE3gRA1iFiJ6CcIiACIiACIiACIlDeBCTwljfRqPK2bt3qPk3LBDfcGWecYevXr4\/Kpc1kJTBs2DDDJWv9VC9zK6Zwb3mneyz5rwpWt2nTpo2b+Eu78dla+srkr3lCa5iUB2P84h6bMGFCUtZPlcohwLjFvRXt\/ApTObky61cCbxzbmw6c1RwOP\/zwYJLbSSedZN27d5fQG0fu5VU0HfrLL79cXsWpnDgQoI24x1g5hQmkOO6xc8891xCq4nBIFbmfBGiX3r17u9VtaK\/ly5e7Ei+88EL1i45Ecv+wLGd2dnZyVzLDa4fs8dVXX9kll1wSyB7cazgm3WcqHgm8cWz5qVOn2ooVK6xv377BUW6++WYXLhdBypWkn\/ImQGeBxomOvbzLVnnlR4B2euedd1ynHu7Eucf4EAwfhCm\/o6mk8iIwefJka9WqlSH0UmadOnXstttus82bN5sXfomXSz4CaAc1diVfu0TXaPv27bZ27Vrr1KlTdFJGb0vgjWPz80U2OvZmzZoFR+HjFGigGKgZsIMEBZKCAG2CxpDO4o033rB27dolRb1UiYIEEJSee+45Gz16dMHESAz3X8TTf5IRoL1oN9qvPKumsuJLAFOGkSNHGq5JkybxPZhK3y8C69ats3379lmjRo32q5x021kCb5xaFMGJVwqYMxTWsSNQ8RQWp8Or2DISoK0YjOfPn2+HHnpoGUvRbhVJgM4ebeHRRx9dkdXQsWMkgCA1ePBgp\/XVA2aM0BKcjfHsxhtvNJQ1mAsl+PA6XCkJ8HZry5YtNmjQoMBOHpt5zIlKWVRaZc8ggTd52k0DcfK0hWqSXgQYmO+\/\/37DpAH7tfQ6u\/Q6G9oK0yGEXB5QMGvggTO9zjI9zgbzPJQ0mAulxxml91n4t1vjxo0LbHj5IixfhsUsJb3Pvuizk8BbNBuliIAIpBiBu+66yxYtWmR33HGHYT6UYtVPXHWT4EgIt7xNYSLNgw8+aAzGTEJMgqqpCiECaAXHjx+veyrEJNmDmA3985\/\/tOOOOy6oKjbzPFw+9dRTxsNmkJBBAQm8FdDY\/umrAg6tQ4pA2hJgGR4m1LBiQ3gSW9qecBqdGO2FRn769OlaqSGJ2hXBiDcmmDHQRklUNVWllAR4wMTEEk19pppTFiXwlhKlskcT8BcXdrx0GtHpRxxxhNWqVSs6WtsiIAKlJMD9xatxCbulBKfsIlACgVWrVrmVhri3\/HquaAlZloxVbFiPFxvsEopRsggkBQEJvHFsBmx1o5+m6Bw++OADO\/300w2hOI6HV9EikPYEEHZZVYPl\/2bNmmXx0UKlPcaEnSDtxcMJ2vjwQYlHOSBFQJhKxYd5Jc6rccxOvMNkiFUasOdlcq9Mhyq+ncI14F7iHsMR9mmEuceYeJipbSaB118NcfB5RUexzG7lYiPMUzG+TyMsJwIiUDYC2Owi7DKphsG5bKVor0QR4CH\/yiuvNDSG4ckztCOCFGnkSVR9dBwRSDcC3D\/cR9xP9Iv+\/LjH6CvD3wXwaaXyUzizBN44Nh5PUa+88oo7QuvWrd3yIGh3Z8yYoQk1jop+RKDsBJhM87e\/\/c1++OEHN+HJv3L1frSGo+xH0p7lSQAtPNr4oUOHuj6R9qJfZIAmrTyPpbJEIBMJcB9xjzHZkPsLxz1Gf5nJigEJvHG+G3ja8jOReSWkV0BxBl6Oxfu2Y8ZrORab7kUl7PzouKNft3KPecd9RxsmrEI6UMwEottO\/WLM6Co8I4oc2mvAgAEVXhdVoGgCuscKspHAW5CJYkRABERABERABERgPwlo92QiIIE3mVpDdREBERABERABERABESh3AhJ4yx2pChSB2AkopwiIgAiIgAiIQPwJSOCNP2MdQQREQAREQAREoHgCShWBuBKQwBtXvCpcBERABERABERABESgoglI4K3oFtDxYyegnCIgAiIgAiIgAiJQBgISeMsATbuIgAiIgAiIQEUS0LFFQARKR0ACb+l4KbcIiIAIiIAIiIAIiECKEZDAm2INFnt1lVMERKA8CQwbNiz4MhhfLirKTZgwwXBt2rQxvgZXnnUobVlbt241vjjn68o5xFLGvHnzgnNNhvOIpc7KIwIiIALFEZDAWxwdpYmACIhALoGTTjrJevbsGTgEwWrVqtl5550XxJHesmXL3D2Sxzv66KONLwZedtllMVWqVatW9sgjj9hvfvObmPInfSZVUAREIOMJSODN+EtAAERABGIh0KNHD7vvvvsCd\/7551vNmjVt0KBBQRzpZ511lvHZVT57zOc9Yyk73nkaNGjgBPPjjz8+pkMddthhduGFF9oxxxwTU35lEgEREIFkJyCBN6eF9CsCIiAC5UYg2qSB7TPOOMMWLVpk3bt3txYtWljr1q1t4sSJtnPnTpsyZYqhQSae9FWrVuWry44dO1xe8mCegP+nP\/3JMFnIlzHGjW+\/\/daGDBliCOSU1759e1c+x4mxCGUTAREQgZQiIIE3pZpLlRUBEUhVAuvWrbPBgwdb06ZNDWEVkwjMDPr16+cE3ptuusluuOEG++STTwxb2x9++MGdKkIt+z388MPWuXNnZ2qAP3XqVLvuuuvM53OZY\/ihPITdJUuW2MCBA23s2LFO2B41apQ9\/vjjZhZDIcoiAiIgAilGQAJvijWYqisCIpCaBHbt2mVXXXWVs6X97W9\/a\/fee69harB+\/Xp75pln3OQyBNE+ffrYF198YdnZ2e5E33rrLfvggw+cMHrnnXc6UwP8SZMm2bJly+ztt992+WL9QXu8YsUKu\/32253AfNFFF7k6derUyU2y+\/7772MtSvlEQAREIGUIlEngTZmzU0VFQAREIEkI1KlTxzAdqFSpkqsR9r9Vq1Y1NL1HHHGEi+OHCWN79uyxvXv32r59+wyB98ADD7S1a9faK6+8EjgE1xo1atj8+fPZLWZXr149ozzMLBCkd+\/ebUy+GzdunD355JMuLebClFEEREAEUoSABN4UaShVUwREICkJxFypSpUqWeXKBbtc4ipVyhGCLfKH8Bnx3P+2bdvsm2++sa+++spuueUWu\/766wP3+9\/\/3jZu3GhoiMnndojhp3Hjxta\/f39buXKlMREPO158hGnZ8MYAUFlEQARSkkDB3jclT0OVFgEREIH0JdCuXTtbvny5E1IRVMOOCW+1a9eO+eQrVapkl19+ub333nv26KOPGmX\/61\/\/coL0tddeW2qb4JgPrIwiIAJpTiC5T08Cb3K3j2onAiKQwQQwe6hfv76tXr3aNm3aVK4kMLG44IILnP3w0qVLjUlzixcvNpZTK9cDqTAREAERSAICEniToBFUBRHIFAI6z9IRwMaX9X4xa8DkAJteX8KHH35oJ554oo0ZM8ZHxeTzFbWTTz45n+0vZhSsvRtTAcokAiIgAilIQAJvCjaaqiwCIpA5BM4880zr0KGDW0kBk4M5c+a4Zc2wwz344IPtd7\/7Xalg\/PKXvzS0xkOHDnXlIEjzwYyRI0c684a2bduWqjxlFgERKBMB7ZRgAhJ4EwxchxMBERCB0hCoW7euPfbYY3bNNdcYWl2WLmMNXoRgljMLr\/AQS7mHH364saQZ2mHKYSLc9OnTrVu3bm6NX0wdYilHeURABEQglQhI4E2l1lJdM4uAzjapCQwYMMDZu7LKgUX9RadFb5O9YcOGzqyAj0+w7V3Hjh0LlIsQyioNLCPGhLWPPvrIfTCitMKuPwYrNbAsGeX48kaMGKElyTwg+SIgAmlHQAJv2jWpTkgEREAEREAE0ouAzkYE9peABN79Jaj9RUAERCDJCWzYsMH++te\/OpOIWKr69ddfuw9cfPLJJ7FkVx4REAERSHoCEniTvolUwdgIKJcIiEBRBD777DMbNmyYTZs2rags+eJXrFjh1uX9y1\/+ki9eGyIgAiKQqgQk8KZqy6neIiACIlACAWx\/n3vuueCDFdH2wkXtjh0xtr041uUtzE65qH0VnwQEVAUREIECBP4PAAD\/\/wc6qFgAAAAGSURBVAMAE430f61JU9gAAAAASUVORK5CYII=","height":337,"width":560}}
%---
%[output:54806142]
%   data: {"dataType":"text","outputData":{"text":"Exact static tip deflection of the distributed wind at the gust peak: 5.510 um\n","truncated":false}}
%---
%[output:91478ab6]
%   data: {"dataType":"text","outputData":{"text":"                             <strong>f_n_Hz<\/strong>    <strong>Q_per_line_load_m<\/strong>    <strong>PeakSway_um<\/strong>    <strong>SteadyPeak_um<\/strong>    <strong>SteadyPeak_over_exact_static<\/strong>\n                             <strong>______<\/strong>    <strong>_________________<\/strong>    <strong>___________<\/strong>    <strong>_____________<\/strong>    <strong>____________________________<\/strong>\n\n    <strong>#1 lumped-model      <\/strong>    17.307            2.5            22.216          14.741                   2.6753           \n    <strong>#2 Ritz phi_t        <\/strong>    17.307         0.9375             8.331          5.5277                   1.0032           \n    <strong>#2 Ritz phi_line     <\/strong>    17.522              1            8.3296           5.527                   1.0031           \n    <strong>Ritz phi_t, K_b - K_g<\/strong>    17.303         0.9375            8.3348          5.5301                   1.0037           \n\n","truncated":false}}
%---
%[output:3dae280a]
%   data: {"dataType":"text","outputData":{"text":"Case A (solid disk, clamped rim) deflection w_A(r):\n","truncated":false}}
%---
%[output:1af315bb]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{p_0 \\,{{\\left(R^2 -r^2 \\right)}}^2 }{64\\,\\textrm{D}}"}}
%---
%[output:7c3b4b46]
%   data: {"dataType":"text","outputData":{"text":"Case B (annulus, clamped hub, free rim) deflection w_B(r):\n","truncated":false}}
%---
%[output:10153e16]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\begin{array}{l}\n\\frac{R^4 \\,p_0 \\,{\\left(\\frac{r^4 }{64\\,R^4 }+\\frac{r^2 \\,\\sigma_3 }{8\\,R^2 }+\\frac{r^2 \\,{\\left(R^4 \\,\\nu +a^4 \\,\\nu +3\\,R^4 -a^4 +2\\,R^2 \\,a^2 -2\\,R^2 \\,a^2 \\,\\nu -4\\,R^2 \\,a^2 \\,\\log \\left(R\\right)+4\\,R^2 \\,a^2 \\,\\log \\left(a\\right)-4\\,R^2 \\,a^2 \\,\\nu \\,\\log \\left(a\\right)+4\\,R^2 \\,a^2 \\,\\nu \\,\\log \\left(R\\right)\\right)}}{32\\,R^4 \\,\\sigma_1 }-\\frac{a^2 \\,{\\left(12\\,R^4 \\,\\sigma_2 +2\\,R^4 \\,\\nu +a^4 \\,\\nu +6\\,R^4 +16\\,R^4 \\,{\\sigma_2 }^2 -a^4 +5\\,R^2 \\,a^2 +4\\,R^2 \\,a^2 \\,\\sigma_2 +16\\,R^4 \\,\\nu \\,{\\sigma_2 }^2 -3\\,R^2 \\,a^2 \\,\\nu +4\\,R^4 \\,\\nu \\,\\sigma_2 +4\\,R^2 \\,a^2 \\,\\nu \\,\\sigma_2 \\right)}}{64\\,R^4 \\,\\sigma_1 }+\\frac{a^2 \\,\\sigma_3 \\,{\\left(a^2 \\,\\nu -R^2 \\,\\nu +R^2 +a^2 +4\\,R^2 \\,\\sigma_2 \\,{\\left(\\nu +1\\right)}\\right)}}{16\\,R^2 \\,\\sigma_1 }\\right)}}{\\textrm{D}}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =R^2 \\,\\nu -a^2 \\,\\nu +R^2 +a^2 \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\log \\left(R\\right)-\\log \\left(a\\right)\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =\\log \\left(R\\right)-\\log \\left(r\\right)\n\\end{array}"}}
%---
%[output:5c4fcaf6]
%   data: {"dataType":"text","outputData":{"text":"Self-weight 441.5 N = 140.5 Pa over the dish; at 0 deg elevation: normal part -0.0 Pa (bends the plate), in-plane part 140.5 Pa (not modelled)\n","truncated":false}}
%---
%[output:28f3b431]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:8170fec3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:62200cae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:964446b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8ad4db5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5628c2f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:83eb0070]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3ae2b2de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:96f52379]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:997e78f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:87ec659a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:59c68b18]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0caae702]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9b5502ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:63bee956]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4f8846f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:93140652]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:21ae7169]
%   data: {"dataType":"text","outputData":{"text":"Case A modal mass:\n","truncated":false}}
%---
%[output:17bcaa4a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:60c92f0a]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{M_{\\textrm{dish}} }{5}"}}
%---
%[output:29584962]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:69412218]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:42a17e7f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2581fbca]
%   data: {"dataType":"text","outputData":{"text":"Case A modal area:\n","truncated":false}}
%---
%[output:243f62ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:71f6998f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{\\pi \\,R^2 }{3}"}}
%---
%[output:3417f6f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4581e01c]
%   data: {"dataType":"text","outputData":{"text":"Case A stiffness:\n","truncated":false}}
%---
%[output:280d53ee]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{64\\,\\pi \\,\\textrm{D}}{3\\,R^2 }"}}
%---
%[output:6e0a9860]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:222a7955]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:53128d67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4b9019a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1b76cd87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:11014414]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0b5d8d26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:94e1bb69]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2dae094e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1983974e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:42ebfabb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:13ead4b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1967a5a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:72d5c4cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6c3c8a5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8040fa78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8f34bbc8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6e146d4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:206000bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:71311b91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:547e7ab3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:29c1aaa4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4e06b12e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6992b9ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5e0f0fed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8b325637]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:78d435e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5162f37e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5871ab63]
%   data: {"dataType":"text","outputData":{"text":"Dish #2 (Case A): M_d = 9 kg, K_d = 7.448e+05 N\/m, A_m = 1.0472 m^2, f_n = 45.785 Hz\n","truncated":false}}
%---
%[output:6a4faa08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:96c63b70]
%   data: {"dataType":"text","outputData":{"text":"Dish #2 (Case B): M_d = 19.58 kg, K_d = 2.719e+05 N\/m, A_m = 1.8678 m^2, f_n = 18.754 Hz\n","truncated":false}}
%---
%[output:4aff5ece]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:045d0012]
%   data: {"dataType":"text","outputData":{"text":"Dish #2 check: peak plate deflection, dynamic vs quasi-static\n","truncated":false}}
%---
%[output:61339828]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9b1f4b9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:549dc5cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7697b005]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:08c1fc54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:214808a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6fa828b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3843f62b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8d6cae1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:11d3d358]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6f22ee26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3b0eddb6]
%   data: {"dataType":"text","outputData":{"text":"            <strong>Plate<\/strong>             <strong>Elevation_deg<\/strong>    <strong>StartupPeak_um<\/strong>    <strong>SteadyPeak_um<\/strong>    <strong>QuasiStaticPeak_um<\/strong>    <strong>Dynamic_over_QuasiStatic<\/strong>\n","truncated":false}}
%---
%[output:71287fb0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:48fa193f]
%   data: {"dataType":"text","outputData":{"text":"    <strong>______________________<\/strong>    <strong>_____________<\/strong>    <strong>______________<\/strong>    <strong>_____________<\/strong>    <strong>__________________<\/strong>    <strong>________________________<\/strong>\n\n","truncated":false}}
%---
%[output:146670de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6ba23b93]
%   data: {"dataType":"text","outputData":{"text":"    \"Case A (clamped rim)\"          0               60.89           40.304              40.303                        1         \n","truncated":false}}
%---
%[output:6cef442b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:994ad0f4]
%   data: {"dataType":"text","outputData":{"text":"    \"Case A (clamped rim)\"         45              182.76            168.2               168.2                        1         \n","truncated":false}}
%---
%[output:6374ad02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7e9c66a3]
%   data: {"dataType":"text","outputData":{"text":"    \"Case A (clamped rim)\"         90              197.57           197.57              197.57                        1         \n","truncated":false}}
%---
%[output:43cac026]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2d52a723]
%   data: {"dataType":"text","outputData":{"text":"    \"Case B (clamped hub)\"          0              299.65           198.69              196.88                   1.0092         \n","truncated":false}}
%---
%[output:29ff9b26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3073f8ae]
%   data: {"dataType":"text","outputData":{"text":"    \"Case B (clamped hub)\"         45              894.34           822.95              821.67                   1.0016         \n","truncated":false}}
%---
%[output:89017ec5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3b22a42b]
%   data: {"dataType":"text","outputData":{"text":"    \"Case B (clamped hub)\"         90              965.13           965.13              965.13                        1         \n\n","truncated":false}}
%---
%[output:8bbc42f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:282d09a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:544dbbc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4a9bae7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:46ac682e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5f74f4ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:05202769]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:56ab5d16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1acf5954]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:56d0c7b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:02d3cfb5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2d9b656a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5ed0aba6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:579ceedb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:63fcb2cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3dc6cb3d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3a91cc4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5bb9428a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:965ca6fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:09815d6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:90180336]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8e9ab67b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:85f9a4e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:64b5005d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2833c2a6]
%   data: {"dataType":"text","outputData":{"text":"Physical-area RMS factors [raw, piston, defocus, residual]:\n","truncated":false}}
%---
%[output:1920200f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:09bb7a2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:35e4d1e8]
%   data: {"dataType":"text","outputData":{"text":"    0.4388    0.3240    0.2863    0.0745\n\n","truncated":false}}
%---
%[output:612eba1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:009a3847]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2c6f6342]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7e616dc5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:45cdd578]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:00e7bf71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9fa129d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:724c27a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7282728d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:40be6060]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:27015193]
%   data: {"dataType":"text","outputData":{"text":"Residual RMS factor per unit amplitude: Case A 0.07449 (centre), Case B 0.04622 (rim)\n","truncated":false}}
%---
%[output:9a13c586]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1f715cf8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4bbd129d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:36b22933]
%   data: {"dataType":"text","outputData":{"text":"             <strong>Model<\/strong>              <strong>f1_Hz<\/strong>     <strong>f1_over_f_gust<\/strong>\n    <strong>________________________<\/strong>    <strong>______<\/strong>    <strong>______________<\/strong>\n\n    \"Tower, lumped-model\"       17.307        34.615    \n    \"Tower, Ritz phi_t\"         17.307        34.615    \n    \"Tower, Ritz phi_line\"      17.522        35.045    \n    \"Tower, coupled (Kb-Kg)\"    17.303        34.607    \n    \"Disk A, one mode\"          45.785         91.57    \n    \"Disk B, one mode\"          18.754        37.509    \n\n","truncated":false}}
%---
%[output:8f17dd5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:46a0aa6b]
%   data: {"dataType":"text","outputData":{"text":"Ruze screening RMS limit at 225 GHz for eta >= 0.60: 75.78 um\n","truncated":false}}
%---
%[output:014a6c03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8b5c8b95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:166ed07f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9cd37f86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:85dcc2ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9d90d5fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8c1ddb14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:879659dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:502f0208]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3a469d4f]
%   data: {"dataType":"text","outputData":{"text":"Residual surface RMS, ","truncated":false}}
%---
%[output:432d1e2d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9c88c1be]
%   data: {"dataType":"text","outputData":{"text":"Disk B: clamped hub, one mode","truncated":false}}
%---
%[output:982e4ef0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2554db69]
%   data: {"dataType":"text","outputData":{"text":", elevation ","truncated":false}}
%---
%[output:0791430f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7abdf61b]
%   data: {"dataType":"text","outputData":{"text":"90","truncated":false}}
%---
%[output:5a74e9ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0cbc117f]
%   data: {"dataType":"text","outputData":{"text":" deg (L_b = ","truncated":false}}
%---
%[output:873516b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7d898731]
%   data: {"dataType":"text","outputData":{"text":"26.8 kg):\n  start state, gravity only (t = 0)          =  44.611 um\n  steady-state peak (t > 2.26 s)              =  44.611 um\n  worst case, gust peak adding to weight      =  44.611 um\n  Ruze limit                                  =   75.78 um\n  peak base-inertia force \/ weight (cos = 1)  = 5.11e-03\n  peak tower-tip slope (pointing error)       = 0.0172 arcmin\n","truncated":false}}
%---
%[output:933e7ed1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7cb153b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:678a484b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9d1e6c88]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:53bd2ecb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:94e2c576]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1fa97fac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3aa67aa0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3dead744]
%   data: {"dataType":"text","outputData":{"text":"(a) Height change of the paraboloid for a focal-length change, in terms of u:\n","truncated":false}}
%---
%[output:31ba9613]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"-\\frac{R^2 \\,\\delta_{\\textrm{focal}} \\,u}{4\\,{{\\textrm{focal}}_{\\textrm{len}} }^2 }"}}
%---
%[output:46acf3ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:13d15e49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:00567bb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:781cd027]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:37fcbea9]
%   data: {"dataType":"text","outputData":{"text":"(b) Projected ring area per unit u (a constant, so equal steps in u are equal areas):\n","truncated":false}}
%---
%[output:7f175c70]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:65688e6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:67210f7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4db80dab]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\pi \\,R^2"}}
%---
%[output:1fb6392a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:32e1ae44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:20eeb69d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0d9e26e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7d894329]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:588609c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9aa30371]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:70a84856]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:589293a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:151ab937]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:891b8cfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0d2a0f7f]
%   data: {"dataType":"text","outputData":{"text":"(c) True-area weight per projected area, in terms of u and beta:\n","truncated":false}}
%---
%[output:8967f26e]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\sqrt{\\beta \\,u+1}"}}
%---
%[output:3ead5be2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9dd974da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8e06af0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1c84ff78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1b4ab0fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:19a14edc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0655788b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7d0aa295]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7de62d5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4bfc4ea5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:84aeac3c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8e2a1e71]
%   data: {"dataType":"text","outputData":{"text":"Case A best fit (flat weighting): piston = ","truncated":false}}
%---
%[output:15acea20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1846de84]
%   data: {"dataType":"text","outputData":{"text":"5\/6","truncated":false}}
%---
%[output:650a4d28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:38462457]
%   data: {"dataType":"text","outputData":{"text":", defocus = ","truncated":false}}
%---
%[output:0c7e67f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:25d6776c]
%   data: {"dataType":"text","outputData":{"text":"-1","truncated":false}}
%---
%[output:5fe63709]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:58cdf1e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:168e6b08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8b10c884]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9ed670e7]
%   data: {"dataType":"text","outputData":{"text":"Case A residual e(u):\n","truncated":false}}
%---
%[output:88859b46]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"u^2 -u+\\frac{1}{6}"}}
%---
%[output:82db0c2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:98c1beb2]
%   data: {"dataType":"text","outputData":{"text":"Case A residual factor, flat-plate assumption (exact):\n","truncated":false}}
%---
%[output:00212b2f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{\\sqrt{5}}{30}"}}
%---
%[output:7f5211b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1e19f30f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2966a5b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8e020f81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1a0e2275]
%   data: {"dataType":"text","outputData":{"text":"Case A raw RMS factor = 5^(1\/2)\/5; refocusing removes 5\/6 = 83.3% of it\n","truncated":false}}
%---
%[output:212f101d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:54371559]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5d1bc094]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1cf9a091]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2b7a8153]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:888efc43]
%   data: {"dataType":"text","outputData":{"text":"Case A residual: e(0) = ","truncated":false}}
%---
%[output:1b7ec761]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7d1ec4e3]
%   data: {"dataType":"text","outputData":{"text":"1\/6","truncated":false}}
%---
%[output:4fade0bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9e89eebb]
%   data: {"dataType":"text","outputData":{"text":", e(1) = ","truncated":false}}
%---
%[output:4093966c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9aead909]
%   data: {"dataType":"text","outputData":{"text":"1\/6","truncated":false}}
%---
%[output:7bbae824]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:45c22e81]
%   data: {"dataType":"text","outputData":{"text":", e(","truncated":false}}
%---
%[output:50cc853a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:73d2e459]
%   data: {"dataType":"text","outputData":{"text":"1\/2","truncated":false}}
%---
%[output:2dcaa3cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:260dbaf7]
%   data: {"dataType":"text","outputData":{"text":") = -1\/12 (stationary point)\n","truncated":false}}
%---
%[output:8670b1e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:67ee1a3d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:20218990]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:50d9dd25]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0265610e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:14790efa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9add2ad8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:16ad56fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:50f34e01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2f0325b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:87e0d49f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:685295c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:826003bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:882420af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:65f0dab6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:32af5974]
%   data: {"dataType":"text","outputData":{"text":"              <strong>Flat_quadrature<\/strong>    <strong>Paraboloid_quadrature<\/strong>    <strong>Flat_mesh<\/strong>    <strong>Paraboloid_mesh<\/strong>\n              <strong>_______________<\/strong>    <strong>_____________________<\/strong>    <strong>_________<\/strong>    <strong>_______________<\/strong>\n\n    <strong>Case A<\/strong>       0.074536              0.074494           0.074536        0.074494    \n    <strong>Case B<\/strong>       0.046529              0.046222           0.046529        0.046222    \n\n","truncated":false}}
%---
%[output:998ba623]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:50cc1cf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:29d5b5ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6e793794]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8c14aaeb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:85a027cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:64c5f01c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:69257ae9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:96d870bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:68069839]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:929141a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:20650835]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:642ea8d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:373ca1fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:18602282]
%   data: {"dataType":"text","outputData":{"text":"Worst-case bending size eta_worst:\n","truncated":false}}
%---
%[output:5b0c1c55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1e4b4919]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{A_m \\,{\\left(g\\,m_A \\,\\sin_{\\textrm{elev},\\textrm{abs}} +C_p \\,k_n \\,p_{\\inf } \\,{\\left(g_f +1\\right)}\\right)}}{K_d }"}}
%---
%[output:9f2e5f2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4c8ead0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4ca67804]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:071c7e02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:351156b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:610f613a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:14f200be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:720065c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:092c5c33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:422bbed5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:321eb212]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:258db12e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:42189739]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9a243612]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:67369cb3]
%   data: {"dataType":"text","outputData":{"text":"Ruze RMS limit, solved from efficiency = eta_min:\n","truncated":false}}
%---
%[output:394e828f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{\\lambda_{\\textrm{rf}} \\,\\sqrt{-\\log \\left(\\eta_{\\min } \\right)}}{4\\,\\pi }"}}
%---
%[output:63454957]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2bfd36a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2361add0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:722846cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9e68fee2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9e8bcdfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:518ff14f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3d26dac0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0826dd97]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:147cd0a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8ef7b83e]
%   data: {"dataType":"text","outputData":{"text":"Surface RMS per reflector surface model, elevation ","truncated":false}}
%---
%[output:50f574b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:11e576da]
%   data: {"dataType":"text","outputData":{"text":"90","truncated":false}}
%---
%[output:08e109e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:77f1b85a]
%   data: {"dataType":"text","outputData":{"text":" deg, wind ","truncated":false}}
%---
%[output:41d8c077]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:56f7cc79]
%   data: {"dataType":"text","outputData":{"text":"6 m\/s, gust peak (Ruze limit 75.78 um):\n","truncated":false}}
%---
%[output:1974face]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8d8d80f2]
%   data: {"dataType":"text","outputData":{"text":"                          <strong>Raw_um<\/strong>    <strong>PistonRemoved_um<\/strong>    <strong>Residual_um<\/strong>    <strong>RuzeEfficiency<\/strong>\n                          <strong>______<\/strong>    <strong>________________<\/strong>    <strong>___________<\/strong>    <strong>______________<\/strong>\n\n    <strong>Case A, flat      <\/strong>    88.354         58.903           14.726           0.9809    \n    <strong>Case A, paraboloid<\/strong>    86.687         58.451           14.717          0.98092    \n    <strong>Case B, flat      <\/strong>    639.91          271.2           44.907          0.83579    \n    <strong>Case B, paraboloid<\/strong>    647.02         269.73           44.611          0.83776    \n\n","truncated":false}}
%---
%[output:24215b1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6fd05a4d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:523a015f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8a73abdf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7ebf4fc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3d19546e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:98adf10a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:610aedc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:381e3a65]
%   data: {"dataType":"text","outputData":{"text":"                          <strong>Raw_um<\/strong>    <strong>Piston_um<\/strong>    <strong>Defocus_um<\/strong>    <strong>Residual_um<\/strong>    <strong>Residual_share_pct<\/strong>\n                          <strong>______<\/strong>    <strong>_________<\/strong>    <strong>__________<\/strong>    <strong>___________<\/strong>    <strong>__________________<\/strong>\n\n    <strong>Case A, flat      <\/strong>    88.354     65.855        57.032        14.726             2.7778      \n    <strong>Case A, paraboloid<\/strong>    86.687     64.017        56.568        14.717             2.8824      \n    <strong>Case B, flat      <\/strong>    639.91      579.6        267.45        44.907            0.49248      \n    <strong>Case B, paraboloid<\/strong>    647.02     588.11        266.01        44.611            0.47539      \n\n","truncated":false}}
%---
%[output:5f2ab6c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3f87dc83]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9ad12b24]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8045a21c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:780e45a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6409e4a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0e251706]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8f6e4863]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4e6b99f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:85702a69]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3e2fa706]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:309bcd49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0cc2cc3b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:63cc2df3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1fd339d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4b20fe49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2c6997e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:40107e6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:21db08d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:07343fc3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:81b5ba64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:965451e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:76ef2bf6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4eb96057]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:21dbdaed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1f9671ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3137f95b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:09e23664]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:07e3db6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:14695a36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:27711343]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0913fea3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:12c12072]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:56543b81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4b2aa982]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:69089d7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:990a597f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0b423e72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0788f97d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:89b1f39e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7db6b46f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:36e63e55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8f3e015a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1463273b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9c78f762]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:91c7f05f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6dc10915]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9b4ebb0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:080920b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3c8fe4a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:22317ff6]
%   data: {"dataType":"text","outputData":{"text":"Worst-case residual RMS sigma_worst:\n","truncated":false}}
%---
%[output:0fa313d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6d8360fe]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{A_m \\,{\\textrm{residual}}_{\\textrm{factor}} \\,{\\left(g\\,m_A \\,\\sin_{\\textrm{elev},\\textrm{abs}} +C_p \\,k_n \\,p_{\\inf } \\,{\\left(g_f +1\\right)}\\right)}}{K_d }"}}
%---
%[output:0ebf9f35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:98ba01c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3e98e5aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6985a994]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:768f3467]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:29fc20f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:92397eac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4dc70225]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4e80adac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0e933e29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:853287b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:964f9592]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9fce0104]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:33c7f5f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:203c694d]
%   data: {"dataType":"text","outputData":{"text":"Residual RMS at elevation 90 deg, wind 6 m\/s (Ruze limit 75.78 um):\n","truncated":false}}
%---
%[output:7234077c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2578dbec]
%   data: {"dataType":"text","outputData":{"text":"                             <strong>SteadyState_um<\/strong>    <strong>Peak_um<\/strong>    <strong>Change_vs_flat_pct<\/strong>    <strong>RuzeEfficiency_at_peak<\/strong>    <strong>Within_Ruze_limit<\/strong>\n                             <strong>______________<\/strong>    <strong>_______<\/strong>    <strong>__________________<\/strong>    <strong>______________________<\/strong>    <strong>_________________<\/strong>\n\n    <strong>Case A, flat plate   <\/strong>        14.726        14.726                 0                 0.9809                  true       \n    <strong>Case A, true parabol…<\/strong>        14.717        14.717         -0.056207                0.98092                  true       \n    <strong>Case B, flat plate   <\/strong>        44.907        44.907                 0                0.83579                  true       \n    <strong>Case B, true parabol…<\/strong>        44.611        44.611          -0.65926                0.83776                  true       \n\n","truncated":false}}
%---
%[output:7690d45b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:16dfc757]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2c75789c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:669fbf34]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2da62fe4]
%   data: {"dataType":"text","outputData":{"text":"Worst-case residual RMS versus wind speed V:\n","truncated":false}}
%---
%[output:3d8db93a]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{A_m \\,{\\textrm{residual}}_{\\textrm{factor}} \\,{\\left(\\frac{C_p \\,k_n \\,\\rho_{\\textrm{air}} \\,{\\left(g_f +1\\right)}\\,{{\\textrm{wind}}_V }^2 }{2}+g\\,m_A \\,\\sin_{\\textrm{elev},\\textrm{abs}} \\right)}}{K_d }"}}
%---
%[output:869e28a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:87960db1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:518e2c61]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:41ee0220]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9268efeb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:485ef88a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7bb895bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7d18e504]
%   data: {"dataType":"text","outputData":{"text":"Gravity part (V = 0):\n","truncated":false}}
%---
%[output:098c5a17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8e6db08f]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{A_m \\,g\\,m_A \\,{\\textrm{residual}}_{\\textrm{factor}} \\,\\sin_{\\textrm{elev},\\textrm{abs}} }{K_d }"}}
%---
%[output:6535cc22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:44cc79d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:338b4cc8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7587017e]
%   data: {"dataType":"text","outputData":{"text":"Wind part:\n","truncated":false}}
%---
%[output:29e87d51]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:57f6b99a]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{A_m \\,C_p \\,k_n \\,{\\textrm{residual}}_{\\textrm{factor}} \\,\\rho_{\\textrm{air}} \\,{{\\textrm{wind}}_V }^2 \\,{\\left(g_f +1\\right)}}{2\\,K_d }"}}
%---
%[output:270d2ed2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9aa94e7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1771c7d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1882bb95]
%   data: {"dataType":"text","outputData":{"text":"Doubling the wind speed multiplies the wind part by 4\n","truncated":false}}
%---
%[output:2f6d9085]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:89effa97]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3b132792]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3d9d8beb]
%   data: {"dataType":"text","outputData":{"text":"Equivalent steady wind = (g_f + 1)^(1\/2) x gust mean speed = 1.140 x for g_f = 0.30\n","truncated":false}}
%---
%[output:98607d6d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:588540fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7185a8fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:17c5017d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:70ae8c77]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:11b255a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:222c4708]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:2bada1e5]
%   data: {"dataType":"text","outputData":{"text":"Critical wind speed V_crit:\n","truncated":false}}
%---
%[output:082e5e86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1e5ee7de]
%   data: {"dataType":"symbolic","outputData":{"name":"","value":"\\frac{\\sqrt{2}\\,\\sqrt{K_d }\\,\\sqrt{\\frac{\\sigma_{\\max ,\\textrm{sym}} -\\frac{A_m \\,g\\,m_A \\,{\\textrm{residual}}_{\\textrm{factor}} \\,\\sin_{\\textrm{elev},\\textrm{abs}} }{K_d }}{k_n }}}{\\sqrt{A_m }\\,\\sqrt{C_p }\\,\\sqrt{{\\textrm{residual}}_{\\textrm{factor}} }\\,\\sqrt{\\rho_{\\textrm{air}} }\\,\\sqrt{g_f +1}}"}}
%---
%[output:02d6fb1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:67dd7d90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1e6474ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8dbe2bf8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8b357688]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:075409db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0e6f3390]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:794085c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3fda3d1c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7aa8c2f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1efcdc49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5d81d09c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:005a3d41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:536ddcb5]
%   data: {"dataType":"text","outputData":{"text":"    <strong>Elevation_deg<\/strong>    <strong>RMS_at_0ms_um<\/strong>    <strong>RMS_at_design_wind_um<\/strong>    <strong>Vcrit_ms<\/strong>\n","truncated":false}}
%---
%[output:6e01a8ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6b0e9ce6]
%   data: {"dataType":"text","outputData":{"text":"    <strong>_____________<\/strong>    <strong>_____________<\/strong>    <strong>_____________________<\/strong>    <strong>________<\/strong>\n\n","truncated":false}}
%---
%[output:6406317a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4769d9db]
%   data: {"dataType":"text","outputData":{"text":"          0                  0               9.1608             17.257 \n","truncated":false}}
%---
%[output:6e14cf74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:64a526e5]
%   data: {"dataType":"text","outputData":{"text":"         45             31.754               38.231             15.643 \n","truncated":false}}
%---
%[output:85437dc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:85e31600]
%   data: {"dataType":"text","outputData":{"text":"         90             44.907               44.907                Inf \n","truncated":false}}
%---
%[output:7ccd7db4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5bb2b523]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:93b82c65]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:105a419c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:1b3d5686]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:979f87b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:43eb49c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5c45bbe8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:3ff92b1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:8b2af639]
%   data: {"dataType":"text","outputData":{"text":"    <strong>Elevation_deg<\/strong>    <strong>Vcrit_flat_ms<\/strong>    <strong>Vcrit_paraboloid_ms<\/strong>\n    <strong>_____________<\/strong>    <strong>_____________<\/strong>    <strong>___________________<\/strong>\n\n          0             17.257              17.314       \n         45             15.643              15.732       \n         90                Inf                 Inf       \n\n","truncated":false}}
%---
%[output:9e4fd169]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:673cc0e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:62233bd6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:79ef6561]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:44e1d1c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:40fe681c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5e8d7843]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:52874a2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:7f342048]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:188c1e91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:06db1b34]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0ad942b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9d614c91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:4b437120]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:60948342]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:97872dfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:72979fed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:511f9c5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6b87809c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:320b45b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:0422c1de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:26d1eee1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:90bafddd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:74f09957]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5657a1a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:5b02579a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:664e5bf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:302186c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:6b3a8940]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
%[output:9e46565a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Unable to locate a personal folder for $documents\\MATLAB."}}
%---
