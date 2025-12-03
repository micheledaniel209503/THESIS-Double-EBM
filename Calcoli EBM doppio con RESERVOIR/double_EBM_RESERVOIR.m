function [Um_ebm_B, Um_ebm_T, Um_res, C_bot, alpha] = double_EBM_RESERVOIR(h, rc, params)
% Compute quantities in the form of "maps", i.e. matrices for values of Um_bot,
% Um_top, C_bot computed for the given domain. Returns a column (fixed rc,
% scroll h)
% INPUTS:
% state = [h, rc]
% params = struct of parameters (everything that is not [h, rc]
% OUTPUTS:
% matrices Um_bot, Um_top, C_bot for the given state and parameters

% initialization
ri = params.ri; % [m]
ro = params.ro; % [m]
h0 = params.h0; % [m]
r1r = params.r1r; % [m]
r2r = params.r2r; % [m]
hr = params.hr; % [m]
R0 = params.R0; % [m]
R1 = params.R1; % [m]
nu = params.nu; % Poisson coefficient
E = params.Y;   % Elastic module
t0 = params.t0; % [m] polymer thickness
t0_res = params.t0_res; % [m] reservoir (VHB) thickness at the unstretched state (before manufacturing reservoir)
epsilon_p = params.epsP; % [F/m] permittivity of polymer
epsilon_o = params.epsO; % [F/m] permittivity of oil
to_left = params.to_res; % [m] residual oil thickness
mu = params.mu; % [Pa] Neo-Hook parameter for VHB 

% initialization
Um_ebm_B = zeros(size(h)); Um_ebm_T = zeros(size(h)); Um_res = zeros(size(h)); C_bot = zeros(size(h));

% fixed rc, scroll columns (h)
alpha = alpha_sol(h, rc, params); % compute alpha (vector) for all h values, fixed rc (pick alpha(i) in the for)
b = sqrt((r2r-r1r)^2 + hr.^2); % compute b (constant) value of initial length of reservoir side CONFIGURATION D
l = b + 8/3*alpha.^2./b - 32/15*alpha.^4./(b).^3; % compute l generic length of reservoir side CONFIGURATION 1
r3 = linspace(ri,ro,40)'; % for integration

% 1/2 LEAP
h_bot = (h0 - h)/2; % 1/2 HALF height of the bottom actuator 
h_top = (h0 + h)/2; % 1/2 HALF height of the top actuator
l0_ebm = ro - ri; % initial lenght, fully zipped (both ebms) CONFIGURATION D

for i = 1:numel(h) % scroll columns
    % BOTTOM ACTUATOR: compute Umbot and Cbot
    % simplified epsilon1 for ebm
    l2_ebm_B = ro-rc;
    l1_ebm_B = sqrt(h_bot(i).^2 + (rc-ri)^2); % h_bot i already 1/2
    l_ebm_B = l2_ebm_B + l1_ebm_B;
    eps1_ebm_B = l_ebm_B./l0_ebm - 1; % meridian
    % compute Umbot
    Um_ebm_B(i) = pi*t0*E/(1-nu^2) * (trapz(r3, r3.*(eps1_ebm_B.^2)));

    % compute Cbot
    % CHECK
    % THISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    % THIS WAS TAKEN FROM SINGLE_EBM_RESERVOIR2
    C_bot(i) = ((epsilon_p*epsilon_o*pi*(ro^2-rc^2)) / (2*epsilon_o*t0 + epsilon_p*to_left));

    % TOP ACTUATOR: compute Um_top
    l_ebm_T = sqrt((ro-ri)^2 + (h_top(i)).^2);
    eps1_ebm_T = l_ebm_T./l0_ebm - 1;
    Um_ebm_T(i) = pi*t0*E/(1-nu^2)*(trapz(r3, (eps1_ebm_T).^2.*r3)); 

    % RESERVOIR: compute Um_res
    lambda_p = R1/R0; % preload factor
    dR0 = (r2r - r1r)/lambda_p; % anulus radius at rest CONFIGURATION RD
    A0_res = pi*((r2r/lambda_p)^2 - (r1r/lambda_p)^2); % area at rest CONFIGURATION RD
    lambda1_res = l(i)/dR0; % [-] meridian stretch of VHB ring
    lambda2_res = lambda_p; % [-] circumferencial stretch = 80/50 (PRELOAD)
    Psi_res = mu/2*(lambda1_res^2 + lambda2_res^2 + (lambda1_res*lambda2_res)^-2 - 3); % [J/m^3] energy density 
    Um_res(i)   = A0_res*t0_res*Psi_res; % [J] energy of VHB reservoir

end

end