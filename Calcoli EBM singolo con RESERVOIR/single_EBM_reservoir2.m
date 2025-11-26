function [Um_ebm, Um_res, C_ebm] = single_EBM_reservoir2(h, rc, params)
% Compute quantities in the form of "maps", i.e. matrices for values of Um_ebm,
% Um_res, C_ebm computed for the given domain. Returns a column (fixed rc,
% scroll h)
% INPUTS:
% state = [h, rc]
% params = struct of parameters (everything that is not [h, rc])
% OUTPUTS:
% matrices Um_ebm, Um_res, C_ebm for the given state and parameters

% initialization
ri = params.ri; % [m]
ro = params.ro; % [m]
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

Um_ebm = zeros(size(h)); Um_res = zeros(size(h)); C_ebm = zeros(size(h));
lambda1_res_v = zeros(size(h)); lambda2_res_v = zeros(size(h));

% fixed rc, scroll columns (h)
alpha = alpha_sol(h, rc, params); % compute alpha (vector) for all h values, fixed rc (pick alpha(i) in the for)
b = sqrt((r2r-r1r)^2 + hr.^2); % compute b (vector) for all h values, fixed rc (pick b(i) in the for)
l = b + 8/3*alpha.^2./b - 32/15*alpha.^4./(b).^3; % compute l
r1 = linspace(ri,rc,20)';   r2 = linspace(rc,ro,20)'; % for integration
r3 = linspace(ri,ro,40)';
l0_ebm = ro-ri; % initial lenght, fully zipped 
% 1/2 LEAP
%uc_ebm = -h.^2/(2*(rc-ri)^2)/(1/(rc-ri)-(rc^2+ro^2)/(rc^2-ro^2)/rc);

for i = 1:numel(h) % scroll columns
    % ACTUATOR: compute Um and C
    % simplified epsilon1 for ebm
    l2_ebm = ro-rc;
    l1_ebm = sqrt(h(i).^2 + (rc-ri)^2);
    l_ebm = l2_ebm + l1_ebm;
    eps1_ebm = l_ebm./l0_ebm - 1; % meridian
    % compute Um_ebm
    %Um_ebm(i) = pi*t0*E/(1-nu^2) * (trapz(r1, r1.*(eps1_1.^2)) + trapz(r2, r2.*(eps1_2.^2)));
    Um_ebm(i) = pi*t0*E/(1-nu^2) * (trapz(r3, r3.*(eps1_ebm.^2)));
    % compute C
    %Cp = 2*pi*epsilon_p/t0 * trapz( r2, (1+eps1_2).*(1+eps2_2)./(1+eps3_2).*r2 ); % polymer capacitance
    C_ebm(i) = ((epsilon_p*epsilon_o*pi*(ro^2-rc^2)) / (2*epsilon_o*t0 + epsilon_p*to_left));
    % if to_left > 0
    %      Co = 2*pi*epsilon_o/to_left * trapz(r2, (1+eps1_2).*(1+eps2_2).*r2); % oil capacitance
    %      C_ebm(i) = Cp*Co/(Cp+Co); % oil capacitance is in series with poly
    % else
    %      C_ebm(i) = Cp;
    % end

    % RESERVOIR
    lambda_p = R1/R0; % preload factor
    dR0 = (r2r - r1r)/lambda_p; % anulus radius at rest
    A0_res = pi*((r2r/lambda_p)^2 - (r1r/lambda_p)^2); % area at rest
    lambda1_res = l(i)/dR0; % [-] meridian stretch
    lambda2_res = lambda_p; % [-] circumferencial stretch = 80/50 (PRELOAD)
    Psi_res = mu/2*(lambda1_res^2 + lambda2_res^2 + (lambda1_res*lambda2_res)^-2 - 3); % [J/m^3] energy density 
    Um_res(i)   = A0_res*t0_res*Psi_res; % [J] energy

    % to workspace
    lambda1_res_v(i) = lambda1_res;
    lambda2_res_v(i) = lambda2_res;

end

end