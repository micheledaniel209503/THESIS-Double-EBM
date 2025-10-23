function [Um_bot, Um_top, C_bot, alpha] = double_EBM(h, rc, params)
% Compute quantities in the form of "maps", i.e. matrices for values of Um_bot,
% Um_top, C_bot computed for the given domain. Returns a column (fixed rc,
% scroll h)
% INPUTS:
% state = [h, rc]
% params = struct of parameters (everything that is not [h, rc]
% OUTPUTS:
% matrices Um_bot, Um_top, C_bot for the given state and parameters

% initialization
h0 = params.h0; % [m] initial height of an actuator (FULL HEIGHT, not 1/2)
ri = params.ri; % [m]
ro = params.ro; % [m]
nu = params.nu; % Poisson coefficient
E = params.Y;   % Elastic module
t0 = params.t0; % polymer thickness
epsilon_p = params.epsP;
epsilon_o = params.epsO;
to_left = params.to_res;

Um_bot = zeros(size(h)); Um_top = zeros(size(h)); C_bot = zeros(size(h));

% fixed rc, scroll columns (h)
alpha = alpha_sol(h, rc, params); % compute alpha (vector) for all h values, fixed rc (pick alpha(i) in the for)
b = sqrt((ro-ri)^2 + ((h+h0)/2).^2); % compute b (vector) for all h values, fixed rc (pick b(i) in the for)
r1 = linspace(ri,rc,20)';   r2 = linspace(rc,ro,20)';  r3 = linspace(ri,ro,40)';% for integration

% 1/2 LEAP
% !!!!!!!!!!!!!! NOTE: MUST SUBSTITUTE h WITH THE PROPER hbot !!!!!!!!!!!!!!!!
h_bot = (h0 - h)/2; % 1/2 HALF height of the bottom actuator 
h_top = (h0 + h)/2; % 1/2 HALF height of the top actuator
uc_bot = -(h_bot).^2/(2*(rc-ri)^2)/(1/(rc-ri)-(rc^2+ro^2)/(rc^2-ro^2)/rc);
%fake_rc_ro = 0.9999999999*ro;
fake_rc_ro = 1.0000000001*ro;
uc_top = -(h_top).^2/(2*(fake_rc_ro-ri)^2)/(1/(fake_rc_ro-ri)-(fake_rc_ro^2+ro^2)/(fake_rc_ro^2-ro^2)/fake_rc_ro);

for i = 1:numel(h) % scroll columns
    % BOTTOM ACTUATOR: compute Umbot and Cbot
    % internal region r1 in [ri, rc]
    eps1_1 = uc_bot(i)/(rc-ri) + (h_bot(i).^2)/(2*(rc-ri)^2);  % meridian
    eps2_1 = 0;                                            % circumferential
    eps3_1 = 0;                                            % normal
    % external region r2 in [rc, ro]
    eps1_2 = -uc_bot(i)*rc/(ro^2-rc^2) * (1 + ro^2./r2.^2);     % meridian
    eps2_2 = 0;                                             % circumferential
    eps3_2 = 0;                                             % normal
    % compute Umbot
    Um_bot(i) = pi*t0*E/(1-nu^2) * (trapz(r1, r1.*(eps1_1.^2 + eps2_1.^2 + 2*nu*eps1_1.*eps2_1)) + trapz(r2, r2.*(eps1_2.^2 + eps2_2.^2 + 2*nu*eps1_2.*eps2_2)) );
    % compute Cbot
    Cp = 2*pi*epsilon_p/t0 * trapz( r2, (1+eps1_2).*(1+eps2_2)./(1+eps3_2).*r2 ); % polymer capacitance
    if to_left > 0
         Co = 2*pi*epsilon_o/to_left * trapz(r2, (1+eps1_2).*(1+eps2_2).*r2); % oil capacitance
         C_bot(i) = Cp*Co/(Cp+Co); % oil capacitance is in series with poly
    else
         C_bot(i) = Cp;
    end

    % TOP ACTUATOR: compute Umtop
    eps1_bulge = 8/3*(alpha(i)/b(i))^2 + 32/15*(alpha(i)/b(i))^4; % strain expression for the bulge
    eps1_cone = uc_top(i)/(fake_rc_ro-ri) + (h_top(i).^2)/(2*(fake_rc_ro-ri)^2);
    eps1_top = eps1_bulge + eps1_cone;
    % compute Umtop
    Um_top(i) = pi*t0*E/(1-nu^2)*(trapz(r3, (eps1_top).^2.*r3));

end

end