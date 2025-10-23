clear; close all; clc

%% PARAMETERS
Y = 1.25e9;              % [Pa]
nu = 0.30;
epsP = 3.9*8.854e-12;    % [F/m] polymer
epsO = 2.7*8.854e-12;    % [F/m] oil
params.Y = Y; params.nu = nu; params.epsP = epsP; params.epsO = epsO;

ri = 8e-3;               % [m]
ro = 15e-3;              % [m]
h0 = 1.5e-3;             % [m] initial height of single actuator (FULL HEIGHT, not 1/2)
t0 = 25.4e-6;            % [m] polymer thickness
to_res = 5e-6/2;         % [m] residual oil film (half)
params.ri = ri; params.ro = ro; params.h0 = h0; params.t0 = t0; params.to_res = 5e-6/2;

% dominio: supponiamo che EE si sposti verso il basso di
% massimo "h0" mm (fino a "fondo corsa")
% supponiamo che rc possa variare da ro = 15 mm fino a 9 mm
lh = 100;
lr = 100;
h_vec = linspace(0.001e-5, h0-0.001e-5, lh)';
rc_vec = linspace(ro*1.00000001, ro - 6e-3, lr);

%% MAPS
% building the "maps" on the domain [h_vec, rc_vec] to solve the equations
% numerically

Um_bot = zeros(lh, lr); Um_top = zeros(lh, lr); C_bot = zeros(lh, lr); alpha = zeros(lh, lr);

for i = 1:lr
    [Um_bot(:, i), Um_top(:, i), C_bot(:, i), alpha(:, i)] = double_EBM(h_vec, rc_vec(i), params);
end

% DOUBLE CONES
Um_bott = 2*Um_bot;
Um_topt = 2*Um_top;
C_bott = 0.5*C_bot; % series
% NOTE: results from double_EBM.m are reasonable in terms of TREND of
% Um_bott(h,rc) and Um_topt(h,rc). Also C(h,rc) makes sense as well as
% alpha(h,rc). The only doubt is the values of Um are correct or not. The
% orders of magnitude seem reasonable when compared to the single EBM model

% midpoints for the midpoint grids
rci = 0.5*(rc_vec(1:end-1) + rc_vec(2:end));
hi  = 0.5*(h_vec(1:end-1)  + h_vec(2:end));

% d/drc at midpoint nodes in rc (fix h, differences in rc)
dUmtop_drc = zeros(lh, lr-1);
dUmbot_drc = zeros(lh, lr-1);
dCbot_drc = zeros(lh, lr-1);
% compute numerical derivatives in rc at the midpoints
for i = 1:lh
    dUmtop_drc(i, :) = diff(Um_topt(i, :))./diff(rc_vec);
    dUmbot_drc(i, :) = diff(Um_bott(i, :))./diff(rc_vec);
    dCbot_drc(i, :) = diff(C_bott(i, :))./diff(rc_vec);
end
% interpolate in h -> hi (midpoint in h) to align d/drc with (hi, rci)
dUmtop_drci = zeros(lh-1, lr-1);
dUmbot_drci = zeros(lh-1, lr-1);
dCbot_drci = zeros(lh-1, lr-1);
for i = 1 : lr-1
    dUmtop_drci(:,i) = interp1(h_vec, dUmtop_drc(:,i), hi, 'linear');
    dUmbot_drci(:,i) = interp1(h_vec, dUmbot_drc(:,i), hi, 'linear');
    dCbot_drci(:,i) = interp1(h_vec, dCbot_drc(:,i), hi, 'linear');
end

% d/dh at midpoint nodes in h (fix rc, differences in h)
dUmtop_dh = zeros(lh-1, lr);
dUmbot_dh = zeros(lh-1, lr);
dCbot_dh = zeros(lh-1, lr);
% compute numerical derivatives in rc at the midpoints
for i = 1:lr
    dUmtop_dh(:, i) = diff(Um_topt(:, i))./diff(h_vec);
    dUmbot_dh(:, i) = diff(Um_bott(:, i))./diff(h_vec);
    dCbot_dh(:, i) = diff(C_bott(:, i))./diff(h_vec);
end
% interpolate in rc -> rci (midpoint in rc) to align d/dh with (hi, rci)
dUmtop_dhi = zeros(lh-1, lr-1);
dUmbot_dhi = zeros(lh-1, lr-1);
dCbot_dhi = zeros(lh-1, lr-1);
% compute numerical derivatives in h at the midpoints
for i = 1:lr-1
    dUmtop_dhi(i, :) = interp1(rc_vec, dUmtop_dh(i, :), rci, 'linear');
    dUmbot_dhi(i, :) = interp1(rc_vec, dUmbot_dh(i, :), rci, 'linear');
    dCbot_dhi(i, :) = interp1(rc_vec, dCbot_dh(i, :), rci, 'linear');
end

%% SOLVE THE EQUATIONS
