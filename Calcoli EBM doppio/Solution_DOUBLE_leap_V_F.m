clear; close all; clc

%% PARAMETERS
Y = 1.25e9;              % [Pa]
nu = 0.30;
epsP = 3.9*8.854e-12;    % [F/m] polymer
epsO = 2.7*8.854e-12;    % [F/m] oil
params.Y = Y; params.nu = nu; params.epsP = epsP; params.epsO = epsO;

ri = 8e-3;               % [m]
ro = 15e-3;              % [m]
% h0 = 1.5e-3;             % [m] initial height of single actuator (FULL HEIGHT, not 1/2)
h0 = 1.5e-3;
t0 = 25.4e-6;            % [m] polymer thickness
to_res = 5e-6/2;         % [m] residual oil film (half)
params.ri = ri; params.ro = ro; params.h0 = h0; params.t0 = t0; params.to_res = to_res;

% dominio: supponiamo che EE si sposti verso il basso di
% massimo "h0" mm (fino a "fondo corsa")
% supponiamo che rc possa variare da ro = 15 mm fino a 9 mm
lh = 100;
lr = 100;
h_vec = linspace(1e-10, h0-1e-10, lh)';
rc_vec = linspace(ro - 6e-3, ro*0.9999999999, lr);

%% MAPS
% building the "maps" on the domain [h_vec, rc_vec] to solve the equations
% numerically

Um_bot = zeros(lh, lr); Um_top = zeros(lh, lr); C_bot = zeros(lh, lr); alpha = zeros(lh, lr);
eps1_cone = zeros(lh, lr);
eps1_bulge = zeros(lh, lr);
eps1_top = zeros(lh, lr);
eps1_single = zeros(lh, lr);

for i = 1:lr
    [Um_bot(:, i), Um_top(:, i), C_bot(:, i), alpha(:, i), eps1_cone(:, i), eps1_bulge(:, i), eps1_top(:, i), eps1_single(:, i)] = double_EBM(h_vec, rc_vec(i), params);
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


%% Solve system for F(h)
% eq_h_map(h,rc) = dU/dh = dUmbot_dhi + dUmtop_dhi  (size: [lh-1, lr-1]),
% on the midpoints basically
eq_h_map = dUmbot_dhi + dUmtop_dhi;

V_min = 0.1e3; V_max = 8e3;
V_list = linspace(V_min, V_max, 5);  lv = numel(V_list);

clr = lines(numel(V_list));

figure(1);
subplot(2,1,1); hold on; grid on; grid minor;
xlabel('h [mm]'); ylabel('F(h) [N]'); title('F(h) (F>0 upwards) at different V levels');

rc_curves = cell(lv,1);

for ivf = 1:numel(V_list)
    V0 = V_list(ivf);

    rc_star_of_h = NaN(size(hi));   % r_c^*(h;V0) on midpoints hi
    F_of_h       = NaN(size(hi));   % F(h;V0) = dU/dh evaluated in (h, r_c^*(h))
    rc_prev = median(rci);          % warm-start

    for ih = 1:numel(hi)
        % residual in rc at h = hi
        dUdr_row = dUmbot_drci(ih,:) + dUmtop_drci(ih,:);   % [N/m]
        dCdr_row = dCbot_drci(ih,:);                        % [F/m]
        Rrc = dUdr_row - 0.5*V0^2 .* dCdr_row;

        xr = rci; ok = isfinite(Rrc);
        if nnz(ok) < 2, continue; end
        Rrc = Rrc(ok); xr = xr(ok);

        % look for a change in sign
        sgn = Rrc(1:end-1).*Rrc(2:end);
        idx = find(sgn <= 0);
        if ~isempty(idx)
            [~,pick] = min(abs(xr(idx) - rc_prev));
            k = idx(pick);
            rc_sol = fzero(@(x) interp1(xr, Rrc, x, 'linear','extrap'), [xr(k), xr(k+1)]);
        else
            % fallback
            [~,kk] = min(abs(Rrc)); rc_sol = xr(kk);
        end
        
        % clamp to domain
        rc_sol = min(max(rc_sol, rci(1)), rci(end));

        % warm start
        rc_star_of_h(ih) = rc_sol;
        rc_prev = rc_sol;
        rc_curves{ivf} = rc_star_of_h;   % store the curve

        % F(h) = dU/dh in (h=hi(ih), rc=rc_sol) compute F at h,rc_sol
        F_of_h(ih) = interp1(rci, eq_h_map(ih,:), rc_sol, 'pchip', 'extrap');  % [N]
    end

    % plot NOTE: THE SIGN OF F IS INVERTED; from the equations F is
    % directed as h, so if we invert it: F positive means force upwards

    plot(hi*1e3, -F_of_h, '-', 'LineWidth', 1.3, 'Color', clr(ivf,:), ...
        'DisplayName', sprintf('V = %.1f kV', V0*1e-3));
    
end

legend('Location','best');

subplot(2,1,2)
hold on; grid on; grid minor;
xlabel('h [mm]'); ylabel('r_c(h) [mm]'); 
title('rc(h) at different V levels');

for ivf = 1:numel(V_list)
    rc_star_of_h = rc_curves{ivf};
    plot(hi*1e3, rc_star_of_h*1e3, '-', 'LineWidth', 1.3, 'Color', clr(ivf,:), ...
        'DisplayName', sprintf('V = %.1f kV', V_list(ivf)/1e3));
end

legend('Location','best');