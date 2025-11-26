clear; close all; clc

%% PARAMETERS
Y = 1.25e9;              % [Pa]
nu = 0.30;
epsP = 3.9*8.854e-12;    % [F/m] polymer
epsO = 2.7*8.854e-12;    % [F/m] oil
ri = 8e-3;               % [m]
ro = 15e-3;              % [m]
h0 = 1.5e-3;             % [m] initial preload height (FULL HEIGHT OF A SINGLE EBM)
r1r = ri;
%r2r = ro + 1.5e-3;
r2r = 23e-3;
hr = 10e-3;             % [m] height of the VHB reservoir (internal cylindrical support)
%hr = 1e-6;
t0 = 25.4e-6;            % [m] keplan thickness
t0_res = 0.5e-3;         % [m] thickness of unstretched VHB disc (reservoir)
to_res = 5e-6/2;         % [m] residual oil film (half)
%mu = 5e5;
mu = 1.83e5;             % [Pa] Neo-Hook's parameter for VHB reservoir @ approx 20°C
%mu = 1e4;
R0 = 50e-3/2;            % [m] initial radius of the unstretched VHB disc
R1 = 80e-3/2;            % [m] final radius of the stretched VHB disc (flat)

params.Y = Y; params.nu = nu; params.epsP = epsP; params.epsO = epsO;
params.ri = ri; params.ro = ro; params.h0 = h0; params.r1r = r1r; params.r2r = r2r; params.hr = hr; params.t0 = t0; params.t0_res = t0_res;
params.to_res = to_res; params.mu = mu; params.R0 = R0; params.R1 = R1;

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

Um_bot = zeros(lh, lr); Um_top = zeros(lh, lr); Um_res = zeros(lh, lr); C_bot = zeros(lh, lr); alpha = zeros(lh, lr);

for i = 1:lr
    [Um_bot(:, i), Um_top(:, i), Um_res(:, i), C_bot(:, i), alpha(:, i)] = double_EBM_RESERVOIR(h_vec, rc_vec(i), params);
end

%%
% DOUBLE CONES
Um_bott = 2*Um_bot;
Um_topt = 2*Um_top;
Um_rest = Um_res; % already complete
C_bott = C_bot; % already complete

% midpoints for the midpoint grids
rci = 0.5*(rc_vec(1:end-1) + rc_vec(2:end));
hi  = 0.5*(h_vec(1:end-1)  + h_vec(2:end));

% d/drc at midpoint nodes in rc (fix h, differences in rc)
dUmtop_drc = zeros(lh, lr-1);
dUmbot_drc = zeros(lh, lr-1);
dUmres_drc = zeros(lh, lr-1);
dCbot_drc = zeros(lh, lr-1);

% compute numerical derivatives in rc at the midpoints
for i = 1:lh
    dUmtop_drc(i, :) = diff(Um_topt(i, :))./diff(rc_vec);
    dUmbot_drc(i, :) = diff(Um_bott(i, :))./diff(rc_vec);
    dUmres_drc(i, :) = diff(Um_rest(i, :))./diff(rc_vec);
    dCbot_drc(i, :) = diff(C_bott(i, :))./diff(rc_vec);
end
% interpolate in h -> hi (midpoint in h) to align d/drc with (hi, rci)
dUmtop_drci = zeros(lh-1, lr-1);
dUmbot_drci = zeros(lh-1, lr-1);
dUmres_drci = zeros(lh-1, lr-1);
dCbot_drci = zeros(lh-1, lr-1);
for i = 1 : lr-1
    dUmtop_drci(:,i) = interp1(h_vec, dUmtop_drc(:,i), hi, 'linear');
    dUmbot_drci(:,i) = interp1(h_vec, dUmbot_drc(:,i), hi, 'linear');
    dUmres_drci(:,i) = interp1(h_vec, dUmres_drc(:,i), hi, 'linear');
    dCbot_drci(:,i) = interp1(h_vec, dCbot_drc(:,i), hi, 'linear');
end

% d/dh at midpoint nodes in h (fix rc, differences in h)
dUmtop_dh = zeros(lh-1, lr);
dUmbot_dh = zeros(lh-1, lr);
dUmres_dh = zeros(lh-1, lr);
dCbot_dh = zeros(lh-1, lr);
% compute numerical derivatives in rc at the midpoints
for i = 1:lr
    dUmtop_dh(:, i) = diff(Um_topt(:, i))./diff(h_vec);
    dUmbot_dh(:, i) = diff(Um_bott(:, i))./diff(h_vec);
    dUmres_dh(:, i) = diff(Um_rest(:, i))./diff(h_vec);
    dCbot_dh(:, i) = diff(C_bott(:, i))./diff(h_vec);
end
% interpolate in rc -> rci (midpoint in rc) to align d/dh with (hi, rci)
dUmtop_dhi = zeros(lh-1, lr-1);
dUmbot_dhi = zeros(lh-1, lr-1);
dUmres_dhi = zeros(lh-1, lr-1);
dCbot_dhi = zeros(lh-1, lr-1);
% compute numerical derivatives in h at the midpoints
for i = 1:lr-1
    dUmtop_dhi(i, :) = interp1(rc_vec, dUmtop_dh(i, :), rci, 'linear');
    dUmbot_dhi(i, :) = interp1(rc_vec, dUmbot_dh(i, :), rci, 'linear');
    dUmres_dhi(i, :) = interp1(rc_vec, dUmres_dh(i, :), rci, 'linear');
    dCbot_dhi(i, :) = interp1(rc_vec, dCbot_dh(i, :), rci, 'linear');
end

%% 1) Solve eq_h(h,rc) for h_eq(rc) (indipendent from V)
% eq_h_map(h,rc) = dU/dh = dUmbot_dhi + dUmtop_dhi  (size: [lh-1, lr-1]),
% on the midpoints basically
% solve eq_h (dh) for h(rc)
eq_h_map = dUmbot_dhi + dUmtop_dhi +dUmres_dhi;
h_eq = NaN(1, lr-1); % h_eq at each rc midpoint rci(k) (solution for each rci)

for k = 1:lr-1 % loop through rci values
    Eh = eq_h_map(:, k); % equation eq_h at rci(k)
    H = hi; % midpoint vector of h (stays the same)

    ok = isfinite(Eh); % check
    if nnz(ok) < 2, continue; end % check (very unfortunate case)
    Eh = Eh(ok); H = H(ok); % check

    % look for a zero (change in sign);
    idx = find(Eh(1:end-1).*Eh(2:end) <= 0, 1, 'first'); % find the change in sign for the function
    if ~isempty(idx)
        h_eq(k) = fzero(@(h) interp1(H, Eh, h, 'pchip'), [H(idx), H(idx+1)]); % solution!
    else % fallback --> take the minimum value of the equation (we're looking for a zero)
        [~,km] = min(abs(Eh)); h_eq(k) = H(km);
    end
end

% smooth the curve
h_eq = smoothdata(h_eq, 'movmean', 3);

%% 2) For each V, solve eq_rc(h_eq(rc), rc, V) = 0. Solve for rc*, then h*(rc*) using h_eq(rc*)
V_min = 1e3; V_max = 8e3;
V_list = linspace(V_min, V_max, round(V_max/V_min)*10);  lv = numel(V_list);

h_star    = NaN(lv,1);
rc_star_pure = NaN(lv,1);
rc_star   = NaN(lv,1);
alpha_star= NaN(lv,1);

% warm-start (helps in case of multiple zeros, not necessary if domains are sensible and problem is well-posed)
rc_prev = median(rci(isfinite(h_eq)));

fallback_eq_rc = NaN(1, lr-1); % check if we need fallback or if the rc_eq is solvable

for iv = 1:lv % for each value of V
    Vt = V_list(iv); % "target", or selected value of V

    % build residue eq_rc(rc) = [dU/dr_c](h_eq(rc),rc) - 0.5*V^2 [dC/dr_c](h_eq(rc),rc)
    Rrc = NaN(1, lr-1); % residue
    for k = 1:lr-1 % for each rc, we have a solution h_eq(rc) that is compatible with eq_h
        if isnan(h_eq(k)), continue; end % very unfortunate case in which eq_h was not solvable (unlikely if domain is reasonable)
        % interpolate along h (column k) to evaluate dU/dr_c and dC/dr_c at point h=h_eq(k)
        % this way we can extract the value of the derivatives at the h_eq point
        dUdr_at = interp1(hi, dUmbot_drci(:,k) + dUmtop_drci(:,k) + dUmres_drci(:,k), h_eq(k), 'linear', 'extrap'); % elastic energy along h, for the selected rc
        dCdr_at = interp1(hi, dCbot_drci(:,k), h_eq(k), 'linear', 'extrap'); % "electric energy" (withot 0.5*V^2) along h, for the selected rc
        Rrc(k)  = dUdr_at - 0.5*Vt^2 * dCdr_at; % build residue (target is zero, so residue will be dUdr_at -0.5*Vt^2*dCdr_at - 0)
    end
    
    % filter out the Residue array, on which we will look for a zero
    xr = rci; ok = isfinite(Rrc);
    if nnz(ok) < 2, continue; end
    Rrc = Rrc(ok); xr = xr(ok); 

    % Look for a zero, if multiples, choose the one that is closest to
    % warm-start rc_prev (unlikely but for safety)
    sgn = Rrc(1:end-1).*Rrc(2:end); % sign vector
    idx = find(sgn <= 0); % find first change in sign in the residue vector
    if ~isempty(idx) % we found a zero --> refine the domain interpolating between the positive and negative sign points
        [~,pick] = min(abs(xr(idx) - rc_prev));
        k = idx(pick);
        rc_sol = fzero(@(x) interp1(xr, Rrc, x, 'linear','extrap'), [xr(k), xr(k+1)]); % most likely solution
    else
        fallback_eq_rc(iv) = 1; % for debugging
        [~,kk] = min(abs(Rrc)); rc_sol = xr(kk);   % fallback: minimum residual
    end
    % clamp on the domain --> THIS IS POSSIBLY WHAT MAKES THE PLOTS
    % STRANGE! Well not really
    rc_sol_pure = rc_sol;
    rc_sol = min(max(rc_sol, rci(1)), rci(end));

    % solution
    rc_star_pure(iv) = rc_sol_pure;
    rc_star(iv) = rc_sol;
    h_star(iv)  = interp1(rci, h_eq, rc_sol, 'pchip', 'extrap');  % h*(V)=h_eq(rc*)
    alpha_star(iv) = alpha_sol(h_star(iv), rc_star(iv), params);

    rc_prev = rc_star(iv);   % warm-start for the next Vi will be the solution of the previous Vi level
end


%% PLOTS
figure;
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; plot(V_list/1e3, h_star*1e3, '-o','LineWidth',1.2); grid on
xlabel('V [kV]'); ylabel('h^* [mm]'); title('h^*(V) via eliminazione eq_h');
nexttile; plot(V_list/1e3, rc_star*1e3,'-o','LineWidth',1.2); grid on
xlabel('V [kV]'); ylabel('r_c^* [mm]'); title('r_c^*(V)');
nexttile; plot(V_list/1e3, alpha_star*1e3,'-o','LineWidth',1.2); grid on
xlabel('V [kV]'); ylabel('\alpha^* [mm]'); title('\alpha^*(V)');

%%
