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


%% 1) Solve eq_h(h,rc) for h_eq(rc) (indipendent from V)
% eq_h_map(h,rc) = dU/dh = dUmbot_dhi + dUmtop_dhi  (size: [lh-1, lr-1]),
% on the midpoints basically
% solve eq_h (dh) for h(rc)
eq_h_map = dUmbot_dhi + dUmtop_dhi;
fallback_eq_h = NaN(1, lr-1); % check if we need fallback or if the h_eq is solvable

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
        fallback_eq_h(k) = 1;
        [~,km] = min(abs(Eh)); h_eq(k) = H(km);
    end
end

% points in which we need a fallback:
count_fallbacks_eq_h = sum(fallback_eq_h == 1);
% NOTE: there's a couple of points, at rc --> ro where eq_h is not zero and
% we need the fallback! these points are actually "clamped" to h = 0, in a
% way. The real solution is...?? It should be zero.
% these points are just 2 and are placed near rc->ro

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
        dUdr_at = interp1(hi, dUmbot_drci(:,k) + dUmtop_drci(:,k), h_eq(k), 'linear', 'extrap'); % elastic energy along h, for the selected rc
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


% points in which we need a fallback:
count_fallbacks_eq_rc = sum(fallback_eq_rc == 1);
% NOTE:
% there are not multiple zeros for any of the two equations with these
% domains of h_vec, rc_vec. This is a good thing I would say.

%% PLOTS
figure;
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; plot(V_list/1e3, h_star*1e3, '-o','LineWidth',1.2); grid on
xlabel('V [kV]'); ylabel('h^* [mm]'); title('h^*(V) via eliminazione eq_h');
nexttile; plot(V_list/1e3, rc_star*1e3,'-o','LineWidth',1.2); grid on
xlabel('V [kV]'); ylabel('r_c^* [mm]'); title('r_c^*(V)');
nexttile; plot(V_list/1e3, alpha_star*1e3,'-o','LineWidth',1.2); grid on
xlabel('V [kV]'); ylabel('\alpha^* [mm]'); title('\alpha^*(V)');

%% NOTE
% there's a weird behaviour up until approx V<2.5 kV, in which the solution
% is almost "clamped". This needs to be investigated but it feels like
% there's a mismatch in terms of energies due to something. Maybe an
% incongruence in the Ueltop and Uelbot or something like that.

%% PLOT OF THE CROSS SECTION

rc = ro - 1e-8;
h = 1e-8;              % [m] stroke (downwards) 
alpha = alpha_sol(h, rc, params);
h_bot = h0 - h;
h_top = h0 + h;

% Checks
assert(ri>0 && ro>ri, 'Rule: 0 < ri < ro.');
rc = min(max(rc, ri), ro);  % clamp rc in [ri, ro]
h_bot  = max(h_bot, 0);

% BOTTOM EBM points
P0 = [0,   h_bot];
P1D = [ri,  h_bot];
P2D = [rc,  h_bot/2];
P3D = [ro,  h_bot/2];
P4D = [ri, 0];

P1L = [-ri, h_bot];
P2L = [-rc, h_bot/2];
P3L = [-ro, h_bot/2];
P4L = [-ri, 0];

% TOP EBM points
% check h_bot + h_top = 2*h0
h_tot = h_bot + h_top;
assert(h_tot == 2*h0, 'Incongruence in: h_tot != 2*h0')
P0t = [0, h_tot];
P1Dt = [ri, h_tot];
P3Dt = [ro, h_bot + h_top/2];
P3Lt = [-ro, h_bot + h_top/2];
P1Lt = [-ri, h_tot];


% figure
figure('Color','w','Name','BOTTOM EBM cross section','NumberTitle','off');
ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');

% function to plot a segment
plotSeg = @(A,B,sty,lw) plot(ax,[ A(1)  B(1)]*1e3,[ A(2)  B(2)]*1e3,sty,'LineWidth',lw);

% initial "rest point"
yline(h0*1e3, 'k--');

% BOTTOM EBM segments
plotSeg(P0, P1D, 'b', 2.0);
plotSeg(P1D, P2D, 'b', 2.0);
plotSeg(P2D, P3D, 'b', 2.0);
plotSeg(P2D, P4D, 'b', 2.0);
plotSeg(P4D, P4L, 'b', 2.0);
plotSeg(P4L, P2L, 'b', 2.0);
plotSeg(P2L, P3L, 'b', 2.0);
plotSeg(P2L, P1L, 'b', 2.0);
plotSeg(P1L, P0, 'b', 2.0);

% TOP EBM segments
plotSeg(P1Lt, P1Dt, 'b', 2.0);
[xarc, yarc] = arc_two_pts(P3Dt, P1Dt, alpha, 10, -1);
plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
[xarc, yarc] = arc_two_pts(P3Dt, P1D, alpha, 10, +1);
plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
[xarc, yarc] = arc_two_pts(P1Lt, P3Lt, alpha, 10, -1);
plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
[xarc, yarc] = arc_two_pts(P3Lt, P1L, alpha, 10, +1);
plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);

% axes limits
xlim(ax,[-ro ro]*1e3);
ylim(ax,[0 h0]*1e3*2.5);

%% PLOT CROSS SECTION ALONG THE SOLUTION

rc_star;
h_star;
alpha_star;

rc = rc_star;
h = h_star;              % [m] stroke (downwards) 
alpha = alpha_star;
h_bot = h0 - h;
h_top = h0 + h;

% Checks
assert(ri>0 && ro>ri, 'Rule: 0 < ri < ro.');
rc = min(max(rc, ri), ro);  % clamp rc in [ri, ro]
h_bot  = max(h_bot, 0);

% figure
figure('Color','w','Name','BOTTOM EBM cross section','NumberTitle','off');
ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');
% axes limits
xlim(ax,[-ro ro]*1e3);
ylim(ax,[0 h0]*1e3*2.5);
% initial "rest point"
yline(h0*1e3, 'k--');

% function to plot a segment
plotSeg = @(A,B,sty,lw) plot(ax,[ A(1)  B(1)]*1e3,[ A(2)  B(2)]*1e3,sty,'LineWidth',lw);

for i = 1:length(rc_star)
    % clean axes at each frame iteration
    cla(ax);
    grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    xlim(ax,[-ro ro]*1e3);
    ylim(ax,[0 h0]*1e3*2.5);
    yline(ax, h0*1e3, 'k--');

    % BOTTOM EBM points
    P0 = [0,   h_bot(i)];
    P1D = [ri,  h_bot(i)];
    P2D = [rc(i),  h_bot(i)/2];
    P3D = [ro,  h_bot(i)/2];
    P4D = [ri, 0];
    
    P1L = [-ri, h_bot(i)];
    P2L = [-rc(i), h_bot(i)/2];
    P3L = [-ro, h_bot(i)/2];
    P4L = [-ri, 0];
    
    % TOP EBM points
    % check h_bot + h_top = 2*h0
    h_tot(i) = h_bot(i) + h_top(i);
    assert(h_tot(i) == 2*h0, 'Incongruence in: h_tot != 2*h0')
    P0t = [0, h_tot(i)];
    P1Dt = [ri, h_tot(i)];
    P3Dt = [ro, h_bot(i) + h_top(i)/2];
    P3Lt = [-ro, h_bot(i) + h_top(i)/2];
    P1Lt = [-ri, h_tot(i)];
    
    % BOTTOM EBM segments
    plotSeg(P0, P1D, 'b', 2.0);
    plotSeg(P1D, P2D, 'b', 2.0);
    plotSeg(P2D, P3D, 'b', 2.0);
    plotSeg(P2D, P4D, 'b', 2.0);
    plotSeg(P4D, P4L, 'b', 2.0);
    plotSeg(P4L, P2L, 'b', 2.0);
    plotSeg(P2L, P3L, 'b', 2.0);
    plotSeg(P2L, P1L, 'b', 2.0);
    plotSeg(P1L, P0, 'b', 2.0);
    
    % TOP EBM segments
    plotSeg(P1Lt, P1Dt, 'b', 2.0);
    [xarc, yarc] = arc_two_pts(P3Dt, P1Dt, alpha(i), 10, -1);
    plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
    [xarc, yarc] = arc_two_pts(P3Dt, P1D, alpha(i), 10, +1);
    plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
    [xarc, yarc] = arc_two_pts(P1Lt, P3Lt, alpha(i), 10, -1);
    plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
    [xarc, yarc] = arc_two_pts(P3Lt, P1L, alpha(i), 10, +1);
    plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);

    title(ax, sprintf('Progress %.0f %%   rc = %.2f mm,  h = %.2f mm,  \\alpha = %.1f mm, V = %.1f kV', ...
        i/numel(rc)*100, rc(i)*1e3, h(i)*1e3, alpha(i)*1e3, V_list(i)*1e-3));

    drawnow limitrate;   % refresh

end

%% CHECK: plot elastic potential energy in the domain, to analyze the path of the solution
Um_tot = Um_topt + Um_bott;   % [J] total elastic potential energy

[RC, H] = meshgrid(rc_vec, h_vec);

% find minimum
[Um_min, idx_min] = min(Um_tot(:));
[ih_min, ir_min] = ind2sub(size(Um_tot), idx_min);
h_min  = h_vec(ih_min);
rc_min = rc_vec(ir_min);

figure('Color','w','Name','U_{m,tot}(h,r_c)','NumberTitle','off');

% 3D plot
subplot(1,2,1);
surf(RC*1e3, H*1e3, Um_tot, 'EdgeColor','none'); 
hold on;
plot3(rc_min*1e3, h_min*1e3, Um_min, 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
hold off;
xlabel('r_c [mm]'); ylabel('h [mm]'); zlabel('U_{m,tot} [J]');
title('Potential surface');
colorbar; grid on; view(35,30);

% contour plot
subplot(1,2,2);
contourf(RC*1e3, H*1e3, Um_tot, 30, 'LineColor','none');
hold on;
plot(rc_min*1e3, h_min*1e3, 'ko', 'MarkerFaceColor','w', 'MarkerSize',6);
hold off;
xlabel('r_c [mm]'); ylabel('h [mm]');
title(sprintf('Contour (minimum: h = %.3f mm, r_c = %.3f mm, U = %.4g J)', ...
    h_min*1e3, rc_min*1e3, Um_min));
axis tight; grid on; colorbar;

% OVERLAY: path of the solution (h*, r_c*)
h_sol  = h_star;
rc_sol = rc_star;
V_sol  = V_list;

% evaluate the elastic potential energy on the points of the solution
F_Umtot = griddedInterpolant({h_vec, rc_vec}, Um_tot, 'linear', 'nearest');
U_sol   = F_Umtot(h_sol, rc_sol);

% SURFACE PLOT
subplot(1,2,1); hold on;
% path
plot3(rc_sol*1e3, h_sol*1e3, U_sol, 'k-', 'LineWidth', 1.6, 'DisplayName','(h^*, r_c^*)');
% points
plot3(rc_sol*1e3, h_sol*1e3, U_sol, 'wo', 'MarkerFaceColor','k', 'MarkerSize', 5, 'HandleVisibility','off');

% labels for the voltage
if numel(V_sol) >= 3
    idx_lab = unique(round(linspace(1, numel(V_sol), 3))); % min, mean, max
    for ii = idx_lab
        text(rc_sol(ii)*1e3, h_sol(ii)*1e3, U_sol(ii), ...
            sprintf('  V=%.1f kV', V_sol(ii)*1e-3), ...
            'Color','w','FontSize',9,'FontWeight','bold');
    end
end
hold off;

% CONTOUR PLOT
subplot(1,2,2); hold on;
plot(rc_sol*1e3, h_sol*1e3, 'k-', 'LineWidth', 1.6, 'DisplayName','(h^*, r_c^*)');
plot(rc_sol*1e3, h_sol*1e3, 'wo', 'MarkerFaceColor','k', 'MarkerSize', 5, 'HandleVisibility','off');

% labels for the voltage
if numel(V_sol) >= 3
    idx_lab = unique(round(linspace(1, numel(V_sol), 3)));
    for ii = idx_lab
        text(rc_sol(ii)*1e3, h_sol(ii)*1e3, ...
            sprintf('  V=%.1f kV', V_sol(ii)*1e-3), ...
            'Color','w','FontSize',9,'FontWeight','bold');
    end
end
hold off;
