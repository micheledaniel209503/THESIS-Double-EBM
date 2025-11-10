%% LEAP – MAPS U(h,rc) e C(h,rc) + F–z curves at fixed V
clear; close all; clc

%% PARAMETERS
Y = 1.25e9;              % [Pa]
nu = 0.30;
epsP = 3.9*8.854e-12;    % [F/m] polymer
epsO = 2.7*8.854e-12;    % [F/m] oil
ri = 8e-3;               % [m]
ro = 15e-3;              % [m]
r1r = ri;
r2r = ro;
hr = 7.5e-3;               % [m] height of the reservoir (internal cylindrical support)
t0 = 25.4e-6;            % [m]
t0_res = 0.5e-3;         % [m] thickness of unstretched VHB disc (reservoir)
to_res = 5e-6/2;         % [m] residual oil film (half)
mu = 4e4;  % 40 KPa??          % [-]
R0 = 50e-3/2;            % [m] initial radius of the unstretched VHB disc

params.Y = Y; params.nu = nu; params.epsP = epsP; params.epsO = epsO;
params.ri = ri; params.ro = ro; params.r1r = r1r; params.r2r = r2r; params.hr = hr; params.t0 = t0; params.t0_res = t0_res;
params.to_res = to_res; params.mu = mu; params.R0 = R0;

% create samples for the maps
lr = 151; lh = 201;
rc_vec = linspace(1.03*ri, 0.999*ro, lr);
h_vec  = linspace(0, 1.5e-3, lh).'; % half of the total height (total max h will be 3.0 mm)

%% BUILD MAPS

% compute terms for HALF the double cone
Um_ebm = zeros(lh, lr); Um_res = zeros(lh, lr); C_ebm = zeros(lh, lr);
for j = 1:lr
    [Um_ebm(:,j), Um_res(:,j), C_ebm(:,j)] = single_EBM_reservoir(h_vec, rc_vec(j), params);
end

% compute terms for ENTIRE double cone
% double cone
Um_ebm = 2*Um_ebm;                 % total energy of ebm
Um_res = Um_res;
C_ebm = 0.5*C_ebm;                 % equivalent capacitance
z = 2*h_vec;               % total displacement

% midpoint grids
rci = 0.5*(rc_vec(1:end-1) + rc_vec(2:end));
hi  = 0.5*(h_vec(1:end-1)  + h_vec(2:end));
zi  = 2*hi;                              % doppio cono

% d/drc at midpoint nodes in rc (fixed rows, diff on columns) 
dUebm_drc = zeros(lh, lr-1);
dUres_drc = zeros(lh, lr-1);
dCebm_drc = zeros(lh, lr-1);

for j = 1:lh  % derivatives wrt rc
    dUebm_drc(j,:) = diff(Um_ebm(j,:))./diff(rc_vec);
    dUres_drc(j,:) = diff(Um_res(j,:))./diff(rc_vec);
    dCebm_drc(j,:) = diff(C_ebm(j,:))./diff(rc_vec);
end

% Interpolate in h -> hi (midpoint in h) to align d/drc with (hi, rci)
dUebm_dri = zeros(lh-1, lr-1);
dUres_dri = zeros(lh-1, lr-1);
dCebm_dri = zeros(lh-1, lr-1);

for i = 1:lr-1
    dUebm_dri(:,i) = interp1(h_vec, dUebm_drc(:,i), hi, 'linear');
    dUres_dri(:,i) = interp1(h_vec, dUres_drc(:,i), hi, 'linear');
    dCebm_dri(:,i) = interp1(h_vec, dCebm_drc(:,i), hi, 'linear');
end

% d/dh at midpoint nodes in h (fixed columns, diff on rows)
dUebm_dh = zeros(lh-1, lr);
dUres_dh = zeros(lh-1, lr);
dCebm_dh = zeros(lh-1, lr);

for i = 1:lr % derivatives wrt h
    dUebm_dh(:,i) = diff(Um_ebm(:,i))./diff(z);
    dUres_dh(:,i) = diff(Um_res(:,i))./diff(z);
    dCebm_dh(:,i) = diff(C_ebm(:,i))./diff(z);
end

% Interpolate in rc -> rci (midpoint in rc) to align d/dh with (hi, rci)
dUebm_dhi = zeros(lh-1, lr-1);
dUres_dhi = zeros(lh-1, lr-1);
dCebm_dhi = zeros(lh-1, lr-1);

for j = 1:lh-1
    dUebm_dhi(j,:) = interp1(rc_vec, dUebm_dh(j,:), rci, 'linear');
    dUres_dhi(j,:) = interp1(rc_vec, dUres_dh(j,:), rci, 'linear');
    dCebm_dhi(j,:) = interp1(rc_vec, dCebm_dh(j,:), rci, 'linear');
end

% V^2/2 on midpoint grid (hi, rci)
tol = 1e-12;
mask = abs(dCebm_dri) > tol;
Vsq_half       = NaN(lh-1, lr-1); % midpoint grid
Vsq_half(mask) = (dUebm_dri(mask) + dUres_dri(mask))./dCebm_dri(mask);
Vi = sqrt( max(2*Vsq_half, 0) );               % Vi = V(h_i, rci_k) MAP
Fi = (dUebm_dhi + dUres_dhi) - Vsq_half.*dCebm_dhi;                  % Fi = F(h_i, rci_k) MAP

%% SOLUTION of the equations
V_min = 1; V_max = 8e3;
V_list = linspace(V_min, V_max, 5);  lv = numel(V_list);

% initialize solution vectors
rc_v = NaN(lh-1, lv);
alpha_v = NaN(lh-1, lv);
F_v  = NaN(lh-1, lv);

for iv = 1:lv
    Vt = V_list(iv);
    for j = 1:lh-1 % for each row (each h)
        Vrow = Vi(j,:);  xr = rci; % extract V vector and rc vector for the specific h

        % extract and clean the row
        ok = isfinite(Vrow);
        if nnz(ok) < 2, continue; end
        Vrow = Vrow(ok);  xr = xr(ok);
        rres = Vrow - Vt; % vector of residuals on which we will look for a zero (root)

        % look for root means look for a change in sign (go across zero)
        k = find(rres(1:end-1).*rres(2:end) <= 0, 1, 'first');
        if ~isempty(k) % if we found a root, refine the interval in between the two points k and k+1 through interpolation
            a = xr(k); b = xr(k+1);
            % find a root on the interpolating function
            rc_sol = fzero(@(x) interp1(xr, Vrow, x, 'linear','extrap') - Vt, [a,b]); % V(h,rc*) - Vtarget = 0 --> rc*
        else
            % fallback: closest point (very unlucky case in which there's
            % no root --> root is on the minimum of the residuals
            [~,kk] = min(abs(rres)); rc_sol = xr(kk);
        end
        
        % for this value of h, we have the solution rc* --> compute F
        rc_v(j,iv) = min(max(rc_sol, rci(1)), rci(end)); % clamp (in case rc_sol is out of the grid)
        alpha_v(j,iv) = alpha_sol(h_vec(j),rc_v(j,iv), params);
        F_v(j,iv)  = interp1(rci, Fi(j,:), rc_v(j,iv), 'linear','extrap');
    end
end

%% PLOTS
sty='-';
leg = arrayfun(@(v) sprintf('V = %.1f kV', v/1e3), V_list, 'UniformOutput', false);
figure(3); hold on; grid on

for iv=1:lv
    plot(1e3*z(1:end-1), F_v(:,iv), sty,'LineWidth',1.3, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([0 14]);
end
xlabel('\itz\rm (mm)'); ylabel('\itF\rm (N)')
legend(leg,'Location','best'); 
title('Force–stroke')

figure(5); hold on; grid on
for iv=1:lv
    plot(1e3*zi, 1e3*rc_v(:,iv), sty, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([8 15]);
end
xlabel('\itz\rm (mm)'); ylabel('\itr_c\rm (mm)')
legend(leg, 'Location', 'best');
title('r_c at equilibrium')

figure(6); hold on; grid on
for iv=1:lv
    plot(1e3*zi, 1e3*alpha_v(:,iv), sty, 'Color',[iv/lv 0 0])
    xlim([0 3]);
end
xlabel('\itz\rm (mm)'); ylabel('\italpha\rm (mm)')
legend(leg, 'Location', 'best');
title('alpha at equilibrium')

%% ANIMATION

%% PLOT OF THE CROSS SECTION

rc = ro - 1e-3;
h = 1.5e-3;              % [m] height of ebm
alpha = alpha_sol(h, rc, params);

% Checks
assert(ri>0 && ro>ri, 'Rule: 0 < ri < ro.');
rc = min(max(rc, ri), ro);  % clamp rc in [ri, ro]

% EBM
P0D = [ri,  h];
P0L = [-ri, h];
P1D = [ro, 0];
P1L = [-ro, 0];
P3D = [rc, 0];
P3L = [-rc, 0];
P2D = [ri, -h];
P2L = [-ri, -h];

% RESERVOIR
P1RD = [r2r, h];
P1RL = [-r2r, h];
P2RD = [r1r, h+hr];
P2RL = [-r1r, h+hr];

% figure
figure('Color','w','Name','SYSTEM CROSS SECTION','NumberTitle','off');
ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');

% function to plot a segment
plotSeg = @(A,B,sty,lw) plot(ax,[ A(1)  B(1)]*1e3,[ A(2)  B(2)]*1e3,sty,'LineWidth',lw);

% EBM segments
plotSeg(P0L, P0D, 'b', 2.0);
plotSeg(P0D, P3D, 'b', 2.0);
plotSeg(P3D, P1D, 'b', 2.0);
plotSeg(P3D, P2D, 'b', 2.0);
plotSeg(P2D, P2L, 'b', 2.0);
plotSeg(P2L, P3L, 'b', 2.0);
plotSeg(P3L, P1L, 'b', 2.0);
plotSeg(P3L, P0L, 'b', 2.0);


% RESERVOIR segments and arcs
plotSeg(P0L, P1RL, 'b', 2.0);
plotSeg(P0D, P1RD, 'b', 2.0);
plotSeg(P2RL, P2RD, 'b', 2.0);

[xarc, yarc] = arc_two_pts(P1RD, P2RD, alpha, 10, 1);
plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
[xarc, yarc] = arc_two_pts(P1RL, P2RL, alpha, 10, 1);
plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);

% axes limits
%xlim(ax,[-ro ro]*1e3);
%ylim(ax,[-h h+hr]*1e3*2.5);

%% ANIMATION

Vt     = 8e3;
[~,iv] = min(abs(V_list - Vt));                    
rc_col = rc_v(:,iv);
h_col  = h_vec(1:end-1);
alpha_col = alpha_v(:,iv);

% figure
figure('Color','w','Name','SYSTEM cross section','NumberTitle','off');
ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');

% function to plot a segment
plotSeg = @(A,B,sty,lw) plot(ax,[ A(1)  B(1)]*1e3,[ A(2)  B(2)]*1e3,sty,'LineWidth',lw);

for i = 1:length(h_col)
    
    h  = h_col(i);
    rc = rc_col(i);
    alpha  = alpha_col(i);

    % clean axes at each frame iteration
    cla(ax);
    grid(ax,'on'); box(ax,'on'); axis(ax,'equal');

    % EBM
    P0D = [ri,  h];
    P0L = [-ri, h];
    P1D = [ro, 0];
    P1L = [-ro, 0];
    P3D = [rc, 0];
    P3L = [-rc, 0];
    P2D = [ri, -h];
    P2L = [-ri, -h];
    
    % RESERVOIR
    P1RD = [r2r, h];
    P1RL = [-r2r, h];
    P2RD = [r1r, h+hr];
    P2RL = [-r1r, h+hr];

    % EBM segments
    plotSeg(P0L, P0D, 'b', 2.0);
    plotSeg(P0D, P3D, 'b', 2.0);
    plotSeg(P3D, P1D, 'b', 2.0);
    plotSeg(P3D, P2D, 'b', 2.0);
    plotSeg(P2D, P2L, 'b', 2.0);
    plotSeg(P2L, P3L, 'b', 2.0);
    plotSeg(P3L, P1L, 'b', 2.0);
    plotSeg(P3L, P0L, 'b', 2.0);
    
    
    % RESERVOIR segments and arcs
    plotSeg(P0L, P1RL, 'b', 2.0);
    plotSeg(P0D, P1RD, 'b', 2.0);
    plotSeg(P2RL, P2RD, 'b', 2.0);
    
    [xarc, yarc] = arc_two_pts(P1RD, P2RD, alpha, 10, 1);
    plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);
    [xarc, yarc] = arc_two_pts(P1RL, P2RL, alpha, 10, 1);
    plot(ax, xarc*1e3, yarc*1e3, 'b', 'LineWidth', 2.0);


    drawnow limitrate;   % refresh
end