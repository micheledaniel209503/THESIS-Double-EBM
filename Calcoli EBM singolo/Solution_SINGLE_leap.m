%% LEAP – MAPS U(h,rc) e C(h,rc) + F–z curves at fixed V
clear; close all; clc

%% PARAMETERS
Y = 1.25e9;              % [Pa]
nu = 0.30;
epsP = 3.9*8.854e-12;    % [F/m] polymer
epsO = 2.7*8.854e-12;    % [F/m] oil

ri = 8e-3;               % [m]
ro = 15e-3;              % [m]
t0 = 25.4e-6;            % [m]
to_res = 5e-6/2;         % [m] residual oil film (half)

alpha = 0;
fo    = 1.0;

% create samples for the maps
lr = 151; lh = 201;
rc_vec = linspace(1.03*ri, 0.999*ro, lr);
h_vec  = linspace(0, 1.5e-3, lh).'; % half of the total height (total max h will be 3.0 mm)

%% BUILD MAPS

% compute terms for HALF the double cone
Uel = zeros(lh,lr);   C = zeros(lh,lr);   Omega = zeros(lh,lr);
for j = 1:lr
    [Uel(:,j), C(:,j), Omega(:,j)] = LEAP_conical(h_vec, rc_vec(j), ...
        ri, ro, t0, Y, nu, epsP, epsO, to_res);
end

Uel_simp = zeros(lh,lr); C_simp = zeros(lh,lr); Omega_simp = zeros(lh,lr);
for k = 1:lr
    [Uel_simp(:,k), C_simp(:,k), Omega_simp(:,k)] = LEAP_conical_simp(h_vec, rc_vec(k), ...
        ri, ro, t0, Y, nu, epsP, epsO, to_res);
end

% compute terms for ENTIRE double cone
% double cone
U = 2*Uel;                 % total energy
U_simp = 2*Uel_simp;
C = 0.5*C;                 % equivalent capacitance
C_simp = 0.5*C_simp;
z = 2*h_vec;               % total displacement

% midpoint grids
rci = 0.5*(rc_vec(1:end-1) + rc_vec(2:end));
hi  = 0.5*(h_vec(1:end-1)  + h_vec(2:end));
zi  = 2*hi;                              % doppio cono

% d/drc at midpoint nodes in rc (fixed rows, diff on columns) 
dUdr = zeros(lh, lr-1);
dCdr = zeros(lh, lr-1);
dUdr_simp = zeros(lh, lr-1);
dCdr_simp = zeros(lh, lr-1);

for j = 1:lh  % derivatives wrt rc
    dUdr(j,:) = diff(U(j,:))./diff(rc_vec);
    dCdr(j,:) = diff(C(j,:))./diff(rc_vec);

    dUdr_simp(j,:) = diff(U_simp(j,:))./diff(rc_vec);
    dCdr_simp(j,:) = diff(C_simp(j,:))./diff(rc_vec);
end

% Interpolate in h -> hi (midpoint in h) to align d/drc with (hi, rci)
dUdri = zeros(lh-1, lr-1);
dCdri = zeros(lh-1, lr-1);
dUdri_simp = zeros(lh-1, lr-1);
dCdri_simp = zeros(lh-1, lr-1);

for i = 1:lr-1
    dUdri(:,i) = interp1(h_vec, dUdr(:,i), hi, 'linear');
    dCdri(:,i) = interp1(h_vec, dCdr(:,i), hi, 'linear');
    dUdri_simp(:,i) = interp1(h_vec, dUdr_simp(:,i), hi, 'linear');
    dCdri_simp(:,i) = interp1(h_vec, dCdr_simp(:,i), hi, 'linear');
end

% d/dh at midpoint nodes in h (fixed columns, diff on rows)
dUdh = zeros(lh-1, lr);
dCdh = zeros(lh-1, lr);
dUdh_simp = zeros(lh-1, lr);
dCdh_simp = zeros(lh-1, lr);

for i = 1:lr % derivatives wrt h
    dUdh(:,i) = diff(U(:,i))./diff(z);
    dCdh(:,i) = diff(C(:,i))./diff(z);
    dUdh_simp(:,i) = diff(U_simp(:,i))./diff(z);
    dCdh_simp(:,i) = diff(C_simp(:,i))./diff(z);
end

% Interpolate in rc -> rci (midpoint in rc) to align d/dh with (hi, rci)
dUdhi = zeros(lh-1, lr-1);
dCdhi = zeros(lh-1, lr-1);
dUdhi_simp = zeros(lh-1, lr-1);
dCdhi_simp = zeros(lh-1, lr-1);

for j = 1:lh-1
    dUdhi(j,:) = interp1(rc_vec, dUdh(j,:), rci, 'linear');
    dCdhi(j,:) = interp1(rc_vec, dCdh(j,:), rci, 'linear');
    dUdhi_simp(j,:) = interp1(rc_vec, dUdh_simp(j,:), rci, 'linear');
    dCdhi_simp(j,:) = interp1(rc_vec, dCdh_simp(j,:), rci, 'linear');
end

% V^2/2 and F on midpoint grid (hi, rci)
tol = 1e-12;
mask = abs(dCdri) > tol;
Vsq_half       = NaN(lh-1, lr-1); % midpoint grid
Vsq_half(mask) = dUdri(mask)./dCdri(mask);
Vi = sqrt( max(2*Vsq_half, 0) );               % Vi = V(h_i, rci_k) MAP
Fi = dUdhi - Vsq_half.*dCdhi;                  % Fi = F(h_i, rci_k) MAP

mask_simp = abs(dCdri_simp) > tol;
Vsq_half_simp       = NaN(lh-1, lr-1); % midpoint grid
Vsq_half_simp(mask_simp) = dUdri_simp(mask_simp)./dCdri_simp(mask_simp);
Vi_simp = sqrt( max(2*Vsq_half_simp, 0) );               % Vi = V(h_i, rci_k) MAP
Fi_simp = dUdhi_simp - Vsq_half_simp.*dCdhi_simp;                  % Fi = F(h_i, rci_k) MAP

%% CURVES @V=0
F0 = pi*(ro^2 - ri^2)*t0*Y/(2*(1-nu^2)*(ro-ri)^4) * h_vec.^3;
alpha0 = atan2(h_vec, ro-ri);  % zipping angle @V=0
% identical for simplified model

%% CURVES F–z for specific V values
V_list = (2e3:2e3:8e3).';  % [V] list of voltages
lv = numel(V_list);

% initialize solution vectors [200*4] [h*V]
rc_v = NaN(lh-1, lv);
F_v  = NaN(lh-1, lv);
alpha_v = NaN(lh-1, lv);
rc_v_simp = NaN(lh-1, lv);
F_v_simp  = NaN(lh-1, lv);
alpha_v_simp = NaN(lh-1, lv);

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
        F_v(j,iv)  = interp1(rci, Fi(j,:), rc_v(j,iv), 'linear','extrap');
        alpha_v(j,iv) = atan2(hi(j), rc_v(j,iv)-ri);
    end
end


% simplified model
for iv = 1:lv
    Vt = V_list(iv);
    for j = 1:lh-1 % for each row (each h)
        Vrow = Vi_simp(j,:);  xr = rci; % extract V vector and rc vector for the specific h

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
            rc_sol_simp = fzero(@(x) interp1(xr, Vrow, x, 'linear','extrap') - Vt, [a,b]); % V(h,rc*) - Vtarget = 0 --> rc*
        else
            % fallback: closest point (very unlucky case in which there's
            % no root --> root is on the minimum of the residuals
            [~,kk] = min(abs(rres)); rc_sol_simp = xr(kk);
        end
        
        % for this value of h, we have the solution rc* --> compute F
        rc_v_simp(j,iv) = min(max(rc_sol_simp, rci(1)), rci(end)); % clamp (in case rc_sol is out of the grid)
        F_v_simp(j,iv)  = interp1(rci, Fi_simp(j,:), rc_v_simp(j,iv), 'linear','extrap');
        alpha_v_simp(j,iv) = atan2(hi(j), rc_v_simp(j,iv)-ri);
    end
end

% SOLUTION will be:
% F_v will be a vector as F(h,V)
% alpha_v will be a vector as alpha(h,V)

%% PLOTS
% FULL MODEL
sty='-';
figure(3); hold on; grid on
plot(1e3*(fo-alpha)/(1-alpha)*z, (fo+alpha)/(1+alpha)*F0, 'k', 'LineWidth',1.2) % F @V=0
leg = {'\itV\rm=0'};

for iv=1:lv
    plot(1e3*(fo-alpha)/(1-alpha)*z(1:end-1), (fo+alpha)/(1+alpha)*F_v(:,iv), sty,'LineWidth',1.3, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([0 14]);
    leg{end+1} = sprintf('\\itV\\rm=%.0f kV', V_list(iv)/1e3);
end
xlabel('\itz\rm (mm)'); ylabel('\itF\rm (N)')

legend(leg,'Location','best'); title('Force–stroke')

figure(5); hold on; grid on
for iv=1:lv
    plot(1e3*zi, 1e3*rc_v(:,iv), sty, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([8 15]);
end
xlabel('\itz\rm (mm)'); ylabel('\itr_c\rm (mm)')

legend(leg(2:end),'Location','best'); title('r_c at equilibrium')

figure(7); hold on; grid on
plot(1e3*z, alpha0, 'k', 'LineWidth',1.2)
for iv=1:lv
    plot(1e3*z(1:end-1), alpha_v(:,iv), sty, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([0 0.4]);
end
xlabel('z (m)'); ylabel('\alpha (rad)')

legend(leg,'Location','best'); title('\alpha vs z')


% SIMPLIFIED (REDUCED) MODEL
figure(4); hold on; grid on
plot(1e3*(fo-alpha)/(1-alpha)*z, (fo+alpha)/(1+alpha)*F0, 'k', 'LineWidth',1.2) % F @V=0
leg = {'\itV\rm=0'};

for iv=1:lv
    plot(1e3*(fo-alpha)/(1-alpha)*z(1:end-1), (fo+alpha)/(1+alpha)*F_v_simp(:,iv), sty,'LineWidth',1.3, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([0 14]);
    leg{end+1} = sprintf('\\itV\\rm=%.0f kV', V_list(iv)/1e3);
end
xlabel('\itz\rm (mm)'); ylabel('\itF\rm (N)')

legend(leg,'Location','best'); title('Force–stroke REDUCED')

figure(6); hold on; grid on
for iv=1:lv
    plot(1e3*zi, 1e3*rc_v_simp(:,iv), sty, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([8 15]);
end
xlabel('\itz\rm (mm)'); ylabel('\itr_c\rm (mm)')

legend(leg(2:end),'Location','best'); title('r_c at equilibrium REDUCED')

figure(8); hold on; grid on
plot(1e3*z, alpha0, 'k', 'LineWidth',1.2)
for iv=1:lv
    plot(1e3*z(1:end-1), alpha_v_simp(:,iv), sty, 'Color',[iv/lv 0 0])
    xlim([0 3]);
    ylim([0 0.4]);
end
xlabel('z (m)'); ylabel('\alpha (rad)')

legend(leg,'Location','best'); title('\alpha vs z REDUCED')




%% LEAP FUNCTIONS
% computes Uel, C, Omega for the HALF cone
function [Uel,C,Omega] = LEAP_conical(h,rc,ri,ro,t0,E,nu,epsilon_p,epsilon_o,to_left)
% 1/2 LEAP
    uc = -h.^2/(2*(rc-ri)^2)/(1/(rc-ri)-(rc^2+ro^2)/(rc^2-ro^2)/rc);

    r1 = linspace(ri,rc,20)';   r2 = linspace(rc,ro,20)';

    Uel = zeros(size(h));  C = zeros(size(h));  Omega = zeros(size(h));
    for i = 1:numel(h)
        % internal region r1 in [ri, rc]
        eps1_1 =  uc(i)/(rc-ri) + (h(i).^2)/(2*(rc-ri)^2);      % meridian
        eps2_1 =  uc(i)/(rc-ri) * (1 - ri./r1);                 % circumferential
        eps3_1 = -nu/(1-nu) * (eps1_1 + eps2_1);                % normal

        % external region r2 in [rc, ro]
        eps1_2 = -uc(i)*rc/(ro^2-rc^2) * (1 + ro^2./r2.^2);                 % meridian
        eps2_2 =  uc(i)*rc/(ro^2-rc^2) * (ro^2./r2.^2 - 1);                 % circumferential
        eps3_2 = -nu/(1-nu) * (eps1_2 + eps2_2);                            % normal


        Uel(i) = pi*t0*E/(1-nu^2) * (trapz(r1, r1.*(eps1_1.^2 + eps2_1.^2 + 2*nu*eps1_1.*eps2_1)) + trapz(r2, r2.*(eps1_2.^2 + eps2_2.^2 + 2*nu*eps1_2.*eps2_2)) );


        w = h(i)/(rc-ri) * (r1 - ri);

        Cp = 2*pi*epsilon_p/t0 * trapz( r2, (1+eps1_2).*(1+eps2_2)./(1+eps3_2).*r2 ); % polymer capacitance

        if to_left > 0
            Co = 2*pi*epsilon_o/to_left * trapz(r2, (1+eps1_2).*(1+eps2_2).*r2); % oil capacitcance
            C(i) = Cp*Co/(Cp+Co); % oil capacitance is in series with poly
        else
            C(i) = Cp;
        end

        Omega(i) = pi/3*h(i)*(rc^2+ri^2+ri*rc) + pi*(ro^2-rc^2)*to_left;
    end
end

%% SIMPLIFIED MODEL LEAP FUNCTION
% computes Uel, C, Omega for the HALF cone
function [Uel,C,Omega] = LEAP_conical_simp(h,rc,ri,ro,t0,E,nu,epsilon_p,epsilon_o,to_left)
% 1/2 LEAP
    uc = -h.^2/(2*(rc-ri)^2)/(1/(rc-ri)-(rc^2+ro^2)/(rc^2-ro^2)/rc);

    r1 = linspace(ri,rc,20)';   r2 = linspace(rc,ro,20)';

    Uel = zeros(size(h));  C = zeros(size(h));  Omega = zeros(size(h));
    for i = 1:numel(h)
        % internal region r1 in [ri, rc]
        eps1_1 = uc(i)/(rc-ri) + (h(i).^2)/(2*(rc-ri)^2);      % meridian
        eps2_1 =  0;                                           % circumferential
        eps3_1 = 0;                                            % normal

        % external region r2 in [rc, ro]
        eps1_2 = -uc(i)*rc/(ro^2-rc^2) * (1 + ro^2./r2.^2);     % meridian
        eps2_2 = 0;                                             % circumferential
        eps3_2 = 0;                                             % normal

        Uel(i) = pi*t0*E/(1-nu^2) * (trapz(r1, r1.*(eps1_1.^2 + eps2_1.^2 + 2*nu*eps1_1.*eps2_1)) + trapz(r2, r2.*(eps1_2.^2 + eps2_2.^2 + 2*nu*eps1_2.*eps2_2)) );

        w = h(i)/(rc-ri) * (r1 - ri);
       
        Cp = 2*pi*epsilon_p/t0 * trapz( r2, (1+eps1_2).*(1+eps2_2)./(1+eps3_2).*r2 ); % polymer capacitance

        if to_left > 0
            Co = 2*pi*epsilon_o/to_left * trapz(r2, (1+eps1_2).*(1+eps2_2).*r2); % oil capacitance
            C(i) = Cp*Co/(Cp+Co); % oil capacitance is in series with poly
        else
            C(i) = Cp;
        end

        Omega(i) = pi/3*h(i)*(rc^2+ri^2+ri*rc) + pi*(ro^2-rc^2)*to_left;
    end
end
