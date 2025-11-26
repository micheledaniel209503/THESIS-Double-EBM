clear all;
%close all;
clc;

% SYSTEM TOPOLOGY: 2 ebms, stacked one over the other + a reservoir on top
% that bulges (alpha) to accomodate the ebm's change in volume due to h0 !=
% 0 and h>0, rc<ro. The reservoir communicates with both ebms.

ri = 8e-3; % internal radius of ebms
ro = 15e-3; % external radius of ebms
r1r = ri + 1e-3; % internal radius of reservoir
r2r = ro + 3e-3; % external radius of reservoir
hr = 7.5e-3; % reservoir height
h0 = 1.5e-3; % [mm] preload (full height of the single ebm)


%% Valutazione alpha(h,rc)

params.r1r = r1r;
params.r2r = r2r;
params.hr = hr;
params.ro = ro;
params.ri = ri;
params.h0 = h0;

% system at configuration D
params.h0 = 0; % (h=0, rc=ro) AND h0=0
alpha_solD = alpha_sol(0, ro, params)

% system at configuration 0 (h=0, rc=ro)
params.h0 = 1.5e-3;
alpha_sol0 = alpha_sol(0, ro, params)

% system at configuratio 1 (generic h>0, rc<ro)
alpha_sol1 = alpha_sol(1e-3, ro-5e-3, params)
%%
% analisi della soluzione: supponiamo che EE si sposti verso il basso di
% massimo h0 (fino a "fondo corsa")
% supponiamo che rc vari da ro = 15 mm fino a 9 mm
% supponiamo di partire dalla configurazione 0, quindi dopo il montaggio
params.h0 = h0;
h_vec = linspace(0, h0, 100)';
rc_vec = linspace(ro, ro - 6e-3, 100)';

% RISULTATI
% scorri colonna per scorrere h
% scorri riga per scorrere rc
alpha_mat = zeros(100, 100);

for i = 1:length(rc_vec)
    rc_val = rc_vec(i);
    alpha_mat(:, i) = alpha_sol(h_vec, rc_val, params); % compute the column
end

% plot alpha(h) per rc fissato
figure;
rc_samples = [15e-3, 13e-3, 11e-3, 9e-3]; % valori campione di rc
hold on; grid on; box on;
for rc_val = rc_samples
    [~, idx] = min(abs(rc_vec - rc_val));
    plot(h_vec*1e3, alpha_mat(:, idx), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('r_c = %.1f mm', rc_val*1e3));
end
xlabel('h [mm]');
ylabel('alpha(h)');
title('alpha(h) per diversi rc fissati');
legend('Location','best');
hold off;

% plot alpha(rc) per h fissato
figure;
h_samples = [0, 0.5e-3, 1e-3, h0]; % valori campione di h
hold on; grid on; box on;
for h_val = h_samples
    [~, idx] = min(abs(h_vec - h_val));
    plot(rc_vec*1e3, alpha_mat(idx, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('h = %.1f mm', h_val*1e3));
end
xlabel('rc [mm]');
ylabel('alpha(rc)');
title('alpha(rc) per diversi h fissati');
legend('Location','best');
hold off;

