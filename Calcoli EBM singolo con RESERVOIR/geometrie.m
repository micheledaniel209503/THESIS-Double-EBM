%% TEST ALPHA_SOL
clear all;
close all;
clc;

ri = 8e-3;
ro = 15e-3;
r1r = 8e-3; % = ri
r2r = 15e-3; % = ro
hr = 10e-3; % 1 cm

params.r1r = r1r;
params.r2r = r2r;
params.hr = hr;
params.ri = ri;
params.ro = ro;

alpha_sol(0 + 2e-3, ro - 2e-3, params);

% Valutazione alpha(h,rc)

h_vec = linspace(0, 1.5e-3, 100)';
rc_vec = linspace(ro, ro - 7e-3, 100)';

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
h_samples = [0, 0.5e-3, 1e-3, 1.5e-3]; % valori campione di h
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

