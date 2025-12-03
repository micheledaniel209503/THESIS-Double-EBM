%% KAPTON TEST for characterization: find Y young modulus
% PARAMETERS
par.nu = 0.30;
par.ri = 8e-3;               % [m]
par.ro = 15e-3;              % [m]
par.t0 = 25.4e-6;
Y_test = 1.25e9;
h_test = linspace(0, 2e-3, 100);
F_test = kapton_force_model(h_test, Y_test,par);

figure
plot(h_test*1e3, F_test);
xlabel('displacement [mm]')
ylabel('force [N]')