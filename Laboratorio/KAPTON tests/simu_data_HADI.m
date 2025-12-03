clc
clear all
%% USING HADI's DATA
load('21.11.2025_3stackEBM_0KV_T=4s(repeating3).mat')
Time = logsOut{1}.Values.Time; % [s]
Potent = logsOut{5}.Values.Data; % [mm]
Potent_gain = 27.425;

Loadcell_gain = 1.52; %N/volt
Loadcell = (logsOut{2}.Values.Data);

laser_gain = 32;
Laser = logsOut{4}.Values.Data;

figure
plot(Time, Loadcell)
xlabel('time [s]')
ylabel('force [N]')
title('RAW DATA')
%%
off_Loadcell = mean(Loadcell(Time<5));
off_laser = mean(Laser(Time<5));
off_Potent = mean(Potent(Time<5));

Potent = (Potent-off_Potent) * Potent_gain;
Loadcell = (Loadcell - off_Loadcell) * Loadcell_gain;
Laser = (Laser - off_laser) * laser_gain;
Loadcell = smooth(Loadcell,5000);
Laser = smooth(Laser,5000);

% SINGLE EBM: divide by 3 the displacement
Laser = Laser./3;


figure(1)
plot(Time, Loadcell)
xlabel('time [s]')
ylabel('force [N]')

figure(2)
plot(Time, Laser)
xlabel('time [s]')
ylabel('displacement [mm]')

figure(3)
plot(Time, Potent)
xlabel('time [s]')
ylabel('displacement [mm]')

figure(4)
plot(Laser, Loadcell)
xlabel('displacement [mm]')
ylabel('force [N]')

%% processing
t_start = 18;
t_end = 144;
idx_start = find(Time == t_start);
idx_end = find(Time == t_end);

Loadcell = Loadcell(idx_start:idx_end);
Laser = Laser(idx_start:idx_end);
%Potent = Potent(idx_start:idx_end);

figure(5)
plot(Laser, Loadcell)
xlabel('displacement [mm]')
ylabel('force [N]')
title('DATA - force vs displacement')

%% KAPTON TEST for characterization: find Y young modulus
% PARAMETERS
par.nu = 0.30;
par.ri = 8e-3;               % [m]
par.ro = 15e-3;              % [m]
par.t0 = 25.4e-6;
%Y_test = 1.25e9;
Y_test = 1e9;
h_test = linspace(0, 2.5e-3, 100);
F_test = kapton_force_model(h_test, Y_test,par);

figure
plot(h_test*1e3, F_test);
xlabel('displacement [mm]')
ylabel('force [N]')

figure(13)
plot(h_test*1e3, F_test)
hold on
plot(Laser, Loadcell -0.1)
xlabel('displacement [mm]')
ylabel('force [N]')
legend('Model', 'Data')
%% FITTING
% force computed by the model given x, mu and the parameters
model_fun = @(mu, x) kapton_force_model(x, mu, par);
% initial estimate of mu (by eye)
Y0 = Y_test;    % [Pa]

% constraints on mu
lb = 2.5e6; % lower bound
ub = 2.5e9; % upper bound

% lsqcurvefit will minimize ||F(x_all, mu, params) - F_all||^2 on mu,
% starting from mu0
Y_est = lsqcurvefit(model_fun, Y0, Laser*1e-3, Loadcell-0.1, lb, ub);

F_model = kapton_force_model(Laser*1e-3, Y_est, par);

figure
plot(Laser, F_model)
hold on
plot(Laser, Loadcell-0.1)
xlabel('displacement [mm]')
ylabel('force [N]')
legend('fit', 'data')

%% OBSERVATIONS
% Using 3 ebms and then dividing by 3 the displacement gives similar
% results: the Young modulus is still roughly 0.8 GPa, even a bit lower
% than the single ebm: this may be due to even more VHB put in series,
% which adds even more compliancy.