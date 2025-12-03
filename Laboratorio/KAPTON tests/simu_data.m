clear all;
close all;

%% DATA ANALYSIS 
% pull an empty single ebm --> characterize the Kapton using the model
% (simplified model, only 1 strain along meridian).
% NOTE: run3 with bias 130 and stroke 30 is the best one. The others are
% not very good.
load("run3_KAPTON_bias130_stroke30.mat")
Time = logsOut{1}.Values.Time;
Loadcell_gain = 1.52; %N/volt
Laser = (logsOut{2}.Values.Data) ;
Potent = (logsOut{1}.Values.Data);
off_Potent = mean(Potent(Time<5));

Loadcell = (logsOut{3}.Values.Data);
off_Loadcell = mean(Loadcell(Time<5));

off_laser = mean(Laser(Time<5));
laser_gain = 32;

Potent_gain = 27.425;

Potent = (Potent-off_Potent) * Potent_gain;
Loadcell = (Loadcell - off_Loadcell) * Loadcell_gain;
Laser = (Laser - off_laser) * laser_gain;
Loadcell = smooth(Loadcell,5000);
Laser = smooth(Laser,5000);



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
t_start = 30.5;
t_end = 71.8;
idx_start = find(Time == t_start);
idx_end = find(Time == t_end);

Loadcell = Loadcell(idx_start:idx_end);
Laser = Laser(idx_start:idx_end);
%Potent = Potent(idx_start:idx_end);

figure(5)
plot(Laser, Loadcell)
xlabel('displacement [mm]')
ylabel('force [N]')

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

%% FITTING
% force computed by the model given x, mu and the parameters
model_fun = @(mu, x) kapton_force_model(x, mu, par);
% initial estimate of mu (by eye)
Y0 = Y_test;    % [Pa]

% constraints on mu
lb = 0; % lower bound
ub = 2.5e9; % upper bound

% lsqcurvefit will minimize ||F(x_all, mu, params) - F_all||^2 on mu,
% starting from mu0
Y_est = lsqcurvefit(model_fun, Y0, Laser*1e-3, Loadcell - 0.1, lb, ub);

F_model = kapton_force_model(Laser*1e-3, Y_est, par);

figure
plot(Laser, F_model)
hold on
plot(Laser, Loadcell - 0.1)
xlabel('displacement [mm]')
ylabel('force [N]')
legend('fit', 'data')

%% OBSERVATIONS
% the resulting value of Young modulus is a bit lower than expected --> may
% due to the presence of VHB in series that acts as a softer spring and
% lowers the overall stiffness.