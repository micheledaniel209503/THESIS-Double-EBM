clear all; 
%close all;
T_sample = 1e-4; % [s] sampling period
slope = 0.5; %
max_pos = 28; % [mm] max stroke of EE after bias
T = 4;
T_v = 4;
T_per = 1/slope;
Potent_gain = 27.425;
Loadcell_gain = 1.52;

v_slope = 0.502; % sec
v_0 = 5; % Kv
Kp = 15; % 20;
Load_pos = 2; % mm
Ki = 1; %1
N_start = 1;
N_finish = ceil(max_pos/(slope*T)) + 3;

Duration_time = (N_finish + 1)*T;%second
omega_filt = 2*pi*100;
Voltage = 5000;

num_periods = 5;
bias = 80;
run_time = 20; %start_time+(3/slope*max_pos)*num_periods+stop_time;

%% DATA ANALYSIS - SECOND VHB DISC - 16mm
load('run3_data_logsout_26_11_15mm_SECOND_VHB.mat');
%logsOut.plot
% offset
% loadcell starts measuring an increase in force at t = 21.8 and stops
% increasing at t = 51.8
time_start = 21.8;
time_end = 44.0;
timeData = logsOut{1}.Values.Time;
idx_start = find(timeData == time_start);
idx_end = find(timeData == time_end);
time = timeData(idx_start:idx_end);

% Extract relevant data from logsOut starting from the identified index
forceData = logsOut{3}.Values.Data; % [V]
forceData = forceData*Loadcell_gain; % [N]
dispData = logsOut{4}.Values.Data; % [mm]
force = forceData(idx_start:idx_end); % [N]
disp = dispData(idx_start:idx_end); % [mm]

% compute offsets
force_offset = mean(forceData(1:idx_start));
force = force - force_offset;
force16 = (-1)*force;
disp_offset = dispData(idx_start);
disp16 = disp - disp_offset;

figure
subplot(2,1,1)
plot(time, force16)
xlim([time_start, time_end]);
xlabel('time [s]')
ylabel('force [N]')
subplot(2,1,2)
plot(time, disp16)
xlim([time_start, time_end]);
xlabel('time [s]')
ylabel('displacement [mm]')

figure
plot(disp16, force16)
xlabel('displacement [mm]')
ylabel('force [N]')

%% DATA ANALYSIS - SECOND VHB DISC - 20mm
load('run4_data_logsout_26_11_20mm_SECOND_VHB.mat');
%logsOut.plot
% offset
% loadcell starts measuring an increase in force at t = 21.8 and stops
% increasing at t = 51.8
time_start = 21.8;
time_end = 51.78;
timeData = logsOut{1}.Values.Time;
idx_start = find(timeData == time_start);
idx_end = find(timeData == time_end);
time = timeData(idx_start:idx_end);

% Extract relevant data from logsOut starting from the identified index
forceData = logsOut{3}.Values.Data; % [V]
forceData = forceData*Loadcell_gain; % [N]
dispData = logsOut{4}.Values.Data; % [mm]
force = forceData(idx_start:idx_end); % [N]
disp = dispData(idx_start:idx_end); % [mm]

% compute offsets
force_offset = mean(forceData(1:idx_start));
force = force - force_offset;
force20 = (-1)*force;
disp_offset = dispData(idx_start);
disp20 = disp - disp_offset;

figure
subplot(2,1,1)
plot(time, force20)
xlim([time_start, time_end]);
xlabel('time [s]')
ylabel('force [N]')
subplot(2,1,2)
plot(time, disp20)
xlim([time_start, time_end]);
xlabel('time [s]')
ylabel('displacement [mm]')

figure
plot(disp20, force20)
xlabel('displacement [mm]')
ylabel('force [N]')

%% DATA ANALYSIS - SECOND VHB DISC - 28mm
load('run5_data_logsout_26_11_28mm_SECOND_VHB.mat');
%logsOut.plot

% offset
% loadcell starts measuring an increase in force at t = 21.8 and stops
% increasing at t = 51.8
time_start = 21.8;
time_end = 68.0;
timeData = logsOut{1}.Values.Time;
idx_start = find(timeData == time_start);
idx_end = find(timeData == time_end);
time = timeData(idx_start:idx_end);

% Extract relevant data from logsOut starting from the identified index
forceData = logsOut{3}.Values.Data; % [V]
forceData = forceData*Loadcell_gain; % [N]
dispData = logsOut{4}.Values.Data; % [mm]
force = forceData(idx_start:idx_end); % [N]
disp = dispData(idx_start:idx_end); % [mm]

% compute offsets
force_offset = mean(forceData(1:idx_start));
force = force - force_offset;
force28 = (-1)*force;
disp_offset = dispData(idx_start);
disp28 = disp - disp_offset;

figure
subplot(2,1,1)
plot(time, force28)
xlim([time_start, time_end]);
xlabel('time [s]')
ylabel('force [N]')
subplot(2,1,2)
plot(time, disp28)
xlim([time_start, time_end]);
xlabel('time [s]')
ylabel('displacement [mm]')

figure
plot(disp28, force28)
xlabel('displacement [mm]')
ylabel('force [N]')

%% FINAL PLOTS
figure
plot(disp16, force16)
hold on
plot(disp20, force20)
plot(disp28, force28)
xlabel('displacement [mm]')
ylabel('force [N]')
xlim([min(disp28) max(disp28)])
ylim([min(force28) max(force28)])
grid on
grid minor
legend('run 3', 'run 4', 'run 5', 'Location','best')
title('Force vs Displacement for pre-stretched VHB')

%% FITTING
% prepare data
% conversion in meters
x16 = disp16 * 1e-3; % [m]
x20 = disp20 * 1e-3; % [m]
x28 = disp28 * 1e-3; % [m]
F16 = force16; % [N]
F20 = force20; % [N]
F28 = force28; % [N]

% arrays of data of the three experiments
x_all = [x16; x20; x28];
F_all = [F16; F20; F28];

% parameters of the model
par.R0       = 0.5*55e-3; % [m]
par.R1       = 0.5*75e-3; % [m]
par.lambda_p = par.R1 / par.R0;
par.r1r      = 0.5*15e-3; % [m]
par.r2r      = par.R0; % [m]
par.t0       = 0.5e-3; % [m]

% force computed by the model given x, mu and the parameters
model_fun = @(mu, x) vhb_force_model(x, mu, par);
% initial estimate of mu (by eye)
mu0 = 0.4e5;    % [Pa]

% constraints on mu
lb = 0; % lower bound
ub = 1e6; % 1 MPa

% lsqcurvefit will minimize ||F(x_all, mu, params) - F_all||^2 on mu,
% starting from mu0
mu_est = lsqcurvefit(model_fun, mu0, x_all, F_all, lb, ub);

%% CHECK FITTING RESULTS
% compute F with the model using the estimated mu, on the experimental
% points
F16_fit = vhb_force_model(x16, mu_est, par);
F20_fit = vhb_force_model(x20, mu_est, par);
F28_fit = vhb_force_model(x28, mu_est, par);

figure; hold on; grid on;
plot(x16*1e3, F16, 'b.',  'DisplayName','run3 - data');
plot(x16*1e3, F16_fit, 'b-', 'DisplayName','run3 - fit');

plot(x20*1e3, F20, 'r.',  'DisplayName','run4 - data');
plot(x20*1e3, F20_fit, 'r-', 'DisplayName','run4 - fit');

plot(x28*1e3, F28, 'k.',  'DisplayName','run5 - data');
plot(x28*1e3, F28_fit, 'k-', 'DisplayName','run5 - fit');

xlabel('displacement [mm]');
ylabel('force [N]');
legend('Location','best');
title(sprintf('Neo-Hooke fit with \\mu = %.3e Pa', mu_est));
