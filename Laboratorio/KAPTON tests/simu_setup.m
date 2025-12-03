clear all; 
%close all;
T_sample = 1e-4; % [s] sampling period
slope = 0.5; %

bias = 130; % [mm]
max_pos = 30; % [mm] max stroke of EE after bias

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
run_time = 20; %start_time+(3/slope*max_pos)*num_periods+stop_time;
