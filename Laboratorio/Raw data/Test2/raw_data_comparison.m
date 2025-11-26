clear all;
close all;
clc;

%% TEST 2
% TEST 2: one dataset is different from the other in terms of VOLTAGE

%% RAW DATA PROCESSING for dataset V = 5 kV
addpath 'C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data';
load("C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data\Test2\13.10.2025_test2(change of voltage)5kv-T=4sec.mat");

rawData_5 = logsOut(1);

potent_5 = rawData_5{1}; % potentiometer measurement
loadcell_5 = rawData_5{2}; % loadcell measurement
voltage_5 = rawData_5{3}; % voltage amplifier measurement
laser_5 = rawData_5{4}; % laser measurement
potent_mm_5 = rawData_5{5};

loadcell_gain = 1.52; % [N/V]
laser_gain = 32; % [mm/V]
voltage_gain = 1500; % [V/V] gain of the high voltage amplifier (from low V to high V)

loadcell_N_5 = (loadcell_5.Values.Data).*loadcell_gain; % [N] measured force
laser_mm_5 = (laser_5.Values.Data).*laser_gain; % [mm] measurement displacement of EE
voltage_V_5 = ( voltage_5.Values.Data).*voltage_gain; % [V] voltage at EBM

%% RAW DATA PROCESSING for dataset V = 4 kV
addpath 'C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data';
load("C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data\Test2\13.10.2025_test2(changeofvoltage)4kv-T=4sec.mat");

rawData_4 = logsOut(1);

potent_4 = rawData_4{1}; % potentiometer measurement
loadcell_4 = rawData_4{2}; % loadcell measurement
voltage_4 = rawData_4{3}; % voltage amplifier measurement
laser_4 = rawData_4{4}; % laser measurement
potent_mm_4 = rawData_4{5};

loadcell_N_4 = (loadcell_4.Values.Data).*loadcell_gain; % [N] measured force
laser_mm_4 = (laser_4.Values.Data).*laser_gain; % [mm] measurement displacement of EE
voltage_V_4 = ( voltage_4.Values.Data).*voltage_gain; % [V] voltage at EBM

%% RAW DATA PROCESSING for dataset V = 6 kV
addpath 'C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data';
load("C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data\Test2\13.10.2025_test2(changeofvoltage)6kv-T=4sec.mat");

rawData_6 = logsOut(1);

potent_6 = rawData_6{1}; % potentiometer measurement
loadcell_6 = rawData_6{2}; % loadcell measurement
voltage_6 = rawData_6{3}; % voltage amplifier measurement
laser_6 = rawData_6{4}; % laser measurement
potent_mm_6 = rawData_6{5};

loadcell_N_6 = (loadcell_6.Values.Data).*loadcell_gain; % [N] measured force
laser_mm_6 = (laser_6.Values.Data).*laser_gain; % [mm] measurement displacement of EE
voltage_V_6 = ( voltage_6.Values.Data).*voltage_gain; % [V] voltage at EBM


%% PLOTS

time = loadcell_4.Values.Time;
time_start = 29; % [s] time at which EE starts moving upwards
idx_start = find(time == time_start);
time_end = 152; % [s] time at which voltage drops to zero
idx_end = find(time == time_end);

figure
subplot(3,1,1)
plot(loadcell_4.Values.Time, loadcell_N_4, 'DisplayName', '4 kV');
hold on
plot(loadcell_5.Values.Time, loadcell_N_5, 'DisplayName', '5 kV');
plot(loadcell_6.Values.Time, loadcell_N_6, 'DisplayName', '6 kV');
xline(time_start, 'k--');
xline(time_end, 'k--');
legend('show')
xlabel('time [s]')
ylabel('force [N]')

subplot(3,1,2)
plot(laser_4.Values.Time, laser_mm_4, 'DisplayName', '4 kV');
hold on
plot(laser_5.Values.Time, laser_mm_5, 'DisplayName', '5 kV');
plot(laser_6.Values.Time, laser_mm_6, 'DisplayName', '6 kV');
xline(time_start, 'k--');
xline(time_end, 'k--');
legend('show')
xlabel('time [s]')
ylabel('displacement [mm]')

subplot(3,1,3)
plot(voltage_4.Values.Time, voltage_V_4./1e3, 'DisplayName', '4 kV');
hold on
plot(voltage_5.Values.Time, voltage_V_5./1e3, 'DisplayName', '5 kV');
plot(voltage_6.Values.Time, voltage_V_6./1e3, 'DisplayName', '6 kV');
xline(time_start, 'k--');
xline(time_end, 'k--');
legend('show')
xlabel('time [s]')
ylabel('voltage [kV]')
sgtitle('RAW DATA')

%% FORCE AND DISPLACEMENT OFFSET
% if we consider that the experiment really begins when the EE starts
% moving, then at that point the force should be zero. RIGHT??
% Calculate the force offset based on the loadcell readings at the start time
forceOffset_4 = loadcell_N_4(idx_start);
forceOffset_5 = loadcell_N_5(idx_start);
forceOffset_6 = loadcell_N_6(idx_start);
% there's a 3 grams difference in the offsets between (4,5) and 6

loadcell_N_4_nooffset = loadcell_N_4 - forceOffset_4;
loadcell_N_5_nooffset = loadcell_N_5 - forceOffset_5;
loadcell_N_6_nooffset = loadcell_N_6 - forceOffset_6;

% the laser offset is the mean for idx < idx_start
laser_offset_mm_4 = mean(laser_mm_4(1:idx_start));
laser_offset_mm_5 = mean(laser_mm_5(1:idx_start));
laser_offset_mm_6 = mean(laser_mm_6(1:idx_start));

laser_mm_4_nooffset = laser_mm_4 - laser_offset_mm_4; % [mm] de-offset laser displacement
laser_mm_5_nooffset = laser_mm_5 - laser_offset_mm_5; % [mm] de-offset laser displacement
laser_mm_6_nooffset = laser_mm_6 - laser_offset_mm_6; % [mm] de-offset laser displacement

figure
subplot(2,1,1)
plot(time, laser_mm_4_nooffset, 'DisplayName', '4 kV')
hold on
plot(time, laser_mm_5_nooffset, 'DisplayName', '5 kV')
plot(time, laser_mm_6_nooffset, 'DisplayName', '6 kV')
xline(time_start, 'k--');
xline(time_end, 'k--');
grid on
legend('show')
xlabel('time [s]')
ylabel('displacement [mm]')

subplot(2,1,2)
plot(time, loadcell_N_4_nooffset, 'DisplayName', '4 kV')
hold on
plot(time, loadcell_N_5_nooffset, 'DisplayName', '5 kV')
plot(time, loadcell_N_6_nooffset, 'DisplayName', '6 kV')
xline(time_start, 'k--');
xline(time_end, 'k--');
grid on
legend('show')
xlabel('time [s]')
ylabel('force [N]')
sgtitle('OFFSET CORRECTION')


%% DATA EXTRACTION
Tsquare = 2; % [s]
Tact = 0.5; % [s] half

time_stamps_up = time_start:Tsquare:time_end; % time values from 29 to time_end with a step of 4 seconds
time_stamps_down = time_start+Tact:Tsquare:time_end+Tact; % minima

idx_up = find(ismember(time, time_stamps_up));
idx_down = find(ismember(time, time_stamps_down));

% extract maxima and minima
laser_mm_4_maxima = laser_mm_4_nooffset(idx_up); % at UP (black lines) --> on
laser_mm_4_minima = laser_mm_4_nooffset(idx_down); % at DOWN (red lines) --> off
laser_mm_5_maxima = laser_mm_5_nooffset(idx_up); % at UP (black lines)
laser_mm_5_minima = laser_mm_5_nooffset(idx_down); % at DOWN (red lines)
laser_mm_6_maxima = laser_mm_6_nooffset(idx_up); % at UP (black lines)
laser_mm_6_minima = laser_mm_6_nooffset(idx_down); % at DOWN (red lines)

% throw away the first maxima and the last minima
laser_mm_4_maxima = laser_mm_4_maxima(2:end); % remove the first maxima
laser_mm_4_minima = laser_mm_4_minima(1:end-1); % remove the last minima
laser_mm_5_maxima = laser_mm_5_maxima(2:end); % remove the first maxima
laser_mm_5_minima = laser_mm_5_minima(1:end-1); % remove the last minima
laser_mm_6_maxima = laser_mm_6_maxima(2:end); % remove the first maxima
laser_mm_6_minima = laser_mm_6_minima(1:end-1); % remove the last minima
idx_up = idx_up(2:end); % update indices for maxima
idx_down = idx_down(1:end-1); % update indices for minima

% compute stroke
% stroke indexes will be between the up and down indexes (it's just an
% assumption)
stroke_4 = laser_mm_4_maxima - laser_mm_4_minima;
stroke_5 = laser_mm_5_maxima - laser_mm_5_minima;
stroke_6 = laser_mm_6_maxima - laser_mm_6_minima;
idx_stroke = (idx_up + idx_down) / 2;

loadcell_N_4_nooffset_og = loadcell_N_4_nooffset;
loadcell_N_5_nooffset_og = loadcell_N_5_nooffset;
loadcell_N_6_nooffset_og = loadcell_N_6_nooffset;

% extract force minima for the force
loadcell_N_04_nooffset = loadcell_N_4_nooffset_og(idx_up);
loadcell_N_05_nooffset = loadcell_N_5_nooffset_og(idx_up);
loadcell_N_06_nooffset = loadcell_N_6_nooffset_og(idx_up);

% extract force maxima for the force
loadcell_N_4_nooffset = loadcell_N_4_nooffset_og(idx_down);
loadcell_N_5_nooffset = loadcell_N_5_nooffset_og(idx_down);
loadcell_N_6_nooffset = loadcell_N_6_nooffset_og(idx_down);

%% STROKE vs FORCE plot
% !!!!!!!!! NOTE: idx_stroke is to be chosen correctly !!!!!!!!!! 
% force at red is the representative one --> when EBM is DOWN (acutated)

max_stroke_4 = max(stroke_4);
max_stroke_5 = max(stroke_5);
max_stroke_6 = max(stroke_6);

figure 
subplot(3,1,1)
plot(loadcell_N_4_nooffset, stroke_4, 'bo', 'MarkerFaceColor', 'b');
yline(max_stroke_4, 'k--');
grid on
grid minor
xlabel('Force [N]')
ylabel('Stroke [mm]')
xlim([0 2.5]);
ylim([0 6]);
title('Stroke vs Force @ V = 4 kV')

subplot(3,1,2)
plot(loadcell_N_5_nooffset, stroke_5, 'bo', 'MarkerFaceColor', 'b');
yline(max_stroke_5, 'k--');
grid on
grid minor
xlabel('Force [N]')
ylabel('Stroke [mm]')
xlim([0 2.5]);
ylim([0 6]);
title('Stroke vs Force @ V = 5 kV')

subplot(3,1,3)
plot(loadcell_N_6_nooffset, stroke_6, 'bo', 'MarkerFaceColor', 'b');
yline(max_stroke_6, 'k--');
grid on
grid minor
xlabel('Force [N]')
ylabel('Stroke [mm]')
xlim([0 2.5]);
ylim([0 6]);
title('Stroke vs Force @ V = 6 kV')


figure
plot(laser_mm_4_maxima, loadcell_N_04_nooffset, 'DisplayName', 'V = 0 kV test 4 kV'); % F(d) @V=0 V
hold on
plot(laser_mm_5_maxima, loadcell_N_05_nooffset, 'DisplayName', 'V = 0 kV test 5 kV'); % F(d) @V=0 V
plot(laser_mm_6_maxima, loadcell_N_06_nooffset, 'DisplayName', 'V = 0 kV test 6 kV'); % F(d) @V=0 V
plot(laser_mm_4_minima, loadcell_N_4_nooffset, 'DisplayName', 'V = 4 kV'); % F(d) @V=4 kV
plot(laser_mm_5_minima, loadcell_N_5_nooffset, 'DisplayName', 'V = 5 kV'); % F(d) @V=5 kV
plot(laser_mm_6_minima, loadcell_N_6_nooffset, 'DisplayName', 'V = 6 kV'); % F(d) @V=6 kV
xlabel('d [mm]')
ylabel('F [N]')
legend('show')
legend('Location', 'best');
grid on
grid minor
title('Force vs displacement at different voltages')

%%
figure
plot(laser_mm_4_minima, loadcell_N_04_nooffset, 'DisplayName', 'V = 0 kV test 4 kV'); % F(d) @V=0 V
hold on
plot(laser_mm_4_maxima, loadcell_N_4_nooffset, 'DisplayName', 'V = 4 kV'); % F(d) @V=4 kV
legend('show')
xlabel('d [mm]')
ylabel('F [N]')
grid on

%% Plot for each dataset

% DISPLACEMENTS
figure
subplot(3, 1, 1)
plot(time(idx_start:idx_end), laser_mm_4_nooffset(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), laser_mm_4_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), laser_mm_4_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('4 kV')

subplot(3,1,2)
plot(time(idx_start:idx_end), laser_mm_5_nooffset(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), laser_mm_5_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), laser_mm_5_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('5 kV')

subplot(3,1,3)
plot(time(idx_start:idx_end), laser_mm_6_nooffset(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), laser_mm_6_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), laser_mm_6_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('6 kV')

% FORCES
figure
subplot(3, 1, 1)
plot(time(idx_start:idx_end), loadcell_N_4_nooffset_og(idx_start:idx_end))
hold on
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), loadcell_N_4_nooffset(i), 'ro', 'MarkerFaceColor', 'r');
end
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), loadcell_N_04_nooffset(i), 'ko', 'MarkerFaceColor', 'k');
end
xlabel('Time [s]')
ylabel('Force [mm]')
title('4 kV')

subplot(3, 1, 2)
plot(time(idx_start:idx_end), loadcell_N_5_nooffset_og(idx_start:idx_end))
hold on
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), loadcell_N_5_nooffset(i), 'ro', 'MarkerFaceColor', 'r');
end
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'r--');
    plot(time(idx_up(i)), loadcell_N_05_nooffset(i), 'ko', 'MarkerFaceColor', 'k');
end
xlabel('Time [s]')
ylabel('Force [mm]')
title('5 kV')

subplot(3, 1, 3)
plot(time(idx_start:idx_end), loadcell_N_6_nooffset_og(idx_start:idx_end))
hold on
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), loadcell_N_6_nooffset(i), 'ro', 'MarkerFaceColor', 'r');
end
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), loadcell_N_06_nooffset(i), 'ko', 'MarkerFaceColor', 'k');
end
xlabel('Time [s]')
ylabel('Force [mm]')
title('6 kV')
