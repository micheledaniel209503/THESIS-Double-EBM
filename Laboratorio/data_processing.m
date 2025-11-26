clear all;
close all;
clc;

%% TEST 2
% TEST 2: one dataset is different from the other in terms of VOLTAGE

%% RAW DATA PROCESSING for dataset V = 5 kV
addpath 'C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data';
load("C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data\Test2\13.10.2025_test2(change of voltage)5kv-T=4sec.mat");

rawData = logsOut(1);

potent = rawData{1}; % potentiometer measurement
loadcell = rawData{2}; % loadcell measurement
voltage = rawData{3}; % voltage amplifier measurement
laser = rawData{4}; % laser measurement
potent_mm = rawData{5};

loadcell_gain = 1.52; % [N/V]
laser_gain = 32; % [mm/V]
voltage_gain = 1500; % [V/V] gain of the high voltage amplifier (from low V to high V)

loadcell_N = (loadcell.Values.Data).*loadcell_gain; % [N] measured force
laser_mm = (laser.Values.Data).*laser_gain; % [mm] measurement displacement of EE
voltage_V = ( voltage.Values.Data).*voltage_gain; % [V] voltage at EBM

%% PLOTS - RAW DATA
figure
plot(potent.Values.Time, potent.Values.Data)
xlabel('time [s]')
ylabel('Voltage [V]')
title('Potentiometer readings')

figure
plot(potent_mm.Values.Time, potent_mm.Values.Data)
xlabel('time [s]')
ylabel('Stroke [mm]')
title('Potentiometer readings')

figure
plot(loadcell.Values.Time, loadcell.Values.Data)
xlabel('time [s]')
ylabel('Output voltage [V]')
title('Loadcell readings')

figure
plot(voltage.Values.Time, voltage.Values.Data)
xlabel('time [s]')
ylabel('EBM voltage [V]')
title('Voltage readings')

figure
plot(laser.Values.Time, laser.Values.Data)
xlabel('time [s]')
ylabel('Output voltage [V]')
title('Laser readings')

figure
plot(laser.Values.Time, laser_mm)
xlabel('time [s]')
ylabel('x  [mm]')
title('Laser measurement: EE displacement')

figure
plot(loadcell.Values.Time, loadcell_N)
xlabel('time [s]')
ylabel('F  [N]')
title('Force at EE')

figure
plot(voltage.Values.Time, voltage_V.*1e-3)
xlabel('time [s]')
ylabel('EBM voltage [kV]')
title('Voltage at EBM poles')

%% STROKE - data processing
close all;

time = loadcell.Values.Time;
time_start = 29; % [s] time at which EE starts moving upwards
idx_start = find(time == time_start);
time_end = 152; % [s] time at which voltage drops to zero
idx_end = find(time == time_end);
laser_offset_mm = laser_mm(idx_start); % [mm] distance of EE when it starts moving
laser_mm_nooffset = laser_mm - laser_offset_mm; % [mm] de-offset laser displacement

figure
subplot(2, 1, 1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
xlabel('time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')
subplot(2, 1, 2)
plot(time(idx_start:idx_end), voltage.Values.Data(idx_start:idx_end));
xlabel('time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

% TO FIND STROKE: max - min at every switch in V, but offsetted by the
% trend
% REMOVE THE TREND: interpolate the minima, interpolate the maxima
%Tsquare = 4; % [s] wave period
Tsquare = 2; % [s] half

figure
subplot(2, 1, 1)
plot(time(idx_start:40e4), laser_mm(idx_start:40e4) - laser_offset_mm)
xlabel('time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')
subplot(2, 1, 2)
plot(time(idx_start:40e4), voltage.Values.Data(idx_start:40e4));
xlabel('time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

%% find minima and maxima
time_stamps_up = time_start:Tsquare:time_end; % time values from 29 to time_end with a step of 4 seconds
time_stamps_down = time_start+1:Tsquare:time_end+1; % minima

idx_up = find(ismember(time, time_stamps_up));
idx_down = find(ismember(time, time_stamps_down));

% extract maxima and minima
laser_mm_maxima = laser_mm_nooffset(idx_up); % at UP (black lines)
laser_mm_minima = laser_mm_nooffset(idx_down); % at DOWN (red lines)

% throw away the first maxima and the last minima
laser_mm_maxima = laser_mm_maxima(2:end); % remove the first maxima
laser_mm_minima = laser_mm_minima(1:end-1); % remove the last minima
idx_up = idx_up(2:end); % update indices for maxima
idx_down = idx_down(1:end-1); % update indices for minima

% compute stroke
% stroke indexes will be between the up and down indexes (it's just an
% assumption)
stroke = laser_mm_maxima - laser_mm_minima;
stroke_V5 = stroke;
idx_stroke = (idx_up + idx_down) / 2;

figure
subplot(3, 1, 1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), laser_mm_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), laser_mm_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')

subplot(3, 1, 2)
plot(time(idx_down), stroke, 'bo', 'MarkerFaceColor', 'b');
grid on
xlabel('Time [s]')
ylabel('Stroke [mm]')
title('Stroke estimate')

subplot(3, 1, 3)
plot(time(idx_start:idx_end), voltage.Values.Data(idx_start:idx_end));
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
end
xlabel('Time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

% plot the stroke points
figure
plot((1:numel(stroke)), stroke, 'ko', 'MarkerFaceColor', 'k')

% the up and down trend is due to the different polarity that induces the
% different peaks: the higher stroke points correspond to NEGATIVE polarity

%% FORCE - data processing
% the weight of the EBM + wire + plastic connector should be removed from
% the loadcell readings
loadcell_offset_V = -0.34; % [V] voltage corresponding to the weight of EBM + wire + plastic connector (pin)
loadcell_nooffset = loadcell.Values.Data - loadcell_offset_V; % [V] debiased loadcell readings

loadcell_nooffset_N = loadcell_nooffset.*loadcell_gain; % [N] force acting on the loadcell, debiased = force on the EE
loadcell_nooffset_N_V5 = loadcell_nooffset_N;

figure
plot(time, loadcell_nooffset_N)
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE (EBM)')

% experiment domain (start = when EE starts moving, end = when EE reaches
% maximum height)
figure
plot(time(idx_start:idx_end), loadcell_nooffset_N(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
end
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE (EBM)')


%% STROKE vs FORCE plot
% !!!!!!!!! NOTE: idx_stroke is to be chosen correctly !!!!!!!!!! 
% force at red is the representative one --> when EBM is DOWN (acutated)

figure 
plot(loadcell_nooffset_N(idx_down), stroke, 'bo', 'MarkerFaceColor', 'b'); % this one should be the correct one
grid on
xlabel('Force [N]')
ylabel('Stroke [mm]')
title('Stroke vs Force')

%% FORCE vs DISPLACEMENT curves
% Building F(h) curves

% the minima of force correspond to curves @V = 0
% the maxima of displacement (idx_up) correspond to curves @V = 0

% the maxima of force correspond to curves @V = 5 kV
% the minima of displacement (idx_down) correspond to curves @V = 5 kV

% extract force minima for the force
loadcell_N_V05 = loadcell_nooffset_N(idx_up);
% extract force maxima for the force
loadcell_N_V5 = loadcell_nooffset_N(idx_down);

figure
subplot(2,1,1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
hold on
for i = 1:length(idx_up)
    plot(time(idx_up(i)), laser_mm_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    plot(time(idx_down(i)), laser_mm_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
grid on
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')

subplot(2,1,2)
plot(time(idx_start:idx_end), loadcell_nooffset_N(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    plot(time(idx_up(i)), loadcell_N_V05(i),  'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    plot(time(idx_down(i)), loadcell_N_V5(i),  'ro', 'MarkerFaceColor', 'r');
end
grid on
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE with offset correction (EBM)')

% plot F(d) where d = displacement, not the real height of the system

figure
plot(laser_mm_maxima, loadcell_N_V05, 'DisplayName', 'V = 0 kV');% F(d) @V=0 V
hold on
plot(laser_mm_minima, loadcell_N_V5, 'DisplayName', 'V = 5 kV'); % F(d) @V=5 kV
xlabel('d [mm]')
ylabel('F [N]')
legend('show')
grid on
grid minor
title('Force vs displacement at different voltages')

% weird pattern due to higher FORCE and higher STROKE when V<0 wrt when V>0

%%







% ----------------------------------------------------







close all;
clc;

%% RAW DATA PROCESSING for dataset V = 4 kV 
addpath 'C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data';
load("C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data\Test2\13.10.2025_test2(changeofvoltage)4kv-T=4sec.mat");

rawData = logsOut(1);

potent = rawData{1}; % potentiometer measurement
loadcell = rawData{2}; % loadcell measurement
voltage = rawData{3}; % voltage amplifier measurement
laser = rawData{4}; % laser measurement
potent_mm = rawData{5};

loadcell_gain = 1.52; % [N/V]
laser_gain = 32; % [mm/V]
voltage_gain = 1500; % [V/V] gain of the high voltage amplifier (from low V to high V)

loadcell_N = (loadcell.Values.Data).*loadcell_gain; % [N] measured force
laser_mm = (laser.Values.Data).*laser_gain; % [mm] measurement displacement of EE
voltage_V = ( voltage.Values.Data).*voltage_gain; % [V] voltage at EBM

%% PLOTS - RAW DATA
figure
plot(potent.Values.Time, potent.Values.Data)
xlabel('time [s]')
ylabel('Voltage [V]')
title('Potentiometer readings')

figure
plot(potent_mm.Values.Time, potent_mm.Values.Data)
xlabel('time [s]')
ylabel('Stroke [mm]')
title('Potentiometer readings')

figure
plot(loadcell.Values.Time, loadcell.Values.Data)
xlabel('time [s]')
ylabel('Output voltage [V]')
title('Loadcell readings')

figure
plot(voltage.Values.Time, voltage.Values.Data)
xlabel('time [s]')
ylabel('EBM voltage [V]')
title('Voltage readings')

figure
plot(laser.Values.Time, laser.Values.Data)
xlabel('time [s]')
ylabel('Output voltage [V]')
title('Laser readings')

figure
plot(laser.Values.Time, laser_mm)
xlabel('time [s]')
ylabel('x  [mm]')
title('Laser measurement: EE displacement')

figure
plot(loadcell.Values.Time, loadcell_N)
xlabel('time [s]')
ylabel('F  [N]')
title('Force at EE')

figure
plot(voltage.Values.Time, voltage_V.*1e-3)
xlabel('time [s]')
ylabel('EBM voltage [kV]')
title('Voltage at EBM poles')


%% STROKE - data processing

time = loadcell.Values.Time;
time_start = 29; % [s] time at which EE starts moving upwards
idx_start = find(time == time_start);
time_end = 152; % [s] time at which voltage drops to zero
idx_end = find(time == time_end);
laser_offset_mm = laser_mm(idx_start); % [mm] distance of EE when it starts moving
laser_mm_nooffset = laser_mm - laser_offset_mm; % [mm] de-offset laser displacement

figure
subplot(2, 1, 1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
xlabel('time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')
subplot(2, 1, 2)
plot(time(idx_start:idx_end), voltage.Values.Data(idx_start:idx_end));
xlabel('time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

% TO FIND STROKE: max - min at every switch in V, but offsetted by the
% trend
% REMOVE THE TREND: interpolate the minima, interpolate the maxima
%Tsquare = 4; % [s] wave period
Tsquare = 2; % [s] half

figure
subplot(2, 1, 1)
plot(time(idx_start:40e4), laser_mm(idx_start:40e4) - laser_offset_mm)
xlabel('time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')
subplot(2, 1, 2)
plot(time(idx_start:40e4), voltage.Values.Data(idx_start:40e4));
xlabel('time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

%% find minima and maxima
time_stamps_up = time_start:Tsquare:time_end; % time values from 29 to time_end with a step of 4 seconds
time_stamps_down = time_start+1:Tsquare:time_end+1; % minima

idx_up = find(ismember(time, time_stamps_up));
idx_down = find(ismember(time, time_stamps_down));

% extract maxima and minima
laser_mm_maxima = laser_mm_nooffset(idx_up); % at UP (black lines)
laser_mm_minima = laser_mm_nooffset(idx_down); % at DOWN (red lines)

% throw away the first maxima and the last minima
laser_mm_maxima = laser_mm_maxima(2:end); % remove the first maxima
laser_mm_minima = laser_mm_minima(1:end-1); % remove the last minima
idx_up = idx_up(2:end); % update indices for maxima
idx_down = idx_down(1:end-1); % update indices for minima

% compute stroke
% stroke indexes will be between the up and down indexes (it's just an
% assumption)
stroke = laser_mm_maxima - laser_mm_minima;
stroke_V4 = stroke;
idx_stroke = (idx_up + idx_down) / 2;

figure
subplot(3, 1, 1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), laser_mm_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), laser_mm_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')

subplot(3, 1, 2)
plot(time(idx_down), stroke, 'bo', 'MarkerFaceColor', 'b');
grid on
xlabel('Time [s]')
ylabel('Stroke [mm]')
title('Stroke estimate')

subplot(3, 1, 3)
plot(time(idx_start:idx_end), voltage.Values.Data(idx_start:idx_end));
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
end
xlabel('Time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

% plot the stroke points
figure
plot((1:numel(stroke)), stroke, 'ko', 'MarkerFaceColor', 'k')

%% FORCE - data processing
% the weight of the EBM + wire + plastic connector should be removed from
% the loadcell readings
loadcell_offset_V = -0.34; % [V] voltage corresponding to the weight of EBM + wire + plastic connector (pin)
loadcell_nooffset = loadcell.Values.Data - loadcell_offset_V; % [V] debiased loadcell readings

loadcell_nooffset_N = loadcell_nooffset.*loadcell_gain; % [N] force acting on the loadcell, debiased = force on the EE
loadcell_nooffset_N_V4 = loadcell_nooffset_N;

figure
plot(time, loadcell_nooffset_N)
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE (EBM)')

% experiment domain (start = when EE starts moving, end = when EE reaches
% maximum height)
figure
plot(time(idx_start:idx_end), loadcell_nooffset_N(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
end
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE (EBM)')

%% STROKE vs FORCE plot
% !!!!!!!!! NOTE: idx_stroke is to be chosen correctly !!!!!!!!!! 
% force at red is the representative one --> when EBM is DOWN (acutated)

figure 
plot(loadcell_nooffset_N(idx_down), stroke_V4, 'bo', 'MarkerFaceColor', 'b'); % this one should be the correct one
grid on
xlabel('Force [N]')
ylabel('Stroke [mm]')
title('Stroke vs Force')

%% FORCE vs DISPLACEMENT curves
% Building F(h) curves

% the minima of force correspond to curves @V = 0
% the maxima of displacement (idx_up) correspond to curves @V = 0

% the maxima of force correspond to curves @V = 4 kV
% the minima of displacement (idx_down) correspond to curves @V = 4 kV

% extract force minima for the force
loadcell_N_V04 = loadcell_nooffset_N(idx_up);
% extract force maxima for the force
loadcell_N_V4 = loadcell_nooffset_N(idx_down);

figure
subplot(2,1,1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
hold on
for i = 1:length(idx_up)
    plot(time(idx_up(i)), laser_mm_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    plot(time(idx_down(i)), laser_mm_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
grid on
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')

subplot(2,1,2)
plot(time(idx_start:idx_end), loadcell_nooffset_N(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    plot(time(idx_up(i)), loadcell_N_V04(i),  'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    plot(time(idx_down(i)), loadcell_N_V4(i),  'ro', 'MarkerFaceColor', 'r');
end
grid on
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE with offset correction (EBM)')

% plot F(d) where d = displacement, not the real height of the system

figure
plot(laser_mm_maxima, loadcell_N_V04, 'DisplayName', 'V = 0 kV');% F(d) @V=0 V
hold on
plot(laser_mm_minima, loadcell_N_V4, 'DisplayName', 'V = 5 kV'); % F(d) @V=4 kV
xlabel('d [mm]')
ylabel('F [N]')
legend('show')
grid on
grid minor
title('Force vs displacement at different voltages')

% weird pattern due to higher FORCE and higher STROKE when V<0 wrt when V>0

%%







% ----------------------------------------------------







close all;
clc;

%% RAW DATA PROCESSING for dataset V = 4 kV 
addpath 'C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data';
load("C:\Users\paola\Desktop\Tesi\Laboratorio\Raw data\Test2\13.10.2025_test2(changeofvoltage)6kv-T=4sec.mat");

rawData = logsOut(1);

potent = rawData{1}; % potentiometer measurement
loadcell = rawData{2}; % loadcell measurement
voltage = rawData{3}; % voltage amplifier measurement
laser = rawData{4}; % laser measurement
potent_mm = rawData{5};

loadcell_gain = 1.52; % [N/V]
laser_gain = 32; % [mm/V]
voltage_gain = 1500; % [V/V] gain of the high voltage amplifier (from low V to high V)

loadcell_N = (loadcell.Values.Data).*loadcell_gain; % [N] measured force
laser_mm = (laser.Values.Data).*laser_gain; % [mm] measurement displacement of EE
voltage_V = ( voltage.Values.Data).*voltage_gain; % [V] voltage at EBM

%% PLOTS - RAW DATA
figure
plot(potent.Values.Time, potent.Values.Data)
xlabel('time [s]')
ylabel('Voltage [V]')
title('Potentiometer readings')

figure
plot(potent_mm.Values.Time, potent_mm.Values.Data)
xlabel('time [s]')
ylabel('Stroke [mm]')
title('Potentiometer readings')

figure
plot(loadcell.Values.Time, loadcell.Values.Data)
xlabel('time [s]')
ylabel('Output voltage [V]')
title('Loadcell readings')

figure
plot(voltage.Values.Time, voltage.Values.Data)
xlabel('time [s]')
ylabel('EBM voltage [V]')
title('Voltage readings')

figure
plot(laser.Values.Time, laser.Values.Data)
xlabel('time [s]')
ylabel('Output voltage [V]')
title('Laser readings')

figure
plot(laser.Values.Time, laser_mm)
xlabel('time [s]')
ylabel('x  [mm]')
title('Laser measurement: EE displacement')

figure
plot(loadcell.Values.Time, loadcell_N)
xlabel('time [s]')
ylabel('F  [N]')
title('Force at EE')

figure
plot(voltage.Values.Time, voltage_V.*1e-3)
xlabel('time [s]')
ylabel('EBM voltage [kV]')
title('Voltage at EBM poles')


%% STROKE - data processing

time = loadcell.Values.Time;
time_start = 29; % [s] time at which EE starts moving upwards
idx_start = find(time == time_start);
time_end = 152; % [s] time at which voltage drops to zero
idx_end = find(time == time_end);
laser_offset_mm = laser_mm(280000); % [mm] distance of EE when it starts moving
laser_mm_nooffset = laser_mm - laser_offset_mm; % [mm] de-offset laser displacement

figure
subplot(2, 1, 1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
xlabel('time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')
subplot(2, 1, 2)
plot(time(idx_start:idx_end), voltage.Values.Data(idx_start:idx_end));
xlabel('time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

% TO FIND STROKE: max - min at every switch in V, but offsetted by the
% trend
% REMOVE THE TREND: interpolate the minima, interpolate the maxima
%Tsquare = 4; % [s] wave period
Tsquare = 2; % [s] half

figure
subplot(2, 1, 1)
plot(time(idx_start:40e4), laser_mm(idx_start:40e4) - laser_offset_mm)
xlabel('time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')
subplot(2, 1, 2)
plot(time(idx_start:40e4), voltage.Values.Data(idx_start:40e4));
xlabel('time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

%% find minima and maxima
time_stamps_up = time_start:Tsquare:time_end; % time values from 29 to time_end with a step of 4 seconds
time_stamps_down = time_start+1:Tsquare:time_end+1; % minima

idx_up = find(ismember(time, time_stamps_up));
idx_down = find(ismember(time, time_stamps_down));

% extract maxima and minima
laser_mm_maxima = laser_mm_nooffset(idx_up); % at UP (black lines)
laser_mm_minima = laser_mm_nooffset(idx_down); % at DOWN (red lines)

% throw away the first maxima and the last minima
laser_mm_maxima = laser_mm_maxima(2:end); % remove the first maxima
laser_mm_minima = laser_mm_minima(1:end-1); % remove the last minima
idx_up = idx_up(2:end); % update indices for maxima
idx_down = idx_down(1:end-1); % update indices for minima

% compute stroke
% stroke indexes will be between the up and down indexes (it's just an
% assumption)
stroke = laser_mm_maxima - laser_mm_minima;
stroke_V6 = stroke;
idx_stroke = (idx_up + idx_down) / 2;

figure
subplot(3, 1, 1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
    plot(time(idx_up(i)), laser_mm_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
    plot(time(idx_down(i)), laser_mm_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')

subplot(3, 1, 2)
plot(time(idx_down), stroke, 'bo', 'MarkerFaceColor', 'b');
grid on
xlabel('Time [s]')
ylabel('Stroke [mm]')
title('Stroke estimate')

subplot(3, 1, 3)
plot(time(idx_start:idx_end), voltage.Values.Data(idx_start:idx_end));
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
end
xlabel('Time [s]');
ylabel('Voltage [V]')
title('EBM voltage')

% plot the stroke points
figure
plot((1:numel(stroke)), stroke, 'ko', 'MarkerFaceColor', 'k')

%% FORCE - data processing
% the weight of the EBM + wire + plastic connector should be removed from
% the loadcell readings
loadcell_offset_V = -0.34; % [V] voltage corresponding to the weight of EBM + wire + plastic connector (pin)
loadcell_nooffset = loadcell.Values.Data - loadcell_offset_V; % [V] debiased loadcell readings

loadcell_nooffset_N = loadcell_nooffset.*loadcell_gain; % [N] force acting on the loadcell, debiased = force on the EE
loadcell_nooffset_N_V6 = loadcell_nooffset_N;

figure
plot(time, loadcell_nooffset_N)
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE (EBM)')

% experiment domain (start = when EE starts moving, end = when EE reaches
% maximum height)
figure
plot(time(idx_start:idx_end), loadcell_nooffset_N(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    xline(time(idx_up(i)), 'k--');
end
for i = 1:length(idx_down)
    xline(time(idx_down(i)), 'r--');
end
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE (EBM)')

%% STROKE vs FORCE plot
% !!!!!!!!! NOTE: idx_stroke is to be chosen correctly !!!!!!!!!! 
% force at red is the representative one --> when EBM is DOWN (acutated)

figure 
plot(loadcell_nooffset_N(idx_down), stroke_V6, 'bo', 'MarkerFaceColor', 'b'); % this one should be the correct one
grid on
xlabel('Force [N]')
ylabel('Stroke [mm]')
title('Stroke vs Force')

%% FORCE vs DISPLACEMENT curves
% Building F(h) curves

% the minima of force correspond to curves @V = 0
% the maxima of displacement (idx_up) correspond to curves @V = 0

% the maxima of force correspond to curves @V = 6 kV
% the minima of displacement (idx_down) correspond to curves @V = 6 kV

% extract force minima for the force
loadcell_N_V06 = loadcell_nooffset_N(idx_up);
% extract force maxima for the force
loadcell_N_V6 = loadcell_nooffset_N(idx_down);

figure
subplot(2,1,1)
plot(time(idx_start:idx_end), laser_mm(idx_start:idx_end) - laser_offset_mm)
hold on
for i = 1:length(idx_up)
    plot(time(idx_up(i)), laser_mm_maxima(i), 'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    plot(time(idx_down(i)), laser_mm_minima(i), 'ro', 'MarkerFaceColor', 'r');
end
grid on
xlabel('Time [s]')
ylabel('Displacement [mm]')
title('Laser measurement: EE displacement with offset correction')

subplot(2,1,2)
plot(time(idx_start:idx_end), loadcell_nooffset_N(idx_start:idx_end))
hold on
for i = 1:length(idx_up)
    plot(time(idx_up(i)), loadcell_N_V06(i),  'ko', 'MarkerFaceColor', 'k');
end
for i = 1:length(idx_down)
    plot(time(idx_down(i)), loadcell_N_V6(i),  'ro', 'MarkerFaceColor', 'r');
end
grid on
xlabel('Time [s]')
ylabel('Force [N]')
title('Force on the EE with offset correction (EBM)')

% plot F(d) where d = displacement, not the real height of the system

figure
plot(laser_mm_maxima, loadcell_N_V06, 'DisplayName', 'V = 0 kV');% F(d) @V=0 V
hold on
plot(laser_mm_minima, loadcell_N_V6, 'DisplayName', 'V = 6 kV'); % F(d) @V=6 kV
xlabel('d [mm]')
ylabel('F [N]')
legend('show')
grid on
grid minor
title('Force vs displacement at different voltages')

% weird pattern due to higher FORCE and higher STROKE when V<0 wrt when V>0


%% COMPARISON
close all;

figure
plot(laser_mm_maxima, loadcell_N_V04, 'DisplayName', 'V = 0 kV test 4'); % F(d) @V=0 V
hold on
plot(laser_mm_maxima, loadcell_N_V05, 'DisplayName', 'V = 0 kV test 5'); % F(d) @V=0 V
plot(laser_mm_maxima, loadcell_N_V06, 'DisplayName', 'V = 0 kV test 6'); % F(d) @V=0 V
plot(laser_mm_minima, loadcell_N_V4, 'DisplayName', 'V = 4 kV'); % F(d) @V=4 kV
plot(laser_mm_minima, loadcell_N_V5, 'DisplayName', 'V = 5 kV'); % F(d) @V=5 kV
plot(laser_mm_minima, loadcell_N_V6, 'DisplayName', 'V = 6 kV'); % F(d) @V=6 kV

xlabel('d [mm]')
ylabel('F [N]')
legend('show')
grid on
grid minor
title('Force vs displacement at different voltages')


figure 
subplot(3,1,1)
plot(loadcell_nooffset_N_V4(idx_down), stroke_V4, 'bo', 'MarkerFaceColor', 'b');
grid on
xlabel('Force [N]')
ylabel('Stroke [mm]')
xlim([0 2.5]);
ylim([0 6]);
title('Stroke vs Force @ V = 4 kV')

subplot(3,1,2)
plot(loadcell_nooffset_N_V5(idx_down), stroke_V5, 'bo', 'MarkerFaceColor', 'b');
grid on
xlabel('Force [N]')
ylabel('Stroke [mm]')
xlim([0 2.5]);
ylim([0 6]);
title('Stroke vs Force @ V = 5 kV')

subplot(3,1,3)
plot(loadcell_nooffset_N_V6(idx_down), stroke_V6, 'bo', 'MarkerFaceColor', 'b');
grid on
xlabel('Force [N]')
ylabel('Stroke [mm]')
xlim([0 2.5]);
ylim([0 6]);
title('Stroke vs Force @ V = 6 kV')