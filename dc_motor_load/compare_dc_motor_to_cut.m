%% Compare Actual cut data to dc motor data
clear all; close all; clc;
startup;
fs = 8000;                              % sampling rate [Hz]
l2l_flag = false;                       % PM motor = true, induction = false
P = 2;                                  % number of pole pairs, 2 for heavy duty
encoder_count = 1800;                   % encoder count
mill = 'heavy duty';                    % 'heavy duty' or 'lighty duty'

useTrueSpeed = false;                   % PM motor = true, induction = false


%% Preprocess data 
filename = 'oct24_46Hz';
replicate_current = true;
data_46Hz = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill,replicate_current);

filename = 'oct22_60Hz';
data_60 = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

filename = 'oct23_1800_30Hz';
data_29 = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

filename = 'oct23_660_30Hz';
data_660 = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

filename = 'oct23_open';
data_open = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

filename = 'oct23_1800_2_5Hz';
data_2_5 = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);


%% Plot Speeds for comparison

speed_interval = 200;

figure(1);

subplot(611);
plot(data_46Hz.Speed(1.5e5:1.5e5+speed_interval));
ylabel('46 Hz')

subplot(612);
plot(data_60.Speed(2e5:2e5+speed_interval));
ylabel('60 Hz')

subplot(613);
plot(data_29.Speed(2e4:2e4+speed_interval));
ylabel('29 Hz')

subplot(614);
plot(data_660.Speed(2e4:2e4+speed_interval));
ylabel('10 Hz')

subplot(615);
plot(data_open.Speed(8e4:8e4+speed_interval));
ylabel('Open')

subplot(616)
plot(data_2_5.Speed(8e4:8e4+speed_interval));
ylabel('2.5Hz')
