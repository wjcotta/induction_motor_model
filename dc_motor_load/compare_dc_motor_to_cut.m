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
filename = '22augcut2';
data_cut = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

filename = 'oct22_60Hz';
data_60 = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

filename = 'oct22_29Hz';
data_29 = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

filename = 'oct22_10Hz';
data_10 = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill);

%% Plot Speeds for comparison

speed_interval = 1000;

figure(1);

subplot(411);
plot(data_cut.Speed(3e5:3e5+speed_interval));

subplot(412);
plot(data_60.Speed(2e5:2e5+speed_interval));

subplot(413);
plot(data_29.Speed(2e5:2e5+speed_interval));

subplot(414);
plot(data_10.Speed(2e5:2e5+speed_interval));