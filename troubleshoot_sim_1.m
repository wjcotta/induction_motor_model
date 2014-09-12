clear all; close all; clc
startup

% Turn on mill motor parameters in 'inparam'.
generate_sim_data
sim_data = load('sim_data');

pulley_ratio = 1;
sim_flag = false;
pulley_ratio = 0.54154*pulley_ratio;    % adjusted for number of flutes
fs = 8000;                              % sampling rate [Hz]
l2l_flag = false;                       % PM motor = true, induction = false
P = 2;                                  % number of pole pairs, 2 for heavy duty
encoder_count = 1800;                   % encoder count
mill = 'heavy duty';                    % 'heavy duty' or 'lighty duty'
avg_time_window = 20;                   % data time window
tune = 1;                               % 0.1 - 100
useTrueSpeed = false;                   % PM motor = true, induction = false

real_data = preprocess_motor_data('22augcut1', fs, l2l_flag, P, encoder_count, mill);

x=[1:1:1001];

% Check sim voltages (using mill motor params) vs preprocess voltages.
figure(1)
plot(x,sim_data.data(9000:10000,1),'r',x,real_data.Vwye(9000:10000,1),'g')

% Check sim currents (using mill motor params) vs preprocess currents.
figure(2)
plot(x,sim_data.data(9000:10000,4),'r',x,real_data.Iabc(9000:10000,1),'g')

