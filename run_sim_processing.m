clear all; close all; clc;

generate_sim_data
clear all;

filename = 'sim_data';
pulley_ratio = 1;
      
sim_flag = true;

[data] = process_data(filename, pulley_ratio, sim_flag);

% Simulation parameters defined in 'indparam'
% Rs = 0.435, Rr = 0.816, Lm = 26.13/377 = 0.069, Lal = 0.754/377 = 0.002

calc_params = data.Parameters
slip = data.wr/(data.we/data.P)