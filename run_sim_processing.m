clear all; close all; clc;

generate_sim_data
clear all;

filename = 'sim_data';
pulley_ratio = 1;
pulley_ratio = 0.5*pulley_ratio;        
sim_flag = true;

[data] = process_data(filename, pulley_ratio, sim_flag);
