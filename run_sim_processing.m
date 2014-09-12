clear all; close all; clc;

generate_sim_data
clear all;

filename = 'sim_data';
pulley_ratio = 1;
sim_flag = true;

[data] = process_data(filename, pulley_ratio, sim_flag);
