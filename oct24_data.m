% Script for testing oct24_data
clear all; close all; clc

filename = 'oct24_constant';
pulley_ratio =.54154;
sim_flag= false;
scan_flag=true;
data = process_data(filename, pulley_ratio, sim_flag,scan_flag);

%%
%{
close all
figure(1)
plotyy(data.Time,data.Speed,data.Time,data.dc_motor_current)
%}