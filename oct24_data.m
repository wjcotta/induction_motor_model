% Script for testing oct24_data
clear all; close all; clc

filename = 'oct24_46Hz';
pulley_ratio =.54154;
sim_flag= false;
scan_flag=false;
process_data(filename, pulley_ratio, sim_flag,scan_flag)