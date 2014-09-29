clear all;close all; clc;
encoder_count = 1800;
mill= 'heavy_duty'
P=2;
l2l_flag = 'False';
fs = 8000;

data8 = preprocess_motor_data('startup8', fs, l2l_flag, P, encoder_count, mill)
data1 = preprocess_motor_data('startup1', fs, l2l_flag, P, encoder_count, mill,1)

plot(data1.Speed)