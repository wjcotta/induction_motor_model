clear all; close all; clc;

load 22augcut1.mat

Ia = data(:,1);
Ib = data(:,2);
Ic = data(:,3);

Ic_ver = -(Ia + Ib);

plot(Ic, 'b')
hold on 
plot(Ic_ver, 'r')