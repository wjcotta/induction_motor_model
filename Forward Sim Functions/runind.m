%  This script runs the AC induction motor simulation during free acceleration
%
%  This software is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY.  It is for educational use only.  Please do not
%  distribute or sell this software, or remove the copyright notice.  
%
%  Copyright, 1995, 1998, 2000  Steven B. Leeb

clear all; clc; close all;
indparam

t_start = 0;
t_final = 45;
fs = 8000;                          
T = 1/fs;                          
tspan = t_start:T:t_final;

state0 = [0 0 0 0 0 0]';

[t, statev] = ode45('ind', tspan, state0);

figure(1)
plot(t, statev(:,5));
xlabel('Time (s)');
ylabel('Electrical Rotor Speed (rad/s)');
axis([0 45 0 250])
title('Motor Spin-up');

figure(2)
plot(t, statev(:,6));
xlabel('Time (s)');
ylabel('Mechanical Shaft Angle (rad)');
axis([0 45 0 6000])
title('Motor Spin-up');

ind1 = find(t,1.000);
ind2 = find(round(statev(:,5)), 150, 'last');
f1 = fit(t(ind1:ind2,:), statev(ind1:ind2,5), 'poly1');
Tc = -1/(f1.p1/f1.p2);
