%  This script loads the machine parameters for a 3hp, 180 V (L-N, peak) AC
%  induction machine into the global environment.
%  Run this code to load sample machine parameters before simulating with
%  ind.
%
%  This software is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY.  It is for educational use only.  Please do not
%  distribute or sell this software, or remove the copyright notice.  
%
%  Copyright, 1995, 1997, 1998, 2000  Steven B. Leeb

global P we Rs Rr Lm Ls Lr J Bl vds vqs vqr vdr torque_load

%% Parameters for mill motor (mill_flag = true):

% P = 2;                              % Number of pole pairs
% we = 377;                           % Base electrical frequency, rad/s (60 Hz)
% 
% Rs = 2.3931/3;                      % Stator resistance
% Rr = 1.5474/3;                      % Rotor resistance
% Lm = 0.1613;                        % Magnetizing inductance
% Ls = 0.0079;                        % Stator leakage
% Lr = 0.0079;                        % Rotor leakage
% 
% J = 0.00170;                        % Rotor inertia
% Bl = 0.0439;                        % Load damping coefficient
% 
% vds = -80;                          % D-axis stator voltage
% vqs = 191;                          % Q-axis stator voltage
% vqr = 0.0;                          % D-axis rotor voltage
% vdr = 0.0;                          % Q-axis rotor voltage
% 
% torque_load = @(t) sin(377*t);      % Load torque

%% Example parameters for a 3 Hp, 180V (L-N, peak) AC induction machine (mill_flag = false):

P = 2;                              % Number of pole pairs
we = 377;                           % Base electrical frequency, rad/s (60 Hz)

Rs = 0.435;                         % Stator resistance
Rr = 0.816;                         % Rotor resistance
Lm = 26.13/we;                      % Magnetizing inductance
Ls = 0.754/we;                      % Stator leakage
Lr = 0.754/we;                      % Rotor leakage

J = 0.089;                          % Rotor inertia
Bl = 0.00;                          % Load damping coefficient (0.02 originally)

vds = 0.0;                          % D axis stator voltage
vqs = 180.0;                        % Q axis stator voltage
vqr = 0.0;                          % D axis rotor voltage
vdr = 0.0;                          % Q axis rotor voltage

%torque_load = @(t) 1000*sin(377/4*t);      % Load torque
torque_load = @(t) 12.5*sin(377/4*t);      % Load torque




