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

global P we Rs Rr Lm Ls Lr J Bl vds vqs vqr vdr torque_load slip

%% Parameters for mill motor (mill_flag = true):

P = 2;                              % Number of pole pairs
we = 377;                           % Base electrical frequency, rad/s (60 Hz)
 
Rs = 2.4480;                      % Stator resistance
Rr = 1.5461;                      % Rotor resistance
Lm = 0.1619;                        % Magnetizing inductance
Ls = 0.0075;                        % Stator leakage
Lr = 0.0075;                        % Rotor leakage

%J = 0.220;
J = 0.2350;                        % Rotor inertia
Bl = 0.0441;
%Bl = 0.04;
%Bl = 0.0439;                        % Load damping coefficient

vds = 0;                          % D-axis stator voltage
vqs = 180;                          % Q-axis stator voltage
vqr = 0.0;                          % D-axis rotor voltage
vdr = 0.0;                          % Q-axis rotor voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIX THE SLIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slip =.9354;                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIX THE SLIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% torque_load = @(t) sin(377*(1/pulley_ratio)*t);      % Load torque

%% Example parameters for a 3 Hp, 180V (L-N, peak) AC induction machine (mill_flag = false):
% 
% P = 2;                              % Number of pole pairs
% we = 377;                           % Base electrical frequency, rad/s (60 Hz)
% 
% Rs = 0.435;                         % Stator resistance
% Rr = 0.816;                         % Rotor resistance
% Lm = 26.13/we;                      % Magnetizing inductance
% Ls = 0.754/we;                      % Stator leakage
% Lr = 0.754/we;                      % Rotor leakage
% 
% J = 0.089;                          % Rotor inertia
% Bl = 0.004;                         % Load damping coefficient 
% 
% vds = 0.0;                          % D axis stator voltage
% vqs = 180.0;                        % Q axis stator voltage
% vqr = 0.0;                          % D axis rotor voltage
% vdr = 0.0;                          % Q axis rotor voltage
% 
% slip = 0.9975;
% % torque_load = @(t) 10*sin(slip*377/2*t);      % Load torque
% 




