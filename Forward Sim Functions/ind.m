function [slopes] = ind(t, statev)
%  This script computes the state variable derivatives for a fifth
%  order model of a balanced, three phase AC induction machine.
%  The state variables are the D and Q stator and rotor fluxes, and
%  the rotor speed (wr).
%  Use this script with ODE45 to simulate the performance of the 
%  induction machine.
%  
%  Run indparam.m to load sample machine parameters before simulating.
%
%  This software is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY.  It is for educational use only.  Please do not
%  distribute or sell this software, or remove the copyright notice.  
%
%  Copyright, 1995, 1998, 2002  Steven B. Leeb
%  The variable w determines the reference frame in which the simulation
%  will be conducted.  With w = 377, the simulation variables will be in
%  a synchronously rotating reference frame.

global P we Rs Rr Lm Ls Lr J Bl vds vqs vqr vdr torque_load slip

w = 377;

Las = Ls + Lm;
Lar = Lr + Lm;
D = Lm*Lm - Las*Lar;

% Turn motor off at 1 sec 
% if t > 1.0
%    vqs = 0.0;
% end

lamqs = statev(1);
lamds = statev(2);
lamqr = statev(3);
lamdr = statev(4);
wr    = statev(5);
th    = statev(6);

% torque_load_total = torque_load(t) + Bl*wr/P;
torque_load = 10*sin(slip*th/2) + Bl*wr/P;

iqs = (Lm*lamqr - Lar*lamqs)/D;
ids = (Lm*lamdr - Lar*lamds)/D;
iqr = (Lm*lamqs - Las*lamqr)/D;
idr = (Lm*lamds - Las*lamdr)/D;

torque_mag = (3/2)*P*(lamqr*idr - lamdr*iqr);

s1 = (vqs - w*lamds - Rs*iqs);
s2 = (vds + w*lamqs - Rs*ids);
s3 = (vqr - (w - wr)*lamdr - Rr*iqr);
s4 = (vdr + (w - wr)*lamqr - Rr*idr);

% s5 = P*(torque_mag - torque_load_total)/J;
s5 = P*(torque_mag - torque_load)/J;
s6 = 377;

slopes = [s1 s2 s3 s4 s5 s6]';



