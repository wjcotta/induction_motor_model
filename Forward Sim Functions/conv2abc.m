function [idq, vabc, iabc, torque_mag, power] = conv2abc(t, statev)
% This script transforms the DQ stator and rotor fluxes computed using
% ind.m and ode45 back into laboratory frame stator currents and voltages,
% e.g., ias and vas for phase a. 
%
% The variable w determines the reference frame in which the simulation
% will be conducted.  With w = 377, the simulation variables will be in
% a synchronously rotating reference frame.
%
% On return, the output matrices m and m2 contain:
%	m = [ids iqs idr iqr];
%	m2 = [T ias ibs ics vas vbs vcs];
% This script also plots the simulated rotor torque versus speed on return.
%
% This software is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY.  It is for educational use only.  Please do not
% distribute or sell this software, or remove the copyright notice.  
%
% Copyright, 1995, 1998, 2002  Steven B. Leeb

global P we Rs Rr Lm Ls Lr J Bl vds vqs vqr vdr torque_load

% Select reference frame (see ind.m).
w = 377; 

Las = Ls + Lm;
Lar = Lr + Lm;
D = Lm*Lm - Las*Lar;

lamqs = statev(:,1);
lamds = statev(:,2);
lamqr = statev(:,3);
lamdr = statev(:,4);
wr    = statev(:,5); 
th    = statev(:,6);

ids   = (Lm*lamdr - Lar*lamds)/D;
iqs   = (Lm*lamqr - Lar*lamqs)/D;
idr   = (Lm*lamds - Las*lamdr)/D;
iqr   = (Lm*lamqs - Las*lamqr)/D;

% Calculate power in dq0 reference frame. 
Pdq = vds*ids + vqs*iqs;

torque_mag = (3/2)*P*(lamqr.*idr - lamdr.*iqr);

vas = cos(th).*vds - sin(th).*vqs;
vbs = cos(th - 2*pi/3).*vds - sin(th - 2*pi/3).*vqs;
vcs = cos(th + 2*pi/3).*vds - sin(th + 2*pi/3).*vqs;
ias = cos(th).*ids - sin(th).*iqs;
ibs = cos(th - 2*pi/3).*ids - sin(th - 2*pi/3).*iqs;
ics = cos(th + 2*pi/3).*ids - sin(th + 2*pi/3).*iqs;

% Multiply by leading factor of sqrt(2/3) for power-invariant dq0 transform. 
f = sqrt(2/3);
vas = f*vas;
vbs = f*vbs;
vcs = f*vcs;
ias = f*ias;
ibs = f*ibs;
ics = f*ics;

% Calculate power in abc reference frame as verification.
Pabc = vas.*ias + vbs.*ibs + vcs.*ics;

idq = [ids iqs idr iqr];
vabc = [vas vbs vcs];
iabc = [ias ibs ics];
power = [Pdq Pabc];
plot(3600*wr/(P*377), torque_mag);
plot(t,power(:,1),t,power(:,2));


