clear all; close all; clc;

generate_sim_data
clear all;

filename = 'sim_data';
pulley_ratio = 1;
sim_flag = true;

[data] = process_data(filename, pulley_ratio, sim_flag);

idq = load('sim_idq');
idq = idq.idq;

V_alpha = real(data.V_clark);
V_beta = imag(data.V_clark);

% idq has the form [ids iqs idr iqr].

%% Check stator current transformation
%{
% Check ids (sim) vs Ids (preprocess)
figure(1)
plot(idq(:,1))
hold on
plot(real(data.Idq),'r')

% Check iqs (sim) vs Iqs (preprocess)
figure(2)
plot(idq(:,2))
hold on
plot(imag(data.Idq),'r')

% Plot Ids DeMod
figure(3)
plot(real(data.IdqDeMod))

% Plot Iqs DemMod
figure(4)
plot(imag(data.IdqDeMod))
%}

%% Check rotor current transformation
%{
% Check idr (sim) rs Idr_check (preprocess)
figure(1)
plot(idq(:,3))
hold on
plot(real(data.Idq_rotor_check),'r')

% Check iqr (sim) vs Iqr_check (preprocess)
figure(2)
plot(idq(:,4))
hold on
plot(imag(data.Idq_rotor_check),'r')

% Plot Idr (after phi_offset)
figure(3)
plot(real(data.Idq_rotor))

% Plot Iqr (after phi_offset)
figure(4)
plot(imag(data.Idq_rotor))
%}

