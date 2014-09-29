clear all; close all; clc
startup
sampling_rate = 8000;
mill='heavy duty';
encoder_count = 1800;
useTrueSpeed=false;
l2l_flag = false; %PM motor true Induction False
pulley_ratio=0.9602;     % heavy duty   
Pole_pairs = 2; % Heavy Duty

%% Calculate b from air cut
% Use this section to recalculate b if necessary
data_air = process_mill_nilm('no_load',1);
data_spindown = Pre_process_three_phase_motor_data('spindown1', sampling_rate, ...
    l2l_flag, Pole_pairs, encoder_count, mill);
close all
%%
close all

speed_air = data_air.Speed_verification;
speed_air = [speed_air',speed_air']';
time_air = [0:1/8000:(size(speed_air,1)-1)/8000]';
f1 = fit(time_air,speed_air,'poly9')
speed_smooth = feval(f1,time_air');

torq_air = data_air.Torque_mag;

speed_diff = [diff(speed_smooth)]./(1/8000);
%{
matches =[];
for i=1:size(speed_diff,1)
    if (abs(speed_diff(i)))<0.1
        matches=[matches,i];
        speed_smooth(i)
    end
end
matches
%}
speed_smooth =speed_smooth(278:781);
speed_diff = speed_diff(277:780);
torq_air = [torq_air',torq_air'];
torq_air = torq_air(278:781)';
x = [0:1:size(speed_diff,1)-1];
plot(x,speed_smooth,x,speed_diff)

size(speed_smooth)
size(speed_diff)
size(torq_air)
%%

J =[0.1600:0.01:0.3];
J_calc = zeros(size(J));
b = zeros(size(speed_smooth,1),size(J,2));
stds = zeros(size(J,2),1);
%J= 0.0012;
for i=1:size(J,2);
    J_temp = J(i);
    
    b_no_j = zeros(size(speed_smooth,1),1);
    for k =1:size(speed_smooth,1);
        b(k,i) = (torq_air(k) - J_temp*speed_diff(k)) / speed_smooth(k);
        b_no_j(k) = torq_air(k)/speed_smooth(k);
    end
    time_air = [0:1/8000:(size(b,1)-1)/8000];

    spindown_index = [2.231e4 5.841e4];

    speed_spindown = data_spindown.Speed(spindown_index(1):spindown_index(2));
    time_spindown = [0:1/8000:(size(speed_spindown,1)-1)/8000]';

    % First order
    f1 = fit(time_spindown,speed_spindown,'poly1');
    p1 = f1.p1/f1.p2;
    Tc1 = -1/p1;

    J_calc(i) = Tc1*mean(b(i));
    
    %plot(time_spindown,speed_spindown)
    %text(0,max(speed_spindown),'J: '+str(J)+'\nb: '+str(b))
    stds(i,1) = std(b(81:401,i));
end
[~,I] = min(stds);
x = [0:1:size(J,2)-1];
figure(1)
subplot(2,1,1)
plot(time_air,b,'r',time_air,b_no_j,'g')
subplot(2,1,2)
plot(x,J,'r',x,J_calc,'g')

matches = [];
for i=1:size(J,2)
    if round(J(i)*1000) ==round(J_calc(i)*1000);
        matches =[matches J(i)];
    end
end
matches




%% Spindown2
%{
data = Pre_process_three_phase_motor_data('spindown2', sampling_rate, ...
    l2l_flag, Pole_pairs, encoder_count, mill);
%%
plot(data.Speed)

b = 0.0439
spindown_index = [2.32e4 5.8e4];

speed_spindown = data.Speed(spindown_index(1):spindown_index(2));
time_spindown = [0:1/8000:(size(speed1,1)-1)/8000]';
plot(time,speed1)

% First order
f1 = fit(time,speed1,'poly9')
p1 = f1.p1/f1.p2;
Tc1 = -1/p1

J = Tc1*b


%}

%% Start up
%{
data = Pre_process_three_phase_motor_data('startup4', sampling_rate, ...
    l2l_flag, Pole_pairs, encoder_count, mill);
%%
x = [0:1/8000:(size(data.Speed,1)-1)/8000]';

plotyy(x,data.Speed,x,data.Power)
torque = data.Power./data.Speed;

%% Startup2

air_index = [4.5e4 6.5e4];
start_index = [2.503e4 2.737e4];
b = 0.0080;
J = 9.444e-4;

%% Startup3
%{
air_index = [4e4 6e4];
start_index = [2.62e4 2.82e4];
b = 0.0079;
J = 9.5003e-4;
%}
%% Startup4
%{
air_index = [4e4 6e4];
start_index = [2.48e4 2.64e4];
b = 0.0078;
J=9.6702e-4;
%}
%% Startup J estimation

speed_air = data.Speed(air_index(1,1):air_index(1,2));
power_air = data.Power(air_index(1,1):air_index(1,2));
torque_air = power_air./speed_air;
b = mean(torque_air./speed_air)

speed1 = data.Speed(start_index(1):start_index(2));
power1 = data.Power(start_index(1):start_index(2));
torque1 = torque(start_index(1):start_index(2));
time = [0:1/8000:(size(torque1,1)-1)/8000]';
exp_eqn= zeros(1,size(time,1));
for i=1:size(time,1)
    exp_eqn(i) = torque1(i)/b - speed1(i);
end
exp_eqn = exp_eqn';

figure(1)
subplot(2,1,1)
plotyy(time,speed1,time,power1)
subplot(2,1,2)
plot(time,exp_eqn)

x= [1:1:size(torque1,1)]';
f1 = fit(time,exp_eqn,'exp1');
f1
%}