function [data] = process_data(filename, pulley_ratio, sim_flag,scan_flag)
%% Input parameters and run processing scripts
if nargin<4;
    scan_flag=false;
end
fs = 8000;                              % sampling rate [Hz]
l2l_flag = false;                       % PM motor = true, induction = false
P = 2;                                  % number of pole pairs, 2 for heavy duty
encoder_count = 1800;                   % encoder count
mill = 'heavy duty';                    % 'heavy duty' or 'lighty duty'
avg_time_window = 20;                   % data time window
tune = 1;                               % 0.1 - 100
useTrueSpeed = false;                   % PM motor = true, induction = false
replicate_current=true;



if sim_flag == false
    data = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill,replicate_current);
else
    data = preprocess_sim_data(filename, fs, l2l_flag, P, encoder_count, mill);
end

data = find_stable_region(data, l2l_flag, avg_time_window);
if pulley_ratio == -1;
    load_wr = -1;
else;
    load_wr = data.wr*pulley_ratio;
end
slip_check(data);


%% Using DC Motor
%load_wr = 2*pi*10;%*pulley_ratio
%% 

%% Inversion Functions
if scan_flag == false;
    data = Decompose_Three_Phase_Motor_Data(data,load_wr,data.indl_schantz)
    
    %data = calc_params_run_sim(data, load_wr, P, tune, useTrueSpeed);
    data = speed_inversion_error_space_sim_fixed(data, load_wr, P, tune, useTrueSpeed);
    data.Time_verification = [0:1/8000:(size(data.Idq_Verification,1)-1)/8000];

elseif scan_flag == true;
    load_wr_range = [10,400];
    % Last frequency test 31.71 rad/s
    d_wr = 0.01;
    data = scan_load_wr(pulley_ratio,load_wr_range,d_wr,data,tune)
    return
end 
    
%% Estimate torque

Lm = data.Parameters(1,3);
Lal = data.Parameters(1,4);
Las = Lm + Lal;
Lar = Lm + Lal;
P = data.P;
D = Lm^2 -(Las*Lar);

lam_s = data.lambdas(:,2) + j*data.lambdas(:,1);

lam_r = data.lambdas(:,4) + j*data.lambdas(:,3);

J = [9.444e-4 9.5003e-4 9.6702e-4];
J = mean(J);
b = [0.0078 0.0080 0.0079];
b = mean(b);

speed = data.Speed_verification;
speed_diff = diff(speed)/(1/8000);
speed_diff = vertcat([0],speed_diff);
idq = data.Idq_Verification;
for i=1:size(idq,1);
    idq(i,1) = -real(idq(i,1)) + 1i*imag(idq(i,1));
end
lam_r(1,:);
for i=1:size(lam_r,1)
torque_mag(i,1) =  (Lm*P)/D * imag(conj(lam_s(i))*lam_r(i));
%torque_mag(i,1) = ((Lm*P)/D)*(imag(lam_r(i))*imag(lam_s(i,1))*(-1));
torque_load(i,1) = torque_mag(i,1) - speed(i,1)*b -J*speed_diff(i,1);
end

data.torque_mag = torque_mag;
data.torque_load = torque_load;
%}
end
