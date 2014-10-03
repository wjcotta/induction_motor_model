close all;clear all;clc
load('dc_motor_data')
fs= 8000;
encoder_count = 1800;
we = 377;
P=2;
slip = 1;
%% Import data, and build data object
data.dc1 = dcmotor1; % corrupted shaft tachometer data
data.dc2 = dcmotor2; 
data.dc3 = dcmotor3; 
data.dc4 = dcmotor4; 
data_fields = fields(data);
data_size = length(data_fields);

% Create a time vector for each set of data
time=0;
for i =1:data_size;
    name = char(data_fields(i));
    temp = getfield(data,name);
    temp = [0:1/fs:(size(temp,1)-1)/fs]';
    time = setfield(time,name,temp);
end

%% Bounds of Steady State Operation
% col1 = begin ss
% col2 = trans off
% col3 = below 10 rad/s
% col4 = end of usable data
bounds.dc1 = [3.045e4 , 5.272e4, 0,0];
bounds.dc2 = [3e4     , 7.095e4, 0,0];
bounds.dc3 = [3.1e4   , 6.537e4, 0,0];
bounds.dc4 = [2.6e4   , 7.836e4, 0,0];

% Make all speeds positive and determine speed to voltage ratio >bounds(3);
for i = 1:data_size;
    name = char(data_fields(i));
    temp_data = getfield(data,name);
    temp_bounds = getfield(bounds,name);
    
    start_mean = mean(temp_data(1:100,1));
    
    % If the motor is rotating backwards, flip the data right side up.
    if mean(temp_data(:,1) - start_mean) < 0;
        temp_data(:,1) = temp_data(:,1) + 2*abs(temp_data(:,1)-start_mean);
    end
    
    ss_mean = mean(temp_data(temp_bounds(1):temp_bounds(2),1));
    
    % Compute coefficient to turn counts to speed
    slope =  (ss_mean-start_mean) / (we*slip/(P));
    temp_data(:,1) = (temp_data(:,1)-start_mean)/slope;
    
    % Find where speed drops below 10 rad/s
    temp_bounds(3) = temp_bounds(2)+find(temp_data(temp_bounds(2):end,1)<10,1);
    
    bounds = setfield(bounds,name,temp_bounds);
    data = setfield(data,name,temp_data);
    
end
%% Plot Verification Section
% Plot data to verify it is formatted correctly
% for i = 1:data_size;
%     name = char(data_fields(i));
%     temp = getfield(data,name);
%     subplot(data_size,1,i);
%     plot(temp(:,1))
% end

% Plot speed and currents
% for i = 1:data_size;
%     name = char(data_fields(i));
%     temp = getfield(data,name);
%     temp_time = getfield(time,name);
%     turn_off = getfield(bounds,name);
%     subplot(data_size,1,i);
%     plot(temp_time(turn_off(2):turn_off(3)), temp(turn_off(2):turn_off(3),1))
% end
%%
fit_linear = 0;
fit_exp = 0;

for i =1:data_size;
    
    name = char(data_fields(i));
    temp_data = getfield(data,name);
    temp_time = getfield(time,name);
    temp_bounds = getfield(bounds,name);
    delete_ind = [size(temp_data,1)];
    
    for k =1:size(temp_data,1);
        if isnan(temp_data(k,1))==1
            delete_ind = [delete_ind,k];
            
        end
    end
    
    temp_bounds(4) = min(delete_ind)-1;
    temp_data = temp_data(temp_bounds(2):temp_bounds(3),1);
    temp_time = temp_time(temp_bounds(2):temp_bounds(3),1);
    fit1 = fit(temp_time,temp_data(:,1),'poly1');
    fit2 = fit(temp_time,temp_data(:,1),'exp1');
    
    fit_linear = setfield(fit_linear,name,fit1);
    fit_exp = setfield(fit_exp,name,fit2);
    bounds = setfield(bounds,name,temp_bounds);
end
close all
%% Plot DC Motor
name = 'dc3'
temp_data = getfield(data,name);
temp_time = getfield(time,name);
temp_bounds = getfield(bounds,name);
temp_linear = getfield(fit_linear,name);
temp_exp = getfield(fit_exp,name);
figure(1);
subplot(2,1,1)
plot(temp_time,temp_data(:,1));
subplot(2,1,2)
hold on
plot(temp_time(temp_bounds(2):temp_bounds(3),:),...
     temp_data(temp_bounds(2):temp_bounds(3),1),'*g');
plot(temp_linear,'b');
plot(temp_exp,'k');
hold off

%% Use Shaft Encoder
t_shift = @(a,t0,fs) ifft(fft(a).*exp(-1i*2*pi*t0*(ifftshift((0:length(a)-1)' -ceil((length(a)-1)/2))/length(a)*fs)));

% Divisor = 8
figure(2)

divisor = 1;
Speed = speed_clean(temp_data(:,3),divisor,encoder_count);
ms = mean(Speed(temp_bounds(1):temp_bounds(2))/2/pi);
dtS = (-1/(ms*encoder_count/divisor));

subplot(2,1,1)
plot(temp_time(temp_bounds(1):temp_bounds(2)),Speed(temp_bounds(1):temp_bounds(2)));
ylim([-10 400]);

subplot(2,1,2)
Speed = t_shift(Speed,dtS,fs)
plot(temp_time,Speed)
ylim([-10 400]);

% Divisor = 1 
figure(3)

divisor = 8;
Speed = speed_clean(temp_data(:,3),divisor,encoder_count);
ms = mean(Speed(temp_bounds(1):temp_bounds(2))/2/pi);
dtS = (-1/(ms*encoder_count/divisor));

subplot(2,1,1)
plot(temp_time(temp_bounds(1):temp_bounds(2)),Speed(temp_bounds(1):temp_bounds(2)));
ylim([-10 400]);

subplot(2,1,2)
Speed = t_shift(Speed(temp_bounds(1):temp_bounds(2)),dtS,fs)
plot(temp_time(temp_bounds(1):temp_bounds(2)),real(Speed))
ylim([-10 400]);