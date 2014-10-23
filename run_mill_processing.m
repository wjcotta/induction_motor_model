clear all; close all; clc;
startup;
%{
filename1 = '22augcut1';
filename2 = '22augcut2';
filename3 = '22augcut3';
filename4 = '22augcut5';
%}
%filename = '22augcut2';
filename = 'oct22_29Hz';
pulley_ratio = 0.54154-0.001 ;    % adjusted for number of flutes
%load_freq =60*2*pi;
%wr = 186.765;
%pulley_ratio = wr/load_freq;
sim_flag = false;
[data] = process_data(filename, pulley_ratio, sim_flag);

%slip = 
%{
[data1] = process_data(filename1, pulley_ratio, sim_flag);
[data2] = process_data(filename2, pulley_ratio, sim_flag);
[data3] = process_data(filename3, pulley_ratio, sim_flag);
[data4] = process_data(filename4, pulley_ratio, sim_flag);
close all
%}
%{
[data1] = process_data('no_load',0.5);
[data2] = process_data('cut2Jun23_damaged_4flute',0.5);
[data3] = process_data('cut3Jun23_intact_4flute',0.5);
close all


num_cols = 4;

%% Make torques the same length
%{
close all
[~,dex] = min(data4.torque_load);
rows = size(data4.torque_load,1);

%data4.torque_load = circshift(data4.torque_load,rows/2 - ind)
intact = fit(data1.Time_verification',data1.torque_load,'poly9','Normalize','on')
damaged = fit(data4.Time_verification',data4.torque_load,'poly9','Normalize','on')
intact

figure()
plot(intact,'b');
figure()
plot(damaged,'r');%data1.Time_verification',data1.torque_load,'g')
%plot(damaged,data4.Time_verification',data4.torque_load,'r')
%}
%%

if num_cols ==4;
    num_rows = min([size(data1.torque_load,1), size(data2.torque_load,1),...
        size(data3.torque_load,1),size(data4.torque_load,1)]);
    torque_matrix=zeros(num_rows,8);
    torque_matrix(:,1) = data1.torque_load(1:num_rows,1);
    torque_matrix(:,2) = data2.torque_load(1:num_rows,1);
    torque_matrix(:,3) = data3.torque_load(1:num_rows,1);
    torque_matrix(:,4) = data4.torque_load(1:num_rows,1);

elseif num_cols ==3;
    num_rows = min([size(data1.torque_load,1), size(data2.torque_load,1),size(data3.torque_load,1)]);
    torque_matrix=zeros(num_rows,6);
    torque_matrix(:,1) = data1.torque_load(1:num_rows,1);
    torque_matrix(:,2) = data2.torque_load(1:num_rows,1);
    torque_matrix(:,3) = data3.torque_load(1:num_rows,1);
end

if num_cols ==4;
    torque_matrix(:,5) = data1.torque_mag(1:num_rows,1);
    torque_matrix(:,6) = data2.torque_mag(1:num_rows,1);
    torque_matrix(:,7) = data3.torque_mag(1:num_rows,1);
    torque_matrix(:,8) = data4.torque_mag(1:num_rows,1);

elseif num_cols ==3;
    torque_matrix(:,4) = data1.torque_mag(1:num_rows,1);
    torque_matrix(:,5) = data2.torque_mag(1:num_rows,1);
    torque_matrix(:,6) = data3.torque_mag(1:num_rows,1);
end

%% Plot Torques
x = [1:1:num_rows]';
for k=[2:num_cols];
    torque1 = torque_matrix(:,k);
    rmse_min = 1000;
    for i=1:size(torque1,1);
        shifted = circshift(torque1,i);
        rmse = rms(shifted - torque_matrix(:,1));
    
        if rmse < rmse_min;
            rmse_min = rmse;
            shift_index = i;
            shift_index;
        end
    end
    torque_matrix(:,k) = circshift(torque1,shift_index);
end
%% 3 cols
if num_cols ==3;
    time = [0:1/8000:(num_rows-1)/8000]';
    intact = fit(time,torque_matrix(:,3),'poly9','Normalize','on')
    damaged = fit(time,torque_matrix(:,2),'poly9','Normalize','on')
    no_load = fit(time,torque_matrix(:,1),'poly9','Normalize','on')

    close all
    figure()
    hold on
    plot(intact,'g')
    plot(damaged,'r')
    plot(no_load,'k')
    axis([0.01 0.12 0.5 1.8])
    legend('intact', 'damaged', 'no load');
    ylabel('Load Torque')
    xlabel('Time (1/8000 s)')
    title('4 Flute')
    hold off
end
%% 4 cols
if num_cols ==4;
    time = [0:1/8000:(num_rows-1)/8000]';
    intact005 = fit(time,torque_matrix(:,1),'poly9','Normalize','on')
    intact010 = fit(time,torque_matrix(:,2),'poly9','Normalize','on')
    intact010v2 = fit(time,torque_matrix(:,3),'poly9','Normalize','on')
    damaged0075 = fit(time,torque_matrix(:,4),'poly9','Normalize','on')

    close all
    figure()
    hold on
    plot(intact005,'g')
    plot(intact010,'b')
    plot(intact010v2,'k')
    plot(damaged0075,'r');
    axis([0 0.0615 0.5 1.8])
    legend('intact 0.05', 'intact 0.10','intact 0.10','damaged 0.075');
    ylabel('Load Torque')
    xlabel('Time (1/8000 s)')
    title('2 Flute')
    hold off
    figure()
    plot(damaged0075,'r');
end
%%

%%
for k=[num_cols+2:2*num_cols];
    torque1 = torque_matrix(:,k);
    rmse_min = 1000;
    for i=1:size(torque1,1);
        shifted = circshift(torque1,i);
        rmse = rms(shifted - torque_matrix(:,num_cols+1));
    
        if rmse < rmse_min;
            rmse_min = rmse;
            shift_index = i;
            shift_index;
        end
    end
    torque_matrix(:,k) = circshift(torque1,shift_index);
end


figure()
if num_cols ==3;
    subplot(2,1,1)
    plot(x,torque_matrix(:,1),'k',x,torque_matrix(:,2),'r',...
        x,torque_matrix(:,3),'g')
    legend('no load','damaged','intact')
    ylabel('Load Torque')
    xlabel('Time (1/8000 s)')
    title('Changing Parameters 4 Flute')
    subplot(2,1,2)
    plot(x,torque_matrix(:,4),'k',x,torque_matrix(:,5),'r',...
    x,torque_matrix(:,6),'g')
    legend('no load','damaged','intact')
    ylabel('Magentic Torque')
    xlabel('Time (1/8000 s)')
    title('Changing Parameters 4 Flute')
    
elseif num_cols ==4;
    subplot(2,1,1)
    plot(x,torque_matrix(:,1),'r',x,torque_matrix(:,2),'k',...
        x,torque_matrix(:,3),'b',x,torque_matrix(:,4),'g')
    legend('1: 0.05','2: 0.10','3: 0.10','4: 0.05')
    ylabel('Load Torque')
    xlabel('Time (1/8000 s)')
    title('Changing Parameters 2 Flute')
    
    subplot(2,1,2)
    plot(x,torque_matrix(:,5),'r',x,torque_matrix(:,6),'k',...
        x,torque_matrix(:,7),'b',x,torque_matrix(:,8),'g')
    legend('1: 0.05','2: 0.10','3: 0.05','4: 0.075 damaged')
    ylabel('Magentic Torque')
    xlabel('Time (1/8000 s)')
    title('Changing Parameters 2 Flute')
end
%%
%}