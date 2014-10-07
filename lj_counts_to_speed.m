function [ speed ] = lj_counts_to_speed( counts,divisor,encoder_counts );
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fs=8000;
clock = 48e6; % Labjack system clock
bins = 2^16;
x = [0:1:size(counts,1)-1];
A = counts;
B = medfilt1(A,5,1e4);
for i=1:size(B,1)-1;
    if abs(B(i,1)-B(i+1,1)) > bins/2;
        B(i+1:end,1) = B(i+1:end,1)+bins;
    end
end

speed=0;
%
%plots = 3;
%figure(1)
%subplot(plots,1,1)
%scatter(x,A);
%subplot(plots,1,2)
%scatter(x,B);
%subplot(plots,1,3)

lj_divisor = 732.42;
speed= (clock*divisor)./(encoder_counts*B); % Speed in frequency
speed_w= 2*pi*speed;
speed_filtered = medfilt1(speed_w,5,1e4);
x = [0:1:size(counts,1)-1]';
%figure(1)
%plot(speed_w)
%pause
time = x(9e4:end,1)/fs + 9e4/fs

speed_fit = speed_w(9e4:end,1)
fit_obj = fit(time,speed_fit,'exp1')

close all

plot(fit_obj);
hold on
plot(time,speed_fit)
hold off
%plot(time,speed_fit);
%hold off
%}
%{
plots = 4;
figure(1)

subplot(plots,1,1);
scatter(x,counts);
ylabel('Raw from labjack')

subplot(plots,1,2);
scatter(x,speed);
ylim([0 60]);
ylabel('Converted to Speed (RPS)')

subplot(plots,1,3);
scatter(x,speed_w);
ylim([0 210]);
ylabel('Converted to Speed (Rad/S)')

subplot(plots,1,4);
scatter(x,speed_filtered);
ylim([0 210]);
ylabel('Median Filtered Speed (Rad/S)')

figure(2)
scatter(x,speed_filtered);
ylim([0 210]);
ylabel('Median Filtered Speed (Rad/S)')
end
%}




end
