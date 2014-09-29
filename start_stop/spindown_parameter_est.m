clear all;close all;clc;

filename = 'spindown2';
% spindown2
bounds = [3.25e4,9.424e4];
stop_8 = [2.299e4, 4.348e4];


encoder_count = 1800;
fs = 8000;
t_shift = @(a,t0,fs)ifft(fft(a).*exp(-1i*2*pi*t0*(ifftshift((0:length(a)-1)' -ceil((length(a)-1)/2))/length(a)*fs)));

load(filename);
temp = data;
clear data



data.raw = temp;
data.time = [bounds(1,1)/fs : 1/fs : bounds(1,2)/fs];
data.current = data.raw(:,1);
data.voltage = data.raw(:,2);
data.speed8 = data.raw(bounds(1,1):bounds(1,2),3);
data.speed1 = data.raw(bounds(1,1):bounds(1,2),4);





divisor = 8;
data.speed8 = speed_clean(data.speed8,divisor,encoder_count);
ms = mean(data.speed8/2/pi);
dtS = (-1/(ms*encoder_count/divisor));
data.speed8 = t_shift(data.speed8,dtS,fs);

divisor = 1;
data.speed1 = speed_clean(data.speed1,divisor,encoder_count);

num_del =0;
temp = data.speed1;
for i=1:length(temp)
    [C,I] = max(temp);
    if C == inf
        max_inf = I;
        temp = [temp(1:I-1,1);temp(I+1:end,1)];
        num_del = num_del + 1;
    end
end
max(temp)
ms = mean(temp)
dtS = (-1/(ms*encoder_count/divisor))
data.speed1 = t_shift(temp,dtS,fs);

%%
% spindown2
stop_1 = [4.348e4+num_del,length(data.speed1)];

speed_stop = [data.speed8(stop_8(1,1):stop_8(1,2)) ; ...
              data.speed1(stop_1(1,1):stop_1(1,2)) ];

          

close all
fc = 1/1800;
Wn = (2/fs)*fc;
b = fir1(20,Wn,'low',kaiser(21,3));
speed_stop = filter(b,1,speed_stop);
speed_stop = speed_stop(50:end);
x = [0:1/fs:(length(speed_stop)-1)/fs]';
plot(x,speed_stop)
%%
f1 = fit(x,speed_stop,'poly1')

p1 = f1.p1/f1.p2;
Tc = -1/p1

h = 1/fs;
plot(diff(speed_stop)/h)


    
          