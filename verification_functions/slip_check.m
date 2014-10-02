function [] = slip_check(speed,elec_freq,Pole_pairs,varagin)
%UNTITLED Summary of this function goes here
%   Calculates the slip of motor based on speed 
mean_speed = 0;
for i = 1:1000:length(speed)-1000;
    if mean(speed(i:i+1000)) > mean_speed ;
        mean_speed = mean(speed(i:i+1000));
    end
end
slip = 1 -  mean_speed/(elec_freq/Pole_pairs)
if slip > 0.10;
    error(strcat('SlipError: Motor slip exceeds 10%, something is probably',...
          ' wrong, as motors are generally unstable above 7% slip'))