function [] = slip_check(data_speed,elec_freq,Pole_pairs,varagin)
% Calculates the slip of motor based on speed and throws error if it
% exceeds 10%. Function can recieve either a data object from the motor
% model inversion code. Or the user can input speed, electrical r 
if nargin ==1;
    mean_speed = data_speed.wr;
    speed = data_speed.Speed;
    elec_freq = data_speed.we_mean;
    Pole_pairs = data_speed.P;
end

if nargin >1;
    mean_speed = 0;
    for i = 1:1000:length(speed)-1000;
        if mean(speed(i:i+1000)) > mean_speed ;
            mean_speed = mean(speed(i:i+1000))
        end
    end
end
slip = 1 -  mean_speed/(elec_freq/Pole_pairs)
if abs(slip) > 0.10;
    error(strcat('SlipError: Motor slip exceeds 10%, something is probably',...
          ' wrong, as motors are generally unstable above 7% slip'))
      
end
end