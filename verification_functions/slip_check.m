function [] = slip_check(iabc,speed,elec_freq,Pole_pairs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h = 100;
phase_rms=zeros(length(speed)-h,3);
%% Find Transients
for i=1:length(speed)-h;
    phase_rms(i,:) = rms(iabc(i:i+h,:));
end

size(phase_rms)

figure(1)
subplot(2,1,1)
plot(iabc)
subplot(2,1,2)
plot(phase_rms)
end


