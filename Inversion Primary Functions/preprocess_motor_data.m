function out = preprocess_motor_data(filename, fs, l2l_flag, P, encoder_count, mill,divisor_flag,varagin)

% This function converts the raw output of the Labjack DAQ device into physical units suitable for calculation, does time-shifting to account for interchannel delay,
% performs Clark and Park transformations (accounting for variable electric frequency over course of experiment), and produces an estimate of instantaneous electric frequency.

if (nargin < 5)
    encoder_count = 1800;
end
if (nargin < 7)
    divisor_flag = false;
end
%% Load data
load(filename)
out.filename = filename;
% L = length(data);

%% Set calibration constants for the three-phase HP power supply and analyzer.
% Schantz box configured so all voltage channels measure the same voltage phase, likewise for the current channels.

%Slopes = [0.000935547853831 0.001010673290357 0.001015083089409 0.091020193288373 0.091089882503841 0.090792761167529];
Slopes = [0.001010673290357  0.000935547853831 0.001015083089409 0.091020193288373 0.091089882503841 0.090792761167529];


%Slopes = [0.00062334 0.00063414 0.0006251 0.062319 0.05787 0.063953];
%Slopes = [0.0016138 0.0016195 0.0016240 0.090972 0.091101 0.090825];
%Slopes = [0.001041415057298 0.001072327641887 0.001014188488812 0.091013485692896 0.091075712704749 0.090784237424945];

% Offsets = [32940.42288  33037.0232  33051.21264 32983.8873  32987.6873  32987.68735];

% {
% if ~l2l_flag
%     Slopes = [0.00101010738971818, 0.00101176539509996, 0.00100682083611568, 0.0914369606121276,  0.0913493710967854,  0.0916521868611356];
%     % These are the mill box constants.
% else
%     Slopes = [0.00101139,  0.00101397, 0.0010135, 1, 1, 1];
%     % Scaling here for AIN0-AIN2 different than above due to inline AA filter.
% end
% }

%% Filter data to remove stator bar harmonics.

Wp = 1200/(fs/2);                   % Passband edge frequency
Ws = 1600/(fs/2);                   % Stopband edge frequency
[n, Wn] = buttord(Wp, Ws, 1, 10);   % Calculates the minimum order and cutoff frequency of a Butterworth filter required to meet desired specifications.
[b, a] = butter(n, Wn);             % Returns filter numberator and denominator coefficients in vectors b and a, respectively.

% Use outputs from butter function to filter the data.
data(:,1) = filtfilt(b,a,data(:,1));
data(:,2) = filtfilt(b,a,data(:,2));
data(:,3) = filtfilt(b,a,data(:,3));
data(:,4) = filtfilt(b,a,data(:,4));
data(:,5) = filtfilt(b,a,data(:,5));
data(:,6) = filtfilt(b,a,data(:,6));

[~, sI, eI] = cycle_mean(data(:,4));
% Finds mean between first maximum and last minimum for a rough estimate of the cycle mean.
% Simply taking mean of entire signal can be misleading if signal is not composed of integer number of periods.
% This ensures we start at a rise and end at a fall
% We now have periodic signal vectors in the long data, which will cause less edge effect with fft-based processing later.

AIN0 = data(sI:eI,1);
AIN1 = data(sI:eI,2);
AIN2 = data(sI:eI,3);
AIN3 = data(sI:eI,4);   % Vab if line to line measurement
AIN4 = data(sI:eI,5);   % Vbc if line to line measurement
AIN5 = data(sI:eI,6);   % Vca if line to line measurement

% AIN3 generally always a voltage measurement.

%% Generate Speed Vector from Encoder Counts

% Speed measurememnt done by counting clock pulses at 48 mHz between 2 pairs of rising or falling edges from a single encoder line.
% Thus calculated speed is most acurate midway between the edge pairs from which it is derived.
% The labjack latches the time between edge pairs after each applicable edge, making the speed derived from the reported time
% always at least 1/5 an edge pair duration old and up to 1 1/2 edge-pair durations olds (depending on when the latch contents
% are read by the data logger), for an average of 1 edge-pair duration time delay in the derived speed.

if size(data,2) == 7
    Counts = data(sI:eI,7);
    divisor = 8;
    if divisor_flag ~= false;
        divisor = divisor_flag
    end
    Speed = speed_clean(Counts,divisor,encoder_count);
    ms = mean(Speed/2/pi);
    dtS = (-1/(ms*encoder_count/divisor));
end

if size(data,2) == 8
    coef = [-0.000152179987790691, 7.44468811677644];
    Counts = data(sI:eI,8);
    divisor = 8;
    Speed = speed_clean(Counts,divisor,encoder_count);
    ms = mean(Speed/2/pi);
    dtS = (-1/(ms*encoder_count/divisor));
    Ac = data(sI:eI,7);
    mA = medfilt1(Ac,5);
    temp = find(abs(Ac-mA)>2*std(Ac));
    Ac(temp) = mA(temp);
    out.Aux_channel = Ac*coef(1) + coef(2);
end

if size(data,2) == 11
    divisor = 32;
    Counts = data(sI:eI,10:11);
    Speed = speed_clean(Counts,divisor,encoder_count);
    ms = mean(Speed/2/pi);
    dtS = (-1/(ms*encoder_count/divisor));
end

if size(data,2) == 10

    nom_Slope = 1.5629e-4;
    nom_Offset = -5.1760;

    out.pu_volts_base = 190/sqrt(3);
    out.pu_current_base = 20.88;

    out.TI_pwmdac_chan0 = ((nom_Slope.*data(sI:eI,7)) + nom_Offset-1.65)/1.65;
    out.TI_pwmdac_chan1 = ((nom_Slope.*data(sI:eI,8)) + nom_Offset-1.65)/1.65;
    out.TI_pwmdac_chan2 = ((nom_Slope.*data(sI:eI,9)) + nom_Offset-1.65)/1.65;
    out.TI_pwmdac_chan3 = ((nom_Slope.*data(sI:eI,10)) + nom_Offset-1.65)/1.65;

end

clear data

%% Convert Labjack bits to currents

Volts_LJ = @(Bits) (Bits/65536)*(5.07+5.18)-5.18;
I_measured = @(Bits) 1000*(Volts_LJ(Bits)/155);
% The resistance is set at 155 ohms using the switches in the nilm box (all resistors turned on).
% 1000 is the conversion for the LA-55-P current transducer from measure current to output current.

currents = [AIN0, AIN1, AIN2];

for i=1:3;
    for k=1:length(currents);
        temp(k,i)=I_measured(currents(k,i));
    end
end

AIN0 = temp(:,1) - mean(temp(:,1));
AIN1 = temp(:,2) - mean(temp(:,2));
AIN2 = temp(:,3) - mean(temp(:,3));

clear temp currents

%% Use Slopes for Voltage

AIN3 = (AIN3 - mean(AIN3))*Slopes(4);
AIN4 = (AIN4 - mean(AIN4))*Slopes(5);
AIN5 = (AIN5 - mean(AIN5))*Slopes(6);

dt = 11.55/1e6;
% Best fit interchannel delay from test runs with 3 consecutive channels reading the same voltage phase.

t_shift = @(a,t0,fs) ifft(fft(a).*exp(-1i*2*pi*t0*(ifftshift((0:length(a)-1)' -ceil((length(a)-1)/2))/length(a)*fs)));
% fftshift and ifftshift switch the left and right half in order to put the zero frequency at the center.
% Makes makes the first and last 50 points weird.

theta_shift = @(a,theta) ifft(fft(a).*exp(-1i*theta*(ifftshift((0:length(a)-1)' -ceil((length(a)-1)/2)))));

AIN0 = t_shift(AIN0,0*dt,fs);
AIN1 = t_shift(AIN1,1*dt,fs);
AIN2 = t_shift(AIN2,2*dt,fs);
AIN3 = t_shift(AIN3,3*dt,fs);
AIN4 = t_shift(AIN4,4*dt,fs);
AIN5 = t_shift(AIN5,5*dt,fs);

if exist('Speed','var') % checks if variables are defined
    Speed = t_shift(Speed,dtS,fs);
end

%% Important things for PM Motors only
%{
if l2l_flag %inverter connected wiring
    Ia = AIN0;
    Ib = AIN1;
    Ic = AIN2;
    Vbc = real(AIN4);
    Vca = real(AIN5);
    Vab = real(AIN5);

    % need to convert form Line to Line voltages to Line to Neutral Voltages
    % following laughman Thesis.
    delta2wye_model = [1 -1 0; 0 1 -1; -1 0 1; 1 1 1];
    Vdelta = [Vab Vbc Vca zeros(size(Vab))];
    Vwye = zeros(length(Vdelta),3);

    delta2wye_Operator = (delta2wye_model'*delta2wye_model)\delta2wye_model';

    for k = 1:length(Vdelta)
        Vwye(k,:) = (delta2wye_Operator*Vdelta(k,:)')';
    end
    SlopesAAF = [0.0916399, 0.0910994, 0.0911048];
    out.V_l2l = [Vab Vbc Vca].*mean(SlopesAAF);
    Va = Vwye(:,1).*SlopesAAF(1);
    Vb = Vwye(:,2).*SlopesAAF(2);
    Vc = Vwye(:,3).*SlopesAAF(3);
else
%}

%% Map Labjack inputs.

if mill=='light duty';
    Ia = AIN0; % For early heavy duty data this is negative.
    Ib = AIN1;
    Ic = AIN2;
    Va = AIN4; % These are remapped due to configuration differences in the NILM box.
    Vb = AIN3;
    Vc = AIN5;
else
    % Post Safety Modifications
%     Ia = AIN0;
%     Ib = -AIN1; % current sensor are on backwards so this is negative
%     Ic = -AIN2; % current sensor are on backwards so this is negative
%     Va = AIN3; % These are remapped due to configuration differences in the NILM box.
%     Vb = AIN4;
%     Vc = AIN5;
    % Pre Safety Modifications
    Ia = AIN1;
    Ib = AIN0; % These are remapped due to configuration differences in the NILM box.
    Ic = AIN2; 
    Va = AIN3; 
    Vb = AIN4;
    Vc = AIN5;
end

%%

out.Speed = real(Speed);
out.Power = real(Va.*Ia + Vb.*Ib + Vc.*Ic);
out.Vwye = real([Va Vb Vc]);
out.Iabc = real([Ia Ib Ic]);

L = length(out.Power);
T = 1/fs;

out.Time = (0:L-1)'*T;
out.T = T;
out.fs = fs;

%% Complex Clarke Transform

a = exp(1i*2*pi/3);

out.V_clark = sqrt(2/3)*(Va + a.*Vb + a^2.*Vc); % Imaginary component is a real number, but stored there for convenience
out.I_clark = sqrt(2/3)*(Ia + a.*Ib + a^2.*Ic);

out.V_gamma = sqrt(2/3)*(1/2)*(Va + Vb + Vc);   % Not used, computed for completeness
out.I_gamma = sqrt(2/3)*(1/2)*(Ia + Ib + Ic);


%% Park Transform

if l2l_flag  % inverter-driven, so very noisy
    intercept_times = Zero_Intercept_Times(real(resample(out.V_clark,1,10)),out.Time(1:10:end)');
else
    intercept_times = Zero_Intercept_Times(real(out.V_clark),out.Time');
end

diff_time = intercept_times(2:end)-intercept_times(1:end-1);

out.we_mean = pi/mean(diff_time);

mixer = exp(-1i*out.we_mean*out.Time);
V_mixed = out.V_clark.*mixer;

% Line below loads precomputed FIR filter coefficients form file much faster than recomputing them each time.
% Contents are in vector called 'b'.
load('LP_FIR_filter_coeff.mat');
N = length(b);

V_mixed_appended = [mean(V_mixed(1:100))*ones(N,1); V_mixed; mean(V_mixed(end-100+1:end))*ones(N,1)];
temp = filtfilt(b,1,V_mixed_appended);

V_filt = temp(N+1:end-N);

t = (0:L-1)'*T;

phase = unwrap(angle(V_filt))+out.we_mean*t;
park_vect = exp(-1i*phase);

Idq = out.I_clark.*park_vect; % Current output from the Park Transform

%%%%%%%%%% DEBUGGING %%%%%%%%%%
Idq = complex(-imag(Idq),real(Idq));

Ids_mean = mean(real(Idq));
Iqs_mean = mean(imag(Idq));

%%%%%%%%%% DEBUGGING %%%%%%%%%% (changed from - to + pi/4)
phi_offset = atan2(Ids_mean, Iqs_mean) + pi/4;   % Rotate currents to have same mean positive values.
% phi_offset = atan2(Ids_mean, Iqs_mean);        % Rotate Park Transform to make D-axis have mean zero value.

park_vect_rot = park_vect.*exp(1i*phi_offset);

out.VdqDeMod = out.V_clark.*park_vect_rot; % Ideal non-zero crossing Park Transform output
out.IdqDeMod = out.I_clark.*park_vect_rot;

out.V_zero = sqrt(3/2)*(sqrt(2)/2)*(Va + Vb + Vc); % Not used, computed for completeness
out.I_zero = sqrt(3/2)*(sqrt(2)/2)*(Ia + Ib + Ic);

out.phase = phase-phi_offset;

temp = (fs/12)*filter([-1 8 0 -8 1],1,out.phase);
out.WeDeMod = [mean(temp(5:50))*ones(4,1);temp(5:end)];

if isfield(out, 'Speed')
    phase_rotor = cumtrapz(out.Time,P*out.Speed);
    park_vect_rotor = exp(-1i*phase_rotor);
    Idq_rotor = out.I_clark.*park_vect_rotor;
    Ids_mean = mean(real(Idq_rotor));
    Iqs_mean = mean(imag(Idq_rotor));
    phi_offset = atan2(Ids_mean, Iqs_mean) + pi/4;   % Rotate currents to have same mean positive values.
    % phi_offset = atan2(Ids_mean, Iqs_mean);        % Rotate Park Transform to make D-axis have zero value.
    park_vect_rotor_shift = park_vect_rotor.*exp(1i*phi_offset);
    out.Vdq_rotor = out.V_clark.*park_vect_rotor_shift;
    out.Idq_rotor = out.I_clark.*park_vect_rotor_shift;
    out.phase_rotor = -1*unwrap(angle(park_vect_rotor_shift));
end

out.P = P;

end
