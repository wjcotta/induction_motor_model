function out = preprocess_sim_data(filename, fs, l2l_flag, P, encoder_count, mill)

% This function converts the raw output of the Labjack DAQ device into physical units suitable for calculation, does time-shifting to account for interchannel delay, 
% performs Clark and Park transformations (accounting for variable electric frequency over course of experiment), and produces an estimate of instantaneous electric frequency.

if (nargin < 5)
    encoder_count = 1800;
end

%% Load and unpack data

load(filename)
out.filename = filename;

Va = data(:,1);
Vb = data(:,2);
Vc = data(:,3);
Ia = data(:,4);
Ib = data(:,5);
Ic = data(:,6);
Speed = data(:,7);

%% Assign outputs

out.Speed = Speed;
out.Power = Va.*Ia + Vb.*Ib + Vc.*Ic;
out.Vwye = [Va Vb Vc];
out.Iabc = [Ia Ib Ic];

L = length(out.Power);
T = 1/fs;

out.Time = (0:L-1)'*T;
out.T = T;
out.fs = fs;

%% Complex Clarke Transform

a = exp(1i*2*pi/3);

%%%%%%%%%% DEBUGGING %%%%%%%%%% (switched b and c)
out.V_clark = sqrt(2/3)*(Va + a.*Vb + a^2.*Vc); % Imaginary component is a real number, but stored there for convenience
out.I_clark = sqrt(2/3)*(Ia + a.*Ib + a^2.*Ic);

% Alternate formula for Clarke Transformation used for verification
out.V_alpha = sqrt(2/3)*(Va - (1/2)*Vb - (1/2)*Vc);
out.V_beta = sqrt(2/3)*(0 + sqrt(3)/2*Vb - sqrt(3)/2*Vc);

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
out.Idq = Idq; 

Ids_mean = mean(real(Idq));
Iqs_mean = mean(imag(Idq));

%%%%%%%%%% DEBUGGING %%%%%%%%%% (changed from - to + pi/4)
phi_offset = atan2(Ids_mean, Iqs_mean) + pi/4;    % Rotate currents to have same mean positive values.
% phi_offset = atan2(Ids_mean, Iqs_mean);         % Rotate Park Transform to make D-axis have mean zero value.

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
    
    %%%%%%%%%% DEBUGGING %%%%%%%%%%
    out.Idq_rotor_check = Idq_rotor;
    
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