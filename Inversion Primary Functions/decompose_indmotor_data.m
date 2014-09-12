function out = decompose_indmotor_data(data, load_wr, focus_ind)

out = data;
fs = data.fs;

% If load_wr = -1, then the load periodicity is equal to the shaft rotational speed.
if nargin < 2
    load_wr = -1;
    focus_ind = data.indl_schantz;
elseif nargin <3;
    focus_ind = data.indl_schantz;
end

mI = mean(out.IdqDeMod);

% Usually the highest energy is in the line synchronous harmonics.
out.I_line_sync = reference_angle_domain_average(out.IdqDeMod(focus_ind), out.WeDeMod(focus_ind), fs);
out.I_remainder1 = out.IdqDeMod(focus_ind(1):focus_ind(end-1))-(out.I_line_sync);

% Next do the load synchronous harmonics.
if load_wr ~= -1;
    out.I_load_sync = direct_sychronous_average_new(out.I_remainder1, load_wr, fs)+mI;
end

if isfield(out,'speed_estimate')
    out.I_rot_sync = reference_angle_domain_average(out.IdqDeMod(focus_ind), out.speed_estimate, fs);
else
    out.I_rot_sync = direct_sychronous_average_new(out.IdqDeMod(focus_ind(1):focus_ind(end-1)), out.wr, fs);
    out.V_rot_sync = direct_sychronous_average_new(out.VdqDeMod(focus_ind(1):focus_ind(end-1)), out.wr, fs);
end

% For strong rotational periodicities, its possible to track changing speed with phase of the current harmonic at 1x shaft speed.  
% Code for this with zero crossings is in older methods.
% out.I_rot_sync = reference_angle_domain_average(out.IdqDeMod(focus_ind), out.slip_speed_est, fs);

if isfield(out,'Speed');
  
    out.Speed_rot_synch_average = direct_sychronous_average_new(out.Speed(focus_ind(1):focus_ind(end-1)), out.wr, fs);
    out.Speed_line_synch_average = reference_angle_domain_average(out.Speed(focus_ind), out.WeDeMod(focus_ind), fs);
    
    if load_wr ~= -1;
        out.Speed_load_synch_average = direct_sychronous_average_new(out.Speed(focus_ind(1):focus_ind(end-1)), load_wr, fs);
        % Speed_average = out.Speed_rot_synch_average+out.Speed_line_synch_average+out.Speed_load_synch_average;
        Speed_average = out.Speed_rot_synch_average+out.Speed_load_synch_average;
        out.Speed_average = Speed_average - mean(Speed_average)/2;
    else
        Speed_average = out.Speed_rot_synch_average+out.Speed_line_synch_average;
        out.Speed_average = Speed_average - mean(Speed_average)/2;
    end
    
end

end

