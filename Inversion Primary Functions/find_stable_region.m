function out = find_stable_region(data, auto_flag, num_sec)

out = data;
clear data;
time = out.Time;

if nargin < 2
    auto_flag = false;
end

if nargin < 3
    num_sec = 5;
end

%%
if auto_flag % PM Motor
    
    window_length = round(num_sec*out.fs);
    rms_we_variation = running_rms(out.WeDeMod-out.we_mean,window_length);
    
    [~, iind_center] = min(rms_we_variation);
    time_ind = iind_center:(iind_center + window_length-1);
    
    % h = figure;
    % plot(time,out.WeDeMod/2/pi,time(time_ind),out.WeDeMod(time_ind)/2/pi,'r')

else % Induction Motor

    for i = 1:length(out.VdqDeMod);
        p_inst(i) = real(out.VdqDeMod(i))*real(out.IdqDeMod(i)) + imag(out.VdqDeMod(i))*imag(out.IdqDeMod(i));
        q_inst(i) = imag(out.VdqDeMod(i))*real(out.IdqDeMod(i)) - real(out.VdqDeMod(i))*imag(out.IdqDeMod(i));
        app_inst(i) = sqrt((p_inst(i))^2+(q_inst(i))^2);
    end
    
    app_mean = mean(app_inst);
    window_length = round(num_sec*out.fs);
    rms_app_variation = running_rms(app_inst-app_mean,window_length);
    [~, iind_center] = min(rms_app_variation);
    time_ind = iind_center:(iind_center + window_length-1);
    %{
    h = figure;
    plot(time, out.WeDeMod/2/pi);
    title('Drive Frequency vs Time.  Please select left then right bounds of the most stable few seconds');
    x = ginput(2);
    if (x(1,1) > x(2,1))
        t = x(2,1);
        x(2,1) = x(1,1);
        x(1,1) = t;
    end
    
    ind_estimate = (time > x(1)) & (time < x(2));
    time = time(ind_estimate);
    time_ind = find((out.Time > time(1)) & (out.Time < time(end)));
    %hold on
    %plot(out.Time(time_ind),out.WeDeMod(time_ind)/2/pi,'r')
    %} 
end

%%
[~, sI, eI] = cycle_mean(real(out.V_clark(time_ind)));  % find indices for whole number of cycles

out.indl_schantz = time_ind(1)-1 + (sI:eI);

% PM motors:
% out.wr = mean(out.WeDeMod(out.indl_schantz)/out.P);

out.we = mean(out.WeDeMod(out.indl_schantz));

% Induction motors:
out.wr = mean(out.Speed(out.indl_schantz));

%pause
%close(h)

end
