function [out single] = reference_angle_domain_average(data, wr_reference, fs)
% this interpolates the data vector into an equal elapsed angle per sample
% space using speed reference angle vector
% performs sychronous averaging in this space, then returns data to
% origional space.
% due to interpolation limits the output is one data point shorter than the
% input.

mw = mean(wr_reference);

L = length(data);

rps = mw/(2*pi);
Period = 1/rps;

new_fs = ceil(Period*fs)/Period;

spp = round(Period*new_fs);

Time = (0:L-1)/fs;

rot_angle_est = cumtrapz(Time,wr_reference);

angle_interp = rot_angle_est(1):(2*pi)/spp:rot_angle_est(end);

max_length = 1e5;

if L <= max_length;
    
    data_equal_angle_interp = interp1(rot_angle_est, data, angle_interp);
    time_equal_angle_interp = interp1(rot_angle_est, Time, angle_interp);
    
else
    
    data_equal_angle_interp = [];
    time_equal_angle_interp = [];
    
    for i = 1:ceil(L/max_length);
        
        indS = (i-1)*max_length+1;
        indE = indS+max_length;
        indE = min(indE, L);
        
        ind_angle_S = find(angle_interp >= rot_angle_est(indS),1,'first');
        ind_angle_E = find(angle_interp < rot_angle_est(indE),1,'last');
        
        data_equal_angle_interp = [data_equal_angle_interp interp1(rot_angle_est(indS:indE), data(indS:indE), angle_interp(ind_angle_S:ind_angle_E))];
        time_equal_angle_interp = [time_equal_angle_interp interp1(rot_angle_est(indS:indE), Time(indS:indE), angle_interp(ind_angle_S:ind_angle_E))];
        
    end
end

[avg_data single] = sychronous_average(data_equal_angle_interp);

L2 = length(avg_data);

if L2 <= max_length;
    out = interp1(time_equal_angle_interp, avg_data, Time);
else
    out = [];
    
    for i = 1:ceil(L2/max_length);
        
        indS = (i-1)*max_length+1;
        indE = indS+max_length;
        indE = min(indE, L2);
        
        ind_time_S = find(Time >= time_equal_angle_interp(indS),1,'first');
        ind_time_E = find(Time < time_equal_angle_interp(indE),1,'last');
        if indE == L2;
            ind_time_E = length(Time);
        end
        out = [out interp1(time_equal_angle_interp(indS:indE), avg_data(indS:indE), Time(ind_time_S:ind_time_E))];
                
    end
end

out = reshape(out,size(data));
out(end) = [];

    function [out h] = sychronous_average(meas)
        C = floor(length(meas)/spp);
        h = reshape(meas(1:C*spp),spp,C);
        h = mean(h,2);
        out = repmat(h,C,1);
        Lm = length(meas);
        out = [out; out(1:(Lm-(C*spp)))];
    end

end
