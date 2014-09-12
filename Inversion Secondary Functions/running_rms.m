function out = running_rms(in, L, detrend_flag)

if nargin < 3
    detrend_flag = false;
end

[num_row num_col] = size(in);

if (num_row == 1) || (num_col == 1)

if detrend_flag
    in = detrend(in);
end

if size(in,1) > size(in,2)
   temp = [0; cumsum(in.^2)];
else
   temp = [0 cumsum(in.^2)];
end

% L is length in data points of window.

out = sqrt((1/L)*(temp((L+1):end)-temp(1:end-(L+1)+1)));

else
    % assume a stack of row vectors
    if detrend_flag
        in = detrend(in')';
    end
    temp = [zeros(num_row,1), cumsum(in.^2,2)];
    out = sqrt((1/L)*(temp(:,(L+1):end)-temp(:,1:end-(L+1)+1)));
end