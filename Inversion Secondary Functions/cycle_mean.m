function [cm start_I end_I] = cycle_mean(dat)
% Finds mean between first maximum and last minimum for a rough estimate of
% the cycle mean.  Simply taking mean of entire signal can be misleading if
% signal is not composed of integer number of periods.

Y = fft(dat-mean(dat));
[~, ii] = max(abs(Y));
L = length(dat);

dL = round(L/1000/2);

cut_window = round(10*L/ii);  

% fft_based_filtering:
tY = zeros(size(Y));
tY(ii-dL:ii+dL) = Y(ii-dL:ii+dL);
tY(end-ii+2-dL:end-ii+2+dL) = Y(end-ii+2-dL:end-ii+2+dL);
clean = ifft(tY);

% need to remove 10 cycle window lengts to reduce edge effects

cdat = clean(cut_window:end-cut_window+1);

% we want to find the first and last positive going zero crossing in the
% dat vector:

mask_up = cdat(1:end-1) < 0 & cdat(2:end) > 0;  % flag all positive going zero crissings;

start_I_temp = find(mask_up == 1,1,'first');
end_I_temp = find(mask_up == 1,1,'last');

temp = 1:length(cdat);

if abs(cdat(start_I_temp)) <= abs(cdat(start_I_temp+1))  % first point closer to zero
    start_I = start_I_temp+cut_window-1;
else
    start_I = start_I_temp+1+cut_window-1;
end

if abs(cdat(end_I_temp)) <= abs(cdat(end_I_temp+1))  % first point closer to zero
    end_I = end_I_temp+cut_window-1;
else
    end_I = end_I_temp+1+cut_window-1;
end

% figure
% plot(temp,[dat(cut_window:end-cut_window+1), cdat],temp(mask_up),cdat(mask_up),'r*',temp([start_I end_I]),cdat([start_I,end_I]),'g*')
% pause

cm = mean(dat(start_I:end_I));

% start_I = start_I - cut_window;
% end_I = end_I + cut_window;

end