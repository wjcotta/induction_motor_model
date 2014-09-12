format short
format compact
rms = @(a)norm(a)./sqrt(length(a));
time_me = @(x,fs)(0:length(x)-1)/fs;
freq_me_w = @(a)ifftshift((0:length(a)-1)' -ceil((length(a)-1)/2));
freq_me_fs = @(a, fs) freq_me_w(a)./length(a)*fs;
xc = @(A,B) ifft(fft(A).*conj(fft(B)));

% Change current directory if still at default. 
if strcmp(pwd,'D:\Joshua\Documents\MATLAB')
     cd 'D:\Joshua\Dropbox\Motor NILM\Induction Motor Model'
end