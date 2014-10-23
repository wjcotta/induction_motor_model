function out = speed_inversion_error_space_sim_fixed(data, wr_load, P, tune, useTrueSpeed)
out = data;
fs = out.fs;
we = mean(out.WeDeMod(out.indl_schantz));

% Starting point for parameters

Rs = 4;
Rr = 4;
Lm =  0.15;
Ls = 0.006;
Lr = 0.02;

Lal = mean([Ls Lr]);

% % Laughman Thesis Parameters
% %Cold Motor Fit to startup transient current:
% Rs = 4.96;
% Rr = 3.24;
% Lm = 53.67/we;  %0.142
% Ls = 1.64/we;   %0.00435
% Lr = 9.23/we;   %0.0245
% %Hot Motor Fit to startup transient current:
% Rs = 6.25;
% Rr = 4.03;
% Lm = 57.75/we;   %0.153
% Ls = 3.14/we;    %0.00833
% Lr = 7.71/we;    %0.0205
% %Cold Motor Fit to startup transient current and torque speed curve:
% Rs = 3.94;
% Rr = 4.12;
% Lm = 56.32/we;   %0.150
% Ls = 1.71/we;    %0.00454
% Lr = 8.86/we;    %0.0235
% %Hot Motor Fit to startup transient current and torque speed curve:
% Rs = 5.24;
% Rr = 4.76;
% Lm = 57.84/we;   %0.153
% Ls = 2.86/we;    %0.00759
% Lr = 7.74/we;    %0.0205

Las = Ls + Lm;
Lar = Lr + Lm;
D = Lm*Lm - Las*Lar;

rms = @(x)norm(x)/sqrt(length(x));

% if wr_load ~= -1 then do no interpolation and do joint Rot and Load inversion.
% if wr_load = -1, then do interpolation and Rot inversion

if wr_load == -1
    
    I_temp = out.I_rot_sync;
    V_temp = out.VdqDeMod(out.indl_schantz);
    Speed_temp = out.Speed_rot_synch_average;
    
    L = length(I_temp);
    mI = mean(I_temp);
    mV = mean(V_temp);
    mS = mean(Speed_temp);
    
    rps = out.wr/(2*pi);
    Period = 1/rps;
    
    new_fs = ceil(Period*fs)/Period
    spp = round(Period*new_fs);
    new_L = round(L/fs*new_fs);
    
    fs = new_fs;
    
    I_long = interpft(I_temp-mI,new_L)+mI;
    V_long = interpft(V_temp-mV,new_L)+mV;
    Speed_long = interpft(Speed_temp-mS,new_L)+mS;
    
    Nrot = spp;
else
    I_long = out.I_rot_sync;
    V_long = out.VdqDeMod(out.indl_schantz);
    Speed_long = out.Speed_average;
    
    period_rot = 2*pi/out.wr;
    period_load = 2*pi/wr_load;
    
    Nrot = round(period_rot*fs);
    Nload = round(period_load*fs);
    
end

Wp = 600/(fs/2);
Ws = 850/(fs/2);
[n Wn] = buttord(Wp, Ws, 1, 10);  % filter to remove stator bar harmonics
[b a] = butter(n, Wn);

I_rot_syncF = filtfilt(b,a,I_long);
Idq_rot = I_rot_syncF(1e4:1e4+Nrot-1);

out.Idq_rot_Averaged = I_long(1e4:1e4+Nrot-1);
out.new_fs = fs;

Vdq_rot = mean(V_long).*ones(Nrot,1);

Idq_sim_compare = Idq_rot;

if wr_load ~= -1
    I_load_syncF = filtfilt(b,a,out.I_load_sync);
    Idq_load = I_load_syncF(1e4:1e4+Nload-1);
    Vdq_load = mean(out.VdqDeMod).*ones(Nload,1);
    
    
    if Nload > Nrot
        s2 = floor(Nload/Nrot);
        r2 = Nload - s2*Nrot;
        
        Idq_sim_compare = Idq_load + [repmat(Idq_rot,s2,1); Idq_rot(1:r2)];
        Idq_sim_compare = Idq_sim_compare - mean(Idq_sim_compare)/2;
    else
        s2 = floor(Nrot/Nload);
        r2 = Nrot - s2*Nload;
        
        Idq_sim_compare = Idq_rot + [repmat(Idq_load,s2,1); Idq_load(1:r2)];
        Idq_sim_compare = Idq_sim_compare - mean(Idq_sim_compare)/2;
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find best parameters 4 param model.
Initial_Parameters = [1 1 1 1];
Scaling_Factors = [Rs Rr Lm Lal];

opt = optimset('MaxFunEvals',5e4,'Display','iter','MaxIter',100, 'TolX',1e-10, 'TolFun',1e-10, 'Jacobian','on');


if wr_load == -1
    Parameters = lsqnonlin(@error_function_mean_speed_with_Jacobian_Rot, Initial_Parameters,[.1 .1 .1 .1], [10 10 10 10], opt);
else
    Parameters = lsqnonlin(@error_function_mean_speed_with_Jacobian_Rot_and_Load, Initial_Parameters,[.1 .1 .1 .1], [10 10 10 10], opt);
end
%{
% FIND_FLAG
Parameters = [2.4819 1.5868 0.1610 0.0080];
Parameters = Parameters./Scaling_Factors;
%}
Rs = Parameters(1)*Scaling_Factors(1);
Rr = Parameters(2)*Scaling_Factors(2);
Lm = Parameters(3)*Scaling_Factors(3);
Lal = Parameters(4)*Scaling_Factors(4);

Las = Lm + Lal;
Lar = Lm + Lal;

D = Lm^2 - Las*Lar;

if wr_load == -1
    [error, junk, Wr_est, sim_state] = error_function_mean_speed_with_Jacobian_Rot(Parameters);
else
    [error, junk, Wr_est, sim_state] = error_function_mean_speed_with_Jacobian_Rot_and_Load(Parameters);
end

out.Wr_est = Wr_est;

out.Parameters = [Rs, Rr, Lm, Lal];

out.error = error;

%%% Find best parameters 5 param model.
% Initial_Parameters = [1 1 1 1 1];
% Scaling_Factors = [Rs Rr Lm Ls Lr];
%
% opt = optimset('MaxFunEvals',5e4,'Display','iter','MaxIter',100, 'TolX',1e-10, 'TolFun',1e-10, 'Jacobian','on');
%
% if wr_load == -1
%     Parameters = lsqnonlin(@error_function_mean_speed_with_Jacobian_Rot_5_param, Initial_Parameters,[.1 .1 .1 .1 .1], [10 10 10 10 10], opt);
% else
% %     Parameters = lsqnonlin(@error_function_mean_speed_with_Jacobian_Rot_and_Load, Initial_Parameters,[.1 .1 .1 .1 .1], [10 10 10 10 10], opt);
% end
%
% Rs = Parameters(1)*Scaling_Factors(1)
% Rr = Parameters(2)*Scaling_Factors(2)
% Lm = Parameters(3)*Scaling_Factors(3)
% Ls = Parameters(4)*Scaling_Factors(4)
% Lr = Parameters(5)*Scaling_Factors(5)
%
% Las = Lm + Ls;
% Lar = Lm + Lr;

% D = Lm^2 - Las*Lar;
%
% if wr_load == -1
%     [error, junk, Wr_est, sim_state] = error_function_mean_speed_with_Jacobian_Rot_5_param(Parameters);
% else
% %     [error, junk, Wr_est] = error_function_mean_speed_with_Jacobian_Rot_and_Load(Parameters);
% end

% out.Wr_est = Wr_est;
%
% out.Parameters = [Rs, Rr, Lm, Ls, Lr];
%
% out.error = error;

if wr_load == -1
    RotTime = (0:Nrot-1)/fs;
    sim_point_count = Nrot;
else
    RotTime = (0:max([Nrot Nload])-1)/fs;
    sim_point_count = max([Nrot Nload]);
end

out.RotTime = RotTime;

sim_time_vector = 0:1/fs:(sim_point_count-1)/fs;

[t y] = ode45(@induction_motor_parameter_fitting_model_arbitrary_load, sim_time_vector,sim_state);

% lamqs = y(end-Nrot+1:end,1);
% lamds = y(end-Nrot+1:end,2);
% lamqr = y(end-Nrot+1:end,3);
% lamdr = y(end-Nrot+1:end,4);

lamqs = y(:,1);
lamds = y(:,2);
lamqr = y(:,3);
lamdr = y(:,4);
out.lambdas = y;






Idq_Estimated = ((Lm*lamdr - Lar*lamds) + 1i*(Lm*lamqr - Lar*lamqs))/D;

out.Idq_Verification = Idq_Estimated;

N = length(out.Wr_est);
% close all;
h = figure;
set(h,'Outerposition',[100,100,1200,700])

uicontrol('Style', 'text', 'String', data.filename, 'Units','normalized', 'Position', [0.45 0.975 0.15 0.025]); 

subplot(3,2,1);
A = real(out.Wr_est)/2/pi;
B = Speed_long(1e4:1e4+N-1)/2/pi;

out.Speed_verification = B;

xc = ifft(fft(A).*conj(fft(B)));
[junk, shift] = max(xc);
% shift = 0;
B = circshift(B,shift);
a1 = plot([A B]);
legend('Estimated','Real')

wr_mean = out.wr;
dd = mean(real(out.Wr_est))-wr_mean;
title(['Mean Speed Diff: ' num2str(dd) ' rps, RMS imag part: ' num2str(rms(imag(out.Wr_est))) '.'])
subplot(3,2,3);
plot(imag(out.Wr_est)/2/pi);
subplot(3,2,5);
plot(out.error);

subplot(2,2,2)
plot(RotTime,real(Idq_sim_compare),'b',RotTime,real(Idq_Estimated),'r')
title('I_d_s Stator Current Comparison');
legend('Measured','Estimated')
xlabel('Time (sec)');
ylabel('Current (amps)');


subplot(2,2,4)
plot(RotTime,imag(Idq_sim_compare),'b',RotTime,imag(Idq_Estimated),'r')
title('I_q_s Stator Current Comparison');
legend('Measured','Estimated')
xlabel('Time (sec)');
ylabel('Current (amps)');
figure(gcf)


%%
    function [combNum, combDem] = recombine(N1, D1, N2, D2, N3, D3)
        if nargin == 2
            combNum = N1;
            combDem = D1;
            return
        elseif nargin == 4;
            L1 = length(N1);
            L2 = length(N2);
            if L1 < L2  % make loongest in N1.
                T = N1;
                N1 = N2;
                N2 = T;
                T = D1;
                D1 = D2;
                D2 = T;
                L1 = length(N1);
                L2 = length(N2);
            end
            s2 = floor(L1/L2);
            r2 = L1 - s2*L2;
            combNum = N1 + [repmat(N2,s2,1); N2(1:r2)];
            combDem = D1 + [repmat(D2,s2,1); D2(1:r2)];
            combNum = combNum - mean(combNum)/2;
            combDem = combDem - mean(combDem)/2;
            return
        elseif nargin == 6; % assuming there will be no equal lengths
            LL = [length(N1); length(N2); length(N3)];
            cN = {N1, N2, N3};
            cD = {D1, D2, D3};
            [B IX] = sort(LL,'descend');
            cN = cN(IX);
            cD = cD(IX);
            LL = LL(IX);
            %%% sorted
            s2 = floor(LL(1)/LL(2));
            r2 = LL(1) - s1*LL(2);
            s3 = floor(LL(1)/LL(3));
            r3 = LL(1) - s1*LL(3);
            combNum = cN{1} + [repmat(cN{2},s2,1); cN{2}(1:r2)] + [repmat(cN{3},s3,1); cN{3}(1:r3)];
            combDem = cD{1} + [repmat(cD{2},s2,1); cD{2}(1:r2)] + [repmat(cD{3},s3,1); cD{3}(1:r3)];
            combNum = combNum - 2*mean(combNum)/3;
            combDem = combDem - 2*mean(combDem)/3;
            return
        end
    end

    function [Num, Dem, VsK, IsK, fK, Lambda_sK, Lambda_rK, IrK, k, Num_m, Dem_m, IsK_m, fK_m, Lambda_sK_m, Lambda_rK_m, IrK_m] = partial_fft_invert(Is, Vs)
        N = length(Is);
        Time = N/fs;
        
        wo = 2*pi/Time;
        if mod(N,2) == 0
            k = wo*[0:N/2 - 1, 0, -N/2  + 1:-1]';
        else
            k = wo*[0:(N-1)/2, -(N-1)/2:-1]';
        end
        
        VsK = fft(Vs);
        IsK = fft(Is);
        if mod(N,2) == 0
            VsK(N/2 + 1) = 0;
            IsK(N/2 + 1) = 0;
        end
        
        fK = VsK - Rs*IsK;
        Lambda_sK = fK./(1i*(k + we));
        Lambda_rK = (D*IsK + Lar*Lambda_sK)/Lm;
        IrK = (Lm*Lambda_sK - Las*Lambda_rK)/D;
        Num = ifft((1i*k).*Lambda_rK + Rr*IrK);
        Dem = ifft(1i*Lambda_rK);
        
        
        IsK_m = fft(mean(Is)*ones(size(Is)));
        
        fK_m = VsK - Rs*IsK_m;
        Lambda_sK_m = fK./(1i*(k + we));
        Lambda_rK_m = (D*IsK_m + Lar*Lambda_sK_m)/Lm;
        IrK_m = (Lm*Lambda_sK_m - Las*Lambda_rK_m)/D;
        Num_m = ifft((1i*k).*Lambda_rK_m + Rr*IrK_m);
        Dem_m = ifft(1i*Lambda_rK_m);
        
        
        
    end


% Nested function that computes the objective function

    function error = error_function_mean_speed_with_Jacobian_Rot_SA(Parameters)
        error = error_function_mean_speed_with_Jacobian_Rot(Parameters);
        error = sum(error);
    end

    function [error, J, Wr_est, sim_state] = error_function_mean_speed_with_Jacobian_Rot(Parameters)
        Rs = Parameters(1)*Scaling_Factors(1);
        Rr = Parameters(2)*Scaling_Factors(2);
        Lm = Parameters(3)*Scaling_Factors(3);
        Lal = Parameters(4)*Scaling_Factors(4);
        
        Las = Lm + Lal;
        Lar = Lm + Lal;
        D = Lm^2 - Las*Lar;
        
        [Num, Dem, VsK, IsK, fK, Lambda_sK, Lambda_rK, IrK, k ...
            Num_m, Dem_m, IsK_m, fK_m, Lambda_sK_m, Lambda_rK_m, IrK_m] = partial_fft_invert(Idq_rot, Vdq_rot);
        
        slip = Num./Dem;
        slip_m = Num_m./Dem_m;
        
        Wr_est = (2/(2*P))*(slip + we);
        Wr_est_m = (2/(2*P))*(slip_m + we);
        wr_mean = out.wr;
        
        
        if ~useTrueSpeed
            error = tune*imag(Wr_est).^2 + (wr_mean-mean(real(Wr_est_m))).^2;
        else
            N = length(Wr_est);
            A = real(Wr_est)/2/pi;
            B = Speed_long(1e4:1e4+N-1)/2/pi;
            xc = ifft(fft(A).*conj(fft(B)));
            [junk, shift] = max(xc);
            shift = 0;
            B = circshift(B,shift);
            error = tune*imag(Wr_est).^2 + (B-A).^2;
        end
        
        if nargout > 1
            
            dfK_dRs = -IsK;
            
            dlam_sK_dRs = dfK_dRs./(1i*(k+we));
            
            dlam_rK_dRs = Lar/Lm*dlam_sK_dRs;
            dlam_rK_dLm = (Lal^2*IsK - Lal*Lambda_sK)/(Lm^2);
            dlam_rK_dLal = -2*Lal*IsK/Lm - 2*IsK + Lambda_sK/Lm;
            
            dIrK_dRs = Lm*dlam_sK_dRs/D - Las*dlam_rK_dRs/D;
            dIrK_dLm = Lambda_rK/(-D) + ((-Lambda_sK*(-D) + Lm*Lambda_sK*2*Lal) + Lm*(dlam_rK_dLm*-D - Lambda_rK*2*Lal) + (Lal*dlam_rK_dLm*(2*Lm*Lal + Lal^2) - Lal*Lambda_rK*2*Lal))/((-D)^2);
            dIrK_dLal = (Lambda_rK)/(-D) + ((Lm*Lambda_sK*(2*Lm + 2*Lal)) + ((Lm*dlam_rK_dLal)*(2*Lm*Lal + Lal^2) - Lm*Lambda_rK*(2*Lm + 2*Lal)) + Lal*(dlam_rK_dLal*(-D) - Lambda_rK*(2*Lm + 2*Lal)))/((-D)^2);
            
            dslip_dRs = (ifft(1i*k.*dlam_rK_dRs + Rr*dIrK_dRs).*ifft(1i*Lambda_rK) - ifft(1i*k.*Lambda_rK + Rr*IrK).*ifft(1i*dlam_rK_dRs))./(ifft(1i*Lambda_rK).^2);
            dslip_dLm = (ifft(1i*k.*dlam_rK_dLm + Rr*dIrK_dLm).*ifft(1i*Lambda_rK) - ifft(1i*k.*Lambda_rK + Rr*IrK).*ifft(1i*dlam_rK_dLm))./(ifft(1i*Lambda_rK).^2);
            dslip_dLal = (ifft(1i*k.*dlam_rK_dLal + Rr*dIrK_dLal).*ifft(1i*Lambda_rK) - ifft(1i*k.*Lambda_rK + Rr*IrK).*ifft(1i*dlam_rK_dLal))./(ifft(1i*Lambda_rK).^2);
            dslip_dRr = ifft(IrK)./ifft(1i*Lambda_rK);
            
            dwr_dRs = (2/(2*P))*dslip_dRs;
            dwr_dRr = (2/(2*P))*dslip_dRr;
            dwr_dLm = (2/(2*P))*dslip_dLm;
            dwr_dLal = (2/(2*P))*dslip_dLal;
            
            %%%%%%%%%%%%%%%5
            dfK_dRs_m = -IsK_m;
            
            dlam_sK_dRs_m = dfK_dRs_m./(1i*(k+we));
            
            dlam_rK_dRs_m = Lar/Lm*dlam_sK_dRs_m;
            dlam_rK_dLm_m = (Lal^2*IsK_m - Lal*Lambda_sK_m)/(Lm^2);
            dlam_rK_dLal_m = -2*Lal*IsK_m/Lm - 2*IsK_m + Lambda_sK_m/Lm;
            
            dIrK_dRs_m = Lm*dlam_sK_dRs_m/D - Las*dlam_rK_dRs_m/D;
            dIrK_dLm_m = Lambda_rK_m/(-D) + ((-Lambda_sK_m*(-D) + Lm*Lambda_sK_m*2*Lal) + Lm*(dlam_rK_dLm_m*-D - Lambda_rK_m*2*Lal) + (Lal*dlam_rK_dLm_m*(2*Lm*Lal + Lal^2) - Lal*Lambda_rK_m*2*Lal))/((-D)^2);
            dIrK_dLal_m = (Lambda_rK_m)/(-D) + ((Lm*Lambda_sK_m*(2*Lm + 2*Lal)) + ((Lm*dlam_rK_dLal_m)*(2*Lm*Lal + Lal^2) - Lm*Lambda_rK_m*(2*Lm + 2*Lal)) + Lal*(dlam_rK_dLal_m*(-D) - Lambda_rK_m*(2*Lm + 2*Lal)))/((-D)^2);
            
            dslip_dRs_m = (ifft(1i*k.*dlam_rK_dRs_m + Rr*dIrK_dRs_m).*ifft(1i*Lambda_rK_m) - ifft(1i*k.*Lambda_rK_m + Rr*IrK_m).*ifft(1i*dlam_rK_dRs_m))./(ifft(1i*Lambda_rK_m).^2);
            dslip_dLm_m = (ifft(1i*k.*dlam_rK_dLm_m + Rr*dIrK_dLm_m).*ifft(1i*Lambda_rK_m) - ifft(1i*k.*Lambda_rK_m + Rr*IrK_m).*ifft(1i*dlam_rK_dLm_m))./(ifft(1i*Lambda_rK_m).^2);
            dslip_dLal_m = (ifft(1i*k.*dlam_rK_dLal_m + Rr*dIrK_dLal_m).*ifft(1i*Lambda_rK_m) - ifft(1i*k.*Lambda_rK_m + Rr*IrK_m).*ifft(1i*dlam_rK_dLal_m))./(ifft(1i*Lambda_rK_m).^2);
            dslip_dRr_m = ifft(IrK_m)./ifft(1i*Lambda_rK_m);
            
            dwr_dRs_m = (2/(2*P))*dslip_dRs_m;
            dwr_dRr_m = (2/(2*P))*dslip_dRr_m;
            dwr_dLm_m = (2/(2*P))*dslip_dLm_m;
            dwr_dLal_m = (2/(2*P))*dslip_dLal_m;
            %%%%%%%%%%%%%%%%%%5
            if ~useTrueSpeed
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dRs_m))))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dRr_m))))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dLm_m))))*Scaling_Factors(3);
                derror_dLal = (2*tune*imag(Wr_est).*imag(dwr_dLal) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dLal_m))))*Scaling_Factors(4);
            else
                N = length(Wr_est);
                B = Speed_long(1e4:1e4+N-1);
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(B - real(Wr_est)).*(-real(dwr_dRs)))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(B - real(Wr_est)).*(-real(dwr_dRr)))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(B - real(Wr_est)).*(-real(dwr_dLm)))*Scaling_Factors(3);
                derror_dLal = (2*tune*imag(Wr_est).*imag(dwr_dLal) + 2*(B - real(Wr_est)).*(-real(dwr_dLal)))*Scaling_Factors(4);
            end
            
            J = [derror_dRs derror_dRr derror_dLm derror_dLal];
        end
        
        if nargout > 3
            
            Lambda_s = ifft(Lambda_sK);
            Lambda_r = ifft(Lambda_rK);
            
            sim_state = [imag(Lambda_s(1)), real(Lambda_s(1)) imag(Lambda_r(1)), real(Lambda_r(1)), 0];
        end
        
        
    end

    function [error, J, Wr_est, sim_state] = error_function_mean_speed_with_Jacobian_Rot_and_Load(Parameters)
        Rs = Parameters(1)*Scaling_Factors(1);
        Rr = Parameters(2)*Scaling_Factors(2);
        Lm = Parameters(3)*Scaling_Factors(3);
        Lal = Parameters(4)*Scaling_Factors(4);
        
        Las = Lm + Lal;
        Lar = Lm + Lal;
        D = Lm^2 - Las*Lar;
        
        % don't actually care the difference between Rot and Load, just
        % make sure Load cariables have the longer vector
        
        if Nload > Nrot
            [NumRot, DemRot, VsRotK, IsRotK, fRotK, LambdaRot_sK, LambdaRot_rK, IrRotK, kRot] = partial_fft_invert(Idq_rot, Vdq_rot);
            [NumLoad, DemLoad, VsLoadK, IsLoadK, fLoadK, LambdaLoad_sK, LambdaLoad_rK, IrLoadK, kLoad] = partial_fft_invert(Idq_load, Vdq_load);
            s2 = floor(Nload/Nrot);
            r2 = Nload - s2*Nrot;
        else
            [NumRot, DemRot, VsRotK, IsRotK, fRotK, LambdaRot_sK, LambdaRot_rK, IrRotK, kRot] = partial_fft_invert(Idq_load, Vdq_load);
            [NumLoad, DemLoad, VsLoadK, IsLoadK, fLoadK, LambdaLoad_sK, LambdaLoad_rK, IrLoadK, kLoad] = partial_fft_invert(Idq_rot, Vdq_rot);
            T = Nload;
            Nload = Nrot;
            Nrot = T;
            s2 = floor(Nload/Nrot);
            r2 = Nload - s2*Nrot;
            
            T = Nload;
            Nload = Nrot;
            Nrot = T;
            
        end
        
        
        combNum = NumLoad + [repmat(NumRot,s2,1); NumRot(1:r2)];
        combDem = DemLoad + [repmat(DemRot,s2,1); DemRot(1:r2)];
        
        combNum = combNum - mean(combNum)/2;
        combDem = combDem - mean(combDem)/2;
        
        slip = combNum./combDem;
        
        Wr_est = (2/(2*P))*(slip + we);
        wr_mean = out.wr;
        
        
        if ~useTrueSpeed
            error = tune*imag(Wr_est).^2 + (wr_mean-mean(real(Wr_est))).^2;
        else
            N = length(Wr_est);
            A = real(Wr_est)/2/pi;
            B = Speed_long(1e4:1e4+N-1)/2/pi;
            xc = ifft(fft(A).*conj(fft(B)));
            [junk, shift] = max(xc);
            shift = 0;
            B = circshift(B,shift);
            
            error = tune*imag(Wr_est).^2 + (B-A).^2;
        end
        
        if nargout > 1
            
            dfRotK_dRs = -IsRotK;
            
            dlamRot_sK_dRs = dfRotK_dRs./(1i*(kRot+we));
            
            dlamRot_rK_dRs = Lar/Lm*dlamRot_sK_dRs;
            dlamRot_rK_dLm = (Lal^2*IsRotK - Lal*LambdaRot_sK)/(Lm^2);
            dlamRot_rK_dLal = -2*Lal*IsRotK/Lm - 2*IsRotK + LambdaRot_sK/Lm;
            
            dIrRotK_dRs = Lm*dlamRot_sK_dRs/D - Las*dlamRot_rK_dRs/D;
            dIrRotK_dLm = LambdaRot_rK/(-D) + ((-LambdaRot_sK*(-D) + Lm*LambdaRot_sK*2*Lal) + Lm*(dlamRot_rK_dLm*-D - LambdaRot_rK*2*Lal) + (Lal*dlamRot_rK_dLm*(2*Lm*Lal + Lal^2) - Lal*LambdaRot_rK*2*Lal))/((-D)^2);
            dIrRotK_dLal = (LambdaRot_rK)/(-D) + ((Lm*LambdaRot_sK*(2*Lm + 2*Lal)) + ((Lm*dlamRot_rK_dLal)*(2*Lm*Lal + Lal^2) - Lm*LambdaRot_rK*(2*Lm + 2*Lal)) + Lal*(dlamRot_rK_dLal*(-D) - LambdaRot_rK*(2*Lm + 2*Lal)))/((-D)^2);
            
            
            dfLoadK_dRs = -IsLoadK;
            
            dlamLoad_sK_dRs = dfLoadK_dRs./(1i*(kLoad+we));
            
            dlamLoad_rK_dRs = Lar/Lm*dlamLoad_sK_dRs;
            dlamLoad_rK_dLm = (Lal^2*IsLoadK - Lal*LambdaLoad_sK)/(Lm^2);
            dlamLoad_rK_dLal = -2*Lal*IsLoadK/Lm - 2*IsLoadK + LambdaLoad_sK/Lm;
            
            dIrLoadK_dRs = Lm*dlamLoad_sK_dRs/D - Las*dlamLoad_rK_dRs/D;
            dIrLoadK_dLm = LambdaLoad_rK/(-D) + ((-LambdaLoad_sK*(-D) + Lm*LambdaLoad_sK*2*Lal) + Lm*(dlamLoad_rK_dLm*-D - LambdaLoad_rK*2*Lal) + (Lal*dlamLoad_rK_dLm*(2*Lm*Lal + Lal^2) - Lal*LambdaLoad_rK*2*Lal))/((-D)^2);
            dIrLoadK_dLal = (LambdaLoad_rK)/(-D) + ((Lm*LambdaLoad_sK*(2*Lm + 2*Lal)) + ((Lm*dlamLoad_rK_dLal)*(2*Lm*Lal + Lal^2) - Lm*LambdaLoad_rK*(2*Lm + 2*Lal)) + Lal*(dlamLoad_rK_dLal*(-D) - LambdaLoad_rK*(2*Lm + 2*Lal)))/((-D)^2);
            
            Num_dRs_temp = ifft(1i*kRot.*dlamRot_rK_dRs + Rr*dIrRotK_dRs);
            Num_dRs = [repmat(Num_dRs_temp,s2,1); Num_dRs_temp(1:r2)] + ifft(1i*kLoad.*dlamLoad_rK_dRs + Rr*dIrLoadK_dRs);
            Num_dRs = Num_dRs - mean(Num_dRs)/2;
            
            Dem_dRs_temp = ifft(1i*dlamRot_rK_dRs);
            Dem_dRs = [repmat(Dem_dRs_temp,s2,1); Dem_dRs_temp(1:r2)] + ifft(1i*dlamLoad_rK_dRs);
            Dem_dRs = Dem_dRs - mean(Dem_dRs)/2;
            
            dslip_dRs = (Num_dRs.*combDem - combNum.*Dem_dRs)./(combDem.^2);
            
            Num_dLm_temp = ifft(1i*kRot.*dlamRot_rK_dLm + Rr*dIrRotK_dLm);
            Num_dLm = [repmat(Num_dLm_temp,s2,1); Num_dLm_temp(1:r2)] + ifft(1i*kLoad.*dlamLoad_rK_dLm + Rr*dIrLoadK_dLm);
            Num_dLm = Num_dLm - mean(Num_dLm)/2;
            
            Dem_dLm_temp = ifft(1i*dlamRot_rK_dLm);
            Dem_dLm = [repmat(Dem_dLm_temp,s2,1); Dem_dLm_temp(1:r2)] + ifft(1i*dlamLoad_rK_dLm);
            Dem_dLm = Dem_dLm - mean(Dem_dLm)/2;
            
            dslip_dLm = (Num_dLm.*combDem - combNum.*Dem_dLm)./(combDem.^2);
            
            Num_dLal_temp = ifft(1i*kRot.*dlamRot_rK_dLal + Rr*dIrRotK_dLal);
            Num_dLal = [repmat(Num_dLal_temp,s2,1); Num_dLal_temp(1:r2)] + ifft(1i*kLoad.*dlamLoad_rK_dLal + Rr*dIrLoadK_dLal);
            Num_dLal = Num_dLal - mean(Num_dLal)/2;
            
            Dem_dLal_temp = ifft(1i*dlamRot_rK_dLal);
            Dem_dLal = [repmat(Dem_dLal_temp,s2,1); Dem_dLal_temp(1:r2)] + ifft(1i*dlamLoad_rK_dLal);
            Dem_dLal = Dem_dLal - mean(Dem_dLal)/2;
            
            dslip_dLal = (Num_dLal.*combDem - combNum.*Dem_dLal)./(combDem.^2);
            
            Num_dRr_temp = ifft(IrRotK);
            Num_dRr = [repmat(Num_dRr_temp,s2,1); Num_dRr_temp(1:r2)] + ifft(IrLoadK);
            Num_dRr = Num_dRr - mean(Num_dRr)/2;
            
            dslip_dRr = Num_dRr./combDem;
            
            dwr_dRs = (2/(2*P))*dslip_dRs;
            dwr_dRr = (2/(2*P))*dslip_dRr;
            dwr_dLm = (2/(2*P))*dslip_dLm;
            dwr_dLal = (2/(2*P))*dslip_dLal;
            
            
            if ~useTrueSpeed
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(wr_mean - mean(real(Wr_est))).*(-mean(real(dwr_dRs))))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(wr_mean - mean(real(Wr_est))).*(-mean(real(dwr_dRr))))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(wr_mean - mean(real(Wr_est))).*(-mean(real(dwr_dLm))))*Scaling_Factors(3);
                derror_dLal = (2*tune*imag(Wr_est).*imag(dwr_dLal) + 2*(wr_mean - mean(real(Wr_est))).*(-mean(real(dwr_dLal))))*Scaling_Factors(4);
            else
                N = length(Wr_est);
                B = Speed_long(1e4:1e4+N-1);
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(B - real(Wr_est)).*(-real(dwr_dRs)))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(B - real(Wr_est)).*(-real(dwr_dRr)))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(B - real(Wr_est)).*(-real(dwr_dLm)))*Scaling_Factors(3);
                derror_dLal = (2*tune*imag(Wr_est).*imag(dwr_dLal) + 2*(B - real(Wr_est)).*(-real(dwr_dLal)))*Scaling_Factors(4);
            end
            
            
            J = [derror_dRs derror_dRr derror_dLm derror_dLal];
        end
        
        
        if nargout > 3
            
            LambdaRot_s = ifft(LambdaRot_sK);
            LambdaLoad_s = ifft(LambdaLoad_sK);
            
            LambdaRot_r = ifft(LambdaRot_rK);
            LambdaLoad_r = ifft(LambdaLoad_rK);
            
            Lambda_s = LambdaLoad_s + [repmat(LambdaRot_s,s2,1); LambdaRot_s(1:r2)];
            Lambda_r = LambdaLoad_r + [repmat(LambdaRot_r,s2,1); LambdaRot_r(1:r2)];
            
            Lambda_s = Lambda_s - mean(Lambda_s)/2;
            Lambda_r = Lambda_r - mean(Lambda_r)/2;
                                
            sim_state = [imag(Lambda_s(1)), real(Lambda_s(1)) imag(Lambda_r(1)), real(Lambda_r(1)), 0];
        end
        
    end

    function [error, J, Wr_est] = error_function_mean_speed_with_Jacobian_Rot_and_Load2(Parameters)
        Rs = Parameters(1)*Scaling_Factors(1);
        Rr = Parameters(2)*Scaling_Factors(2);
        Lm = Parameters(3)*Scaling_Factors(3);
        Lal = Parameters(4)*Scaling_Factors(4);
        
        Las = Lm + Lal;
        Lar = Lm + Lal;
        D = Lm^2 - Las*Lar;
        
        % don't actually care the difference between Rot and Load, just
        % make sure Load cariables have the longer vector
        
        if Nload > Nrot
            [NumRot, DemRot, VsRotK, IsRotK, fRotK, LambdaRot_sK, LambdaRot_rK, IrRotK, kRot] = partial_fft_invert(Idq_rot, Vdq_rot);
            [NumLoad, DemLoad, VsLoadK, IsLoadK, fLoadK, LambdaLoad_sK, LambdaLoad_rK, IrLoadK, kLoad] = partial_fft_invert(Idq_load, Vdq_load);
            s2 = floor(Nload/Nrot);
            r2 = Nload - s2*Nrot;
        else
            [NumRot, DemRot, VsRotK, IsRotK, fRotK, LambdaRot_sK, LambdaRot_rK, IrRotK, kRot] = partial_fft_invert(Idq_load, Vdq_load);
            [NumLoad, DemLoad, VsLoadK, IsLoadK, fLoadK, LambdaLoad_sK, LambdaLoad_rK, IrLoadK, kLoad] = partial_fft_invert(Idq_rot, Vdq_rot);
            T = Nload;
            Nload = Nrot;
            Nrot = T;
            s2 = floor(Nload/Nrot);
            r2 = Nload - s2*Nrot;
            
            T = Nload;
            Nload = Nrot;
            Nrot = T;
            
        end
        
        
        combNum = NumLoad + [repmat(NumRot,s2,1); NumRot(1:r2)];
        combDem = DemLoad + [repmat(DemRot,s2,1); DemRot(1:r2)];
        
        combNum = combNum - mean(combNum)/2;
        combDem = combDem - mean(combDem)/2;
        
        combSlip = combNum./combDem;
        
        RotSlip = NumRot./DemRot;
        LoadSlip = NumLoad./DemLoad;
        
        Wr_Rot_mean = mean(real((2/(2*P))*(RotSlip + we)));
        Wr_Load_mean = mean(real((2/(2*P))*(LoadSlip + we)));
        
        Wr_est = (2/(2*P))*(combSlip + we);
        wr_mean = out.wr;
        
        
        if ~useTrueSpeed
            error = tune*imag(Wr_est).^2 + (wr_mean-Wr_Rot_mean).^2 + (wr_mean-Wr_Load_mean).^2;
        else
            N = length(Wr_est);
            A = real(Wr_est)/2/pi;
            B = Speed_long(1e4:1e4+N-1)/2/pi;
            xc = ifft(fft(A).*conj(fft(B)));
            [junk, shift] = max(xc);
            shift = 0;
            B = circshift(B,shift);
            
            error = tune*imag(Wr_est).^2 + (B-A).^2;
        end
        
        if nargout > 1
            
            dfRotK_dRs = -IsRotK;
            
            dlamRot_sK_dRs = dfRotK_dRs./(1i*(kRot+we));
            
            dlamRot_rK_dRs = Lar/Lm*dlamRot_sK_dRs;
            dlamRot_rK_dLm = (Lal^2*IsRotK - Lal*LambdaRot_sK)/(Lm^2);
            dlamRot_rK_dLal = -2*Lal*IsRotK/Lm - 2*IsRotK + LambdaRot_sK/Lm;
            
            dIrRotK_dRs = Lm*dlamRot_sK_dRs/D - Las*dlamRot_rK_dRs/D;
            dIrRotK_dLm = LambdaRot_rK/(-D) + ((-LambdaRot_sK*(-D) + Lm*LambdaRot_sK*2*Lal) + Lm*(dlamRot_rK_dLm*-D - LambdaRot_rK*2*Lal) + (Lal*dlamRot_rK_dLm*(2*Lm*Lal + Lal^2) - Lal*LambdaRot_rK*2*Lal))/((-D)^2);
            dIrRotK_dLal = (LambdaRot_rK)/(-D) + ((Lm*LambdaRot_sK*(2*Lm + 2*Lal)) + ((Lm*dlamRot_rK_dLal)*(2*Lm*Lal + Lal^2) - Lm*LambdaRot_rK*(2*Lm + 2*Lal)) + Lal*(dlamRot_rK_dLal*(-D) - LambdaRot_rK*(2*Lm + 2*Lal)))/((-D)^2);
            
            dslipRot_dRs = (ifft(1i*kRot.*dlamRot_rK_dRs + Rr*dIrRotK_dRs).*ifft(1i*LambdaRot_rK) - ifft(1i*kRot.*LambdaRot_rK + Rr*IrRotK).*ifft(1i*dlamRot_rK_dRs))./(ifft(1i*LambdaRot_rK).^2);
            dslipRot_dLm = (ifft(1i*kRot.*dlamRot_rK_dLm + Rr*dIrRotK_dLm).*ifft(1i*LambdaRot_rK) - ifft(1i*kRot.*LambdaRot_rK + Rr*IrRotK).*ifft(1i*dlamRot_rK_dLm))./(ifft(1i*LambdaRot_rK).^2);
            dslipRot_dLal = (ifft(1i*kRot.*dlamRot_rK_dLal + Rr*dIrRotK_dLal).*ifft(1i*LambdaRot_rK) - ifft(1i*kRot.*LambdaRot_rK + Rr*IrRotK).*ifft(1i*dlamRot_rK_dLal))./(ifft(1i*LambdaRot_rK).^2);
            dslipRot_dRr = ifft(IrRotK)./ifft(1i*LambdaRot_rK);
            
            dwrRot_dRs = (2/(2*P))*dslipRot_dRs;
            dwrRot_dRr = (2/(2*P))*dslipRot_dRr;
            dwrRot_dLm = (2/(2*P))*dslipRot_dLm;
            dwrRot_dLal = (2/(2*P))*dslipRot_dLal;
            
            
            dfLoadK_dRs = -IsLoadK;
            
            dlamLoad_sK_dRs = dfLoadK_dRs./(1i*(kLoad+we));
            
            dlamLoad_rK_dRs = Lar/Lm*dlamLoad_sK_dRs;
            dlamLoad_rK_dLm = (Lal^2*IsLoadK - Lal*LambdaLoad_sK)/(Lm^2);
            dlamLoad_rK_dLal = -2*Lal*IsLoadK/Lm - 2*IsLoadK + LambdaLoad_sK/Lm;
            
            dIrLoadK_dRs = Lm*dlamLoad_sK_dRs/D - Las*dlamLoad_rK_dRs/D;
            dIrLoadK_dLm = LambdaLoad_rK/(-D) + ((-LambdaLoad_sK*(-D) + Lm*LambdaLoad_sK*2*Lal) + Lm*(dlamLoad_rK_dLm*-D - LambdaLoad_rK*2*Lal) + (Lal*dlamLoad_rK_dLm*(2*Lm*Lal + Lal^2) - Lal*LambdaLoad_rK*2*Lal))/((-D)^2);
            dIrLoadK_dLal = (LambdaLoad_rK)/(-D) + ((Lm*LambdaLoad_sK*(2*Lm + 2*Lal)) + ((Lm*dlamLoad_rK_dLal)*(2*Lm*Lal + Lal^2) - Lm*LambdaLoad_rK*(2*Lm + 2*Lal)) + Lal*(dlamLoad_rK_dLal*(-D) - LambdaLoad_rK*(2*Lm + 2*Lal)))/((-D)^2);
            
            dslipLoad_dRs = (ifft(1i*kLoad.*dlamLoad_rK_dRs + Rr*dIrLoadK_dRs).*ifft(1i*LambdaLoad_rK) - ifft(1i*kLoad.*LambdaLoad_rK + Rr*IrLoadK).*ifft(1i*dlamLoad_rK_dRs))./(ifft(1i*LambdaLoad_rK).^2);
            dslipLoad_dLm = (ifft(1i*kLoad.*dlamLoad_rK_dLm + Rr*dIrLoadK_dLm).*ifft(1i*LambdaLoad_rK) - ifft(1i*kLoad.*LambdaLoad_rK + Rr*IrLoadK).*ifft(1i*dlamLoad_rK_dLm))./(ifft(1i*LambdaLoad_rK).^2);
            dslipLoad_dLal = (ifft(1i*kLoad.*dlamLoad_rK_dLal + Rr*dIrLoadK_dLal).*ifft(1i*LambdaLoad_rK) - ifft(1i*kLoad.*LambdaLoad_rK + Rr*IrLoadK).*ifft(1i*dlamLoad_rK_dLal))./(ifft(1i*LambdaLoad_rK).^2);
            dslipLoad_dRr = ifft(IrLoadK)./ifft(1i*LambdaLoad_rK);
            
            dwrLoad_dRs = (2/(2*P))*dslipLoad_dRs;
            dwrLoad_dRr = (2/(2*P))*dslipLoad_dRr;
            dwrLoad_dLm = (2/(2*P))*dslipLoad_dLm;
            dwrLoad_dLal = (2/(2*P))*dslipLoad_dLal;
            
            Num_dRs_temp = ifft(1i*kRot.*dlamRot_rK_dRs + Rr*dIrRotK_dRs);
            Num_dRs = [repmat(Num_dRs_temp,s2,1); Num_dRs_temp(1:r2)] + ifft(1i*kLoad.*dlamLoad_rK_dRs + Rr*dIrLoadK_dRs);
            Num_dRs = Num_dRs - mean(Num_dRs)/2;
            
            Dem_dRs_temp = ifft(1i*dlamRot_rK_dRs);
            Dem_dRs = [repmat(Dem_dRs_temp,s2,1); Dem_dRs_temp(1:r2)] + ifft(1i*dlamLoad_rK_dRs);
            Dem_dRs = Dem_dRs - mean(Dem_dRs)/2;
            
            dslip_dRs = (Num_dRs.*combDem - combNum.*Dem_dRs)./(combDem.^2);
            
            Num_dLm_temp = ifft(1i*kRot.*dlamRot_rK_dLm + Rr*dIrRotK_dLm);
            Num_dLm = [repmat(Num_dLm_temp,s2,1); Num_dLm_temp(1:r2)] + ifft(1i*kLoad.*dlamLoad_rK_dLm + Rr*dIrLoadK_dLm);
            Num_dLm = Num_dLm - mean(Num_dLm)/2;
            
            Dem_dLm_temp = ifft(1i*dlamRot_rK_dLm);
            Dem_dLm = [repmat(Dem_dLm_temp,s2,1); Dem_dLm_temp(1:r2)] + ifft(1i*dlamLoad_rK_dLm);
            Dem_dLm = Dem_dLm - mean(Dem_dLm)/2;
            
            dslip_dLm = (Num_dLm.*combDem - combNum.*Dem_dLm)./(combDem.^2);
            
            Num_dLal_temp = ifft(1i*kRot.*dlamRot_rK_dLal + Rr*dIrRotK_dLal);
            Num_dLal = [repmat(Num_dLal_temp,s2,1); Num_dLal_temp(1:r2)] + ifft(1i*kLoad.*dlamLoad_rK_dLal + Rr*dIrLoadK_dLal);
            Num_dLal = Num_dLal - mean(Num_dLal)/2;
            
            Dem_dLal_temp = ifft(1i*dlamRot_rK_dLal);
            Dem_dLal = [repmat(Dem_dLal_temp,s2,1); Dem_dLal_temp(1:r2)] + ifft(1i*dlamLoad_rK_dLal);
            Dem_dLal = Dem_dLal - mean(Dem_dLal)/2;
            
            dslip_dLal = (Num_dLal.*combDem - combNum.*Dem_dLal)./(combDem.^2);
            
            Num_dRr_temp = ifft(IrRotK);
            Num_dRr = [repmat(Num_dRr_temp,s2,1); Num_dRr_temp(1:r2)] + ifft(IrLoadK);
            Num_dRr = Num_dRr - mean(Num_dRr)/2;
            
            dslip_dRr = Num_dRr./combDem;
            
            dwr_dRs = (2/(2*P))*dslip_dRs;
            dwr_dRr = (2/(2*P))*dslip_dRr;
            dwr_dLm = (2/(2*P))*dslip_dLm;
            dwr_dLal = (2/(2*P))*dslip_dLal;
            
            
            if ~useTrueSpeed
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(wr_mean - Wr_Rot_mean).*(-mean(real(dwrRot_dRs))) + 2*(wr_mean - Wr_Load_mean).*(-mean(real(dwrLoad_dRs))))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(wr_mean - Wr_Rot_mean).*(-mean(real(dwrRot_dRr))) + 2*(wr_mean - Wr_Load_mean).*(-mean(real(dwrLoad_dRr))))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(wr_mean - Wr_Rot_mean).*(-mean(real(dwrRot_dLm))) + 2*(wr_mean - Wr_Load_mean).*(-mean(real(dwrLoad_dLm))))*Scaling_Factors(3);
                derror_dLal = (2*tune*imag(Wr_est).*imag(dwr_dLal) + 2*(wr_mean - Wr_Rot_mean).*(-mean(real(dwrRot_dLal))) + 2*(wr_mean - Wr_Load_mean).*(-mean(real(dwrLoad_dLal))))*Scaling_Factors(4);
            else
                N = length(Wr_est);
                B = Speed_long(1e4:1e4+N-1);
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(B - real(Wr_est)).*(-real(dwr_dRs)))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(B - real(Wr_est)).*(-real(dwr_dRr)))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(B - real(Wr_est)).*(-real(dwr_dLm)))*Scaling_Factors(3);
                derror_dLal = (2*tune*imag(Wr_est).*imag(dwr_dLal) + 2*(B - real(Wr_est)).*(-real(dwr_dLal)))*Scaling_Factors(4);
            end
            
            
            J = [derror_dRs derror_dRr derror_dLm derror_dLal];
        end
        
    end

    function dx = induction_motor_parameter_fitting_model_arbitrary_load(t, x)
        %  HI BUNNY I AM MOST AWESUM CODE IN DA UNIVERSE
        
        lamqs = x(1);
        lamds = x(2);
        lamqr = x(3);
        lamdr = x(4);
        
        VDS = real(Vdq_rot(1)); %%% these values are constant
        VQS = imag(Vdq_rot(1));
        
        if t > RotTime(end)
            time = mod(t,RotTime(end));
        else
            time = t;
        end
        
        Wr = interp1(RotTime,real(out.Wr_est),time);
        
        idr   = (Lm*lamds - Las*lamdr)/D;
        iqr   = (Lm*lamqs - Las*lamqr)/D;
        iqs   = (Lm*lamqr - Lar*lamqs)/D;
        ids   = (Lm*lamdr - Lar*lamds)/D;
        
        s1 = (VQS - we*lamds - Rs*iqs);
        s2 = (VDS + we*lamqs - Rs*ids);
        s3 = (P*Wr - we)*lamdr - Rr*iqr;
        s4 = (we - P*Wr)*lamqr - Rr*idr;
        
        dx = [s1 s2 s3 s4 we]';
    end


    function [error, J, Wr_est, sim_state] = error_function_mean_speed_with_Jacobian_Rot_5_param(Parameters)
        Rs = Parameters(1)*Scaling_Factors(1);
        Rr = Parameters(2)*Scaling_Factors(2);
        Lm = Parameters(3)*Scaling_Factors(3);
        Ls = Parameters(4)*Scaling_Factors(4);
        Lr = Parameters(5)*Scaling_Factors(5);
        
        Las = Lm + Ls;
        Lar = Lm + Lr;
        D = Lm^2 - Las*Lar;
        
        [Num, Dem, VsK, IsK, fK, Lambda_sK, Lambda_rK, IrK, k ...
            Num_m, Dem_m, IsK_m, fK_m, Lambda_sK_m, Lambda_rK_m, IrK_m] = partial_fft_invert(Idq_rot, Vdq_rot);
        
        slip = Num./Dem;
        slip_m = Num_m./Dem_m;
        
        Wr_est = (2/(2*P))*(slip + we);
        Wr_est_m = (2/(2*P))*(slip_m + we);
        wr_mean = out.wr;
        
        
        if ~useTrueSpeed
            error = tune*imag(Wr_est).^2 + (wr_mean-mean(real(Wr_est_m))).^2;
        else
            N = length(Wr_est);
            A = real(Wr_est)/2/pi;
            B = Speed_long(1e4:1e4+N-1)/2/pi;
            xc = ifft(fft(A).*conj(fft(B)));
            [junk, shift] = max(xc);
            shift = 0;
            B = circshift(B,shift);
            error = tune*imag(Wr_est).^2 + (B-A).^2;
        end
        
        if nargout > 1
            
            dfK_dRs = -IsK;
            
            dlam_sK_dRs = dfK_dRs./(1i*(k+we));
            
            dlam_rK_dRs = Lar/Lm*dlam_sK_dRs;                   % converted
            dlam_rK_dLm = (Lr*Ls*IsK - Lr*Lambda_sK)/(Lm^2);    % converted
            dlam_rK_dLs = -IsK - Lr*IsK/Lm;                     % converted
            dlam_rK_dLr = -IsK - Ls*IsK/Lm + Lambda_sK/Lm;      % converted
            
            dIrK_dRs = Lm*dlam_sK_dRs/D - Las*dlam_rK_dRs/D;    % converted
            dIrK_dLm = Lambda_rK/(-D) + ((-Lambda_sK*(-D) + Lm*Lambda_sK*(Ls+Lr)) + Lm*(dlam_rK_dLm*-D - Lambda_rK*(Lr+Ls)) + (Ls*dlam_rK_dLm*(-D) - Ls*Lambda_rK*(Ls+Lr)))/((-D)^2);       % converted
            dIrK_dLs = (Lambda_rK)/(-D) + ((Lm*Lambda_sK*(Lm + Lr)) + ((Lm*dlam_rK_dLs)*(-D) - Lm*Lambda_rK*(Lm + Lr)) + Ls*(dlam_rK_dLs*(-D) - Lambda_rK*(Lm + Lr)))/((-D)^2);             % converted
            dIrK_dLr = ((Lm*Lambda_sK*(Lm + Ls)) + ((Lm*dlam_rK_dLr)*(-D) - Lm*Lambda_rK*(Lm + Ls)) + Ls*(dlam_rK_dLr*(-D) - Lambda_rK*(Lm + Ls)))/((-D)^2);                                % converted
            
            dslip_dRs = (ifft(1i*k.*dlam_rK_dRs + Rr*dIrK_dRs).*ifft(1i*Lambda_rK) - ifft(1i*k.*Lambda_rK + Rr*IrK).*ifft(1i*dlam_rK_dRs))./(ifft(1i*Lambda_rK).^2);    % converted
            dslip_dLm = (ifft(1i*k.*dlam_rK_dLm + Rr*dIrK_dLm).*ifft(1i*Lambda_rK) - ifft(1i*k.*Lambda_rK + Rr*IrK).*ifft(1i*dlam_rK_dLm))./(ifft(1i*Lambda_rK).^2);    % converted
            dslip_dLs = (ifft(1i*k.*dlam_rK_dLs + Rr*dIrK_dLs).*ifft(1i*Lambda_rK) - ifft(1i*k.*Lambda_rK + Rr*IrK).*ifft(1i*dlam_rK_dLs))./(ifft(1i*Lambda_rK).^2);    % converted
            dslip_dLr = (ifft(1i*k.*dlam_rK_dLr + Rr*dIrK_dLr).*ifft(1i*Lambda_rK) - ifft(1i*k.*Lambda_rK + Rr*IrK).*ifft(1i*dlam_rK_dLr))./(ifft(1i*Lambda_rK).^2);    % converted
            dslip_dRr = ifft(IrK)./ifft(1i*Lambda_rK);          % converted
            
            dwr_dRs = (2/(2*P))*dslip_dRs;                      % converted
            dwr_dRr = (2/(2*P))*dslip_dRr;                      % converted
            dwr_dLm = (2/(2*P))*dslip_dLm;                      % converted
            dwr_dLs = (2/(2*P))*dslip_dLs;                      % converted
            dwr_dLr = (2/(2*P))*dslip_dLr;                      % converted
            %%%%%%%%%%%%%%%5
            if ~useTrueSpeed  % don't calculate mean based partials if using true speed.
                
                dfK_dRs_m = -IsK_m;
                
                dlam_sK_dRs_m = dfK_dRs_m./(1i*(k+we));
                
                dlam_rK_dRs_m = Lar/Lm*dlam_sK_dRs_m;
                dlam_rK_dLm_m = (Lr*Ls*IsK_m - Lr*Lambda_sK_m)/(Lm^2);
                dlam_rK_dLs_m = -IsK_m - Lr*IsK_m/Lm;
                dlam_rK_dLr_m = -IsK_m - Ls*IsK_m/Lm + Lambda_sK_m/Lm;
                
                dIrK_dRs_m = Lm*dlam_sK_dRs_m/D - Las*dlam_rK_dRs_m/D;    % converted
                dIrK_dLm_m = Lambda_rK_m/(-D) + ((-Lambda_sK_m*(-D) + Lm*Lambda_sK_m*(Ls+Lr)) + Lm*(dlam_rK_dLm_m*-D - Lambda_rK_m*(Lr+Ls)) + (Ls*dlam_rK_dLm_m*(-D) - Ls*Lambda_rK_m*(Ls+Lr)))/((-D)^2);       % converted
                dIrK_dLs_m = (Lambda_rK_m)/(-D) + ((Lm*Lambda_sK_m*(Lm + Lr)) + ((Lm*dlam_rK_dLs_m)*(-D) - Lm*Lambda_rK_m*(Lm + Lr)) + Ls*(dlam_rK_dLs_m*(-D) - Lambda_rK_m*(Lm + Lr)))/((-D)^2);             % converted
                dIrK_dLr_m = ((Lm*Lambda_sK_m*(Lm + Ls)) + ((Lm*dlam_rK_dLr_m)*(-D) - Lm*Lambda_rK_m*(Lm + Ls)) + Ls*(dlam_rK_dLr_m*(-D) - Lambda_rK_m*(Lm + Ls)))/((-D)^2);                                % converted
                
                dslip_dRs_m = (ifft(1i*k.*dlam_rK_dRs_m + Rr*dIrK_dRs_m).*ifft(1i*Lambda_rK_m) - ifft(1i*k.*Lambda_rK_m + Rr*IrK_m).*ifft(1i*dlam_rK_dRs_m))./(ifft(1i*Lambda_rK_m).^2);
                dslip_dLm_m = (ifft(1i*k.*dlam_rK_dLm_m + Rr*dIrK_dLm_m).*ifft(1i*Lambda_rK_m) - ifft(1i*k.*Lambda_rK_m + Rr*IrK_m).*ifft(1i*dlam_rK_dLm_m))./(ifft(1i*Lambda_rK_m).^2);
                dslip_dLs_m = (ifft(1i*k.*dlam_rK_dLs_m + Rr*dIrK_dLs_m).*ifft(1i*Lambda_rK_m) - ifft(1i*k.*Lambda_rK_m + Rr*IrK_m).*ifft(1i*dlam_rK_dLs_m))./(ifft(1i*Lambda_rK_m).^2);
                dslip_dLr_m = (ifft(1i*k.*dlam_rK_dLr_m + Rr*dIrK_dLr_m).*ifft(1i*Lambda_rK_m) - ifft(1i*k.*Lambda_rK_m + Rr*IrK_m).*ifft(1i*dlam_rK_dLr_m))./(ifft(1i*Lambda_rK_m).^2);
                dslip_dRr_m = ifft(IrK_m)./ifft(1i*Lambda_rK_m);
                
                dwr_dRs_m = (2/(2*P))*dslip_dRs_m;
                dwr_dRr_m = (2/(2*P))*dslip_dRr_m;
                dwr_dLm_m = (2/(2*P))*dslip_dLm_m;
                dwr_dLs_m = (2/(2*P))*dslip_dLs_m;
                dwr_dLr_m = (2/(2*P))*dslip_dLr_m;
            end
            %%%%%%%%%%%%%%%%%%5
            if ~useTrueSpeed
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dRs_m))))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dRr_m))))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dLm_m))))*Scaling_Factors(3);
                derror_dLs = (2*tune*imag(Wr_est).*imag(dwr_dLs) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dLs_m))))*Scaling_Factors(4);
                derror_dLr = (2*tune*imag(Wr_est).*imag(dwr_dLr) + 2*(wr_mean - mean(real(Wr_est_m))).*(-mean(real(dwr_dLr_m))))*Scaling_Factors(5);
            else
                N = length(Wr_est);
                B = Speed_long(1e4:1e4+N-1);
                derror_dRs = (2*tune*imag(Wr_est).*imag(dwr_dRs) + 2*(B - real(Wr_est)).*(-real(dwr_dRs)))*Scaling_Factors(1);
                derror_dRr = (2*tune*imag(Wr_est).*imag(dwr_dRr) + 2*(B - real(Wr_est)).*(-real(dwr_dRr)))*Scaling_Factors(2);
                derror_dLm = (2*tune*imag(Wr_est).*imag(dwr_dLm) + 2*(B - real(Wr_est)).*(-real(dwr_dLm)))*Scaling_Factors(3);
                derror_dLs = (2*tune*imag(Wr_est).*imag(dwr_dLs) + 2*(B - real(Wr_est)).*(-real(dwr_dLs)))*Scaling_Factors(4);
                derror_dLr = (2*tune*imag(Wr_est).*imag(dwr_dLr) + 2*(B - real(Wr_est)).*(-real(dwr_dLr)))*Scaling_Factors(5);
            end
            
            J = [derror_dRs derror_dRr derror_dLm derror_dLs derror_dLr];
        end
        
        if nargout > 3
            
            Lambda_s = ifft(Lambda_sK);
            Lambda_r = ifft(Lambda_rK);
            
            sim_state = [imag(Lambda_s(1)), real(Lambda_s(1)) imag(Lambda_r(1)), real(Lambda_r(1)), 0];
        end
        
        
    end



end
