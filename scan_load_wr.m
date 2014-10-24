scan_load_wr(pulley_ratio,load_wr_range,d_wr,data,tune,useTrueSpeed)
if nargin <6
    useTrueSpeed=false;
end
P=data.P;
load_wr = data.wr*pulley_ratio;
wr_candidates = [load_wr_range(1):dt:load_wr]+[load_wr+dt:dt:load_wr_range(2)]
data = speed_inversion_error_space_sim_fixed( data, load_wr, P, tune, useTrueSpeed);