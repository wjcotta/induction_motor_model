function [data] = scan_load_wr(pulley_ratio,load_wr_range,d_wr,data,tune,useTrueSpeed)
%%% This function scans through a series user defined load speeds running
%%% Chris Schantz's inversion code to search for a best fit.


if nargin <6
    useTrueSpeed=false;
end

P=data.P;
load_wr = data.wr*pulley_ratio;
wr_candidates = [load_wr_range(1) : d_wr : load_wr]+...
                [load_wr+d_wr : d_wr : load_wr_range(2)]


%data = decompose_indmotor_data(data, load_wr);
%data = speed_inversion_error_space_sim_fixed( data, load_wr, P, tune, useTrueSpeed);
end