pulley_ratio = .54154;
sim_flag     = false;
scan_flag    = true;
%{
filename = 'oct23_noload';
data_no_load = process_data(filename, pulley_ratio, sim_flag)
pause
pulley_ratio = 2*pulley_ratio;
filename = 'oct23_alum_damaged';
data_alum_damaged = process_data(filename, pulley_ratio, sim_flag)
pause
filename = 'oct23_alum_intact';
data_alum_intact = process_data(filename, pulley_ratio, sim_flag)
pause

filename = 'oct23_wood_intact';
data_wood_intact = process_data(filename, pulley_ratio, sim_flag)
pause

filename = 'oct23_wood_damaged';
data_wood_damaged = process_data(filename, pulley_ratio, sim_flag)
pause

%}

filename = 'oct23_alum_damaged';
data = process_data(filename, pulley_ratio, sim_flag, scan_flag);