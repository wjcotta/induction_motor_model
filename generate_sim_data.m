% This file generates sim data (currents, voltages, and speed) using the
% induction motor model. 

% Run motor simulation using parameters defined in 'indparam'.
% 'runind' provides the time, 'tout', and the state, 'yout'. 
runind;

% Convert from dq0 frame back to abc frame.
% 'data' contains 6 columns: Va, Vb, Vc, Ia, Ib, Ic
[idq, vabc, iabc, ~, ~] = conv2abc(t, statev);
data = [vabc iabc];
close all;

% Append data with the motor speed from the solution to runind. 
Speed = statev(:,5)/P;
data = [data, Speed];

plot(t, data(:,7))
xlabel('Time (s)')
ylabel('Mechanical Rotor Speed (rad/s)')
title('Motor Spin-up')
axis([0 45 0 125])

% We now have 7 columns of data, saved in 'sim_data'.
save('sim_data', 'data')

% Save dq0 currents for debugging.
save('sim_idq', 'idq')

