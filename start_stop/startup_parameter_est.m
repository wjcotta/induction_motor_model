load('startup8')


encoder_count = 1800;
mill= 'heavy_duty'
P=2;
l2l_flag = 'False';
fs = 8000;

data8 = Pre_process_three_phase_motor_data('startup8', fs, l2l_flag, P, encoder_count, mill)
data1 = Pre_process_three_phase_motor_data('startup1', fs, l2l_flag, P, encoder_count, mill,1)

