%readout for XFCT measurements

clf, clear, close all
filepath = '~/Dropbox/CT_MATHEMATICS/MATLAB_CT/XFCT_DATA/';

angle = '90';
tube_amps = '2';
tube_voltage_vec = 80:100;

delta_E = 1;

alpha_noise = zeros(1,length(tube_voltage_vec));
beta_noise = zeros(1,length(tube_voltage_vec));

i = 1;
for tube_voltage = tube_voltage_vec
    
    load(['XFCT_data' angle '_' tube_amps '_' num2str(tube_voltage)];
    
    beta_index = find(Es > 78,1,'first');
    
    alpha_index = find(Es > 68.8,1,'first');
    
    alpha_noise(i) = sum(data(alpha_index-delta_E:alpha_index+delta_E));
    beta_noise(i) = sum(data(beta_index-delta_E:beta_index+delta_E));
    i = i +1;
end

save('XFCT_noise')

%%
clf, clear, close all
load('XFCT_noise')
plot(tube_voltage_vec,alpha_noise)

plot(tube_voltage_vec,beta_noise)

