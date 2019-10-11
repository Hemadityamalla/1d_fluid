%Order of convergence checking- too lazy to be able to read files in python

clear; clc;format long
skip_line = 1; % The first line has the field names
h_data = importdata('200_dx.txt', ' ', skip_line);
h2_data = importdata('200_dx_half.txt', ' ', skip_line);
h4_data = importdata('200_dx_quarter.txt', ' ', skip_line);
h8_data = importdata('200_dx_eighth.txt', ' ', skip_line);


%For the electron density (sigma)
order1 = log2(norm(abs(h_data.data(:,2) - h2_data.data(1:2:end,2)))/norm(abs(h2_data.data(:,2) - h4_data.data(1:2:end,2))))

order2 = log2(norm(abs(h2_data.data(:,2) - h4_data.data(1:2:end,2)))/norm(abs(h4_data.data(:,2) - h8_data.data(1:2:end,2))))