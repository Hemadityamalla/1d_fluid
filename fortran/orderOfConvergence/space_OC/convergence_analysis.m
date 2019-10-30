%Order of convergence checking- too lazy to be able to read files in python

clear; clc;format long
skip_line = 1; % The first line has the field names
field = 2; %2,3,4

%test for space convergence analysis- had to take second order time to
%avoid really small time steps


dx = importdata('0.5.txt',' ',skip_line);
dx_2 = importdata('0.25.txt',' ',skip_line);
dx_4 = importdata('0.125.txt',' ',skip_line);
dx_8 = importdata('0.0625.txt',' ',skip_line);
dx_16 = importdata('0.03125.txt',' ',skip_line);

order1 = log2(norm(dx.data(:,field) - dx_2.data(1:2:end,field),2)/norm(dx_2.data(:,field) - dx_4.data(1:2:end,field),2))

order2 = log2(norm(dx_2.data(:,field) - dx_4.data(1:2:end,field),2)/norm(dx_4.data(:,field) - dx_8.data(1:2:end,field),2))


order3 = log2((norm(dx_4.data(:,field) - dx_8.data(1:2:end,field),2)/norm(dx_8.data(:,field) - dx_16.data(1:2:end,field),2)))