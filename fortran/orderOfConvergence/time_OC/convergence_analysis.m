%Order of convergence checking- too lazy to be able to read files in python

clear; clc;format long
skip_line = 1; % The first line has the field names
field = 2; %2,3,4

%test for time convergence analysis
% dt = importdata('dt_first.txt',' ',skip_line); %dt = 0.1
% dthalf = importdata('dt_2_first.txt',' ',skip_line); %dt = 0.05
% dtfourth = importdata('dt_4_first.txt',' ',skip_line); %dt = 0.025
% dteighth = importdata('dt_8_first.txt',' ',skip_line); %dt = 0.0125
% dtsixteen = importdata('dt_16_first.txt',' ', skip_line); %dt = 0.00625
% dt32 = importdata('dt_32_first.txt',' ', skip_line); %dt = 0.003125
% 
% 
% 
% order2 = log2(norm(dt.data(:,field) - dthalf.data(:,field))/norm(dthalf.data(:,field) - dtfourth.data(:,field)))
% 
% order3 = log2(norm(dthalf.data(:,field) - dtfourth.data(:,field))/norm(dtfourth.data(:,field) - dteighth.data(:,field)))
% 
% order4 = log2(norm(dtfourth.data(:,field) - dteighth.data(:,field))/norm(dteighth.data(:,field) - dtsixteen.data(:,field)))
% 
% order5 = log2(norm(dteighth.data(:,field) - dtsixteen.data(:,field))/norm(dt32.data(:,field) - dtsixteen.data(:,field)))

%test for space convergence analysis- had to take second order time to
%avoid really small time steps


dx = importdata('0.25',' ',skip_line);
dx_2 = importdata('0.125',' ',skip_line);
dx_4 = importdata('0.0625',' ',skip_line);
dx_8 = importdata('0.03125',' ',skip_line);

order1 = log2(norm(dx.data(:,field) - dx_2.data(:,field))/norm(dx_2.data(:,field) - dx_4.data(:,field)))

order2 = log2(norm(dx_2.data(:,field) - dx_4.data(:,field))/norm(dx_4.data(:,field) - dx_8.data(:,field)))