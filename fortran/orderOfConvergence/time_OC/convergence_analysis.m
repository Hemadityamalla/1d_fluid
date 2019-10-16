%Order of convergence checking- too lazy to be able to read files in python

clear; clc;format long
skip_line = 1; % The first line has the field names

%test for time convergence analysis
%Exact solution dx = 0.25, dt = 0.001
%exact = importdata('exact.txt',' ',skip_line); 
dt = importdata('dt_first.txt',' ',skip_line); %dt = 0.1
dthalf = importdata('dt_2_first.txt',' ',skip_line); %dt = 0.05
dtfourth = importdata('dt_4_first.txt',' ',skip_line); %dt = 0.025
dteighth = importdata('dt_8_first.txt',' ',skip_line); %dt = 0.0125
dtsixteen = importdata('dt_16_first.txt',' ', skip_line); %dt = 0.00625
dt32 = importdata('dt_32_first.txt',' ', skip_line); %dt = 0.003125



%For the electron density (sigma)
%order1 = log2(norm(dt.data(:,2) - exact.data(:,2))/norm(dthalf.data(:,2) - exact.data(:,2)))
field = 2; %2,3,4

order2 = log2(norm(dt.data(:,field) - dthalf.data(:,field))/norm(dthalf.data(:,field) - dtfourth.data(:,field)))

order3 = log2(norm(dthalf.data(:,field) - dtfourth.data(:,field))/norm(dtfourth.data(:,field) - dteighth.data(:,field)))

order4 = log2(norm(dtfourth.data(:,field) - dteighth.data(:,field))/norm(dteighth.data(:,field) - dtsixteen.data(:,field)))

order5 = log2(norm(dteighth.data(:,field) - dtsixteen.data(:,field))/norm(dt32.data(:,field) - dtsixteen.data(:,field)))