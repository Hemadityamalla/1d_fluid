%Order of convergence checking- too lazy to be able to read files in python

clear; clc;format long
skip_line = 1; % The first line has the field names

%test for time convergence analysis
%Exact solution dx = 0.25, dt = 0.001
exact = importdata('exact.txt',' ',skip_line); 
dt = importdata('dt.txt',' ',skip_line); %dt = 0.1
dthalf = importdata('dt_2.txt',' ',skip_line); %dt = 0.05
dtfourth = importdata('dt_4.txt',' ',skip_line); %dt = 0.025
dteighth = importdata('dt_eighth.txt',' ',skip_line); %dt = 0.0125


%For the electron density (sigma)
order1 = log2(norm(dt.data(:,2) - exact.data(:,2))/norm(dthalf.data(:,2) - exact.data(:,2)))

order2 = log2(norm(dt.data(:,2) - dthalf.data(:,2))/norm(dthalf.data(:,2) - dtfourth.data(:,2)))

order3 = log2(norm(dthalf.data(:,2) - dtfourth.data(:,2))/norm(dtfourth.data(:,2) - dteighth.data(:,2)))