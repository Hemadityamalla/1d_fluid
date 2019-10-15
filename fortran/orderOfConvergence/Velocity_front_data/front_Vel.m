%Script to calculate and compare the front propagation velocity


clear; clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

xpos = 500;ypos = 500; width = 1200; height = 800;


skip_line = 1; % The first line has the field names
h_data = importdata('frontPos.dat', ' ', skip_line);
h_data_small = importdata('frontPos_small.dat', ' ', skip_line);


%Asymptoic front velocity
E_b = 1.0;
D = 0.1;
v_exact = E_b + 2.0*sqrt(D*E_b*exp(-1.0/E_b));

%Theoretical (todo)
%zb = 31.0;
%ztheoretical = @(tau) zb + E_b*tau + sqrt*(E_b + sqrt(4.0*D*E_b*exp(-1.0/E_b))); 

%Simulation
dt = 0.01;
times = h_data.data(:,1); 
[fpos,ia,ic] = unique(h_data.data(:,2));
tpos = times(ia);
vfront = diff(fpos)./diff(tpos);

times_small = h_data_small.data(:,1); 
[fpos_small,ia,ic] = unique(h_data_small.data(:,2));
tpos_small = times(ia);
vfront_small = diff(fpos_small)./diff(tpos_small);
figure(1)
plot(tpos(2:5:end),vfront(1:5:end),'ro-',tpos_small(2:5:end),vfront_small(1:5:end),'go-', times, v_exact*ones(length(times),1),'b--');
xlabel('\tau');ylabel('v_f');
legend('Numerical','Numerical small dx','asymptotic');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;