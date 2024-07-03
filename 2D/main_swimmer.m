%main_swimmer.m 

%swimmmer on which the velocity is prescribed and then computing the
%velocity in the fluid surrounding the swimmer assuming net force and net
%torque is zero. 

%Developed by Shilpa Khatri and Ricardo Cortez  
%July 2024 

clear all
close all


%% Parameters to set 

%viscosity 
mu = 1; 

%wavenumber of flagella; 
k = 2*pi;

%frequency 
w = 2*pi; 

%final parameterization of flagella
slength = 1; 

%number of points on flagella 
Ns = 40; 

%radius of body
a = 0.25;

%domain on which velocity will be computed 
x1min = -1; 
x1max = 2; 
x2min = -1; 
x2max = 1; 

%resolution for velocity grid in fluid 
Nx1 = 60; 
Nx2 = 40; 

%final time 
tfinal = 1; 

%number of timesteps 
tsteps = 20; 

%% Setting velocity of swimmer, computing force on swimmer 
% and then computing velocity in fluid domain

%discretization of flagella - discretization is set same for body 
ds = slength/Ns; 
s = 0:ds:slength;  %parameterization
s = s';

%discretization is set the same for the body 
bodycirc = 2*pi*a; 
N = floor(bodycirc/ds); 
dt = 2*pi/N; 
t = dt:dt:2*pi-dt/2; %do not want a point at 0 since on flagella 
t = t'; 

%location of body 
zc1 = -a; 
z1 = zc1 + a*cos(t);
z2 = a*sin(t); 

%velocity of body 
vz1 = 0*z1; 
vz2 = 0*z2; 

%regularization parameter
ep = ds/4; 

%points on which velocity will be computed from forces
xx1 = linspace(x1min,x1max,Nx1);
xx2 = linspace(x2min,x2max,Nx2); 
[x1m,x2m] = meshgrid(xx1,xx2); 
x1 = x1m(:);
x2 = x2m(:);

%timestep 
dtstep = tfinal/tsteps; 
counter = 1; 

%for plots 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

%looping over time
for t = 0:dtstep:tfinal
    %if want to run to the end of one period run to tfinal - dt 

    %location of flagella 
    y1 = s; 
    y2 = 0.2*s.*sin(k*s-w*t);

    %velocity of flagella (derivative wrt time of location) 
    vy1 = 0*s; 
    vy2 = -0.2*w*s.*cos(k*s-w*t); 

    %concatenate location of swimmer and flagella 
    sw1 = [y1;z1];
    sw2 = [y2;z2]; 

    %contcateante velocity of swimmer and flagella 
    v1 = [vy1; vz1];
    v2 = [vy2; vz2];

    %computing forces on swimmer 
    [f,uo,omega] = RegStokeslets2D_velocitytoforce_augmented([sw1 sw2],[sw1 sw2],[v1 v2],ep,mu);

    %%computing velocity on swimmer from forces - just to check 
    uswimmer = RegStokeslets2D_forcetovelocity([sw1 sw2],f,[sw1 sw2],ep,mu);
    %%adding in background velocities to swimmer velocity - checks that we are getting what we expect 
    %uswimmer_b = [uswimmer(:,1)-uo(1)+sw2*omega uswimmer(:,2)-uo(2)-sw1*omega];
    %max(max([abs(uswimmer_b(:,1)-v1) abs(uswimmer_b(:,2)-v2)]))

    %computing velocity on fluid grid 
    u = RegStokeslets2D_forcetovelocity([sw1 sw2],f,[x1 x2],ep,mu);
    u1 = u(:,1);
    u2 = u(:,2); 
    %changing frame of reference to swimmer 
    u1_b = u(:,1) - uo(1) + x2*omega;
    u2_b = u(:,2) - uo(2) - x1*omega;
    u1m = reshape(u1,size(xx2,2),size(xx1,2)); 
    u2m = reshape(u2,size(xx2,2),size(xx1,2));
    %changing frame of reference to swimmer  
    u1m_b = u1m - uo(1) + x2m*omega;
    u2m_b = u2m - uo(2) - x1m*omega; 

    %keeping track of uo and omega in time 
    uot(counter,:) = uo'; 
    omegat(counter) = omega; 

    %update counter 
    counter = counter + 1; 

    figure(1)
    hold off
    quiver([sw1;x1],[sw2;x2],[v1;u1_b],[v2;u2_b],'k')
    hold on
    plot(y1,y2,'.-k')
    plot(z1,z2,'.-k')
    axis equal
  
    figure(2)
    hold off 
    pcolor(x1m,x2m,sqrt(u1m_b.^2+u2m_b.^2) )
    shading interp
    colorbar
    clim([0,1])
    hold on
    h = streamslice(x1m,x2m,u1m_b,u2m_b);
    set( h, 'Color', [0.7 0.7 0.7] )
    set( h, 'LineWidth', 1)
    plot(y1,y2,'.-k')
    plot(z1,z2,'.-k')
    axis equal 
    
    pause(0.2)

end

