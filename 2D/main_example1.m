%main_example1.m 

%Example 1 Section 3.1 of Cortez, SIAM J. Sci Comput. 2001 
%Flow past a cylinder of radius a 
%Setting velocity on cylinder to be (1,0) and then computing the velocity
%in the fluid surrounding the cylinder 

%Developed by Shilpa Khatri and Ricardo Cortez  
%July 2024 

clear all
close all 

%% Parameters to set 

%setting the viscosity 
mu = 1; 

%radius of cylinder  
a = 1; 

%domain on which velocity is computed 
x1min = -2; 
x1max = 2; 
x2min = -2; 
x2max = 2; 

%number of points on boundary where velocity is set 
N = 200;  

%resolution for velocity grid in fluid 
Nx1 = 80; 
Nx2 = 80; 

%% Setting velocity on cylinder boundary, computing force on cylinder boundary, 
% and then computing velocity in fluid domain

%discretization of cylinder boundary 
dt = 2*pi/N; 
t = 0:dt:2*pi-dt/2;
t = t';

%regularization parameter
ep = dt/4; 

%cylinder boundary on which forces will be computed from setting velocity 
y1 = a*cos(t); 
y2 = a*sin(t);

%velocity on cylinder boundary 
v1 = 1 + 0*t; 
v2 = 0*t; 

%computing forces on cylinder boundary 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[v1,v2],ep,mu);
f1 = f(:,1); 
f2 = f(:,2);

%points on which velocity will be computed from forces on circle 
xx1 = linspace(x1min,x1max,Nx1);
xx2 = linspace(x2min,x2max,Nx2); 
[x1m,x2m] = ndgrid(xx1,xx2); 
x1 = x1m(:);
x2 = x2m(:);

%computing velocity on grid 
u = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu);
u1 = u(:,1); 
u2 = u(:,2);
u1m = reshape(u1,size(xx1,2),size(xx2,2)); 
u2m = reshape(u2,size(xx1,2),size(xx2,2));

%% Computing error 

%exact solution 
fo1 = 8*pi/(1-2*log(a));

for i = 1:length(xx1)

    for j = 1:length(xx2)

        r = sqrt(x1m(i,j).^2 + x2m(i,j).^2); %radius 
        fodotx = fo1*x1m(i,j); %fo dot x

        %computing error 
        if (r >= 1)
            
            uexact1(i,j) = -fo1*(2*log(r) - a^2/(r^2))/8/pi + fodotx*x1m(i,j)*(1-a^2/(r^2))/4/pi/(r^2); 
            uexact2(i,j) = fodotx*x2m(i,j)*(1-a^2/(r^2))/4/pi/(r^2);

            error1(i,j) = abs(u1m(i,j)-uexact1(i,j));
            error2(i,j) = abs(u2m(i,j)-uexact2(i,j)); 

        else
        %setting error interior to cylinder = 0
            
            error1(i,j) = 0; 
            error2(i,j) = 0; 

        end

    end
end

errormag = sqrt(error1.^2 + error2.^2);

%prints the max error 
fprintf('maximum error in u1: %d \n',max(max(error1)));
fprintf('maximum error in u2: %d \n',max(max(error2)));

%% Plotting figures 
skip = 4; %for quiver plots 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

figure(1) 
plot(y1,y2,'k.-')
hold on 
quiver(y1,y2,v1,v2,'r')
axis equal 
quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end)-sum(f1)/8/pi/mu,u2m(1:skip:end,1:skip:end),'k')
xlim([x1min,x1max])
ylim([x2min,x2max])
title('Velocity set and Computed Fluid Velocity')

figure(2)
plot(y1,y2,'k.-')
hold on 
quiver(y1,y2,v1,v2,'r')
axis equal 
quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),uexact1(1:skip:end,1:skip:end),uexact2(1:skip:end,1:skip:end),'k')
xlim([x1min,x1max])
ylim([x2min,x2max])
title('Velocity set and Exact Fluid Velocity')

figure(3) 
plot(y1,y2,'k.-')
hold on 
quiver(y1,y2,f1,f2,'r')
axis equal 
xlim([x1min,x1max])
ylim([x2min,x2max])
title('Forces Computed on Cylinder')

figure(4)  
pcolor(x1m,x2m,log10(errormag)+eps)
shading interp
colorbar
caxis( [-5 -2] )
hold on
plot(y1,y2,'k.-') 
axis equal 
xlim([x1min,x1max])
ylim([x2min,x2max])
title('Error Magnitude in the Fluid Velocity')


