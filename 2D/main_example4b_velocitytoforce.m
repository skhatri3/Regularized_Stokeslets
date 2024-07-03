%main_example4b_velocitytoforce

%Example 4b Section 4.2 of Cortez, SIAM J. Sci Comput. 2001 
%Setting velocity on boundary of a cylinder of radius 1 and computing the
%forces at those same points 

%Note: you have to be careful as this solution is nonunique - for any
%solution on a circle you can add constant magnitude normal forces which
%does not add any velocity (this means as we increase N, the problem should
%be more ill posed). 

%Developed by Shilpa Khatri and Ricardo Cortez  
%July 2024 

clear all 
close all

%% Parameters to set 

%setting the viscosity
mu = 1; 

%number of points on boundary where velocity is set and force is computed 
N = 100;    

%% Setting forces and computing velocity 

%discretization of cylinder boundary 
dt = 2*pi/N; 
t = 0:dt:2*pi-dt/2;
t = t';

%regularization parameter
ep = dt/4; 

%cylinder on which velocity is set - cylinder is of radius 1 
y1 = cos(t);
y2 = sin(t); 

%tangential of cylinder boundary - needed for exact solution of forces 
yp1 = -sin(t); 
yp2 = cos(t); 

%velocity on cylinder boundary 
u1 = -cos(2*t)/8 + 5*cos(4*t)/16 - cos(4*t)/4; 
u2 = sin(2*t)/8 + 5*sin(4*t)/16 - sin(4*t)/4; 

%computing the force 
f = RegStokeslets2D_velocitytoforce([y1,y2],[y1,y2],[u1,u2],ep,mu);

%f is a force and to compare to exact solution need force density - divide
%by radius*dt 
f1 = f(:,1)/dt;
f2 = f(:,2)/dt; 

%% Computing error 

%exact solution 
exactf1 = 2*sin(3*t).*yp1; 
exactf2 = 2*sin(3*t).*yp2; 

%computing error 
error1 = abs(f1-exactf1); 
error2 = abs(f2-exactf2); 

%prints the max error 
fprintf('maximum error in f1: %d \n',max(max(error1)));
fprintf('maximum error in f2: %d \n',max(max(error2)));

%% Plotting figures 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

figure(1) 
plot(t,f1,'k.-')
hold on
plot(t,exactf1,'r--')
plot(t,f2,'b.-')
plot(t,exactf2,'g--')
title('Numerical and Exact Solutions')
legend('f1','fexact1','f2','fexact2')
xlabel('t: angle around cylinder')

figure(2) 
plot(t,error1,'k.-')
hold on 
plot(t,error2,'b.-')
title('Error')
legend('f1 error','f2 error')
xlabel('t: angle around cylinder')
