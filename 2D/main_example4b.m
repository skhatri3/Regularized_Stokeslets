%main_example4b  

%Example 4b Section 4.2 of Cortez, SIAM J. Sci Comput. 2001 
%Setting tangential forces on a cylinder of radius 1 and computing the
%velocity at given points 

%Developed by Shilpa Khatri and Ricardo Cortez  
%July 2024 

clear all
close all 

%% Parameters to set 

%setting the viscosity
mu = 1; 

%domain on which velocity is computed 
x1min = 0; 
x1max = 2; 
%x2min and x2max to compare with plots in paper 
%x2min = 3/10; 
%x2max = 3/10; 
%x2min and x2max to see a 2d domain 
x2min = 0; 
x2max = 2; 

%number of points on boundary where force is applied 
N = 100;  

%resolution for velocity grid
Nx1 = 80;  
%Nx2 = 1; %corresponds to compare with plots in paper above 
Nx2 = 80; %for a 2d domain 

%% Setting forces and computing velocity 

%discretization of cylinder boundary 
dt = 2*pi/N; 
t = 0:dt:2*pi-dt/2;
t = t';

%regularization parameter
ep = dt/4; 

%cylinder on which forces are applied - cylinder is of radius 1 
y1 = cos(t); 
y2 = sin(t); 

%tangential of cylinder boundary
yp1 = -sin(t);
yp2 = cos(t);

%forces on cylinder boundary 
%note that the force density is given in the paper - multiply by radius*dt  
f1 = 2*sin(3*t).*yp1*dt; 
f2 = 2*sin(3*t).*yp2*dt;

%points on which velocity will be computed 
xx1 = linspace(x1min,x1max,Nx1);
xx2 = linspace(x2min,x2max,Nx2); 
[x1m,x2m] = ndgrid(xx1,xx2); 
x1 = x1m(:);
x2 = x2m(:);

%computing velocity 
u = RegStokeslets2D_forcetovelocity([y1,y2],[f1,f2],[x1,x2],ep,mu);
u1 = u(:,1);
u2 = u(:,2); 
u1m = reshape(u1,size(xx1,2),size(xx2,2)); 
u2m = reshape(u2,size(xx1,2),size(xx2,2));

%% Computing error 

%exact solution 

for i = 1:length(xx1)

    for j = 1:length(xx2)

        r = sqrt(x1m(i,j).^2 + x2m(i,j).^2); %radius 
        s = atan2(x2m(i,j),x1m(i,j)); %angle

        if (r < 1)
    
            uexact1(i,j) = cos(2*s)*r^2/8 + cos(4*s)*r^4/16 - cos(2*s)*r^4/4; 
            uexact2(i,j) = -sin(2*s)*r^2/8 + sin(4*s)*r^4/16 + sin(2*s)*r^4/4; 

        else

           uexact1(i,j) = -cos(2*s)/(r^2)/8 + 5*cos(4*s)/(r^4)/16 - cos(4*s)/(r^2)/4; 
           uexact2(i,j) = sin(2*s)/(r^2)/8 + 5*sin(4*s)/(r^4)/16 - sin(4*s)/(r^2)/4; 

        end

    end

end

%computing error 
error1 = abs(u1m-uexact1);
error2 = abs(u2m-uexact2); 
errormag = sqrt(error1.^2 + error2.^2); 

%prints the max error 
fprintf('maximum error in u1: %d \n',max(max(error1)));
fprintf('maximum error in u2: %d \n',max(max(error2)));

%% Plotting figures 
skip = 4; %for quiver plots 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

%plots as in paper, velocity is computed on a line, x2 = constant 
if (Nx2 == 1) 

    figure(1) 
    plot(x1,u1,'k.-')
    hold on
    plot(x1,uexact1,'r--')
    plot(x1,u2,'b.-')
    plot(x1,uexact2,'g--')
    title('Numerical and Exact Solutions')
    legend('u1','uexact1','u2','uexact2')
    xlabel('x1')

    figure(2) 
    plot(x1,error1,'k.-')
    hold on 
    plot(x1,error2,'b.-')
    title('Error')
    legend('u1 error','u2 error')
    xlabel('x1')
    
%plots where the velocity is computed in a 2d domain 
else 

    figure(1) 
    plot(y1,y2,'k.-')
    hold on 
    quiver(y1,y2,f1,f2,'r')
    axis equal 
    quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end),u2m(1:skip:end,1:skip:end),'k')
    xlim([x1min,x1max])
    ylim([x2min,x2max])
    title('Forces and Computed Velocity')
    
    figure(2) 
    plot(y1,y2,'k.-')
    hold on 
    quiver(y1,y2,f1,f2,'r')
    axis equal 
    quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),uexact1(1:skip:end,1:skip:end),uexact2(1:skip:end,1:skip:end),'k')
    xlim([x1min,x1max])
    ylim([x2min,x2max])
    title('Forces and Exact Velocity')
    
    figure(3) 
    plot(y1,y2,'k.-')
    hold on 
    quiver(y1,y2,f1,f2,'r')
    axis equal 
    quiver(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),error1(1:skip:end,1:skip:end),error2(1:skip:end,1:skip:end),5,'k')
    xlim([x1min,x1max])
    ylim([x2min,x2max])
    title('Forces and Error in the Velocity')
    
    figure(4)  
    pcolor(x1m,x2m,log10(errormag)+eps)
    shading interp
    colorbar
    caxis( [-5 -3] )
    hold on
    plot(y1,y2,'k.-') 
    quiver(y1,y2,f1,f2,'r')
    axis equal 
    xlim([x1min,x1max])
    ylim([x2min,x2max])
    title('Forces and Error Magnitude in the Velocity')

end



