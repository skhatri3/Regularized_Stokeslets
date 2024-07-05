%main_example1_velocitytoforce

% Example 1 Section 4 of Cortez, Fauci, Medovikov, Physics of Fluids 2005 
% Setting forces on a sphere of radius a and computing the
% velocity at given points 
% Velocity on sphere is set to (0,0,1) and the force is computed 

% Developed by Shilpa Khatri and Ricardo Cortez
% July 2024 

clear all
close all 

%% Parameters to set 

%setting the viscosity
mu = 1; 

%radius of sphere 
a = 2; 
%surface area of sphere
sa = 4*pi*a*a;

%sphere will be discretized using an N x N uniform grid on the face of a
%cube in which the sphere is embedded and then mapping the points to the
%surface of the sphere (see function sphere_surface)
%number of points N for the N x N grid on each of the six sphere patches
N = 24;

%velocity of translating sphere 
uo1 = 0; 
uo2 = 0; 
uo3 = 1; 

% 2d surface (x1-x3 plane) on which velocity is computed 
x1min = -5; 
x1max = 5; 
x2fixed = 0; 
x3min = -5; 
x3max = 5; 

%resolution for velocity grid
Nx1 = 200;  
Nx3 = 200; 

%% Setting forces and computing velocity 

%sphere points on which forces are applied - discretization uses an 
% N x N uniform grid on the face of a cube in which the unit sphere is embedded 
% and then maps the points to the surface of the sphere 
% darea - discretization (grid size) on sphere surface 
[y1, y2, y3, darea] = six_patch_sphere_surface(N,a); 

%regularization parameter
ep = mean(sqrt(darea))/2; 

%total  number of points on surface 
Npts = 6*N*N; 

%velocity on surface of sphere 
us1 = uo1*ones(Npts,1); 
us2 = uo2*ones(Npts,1); 
us3 = uo3*ones(Npts,1); 

%computing velocity on surface of sphere 
fs = RegStokeslets3D_velocitytoforce([y1,y2,y3],[y1,y2,y3],[us1,us2,us3],ep,mu);
%fs = M\[us1;us2;us3];
%f1 = fs(1:Npts); 
%f2 = fs(Npts+1:2*Npts);
%f3 = fs(2*Npts+1:end); 
f1 = fs(:,1);
f2 = fs(:,2); 
f3 = fs(:,3); 

%points on which velocity will be computed 
xx1 = linspace(x1min,x1max,Nx1);
xx3 = linspace(x3min,x3max,Nx3); 
[x1m,x3m] = ndgrid(xx1,xx3); 
x2m = x1m*0 + x2fixed; 
x1 = x1m(:);
x2 = x2m(:);
x3 = x3m(:);

%discretization of plane on which velocity will be computed 
h1 = (x1max-x1min)/Nx1; 
h3 = (x3max-x3min)/Nx3; 

%area of plane on which velocity will be computed 
pa = (x3max-x3min)*(x1max-x1min);

%computing velocity 
u = RegStokeslets3D_forcetovelocity([y1,y2,y3],[f1,f2,f3],[x1,x2,x3],ep,mu);
u1 = u(:,1);
u2 = u(:,2); 
u3 = u(:,3); 
u1m = reshape(u1,size(xx1,2),size(xx3,2)); 
u2m = reshape(u2,size(xx1,2),size(xx3,2));
u3m = reshape(u3,size(xx1,2),size(xx3,2));


%% Computing error of forces on surface of sphere 

%solution on surface 
%forces on fluid from sphere boundary dependent on translating velocity 
%note that the force density (traction) is given in the paper - multiply by darea  
fexact1 = 3*mu*uo1*ones(Npts,1).*darea/2/a; 
fexact2 = 3*mu*uo2*ones(Npts,1).*darea/2/a; 
fexact3 = 3*mu*uo3*ones(Npts,1).*darea/2/a;  

%computing error 
errors1 = abs(f1-fexact1);
errors2 = abs(f2-fexact2); 
errors3 = abs(f3-fexact3); 
errorsmag = sqrt(errors1.^2 + errors2.^2 + errors3.^2); 
l2errors1 = sqrt(sum(errors1.^2.*darea)/sa); 
l2errors2 = sqrt(sum(errors2.^2.*darea)/sa); 
l2errors3 = sqrt(sum(errors3.^2.*darea)/sa); 

fprintf('Error in force on surface of sphere \n')
%prints the max error 
fprintf('maximum error in f1: %d \n',max(errors1));
fprintf('maximum error in f2: %d \n',max(errors2));
fprintf('maximum error in f3: %d \n',max(errors3));

%prints the l2 error 
fprintf('l2 error in f1: %d \n',l2errors1);
fprintf('l2 error in f2: %d \n',l2errors2);
fprintf('l2 error in f3: %d \n',l2errors3);


%% Computing error on x1-x3 plane

%exact solution on x1-x3 plane - assumes velocity on sphere boundary is
%(0,0,uo3)

for i = 1:length(xx1)

    for j = 1:length(xx3)

      r = sqrt(x1m(i,j).^2 + x2m(i,j).^2 + x3m(i,j).^2); %radius 

        if (r < a)
    
            uexact1(i,j) = uo1; 
            uexact2(i,j) = uo2;
            uexact3(i,j) = uo3; 

        else

           uexactc =  (3*a*uo3/4) * ( (1/(r^3)) - ((a^2)/(r^5)) );
           uexactc2 = (a*uo3/4/r) * ( 3 + (a^2)/(r^2) );
           uexact1(i,j) = uexactc*x1m(i,j)*x3m(i,j); 
           uexact2(i,j) = uexactc*x2m(i,j)*x3m(i,j); 
           uexact3(i,j) = uexactc*x3m(i,j)*x3m(i,j) + uexactc2; 

        end

    end

end

%computing error 
error1 = abs(u1m-uexact1);
error2 = abs(u2m-uexact2); 
error3 = abs(u3m-uexact3);
errormag = sqrt(error1.^2 + error2.^2 + error3.^2); 
l2error1 = sqrt(sum(sum(error1.^2*h1*h3))/pa); 
l2error2 = sqrt(sum(sum(error2.^2*h1*h3))/pa); 
l2error3 = sqrt(sum(sum(error3.^2*h1*h3))/pa); 

fprintf('Error in velocity on 2D plane \n')
%prints the max error 
fprintf('maximum error in u1: %d \n',max(max(error1)));
fprintf('maximum error in u2: %d \n',max(max(error2)));
fprintf('maximum error in u3: %d \n',max(max(error3)));

%prints the l2 error 
fprintf('l2 error in u1: %d \n',l2error1);
fprintf('l2 error in u2: %d \n',l2error2);
fprintf('l2 error in u3: %d \n',l2error3);

%% Plotting figures 
skip = 8; %for quiver plots 
set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',2.0,...
      'defaultlinelinewidth',2.0,'defaultlinemarkersize',10.0)

%plots on surface points 

figure(1)
plot(f1,'k-')
hold on 
plot(f2,'b-')
plot(f3,'r-') 
plot(fexact1,'k--') 
plot(fexact2,'b--')
plot(fexact3,'r--') 
xlabel('surface points')
legend('f1','f2','f3','f1exact','f2exact','f3exact')
title('Numerical and Exact Solutions on Surface Points')

figure(2)
plot(errors1,'k-')
hold on 
plot(errors2,'b-')
plot(errors3,'r-') 
xlabel('surface points')
legend('f1 error','f2 error','f3 error')
title('Error on Surface Points')

%plotting on 2D surface in domain 

figure(3)
plot3(y1,y2,y3,'.')
hold on 
quiver3(y1,y2,y3,f1,f2,f3)
axis equal 
splot = surf(x1m,x2m,x3m,uexact3 - uo3);
splot.EdgeColor = 'none';
colorbar
quiver3(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),x3m(1:skip:end,1:skip:end),uexact1(1:skip:end,1:skip:end)-uo1 ,uexact2(1:skip:end,1:skip:end)-uo2,uexact3(1:skip:end,1:skip:end)-uo3,'k')
view(0,0)
title('Exact Solution')

figure(4)
plot3(y1,y2,y3,'.')
hold on 
quiver3(y1,y2,y3,f1,f2,f3)
axis equal 
splot = surf(x1m,x2m,x3m,u3m - uo3);
splot.EdgeColor = 'none';
colorbar
quiver3(x1m(1:skip:end,1:skip:end),x2m(1:skip:end,1:skip:end),x3m(1:skip:end,1:skip:end),u1m(1:skip:end,1:skip:end)-uo1,u2m(1:skip:end,1:skip:end)-uo2,u3m(1:skip:end,1:skip:end)-uo3,'k')
view(0,0)
title('Numerical Solution')

figure(5)
plot3(y1,y2,y3,'.')
hold on 
quiver3(y1,y2,y3,f1,f2,f3)
axis equal 
splot = surf(x1m,x2m,x3m,log10(errormag)+eps);
splot.EdgeColor = 'none';
colorbar
view(0,0)
title('Error')


