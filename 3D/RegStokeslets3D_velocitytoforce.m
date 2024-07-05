function [f] = RegStokeslets3D_velocitytoforce(y,x,u,ep,mu)

% Computes forces at set of points when given a set of points and
% velocities at those points using the Method of Regularized Stokeslets 
% Based on Cortez, Fauci, Medovikov, Physics of Fluids 2005  

%Developed by Shilpa Khatri and Ricardo Cortez 
%July 2024 

%y = (y1,y2,y3) source points
%f = (f1,f2,f3) forces at those source points (not force density)
%x = (x1,x2,x3) target points 
%u = (u1,u2,u3) velocity evaluated at those target points 
%mu is the viscosity 
%ep is the width of the regularization 

N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);
u1 = u(:,1);
u2 = u(:,2);
u3 = u(:,3); 
x1 = x(:,1); 
x2 = x(:,2); 
x3 = x(:,3); 

%initializing the forces
f1 = zeros(N,1);
f2 = zeros(N,1); 
f3 = zeros(N,1); 

%building matrix 

%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1(:,k) = x1(:) - y1(k); 
    XY2(:,k) = x2(:) - y2(k);  
    XY3(:,k) = x3(:) - y3(k); 

end

R2 = XY1.^2 + XY2.^2 + XY3.^2 + ep^2; 
R = sqrt( R2 ); 

[A, B] = reg_fncs(ep,R); 

M11 = A + B.*XY1.*XY1; 
M22 = A + B.*XY2.*XY2; 
M33 = A + B.*XY3.*XY3; 
M12 = B.*XY1.*XY2; 
M13 = B.*XY1.*XY3; 
M23 = B.*XY2.*XY3; 

M = [M11 M12 M13; M12 M22 M23; M13 M23 M33]/8/pi/mu;

%solving for force when given a velocity 
uu = [u1;u2;u3];
ff = M\uu; 
f1 = ff(1:N); 
f2 = ff(N+1:2*N);
f3 = ff(2*N+1:end); 

%repacking output 
f = [f1,f2,f3]; 

