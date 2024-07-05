function [f,ut,omega] = RegStokeslets2D_velocitytoforce_augmented(y,x,u,ep,mu)

% Computes forces at set of points and constant translational and rotational 
% background velocities when given a set of points and velocities at those 
% points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001

%The matrix is augmented so that the sum of the forces is zero and the sum
%of the torque is 0. 

% Developed by Shilpa Khatri and Ricardo Cortez 
% July 2024 

%y = (y1,y2) source points
%f = (f1,f2) forces at those source points (not force density) 
%x = (x1,x2) target points 
%u = (u1,u2) velocity at those target points 
%mu is the viscosity 
%ep is the width of the regularization 

N = size(y,1); %number of source points 
Nt = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
u1 = u(:,1);
u2 = u(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 

%initializing the forces
f1 = zeros(N,1);
f2 = zeros(N,1); 

%building matrix 

%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1(:,k) = x1(:) - y1(k); 
    XY2(:,k) = x2(:) - y2(k);  

end

R2 = XY1.^2 + XY2.^2 + ep^2; 
R = sqrt( R2 ); 

[A, B] = reg_fncs(ep,R); 

M11 = A + B.*XY1.*XY1; 
M22 = A + B.*XY2.*XY2; 
M12 = B.*XY1.*XY2; 
%M21 = B.*XY2.*XY1; 

M = [M11 M12; M12 M22]/4/pi/mu;

%Augmenting matrix 

%sum of forces = 0 
M3(1,1:N) = 1; 
M3(1,N+1:2*N) = 0; 
M3(2,1:N) = 0; 
M3(2,N+1:2*N) = 1; 

%sum of torque = 0 
M4(1:N) = -y2'; 
M4(N+1:2*N) = y1'; 

M = [M;M3;M4]; 

%to incorporate translational rigid body motion 
M5(1:Nt,1) = -1; 
M5(Nt+1:2*Nt,1) = 0;
M5(1:Nt,2) = 0; 
M5(Nt+1:2*Nt,2) = -1; 
M5 = [M5;0 0;0 0;0 0];

%to incorporate rotational rigid body motion 
M6 = [x2; -x1; 0; 0; 0];

M = [M M5 M6]; 

%solving for force and translational and rotational velocity 
%when given a velocity 
uu = [u1;u2;0;0;0];
ff = M\uu; 
f1 = ff(1:N); 
f2 = ff(N+1:2*N);
ut = ff(2*N+1:2*N+2);
omega = ff(end); 

%repacking force output 
f = [f1,f2]; 


