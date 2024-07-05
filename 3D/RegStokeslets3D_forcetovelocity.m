function [u] = RegStokeslets3D_forcetovelocity(y,f,x,ep,mu)

% Computes velocities at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
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
f1 = f(:,1);
f2 = f(:,2); 
f3 = f(:,3);
x1 = x(:,1); 
x2 = x(:,2);
x3 = x(:,3); 

%initializing the velocity 
u1 = zeros(M,1);
u2 = zeros(M,1); 
u3 = zeros(M,1);

%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2 = x2(:) - y2(k);  
    XY3 = x3(:) - y3(k); 
    R2 = XY1.^2 + XY2.^2 + XY3.^2 + ep^2; 
    R = sqrt( R2 ); 

    %computing the velocity 
    [A, B] = reg_fncs(ep,R);

    fdotXY = f1(k)*XY1 + f2(k)*XY2 + f3(k)*XY3;

    u1(:) = u1(:) + f1(k)*A + fdotXY.*B.*XY1;  
    u2(:) = u2(:) + f2(k)*A + fdotXY.*B.*XY2;
    u3(:) = u3(:) + f3(k)*A + fdotXY.*B.*XY3;

end

u1 = u1/8/pi/mu; 
u2 = u2/8/pi/mu; 
u3 = u3/8/pi/mu; 

%repacking output 
u = [u1 u2 u3]; 


