function [u] = RegStokeslets2D_forcetovelocity(y,f,x,ep,mu)

% Computes velocities at set of points when given a set of points and
% forces at those points using the Method of Regularized Stokeslets 
% Based on Cortez, SIAM J. Sci Comput. 2001
% Constant flow is set to 0 here (Uo in Eqn 9)

% Developed by Shilpa Khatri and Ricardo Cortez 
% July 2024 

%y = (y1,y2) source points
%f = (f1,f2) forces at those source points (not force density)
%x = (x1,x2) target points 
%u = (u1,u2) velocity evaluated at those target points 
%mu is the viscosity 
%ep is the width of the regularization 

N = size(y,1); %number of source points 
M = size(x,1); %number of target points 

%unpacking the inputs 
y1 = y(:,1);
y2 = y(:,2);
f1 = f(:,1);
f2 = f(:,2); 
x1 = x(:,1); 
x2 = x(:,2); 

%initializing the velocity 
u1 = zeros(M,1);
u2 = zeros(M,1); 

%loop over source points    
for k = 1:N 

    %distance between target and source points 
    XY1 = x1(:) - y1(k); 
    XY2 = x2(:) - y2(k);  
    R2 = XY1.^2 + XY2.^2 + ep^2; 
    R = sqrt( R2 ); 

    %computing the velocity 
    [A, B] = reg_fncs(ep,R);

    fdotXY = f1(k)*XY1 + f2(k)*XY2; 

    u1(:) = u1(:) + f1(k)*A + fdotXY.*B.*XY1;  
    u2(:) = u2(:) + f2(k)*A + fdotXY.*B.*XY2;

end

u1 = u1/4/pi/mu; 
u2 = u2/4/pi/mu; 

%repacking output 
u = [u1 u2]; 


