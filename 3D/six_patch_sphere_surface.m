function [Px,Py,Pz,dArea]=six_patch_sphere_surface(N,r)

% Discretizes the surface of a sphere of radius r using the six-patch structure
% the discretization uses an N x N uniform grid on the face of a cube 
% in which the sphere is embedded and then maps the points to the 
% surface of the sphere 
% Based on description given in Cortez, Fauci, Medovikov, Physics of Fluids 2005 

% Developed by Ricardo Cortez and Shilpa Khatri 
% July 2024 

%r: radius of sphere 
%N: number of points in one of the patches of the six (the patch is
%discretized as N x N 
%[Px, Py, Pz]: x, y, and z coordinates of points on the surface of the
%sphere 
%darea: gridsize of each point 

%discretization in 1D for a unit sphere 
h = 2/N; 

%When we discretize the cube to find the sphere discretization,
%it is better to have the points on the faces of the cube at the
%center of the grid cells.
x = (-1+h/2 : h : 1-h/2)';
[a1,a2] = meshgrid(x,x); 
a3 = ones(size(a1));

%endpoints of each discretized square 
am = a1(:)-h/2; bm = a2(:)-h/2;
ap = a1(:)+h/2; bp = a2(:)+h/2;

%the exact surface area of each projected grid cell on unit sphere
FF1 = 1/2*atan2(2*ap.*bp.*sqrt(1+ap.^2+bp.^2),1+bp.^2-ap.^2.*(-1+bp.^2));
FF2 = 1/2*atan2(2*am.*bp.*sqrt(1+am.^2+bp.^2),1+bp.^2-am.^2.*(-1+bp.^2));
FF3 = 1/2*atan2(2*ap.*bm.*sqrt(1+ap.^2+bm.^2),1+bm.^2-ap.^2.*(-1+bm.^2));
FF4 = 1/2*atan2(2*am.*bm.*sqrt(1+am.^2+bm.^2),1+bm.^2-am.^2.*(-1+bm.^2));
ddArea = (FF1-FF2)-(FF3-FF4);

%six patches of sphere concatenated 
xb=[a3(:);-a3(:); a1(:);a1(:); a1(:);a1(:)];
yb=[a1(:); a1(:);-a3(:);a3(:); a2(:);a2(:)];
zb=[a2(:); a2(:); a2(:);a2(:);-a3(:);a3(:)];
dArea = [ddArea;ddArea;ddArea;ddArea;ddArea;ddArea];

%rescaling points on surface of unit sphere 
bmag = sqrt( xb.^2 + yb.^2 + zb.^2 );
Px = (xb./bmag);
Py = (yb./bmag);
Pz = (zb./bmag);

%rescale points and area for sphere of radius r 
Px = r*Px; 
Py = r*Py; 
Pz = r*Pz; 
dArea = r*r*dArea; 

