function [A, B] = reg_fncs(ep,R)

% Computes two components of blobs (aka H1 and H2 in RC's notes) 
% for the Method of Regularized Stokeslets in 3D 
% Based on Cortez, Fauci, Medovikov, Physics of Fluids 2005 

% Developed by Shilpa Khatri and Ricardo Cortez 
% June 2024 

%ep: blob width (regularization parameter)
%R: distance between targe and source point + regularization 
%   (R = sqrt(|x-y|^2 + ep^2))

% Blob as given in Cortez, Fauci, Medovikov, Physics of Fluids 2005, Eqn 9 
A = (R.^2 + ep.^2) ./ (R.^3);
B = 1./ (R.^3);

% Another blob in 3D 
%R2 = R.^2; 
%R3 = R.^3; 
%d2 = ep.^2; 

%B = (2*R3.^2 +3*d2*R2.^2+35*d2^3)./(2*R3.^3) ;
%A = (R2 + d2).*B -2*d2*(15*d2^2+R2.^2)./(R3.*R2.^2) ;