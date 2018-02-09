function [phi, theta] = torus_inversion(x,y,z,r,R)
%TORUS_INVERSION
% [phi, theta] = torus_inversion(x,y,z,r,R)
% For the three input vectors x,y,z, returns the corresponding values of
% angles  on the torus parametrized as
%     [ (R+r*cos(theta)) *  cos(phi); 
%       (R+r*cos(theta)) * sin(phi); 
%               r*sin(theta) ]

Rtilda = sqrt(x.^2 + y.^2);
theta = atan2(z, (Rtilda - R)/r);
phi = atan2(y./Rtilda, x./Rtilda);