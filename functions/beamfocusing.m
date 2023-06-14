function [a] = beamfocusing(para, r, theta)
%Near-field array response vector
%   [a] = beamfocusing(para, r, theta)
%Inputs:
%   para: structure of the initial parameters
%   r: distance
%   theta: direction
%Outputs:
%   a: array response vector
%Date: 14/06/2023
%Author: Zhaolin Wang



n = (-(para.N-1)/2 : (para.N-1)/2) * para.d;
n = n';

% distance
% r = -n*cos(theta) + n.^2*(sin(theta))^2/(2*r);
r = sqrt(r^2 + n.^2 - 2*r*n*cos(theta)) - r;

% beamfocusing vector
a = exp( -1i * 2 * pi * para.f/para.c * r );

end