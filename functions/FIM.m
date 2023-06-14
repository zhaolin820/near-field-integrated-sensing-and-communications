function [J_11, J_12, J_22] = FIM(para, Rx, beta, scale)
%Calculate the FIM matrix
%  [J_11, J_12, J_22] = FIM(para, Rx, beta, scale)
%Inputs:
%   para: structure of the initial parameters
%   Rx: covariance matrix of transmit signal
%   beta: complex channel gain
%   scale: a scaling factor to facilitate optimization
%Outputs:
%   J_11, J_12, J_22: elements of FIM
%Date: 14/06/2023
%Author: Zhaolin Wang


lambda = para.c/para.f;
r = para.r_s;
theta = para.theta_s;
a = beamfocusing(para, r, theta);
G = a*a.';

%% Calculate derivative of G
n = (-(para.N-1)/2 : (para.N-1)/2) * para.d;
n = n';
r_n = sqrt(r^2 + n.^2 - 2*r*n*cos(theta));
a_r = -1i*2*pi/lambda * ((r - n*cos(theta))./ r_n - 1) .*  a;
a_theta = -1i*2*pi/lambda * (r*n*sin(theta)) ./ r_n .* a;
G_r = 2*a_r*a.';
G_theta = 2*a_theta*a.';

%% Calculate elements of FIM
J_11 = scale * 2*abs(beta)^2 * real([ trace(G_r*Rx*G_r'), trace(G_r*Rx*G_theta'); ...
    trace(G_r*Rx*G_theta'), trace(G_theta*Rx*G_theta') ]);

J_12 = scale * 2*conj(beta) * real([ trace(G*Rx*G_r'), 1i*trace(G*Rx*G_r'); ...
    trace(G*Rx*G_theta'), 1i*trace(G*Rx*G_theta')]);

J_22 =scale * eye(2) * 2 * real(trace( G*Rx*G' ));
end
