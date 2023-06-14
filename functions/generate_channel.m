function [H, G, beta_s, r, theta, r_s, theta_s] = generate_channel(para)
%Generate communication and sensing channels
%  [H, G, beta_s, r, theta, r_s, theta_s] = generate_channel(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   H: communication channels
%   G: target response matrix
%Date: 14/06/2023
%Author: Zhaolin Wang

%% communication channels
% generate user location
lambda = para.c/para.f;
Rayleigh_distance = 2*para.D^2/lambda;
r = rand(para.K, 1) * Rayleigh_distance;
theta = rand(para.K,1)*pi;

% generate channels
H = zeros(para.N, para.K);
for k = 1:para.K
    beta = sqrt(1/para.noise) * sqrt(para.rho_0)/r(k) * exp(-1i*2 * pi * para.f/para.c * r(k));
    H(:,k) = beta * beamfocusing(para, r(k), theta(k));
end

%% target response matrix
r_s = para.r_s;
theta_s = para.theta_s;

beta_reflection = sqrt(1/2) * (randn(1) + 1i*randn(1));
beta_s = sqrt(para.rho_0)/(2*r_s) * exp(-1i*2 * pi * para.f/para.c * 2*r_s) * beta_reflection;
a = beamfocusing(para, r_s, theta_s);
G = beta_s * (a*a.');

end