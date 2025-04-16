function [rate] = rate_calculator(para, H, Rx, f)
%MUSIC algorithm
%   [rate] = rate_calculator(para, H, Rx, f)
%Inputs:
%   para: structure of the initial parameters
%   H: communication channels
%   Rx: covariance matrix of transmit signal
%   f: beamformers of communication signals
%Outputs:
%   rate: communication rates for all users
%Date: 16/04/2025
%Author: Zhaolin Wang

rate = zeros(para.K, 1);
for k = 1:para.K
    hk = H(:,k);
    fk = f(:,k);
    SINR_k = abs(hk.'*fk)^2 / ( real(hk.'*Rx*conj(hk)) -  abs(hk.'*fk)^2 + 1);
    rate(k) = log2(1 + SINR_k);
end


end

