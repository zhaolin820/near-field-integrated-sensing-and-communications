function [Rx, f] = SDR_fully_digital(para, H, beta_s, scale)
%MUSIC algorithm
%   [spectrum, X, Y] = MUSIC_estimation(para, Rx, f, G)
%Inputs:
%   para: structure of the initial parameters
%   H: communication channels
%   beta_s: complex sensing channel gain
%   scale: a scaling factor to facilitate optimization
%Outputs:
%   Rx: covariance matrix of transmit signal
%   f: beamformers of communication signals
%Date: 14/06/2023
%Author: Zhaolin Wang


cvx_begin
    % optimization variables
    variable F(para.N, para.N, para.K) complex
    variable U(2, 2) semidefinite
    variable Rx(para.N, para.N) complex semidefinite

    [J_11, J_12, J_22] = FIM(para, Rx, beta_s, scale);
    [J_11 - U, J_12; J_12', J_22] == hermitian_semidefinite(4);

    for k = 1:para.K
        hk = H(:,k);
        Fk = F(:,:,k);
        Fk == hermitian_semidefinite(para.N);
        quad_form(conj(hk), Fk) >= (2^para.Rmin - 1) *(quad_form(conj(hk), Rx) - quad_form(conj(hk), Fk) + 1);
    end

    real(trace(Rx)) <= para.Pt;

    Rx - sum(F,3) == hermitian_semidefinite(para.N);
    
    obj = trace_inv(U);
    minimize(obj);  
cvx_end   

% construct rank-one solution
f = zeros(para.N, para.K);
for k = 1:para.K
    hk = H(:,k);
    Fk = F(:,:,k);
    f(:,k) = (hk.'*Fk*conj(hk))^(-1/2) * Fk* conj(hk);
end

end

