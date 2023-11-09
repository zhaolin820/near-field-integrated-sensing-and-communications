function [spectrum, X, Y] = MUSIC_estimation(para, Rx, f, G)
%MUSIC algorithm
%   [spectrum, X, Y] = MUSIC_estimation(para, Rx, f, G)
%Inputs:
%   para: structure of the initial parameters
%   Rx: covariance matrix of transmit signal
%   f: beamformers of communication signals
%   G: traget response matrix
%Outputs:
%   spectrum: MUSIC spectrum
%   X, Y: coordinations
%Date: 14/06/2023
%Author: Zhaolin Wang


%% Generate signals
Rs = Rx - f*f';
A = matrix_decomposition(Rs);
L = chol(A*A');
s = L'*(randn(para.N,para.T)+1i*randn(para.N,para.T))/sqrt(2); % dedicated sensing signal
c = (randn(para.K, para.T)+1i*randn(para.K, para.T))/sqrt(2); % communication signal
N_s = sqrt(1/2*para.noise) * ( randn(para.N, para.T) + 1i*randn(para.N, para.T) ); % noise

X = f*c + s; % transmit signal
Y_s = G*X + N_s; % receives signal


%% MUSIC algorithm
R = 1/para.T * (Y_s*Y_s');
[U,D] = eig(R);
[~, ind] = sort(diag(D), 'descend');
U = U(:, ind);
Uz = U(:,2:end);
U = Uz*Uz';

m = 500;
X = linspace(0, 40, m);
Y = linspace(0, 40, m);
[X, Y] = meshgrid(X, Y);
[theta_all, r_all] = cart2pol(X, Y);
spectrum = zeros(length(r_all), length(theta_all));
for i = 1:m
    parfor j = 1:m
        aa = beamfocusing(para, r_all(i,j), theta_all(i,j));
        spectrum(i,j) = 1/real( aa' * U * aa );
    end
end
spectrum = spectrum ./ max(max(spectrum));

end

function [U_de] = matrix_decomposition(U)
    [P,C] = eig(U) ;
    U_de = P * (C.^(1/2)) * P';
end
