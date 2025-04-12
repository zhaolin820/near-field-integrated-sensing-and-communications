clc
clear all
close all

addpath('./functions');
cvx_solver mosek


para = para_init();

[H, G, beta_s, r, theta, r_s, theta_s] = generate_channel(para);
scale = 1e2; % if the optimization is failed, try to adjust this scale factor

% Optimize the transmit waveform
[Rx, f] = SDR_fully_digital(para, H, beta_s, scale);


% Calculate the CRB matrix
[J_11, J_12, J_22] = FIM(para, Rx, beta_s, scale);
scale = scale*para.noise/para.T;
J_11 = J_11/scale;
J_12 = J_12/scale;
J_22 = J_22/scale;
CRB = real(inv(J_11 - J_12*inv(J_22)*J_12'));

% MUSIC algorithm
[spectrum, X, Y] = MUSIC_estimation(para, Rx, f, G);

figure; colormap jet;
mesh(X,Y,10*log10(spectrum)); 
view([60,60]);
xlim([0,40]);
ylim([0,40]);
xlabel('x (m)'); ylabel('y (m)'); zlabel('Spectrum (dB)');



