function [N_theta] = N_cov(theta, T, r, SNR_1000)
%
% This function computes the Measurement Noise Covariance Matrix for an
% input waveform library.
%
% INPUTS:
% theta: The waveform library vector
% T:    Measurement transition matrix
% r:  	Target range (m);
% SNR_1000: Target RX SNR at 100om (dB)
% c:    Signal velocity in TX medium
% wc:   Carrier frequency (rad/s)

% OUTPUTS:
% N_theta: Vector of Noise Covariance Matrices for the waveform library



% %%  Compute Noise Covariance for Waveforms
% SNR = 0; % 0dB SNR
% c = 1500; % speed of signal in medium (m/s) w
% wc = 2+pi*25E3; % 25KHz TX centre frequency
% T = [c/2 0; 0 c/(2*wc)]; % transformation matrix b/w rx estimator and tracking system measurement vector

eta_1000 = 10^(SNR_1000/10);
eta = (1000/r)^4 * eta_1000;
SNR = 10*log10(eta);

K = size(theta,1);      % number of unique waveforms
N_theta = zeros(2,2,K);       % initialise waveform noise covariance matrix

waveform_class = theta(:,1);
lambda = theta(:,2);
chirp_rate = theta(:,3);


parfor k = 1:K
    
    wc = waveform_class(k);
    l = lambda(k);
    b = chirp_rate(k);
    
    J = fisher(wc, l, b, SNR);
    N_theta(:,:,k) = T/J*T';             % calc noise covariance matrices
end

