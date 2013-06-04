function [xest_updated, P_updated, K, S] = kf_corr (xest, y, N, P, H)

% x_hat = state estimate
% K - Kalman Gain
% S - Innovation covariance
% N - Measurement Noise covariance
% P - Estimate covariance



S = H*P*H' + N;
K = P*H'*inv(S);          % compute Kalman Gain
xest_updated = xest+ K*(y-H*xest);
P_updated = P - K*S*K';
