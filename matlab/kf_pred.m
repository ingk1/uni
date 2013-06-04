function [xest_pred, Ppred] = kf_pred (xest, P, F, Q, G)

% x_hat = state estimate
% K - Kalman Gain
% S - Innovation covariance
% N - Measurement Noise covariance
% P - Estimate covariance

xest_pred = F*xest;
Ppred = F*P*F' + G*Q*G';        % only acceleration noise is 

