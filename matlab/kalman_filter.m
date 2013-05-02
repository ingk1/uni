clear all;
close all;
clc;

c = 1500;   % speed of sound (m/s)
wc = 25E3*2*pi;       % TX centre frequnecy (f*2pi)

nk = 500; % max num iterations

% Kalman Filter

dt = 2;       % model update rate (s)
t = 0:dt:dt*(nk-1);  % time vector;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Target Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% F is the target state vector, for calculating x which includes range, range rate and radial
% acceleration.

a_corr = 0.833;       % acceleration correction coefficient
% a_corr = 1;

F = [1 dt dt^2/2; 0 1 dt; 0 0 a_corr];      % target motion model matrix (s, v, a)'
G = [0 0 1];
H = [1 0 0; 0 1 0];          % observation matrix


% Generate target flight (with constant acceleration)
x = zeros(3, nk);        % target state vector (range; range-rate; acceleration)
y = zeros(2, nk);        % measurement vector  (range; range-rate)


% Generate process noise
w = zeros(3, nk);                    % target trajectory noise
sigma = sqrt(0.01);                  % target model process noise(stdev of target accl (m/s^2))
w(3, :) = sigma*randn(1,nk);         % generate acceleration noise ~N(0,sigma^2)
w(1:2,:) = [dt^2/2; dt]*w(3,:);     % calc range and range rate noise



B = F(:,3);                         % system force (acceleration)
Q = B*B'*sigma^2;                   % process noise covariance Q = G*G'*sigma^2 such that w~N(0,Q)

% Generate measurement noise
sigma = [300 0.5];                   % measurement model noise [phase noise, frequency shift noise]
n = diag(sigma)*randn(2,nk);         % measurement noise vector


N = diag(sigma.^2);                 % measurement noise covariance N = diag(sigma1^2, sigma2^2) i.e. signals are independant



% T = diag(c/2, c/2*wc);
% N = 1/n * T * 1/U *T;

% Initialise target state vector x(k) for k=1
range = 500;         % range (m)
vel   = 5.555;       % radial velocity (m/s)
accl  = 0.5;         % radial acceleration (m/s^2)
u = [range vel accl]';              % initial target state

% Target State Model
x(:,1) = u + G*w(:,1);          % target state vector (r, rr, ra)
y(:,1) = H*x(:,1) + n(:,1);       % measurement vector of target state (r, rr)
    
for k = 1:nk-1    
    x(:,k+1) = F*x(:,k) + G*w(:,k);     % target state vector (r, rr, ra)
    y(:,k+1) = H*x(:,k) + n(:,k);       % measurement vector of target state (r, rr)
end




% Now for the Kalman Filter

% x_hat = state estimate
% K - Kalman Gain
% S - Innovation covariance
% N - Measurement Noise covariance
% P - Estimate covariance

% Q calc covariance of process noise (w) from initial value of w

x_hat = zeros(3, nk);


x_hat(:,1) = [y(:,1);a_corr];

% x_hat(1:3,1) = [1000;0;0];
% x_hat(3,1) = 1;

P = zeros(3,3);
mse = zeros(1,nk);      % Mean-Square Track Error
% 
for k = 1:nk-1
%     


S = H*P*H' + N;
K = P*H'*1/S;
x_hat(:,k) = x_hat(:,k)+ K*(y(:,k)-H*x_hat(:,k));
P = P -K*S*K';
% 
mse(k) = trace(P);             % calc the mean-square track error
x_hat(:,k+1) = F*x_hat(:,k);
P = F*P*F' + G*Q*G';        % only acceleration noise is 


%     
%     5
%     6
%     
%     1
%     2
%     3
%     4
end

t = t(1:nk-1);

figure(1)
range_target = x(1, 1:nk-1);
range_meas = y(1,1:nk-1);
range_est = x_hat(1, 2:nk);
plot(t, range_target); hold on
plot(t, range_meas, 'r'); hold on
plot(t, range_est, 'k'); hold on
legend('Target Range', 'Measured Range', 'Estimated Range')

figure(2)
vel_target = x(2,1:nk-1);
vel_meas = y(2,1:nk-1);
vel_est = x_hat(2, 2:nk);
plot(t, vel_target); hold on
plot(t, vel_meas, 'r'); hold on
plot(t, vel_est, 'k'); hold on
legend('Target Velocity', 'Measured Velocity', 'Estimated Velocity')

figure(3)
plot(t, mse(1:nk-1))

% figure(4)
% range_err  = abs(range_target(1:nk-1)-range_est(2:nk));
% 
% plot(t(1:nk-1),range_err)
% 
% 
% subfunction target_model (x1, nk)


figure(5)
range_err  = abs(range_target-range_est);

plot(t,range_err)

    