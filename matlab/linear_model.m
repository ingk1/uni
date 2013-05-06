function [x, y] = linear_model(x1, F, G, H, w, n, N)

% x1 is the initial state [range; velocity; acceleration]
% F is the target motion model matrix e.g. [speed; velocity; acceleration]
% G is the target motion model noise matrix 
% H is the observation matrix
% n is the 
% w is the measurement noise (standard deviation)
% N is the number of iterations

% Generate target flight (with constant acceleration)
x = zeros(3, N);        % target state vector (range; range-rate; acceleration)
y = zeros(2, N);        % measurement vector  (range; range-rate)


% Generate process noise

sigma = sqrt(w);                    % target model process noise(stdev of target accl (m/s^2))
% w = zeros(3, N);                    % target trajectory noise
w = diag(sigma)*randn(3,N);              % generate acceleration noise ~N(0,sigma^2) (iid)
% w(1:2,:) = [dt^2/2; dt]*w(3,:);     % calc range and range rate noise

% Generate measurement noise
sigma = n;                          % measurement model noise [phase noise, frequency shift noise]
n = diag(sigma)*randn(2,N);         % measurement noise vector

% Target State Model
x(:,1) = x1 + G*w(:,1);             % target state vector (r, rr, ra)
y(:,1) = H*x(:,1) + n(:,1);         % measurement vector of target state (r, rr)
    
for k = 1:N-1    
    x(:,k+1) = F*x(:,k) + G*w(:,k); % target state vector (r, rr, ra)
    y(:,k+1) = H*x(:,k) + n(:,k);   % measurement vector of target state (r, rr)
end

