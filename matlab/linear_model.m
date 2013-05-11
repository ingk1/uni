function [x] = linear_model(x1, P0, F, Q, G, N) 

% x1 is the initial state [range; velocity; acceleration]
% F is the target motion model matrix e.g. [speed; velocity; acceleration]
% G is the target motion model noise matrix 
% H is the observation matrix
% n is the 
% w is the measurement noise (standard deviation)
% N is the number of iterations

% Generate target flight (with constant acceleration)
x = zeros(3, N+1);        % target state vector (range; range-rate; acceleration)
% y = zeros(2, N);        % measurement vector  (range; range-rate)


% Generate process noise

%sigma = sqrt(w);                    % target model process noise(stdev of target accl (m/s^2))
% w = zeros(3, N);                    % target trajectory noise
w = (G*chol(Q))'*randn(3,N);
% w = chol(Q)'*randn(3,N);              % generate acceleration noise ~N(0,sigma^2) (iid)
% w(1:2,:) = [dt^2/2; dt]*w(3,:);     % calc range and range rate noise

% Generate measurement noise
%sigma = n;                          % measurement model noise [phase noise, frequency shift noise]
%n = diag(sigma)*randn(2,N);         % measurement noise vector
% n = chol(R)'*randn*2,N);

% Target State Model
x(:,1) = x1 + chol(P0)'*randn(3,1);%w(:,1);             % target state vector (r, rr, ra)
    
for k = 2:N+1    
    x(:,k) = F*x(:,k-1) + w(:,k-1); % target state vector (r, rr, ra)
%     y(:,k-1) = H*x(:,k) + n(:,k-1);   % measurement vector of target state (r, rr)
end

