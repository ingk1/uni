clear all;
close all;
clc;

c = 1500;               % speed of signal in medium (m/s)
wc = 25E3*2*pi;         % 25KHz TX centre frequency
nk = 31;                % max num iterations
dt = 2;                 % model update rate (s)
SNR_1000 = 0;           % Target RX reference SNR at 1000m (dB)

% Define the Waveform Library
waveform_class = [1 2 3 1 2];
upw = [0.01 0.3 250;
    0.01/7.4338 0.03/7.4338 50;
    0.01/7.4338 0.06/7.4338 50;
    11.83E-3/2 11.83E-3/2 1;
    12E-3/7.4338 12E-3/7.4338 1];
sweep_rate = [0 0 1;
    0 0 1;
    -15E3*2*pi 0 50;
    0 0 1;
    0 0 1];
debug = 0;

theta = waveform_lib(waveform_class, upw, sweep_rate, debug);
N = sum(upw(:,3).*sweep_rate(:,3));     % size of waveform library

% OPSW 1 = Minimisation of Mean Squre Tracking Error
% OPSW 2 = Minimisation of Validation Gate Volume
% OPSW 3 = Minimisation of Validation Gate Volume (closed form solution)
opsw = 3;                       %   optimal wavefrom selection technique

n_con = 2;                      % number of predefined conventional waveforms in the library
switch opsw
    case 1
        N = upw(1,3)+n_con;
        con_wfm = N-1;                       % conventional waveform
        theta = [theta(1:N-n_con,:); theta(end-n_con+1:end,:)];
        opt = 1:upw(1,3);                    % optimal waveform library
    case 2
        con_wfm = N;                       % conventional waveform
        opt = 1:N-n_con;                % optimal waveform library
    case 3
        con_wfm = N;                       % conventional waveform
        N = nk+n_con;                   % define library size
        theta = [theta(con_wfm,:);zeros(nk-1,3); theta(end-n_con+1:end,:)];
        con_wfm = N;                     % update conventional waveform with smaller library
        opt = 1:nk;                        % optimal waveform library
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Target Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F is the target state vector, for calculating x which includes range, range rate and radial
% acceleration.
a_corr = 0.833;       % acceleration correction coefficient
% a_corr = 1;

F = [1 dt dt^2/2; 0 1 dt; 0 0 a_corr];      % target motion model matrix (s, v, a)'
G = [0 0 1]';
H = [1 0 0; 0 1 0];          % observation matrix
T = [c/2 0; 0 c/(2*wc)];     % transformation matrix b/w rx estimator and tracking system measurement vector

% Initialise target state vector x(k) for k=1
range = 500;             % range (m)
vel   = 5.555;          % radial velocity (m/s)
accl  = 0.5;            % radial acceleration (m/s^2)  (0.05g)
x1 = [range vel accl]'; % initial target state
q = 0.005;               % 0.005 m/s^2 process (acceleration) noise
Q  = q*dt;               % 0.01 process noise


% define x0 based on the initial conditions
x0 = F\x1;      % x0 = inv(F)*x1



% Define the initial state process noise
P0 = diag(10*[10 0.5 0.05]);
x = linear_model(x1, P0, F, Q, G, nk);


xest_con = zeros(3, nk+1);
xest_con(:,1) = x1;
xest_con(:,2) = F*x1;

% Initialise optimal waveform tracker output
xest_opt = xest_con;

%Then calc the diffb/w x1 and x0 to determine what P0 should be

P_con = zeros(3,3,nk+1);
P_con(:,:,1) = P0;
P_con(:,:,2) = F*P0*F'+Q;
P_opt = P_con;

S_con = zeros(2,2,nk+1);
S_opt = S_con;

mse_con = zeros(1,nk);      % Mean-Square Track Error
mse_opt = mse_con;

vgv_con = zeros(1,nk);      % Validation Gate Volume
vgv_opt = vgv_con;

Kalgain_con = zeros(3,2,nk);
Kalgain_opt = Kalgain_con;

opt_wfm = zeros(1, nk+1);
opt_wfm(1:2) = [con_wfm con_wfm];            % select the initial optimal waveforms

k = 1;
% matlabpool(4)

% Compute the Noise Covariance Matrices for the waveform library
rk = x(1,k);
N_theta = N_cov(theta, T, rk, SNR_1000);

% Select the noise covariance matrix for the initial waveform parameters
Nk_con = N_theta(:,:,con_wfm);                       % measurement noise covariance
y_cov(:,k) = H*x(:,k) + chol(Nk_con)'*randn(2,1);    % measurement vector of target state (r, rr)

Nk_opt = N_theta(:,:,opt_wfm(1));                       % measurement noise covariance
y_opt(:,k) = H*x(:,k) + chol(Nk_opt)'*randn(2,1);    % measurement vector of target state (r, rr)

mse_con(k) = trace(P_con(:,:,k));
mse_opt(k) = trace(P_opt(:,:,k));


for k = 2:nk
    
    % Calculate the measurement
    rk = x(1,k);     % target range (m)
    N_theta = N_cov(theta, T, rk, SNR_1000);
    Nk_act = N_theta(:,:,con_wfm);
    y_cov(:,k) = H*x(:,k) + chol(Nk_act)'*randn(2,1);   % measurement vector of target state (r, rr)
    yk = y_cov(:,k);    
    Pk = P_con(:,:,k);
    %%%%%%%%%%%%% Perform conventional tracking  %%%%%%%%%%%%%%%%%%
    % Select measurement noise covariance for convectional track

    % Recalculated measurement noise covariance for estimated range
    
    xestk = xest_con(:,k);
    Nk_con = N_cov(theta(con_wfm,:), T, xestk(1), SNR_1000); % measurement noise covariance
    % Perform Kalman filter update
    
    [xest_updated, P_updated, Kk, Sk] = kf_corr (xestk, yk, Nk_con, Pk, H);
    xest_con(:,k) = xest_updated;
    P_con(:,:,k) = P_updated;
    S_con(:,:,k) = Sk;
    
    % Perform Kalman filter prediction
    [xest_pred, Ppred] = kf_pred (xest_updated, P_updated, F, Q, G);
    
    xest_con(:,k+1) = xest_pred;
    P_con(:,:,k+1) = Ppred;
    
    Kalgain_con(:,:,k) = Kk;
    Pk = P_updated;
    mse_con(k) = trace(Pk);
    vgv_con(k) = det(Sk);
    
    %%%%%%%%%%%%% Perform optimal waveform selection tracking  %%%%%%%%%%%%%%%%%%
    % Calculated measurement noise
    Nk_act = N_theta(:,:,opt_wfm(k));
    y_opt(:,k) = H*x(:,k) + chol(Nk_act)'*randn(2,1);   % measurement vector of target state (r, rr)
    yk = y_opt(:,k);
    Pk = P_opt(:,:,k);
    
    % Perform Kalman filter update
    
    % Recalculated measurement noise covariance for estimated range
    
    xestk = xest_opt(:,k);

    Nk_opt = N_cov(theta(opt_wfm(k),:), T, xestk(1), SNR_1000); % measurement noise covariance
    
    
    [xest_updated, P_updated, Kk, Sk] = kf_corr (xestk, yk, Nk_opt, Pk, H);
    xest_opt(:,k) = xest_updated;
    P_opt(:,:,k) = P_updated;
    
    % Perform Kalman filter prediction
    [xest_pred, Ppred] = kf_pred (xest_updated, P_updated, F, Q, G);
    
    xest_opt(:,k+1) = xest_pred;
    P_opt(:,:,k+1) = Ppred;
    S_opt(:,:,k) = Sk;
    
    Kalgain_opt(:,:,k) = Kk;
    Pk = P_updated;
    mse_opt(k) = trace(Pk);
    vgv_opt(k) = det(Sk);
    
    % Compute predicted measurement noise covariance matrix
    xk_1 = F*xest_opt(:,k);
    rk_1 = xk_1(1);         % next predicted range
   
    % Compute predicted track error covariance matrix Pk+1|k+1
    N_theta_pred = N_cov(theta, T, rk_1, SNR_1000);
        
    if (opsw == 1) ||( opsw == 2)
        
        parfor n = 1:N
            % Calculate the Predicted Measurement Noise for the entire waveform library
            [~, P_theta(:,:,n), ~, S_theta(:,:,n)] = kf_corr (0, 0, N_theta_pred(:,:,n), Ppred, H);
            
            switch opsw
                case 1 % Minimise Mean Square Track Error
                    mse_theta(n) = trace(P_theta(:,:,n));
                case 2 % Minimise Validation Gate Volume
                    vgv_theta(n) = det(S_theta(:,:,n));
            end
        end
        
        switch opsw
            case 1 % Minimise Mean Square Track Error
                [~, idx] = min(mse_theta(opt));
            case 2 % Minimise Validation Gate Volume
                [~, idx] = min(vgv_theta(opt));
        end
    % Select the next waveform    
    opt_wfm(k+1) = opt(idx);        
        
    end
    
    
    if opsw == 3 % Minimise Validation Gate Volume (Closed Form Solution)
        % Calculate the Optimal waveform parameters for each waveform type
        P = Ppred;
        % Triangular
        l_tri = (30*P(1,1) / (wc^2*P(2,2)))^(1/4);
        % Gaussian
        l_gau = (P(1,1) / (wc^2*P(2,2)))^(1/4);
        % Gaussian with LFM
        b_gau_lfm = (-wc*P(1,2)) / (2*P(1,1));
        l_gau_lfm = (P(1, 1)^2/(wc^2*(P(1, 1)*P(2, 2)-P(1, 2)^2)))^(1/4);
        % Generate a temporary waveform library 
        theta_vgv = [
            1, l_tri    ,0;
            2, l_gau    , 0;
            3, l_gau_lfm, b_gau_lfm];
        % Calculate the Predicted Measurement Noise for each of waveform types
        N_theta_pred = N_cov(theta_vgv, T, rk_1, SNR_1000);
        % Calculate the VGV for each waveform
        parfor n = 1:3
            [~, P_theta(:,:,n), ~, S_theta(:,:,n)] = kf_corr (0, 0, N_theta_pred(:,:,n), Ppred, H);
            vgv_theta(n) = det(S_theta(:,:,n));
        end
        % Minimise Validation Gate Volume
        [~, idx] = min(vgv_theta);
        
        % Save the Optimal Waveform into the Waveform Library
        theta(k,:) = theta_vgv(idx,:);
        % Select the next waveform
        opt_wfm(k+1) = opt(k);
        
    end
    
end

% matlabpool close

opt_wfm = opt_wfm(1:nk);

t = 0:dt:dt*(nk-1);  % time vector;
x = x(:,1:nk);
xest_con = xest_con(:,1:nk);


figure(1)
range_target = x(1, 1:nk);
range_meas = y_cov(1,1:nk);
range_est = xest_con(1, 1:nk);
plot(t, range_target); hold on
plot(t, range_meas, 'r'); hold on
plot(t, range_est, 'k'); hold on
plot(t, xest_opt(1, 1:nk), 'g'); hold on
legend('Target Range', 'Measured Range', 'Con Est', 'Opt Est')

figure(2)
vel_target = x(2,1:nk);
vel_meas = y_cov(2,1:nk);
vel_est = xest_con(2, 1:nk);
plot(t, vel_target); hold on
plot(t, vel_meas, 'r'); hold on
plot(t, vel_est, 'k'); hold on
plot(t, xest_opt(2, 1:nk), 'g'); hold on
legend('Target Velocity', 'Measured Velocity', 'Conventional Est', 'Optimal Est')

figure(3)
switch opsw
    case 1
        y1 = mse_con';
        y2 = mse_opt';
        ylab = 'Total Mean Square Error';
        plot([t', t'], [y1, y2], 'Marker', '.');
        case {2,3}
        y1 = vgv_con';
        y2 = vgv_opt';
        ylab = 'Validation Gate Volume';
        semilogy([t', t'], [y1, y2], 'Marker', '.');
end

xlabel('Elasped Track Time (sec)')
ylabel(ylab)
legend('Conventional', 'Optimal')

pulse_length = theta(opt_wfm, 2);
figure(4)
plot(t, pulse_length', 'Marker', '.')
xlabel('Elasped Track Time (sec)')
ylabel('Effective Pulse Length  \lambda (sec)')

wfm_class = theta(opt_wfm, 1);
figure(5)
plot(t, wfm_class', 'Marker', '.')
xlabel('Elasped Track Time (sec)')
ylabel('Waveform Class')

fm_sweep = theta(opt_wfm, 3)/2/pi/1000;
figure(6)
plot(t, fm_sweep', 'Marker', '.')
xlabel('Elasped Track Time (sec)')
ylabel('FM Sweep Rate (kHz/sec)')


pl_scale = zeros(nk, 1);
for k = 1:3;
    
    idx = find(wfm_class == k);
    
    switch k
        case 1
            pl_scale(idx) = 2;
            
        case {2,3}
            pl_scale(idx) = 7.4338;
    end
end

total_pulse_length = theta(opt_wfm, 2).* pl_scale;
figure(7)
plot(t, total_pulse_length', 'Marker', '.')
xlabel('Elasped Track Time (sec)')
ylabel('Total Pulse Length (sec)')
    