function [theta] = waveform_lib(waveform_class, upw, sweep_rate, debug)
%
% This function returns a waveform parameter vector; theta. Theta is a Nx3
% matrix that contains the waveform parameters for each pulse type:
% 
% Column 1: Waveform class (1 = Tri, 2 = Gaussian, 3 = Gaussian with LFM
% chirp)
% Column 2: Gaussian waveform pulse length parameter, lambda (s)
% Note, The effective pulse length is 7.4338*lambda for a Gaussian pulse
% Column 3: The LFM sweep rate (Hz/s)
%
% 
% Inputs:
% waveform class   - This is just a numerical value that is assigned to  
%                    allow the pulse type to be uniquely catergorised
% upw  (s)         - Vector of length 3 that contains the minimum and
%                    maximum effective pulse length of the waveform and the number
%                    of linearly spaced values
% sweep_rate (Hz/s)- Vector of length 3 that contains the minimum and
%                    maximum LFM sweep rate of the waveform and the number
%                    of linearly spaced values


% Waveform Parameters
% waveform_class = [1 2 3 1 1];
% upw = [0.01 0.3 5;
%     0.01 0.03 5;
%     0.01 0.06 5;
%     11.83E-3/2 11.83E-3/2 1;
%     12E-3/7.4338 12E-3/7.4338 1];
% sweep_rate = [0 0 1;
%     0 0 1;
%     -15E3*2*pi 0 5;
%     0 0 1;
%     0 0 1];

% debug = 0;


num_waveforms = upw(:,3).*sweep_rate(:,3);

K = length(waveform_class);             % number of different types of waveform classes
N_waveforms = sum(num_waveforms);       % total number of waveform in the library

theta = zeros(N_waveforms, 3);          % the waveform parameter library
idx = 1;

for k = 1:K
    
    Nk = num_waveforms(k);
    
    % Generate the pulselength values
    lambdak = linspace(upw(k,1), upw(k,2), upw(k,3));
    LAMBDAk = repmat(lambdak,sweep_rate(k,3),1);
    LAMBDAk = reshape(LAMBDAk,1,Nk)';
    
    % Generate the LFM sweeprate values
    bk = linspace(sweep_rate(k,1), sweep_rate(k,2), sweep_rate(k,3));
    Bk = repmat(bk,1,upw(k,3))';
    
    % Compile the Waveform Library
    theta(idx:idx+num_waveforms(k)-1, 1) = waveform_class(k)*ones(1,Nk);
    theta(idx:idx+num_waveforms(k)-1, 2) = LAMBDAk;
    theta(idx:idx+num_waveforms(k)-1, 3) = Bk;
    
    idx = idx+num_waveforms(k);
end

 
% Setup debug flags if enabled
if (~isempty(debug)) && (debug ~= 0);
    k_debug = debug;    % sets the waveform to plot
    debug = 1;      % allows plotting of pulse envelope in time and frequency domains
    
end


% theta = [];
% %% Triangular pulse waveform
% waveform_class = 1;
% upw = [0.01 0.3]; % triangular pulse min and max length (s)
% sweep_rate = [0 0];
% N = 5; % number of different pulse types for Triangular waveform class
% theta = [theta; waveform_param(waveform_class, upw, sweep_rate, N)];
% %% Gaussian Pulse Waveform
% waveform_class = 2;
% upw = [0.01 0.3]; % triangular pulse min and max length (s)
% sweep_rate = [0 0];
% N = 5; % number of different pulse types for Triangular waveform class
% theta = [theta; waveform_param(waveform_class, upw, sweep_rate, N)];
% %% Gaussian LFM Pulse Waveform
% waveform_class = 3;
% upw = [0.01 0.3]; % triangular pulse min and max length (s)
% sweep_rate = [-15E3 0];
% N = 5; % number of different pulse types for Triangular waveform class
% theta = [theta; waveform_param(waveform_class, upw, sweep_rate, N)];



if debug
    
    k = k_debug;    
    waveform_class = theta(k,1);
    lambda = theta(k,2);
    chirp_rate = theta(k,3);
    epl = 7.4338*lambda;  % effective pulse length is = 7.4338*lambda
    
    fs = 1E3;      % 1KHz sampling rate
    ts = 1/fs;
    t = -lambda:ts:lambda;
    N = length(t);
    f = [0:N-1]*fs/(N-1);
    w = 2*pi*f;
    

    
    switch waveform_class
        
        case 1      % Class 1:Triangular Waveform
            t = -lambda:ts:lambda;
            st = sqrt(2/lambda)*(1-abs(t)/lambda);       % pulse signal (time domain)
            Sw = sqrt(2/lambda)*(sinc(w*lambda/2)).^2;
            
        case 2 % Class 2: Gaussian AM waveform
            
            t = -5*lambda:ts:5*lambda;
            st = (1/(pi*lambda^2))^(1/4)*exp(-t.^2/(2*lambda^2));
            Sw = (4*pi*lambda^2)^(1/4)*exp(-w.^2*lambda^2/2);
            
            
        case 3  % Class 3: Gaussian LFM waveform_class - need to put the right parameters in here
            t = -5*lambda:ts:5*lambda;
            st = (1/(pi*lambda^2))^(1/4)*exp(-t.^2/(2*lambda^2));
            Sw = (4*pi*lambda^2)^(1/4)*exp(-w.^2*lambda^2/2);
    end
    
    figure (55)
    subplot(2,1,1), plot(t, st);        % time domain plot
    subplot(2,1,2), plot(w, Sw),        % frequency domain plot
    xlim([0,fs/2]);
end

