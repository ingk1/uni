function J = fisher(pulse_class, lambda, chirp_rate, SNR)

% This function computes the (2x2) Fisher Information matrix J(THETA) as a
% function of the pulse waveform parameters. It assumes a zero mean time
% and frequency (t and w). Pulse amplitudesare normalised such that the
% energy s(t)s*(t) = 1
%
% INPUTS:
% pulse_class - 1: AM Triangular Pulse
%               2: AM Gaussian Pulse
%               3: AM Gaussian with LFM chirp
%
% lambda - the Gaussian waveform pulse length parameter (s)
% The effective pulse length is 7.4338*lambda
%
% chirp_rate - The LFM sweep rate (rad/s^2)
% SNR - The signal to noise ratio of the RX signal (dB)
% This is normally eta = 2*Er/No




% Assume that the mean frequency and time are both zero (the pulses are
% centred about the origin)
w = 0;
t = 0;

l = lambda;
b = chirp_rate;

% Compute the mean w^2, mean w^2 and wt for each waveform type
switch pulse_class
    
    case 0  % No pulse
        w_2 = 0;
        t_2 = 0;
        wt = 0; 
    
    case 1  % AM Triangular Pulse
        
        w_2 = 3/l^2;
        t_2 = l^2/10;
        wt = 0; 
        
    case 2  % AM Gaussian Pulse
        
        w_2 = 1/(2*l^2);
        t_2 = l^2/2;
        wt = 0;        
        
    case 3      % AM Gaussian Pulse with LFM chirp
        
        w_2 = (1/(2*l^2)) + 2*b^2*l^2;
        t_2 = l^2/2;
        wt = b*l^2;        
        
    otherwise
        
        error('Pulse case not valid')
        
end


% convert SNR from DB to a linear ratio
eta = 10^(SNR/10);      % SNR (ratio)

% Define symmetric matrix
% note that w denotes w bar, the mean freqency and t denotes t bar, the
% mean time.

u11 = w_2 - w^2;        %(w^2) - (w)^2 where w is w bar
u12 = wt - w*t;         % (wt) - w*t
u21 = wt - w*t;         % (wt) - w*t
u22 = t_2 - t^2;         % (t^2) - (t)^2

U = [u11 u12; u21 u22];
J = eta*U;              % calc Fisher information Matrix



