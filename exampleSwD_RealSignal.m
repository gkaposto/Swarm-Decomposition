%
% This script realizes the real signal decomposition using SwD that presented in the paper: 
% "Swarm decomposition: A novel signal analysis using swarm intelligence"
% by G.Apostolidis and L.Hadjileontiadis, published in Elsevier Signal Processing, Volume 132, March 2017, Pages 40–50
%
% The real signal is an EEG signal found in https://physionet.org/physiobank/database/
%
%% ******************* Load and pre-process the singal *******************
load 'chb01_01_edfm.mat';

signal_idx = 20;
x    = val(signal_idx, :);
x    = (x - mean(x)) / max(abs(x));
fs   = 256;
dt   = 1 / fs;
L    = length(x);
t    = 0:dt:(L - 1) * dt;
nfft = 2 ^ nextpow2(L);
X    = abs(fft(x, nfft)) / L;
f    = (fs / 2) * linspace(0, 1, nfft/2);

figure;
subplot(1, 2, 1); plot(t, x); 
subplot(1, 2, 2); plot(f, X(1:end/2));

% bandpass filtering
[b, a] = butter(5, [7 31.5] / (fs / 2), 'bandpass');
x_filt = filtfilt(b, a, x);
nfft   = 2 ^ nextpow2(length(x_filt));
X_filt = abs(fft(x_filt, nfft)) / nfft;
f      = (fs / 2) * linspace(0, 1, nfft/2);

figure; 
subplot(1, 2, 1); plot(t, x_filt); 
subplot(1, 2, 2); plot(f, X_filt(1:end/2));

% downsamping
q     = 4; 
fs2   = fs / q;
dt2   = 1 / fs2;
x_dec = resample(x_filt, 1, q); 
L_dec = length(x_dec);
nfft  = 2 ^ nextpow2(L_dec);
X_dec = abs(fft(x_dec, nfft)) / L_dec; 
t2    = 0:dt2:(length(x_dec) - 1) * dt2;
f2    = (fs2 / 2) * linspace(0, 1, nfft/2);

figure; 
subplot(1, 2, 1); plot(t2, x_dec); 
subplot(1, 2, 2); plot(f2, X_dec(1:end/2));
%% ******************* SWD execution *******************
L1    = length(x_dec);
s_SWD = [x_dec, zeros(1,100)];  % zero-padding
L2    = length(s_SWD);

% parameter setting
% ------------------------
% How coarse/fine the decomposition will be depends on the smoothing of the FFT spectrum using Savitsky-Golay filter, according to the paper.
% Here a better approach is adopted calculating the Welch spectrum and we configure the smoothness 
% via the Welch spectrum parameters which are the Welch window length and the overlapping (https://en.wikipedia.org/wiki/Welch's_method).
% If the Welch window is equal to the length of the input signal and the overlap is 0% the Welch spectrum is identical to FFT-based power spectrum.
% More segments (i.e. "small" windows) means more smoothing so more coarse decomposition, whereas, 
% less segments (i.e. "big" windows) means less smoothing so finer decomposition. 
% In short, the smoothing process is replaced from Savitsky-Golay filtering of FFT to calculating Welch spectrum. 
% ------------------------
welch_window      = round(L2 / 16);
welch_no_overlap  = round(welch_window / 2);
welch_nfft        = 2^nextpow2(L2);
param_struct      = struct('P_th', 0.2, ...
                           'StD_th', 0.2, ...
                           'Welch_window', welch_window, ...
                           'Welch_no_overlap', welch_no_overlap, ...
                           'Welch_nfft', welch_nfft);

y_SWD_res = SwD(s_SWD, param_struct); disp('end')
y_SWD = y_SWD_res.';

%% ******************* Plot results *******************
no_comp_SWD = size(y_SWD, 1);

figure; 
for i = 1:1:no_comp_SWD
    subplot(no_comp_SWD, 1, i); plot(t2, y_SWD(i, 1:L1));   
end

figure; 
for i = 1:1:no_comp_SWD
    nfft = 2 ^ nextpow2(length(y_SWD(i, :)));
    Y    = abs(fft(y_SWD(i,:), nfft)) / nfft;
    ff   = (fs2 / 2) * linspace(0,1, nfft / 2);
    
    subplot(no_comp_SWD, 1, i); plot(ff, Y(1:end/2));
end
