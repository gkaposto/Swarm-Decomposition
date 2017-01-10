%
% This script realizes the AM decomposition using SwD that presented in the paper: 
% "Swarm decomposition: A novel signal analysis using swarm intelligence"
% by G.Apostolidis and L.Hadjileontiadis, published in Elsevier Signal Processing, Volume 132, March 2017, Pages 40–50
%
%% ******************* AM signal generation *******************
Q     = 339;
A     = [2 1 0.9 3 2.5];
a     = [0.0042 0.0037 0.0029 0.0037 0.0033];
F     = [100 140 190 230 270];
theta = [0 pi pi/2 0 pi/3];
Fs    = 6400;

[x, x_c] = AMGenerator(Q, A, a, F, theta, Fs);
L = size(x_c, 1);
n = 0:1:Q-1;

figure; 
subplot(L+1, 1, 1);plot(n, x, 'b');
for i = 1:1:L
    subplot(L+1, 1, i+1); plot(n, x_c(i, :), 'b');
end

figure; 
for i = 1:1:L
    X  = fftshift(abs(fft(x_c(i, :))));
    ff = (Fs / 2) * linspace(-1, 1, length(X));
    
    subplot(L, 1, i); plot(ff, X, 'b');
end

%% ******************** SwD execution **************************

% decimation to avoid the existance of very small oscillatory modes.
decimation_factor = 4; 
zero_padding      = 100;
x_dec             = [resample(x, 1, decimation_factor), zeros(1, zero_padding)];

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

welch_window      = round(length(x_dec) / 2);
welch_no_overlap  = round(welch_window / 2);
welch_nfft        = 2^nextpow2(length(x_dec));
param_struct      = struct('P_th', 0.05, ...
                           'StD_th', 0.05, ...
                           'Welch_window', welch_window, ...
                           'Welch_no_overlap', welch_no_overlap, ...
                           'Welch_nfft', welch_nfft);

% SwD execution
y_res = SwD(x_dec, param_struct);

% interpolation to reverse the decimation
y = resample(y_res, p, 1);
y = y(1:length(x), :);
