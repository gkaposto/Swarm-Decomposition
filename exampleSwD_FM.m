%
% This script realizes the FM decomposition using SwD that presented in the paper: 
% "Swarm decomposition: A novel signal analysis using swarm intelligence"
% by G.Apostolidis and L.Hadjileontiadis, published in Elsevier Signal Processing, Volume 132, March 2017, Pages 40–50
%
% ps: It may take 1-2 minutes to execute.
%% ******************* FM signal generation *******************
Q     = 299;
Fs    = 20000;
A     = [3 2 1];
F     = [700 1000 1300];
theta = [0 pi 0];
beta  = [0.9/Q 0.5/Q 0.1/Q];

[x ,x_c] = FMGenerator(Q, A, F, beta, theta, Fs);

n = 0:1:Q-1;
figure; plot(n, x);

no_comp = size(x_c, 1);
nfft    = 2^nextpow2(length(x));
f_axis  = (Fs / 2) * linspace(-1, 1, nfft);

for i = 1:no_comp
    
    X_c = fftshift(abs(fft(x_c(i, :), nfft)));
    
    figure; 
    subplot(2, 1, 1); plot(n, x_c(i, :));
    subplot(2, 1, 2); plot(f_axis, X_c);
end

%% ******************** SwD execution **************************
decimation_factor = 2;
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

welch_window      = round(length(x) / 4);
welch_no_overlap  = round(welch_window / 2);
welch_nfft        = 2^nextpow2(length(x_dec));
param_struct      = struct('P_th', 0.05, ...
                           'StD_th', 0.01, ...
                           'Welch_window', welch_window, ...
                           'Welch_no_overlap', welch_no_overlap, ...
                           'Welch_nfft', welch_nfft);

% SwD execution
y_res = SwD(x, param_struct);

% interpolation to reverse the decimation
y = resample(y_res, p, 1);
y = y(1:length(x), :);
