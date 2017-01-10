%
% This script realizes the synthetic non-stationary signal decomposition using SwD that presented in the paper: 
% "Swarm decomposition: A novel signal analysis using swarm intelligence"
% by G.Apostolidis and L.Hadjileontiadis, published in Elsevier Signal Processing, Volume 132, March 2017, Pages 40–50
%
%% ******************* Generate synthetic non-stationary singal *******************
Q     = 499;
coord = [250, 0.2*pi, 300, 0.7;...
         125, 0.6*pi, 125, 1.5; ...
         375, 0.7*pi, 100, 1.9];
SNR   = 10;

[x, x_c, sig] = EvaluationSignalGenerator(Q, coord, SNR);

x   = x.';
x_c = x_c.';
n   = 0:1:Q - 1;
L   = size(x_c, 1);

% Plotting
figure; plot(n, x, 'b');

figure; 
subplot(L+2, 1, 1); plot(n, x, 'b');
for i = 1:1:L
    subplot(L+2, 1, i+1); plot(n, x_c(i, :), 'b');
end
subplot(L+2, 1, L+2); plot(n, x-sig.', 'b');

figure; 
for i = 1:1:L
    X  = fftshift(abs(fft(x_c(i,:)))) / Q;
    ff = linspace(-1, 1, length(X));
    
    subplot(L, 1, i); plot(ff, X, 'b');
end
%% ******************* SWD execution *******************

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
welch_window      = round(Q / 8);
welch_no_overlap  = round(welch_window / 2);
welch_nfft        = 2^nextpow2(Q);
param_struct = struct('P_th', 0.05, ...
                      'StD_th', 0.05, ...
                      'Welch_window', welch_window, ...
                      'Welch_no_overlap', welch_no_overlap, ...
                      'Welch_nfft', welch_nfft);

y_SWD_res = SwD(x, param_struct); disp('end')

% Plotting
no_components = size(y_SWD_res, 2);
figure; 
for i = 1:1:no_components
    subplot(no_components, 1, i); plot(y_SWD_res(:, i));
end

figure; 
for i = 1:1:no_components
    Y = fftshift(abs(fft(y_SWD_res(:, i)))) / Q;
    ff = linspace(-1, 1, length(Y)); 
    
    subplot(no_components, 1, i); plot(ff, Y);
end