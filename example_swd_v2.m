%clear; clc; close all hidden;
%% Input signal
% Load the signal 
[s, fs] = bat;       % obtained from https://ltfat.github.io/doc/signals/bat.html

% Normalize the signal
x       = (s - mean(s)) / std(s);
[Px, w] = pwelch(x);

%% SWD
cmps_thresh  = 0.025;
detail_depth = 1e-5;
cmps         = SwD_v2(x, cmps_thresh, detail_depth);

%% Reconstruction
reconstructed_x = sum(cmps, 2);
[Px_recon, ~]   = pwelch(reconstructed_x);

%% Plot
figure;
subplot(2, 1, 1); title('Time-Domain'); xlabel('samples');
hold on; box on; grid on;
plot(x); plot(reconstructed_x, '--'); legend('x', 'x recon');
subplot(2, 1, 2); title('Frequency Domain'); xlabel('frequency (rad/{\pi})');
hold on; box on; grid on;
plot(w/pi, Px); plot(w/pi, Px_recon, '--'); legend('x', 'x recon');

disp(['Number of extracted components = ', num2str(size(cmps, 2))]);

energy_cmps        = sum(cmps.^2, 1);
energy_total       = sum(energy_cmps);
plot_energy_thresh = 0.000;
idx                = find((energy_cmps / energy_total) <= plot_energy_thresh);
idx2               = find((energy_cmps / energy_total) > plot_energy_thresh);
if (isempty(idx))
    cmps_new = cmps(:, idx2);
else
    cmps_new = [cmps(:, idx2), sum(cmps(:, idx), 2)];
end
nocmps = size(cmps_new, 2);

figure;
for i = 1:nocmps-1
    
    subplot(nocmps-1, 1, i); plot(cmps_new(:, i));
    
end
