%% Step 1: Initialize parameters
clear; close all; clc;

% Physical constants
c = 343;                % speed of sound (m/s) at 20°C
fs = 16000;             % sampling frequency (Hz)
duration = 5;           % seconds
nSamples = fs * duration;

% Room dimensions (from thesis: 1.2 x 1.5 x 5 m)
room = [1.2, 1.5, 5];

% Microphone array: 4-element linear array along y-axis, 5 cm spacing
M = 4;                  % number of microphones
d = 0.05;               % spacing (m)
array_center = [0.6, 0.725, 0.5];   % center of array (x, y, z)

% Compute individual microphone positions
mic_pos = zeros(M, 3);
for m = 1:M
    % Linear along y: equally spaced around center
    mic_pos(m, :) = [array_center(1), ...
                     array_center(2) - (M-1)*d/2 + (m-1)*d, ...
                     array_center(3)];
end

% Source positions (from thesis Table 3.1 / Figure 3.1)
target_pos = [0.125, 0.725, 0.5];   % target speech
noise_pos  = [0.6, 1.2, 0.5];       % interference noise (example)

% Distances from array center
target_dist = norm(target_pos - array_center);
noise_dist  = norm(noise_pos - array_center);

fprintf('Target distance: %.2f m\n', target_dist);
fprintf('Noise distance:  %.2f m\n', noise_dist);