

clear; close all; clc;

%% -------------------- STEP 1: PARAMETERS --------------------------------
c = 343;                    % speed of sound (m/s)
fs = 16000;                 % sampling frequency (Hz)
duration = 3;                % signal duration (seconds)
nSamples = fs * duration;

% Array geometry: 8-element linear array, 4 cm spacing
M_full = 8;                  % total microphones
d = 0.04;                    % spacing (m)

% Subarray configuration
subarray_size = 4;           % microphones per subarray
subarray_overlap = 2;        % overlapping microphones
n_subarrays = M_full - subarray_size + 1;  % number of subarrays

% Calculate microphone positions
array_center = [0, 0, 0];    % center at origin
mic_pos = zeros(M_full, 3);
for m = 1:M_full
    mic_pos(m,:) = [0, (m - (M_full+1)/2) * d, 0];
end

% Display array configuration
fprintf('=== Array Configuration ===\n');
fprintf('Full array: %d microphones, spacing = %.2f cm\n', M_full, d*100);
fprintf('Subarrays: %d subarrays of %d mics with %d overlap\n', ...
        n_subarrays, subarray_size, subarray_overlap);
fprintf('Total degrees of freedom: %d\n', n_subarrays * subarray_size);

%% -------------------- STEP 2: CREATE SUBARRAYS -------------------------
% Define subarray indices
subarray_indices = cell(n_subarrays, 1);
for s = 1:n_subarrays
    start_idx = s;
    end_idx = s + subarray_size - 1;
    subarray_indices{s} = start_idx:end_idx;
end

% Display subarray composition
fprintf('\n=== Subarray Composition ===\n');
for s = 1:n_subarrays
    fprintf('Subarray %d: Mics ', s);
    fprintf('%d ', subarray_indices{s});
    fprintf('\n');
end

%% -------------------- STEP 3: DESIGN FILTERS ---------------------------
% We'll design filters for each subarray using frequency-domain method

% Filter parameters
L = 64;                     % filter length
freq_range = [100, 6000];    % frequency range of interest
N_freq = 512;                % number of frequency points
f = linspace(0, fs/2, N_freq);  % frequency vector

% Target frequency response (desired beampattern)
steering_angle = 90;         % desired look direction (degrees)

% Pre-allocate filter coefficients for each subarray
subarray_filters = zeros(L, n_subarrays);

% Design filters for each subarray
for s = 1:n_subarrays
    % Get microphone positions for this subarray
    sub_mic_pos = mic_pos(subarray_indices{s}, :);
    
    % Design superdirective beamformer for this subarray
    % (This is a simplified version - in practice, use MVDR or similar)
    H_freq = zeros(N_freq, subarray_size);
    
    for fi = 1:N_freq
        freq = f(fi);
        if freq < 100 || freq > 6000
            H_freq(fi,:) = 0;
            continue;
        end
        
        % Steering vector for look direction
        d_steer = steering_vector(sub_mic_pos, steering_angle, freq, c);
        
        % Simple superdirective: pseudo-inverse of steering matrix
        H_freq(fi,:) = d_steer' / (d_steer * d_steer' + 0.01);
    end
    
    % Convert to time-domain FIR filter
    % (This is simplified - use proper filter design)
    for m = 1:subarray_size
        % IFFT to get impulse response
        h = real(ifft([H_freq(:,m); conj(flipud(H_freq(2:end-1,m)))]));
        
        % Store in filter bank (simplified - truncate to length L)
        subarray_filters(:, (s-1)*subarray_size + m) = h(1:L);
    end
end

fprintf('\n=== Filter Design Complete ===\n');
fprintf('Filter length: %d taps\n', L);
fprintf('Number of filters: %d\n', n_subarrays * subarray_size);

%% -------------------- STEP 4: GENERATE TEST SIGNALS --------------------
% Create time vector
t = (0:nSamples-1)' / fs;

% Target signal (speech-like)
f0 = 150;
target = sin(2*pi*f0*t) + 0.5*sin(2*pi*2*f0*t) + 0.3*sin(2*pi*3*f0*t);
target = target .* (0.8 + 0.2*sin(2*pi*4*t));
target = target / max(abs(target));

% Interference from different directions
interf1 = 0.7 * sin(2*pi*200*t + randn*0.1) .* (0.5 + 0.5*sin(2*pi*5*t));
interf2 = 0.5 * randn(nSamples, 1);
[b, a] = butter(4, [500, 3000]/(fs/2), 'bandpass');
interf2 = filter(b, a, interf2);
interf2 = interf2 / max(abs(interf2)) * 0.5;

% Source positions
target_angle = 90;          % desired direction
interf1_angle = 45;         % interference 1 direction
interf2_angle = 135;        % interference 2 direction

% Distances
target_dist = 2.0;
interf1_dist = 2.5;
interf2_dist = 3.0;

fprintf('\n=== Source Configuration ===\n');
fprintf('Target: angle = %d°, distance = %.1f m\n', target_angle, target_dist);
fprintf('Interf1: angle = %d°, distance = %.1f m\n', interf1_angle, interf1_dist);
fprintf('Interf2: angle = %d°, distance = %.1f m\n', interf2_angle, interf2_dist);

%% -------------------- STEP 5: GENERATE RECEIVED SIGNALS ----------------
% Calculate delays for each source to each microphone
target_delays = zeros(M_full, 1);
interf1_delays = zeros(M_full, 1);
interf2_delays = zeros(M_full, 1);

for m = 1:M_full
    % Target
    src_pos = [target_dist * cosd(target_angle), ...
               target_dist * sind(target_angle), 0];
    target_delays(m) = norm(src_pos - mic_pos(m,:)) / c;
    
    % Interference 1
    src_pos = [interf1_dist * cosd(interf1_angle), ...
               interf1_dist * sind(interf1_angle), 0];
    interf1_delays(m) = norm(src_pos - mic_pos(m,:)) / c;
    
    % Interference 2
    src_pos = [interf2_dist * cosd(interf2_angle), ...
               interf2_dist * sind(interf2_angle), 0];
    interf2_delays(m) = norm(src_pos - mic_pos(m,:)) / c;
end

% Convert to samples and generate received signals
target_delay_samps = round(target_delays * fs);
interf1_delay_samps = round(interf1_delays * fs);
interf2_delay_samps = round(interf2_delays * fs);

% Attenuation (1/r)
target_atten = 1 ./ target_dist;
interf1_atten = 1 ./ interf1_dist;
interf2_atten = 1 ./ interf2_dist;

% Generate multi-channel received signals
received = zeros(nSamples, M_full);

for m = 1:M_full
    % Add target
    if target_delay_samps(m) > 0
        sig = [zeros(target_delay_samps(m),1); ...
               target(1:end-target_delay_samps(m))];
    else
        sig = target;
    end
    received(:,m) = received(:,m) + target_atten * sig;
    
    % Add interference 1
    if interf1_delay_samps(m) > 0
        sig = [zeros(interf1_delay_samps(m),1); ...
               interf1(1:end-interf1_delay_samps(m))];
    else
        sig = interf1;
    end
    received(:,m) = received(:,m) + 0.7 * interf1_atten * sig;
    
    % Add interference 2
    if interf2_delay_samps(m) > 0
        sig = [zeros(interf2_delay_samps(m),1); ...
               interf2(1:end-interf2_delay_samps(m))];
    else
        sig = interf2;
    end
    received(:,m) = received(:,m) + 0.5 * interf2_atten * sig;
end

% Add sensor noise
received = received + 0.005 * randn(size(received));

%% -------------------- STEP 6: APPLY FILTER-AND-SUM BEAMFORMING ---------
% Method 1: Direct filter-and-sum on full array
fprintf('\n=== Applying Beamforming ===\n');

% We'll use a simple frequency-domain implementation
beamformed_fs = zeros(nSamples, 1);

% STFT parameters
win_len = 512;
hop = 256;
win = hann(win_len);

% STFT of received signals
X = zeros(win_len/2+1, floor((nSamples-win_len)/hop)+1, M_full);
for m = 1:M_full
    X(:,:,m) = spectrogram(received(:,m), win, win_len-hop, win_len, fs);
end

% Process each frequency bin
Y = zeros(size(X,1), size(X,2));
for fi = 1:size(X,1)
    freq = (fi-1) * fs / win_len;
    
    if freq >= 100 && freq <= 6000
        % Steering vector for look direction
        d = steering_vector(mic_pos, target_angle, freq, c);
        
        % Simple superdirective weights
        w = d' / (d * d' + 0.01);
        
        % Apply weights
        for ti = 1:size(X,2)
            x_vec = squeeze(X(fi,ti,:));
            Y(fi,ti) = w * x_vec;
        end
    else
        Y(fi,:) = 0;
    end
end

% ISTFT to get time-domain signal
beamformed_fs = istft(Y, win, hop, nSamples);

%% -------------------- STEP 7: SUBARRAY PROCESSING ---------------------
% Method 2: Subarray beamforming
subarray_outputs = zeros(nSamples, n_subarrays);

for s = 1:n_subarrays
    % Get signals for this subarray
    sub_mics = subarray_indices{s};
    sub_received = received(:, sub_mics);
    
    % Apply subarray beamforming (simpler delay-sum for subarray)
    sub_delays = target_delay_samps(sub_mics);
    min_delay = min(sub_delays);
    align_delays = sub_delays - min_delay;
    
    aligned = zeros(nSamples, subarray_size);
    for m = 1:subarray_size
        if align_delays(m) > 0
            aligned(1:end-align_delays(m), m) = ...
                sub_received(1+align_delays(m):end, m);
        else
            aligned(:, m) = sub_received(:, m);
        end
    end
    
    subarray_outputs(:, s) = sum(aligned, 2) / subarray_size;
end

% Combine subarray outputs (simple averaging)
beamformed_subarray = mean(subarray_outputs, 2);

%% -------------------- STEP 8: RESULTS VISUALIZATION -------------------
figure('Position', [50 50 1400 900]);

% Time-domain comparison (zoomed in)
subplot(3,2,1);
plot(t(1:2000), target(1:2000), 'b', 'LineWidth', 1);
hold on;
plot(t(1:2000), received(1:2000,1), 'r', 'LineWidth', 0.5);
xlabel('Time (s)'); ylabel('Amplitude');
title('Clean Target vs Noisy Microphone 1');
legend('Target', 'Mic 1');
xlim([0 0.1]);
grid on;

subplot(3,2,2);
plot(t(1:2000), beamformed_fs(1:2000), 'g', 'LineWidth', 1);
hold on;
plot(t(1:2000), beamformed_subarray(1:2000), 'm', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amplitude');
title('Beamformed Outputs');
legend('Full FS', 'Subarray');
xlim([0 0.1]);
grid on;

% Spectrograms
subplot(3,2,3);
spectrogram(target, hann(512), 256, 512, fs, 'yaxis');
title('Clean Target Spectrogram');
caxis([-100 -20]);

subplot(3,2,4);
spectrogram(beamformed_fs, hann(512), 256, 512, fs, 'yaxis');
title('Full Filter-Sum Output');
caxis([-100 -20]);

subplot(3,2,5);
spectrogram(beamformed_subarray, hann(512), 256, 512, fs, 'yaxis');
title('Subarray Output');
caxis([-100 -20]);

subplot(3,2,6);
% Residual noise comparison
residual_fs = beamformed_fs - [target; zeros(length(beamformed_fs)-length(target),1)];
residual_sub = beamformed_subarray - [target; zeros(length(beamformed_subarray)-length(target),1)];

plot(t(1:2000), residual_fs(1:2000), 'r', t(1:2000), residual_sub(1:2000), 'b');
xlabel('Time (s)'); ylabel('Amplitude');
title('Residual Noise Comparison');
legend('Full FS', 'Subarray');
xlim([0 0.1]);
grid on;

sgtitle('Filter-and-Sum vs Subarray Beamforming Comparison');

%% -------------------- STEP 9: BEAMPATTERN COMPARISON ------------------
figure('Position', [50 50 1200 400]);

frequencies = [500, 2000, 4000];
angles = 0:359;

for fi = 1:length(frequencies)
    freq = frequencies(fi);
    
    % Full array beampattern
    pattern_full = zeros(size(angles));
    for a = 1:length(angles)
        d = steering_vector(mic_pos, angles(a), freq, c);
        w = steering_vector(mic_pos, target_angle, freq, c)';
        pattern_full(a) = abs(w * d) / M_full;
    end
    
    % Subarray beampattern (average of all subarrays)
    pattern_sub = zeros(size(angles));
    for s = 1:n_subarrays
        sub_mics = subarray_indices{s};
        sub_pos = mic_pos(sub_mics, :);
        for a = 1:length(angles)
            d = steering_vector(sub_pos, angles(a), freq, c);
            w = steering_vector(sub_pos, target_angle, freq, c)';
            pattern_sub(a) = pattern_sub(a) + abs(w * d) / subarray_size;
        end
    end
    pattern_sub = pattern_sub / n_subarrays;
    
    subplot(1,3,fi);
    polarplot(deg2rad(angles), pattern_full, 'b-', 'LineWidth', 2);
    hold on;
    polarplot(deg2rad(angles), pattern_sub, 'r--', 'LineWidth', 2);
    title(sprintf('Beampattern at %d Hz', freq));
    legend('Full', 'Subarray');
    rlim([0 1]);
end
sgtitle('Beampattern Comparison: Full Array vs Subarray Averaging');

%% -------------------- STEP 10: PERFORMANCE METRICS --------------------
% Compute SNR improvement
input_power = mean(target.^2);
noise_power = mean((received(:,1) - target).^2);
input_SNR = 10*log10(input_power / noise_power);

output_power_fs = mean(beamformed_fs.^2);
residual_fs = beamformed_fs - target(1:length(beamformed_fs));
output_SNR_fs = 10*log10(output_power_fs / mean(residual_fs.^2));

output_power_sub = mean(beamformed_subarray.^2);
residual_sub = beamformed_subarray - target(1:length(beamformed_subarray));
output_SNR_sub = 10*log10(output_power_sub / mean(residual_sub.^2));

fprintf('\n=== Performance Metrics ===\n');
fprintf('Input SNR: %.2f dB\n', input_SNR);
fprintf('Full Filter-Sum Output SNR: %.2f dB (improvement: %.2f dB)\n', ...
        output_SNR_fs, output_SNR_fs - input_SNR);
fprintf('Subarray Output SNR: %.2f dB (improvement: %.2f dB)\n', ...
        output_SNR_sub, output_SNR_sub - input_SNR);

%% -------------------- HELPER FUNCTIONS --------------------------------
function d = steering_vector(mic_pos, angle, freq, c)
    % Calculate steering vector for given direction
    M = size(mic_pos, 1);
    array_center = mean(mic_pos, 1);
    
    % Source position in far-field
    source_pos = array_center + 10 * [cosd(angle), sind(angle), 0];
    
    % Calculate delays
    delays = zeros(M, 1);
    for m = 1:M
        delays(m) = norm(source_pos - mic_pos(m,:)) / c;
    end
    min_delay = min(delays);
    delays = delays - min_delay;
    
    % Steering vector
    d = exp(-1j * 2 * pi * freq * delays);
end

function x = istft(S, win, hop, nSamples)
    % Simple inverse STFT
    [nFreq, nFrames] = size(S);
    win_len = length(win);
    
    x = zeros(nSamples, 1);
    for t = 1:nFrames
        % IFFT
        frame = real(ifft([S(:,t); conj(flipud(S(2:end-1,t)))]));
        
        % Overlap-add
        idx = (t-1)*hop + (1:win_len);
        if max(idx) <= nSamples
            x(idx) = x(idx) + frame .* win;
        end
    end
end