%% MVDR BEAMFORMING EXAMPLE - COMPLETE CORRECTED CODE
% For a uniform linear array (ULA) with M microphones.
% Simulates a target source at 60° and an interferer at 120°.
% Compares MVDR output with Delay‑Sum beamformer.

clear; close all; clc;

%% Fix random seed for reproducibility
rng(42);
fprintf('Random seed fixed for reproducibility\n');

%% -------------------- PARAMETERS ----------------------------------------
fprintf('\n=== STEP 1: PARAMETER INITIALIZATION ===\n');

c = 343;                    % speed of sound (m/s)
fs = 16000;                 % sampling frequency (Hz)
duration = 3;                % signal duration (seconds)
nSamples = fs * duration;    % total number of samples

fprintf('fs = %d Hz\n', fs);
fprintf('duration = %.2f s\n', duration);
fprintf('nSamples = %d\n', nSamples);

% Array geometry
M = 6;                      % number of microphones
d = 0.05;                   % spacing (m)   (5 cm)
mic_pos = (0:M-1)' * d;      % positions along x-axis

% Source directions (azimuth angles, in degrees)
target_angle = 60;           % desired direction
interf_angle = 120;          % interference direction

% Source distances (for attenuation)
target_dist = 2.0;
interf_dist = 2.5;

fprintf('Array: %d microphones, spacing = %.2f cm\n', M, d*100);
fprintf('Target direction: %d°, Interferer direction: %d°\n', target_angle, interf_angle);

%% -------------------- GENERATE TEST SIGNALS ----------------------------
fprintf('\n=== STEP 2: SIGNAL GENERATION ===\n');

t = (0:nSamples-1)' / fs;
fprintf('Time vector created: length(t) = %d\n', length(t));
fprintf('Time range: %.3f s to %.3f s\n', t(1), t(end));

% PROTECT t FROM BEING OVERWRITTEN - Create a safe copy
t_safe = t;
fprintf('Created t_safe backup with length %d\n', length(t_safe));

% Target: speech-like signal (harmonic complex)
f0 = 150;
target = sin(2*pi*f0*t) + 0.5*sin(2*pi*2*f0*t) + 0.3*sin(2*pi*3*f0*t);
target = target .* (0.8 + 0.2*sin(2*pi*4*t));   % amplitude modulation
target = target / max(abs(target));
fprintf('Target signal created: length = %d\n', length(target));

% Interference: bandpass filtered noise - MAKE IT STRONGER
interf = randn(nSamples, 1);
[b, a] = butter(4, [500 3000]/(fs/2), 'bandpass');
interf = filter(b, a, interf);
interf = interf / max(abs(interf)) * 1.5;  % Increased from 0.7 to 1.5
fprintf('Interference scaled to 1.50\n');
fprintf('Interference signal created: length = %d\n', length(interf));

%% -------------------- VERIFY SIGNAL CREATION ---------------------------
fprintf('\n=== VERIFYING SIGNAL CREATION ===\n');
fprintf('t created: length = %d\n', length(t));
fprintf('t_safe created: length = %d\n', length(t_safe));
fprintf('target created: length = %d\n', length(target));
fprintf('Sample values - t(1:5):\n');
disp(t(1:5));
fprintf('Sample values - target(1:5):\n');
disp(target(1:5));

%% -------------------- PROPAGATION DELAYS & ATTENUATION ----------------
fprintf('\n=== STEP 3: DELAY CALCULATION ===\n');

% Compute time delays from each source to each microphone.
% Far‑field assumption: plane wave, delays depend only on angle.
target_delays = (mic_pos * cosd(target_angle)) / c;
interf_delays = (mic_pos * cosd(interf_angle)) / c;

% Convert to samples
target_delay_samps = round(target_delays * fs);
interf_delay_samps = round(interf_delays * fs);

% Attenuation (1/r)
target_atten = 1 / target_dist;
interf_atten = 1 / interf_dist;

fprintf('Target delays (samples): '); fprintf('%d ', target_delay_samps); fprintf('\n');
fprintf('Interferer delays (samples): '); fprintf('%d ', interf_delay_samps); fprintf('\n');

%% -------------------- GENERATE MICROPHONE SIGNALS ----------------------
fprintf('\n=== STEP 4: MICROPHONE SIGNAL GENERATION ===\n');

received = zeros(nSamples, M);

for m = 1:M
    % Target contribution
    dT = target_delay_samps(m);
    if dT > 0
        sigT = [zeros(dT,1); target(1:end-dT)];
    elseif dT < 0
        sigT = target(-dT+1:end);
        sigT = [sigT; zeros(-dT,1)];
    else
        sigT = target;
    end
    received(:,m) = received(:,m) + target_atten * sigT;
    
    % Interference contribution
    dI = interf_delay_samps(m);
    if dI > 0
        sigI = [zeros(dI,1); interf(1:end-dI)];
    elseif dI < 0
        sigI = interf(-dI+1:end);
        sigI = [sigI; zeros(-dI,1)];
    else
        sigI = interf;
    end
    received(:,m) = received(:,m) + interf_atten * sigI;
end

% Add a small amount of sensor noise
received = received + 0.005 * randn(size(received));
fprintf('Microphone signals created: size(received) = [%d, %d]\n', size(received,1), size(received,2));

%% -------------------- STEERING VECTOR FUNCTION -------------------------
function d = steering_vector(mic_pos, angle, freq, c)
    % CORRECTED steering vector
    % mic_pos: column vector of microphone positions (M x 1)
    % angle:   direction in degrees
    % freq:    frequency in Hz
    % returns: steering vector (M x 1)
    d = exp(1j * 2 * pi * freq * (mic_pos * cosd(angle)) / c);
end

fprintf('\n=== STEP 5: STEERING VECTOR FUNCTION DEFINED ===\n');

%% -------------------- CORRECTED MVDR IMPLEMENTATION -------------------
fprintf('\n=== STEP 6: CORRECTED MVDR BEAMFORMING ===\n');

% STFT parameters
win_len = 512;               % window length (samples)
hop = 256;                   % hop size
win = hann(win_len);

% Compute STFT of each microphone signal
nFreq = win_len/2 + 1;
nFrames = floor((nSamples - win_len)/hop) + 1;
fprintf('STFT parameters: %d frequency bins, %d frames\n', nFreq, nFrames);

X = zeros(nFreq, nFrames, M);
for m = 1:M
    X(:,:,m) = spectrogram(received(:,m), win, win_len-hop, win_len, fs);
end
fprintf('STFT computed for all %d channels\n', M);

% Pre-allocate output STFT
Y_mvdr = zeros(nFreq, nFrames);
Y_ds = zeros(nFreq, nFrames);

% Regularization factor
reg = 0.1;  % Strong regularization
fprintf('Regularization factor: %.4f\n', reg);

% Process each frequency bin
for f = 1:nFreq
    freq = (f-1) * fs / win_len;
    
    % Skip extreme frequencies
    if freq < 100 || freq > 7000
        continue;
    end
    
    % Get data for this frequency (all frames, all microphones)
    X_f = squeeze(X(f,:,:));  % nFrames x M
    
    % Steering vector for target direction
    d_target = steering_vector(mic_pos, target_angle, freq, c);
    
    % ----- CORRECTED DELAY-SUM WEIGHTS -----
    % For a signal arriving from angle theta, the steering vector is exp(1j*phi)
    % To align signals, we need to apply conjugate weights
    w_ds = conj(d_target) / M;
    
    % ----- MVDR WEIGHTS -----
    % Estimate covariance matrix
    R = (X_f' * X_f) / nFrames;
    
    % Add regularization (diagonal loading) - adaptive
    R_reg = R + reg * trace(R)/M * eye(M);
    
    % Compute MVDR weights with numerical stability
    try
        % Use Cholesky decomposition for better numerical stability
        [L, p] = chol(R_reg, 'lower');
        if p == 0
            % R_reg is positive definite, use Cholesky
            R_inv_d = L' \ (L \ d_target);
        else
            % Fall back to pseudoinverse
            R_inv_d = pinv(R_reg) * d_target;
        end
        w_mvdr = R_inv_d / (d_target' * R_inv_d);
        
        % Check for numerical issues
        if any(isnan(w_mvdr)) || any(isinf(w_mvdr)) || norm(w_mvdr) > 100
            error('Numerical issues detected');
        end
    catch
        % If MVDR fails, use delay-sum weights
        w_mvdr = w_ds;
    end
    
    % Apply weights to each frame - CORRECTED: X_f * conj(w)
    Y_ds(f,:) = X_f * conj(w_ds);
    Y_mvdr(f,:) = X_f * conj(w_mvdr);
end

fprintf('Frequency-domain processing complete\n');

%% -------------------- INVERSE STFT FUNCTION ----------------------------
function x = istft_mvdr(S, win, hop, nSamples)
    % Corrected inverse STFT with proper scaling
    [nFreq, nFrames] = size(S);
    win_len = length(win);
    x = zeros(nSamples, 1);
    
    for frame = 1:nFrames
        % Reconstruct full spectrum (positive and negative frequencies)
        full_fft = [S(:,frame); conj(flipud(S(2:end-1,frame)))];
        
        % IFFT to get time domain
        frame_signal = real(ifft(full_fft)) * (win_len / hop);  % Scaling factor
        
        % Apply synthesis window
        frame_signal = frame_signal .* win;
        
        % Overlap-add
        idx = (frame-1)*hop + (1:win_len);
        if max(idx) <= nSamples
            x(idx) = x(idx) + frame_signal;
        end
    end
    
    % Normalize
    x = x / (win_len / hop);
end

% Inverse STFT to obtain time‑domain signals
fprintf('\n=== STEP 7: TIME-DOMAIN RECONSTRUCTION ===\n');
beamformed_ds = istft_mvdr(Y_ds, win, hop, nSamples);
beamformed_mvdr = istft_mvdr(Y_mvdr, win, hop, nSamples);
fprintf('Beamformed signals created: length(beamformed_ds) = %d\n', length(beamformed_ds));
fprintf('Beamformed signals created: length(beamformed_mvdr) = %d\n', length(beamformed_mvdr));

%% -------------------- VERIFY FREQUENCY-DOMAIN DS ----------------------
fprintf('\n=== VERIFYING FREQUENCY-DOMAIN DS ===\n');

% Create correct time-domain delay-sum for verification
correct_ds = zeros(nSamples, 1);
min_delay = min(target_delay_samps);
for m = 1:M
    shift = target_delay_samps(m) - min_delay;
    if shift > 0
        % Correct alignment: delay OTHER channels, not this one
        correct_ds(shift+1:end) = correct_ds(shift+1:end) + received(1:end-shift, m);
    else
        correct_ds = correct_ds + received(:, m);
    end
end
correct_ds = correct_ds / M;

% Compare with frequency-domain DS
corr_val = corrcoef(correct_ds(1000:9000), beamformed_ds(1000:9000));
fprintf('Correlation between correct DS and freq-domain DS: %.3f\n', corr_val(1,2));

if corr_val(1,2) < 0.95
    fprintf('ERROR: Frequency-domain DS still wrong! Check implementation.\n');
else
    fprintf('✓ Frequency-domain DS is correct!\n');
end

%% -------------------- DIAGNOSTICS - CHECK BEAMFORMER OUTPUT -----------
fprintf('\n=== DIAGNOSTICS - CHECKING BEAMFORMER OUTPUT ===\n');

% Signal statistics
fprintf('\nSignal Statistics:\n');
fprintf('  Target: mean=%.3f, std=%.3f, range=[%.2f, %.2f]\n', ...
        mean(target), std(target), min(target), max(target));
fprintf('  Mic 1: mean=%.3f, std=%.3f, range=[%.2f, %.2f]\n', ...
        mean(received(:,1)), std(received(:,1)), min(received(:,1)), max(received(:,1)));
fprintf('  DS out: mean=%.3f, std=%.3f, range=[%.2f, %.2f]\n', ...
        mean(beamformed_ds), std(beamformed_ds), min(beamformed_ds), max(beamformed_ds));
fprintf('  MVDR out: mean=%.3f, std=%.3f, range=[%.2f, %.2f]\n', ...
        mean(beamformed_mvdr), std(beamformed_mvdr), min(beamformed_mvdr), max(beamformed_mvdr));

% Correlation with target
[c_ds, l_ds] = xcorr(beamformed_ds, target);
[c_mvdr, l_mvdr] = xcorr(beamformed_mvdr, target);
[peak_ds, idx_ds] = max(abs(c_ds));
[peak_mvdr, idx_mvdr] = max(abs(c_mvdr));
fprintf('\nCorrelation with target:\n');
fprintf('  DS: peak=%.3f at lag=%d samples\n', peak_ds, l_ds(idx_ds));
fprintf('  MVDR: peak=%.3f at lag=%d samples\n', peak_mvdr, l_mvdr(idx_mvdr));

%% -------------------- PROTECT TIME VECTOR CHECK -----------------------
fprintf('\n=== PROTECT TIME VECTOR CHECK ===\n');
if exist('t_safe', 'var')
    fprintf('t_safe exists with length %d\n', length(t_safe));
else
    error('t_safe does not exist!');
end

%% -------------------- PLOT RESULTS -------------------------------------
fprintf('\n=== STEP 8: PLOTTING RESULTS ===\n');

% Check variables
fprintf('Variable check:\n');
fprintf('  length(t_safe): %d\n', length(t_safe));
fprintf('  length(target): %d\n', length(target));
fprintf('  size(received): [%d, %d]\n', size(received,1), size(received,2));
fprintf('  length(beamformed_ds): %d\n', length(beamformed_ds));
fprintf('  length(beamformed_mvdr): %d\n', length(beamformed_mvdr));

% Determine safe plotting length
plot_len = min(2000, length(t_safe));
fprintf('Plotting %d samples (%.3f seconds)\n', plot_len, plot_len/fs);

t_plot = t_safe(1:plot_len);

figure('Position', [100 100 1200 800]);

% Subplot 1: Clean target vs noisy microphone
subplot(3,1,1);
plot(t_plot, target(1:plot_len), 'b', 'LineWidth', 1.5);
hold on;
plot(t_plot, received(1:plot_len,1), 'r', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Amplitude');
title('Clean Target vs Noisy Microphone 1');
legend('Target', 'Mic 1');
grid on;
xlim([t_plot(1), t_plot(end)]);
ylim([-1.5 1.5]);

% Subplot 2: Beamformed outputs comparison
subplot(3,1,2);
plot(t_plot, beamformed_ds(1:plot_len), 'g', 'LineWidth', 1.2);
hold on;
plot(t_plot, beamformed_mvdr(1:plot_len), 'm', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Amplitude');
title('Beamformed Outputs Comparison');
legend('Delay‑Sum', 'MVDR');
grid on;
xlim([t_plot(1), t_plot(end)]);
ylim([-1.5 1.5]);

% Subplot 3: MVDR spectrogram
subplot(3,1,3);
spectrogram(beamformed_mvdr, hann(512), 256, 512, fs, 'yaxis');
title('MVDR Output Spectrogram');
caxis([-100 -20]);
colorbar;

sgtitle('MVDR Beamforming Results');
fprintf('Plotting complete\n');

%% -------------------- BEAMPATTERN COMPARISON ---------------------------
fprintf('\n=== STEP 9: BEAMPATTERN ANALYSIS ===\n');

figure('Position', [100 100 1200 400]);
freqs_to_plot = [500, 2000, 4000];
angles = 0:359;

for idx = 1:3
    freq = freqs_to_plot(idx);
    d_target = steering_vector(mic_pos, target_angle, freq, c);
    
    % Delay‑Sum beampattern
    w_ds = conj(d_target) / M;
    pat_ds = zeros(size(angles));
    for a = 1:length(angles)
        d_angle = steering_vector(mic_pos, angles(a), freq, c);
        pat_ds(a) = abs(w_ds' * d_angle);
    end
    
    % MVDR beampattern (using theoretical covariance with interferer)
    d_interf = steering_vector(mic_pos, interf_angle, freq, c);
    % Assume interference power = 1, noise power = 0.1
    R_theory = d_interf * d_interf' + 0.1 * eye(M);
    w_mvdr = (R_theory \ d_target) / (d_target' * (R_theory \ d_target));
    pat_mvdr = zeros(size(angles));
    for a = 1:length(angles)
        d_angle = steering_vector(mic_pos, angles(a), freq, c);
        pat_mvdr(a) = abs(w_mvdr' * d_angle);
    end
    
    subplot(1,3,idx);
    polarplot(deg2rad(angles), pat_ds, 'b-', 'LineWidth', 2);
    hold on;
    polarplot(deg2rad(angles), pat_mvdr, 'r--', 'LineWidth', 2);
    title(sprintf('%d Hz', freq));
    legend('DS', 'MVDR', 'Location', 'southoutside');
    rlim([0 1]);
    grid on;
end
sgtitle('Beampattern Comparison: Delay-Sum (blue) vs MVDR (red)');

%% -------------------- PERFORMANCE METRICS ------------------------------
fprintf('\n=== STEP 10: PERFORMANCE METRICS ===\n');

% Compute input SNR (using first microphone as reference)
% First, align the target with the microphone signal
[c_target, l_target] = xcorr(received(:,1), target);
[~, idx_target] = max(abs(c_target));
target_delay = l_target(idx_target);

if target_delay > 0
    target_aligned_input = [zeros(target_delay,1); target(1:end-target_delay)];
elseif target_delay < 0
    target_aligned_input = target(-target_delay+1:end);
    target_aligned_input = [target_aligned_input; zeros(-target_delay,1)];
else
    target_aligned_input = target;
end
len_input = min(length(target_aligned_input), length(received(:,1)));
target_aligned_input = target_aligned_input(1:len_input);
mic1_trim = received(1:len_input,1);

noise_in_mic1 = mic1_trim - target_aligned_input;
input_SNR = 10*log10(mean(target_aligned_input.^2) / mean(noise_in_mic1.^2));
fprintf('Input SNR (Mic 1): %.2f dB\n', input_SNR);

% Align outputs with target using cross-correlation
fprintf('\nAligning beamformed outputs with target...\n');

% Delay-Sum alignment
[~, idx_ds] = max(abs(c_ds));
delay_ds = l_ds(idx_ds);
fprintf('DS optimal lag: %d samples\n', delay_ds);

if delay_ds > 0
    target_aligned_ds = [zeros(delay_ds,1); target(1:end-delay_ds)];
elseif delay_ds < 0
    target_aligned_ds = target(-delay_ds+1:end);
    target_aligned_ds = [target_aligned_ds; zeros(-delay_ds,1)];
else
    target_aligned_ds = target;
end
len_ds = min(length(target_aligned_ds), length(beamformed_ds));
target_aligned_ds = target_aligned_ds(1:len_ds);
bf_ds_trim = beamformed_ds(1:len_ds);

% MVDR alignment
[~, idx_mvdr] = max(abs(c_mvdr));
delay_mvdr = l_mvdr(idx_mvdr);
fprintf('MVDR optimal lag: %d samples\n', delay_mvdr);

if delay_mvdr > 0
    target_aligned_mvdr = [zeros(delay_mvdr,1); target(1:end-delay_mvdr)];
elseif delay_mvdr < 0
    target_aligned_mvdr = target(-delay_mvdr+1:end);
    target_aligned_mvdr = [target_aligned_mvdr; zeros(-delay_mvdr,1)];
else
    target_aligned_mvdr = target;
end
len_mvdr = min(length(target_aligned_mvdr), length(beamformed_mvdr));
target_aligned_mvdr = target_aligned_mvdr(1:len_mvdr);
bf_mvdr_trim = beamformed_mvdr(1:len_mvdr);

% Compute residuals and output SNR
residual_ds = bf_ds_trim - target_aligned_ds;
residual_mvdr = bf_mvdr_trim - target_aligned_mvdr;

output_SNR_ds = 10*log10(mean(target_aligned_ds.^2) / mean(residual_ds.^2));
output_SNR_mvdr = 10*log10(mean(target_aligned_mvdr.^2) / mean(residual_mvdr.^2));

fprintf('\n===== FINAL PERFORMANCE RESULTS =====\n');
fprintf('Delay‑Sum Output SNR: %.2f dB (gain: %.2f dB)\n', ...
        output_SNR_ds, output_SNR_ds - input_SNR);
fprintf('MVDR Output SNR:      %.2f dB (gain: %.2f dB)\n', ...
        output_SNR_mvdr, output_SNR_mvdr - input_SNR);

% Plot residual comparison
figure('Position', [100 100 800 400]);
plot(t_plot, residual_ds(1:plot_len), 'r', 'LineWidth', 0.8);
hold on;
plot(t_plot, residual_mvdr(1:plot_len), 'b', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Amplitude');
title('Residual Noise Comparison (after alignment)');
legend('Delay-Sum Residual', 'MVDR Residual');
xlim([t_plot(1), t_plot(end)]);
grid on;

%% -------------------- FINAL DIAGNOSTIC --------------------------------
fprintf('\n=== FINAL DIAGNOSTIC ===\n');

% Check if beamformers are preserving the target
target_energy = sum(target.^2);
ds_energy = sum(beamformed_ds.^2);
mvdr_energy = sum(beamformed_mvdr.^2);

fprintf('Target energy: %.2f\n', target_energy);
fprintf('DS energy: %.2f (ratio: %.2f)\n', ds_energy, ds_energy/target_energy);
fprintf('MVDR energy: %.2f (ratio: %.2f)\n', mvdr_energy, mvdr_energy/target_energy);

if ds_energy/target_energy < 0.5
    fprintf('WARNING: DS is attenuating the target too much!\n');
end
if mvdr_energy/target_energy < 0.3
    fprintf('WARNING: MVDR is attenuating the target too much!\n');
end

%% -------------------- FINAL VARIABLE CHECK ----------------------------
fprintf('\n=== FINAL VARIABLE CHECK ===\n');
fprintf('At end of script:\n');
fprintf('  exist(''t_safe''): %d\n', exist('t_safe', 'var'));
fprintf('  length(t_safe): %d\n', length(t_safe));
fprintf('  t_safe(1): %.6f, t_safe(end): %.6f\n', t_safe(1), t_safe(end));
fprintf('Script completed at: %s\n', datestr(now));
fprintf('\n=== SCRIPT COMPLETED SUCCESSFULLY ===\n');