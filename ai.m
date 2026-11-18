% =========================================================================
% Elliptic Filter Audio Processing with Waveforms & Spectrograms
% Target: 48 kHz Mono WAV | Cutoff Frequency: 22 kHz
% =========================================================================

clear; close all; clc;

% 1. File Configuration
inputFileName = 'untitled.wav';   % Replace with your input file name
outputFileName = 'filtered_audio.wav';

% 2. Read Audio File
[audioIn, fs] = audioread(inputFileName);

% Verify sampling rate constraint
if fs ~= 48000
    warning('Sampling rate is %d Hz, expected 48000 Hz.', fs);
end

% Ensure mono audio
if size(audioIn, 2) > 1
    audioIn = mean(audioIn, 2);
    disp('Converted stereo audio to mono.');
end

% 3. Filter Specifications (Updated for 22 kHz Cutoff)
fp = 22000;               % Passband edge frequency (22 kHz)
fs_stop = 23500;          % Stopband edge frequency (23.5 kHz)
rp = 0.5;                 % Passband ripple in dB
rs = 60;                  % Stopband attenuation in dB

% Normalize frequencies relative to Nyquist frequency (fs / 2 = 24000 Hz)
Wp = fp / (fs / 2);
Ws = fs_stop / (fs / 2);

% 4. Design Elliptic Filter
[n, Wn] = ellipord(Wp, Ws, rp, rs);
fprintf('Designed Elliptic Filter Order: %d\n', n);

[b, a] = ellip(n, rp, rs, Wn, 'low');

% 5. Apply Zero-Phase Filter
audioOut = filtfilt(b, a, audioIn);

% Normalize output to prevent clipping
audioOut = audioOut / max(abs(audioOut)) * 0.95;

% 6. Save Filtered Audio
audiowrite(outputFileName, audioOut, fs);
fprintf('Filtered audio saved to: %s\n', outputFileName);

% 7. Time Vector for Waveforms
t = (0:length(audioIn)-1) / fs;

% 8. Comprehensive Visualization (Waveforms + Spectrograms)
figure('Color', 'white', 'Position', [100, 50, 1100, 850]);

% --- Subplot 1: Time-Domain Waveforms ---
subplot(3, 1, 1);
plot(t, audioIn, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hold on;
plot(t, audioOut, 'Color', [0 0.4 0.8], 'LineWidth', 1.0);
grid on;
title('Time-Domain Waveforms (Original vs. Elliptic Filtered @ 22 kHz)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original', 'Filtered', 'Location', 'northeast');
xlim([0, t(end)]); 
hold off;

% Spectrogram Parameters
windowLen = 1024;
overlapLen = 512;
nfft = 2048;

% --- Subplot 2: Original Audio Spectrogram ---
subplot(3, 1, 2);
spectrogram(audioIn, windowLen, overlapLen, nfft, fs, 'yaxis');
title('Original Audio Spectrogram (48 kHz Mono)');
ylim([0, fs/2/1000]);
ylabel('Freq (kHz)');
colorbar;
clim([-120, 0]);

% --- Subplot 3: Filtered Audio Spectrogram ---
subplot(3, 1, 3);
spectrogram(audioOut, windowLen, overlapLen, nfft, fs, 'yaxis');
title('Filtered Audio Spectrogram (Cutoff: 22 kHz)');
ylim([0, fs/2/1000]);
ylabel('Freq (kHz)');
colorbar;
clim([-120, 0]);

% Highlight passband and stopband edges on the filtered spectrogram
hold on;
yline(fp/1000, '--r', 'LineWidth', 1.2);
yline(fs_stop/1000, '--w', 'LineWidth', 1.2);
hold off;