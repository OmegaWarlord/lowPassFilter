clear; close all; clc;

% read the audio
inputFileName = 'input_audio.wav';
noisyAudioFileName = 'noisy_audio.wav';
outputFileName = 'filtered_audio.wav';
[audioIn, fs] = audioread(inputFileName);

noiseLevel = 0.1;  % adjust as needed
hfNoise = randn(size(audioIn));
% uses a butter filter just to get the high frequency noise
[bNoise, aNoise] = butter(4, 22000 / (fs / 2), 'high');
hfNoise = filter(bNoise, aNoise, hfNoise);

audioInNoisy = audioIn + noiseLevel * hfNoise;
audiowrite(noisyAudioFileName, audioInNoisy, fs);

% project specifications
fp = 22000; 
fs_stop = 23500;
rp = 1;
rs = 17;
order = 20;


% normalized passband for ellip
Wp = fp / (fs / 2);
Ws = fs_stop / (fs / 2);

[z, p, k] = ellip(order, rp, rs, Wp, 'low');
[sos, g] = zp2sos(z,p,k);

audioOut = g * sosfilt(sos, audioInNoisy);

peak = max(abs(audioOut), [], 'all');
if peak > 0
    audioOut = 0.95 * audioOut / peak;
end

audiowrite(outputFileName, audioOut, fs);
fprintf('Filtered audio saved to: %s\n', outputFileName);

t = (0:length(audioInNoisy)-1) / fs;

figure('Color', 'white', 'Position', [100, 50, 1100, 850]);

subplot(3, 1, 1);
plot(t, audioInNoisy, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hold on;
plot(t, audioOut, 'Color', [0 0.4 0.8], 'LineWidth', 1.0);
grid on;
title('Time-Domain Waveforms (Noisy Input vs. Elliptic Filtered @ 22 kHz)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noisy input', 'Filtered', 'Location', 'northeast');
xlim([0, t(end)]);
hold off;

windowLen = 1024;
overlapLen = 512;
nfft = 2048;

subplot(3, 1, 2);
spectrogram(audioInNoisy, windowLen, overlapLen, nfft, fs, 'yaxis');
title('Noisy Audio Spectrogram');
ylim([0, fs/2/1000]);
ylabel('Freq (kHz)');
colorbar;
clim([-120, 0]);

subplot(3, 1, 3);
spectrogram(audioOut, windowLen, overlapLen, nfft, fs, 'yaxis');
title('Filtered Audio Spectrogram (Cutoff: 22 kHz)');
ylim([0, fs/2/1000]);
ylabel('Freq (kHz)');
colorbar;
clim([-120, 0]);

hold on;
yline(fp/1000, '--r', 'LineWidth', 1.2);
hold off;
