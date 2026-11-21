clear; close all; clc;

%garrett was here

% read the audio
inputFileName = 'input_audio.wav';
outputFileName = 'filtered_audio.wav';
[audioIn, fs] = audioread(inputFileName);

% project specifications
fp = 22000; 
rp = 1;
rs = 17;
order = 4;


% normalized passband for ellip
%Wp = fp / (fs / 2);

Wp = 2 * pi * fp;

[z, p, k] = ellip(order, rp, rs, Wp, 'low', 's');

[b, a] = zp2tf(z, p, k);

[bz,az] = impinvar(b,a,fs, 1e-4);

[sos, g] = tf2sos(bz, az);

audioOut = g * sosfilt(sos, audioIn, 1);

display(audioOut);

peakVal = max(abs(audioOut(:)));
if peakVal > 1
    audioOut = audioOut / peakVal;
end

audiowrite(outputFileName, audioOut, fs);
fprintf('Filtered audio saved to: %s\n', outputFileName);


t = (0:length(audioIn)-1) / fs;

figure('Color', 'white', 'Position', [100, 50, 1100, 850]);

subplot(3, 1, 1);
plot(t, audioIn, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
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
spectrogram(audioIn, windowLen, overlapLen, nfft, fs, 'yaxis');
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
