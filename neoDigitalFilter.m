clear; close all; clc;

%read audio data
[audioIn, fs] = audioread('input_audio.wav');

%filter specs
fp = 22000;
rp = 1;
rs = 17;

%analog elliptic prototype
Wp = 2*pi*fp/fs;                                        %normalize fp
[z, p, k] = ellip(20, rp, rs, Wp, 'low', 's');          %generate an array of the poles and zeroes

%compute from zeros & poles
np = numel(p);                                          %size of array of poles
r = zeros(np, 1);                                       %what
for i = 1:np                                            %for loop over np
    others = p([1:i-1, i+1:np]);                        %what
    %display(others);
    r(i) = k*prod(p(i)-z)/prod(p(i)-others);            %what
end


used = false(np, 1); %what
B = []; A = [];                                         %crate empty matrices for B and A
for i = 1:np                                            %for loop over np
    if used(i), continue; end                           %what
    if abs(imag(p(i))) < 1e-12 * max(1, abs(p(i)))      %check if pole is real
        B = [B; real(r(i)) 0 0];                        %what
        A = [A; 1 -exp(real(p(i))) 0];                  %WHAT
        used(i) = true;                                 %what
    else
        j = find(~used & abs(p - conj(p(i))) < 1e-9*abs(p(i)), 1);
        used([i j]) = true;
        ez = exp(p(i));
        B = [B; 2*real(r(i)), -2*real(r(i)*conj(ez)), 0];
        A = [A; 1, -2*real(ez), abs(ez)^2];
    end
end

nsec = size(B,1);

%parallel sum of biquads + feedthrough
audioOut = k * audioIn;
for s = 1:nsec
    audioOut = audioOut + filter(B(s, :), A(s, :), audioIn,  [], 1);
end


%check for audio cliping
peakVal = max(abs(audioOut(:)));
if peakVal > 1, audioOut = audioOut / peakVal; end


%write audio
audiowrite('output_audio.wav', audioOut, fs);

%plots check later
t = (0:size(audioIn,1)-1) / fs;

figure('Color', 'white', 'Position', [100, 50, 1100, 850]);

subplot(3, 1, 1);
plot(t, audioIn(:,1), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hold on;
plot(t, audioOut(:,1), 'Color', [0 0.4 0.8], 'LineWidth', 1.0);
grid on;
title(sprintf('Time-Domain Waveforms (Noisy Input vs. Order 20 Elliptic, Impulse Invariance @ %g kHz)', fp/1000));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Noisy input', 'Filtered', 'Location', 'northeast');
xlim([0, t(end)]);
hold off;

windowLen = 1024; overlapLen = 512; nfft = 2048;

subplot(3, 1, 2);
spectrogram(audioIn(:,1), windowLen, overlapLen, nfft, fs, 'yaxis');
title('Noisy Audio Spectrogram');
ylim([0, fs/2/1000]); ylabel('Freq (kHz)');
colorbar; clim([-120, 0]);

subplot(3, 1, 3);
spectrogram(audioOut(:,1), windowLen, overlapLen, nfft, fs, 'yaxis');
title(sprintf('Filtered Audio Spectrogram (Cutoff: %g kHz)', fp/1000));
ylim([0, fs/2/1000]); ylabel('Freq (kHz)');
colorbar; clim([-120, 0]);
hold on; yline(fp/1000, '--r', 'LineWidth', 1.2); hold off;

%% Frequency response of the designed digital filter
nf = 8192;
w  = linspace(0, pi, nf);
H  = k * ones(1, nf);
for s = 1:nsec
    H = H + freqz(B(s,:), A(s,:), w);
end
figure('Color','white');
plot(w/pi*fs/2/1000, 20*log10(abs(H)+eps), 'LineWidth', 1.2); grid on;
xline(fp/1000, '--r');
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('Impulse-Invariant Elliptic Filter Response');
ylim([-20 5]);

        