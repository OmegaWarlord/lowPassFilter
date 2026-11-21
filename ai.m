clear; close all; clc;

% read the audio
inputFileName  = 'input_audio.wav';
outputFileName = 'filtered_audio.wav';
[audioIn, fs] = audioread(inputFileName);

% project specifications
fp    = 22000;   % passband edge (Hz)
rp    = 1;       % passband ripple (dB)
rs    = 17;      % stopband attenuation (dB)
order = 20;

if fp >= fs/2
    error('fp (%g Hz) must be below Nyquist (%g Hz).', fp, fs/2);
end

%% 1) Analog elliptic prototype (normalized time: T = 1 sample, rad/sample)
% Working in rad/sample keeps all numbers near unity. Building huge
% polynomials with Wp = 2*pi*22000 and order 20 is numerically useless.
Wp = 2*pi*fp/fs;
[z, p, k] = ellip(order, rp, rs, Wp, 'low', 's');

%% 2) Impulse invariance computed directly from zeros/poles
% H(s) = k*prod(s-z)/prod(s-p) = d + sum_i r_i/(s-p_i)
% Residues from zpk (no polynomial expansion, so no ill-conditioning):
np = numel(p);
r  = zeros(np,1);
for i = 1:np
    others = p([1:i-1, i+1:np]);
    r(i) = k * prod(p(i) - z) / prod(p(i) - others);
end
% Equal numerator/denominator degree (even-order elliptic) gives a
% direct feedthrough term d = k (leading coefficient ratio).
d = k;

% Map each pole: r/(s-p) -> r/(1 - e^p z^-1)   (T = 1, h[n] = h_c(n))
% Combine complex-conjugate pairs into real second-order sections.
used = false(np,1);
B = []; A = [];
for i = 1:np
    if used(i), continue; end
    if abs(imag(p(i))) < 1e-12 * max(1,abs(p(i)))    % real pole
        B = [B; real(r(i)) 0 0];
        A = [A; 1 -exp(real(p(i))) 0];
        used(i) = true;
    else
        j = find(~used & abs(p - conj(p(i))) < 1e-9*abs(p(i)), 1);
        used([i j]) = true;
        ez = exp(p(i));
        B = [B; 2*real(r(i)), -2*real(r(i)*conj(ez)), 0];
        A = [A; 1, -2*real(ez), abs(ez)^2];
    end
end
nsec = size(B,1);

% Stability check (impulse invariance maps LHP poles inside unit circle)
fprintf('Order %d, %d parallel sections, max |pole| = %.6f\n', ...
        order, nsec, max(abs(exp(p))));

%% 3) Filter: parallel sum of biquads + feedthrough
audioOut = d * audioIn;
for s = 1:nsec
    audioOut = audioOut + filter(B(s,:), A(s,:), audioIn, [], 1);
end

peakVal = max(abs(audioOut(:)));
if peakVal > 1
    audioOut = audioOut / peakVal;
end

audiowrite(outputFileName, audioOut, fs);
fprintf('Filtered audio saved to: %s\n', outputFileName);

%% Plots
t = (0:size(audioIn,1)-1) / fs;

figure('Color', 'white', 'Position', [100, 50, 1100, 850]);

subplot(3, 1, 1);
plot(t, audioIn(:,1), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hold on;
plot(t, audioOut(:,1), 'Color', [0 0.4 0.8], 'LineWidth', 1.0);
grid on;
title(sprintf('Time-Domain Waveforms (Noisy Input vs. Order-%d Elliptic, Impulse Invariance @ %g kHz)', order, fp/1000));
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
H  = d * ones(1, nf);
for s = 1:nsec
    H = H + freqz(B(s,:), A(s,:), w);
end
figure('Color','white');
plot(w/pi*fs/2/1000, 20*log10(abs(H)+eps), 'LineWidth', 1.2); grid on;
xline(fp/1000, '--r');
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('Impulse-Invariant Elliptic Filter Response');
ylim([-100 5]);