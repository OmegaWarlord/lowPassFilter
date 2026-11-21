clear; close all; clc;

%read audio data
[audioIn, fs] = audioread('input_audio.wav');

%filter specs
fp = 22000;
rp = 1;
rs = 17;

%analog elliptic prototype
Wp = 2*pi*fp/fs;                                                            %normalize fp
[zeroes, poles, gain] = ellip(20, rp, rs, Wp, 'low', 's');                  %generate an array of the poles and zeroes

%Transfer function
%H(s) = K · prod(s − z_j) / prod(s − p_k)
np = numel(poles);                                                          %size of array of poles
residue = zeros(np);                                                        %creates an array of scalar 0s of size np
for i = 1:np                                                                %for loop over np
    p = poles(i)-poles([1:i-1, i+1:np]);                                    %calculate the difference pole at i and other poles
    residue(i) = gain*prod(poles(i)-zeroes)/prod(p);                        %calculate the residue of the partial fraction expansion
end


used = false(np);                                                           %creates an array of logical 0s of size np to check if pole has been visited?
B = []; A = [];                                                             %crate empty matrices for coefficients B and A
for i = 1:np                                                                %for loop over np
    if used(i), continue; end                                               %check if pole was already checked
    if abs(imag(poles(i))) < 1e-12 * max(abs(poles(i)))                     %check if pole is real
        B = [B; real(residue(i))];                                          %add coefficient to B array
        A = [A; 1-exp(real(poles(i)))];                                     %add coefficient to A array                                                   %Check pole to used
    else
        ez = exp(poles(i));                                                 %exponential of pole i         
        B = [B; 2*real(residue(i)), -2*real(residue(i)*conj(ez)), 0];       %coefficient of b
        A = [A; 1, -2*real(ez), abs(ez)^2];                                 %coefficient of a
    end
    j = find(~used & abs(poles-conj(poles(i)))<1e-9*abs(poles(i)), 1);      %finds the location of the complex conjugate of the pole
    used([i j]) = true;                                                     %sets used to that to true
end


%parallel sum of biquads + feedthrough
audioOut = gain * audioIn;                                                  %start audioOut
for s = 1:size(B,1)
    audioOut = audioOut + filter(B(s, :), A(s, :), audioIn);                %add filtered value from B & A to audio out
end


%check for audio cliping
peakVal = max(abs(audioOut(:)));
if peakVal > 1, audioOut = audioOut / peakVal; end


%write audio
audiowrite('output_audio.wav', audioOut, fs);

t = (0:size(audioIn,1)-1) / fs;

figure('Color', 'white', 'Position', [100, 50, 1100, 850]);

subplot(3, 1, 1);
plot(t, audioIn(:,1), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); 
hold on;
plot(t, audioOut(:,1), 'Color', [0 0.4 0.8], 'LineWidth', 1.0); 
title(sprintf('Time-Domain Waveforms (Noisy Input vs. Order 20 Elliptic, Impulse Invariance @ %g kHz)', fp/1000));
xlabel('Time (s)'); 
ylabel('Amplitude');
legend('Noisy input', 'Filtered', 'Location', 'northeast');
xlim([0, t(end)]);
hold off;
grid on;

windowLen = 1024; overlapLen = 512; nfft = 2048;

subplot(3, 1, 2);
spectrogram(audioIn(:,1), windowLen, overlapLen, nfft, fs, 'yaxis');
title('Noisy Audio Spectrogram');
ylim([0, fs/2/1000]); 
ylabel('Freq (kHz)');
colorbar; 
clim([-120, 0]);

subplot(3, 1, 3);
spectrogram(audioOut(:,1), windowLen, overlapLen, nfft, fs, 'yaxis');
title(sprintf('Filtered Audio Spectrogram (Cutoff: %g kHz)', fp/1000));
ylim([0, fs/2/1000]); 
ylabel('Freq (kHz)');
colorbar; 
clim([-120, 0]);
hold on; 
yline(fp/1000, '--r', 'LineWidth', 1.2); 
hold off;

%% Frequency response of the designed digital filter
nf = 8192;
w  = linspace(0, pi, nf);
H  = gain * ones(1, nf);
for s = 1:size(B,1)
    H = H + freqz(B(s,:), A(s,:), w);
end
figure('Color','white');
plot(w/pi*fs/2/1000, 20*log10(abs(H)+eps), 'LineWidth', 1.2); 
xline(fp/1000, '--r');
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('Impulse-Invariant Elliptic Filter Response');
ylim([-20 5]);
grid on;

        
%% DFT before / after filtering, with filter response
function plot_dft(audioIn, audioOut, fs, B, A, gain, scale, fp)

if nargin < 7 || isempty(scale)
    scale = 'db';
end
if nargin < 8 || isempty(fp)
    fp = 22000;
end

x = double(audioIn(:,1));
y = double(audioOut(:,1));

segLen = min(numel(x), 65536);
startIdx = floor((numel(x) - segLen) / 2) + 1;
idx = startIdx:(startIdx + segLen - 1);

if segLen > 1
    n = (0:segLen-1).';
    win = 0.5 - 0.5*cos(2*pi*n/(segLen - 1));
else
    win = 1;
end

xSeg = x(idx) .* win;
ySeg = y(idx) .* win;

nfft = 2^nextpow2(segLen);

X = fftshift(fft(xSeg, nfft));
Y = fftshift(fft(ySeg, nfft));
f = (-nfft/2:nfft/2-1).' * fs / nfft;

Xmag = abs(X);
Ymag = abs(Y);

ref = max(Xmag);
if ref == 0
    ref = 1;
end

Xmag = Xmag / ref;
Ymag = Ymag / ref;

w = (0:nfft-1).' * (2*pi/nfft);
H = gain * ones(nfft, 1);
for s = 1:size(B, 1)
    H = H + freqz(B(s,:), A(s,:), w);
end
H = fftshift(H);
Hmag = abs(H);

if strcmpi(scale, 'db')
    Xplot = 20*log10(Xmag + eps);
    Yplot = 20*log10(Ymag + eps);
    Hplot = 20*log10(Hmag + eps);
    yLabelText = 'Magnitude (dB)';
    yLimits = [-120 10];
else
    Xplot = Xmag;
    Yplot = Ymag;
    Hplot = Hmag;
    yLabelText = 'Magnitude';
    yLimits = [0 1.1];
end

figure('Color', 'white', 'Position', [100 100 1000 650]);

subplot(2,1,1);
plot(f/1000, Xplot, 'Color', [0.2 0.45 0.75], 'LineWidth', 1.1);
hold on;
plot(f/1000, Hplot, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
xline(fp/1000, '--r', 'LineWidth', 1.2);
xline(-fp/1000, '--r', 'LineWidth', 1.2);
hold off;
title('Two-Sided DFT of Original Signal');
xlabel('Frequency (kHz)');
ylabel(yLabelText);
xlim([-fs/2 fs/2]/1000);
ylim(yLimits);
grid on;
legend('Input spectrum', 'Filter response', 'Cutoff', 'Location', 'northeast');

subplot(2,1,2);
plot(f/1000, Yplot, 'Color', [0.2 0.45 0.75], 'LineWidth', 1.1);
hold on;
plot(f/1000, Hplot, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
xline(fp/1000, '--r', 'LineWidth', 1.2);
xline(-fp/1000, '--r', 'LineWidth', 1.2);
hold off;
title('Two-Sided DFT of Filtered Output');
xlabel('Frequency (kHz)');
ylabel(yLabelText);
xlim([-fs/2 fs/2]/1000);
ylim(yLimits);
grid on;
legend('Output spectrum', 'Filter response', 'Cutoff', 'Location', 'northeast');
end

plot_dft(audioIn, audioOut, fs, B, A, gain, 'db', fp);
