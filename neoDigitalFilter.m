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
B = []; A = [];                                                             %crate empty matrices for coeficcients B and A
for i = 1:np                                                                %for loop over np
    if used(i), continue; end                                               %check if pole was already checked
    if abs(imag(poles(i))) < 1e-12 * max(1, abs(poles(i)))                  %check if pole is real
        B = [B; real(residue(i))];                                          %add coefficient to B array
        A = [A; 1-exp(real(poles(i)))];                                     %add coefficient to A array
        used(i) = true;                                                     %Check pole to used
    else
        j = find(~used & abs(poles-conj(poles(i)))<1e-9*abs(poles(i)), 1);  %finds the location of the complex conjugate of the pole
        used([i j]) = true;                                                 %sets used to that to true
        ez = exp(poles(i));                                                 %exponential of pole i         
        B = [B; 2*real(residue(i)), -2*real(residue(i)*conj(ez)), 0];       %coefficient of b
        A = [A; 1, -2*real(ez), abs(ez)^2];                                 %coefficient of a
    end
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

%plots check later
t = (0:size(audioIn,1)-1) / fs;

figure('Color', 'white', 'Position', [100, 50, 1100, 850]);

subplot(3, 1, 1);
plot(t, audioIn(:,1), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); 
hold on;
plot(t, audioOut(:,1), 'Color', [0 0.4 0.8], 'LineWidth', 1.0); 
grid on;
title(sprintf('Time-Domain Waveforms (Noisy Input vs. Order 20 Elliptic, Impulse Invariance @ %g kHz)', fp/1000));
xlabel('Time (s)'); 
ylabel('Amplitude');
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
H  = gain * ones(1, nf);
for s = 1:size(B,1)
    H = H + freqz(B(s,:), A(s,:), w);
end
figure('Color','white');
plot(w/pi*fs/2/1000, 20*log10(abs(H)+eps), 'LineWidth', 1.2); grid on;
xline(fp/1000, '--r');
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('Impulse-Invariant Elliptic Filter Response');
ylim([-20 5]);

        