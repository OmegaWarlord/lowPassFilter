[y1, fs] = audioread('/untitled2.wav');
y2 = filter(filter2filt,y1);

audiowrite('filtered.wav', y2, fs);


% If stereo, convert to mono by taking the first channel
if size(y1, 2) > 1, y1 = y1(:, 1); end
if size(y2, 2) > 1, y2 = y2(:, 1); end

% 2. Create time vectors for both signals
t1 = (0:length(y1)-1) / fs;
t2 = (0:length(y2)-1) / fs;

% 3. Create a figure with 4 subplots (2 waveforms + 2 spectrograms)
figure('Color', 'white');

% --- Audio 1: Waveform ---
subplot(2, 2, 1);
plot(t1, y1);
title('Audio 1: Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- Audio 1: Spectrogram ---
subplot(2, 2, 2);
window = 512;
noverlap = 256;
nfft = 512;
spectrogram(y1, window, noverlap, nfft, fs, 'yaxis');
title('Audio 1: Spectrogram');

% --- Audio 2: Waveform ---
subplot(2, 2, 3);
plot(t2, y2, 'r'); % 'r' for red color differentiation
title('Audio 2: Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- Audio 2: Spectrogram ---
subplot(2, 2, 4);
spectrogram(y2, window, noverlap, nfft, fs, 'yaxis');
title('Audio 2: Spectrogram');