Fs = 48000;

audioSignals = cell(50, 1);

for i = 1:50
    audioSignals{i} = preProcess(i-1);
end
audioSignals;

% Each index of audioSignals contains the audio signals of each one of the digits from that sample
digit5 = cell(50, 1);
for i = 1:50
    digit5{i} = audioSignals{i}{6};
end

fftData = fft(digit5{2});
fftData = 1:floor(length(fftData)/2);
figure;
ampSpectrum = abs(fftData);
ampSpectrum = ampSpectrum/length(fftData);
plot(ampSpectrum(1:3000), 'Color', 'red');


% figure;
% for i = 1:10
%     subplot(5,2,i);
%     fftData = fft(audioSignals{i});
%     fftData = fftData(1:floor(length(fftData)/2)); % only take the first half of the data because the second half is just a mirror image

%     ampSpectrum = abs(fftData); % get the magnitude of the complex fourier series coefficients
%     ampSpectrum = ampSpectrum/length(fftData); % normalize by the number of samples

%     plot(ampSpectrum(1:3000), 'Color', 'red'); % set the color to yellow
%     title(['Amplitude Spectrum of Digit ', num2str(i)]);
%     xlabel('Frequency (Hz)');
%     ylabel('Amplitude');
% end
