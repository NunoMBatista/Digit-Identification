Fs = 48000;

audioSignals = preProcess(1)

figure;
for i = 1:10
    subplot(5,2,i);
    fftData = fft(audioSignals{i});
    fftData = fftData(1:floor(length(fftData)/2)); % only take the first half of the data because the second half is just a mirror image

    ampSpectrum = abs(fftData); % get the magnitude of the complex fourier series coefficients
    ampSpectrum = ampSpectrum/length(fftData); % normalize by the number of samples

    plot(ampSpectrum(1:3000), 'Color', 'red'); % set the color to yellow
    title(['Amplitude Spectrum of Digit ', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
end
