audioSignals = cell(50, 1); % Each one of the 50 samples will contain the audio signals of each one of the digits
Fs = 48000;

for i = 1:50
    audioSignals{i} = preProcess(i - 1); % get the audio signals of each one of the digits from the sample
end

audioSignalMedians = cell(10, 1);

% Spectrogram parameters
windowSize = round(0.0032 * Fs); % Window size of 0.0032s
overlap = round(0.0016 * Fs); % Overlap should be half the window size
nfft = 2^nextpow2(windowSize); % Number of points for the FFT

% Start plotting the spectrograms
figure;
for i = 1:10
    % Get the audio signal
    audioSignal = audioSignals{1}{i}; % Get the first example of each digit
    
    % Find the last non-zero element
    lastNonZero = find(audioSignal ~= 0, 1, 'last');
    
    % Trim the audio signal
    trimmedAudioSignal = audioSignal(1:lastNonZero);
    
    subplot(5, 2, i);
    [s, f, t] = spectrogram(trimmedAudioSignal, hamming(windowSize), overlap, nfft, Fs, 'yaxis'); % Get the spectrogram data

    % Plot the spectrogram in logarithmic scale
    imagesc(t*1000, f, 10*log10(abs(s)));

    axis xy;
    title(['Digit ' num2str(i - 1)]);

    % Label the axes
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    
    % Colorbar label
    cb = colorbar;
    ylabel(cb, 'Power/Frequency (dB/Hz)');
end

spectrogramFeatures = containers.Map();


meanPowerFreqBandDigit = cell(10, 1);
meanPowerTimeBandDigit = cell(10, 1);
peakPowerTimeBandDigit = cell(10, 1);
spectralFlatnessTimeBandDigit = cell(10, 1);

everyPeakPowerTimeBand = [];

maxTimeWindows = 530;
maxFreqBands = 257;
for digit = 1:10
    % These matrices will have 50 rows, one for each sample, and f or t columns, one for each frequency band or time window
    % We will then calculate the mean of each column to get the mean of each frequency band or time window over the 50 samples
    MPFB = []; % Mean Power Frequency Band
    MPTB = []; % Mean Power Time Band
    PPTB = []; % Peak Power Time Band
    SFTB = []; % Spectral Flatness Time Band

    for sample = 1:50
        % Find the last non-zero element
        lastNonZero = find(audioSignals{sample}{digit} ~= 0, 1, 'last');
        
        % Trim the audio signal
        trimmedAudioSignal = audioSignals{sample}{digit}(1:lastNonZero);
        
        % Get the STFT from each sample
        % s is a matrix containing the STFT of the signal (rows -> frequencies, columns -> time, values -> power) 
        % f is a vector containing the frequency values
        % t is a vector containing the time values
        [s, f, t] = spectrogram(trimmedAudioSignal, hamming(windowSize), overlap, nfft, Fs, 'yaxis');
        
         
        powerSpectrum = abs(s) .^ 2;
        

        
        % Calculate the mean power per frequency band
        curMPFB = mean(powerSpectrum, 2)'; % Column Matrix, transposed to be a row matrix
        curMPTB = mean(powerSpectrum, 1); % Row Matrix
        curPPTB = max(powerSpectrum, [], 1); % Row Matrix
        curSFTB = geomean(powerSpectrum) ./ mean(powerSpectrum); % Row Matrix

        % size(mean(powerSpectrum, 1))
        
        MPFB = [MPFB; curMPFB];
        
        % If the number of time windows is lower than the maximum, add silence to the end
        if size(curMPTB, 1) < maxTimeWindows
            curMPTB = [curMPTB zeros(1, maxTimeWindows - size(curMPTB, 2))];
        end
        MPTB = [MPTB; curMPTB];

        % If the number of time windows is lower than the maximum, add silence to the end
        if size(curPPTB, 1) < maxTimeWindows
            curPPTB = [curPPTB zeros(1, maxTimeWindows - size(curPPTB, 2))];
        end
        PPTB = [PPTB; curPPTB];

        % If the number of time windows is lower than the maximum, add silence to the end
        if size(curSFTB, 1) < maxTimeWindows
            curSFTB = [curSFTB zeros(1, maxTimeWindows - size(curSFTB, 2))];
        end
        SFTB = [SFTB; curSFTB];
        

    end

    everyPeakPowerTimeBand = [everyPeakPowerTimeBand; PPTB];

    % Turn the matrix into a vector by getting the mean of each column
    MPFB = mean(MPFB, 1);
    MPTB = mean(MPTB, 1);
    PPTB = mean(PPTB, 1);
    SFTB = mean(SFTB, 1);

    meanPowerFreqBandDigit{digit} = MPFB;
    meanPowerTimeBandDigit{digit} = MPTB;
    peakPowerTimeBandDigit{digit} = PPTB;
    spectralFlatnessTimeBandDigit{digit} = SFTB;    


end

% figure;
% plot3DScatterPlot(everyPeakPowerTimeBand, 'Peak Power Time Band');


for i = 1:10
    % Get the mean power per time band for only one digit
    lastNonZero = find(audioSignals{2}{i} ~= 0, 1, 'last');
    trimmedAudioSignal = audioSignals{2}{i}(1:lastNonZero);
    [s, f, t] = spectrogram(trimmedAudioSignal, hamming(windowSize), overlap, nfft, Fs, 'yaxis');
    powerSpectrum = abs(s) .^ 2;

%    Calculate the mean power per time band
    newSFTB = geomean(powerSpectrum) ./ mean(powerSpectrum);
    if size(newSFTB, 2) < maxTimeWindows
        newSFTB = [newSFTB zeros(1, maxTimeWindows - size(newSFTB, 2))];
    end

    % newMPTB = mean(powerSpectrum, 1);
    % if size(newMPTB, 2) < maxTimeWindows
    %     newMPTB = [newMPTB zeros(1, maxTimeWindows - size(newMPTB, 2))];
    % end



    maxCoef = 0;
    for j = 1:10
        r = corrcoef(newSFTB, meanPowerTimeBandDigit{j});
        if r(1, 2) > maxCoef
            maxCoef = r(1, 2);
            maxDigit = j;
        end
        %disp(['Correlation between digit 0 and digit ' num2str(i) ': ' num2str(r(1, 2))]);
    end

    disp(['The digit with the highest correlation with digit ' num2str(i) ' is digit ' num2str(maxDigit) ' with a correlation of ' num2str(maxCoef)]);
end

figure;
for i = 1:10   
    data = meanPowerFreqBandDigit{i}(:);

    % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
    plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
    hold on;
end
legend('Location', 'Best');

figure;
for i = 1:10   
    data = meanPowerTimeBandDigit{i}(:);

    % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
    plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
    hold on;
end
legend('Location', 'Best');

figure;
for i = 1:10   
    data = peakPowerTimeBandDigit{i}(:);

    % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
    plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
    hold on;
end
legend('Location', 'Best');

figure; 
for i = 1:10   
    data = spectralFlatnessTimeBandDigit{i}(:);

    % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
    plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
    hold on;
end
legend('Location', 'Best');









function plot3DScatterPlot(data, plotTitle)
    % Get the number of rows and columns
    [nRows, nCols] = size(data);

    % Create a meshgrid for the row and column indices
    [X, Y] = meshgrid(1:nCols, 1:nRows);

    % Flatten the matrices
    X = X(:);
    Y = Y(:);
    Z = data(:);

    colors = jet(10);

    colorIndices = ceil(Y / 50);

    % Create the 3D scatter plot
    scatter3(X, Y, Z, 10, colors(colorIndices, :), 'filled');
    xlabel('Row Index');
    ylabel('Column Index');
    zlabel('Value');
    title(plotTitle);
    grid on;

    hold on;
    for i = 1:10
        scatter3(nan, nan, nan, 10, colors(i, :), 'filled', 'DisplayName', ['Digit ' num2str(i - 1)]);
    end
    hold off;

    legend('Location', 'Best');
end