audioSignals = cell(50, 1); % Each one of the 50 samples will contain the audio signals of each one of the digits
Fs = 48000;

for i = 1:50
    audioSignals{i} = preProcess(i - 1); % get the audio signals of each one of the digits from the sample
end

audioSignalMedians = cell(10, 1);

% Spectrogram parameters
windowSize = round(0.0032 * Fs); % Window size of 0.0032ms
overlap = round(0.0016 * Fs); % Overlap should be half the window size
nfft = 2^nextpow2(windowSize); % Number of points for the FFT

% Start plotting the spectrograms
% figure;
% for i = 1:10
%     % Get the audio signal
%     audioSignal = audioSignals{1}{i}; % Get the first example of each digit
    
%     % Find the last non-zero element
%     lastNonZero = find(audioSignal ~= 0, 1, 'last');
    
%     % Trim the audio signal
%     trimmedAudioSignal = audioSignal(1:lastNonZero);
    
%     subplot(5, 2, i);
%     [s, f, t] = spectrogram(trimmedAudioSignal, hamming(windowSize), overlap, nfft, Fs, 'yaxis'); % Get the spectrogram data

%     % Plot the spectrogram in logarithmic scale
%     imagesc(t*1000, f, 10*log10(abs(s)));

%     axis xy;
%     title(['Digit ' num2str(i - 1)]);

%     % Label the axes
%     xlabel('Time (ms)');
%     ylabel('Frequency (Hz)');
    
%     % Colorbar label
%     cb = colorbar;
%     ylabel(cb, 'Power/Frequency (dB/Hz)');
% end

spectrogramFeatures = containers.Map();

maxTimeWindows = 200;
maxFreqBands = 129;

meanPowerFreqBandDigit = cell(10, 1);
meanPowerTimeBandDigit = cell(10, 1);
peakPowerTimeBandDigit = cell(10, 1);
spectralFlatnessTimeBandDigit = cell(10, 1);
powerSTDTimeBandDigit = cell(10, 1);
spectralFluxTimeBandDigit = cell(10, 1);
spectralRollOffTimeBandDigit = cell(10, 1);
powerSTDFreqBandDigit = cell(10, 1);

spectralCentroidDigit = zeros(10, 50);

% These matrices will have 50 * 10 rows, one for each sample of each digit and f or t columns, one for each frequency band or time window
% They will later be used to get 3d scatter plots
everyPeakPowerTimeBand = [];
everySpectralFluxTimeBand = [];

everySpectralRollOffTimeBand = [];

for digit = 1:10
    % These matrices will have 50 rows, one for each sample, and f or t columns, one for each frequency band or time window
    % We will then calculate the mean of each column to get the mean of each frequency band or time window over the 50 samples
    % Preallocate matrices
    MPFB = zeros(50, maxFreqBands); % Mean Power per Frequency Band
    MPTB = zeros(50, maxTimeWindows); % Mean Power per Time Band
    PPTB = zeros(50, maxTimeWindows); % Peak Power per Time Band
    SFTB = zeros(50, maxTimeWindows); % Spectral Flatness per Time Band
    PSTDTB = zeros(50, maxTimeWindows); % Power STD per Time Band
    SFluxTB = zeros(50, maxTimeWindows); % Spectral Flux per Time Band
    SROTB = zeros(50, maxTimeWindows); % Spectral Roll Off per Time Band
    PSTDFB = zeros(50, maxFreqBands); % Power STD per Frequency Band

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
        
        % Apply logarithmic scale to the power values
        powerSpectrum = abs(s) .^ 2;
        powerSpectrum = abs(10 * log10(powerSpectrum)); 
        
        % When appending to the overall matrices, we need to make sure that the number of time windows is the same for all samples
        
        % Calculate the mean power per frequency band of the current sample
        curMPFB = mean(powerSpectrum, 2)'; % Column Matrix, transposed to be a row matrix
        MPFB(sample, :) = curMPFB;
        
        % Calculate the mean power per time band of the current sample
        curMPTB = mean(powerSpectrum, 1); % Row Matrix

        MPTB(sample, :) = addSilence(curMPTB, maxTimeWindows);

        % Calculate the peak power per time band of the current sample
        curPPTB = max(powerSpectrum, [], 1); % Row Matrix
        PPTB(sample, :) = addSilence(curPPTB, maxTimeWindows);
        
        % Calculate the spectral flatness per time band of the current sample
        curSFTB = geomean(powerSpectrum) ./ mean(powerSpectrum); % Row Matrix
        SFTB(sample, :) = addSilence(curSFTB, maxTimeWindows);
        
        % Calculate the power STD per time band of the current sample
        curPSTDTB = std(powerSpectrum, 0, 1); % Row Matrix
        PSTDTB(sample, :) = addSilence(curPSTDTB, maxTimeWindows);
        
        % Calculate the spectral centroid of the sample
        spectralCentroidDigit(digit, sample) = sum(f .* mean(abs(s), 2)) / sum(mean(abs(s), 2));
        
        % Calculate the Spectral Flux
        curSFluxTB = zeros(1, size(powerSpectrum, 2));
        curSFluxTB(1) = 0;
        for frame = 2:size(powerSpectrum, 2)
            curSFluxTB(frame) = sum((powerSpectrum(:, frame) - powerSpectrum(:, frame - 1)).^2);
        end
        SFluxTB(sample, :) = addSilence(curSFluxTB, maxTimeWindows);

        % Calculate the Spectral Roll Off
        curSROTB = zeros(1, size(powerSpectrum, 2));
        cumulativeSum = cumsum(powerSpectrum, 1);
        totalSum = sum(powerSpectrum, 1);
        for frame = 1:size(powerSpectrum, 2)
            [~, rollOffIndex] = min(abs(cumulativeSum(:, frame) - 0.85 * totalSum(frame)));
            curSROTB(frame) = f(rollOffIndex);
        end
        SROTB(sample, :) = addSilence(curSROTB, maxTimeWindows);

        % Calculate the Power STD per Frequency Band
        curPSTDFB = std(powerSpectrum, 0, 2)'; % Column Matrix, transposed to be a row matrix
        PSTDFB(sample, :) = curPSTDFB;

    end
    

    everyPeakPowerTimeBand = [everyPeakPowerTimeBand; PPTB];
    everySpectralFluxTimeBand = [everySpectralFluxTimeBand; SFluxTB];
    everySpectralRollOffTimeBand = [everySpectralRollOffTimeBand; SROTB];

    % Turn the matrix into a vector by getting the mean of each column
    MPFB = mean(MPFB, 1);
    MPTB = mean(MPTB, 1);
    PPTB = mean(PPTB, 1);
    SFTB = mean(SFTB, 1);
    PSTDTB = mean(PSTDTB, 1);
    SFluxTB = mean(SFluxTB, 1);
    SROTB = mean(SROTB, 1);
    PSTDFB = mean(PSTDFB, 1);

    meanPowerFreqBandDigit{digit} = MPFB;
    meanPowerTimeBandDigit{digit} = MPTB;
    peakPowerTimeBandDigit{digit} = PPTB;
    spectralFlatnessTimeBandDigit{digit} = SFTB;    
    powerSTDTimeBandDigit{digit} = mean(PSTDTB, 1);
    spectralFluxTimeBandDigit{digit} = SFluxTB;
    spectralRollOffTimeBandDigit{digit} = SROTB;
    powerSTDFreqBandDigit{digit} = PSTDFB;

end

spectrogramFeatures('Mean Power per Frequency Band') = meanPowerFreqBandDigit;
spectrogramFeatures('Mean Power per Time Band') = meanPowerTimeBandDigit;
spectrogramFeatures('Peak Power per Time Band') = peakPowerTimeBandDigit;
spectrogramFeatures('Spectral Flatness per Time Band') = spectralFlatnessTimeBandDigit;
spectrogramFeatures('Power STD per Time Band') = powerSTDTimeBandDigit;
spectrogramFeatures('Spectral Flux per Time Band') = spectralFluxTimeBandDigit;
spectrogramFeatures('Spectral Roll-Off per Time Band') = spectralRollOffTimeBandDigit;

spectrogramFeatures('Spectral Centroid') = spectralCentroidDigit;


spectrogramFeatures('Power STD per Frequency Band') = powerSTDFreqBandDigit;

featuresStrings = {
    'Mean Power per Frequency Band', 
    'Mean Power per Time Band', 
    'Peak Power per Time Band', 
    'Spectral Flatness per Time Band', 
    'Power STD per Time Band', 
    'Spectral Flux per Time Band', 
    'Spectral Roll-Off per Time Band',
    'Power STD per Frequency Band'
};

reducedMeanPowerTimeBandDigit = extractBestWindow(everyPeakPowerTimeBand, 5, 10, 1);

reducedFeature = cell2mat(reducedMeanPowerTimeBandDigit)
size(reducedFeature)

figure;
allFeatureValues = []
digitLabels = []

for digit = 1:10
    currentFeatureValues = reducedFeature(i);

    allFeatureValues = [allFeatureValues; currentFeatureValues];

    digitLabels = [digitLabels; repmat(digit, length(currentFeatureValues), 1)]
end    

boxplot(allFeatureValues, digitLabels)

xticklabels(0:9)

title(plotTitle)
xlabel('Digit')
ylabel




% for i = 1:length(featuresStrings)
%     figure;
%     curFeat = spectrogramFeatures(featuresStrings{i});
%     for j = 1:10   
%         data = curFeat{j}(:);
%         % data = meanPowerFreqBandDigit{i}(:);
    
%         % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
%         plot(data, 'DisplayName', ['Digit ' num2str(j - 1)]);
%         hold on;
%     end
%     title(featuresStrings{i});
%     colormap(jet(10));
%     legend('Location', 'Best');
% end

% figure;
% plot3DScatterPlot(everyPeakPowerTimeBand, 'Peak Power Time Band');

% figure;
% plot3DScatterPlot(everySpectralFluxTimeBand, 'Spectral Flux Time Band');

% figure;
% plot3DScatterPlot(everySpectralRollOffTimeBand, 'Spectral Roll Off Time Band');

% for i = 1:10
%     % Get the mean power per time band for only one digit
%     lastNonZero = find(audioSignals{10}{i} ~= 0, 1, 'last');
%     trimmedAudioSignal = audioSignals{10}{i}(1:lastNonZero);
%     [s, f, t] = spectrogram(trimmedAudioSignal, hamming(windowSize), overlap, nfft, Fs, 'yaxis');
%     powerSpectrum = abs(s) .^ 2;
%     powerSpectrum = abs(10 * log10(powerSpectrum));

%     % Calculate the mean power per time band
%     newMPTB = mean(powerSpectrum, 1);
    
%     % Calculate the spectral flatness per time band
%     newSFTB = geomean(powerSpectrum) ./ mean(powerSpectrum);
%     newSFTB = addSilence(newSFTB, maxTimeWindows);
 
%     maxCoef = 0;
%     for j = 1:10
%         r = corrcoef(newSFTB, meanPowerTimeBandDigit{j});
%         if r(1, 2) > maxCoef
%             maxCoef = r(1, 2);
%             maxDigit = j;
%         end
%         %disp(['Correlation between digit 0 and digit ' num2str(i) ': ' num2str(r(1, 2))]);
%     end

%     disp(['The digit with the highest correlation with digit ' num2str(i) ' is digit ' num2str(maxDigit) ' with a correlation of ' num2str(maxCoef)]);
% end

% figure;
% for i = 1:10   
%     data = meanPowerFreqBandDigit{i}(:);

%     % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
%     plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
%     hold on;
% end
% legend('Location', 'Best');

% figure;
% for i = 1:10   
%     data = meanPowerTimeBandDigit{i}(:);

%     % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
%     plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
%     hold on;
% end
% legend('Location', 'Best');

% figure;
% for i = 1:10   
%     data = peakPowerTimeBandDigit{i}(:);

%     % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
%     plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
%     hold on;
% end
% legend('Location', 'Best');

% figure; 
% for i = 1:10   
%     data = spectralFlatnessTimeBandDigit{i}(:);

%     % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
%     plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
%     hold on;
% end
% legend('Location', 'Best');

% figure; 
% for i = 1:10   
%     data = powerSTDTimeBandDigit{i}(:);

%     % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
%     plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
%     hold on;
% end
% legend('Location', 'Best');

% figure;
% for i = 1:10   
%     data = spectralFluxTimeBandDigit{i}(:);

%     % Plot meanPowerFreqBandDigit for each digit in the same 2d plot, diferentiating digits by color
%     plot(data, 'DisplayName', ['Digit ' num2str(i - 1)]);
%     hold on;
% end
% legend('Location', 'Best');




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




function [newData] = addSilence(data, maxTimeWindows)
    if size(data, 2) < maxTimeWindows
        newData = [data zeros(1, maxTimeWindows - size(data, 2))];
    else
        newData = data(1:maxTimeWindows);
    end
end


function [reducedFeatureDigit] = extractBestWindow(everyWindowFeatureDigit, low, high, reductionType)
    reducedFeatureDigit = cell(10, 50);
    % Reduction types
    % 1 -> Mean
    % 2 -> Median
    % 3 -> Max
    % 4 -> Standard Deviation

    if reductionType == 1
        for dig = 1:10
            for samp = 1:50
                reducedFeatureDigit{dig, samp} = mean(everyWindowFeatureDigit(dig*samp, low:high));
            end
        end
    end
end
