freqNum = 3000; % number of frequencies to plot
medianAmpSpectrumMat = zeros(freqNum, 10); % matrix to store the median amplitude spectrum of each digit
Q25AmpSpectrumMat = zeros(freqNum, 10); % matrix to store the 25th percentile amplitude spectrum of each digit
Q75AmpSpectrumMat = zeros(freqNum, 10); % matrix to store the 75th percentile amplitude spectrum of each digit
everyAmpSpectrumMat = cell(10, 1); % matrix to store the amplitude spectrum of each digit

audioSignals = cell(50, 1); % Each one of the 50 samples will contain the audio signals of each one of the digits
for i = 1:50
    audioSignals{i} = preProcess(i-1); % get the audio signals of each one of the digits from the sample
end

meanAmpSpectrumMat = zeros(freqNum, 10); % matrix to store the mean amplitude spectrum of each digit
[medianAmpSpectrumMat, Q25AmpSpectrumMat, Q75AmpSpectrumMat, meanAmpSpectrumMat, everyAmpSpectrumMat] = digitsAmpSpectrums(audioSignals, 'default');

% digitsAmpSpectrums('hamming');
% digitsAmpSpectrums('hann');
% digitsAmpSpectrums('blackman');

spectralFeaturesStrings = {'Spectral Mean'};

spectralFeatures = containers.Map(); % Hashmap containing the spectral features of each digit

% Get spectral mean of each digit
spectralMean = cell(10, 1);
curSampleSpectralMean = cell(50, 1);
for digit = 1:10
    for sample = 1:50
        curSampleSpectralMean{sample} = mean(everyAmpSpectrumMat{digit}(:, sample), 1);
    end
    spectralMean{digit} = curSampleSpectralMean;
end
spectralFeatures('Spectral Mean') = spectralMean;


% Make a boxplot of the spectral mean of each digit
boxplotFeature(spectralFeatures('Spectral Mean'), 'Spectral Mean', 'Boxplot of Spectral Means for Each Digit');
