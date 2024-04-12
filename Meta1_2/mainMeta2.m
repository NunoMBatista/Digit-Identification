freqNum = 3000; % number of frequencies to plot
medianAmpSpectrumMat = zeros(freqNum, 10); % matrix to store the median amplitude spectrum of each digit
Q25AmpSpectrumMat = zeros(freqNum, 10); % matrix to store the 25th percentile amplitude spectrum of each digit
Q75AmpSpectrumMat = zeros(freqNum, 10); % matrix to store the 75th percentile amplitude spectrum of each digit
meanAmpSpectrumMat = zeros(freqNum, 10); % matrix to store the mean amplitude spectrum of each digit
everyAmpSpectrumMat = cell(10, 1); % matrix to store the amplitude spectrum of each digit

audioSignals = cell(50, 1); % Each one of the 50 samples will contain the audio signals of each one of the digits
for i = 1:50
    audioSignals{i} = preProcess(i-1); % get the audio signals of each one of the digits from the sample
end

[medianAmpSpectrumMat, Q25AmpSpectrumMat, Q75AmpSpectrumMat, meanAmpSpectrumMat, everyAmpSpectrumMat] = digitsAmpSpectrums(audioSignals, 'default');

% digitsAmpSpectrums('hamming');
% digitsAmpSpectrums('hann');
% digitsAmpSpectrums('blackman');

spectralFeaturesStrings = {'Spectral Mean'};
spectralFeatures = containers.Map(); % Hashmap containing the spectral features of each digit


% Array with the spectral mean of each digit
spectralMean = cell(10, 1);
% Array with the spectral flux of each digit
spectralFlux = cell(10, 1);
% Array with the spectral spread of each digit
spectralSpread = cell(10, 1);
% Array with the spectral max absolute value of each digit
spectralMaxAbs = cell(10, 1);
% Array with the spectral flatness of each digit
spectralFlatness = cell(10, 1);
% Array with the spectral variance of each digit
spectralVariance = cell(10, 1);


curSampleSpectralMean = cell(50, 1);
curSampleSpectralFlux = cell(50, 1);
curSampleSpectralSpread = cell(50, 1);
curSampleSpectralMaxAbs = cell(50, 1);
curSampleSpectralFlatness = cell(50, 1);
curSampleSpectralVariance = cell(50, 1);


for digit = 1:10
    for sample = 1:50
        curSampleSpectralMean{sample} = mean(everyAmpSpectrumMat{digit}(:, sample), 1); 
        curSampleSpectralFlux{sample} = sum(abs(diff(everyAmpSpectrumMat{digit}(:, sample))), 1);
        curSampleSpectralSpread{sample} = std(everyAmpSpectrumMat{digit}(:, sample), 1);
        curSampleSpectralMaxAbs{sample} = max(abs(everyAmpSpectrumMat{digit}(:, sample)), [], 1);
        curSampleSpectralFlatness{sample} = geomean(everyAmpSpectrumMat{digit}(:, sample)) / mean(everyAmpSpectrumMat{digit}(:, sample));
        curSampleSpectralVariance{sample} = var(everyAmpSpectrumMat{digit}(:, sample), 1);
    end
    spectralMean{digit} = curSampleSpectralMean;
    spectralFlux{digit} = curSampleSpectralFlux;
    spectralSpread{digit} = curSampleSpectralSpread;
    spectralMaxAbs{digit} = curSampleSpectralMaxAbs;
    spectralFlatness{digit} = curSampleSpectralFlatness;
    spectralVariance{digit} = curSampleSpectralVariance;
end

spectralFeatures('Spectral Mean') = spectralMean;
spectralFeatures('Spectral Flux') = spectralFlux;
spectralFeatures('Spectral Spread') = spectralSpread;
spectralFeatures('Spectral Max Abs') = spectralMaxAbs;
spectralFeatures('Spectral Variance') = spectralVariance;
spectralFeatures('Spectral Flatness') = spectralFlatness;

% Get the spectral 
figure;
% Make a boxplot of the spectral mean of each digit
subplot(2, 3, 1);
boxplotFeature(spectralFeatures('Spectral Mean'), 'Spectral Mean', 'Boxplot of Spectral Means for Each Digit');

% Make a boxplot of the spectral centroid of each digit
subplot(2, 3, 2);
boxplotFeature(spectralFeatures('Spectral Flux'), 'Spectral Flux', 'Boxplot of Spectral Flux for Each Digit');

% Make a boxplot of the spectral spread of each digit
subplot(2, 3, 3);
boxplotFeature(spectralFeatures('Spectral Spread'), 'Spectral Spread', 'Boxplot of Spectral Spread for Each Digit');

% Make a boxplot of the spectral max absolute value of each digit
subplot(2, 3, 4);
boxplotFeature(spectralFeatures('Spectral Max Abs'), 'Spectral Max Abs', 'Boxplot of Spectral Max Abs for Each Digit');

% Make a boxplot of the spectral flatness of each digit
subplot(2, 3, 5); 
boxplotFeature(spectralFeatures('Spectral Flatness'), 'Spectral Flatness', 'Boxplot of Spectral Flatness for Each Digit');

% Make a boxplot of the spectral variance of each digit
subplot(2, 3, 6);
boxplotFeature(spectralFeatures('Spectral Variance'), 'Spectral Variance', 'Boxplot of Spectral Variance for Each Digit');

% Make a boxplot of the spectral flatness of each digit
subplot(2, 3, 5); 
boxplotFeature(spectralFeatures('Spectral Flatness'), 'Spectral Flatness', 'Boxplot of Spectral Flatness for Each Digit');