audioFeatures = zeros(50, 10, 5);

% Start by plotting the first sample of each digit (id = 0)
figure(1);
getSpectrum(1, 1);
sgtitle('Signal examples');
% i + 1 because MATLAB uses 1-based indexing but the audio file names use 0-based indexing

