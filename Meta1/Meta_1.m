% Duration of the longest audio
maxRows = 0;
for i = 0:9
    [y, Fs] = audioread(sprintf("digits/%d_16_%d.wav", i, 0));
    [rows, cols] = size(y);
    if(maxRows < rows)
        maxRows = rows;
    end
end
maxDuration = maxRows/Fs;

% Create 9x3 matrix to store each one of the selected audio features 
% Energy(1), average pitch(2) and standard deviation(3)
audioFeatures = zeros(9, 3);

for i = 0:9
    % y is the audio signal
    % Fs is the sampling frequency
    [y, Fs] = audioread(sprintf("Samples/%d_16_%d.wav", i, 0));
    
    % Normalize the signal based on the maximum amplitude
    y = y / max(abs(y));

    % Ts is the sampling period
    Ts = 1 / Fs;

    % [rows, cols] are the original dimensions of y (there are #rows samples)
    [rows, cols] = size(y);

    % Calculate frame energy
    % Frame size in seconds
    frameSize = 0.01;
    % Get number of samples per frame
    frameSamples = round(frameSize * Fs);
    % Calculate number of frames
    numFrames = floor(rows / frameSamples);
    
    % Calculate frame energy, the sum of squares of the samples in the frame (not the integral because we are in the discrete domain)
    frameEnergy = zeros(numFrames, 1);
    % Iterate through every frame
    for j = 1:numFrames
        % The frame range is from the (j-1)th*frameSamples frame to the (j+1)th*frameSamples frame
        frame = y((j - 1)*frameSamples + 1:j*frameSamples);
        % Get the frame's energy
        frameEnergy(j) = sum(frame .^ 2);
    end
    % End frame energy calculation

    % Get total energy
    audioFeatures(i+1, 1) = sum(frameEnergy);

    % Get pitch
    pitchVector = pitch(y, Fs);
    audioFeatures(i+1, 2) = mean(pitchVector);
    
    % Get standard deviation (multiplied by 1000 to get a more readable value, even after normalization)
    audioFeatures(i+1, 3) = std(y)*1000;

    % Find first index of the first frame with energy above threshold
    energyThreshold = 0.5;
    startFrame = find(frameEnergy > energyThreshold, 1);
    
    % Get the TIME index of the first frame with enery above threshold
    startSample = (startFrame - 1) * frameSamples + 1;

    % Trim y in order to remove low energy values
    y = y(startSample:end);
    % Add silence to the end of y
    [curRows, curCols] = size(y);
    concatY = zeros(maxRows - curRows, 1);
    y = [y; concatY];
        
    
    [curRows, curCols] = size(y);
    % t is the time vector
    t = (0:maxRows-1);
    % t is the time vector in seconds (starts at the first frame that excedes the energy threshold)
    t = (t .* Ts);

    subplot(5, 2, i + 1);
    plot(t, y');
    xlabel('Time (s)');
    ylabel('Amplitude');

    label = sprintf("%d", i);
    title(label);


end
% Make a 3d scatter plot of the audio features
figure;
scatter3(audioFeatures(:, 1), audioFeatures(:, 2), audioFeatures(:, 3), [], 1:size(audioFeatures, 1), 'filled');
xlabel('Energy');
ylabel('Average Pitch');
zlabel('Standard Deviation');
title('Audio Features');
grid on;
% Label each color with the corresponding digit
digits = 0:9;
for i = 1:size(audioFeatures, 1)
    text(audioFeatures(i, 1), audioFeatures(i, 2), audioFeatures(i, 3), num2str(digits(i)));
end