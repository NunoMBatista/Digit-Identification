%cd Samples
clear;
for i = 0:9
    % y is the audio signal
    % Fs is the sampling frequency
    [y, Fs] = audioread(sprintf("Samples/%d_16_%d.wav", i, 0));
    
    % Ts is the sampling period
    Ts = 1 / Fs;
    % [rows, cols] are the dimensions of y
    [rows, cols] = size(y);
    % t is the time vector
    t = 0:rows - 1;
    % t is the time vector in seconds
    t = t .* Ts;
    
    % Frame size in seconds
    frameSize = 0.02;
    % Calculate frame energy
    frameSamples = round(frameSize * Fs);

    % Calculate frame energy, the sum of squares of the samples in the frame (discrete domain)
    numFrames = floor(ws / frameSamples);
    frameEnergy = zeros(numFrames, 1);
    for j = 1:numFrames
        frame = y((j-1)*frameSamples + 1:j*frameSamples);
        frameEnergy(j) = sum(frame.^2);
    end


    % Find first frame with energy above threshold
    energyThreshold = 0.001; % Adjust this threshold as needed
    startFrame = find(frameEnergy > energyThreshold, 1);

    startSample = (startFrame - 1) * frameSamples + 1;

    % Trim audio signal from startFrame onwards
    y = y(startSample:end);
    t = t(startSample:end);
    
    subplot(5, 2, i + 1);
    plot(t, y');
    label = sprintf("%d", i);
    title(label);
end