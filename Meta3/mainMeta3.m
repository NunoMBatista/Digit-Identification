audioSignals = cell(50, 1); % Each one of the 50 samples will contain the audio signals of each one of the digits
Fs = 48000;

for i = 1:50
    audioSignals{i} = preProcess(i - 1); % get the audio signals of each one of the digits from the sample
end

audioSignalMedians = cell(10, 1);

% for i = 1:10
%     curDigitSignals = [];
%     % Get max length of current digit
%     maxLength = 0;
%     for j = 1:50
%         if length(audioSignals{j}{i}) > maxLength
%             maxLength = length(audioSignals{j}{i});
%         end
%     end

%     for j = 1:50
%         % Add silence to the end of y
%         concatY = zeros(maxLength - length(audioSignals{j}{i}), 1);
%         audioSignals{j}{i} = [audioSignals{j}{i}; concatY];

%         curDigitSignals = [curDigitSignals; audioSignals{j}{i}'];
%     end
%     %audioSignalMedians{i} = mean(curDigitSignals, 1);
%     audioSignalMedians{i} = curDigitSignals(1, :) / max(abs(curDigitSignals(1, :)));

%     % Uncomment the following line to plot the median spectrogram
%     %audioSignalMedians{i} = audioSignalMedians{i} / max(abs(audioSignalMedians{i}));
%end

% Spectrogram parameters
windowSize = round(0.0064 * Fs); % Window size of 0.0064s
overlap = round(0.0032 * Fs); % Overlap of 0.0032s
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

    % Plot the spectrogram
    %imagesc(t*1000, f, 10*log10(abs(s)));
    imagesc(t*1000, f, log10(abs(s)));
    %imagesc(t*1000, f, abs(s));

    axis xy;
    title(['Digit ' num2str(i - 1)]);

    % Label the axes
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    
    % Colorbar label
    cb = colorbar;
    ylabel(cb, 'Power/Frequency (dB/Hz)');
end

spectralTimeFeatures = containers.Map();

% Potência média por banda de frequência
% Potência média por janela de tempo 
% Potência total por banda de frequência

for digit = 1:10
    for sample = 1:50
        currentDigitSignal = audioSignals{sample}{digit};

        % Get the stft from each sample
        [s, f, t] = spectrogram(audioSignals{sample}{digit}, hamming(windowSize), overlap, nfft, Fs, 'yaxis');
        
        % Get the time/frequency features
        spectralTimeFeatures(digit, sample) = mean(abs(s), 2);
    end
end