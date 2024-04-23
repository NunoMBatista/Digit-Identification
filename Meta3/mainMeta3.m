audioSignals = cell(50, 1); % Each one of the 50 samples will contain the audio signals of each one of the digits

for i = 1:50
    audioSignals{i} = preProcess(i - 1); % get the audio signals of each one of the digits from the sample
end

disp('olaaaaaaaaaaaaaaaaaaaa\n');

audioSignalMedians = cell(10, 1);

% Get the median of each row in the audioSignals matrix


% Plot the median audio signals
figure;
for i = 1:10
    subplot(2, 5, i);
    plot(audioSignalMedians{i});
    title(sprintf("Digit %d", i - 1));
end