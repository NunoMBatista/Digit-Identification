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