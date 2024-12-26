% Assume emgMatrix is the input matrix with size 20 x 50000 (20 signals, 50000 timestamps)
[rows, cols] = size(emgMatrix);

% Create a figure for the plot
figure;
hold on;  % Hold the plot to overlay multiple signals

% Define colors for the first 10 signals
color = [0 0 1];;  % Get 10 distinct colors from the 'lines' colormap

% Loop through each EMG signal (each row of the matrix)
for i = 1
    % Extract the ith EMG signal (row of emgMatrix)
    signal = emgMatrix(i, :);
    
    % Plot the first 10 signals in different colors
    
    plot(signal, 'Color', color, 'LineWidth', 1.5, 'DisplayName', ['EMG ' num2str(i)]);
    
end

% Add title and labels
title('20 EMG Signals');
xlabel('Time (samples)');
ylabel('Amplitude');

% Add a legend for the first 10 signals
legend([arrayfun(@(i) ['EMG ' num2str(i)], 1:10, 'UniformOutput', false), {'Last 10 EMGs'}]);

% Release the hold on the plot
hold off;
