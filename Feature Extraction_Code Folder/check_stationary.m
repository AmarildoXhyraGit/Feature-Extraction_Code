% Assuming signal is the input EMG signal
signal = emgMatrix(1,:)
window_size = 20000;  % Define a window size
num_windows = floor(length(signal) / window_size);

means = zeros(1, num_windows);
variances = zeros(1, num_windows);

for i = 1:num_windows
    % Extract the i-th window from the signal
    window = signal((i-1)*window_size + 1 : i*window_size);
    means(i) = mean(window);
    variances(i) = var(window);
end

% Plot the mean and variance over time
figure;
subplot(2,1,1);
plot(means);
title('Mean over time');
xlabel('Window');
ylabel('Mean');

subplot(2,1,2);
plot(variances);
title('Variance over time');
xlabel('Window');
ylabel('Variance');
