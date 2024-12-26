% Assuming 'EMG' is a matrix with 20 rows (signals) and 50000 columns (time points)

% Initialize the feature matrix to store the means of wavelet coefficients
num_signals = size(emgMatrix, 1);
wavelet_means = zeros(num_signals, 1);

% Loop through each signal and calculate the continuous wavelet transform (CWT)
for i = 1:num_signals
    % Compute the continuous wavelet transform (CWT)
    [wt, ~] = cwt(emgMatrix(i, :), 'amor'); % 'amor' is the Morlet wavelet (adjust if needed)
    
    % Calculate the mean of the absolute values of the wavelet coefficients
    wavelet_means(i) = mean(abs(wt(:)));
end
