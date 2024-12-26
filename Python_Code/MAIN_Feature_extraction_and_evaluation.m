clear variables;

N = 50000;  
numOfSignals = 20;
samplingFrequency = 10000;
t = 0:1/samplingFrequency:1;
emgMatrix = zeros(numOfSignals,N);

% creating the EMG matrix of all the signals
for i = 1:numOfSignals
    filename = sprintf('EMG%d.mat', i);
    data = load(filename);
    row = data.EMG;
    emgMatrix(i, :) = row;
end

                                     % Features
%| Energy
  energy_per_row = sum(abs(emgMatrix).^2, 2);

%| Power
  fs = 10000;
  num_samples = size(emgMatrix, 2);
  T = num_samples / fs; 
  power_per_row = energy_per_row / T;

%| Mean 
  mean_values = mean(emgMatrix ,2);  
    
%| Variance
variance_values = var(emgMatrix,0, 2);

%| Standard Deviation 
standard_deviation_value = std(emgMatrix, 0, 2);

%| Skewness
skewness_value = skewness(emgMatrix, 0, 2);

%| Kurtosis
kurtosis_value= kurtosis(emgMatrix, 0, 2);

%| Peak Amplitude
peak_amplitude_values = max(emgMatrix, [], 2);  

%| Autocorrelation
autocorr_values = zeros(20, N); 
for i = 1:20
    [acor, lag] = xcorr(emgMatrix(i, :), 'coeff');  
    autocorr_values(i, :) = [acor(N:end)];
end
auc_scores = mean(autocorr_values, 2);  

%| Moment Coefficient of Skewness
moment_skewness_values = zeros(20, 1);  
for i = 1:20 
    E_X = mean(emgMatrix(i, :));             
    E_X2 = mean(emgMatrix(i, :).^2);           
    E_X3 = mean(emgMatrix(i, :).^3);          
    moment_skewness_values(i) = E_X3 / (E_X2^(3/2));
end

%| Moment Coefficient of Kurtosis
moment_kurtosis_values = zeros(20, 1);  
for i = 1:20 
    E_X2 = mean(emgMatrix(i, :).^2);        
    E_X4 = mean(emgMatrix(i, :).^4);       
    moment_kurtosis_values(i) = E_X4 / (E_X2^2);
end

%| Entropy
emgMatrixNonZero = abs(emgMatrix);
normalizedMatrix = emgMatrixNonZero ./ (sum(emgMatrixNonZero, 2) + eps);
entropyValues = -sum(normalizedMatrix .* log2(normalizedMatrix + eps), 2);

%| High/Low ratio
lowFreq = [0, 100];  
highFreq = [100, 500];
hlRatio = zeros(numOfSignals, 1);
for i = 1:numOfSignals
    powerSpectrum = abs(fft(emgMatrix(i, :))).^2 / N;
    freq = (0:N-1) * (samplingFrequency / N);
    lowFreqPower = sum(powerSpectrum(freq >= lowFreq(1) & freq <= lowFreq(2)));
    highFreqPower = sum(powerSpectrum(freq > highFreq(1) & freq <= highFreq(2)));
    hlRatio(i) = highFreqPower / (lowFreqPower + eps);
end

%| Median frequency
medianFrequencies = zeros(numOfSignals, 1);
for i = 1:numOfSignals
    P = abs(fft(emgMatrix(i, :))).^2;
    P = P(1:floor(end/2));
    f = (0:length(P)-1) * (samplingFrequency / length(emgMatrix(i, :))); 
    P_total = sum(P);
    cumulativePower = cumsum(P);   
    medianFrequencies(i) = f(find(cumulativePower >= P_total / 2, 1));
end

%| Relative energy per frequency band
frequencyBands = [0 4;   
                  4 8;   
                  8 12;  
                  12 30; 
                  30 50]; 
numOfBands = size(frequencyBands, 1);
relativeEnergy = zeros(numOfSignals, numOfBands);
for i = 1:numOfSignals
    Y = fft(emgMatrix(i, :));
    P = abs(Y).^2 / N; 
    f = (0:(N/2)) * (samplingFrequency / N);  
    P = P(1:length(f));
    totalPower = sum(P);   
    for j = 1:numOfBands
        bandIndices = find(f >= frequencyBands(j, 1) & f < frequencyBands(j, 2));
        bandEnergy = sum(P(bandIndices));
        relativeEnergy(i, j) = bandEnergy / totalPower;
    end
end
avg = mean(relativeEnergy, 2);

%| Root mean square
rmsValues = sqrt(mean(emgMatrix.^2, 2));

%| Time reversibility
t1 = 1;
timeReversibility = zeros(numOfSignals, 1);
for i = 1:numOfSignals
    signal = emgMatrix(i, :);
    sum_tr = 0;
    for n = (t1+1):N
        sum_tr = sum_tr + (signal(n) - signal(n-t1))^3;
    end
    timeReversibility(i) = sum_tr / (N-t1);
end

%| zScore
eps1 = 1e-10;
zScoreMatrix = (emgMatrix - mean(emgMatrix, 2)) ./ (std(emgMatrix, 0, 2) + eps1);
meanZScore = mean(zScoreMatrix, 2); 

% Variables for PSD Features
fs = 10000;  
[rows, cols] = size(emgMatrix);
noverlap = 250;         
nfft = 1024;   
window = hamming(500);  

%| Average PSD
averages_psd = zeros(rows, 1);
for i = 1:rows
    signal = emgMatrix(i, :);
    [pxx, f] = pwelch(signal, window, noverlap, nfft, fs);
    averages_psd(i) = mean(pxx);
end

%| Variance PSD
variance_psd = zeros(rows, 1);
for i = 1:rows
    signal = emgMatrix(i, :);
    [pxx, f] = pwelch(signal, window, noverlap, nfft, fs);
    variance_psd(i) = var(pxx);
end

%| Kurtosis PSD
kurtosis_psd = zeros(rows, 1);
for i = 1:rows
    signal = emgMatrix(i, :);
    [pxx, f] = pwelch(signal, window, noverlap, nfft, fs);
    kurtosis_psd(i) = kurtosis(pxx);
end

%| Entropy PSD
entropies = zeros(rows, 1);
for i = 1:rows
    signal = emgMatrix(i, :);
    [pxx, f] = pwelch(signal, window, noverlap, nfft, fs);
    pxx_norm = pxx / sum(pxx); 
    pxx_norm(pxx_norm == 0) = eps;
    entropies(i) = -sum(pxx_norm .* log2(pxx_norm));
end


%| To display features for Boxplot and AUC ,we have to write the following lines of code ,and change only the name of required features

% Display Moment Coefficient of Kurtosis using your boxplot function
%boxplot_two_groups(moment_kurtosis_values);

% Compute AUC for Moment Coefficient of Kurtosis
compute_auc(moment_kurtosis_values);