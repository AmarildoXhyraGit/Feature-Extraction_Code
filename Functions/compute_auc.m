function auc = compute_auc(features)
    % Compute AUC for two groups of features classified as two classes.
    % Inputs:
    %   features - A 20xN matrix where each row is a signal (N is the number of samples per signal)
    % Outputs:
    %   auc - Area Under the Curve (AUC) for the classification task based on the features
    
    % Validate input matrix has exactly 20 rows
    if size(features, 1) ~= 20
        error('Input matrix must have exactly 20 rows.');
    end
    
    % Split the matrix into two groups of 10 rows each
    feature_group1 = features(1:10, :);  % First 10 rows (Class 1)
    feature_group2 = features(11:20, :); % Last 10 rows (Class 2)
    
    
    % Combine the features and create labels for the classes
    features = [feature_group1; feature_group2];  % 20 feature values (10 for each group)
    labels = [ones(10, 1); zeros(10, 1)];         % Class 1 for group1, Class 0 for group2
    
    % Use the perfcurve function to compute ROC and AUC
    [X, Y, T, auc] = perfcurve(labels, features, 0); % AUC for class 0 vs class 1
    specificity = 1 - X; % Specificity is 1 - FPR

    % CHOOSE BEST THRESHOLD USING YOUDEN'S J STATISTIC
    % Youden's J statistic
    J = Y + specificity - 1;
    
    % Find the threshold that maximizes Youden's J statistic
    [~, optimalIndex] = max(J);
    optimalThreshold = T(optimalIndex);
    
    disp(['Optimal Threshold (Max Youden''s J): ', num2str(optimalThreshold)]);

    
    % Step 1: Apply the optimal threshold to get predicted labels
    predicted_labels = features < optimalThreshold;
    
    % Step 2: Calculate confusion matrix components
    TP = sum((predicted_labels == 0) & (labels == 0)); % True Positives
    TN = sum((predicted_labels == 1) & (labels == 1)); % True Negatives
    FP = sum((predicted_labels == 0) & (labels == 1)); % False Positives
    FN = sum((predicted_labels == 1) & (labels == 0)); % False Negatives
    
    % Step 3: Calculate Accuracy, Precision, Recall, and F1 Score
    accuracy = (TP + TN) / (TP + TN + FP + FN); 
    precision = TP / (TP + FP); 
    recall = TP / (TP + FN); 
    f1_score = 2 * (precision * recall) / (precision + recall);

    % Display the results
    disp(['Accuracy: ', num2str(accuracy)]);
    disp(['Precision: ', num2str(precision)]);
    disp(['Recall (Sensitivity): ', num2str(recall)]);
    disp(['F1 Score: ', num2str(f1_score)]);


    
    % Plot the ROC curve
    figure;
    plot(X, Y);
    xlabel('False positive rate');
    ylabel('True positive rate');
    title(['ROC Curve (AUC = ' num2str(auc) ')']);

end
