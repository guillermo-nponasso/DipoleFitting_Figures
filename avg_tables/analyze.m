data_table = readtable('combined_avg.csv');
data_table = data_table(:,2:end);
disp(data_table)

numeric_data = table2array(data_table);

std_vec = std(numeric_data,0,1);
normalization_matrix = (std_vec' * std_vec).^(-1);
covariance_matrix = cov(numeric_data);
fprintf("Error - tissue thickness correlation matrix:\n");
correlation_matrix = normalization_matrix.*covariance_matrix