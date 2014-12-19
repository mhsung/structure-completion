close all;
clear all;

train_path = 'dataset_1/train/';
num_train_files = 300;

% for k = 1:num_train_files
%     disp(['File ID = ', num2str(k)]);
%     generate_dataset_1(train_path, k, false);
% end
% pause;

test_path = 'dataset_1/test/';
num_test_files = 300;

for k = 1:num_test_files
    disp(['File ID = ', num2str(k)]);
    generate_dataset_1(test_path, k, true);
end