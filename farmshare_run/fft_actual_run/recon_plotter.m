% plotter for fft_actual_run
clear all; close all;

%% load data

folder_indices_with_results = [0:23];
num_images = numel(folder_indices_with_results);
runtimesSeconds = zeros(num_images, 1);
filenames = cell(num_images, 1);

for i = 1:num_images
    % loads data into a struct called reconstructionResults
    load(sprintf('run%d/reconstructionResults.mat', folder_indices_with_results(i)));
    filenames{i} = reconstructionResults{1}.filename;
   
    % lazy init
    if(i == 1)
        measurements_swept = numel(reconstructionResults);
        
        SNR_array = zeros(num_images, measurements_swept);
        percent_num_measurements = zeros(measurements_swept, 1);
        num_angular_images = reconstructionResults{1}.numAngularViews;
        display(sprintf('angular dimensions are %d %d', reconstructionResults{1}.angularImageWidth, ...
                    reconstructionResults{1}.angularImageHeight));
                
        % read out number of measurements
        for k = 1: measurements_swept
            percent_num_measurements(k) = reconstructionResults{k}.fractionOfMeasurements;
        end
        
    end
    
    for k = 1: measurements_swept
        SNR_array(i, k) = real(reconstructionResults{k}.SNR);
        
        runtimesSeconds(i) = runtimesSeconds(i) + reconstructionResults{k}.reconstructionTime;
    end
end

%% find index of worst and best
[worst_snr worst_index] = min(squeeze(SNR_array(:, 1)));
worst_filename = filenames{worst_index}

% print the filenames of those
[best_snr best_index] = max(squeeze(SNR_array(:, 1)));
best_filename = filenames{best_index}

%% plot data
close all;

SNR_averages = mean(SNR_array, 1);
SNR_stdev = sqrt(var(SNR_array));
figure;
errorbar(percent_num_measurements, SNR_averages, SNR_stdev);
title(['SNR Average for ' char(reconstructionResults{1}.reconBasis)])
ylim([0 40])
grid on;
