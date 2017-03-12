% plotter for fft_actual_run
clear all; close all;

%% load data

folder_indices_with_results = [0:1];
num_images = numel(folder_indices_with_results);

for i = 1:num_images
    % loads data into a struct called reconstructionResults
    load(sprintf('run%d/reconstructionResults.mat', folder_indices_with_results(i)));
    
    % lazy init
    if(i == 1)
        measurements_swept = numel(reconstructionResults);
        
        SNR_array = zeros(num_images, measurements_swept);
        percent_num_measurements = zeros(measurements_swept, 1);
        num_angular_images = reconstructionResults{1}.numAngularViews;
       
        % read out number of measurements
        for k = 1: measurements_swept
            percent_num_measurements(k) = reconstructionResults{k}.fractionOfMeasurements;
        end
    end
    
    for k = 1: measurements_swept
        SNR_array(i, k) = real(reconstructionResults{k}.SNR);
    end
end

%% plot data
close all;

SNR_averages = mean(SNR_array, 1);
SNR_stdev = sqrt(var(SNR_array));
figure;
errorbar(percent_num_measurements, SNR_averages, SNR_stdev);
title('SNR Average for fft')
ylim([0 40])
grid on;
