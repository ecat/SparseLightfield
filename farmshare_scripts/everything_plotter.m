% plotter for fft_actual_run
%clear all; close all;

%% load data
experiments = {'dct_actual_run', 'fft_actual_run', 'tv_actual_lambda_0point1'};
% create structo to hold all results
all_results = containers.Map();

for k = 1:numel(experiments)
    experiment = experiments{k}
    
    folder_indices_with_results = [0:23];
    num_images = numel(folder_indices_with_results);
    runtimesSeconds = zeros(num_images, 1);
    filenames = cell(num_images, 1);

    for i = 1:num_images
        % loads data into a struct called reconstructionResults
        load(sprintf('farmshare_run/%s/run%d/reconstructionResults.mat', experiment, folder_indices_with_results(i)));
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
    worst_filename = filenames{worst_index};

    % print the filenames of those
    [best_snr best_index] = max(squeeze(SNR_array(:, 1)));
    best_filename = filenames{best_index};
    
    
    % compile all results for this basis into a struct
    single_basis_results = struct();
    single_basis_results.worst_filename = worst_filename;
    single_basis_results.best_filename = best_filename;
    single_basis_results.SNR_array = SNR_array;
    single_basis_results.percent_num_measurements = percent_num_measurements;
    single_basis_results.num_angular_images = num_angular_images;
    single_basis_results.worst_index = worst_index;
    single_basis_results.best_index = best_index;
    
    all_results(experiment) = single_basis_results;
    
end

%% plot SNR averages data
figure;

for k = 1:numel(experiments)
    experiment = experiments{k};
    
    case_to_plot = all_results(experiment);
    case_to_plot.best_filename
    case_to_plot.worst_filename
    
    SNR_averages = mean(case_to_plot.SNR_array, 1);
    SNR_stdev = sqrt(var(case_to_plot.SNR_array));
    x_values = case_to_plot.percent_num_measurements;

    errorbar(x_values, SNR_averages, SNR_stdev);
    hold on;
    title('SNR Average')
    ylim([0 40])
    grid on;
end

legend(experiments);

%% plot composite images figure, first plot the worst image at index 22
% best image is at index 6
% second worse image is at index 12
for image_index = [6 12 22]
    figure;
    run_index = image_index - 1;
    u_index = 5;
    v_index = 5;
    for basis_index = 1:numel(experiments)
        experiment = experiments{basis_index};
        case_to_plot = all_results(experiment);

        % load the experiment
        load(sprintf('farmshare_run/%s/run%d/reconstructionResults.mat', experiment, run_index));

        for compression_index = 1:4
            subplot_index = compression_index + 4 * (basis_index-1);
            subplot(numel(experiments), 4, subplot_index)
            imshow(squeeze( ...
                real(reconstructionResults{compression_index}.recoveredLightField(:, :, v_index, u_index, :))));

        end
    end
end

%% show sparsity in angular domain
%% load image
parameters.filename = 'occlusions_15_eslf.png';
parameters.angularLightFieldSize = 10;
parameters.angularViewResizeFactor = 6;
parameters.brightnessScale = 4;

lightFieldImage = LightFieldImage(parameters);

% put x indices in normalized coordinates
x_indices = [0.5:0.02:0.6];
y_indices = [0.5];
for image_index = [13]
    figure;
    run_index = image_index - 1;

    for basis_index = 1:numel(experiments) + 1

        % load the experiment
        if(basis_index > numel(experiments))
            recoveredLightField = lightFieldImage.lightField;
        else
            experiment = experiments{basis_index};
            load(sprintf('farmshare_run/%s/run%d/reconstructionResults.mat', experiment, run_index));
            compression_index = 1;
            display(reconstructionResults{1}.reconBasis);
            display(reconstructionResults{1}.SNR);
            
            recoveredLightField = reconstructionResults{compression_index}.recoveredLightField;
        end

        for k = 1:numel(x_indices)
            lfSize = size(recoveredLightField);
            image_width = lfSize(2);
            image_height = lfSize(1);
            x_index = round(x_indices(k) * image_width);
            y_index = round(y_indices(k) * image_height);
            
            subplot_index = k + numel(x_indices) * (basis_index-1);
            subplot(numel(experiments) + 1, numel(x_indices), subplot_index)
            imshow(squeeze( ...
                real(recoveredLightField(x_index, y_index, :, :, :))));

        end
    end
end