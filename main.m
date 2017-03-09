%clear all; close all;

%% load bpdn solver
addpath('spgl1-1.9')

%% perform reconstruction over a lightfield image
[originalLF, reconstructionResults] = lightfield_reconstruction('flowers_plants_1_eslf.png')


%% compare results
u_index = 3; v_index = 3;

figure;
subplot(2, 1, 1)
imshow(originalLF.getImage(v_index, u_index));
title(['original image  v = ' num2str(v_index) ' u = ' num2str(u_index)])
subplot(2, 1, 2)
imshow(squeeze(reconstructionResults{1}.recoveredLightField(:, :, v_index, u_index, :)));
title(sprintf('recovered image SNR = %.2f dB compression = %.2f', reconstructionResults{1}.SNR, ...
    reconstructionResults{1}.fractionOfMeasurements))

figure;
subplot(2, 1, 1)
imshow(squeeze(reconstructionResults{1}.recoveredLightField(:, :, 1, 1, :)));
title(sprintf('View 1 1'))
subplot(2, 1, 2)
imshow(squeeze(reconstructionResults{1}.recoveredLightField(:, :, 6, 6, :)));
title(sprintf('View 6 6'))



