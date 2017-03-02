clear all; close all;

%% load bpdn solver
addpath('spgl1-1.9')


%% load image
parameters.filename = 'flowers_plants_1_eslf.png';
parameters.angularLightFieldSize = 12;
parameters.angularViewResizeFactor = 2;
parameters.brightnessScale = 4;

lightFieldImage = LightFieldImage(parameters)

%% perform cs reconstruction
reconParams.numMeasurements = 16;
reconParams.reconBasis = ReconstructionBasis.FFT;
t_reconstruction = tic;
[recoveredLightField] = cs_reconstruction(lightFieldImage, reconParams);
toc(t_reconstruction)

%% compare results
u_index = 3; v_index = 3;

diff = lightFieldImage.lightField - recoveredLightField;
mse = mean(diff(:).^2)
SNR = 20 * log10(mean(lightFieldImage.lightField(:).^2)/mse);

figure;
subplot(2, 1, 1)
imshow(lightFieldImage.getImage(v_index, u_index));
title(['original image  v = ' num2str(v_index) ' u = ' num2str(u_index)])
subplot(2, 1, 2)
imshow(squeeze(recoveredLightField(:, :, v_index, u_index, :)));
title(sprintf('recovered image MSE = %.2f', mse))




