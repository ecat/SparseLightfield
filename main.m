clear all; close all;

%% load bpdn solver
addpath('spgl1-1.9')


%% load image
parameters.filename = 'flowers_plants_1_eslf.png';
parameters.angularLightFieldSize = 8;
parameters.angularViewResizeFactor = 2;
parameters.brightnessScale = 4;

lightFieldImage = LightFieldImage(parameters)

%% perform cs reconstruction
[recoveredLightField] = cs_reconstruction(lightFieldImage);


%% compare results
u = 3; v = 3;


diff = lightFieldImage.lightField - recoveredLightField;
mse = sum(diff(:).^2)

figure;
subplot(2, 1, 1)
imshow(lightFieldImage.getImage(v, u));
title(['original image  v = ' num2str(v) ' u = ' num2str(u)])
subplot(2, 1, 2)
imshow(squeeze(recoveredLightField(:, :, v, u, :)));
title(sprintf('recovered image MSE = %.2f', mse))

