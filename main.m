clear all; close all;

%% load bpdn solver
addpath('spgl1-1.9')


%% load image
parameters.filename = 'flowers_plants_1_eslf.png';
parameters.angularLightFieldSize = 4;
%parameters.angularViewResizeFactor = 1;

lightFieldImage = LightFieldImage(parameters)

%imshow(lightFieldImage.getFFTRawImage())

%% perform cs reconstruction
cs_reconstruction(lightFieldImage);


