addpath('spgl1-1.9')

parameters.filename = 'flowers_plants_1_eslf.png';
parameters.downsampledAngularLightFieldSize = 4;

lightFieldImage = LightFieldImage(parameters)

imshow(lightFieldImage.getFFTRawImage())


