function [ originalLightFieldImage reconstruction_results ] = lightfield_reconstruction( filename )
%LIGHTFIELD_RECONSTRUCTION Sweeps over different reconstruction parameteres
%for a given filename

    reconstruction_results = {};

    n = 1;
    %% load image
    parameters.filename = filename;
    parameters.angularLightFieldSize = 6;
    parameters.angularViewResizeFactor = 4;
    parameters.brightnessScale = 4;

    lightFieldImage = LightFieldImage(parameters);

    %% perform cs reconstruction over different parameters
    sweepTimer = tic;
    for numMeasurements = 2%[2 4 8 10 16]
        display(sprintf('Performing reconstruction using %d measurements.', numMeasurements))
        reconParams.numMeasurements = numMeasurements;
        reconParams.reconBasis = ReconstructionBasis.FFT;

        [recoveredLightFieldResults] = cs_reconstruction(lightFieldImage, reconParams);
        reconstruction_results{n} = recoveredLightFieldResults;
        
        n = n + 1;
    end
    display('Total time taken: ')
    toc(sweepTimer)

    originalLightFieldImage = lightFieldImage;
end

