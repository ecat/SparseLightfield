function [ originalLightFieldImage reconstruction_results ] = lightfield_reconstruction( filename )
%LIGHTFIELD_RECONSTRUCTION Sweeps over different reconstruction parameteres
%for a given filename

    %% load image
    parameters.filename = filename;
    parameters.angularLightFieldSize = 10;
    parameters.angularViewResizeFactor = 4;
    parameters.brightnessScale = 4;

    lightFieldImage = LightFieldImage(parameters);

    %% perform cs reconstruction over different parameters
    sweepTimer = tic;
    numMeasurements = [2 4 8 16];
    
    reconstruction_results = cell(numel(numMeasurements), 1);
    poolobj = parpool(4);
    
    parfor k = 1: numel(numMeasurements);
        display(sprintf('Performing reconstruction using %d measurements.', numMeasurements))
        reconParams = struct();
        reconParams.numMeasurements = numMeasurements;
        reconParams.reconBasis = ReconstructionBasis.FFT;

        [recoveredLightFieldResults] = cs_reconstruction(lightFieldImage, reconParams);
        reconstruction_results{k} = recoveredLightFieldResults;
    end
    
    delete(poolobj)
    
    display('Total time taken: ')
    toc(sweepTimer)

    originalLightFieldImage = lightFieldImage;
end

