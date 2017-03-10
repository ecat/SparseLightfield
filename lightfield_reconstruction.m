function [ originalLightFieldImage, reconstruction_results ] = lightfield_reconstruction( filename )
%LIGHTFIELD_RECONSTRUCTION Sweeps over different reconstruction parameteres
%for a given filename

    %% load image
    parameters.filename = filename;
    parameters.angularLightFieldSize = 10;
    parameters.angularViewResizeFactor = 8;
    parameters.brightnessScale = 4;

    lightFieldImage = LightFieldImage(parameters);

    %% perform cs reconstruction over different parameters
    sweepTimer = tic;

    numMeasurements = 8;%[2 4 8 12 16 22 24 26 28 32 36 40];

    
    reconstruction_results = cell(numel(numMeasurements), 1);
    
    % on corn default is 12 threads
    parfor k = 1: numel(numMeasurements);
    %for k = 1: numel(numMeasurements)
        %display(sprintf('Performing reconstruction using %d measurements.', numMeasurements))
        reconParams = struct();
        reconParams.numMeasurements = numMeasurements(k);
        reconParams.reconBasis = ReconstructionBasis.FFT;
        %reconParams.reconBasis = ReconstructionBasis.HAAR;
        %reconParams.reconBasis = ReconstructionBasis.DCT;
        reconParams.reconBasis = ReconstructionBasis.TV_PRIOR;


        [recoveredLightFieldResults] = cs_reconstruction(lightFieldImage, reconParams);
        recoveredLightFieldResults.filename = filename;
        reconstruction_results{k} = recoveredLightFieldResults;
    end
    

    
    display('Total time taken: ')
    toc(sweepTimer)

    originalLightFieldImage = lightFieldImage;
end

