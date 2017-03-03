function [recoveredLightFieldResults] = cs_reconstruction(lightFieldImage, reconParams)

    %% create default parameters
    M = 2;
    reconBasis = ReconstructionBasis.FFT;
    

    %% read in parameters structure
    if(isfield(reconParams, 'numMeasurements'))
        % number of simulated measurements 
        assert(reconParams.numMeasurements > 1) % otherwise there is a bug with squeeze
        M = reconParams.numMeasurements;    
    elseif(isfield(reconParams, 'reconBasis'))
        reconBasis = reconParams.reconBasis;
    end   
    
    
    %% perform CS reconstruction    
    t_reconstruction = tic;

    % generate masks
    masks = rand(lightFieldImage.angularLightFieldSize, ...
        lightFieldImage.angularLightFieldSize, M);
    
    % create array for the recovered lightfield
    recoveredLightField = zeros(size(lightFieldImage.lightField));
        
    channels = 1:3;
    for c = channels
        % apply masks to the raw data to simulate measurements
        lightFieldImageSingleChannel = lightFieldImage.lightField(:, :, :, :, c);           
        Y = applyMasks(lightFieldImageSingleChannel, masks);

        % specify parameters for bpdn
        sigma = 0.1;            
        bpdnOptions = struct();
        bpdnOptions.iterations = 400;
        
        % solve reconstruction
        [x r g info] = spg_bpdn(@AReconFourierBasis, vectorizeLightField(Y), sigma, bpdnOptions);

        % reformat x into angular light fields
        recoveredLightFieldSingleChannel = reshape(x, [lightFieldImage.imageHeight, lightFieldImage.imageWidth, ...
            lightFieldImage.angularLightFieldSize, lightFieldImage.angularLightFieldSize]);
        
        % take inverse fft of recovered solution
        recoveredLightFieldSingleChannel = inverseBasisOperator(recoveredLightFieldSingleChannel);
        
        % save recovered solution into output
        recoveredLightField(:, :, :, :, c) = recoveredLightFieldSingleChannel;        
    end    

    % compute SNR   
    diff = lightFieldImage.lightField - recoveredLightField;
    mse = mean(diff(:).^2);
    SNR = 20 * log10(mean(lightFieldImage.lightField(:).^2)/mse);
    
    %% assemble results into structure
    recoveredLightFieldResults = struct();
    recoveredLightFieldResults.recoveredLightField = recoveredLightField;
    recoveredLightFieldResults.SNR = SNR;
    recoveredLightFieldResults.reconstructionTime = toc(t_reconstruction);
    recoveredLightFieldResults.numMeasurements = M;
    recoveredLightFieldResults.numAngularViews = lightFieldImage.angularLightFieldSize.^2;
    recoveredLightFieldResults.reconBasis = reconBasis;
    recoveredLightFieldResults.fractionOfMeasurements = M/(lightFieldImage.angularLightFieldSize.^2);
    
    return
    
    %%%%%%% Declare private functions %%%%%%%%%%%%%%%%    
    function y = forwardBasisOperator(lightFieldSingleChannel)        
        if(reconBasis == ReconstructionBasis.FFT)
            y = fft(fft(lightFieldSingleChannel, [],3), [], 4);
        else
            assert(false, 'invalid reconBasis');              
        end
    end

    function y = inverseBasisOperator(lightFieldSingleChannel)
        if(reconBasis == ReconstructionBasis.FFT)
            y = ifft(ifft(lightFieldSingleChannel, [], 3), [], 4);
        else            
            assert(false, 'invalid reconBasis');  
        end
    end
    

    %% function provided to BPDN solver that does reconstruction in Fourier basis
    function v  = AReconFourierBasis( w, mode )
    %ARECONFOURIERBASIS 
        imageWidth = lightFieldImage.imageWidth;
        imageHeight = lightFieldImage.imageHeight;

        if(mode == 1)
            % w is in the sparse basis and is the same size as the
            % lightfield we want to recover
            w1 = reshape(w, [imageHeight imageWidth ...
                lightFieldImage.angularLightFieldSize lightFieldImage.angularLightFieldSize]);
            
            lightFieldSpatialDomain = inverseBasisOperator(w1);

            
            % apply masks            
            y = applyMasks(lightFieldSpatialDomain, masks);
            
            % vectorize
            v = vectorizeLightField(y);
            
        elseif(mode == 2)            
            % the provided w is the same size as the measured data
            w2 = reshape(w, [imageHeight imageWidth M]);
            
            % apply masks adjoint operation, creates linear combination of
            % masks
            adjointOutput = zeros(imageHeight, imageWidth, ...
                lightFieldImage.angularLightFieldSize, lightFieldImage.angularLightFieldSize);

            for v = 1:lightFieldImage.angularLightFieldSize
               for u = 1:lightFieldImage.angularLightFieldSize                   
                   % sum linear combinations of the masks
                   for k = 1: M
                      adjointOutput(:, :, v, u) = adjointOutput(:, :, v, u) + ...
                          squeeze(masks(v, u, k)) * w2(:, :, k);
                   end
               end
            end                
            
            % apply sparse basis adjoint operation
            lightFieldFourierDomain = forwardBasisOperator(adjointOutput);
            
            % vectorize output
            v = lightFieldFourierDomain(:);

        end
    end
    
    %% input simulated light field data, size imageHeight x imageWidth x M
    % outputs vectorized data in images ordered
    function Y2 = vectorizeLightField(X)
        assert(numel(size(X)) == 3)
        imageWidth = lightFieldImage.imageWidth;
        imageHeight = lightFieldImage.imageHeight;
        numMeasurements = M;
        
        Y1 = reshape(X, [imageHeight*imageWidth numMeasurements]);
        Y2 = reshape(Y1, [numel(Y1), 1]);      
    end

    %% undo vectorizeLightField
    function Y2 = unvectorizeLightField(X)
        imageWidth = lightFieldImage.imageWidth;
        imageHeight = lightFieldImage.imageHeight;
        numMeasurements = M;
        assert(numel(X) == imageWidth * imageHeight * numMeasurements)
        Y1 = reshape(X,[imageHeight.*imageWidth numMeasurements]);
        Y2 = reshape(Y1,[imageHeight imageWidth numMeasurements]);
    end
    
    %% applies mask to a single channel light field 
    function Y = applyMasks(angularLightFields, measurement_masks)
        
        % create array for our measurements
        Y = zeros(lightFieldImage.imageHeight, ...
                lightFieldImage.imageWidth, M);

        for m = 1:M
            for u = 1: lightFieldImage.angularLightFieldSize        
                for v = 1: lightFieldImage.angularLightFieldSize
                    Y(:, :, m) = Y(:, :, m) + measurement_masks(v, u, m) .* ...
                        angularLightFields(:, :, v, u);  
                end
            end       
            
            % add noise to each measurement, disabled for now
            % stdev = 0.25;
            % Y(:, :, m) = Y(:, :, m) + stdev * randn(size(Y(:, :, m)));
            
        end
    end
end