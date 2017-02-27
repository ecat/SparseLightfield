function cs_reconstruction(lightFieldImage)

    %% perform CS reconstruction

    % number of measurements
    M = 8;

    % generate masks
    masks = rand(lightFieldImage.angularLightFieldSize, ...
        lightFieldImage.angularLightFieldSize, M);
    % create array for our measurements
    Y = zeros(lightFieldImage.imageHeight, lightFieldImage.imageWidth , M);
      
    %% applies mask to a light field single channel
    function Y = applyMasks(angularLightFields, measurement_masks)
        % add some asserts
        Y = zeros(lightFieldImage.imageHeight, ...
                lightFieldImage.imageWidth, M);

        for m = 1:M
            for u = 1: lightFieldImage.angularLightFieldSize        
                for v = 1: lightFieldImage.angularLightFieldSize
                    Y(:, :, m) = Y(:, :, m) + measurement_masks(v, u, m) .* ...
                        angularLightFields(:, :, v, u);  
                end
            end        
        end
    end
    
    for c = 1
        % apply masks to the raw data to simulate measurements
        lightFieldImageSingleChannel = lightFieldImage.lightField(:, :, :, :, c);           
        Y = applyMasks(lightFieldImageSingleChannel, masks);

        % solve the reconstruction problem
        sigma = 0.1;    
        
        spg_bpdn(AReconFourierBasis, vectorizeLightField(Y), sigma)

        % calculate MSE between recovered image and original image
        % lightFieldImage.lightField   %% Change this
    end

    %% function provided to BPDN solver
    function v  = AReconFourierBasis( w, mode )
    %ARECONFOURIERBASIS 
        imageWidth = lightFieldImage.imageWidth;
        imageHeight = lightFieldImage.imageHeight;

        % w is in the sparse basis
        if(mode == 1)
            w1 = reshape(w, [imageHeight imageWidth ...
                lightFieldImage.angularLightFieldSize lightFieldImage.angularLightFieldSize]);
            
            % apply ifft in the u and v axes
            lightFieldSpatialDomain = ifft(ifft(w1, [], 3), [], 4);
            
            % apply masks            
            y = applyMasks(lightFieldSpatialDomain, masks);
            
            % vectorize
            v = vectorizeLightField(y);
            
        elseif(mode == 2)
            % apply masks adjoint operation, creates linear combination of
            % masks
            w2 = reshape(w, [imageHeight imageWidth ...
                lightFieldImage.angularLightFieldSize lightFieldImage.angularLightFieldSize]);
            aa = reshape(w2(lightFieldImage.angularLightFieldSize*lightFieldImage.angularLightFieldSize, M))
            % apply sparse basis adjoint operation
            % fft2
            lightFieldFourierDomain = fft(fft(w2, [],3), [], 4);
            %v = sum(repmat(reshape(w, [1 1 M]), ...
            
            %    [angularViewResizeFactor, angularViewResizeFactor, 1]).*masks, 3);
        end
    end
    
    %% input measured light field data, size imageHeight x imageWidth x M
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
    

end