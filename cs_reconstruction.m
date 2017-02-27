function cs_reconstruction(lightFieldImage)

    %% perform CS reconstruction

    % number of measurements
    M = 8;

    % generate masks
    masks = rand(lightFieldImage.angularLightFieldSize, ...
        lightFieldImage.angularLightFieldSize, M);
    % create array for our measurements
    Y = zeros(lightFieldImage.imageHeight, lightFieldImage.imageWidth , M);
  
    
    for c = 1:3
    
        % apply masks to the raw data to simulate measurements
        for m = 1:M
            lightFieldImageSingleChannel = lightFieldImage.lightField(:, :, :, :, c);
            lightFieldMeasurement = zeros(size(lightFieldImage.imageHeight, ...
                lightFieldImage.imageWidth));

            for u = 1: lightFieldImage.angularLightFieldSize        
                for v = 1: lightFieldImage.angularLightFieldSize
                    lightFieldMeasurement = lightFieldMeasurement + masks(v, u, m) .* ...
                        lightFieldImageSingleChannel(:, :, v, u);  
                end
            end
            
            Y(:, :, m) = lightFieldMeasurement;
        end
        

        % solve the reconstruction problem
        sigma = 0.1;    
        
        
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
    
    %B = measurements vectorized

    %spg_bpdn(AReconFourierBasis, B, sigma)  %% Change this        

% 
%     % calculate MSE between recovered image and original image
%     % lightFieldImage.lightField   %% Change this
%   
%     function v  = AReconFourierBasis( w, mode )
%     %ARECONFOURIERBASIS 
%         if mode ==1
%         return   
%         elseif mode == 2
%         return  v = sum(repmat(reshape(w, [1 1 M]), ...
%             [angularViewResizeFactor, angularViewResizeFactor, 1]).*masks, 3);
%         end
%         
%     
%         % w is in the sparse basis
%         if(mode == 1)
%             % apply ifft2 in the u and v axes
%             temp1 = repmat(ifft2(w), [1 1 M]);
%             % apply masks
%             v = squeeze(sum(sum(masks .* temp1,1),2));
%             
%         elseif(mode == 2)
%             % apply masks adjoint operation, creates linear combination of
%             % masks
% 
% 
%             % apply sparse basis adjoint operation
%             % fft2
% 
%         end
%     end


end