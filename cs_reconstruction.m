function [recoveredLightFieldResults] = cs_reconstruction(lightFieldImage, reconParams)

    %% create default parameters
    M = 2;
    reconBasis = ReconstructionBasis.FFT;
    
    %% read in parameters structure
    if(isfield(reconParams, 'numMeasurements'))
        % number of simulated measurements 
        assert(reconParams.numMeasurements > 1) % otherwise there is a bug with squeeze
        M = reconParams.numMeasurements;    
    end
    
    if(isfield(reconParams, 'reconBasis'))
        reconBasis = reconParams.reconBasis;
    end   
    
    display(sprintf('Performing reconstruction for %d measurements and basis', M) );
    display(reconBasis);
    
    %% perform CS reconstruction    
    t_reconstruction = tic;

    % generate masks
    masks = rand(lightFieldImage.angularLightFieldSize, ...
        lightFieldImage.angularLightFieldSize, M);
    
    % create array for the recovered lightfield
    recoveredLightField = zeros(size(lightFieldImage.lightField));
        
    channels = [1:3];
    for c = channels
        display(['START M = ' num2str(M) ' channel ' num2str(c)])
        % apply masks to the raw data to simulate measurements
        lightFieldImageSingleChannel = lightFieldImage.lightField(:, :, :, :, c);           
        Y = applyMasks(lightFieldImageSingleChannel, masks);

        % specify parameters for bpdn
        bpOptions = struct();
        bpOptions.iterations = 1500;
        bpOptions.verbosity = 1;
        
        % solve reconstruction
        sigma = 0.0;
        
        if(reconBasis == ReconstructionBasis.FFT || reconBasis == ReconstructionBasis.DCT)
            [x_rec r g info] = spg_bpdn(@AReconFourierBasis, vectorizeLightField(Y), sigma, bpOptions);
        elseif(reconBasis == ReconstructionBasis.HAAR)
            [x_rec r g info] = spg_bpdn(@AReconWaveletBasis, vectorizeLightField(Y), sigma, bpOptions);
        elseif(reconBasis == ReconstructionBasis.TV_PRIOR)
            lambda = 0.1;
            [recoveredLightFieldSingleChannel] = reconSparseGradients(Y, lambda);
        else
            assert(false)
        end
        
        if(reconBasis == ReconstructionBasis.FFT || reconBasis == ReconstructionBasis.DCT)
            % reformat x into angular light fields
            recoveredLightFieldSingleChannel = reshape(x_rec, [lightFieldImage.imageHeight, lightFieldImage.imageWidth, ...
                lightFieldImage.angularLightFieldSize, lightFieldImage.angularLightFieldSize]);
            
            % take inverse fft of recovered solution
            recoveredLightFieldSingleChannel = inverseBasisOperator(recoveredLightFieldSingleChannel);
        elseif(reconBasis == ReconstructionBasis.HAAR)
            recoveredLightFieldSingleChannelSparse = reshape(x_rec, [lightFieldImage.imageHeight, lightFieldImage.imageWidth, ...
                lightFieldImage.angularLightFieldSize/2, lightFieldImage.angularLightFieldSize/2, 4]);           
            recoveredLightFieldSingleChannel = backwardWaveletOperator(recoveredLightFieldSingleChannelSparse);
        elseif(reconBasis == ReconstructionBasis.TV_PRIOR)
            % do nothing
        else
            assert(false)
        end
        
        % save recovered solution into output
        recoveredLightField(:, :, :, :, c) = recoveredLightFieldSingleChannel;        
    end    

    % compute SNR   
    diff = lightFieldImage.lightField - recoveredLightField;
    mse = mean(diff(:).^2);
    SNR = 10 * log10(mean(lightFieldImage.lightField(:).^2)/mse);
    
    %% assemble results into structure
    recoveredLightFieldResults = struct();
    recoveredLightFieldResults.recoveredLightField = recoveredLightField;
    recoveredLightFieldResults.SNR = SNR;
    recoveredLightFieldResults.reconstructionTime = toc(t_reconstruction);
    recoveredLightFieldResults.numMeasurements = M;
    recoveredLightFieldResults.numAngularViews = lightFieldImage.angularLightFieldSize.^2;
    recoveredLightFieldResults.reconBasis = reconBasis;
    recoveredLightFieldResults.fractionOfMeasurements = M/(lightFieldImage.angularLightFieldSize.^2);
    recoveredLightFieldResults.angularImageWidth = lightFieldImage.imageWidth;
    recoveredLightFieldResults.angularImageHeight = lightFieldImage.imageHeight;
    
    return
    
    %%%%%%% Declare private functions %%%%%%%%%%%%%%%%    
    function y = reconSparseGradients(measuredLightField, lambda)        
        rho = 10;
        kappa = lambda/rho;
        option_fs.filterstrength = sqrt(kappa);
        size_y = lightFieldImage.imageHeight;
        size_x = lightFieldImage.imageWidth;
        size_v = lightFieldImage.angularLightFieldSize;
        size_u = lightFieldImage.angularLightFieldSize;
        imageResolution = [size_v size_u];
        
        % define function handle for image formation
        % dimension: size_y*size_x to M*1
        Afun    = @(x) squeeze(sum(sum(masks .* repmat(x, [1 1 M]),1),2) );

        % dimension: M*1 to size_y*size_x
        Atfun   = @(x) sum(repmat(reshape(x, [1 1 M]), [size_v size_u 1]).*masks, 3);

        % dimension: size_y*size_x to size_y*size_x        
        AtAfun = @(x) Atfun(Afun(reshape(x,imageResolution)));

        % dimension: size_y*size_x to size_y*size_x        
        optDtDx = @(x) opDtx(opDx(reshape(x,imageResolution)));

        
        %initialize x, z, u;
        x = zeros(imageResolution);
        z = zeros(size_v,size_u,2);
        u = zeros(size_v,size_u,2);
        y = zeros([size_y, size_x, size_v, size_u]);
        
        
        for dy = 1:size_y
            for dx = 1:size_x
                b = squeeze(measuredLightField(dy,dx,:,:));
                for iterations = 1:25
                    % update x
                    v = z-u;
                    v1 = v(:,:,1);
                    v2 = v(:,:,2);
                    imageSize = size_v * size_u;
                    pcg_term1 = @(x) reshape(AtAfun(x) + rho.*optDtDx(x),[imageSize 1]);
                    
                    pcg_term2 = @(v) reshape(Atfun(b) + rho.*opDtx(v), [imageSize 1]);
                    
                    [x, trash] = pcg(pcg_term1,pcg_term2(v),1e-10,50); % trash is to suppress output
                    x = reshape(x,imageResolution);
                    
                    % update z
                    v = opDx(x) + u;
                    
                    first_case = v > kappa;
                    z1 = first_case.*(v-kappa);
                    second_case = abs(v) <= kappa;
                    z2 = second_case.*0;
                    third_case = v < -kappa;
                    z3 = third_case.*(v+kappa);
                    z = z1+z2+z3;
                    
                    % update u
                    u = u + opDx(x) -z;
                    
                end
                y(dy,dx,:,:) = x;
            end
        end
    end

    
    function y = forwardWaveletOperator(lightFieldSingleChannel)
        imageWidth = lightFieldImage.imageWidth;
        imageHeight = lightFieldImage.imageHeight;
        
        lightFieldWaveletDomain = zeros(imageHeight, imageWidth, ...
            lightFieldImage.angularLightFieldSize/2, lightFieldImage.angularLightFieldSize/2, 4);
        
        % apply dwt2 to data
        for dy = 1:imageHeight
            for dx = 1:imageWidth
                [cA, cH, cV, cD] = dwt2(squeeze(lightFieldSingleChannel(dy, dx, :, :)), 'haar');
                lightFieldWaveletDomain(dy, dx, :, :, 1) = cA;
                lightFieldWaveletDomain(dy, dx, :, :, 2) = cH;
                lightFieldWaveletDomain(dy, dx, :, :, 3) = cV;
                lightFieldWaveletDomain(dy, dx, :, :, 4) = cD;
            end
        end
        
        y = lightFieldWaveletDomain;
    end
    
    function y = backwardWaveletOperator(lightFieldSingleChannel)  
        imageWidth = lightFieldImage.imageWidth;
        imageHeight = lightFieldImage.imageHeight;   
        
        assert(numel(size(lightFieldSingleChannel)) == 5);
        
        cA = squeeze(lightFieldSingleChannel(:, :, :, :, 1));
        cH = squeeze(lightFieldSingleChannel(:, :, :, :, 2));
        cV = squeeze(lightFieldSingleChannel(:, :, :, :, 3));
        cD = squeeze(lightFieldSingleChannel(:, :, :, :, 4));
        
        lightFieldSpatialDomain = zeros(size(squeeze(lightFieldImage.lightField(:, :, :, :, 1))));
        
        for dy = 1: imageHeight
            for dx = 1:imageWidth
                cA_yx = squeeze(cA(dy, dx, :, :));
                cH_yx = squeeze(cH(dy, dx, :, :));
                cV_yx = squeeze(cV(dy, dx, :, :));
                cD_yx = squeeze(cD(dy, dx, :, :));
                
                lightFieldSpatialDomain(dy, dx, :, :) = ...
                    idwt2(cA_yx, cH_yx, cV_yx, cD_yx, 'haar');
            end
        end
        
        y = lightFieldSpatialDomain;
    end
    
    
    function y = forwardBasisOperator(lightFieldSingleChannel)  
        if(reconBasis == ReconstructionBasis.FFT)
            y = fft(fft(lightFieldSingleChannel, [], 3), [], 4);
        elseif(reconBasis == ReconstructionBasis.DCT)
            
            y = zeros(size(lightFieldImageSingleChannel));
            imageWidth = lightFieldImage.imageWidth;
            imageHeight = lightFieldImage.imageHeight;
            
            temp = permute(lightFieldSingleChannel,[3 4 1 2]);
            for dy = 1:imageHeight
                for dx = 1:imageWidth
                    y(dy,dx,:,:) = dct2(squeeze(temp(:,:,dy,dx)));
                end
            end
        else
            assert(false, 'invalid reconBasis');
        end
    end

    function y = inverseBasisOperator(lightFieldSingleChannel)

        if(reconBasis == ReconstructionBasis.FFT)
            y = ifft(ifft(lightFieldSingleChannel, [], 3), [], 4);
        elseif(reconBasis == ReconstructionBasis.DCT)
            imageWidth = lightFieldImage.imageWidth;
            imageHeight = lightFieldImage.imageHeight; 
            y = zeros(size(lightFieldImageSingleChannel));            
            
            temp = permute(lightFieldSingleChannel,[3 4 1 2]);
            for dy =  1:imageHeight
                for dx = 1:imageWidth
                    y(dy,dx,:,:) = idct2(squeeze(temp(:,:,dy,dx))); 
                end
            end
        else
            assert(false, 'invalid reconBasis');  
        end
    end
    
    function v = AReconWaveletBasis(w, mode)
        imageWidth = lightFieldImage.imageWidth;
        imageHeight = lightFieldImage.imageHeight;
        
        if(mode == 1)
            % w is in the sparse basis
            w1 = reshape(w, imageHeight, imageWidth, ...
                lightFieldImage.angularLightFieldSize/2, lightFieldImage.angularLightFieldSize/2, 4);
            
            lightFieldSpatialDomain = backwardWaveletOperator(w1);
            
            % apply masks
            y = applyMasks(lightFieldSpatialDomain, masks);
            
            % vectorize
            v = vectorizeLightField(y);
        elseif(mode == 2)
            % the provided w is the same size as the measured data
            w2 = reshape(w, [imageHeight imageWidth M]);
            
            % apply masks adjoint operation, creates linear combination of
            % masks                   
            w3 = repmat(w2, [1, 1, 1, ...
                lightFieldImage.angularLightFieldSize, lightFieldImage.angularLightFieldSize]);
            
            w4 = permute(w3, [1 2 4 5 3]);
            
            masks_repeated = repmat(masks, [1, 1, 1, lightFieldImage.imageHeight, lightFieldImage.imageWidth]);
            masks_repeated = permute(masks_repeated, [4 5 1 2 3]);
            
            adjointOutput = sum(w4 .* masks_repeated, 5);
            
            lightFieldWaveletDomain = forwardWaveletOperator(adjointOutput);
                        
            v = lightFieldWaveletDomain(:);
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
            w3 = repmat(w2, [1, 1, 1, ...
                lightFieldImage.angularLightFieldSize, lightFieldImage.angularLightFieldSize]);
            
            w4 = permute(w3, [1 2 4 5 3]);
            
            masks_repeated = repmat(masks, [1, 1, 1, lightFieldImage.imageHeight, lightFieldImage.imageWidth]);
            masks_repeated = permute(masks_repeated, [4 5 1 2 3]);
            
            adjointOutput = sum(w4 .* masks_repeated, 5);
            
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