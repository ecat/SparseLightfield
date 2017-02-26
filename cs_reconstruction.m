function cs_reconstruction(lightFieldImage)

    %% perform CS reconstruction

    % number of measurements
    M = 8;

    % generate masks
    % uniform masks, size (angularViewResizeFactor x angularViewResizeFactor x M)
    masks = zeros(3, 3); %% Change this

    % apply masks to the raw data
    % Y = masks * raw data;  %% Change this

    % solve the reconstruction problem
    sigma = 0.1;
    % B = measurements vectorized

    % perform once for each colour channel
    for c = 1 :3
        %spg_bpdn(AReconFourierBasis, B, sigma)  %% Change this        
    end

    % calculate MSE between recovered image and original image
    % lightFieldImage.lightField   %% Change this

    function v  = AReconFourierBasis( w, mode )
    %ARECONFOURIERBASIS 
    %    returns  v = A *w  if mode == 1;
    %             v = A'*w  if mode == 2. 

        masks
    
        % w is in the sparse basis
        if(mode == 1)
            % apply ifft2 in the u and v axes

            % apply masks
        elseif(mode == 2)
            % apply masks adjoint operation, creates linear combination of
            % masks


            % apply sparse basis adjoint operation
            % fft2

        end
    end


end