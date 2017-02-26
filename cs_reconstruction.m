function cs_reconstruction(lightFieldImage)

    %% perform CS reconstruction

    % number of measurements
    M = 8;

    % generate masks
       masks = rand(angularViewResizeFactor,angularViewResizeFactor,M);
  
%     % uniform masks, size (angularViewResizeFactor x angularViewResizeFactor x M)
%     masks = zeros(3, 3); %% Change this

    % apply masks to the raw data
      maskedImage = masks .* lightFieldImage ;  

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
        if mode ==1
        return   
        elseif mode == 2
        return  v = sum(repmat(reshape(w, [1 1 M]), ...
            [angularViewResizeFactor, angularViewResizeFactor, 1]).*masks, 3);
        end
        
    
        % w is in the sparse basis
        if(mode == 1)
            % apply ifft2 in the u and v axes
            temp1 = repmat(ifft2(w), [1 1 M]);
            % apply masks
            v = squeeze(sum(sum(masks .* temp1,1),2));
            
        elseif(mode == 2)
            % apply masks adjoint operation, creates linear combination of
            % masks


            % apply sparse basis adjoint operation
            % fft2

        end
    end


end