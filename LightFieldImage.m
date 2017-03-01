classdef LightFieldImage
    %LIGHTFIELDIMAGE Manages accessing and loading light field
    %
    
    properties        
        lightField % 4D buffer, (height, width, v, u) 
        imageWidth, imageHeight % width and height of image        
        microlensSize % Number of pixels underneath a microlens (assumes squares)
        LF_microlens % raw buffer
        
        % default is 14 x 14, if we want to reduce the number of angular light field
        % views, we can specify a smaller number that copies only the
        % center angular views denoted by lightFieldSize
        angularLightFieldSize 
        % resize each of the angular views, 2x makes the image half the
        % size in each dimension
        angularViewResizeFactor
        
    end
    
    methods
        %% Constructor
        % parameters is a that contains:
        % filename
        % (optional)
        % angularLightFieldSize
        % angularViewResizeFactor
        % brightnessScale
        function obj = LightFieldImage(params)    
            
            filename = params.filename;
            
            obj.microlensSize = 14;
            
            % parse parameters
            if(isfield(params, 'angularLightFieldSize'))
                assert( mod(params.angularLightFieldSize, 2) == 0)
                obj.angularLightFieldSize = params.angularLightFieldSize;                
            else 
                obj.angularLightFieldSize = obj.microlensSize;                
            end
            
            if(isfield(params, 'angularViewResizeFactor'))
                obj.angularViewResizeFactor = params.angularViewResizeFactor;
            else
                obj.angularViewResizeFactor = 1; 
            end
            
            if(isfield(params, 'brightnessScale'))
                brightnessScale = params.brightnessScale;               
            else                
                brightnessScale = 1;
            end
                      
            % load lightField image
            obj.LF_microlens = im2double(imread(filename)) * brightnessScale;
            % size of the raw light field
            num_y = size(obj.LF_microlens, 1); 
            num_x = size(obj.LF_microlens, 2); 
            mImageWidth = num_x / obj.microlensSize;
            mImageHeight = num_y / obj.microlensSize;
            
            % arrange views into single images                             
            mLightField = zeros(mImageHeight, mImageWidth, ...
                obj.microlensSize, obj.microlensSize, 3);
            
            for ky = 1:obj.microlensSize
                for kx = 1:obj.microlensSize 
                    mLightField(:, : ,ky, kx, :) = ...
                        obj.LF_microlens(ky:obj.microlensSize:end, kx:obj.microlensSize:end, :);
                end 
            end
            
            % center crop to downsample the views            
            if(obj.angularLightFieldSize ~= obj.microlensSize || ...
                    obj.angularViewResizeFactor ~= 1)
                mImageHeight = round(mImageHeight/obj.angularViewResizeFactor);
                mImageWidth = round(mImageWidth/obj.angularViewResizeFactor);
                downsampledLightField = zeros(mImageHeight, mImageWidth, ...
                    obj.angularLightFieldSize, obj.angularLightFieldSize, ...
                    3);
                
                for ky = 1:obj.angularLightFieldSize
                    for kx = 1:obj.angularLightFieldSize
                        offsetX = (obj.microlensSize -  obj.angularLightFieldSize)/2;
                        offsetY = (obj.microlensSize -  obj.angularLightFieldSize)/2;                
                        originalImage = mLightField(:, :, ky + offsetY, kx + offsetX, :);

                        resizedImage = imresize(originalImage, [mImageHeight, mImageWidth], 'bilinear');
                        downsampledLightField(:, :, ky, kx, :) = ...
                            resizedImage;
                    end                    
                end
                
                
                mLightField = downsampledLightField;
            end
            
            obj.lightField = mLightField;
            
            obj.imageHeight = mImageHeight;
            obj.imageWidth = mImageWidth;
        end
        
        %% See lightfield from a certain perspective
        function outputImage = getImage(obj, v, u)
            if( u > obj.angularLightFieldSize || v > obj.angularLightFieldSize)
                error('u or v exceeds bounds')
            end
            
            outputImage = squeeze(obj.lightField(:, :, v, u, :));            
        end
        
        %% tiles the views in u v spatial domain
        function outputImage = getTiledImage(obj)
            width = obj.imageWidth;
            height = obj.imageHeight;
            
           outputImage = zeros(height * obj.angularLightFieldSize, ...
               width * obj.angularLightFieldSize, 3);
           
           for ky = 1:obj.angularLightFieldSize
              for kx = 1:obj.angularLightFieldSize
                    outputImage( (ky -1) * height + 1: ky * height, ...
                    (kx -1) * width + 1: kx * width, :) = obj.getImage(ky, kx);
              end
           end
            
        end
        
        %% returns the log magnitude FFT of the average of the channels
        %% of the raw lightfield image
        function outputImage = getFFTRawImage(obj)
            F = fft2(squeeze(mean(obj.LF_microlens(:, :, :), 3)));
            F = fftshift(F); % Center FFT
            F = abs(F); % Get the magnitude
            F = log(F+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
            F = mat2gray(F); % Use mat2gray to scale the image between 0 and 1            
            outputImage = F;
        end

    end
    
end

