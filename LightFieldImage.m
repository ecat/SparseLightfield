classdef LightFieldImage
    %LIGHTFIELDIMAGE Manages accessing and loading light field
    %
    
    properties        
        lightField % 4D buffer, (height, width, v, u)        
        imageWidth, imageHeight % width and height of image        
        microlensSize % Number of pixels underneath a microlens (assumes squares)
    end
    
    methods
        %% Constructor
        function obj = LightFieldImage(filename)    
            obj.microlensSize = 14;
            
            % load lightField image
            LF_microlens = im2double(imread(filename)) .*5;
            % size of the raw light field
            num_y = size(LF_microlens, 1); 
            num_x = size(LF_microlens, 2); 
            obj.imageWidth = num_x / obj.microlensSize;
            obj.imageHeight = num_y / obj.microlensSize;
            
            % arrange views into single images                             
            obj.lightField = zeros(obj.imageHeight, obj.imageWidth, ...
                obj.microlensSize, obj.microlensSize, 3);
            for ky = 1:obj.microlensSize
                for kx = 1:obj.microlensSize 
                    obj.lightField(:, : ,ky, kx, :) = ...
                        LF_microlens(ky:obj.microlensSize:end, kx:obj.microlensSize:end, :);
                end 
            end
        end
        
        %% See lightfield from a certain perspective
        function outputImage = getImage(obj, u, v)
            if( u > obj.microlensSize || v > obj.microlensSize)
                error('u or v exceeds bounds')
            end
            
            outputImage = squeeze(obj.lightField(:, :, v, u, :));
            
        end
        
        function outputImage = getTiledImage(obj)
            width = obj.imageWidth;
            height = obj.imageHeight;
            
           outputImage = zeros(height * obj.microlensSize, ...
               width * obj.microlensSize, 3);
           
           for ky = 1:obj.microlensSize
              for kx = 1:obj.microlensSize
                    outputImage( (ky -1) * height + 1: ky * height, ...
                    (kx -1) * width + 1: kx * width, :) = obj.getImage(ky, kx);
              end
           end
            
        end
    end
    
end

