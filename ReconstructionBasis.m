classdef ReconstructionBasis
    %ReconstructionBasis enum to flip between different sparse basis in CS
    % reconstruction
    %   Detailed explanation goes here
    
    enumeration
        FFT
        HAAR
        DCT
        TV_PRIOR
        %WAVELET = {'HAAR','MORLET','SYMLET','MEXHAT'}   % 'haar','morl','sym','mexh'
    end
    
end

