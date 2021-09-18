%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WSINDy_PDE: compute ND convolution utilizing separability over 'valid' points, 
%%%%%%%%%%%% i.e. those that do not require zero-padding. Uses FFT to
%%%%%%%%%%%% compute each 1D convolution.
%%%%%%%%%%%%
%%%%%%%%%%%% Copyright 2020, All Rights Reserved
%%%%%%%%%%%% Code by Daniel A. Messenger
%%%%%%%%%%%% For Paper, "Weak SINDy for Partial Differential Equations"
%%%%%%%%%%%% by D. A. Messenger and D. M. Bortz

function X = convNDfft(X,cols,sub_inds,ver)
    Ns = size(X);
    dim = length(Ns);
    for k=1:dim
        if ver==1
            col = cols{k}(:);
            n = length(col);
%            col_ifft = Ns(k)*ifft([flipud(col);zeros(Ns(k)-n,1)]);
            col_ifft = fft([zeros(Ns(k)-n,1);col]);
        else
            col_ifft = cols{k}(:);
        end
        
        if dim ==1
            shift = [1 2];
            shift_back = shift;
        else
            shift = circshift(1:dim,1-k);
            shift_back=circshift(1:dim,k-1);
        end
        
        X = ifft(col_ifft.*fft(permute(X,shift)));
        inds = sub_inds(k);
        for j=2:dim
            inds{j} = ':';
        end
        X = X(inds{:});
        X = permute(X,shift_back);
                
    end
    X = real(X);
end