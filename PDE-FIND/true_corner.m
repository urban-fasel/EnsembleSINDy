function [tc,err] = true_corner(Ufft,sigma,c)
    dims = size(Ufft);
    dim = length(dims);
    sigmafft = sqrt(dims(1))*sigma;
    Ufft = abs(fftshift(fft(Ufft)));
    M= prod(dims(2:end));
    if dim>2
        Ufft = reshape(Ufft,[dims(1),M]);
    end
    inds = zeros(1,M);
    for j=1:M
        inds_temp = find(Ufft(floor(dims(1)/2):end,j)>sigmafft,1,'last');
        if isempty(inds_temp)
            inds(j) = 1;
        else
            inds(j) = inds_temp;
        end    
    end
    tc = mean(inds);
    err = (tc-c)/dims(1)*2;
end