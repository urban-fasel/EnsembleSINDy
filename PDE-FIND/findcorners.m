function [mx,mt,px,pt,sig_est,corners_all] = findcorners(U_obs,xs,tau,tauhat,max_dx,max_dt,opt)

mx = [];
mt = [];
L = length(xs{1});
T = length(xs{end});
corners_all = cell(length(U_obs),1);

if opt==1
    l = @(m,k,N) log((2*m-1)./m.^2).*(4*pi^2*k^2*m.^2-3*N^2*tauhat^2)-2*N^2*tauhat^2*log(tau);
elseif opt ==2
    px = [];
    pt = [];
end

for n= 1:length(U_obs)
    [corners,sig_est] = findcornerpts(U_obs{n},xs);
    for d = 1:length(xs)-1
        N = length(xs{d});
        k = corners{d}(2);
        if opt==2
            px = [px 2*pi*k/tauhat/N];
            mnew = 1+N*tauhat/2/pi/k*sqrt(-2*log(tau));
        elseif opt==1
            mstar1 = sqrt(3)/pi*N/2/k*tauhat;
            mstar2 = 1/pi*tauhat*(N/2)/k*sqrt(log(exp(1)^3/tau^8));
            mnew = fzero(@(m)l(m,k,N), [mstar1 mstar2]); 
            if mnew>N/2-1
                mnew = N/2/k;
            end
        end
        mx = [mx mnew];
        L = min(L,N);
    end
    k = corners{end}(2);
    if opt == 2
        mnew = 1+T*tauhat/2/pi/k*sqrt(-2*log(tau));
        pt = [pt 2*pi*corners{end}(2)/tauhat/length(xs{end})];
    elseif opt ==1
        mnew = fzero(@(m)l(m,k,T), [1 2/sqrt(tau)]); % provably a bracket for [1 (1+sqrt(1-tol))/tol] for all j < certain value
        if mnew>T/2-1
            mnew = T/2/k;
        end
    end
    mt = [mt mnew];
    corners_all{n}=corners;
end

mx = min(floor((L-1)/2),ceil(mean(mx)));
mt = min(floor((T-1)/2),ceil(mean(mt)));

if opt ==1
    px = max(max_dx+2,floor(log(tau)/log(1-(1-1/mx)^2)));
    pt = max(max_dt+2,floor(log(tau)/log(1-(1-1/mt)^2)));
elseif opt==2
    px = mean(px);
    pt = mean(pt);
end
end

function [corners,sig_est] = findcornerpts(U_obs,xs)
    dims = size(U_obs);
    dim = length(dims);
    corners = cell(dim,1);
    for d=1:dim
        if dim ==1
            shift = [1 2];
        else
            shift = circshift(1:dim,1-d);
        end
        dim_perm = dims(shift);
        x = xs{d}(:);
        L = length(x);
        wn = ((0:L-1)-floor(L/2))'*(2*pi)/range(x);
        xx = wn(1:ceil(end/2));
        NN = length(xx);
        Ufft = abs(fftshift(fft(permute(U_obs,shift))));       
        if dim>2
            Ufft = reshape(Ufft,[L,prod(dim_perm(2:end))]);
        end
        Ufft = mean(Ufft,2);
        Ufft = cumsum(Ufft);            
        Ufft = Ufft(1:ceil(L/2),:);
        errs = zeros(NN-6,1);
        for k=4:NN-3
           subinds1 = 1:k;
           subinds2 = k:NN;
           Ufft_av1 = Ufft(subinds1);
           Ufft_av2 = Ufft(subinds2);
           m1 = range(Ufft_av1)/range(xx(subinds1));
           m2 = range(Ufft_av2)/range(xx(subinds2));
           L1 = min(Ufft_av1)+m1*(xx(subinds1)-xx(1));
           L2 = max(Ufft_av2)+m2*(xx(subinds2)-xx(end));
           errs(k-3) = sqrt(sum(((L1-Ufft_av1)./Ufft_av1).^2) + sum(((L2-Ufft_av2)./Ufft_av2).^2)); % relative l2 
        end
        [~,tstarind] = min(errs);
        tstar = -xx(tstarind);
        corners{d} = [tstar NN-tstarind-3];
    end
    Ufft = abs(fftshift(fft(permute(U_obs,shift))));
    Ufft = Ufft([1:tstarind-3 2*NN-tstarind+2:end],:);
    sig_est = sqrt(mean(Ufft(:).^2)/L);
end