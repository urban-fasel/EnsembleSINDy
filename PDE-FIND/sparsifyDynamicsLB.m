function [Xi,its,thrs_EL,XiDBs] = sparsifyDynamicsLB(Theta_0,dXdt,lambda,gamma,M_0,LBp)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%
% modified by Daniel A. Messenger, 2020 to prevent return of zero vector
% and include regularization
%
% compute Sparse regression: sequential least squares
%
% library bagging ensemble SINDy

nE1 = LBp.nE1;
nE2 = LBp.nE2;
ensT = LBp.ensT;
nEnsemblesDD = LBp.nE3;
ensembleT = LBp.ensT2;

n = min(size(dXdt));

nEnsemble1 = round(nE1*size(Theta_0,2));
mOutBS = zeros(nEnsemble1,n,nE2);
libOutBS = zeros(nEnsemble1,nE2);

% nn = size(Theta,2);
nn = size(Theta_0,2);
% nn = nEnsemble1;

for iii = 1:nE2
    rs = RandStream('mlfg6331_64','Seed',iii); 
    libOutBS(:,iii) = datasample(rs,1:size(Theta_0,2),nEnsemble1,'Replace',false)';

    Theta = Theta_0(:,libOutBS(:,iii));
    M = M_0(libOutBS(:,iii));
    
    if  gamma ~= 0
        Theta = [Theta;gamma*eye(nn)];
        dXdt = [dXdt;zeros(nn,n)];
    end

    Xi = Theta \ dXdt;  % initial guess: Least-squares
    if ~isempty(M)
        Xi = M.*Xi;
    end

    if isempty(M)
        thrs_EL = lambda;
    else
        bnds = norm(dXdt)./vecnorm(Theta)'.*M; 
        LBs = lambda*max(1,bnds);
        UBs = 1/lambda*min(1,bnds);
        thrs_EL = [LBs bnds UBs];
    end

    smallinds = 0*Xi;
    for j=1:nn
        if ~isempty(M)
            smallinds_new = or(abs(Xi)<LBs,abs(Xi)>UBs);
            if all(smallinds_new(:)==smallinds(:))
                its = j;
%                 return
                continue
            else
                smallinds = smallinds_new;
                Xi(smallinds)=0;    
                for ind=1:n
                    Xi(~smallinds,ind) = M(~smallinds).*(Theta(:,~smallinds)\dXdt(:,ind));
                end
            end
        else
            smallinds_new = (abs(Xi)<lambda);
            if all(smallinds_new(:)==smallinds(:))
                its = j;
%                 return
                continue
            else
                smallinds = smallinds_new;
                Xi(smallinds)=0;
                for ind = 1:n        
                    biginds = ~smallinds(:,ind);
                    Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
                end
            end
        end
    end
    its = nn;

    mOutBS(:,:,iii) = Xi;

end


inclProbBS = zeros(size(Theta_0,2),n);
for iii = 1:nE2
    for jjj = 1:n
        for kkk = 1:nEnsemble1
            if mOutBS(kkk,jjj,iii) ~= 0
                inclProbBS(libOutBS(kkk,iii),jjj) = inclProbBS(libOutBS(kkk,iii),jjj) + 1;
            end
        end
    end
end
inclProbBS = inclProbBS/nE2*size(Theta_0,2)/nEnsemble1;

Xi = zeros(size(Theta_0,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;
%     XiBias = sparsifyDynamics(Theta_0(:,libEntry),dXdt(:,iii),lambda,gamma,M);
    [XiBias,its,thrs_EL] = sparsifyDynamics(Theta_0(:,libEntry),dXdt(:,iii),lambda,gamma,M_0(libEntry));
    Xi(libEntry,iii) = XiBias;
end

if LBp.DB
    XiDB = zeros(size(Theta_0,2),n);
    XiDBmed = zeros(size(Theta_0,2),n);
    XiDBs = zeros(size(Theta_0,2),n);
    for iii = 1:n
        libEntry = inclProbBS(:,iii)>ensT;
        
        rng(iii,'twister')
        
        MIn = M_0(libEntry);
        bootstatDD = bootstrp(nEnsemblesDD,@(ThetaIn,dxIn)sparsifyDynamics(ThetaIn,dxIn,lambda,gamma,MIn),Theta_0(:,libEntry),dXdt(:,iii)); 

        XiDBe = [];
        XiDBnz = [];
        for iE = 1:nEnsemblesDD
            XiDBe(:,iE) = reshape(bootstatDD(iE,:),size(Theta_0(:,libEntry),2),1);
            XiDBnz(:,iE) = XiDBe(:,iE)~=0;
        end

        % Thresholded bootstrap aggregating (bagging, from bootstrap aggregating)
        XiDBnzM = mean(XiDBnz,2); % mean of non-zero values in ensemble
        XiDBnzM(XiDBnzM<ensembleT) = 0; % threshold: set all parameters that have an inclusion probability below threshold to zero

        XiDBmean = mean(XiDBe,2);
        XiDBmedian = median(XiDBe,2);

        XiDBmean(XiDBnzM==0)=0; 
        XiDBmedian(XiDBnzM==0)=0; 

        XiDB(libEntry,iii) = XiDBmean;
        XiDBmed(libEntry,iii) = XiDBmedian;
        
        XiDBstd = std(XiDBe')';
        XiDBstd(XiDBnzM==0)=0; 
        XiDBs(libEntry,iii) = XiDBstd;
                        
    end
    Xi = XiDB;
%     Xi = XiDBmed;
end

%    bnds = norm(dXdt_reg).^2./abs(dXdt_reg'*Theta_reg)'.*M; 


%             while size(find(smallinds)) == size(Xi(:))
%                 lambda = lambda/2;
%                 smallinds = (abs(Xi)<lambda);
%             end
