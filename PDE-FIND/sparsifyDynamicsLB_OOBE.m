function [Xi,its,thrs_EL,XiDBs] = sparsifyDynamicsLB_OOBE(Theta_0,dXdt,lambda,gamma,M_0,LBp)
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
% and out of bag error method

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
    if ~isempty(M_0)
        M = M_0(libOutBS(:,iii));
    else
        M = [];
    end
    
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
    if ~isempty(M_0)
        [XiBias,its,thrs_EL] = sparsifyDynamics(Theta_0(:,libEntry),dXdt(:,iii),lambda,gamma,M_0(libEntry));
    else
        [XiBias,its,thrs_EL] = sparsifyDynamics(Theta_0(:,libEntry),dXdt(:,iii),lambda,gamma,[]);
    end
    Xi(libEntry,iii) = XiBias;
end

if ~LBp.DB
    XiDBs = zeros(size(Theta_0,2),n);
end

if LBp.DB
    XiDB = zeros(size(Theta_0,2),n);
    XiDBmed = zeros(size(Theta_0,2),n);
    XiDBs = zeros(size(Theta_0,2),n);
    for iii = 1:n
        libEntry = inclProbBS(:,iii)>ensT;
                
        rng(iii,'twister')
        
        if ~isempty(M_0)
            MIn = M_0(libEntry);
        else
            MIn =[];
        end
%         MIn = [];
%         bootstatDD = bootstrp(nEnsemblesDD,@(ThetaIn,dxIn)sparsifyDynamics(ThetaIn,dxIn,lambda,gamma,MIn),Theta_0(:,libEntry),dXdt(:,iii)); 
        [bootstatDD,bootstatDDn] = bootstrp(nEnsemblesDD,@(ThetaIn,dxIn)sparsifyDynamics(ThetaIn,dxIn,lambda,gamma,MIn),Theta_0(:,libEntry),dXdt(:,iii)); 
%         [bootstatDD,bootstatDDn] = bootstrp(nEnsemblesDD,@(ThetaIn,dxIn)sparsifyDynamics(ThetaIn,dxIn,lambda,gamma,M_0),Theta_0,dXdt(:,iii)); 


%         nEnsemble3 = round(0.6*size(Theta_0,1));
%         bootstatDD = zeros(nEnsemblesDD,sum(libEntry));
%         bootstatDDn = zeros(nEnsemble3,nEnsemblesDD);
% %         MIn = [];
%         for jjj = 1:nEnsemblesDD
%             rng(jjj,'twister')
%             rs = RandStream('mlfg6331_64','Seed',jjj); 
%             bootstatDDn(:,jjj) = datasample(rs,1:size(Theta_0,1),nEnsemble3,'Replace',false)';
%             bootstatDD(jjj,:) = sparsifyDynamics(Theta_0(bootstatDDn(:,jjj),libEntry),dXdt(bootstatDDn(:,jjj),iii),lambda,gamma,MIn);
%         end

%         true_nz_weightsOut ID 5 and 10
        
        % Find the unique values
        for jj = 1:nEnsemblesDD
            XX = [1:size(Theta_0,1), bootstatDDn(:,jj)'];
            nUnique = histc(XX, unique(XX));
            uniqueVals = find(nUnique == 1);
        
%             OOSeps(jj) = abs(1 - mean(abs((Theta_0(uniqueVals,libEntry)*bootstatDD(jj,:)')./dXdt(uniqueVals,iii))));
%             OOSeps(jj) = mean(((Theta_0(uniqueVals,libEntry)*bootstatDD(jj,:)') - dXdt(uniqueVals,iii))./dXdt(uniqueVals,iii));
%             OOSeps2(jj) = mean(((Theta_0(:,libEntry)*bootstatDD(jj,:)') - dXdt(:,iii))./dXdt(:,iii));
%             OOSeps(jj) = abs(1-mean(abs((Theta_0(uniqueVals,libEntry)*bootstatDD(jj,:)')./dXdt(uniqueVals,iii))));
            
            OOSeps(jj) = norm(Theta_0(uniqueVals,libEntry)*(bootstatDD(jj,:))' - dXdt(uniqueVals,iii))/norm(dXdt(uniqueVals,iii));
            OOSeps2(jj) = norm(Theta_0(bootstatDDn(:,jj),libEntry)*(bootstatDD(jj,:))' - dXdt(bootstatDDn(:,jj),iii))/norm(dXdt(bootstatDDn(:,jj),iii));
%             OOSeps3(jj) = norm(Theta_0(bootstatDDn(:,jj),libEntry)*(bootstatDD(jj,:)./MIn)' - dXdt(bootstatDDn(:,jj),iii))/norm(dXdt(bootstatDDn(:,jj),iii));
            
%             OOSeps2(jj) = mean((Theta_0(:,libEntry)*bootstatDD(jj,:)')./dXdt(:,iii));
%             OOSeps(jj) = abs(1 - mean(abs((Theta_0(uniqueVals,:)*bootstatDD(jj,:)')./dXdt(uniqueVals,iii))));
        end
        
%         warning('its a mess')
        
%         Theta_0(5,libEntry)*bootstatDD(jj,:)'
%         dXdt(5,iii)
        
%         meanOOSeps = mean(OOSeps)
%         meanOOSeps2 = mean(OOSeps2)

%         test = dXdt(uniqueVals,iii);
%         test = mean(((Theta_0(uniqueVals,libEntry)*bootstatDD(jj,:)') - dXdt(uniqueVals,iii))./dXdt(uniqueVals,iii))
%         test = mean(((Theta_0(:,libEntry)*bootstatDD(jj,:)') - dXdt(:,iii))./dXdt(:,iii))
        
%         figure
%         plot(test,'.')
%         
%         epsMax = 0.25;
%         epsMaxDelta = 0.01;
%         OOSsmall = [];
%         while isempty(OOSsmall)
%             OOSsmall = find(OOSeps < epsMax);
%             epsMax = epsMax + epsMaxDelta;
%         end

%         figure
%         plot(OOSeps,'.'); hold on
% %         plot([0 nEnsemblesDD], [epsMax epsMax])
%         
        [~,OOSsmall] = mink(OOSeps,round(0.1*nEnsemblesDD)); % choose 30% best models
%         [~,OOSlarge] = maxk(OOSeps,round(0.3*nEnsemblesDD)); % choose 30% best models
%         OOSsmall = OOSlarge;
        
%         figure
%         plot(OOSeps(OOSsmall),'b.'); hold on
%         plot(OOSeps(OOSlarge),'r.'); hold on
        
                 
        XiDBe = [];
        XiDBnz = [];
        for iE = 1:nEnsemblesDD
            XiDBe(:,iE) = reshape(bootstatDD(iE,:),size(Theta_0(:,libEntry),2),1);
%             XiDBe(:,iE) = reshape(bootstatDD(iE,:),size(Theta_0,2),1);
            XiDBnz(:,iE) = XiDBe(:,iE)~=0;
        end

        % Thresholded bootstrap aggregating (bagging, from bootstrap aggregating)
        XiDBnzM = mean(XiDBnz,2); % mean of non-zero values in ensemble
        XiDBnzM(XiDBnzM<ensembleT) = 0; % threshold: set all parameters that have an inclusion probability below threshold to zero

        XiDBmean = mean(XiDBe,2);
        XiDBmedian = median(XiDBe,2);
        XiDBmean2 = mean(XiDBe(:,OOSsmall),2);

        XiDBmean(XiDBnzM==0)=0; 
        XiDBmedian(XiDBnzM==0)=0; 
        XiDBmean2(XiDBnzM==0)=0; 
        
%         XiDB(libEntry,iii) = XiDBmean;
        XiDB(libEntry,iii) = XiDBmean2;
        XiDBmed(libEntry,iii) = XiDBmedian;
%         XiDB(:,iii) = XiDBmean2;
%         XiDBmed(:,iii) = XiDBmedian;
        
        XiDBstd = std(XiDBe')';
        XiDBstd(XiDBnzM==0)=0; 
        XiDBs(libEntry,iii) = XiDBstd;
%         XiDBs(:,iii) = XiDBstd;
                        
    end
    Xi = XiDB;
%     Xi = XiDBmed;
end

%    bnds = norm(dXdt_reg).^2./abs(dXdt_reg'*Theta_reg)'.*M; 


%             while size(find(smallinds)) == size(Xi(:))
%                 lambda = lambda/2;
%                 smallinds = (abs(Xi)<lambda);
%             end
