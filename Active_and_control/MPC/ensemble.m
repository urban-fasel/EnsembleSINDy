function [Xi, Xis, XiE] = ensemble(Theta_0,dxobs_0,lambda)

% ensemble hyperparameters
% data ensembling
nEnsembles = 100; % number of bootstraps (SINDy models using sampled data) in ensemble
ensembleT = 0.6; % Threshold model coefficient inclusion probability: set ensemble SINDy model coefficient to zero if inclusion probability is below ensembleT

% library
nEnsemble1P = 0.9; % percentage of full library that is sampled without replacement for library bagging
nEnsemble2 = 100; % number of bootstraps (SINDy models using sampled library terms) in ensemble
ensT = 0.4; % Threshold library term inclusion probabilty: cut library entries that occur less than ensT

% double bagging
nEnsemblesDD = 100; % number of models in ensemble for data bagging after library bagging

n = size(dxobs_0,2);

%% double bagging SINDy
            
%% Bagging SINDy library
% randomly sample library terms without replacement and throw away terms
% with low inclusion probability
nEnsemble1 = round(nEnsemble1P*size(Theta_0,2));
mOutBS = zeros(nEnsemble1,n,nEnsemble2);
libOutBS = zeros(nEnsemble1,nEnsemble2);
for iii = 1:nEnsemble2
    rs = RandStream('mlfg6331_64','Seed',iii); 
    libOutBS(:,iii) = datasample(rs,1:size(Theta_0,2),nEnsemble1,'Replace',false)';
%     mOutBS(:,:,iii) = sparsifyDynamics(Theta_0(:,libOutBS(:,iii)),dxobs_0,lambda,n,gamma);
    mOutBS(:,:,iii) = sparsifyDynamics(Theta_0(:,libOutBS(:,iii)),dxobs_0,lambda,n);
end

inclProbBS = zeros(size(Theta_0,2),n);
for iii = 1:nEnsemble2
    for jjj = 1:n
        for kkk = 1:nEnsemble1
            if mOutBS(kkk,jjj,iii) ~= 0
                inclProbBS(libOutBS(kkk,iii),jjj) = inclProbBS(libOutBS(kkk,iii),jjj) + 1;
            end
        end
    end
end
inclProbBS = inclProbBS/nEnsemble2*size(Theta_0,2)/nEnsemble1;

XiD = zeros(size(Theta_0,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;
%     XiBias = sparsifyDynamics(Theta_0(:,libEntry),dxobs_0(:,iii),lambda,1,gamma);
    XiBias = sparsifyDynamics(Theta_0(:,libEntry),dxobs_0(:,iii),lambda,1);
    XiD(libEntry,iii) = XiBias;
end

                
%% Double bagging SINDy 
% randomly sample library terms without replacement and throw away terms
% with low inclusion probability
% then on smaller library, do bagging

XiDB = zeros(size(Theta_0,2),n);
XiDBmed = zeros(size(Theta_0,2),n);
XiDBs = zeros(size(Theta_0,2),n);
XiDBeOut = zeros(size(Theta_0,2),n,nEnsemblesDD);
inclProbDB = zeros(size(Theta_0,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;

%     bootstatDD = bootstrp(nEnsemblesDD,@(Theta,dx)sparsifyDynamics(Theta,dx,lambda,1,gamma),Theta_0(:,libEntry),dxobs_0(:,iii)); 
    bootstatDD = bootstrp(nEnsemblesDD,@(Theta,dx)sparsifyDynamics(Theta,dx,lambda,1),Theta_0(:,libEntry),dxobs_0(:,iii)); 
    
    XiDBe = [];
    XiDBnz = [];
    for iE = 1:nEnsemblesDD
        XiDBe(:,iE) = reshape(bootstatDD(iE,:),size(Theta_0(:,libEntry),2),1);
        XiDBnz(:,iE) = XiDBe(:,iE)~=0;
        
        XiDBeOut(libEntry,iii,iE) = XiDBe(:,iE);
    end

    % Thresholded bootstrap aggregating (bagging, from bootstrap aggregating)
    XiDBnzM = mean(XiDBnz,2); % mean of non-zero values in ensemble
    inclProbDB(libEntry,iii) = XiDBnzM;
    XiDBnzM(XiDBnzM<ensembleT) = 0; % threshold: set all parameters that have an inclusion probability below threshold to zero

    XiDBmean = mean(XiDBe,2);
    XiDBmedian = median(XiDBe,2);
    XiDBstd = std(XiDBe')';

    XiDBmean(XiDBnzM==0)=0; 
    XiDBmedian(XiDBnzM==0)=0; 
    XiDBstd(XiDBnzM==0)=0; 
    
    XiDB(libEntry,iii) = XiDBmean;
    XiDBmed(libEntry,iii) = XiDBmedian;
    XiDBs(libEntry,iii) = XiDBstd;
    
end

Xi = XiDBmed;
Xis = XiDBs;
XiE = XiDBeOut;

