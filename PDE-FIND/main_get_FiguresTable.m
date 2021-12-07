
%% plot figures for table in paper

%% Load data

clc
clear all
close all

for pde_num = 1:6

    pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','rxn_diff','rxn_diff'};
    pde_Plot = {'burgers','KdV','KS','NLS','rxn_diff_u','rxn_diff_v'};


    if pde_num < 5
        load(['datasets/',pde_names{pde_num}])
    else
%         load('datasets/rxn_diff.mat')
        load('C:\Users\ufase\OneDrive - UW\Documents\GitHub\LargeData\rxn_diff.mat')
        tPlot = 100;
        U_exact = [];
        if pde_num == 5
            U_exact{1} = u(:,:,tPlot);
        else
            U_exact{1} = v(:,:,tPlot);
        end
        xs{1} = x';
        xs{2} = y';
    end
    
    for nP = [1 0]
        coarsen_data = [[1 1];[1 1];[1 1]];
        if nP
            if pde_num == 4
                sigma_NR = 0.5; 
            elseif pde_num == 5 || pde_num == 6
                sigma_NR = 0.2;
            else
                sigma_NR = 1.0;
            end
        else
            sigma_NR = 0;
        end

        noise_dist = 0; 
        noise_alg = 0;

        nNoise = 1;

        iN = 1;
        rng(iN,'twister')
        rng_seed = rng().Seed; 

        nstates = length(U_exact);
        dim = length(size(U_exact{1}));
        inds = cell(1,dim);
        for j=1:dim
            N = length(xs{j});
            inds{j} = 1:coarsen_data(j,1):floor(N/coarsen_data(j,2));
            xs{j} = xs{j}(inds{j});
        end
        for j=1:nstates
            U_exact{j} = U_exact{j}(inds{:});
        end

        rng(rng_seed);
        [U_obs,noise,snr,sigma] = gen_noise(U_exact,sigma_NR,noise_dist,noise_alg,rng_seed,0);
        dims = size(U_obs{1});
        dim = length(dims);


        C1 = [0 119 187]/255;
        C2 = [51 187 238]/255;
        C3 = [0 153 136]/255;
        C4 = [238 119 51]/255;
        C5 = [204 51 17]/255;
        C6 = [238 51 119]/255;
        C7 = [187 187 187]/255;
        C8 = [80 80 80]/255;
        C9 = [140 140 140]/255;
        C10 = [0 128 255]/255;

        yellow = [255 255 153]/255;

        mymap = [C2
            C10
            0 0 0
            C5
            yellow];

        interpMap = interp1(1:5,mymap,linspace(1,5,200));

        labelS = 24;
        sizeX = 800;
        if pde_num == 5 || pde_num == 6
            sizeY = 300;
        else
            sizeY = 600;
        end
        f = figure('Position', [10 10 sizeX sizeY]);
        set(gcf,'units','normalized')
        set(gca,'LooseInset',get(gca,'TightInset'));
%         surf(xs{1},xs{2},U_obs{1}', 'EdgeColor','none')
        surf(xs{1},xs{2},U_obs{1}', 'EdgeColor','none','FaceColor','interp')
        xlim([min(xs{1}) max(xs{1})])
        ylim([min(xs{2}) max(xs{2})])
        if nP == 1
            zlimI = [min(min(U_obs{1}))*1.1 max(max(U_obs{1}))*1.1];
        end
        zlim(zlimI)
        colormap(autumn)
        colormap(interpMap)
        view([15 55])
        % view([0 90])
        xlabel('$x$','interpreter','latex','fontsize',labelS)
        ylabel('$t$','interpreter','latex','fontsize',labelS)
        % set(gca, 'TickLabelInterpreter','latex','fontsize',labelS)
        xticks([])
        yticks([])
        zticks([])
        axis off
%         if pde_num == 5 || pde_num == 6
%             view([15 15])
%             axis equal
%         end
        saveas(gcf,['PlotsPaper/',pde_Plot{pde_num}, num2str(nP),'3D','.png'])
    end
end
