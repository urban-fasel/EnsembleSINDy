function plotUQ_Lorenz(XiE,XiTrue,Xi,lib)

C1 = [0 119 187]/255;
C2 = [51 187 238]/255;
C3 = [0 153 136]/255; % green
C4 = [238 119 51]/255;
C5 = [204 51 17]/255;
C6 = [238 51 119]/255;
C7 = [187 187 187]/255;
C8 = [80 80 80]/255;
C9 = [140 140 140]/255;
C10 = [0 128 255]/255;

lw = 1.0;
lwlb = 1.0;
lwlbA = 1.0;

lss = '-'; % line style sindy
lslb = '-'; % line style library bagging
lslb1 = '-';

Ct = C8;%C7;%C9; % 'r'
Cs = C5;
Clb = C10;
ClbA = Clb;
Clb2 = C7;

colors = viridis(100);
% Cs = colors(98,:); 
% Clb = colors(70,:); 

fos = 10; % fontsize
fosT = 14; % fontsize

plotAllDist = 1; % plot all or only thresholded distributions 

sizeX = 600;
sizeY = 500;

figure('Position', [10 10 sizeX sizeY])

nbins = 10;
yBB = 400; % y blue background
nP1 = size(XiE,1);
nP2 = size(XiE,2);

% xmin = min(XiOut(:,:,iP))-1;%0.2*abs(min(XiOut(:,:,iP)));
% xmax = max(XiOut(:,:,iP))+1;%0.2*abs(max(XiOut(:,:,iP)));
xmin = [-12 -3 -5];%0.2*abs(min(XiOut(:,:,iP)));
xmax = [12 30 3];%0.2*abs(max(XiOut(:,:,iP)));
nn = 1;
for ip = 1:nP1
    for jp = 1:nP2
        xPtest = XiE(ip,jp,:);
        xPtest = xPtest(:);
        
        subplot(nP1,nP2,nn)
        if (XiTrue(ip,jp)~=0)
%             axis on
            patch([xmin(jp) xmin(jp) xmax(jp) xmax(jp)],[0 yBB yBB 0],Clb,'FaceAlpha',0.1); hold on
        end
        if plotAllDist
            if sum(xPtest~=0) > 0
                h = histfit(xPtest,nbins,'normal'); hold on
                h(1).FaceColor = Clb2;
                h(1).FaceAlpha = 0;%0.5;
                h(1).EdgeColor = Clb2;
                h(1).EdgeAlpha = 0;
                h(2).Color = Clb2;
                h(2).LineWidth = 1;
                xlim([xmin(jp) xmax(jp)])
                mY = max(h(2).YData);
                area(h(2).XData,h(2).YData,'FaceColor',Clb,'FaceAlpha',0.4); hold on
                ylim([0 mY])
                if Xi(ip,jp) ~= 0
                    plot(XiTrue(ip,jp),0,'x','Color',Ct,'Linewidth',lw,'MarkerSize',8); hold on
                    plot([Xi(ip,jp) Xi(ip,jp)],[0 mY],'-','Color',Clb,'Linewidth',lw); hold on
                else
                    plot([mean(xPtest) mean(xPtest)],[0 mY],'-','Color',Clb,'Linewidth',lw); hold on
                end
            else
                plot([0 0],[0 1],'-','Color',Clb,'Linewidth',lw)
                ylim([0 1])
            end
%             axis off
        else
             if Xi(ip,jp) ~= 0
    %             h = histfit(xPtest,nbins,'kernel');
                h = histfit(xPtest,nbins,'normal');
                h(1).FaceColor = Clb2;
                h(1).FaceAlpha = 0;%0.5;
                h(1).EdgeColor = Clb2;
                h(1).EdgeAlpha = 0;
                h(2).Color = Clb2;
                mY = max(h(2).YData);
                ylim([0 mY])
                if (XiTrue(ip,jp)~=0)
                    plot(XiTrue(ip,jp),0,'x','Color',Ct,'Linewidth',lw); hold on
                    plot([Xi(ip,jp) Xi(ip,jp)],[0 mY],'-','Color',Clb,'Linewidth',lw); hold on
                end
             else
                plot([0 0],[0 1],'-','Color',Clb,'Linewidth',lw)
                ylim([0 1])
             end
        end
        xlim([xmin(jp) xmax(jp)])
      
        yticks([])
        if ip~=nP1
            xticks([])
        end
        if ip==nP1
            if jp == 1
                xlabel('Coefficient value (xDot)','interpreter','latex','FontSize',fos)
            elseif jp == 2
                xlabel('Coefficient value (yDot)','interpreter','latex','FontSize',fos)
            else
                xlabel('Coefficient value (zDot)','interpreter','latex','FontSize',fos)
            end
            set(gca,'ticklabelinterpreter','latex','FontSize',fos)
        end
        
        if ip == 1
        	if jp == 1
                title('xDot','interpreter','latex','FontSize',fos)
            elseif jp == 2
                title('yDot','interpreter','latex','FontSize',fos)
            else 
                title('zDot','interpreter','latex','FontSize',fos)
            end
        end
        
        if jp == 1
            ylabel(lib{ip,1},'interpreter','latex','FontSize',fos)
        end
        
        nn = nn + 1;
    end
end

sgtitle('Uncertainty in model parameters','interpreter','latex','FontSize',fosT)

