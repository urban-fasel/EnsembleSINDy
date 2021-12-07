function XL = plotUQ_LorenzActive(XiE,XiTrue,Xi,lib,iLoop,XL)

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

Ct = C8;
Cs = C5;
Clb = C10;
ClbA = Clb;
Clb2 = C7;
C7b = 220*[1 1 1]/255;

%% colors figure1
% blue: x
blue = [44,127,184]./255;
% green: y
green = [49,163,84]./255;
% orange: z
orange = [240,59,32]./255;

colorsNew = [blue; green; orange];


fos = 14; 
fosT = 14;

plotAllDist = 1; % plot all or only thresholded distributions 

sizeX = 500;
sizeY = 600;

if iLoop == 1
    figure('Position', [10 10 sizeX sizeY])
end

nbins = 10;
yBB = 400; % y blue background
nP1 = size(XiE,1);
nP2 = size(XiE,2);

xmin = [-12 -3 -5];
xmax = [12 30 3];
nn = 1;
for ip = 1:nP1
    for jp = 1:nP2
        xPtest = XiE(ip,jp,:);
        xPtest = xPtest(:);
        
        Clb = colorsNew(jp,:);
        
        subplot(nP1,nP2,nn)
        if plotAllDist
            if sum(xPtest~=0) > 0
                h = histfit(xPtest,nbins,'kernel'); hold on
                h(1).FaceColor = Clb2;
                h(1).FaceAlpha = 0;
                h(1).EdgeColor = Clb2;
                h(1).EdgeAlpha = 0;
                if iLoop == 1
                    h(2).Color = C7b;
                    h(2).LineWidth = 0.5;
                else
                    h(2).Color = Clb2;
                    h(2).LineWidth = 1.5;
                end
                mY = 1;
                h(2).YData = h(2).YData/max(h(2).YData);
                if iLoop == 1
                    area(h(2).XData,h(2).YData,'FaceColor',Clb,'FaceAlpha',0.1); hold on
                else
                    area(h(2).XData,h(2).YData,'FaceColor',Clb,'FaceAlpha',0.4); hold on
                end
                ylim([0 mY])
                if iLoop == 1
                    XL(ip,jp,:) = [min(h(2).XData) max(h(2).XData)];
                else
                    xlim([XL(ip,jp,1) XL(ip,jp,2)]);
                end
                if Xi(ip,jp) ~= 0
                    plot(XiTrue(ip,jp),0,'x','Color',Ct,'Linewidth',lw,'MarkerSize',14); hold on
                    if iLoop == 1
                        plot([Xi(ip,jp) Xi(ip,jp)],[0 mY],':','Color',Clb,'Linewidth',lw); hold on
                    else
                        plot([Xi(ip,jp) Xi(ip,jp)],[0 mY],'-','Color',Clb,'Linewidth',lw); hold on
                    end
                else
                    plot([mean(xPtest) mean(xPtest)],[0 mY],'-','Color',Clb,'Linewidth',lw); hold on
                end
            else
                plot([0 0],[0 1],'-','Color',Clb,'Linewidth',lw)
                ylim([0 1])
            end
        else
             if Xi(ip,jp) ~= 0
                h = histfit(xPtest,nbins,'kernel');
                h(1).FaceColor = Clb2;
                h(1).FaceAlpha = 0;
                h(1).EdgeColor = Clb2;
                h(1).EdgeAlpha = 0;
                h(2).Color = Clb2;
                mY = 1;
                h(2).YData = h(2).YData/max(h(2).YData);
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
      
        yticks([])
        if ip~=nP1
            xticks([])
        else
            xticklabels({' ',' ',' '})
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',2)
        end
        if ip==nP1
            if jp == 1
            elseif jp == 2
                set(gca,'ticklabelinterpreter','latex','FontSize',11)
                xlabel('Identified model coefficient PDFs','interpreter','latex','FontSize',fos)
            end
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
            ylabel(lib{ip,1},'interpreter','latex','FontSize',fos,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
        end
        
        nn = nn + 1;
    end
end

