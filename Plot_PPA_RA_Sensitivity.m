clear;
close all;
clc;

addpath('PPA_RA_Tests');
load('RAgTest_Name.mat','testName');
NumTests=length(testName);

load('PCR_mapping.mat','VmI','VmS','t1','t2');
PCR_Map.VmI=VmI;
PCR_Map.VmS=VmS;
PCR_Map.t1=t1;
PCR_Map.t2=t2;

load('Peak_Infection.mat','mmv','tsv');
ts=[3.1 4.4 5.723];

for ii=1:NumTests
    figure('units','normalized','outerposition',[0 0 1 1]);
    for jj=1:3
        
        [MLE_RTPCR,MLE_Ag,U_RTPCR,U_Ag,MLE_PPA,U_PPA,Dt,totalpos,truepos,w,t] = Sensitivity_for_Plotting(testName{ii},ts(jj),pchip(tsv,mmv,ts(jj)),PCR_Map);
        switch jj
            case 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                % PPA
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                subplot('Position',[0.061491596638655,0.562256332320174,0.424600840336134,0.390000000000001]);
                
                patch([t flip(t)],[prctile(U_PPA,2.5) flip(prctile(U_PPA,97.5))],hex2rgb('#004445'),'LineStyle','none','Facealpha',0.2); hold on;
                p1=plot(t,MLE_PPA,'-','color',hex2rgb('#004445'),'LineWidth',2); hold on;
                p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,'o','filled','MarkerEdgeColor',hex2rgb('#004445'),'MarkerFaceColor',hex2rgb('#004445'));
                scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,'o','LineWidth',2,'MarkerEdgeColor',hex2rgb('#004445'));

                set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
                ytickformat('percentage')
                xlabel('Days since symptom onset','Fontsize',18);
                ylabel({'Percent','positive agreement'},'Fontsize',18,'Position',[-2.94470530625081,50,-0.999999999999986]);
                
                text(-5.697,98.8,'A','Fontsize',30,'FontWeight','bold');
                text(42.33621755253399,108.051948051948,AdjustedNames_Plotting(testName{ii}),'Fontsize',28,'FontWeight','bold','HorizontalAlignment','center');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                % Omicron
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                subplot('Position',[0.561491596638655,0.562256332320174,0.424600840336134,0.390000000000001]);
            case 2
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                % Delta
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                subplot('Position',[0.061491596638655,0.084,0.424600840336134,0.390000000000001]);
            case 3
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                % Original
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                subplot('Position',[0.561491596638655,0.084,0.424600840336134,0.390000000000001]);
        end
        
        patch([t flip(t)],[prctile(U_RTPCR,2.5) flip(prctile(U_RTPCR,97.5))],'k','LineStyle','none','Facealpha',0.2); hold on;
        patch([t flip(t)],[prctile(U_Ag,2.5) flip(prctile(U_Ag,97.5))],hex2rgb('#004445'),'LineStyle','none','Facealpha',0.2); 
        p1=plot(t,MLE_RTPCR,'-.','color','k','LineWidth',2); 
        p2=plot(t,MLE_Ag,'-','color',hex2rgb('#004445'),'LineWidth',2);
        plot(ts(jj).*ones(101,1),linspace(0,1,101),'-.','LineWidth',1.5,'color',[0.75 0.75 0.75]);
        set(gca,'LineWidth',2,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:0.1:1],'Ylim',[0 1]);
        ylabel({'Diagnostic sensitivity'},'Fontsize',18);
        xlabel('Days post-infection','Fontsize',18);
        text(-6.54739985398804,0.988,char(65+jj),'Fontsize',30,'FontWeight','bold');
        legend([p1 p2],{'RT-PCR','Rapid antigen test'},'Fontsize',18);
        legend boxoff;
    end
print(gcf,['PPA_Sensitivity_' testName{ii} '.png'],'-dpng','-r300');
end

rmpath('PPA_RA_Tests');