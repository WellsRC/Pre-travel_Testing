close all;
clear;
clc;
% R0 for Delta (The reproductive number of the Delta variant of SARS-CoV-2
% is far higher compared to the ancestral SARS-CoV-2 virus)

R0v=[5.08*(1+0.365) 5.08 2.79];
VOCv={'Omicron','Delta','Original'};

pA=0.351; % Proportion of asymptomatic
for vv=1:3
    R0VOC=R0v(vv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % No test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_No_Testing_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RNoTest=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RTPCR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_RTPCR_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RRTPCR=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BinaxNOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_BinaxNow_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RBinax=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BDVeritor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_BD Veritor_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RBDVeritor=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sofia
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_Sofia (FDA)_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RSofia=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lumarax
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_LumiraDX (Anterior Nasal Swab)_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RLumarax=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CareStart
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_CareStart (Anterior Nasal Swab - FDA)_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RCareStart=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);

    Lcolour=[hex2rgb('#9A9EAB');hex2rgb('#000B29'); hex2rgb('#807dba'); hex2rgb('#CF3721'); hex2rgb('#4CB5F5'); hex2rgb('#F5BE41');  hex2rgb('#A2C523');];
    Lstyle={':','-','-.','-.','-.','-.','-.'};

    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);

    subplot('Position',[0.111310190369541,0.136779071547693,0.399999999999999,0.81]);

    pdt=plot(w,100.*(1-RRTPCR./RNoTest),w,100.*(1-RBDVeritor./RNoTest),w,100.*(1-RBinax./RNoTest),w,100.*(1-RCareStart./RNoTest),w,100.*(1-RLumarax./RNoTest),w,100.*(1-RSofia./RNoTest),'LineWidth',2);
    box off;
    for ii=1:6
       pdt(ii).Color=Lcolour(ii+1,:); 
       pdt(ii).LineStyle=Lstyle{ii+1}; 
    end
    set(gca,'LineWidth',2,'tickdir','out','XTick',[0:7],'XminorTick','on','YTick',[0:10:100],'Yminortick','on','Fontsize',22, 'xdir', 'reverse' );
    ytickformat('percentage');
    xlim([0 7]);
    ylim([0 50])
    text(8.83660789473684,83.07199999999997./80.*50,'A','Fontsize',34','FontWeight','bold');
    legend({'RT-PCR','BD Veritor','BinaxNOW','CareStart','LumiraDx','Sofia'},'Fontsize',22,'Location','NorthWest')
    legend boxoff;

    xlabel('Days pre-arrival test conducted','Fontsize',28);
    ylabel({'Reduction in post-arrival','transmission from no testing'},'Fontsize',28);
    % print(gcf,['PreDeparture_Expected_Transmission_NatComm.png'],'-dpng','-r600');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Probability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Risk=1;

    % figure('units','normalized','outerposition',[0.1 0.1 0.4 0.8]);

    subplot('Position',[0.592236842105263,0.136779071547693,0.4,0.81]);

    pdt=plot(w,100.*Probability_Onward(RNoTest,Risk),w,100.*Probability_Onward(RRTPCR,Risk),w,100.*Probability_Onward(RBDVeritor,Risk),w,100.*Probability_Onward(RBinax,Risk),w,100.*Probability_Onward(RCareStart,Risk),w,100.*Probability_Onward(RLumarax,Risk),w,100.*Probability_Onward(RSofia,Risk),'LineWidth',2);
    box off;
    for ii=1:7
       pdt(ii).Color=Lcolour(ii,:); 
       pdt(ii).LineStyle=Lstyle{ii}; 
    end
    set(gca,'LineWidth',2,'tickdir','out','XTick',[0:7],'XminorTick','on','YTick',[0:5:45],'Yminortick','on','Fontsize',22, 'xdir', 'reverse');
    ylim([0 42])
    ytickformat('percentage');
    xlim([0 7]);
    xlabel('Days pre-arrival test conducted','Fontsize',28);
    ylabel({'Probability of post-arrival transmission'},'Fontsize',28);
    legend({'No test','RT-PCR','BD Veritor','BinaxNOW','CareStart','LumiraDx','Sofia'},'Fontsize',22,'Location','SouthWest')
    legend boxoff;
    text(8.216697819314641,43.5376,'B','Fontsize',34','FontWeight','bold');
    print(gcf,['PreDeparture_' VOCv{vv} '.png'],'-dpng','-r600');
end