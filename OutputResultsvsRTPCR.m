clear;
clc;
R0v=[5.08*(1+0.365) 5.08 2.79];
VOCv={'Omicron','Delta','Original'};
rng default;
pA=0.351; % Proportion of asymptomatic
par=fminbnd(@(x)((betainv(0.025,x,x*(1-pA)./pA)-0.307).^2+(betainv(0.975,x,x*(1-pA)./pA)-0.399).^2),[0],[400]);
pAv=betarnd(par,par*(1-pA)./pA,1000,1);

load('RAgTest_Name_Final.mat','testName');
NumTests=length(testName);

MLE_T=zeros(NumTests,1);
All_T=zeros(NumTests,1000);

vv=1;
w_test=0.5;

    R0VOC=R0v(vv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % No test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(['Pre_Testing_RTPCR_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
    RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
    RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

    RS(w==0)=PretestIS(w==0);
    RA(w==0)=PretestIA(w==0);
    
    RS(w==ts)=PosttestIS(w==ts);
    RA(w==td)=PosttestIA(w==td);

    RRTPCR=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);
    RRTPCR=RRTPCR(w==w_test);
   
    
    load(['Pre_Testing_RTPCR_' VOCv{vv} '_Uncertainty.mat'],'PosttestISv','PosttestIAv','PretestISv','PretestIAv');
    w=repmat([0 0.5 1 2 3 7],1000,1);
    RS=(w./ts).*PosttestISv+((ts-w)./ts).*PretestISv;
    RA=(w./td).*PosttestIAv+((td-w)./td).*PretestIAv;

    RS(w==0)=PretestISv(w==0);
    RA(w==0)=PretestIAv(w==0);
    
    
    RS(w==ts)=PosttestISv(w==ts);
    RA(w==td)=PosttestIAv(w==td);

    RRTPCRv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);
    
    RRTPCRv=RRTPCRv(w==w_test);
    
    for jj=1:NumTests
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RA TEST
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(['Pre_Testing_' testName{jj} '_' VOCv{vv} '.mat'],'w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
        RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
        RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

        RS(w==0)=PretestIS(w==0);
        RA(w==0)=PretestIA(w==0);


        RS(w==ts)=PosttestIS(w==ts);
        RA(w==td)=PosttestIA(w==td);

        R_RATest=pA.*RA.*(R0VOC./R0A)+(1-pA).*RS.*(R0VOC./R0S);
        if(R_RATest(w==0)>RRTPCR)
            MLE_T(jj)=1; 
        end
        
        load(['Pre_Testing_' testName{jj} '_' VOCv{vv} '_Uncertainty.mat'],'PosttestISv','PosttestIAv','PretestISv','PretestIAv');
        w=repmat([0 0.5 1 2 3 7],1000,1);
        RS=(w./ts).*PosttestISv+((ts-w)./ts).*PretestISv;
        RA=(w./td).*PosttestIAv+((td-w)./td).*PretestIAv;

        RS(w==0)=PretestISv(w==0);
        RA(w==0)=PretestIAv(w==0);


        RS(w==ts)=PosttestISv(w==ts);
        RA(w==td)=PosttestIAv(w==td);

        R_RATestv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);

        R_temp=R_RATestv-RRTPCRv; 


        All_T(jj,:)=R_temp(w==0);
    end
    
    All_T(All_T>0)=1;
    All_T(All_T<=0)=0;