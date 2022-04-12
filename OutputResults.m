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

MLE_Red=zeros(NumTests,1);
All_Red=zeros(NumTests,1000);

MLE_PT=zeros(NumTests,1);
All_PT=zeros(NumTests,1000);

vv=1;


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
    RNoTest=RNoTest(w==0);
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
    
    RA=repmat(RA(w==0),1000,1);
    RS=repmat(RS(w==0),1000,1);

    RNoTestv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);
    
    
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


        Risk=1;
        MLE_Red(jj,:)=1-R_RATest(w==0)./RNoTest;   
        MLE_PT(jj,:)=[100.*Probability_Onward(R_RATest(w==0),Risk);];
        
        load(['Pre_Testing_' testName{jj} '_' VOCv{vv} '_Uncertainty.mat'],'PosttestISv','PosttestIAv','PretestISv','PretestIAv');
        w=repmat([0 0.5 1 2 3],1000,1);
        RS=(w./ts).*PosttestISv+((ts-w)./ts).*PretestISv;
        RA=(w./td).*PosttestIAv+((td-w)./td).*PretestIAv;

        RS(w==0)=PretestISv(w==0);
        RA(w==0)=PretestIAv(w==0);


        RS(w==ts)=PosttestISv(w==ts);
        RA(w==td)=PosttestIAv(w==td);

        R_RATestv=pAv.*RA.*(R0VOC./R0A)+(1-pAv).*RS.*(R0VOC./R0S);

        R_temp=1-R_RATestv./RNoTestv; 


        All_Red(jj,:)=R_temp(w==0);

        All_PT(jj,:)=100.*Probability_Onward(R_RATestv(w==0),Risk);
    end