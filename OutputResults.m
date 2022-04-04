clear;
clc;


R0Delta=5.08;

pA=0.308; % Proportion of asymptomatic


Risk=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No test (NatComm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_NoTest_NatComm.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RNoTest=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);

fprintf('=================================\n');
fprintf('No test \n');
fprintf('=================================\n');
fprintf('Expected transmission No Test: %3.2f \n',RNoTest(1));
fprintf('Probability of transmission No Test: %3.0f \n',100.*Probability_Onward(RNoTest(1),Risk));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-departure RT-PCR testing. (NatComm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RTPCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_RTPCR_NatComm.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RRTPCR=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BinaxNOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_BinaxNow__NatComm.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RBinax=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDVeritor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_BD Veritor__NatComm.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RBDVeritor=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sofia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_Sofia (FDA)__NatComm.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RSofia=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lumarax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_LumiraDX (Anterior Nasal Swab)__NatComm.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RLumarax=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CareStart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Pre_Testing_CareStart (Anterior Nasal Swab - FDA)__NatComm.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RCareStart=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);



fprintf('=================================\n');
fprintf('RT-PCR test \n');
fprintf('=================================\n');
fprintf('Expected transmission RT-PCR 24 hr: %3.2f \n',RRTPCR(w==1));
fprintf('Expected transmission RT-PCR 12 hr: %3.2f \n',RRTPCR(w==0.5));
fprintf('Expected transmission RT-PCR 72 hr: %3.2f \n',RRTPCR(w==3));
fprintf('Expected transmission RT-PCR 1 week: %3.2f \n',RRTPCR(w==7));


fprintf('Probability transmission RT-PCR 24 hr: %3.1f \n',100.*Probability_Onward(RRTPCR(w==1),Risk));
fprintf('Probability transmission RT-PCR 1 week: %3.1f \n',100.*Probability_Onward(RRTPCR(w==7),Risk));


fprintf('=================================\n');
fprintf('Rapid Antigen test \n');
fprintf('=================================\n');

load('Table_NatComm.mat','EPT','ProbTrans');

fprintf('Range Expected transmission 5 popular RA test: %3.2f  to %3.2f \n',[min([RLumarax(w==0) RSofia(w==0) RBinax(w==0) RBDVeritor(w==0) RCareStart(w==0)]) max([RLumarax(w==0) RSofia(w==0) RBinax(w==0) RBDVeritor(w==0) RCareStart(w==0)])] );

fprintf('Range Probability transmission 5 popular RA test: %3.0f  to %3.0f \n',100.*Probability_Onward([min([RLumarax(w==0) RSofia(w==0) RBinax(w==0) RBDVeritor(w==0) RCareStart(w==0)]) max([RLumarax(w==0) RSofia(w==0) RBinax(w==0) RBDVeritor(w==0) RCareStart(w==0)])],Risk) );

fprintf('Median (Range) Expected transmission ALL RA test: %3.2f (%3.2f  to %3.2f) \n',[ median(EPT(2:end,end)) min(EPT(2:end,end)) max(EPT(2:end,end))]);
fprintf('Median probability transmission ALL RA test: %3.0f (%3.0f  to %3.0f) \n',100.*[ median(ProbTrans(2:end,end)) min(ProbTrans(2:end,end)) max(ProbTrans(2:end,end))]);


fRAp=find(EPT(2:end,end)>EPT(1,end-1));

fprintf('Number of RA tests that do worse than 24 h RT-PCR: %2.0f \n',length(fRAp));

fprintf('Time RT-PCR conducted to do better than 50%% of RA: %3.1f \n',24*w(max(find(RRTPCR<median(EPT(2:end,end))))))


RedC=1-RRTPCR(w==0)./RRTPCR(w==7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RTPCR (Hellewell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf('=================================\n');
fprintf('Hellewell \n');
fprintf('=================================\n');

load('Pre_Testing_RTPCR_Hellewell.mat','w','ts','td','PosttestIS','PosttestIA','PretestIS','PretestIA','R0A','R0S');
RS=(w./ts).*PosttestIS+((ts-w)./ts).*PretestIS;
RA=(w./td).*PosttestIA+((td-w)./td).*PretestIA;

RS(w==0)=PretestIS(w==0);
RA(w==0)=PretestIA(w==0);

RRTPCRH=pA.*RA.*(R0Delta./R0A)+(1-pA).*RS.*(R0Delta./R0S);


RedH=1-RRTPCRH(w==0)./RRTPCRH(w==7);

fprintf('Percent reduction Hellewell vs Wells: %3.0f vs %3.0f  \n',100.*[RedH RedC]);


load('Table_Hellewell.mat','EPT','ProbTrans');


fRAp=find(EPT(2:end,end)>EPT(1,end-1));

fprintf('Number of RA tests that do worse than 24 h RT-PCR for Hellewell: %2.0f \n',length(fRAp));


fprintf('Median probability transmission ALL RA test Hellewell: %3.0f (%3.0f  to %3.0f) \n',100.*[ median(ProbTrans(2:end,end)) min(ProbTrans(2:end,end)) max(ProbTrans(2:end,end))]);

fprintf('Time RT-PCR conducted to do better than 50%% of RA for Hellewell: %3.1f \n',24*w(max(find(RRTPCRH<median(EPT(2:end,end))))))