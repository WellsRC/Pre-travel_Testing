clear;
clc;

load('RAgTest_Name_Final.mat','testName');
NumTests=length(testName);

Name=cell(NumTests,1);
beta_0=cell(NumTests,1);
beta_1=cell(NumTests,1);
for jj=1:NumTests
    Name{jj}=AdjustedNames_Plotting(testName{jj});
    load(['Pre_Testing_' testName{jj} '_Omicron.mat'],'betaAg');
    betaAg_MLE=betaAg;
    load([testName{jj} '_Parameter_Sample_Uncertainty.mat'],'betaAgv');
    betaAg_LB_UB_0=prctile(betaAgv(:,1),[2.5 97.5]);
    betaAg_LB_UB_1=prctile(betaAgv(:,2),[2.5 97.5]);
    beta_0{jj}=[num2str(betaAg_MLE(1),'%3.2f') ' (' num2str(betaAg_LB_UB_0(1),'%3.2f') char(8211) num2str(betaAg_LB_UB_0(2),'%3.2f') ')'];               
    beta_1{jj}=[num2str(betaAg_MLE(2),'%3.2e') ' (' num2str(betaAg_LB_UB_1(1),'%3.2e') char(8211) num2str(betaAg_LB_UB_1(2),'%3.2e') ')']; 
end

writetable(table(Name,beta_0,beta_1),['beta_Estimate_PPA.csv']);