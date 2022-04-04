function [betaRTPCR]=ParameterRTPCRTest(MLEv)
if(MLEv~=0)
            load(['MLE-Estimate-RTPCR.mat'],'beta')
    betaRTPCR=beta;
else
    load(['Uncertainty-Beta-Estimate-RTPCR.mat'],'L','betaRTU');
    w=cumsum(exp(L)./sum(exp(L)));
    r=rand(1);
    findx=find(w>=r,1);
    betaRTPCR=betaRTU(findx,:);
end
end