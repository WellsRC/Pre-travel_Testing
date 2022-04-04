function [betaAg]=ParameterCOVIDTest(testName,MLEv)
betaAg=[];
if(~isempty(testName))
    if(MLEv~=0)
        load([testName '_LR_Parameters.mat'],'beta');
        betaAg=beta;
    else
        load([testName '_LR_Uncertainty.mat'],'L','betaS','beta');
        w=cumsum(exp(L)./sum(exp(L)));
        r=rand(1);
        findx=find(w>=r,1);
        betaAg=beta.*(1+betaS(findx,:));
    end
end

end