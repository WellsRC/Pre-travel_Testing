function p = PCRSens(t,beta,tL)
 V = ViralShedding_Symptomatic(t,tL,inf);
  p= V.^beta(1)./(V.^beta(1)+beta(2));
end

