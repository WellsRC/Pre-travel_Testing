function [f,F] = DistIncubation(Inc)
% https://www.eurosurveillance.org/docserver/fulltext/eurosurveillance/26/50/eurosurv-26-50-2.pdf?expires=1647872499&id=id&accname=guest&checksum=19A9B093BCF9B4114BB7835B8A4FACA9
% m=3.1;
% v=2.6^2;
% p=[log(m^2/sqrt(v+m^2)) sqrt(log(v/m^2+1))];
% [pp,fval]=fmincon(@(x)[(logninv(0.5,x(1),x(2))-3).^2+ (logninv(0.75,x(1),x(2))-4).^2+ (logninv(0.975,x(1),x(2))-8).^2],p)

p=[1.06453632885258,0.516594915300738];
 
 f=lognpdf(Inc,p(1),p(2));
 F=logncdf(Inc,p(1),p(2));

end

