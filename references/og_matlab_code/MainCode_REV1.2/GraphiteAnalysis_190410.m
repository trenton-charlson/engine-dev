T_wg_min = 290; %Min wall temp considered
T_wg_max = 1500; %Max wall temp considered

[r, ~] = size(engineProps);

nStat = 125;

xVec = (nStat:1:r)';

qconv_minTVec = zeros(r,1);
qconv_maxTVec = zeros(r,1);

for i = 1:r
[~,qconv_minT,~] = bartzCorrelation(engineProps,T_wg_min,R_t,P_c,C_star,R_tCurve,i);
qconv_minTVec(i) = qconv_minT;
[~,qconv_maxT,~] = bartzCorrelation(engineProps,T_wg_max,R_t,P_c,C_star,R_tCurve,i);
qconv_maxTVec(i) = qconv_maxT;
end

qconv_minTVec = qconv_minTVec(nStat:end);
qconv_maxTVec = qconv_maxTVec(nStat:end);
qconv_avg = mean([qconv_minTVec';qconv_maxTVec'])';

plot(xVec,qconv_minTVec,'r')
hold on
plot(xVec,qconv_maxTVec,'b')
hold on
plot(xVec,qconv_avg,'k')

ALocal = [];
for i = nStat:r
    R = engineProps(i,1);
    dZ = engineProps(i,2) - engineProps(i-1);
    A = pi()*(R^2)*dZ;
    ALocal(end+1,1) = A;
end
E = sum(ALocal.*qconv_avg)
cpGraph = 0.71*1000; %J/kgK
mGraphite = 1.1; %kg
fireTime = 2; %Seconds - worst possible case

dT = (E*fireTime)/(cpGraph*mGraphite)

L_graph = 0.11688; %m
CTE_graph = 8*10^-6; %worst case - m/mK
R_graph = 0.09871; %m
dR_graph = R_graph*CTE_graph*dT

dL_graph = L_graph*CTE_graph*dT





