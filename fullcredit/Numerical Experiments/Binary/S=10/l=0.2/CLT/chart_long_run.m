hold on

xMin = 0;
xMax = 35350;
yMin = 0.0007;
yMax = 0.0014;
trueAns =  0.00108;
trueAnsCLT = 0.0010464342064605502112539792136658;
X = 1:xMax;

TwoLvlISX = [33222.1406,33222.1406,33222.1406];
        
TwoLvlISY = [0.0010854822534959573249180309062467, 0.0010829015699632953562830950389184, 0.0010791416233989011182192729876306];

TwoLvlISYMean = arrayfun(@(i) mean(TwoLvlISY((3*i-2):(3*i))), 1:(length(TwoLvlISY)/3));

OneLvlISCLTX = [3935.2656,3935.2656,3935.2656];
        
OneLvlISCLTY = [0.0010493684293107112884269049857267, 0.0010495767203081316984447646234457, 0.0010486001322307931648419865311439];

OneLvlISCLTYMean = arrayfun(@(i) mean(OneLvlISCLTY((3*i-2):(3*i))), 1:(length(OneLvlISCLTY)/3));


s3 = scatter(TwoLvlISX,TwoLvlISY,'b', 'filled');
s4 = scatter(OneLvlISCLTX,OneLvlISCLTY,'m', 'filled');

axis([xMin,xMax,yMin,yMax]);

r1 = refline(0,trueAns);
r1.Color = 'r';
r1.LineStyle = '--';
r2 = refline(0,trueAnsCLT);
r2.Color = 'm';
r2.LineStyle = '--';

legend('2LvlIS','1LvlISCLT','True Ans', 'True Ans CLT');
%title('S:10, l:0.20')
xlabel('Runtime (Seconds)');
ylabel('P(L > l)');


hold off
