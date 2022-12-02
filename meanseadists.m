clear all; close all;
load SSPrises
SSP19=(SSP19-6)/100; %minus current mean sea level (6cm ) and into m (/100)
SSP26=(SSP26-6)/100; 
SSP45=(SSP45-6)/100; 
SSP7=(SSP7-6)/100; 
SSP85=(SSP85-6)/100; 
SSP26low=(SSP26low-6)/100; 
SSP85low=(SSP85low-6)/100; 

global target
for YY=1:13
target=SSP19(YY,:);
startingparms =[ target(3)  (target(5)-target(1))/3.3   0.219171746062695];  %first guesses for distribution parameters (location,scale and shape) 
% Use fminsearch to find the best parameters; the error function is below.
% Note that the target values you are searching for are specified inside the error function.
finalparms = fminsearch(@SkewNormErr,startingparms,optimset('MaxFunEvals',1e4,'Maxiter',1e4));
% Create a distribution with the final best parameters that fminsearch could find:
eval(['skn19_',int2str(YY),' = SkewNor(finalparms(1),finalparms(2),finalparms(3));'])
% Check that it gives the desired mean & percentile values.
% These should match the "Known" values in the error function.
eval(['diffSN19(YY,1)=skn19_',int2str(YY),'.InverseCDF(.05)-target(1);']);
eval(['diffSN19(YY,2)=skn19_',int2str(YY),'.InverseCDF(.17)-target(2);']);
eval(['diffSN19(YY,3)=skn19_',int2str(YY),'.InverseCDF(.50)-target(3);']);
eval(['diffSN19(YY,4)=skn19_',int2str(YY),'.InverseCDF(.83)-target(4);']);
eval(['diffSN19(YY,5)=skn19_',int2str(YY),'.InverseCDF(.95)-target(5);']);
end

for YY=1:13
target=SSP26(YY,:);
startingparms =[ target(3)  (target(5)-target(1))/3.3   0.219171746062695];
% Use fminsearch to find the best parameters; the error function is below.
% Note that the target values you are searching for are specified inside the error function.
finalparms = fminsearch(@SkewNormErr,startingparms,optimset('MaxFunEvals',1e4,'Maxiter',1e4));
% Create a distribution with the final best parameters that fminsearch could find:
eval(['skn26_',int2str(YY),' = SkewNor(finalparms(1),finalparms(2),finalparms(3));'])
% Check that it gives the desired mean & percentile values.
% These should match the "Known" values in the error function.
eval(['diffSN26(YY,1)=skn26_',int2str(YY),'.InverseCDF(.05)-target(1);']);
eval(['diffSN26(YY,2)=skn26_',int2str(YY),'.InverseCDF(.17)-target(2);']);
eval(['diffSN26(YY,3)=skn26_',int2str(YY),'.InverseCDF(.50)-target(3);']);
eval(['diffSN26(YY,4)=skn26_',int2str(YY),'.InverseCDF(.83)-target(4);']);
eval(['diffSN26(YY,5)=skn26_',int2str(YY),'.InverseCDF(.95)-target(5);']);
end


for YY=1:13
target=SSP45(YY,:);
startingparms =[ target(3)  (target(5)-target(1))/3.3   0.219171746062695];
% Use fminsearch to find the best parameters; the error function is below.
% Note that the target values you are searching for are specified inside the error function.
finalparms = fminsearch(@SkewNormErr,startingparms,optimset('MaxFunEvals',1e5,'Maxiter',1e5));
% Create a distribution with the final best parameters that fminsearch could find:
eval(['skn45_',int2str(YY),' = SkewNor(finalparms(1),finalparms(2),finalparms(3));'])
% Check that it gives the desired mean & percentile values.
% These should match the "Known" values in the error function.
eval(['diffSN45(YY,1)=skn45_',int2str(YY),'.InverseCDF(.05)-target(1);']);
eval(['diffSN45(YY,2)=skn45_',int2str(YY),'.InverseCDF(.17)-target(2);']);
eval(['diffSN45(YY,3)=skn45_',int2str(YY),'.InverseCDF(.50)-target(3);']);
eval(['diffSN45(YY,4)=skn45_',int2str(YY),'.InverseCDF(.83)-target(4);']);
eval(['diffSN45(YY,5)=skn45_',int2str(YY),'.InverseCDF(.95)-target(5);']);
end

for YY=1:13
target=SSP7(YY,:);
startingparms =[ target(3)  (target(5)-target(1))/3.3   0.219171746062695];
% Use fminsearch to find the best parameters; the error function is below.
% Note that the target values you are searching for are specified inside the error function.
finalparms = fminsearch(@SkewNormErr,startingparms,optimset('MaxFunEvals',1e5,'Maxiter',1e5));
% Create a distribution with the final best parameters that fminsearch could find:
eval(['skn7_',int2str(YY),' = SkewNor(finalparms(1),finalparms(2),finalparms(3));'])
% Check that it gives the desired mean & percentile values.
% These should match the "Known" values in the error function.
eval(['diffSN7(YY,1)=skn7_',int2str(YY),'.InverseCDF(.05)-target(1);']);
eval(['diffSN7(YY,2)=skn7_',int2str(YY),'.InverseCDF(.17)-target(2);']);
eval(['diffSN7(YY,3)=skn7_',int2str(YY),'.InverseCDF(.50)-target(3);']);
eval(['diffSN7(YY,4)=skn7_',int2str(YY),'.InverseCDF(.83)-target(4);']);
eval(['diffSN7(YY,5)=skn7_',int2str(YY),'.InverseCDF(.95)-target(5);']);
end

for YY=1:13
target=SSP85(YY,:);
startingparms =[ target(3)  (target(5)-target(1))/3.3   0.219171746062695];
% Use fminsearch to find the best parameters; the error function is below.
% Note that the target values you are searching for are specified inside the error function.
finalparms = fminsearch(@SkewNormErr,startingparms,optimset('MaxFunEvals',1e5,'Maxiter',1e5));
% Create a distribution with the final best parameters that fminsearch could find:
eval(['skn85_',int2str(YY),' = SkewNor(finalparms(1),finalparms(2),finalparms(3));'])
% Check that it gives the desired mean & percentile values.
% These should match the "Known" values in the error function.
eval(['diffSN85(YY,1)=skn85_',int2str(YY),'.InverseCDF(.05)-target(1);']);
eval(['diffSN85(YY,2)=skn85_',int2str(YY),'.InverseCDF(.17)-target(2);']);
eval(['diffSN85(YY,3)=skn85_',int2str(YY),'.InverseCDF(.50)-target(3);']);
eval(['diffSN85(YY,4)=skn85_',int2str(YY),'.InverseCDF(.83)-target(4);']);
eval(['diffSN85(YY,5)=skn85_',int2str(YY),'.InverseCDF(.95)-target(5);']);
end

for YY=1:13
target=SSP85low(YY,:);
startingparms =[ target(3)/2  (target(5)-target(1))/3.3   1.219171746062695];
% Use fminsearch to find the best parameters; the error function is below.
% Note that the target values you are searching for are specified inside the error function.
if YY<=7   %A fix for the year 2130-2150 where the fit was not good
    finalparms = fminsearch(@SkewNormErr,startingparms,optimset('MaxFunEvals',1e5,'Maxiter',1e5));
else
    finalparms = fminsearch(@SkewNormErr_3t,startingparms,optimset('MaxFunEvals',1e5,'Maxiter',1e5));
end
% Create a distribution with the final best parameters that fminsearch could find:
eval(['skn85low_',int2str(YY),' = SkewNor(finalparms(1),finalparms(2),finalparms(3));'])
% Check that it gives the desired mean & percentile values.
% These should match the "Known" values in the error function.
eval(['diffSN85low(YY,1)=skn85low_',int2str(YY),'.InverseCDF(.05)-target(1);']);
eval(['diffSN85low(YY,2)=skn85low_',int2str(YY),'.InverseCDF(.17)-target(2);']);
eval(['diffSN85low(YY,3)=skn85low_',int2str(YY),'.InverseCDF(.50)-target(3);']);
eval(['diffSN85low(YY,4)=skn85low_',int2str(YY),'.InverseCDF(.83)-target(4);']);
eval(['diffSN85low(YY,5)=skn85low_',int2str(YY),'.InverseCDF(.95)-target(5);']);
end

for YY=1:13
target=SSP26low(YY,:);
startingparms =[ target(3)  (target(5)-target(1))/3.3   0.219171746062695];
% Use fminsearch to find the best parameters; the error function is below.
% Note that the target values you are searching for are specified inside the error function.
finalparms = fminsearch(@SkewNormErr,startingparms,optimset('MaxFunEvals',1e5,'Maxiter',1e5));
% Create a distribution with the final best parameters that fminsearch could find:
eval(['skn26low_',int2str(YY),' = SkewNor(finalparms(1),finalparms(2),finalparms(3));'])
% Check that it gives the desired mean & percentile values.
% These should match the "Known" values in the error function.
eval(['diffSN26low(YY,1)=skn26low_',int2str(YY),'.InverseCDF(.05)-target(1);']);
eval(['diffSN26low(YY,2)=skn26low_',int2str(YY),'.InverseCDF(.17)-target(2);']);
eval(['diffSN26low(YY,3)=skn26low_',int2str(YY),'.InverseCDF(.50)-target(3);']);
eval(['diffSN26low(YY,4)=skn26low_',int2str(YY),'.InverseCDF(.83)-target(4);']);
eval(['diffSN26low(YY,5)=skn26low_',int2str(YY),'.InverseCDF(.95)-target(5);']);
end
