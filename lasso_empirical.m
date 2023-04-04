function [realdata_all] = lasso_empirical(ReachS,Bin10smooth)

% lasso regression model

% clip and concatenate reaches
ReachS(1).kin_cat_10ms = [];
for i=1:length(ReachS)
    if ReachS(i).exclude == 0 
        ReachS(1).kin_cat_10ms = vertcat(ReachS(1).kin_cat_10ms, ReachS(i).kin_10ms(300:500,:));
    end
end

% pull out predictors 
pos_x = ReachS(1).kin_cat_10ms(:,2);
pos_y = ReachS(1).kin_cat_10ms(:,3);
pos_z = ReachS(1).kin_cat_10ms(:,4);
vel = ReachS(1).kin_cat_10ms(:,5);

vel_x = ReachS(1).kin_cat_10ms(:,6);
vel_x_up_la = ReachS(1).kin_cat_10ms(:,6) >= 0;
vel_x_up = ReachS(1).kin_cat_10ms(:,6).*vel_x_up_la;
vel_x_down_la = ReachS(1).kin_cat_10ms(:,6) < 0;
vel_x_down = abs(ReachS(1).kin_cat_10ms(:,6).*vel_x_down_la);

clear vel_x_up_la
clear vel_x_down_la

vel_y = ReachS(1).kin_cat_10ms(:,7);
vel_y_up_la = ReachS(1).kin_cat_10ms(:,7) >= 0;
vel_y_up = ReachS(1).kin_cat_10ms(:,7).*vel_y_up_la;
vel_y_down_la = ReachS(1).kin_cat_10ms(:,7) < 0;
vel_y_down = abs(ReachS(1).kin_cat_10ms(:,7).*vel_y_down_la);

clear vel_y_up_la
clear vel_y_down_la

vel_z = ReachS(1).kin_cat_10ms(:,8);
vel_z_up_la = ReachS(1).kin_cat_10ms(:,8) >= 0;
vel_z_up = ReachS(1).kin_cat_10ms(:,8).*vel_z_up_la;
vel_z_down_la = ReachS(1).kin_cat_10ms(:,8) < 0;
vel_z_down = abs(ReachS(1).kin_cat_10ms(:,8).*vel_z_down_la);

clear vel_z_up_la
clear vel_z_down_la

acc = ReachS(1).kin_cat_10ms(:,9);

acc_x = ReachS(1).kin_cat_10ms(:,10);
acc_x_up_la = ReachS(1).kin_cat_10ms(:,10) >= 0;
acc_x_up = ReachS(1).kin_cat_10ms(:,10).*acc_x_up_la;
acc_x_down_la = ReachS(1).kin_cat_10ms(:,10) < 0;
acc_x_down = abs(ReachS(1).kin_cat_10ms(:,10).*acc_x_down_la);

clear acc_x_up_la
clear acc_x_down_la

acc_y = ReachS(1).kin_cat_10ms(:,11);
acc_y_up_la = ReachS(1).kin_cat_10ms(:,11) >= 0;
acc_y_up = ReachS(1).kin_cat_10ms(:,11).*acc_y_up_la;
acc_y_down_la = ReachS(1).kin_cat_10ms(:,11) < 0;
acc_y_down = abs(ReachS(1).kin_cat_10ms(:,11).*acc_y_down_la);

clear acc_y_up_la
clear acc_y_down_la

acc_z = ReachS(1).kin_cat_10ms(:,12);
acc_z_up_la = ReachS(1).kin_cat_10ms(:,12) >= 0;
acc_z_up = ReachS(1).kin_cat_10ms(:,12).*acc_z_up_la;
acc_z_down_la = ReachS(1).kin_cat_10ms(:,12) < 0;
acc_z_down = abs(ReachS(1).kin_cat_10ms(:,12).*acc_z_down_la);

clear acc_z_up_la
clear acc_z_down_la

Predictors = [pos_x pos_y pos_z vel vel_x vel_x_up vel_x_down vel_y vel_y_up vel_y_down vel_z vel_z_up vel_z_down acc acc_x acc_x_up acc_x_down acc_y acc_y_up acc_y_down acc_z acc_z_up acc_z_down];

clear pos_x pos_y pos_z vel vel_x vel_x_up vel_x_down vel_y vel_y_up vel_y_down vel_z vel_z_up vel_z_down acc acc_x acc_x_up acc_x_down acc_y acc_y_up acc_y_down acc_z acc_z_up acc_z_down

%% Regression

% find aligned FR start and end points
for i=1:length(ReachS)
    if ReachS(i).exclude == 0 
        [~,idx1(i,1)]=min(abs(Bin10smooth(:,1)-ReachS(i).kin_10ms(300,1)));
        idx2(i,1)=idx1(i,1)+200;
    end
end

% define lag and steps to use
lag_t = 30; % in 10's of ms
step_size = 1;
lag = lag_t/step_size;

% concatenate FRs at lag
a=1;
for ii=-lag:0
Responses = [];
for i=1:length(ReachS)
    if ReachS(i).exclude == 0
        Responses = vertcat(Responses, Bin10smooth(idx1(i,1)+(ii*step_size):idx2(i,1)+(ii*step_size),2));
    end
end

% regression
[B,FitInfo] = lasso(Predictors,Responses,'Standardize',true,'CV',10,'PredictorNames',{'pos_x' 'pos_y' 'pos_z' 'vel' 'vel_x' 'vel_x_up' 'vel_x_down' 'vel_y' 'vel_y_up' 'vel_y_down' 'vel_z' 'vel_z_up' 'vel_z_down' 'acc' 'acc_x' 'acc_x_up' 'acc_x_down' 'acc_y' 'acc_y_up' 'acc_y_down' 'acc_z' 'acc_z_up' 'acc_z_down'});

% put regression results into struct
Regression(a).B=B;
Regression(a).FitInfo=FitInfo;

% find lambda min that minimized MSE
[~,idx4] = min(Regression(a).FitInfo.MSE);

% get coefficients at min lag min lambda
X = Regression(a).B(:,idx4);
coef0 = Regression(a).FitInfo.Intercept(1,idx4);

yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
     Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
     Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
     Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
     Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
     Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
     Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
     Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;

% calculate R2 
ybar = mean(Responses);
SSresid = sum((Responses-yhat).^2);
SStotal = sum((Responses-ybar).^2);
rsq = 1 - (SSresid/SStotal);
RSQ(a) = rsq;

a=a+1;
end

%% Reporting 

% find lag that maximized R2
[rsq, idx] = max(RSQ);
lagmax = -lag+(idx-1);

% recompute values for that lag
[~,idx4] = min(Regression(idx).FitInfo.MSE);

X = Regression(idx).B(:,idx4);

coef0 = Regression(idx).FitInfo.Intercept(1,idx4);

yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
     Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
     Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
     Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
     Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
     Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
     Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
     Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;
 
ii = lagmax;
Responses = [];
for i=1:length(ReachS)
    if ReachS(i).exclude == 0
        Responses = vertcat(Responses, Bin10smooth(idx1(i,1)+(ii*step_size):idx2(i,1)+(ii*step_size),2));
    end
end

% report data
realdata_all.MSEmin.Lambda = Regression(idx).FitInfo.Lambda;
realdata_all.MSEmin.DF = Regression(idx).FitInfo.DF;
realdata_all.MSEmin.MSE = Regression(idx).FitInfo.MSE;
realdata_all.MSEmin.rsq = rsq;
realdata_all.MSEmin.lagmax = lagmax*10; % binned at 10 ms
realdata_all.MSEmin.predictors = X;
realdata_all.MSEmin.predictornames = Regression(1).FitInfo.PredictorNames; % doesn't matter which regression you get names from because they are all the same
realdata_all.MSEmin.predictedFR = yhat;
realdata_all.MSEmin.actualFR = Responses;


%% Get MSEs from regression with all at Lambda1SE

for i = 1:length(Regression)
    % find lambda min that minimized MSE
    idx4 = Regression(i).FitInfo.Index1SE;

    % get coefficients at min lag min lambda
    X = Regression(i).B(:,idx4);
    coef0 = Regression(i).FitInfo.Intercept(1,idx4);

    yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
         Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
         Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
         Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
         Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
         Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
         Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
         Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;

    % calculate R2 
    ybar = mean(Responses);
    SSresid = sum((Responses-yhat).^2);
    SStotal = sum((Responses-ybar).^2);
    rsq = 1 - (SSresid/SStotal);
    RSQ_1SE(i) = rsq;
end

%% Reporting 

% find lag that maximized R2
[rsq, idx] = max(RSQ_1SE);
lagmax = -lag+(idx-1);

% recompute values for that la
idx4 = Regression(idx).FitInfo.Index1SE;

X = Regression(idx).B(:,idx4);

coef0 = Regression(idx).FitInfo.Intercept(1,idx4);

yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
     Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
     Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
     Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
     Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
     Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
     Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
     Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;
 
ii = lagmax;
Responses = [];
for i=1:length(ReachS)
    if ReachS(i).exclude == 0
        Responses = vertcat(Responses, Bin10smooth(idx1(i,1)+(ii*step_size):idx2(i,1)+(ii*step_size),2));
    end
end

% report data
realdata_all.MSE1SE.Lambda = Regression(idx).FitInfo.Lambda;
realdata_all.MSE1SE.DF = Regression(idx).FitInfo.DF;
realdata_all.MSE1SE.MSE = Regression(idx).FitInfo.MSE;
realdata_all.MSE1SE.rsq = rsq;
realdata_all.MSE1SE.lagmax = lagmax*10; % binned at 10 ms
realdata_all.MSE1SE.predictors = X;
realdata_all.MSE1SE.predictornames = Regression(1).FitInfo.PredictorNames; % doesn't matter which regression you get names from because they are all the same
realdata_all.MSE1SE.predictedFR = yhat;
realdata_all.MSE1SE.actualFR = Responses;


%% Get MSEs from regression with lag at DF3

for i = 1:length(Regression)
    % find lambda that minimized MSE with 3 variables
    for ii = 1:length(Regression(i).FitInfo.DF)
        if Regression(i).FitInfo.DF(ii) == 3
            idx4 = ii;
            break
        end
    end

    % get coefficients at min lag min lambda
    X = Regression(i).B(:,idx4);
    coef0 = Regression(i).FitInfo.Intercept(1,idx4);

    yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
         Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
         Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
         Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
         Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
         Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
         Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
         Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;

    % calculate R2 
    ybar = mean(Responses);
    SSresid = sum((Responses-yhat).^2);
    SStotal = sum((Responses-ybar).^2);
    rsq = 1 - (SSresid/SStotal);
    RSQ_DF3(i) = rsq;
end

%% Reporting 

% find lag that maximized R2
[rsq, idx] = max(RSQ_DF3);
lagmax = -lag+(idx-1);

% recompute values for that la
% find lambda that minimized MSE with 3 variables
for ii = 1:length(Regression(idx).FitInfo.DF)
    if Regression(idx).FitInfo.DF == 3
        idx4 = ii;
        break
    end
end
    
X = Regression(idx).B(:,idx4);

coef0 = Regression(idx).FitInfo.Intercept(1,idx4);

yhat = Predictors(:,1)*X(1) + Predictors(:,2)*X(2) + Predictors(:,3)*X(3) + ...
     Predictors(:,4)*X(4) + Predictors(:,5)*X(5) + Predictors(:,6)*X(6) +...
     Predictors(:,7)*X(7) + Predictors(:,8)*X(8) + Predictors(:,9)*X(9) + ...
     Predictors(:,10)*X(10) + Predictors(:,11)*X(11) + Predictors(:,12)*X(12) + ...
     Predictors(:,13)*X(13) + Predictors(:,14)*X(14) + Predictors(:,15)*X(15) + ...
     Predictors(:,16)*X(16) + Predictors(:,17)*X(17) + Predictors(:,18)*X(18) + ...
     Predictors(:,19)*X(19) + Predictors(:,20)*X(20) + Predictors(:,21)*X(21) + ...
     Predictors(:,22)*X(22) + Predictors(:,23)*X(23) + coef0;
 
ii = lagmax;
Responses = [];
for i=1:length(ReachS)
    if ReachS(i).exclude == 0
        Responses = vertcat(Responses, Bin10smooth(idx1(i,1)+(ii*step_size):idx2(i,1)+(ii*step_size),2));
    end
end

% report data
realdata_all.MSEDF3.Lambda = Regression(idx).FitInfo.Lambda;
realdata_all.MSEDF3.DF = Regression(idx).FitInfo.DF;
realdata_all.MSEDF3.MSE = Regression(idx).FitInfo.MSE;
realdata_all.MSEDF3.rsq = rsq;
realdata_all.MSEDF3.lagmax = lagmax*10; % binned at 10 ms
realdata_all.MSEDF3.predictors = X;
realdata_all.MSEDF3.predictornames = Regression(1).FitInfo.PredictorNames ;% doesn't matter which regression you get names from because they are all the same
realdata_all.MSEDF3.predictedFR = yhat;
realdata_all.MSEDF3.actualFR = Responses;

