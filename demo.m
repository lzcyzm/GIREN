clear all;close all;clc;

load toy_data;
% k-fold cross validation
N = 3;
%for alpha_index = 1:9
%    for lambda_index = 1:9
for j = 1:9
indices = crossvalind('Kfold',y,N);
for i = 1:N
% call giren.m
    test = (indices == i); train = ~test;
    x_train = drug_fc(:,train);
    y_train = y(train);
    x_test = drug_fc(:,test);
    y_test = y(test);

    [G S] = size(x_train);
    temp = rand(G,1);
    parameters.beta = temp/sum(temp);
    parameters.beta0 = mean(y_train-x_train'*parameters.beta);
    parameters.alpha = 0.4;%alpha_id*0.1;             % elastic netÖÐalphaÎª£¨0,1£©
    parameters.lambda = 0.4;%lambda_id*0.1;             % penalty factor of L1 and L2 norm
    parameters.gamma =  parameters.lambda * 0.2;%gamma_id*0.1;   % penalty factor of gene-gene interaction network, gamma/lambda < 0.5
    parameters.crit = 1e-6;
    parameters.num_its = 2000;

    % GIREN method
    [beta0,beta,J,beta_rec] = giren(x_train,y_train,ppi_network,parameters);
    y_test_hat = x_test'*beta;
    y_test_giren = exp(y_test_hat)./(1+exp(y_test_hat));
    [X_giren,Y_giren,T_giren,AUC_giren] = perfcurve(y_test,y_test_giren,1);
    
    % elastic net method
    [b fitinfo] = lasso(x_train',y_train);
    y_test_hat = x_test'*b(:,end);
    y_test_elasticnet = exp(y_test_hat)./(1+exp(y_test_hat));
    [X_elasticnet,Y_elasticnet,T_elasticnet,AUC_elasticnet] = perfcurve(y_test,y_test_elasticnet,1);
   
    file_name = strcat('result-',num2str(i*10+j));
    save(strcat(file_name,'.mat'),'J','beta_rec','AUC_giren','b','AUC_elasticnet');
end
end
%    end
%end

