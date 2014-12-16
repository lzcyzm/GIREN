% function for gene interaction network regularized elastic net algorithm
% Reference:
% Input:
%       x:  input miroarray matrix, rows are genes while cols are samples
%           size of G*S
%       y:  output response, clinical outcome or subtype class label,
%           length of S
%       ppi_mat:    protein-protein interaction network matrix
%       parameters: initialized parameters
%       beta0,beta: initialized parameters
%       alpha,gamma,lambda: initialized parameters
%       crit:    criticized threshold for giren algorithm
%       num_its: maximized number of iterations
% Output:
%       beta,beta0

function [beta0,beta,J,beta_rec]=giren(x_train,y_train,ppi_mat,parameters)

    beta0 = parameters.beta0;
    beta = parameters.beta;
    alpha = parameters.alpha;
    gamma = parameters.gamma;
    lambda = parameters.lambda;
    crit0 = parameters.crit;
    num_its0 = parameters.num_its;
    crit = 1;
    num_its = 0;
    
    [G,S] = size(x_train);
    while (crit > crit0 || num_its <= num_its0)
    		num_its = num_its + 1;
    		obj_ini = .5/S*(y_train-beta0-x_train'*beta)'*(y_train-beta0-x_train'*beta)+lambda*(.5*(1-alpha)*(beta'*beta)+alpha*sum(abs(beta)))-gamma*(beta'*ppi_mat*beta);
    		% ¸üĞÂbeta
            for j = 1:G
                e = 2*sum(x_train(j,:)*(y_train-beta0));
                a = 2*sum(sum(x_train.^2)) - gamma*sum(ppi_mat(:,j)) + lambda*(1-alpha);
                if lambda*alpha >= abs(e)
                    beta_new(j) = 0;
                    %disp('beta value set to 0!');
                else
                    if a > 0
                        if e > 0
                            beta_new(j) = (e - lambda*alpha)/a;
                        else
                            beta_new(j) = (e + lambda*alpha)/a;
                        end
                    else
                        if e > 0
                            beta_new(j) = (e + lambda*alpha)/a;
                        else
                            beta_new(j) = (e - lambda*alpha)/a;
                        end
                    end
                end
            end
            beta = beta_new';
    		beta0 = mean(y_train-x_train'*beta);
    		obj_now = .5/S*sum((y_train-beta0-x_train'*beta).^2)+lambda*(sum(.5*(1-alpha)*(beta.^2)+alpha*abs(beta)))-gamma*(beta'*ppi_mat*beta);
    		crit = abs((obj_now - obj_ini)/(obj_ini));
            J(num_its) = obj_now;
            % its_rec(num_its,:) = [obj_ini obj_now crit];
            beta_rec(num_its,:) = [beta0;beta];
    end
return