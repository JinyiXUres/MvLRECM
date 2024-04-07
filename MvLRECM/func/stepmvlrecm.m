function [mvU_new, mvcenter, weight, obj_value] = stepmvlrecm(data, mvU, weight, S, cluster_n, alpha, delta2, theta, eta, rho, FOC)
%STEPFCM One step in fuzzy c-mean clustering.
%   [mvU_NEW, mvCENTER, obj_value] = STEPFCM(DATA, mvU, CLUSTER_N, ALPHA, DELTA2)
%   performs one iteration of Multi-view LowRank evidential version of fuzzy c-mean clustering, where
%
%   DATA: 1*view_n cell, data to be clustered. (data_n * dimen_n)
%   mvU: 1*view_n cell, 2^cluster_n by data_n (U{m}(i,j) is the MF value of data j in cluster i in mth view.)
%   weight: weight for views sum(weight)=1
%   CLUSTER_N: number of class, all 2^CLUSTER_N clusters .
%   APLHA:exponent to control the degree of penalization of the cardinal A.  
%   DELTA2: the parameter control the amount of data considered as outliers
%   THETA: the parameter control the influence of Low Rand Constraints
%   ETA, RHO: the parameter control the influence of weight entropy

%   mvU_NEW: new cell.
%   mvCENTER: 1*view_n cell, centers of clusters. (2^cluster_n * dimen_n)
%   obj_value: objective function after update.

mvcenter=mvlrecm_center(data, mvU, S, alpha);

[mvU, mveudis2]=mvlrecm_U(data, mvcenter, weight, S, alpha, delta2, theta, FOC);

weight=mvlrecm_weight(data, mvU, mveudis2, S, alpha, delta2, eta);

mvU_new=mvU;

[mnZ]=mvrlecm_Z(data, mvU, rho);

data_n=size(data{1},1);
view_n = length(data);
%Z: 1*data_n cell (2^cluster_n * view_n) 

mvU_new=cell(1,view_n);
%1*view_n cell 2^cluster_n * data_n
mnU=cell(1,data_n);
for i=1:data_n
    for m=1:view_n
        mnU{i}(:,m)=mnZ{i}(:,m)./sum(mnZ{i}(:,m));
        mvU_new{m}(:,i)=mnU{i}(:,m);
    end
end

obj_value = mvlrecm_obj(weight, mnU, mvU_new, mveudis2, S, alpha, delta2, theta, eta, rho);  % objective function
% obj_value=0;
end
