function obj_value = mvlrecm_obj(weight, mnU, mvU, mveudis2, S, alpha, delta2, theta, eta, rho)  % objective function
%Calculate the objective function
%   weight: weight of views
%   mnU: 1*data_n cell, 2^cluster_n * view_n
%   mvU: 1*view_n cell, 2^cluster_n * data_n (U{m}(i,j) is the MF value of data j in cluster i in mth view.)
%   EUDIS2: size (2^cluster_n*data_n)
%   S: 2^CLUSTER_N * CLUSTER_N, S_ij=1 if class_i in cluster_j, 
%                               otherwise, S_ij=0.
%   CLUSTER_N: number of class, all 2^CLUSTER_N clusters .
%   APLHA:exponent to control the degree of penalization of the cardinal A.  
%   DELTA2: the parameter control the amount of data considered as outliers
%   THETA, RHO: the parameter control the influence of Low Rand Constraints
%   ETA: the parameter control the influence of weight entropy
data_n = size(mvU{1}, 2); 
c=sum(S,2);
c_rep=repmat((c(1:end-1,:).^alpha),1,data_n);
obj_value_part1=0; % W*U
for m=1:length(weight)
    value1=c_rep .* (mvU{m}(1:end-1,:).^2) .* mveudis2{m}(1:end-1,:);
    value2=delta2*sum(mvU{m}(end,:).^2);
    obj_value_part1=obj_value_part1+weight(m)*(sum(sum(value1))+value2);
end
obj_value_part2=0; %LowRank(U)
for i=1:data_n;
    obj_value_part2=obj_value_part2+theta*rho*rank(mnU{i});
end

obj_value_part3=0; %weight entropy
for m=1:length(weight)
    obj_value_part3=obj_value_part3+eta*(weight(m)*log(weight(m)));
end

obj_value=obj_value_part1+obj_value_part2+obj_value_part3;
end