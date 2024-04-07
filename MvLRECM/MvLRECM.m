function [index, sumU] = MvLRECM(data, y0, FOC, options )
%ECM Data set clustering using evidential verson of fuzzy c-means clustering.
%
%   [CENTER, U, OBJ_FCN] = ECM(DATA, N_CLUSTER) finds 2^N_CLUSTER number of
%   clusters in the data set DATA. DATA is size N-by-M, where N is the number of
%   data points and M is the number of coordinates for each data point. The
%   coordinates for each cluster center are returned in the rows of the matrix
%   CENTER. The membership function matrix U contains the grade of membership of
%   each DATA point in each cluster. The values 0 and 1 indicate no membership
%   and full membership respectively. Grades between 0 and 1 indicate that the
%   data point has partial membership in a cluster. At each iteration, an
%   objective function is minimized to find the best location for the clusters
%   and its values are returned in OBJ_FCN.
%


if nargin ~= 2 && nargin ~= 3 && nargin ~= 4,
	error('Too many or too few input arguments!');
end

% Change the following to set default options
default_options = [2;   % alpha, exponent to control the degree of penalization.  (default: 2.0)
                   20;  % delta^2, the parameter control the amount of data considered as outliers  (default: 20)
                   1; % theta, the parameter control the influence of Low Rand Constraints
                   1; % eta, the parameter control the influence of weight entropy
                   2^(-length(unique(y0))/2); %rho
                   100;	% max. number of iteration
                   1e-5;	% min. amount of improvement
                   1];	% info display during iteration 


if nargin == 2
    FOC=0;
end

if nargin < 4
    options=default_options;
else
	% If "options" is not fully specified, pad it with default values.
	if length(options) < 8,
		tmp = default_options;
		tmp(1:length(options)) = options;
		options = tmp;
	end
	% If some entries of "options" are nan's, replace them with defaults.
	nan_index = find(isnan(options)==1);
	options(nan_index) = default_options(nan_index);
end

alpha = options(1);		% Exponent for C
delta2 = options(2);		% Parameter for outliers
theta = options(3);
eta = options(4);
rho = options(5);
max_iter = options(6);		% Max. iteration
min_impro = options(7);		% Min. improvement
display = options(8);		% Display info or not

obj_value = zeros(max_iter, 1);	% Array for objective function

data_n=length(y0);
cluster_n=length(unique(y0));
view_n=length(data);
[mvU, weight]= initlrecm(cluster_n, data_n, view_n);			% Initial fuzzy partition

S=zeros(2^cluster_n,cluster_n); % S_ij=1 if w_k \in A_j
for i=1:cluster_n
    comb=nchoosek(1:cluster_n,i);
    zero_index=find(sum(S,2)==0);
    first_zero_index=zero_index(1);
    for j=1:size(comb,1)
        S(first_zero_index-1+j,comb(j,:))=1;
    end  
end

% Main loop
for i = 1:max_iter,
	[mvU, mvcenter, weight, obj_value(i)] = stepmvlrecm(data, mvU, weight, S, cluster_n, alpha, delta2, theta, eta, rho, FOC);
	if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_value(i));
	end
	% check termination condition
	if i > 1,
		if abs(obj_value(i) - obj_value(i-1)) < min_impro, break; end,
	end
end

iter_n = i;	% Actual number of iterations 
obj_value(iter_n+1:max_iter) = [];
sumU=0;

for view_i=1:view_n
    sumU=sumU+mvU{view_i}*weight(view_i);
end
maxU = max(sumU);
index=zeros(data_n,2^cluster_n);
for i=1:2^cluster_n
    x=find(sumU(i,:) == maxU);
    index(x,i)=1;  
end

for i=1:data_n
    if sum(index(i,:))>1
        xS=sum(S(find(index(i,:)>0),:));
        for s=1:2^cluster_n
            if all(xS==S(s,:))
                index(i,:)=0;
                index(i,s)=1;
                break
            end
        end
    end
end

end
