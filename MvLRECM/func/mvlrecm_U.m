function [mvU, mveudis2]=mvlrecm_U(data, mvcenter, weight, S, alpha, delta2, theta, FOC)
%Update Clustering Center
% calculate the Euclidean Distance
% upload U
%   weight: weight of views
%   DATA: matrix of data to be clustered. (Each row is a data point.)
%   mvCENTER: center of clusters in each view. (Each row is a center.)
%   S: 2^CLUSTER_N * CLUSTER_N, S_ij=1 if class_i in cluster_j, 
%                               otherwise, S_ij=0.
%   APLHA: exponent to control the degree of penalization of the cardinal A. 
%   mvU: 1*view_n cell, 2^cluster_n * data_n

cluster_n = size(S,2);
view_n=length(data);
data_n = size(data{1}, 1); 


c=sum(S,2);
c_rep=repmat((c(1:end-1,:).^alpha),1,data_n);

mveudis2=cell(1,view_n);
%Phi=zeros(1,data_n); % ¡ø
mvU=cell(1,view_n);
for m=1:view_n
    eudis2=zeros(2^cluster_n,data_n);
    for i=1:data_n
        for j=1:2^cluster_n
            eudis2(j,i)=EuDist2(mvcenter{m}(j,:),data{m}(i,:),0);
        end
    end
    mveudis2{m}=eudis2;
% U(j,i)=(c(j)^(alpha*betafun)*eudis2(j,i)^betafun) ...
% /(sum((c(1:end-1).^(alpha*betafun)).*(eudis2(1:end-1,i).^betafun))+delta2^betafun);
    
%     value1=weight(m)*c_rep .* eudis2(1:end-1,:)+theta;
%     value2=weight(m)*delta2+theta;   
    value1=c_rep .* eudis2(1:end-1,:)+theta;
    value2=delta2+theta;   
    Phi=(sum(value1.^-1)+value2.^-1).^-1;
    
    U=zeros(2^cluster_n, data_n);
    U(1:end-1,:)=repmat(Phi,2^cluster_n-1,1)./value1;
    
    if FOC==1
        for n=1:data_n
        thred=1-rand(1)*0.01-1/2^cluster_n;    
        sumU=sum(U(1:end-1,n));
            if sumU<thred
                U(1:end-1,n)=U(1:end-1,n).*(thred/sumU);
            end
        end
    end
    U(end,:)=1-sum(U(1:end-1,:)); 
%     
    mvU{m}=U;
end
end