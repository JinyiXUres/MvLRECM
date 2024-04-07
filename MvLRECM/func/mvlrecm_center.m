function [ mvcenter ] = mvlrecm_center( data, mvU, S, alpha )
%Update Clustering Center
% upload class (singleton cluster) center first
% calculate the meta-cluster centers
%   DATA: matrix of data to be clustered. (Each row is a data point.)
%   U: partition matrix. (U(i,j) is the MF value of data j in cluster i.)
%   S: 2^CLUSTER_N * CLUSTER_N, S_ij=1 if class_i in cluster_j, 
%                               otherwise, S_ij=0.
%   CLUSTER_N: number of class, 2^CLUSTER_N clusters .
%   BETA: exponent (> 1) for the partition matrix.
%   APLHA:exponent to control the degree of penalization of the cardinal A.  
%   B: a matrix of size(CLUSTER_N,IN_N)
%   H: a matrix of size(CLUSTER_N,CLUSTER_N)
%   HV=B, V is the singleton cluster center (Each row is a center.)
%   CENTER: center of clusters. (Each row is a center.)

data_n = size(data{1}, 1); 
cluster_n = size(S,2);
view_n=length(data);
mvcenter=cell(1,view_n);
for m=1:view_n

    in_n = size(data{m}, 2);
    U=mvU{m};
    x=data{m};
%     B=zeros(cluster_n,in_n);
%     for l=1:cluster_n
%         for q=1:in_n
%             index_c=find(S(:,l)==1);
%             c=sum(S(index_c,:),2);
%             c_rep=repmat((c.^(alpha-1)),1,data_n);
%             B(l,q)=sum(x(:,q)'.*sum(c_rep.*(U(index_c,:).^2)));
%         end
%     end
%     H=zeros(cluster_n,cluster_n);
%     for l=1:cluster_n
%         for k=1:cluster_n
%             index_c=find(sum(S(:,[l,k]),2)==2);
%             c=sum(S(index_c,:),2);
%             H(l,k)=sum(sum(U(index_c,:).^2,2).*(c.^(alpha-2)));
%         end
%     end
%     V=H^-1*B;
%     center=zeros(2^cluster_n,in_n);
%     center(1:cluster_n,:)=V;
%     for i=cluster_n+1:2^cluster_n
%         center(i,:)=sum(V(find(S(i,:)==1),:))/sum(S(i,:));
%     end
    [ center ] = ecm_center( x, U, S, cluster_n, 2, alpha );
    mvcenter{m}=center;
end
end

