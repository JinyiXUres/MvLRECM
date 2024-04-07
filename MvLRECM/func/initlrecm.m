function [mvU, weight]= initecm(cluster_n, data_n, view_n)
%INITFCM Generate initial evidential partition matrix for ECM clustering.
%   U = INITECM(CLUSTER_N, DATA_N, VIEW_N) randomly generates VIEW_N evidential partition
%   matrix U, each U is 2^CLUSTER_N by DATA_N, where 2^CLUSTER_N is number of
%   clusters and DATA_N is number of data points. The summation of each
%   column of the generated U is equal to unity, as required by ecm clustering.
mvU=cell(1,view_n);
for i=1:view_n
    U_sing = rand(cluster_n, data_n);
    col_sum = sum(U_sing);
    U_sing = U_sing./col_sum(ones(cluster_n, 1), :);
    U=zeros(2^cluster_n, data_n);
    U(1:size(U_sing,1), 1:size(U_sing,2))=U_sing;
    mvU{i}=U;
end
weight=ones(1,view_n)./(view_n);
end
