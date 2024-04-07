function [mnZ]=mvrlecm_Z(data, mvU, rho)
view_n = length(data);
data_n=size(data{1},1);
%mvU: 1*view_n cell 2^cluster_n * data_n
mnU=cell(1,data_n);
%1*data_n cell (2^cluster_n * view_n) 
mnZ=cell(1,data_n);
for i=1:data_n
    for m=1:view_n
        mnU{i}(:,m)=mvU{m}(:,i);
    end
    [U,S,V] = svd(mnU{i});
    zS=S-rho/2;
    zS(find(zS<0))=0;
    Z=U*zS*V';
    mnZ{i}=Z;
end
end