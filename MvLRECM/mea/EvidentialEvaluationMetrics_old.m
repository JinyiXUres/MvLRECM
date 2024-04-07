function [res, BestMapY]=  EvidentialEvaluationMetrics_old(y0, y)
%res=[ACC, NMI, purity, F1score, Precision, Recall, RI, Uncertain]
if size(y0,2) ~= 1
    y0 = y0';
end;
data_n=length(y0);
Uy = unique(y0);
cluster_n=length(Uy);
proY0 = zeros(data_n,1);
if cluster_n ~= max(y)
    for i = 1:cluster_n
        proY0(find(y0 == Uy(i))) = i;
    end;
    y0 = proY0;
end;


if size(y,2)==data_n
    y = y';
end;

proY=zeros(data_n,cluster_n);
if size(y,2)==1  % normal MVC
    for i=1:data_n
        proY(i,y(i))=1;
    end
elseif size(y,2)==cluster_n
    proY=y;
elseif size(y,2)>cluster_n
    S=zeros(2^cluster_n,cluster_n); % S_ij=1 if w_k \in A_j
    T=[];
    for i=1:cluster_n
        comb=nchoosek(1:cluster_n,i);
        zero_index=find(sum(S,2)==0);
        first_zero_index=zero_index(1);
        for j=1:size(comb,1)
            S(first_zero_index-1+j,comb(j,:))=1;
        end  
        T(i)=first_zero_index-1+j;
    end
    t=find(T==size(y,2));
    if size(y,2)==2^cluster_n || t
        for i=1:data_n
            index1=find(y(i,:)==1);
            index2=find(S(index1(1),:)==1);
            proY(i,index2)=1/length(index2);
        end
%     elseif t
%         S=S(1:T(t),:);
%         for i=1:data_n
%             index1=find(y(i,:)==1);
%             index2=find(S(index1,:)==1);
%             proY(i,index2)=1/length(index2);
%         end
    else
        fprintf('Y is Wrong£¡')
    end
else
    fprintf('Y is Wrong£¡')
end

[BestMapY, G] = BestMapping(y0,proY);

%% Purity
profullY=proY>0;
if size(y,2)==1 
    purity=sum(max(G))/data_n;
else
    for i=1:log2(size(y,2))
        index1=find(profullY(:,i)==1);
        for j=1:cluster_n
            fullG(i,j)=sum(y0(index1)==j);
        end
    end
    purity=sum(max(fullG))/data_n;
end


%% ACC
proY0=zeros(data_n,cluster_n);
for i=1:data_n
    proY0(i,y0(i))=1;
end
ACC=sum(sum(proY0.*(BestMapY>0)))/data_n;

%% Precision, Recall, F1score
TP=0;
FP=0;
TN=0;
FN=0;
for i=1:data_n
    for j=i+1:data_n
        if y0(i)==y0(j)
            if sum(BestMapY(i,:).*BestMapY(j,:))~=0
                TP=TP+1;
            else
                FN=FN+1;
            end
        else
            if sum((BestMapY(i,:)+BestMapY(j,:))>0)==1 && sum(BestMapY(i,:).*BestMapY(j,:))~=0
                FP=FP+1;
            else
                TN=TN+1;
            end
        end
    end
end
%Precision
Precision=TP/(TP+FP);
%Recall
Recall=TP/(TP+FN);
% F1score
F1score=2*(Precision*Recall/(Precision+Recall));

%% Uncertain
Uncertain=sum(sum(y(:,cluster_n+1:end-1)))/data_n;

%% RI
RI=(TP+TN)/(TP+FP+TN+FN);

%% NMI
G = G+eps;
sumG = sum(G(:));
Py0 = sum(G,2)/sumG;
Hy0=0;
for i=1:numel(Py0)
    Hy0=Hy0-Py0(i).*log2(Py0(i));
end

Py = sum(G,1)/sumG;
Hy=0;
for i=1:numel(Py)
    Hy=Hy-Py(i).*log2(Py(i));
end

Py0y = G./repmat(sum(G,1),cluster_n,1);
Hy0y=0;
for i=1:size(Py0y,2)
    tep=0;
    for j=1:size(Py0y,1)
        tep=tep-Py0y(j,i).*log2(Py0y(j,i));
    end
    Hy0y=Hy0y+Py(i)*tep;
end
I=Hy0-Hy0y;
NMI = 2*I / (Hy0+Hy);

res=[ACC, NMI, purity, F1score, Precision, Recall, RI, Uncertain];
end

