function [NewLabel, NewG] = BestMapping(La1,La2)

%La1: truelabel
%La2: cluster result (data_n*cluster_n) or data_n*1
%NewLabel: result
%G: right pairs matrix

cluster_n=length(unique(La1));

preLa2=La2;
if size(La2,2)==1 && cluster_n>=length(unique(La2)) % normal MVC
    proLa2=zeros(size(La2,1),cluster_n);
    for i=1:size(La2,1)
        proLa2(i,La2(i))=1;
    end
    La2=proLa2;
% elseif cluster_n<size(La2,2)
%     S=zeros(2^cluster_n,cluster_n); % S_ij=1 if w_k \in A_j
%     for i=1:cluster_n
%         comb=nchoosek(1:cluster_n,i);
%         zero_index=find(sum(S,2)==0);
%         first_zero_index=zero_index(1);
%         for j=1:size(comb,1)
%             S(first_zero_index-1+j,comb(j,:))=1;
%         end  
%     end
%     proLa2=zeros(size(La2,1),cluster_n);
%     for i=1:size(La2,1)
%         index1=find(preLa2(i,:)==1);
%         index2=find(S(index1,:)==1);
%         proLa2(i,index2)=1/length(index2);
%     end
%     La2=proLa2;
end


Label1=unique(La1);
L1=length(Label1);

L2=size(La2,2);
Label2=1:L2;


%构建计算两种分类标签重复度的矩阵G
G = zeros(max(L1,L2),max(L1,L2));
for i=1:L1
    index1= La1==Label1(i);
    for j=1:L2
        index2= La2(:,Label2(j))>0;
        if ~isempty(find(index1.*index2==1))
            G(i,j)=sum(La2(find(index1.*index2==1),Label2(j)));
        else
            G(i,j)=0;
        end
    end
end

%利用匈牙利算法计算出映射重排后的矩阵
[index]=munkres(-G);
%将映射重排结果转换为一个存储有映射重排后标签顺序的行向量
[temp]=MarkReplace(index);
%生成映射重排后的标签NewLabel
NewLabel=zeros(size(La2));
for i=1:L2
    NewLabel(:,temp(i))=La2(:,i);
end
NewG=zeros(size(G));
for i=1:L2
    NewG(:,temp(i))=G(:,i);
end

if size(preLa2,2)==1
    proLa2=zeros(size(preLa2,1),1);
    for i=1:size(preLa2,1)
        proLa2(i)=find(NewLabel(i,:)==1);
    end
    NewLabel=proLa2;
end
% 
% NewSinLabel=zeros(size(La1));
% if cluster_n<size(preLa2,2)
%     NewLabel=NewLabel>0;
%     NewSinLabel=sum(NewLabel,2);
%     index=find(NewSinLabel>1);
%     for i=index'
%         NewSinLabel(i)=find((sum(repmat(NewLabel(i,:),2^cluster_n,1).*S,2)==sum(NewLabel(i,:))).*(sum(NewLabel(i,:))==sum(S,2)));
%     end
% end
% 
end
