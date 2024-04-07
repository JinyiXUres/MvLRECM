%将存储标签顺序的空间矩阵转换为一个行向量
function [assignment] = MarkReplace(MarkMat)

[rows,cols]=size(MarkMat);

assignment=zeros(1,cols);

for i=1:rows
    for j=1:cols
        if MarkMat(i,j)==1
            assignment(1,j)=i;
        end
    end
end
end
