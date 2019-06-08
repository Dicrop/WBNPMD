function [S_F]=WBNPMD(rFunctionalArray,dWeightArray,rWeightArray,dSemanticArray,interaction,alpha,beta)

%% 初始化miRNA-疾病关系矩阵
m=495;
n=383;
S_r=zeros(m,m);
S_d=zeros(n,n);

%% 高斯核相互作用
[GaussR,GaussD] = gaussiansimilarity(interaction,m,n);

for i=1:n
    for j=1:n
        if (dWeightArray(i,j) == 1) S_d(i,j)=dSemanticArray(i,j);
        else if (dWeightArray(i,j) == 0) S_d(i,j)=GaussD(i,j);
            end
        end
    end
end


for i=1:m
    for j=1:m
        if (rWeightArray(i,j) == 1) S_r(i,j)=rFunctionalArray(i,j);
        else if (rWeightArray(i,j) == 0) S_r(i,j)=GaussR(i,j);
            end
        end
    end
end

test=interaction';

%% 构建初始权重
 
y1=S_r*interaction./sum(S_r,2);%.*interaction; %miRNA权重计算

y2=S_d*test./sum(S_d,2); %疾病权重计算

%% 第一次 miRNA (邻接矩阵495 by 383)

W1=zeros(m,m);
ini1=y1';

for i=1:m
    gmq=bsxfun(@times,ini1,ini1(:,i));
    a=sum(ini1,2);
    q=bsxfun(@rdivide,gmq,a);
    q(find(isnan(q)==1))=0;
    W1(i,:)=sum(q)./sum(ini1);
    W1(find(isnan(W1)==1))=0;
end

% WT=W1+(-alpha)*W1*W1;
% S_M=W1*y1;

mod1=interaction.*(sum(interaction)).^(-beta);%初始分数测试
S_M=W1*mod1;

%% 第二次 疾病 （邻接矩阵383 by 495)

W2=zeros(n,n);
ini_2=y2';

for i=1:n
    gmq=bsxfun(@times,ini_2,ini_2(:,i));
    a=sum(ini_2,2);
    q=bsxfun(@rdivide,gmq,a);
    q(find(isnan(q)==1))=0;
    W2(i,:)=sum(q)./sum(ini_2);
    W2(find(isnan(W2)==1))=0;
end

% WT=W2+(-alpha)*W2*W2;
% S_D=W2*y2;

mod2=test.*(sum(test)).^(-beta);%初始分数测试
S_D=W2*mod2;

%% 计算最终得分
S_F=(S_M+S_D')/2;
% S_F=S_M;
% S_F=S_D';

toc
end
