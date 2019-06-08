function AUC=calAUC(inputG,outputG)
%Alternative calculation method for AUC

%load input
%load output
%load auc_matrix2
%outputG=auc_matrix;

%[n,m]=size(inputG);
%ini=zeros(n,n);
%p=4796;

A=find(inputG(:,:)==1);
B=find(inputG(:,:)==0);

[n,~]=size(A); %size of association matrix
[m,~]=size(B); %size of non-association matrix

%A=A(randperm(numel(A))); %randomize asociation matrix
%B=B(randperm(numel(B))); %randomize non-association matrix

%n=p;m=p;

out=0;
for i=1:n
    i;
    for j=1:m
        if(outputG(A(i))>outputG(B(j)))
            out=out+1;
        end
        if(outputG(A(i))==outputG(B(j)))
            out=out+0.5;
        end
    end
end

AUC=out/(m*n)

end
 