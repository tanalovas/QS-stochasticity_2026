
function [bigM, bigB]=matSecMoms(M,B,n)

nMoments=n*(n+1)/2;
bigM=zeros(nMoments);
bigB=zeros(nMoments,1);

for K=1:n
    for L=K:n
        indexRow=(K-1)*n+1-(K-2)*(K-1)/2+L-K;
        for j=1:n
            if j>=L
                indexColumn1 = (L-1)*n+1-(L-2)*(L-1)/2+j-L;
            else
                indexColumn1 = (j-1)*n+1-(j-2)*(j-1)/2+L-j;
            end

            bigM(indexRow,indexColumn1)= bigM(indexRow,indexColumn1)+M(K,j);

            if j>=K
                indexColumn2 = (K-1)*n+1-(K-2)*(K-1)/2+j-K;
            else
                indexColumn2 = (j-1)*n+1-(j-2)*(j-1)/2+K-j;
            end

            bigM(indexRow,indexColumn2)= bigM(indexRow,indexColumn2)+M(L,j);
            
        end

         bigB(indexRow)=B(K,L);
        
    end
end