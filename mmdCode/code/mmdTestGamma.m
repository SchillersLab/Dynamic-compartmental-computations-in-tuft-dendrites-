%This function implements the MMD two-sample test using a Gamma approximation
%to the test threshold

%Arthur Gretton
%07/12/08

%Inputs: 
%        X contains dx columns, m rows. Each row is an i.i.d sample
%        Y contains dy columns, m rows. Each row is an i.i.d sample
%        alpha is the level of the test
%        params.sig is kernel size. If -1, use median distance heuristic.


%Outputs: 
%        thresh: test threshold for level alpha test
%        testStat: test statistic m*MMD (biased)


function [testStat,thresh,params] = mmdTestGamma(X,Y,alpha,params);

    
m=size(X,1);



%Set kernel size to median distance between points IN AGGREGATE SAMPLE
if params.sig == -1
  Z = [X;Y];  %aggregate the sample
  size1=size(Z,1);
    if size1>100
      Zmed = Z(1:100,:);
      size1 = 100;
    else
      Zmed = Z;
    end
    G = sum((Zmed.*Zmed),2);
    Q = repmat(G,1,size1);
    R = repmat(G',size1,1);
    dists = Q + R - 2*Zmed*Zmed';
    dists = dists-tril(dists);
    dists=reshape(dists,size1^2,1);
    params.sig = sqrt(0.5*median(dists(dists>0)));  %rbf_dot has factor of two in kernel
end


%MMD statistic. Here we use biased 
%statistic and get Gamma fit.

K = rbf_dot(X,X,params.sig);
L = rbf_dot(Y,Y,params.sig);
KL = rbf_dot(X,Y,params.sig);


testStat = 1/m^2 * sum(sum(K + L - KL - KL'));
testStat = testStat * m; %null distirbution on m*MMD_b

%mean under H0 of MMD  
meanMMD  = 2/m * ( 1  - 1/m*sum(diag(KL))  );  

K = K - diag(diag(K));
L = L - diag(diag(L));
KL = KL - diag(diag(KL));
  
  
%Variance under H0 of MMD
varMMD = 2/m/(m-1) * 1/m/(m-1) * sum(sum( (K + L - KL - KL').^2 ));


al = meanMMD^2 / varMMD;
bet = varMMD*m / meanMMD;   

thresh = icdf('gam',1-alpha,al,bet);    %%%% TEST THRESHOLD
