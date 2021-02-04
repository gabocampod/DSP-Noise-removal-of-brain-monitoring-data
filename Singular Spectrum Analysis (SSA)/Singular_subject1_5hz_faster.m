[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
%%
%%STEP 1: set parameters and define input X
load('Test_1mA_5Hz_ERP_SUBJ1.mat')
RAW_SUBJECT_DATA = EEG.Data(1,:);

subject1_tacs_t1 = 15025;
subject1_tacs_t2 = subject1_tacs_t1 + 29999 ;

tic;
X_init = RAW_SUBJECT_DATA(1,subject1_tacs_t1: subject1_tacs_t2);  %FOR ORIGINAL TEST 
X = X_init - mean(X_init);
M = 200;    % window length = embedding dimension
N = length(X);
t = (1:N)';
%%
%%step 2: Calculate Covariance matrix C

Y=zeros(N-M+1,M);

for m=1:M
  Y(:,m) = X((1:N-M+1)+m-1);
end;
Cemb=Y'*Y / (N-M+1);
C = Cemb;
%%
%%STEP 3: Calculate eigenvalues (Lambda) and eigen vectors (RHO)

[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors

%%
%%STEP 4: Calculate prncipal components (PC) PC = Y*RHO;
PC = Y*RHO;

%%
%%STEP 5: Calculate reconstructed components  RC = PC*RHO
RC=zeros(N,M);
for m=1:2
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;

summed = sum(RC(:,1:2),2);
free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));
toc;
%%
%%STEP 6: Compare reconstructed to original
load('Test_Sham_ERP_SUBJ1.mat') 

SHAM_ORIGINAL_BEST=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, 15015:45014)));


free_single_5HZ_1500 = free_filtered;

figure(7);  %COMPARE SHAM TO ARTIFACT FREE SIGNAL 
plot(SHAM_ORIGINAL_BEST)
hold on
plot(free_filtered)
hold off
title('COMPARE SHAM TO ARTIFACT FREE SIGNAL ')

R = corrcoef(SHAM_ORIGINAL_BEST(1, 5000:25000), free_filtered(1, 5000:25000));
similarity_best = R(1,2);





