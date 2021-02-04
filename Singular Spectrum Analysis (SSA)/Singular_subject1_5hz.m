[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
%%
%%STEP 1: set parameters and define input X

load('Test_1mA_5Hz_ERP_SUBJ1.mat')
RAW_SUBJECT_DATA = EEG.Data(1,:);

subject1_tacs_t1 = 15025;
subject1_tacs_t2 = subject1_tacs_t1 + 29999 ;

X_init = RAW_SUBJECT_DATA(1,subject1_tacs_t1: subject1_tacs_t2);  %FOR ORIGINAL TEST 
X = X_init - mean(X_init);
M = 30;    % window length = embedding dimension
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

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:3
  subplot(2,2,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
end;

%%
%%STEP 6: Compare reconstructed to original
load('Test_Sham_ERP_SUBJ1.mat') 

SHAM_ORIGINAL_BEST=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, 15015:45014)));

summed = sum(RC(:,1:2),2);
free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));

figure(6);  %COMPARE ORIGINAL (x) TO reconstructed (sum of ALL RC)
set(gcf,'name','Original time series X and reconstruction RC')
clf;
plot(t,X,'b-',t,sum(RC(:,:),2),'r-');
legend('Original','Complete reconstruction');
title('ORIGINAL VS RECONSTRUCTED')

figure(7);  %COMPARE SHAM TO ARTIFACT FREE SIGNAL 
plot(SHAM_ORIGINAL_BEST)
hold on
plot(free_filtered)
hold off
title('COMPARE SHAM TO ARTIFACT FREE SIGNAL ')

R = corrcoef(SHAM_ORIGINAL_BEST(1, 3000:27000), free_filtered(1, 3000:27000));
similarity_best = R(1,2);

similarity_best_loop2 = 0;

for i= 0:1:25000
        
    low = i+1;
    
    high = 5000+i;
    R = corrcoef(free_filtered(1, low:high),  SHAM_ORIGINAL_BEST(1,low:high));
    similarity = R(1,2);
    
    if similarity > similarity_best_loop2
        similarity_best_loop2 = similarity;
        low_best_loop2 = low;
        high_best_loop2 = high;
    end         

end


%[Cxy,f] = mscohere(SHAM_filtered,free_filtered,[],[],[],500);

%figure(8);  
%plot(f, Cxy)



