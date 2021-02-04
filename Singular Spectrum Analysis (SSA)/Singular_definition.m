load('ERP_sham_test1_test2_cc')

[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
%%
%%STEP 1: set parameters and define input X

X_init = ERPTEST1.E1(1,5000:24999);  %FOR ORIGINAL TEST 
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

figure(2);
set(gcf,'name','Covariance matrix');
clf;
imagesc(Cemb);
axis square
set(gca,'clim',[-1 1]);
colorbar

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
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:10
  subplot(5,2,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
end;

%%
%%STEP 6: Compare reconstructed to original

SHAM = ERPSHAM.ES(1,5000:24999);
SHAM_filtered = filtfilt(b,a,filtfilt(d,c,SHAM));

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
plot(SHAM_filtered)
hold on
plot(free_filtered)
hold off
title('COMPARE SHAM TO ARTIFACT FREE SIGNAL ')

R = corrcoef(SHAM_filtered, free_filtered);

similarity = R(1,2);

[Cxy,f] = mscohere(SHAM_filtered,free_filtered,[],[],[],500);

%figure(8);  
%plot(f, Cxy)



