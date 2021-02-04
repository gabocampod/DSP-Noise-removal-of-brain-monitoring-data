load('CEEMDAN_FAST_20000_50_500_002_RESULTS')
load('ERP_sham_test1_test2_cc')
[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
%%
%%STEP 1: set parameters and define input X

modes = CEEMDAN_FAST_20000_50_500_002;

 figure(10);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:))
end

X_init = modes(4,:);  %FOR RECONSTRUCTED EMD
%X = filtfilt(b,a,X_init);
X = X_init - mean(X_init);
%X = X_init;  
M = 1500;    % window length = embedding dimension
N = length(X);
t = (1:N)';

figure(1);
subplot(2,1,1)
plot(X_init(1,:))
subplot(2,1,2)
plot(X(1,:))
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

figure(3);
set(gcf,'name','Eigenvectors RHO and eigenvalues LAMBDA')
clf;
subplot(3,1,1);
plot(LAMBDA,'o-');
subplot(3,1,2);
plot(RHO(:,1:2), '-');
legend('1', '2');
subplot(3,1,3);
plot(RHO(:,3:4), '-');
legend('3', '4');

%%
%%STEP 4: Calculate prncipal components (PC) PC = Y*RHO;
PC = Y*RHO;

for m=1:4
  subplot(4,1,m);
  plot(t(1:N-M+1),PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  %ylim([-10 10]);
end;

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
for m=1:4
  subplot(4,1,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  %ylim([-1 1]);
end;

%%
%%STEP 6: Compare reconstructed to original

SHAM = ERPSHAM.ES(1,5000:24999);
SHAM_filtered = filtfilt(b,a,filtfilt(d,c,SHAM));

summed = sum(RC(:,1:2),2);
modes(4,:) = X - summed';

reconstructed = zeros(1,20000);

for i= 1:1:11
    reconstructed(1,:) = reconstructed(1,:) + modes(i,:);
end

free_reconstructed_single = filtfilt(b,a,filtfilt(d,c,reconstructed));

figure(6);
set(gcf,'name','Original time series X and reconstruction RC')
clf;
plot(t,X,'b-',t,sum(RC(:,:),2),'r-');
legend('Original','Complete reconstruction');

figure(7);
plot(SHAM_filtered)
hold on
plot(free_reconstructed_single)
hold off

R = corrcoef(SHAM_filtered, free_reconstructed_single);

similarity = R(1,2);

%%No longer spikes at start and end, meaning some components of mode 4 are making so that 
%%signal levels to correct value






