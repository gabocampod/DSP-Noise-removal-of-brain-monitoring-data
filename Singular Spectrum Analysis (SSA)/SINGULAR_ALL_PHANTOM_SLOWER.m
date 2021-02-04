[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5001;
upper = 25000;
%%
%%STEP 1: set parameters and define input X

load('Test_Sham_ALPHA_PH.mat') 
sham_tacs_t1 = 15048;   %%CHANGE sham boundary !!!!!!!!!!
sham_tacs_t2 = sham_tacs_t1 + 29999;
ERP_Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));

load('Test_1mA_5Hz_Alpha_PH.mat') %%CHANGE data !!!!!!!!!!
tacs_t1 = 15031;         %%CHANGE data boundary !!!!!!!!!!
tacs_t2 = tacs_t1 + 29999;

X_init = EEG.Data(1,tacs_t1:tacs_t2);  %FOR ORIGINAL TEST 
X = X_init - mean(X_init);
M = 1000;    % window length = embedding dimension
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
for m=1:10
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m= 1:10
  subplot(5,2,(m));
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  %ylim([-1 1]);
end;

%%
%%STEP 6: Compare reconstructed to original
summed = sum(RC(:,1:2),2);
free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));

figure(7);
plot(ERP_Sham_SMA)
hold on
plot(free_filtered)
hold off

R = corrcoef(ERP_Sham_SMA(1,lower:upper), free_filtered(1,lower:upper));

similarity = R(1,2);


noise_SUBJECT2 = (ERP_Sham_SMA(1,lower:upper) - free_filtered(1,lower:upper));
SNR_SUBJECT2= (rms(free_filtered(1,lower:upper) / rms(noise_SUBJECT2) ))^2;
SNR_SUBJECT2_db = db(SNR_SUBJECT2, 'power');



