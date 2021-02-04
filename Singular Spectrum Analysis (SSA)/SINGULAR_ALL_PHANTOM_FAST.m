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

load('Test_1mA_5Hz_ALPHA_PH.mat') %%CHANGE data !!!!!!!!!!
tacs_t1 = 15031;        %%CHANGE data boundary !!!!!!!!!!
tacs_t2 = tacs_t1 + 29999;

X_init = EEG.Data(1,tacs_t1:tacs_t2);  %FOR ORIGINAL TEST 
X = X_init - mean(X_init);
M = 1000;    % window length = embedding dimension
N = length(X);
t = (1:N)';

%%
%%step 2: Calculate Covariance matrix C
tic;
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
for m=1:4
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
  plot(t,RC(:,m)/10000,'r-');
  ylabel(sprintf('RC %d',m));
end;
xlabel('samples');  

%%
%%STEP 6: Compare reconstructed to original
summed = sum(RC(:,1:2),2);
free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));
toc;

Data_to_save = free_filtered;

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




