[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

%%
%%STEP 1: set parameters and define input X
%start = 18921 (0-5) 168920
%start 2 = 168920 (5-10) 183920
load('CoStim_Result_Subjectfdz1_Sess1_Stim2000mA.mat')
sham_tacs_t1 = 320000;
sham_tacs_t2 = sham_tacs_t1 + 149999;
Sham_SMA=filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));


%top lim = 1218870
load('CoStim_Result_Subjectfdz2_Sess2_Stim2000mA.mat')
tacs_t1 = 168920;        
tacs_t2 = tacs_t1 + 149999; %5 MINUTES (5*60*500)

X_init = CoStim_Result.EEG.Data(1,tacs_t1:tacs_t2);  %FOR ORIGINAL TEST 
X = X_init - mean(X_init);
M = 200;    % window length = embedding dimension
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
for m=1:20
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
end;

%%
%%STEP 6: Compare reconstructed to original
summed = sum(RC(:,1:2),2);
free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));
toc;

figure(6);
plot(Sham_SMA, 'k');
hold on
plot(free_filtered, 'r');




 



