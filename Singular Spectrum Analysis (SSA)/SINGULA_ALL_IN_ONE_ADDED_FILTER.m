[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

lower = 5001;
upper = 25000;
%%
%%STEP 1: set parameters and define input X

load('Test_Sham_ALPHA_PH.mat') 
sham_tacs_t1 = 15048;   %%CHANGE sham boundary !!!!!!!!!!
sham_tacs_t2 = sham_tacs_t1 + 29999;
ERP_Sham_SMA = filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));


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
for m=1:2
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;

%%
%%STEP 6: Compare reconstructed to original
summed = sum(RC(:,1:2),2);
free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));
%free_filtered = filtfilt(da,ca, free_filtered);
toc;

SINGULAR_RESULT_SUBJECT2 = free_filtered;

%%
%%COEFFICIENT
R1 = corrcoef(ERP_Sham_SMA(1,lower:upper), SINGULAR_RESULT_SUBJECT2(1,lower:upper));
Coefficient_SUBJECT2 = R1(1,2);

%%
%%RMS
RMS_SUBJECT2 = 10^-3.*rms(ERP_Sham_SMA(1,lower:upper) - SINGULAR_RESULT_SUBJECT2(1,lower:upper));

%%
%%SNR 
noise_SUBJECT2 = (ERP_Sham_SMA(1,lower:upper) - SINGULAR_RESULT_SUBJECT2(1,lower:upper));
SNR_SUBJECT2= (rms(SINGULAR_RESULT_SUBJECT2(1,lower:upper) / rms(noise_SUBJECT2) ))^2;
SNR_SUBJECT2_db = db(SNR_SUBJECT2, 'power');
%%
%%Coherence
[cxy_SUBJECT1,f] = mscohere(ERP_Sham_SMA(1,lower:upper), SINGULAR_RESULT_SUBJECT2(1,lower:upper), [],[],[],500);

figure(1)
plot(f, cxy_SUBJECT1);

%%
%%COHERENCE VALUE

Fs = 500;            % Sampling frequency                         
L = 20000;             % Length of signal
FFT_SHAM_FILTERED = fft(ERP_Sham_SMA(1,lower:upper));
FFT_FREE_FILTERED = fft(SINGULAR_RESULT_SUBJECT2(1,lower:upper));

P2 = abs(FFT_SHAM_FILTERED/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P2_2 = abs(FFT_FREE_FILTERED/L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);

R_FFT = corrcoef(P2, P2_2);
Coherence_value = R_FFT(1,2);

%%
%%RUN IN ARTIFACT
low_in = 1;
high_in = 500;
similarity_in = 0;
treshhold = 5; %db
while similarity_in < treshhold
    
    noise_in = (ERP_Sham_SMA(1,low_in:high_in) - SINGULAR_RESULT_SUBJECT2(1,low_in:high_in));
    SNR_in= (rms(SINGULAR_RESULT_SUBJECT2(1,low_in:high_in) / rms(noise_in) ))^2;
    SNR_in_db = db(SNR_in, 'power');
    similarity_in = SNR_in_db;
    
    low_in = low_in +1;
    high_in = high_in +1;
end

RUN_IN_ARTIFACT = (low_in - 1)/500;


%%
%%RUN OUT ARTIFACT
low_out = 29501;
high_out = 30000;
similarity_out = 0;

while similarity_out < treshhold
    
    noise_out = (ERP_Sham_SMA(1,low_out:high_out) - SINGULAR_RESULT_SUBJECT2(1,low_out:high_out));
    SNR_out= (rms(SINGULAR_RESULT_SUBJECT2(1,low_out:high_out) / rms(noise_out) ))^2;
    SNR_out_db = db(SNR_out, 'power');
    similarity_out = SNR_out_db;
    
    low_out = low_out - 1;
    high_out = high_out - 1;
end
high_out_tt = high_out + 1;
high_out = 30000 - high_out + 1;

RUN_OUT_ARTIFACT = (high_out - 1)/500;



figure(2);
plot(ERP_Sham_SMA(1,:))
hold on
plot(SINGULAR_RESULT_SUBJECT2(1,:))





