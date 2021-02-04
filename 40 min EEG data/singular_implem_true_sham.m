[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
[da,ca]= butter(3,[38/250 40/250], 'stop');


load('session1_SHAM.mat')      %%SHAM DEFINITION
sham_tacs_t1 = 320000;
sham_tacs_t2 = sham_tacs_t1 + 149999;
Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG_mine(1, sham_tacs_t1:sham_tacs_t2)));
Sham_SMA = filtfilt(da,ca,Sham_SMA);

%load('Test_Sham_ALPHA_PH.mat') 
%sham_tacs_t1 = 15032;   %%CHANGE sham boundary
%sham_tacs_t2 = sham_tacs_t1 + 29999;
%Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));


load('session2_DATA.mat')     %%RAW DATA DEFINITION
tacs_t1 = 90000;        
tacs_t2 = tacs_t1 + 149999; %5 MINUTES (5*60*500)

X_init = EEG_mine_ses2(4,tacs_t1:tacs_t2);  %FOR ORIGINAL TEST 

RAW_FILTERED = filtfilt(b,a,filtfilt(d,c,X_init));

Fs = 500;                       %FFT OF RAW FILTERED                          
L = 150000;            
Y = fft(RAW_FILTERED);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Y = fft(Sham_SMA);             %FFT OF SHAM
P22 = abs(Y/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);
f = Fs*(0:(150000/2))/150000;
figure(2);
plot(f,P1, 'k') 
hold on
plot(f,P12, 'b') 
title('filtered raw DATA and filtered SHAM')
xlabel('f(Hz)')
ylabel('|P1(f)|')
ylim([0 600000])


%%
%%STEP 1: set parameters and define input X

X = X_init - mean(X_init);
M = 200;    % window length = embedding dimension
N = length(X);
t = (1:N)';

figure(1);
plot(X)
%%
%%step 2: Calculate Covariance matrix C
tic;
Y=zeros(N-M+1,M);

for m=1:M
  Y(:,m) = X((1:N-M+1)+m-1);
end
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
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
end



RC_TO_QUIT = zeros(150000,1000);

inner_col = 1;
for i=1:1:length(RC(1,:))
    if max(abs(RC(:,i))) > 100000
        RC_TO_QUIT(:,inner_col) =  RC(:,i);
        inner_col = inner_col + 1;
    end
end

summed = sum(RC_TO_QUIT(:,1:100),2);

free_single = X - summed';

free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));
toc;
%%
%%STEP 6: Compare reconstructed to original


figure(6);
subplot(2,1,1)
plot(Sham_SMA, 'k');
subplot(2,1,2)
plot(free_filtered, 'r');

figure(7);
plot(free_filtered, 'r');
hold on
plot(Sham_SMA, 'k');


Y = fft(Sham_SMA);                %%FFT OF SHAM
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

L= 150000;
Y = fft(free_filtered);            %%FFT OF FREE FILTERED       
P22 = abs(Y/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);
figure(10);
plot(f,P1, 'k') 
hold on
plot(f,P12, 'r') 
title('sham (black) and filtered free (red)')
xlabel('f(Hz)')
ylabel('|P1(f)|')




 



