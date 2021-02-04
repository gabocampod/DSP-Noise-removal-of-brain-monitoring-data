[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
[da,ca]= butter(3,[38/250 40/250], 'stop');


load('session1_SHAM.mat')      %%SHAM DEFINITION
sham_tacs_t1 = 320000;
sham_tacs_t2 = sham_tacs_t1 + 149999;
Sham_SMA=filtfilt(b,a,filtfilt(d,c,EEG_mine(1, sham_tacs_t1:sham_tacs_t2)));
%Sham_SMA = filtfilt(da,ca,Sham_SMA);

load('session2_DATA.mat')     %%RAW DATA DEFINITION
tacs_t1 = 90000;        
tacs_t2 = tacs_t1 + 149999; %5 MINUTES (5*60*500)

X_init = EEG_mine_ses2(4,tacs_t1:tacs_t2);  %FOR ORIGINAL TEST 
X = X_init - mean(X_init);
M = 1000;    
N = length(X);
t = (1:N)';

RAW_FILTERED = filtfilt(b,a,filtfilt(d,c,X_init));

load('RC_1000_FROM90K.MAT') 

figure(29)
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
end

figure(120);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 21;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end

figure(140);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 41;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end

figure(160);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 61;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end

figure(180);
set(gcf,'name','Reconstructed components RCs')
clf;
k = 81;
for m=1:20
  subplot(10,2,m);
  plot(t,RC(:,k),'r-');
  ylabel(sprintf('RC %d',k));
  k = k+1;
end
%}

RC_TO_QUIT = zeros(150000,1000);
val = zeros(1000);
R_best = 0;
tresh_best = 0;

%FFT SHAM
f_samp_b = 500;
eeg_b = Sham_SMA;

% FFT settings

epoch_duration = 1; n_fft = 2^15; windows = f_samp_b * epoch_duration; overlap = 0.1 * f_samp_b;
eeg_b_to_analyse = eeg_b(1, 40*f_samp_b:end);

[pxx_b_sham, f_b_sham] = pwelch(eeg_b_to_analyse,windows,overlap,n_fft,f_samp_b);

for tresh = 50000:-1000:1000
    
    inner_col = 1;
    for i=1:1:length(RC(1,:))
        if max(abs(RC(:,i))) > tresh
            RC_TO_QUIT(:,inner_col) =  RC(:,i);
            inner_col = inner_col + 1;
            val(inner_col) = i;
        end
    end
    
   summed = sum(RC_TO_QUIT(:,:),2);
   free_single = X - summed';
   
   free_filtered = filtfilt(b,a,filtfilt(d,c,free_single));
   
    %FFT FILTERED RESULTS

    f_samp_b = 500;
    eeg_b = free_filtered;

    % FFT settings

    epoch_duration = 1; n_fft = 2^15; windows = f_samp_b * epoch_duration; overlap = 0.1 * f_samp_b;

    eeg_b_to_analyse = eeg_b(1, 40*f_samp_b:end);

    [pxx_b, f_b] = pwelch(eeg_b_to_analyse,windows,overlap,n_fft,f_samp_b);

    noise_1 = (pxx_b_sham - pxx_b);
    SNR_SUBJECT1= (rms(pxx_b) / rms(noise_1))^2;
    R = db(SNR_SUBJECT1, 'power');
   
   if (R > R_best)
       R_best = R;
       tresh_best = tresh;
       free_filtered_best= free_filtered;
   end
   
end

%FFT FILTERED RESULTS

f_samp_b = 500;
eeg_b = free_filtered_best;

% FFT settings

epoch_duration = 1; n_fft = 2^15; windows = f_samp_b * epoch_duration; overlap = 0.1 * f_samp_b;
eeg_b_to_analyse = eeg_b(1, 40*f_samp_b:end);

[pxx_b, f_b] = pwelch(eeg_b_to_analyse,windows,overlap,n_fft,f_samp_b);


figure(6);
subplot(2,1,1)
plot(Sham_SMA, 'k');
subplot(2,1,2)
plot(free_filtered_best, 'r');
title('SHAM (black) vs FREE OF ARTIFACT (red)')

figure(7);
plot(free_filtered_best, 'r');
hold on
plot(Sham_SMA, 'k');
title('SHAM (black) and FREE OF ARTIFACT (red)')

%%FFT RAW FILTERED
f_samp_b = 500;
eeg_b = X_init;

%FFT settings

epoch_duration = 1; n_fft = 2^15; windows = f_samp_b * epoch_duration; overlap = 0.1 * f_samp_b;

eeg_b_to_analyse = eeg_b(1, 40*f_samp_b:end);

[pxx_b_raw, f_b_raw] = pwelch(eeg_b_to_analyse,windows,overlap,n_fft,f_samp_b);

figure(12);
semilogx(f_b_sham,10*log10(pxx_b_sham),'LineWidth',2,'Color','b'); hold all
semilogx(f_b,10*log10(pxx_b),'LineWidth',2,'Color','r'); hold all
semilogx(f_b_raw,10*log10(pxx_b_raw),'LineWidth',2,'Color','k'); hold all
title('FREE OF ARTIFACT FILTERED (red) vs sham (red) vs raw (black)')


load('Test_Sham_ALPHA_PH.mat') 
sham_tacs_t1 = 15030;   %%CHANGE sham boundary
sham_tacs_t2 = sham_tacs_t1 + 29999;
Sham_SMA_ALPHA=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));


load('Test_Sham_ERP_PH.mat') 
sham_tacs_t1 = 15032;   %%CHANGE sham boundary
sham_tacs_t2 = sham_tacs_t1 + 29999;
Sham_SMA_ERP=filtfilt(b,a,filtfilt(d,c,EEG.Data(1, sham_tacs_t1:sham_tacs_t2)));





