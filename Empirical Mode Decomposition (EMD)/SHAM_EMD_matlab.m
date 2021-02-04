load('ERP_sham_test1_test2_cc')
load('CEEMDAN_FAST_20000_50_500_002_RESULTS')
load('CEEMDAN_FAST_20000points_50_500_02_RESULTS')
load('CEEMDAN_of_mode_4_results')
load('CEEMDAN_FAST_TEST2_50_500_002_results.mat')

[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

X = ERPSHAM.ES(1,5000:24999);
Y = ERPTEST1.E1(1,5000:24999);
%% ERP DATA MODES PLOT

modes = CEEMDAN_FAST_TEST2_50_500_002_results;

figure(1);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:))
end 

reconstructed_mode4 = zeros(1,20000);

for i= 1:1:11
    reconstructed_mode4(1,:) = reconstructed_mode4(1,:) + modes_of_mode(i,:);
end

RR = corrcoef(reconstructed_mode4(1,:), modes(4,:));

similarityR = RR(1,2);

figure(21);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes_of_mode(i,:))
end 


%modes(4,:)= modes(4,:) - modes_of_mode(4,:);

reconstructed= zeros(1,20000);

for i= 1:1:11
    reconstructed(1,:) = reconstructed(1,:) + modes(i,:);
end


%% SHAM MODES PLOT

%modes_sham = CEEMDAN_FAST__SHAM_20000_50_500_002; reconstructed_sham =

%% Create reconstructed_free and filters
reconstructed_artif_free = reconstructed(1,:) - modes(4,:);

SHAM_after_filter = filtfilt(b,a,filtfilt(d,c,X));
FREE_after_filter = filtfilt(b,a,filtfilt(d,c,reconstructed_artif_free));

SHAM_only_high_pass = filtfilt(b,a,X);
FREE_only_high_pass = filtfilt(b,a,reconstructed_artif_free);

%% Plots for filters

figure(3);
plot(SHAM_after_filter)
hold on
plot(FREE_after_filter)
hold off
title('AFTER FILTER')

similarity_best = 0;

for i= 0:1:15000
    low = 1 + i;
    high = 5000 + i;
    
    R = corrcoef(SHAM_after_filter(1, low: high), FREE_after_filter(1, low: high));
    similarity = R(1,2);
    
    if similarity > similarity_best
        similarity_best = similarity;
        low_best = low;
        high_best = high;
    end
end
 
R_check = corrcoef(SHAM_after_filter(1, 6600: 11600), FREE_after_filter(1, 6600: 11600));
similarity_check = R_check(1,2);

figure(3);
plot(SHAM_after_filter(1, low_best: high_best))
hold on
plot(FREE_after_filter(1, low_best: high_best))
hold off
xlim([0 5000])
title('AFTER FILTER')

%figure(4);
%plot(X)
%hold on
%plot(reconstructed_artif_free)
%hold off
%title('BEFORE FILTER')

%figure(5);
%plot(SHAM_only_high_pass)
%hold on
%plot(FREE_only_high_pass)
%hold off
%title('Only high pass filter')
