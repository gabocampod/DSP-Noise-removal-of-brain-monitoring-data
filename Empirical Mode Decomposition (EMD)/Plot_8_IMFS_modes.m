load('CEEMDAN_FAST_TEST2_5TO25_50_500_002_results')
load('ERP_sham_test1_test2_cc')

[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

ERP_SHAM = filtfilt(b,a,filtfilt(d,c,ERPSHAM.ES(1,:)));

modes = CEEMDAN_FAST_TEST2_5TO25_50_500_002_results;

t = 5:1/500:25/500;

figure(10);
for i = 1:1:11
    subplot(6,2,i)
    plot(modes(i,:)/10000)
end 


sum = zeros(1,max(size(modes)));

for i = 1:1:7
    sum(1,:) = sum(1,:) + modes(i,:);
end

sum_plot = sum(1,:) - modes(4,:) - modes(5,:);

sum_plot_filtered = filtfilt(b,a,filtfilt(d,c, sum_plot));

figure(9)
subplot(3,1,1)
plot(ERP_SHAM(1,24500:24999))
xlim([50 450])

subplot(3,1,2)
plot(sum_plot_filtered(1,:))
xlim([50 450])

subplot(3,1,3)
plot(sum_plot(1,:))
xlim([50 450])




