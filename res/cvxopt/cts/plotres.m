el1sparse = csvread('el1sparse_det.csv'); 
l1ls = csvread('matlab_res.csv'); 
el11 = csvread('el11.csv');
loqo_el1sparse = csvread('res_loqo_el1sparse.csv'); 
loqo_el1 = csvread('res_loqo_el1.csv'); 

spar = [2*(1:50),125,150];
spar2 = [2,10,20,30,40,50,60,70,80,90,100,125,150];
spar3 = [2,50,100,125,150];
len1 = length(spar);
len2 = length(spar2);
len3 = length(spar3);

avg_time_sparse = zeros(len1,1); %param matrix 
avg_maxerror_sparse = zeros(len1,1);
avg_l1error_sparse = zeros(len1,1);
sd_time_sparse = zeros(len1,1);
sd_maxerror_sparse = zeros(len1,1);
sd_l1error_sparse = zeros(len1,1);
avg_time = zeros(len1,1); %param vector 
avg_maxerror = zeros(len1,1);
avg_l1error = zeros(len1,1);
sd_time = zeros(len1,1);
sd_maxerror = zeros(len1,1);
sd_l1error = zeros(len1,1);
avg_time_l1ls = zeros(len2,1); %l1ls 
avg_maxerror_l1ls = zeros(len2,1); 
avg_l1error_l1ls = zeros(len2,1);
sd_time_l1ls = zeros(len2,1);
sd_maxerror_l1ls = zeros(len2,1);
sd_l1error_l1ls = zeros(len2,1);
avg_time_loqosparse = zeros(len2,1); %loqo el1sparse
avg_l1error_loqosparse = zeros(len2,1);
avg_maxerror_loqosparse = zeros(len2,1);
sd_time_loqosparse = zeros(len2,1);
sd_l1error_loqosparse = zeros(len2,1);
sd_maxerror_loqosparse = zeros(len2,1);
avg_time_loqo = zeros(len3,1); %loqo el1
avg_l1error_loqo = zeros(len3,1);
avg_maxerror_loqo = zeros(len3,1);
sd_time_loqo = zeros(len3,1);
sd_l1error_loqo = zeros(len3,1);
sd_maxerror_loqo = zeros(len3,1);

counter_10 = 1;
for i=1:50
    avg_time_sparse(i) = mean(el1sparse(counter_10:(counter_10+9),2));
    avg_maxerror_sparse(i) = mean(el1sparse(counter_10:(counter_10+9),3));
    avg_l1error_sparse(i) = mean(el1sparse(counter_10:(counter_10+9),4));
    sd_time_sparse(i) = std(el1sparse(counter_10:(counter_10+9),2));
    sd_maxerror_sparse(i) = std(el1sparse(counter_10:(counter_10+9),3));
    sd_l1error_sparse(i) = std(el1sparse(counter_10:(counter_10+9),4));
       
    avg_time(i) = mean(el1(counter_10:(counter_10+9),2));
    avg_maxerror(i) = mean(el1(counter_10:(counter_10+9),3));
    avg_l1error(i) = mean(el1(counter_10:(counter_10+9),4));
    sd_time(i) = std(el1(counter_10:(counter_10+9),2));
    sd_maxerror(i) = std(el1(counter_10:(counter_10+9),3));
    sd_l1error(i) = std(el1(counter_10:(counter_10+9),4));

    counter_10 = counter_10+10;
end

counter_2 = counter_10;
for i=51:52
    avg_time_sparse(i) = mean(el1sparse(counter_2:(counter_2+1),2));
    avg_maxerror_sparse(i) = mean(el1sparse(counter_2:(counter_2+1),3));
    avg_l1error_sparse(i) = mean(el1sparse(counter_2:(counter_2+1),4));
    sd_time_sparse(i) = std(el1sparse(counter_2:(counter_2+1),2));
    sd_maxerror_sparse(i) = std(el1sparse(counter_2:(counter_2+1),3));
    sd_l1error_sparse(i) = std(el1sparse(counter_2:(counter_2+1),4));
    
    avg_time(i) = mean(el1(counter_2:(counter_2+1),2));
    avg_maxerror(i) = mean(el1(counter_2:(counter_2+1),3));
    avg_l1error(i) = mean(el1(counter_2:(counter_2+1),4));
    sd_time(i) = std(el1(counter_2:(counter_2+1),2));
    sd_maxerror(i) = std(el1(counter_2:(counter_2+1),3));
    sd_l1error(i) = std(el1(counter_2:(counter_2+1),4));
    
    counter_2 = counter_2+2;
end

counter_2 = 1;
for i=1:len2
    avg_time_loqosparse(i) = mean(loqo_el1sparse(counter_2:(counter_2+1),2));
    avg_maxerror_loqosparse(i) = mean(loqo_el1sparse(counter_2:(counter_2+1),3));
    avg_l1error_loqosparse(i) = mean(loqo_el1sparse(counter_2:(counter_2+1),4));
    sd_time_loqosparse(i) = std(loqo_el1sparse(counter_2:(counter_2+1),2));
    sd_maxerror_loqosparse(i) = std(loqo_el1sparse(counter_2:(counter_2+1),3));
    sd_l1error_loqosparse(i) = std(loqo_el1sparse(counter_2:(counter_2+1),4));
    
    avg_time_l1ls(i) = mean(l1ls(counter_2:(counter_2+1),2));
    avg_maxerror_l1ls(i) = mean(l1ls(counter_2:(counter_2+1),3));
    avg_l1error_l1ls(i) = mean(l1ls(counter_2:(counter_2+1),4));
    sd_time_l1ls(i) = std(l1ls(counter_2:(counter_2+1),2));
    sd_maxerror_l1ls(i) = std(l1ls(counter_2:(counter_2+1),3));
    sd_l1error_l1ls(i) = std(l1ls(counter_2:(counter_2+1),4));
    
    counter_2 = counter_2+2;
end

counter_2 = 1;
for i=1:len3
    avg_time_loqo(i) = mean(loqo_el1(counter_2:(counter_2+1),2));
    avg_maxerror_loqo(i) = mean(loqo_el1(counter_2:(counter_2+1),3));
    avg_l1error_loqo(i) = mean(loqo_el1(counter_2:(counter_2+1),4));
    sd_time_loqo(i) = std(loqo_el1(counter_2:(counter_2+1),2));
    sd_maxerror_loqo(i) = std(loqo_el1(counter_2:(counter_2+1),3));
    sd_l1error_loqo(i) = std(loqo_el1(counter_2:(counter_2+1),4));
   
    counter_2 = counter_2+2;
end
%%


fig = figure(1);
set(fig, 'Position', [1 1 900 600]);
semilogy(spar,avg_time_sparse,'bv',...
    spar,avg_time,'gv',...
    spar2,avg_time_l1ls,'rv',...
    spar2,avg_time_loqosparse,'kv',spar3,avg_time_loqo,'k^',...
    'LineWidth',2);
hold on;
semilogy(spar,avg_time_sparse,'b-',...
    spar,avg_time,'g-',...
    spar2,avg_time_l1ls,'r-',......
    spar2,avg_time_loqosparse,'k-',spar3,avg_time_loqo,'k-',...
    'LineWidth',2);
errorbar(spar,avg_time_sparse,sd_time_sparse,'b');
errorbar(spar,avg_time,sd_time,'g');
errorbar(spar2,avg_time_l1ls,sd_time_l1ls,'r');
errorbar(spar2,avg_time_loqosparse,sd_time_loqosparse,'k');
errorbar(spar3,avg_time_loqo,sd_time_loqo,'k');
title('Solution Time vs. Sparsity');
xlabel('Number of nonzeros in signal');
ylabel('Solution Time (seconds)');
legend('Simplex Matrix',...
     'Simplex Vector',...
    'L_1\_L_S',...
    'Loqo Matrix',...
    'Loqo Vector',...
		'Location','SouthEast');
hold off;
print(fig,'-dpng','plot_time.png');
close;

%%
fig = figure(1);
set(fig, 'Position', [1 1 900 600]);
plot(spar,avg_maxerror_sparse,'bv',...
    spar,avg_maxerror,'gv',...
    spar2,avg_maxerror_l1ls,'rv',...
    spar2,avg_maxerror_loqosparse,'kv',spar3,avg_maxerror_loqo,'k^',...
    'LineWidth',2);
hold on;
plot(spar,avg_maxerror_sparse,'b-',...
    spar,avg_maxerror,'g-',...
    spar2,avg_maxerror_l1ls,'r-',...
    spar2,avg_maxerror_loqosparse,'k-',spar3,avg_maxerror_loqo,'k-',...
    'LineWidth',2);
errorbar(spar,avg_maxerror_sparse,sd_maxerror_sparse,'b');
errorbar(spar,avg_maxerror,sd_maxerror,'g');
errorbar(spar2,avg_maxerror_l1ls,sd_maxerror_l1ls,'r');
errorbar(spar2,avg_maxerror_loqosparse,sd_maxerror_loqosparse,'k');
errorbar(spar3,avg_maxerror_loqo,sd_maxerror_loqo,'k');
title('Maximum Error vs. Sparsity');
xlabel('Number of nonzeros in signal');
ylabel('Maximum Error');
legend('Simplex Matrix',...
     'Simplex Vector',...
    'L_1\_L_S',...
    'Loqo Matrix',...
    'Loqo Vector',...
		'Location','NorthWest');
hold off;
print(fig,'-dpng','plot_maxerror.png');
close;

%%
fig = figure(1);
set(fig, 'Position', [1 1 900 600]);
plot(spar,avg_l1error_sparse,'bv',...
    spar,avg_l1error,'gv',...
    spar2,avg_l1error_l1ls,'rv',...
    spar2,avg_l1error_loqosparse,'kv',spar3,avg_l1error_loqo,'k^',...
    'LineWidth',2);
hold on;
plot(spar,avg_l1error_sparse,'b-',...
    spar,avg_l1error,'g-',...
    spar2,avg_l1error_l1ls,'r-',...
    spar2,avg_l1error_loqosparse,'k-',spar3,avg_l1error_loqo,'k-',...
    'LineWidth',2);
errorbar(spar,avg_l1error_sparse,sd_l1error_sparse,'b');
errorbar(spar,avg_l1error,sd_l1error,'g');
errorbar(spar2,avg_l1error_l1ls,sd_l1error_l1ls,'r');
errorbar(spar2,avg_l1error_loqosparse,sd_l1error_loqosparse,'k');
errorbar(spar3,avg_l1error_loqo,sd_l1error_loqo,'k');
title('L1 Error vs. Sparsity');
xlabel('Number of nonzeros in signal');
ylabel('L1 Error');
legend('Simplex Matrix',...
     'Simplex Vector',...
    'L_1\_L_S',...
    'Loqo Matrix',...
    'Loqo Vector',...
		'Location','NorthWest');
    hold off;
print(fig,'-dpng','plot_l1error.png');
close;