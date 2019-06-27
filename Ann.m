close all;
clear;

%% 训练参数调试
global n_max Group_Span Intra_Group_Power_Step Inter_Group_Power_Step Hide_Layer_Cells x_bound_factor_P2 y_bound_factor_P2 x_bound_factor_P1
Group_Span = 3;
n_max = 20;        % learning rate
Intra_Group_Power_Step = 20;        % mW
Inter_Group_Power_Step = 2;        % mW
Hide_Layer_Cells = 4;
x_bound_factor_P2 = 5.5;
x_bound_factor_P1 = 4.4;
y_bound_factor_P2 = 1;

%% 双环仿真
if exist('PM1.mat')
    load('PM1.mat');
    monitor = PM1;
else
t1=0.801; %耦合强度(传输系数
t2=0.9764;
t3=0.9899;
k1=sqrt(1-t1^2); %由t推出的交叉耦合系数
k2=sqrt(1-t2^2);
k3 = sqrt(1-t3^2);

R=10e-6; %环半径
phase_mis1 = 0.0*pi; 

% w_min=-2e-9; w_step=0.001e-9; w_max=8e-9;
% misalignment = w_min:w_step:w_max;

% phase_shift = (-0.5:0.001:0.5)*pi;
phase_range = (-1.5:0.001:1.5)*pi; %一个周期，精度

alpha = 1.98357e2;          % transmission loss factor, equal to 17.23 dB/cm


a = exp(-alpha*(2*pi*R));


f1 = t1*t2*t3*a; f2 = t1*k2*t3*a;
f3 = -f2^2/(1-f1)+f1;
r_max = (f1-f3)/(1-f3);
r_min = (f1-f3)/(1+f3);
ph1 = exp(1i*(phase_range));
ph2 = exp(1i*(phase_range + phase_mis1));
M1=zeros(length(ph1), length(ph2));
M2=M1; ph1s = M1; ph2s = M1;
for index1 = 1:length(ph1)
    for index2 = 1:length(ph2)
        ph1s(index1, index2) = phase_range(index1);
        ph2s(index1, index2) = phase_range(index2) + phase_mis1;
        M1(index1,index2) = (k1*k3)^2*a^0.5* abs( (ph2(index2)-f1)./((ph1(index1)-f1).*(ph2(index2)-f1)+f2^2) ).^2;  %Monitor1的公式
        M2(index1,index2) = (k1*k3)^2*a^0.5* abs( (f2)./((ph1(index1)-f1).*(ph2(index2)-f1)+f2^2) ).^2;              %Monitor2的公式
    end
end

P1 = [0 0.5 1 1 1 1];         % mW
P2 = [1 1 1 0.5 0.2 0];         % mW

for k = 6
PM1 = M1*P1(k)+M2*P2(k);
% PM1 = M2*P2(k);
PM1 = PM1./max(max(PM1));
PM1_dB = 10*log10(PM1);
min(min(PM1))
% figure;
% pcolor(ph1s./pi,ph2s./pi,PM1);
% shading interp;
% xlabel('round trip phase1 (pi)');
% ylabel('round trip phase2 (pi)');
% title(['P1 = ' num2str(P1(k)) ' mW; ' 'P2 = ' num2str(P2(k)) ' mW'])
% set(gca, 'CLim', [0.01 1]);
% colorbar;
% colormap(jet);
end
    save('PM1.mat', 'PM1', 'phase_range');
monitor = PM1; %接口函数
end

%% 双环采样，构造ADC采样集
global ADC_Set DAC_Steps

DAC_Steps = 256;
% Max_power = 128;        % mW, 0.5 mW/step
% Intra_Group_Power_Step = 20;        % mW % ? 
% %Inter_Group_Power_Step = 2;        % mW
% Inter_Group_Power_Step = 2;   
% %Group_Data_Count = 3;
% Group_Data_Count = 20;
%PM1=PM1';
%[through, drop, original_monitor]=Ring_Script(3); % 是3或者5还不确定，都试一下
%[through, drop, original_monitor]=Ring_Script(3);

%monitor = monitor / max(monitor) + (rand(size(monitor)) - 0.5)*0e-3/3.3;     % 添加噪声，about 3-mVpp noise
% monitor(monitor < 0) = 0; 
% monitor(monitor > 1) = 1;
monitor = round(monitor * 2^16);     % 11-bit ADC
DAC_Set_1D = round(sqrt((0:DAC_Steps-1) / DAC_Steps) * 2 ^ 12).';     % 12-bit DAC
%DAC_Set_2D_1 = 
index = round((DAC_Set_1D/2^12).^2 * (length(monitor)-1)) + 1;
ADC_Set = 10*log10(monitor(index,index));
%ADC_Set = monitor(index,index);
% for k = 6

figure;
% ph1s = ph1s(index,index);
% ph2s = ph2s(index,index);
%pcolor(ph1s./pi,ph2s./pi,PM1);
ADC_Set = ADC_Set./max(max(ADC_Set));
% pcolor(ph1s./pi,ph2s./pi,ADC_Set);
pcolor(phase_range(index) / pi,phase_range(index) / pi,ADC_Set.');
shading interp;
xlabel('round trip phase1 (pi)');
ylabel('round trip phase2 (pi)');
% title(['P1 = ' num2str(P1(k)) ' mW; ' 'P2 = ' num2str(P2(k)) ' mW'])
%set(gca, 'CLim', [0.01 1]);
colorbar;
colormap(jet);
% end

%% 鞍点搜索—寻找最大值的坐位置

index_P1=randi(256);
index_P2=index_P1 + randi(round([-256 256]*0.2));
if index_P2 < 1
    index_P2 = 1;
end
if index_P2 > 256
    index_P2 = 256;
end
fprintf('initial_index_P1 = %u, initial_index_P2 = %u\n', index_P1, index_P2);

P_max = 256;
delt_p1=32;
delt_p2=8;
delt_p3=2;
measure_value1 = Get_ADC(index_P1,index_P2);
PM_max=measure_value1;
fprintf('initial_ADC = %f\n', PM_max);
index_searchP1=index_P1;
index_searchP2=index_P2;
deltP1P2 = 6;
while deltP1P2>=5
%sub process i
for  k = 0:(P_max/delt_p1)
   index_P1=index_P1+delt_p1*k;
   index_P1= mod(index_P1, P_max)+1;
   index_P2=index_P2 +delt_p1*k ;
   index_P2= mod(index_P2, P_max)+1;
   measure_value2 = Get_ADC(index_P1,index_P2);
   if measure_value2 > PM_max
       PM_max=measure_value2;
       index_searchP1=index_P1;
       index_searchP2=index_P2;
   end
end

index_P1=index_searchP1;
index_P2=index_searchP2;
fprintf('sub1_index_P1 = %u, sub1_index_P2 = %u\n', index_P1, index_P2);
fprintf('sub1_ADC = %f\n', PM_max);

%sub process ii
P_start=index_P2 - P_max/8;
P_end=index_P2 + P_max/8;
measure_value1 = Get_ADC(index_P1,P_start);
PM_min=measure_value1;
%fprintf('ADC = %f\n', PM_min);
index_searchP2 = P_start;
index_P2=P_start;
% for  k = 0:(P_max/delt_p2)
while index_P2 < P_end
   index_P2=index_P2 +delt_p2;
   %index_P2= mod(index_P2, P_max)+1;
   measure_value2 = Get_ADC(index_P1,index_P2);
   if measure_value2 < PM_min
       PM_min=measure_value2;
       index_searchP2=index_P2;
   end
end

index_P2=index_searchP2;
fprintf('sub2_index_P1 = %u, sub2_index_P2 = %u\n', index_P1, index_P2);
fprintf('sub2_ADC = %f\n', PM_min);
fprintf('%%%%%%\n');
%sub process iii
P_start=index_P1 - P_max/32;
P_end=index_P1 + P_max/32;
measure_value1 = Get_ADC(P_start,index_P2);
PM_max=measure_value1;
% fprintf('ADC = %f\n', PM_max);
index_searchP1 = P_start;
index_P1=P_start;
% for  k = 0:(P_max/delt_p2)
while index_P1 < P_end
   index_P1=index_P1 +delt_p3;
   %index_P2= mod(index_P2, P_max)+1;
   measure_value2 = Get_ADC(index_P1,index_P2);
   if measure_value2 > PM_max
       PM_max=measure_value2;
       index_searchP1=index_P1;
   end
end

index_P1=index_searchP1;
fprintf('sub3_index_P1 = %u,sub3_ index_P2 = %u\n', index_P1, index_P2);
fprintf('sub3_ADC = %f\n', PM_max);

%sub process iv
P_start=index_P2 - P_max/64;
P_end=index_P2 + P_max/64;
measure_value1 = Get_ADC(index_P1,P_start);
PM_min=measure_value1;
%fprintf('ADC = %f\n', PM_min);
index_searchP2 = P_start;
index_P2=P_start;
% for  k = 0:(P_max/delt_p2)
while index_P2 < P_end
   index_P2=index_P2 +delt_p3;
   %index_P2= mod(index_P2, P_max)+1;
   measure_value2 = Get_ADC(index_P1,index_P2);
   if measure_value2 < PM_min
       PM_min=measure_value2;
       index_searchP2=index_P2;
   end
end

index_P2=index_searchP2;
deltP1P2 = abs(index_P1- index_P2);
end
fprintf('sub4_index_P1 = %u, sub4_index_P2 = %u\n', index_P1, index_P2);
fprintf('sub4_ADC = %f\n', PM_min);
target_set_index_P1= index_P1;
target_set_index_P2= index_P2;
%% 构造训练集
 DAC_Steps = 256;
 Max_power = 128;        % mW, 0.5 mW/step
intra_step = Intra_Group_Power_Step/Max_power*DAC_Steps; %40
inter_step = Inter_Group_Power_Step/Max_power*DAC_Steps; % 4
Group_Count_1D = 21;
CenterPoint_Cell_1 = cell(Group_Count_1D,Group_Count_1D);
FivePoint_Cell_2 = cell(1,5);
IndexAndValue_Cell_3 = cell(1,3);
%先把中心点矩阵找到
start_index_P1 = index_P1 - inter_step * (Group_Count_1D-1)/2;
start_index_P2 = index_P2 - inter_step * (Group_Count_1D-1)/2;
for i = 1:Group_Count_1D
    for j = 1:Group_Count_1D
        
       IndexAndValue_Cell_3={start_index_P1-intra_step,start_index_P2-intra_step,...
           Get_ADC(start_index_P1-intra_step,start_index_P2-intra_step)};
       FivePoint_Cell_2(1,1)={IndexAndValue_Cell_3};
       IndexAndValue_Cell_3={start_index_P1-intra_step,start_index_P2+intra_step,...
           Get_ADC(start_index_P1-intra_step,start_index_P2+intra_step)};
       FivePoint_Cell_2(1,2)={IndexAndValue_Cell_3}; 
       IndexAndValue_Cell_3={start_index_P1,start_index_P2,...
           Get_ADC(start_index_P1,start_index_P2)};
       FivePoint_Cell_2(1,3)={IndexAndValue_Cell_3}; 
       IndexAndValue_Cell_3={start_index_P1+intra_step,start_index_P2-intra_step,...
           Get_ADC(start_index_P1+intra_step,start_index_P2-intra_step)};
       FivePoint_Cell_2(1,4)={IndexAndValue_Cell_3}; 
       IndexAndValue_Cell_3={start_index_P1+intra_step,start_index_P2+intra_step,...
           Get_ADC(start_index_P1+intra_step,start_index_P2+intra_step)};
       FivePoint_Cell_2(1,5)={IndexAndValue_Cell_3}; 

       CenterPoint_Cell_1(i,j) = {FivePoint_Cell_2};
       j=j+1;
       start_index_P1 = start_index_P1 + inter_step;
    end
    start_index_P1 = index_P1 - inter_step * (Group_Count_1D-1)/2;
    start_index_P2 = start_index_P2 + inter_step;     
    i=i+1;
    
end
%% 差模方法，对P1来说input_data441x5不改变，对P2来说input_data要改变为400X3的数据，去掉第一行和最后一行
%% 差模方法（更新后）：对P2来说input_data要改变为420X5的数据，即每个5点Group和它垂直方向上相邻的5点Group做差值
ADC_Groups_P2 = zeros(Group_Count_1D*(Group_Count_1D-Group_Span), 5); %420组5点数据
k=1;
for i = 1:Group_Count_1D-1
    for j = 1:Group_Count_1D
        for m = 1:5
            ADC_Groups_P2(k,m) = CenterPoint_Cell_1{i,j}{1,m}{1,3} - CenterPoint_Cell_1{i+1,j}{1,m}{1,3};
        end
%         ADC_Groups_P2(k,1) = CenterPoint_Cell_1{i,j+1}{1,4}{1,3}- CenterPoint_Cell_1{i,j}{1,1}{1,3};
%         ADC_Groups_P2(k,2) = (CenterPoint_Cell_1{i,j+1}{1,3}{1,3}- (CenterPoint_Cell_1{i,j}{1,1}{1,3}...
%             +CenterPoint_Cell_1{i,j}{1,2}{1,3})*0.5)*0.5;
%         ADC_Groups_P2(k,3) = CenterPoint_Cell_1{i,j+1}{1,5}{1,3}- CenterPoint_Cell_1{i,j}{1,2}{1,3};
        k = k+1;
        j = j+1;
    end
    i = i+1;
end

x_bound2 = [min(ADC_Groups_P2(:,3)), max(ADC_Groups_P2(:,1))];
input_data_P2 = (ADC_Groups_P2 - (x_bound2(1) + x_bound2(2))/2) / (x_bound2(2) - x_bound2(1))*x_bound_factor_P2;

ADC_Groups_P1 = zeros(Group_Count_1D^2,5); %400组5点数据
k=1;
for i = 1:Group_Count_1D
    for j = 1:Group_Count_1D
        for m = 1:5
            ADC_Groups_P1(k,m) = CenterPoint_Cell_1{i,j}{1,m}{1,3};
        end
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
x_bound1 = [min(ADC_Groups_P1(:,5)), max(ADC_Groups_P1(:,1))];
input_data_P1 = (ADC_Groups_P1 - (x_bound1(1) + x_bound1(2))/2) / (x_bound1(2) - x_bound1(1)) *x_bound_factor_P1;

k=1;
TO_Groups_P1 = zeros(Group_Count_1D^2,1);
for i = 1:Group_Count_1D
    for j = 1:Group_Count_1D
        TO_Groups_P1(k,1) = CenterPoint_Cell_1{i,j}{1,3}{1,1} - index_P1;
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
y_bound_P1 = [min(TO_Groups_P1), max(TO_Groups_P1)];
output_data_P1 = (TO_Groups_P1 - y_bound_P1(1)) / (y_bound_P1(2) - y_bound_P1(1));

%outputP2
k=1;
TO_Groups_P2 = zeros(Group_Count_1D*(Group_Count_1D-Group_Span),1);
for i = 1:Group_Count_1D-1
    for j = 1:Group_Count_1D
%         TO_Groups_P2(k,1) = CenterPoint_Cell_1{i,j+1}{1,3}{1,2} - index_P2;
        TO_Groups_P2(k,1) = CenterPoint_Cell_1{i,j}{1,3}{1,2} - index_P2;
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
y_bound_P2 = [min(TO_Groups_P2), max(TO_Groups_P2)] / y_bound_factor_P2;
output_data_P2 = (TO_Groups_P2 - y_bound_P2(1)) / (y_bound_P2(2) - y_bound_P2(1));

figure; hold on;
for k = 1:5
h = plot(ADC_Groups_P1(:,k), TO_Groups_P1,'o', 'markersize', 6);
set(h, 'MarkerFaceColor', h.Color);
xlabel('ADC Groups');
ylabel('TO Groups P1');
end
%legend({'TO Groups P1'});
%legend({'1', '2', '3'});
hold off;
figure; hold on;
for k = 1:5
h = plot(input_data_P1(:,k), output_data_P1,'o', 'markersize', 6);
set(h, 'MarkerFaceColor', h.Color);
xlabel('input_data_P1');
ylabel('output_data_P1');
end
%legend({'TO Groups P1'});
%legend({'1', '2', '3'});
hold off;
figure; hold on;
for k = 1:5
h = plot(ADC_Groups_P2(:,k), TO_Groups_P2,'o', 'markersize', 6);
% h = plot(input_data_P2(:,k), output_data_P2,'o', 'markersize', 6);
set(h, 'MarkerFaceColor', h.Color);
xlabel('ADC Groups P2');
ylabel('TO Groups P2');
end
figure; hold on;
for k = 1:5
% h = plot(ADC_Groups_P2(:,k), TO_Groups_P2,'o', 'markersize', 6);
h = plot(input_data_P2(:,k), output_data_P2,'o', 'markersize', 6);
set(h, 'MarkerFaceColor', h.Color);
xlabel('input data P2');
ylabel('output data P2');
end
%legend({'TO Groups P1'});
%legend({'1', '2', '3'});
hold off;


%% 训练
Group_Data_Count = 5;
iter = 1;
avg_error = 0;
avg_error_P1 = 0;
avg_error_P2 = 0;
for k = 1:iter
    Epoch = 150;
%     [errors, para_monitor] = BP_test_v2(input_data_P1, output_data_P1, Epoch);
%     errors_P1 = errors(:,1);
%     errors_P2 = errors(:,2);
%     para_monitor_P1 = para_monitor(:,1);
%     para_monitor_P2 = para_monitor(:,2);
    [errors_P1, para_monitor_P1] = BP_test_v2(input_data_P1, output_data_P1, Epoch);
    estimated_output_data_P1 = Neural_Network_v2(input_data_P1);
    [errors_P2, para_monitor_P2] = BP_test_v2(input_data_P2, output_data_P2, Epoch);
    estimated_output_data_P2 = Neural_Network_v2(input_data_P2);
    

    errors_P1 = sqrt(errors_P1) * (y_bound_P1(2) - y_bound_P1(1)) / DAC_Steps * Max_power;
    errors_P2 = sqrt(errors_P2) * (y_bound_P2(2) - y_bound_P2(1)) / DAC_Steps * Max_power;
    max_errors_P1 = para_monitor_P1;
    max_errors_P2 = para_monitor_P2;
    max_errors_P1 = sqrt(max_errors_P1) * (y_bound_P1(2) - y_bound_P1(1)) / DAC_Steps * Max_power;
    max_errors_P2 = sqrt(max_errors_P2) * (y_bound_P2(2) - y_bound_P2(1)) / DAC_Steps * Max_power;
    error_30_epoch_P1 = errors_P1(30); error_end_P1 = errors_P1(end);
    error_30_epoch_P2 = errors_P2(30); error_end_P2 = errors_P2(end);

    if k == 1        
        figure; hold on;
        for n = 1:Group_Data_Count
        h = plot(input_data_P1(:,n), estimated_output_data_P1,'o', 'markersize', 6);
        set(h, 'MarkerFaceColor', h.Color);
        end
        %set(gca, 'ylim', [0, 1]);
        %legend({'p1', 'p2', 'p3'});
        hold off;
        xlabel('Monitored power (a. u.)');
        ylabel('TO tuning power (a. u.)');
        set(gca,'linewidth', 1, 'fontname', 'Microsoft Yahei', 'fontsize', 18, 'GridLineStyle', '--');
        
        figure; hold on;
        for n = 1:5
        h = plot(input_data_P2(:,n), estimated_output_data_P2,'o', 'markersize', 6);
        set(h, 'MarkerFaceColor', h.Color);
        end
        %set(gca, 'ylim', [0, 1]);
        %legend({'p1', 'p2', 'p3'});
        hold off;
        xlabel('Monitored power (a. u.)');
        ylabel('TO tuning power (a. u.)');
        set(gca,'linewidth', 1, 'fontname', 'Microsoft Yahei', 'fontsize', 18, 'GridLineStyle', '--');
        
    
        figure;
        semilogy(1:Epoch, errors_P1, 1:Epoch, max_errors_P1, 'linewidth', 2);

        title([num2str(error_30_epoch_P1) ' mW, ' num2str(error_end_P1) ' mW, ' num2str(max_errors_P1(end)) ' mW']);
        set(gca, 'ylim', [0.2 40]);
        xlabel('Epoch');
        ylabel('Estimation error (mW)');
        legend({'Average Error P1', 'Maximal Error P1'});
        set(gca,'linewidth', 1, 'fontname', 'Microsoft Yahei', 'fontsize', 18, 'GridLineStyle', '--');
        width=800;%宽度，像素数
        height=500;%高度
        left=200;%距屏幕左下角水平距离
        bottem=100;%距屏幕左下角垂直距离
        set(gcf,'position',[left,bottem,width,height])     
        
         figure;
        semilogy(1:Epoch, errors_P2, 1:Epoch, max_errors_P2, 'linewidth', 2);

        title([num2str(error_30_epoch_P2) ' mW, ' num2str(error_end_P2) ' mW, ' num2str(max_errors_P2(end)) ' mW']);
        set(gca, 'ylim', [0.2 40]);
        xlabel('Epoch');
        ylabel('Estimation error (mW)');
        legend({'Average Error P2', 'Maximal Error P2'});
        set(gca,'linewidth', 1, 'fontname', 'Microsoft Yahei', 'fontsize', 18, 'GridLineStyle', '--');
        width=800;%宽度，像素数
        height=500;%高度
        left=200;%距屏幕左下角水平距离
        bottem=100;%距屏幕左下角垂直距离
        set(gcf,'position',[left,bottem,width,height])   
    end
    avg_error_P1 = avg_error_P1 + error_end_P1;
    avg_error_P2 = avg_error_P2 + error_end_P2;

end
    fprintf('Average Error_P1 = %.2f\n\n\n',avg_error_P1/iter);
    fprintf('Average Error_P2 = %.2f\n\n\n',avg_error_P2/iter);
    
function [errors, para_monitor] = BP_test_v2(input_data, output_data, Epoch) %para_monitor = max_errors;
    global v w gama  theta Hide_Layer_Cells Output_Data_Count n_max
    Input_Data_Count = 5;
    Output_Data_Count = 1;

    v = rand(Input_Data_Count, Hide_Layer_Cells);%生成5x5的0到1之间的随机数
    w = rand(Hide_Layer_Cells, Output_Data_Count);
    gama = rand(Hide_Layer_Cells, 1);
    theta = rand(Output_Data_Count, 1);
    b = zeros(Hide_Layer_Cells, 1);
    e = zeros(Hide_Layer_Cells, 1);
    
    %%  Neural Network
    Group_Count = length(output_data); %暂定为26个数
    errors = zeros(Epoch,Output_Data_Count);
    max_errors = zeros(Epoch,Output_Data_Count);
        random_indices = randperm(Group_Count);  % 随机顺序测试这26个点
    for epoch = 1:Epoch
        n = n_max/epoch^(1/2);
        for k = 1:Group_Count
            %% 如何构建随机数
%             data_index = mod(k*7,Group_Count)+1;
            data_index = random_indices(k); %随机选取某组5点数据
            x = input_data(data_index,:).';
            y = output_data(data_index, :).';
            y_ = zeros(size(y));
            g = zeros(size(y));
            
            for h = 1:Hide_Layer_Cells
                b(h) = sigmoid(sum(x .* v(:,h)) - gama(h));
            end
            for j = 1:Output_Data_Count
                y_(j) = sigmoid(sum(b .* w(:,j)) - theta(j));
                %更新权重
                g(j) = y_(j)*(1-y_(j))*(y(j)-y_(j));
            end
            
            for h = 1:Hide_Layer_Cells
                e(h) = b(h)*(1-b(h))*sum(w(h,:).'.*g);
            end
            
            for j = 1:Output_Data_Count
                w(:, j) = w(:, j) + n*g(j)*b;
            end
            
            for h = 1:Hide_Layer_Cells
                v(:,h) = v(:,h) + n*e(h)*x;
                gama(h) = gama(h)-n*e(h);
            end
            theta = theta-n*g;
        end
        estimated_output_data = Neural_Network_v2(input_data);
        for j = 1:Output_Data_Count
            errors(epoch, j) = mean((output_data(:, j)-estimated_output_data(:, j)).^2);
            max_errors(epoch, j) = max((output_data(:, j)-estimated_output_data(:, j)).^2);
        end
    end
    
    para_monitor = max_errors;
end

function [errors, para_monitor] = BP_test_v3(input_data, output_data, Epoch) %para_monitor = max_errors;
    global v w gama  theta Hide_Layer_Cells Output_Data_Count
    %Group_Data_Count = 3;
    Input_Data_Count = 5;
    Output_Data_Count = 1;
    
    %Hide_Layer_Cells = 3;
    Hide_Layer_Cells = 5;
    n_max = 10;        % learning rate

    v = rand(Input_Data_Count, Hide_Layer_Cells);%生成5x5的0到1之间的随机数
    w = rand(Hide_Layer_Cells, Output_Data_Count);
    gama = rand(Hide_Layer_Cells, 1);
    theta = rand(Output_Data_Count, 1);
    b = zeros(Hide_Layer_Cells, 1);
    e = zeros(Hide_Layer_Cells, 1);
  %%  Neural Network
    Group_Count = length(output_data); %暂定为26个数
    errors = zeros(Epoch,Output_Data_Count);
    max_errors = zeros(Epoch,Output_Data_Count);
        random_indices = randperm(Group_Count);  % 随机顺序测试这些点
    for epoch = 1:Epoch
        n = n_max/epoch^(1/2);
        for k = 1:Group_Count
            %% 如何构建随机数
%             data_index = mod(k*7,Group_Count)+1;
            data_index = random_indices(k); %随机选取某组5点数据
            x = input_data(data_index,:).';
            y = output_data(data_index, :).';
            y_ = zeros(size(y));
            g = zeros(size(y));
            
            for h = 1:Hide_Layer_Cells
                b(h) = sigmoid(sum(x .* v(:,h)) - gama(h));
            end
            for j = 1:Output_Data_Count
                y_(j) = sigmoid(sum(b .* w(:,j)) - theta(j));
                %更新权重
                g(j) = y_(j)*(1-y_(j))*(y(j)-y_(j));
            end
            
            for h = 1:Hide_Layer_Cells
                e(h) = b(h)*(1-b(h))*sum(w(h,:).'.*g);
            end
            
            for j = 1:Output_Data_Count
                w(:, j) = w(:, j) + n*g(j)*b;
            end
            
            for h = 1:Hide_Layer_Cells
                v(:,h) = v(:,h) + n*e(h)*x;
                gama(h) = gama(h)-n*e(h);
            end
            theta = theta-n*g;
        end
        estimated_output_data = Neural_Network_v2(input_data);
        for j = 1:Output_Data_Count
            errors(epoch, j) = mean((output_data(:, j)-estimated_output_data(:, j)).^2);
            max_errors(epoch, j) = max((output_data(:, j)-estimated_output_data(:, j)).^2);
        end
    end
    
    para_monitor = max_errors;
end

function output_data = Neural_Network_v2(input_data)
    global v w gama  theta Hide_Layer_Cells Output_Data_Count
    b = zeros(Hide_Layer_Cells, 1);
    Group_Count = max(size(input_data));
    output_data = zeros(Group_Count,Output_Data_Count);
    for k = 1:Group_Count
        x = input_data(k,:).';

        for h = 1:Hide_Layer_Cells
            b(h) = sigmoid(sum(x .* v(:,h)) - gama(h));
        end
        for j = 1:Output_Data_Count
            output_data(k, j) = sigmoid(sum(b .* w(:,j)) - theta(j));
        end
    end
end

function y = sigmoid(x)
    y = 1./(1+exp(-x));
end

function output = Get_ADC(index_P1,index_P2)
    global ADC_Set DAC_Steps
    
    rand_span = randi(floor([DAC_Steps / 3, DAC_Steps / 3 * 2]));
    if index_P1 < 1
        index_P1 = index_P1 + rand_span;
    elseif index_P1 > DAC_Steps
        index_P1 = index_P1 - rand_span;
    end
    
    rand_span = randi(floor([DAC_Steps / 3, DAC_Steps / 3 * 2]));
    if index_P2 < 1
        index_P2 = index_P2 + rand_span;
    elseif index_P2 > DAC_Steps
        index_P2 = index_P2 - rand_span;
    end
    
    output = ADC_Set(index_P1,index_P2);
end

