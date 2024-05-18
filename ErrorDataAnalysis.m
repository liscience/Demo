%%遍历绘图
%%机动执行阈下轨道偏差演化终端绘图程序
clc
clear 
close all
% 设置文件夹路径
folder_path = 'ErrorData';
figure_folder = 'ErrorDataFigure';
Processfolder_path = 'ErrorDataProcess';
%数据格式：
%1-6列为终端位置速度误差（EJ2000）
%7-12列为第1次修正时刻的位置速度误差（EJ2000）
%13-18列为第2次修正时刻的位置速度误差（EJ2000）
%第19-21列为第1次修正脉冲矢量
%第22-24列为第2次修正脉冲矢量
%第25列为第1次修正是否施加的标识
%第26列为第2次修正是否施加的标识
%27-32列为终端位置速度误差（MJ2000）
%33-38列为终端位置速度误差（月心LVLH）

% 创建图像保存文件夹
if ~exist(figure_folder, 'dir')
    mkdir(figure_folder);
end

% 获取文件夹中所有txt文件列表
file_list = dir(fullfile(folder_path, '*.txt'));

% 循环处理每个txt文件
for i = 1:length(file_list)
    % 读取txt文件
    i=20;
    file_name = file_list(i).name;
    file_path = fullfile(folder_path, file_name);
    data = load(file_path);
    % 绘制三维散点图
    threshold = 10000;
    rows_to_remove1 = any(vecnorm(data(:,22:24),2,2)>threshold, 2);

    data = data(~rows_to_remove1, :);
    %%
%     % 假设要去除的行号存储在 remove_indices 中
% remove_indices = [362, 800, 2609]; % 示例中假设要去除第5、10和15行的数据
% 
% % 使用逻辑索引删除指定行
% data(remove_indices, :) = [];
    % 获取文件名的后四个元素
    [~, file_title, ~] = fileparts(file_name);
    file_elements = strsplit(file_title, '-');
    file_suffix = strjoin(file_elements(end-3:end), '-');
    
    % 保存为对应名称的dat文件
    dat_file_name = fullfile(Processfolder_path, strcat(file_title, '.dat'));
    save(dat_file_name, 'data', '-ascii');
    %%
    % 绘制三维散点图
    figure;
    scatter3(data(:,33), data(:,34), data(:,35), '.');
    title(['终端位置偏差分布(月心LVLH)', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');   
%     hold on
%     mu = mean(data(:,33:35));   
%     sigma = cov(data(:,33:35));
%     [xx,yy,zz] = plot3sigma3D(mu,sigma,100);%绘制椭球
    % 保存三维散点图
    fig_name = fullfile(figure_folder, [file_title, '_3D.fig']);
    saveas(gcf, fig_name);
    %%
    % 绘制三个二维投影图
    figure;
    subplot(1, 3, 1);
    scatter(data(:,33), data(:,34), '.');
    title(['X-Y平面投影', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    
    subplot(1, 3, 2);
    scatter(data(:,33), data(:,35), '.');
    title(['X-Z平面投影', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Z/m');
    
    subplot(1, 3, 3);
    scatter(data(:,34), data(:,35), '.');
    title(['Y-Z平面投影', '-', file_suffix]);
    xlabel('Y/m');
    ylabel('Z/m');
    
    % 保存二维投影图
    fig_name = fullfile(figure_folder, [file_title, '_2D.fig']);
    saveas(gcf, fig_name);
    %%
    % 绘制第1次修正时刻三维散点图
    figure;
    scatter3(data(:,7), data(:,8), data(:,9), '.');
    title(['第1次修正时刻位置偏差分布', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');
%%
    % 绘制第2次修正时刻三维散点图
    figure;
    scatter3(data(:,13), data(:,14), data(:,15), '.');
    title(['第2次修正时刻位置偏差分布', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');
 %%
% 绘制机动跟踪图
positions1 = data(:, 7:9); % 第1次修正时刻的位置
positions2 = data(:, 13:15); % 第2次修正时刻的位置

labels1 = data(:, 25); % 第1次修正是否执行的标志
labels2 = data(:, 26); % 第2次修正是否执行的标志
figure;
% 绘制第一个子图窗，显示第一组位置点
subplot(1, 2, 1);
hold on;
% 绘制第一组位置点
red_indices = labels1 == 1; % 第一组标识为1的索引
blue_indices = ~red_indices; % 第一组标识为0的索引
scatter3(positions1(red_indices, 1), positions1(red_indices, 2), positions1(red_indices, 3), 'r.'); % 标记为红色
scatter3(positions1(blue_indices, 1), positions1(blue_indices, 2), positions1(blue_indices, 3), 'b.'); % 标记为蓝色
hold off;
% 设置图例和标签
legend('执行机动（红色）', '不执行机动（蓝色）');
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
title('第1次修正时刻机动执行情况分布图');
view(3); % 显示三维视图

% 绘制第二个子图窗，显示第二组位置点
subplot(1, 2, 2);
hold on;
% 绘制第二组位置点
red_indices = labels2 == 1; % 第二组标识为1的索引
blue_indices = ~red_indices; % 第二组标识为0的索引
scatter3(positions2(red_indices, 1), positions2(red_indices, 2), positions2(red_indices, 3), 'r.'); % 标记为红色
scatter3(positions2(blue_indices, 1), positions2(blue_indices, 2), positions2(blue_indices, 3), 'b.'); % 标记为蓝色
hold off;
% 设置图例和标签
legend('执行机动（红色）', '不执行机动（蓝色）');
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
title('第2次修正时刻机动执行情况分布图');
view(3); % 显示三维视图
    
%%
% 绘制终端时刻空间点分布特性
figure;
% 创建空间点坐标向量
x = data(:, 33);
y = data(:, 34);
z = data(:, 35);

% 创建颜色向量
% 初始化颜色向量为全蓝色
colors = repmat([0, 0, 1], size(data, 1), 1); 
red_indices = data(:, 25) == 1 & data(:, 26) == 1;%对应两次机动都执行
green_indices = data(:, 25) == 1 & data(:, 26) == 0;%对应第1次机动执行，第2次机动不执行
yellow_indices = data(:, 25) == 0 & data(:, 26) == 1;%对应第1次机动不，第2次机动执行

% 设置不同颜色的空间点
colors(red_indices, :) = repmat([1, 0, 0], sum(red_indices), 1); % 红色
colors(green_indices, :) = repmat([0, 1, 0], sum(green_indices), 1); % 绿色
colors(yellow_indices, :) = repmat([1, 1, 0], sum(yellow_indices), 1); % 黄色

% 绘制三维散点图
scatter3(x, y, z, 7, colors, 'filled');
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
title('终端位置误差与机动分布');
% 添加图例
hold on; % 保持图形，以便添加多个图例项
scatter3(nan, nan, nan, 20, [1 0 0], 'filled'); % 添加红色点的图例项
scatter3(nan, nan, nan, 20, [0 1 0], 'filled'); % 添加绿色点的图例项
scatter3(nan, nan, nan, 20, [1 1 0], 'filled'); % 添加黄色点的图例项
hold off; % 释放图形

% 设置图例
% legend('红色', '绿色', '黄色', '蓝色');
% 设置图例
% legend('两次机动均执行', '仅执行第1次机动', '仅执行第2次机动', '两次机动均不执行');
%%
data1 = data(:,19:21);
data2 = data(:,7:9);
data3 = vecnorm(data1,2,2);
data4 = data(:,22:24);
data5 = data(:,13:15);
data6 = vecnorm(data4,2,2);
data3(data3 == 0) = NaN;
data6(data6 == 0) = NaN;
threshold1 = 10;threshold2 = 10; % 假设阈值
points_below_threshold = data2(data3 < threshold1, :);
alpha_shape = alphaShape(points_below_threshold(:,1), points_below_threshold(:,2), points_below_threshold(:,3));
points_below_threshold1 = data5(data6 < threshold2, :);
alpha_shape1 = alphaShape(points_below_threshold1(:,1), points_below_threshold1(:,2), points_below_threshold1(:,3));

% 绘制等高线图
figure;
scatter3(data2(:,1), data2(:,2), data2(:,3), 7, data3, 'filled','DisplayName', '机动修正量大小/(m/s)'); % 绘制散点图
hold on;
plot(alpha_shape,  'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.7,'DisplayName', '机动执行阈包络面'); % 绘制包络面
colorbar; % 添加颜色条
xlabel('X轴方向位置/m');
ylabel('Y轴方向位置/m');
zlabel('Z轴方向位置/m');
legend('Location', 'best'); % 添加图例说明，位置自动选择最佳位置
figure;
scatter3(data5(:,1), data5(:,2), data5(:,3), 10, data6, 'filled','DisplayName', '机动修正量大小/(m/s)'); % 绘制散点图
hold on;
plot(alpha_shape1, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.7,'DisplayName', '机动执行阈包络面'); % 绘制包络面
colorbar; % 添加颜色条
xlabel('X轴方向位置/m');
ylabel('Y轴方向位置/m');
zlabel('Z轴方向位置/m');
legend('Location', 'best'); % 添加图例说明，位置自动选择最佳位置
    
    close all
end
