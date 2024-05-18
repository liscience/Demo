%%������ͼ
%%����ִ�����¹��ƫ���ݻ��ն˻�ͼ����
clc
clear 
close all
% �����ļ���·��
folder_path = 'ErrorData';
figure_folder = 'ErrorDataFigure';
Processfolder_path = 'ErrorDataProcess';
%���ݸ�ʽ��
%1-6��Ϊ�ն�λ���ٶ���EJ2000��
%7-12��Ϊ��1������ʱ�̵�λ���ٶ���EJ2000��
%13-18��Ϊ��2������ʱ�̵�λ���ٶ���EJ2000��
%��19-21��Ϊ��1����������ʸ��
%��22-24��Ϊ��2����������ʸ��
%��25��Ϊ��1�������Ƿ�ʩ�ӵı�ʶ
%��26��Ϊ��2�������Ƿ�ʩ�ӵı�ʶ
%27-32��Ϊ�ն�λ���ٶ���MJ2000��
%33-38��Ϊ�ն�λ���ٶ�������LVLH��

% ����ͼ�񱣴��ļ���
if ~exist(figure_folder, 'dir')
    mkdir(figure_folder);
end

% ��ȡ�ļ���������txt�ļ��б�
file_list = dir(fullfile(folder_path, '*.txt'));

% ѭ������ÿ��txt�ļ�
for i = 1:length(file_list)
    % ��ȡtxt�ļ�
    i=20;
    file_name = file_list(i).name;
    file_path = fullfile(folder_path, file_name);
    data = load(file_path);
    % ������άɢ��ͼ
    threshold = 10000;
    rows_to_remove1 = any(vecnorm(data(:,22:24),2,2)>threshold, 2);

    data = data(~rows_to_remove1, :);
    %%
%     % ����Ҫȥ�����кŴ洢�� remove_indices ��
% remove_indices = [362, 800, 2609]; % ʾ���м���Ҫȥ����5��10��15�е�����
% 
% % ʹ���߼�����ɾ��ָ����
% data(remove_indices, :) = [];
    % ��ȡ�ļ����ĺ��ĸ�Ԫ��
    [~, file_title, ~] = fileparts(file_name);
    file_elements = strsplit(file_title, '-');
    file_suffix = strjoin(file_elements(end-3:end), '-');
    
    % ����Ϊ��Ӧ���Ƶ�dat�ļ�
    dat_file_name = fullfile(Processfolder_path, strcat(file_title, '.dat'));
    save(dat_file_name, 'data', '-ascii');
    %%
    % ������άɢ��ͼ
    figure;
    scatter3(data(:,33), data(:,34), data(:,35), '.');
    title(['�ն�λ��ƫ��ֲ�(����LVLH)', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');   
%     hold on
%     mu = mean(data(:,33:35));   
%     sigma = cov(data(:,33:35));
%     [xx,yy,zz] = plot3sigma3D(mu,sigma,100);%��������
    % ������άɢ��ͼ
    fig_name = fullfile(figure_folder, [file_title, '_3D.fig']);
    saveas(gcf, fig_name);
    %%
    % ����������άͶӰͼ
    figure;
    subplot(1, 3, 1);
    scatter(data(:,33), data(:,34), '.');
    title(['X-Yƽ��ͶӰ', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    
    subplot(1, 3, 2);
    scatter(data(:,33), data(:,35), '.');
    title(['X-Zƽ��ͶӰ', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Z/m');
    
    subplot(1, 3, 3);
    scatter(data(:,34), data(:,35), '.');
    title(['Y-Zƽ��ͶӰ', '-', file_suffix]);
    xlabel('Y/m');
    ylabel('Z/m');
    
    % �����άͶӰͼ
    fig_name = fullfile(figure_folder, [file_title, '_2D.fig']);
    saveas(gcf, fig_name);
    %%
    % ���Ƶ�1������ʱ����άɢ��ͼ
    figure;
    scatter3(data(:,7), data(:,8), data(:,9), '.');
    title(['��1������ʱ��λ��ƫ��ֲ�', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');
%%
    % ���Ƶ�2������ʱ����άɢ��ͼ
    figure;
    scatter3(data(:,13), data(:,14), data(:,15), '.');
    title(['��2������ʱ��λ��ƫ��ֲ�', '-', file_suffix]);
    xlabel('X/m');
    ylabel('Y/m');
    zlabel('Z/m');
 %%
% ���ƻ�������ͼ
positions1 = data(:, 7:9); % ��1������ʱ�̵�λ��
positions2 = data(:, 13:15); % ��2������ʱ�̵�λ��

labels1 = data(:, 25); % ��1�������Ƿ�ִ�еı�־
labels2 = data(:, 26); % ��2�������Ƿ�ִ�еı�־
figure;
% ���Ƶ�һ����ͼ������ʾ��һ��λ�õ�
subplot(1, 2, 1);
hold on;
% ���Ƶ�һ��λ�õ�
red_indices = labels1 == 1; % ��һ���ʶΪ1������
blue_indices = ~red_indices; % ��һ���ʶΪ0������
scatter3(positions1(red_indices, 1), positions1(red_indices, 2), positions1(red_indices, 3), 'r.'); % ���Ϊ��ɫ
scatter3(positions1(blue_indices, 1), positions1(blue_indices, 2), positions1(blue_indices, 3), 'b.'); % ���Ϊ��ɫ
hold off;
% ����ͼ���ͱ�ǩ
legend('ִ�л�������ɫ��', '��ִ�л�������ɫ��');
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
title('��1������ʱ�̻���ִ������ֲ�ͼ');
view(3); % ��ʾ��ά��ͼ

% ���Ƶڶ�����ͼ������ʾ�ڶ���λ�õ�
subplot(1, 2, 2);
hold on;
% ���Ƶڶ���λ�õ�
red_indices = labels2 == 1; % �ڶ����ʶΪ1������
blue_indices = ~red_indices; % �ڶ����ʶΪ0������
scatter3(positions2(red_indices, 1), positions2(red_indices, 2), positions2(red_indices, 3), 'r.'); % ���Ϊ��ɫ
scatter3(positions2(blue_indices, 1), positions2(blue_indices, 2), positions2(blue_indices, 3), 'b.'); % ���Ϊ��ɫ
hold off;
% ����ͼ���ͱ�ǩ
legend('ִ�л�������ɫ��', '��ִ�л�������ɫ��');
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
title('��2������ʱ�̻���ִ������ֲ�ͼ');
view(3); % ��ʾ��ά��ͼ
    
%%
% �����ն�ʱ�̿ռ��ֲ�����
figure;
% �����ռ����������
x = data(:, 33);
y = data(:, 34);
z = data(:, 35);

% ������ɫ����
% ��ʼ����ɫ����Ϊȫ��ɫ
colors = repmat([0, 0, 1], size(data, 1), 1); 
red_indices = data(:, 25) == 1 & data(:, 26) == 1;%��Ӧ���λ�����ִ��
green_indices = data(:, 25) == 1 & data(:, 26) == 0;%��Ӧ��1�λ���ִ�У���2�λ�����ִ��
yellow_indices = data(:, 25) == 0 & data(:, 26) == 1;%��Ӧ��1�λ���������2�λ���ִ��

% ���ò�ͬ��ɫ�Ŀռ��
colors(red_indices, :) = repmat([1, 0, 0], sum(red_indices), 1); % ��ɫ
colors(green_indices, :) = repmat([0, 1, 0], sum(green_indices), 1); % ��ɫ
colors(yellow_indices, :) = repmat([1, 1, 0], sum(yellow_indices), 1); % ��ɫ

% ������άɢ��ͼ
scatter3(x, y, z, 7, colors, 'filled');
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
title('�ն�λ�����������ֲ�');
% ���ͼ��
hold on; % ����ͼ�Σ��Ա���Ӷ��ͼ����
scatter3(nan, nan, nan, 20, [1 0 0], 'filled'); % ��Ӻ�ɫ���ͼ����
scatter3(nan, nan, nan, 20, [0 1 0], 'filled'); % �����ɫ���ͼ����
scatter3(nan, nan, nan, 20, [1 1 0], 'filled'); % ��ӻ�ɫ���ͼ����
hold off; % �ͷ�ͼ��

% ����ͼ��
% legend('��ɫ', '��ɫ', '��ɫ', '��ɫ');
% ����ͼ��
% legend('���λ�����ִ��', '��ִ�е�1�λ���', '��ִ�е�2�λ���', '���λ�������ִ��');
%%
data1 = data(:,19:21);
data2 = data(:,7:9);
data3 = vecnorm(data1,2,2);
data4 = data(:,22:24);
data5 = data(:,13:15);
data6 = vecnorm(data4,2,2);
data3(data3 == 0) = NaN;
data6(data6 == 0) = NaN;
threshold1 = 10;threshold2 = 10; % ������ֵ
points_below_threshold = data2(data3 < threshold1, :);
alpha_shape = alphaShape(points_below_threshold(:,1), points_below_threshold(:,2), points_below_threshold(:,3));
points_below_threshold1 = data5(data6 < threshold2, :);
alpha_shape1 = alphaShape(points_below_threshold1(:,1), points_below_threshold1(:,2), points_below_threshold1(:,3));

% ���Ƶȸ���ͼ
figure;
scatter3(data2(:,1), data2(:,2), data2(:,3), 7, data3, 'filled','DisplayName', '������������С/(m/s)'); % ����ɢ��ͼ
hold on;
plot(alpha_shape,  'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.7,'DisplayName', '����ִ���а�����'); % ���ư�����
colorbar; % �����ɫ��
xlabel('X�᷽��λ��/m');
ylabel('Y�᷽��λ��/m');
zlabel('Z�᷽��λ��/m');
legend('Location', 'best'); % ���ͼ��˵����λ���Զ�ѡ�����λ��
figure;
scatter3(data5(:,1), data5(:,2), data5(:,3), 10, data6, 'filled','DisplayName', '������������С/(m/s)'); % ����ɢ��ͼ
hold on;
plot(alpha_shape1, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.7,'DisplayName', '����ִ���а�����'); % ���ư�����
colorbar; % �����ɫ��
xlabel('X�᷽��λ��/m');
ylabel('Y�᷽��λ��/m');
zlabel('Z�᷽��λ��/m');
legend('Location', 'best'); % ���ͼ��˵����λ���Զ�ѡ�����λ��
    
    close all
end
