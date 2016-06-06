%% Save all the data in one matrix
clear
data_dirs = dir('.');
count = 0;
folder = {};
for i = 3:length(data_dirs)
    if data_dirs(i).isdir == 1
        count = count + 1;
        folder{count} = data_dirs(i).name;
        folder_name = char(folder(count));
        F = dir(strcat(folder_name, '/*.dat'));
        num_files = length(F);
        a = length(csvread(strcat(folder_name, '/', F(1).name)));
        for j = 1:num_files
            data(((count-1)*a+1):count*a,j) = csvread(strcat(folder_name, '/', F(j).name));
        end
    end
end

%% Plot four subplots, each one is the comparision of different size of lattice
subplot(2,2,1);
for i = 1:count
    scatter(data(1,:), data(3 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel('Average Energy');
hold off;

subplot(2,2,2);
for i = 1:count
    scatter(data(1,:), data(4 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel();
hold off;

subplot(2,2,3);
for i = 1:count
    scatter(data(1,:), data(5 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel('Average Magnetization');
hold off;

subplot(2,2,4);
for i = 1:count
    scatter(data(1,:), data(6 + a*(i-1),:));
    hold on;
end
legend(folder);
xlabel('Temperature');
ylabel();
hold off;
saveas(gcf,'plot','png')