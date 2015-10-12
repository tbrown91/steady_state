best = fopen('best_values.txt','r');
data = fscanf(best,'%f %f %f',[3 Inf]);
data = data';

%Plot histogram of best values
figure();
subplot(3,1,1);
hist(data(:,1),100);
title('alpha/beta');

subplot(3,1,2);
hist(data(:,2),100);
title('alpha/gamma');

subplot(3,1,3);
hist(data(:,3),100);
title('alpha/delta');

%Plot all values against each other
figure();
plot3(data(:,1),data(:,2),data(:,3),'b.');
xlabel('alpha/beta');
ylabel('alpha/gamma');
zlabel('alpha/delta');
hold on;
grid('on');