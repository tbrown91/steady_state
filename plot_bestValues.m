best = fopen('best_values.txt','r');
data = fscanf(best,'%f %f %f',[3 Inf]);
data = data';

%Plot histogram of best values
figure();
subplot(3,1,1);
hist(data(:,1),500);
title('alpha/beta');

subplot(3,1,2);
hist(data(:,2),500);
title('alpha/gamma');

subplot(3,1,3);
hist(data(:,3),500);
title('alpha/delta');

%Plot all values against each other
figure();
plot3(1./data(:,1),1./data(:,2),1./data(:,3),'b.');
xlabel('alpha/beta');
ylabel('alpha/gamma');
zlabel('alpha/delta');
hold on;
grid('on');

figure();
plot(data(:,1),data(:,2),'b.');
xlabel('beta+gamma');
ylabel('beta*gamma');