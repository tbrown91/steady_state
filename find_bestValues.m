file = 'test3';

data = importdata(file,'\t',1);

appended_data = data.data;

appended_data(any(isnan(appended_data),2),:) = [];

unique_data = unique(appended_data,'rows');

sorted_data = sortrows(unique_data,1);

valuesFile = fopen('best_values.txt','w');

for i = 1:40000
    for j = 2:3
       fprintf(valuesFile,'%f\t',sorted_data(i,j));
    end
    fprintf(valuesFile,'%f\n',sorted_data(i,4));
end
fclose('all');