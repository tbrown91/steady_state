nuc_simFile = 'simulations_nuc.txt';
cyt_simFile = 'simulations_cyt.txt';
nuc_ssFile = 'steady_statesNuc.txt';
cyt_ssFile = 'steady_statesCyt.txt';

nuc_simData = importdata(nuc_simFile,'\t',0);
cyt_simData = importdata(cyt_simFile,'\t',0);
nuc_ssData = importdata(nuc_ssFile,'\t',0);
cyt_ssData = importdata(cyt_ssFile,'\t',0);

figure();
for i =1:length(nuc_simData)
    subplot(1,2,1);
    plot(0:49,nuc_simData(i,:),'b',0:49,nuc_ssData(i,:),'r');
    xlim([0,10]);
    subplot(1,2,2);
    plot(0:49,cyt_simData(i,:),'b',0:49,cyt_ssData(i,:),'r');
    xlim([0,10]);
    pause(0.1);
end
