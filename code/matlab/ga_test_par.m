%% STDP Curve Fitting via Genetic Algorithm
function ga_test_par(numpars, index,expname, seed, params)


%% Importing Experimental Variables
gadata=exp_data(numpars,index,expname);

%% Plotting 
fitting_fig = figure('NumberTitle','off','doublebuffer','on',...
       'Name','Fitting STDP Curve',...
       'Units','normalized','Position',[0.4 0.2 0.55 0.7]);
time_w = min(gadata.stdp_time):0.1:max(gadata.stdp_time);

integr=integration( ...
  time_w, ...
  params(1),params(2),params(3), ...
  params(4),params(5),params(6), ...
  params(7),params(8),params(9), ...
  params(10),params(11)); 


plot(time_w, integr);
hold on;
plot(gadata.stdp_time,gadata.stdp_fun,'*r');
set(gca,'XLim',[min(gadata.stdp_time), max(gadata.stdp_time)]);
hold on;
line([min(gadata.stdp_time) max(gadata.stdp_time)],mean(gadata.stdp_fun)*ones(1,2),'Color',[0 1 0]);
pause(1)
fname = ['curve-',expname '-' num2str(numpars) '-' num2str(index)];
print(gcf,'-dpng', fname);
fname = ['curve-',expname '-' num2str(numpars) '-' num2str(index) '-' num2str(seed)  '_params'];
params_file= fopen(fname,'w');
fprintf(params_file, '%params');
for t = 1:length(params) 
  fprintf(params_file, '%e ' , params(t) );
end
fprintf(params_file, '\n');
fclose(params_file);
fname = ['curve-',expname '-' num2str(numpars) '-' num2str(index) '-' num2str(seed)  '_data'];
data_file= fopen(fname,'w');
for t = 1:length(time_w) 
  fprintf(data_file, '%e %e\n' , time_w(t), integr(t) );
end
fprintf(data_file, '\n');
fclose(data_file);

fname = ['curve-',expname '-' num2str(numpars) '-' num2str(index) '-' num2str(seed)  '_orig_data'];
data_file= fopen(fname,'w');
for t = 1:length(gadata.stdp_time) 
  fprintf(data_file, '%e %e\n' ,gadata.stdp_time(t), gadata.stdp_fun(t) );
end
fprintf(data_file, '\n');
fclose(data_file);

