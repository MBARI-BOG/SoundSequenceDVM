% Look at short wave radiation over duration of CN19S experiment
% 05-29-19 to 06-06-19
% KP

%% Set parameters
%fname = '../data/asimet/m1_asimetlwr_20180806.nc'
fname = '/Users/kpitz/github/MBARI-BOG/CN19S_12S/data/asimet/m1_asimetswr_20180806.nc'
% Plots/Formatted data save location:
plot_dir = '/Users/kpitz/github/MBARI-BOG/CN19S_12S/data/asimet/'

%% look at the netcdf file
ncdisp(fname)
% esecs
% units     = 'seconds since 1970-01-01 00:00:00'
% SW_flux
% units         = 'watts per square meter'

%% Import data
% Get information about what is in the nc file
tmp = ncinfo([fname]);
% Varable names
vars = {tmp.Variables.Name};
for i = 1:length(vars)
    % current variable
    cvar = char(vars(i));
    M1.(cvar) = ncread([fname],cvar);
end
M1.dt = datetime(1970,1,1) + seconds(M1.esecs);
% Check date makes sense
datestr(M1.dt);

%M1.month =  month( M1.dt )
%M1.year =  year( M1.dt )

%plot(M1.dt, M1.SW_flux)
%M = mean(M1, year) 

%% More broadly limit at first before interpolation

%find indices of data in time frame
%ix=find([M1.dt]>=datetime(2019,5,29, 0, 0, 0) & [M1.dt]<=datetime(2019,6,7, 0, 0, 0));
% GMT -7 = PDT
% midnight PDT = 7am (07:00) GMT
ix=find([M1.dt]>=datetime(2019,5,29, 7, 0, 0) & [M1.dt]<=datetime(2019,6,7, 7, 0, 0));

%remove a dimension from flux data
flux_lim = squeeze(M1.SW_flux(:,:,:,:,:,:,1,:));

%create new structure with limited data:
M1lim = struct('dt', M1.dt(ix), 'SW_flux', flux_lim(ix), 'esecs', M1.esecs(ix))

%shift time zone (subtract 7 hours GMT to PDT)
M1lim.dt = M1lim.dt - hours(7)


% get time of day
% T = timeofday(DT) returns a duration array whose values equal the elapsed time since midnight for each element in DT.
M1lim.timeofday = timeofday(M1lim.dt)
%M1lim.day = day(M1lim.dt)


% plot through time
scatter(M1lim.dt, M1lim.SW_flux)
hold on 
plot(M1lim.dt, M1lim.SW_flux)
% plot by time of day
%scatter(M1lim.timeofday, M1lim.SW_flux)
% color by time (seconds since 1970)
scatter(M1lim.timeofday, M1lim.SW_flux,25, M1lim.esecs, 'filled')
colorbar


% % Try with timetable
% TT = timetable(M1lim.dt, M1lim.SW_flux)
% TT2 = retime(TT, 'minutely', 'spline')
% dt = minutes(15);
% TT3 = retime(TT,'regular','linear','TimeStep',dt)
% TT4 = retime(TT2,'regular','linear','TimeStep',dt)
% 
% TT2.timeofday  = timeofday(TT2.Time)
% 
% % plot differences
% scatter(TT, "Time", "Var1", 'filled', 'MarkerFaceColor','red')
% hold on
% scatter(TT2, "Time", "Var1", "filled",'MarkerFaceColor', 'black')
% hold on
% scatter(TT3, "Time", "Var1", "filled",'MarkerFaceColor', 'green')
% hold on
% scatter(TT4, "Time", "Var1", "filled",'MarkerFaceColor', 'blue')
% 
% scatter(TT2, "timeofday", "Var1")

%M1lim has 648 points

%Interpolate not assuming uniform times, but with Signal Processing Toolbox
%M1lim.interp = resample(M1lim.SW_flux, M1lim.dt, 1/(15*60));   %1/(15*60) Hz -> 15 minutes apart
[y,ty] = resample(M1lim.SW_flux, M1lim.dt, 1/(15*60));   %1/(15*60) Hz -> 15 minutes apart


scatter(ty,y)
hold on
scatter(M1lim.dt, M1lim.SW_flux,25, M1lim.esecs, 'filled')

%create new structure with interpolated data:
M1interp = struct('dt', ty, 'SW_flux', y)
M1interp.timeofday = timeofday(M1interp.dt)
M1interp.day = day(M1interp.dt)

% plot by time of day
scatter(M1interp.timeofday, M1interp.SW_flux)


T = struct2table(M1interp);
% take mean over time of day, and SD
tblstats1 = grpstats(T,"timeofday", ["mean", "std"])

tblstats1.Properties.VariableNames
Tnew= removevars(tblstats1,{'mean_dt','std_dt','mean_day','std_day' })
plot(Tnew, "timeofday", "mean_SW_flux")

% Add in standard deviation to plot:
curve1 = Tnew.mean_SW_flux + Tnew.std_SW_flux
curve2 = Tnew.mean_SW_flux - Tnew.std_SW_flux

plot(Tnew, "timeofday", "mean_SW_flux")
hold on
scatter(Tnew.timeofday, curve1,25,'filled', 'red')
hold on 
scatter(Tnew.timeofday, curve2,25,'filled', 'black')

% Export data to csv file for R plotting
writetable(Tnew, '/Users/kpitz/github/MBARI-BOG/CN19S_12S/data/asimet/Interp_Mean_Values.csv')


%% Combine plots

% plot by time of day
scatter(M1lim.timeofday, M1lim.SW_flux)
hold on
% color by time (seconds since 1970)
scatter(M1lim.timeofday, M1lim.SW_flux,25, M1lim.esecs, 'filled')
hold on
colorbar


plot(Tnew, "timeofday", "mean_SW_flux")
hold on
scatter(Tnew.timeofday, curve1,25,'filled', 'red')
hold on 
scatter(Tnew.timeofday, curve2,25,'filled', 'black')
hold on
scatter(M1lim.timeofday, M1lim.SW_flux,25, M1lim.esecs, 'filled')
saveas(gcf,'/Users/kpitz/github/MBARI-BOG/CN19S_12S/data/asimet/Interpolated_mean_chart.png')