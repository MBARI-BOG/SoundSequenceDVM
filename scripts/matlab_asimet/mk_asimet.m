
function mk_asimet(deployment,d_file)
% Create asimet qc graphs for M1
% deployment='202008'; d_file='20200825';
% mk_asimet(deployment,d_file)
close all; 
[yr,mo,da]=datevec(now);
disp(['Begin M1 Asimet processing ',datestr(now)])
clear mo da

% for running on the virtual machine bogdata:
path_prefix = '/mbari'
% for running locally on a mac (must have loaded ssdsdata, Diatom):
% path_prefix = '/Volumes'

% Plot save location:
plot_dir = [path_prefix,'/Diatom/users/kpitz/mooring_test']

try
    %fname=['\\atlas\ssdsdata\deployments\m1\',deployment,'\m1_asimetlwr_',d_file,'.nc'];
    fname=[path_prefix,'/ssdsdata/deployments/m1/202207/m1_asimetlwr_20220718.nc'];
    fid=fopen(fname);
    if fid>-1
        %set figure filename
        filename = [plot_dir, '/asimetlwr.jpg'];
        dat=ncread(fname,'LW_flux');
        utc=ncread(fname,'esecs'); 
        sdn=utc2sdn(utc);
        LW=dat(:);
        idx=find(LW==-99999);
        LW(idx)=NaN;
        % create subplots week and 24 hr
        subplot(2,1,2)
        idx=find(sdn>floor(now)-7);
        plot(sdn(idx),LW(idx),'b-');
        grid on; box on;
        set(gca,'YLim',[250 450],'fontsize',6);
        datetick('x',6,'keeplimits','keepticks');
        title(['7 Day M1 Asimet LW Flux ',num2str(yr)],'fontsize',8);
        %     ylabel('Long Wave Flux (watts m^-^2)','fontsize',8);
        ylabel('Long Wave Flux (watts m^-^2)','fontsize',8);
        xlabel('Date GMT','fontsize',8);
        % 24 hr
        subplot(2,1,1)
        idx=find(sdn>max(sdn)-1);
        plot(sdn(idx),LW(idx),'r-');
        grid on; box on;
        set(gca,'YLim',[250 450]);
        datetick('x',15,'keeplimits','keepticks');
        title(['24 Hour M1 Asimet LW Flux ',datestr(floor(max(sdn(idx))))],'fontsize',8);
        %     ylabel('Long Wave Flux (watts m^-^2)','fontsize',8);
        ylabel('Long Wave Flux (watts m^-^2)','fontsize',8);
        xlabel('Time Pacific Standard Time','fontsize',8);
        set(gca,'fontsize',6);
        set(gcf,'Position',[304 122 762 736],'PaperPosition',[0.25 0.25 6 4]);
        % print -djpeg95 -painters -r125 \\atlas\diatom\oasis\tests\asimetlwr.jpg
        % print -djpeg95 -painters -r125 /Users/kpitz/github/MBARI-BOG/mooring_processing/matlab/asimet/test_outputs/asimetlwr.jpg
        print('-djpeg95', '-painters', '-r125',sprintf(filename))
    else
        disp('Error netcdf file not opened');
        
    end
    
    close all; %clear all;
    [yr,mo,da]=datevec(now);
    clear mo da
        %fname=['\\atlas\ssdsdata\deployments\m1\',deployment,'\m1_asimetswr_',d_file,'.nc'];
        fname=[path_prefix,'/ssdsdata/deployments/m1/202207/m1_asimetswr_20220718.nc'];
    fid=fopen(fname);
    if fid>-1
        fclose(fid);
        %set figure filename
        filename = [plot_dir, '/asimetswr.jpg'];
        %     ncid=netcdf('\\atlas\ssdsdata\deployments\m1\current_netcdfs\asimetswr.nc','read');
        utc=ncread(fname,'esecs'); %ncid{'esecs'};
        dat=ncread(fname,'SW_flux'); %ncid{'esecs'};
        sdn=utc2sdn(utc);
        SW=dat(:);
        SW(SW==-99999)=NaN;
        % create subplots week and 24 hr
        subplot(2,1,2)
        idx=find(sdn>floor(now)-7);
        plot(sdn(idx),SW(idx),'b-');
        grid on; box on;
        %     set(gca,'YLim',[250 450]);
        datetick('x',6,'keeplimits','keepticks');
        title(['7 Day M1 Asimet SW Flux ',num2str(yr)],'fontsize',8);
        ylabel('Short Wave Flux (watts m^-^2)','fontsize',8);
        xlabel('Date GMT','fontsize',8);
        set(gca,'fontsize',6);
        % 24 hr
        subplot(2,1,1)
        idx=find(sdn>max(sdn)-1);
        plot(sdn(idx),SW(idx),'r-');
        grid on; box on;
        %     set(gca,'YLim',[250 450]);
        datetick('x',15,'keeplimits','keepticks');
        title(['24 Hour M1 Asimet SW Flux ',datestr(floor(max(sdn(idx))))],'fontsize',8);
        ylabel('Short Wave Flux (watts m^-^2)','fontsize',8);
        xlabel('Time Pacific Standard Time','fontsize',8);
        set(gca,'fontsize',6);
        set(gcf,'Position',[304 122 762 736],'PaperPosition',[0.25 0.25 6 4]);
        % print -djpeg95 -painters -r125 /Users/kpitz/github/MBARI-BOG/mooring_processing/matlab/asimet/test_outputs/asimetswr.jpg
        print('-djpeg95', '-painters', '-r125',sprintf(filename))
    else
        disp('Error netcdf file not opened');
        
    end
    disp(['Completed M1 Asimet processing ',datestr(now)])
    
    close all; clear all;
    % exit
catch ME
    disp(ME.identifier);
    %     rethrow(ME);
    % exit
end %try