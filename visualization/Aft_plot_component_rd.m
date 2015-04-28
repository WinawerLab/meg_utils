function Aft_plot_component_rd(ft_comp,channels,layout,trl, shift, fs,Nsec,Position)
%
%  Aft_plot_component_rd(ft_comp,channels,layout,trl, shift, fs,Nsec,Position)
%
% Aft_plot_component will plot the first Nsec seconds, power spectra , event related average 
% and topography of specified components (channels) from a fieldtrip structure returned by ft_componentanalysis.
% 
% This function needs the trl array data for trial structure and assumes
% the component data continuous from sample 1 of trial 1 through the last
% sample of the last trial
% 
% shift is the difference between the data and trl events (if you didn't
% change the data start position use shift =0)
%
%
if nargin<6, fs = 1000; end
if nargin<7, Nsec = 3; end

    run_this = ft_comp.trial{1}(channels,:);
    [M,N] = size(run_this);
    
    this_topo = ft_comp.topo;
    
    chanX = layout.pos(1:size(this_topo,1),1);
    chanY = layout.pos(1:size(this_topo,1),2);
   
    
    this_size = (length(run_this)-trl(end,1)+shift);
    EventRelatedX_axis = trl(1,3):this_size-abs(trl(1,3))-1;
    avg = zeros(M,this_size);

    for i = 1:length(trl)
         avg = avg+run_this(:,trl(i,1)-shift:trl(i,1)-shift+this_size-1);
     end
     avg = avg./i;
     
    
     Pxx = zeros(1001,M);
     MaxFview = 50;
     
     if std(run_this(1,1:round(Nsec*fs)))<1
         norm = 1;
         warning('Must normalize data to units of 1e-13');
     else
         norm = 0;
     end
     
     cnt = 1;
     for i = 1:M
%        [Pxx,F] = pwelch(run_this(i,:),[],[],2000,fs);
       [Pxx,F] = pwelch(run_this(i,:),[],[],[],fs); % rd
       Psel = 1:find(F<MaxFview,1,'last');
     
       subplot(M,3,cnt);cnt  = cnt+1;
       
       if norm
            %plot(run_this(i,1:round(Nsec*fs))./max(run_this(i,1:round(Nsec*fs)))*100);
            plot(run_this(i,1:round(Nsec*fs))./1e-13);
       else
           plot(run_this(i,1:round(Nsec*fs)));
       end
       title(ft_comp.label(channels(i)));
       subplot(M,3,cnt);cnt  = cnt+1;
       
       if norm
            plot(EventRelatedX_axis,avg(i,:)./1e-13);
            %plot(trl(1,3):trl(1,2)-trl(1,1)+trl(1,3),avg(i,:)./max(abs(avg(i,:)))*100);
       else
            plot(EventRelatedX_axis,avg(i,:));
       end
       xlim([EventRelatedX_axis(1) EventRelatedX_axis(end)]);
       
       subplot(M,3,cnt);cnt  = cnt+1; 
       plot(F(Psel),Pxx(Psel)./max(Pxx(Psel))*100);
      
       chanZ = this_topo(:,channels(i));

     ax = axis;
     hpos = round(ax(2)*2/3);
     vpos = round(ax(4)*2/3);
     width = hpos/1.5;
     height = vpos/1.5;
        ft_plot_topo(chanX, chanY, chanZ, 'mask', ...
            layout.mask, 'interplim', 'mask', 'outline', ...
            layout.outline, 'tag', 'topography',...
            'hpos',hpos,...
            'vpos',vpos,...
            'width',width,...
            'height',height);
     end
     if exist('Position')
         set(gcf,'Position',Position);
     end
end

