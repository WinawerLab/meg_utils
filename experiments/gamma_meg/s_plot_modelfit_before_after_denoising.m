% Visualizes summary statistics computed in s_GAMMA_MEG_analysis.m

project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';


fs                            = 1000;
intertrial_trigger_num        = 10;
which_data_to_visualize       = 9;
save_images                   = true;

meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

condition_names               = {   ...
    'White Noise' ...
    'Binarized White Noise' ...
    'Pink Noise' ...
    'Brown Noise' ...
    'Gratings(0.36 cpd)' ...
    'Gratings(0.73 cpd)' ...
    'Gratings(1.45 cpd)' ...
    'Gratings(2.90 cpd)' ...
    'Plaid'...
    'Blank'};

%% Loop over data sets
cmap = jet(9);
for subject_num = which_data_to_visualize
    
    load_pth    = fullfile(project_pth, subj_pths{subject_num}, 'processed');
    datasets        =  dir(fullfile(load_pth, '*boot*'));
    before         = load(fullfile(load_pth, datasets(1).name));
    if subject_num < 6
        after         = load(fullfile(load_pth, datasets(3).name));
    else 
         after         = load(fullfile(load_pth, datasets(2).name));
    end
         
    denoisedData = dir(fullfile(load_pth, sprintf('s0%d_denoisedData*',subject_num+1)));
    denoisedData = load(fullfile(load_pth,denoisedData(1).name));
    
    
    model_fit_before = nanmedian(before.fit_f2,4);
    
   
    t = (1:500)/1000;
    f = (0:length(t)-1)/max(t);
    for chan =  1
        figure; 
        for ii = 1:9
        plot(f,model_fit_before(ii,:,chan), 'color', cmap(ii,:,:), 'LineWidth',2); hold on;
        end
         plot(f,model_fit_before(10,:,chan),'k--','LineWidth',2);
        set(gca, 'YScale','log','XScale','log','LineWidth',2)
        xlim([10 200])
        ylim([10.^-.2, 10.^.2])
        xlabel('Frequency (Hz)','FontSize',18)
        ylabel('Power','FontSize',18)
        legend(condition_names{:}, 'FontSize',18)
        box off        set(gcf, 'color','w')

        title('Model fit, subject 6, chan 1, before denoising','FontSize',18)
    end
    
    model_fit_after = nanmedian(after.fit_f2,4);
    
    t = (1:800)/1000;
    f = (0:length(t)-1)/max(t);
    for chan =  1
        figure; 
        for ii = 1:9
        plot(f,model_fit_after(ii,:,chan), 'color', cmap(ii,:,:), 'LineWidth',2); hold on;
        end
         plot(f,model_fit_after(10,:,chan),'k--','LineWidth',2);
        set(gca, 'YScale','log','XScale','log','LineWidth',2)
        xlim([10 150])
         ylim([10.^-.2, 10.^.2])
        xlabel('Frequency (Hz)','FontSize',18)
        ylabel('Power','FontSize',18)
        legend(condition_names{:}, 'FontSize',18)
        title('Model fit, subject 6, chan 1, after denoising','FontSize',18)

    end
    
   
    
end