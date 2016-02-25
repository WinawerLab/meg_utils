%% Visualize Neuromag306 data
%
% http://www.megwiki.org/images/f/fd/Elekta_Neuromag_User's_Manual_-_System_Hardware.pdf
%
% The 306 sensors comprise 102 magnetometers that measure the Bz
% component which is perpendicular to the surface of the detector of the
% field, 102 planar gradiometers measuring the gradient ?Bz/?x, and 102
% planar gradiometers measuring ?Bz/?y-component of the gradient. The
% sensors are arranged in triple sensor elements each comprising two
% orthogonal planar gradiometers and one magnetomer in the same plane
% as the planar gradiometers.
%
% Data can be viewed in several ways: 
%   - 204 gradiometers ('neuromag306planar') 
%   - 102 planar gradiometers by combining the x- and y (neuromag306cmb)
%   - 306 all channels ('neuromag306all')
%   - 102 magnetometers ('neuromag306mag')

layouts = {'neuromag306cmb' 'neuromag306planar' 'neuromag306all' 'neuromag306mag'};
n = [102 204 306 102];

figure
for l = 1:length(layouts)
    
    subplot(2,2,l)
    cfg = [];
    cfg.marker    = 'on';
    cfg.layout    = layouts{l};
    cfg.electrodes = 'numbers';
    cfg.style='straight';
    
    topoplot(cfg, (1:n(l))'); 
    title(strsplit(layouts{l}, 'neuromag306'))
end

