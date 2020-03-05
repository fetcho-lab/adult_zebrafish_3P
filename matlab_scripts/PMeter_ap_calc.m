function P_real=PMeter_ap_calc(P0,aperture_size)
% this function calculates the effective power of a reduced size aperture
% the parameters are the P0 measured power (which is the power that is
% measured by the refeence PD assuming full aperture), and the aperture
% diameter.
% P0 - measured power (for full diameter)
% aperture_size - Aperture diameter

full_ap=15.12;%
beam_2W=21.5;%
if nargin>1
    P_real=P0./(1-exp(-2*full_ap.^2/beam_2W^2))*(1-exp(-2*aperture_size.^2/beam_2W.^2));
end