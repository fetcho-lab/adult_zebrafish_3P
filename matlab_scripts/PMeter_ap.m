function signal=PMeter_ap(signal,cal_a,cal_b,aperature_size)
full_ap=15.12;
beam_2W=21.5;
signal = signal*cal_a + cal_b;
signal = signal./(1-exp(-2*full_ap.^2/beam_2W^2)) * (1 - exp(-2 * aperature_size.^2/beam_2W.^2));