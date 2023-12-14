function alginfo = alg_info() %<<<1
% Part of QWTB. Info script for algorithm PJVS-DIFF
%
% See also qwtb

    alginfo.id = 'PJVS-DIFF';
    alginfo.name = 'Amplitude from PJVS differential measurement';
    alginfo.desc = 'Calculates amplitude of main signal component from differentiall measurements using Programmable Josephson Voltage Standard.';
    alginfo.citation = 'Implementation: Martin Sira';
    alginfo.remarks = 'Sampling, PJVS and DUT signals has to be coherent.';
    alginfo.license = 'MIT';
    alginfo.providesGUF = 0;
    alginfo.providesMCM = 0;  

    % --- input quantities %<<<2
    pid = 1;

    % --- common inputs %<<<3
    % sample data
    alginfo.inputs(pid).name = 'fs';
    alginfo.inputs(pid).desc = 'Sampling frequency';
    alginfo.inputs(pid).alternative = 1;
    alginfo.inputs(pid).optional = 0;
    alginfo.inputs(pid).parameter = 0;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'Ts';
    alginfo.inputs(pid).desc = 'Sampling time';
    alginfo.inputs(pid).alternative = 1;
    alginfo.inputs(pid).optional = 0;
    alginfo.inputs(pid).parameter = 0;
    pid = pid + 1;

    alginfo.inputs(pid).name = 't';
    alginfo.inputs(pid).desc = 'Time series';
    alginfo.inputs(pid).alternative = 1;
    alginfo.inputs(pid).optional = 0;
    alginfo.inputs(pid).parameter = 0;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'f';
    alginfo.inputs(pid).desc = 'Signal frequency';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 0;
    alginfo.inputs(pid).parameter = 0;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'y';
    alginfo.inputs(pid).desc = 'Sampled voltage';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 0;
    alginfo.inputs(pid).parameter = 0;
    pid = pid + 1;

    % --- pjvs data inputs %<<<3
    alginfo.inputs(pid).name = 'Uref';
    alginfo.inputs(pid).desc = 'PJVS reference values (V)';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 0;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'fseg';
    alginfo.inputs(pid).desc = 'Frequency of PJVS segments (Hz)';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 0;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'Rs';
    alginfo.inputs(pid).desc = 'PJVS remove before PJVS step (samples)';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 1;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'Re';
    alginfo.inputs(pid).desc = 'PJVS remove after PJVS step (samples)';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 1;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'Ms';
    alginfo.inputs(pid).desc = 'PJVS remove at record beginning (samples)';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 1;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'Me';
    alginfo.inputs(pid).desc = 'PJVS remove at record end (samples)';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 1;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'plots';
    alginfo.inputs(pid).desc = 'Algorithm will plot figures if nonzero. Values are 1 (only basic figures) to 4 (plot all).';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 1;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    alginfo.inputs(pid).name = 'data_folder';
    alginfo.inputs(pid).desc = 'Path where the plots will be saved.';
    alginfo.inputs(pid).alternative = 0;
    alginfo.inputs(pid).optional = 1;
    alginfo.inputs(pid).parameter = 1;
    pid = pid + 1;

    % outputs %<<<2
    pid = 1;
    alginfo.outputs(pid).name = 'U';
    alginfo.outputs(pid).desc = 'RMS voltage of whole signal';
    pid = pid + 1;
    
    alginfo.outputs(pid).name = 'U_t';
    alginfo.outputs(pid).desc = 'RMS voltage per signal period';
    pid = pid + 1;
    
    alginfo.outputs(pid).name = 'U_fft';
    alginfo.outputs(pid).desc = 'Amplitude from FFT calculation';
    pid = pid + 1;
    
    alginfo.outputs(pid).name = 'U_fft_l';
    alginfo.outputs(pid).desc = 'Amplitude of noise 1 frequency bin to the left of the main signal.';
    pid = pid + 1;
    
    alginfo.outputs(pid).name = 'U_fft_r';
    alginfo.outputs(pid).desc = 'Amplitude of noise 1 frequency bin to the right of the main signal.';
    pid = pid + 1;

end % function alginfo = alg_info() %<<<1

% vim settings modeline: vim: foldmarker=%<<<,%>>> fdm=marker fen ft=octave textwidth=80 tabstop=4 shiftwidth=4
