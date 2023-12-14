function dataout = alg_wrapper(datain, calcset)
% Part of QWTB. Wrapper script for algorithm SFDR.
%
% See also qwtb

% Format input data --------------------------- %<<<1
% alternatives
if isfield(datain, 'fs')
    sigconfig.fs = datain.fs.v;
elseif isfield(datain, 'Ts')
    sigconfig.fs = 1/datain.Ts.v;
    if calcset.verbose
        disp('QWTB: SP-WFFT wrapper: sampling frequency was calculated from sampling time')
    end
else
    sigconfig.fs = 1/mean(diff(datain.t.v));
    if calcset.verbose
        disp('QWTB: SP-WFFT wrapper: sampling frequency was calculated from time series')
    end
end

sigconfig.PRs = datain.Rs.v;
sigconfig.PRe = datain.Re.v;
sigconfig.MRs = datain.Ms.v;
sigconfig.MRe = datain.Me.v;
sigconfig.fseg = datain.fseg.v;
sigconfig.f = datain.f.v;

% debug structure:
dbg = check_gen_dbg(datain.plots.v);
dbg.plotpath = datain.data_folder.v;

% XXX ENSURE data_folder exists!

% Call algorithm ---------------------------  %<<<1
% algorithm definition:
% function [A_rms_total, A_rms_t, A_fft, A_fft_l, A_fft_r] = diff_process(sigconfig, y, Uref1period, dbg) %<<<1
%
[dataout.U.v, dataout.U_t.v, dataout.U_fft.v, dataout.U_fft_l.v, dataout.U_fft_r.v] = diff_process(...
    sigconfig, ...
    datain.y.v, ...
    datain.Uref.v, ...
    dbg);

% Format output data:  --------------------------- %<<<1

end % function

% vim settings modeline: vim: foldmarker=%<<<,%>>> fdm=marker fen ft=octave textwidth=80 tabstop=4 shiftwidth=4
