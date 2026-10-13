function mucalib_depthcalc();

% Function calculating depth profile muon production parameters without sigma0 and fstar.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2017-2026 (jakob.heyman@gu.se)

clear all;
tic();

% run and load expage constants
make_consts_expage;
load consts_expage;

% Calibration sites for production from muons
site = {'Beacon','LeymonHigh','LeymonLow','LaCiotat','LeymonHigh14'};

for i = 1:numel(site);
    % fix nucl10, nucl26, nucl14
    if strcmp(site{i},'LeymonHigh14');
        nucl10 = 0; nucl26 = 0; nucl14 = 1;
    else;
        nucl10 = 1; nucl26 = 1; nucl14 = 0;
    end;

    % display site
    fprintf(1,'Site: %s\n',site{i});

    % load site data from mucalib_data
    mucal = mucalib_data(site{i});

    % calculate atmospheric pressure
    atm = ERA40atm(mucal.lat,mucal.lon,mucal.elv);

    % Fix RcEst for muon calculation
    LSDfix = LSD_fix(mucal.lat,mucal.lon,1E7,mucal.samplingyr,consts);

    % muon production (without sigma0 and fstar!)
    fprintf(1,'calculating muon P (up to %.0f)\n',mucal.dz(end));
    Pmud = mucalib_Pmu(mucal.dz,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,nucl14,consts);
    fprintf(1,'\ndone!\n');

    % save parameters in mucalpre.site.parameter
    mucalpre.(site{i}).dz = mucal.dz';
    if nucl10 == 1;
        mucalpre.(site{i}).Pfast10d = Pmud.P_fast10';
        mucalpre.(site{i}).Pneg10d = Pmud.P_neg10';
    end;
    if nucl26 == 1;
        mucalpre.(site{i}).Pfast26d = Pmud.P_fast26';
        mucalpre.(site{i}).Pneg26d = Pmud.P_neg26';
    end;
    if nucl14 == 1;
        mucalpre.(site{i}).Pfast14d = Pmud.P_fast14';
        mucalpre.(site{i}).Pneg14d = Pmud.P_neg14';
    end;
end;

% save mucalpre as mucalib_depthcalc_precalc.mat
fprintf(1,'saving mucalib_depthcalc_precalc.mat');
save('-v7','-zip','mucalib_depthcalc_precalc.mat','mucalpre');
fprintf(1,' - done!\n');

toc();
