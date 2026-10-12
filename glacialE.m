function glacialE();

% 10Be/26Al/14C glacial erosion calculator.
% Function for quantification/investigation of glacial erosion based on 10Be/26Al/14C conc.
% Interpolates the best fitting glacial erosion to yield the 10Be/26Al/14C conc given as input
% based on the time-dependent cosmogenic nuclide production rate.
%
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman (jakob.heyman@gu.se) 2018-2026

clear; close all; tic();

% What version is this?
ver = '202610';

% =============== MAKE CHOICES HERE ================================================================
% max time (determines the start of the simulation - max 1E7)
mt = 1E7;

% calculate glacial erosion rate and/or incremental erosion depth steps? 1 = yes, 0 = no
glacErate = 1;
glacEstep = 1;

% Quasi-Monte Carlo sample size for uncertainty (rank-1 lattice rule: qmc should be a power of 2)
qmc = 2^10;

% calculate combined 10Be-26Al erosion? 1 = yes, 0 = no
nucl1026 = 0;

% plotting choices
plch.depth = 1;         % plot sample depth
plch.P = 0;             % plot sample production
plch.N = 0;             % plot nuclide build-up
plch.Pprcnt = 0;        % plot sample production as percent of P at time 0
plch.Nprcnt = 0;        % plot sample conc as percent of N at time 0
plch.combined = 0;      % plot 10Be and 26Al data together? (requires combined_full ~= 1)
plch.combined_full = 0; % plot 10Be, 26Al and potential combined 1026 data together?
plch.uncline = 0;       % plot uncertainties as lines
plch.uncfill = 0;       % plot uncertainties as filled areas
plch.probdensmap = 0;   % plot summed probability density map (pdm) of Latin Hypercube output
  plch.pdmX = 201;      % number of x points in probability density map
  plch.pdmY = 151;      % number of x points in probability density map
  plch.pdm_save = 1;    % save probability density map as txt file
plch.pct_area = 1;      % plot summed percentiles
  plch.pct_dt = 1000;    % time step (yr) for percentile plot
  plch.pct_areaN = 20; % number of percentile groups
  plch.colorbar = 0;    % plot percentile area colorbar
plch.pct_lines = [25 50 75];    % vector for plotting summed percentile lines
plch.pct_unc = 0;       % calculate percentiles based on uncertainty data
plch.banana = 0;        % plot 26/10 banana - only possible when nucl1026 = 1
plch.maxt = 120000;     % max time for plotting (yr)
plch.maxt = plch.maxt.*1E-3; % change to ka
plch.maxd = 15;         % max depth for plotting (m)
plch.cutuncmaxd = 1;    % cut uncertainties along maxd (for vector handling in other software)
plch.cutmaxxy = 1;      % cut lines at plot edges
plch.splitdd = 5;       % split sample depth lines for samples with max erosion at depth in m
plch.clr10 = [0 0 0];   % color for 10Be: red
plch.clr26 = [0 0 1];   % color for 26Al: blue
plch.clr1026 = [0 0 0]; % color for combined 10Be+26Al: black
plch.clr14 = [0 1 0];   % color for 14C: green

% estimate sensitivity of individual parameters? 1 = yes, 0 = no
uncsens.yes = 0; % if 0: no sensitivity calculations
uncsens.N = 1;
uncsens.Pref = 1;
uncsens.l = 0;
uncsens.dens = 1;
uncsens.erosion = 1;
uncsens.glacE = 1;
uncsens.simt = 1;
uncsens.icevalue = 1;
uncsens.deglac = 1;
uncsens.isostsubm = 1;
uncsens.full = 1; % display full uncertainty in plot?
% uncertainty sensitivity colors for plotting
sensclr.N = [1 0.1 0.1];
sensclr.Pref = [0.6 0.6 0.6];
sensclr.l = [0.6 0.6 1];
sensclr.dens = [1 0.6 0];
sensclr.erosion = [0 0 0];
sensclr.glacE = [0 0 0];
sensclr.simt = [0.5 0.5 0];
sensclr.icevalue = [0 0 0.8];
sensclr.deglac = [1 0 1];
sensclr.isostsubm = [0 1 0];
sensclr.full = [0.8 0.8 0.8];

% time vector (yr) for writing sample depth (m) in output
tdv = [120000,2588000];
% tdv = [];
% time vector (yr) for writing glacial erosion depth (m) in output
gdv = [120000,2588000];
gdv = [];

% ice cover proxy file - default is Lisiecki and Raymo (2005) d18O ice volume record
% column 1: year   column 2: proxy value
iceproxy = 'LR04.txt';
iceproxy_yr0 = 1950; % start year (= 0) in iceproxy file

% glacial (or non-glacial) erosion vector
Etestv = [(0:0.01:0.99) logspace(0,4,900)]; % mm/ka or cm/glac
% ==================================================================================================

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','isostsubm','isostP','glacErv','glacErv_tv',...
    'glacEiv','glacEiv_tv','erosionv','erosionv_tv','burialdepthv','burialdepthv_tv','shieldv',...
    'shieldv_tv','icevaluev','icevaluev_tv','ice_tv','noice_tv','lat','long','elv','thick',...
    'dens','densunc','shield','erosion','erosionunc','N10','N10unc','N26','N26unc','N14',...
    'N14unc','samplingyr','pressure','deglac','deglacunc','deglac_yr0','icevalue','icevalueunc',...
    'glacE','glacEunc','isostsubmunc','burialdepth','simt','simtunc','simt_yr0','densmin','densmax',...
    'erosionmin','erosionmax','deglacmin','deglacmax','icevaluemin','icevaluemax','glacEmin',...
    'glacEmax','isostsubmmin','isostsubmmax','simtmin','simtmax','isostPmod','isostsubmmod',...
    'glacErv_tv0','glacEiv_tv0','erosionv_tv0','burialdepthv_tv0','shieldv_tv0','icevaluev_tv0',...
    'ice_tv0','noice_tv0'};
vartypes = {'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s',...
    '%s','%s','%s','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n','%n'};
% optional variables for pars
opt_pars = {'burialdepth','glacErv','glacErv_tv','glacErv_tv0','glacEiv','glacEiv_tv',...
    'glacEiv_tv0','erosionv','erosionv_tv','erosionv_tv0','burialdepthv','burialdepthv_tv',...
    'burialdepthv_tv0','shieldv','shieldv_tv','shieldv_tv0','icevaluev','icevaluev_tv',...
    'icevaluev_tv0','ice_tv','ice_tv0','noice_tv','noice_tv0'};
% read input file
fid = fopen('input.txt');
varsin = strsplit(fgetl(fid)); % read first line
if(ismember(varsin,varnames)); % if first line contain only variable names
    [testi,vari] = ismember(varsin,varnames); % find index of varnames
    typestr = vartypes{vari(1)}; % fix type string
    for i = 2:numel(vari); % fix type string
        typestr = [typestr ' ' vartypes{vari(i)}];
    end;
else; % if no variable names in first line
    fprintf(1,'ERROR! Something is wrong with the input\n');
    return;
end;
indata = textscan(fid,typestr,'CommentStyle','%','TreatAsEmpty','-'); % scan data
for i = 1:numel(varsin); % fix variables
    samplein.(varsin{i}) = indata{i};
end;
fclose(fid);
% ==================================================================================================

% fix for sensitivity calculations
if uncsens.yes == 0;
    sensfields = fieldnames(uncsens);
    for i = 1:numel(sensfields); uncsens.(sensfields{i}) = 0; end;
end;

% run and load expage constants
make_consts_expage;
load consts_expage;

% number of input rows
nn = numel(samplein.sample);

% fix Prefs, decay constants, and fill empty N with 0
n = {'10','26','14'};
for i = 1:numel(n);
    if isfield(samplein,['N',n{i}]);
        % fix Pref and Prefunc
        if isfield(samplein,['Pref',n{i}]) == 0;
            samplein.(['Pref',n{i}])(1:nn,1) = NaN(nn,1);
        end;
        if isfield(samplein,['Pref',n{i},'unc']) == 0;
            samplein.(['Pref',n{i},'unc'])(1:nn,1) = NaN(nn,1);
        end;
        if isfield(samplein,'isostP'); % if using isostatic P
            isoPv = ((strcmp(samplein.isostP,'-') == 0) & isnan(samplein.(['Pref',n{i}])));
            isoPuncv = ((strcmp(samplein.isostP,'-') == 0) & isnan(samplein.(['Pref',n{i},'unc'])));
            samplein.(['Pref',n{i}])(isoPv) = consts.(['Pref',n{i},'iso']);
            samplein.(['Pref',n{i},'unc'])(isoPuncv) = consts.(['Pref',n{i},'isounc']);
        end;
        samplein.(['Pref',n{i}])(isnan(samplein.(['Pref',n{i}]))) = consts.(['Pref',n{i}]);
        samplein.(['Pref',n{i},'unc'])(isnan(samplein.(['Pref',n{i},'unc']))) = ...
            consts.(['Pref',n{i},'unc']);
    else;
        % if no N: fill with 0
        samplein.(['N',n{i}])(1:nn,1) = 0;
    end;
    % if no Nunc: fill with 0
    if isfield(samplein,['N',n{i},'unc']) == 0; samplein.(['N',n{i},'unc'])(1:nn,1) = 0; end;
    % if there is NaN in N10/N26/N14: replace with 0
    samplein.(['N',n{i}])(isnan(samplein.(['N',n{i}]))) = 0;
    samplein.(['N',n{i},'unc'])(isnan(samplein.(['N',n{i},'unc']))) = 0;
end;

% convert 10Be concentrations according to standards
if isfield(samplein,'std10') == 0; samplein.std10(1:nn,1) = {'0'}; end;
[testi,stdi] = ismember(samplein.std10,consts.std10); % find index of standard conversion factors
if any(testi == 0); % unknown standard name -> friendly error instead of hitting std10_cf(0)
    badstd = unique(samplein.std10(testi == 0));
    errmsg = sprintf('%s, ', badstd{:});
    error(['glacialE: unknown 10Be standard name(s): ' errmsg(1:end-2) ...
        ' (accepted standards are defined in make_consts_expage.m)']);
end;
mult10 = consts.std10_cf(stdi); % pick out conversion factor
samplein.N10 = samplein.N10 .* mult10;
samplein.N10unc = samplein.N10unc .* mult10;

% convert 26Al concentrations according to standards
if isfield(samplein,'std26') == 0; samplein.std26(1:numel(samplein.sample),1) = {'0'}; end;
[testi,stdi] = ismember(samplein.std26,consts.std26); % find index of standard conversion factors
if any(testi == 0); % unknown standard name -> friendly error instead of hitting std26_cf(0)
    badstd = unique(samplein.std26(testi == 0));
    errmsg = sprintf('%s, ', badstd{:});
    error(['glacialE: unknown 26Al standard name(s): ' errmsg(1:end-2) ...
        ' (accepted standards are defined in make_consts_expage.m)']);
end;
mult26 = consts.std26_cf(stdi); % pick out conversion factor
samplein.N26 = samplein.N26 .* mult26;
samplein.N26unc = samplein.N26unc .* mult26;

% fix longitude values
samplein.long(find(samplein.long < 0)) = samplein.long(find(samplein.long < 0)) + 360;

% fix sample pressure
if isfield(samplein,'pressure') == 0;
    % if there is no pressure flag: use std
    if isfield(samplein,'Pflag') == 0; samplein.Pflag(1:nn,1) = {'std'}; end;
    stdv = strcmp(samplein.Pflag,'std');
    antv = strcmp(samplein.Pflag,'ant');
    prev = strcmp(samplein.Pflag,'pre');
    samplein.pressure(stdv) = ERA40atm(samplein.lat(stdv),samplein.long(stdv),samplein.elv(stdv));
    samplein.pressure(antv) = antatm(samplein.elv(antv));
    samplein.pressure(prev) = samplein.elv(prev);
end;

% read iceproxy file
fid = fopen(iceproxy);
inproxy = textscan(fid,'%n %n','CommentStyle','%');
iceproxy_tv = inproxy{1};
iceproxyin = inproxy{2};
fclose(fid);

% fix isostatic adjustment input
if isfield(samplein,'isostP');
    samplein.isostP = isost_data_load(samplein.isostP);
end;
if isfield(samplein,'isostsubm');
    samplein.isostsubm = isost_data_load(samplein.isostsubm);
end;

% Muon shielding depth vector to ~200,000 g/cm2
z_mu = [(0:3:12) logspace(1.18,5.3,50)];

% glacial erosion matrices
glacErm = [];
glacEim = [];

% number of elements in Etestv (for precalc)
Etestvn = numel(Etestv);

% check if input contains both erosion and glacE
if isfield(samplein,'erosion') && isfield(samplein,'glacE');
    fprintf(1,'Input contains both glacial and nonglacial erosion - ');
    fprintf(1,'there is nothing to calculate!\n');
    return;
end;

% check if erosion calculation should be done for glacial erosion or non-glacial erosion
if isfield(samplein,'glacE');
    glaccalc = 0; % calculate nonglacial erosion
else;
    glaccalc = 1; % calculate glacial erosion
    % if there is no erosion in input: assume zero erosion
    if isfield(samplein,'erosion') == 0 && ...
            (isfield(samplein,'erosionmin') == 0 || isfield(samplein,'erosionmax') == 0);
        samplein.erosion(1:nn,1) = 0;
    end;
end;

% fix for output and plotting
output(1,1) = {'sample'};
outcol.sample = 1;
[output,outcol,sensstr] = fixoutput1(output,outcol,samplein,glacErate,glacEstep,glaccalc,uncsens,...
    tdv,gdv,nucl1026);

% for plotting
pl = plfix(plch,samplein,nucl1026,glacErate,glacEstep,sensstr);

% fix l vector for uncertainty: deterministic rank-1 lattice low-discrepancy sampling
l10v = norminv(latt_u(qmc,3),consts.l10,consts.l10unc)';
l26v = norminv(latt_u(qmc,3),consts.l26,consts.l26unc)';
l14v = norminv(latt_u(qmc,3),consts.l14,consts.l14unc)';

% parameters in input
sampleinfields = fieldnames(samplein);

% pick out samples one by one
for i = 1:numel(samplein.lat);
    % pick out sample data
    for j = 1:numel(sampleinfields);
        sample.(sampleinfields{j}) = samplein.(sampleinfields{j})(i);
        % remove fields with input '-' or numeric fields with NaN input
        % (guard ischar/iscell so strcmp never sees a struct/numeric in MATLAB)
        sval = sample.(sampleinfields{j});
        if ((ischar(sval) || (iscell(sval) && numel(sval) == 1 && ischar(sval{1}))) && strcmp(sval,'-')) || ...
                (isnumeric(sval) && isnan(sval));
            sample = rmfield(sample,sampleinfields{j});
        end;
    end;

    % write sample name to output
    output(i+1,1) = sample.sample;

    % Set nucl to 0 for 10/26/14 and check if there is N10/N26/N14
    nucl10 = 0; nucl26 = 0; nucl14 = 0;
    if sample.N10 > 0; nucl10 = 1; end;
    if sample.N26 > 0; nucl26 = 1; end;
    if sample.N14 > 0; nucl14 = 1; end;

    % if no 10Be/26Al/14C: move on
    if nucl10+nucl26+nucl14 == 0;
        clear sample;
        continue;
    end;

    % remove potential isostatic adjustment if having atmospheric pressure as input
    if isfield(sample,'isostP') && strcmp(sample.Pflag,'pre');
        sample = rmfield(sample,'isostP');
    end;
    
    % fix for min/max parameters without central value
    minmaxp = {'dens','erosion','deglac','icevalue','glacE','isostsubm','simt'};
    for j = 1:numel(minmaxp);
        if isfield(sample,[minmaxp{j} 'min']) && isfield(sample,[minmaxp{j} 'max']) && ...
            isfield(sample,minmaxp{j})==0;
            sample.(minmaxp{j}) = mean([sample.([minmaxp{j} 'min']) sample.([minmaxp{j} 'max'])]);
        end;
    end;

    % fix for using simulation zero year
    if isfield(sample,'simt_yr0');
        sample.samplingyr_in = sample.samplingyr;
        sample.samplingyr = sample.simt_yr0;
    end;

    % fix for deglac_yr0
    if isfield(sample,'deglac_yr0') && isfield(sample,'deglac');
        sample.deglac = sample.deglac + sample.samplingyr - sample.deglac_yr0;
        if isfield(sample,'deglacmin') && isfield(sample,'deglacmax');
            sample.deglacmin = sample.deglacmin + sample.samplingyr - sample.deglac_yr0;
            sample.deglacmax = sample.deglacmax + sample.samplingyr - sample.deglac_yr0;
        end;
    end;

    % find out if calculating glacial or non-glacial erosion
    if isfield(sample,'glacE') == 0;
        sample.glacE = Etestv;
    else;
        sample.erosion = Etestv;
    end;
    
    % Prefs
    if isfield(sample,'Pref10');
        Pref10 = sample.Pref10; Pref10unc = sample.Pref10unc;
        Pref10v = norminv(latt_u(qmc,1),Pref10,Pref10unc)';
    end;
    if isfield(sample,'Pref26');
        Pref26 = sample.Pref26; Pref26unc = sample.Pref26unc;
        Pref26v = norminv(latt_u(qmc,1),Pref26,Pref26unc)';
    end;
    if isfield(sample,'Pref14');
        Pref14 = sample.Pref14; Pref14unc = sample.Pref14unc;
        Pref14v = norminv(latt_u(qmc,1),Pref14,Pref14unc)';
    end;

    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample{1});
    
    % set simulation max time
    dimn = 4; % dimension number for lattice points generation
    if isfield(sample,'simtmin') && isfield(sample,'simtmax');
        sample.simtqmc = sample.simtmin + (sample.simtmax-sample.simtmin)*latt_u(qmc,dimn)';
        sample.mt = min([max(sample.simtqmc) 1E7]);
        dimn = 5;
    end;
    if isfield(sample,'simt') && isfield(sample,'simtunc');
        if isfield(sample,'simtqmc');
            fprintf(1,'\nINPUT ERROR! This sample has both min/max and +/- uncertainty for parameter simt');
            dimn = 4;
        end;
        sample.simtqmc = norminv(latt_u(qmc,dimn),sample.simt,sample.simtunc)';
        sample.mt = min([max(sample.simtqmc) 1E7]);
        dimn = 5;
    end;
    if isfield(sample,'simt') && isfield(sample,'mt')==0;
        sample.mt = sample.simt;
    end;
    if isfield(sample,'mt') == 0;
        sample.mt = mt;
    end;
    
    % Age Relative to t0=2010 - LSD tv from LSD_fix
    % tv = [0:10:50 60:100:3060 3260:200:75060 76060:1000:800060 802060:2000:2000060 logspace(log10(2002060),7,200)];

    % Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
    LSDfix = LSD_fix(sample.lat,sample.long,sample.mt,sample.samplingyr,consts);

    % add time points for ice cover start/end and potential submergence
    [sample,LSDfix,pars] = add_tv_points(sample,LSDfix,iceproxy_tv,iceproxyin,iceproxy_yr0);

    % Production from muons along shielding depth vector z_mu
    Pmu_z = P_mu_expage(z_mu,sample.pressure,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,nucl14,...
        consts,'no');
    
    % fix atmospheric pressure if using isostatic adjustment
    if isfield(sample,'isostP');
        % calculate elevation development
        sample.elvv = isost_elv(sample.isostP{1},sample);
        if isfield(sample,'isostPmod');
            sample.elvv = (sample.elvv - sample.elv) .* sample.isostPmod + sample.elv;
        end;
        % calculate atmospheric pressure
        if strcmp(sample.Pflag,'std');
            sample.pressure = ERA40atm(sample.lat,sample.long,sample.elvv);
        elseif strcmp(sample.Pflag,'ant');
            sample.pressure = antatm(sample.elvv);
        end;
    end;
    
    % spallation production scaling
    Psp = P_sp_expage(sample.pressure,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,nucl10,nucl26,nucl14);
    
    % interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
    Lsp = rawattenuationlength(sample.pressure,LSDfix.Rc);
    
    % fix parameters for calculation preparations
    pars.thick = sample.thick;
    pars.Lsp = Lsp(:);
    pars.samplingyr = sample.samplingyr;
    pars.tv = sample.tv;
    pars.z_mu = z_mu;
    % fix optional variables for pars
    samplefields = fieldnames(sample);
    optvars = samplefields(ismember(samplefields,opt_pars));
    for j = 1:numel(optvars);
        if isnumeric(sample.(optvars{j}));
            pars.(optvars{j}) = sample.(optvars{j});
        else;
            pars.(optvars{j}) = str2num(sample.(optvars{j}){1});
        end;
    end;
    % fix shielding vector if using time-varying shielding
    if isfield(pars,'shieldv') && isfield(pars,'shieldv_tv');
        sample.shield(1:numel(sample.tv),1) = sample.shield;
        minmaxv = fix_minmaxv(pars,'shieldv');
        for j = 1:size(minmaxv,1);
            sample.shield(minmaxv(j,1):minmaxv(j,2),1) = pars.shieldv(j);
        end;
    end;
    pars.shield = sample.shield;
    
    % fix variable parameters
    pars.glacE = sample.glacE;
    pars.icevalue = sample.icevalue;
    pars.erosion = sample.erosion;
    pars.dens = sample.dens;
    if isfield(sample,'deglac'); pars.deglac = sample.deglac; end;
    if isfield(sample,'delv'); pars.delv = sample.delv; pars.elv = sample.elv; end;
    if isfield(sample,'simt'); pars.simt = sample.simt; else; pars.simt = sample.mt; end;
    
    % calculate ice cover, erosion, and depth matrices
    prp = precalc(pars,Etestvn);
    if isfield(prp,'waterdm'); pars.waterdm = prp.waterdm; end;
    prp2 = getdepth(prp,pars,glacErate,glacEstep,Etestvn);

    % nuclide-specific erosion calculations
    if nucl10 == 1; % 10Be
        % fix for using simulation zero year
        if isfield(sample,'simt_yr0') && prp.icev(1)==0;
            sample.N10 = adjustN(sample,pars,sample.N10,Pref10,Psp.sp10,Pmu_z.mu10,consts.l10);
        end;
        % fix nuclide specific parameters
        np10.Psp = Psp.sp10;
        np10.Pmu_z = Pmu_z.mu10;
        np10.N = sample.N10; np10.Nunc = sample.N10unc;
        np10.Nv = norminv(latt_u(qmc,2),sample.N10,sample.N10unc)';
        np10.Nv(np10.Nv<0) = 0; % don't allow negative N
        np10.Pref = Pref10; np10.Prefunc = Pref10unc;
        np10.l = consts.l10;
        % calculate erosion and internal unc
        N10E = nuclE(pars,prp2,np10,Etestv,glaccalc);
    end;
    if nucl26 == 1; % 26Al
        % fix for using simulation zero year
        if isfield(sample,'simt_yr0') && prp.icev(1)==0;
            sample.N26 = adjustN(sample,pars,sample.N26,Pref26,Psp.sp26,Pmu_z.mu26,consts.l26);
        end;
        % fix nuclide specific parameters
        np26.Psp = Psp.sp26;
        np26.Pmu_z = Pmu_z.mu26;
        np26.N = sample.N26; np26.Nunc = sample.N26unc;
        np26.Nv = norminv(latt_u(qmc,2),sample.N26,sample.N26unc)';
        np26.Nv(np26.Nv<0) = 0; % don't allow negative N
        np26.Pref = Pref26; np26.Prefunc = Pref26unc;
        np26.l = consts.l26;
        % calculate erosion and internal unc
        N26E = nuclE(pars,prp2,np26,Etestv,glaccalc);
    end;
    if nucl14 == 1; % 14C
        % fix for using simulation zero year
        if isfield(sample,'simt_yr0') && prp.icev(1)==0;
            sample.N14 = adjustN(sample,pars,sample.N14,Pref14,Psp.sp14,Pmu_z.mu14,consts.l14);
        end;
        % fix nuclide specific parameters
        np14.Psp = Psp.sp14;
        np14.Pmu_z = Pmu_z.mu14;
        np14.N = sample.N14; np14.Nunc = sample.N14unc;
        np14.Nv = norminv(latt_u(qmc,2),sample.N14,sample.N14unc)';
        np14.Nv(np14.Nv<0) = 0; % don't allow negative N
        np14.Pref = Pref14; np14.Prefunc = Pref14unc;
        np14.l = consts.l14;
        % calculate erosion and internal unc
        N14E = nuclE(pars,prp2,np14,Etestv,glaccalc);
    end;
    
    % generate uncertainty parameters for qmc iterations
    parsu = uncqmc(pars,sample,qmc,dimn);
    if isfield(parsu,'waterdm'); parsu = rmfield(parsu,'waterdm'); end;
    
    % do precalc for uncertainty estimation
    prpunc = precalc(parsu,qmc);
    if isfield(prpunc,'waterdm'); parsu.waterdm = prpunc.waterdm; end;
    
    % fix for waterdm
    if isfield(prp,'waterdm'); pars.waterdm = prp.waterdm(:,1); end;
    
    % do precalc for uncertainty sensitivity and/or 10Be+26Al calc
    prpuncsens = struct;
    if numel(sensstr)>0 || nucl10*nucl26*nucl1026==1;
        prpuncsens = precalc(pars,qmc);
    end;
    
    % do precalc for depthPconc/tdv/1026
    prp1 = precalc(pars,1);
    
    % nuclide-specific uncertainty calculations
    if nucl10 == 1; % 10Be
        [output,pl,unc10] = uncfix(output,pl,outcol,plch,uncsens,sensstr,sample,pars,parsu,prp1,...
            prpunc,prpuncsens,Etestv,tdv,gdv,N10E,np10,Pref10v,l10v,glaccalc,i,qmc,10);
    end;
    if nucl26 == 1; % 26Al
        [output,pl,unc26] = uncfix(output,pl,outcol,plch,uncsens,sensstr,sample,pars,parsu,prp1,...
            prpunc,prpuncsens,Etestv,tdv,gdv,N26E,np26,Pref26v,l26v,glaccalc,i,qmc,26);
    end;
    if nucl14 == 1; % 14C
        [output,pl,unc14] = uncfix(output,pl,outcol,plch,uncsens,sensstr,sample,pars,parsu,prp1,...
            prpunc,prpuncsens,Etestv,tdv,gdv,N14E,np14,Pref14v,l14v,glaccalc,i,qmc,14);
    end;
    % 10Be + 26Al
    if nucl10==1 && nucl26==1 && nucl1026==1;
        % calculate depth P for erosion line and mean P for sample data
        if plch.banana == 1 || (isfield(N10E,'Nendr') && isfield(N26E,'Nendr')) || ...
                (isfield(N10E,'Nendi') && isfield(N26E,'Nendi'));
            [pl,P10u,P26u] = P1026fix(pl,sample,LSDfix,Pmu_z,z_mu,np10,np26);
        end;
        % define E1026 as structure
        E1026 = struct;
        % calculate 1026 erosion
        if isfield(N10E,'Nendr') && isfield(N26E,'Nendr');
            [output,pl,E1026] = Ecalc1026(output,pl,E1026,outcol,plch,sample,pars,parsu,prp1,...
                prpunc,prpuncsens,Etestv,N10E,N26E,unc10,unc26,np10,np26,Pref10v,Pref26v,l10v,...
                l26v,P10u,P26u,tdv,gdv,glaccalc,qmc,i,'r');
        end;
        if isfield(N10E,'Nendi') && isfield(N26E,'Nendi');
            [output,pl,E1026] = Ecalc1026(output,pl,E1026,outcol,plch,sample,pars,parsu,prp1,...
                prpunc,prpuncsens,Etestv,N10E,N26E,unc10,unc26,np10,np26,Pref10v,Pref26v,l10v,...
                l26v,P10u,P26u,tdv,gdv,glaccalc,qmc,i,'i');
        end;
        % fill pl.banana with N
        if plch.banana == 1;
            pl = bananaN(pl,E1026,np10,np26,P10u,P26u);
        end;
    end;
    
    % fix plot structure
    if plch.depth==1 || plch.P==1 || plch.N==1;
        pl.tm(1:numel(sample.tv),end+1) = sample.tv.*1E-3;
        pl.tmunc(1:numel(sample.tv)*2,end+1) = [flip(sample.tv.*1E-3);sample.tv.*1E-3];
        pl.mm(1:numel(sample.tv),end+1) = 1;
        pl.mmunc(1:numel(sample.tv)*2,end+1) = 1;
    end;
    
    % new line
    fprintf(1,'\n');
    
    % clear parameters
    clear sample; clear pars; clear parsu;
    clear prp; clear prp2; clear prpunc; clear prp1; clear prpuncsens;
    clear np10; clear np26; clear np14;
    clear N10E; clear N26E; clear N14E; clear E1026;
    clear unc10; clear unc26; clear unc14;
end;

% fix and save output ======================================
if sum(samplein.N10 + samplein.N26 + samplein.N14)>0;
    % fix output string
    outstr = '%s';
    for j = 1:size(output,2)-1;
        outstr = strcat(outstr,'\t%s');
    end;
    outstr = strcat(outstr,'\n');

    % fill empty cells with '-'
    nullidx = cellfun(@isempty,output);
    output(nullidx) = {'-'};

    % write out-glacialE.txt
    out = fopen('out-glacialE.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    fclose(out);
end;
% ==========================================================

% plotting ======================================================================================
% plot dPN
dPNstr = {'d','P','N'};
nuclstr = {'10','26','1026','14'};
ristr = {'r','i'};
for a = 1:numel(dPNstr);
    for b = 1:numel(nuclstr);
        for c = 1:numel(ristr);
            pl = plot_dPN(pl,plch,dPNstr{a},nuclstr{b},ristr{c});
        end;
    end;
end;
% uncertainty sensitivity plot
nuclstr = {'10','10','26','26','14','14'};
ristr = {'r','i','r','i','r','i'};
nuclv = [10 10 26 26 14 14];
for j = 1:6;
    if isfield(pl,['sens' nuclstr{j} ristr{j}]);
        plot_uncsens(pl.(['sens' nuclstr{j} ristr{j}]),sensstr,sensclr,glaccalc,ristr{j},nuclv(j));
    end;
end;
% banana plot
if isfield(pl,'banana');
    plot_banana(pl.banana,consts,z_mu);
end;
% ===============================================================================================

toc()
clear;
% end glacialE function ============================================================================


% subfunction precalc ==============================================================================
function out = precalc(pars,numcols);
% calculate ice cover, erosion, and depth matrices
    % fix ice cover record
    icem = (pars.iceproxy >= pars.icevalue);
    icem = double(icem); % fix for matlab

    % fix deglaciation based on sample.deglac
    if isfield(pars,'deglac');
        icem = fix_deglac(pars.deglac,icem,pars.tv);
    end;

    % fix specified ice-cover/noice-cover periods
    if isfield(pars,'ice_tv') || isfield(pars,'noice_tv');
        icem = fixicem(pars,icem);
    end;
    
    % fix icetest matrix
    out.icetest = icem;
    
    % fix water depth matrix - water is assumed to have a density of 1 g/cm3
    if isfield(pars,'delv');
        if isfield(pars,'waterdm');
            out.waterdm = pars.waterdm;
        else;
            out.waterdm = get_waterdepth(pars.elv,pars.delv,icem,pars.tv);
        end;
        if size(out.waterdm,2) == 1;
            out.waterdm = repmat(out.waterdm,1,numcols);
        end;
    else;
        out.waterdm = 0;
    end;
    
    % save icev for plotting and fix ice matrix
    if size(icem,2) == 1;
        out.icev = icem;
        icem = repmat(icem(:),1,numcols);
    end;
    
    % save icem and noicem to output
    out.icem = icem;
    out.noicem = (icem == 0);
% end subfunction precalc ==========================================================================


% subfunction adjustN ==============================================================================
function out = adjustN(sample,pars,N,Pref,Psp,Pmu_z,l);
    % adjust sample N to simt_yr0
    nyr = sample.samplingyr_in - sample.simt_yr0; % number of years to adjust for
    dcf = exp(-nyr.*l); % decay factor
    dpfs = exp(-nyr.*pars.erosion(1).*1E-4.*sample.dens./pars.Lsp(1)); % erosion factor
    if sample.thick > 0;
        thickSF = (pars.Lsp(1)./(sample.dens.*sample.thick)).*...
            (1 - exp(((-1.*sample.dens.*sample.thick)./pars.Lsp(1)))); % thickness scaling factor
    else; thickSF = 1; end;
    Psp_adj = Psp(1).*Pref.*thickSF.*sample.shield(1); % spallation P
    Pmu_adj = interp1(pars.z_mu,Pmu_z,(sample.thick+nyr.*pars.erosion(1).*1E-4).*sample.dens./2,...
        'pchip') .* sample.shield(1); % muon P
    out = (N - nyr.*(Psp_adj.*dpfs.*dcf + Pmu_adj.*dcf)) ./ dcf; % adjusted N
% end subfunction adjustN ==========================================================================


% subfunction getdepth =============================================================================
function out = getdepth(prp,pars,glacErate,glacEstep,numcols);
    % save icem and noicem in output
    out.icem = prp.icem;
    out.noicem = prp.noicem;
    
    % fix non-glacial erosion matrix
    nonglacEm = repmat(pars.erosion,numel(pars.tv),numcols./numel(pars.erosion));
    
    % fix glacial erosion matrices
    if glacErate == 1;
        glacEm.rr = repmat(pars.glacE,numel(pars.tv),numcols./numel(pars.glacE));
        out.icetestr = prp.icetest;
    end;
    if glacEstep == 1;
        glacEm.ii = repmat(pars.glacE,numel(pars.tv),numcols./numel(pars.glacE));
        out.icetesti = prp.icetest;
    end;
    
    % fix non-glacial and glacial erosion in specific time periods
    if isfield(pars,'glacErv') || isfield(pars,'glacEiv') || isfield(pars,'erosionv');
        [nonglacEm,glacEm,out] = fixE(pars,nonglacEm,glacEm,out);
    end;
    
    % fix dens matrix and save to output
    densm = repmat(pars.dens,numel(pars.tv),numcols./numel(pars.dens)); % fix dens matrix
    out.densm = densm;
    
    % fix water mask
    watermask = (prp.waterdm==0);
    
    % calculate non-glacial erosion depth
    nonglacdm = cumtrapz(pars.tv,out.noicem.*nonglacEm.*1E-4.*watermask).*densm; % g/cm2

    % add potential burial depth
    if isfield(pars,'burialdepth') || ...
            (isfield(pars,'burialdepthv') && isfield(pars,'burialdepthv_tv'));
        if isfield(pars,'burialdepth') == 0; pars.burialdepth = 0; end;
        pars.burialdepth(1:numel(pars.tv),1:numcols) = pars.burialdepth;
        if isfield(pars,'burialdepthv') && isfield(pars,'burialdepthv_tv')
            minmaxv = fix_minmaxv(pars,'burialdepthv');
            for j = 1:size(minmaxv,1);
                pars.burialdepth(minmaxv(j,1):minmaxv(j,2),:) = pars.burialdepthv(j);
            end;
        end;
        nonglacdm = nonglacdm + pars.burialdepth.*densm; % add burial depth (g/cm2)
    end;
    
    % calculate glacial erosion depth matrix
    if isfield(glacEm,'rr');
        glacdrm = cumtrapz(pars.tv,out.icem.*glacEm.rr.*1E-4).*densm; % g/cm2
    end;
    if isfield(glacEm,'ri');
        % find deglaciation events
        deglacm = (filter([1 2],1,out.icem) == 1);
        glacdrm = glacdrm + cumsum(deglacm.*glacEm.ri).*densm; % g/cm2
    end;
    if isfield(glacEm,'ii');
        % find deglaciation events
        deglacm = (filter([1 2],1,out.icem) == 1);
        glacdim = cumsum(deglacm.*glacEm.ii).*densm; % g/cm2
    end;
    if isfield(glacEm,'ir');
        glacdim = glacdim + cumtrapz(pars.tv,out.icem.*glacEm.ir.*1E-4).*densm; % g/cm2
    end;
    
    % calculate full depth matrix (glac and nonglac plus water depth) and find deglac index
    if exist('glacdrm');
        out.fulldmr = glacdrm + nonglacdm + prp.waterdm;
        out.dmmr = (out.fulldmr-prp.waterdm)./densm.*1E-2; % depth matrix (m)
        out.gmmr = glacdrm./densm.*1E-2; % glacial erosion depth matrix (m)
    end;
    if exist('glacdim');
        out.fulldmi = glacdim + nonglacdm + prp.waterdm;
        out.dmmi = (out.fulldmi-prp.waterdm)./densm.*1E-2; % depth matrix (m)
        out.gmmi = glacdim./densm.*1E-2; % glacial erosion depth matrix (m)
    end;
% end subfunction getdepth =========================================================================


% subfunction fix_deglac ===========================================================================
function out = fix_deglac(deglac,icem,tv);
    % fix deglaciation in icem
    if numel(deglac) > 1;
        deglac = repmat(deglac(:)',size(tv));
        tv = repmat(tv,1,size(deglac,2));
        icem = repmat(icem,1,size(deglac,2)./size(icem,2));
    end;
    predeglac = (tv>=deglac);
    predeglacin = (cumsum(icem)>=1);
    newice = (predeglac-predeglacin == 1);
    out = icem.*predeglac + newice;
    % set ice at deglac (for cases where the last ice cover in icem starts later than deglac)
    [r,c] = find(predeglac==1);
    [cc,cidx] = unique(c);
    out(sub2ind(size(out),r(cidx),cc)) = 1;
% end subfunction fix_deglac =======================================================================


% subfunction fix_minmaxv ==========================================================================
function minmaxv = fix_minmaxv(pars,str);
    % refer tv to samplingyr if tv0 exists
    if isfield(pars,[str '_tv0']);
        pars.([str '_tv']) = pars.([str '_tv']) + pars.samplingyr - pars.([str '_tv0']);
    end;
    for j = 1:size(pars.([str '_tv']),1);
        if pars.([str '_tv'])(j,1) <= pars.tv(end) && pars.([str '_tv'])(j,2) >= pars.tv(1);
            minmaxv(j,1) = find(pars.tv >= pars.([str '_tv'])(j,1),1,'first');
            minmaxv(j,2) = find(pars.tv <= pars.([str '_tv'])(j,2),1,'last');
        end;
    end;
% end subfunction fix_minmaxv ======================================================================


% subfunction fixicem ==============================================================================
function icem = fixicem(pars,icem);
    % fix specified period(s) with ice cover
    if isfield(pars,'ice_tv');
        minmaxv = fix_minmaxv(pars,'ice');
        for j = 1:size(minmaxv,1);
            icem(minmaxv(j,1):minmaxv(j,2),:) = 1;
        end;
    end;
    % fix specified period(s) with no ice cover
    if isfield(pars,'noice_tv');
        minmaxv = fix_minmaxv(pars,'noice');
        for j = 1:size(minmaxv,1);
            icem(minmaxv(j,1):minmaxv(j,2),:) = 0;
        end;
    end;
% end subfunction fixicem ==========================================================================


% subfunction fixE =================================================================================
function [nonglacEm,glacEm,out] = fixE(pars,nonglacEm,glacEm,out);
    % fix specific nonglacial erosion in specified period(s)
    if isfield(pars,'erosionv') && isfield(pars,'erosionv_tv');
        minmaxv = fix_minmaxv(pars,'erosionv');
        for j = 1:size(minmaxv,1);
            nonglacEm(minmaxv(j,1):minmaxv(j,2),:) = pars.erosionv(j);
        end;
    end;
    % fix specific glacial erosion rate in specified period(s)
    if isfield(pars,'glacErv') && isfield(pars,'glacErv_tv');
        if isfield(glacEm,'ii'); glacEm.ir = zeros(size(glacEm.ii)); end;
        minmaxv = fix_minmaxv(pars,'glacErv');
        for j = 1:size(minmaxv,1);
            if isfield(glacEm,'rr'); glacEm.rr(minmaxv(j,1):minmaxv(j,2),:) = pars.glacErv(j); end;
            if isfield(glacEm,'ii'); glacEm.ir(minmaxv(j,1):minmaxv(j,2),:) = pars.glacErv(j); end;
            out.icetestr(minmaxv(j,1):minmaxv(j,2),:) = 0;
        end;
    end;
    % fix specific incremental glacial erosion step in specified period(s)
    if isfield(pars,'glacEiv') && isfield(pars,'glacEiv_tv');
        if isfield(glacEm,'rr'); glacEm.ri = zeros(size(glacEm.rr)); end;
        minmaxv = fix_minmaxv(pars,'glacEiv');
        for j = 1:size(minmaxv,1);
            if isfield(glacEm,'ii'); glacEm.ii(minmaxv(j,1):minmaxv(j,2),:) = pars.glacEiv(j); end;
            if isfield(glacEm,'rr'); glacEm.ri(minmaxv(j,1):minmaxv(j,2),:) = pars.glacEiv(j); end;
            out.icetesti(minmaxv(j,1):minmaxv(j,2),:) = 0;
        end;
    end;
% end subfunction fixE =============================================================================


% subfunction nuclE ================================================================================
function out = nuclE(pars,prp2,np,Evect,glaccalc);
    % find deglacial idx for minimum erosion estimation
    deglac_idx = find(prp2.icem(:,1)==1,1,'first');
    % define out structure
    out = struct;
    if isfield(prp2,'fulldmr'); % glacial erosion rate
        % calculate end conc and interpolate erosion
        Pm = getPm(pars,prp2,np.Psp,np.Pmu_z,np.Pref,prp2.fulldmr);
        out.Nendr = Ncalc(Pm,pars.tv,np.l);
        out.r = Ecalc(out.Nendr,Evect,np.N);
        % calculate deglacial N with 1% margin for estimating minimum erosion
        Ndegl = Ncalc(Pm(1:deglac_idx,1),pars.tv(1:deglac_idx),np.l);
        out.rmin = Ecalc(out.Nendr,Evect,Ndegl*1.01);
        if numel(deglac_idx) == 0; out.rmin = out.r; end;
    end;
    if isfield(prp2,'fulldmi'); % incremental glacial erosion steps
        % calculate end conc and interpolate erosion
        Pm = getPm(pars,prp2,np.Psp,np.Pmu_z,np.Pref,prp2.fulldmi);
        out.Nendi = Ncalc(Pm,pars.tv,np.l);
        out.i = Ecalc(out.Nendi,Evect,np.N);
        % calculate deglacial N with 1% margin for estimating minimum erosion
        Ndegl = Ncalc(Pm(1:deglac_idx,1),pars.tv(1:deglac_idx),np.l);
        out.imin = Ecalc(out.Nendi,Evect,Ndegl*1.01);
        if numel(deglac_idx) == 0; out.imin = out.i; end;
    end;
    out.Evect = Evect;
% end subfunction nuclE ============================================================================


% subfunction getPmc ===============================================================================
function pmc = getPmc(pars,prp2,Psp,Pmu_z,Pref);
% build all loop-invariant production terms for getPm; computed once per (pars,prp2,nuclide) pair
% and reused for every depth matrix (e.g. the QMC while-loop iterations in Eunccalc)
    ncol = size(prp2.noicem,2);
    Lspm = repmat(pars.Lsp,1,ncol); % fix Lsp matrix
    % fix thickness scaling factor
    if pars.thick > 0;
        thickSF = (Lspm./(prp2.densm.*pars.thick)) .* ...
            (1 - exp(((-1.*prp2.densm.*pars.thick)./Lspm)));
    else; thickSF = 1; end;
    % fix shielding matrix if having a time-varying shielding factor
    if numel(pars.shield) > 1; shieldm = repmat(pars.shield,1,ncol); else; shieldm = pars.shield; end;
    % fold spallation P, Pref, thickSF, shield and noicem into one constant matrix
    Prefm = repmat(Pref,numel(pars.tv),ncol./numel(Pref));
    PspPrf = repmat(Psp(:),1,ncol) .* Prefm .* thickSF .* shieldm .* prp2.noicem;
    % time with production
    simp = (repmat(pars.tv,[1 ncol]) <= ...
        repmat(pars.simt,[numel(pars.tv) ncol/numel(pars.simt)]));
    % muon P: cache the pchip pp-form on the z_mu table; getPm evaluates it directly with the same
    % lookup + Horner arithmetic as Octave's ppval (= original interp1(...,'pchip')), reusing the
    % pp-form instead of rebuilding it on every depth matrix
    pp = pchip(pars.z_mu',Pmu_z');
    pmc.brk = pp.breaks(:); % force column (Octave pchip returns a row) to keep idx/h column-shaped
    pmc.coefs = pp.coefs; % n x 4, highest power first (as in ppval)
    pmc.zmu_hi = pars.z_mu(end); % beyond this the original interp1 gave NaN (->0 after NaN fix)
    pmc.Lspm = Lspm;
    pmc.PspPrf = PspPrf;
    pmc.simp = simp;
    pmc.munmask = shieldm .* prp2.noicem .* simp; % mask for muon P
    pmc.thickdens = pars.thick .* prp2.densm ./ 2; % depth offset for muon P (constant per run)
% end subfunction getPmc ===========================================================================

% subfunction getPm ================================================================================
function Pm = getPm(pars,prp2,Psp,Pmu_z,Pref,depthm,pmc);
    % build production cache of loop-invariant terms if not passed in
    if nargin < 7; pmc = getPmc(pars,prp2,Psp,Pmu_z,Pref); end;
    dpfs = exp(-depthm./pmc.Lspm); % spal depth dependence matrix
    Pspm = pmc.PspPrf .* pmc.simp .* dpfs; % spal prod matrix
    % muon P: evaluate the cached pchip pp-form at the muon center depth with exactly the same
    % arithmetic as Octave's ppval (lookup + Horner), but without ppval's reshape/shiftdim
    % overhead - results are identical to the original interp1(...,'pchip') muon term
    dkm = depthm + pmc.thickdens; % muon P center depth
    oob = dkm > pmc.zmu_hi; % beyond the z_mu table: original code got NaN, then set it to 0
    idx = lookup(pmc.brk, dkm, 'lr'); % same interval search as ppval, matrix-shaped like dkm
    h = (dkm - pmc.brk(idx))(:); % work on columns; per-element arithmetic identical to ppval
    Pc = pmc.coefs(idx(:),:); % pchip coefficients (highest power first, as in ppval)
    mu = Pc(:,1);
    mu = mu .* h + Pc(:,2);
    mu = mu .* h + Pc(:,3);
    mu = mu .* h + Pc(:,4);
    Pmum = reshape(mu,size(dkm));
    Pmum(oob) = 0; % replicate the original NaN->0 fix
    Pmum = Pmum .* pmc.munmask; % muon prod matrix
    Pm = Pspm + Pmum; % production matrix
% end subfunction getPm ============================================================================


% subfunction Ncalc ================================================================================
function Nend = Ncalc(P,tv,l);
    % Calculate N(end) including decay and erosion
    lm = repmat(l,numel(tv),size(P,2)./numel(l)); % lambda matrix
    tm = repmat(tv,1,size(P,2)); % time matrix
    dcf = exp(-tm.*lm); % decay factor
    Nend = trapz(tv,P.*dcf);
% end subfunction Ncalc ============================================================================


% subfunction Ecalc ================================================================================
function NE = Ecalc(Nend,Evect,N);
    % remove points with same N from Nend and Evect
    rmidx = find(diff(Nend)==0,1,'first');
    if numel(rmidx)>0; Nend(rmidx+1:end) = []; Evect(rmidx+1:end) = []; end;
    % set too low and too high N to limits of Nend
    N(N<min(Nend)) = min(Nend);
    N(N>max(Nend)) = max(Nend);
    % interpolate erosion
    if numel(Nend)>1;
        NE = interp1(Nend,Evect,N,'pchip');
    else; % sample is saturated and we set the erosion to zero
        NE = 0;
    end;
% end subfunction Ecalc ============================================================================


% subfunction uncqmc ===============================================================================
function pars = uncqmc(pars,sample,qmc,dimn);
    % uncertainty parameter string
    ustr = {'icevalue','erosion','glacE','isostsubm','deglac','dens'};
    pars.ustr = {};
    sample.isostsubm = 1; % fix for isostsubm
    % loop for uncertainties
    for i = 1:numel(ustr);
        minmax = 0;
        % fix parameters with min/max limits
        if isfield(sample,[ustr{i} 'min']) && isfield(sample,[ustr{i} 'max']);
            % evenly spaced parameters via rank-1 lattice
            pars.(ustr{i}) = sample.([ustr{i} 'min']) + ...
                (sample.([ustr{i} 'max'])-sample.([ustr{i} 'min']))*latt_u(qmc,dimn)';
            pars.ustr(end+1) = ustr(i);
            dimn = dimn + 1;
            minmax = 1;
        end;
        % fix parameters with normal uncertainties
        if isfield(sample,[ustr{i} 'unc']);
            % test if both min/max and +/- unc is present in input
            if minmax == 1;
                fprintf(1,'\nINPUT ERROR! This sample has both min/max and +/- uncertainty for parameter %s',ustr{i});
                dimn = dimn - 1;
            end;
            sz = size(sample.(ustr{i}),1);
            normv = norminv(latt_u(qmc,dimn),0,1);
            % deterministic rank-1 lattice low-discrepancy sampling
            uncv = normv';
            unc_qmcm = repmat(uncv.*sample.([ustr{i} 'unc']),sz,1);
            pars.(ustr{i}) = repmat(sample.(ustr{i}),1,qmc) + unc_qmcm;
            pars.(ustr{i})(pars.(ustr{i})<0) = 0; % don't allow negative values
            if any(strcmp(ustr{i},pars.ustr))==0;
                pars.ustr(end+1) = ustr(i);
            end;
            dimn = dimn + 1;
        end;
    end;

    for i = 1:numel(ustr);

    end;
    % fix for isostsubm
    if isfield(pars,'isostsubm');
        pars.delv = repmat(pars.delv,1,qmc) .* repmat(pars.isostsubm,size(pars.delv));
    end;
    % fix for simt
    if isfield(sample,'simtqmc');
        pars.simt = sample.simtqmc;
        pars.ustr(end+1) = 'simt';
    end;
% end subfunction uncqmc ===========================================================================


% subfunction latt_u ===============================================================================
function u = latt_u(n, dim);
% deterministic rank-1 (Korobov) lattice points in (0,1) for a given dimension.
% u_k = (mod((k-1)*z_dim, n) + 0.5) / n  with  z_dim = a^(dim-1) mod n.
% The multiplier 'a' is taken from an offline-optimized table (minimises the worst pairwise cross-
% dimension rank correlation over parameter dims 1..27), so every parameter margin is exactly one
% point per stratum and parameters are decorrelated in a deterministic, reproducible QMC fashion.
% Developed with the help of OpenCode Big Pickle model, September 2026.
    ns = 2 .^ (6:16);
    atab = [3 5 51 41 837 583 3729 2183 8381 16759 33087];
    if nargin < 2; dim = 1; end;
    if dim < 1; dim = 1; end;
    k = round(log2(n));
    ki = k - 5;
    if ki >= 1 && ki <= numel(ns);
        n0 = ns(ki);
        a = atab(ki);
    else;
        [n0, a] = latt_fallback(n, ns, atab);
    end;
    if n0 ~= n;
        warning('glacialE:qmcNAS',['qmc=' num2str(n) ' is not a power of two; using nearest table point ' num2str(n0)]);
    end;
    % z = a^(dim-1) mod n0
    z = 1;
    for j = 2:dim;
        z = mod(z*a, n0);
    end;
    if mod(z,2) == 0; z = z + 1; end;  % keep odd (coprime with power of 2)
    z = mod(z, n0); if z == 0; z = 1; end;
    u = (mod((0:n-1)'*z, n0) + 0.5) / n0;
% end subfunction latt_u ===========================================================================


% subfunction latt_fallback ========================================================================
function [n0, a] = latt_fallback(n, ns, atab);
% choose the nearest power of two in the table and its multiplier
    [~, ii] = min(abs(ns - n));
    n0 = ns(ii);
    a = atab(ii);
% end subfunction latt_fallback ====================================================================


% subfunction uncfix ===============================================================================
function [output,pl,unc] = uncfix(output,pl,outcol,plch,uncsens,sensstr,sample,pars,parsu,prp1,...
    prpunc,prpuncsens,Etestv,tdv,gdv,NE,np,Prefv,lv,glaccalc,i,qmc,nucl);
    % fix strings
    if nucl == 10; nstr = '10'; nstr2 = '10Be';
    elseif nucl == 26; nstr = '26'; nstr2 = '26Al';
    elseif nucl == 14; nstr = '14'; nstr2 = '14C'; end;
    % fix Pref and l
    parsu.Pref = Prefv; parsu.l = lv;
    % define unc structure
    unc = struct;
    % calculate uncertainty for glacial erosion rate case
    if isfield(NE,'r');
        unc = Eunccalc(unc,NE,prpunc,parsu,np,'r',glaccalc,qmc);
    end;
    % calculate uncertainty for incremental glacial erosion case
    if isfield(NE,'i');
        unc = Eunccalc(unc,NE,prpunc,parsu,np,'i',glaccalc,qmc);
    end;
    % display output and fill for saving
    if isfield(NE,'r');
        [output,pl] = fixoutplot(output,pl,sample,NE,unc,nucl,'r',glaccalc,i,outcol,tdv,gdv,plch,...
            pars,parsu,prp1,prpunc,Etestv,np);
    end
    if isfield(NE,'i');
        [output,pl] = fixoutplot(output,pl,sample,NE,unc,nucl,'i',glaccalc,i,outcol,tdv,gdv,plch,...
            pars,parsu,prp1,prpunc,Etestv,np);
    end;
    % estimate uncertainty sensitivity
    if numel(sensstr) > 0;        
        if any(strcmp('full',sensstr)); snum = numel(sensstr)-1; else; snum = numel(sensstr); end;
        fprintf(1,'\n%s sensitivity for %.0f parameters:',nstr2,snum);
        sensn = Euncsens(pars,parsu,prpuncsens,unc,np,NE,uncsens,sample,qmc,glaccalc);
        if isfield(NE,'r');
            [output,pl.(['sens' nstr 'r'])] = Euncsens_outpl(output,outcol.(['E' nstr 'rsens1']),...
                outcol.(['E' nstr 'rsens2']),pl.(['sens' nstr 'r']),sensn,sensstr,NE,unc,'r');
        end;
        if isfield(NE,'i');
            [output,pl.(['sens' nstr 'i'])] = Euncsens_outpl(output,outcol.(['E' nstr 'isens1']),...
                outcol.(['E' nstr 'isens2']),pl.(['sens' nstr 'i']),sensn,sensstr,NE,unc,'i');
        end;
    end;
% end subfunction uncfix ===========================================================================


% subfunction Eunccalc =============================================================================
function out = Eunccalc(out,NE,prpunc,parsu,np,ri,glaccalc,qmc);
% calculate erosion uncertainty based on input parameters in parsu (np for N)
    % fix glacErate/glacEstep
    if strcmp(ri,'r'); glacErate = 1; glacEstep = 0;
    elseif strcmp(ri,'i'); glacErate = 0; glacEstep = 1; end;
    
    % fix erosion parameter to test
    if glaccalc == 1; parsu.glacE = NE.(ri); else; parsu.erosion = NE.(ri); end;
    
    % calculate depth matrix
    prpunc2 = getdepth(prpunc,parsu,glacErate,glacEstep,qmc);
    
    % fill Em
    Em = repmat(NE.(ri),1,qmc);
    
    % calculate Nend
    % build production cache (invariant across the depth-matrix updates below, reused by getPm)
    pmc = getPmc(parsu,prpunc2,np.Psp,np.Pmu_z,parsu.Pref);
    Pm = getPm(parsu,prpunc2,np.Psp,np.Pmu_z,parsu.Pref,prpunc2.(['fulldm' ri]),pmc);
    
    % calculate decay factor matrix
    lm = repmat(parsu.l,numel(parsu.tv),qmc./numel(parsu.l)); % lambda matrix
    tm = repmat(parsu.tv,1,qmc); % time matrix
    dcf = exp(-tm.*lm); % decay factor
    
    % calculate Nend including decay and erosion
    Nend = trapz(parsu.tv,Pm.*dcf);

    % set erosion number two based on mid E (Em) and ratio between calculated and expected N
    Em(end+1,:) = Em .* Nend./np.Nv;
    
    % fix erosion parameter to test
    if glaccalc == 1; parsu.glacE = Em(end,:); else; parsu.erosion = Em(end,:); end;
    % calculate depth matrix
    prpunc2 = getdepth(prpunc,parsu,glacErate,glacEstep,qmc);
    % calculate production and Nend
    Pm = getPm(parsu,prpunc2,np.Psp,np.Pmu_z,parsu.Pref,prpunc2.(['fulldm' ri]),pmc);
    Nend(end+1,:) = trapz(parsu.tv,Pm.*dcf);
    % mean deviation from N used as test parameter
    Ntest = mean(abs(Nend(end,:)./np.Nv-1));
    
    % while loop for linear interpolation of E
    j = 0; % counter for while loop
    while Ntest(end) > 0.005;
        j = j+1; % counter
        % estimate erosion based on linear interpolation from two last calculated Em and Nend
        Em(end+1,:) = min(Em(end-1,:) + (np.Nv-Nend(end-1,:)).*(Em(end,:)-Em(end-1,:))./...
            (Nend(end,:)-Nend(end-1,:)),max(NE.Evect));
        % don't allow too large change of E per step
        Etemp1 = Em(end-1,:); Etemp2 = Em(end,:);
        Etemp2(Etemp2>Etemp1.*2) = min(Etemp1(Etemp2>Etemp1.*2).*2,max(NE.Evect));
        Etemp2(Etemp2<Etemp1./2) = min(Etemp1(Etemp2<Etemp1./2)./2,max(NE.Evect));
        Etemp2(isnan(Etemp2)) = NE.(ri);
        Em(end,:) = Etemp2;
        % fix erosion parameter to test
        if glaccalc == 1; parsu.glacE = Em(end,:); else; parsu.erosion = Em(end,:); end;
        % calculate depth matrix
        prpunc2 = getdepth(prpunc,parsu,glacErate,glacEstep,qmc);
        % calculate production and Nend
        Pm = getPm(parsu,prpunc2,np.Psp,np.Pmu_z,parsu.Pref,prpunc2.(['fulldm' ri]),pmc);
        Nend(end+1,:) = trapz(parsu.tv,Pm.*dcf);
        % mean deviation from N used as test parameter (max 0.5%)
        Ntest(end+1) = mean(abs(Nend(end,:)./np.Nv-1));
        % after 5 loops: pick E yielding the lowest Ntest
        if j == 5; [mintest idx] = min(Ntest); Em = Em(idx+1,:); Ntest = 0; end;
    end;
    
    % pick out best fit E vector
    Ev = Em(end,:);
    Ev(isnan(Ev)) = NE.(ri); % set NaN E to midpoint E
    
    % full-sample 68.27% two-sided percentiles (84.135 / 15.865) of Ev
    out.([ri 'pos']) = prctile(Ev,84.135)-NE.(ri);
    out.([ri 'neg']) = NE.(ri)-min(prctile(Ev,15.865),NE.([ri 'min']));

    % save Ev in output
    out.(['Ev' ri]) = Ev;

    % don't allow too large positive unc
    if NE.(ri)+out.([ri 'pos']) >= max(NE.Evect);
        out.([ri 'pos']) = max(NE.Evect);
    end;
% end subfunction Eunccalc =========================================================================


% subfunction Euncsens =============================================================================
function sens = Euncsens(pars,parsu,prpuncsens,uncp,np,NE,uncsens,sample,qmc,glaccalc);
% function for estimating the uncertainty sensitivity to various parameters
    % fix for Pref and l
    pars.Pref = np.Pref;
    pars.l = np.l;
    % fix for waterdm
    if isfield(prpuncsens,'waterdm'); pars.waterdm = prpuncsens.waterdm; end;
    % strings for calculations
    str1 = {'N','Pref','l'};
    str2 = {'dens','erosion','glacE','simt','icevalue','deglac','isostsubm'};
    parn = 0; % counter
    % calculate N/Pref/l sensitivity
    for j = 1:numel(str1);
        if uncsens.(str1{j}) == 1;
            parn = parn+1; fprintf(1,' %.0f',parn);
            sens.(str1{j}) = Euncsens2(pars,parsu,prpuncsens,uncp,np,NE,qmc,glaccalc,str1{j});
        end;
    end;
    % calculate sensitivity for all other parameters
    for j = 1:numel(str2);
        if uncsens.(str2{j}) == 1 && any(strcmp(str2{j},parsu.ustr));
            parn = parn+1; fprintf(1,' %.0f',parn);
            sens.(str2{j}) = Euncsens2(pars,parsu,prpuncsens,uncp,np,NE,qmc,glaccalc,str2{j});
        end;
    end;
    fprintf(1,' done!');
% end subfunction Euncsens =========================================================================


% subfunction Euncsens2 ============================================================================
function out = Euncsens2(pars,parsu,prpuncsens,uncp,np,NE,qmc,glaccalc,parstr);
    % fix for isostsubm
    if strcmp(parstr,'isostsubm'); parstr = 'delv'; end;
    % fix for N sensitivity
    if strcmp(parstr,'N')==0; np.Nv = np.N; end;
    % fix sensitivity parameter
    if any(strcmp(parstr,{'Pref','l','dens','erosion','glacE','simt','icevalue','deglac','delv'}));
        pars.(parstr) = parsu.(parstr);
    end;
    % fix precalc for sensitivity with varying ice cover or isostatic uplift
    if any(strcmp(parstr,{'deglac','icevalue','delv'}));
        if isfield(pars,'waterdm'); pars = rmfield(pars,'waterdm'); end;
        prpuncsens = precalc(pars,qmc);
    end;
    % define out structure
    out = struct;
    % erosion rate case
    if isfield(NE,'r'); out = Eunccalc(out,NE,prpuncsens,pars,np,'r',glaccalc,qmc); end;
    % incremental depth step case
    if isfield(NE,'i'); out = Eunccalc(out,NE,prpuncsens,pars,np,'i',glaccalc,qmc); end;
% end subfunction Euncsens2 ========================================================================


% subfunction Euncsens_outpl =======================================================================
function [output,plsens] = Euncsens_outpl(output,col1,col2,plsens,sens,sensstr,NE,unc,ri);
% fix output for uncertainty sensitivity output
    outp = []; outn = []; outstr = {};
    if strcmp(sensstr(end),'full'); jend = numel(sensstr)-1; else; jend = numel(sensstr); end;
    for j = 1:jend;
        outp(end+1) = sens.(sensstr{j}).([ri 'pos']);
        outn(end+1) = sens.(sensstr{j}).([ri 'neg']);
        outstr(end+1) = {num2str(outp(end),'%.2f')};
        outstr(end+1) = {num2str(outn(end),'%.2f')};
    end;
    if strcmp(sensstr(end),'full');
        outp(end+1) = unc.([ri 'pos']);
        outn(end+1) = unc.([ri 'neg']);
    end;
    if NE.(ri) > 0;
        plsens.pos(end+1,:) = 1 + outp./NE.(ri);
        plsens.neg(end+1,:) = 1 - outn./NE.(ri);
    else;
        plsens.pos(end+1,:) = 100;
        plsens.neg(end+1,:) = 0;
    end;
    output(end,col1:col2) = outstr;
% end subfunction Euncsens_outpl ===================================================================


% subfunction P1026fix =============================================================================
function [pl,P10u,P26u] = P1026fix(pl,sample,LSDfix,Pmu_z,z,np10,np26);
% calculate normalized P10 and P26
    % interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
    Lsp = rawattenuationlength(sample.pressure(1),LSDfix.RcEst);
    % fix thickness scaling factor
    if sample.thick > 0;
        thickSF = Lsp./(sample.dens.*sample.thick) .* (1-exp(-1.*sample.dens.*sample.thick./Lsp));
    else; thickSF = 1; end;
    % fix shielding factor
    shield = sample.shield;
    if numel(shield) > 1; shield = trapz(sample.tv,sample.shield)./sample.tv(end); end;
    % sample production from spallation
    P0sp10 = np10.Psp(1).*np10.Pref.*thickSF.*shield;
    P0sp26 = np26.Psp(1).*np26.Pref.*thickSF.*shield;
    % sample production from muons
    P0mu10 = interp1(z,Pmu_z.mu10,sample.thick.*sample.dens./2,'pchip') .* shield;
    P0mu26 = interp1(z,Pmu_z.mu26,sample.thick.*sample.dens./2,'pchip') .* shield;
    % mean sample P
    P10u = P0sp10 + P0mu10;
    P26u = P0sp26 + P0mu26;
    % depth P saved in pl
    if isfield(pl,'banana');
        pl.banana.P10z(end+1,:) = (np10.Psp(1).*np10.Pref.*exp(-z./Lsp) + Pmu_z.mu10) .* shield;
        pl.banana.P26z(end+1,:) = (np26.Psp(1).*np26.Pref.*exp(-z./Lsp) + Pmu_z.mu26) .* shield;
    end;
% end subfunction P1026fix =========================================================================


% subfunction Ecalc1026 ============================================================================
function [output,pl,E1026] = Ecalc1026(output,pl,E1026,outcol,plch,sample,pars,parsu,prp1,prpunc,...
    prpuncsens,Etestv,N10E,N26E,unc10,unc26,np10,np26,Pref10v,Pref26v,l10v,l26v,P10u,P26u,tdv,...
    gdv,glaccalc,qmc,i,ri);
    % calculate 10Be and 26Al uncertainty with N, Pref, and l uncertainties
    pars.Pref = Pref10v; pars.l = l10v;
    Eunc10 = Eunccalc(E1026,N10E,prpuncsens,pars,np10,ri,glaccalc,qmc);
    pars.Pref = Pref26v; pars.l = l26v;
    Eunc26 = Eunccalc(E1026,N26E,prpuncsens,pars,np26,ri,glaccalc,qmc);
    % calculate weighted E using internal uncertainties and the EVM method
    [E1026.(ri),intp,intn] = evm([N10E.(ri),N26E.(ri)],...
        [Eunc10.([ri 'pos']),Eunc26.([ri 'pos'])],[Eunc10.([ri 'neg']),Eunc26.([ri 'neg'])]);
    % calculate mean E for external uncertainty estimation - use arithmetic mean of 10/26 values
    E1026.(['Ev' ri]) = mean([unc10.(['Ev' ri]);unc26.(['Ev' ri])]);
    % full-sample 68.27% two-sided percentiles (84.135 / 15.865) of Ev
    E1026.([ri 'pos']) = prctile(E1026.(['Ev' ri]),84.135)-E1026.(ri);
    E1026.([ri 'neg']) = E1026.(ri)-prctile(E1026.(['Ev' ri]),15.865);
    % don't allow too large positive unc
    if E1026.(ri)+E1026.([ri 'pos']) >= max(Etestv);
        E1026.([ri 'pos']) = max(Etestv);
    end;
    % display output and fill for saving
    np1026 = struct; % fix for fixoutplot
    [output,pl] = fixoutplot(output,pl,sample,E1026,E1026,1026,ri,glaccalc,i,outcol,tdv,gdv,plch,...
        pars,parsu,prp1,prpunc,Etestv,np1026);
    % test overlap and save in output
    if N10E.(ri)-Eunc10.([ri 'neg']) <= N26E.(ri)+Eunc26.([ri 'pos']) && ...
            N10E.(ri)+Eunc10.([ri 'pos']) >= N26E.(ri)-Eunc26.([ri 'neg']) && ...
            (E1026.(ri)+E1026.([ri 'pos'])+E1026.([ri 'neg'])) > 0;
        output(i+1,outcol.(['E1026' ri 'overlap'])) = {'yes'};
        E1026.(['overlap' ri]) = 1; str = '';
    else;
        output(i+1,outcol.(['E1026' ri 'overlap'])) = {'no'};
        fprintf(1,'   no 10Be/26Al overlap');
        E1026.(['overlap' ri]) = 0; str = '_out';
    end;
    % calculate Npath for banana
    if plch.banana==1 && isfield(E1026,ri);
        % fix mid-point erosion parameter
        if glaccalc == 1; pars.glacE = E1026.(ri); else; pars.erosion = E1026.(ri); end;
        % calculate depth for mid-point E
        if strcmp(ri,'r'); prp2 = getdepth(prp1,pars,1,0,1); end;
        if strcmp(ri,'i'); prp2 = getdepth(prp1,pars,0,1,1); end;
        % calculate P and Npath and fill pl
        Pmid10 = getPm(pars,prp2,np10.Psp,np10.Pmu_z,np10.Pref,prp2.(['fulldm' ri]));
        Pmid26 = getPm(pars,prp2,np26.Psp,np26.Pmu_z,np26.Pref,prp2.(['fulldm' ri]));
        pl.banana.(['Npath10' ri str])(1:numel(Pmid10),end+1) = ...
            Npathcalc(pars.tv,Pmid10,np10.l)./P10u;
        pl.banana.(['Npath26' ri str])(1:numel(Pmid26),end+1) = ...
            Npathcalc(pars.tv,Pmid26,np26.l)./P26u;
    end
% end subfunction Ecalc1026 ========================================================================


% subfunction depthPconc ===========================================================================
function [output,pl] = depthPconc(output,pl,pars,parsu,prp1,prpunc,Emid,uncp,uncn,Ev,Etestv,np,...
    plch,tdv,gdv,glaccalc,i,outcol,ri,nstr);
% fix and calculate sample depth/production/concentration over simulation time for plotting/output
    % fix glacErate/glacEstep
    if strcmp(ri,'r'); glacErate = 1; glacEstep = 0;
    elseif strcmp(ri,'i'); glacErate = 0; glacEstep = 1; end;
    % fix mid-point erosion parameter for depth matrix
    if glaccalc == 1; pars.glacE = Emid; else; pars.erosion = Emid; end;
    % calculate depth for mid-point E
    prp2 = getdepth(prp1,pars,glacErate,glacEstep,1);
    deglidx = find(prp2.(['icetest' ri])==1,1,'first');
    Ev(Ev<0) = 0; % set negative E to 0
    % fix erosion parameter for uncertainty depth matrix
    if glaccalc == 1; parsu.glacE = Ev; else; parsu.erosion = Ev; end;
    % calculate depth matrix
    prp2u = getdepth(prpunc,parsu,glacErate,glacEstep,numel(Ev));
    deglidxunc = find(sum(prp2u.(['icetest' ri])')>0,1,'first');
    % plot depth or get depth at specific time? ===========================
    if plch.depth == 1 || numel(tdv) > 0;
        % fix depth matrices
        dvm = prp2.(['dmm' ri]);
        dmm = prp2u.(['dmm' ri]);

        % calculate total depths
        [output,dvm,depthp,depthn] = tdv_fix(output,dvm,dmm,Emid,uncp,uncn,Etestv,deglidx,...
            deglidxunc,tdv,pars,outcol,nstr,ri,i,'d');

        % save depth and pos and neg uncertainty for plotting
        if plch.depth == 1;
            pl.([ri nstr 'd']).y(1:numel(dvm),end+1) = dvm;
            pl.([ri nstr 'd']).yunc(1:numel(dvm)*2,end+1) = [flip(depthn);depthp];
            % save if sample has maximum possible erosion (for splitting depth lines)
            pl.([ri nstr 'd']).emax(end+1) = (Emid == max(Etestv));
            pl.([ri nstr 'd']).t(1:numel(pars.tv),end+1) = pars.tv.*1E-3;
            pl.([ri nstr 'd']).tunc(1:numel(pars.tv)*2,end+1) = [flip(pars.tv.*1E-3);pars.tv.*1E-3];
            pl.([ri nstr 'd']).m(1:numel(pars.tv),end+1) = 1;
            pl.([ri nstr 'd']).munc(1:numel(pars.tv)*2,end+1) = 1;
            if plch.probdensmap == 1;
                pl.([ri nstr 'd']) = pdm_fix(pl.([ri nstr 'd']),dmm,pars.tv,plch,plch.maxd);
            end;
            if plch.pct_area == 1 || any(plch.pct_lines > 0);
                % interpret mm to plotting time
                if plch.pct_unc == 1;
                    mm_plot = interp1(pars.tv,dmm,pl.([ri nstr 'd']).prcnttv);
                else;
                    mm_plot = interp1(pars.tv,dvm(:),pl.([ri nstr 'd']).prcnttv);
                end;
                mm_plot(mm_plot>plch.maxd) = plch.maxd;
                if plch.pct_unc == 1;
                    pl.([ri nstr 'd']).prcntm(end+1:end+size(dmm,2),:) = mm_plot';
                else;
                    pl.([ri nstr 'd']).prcntm(end+1,:) = mm_plot';
                end;
            end;
        end;
    end;
    % get glacial erosion depth at specific time? ===========================
    if numel(gdv) > 0 && glaccalc == 1;
        % fix glacpars
        glacparcu = pars;
        glacparcu.glacE = Ev;
        % calculate depth matrix
        glacprp2u = getdepth(prp1,glacparcu,glacErate,glacEstep,numel(Ev));
        deglidxunc = find(sum(glacprp2u.(['icetest' ri])')>0,1,'first');
        % fix glacial depth matrices
        gvm = prp2.(['gmm' ri]);
        gmm = glacprp2u.(['gmm' ri]);

        % calculate glacial erosion depths
        [output,gvm,depthp,depthn] = tdv_fix(output,gvm,gmm,Emid,uncp,uncn,Etestv,deglidx,...
            deglidxunc,gdv,pars,outcol,nstr,ri,i,'gd');
    end;


    % fix for P and N =====================================================
    if any(strcmp(nstr,{'10','26','14'})) && (plch.P == 1 || plch.N);
        % calculate mid and unc P
        Pmid = getPm(pars,prp2,np.Psp,np.Pmu_z,np.Pref,prp2.(['fulldm' ri]));
        Pmunc = getPm(parsu,prp2u,np.Psp,np.Pmu_z,parsu.Pref,prp2u.(['fulldm' ri]));
    end;
    % plot P ==============================================================
    if any(strcmp(nstr,{'10','26','14'})) && plch.P==1;
        % calculate pos and neg uncertainties for P
        Pmidm = repmat(Pmid',size(Pmunc,2),1);
        Pmuncp = Pmunc'; Pmuncp(Pmuncp<=Pmidm) = NaN;
        Pmuncn = Pmunc'; Pmuncn(Pmuncn>Pmidm) = NaN;
        Pp = prctile(Pmuncp,68.27)';
        Pn = prctile(Pmuncn,31.73)';

        % normalize to P at t = 0
        if plch.Pprcnt == 1; P0 = Pmid(1)/100; else; P0 = 1; end;

        % save P and pos and neg uncertainty for plotting
        pl.([ri nstr 'P']).y(1:numel(Pmid),end+1) = Pmid./P0;
        pl.([ri nstr 'P']).yunc(1:numel(Pmid)*2,end+1) = [flip(Pp);Pn]./P0;
        pl.([ri nstr 'P']).t(1:numel(pars.tv),end+1) = pars.tv.*1E-3;
        pl.([ri nstr 'P']).tunc(1:numel(pars.tv)*2,end+1) = [flip(pars.tv.*1E-3);pars.tv.*1E-3];
        pl.([ri nstr 'P']).m(1:numel(pars.tv),end+1) = 1;
        pl.([ri nstr 'P']).munc(1:numel(pars.tv)*2,end+1) = 1;
        if plch.probdensmap == 1;
            pl.([ri nstr 'P']) = pdm_fix(pl.([ri nstr 'P']),Pmunc,pars.tv,plch,max(max(Pmunc)));
        end;
        if plch.pct_area == 1 || any(plch.pct_lines > 0);
            % interpret mm to plotting time
            if plch.pct_unc == 1;
                mm_plot = interp1(pars.tv,Pmunc,pl.([ri nstr 'P']).prcnttv);
                pl.([ri nstr 'P']).prcntm(end+1:end+size(Pmunc,2),:) = mm_plot';
            else;
                mm_plot = interp1(pars.tv,Pmid(:),pl.([ri nstr 'P']).prcnttv);
                pl.([ri nstr 'P']).prcntm(end+1,:) = mm_plot';
            end;
        end;
    end;
    % plot N ==============================================================
    if any(strcmp(nstr,{'10','26','14'})) && plch.N==1;
        % normalize to N at t = 0
        if plch.Nprcnt == 1; N0 = Nmid(1)/100; else; N0 = 1; end;

        % calculate N and uncertainty through simulation tv
        Nmid = Npathcalc(pars.tv,Pmid,np.l)';
        Nmunc = Npathcalc(pars.tv,Pmunc,np.l)';
        Nmidm = repmat(Nmid,size(Nmunc,1),1);
        Nmuncp = Nmunc; Nmuncp(Nmuncp<Nmidm) = NaN;
        Nmuncn = Nmunc; Nmuncn(Nmuncn>=Nmidm) = NaN;
        Np = prctile(Nmuncp,68.27)' ./ N0;
        Nn = prctile(Nmuncn,31.73)' ./ N0;
        
        % save N and pos and neg uncertainty for plotting
        pl.([ri nstr 'N']).y(1:numel(Nmid),end+1) = Nmid ./ N0;
        pl.([ri nstr 'N']).yunc(1:numel(Pmid)*2,end+1) = [flip(Np);Nn];
        pl.([ri nstr 'N']).t(1:numel(pars.tv),end+1) = pars.tv.*1E-3;
        pl.([ri nstr 'N']).tunc(1:numel(pars.tv)*2,end+1) = [flip(pars.tv.*1E-3);pars.tv.*1E-3];
        pl.([ri nstr 'N']).m(1:numel(pars.tv),end+1) = 1;
        pl.([ri nstr 'N']).munc(1:numel(pars.tv)*2,end+1) = 1;
        if plch.probdensmap == 1;
            pl.([ri nstr 'N']) = pdm_fix(pl.([ri nstr 'N']),Nmunc',pars.tv,plch,max(max(Nmunc)));
        end;
        if plch.pct_area == 1 || any(plch.pct_lines > 0);
            % interpret mm to plotting time
            if plch.pct_unc == 1;
                mm_plot = interp1(pars.tv,Nmunc',pl.([ri nstr 'N']).prcnttv);
                pl.([ri nstr 'N']).prcntm(end+1:end+size(Nmunc,1),:) = mm_plot';
            else;
                mm_plot = interp1(pars.tv,Nmid(:),pl.([ri nstr 'N']).prcnttv);
                pl.([ri nstr 'N']).prcntm(end+1,:) = mm_plot';
            end;
        end;
    end;
% end subfunction depthPconc =======================================================================


% subfunction tdv_fix ==============================================================================
function [output,dvm,depthp,depthn] = tdv_fix(output,dvm,dmm,Emid,uncp,uncn,Etestv,deglidx,deglidxunc,tdv,pars,...
    outcol,nstr,ri,i,gd);
    % fix for unconstrained erosion with constant depth at 1 km
    if Emid == max(Etestv); dvm(deglidx:end) = 1E3; end;

    % calculate pos and neg uncertainties for depth history
    dvmm = repmat(dvm',size(dmm,2),1);
    dmmp = dmm'; dmmp(dmmp<dvmm) = NaN;
    dmmn = dmm'; dmmn(dmmn>=dvmm) = NaN;
    dvmp = prctile(dmmp,68.27) - dvm';
    dvmn = dvm' - prctile(dmmn,31.73);

    % fix for negative uncertainty
    if uncn == 0; dvmn(deglidxunc:end) = 0; end;
    depthn = flip(cummin(flip(dvm-dvmn'))); % cummin to avoid wobbling uncertainty

    % fix for unconstrained erosion with constant depth at 1 km
    if uncp == max(Etestv); dvmp(deglidxunc:end) = 1E3; end;
    depthp = cummax(dvm+dvmp'); % cummax to avoid wobbling uncertainty

    % if saving time-depth data
    if numel(tdv) > 0;
        % fix time-depth strings and save in output
        for j = 1:numel(tdv);
            tdstr(j*3-2) = {num2str(interp1(pars.tv,dvm,tdv(j)),'%.2f')};
            tdstr(j*3-1) = {num2str(interp1(pars.tv,depthp-dvm,tdv(j)),'%.2f')};
            tdstr(j*3) = {num2str(interp1(pars.tv,dvm-depthn,tdv(j)),'%.2f')};
        end;
        output(i+1,outcol.(['E' nstr ri gd '1']):outcol.(['E' nstr ri gd '2'])) = tdstr;
    end;
% end subfunction tdv_fix ==========================================================================


% subfunction pdm_fix ==============================================================================
function inpdm = pdm_fix(inpdm,mm,tv,plch,ymax);
    % calculate t and m step
    maxt = plch.maxt .* 1E3; % convert to yr
    xstep = maxt/(plch.pdmX-1);
    ystep = ymax/(plch.pdmY-1);
    % interpret mm to plotting time
    tv_plot = (0:xstep:maxt)';
    mm_plot = interp1(tv,mm,tv_plot);

    % fix time and depth vectors
    tvv = repmat(tv_plot,size(mm_plot,2),1);
    mmv = mm_plot(:);

    % fix first pdm
    [pdm1,centers] = hist3([mmv,tvv],'ctrs',{[(0:ystep:ymax) ymax+ystep],[(0:xstep:maxt)]});

    % remove last row (number of points outside domain)
    inpdm.pdm = inpdm.pdm + pdm1(1:end-1,1:end);
    inpdm.xv = (0:xstep*1E-3:maxt*1E-3);
    inpdm.yv = (0:ystep:ymax);
    %[inpdm.xm,inpdm.ym] = meshgrid((0:xstep*1E-3:maxt*1E-3),(0:ystep:ymax));
% end subfunction pdm_fix ==========================================================================


% subfunction plfix ================================================================================
function pl = plfix(plch,sample,nucl1026,glacErate,glacEstep,sensstr);
    pl = struct;
    num = 0; % figure number
    % fix for dPN figures
    dPN = {}; leg = {};
    if plch.depth == 1; dPN(end+1) = {'d'}; leg(end+1) = {'Sample depth (m)'}; end;
    if plch.P == 1;
        dPN(end+1) = {'P'};
        if plch.Pprcnt == 1; leg(end+1) = {'P (%)'}; else; leg(end+1) = {'P (atoms/g/yr)'}; end;
    end;
    if plch.N == 1;
        dPN(end+1) = {'N'};
        if plch.Nprcnt == 1; leg(end+1) = {'N (%)'}; else; leg(end+1) = {'N (atoms/g)'}; end;
    end;
    for j = 1:numel(dPN);
        if sum(sample.N10) > 0;
            if glacErate == 1; [pl,num] = plfixdPN(pl,num,plch,'r','10',dPN{j},leg{j}); end;
            if glacEstep == 1; [pl,num] = plfixdPN(pl,num,plch,'i','10',dPN{j},leg{j}); end;
        end;
        if sum(sample.N26) > 0;
            if glacErate == 1; [pl,num] = plfixdPN(pl,num,plch,'r','26',dPN{j},leg{j}); end;
            if glacEstep == 1; [pl,num] = plfixdPN(pl,num,plch,'i','26',dPN{j},leg{j}); end;
        end;
        if sum(sample.N14) > 0;
            if glacErate == 1; [pl,num] = plfixdPN(pl,num,plch,'r','14',dPN{j},leg{j}); end;
            if glacEstep == 1; [pl,num] = plfixdPN(pl,num,plch,'i','14',dPN{j},leg{j}); end;
        end;
        if strcmp(dPN{j},'d') && nucl1026==1 && sum(sample.N10.*sample.N26)>0;
            if glacErate == 1; [pl,num] = plfixdPN(pl,num,plch,'r','1026','d',leg{j}); end;
            if glacEstep == 1; [pl,num] = plfixdPN(pl,num,plch,'i','1026','d',leg{j}); end;
        end;
    end;
    % fix for sensitivity figures
    if numel(sensstr) > 0;
        ratestepv = [glacErate glacEstep glacErate glacEstep glacErate glacEstep];
        riv = {'r','i','r','i','r','i'};
        nuclv = {'10','10','26','26','14','14'};
        for j = 1:6;
            if sum(sample.(['N' nuclv{j}]))>0 && ratestepv(j)==1;
                pl.(['sens' nuclv{j} riv{j}]).pos = [];
                pl.(['sens' nuclv{j} riv{j}]).neg = [];
                num = num + 1;
                pl.(['sens' nuclv{j} riv{j}]).num = num;
            end;
        end;
    end;
    % fix for banana figure
    if plch.banana==1 && nucl1026==1 && sum(sample.N10.*sample.N26)>0;
        num = num + 1;
        pl.banana.num = num;
        pl.banana.P10z = []; pl.banana.P26z = [];
        if glacErate == 1;
            pl.banana.Npath10r = []; pl.banana.Npath26r = [];
            pl.banana.Npath10r_out = []; pl.banana.Npath26r_out = [];
        end;
        if glacEstep == 1;
            pl.banana.Npath10i = []; pl.banana.Npath26i = [];
            pl.banana.Npath10i_out = []; pl.banana.Npath26i_out = [];
        end;
        pl.banana.N10 = []; pl.banana.N26 = [];
        pl.banana.Nuncx = []; pl.banana.Nuncy = []; pl.banana.Nuncm = [];
        pl.banana.N10_out = []; pl.banana.N26_out = [];
        pl.banana.Nuncx_out = []; pl.banana.Nuncy_out = []; pl.banana.Nuncm_out = [];
    end;
% end subfunction plfix ============================================================================


% subfunction plfixdPN =============================================================================
function [pl,num] = plfixdPN(pl,num,plch,ri,nucl,dPN,leg);
    % fix plot number
    if (plch.combined==1 && isfield(pl,[ri '10' dPN]) && strcmp(nucl,'1026')==0) || ...
            (plch.combined_full==1 && isfield(pl,[ri '10' dPN]));
        nn = pl.([ri '10' dPN]).num;
    else;
        num = num + 1; nn = num;
    end;
    pl.([ri nucl dPN]).num = nn;
    pl.([ri nucl dPN]).t = [];
    pl.([ri nucl dPN]).tunc = [];
    pl.([ri nucl dPN]).m = [];
    pl.([ri nucl dPN]).munc = [];
    pl.([ri nucl dPN]).y = [];
    pl.([ri nucl dPN]).yunc = [];
    pl.([ri nucl dPN]).emax = [];
    pl.([ri nucl dPN]).clr = plch.(['clr' nucl]);
    pl.([ri nucl dPN]).yleg = leg;
    pl.([ri nucl dPN]).pdm = zeros(plch.pdmY,plch.pdmX);
    pl.([ri nucl dPN]).prcntm = [];
    pl.([ri nucl dPN]).prcnttv = (0:plch.pct_dt:plch.maxt.*1E3);
    % fix tm, tmunc, mm, mmunc
    pl.tm = []; pl.tmunc = []; pl.mm = []; pl.mmunc = [];
% end subfunction plfixdPN =========================================================================


% subfunction fixoutput1 ===========================================================================
function [output,outcol,sensstr] = fixoutput1(output,outcol,sample,glacErate,glacEstep,glaccalc,...
    uncsens,tdv,gdv,nucl1026);
% function for fixing output header and outcol numbers
    % fix sensstr for sensitivity output ==================================
    sensstr = {};
    if uncsens.N==1; sensstr(end+1) = {'N'}; end;
    if uncsens.Pref==1; sensstr(end+1) = {'Pref'}; end;
    if uncsens.l==1; sensstr(end+1) = {'l'}; end;
    if uncsens.dens==1 && (isfield(sample,'densunc') || ...
        (isfield(sample,'densmin') && isfield(sample,'densmax'))); sensstr(end+1) = {'dens'}; end;
    if uncsens.erosion==1 && (isfield(sample,'erosionunc') || ...
        (isfield(sample,'erosionmin') && isfield(sample,'erosionmax'))) && glaccalc==1;
        sensstr(end+1) = {'erosion'};
    end;
    if uncsens.glacE==1 && (isfield(sample,'glacEunc') || ...
        (isfield(sample,'glacEmin') && isfield(sample,'glacEmax'))) && glaccalc==0;
        sensstr(end+1) = {'glacE'};
    end;
    parstr = {'simt','icevalue','deglac','isostsubm'};
    for j = 1:numel(parstr);
        if uncsens.(parstr{j})==1 && (isfield(sample,[parstr{j} 'unc']) || ...
            (isfield(sample,[parstr{j} 'min']) && isfield(sample,[parstr{j} 'max'])));
            sensstr(end+1) = parstr(j);
        end;
    end;
    % fix vectors for for loop ============================================
    if glacErate == 1; r = 'r'; else; r = '0'; end;
    if glacEstep == 1; i = 'i'; else; i = '0'; end;
    riv = {r,i,r,i,r,i,r,i};
    test10 = 0; test26 = 0; test14 = 0; test1026 = 0;
    if sum(sample.N10) > 0; test10 = 1; end;
    if sum(sample.N26) > 0; test26 = 1; end;
    if sum(sample.N14) > 0; test14 = 1; end;
    if sum(sample.N10.*sample.N26)>0 && nucl1026==1; test1026 = 1; end;
    testNv = [test10 test10 test26 test26 test14 test14 test1026 test1026];
    nstrv = {'10','10','26','26','14','14','1026','1026'};
    % for loop for nuclide and rate/step case =============================
    for j = 1:8;
        % check for input
        if strcmp(riv{j},'0') || testNv(j) == 0; continue; end;
        % fix ri and nstr
        ri = riv{j}; nstr = nstrv{j};
        % fix unit string
        if strcmp(riv{j},'i') && glaccalc==1; unit = '(cm/glac)'; else; unit = '(mm/ka)'; end;
        % fix i for nonglacial erosion calculation with incremental glacial erosion
        if strcmp(riv{j},'i') && glaccalc==0; ii = 'i'; else; ii = ''; end;
        % fix output and outcol for E and unc
        output(1,end+1:end+3) = {['E' nstr ii unit],['E' nstr ii '+' unit],['E' nstr ii '-' unit]};
        outcol.(['E' nstr ri]) = max(cell2mat(struct2cell(outcol)))+1;
        outcol.(['E' nstr ri 'p']) = outcol.(['E' nstr ri])+1;
        outcol.(['E' nstr ri 'n']) = outcol.(['E' nstr ri 'p'])+1;
        % fix for 1026 overlap header
        if strcmp(nstr,'1026');
            output(1,end+1) = {['E1026' ii unit 'overlap']};
            outcol.(['E1026' ri 'overlap']) = max(cell2mat(struct2cell(outcol)))+1;
        end;
        % fix uncertainty sensitivity output
        if numel(sensstr)>0 && strcmp(nstr,'1026')==0;
            for i = 1:numel(sensstr);
                output(1,end+1:end+2) = {[sensstr{i} '+' unit],[sensstr{i} '-' unit]};
            end;
            outcol.(['E' nstr ri 'sens1']) = max(cell2mat(struct2cell(outcol)))+1;
            outcol.(['E' nstr ri 'sens2']) = max(cell2mat(struct2cell(outcol)))+numel(sensstr)*2-1;
        end;
        % fix for time depth points headers
        for k = 1:numel(tdv);
            if tdv(k) < 1E3;
                tstr = [num2str(tdv(k),'%.0f') 'yr'];
            elseif tdv(k) < 1E6;
                tstr = [num2str(tdv(k).*1E-3,'%.0f') 'ka'];
            elseif tdv(k).*1E-6 == floor(tdv(k).*1E-6);
                tstr = [num2str(tdv(k).*1E-6,'%.0f') 'Ma'];
            else;
                tstr = [num2str(tdv(k).*1E-6,'%.1f') 'Ma'];
            end;
            tdheader(k*3-2:k*3) = {[tstr '(m)'],[tstr '+(m)'],[tstr '-(m)']};
        end;
        if numel(tdv)>0;
            output(1,end+1:end+numel(tdheader)) = tdheader;
            outcol.(['E' nstr ri 'd1']) = max(cell2mat(struct2cell(outcol)))+1;
            outcol.(['E' nstr ri 'd2']) = max(cell2mat(struct2cell(outcol)))+numel(tdheader)-1;
        end;
        % fix for time glacial depth points headers
        for k = 1:numel(gdv);
            if gdv(k) < 1E3;
                tstr = ['glacE' num2str(gdv(k),'%.0f') 'yr'];
            elseif gdv(k) < 1E6;
                tstr = ['glacE' num2str(gdv(k).*1E-3,'%.0f') 'ka'];
            elseif gdv(k).*1E-6 == floor(gdv(k).*1E-6);
                tstr = ['glacE' num2str(gdv(k).*1E-6,'%.0f') 'Ma'];
            else;
                tstr = ['glacE' num2str(gdv(k).*1E-6,'%.1f') 'Ma'];
            end;
            gdheader(k*3-2:k*3) = {[tstr '(m)'],[tstr '+(m)'],[tstr '-(m)']};
        end;
        if numel(gdv)>0;
            output(1,end+1:end+numel(gdheader)) = gdheader;
            outcol.(['E' nstr ri 'gd1']) = max(cell2mat(struct2cell(outcol)))+1;
            outcol.(['E' nstr ri 'gd2']) = max(cell2mat(struct2cell(outcol)))+numel(gdheader)-1;
        end;
    end;
    % fix for sensstr
    if numel(sensstr)>0 && uncsens.full==1; % fix for plotting full uncertainty
        sensstr(end+1) = {'full'};
    end;
% end subfunction fixoutput1 =======================================================================


% subfunction fixoutplot ===========================================================================
function [output,pl] = fixoutplot(output,pl,sample,NE,unc,nucl,ri,glaccalc,i,outcol,tdv,gdv,plch,...
    pars,parsu,prp1,prpunc,Etestv,np);
% display erosion and fix for output and plotting
    % fix nuclide strings
    if nucl == 10; nstr = '10'; nstr1 = '10Be';
    elseif nucl == 26; nstr = '26'; nstr1 = '26Al';
    elseif nucl == 14; nstr = '14'; nstr1 = '14C';
    elseif nucl == 1026; nstr = '1026'; nstr1 = '10Be+26Al'; end;
    % fix unit string
    if glaccalc==1 && strcmp(ri,'i'); unit = 'cm/glac';
    else; unit = 'mm/ka'; end;
    % display erosion
    fprintf(1,['\n' nstr1 ': %.2f +/- [%.2f / %.2f] '],NE.(ri),unc.([ri 'pos']),unc.([ri 'neg']));
    fprintf(1,unit);
    output(i+1,outcol.(['E' nstr ri])) = {num2str(NE.(ri),'%.2f')};
    output(i+1,outcol.(['E' nstr ri 'p'])) = {num2str(unc.([ri 'pos']),'%.2f')};
    output(i+1,outcol.(['E' nstr ri 'n'])) = {num2str(unc.([ri 'neg']),'%.2f')};
    % fix for plotting and depth at time points
    if numel(tdv) > 0 || numel(gdv) > 0 || sum([plch.depth plch.P plch.N]) > 0;
        [output,pl] = depthPconc(output,pl,pars,parsu,prp1,prpunc,NE.(ri),unc.([ri 'pos']),...
            unc.([ri 'neg']),unc.(['Ev' ri]),Etestv,np,plch,tdv,gdv,glaccalc,i,outcol,ri,nstr);
    end;
% end subfunction fixoutplot =======================================================================


% subfunction plot_dPN =============================================================================
function pl = plot_dPN(pl,plch,dPN,nucl,ri);
    if isfield(pl,[ri nucl dPN]) == 0; return; end;
    % fix for figure name
    figname = [num2str(pl.([ri nucl dPN]).num) ': '];
    if strcmp(ri,'r'); ratestep = 'Erate'; else; ratestep = 'Estep'; end;
    if strcmp(nucl,'10') == 0 && isfield(pl,[ri '10' dPN]) && ...
            pl.([ri nucl dPN]).num == pl.([ri '10' dPN]).num;
        figname = [figname ratestep '-' dPN];
    else;
        figname = [figname ratestep nucl dPN];
    end;
    % fix for legend
    if strcmp(nucl,'10'); nstr = '^{10}Be';
    elseif strcmp(nucl,'26'); nstr = '^{26}Al';
    elseif strcmp(nucl,'14'); nstr = '^{14}C';
    elseif strcmp(nucl,'1026'); nstr = '^{10}Be + ^{26}Al'; end;
    figure(pl.([ri nucl dPN]).num,'name',figname,'NumberTitle','off'); hold on; box on;
    % plot probability density map
    if plch.probdensmap == 1;
        % find max value in last column of pdm and use for scaling the colors
        lastcolmax = max(pl.([ri nucl dPN]).pdm(:,end));
        levels = linspace(0,lastcolmax*3,200);
        % make meshgrids
        [xm,ym] = meshgrid(pl.([ri nucl dPN]).xv,pl.([ri nucl dPN]).yv);
        % draw probability density map with contourf
        contourf(xm,ym,pl.([ri nucl dPN]).pdm,levels,'linestyle','none');
        % fix colors
        pdmclr = customcolormap([0 0.5 0.75 0.9 1],[0.5 0 0;1 0 0;1 0.9 0;0.6 0.7 1;1 1 1],200);
        colormap(pdmclr);
        plch.(['clr' nucl]) = 'black'; % use black for mid lines
        % save txt file
        if plch.pdm_save == 1;
            fid = fopen(['glacE-' ri nucl dPN '-pdm.txt'], 'w');
            for i=1:size(pl.([ri nucl dPN]).pdm,1);
                fprintf(fid,'%f ',pl.([ri nucl dPN]).pdm(i,:));
                fprintf(fid,'\n');
            end;
            fclose(fid);
        end;
    end;
    % plot full percentiles
    if plch.pct_area == 1;
        % fix colors
        values = [0 0.1 0.35 0.5 0.65 0.9 1]';
        clrv = [1 1 0.95;1 1 0;1 0 0;0.3 0 0;1 0 0;1 1 0;1 1 0.95];
        step = 1/plch.pct_areaN;
        values2 = (step/2:step:1-step/2);
        prcntclr = interp1(values,clrv,values2);
        % interpolate dPN
        prcntm = prctile(pl.([ri nucl dPN]).prcntm,linspace(0,100,plch.pct_areaN+1));
        prcntmf = flip(prcntm')';
        prcttv = [pl.([ri nucl dPN]).prcnttv,flip(pl.([ri nucl dPN]).prcnttv)] .* 1E-3;
        % plot each percentile interval
        for i=1:plch.pct_areaN;
            patch('XData',prcttv,'YData',[prcntm(i,:),prcntmf(i+1,:)],'FaceColor',prcntclr(i,:),...
                'EdgeColor','none','FaceAlpha',0.7);
        end;
        plch.(['clr' nucl]) = 'black'; % use black for mid lines
        if plch.colorbar == 1;
            colormap(prcntclr);
            colorbar();
        end;
    end;
    % fix for samples with varying tv
    pl.([ri nucl dPN]).y(pl.([ri nucl dPN]).m==0) = NaN;
    pl.([ri nucl dPN]).yunc(pl.([ri nucl dPN]).munc==0) = NaN;
    pl.([ri nucl dPN]).t(pl.([ri nucl dPN]).m==0) = NaN;
    pl.([ri nucl dPN]).tunc(pl.([ri nucl dPN]).munc==0) = NaN;
    % plot mid line
    [t,y] = cutplot(pl.([ri nucl dPN]).t,pl.([ri nucl dPN]).y,plch,dPN);
    td = []; yd = [];
    if strcmp(dPN,'d') && plch.splitdd > 0;
        % split depth lines at plch.splitdd for samples with maximum possible erosion,
        % deeper part plotted dashed (other samples are plotted as full solid lines)
        [ts,ys] = cut_half(t,y,plch.splitdd,2,1);   % shallow part
        [td,yd] = cut_half(t,y,plch.splitdd,2,0);   % deeper part
        notmax = find(~pl.([ri nucl dPN]).emax);
        if numel(notmax) > 0;
            n = size(t,1);
            ts(1:n,notmax) = t(:,notmax); ts(n+1:end,notmax) = NaN;
            ys(1:n,notmax) = y(:,notmax); ys(n+1:end,notmax) = NaN;
            td(:,notmax) = NaN; yd(:,notmax) = NaN;
        end;
        t = ts; y = ys;
    end;
    leg = plot(t,y,'color',plch.(['clr' nucl]),'linewidth',0.1);
    if any(~isnan(td(:)));
        plot(td,yd,'color',plch.(['clr' nucl]),'linewidth',0.1,'linestyle','--');
    end;
    pl.([ri nucl dPN]).leg = leg(1); % fix for legend
    pl.([ri nucl dPN]).legin = {nstr};
    legend(pl.([ri nucl dPN]).leg,pl.([ri nucl dPN]).legin,'location','northwest','AutoUpdate',...
        'off');
    if (strcmp(nucl,'26') || strcmp(nucl,'1026') || strcmp(nucl,'14')) && isfield(pl,[ri '10' dPN]);
        if pl.([ri nucl dPN]).num == pl.([ri '10' dPN]).num;
            pl.([ri '10' dPN]).leg(end+1) = pl.([ri nucl dPN]).leg;
            pl.([ri '10' dPN]).legin(end+1) = pl.([ri nucl dPN]).legin;
            legend(pl.([ri '10' dPN]).leg,pl.([ri '10' dPN]).legin,'location','northwest',...
                'AutoUpdate','off');
        end;
    end;
    % plot uncertainties
    if plch.uncline == 1;
        [tunc,yunc] = cutplot(pl.([ri nucl dPN]).tunc,pl.([ri nucl dPN]).yunc,plch,dPN);
        plot(tunc,yunc,'color',...
            min([plch.(['clr' nucl])+0.7;1 1 1]));
    end;
    % cut uncertainties at maxd
    if strcmp(dPN,'d') && plch.cutuncmaxd == 1;
        pl.([ri nucl dPN]).yunc(pl.([ri nucl dPN]).yunc > plch.maxd) = plch.maxd;
    end;
    if plch.uncfill == 1;
        % remove yunc nan points
        notnan_idx = ~isnan(pl.([ri nucl dPN]).yunc);
        pl.([ri nucl dPN]).tunc = pl.([ri nucl dPN]).tunc(notnan_idx);
        pl.([ri nucl dPN]).yunc = pl.([ri nucl dPN]).yunc(notnan_idx);
        patch('XData',pl.([ri nucl dPN]).tunc,'YData',pl.([ri nucl dPN]).yunc,'FaceColor',...
            plch.(['clr' nucl]),'EdgeColor','none','FaceAlpha',0.1);
    end;
    % plot percentile lines
    if numel(plch.pct_lines) > 0;
        % interpolate dPN
        prcntm = prctile(pl.([ri nucl dPN]).prcntm,plch.pct_lines);
        prcttv = pl.([ri nucl dPN]).prcnttv .* 1E-3;
        % plot percentile lines
        for i=1:numel(plch.pct_lines);
            [tp,yp] = cutplot(prcttv,prcntm(i,:),plch,dPN);
            plot(tp,yp,'color',[1 0 1],'linewidth',0.3);
        end;
    end;
    set(gca,'layer','top'); % plot axis on top
    set(gca,'TickDir','out');
    xlim([0 plch.maxt]); % fix x-axis
    if strcmp(dPN,'d'); ylim([0 plch.maxd]); axis('ij'); end; % fix reversed y-axis for depth plot
    xlim([0 plch.maxt]); set(gca(),'xdir','reverse'); % fix x-axis
    xlabel('Time (ka)'); ylabel(pl.([ri nucl dPN]).yleg); % fix axis labels
    hold off;
% end subfunction plot_dPN =========================================================================


% subfunction cutplot =============================================================================
function [t,y] = cutplot(t,y,plch,dPN);
% cut plotted lines at plch.maxt and, for depth, at plch.maxd (if plch.cutmaxxy == 1)
    if plch.cutmaxxy ~= 1; return; end;
    [t,y] = cut_half(t,y,plch.maxt,1,1);
    if strcmp(dPN,'d'); [t,y] = cut_half(t,y,plch.maxd,2,1); end;
% end subfunction cutplot =========================================================================


% subfunction cut_half ============================================================================
function [t,y] = cut_half(t,y,bnd,ax,keep);
% cut lines at x = bnd (ax = 1) or y = bnd (ax = 2), keeping points with values <= bnd
% (keep == 1) or >= bnd (keep == 0) and setting all other points to NaN. The exact crossing
% point with the cut line is inserted, so the lines are cut exactly at bnd.
    if ax == 1; v = t; else; v = y; end;
    rowvec = size(t,1) == 1 && size(t,2) > 1; % fix row vector input
    if rowvec; t = t'; y = y'; v = v'; end;
    [n,nc] = size(t);
    tn = nan(2*n,nc); yn = nan(2*n,nc); % fix output (crossing points are added)
    for j = 1:nc;
        vj = v(:,j); tj = t(:,j); yj = y(:,j);
        fin = ~isnan(tj) & ~isnan(yj);
        if keep == 1; in = vj <= bnd & fin; else; in = vj >= bnd & fin; end;
        k = 0;
        for i = 1:n;
            if i > 1 && fin(i) && fin(i-1) && in(i) ~= in(i-1);
                % insert exact crossing point with the cut line
                f = (bnd - vj(i-1)) ./ (vj(i) - vj(i-1));
                k = k + 1;
                tn(k,j) = tj(i-1) + f .* (tj(i) - tj(i-1));
                yn(k,j) = yj(i-1) + f .* (yj(i) - yj(i-1));
            end;
            if in(i);
                if i > 1 && ~fin(i-1);
                    % keep breaks caused by NaN points in the input
                    k = k + 1;
                end;
                k = k + 1;
                tn(k,j) = tj(i);
                yn(k,j) = yj(i);
            end;
        end;
    end;
    t = tn; y = yn;
    if rowvec; t = t'; y = y'; end;
% end subfunction cut_half =======================================================================


% subfunction plot_uncsens =========================================================================
function plot_uncsens(plsens,sensstr,sensclr,glaccalc,ri,nucl);
    figname = [num2str(length(findobj('type','figure')) + 1) ': '];
    if strcmp(ri,'r'); ratestep = 'Erate'; else; ratestep = 'Estep'; end;
    if nucl == 10; nuclstr = '^{10}Be'; figname = [figname ratestep '10'];
    elseif nucl == 26; nuclstr = '^{26}Al'; figname = [figname ratestep '26'];
    elseif nucl == 14; nuclstr = '^{14}C'; figname = [figname ratestep '14']; end;
    figure('name',figname,'NumberTitle','off'); hold on; box on;
    xm = [(1:1:size(plsens.pos,1));(1:1:size(plsens.pos,1))];
    for i = 1:numel(sensstr);
        legtemp = plot(xm,[plsens.pos(:,i)';plsens.neg(:,i)'],'color',sensclr.(sensstr{i}));
        leg(i) = legtemp(1);
        xm = xm + size(plsens.pos,1) + 1;
    end;
    legend(leg,sensstr,'AutoUpdate','off');
    xlim([0 size(plsens.pos,1)*numel(sensstr)+numel(sensstr)]);
    set(gca,'xtick',[]); set(gca,'xticklabel',[]);
    if glaccalc == 1; glstr = 'glacial '; else; glstr = ''; end;
    ylabel(['Relative ' glstr 'erosion uncertainties']);
    if glaccalc==1 && strcmp(ri,'i');
        etypestr = ' incremental erosion';
    else;
        etypestr = ' erosion rate';
    end;
    title([nuclstr etypestr]);
    hold off;
% end subfunction plot_uncsens =====================================================================


% subfunction bananaN ==============================================================================
function pl = bananaN(pl,E1026,np10,np26,P10u,P26u);
    % fix for samples with overlapping 1026 E and samples with no 1026 E
    if (isfield(E1026,'r') && E1026.overlapr == 1) || (isfield(E1026,'i') && E1026.overlapi == 1);
        str = '';
    else;
        str = '_out';
    end;
    % fill pl with normalized N
    pl.banana.(['N10' str])(end+1,1) = np10.N/P10u;
    pl.banana.(['N26' str])(end+1,1) = np26.N/P26u;
    % add Pref uncertainty to N uncertainty
    N10unc = sqrt(np10.Nunc^2 + (np10.N*np10.Prefunc/np10.Pref)^2);
    N26unc = sqrt(np26.Nunc^2 + (np26.N*np26.Prefunc/np26.Pref)^2);
    % calculate uncertainty line for banana plot
    [xx,yy] = bananaNunc(np10.N/P10u,np26.N/P26u,N10unc/P10u,N26unc/P26u);
    % fill pl with x and y values
    pl.banana.(['Nuncx' str])(1:numel(xx),end+1) = xx;
    pl.banana.(['Nuncy' str])(1:numel(yy),end+1) = yy;
    pl.banana.(['Nuncm' str])(1:numel(xx),end+1) = 1;
% end subfunction bananaN ==========================================================================


% subfunction bananaNunc ===========================================================================
function [outx,outy] = bananaNunc(N10n,N26n,N10uncn,N26uncn);
    % estimate range and create mesh
    R = (N26n/N10n);
    delR = sqrt((N26uncn/N26n)^2 + (N10uncn/N10n)^2);
    [x,y] = meshgrid((N10n-4*N10uncn):(0.1*N10uncn):(N10n+4*N10uncn),...
        (R*(1-4*delR)):(0.1*R*delR):(R*(1+4*delR)));
    % calculate PDF
    Prob = x.*exp(-0.5.*((((y.*x) - N26n)./N26uncn).^2 + ((x - N10n)./N10uncn).^2));
    % Normalize to volume = 1
    normP = Prob ./ sum(sum(Prob));
    % Now we need to figure out cumulative probabilities.
    % multiply by 10000 to achieve manageable number of values, round to get integers:
    normP = normP * 10000;
    intP = round(normP);
    for a = 1:max(max(intP));
        cumprob(a) = sum(intP(find(intP >= a)))/10000;
    end;
    probs = 1:max(max(intP));
    sigma1 = find(abs(cumprob - 0.6827) == min(abs(cumprob - 0.6827)));
    sigma2 = find(abs(cumprob - 0.9545) == min(abs(cumprob - 0.9545)));
    % weed out cases where adjacent probs are same -- rounding error
    if length(sigma1) ~= 1; sigma1 = min(sigma1); end;
    if length(sigma2) ~= 1; sigma2 = min(sigma2); end;
    % Now draw the contours.
    cmat = contourc(x(1,:),y(:,1),normP,[sigma1 sigma2]);
    % Sometimes contourc returns several contours for one level - grid size issue?
    % This is spurious, so plot only the major one.
    cntStarts = find(cmat(1,:) == sigma1);
    cntSizes = cmat(2,cntStarts);
    cntToPlot = find(cntSizes == max(cntSizes));
    outx = cmat(1,(cntStarts(cntToPlot)+1):(cntStarts(cntToPlot) + cntSizes(cntToPlot)))';
    outy = cmat(2,(cntStarts(cntToPlot)+1):(cntStarts(cntToPlot) + cntSizes(cntToPlot)))';
% end subfunction bananaNunc =======================================================================


% subfunction Npathcalc ============================================================================
function Nm = Npathcalc(tv,P,l);
    % fix flipped full prod matrices and tv matrix
    Pfl = flip(P);
    tvfl = flip(tv);
    % fix output N matric
    Nm(1,size(Pfl,2)) = 0;
    % calculate N through time
    for i = 1:numel(tvfl)-1;
        Nm(i+1,:) = (Nm(end,:) + Pfl(i,:).*(tvfl(i)-tvfl(i+1))) .* exp(-(tvfl(i)-tvfl(i+1)).*l);
    end;
    Nm = flip(Nm);
% end subfunction Npathcalc ========================================================================


% subfunction plot_banana ==========================================================================
function plot_banana(pl,consts,z);
    figname = [num2str(length(findobj('type','figure')) + 1) ': 26/10-banana'];
    figure('name',figname,'NumberTitle','off'); hold on; box on;
    %figure(pl.num); hold on; box on;
    % create data for the simple exposure line
    tempt = [1 (100:100:900) logspace(3,7,60)]';
    be = (1/consts.l10)*(1-exp(-consts.l10*tempt));
    al = (1/consts.l26)*(1-exp(-consts.l26*tempt));
    % make data for simple burial lines
    burv = [5E5 1E6 1.5E6 2E6 2.5E6 3E6];
    bu_be = []; bu_al = [];
    for j = 1:numel(burv)
        bu_be(:,end+1) = be.*exp(-consts.l10*burv(j));
        bu_al(:,end+1) = al.*exp(-consts.l26*burv(j));
    end;
    % plot simple burial lines
    semilogx(bu_be,bu_al./bu_be,'linestyle','--','color','black');
    % make data for burial pathways
    pathv = [1E4 1E5 1E6 1E7];
    burv = (0:1E5:1E7)';
    path_be = []; path_al = [];
    for j = 1:numel(pathv);
        path_be(:,end+1) = (1/consts.l10)*(1-exp(-consts.l10*pathv(j))).*exp(-consts.l10.*burv);
        path_al(:,end+1) = (1/consts.l26)*(1-exp(-consts.l26*pathv(j))).*exp(-consts.l26.*burv);
    end;
    % plot burial pathways
    semilogx(path_be,path_al./path_be,'linestyle','--','color','black');
    % fix data for erosion line
    P10z = mean(pl.P10z)'; P26z = mean(pl.P26z)';
    tempe = logspace(-5,1,100);
    tvm = repmat(z',size(tempe))./repmat(tempe,size(z'));
    dcf10 = exp(-tvm.*consts.l10); dcf26 = exp(-tvm.*consts.l26);
    N10tv = trapz_m(tvm,repmat(P10z,size(tempe)).*dcf10);
    N26tv = trapz_m(tvm,repmat(P26z,size(tempe)).*dcf26);
    bee = N10tv./P10z(1); ale = N26tv./P26z(1);
    bee(1,2:end+1) = bee; ale(1,2:end+1) = ale;
    bee(1,1) = be(end); ale(1,1) = al(end);
    % plot erosion end-point line
    semilogx(bee,ale./bee,'color','black');
    % plot simple exposure line
    semilogx(be,al./be,'color','black');
    % plot colors
    Nin = [1 0.35 0]; Nout = [0.7 0.7 0.7]; pathr = [0.66 0.37 0]; pathi = [1 0.7 0.15];
    % plot sample N uncertainty lines for samples without 1026 overlap
    if numel(pl.Nuncx_out) > 0;
        pl.Nuncx_out(pl.Nuncm_out==0) = NaN; pl.Nuncy_out(pl.Nuncm_out==0) = NaN;
        semilogx(pl.Nuncx_out,pl.Nuncy_out,'color',Nout);
    end;
    % plot sample N uncertainty lines for samples with 1026 overlap
    if numel(pl.Nuncx) > 0;
        pl.Nuncx(pl.Nuncm==0) = NaN; pl.Nuncy(pl.Nuncm==0) = NaN;
        semilogx(pl.Nuncx,pl.Nuncy,'color',Nin);
    end;
    % plot sample N
    if numel(pl.N10_out) > 0;
        semilogx(pl.N10_out,pl.N26_out./pl.N10_out,'o','markersize',4,'color',Nout);
    end;
    if numel(pl.N10) > 0;
        semilogx(pl.N10,pl.N26./pl.N10,'o','markersize',4,'color',Nin);
    end;
    % plot sample N pathways
    leg = []; legin = {};
    if isfield(pl,'Npath10r') && numel(pl.Npath10r)>0;
        if isfield(pl,'Npath10i')==0; pathr = Nin; end;
        pl.Npath10r(pl.Npath10r<=0) = NaN; pl.Npath26r(pl.Npath26r<=0) = NaN;
        leg1 = semilogx(pl.Npath10r,pl.Npath26r./pl.Npath10r,'color',pathr);
        leg(end+1) = leg1(1);
        legin(end+1) = 'Erate';
    end;
    if isfield(pl,'Npath10i') && numel(pl.Npath10i)>0;
        if isfield(pl,'Npath10r')==0; pathi = Nin; end;
        pl.Npath10i(pl.Npath10i<=0) = NaN; pl.Npath26i(pl.Npath26i<=0) = NaN;
        leg1 = semilogx(pl.Npath10i,pl.Npath26i./pl.Npath10i,'color',pathi);
        leg(end+1) = leg1(1);
        legin(end+1) = 'Estep';
    end;
    axis([1E3 3E6 0.2 1.2]);
    xlabel('[^{10}Be]*');
    ylabel('[^{26}Al]*/[^{10}Be]*');
    legend(leg,legin,'location','northeast','AutoUpdate','off');
    set(gca,'layer','top'); % plot axis on top
    set(gca,'XScale','log'); % fix for matlab
    set(gca,'TickDir','out');
    hold off;
% end subfunction plot_banana ======================================================================


% subfunction trapz_m ==============================================================================
function z = trapz_m(x,y);
% trapz function for matrices
    dim = 1;
    nd = ndims(y);
    sz = size(y);
    n = sz(dim);
    idx1 = repmat ({':'}, [nd, 1]);
    idx2 = idx1;
    idx1{dim} = 2 : n;
    idx2{dim} = 1 : (n - 1);
    z = 0.5 * sum (diff (x, 1, dim) .* (y(idx1{:}) + y(idx2{:})), dim);
% end subfunction trapz_m ==========================================================================


% subfunction get_waterdepth =======================================================================
function [waterdm,delvmfull] = get_waterdepth(elv,delv,icem,tv);
    % fix ice matrix
    if size(delv,2)>1 && size(icem,2)==1;
        icem = repmat(icem,1,size(delv,2));
    end;
    
    % fix delv matrix
    if size(icem,2)>1 && size(delv,2)==1;
        delv = repmat(delv,1,size(icem,2));
    end;
    
    % fix noice matrix
    noicem = (icem == 0);
    
    % fix postglacial delv matrix
    delvm = delv .* (cumsum(icem)==0);
    delvm_min = min(delvm); % max elevation change (assumed to be at deglaciation)
    
    % flipped ice/noice matrices and tv
    flipicem = flip(icem);
    flipnoicem = flip(noicem);
    flipdtv = flip([diff(tv);0]);
    % uplift matrix
    uplm(1,size(icem,2)) = 0;
    % uplift constant
    uplconst = 0.0002;
    % max uplift vector
    maxuplv = repmat(30000,1,size(icem,2));
    
    % calculate uplift through simulation time
    for i = 2:numel(flipdtv);
        uplm(i,:) = min([uplm(i-1,:)+flipicem(i,:).*flipdtv(i);maxuplv]) .* ...
            exp(-uplconst.*flipdtv(i).*flipnoicem(i,:));
    end;
    
    % flip back
    uplm = flip(uplm);
    
    % calculated deglac and present-day uplift
    uplm_degl = repmat(max(uplm.*(cumsum(icem)==0)),size(tv));
    uplm_present = repmat(uplm(1,:),size(tv));
    
    % scale uplm based on delvm_min to get simulated delv
    delvmsim = (uplm-uplm_present)./uplm_degl.*delvm_min.*(cumsum(icem)>0);
    
    % combine postglacial and pre-LGM uplifts and remove uplift in ice cover periods
    delvmfull = (delvm + delvmsim) .* noicem;
    
    % water depth matrix (cm or g/cm2 if assuming a water density of 1)
    waterdm = (-delvmfull>elv).*(-delvmfull-elv).*1E2;
% end subfunction get_waterdepth ===================================================================


% subfunction add_tv_points ========================================================================
function [sample,LSDfix,pars] = add_tv_points(sample,LSDfix,iceproxytv,iceproxy,iceproxy_yr0);
% add time points for ice cover start/end and potential submergence
    % make iceproxy and iceproxytv row vectors
    iceproxy = iceproxy(:)'; iceproxytv = iceproxytv(:)';

    % parse time-varying icevalue and refer its time bands to samplingyr
    if isfield(sample,'icevaluev') && isfield(sample,'icevaluev_tv');
        sample.icevaluev = str2num(sample.icevaluev{1});
        sample.icevaluev_tv = str2num(sample.icevaluev_tv{1});
        if isfield(sample,'icevaluev_tv0'); % fix year 0
            sample.icevaluev_tv = sample.icevaluev_tv + sample.samplingyr - sample.icevaluev_tv0;
            sample = rmfield(sample,'icevaluev_tv0'); % already referred to samplingyr
        end;
    end;

    % adjust iceproxy and iceproxytv to LSDfix.tv
    iceproxytv = iceproxytv + sample.samplingyr - iceproxy_yr0;
    % first add one value outside each end to enable extrapolation of the end values
    iceproxytvpad = [iceproxytv(1)-1,iceproxytv,iceproxytv(end)+1];
    iceproxypad = [iceproxy(1),iceproxy,iceproxy(end)];
    % adjust to LSDfix.tv - based on the iceproxy ages only, not the added end values
    keepidx = (iceproxytv > LSDfix.tv(1)) & (iceproxytv < LSDfix.tv(end));
    proxytv = [LSDfix.tv(1),iceproxytv(keepidx),LSDfix.tv(end)];
    % add the icevaluev time bands to avoid quantising them to the iceproxy time step
    if isfield(sample,'icevaluev') && isfield(sample,'icevaluev_tv');
        edgeidx = (sample.icevaluev_tv > LSDfix.tv(1)) & (sample.icevaluev_tv < LSDfix.tv(end));
        proxytv = unique(sort([proxytv,sample.icevaluev_tv(edgeidx)'])); % drop duplicate bands
    end;
    proxy = interp1(iceproxytvpad,iceproxypad,proxytv,'linear','extrap');
    
    % fix icevalue vector if using time-varying icevalue
    if isfield(sample,'icevaluev') && isfield(sample,'icevaluev_tv');
        if isfield(sample,'icevalue')==0; sample.icevalue = 9; end; % fix for no icevalue input
        sample.icevalue(1,1:numel(proxytv)) = sample.icevalue;
        sample.tv = proxytv;
        % icevaluev_tv is already referred to samplingyr, so fix_minmaxv shifts it no further
        minmaxv = fix_minmaxv(sample,'icevaluev');
        for j = 1:size(minmaxv,1);
            sample.icevalue(1,minmaxv(j,1):minmaxv(j,2)) = sample.icevaluev(j);
        end;
    end;
    icevalue = sample.icevalue;

    % find break points based on proxy and icevalue
    gl1 = find(diff(proxy >= icevalue) == -1);
    degl1 = find(diff(proxy >= icevalue) == 1);
    gl2 = gl1 + 1;
    degl2 = degl1 + 1;

    % fix icevalue vectors for linear interpolation below
    if numel(sample.icevalue) == 1;
        ivg1 = ivg2 = ivd1 = ivd2 = icevalue;
    else;
        ivg1 = icevalue(gl1); ivg2 = icevalue(gl2);
        ivd1 = icevalue(degl1); ivd2 = icevalue(degl2);
    end;

    % find exact time points based on linear interpolation
    glactv = proxytv(gl1) + (ivg1-proxy(gl1)) .* (proxytv(gl2)-proxytv(gl1)) ./ ...
             (proxy(gl2)-proxy(gl1)-ivg2+ivg1);
    deglactv = proxytv(degl1) + (ivd1-proxy(degl1)) .* (proxytv(degl2)-proxytv(degl1)) ./ ...
               (proxy(degl2)-proxy(degl1)-ivd2+ivd1);

    % add potential deglac point to deglactv
    if isfield(sample,'deglac');
        deglactv = [deglactv,sample.deglac];
    end;
    
    % fix specified ice periods if using ice_tv
    if isfield(sample,'ice_tv');
        if isfield(sample,'ice_tv0'); % fix year 0
            icetv = str2num(sample.ice_tv{1}) + sample.samplingyr - sample.ice_tv0;
        else; icetv = str2num(sample.ice_tv{1}); end;
        for j = 1:size(icetv,1);
            % remove points within icetv periods
            glactv((glactv<=icetv(j,2))==(glactv>=icetv(j,1))) = [];
            deglactv((deglactv<=icetv(j,2))==(deglactv>=icetv(j,1))) = [];
            % add points from icetv to glactv and deglactv
            if (min(glactv(glactv>icetv(j,2)))>=min(deglactv(deglactv>icetv(j,2)))) || ...
                    (max(glactv)<icetv(j,2) && (max(glactv)>max(deglactv) || numel(deglactv)==0))...
                    || (numel(glactv)==0 && numel(deglactv)==0);
                glactv = [glactv,icetv(j,2)];
            end;
            if (min(glactv(glactv>icetv(j,1)))>=min(deglactv(deglactv>icetv(j,1)))) || ...
                    (max(glactv)<icetv(j,1) && (max(glactv)>max(deglactv) || numel(deglactv)==0))...
                    || (numel(glactv)==0 && numel(deglactv)==0);
                deglactv = [deglactv,icetv(j,1)];
            end;
        end;
    end;
    
    % fix specified ice free periods if using noice_tv
    if isfield(sample,'noice_tv');
        if isfield(sample,'noice_tv0'); % fix year 0
            noicetv = str2num(sample.noice_tv{1}) + sample.samplingyr - sample.noice_tv0;
        else; noicetv = str2num(sample.noice_tv{1}); end;
        for j = 1:size(noicetv,1);
            % remove points within noicetv periods
            glactv((glactv<=noicetv(j,2))==(glactv>=noicetv(j,1))) = [];
            deglactv((deglactv<=noicetv(j,2))==(deglactv>=noicetv(j,1))) = [];
            % add points from noicetv to glactv and deglactv
            if (min(glactv(glactv>noicetv(j,2)))>=min(deglactv(deglactv>noicetv(j,2)))) || ...
                    (max(deglactv)<noicetv(j,2) && (max(deglactv)>max(glactv) || ...
                    numel(glactv)==0)) || (max(glactv)>noicetv(j,2) && numel(deglactv)==0);
                glactv = [glactv,noicetv(j,2)];
            end;
            if (min(glactv(glactv>noicetv(j,1)))>=min(deglactv(deglactv>noicetv(j,1)))) || ...
                    (max(deglactv)<noicetv(j,1) && (max(deglactv)>max(glactv) || ...
                    numel(glactv)==0)) || (max(glactv)>noicetv(j,1) && numel(deglactv)==0);
                deglactv = [deglactv,noicetv(j,1)];
            end;
        end;
    end;

    % add +/- 1 yr for each point in glactv/deglactv
    sample.tv = sort([LSDfix.tv,glactv-1,glactv+1,deglactv-1,deglactv+1]);

    % remove potential duplicate time points
    sample.tv(diff(sample.tv)==0) = [];

    % interpolate Rc/SPhi/iceproxy to new tv
    LSDfix.Rc = interp1(LSDfix.tv,LSDfix.Rc,sample.tv);
    LSDfix.SPhi = interp1(LSDfix.tv,LSDfix.SPhi,sample.tv);
    pars.iceproxy = interp1(proxytv,proxy,sample.tv(:));
    if numel(sample.icevalue) > 1;
        sample.icevalue = interp1(proxytv,sample.icevalue,sample.tv(:));
    end;
    sample.tv = sample.tv(:);

    % fix submergence vectors if using isostatic adjustment
    if isfield(sample,'isostsubm');
        sample.submelv = isost_elv(sample.isostsubm{1},sample);
        sample.delv = sample.submelv' - sample.elv;
        if isfield(sample,'isostsubmmod') && sample.isostsubmmod>0;
            sample.delv = sample.delv .* sample.isostsubmmod;
        end;
        % fix ice cover record
        icem = (pars.iceproxy >= sample.icevalue);
        icem = double(icem); % fix for matlab
        % fix deglaciation based on sample.deglac
        if isfield(sample,'deglac');
            icem = fix_deglac(sample.deglac,icem,sample.tv);
        end;
        % fix specified ice-cover/noice-cover periods
        if isfield(sample,'ice_tv') || isfield(sample,'noice_tv');
            icem = fixicem(sample,icem);
        end;
        % calculate water depth
        [pars.waterdm,delvmfull] = get_waterdepth(sample.elv,sample.delv,icem,sample.tv);
        % pick out boundary points for water/air change during ice-free periods
        bp1 = find(diff(-delvmfull > sample.elv) ~= 0); % boundary points 1
        bp1(icem(bp1+1) == 1) = []; % remove points associated with deglaciation
        bp2 = bp1 + 1; % boundary points 2
        % check if any elevation boundary point is exactly 0 and remove
        rmv = find(-delvmfull(bp1) == sample.elv); bp1(rmv) = []; bp2(rmv) = [];
        rmv = find(-delvmfull(bp2) == sample.elv); bp1(rmv) = []; bp2(rmv) = [];
        % calculate ratio of step between sample.delv(bp1) and sample.delv(bp2)
        stepratio = (-delvmfull(bp1)-sample.elv)./(delvmfull(bp2)-delvmfull(bp1));
        % fix updated tv with added points for water/air boundaries
        newtv = sort([sample.tv;sample.tv(bp1)+(sample.tv(bp2)-sample.tv(bp1)).*stepratio]);
        % interpolate Rc, SPhi, iceproxy, waterdm, delv to newtv
        LSDfix.Rc = interp1(sample.tv,LSDfix.Rc,newtv);
        LSDfix.SPhi = interp1(sample.tv,LSDfix.SPhi,newtv);
        pars.iceproxy = interp1(sample.tv,pars.iceproxy,newtv);
        pars.waterdm = interp1(sample.tv,pars.waterdm,newtv);
        sample.delv = interp1(sample.tv,sample.delv,newtv);
        if numel(sample.icevalue) > 1;
            sample.icevalue = interp1(sample.tv,sample.icevalue,newtv);
        end;
        sample.tv = newtv;
    end;
% end subfunction add_tv_points ====================================================================


% subfunction evm ==================================================================================
function [mid,uncp,uncn] = evm(midv,uncpv,uncnv);
    % average and uncertainty estimation using the "expected value method" (Birch and Singh 2014)
    num = numel(midv); % number of input samples
    
    % fix data (make matrices)
    midm = repmat(midv,num,1);
    uncpm = repmat(uncpv,num,1);
    uncnm = repmat(uncnv,num,1);
    xm = repmat(midv',1,num);
    
    % calculate probability matrix
    Mui = sum(sqrt(2./(pi.*(uncpm+uncnm).^2)).*exp(-(xm-midm).^2./(2.*uncnm.^2)).*(xm<=midm) + ...
        sqrt(2./(pi.*(uncpm+uncnm).^2)).*exp(-(xm-midm).^2./(2.*uncpm.^2)).*(xm>midm)) ./ num;
    
    % calculate summed probability for all samples
    Muj = sum(Mui);
    
    % calculate sample weights
    wi = Mui./Muj;
    
    % calculate weighted mean mid
    mid = sum(wi.*midv);
    
    % uncertainty estimation
    uncintp = sqrt(sum(wi.^2.*uncpv.^2));
    uncintn = sqrt(sum(wi.^2.*uncnv.^2));
    uncext = sqrt(sum(wi.*(midv-mid).^2));
    uncp = max(uncintp,uncext);
    uncn = max(uncintn,uncext);
% end subfunction evm ==============================================================================


% subfunction customcolormap =======================================================================
function J = customcolormap(positions,colors,m);
    % Compute positions along the samples
    color_samples = round((m-1)*positions)+1;
    % Make the gradients among colors
    J = zeros(m,3);
    J(color_samples,:) = colors;
    diff_samples = diff(color_samples)-1;
    for d = 1:1:length(diff_samples)
        if diff_samples(d)~=0
            color1 = colors(d,:);
            color2 = colors(d+1,:);
            G = zeros(diff_samples(d),3);
            for idx_rgb = 1:3
                g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                g([1, length(g)]) = [];
                G(:,idx_rgb) = g';
            end
            J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
        end
    end
    J = flipud(J);
% end subfunction customcolormap ===================================================================
