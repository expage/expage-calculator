function out = make_consts_expage();

% This function creates and saves a structure with relevant constants and external data for the
% expage calculator.
% Based on make_al_be_consts_v22.m (Balco et al. 2008) and make_consts_LSD.m Lifton et al. (2014),
% created and modified by Greg Balco, Brent Goehring, Nat Lifton.
%
% Paleomagnetic records updated to Lifton (2016) data: PMag_Sep12.mat changed to PavonPmag.mat.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2018-2026 (jakob.heyman@gu.se)

consts.version = '202610';
consts.prepdate = fix(clock);

% Be-10 decay constant -- new Euro value
consts.l10 = -log(0.5)./1.387e6; % Chmeleff/Korschinek value
dldt = -log(0.5).*(1.387e6^-2);
consts.l10unc = sqrt((dldt.*0.012e6)^2); % Chmeleff/Korschinek value

% Al-26 decay constant -- value compatible with Nishiizumi standards
% lambda = 9.83e-7 --> t(1/2) = 0.705e6 yr
% See Nishiizumi (2004) for details.
consts.l26 = 9.83e-7; 
consts.l26unc = 2.5e-8;

% C-14 decay constant -- t(1/2) = 5700 ± 30 yr (Kutschera 2013)
consts.l14 = -log(0.5)./5700;
dldt = -log(0.5).*(5700^-2);
consts.l14unc = sqrt((dldt.*30)^2);

% Effective attenuation length for spallation in rock
% Commonly accepted value: 160 g/cm2 (Gosse and Phillips 2001)
% Only used for simple age calculation - else rawattenuationlength is used.
consts.Lsp = 160;

% Be-10 standardization info.
% Standards comparison/conversion lookup table. A zero placeholder is allowed.
consts.std10 = {'07KNSTD','KNSTD','NIST_Certified','NIST_30000','NIST_30200','NIST_30300',...
    'NIST_30500','NIST_30600','NIST_27900','LLNL31000','LLNL10000','LLNL3000','LLNL1000',...
    'LLNL300','S555','S2007','BEST433','S555N','S2007N','S2010N','BEST433N','STD11','SMDBe12',...
    '0','-'};
consts.std10_cf = [1 0.9042 1.0425 0.9313 0.9251 0.9221 0.9157 0.9130 1 0.8761 0.9042 0.8644 ...
    0.9313 0.8562 0.9124 0.9124 0.9124 1 1 1 1 1 1 1 1]';

% Same for Al-26. A zero placeholder is allowed.
consts.std26 = {'KNSTD','ZAL94','ZAL94N','SMAL11','ASTER','Z92-0222','SMDAl11','0','-'};
consts.std26_cf = [1 0.9134 1 1.021 1.021 1 1.018 1 1]';

% Reference production rates at SLHL for spallation. Letter codes sp and nu refer to the LSD
% spallation scaling and the LSD nuclide-specific scaling of Lifton et al. (2014), respectively.
% 10Be production rates are referenced to 07KNSTD.
% 26Al production rates are referenced to KNSTD.

% Be-10 production rates - global expage-202610 ref prod rate
consts.Pref10 = 4.17;
consts.Pref10unc = 0.25;
consts.Pref10iso = 4.18;
consts.Pref10isounc = 0.25;
consts.P10_ref_sp = 4; % not used and not calibrated!
consts.delP10_ref_sp = 0.3; % not used and not calibrated!

% Al-26 production rates - global expage-202610 ref prod rate
consts.Pref26 = 28.48;
consts.Pref26unc = 2.36;
consts.Pref26iso = 28.36;
consts.Pref26isounc = 2.49;
consts.P26_ref_sp = 28; % not used and not calibrated!
consts.delP26_ref_sp = 3; % not used and not calibrated!

% C-14 production rates - global expage-202610 ref prod rate
consts.Pref14 = 12.75;
consts.Pref14unc = 1.46;
consts.Pref14iso = 12.60;
consts.Pref14isounc = 1.40;

% % GLOBAL LOCKNE DATA =================================================
% % Be-10 production rates - global expage-202610 ref prod rate
% consts.Pref10 = 4.18;
% consts.Pref10unc = 0.25;
% consts.Pref10iso = 4.20;
% consts.Pref10isounc = 0.26;
% consts.P10_ref_sp = 4; % not used and not calibrated!
% consts.delP10_ref_sp = 0.3; % not used and not calibrated!
%
% % Al-26 production rates - global expage-202610 ref prod rate
% consts.Pref26 = 28.80;
% consts.Pref26unc = 2.46;
% consts.Pref26iso = 28.93;
% consts.Pref26isounc = 2.80;
% consts.P26_ref_sp = 28; % not used and not calibrated!
% consts.delP26_ref_sp = 3; % not used and not calibrated!
%
% % C-14 production rates - global expage-202610 ref prod rate
% consts.Pref14 = 12.89;
% consts.Pref14unc = 1.45;
% consts.Pref14iso = 12.90;
% consts.Pref14isounc = 1.46;
% % ====================================================================

% Atomic number densities (atoms target/g mineral)
% (Moles target * Avogadro's Number / Formula weight)
% Quartz (SiO2)
consts.NatomsQtzO = 2.00456e22;
consts.NatomsQtzSi = 1.00228e22;

% Muon interaction cross-sections. All follow Heisinger (2002a,b).
consts.Natoms10 = 2.006e22;
consts.Natoms26 = 1.003e22;
consts.Natoms14 = 2.006e22;
consts.Natoms3 = 2.006e22; % not used in present version

% sigma0 and fstar calibrated for the expage calculator using mucalib.m (mucalib14.m for 14C)
% (ref prodrate expage-202610)
consts.sigma0_10nu = 2.536e-31;
consts.delsigma0_10nu = 0.562e-31;
consts.fstar10nu = 1.755e-3;
consts.delfstar10nu = 0.378e-3;
consts.sigma0_26nu = 3.748e-30;
consts.delsigma0_26nu = 1.559e-30;
consts.fstar26nu = 1.283e-2;
consts.delfstar26nu = 0.536e-2;
consts.sigma0_14nu = 0.769e-31;
consts.delsigma0_14nu = 123.793e-31;
consts.fstar14nu = 1.321e-1;
consts.delfstar14nu = 1.021e-1;

% % GLOBAL LOCKNE DATA =================================================
% % sigma0 and fstar calibrated for the expage calculator using mucalib.m (mucalib14.m for 14C)
% consts.sigma0_10nu = 2.538e-31;
% consts.delsigma0_10nu = 0.564e-31;
% consts.fstar10nu = 1.754e-3;
% consts.delfstar10nu = 0.381e-3;
% consts.sigma0_26nu = 3.733e-30;
% consts.delsigma0_26nu = 1.557e-30;
% consts.fstar26nu = 1.292e-2;
% consts.delfstar26nu = 0.536e-2;
% consts.sigma0_14nu = 0.739e-31;
% consts.delsigma0_14nu = 123.833e-31;
% consts.fstar14nu = 1.314e-1;
% consts.delfstar14nu = 1.021e-1;
% % ====================================================================

% sigma0 and fstar from CRONUScalc calculator (Marrero et al. 2016; Phillips et al. 2016)
% not used in present version
consts.sigma0_10sp = 0.252E-30;
consts.sigma0_26sp = 4.10E-30;
consts.delsigma0_10sp = 0.015e-30;
consts.delsigma0_26sp = 0.73e-30;
consts.sigma0_14sp = 8.79321E-30;
consts.fstar10sp = 1.89e-3;
consts.fstar26sp = 11.7e-3;
consts.delfstar10sp = 0.18e-3;
consts.delfstar26sp = 2.6e-3;
consts.fstar14sp = 0.137;

% k_negpartial from Balco 2017
consts.k_negpartial_10 = (0.704 * 0.1828); % per Beacon Core Fit (Balco, 2017) 8/8/17 NL
consts.delk_negpartial_10 = (0.704 * 0.1828 * 0.0003); % per Beacon Core Fit (Balco, 2017) 8/8/17 NL
% consts.k_negpartial_10 = (0.704 * 0.1828)./1.106; % per CRONUScalc (Marrero et al., 2016)
% consts.delk_negpartial_10 = (0.704 * 0.1828 * 0.0003)./1.106; % per CRONUScalc (Marrero et al., 2016)
consts.k_negpartial_14 = 0.704 * 0.1828;
consts.delk_negpartial_14 = 0.704 * 0.1828 * 0.0011;
consts.k_negpartial_26 = 0.296 * 0.6559;
consts.delk_negpartial_26 = 0.296 * 0.6559 * 0.002;

% Spallogenic Nuclide Production Cross-Sections (n & p) from Bob Reedy 9/2010
% High-resolution cross-sections (800 energy bins instead of original 200)
load XSectionsHR;

consts.E = E;

consts.O16nxBe10 = Onx10Be;
consts.O16pxBe10 = Opx10Be;
consts.SinxBe10 = Sinx10Be;
consts.SipxBe10 = Sipx10Be;
consts.O16nn2pC14 = Onn2p14C;
consts.O16pxC14 = Opx14C;
consts.SinxC14 = Sinx14C;
consts.SipxC14 = Sipx14C;
consts.Aln2nAl26 = Aln2n26Al;
consts.AlppnAl26 = Alppn26Al;
consts.SinxAl26 = Sinx26Al;
consts.SipxAl26 = Sipx26Al;

% Paleomagnetic records for use in time-dependent production rate schemes
% Derived from Nat Lifton's compilation of paleomagnetic data from various sources.
% Load the magnetic field data - Lifton (2016) version:
% 2010-1950: Rc from DGRFs
% 0-14 ka:   SHA.DIF.14k
% 14-75 ka:  GLOPIS-75
% >75 ka:    PADM2M
load PavonPmag %Pavon to 14kyr, GLOPIS-75 to 75 ka, PADM2M to 2Ma, 1950-2010 Rc. Same S and SPhi temporal spacing as PMag_Aug13. 7/11/14

% Relative dipole moment and time vector - Pavon 7/11/14
consts.GLOPADMM0 = GLOPADMM0;
consts.tM = tM;

consts.PavonRc = PavonRc; % gridded Rc data block - 7/11/14

consts.lat_Rc = lat_Rc; % lat and lon indices for Rc data block
consts.lon_Rc = lon_Rc;
consts.tRc = tRc; % time vector for Rc data block

% These start at 14000 yr -- time slices are 100-yr from 0 to 3000, 200-yr from 3000 to 75000,
% 1000-yr from 75-ka to 800-ka, 2000-yr from 800-ka to 2-Ma, and 10-Ma at the end. The interval
% spacing is adjusted from that of Lifton et al. (2014) to improve computational speed while
% remaining faithful to the temporal variability in each model.

% Cutoff rigidity blocks for past 14000 yr. 
% PavonRc has lat x lon x time blocks of Rc values for the past 14000 years.
% Cutoff rigidity obtained by Nat Lifton using trajectory tracing from the SHA.DIF.14k magnetic
% field model of Pavón-Carrasco et al. (2014).  
consts.PavonRc = PavonRc; % data block
consts.lat_Rc = lat_Rc; % lat and lon indices for Rc data block
consts.lon_Rc = lon_Rc;
consts.tRc = tRc; % time vector for Rc data block

% Effective pole positions and field strengths inferred from Pavon field
% reconstructions for last 14000 yr. These are used in the
% paleomagnetically-corrected implementation of the Lal SF. They are for
% the same times as the RC slices in the data block above. Generated by Nat Lifton
consts.MM0_Pavon = MM0_Pavon;
consts.lat_pp_Pavon = lat_pp_Pavon;
consts.lon_pp_Pavon = lon_pp_Pavon;

% Solar variability from Usoskin et al. 2011
% 0-11260 yr

% Per Tatsuhiko Sato, personal communication, 2013, convert annually averaged Usoskin et al. (2011)
% solar modulation potential to Sato Force Field Potential due to different assumed Local
% Interstellar Spectrum and other factors
SPhi = 1.1381076.*SPhi - 1.2738468e-4.*SPhi.^2;

consts.SPhi = SPhi;
consts.SPhiInf = trapz(tRc(1:78),SPhi)/tRc(78); % time-averaged SPhi

%Calculated reference PRs using hi-resolution cross-sections (800 energy
%bins instead of original 200, per Argento et al., 2015) 7/20/15
load ReferenceHR
%Reference values for scaling via Sato et al. (2008) spectra
% consts.E = E;

consts.P10nRef_q = P10nRef_q; %10Be neutron reference production in SiO2
consts.P10pRef_q = P10pRef_q; %10Be proton reference production in SiO2
consts.P14nRef_q = P14nRef_q; %14C neutron reference production in SiO2
consts.P14pRef_q = P14pRef_q; %14C proton reference production in SiO2
consts.P26nRef_q = P26nRef_q; %26Al neutron reference production in SiO2
consts.P26pRef_q = P26pRef_q; %26Al proton reference production in SiO2
consts.nfluxRef = nfluxRef; %integral reference neutron flux
consts.pfluxRef = pfluxRef; %integral reference proton flux
consts.ethfluxRef = ethfluxRef; %integral reference epithermal flux
consts.thfluxRef = thfluxRef; %integral reference thermal neutron flux
consts.mfluxRef = mfluxRef; % reference muon flux components

% Finish up
save consts_expage consts
disp(['consts_expage v. ' consts.version ' saved']);
