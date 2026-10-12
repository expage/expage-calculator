function out = LSD_fix(lat,lon,age,samplingyr,consts)

% Fix w, Rc, and SPhi to be used by P_mu_LSD and LSDspal.
% syntax : LSD_fix(lat,lon,age,w,consts);

% lat = sample latitude in deg N (negative values for S hemisphere)
% lon = sample longitude in deg E (negative values for W longitudes, or 0-360 degrees E)
% age = age for production scaling
% w = gravimetric fractional water content - 0.066 is default typically about 14% volumetric per
%     Fred Phillips. -1 gives default value
% samplingyr = year of sampling
% consts = consts_LSD (LSD constants structure)
%
% Input values as scalars
%
% Based on code written by Greg Balco -- Berkeley Geochronology Center
% balcs@bgc.org
% 
% Modified by Brent Goehring and  Nat Lifton -- Purdue University
% nlifton@purdue.edu, bgoehrin@purdue.edu
%
% Modified by Jakob Heyman (jakob.heyman@gu.se) 2015-2026
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 3, as published by the Free Software Foundation (www.fsf.org).

% Make the time vector
% Age Relative to t0=2010
tv1 = [0:10:50 60:100:3060 3260:200:75060 76060:1000:800060 802060:2000:2000060 logspace(log10(2002060),7,200)];

Rc = zeros(1,length(tv1));

% Need solar modulation parameter
SPhi = zeros(size(tv1)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
SPhi(1:78) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)

out.w = 0.066; % water content for Sato & Niita (2006), gravimetric, about 14% volumetric per Fred Phillips

% interpolate an M for tv > 0 for dipole field. 10/1/15 Pavon et al. SHA.DIF.14k
temp_M1 = interp1(consts.tRc,consts.MM0_Pavon,tv1(1:91));
temp_M2 = interp1(consts.tM,consts.GLOPADMM0,tv1(92:end));
temp_M = cat(2,temp_M1,temp_M2);

% TIME-VARYING GAD FIELD
%   Lifton et al. 2014 -  Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average.
GDRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];

% catch for negative longitudes before Rc interpolation
if lon < 0; lon = lon + 360; end;

% Make up the Rc vectors.
% Geocentric dipole 0-14 ka (Pavon-Carrasco et al., 2014),
% similar to University of Washington online calculator version 3
LDLat = abs(90-angdist(lat,lon,consts.lat_pp_Pavon,consts.lon_pp_Pavon));
Rc(1:91) = temp_M(1:91).*polyval(GDRc,cos(d2r(LDLat)));
Rc(92:end) = temp_M(92:end).*polyval(GDRc,cos(d2r(lat)));

zIndex = Rc(:) < 0;
Rc(zIndex) = 0;

% Time-average constant RcEst for muon production
out.RcEst = trapz(tv1,Rc)/tv1(end);

% adjust tv, tv1, Rc, and SPhi to sampling year
if samplingyr < 2010;
    clipidx = find(tv1 > 2010-samplingyr,1,'first');
    tv = [0 (tv1(clipidx:end)-2010+samplingyr)];
    tv1 = [(tv1-2010+samplingyr) tv1(end)];
    Rc = [Rc Rc(end)];
    SPhi = [SPhi SPhi(end)];
elseif samplingyr > 2010; % assume 2010 Rc and SPhi value for all years >2010
    tv = [0 (tv1+samplingyr-2010)];
    tv1 = tv;
    Rc = [Rc(1) Rc];
    SPhi = [SPhi(1) SPhi];
elseif samplingyr == 2010;
    tv = tv1;
end;

% chop off tv at age
clipidx = find(tv < age,1,'last');
out.tv = [tv(1:clipidx) age];

% shorten Rc and SPhi commensurately
out.Rc = interp1(tv1,Rc,out.tv);
out.SPhi = interp1(tv1,SPhi,out.tv);


% subfunction angdist ==============================================================================
function out = angdist(lat1,long1,lat2,long2);
    lata = d2r(lat1);longa = d2r(long1);
    latb = d2r(lat2);longb = d2r(long2);

    theta = acos( (cos(lata).*cos(latb).*((cos(longa).*cos(longb)) + (sin(longa).*sin(longb)))) + (sin(lata).*sin(latb)));
    out_1 = r2d(theta);

    % make vector horizontal

    if size(out_1,1) > size(out_1,2);
        out = out_1';
    else;
        out = out_1;
    end;
% end subfunction angdist ==========================================================================


% subfunction d2r ==================================================================================
function out = d2r(in);
    out = (in./360).*2.*pi;
% end subfunction d2r ==============================================================================


% subfunction r2d ==================================================================================
function out = r2d(in);
    out = (in*360)/(2*pi);
% end subfunction r2d ==============================================================================