## Copyright (C) 2001 Paulo Neis <p_neis@yahoo.com.br>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING. If not, see
## <https://www.gnu.org/licenses/>.

## LOCAL SHADOW of signal-package ncauer.m.
##
## Two changes from upstream
##   https://github.com/gnu-octave/octave-signal/blob/main/inst/ncauer.m:
##
##   1. __ellip_ws invokes fminbnd with an explicit TolX of 1e-15 instead of
##      the ~1e-4 Octave default.  The loose default propagated into analog pole
##      locations at the ~1e-6 level, corrupting the reference.
##
##   2. Step comments added to ncauer() for readability.
##
## To verify this file against upstream, diff it against the URL above.  The
## only algorithmic difference should be the fminbnd options struct inside
## __ellip_ws.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{z}, @var{p}, @var{g}] =} cauer(@var{Rp}, @var{Rs}, @var{n})
## Analog prototype for Cauer filter.
##
## @table @asis
## @item Rp
## Passband ripple
## @item Rs
## Stopband ripple
## @item n
## Desired order
## @item z
## complex vector of zeros for the model.
## @item p
## complex vector of poles for the model.
## @item g
## gain value.
## @end table
##
## References:
##
## - Serra, Celso Penteado, Teoria e Projeto de Filtros, Campinas: CARTGRAF,
##   1983.
##
## - Lamar, Marcus Vinicius, Notas de aula da disciplina TE 456 - Circuitos
##   Analogicos II, UFPR, 2001/2002.
## @end deftypefn

function [zer, pol, T0] = ncauer (Rp, Rs, n)

  ## Step 1 & 2 — design modulus k and nome q.
  ##   ncauer finds ws first via an iterative degree-equation solver, then
  ##   derives k = 1/ws and the nome q from k.  (constfilt takes a closed-form
  ##   path instead; only this step differs between the two implementations.)
  wp=1;
  ws=__ellip_ws(n, Rp, Rs);
  k=wp/ws;
  k1=sqrt(1-k^2);
  q0=(1/2)*((1-sqrt(k1))/(1+sqrt(k1)));
  q= q0 + 2*q0^5 + 15*q0^9 + 150*q0^13; #(....)
  D=(10^(0.1*Rs)-1)/(10^(0.1*Rp)-1);

  ## Filter order maybe this, but not used now:
  ## n=ceil(log10(16*D)/log10(1/q))

  ## Step 3 — pole-shift parameter sig0.
  ##   l encodes the passband-ripple spec as a hyperbolic angle, divided by n
  ##   to spread it across the n pole positions.  sig01/sig02 are theta-function
  ##   series that evaluate the Jacobi elliptic function at that angle.
  l=(1/(2*n))*log((10^(0.05*Rp)+1)/(10^(0.05*Rp)-1));
  sig01=0; sig02=0;
  for m=0 : 30
    sig01=sig01+(-1)^m * q^(m*(m+1)) * sinh((2*m+1)*l);
  endfor
  for m=1 : 30
    sig02=sig02+(-1)^m * q^(m^2) * cosh(2*m*l);
  endfor
  sig0=abs((2*q^(1/4)*sig01)/(1+2*sig02));

  ## Step 4 — zero positions wi (one per conjugate pair).
  ##   wi = sn(i*K(k)/n, k) via theta-function series; mu shifts the index for
  ##   even orders so zeros land at half-integer multiples of K/n.
  w=sqrt((1+k*sig0^2)*(1+sig0^2/k));
  if rem(n,2)
    r=(n-1)/2;
  else
    r=n/2;
  endif
  wi=zeros(1,r);
  for ii=1 : r
    if rem(n,2)
      mu=ii;
    else
      mu=ii-1/2;
    endif
    soma1=0;
    for m=0 : 30
      soma1 = soma1 + 2*q^(1/4) * ((-1)^m * q^(m*(m+1)) * sin(((2*m+1)*pi*mu)/n));
    endfor
    soma2=0;
    for m=1 : 30
      soma2 = soma2 + 2*((-1)^m * q^(m^2) * cos((2*m*pi*mu)/n));
    endfor
    wi(ii)=(soma1/(1+soma2));
  endfor

  ## Step 5 — pole and zero placement.
  ##   Vi = cn(ui,k)*dn(ui,k) via Pythagorean identities on wi = sn(ui,k).
  ##   B0i/B1i are the quadratic and linear denominator coefficients per pair.
  ##   sqrt(ws) scaling puts the stopband edge at ws rad/s.
  Vi=sqrt((1-(k.*(wi.^2))).*(1-(wi.^2)/k));
  A0i=1./(wi.^2);
  sqrA0i=1./(wi);
  B0i=((sig0.*Vi).^2 + (w.*wi).^2)./((1+sig0^2.*wi.^2).^2);
  B1i=(2 * sig0.*Vi)./(1 + sig0^2 * wi.^2);

  ## Step 6 — gain normalization.
  ##   odd n:  H(0) = 1  (DC gain exactly 1)
  ##   even n: H(0) = Gp = 10^(-0.05*Rp)  (DC gain at passband floor)
  if rem(n,2)
    T0=sig0*prod(B0i./A0i)*sqrt(ws);
  else
    T0=10^(-0.05*Rp)*prod(B0i./A0i);
  endif

  ## Step 7 — zeros, poles, stopband scaling.
  zer=[i*sqrA0i, -i*sqrA0i];
  pol=[(-2*sig0*Vi+2*i*wi.*w)./(2*(1+sig0^2*wi.^2)), (-2*sig0*Vi-2*i*wi.*w)./(2*(1+sig0^2*wi.^2))];

  ## If n odd, there is a real pole at -sig0:
  if rem(n,2)
    pol=[pol, -sig0];
  endif

  pol=(sqrt(ws)).*pol;
  zer=(sqrt(ws)).*zer;

endfunction

## usage: ws = __ellip_ws(n, rp, rs)
## Calculate the stop band edge for the Cauer filter.

function ws=__ellip_ws(n, rp, rs)

  kl0 = ((10^(0.1*rp)-1)/(10^(0.1*rs)-1));
  k0  = (1-kl0);
  int = ellipke([kl0 ; k0]);
  ql0 = int(1);
  q0  = int(2);
  x   = n*ql0/q0;
  ## PATCH (one of two changes vs upstream
  ##   https://github.com/gnu-octave/octave-signal/blob/main/inst/ncauer.m):
  ## upstream passes no options; we force TolX = 1e-15 so fminbnd converges
  ## to kl at machine precision instead of the default ~1e-4 that
  ## contaminates downstream pole locations.
  kl  = fminbnd(@(y) __ellip_ws_min(y,x), eps, 1-eps, ...
                optimset('TolX', 1e-15, 'MaxIter', 10000));
  ws  = sqrt(1/kl);

endfunction

## usage: err = __ellip_ws_min(kl, x)

function err=__ellip_ws_min(kl, x)

  int=ellipke([kl; 1-kl]);
  ql=int(1);
  q=int(2);
  err=abs((ql/q)-x);

endfunction
