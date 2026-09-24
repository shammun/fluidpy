"""Build the Chapter 7 teaching notebook: ``notebooks/ch07_gravity_waves.ipynb``.

Source of truth: ``analysis/ch07_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract
— every call goes to ``fluidpy.ch07_gravity_waves``, imported as ``ch07``, which re-exports the new core module
``core.waves`` (imported as ``W``)), Part E (prerequisite ledger → 14 primers P165–P178 and one-line reminders of earlier
primers), Part F (the 37 derivations D01–D37, one move per step) and ``analysis/ch07_curation.md`` (IDs, depths,
section coverage). Content decisions from ``reports/ch07_verification.md`` override the design where they differ:
the exact-condition residuals of the linear solution scale as (ka)·aω — relative O(ka) — (O3, following D03, not
analysis V25); Stokes' (7.82)–(7.83) is presented with the consistent third-order expansion (γ = 1) and the reason the
literal truncated Exercise-7.2 set-up gives γ = 3/8 (O4); the interfacial potential energy is ¼Δρga²; the internal-wave
group velocity for k < 0 is ∇_K ω (the printed (7.145) holds for k > 0 only); the explicit dispersion approximations are
quoted with the primary bounds (1.7 % and 0.8 %, O2). Physics lives in ``fluidpy``; cells only call it.

**Derivations are read from Part F at build time** (``part_f()`` below, the ch03–ch06 parser, extended for the
explainer-only step fields ``*live (E1):*``, ``*set (E4):*``, ``*watch (E2):*``, which are dropped): goal, start, plan,
tools, assumptions, every step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word.
The four ★★★ sympy checks (D20, D29, D30, D33) are written here, every line commented, and re-run the construction
(D20's Gaussian integral is done once in closed form instead of the design sketch's complex ``integrate``, which took
33 s).

Book numbers never printed (rule 9): our own worked numbers only; (7.59) is shown with our IAPWS values.

Run:  .venv/Scripts/python.exe notebooks/build_ch07.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch07.py --dump     (prints the parsed Part F derivations only)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch07")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# Ch. 7 from analysis/ch07.md §2 (transcribed from the rendered pages p281–p331), corrected forms for the printed slips
# (7.66) ½Δω t, (7.105) e^{i(kx−ωt)}, (7.138)/(7.145) for k > 0; earlier chapters from their notebooks.
# ---------------------------------------------------------------------------------------------------------------------
EQ = {
    "7.1": r"\eta(x,t)=a\cos\Big[\frac{2\pi}{\lambda}(x-ct)\Big]",
    "7.2": r"\eta(x,t)=a\cos[kx-\omega t]",
    "7.3": r"\frac{2\pi}{\lambda}(x_{crest}-ct)=2n\pi=kx_{crest}-\omega t",
    "7.4": r"c=\omega/k=\lambda\nu",
    "7.5": r"\eta=a\cos(kx+ly+mz-\omega t)=a\cos(\mathbf K\cdot\mathbf x-\omega t)",
    "7.6": r"K^2=k^2+l^2+m^2",
    "7.7": r"\lambda=2\pi/K",
    "7.8": r"\mathbf c=(\omega/K)\,\mathbf e_K,\ \mathbf e_K=\mathbf K/K",
    "7.9": r"\omega_0=\omega+\mathbf U\cdot\mathbf K",
    "7.10": r"u=\partial\phi/\partial x,\ w=\partial\phi/\partial z",
    "7.11": r"\partial^2\phi/\partial x^2+\partial^2\phi/\partial z^2=0",
    "7.12": r"w=\partial\phi/\partial z=0\ \text{on}\ z=-H",
    "7.13": r"(\mathbf n\cdot\mathbf u)_{z=\eta}=\mathbf n\cdot\mathbf U_s",
    "7.14": r"\mathbf n=\big(-(\partial\eta/\partial x)\mathbf e_x+\mathbf e_z\big)\big/\sqrt{(\partial\eta/\partial x)^2+1}",
    "7.15": r"\mathbf U_s=(\partial\eta/\partial t)\,\mathbf e_z",
    "7.16": r"\Big(\frac{\partial\phi}{\partial z}\Big)_{z=\eta}=\frac{\partial\eta}{\partial t}+\frac{\partial\eta}{\partial x}\Big(\frac{\partial\phi}{\partial x}\Big)_{z=\eta}",
    "7.17": r"\Big(\frac{\partial\phi}{\partial z}\Big)_{z=\eta}\cong\frac{\partial\eta}{\partial t}",
    "7.18": r"\Big(\frac{\partial\phi}{\partial z}\Big)_{z=0}\cong\frac{\partial\eta}{\partial t}",
    "7.19": r"(p)_{z=\eta}=0",
    "7.20": r"\frac{\partial\phi}{\partial t}+\frac p\rho+gz\cong0",
    "7.21": r"\Big(\frac{\partial\phi}{\partial t}\Big)_{z=\eta}\cong\Big(\frac{\partial\phi}{\partial t}\Big)_{z=0}\cong-g\eta",
    "7.22": r"\phi(x,z,t)=f(z)\sin(kx-\omega(k)t)",
    "7.23": r"\frac{d^2f}{dz^2}-k^2f=0,\ f=Ae^{kz}+Be^{-kz}",
    "7.24": r"k(Ae^{-kH}-Be^{+kH})\sin(kx-\omega t)=0\ \text{or}\ B=Ae^{-2kH}",
    "7.25": r"k(A-B)=\omega a",
    "7.26": r"\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)",
    "7.27": r"u=a\omega\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t),\ w=a\omega\frac{\sinh k(z+H)}{\sinh kH}\sin(kx-\omega t)",
    "7.28": r"\omega=\sqrt{gk\tanh kH}",
    "7.29": r"c=\sqrt{\frac gk\tanh kH}",
    "7.30": r"p'\equiv p+\rho gz",
    "7.31": r"p'=\rho ga\frac{\cosh k(z+H)}{\cosh kH}\cos(kx-\omega t)",
    "7.32": r"\frac{dx_p}{dt}=u(x_p,z_p,t),\ \frac{dz_p}{dt}=w(x_p,z_p,t)",
    "7.33": r"\frac{dx_p}{dt}=a\omega\frac{\cosh k(z_p+H)}{\sinh kH}\cos(kx_p-\omega t),\ \frac{dz_p}{dt}=a\omega\frac{\sinh k(z_p+H)}{\sinh kH}\sin(kx_p-\omega t)",
    "7.34a": r"\frac{d\xi}{dt}\cong a\omega\frac{\cosh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)",
    "7.34b": r"\frac{d\zeta}{dt}\cong a\omega\frac{\sinh k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)",
    "7.35a": r"\xi\cong-a\frac{\cosh k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)",
    "7.35b": r"\zeta\cong a\frac{\sinh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)",
    "7.36": r"\frac{\xi^2}{A^2}+\frac{\zeta^2}{B^2}=1,\ A=a\frac{\cosh k(z_0+H)}{\sinh kH},\ B=a\frac{\sinh k(z_0+H)}{\sinh kH}",
    "7.37": r"\psi=\frac{a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}\cos(kx-\omega t)",
    "7.38": r"E_k=\frac{\rho}{2\lambda}\int_0^\lambda\!\!\int_{-H}^0(u^2+w^2)\,dz\,dx",
    "7.39": r"E_k=\tfrac12\rho g\overline{\eta^2}",
    "7.40": r"E_p=\frac{\rho g}{2\lambda}\int_0^\lambda\eta^2\,dx",
    "7.41": r"E_p=\tfrac12\rho g\overline{\eta^2}",
    "7.42": r"E=E_p+E_k=\rho g\overline{\eta^2}=\tfrac12\rho ga^2",
    "7.43": r"F=\Big\langle\int_{-H}^0p'u\,dz\Big\rangle",
    "7.44": r"F=\Big[\tfrac12\rho ga^2\Big]\Big[\frac c2\Big(1+\frac{2kH}{\sinh2kH}\Big)\Big]",
    "7.45": r"c=\sqrt{g/k}=\sqrt{g\lambda/2\pi}",
    "7.46": r"\xi\cong-ae^{kz_0}\sin(kx_0-\omega t),\ \zeta\cong ae^{kz_0}\cos(kx_0-\omega t)",
    "7.47": r"u=a\omega e^{kz}\cos(kx-\omega t),\ w=a\omega e^{kz}\sin(kx-\omega t)",
    "7.48": r"p'=\rho gae^{kz}\cos(kx-\omega t)",
    "7.49": r"c=\sqrt{gH}",
    "7.50": r"\xi\cong-\frac a{kH}\sin(kx_0-\omega t),\ \zeta\cong a\Big(1+\frac {z_0}H\Big)\cos(kx_0-\omega t)",
    "7.51": r"u=\frac{a\omega}{kH}\cos(kx-\omega t),\ w=a\omega\Big(1+\frac zH\Big)\sin(kx-\omega t)",
    "7.52": r"p'=\rho ga\cos(kx-\omega t)=\rho g\eta",
    "7.53": r"p_a-(p)_{z=\eta}=\sigma\frac1R=\sigma\frac{\partial^2\eta/\partial x^2}{[1+(\partial\eta/\partial x)^2]^{3/2}}\cong\sigma\frac{\partial^2\eta}{\partial x^2}",
    "7.54": r"(p)_{z=\eta}=-\sigma\frac{\partial^2\eta}{\partial x^2}",
    "7.55": r"\Big(\frac{\partial\phi}{\partial t}\Big)_{z=0}=\frac\sigma\rho\frac{\partial^2\eta}{\partial x^2}-g\eta",
    "7.56": r"\omega=\sqrt{k\Big(g+\frac{\sigma k^2}{\rho}\Big)\tanh kH}",
    "7.57": r"c=\sqrt{\Big(\frac gk+\frac{\sigma k}{\rho}\Big)\tanh kH}",
    "7.58": r"c_{min}=\Big[\frac{4g\sigma}{\rho}\Big]^{1/4}\ \text{at}\ \lambda_m=2\pi\sqrt{\frac{\sigma}{\rho g}}",
    "7.59": r"c_{min}\approx23.1\ \text{cm/s at}\ \lambda_m\approx1.71\ \text{cm (air–water, our IAPWS numbers)}",
    "7.60": r"c=\sqrt{\frac{2\pi\sigma}{\rho\lambda}}",
    "7.61": r"\eta(x,t)=a\cos[kx+\omega t]",
    "7.62": r"\psi=\frac{2a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}\sin kx\sin\omega t",
    "7.63": r"u=2a\omega\frac{\cosh k(z+H)}{\sinh kH}\sin kx\sin\omega t",
    "7.64": r"kL=(n+1)\pi,\ \lambda=\frac{2L}{n+1}",
    "7.65": r"\omega=\sqrt{\frac{\pi g(n+1)}{L}\tanh\Big[\frac{(n+1)\pi H}{L}\Big]}",
    "7.66": r"\eta=2a\cos\big(\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,t\big)\cos(kx-\omega t)",
    "7.67": r"c_g=\Delta\omega/\Delta k=d\omega/dk",
    "7.68": r"\eta=a(x-c_gt)\cos(kx-\omega t)",
    "7.69": r"c_g=\frac c2\Big[1+\frac{2kH}{\sinh(2kH)}\Big]",
    "7.70": r"c_g=c/2\ \text{(deep)},\ c_g=c\ \text{(shallow)}",
    "7.71": r"F=Ec_g,\ E=\tfrac12\rho ga^2",
    "7.72": r"\eta=a(x,t)\cos[\theta(x,t)]",
    "7.73": r"k\equiv\partial\theta/\partial x,\ \omega\equiv-\partial\theta/\partial t",
    "7.74": r"\partial k/\partial t+\partial\omega/\partial x=0",
    "7.75": r"\frac{\partial k}{\partial t}+c_g\frac{\partial k}{\partial x}=0",
    "7.76": r"\omega=\sqrt{gk\tanh[kH(x)]}=\omega(k,x)",
    "7.77": r"\partial\omega(k,x)/\partial k=c_g",
    "7.78": r"c_g\frac{\partial k}{\partial t}=\frac{\partial\omega}{\partial t}",
    "7.79": r"\frac{\partial\omega}{\partial t}+c_g\frac{\partial\omega}{\partial x}=0",
    "7.80": r"Q^2\Big(\frac1{H_2}-\frac1{H_1}\Big)=\tfrac12g(H_1^2-H_2^2)",
    "7.81": r"\frac{H_2}{H_1}=\tfrac12\big(-1+\sqrt{1+8\mathrm{Fr}_1^2}\big)",
    "7.82": r"\eta=a\cos k(x-ct)+\tfrac12ka^2\cos2k(x-ct)+\tfrac38k^2a^3\cos3k(x-ct)+\ldots",
    "7.83": r"c=\sqrt{\frac gk(1+k^2a^2+\ldots)}",
    "7.84a": r"\frac{dx_p}{dt}=u(x_0,z_0,t)+\xi\Big(\frac{\partial u}{\partial x}\Big)_0+\zeta\Big(\frac{\partial u}{\partial z}\Big)_0+\ldots",
    "7.84b": r"\frac{dz_p}{dt}=w(x_0,z_0,t)+\xi\Big(\frac{\partial w}{\partial x}\Big)_0+\zeta\Big(\frac{\partial w}{\partial z}\Big)_0+\ldots",
    "7.85": r"\bar u_L=a^2\omega ke^{2kz_0}",
    "7.86": r"\bar u_L=a^2\omega k\frac{\cosh 2k(z_0+H)}{2\sinh^2kH}",
    "7.87": r"\frac{\partial\eta}{\partial t}+c_0\frac{\partial\eta}{\partial x}+\frac32c_0\frac\eta H\frac{\partial\eta}{\partial x}+\frac16c_0H^2\frac{\partial^3\eta}{\partial x^3}=0",
    "7.88": r"\eta=a\,\mathrm{sech}^2\Big[\Big(\frac{3a}{4H^3}\Big)^{1/2}(x-ct)\Big],\ c=c_0\Big(1+\frac a{2H}\Big)",
    "7.89": r"\zeta(x,t)=a\exp[i(kx-\omega t)]",
    "7.90": r"\nabla^2\phi_1=0,\ \nabla^2\phi_2=0",
    "7.91": r"\phi_1\to0\ \text{as}\ z\to\infty",
    "7.92": r"\phi_2\to0\ \text{as}\ z\to-\infty",
    "7.93": r"\frac{\partial\phi_1}{\partial z}=\frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}{\partial t}\ \text{at}\ z=0",
    "7.94": r"\rho_1\frac{\partial\phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+\rho_2g\zeta\ \text{at}\ z=0",
    "7.95": r"\omega=\sqrt{gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}}=\varepsilon\sqrt{gk}",
    "7.96": r"E=\tfrac12(\rho_2-\rho_1)ga^2",
    "7.97": r"\phi_2\to0\ \text{at}\ z\to-\infty",
    "7.98": r"\frac{\partial\phi_1}{\partial z}=\frac{\partial\eta}{\partial t}\ \text{at}\ z=0",
    "7.99": r"\frac{\partial\phi_1}{\partial t}+g\eta=0\ \text{at}\ z=0",
    "7.100": r"\frac{\partial\phi_1}{\partial z}=\frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}{\partial t}\ \text{at}\ z=-H",
    "7.101": r"\rho_1\frac{\partial\phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+\rho_2g\zeta\ \text{at}\ z=-H",
    "7.102": r"\eta=ae^{i(kx-\omega t)}",
    "7.103": r"\zeta=be^{i(kx-\omega t)}",
    "7.104": r"\phi_1=(Ae^{kz}+Be^{-kz})e^{i(kx-\omega t)}",
    "7.105": r"\phi_2=Ce^{kz}e^{i(kx-\omega t)}",
    "7.106": r"A=-\frac{ia}2\Big(\frac\omega k+\frac g\omega\Big)",
    "7.107": r"B=\frac{ia}2\Big(\frac\omega k-\frac g\omega\Big)",
    "7.108": r"C=-\frac{ia}2\Big(\frac\omega k+\frac g\omega\Big)-\frac{ia}2\Big(\frac\omega k-\frac g\omega\Big)e^{2kH}",
    "7.109": r"b=\frac a2\Big(1+\frac{gk}{\omega^2}\Big)e^{-kH}+\frac a2\Big(1-\frac{gk}{\omega^2}\Big)e^{kH}",
    "7.110": r"\Big(\frac{\omega^2}{gk}-1\Big)\Big\{\frac{\omega^2}{gk}[\rho_1\sinh kH+\rho_2\cosh kH]-(\rho_2-\rho_1)\sinh kH\Big\}=0",
    "7.111": r"\omega^2=gk",
    "7.112": r"b=ae^{-kH}",
    "7.113": r"\omega^2=\frac{gk(\rho_2-\rho_1)\sinh kH}{\rho_2\cosh kH+\rho_1\sinh kH}",
    "7.114": r"\eta=-\zeta\Big(\frac{\rho_2-\rho_1}{\rho_1}\Big)e^{-kH}",
    "7.115": r"\omega^2=kg\Big(\frac{\rho_2-\rho_1}{\rho_2}\Big)kH",
    "7.116": r"c=[g'H]^{1/2}",
    "7.117": r"g'=g\Big(\frac{\rho_2-\rho_1}{\rho_2}\Big)",
    "7.118": r"\eta=-\zeta\Big(\frac{\rho_2-\rho_1}{\rho_1}\Big)",
    "7.119": r"p'=-\rho_1\frac{\partial\phi_1}{\partial t}=\rho_1g\eta",
    "7.120": r"\frac{\partial u}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial x}",
    "7.121": r"\frac{\partial v}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial y}",
    "7.122": r"\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial z}-\frac{\rho g}{\rho_0}",
    "7.123": r"0=-\frac1{\rho_0}\frac{d\bar p}{dz}-\frac{\bar\rho g}{\rho_0}",
    "7.124": r"p=\bar p(z)+p',\ \rho=\bar\rho(z)+\rho'",
    "7.125": r"\frac{\partial}{\partial t}(\bar\rho+\rho')+\mathbf u\cdot\nabla(\bar\rho+\rho')=0",
    "7.126": r"\frac{\partial\rho'}{\partial t}+w\frac{d\bar\rho}{dz}=0",
    "7.127": r"N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}",
    "7.128": r"\frac{\partial u}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial x}",
    "7.129": r"\frac{\partial v}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial y}",
    "7.130": r"\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial z}-\frac{\rho'g}{\rho_0}",
    "7.131": r"\frac{\partial\rho'}{\partial t}-\frac{N^2\rho_0}gw=0",
    "7.132": r"\frac1{\rho_0}\nabla_H^2p'=\frac{\partial^2w}{\partial z\,\partial t}",
    "7.133": r"\frac1{\rho_0}\frac{\partial^2p'}{\partial t\,\partial z}=-\frac{\partial^2w}{\partial t^2}-N^2w",
    "7.134": r"\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0",
    "7.135": r"\omega=\omega(k,l,m)=\omega(\mathbf K)",
    "7.136": r"w=w_0e^{i(kx+ly+mz-\omega t)}",
    "7.137": r"\omega^2=\frac{k^2+l^2}{k^2+l^2+m^2}N^2",
    "7.138": r"\omega=\frac{kN}{\sqrt{k^2+m^2}}=\frac{kN}K\ (k>0)",
    "7.139": r"\omega=N\cos\theta",
    "7.140": r"u=u_0e^{i(kx+ly+mz-\omega t)}",
    "7.141": r"\mathbf K\cdot\mathbf u=0",
    "7.142": r"\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0",
    "7.143": r"\mathbf c_g=\mathbf e_x\frac{\partial\omega}{\partial k}+\mathbf e_y\frac{\partial\omega}{\partial l}+\mathbf e_z\frac{\partial\omega}{\partial m}",
    "7.144": r"\mathbf c=\frac\omega{K^2}(k\mathbf e_x+m\mathbf e_z)",
    "7.145": r"\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)\ (k>0)",
    "7.146": r"\mathbf c_g\cdot\mathbf c=0",
    "7.147": r"\frac{\partial}{\partial t}\Big[\tfrac12\rho_0(u^2+v^2+w^2)\Big]+g\rho'w+\nabla\cdot(p'\mathbf u)=0",
    "7.148": r"\frac{\partial E_p}{\partial t}=g\rho'w=\frac{\partial}{\partial t}\Big[\frac{g^2\rho'^2}{2\rho_0N^2}\Big]",
    "7.149": r"\rho'=\frac{N^2\rho_0\zeta}{g}",
    "7.150": r"E_p=\frac{g^2\rho'^2}{2\rho_0N^2}=\tfrac12N^2\rho_0\zeta^2",
    "7.151": r"E_p=\tfrac14(\rho_2-\rho_1)ga^2",
    "7.152": r"N^2=\frac g{\rho_0}(\rho_2-\rho_1)\delta(z)",
    "7.153": r"p'=-\frac{\omega m\rho_0}{k^2}\hat we^{i\theta},\ \rho'=\frac{iN^2\rho_0}{\omega g}\hat we^{i\theta},\ u=-\frac mk\hat we^{i\theta}",
    "7.154": r"E_k=\tfrac14\rho_0\Big(\frac{m^2}{k^2}+1\Big)\hat w^2",
    "7.155": r"E_p=\frac{N^2\rho_0}{4\omega^2}\hat w^2",
    "7.156": r"E_k=E_p",
    "7.157": r"E=\tfrac12\rho_0\Big(\frac{m^2}{k^2}+1\Big)\hat w^2",
    "7.158": r"\mathbf F=\overline{p'\mathbf u}=\frac{\rho_0\omega m\hat w^2}{2k^2}\Big(\mathbf e_x\frac mk-\mathbf e_z\Big)",
    "7.159": r"\mathbf F=\mathbf c_gE",
    # earlier chapters (from their notebooks)
    "1.5": r"\Delta p=\sigma\Big(\frac1{R_1}+\frac1{R_2}\Big)",
    "1.29": r"N^2=-\frac{g}{\rho}\Big(\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\Big)",
    "3.8": r"d\mathbf r/dt=\mathbf u(\mathbf r,t)",
    "4.9": r"\frac{D\rho}{Dt}\equiv\frac{\partial\rho}{\partial t}+\mathbf u\cdot\nabla\rho=0",
    "4.10": r"\nabla\cdot\mathbf u=0",
    "4.17": r"\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA",
    "4.56": r"\rho\frac{D}{Dt}\big(\tfrac12u_j^2\big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}",
    "4.83": r"\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+\frac p\rho=\text{const}",
    "4.90": r"\frac{\partial\eta}{\partial t}+(\mathbf u_s\cdot\nabla)\eta=0\ \text{on}\ \eta=0",
    "4.91": r"\frac{D\eta}{Dt}=\frac{\partial\eta}{\partial t}+(\mathbf u\cdot\nabla)\eta=0\ \text{on}\ \eta=0",
    "4.104": r"\mathrm{Fr}=\frac{U}{\sqrt{gl}}",
    "6.16": r"\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{surface}}",
}
EQ["7.34"] = EQ["7.34a"] + r",\ " + EQ["7.34b"]
EQ["7.35"] = EQ["7.35a"] + r",\ " + EQ["7.35b"]
EQ["7.84"] = EQ["7.84a"]


def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ *({num})*"


def EE(*nums: str) -> str:
    """Several equations, each with its number, joined for prose."""
    return " · ".join(E(n) for n in nums)


_EQ_REF = re.compile(r"\(((?:\d)\.\d+[a-d]?)\)")


def _plain_eq_follows(after: str) -> bool:
    """Is the equation already written in plain symbols right after its number ("(7.10) u = ∂φ/∂x", "(7.28): ω² = gk
    tanh kH")? Then inserting its LaTeX again would only duplicate it (lesson review round 1, Should-fix 2)."""
    a = after.lstrip(" ,:*")[:40]
    pos = a.find("=")
    if pos <= 0 or pos > 14:
        return False
    head, rhs = a[:pos], a[pos + 1:].lstrip()
    if re.search(r"[A-Za-z]{4,}|[.;()→–]", head) or rhs[:1].isdigit():   # prose, a range, or a computed number
        return False
    stop = {"is", "of", "to", "in", "at", "as", "be", "by", "the", "and", "so", "it", "on", "for", "with"}
    return not any(w.lower() in stop for w in head.split())


def _maths_before(before: str) -> bool:
    """Does maths end right before the number (possibly across a line break or an italic star)?"""
    return before.rstrip(" *\n").endswith("$")


def show_eqs(text: str) -> str:
    """Write the equation next to the first mention of a book equation number "(7.28)" in a text that names it without
    writing it — the house rule "show the equation, not just its number"."""
    done: set[str] = set()

    def rep(m):
        n = m.group(1)
        nxt = text[m.end():m.end() + 16]
        if (n in done or n not in EQ or EQ[n] in text or nxt.lstrip(", :").startswith("$")
                or _plain_eq_follows(text[m.end():]) or _maths_before(text[:m.start()])
                or nxt.startswith(" — the line above") or nxt.startswith(" — the result")
                or nxt[:1] in "/–-)" or text[max(0, m.start() - 1):m.start()] in ("(", "–", "-")):
            return m.group(0)
        done.add(n)
        return f"({n}), ${EQ[n]}$"
    return _EQ_REF.sub(rep, text)


# ---------------------------------------------------------------------------------------------------------------------
# Part F reader: the derivations, word for word (the ch03–ch06 parser; step fields may carry an explainer qualifier)
# ---------------------------------------------------------------------------------------------------------------------
def _join(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(x.strip() for x in lines)).strip()


def _abs_plain(s: str) -> str:
    r"""Part F writes absolute values as \lvert … \rvert everywhere; outside $…$ they become plain |…|."""
    parts = re.split(r"(\$[^$]+\$)", s)
    return "".join(q if q.startswith("$") else re.sub(r"\\rvert", "|", re.sub(r"\\lvert\s?", "|", q)) for q in parts)


def _tex_escape(s: str) -> str:
    return (s.replace("\\", r"\backslash ").replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")
            .replace("_", r"\_").replace("{", r"\{").replace("}", r"\}"))


def display_tex(s: str, width: int = 78) -> str:
    """A Part F line (maths in `$…$`, possibly with words around it) as the body of one display-math block."""
    s = s.strip()
    parts = re.split(r"(\$[^$]+\$)", s)
    parts = [p for p in parts if p.strip()]
    if len(parts) == 1 and parts[0].startswith("$"):
        return parts[0][1:-1]
    tokens: list[tuple[str, str]] = []
    for p in parts:
        if p.startswith("$"):
            tokens.append(("m", p[1:-1]))
        else:   # the antiphase mode
            tokens += [("w", w.replace("**", "").replace("*", "").replace("`", "").replace("\\|", "|")) for w in p.split()]
    lines, cur, n = [], [], 0
    for kind, val in tokens:
        size = len(val) if kind == "w" else max(4, len(val) // 3)
        if cur and n + size > width:
            lines.append(cur)
            cur, n = [], 0
        cur.append((kind, val))
        n += size + 1
    if cur:
        lines.append(cur)

    def render(line):
        out, words = [], []
        for kind, val in line:
            if kind == "w":
                words.append(val)
                continue
            if words:
                out.append(r"\text{" + _tex_escape(" ".join(words)) + r"}")
                words = []
            out.append(val)
        if words:
            out.append(r"\text{" + _tex_escape(" ".join(words)) + r"}")
        return r"\ ".join(out)

    if len(lines) == 1:
        return render(lines[0])
    return r"\begin{array}{l}" + r" \\ ".join(render(ln) for ln in lines) + r"\end{array}"


def _split_words(s: str) -> tuple[str, str]:
    """'tex — *in words:* plain' → (tex, plain)."""
    if "*in words:*" in s:
        a, b = s.split("*in words:*", 1)
        return a.strip().rstrip("—").strip(), b.strip()
    return s.strip(), ""


_SELF = re.compile(r"(This is (?:the book's )?|this is the book's |is the book's |it is the book's |This is exactly )\(((?:\d)\.\d+[a-d]?)\)")


def _self_ref(text: str) -> str:
    """"this is the book's (7.133)" inside a step's *why* points at the line just displayed: say so."""
    return _SELF.sub(lambda m: f"{m.group(1)}({m.group(2)}) — the line above", text)


def eq_display(nums: list[str], words: str = "") -> str:
    """Display-math body writing several book equations with their numbers (for Starts and Results that only name
    numbers in Part F)."""
    body = r",\qquad ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in nums)
    if len(nums) > 2:
        body = r"\begin{array}{l}" + r" \\ ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in nums) + r"\end{array}"
    return body + (r"\ \ \text{" + _tex_escape(words) + "}" if words else "")


def part_f() -> dict[str, dict]:
    """Parse Part F of the design into {D01: dict(title, goal, start, plan, tools, assumptions, steps, result, check,
    meaning, traps)}."""
    text = (ROOT / "analysis" / "ch07_design.md").read_text(encoding="utf-8")
    part = text.split("## Part F", 1)[1]
    chunks = re.split(r"^### (D\d\d) · ", part, flags=re.M)
    out: dict[str, dict] = {}
    for key, body in zip(chunks[1::2], chunks[2::2]):
        lines = body.splitlines()
        title = _abs_plain(re.split(r" — ★", lines[0])[0].strip())
        fields: dict[str, list[str]] = {}
        steps: list[list[str]] = []
        cur = None
        for ln in lines[1:]:
            if ln.startswith("---") or ln.startswith("## "):
                break
            m = re.match(r"^- \*\*(.+?)\*\*\s*(.*)$", ln)
            if m:
                cur = m.group(1).strip().rstrip(".").strip()
                fields[cur] = [m.group(2)]
                continue
            if cur == "Steps":
                ms = re.match(r"^\s{2}(\d+)\. (.*)$", ln)
                if ms:
                    steps.append([ms.group(2)])
                    continue
                if steps and ln.strip():
                    steps[-1].append(ln)
                    continue
            if cur and ln.strip():
                fields[cur].append(ln)
        f = {k: _abs_plain(_join(v)) for k, v in fields.items() if not k.startswith("sympy")}
        parsed = []
        for st in steps:
            s = _abs_plain(_join(st))
            bits = re.split(r"(?:^|\s·\s)\*(did|tex|why|plain|live|set|watch)(?:\s*\([^)]*\))?:\*\s*", s)
            d = {name: val.strip() for name, val in zip(bits[1::2], bits[2::2])}
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=_self_ref(d["why"]),
                               plain=_self_ref(d["plain"])))
        start_tex, start_plain = _split_words(f["Start"])
        res_tex, res_plain = _split_words(f["Result"])
        res_tex = re.sub(r"\s*\*?\((?:\d)\.\d+[a-d]?\)\*?\.?\s*$", "", res_tex) if res_tex.count("$") >= 2 else res_tex
        plan = [p.strip() for p in re.split(r"\(\d+\)\s*", f["Plan"]) if p.strip()]
        tools = [t.strip() for t in f["Tools"].split(" · ") if t.strip()]
        out[key] = dict(title=title, goal=f["Goal"], start=(display_tex(start_tex), start_plain),
                        plan=plan, tools=tools, assumptions=f.get("Assumptions", ""), steps=parsed,
                        result=(display_tex(res_tex), res_plain), check=f.get("Check", ""),
                        meaning=f.get("What it means", ""), traps=f.get("Traps", ""))
    return out


PF = part_f()


def pf_sub(key: str, field: str, old: str, new: str) -> None:
    """Edit one Part F field after parsing (pointers to the cells that run a check, verification decisions); fails
    loudly if the design text changed. ``field`` = goal/check/meaning/traps/assumptions, ``tools``, or
    ``stepN.did|why|plain``."""
    d = PF[key]
    if field.startswith("step"):
        i, part = field[4:].split(".")
        st = d["steps"][int(i) - 1]
        assert old in st[part], (key, field, old)
        st[part] = st[part].replace(old, new)
        return
    if field == "tools":
        hit = [j for j, t in enumerate(d["tools"]) if old in t]
        assert hit, (key, field, old)
        d["tools"][hit[0]] = d["tools"][hit[0]].replace(old, new)
        return
    assert old in d[field], (key, field, old)
    d[field] = d[field].replace(old, new)


if "--dump" in sys.argv:
    for k, d in PF.items():
        print(f"=== {k} {d['title']}  ({len(d['steps'])} steps)")
        print("  START", d["start"][0][:150], "|", d["start"][1][:80])
        for i, s in enumerate(d["steps"], 1):
            print(f"  {i}. {s['did']} :: {s['tex'][:110]}  || why={len(s['why'].split())}w")
        print("  RESULT", d["result"][0][:150], "|", d["result"][1][:80])
        print("  CHECK", d["check"])
        print("  TOOLS", d["tools"])
    sys.exit(0)

# ---------------------------------------------------------------------------------------------------------------------
# Every "check" names a cell of this notebook that really runs it; Starts/Results that only name numbers show the
# equations; the verification decisions (reports/ch07_verification.md) are applied.
# ---------------------------------------------------------------------------------------------------------------------
pf_sub("D02", "check", "(C02 code).", "(the code cell after the tiny example of C02 prints it: about 5×10⁻³ m/s for a = 5 cm, k = 1 rad/m).")
pf_sub("D08", "check", "✓ (E1 preset)", "✓ (`ch07.depth_regime` returns this error)")
pf_sub("D11", "check", "✓ (C05 code)", "✓ (the code cell after the tiny example of C05)")
pf_sub("D13", "check", "(C06 code with `quad`)", "(the code cell after the tiny example of C06, also with `quad`)")
pf_sub("D17", "check", "(C08 code)", "(the code cell after the tiny example of C08)")
pf_sub("D22", "check", "(C10 code)", "(the code cell after the tiny example of C10)")
pf_sub("D27", "check", "drift at (7.85) within 2 % at ka = 0.05 ✓",
       "drift at (7.85) within about 1 % at ka = 0.02 and 3–5 % at ka = 0.05 — the gap is itself O(ka), the size of the "
       "terms the derivation dropped ✓")
pf_sub("D16", "meaning", "(≈ 17.8 cm/s)", "(0.1776 m/s at λ ≈ 4.35 cm for clean water, computed with `min_group_velocity` in C09)")
pf_sub("D24", "step1.why", "Hamilton's equations for rays)",
       r"Hamilton's equations for rays — the same pair as Hamilton's equations of mechanics, $dx/dt=\partial\mathcal H/\partial p$, "
       r"$dp/dt=-\partial\mathcal H/\partial x$, with the ray's point x as the position, the wavenumber k as the momentum and "
       r"ω as the energy, the Hamiltonian $\mathcal H$ (a script letter, not the depth H))")
pf_sub("D34", "check", "(C15 code: three equal values)", "(the code cell after the tiny example of C15: three equal values)")
PF["D09"]["start"] = (eq_display(["7.29", "7.35a", "7.35b", "7.27", "7.31"]), "the speed, orbits, velocities and pressure of the general-depth wave")
PF["D10"]["start"] = (eq_display(["7.32", "7.27"]), "each parcel moves with the local velocity of the linear wave")
PF["D11"]["start"] = (eq_display(["7.35a", "7.35b"]), "the linear excursions of D10")
PF["D15"]["result"] = (eq_display(["7.53", "7.54", "7.55", "7.56", "7.57"]), "tension adds σk²/ρ to g")
PF["D17"]["result"] = (eq_display(["7.62", "7.63"], "and the surface 2a cos kx cos ωt"), "a standing wave, its streamlines and its horizontal velocity")
PF["D19"]["result"] = (eq_display(["7.66", "7.67"]), "fast carrier at c, slow envelope at the group velocity")
PF["D29"]["result"] = (eq_display(["7.106", "7.107", "7.108", "7.109"]), "every constant in terms of a and ω")
PF["D30"]["result"] = (eq_display(["7.110"]), "two factors: two modes")
PF["D31"]["start"] = (eq_display(["7.110", "7.109"]), "the two-layer dispersion relation and the interface amplitude")
PF["D32"]["result"] = (eq_display(["7.126", "7.128", "7.129", "7.130", "7.131", "4.10"]), "five linear equations for u, v, w, p′, ρ′")
PF["D33"]["start"] = (eq_display(["4.10", "7.128", "7.129", "7.130", "7.131"]), "the linear internal-wave equations of D32")
PF["D37"]["result"] = (eq_display(["7.153", "7.158", "7.159"]), "every field through w; the energy flux is group velocity times energy")
# Starts and steps that only named equation numbers now write the equations (lesson review round 1, Must-fix 6)
PF["D30"]["start"] = (eq_display(["7.101", "7.106", "7.107", "7.108", "7.109"]),
                      "the interface pressure condition, with the constants of D29 in terms of a and ω")
PF["D32"]["start"] = (r"\begin{array}{l}" + r" \\ ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in ("4.9", "4.10", "7.120", "7.121", "7.122", "7.123"))
                      + r" \\ p=\bar p(z)+p',\quad\rho=\bar\rho(z)+\rho'\ \text{(7.124)}\end{array}",
                      "mass conservation, incompressibility, the momentum equations, the resting (hydrostatic) state and the split into base plus perturbation")
PF["D32"]["steps"][8]["tex"] = (r"\begin{array}{l}\mathbf u\cdot\nabla\mathbf u\ \text{already dropped in} \\ "
                                + r" \\ ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in ("7.120", "7.121", "7.122")) + r"\end{array}")
PF["D32"]["steps"][9]["tex"] = (r"\begin{array}{l}" + r" \\ ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in ("4.10", "7.128", "7.129", "7.130", "7.131"))
                                + r" \\ \text{5 equations for}\ u,v,w,p',\rho'\end{array}")
PF["D35"]["start"] = (eq_display(["7.140", "4.10"]), "a plane wave in u (and the same form for v and w), and incompressibility")
PF["D37"]["start"] = (r"\begin{array}{l}[u,w,p',\rho']=[\hat u,\hat w,\hat p,\hat\rho]e^{i(kx+mz-\omega t)}\ \text{with} \\ "
                      + r" \\ ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in ("7.128", "7.131", "7.132")) + r"\end{array}",
                      "plane-wave fields in the horizontal momentum, density and pressure equations of D32–D33")
# Results whose exponentials sat inside \text{} (review Should-fix 4) and D28's start written out in full
PF["D08"]["result"] = (r"\begin{array}{l}c=\sqrt{g/k}\ \text{(7.45)},\qquad " + EQ["7.46"] + r"\ \text{(7.46): circles of radius}\ ae^{kz_0} \\ "
                       + EQ["7.47"] + r"\ \text{(7.47)},\qquad " + EQ["7.48"] + r"\ \text{(7.48)}\end{array}", PF["D08"]["result"][1])
PF["D28"]["start"] = (eq_display(["7.89", "7.90", "7.91", "7.92", "7.93", "7.94"]), PF["D28"]["start"][1] or "the interface problem: a wave on the interface, Laplace in each fluid, decay far away, and the kinematic and pressure conditions at z = 0")
PF["D28"]["result"] = (r"\omega=\varepsilon\sqrt{gk},\ \ \varepsilon^2=\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\ \text{(7.95)};\qquad u_1=-\omega ae^{-kz}e^{i\theta},\ \ u_2=\omega ae^{kz}e^{i\theta}\ \text{(note N107)}",
                       PF["D28"]["result"][1])
PF["D31"]["result"] = (r"\begin{array}{l}\text{barotropic:}\ \omega^2=gk,\ \ b=ae^{-kH} \\ \text{baroclinic:}\ " + EQ["7.113"]
                       + r"\ \text{(7.113)},\ \ \frac\eta\zeta=-\frac{\rho_2-\rho_1}{\rho_1}e^{-kH} \\ \text{long waves:}\ c=\sqrt{g'H}\end{array}",
                       PF["D31"]["result"][1])
# D14 step 5 names the first form of the pressure: show that form (review Should-fix 3)
PF["D14"]["steps"][4]["why"] = (r"The first form of the pressure, $p'=-\rho\frac{\partial\phi}{\partial t}=\rho\frac{a\omega^2}{k}\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t)$ "
                                r"(the linearised Bernoulli equation applied to the potential of D05), and u from (7.27) (the book's 'u from (7.28)' means (7.27)).")
# D07 step 4: say that the bracket is c² up to a positive constant (review Should-fix 7)
pf_sub("D07", "step4.why", "Product and chain rules:",
       "By $c=\\sqrt{\\frac gk\\tanh kH}$ (7.29) with k = 2π/λ, c² = (g/2π)·λ tanh(2πH/λ): the bracket is c² up to the positive "
       "constant g/2π, which cannot change the sign of the slope. Product and chain rules:")
# D26 step 7: no pipeline jargon (review Should-fix 8)
pf_sub("D26", "step7.why", " (The extracted text drops the minus sign; the page has it.)", " (Some copies lose this minus sign; the energy argument below needs it.)")
# D03/D04 size the terms from η alone, so C03's formulas are not needed yet (review Should-fix 5)
pf_sub("D03", "step1.why", "velocities are aω (C03's (7.27))",
       "the water at the surface has to keep up with a surface that rises at ∂η/∂t ∼ aω, so the velocities are of size aω "
       "too (no formula from C03 needed)")
pf_sub("D04", "step3.why", "velocities are aω and φ ∼ aω/k (7.26)",
       "velocities are aω (the surface rises at that rate), and a potential whose gradient aω changes over a distance "
       "1/k has size φ ∼ aω/k")
PF["D31"]["title"] = r"Roots, mode shapes and long waves: $c=\sqrt{g'H}$"   # the range of numbers is shown in the heading equation
# lesson review round 2: no exercise labels or book quotes in prose (rule 9), one copy of each equation, the second
# particle equation written in maths (not \text{z\_p}), and the Hamiltonian letter distinct from the depth H
pf_sub("D21", "check", " (Exercise 7.9, `group_velocity(g=0)`)", " (`group_velocity(g=0)` returns it)")
pf_sub("D30", "goal", ' (The book: "after some algebraic manipulations", Exercise 7.19; written out here.)',
       " The book skips this algebra; we write it out move by move.")
pf_sub("D31", "step2.why", "In (7.109) gk/ω² = 1: 1 + 1/s = 2 and 1 − 1/s = 0.",
       f"Put gk/ω² = 1 into the interface amplitude ${EQ['7.109']}$ (7.109): the first bracket becomes 1 + 1 = 2 and the "
       "second 1 − 1 = 0.")
PF["D10"]["steps"][0]["tex"] = (r"\begin{array}{l}\frac{dx_p}{dt}=a\omega\frac{\cosh k(z_p+H)}{\sinh kH}\cos(kx_p-\omega t) \\ "
                                r"\frac{dz_p}{dt}=a\omega\frac{\sinh k(z_p+H)}{\sinh kH}\sin(kx_p-\omega t)\ \text{(7.33)}\end{array}")


# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch06.py)
# ---------------------------------------------------------------------------------------------------------------------
def core(cid: str, title: str, question: str, eqs: tuple[str, ...] = ()) -> None:
    """A CORE block heading with its id; ``eqs`` = the book equations named in the title, written out under it."""
    q = textwrap.dedent(question).strip()
    if eqs:
        q += "\n\n*In one line:* " + " · ".join(E(n) for n in eqs)
    nb.core(cid, f"{title} `{cid}`", question=q)


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` is exactly the Part E concept text (the ledger check matches it)."""
    nb.primer(term, f"`{pid}` · {textwrap.dedent(text).strip()}", code=textwrap.dedent(code).strip("\n") if code else None)


def remind(items: list[tuple[str, str]], lead: str = "") -> None:
    """One 🔁 cell reminding several tools primed in earlier chapters (knowledge/primers.md) — one sentence each."""
    body = "\n".join(f"- **{c}** — {textwrap.dedent(t).strip()}" for c, t in items)
    head = "> 🔁 **Tools from earlier chapters used here** (see `knowledge/primers.md`)"
    nb.md(f"{head}{(' — ' + lead) if lead else ''}\n\n{body}", tags=["primer"])
    for c, _ in items:
        nb.primers.append(f"{c} (reminder)")


def note(nid: str, title: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note with its curation id(s)."""
    body = textwrap.dedent(text).strip()
    sep = " " if body.startswith("—") or not body else " — "
    nb.note(f"**{title}** `{nid}`{sep}{body}", equation=equation, ref=ref)


def _tidy_check(s: str) -> str:
    return re.sub(r"\s*`check_src` sketch.*$", "", s, flags=re.S)


def _title_no_numbers(title: str) -> str:
    """Drop bare equation numbers from a derivation title (the heading shows the key equation itself instead)."""
    num = r"\(\d\.\d+[a-d]?(?:,\s*(?:\d\.\d+)?[a-d]?)*\)"
    t_ = re.sub(r"\s*(?:—\s*)?" + num + r"(?:\s*(?:→|–|,|and|or)\s*" + num + r")*", "", title)
    t_ = re.sub(r"\s+([,:)])", r"\1", t_)
    t_ = re.sub(r"\s*[—:,]\s*$", "", t_)
    t_ = re.sub(r":\s*→\s*", ": ", t_)
    return re.sub(r"\s{2,}", " ", t_).strip()


def D(key: str, ref: str = "", check_src: str | None = None, extra_check: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then its traps as a ⚠️ callout."""
    d = PF[key]
    if ref in EQ:
        ref = f"{ref}: ${EQ[ref]}$"                     # the heading shows the key equation, not only its number
    S = show_eqs
    check = S(_tidy_check(d["check"])) + (f" {extra_check}" if extra_check else "")
    goal = S(d["goal"]) + (f"\n\n**Assumptions.** {S(d['assumptions'])}" if d["assumptions"] else "")
    steps = [dict(st, why=S(st["why"]), plain=S(st["plain"])) for st in d["steps"]]
    nb.derivation(key, f"{_title_no_numbers(d['title'])} `{key}`", ref=ref, goal=goal, start=(d["start"][0], S(d["start"][1])),
                  plan=[S(x) for x in d["plan"]], uses=[S(x) for x in d["tools"]], steps=steps,
                  result=(d["result"][0], S(d["result"][1])), interpret=S(d["meaning"]), check=check,
                  check_src=textwrap.dedent(check_src).strip("\n") if check_src else None)
    if d["traps"]:
        nb.md(f"> ⚠️ **Common confusion (traps in `{key}`):** {S(d['traps'])}")


def see_read_change(see: str, read: str, change: str) -> None:
    nb.figure_notes(see, read, change)


def whatif(text: str) -> None:
    nb.md(f"**What would change if…** {textwrap.dedent(text).strip()}")


def idea(sketch: str, words: str = "") -> None:
    """The idea block: an ASCII sketch plus one or two sentences; the book equations the sketch names by number are
    written out underneath it (the sketch itself is plain text)."""
    sk = textwrap.dedent(sketch).strip("\n")
    nums = [n for n in dict.fromkeys(re.findall(r"\((\d\.\d+[a-d]?)\)", sk)) if n in EQ]
    eqs = ("\n\n*The numbered equations in the sketch:* " + " · ".join(E(n) for n in nums)) if nums else ""
    nb.md("#### The idea\n\n```\n" + sk + "\n```" + eqs + (f"\n\n{textwrap.dedent(words).strip()}" if words else ""))


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + textwrap.dedent(text).strip())


def explainer(slug: str, heading: str, why_static: str, body: str, tries: list[str]) -> None:
    """An embedded explainer, opening with the "why interactive rather than a static figure" sentence (ch04 lesson)."""
    why = f"**Why interactive rather than a static figure:** {textwrap.dedent(why_static).strip()}"
    if body:
        why += "\n\n" + textwrap.dedent(body).strip()
    nb.explainer(slug, heading=heading, why=show_eqs(why), tries=[show_eqs(t) for t in tries])


def confusion(text: str) -> None:
    nb.md("> ⚠️ **Common confusion:** " + textwrap.dedent(text).strip())


# ---------------------------------------------------------------------------------------------------------------------
# Final passes: every book equation named by number in prose gets the equation written next to it.
# Only parenthesised groups count as equation mentions — "(7.28)", "(7.116, 7.117)", "(7.34a, b)", "(Eq. 7.4)" — so a
# number such as "7.8 m/s" is never mistaken for Eq. (7.8).
# ---------------------------------------------------------------------------------------------------------------------
_PROSE_SPLIT = re.compile(r"(```.*?```|\$\$.*?\$\$|\$[^$]+\$|`[^`\n]*`)", re.S)
_GROUP = re.compile(r"\((?:Eqs?\.\s*)?(\d\.\d+[a-d]?(?:\s*(?:,|–|-|and)\s*(?:\d\.\d+[a-d]?|[a-d]))*)\)")


def _labels(group: str) -> list[str]:
    out, base = [], ""
    for tok in re.split(r"\s*(?:,|–|-|and)\s*", group):
        if re.fullmatch(r"\d\.\d+[a-d]?", tok):
            out.append(tok)
            base = re.sub(r"[a-d]$", "", tok)
        elif re.fullmatch(r"[a-d]", tok) and base:
            out.append(base + tok)
    return [n for n in out if n in EQ]


def _shown(n: str, unit: str) -> bool:
    """Is equation n written in this unit (its LaTeX, or a display tagged with its number)?"""
    return EQ[n] in unit or ("\\text{(" + n + ")}") in unit or ("(" + n + ")}") in unit


def _mentions(unit: str) -> list[tuple[int, int, list[str], bool]]:
    """(start, end, labels, adjacent_to_maths) of every parenthesised equation group in the prose of a unit (outside
    maths, code, headings and derivation step titles)."""
    out, pos = [], 0
    parts = _PROSE_SPLIT.split(unit)
    for k, seg in enumerate(parts):
        if k % 2 == 0:
            for m in _GROUP.finditer(seg):
                labs = _labels(m.group(1))
                if not labs:
                    continue
                before = seg[:m.start()]
                line = before[before.rfind("\n") + 1:]
                if line.lstrip().startswith("#") or line.startswith("**Step "):
                    continue
                s0, e0 = pos + m.start(), pos + m.end()
                adj = _maths_before(unit[:s0]) or unit[e0:].lstrip(" *,:").startswith("$")
                plain = len(labs) == 1 and _plain_eq_follows(unit[e0:])    # written in plain symbols right after
                out.append((s0, e0, labs, adj, plain))
        pos += len(seg)
    merged: list = []                                    # a range "(7.120)–(7.122)" becomes one mention of 7.120…7.122
    for s0, e0, labs, adj, plain in out:
        if merged and s0 == merged[-1][1] + 1 and unit[merged[-1][1]] in "–-" and len(merged[-1][2]) == 1 and len(labs) == 1:
            a_, b_ = merged[-1][2][0], labs[0]
            ch_a, n_a = a_.split("."); ch_b, n_b = b_.split(".")
            span = [a_, b_]
            if ch_a == ch_b and n_a.isdigit() and n_b.isdigit() and 0 < int(n_b) - int(n_a) <= 4:
                span = [f"{ch_a}.{j}" for j in range(int(n_a), int(n_b) + 1) if f"{ch_a}.{j}" in EQ]
            merged[-1] = (merged[-1][0], e0, span, merged[-1][3] or adj, False)   # a range is never "written in plain"
        else:
            merged.append((s0, e0, labs, adj, plain))
    return [(a_, b_, c_, d_ or e_) for a_, b_, c_, d_, e_ in merged]


def finalize_equations() -> int:
    """Every equation group still named without its equation gets the equation(s) written right after it (first mention
    in each unit; a derivation's steps are separate units)."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        units = c.source.split("\n---\n")
        for u, src in enumerate(units):
            done: set[str] = set()
            new, last = [], 0
            for s, e, labs, adj in _mentions(src):
                miss = [n for n in labs if not _shown(n, src) and n not in done]
                done.update(labs)
                if adj or not miss:
                    continue
                new.append(src[last:e] + ", $" + r",\ ".join(EQ[n] for n in miss) + "$")
                last = e
            if new:
                units[u] = "".join(new) + src[last:]
        out = "\n---\n".join(units)
        if out != c.source:
            c.source = out
            changed += 1
    return changed


_UNI_TEX = {"−": "-", "₀": "_0", "₁": "_1", "₂": "_2", "θ": r"\theta ", "Θ": r"\Theta ", "λ": r"\lambda ",
            "ω": r"\omega ", "±": r"\pm ", "∓": r"\mp ", "ρ": r"\rho ", "σ": r"\sigma "}


def _exp_to_maths(m: re.Match) -> str:
    """'e^{−2kH}' written in plain text (Part F why/in-words lines) → the maths $e^{-2kH}$."""
    body = "".join(_UNI_TEX.get(ch, ch) for ch in m.group(0))        # a run such as e^{−2kH}e^{−kz} stays one maths span
    return "$" + re.sub(r"\s+\}", "}", body) + "$"


def tidy_raw_tex() -> int:
    """Plain-text exponents copied from Part F ('e^{kz₀}', '(4gσ/ρ)^{1/4}') showed their braces on the page (lesson
    review round 2, Should-fix 2): an exponential becomes inline maths, a bracket's power becomes '^(1/4)'."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        parts = _PROSE_SPLIT.split(c.source)
        for k in range(0, len(parts), 2):
            seg = re.sub(r"(?:(?<![A-Za-z\\])[A-Za-z0-9])?(?<!\\)(?:e\^\{[^{}$]*\})+", _exp_to_maths, parts[k])   # with a one-letter factor (Ae^{kz})
            parts[k] = re.sub(r"\)\^\{([^{}$]*)\}", r")^(\1)", seg)
            if k and parts[k].startswith("$") and parts[k - 1].endswith("$"):
                parts[k] = " " + parts[k]                  # never glue two inline maths into a "$$"
            if k + 1 < len(parts) and parts[k].endswith("$") and parts[k + 1].startswith("$"):
                parts[k] += " "
        out = "".join(parts)
        assert out.count("$$") == c.source.count("$$"), c.source[:120]   # no accidental display maths
        if out != c.source:
            c.source = out
            changed += 1
    return changed


def self_check_prose() -> list[str]:
    """No TeX command or TeX exponent outside maths or code, no escaped underscore inside \\text{}, and no equation
    group named without being shown."""
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "markdown":
            continue
        prose = "".join(seg for k, seg in enumerate(_PROSE_SPLIT.split(c.source)) if k % 2 == 0)
        cmds = re.findall(r"\\[A-Za-z]+|[\^_]\{", prose) + re.findall(r"\\text\{[^}]*\\_", c.source)
        if cmds:
            bad.append(f"cell {i}: TeX outside maths {sorted(set(cmds))[:5]}: {c.source[:80]!r}")
        real = []
        for unit in c.source.split("\n---\n"):
            seen: set[str] = set()
            for _s, _e, labs, adj in _mentions(unit):
                for n in labs:
                    if not (n in seen or adj or _shown(n, unit)):
                        real.append(n)
                    seen.add(n)
        if real:
            bad.append(f"cell {i}: equations named but not shown {sorted(set(real))}: {c.source[:80]!r}")
    return bad


def self_check_numbers() -> list[str]:
    """The coverage tool's rule 8, run at build time: a markdown cell that names a book equation number must contain
    maths. Returns the offending cells (the builder fails on any)."""
    from coverage_check import EQ_REF
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type == "markdown":
            nums = sorted(set(EQ_REF.findall(c.source)))
            if nums and "$" not in c.source and "\\begin{" not in c.source:
                bad.append(f"cell {i}: {nums}: {c.source[:120]!r}")
    return bad


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=r"""
Drop a stone in a pond, watch swell roll onto a beach, or feel a lake slosh after a storm: the water surface wants to be
flat, gravity pulls every bump back, and the water's inertia overshoots — so the surface oscillates and the oscillation
travels. This chapter turns that picture into one formula, the dispersion relation $\omega=\sqrt{gk\tanh kH}$
*(Eq. 7.28)*, and then reads everything off it: how fast crests move and why long waves outrun short ones, how the water
under a wave goes round in orbits, how much energy a wave holds and at which speed that energy travels (the group
velocity $c_g=d\omega/dk$ *(Eq. 7.67)*, not the crest speed), why waves turn toward a beach, what happens when they
steepen into a hydraulic jump, and how a small density step or a continuous stratification carries slow, huge internal
waves whose energy leaves at right angles to their crests, $\mathbf c_g\cdot\mathbf c=0$ *(Eq. 7.146)*. Ocean swell,
tides, tsunamis, thermocline waves and the internal waves of the atmosphere are all here.
""",
    roadmap=[
        r"§7.1 the sinusoidal wave and its vocabulary, $c=\omega/k=\lambda\nu$ *(7.4)* (C01)",
        r"§7.2 the linear free-surface problem, conditions moved to z = 0 (C02); the dispersion relation $\omega=\sqrt{gk\tanh kH}$ *(7.28)* (C03); phase speed, deep and shallow water (C04); particle orbits (C05); wave energy $E=\tfrac12\rho ga^2$ *(7.42)* and its flux (C06)",
        r"§7.3 capillary–gravity waves and the slowest ripple, $c_{min}=(4g\sigma/\rho)^{1/4}$ *(7.58)* (C07)",
        r"§7.4 standing waves and seiches (C08)",
        r"§7.5 group velocity $c_g=d\omega/dk$ *(7.67)* and $F=Ec_g$ *(7.71)* (C09); rays and refraction (C10)",
        r"§7.6 the hydraulic jump (and solitons) (C11); Stokes waves and Stokes drift $\bar u_L=a^2\omega ke^{2kz_0}$ *(7.85)* (C12)",
        r"§7.7 waves on a density interface, $\omega=\varepsilon\sqrt{gk}$ *(7.95)* (C13); barotropic and baroclinic modes, reduced gravity $c=\sqrt{g'H}$ *(7.116)* (C14)",
        r"§7.8 internal waves, $\omega=N\cos\theta$ *(7.139)* (C15); $\mathbf c\perp\mathbf c_g$, beams and $\mathbf F=\mathbf c_gE$ *(7.159)* (C16)",
    ],
    prerequisites=[
        "velocity potential and Laplace's equation (Ch. 6 §6.2)",
        "unsteady Bernoulli (Ch. 4 §4.9)",
        "kinematic and dynamic conditions at a moving surface, surface tension (Ch. 4 §4.10, Ch. 1 §1.6)",
        "control-volume momentum (Ch. 4 §4.4)",
        "path lines (Ch. 3 §3.2)",
        "Boussinesq equations and the buoyancy frequency N (Ch. 4 §4.9, Ch. 1 §1.10)",
        "complex exponentials (Ch. 1, Ch. 6)",
    ],
)
for _c in nb.cells:                                     # the title cell's promise: equations are shown, not only cited
    _c.source = _c.source.replace("equations are cited by their numbers so you can follow along in your copy",
                                  "every equation is shown in full together with its number, so you can follow along in your copy")
nb.explainer_index([
    ("dispersion_relation", "Why do long waves outrun short ones?",
     "C03 C04: ω(k), c(λ), the deep and shallow limits and the pressure under a wave on one clock"),
    ("particle_orbits", "What does the water under a wave actually do?",
     "C05 C12: orbits shrink with depth, circles become flat ellipses, and at finite amplitude they fail to close — Stokes drift"),
    ("capillary_gravity_waves", "Why is there a slowest ripple?",
     "C07: gravity for long waves, surface tension for short ones, a minimum speed of about 23 cm/s in between"),
    ("seiche_standing_waves", "Which waves can live in a lake?",
     "C08: two opposite waves make fixed nodes; the walls pick the periods"),
    ("group_velocity_packets", "Where does a wave's energy go?",
     "C09: crests move at c, the envelope and the energy at c_g = dω/dk"),
    ("wave_rays_refraction", "Why do waves arrive parallel to the beach?",
     "C10: frequency is fixed along a ray, the wavelength shrinks with depth, the crest swings round"),
    ("hydraulic_jump", "How high does the water jump?",
     "C11: momentum sets the height, energy sets the direction"),
    ("two_layer_modes", "What is a baroclinic mode?",
     "C13 C14: a layer over deep water rings in two ways — surface and interface in phase, or in antiphase and slow"),
    ("internal_wave_beams", "Why does energy leave at right angles to the crests?",
     "C15 C16: ω = N cos θ, c ⟂ c_g, the St Andrew's cross"),
])
nb.setup()
nb.code(r"""
import numpy as np                                       # arrays and maths (Ch. 1 primer P03)
import sympy as sp                                       # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                          # static figures (Ch. 1 primer P01)
import plotly.graph_objects as go                        # extra plotly traces for page-surviving figures (Ch. 1 P41)
import warnings                                          # to catch (and show) one deliberate physics warning
from fluidpy import ch07_gravity_waves as ch07           # the tested chapter-7 module; it re-exports core.waves
from fluidpy.core import waves as W                      # the new wave toolkit: dispersion relations, packets, rays
from fluidpy import ch01_introduction as ch01            # Ch. 1: surface tension of water (IAPWS)
from fluidpy import ch03_kinematics as ch03              # Ch. 3: path lines (Example 3.1)
from fluidpy import ch04_conservation_laws as ch04       # Ch. 4: the wave test field and the bore speed
from fluidpy import ch05_vorticity_dynamics as ch05      # Ch. 5: vortex sheets
from fluidpy.core import interfaces, bernoulli, potential, streamfunction, similarity, stratification   # earlier tools
from fluidpy.core.interact import slider_figure, animate_figure, live   # plotly sliders (P17), time players, widgets (P47)
from fluidpy.core.anim import animate, ffmpeg_path       # matplotlib animations (P16); show_animation came with setup
from fluidpy.core.style import COLORS, savefig           # the house palette and a helper that saves PNGs to outputs/ch07
from tools.convergence import observed_order             # slope of log(error) vs log(step) (P13)
from scripts.ch07_drawings import tank, wave_surface, orbit_ghosts, cv_box, beam_labels, two_layer_sketch, NAVY  # drawing only
ffmpeg_path()                                            # find ffmpeg once, for the MP4 animations
G = ch07.G_BOOK                                          # g = 9.81 m/s², the book's value used by every ch07 function
plt.rcParams["contour.negative_linestyle"] = "solid"     # negative ψ levels are streamlines too: draw them solid
import logging                                           # Python's message system (matplotlib reports font fallbacks through it)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless "font weight not found" notes


def recolor(fig, colors, dashes=None):                   # give plotly traces the notebook's colours by trace name
    for tr in fig.data:                                  # every trace of every slider step
        if tr.name in colors:                            # a name we assigned a colour to
            tr.line.color = colors[tr.name]              # same colour meaning as in the matplotlib figures
            tr.marker.color = colors[tr.name]            # markers too (for "markers" traces)
        if dashes and tr.name in dashes:                 # optional dash pattern ("dash", "dot")
            tr.line.dash = dashes[tr.name]
    return fig                                           # the same figure, restyled


def logaxes(fig, x=None, y=None):                        # switch plotly axes to log scale with fixed ranges [lo, hi]
    if x is not None:                                    # x range in data units
        fig.update_xaxes(type="log", range=[np.log10(x[0]), np.log10(x[1])])   # plotly wants log10 of the limits
    if y is not None:                                    # y range in data units
        fig.update_yaxes(type="log", range=[np.log10(y[0]), np.log10(y[1])])
    return fig                                           # the same figure


print(len([n for n in dir(ch07) if not n.startswith("_")]), "public names in fluidpy.ch07_gravity_waves; g =", G, "m/s²")
""", explain=r"""
1. Numerical, symbolic and plotting libraries (all primed in Ch. 1).
2. `ch07` is the chapter module. It **re-exports the new wave toolkit `core.waves`**, so `ch07.omega_gravity` and
   `W.omega_gravity` are the same function; every function cites its § and equation and is tested in `tests/test_ch07.py`.
3. Earlier chapters' modules are used in the recaps and for parity checks (the Ch. 4 wave test field, the Ch. 5 vortex
   sheet, the Ch. 1 surface tension, the Ch. 3 path line).
4. `scripts/ch07_drawings.py` only draws (tank outlines, surfaces, orbit ghosts, the jump's control volume, beam labels,
   the two-layer sketch); it computes no physics.
5. `G = 9.81` m/s²: every ch07 function uses the book's g (Ch. 1's `G0 = 9.80665` appears only where we call Ch. 1–4
   functions, and we pass g explicitly there).
6. `recolor` and `logaxes` only restyle plotly figures.
""")
nb.md(r"""
### ⚠️ Conventions in this chapter (read once; each is repeated where it bites)

| Symbol | Meaning here | Before (and where it bites) |
|---|---|---|
| ω | angular **frequency** of a wave [rad/s] | the vorticity in Ch. 2–6 |
| η(x, t) | **elevation** of the free surface above z = 0 [m] | Ch. 4's $\frac{D\eta}{Dt}=0$ *(4.91)* used η for the surface *function*; here that function is f = z − η |
| z, H | z **up**, still surface z = 0, flat bottom z = −H; H = ∞ means deep water (`H=np.inf` in the code) | §7.7's two-layer problem: origin at the free surface, interface at z = −H |
| p, p′ | p is **gauge** pressure; the wave part $p'\equiv p+\rho gz$ *(7.30)* in §7.2 | in §7.8, $p=\bar p(z)+p',\ \rho=\bar\rho(z)+\rho'$ *(7.124)*: p′ is measured from the stratified background |
| ζ | a particle's vertical excursion (§7.2, §7.6, §7.8) | the **interface** displacement in §7.7 |
| θ | the local phase θ(x, t) in $\eta=a(x,t)\cos[\theta(x,t)]$ *(7.72)* | the angle of **K** above the horizontal in $\omega=N\cos\theta$ *(7.139)* |
| g′ | the book's reduced gravity $g'=g\big(\frac{\rho_2-\rho_1}{\rho_2}\big)$ *(7.117)* (lower density below the line) | Ch. 4's version divides by ρ₁: 0.2 % apart in the ocean, 25 % for oil over water |
| E, F | surface and interface waves: E per unit **area** [J/m²], F per metre of **crest** [W/m] | internal waves (§7.8): E per unit **volume** [J/m³], F per unit area [W/m²] |
| g | 9.81 m/s² in every ch07 function (`ch07.G_BOOK`) | Ch. 1's `G0` = 9.80665 m/s² elsewhere in `fluidpy.core` |

We compute with the book's conventions and name the other one wherever it differs.

**Colours** (one meaning each, also in the explainers): free surface **blue** · particles and orbits **teal** · phase
(crests, c) **orange** · group and energy (envelope, $c_g$, F) **purple** · pressure **amber** · surface tension **rose** ·
lower (denser) layer **navy** · deep-water limit grey dashed · shallow-water limit grey dotted · residuals black ·
a printed slip shown as a ghost **rose dashed**.
""")
nb.md(r"""
**Book slips we correct in this notebook** (each is explained where it is used):

| The book prints | We use |
|---|---|
| (7.66) with $\tfrac12\Delta\omega\,x$ in the envelope | $\eta=2a\cos\big(\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,t\big)\cos(kx-\omega t)$ — with x the envelope could not move (D19) |
| "u from (7.28)" before the energy flux | u from (7.27), $u=a\omega\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t)$; (7.28) is the dispersion relation (D14) |
| (7.105) with $e^{i(kz-\omega t)}$ | $\phi_2=Ce^{kz}e^{i(kx-\omega t)}$ — the printed form fails Laplace's equation (D29) |
| (7.98) with ∂φ₁/dz | ∂φ₁/∂z: $\frac{\partial\phi_1}{\partial z}=\frac{\partial\eta}{\partial t}$ at z = 0 |
| §7.1 "y = 0" for the still surface | z = 0 |
| the Ursell remark citing (7.88) | the terms of the KdV equation (7.87), $\frac{\partial\eta}{\partial t}+c_0\frac{\partial\eta}{\partial x}+\frac32c_0\frac\eta H\frac{\partial\eta}{\partial x}+\frac16c_0H^2\frac{\partial^3\eta}{\partial x^3}=0$ |
| $\omega=kN/K$ *(7.138)* and $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$ *(7.145)* while Fig. 7.29 draws k < 0 | $\omega=N\lvert k\rvert/K$ and $\mathbf c_g=\nabla_{\mathbf K}\omega$ for either sign (D36) |
""")
nb.md(r"""
> 🔁 **Tools from earlier chapters used in this one** (one line each where first used; the full primers are in
> `knowledge/primers.md`): partial derivative (P25), first-order Taylor (P26) and its multivariable form (P98), definite
> integral (P27), `solve_ivp` (P31, P94), RK4 by hand (P95), sympy (P40) and `expand`/`series` (P117), the linear
> second-order ODE and the $e^{\lambda z}$ trial (P44), Euler's formula (P45), the chain rule (P49, P91), orders of
> smallness (P68), level sets and the normal (P75), `np.meshgrid` (P76), broadcasting (P77), contour and quiver plots
> (P78), eigenvalues — the idea only (P80), the complex conjugate (P81), `quad`/`dblquad` (P87), parametric curves (P92),
> `np.expm1` and cancellation (P107), `brentq` (P108), momentum flux (P114), Schwarz's theorem (P121), power of a force
> (P127), reduced gravity (P131), moving level sets (P132), Fourier modes and the FFT (P142), the Dirac delta (P149),
> integrals of sines and cosines over a period (P151), the complex plane in numpy (P153), sympy for complex algebra
> (P158), Vieta (P71), `animate` (P16), `slider_figure` (P17), `show_viz` (P18), live widgets (P47), `assert np.allclose`
> (P15), log–log slopes (P13). New primers of this chapter: P165–P178.
""")

# =====================================================================================================================
# A.1 §7.1 — C01 (N01–N11, N167, D01, P165)
# =====================================================================================================================
nb.section("7.1", "Introduction", intro=r"""
**What is this section about?** Before any fluid mechanics, the words: a wave's amplitude, wavelength, period,
frequency, wavenumber and phase, the speed at which a crest moves, the wavenumber vector of a wave travelling in any
direction, and what an observer on a current measures. Everything later in the chapter is written in these words.
""")
core("C01", r"The sinusoidal travelling wave and its vocabulary", r"""
What exactly moves at "the wave speed" — and how do λ, T, k, ω and c fit together?
""", eqs=("7.1", "7.4"))
problem(r"""
Stand on a pier and watch swell arrive: a buoy tied below you rises and falls every eight seconds, but it does not travel
toward the beach. What travels is the *shape* — the crest. We want exact words for that shape and one clean definition of
how fast it moves, because every later result (how fast a tsunami crosses an ocean, how fast its energy arrives) is a
statement about those speeds.
""")
note("N01", "Three families of waves", r"""
Waves on an interface are restored by gravity and surface tension (this chapter's first half); internal waves inside a
stratified fluid by buoyancy (its end); compression waves by compressibility (Ch. 15). Water waves are neither purely
longitudinal nor transverse: the water goes round in loops (C05). The whole chapter assumes small amplitude (so
everything is linear) and waves fast compared with the Earth's rotation (rotation joins in Ch. 13: Poincaré, Kelvin and
Rossby waves).
""")
nb.md(r"""
#### The idea

```
snapshot at t:       η(x) = a cos(kx − ωt)      crest where the phase kx − ωt = 0, 2π, 4π, …
a moment Δt later:   the same cosine, slid right by c·Δt
phase fixed:  kx − ωt = const  ⇒  x = (ω/k)t + const  ⇒  c = ω/k
```

| symbol | name | unit | from the others |
|---|---|---|---|
| a | amplitude (crest height above the still level) | m | — |
| λ | wavelength (crest to crest) | m | λ = 2π/k |
| k | wavenumber (radians of phase per metre) | rad/m | k = 2π/λ |
| T | period (crest to crest at a fixed point) | s | T = λ/c |
| ν | cyclic frequency | Hz | ν = 1/T |
| ω | angular frequency (radians of phase per second) | rad/s | ω = 2πν |
| c | phase speed (speed of a crest) | m/s | c = ω/k = λν |

In these letters the wave is $\eta(x,t)=a\cos\big[\frac{2\pi}{\lambda}(x-ct)\big]$ *(7.1)*. **A wave is a moving phase, not
moving water.**
""")
confusion(r"""
in Ch. 3–6 ω was the vorticity (a spin rate, 1/s). From here on ω is the wave's angular frequency (rad/s): how many
radians of phase pass a fixed point per second, $\omega=2\pi/T$.
""")
P("P165", "phase of a wave", r"""
The phase is the argument of the cosine, kx − ωt, measured in radians: 0 at a crest, π at a trough, 2π at the next
crest. Radians per metre (k) and radians per second (ω) are the natural units; cycles (1/λ, ν) differ by 2π. A point that
keeps its phase fixed moves to the right when ω/k > 0.
""", code=r"""
import numpy as np                                # numbers
k, w = 2*np.pi/100, 2*np.pi/8                     # λ = 100 m, T = 8 s → k [rad/m], ω [rad/s]
x = np.array([0.0, 25.0, 50.0])                   # three positions [m]
print(k*x - w*0.0, np.cos(k*x))                   # phases [0, 1.571, 3.142] rad → cos [1, 0, −1]
""")
note("N03", "The same wave in k and ω", r"""
With $k=2\pi/\lambda$ and $\omega=2\pi/T$ this is $\eta(x,t)=a\cos\big[\frac{2\pi}{\lambda}(x-ct)\big]$ *(7.1)* with
$c=\omega/k$. Number: λ = 100 m, T = 8 s → k = 0.0628 rad/m, ω = 0.785 rad/s, c = 12.5 m/s.
""", equation=EQ["7.2"], ref="7.2")
D("D01", ref="7.4")
note("N04, N05", "Where the crest condition and the speed sit", r"""
The crest condition $\frac{2\pi}{\lambda}(x_{crest}-ct)=2n\pi=kx_{crest}-\omega t$ *(7.3)* is D01's step 2; the integer n
only labels which crest. The speed $c=\omega/k=\lambda\nu$ *(7.4)* is D01's result: a crest's Δx/Δt.
""")
nb.worked_example("a swell with λ = 100 m and T = 8 s", r"""
1. $k=2\pi/\lambda=6.2832/100=0.0628$ rad/m.
2. $\omega=2\pi/T=6.2832/8=0.785$ rad/s.
3. $\nu=1/T=0.125$ Hz.
4. $c=\omega/k=0.785/0.0628=12.5$ m/s — the same as $\lambda\nu=100\times0.125$.
5. In 2 s a crest moves $c\,\Delta t=25$ m: a quarter wavelength.
""")
remind([
    ("Python dictionaries (fluidpy results)", "fluidpy returns several named results at once as a dictionary: `p[\"c\"]` (Ch. 1 P23)."),
    ("f-strings", "`f\"{x:.3g} m/s\"` puts a number with 3 significant figures into a sentence (Ch. 1 P04)."),
    ("np.linspace and np.logspace", "evenly spaced numbers on a linear axis, or with a constant factor on a log axis (Ch. 1 P06)."),
    ("matplotlib figures", "`fig, ax = plt.subplots()`, `ax.plot`, labels with units; every figure is followed by What you see / How to read it / What would change if (Ch. 1 P01)."),
    ("assert np.allclose", "stops the notebook if two results differ beyond a tolerance — our proof that a from-scratch version matches the library (Ch. 1 P15)."),
])
nb.code(r"""
p = ch07.wave_parameters(lam=100.0, T=8.0)               # any consistent pair → all six descriptors, (7.1)–(7.4)
print({name: round(v, 5) for name, v in p.items()})      # k [rad/m], lam [m], omega [rad/s], T [s], nu [Hz], c [m/s]
t = np.array([0.0, 1.0, 2.0])                            # three times [s]
print("crest n = 0 at x =", ch07.crest_positions(t, p["k"], p["omega"], n=0), "m")   # (7.3) solved for x_crest
""", explain=r"""
1. `wave_parameters` accepts any consistent pair (here λ and T) and returns all six quantities; it refuses an
   inconsistent or under-determined input.
2. `crest_positions` solves the crest condition $kx_{crest}-\omega t=2n\pi$ *(7.3)* for $x_{crest}$ at each time: the crest
   advances 12.5 m every second — c = 12.5 m/s.
""")
nb.md(r"""
**From scratch — track the crest on a fine grid.** `np.argmax` returns the index of the largest value, so
`x[np.argmax(eta)]` is where the surface is highest; `np.polyfit(ts, xs, 1)[0]` is the slope of the best straight line
through the points (ts, xs) — the crest's speed.
""")
nb.check_agree(r"""
x = np.linspace(-10, 110, 120001)                        # a 1 mm grid [m]
ts = np.array([0.0, 0.5, 1.0, 1.5, 2.0])                 # five times [s]
xs = []                                                   # the crest positions we find [m]
for t_ in ts:                                             # at each time …
    eta = W.sinusoid(x, t_, 1.0, p["k"], p["omega"])     # (7.2): η = a cos(kx − ωt) with a = 1 m
    xs.append(x[np.argmax(eta * (x < 50))])              # … the highest point left of 50 m: the crest that began at 0
slope = np.polyfit(ts, xs, 1)[0]                          # speed of the tracked crest [m/s]
print(f"tracked crest speed = {slope:.3f} m/s, ω/k = {p['c']:.3f} m/s")
assert np.allclose(slope, p["c"], rtol=1e-3)             # same numbers → the crest really moves at c = ω/k
""")
nb.md("Tracking the highest point of the cosine on a 1 mm grid gives 12.500 m/s — the crest really moves at ω/k.")
nb.figure(r"""
lam, a, T = 100.0, 1.0, 8.0                               # the swell of the tiny example [m, m, s]
k, w = 2*np.pi/lam, 2*np.pi/T                             # k [rad/m], ω [rad/s]
x = np.linspace(0, 200, 801)                              # two wavelengths [m]
fig, ax = plt.subplots(figsize=(7, 3))                   # one panel
ax.plot(x, W.sinusoid(x, 0.0, a, k, w), color=COLORS["blue"], lw=2, label="t = 0 s")          # (7.2) at t = 0
ax.plot(x, W.sinusoid(x, 2.0, a, k, w), color=COLORS["blue"], lw=2, ls="--", label="t = 2 s")  # the same wave 2 s later
xc0, xc2 = ch07.crest_positions(0.0, k, w, 1), ch07.crest_positions(2.0, k, w, 1)   # crest n = 1 at both times [m]
ax.plot([xc0, xc2], [a, a], "o", color=COLORS["orange"], ms=8)                        # the tracked crest (orange)
ax.annotate("", xy=(xc2, 1.25), xytext=(xc0, 1.25), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2))
ax.text(0.5*(xc0 + xc2), 1.35, "cΔt = 25 m", color=COLORS["orange"], ha="center", fontsize=9)
xb = 60.0                                                  # a buoy moored at x = 60 m
for tt, mk in ((0.0, "s"), (2.0, "D")):                    # its height at the two times
    ax.plot(xb, W.sinusoid(xb, tt, a, k, w), mk, color=COLORS["teal"], ms=8, mfc="white" if tt else COLORS["teal"])
ax.text(xb + 3, -1.35, "buoy: moves up and down only", color=COLORS["teal"], fontsize=8)
ax.set_xlabel("x [m]"); ax.set_ylabel("η [m]"); ax.set_ylim(-1.6, 1.7)   # axis labels with units
ax.set_title("The shape moves; the water at the buoy only goes up and down"); ax.legend(loc="lower right", fontsize=8)
savefig(fig, "ch07", "c01_travelling_wave"); plt.show()   # save the PNG to outputs/ch07, then draw
""", see=r"""Two snapshots of the same cosine, the dashed one shifted right by 25 m; the orange crest moved along with it;
the buoy (teal square at t = 0, open diamond at t = 2 s) changed height, not position.""",
    read=r"""Measure the shift of any crest and divide by the time: $25\ \text{m}/2\ \text{s}=12.5$ m/s $=\omega/k$.""",
    change=r"""…T were 16 s with the same λ: ω halves and the crest moves only 12.5 m in 2 s — half as fast
($c=\lambda/T$).""")
note("N02", "Why sinusoids", r"""
For linear waves any surface shape is a sum of sinusoids,
$\eta=\sum_k\hat\eta_ke^{i(kx-\omega(k)t)}+\text{c.c.}$, and each one moves with its own ω(k) (our
`W.linear_evolve` does the general case with the FFT). So one sinusoid is the building block.
""")
remind([
    ("Fourier modes and the FFT", r"on a periodic domain `np.fft.fft` splits a signal into modes $e^{ikx}$, ∂/∂x becomes multiplication by ik, and `np.fft.ifft` adds them back (Ch. 5 P142)."),
])
nb.code(r"""
L = 2*np.pi/0.02                                         # a periodic domain 314.16 m long [m]
x = np.linspace(0, L, 4096, endpoint=False)              # grid [m]
k1, k2 = 5*0.02, 20*0.02                                 # two deep-water modes that fit the domain [rad/m]
eta0 = np.cos(k1*x) + np.cos(k2*x)                       # the initial surface: a sum of two sinusoids [m]
eta1 = W.linear_evolve(eta0, x, 1.0, direction=+1, H=np.inf)   # every Fourier mode moved with its own ω(k) for 1 s
E0, E1 = np.fft.fft(eta0), np.fft.fft(eta1)              # the two spectra
for j, kk in ((5, k1), (20, k2)):                        # mode index j has k = j·2π/L
    shift = (np.angle(E0[j]) - np.angle(E1[j]))/kk       # phase lost by the mode / k = distance moved [m]
    print(f"k = {kk:.2f} rad/m: moved {shift:.2f} m in 1 s; deep-water c = √(g/k) = {np.sqrt(G/kk):.2f} m/s")
""", explain=r"""
1. `linear_evolve` takes the FFT of η₀, multiplies each mode by $e^{-i\omega(k)t}$ (one-way, to +x) with the deep-water
   ω(k) $=\sqrt{gk}$ (derived in C03–C04; take it on trust for now), and adds the modes back with the inverse FFT.
2. The phase each mode lost, divided by its k, is how far it moved: the long mode (k = 0.1 rad/m) moves 9.90 m in 1 s,
   the short one (k = 0.4 rad/m) 4.95 m — the long one twice as fast. That is dispersion, explained in C04.
""")
note("N06–N09", "Waves in any direction", r"""
$\eta=a\cos(kx+ly+mz-\omega t)=a\cos(\mathbf K\cdot\mathbf x-\omega t)$ *(7.5)* with the wavenumber vector
$\mathbf K=(k,l,m)$ and $K^2=k^2+l^2+m^2$ *(7.6)*. Crests are the planes $\mathbf K\cdot\mathbf x-\omega t=2n\pi$,
perpendicular to K and $\lambda=2\pi/K$ *(7.7)* apart along K; they advance along K at
$\mathbf c=(\omega/K)\,\mathbf e_K,\ \mathbf e_K=\mathbf K/K$ *(7.8)* (met again in C16, where c and the energy velocity
are perpendicular).
""")
remind([
    ("dot product K·x", r"$\mathbf K\cdot\mathbf x=kx+ly+mz$, the sum of the products of components (Ch. 2 §2.2)."),
    ("level sets and the normal to a surface", "a surface where a function is constant has the function's gradient as its normal (Ch. 2 P75): here the gradient of K·x is K itself, so crests are perpendicular to K."),
])
nb.code(r"""
K = np.array([1.0, 1.0])                                  # a wavenumber vector (k, l) [rad/m]
print("|K| =", np.hypot(*K), "rad/m;  λ = 2π/K =", 2*np.pi/np.hypot(*K), "m")   # (7.6), (7.7)
print("c =", ch07.phase_velocity_vector(K, 1.0), "m/s")  # (7.8) with ω = 1 rad/s: along K, size ω/K
print("trace speeds (ω/k, ω/l) =", ch07.trace_velocities(K, 1.0), "m/s")   # how fast crests cut each axis
""")
remind([
    ("np.meshgrid and the grid layout", "`X, Y = np.meshgrid(x, y)` gives every (x, y) pair; rows run along y, columns along x (Ch. 2 P76)."),
    ("contour and quiver plots", "`ax.contour(X, Y, F, levels=[…])` draws the lines where F takes the given values; `quiver` draws arrows (Ch. 2 P78)."),
    ("numpy broadcasting", "an operation between arrays of compatible shapes acts at every point at once (Ch. 2 P77)."),
])
P("P179", "np.stack and np.c_", r"""
Two ways to glue arrays together. `np.stack([X, Y], axis=-1)` puts equal-shaped arrays side by side on a **new last
axis**, so every grid point carries its (x, y) pair — the layout `plane_wave` expects for positions. `np.c_[a, b]`
stacks 1-D arrays as **columns**: two points (x₀, y₀) and (x₁, y₁) given as the columns a and b become one row of x's
and one row of y's, which is what `ax.plot(*…)` needs to draw the segment between them.""", code=r"""
import numpy as np                                        # arrays
X, Y = np.meshgrid([0, 1], [0, 1])                        # a 2 × 2 grid of x and y values
print(np.stack([X, Y], axis=-1).shape)                    # (2, 2, 2): the last axis holds (x, y) at each point
print(np.c_[[0, 0], [3, 4]])                              # columns (0,0) and (3,4) → rows [0 3] (the x's) and [0 4] (the y's)
""")
nb.figure(r"""
xg = np.linspace(0, 20, 201)                              # a 20 m × 20 m patch [m]
XX, YY = np.meshgrid(xg, xg)                              # every (x, y) pair
K = np.array([1.0, 1.0])                                  # the wavenumber vector [rad/m]
eta = ch07.plane_wave(np.stack([XX, YY], axis=-1), K, 1.0, a=1.0, t=0.0)   # (7.5) at t = 0 (components on the last axis)
fig, ax = plt.subplots(figsize=(5, 5))                    # square panel so angles look right
ax.contour(XX, YY, eta, levels=[0.999], colors=[COLORS["orange"]], linewidths=1.5)   # crest lines η = a
ax.annotate("", xy=(17, 5), xytext=(14, 2), arrowprops=dict(arrowstyle="->", color=COLORS["ink"], lw=2))   # K direction
ax.text(15.2, 2.2, "K", fontsize=12)
lamK = 2*np.pi/np.hypot(*K)                               # spacing of crests along K [m]
ax.plot([2*np.pi, 4*np.pi], [0.3, 0.3], color=COLORS["muted"], lw=3)   # spacing along x: 2π/k = 6.28 m
ax.text(3*np.pi, 0.8, "2π/k = 6.28 m along x", ha="center", fontsize=8, color=COLORS["muted"])
c0 = np.array([np.pi, np.pi])                             # a point on a crest line (K·x = 2π)
ax.plot(*np.c_[c0, c0 + lamK*K/np.hypot(*K)], color=COLORS["accent"], lw=3)   # spacing along K: λ = 4.44 m
ax.text(c0[0] + 3.6, c0[1] + 1.2, f"λ = 2π/K = {lamK:.2f} m", color=COLORS["accent"], fontsize=8,
        bbox=dict(facecolor="white", edgecolor="none", pad=1))
ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
ax.set_title("Crests perpendicular to K: further apart along the axes", fontsize=10)
savefig(fig, "ch07", "c01_plane_wave"); plt.show()
""", see=r"""Straight crest lines (orange) at 45°, $\lambda=2\pi/K=4.44$ m apart along K (purple bar) but 6.28 m apart along
x (grey bar) — our remake of Fig. 7.1 `N167`.""",
    read=r"""Along an axis the crests are further apart and pass a fixed point at the same rate ω, so they *seem* to move
faster: $c_x=\omega/k=1.0$ m/s against $c=\omega/K=0.707$ m/s along K.""",
    change=r"""…K turned to (1, 0): the crests stand perpendicular to x, $c_x=c$, and $c_y=\infty$ (the crest lines never cross
the y-axis).""")
confusion(r"""
`N10` the trace speeds $c_x=\omega/k$, $c_y=\omega/l$, $c_z=\omega/m$ are each ≥ $c=\omega/K$, so they are **not** the
components of c (components would be smaller). Reciprocals do not add like components: $c_x\mathbf e_x+c_y\mathbf e_y\ne
\mathbf c$ — its length is 1.41 m/s here against 0.71 m/s for c.
""")
note("N11", "A current shifts the frequency you see", r"""
In water moving at uniform U the crests are carried along, so a fixed probe sees the observed frequency ω₀ below; ω is the
intrinsic one, measured drifting with the water — one substitution $x'=x-Ut$ in the phase (the frames of Ch. 3 P96). A
frozen pattern (ω = 0) swept past at U gives $\omega_0=Uk$. From here on every ω is intrinsic.
""", equation=EQ["7.9"], ref="7.9")
remind([
    ("change of frame x′ = x − Ut", "two frames moving apart at constant U see velocities that differ by U (Galilean frames, Ch. 3 P96); writing the phase with x′ = x − Ut gives (7.9)."),
])
nb.code(r"""
w = 2*np.pi/8; k = w**2/G                               # an 8 s deep-water swell: ω [rad/s], k = ω²/g [rad/m] (C04)
for U in (-1.0, 0.0, 1.0):                               # a current against, none, and with the swell [m/s]
    w0 = ch07.doppler_frequency(w, np.array([U, 0.0]), np.array([k, 0.0]))   # (7.9): ω₀ = ω + U·K
    print(f"U = {U:+.0f} m/s: ω₀ = {w0:.4f} rad/s, observed period {2*np.pi/w0:.2f} s")
""", explain=r"""
The first line uses the deep-water rule $k=\omega^2/g$ (derived in C03–C04; take it on trust for now). Against a 1 m/s
current the probe sees the 8 s swell as an 8.70 s wave, with the current as a 7.41 s wave. Ch. 13 uses
exactly this for mountain waves and Rossby waves in a wind.
""")
whatif(r"""
…the medium itself set a rule between ω and k? For water it does — $\omega=\sqrt{gk\tanh kH}$ *(7.28)* — and that single
rule makes long waves faster than short ones (C03–C04). The vocabulary stays; only the relation changes from wave to wave
(capillary C07, interfacial C13, internal C15).
""")

# =====================================================================================================================
# A.2 §7.2 — R01–R06, C02 (N12–N17, N168, D02–D04, P166), C03 (N18–N25, D05, D06, P167, P168)
# =====================================================================================================================
nb.section("7.2", "Linear Liquid-Surface Gravity Waves", intro=r"""
**What is this section about?** The heart of the chapter. We set up the flow under small waves on a layer of depth H —
Laplace's equation with conditions at the bottom and at the moving surface — make it linear, solve it, and read off the
dispersion relation, the speed of crests in deep and shallow water, the pressure under a wave, the orbits of the water
and the energy a wave carries.
""")
nb.recap("R01", "The velocity potential", r"""
Where the flow has no vorticity the velocity is the gradient of one scalar, $u=\partial\phi/\partial x,\ w=\partial\phi/
\partial z$ *(7.10)*. Waves started from rest by gravity stay irrotational (Kelvin's theorem, Ch. 5), so this is exact here.
""", where="Ch. 6 §6.2")
nb.recap("R02", "Laplace's equation", r"""
Incompressibility $\partial u/\partial x+\partial w/\partial z=0$ with (7.10) gives
$\partial^2\phi/\partial x^2+\partial^2\phi/\partial z^2=0$ *(7.11)* — the same Laplacian as Ch. 6 with y renamed z. A
function with zero Laplacian is *harmonic*; below, Ch. 6's five-point stencil (`potential.laplacian_residual`, y ↦ z)
tests a deep-water-shaped potential.
""", where="Ch. 6 §6.2")
remind([
    ("functions as arguments and lambda", "`lambda x, z: …` is a one-line function; fluidpy tools take fields as functions (Ch. 1 P29)."),
])
nb.code(r"""
phi = lambda x, z: np.exp(z) * np.sin(x)                # a deep-water-shaped potential with k = 1 rad/m [m²/s]
print(potential.laplacian_residual(phi, 0.3, -0.5))      # ∇²φ at (x, z) = (0.3, −0.5) m: stencil round-off, ≈ 0
""", explain=r"""
$\phi=e^z\sin x$ is harmonic ($\partial_x^2\phi=-\phi$ and $\partial_z^2\phi=+\phi$ cancel), so Ch. 6's five-point stencil
returns only round-off: this is the shape of potential C03 will find under deep-water waves.
""")
nb.recap("R03", "No flow through the bottom", r"""
A solid flat bottom at z = −H lets no water through: $w=\partial\phi/\partial z=0\ \text{on}\ z=-H$ *(7.12)* — the
wall-is-a-streamline condition of Ch. 6, $\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{surface}}$ *(6.16)*,
for a wall at rest.
""", where="Ch. 6 §6.2")
nb.recap("R04", "The kinematic condition at a moving surface", r"""
Water on the surface stays on it: the fluid velocity normal to the surface equals the surface's own normal speed,
$(\mathbf n\cdot\mathbf u)_{z=\eta}=\mathbf n\cdot\mathbf U_s$ *(7.13)*. Ch. 4 wrote it as
$\frac{D\eta}{Dt}=\frac{\partial\eta}{\partial t}+(\mathbf u\cdot\nabla)\eta=0$ *(4.91)*, on the surface. ⚠️ In Ch. 4 η was
that surface *function*; here η(x, t) is the *elevation* and the function is f = z − η. Below, Ch. 4's test wave (this
chapter derives it in C03) is put into the exact condition on its own surface.
""", where="Ch. 4 §4.10")
nb.code(r"""
f_imp, u_fn = ch04.linear_wave_fields(a=0.05, k=1.0, H=2.0)   # Ch. 4 test field: level-set function z − η and velocity
X = np.array([[0.3], [0.05*np.cos(0.3)]])                     # one point on the surface at t = 0: (x, η(x)) [m]
print(interfaces.kinematic_bc_residual(f_imp, u_fn, X, 0.0))   # D(z − η)/Dt there [m/s]: about 5e-3, not 0
print(ch04.wave_kinematic_residual(0.3, 0.0, a=0.05, k=1.0, H=2.0))   # the same number through Ch. 4's wrapper
w_test = np.sqrt(G*1.0*np.tanh(1.0*2.0))                      # ω of the test wave, √(gk tanh kH) with k = 1, H = 2 [rad/s] (C03)
print(f"aω = {0.05*w_test:.3f} m/s;  ka·aω = {1.0*0.05*0.05*w_test:.1e} m/s")   # the velocity scale and the size of the dropped terms
""", explain=r"""
The *exact* condition is not met by the linear solution: it misses by about $5\times10^{-3}$ m/s. The last line prints
the scales: with k = 1 rad/m, a = 5 cm and H = 2 m, $\omega=\sqrt{gk\tanh kH}=3.08$ rad/s (the formula C03 derives —
take it on trust for now), so the velocity scale is $a\omega\approx0.15$ m/s and $ka\cdot a\omega\approx8\times10^{-3}$
m/s. The residual $4.7\times10^{-3}$ m/s is of that order — the gap C02 is about.
""")
nb.recap("R05", "The free surface feels only the air", r"""
With no surface tension and no shear, stress continuity (Ch. 4 §4.10) says the water pressure just below the surface
equals the air's: $(p)_{z=\eta}=0$ *(7.19)* in **gauge** pressure (atmospheric = 0).
""", where="Ch. 4 §4.10")
nb.recap("R06", "Unsteady Bernoulli", r"""
For irrotational flow $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+\frac p\rho=\text{const}$
*(4.83)*; dropping the small square gives $\frac{\partial\phi}{\partial t}+\frac p\rho+gz\cong0$ *(7.20)* (the constant is
zero on the still surface far away). Below, Ch. 4's function and ours agree for $\partial\phi/\partial t=-0.5$ m²/s² at
z = −1 m.
""", where="Ch. 4 §4.9")
nb.code(r"""
print(bernoulli.unsteady_bernoulli_pressure(dphi_dt=-0.5, speed=0.0, z=-1.0, rho=1000.0, g=9.81))   # Ch. 4: p [Pa]
print(ch07.linear_bernoulli_pressure(-0.5, -1.0))        # (7.20) and p′ = p + ρgz (7.30) [Pa]
""", explain=r"""
Both give the gauge pressure $p=-\rho(\phi_t+gz)=-1000(-0.5-9.81)=10\,310$ Pa; `ch07` also returns the wave part
$p'=p+\rho gz=500$ Pa.
""")

# ---- C02 ------------------------------------------------------------------------------------------------------------
core("C02", r"The linear free-surface problem: move the conditions to z = 0", r"""
The surface is itself an unknown that moves. How can we impose conditions on it — and what do we give up?
""", eqs=("7.18",))
problem(r"""
A leaf floating on swell rides up and down with the surface and never leaves it; the air above pushes on the water with
the same pressure everywhere. Those are our two surface conditions. The trouble: they hold on the surface z = η(x, t),
which is part of the answer. For small, gentle waves we can apply them on the flat level z = 0 instead — and we want to
know exactly what that costs.
""")
note("N12", "The set-up", r"""
A liquid of constant density and uniform depth H, small slope a/λ ≪ 1 and small amplitude a/H ≪ 1, no surface tension
(C07 adds it), the air ignored, irrotational motion caused only by the waves; x along the propagation, z up. (The book
says "y = 0" once for the still surface; it means z = 0.)
""")
nb.figure(r"""
H, L, a, k = 1.0, 4.0, 0.15, 2*np.pi/2.0                  # a drawing: depth 1 m, 4 m of tank, 15 cm wave, λ = 2 m
fig, ax = plt.subplots(figsize=(7, 3))                    # one panel
tank(ax, H, L)                                            # bottom z = −H (hatched), dashed still level z = 0, z-axis
x = np.linspace(0, L, 400)                                # positions [m]
wave_surface(ax, x, a*np.cos(k*x), fill_to=-H)            # the free surface η = a cos kx (blue) with the water shaded
ax.annotate("", xy=(2.0, a), xytext=(2.0, 0), arrowprops=dict(arrowstyle="<->", color=COLORS["ink"]))   # amplitude a
ax.text(2.05, 0.04, "a", fontsize=11)
ax.annotate("", xy=(3.3, -H), xytext=(3.3, 0), arrowprops=dict(arrowstyle="<->", color=COLORS["ink"]))  # depth H
ax.text(3.35, -0.5, "H", fontsize=11)
ax.text(0.55, 0.22, "η(x, t)", color=COLORS["blue"], fontsize=10)
ax.annotate("g", xy=(3.8, -0.35), xytext=(3.8, 0.0), arrowprops=dict(arrowstyle="->", color=COLORS["ink"]), fontsize=10)
ax.set_ylim(-1.15, 0.35); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")
ax.set_title("Where the conditions live: on the blue curve, moved to the dashed line", fontsize=10)
savefig(fig, "ch07", "c02_geometry"); plt.show()
""", see=r"""The geometry of every problem in §7.2–§7.6 (our remake of Fig. 7.2 `N168`): the still level z = 0 (dashed),
the surface η = a cos kx (blue), the bottom z = −H, and gravity pointing down.""",
    read=r"""z = 0 is where the surface would be at rest; the conditions of this block are first written on the blue curve,
then moved to the dashed line.""",
    change=r"""…H → ∞: the bottom disappears and $w=\partial\phi/\partial z=0\ \text{on}\ z=-H$ *(7.12)* becomes "φ → 0 as
z → −∞" (deep water, D08).""")
nb.md(r"""
#### The idea

Any smooth quantity evaluated on the surface can be Taylor-expanded about the flat level:

```
condition lives on z = η(x, t)            (unknown, moving)
F(η) = F(0) + η F′(0) + …                 (Taylor about z = 0)
keep F(0):  error ≈ η F′ ≈ ka × (kept term)  ⇒  negligible when ka ≪ 1
```

**Linear = drop products of small quantities, then apply the conditions on z = 0.**
""")
P("P166", "Taylor transfer of a boundary condition", r"""
A condition on a moving boundary z = η can be written on the fixed level z = 0 by expanding each quantity F about
z = 0: F(η) ≈ F(0) + ηF′(0) (Ch. 1 P26). If F varies on the scale 1/k and η ~ a, the correction is about ka times F — so
for gentle waves we keep only F(0), and both of the book's approximations, $\big(\frac{\partial\phi}{\partial
z}\big)_{z=\eta}\cong\frac{\partial\eta}{\partial t}$ *(7.17)* and $\big(\frac{\partial\phi}{\partial z}\big)_{z=0}\cong
\frac{\partial\eta}{\partial t}$ *(7.18)*, make errors of the same order ka.
""", code=r"""
import numpy as np                                # numbers
k, eta = 1.0, 0.05                                # decay rate [1/m] and a small surface height [m]
F = lambda z: np.exp(k*z)                         # a quantity that varies on the scale 1/k
print(F(eta), F(0) + eta*k*F(0), F(0))            # 1.05127, 1.05, 1.0: dropping ηF′ costs 5 % = ka
""")
remind([
    ("moving level set and its normal speed", "a surface f = 0 has normal ∇f/|∇f| and a point riding on it keeps f = 0, so Df/Dt = 0 there (Ch. 4 P132)."),
    ("chain rule along a path", r"$df/dt=f_x\dot x+f_z\dot z+f_t$ along a moving point (Ch. 3 P91)."),
    ("level sets and the normal to a surface", "the gradient of a function is perpendicular to its level sets (Ch. 2 P75)."),
], lead="the tools of D02")
D("D02", ref="7.16")
note("N13–N15", "The pieces of D02", r"""
The upward unit normal of the surface f = z − η(x, t) = 0 is
$\mathbf n=\big(-(\partial\eta/\partial x)\mathbf e_x+\mathbf e_z\big)\big/\sqrt{(\partial\eta/\partial x)^2+1}$ *(7.14)*
(D02 step 3; for a slope of 0.1 it is (−0.0995, 0.9950), printed below). ⚠️ In Ch. 4 the level-set function was called
η; here η is the height and the function is z − η. Only the normal part of the surface velocity matters, so we may follow
the surface point straight above x: $\mathbf U_s=(\partial\eta/\partial t)\,\mathbf e_z$ *(7.15)*. The exact kinematic
condition (7.16) is D02's result — it is $D(z-\eta)/Dt=0$.
""")
nb.code(r"""
print(ch07.surface_normal(0.1))                          # (7.14) for a surface slope ∂η/∂x = 0.1: (n_x, n_z) [–]
""")
remind([
    ("limits and orders of smallness", "a product of two small quantities is smaller than either by one more factor of smallness (Ch. 2 P68)."),
    ("order-of-magnitude scaling", "replace each quantity by its typical size — ∂η/∂t ~ aω, ∂η/∂x ~ ka — and compare terms (Ch. 4 P130)."),
    ("first-order Taylor expansion", r"$f(z+dz)\approx f(z)+f'(z)\,dz$, error ∝ dz² (Ch. 1 P26)."),
    ("central differences", r"$f'(x)\approx[f(x+h)-f(x-h)]/2h$, error ∝ h² (Ch. 1 P21) — used in the from-scratch cell below."),
    ("power laws, log–log slopes, observed order", r"$y=Cx^p$ is a straight line of slope p on log–log axes; `observed_order` fits p (Ch. 1 P13) — used in the residual figure below."),
], lead="the tools of D03 and D04")
D("D03", ref="7.18")
note("N16", "The kinematic approximation", r"""
Dropping the slope × velocity product gives $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}\cong\frac{\partial\eta}
{\partial t}$ *(7.17)* (D03 step 3): the dropped term is about $ka\cdot a\omega$, the kept ones $a\omega$ — a relative
error of order ka.
""")
D("D04", ref="7.21")
note("N17", "The linearised dynamic condition", r"""
$\big(\frac{\partial\phi}{\partial t}\big)_{z=\eta}\cong\big(\frac{\partial\phi}{\partial t}\big)_{z=0}\cong-g\eta$
*(7.21)* is D04's result.
""")
nb.md(r"""
#### The linear problem in one box (the whole of C03 starts from it)

$$\nabla^2\phi=0\ \text{(7.11)},\qquad \frac{\partial\phi}{\partial z}=0\ \text{on } z=-H\ \text{(7.12)},\qquad
\Big(\frac{\partial\phi}{\partial z}\Big)_{z=0}=\frac{\partial\eta}{\partial t}\ \text{(7.18)},\qquad
\Big(\frac{\partial\phi}{\partial t}\Big)_{z=0}=-g\eta\ \text{(7.21)}.$$

Four linear equations on fixed boundaries: sums of solutions are solutions (note `N02`).
""")
nb.worked_example("how much do we throw away?", r"""
A swell with a = 0.5 m, λ = 50 m.

1. $k=2\pi/50=0.126$ rad/m.
2. $ka=0.063$.
3. The kept terms of the kinematic condition are about $a\omega$; the dropped slope × velocity term is about
   $ka\cdot a\omega$ — 6 % of the kept one.
4. A steep wave with ka = 0.3 (close to breaking, C12) loses 30 %: linear theory is then only a first guess.
""")
nb.code(r"""
for a in (0.05, 0.3):                                    # a gentle (ka = 0.05) and a steep (ka = 0.3) wave, k = 1 rad/m, H = 2 m
    r = ch07.free_surface_residuals(0.3, 0.0, a, 1.0, H=2.0)   # residual of every condition at x = 0.3 m, t = 0
    print(f"a = {a}:", {name: f"{v:.2e}" for name, v in r.items() if name != "scales"})   # m/s or m²/s²
""", explain=r"""
1. For the linear wave of C03 (k = 1 rad/m, H = 2 m) the function evaluates the residual of each condition of the box:
   `laplace` (7.11), `bottom` (7.12), `kinematic_linear` (7.18), `dynamic_linear` (7.21) — all at round-off: **the
   linear solution solves the linear problem exactly**.
2. `kinematic_exact` is the exact condition (7.16) on the true surface z = η, `kinematic_17` the partial step (7.17),
   `dynamic_exact` the full Bernoulli equation (4.83) with p = 0 on z = η. They are of order $10^{-3}$–$10^{-2}$ at
   a = 5 cm and much larger at a = 30 cm — the kinematic one 43 times: six times the amplitude, squared, plus the
   growing higher-order terms (at a single point the dynamic one also depends on where in the wave x sits).
""")
nb.md(r"""
**From scratch — the exact kinematic residual (7.16) at one point**, by central differences in x and t
(Ch. 1 P21: $f'(x)\approx[f(x+h)-f(x-h)]/2h$):
""")
nb.check_agree(r"""
h, x, t = 1e-6, 0.3, 0.0                                  # difference step, point [m] and time [s]
F = lambda X, Z, T: ch07.wave_fields(X, Z, T, 0.05, 1.0, H=2.0)   # the linear wave: η, u, w, … at (x, z, t)
eta = F(x, 0.0, t)["eta"]                                 # surface height at x [m]
eta_t = (F(x, 0.0, t + h)["eta"] - F(x, 0.0, t - h)["eta"])/(2*h)   # ∂η/∂t [m/s]
eta_x = (F(x + h, 0.0, t)["eta"] - F(x - h, 0.0, t)["eta"])/(2*h)   # ∂η/∂x [–]
res = F(x, eta, t)["w"] - eta_t - eta_x*F(x, eta, t)["u"]           # (7.16): w − η_t − η_x u, all on z = η [m/s]
print(f"our residual {res:.6e} m/s")
assert np.allclose(res, ch07.free_surface_residuals(x, t, 0.05, 1.0, H=2.0)["kinematic_exact"], rtol=1e-5)
""")
nb.md("Our three-line evaluation of (7.16) on the true surface matches the library's residual.")
nb.figure(r"""
ka = np.logspace(-3, -0.5, 12 if not FAST else 6)         # steepness values from 0.001 to 0.32
s = ch07.free_surface_residual_scan(ka, kH=1.0)           # max residuals over a wavelength, k = 1 rad/m fixed
fig, ax = plt.subplots(figsize=(6.5, 3.4))                # one log–log panel
ax.loglog(ka, s["kinematic_abs"], "k-o", ms=3, label="kinematic (7.16), absolute [m/s]")   # misfit of (7.16) on z = η
ax.loglog(ka, s["dynamic_abs"], "k--s", ms=3, label="dynamic (4.83) with p = 0, absolute [m²/s²]")   # misfit of full Bernoulli with p = 0 on z = η
ax.loglog(ka, s["kinematic_rel"], "-^", color=COLORS["muted"], ms=3, label="kinematic ÷ aω (relative)")   # the same misfit ÷ the kept term aω
ax.loglog(ka, 0.6*ka**2, ":", color=COLORS["ink"], lw=1); ax.text(3e-3, 1.2e-5, "slope 2", fontsize=8)   # guides
ax.loglog(ka, 0.8*ka, ":", color=COLORS["muted"], lw=1); ax.text(2e-3, 3e-3, "slope 1", fontsize=8, color=COLORS["muted"])   # slope-1 guide
for kv, lab in ((0.063, "swell ka = 0.063"), (0.3, "steep ka = 0.3")):   # two markers
    ax.axvline(kv, color=COLORS["orange"], ls="--", lw=1); ax.text(kv*1.05, 2e-6, lab, rotation=90, fontsize=7, color=COLORS["orange"])   # marker line and its label
ax.set_xlabel("steepness ka [–]"); ax.set_ylabel("residual of the exact condition")   # axis labels
ax.set_title("What linearisation drops grows like (ka)² — (ka)·aω", fontsize=10); ax.legend(fontsize=7, loc="upper left")   # the message as title
savefig(fig, "ch07", "c02_residual_scan"); plt.show()   # save the PNG to outputs/ch07, then draw
print(f"slopes: kinematic absolute {s['slope_abs']:.2f}, dynamic absolute {s['slope_abs_dynamic']:.2f}, "
      f"kinematic relative {s['slope_rel']:.2f}")   # observed log–log slopes (P13)
""", see=r"""Two black lines of slope 2 (the absolute misfit of the exact conditions, at fixed k) and one grey line of slope 1
(the misfit divided by the kept term aω); the fitted slopes are printed under the figure (2.02, 2.01 and 1.02).""",
    read=r"""The dropped terms are of order $(ka)\cdot a\omega$: **relative** to the kept terms the error is O(ka) — about 6 %
at the swell marker and about 30 % at the steep-wave marker. The approximation is as good as ka is small.""",
    change=r"""…the waves were in shallow water (kH = 0.2): the same slopes, but a/H ≪ 1 becomes the stricter condition (the
bottom term of D03 step 8).""")
confusion(r"""
"linear" does not mean "the flow is slow". It means every product of two wave quantities (slope × velocity, velocity²)
is dropped **and** the conditions are moved to z = 0; both errors are of the same order ka, so doing one without the
other is inconsistent.
""")
whatif(r"""
…the surface had surface tension? The dynamic condition gains a curvature term,
$\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=\frac\sigma\rho\frac{\partial^2\eta}{\partial x^2}-g\eta$ *(7.55)* (C07).
…there were a second fluid above? Both conditions become two-sided, (7.93)–(7.94) (C13) — and with a shear flow they give
the Kelvin–Helmholtz instability of Ch. 11.
""")

# ---- C03 ------------------------------------------------------------------------------------------------------------
core("C03", r"The dispersion relation", r"""
Given a wavelength and a depth, how fast must the surface oscillate — and why is there only one answer?
""", eqs=("7.28",))
problem(r"""
In a wave tank, push the paddle slowly and long waves come out; push it fast and short ones do. The water decides the
wavelength that goes with each frequency. We want that rule — the dispersion relation — because it is the chapter's hub:
speeds, orbits, energy flux, seiches, refraction and group velocity are all read off it. Ch. 4 already *used*
$\omega^2=gk\tanh kH$ as a test field (`ch04.linear_wave_surface`); here we derive it for the first time.
""")
idea(r"""
try   φ = f(z) · sin(kx − ωt)        (the surface is a cosine, so φ must be a sine)
Laplace (7.11)          ⇒  f″ = k² f      ⇒  f = A e^{kz} + B e^{−kz}
bottom (7.12)           ⇒  B/A fixed      ⇒  f ∝ cosh k(z + H)
kinematic (7.18)        ⇒  size fixed by a and ω
dynamic (7.21)          ⇒  only one ω for each k   ←  the dispersion relation
""")
P("P167", "separation of variables for a PDE", r"""
When a linear PDE and its boundaries are straight in x and z, try a product: one function of z times a wave in x.
Substituting turns ∂²/∂x² into −k², so the PDE becomes an ordinary differential equation for the z-part (Ch. 1 P44
solved those with an $e^{\lambda z}$ trial). The x-part is chosen by the boundary conditions — here the surface cosine.
(Not the same move as Ch. 1's "separation of variables" for an ODE, P42, which put each variable on its own side.)
""", code=r"""
import sympy as sp                                # symbolic algebra
x, z, k = sp.symbols('x z k', positive=True)     # coordinates and the wavenumber
f = sp.Function('f')                              # the unknown depth profile
phi = f(z)*sp.sin(k*x)                            # the product trial
print(sp.simplify((sp.diff(phi, x, 2) + sp.diff(phi, z, 2))/sp.sin(k*x)))   # -k**2*f(z) + f''(z): an ODE in z
""")
P("P168", "hyperbolic functions cosh, sinh, tanh", r"""
$\cosh x=(e^x+e^{-x})/2$, $\sinh x=(e^x-e^{-x})/2$, $\tanh=\sinh/\cosh$, $\coth=1/\tanh$, $\mathrm{sech}=1/\cosh$. Rules
used in this chapter: $\cosh^2-\sinh^2=1$, $\cosh^2+\sinh^2=\cosh2x$, $2\sinh x\cosh x=\sinh2x$, $(\cosh)'=\sinh$,
$(\sinh)'=\cosh$, $(\tanh)'=\mathrm{sech}^2$. Small x: cosh ≈ 1, sinh ≈ tanh ≈ x. Large x: cosh ≈ sinh ≈ $e^x/2$,
tanh → 1 (tanh 2 = 0.964). Here x = kH: tanh kH switches between shallow and deep water.
""", code=r"""
import numpy as np                                # numbers
x = np.array([0.1, 1.0, 2.0, 5.0])               # four arguments
print(np.cosh(x), np.sinh(x), np.tanh(x))         # tanh: 0.0997 0.762 0.964 1.000
print(np.cosh(x)**2 - np.sinh(x)**2)              # [1. 1. 1. 1.] — the identity
""")
note("N18", "The initial shape", r"""
The linear problem also needs an initial surface shape; the book picks η(x, 0) = a cos kx (the general case is note N02's
Fourier sum).
""")
note("N19", "The trial", "D05's step 1:", equation=EQ["7.22"], ref="7.22")
remind([
    ("linear second-order ODE and the e^{λz} trial", r"try $f=e^{\lambda z}$: $f''-k^2f=0$ gives $\lambda^2=k^2$, so λ = ±k — real roots, exponentials (Ch. 1 P44)."),
    ("exponent rules", r"$e^{a}e^{b}=e^{a+b}$; factoring $e^{-kH}$ out of $e^{kz}+e^{-2kH}e^{-kz}$ uses it (Ch. 1 P43)."),
], lead="two linear equations in two unknowns are solved in D05 steps 7–8; matching coefficients of cos(kx − ωt) is explained in D06 step 4")
D("D05", ref="7.26")
note("N20–N23", "Where D05 passes", r"""
D05 passes through $\frac{d^2f}{dz^2}-k^2f=0,\ f=Ae^{kz}+Be^{-kz}$ *(7.23)* (steps 2–4), the bottom condition
$k(Ae^{-kH}-Be^{+kH})\sin(kx-\omega t)=0$, i.e. $B=Ae^{-2kH}$ *(7.24)* (step 5), the kinematic condition
$k(A-B)=\omega a$ *(7.25)* (step 6), and ends at $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$
*(7.26)*.
""")
note("N24", "The velocities", r"""
One derivative each gives the velocities below (∂/∂x of the sine gives k cos; d/dz of cosh k(z + H) gives
k sinh k(z + H)). They came from the kinematic conditions alone: the pressure follows afterwards from Bernoulli. The code
cell checks that Ch. 4's test field is exactly this derived solution.
""", equation=EQ["7.27"], ref="7.27")
nb.code(r"""
x, z, t = 0.7, -0.4, 1.3                                  # a point under the wave [m] and a time [s]
f7 = ch07.wave_fields(x, z, t, 0.05, 1.0, H=2.0)          # this chapter's (7.26)–(7.27): a = 5 cm, k = 1 rad/m, H = 2 m
f4 = ch04.linear_wave_surface(x, z, t, a=0.05, k=1.0, H=2.0)   # Ch. 4's test field at the same point
print(f7["u"], f7["w"], "|", f4["u"], f4["w"])            # u, w [m/s] from both
assert np.allclose([f7["u"], f7["w"]], [f4["u"], f4["w"]], rtol=1e-12)   # identical
""")
D("D06", ref="7.28")
note("N25", "The line that produced (7.28)", r"""
$\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=-\frac{a\omega^2}{k}\frac{\cosh kH}{\sinh kH}\cos(kx-\omega t)\cong-g\eta
=-ag\cos(kx-\omega t)$ (D06 steps 2–4).
""")
nb.worked_example("H = 1 m, k = 1 rad/m", r"""
1. $kH=1$, $\tanh 1=0.762$.
2. $\omega^2=gk\tanh kH=9.81\times1\times0.762=7.47$ s⁻².
3. $\omega=2.73$ rad/s.
4. $T=2\pi/\omega=2.30$ s.
5. $\lambda=2\pi/k=6.28$ m and $c=\omega/k=2.73$ m/s.

Check with deep water: $\sqrt{gk}=3.13$ rad/s would be 15 % faster — the 1 m bottom slows this 6 m wave.
""")
remind([
    ("sympy", "`sp.symbols`, `sp.diff`, `sp.simplify`; an expression that simplifies to 0 is an identity (Ch. 1 P40)."),
    ("scipy.optimize.brentq", "the root of f(x) = 0 inside a bracket [a, b] where f changes sign (Ch. 3 P108)."),
    ("np.expm1 and cancellation near zero", "subtracting nearly equal numbers loses digits; rewriting the formula avoids it (Ch. 3 P107)."),
    ("overflow of cosh/sinh and the stable ratio", r"$\cosh x$ exceeds the largest float near x ≈ 710; a ratio such as $\cosh k(z+H)/\sinh kH$ is rewritten with decaying exponentials so that it never overflows (explained in the cell below)."),
])
nb.code(r"""
print("ω(k = 1, H = 1) =", ch07.omega_gravity(1.0, 1.0), "rad/s")      # (7.28): the tiny example
print("T(λ = 50 m, H = 10 m) =", ch07.period_from_wavelength(50.0, 10.0), "s;  T(λ = 156 m, deep) =",
      ch07.period_from_wavelength(156.0), "s")                          # (7.28) in its T–λ form
k = ch07.wavenumber_from_omega(2*np.pi/10, 30.0)                        # invert (7.28): T = 10 s on a 30 m shelf
print(f"k = {k:.6f} rad/m, λ = {2*np.pi/k:.1f} m")
s = ch07.surface_wave_sympy()                                           # D05–D06 redone in sympy (cached)
print(s["dispersion"], s["laplace_residual"], s["bottom_residual"], s["kinematic_residual"], s["regroup_identity"])
""", explain=r"""
1. `omega_gravity(k, H)` is the forward relation (7.28): 2.733 rad/s for the tiny example.
2. `period_from_wavelength` is its T–λ form $T=\sqrt{\frac{2\pi\lambda}{g}\coth\frac{2\pi H}{\lambda}}$ (D06 step 7): 6.14 s
   for λ = 50 m on H = 10 m (kH = 1.26), and 9.996 s for a 156 m wave in deep water.
3. `wavenumber_from_omega` inverts (7.28) with `brentq` between a deep-water and a shallow-water root that bracket the
   answer: T = 10 s on a 30 m shelf has λ = 137.3 m.
4. `surface_wave_sympy` redoes D05–D06 symbolically: it prints the dispersion relation and four zeros (Laplace, bottom,
   kinematic condition, and the regrouping of D05 step 9).
5. Numerics: cosh and sinh overflow near kH ≈ 710, so the library writes the ratio
   $\frac{\cosh k(z+H)}{\sinh kH}=e^{kz}\frac{1+e^{-2k(z+H)}}{1-e^{-2kH}}$ (D08 step 4) and uses `np.tanh`, which saturates
   safely at 1.
""")
remind([
    ("Newton's method", r"improve a root guess by $k\leftarrow k-f(k)/f'(k)$; it converges very fast near the root (Ch. 6 P152)."),
])
nb.check_agree(r"""
w0, H = 2*np.pi/10, 30.0                                  # T = 10 s swell on a 30 m shelf: ω [rad/s], depth [m]
k = w0**2/G                                               # start from the deep-water guess k = ω²/g [rad/m]
for i in range(5):                                        # Newton on f(k) = ω² − gk tanh kH
    f = w0**2 - G*k*np.tanh(k*H)                          # the residual of (7.28) [1/s²]
    fp = -(G*np.tanh(k*H) + G*k*H/np.cosh(k*H)**2)        # f′(k) by the product and chain rules
    k -= f/fp                                             # Newton step
    print(i, f"{k:.9f}")
assert np.allclose(k, ch07.wavenumber_from_omega(w0, H), rtol=1e-10)   # same root as the library's brentq
""")
nb.md(r"""
Newton needs three iterations for ten digits. **V5 cross-check** against two published explicit approximations of (7.28)
(Fenton & McKee 1990, error ≤ 1.7 % in wavelength; Guo 2002, ≤ 0.8 %), both exact in the deep and shallow limits:
""")
nb.code(r"""
w = 2*np.pi/10                                            # a 10 s swell [rad/s]
for H in (2.0, 30.0, 300.0):                              # shallow, intermediate, deep [m]
    exact = ch07.wavenumber_from_omega(w, H)*H           # kH from (7.28) by brentq
    fm, guo = ch07.fenton_mckee_kh(w, H), ch07.guo_kh(w, H)   # the two explicit formulas
    print(f"H = {H:5.0f} m: kH = {exact:.5f};  Fenton–McKee {fm/exact - 1:+.2%}, Guo {guo/exact - 1:+.2%}")
""", explain=r"""
Both formulas stay within their published bounds (Fenton–McKee within 0.9 %, Guo within 0.5 % here) and become exact in deep water — independent
evidence that our inverse of (7.28) is right.
""")
nb.figure(r"""
k = np.logspace(-3, 0, 300)                               # wavenumbers [rad/m]
fig, ax = plt.subplots(figsize=(6.5, 3.8))                # one log–log panel
blues = ["#93c5fd", "#3b82f6", "#1e3a8a"]                 # three shades of blue for three depths
for H, cl in zip((1.0, 10.0, 100.0), blues):              # three depths [m]
    ax.loglog(k, ch07.omega_gravity(k, H), color=cl, lw=2, label=f"H = {H:g} m")    # (7.28)
    ax.loglog(k, k*np.sqrt(G*H), ":", color=COLORS["muted"], lw=1)                  # shallow line k√(gH)
    ax.plot(1/H, ch07.omega_gravity(1/H, H), "o", color=cl, ms=5)                   # kH = 1 on each curve
ax.loglog(k, np.sqrt(G*k), "--", color=COLORS["muted"], lw=1.5, label="deep √(gk) (7.45)")   # deep asymptote
ax.set_ylim(1e-3, 5); ax.set_xlabel("wavenumber k [rad/m]"); ax.set_ylabel("ω [rad/s]")
ax.set_title("Each depth follows k√(gH) for long waves, then joins √(gk)", fontsize=10); ax.legend(fontsize=8)
savefig(fig, "ch07", "c03_dispersion"); plt.show()
""", see=r"""Three curves that start on straight lines of slope 1 (dotted, $k\sqrt{gH}$) and bend onto the common
slope-½ line (dashed, $\sqrt{gk}$); a dot marks kH = 1 on each.""",
    read=r"""Left of kH ≈ 0.3 the wave feels the bottom ($\omega\propto k$: all long waves move together); right of kH ≈ 3
it does not ($\omega\propto\sqrt k$).""",
    change=r"""…the depth were halved: each curve's bend moves to twice the k; the deep line does not move at all.""")
whatif(r"""
…the surface had tension σ? Every g in D06 becomes $g+\sigma k^2/\rho$,
$\omega=\sqrt{k\big(g+\frac{\sigma k^2}{\rho}\big)\tanh kH}$ *(7.56)* (C07). …there were a lighter fluid above instead of
air? gk is multiplied by $(\rho_2-\rho_1)/(\rho_2+\rho_1)$, $\omega=\varepsilon\sqrt{gk}$ *(7.95)* (C13).
""")

# ---- C04 ------------------------------------------------------------------------------------------------------------
nb.recap("R07", "Perturbation pressure", r"""
Split the pressure into its still-water part and the wave's part: $p'\equiv p+\rho gz$ *(7.30)* (in gauge pressure the
still water has p = −ρgz). Ch. 4 made the same split, $p=p_s+p'$, for the Boussinesq equations.
""", where="Ch. 4 §4.9")
core("C04", r"Phase speed: dispersion, deep water and shallow water", r"""
When does the depth matter — and why do long waves all travel together?
""", eqs=("7.29",))
problem(r"""
A tsunami crosses the Pacific in less than a day; the swell from a storm in the Southern Ocean needs about a week to reach
California. Both are gravity waves on the same water. The difference is the wavelength compared with the depth — we want
the rule, its two limits and the numbers.
""")
nb.md(r"""
#### The idea

| regime | condition (book) | tanh kH ≈ | phase speed | pressure below | dispersive? |
|---|---|---|---|---|---|
| deep | kH > 2 (H > 0.32λ) | 1 | $\sqrt{g/k}$ *(7.45)* | decays like $e^{kz}$ *(7.48)* | yes: c ∝ √λ |
| intermediate | 0.44 < kH < 2 | tanh kH | $\sqrt{\frac gk\tanh kH}$ *(7.29)* | $\frac{\cosh k(z+H)}{\cosh kH}$ *(7.31)* | yes |
| shallow | H < 0.07λ (kH < 0.44) | kH | $\sqrt{gH}$ *(7.49)* | hydrostatic $\rho g\eta$ *(7.52)* | no |

**Depth is invisible to a wave shorter than about three times the depth (H > 0.32λ); for longer waves it sets the speed.**
""")
note("N37", "The two limits of the hyperbolic functions", r"""
(primer P168) tanh x → 1 for large x ($\tanh2=0.964$ — the deep-water threshold); tanh x ≈ sinh x ≈ x and cosh x ≈ 1
for small x. The figure is our remake of Fig. 7.7 `N173`.
""")
nb.figure(r"""
xh = np.linspace(0, 2.3, 200)                             # the argument x = kH [–]
fig, ax = plt.subplots(figsize=(5.5, 3.2))                # one panel
ax.plot(xh, np.cosh(xh), color=COLORS["amber"], lw=2, label="cosh x")   # cosh starts at 1
ax.plot(xh, np.sinh(xh), color=COLORS["teal"], lw=2, label="sinh x")   # sinh starts at 0 with slope 1
ax.plot(xh, np.tanh(xh), color=COLORS["blue"], lw=2, label="tanh x")   # tanh flattens under 1
ax.plot(xh, xh, ":", color=COLORS["muted"], lw=1, label="y = x (small-x form)")   # the small-argument form
ax.axhline(1, color=COLORS["muted"], ls="--", lw=1)       # tanh's ceiling
ax.plot(2, np.tanh(2), "o", color=COLORS["orange"]); ax.text(1.55, 1.12, "tanh 2 = 0.964", color=COLORS["orange"], fontsize=8)   # the deep-water threshold kH = 2
ax.set_ylim(0, 3); ax.set_xlabel("x = kH [–]"); ax.set_ylabel("value [–]")   # axes
ax.set_title("Small x: sinh ≈ tanh ≈ x; large x: tanh → 1", fontsize=10); ax.legend(fontsize=8)   # title and legend
savefig(fig, "ch07", "c04_hyperbolic"); plt.show()   # save, then draw
""", see=r"""cosh starts at 1, sinh and tanh at 0 with slope 1; tanh flattens under the dashed line 1.""",
    read=r"""Left: sinh ≈ tanh ≈ x (shallow water); right: tanh ≈ 1 (deep water) — already 0.964 at x = 2.""",
    change=r"""…we plotted to x = 5: cosh and sinh become indistinguishable (both ≈ $e^x/2$).""")
remind([
    ("product rule", r"$d(uv)=u\,dv+v\,du$ (Ch. 1 P38); with the chain rule (P49) it differentiates $\lambda\tanh(2\pi H/\lambda)$ in D07."),
])
D("D07", ref="7.29")
nb.code(r"""
for lam in (50.0, 500.0):                                 # two wavelengths on 10 m of water [m]
    print(f"λ = {lam:.0f} m: c = {ch07.phase_speed(2*np.pi/lam, 10.0):.2f} m/s")   # (7.29): longer is faster
print(f"√(gH) = {np.sqrt(G*10.0):.2f} m/s (never exceeded)")                        # the shallow-water ceiling
""")
D("D08", ref="7.45")
note("N38, N42", "Deep water", r"""
$c=\sqrt{g/k}=\sqrt{g\lambda/2\pi}$ *(7.45)*, within 2 % once kH > 2 (H > 0.32λ, "deeper than a third of a
wavelength"); still dispersive — c grows like √λ. The pressure $p'=\rho gae^{kz}\cos(kx-\omega t)$ *(7.48)* has only
$e^{-\pi}=4.3\,\%$ left at depth λ/2, so a bottom pressure gauge in deep water is a low-pass filter — it sees swell only
when the water is shallow compared with λ.
""")
remind([
    ("Taylor series of tanh and of a square root", r"$\tanh x\approx x-x^3/3$ and $\sqrt{1-\varepsilon}\approx1-\varepsilon/2$; sympy's `series` produces them (Ch. 4 P117, Ch. 1 P26)."),
])
D("D09", ref="7.49")
note("N43", "Shallow water", r"""
$c=\sqrt{gH}$ *(7.49)* — the same speed for every long wave (nondispersive), within 3 % when H < 0.07λ; it matches the
control-volume bore speed of Ch. 4 Example 4.3 in the small-step limit (Ch. 4 §4.4). Below: a 200 km tsunami on 4 km of
ocean.
""")
nb.code(r"""
print(ch07.phase_speed(2*np.pi/2e5, 4000.0), "m/s from (7.29)")          # λ = 200 km, H = 4 km
print(np.sqrt(G*4000.0), "m/s = √(gH), (7.49)")                         # the shallow-water limit
print(ch04.bore_speed(4000.0, 4000.001, g=G), "m/s: Ch. 4 bore speed for a 1 mm step")   # Example 4.3, vanishing step
""", explain=r"""
(7.29) gives 197.57 m/s, the shallow limit 198.09 m/s (0.3 % faster) and Ch. 4's bore speed for a vanishing step
198.09 m/s: the momentum balance of Ch. 4 and the dispersion relation agree for long, small waves.
""")
note("N46", "Shallow waves are hydrostatic", r"""
In shallow water the pressure at every depth is just the weight of the extra water above — the hydrostatic shallow-water
model of Ch. 13 (tides, tsunamis, Kelvin waves).
""", equation=EQ["7.52"], ref="7.52")
note("N26", "The pressure under a wave", r"""
From $\frac{\partial\phi}{\partial t}+\frac p\rho+gz\cong0$ *(7.20)* with $p'\equiv p+\rho gz$ *(7.30)* and the potential
$\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ *(7.26)*, then
$\omega=\sqrt{gk\tanh kH}$ *(7.28)* to remove ω² (two substitutions; the book writes them in one line):
$p'=-\rho\frac{\partial\phi}{\partial t}=\rho\frac{a\omega^2}{k}\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t)$, which is
(7.31) below. It decays with depth at a rate set by k: short waves are felt only near the surface.
""", equation=EQ["7.31"], ref="7.31")
nb.worked_example("swell, a tsunami and a bottom gauge", r"""
1. T = 10 s in the open ocean (deep): $\lambda=gT^2/2\pi=9.81\times100/6.283=156$ m and $c=\lambda/T=15.6$ m/s.
2. The same swell on a 30 m shelf: solve (7.28) → λ = 137 m, c = 13.7 m/s (kH = 1.37: it has started to feel the bottom).
3. A tsunami, H = 4 km: $c=\sqrt{gH}=\sqrt{39\,240}=198$ m/s ≈ 713 km/h — the Pacific (≈ 15 000 km) in ≈ 21 h.
4. λ = 50 m, a = 1 m on H = 10 m: $\rho ga=9.81$ kPa at z = 0; at the bottom divide by $\cosh kH=\cosh1.257=1.90$ →
   5.17 kPa.
""")
note("N39", "Ocean scales", r"""
Wind waves of T ≈ 10 s are ≈ 156 m long, so over a 100 m shelf and a 4 km ocean they are deep-water waves until close to
shore; tides and tsunamis (hundreds of km) are always shallow-water waves. The code solves $\omega=\sqrt{gk\tanh kH}$
*(7.28)* for λ at fixed T = 10 s and classifies each case with the book's thresholds, printing both limit errors.
""")
nb.code(r"""
for H in (np.inf, 4000.0, 100.0, 30.0, 10.0):             # the same 10 s swell over five depths [m]
    lam = ch07.wavelength_from_period(10.0, H)            # λ from (7.28) with T fixed [m]
    r = ch07.depth_regime(2*np.pi/lam, H)                 # regime and the errors of √(g/k) and √(gH)
    print(f"H = {H:>6} m: λ = {lam:6.1f} m, kH = {r['kH']:7.2f}, {r['regime']:12s} "
          f"error of √(g/k) {r['deep_error']:.2%}, of √(gH) {r['shallow_error']:.2%}")
for lam in (5.0, 50.0, 500.0):                            # three wavelengths on H = 10 m
    k = 2*np.pi/lam                                       # wavenumber [rad/m]
    print(f"λ = {lam:5.0f} m: p′/ρga at z = 0, −5, −10 m:", ch07.pressure_response(k, np.array([0.0, -5.0, -10.0]), 10.0))
""", explain=r"""
1. λ = 156.1, 156.1, 156.0, 137.3 and 92.4 m: deep, deep, deep (kH = 4.03, error 0.03 %), intermediate (kH = 1.37) and
   intermediate (kH = 0.68). A 10 s swell only starts to feel the bottom within a few tens of metres of depth.
2. The cosh ratio of (7.31) at three depths: a 5 m wave leaves 0.2 % at 5 m depth and nothing at the bottom; a 50 m wave
   still presses 53 % at the bottom; a 500 m wave presses 99 % at every depth — hydrostatic.
""")
nb.figure(r"""
kH = np.logspace(-2, 1, 400)                              # depth ratio [–]
ratio = np.sqrt(np.tanh(kH)/kH)                           # c/√(gH) from (7.29)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6))      # speed | pressure
a1.semilogx(kH, ratio, color=COLORS["ink"], lw=2, label="(7.29)")
a1.semilogx(kH, 1/np.sqrt(kH), "--", color=COLORS["muted"], label="deep limit 1/√(kH)")   # (7.45)/√(gH)
a1.axhline(1, ls=":", color=COLORS["muted"], label="shallow limit 1")                     # (7.49)
a1.axvspan(2, 10, color=COLORS["blue"], alpha=0.08); a1.axvspan(0.01, 0.44, color=COLORS["amber"], alpha=0.1)
for xv in (0.44, 2.0):                                    # the book's thresholds
    a1.axvline(xv, color=COLORS["orange"], lw=1); a1.text(xv*1.05, 1.55, f"kH = {xv}", fontsize=7, color=COLORS["orange"])
a1.set_ylim(0, 1.7); a1.set_xlabel("kH [–]"); a1.set_ylabel(r"$c/\sqrt{gH}$ [–]"); a1.legend(fontsize=7, loc="lower left")
a1.set_title("speed: shallow band (sand), deep band (blue)", fontsize=9)
zH = np.linspace(-1, 0, 200)                              # z/H [–]
for kh, cl in zip((0.3, 1.0, 3.0), ("#fcd34d", COLORS["amber"], "#92400e")):   # three depth ratios (amber shades)
    a2.plot(ch07.pressure_response(kh, zH, 1.0), zH, color=cl, lw=2, label=f"kH = {kh}")   # (7.31) with k = kH, H = 1
    a2.plot(np.exp(kh*zH), zH, "--", color=cl, lw=1)      # e^{kz} (7.48) for the same k
a2.axvline(1, ls=":", color=COLORS["muted"], label="hydrostatic (7.52)")
a2.set_xlabel(r"$p'/\rho g a$ under a crest [–]"); a2.set_ylabel("z/H [–]"); a2.legend(fontsize=7)
a2.set_title("pressure: from vertical lines to fast decay", fontsize=9)
fig.suptitle("Where each limit holds", fontweight="bold")
savefig(fig, "ch07", "c04_limits"); plt.show()
print("error of √(g/k) at kH = 2:", f"{1 - np.sqrt(np.tanh(2.0)):.2%};", "error of √(gH) at kH = 0.44:",
      f"{1 - np.sqrt(np.tanh(0.44)/0.44):.2%}")           # the book's '2 %' and '3 %'
""", see=r"""Left: the speed curve leaves the shallow value near kH ≈ 0.4 and joins the deep curve near kH ≈ 2; right: the
pressure profiles go from nearly vertical lines (kH = 0.3) to a fast decay (kH = 3), each close to its dashed
$e^{kz}$ ghost in deep water.""",
    read=r"""Read the error of a limit as the gap between the solid curve and the ghost: 3 % at kH = 0.44, 1.8 % at kH = 2
(printed under the figure).""",
    change=r"""…the depth doubled at fixed λ: every point slides right by log 2 — the wave becomes "deeper".""")
remind([
    ("slider_figure", "computes every slider position in Python up front, so the plotly figure keeps working on the published page (Ch. 1 P17)."),
    ("live widgets", "ipywidgets sliders that re-run a function while a kernel runs; frozen on the page (Ch. 1 P47)."),
    ("show_viz", "embeds a chapter explainer in the notebook (Ch. 1 P18)."),
    ("explainer tabs (Walkthrough, Explore, Explain, Derivation, Equations, Code, Check)", "start with the Walkthrough, play in Explore, read every number worked out in Explain, step through the Derivation, finish with Check yourself (Ch. 1 P18)."),
])
nb.plotly(r"""
lam = np.logspace(-1, 6, 150)                             # wavelengths 0.1 m … 1000 km [m]
fig = slider_figure(lambda H: {"finite depth (7.29)": (lam, ch07.phase_speed(2*np.pi/lam, H)),   # for each H: c(λ) from (7.29) …
                               "deep limit (7.45)": (lam, np.sqrt(G*lam/(2*np.pi))),   # … the deep line √(gλ/2π) …
                               "shallow limit (7.49)": (lam, np.sqrt(G*H) + 0*lam)},   # … and the shallow plateau √(gH)
                    "H", np.logspace(np.log10(0.5), np.log10(5000), 20 if not FAST else 10), unit="m",   # slider: depth 0.5 m … 5 km
                    xlabel="wavelength λ [m]", ylabel="c [m/s]", title="Depth only matters for waves longer than about 3H")   # labels and message
logaxes(fig, x=(0.1, 1e6), y=(0.3, 400))                  # both axes logarithmic, fixed
recolor(fig, {"finite depth (7.29)": COLORS["ink"], "deep limit (7.45)": COLORS["muted"],   # colours by meaning …
              "shallow limit (7.49)": COLORS["muted"]}, dashes={"deep limit (7.45)": "dash", "shallow limit (7.49)": "dot"})   # … deep dashed, shallow dotted
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** Drag H: the black curve follows the dashed deep line for short waves and the dotted shallow plateau $\sqrt{gH}$ for long
ones; the corner between them sits where kH is between about 0.44 and 2, i.e. λ between about 3H (H > 0.32λ: deep)
and 14H (H < 0.07λ: shallow). Every slider position was computed up front, so the figure works on the page.

**What would change if…** g were smaller (a lighter planet): every curve would drop by the same factor $\sqrt g$, but the
corner would stay at λ ≈ 3H–14H, because it depends only on kH.
""")
nb.plotly(r"""
z = np.linspace(-20, 0, 80)                              # depth below the still surface, H = 20 m [m]
fig = slider_figure(lambda lam: {"cosh ratio (7.31)": (ch07.pressure_response(2*np.pi/lam, z, 20.0), z),   # for each λ: p′/ρga from (7.31) …
                                 "hydrostatic (7.52)": (1 + 0*z, z),   # … the hydrostatic line …
                                 "deep e^{kz} (7.48)": (np.exp(2*np.pi/lam*z), z)},   # … and the deep-water e^{kz}
                    "λ", np.logspace(0, 3, 20 if not FAST else 10), unit="m", xlabel="p′/ρga [–]", ylabel="z [m]",   # slider: λ from 1 m to 1 km
                    title="When a bottom gauge sees the wave (H = 20 m)", xrange=(0, 1.05), yrange=(-20, 0))   # fixed axes
recolor(fig, {"cosh ratio (7.31)": COLORS["amber"], "hydrostatic (7.52)": COLORS["muted"], "deep e^{kz} (7.48)": COLORS["muted"]},   # pressure amber, limits grey
        dashes={"hydrostatic (7.52)": "dot", "deep e^{kz} (7.48)": "dash"})   # limit styles
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** Short waves (left end of the slider) leave the bottom untouched; from λ ≈ 100 m on, the amber profile straightens toward
the hydrostatic line: the whole column feels the wave.

**What would change if…** the gauge sat on a 100 m deep bottom: the amber profile would need waves about five times longer before it
straightened, so the gauge would miss most wind waves and see only swell and tides.
""")
nb.live(r"""
_ = live(lambda lam=100.0, H=20.0: print({k_: (round(v, 4) if isinstance(v, float) else v)            # every number of E1
                                     for k_, v in ch07.dispersion_state(lam, H).items()}),
     lam=(1.0, 2000.0, 1.0), H=(0.5, 4000.0, 0.5))       # sliders: wavelength [m], depth [m]
""", explain=r"""
`dispersion_state(λ, H)` returns k, kH, ω, T, c, c_g, both limits and their errors, the regime, the fraction of the
surface pressure felt at the bottom, and how long crests and energy need to cross 1000 km. The plotly slider above
carries the same idea on the page.
""")
explainer("dispersion_relation", "Why long waves outrun short ones",
          r"""Three linked pictures of one number: the tank with a crest marker moving at c and the pressure fading with
depth, the c(λ) curve with its two limits, and the pressure profile. Dragging λ across kH ≈ 0.3…3 moves all three at
once — a static figure can show only one.""", "",
          ["Start at the *wind swell* preset, then drag H down to 10 m: watch λ and c shrink while T stays 10 s.",
           "Choose *tsunami*: the status says shallow — $c=\\sqrt{gH}$ *(7.49)* whatever λ you pick.",
           "Click the tank at the bottom for the pressure arithmetic; compare λ = 50 m and λ = 500 m.",
           "Open the Derivation tab (D06) and step to 'cancel the cosine'."])
confusion(r"""
"deep water" is not a depth in metres. A 10 m pond is deep for 1 m ripples and shallow for a 1 km seiche; the ocean is
deep for swell and shallow for tides. Always compare H with λ.
""")
whatif(r"""
…the Earth rotated fast enough to matter (periods of hours)? Long waves in Ch. 13 add the Coriolis frequency f:
$\omega^2=f^2+gHk^2$ (Poincaré waves) — this section's $c=\sqrt{gH}$ *(7.49)* is the limit f → 0.
""")

# ---- C05 ------------------------------------------------------------------------------------------------------------
nb.recap("R08", "Path lines", r"""
A fluid particle's path $x_p(t)$ obeys $\frac{dx_p}{dt}=u(x_p,z_p,t),\ \frac{dz_p}{dt}=w(x_p,z_p,t)$ *(7.32)* — the
path-line equation $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ *(3.8)*. Ch. 3 Example 3.1 followed a particle in an oscillating
flow and found a closed circle, as we are about to under a wave.
""", where="Ch. 3 §3.2")
nb.code(r"""
r = ch03.example_3_1(0.5)                                 # Ch. 3 Example 3.1 at t′ = 0.5 s (ξ₀ = 1 m, ω = 1 rad/s)
print("path line of Example 3.1: a circle of radius", r["radius"], "m")   # oscillating flow → closed loop
""", explain=r"""
Ch. 3's example of a flow that oscillates in time: its path lines are closed loops, although the streamlines at any
instant are straight — the same distinction C05 needs for the water under a wave.
""")
core("C05", r"Particle orbits: closed ellipses", r"""
If the wave moves but the water does not travel, what path does a water parcel follow?
""", eqs=("7.36",))
problem(r"""
Throw a cork on swell and release a neutrally buoyant float a few metres down: both go round and come back. Divers feel
the surge forward under a crest and backward under a trough. We want the shape of those loops at every depth — the
picture behind Stokes drift (C12), ocean floats (Ch. 13) and sediment moved to and fro on the sea bed.
""")
nb.md(r"""
#### The idea

```
deep water:     circles of radius a e^{kz₀}, shrinking fast with depth
intermediate:   ellipses, flatter as you go down
shallow water:  flat ellipses of the same width a/kH at every depth, height → 0 at the bottom
all parcels of one vertical column are at the same point of their loops (in phase); clockwise for a wave to the right
```
""")
note("N27", "The path lines under the wave", r"""
With the velocities $u=a\omega\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t),\ w=a\omega\frac{\sinh k(z+H)}{\sinh kH}
\sin(kx-\omega t)$ *(7.27)* inserted, the path lines are nonlinear in the unknown position:
""", equation=EQ["7.33"], ref="7.33")
remind([
    ("multivariable Taylor expansion", r"$u(x_0+\xi,z_0+\zeta)\approx u(x_0,z_0)+\xi u_x+\zeta u_z$, error ∝ |(ξ, ζ)|² (Ch. 3 P98)."),
], lead="antiderivatives of sin and cos are worked in D10 steps 6–7")
D("D10", ref="7.35a")
note("N28, N29", "Frozen arguments and the excursions", r"""
Freezing the right sides at the mean position gives $\frac{d\xi}{dt}\cong a\omega\frac{\cosh k(z_0+H)}{\sinh kH}
\cos(kx_0-\omega t)$ *(7.34a)* and $\frac{d\zeta}{dt}\cong a\omega\frac{\sinh k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)$
*(7.34b)* (D10 steps 4–5), and one integration gives the excursions $\xi\cong-a\frac{\cosh k(z_0+H)}{\sinh kH}\sin(kx_0-
\omega t)$ *(7.35a)*, $\zeta\cong a\frac{\sinh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)$ *(7.35b)* — purely oscillatory, so
the assumption "the parcel stays near $(x_0,z_0)$" is self-consistent.
""")
remind([
    ("ellipse, semi-axes, foci", r"an ellipse $x^2/A^2+y^2/B^2=1$ with A ≥ B has its foci $\sqrt{A^2-B^2}$ from the centre; a linear map of a circle is an ellipse (Ch. 3 P104)."),
], lead="sin² + cos² = 1 is D11 step 3")
D("D11", ref="7.36")
note("N40, N41", "Deep water", r"""
(kH > 2): $\xi\cong-ae^{kz_0}\sin(kx_0-\omega t),\ \zeta\cong ae^{kz_0}\cos(kx_0-\omega t)$ *(7.46)* — circles whose radius
halves every $\ln2/k=0.11\lambda$ of depth; the velocity $u=a\omega e^{kz}\cos(kx-\omega t),\ w=a\omega e^{kz}\sin(kx-
\omega t)$ *(7.47)* turns clockwise at ω with constant size $a\omega e^{kz}$.
""")
note("N44, N45", "Shallow water", r"""
(kH ≪ 1): $\xi\cong-\frac a{kH}\sin(kx_0-\omega t),\ \zeta\cong a\big(1+\frac{z_0}H\big)\cos(kx_0-\omega t)$ *(7.50)* —
flat ellipses of width a/kH at every depth; $u=\frac{a\omega}{kH}\cos(kx-\omega t),\ w=a\omega\big(1+\frac zH\big)
\sin(kx-\omega t)$ *(7.51)*, so $w/u\lesssim kH$ (0.2 for kH = 0.2).
""")
nb.worked_example("orbits under a 20 m wave in 50 m of water", r"""
1. $k=2\pi/20=0.314$ rad/m, $kH=15.7$: deep. Take a = 0.1 m, so ka = 0.031 — small, as linear theory needs.
2. Radius $ae^{kz_0}$: 10 cm at the surface.
3. At $z_0=-2$ m: $e^{-0.628}=0.533$ → 5.3 cm.
4. At $z_0=-10$ m (half a wavelength): $e^{-\pi}=0.043$ → 4 mm.
5. The bottom (−50 m) does not feel the wave: $e^{-15.7}\approx10^{-7}$.

Every radius is proportional to a; the code below uses a = 1 m as a unit amplitude, so its lengths are "per metre of a".
""")
remind([
    ("scipy.integrate.solve_ivp", "adaptive Runge–Kutta that picks its own steps to meet a tolerance (Ch. 1 P31)."),
    ("solve_ivp options and DOP853 tolerances", "`method=\"DOP853\"` is an 8th-order scheme, `rtol`/`atol` set the accuracy and `t_eval` the output times (Ch. 3 P94); tight tolerances make a tiny drift trustworthy."),
])
nb.code(r"""
for kH in (3.0, 1.0, 0.3):                                # k = 1 rad/m, so H = kH [m]; a = 1 m
    ax_ = ch07.orbit_semi_axes(np.array([0.0, -kH/2]), 1.0, 1.0, kH)   # (7.36) at the surface and at mid-depth
    print(f"kH = {kH}: A = {ax_['A']}, B = {ax_['B']} m; focal half-distance a/sinh kH = {ax_['focal_half']:.4f} m")
T1 = 2*np.pi/ch07.omega_gravity(1.0, 1.0)                 # period for k = 1 rad/m, H = 1 m [s]
tt = np.linspace(0, 8*T1, 801)                            # eight whole periods [s]
p_lin = ch07.particle_path(0.0, -0.5, tt, 0.05, 1.0, 1.0, model="linear")   # (7.34): frozen at the mean position
p_ex = ch07.particle_path(0.0, -0.5, tt, 0.05, 1.0, 1.0, model="exact")     # (7.33): the true path line
print(f"after 8 periods: linear path moved {p_lin['x'][-1] - p_lin['x'][0]:.1e} m, exact path {p_ex['x'][-1] - p_ex['x'][0]:.4f} m")
""", explain=r"""
1. The semi-axes of (7.36) at the surface and at mid-depth: for kH = 3, A = 1.005 and B = 1.000 m at the surface (almost
   a circle), 0.235 and 0.213 m at mid-depth; kH = 1: 1.313/1.000 and 0.960/0.443; kH = 0.3: 3.433/1.000 and 3.321/0.494.
   The focal half-distance a/sinh kH — 0.0998, 0.851, 3.284 m — is the same at every depth.
2. `particle_path` integrates the frozen-argument equations (7.34) or the exact path lines (7.33) with `solve_ivp`
   (DOP853 and tight tolerances, so the tiny drift is not numerical error). After eight whole periods the linear parcel is
   back where it started; the exact one has crept forward by several centimetres — the Stokes drift of C12.
""")
remind([
    ("RK4 by hand", "four slope samples per step weighted 1-2-2-1; error ∝ Δt⁴ (Ch. 3 P95)."),
])
nb.check_agree(r"""
a, k, H, x0, z0 = 0.05, 1.0, 1.0, 0.0, -0.5               # wave and mean position [m, rad/m, m, m, m]
w = ch07.omega_gravity(k, H); T = 2*np.pi/w                # ω [rad/s] and period [s]
Ca = np.cosh(k*(z0 + H))/np.sinh(k*H); Sa = np.sinh(k*(z0 + H))/np.sinh(k*H)   # the depth ratios of (7.34)
rhs = lambda t: np.array([a*w*Ca*np.cos(k*x0 - w*t), a*w*Sa*np.sin(k*x0 - w*t)])   # (7.34a, b): dξ/dt, dζ/dt [m/s]
dt, n = T/400, 400                                        # 400 RK4 steps over one period
y = np.array(ch07.orbit_linear(x0, z0, 0.0, a, k, H))     # start on the orbit: (ξ, ζ) at t = 0 [m]
mine = [y.copy()]                                          # the RK4 path
for i in range(n):                                         # classic RK4 (1-2-2-1)
    t = i*dt
    k1 = rhs(t); k2 = rhs(t + dt/2); k3 = rhs(t + dt/2); k4 = rhs(t + dt)   # slopes (the right side depends on t only)
    y = y + dt*(k1 + 2*k2 + 2*k3 + k4)/6
    mine.append(y.copy())
mine = np.array(mine)                                      # (n + 1, 2)
lib = np.array(ch07.orbit_linear(x0, z0, np.arange(n + 1)*dt, a, k, H)).T   # (7.35) at the same times
print("largest difference:", np.abs(mine - lib).max(), "m")
assert np.allclose(mine, lib, atol=1e-6)                   # RK4 draws the same ellipse as the closed form
""")
nb.md("RK4 on the frozen-argument equations (7.34) draws the same ellipse as the closed form (7.35).")
nb.figure(r"""
fig, axs = plt.subplots(1, 4, figsize=(10, 3.4), gridspec_kw=dict(width_ratios=[0.8, 1, 1, 1]))
s = np.linspace(0, 2*np.pi, 200)                          # parameter round an ellipse
A0, B0 = 1.0, 0.6                                         # (a) a schematic orbit
axs[0].plot(A0*np.cos(s), B0*np.sin(s), color=COLORS["teal"], lw=2)
axs[0].plot(0, 0, "k+", ms=10); axs[0].text(-0.75, 0.08, "(x₀, z₀)", fontsize=8)
ph = 0.9                                                  # a phase along the orbit
px, pz = A0*np.cos(ph), B0*np.sin(ph)                     # the parcel's position relative to the mean
axs[0].annotate("", xy=(px, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["ink"]))   # ξ
axs[0].annotate("", xy=(px, pz), xytext=(px, 0), arrowprops=dict(arrowstyle="->", color=COLORS["ink"]))  # ζ
axs[0].text(px/2, -0.2, "ξ", fontsize=10); axs[0].text(px + 0.05, pz/2, "ζ", fontsize=10)
axs[0].annotate("", xy=(0.25, B0), xytext=(-0.25, B0), arrowprops=dict(arrowstyle="->", color=COLORS["teal"], lw=2))   # clockwise
axs[0].set_title("(a) one orbit, clockwise", fontsize=9); axs[0].set_aspect("equal"); axs[0].axis("off")
for axp, kH, tag in zip(axs[1:], (3.0, 1.0, 0.3), "bcd"):  # three depth ratios, k = 1 rad/m
    H = kH; a = 0.08*H                                     # depth [m] and an amplitude drawn at 8 % of the depth [m]
    xw = np.linspace(-1.2*H - 0.4, 1.2*H + 0.4, 200)        # window [m]
    axp.plot(xw, a*np.cos(xw), color=COLORS["blue"], lw=2)  # the surface at t = 0
    axp.fill_between(xw, -1.06*H, -H, color="#c8b99a")     # bottom
    z0s = np.linspace(-0.9*H, 0.0, 5)                      # five mean depths [m]
    ax_ = ch07.orbit_semi_axes(z0s, a, 1.0, H)             # (7.36) semi-axes at those depths
    orbit_ghosts(axp, 0.0, z0s, ax_["A"], ax_["B"])        # the orbits (teal)
    axp.set_aspect("equal"); axp.set_ylim(-1.08*H, 0.25*H); axp.set_xlim(xw[0], xw[-1])
    axp.set_title(f"({tag}) kH = {kH}", fontsize=9); axp.set_xlabel("x [m]")
axs[1].set_ylabel("z [m]")
fig.suptitle("Circles in deep water, flat ellipses in shallow water", fontweight="bold", y=0.9)
savefig(fig, "ch07", "c05_orbits"); plt.show()
""", see=r"""(a) one orbit with the mean position $(x_0,z_0)$ and the excursion (ξ, ζ), travelled clockwise; (b)–(d) the
orbits at five depths for kH = 3, 1, 0.3 (our remake of Figs. 7.3 and 7.4 `N169`, `N170`): circles that shrink fast
with depth, ellipses, flat ellipses reaching the bottom.""",
    read=r"""In (d) the horizontal axis is almost the same at every depth (a/kH) and the vertical axis falls linearly to zero
at the bottom, as $\zeta\cong a\big(1+\frac{z_0}H\big)\cos(kx_0-\omega t)$ *(7.50)* says.""",
    change=r"""…a doubled: every orbit doubles, the shapes do not change (linear theory) — until ka ~ 0.3, when the loops open
(C12).""")
note("N30", "The stream function of the wave", r"""
Integrating $u=\partial\psi/\partial z$ in z (and checking $-\partial\psi/\partial x=w$) gives the stream function below
(the convention of Ch. 4 §4.3 with y → z): ψ = 0 on the bottom and on the vertical lines below the zeros of
η. Under a crest the water moves forward at every depth, under a trough backward. The code differentiates ψ numerically
with Ch. 4's tool and gets the velocities (7.27) back.
""", equation=EQ["7.37"], ref="7.37")
nb.code(r"""
psi = lambda x, z: ch07.wave_fields(x, z, 0.0, 0.05, 1.0, H=2.0)["psi"]   # (7.37) at t = 0 [m²/s]
u, w = streamfunction.velocity_from_streamfunction_2d(psi, 0.4, -0.7)     # Ch. 4: u = ∂ψ/∂z, w = −∂ψ/∂x (y ↦ z)
f = ch07.wave_fields(0.4, -0.7, 0.0, 0.05, 1.0, H=2.0)                    # (7.27) directly
print(u, w, "|", f["u"], f["w"])                                          # equal to ~1e-10 m/s
""")
nb.md(r"**What does this show?** Differentiating the stream function " + E("7.37") + r" numerically (u = ∂ψ/∂z, w = −∂ψ/∂x, Ch. 4) gives u ≈ 0.0770 m/s and w ≈ 0.0280 m/s at (x, z) = (0.4, −0.7) m — the same numbers as the velocity formulas " + E("7.27") + r" to about $10^{-12}$ m/s. So ψ really is the stream function of this wave.")
nb.figure(r"""
k, H, a = 1.0, 1.0, 0.05                                  # kH = 1, amplitude 5 cm [rad/m, m, m]
xg, zg = np.linspace(0, 2*np.pi, 161 if not FAST else 81), np.linspace(-H, 0, 81 if not FAST else 41)
XX, ZZ = np.meshgrid(xg, zg)                              # grid over one wavelength and the depth
f = ch07.wave_fields(XX, ZZ, 0.0, a, k, H)                # ψ, u, w at t = 0
fig, ax = plt.subplots(figsize=(7, 3.2))                  # one panel
levels = np.linspace(-np.abs(f["psi"]).max(), np.abs(f["psi"]).max(), 13)   # ψ levels [m²/s]
ax.contour(XX, ZZ, f["psi"], levels=levels, colors=[COLORS["teal"]], linewidths=1)   # streamlines
ax.contour(XX, ZZ, f["psi"], levels=[0.0], colors=[COLORS["ink"]], linewidths=2.5)   # ψ = 0 (bold)
sl = (slice(None, None, 8), slice(None, None, 16))        # a sparse set of arrows
ax.quiver(XX[sl], ZZ[sl], f["u"][sl], f["w"][sl], color=COLORS["muted"], width=0.004)   # (u, w)
ax.plot(xg, 3*f["eta"][-1], color=COLORS["blue"], lw=2)   # the surface (height ×3 for visibility)
ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.set_ylim(-H, 0.2)
ax.set_title("Streamlines of a progressive wave (kH = 1, t = 0)", fontsize=10)
savefig(fig, "ch07", "c05_streamlines"); plt.show()
""", see=r"""Loops of ψ that start and end on the surface, forward flow under the crest (x = 0, 2π) and backward flow
under the trough (x = π); bold ψ = 0 lines along the bottom and below the zeros of η (our remake of Fig. 7.5 `N171`).""",
    read=r"""Streamlines are a snapshot of the velocity field; the particles themselves go round the orbits of the previous
figure.""",
    change=r"""…t advanced by T/4: the whole pattern slides a quarter wavelength right (it moves with the wave).""")
remind([
    ("animate and show_animation", "`update(i)` moves the artists already drawn for frame i; `player=\"video\"` plays smoothly (Ch. 1 P16)."),
])
nb.animation(r"""
fig, axs = plt.subplots(1, 3, figsize=(10, 3.2))          # kH = 3, 1, 0.3 side by side
k = 1.0                                                   # wavenumber [rad/m]; λ = 2π m
arts = []                                                 # (surface, dots, column line, parameters) per panel
for axp, kH in zip(axs, (3.0, 1.0, 0.3)):
    H = kH; a = 0.06*H                                     # depth and amplitude [m]
    w = ch07.omega_gravity(k, H)                           # ω [rad/s]
    X0, Z0 = np.meshgrid(np.linspace(0.3, 2*np.pi - 0.3, 7), np.linspace(-0.85*H, -0.02*H, 5))   # 5 × 7 mean positions
    ax_ = ch07.orbit_semi_axes(Z0, a, k, H)                # semi-axes (7.36)
    orbit_ghosts(axp, X0.ravel(), Z0.ravel(), ax_["A"].ravel(), ax_["B"].ravel(), color=COLORS["teal"], lw=0.4, alpha=0.4)
    xs = np.linspace(0, 2*np.pi, 200)                      # surface points [m]
    (surf,) = axp.plot(xs, a*np.cos(k*xs), color=COLORS["blue"], lw=2)
    dots = axp.scatter(X0.ravel(), Z0.ravel(), s=10, color=COLORS["teal"], zorder=3)
    (col,) = axp.plot(X0[:, 3], Z0[:, 3], color=COLORS["ink"], lw=1)   # one vertical column joined by a line
    axp.fill_between(xs, -1.06*H, -H, color="#c8b99a")
    axp.set_xlim(0, 2*np.pi); axp.set_ylim(-1.08*H, 0.3*H); axp.set_title(f"kH = {kH}", fontsize=9); axp.set_xticks([])
    arts.append((surf, dots, col, X0, Z0, a, H, w, xs))
nf = 36 if not FAST else 18                                # frames over one period of each panel


def update(i):                                             # frame i: advance every panel by i/nf of its own period
    for surf, dots, col, X0, Z0, a, H, w, xs in arts:
        t = i/nf*2*np.pi/w                                 # time in this panel [s]
        surf.set_ydata(a*np.cos(k*xs - w*t))               # (7.2): the surface
        xi, ze = ch07.orbit_linear(X0, Z0, t, a, k, H)     # (7.35): every parcel's excursion
        dots.set_offsets(np.c_[(X0 + xi).ravel(), (Z0 + ze).ravel()])
        col.set_data(X0[:, 3] + xi[:, 3], Z0[:, 3] + ze[:, 3])   # the column moves as a whole
    return []   # nothing to blit


fig.suptitle("Parcels go round their orbits in phase down each column", fontsize=10)
show_animation(animate(update, frames=nf, fig=fig, interval=60))   # MP4 in the page (or a frame player)
""")
see_read_change(r"""The surface (blue) moves right; every parcel (teal) goes round its faint orbit ghost; the black line joins one
column of parcels — the `A1` animation.""",
                r"""The column stays straight and swings as a whole: all its parcels are at the same point of their loops (in
phase). In shallow water (right) the swing is almost purely horizontal.""",
                r"""…the amplitude were not small? The loops would no longer close — the column would lean forward period by
period (C12).""")
whatif(r"""
…the amplitude were not small? The parcel is slightly further forward at the top of its loop, where the forward velocity
is larger, than it is backward at the bottom: the loop does not close. That second-order creep is the Stokes drift of C12.
""")

# ---- C06 ------------------------------------------------------------------------------------------------------------
core("C06", r"Wave energy, equipartition and the energy flux", r"""
How much energy does a wave hold per square metre of sea — and how fast is it delivered to a beach?
""", eqs=("7.42", "7.44"))
problem(r"""
A storm thousands of kilometres away sends swell that breaks on a beach and keeps a surfer moving. The energy crossed the
ocean without the water crossing it. We want how much energy a wave stores (moving water and lifted water) and the power
it carries per metre of crest — the numbers wave-energy engineers and swell forecasters use.
""")
nb.md(r"""
#### The idea

Kinetic energy lives in the orbiting water under the surface; potential energy in the humps and hollows of the surface
itself. Average both over a wavelength and a period:

```
E_k = ½ρg·mean(η²)   (moving water)        E_p = ½ρg·mean(η²)   (lifted water)
E   = ρg·mean(η²) = ½ρga²                  F = E × (a speed)  —  which speed?
```

Two kinds of average appear: the overbar $\overline{\eta^2}=\frac1\lambda\int_0^\lambda\eta^2dx$ over a wavelength and
⟨·⟩ over a period; for a sinusoid both give ⟨cos²⟩ = ½ (Ch. 6 P151).
""")
note("N31", "The kinetic energy per unit horizontal area", r"""
Depth-integrated and wavelength-averaged, it is (7.38); with (7.27) it splits into an x-part and a z-part (D12 step 3).
The z-integral stops at 0, not at η: the slab in between changes E only at third order in a.
""", equation=EQ["7.38"], ref="7.38")
remind([
    ("integrals of sin² and cos² over a period", "over a full period sin and cos average to 0, sin² and cos² to ½ (Ch. 6 P151)."),
    ("definite integral", r"$\int_a^b f\,dz$ adds pieces f dz between the limits (Ch. 1 P27); the integrals of cosh² and sinh² are done in D12 steps 5–6."),
])
D("D12", ref="7.39")
note("N32", "Kinetic energy", r"Gravity appears through the dispersion relation $\omega^2=gk\tanh kH$ *(7.28)*, which trades the velocities for g:", equation=EQ["7.39"], ref="7.39")
note("N33", "The potential energy", r"""
is the work to deform the flat surface:
$E_p=\frac{\rho g}{\lambda}\int_0^\lambda\!\!\int_{-H}^\eta z\,dz\,dx-\frac{\rho g}{\lambda}\int_0^\lambda\!\!\int_{-H}^0
z\,dz\,dx=\frac{\rho g}{\lambda}\int_0^\lambda\!\!\int_0^\eta z\,dz\,dx$, which is (7.40). The sketch below (our version of
Fig. 7.6 `N172`) shows the column swap behind it.
""", equation=EQ["7.40"], ref="7.40")
nb.figure(r"""
x = np.linspace(0, 2*np.pi, 300)                          # one wavelength, k = 1 rad/m [m]
fig, ax = plt.subplots(figsize=(6.5, 2.8))                # one panel
ax.plot(x, 0.5*np.cos(x), color=COLORS["blue"], lw=2)     # the surface η = a cos kx, a = 0.5 m
ax.axhline(0, color=COLORS["muted"], ls="--", lw=1)       # still level
xa, xb, dx = np.pi, 0.0 + 0.25, 0.25                      # column A at the trough, column B on the crest side [m]
ax.fill_between([xa - dx/2, xa + dx/2], [0.5*np.cos(xa)]*2, [0, 0], color=COLORS["rose"], alpha=0.4, hatch="//")   # column A: the missing water under the trough
ax.fill_between([xb - dx/2, xb + dx/2], [0, 0], [0.5*np.cos(xb)]*2, color=COLORS["teal"], alpha=0.4, hatch="\\\\")   # column B: the lifted water on the crest
ax.annotate("", xy=(xb + 0.2, 0.3), xytext=(xa - 0.2, -0.3), arrowprops=dict(arrowstyle="->", color=COLORS["ink"], connectionstyle="arc3,rad=-0.3"))
ax.text(xa + 0.2, -0.3, "A: water removed", fontsize=8, color=COLORS["rose"])   # label A
ax.text(1.5, 0.42, "B: lifted by η onto the crest", fontsize=8, color=COLORS["teal"])   # label B
ax.text(xa - 0.08, 0.04, "dx", fontsize=8)   # the column width
ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.set_ylim(-0.7, 0.7)   # axes
ax.set_title("Building the surface: lift the column under a trough onto a crest", fontsize=10)   # the message
savefig(fig, "ch07", "c06_column_swap"); plt.show()   # save, then draw
""", see=r"""The column A of height |η| under a trough (rose) is removed and put on the crest as column B (teal).""",
    read=r"""Lifting a column of mass ρη dx (per unit width) through a height η costs $\rho g\eta^2dx$; over half a wavelength
the swaps build a whole wavelength of surface — the $\frac{\rho g}{2\lambda}\int\eta^2dx$ of (7.40).""",
    change=r"""…the fluid above were not air but a lighter liquid ρ₁: each swap also moves ρ₁ down — the work uses ρ₂ − ρ₁
(note N105, C13).""")
D("D13", ref="7.42")
note("N34", "Equipartition", r"""
$E_p=\tfrac12\rho g\overline{\eta^2}$ *(7.41)* equals $E_k=\tfrac12\rho g\overline{\eta^2}$ *(7.39)*: small oscillations of a
conservative system share their energy equally between motion and height. It fails once Coriolis forces act (Ch. 13:
geostrophic adjustment leaves potential energy locked in).
""")
remind([
    ("power of a force", "a force does work at the rate F·u; per unit area of a surface, a pressure does p·u·n (Ch. 4 P127)."),
])
note("N35", "The energy flux across a vertical line", r"""
is the pressure work: $F=\big\langle\int_{-H}^0pu\,dz\big\rangle=\big\langle\int_{-H}^0p'u\,dz\big\rangle-\rho g\langle
u\rangle\int_{-H}^0z\,dz$, which is (7.43) because ⟨u⟩ = 0 at every depth.
""", equation=EQ["7.43"], ref="7.43")
D("D14", ref="7.44")
note("N36", "Energy times a speed that is not c", r"""
$F=\big[\tfrac12\rho ga^2\big]\big[\frac c2\big(1+\frac{2kH}{\sinh2kH}\big)\big]$ *(7.44)*. The book calls the second factor
the group speed and returns to it in §7.5 (C09: $c_g=\frac c2\big[1+\frac{2kH}{\sinh(2kH)}\big]$ *(7.69)*). Where the book
says "u from (7.28)" before (7.44) it means (7.27): (7.28) is the dispersion relation.
""")
nb.worked_example("energy and power of a metre of swell", r"""
1. a = 1 m: $E=\tfrac12\rho ga^2=\tfrac12\times1000\times9.81\times1=4905$ J/m², half of it kinetic.
2. T = 10 s deep water: c = 15.6 m/s, and the second factor of (7.44) is c/2 = 7.81 m/s.
3. $F=4905\times7.81=38\,300$ W per metre of crest ≈ 38 kW/m.
4. If the energy moved at c we would have claimed 76.6 kW/m — twice too much.
""")
remind([
    ("scipy.integrate.quad and dblquad", "adaptive integration of a function over an interval or an area (inner limits may be functions) (Ch. 3 P87)."),
    ("midpoint double sums", "an area integral as Σ f × cell area at cell centres, error ∝ h² (Ch. 2 P83)."),
])
nb.code(r"""
for kH in (0.5, 1.0, 3.0):                                # k = 1 rad/m, so H = kH [m]; a = 1 m
    c_ = ch07.wave_energy(1.0, 1.0, kH, method="closed")  # (7.39), (7.41), (7.42)
    q_ = ch07.wave_energy(1.0, 1.0, kH, method="quad")    # dblquad of the integrals (7.38) and (7.40)
    print(f"kH = {kH}: closed {c_}, quad Ek = {q_['Ek']:.6f}, Ep = {q_['Ep']:.6f} J/m²")
k = (2*np.pi/10)**2/G                                     # a 10 s deep-water swell [rad/m]
print(f"F = {ch07.energy_flux(1.0, k):.1f} W/m (7.44); by averaging p′u over depth and period: "
      f"{ch07.energy_flux(1.0, k, method='quad'):.1f} W/m")
""", explain=r"""
1. `wave_energy(method="closed")` returns $E_k=E_p=2452.5$ J/m² and E = 4905 J/m² for a = 1 m, whatever kH is;
   `method="quad"` integrates (7.38) and (7.40) with `dblquad` over one wavelength and the depth — the same numbers to
   about $10^{-12}$.
2. `energy_flux` gives (7.44) = 38 291 W/m for the 10 s swell; the `quad` version averages p′u (from `wave_fields`) over
   the depth and one period: the same number.
""")
nb.check_agree(r"""
k, H, a, rho = 1.0, 1.0, 1.0, 1000.0                     # kH = 1 wave [rad/m, m, m, kg/m³]
nx, nz = 400, 200                                        # midpoint grid over one wavelength and the depth
xm = (np.arange(nx) + 0.5)*(2*np.pi/k)/nx                # cell centres in x [m]
zm = -H + (np.arange(nz) + 0.5)*H/nz                     # cell centres in z [m]
XX, ZZ = np.meshgrid(xm, zm)                             # every cell centre
f = ch07.wave_fields(XX, ZZ, 0.0, a, k, H)               # u, w at t = 0 [m/s]
Ek_mine = np.sum(0.5*rho*(f["u"]**2 + f["w"]**2))*(2*np.pi/k/nx)*(H/nz)/(2*np.pi/k)   # Σ ½ρ|u|² dx dz / λ [J/m²]
print(f"midpoint sum: E_k = {Ek_mine:.3f} J/m²")
assert np.allclose(Ek_mine, ch07.wave_energy(a, k, H)["Ek"], rtol=1e-4)   # same number as ¼ρga² = 2452.5 J/m²
""")
nb.md(r"**What does this show?** The midpoint sum of ½ρ(u² + w²) over one wavelength and the whole depth, the integral " + E("7.38") + r", gives 2452.49 J/m²; the closed form " + E("7.39") + r" = ¼ρga² = ¼ × 1000 × 9.81 × 1² = 2452.5 J/m². They agree to about 4 parts in a million, so the integral really does collapse to a quarter of ρga².")
nb.figure(r"""
k, H, a, rho = 1.0, 1.0, 1.0, 1000.0                      # kH = 1 [rad/m, m, m, kg/m³]
xs = np.linspace(0, 2*np.pi, 200)                         # one wavelength [m]
zs = np.linspace(-H, 0, 201)                              # depth grid for the column integral [m]
XX, ZZ = np.meshgrid(xs, zs)                              # (z, x) grid
f = ch07.wave_fields(XX, ZZ, 0.0, a, k, H)                # u, w at t = 0
ke_col = np.trapezoid(0.5*rho*(f["u"]**2 + f["w"]**2), zs, axis=0)   # ∫½ρ|u|² dz at each x [J/m²]
pe_col = 0.5*rho*G*(a*np.cos(k*xs))**2                    # ½ρgη² at each x [J/m²]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4))      # energy along a wavelength | c_g/c
a1.plot(xs, ke_col, color=COLORS["teal"], lw=2, label="kinetic, ∫½ρ|u|²dz")
a1.plot(xs, pe_col, color=COLORS["blue"], lw=2, label="potential, ½ρgη²")
a1.axhline(ke_col.mean(), color=COLORS["teal"], ls="--", lw=1); a1.axhline(pe_col.mean(), color=COLORS["blue"], ls=":", lw=1)
a1.set_xlabel("x [m]"); a1.set_ylabel("energy per area [J/m²]"); a1.legend(fontsize=7); a1.set_title("equal averages (dashed)", fontsize=9)
kHs = np.logspace(-1.5, 1.2, 200)                         # depth ratios [–]
a2.semilogx(kHs, 0.5*(1 + 2*kHs/np.sinh(2*kHs)), color=COLORS["accent"], lw=2)   # F/(Ec) = c_g/c from (7.44)
for Hc, lab in ((100.0, "swell, H = 100 m"), (30.0, "30 m shelf")):   # the 10 s swell of the tiny examples on two depths [m]
    kv = ch07.wavenumber_from_omega(2*np.pi/10, Hc)*Hc     # its kH on that depth (solves ω² = gk tanh kH) [–]
    ratio = 0.5*(1 + 2*kv/np.sinh(2*kv))                  # F/(Ec) from (7.44) [–]
    print(f"{lab}: kH = {kv:.2f}, F/(Ec) = {ratio:.2f}")  # the numbers the dots mark
    a2.plot(kv, ratio, "o", color=COLORS["orange"]); a2.text(kv*0.6, ratio + 0.04, lab, fontsize=7)
a2.set_ylim(0.45, 1.05); a2.set_xlabel("kH [–]"); a2.set_ylabel("F/(Ec) [–]"); a2.set_title("the energy speed ÷ crest speed", fontsize=9)
fig.suptitle("Energy is shared equally and travels slower than the crests", fontweight="bold")
savefig(fig, "ch07", "c06_energy"); plt.show()
""", see=r"""Left: the depth-integrated kinetic energy (teal) and ½ρgη² (blue) along one wavelength — different shapes, the
same average. Right: the ratio F/(Ec) falls from 1 (shallow) to ½ (deep).""",
    read=r"""In deep water energy arrives at half the crest speed (C09 explains why); the 10 s swell on the 30 m shelf
(kH = 1.37, printed above the figure) carries it at 0.68 of the crest speed.""",
    change=r"""…the wave were in 2 m of water: kH ≈ 0.2 and the ratio ≈ 0.99 — energy and crests travel together.""")
confusion(r"""
surface-wave E is energy per square metre of sea surface (the whole water column below it), and F is power per metre of
crest. For internal waves in C16 the same letters mean per cubic metre and per square metre.
""")
whatif(r"""
…the Earth's rotation joined in (long waves, Ch. 13)? Kinetic and potential energy stop being equal: geostrophic
adjustment keeps some potential energy in a balanced current.
""")
nb.pointer("**Refraction** (the book's Figs. 7.8–7.9: waves turning toward a beach and wrapping round an island) belongs to "
           "this section of the book; we teach it in C10 (§7.5), where the rays that explain it are derived (notes N47, "
           "N48, N174, N175).")

# =====================================================================================================================
# A.3 §7.3 — R09, C07 (N49–N57, N176, D15, D16, P169, P170, E3)
# =====================================================================================================================
nb.section("7.3", "Influence of Surface Tension", intro=r"""
**What is this section about?** A curved surface pulls itself flat, like a stretched membrane — a second restoring force
besides gravity. It only matters for short waves, but there it wins, and together the two forces leave a slowest possible
wave on water.
""")
nb.recap("R09", "The Laplace pressure jump", r"""
Across a curved interface the pressure jumps by $\Delta p=\sigma\Big(\frac1{R_1}+\frac1{R_2}\Big)$ *(1.5)*, higher on the
side of the centres of curvature; Ch. 4 re-derived it as a force balance on a small patch (`interfaces.laplace_jump_from_balance`).
Below: a surface curved with a 1 cm radius in one direction and straight in the other (a 10⁹ m radius), as under a
straight crest.
""", where="Ch. 1 §1.6, Ch. 4 §4.10")
nb.code(r"""
print(interfaces.laplace_jump_from_balance(0.0727, 0.01, 1e9), "Pa")   # σ = 0.0727 N/m, R₁ = 1 cm, R₂ ≈ ∞ → σ/R₁ = 7.27 Pa
""")
core("C07", r"Capillary–gravity waves and the minimum speed", r"""
Why do the smallest waves on a pond run faster than slightly longer ones — and is there a slowest wave?
""", eqs=("7.57", "7.58"))
note("N49", "A second restoring force", r"""
Surface tension acts like a stretched membrane. Pure capillary waves (no gravity) are rare on Earth, so we extend §7.2
rather than replace it.
""")
problem(r"""
Rain on a lake makes tiny ripples that race ahead of the wind waves; a stone makes a ring of ripples with a calm patch in
the middle. Short waves feel the surface's tension, long ones gravity. We want the speed with both forces and the
wavelength where the speed is smallest.
""")
nb.md(r"""
#### The idea

| restoring force | stronger for | term in c² (deep water) | branch |
|---|---|---|---|
| gravity (weight of the hump) | long waves | g/k | $c=\sqrt{g\lambda/2\pi}$ grows with λ |
| surface tension (curvature) | short waves | σk/ρ | $c=\sqrt{2\pi\sigma/\rho\lambda}$ grows as λ shrinks |

**One force favours long waves, the other short ones — so the speed has a minimum in between.**
""")
P("P169", "curvature of a plane curve", r"""
For a curve z = η(x) the curvature (1/radius of the circle that fits the curve best there) is
$1/R=\eta_{xx}/(1+\eta_x^2)^{3/2}$; for gentle slopes $1/R\approx\eta_{xx}$. Sign: under a crest $\eta_{xx}<0$ and the
centre of curvature lies below, inside the water.
""", code=r"""
import numpy as np                                # numbers
a, k = 0.001, 2*np.pi/0.01                        # a 1 mm ripple with λ = 1 cm [m, rad/m]
etaxx_crest = -a*k**2                             # η = a cos kx at x = 0: η'' = −ak² [1/m]
print(etaxx_crest, 1/abs(etaxx_crest))            # −394.8 1/m: radius 2.5 mm at the crest
""")
D("D15", ref="7.56")
note("N50–N53", "Where D15 passes", r"""
D15 passes through the curvature condition
$p_a-(p)_{z=\eta}=\sigma\frac1R=\sigma\frac{\partial^2\eta/\partial x^2}{[1+(\partial\eta/\partial x)^2]^{3/2}}\cong\sigma
\frac{\partial^2\eta}{\partial x^2}$ *(7.53)*, its gauge form $(p)_{z=\eta}=-\sigma\frac{\partial^2\eta}{\partial x^2}$
*(7.54)*, the new dynamic condition $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=\frac\sigma\rho\frac{\partial^2\eta}
{\partial x^2}-g\eta$ *(7.55)* and ends at $\omega=\sqrt{k\big(g+\frac{\sigma k^2}{\rho}\big)\tanh kH}$ *(7.56)*: only the
dispersion relation changes, φ keeps its form (7.26).
""")
remind([
    ("slope zero at a minimum", "a smooth function has zero slope at its minimum and a positive second derivative there (Ch. 1 P19)."),
])
D("D16", ref="7.58")
note("N54", "The minimum", r"""
For λ < λ_m surface tension dominates.
""", equation=EQ["7.58"], ref="7.58")
nb.worked_example("clean water at 20 °C", r"""
1. σ = 72.74 mN/m (IAPWS, Ch. 1), ρ = 998.2 kg/m³, g = 9.81 m/s².
2. $\lambda_m=2\pi\sqrt{0.07274/(998.2\times9.81)}=2\pi\times2.73\times10^{-3}=0.0171$ m = 1.71 cm.
3. $c_{min}=(4\times9.81\times0.07274/998.2)^{1/4}=(2.86\times10^{-3})^{1/4}=0.231$ m/s = 23.1 cm/s.
4. Under the crest of a 1 mm ripple with λ = 1 cm: $k=2\pi/0.01=628.3$ rad/m, $k^2=3.95\times10^5$ m⁻², so
   $p=-\sigma\eta_{xx}=\sigma ak^2=0.07274\times0.001\times3.95\times10^5=+28.7$ Pa, pushing the crest down.
""")
note("N55–N57", "Air–water numbers and pure capillary waves", r"""
With our numbers above, (7.59) reads $c_{min}\approx23.1$ cm/s at $\lambda_m\approx1.71$ cm (the book's printed values
are kept in the private test file). `N56` Only ripples shorter than about 7 cm feel σ noticeably, and below about 4 mm
gravity hardly matters; soap and oil films lower σ and change the ripples. `N57` Dropping g in deep water gives pure
capillary waves — shorter is faster:
""", equation=EQ["7.60"], ref="7.60")
P("P170", "scipy.optimize.minimize_scalar", r"""
Finds the minimum of a function of one variable inside a bracket (`method="bounded"`); we use it to check a formula for a
minimum by brute force.
""", code=r"""
from scipy.optimize import minimize_scalar        # 1-D minimiser
r = minimize_scalar(lambda x: (x - 2)**2 + 1, bounds=(0, 5), method="bounded")   # a parabola with its minimum at x = 2
print(r.x, r.fun)                                 # 2.0 1.0
""")
nb.code(r"""
sig, rho = ch01.surface_tension_water(293.15), 998.2     # Ch. 1: σ of clean water at 20 °C [N/m]; density [kg/m³]
print(f"σ = {sig:.5f} N/m;", ch07.capillary_minimum(sig, rho))   # (7.58): c_min [m/s], λ_m [m], k_m [rad/m]
for lam in (0.004, 0.0171, 0.07, 1.0):                   # four wavelengths [m]
    k = 2*np.pi/lam                                      # wavenumber [rad/m]
    print(f"λ = {lam*100:5.2f} cm: c = {ch07.phase_speed(k, np.inf, sigma=sig, rho=rho):.3f} m/s with σ, "
          f"{ch07.phase_speed(k, np.inf):.3f} m/s without")   # (7.57) vs (7.45)
print(f"crest pressure of a 1 mm, 1 cm ripple: {ch07.capillary_surface_pressure(-0.001*(2*np.pi/0.01)**2, sig):.1f} Pa;",
      "curvature (linear form):", ch07.curvature(0.0, -394.8, linear=True), "1/m")   # (7.54), (7.53)
""", explain=r"""
1. `surface_tension_water` (the IAPWS fit of Ch. 1) gives σ = 0.07274 N/m; `capillary_minimum` returns
   $c_{min}=0.2312$ m/s at $\lambda_m=1.712$ cm.
2. At λ = 4 mm tension wins (0.347 m/s with σ, only 0.079 without); at 1.71 cm the speed is the minimum; at 7 cm the two
   differ by 3 %; at 1 m surface tension is invisible (1.250 m/s both ways).
3. The liquid pressure under the crest of a 1 mm ripple, $(p)_{z=\eta}=-\sigma\frac{\partial^2\eta}{\partial x^2}$ *(7.54)*,
   is +28.7 Pa — the crest is pushed down, like gravity.
""")
nb.check_agree(r"""
from scipy.optimize import minimize_scalar                 # the 1-D minimiser of primer P170
lam = np.logspace(-3.5, 0, 4001)                           # a log grid of wavelengths [m]
c = ch07.phase_speed(2*np.pi/lam, np.inf, sigma=sig, rho=rho)   # (7.57) in deep water [m/s]
i = np.argmin(c)                                           # the smallest sampled speed
r = minimize_scalar(lambda kk: G/kk + sig*kk/rho, bounds=(2*np.pi/lam[i + 1], 2*np.pi/lam[i - 1]), method="bounded",
                    options=dict(xatol=1e-10))             # refine: minimise c² = g/k + σk/ρ in k
cmin_mine, lam_mine = np.sqrt(r.fun), 2*np.pi/r.x          # [m/s], [m]
d = ch07.capillary_minimum(sig, rho)                       # the closed form (7.58)
print(f"brute force: c_min = {cmin_mine:.6f} m/s at λ = {lam_mine*100:.5f} cm")
assert np.allclose([cmin_mine, lam_mine], [d["c_min"], d["lam_m"]], rtol=1e-6)   # same numbers as (7.58)
""")
nb.figure(r"""
lam = np.logspace(-3, 3, 400)                              # wavelengths 1 mm … 1 km [m]
k = 2*np.pi/lam                                            # wavenumbers [rad/m]
sg, rh, H = 0.0727, 1000.0, 1.0                            # σ [N/m], ρ [kg/m³], a 1 m deep tank [m]
fig, ax = plt.subplots(figsize=(6.5, 4))                   # one log–log panel
ax.loglog(lam, ch07.phase_speed(k, H, sigma=sg, rho=rh), color=COLORS["ink"], lw=2, label="(7.57), H = 1 m")   # c with both forces, (7.57)
ax.loglog(lam, np.sqrt(2*np.pi*sg/(rh*lam)), "--", color=COLORS["rose"], label="capillary branch (7.60)")   # tension only
ax.loglog(lam, np.sqrt(G*lam/(2*np.pi)), "--", color=COLORS["blue"], label="deep gravity branch (7.45)")   # gravity only, deep
ax.axhline(np.sqrt(G*H), ls=":", color=COLORS["muted"], label="shallow plateau √(gH)")   # very long waves on 1 m
m = ch07.capillary_minimum(sg, rh, G)                      # (7.58)
ax.plot(m["lam_m"], m["c_min"], "o", color=COLORS["orange"], ms=7)   # the minimum (7.58)
ax.text(m["lam_m"]*1.3, m["c_min"]*0.8, f"c_min = {m['c_min']*100:.1f} cm/s at λ_m = {m['lam_m']*100:.2f} cm", fontsize=8, color=COLORS["orange"])
ax.axvspan(1e-3, 0.004, color=COLORS["rose"], alpha=0.08); ax.text(1.1e-3, 3.2, "capillary\n< 4 mm", fontsize=7, color=COLORS["rose"])
ax.axvspan(0.004, 0.07, color=COLORS["amber"], alpha=0.06); ax.text(8e-3, 3.2, "σ matters\n< 7 cm", fontsize=7, color=COLORS["amber"])
ax.set_ylim(0.1, 5); ax.set_xlabel("wavelength λ [m]"); ax.set_ylabel("c [m/s]")   # axes with units
ax.set_title("Two restoring forces leave a slowest wave", fontsize=10); ax.legend(fontsize=7, loc="lower right")   # message and legend
savefig(fig, "ch07", "c07_capillary"); plt.show()   # save, then draw
""", see=r"""A valley at 1.7 cm between the rising capillary branch (rose) and the gravity branch (blue), then the $\sqrt{gH}$
plateau for very long waves (our remake of Fig. 7.10 `N176`).""",
    read=r"""At any λ the square of the black curve is the sum of the squares of the two dashed branches, times tanh kH.""",
    change=r"""…soap halved σ: the valley moves to 1.2 cm and 19 cm/s ($\lambda_m\propto\sigma^{1/2}$, $c_{min}\propto\sigma^{1/4}$).""")
nb.plotly(r"""
lam = np.logspace(-3.5, 0.5, 300)                          # wavelengths [m]
def curves(sg):                                            # c(λ) and the minimum for one surface tension σ [N/m]
    c = ch07.phase_speed(2*np.pi/lam, np.inf, sigma=sg, rho=1000.0)          # (7.57), deep water
    m = ch07.capillary_minimum(sg, 1000.0, G) if sg > 0 else {"lam_m": np.nan, "c_min": np.nan}   # (7.58)
    return {"c(λ) (7.57)": (lam, c), "minimum (7.58)": ([m["lam_m"]], [m["c_min"]])}
fig = slider_figure(curves, "σ", np.linspace(0.0, 0.5, 26 if not FAST else 11), unit="N/m", xlabel="λ [m]",
                    ylabel="c [m/s]", title="The minimum slides with σ", modes={"minimum (7.58)": "markers"})
logaxes(fig, x=(10**-3.5, 10**0.5), y=(0.05, 3))           # log axes, fixed
recolor(fig, {"c(λ) (7.57)": COLORS["ink"], "minimum (7.58)": COLORS["orange"]})
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** At σ = 0 there is no minimum (pure gravity waves get ever slower as λ shrinks); raising σ lifts the short-wave branch and
pushes the minimum to longer, faster waves.

**What would change if…** the liquid were denser at the same σ (mercury-like): σ/ρ is what matters, so the minimum would move back
toward shorter, slower waves.
""")
explainer("capillary_gravity_waves", "Why is there a slowest ripple?",
          r"""Sweep λ over four decades and change the liquid: the gravity and surface-tension term bars trade places, the
minimum moves, and the status names the branch. Only interaction shows the two forces competing.""", "",
          ["Press *λ = λ_m*: the two bars are equal.",
           "Switch the liquid to mercury: the minimum moves to about 1.2 cm and 19 cm/s.",
           "Drag λ to 5 mm and read the crest pressure in the Explain tab.",
           "Step the Derivation (D16) to 'set the slope to zero'."])
whatif(r"""
…the group speed were asked for instead of the phase speed? It also has a minimum — lower than $c_{min}$ and at a longer
wavelength; C09 computes it with `min_group_velocity` (0.1776 m/s at λ ≈ 4.35 cm for clean water) — which is why the ring
of ripples from a stone leaves a calm centre (C09, note N72).
""")

# =====================================================================================================================
# A.4 §7.4 — C08 (N58–N63, N177, N178, D17, D18, P171, A3, E4)
# =====================================================================================================================
nb.section("7.4", "Standing Waves", intro=r"""
**What is this section about?** Two equal waves travelling in opposite directions add into a pattern that does not
travel: fixed nodes, sloshing in between. Walls allow only the patterns whose nodes fit — so a lake, a harbour or a bathtub
rings at its own set of periods (seiches).
""")
core("C08", r"Standing waves and seiches", r"""
Which waves can live in a closed basin, and with which periods?
""", eqs=("7.64", "7.65"))
problem(r"""
Push the water in a bathtub once and it sloshes back and forth at a fixed rhythm. Lake Geneva does the same after a storm,
with a period of about an hour; harbours do it and ships break their moorings. We want to know which sloshing patterns a
basin allows and their periods.
""")
nb.md(r"""
#### The idea

```
right-going  a cos(kx − ωt)  +  left-going  a cos(kx + ωt)  =  2a cos kx · cos ωt
the shape cos kx never moves; it only breathes in time (cos ωt)
nodes of η at kx = π/2, 3π/2, …   — walls must sit where u = 0 (antinodes of η)
```
""")
P("P171", "sum-to-product identities", r"""
cos A + cos B = 2 cos((A − B)/2) cos((A + B)/2) and cos A − cos B = −2 sin((A + B)/2) sin((A − B)/2): a sum of two waves
becomes a product of a slow and a fast factor. They turn two travelling waves into a standing wave here and into beats in
C09.
""", code=r"""
import numpy as np                                # numbers
A, B = 1.1, 0.3                                   # any two angles [rad]
print(np.cos(A) + np.cos(B), 2*np.cos((A - B)/2)*np.cos((A + B)/2))   # both 1.4089
""")
note("N58", "The left-going wave", r"""
It is $\eta(x,t)=a\cos[kx-\omega t]$ *(7.2)* with ω → −ω: its crests move to −x.
""", equation=EQ["7.61"], ref="7.61")
nb.code(r"""
print(ch07.sinusoid(0.0, np.array([0.0, 1.0]), 1.0, 1.0, 1.0, direction=-1))   # (7.61) at x = 0, t = 0 and 1 s
print(ch07.sinusoid(np.array([0.0, -1.0]), 1.0, 1.0, 1.0, 1.0, direction=-1))  # at t = 1 s the crest is at x = −1 m
""", explain=r"""
With `direction=-1` the sinusoid is $a\cos(kx+\omega t)$: at x = 0 the surface falls from 1 to cos 1 = 0.540 m in 1 s,
and at t = 1 s the value 1 (the crest) sits at x = −1 m — the wave has moved one metre to the left.
""")
D("D17", ref="7.63")
note("N59, N61", "The standing wave's stream function and velocity", r"""
$\psi=\frac{a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}[\cos(kx-\omega t)-\cos(kx+\omega t)]=\frac{2a\omega}{k}\frac{\sinh
k(z+H)}{\sinh kH}\sin kx\sin\omega t$ *(7.62)* and $u=2a\omega\frac{\cosh k(z+H)}{\sinh kH}\sin kx\sin\omega t$ *(7.63)*:
where η has a node, u is largest.
""")
nb.figure(r"""
k, H, a = 1.0, 1.0, 0.1                                    # kH = 1; amplitude 10 cm [rad/m, m, m]
w = ch07.omega_gravity(k, H)                               # ω [rad/s]
xg, zg = np.linspace(0, 2*np.pi, 161 if not FAST else 81), np.linspace(-H, 0, 81 if not FAST else 41)
XX, ZZ = np.meshgrid(xg, zg)                               # grid over one wavelength and the depth
f = ch07.standing_wave_fields(XX, ZZ, (np.pi/2)/w, a, k, H)   # ωt = π/2: largest flow, η = 0
fig, ax = plt.subplots(figsize=(7, 3.2))                   # one panel
lv = np.linspace(-np.abs(f["psi"]).max(), np.abs(f["psi"]).max(), 13)
ax.contour(XX, ZZ, f["psi"], levels=lv, colors=[COLORS["teal"]], linewidths=1)   # streamlines of (7.62)
ax.contour(XX, ZZ, f["psi"], levels=[0.0], colors=[COLORS["ink"]], linewidths=2.5)   # ψ = 0
eta_pi = ch07.standing_wave_fields(xg, 0.0, np.pi/w, a, k, H)["eta"]   # the surface at ωt = π (drawn ×1)
ax.plot(xg, eta_pi, color=COLORS["blue"], lw=2, label="surface at ωt = π")
ax.set_ylim(-H, 0.25); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.legend(fontsize=7, loc="lower left")
ax.set_title("Standing wave: streamlines at ωt = π/2 (kH = 1)", fontsize=10)
savefig(fig, "ch07", "c08_standing_psi"); plt.show()
""", see=r"""Streamlines leaving the region under one crest and entering the region under the next trough, with vertical
ψ = 0 lines at the crests and troughs (our remake of Fig. 7.11 `N177`).""",
    read=r"""The vertical ψ = 0 lines at kx = 0, π, 2π are where a wall can stand: u = 0 there at all times.""",
    change=r"""…the snapshot were at ωt = 0: η is largest and the flow is zero everywhere — the energy is all potential at
that instant.""")
remind([
    ("seiche and the eigenvalue idea", "in A·b = λb only special λ allow a non-zero b (Ch. 2 P80); here the walls allow only special wavelengths — the basin's own modes."),
])
note("N60", "A seiche", r"""
is a standing oscillation of a closed basin (length L, depth H, vertical walls): only patterns with u = 0 at both walls
survive — an eigenvalue problem (the idea of Ch. 2 P80: the boundary picks discrete solutions).
""")
D("D18", ref="7.65")
note("N62, N63", "Allowed wavelengths and the lake's frequencies", r"""
Walls at x = 0 and L allow $kL=(n+1)\pi,\ \lambda=\frac{2L}{n+1}$ *(7.64)* — the longest is 2L, then L, 2L/3 … — and the
natural frequencies (7.65). A rectangular basin L × b rings with $k^2=(m\pi/L)^2+(n\pi/b)^2$ (
`ch07.basin_modes`).
""", equation=EQ["7.65"], ref="7.65")
nb.worked_example("a bathtub and a lake", r"""
1. Bathtub L = 1.5 m, H = 0.2 m, n = 0: $k=\pi/L=2.094$ rad/m, kH = 0.419, tanh = 0.396.
2. $\omega^2=9.81\times2.094\times0.396=8.14$ s⁻² → ω = 2.85 rad/s, T = 2.20 s.
3. Shallow estimate $T=2L/\sqrt{gH}=3/1.40=2.14$ s (3 % short).
4. A lake L = 50 km, H = 100 m: kH = 0.0063 — shallow; $T_0=2L/\sqrt{gH}=100\,000/31.3=3193$ s ≈ 53 min.
""")
nb.code(r"""
print(ch07.seiche_modes(1.5, 0.2, 0))                      # bathtub, gravest mode: k, λ, ω, T
print(ch07.seiche_modes(50e3, 100.0, 0))                   # a 50 km lake, 100 m deep
print("lake periods [min]:", [round(ch07.seiche_modes(50e3, 100.0, n)["T"]/60, 1) for n in range(4)])   # n = 0…3
print(ch07.basin_modes(1.5, 0.7, 0.2, 1, 1))               # rectangular 1.5 m × 0.7 m tub, mode (1, 1)
x = np.linspace(0, 2, 5)                                   # five points [m]
s = ch07.standing_wave_fields(x, -0.1, 0.4, 0.01, np.pi, 0.2)                   # the standing wave (D17)
r = ch07.wave_fields(x, -0.1, 0.4, 0.01, np.pi, 0.2)                            # right-going wave (7.2)
l = ch07.wave_fields(x, -0.1, 0.4, 0.01, np.pi, 0.2, direction=-1)              # left-going wave (7.61)
assert np.allclose(s["eta"], r["eta"] + l["eta"]) and np.allclose(s["u"], r["u"] + l["u"])   # sum of the two
print("standing wave = sum of the two travelling waves ✓")
""", explain=r"""
1. The bathtub rings at T = 2.203 s; the lake at 3193 s = 53.2 min, then 26.6, 17.7 and 13.3 min — in a shallow basin
   the period falls like 1/(n + 1).
2. `basin_modes` adds a cross-basin wavenumber: $k^2=(m\pi/L)^2+(n\pi/b)^2$.
3. The assertion: the standing wave is exactly the sum of the two travelling waves of C03, both in η and in u.
""")
nb.check_agree(r"""
L, H = 1.5, 0.2                                            # the bathtub [m]
Ts = [2*np.pi/np.sqrt(G*(n + 1)*np.pi/L*np.tanh((n + 1)*np.pi*H/L)) for n in range(5)]   # (7.65) by hand [s]
print(np.round(Ts, 4))
assert np.allclose(Ts, [ch07.seiche_modes(L, H, n)["T"] for n in range(5)], rtol=1e-12)   # same periods
""")
nb.md(r"**What does this show?** The first five periods of a 1.5 m tub holding 0.2 m of water: 2.20 s for the gravest mode (n = 0), then 1.18, 0.87, 0.72 and 0.63 s. Our line of " + E("7.65") + r" (with T = 2π/ω) and `seiche_modes` agree to $10^{-12}$. The periods fall more slowly than 1/(n + 1): the shorter modes feel the bottom less and run slower than $\sqrt{gH}$.")
nb.figure(r"""
L, H = 1.0, 0.1                                            # a unit basin [m]
x = np.linspace(0, L, 300)                                 # positions [m]
fig, axs = plt.subplots(2, 1, figsize=(7, 3.8), sharex=True)   # η above, u below
for n, cl in ((0, COLORS["blue"]), (1, COLORS["accent"])):  # the two gravest modes
    k = (n + 1)*np.pi/L                                    # (7.64)
    axs[0].plot(x, np.cos(k*x), color=cl, lw=1, label=f"n = {n}")          # η ∝ cos kx (at ωt = 0)
    axs[1].plot(x, np.sin(k*x), color=cl, lw=2, label=f"n = {n}")          # u ∝ sin kx (at ωt = π/2), (7.63)
    for xn in x[np.abs(np.cos(k*x)) < 0.012]:              # nodes of η
        axs[0].plot(xn, 0, "|", color=cl, ms=10)
for axp in axs:
    axp.axvline(0, color=COLORS["ink"], lw=3); axp.axvline(L, color=COLORS["ink"], lw=3); axp.axhline(0, color=COLORS["muted"], lw=0.5)
axs[0].set_ylabel("η (scaled)"); axs[1].set_ylabel("u (scaled)"); axs[1].set_xlabel("x/L [–]"); axs[0].legend(fontsize=7, ncol=2)
axs[0].set_title("Seiche modes: u = 0 at the walls, largest where η has its node", fontsize=10)
savefig(fig, "ch07", "c08_seiche_modes"); plt.show()
""", see=r"""Top, the surface shapes of modes n = 0 and 1 with their nodes marked; bottom, the horizontal velocity: zero at
both walls, one bulge for n = 0, two opposite bulges for n = 1 (our remake of Fig. 7.12 `N178`).""",
    read=r"""Where u is largest η has its node: water rushes through the node from one side of the basin to the other.""",
    change=r"""…n = 2: three bulges of u and λ = 2L/3; in a shallow lake the period shortens to a third of the gravest.""")
nb.plotly(r"""
xb = np.linspace(0, 1, 200)                                # a unit basin, x/L [–]
fig = slider_figure(lambda n: {"η (7.64)": (xb, np.cos((n + 1)*np.pi*xb)),        # surface shape of mode n
                               "u (7.63)": (xb, np.sin((n + 1)*np.pi*xb))},       # velocity shape of mode n
                    "n", range(6), xlabel="x/L [–]", ylabel="shape [–]", title="Seiche mode n: kL = (n + 1)π",   # slider over the mode number n = 0…5
                    yrange=(-1.1, 1.1))   # fixed vertical range
recolor(fig, {"η (7.64)": COLORS["blue"], "u (7.63)": COLORS["accent"]})   # surface blue, velocity purple
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** Step n from 0 to 5: the surface always has antinodes at the walls, the velocity always vanishes there, and one more half
wavelength fits each time.

**What would change if…** one end of the basin were open to the sea (a bay): the open end would be a node of η instead of an antinode,
and a quarter wavelength (not a half) would fit in the longest mode.
""")
nb.animation(r"""
k, H, a = 1.0, 1.0, 0.12                                    # kH = 1, amplitude 12 cm (exaggerated) [rad/m, m, m]
L = 2*np.pi                                                 # basin length for mode n = 1: λ = L [m]
w = ch07.omega_gravity(k, H)                                # ω [rad/s]
ng = 121 if not FAST else 81                                # contour grid size
xg, zg = np.linspace(0, L, ng), np.linspace(-H, 0, ng//2)   # grid [m]
XX, ZZ = np.meshgrid(xg, zg)
fig, ax = plt.subplots(figsize=(7, 3.2))
ax.axvline(0, color=COLORS["ink"], lw=3); ax.axvline(L, color=COLORS["ink"], lw=3)   # the walls
(ghost_r,) = ax.plot(xg, 0*xg, color=COLORS["orange"], lw=1, alpha=0.5)   # right-going component (ghost)
(ghost_l,) = ax.plot(xg, 0*xg, color=COLORS["orange"], lw=1, alpha=0.5, ls="--")   # left-going component
(surf,) = ax.plot(xg, 0*xg, color=COLORS["blue"], lw=2)    # their sum: the standing wave
ax.set_xlim(0, L); ax.set_ylim(-H, 0.3); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")
held = []                                                   # the current streamline set (redrawn every frame)
nf = 40 if not FAST else 20                                 # frames over one period


def update(i):   # draw frame i
    t = i/nf*2*np.pi/w                                      # time [s]
    ghost_r.set_ydata(a*np.cos(k*xg - w*t)); ghost_l.set_ydata(a*np.cos(k*xg + w*t))   # (7.2) and (7.61)
    surf.set_ydata(ch07.standing_wave_fields(xg, 0.0, t, a, k, H)["eta"])             # 2a cos kx cos ωt
    for c_ in held:                                          # remove last frame's streamlines
        c_.remove()
    held.clear()
    psi = ch07.standing_wave_fields(XX, ZZ, t, a, k, H)["psi"]   # (7.62)
    if np.abs(psi).max() > 1e-6:                              # skip the instants with no flow
        held.append(ax.contour(XX, ZZ, psi, levels=np.linspace(-0.6, 0.6, 13), colors=[COLORS["teal"]], linewidths=0.8))
    ax.set_title(f"standing wave in a basin, ωt = {w*t:.2f} rad", fontsize=10)
    return []   # nothing to blit


show_animation(animate(update, frames=nf, fig=fig, interval=70))
""")
see_read_change(r"""The two travelling components (orange ghosts) slide through each other; their sum (blue) never moves
sideways; the streamlines (teal) reverse every half period — the `A3` animation.""",
                r"""The walls stand where u = 0 at every instant; the water sloshes from under one crest to the other side,
passing through the node in the middle.""",
                r"""…one ghost were removed: the blue curve would travel — a single progressive wave cannot live in a closed
basin.""")
explainer("seiche_standing_waves", "Which waves can live in a lake?",
          r"""Toggle the two travelling waves and watch their sum pin its nodes; then change the basin and the mode and
compare the period with real basins in the table.""", "",
          ["Switch off the left-going wave: the pattern starts to travel.",
           "Choose the 50 km lake and step n from 0 to 3: the period drops as 1/(n + 1).",
           "Drag H up in the bathtub until the shallow estimate fails by 10 %."])
whatif(r"""
…the basin were a bay open to the sea? Its open end becomes a node of η, so the longest wave is 4L (a quarter wave) — the
harbour-resonance rule; tides in long gulfs resonate this way (Ch. 13).
""")

# =====================================================================================================================
# A.5 §7.5 — C09 (N64–N73, N179–N182, D19–D21, P172, P173, A2, E5), C10 (N47, N48, N74–N80, N174, N175, N183, N184,
#            D22–D24, P174, P175, E6)
# =====================================================================================================================
nb.section("7.5", "Group Velocity, Energy Flux, and Dispersion", intro=r"""
**What is this section about?** When the speed depends on the wavelength, a group of waves moves at a different speed
from its crests. That group speed, dω/dk, is the speed of the energy — the answer to C06's cliff-hanger. Then we let the
depth change slowly: frequency stays constant along rays while wavelength and direction change, which is why waves turn
toward beaches.
""")
core("C09", r"Group velocity and the energy flux", r"""
Why do crests appear at the back of a group of swell and vanish at its front — and which speed carries the energy?
""", eqs=("7.67", "7.71"))
note("N64", "Dispersion is common", r"""
Speeds that depend on wavelength are common for waves on interfaces between materials (Rayleigh and Stoneley waves in
solids, liquid–liquid interfaces); we stay with water.
""")
problem(r"""
Surfers wait for "sets": swell arrives in groups of a few big waves, then a lull. Watch one crest in a group and it runs
forward through the group and fades at the front, while new crests grow at the back. The group — and the energy — is
slower than the crests. We want that speed and the proof that energy rides with it.
""")
nb.md(r"""
#### The idea

```
ω(k) curve:   slope of the chord from the origin to (k, ω)   = ω/k   = c    (crests)
              slope of the tangent at (k, ω)                 = dω/dk = c_g  (groups, energy)
deep water:   ω = √(gk) bends down  ⇒  tangent flatter than chord  ⇒  c_g = c/2
```
""")
remind([
    ("derivative as the limit of a difference quotient", r"$d\omega/dk=\lim_{\Delta k\to0}\Delta\omega/\Delta k$: the chord slope becomes the tangent slope (Ch. 1 P19)."),
])
D("D19", ref="7.67")
note("N65", "Beats", r"""
Two nearby waves make the beat pattern (7.66). The book prints ½Δω **x** in (7.66); it must be ½Δω **t** — the next line
of the text and the figure use t, and with x the envelope could not move. The code compares the two forms: the node of
the correct envelope moves by $\frac{\Delta\omega}{\Delta k}\times10$ s, the printed one does not move.
""", equation=EQ["7.66"], ref="7.66")
P("P180", "np.sign and np.nonzero: finding where a curve crosses zero", r"""
`np.sign(e)` turns every value into −1, 0 or +1. Where two neighbours `e[:-1]` and `e[1:]` have different signs the curve
crossed zero between them. `np.nonzero(mask)` returns the indices where a True/False array is True (as a tuple, one array
per axis — hence the `[0]`), so `np.nonzero(...)[0][0]` is the first crossing; a straight line between the two
neighbours then places the zero precisely.""", code=r"""
import numpy as np                                        # arrays
e = np.array([2.0, 1.0, -1.0, -2.0])                      # a curve sampled at x = 0, 1, 2, 3
i = np.nonzero(np.sign(e[:-1]) != np.sign(e[1:]))[0][0]   # first index before a sign change → 1
print(i, i - e[i]*(1.0)/(e[i + 1] - e[i]))                 # 1 and the interpolated zero at x = 1.5
""")
nb.code(r"""
x = np.linspace(-40, 80, 24001)                           # positions [m]
env = lambda t, printed: ch07.beat_wave(x, t, 0.9, 1.1, printed=printed)["envelope"]   # (7.66), deep water
b = ch07.beat_wave(x, 0.0, 0.9, 1.1)                      # k₁ = 0.9, k₂ = 1.1 rad/m
def node(e):                                              # first zero of the envelope beyond x = 5 m [m]
    i = np.nonzero((x[:-1] > 5) & (np.sign(e[:-1]) != np.sign(e[1:])))[0][0]   # first sign change
    return x[i] - e[i]*(x[i + 1] - x[i])/(e[i + 1] - e[i])   # linear interpolation to the zero
for pr in (False, True):                                  # correct ½Δω t, then the printed ½Δω x
    print(f"printed={pr}: envelope node at t = 0: {node(env(0.0, pr)):.3f} m, at t = 10 s: {node(env(10.0, pr)):.3f} m")
print(f"Δω/Δk = {b['cg_finite']:.4f} m/s, dω/dk at k = 1 = {ch07.group_velocity(1.0):.4f} m/s, c = {b['c']:.4f} m/s")
""", explain=r"""
With ½Δω t the node moves from 15.7 m to 31.4 m in 10 s — 1.568 m/s, close to the tangent slope $d\omega/dk=\tfrac12\sqrt{g/k}=1.566$ m/s
at k = 1 rad/m, while the crests move at c = 3.13 m/s. With the printed ½Δω x the node stays where it is: that form
describes a frozen pattern, not a wave group.
""")
nb.figure(r"""
x = np.linspace(-10, 90, 3000)                            # positions [m]
fig, ax = plt.subplots(figsize=(7, 3))                    # one panel
b0, b1 = ch07.beat_wave(x, 0.0, 0.9, 1.1), ch07.beat_wave(x, 6.0, 0.9, 1.1)   # t = 0 and t = 6 s
ax.plot(x, b1["eta"], color=COLORS["blue"], lw=1)         # the beats at t = 6 s
ax.plot(x, np.abs(b1["envelope"]), color=COLORS["accent"], lw=2); ax.plot(x, -np.abs(b1["envelope"]), color=COLORS["accent"], lw=2)   # ± the envelope
ax.plot(x, np.abs(ch07.beat_wave(x, 6.0, 0.9, 1.1, printed=True)["envelope"]), ls="--", color=COLORS["rose"], lw=1,   # the printed slip's envelope
        label="printed ½Δω x: frozen envelope")          # the slip's ghost (it does not move)
xn0 = np.pi/b0["dk"]                                      # a node of the correct envelope at t = 0 [m]
xn1 = xn0 + b0["cg_finite"]*6.0                           # it moves at Δω/Δk [m]
ax.plot([xn0, xn1], [0, 0], "ko"); ax.annotate("", xy=(xn1, -2.3), xytext=(xn0, -2.3), arrowprops=dict(arrowstyle="->", color=COLORS["accent"]))   # the node then and now, and its arrow
ax.text(xn0, -2.8, "node: c_g t", fontsize=8, color=COLORS["accent"])   # label
xc0 = 0.0; xc1 = xc0 + b0["c"]*6.0                        # a crest moves at c [m]
ax.plot([xc1], [np.abs(b1["envelope"][np.argmin(np.abs(x - xc1))])], "o", color=COLORS["orange"])   # the tracked crest at t = 6 s
ax.annotate("", xy=(xc1, 2.4), xytext=(xc0, 2.4), arrowprops=dict(arrowstyle="->", color=COLORS["orange"]))   # how far the crest moved
ax.text(xc0, 2.55, "crest: c t", fontsize=8, color=COLORS["orange"])   # label
ax.set_ylim(-3.1, 3.1); ax.set_xlabel("x [m]"); ax.set_ylabel("η [m]"); ax.legend(fontsize=7, loc="lower right")   # axes and legend
ax.set_title("Beats: fast crests inside slow bulges (t = 6 s)", fontsize=10)   # the message
savefig(fig, "ch07", "c09_beats"); plt.show()   # save, then draw
""", see=r"""Fast crests (blue) inside slow bulges bounded by the purple envelope; in 6 s a crest (orange) moved 18.8 m but a node
(black) only 9.4 m; the rose dashed ghost is the printed envelope, which stays put (our remake of Fig. 7.13 `N179`).""",
    read=r"""Nodes travel at Δω/Δk: no energy crosses a node, so the energy travels at the node speed.""",
    change=r"""…Δk were halved: the groups become twice as long and the node speed tends to dω/dk.""")
nb.figure(r"""
k = np.linspace(0, 2, 300)                                # wavenumbers [rad/m]
fig, ax = plt.subplots(figsize=(5.5, 3.5))                # one panel
ax.plot(k, ch07.omega_gravity(k), color=COLORS["ink"], lw=2, label="ω = √(gk) (deep)")   # (7.28) deep
k0, w0, cg0 = 1.0, ch07.omega_gravity(1.0), ch07.group_velocity(1.0)   # the point k = 1 rad/m
ax.plot(k, w0/k0*k, color=COLORS["orange"], lw=1.5, label=f"chord: slope c = {w0/k0:.2f} m/s")        # ω/k
ax.plot(k, w0 + cg0*(k - k0), color=COLORS["accent"], lw=1.5, label=f"tangent: slope c_g = {cg0:.2f} m/s")   # dω/dk
ax.plot(k0, w0, "o", color=COLORS["ink"])
ax.set_ylim(0, 5); ax.set_xlabel("k [rad/m]"); ax.set_ylabel("ω [rad/s]"); ax.legend(fontsize=8)
ax.set_title("Crest speed = chord slope; group speed = tangent slope", fontsize=10)
savefig(fig, "ch07", "c09_chord_tangent"); plt.show()
""", see=r"""The deep-water ω(k) with, at k = 1 rad/m, the chord from the origin (orange, slope c = 3.13 m/s) and the tangent
(purple, slope $c_g=1.57$ m/s) — our remake of Fig. 7.14 `N180`.""",
    read=r"""$c_g<c$ wherever ω(k) bends down (gravity waves); $c_g>c$ where it bends up (capillary waves).""",
    change=r"""…shallow water: ω = k√(gH) is a straight line through the origin — chord and tangent coincide, $c_g=c$.""")
P("P172", "Fourier integral and a packet's spectrum", r"""
A real wave group is a continuous sum of sinusoids, $\eta(x,t)=\int A(k)e^{i(kx-\omega(k)t)}dk$ (real part understood),
where A(k) is its spectrum. A narrow spectrum of width δk around k₀ means a long group of length ~1/δk, and near k₀ we may
Taylor-expand ω(k) (Ch. 1 P26). This is Ch. 5 P142's discrete Fourier modes with the sum turned into an integral.
""", code=r"""
import numpy as np                                # numbers
x = np.linspace(-200, 200, 4001)                  # positions [m]
k = np.linspace(0.9, 1.1, 201); A = np.exp(-(k - 1)**2/(2*0.02**2))   # a narrow spectrum around k₀ = 1 rad/m
eta = (A[:, None]*np.cos(k[:, None]*x)).sum(0)*(k[1] - k[0])          # the integral as a sum over k
print(x[np.argmax(eta)], np.abs(eta[np.abs(x) > 150]).max()/eta.max())   # peak at 0; tails < 1 %: ~1/0.02 = 50 m long
""")
note("N66", "A wave packet", r"""
All wavenumbers in a narrow band δk around k; in space a nearly sinusoidal wave whose amplitude dies away over a length
∝ 1/δk. The figure is our remake of Fig. 7.15 `N181`.
""")
P("P181", "np.fft.rfft and np.fft.rfftfreq", r"""
For a **real** signal the FFT's negative-wavenumber half mirrors the positive half, so `np.fft.rfft(eta)` returns only the
modes with k ≥ 0 (a reminder of the FFT itself: Ch. 5 P142). `np.fft.rfftfreq(n, d)` gives the matching frequencies in
cycles per unit length for n samples spaced d apart; multiplying by 2π turns them into wavenumbers in rad/m. The size
`np.abs(...)` of each complex mode is the height of the spectrum.""", code=r"""
import numpy as np                                        # arrays
x = np.arange(8)*0.5                                      # 8 samples, 0.5 m apart (a 4 m periodic domain)
F = np.fft.rfft(np.cos(2*np.pi*x/2.0))                    # a cosine of wavelength 2 m
print(2*np.pi*np.fft.rfftfreq(8, d=0.5), np.abs(F).round(3))   # k = 0, 1.57, 3.14, 4.71, 6.28 rad/m; peak at k = π rad/m
""")
nb.figure(r"""
x = np.linspace(-150, 150, 4096, endpoint=False)          # a periodic domain [m]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.2))      # packet | spectrum
for sx, cl in ((30.0, COLORS["blue"]), (15.0, COLORS["teal"])):   # a long and a short packet [m]
    p = ch07.gaussian_packet(x, 0.0, 1.0, 1.0, sx)        # η = e^{−x²/2σx²} cos k₀x at t = 0
    if sx == 30.0:
        a1.plot(x, p["eta"], color=cl, lw=1); a1.plot(x, p["envelope"], color=COLORS["accent"], lw=2)
    kk = 2*np.pi*np.fft.rfftfreq(x.size, d=x[1] - x[0])   # wavenumbers of the FFT [rad/m]
    a2.plot(kk, np.abs(np.fft.rfft(p["eta"])), color=cl, lw=2, label=f"σ_x = {sx:.0f} m")   # |spectrum|
a1.set_xlabel("x [m]"); a1.set_ylabel("η [m]"); a1.set_title("a packet and its envelope (σ_x = 30 m)", fontsize=9)
a2.set_xlim(0.7, 1.3); a2.set_xlabel("k [rad/m]"); a2.set_ylabel("|FFT| [–]"); a2.legend(fontsize=7)
a2.set_title("a shorter packet has a wider spectrum", fontsize=9)
savefig(fig, "ch07", "c09_packet_spectrum"); plt.show()
""", see=r"""Left, a Gaussian packet (blue) inside its envelope (purple); right, the size of its Fourier spectrum around
k₀ = 1 rad/m, and that of a packet half as long (teal), twice as wide.""",
    read=r"""Length × width ≈ constant: for a Gaussian, σ_x σ_k = 1, so 30 m ↔ 0.033 rad/m and 15 m ↔ 0.067 rad/m.""",
    change=r"""…the packet were halved again: its spectrum doubles again — and it spreads faster in time (D20's dropped term).""")
remind([
    ("substitution in an integral", r"replace k by $k_0+\kappa$ and dk by dκ (Ch. 3 P106)."),
    ("complex amplitudes", "Re is taken at the end, so $e^{i\\theta}$ factors can be pulled out of integrals (primer P176 in C13 develops this; here only $\\mathrm{Re}[a\\,e^{i\\theta}]=a\\cos\\theta$ for real a is used, Ch. 1 P45)."),
])
D("D20", ref="7.68", check_src=r"""
import sympy as sp                                              # symbolic algebra
k0, kap, x, t, w0, cg, w2 = sp.symbols('k0 kappa x t omega0 c_g omega2', real=True)   # carrier, offset, position, time, ω₀, c_g, ω″
xi = sp.symbols('xi', real=True)                                # the moving coordinate ξ = x − c_g t
omega_taylor = w0 + cg*kap + w2*kap**2/2                        # step 5: ω(k₀ + κ) to second order
phase = (k0 + kap)*x - omega_taylor*t                           # steps 3 and 6: the phase kx − ωt with k = k₀ + κ
split = (k0*x - w0*t) + kap*(x - cg*t) - w2*kap**2*t/2          # step 6: carrier + envelope + small rest
print("step 6 regrouping:", sp.expand(phase - split))           # 0: the regrouping is exact
al, bt = sp.symbols('alpha beta', positive=True)                # a Gaussian spectrum e^{−ακ²} and a position β
gauss = sp.integrate(sp.exp(-al*kap**2 + sp.I*bt*kap), (kap, -sp.oo, sp.oo))    # ∫ e^{−ακ² + iβκ} dκ
print("Gaussian integral:", sp.simplify(gauss - sp.sqrt(sp.pi/al)*sp.exp(-bt**2/(4*al))))   # 0: = √(π/α) e^{−β²/4α}
s, tt = sp.symbols('s tau', positive=True)                      # spectral width δk and a (positive) time τ
Gauss = lambda A_, B_: sp.sqrt(sp.pi/A_)*sp.exp(-B_**2/(4*A_))  # the closed form just verified, ∫e^{−Aκ² + iBκ}dκ
a0 = Gauss(1/(2*s**2), x)                                       # step 4: a(x) for the spectrum A(k₀ + κ) = e^{−κ²/2s²}
env1 = Gauss(1/(2*s**2), x - cg*tt)                             # step 9 with ω″ dropped: the integral with κ(x − c_g t)
print("steps 9-10: envelope − a(x − c_g t) =", sp.simplify(env1 - a0.subs(x, x - cg*tt)),   # 0: it IS the initial shape …
      "; ∂/∂t + c_g ∂/∂x of it =", sp.simplify(sp.diff(env1, tt) + cg*sp.diff(env1, x)))   # … carried rigidly at c_g
alpha2 = 1/(2*s**2) + sp.I*w2*tt/2                              # keep ω″: the exponent −κ²/2s² − iω″κ²t/2, α complex
m2 = sp.simplify(sp.pi/sp.Abs(alpha2)*sp.exp(-xi**2*sp.re(1/alpha2)/2))   # |envelope|² = |√(π/α)|² |e^{−ξ²/4α}|²
print("|envelope|² =", m2)
print("slope of |envelope|² at ξ = 0:", sp.simplify(sp.diff(m2, xi).subs(xi, 0)))   # 0: the peak rides at x = c_g t (step 12)
print("peak ratio:", sp.simplify(m2.subs(xi, 0)/m2.subs({xi: 0, tt: 0})))   # 1/√(1 + (ω″s²t)²): spreading once t ~ 1/(ω″δk²), step 7
""")
note("N67", "The packet keeps its shape for a while", r"""
For short times the packet's envelope moves at $c_g$: $\eta=a(x-c_gt)\cos(kx-\omega t)$ *(7.68)* (the book cites Phillips;
D20 writes it out). The nodes of the envelope move at $c_g$ and no energy crosses them.
""")
remind([
    ("chain rule", r"$d\,\omega^2/dk=2\omega\,d\omega/dk$ (Ch. 1 P49); the product rule differentiates $k\tanh kH$ (Ch. 1 P38)."),
])
D("D21", ref="7.69")
note("N68–N70", "Group velocity of water waves", r"""
$c_g=\frac c2\Big[1+\frac{2kH}{\sinh(2kH)}\Big]$ *(7.69)*; $c_g=c/2\ \text{(deep water)},\ c_g=c\ \text{(shallow water)}$ *(7.70)*;
pure capillary waves have $c_g=3c/2$ (checked below with `group_velocity(g=0)`); and
$F=E\frac c2\big[1+\frac{2kH}{\sinh(2kH)}\big]=Ec_g$ *(7.71)*, $E=\rho ga^2/2$ — **the rate of energy transmission is energy
times group velocity**.
""")
note("N71", "In three dimensions", r"""
$c_{gi}=\partial\omega/\partial K_i$: the group velocity is the gradient of ω in wavenumber space (used for internal waves
in C16). For deep water with K = (0.6, 0.8) rad/m (|K| = 1) it points along K with size $\tfrac12\sqrt{g/K}=1.566$ m/s:
""", equation=r"\mathbf c_g=\nabla_{\mathbf K}\,\omega(\mathbf K)")
nb.code(r"""
print(ch07.group_velocity_vector(lambda K: np.sqrt(G*np.hypot(K[0], K[1])), np.array([0.6, 0.8])), "m/s")   # ∂ω/∂K_i
""")
nb.worked_example("swell from a distant storm", r"""
1. Deep water, λ = 156 m (T = 10 s): $c=\sqrt{g\lambda/2\pi}=15.6$ m/s.
2. $c_g=c/2=7.8$ m/s.
3. A storm 1000 km away: the first energy arrives after $10^6/7.8=128\,000$ s ≈ 35.6 h — not 17.8 h, as the crest speed
   would suggest.
4. Longer swell (T = 20 s) arrives first: $c_g=15.6$ m/s, 17.8 h — a forecaster reads the storm's distance from how fast
   the period falls.
""")
nb.code(r"""
k = 2*np.pi/156.0                                          # the swell of the tiny example [rad/m]
print(f"deep: c = {ch07.phase_speed(k):.2f}, c_g = {ch07.group_velocity(k):.2f} m/s")          # (7.29), (7.69)
print(f"on 30 m: c = {ch07.phase_speed(k, 30.0):.2f}, c_g = {ch07.group_velocity(k, 30.0):.2f} m/s")   # kH = 1.21
print(f"H = 10 m, λ = 50 m: c = {ch07.phase_speed(2*np.pi/50, 10.0):.2f}, c_g = {ch07.group_velocity(2*np.pi/50, 10.0):.2f} m/s")
print(f"complex-step dω/dk = {ch07.group_velocity_numeric(None, k):.10f} m/s")                 # numerical (7.67)
kc = 2*np.pi/0.005                                          # a 5 mm ripple [rad/m]
print("pure capillary c_g/c =", ch07.group_velocity(kc, sigma=0.0727, g=0.0)/ch07.phase_speed(kc, sigma=0.0727, g=0.0))
ps = ch07.packet_state(k, distance=1e6)                    # everything about this group, storm 1000 km away
print({n: (round(v, 3) if isinstance(v, float) else v) for n, v in ps.items()})
assert np.isclose(ch07.energy_flux(1.0, k), ch07.wave_energy_density(1.0)*ch07.group_velocity(k))   # (7.44) = E c_g, (7.71)
""", explain=r"""
1. Deep water: c = 15.61 and $c_g=7.80$ m/s (half). On a 30 m shelf (kH = 1.21) c = 14.27 m/s but $c_g=10.24$ m/s — larger
   than in deep water: in intermediate depth $c_g/c$ rises toward 1. For H = 10 m, λ = 50 m: 8.15 and 5.74 m/s (D21's check).
2. The complex-step derivative $\mathrm{Im}\,\omega(k+ih)/h$ has no subtraction, so no cancellation error (Ch. 3 P107):
   it equals (7.69) to ten digits.
3. Pure capillary waves: $c_g/c=1.5$ — groups outrun crests.
4. `packet_state` bundles c, $c_g$, E, F = E c_g and the arrival time from a storm 1000 km away (35.6 h); the assertion
   is (7.71): the flux of C06 is energy times group velocity.
""")
P("P173", "envelope with scipy.signal.hilbert", r"""
The modulus of the "analytic signal" `np.abs(scipy.signal.hilbert(eta))` follows the slowly varying amplitude of a fast
oscillation — the envelope — without fitting. We use it to track a packet's peak (`ch07.envelope` wraps it).
""", code=r"""
import numpy as np, scipy.signal as ss            # numbers and signal tools
x = np.linspace(0, 100, 2001); eta = np.exp(-((x - 50)/10)**2)*np.cos(3*x)   # a packet
env = np.abs(ss.hilbert(eta))                     # its envelope
print(x[np.argmax(env)], env.max())               # 50.0, ≈ 1.0
""")
remind([
    ("periodic FFT domain and wrap-around", "the FFT treats the domain as periodic: whatever leaves at the right re-enters at the left, so the domain must be long enough that the packet never reaches its end (Ch. 5 P142)."),
])
nb.md(r"""
**From scratch — two checks of $c_g$.** (1) A central difference of ω(k). (2) Evolve a Gaussian packet (k₀ = 1 rad/m,
σ_x = 30 m, deep water) with the FFT (`linear_evolve`, each mode with its own ω(k)) on a periodic domain long enough that
nothing wraps round, and fit the speed of the envelope's peak. The run is kept (`packet_run`) for the x–t diagram of C10.
""")
nb.check_agree(r"""
h = 1e-6                                                   # difference step [rad/m] (k = 0.040 rad/m here)
cg_mine =(ch07.omega_gravity(k + h) - ch07.omega_gravity(k - h))/(2*h)   # (1) central difference of (7.28), deep
assert np.allclose(cg_mine, ch07.group_velocity(k), rtol=1e-8)
N = 4096 if not FAST else 1024                             # FFT points (spectral: k₀ = 1 is resolved either way)
xp = np.linspace(0, 2000, N, endpoint=False)               # periodic domain [m]
eta0 = np.exp(-(xp - 300)**2/(2*30**2))*np.cos(1.0*(xp - 300))   # the packet at t = 0 [m]
tp = np.linspace(0, 200, 11)                               # eleven times [s]
packet_run = ch07.linear_evolve(eta0, xp, tp, direction=+1, H=np.inf)   # (7.28) mode by mode, one-way to +x
peaks = [xp[np.argmax(ch07.envelope(e))] for e in packet_run]   # envelope peak at each time [m]
slope = np.polyfit(tp, peaks, 1)[0]                        # its speed [m/s]
print(f"(1) central difference: {cg_mine:.6f} m/s;  (2) envelope peak speed: {slope:.4f} m/s;  c_g = {ch07.group_velocity(1.0):.4f} m/s")
assert np.allclose(slope, ch07.group_velocity(1.0), rtol=1e-2)   # the envelope moves at c_g, not at c = 3.13 m/s
""")
P("P182", "np.interp: reading a curve between its samples", r"""
`np.interp(x_new, x, y)` returns the value of the sampled curve (x, y) at the point x_new by drawing a straight line
between the two neighbouring samples (x must be increasing). Below it reads the envelope's height at the tracked crest,
which sits between grid points.""", code=r"""
import numpy as np                                        # arrays
print(np.interp(1.5, [0.0, 1.0, 2.0], [0.0, 10.0, 30.0])) # halfway between 10 and 30 → 20.0
""")
nb.animation(r"""
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 4.6))       # packet (top) | stone in a pond (bottom)
xa = np.linspace(-60, 260, 1200)                           # positions for the packet [m]
p0 = ch07.gaussian_packet(xa, 0.0, 1.0, 1.0, 20.0)         # deep-water packet, k₀ = 1 rad/m, σ_x = 20 m
c, cg = p0["c"], p0["cg"]                                  # 3.13 and 1.57 m/s
(car,) = a1.plot(xa, p0["eta"], color=COLORS["blue"], lw=1)                     # carrier
(env_u,) = a1.plot(xa, p0["envelope"], color=COLORS["accent"], lw=2)             # envelope
(dot,) = a1.plot([], [], "o", color=COLORS["orange"], ms=7)                      # one tracked crest
xcrest0 = -2*np.pi*6                                       # a crest that starts six wavelengths behind the centre [m]
a1.set_xlim(xa[0], xa[-1]); a1.set_ylim(-1.1, 1.1); a1.set_ylabel("η [m]")
txt = a1.text(-55, 0.85, "", fontsize=8)
xs_ = np.linspace(0, 1.0, 800)                             # pond: distance from the stone [m]
(pond,) = a2.plot(xs_, 0*xs_, color=COLORS["blue"], lw=1)
cgmin = ch07.min_group_velocity(0.07274, 998.2)["cg_min"]  # minimum group velocity of ripples [m/s]
calm = a2.axvline(0, color=COLORS["accent"], ls="--", lw=1)   # the calm-centre edge at c_g,min t
a2.set_xlim(0, 1); a2.set_ylim(-0.012, 0.012); a2.set_xlabel("x [m] (packet above: x in m too)"); a2.set_ylabel("η [m]")
nf = 48 if not FAST else 24                                # frames


def update(i):   # draw frame i
    t = i/nf*120.0                                         # packet time 0…120 s
    p = ch07.gaussian_packet(xa, t, 1.0, 1.0, 20.0)        # the packet at t (spreads slowly: ω″ kept)
    car.set_ydata(p["eta"]); env_u.set_ydata(p["envelope"])
    xc = xcrest0 + c*t                                     # the tracked crest moves at c
    amp = np.interp(xc, xa, p["envelope"])                 # the envelope height where it is
    dot.set_data([xc], [amp] if amp > 0.03 else [np.nan])  # visible only inside the group
    txt.set_text(f"t = {t:5.1f} s   crest at c = {c:.2f} m/s, group at c_g = {cg:.2f} m/s")
    tp_ = 0.02 + i/nf*2.0                                  # pond time 0.02…2 s
    pond.set_ydata(ch07.pond_ripples(xs_, tp_, width=0.01, a=0.01, rho=998.2, sigma=0.07274))   # stone-in-a-pond train
    calm.set_xdata([cgmin*tp_, cgmin*tp_])
    a2.set_title(f"stone in a pond, t = {tp_:.2f} s: calm inside c_g,min·t (purple)", fontsize=9)
    return []   # nothing to blit


show_animation(animate(update, frames=nf, fig=fig, interval=60))
""")
see_read_change(r"""Top: a deep-water packet; the orange dot rides one crest, enters at the rear of the group, runs through it
and disappears at the front while the purple envelope moves at half its speed. Bottom: the ripple train from a stone
(released hump 1 cm wide) spreading, with the calm-centre edge (purple dashed) — the `A2` animation.""",
                r"""A crest lives about as long as it takes to cross the group at the speed difference c − c_g. In the pond the
waves sort themselves by group velocity; nothing is left inside the dashed line, which moves at the minimum group
velocity.""",
                r"""…the water were shallow (ω = k√(gH)): crest and envelope would move together; the dot would stay at the
same place in the group forever.""")
remind([
    ("animate_figure (plotly time slider)", "like `slider_figure` (Ch. 1 P17) but the slider is time and a ▶ button plays it; every frame is computed up front, so it works on the page."),
])
nb.plotly(r"""
xq = np.linspace(-50, 400, 400)                           # positions [m]
p0 = ch07.gaussian_packet(xq, 0.0, 1.0, 1.0, 20.0)         # deep packet, k₀ = 1 rad/m
def frame(t):                                              # everything drawn at time t [s]
    p = ch07.gaussian_packet(xq, t, 1.0, 1.0, 20.0, order=2)   # carrier and envelope (ω″ kept)
    xg = p0["cg"]*t                                        # the group centre moves at c_g [m]
    xc = xg + (p0["c"] - p0["cg"])*t % 60 - 30             # a crest marker: speed c, re-entering the group at the rear [m]
    return {"η": (xq, p["eta"]), "envelope": (xq, p["envelope"]),
            "marker at c_g": ([xg, xg], [-1.2, 1.2]), "marker at c": ([xc, xc], [-1.2, 1.2])}
fig = animate_figure(frame, np.linspace(0, 200, 21 if not FAST else 11), xlabel="x [m]", ylabel="η [m]",
                     title="Crests at c, the group at c_g (deep water, k₀ = 1 rad/m)", yrange=(-1.3, 1.3))
for tr in list(fig.data) + [d for fr in fig.frames for d in fr.data]:   # colours by meaning
    tr.line.color = {"η": COLORS["blue"], "envelope": COLORS["accent"], "marker at c_g": COLORS["accent"],
                     "marker at c": COLORS["orange"]}[tr.name]
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** Press ▶: the purple bar stays at the centre of the envelope (speed $c_g=1.57$ m/s); the orange bar runs through the group
at c = 3.13 m/s relative to the ground (it re-enters at the rear each time it leaves the front, like a new crest being
born). This player works on the page without a kernel.

**What would change if…** the water were shallow: both bars would move together at √(gH) and the orange bar would never leave the
group.
""")
note("N72", "A stone in a pond", r"""
(the book's Fig. 7.16) The impact contains all wavelengths; each travels at its own $c_g$, so the train sorts itself —
longest gravity waves in front — and because capillary–gravity waves have a **minimum group velocity**
$c_{g,min}$ (0.1776 m/s at λ ≈ 4.35 cm for clean water at 20 °C, σ = 72.74 mN/m — printed by the next cell), nothing is
left inside a circle of radius $c_{g,min}t$: a calm centre. Heights fall as the train lengthens; viscosity finally damps
each wave like $a_0e^{-2\nu k^2t}$ (Ch. 8's tools).
""")
nb.code(r"""
print(ch07.min_group_velocity(0.07274, 998.2))            # c_g,min [m/s] and its wavenumber and wavelength
print("a 1 cm ripple after 10 s keeps", ch07.viscous_decay(1.0, 2*np.pi/0.01, 1e-6, 10.0), "of its height")   # e^{−2νk²t}
""", explain=r"""
$c_{g,min}=0.1776$ m/s at λ = 4.35 cm; the viscous factor for a 1 cm ripple over 10 s is $e^{-2\times10^{-6}\times394\,784
\times10}=e^{-7.9}=3.7\times10^{-4}$ — short ripples die fast.
""")
nb.figure(r"""
xs_ = np.linspace(0, 1.0, 1000)                            # distance from the stone [m]
cgmin = ch07.min_group_velocity(0.07274, 998.2)["cg_min"]  # [m/s]
fig, axs = plt.subplots(3, 1, figsize=(7, 4.2), sharex=True)   # three times, stacked
for axp, t in zip(axs, (0.5, 1.0, 2.0)):                   # times after the impact [s]
    axp.plot(xs_, 1e3*ch07.pond_ripples(xs_, t, width=0.01, a=0.01, rho=998.2, sigma=0.07274), color=COLORS["blue"], lw=1)   # the train at time t, in mm
    axp.axvline(cgmin*t, color=COLORS["accent"], ls="--", lw=1.5)   # calm-centre edge
    axp.text(0.8, 0.5, f"t = {t} s", transform=axp.transAxes, fontsize=8); axp.set_ylabel("η [mm]")   # time label and unit
    axp.set_ylim(-4, 4)   # same scale in every panel
axs[-1].set_xlabel("distance from the impact [m]")   # shared x label
axs[0].set_title("A stone in a pond: the train stretches, the centre stays calm", fontsize=10)   # the message
savefig(fig, "ch07", "c09_pond"); plt.show()   # save, then draw
""", see=r"""The ripple train from a released 1 cm hump at 0.5, 1 and 2 s (our remake of Fig. 7.16 `N182`); the purple dashed line
is at $c_{g,min}t$.""",
    read=r"""The longest waves lead; inside the dashed line the surface stays flat — its edge moves at the minimum group
velocity (0.1776 m/s, printed two cells above), not at the minimum phase speed $c_{min}=23.1$ cm/s of C07.""",
    change=r"""…σ = 0 (no surface tension): no minimum group velocity, and ever shorter waves fill the centre.""")
note("N73", "Following one crest", r"""
in a deep-water group: it runs from the rear to the front at twice the group speed, its wavelength grows as it goes, and it
dies at the front while new crests are born at the rear (seen in the A2 animation and in the explainer).
""")
explainer("group_velocity_packets", "Where does a wave's energy go?",
          r"""Follow one crest with a marker while the envelope moves at half its speed; switch the dispersion to shallow or
capillary and see $c_g=c$ or $c_g>c$; the chord and tangent of ω(k) move with the dot. Crests passing through a group can
only be seen moving.""", "",
          ["Play the *swell set*: count how many crests pass through one group.",
           "Choose *shallow*: the packet keeps its shape — no dispersion.",
           "Choose *capillary*: crests now appear at the front and vanish at the back.",
           "Open the Derivation (D20) and toggle the dropped ω″ term: the packet spreads."])
confusion(r"""
"the wave speed" is two speeds. $c=\omega/k$ is how fast crests move; $c_g=d\omega/dk$ is how fast the envelope, the energy
and any signal move. For deep-water gravity waves $c_g=c/2$; for capillary waves $c_g=3c/2$; only non-dispersive waves
(shallow water, sound) have $c_g=c$.
""")
whatif(r"""
…the medium changed slowly along the way (depth decreasing toward a beach)? k and ω can no longer both stay constant: C10
shows that ω is the one that is carried unchanged along a path moving at $c_g$.
""")

# ---- C10 ------------------------------------------------------------------------------------------------------------
core("C10", r"Kinematic wave theory, rays and refraction", r"""
How do the wavenumber and the frequency of a wave train change when the depth changes slowly — and why do waves turn
toward the shore?
""", eqs=("7.79",))
problem(r"""
Walk along a beach on any day: the breakers arrive nearly parallel to the shore, whatever direction the wind blew offshore.
Waves round an island even bend toward its sides. We want the rules that carry a wave train's wavelength, frequency and
direction across slowly changing depth — **ray theory**, also called **WKB theory** (after Wentzel, Kramers and
Brillouin): the approximation for a medium that changes little over one wavelength, so that locally the wave looks like a
plain sinusoid whose k and ω drift slowly. D22–D23 below are an example of it, and Ch. 13 uses the same theory for
internal and planetary waves.
""")
nb.md(r"""
#### The idea

```
phase θ(x,t):  k = ∂θ/∂x (crests per metre),  ω = −∂θ/∂t (crests per second passing a point)
crests are conserved  ⇒  ∂k/∂t + ∂ω/∂x = 0
uniform depth:  k is carried at c_g (rays are straight lines in the x–t plane)
depth H(x):     ω is carried at c_g, k changes;  the shoreward end of a crest slows  ⇒  the crest turns
```
""")
note("N74, N75", "Local wavenumber and frequency", r"""
A slowly varying train $\eta=a(x,t)\cos[\theta(x,t)]$ *(7.72)* has a local wavenumber and frequency
$k\equiv\partial\theta/\partial x,\ \omega\equiv-\partial\theta/\partial t$ *(7.73)* (for θ = kx − ωt they are the usual k
and ω). The code differentiates a chirped phase θ = 0.5x + 0.002x² − 2t numerically: k = 0.54 rad/m at x = 10 m, ω = 2
rad/s.
""")
nb.code(r"""
theta = lambda x, t: 0.5*x + 0.002*x**2 - 2.0*t           # a chirped wave train's phase [rad]
print(ch07.local_wavenumber_frequency(theta, 10.0, 0.0))  # (7.73): (k [rad/m], ω [rad/s]) at x = 10 m, t = 0
""")
P("P174", "first-order wave equation and characteristics", r"""
The equation ∂q/∂t + c ∂q/∂x = 0 says q does not change for an observer moving at speed c: along each line dx/dt = c in
the x–t plane (a characteristic) q is constant. It is the material derivative of Ch. 3 with the velocity replaced by c.
""", code=r"""
import numpy as np                                # numbers
q0 = lambda x: np.exp(-x**2)                      # a bump at t = 0
c, t = 2.0, 3.0                                   # speed [m/s] and time [s]
q = lambda x, t: q0(x - c*t)                      # the solution: the bump shifted by ct
print(q(6.0, 3.0), q0(0.0))                       # 1.0 1.0 — the value rode along x = ct
""")
remind([
    ("mixed partial derivatives commute", r"$\partial^2\theta/\partial t\,\partial x=\partial^2\theta/\partial x\,\partial t$ for smooth θ — Schwarz's theorem (Ch. 4 P121)."),
    ("chain rule along a path", "ω depends on x only through k, so ∂ω/∂x = (dω/dk)(∂k/∂x) (Ch. 3 P91)."),
])
D("D22", ref="7.75")
note("N76, N77", "Wavenumbers ride at the group velocity", r"""
Crest conservation $\partial k/\partial t+\partial\omega/\partial x=0$ *(7.74)* and, with a dispersion relation,
$\frac{\partial k}{\partial t}+c_g\frac{\partial k}{\partial x}=0$ *(7.75)*: **wavenumbers are carried at the group
velocity** — an observer moving at $c_g$ always sees the same wavelength. The x–t diagram below (our remake of Fig. 7.17
`N183`) uses the FFT packet of C09.
""")
nb.figure(r"""
tq = np.linspace(0, 200, 121 if not FAST else 61)          # times [s]
run = ch07.linear_evolve(eta0, xp, tq, direction=+1, H=np.inf)   # the C09 packet, now at many times
sel = (xp > 150) & (xp < 750)                              # the part of the domain the packet visits [m]
fig, ax = plt.subplots(figsize=(7, 3.8))                   # x–t diagram
env_xt = ch07.envelope(run, axis=-1)                       # the packet's envelope at every time [m]
ax.pcolormesh(xp[sel], tq, env_xt[:, sel], cmap="Greys", shading="auto", vmin=0, vmax=1.2)   # the group (grey band)
c0, cg0 = ch07.phase_speed(1.0), ch07.group_velocity(1.0)  # 3.13 and 1.57 m/s
crest_x, crest_t = [], []                                  # one crest followed through the group
x_prev, lost = None, False                                 # its position at the previous time [m]; True once it has faded
for j, t_ in enumerate(tq):                                # crest positions at each time: local maxima inside the group
    e = run[j]
    ic = np.nonzero((e[1:-1] > e[:-2]) & (e[1:-1] >= e[2:]) & (e[1:-1] > 0.15))[0] + 1   # local maxima (P180)
    if ic.size and not lost:                               # follow the crest that starts about 12 m behind the centre
        guess = 288.0 if x_prev is None else x_prev + c0*(tq[1] - tq[0])   # where it should be now: it moves at c [m]
        jn = ic[np.argmin(np.abs(xp[ic] - guess))]         # the nearest crest
        if abs(xp[jn] - guess) < 2.0:                      # still the same crest
            crest_x.append(xp[jn]); crest_t.append(t_); x_prev = xp[jn]   # record it
        else:
            lost = True                                    # it has faded at the front of the group: stop following
    elif ic.size == 0 and x_prev is not None:
        lost = True                                        # no crest left above the threshold
for n_ in range(-2, 64, 2):                                # every 2nd crest of the carrier cos(x − 300): crest at 300 − 2πn at t = 0
    xl = 300 - 2*np.pi*n_ + c0*tq                          # its path if it moves at c [m]
    amp = np.array([np.interp(xl[j], xp, env_xt[j]) for j in range(tq.size)])   # the group's height along that path (P182)
    ax.plot(np.where(amp > 0.15, xl, np.nan), tq, color=COLORS["orange"], lw=0.8)   # drawn only while inside the group
ax.plot(crest_x, crest_t, color=COLORS["orange"], lw=3, label="one crest found in the computed surface")   # one crest's whole life
ax.plot([], [], color=COLORS["orange"], lw=0.8, label="lines of slope c inside the group (every 2nd crest)")   # legend entry
for x0 in (240, 300, 360):                                 # rays: lines of slope c_g (purple)
    ax.plot(x0 + cg0*tq, tq, color=COLORS["accent"], lw=2)
ax.plot(250 + c0*tq, tq, color=COLORS["orange"], lw=1, ls="--")   # a guide of slope c for comparison
ax.set_xlim(150, 750); ax.set_ylim(0, 200); ax.set_xlabel("x [m]"); ax.set_ylabel("t [s]")
ax.legend(fontsize=7, loc="upper left")                   # the followed crest
ax.set_title("x–t diagram: crest lines (slope c) cross the rays (slope c_g)", fontsize=10)
savefig(fig, "ch07", "c10_xt_diagram"); plt.show()
""", see=r"""The packet is the grey band. The thick orange line is one crest found in the computed surface (local maxima,
time level by time level) from near the rear of the group until it fades at the front; the thin orange segments are lines
of slope c drawn through every second crest, only where the group is present. All are flatter than the band (speed
c = 3.13 m/s, like the dashed guide), so crests enter at the rear of the band and leave at its front, where they fade. The band itself
follows the thick purple rays (speed $c_g=1.57$ m/s) and widens slowly with time.""",
    read=r"""Along a purple line the wavenumber and frequency stay fixed; along an orange line a crest keeps its phase and
ends when it leaves the band at the front.""",
    change=r"""…shallow water: both families have the same slope — crests and groups move together.""")
note("N78–N80", "A slowly varying depth", r"""
In slowly varying depth $\omega=\sqrt{gk\tanh[kH(x)]}$ has the form $\omega=\omega(k,x)$ *(7.76)*, a local group velocity is
$\partial\omega(k,x)/\partial k=c_g$ *(7.77)*, and multiplying by ∂k/∂t gives $c_g\frac{\partial k}{\partial t}=\frac
{\partial\omega}{\partial t}$ *(7.78)*.
""")
D("D23", ref="7.79")
P("P175", "Snell's law for waves", r"""
When a wave crosses a medium whose properties change only in one direction (say x), the wavenumber component along the
other direction (y) cannot change: the crests must match along every line of constant x. With |k| changing, the angle α
between k and the x-direction obeys |k| sin α = constant — the optics rule n₁ sin α₁ = n₂ sin α₂.
""", code=r"""
import numpy as np                                # numbers
k1, a1 = 0.07, np.radians(30)                     # offshore wavenumber [rad/m] and angle
k2 = 0.18                                         # a larger k closer to shore (shallower water) [rad/m]
print(np.degrees(np.arcsin(k1*np.sin(a1)/k2)))    # 11.2° — the crest has turned toward the shore
""")
D("D24")                                                   # Snell's law: no book number (the book argues in words)
note("N47, N48", "Refraction on a beach and round an island", r"""
On a sloping beach (the book's Fig. 7.8) ω is fixed along the path, while $c=\sqrt{gH}$ and λ shrink as H falls, so
$k\sin\alpha$ constant (D24) forces α toward zero — crests end up parallel to the depth contours. Round a circular island
with a gently sloping beach (Fig. 7.9) the same bending brings crests onto the shadow side too — optics' lens, made of
water depth.
""")
nb.worked_example("swell turning on a 1:50 beach", r"""
T = 8 s, arriving at 30° to the shore normal in 20 m of water.

1. $\omega=2\pi/8=0.785$ rad/s (fixed along the ray).
2. At H = 20 m solve (7.28): k = 0.0708 rad/m (λ = 88.8 m).
3. At H = 2 m: k = 0.181 rad/m (λ = 34.7 m).
4. Snell: $\sin\alpha_2=0.0708\times\sin30°/0.181=0.195$ → α₂ = 11.3°.
5. The crest turned by almost 19° while crossing 900 m of beach (from 20 m to 2 m depth on a 1:50 slope).
""")
remind([
    ("np.arcsin, np.radians, np.degrees", "inverse sine in radians, and conversions degrees ↔ radians (numpy works in radians; Ch. 2 P66)."),
])
nb.code(r"""
print(ch07.crest_conservation_residual(theta, 10.0, 0.0))          # (7.74) for the chirp: residual ≈ 1e-10
r = ch07.snell_ray_plane_beach(np.linspace(1000, 100, 10), np.radians(30), 1000.0, 8.0, 0.02)   # T = 8 s, 30°, slope 1:50
print("α along the ray [deg]:", np.round(np.degrees(r["alpha"]), 2))
print("k sin α [rad/m]:", r["snell"])                              # Snell's invariant: constant
s = ch07.refraction_state(8.0, np.radians(30), 20.0, 2.0)          # one point: H₀ = 20 m → H = 2 m
print(f"at H = 2 m: α = {s['alpha_deg']:.2f}°, λ = {s['lam']:.1f} m")
ray = ch07.ray_trace(None, (1000.0, 0.0), (-0.0708*np.cos(np.radians(30)), 0.0708*np.sin(np.radians(30))), (0, 250),
                     H_fn=lambda X: 0.02*X[0])                     # Hamilton's ray equations over the same beach
print("relative drift of ω along the traced ray:", ray["omega_drift"])
""", explain=r"""
1. The chirp satisfies crest conservation to about $10^{-10}$ — it holds for any smooth phase.
2. `snell_ray_plane_beach` fixes ω = 2π/8 along the ray, gets k(x) from (7.28) with `brentq` and α from Snell's law: α
   falls from 30° at 20 m depth to 11.3° at 2 m, while k sin α stays 0.03538 rad/m.
3. `ray_trace` integrates the ray equations $d\mathbf x/dt=\nabla_{\mathbf k}\omega$, $d\mathbf k/dt=-\nabla_{\mathbf x}
   \omega$ (D23 step 8, D24 step 1) with `solve_ivp`: ω drifts by $2\times10^{-9}$ relative over the whole beach (7.79).
""")
nb.check_agree(r"""
wf = lambda kx, ky, x: np.sqrt(G*np.hypot(kx, ky)*np.tanh(np.hypot(kx, ky)*0.02*x))   # ω(k, x) on the 1:50 beach [rad/s]
def rhs(y):                                                # Hamilton's ray equations: y = (x, y, kx, ky)
    x, yy, kx, ky = y
    hk, hx = 1e-7, 1e-2                                    # difference steps in k [rad/m] and x [m]
    dwdkx = (wf(kx + hk, ky, x) - wf(kx - hk, ky, x))/(2*hk)   # ∂ω/∂k_x = c_g,x
    dwdky = (wf(kx, ky + hk, x) - wf(kx, ky - hk, x))/(2*hk)   # ∂ω/∂k_y = c_g,y
    dwdx = (wf(kx, ky, x + hx) - wf(kx, ky, x - hx))/(2*hx)    # ∂ω/∂x (the depth changes with x)
    return np.array([dwdkx, dwdky, -dwdx, 0.0])            # dx/dt, dy/dt, dk_x/dt, dk_y/dt
k0 = ch07.wavenumber_from_omega(2*np.pi/8, 20.0)           # start: T = 8 s at H = 20 m [rad/m]
y = np.array([1000.0, 0.0, -k0*np.cos(np.radians(30)), k0*np.sin(np.radians(30))])   # heading shoreward at 30°
dt = 0.5                                                   # RK4 step [s]
while y[0] > 110.0:                                        # march until H ≈ 2.2 m
    s1 = rhs(y); s2 = rhs(y + dt/2*s1); s3 = rhs(y + dt/2*s2); s4 = rhs(y + dt*s3)   # RK4 (P95)
    y = y + dt*(s1 + 2*s2 + 2*s3 + s4)/6
alpha_mine = np.arcsin(y[3]/np.hypot(y[2], y[3]))          # angle to the shore normal at the end [rad]
r_alpha = float(np.ravel(ch07.snell_ray_plane_beach(np.array([y[0]]), np.radians(30), 1000.0, 8.0, 0.02)["alpha"])[0])   # closed form there
print(f"x = {y[0]:.1f} m: our α = {np.degrees(alpha_mine):.6f}°, Snell {np.degrees(r_alpha):.6f}°;",
      f"ω drift {abs(wf(y[2], y[3], y[0])/(2*np.pi/8) - 1):.1e}")
assert np.allclose(alpha_mine, r_alpha, rtol=1e-6)         # our RK4 ray bends exactly as Snell's law says
""")
nb.figure(r"""
fig, axs = plt.subplots(1, 3, figsize=(10, 3.6))           # beach | island | x–t
# (a) plane beach, 1:50, rays from 1000 m offshore at 30°
xr = np.linspace(1000, 40, 60)                             # offshore distance along the ray [m]
ray0 = ch07.snell_ray_plane_beach(xr, np.radians(30), 1000.0, 8.0, 0.02)   # one ray (the others are shifted copies)
for y0 in np.linspace(-1500, 300, 7):                      # seven rays, shifted along the shore [m]
    axs[0].plot(y0 + ray0["y"], xr, color=COLORS["accent"], lw=1)
for H in (2, 5, 10, 15, 20):                               # depth contours [m]
    axs[0].axhline(H/0.02, color=COLORS["muted"], lw=0.5); axs[0].text(-1580, H/0.02 + 10, f"{H} m", fontsize=6, color=COLORS["muted"])
yy = np.linspace(-1600, 900, 250)                          # along-shore coordinate [m]
k_x = np.sqrt(ray0["k"]**2 - ray0["l"]**2)                 # shoreward wavenumber component on the ray [rad/m]
phase_x = np.concatenate([[0.0], np.cumsum(0.5*(k_x[1:] + k_x[:-1])*(xr[:-1] - xr[1:]))])   # ∫k_x dx toward the shore
PH = ray0["l"]*yy[None, :] + phase_x[:, None]              # phase θ = l y + ∫k_x dx on the grid
lv0 = 4*np.pi*np.floor(PH.min()/(4*np.pi))                 # the first crest level below the smallest phase
axs[0].contour(yy, xr, PH, levels=np.arange(lv0, PH.max(), 4*np.pi), colors=[COLORS["orange"]], linewidths=0.6)   # every 2nd crest
axs[0].set_aspect("equal"); axs[0].set_xlim(-1600, 900); axs[0].set_ylim(40, 1000)   # true angles
axs[0].set_xlabel("along the shore y [m]"); axs[0].set_ylabel("offshore x [m]"); axs[0].set_title("(a) 1:50 beach, T = 8 s", fontsize=9)
# (b) circular island: H rises from 0 at r = 2 km to 30 m at r = 5 km (our bathymetry)
R0 = 2000.0                                                # island radius [m]
Hf = lambda X: float(np.clip(30*(np.hypot(X[0], X[1]) - R0)/3000.0, 0.0, 30.0))   # depth [m]
kI = ch07.wavenumber_from_omega(2*np.pi/12, 30.0)          # T = 12 s on 30 m [rad/m]
th_ = np.linspace(0, 2*np.pi, 200); axs[1].fill(R0*np.cos(th_)/1e3, R0*np.sin(th_)/1e3, color="#c8b99a")
axs[1].plot((R0 + 3000)*np.cos(th_)/1e3, (R0 + 3000)*np.sin(th_)/1e3, color=COLORS["muted"], lw=0.5, ls="--")
n_lee = 0                                                  # rays that land on the island's east (lee) half
y_start = np.concatenate([[-6000, -5500], np.linspace(-4600, 4600, 15 if not FAST else 9), [5500, 6000]])   # across the shelf, plus 2 beyond its edge on each side [m]
for yI in y_start:                                         # parallel rays from the west
    rI = ch07.ray_trace(None, (-9000.0, yI), (kI, 0.0), (0, 1800), H_fn=Hf, n_out=300)   # ray equations (stop at H ≈ 0)
    xe, ye = rI["x"][0][-1], rI["x"][1][-1]                # where the ray ends [m]
    lee = np.hypot(xe, ye) < R0 + 50 and xe > 0            # it reached the shore on the side facing away from the swell
    n_lee += lee                                           # count it
    axs[1].plot(rI["x"][0]/1e3, rI["x"][1]/1e3, color=COLORS["orange"] if lee else COLORS["accent"], lw=2 if lee else 1)
print(f"(b) rays that wrap round onto the lee (east) half of the island: {n_lee}")
axs[1].set_aspect("equal"); axs[1].set_xlim(-9, 9); axs[1].set_ylim(-7, 7)
axs[1].set_xlabel("x [km]"); axs[1].set_ylabel("y [km]"); axs[1].set_title("(b) island, T = 12 s", fontsize=9)
# (c) x–t of four rays straight onto the beach, different periods
for T_, cl in zip((6.0, 8.0, 10.0, 12.0), ("#c4b5fd", "#8b5cf6", "#6d28d9", "#4c1d95")):
    kT = ch07.wavenumber_from_omega(2*np.pi/T_, 20.0)      # start at 20 m depth [rad/m]
    rT = ch07.ray_trace(None, (1000.0, 0.0), (-kT, 0.0), (0, 600), H_fn=lambda X: 0.02*X[0], n_out=300)
    axs[2].plot(rT["t"], rT["x"][0], color=cl, lw=2, label=f"T = {T_:g} s, ω drift {rT['omega_drift']:.0e}")
axs[2].set_xlabel("t [s]"); axs[2].set_ylabel("offshore x [m]"); axs[2].legend(fontsize=6); axs[2].set_title("(c) rays in x–t", fontsize=9)
fig.suptitle("Frequency is carried along rays; the wavelength and direction change", fontweight="bold")
savefig(fig, "ch07", "c10_refraction"); plt.show()
""", see=r"""(a) Rays (purple) coming from offshore at 30° bend toward the shore normal as the depth contours (grey) get
shallower, and the crest lines (orange) turn parallel to the shore — our remake of Fig. 7.8 `N174`. (b) Parallel rays
from the west bending toward a circular island from every side: the central rays land on its west face, the outer ones
curve round its north and south points, and the rays launched about 3.5–4.5 km off the axis (thick orange, counted in the
printed line) wrap far enough to land on its east, sheltered half — Fig. 7.9 `N175` with our bathymetry. The outermost
shelf rays (launched ±4.6 km off the axis, thin purple) bend the most of all: they swing behind the island and cross each
other on its east side, but they never come closer than about 3.8 km to the centre, so they never reach the shore and are
not counted — the traced 30 minutes end while they are still turning (traced longer, they head back out to sea). The four rays launched at ±5.5 and ±6 km never cross the shelf edge (dashed circle): the
depth along them stays 30 m, so they pass by dead straight. (c) Four rays straight
onto the beach in the x–t plane: curved, because they slow down as the water shallows, each keeping its ω to about
$10^{-9}$ — Fig. 7.18 `N184`.""",
    read=r"""Each ray keeps its ω (printed drift); k grows as the depth falls, so the ray bends toward the shallow side. Longer
periods (darker) travel faster in deep water and slow down more near the shore.""",
    change=r"""…a longer period (T = 14 s): the swell feels the bottom sooner, so the turning starts further offshore.""")
nb.plotly(r"""
xr = np.linspace(1000, 40, 40)                             # offshore distance [m]
def rays(a0):                                              # rays over the 1:50 beach for incidence a0 [deg]
    r_ = ch07.snell_ray_plane_beach(xr, np.radians(a0), 1000.0, 8.0, 0.02)   # Snell with (7.28)
    return {f"ray {j}": (y0 + r_["y"], xr) for j, y0 in enumerate(np.linspace(-1500, 300, 4))}   # four rays: shifted copies along the shore
fig = slider_figure(rays, "α₀", np.linspace(0, 60, 13 if not FAST else 7), unit="°", xlabel="along the shore y [m]",   # slider over the incidence angle
                    ylabel="offshore x [m]", title="Rays over a 1:50 beach (T = 8 s): all end nearly normal to the shore",   # labels and message
                    xrange=(-1600, 1400), yrange=(0, 1000))   # fixed axes
recolor(fig, {f"ray {j}": COLORS["accent"] for j in range(4)})   # rays purple (group/energy colour)
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** Whatever the incidence angle α₀, the rays end nearly perpendicular to the shoreline (the bottom of the plot): refraction
straightens them.

**What would change if…** the beach were steeper (1:10): the turning would happen over a fifth of the distance, so the rays would bend
more sharply but end at nearly the same angle, set by Snell's law and the final depth.
""")
explainer("wave_rays_refraction", "Why do waves arrive parallel to the beach?",
          r"""Set the angle and the period and watch rays bend over a beach, a ridge or round an island; the side view prints
ω flat along the selected ray while k, c and $c_g$ change. The conservation is seen, not asserted.""", "",
          ["Preset *swell at 30°*: read α at the 2 m contour (11°).",
           "Switch to *island*: find a ray that reaches the shadow side.",
           "Choose *homogeneous*: the rays are straight — nothing to refract.",
           "Step the Derivation (D23) to 'multiply by c_g'."])
whatif(r"""
…the medium varied in time too (a tide changing the depth)? Then ω is no longer conserved along a ray: dω/dt = ∂ω/∂t at
fixed k — the Doppler-like frequency shift used for waves on tidal currents and for internal waves in Ch. 13.
""")

# =====================================================================================================================
# A.6 §7.6 — R10, C11 (N81–N85, N93–N97, N185, N186, N189, D25, D26, A5, E7), C12 (N86–N92, N187, N188, D27, A4, E2)
# =====================================================================================================================
nb.section("7.6", "Nonlinear Waves in Shallow and Deep Water", intro=r"""
**What is this section about?** Everything so far assumed ka ≪ 1. Now the amplitude matters. In shallow water crests
outrun troughs, fronts steepen and break into a hydraulic jump — unless dispersion balances the steepening and a solitary
wave travels unchanged. In deep water the crests sharpen (Stokes waves) and the orbits of C05 fail to close: the water
drifts forward (Stokes drift).
""")
nb.recap("R10", "The Froude number", r"""
$\mathrm{Fr}=\frac{U}{\sqrt{gl}}$ *(4.104)* compares a flow speed with the shallow-water wave speed of C04,
$c=\sqrt{gH}$ *(7.49)*: Fr > 1 supercritical (waves cannot travel upstream), Fr < 1 subcritical.
""", where="Ch. 4 §4.11")
nb.code(r"""
print(similarity.froude_number(2.97, 0.1, g=9.81))        # 2.97 m/s in 10 cm of water: Fr ≈ 3, supercritical
""", explain=r"""
$\mathrm{Fr}=u/\sqrt{gH}$ with u = 2.97 m/s, g = 9.81 m/s² and H = 0.1 m: 2.97/0.990 = 3.0 — the stream runs three times faster than a long wave can travel on
it, so no signal can move upstream — the condition for a jump.
""")
core("C11", r"The hydraulic jump — and the other fate, the solitary wave", r"""
Given the incoming depth and speed, how high is the jump, where does the energy go — and when does a steep shallow-water
wave not break at all?
""", eqs=("7.81",))
problem(r"""
At the foot of a dam spillway a thin, fast sheet of water suddenly rises into a deep, slow, churning roller; a tidal bore
climbs up a river mouth as a moving step. Both are hydraulic jumps. We want the downstream depth from the upstream depth
and speed, and the price paid in energy. Some steep waves, though, never break: a single hump can travel for kilometres
unchanged (Russell's solitary wave, 1844).
""")
note("N81", "Why shallow waves steepen", r"""
Relative to water moving at u, a small wavelet travels at the local $c'=\sqrt{gH'}$; in the fixed frame at $c=c'+u$.
Under a crest H′ and u are both larger, so crests outrun troughs and the front steepens (the book's Fig. 7.19). The book
stays qualitative; our labelled extension uses the exact simple wave of the shallow-water equations,
$c=3\sqrt{g(H+\eta)}-2\sqrt{gH}$ (a Riemann invariant of the shallow-water equations, Ch. 13/15; stated, not derived),
which breaks at $t_b=-1/\min(\partial c/\partial x)$ at t = 0 (characteristics, primer
P174; Ch. 15 develops them for gases).
""")
remind([
    ("simple wave and breaking time (our extension)", r"each surface point keeps its height and moves at its own speed c(η) along a characteristic (primer P174); where faster points catch slower ones the profile would fold — at $t_b$."),
])
nb.code(r"""
x = np.linspace(0, 200, 2001)                              # positions [m]
eta0 = 0.2*np.exp(-((x - 50)/10)**2)                      # a 20 cm hump on H = 1 m [m]
s = ch07.simple_wave_evolve(eta0, x, np.linspace(0, 12, 7), 1.0)   # every surface point rides its characteristic
print(f"breaking time t_b = {s['t_break']:.2f} s")
eta_test = np.array([-0.1, 0.0, 0.1])                      # a trough, the still level and a crest [m]
print("book c′ + u:", ch07.nonlinear_wavelet_speed(eta_test, 1.0, model="book"), "m/s")    # √(g(H + η)) + c₀η/H
print("simple wave:", ch07.nonlinear_wavelet_speed(eta_test, 1.0), "m/s")                  # 3√(g(H + η)) − 2√(gH)
""", explain=r"""
1. `simple_wave_evolve` moves every point of the hump at its own speed $c(\eta)$ and reports when the profile first
   becomes vertical: $t_b\approx13$ s for this 20 cm hump on 1 m of water.
2. The two wavelet-speed models agree to first order in η/H (2.66/2.65, 3.13, 3.60/3.59 m/s): both are
   $c_0(1+\tfrac32\eta/H)$ plus smaller terms — the coefficient $\tfrac32$ reappears in the KdV equation (7.87) below.
""")
nb.figure(r"""
fig, ax = plt.subplots(figsize=(7, 3))                     # one panel
tb = s["t_break"]                                          # breaking time [s]
for t_, cl in zip((0.0, 4.0, 8.0, tb), ("#bfdbfe", "#60a5fa", "#2563eb", "#1e3a8a")):
    xs_ = ch07.simple_wave_evolve(eta0, x, t_, 1.0)["x_points"]   # positions of the surface points at t_ [m]
    ax.plot(xs_, eta0, color=cl, lw=2, label=f"t = {t_:.1f} s")
ax.set_xlim(20, 110); ax.set_ylim(-0.03, 0.23); ax.set_xlabel("x [m]"); ax.set_ylabel("η [m]"); ax.legend(fontsize=7)
ax.set_title("Crests outrun troughs: the front becomes vertical at t_b", fontsize=10)
savefig(fig, "ch07", "c11_steepening"); plt.show()
""", see=r"""The same hump at four times: its front steepens and becomes vertical at $t_b\approx13$ s (our remake of Fig. 7.19
`N185`). Every surface point keeps its height and moves at its own speed c(η), so the top of the hump gains on the
foot of its front.""",
    read=r"""Points higher up move faster: after $t_b$ the curve would fold over — real water breaks instead.""",
    change=r"""…the hump were halved: $t_b$ roughly doubles (∝ 1/amplitude).""")
note("N82", "The hydraulic jump", r"""
After breaking, the front becomes a hydraulic jump: a nearly steady step joining a shallow supercritical stream (depth H₁,
wave speed $\sqrt{gH_1}<u_1$) to a deeper subcritical one. It may stand still (spillway, a kitchen sink's circular jump) or
move (a tidal bore). In the frame of a moving bore the same relations hold (the figure is our Fig. 7.20 `N186`).
""")
nb.figure(r"""
fig, ax = plt.subplots(figsize=(7, 3))                     # one panel
cv_box(ax, 0.1, 0.377)                                      # channel, jump roller and the dashed control volume
ax.annotate("", xy=(0.28, 0.05), xytext=(0.1, 0.05), arrowprops=dict(arrowstyle="->", color=COLORS["ink"], lw=2))   # fast inflow
ax.text(0.12, 0.07, "u₁", fontsize=10)   # label
ax.annotate("", xy=(0.9, 0.19), xytext=(0.8, 0.19), arrowprops=dict(arrowstyle="->", color=COLORS["ink"], lw=1.2))   # slow outflow
ax.text(0.8, 0.21, "u₂", fontsize=10)   # label
ax.annotate("", xy=(0.11, 0.033), xytext=(-0.04, 0.033), arrowprops=dict(arrowstyle="->", color=COLORS["amber"], lw=2))   # pressure force on face 1
ax.text(-0.04, 0.0, "½ρgH₁²", fontsize=8, color=COLORS["amber"])   # label
ax.annotate("", xy=(0.95, 0.125), xytext=(1.12, 0.125), arrowprops=dict(arrowstyle="->", color=COLORS["amber"], lw=3))   # pressure force on face 2
ax.text(0.98, 0.14, "½ρgH₂²", fontsize=8, color=COLORS["amber"])   # label
ax.set_xlim(-0.08, 1.15); ax.set_ylim(-0.02, 0.52); ax.set_xlabel("x (schematic)"); ax.set_ylabel("depth [m]")   # axes
ax.set_title("The control volume round a jump: H₁ = 0.1 m → H₂ = 0.377 m", fontsize=10)   # the message
savefig(fig, "ch07", "c11_jump_cv"); plt.show()   # save, then draw
""", see=r"""A control volume (dashed) round the jump with two vertical faces: fast shallow water in, slow deep water out, the
hydrostatic pushes ½ρgH² on the faces (amber).""",
    read=r"""Only four things enter the momentum account: two momentum fluxes and two pressure forces (the bed is frictionless,
the roller stays inside).""",
    change=r"""…the observer ran upstream at the bore speed: a moving bore becomes this stationary jump (Galilean frames, Ch. 3).""")
remind([
    ("control-volume momentum balance", r"momentum outflow − inflow = net force, the steady form of $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ *(4.17)*; a face carries the momentum flux ρu(u·n) (Ch. 4 §4.4, P114)."),
    ("quadratic formula, np.roots, the physical root", r"$r^2+r+c=0$ has $r=\tfrac12(-1\pm\sqrt{1-4c})$; the roots multiply to c (Vieta, Ch. 2 P71), so for c < 0 exactly one is positive. `np.roots([1, 1, c])` returns both."),
])
D("D25", ref="7.81")
note("N83", "The momentum balance in the book's form", "", equation=EQ["7.80"], ref="7.80")
D("D26")                                                   # the energy line after (7.81) has no number
note("N84", "Only upward jumps", r"""
The mechanical energy per unit mass of a surface particle drops across the jump,
$E_2-E_1=-(H_2-H_1)\frac{g(H_2-H_1)^2}{4H_1H_2}<0$ whenever H₂ > H₁; a "jump down" (Fr₁ < 1) would create mechanical energy
— the second law forbids it. The lost energy goes into the turbulent roller and ends as heat.
""")
nb.worked_example("a spillway jump", r"""
H₁ = 0.1 m, Fr₁ = 3.

1. $u_1=\mathrm{Fr}_1\sqrt{gH_1}=3\times0.990=2.97$ m/s, $Q=u_1H_1=0.297$ m²/s.
2. $H_2/H_1=\tfrac12(-1+\sqrt{1+72})=\tfrac12(-1+8.544)=3.772$ → H₂ = 0.377 m.
3. $u_2=Q/H_2=0.788$ m/s, $\mathrm{Fr}_2=0.788/\sqrt{9.81\times0.377}=0.41$ (subcritical).
4. Head loss $\frac{(H_2-H_1)^3}{4H_1H_2}=\frac{0.2772^3}{4\times0.1\times0.377}=0.141$ m — the energy of 14 cm of fall,
   ≈ ρgQ × 0.141 = 412 W per metre of width turned into heat.
""")
nb.code(r"""
j = ch07.hydraulic_jump(0.1, Fr1=3.0)                      # (7.81) and the energy loss
print({n: (round(float(v), 4) if not isinstance(v, bool) else v) for n, v in j.items()})   # every result, rounded
print("momentum residual of (7.80):", ch07.jump_momentum_residual(0.1, j["H2"], j["Q"]), "m³/s²")   # ≈ 0
st = ch07.jump_state(0.1, 3.0)                             # the momentum budget per metre of width
print(f"momentum in {st['mom_in']:.1f}, out {st['mom_out']:.1f} N/m; pressure force in {st['p_in']:.2f}, "   # the four terms of the budget …
      f"out {st['p_out']:.1f} N/m; balance {st['residual']:.1e} N/m; power lost {st['power_loss']:.0f} W/m")   # … their balance and the dissipated power
jd = ch07.hydraulic_jump(1.0, Fr1=0.7)                     # a subcritical 'jump'
print(f"Fr₁ = 0.7: H₂/H₁ = {jd['ratio']:.3f}, E₂ − E₁ = {jd['dE']:+.3f} J/kg, allowed = {jd['allowed']}")   # energy would increase: forbidden
print("bore 2 m deep into still water 1 m deep:", ch07.jump_state(1.0, 3.0**0.5, frame="moving")["bore_speed"], "m/s;",   # the moving-bore frame …
      "Ch. 4:", ch04.bore_speed(1.0, 2.0, g=G), "m/s")    # Fr₁ = √3 ⇔ H₂/H₁ = 2; ch04.bore_speed(h_still, h_behind)
""", explain=r"""
1. H₂ = 0.3772 m, u₂ = 0.788 m/s, Fr₂ = 0.410, $E_2-E_1=-1.385$ J/kg, head loss 0.1412 m.
2. The momentum budget balances: 882.9 − 234.1 N/m of momentum flux against 697.9 − 49.05 N/m of pressure force (to
   $10^{-13}$); 412 W per metre of width is dissipated.
3. For Fr₁ = 0.7 the formula still gives a number (H₂/H₁ = 0.609) but the energy *increases*: `allowed = False`.
4. A bore of depth 2 m running into still water 1 m deep (H₂/H₁ = 2 ⇔ Fr₁ = √3) moves at 5.42 m/s — the same as Ch. 4's
   bore speed (`ch04.bore_speed(h_in, h_out)`, still depth first).
""")
nb.check_agree(r"""
r = np.roots([1.0, 1.0, -2*3.0**2])                        # D25 step 9: r² + r − 2Fr₁² = 0 with Fr₁ = 3
print("both roots:", r)                                    # 3.772 and −4.772
r_pos = r[r > 0][0]                                        # physics keeps the positive depth ratio
assert np.allclose(r_pos, ch07.hydraulic_jump(0.1, Fr1=3.0)["ratio"], rtol=1e-12)   # = (7.81)
""")
nb.md(r"**What does this show?** The quadratic r² + r − 2Fr₁² = 0 of D25 step 9 (Fr₁ = 3) has the roots 3.772 and −4.772. A depth ratio cannot be negative, so the jump keeps 3.772 — exactly the ratio " + E("7.81") + r" that `hydraulic_jump` returns.")
nb.figure(r"""
Fr = np.linspace(0.3, 6, 300)                              # upstream Froude numbers [–]
res = [ch07.hydraulic_jump(1.0, Fr1=f) for f in Fr]        # H₁ = 1 m
ratio = np.array([q["ratio"] for q in res]); loss = np.array([q["head_loss"] for q in res])   # H₂/H₁ and loss/H₁
fig, ax = plt.subplots(figsize=(7, 3.4)); ax2 = ax.twinx()  # two y-axes
ax.plot(Fr, ratio, color=COLORS["ink"], lw=2, label="H₂/H₁ (7.81)")   # conjugate depth ratio
ax2.plot(Fr, loss, color=COLORS["rose"], lw=2, label="head loss / H₁")   # energy lost, in metres of fall per metre of H₁
ax.axvspan(0.3, 1.0, color=COLORS["muted"], alpha=0.15); ax.text(0.35, 6.3, "would create\nenergy (forbidden)", fontsize=7)   # the forbidden band
for f, lab in ((1.3, "undular"), (3.0, "our example"), (5.0, "spillway")):   # three marked cases
    q = ch07.hydraulic_jump(1.0, Fr1=f); ax.plot(f, q["ratio"], "o", color=COLORS["orange"]); ax.text(f + 0.1, q["ratio"] - 0.5, lab, fontsize=7)
ax.set_xlabel("Fr₁ [–]"); ax.set_ylabel("H₂/H₁ [–]"); ax2.set_ylabel("head loss / H₁ [–]", color=COLORS["rose"])   # axes
ax.set_ylim(0, 8); ax2.set_ylim(-0.5, 8)   # ranges
ax.set_title("Momentum sets the height, energy sets the direction", fontsize=10)   # the message
savefig(fig, "ch07", "c11_jump_curves"); plt.show()   # save, then draw
""", see=r"""The depth ratio (black) grows roughly like √2·Fr₁; the head loss (rose, right axis) rises steeply; the grey band
Fr₁ < 1 is forbidden.""",
    read=r"""At Fr₁ = 1 nothing happens (H₂ = H₁, no loss); below 1 the formula still gives a number but the loss would be
negative.""",
    change=r"""…H₁ were doubled at the same Fr₁: every ratio stays, the head loss in metres doubles (it scales with H₁).""")
note("N85", "Dispersion against steepening", r"""
In a dispersive medium the different wavelengths separate, which works against steepening; the two can balance and give
waves of permanent form.
""")
note("N93", "The Korteweg–de Vries equation", r"""
For long waves that are only weakly dispersive (λ several times H, but not so long that nonlinearity alone rules),
Korteweg and de Vries (1895) showed (7.87), with $c_0=\sqrt{gH}$: term by term,
shallow-water propagation, nonlinear steepening (note N81, the ³⁄₂ coefficient) and the first correction for the depth
not being quite shallow (weak dispersion).
""", equation=EQ["7.87"], ref="7.87")
note("N94", "Its linear speed", r"""
Without the nonlinear term, η = a cos(kx − ωt) gives $c=c_0(1-\tfrac16k^2H^2)$ — the first two Taylor terms of
$c=\sqrt{\frac gk\tanh kH}$ *(7.29)*. The error grows like (kH)⁴:
""")
nb.figure(r"""
kH = np.logspace(-2, 0, 30)                                # depth ratios [–] (H = 1 m)
err = np.abs(ch07.kdv_linear_phase_speed(kH, 1.0)/ch07.phase_speed(kH, 1.0) - 1)   # |c_KdV/c − 1|
fig, ax = plt.subplots(figsize=(5, 3.4))   # one panel
ax.loglog(kH, err, "o-", color=COLORS["ink"], ms=3, label="|c_KdV/c − 1|")   # the error of KdV's speed
ax.loglog(kH, 19/360*kH**4, ":", color=COLORS["muted"], label="(19/360)(kH)⁴")   # the next Taylor term
ax.set_xlabel("kH [–]"); ax.set_ylabel("relative error [–]"); ax.legend(fontsize=7)   # axes and legend
ax.set_title("KdV's dispersion is (7.29) to order (kH)²", fontsize=10)   # the message
savefig(fig, "ch07", "c11_kdv_speed"); plt.show()   # save, then draw
print(f"fitted slope: {observed_order(kH[:20], err[:20]):.2f}")   # expect 4
""", see=r"""A straight line of slope 4 (fitted slope printed below the figure).""",
    read=r"""KdV reproduces the true phase speed to about 5 % of (kH)⁴: 0.05 % at kH = 0.3, 5 % at kH = 1.""",
    change=r"""…the dispersive term were dropped: the error would be the full $(kH)^2/6$ — slope 2.""")
note("N95", "Which fate wins?", r"""
The ratio of the nonlinear to the dispersive term of (7.87) is $\frac\eta H\frac{\partial\eta}{\partial x}\Big/H^2\frac{\partial^3\eta}
{\partial x^3}\sim\frac{a\lambda^2}{H^3}$ (the Ursell number; the book's text says (7.88) but means the terms of (7.87)).
When it is large the front steepens into a jump; when it is of order one or smaller, dispersion can balance it. Where
exactly the crossover sits depends on how the slopes are estimated: with sinusoidal slopes (∂/∂x → k) the two KdV terms
are equal when $\frac{9a}{k^2H^3}=1$, i.e. at $a\lambda^2/H^3=4\pi^2/9\approx4.4$ (printed below).
""")
nb.code(r"""
for a, lam in ((0.1, 10.0), (0.2, 10.0), (0.05, 30.0)):   # amplitude and wavelength on H = 1 m [m]
    print(f"a = {a} m, λ = {lam:.0f} m: Ursell number aλ²/H³ = {ch07.ursell_number(a, lam, 1.0):.0f}")
print(f"equal KdV terms (sinusoidal slopes) at aλ²/H³ = 4π²/9 = {4*np.pi**2/9:.2f}")   # 9a/(k²H³) = 1 with k = 2π/λ
""", explain=r"""
1. The Ursell number $a\lambda^2/H^3$ for three waves on H = 1 m: 10, 20 and 45 — each far above the crossover, so these
   fronts steepen before dispersion can stop them.
2. The last line is the crossover estimate: the nonlinear term $\frac32c_0\frac aH\,ka$ equals the dispersive term
   $\frac16c_0H^2k^3a$ when $9a/(k^2H^3)=1$.
""")
remind([
    ("solitary wave (7.88), sech² and its derivatives", r"$\mathrm{sech}=1/\cosh$ (primer P168); $\frac{d}{dx}\tanh x=\mathrm{sech}^2x$ and $\frac{d}{dx}\mathrm{sech}\,x=-\mathrm{sech}\,x\tanh x$ — all the calculus the soliton check needs."),
])
note("N96, N97", "Cnoidal and solitary waves", r"""
The periodic permanent-form solutions of (7.87) are cnoidal waves (Jacobi's cn functions, drawn with `ch07.cnoidal_wave`).
The single-hump solution is the solitary wave (7.88): taller is faster (checked by substitution below; the
derivatives used are $\tanh'=\mathrm{sech}^2$ and $\mathrm{sech}'=-\mathrm{sech}\tanh$).
""", equation=EQ["7.88"], ref="7.88")
nb.code(r"""
print(ch07.kdv_residual_sympy())                            # (7.88) put into (7.87): the residual is exactly 0
print(f"a = 0.2 m on H = 1 m: c = c₀(1 + a/2H) = {np.sqrt(G*1.0)*(1 + 0.2/2):.3f} m/s")   # (7.88)
""")
nb.figure(r"""
x = np.linspace(-40, 40, 800)                              # positions [m]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3))         # cnoidal | solitary
cn = ch07.cnoidal_wave(x, 0.0, H=1.0, height=0.2, m=0.99)  # a cnoidal train, crest-to-trough 0.2 m
a1.plot(x, cn["eta"], color=COLORS["blue"], lw=2); a1.axhline(0, color=COLORS["muted"], lw=0.5)   # the cnoidal surface and the still level
a1.set_title(f"cnoidal wave: λ = {cn['wavelength']:.1f} m, c = {cn['c']:.2f} m/s", fontsize=9)   # its wavelength and speed
a2.plot(x, ch07.solitary_wave(x, 0.0, 0.2, 1.0), color=COLORS["blue"], lw=2); a2.axhline(0, color=COLORS["muted"], lw=0.5)   # (7.88)
a2.set_title(f"solitary wave (7.88): a = 0.2 m, c = {ch07.solitary_wave_speed(0.2, 1.0):.2f} m/s", fontsize=9)   # its speed c₀(1 + a/2H)
for axp in (a1, a2):   # both panels:
    axp.set_xlabel("x [m]"); axp.set_ylim(-0.1, 0.25)   # same axes
a1.set_ylabel("η [m] (H = 1 m)")   # unit
savefig(fig, "ch07", "c11_cnoidal_solitary"); plt.show()   # save, then draw
""", see=r"""Left, a cnoidal train with sharp crests and long flat troughs; right, a single hump (our remake of Fig. 7.23
`N189`).""",
    read=r"""Steepening (the $\eta\,\partial\eta/\partial x$ term) and spreading (the $\partial^3\eta/\partial x^3$ term) cancel
exactly for these shapes, so they travel unchanged.""",
    change=r"""…a were doubled: the hump narrows by √2 and speeds up to $c_0(1+a/H)$.""")
nb.animation(r"""
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 4.4), layout="none")   # (a) no dispersion | (b) KdV; manual spacing below
xs_ = np.linspace(0, 200, 2001); e0 = 0.2*np.exp(-((xs_ - 50)/10)**2)   # the hump of N81
sw = ch07.simple_wave_evolve(e0, xs_, 0.0, 1.0)            # breaking time
(l1,) = a1.plot(xs_, e0, color=COLORS["blue"], lw=2); a1.set_xlim(20, 110); a1.set_ylim(-0.03, 0.25)   # panel (a): the hump
a1.set_ylabel("η [m]")   # unit
N = 1024 if not FAST else 512                              # KdV grid (spectral)
xk = np.linspace(-50, 350, N, endpoint=False)              # periodic domain [m]
ek0 = 0.1*np.exp(-(xk/10)**2)                              # a 10 cm Gaussian hump, 10 m wide, H = 1 m
tk = np.linspace(0, 100, 16 if not FAST else 9)           # output times [s]
kd = ch07.kdv_solve(ek0, xk, tk, 1.0, cache=str(ROOT/"outputs"/"ch07"/f"kdv_run_{N}.npz"))   # (7.87), cached
(l2,) = a2.plot(xk, ek0, color=COLORS["blue"], lw=1.5); a2.set_xlim(-50, 350); a2.set_ylim(-0.04, 0.2)   # panel (b): the KdV hump
a2.set_xlabel("x [m]"); a2.set_ylabel("η [m]")   # units


def update(i):   # draw frame i
    t1 = i/(len(tk) - 1)*sw["t_break"]                     # (a) runs up to the breaking time
    l1.set_xdata(ch07.simple_wave_evolve(e0, xs_, t1, 1.0)["x_points"])   # each surface point moved along its characteristic
    a1.set_title(f"(a) no dispersion: steepening, t = {t1:.1f} s (breaks at {sw['t_break']:.1f} s)", fontsize=9)   # time in the title
    l2.set_ydata(kd["eta"][i])                              # (b) the KdV solution at tk[i]
    a2.set_title(f"(b) KdV (7.87): the hump sorts itself into solitons, t = {tk[i]:.0f} s", fontsize=9)   # time in the title
    return []   # nothing to blit


fig.subplots_adjust(left=0.1, right=0.97, bottom=0.1, top=0.92, hspace=0.5)   # margins and room for the two titles
show_animation(animate(update, frames=len(tk), fig=fig, interval=200), player="frames")   # a frame player: step with ◀ ▶
print(f"KdV invariants: mass drift {abs(kd['mass'][-1]/kd['mass'][0] - 1):.1e}, "   # conserved quantities of (7.87) …
      f"energy drift {abs(kd['energy'][-1]/kd['energy'][0] - 1):.1e} (relative); {kd['steps']} steps, cached = {kd['cached']}")   # … should not drift
""")
see_read_change(r"""Step through the frames (◀ ▶): (a) without dispersion the hump's front becomes vertical at $t_b$; (b) with the KdV
dispersive term the same kind of hump splits into a train of solitary waves, the tallest (fastest) in front, followed by a
small dispersive tail — the `A5` player.""",
                r"""The printed invariants (mass and the KdV energy) stay constant to many digits: the solver is accurate.
Taller solitons are faster, $c=c_0\big(1+\frac a{2H}\big)$ *(7.88)*, so the train spreads out in order of height.""",
                r"""…the hump were wider (smaller Ursell number): fewer solitons, eventually just one.""")
explainer("hydraulic_jump", "How high does the water jump — and why only upward?",
          r"""The Fr₁ slider moves the downstream depth, the momentum bars and the energy bar together; pushing Fr₁ below 1
turns the energy change positive — the second law choosing the root is felt, not just stated.""", "",
          ["Preset *spillway*: read H₂ and the head loss.",
           "Drag Fr₁ to 0.7: the status turns red — the energy bar is positive.",
           "Switch to the moving-bore frame: the same depths, a different observer.",
           "Choose the solitary-wave fate and raise a/H: watch the Ursell number grow and the front steepen."])
whatif(r"""
…the layer were a denser fluid under a lighter one (cold air under warm air, salty water under fresh)? g becomes the
reduced gravity g′ of C14 and the same jump appears as an internal hydraulic jump in the lee of mountains (Ch. 13).
""")

# ---- C12 ------------------------------------------------------------------------------------------------------------
core("C12", r"Stokes waves and Stokes drift", r"""
If every orbit is a closed circle, why does floating debris travel with the waves?
""", eqs=("7.85",))
problem(r"""
A swimmer floating beyond the surf slowly drifts toward the beach although "waves only move up and down"; oil, plastic and
larvae at sea are carried downwind by the waves themselves. We want that mean drift: its size, how fast it dies with
depth, and why a current meter fixed in place measures no mean current at all.
""")
note("N86, N87", "Stokes waves", r"""
In 1847 Stokes found periodic waves of finite amplitude in deep water, (7.82), with
$c=\sqrt{\frac gk(1+k^2a^2+\ldots)}$ *(7.83)*: the harmonics sharpen the crests and flatten the troughs, and the speed grows
with amplitude (4.4 % at ka = 0.3). We state them; the code cell below re-derives the coefficients in sympy.
""", equation=EQ["7.82"], ref="7.82")
nb.md(r"""
**How the coefficients are found — and a trap** (verification note O4). Write, with ε = ka and k = g = 1,
$\eta=\varepsilon\cos\theta+\alpha\varepsilon^2\cos2\theta+\delta\varepsilon^3\cos3\theta$,
$\phi=\varepsilon\omega(1+\beta\varepsilon^2)e^{z}\sin\theta$ and $\omega^2=1+\gamma\varepsilon^2$, put them into the exact
kinematic condition (7.16) and the full Bernoulli equation with p = 0, both on z = η, and collect powers of ε and
harmonics:

| order | condition | gives |
|---|---|---|
| ε² | kinematic, sin 2θ | α = ½ — the second harmonic of (7.82) |
| ε³ | kinematic, sin θ and sin 3θ | β = −5/8, δ = 3/8 |
| ε³ | dynamic, cos θ | γ = 1 — the speed correction of (7.83) |

With β = −5/8 and $\omega\approx1+\varepsilon^2/2$, the potential's first harmonic is $\varepsilon(1-\varepsilon^2/8)$: the
amplitude of φ is corrected at second order. **If instead the potential is kept at aω/k (β = 0) and the kinematic
condition only to first order** — a tempting shortcut — the ε³ balance gives γ = 3/8, not 1. The
consistent third-order solution (the coded one, and an independent expansion in the tests) gives γ = 1, i.e. (7.83).
""")
remind([
    ("perturbation series checked in sympy", "`sp.series(expr, eps, 0, 4).removeO()` keeps powers of ε up to ε³, `.coeff(eps, n)` reads one order, and a product of sines is turned into harmonics before reading a coefficient (Ch. 4 P117)."),
])
nb.code(r"""
st = ch07.stokes_expansion_sympy()                        # the expansion above, done in sympy (cached)
print("consistent third order:", st["third_order"])        # alpha = 1/2, beta = −5/8, delta = 3/8, gamma = 1
print("truncated set-up (β frozen at 0):", st["exercise_literal"])   # alpha = 1/2 but gamma = 3/8: β was dropped
print(f"c(ka = 0.3)/c(ka → 0) = {ch07.stokes_wave_speed(1.0, 0.3)/ch07.phase_speed(1.0):.4f}")   # (7.83): 4.4 % faster
""", explain=r"""
1. `stokes_expansion_sympy` projects each condition on its harmonics order by order: α = ½ and δ = ³⁄₈ (the coefficients
   of (7.82)) and γ = 1 (the ka² of (7.83)) — plus β = −⁵⁄₈, the correction of φ's first harmonic.
2. The same machinery with φ frozen at aω/k gives γ = ³⁄₈: dropping a term of the same order as the one kept changes the
   answer.
3. At ka = 0.3 the Stokes wave is 4.4 % faster than the linear one.
""")
note("N88", "The steepest wave", r"""
The steepest possible Stokes wave has a sharp crest of angle 120° (Stokes 1880) and crest-to-trough height/wavelength
of about one-seventh (Schwartz and Fenton's 1982 review); the digits 0.1411 in `ch07.STOKES_LIMIT_STEEPNESS` are from the
computation of Dyachenko, Lushnikov and Korotkevich (2016). Half of it, measured from the mean level, is the largest
amplitude, a ≈ 0.071λ (printed below); beyond it the crest spills (white caps). The truncated series
(7.82) cannot reach that corner; we show profiles only up to ka = 0.3.
""")
nb.code(r"""
Hmax = ch07.STOKES_LIMIT_STEEPNESS                         # crest-to-trough height / wavelength of the steepest wave [–]
print(f"H/λ = {Hmax:.4f} ≈ 1/{1/Hmax:.2f};  largest amplitude a = H/2 ≈ {Hmax/2:.3f} λ;  ka = π·H/λ ≈ {np.pi*Hmax:.2f}")   # a from the mean level
""")
nb.figure(r"""
x = np.linspace(0, 4*np.pi, 800)                           # two wavelengths, k = 1 rad/m [m]
fig, ax = plt.subplots(figsize=(7, 2.8))                   # one panel
ax.plot(x, ch07.stokes_wave_profile(x, 0.0, 0.3, 1.0), color=COLORS["ink"], lw=2, label="Stokes (7.82), ka = 0.3")   # three harmonics
ax.plot(x, 0.3*np.cos(x), "--", color=COLORS["blue"], lw=1.5, label="linear a cos kx")   # first harmonic only
ax.axhline(0, color=COLORS["muted"], lw=0.5); ax.set_xlabel("x [m]"); ax.set_ylabel("η [m]"); ax.legend(fontsize=7)   # still level and axes
ax.set_title("Finite amplitude: sharper crests, flatter troughs", fontsize=10)   # the message
savefig(fig, "ch07", "c12_stokes_profile"); plt.show()   # save, then draw
""", see=r"""The Stokes profile (black) at ka = 0.3 against the linear cosine (blue dashed) — our remake of Fig. 7.21 `N187`.""",
    read=r"""The second harmonic $\tfrac12ka^2\cos2k(x-ct)$ adds at the crests and subtracts at the troughs: crests 0.355 m
high, troughs only 0.265 m deep.""",
    change=r"""…ka = 0.1: the two curves almost coincide (the correction is ∝ ka).""")
note("N89", "Stokes drift", r"""
At finite amplitude a particle is a little further forward at the top of its loop, where the forward velocity is larger,
than it is backward at the bottom, where it is smaller: the loop does not close and the particle creeps forward (the
book's Fig. 7.22). **Lagrangian vs Eulerian mean**: the average of the velocity *following a particle* ($\bar u_L$) versus
the average *at a fixed point* (ū); they differ whenever the particle samples a non-uniform field (Ch. 3 §3.2).
""")
note("N90", "The velocity a parcel really feels", r"""
Keep the first Taylor terms of u about the mean position (Ch. 3 P98): (7.84a), and the same for w, $\frac{dz_p}{dt}=
w(x_0,z_0,t)+\xi\big(\frac{\partial w}{\partial x}\big)_0+\zeta\big(\frac{\partial w}{\partial z}\big)_0+\ldots$ *(7.84b)*.
""", equation=EQ["7.84a"], ref="7.84a")
D("D27", ref="7.85")
remind([
    ("fundamental theorem of calculus", r"$f(z)=f(z_0)+\int_{z_0}^z f'(s)\,ds$ (Ch. 2 P84); with $\partial u/\partial z=\partial w/\partial x$ (no vorticity) it rebuilds u from the bottom up in the note below."),
])
note("N91, N92", "Any depth, and the Eulerian mean", r"""
At any depth the same moves with cosh and sinh give (7.86); the vertical drift is zero.
At a point that stays under water the mean velocity is exactly zero: irrotationality gives
$u=u\vert_{z=-H}+\int_{-H}^z\frac{\partial w}{\partial x}dz$, whose wavelength average vanishes for a periodic wave. All the
transport is Lagrangian.
""", equation=EQ["7.86"], ref="7.86")
nb.code(r"""
print(f"Eulerian mean at z = −0.5 m: {ch07.eulerian_mean_u(-0.5, 0.05, 1.0, 5.0):.1e} m/s;",
      f"Stokes drift there: {ch07.stokes_drift(-0.5, 0.05, 1.0, 5.0):.2e} m/s")   # a = 5 cm, k = 1 rad/m, H = 5 m
""")
nb.worked_example("drift under a 1 m swell", r"""
a = 1 m, T = 8 s, deep water.

1. ω = 0.785 rad/s, k = ω²/g = 0.0629 rad/m.
2. Surface drift $a^2\omega k=1\times0.785\times0.0629=0.0494$ m/s ≈ 4.9 cm/s ≈ 4 km per day.
3. It decays like $e^{2kz_0}$ — twice as fast as the orbits: halved every $\ln2/2k=5.5$ m.
4. Compare the orbital speed aω = 0.785 m/s: the drift is ka = 6 % of it — a second-order effect.
""")
nb.code(r"""
k, a = 1.0, 0.05                                           # deep water, ka = 0.05 [rad/m, m]
for z0 in (0.0, -0.5, -1.0):                               # three mean depths [m]
    an = ch07.stokes_drift(z0, a, k)                       # (7.85) [m/s]
    nu_ = ch07.stokes_drift_numeric(z0, a, k, periods=20 if not FAST else 10)   # exact path lines, drift per time
    print(f"z₀ = {z0:4.1f} m: (7.85) {an:.3e} m/s, exact path lines {nu_:.3e} m/s ({nu_/an - 1:+.1%})")
print(f"1 m, 8 s swell at the surface: {ch07.stokes_drift(0.0, 1.0, (2*np.pi/8)**2/G):.4f} m/s")
""", explain=r"""
1. (7.85) gives 7.83, 2.88 and 1.06 mm/s at 0, 0.5 and 1 m depth. The drift measured from exact path lines
   (`solve_ivp`, DOP853, rtol = 1e-10: the drift is a few % of the orbital speed, so the integrator must be far more
   accurate than that, Ch. 3 P94) agrees within a few per cent — the difference is itself O(ka), the order of the terms
   the derivation dropped.
2. The 1 m swell of the tiny example drifts at 0.0494 m/s at the surface.
""")
nb.check_agree(r"""
a, k, z0 = 0.02, 1.0, -0.5                                 # a gentler wave, ka = 0.02 [m, rad/m, m]
T = 2*np.pi/ch07.omega_gravity(k)                          # period [s]
tt = np.arange(0, 21)*T                                    # 21 whole-period instants [s]
p = ch07.particle_path(0.0, z0, tt, a, k, model="exact")   # the exact path line (7.33), sampled once per period
slope = np.polyfit(tt, p["x"], 1)[0]                       # period-averaged forward speed [m/s]
print(f"measured drift {slope:.4e} m/s, (7.85) {ch07.stokes_drift(z0, a, k):.4e} m/s")
assert np.allclose(slope, ch07.stokes_drift(z0, a, k), rtol=3e-2)   # agree to O(ka) = 2 %: the drift is the mean of the exact motion
""")
nb.md(r"""
The particles really do creep forward at $\bar u_L=a^2\omega ke^{2kz_0}$ *(7.85)* — the analytic drift is the average of the
exact motion, not an extra assumption (at ka = 0.02 the two differ by about 1 %, as expected from the dropped O(ka)
terms).
""")
nb.figure(r"""
k, a = 1.0, 0.05                                            # ka = 0.05 [rad/m, m]
z0 = np.linspace(-3, 0, 200)                                # mean depths [m]
fig, ax = plt.subplots(figsize=(6, 3.6))                    # drift profile
ax.plot(ch07.stokes_drift(z0, a, k)*1e3, z0, color=COLORS["accent"], lw=2, label="deep (7.85)")   # deep water [mm/s]
zH = z0[z0 >= -1.0]                                         # only above the bottom of the H = 1 m case [m]
ax.plot(ch07.stokes_drift(zH, a, k, 1.0)*1e3, zH, "--", color=COLORS["accent"], lw=1.5, label="kH = 1 (7.86)")   # finite depth H = 1 m
ax.axhline(-1.0, color="#c8b99a", lw=2); ax.text(8, -0.95, "bottom for kH = 1", fontsize=7)   # where the H = 1 m case ends
zd = np.linspace(-2.5, 0, 6)                                # six depths for the numerical drift [m]
ax.plot([ch07.stokes_drift_numeric(z, a, k, periods=20 if not FAST else 10)*1e3 for z in zd], zd, "o",   # drift measured on exact path lines
        color=COLORS["accent"], label="exact path lines")   # … as dots
ax.plot(ch07.stokes_drift(0.0, a, k)*1e3*np.exp(k*z0), z0, color=COLORS["teal"], lw=1.5,   # orbit size ∝ e^{kz₀}, scaled
        label=r"orbit radius $\propto e^{kz_0}$ (scaled)")           # decay of the orbits, for comparison
ax.set_xlabel("ū_L [mm/s]"); ax.set_ylabel("mean depth z₀ [m]"); ax.legend(fontsize=7)   # axes and legend
ax.set_title("The drift lives close to the surface", fontsize=10)   # the message
savefig(fig, "ch07", "c12_drift_profile"); plt.show()   # save, then draw
""", see=r"""The drift (purple) falls twice as fast with depth as the orbit size (teal, scaled to the same surface value); the
numerical dots from exact path lines sit on the deep-water curve; the dashed curve is the finite-depth (7.86) with kH = 1.""",
    read=r"""$e^{2kz_0}$ vs $e^{kz_0}$: one wavelength down ($z_0=-\lambda$) the orbit keeps $e^{-2\pi}\approx0.2$ % and the
drift $e^{-4\pi}\approx4\times10^{-6}$.""",
    change=r"""…the depth were finite (kH = 1): the drift does not vanish at the bottom — cosh 2k(z₀ + H) ≥ 1 (dashed).""")
nb.animation(r"""
k, a = 1.0, 0.25                                            # a steep deep-water wave, ka = 0.25 [rad/m, m]
w = ch07.omega_gravity(k); T = 2*np.pi/w                    # ω [rad/s], period [s]
z0s = np.linspace(-2.0, 0.0, 9)                             # nine particles on a vertical line at x = 0 [m]
nf = 50 if not FAST else 25                                 # frames over five periods
tt = np.linspace(0, 5*T, nf)                                # frame times [s]
dl = ch07.dyed_line(z0s, tt, a, k)                          # exact path lines of all particles, all times (one solve_ivp)
fig, ax = plt.subplots(figsize=(7, 3.4))
xs = np.linspace(-np.pi, 2*np.pi, 300)                      # surface points [m]
(surf,) = ax.plot(xs, a*np.cos(k*xs), color=COLORS["blue"], lw=2)
trails = [ax.plot([], [], color=COLORS["teal"], lw=0.8)[0] for _ in z0s]   # each particle's open loop so far
(line,) = ax.plot(dl["x"][0], dl["z"][0], color=COLORS["rose"], lw=2)      # the dyed line
dots = ax.scatter(dl["x"][0], dl["z"][0], s=14, color=COLORS["teal"], zorder=3)
ax.set_xlim(-np.pi, 2*np.pi); ax.set_ylim(-2.4, 0.5); ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")


def update(i):   # draw frame i
    surf.set_ydata(a*np.cos(k*xs - w*tt[i]))                 # the linear surface (the particles follow (7.33))
    for j, tr in enumerate(trails):                          # the loops traced up to now
        tr.set_data(dl["x"][:i + 1, j], dl["z"][:i + 1, j])
    line.set_data(dl["x"][i], dl["z"][i]); dots.set_offsets(np.c_[dl["x"][i], dl["z"][i]])
    ax.set_title(f"exact path lines, ka = 0.25: after {tt[i]/T:.1f} periods the dye line leans forward", fontsize=9)
    return []   # nothing to blit


show_animation(animate(update, frames=nf, fig=fig, interval=60))
""")
see_read_change(r"""Nine particles released on a vertical line under a steep wave: each goes round an open loop (teal trail) and the
rose dye line through them tilts forward period by period — our version of Fig. 7.22 `N188`, the `A4` animation.""",
                r"""The top of the dye line runs ahead; after a few periods the line is curved like $e^{2kz_0}$ — the drift profile
of the previous figure.""",
                r"""…ka were halved: the drift falls by four (∝ a²) while the loops shrink only by two.""")
explainer("particle_orbits", "What does the water under a wave actually do?",
          r"""Tracer particles loop under the moving surface; drag the depth to squash circles into ellipses, then switch the
model from linear to exact and watch the loops open and the dyed column lean forward. Motion and one model toggle are the
whole lesson.""", "",
          ["Start in *deep water*, click a particle at 1 m depth and read its A, B in the inspector.",
           "Drag kH to 0.3: the circles become flat ellipses; the bottom particle only slides.",
           "Choose *steep deep wave (exact)*: after five periods read the drift in the end-of-run card and compare with $a^2\\omega k$.",
           "Open the Derivation (D27) at 'average the two terms'."])
confusion(r"""
"the mean current is zero" and "the water drifts" are both true. A fixed meter averages over whatever water passes it
(zero); a float averages following one parcel ($a^2\omega ke^{2kz_0}$). Ch. 13 adds the Earth's rotation, and the Stokes
drift then drives the Stokes–Coriolis force and Langmuir circulation in the ocean's surface layer.
""")
whatif(r"""
…the waves broke? Breaking transfers momentum to the water as a real (Eulerian) current — the longshore currents and
undertow of beaches — which the irrotational theory here cannot describe.
""")

# =====================================================================================================================
# A.7 §7.7 — R11, C13 (N98–N108, N190–N192, D28, P176), C14 (N109–N131, N193, N194, D29–D31, E8)
# =====================================================================================================================
nb.section("7.7", "Waves on a Density Interface", intro=r"""
**What is this section about?** Replace the air above the water by a slightly lighter liquid — fresh water over salt
water, warm water over cold. Gravity still restores the interface, but only through the small density difference: the
waves become slow and tall. With a free surface on top, a layer of fluid can oscillate in two ways — the surface way and
the interface way — the barotropic and baroclinic modes that layered ocean and atmosphere models are built on.
""")
nb.recap("R11", "Laplace in each fluid", r"""
Each fluid is irrotational and incompressible on its own side, so $\nabla^2\phi_1=0,\ \nabla^2\phi_2=0$ *(7.90)*: one
Laplace equation per layer (Ch. 6 §6.2, recap R02).
""", where="Ch. 6 §6.2")
core("C13", r"Interfacial waves", r"""
How does a wave on a weak density step differ from a wave on the sea surface?
""", eqs=("7.95",))
note("N98", "The set-up", r"""
A lighter fluid ρ₁ above a heavier one ρ₂, both infinitely deep, small slopes, no interfacial tension (the book's Fig.
7.24) — a fjord's fresh upper layer, a sun-warmed surface layer, oil on water.
""")
problem(r"""
In a fjord, a layer of river water floats on sea water; in summer a lake's warm surface layer floats on cold water. Their
interface carries waves you cannot see from a boat but a thermometer string can: slow, and metres tall. We want their
speed and why they are so different from surface waves.
""")
nb.md(r"""
#### The idea

Gravity pulls a raised piece of interface down with the weight of the *difference* ρ₂ − ρ₁, but the inertia to be moved is
both fluids, ρ₂ + ρ₁:

```
surface wave:    ω² = gk · (ρ − 0)/(ρ + 0)      = gk
interface wave:  ω² = gk · (ρ₂ − ρ₁)/(ρ₂ + ρ₁)  = ε² gk    ε ≈ 0.03 for Δρ/ρ = 0.002
```

**A small density step makes a slow wave — and the two fluids slide past each other.**
""")
P("P176", "complex amplitudes", r"""
Write a real wave as the real part of a complex one, $\zeta=\mathrm{Re}\{a\,e^{i(kx-\omega t)}\}$, and drop the Re during
linear algebra: ∂/∂x becomes multiplication by ik and ∂/∂t by −iω, and a complex amplitude $b=\lvert b\rvert e^{i\varphi}$
carries both size and phase shift (Ch. 1 P45 Euler's formula; Ch. 6 P153). ⚠️ Take real parts before multiplying two
fields (energy, flux: C16).
""", code=r"""
import numpy as np                                # numbers
k, w, t = 1.0, 2.0, 0.3                           # wavenumber [rad/m], frequency [rad/s], time [s]
x = np.linspace(0, 1, 3)                          # three positions [m]
zeta = np.real(1.0*np.exp(1j*(k*x - w*t)))        # a cos(kx − ωt) with a = 1 m
dzdt = np.real(-1j*w*1.0*np.exp(1j*(k*x - w*t)))  # ∂/∂t ↦ −iω
print(np.allclose(dzdt, w*np.sin(k*x - w*t)))     # True: the same as differentiating the cosine
""")
note("N99, N100", "Complex notation from here on", r"""
The book writes the interface as (7.89) with Re{} dropped; complex constants A, B, C, b carry phase shifts. `W.real_field(amp,
phase)` returns the real part when we plot.
""", equation=EQ["7.89"], ref="7.89")
note("N101–N104", "The two-fluid problem", r"""
Laplace in each fluid (7.90); $\phi_1\to0\ \text{as}\ z\to\infty$ *(7.91)*; $\phi_2\to0\ \text{as}\ z\to-\infty$ *(7.92)*;
both fluids move with the interface, $\frac{\partial\phi_1}{\partial z}=\frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}
{\partial t}\ \text{at}\ z=0$ *(7.93)*; and the pressure is continuous (no interfacial tension) — C02's conditions, now two-sided:
""", equation=EQ["7.94"], ref="7.94")
remind([
    ("Euler's formula and i² = −1", r"$e^{i\theta}=\cos\theta+i\sin\theta$, so $(-i\omega)(i\omega a/k)=\omega^2a/k$ (Ch. 1 P45)."),
    ("the complex plane in numpy", "`1j` is i; `np.real`, `np.abs`, `np.angle` read a complex number (Ch. 6 P153)."),
    ("np.linalg.solve with complex numbers", "`np.linalg.solve(A, b)` solves A·x = b; complex entries are fine (Ch. 2 P57)."),
])
D("D28", ref="7.95")
nb.worked_example("a thermocline wave", r"""
ρ₁ = 1000, ρ₂ = 1002 kg/m³, λ = 100 m.

1. $\varepsilon^2=2/2002=9.99\times10^{-4}$, ε = 0.0316.
2. A surface wave of this length: $\omega=\sqrt{gk}=0.785$ rad/s, T = 8.0 s.
3. The interfacial wave: ω = 0.0316 × 0.785 = 0.0248 rad/s, T = 253 s ≈ 4.2 min — 32 times slower.
4. The same energy as a 1 m surface wave needs $a=\sqrt{\rho/\Delta\rho}\times1\ \text{m}=\sqrt{500}=22$ m.
""")
P("P183", "warnings.catch_warnings(record=True)", r"""
A Python *warning* is a message a function raises without stopping ("this input makes no physical sense; returning NaN").
Inside `with warnings.catch_warnings(record=True) as caught:` every warning is collected into the list `caught` instead
of being printed, and `warnings.simplefilter("always")` makes sure none is suppressed as a repeat. We then print the
message ourselves, next to the result it explains.""", code=r"""
import warnings                                           # Python's warning system
with warnings.catch_warnings(record=True) as caught:      # collect warnings in the list `caught`
    warnings.simplefilter("always")                       # do not hide repeats
    warnings.warn("heavier fluid on top")                 # a function would raise this
print(len(caught), caught[0].message)                     # 1 heavier fluid on top
""")
nb.code(r"""
k = 2*np.pi/100                                            # λ = 100 m [rad/m]
print(f"ω interface = {ch07.interface_omega(k, 1000.0, 1002.0):.5f} rad/s, surface √(gk) = {ch07.omega_gravity(k):.4f} rad/s")   # (7.95) vs (7.45)
f = ch07.interface_fields(np.array([0.0, 25.0]), np.array([0.0, 0.0]), 0.0, 1.0, k, 1000.0, 1002.0)   # a = 1 m, at a crest and a node
print("u₁ (above) =", f["u1"], " u₂ (below) =", f["u2"], " sheet strength u₂ − u₁ =", f["gamma_sheet"], "m/s")   # opposite velocities
print({n: f"{np.max(np.abs(v)):.1e}" for n, v in ch07.interface_residuals(10.0, 3.0, 1.0, k, 1000.0, 1002.0).items()})   # residuals of (7.90)–(7.94)
print(ch07.interface_energy(1.0, k, 1000.0, 1002.0), ch07.interface_energy(1.0, k, 1000.0, 1002.0, method="quad"))   # closed form and dblquad
with warnings.catch_warnings(record=True) as caught:      # heavy fluid on top: the code refuses a real ω
    warnings.simplefilter("always")   # show every warning
    print("ρ₁ = 1002 over ρ₂ = 1000:", ch07.interface_omega(k, 1002.0, 1000.0), "—", str(caught[0].message)[:60], "…")   # NaN and the reason
""", explain=r"""
1. ω = 0.02481 rad/s against 0.7851 rad/s at the surface: 32 times slower.
2. Under the crest (x = 0) the upper fluid moves backward and the lower one forward at 0.0248 m/s: the interface is a
   vortex sheet of strength $u_2-u_1=0.0496$ m/s (note N107).
3. The residuals of (7.90)–(7.94) are at round-off.
4. $E_k=E_p=4.905$ J/m² and E = 9.81 J/m² for a = 1 m — a quarter of Δρga² each — both in closed form and by `dblquad`.
5. Heavy fluid on top: ω² < 0, so the code returns NaN with a warning — the Rayleigh–Taylor instability of Ch. 11.
""")
nb.check_agree(r"""
a, rho1, rho2 = 1.0, 1000.0, 1002.0                       # amplitude [m]; upper, lower density [kg/m³]
al, be = np.linalg.solve(np.array([[-k, 0], [0, k]], dtype=complex), np.array([-1j*a, -1j*a]))   # (7.93) per unit ω: −kA = kB = −iωa
print("A/ω =", al, " B/ω =", be, " (= ±ia/k)")            # so A = iωa/k, B = −iωa/k (D28 steps 4–5)
# (7.94) with ∂/∂t → −iω and A = ω·al, B = ω·be:  ω²(−iρ₁al + iρ₂be) = (ρ₂ − ρ₁) g a  — linear in ω²
w2_mine = (rho2 - rho1)*G*a/(-1j*rho1*al + 1j*rho2*be)     # ω² [1/s²] (a real number)
w_mine = np.sqrt(w2_mine.real)                             # ω [rad/s]
print(f"ω from our two solves: {w_mine:.6f} rad/s")
assert np.allclose(w_mine, ch07.interface_omega(k, rho1, rho2), rtol=1e-12)   # D28 step 9: ω = ε√(gk)
""")
note("N105, N106", "Interfacial energy", r"""
The kinetic and potential energies are equal, $E_k=E_p=\tfrac14(\rho_2-\rho_1)ga^2$ (the column swap of C06
now lifts ρ₂ and lowers ρ₁ — the book's Fig. 7.25 `N191` is that sketch with a second fluid, not redrawn), so the total is
(7.96): for the same energy an interfacial wave is $\sqrt{\rho_2/(\rho_2-\rho_1)}$ times taller than a surface wave (22×
for the thermocline above).
""", equation=EQ["7.96"], ref="7.96")
note("N107", "The interface is a vortex sheet", r"""
The horizontal velocities $u_1=\partial\phi_1/\partial x=-\omega ae^{-kz}e^{i(kx-\omega t)}$ and $u_2=\omega ae^{kz}e^{i(kx-
\omega t)}$ are opposite: the interface is a **vortex sheet** of strength $u_2-u_1=2\omega a\cos(kx-\omega t)$ (the jump of
Ch. 5's vortex sheets) — the seed of the Kelvin–Helmholtz instability (Ch. 11). A continuous stratification spreads this
vorticity through the fluid, so internal waves (C15) are not irrotational. Ch. 5's function (counterclockwise convention,
$\gamma=u_{below}-u_{above}$) gives the same number:
""")
nb.code(r"""
print(ch05.vortex_sheet_strength(u_above=f["u1"][0], u_below=f["u2"][0]), "=", f["gamma_sheet"][0], "m/s")   # at the crest
""")
nb.figure(r"""
k, a = 2*np.pi/100, 1.0                                    # λ = 100 m, a = 1 m
xs = np.linspace(0, 200, 400)                              # two wavelengths [m]
zs = np.array([-60, -45, -30, -15, 15, 30, 45, 60.0])      # heights of the arrow rows, away from the interface [m]
XX, ZZ = np.meshgrid(np.linspace(0, 200, 17), zs)          # arrow grid
f = ch07.interface_fields(XX, ZZ, 0.0, a, k, 1000.0, 1002.0)   # u, w in the fluid on each side (real parts)
fig, ax = plt.subplots(figsize=(7, 3.2))
ax.fill_between(xs, -70, 10*np.cos(k*xs), color=NAVY, alpha=0.25)   # lower (denser) fluid, interface ×10 for visibility
ax.plot(xs, 10*np.cos(k*xs), color=NAVY, lw=2, label=r"interface $\zeta$ (×10)")
up = ZZ > 0                                                 # arrows above (light) and below (navy)
ax.quiver(XX[up], ZZ[up], f["u"][up], f["w"][up], color=COLORS["muted"], scale=0.4, width=0.004)
ax.quiver(XX[~up], ZZ[~up], f["u"][~up], f["w"][~up], color=NAVY, scale=0.4, width=0.004)
for xc in (8, 100, 192): ax.text(xc, 64, "γ > 0", fontsize=7, ha="center", bbox=dict(facecolor="white", edgecolor="none", pad=1))   # sheet strength positive above crests
for xc in (50, 150): ax.text(xc, 64, "γ < 0", fontsize=7, ha="center", bbox=dict(facecolor="white", edgecolor="none", pad=1))       # and negative above troughs
ax.set_xlim(0, 200); ax.set_ylim(-70, 75); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]"); ax.legend(fontsize=7, loc="lower right")
ax.set_title("Above and below the interface the fluids move in opposite directions", fontsize=10)
savefig(fig, "ch07", "c13_interface"); plt.show()
""", see=r"""Above and below the interface (navy line, height ×10) the arrows point in opposite directions; the sheet strength
changes sign between crests and troughs (our remake of Fig. 7.24 `N190`).""",
    read=r"""The arrows shrink like $e^{-k\lvert z\rvert}$: the motion is confined to about a wavelength around the interface.""",
    change=r"""…ρ₁ → 0 (air over water): the upper arrows vanish and (7.95) returns $\omega=\sqrt{gk}$, the deep-water wave of C04.""")
note("N108, N192", "Dead water", r"""
A ship in a fjord with a thin fresh layer makes interfacial waves and feels an unexpected drag — Bjerknes' explanation of
the old sailors' puzzle (the book's Fig. 7.26; two-layer flows return in Ch. 13).
""")
whatif(r"""
…the upper layer had a free surface a finite distance above? A second wave (on the surface) enters and the two interact —
C14's barotropic and baroclinic modes.
""")

# ---- C14 ------------------------------------------------------------------------------------------------------------
core("C14", r"Barotropic and baroclinic modes; reduced gravity", r"""
What are the two ways a layer over deep water can oscillate — and why is one of them almost invisible at the surface?
""", eqs=("7.116", "7.117"))
note("N109", "The set-up", r"""
(the book's Fig. 7.27) An upper layer of thickness H with a free surface, over an infinitely deep lower layer; the origin
at the mean free surface, the interface at z = −H. Two modes: surface and interface in phase, or in antiphase.
""")
problem(r"""
The summer thermocline of a lake sits 10–50 m down; tides and winds move both the surface and the thermocline. Ocean
models are built from exactly two kinds of motion: the barotropic one, which moves the whole column and shows at the
surface, and the baroclinic one, which heaves the thermocline by metres while the surface barely twitches. We want both
from one calculation.
""")
nb.md(r"""
#### The idea

```
barotropic:  surface ~~~   interface ~ (in phase, smaller by e^{−kH})     ω² = gk          fast
baroclinic:  surface ~ (tiny, antiphase)   interface ~~~~~ (large)        ω² ≈ g′kH·k      slow
```
""")
note("N110–N118", "The problem", r"""
$\phi_2\to0\ \text{at}\ z\to-\infty$ *(7.97)*; $\frac{\partial\phi_1}{\partial z}=\frac{\partial\eta}{\partial t}\ \text{at}\ z=0$
*(7.98)* (the book prints ∂φ₁/dz; read ∂z); $\frac{\partial\phi_1}{\partial t}+g\eta=0\ \text{at}\ z=0$ *(7.99)*;
$\frac{\partial\phi_1}{\partial z}=\frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}{\partial t}\ \text{at}\ z=-H$ *(7.100)*;
$\rho_1\frac{\partial\phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+\rho_2g\zeta\ \text{at}\ z=-H$
*(7.101)*; with $\eta=ae^{i(kx-\omega t)}$ *(7.102)* (a real), $\zeta=be^{i(kx-\omega t)}$ *(7.103)* (b complex: a phase
difference is allowed), $\phi_1=(Ae^{kz}+Be^{-kz})e^{i(kx-\omega t)}$ *(7.104)* and $\phi_2=Ce^{kz}e^{i(kx-\omega t)}$
*(7.105)*. The book prints (7.105) with $e^{i(kz-\omega t)}$; that cannot satisfy Laplace's equation or (7.100) for all x —
`ch07.two_layer_sympy()["printed_7105_residual"]` ≠ 0 shows it (printed in D29's check).
""")
remind([
    ("complex 2 × 2 systems and factorising in sympy", "`sp.solve` handles complex unknowns, `sp.I` is i, `sp.factor` exposes common factors (Ch. 6 P158)."),
])
D("D29", ref="7.109", check_src=r"""
import sympy as sp                                         # symbolic algebra with complex numbers
x, z, t = sp.symbols('x z t', real=True)                   # coordinates and time
k, w, g, H, a = sp.symbols('k omega g H a', positive=True) # wavenumber, frequency, gravity, layer depth, surface amplitude
A, B, C, b = sp.symbols('A B C b')                         # unknown (complex) constants
E_ = sp.exp(sp.I*(k*x - w*t))                              # the common wave factor e^{iθ}
phi1 = (A*sp.exp(k*z) + B*sp.exp(-k*z))*E_                 # (7.104) upper-layer potential
phi2 = C*sp.exp(k*z)*E_                                    # (7.105) with e^{i(kx−ωt)} (the printed kz is a slip)
eta, zeta = a*E_, b*E_                                     # (7.102), (7.103)
eqs = [sp.diff(phi1, z).subs(z, 0) - sp.diff(eta, t),      # (7.98) surface kinematics — D29 step 3
       sp.diff(phi1, t).subs(z, 0) + g*eta,                # (7.99) surface dynamics — step 4
       (sp.diff(phi1, z) - sp.diff(phi2, z)).subs(z, -H),  # (7.100) first equality — step 7
       sp.diff(phi1, z).subs(z, -H) - sp.diff(zeta, t)]    # (7.100) second equality — step 9
sol = sp.solve([sp.simplify(e/E_) for e in eqs], [A, B, C, b], dict=True)[0]   # steps 5, 6, 8, 10
book = {A: -sp.I*a/2*(w/k + g/w), B: sp.I*a/2*(w/k - g/w),
        C: -sp.I*a/2*(w/k + g/w) - sp.I*a/2*(w/k - g/w)*sp.exp(2*k*H),
        b: a/2*(1 + g*k/w**2)*sp.exp(-k*H) + a/2*(1 - g*k/w**2)*sp.exp(k*H)}   # (7.106)–(7.109)
print("solved − book:", [sp.simplify(sol[s_] - book[s_]) for s_ in (A, B, C, b)])   # [0, 0, 0, 0]
phi2_printed = C*sp.exp(k*z)*sp.exp(sp.I*(k*z - w*t))      # the printed (7.105)
r = (sp.diff(phi1, z) - sp.diff(phi2_printed, z)).subs(z, -H).subs(sol)   # (7.100) with it, after steps 5–10
print("printed (7.105) satisfies (7.100) for every x:", sp.simplify(sp.diff(r, x)) == 0)   # False
print("library:", ch07.two_layer_sympy()["residuals"], "| printed Laplace residual:", ch07.two_layer_sympy()["printed_7105_residual"])
""")
note("N119–N122", "The constants", r"""
$A=-\frac{ia}2\big(\frac\omega k+\frac g\omega\big)$ *(7.106)*, $B=\frac{ia}2\big(\frac\omega k-\frac g\omega\big)$
*(7.107)*, $C=-\frac{ia}2\big(\frac\omega k+\frac g\omega\big)-\frac{ia}2\big(\frac\omega k-\frac g\omega\big)e^{2kH}$
*(7.108)* and (7.109) below.
""", equation=EQ["7.109"], ref="7.109")
D("D30", ref="7.110", check_src=r"""
import sympy as sp                                              # symbolic algebra
k, w, g, H, a, r1, r2 = sp.symbols('k omega g H a rho1 rho2', positive=True)   # wave, layer and densities
A = -sp.I*a/2*(w/k + g/w); B = sp.I*a/2*(w/k - g/w)             # (7.106), (7.107)
C = A - B*sp.exp(2*k*H)                                         # (7.108), D29 step 8
b = sp.I*k/w*(A*sp.exp(-k*H) - B*sp.exp(k*H))                   # D29 step 10
lhs = -sp.I*w*(r1*(A*sp.exp(-k*H) + B*sp.exp(k*H)) - r2*C*sp.exp(-k*H))   # D30 step 3, left side
res = lhs - (r2 - r1)*g*b                                       # residual of (7.101) after cancelling e^{iθ}
s = w**2/(g*k)                                                  # step 7: s = ω²/gk
book = (s - 1)*(s*(r1*sp.sinh(k*H) + r2*sp.cosh(k*H)) - (r2 - r1)*sp.sinh(k*H))   # (7.110)
print("residual − (ag²k/ω²)·(7.110):", sp.simplify((res - a*g**2*k/w**2*book).rewrite(sp.exp)))   # 0 (step 12)
print("factored residual:", sp.factor(sp.simplify((res*w**2/(a*g**2*k)).rewrite(sp.exp))))   # shows (gk − ω²) (step 10)
print("library ratio (7.101)/(7.110):", ch07.two_layer_sympy()["common_factor"])   # a g² k/ω²
""")
note("N123", "Two factors, two modes", r"""
The product (7.110) vanishes when either factor does.
""", equation=EQ["7.110"], ref="7.110")
D("D31", ref="7.117")
note("N124–N130", "The two modes", r"""
First root $\omega^2=gk$ *(7.111)* — a deep-water surface wave — with $b=ae^{-kH}$ *(7.112)*: the interface moves in phase,
smaller by $e^{-kH}$: the **barotropic** mode (surfaces of constant pressure and density coincide). Second root
$\omega^2=\frac{gk(\rho_2-\rho_1)\sinh kH}{\rho_2\cosh kH+\rho_1\sinh kH}$ *(7.113)*, which becomes (7.95) as kH → ∞, with
$\eta=-\zeta\big(\frac{\rho_2-\rho_1}{\rho_1}\big)e^{-kH}$ *(7.114)*: surface and interface in antiphase and the surface
tiny — the **baroclinic** (internal) mode; u reverses across the interface. Long waves (kH ≪ 1):
$\omega^2=kg\big(\frac{\rho_2-\rho_1}{\rho_2}\big)kH$ *(7.115)*, so $c=[g'H]^{1/2}$ *(7.116)* with the **reduced gravity**
$g'=g\big(\frac{\rho_2-\rho_1}{\rho_2}\big)$ *(7.117)* — the shallow-water speed with g replaced by g′;
$\eta=-\zeta\big(\frac{\rho_2-\rho_1}{\rho_1}\big)$ *(7.118)*; and $p'=-\rho_1\frac{\partial\phi_1}{\partial t}=\rho_1g\eta$
*(7.119)*: the upper layer is hydrostatic for long waves — shallow means hydrostatic again (C04).
""")
remind([
    ("reduced gravity and its two conventions", "a parcel lighter by Δρ feels g′ = gΔρ/ρ_ref (Ch. 4 P131); which ρ_ref is used differs between books — the callout below gives the numbers."),
])
confusion(r"""
**Which ρ in g′?** The book's $g'=g\big(\frac{\rho_2-\rho_1}{\rho_2}\big)$ *(7.117)* divides by the **lower** density ρ₂;
Ch. 4's `reduced_gravity` (P131) divided by ρ₁. Ocean (1025 over 1027 kg/m³): 0.01910 vs 0.01914 m/s² (0.2 % apart); oil
(800) over water (1000): 1.96 vs 2.45 m/s² (25 %). Our code: `W.reduced_gravity_book(rho1, rho2, ref="lower")` for the
book, `ref="upper"` for Ch. 4's. Ch. 13 uses a reference ρ₀ — for the ocean all three agree to 0.2 %.
""")
nb.code(r"""
for r1, r2 in ((1025.0, 1027.0), (800.0, 1000.0)):         # ocean thermocline; oil over water [kg/m³]
    gl = W.reduced_gravity_book(r1, r2, ref="lower")       # (7.117): divide by ρ₂
    gu = W.reduced_gravity_book(r1, r2, ref="upper")       # Ch. 4: divide by ρ₁
    print(f"ρ₁ = {r1:.0f}, ρ₂ = {r2:.0f}: g′ = {gl:.5f} (book) vs {gu:.5f} m/s² (Ch. 4), "
          f"Ch. 4's own function {similarity.reduced_gravity(r1, r2, g=9.81):.5f}; ratio {gl/gu:.4f}")
""")
nb.worked_example("a thermocline 50 m down", r"""
H = 50 m, ρ₁ = 1000, ρ₂ = 1002 kg/m³.

1. $g'=9.81\times2/1002=0.0196$ m/s².
2. $c=\sqrt{g'H}=\sqrt{0.98}=0.99$ m/s — against $\sqrt{gH}=22.1$ m/s for a long surface wave on the same 50 m.
3. A 10 m heave of the thermocline shows at the surface as $\eta=-10\times2/1000=-0.02$ m: 2 cm, in antiphase.
4. For λ = 1 km (kH = 0.31, not quite long) (7.113) gives T = 19.5 min, c = 0.85 m/s.
""")
nb.code(r"""
for lam in (100.0, 1000.0, 5000.0):                        # three wavelengths [m]
    k = 2*np.pi/lam                                        # [rad/m]
    bt, bc = ch07.two_layer_free_surface_omega(k, 50.0, 1000.0, 1002.0)   # (7.111), (7.113)
    st = ch07.two_layer_state(k, 50.0, 1000.0, 1002.0)     # everything for the baroclinic mode
    print(f"λ = {lam:5.0f} m: barotropic T = {2*np.pi/bt:5.1f} s, baroclinic T = {2*np.pi/bc/60:5.1f} min, "
          f"c_bc = {st['c_bc']:.3f} m/s, η/ζ = {st['eta_over_zeta']:.5f}")
print("long-wave c = √(g′H) =", ch07.two_layer_long_wave_speed(50.0, 1000.0, 1002.0), "m/s")   # (7.116)
m = ch07.two_layer_modes(2*np.pi/1000, 50.0, 1000.0, 1002.0, mode="baroclinic")
print("baroclinic η/ζ at λ = 1 km:", m["eta_over_zeta"], "| book form (7.114):", m["eta_over_zeta_book"])
print("rigid lids, h₁ = 50 m over h₂ = 200 m (rigid-lid form):", ch07.two_layer_rigid_lid_omega(2*np.pi/1000, 50.0, 200.0, 1000.0, 1002.0), "rad/s")
""", explain=r"""
1. Baroclinic periods 4.2, 19.5 and 86.9 min against barotropic 8.0, 25.3 and 56.6 s: the internal mode is 30–90 times
   slower.
2. For λ = 1 km a 1 m heave of the interface lifts the surface by only −1.5 mm (antiphase), approaching −2 mm per metre
   (−Δρ/ρ₁) for long waves.
3. `two_layer_rigid_lid_omega` is the two-layer form between rigid lids (note N131): here 0.0053 rad/s.
""")
nb.check_agree(r"""
k, H, r1, r2 = 2*np.pi/1000, 50.0, 1000.0, 1002.0          # λ = 1 km on the thermocline of the tiny example
S, Cc = np.sinh(k*H), np.cosh(k*H)                         # sinh kH, cosh kH
# (7.110) expanded as a quadratic in s = ω²/gk:
coeffs = [r1*S + r2*Cc, -(r1*S + r2*Cc + (r2 - r1)*S), (r2 - r1)*S]
s_roots = np.roots(coeffs)                                 # the two roots s
w_mine = np.sqrt(G*k*s_roots)                              # ω = √(gk s) [rad/s]
print("s roots:", s_roots, "→ ω:", w_mine)
assert np.allclose(sorted(w_mine), sorted(ch07.two_layer_free_surface_omega(k, H, r1, r2)), rtol=1e-10)   # both modes
""")
nb.md(r"**What does this show?** The quadratic in s = ω²/(gk) has two roots. s = 1 gives ω = 0.248 rad/s, the barotropic mode ($\omega^2=gk$, as for one deep fluid); s ≈ 5 × 10⁻⁴ gives ω ≈ 0.0054 rad/s, the baroclinic mode, about 46 times slower. Both are the zeros of the two factors of " + E("7.110") + r", as `two_layer_free_surface_omega` finds them.")
nb.figure(r"""
k, H, r1, r2 = 2*np.pi/500, 50.0, 1000.0, 1002.0           # λ = 500 m
xs = np.linspace(0, 1000, 300)                             # two wavelengths [m]
fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)   # barotropic | baroclinic
for axp, mode in zip(axs, ("barotropic", "baroclinic")):   # one panel per mode
    md = ch07.two_layer_modes(k, H, r1, r2, mode=mode)     # constants and the amplitude ratio
    ph = k*xs                                               # phase at t = 0
    zeta_true = np.real(md["b"]*np.exp(1j*ph))              # interface displacement for a surface amplitude a = 1 m [m]
    if mode == "barotropic":   # the in-phase mode
        eta_d, zeta_d = 10*np.cos(ph), 10*zeta_true         # both drawn ×10 (the interface is the smaller one)
    else:   # the antiphase mode
        sc = 10/abs(md["b"])                                # draw the interface with a 10 m amplitude …
        zeta_d, eta_d = sc*zeta_true, 50*sc*np.cos(ph)      # … and the surface 50 times larger again
    two_layer_sketch(axp, H, 1000, eta=eta_d, zeta=zeta_d, x=xs)   # light upper layer over navy lower layer
    zz = np.linspace(-2.5*H, 0, 6)                          # arrow heights [m]
    for z_ in zz:                                           # u at x = 0 (crest) and x = 250 (trough), sign only
        up = z_ > -H   # True in the upper layer
        amp = np.real(1j*k*((md["A"]*np.exp(k*z_) + md["B"]*np.exp(-k*z_)) if up else md["C"]*np.exp(k*z_)))   # u = ∂φ/∂x at a crest of η
        for x_, sg in ((500, 1), (250, -1)):               # under a crest of η (u ∝ +amp) and a trough (−amp)
            axp.annotate("", xy=(x_ + 40*np.sign(amp)*sg, z_), xytext=(x_, z_),   # an arrow showing the direction of u
                         arrowprops=dict(arrowstyle="->", color=COLORS["ink"], lw=1.2))   # arrow style
    axp.set_title(f"{mode}: η/ζ = {md['eta_over_zeta']:+.4f}" + (" (surface ×50 extra)" if mode == "baroclinic" else ""), fontsize=9)   # the true amplitude ratio
    axp.set_xlabel("x [m]"); axp.set_ylim(-2.6*H, 25)   # axes
axs[0].set_ylabel("z [m]")   # unit
fig.suptitle("Two ways to oscillate: in phase (barotropic) and in antiphase (baroclinic)", fontweight="bold")   # the message
savefig(fig, "ch07", "c14_modes"); plt.show()   # save, then draw
""", see=r"""Left, the barotropic mode: surface and interface in phase, the arrows point the same way in both layers; right,
the baroclinic mode: a tiny surface signal in antiphase (drawn 50 times larger again) and a large interface wave, with u
reversing across the interface (our remake of Fig. 7.27 `N193`).""",
    read=r"""The printed ratios: η/ζ = +e^{kH} ≈ 1.9 for the barotropic mode (the interface is the smaller one) and −0.0011 for
the baroclinic one.""",
    change=r"""…Δρ → 0: the baroclinic surface signal vanishes and its speed goes to zero like √Δρ.""")
note("N131", "Two thin layers", r"""
When both layers are thin compared with λ and Δρ/ρ ≪ 1 (Boussinesq), the barotropic mode has the same u at every depth and
the baroclinic mode has uniform, opposite u in the two layers (two layers between rigid lids). The profiles
below are $u=\partial\phi/\partial x$ of the two long-wave modes (λ = 5 km, H = 50 m), each scaled to its largest value —
our version of the book's Fig. 7.28 `N194`.
""")
nb.figure(r"""
k, H, r1, r2 = 2*np.pi/5000, 50.0, 1000.0, 1002.0          # a long wave, kH = 0.063
z1, z2 = np.linspace(-H, 0, 50), np.linspace(-4*H, -H, 80) # upper and lower layer depths [m]
fig, ax = plt.subplots(figsize=(5, 3.6))   # one panel
for mode, cl, ls, mk in (("barotropic", COLORS["blue"], "-", None), ("baroclinic", COLORS["accent"], "--", "o")):   # both modes
    md = ch07.two_layer_modes(k, H, r1, r2, mode=mode)   # its constants A, B, C (a = 1 m)
    u1 = np.real(1j*k*(md["A"]*np.exp(k*z1) + md["B"]*np.exp(-k*z1)))   # u = ∂φ₁/∂x at the crest of η
    u2 = np.real(1j*k*md["C"]*np.exp(k*z2))                              # u = ∂φ₂/∂x
    sc = max(np.abs(u1).max(), np.abs(u2).max())   # scale to the largest |u|
    ax.plot(u1/sc, z1, color=cl, lw=2 if mk is None else 1.5, ls=ls, marker=mk, markevery=6, ms=4, label=mode)   # upper layer
    ax.plot(u2/sc, z2, color=cl, lw=2 if mk is None else 1.5, ls=ls, marker=mk, markevery=10, ms=4)   # lower layer
ax.text(0.55, -H/2, "both ≈ 1 in\nthe upper layer", fontsize=7, ha="center", color=COLORS["ink"])   # the curves overlap here
ax.axhline(-H, color=NAVY, lw=1, ls="--"); ax.text(-0.95, -H + 5, "interface", fontsize=7, color=NAVY)   # the interface
ax.axvline(0, color=COLORS["muted"], lw=0.5)   # u = 0
ax.set_xlabel("u / max|u| [–]"); ax.set_ylabel("z [m]"); ax.legend(fontsize=7); ax.set_xlim(-1.1, 1.1)   # axes
ax.set_title("Long waves: uniform u, or a return flow below the interface", fontsize=10)   # the message
savefig(fig, "ch07", "c14_long_wave_u"); plt.show()   # save, then draw
""", see=r"""The barotropic u (blue) is nearly the same at every depth (it decays only like $e^{kz}$ over 1/k ≈ 800 m); the
baroclinic u (purple, dashed with dots) is uniform in the upper layer — on top of the blue line there — and reverses across
the interface into a weak return flow below.""",
    read=r"""In the baroclinic mode the upper layer's flow is balanced by the opposite flow below: here the lower layer is
deep, so the return flow is weak (≈ kH = 6 % of the upper one) but spread over ~1/k; between rigid lids with a thin lower
layer (note N131) it would be as strong as the upper one. Almost no net transport — the signature of internal motion.""",
    change=r"""…the upper layer were much thicker than the wavelength: the profiles would curve like $e^{\pm kz}$ (deep water,
C13).""")
nb.figure(r"""
k = np.logspace(-4, 0, 300)                                # wavenumbers [rad/m]
H, r1, r2 = 50.0, 1000.0, 1002.0                           # thermocline of the tiny example
bt, bc = ch07.two_layer_free_surface_omega(k, H, r1, r2)   # (7.111), (7.113)
fig, ax = plt.subplots(figsize=(6, 3.6))   # one panel
ax.loglog(k, bt, color=COLORS["blue"], lw=2, label="barotropic ω² = gk (7.111)")   # fast branch
ax.loglog(k, bc, color=COLORS["accent"], lw=2, label="baroclinic (7.113)")   # slow branch
ax.loglog(k, ch07.interface_omega(k, r1, r2), "--", color=COLORS["muted"], label="ε√(gk) (7.95)")   # short-wave limit
ax.loglog(k, np.sqrt(W.reduced_gravity_book(r1, r2)*H)*k, ":", color=COLORS["muted"], label="√(g′H)·k (7.116)")   # long-wave limit
ax.set_xlabel("k [rad/m]"); ax.set_ylabel("ω [rad/s]"); ax.legend(fontsize=7)   # axes
ax.set_title("Two branches, far apart (H = 50 m, Δρ = 2 kg/m³)", fontsize=10)   # the message
savefig(fig, "ch07", "c14_two_branches"); plt.show()   # save, then draw
""", see=r"""Two branches of ω(k) about a factor 30 apart; the baroclinic one follows √(g′H)·k for long waves and joins ε√(gk)
for short ones.""",
    read=r"""For kH ≫ 1 the surface is too far away to matter: the interface forgets it, (7.113) → (7.95).""",
    change=r"""…H were doubled: the baroclinic long-wave line moves up by √2.""")
nb.plotly(r"""
k = np.logspace(-4, 0, 120)                                # wavenumbers [rad/m]
def branches(r):                                           # r = Δρ/ρ₂ [–]
    r2 = 1000.0; r1 = r2*(1 - r)                           # densities [kg/m³]
    bt, bc = ch07.two_layer_free_surface_omega(k, 50.0, r1, r2)   # (7.111), (7.113)
    return {"barotropic (7.111)": (k, bt), "baroclinic (7.113)": (k, bc), "ε√(gk) (7.95)": (k, ch07.interface_omega(k, r1, r2))}   # three curves per position
fig = slider_figure(branches, "Δρ/ρ₂", np.logspace(-4, np.log10(0.5), 16 if not FAST else 9), xlabel="k [rad/m]",   # slider: density contrast (log)
                    ylabel="ω [rad/s]", title="The baroclinic branch rises with the density contrast (H = 50 m)")   # labels and message
logaxes(fig, x=(1e-4, 1), y=(1e-6, 5))   # log axes, fixed
recolor(fig, {"barotropic (7.111)": COLORS["blue"], "baroclinic (7.113)": COLORS["accent"], "ε√(gk) (7.95)": COLORS["muted"]},   # colours by meaning
        dashes={"ε√(gk) (7.95)": "dash"})   # limit dashed
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** Drag Δρ/ρ₂: the barotropic branch does not move at all; the baroclinic branch rises like $\sqrt{\Delta\rho}$ and always
joins the two-deep-fluids curve (7.95) for short waves.

**What would change if…** the upper layer were thicker (H = 200 m): the baroclinic long-wave line $\sqrt{g'H}\,k$ would move up by a
factor 2 and bend over into the short-wave curve at a smaller k.
""")
explainer("two_layer_modes", "What is a baroclinic mode?",
          r"""The density-contrast and thickness sliders move both dispersion curves and both mode shapes at once; the reader
sees the interface amplitude dwarf the surface one as Δρ → 0 and recovers the two-deep-fluids case as kH → ∞.""", "",
          ["Preset *ocean thermocline*: compare the two periods in the Explain tab.",
           "Switch between the modes: watch u reverse across the interface in the baroclinic one.",
           "Drag Δρ/ρ up to 0.2 (oil over water) and read how different the two g′ conventions are.",
           "Step the Derivation (D30) to 'factor out (s − 1)'."])
whatif(r"""
…the Earth rotated? The baroclinic long-wave speed $\sqrt{g'H}\approx1$ m/s divided by the Coriolis parameter
f ≈ 10⁻⁴ s⁻¹ gives the internal (baroclinic) Rossby radius ≈ 10 km — the size of ocean eddies (Ch. 13).
""")

# =====================================================================================================================
# A.8 §7.8 — R12–R20, C15 (N132–N142, N145–N147, N196, D32–D34, P177), C16 (N143, N144, N148–N166, N195, N197–N199,
#            D35–D37, P178, A6, E9), S01
# =====================================================================================================================
nb.section("7.8", "Internal Waves in a Continuously Stratified Fluid", intro=r"""
**What is this section about?** Now the density changes smoothly with height, as in most of the ocean and atmosphere.
Displaced parcels bob at the buoyancy frequency N, and waves can travel in any direction — with astonishing rules: the
frequency depends only on the direction of the wave, never on its wavelength; the water moves along the crests; and the
energy travels at right angles to the crests.
""")
nb.recap("R12", "Assumptions", r"""
Boussinesq (density differences matter only in the weight), inviscid, small amplitude (advection dropped), and frequencies
well above the Coriolis frequency (rotation joins in Ch. 13). The momentum equations are Ch. 4's Boussinesq form
$\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ without advection and
viscosity.
""", where="Ch. 4 §4.9")
nb.recap("R13", "Density and continuity equations", r"""
$\frac{D\rho}{Dt}\equiv\frac{\partial\rho}{\partial t}+\mathbf u\cdot\nabla\rho=0$ *(4.9)* (each parcel keeps its density)
and $\nabla\cdot\mathbf u=0$ *(4.10)* (incompressible flow). Below: the internal wave of C16 is divergence-free to
round-off.
""", where="Ch. 4 §4.2")
nb.code(r"""
r = ch07.internal_wave_fields(0.3, -0.2, 1.0, 1.0, 1.0, 0.01, 1e-3, residuals=True)["residuals"]   # k = m = 1 rad/m, N = 0.01 rad/s
print({n: f"{v:.1e}" for n, v in r.items()})               # (4.10), (7.128), (7.130), (7.131) by central differences
""", explain=r"""
The plane internal wave of C15 is put into continuity $\nabla\cdot\mathbf u=0$ *(4.10)* and the linear equations
$\frac{\partial u}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial x}$ *(7.128)*,
$\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial z}-\frac{\rho'g}{\rho_0}$ *(7.130)* and
$\frac{\partial\rho'}{\partial t}-\frac{N^2\rho_0}gw=0$ *(7.131)*; every residual is at the size of the difference-stencil
error — the fields satisfy the equations D32 derives below.
""")
nb.recap("R14", "Linear Boussinesq momentum", r"""
$\frac{\partial u}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial x}$ *(7.120)*,
$\frac{\partial v}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial y}$ *(7.121)*,
$\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial z}-\frac{\rho g}{\rho_0}$ *(7.122)*: ρ₀ in the
inertia, the true ρ only in the weight.
""", where="Ch. 4 §4.9")
nb.recap("R15", "Why Dρ/Dt = 0", r"""
Heat and salt do not diffuse on wave time scales, and the density depends on temperature and salinity but not on pressure
(δρ/ρ = −αδT + βδS), so a parcel keeps its density: $\frac{D\rho}{Dt}=0$ *(4.9)*.
""", where="Ch. 1 §1.10, Ch. 4 §4.2")
nb.recap("R16", "The resting state", r"""
Before the waves, the fluid is hydrostatic: $0=-\frac1{\rho_0}\frac{d\bar p}{dz}-\frac{\bar\rho g}{\rho_0}$ *(7.123)*.
""", where="Ch. 1 §1.7")
nb.recap("R17", "Background plus perturbation", r"""
$p=\bar p(z)+p',\ \rho=\bar\rho(z)+\rho'$ *(7.124)* — ⚠️ here p′ is measured from the stratified background $\bar p(z)$,
not from −ρgz as in $p'\equiv p+\rho gz$ *(7.30)*.
""", where="Ch. 4 §4.9")
nb.recap("R18", "The buoyancy frequency", r"""
$N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ *(7.127)* — $N^2=-\frac{g}{\rho}\big(\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\big)$
*(1.29)* with no adiabatic gradient: a displaced parcel oscillates at N.
""", where="Ch. 1 §1.10")
nb.code(r"""
N2 = stratification.brunt_vaisala_sq(1025.0, -1.045e-2, 0.0, g=9.81)   # (1.29): ρ = 1025 kg/m³, dρ/dz = −0.01045 kg/m⁴
print(f"N² = {N2:.3e} 1/s², N = {np.sqrt(N2):.4f} rad/s, buoyancy period {2*np.pi/np.sqrt(N2)/60:.1f} min")   # a thermocline
""")
nb.recap("R19", "Potential density", r"""
Compressibility is handled by using the potential density (the density brought adiabatically to a reference pressure,
oceanographers' "sigma-theta") in N, so the incompressible equations still apply.
""", where="Ch. 1 §1.10")
nb.recap("R20", "Internal waves are rotational", r"""
Where density surfaces tilt against pressure surfaces the baroclinic torque ∇ρ × ∇p/ρ² makes vorticity (Ch. 5 §5.4), so
internal waves have no velocity potential — C13's vortex sheet is the thin-layer limit.
""", where="Ch. 5 §5.4")

# ---- C15 ------------------------------------------------------------------------------------------------------------
core("C15", r"Internal waves: frequency set by direction", r"""
If the buoyancy frequency N is the only frequency the fluid knows, what sets the frequency of a wave?
""", eqs=("7.139",))
problem(r"""
Fill a tank with salt water that gets fresher toward the top and wiggle a rod: the disturbance does not spread in rings
like on a pond but in slanted beams. In the ocean, tides flowing over ridges radiate such beams through the thermocline; in
the atmosphere, wind over mountains makes lee waves. We want the rule that fixes their frequency.
""")
nb.md(r"""
#### The idea

```
parcel pushed straight up:           it bobs at N (Ch. 1)
parcel pushed along a slope at θ:    only the vertical part of gravity's pull restores it
                                     ⇒ it bobs at N cos θ, whatever the wavelength
so:  ω = N cos θ,   0 < ω < N,   θ = angle of K above the horizontal = angle of the motion from the vertical
```
""")
note("N132", "Put the split into the density equation", r"""
**Linearisation about a base state** keeps products of one small quantity with a background quantity ($w\,d\bar\rho/dz$) and
drops products of two small quantities ($u\,\partial\rho'/\partial x$) — Ch. 2 P68's orders of smallness.
""", equation=EQ["7.125"], ref="7.125")
remind([
    ("material derivative", r"$\frac{D}{Dt}=\frac{\partial}{\partial t}+\mathbf u\cdot\nabla$, the rate of change following a parcel (Ch. 3 §3.4)."),
])
D("D32", ref="7.131")
note("N133–N135", "The linear equations", r"""
$\frac{\partial\rho'}{\partial t}+w\frac{d\bar\rho}{dz}=0$ *(7.126)*: density changes at a point only because the vertical
motion carries the background up or down; the perturbation momentum equations
$\frac{\partial u}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial x}$ *(7.128)*,
$\frac{\partial v}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial y}$ *(7.129)*,
$\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial z}-\frac{\rho'g}{\rho_0}$ *(7.130)* look like
(7.120)–(7.122) with ρ′, p′; and $\frac{\partial\rho'}{\partial t}-\frac{N^2\rho_0}gw=0$ *(7.131)*.
""")
P("P177", "operator elimination for linear PDEs", r"""
With several linear PDEs in several unknowns, apply derivatives (∂/∂t, ∂/∂z, the horizontal Laplacian $\nabla_H^2$) to whole
equations and substitute one into another until a single unknown is left — like eliminating variables from linear
equations. Allowed because the derivatives of smooth fields commute (Schwarz, Ch. 4 P121) and the coefficients do not
depend on the variable being differentiated.
""", code=r"""
import sympy as sp                                # symbolic check of one elimination move
x, t = sp.symbols('x t'); u, p = sp.Function('u')(x, t), sp.Function('p')(x, t)   # two unknown fields
eq1 = sp.Eq(sp.diff(u, t), -sp.diff(p, x))        # ∂u/∂t = −∂p/∂x
print(sp.simplify(sp.diff(eq1.lhs, x) - sp.diff(sp.diff(u, x), t)))   # 0: ∂/∂x of ∂u/∂t is ∂/∂t of ∂u/∂x
""")
D("D33", ref="7.134", check_src=r"""
import sympy as sp                                        # symbolic algebra
x, y, z, t = sp.symbols('x y z t', real=True)             # coordinates and time
rho0, g = sp.symbols('rho0 g', positive=True)             # reference density, gravity
N = sp.Function('N')(z)                                   # buoyancy frequency, allowed to depend on z
p = sp.Function('p')(x, y, z, t); w = sp.Function('w')(x, y, z, t)   # perturbation pressure and vertical velocity
lap_H = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2)     # horizontal Laplacian ∇_H²
eq132 = lap_H(p)/rho0 - sp.diff(w, z, t)                  # (7.132) as "left − right": ∇_H²p′/ρ₀ − w_zt (steps 1–4)
eq133 = sp.diff(p, t, z)/rho0 + sp.diff(w, t, 2) + N**2*w # (7.133) as "left − right": p′_tz/ρ₀ + w_tt + N²w (steps 5–7)
combo = sp.diff(eq132, t, z) - lap_H(eq133)               # steps 8–10: ∂_t∂_z of (7.132) minus ∇_H² of (7.133)
w134 = sp.diff(w, t, 2, z, 2) + lap_H(sp.diff(w, t, 2)) + N**2*lap_H(w)   # (7.134): ∂_tt∇²w + N²∇_H²w (steps 11–12)
print("p′ eliminated, even for N(z):", sp.simplify(combo + w134))          # 0: the elimination is exact
print("library, every step's residual:", {k_: v for k_, v in ch07.boussinesq_linear_sympy().items() if k_.startswith("r")})
""")
note("N136–N138", "Where D33 passes", r"""
D33 passes through $\frac1{\rho_0}\nabla_H^2p'=\frac{\partial^2w}{\partial z\,\partial t}$ *(7.132)*
($\nabla_H^2\equiv\partial^2/\partial x^2+\partial^2/\partial y^2$), $\frac1{\rho_0}\frac{\partial^2p'}{\partial t\,\partial
z}=-\frac{\partial^2w}{\partial t^2}-N^2w$ *(7.133)*, and ends at the w-equation (7.134) — which holds for N(z) too.
""", equation=EQ["7.134"], ref="7.134")
note("N139, N140", "Anisotropy and the plane-wave trial", r"""
The medium has no preferred horizontal direction but a preferred vertical, so ω depends on all of K:
$\omega=\omega(k,l,m)=\omega(\mathbf K)$ *(7.135)* — anisotropic. Try a plane wave $w=w_0e^{i(kx+ly+mz-\omega t)}$
*(7.136)*. **Angle of a vector from its components**: cos θ = |k|/K, computed with `np.arccos`, degrees with
`np.degrees`.
""")
remind([
    ("plane-wave trial (7.136), ∂ → ik and −iω", r"for a field $\propto e^{i(\mathbf K\cdot\mathbf x-\omega t)}$ every ∂/∂x brings ik, ∂/∂z brings im and ∂/∂t brings −iω (primer P176, C13)."),
])
D("D34", ref="7.139")
note("N141, N142", "The internal-wave dispersion relation", r"""
$\omega^2=\frac{k^2+l^2}{k^2+l^2+m^2}N^2$ *(7.137)*; in the x–z plane $\omega=\frac{kN}{\sqrt{k^2+m^2}}=\frac{kN}K$
*(7.138)* — written for k > 0; for a wave travelling toward −x use |k| (our code does).
""")
note("N145", "The two limits", r"""
θ = 0 (K horizontal, m = 0) gives ω = N — vertical columns bobbing, Ch. 1's parcel; θ → π/2 gives ω → 0, where the wave
solution says nothing and the full equations must be used (next note). The code checks the ω = N limit against Ch. 1's
parcel oscillation.
""")
nb.code(r"""
t = np.linspace(0, 1000, 5)                                # times [s]
f = ch07.internal_wave_fields(0.0, 0.0, t, 0.01, 0.0, 0.01, 1e-3)   # m = 0: K horizontal; N = 0.01 rad/s; ŵ = 1 mm/s
z = stratification.parcel_displacement(t, 0.0, 1e-4, w0=1e-3)       # Ch. 1: ζ″ + N²ζ = 0 from ζ = 0 with speed ŵ
print(f["zeta_particle"], "|", z)                          # the particle's vertical displacement both ways [m]
assert np.allclose(f["zeta_particle"], z, atol=1e-12)      # the ω = N wave is Ch. 1's bobbing parcel
""")
note("N146, N147, N196", "Layered flows and blocking", r"""
A possible steady solution of (4.10) and (7.128)–(7.131): w = p′ = ρ′ = 0 with u, v any horizontally non-divergent field,
$\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$ *(7.142)*, at each level separately — strong stratification
allows flat, layered flows (pancake eddies, cloud sheets seen from aircraft, Ch. 13). In strong stratification a 2-D body
blocks a whole horizontal layer of fluid ahead of it (the book's Fig. 7.30: the fluid cannot go over or under — orographic
blocking in Ch. 13).
""")
nb.code(r"""
print(ch07.layered_flow_check(lambda x, y, z: -np.sin(y), lambda x, y, z: np.sin(x), 0.3, 0.7))   # (7.142) and the rest: 0
""")
nb.worked_example("a thermocline internal wave", r"""
N = 0.01 rad/s.

1. Buoyancy period 2π/N = 628 s = 10.5 min — no internal wave can be faster.
2. A wave with K at θ = 45° to the horizontal: ω = N cos 45° = 0.00707 rad/s, T = 14.8 min.
3. The same direction with a 10× longer wavelength: still 14.8 min — only the direction counts.
4. For T = 12 h (a tide): cos θ = (2π/43 200)/0.01 = 0.0145, θ = 89.2° — K almost vertical, the energy beams almost
   horizontal.
""")
nb.code(r"""
print(ch07.internal_wave_omega(1.0, 1.0, 0.01), ch07.internal_wave_omega(0.1, 0.1, 0.01),
      ch07.internal_wave_omega(-1.0, 1.0, 0.01), "rad/s")   # (7.138): size and the sign of k do not matter
print(f"θ for ω/N = 0.707: {np.degrees(ch07.beam_angle(0.00707, 0.01)):.2f}°")   # (7.139) inverted
s = ch07.boussinesq_linear_sympy()                         # D32–D34 in sympy (cached)
print({n: v for n, v in s.items() if n.startswith("r")})   # every residual 0
""", explain=r"""
1. Three wavenumber vectors along the same 45° direction or its mirror — ω = 0.0070711 rad/s every time: the frequency
   knows only the direction of K.
2. `beam_angle` inverts $\omega=N\cos\theta$ *(7.139)*: 45.0°.
3. `boussinesq_linear_sympy` checks every equation of D32–D34 symbolically, including $\omega^2=\frac{k^2+l^2}{k^2+l^2+m^2}N^2$
   *(7.137)*.
""")
nb.check_agree(r"""
Nf, kk, mm = 1.0, 1.0, 1.0                                 # N [rad/s], k, m [rad/m] (unit scales: round-off stays small)
w_right = ch07.internal_wave_omega(kk, mm, Nf)             # (7.138): 0.707 rad/s
def res(om, x=0.3, z=-0.2, t=1.0, h=1e-2):                 # residual of (7.134) for w = cos(kx + mz − ωt), by hand
    wf = lambda X, Z, T: np.cos(kk*X + mm*Z - om*T)        # the trial plane wave
    d2 = lambda f, X, Z, T, dx, dz, dt: (f(X + dx, Z + dz, T + dt) - 2*f(X, Z, T) + f(X - dx, Z - dz, T - dt))/h**2
    lap = lambda X, Z, T: d2(wf, X, Z, T, h, 0, 0) + d2(wf, X, Z, T, 0, h, 0)   # ∇²w
    return d2(lap, x, z, t, 0, 0, h) + Nf**2*d2(wf, x, z, t, h, 0, 0)            # ∂_tt∇²w + N²∂_xx w
res_right, res_wrong = res(w_right), res(Nf)               # with ω from (7.138), and with the wrong ω = N
print(f"residual with the right ω: {res_right:.1e}; with ω = N: {res_wrong:.2f}")
assert abs(res_right) < 1e-4*abs(res_wrong)                # only (7.138) solves (7.134) (what is left is the O(h²) stencil error)
assert np.allclose(res_right, ch07.w_equation_residual(lambda X, Z, T: np.cos(kk*X + mm*Z - w_right*T), 0.3, -0.2, 1.0, Nf), atol=1e-8)
""")
P("P184", "np.errstate(invalid=\"ignore\")", r"""
numpy normally prints a RuntimeWarning when an operation has no valid answer, such as 0/0 → NaN. Inside
`with np.errstate(invalid="ignore"):` that particular complaint is switched off for the lines in the block only. We use
it where a NaN is expected and harmless: at K = 0 the direction of K, and so ω = N cos θ, is undefined.""", code=r"""
import numpy as np                                        # arrays
with np.errstate(invalid="ignore"):                       # 0/0 is expected here: stay quiet
    print(np.array([0.0, 1.0])/np.array([0.0, 2.0]))      # [nan 0.5]
""")
nb.figure(r"""
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6))       # ω/N vs θ | the (k, m) plane
th = np.linspace(0, np.pi/2, 200)                          # angle of K above the horizontal [rad]
a1.plot(np.degrees(th), np.cos(th), color=COLORS["accent"], lw=2)   # (7.139)
for d in (0, 45, 89):                                      # three directions [deg]
    a1.plot(d, np.cos(np.radians(d)), "o", color=COLORS["orange"]); a1.text(d + 2, np.cos(np.radians(d)) + 0.03, f"{d}°", fontsize=8)
a1.set_xlabel("θ, angle of K above the horizontal [deg]"); a1.set_ylabel("ω/N [–]"); a1.set_title("ω = N cos θ (7.139)", fontsize=9)
kg = np.linspace(-2, 2, 201); KX, MZ = np.meshgrid(kg, kg) # the wavenumber plane [rad/m]
with np.errstate(invalid="ignore"):                        # K = 0 (the centre) has no direction: a quiet NaN there
    WN = ch07.internal_wave_omega(KX, MZ, 1.0)             # ω/N on the plane (N = 1)
cs = a2.contour(KX, MZ, WN, levels=[0.2, 0.5, 0.8], colors=[COLORS["accent"]], linewidths=1)   # straight rays
a2.clabel(cs, fontsize=7)
for L_ in (0.6, 1.2, 1.8):                                 # three K of different length on the 45° ray
    a2.annotate("", xy=(L_*np.cos(np.pi/4), L_*np.sin(np.pi/4)), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["ink"]))
a2.set_aspect("equal"); a2.set_xlabel("k [rad/m]"); a2.set_ylabel("m [rad/m]"); a2.set_title("contours of ω/N: straight rays", fontsize=9)
fig.suptitle("Frequency is set by direction, not wavelength", fontweight="bold")
savefig(fig, "ch07", "c15_direction"); plt.show()
""", see=r"""Left, ω/N = cos θ with dots at 0°, 45° and 89°; right, the (k, m) plane with contours of ω/N — straight lines through
the origin — and three K arrows of different length on the 45° line.""",
    read=r"""Moving along a ray (changing |K|) keeps ω; turning K changes it. The three arrows share one frequency.""",
    change=r"""…we drew a surface wave's ω(K) instead: the contours would be circles (ω depends only on |K|) — the opposite
geometry.""")
nb.plotly(r"""
def cross(r):                                              # r = ω/N: four beams and the c, c_g arrows of each
    bx, bz, cx_, cz_, gx, gz = [], [], [], [], [], []
    for ks in (1.0, -1.0):                                 # K pointing right or left
        for ms in (1.0, -1.0):                             # K pointing up or down
            s_ = ch07.internal_wave_state(r, 1.0, 1.0, k_sign=ks, m_sign=ms)   # (7.139), (7.144)–(7.146)
            eb = np.array([s_["cgx"], s_["cgz"]])/max(s_["cg"], 1e-12)          # beam direction = c_g direction
            bx += [0, 2*eb[0], np.nan]; bz += [0, 2*eb[1], np.nan]              # the beam axis (length 2)
            p0 = 1.2*eb                                     # arrows start on the beam
            cx_ += [p0[0], p0[0] + s_["cx"], np.nan]; cz_ += [p0[1], p0[1] + s_["cz"], np.nan]   # c (orange)
            gx += [p0[0], p0[0] + s_["cgx"], np.nan]; gz += [p0[1], p0[1] + s_["cgz"], np.nan]   # c_g (purple)
    return {"beams": (bx, bz), "c (7.144)": (cx_, cz_), "c_g (7.145)": (gx, gz)}
fig = slider_figure(cross, "ω/N", np.linspace(0.05, 0.95, 19 if not FAST else 10), xlabel="x [m]", ylabel="z [m]",
                    title="The angle is set by the frequency: beams at arccos(ω/N) from the vertical",
                    xrange=(-2.5, 2.5), yrange=(-2.5, 2.5), height=520)
recolor(fig, {"beams": COLORS["muted"], "c (7.144)": COLORS["orange"], "c_g (7.145)": COLORS["accent"]}, dashes={"beams": "dash"})
fig.update_yaxes(scaleanchor="x", scaleratio=1)            # equal scales so angles are true
fig.show()   # draw (works on the page: no kernel needed)
""", explain=r"""
**What you see.** Slide ω/N toward 1: the beams stand up toward the vertical and the purple group velocities shrink; toward 0 the beams lie
down. In every beam the orange phase velocity is perpendicular to the purple group velocity.

**What would change if…** N doubled with the same forcing frequency: ω/N would halve, so the beams would lie closer to the horizontal
and energy would spread more sideways than downward.
""")
confusion(r"""
θ in $\omega=N\cos\theta$ *(7.139)* is the angle of **K** above the horizontal. The same θ is the angle of the particle
motion and of the energy beams from the **vertical**. The book's caption of its St Andrew's-cross photograph says "45° with
the horizontal" — true for the beams only because 45° is symmetric. We always print both.
""")
whatif(r"""
…the fluid rotated at f (Ch. 13)? Then $\omega^2=N^2\cos^2\theta+f^2\sin^2\theta$: internal waves live between f and N — the
inertia–gravity waves of the ocean and atmosphere.
""")

# ---- C16 ------------------------------------------------------------------------------------------------------------
core("C16", r"Transverse waves with c perpendicular to c_g: beams and F = c_g E", r"""
Why does an internal wave's energy leave at right angles to the way its crests move?
""", eqs=("7.146", "7.159"))
problem(r"""
Oscillate a cylinder up and down in a stratified tank and four beams appear — the St Andrew's cross; inside each beam the
crests (dark and light stripes) slide *across* the beam, not along it. In the ocean, energy from tides over a ridge goes
down into the abyss in such beams while the phase lines move up. We want to understand this and prove where the energy
goes.
""")
nb.md(r"""
#### The idea

```
incompressible plane wave   ⇒  K·u = 0          the water moves along the crests (transverse)
ω = N|k|/K depends on direction only  ⇒  ∇_K ω ⊥ K   ⇒  c_g ⊥ c
phase up  ⇔  energy down;   horizontal parts of c and c_g point the same way
```
""")
note("N143", "The velocity as a plane wave", "", equation=EQ["7.140"] + r",\ \text{and the same for } v,\ w", ref="7.140")
D("D35", ref="7.141")
note("N144", "Transverse waves", r"""
$\mathbf K\cdot\mathbf u=0$ *(7.141)*: particle motion is perpendicular to K, i.e. along the crests — a transverse (shear)
wave. Surface waves are not transverse: their fields decay in z instead of oscillating.
""")
note("N148", "Group velocity as a gradient", r"""
The group velocity is the gradient of ω in wavenumber space (note N71), component by component, next to the phase velocity
$\mathbf c=(\omega/K)\,\mathbf e_K$ *(7.8)*. **Gradient in wavenumber space**: treat (k, l, m) as coordinates and take
partial derivatives (Ch. 1 P25).
""", equation=EQ["7.143"], ref="7.143")
remind([
    ("gradient in wavenumber space", r"treat (k, l, m) as coordinates: $\nabla_{\mathbf K}\omega=(\partial\omega/\partial k,\partial\omega/\partial l,\partial\omega/\partial m)$, partial derivatives as in Ch. 1 P25."),
    ("quotient rule", r"$\frac{d}{dk}\frac{u}{v}=\frac{u'v-uv'}{v^2}$ — or the product rule on $uv^{-1}$ (Ch. 1 P38), as D36 step 3 does."),
])
D("D36", ref="7.146")
note("N149", "Phase and group velocity of internal waves", r"""
$\mathbf c=\frac\omega{K^2}(k\mathbf e_x+m\mathbf e_z)$ *(7.144)* and $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$
*(7.145)* for k > 0; for k < 0, as drawn in the book's Fig. 7.29, the code uses $\nabla_{\mathbf K}$ of $N\lvert k\rvert/K$
(D36 step 7).
""")
nb.figure(r"""
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.8))       # the book's geometry | the right triangle
v = ch07.internal_wave_velocities(-0.7071, 0.7071, 1.0)    # k < 0: K up-left at 45°, N = 1 rad/s
for vec, cl, lab in ((v["c"], COLORS["orange"], "c"), (v["cg"], COLORS["accent"], "c_g")):   # the two velocities
    a1.annotate("", xy=vec, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=cl, lw=2.5))   # arrow from the origin
    a1.text(vec[0]*1.15, vec[1]*1.15, lab, color=cl, fontsize=11)   # label
s_ = np.linspace(-1, 1, 2)                                 # crest lines ⟂ K through a few points
for off in np.linspace(-0.8, 0.8, 5):   # five crest lines
    a1.plot(off*(-0.7071) + s_*0.7071, off*0.7071 + s_*0.7071, color=COLORS["muted"], lw=0.8)   # direction ⟂ K = (1, 1)/√2
a1.annotate("", xy=(0.35, 0.35), xytext=(-0.35, -0.35), arrowprops=dict(arrowstyle="<->", color=COLORS["teal"], lw=1.5))   # motion along the crests
a1.text(0.38, 0.28, "particle motion", color=COLORS["teal"], fontsize=8)   # label
a1.annotate("", xy=(-0.7071*0.9, 0.7071*0.9), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["ink"], lw=1))   # K (up-left)
a1.text(-0.85, 0.72, "K", fontsize=10)   # label
a1.set_xlim(-1, 1); a1.set_ylim(-1, 1); a1.set_aspect("equal"); a1.set_title("k < 0: c up-left, c_g down-left", fontsize=9)   # true angles
a1.set_xlabel("x"); a1.set_ylabel("z")   # axes
for r, cl in ((0.3, "#c4b5fd"), (0.71, "#8b5cf6"), (0.9, "#4c1d95")):   # ω/N: the right triangles
    st = ch07.internal_wave_state(r, 1.0)                   # k > 0, m > 0
    a2.plot([0, st["cx"], st["cx"] + st["cgx"], 0], [0, st["cz"], st["cz"] + st["cgz"], 0], color=cl, lw=2,   # c, then c_g, then back
            label=f"ω/N = {r}: c ⊥ c_g")   # legend
a2.axhline(0, color=COLORS["muted"], lw=0.5); a2.set_aspect("equal"); a2.legend(fontsize=7, loc="lower left")   # horizontal reference
a2.set_xlabel("x-velocity [m/s]"); a2.set_ylabel("z-velocity [m/s]")   # axes
a2.set_title("c then c_g: legs of a right triangle, horizontal hypotenuse", fontsize=9)   # the message
fig.suptitle("Phase up, energy down — at right angles", fontweight="bold")   # overall message
savefig(fig, "ch07", "c16_geometry"); plt.show()   # save, then draw
""", see=r"""Left, the book's geometry with k < 0 (our remake of Figs. 7.29 and 7.31 `N195`, `N197`): K and c up-left, $c_g$
down-left at a right angle, the crest lines (grey) and the particle motion (teal) along them. Right, for three frequencies,
c followed by $c_g$: the two legs of a right triangle whose closing side is horizontal.""",
    read=r"""The hypotenuse is horizontal: $c_x$ and $c_{gx}$ have the same sign, $c_z=-c_{gz}$.""",
    change=r"""…ω → N: θ → 0, K horizontal, $c_g\to0$ — vertical columns bob and nothing propagates.""")
note("N150", "St Andrew's cross", r"""
A source oscillating at ω < N can only radiate waves with $\cos\theta=\omega/N$, so its energy leaves along four beams at θ
from the vertical (for example ω = 0.71N gives θ = arccos 0.71 ≈ 44.8°, the case of our figure), phase lines moving
across each beam. The figure is **our
superposition of four Gaussian beams — an illustration of the geometry, not the photograph** (`N199`).
""")
nb.figure(r"""
n = 400 if not FAST else 200                               # grid size
xg = np.linspace(-8, 8, n)                                 # [m]
det = ch07.st_andrews_cross(xg, xg, 0.0, omega=0.71, N=1.0, width=1.0, detail=True)   # the illustrative field and its beams
fig, ax = plt.subplots(figsize=(5.5, 5.5))   # square panel
ax.imshow(det["field"], extent=(-8, 8, -8, 8), origin="lower", cmap="RdBu_r", vmin=-1, vmax=1)   # the field: crests red, troughs blue
beam_labels(ax, det["theta"], r=7.5)                        # dashed beam axes at θ from the vertical
for bm in det["beams"]:                                     # c (orange) and c_g (purple) on each beam
    p0 = 4.5*bm["e_beam"]   # a point on the beam axis
    ax.annotate("", xy=p0 + 2.5*bm["cg"]/np.linalg.norm(bm["cg"]), xytext=p0, arrowprops=dict(arrowstyle="->", color=COLORS["accent"], lw=2))   # c_g: along the beam
    ax.annotate("", xy=p0 + 2.5*bm["c"]/np.linalg.norm(bm["c"]), xytext=p0, arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2))   # c: across the beam
ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")   # axes
ax.set_title(f"ω/N = 0.71: beams {np.degrees(det['theta']):.1f}° from the vertical, {90 - np.degrees(det['theta']):.1f}° from the horizontal", fontsize=9)   # both angles
savefig(fig, "ch07", "c16_st_andrews"); plt.show()   # save, then draw
""", see=r"""An X of striped beams (red and blue: crests and troughs of the density perturbation) with the purple $c_g$ arrows
pointing away from the source along each beam and the orange c arrows across it.""",
    read=r"""The stripes are phase lines (constant ρ′); energy runs along the beam, phase across it.""",
    change=r"""…ω/N = 0.3: the beams tilt toward the horizontal (72.5° from the vertical).""")
note("N151", "Real oceans", r"""
In the real ocean N depends on z (below 0.01 rad/s); the plane-wave results hold locally where N changes little over a
vertical wavelength 2π/m, and rays (C10's method) bend as N changes — WKB theory in Ch. 13.
""")
note("N152", "The energy equation", r"""
Multiply (7.128)–(7.130) by ρ₀u, ρ₀v, ρ₀w and add — the Ch. 4 energy move,
$\rho\frac{D}{Dt}\big(\tfrac12u_j^2\big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}$
*(4.56)* — to get (7.147): kinetic energy changes by buoyancy work gρ′w and by the divergence of the energy flux p′u. The
code evaluates it point by point for a plane wave.
""", equation=EQ["7.147"], ref="7.147")
nb.code(r"""
print(ch07.internal_energy_budget_residual(0.3, -0.2, 1.0, 1.0, 1.0, 0.01, 1e-3), "W/m³")   # (7.147): ≈ 0
""")
note("N153–N155", "Available potential energy", r"""
With (7.131) the buoyancy work is a time derivative, $\frac{\partial E_p}{\partial t}=g\rho'w=\frac{\partial}{\partial t}
\Big[\frac{g^2\rho'^2}{2\rho_0N^2}\Big]$ *(7.148)*; with w = ∂ζ/∂t, $\rho'=\frac{N^2\rho_0\zeta}{g}$ *(7.149)*, so the energy
per unit **volume** is (7.150) — the available potential energy of Ch. 13. Number: N = 0.01 rad/s, ζ = 10 m → 5 J/m³.
""", equation=EQ["7.150"], ref="7.150")
remind([
    ("density jump as a delta function", r"a step of height Δρ has a derivative that is zero except at the step, with total Δρ: $-d\bar\rho/dz=\Delta\rho\,\delta(z)$ — the 1-D version of Ch. 6 P149's delta. The next note uses it to write $N^2=\frac g{\rho_0}(\rho_2-\rho_1)\delta(z)$ *(7.152)* and recover $E_p=\tfrac14(\rho_2-\rho_1)ga^2$ *(7.151)*."),
])
note("N156, N157", "Consistency with C13", r"""
For two deep fluids the potential energy per unit area was $E_p=\tfrac14(\rho_2-\rho_1)ga^2$ *(7.151)*; writing a density
jump as $N^2=\frac g{\rho_0}(\rho_2-\rho_1)\delta(z)$ *(7.152)* and integrating (7.150) over z recovers it (a 1-D Dirac
delta is the limit of a narrow step's derivative — Ch. 6 P149 in one dimension). The code smooths the jump over a width ε
and lets ε shrink:
""")
nb.code(r"""
eps_list = [1.0, 0.3, 0.1, 0.03]                           # widths of the smoothed density step [m]
errs = []                                                   # relative errors against ¼Δρga²
for eps in eps_list:
    d = ch07.internal_pe_interface_limit(1.0, 1000.0, 1002.0, eps=eps, detail=True)   # ∫½N²ρ₀⟨ζ²⟩dz [J/m²]
    errs.append(abs(d["rel_error"]))
    print(f"ε = {eps:4.2f} m: {d['Ep_column']:.5f} J/m² vs ¼Δρga² = {d['target']:.4f} ({d['rel_error']:+.2%})")
print(f"observed order in ε: {observed_order(eps_list, errs):.2f}")   # ≈ 1: the error shrinks like ε
""")
remind([
    ("complex conjugate", r"$z^*=a-ib$ for $z=a+ib$; $zz^*=\lvert z\rvert^2$ (Ch. 2 P81)."),
])
P("P178", "mean of a product of real parts", r"""
For two fields written as real parts of complex amplitudes, the average over a period of their product is
$\langle\mathrm{Re}(Ae^{i\theta})\,\mathrm{Re}(Be^{i\theta})\rangle=\tfrac12\mathrm{Re}(AB^*)$ (B* the complex conjugate,
Ch. 2 P81). Two fields 90° out of phase (B = iA) give zero: they carry no mean product.
""", code=r"""
import numpy as np                                # numbers
A, B = 2.0, 3.0*np.exp(1j*0.5)                    # amplitudes with a phase shift
th = np.linspace(0, 2*np.pi, 100001)              # one period of the phase
mean = np.mean(np.real(A*np.exp(1j*th))*np.real(B*np.exp(1j*th)))   # the average of the product of real parts
print(mean, 0.5*np.real(A*np.conj(B)))            # both 2.633
""")
note("N158", "The plane-wave fields", r"""
Assume $[u,w,p',\rho']=[\hat u,\hat w,\hat p,\hat\rho]e^{i(kx+mz-\omega t)}$ with real ŵ (D37 step 1).
""")
note("N160–N163", "The energy of an internal wave (used in D37 step 11)", r"""
Per unit volume and averaged over a period, the kinetic energy is $\tfrac12\rho_0\langle u^2+w^2\rangle$ and the potential
energy is $\tfrac12\frac{g^2}{\rho_0N^2}\langle\rho'^2\rangle$ (the parcel-displacement energy of Ch. 1's N, written with ρ′).
With the ½Re(AB*) rule (P178 above) each mean square is half the squared amplitude, and the amplitudes are the plane-wave
ones D37 finds in steps 3–4 (derived line by line just below): $\lvert\hat u\rvert=\frac mk\hat w$ and $\lvert\hat\rho\rvert=\frac{N^2\rho_0}{\omega g}\hat w$. That gives
$E_k=\tfrac14\rho_0\big(\frac{m^2}{k^2}+1\big)\hat w^2$ *(7.154)* and $E_p=\frac{N^2\rho_0}{4\omega^2}\hat w^2$ *(7.155)*.
They are equal, $E_k=E_p$ *(7.156)*, because the dispersion relation $\omega^2=k^2N^2/(k^2+m^2)$ turns $N^2/\omega^2$ into
$(k^2+m^2)/k^2=\frac{m^2}{k^2}+1$. So the total is $E=\tfrac12\rho_0\big(\frac{m^2}{k^2}+1\big)\hat w^2$ *(7.157)* — the E
that D37 multiplies by $\mathbf c_g$ in step 11. (Stated, not derived line by line; the code after D37 checks
$E_k=E_p$ to round-off.)
""", equation=EQ["7.157"], ref="7.157")
D("D37", ref="7.159")
note("N159, N164–N166", "Polarization and flux", r"""
The polarization relations $p'=-\frac{\omega m\rho_0}{k^2}\hat we^{i\theta},\ \rho'=\frac{iN^2\rho_0}{\omega g}\hat we^{i\theta},\ u=-\frac mk\hat we^{i\theta}$
*(7.153)*: ρ′ is 90° out of phase with w (the factor i), so ρ′ peaks a quarter period after the upward velocity does. The mean
flux $\mathbf F=\frac{\rho_0\omega m\hat w^2}{2k^2}\big(\mathbf e_x\frac mk-\mathbf e_z\big)$ *(7.158)* and $\mathbf c_gE$
reduce to each other with $\omega=kN/K$ *(7.138)*: $\mathbf F=\mathbf c_gE$ *(7.159)* — the same law as $F=Ec_g$ *(7.71)*.
""")
confusion(r"""
**units**: in $\mathbf F=\mathbf c_gE$ *(7.159)* F is per unit **area** and E per unit **volume**; in $F=Ec_g$ *(7.71)* F was
per metre of crest (the whole depth) and E per unit horizontal area. Same law, different bookkeeping.
""")
nb.worked_example("c and c_g at 45°", r"""
N = 1 rad/s, K = 1 rad/m at θ = 45° (k = m = 0.707 rad/m).

1. ω = N cos 45° = 0.707 rad/s.
2. $\mathbf c=\frac{\omega}{K^2}(k,m)=0.707\times(0.707,0.707)=(0.5,0.5)$ m/s: up and to the right.
3. $\mathbf c_g=\frac{Nm}{K^3}(m,-k)=0.707\times(0.707,-0.707)=(0.5,-0.5)$ m/s: down and to the right.
4. $\mathbf c\cdot\mathbf c_g=0.25-0.25=0$: a right angle; $\lvert\mathbf c\rvert=\lvert\mathbf c_g\rvert=0.707$ m/s at this
   angle only.
5. ω/N = 0.71 → θ = arccos 0.71 = 44.8°.
""")
nb.code(r"""
print(ch07.internal_wave_velocities(0.7071, 0.7071, 1.0))  # k > 0: c = (0.5, 0.5), c_g = (0.5, −0.5), dot 0
print(ch07.internal_wave_velocities(-0.7071, 0.7071, 1.0)) # k < 0 (the book's Fig. 7.29): c up-left, c_g down-left
print("printed (7.145) for k < 0:", ch07.internal_wave_velocities(-0.7071, 0.7071, 1.0, printed=True)["cg"])   # wrong way
e = ch07.internal_wave_energy(0.7071, 0.7071, 1.0, 0.01)   # ŵ = 1 cm/s, ρ₀ = 1000 kg/m³
print({n: v for n, v in e.items() if n in ("Ek", "Ep", "E", "F", "cgE")})
assert np.allclose(e["F"], e["cgE"]) and np.isclose(e["Ek"], e["Ep"])   # (7.159) and (7.156)
st = ch07.internal_wave_state(0.71, 1.0)
print(f"ω/N = 0.71: beam {st['beam_from_vertical_deg']:.2f}° from the vertical, {st['beam_from_horizontal_deg']:.2f}° from the horizontal")
""", explain=r"""
1. For k > 0, c = (0.5, 0.5) and $\mathbf c_g=(0.5,-0.5)$ m/s: perpendicular.
2. For k < 0 (K up-left, as the book draws it) c = (−0.5, 0.5) and the sign-safe $\mathbf c_g=\nabla_{\mathbf K}\omega=
   (-0.5,-0.5)$ m/s — down-left, as drawn. The printed (7.145) would give (0.5, 0.5), up-right: wrong for k < 0.
3. $E_k=E_p=0.05$ J/m³ for ŵ = 1 cm/s, and the mean flux F equals $\mathbf c_gE$ component by component.
4. The beam of a source at ω = 0.71N makes 44.77° with the vertical and 45.23° with the horizontal.
""")
nb.check_agree(r"""
wK = lambda k_, m_: ch07.internal_wave_omega(k_, m_, 1.0)  # ω(k, m) = N|k|/K with N = 1 rad/s
k0, m0, h = -0.7071, 0.7071, 1e-6                          # the book's k < 0 case; difference step [rad/m]
cg_mine = np.array([(wK(k0 + h, m0) - wK(k0 - h, m0))/(2*h),   # ∂ω/∂k
                    (wK(k0, m0 + h) - wK(k0, m0 - h))/(2*h)])  # ∂ω/∂m
c = ch07.internal_wave_velocities(k0, m0, 1.0)["c"]         # phase velocity (7.144)
print("our gradient:", cg_mine, " c·c_g =", np.dot(cg_mine, c))
assert np.allclose(cg_mine, ch07.internal_wave_velocities(k0, m0, 1.0)["cg"], rtol=1e-7)   # c_g = ∇_K ω
assert abs(np.dot(cg_mine, c)) < 1e-9                      # perpendicular to c (to the difference error)
""")
nb.md(r"**What does this show?** Our centred differences of ω(k, m) at k = −0.7071, m = 0.7071 rad/m (N = 1 rad/s) give $\mathbf c_g=\nabla_{\mathbf K}\,\omega=(-0.5,-0.5)$ m/s — energy heading down and to the left — the same vector as the sign-safe `internal_wave_velocities`, and $\mathbf c\cdot\mathbf c_g=0$ to round-off: the energy moves at right angles to the phase.")
nb.animation(r"""
n = 256 if not FAST else 128                               # doubly periodic grid
xg = np.linspace(0, 100, n, endpoint=False); zg = np.linspace(0, 100, n, endpoint=False)   # [m]
X, Z = np.meshgrid(xg, zg)                                 # (z, x) layout
k0 = m0 = 1/np.sqrt(2)                                     # K₀ = 1 rad/m at 45° up-right
env0 = np.exp(-((X - 35)**2 + (Z - 65)**2)/(2*8.0**2))     # Gaussian envelope, 8 m wide, centred at (35, 65) m
field0 = env0*np.exp(1j*(k0*X + m0*Z))                     # the analytic packet (complex)
nf = 40 if not FAST else 20                                # frames
tt = np.linspace(0, 60, nf)                                # times [s]
frames = ch07.linear_evolve_2d(field0, xg, zg, tt, lambda k_, m_: np.abs(k_)/np.hypot(k_, m_))   # each mode at N|k|/K, N = 1
v = ch07.internal_wave_velocities(k0, m0, 1.0)             # c = (0.5, 0.5), c_g = (0.5, −0.5) m/s
fig, ax = plt.subplots(figsize=(5.2, 5.2))
im = ax.imshow(frames[0], extent=(0, 100, 0, 100), origin="lower", cmap="RdBu_r", vmin=-1, vmax=1)
ax.annotate("", xy=(35 + 16*v["cg"][0]/0.707, 65 + 16*v["cg"][1]/0.707), xytext=(35, 65), arrowprops=dict(arrowstyle="->", color=COLORS["accent"], lw=2))
ax.annotate("", xy=(35 + 16*v["c"][0]/0.707, 65 + 16*v["c"][1]/0.707), xytext=(35, 65), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2))
ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")   # axes


def update(i):   # draw frame i
    im.set_data(frames[i])                                 # the packet at time tt[i]
    ax.set_title(f"t = {tt[i]:.0f} s: phase ↗ (orange c), packet ↘ (purple c_g)", fontsize=9)
    return []   # nothing to blit


show_animation(animate(update, frames=nf, fig=fig, interval=70))
""")
see_read_change(r"""A packet of internal waves (N = 1 rad/s, K at 45° up-right): the stripes (crests) move up-right while the whole
packet slides down-right, along its own crests — the `A6` animation, our version of the book's Fig. 7.32 `N198`.""",
                r"""The packet moves at $\mathbf c_g=(0.5,-0.5)$ m/s: after 60 s its centre has moved 30 m right and 30 m down,
while each stripe moved 30 m right and 30 m up.""",
                r"""…ω → N (K nearly horizontal): $c_g\to0$ — the packet would hardly move while its stripes still slide.""")
explainer("internal_wave_beams", "Why does energy leave at right angles to the crests?",
          r"""The ω/N slider tilts the beams while the crests visibly slide across them; drag K in the wavenumber plane: ω
stays fixed when only |K| changes, and $c_g$ stays perpendicular to c.""", "",
          ["Preset *45° cross*: read both angles in the status.",
           "Drag the K arrow outward along its ray: ω does not change.",
           "Choose *near N*: the beams stand almost vertical and the energy hardly moves.",
           "Set the book's Fig. 7.29 geometry (k < 0) and check that $c_g$ points down-left."])
whatif(r"""
…N varied with depth, as in the real ocean (N largest in the thermocline)? The beams curve, like the rays of C10, and can
reflect where ω = N — internal waves get trapped in the thermocline (Ch. 13).
""")
nb.pointer("**Exercises 7.1–7.20** are not reproduced (S01). Ideas from them used in this notebook: 7.2 → C12 (Stokes "
           "expansion and its trap); 7.3 → N02 (any initial shape); 7.4 → N30 (ψ); 7.5–7.6 → N63 (basin modes); 7.9–7.10 → "
           "N69, N72 (capillary c_g, minimum c_g); 7.11 → N72 (viscous decay); 7.12 → D20 (Gaussian packet); 7.13–7.14 → "
           "N91 (drift at any depth); 7.15 → N97 (soliton check); 7.18 → N105 (interfacial energies); 7.19 → D30 (two-layer "
           "dispersion); 7.20 → N131 (rigid lids).")
nb.summary(
    clicked=[
        r"**C01** A wave is a moving phase: crests move at $c=\omega/k=\lambda\nu$ *(7.4)* while the water only oscillates.",
        r"**C02** For gentle waves the surface conditions can be moved to the flat level z = 0; what is dropped is of order $(ka)\cdot a\omega$ — relative error O(ka).",
        r"**C03** Laplace plus the bottom, kinematic and dynamic conditions leave one frequency per wavenumber: $\omega=\sqrt{gk\tanh kH}$ *(7.28)*.",
        r"**C04** Long waves outrun short ones; depth matters only for waves longer than about three times the depth, and all long waves travel at $\sqrt{gH}$ *(7.49)*.",
        r"**C05** Under a wave the water goes round closed ellipses — circles in deep water, flat ellipses in shallow water.",
        r"**C06** A wave holds $E=\tfrac12\rho ga^2$ *(7.42)* per square metre, half kinetic and half potential, and carries it at a speed that is not the crest speed.",
        r"**C07** Surface tension adds a second restoring force; together with gravity it makes a slowest wave, 23 cm/s at 1.7 cm.",
        r"**C08** Two opposite waves make a standing wave; walls keep only the patterns whose nodes fit, so a lake rings at its own periods.",
        r"**C09** Groups and energy travel at $c_g=d\omega/dk$ *(7.67)* — half the crest speed in deep water — and $F=Ec_g$ *(7.71)*.",
        r"**C10** Along a ray moving at $c_g$ the frequency never changes; the wavelength shrinks in shallow water, so crests turn parallel to the shore.",
        r"**C11** A hydraulic jump's height comes from momentum alone; energy is lost, which is why only fast shallow flow can jump — or dispersion balances steepening in a soliton.",
        r"**C12** At finite amplitude the orbits do not close: the water drifts forward at $a^2\omega ke^{2kz_0}$ while the mean at a fixed point stays zero; Stokes' speed correction needs a consistent third-order expansion.",
        r"**C13** A small density step gives slow, tall waves, $\omega=\varepsilon\sqrt{gk}$ *(7.95)*, with the two fluids sliding past each other.",
        r"**C14** A layer over deep water rings in two modes: barotropic (fast, at the surface) and baroclinic (slow, $\sqrt{g'H}$, almost invisible at the surface).",
        r"**C15** In a stratified fluid the frequency is set by direction only: $\omega=N\cos\theta$ *(7.139)*, never above N.",
        r"**C16** Internal waves are transverse and their energy moves at right angles to the crests: phase up, energy down, $\mathbf F=\mathbf c_gE$ *(7.159)*.",
    ],
    feeds_forward=[
        "Ch. 11: the interfacial vortex sheet with a shear flow is the Kelvin–Helmholtz instability; the σ-term of (7.56) stabilises short waves.",
        "Ch. 13: √(gH) and hydrostatic shallow water carry tides, tsunamis and Kelvin waves; rotation adds f to ω² (Poincaré and inertia–gravity waves); √(g′H)/f is the baroclinic Rossby radius; rays and WKB follow waves through N(z); equipartition fails in geostrophic adjustment.",
        "Ch. 15: characteristics and shocks are the compressible cousins of the simple wave (note N81) and the hydraulic jump.",
    ],
    left_out=[
        "Exercise text (S01).",
        "Nonlinear internal waves and wave breaking (research literature).",
        "Viscous damping in detail (Ch. 8 tools; named in note N72).",
    ],
)

# ---------------------------------------------------------------------------------------------------------------------
# final pass: every equation named by number is written out; save (nbkit's coverage checks run here)
# ---------------------------------------------------------------------------------------------------------------------
if "--partial" not in sys.argv:
    n_changed = finalize_equations()
    n_tex = tidy_raw_tex()
    print(f"plain-text exponents turned into maths in {n_tex} cells")
    bad = self_check_numbers() + self_check_prose()
    if bad:
        raise SystemExit("markdown cells cite equation numbers without maths:\n  " + "\n  ".join(bad))
    out = nb.save()
    print(f"wrote {out.relative_to(ROOT)} ({len(nb.cells)} cells; equations written out in {n_changed} cells)")


if "--partial" in sys.argv:                                  # development aid: write what exists so far, unchecked
    import nbformat as _nbf
    finalize_equations()
    for _m in self_check_numbers() + self_check_prose():
        print("CHECK", _m)
    _nb2 = _nbf.v4.new_notebook(cells=list(nb.cells))
    _nb2.metadata["kernelspec"] = {"name": "fluidpy-venv", "display_name": "fluidpy-venv", "language": "python"}
    _out = ROOT / "outputs" / "ch07" / "partial.ipynb"
    _nbf.write(_nb2, _out)
    print("wrote", _out, len(nb.cells), "cells")
    sys.exit(0)
