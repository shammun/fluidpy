"""Build the Chapter 11 teaching notebook: ``notebooks/ch11_instability.ipynb``.

Source of truth: ``analysis/ch11_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract —
``fluidpy.ch11_instability`` imported as ``ch11`` re-exports the eigen-solver toolkit ``core.stability`` = ``ST``), Part D
(runtime budget and FAST sizes), Part E (prerequisite ledger → primers P255–P279 and one-line reminders of earlier primers), Part F
(the 25 derivations D01–D25, one move per step) and ``analysis/ch11_curation.md`` (IDs, depths, section coverage), with the
facts of ``reports/ch11_verification.md`` written into the text:

* heavy results (the Taylor–Goldstein growth map, neutral curves, critical points, the Taylor table) are read from the public
  tables in ``reference/ch11/`` through the cached paths of ``ch11``; one live recomputation per table shows the method running;
* ``tg_growth`` / ``tg_growth_map`` use the box ``decay_box(k)`` and the node rule ``decay_map_scale(k)``; a 0 within 0.006 (in J)
  of the neutral curve J = k(1 − k) means kc_i < 0.004, not "stable" — said where the map is shown;
* the integral identities (11.65), (11.69), (11.70) are checked with the core solver at map scale 0.5, away from the neutral curve;
* ``benard_marginal_Ra_det`` is called inside its scan range only; 4-s.f. cached tables are compared with rtol ≥ 6e-4;
* numbers the design computed with scratch scripts and that the tested functions give slightly differently (the double-diffusive
  cubic roots, the finger layer with g = 9.80665) are printed by the cells, never typed into prose;
* the lapse-rate sign: code takes ``dT`` = T_bottom − T_top; Kundu's Γ ≡ dT/dz and the meteorological Γ ≡ −dT/dz are both shown;
* the 12 printed slips are taught in corrected form with a short box (``slip #k``).

**Derivations are read from Part F at build time** (``part_f()`` below, the ch07–ch10 parser): goal, start, plan, tools,
assumptions, every step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word; Starts and Results
that Part F writes only as equation numbers are rewritten with the equations shown; the ★★★ D10, D14, D18, D19, D24 carry the sympy
check cells of Part F, and D03, D08, D11, D21, D23 optional ones.

Book numbers never printed (rule 9): our own inputs everywhere (design convention 11); published benchmarks are read from
``reference/ch11/benchmarks.json`` and cited.

Run:  .venv/Scripts/python.exe notebooks/build_ch11.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch11.py --dump     (prints the parsed Part F derivations only)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch11")
DESIGN = (ROOT / "analysis" / "ch11_design.md").read_text(encoding="utf-8")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# Ch. 11 from the design (each re-read on the rendered pages chapters/pages/ch11, list in the design header; (11.18),
# (11.45)–(11.46), (11.79)–(11.80), (11.88), (11.91) re-read again for this builder), in the CORRECTED form wherever the
# book prints a slip (#1 … #12).  Earlier chapters' equations as their own notebooks show them.
# ---------------------------------------------------------------------------------------------------------------------
_D2 = r"\frac{d^2}{dz^2}"
EQ: dict[str, str] = {
    "11.1": r"u=\hat u(z)\exp\{\mathrm ikx+\mathrm imy+\sigma t\}=\hat u(z)\exp\{\mathrm i\lvert\mathbf K\rvert(\mathbf e_K\cdot\mathbf x-ct)\}",
    "11.2": r"\tilde\phi_1=U_1x+\phi_1,\ \tilde\phi_2=U_2x+\phi_2",
    "11.3": r"\nabla^2\phi_1=0,\ \nabla^2\phi_2=0",
    "11.4": r"\phi_1\to0\ \text{as}\ z\to+\infty",
    "11.5": r"\phi_2\to0\ \text{as}\ z\to-\infty",
    "11.6": r"\mathbf n\cdot\nabla\tilde\phi_1=\mathbf n\cdot\mathbf U_s=\mathbf n\cdot\nabla\tilde\phi_2\ \text{on}\ z=\zeta",
    "11.7": r"p_1=p_2\ \text{on}\ z=\zeta",
    "11.8": r"\mathbf n\cdot\Big\{\frac{\partial\tilde\phi_1}{\partial x}\mathbf e_x+\frac{\partial\tilde\phi_1}{\partial z}\mathbf e_z\Big\}=\mathbf n\cdot\Big\{\frac{\partial\zeta}{\partial t}\mathbf e_z\Big\}=\mathbf n\cdot\Big\{\frac{\partial\tilde\phi_2}{\partial x}\mathbf e_x+\frac{\partial\tilde\phi_2}{\partial z}\mathbf e_z\Big\}",
    "11.9": r"-U_1\frac{\partial\zeta}{\partial x}+\frac{\partial\phi_1}{\partial z}=\frac{\partial\zeta}{\partial t}=-U_2\frac{\partial\zeta}{\partial x}+\frac{\partial\phi_2}{\partial z}\ \text{on}\ z=0",
    "11.10": r"\frac{\partial\tilde\phi_j}{\partial t}+\tfrac12\lvert\nabla\tilde\phi_j\rvert^2+\frac{p_j}{\rho_j}+gz=C_j",
    "11.11": r"\rho_1\big(C_1-\partial_t\tilde\phi_1-\tfrac12\lvert\nabla\tilde\phi_1\rvert^2-gz\big)=\rho_2\big(C_2-\partial_t\tilde\phi_2-\tfrac12\lvert\nabla\tilde\phi_2\rvert^2-gz\big)\ \text{on}\ z=\zeta",
    "11.12": r"\rho_1(C_1-\tfrac12U_1^2)=\rho_2(C_2-\tfrac12U_2^2)",
    "11.13": r"\rho_1\Big(\frac{\partial\phi_1}{\partial t}+U_1\frac{\partial\phi_1}{\partial x}+g\zeta\Big)=\rho_2\Big(\frac{\partial\phi_2}{\partial t}+U_2\frac{\partial\phi_2}{\partial x}+g\zeta\Big)\ \text{on}\ z=0",
    "11.14": r"\phi_j=A_j(z)\exp\{\mathrm ik(x-ct)\}",
    "11.15": r"\phi_1=A_-\exp\{\mathrm ik(x-ct)-kz\},\ \phi_2=A_+\exp\{\mathrm ik(x-ct)+kz\}",
    "11.16": r"-\mathrm iU_1k\zeta_o-kA_-=-\mathrm ikc\zeta_o=-\mathrm iU_2k\zeta_o+kA_+",
    "11.17": r"\rho_1(-\mathrm ikcA_-+\mathrm ikU_1A_-+g\zeta_o)=\rho_2(-\mathrm ikcA_++\mathrm ikU_2A_++g\zeta_o)",
    "11.18": r"c=\frac{\rho_2U_2+\rho_1U_1}{\rho_2+\rho_1}\pm\Big[\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk-\frac{\rho_2\rho_1}{(\rho_2+\rho_1)^2}(U_2-U_1)^2\Big]^{1/2}",
    "11.19": r"c=\pm\Big[\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk\Big]^{1/2}",
    "11.20": r"c=\frac{U_2+U_1}2\pm\mathrm i\frac{U_2-U_1}2",
    "11.21": r"\mathrm{Ra}=\frac{g\alpha\Gamma d^4}{\kappa\nu}",
    "11.22": r"\tilde{\mathbf u}=0+\mathbf u(\mathbf x,t),\ \tilde T=\bar T(z)+T'(\mathbf x,t),\ \tilde p=P(z)+p(\mathbf x,t)",
    "11.23": r"0=-\frac1{\rho_0}\nabla P-g[1-\alpha(\bar T-T_0)]\mathbf e_z,\ \ 0=\kappa\frac{\partial^2\bar T}{\partial z^2}",
    "11.24": r"\bar T(z)=T_0-\tfrac12\Delta T-\Gamma z,\ \ \Gamma\equiv\Delta T/d",
    "11.25": r"\nabla\cdot\mathbf u=0",
    "11.26": r"\frac{\partial\mathbf u}{\partial t}=-\frac1{\rho_0}\nabla p+g\alpha T'\mathbf e_z+\nu\nabla^2\mathbf u",
    "11.27": r"\frac{\partial T'}{\partial t}-w\Gamma=\kappa\nabla^2T'",
    "11.28": r"\frac{\partial}{\partial t}\nabla^2w=-\frac1{\rho_0}\nabla^2\frac{\partial p}{\partial z}+g\alpha\nabla^2T'+\nu\nabla^4w",
    "11.29": r"\frac{\partial}{\partial t}\nabla^2w=+g\alpha\nabla_H^2T'+\nu\nabla^4w",
    "11.30": r"w=\frac{\partial w}{\partial z}=T'=0\ \text{on}\ z=\pm d/2",
    "11.31": r"\Big(\frac{\partial}{\partial t}-\nabla^2\Big)T'=\frac{\Gamma d^2}{\kappa}w",
    "11.32": r"\Big(\frac1\Pr\frac{\partial}{\partial t}-\nabla^2\Big)\nabla^2w=\frac{g\alpha d^2}\nu\nabla_H^2T'",
    "11.33": r"w=\frac{\partial w}{\partial z}=T'=0\ \text{on}\ z=\pm\tfrac12",
    "11.34": r"\Big(\sigma+K^2-" + _D2 + r"\Big)\hat T=\frac{\Gamma d^2}\kappa\hat w",
    "11.35": r"\Big(\frac\sigma\Pr+K^2-" + _D2 + r"\Big)\Big(" + _D2 + r"-K^2\Big)\hat w=-\frac{g\alpha d^2K^2}\nu\hat T",
    "11.36": r"\Big(\sigma+K^2-" + _D2 + r"\Big)\hat T=W",
    "11.37": r"\Big(\frac\sigma\Pr+K^2-" + _D2 + r"\Big)\Big(" + _D2 + r"-K^2\Big)W=-\mathrm{Ra}K^2\hat T",
    "11.38": r"W=\frac{dW}{dz}=\hat T=0\ \text{on}\ z=\pm\tfrac12",
    "11.39": r"\Big(" + _D2 + r"-K^2\Big)\hat T=-W,\ \ \Big(" + _D2 + r"-K^2\Big)^2W=\mathrm{Ra}K^2\hat T",
    "11.40": r"\Big(" + _D2 + r"-K^2\Big)^3W=-\mathrm{Ra}K^2W",
    "11.41": r"W=\frac{dW}{dz}=\Big(" + _D2 + r"-K^2\Big)^2W=0\ \text{on}\ z=\pm\tfrac12",
    "11.42": r"q^2=-K^2\Big[\Big(\frac{\mathrm{Ra}}{K^4}\Big)^{1/3}-1\Big],\ \ q^2=K^2\Big[1+\frac12\Big(\frac{\mathrm{Ra}}{K^4}\Big)^{1/3}(1\pm\mathrm i\sqrt3)\Big]",
    "11.43": r"W=\frac{d^2W}{dz^2}=\frac{d^4W}{dz^4}=0\ \text{on}\ z=\pm\tfrac12",
    "11.44": r"\mathrm{Ra}=\frac{(n^2\pi^2+K^2)^3}{K^2}",
    "11.45": r"\Big(" + _D2 + r"-K^2\Big)\hat T=-W,\ \ \frac{\kappa_s}\kappa\Big(" + _D2 + r"-K^2\Big)\hat s=-W,\ \ \Big(" + _D2 + r"-K^2\Big)^2W=-\mathrm{Ra}K^2\hat T+\mathrm{Rs}'K^2\hat s",
    "11.46": r"\frac{gd^4}\nu\Big[\frac\beta{\kappa_s}\frac{dS}{dz}-\frac\alpha\kappa\frac{d\bar T}{dz}\Big]=657",
    "11.47": r"\frac{D\tilde u_R}{Dt}-\frac{\tilde u_\varphi^2}R=-\frac1\rho\frac{\partial\tilde p}{\partial R}+\nu\Big(\nabla^2\tilde u_R-\frac{\tilde u_R}{R^2}\Big),\ \ \frac{D\tilde u_\varphi}{Dt}+\frac{\tilde u_R\tilde u_\varphi}R=\nu\Big(\nabla^2\tilde u_\varphi-\frac{\tilde u_\varphi}{R^2}\Big),\ \ \frac{D\tilde u_z}{Dt}=-\frac1\rho\frac{\partial\tilde p}{\partial z}+\nu\nabla^2\tilde u_z,\ \ \frac1R\frac{\partial}{\partial R}(R\tilde u_R)+\frac{\partial\tilde u_z}{\partial z}=0",
    "11.48": r"\tilde{\mathbf u}=\mathbf U+\mathbf u,\ \tilde p=P+p",
    "11.49": r"U_R=U_z=0,\ \ U_\varphi=AR+\frac BR,\ \ \frac1\rho\frac{dP}{dR}=\frac{U_\varphi^2}R",
    "11.50": r"\frac{\partial u_R}{\partial t}-\frac{2U_\varphi u_\varphi}R=-\frac1\rho\frac{\partial p}{\partial R}+\nu\Big(\nabla^2u_R-\frac{u_R}{R^2}\Big),\ \ \frac{\partial u_\varphi}{\partial t}+\Big(\frac{dU_\varphi}{dR}+\frac{U_\varphi}R\Big)u_R=\nu\Big(\nabla^2u_\varphi-\frac{u_\varphi}{R^2}\Big),\ \ \frac{\partial u_z}{\partial t}=-\frac1\rho\frac{\partial p}{\partial z}+\nu\nabla^2u_z,\ \ \frac1R\frac{\partial(Ru_R)}{\partial R}+\frac{\partial u_z}{\partial z}=0",
    "11.51": r"\Big(\frac{d^2}{dR^2}-k^2-\sigma\Big)\Big(\frac{d^2}{dR^2}-k^2\Big)\hat u_R=(1+\alpha x)\hat u_\varphi,\ \ \Big(\frac{d^2}{dR^2}-k^2-\sigma\Big)\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R",
    "11.52": r"\mathrm{Ta}\equiv4\Big(\frac{\Omega_1R_1^2-\Omega_2R_2^2}{R_2^2-R_1^2}\Big)\frac{\Omega_1d^4}{\nu^2}",
    "11.53": r"\hat u_R=\frac{d\hat u_R}{dR}=\hat u_\varphi=0\ \text{at}\ x=0,1",
    "11.54": r"\mathrm{Ta}_{cr}=\frac{1708}{\tfrac12(1+\Omega_2/\Omega_1)}",
    "11.55": r"\frac{\partial u}{\partial t}+w\frac{dU}{dz}+U\frac{\partial u}{\partial x}=-\frac1{\rho_0}\frac{\partial p}{\partial x},\ \ \frac{\partial w}{\partial t}+U\frac{\partial w}{\partial x}=-\frac1{\rho_0}\frac{\partial p}{\partial z}-g\frac\rho{\rho_0}",
    "11.56": r"\frac{\partial\rho}{\partial t}+U\frac{\partial\rho}{\partial x}-\frac{\rho_0N^2w}g=0",
    "11.57": r"u=\frac{\partial\psi}{\partial z},\ \ w=-\frac{\partial\psi}{\partial x}",
    "11.58": r"(U-c)\hat\psi'-U'\hat\psi=-\frac{\hat p}{\rho_0}",
    "11.59": r"k^2(U-c)\hat\psi=-\frac{g\hat\rho}{\rho_0}-\frac{\hat p'}{\rho_0}",
    "11.60": r"(U-c)\hat\rho+\frac{\rho_0N^2}g\hat\psi=0",
    "11.61": r"(U-c)\Big(" + _D2 + r"-k^2\Big)\hat\psi-\frac{d^2U}{dz^2}\hat\psi+\frac{N^2}{U-c}\hat\psi=0",
    "11.62": r"\hat\psi(0)=\hat\psi(d)=0",
    "11.63": r"\phi\equiv\frac{\hat\psi}{(U-c)^{1/2}}",
    "11.64": r"\frac d{dz}\Big[(U-c)\frac{d\phi}{dz}\Big]-\Big\{k^2(U-c)+\frac12\frac{d^2U}{dz^2}+\frac{\tfrac14(dU/dz)^2-N^2}{U-c}\Big\}\phi=0",
    "11.65": r"\int\frac{N^2-\tfrac14(dU/dz)^2}{U-c}\lvert\phi\rvert^2dz=\int(U-c)\Big\{\Big\lvert\frac{d\phi}{dz}\Big\rvert^2+k^2\lvert\phi\rvert^2\Big\}dz+\frac12\int\frac{d^2U}{dz^2}\lvert\phi\rvert^2dz",
    "11.66": r"\mathrm{Ri}(z)\equiv\frac{N^2}{(dU/dz)^2}",
    "11.67": r"\mathrm{Ri}>\tfrac14\ \text{everywhere}\ \Rightarrow\ \text{stable}",
    "11.68": r"F\equiv\frac{\hat\psi}{U-c}",
    "11.69": r"\int\big[(U-c_r)^2-c_i^2\big]Q\,dz=\int N^2\lvert F\rvert^2dz",
    "11.70": r"c_i\int(U-c_r)Q\,dz=0",
    "11.71": r"U_{\min}<c_r<U_{\max}",
    "11.72": r"\int\big[U^2-c_r^2-c_i^2\big]Q\,dz>0",
    "11.73": r"\frac{\partial u}{\partial t}+(U+u)\frac{\partial}{\partial x}(U+u)+v\frac{\partial}{\partial y}(U+u)=-\frac{\partial}{\partial x}(P+p)+\frac1{\mathrm{Re}}\nabla^2(U+u)",
    "11.74": r"\frac{\partial u}{\partial t}+U\frac{\partial u}{\partial x}+v\frac{\partial U}{\partial y}=-\frac{\partial p}{\partial x}+\frac1{\mathrm{Re}}\nabla^2u",
    "11.75": r"\frac{\partial v}{\partial t}+U\frac{\partial v}{\partial x}=-\frac{\partial p}{\partial y}+\frac{\nabla^2v}{\mathrm{Re}},\ \ \frac{\partial w}{\partial t}+U\frac{\partial w}{\partial x}=-\frac{\partial p}{\partial z}+\frac{\nabla^2w}{\mathrm{Re}},\ \ \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}=0",
    "11.76": r"[\mathbf u,p]=[\hat{\mathbf u}(y),\hat p(y)]\exp\{\mathrm i(kx+mz-kct)\}",
    "11.77": r"\mathrm ik(U-c)\hat u+\hat vU'=-\mathrm ik\hat p+\frac{\hat u''-(k^2+m^2)\hat u}{\mathrm{Re}},\ \ \mathrm ik(U-c)\hat v=-\hat p'+\frac{\hat v''-(k^2+m^2)\hat v}{\mathrm{Re}},\ \ \mathrm ik(U-c)\hat w=-\mathrm im\hat p+\frac{\hat w''-(k^2+m^2)\hat w}{\mathrm{Re}},\ \ \mathrm ik\hat u+\hat v'+\mathrm im\hat w=0",
    "11.78": r"\bar k=\sqrt{k^2+m^2},\ \bar c=c,\ \bar k\bar u=k\hat u+m\hat w,\ \bar v=\hat v,\ \bar p/\bar k=\hat p/k,\ \bar k\,\overline{\mathrm{Re}}=k\,\mathrm{Re}",
    "11.79": r"(U-c)\Big(\frac{d^2\phi}{dy^2}-k^2\phi\Big)-\frac{d^2U}{dy^2}\phi=\frac1{\mathrm ik\mathrm{Re}}\Big(\frac{d^4\phi}{dy^4}-2k^2\frac{d^2\phi}{dy^2}+k^4\phi\Big)",
    "11.80": r"\phi=\frac{d\phi}{dy}=0\ \text{at}\ y=y_1,\ y_2",
    "11.81": r"(U-c)\Big(\frac{d^2\phi}{dy^2}-k^2\phi\Big)-\frac{d^2U}{dy^2}\phi=0",
    "11.82": r"\phi=0\ \text{at}\ y=y_1,\ y_2",
    "11.83": r"\int\big(\lvert\phi'\rvert^2+k^2\lvert\phi\rvert^2\big)dy+\int\frac1{U-c}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=0",
    "11.84": r"c_i\int\frac1{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=0",
    "11.85": r"\int\frac{U-c_r}{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=-\int\big(\lvert\phi'\rvert^2+k^2\lvert\phi\rvert^2\big)dy<0",
    "11.86": r"(c_r-U_I)\int\frac1{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=0",
    "11.87": r"\hat\psi=\int(U-c)\,dy+A\phi(y)\exp\{\mathrm ikx\}",
    "11.88": r"\frac d{dt}\int\tfrac12u_i^2\,dV=-\int u_iu_j\frac{\partial U_i}{\partial x_j}dV-\Lambda",
    "11.89": r"\dot X=Y,\ \ \dot Y=-(g/l)\sin X",
    "11.90": r"\psi\propto X(t)\cos(\pi z)\sin(kx),\ \ T'\propto Y(t)\cos(\pi z)\cos(kx)+Z(t)\sin(2\pi z)",
    "11.91": r"\dot X=\Pr(Y-X),\ \ \dot Y=-XZ+rX-Y,\ \ \dot Z=XY-bZ",
    "11.92": r"\Big(\frac{d^2}{dR^2}-k^2\Big)^2\hat u_R=(1+\alpha x)\hat u_\varphi",
    "11.93": r"\Big(\frac{d^2}{dR^2}-k^2\Big)\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R",
    "11.94": r"\hat u_\varphi=\sum_{m=1}^\infty C_m\sin(m\pi x)",
    "11.95": r"(U-c)\Big(\frac{d^2\hat v}{dy^2}-k^2\hat v\Big)-\frac{d^2U}{dy^2}\hat v=0",
    "11.96": r"\frac{\partial}{\partial t}(U_i+u_i)+(U_j+u_j)\frac{\partial}{\partial x_j}(U_i+u_i)=-\frac1\rho\frac{\partial}{\partial x_i}(P+p)+\nu\frac{\partial^2}{\partial x_j\partial x_j}(U_i+u_i)",
    # earlier chapters (as their notebooks show them)
    "1.29": r"N^2=-\frac{g}{\rho}\Big(\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\Big)",
    "3.5": r"\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F",
    "4.10": r"\nabla\cdot\mathbf u=0",
    "4.75": r"\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const}",
    "4.86": r"\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u",
    "4.89": r"\frac{DT}{Dt}=\kappa\nabla^2T",
    "5.8": r"\frac{D\Gamma}{Dt}=0",
    "7.127": r"N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}",
    "7.18": r"\Big(\frac{\partial\phi}{\partial z}\Big)_{z=0}\cong\frac{\partial\eta}{\partial t}",
    "7.21": r"\Big(\frac{\partial\phi}{\partial t}\Big)_{z=0}\cong-g\eta",
    "7.95": r"\omega^2=gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}",
    "7.96": r"E=\tfrac12(\rho_2-\rho_1)ga^2",
    "8.9": r"u_\varphi(R)=AR+\frac BR",
    "8.10": r"u_\varphi=\frac{1}{R_2^2-R_1^2}\Big\{\big[\Omega_2R_2^2-\Omega_1R_1^2\big]R-\big[\Omega_2-\Omega_1\big]\frac{R_1^2R_2^2}{R}\Big\}",
    "9.36": r"f'''+\frac{n+1}2ff''-nf'^2+n=0",
    "9.51": r"\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}<0\quad(dp/dx<0)",
    "9.52": r"\Big(\frac{\partial^2u}{\partial y^2}\Big)_{wall}>0\quad(dp/dx>0)",
    "9.71": r"u=u_0\,\mathrm{sech}^2\Big(\frac\eta{\sqrt6}\Big)",
}
CHAPTER_PREFIX = r"11\."

def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ *({num})*"


def EE(*nums: str) -> str:
    """Several equations, each with its number, joined for prose."""
    return " · ".join(E(n) for n in nums)


_EQ_REF = re.compile(r"\((\d{1,2}\.\d+[a-d]?)\)")


def _plain_eq_follows(after: str) -> bool:
    """Is the equation already written in plain symbols right after its number? Then inserting its LaTeX again would only
    duplicate it."""
    a = after.lstrip(" ,:*")[:40]
    pos = a.find("=")
    if pos <= 0 or pos > 14:
        return False
    head, rhs = a[:pos], a[pos + 1:].lstrip()
    if re.search(r"[A-Za-z]{4,}|[.;()→–]", head) or rhs[:1].isdigit():
        return False
    stop = {"is", "of", "to", "in", "at", "as", "be", "by", "the", "and", "so", "it", "on", "for", "with"}
    return not any(w.lower() in stop for w in head.split())


def _maths_before(before: str) -> bool:
    """Does maths end right before the number (possibly across a line break or an italic star)?"""
    return before.rstrip(" *\n").endswith("$")


_PRINTED = re.compile(r"printed|book prints|The text|the text (?:calls|says)|The book says|book writes", re.I)


def _names_printed(before: str) -> bool:
    """Does the text just before an equation number talk about the book's PRINTED (slipped) form? Then the corrected
    equation must not be inserted after the number; the printed form is written out by hand where it matters."""
    return bool(_PRINTED.search(before[-32:]))


def show_eqs(text: str) -> str:
    """Write the equation next to the first mention of a book equation number "(10.25)" in a text that names it without
    writing it — the house rule "show the equation, not just its number"."""
    done: set[str] = set()

    def rep(m):
        n = m.group(1)
        nxt = text[m.end():m.end() + 16]
        if (n in done or n not in EQ or EQ[n] in text or nxt.lstrip(", :").startswith("$")
                or _plain_eq_follows(text[m.end():]) or _maths_before(text[:m.start()])
                or _names_printed(text[:m.start()])
                or nxt[:1] in "/–-)" or text[max(0, m.start() - 1):m.start()] in ("(", "–", "-")):
            return m.group(0)
        done.add(n)
        return f"({n}), ${EQ[n]}$"
    return _EQ_REF.sub(rep, text)


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
        else:
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
    if "*in" in s and "words:*" in s:
        a, b = re.split(r"\*in\s+words:\*", s, maxsplit=1)
        return a.strip().rstrip("—").strip(), b.strip()
    return s.strip(), ""


def eq_display(nums: list[str], words: str = "") -> str:
    """Display-math body writing several book equations with their numbers (for Starts and Results that only name
    numbers in Part F)."""
    body = r",\qquad ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in nums)
    if len(nums) > 2:
        body = r"\begin{array}{l}" + r" \\ ".join(f"{EQ[n]}\\ \\text{{({n})}}" for n in nums) + r"\end{array}"
    return body + (r"\ \ \text{" + _tex_escape(words) + "}" if words else "")



def part_f() -> dict[str, dict]:
    """Parse Part F of the design into {D01: dict(title, stars, goal, start, plan, tools, assumptions, steps, result,
    check, meaning, traps, check_src)} — every field word for word; ``check_src`` is the fenced sympy cell when Part F
    gives one (the ★★★ derivations)."""
    part = DESIGN.split("## Part F", 1)[1]
    chunks = re.split(r"^### (D\d\d) · ", part, flags=re.M)
    out: dict[str, dict] = {}
    for key, body in zip(chunks[1::2], chunks[2::2]):
        code = re.search(r"```python\n(.*?)```", body, flags=re.S)
        check_src = textwrap.dedent(code.group(1)).strip("\n") if code else None
        body = re.sub(r"```python\n.*?```", "", body, flags=re.S)
        lines = body.splitlines()
        title = _abs_plain(re.split(r" — ★", lines[0])[0].strip())
        stars = lines[0].count("★")
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
        f = {k: _abs_plain(_join(v)) for k, v in fields.items() if "sympy" not in k.lower()}
        parsed = []
        for st in steps:
            s = _abs_plain(_join(st))
            bits = re.split(r"(?:^|\s·\s)\*(did|tex|why|plain|live|set|watch)(?:\s*\([^)]*\))?:\*\s*", s)
            d = {name: val.strip() for name, val in zip(bits[1::2], bits[2::2])}
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=d["why"], plain=d["plain"]))
        start_tex, start_plain = _split_words(f["Start"])
        res_tex, res_plain = _split_words(f["Result"])
        plan = [p.strip() for p in re.split(r"\(\d+\)\s*", f["Plan"]) if p.strip()]
        tools = [t.strip().rstrip(".") for t in f["Tools"].split(";") if t.strip()]
        notation = f.get("Steps", "").strip()
        assumptions = f.get("Assumptions", "") + (f" Notation in the steps: {notation.strip('()')}." if notation else "")
        out[key] = dict(title=title, stars=stars, goal=f["Goal"], start=(display_tex(start_tex), start_plain),
                        plan=plan, tools=tools, assumptions=assumptions, steps=parsed,
                        result=(display_tex(res_tex), res_plain), check=f.get("Check", ""),
                        meaning=f.get("What it means", ""), traps=f.get("Traps", ""), check_src=check_src)
    return out


PF = part_f()
assert len(PF) == 25 and sum(len(d["steps"]) for d in PF.values()) == 257, (len(PF), sum(len(d["steps"]) for d in PF.values()))

def pf_sub(key: str, field: str, old: str, new: str) -> None:
    """Edit one Part F field after parsing (pointers to the cells that run a check, verification decisions); fails
    loudly if the design text changed. ``field`` = goal/check/meaning/traps/assumptions/title, ``tools``, or
    ``stepN.did|tex|why|plain``."""
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


def pf_insert(key: str, before: int, did: str, tex: str, why: str, plain: str) -> None:
    """Insert a step that Part F skipped, so that every step follows from the one above by one stated move."""
    PF[key]["steps"].insert(before - 1, dict(did=did, tex=tex, why=why, plain=plain))


def renumber(key: str, start: int, by: int = 1) -> None:
    """After a step was inserted at position ``start``, bump every pointer 'step N' / 'steps N–M' with N ≥ start."""
    def bump(m: re.Match) -> str:
        a = int(m.group(2))
        a2 = a + by if a >= start else a
        b = m.group(3)
        tail = f"–{int(b) + by if int(b) >= start else int(b)}" if b else ""
        return f"{m.group(1)}{a2}{tail}"
    pat = re.compile(r"(steps? )(\d+)(?:–(\d+))?")
    d = PF[key]
    for fld in ("goal", "assumptions", "check", "meaning", "traps"):
        d[fld] = pat.sub(bump, d[fld])
    for st in d["steps"]:
        for part in ("did", "tex", "why", "plain"):
            st[part] = pat.sub(bump, st[part])


def pf_why_add(key: str, step: int, extra: str) -> None:
    """Complete a step whose Part F *why* is only a pointer: say why the move is allowed and why we make it (nbkit rule)."""
    PF[key]["steps"][step - 1]["why"] = PF[key]["steps"][step - 1]["why"].rstrip() + " " + extra


def eq_rows(n: str) -> list[str]:
    """One book equation as display rows: a long multi-part equation is split at its part separators; the number tags
    the last row."""
    rows = [r.strip() for r in EQ[n].split(r",\ \ ")]
    rows[-1] += r"\qquad\text{(" + n + ")}"
    return rows


def eq_display(nums: list[str], *extra: str) -> str:  # noqa: F811  (replaces the ch10 helper: long equations are stacked)
    """Display-math body writing several book equations with their numbers (for Starts and Results that only name
    numbers in Part F); ``extra`` rows (already LaTeX) are appended."""
    rows = [r for n in nums for r in eq_rows(n)] + list(extra)
    if len(rows) == 1:
        return rows[0]
    return r"\begin{array}{l}" + r" \\[4pt] ".join(rows) + r"\end{array}"


# ---------------------------------------------------------------------------------------------------------------------
# Starts and Results that Part F writes only as equation numbers (or as words): shown as maths here.
# ---------------------------------------------------------------------------------------------------------------------
def _set(key: str, field: str, tex: str) -> None:
    PF[key][field] = (tex, PF[key][field][1])


_set("D02", "result", eq_display(["11.9", "11.13"]))
_set("D03", "start", eq_display(["11.3", "11.4", "11.5", "11.9", "11.13", "11.14"]))
_set("D04", "start", eq_display(["11.15", "11.16", "11.17"],
                                r"\text{and now a wall under the lower layer: }\frac{\partial\phi_2}{\partial z}=0\ \text{at}\ z=-h"))
_set("D04", "result", r"\begin{array}{l}c=\frac{\rho_1U_1+\rho_2U_2\coth kh}{\rho_1+\rho_2\coth kh}\pm\Big[\frac{\frac gk(\rho_2-\rho_1)+\sigma_sk}"
                      r"{\rho_1+\rho_2\coth kh}-\frac{\rho_1\rho_2(U_1-U_2)^2\coth kh}{(\rho_1+\rho_2\coth kh)^2}\Big]^{1/2} \\[4pt] "
                      r"\Delta U_{\min}^2=\frac{2\sqrt{g\Delta\rho\,\sigma_s}(\rho_1+\rho_2)}{\rho_1\rho_2}\ \ \text{at}\ \ k_*=\sqrt{g\Delta\rho/\sigma_s}"
                      r"\quad(\text{deep water})\end{array}")
PF["D04"]["result"] = (PF["D04"]["result"][0], "air over water with surface tension 0.074 N/m: the minimum wind is 6.70 m/s at a wavelength of "
                       "1.73 cm — gravity guards the long waves, surface tension the short ones; between them a minimum shear exists.")
_set("D05", "result", eq_display(["11.25", "11.26", "11.27"]))
_set("D06", "result", eq_display(["11.29", "11.27"]))
_set("D07", "start", eq_display(["11.27", "11.29", "11.30"]))
_set("D07", "result", eq_display(["11.36", "11.37", "11.38"]))
_set("D08", "start", eq_display(["11.36", "11.37", "11.38"], r"\text{with a complex growth rate allowed: }\sigma=\sigma_r+\mathrm i\sigma_i"))
_set("D08", "result", r"\sigma_i=0\quad\text{for every Bénard mode when}\ \mathrm{Ra}>0")
_set("D09", "start", eq_display(["11.36", "11.37", "11.38"], r"\text{at the margin }\sigma=0\ \text{(D08)}"))
_set("D09", "result", eq_display(["11.40", "11.41"]))
_set("D10", "start", eq_display(["11.40", "11.41"]))
_set("D11", "start", eq_display(["11.40"], r"w=T'=\mu\Big(\frac{\partial u}{\partial z}+\frac{\partial w}{\partial x}\Big)=\mu\Big(\frac{\partial v}{\partial z}+"
                                           r"\frac{\partial w}{\partial y}\Big)=0\ \text{on the walls (p. 489)}"))
_set("D12", "start", eq_display(["11.36", "11.37"], r"\text{with stress-free isothermal walls: }W=W''=\hat T=0\ \text{on}\ z=\pm\tfrac12"))
_set("D13", "result", eq_display(["11.45", "11.46"]))
_set("D14", "start", eq_display(["11.47", "11.49"]))
_set("D14", "result", eq_display(["11.51", "11.52", "11.53"]))
_set("D15", "start", eq_display(["11.52", "11.51"]))
_set("D16", "result", eq_display(["11.55", "11.56", "11.57"]))
_set("D17", "start", eq_display(["11.55", "11.56", "11.57"], r"[\rho,p,\psi]=[\hat\rho,\hat p,\hat\psi](z)\,e^{\mathrm ik(x-ct)}\ \ \text{(one normal mode)}"))
_set("D17", "result", eq_display(["11.61", "11.62"]))
_set("D18", "start", eq_display(["11.61", "11.62"]))
_set("D19", "start", eq_display(["11.61", "11.62"], r"\text{with a growing mode, }c_i\neq0,\ \text{and stable stratification, }N^2\ge0"))
_set("D19", "result", r"\Big[c_r-\tfrac12(U_{\max}+U_{\min})\Big]^2+c_i^2\le\Big[\tfrac12(U_{\max}-U_{\min})\Big]^2\ \ \text{(p. 507)},\qquad "
                      r"kc_i\le\frac k2(U_{\max}-U_{\min})")
_set("D20", "result", eq_display(["11.77"]))
_set("D21", "result", eq_display(["11.79", "11.80"]))
_set("D22", "start", eq_display(["11.81", "11.82"]))
_set("D22", "result", r"c_i\neq0\ \Rightarrow\ \int\frac{\lvert\phi\rvert^2}{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}dy=0\ \Rightarrow\ "
                      r"\frac{d^2U}{dy^2}\ \text{changes sign somewhere in}\ (y_1,y_2)")
_set("D23", "result", eq_display(["11.88"], r"\text{2-D shear flow: }\frac d{dt}\int\tfrac12(u^2+v^2)\,dV=-\int uv\frac{\partial U}{\partial y}dV-\Lambda,\ \ "
                                            r"\Lambda=\nu\int\Big(\frac{\partial u_i}{\partial x_j}\Big)^2dV\ge0"))
_set("D24", "start", r"\begin{array}{l}u=-\frac{\partial\psi}{\partial z},\ \ w=\frac{\partial\psi}{\partial x}\quad\text{(2-D rolls between stress-free "
                     r"isothermal walls, }-\tfrac12\le z\le\tfrac12) \\[4pt] " + EQ["11.90"] + r"\qquad\text{(11.90)}\end{array}")
_set("D24", "result", eq_display(["11.91"], r"r=\frac{\mathrm{Ra}\,k^2}{(\pi^2+k^2)^3}=\frac{\mathrm{Ra}}{\mathrm{Ra}_{cr}},\qquad b=\frac{4\pi^2}{\pi^2+k^2}"))
_set("D25", "result", r"\begin{array}{l}(0,0,0)\ \text{unstable for}\ r>1;\qquad C_\pm=\big(\pm\sqrt{b(r-1)},\ \pm\sqrt{b(r-1)},\ r-1\big) \\[4pt] "
                      r"C_\pm\ \text{unstable for}\ r>r_H=\frac{\Pr(\Pr+b+3)}{\Pr-b-1}=24.74\quad(\Pr=10,\ b=8/3)\end{array}")

# ---------------------------------------------------------------------------------------------------------------------
# Steps whose Part F *why* is only a pointer: say why the move is allowed AND why we make it (nbkit rule).
# ---------------------------------------------------------------------------------------------------------------------
pf_why_add("D03", 11, "A quadratic a c² + b c + d = 0 always has the two roots (−b ± √(b² − 4ad))/(2a); the factor 2 of b cancels the 2 below.")
pf_why_add("D06", 7, "Nothing new is done here: we only write the outcome of the subtraction with the horizontal Laplacian named.")
pf_why_add("D07", 4, "The conditions themselves do not change; only the place where they hold is written in the new unit of length.")
pf_why_add("D07", 9, "The walls are the same for every horizontal wave, so the conditions pass unchanged to the amplitudes; W is a multiple of the velocity amplitude.")
pf_why_add("D08", 4, "Sorting by σ separates the part that carries the unknown growth rate from the part that does not.")
pf_why_add("D09", 5, "Two of the three conditions of each wall already involve only W, so they can be kept as they are.")
pf_why_add("D12", 3, "The sine is never zero inside the layer, so it can be divided out, leaving a relation between the two amplitudes.")
pf_why_add("D12", 4, "The same division by the common sine leaves a second relation between the two amplitudes.")
pf_why_add("D14", 7, "Continuity is the only equation without a time derivative, so it gives the axial velocity directly; dividing by ik is allowed for k ≠ 0.")
pf_why_add("D15", 1, "We want the simplest case first — the one of the worked example below.")
pf_why_add("D16", 2, "A derivative of U along x is zero, so only the vertical velocity can change the momentum the fluid carries.")
pf_why_add("D16", 4, "The base flow has no vertical velocity, so no shear term appears here; gravity acts on the density disturbance.")
pf_why_add("D17", 1, "Differentiating an exponential multiplies it by the coefficient in its exponent, so each derivative becomes a number.")
pf_why_add("D17", 8, "Inserting the density removes the last unknown except the stream-function amplitude.")
pf_why_add("D18", 8, "A total derivative is what integration by parts needs in the next step.")
pf_why_add("D19", 2, "We need these two derivatives to rewrite every term of the Taylor–Goldstein equation in the new unknown.")
pf_why_add("D20", 6, "Each coefficient of the equations is independent of x, z and t, so the exponential passes through unchanged.")
pf_why_add("D21", 8, "Dividing a whole equation by the same nonzero number keeps it true, and it puts the equation in the book's form.")
pf_why_add("D22", 2, "Dividing isolates the highest derivative, so that multiplying by the conjugate gives a real, positive integral.")
pf_why_add("D22", 8, "A product is zero only if one factor is; the first factor is not, so the integral must vanish.")
pf_why_add("D23", 2, "Multiplying a momentum equation by the velocity turns it into an energy equation — the same move as for Bernoulli in Ch. 4.")
pf_why_add("D23", 6, "The first piece is a divergence (it will integrate to zero); the second is a sum of squares with a minus sign — a pure loss.")
pf_why_add("D23", 7, "Every term has now been sorted into three boxes: a source, a sink, and a divergence that only moves energy around.")
pf_why_add("D23", 11, "The surface integral has vanished, and the remaining viscous integral is a sum of squares, hence never negative.")
pf_why_add("D25", 1, "A product is zero only if a factor is zero, and the Prandtl number is a positive material constant.")
pf_why_add("D25", 4, "Either the first factor vanishes (no motion) or the bracket does; the bracket has real roots only for r > 1.")

# D12: the Result (the two growth rates) is never reached by a step in Part F — add it as the last step, so the numbering
# of the steps shared with the explainer's Derivation tab does not change.
PF["D12"]["steps"].append(dict(
    did="Solve the quadratic for σ",
    tex=r"\sigma_\pm=\frac{\Pr}2\Big[-a^2\big(1+\tfrac1\Pr\big)\pm\sqrt{a^4\big(1-\tfrac1\Pr\big)^2+\frac{4\mathrm{Ra}K^2}{\Pr a^2}}\Big]",
    why="The quadratic formula (Ch. 6 P159) applied to step 6, whose coefficients are 1/Pr, a²(1 + 1/Pr) and a⁴ − RaK²/a²; "
        "under the root is exactly the positive discriminant of step 7. We want the growth rate itself, not only its zero.",
    plain="two real growth rates: the plus root is the one that can become positive."))
PF["D12"]["goal"] += " (The last step — the quadratic formula written out — is our addition, so that the two growth rates are reached by a step.)"

# Cross-chapter numbers that this project's earlier notebooks number differently: cite the section instead.
pf_sub("D16", "step6.why", "Dρ̃/Dt (Ch. 3 (3.5))", r"Dρ̃/Dt, the material derivative $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ (3.5) of Ch. 3,")
pf_sub("D16", "tools", "the material derivative (Ch. 3 (3.5))", r"the material derivative $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ (3.5) of Ch. 3")
pf_sub("D02", "check", "Ch. 7's (7.18) and (7.21)",
       r"Ch. 7's linearised surface conditions $(\partial\phi/\partial z)_{z=0}\cong\partial\eta/\partial t$ (7.18) and $(\partial\phi/\partial t)_{z=0}\cong-g\eta$ (7.21)")
# the measured numbers of the tested functions (the design's scratch values differ in the last digits)
pf_sub("D13", "check", "Thermocline numbers (C06): lhs = 61 254 ✓.", "Thermocline numbers (C06's code cell prints them): the left side is about 6.1 × 10⁴, far above 657 ✓.")
pf_sub("D19", "check", "c = 0.4267i at k = 0.445", "c ≈ 0.426i at k = 0.445")
pf_sub("D07", "check", "every term of (11.36)–(11.37) is a pure number", "every term of the two amplitude equations is a pure number")

# pipeline words → reader words
pf_sub("D06", "traps", "— the planted mutant", "— the mistake the check cell below plants on purpose")
pf_sub("D06", "check", "Planted variant ∇² instead of ∇_H² fails", "A deliberately wrong variant (∇² instead of ∇_H²) fails")
pf_sub("D21", "check", "planted v̂ = +ikφ changes", "a deliberately wrong sign, v̂ = +ikφ, changes")

if "--dump" in sys.argv:
    for k, d in PF.items():
        print(f"=== {k} {d['title']}  ({len(d['steps'])} steps, {d['stars']} stars, check_src={bool(d['check_src'])})")
        print("  START", d["start"][0][:400], "|", d["start"][1][:80])
        for i, s in enumerate(d["steps"], 1):
            print(f"  {i}. {s['did']} :: {s['tex'][:90]}  || why={len(s['why'].split())}w")
        print("  RESULT", d["result"][0][:400], "|", d["result"][1][:80])
        print("  TOOLS", d["tools"])
    sys.exit(0)

# ---------------------------------------------------------------------------------------------------------------------
# Lesson review, round 1: back-pointers replaced by the actual line, two terms explained, one step made unambiguous.
# ---------------------------------------------------------------------------------------------------------------------
pf_sub("D13", "step1.why", "as D05 steps 4–8, with conduction and diffusion profiles T̄(z), S̄(z) of constant gradients.",
       r"insert $\tilde T=\bar T(z)+T'$ and $\tilde s=\bar S(z)+s'$ into the heat equation and the salt equation "
       r"$\frac{\partial\tilde s}{\partial t}+(\tilde{\mathbf u}\cdot\nabla)\tilde s=\kappa_s\nabla^2\tilde s$. The two profiles have constant "
       r"gradients, so their Laplacians vanish; $(\mathbf u\cdot\nabla)\bar T=w\,d\bar T/dz$ because $\bar T$ depends on $z$ only; products of "
       r"disturbances are dropped (the same moves as in D05).")
pf_sub("D13", "step3.why", "the D06 moves with the new buoyancy.",
       r"start from the vertical momentum equation $\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p}{\partial z}+g(\alpha T'-\beta s')+\nu\nabla^2w$; "
       r"take its Laplacian; the divergence of the momentum equation gives $\nabla^2p=\rho_0g\frac{\partial}{\partial z}(\alpha T'-\beta s')$; "
       r"differentiate that in $z$ and subtract, so the pressure cancels and $\nabla^2-\partial_z^2=\nabla_H^2$ remains (the moves of D06).")
pf_sub("D13", "step4.why", "as D07 step 2, now also making w dimensionless, w = (κ/d)W; heat diffusion sets the time unit.",
       r"put $t\to(d^2/\kappa)t$, lengths $\to d$ and $w=(\kappa/d)W$. Every term of the heat equation then carries $\kappa/d^2$ except "
       r"$w\,d\bar T/dz=(\kappa/d)W\,d\bar T/dz$; dividing by $\kappa/d^2$ leaves the factor $d$. The salt equation is divided by the same "
       r"$\kappa/d^2$ (heat diffusion sets the time unit), which leaves $\kappa_s/\kappa$ on its Laplacian.")
pf_sub("D13", "step6.why", "D07 steps 5–6 at the margin.",
       r"on a mode $e^{\mathrm i(kx+ly)+\sigma t}$ replace $\partial_t\to\sigma=0$, $\nabla_H^2\to-K^2$, $\nabla^2\to D^2-K^2$ in the two "
       r"equations of step 4 and in step 3 scaled the same way (divide it by $\nu\kappa/d^5$).")
pf_sub("D04", "step7.why", "D03 steps 10–12 with ρ₂ replaced by ρ₂ coth kh (same identity).",
       r"write $\rho_2'\equiv\rho_2\coth kh$ and expand step 6 in powers of $c$: "
       r"$(\rho_1+\rho_2')c^2-2(\rho_1U_1+\rho_2'U_2)c+\rho_1U_1^2+\rho_2'U_2^2-\frac gk(\rho_2-\rho_1)-\sigma_sk=0$. The quadratic formula and the "
       r"identity $(\rho_1U_1+\rho_2'U_2)^2-(\rho_1+\rho_2')(\rho_1U_1^2+\rho_2'U_2^2)=-\rho_1\rho_2'(U_1-U_2)^2$ (the one of D03 step 12) give "
       r"the line above.")
PF["D14"]["steps"][13]["did"] = "Rescale the swirl amplitude"
PF["D14"]["steps"][13]["tex"] = r"\hat u_\varphi^{\,\text{new}}=\frac{2k^2d^2\Omega_1}{\nu}\,\hat u_\varphi^{\,\text{old}}"
pf_why_add("D14", 14, "Read it as: the new swirl amplitude is the old one multiplied by the constant; from here on the plain symbol means the new one.")
pf_sub("D15", "check", "(2.4 % = O(d/R₁))", "(shortcut/exact = 1.025: the shortcut is 2.5 % higher, a difference of order d/R₁)")
pf_why_add("D25", 7, "This is called a *pitchfork* bifurcation: one fixed point (the origin) loses its stability and two symmetric new ones "
           "(C₊ and C₋) appear — on a plot of the steady X against r, one line splits into three prongs.")
pf_sub("D22", "check", "(`rayleigh_identity_check`)", "(`ch11.rayleigh_identity_check`, run in the cell below)")
pf_sub("D25", "check", "(`ch11.lorenz_eigs`, `lorenz_hopf_r`)", "(`ch11.lorenz_eigs`, `ch11.lorenz_hopf_r` — all three printed by the code cell below)")
pf_sub("D24", "check", "`ch11.lorenz_sympy()` residuals 0.", "the check cell below projects both equations and prints the residuals of the three rescaled equations: 0, 0, 0.")

# D24 (★★★): Part F's check covers steps 9–10 and b; this adds step 7 and the rescaling of steps 11–15
D24_EXTRA = r"""
vort = sp.diff(lap(psi), t)/Pr - Ra*sp.diff(th, x) - lap(lap(psi))       # step 5, vorticity residual (J(psi, lap psi) = 0, step 4)
eqA = sp.solve(proj(vort, sp.cos(sp.pi*z)*sp.sin(k*x)), sp.diff(A, t))[0]   # step 7: project on the roll's own shape
print(sp.simplify(eqA - Pr*(Ra*k/a2*B - a2*A)))                   # -> 0: dA/dt = Pr (Ra k B/a^2 - a^2 A)
r_ = Ra*k**2/a2**3                                                # step 14: r = Ra/Ra_c(k)
aL = sp.pi*k/(sp.sqrt(2)*a2); bL = aL*Ra*k/a2**2; gL = sp.pi*r_   # step 13: the three scale factors alpha_L, beta_L, gamma_L
Xs, Ys, Zs = sp.symbols("X Y Z")                                  # the rescaled amplitudes X = alpha_L A, Y = beta_L B, Z = gamma_L C
back = {A: Xs/aL, B: Ys/bL, C: Zs/gL}                             # A, B, C written with X, Y, Z
dX = (aL*eqA/a2).subs(back)                                       # steps 11-12: dX/dtau = alpha_L (dA/dt)/a^2, with tau = a^2 t
dY = (bL*eqB/a2).subs(back)                                       # ... dY/dtau
dZ = (gL*eqC/a2).subs(back)                                       # ... dZ/dtau
bpar = 4*sp.pi**2/a2                                              # step 15: b = 4 pi^2/(pi^2 + k^2)
print(sp.simplify(dX - Pr*(Ys - Xs)), sp.simplify(dY - (-Xs*Zs + r_*Xs - Ys)), sp.simplify(dZ - (Xs*Ys - bpar*Zs)))   # -> 0 0 0: (11.91)
print(ch11.lorenz_sympy()["scaled_residuals"])                    # the library's own engine does the same projection and rescaling: [0, 0, 0]
"""


if "(11.61)" in PF["D19"]["steps"][2]["did"]:               # same step title as the explainers' Derivation tabs; the equation stays in the step
    PF["D19"]["steps"][2]["did"] = "Insert into the Taylor–Goldstein equation"
    PF["D19"]["steps"][2]["why"] = "The equation is the Taylor–Goldstein equation (11.61). " + PF["D19"]["steps"][2]["why"][0].upper() + PF["D19"]["steps"][2]["why"][1:]
# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch10.py)
# ---------------------------------------------------------------------------------------------------------------------
def core(cid: str, title: str, question: str) -> None:
    """A CORE block heading with its id (the title carries the block's key equation in LaTeX)."""
    nb.core(cid, f"{title} `{cid}`", question=textwrap.dedent(question).strip())


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` is exactly the Part E concept text (the ledger check matches it)."""
    nb.primer(term, f"`{pid}` · {textwrap.dedent(text).strip()}", code=textwrap.dedent(code).strip("\n") if code else None)


def _ledger_rows() -> list[list[str]]:
    part = DESIGN.split("## Part E", 1)[1].split("## Part F", 1)[0]
    rows = [[x.strip() for x in ln.strip().strip("|").split("|")] for ln in part.splitlines() if ln.strip().startswith("|")]
    return [r for r in rows if len(r) >= 3 and not set("".join(r)) <= set("-: ") and r[0].lower() != "concept"]


LEDGER = _ledger_rows()
_PRIMERS_MD = (ROOT / "knowledge" / "primers.md").read_text(encoding="utf-8").splitlines()


def _gist(pid: str) -> str:
    """The one-line gist of an earlier primer from knowledge/primers.md (our words), tidied for a reminder line."""
    for ln in _PRIMERS_MD:
        if not ln.startswith("|"):
            continue
        cells = [c.strip() for c in re.split(r"(?<!\\)\|", ln.strip().strip("|"))]
        if len(cells) >= 3 and re.search(r"\(" + re.escape(pid) + r"[,)]", cells[0]):
            g = cells[2].replace("\\|", "|")
            g = re.sub(r"\s*\((?:Eqs?\.\s*)?\(?\d{1,2}\.\d+[a-z]?\)?(?:[–,-]\s*\(?\d{1,2}\.\d+[a-z]?\)?)*\)", "", g)   # no bare equation numbers
            g = re.sub(r"\bEqs?\.\s*\(?\d{1,2}\.\d+[a-z]?\)?", "", g)
            g = re.sub(r"\(\d{1,2}\.\d+[a-z]?\)", "", g)
            g = g.replace("$", "").replace("\\", "")
            g = re.sub(r"_\{([^{}]*)\}", r"_\1", g)
            cut = re.split(r"(?<=[^0-9])[;.] (?=[A-Za-z⚠(`])", g)                # first clause or two
            out = cut[0]
            for extra in cut[1:]:
                if len(out) + len(extra) > 230:
                    break
                out += "; " + extra
            return out.strip().rstrip(";.") + "."
    return ""


def remind(cid: str, lead: str = "") -> None:
    """One 🔁 cell reminding the tools primed in earlier chapters that Part E lists as first used in ``cid`` — one
    sentence each (the gist recorded in knowledge/primers.md), the concept text verbatim so the ledger check matches."""
    items = [(r[0], r[2]) for r in LEDGER if "primers.md" in r[2] and r[1] == cid]
    if not items:
        return
    lines = []
    for concept, by in items:
        refs = re.findall(r"\(ch(\d\d) (P\d+a?)\)", by)
        where = ", ".join(f"Ch. {int(ch)}, {p}" for ch, p in refs)
        if concept in REMIND:                        # our own sentence, written for the way the tool is used in this chapter
            lines.append(f"- **{concept}** — {REMIND[concept]} *({where})*")
        else:                                        # fall back to the gist recorded in knowledge/primers.md
            gists = [f"{_gist(p) or 'see the earlier primer.'}" for ch, p in refs]
            lines.append(f"- **{concept}** — " + " ".join(gists) + f" *({where})*")
        nb.primers.append(f"{concept} (reminder)")
    head = "> 🔁 **Tools from earlier chapters used in this block** (one line each; the full primers are in the chapters named)"
    nb.md(f"{head}{(' — ' + lead) if lead else ''}\n\n" + "\n".join(lines), tags=["primer"])


def note(nid: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note whose text opens with its curation id(s) in bold ("**N09 [B]**")."""
    nb.note(f"**{nid}** {textwrap.dedent(text).strip()}", equation=equation, ref=ref)


def _title_no_numbers(title: str) -> str:
    """Drop bare equation numbers from a derivation title (the heading shows the key equation itself instead)."""
    num = r"\(\d{1,2}\.\d+[a-d]?(?:,\s*(?:\d{1,2}\.\d+)?[a-d]?)*\)"
    t_ = re.sub(r"\s*(?:—\s*)?" + num + r"(?:\s*(?:→|–|,|and|or)\s*" + num + r")*", "", title)
    t_ = re.sub(r"\(from\s*\)", "", t_)
    t_ = re.sub(r"\s+([,:)])", r"\1", t_)
    t_ = re.sub(r"\s+from(?=\s*[,;)]|\s*$)", "", t_)
    t_ = re.sub(r"\s+(?:and|from)\s+(?=and\b|from\b)", " ", t_)
    t_ = re.sub(r"\s*[—:,]\s*$", "", t_)
    t_ = re.sub(r":\s*→\s*", ": ", t_)
    t_ = re.sub(r"(?:\s*→\s*)+", " → ", t_)
    t_ = re.sub(r"^\s*→\s*|\s*→\s*$", "", t_)
    t_ = re.sub(r"\b(?:from|and)\s*$", "", t_).strip()
    t_ = re.sub(r"\s+(?:from|and)\s+(?=[,;:]|$)", "", t_)
    return re.sub(r"\s{2,}", " ", t_).strip(" ,:")


STARS = {1: "★", 2: "★★", 3: "★★★"}


def D(key: str, ref: str = "", check_src: str | None = None, extra_check: str = "", after: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then its traps as a ⚠️ callout. ``check_src`` defaults
    to the sympy cell Part F gives (★★★); ``after`` = a slip box placed between the derivation and its traps."""
    d = PF[key]
    if ref in EQ:
        short = EQ[ref] if len(EQ[ref]) < 230 else EQ[ref].split(r",\ \ ")[0] + r",\ \dots"
        ref = f"{ref}: ${short}$"                     # the heading shows the key equation, not only its number
    def S(text: str) -> str:                          # working notes of the pipeline are not for the reader; equations are
        text = re.sub(r"\s*\((?:verifier|analyst) note [ivx]+\)", "", text)
        text = re.sub(r"\(optional `check_src`: ", "(the check cell below: ", text)
        text = re.sub(r"Optional `check_src`: ", "The check cell below: ", text)
        return show_eqs(text)
    src = check_src if check_src is not None else d["check_src"]
    check = S(d["check"]) + (f" {extra_check}" if extra_check else "")
    goal = S(d["goal"]) + (f"\n\n**Assumptions.** {S(d['assumptions'])}" if d["assumptions"] else "")
    steps = [dict(st, why=S(st["why"]), plain=S(st["plain"])) for st in d["steps"]]
    nb.derivation(key, f"{_title_no_numbers(d['title'])} `{key}` {STARS[d['stars']]}", ref=ref, goal=goal,
                  start=(d["start"][0], S(d["start"][1])), plan=[S(x) for x in d["plan"]], uses=[S(x) for x in d["tools"]],
                  steps=steps, result=(d["result"][0], S(d["result"][1])), interpret=S(d["meaning"]), check=check,
                  check_src=textwrap.dedent(src).strip("\n") if src else None)
    if key in NOTES_IN:
        nb.md(f"📝 **Notes taught inside `{key}`:** {NOTES_IN[key]}.")
    if after:
        nb.md(textwrap.dedent(after).strip())
    if d["traps"]:
        nb.md(f"> ⚠️ **Common confusion (traps in `{key}`):** {S(d['traps'])}")


def see_read_change(see: str, read: str, change: str) -> None:
    nb.figure_notes(see, read, change)


def whatif(text: str) -> None:
    nb.md(f"**What would change if…** {textwrap.dedent(text).strip()}")


def idea(sketch: str = "", words: str = "") -> None:
    """The idea block: an ASCII sketch plus one or two sentences."""
    sk = textwrap.dedent(sketch).strip("\n")
    body = ("```\n" + sk + "\n```") if sk else ""
    nb.md("#### The idea\n\n" + body + (("\n\n" if body else "") + textwrap.dedent(words).strip() if words else ""))


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + textwrap.dedent(text).strip())


def explainer(slug: str, heading: str, why_static: str, tries: list[str]) -> None:
    """An embedded explainer, opening with the "why interactive rather than a static figure" sentence (ch04 lesson)."""
    why = (f"**Why interactive rather than a static figure:** {textwrap.dedent(why_static).strip()}\n\n"
           "Start with the **Walkthrough**; the **Explain** tab works every number out with your settings, and the "
           "**Derivation** tab steps through the same derivation as above.")
    nb.explainer(slug, heading=heading, why=why, tries=tries)


def confusion(text: str) -> None:
    nb.md("> ⚠️ **Common confusion:** " + textwrap.dedent(text).strip())


def slip(k: int, printed: str, correct: str) -> None:
    """A printed-vs-correct box (the book's printed slip #k, taught in corrected form)."""
    nb.md(f"> ⚠️ **slip #{k} — the book prints** {textwrap.dedent(printed).strip()} **; the correct form is** {textwrap.dedent(correct).strip()}")


def scratch(src: str, explain: str) -> None:
    """A from-scratch cell (hand-written version next to the tested function, then an assert) with its explanation."""
    nb.code(src, explain=explain, tags=["from-scratch"])


def gloss(title: str, text: str) -> None:
    """A one-paragraph gloss of a small tool used right here (no primer needed: a sentence says it)."""
    nb.md(f"*{title}.* {textwrap.dedent(text).strip()}")


def code(src: str, explain: str | None = None) -> None:
    nb.code(src, explain=explain)


def fig(src: str, see: str, read: str, change: str, explain: str | None = None) -> None:
    nb.figure(src, see=see, read=read, change=change, explain=explain)


_PROSE_SPLIT = re.compile(r"(```.*?```|\$\$.*?\$\$|\$[^$]+\$|`[^`\n]*`)", re.S)
_GROUP = re.compile(r"\((?:Eqs?\.\s*)?(\d{1,2}\.\d+[a-d]?(?:\s*(?:,|–|-|and)\s*(?:\d{1,2}\.\d+[a-d]?|[a-d]))*)\)")


def _labels(group: str) -> list[str]:
    out, base = [], ""
    for tok in re.split(r"\s*(?:,|–|-|and)\s*", group):
        if re.fullmatch(r"\d{1,2}\.\d+[a-d]?", tok):
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
                adj = (_maths_before(unit[:s0]) or unit[e0:].lstrip(" *,:").startswith("$")
                       or _names_printed(unit[:s0]))
                plain = len(labs) == 1 and _plain_eq_follows(unit[e0:])
                out.append((s0, e0, labs, adj, plain))
        pos += len(seg)
    merged: list = []
    for s0, e0, labs, adj, plain in out:
        if merged and s0 == merged[-1][1] + 1 and unit[merged[-1][1]] in "–-" and len(merged[-1][2]) == 1 and len(labs) == 1:
            a_, b_ = merged[-1][2][0], labs[0]
            ch_a, n_a = a_.split(".")
            ch_b, n_b = b_.split(".")
            span = [a_, b_]
            if ch_a == ch_b and n_a.isdigit() and n_b.isdigit() and 0 < int(n_b) - int(n_a) <= 4:
                span = [f"{ch_a}.{j}" for j in range(int(n_a), int(n_b) + 1) if f"{ch_a}.{j}" in EQ]
            merged[-1] = (merged[-1][0], e0, span, merged[-1][3] or adj, False)
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
                e2 = e + 2 if src[e:e + 2] == "**" else (e + 1 if src[e:e + 1] == "*" else e)   # after a closing bold/italic mark
                new.append(src[last:e2] + ", $" + r",\ ".join(EQ[n] for n in miss) + "$")
                last = e2
            if new:
                units[u] = "".join(new) + src[last:]
        out = "\n---\n".join(units)
        if out != c.source:
            c.source = out
            changed += 1
    return changed


_UNI_TEX = {"−": "-", "₀": "_0", "₁": "_1", "₂": "_2", "θ": r"\theta ", "λ": r"\lambda ", "ω": r"\omega ", "±": r"\pm ",
            "ρ": r"\rho ", "σ": r"\sigma ", "η": r"\eta ", "ν": r"\nu ", "δ": r"\delta ", "√": r"\sqrt ", "ε": r"\varepsilon ",
            "π": r"\pi "}


def _exp_to_maths(m: re.Match) -> str:
    """'e^{−iθ}' written in plain text (Part F why/in-words lines) → the maths $e^{-i\\theta}$."""
    body = "".join(_UNI_TEX.get(ch, ch) for ch in m.group(0)).replace("²", "^2").replace("³", "^3")
    return "$" + re.sub(r"\s+\}", "}", body) + "$"


def tidy_raw_tex() -> int:
    """Plain-text exponents copied from Part F ('e^{−iθ}', 'x^{−1/4}') would show their braces on the page: an
    exponential becomes inline maths, any other braced power becomes '^(…)'."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        parts = _PROSE_SPLIT.split(c.source)
        for k in range(0, len(parts), 2):
            seg = re.sub(r"(?:(?<![A-Za-z\\])[A-Za-z0-9])?(?<!\\)(?:e\^\{[^{}$]*\})+", _exp_to_maths, parts[k])
            sub = re.split(r"(\$[^$]+\$)", seg)
            seg = "".join(q if q.startswith("$") else re.sub(r"\^\{([^{}$]*)\}", r"^(\1)", q) for q in sub)
            parts[k] = seg
            if k and parts[k].startswith("$") and parts[k - 1].endswith("$"):
                parts[k] = " " + parts[k]
            if k + 1 < len(parts) and parts[k].endswith("$") and parts[k + 1].startswith("$"):
                parts[k] += " "
        out = "".join(parts)
        assert out.count("$$") == c.source.count("$$"), c.source[:120]
        if out != c.source:
            c.source = out
            changed += 1
    return changed


def self_check_prose() -> list[str]:
    """No TeX command or TeX exponent outside maths or code, and no equation group named without being shown."""
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "markdown":
            continue
        prose = "".join(seg for k, seg in enumerate(_PROSE_SPLIT.split(c.source)) if k % 2 == 0)
        cmds = re.findall(r"\\[A-Za-z]+|[\^_]\{", prose) + re.findall(r"\\text\{[^}]*\\_", c.source)
        if cmds:
            bad.append(f"cell {i}: TeX outside maths {sorted(set(cmds))[:5]}: {c.source[:80]!r}")
        maths = [seg for k, seg in enumerate(_PROSE_SPLIT.split(c.source)) if k % 2 and seg.startswith("$")]
        garbled = [m_ for m_ in maths if "^(" in m_]
        if garbled:
            bad.append(f"cell {i}: garbled exponent in maths {garbled[:2]}")
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


def self_check_near(window: int = 160) -> list[str]:
    """Stricter than the coverage tool: EVERY equation number named in prose (not just the first per cell) must have
    maths within `window` characters of it — the equation next to its number, as the lesson reviewer reads it — and must
    be an equation the builder knows (a number missing from EQ would silently escape the other checks)."""
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "markdown":
            continue
        src, pos = c.source, 0
        eqs, p_ = [], 0                                            # (start, end) of maths segments that are equations
        for k, seg in enumerate(_PROSE_SPLIT.split(src)):
            if k % 2 and seg.startswith("$") and re.search(r"=|\\le|\\ge|<|>|\\approx|\\in\b|\\sim|\\to|\\equiv|\\propto|\\Rightarrow|\\cong", seg):
                eqs.append((p_, p_ + len(seg)))
            p_ += len(seg)
        for k, seg in enumerate(_PROSE_SPLIT.split(src)):
            if k % 2 == 0:
                for m in _GROUP.finditer(seg):
                    s0, e0 = pos + m.start(), pos + m.end()
                    if not re.match(r"11\.", m.group(1)):
                        continue                                   # other chapters' numbers: the recap rule handles them
                    toks = re.findall(r"\d{1,2}\.\d+", m.group(1))
                    unknown = [t for t in toks if t not in EQ]
                    near = any(b > s0 - window and a < e0 + window for a, b in eqs)   # an equation (not a lone symbol) nearby
                    near = near or bool(re.match(r"[\s*,:]*(?:as\*\*|\*\*)?\s*\$[^$]{12,}\$", src[e0:]))   # a printed form written right after
                    if unknown or not near:
                        bad.append(f"cell {i}: ({m.group(1)}){' unknown ' + str(unknown) if unknown else ''}"
                                   f"{'' if near else ' no maths nearby'}: …{src[max(0, s0 - 60):e0 + 30]!r}")
            pos += len(seg)
    return bad


def self_check_ctrl() -> list[str]:
    """No control character other than newline in any cell: a TeX command in a non-raw string loses its backslash
    ('\\to' → TAB + 'o', '\\frac' → form feed + 'rac', '\\beta' → backspace + 'eta')."""
    bad = []
    for i, c in enumerate(nb.cells):
        hits = sorted({hex(ord(ch)) for ch in c.source if ord(ch) < 32 and ch != "\n"})
        if hits:
            k = next(j for j, ch in enumerate(c.source) if ord(ch) < 32 and ch != "\n")
            bad.append(f"cell {i}: control characters {hits}: …{c.source[max(0, k - 40):k + 20]!r}")
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


def self_check_ledger() -> list[str]:
    """Every Part E row explained by a primer (new or reminded) has a matching primer entry — the coverage tool's rule
    5, run at build time with the same 18-character match."""
    primers = [p.lower() for p in nb.primers]
    miss = []
    for concept, _first, by, *_ in LEDGER:
        if "primer" in by.lower():
            key = re.sub(r"[`*$\\]", "", concept).lower().split("(")[0].strip()
            if key and not any(key[:18] in p or p[:18] in key for p in primers):
                miss.append(concept)
    return miss


# =====================================================================================================================

# one-sentence reminders of tools primed in earlier chapters (Part E rows "knowledge/primers.md: … — reminder"), written for
# the way each tool is used in THIS chapter; a concept missing here falls back to the gist stored in knowledge/primers.md
REMIND = {
    "complex numbers and Euler's formula": r"$\mathrm i^2=-1$ and $e^{\mathrm i\theta}=\cos\theta+\mathrm i\sin\theta$: a complex exponent is a wave, and a real part in the exponent is growth or decay.",
    "complex amplitudes and Re{·}": r"we write a wave as $\hat u\,e^{\mathrm i(kx-\omega t)}$ and mean its real part; the complex number $\hat u$ carries size and phase, and $\partial/\partial x\to\mathrm ik$.",
    "complex conjugate": r"$z^*=a-\mathrm ib$ flips the sign of the imaginary part; $zz^*=\lvert z\rvert^2\ge0$ — the trick behind every integral theorem of this chapter.",
    "the complex plane in numpy (1j, np.real, np.imag)": "`1j` is the imaginary unit; `np.real(z)`, `np.imag(z)`, `np.abs(z)`, `np.conj(z)` work on whole arrays.",
    "Fourier modes and amplification factors (Ch. 10 bridge)": r"any disturbance is a sum of waves $e^{\mathrm ikx}$, and a linear rule acts on each wave separately — in Ch. 10 it multiplied each mode by a factor per time step, here by $e^{\sigma\Delta t}$.",
    "linear stability by eigenvalues (perturb, linearise)": r"displace a steady state slightly, keep the terms linear in the displacement, and look for solutions $e^{\lambda t}$: any eigenvalue with positive real part means unstable.",
    "solve_ivp": "`scipy.integrate.solve_ivp(f, (t0, t1), y0)` integrates a system of ODEs with steps it chooses itself to meet a tolerance.",
    "animate and show_animation (frames player)": "`animate(update, frames)` builds a matplotlib animation from a function that moves the drawn objects; `show_animation(..., player=\"frames\")` gives step buttons, `\"video\"` plays smoothly.",
    "show_viz": "`show_viz(\"ch11\", slug)` embeds an interactive explainer in the notebook (its **Explain** tab works every number out with your settings).",
    "assert np.allclose": "stops with an error if two results differ beyond a tolerance — our proof that a from-scratch version matches the library.",
    "matplotlib figures": "`fig, ax = plt.subplots()`, `ax.plot(...)`, labels with units; every figure is followed by *What you see / How to read it / What would change if*.",
    "numpy arrays and broadcasting": "arithmetic between arrays (or an array and a number) acts on every element at once — no loops.",
    "f-strings": "`f\"{x:.3f}\"` inserts a value into text with three decimals.",
    "level sets and their normals, (11.8)": r"a surface $f=0$ has the normal $\nabla f/\lvert\nabla f\rvert$; for the interface $f=z-\zeta(x,t)$.",
    "A″ = k²A and exponential solutions": r"a linear ODE with constant coefficients is solved by trying $e^{mz}$; here $m^2=k^2$, so $e^{\pm kz}$.",
    "separation of variables": r"try (a function of $z$) × (a wave in $x$ and $t$): the PDE becomes an ODE for the function of $z$.",
    "Taylor transfer of a boundary condition to z = 0": r"a condition on the moving surface $z=\zeta$ is written on $z=0$ with $F(\zeta)\approx F(0)+\zeta F'(0)$; the correction is of second order.",
    "orders of smallness, dropping products": "a product of two small quantities is much smaller than either; linear theory keeps first-order terms only.",
    "quadratic formula with complex roots": r"$ac^2+bc+d=0$ has $c=\big(-b\pm\sqrt{b^2-4ad}\big)/(2a)$; a negative number under the root gives a complex-conjugate pair.",
    "inequalities under a sign change": "multiplying an inequality by a positive number keeps its direction, by a negative number reverses it.",
    "moving frame, c_r = mean": "seen from a frame moving at constant velocity, every velocity is reduced by the frame's velocity; the physics is unchanged.",
    "hyperbolic functions cosh, sinh, coth": r"$\cosh x=(e^x+e^{-x})/2$, $\sinh x=(e^x-e^{-x})/2$, $\coth x=\cosh x/\sinh x\to1$ for large $x$.",
    "Laplace pressure jump σ_s ∂²ζ/∂x²": r"a curved interface with surface tension $\sigma_s$ carries a pressure jump $\sigma_s\times$ curvature; for a gentle slope the curvature is $\partial^2\zeta/\partial x^2$.",
    "sympy": "algebra with symbols: `sp.symbols`, `sp.diff`, `sp.simplify`; an expression that simplifies to 0 is an identity.",
    "sympy expand, collect, subs": "`sp.expand` multiplies out, `sp.collect` groups by powers, `.subs` substitutes.",
    "scipy.integrate.quad (mixing energy cross-check)": "`quad(f, a, b)` returns the integral of a function and an error estimate.",
    "slider_figure": "a plotly figure whose curves follow a slider; every slider position is computed in advance, so it also works on the web page.",
    "np.logspace": "`np.logspace(a, b, n)` gives n numbers from 10^a to 10^b with a constant ratio between neighbours.",
    "power laws and log axes": "on a logarithmic axis equal distances are equal ratios — the right axis when a quantity spans several decades.",
    "Python dictionaries (returned results)": "`r[\"Ra_c\"]`: fluidpy returns several named results in a dict.",
    "functions as arguments and lambda": "`lambda y: 1 - y**2` is a one-line function; the eigen-solvers take the velocity profile as such a function.",
    "tuple unpacking": "`cp, cm = f()` takes a function's two return values apart in one line.",
    "sign of Γ: (11.21) vs Kundu ch01 vs meteorology (slip #10)": "to change convention, negate the number and flip the inequality.",
    "scaling w ~ κ/d, buoyancy/viscous ~ Ra": "replace each term by its typical size to see which terms balance and which dimensionless number measures their ratio.",
    "divergence of a vector equation; commuting constant-coefficient operators": r"for smooth fields derivatives can be taken in any order, so $\nabla^2$, $\partial_t$ and $\partial_z$ commute when the coefficients are constants.",
    "scaled variables and the chain rule": r"writing $z=d\,z^*$ turns $\partial/\partial z$ into $(1/d)\,\partial/\partial z^*$; each derivative brings one factor of the scale.",
    "operator substitution ∂_t → σ, ∇_H² → −K²": r"on a mode $e^{\mathrm i(kx+ly)+\sigma t}$ every derivative is a multiplication, so operators can be handled like numbers.",
    "eigenvalues and eigenvectors": r"$A\mathbf v=\lambda\mathbf v$ with $\mathbf v\neq0$: the special numbers $\lambda$ for which a nonzero solution exists.",
    "generalised symmetric eigenproblem (contrast)": r"$A\mathbf v=\lambda B\mathbf v$; with symmetric matrices (`eigh`) the eigenvalues are real — here the matrices are not symmetric and we need `eig`.",
    "collocation": "require the equation to hold exactly at a set of points; the unknowns are the values at those points.",
    "integration by parts with boundary terms": r"$\int f\,g'\,dz=[fg]-\int f'\,g\,dz$; the bracket vanishes when $f$ or $g$ is zero at both ends.",
    "matrices and determinants; nonzero solution ⇔ det = 0": r"a homogeneous system $M\mathbf a=0$ has a nonzero solution only if $\det M=0$.",
    "np.linalg.det": "`np.linalg.det(M)` returns the determinant of a square matrix (complex entries allowed).",
    "scipy.optimize.brentq": "`brentq(f, a, b)` finds a root of f between two points where f has opposite signs.",
    "scipy.optimize.minimize_scalar": "`minimize_scalar(f, bounds=(a, b), method=\"bounded\")` finds the minimum of a function of one variable.",
    "np.meshgrid": "`X, Y = np.meshgrid(x, y)` builds the grid of all (x, y) pairs for evaluating a field.",
    "contour and streamline plots": "`ax.contour` draws level lines of a field (streamlines when the field is a stream function), `ax.pcolormesh` colours it.",
    "live widgets": "an ipywidgets slider that re-runs a function while a Python kernel is running (frozen on the web page).",
    "integrals of sines over a period, orthogonality": r"$\int\sin m\pi z\,\sin n\pi z\,dz=0$ over the layer for $m\neq n$: different sine modes do not mix.",
    "Laplacian in cylindrical coordinates": r"$\nabla^2=\partial_R^2+\frac1R\partial_R+\frac1{R^2}\partial_\varphi^2+\partial_z^2$; the extra $1/R$ terms are curvature terms.",
    "anisotropic scaling with two length scales": "when one length (the gap) is much smaller than another (the radius), derivatives across the small length dominate.",
    "complex powers (U − c)^{1/2} and branches": r"a complex number has two square roots; any fixed choice works as long as the number never passes through zero.",
    "product rule": r"$(fg)'=f'g+fg'$.",
    "chain rule": r"$\frac{d}{dz}f(g(z))=f'(g)\,g'(z)$.",
    "Simpson and trapezoid (contrast)": "`np.trapezoid(y, x)` integrates sampled data with straight pieces — fine on a dense uniform grid, poor on a few clustered points.",
    "completing the square (tensors, recalled)": r"$x^2-2ax=(x-a)^2-a^2$.",
    "six Reynolds numbers (which length?)": "a Reynolds number means nothing without its length and velocity scales — half-width or diameter, centreline or mean speed.",
    "reading reference data from a file": "published benchmark numbers are read from `reference/ch11/` with their sources, never typed into the text.",
    "verification versus validation": "verification = are we solving the equations right (benchmarks, convergence); validation = are they the right equations (experiments).",
    "inflection point": r"where the curvature $U''$ changes sign.",
    "sech and (tanh)′ = sech²": r"$\mathrm{sech}\,x=1/\cosh x$ and $\frac{d}{dx}\tanh x=\mathrm{sech}^2x$.",
    "np.sign and np.nonzero": "`np.nonzero(np.diff(np.sign(f)))` gives the places where a sampled curve crosses zero.",
    "caching expensive runs": "results that take minutes are computed once and stored in a small file; the notebook reads them and recomputes one value to show the method.",
    "fundamental theorem of calculus": r"$\int_a^bF'(x)\,dx=F(b)-F(a)$.",
    "mean of a product of real parts": r"for two waves with complex amplitudes $\hat a,\hat b$ the average of their product over a wavelength is $\tfrac12\mathrm{Re}(\hat a\hat b^*)$ — zero when they are 90° out of phase.",
    "np.fft.rfft spectrum": "`np.fft.rfft(signal)` gives the amplitude at each frequency `np.fft.rfftfreq(n, dt)`.",
    "multivariable first-order Taylor expansion": r"$\mathbf f(\mathbf s^*+\boldsymbol\delta)\approx\mathbf f(\mathbf s^*)+J\boldsymbol\delta$ with $J$ the matrix of partial derivatives.",
    "RK4 by hand": "four slope evaluations per step, combined with weights 1, 2, 2, 1 — fourth-order accurate.",
    "solve_ivp options (DOP853, rtol, t_eval)": "`method=\"DOP853\"` is a high-order solver, `rtol` the relative tolerance, `t_eval` the output times.",
    "plotly 3-D (explainers, alternatives)": "interactive 3-D figures (used by the explainers); the static 3-D figure here uses matplotlib instead.",
}

# B/C notes of the curation that are taught INSIDE a derivation (design Part A): listed under the derivation with their ids,
# so that every note id has a visible home (no equation numbers here — the equations are in the steps named)
NOTES_IN = {
    "D02": "**N08 [B]** the exact kinematic condition, where the square-root factor cancels (step 5) · **N09 [B]** its linearised "
           "form on the flat plane (step 7) · **N10 [B]** the pressure matching (step 8) · **N11 [B]** the undisturbed balance "
           "(step 9) · **N12 [B]** the linearised dynamic condition (step 12)",
    "D03": "**N14 [B]** the remnants of the two conditions and the amplitudes they fix (steps 5–7) · **N15 [B]** the quadratic for "
           "the wave speed (steps 9–10)",
    "D05": "**N29 [B]** the exact (nonlinear) disturbance equations, with the heating term from vertical motion (step 7) · "
           "**N30 [B]** the linearised set (step 8)",
    "D06": "**N32 [B]** the Laplacian of the vertical momentum equation (step 1) · **N33 [B]** the pressure Poisson equation and "
           "its vertical derivative (steps 3–4) · **N34 [B]** the pressure-free equation for the vertical velocity (step 7)",
    "D07": "**N36 [B]** the non-dimensional equations and the Prandtl number (steps 2–4) · **N37 [B]** derivatives of a normal "
           "mode become numbers (step 6) · **N38 [B]** the amplitude equations before the velocity is rescaled (step 7)",
    "D09": "**N41 [B]** the marginal pair (steps 1–2) · **N42 [B]** the sixth-order equation (step 4) · **N43 [B]** the rigid "
           "conditions written with the velocity amplitude only (step 6)",
    "D10": "**N45 [B]** the characteristic equation, its three roots and the real exponent (steps 1–5) · **N46 [B]** the even "
           "solution and the 3 × 3 determinant, with slip #1 (steps 7–11)",
    "D11": "**N47 [B]** the stress-free conditions (steps 1–5) · **N48 [B]** the sine modes, with slip #2 (step 7) · **N49 [B]** "
           "the minimum over the wavenumber, with slip #3 (steps 9–11)",
    "D13": "**N55 [B]** the double-diffusive marginal equations and the second sign convention of the Rayleigh number (steps 5–8) "
           "· **N56 [B]** temperature and salt amplitudes locked together, and the threshold that follows (steps 9–12)",
    "D14": "**N60 [B]** the decomposition (step 1) · **N62 [B]** the linearised equations (steps 3–4) · **N63 [B]** normal modes "
           "and the narrow gap (steps 5, 11–12) · **N64 [B]** the narrow-gap pair (steps 13–15)",
    "D15": "**N65 [B]** the Taylor number and its inner-cylinder, narrow-gap form (steps 1–3)",
    "D16": "**N71 [B]** the linearised momentum equations, with slip #4 (steps 3–5)",
    "D17": "**N72 [B]** the three normal-mode equations (steps 2–4)",
    "D18": "**N75 [B]** Howard's change of variable and its two derivatives (steps 2–4) · **N76 [B]** the self-adjoint form — "
           "the book's 'after some rearrangement' written out (steps 5–8) · **N77 [B]** the integral identity and its "
           "imaginary part (steps 9–12)",
    "D19": "**N79 [B]** the second change of variable and its derivatives (steps 1–2) · **N80 [B]** the divergence form, with "
           "slip #11 (step 5) · **N81 [B]** the real and imaginary parts of the integral identity (steps 9–10) · **N82 [B]** "
           "the wave speed lies inside the velocity range (step 9)",
    "D20": "**N84 [B]** the linear streamwise momentum equation (step 3) · **N85 [B]** the other two momentum equations and "
           "continuity (step 4) · **N86 [B]** the three-dimensional normal mode (step 5) · **N87 [B]** the normal-mode "
           "equations (step 7)",
    "D25": "**N116 [B]** the fixed points, the pitchfork at r = 1 and the Hopf point (steps 4, 7, 11)",
}

# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=r"""
Every flow we solved so far was a *possible* flow — it satisfied the equations. Nature only shows us the possible flows that
*survive a nudge*. A pencil can balance on its tip in a textbook, never on a desk. This chapter asks of each steady flow the same
question: put a tiny wave-shaped disturbance on it — does the wave die or grow? Because the disturbance is tiny, the equations
become linear, a wave of one wavelength evolves on its own, and the question becomes an eigenvalue problem: a growth rate σ (or a
complex wave speed c) for every wavenumber. We meet the classic mechanisms one at a time — shear between two streams
(Kelvin–Helmholtz billows, the start of wind waves), heating from below (Bénard convection, the atmosphere's shallow convection),
two diffusivities (salt fingers in the subtropical ocean), rotation (Taylor vortices), shear against stratification (the
Richardson number ¼ behind every ocean and atmosphere mixing scheme), and plain viscous shear (Tollmien–Schlichting waves) —
learn the theorems that decide stability without solving anything (Rayleigh, Fjørtoft, Howard, Miles–Howard, Squire), compute
every number with one numerical tool (Chebyshev collocation), and end with three equations of convection that never settle down:
Lorenz's chaos, the reason weather has a predictability limit.""",
    roadmap=[
        "**Normal modes and the stability vocabulary** — one growth rate per wavelength decides everything (`C01`).",
        "**Kelvin–Helmholtz** — shear against gravity under a square root (`C02`).",
        "**The Bénard eigenproblem** — heating from below, two ODEs, a growth rate σ(K, Ra, Pr) (`C03`).",
        "**The rigid–rigid neutral curve** — convection starts at the bottom of a valley, Ra ≈ 1708 (`C04`).",
        "**Free walls** — the neutral curve by hand, 27π⁴/4 (`C05`).",
        "**Salt fingers** — diffusion that destabilises (`C06`).",
        "**Taylor vortices** between rotating cylinders — rotation plays gravity (`C07`).",
        "**The Taylor–Goldstein equation** — one ODE for stratified shear flow (`C08`).",
        "**The Richardson number ¼** — the Miles–Howard guarantee (`C09`).",
        "**Howard's semicircle** — where unstable wave speeds can be (`C10`).",
        "**Squire's theorem and the Orr–Sommerfeld equation** — viscous parallel flows (`C11`).",
        "**Rayleigh's inflection-point theorem** (and Fjørtoft's) — stability from the shape of a profile (`C12`).",
        "**Plane Poiseuille flow** — viscosity destabilises above Re = 5772 (`C13`).",
        "**The disturbance-energy budget** — Reynolds-stress production against dissipation (`C14`).",
        "**The Lorenz system** — deterministic chaos (`C15`).",
    ],
    prerequisites=[
        "potential flow and Laplace's equation (Ch. 6); interface waves (Ch. 7 §7.7)",
        "vortex sheets and their roll-up (Ch. 5 §5.8); unsteady Bernoulli (Ch. 4 §4.9)",
        "the Boussinesq equations (Ch. 4 §4.9); the buoyancy frequency and stratification (Ch. 1 §1.10, Ch. 7 §7.8)",
        "circular Couette flow (Ch. 8 §8.2); Blasius, Falkner–Skan and jets (Ch. 9)",
        "perturb–linearise–eigenvalues (Ch. 9); Fourier modes and amplification factors (Ch. 10 §10.2); the generalised eigenproblem (Ch. 10)",
    ])
nb.explainer_index([
    ("normal_mode_growth", "When is a flow 'unstable'?",
     "a disturbance is a sum of waves that grow or decay like e^{σt}; the flow is unstable if the growth rate is positive for any wavenumber, and onset is where the growth curve first touches zero"),
    ("kelvin_helmholtz_boundary", "Why does wind raise waves only above a threshold?",
     "gravity holds long waves, surface tension short ones, shear pushes all; two real wave speeds collide and become a growing/decaying pair"),
    ("benard_neutral_curve", "Why does convection start near Ra ≈ 1708?",
     "each cell width has its own marginal Ra; onset is the bottom of that valley"),
    ("salt_fingers", "How can a stably stratified column overturn?",
     "heat diffuses 100× faster than salt, so a displaced parcel keeps its salt and keeps sinking"),
    ("taylor_couette_onset", "Why do stacked vortices appear between spinning cylinders?",
     "rings swap and release energy when the squared circulation falls outward; viscosity adds the same 1708-type threshold"),
    ("richardson_shear_instability", "Why does Ri = ¼ decide whether a shear layer billows?",
     "Ri > ¼ everywhere guarantees stability; below ¼ instability becomes possible, not certain"),
    ("inviscid_shear_criteria", "Which profiles can be unstable without viscosity?",
     "an inflection point with a vorticity maximum, and every unstable wave speed inside Howard's semicircle"),
    ("orr_sommerfeld_neutral_curve", "How can viscosity make a channel flow unstable?",
     "the viscous layer shifts v against u so the Reynolds stress draws energy from the shear: a thumb-shaped neutral curve above Re = 5772"),
    ("lorenz_attractor", "How can three exact equations be unpredictable?",
     "past r ≈ 24.74 every steady state is unstable; nearby starts separate exponentially on a strange attractor"),
])
nb.setup()
nb.code(r"""
import json, time, logging                               # read small reference files; a wall-clock timer; Python's message system
import numpy as np                                       # arrays and maths (Ch. 1 primer P03)
import sympy as sp                                       # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                          # static figures (Ch. 1 primer P01)
from scipy import linalg                                 # dense linear algebra: linalg.eig(A, B) solves A v = c B v (primer P259 below)
from scipy.optimize import brentq, minimize_scalar       # root of a function (P108) and minimum of a function (P170)
from scipy.integrate import solve_ivp, quad              # ODE integrator (P31) and 1-D integral (P87)
from fluidpy import ch11_instability as ch11             # the tested chapter-11 module (re-exports the eigen-solver toolkit)
from fluidpy.core import stability as ST                 # the new toolkit itself: Chebyshev matrices, eigen-solvers, criteria
from fluidpy.core import waves as WV, laminar as LAM, boundary_layer as BL, jets as JET   # Ch. 7, 8, 9 functions reused here
from fluidpy.core import stratification as STRAT, kinematics as KIN                       # Ch. 1 N^2 and Ch. 3 frames of reference
from fluidpy import ch05_vorticity_dynamics as ch05, ch07_gravity_waves as ch07           # Ch. 5 vortex-sheet roll-up, Ch. 7 waves
from fluidpy import ch10_computational_fluid_dynamics as ch10                             # Ch. 10: dominant_frequency of a signal
from fluidpy.core.interact import slider_figure, live    # plotly sliders (P17) and live widgets (P47)
from fluidpy.core.anim import animate                    # matplotlib animations (P16); show_animation came with setup
from fluidpy.core.style import COLORS, savefig           # the house palette and a helper that saves PNGs to outputs/ch11
from fluidpy.core.thermo import G0                       # standard gravity 9.80665 m/s^2 (the default g of every function)
sys.path.insert(0, str(ROOT / "scripts"))                # the chapter's drawing helpers live in scripts/
from ch11_drawings import kh_sketch, benard_sketch, taylor_sketch, cv_sketch, semicircle_sketch, six_profiles   # drawing only
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless "font not found" notes
C_GROW, C_DECAY, C_NEUT, C_BASE = COLORS["rose"], COLORS["teal"], COLORS["accent"], COLORS["muted"]   # growing, decaying, neutral, base state
C_BUOY, C_SHEAR, C_SALT = COLORS["blue"], COLORS["orange"], COLORS["amber"]                           # buoyancy, shear / kinetic energy, salt
BENCH = json.loads((ROOT / "reference" / "ch11" / "benchmarks.json").read_text(encoding="utf-8"))     # published benchmarks, with sources
print(len([n for n in dir(ch11) if not n.startswith("_")]), "public names in fluidpy.ch11_instability")
print("ch11.orr_sommerfeld_eigs is ST.orr_sommerfeld_eigs:", ch11.orr_sommerfeld_eigs is ST.orr_sommerfeld_eigs)   # one function, two names
""", explain=r"""
1. Numerical, symbolic and plotting libraries (primed in Ch. 1). `linalg.eig`, `brentq`, `minimize_scalar` and `solve_ivp` are
   reminded where they are first used.
2. `ch11` is the chapter module. It **re-exports the new eigen-solver toolkit** `core/stability.py` (`ST`: Chebyshev matrices,
   the Rayleigh, Taylor–Goldstein and Orr–Sommerfeld solvers, the stability criteria), so `ch11.orr_sommerfeld_eigs` and
   `ST.orr_sommerfeld_eigs` are the same function (the last line prints `True`). Every function cites its § and equation and is
   tested in `tests/test_ch11.py`.
3. Functions from earlier chapters are reused, not rewritten: interface waves (`WV`), circular Couette flow (`LAM`), Blasius
   (`BL`), the plane jet (`JET`), the buoyancy frequency (`STRAT`), moving frames (`KIN`), the vortex-sheet roll-up (`ch05`).
4. `scripts/ch11_drawings.py` only draws (two streams, a heated layer, a gap with vortices, a control volume) — no physics lives there.
5. `C_GROW, C_DECAY, …` fix one colour per meaning for the whole chapter: growing rose, decaying teal, neutral/marginal purple,
   base state grey, buoyancy blue, shear orange, salt amber. `BENCH` holds the published benchmark numbers (Chandrasekhar 1961,
   Orszag 1971, …) read from `reference/ch11/benchmarks.json`; they are never typed into the text by hand.""")
nb.md(r"""
## ⚠️ Conventions in this chapter (read once, come back when a symbol surprises you)

**One letter, several meanings.** The book reuses symbols between sections. This is what each means here and what we write.

| Book symbol | Meanings in the book | What we write |
|---|---|---|
| $K$, $k$, $\kappa$ | Bénard horizontal wavenumber (§11.4); streamwise or axial wavenumber (§11.3, §11.6–11.10); thermal diffusivity | $K$ for Bénard, $k$ elsewhere; $\kappa$ heat, $\kappa_s$ salt |
| $\alpha$ | thermal expansion (§11.4–11.5); $\Omega_2/\Omega_1-1$ (§11.6) | $\alpha$ thermal; $\mu\equiv\Omega_2/\Omega_1$ for the rotation ratio |
| $\sigma$ | growth rate; surface tension (Exercises 11.1, 11.2) | $\sigma$ growth rate, $\sigma_s$ surface tension |
| $\Gamma$ | $-d\bar T/dz$ (§11.4); circulation $2\pi rU_\theta$ (§11.6); Ch. 1's lapse rate $dT/dz$ | see the Γ table below; code takes `dT` = T_bottom − T_top |
| $U_1$ | upper-stream speed (§11.3); $U$ at the inflection point (§11.9) | $U_1,U_2$ streams; $U_I$ at the inflection point |
| $\phi$ | velocity potential (§11.3); $\hat\psi/(U-c)^{1/2}$ (§11.7); stream-function amplitude (§11.8–11.9) | each block says which |
| $\psi$ | three sign conventions (see slip #9) | each block states its own |
| $c$, $\sigma$ | complex wave speed; growth rate | linked by $\sigma=-\mathrm i\lvert\mathbf K\rvert c$ (D01) |
| Ra | $>0$ heated from below (§11.4); $<0$ heated from below (§11.5) | Ra with the §11.4 sign; C06 says when the §11.5 sign is used |

**Thirteen printed slips, taught in corrected form.** Each has a short ⚠️ box where it is used.

| slip | the book prints | correct | where we fix it |
|---|---|---|---|
| #1 | third term of $(d^2/dz^2-K^2)^2W$ with coefficient $B$ | coefficient $C$: $C(q^{*2}-K^2)^2\cosh q^*z$ | C04, D10 step 9 |
| #2 | $W=A\sin(n\pi z)$ | $W=A\sin n\pi(z+\tfrac12)$ | C05, D11 step 7 |
| #3 | $\frac{d\mathrm{Ra}}{dK^2}=\frac{3(\pi^2+K^2)^2}{K^2}-\frac{3(\pi^2+K^2)^3}{K^4}$ | no 3 on the second term | C05, D11 step 9 |
| #4 | $-\frac1{\rho_0}\frac{\partial p}{\partial x}$ in the $w$-equation of the stratified set | $-\frac1{\rho_0}\frac{\partial p}{\partial z}$ | C08, D16 step 5 |
| #5 | "see (7.96)", $E=\tfrac12(\rho_2-\rho_1)ga^2$, for the static interface | Ch. 7's dispersion relation $\omega^2=gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}$ (7.95) | recap R09 |
| #6 | Tollmien's middle branch $1-b[1-(y/\delta)^2]$ | $1-b(1-y/\delta)^2$ | C13 (continued), N105 |
| #7 | $(d^2/dR^2-k^2)^2\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R$ in Exercise 11.9 | first power of the operator | C07, N124 |
| #8 | $\int[U_{\min}-U][U_{\max}-U]\,dz\le0$ without a weight | the weight $Q\ge0$ inside the integral from the start | C10, D19 step 12 |
| #9 | three stream-function signs: $u=\partial\psi/\partial z$ (§11.7), $u=\partial\psi/\partial y$ (§11.8), $u=-\partial\psi/\partial z$ (§11.14) | not an error — each block states its own | C08, C11, C15 |
| #10 | $\Gamma=-d\bar T/dz$ in §11.4, then Ra redefined with $+d\bar T/dz$ in §11.5 | both conventions side by side | C03, C06 |
| #11 | $-k^2(U-c)F$ in Howard's divergence form | $-k^2(U-c)^2F$ | C10, D19 step 5 |
| #12 | continuity $\frac{\partial}{\partial R}(R\tilde u_R)+\frac{\partial\tilde u_z}{\partial z}=0$ | $\frac1R\frac{\partial}{\partial R}(R\tilde u_R)+\frac{\partial\tilde u_z}{\partial z}=0$ | recap R15, D14 step 4 |
| #13 | the cross-reference "(7.128)" beside $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ in §11.7 | Ch. 7 defines the buoyancy frequency as its equation (7.127); (7.128) there is the momentum equation that follows | recap R18 |

**Four ways this chapter makes things dimensionless.**

| Sections | Length | Time | Notes |
|---|---|---|---|
| §11.4 (Bénard) | depth $d$ | $d^2/\kappa$ | $z\in[-\tfrac12,\tfrac12]$; $w$ stays dimensional until $W\equiv(\Gamma d^2/\kappa)\hat w$ |
| §11.6 (Taylor) | gap $d$ | $d^2/\nu$ | $x=(R-R_1)/d\in[0,1]$; "$d/dR$" in the narrow-gap equations means $d/dx$ |
| §11.8–11.11 (viscous shear flows) | $L$ per flow | $L/U_0$ | Poiseuille: half-width and centreline speed; Blasius: $\delta^*$ and $U_\infty$ |
| §11.14 (Lorenz) | depth $d$ | $d^2/\kappa$, then × $(\pi^2+k^2)$ | rescaled amplitudes $X,Y,Z$ |

**Which Γ? Three conventions for one temperature gradient** (a water layer $d=5$ mm, bottom 2 K warmer than the top):

| Convention | Definition | Value | "Heated from below" means |
|---|---|---|---|
| this chapter, with $\mathrm{Ra}=g\alpha\Gamma d^4/\kappa\nu$ (11.21) | $\Gamma\equiv-d\bar T/dz$ | $+400$ K/m | $\Gamma>0$ |
| Ch. 1, Kundu's lapse rate | $\Gamma\equiv dT/dz$ | $-400$ K/m | $\Gamma<0$ |
| meteorology | $\Gamma_{met}\equiv-dT/dz$ | $+400$ K/m | $\Gamma_{met}>0$ |

To convert: negate the number and flip the inequality (Ch. 1, P48). The code never takes a Γ at all — it takes `dT` =
T_bottom − T_top, so the sign cannot be mistyped.""")
nb.md(r"""
> 🔁 **Tools from earlier chapters used in this one.** Each block reminds the ones it uses in one line where they first appear.
> In short: complex numbers, Euler's formula and complex amplitudes (Ch. 1, 2, 6, 7) carry every wave; the exponential trial for
> a linear ODE (Ch. 1) and separation of variables (Ch. 7) give the depth structure; orders of smallness and the Taylor transfer of
> a boundary condition (Ch. 2, 7) linearise; integration by parts (Ch. 9) proves the theorems; eigenvalues (Ch. 2), collocation
> (Ch. 6) and the generalised eigenproblem (Ch. 10) compute the numbers; `brentq`, `minimize_scalar` and `solve_ivp` (Ch. 1, 3, 7)
> search and integrate; `slider_figure`, `animate` and `show_viz` (Ch. 1) draw.""")

# =====================================================================================================================
# A.1 §11.1 Introduction — N01 N02 N03
# =====================================================================================================================
nb.section("11.1", "Introduction", intro=r"""
**What is this section about?** A steady flow that satisfies the equations of motion is not necessarily a flow you will ever see.
If a tiny disturbance — a ripple, a puff of warmer air, the vibration of a pump — grows instead of dying, the steady flow is
*unstable* and something else takes its place: billows, convection cells, vortices, turbulence. This section sets up the
question; §11.2 gives the method used for the rest of the chapter.""")
note("N01 [B]", r"""
**Basic state and disturbance.** The *basic* (background) state is the steady flow we already know — two streams, a heated layer
at rest, Couette flow, Poiseuille flow. We add a *disturbance* so small that products of disturbances can be dropped (the same
move as Ch. 1's parcel argument for the buoyancy frequency and Ch. 9's perturb–linearise–eigenvalues). *Linear stability* asks
whether infinitesimal disturbances grow; *finite-amplitude* (nonlinear) instability is when only large enough kicks grow. Both can
happen to the same flow (the dimple below).""")
note("N02 [B]", r"""
**Four mechanical pictures** (animated below): a ball in a bowl (every kick dies: stable), on an upturned bowl (every kick grows:
unstable), on a flat table (it stays wherever you put it: neutral), and in a small dimple on a hilltop (small kicks die, a large
kick throws it out: stable to small, unstable to large disturbances).""")
nb.md(r"""
**Animation A1 — four balls, two kicks each.** `ch11.potential_well_demo(shape, x0)` integrates a damped ball in a potential
$V(x)$, $\ddot x=-V'(x)-\gamma\dot x$, with `solve_ivp` (Ch. 1, P31: an adaptive ODE integrator). Our four potentials: bowl
$\tfrac12x^2$, cap $-\tfrac12x^2$, plane $0$, dimple $\tfrac12x^2-\tfrac14x^4$ (its rim is at $\lvert x\rvert=1$). Use the ◀ ▶
buttons to step.""")
nb.animation(r"""
shapes = ("bowl", "cap", "plane", "dimple")                       # the four potentials of our Fig. 11.1 analogue
nf = 16 if FAST else 22                                           # number of frames
runs = {(s, x0): ch11.potential_well_demo(s, x0, damping=0.3, t_end=12.0, n=nf)   # damped ball, displaced by x0, released at rest
        for s in shapes for x0 in (0.3, 1.2)}                     # a small kick (0.3) and a large kick (1.2) for every shape
print({s: (bool(runs[(s, 0.3)]["escaped"]), bool(runs[(s, 1.2)]["escaped"])) for s in shapes})   # did the ball leave? (small, large)
fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.4))
xg = np.linspace(-2.0, 2.0, 200)                                  # positions at which the potential is drawn [-]
dots = {}                                                         # the moving balls, one per (shape, kick)
for ax, s in zip(axs.ravel(), shapes):
    V = runs[(s, 0.3)]["V"]                                       # the potential V(x) of this shape (a function)
    ax.plot(xg, V(xg), color=C_BASE)                              # the landscape the ball rolls on
    for x0, mk in ((0.3, "o"), (1.2, "s")):                       # circle = small kick, square = large kick
        col = C_GROW if runs[(s, x0)]["escaped"] else C_DECAY     # rose if the ball escapes, teal if it stays
        (dots[(s, x0)],) = ax.plot([x0], [V(x0)], mk, color=col, ms=9)
    ax.set(title=s, xlim=(-2, 2), ylim=(-1.2, 1.2), xticks=[], yticks=[])
fig.suptitle("small kick (circle) and large kick (square): teal returns or stays, rose escapes", fontsize=10)

def update(i):                                                    # frame i: move every ball to its position at output time i
    for (s, x0), d in dots.items():
        r = runs[(s, x0)]                                         # this ball's trajectory
        j = min(i, len(r["x"]) - 1)                               # an escaped ball stops being integrated: hold its last point
        xb = float(np.clip(r["x"][j], -2.0, 2.0))                 # keep the marker inside the panel
        d.set_data([xb], [max(float(r["V"](xb)), -1.15)])         # the ball sits on the curve V(x); one that has left the panel is held at its edge
    return list(dots.values())

show_animation(animate(update, frames=nf, fig=fig, interval=150), player="frames", dpi=60)   # a frame player: step through with the buttons
""", explain=r"""
1. `potential_well_demo` returns the times `t`, positions `x`, velocities `v`, the potential `V` (a function) and `escaped`
   (True once the ball has run away; the integration then stops, so its arrays are shorter).
2. The printed dictionary is the verdict for (small kick, large kick): the bowl keeps both, the cap loses both, the plane keeps
   both where they were put, the dimple keeps the small kick and loses the large one.
3. `update(i)` moves each marker along its potential curve (a ball that has run out of the panel — both balls on the cap — is
   held at the panel's edge, so you still see which way it went); `player="frames"` gives buttons instead of a video, because the point
   is to compare the four panels at the same instant.""")
nb.figure_notes(
    see="Four landscapes, each with two balls. In the bowl both roll back and settle; on the cap both run away; on the plane both "
        "stay put; in the dimple the circle settles and the square rolls over the rim and away.",
    read="Stability is a statement about what happens *after* a disturbance. The dimple shows that the answer can depend on the "
         "size of the disturbance: linearly stable (small kicks die), nonlinearly unstable (a large kick escapes).",
    change="…we remove the damping (`damping=0`): the bowl's ball oscillates forever — it neither grows nor decays. That is the "
           "difference between *neutral* and *stable* that §11.2 makes precise.")
note("N03 [C]", r"""
**What this chapter is for, and what it leaves out.** Linear theory predicts *when* a flow first becomes unstable and *which
wavelength* appears; it cannot say what the flow becomes afterwards (§11.12–11.14). We follow disturbances growing in *time*
everywhere at once (temporal instability); disturbances that grow as they travel downstream (spatial instability, Huerre &
Monkewitz 1990) are only named. No rotation of the frame here: the Coriolis force and the baroclinic instability of the
atmosphere and ocean come in Ch. 13 §13.17.""")

# =====================================================================================================================
# A.2 §11.2 Method of Normal Modes — C01
# =====================================================================================================================
nb.section("11.2", "Method of Normal Modes", intro=r"""
**What is this section about?** One recipe used in every section of the chapter: write the disturbance as a sum of waves, follow
each wave alone, and read off its growth rate. This section gives the recipe and the words — stable, neutral, unstable, marginal,
stationary, oscillatory — that the rest of the chapter uses.""")
core("C01", r"The normal mode $u=\hat u(z)\exp\{\mathrm ikx+\mathrm imy+\sigma t\}=\hat u(z)\exp\{\mathrm i\lvert\mathbf K\rvert(\mathbf e_K\cdot\mathbf x-ct)\}$ (11.1) and the stability vocabulary",
     "How can one number per wavelength decide whether a whole flow is stable?")
problem(r"""
A light breeze blows over a pond. The surface is never perfectly flat: there are ripples of every length, from millimetres to
metres, all at once. Some grow into waves, others die. To predict which, we would like to follow each ripple length separately —
and for small ripples we can, because small disturbances obey *linear* equations, where a sum of solutions is a solution. A
weather forecaster asks the same of a jet stream: which wavelength of meander grows fastest?""")
idea(r"""
any small disturbance  =  Σ over k   û(z) · e^{ikx} · e^{σ(k) t}          (Fourier, Ch. 10)
                                      ─────   ──────   ─────────
                                      shape   wave in x  grows or decays in time
linear equations ⇒ each k evolves alone ⇒ one number σ(k) = σ_r + iσ_i per wavelength
  σ_r < 0 for every k   → stable        (every ripple dies)
  σ_r > 0 for some k    → unstable      (that wavelength grows like e^{σ_r t})
  σ_i ≠ 0               → the pattern travels/oscillates;  σ_i = 0 → it grows in place (cells)
""", r"""
Stability is not a property of one wave but of the whole curve $\sigma_r(k)$: the flow is stable only if the curve stays below
zero everywhere, and as a control parameter rises, the wavelength where the curve first touches zero is the pattern you see.
Symbols: $k,m$ wavenumbers in $x,y$ [1/m]; $\sigma=\sigma_r+\mathrm i\sigma_i$ growth rate [1/s]; $c=c_r+\mathrm ic_i$ complex
wave speed [m/s]; $\hat u(z)$ the complex amplitude, a function of the remaining coordinate.""")
remind("C01")
P("P255", "necessary vs sufficient conditions, \"for every k\", and proof by contradiction", r"""
Three pieces of logic the chapter's theorems use. *Sufficient*: "if A then B" — A guarantees B (a Richardson number above ¼
everywhere guarantees stability). *Necessary*: "B only if A" — without A no B, but A alone does not force B (an inflection point
is needed for inviscid instability, but is not enough). "Stable" means $\sigma_r\le0$ **for every** $k$; "unstable" needs **one**
$k$ with $\sigma_r>0$. *Proof by contradiction*: assume the opposite (a growing mode exists), derive something impossible (a
positive number equal to a negative one), so the assumption was false. In numpy, `np.all(cond)` is True only if every entry of
`cond` is True; `np.any(cond)` is True if at least one is.""", code=r"""
import numpy as np                                # arrays
sig_r = np.array([-0.3, -0.1, 0.02, -0.5])        # growth rates at four wavenumbers [1/s]
print(np.all(sig_r < 0))                          # False: 'stable' needs EVERY k decaying
print(np.any(sig_r > 0))                          # True: one growing k is enough for 'unstable'
""")
note("N05 [B]", r"""
**Why one wave at a time is allowed.** The coefficients of the linearised equations do not depend on $x$, $y$ or $t$ (the basic
state is uniform in those directions), so $e^{\mathrm i(kx+my)+\sigma t}$ passes through every derivative unchanged:
$\partial/\partial x\to\mathrm ik$, $\partial/\partial t\to\sigma$ — the same trick as Ch. 10's Fourier error modes, where one
mode was multiplied by an amplification factor each step and here by $e^{\sigma\Delta t}$. Linearity lets a sum of modes evolve as
the sum of their evolutions, so any disturbance is covered. What is left is an ODE in $z$ with $\sigma$ as an unknown — an
**eigenvalue problem**.""")
code(r"""
x = np.linspace(0, 2*np.pi, 400)                           # positions along the interface [m]
ks = np.array([1.0, 2.0, 3.0])                             # three wavenumbers [1/m]
sigmas = np.array([-0.2, 0.1, -0.5])                       # their growth rates [1/s]: only k = 2 grows
t_end = 10.0                                               # look after 10 s
parts = [ST.normal_mode(x, 0.0, t_end, 1.0, k, sigma=s) for k, s in zip(ks, sigmas)]   # each mode alone, Re{e^{ikx + sigma t}} (11.1)
total = sum(parts)                                         # linear problem: the sum of the modes is the disturbance
for k, s, p in zip(ks, sigmas, parts):
    print(f"k = {k:.0f}: amplitude at t = 10 s is {np.max(np.abs(p)):.4f}   (e^(sigma t) = {np.exp(s*t_end):.4f})")
print(f"sum of all three: max = {np.max(np.abs(total)):.4f} - almost all of it is the k = 2 mode")
""", explain=r"""
1. `ST.normal_mode(x, y, t, uhat, k, sigma=…)` evaluates the real part of $\hat u\,e^{\mathrm ikx+\sigma t}$, the first form of (11.1).
2. Each mode's amplitude after 10 s is $e^{\sigma t}$: $e^{-2}=0.135$, $e^{1}=2.718$, $e^{-5}=0.0067$.
3. The sum is dominated by the one growing mode — after a while only the growing wavelength is visible, which is why "one $k$ with
   $\sigma_r>0$" is enough to call a flow unstable.""")
D("D01", ref="11.1")
note("N04 [B]", r"""
**Stable, neutral, unstable** (one mode): $\sigma_r<0$ (or $c_i<0$) is stable; $\sigma_r=0$ (or $c_i=0$) is neutrally stable;
$\sigma_r>0$ (or $c_i>0$) is unstable. For the *flow*: stable if every $k$ is stable, unstable if any $k$ is (P255).""")
code(r"""
print([ST.stability_class(sigma=s) for s in (-0.2, 0.0, 0.3 + 1j)])        # one mode each: decays, stays, grows (and travels)
print(ST.stability_verdict(np.array([-0.3, -0.1, 0.02, -0.5])))            # the flow: a whole array of growth rates over k
print(ST.marginal_type(0.0), ST.marginal_type(0.0 + 2.5j))                 # sigma = 0 exactly, or sigma = 2.5 i: how the onset looks
""", explain=r"""
1. `stability_class` reads the sign of $\sigma_r$ of one mode: `'stable'`, `'neutral'`, `'unstable'`.
2. `stability_verdict` applies "for every $k$" to an array: one positive entry (index 2) makes the flow `'unstable'`.
3. `marginal_type` looks at $\sigma_i$ where $\sigma_r=0$: `'stationary'` (the pattern grows in place) or `'oscillatory'`.""")
note("N06 [B]", r"""
**Marginal state; stationary or oscillatory onset.** The *marginal* state sits on the border between stable and unstable:
$\sigma_r=c_i=0$ **and** a small change of the control parameter (Ra, Re, Ta) makes $\sigma_r>0$. (A *neutral* mode has
$\sigma_r=0$ but need not sit on such a border — a stable interface wave is neutral for every parameter value.)

| onset | $\sigma$ at the margin | what you see | examples |
|---|---|---|---|
| **stationary** | $\sigma_r=0$ and $\sigma_i=0$ | a pattern of cells that grows in place | Bénard cells (C03–C05), Taylor vortices (C07) |
| **oscillatory** (older name: overstability) | $\sigma_r=0$, $\sigma_i\neq0$ | oscillations of growing amplitude | the diffusive regime of double diffusion (C06) |""")
nb.worked_example("both forms of the normal mode for one wave", r"""
Take $\mathbf K=(k,m)=(2,0)$ m⁻¹ and $c=1+0.5\mathrm i$ m/s.

1. Second form: $e^{\mathrm iK(x-ct)}=e^{\mathrm i2x}e^{-\mathrm i2(1+0.5\mathrm i)t}=e^{\mathrm i2x}e^{(1-2\mathrm i)t}$.
2. Compare with $e^{\mathrm ikx+\sigma t}$: $\sigma=1-2\mathrm i$ s⁻¹, i.e. $\sigma=-\mathrm iKc$ (D01).
3. Growth rate $\sigma_r=1$ s⁻¹ $=Kc_i=2\times0.5$ ✓: the amplitude multiplies by $e\approx2.72$ every second.
4. $\sigma_i=-2=-Kc_r$: the crests move at $c_r=1$ m/s toward $+x$ (a crest at $kx-2t=0$ moves right).
5. Verdict: unstable, and travelling (oscillatory at a fixed point).""")
code(r"""
sig = ST.sigma_from_c(2.0, 1 + 0.5j)                       # D01: sigma = -i K c for K = 2 1/m, c = 1 + 0.5i m/s
print("sigma =", sig, "| back to c:", ST.c_from_sigma(2.0, sig))            # (1-2j) and (1+0.5j)
print(ST.stability_class(sigma=sig), "| travelling" if abs(sig.imag) > 0 else "| growing in place")   # 'unstable'; sigma_i != 0 -> it travels
x = np.linspace(0, 2*np.pi, 200)                           # one stretch of the x axis [m]
u_a = ST.normal_mode(x, 0.0, 1.0, 1.0, 2.0, sigma=sig)     # first form of (11.1) at t = 1 s: Re{e^{ikx + sigma t}}
u_b = ST.normal_mode(x, 0.0, 1.0, 1.0, 2.0, c=1 + 0.5j)    # second form of (11.1): Re{e^{iK(x - ct)}}
print("largest difference between the two forms:", np.max(np.abs(u_a - u_b)))   # ~1e-16: the same wave
""", explain=r"""
1. `sigma_from_c` and `c_from_sigma` are the two directions of D01's result $\sigma=-\mathrm i\lvert\mathbf K\rvert c$.
2. The verdict of N04 for this mode, and whether it travels ($\sigma_i\neq0$).
3. The two forms of (11.1) evaluated on a grid agree to round-off.""")
scratch(r"""
u_mine = np.real(np.exp(1j*(2.0*x) + sig*1.0))                     # first form typed out: Re{exp(i k x + sigma t)} at t = 1 s
u_mine2 = np.real(np.exp(1j*2.0*(x - (1 + 0.5j)*1.0)))             # second form typed out: Re{exp(i K (x - c t))}
assert np.allclose(u_mine, u_a) and np.allclose(u_mine2, u_b)      # same numbers -> the library does exactly this
assert np.allclose(u_mine, WV.real_field(np.exp(sig.real*1.0), 2.0*x + sig.imag*1.0))   # Ch. 7's Re{a e^{i phase}} with a = e^{sigma_r t}
print("from-scratch normal mode = ST.normal_mode = Ch. 7's real_field")
""", r"""
1. The two forms of (11.1) written with `np.exp` and `np.real` only.
2. `assert np.allclose(...)` (Ch. 1, P15) stops with an error if they differ from the tested function.
3. The same wave through Ch. 7's `WV.real_field(amp, phase)`: amplitude $e^{\sigma_rt}$, phase $kx+\sigma_it$.""")
fig(r"""
k = np.linspace(20.0, 600.0, 300)                                  # wavenumbers [1/m] (wavelengths from 31 cm down to 1 cm)
grow = np.array([ch11.normal_mode_growth("interface", ki, rho1=1000.0, rho2=1.2, surface_tension=0.074).real
                 for ki in k])                                     # sigma_r(k) [1/s]: water (1000) ABOVE air (1.2), with surface tension
grow0 = np.array([ch11.normal_mode_growth("interface", ki, rho1=1000.0, rho2=1.2, surface_tension=0.0).real for ki in k])   # no surface tension
k_cut = 2*np.pi/ch11.rayleigh_taylor_cutoff(0.074, 1000.0, 1.2)    # the wavenumber where growth stops [1/m]
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.5))
a.plot(k, grow0, "--", color=C_BASE, label="no surface tension")
a.plot(k, grow, color=C_GROW, label=r"with $\sigma_s$ = 0.074 N/m")
a.axvline(k_cut, color=C_NEUT, ls=":")
a.text(k_cut + 8, 6, f"cut-off\n$k$ = {k_cut:.0f} m$^{{-1}}$", color=C_NEUT, fontsize=8)
a.set(xlabel="wavenumber $k$ [m$^{-1}$]", ylabel=r"growth rate $\sigma_r$ [s$^{-1}$]")
a.set_title("(a) water over air: long waves grow", fontsize=10)
a.legend(fontsize=8)
kk = np.array([60.0, 210.0, 340.0])                               # three modes [1/m]
ss = [ch11.normal_mode_growth("interface", ki, rho1=1000.0, rho2=1.2, surface_tension=0.074) for ki in kk]   # their complex sigma
xs = np.linspace(0.0, 0.25, 300)                                   # 25 cm of interface [m]
ts = np.linspace(0.0, 0.30, 120)                                   # the first 0.3 s [s]
eta = np.array([sum(ST.normal_mode(xs, 0.0, t, 1.0, ki, sigma=si) for ki, si in zip(kk, ss)) for t in ts])   # the three modes summed
eta = eta/np.max(np.abs(eta), axis=1, keepdims=True)               # each time row scaled by its own maximum (so the shape is visible)
b.pcolormesh(xs*100, ts, eta, cmap="RdBu_r", shading="auto")
b.set(xlabel="$x$ [cm]", ylabel="time $t$ [s]")
b.set_title("(b) three modes summed: the fastest wins", fontsize=10)
savefig(fig, "ch11", "c01_growth_curve")
print(f"fastest of the three: k = {kk[int(np.argmax([s.real for s in ss]))]:.0f} 1/m with sigma_r = {max(s.real for s in ss):.1f} 1/s")
""",
    see="(a) A hump of growth rates: positive for long waves, falling to zero at the cut-off (purple dotted line) and zero beyond "
        "it; without surface tension (grey dashed) the curve keeps rising. (b) A space–time picture of three waves added "
        "together: at the bottom (t = 0) a mixed pattern, at the top a single clean wavelength.",
    read="The top of the hump is the wavelength you will see; the zero crossing is the shortest wave that can grow (about "
         "1.7 cm here — N121 in C02 derives it). In (b) the pattern at late times is the mode with the largest "
         "$\\sigma_r$: growth is a property of the whole curve, and its maximum wins.",
    change="…we remove surface tension: the hump never comes back down — every shorter wave grows faster (the 'short waves "
           "always lose' of C02).",
    explain=r"""
1. `ch11.normal_mode_growth("interface", k, rho1, rho2, surface_tension)` returns the complex growth rate $\sigma=-\mathrm ikc_+$
   of a static interface (upper density `rho1`, lower `rho2`). With the heavy fluid on top it is real and positive for long waves
   (Rayleigh–Taylor instability); it takes one wavenumber at a time, hence the list comprehension.
2. Panel (b) adds three modes with `ST.normal_mode` and rescales each time row, so you watch the *shape* change, not the size.""")
explainer("normal_mode_growth", "When is a flow 'unstable'?", r"""
Dragging $k$ along the growth-rate curve makes the interface grow, decay or travel at once, and raising the control parameter
lifts the whole curve through zero — the wavelength that goes first appears in front of you; a static plot hides the link between
one point on the curve and the motion it means.""", [
    "Pick 'upside-down interface' and drag k from long to short waves: watch the growth rate peak, then die where surface tension wins.",
    "Switch to 'Bénard (free walls)' and raise Ra slowly from 500: note the K where the growth rate first crosses zero (K ≈ 2.22 at Ra ≈ 657.5).",
    "Open the complex σ-plane: at the onset the root crosses the imaginary axis — on the real axis for Bénard (stationary), off it for a travelling interface wave.",
])
whatif(r"""
…the basic state moved — two streams sliding past each other? Then $c$ is no longer purely real or purely imaginary: it has a
travelling part and a growing part at once, and the growth comes from the shear itself. That is the Kelvin–Helmholtz problem
(C02), the first complete normal-mode calculation.""")
# =====================================================================================================================
# A.3 §11.3 Kelvin–Helmholtz Instability — R01–R11, C02, S01
# =====================================================================================================================
nb.section("11.3", "Kelvin-Helmholtz Instability", intro=r"""
**What is this section about?** Two layers of fluid slide past each other: wind over the sea, a fast current over a slow one in
the ocean, the warm air of a front over cold air. Gravity wants the interface flat (if the heavy fluid is below); the velocity
jump wants to roll it up. We put a small wave on the interface, apply the conditions you met in Ch. 7, and get a formula for the
wave speed whose square root turns imaginary when shear wins. The result explains billow clouds, wind-driven ripples and the start
of mixing in the thermocline.""")
nb.recap("R01", "Total potentials and Kelvin's theorem", r"""
Each layer is irrotational: far upstream the flow is uniform, and Kelvin's theorem (Ch. 5 §5.2, $\frac{D\Gamma}{Dt}=0$ (5.8): the
circulation of a material loop is constant in inviscid, barotropic flow with conservative forces) keeps it irrotational, so each
layer has a potential, $\tilde\phi_1=U_1x+\phi_1,\ \tilde\phi_2=U_2x+\phi_2$ (11.2), with $\phi_j$ the small disturbance
(subscript 1 = upper layer, 2 = lower layer).""", where="Ch. 5 §5.2, Ch. 6")
nb.recap("R02", "Laplace's equation for the disturbances", r"""
Incompressible + irrotational ⇒ $\nabla^2\phi_1=0,\ \nabla^2\phi_2=0$ (11.3) — exactly Ch. 6's potential flow and Ch. 7 §7.1.""",
         where="Ch. 6, Ch. 7 §7.1")
nb.recap("R03", "Decay far away", r"""
$\phi_1\to0$ as $z\to+\infty$ (11.4) and $\phi_2\to0$ as $z\to-\infty$ (11.5): the disturbance lives near the interface, as for
Ch. 7's interface waves; this decides the signs $e^{\mp kz}$ in D03.""", where="Ch. 7 §7.7")
nb.recap("R04", "Kinematic and dynamic interface conditions", r"""
Fluid on either side moves with the interface, $\mathbf n\cdot\nabla\tilde\phi_1=\mathbf n\cdot\mathbf U_s=\mathbf
n\cdot\nabla\tilde\phi_2$ on $z=\zeta$ (11.6) ($\mathbf n$ the unit normal, $\mathbf U_s$ the velocity of the interface), and —
without surface tension — the pressure is continuous, $p_1=p_2$ on $z=\zeta$ (11.7). These are the two-fluid versions of Ch. 7's
free-surface conditions; they are the starting line of D02.""", where="Ch. 4 §4.10, Ch. 7 §7.1")
nb.recap("R05", "The normal of a level set", r"""
The interface is the level set $f=z-\zeta(x,t)=0$, so its unit normal is $\mathbf n=\nabla f/\lvert\nabla f\rvert=
[-(\partial\zeta/\partial x)\mathbf e_x+\mathbf e_z]/\sqrt{1+(\partial\zeta/\partial x)^2}$, which turns the kinematic condition
into $\mathbf n\cdot\{\partial_x\tilde\phi_1\mathbf e_x+\partial_z\tilde\phi_1\mathbf e_z\}=\mathbf n\cdot\{\partial_t\zeta\,\mathbf
e_z\}=\mathbf n\cdot\{\partial_x\tilde\phi_2\mathbf e_x+\partial_z\tilde\phi_2\mathbf e_z\}$ (11.8) — the gradient of a level-set
function is normal to the level set (Ch. 2, P75).""", where="Ch. 2 P75, Ch. 7 §7.1")
nb.recap("R06", "Unsteady Bernoulli in each layer", r"""
Irrotational, unsteady: $\frac{\partial\tilde\phi_j}{\partial t}+\tfrac12\lvert\nabla\tilde\phi_j\rvert^2+\frac{p_j}{\rho_j}+gz=C_j$
(11.10), $j=1,2$ — Ch. 4's unsteady Bernoulli equation with constant density in each layer; $C_j$ is a constant, one per layer.""",
         where="Ch. 4 §4.9")
nb.recap("R07", "A″ = k²A and its exponentials", r"""
A constant-coefficient second-order ODE is solved by the trial $e^{mz}$ (Ch. 1, P44): $A''=k^2A$ gives $m^2=k^2$, so
$A=A_\pm e^{\pm kz}$ (the depth structure after separating variables, Ch. 7 P167).""", where="Ch. 1 P44, Ch. 7 P167")
nb.recap("R08", "The decaying solutions", r"""
With the decay conditions: $\phi_1=A_-\exp\{\mathrm ik(x-ct)-kz\}$, $\phi_2=A_+\exp\{\mathrm ik(x-ct)+kz\}$ (11.15) — the same
potentials as Ch. 7's interface waves (`ch07.interface_fields`).""", where="Ch. 7 §7.7")
nb.recap("R09", "The static interface (U₁ = U₂ = 0)", r"""
Then $c=\pm\big[\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk\big]^{1/2}$ (11.19): the interface waves of Ch. 7,
$\omega^2=gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}$ (7.95) with $\omega=kc$; neutral if the heavy fluid is below ($\rho_2>\rho_1$),
**Rayleigh–Taylor** unstable if it is above (the square root is then imaginary). A parity cell in C02 checks `ch11.kh_phase_speed`
against Ch. 7's `WV.interface_omega`.""", where="Ch. 7 §7.7")
slip(5, r"""a pointer "see (7.96)", which is the energy $E=\tfrac12(\rho_2-\rho_1)ga^2$, for the static-interface limit""",
     r"""the interface dispersion relation of Ch. 7, $\omega^2=gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}$ (7.95).""")
nb.recap("R10", "The vortex sheet", r"""
Equal densities $\rho_1=\rho_2$: the interface is a vortex sheet of strength $\gamma=U_2-U_1$ (Ch. 5 §5.8). **New here:** the
dispersion relation becomes $c=\frac{U_2+U_1}2\pm\mathrm i\frac{U_2-U_1}2$ (11.20), so the growth rate $kc_i=k\lvert U_2-U_1\rvert/2$
is positive at **every** $k$ — a vortex sheet is unstable to all wavelengths, fastest for the shortest.""", where="Ch. 5 §5.8")
nb.recap("R11", "Roll-up of a perturbed vortex sheet", r"""
Ch. 5's `sheet_rollup` (a row of point vortices with a smoothing parameter $\delta$) follows the sheet past the linear stage into
cat's-eye rolls. Animation A2 below puts the linear growth law next to it.""", where="Ch. 5 §5.8")

core("C02", r"The Kelvin–Helmholtz dispersion relation $c=\frac{\rho_2U_2+\rho_1U_1}{\rho_2+\rho_1}\pm\Big[\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk-\frac{\rho_2\rho_1}{(\rho_2+\rho_1)^2}(U_2-U_1)^2\Big]^{1/2}$ (11.18)",
     "When does a velocity jump across an interface make waves grow, and which ones?")
problem(r"""
Blow gently across a glass of water: nothing. Blow harder: ripples appear and grow. In the sky, a fast layer of air over a slow,
denser one sometimes paints a row of breaking-wave "billow clouds". In the ocean thermocline, dye released by divers rolls up
into the same shapes. All are one instability: shear across an interface. We want the threshold — how fast must the wind be? —
and the wavelength that grows.""")
idea(r"""
      U₁, ρ₁ (light)   ───────────►            over a crest the upper stream is squeezed → faster → lower pressure
   ~~~~~~~~~~~~~~~~~ interface ζ = ζ₀ e^{ik(x−ct)} ~~~~   ⇒ the crest is SUCKED up (Bernoulli): shear destabilises
      U₂, ρ₂ (heavy)   ──►                     gravity pulls the crest back down: stratification stabilises
   c = mean speed ± √( gravity term − shear term )   →  √(negative) = imaginary ⇒ one root grows
""", r"""
Everything is decided by the sign under the square root: restoring (gravity, later surface tension) minus driving (the kinetic
energy of the shear).""")
remind("C02")
note("N07 [B]", r"""
**The set-up** (figure below): upper fluid speed $U_1$ [m/s], density $\rho_1$ [kg/m³]; lower $U_2$, $\rho_2$; interface
$z=\zeta(x,t)$ [m]; gravity $g$ [m/s²]; inviscid, irrotational in each layer, no surface tension (added in N120).""")
fig(r"""
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.4))
kh_sketch(a)                                                       # drawing only: two streams and a wavy interface
a.set_title("(a) two streams, one interface", fontsize=10)
X, Z = np.meshgrid(np.linspace(0, 4*np.pi, 160), np.linspace(-1.2, 1.2, 81))   # a grid of (x, z) points [m]
U1, U2 = 1.0, -0.6                                                 # the two stream speeds [m/s]
f = ch11.kh_fields(X, Z, 0.0, 1.0, U1, U2, 1.0, 3.0, g=10.0, zeta0=0.08)   # one neutral mode: k = 1 1/m, densities 1 over 3 kg/m^3
c = f["c"].real                                                    # its (real) wave speed [m/s]
rel = np.where(Z > 0, np.sign(U1 - c), np.sign(U2 - c))*f["u"]     # change of flow SPEED seen by an observer riding on the wave [m/s]
xs = np.linspace(0, 4*np.pi, 300)                                  # fine x for the interface line [m]
zeta = ch11.kh_fields(xs, 0.0, 0.0, 1.0, U1, U2, 1.0, 3.0, g=10.0, zeta0=0.08)["zeta"]   # interface height [m]
pc = b.pcolormesh(X, Z, rel, cmap="RdBu_r", vmin=-0.25, vmax=0.25, shading="auto")
b.plot(xs, zeta, color=COLORS["ink"])
fig.colorbar(pc, ax=b, label="speed change [m/s]")
b.set(xlabel="$x$ [m]", ylabel="$z$ [m]", ylim=(-1.2, 1.2))
b.set_title("(b) riding on the wave: faster (red), slower (blue)", fontsize=10)
plt.show()
print(f"wave speed of this mode c = {c:.3f} m/s (real: a neutral wave)")
""",
    see="(a) Two streams (orange arrows) and a wavy interface. (b) For one normal mode, how much faster (red) or slower (blue) "
        "the fluid moves past an observer who rides on the wave; the black line is the interface. The pattern is strongest at "
        "the interface and fades away from it like $e^{-k\\lvert z\\rvert}$.",
    read="Above a crest the upper stream is squeezed and speeds up (red); under the same crest the lower layer is thicker and "
         "slows down (blue). By Bernoulli the pressure is lower above the crest and higher below it: both push the crest "
         "further up. That is the destabilising effect of shear; gravity pulls the other way.",
    change="…$U_1=U_2$: both layers move past the wave at the same speed, the two pressure changes are those of an ordinary "
           "interface wave of Ch. 7, and nothing can grow.",
    explain=r"""
1. `kh_sketch(ax)` is a drawing helper (no physics).
2. `ch11.kh_fields(x, z, t, k, U1, U2, rho1, rho2, g, zeta0)` returns the real fields of one normal mode — the potentials of
   (11.15) with the amplitudes that D03 derives — as a dictionary: `zeta`, `u`, `w` (the disturbance velocity, taken from the
   upper potential for $z>0$ and from the lower for $z<0$) and the wave speed `c`.
3. Riding on the wave, the upper stream passes at $U_1-c$ and the lower at $U_2-c$; multiplying the disturbance $u$ by the sign
   of that relative velocity gives the change of *speed* (`np.where` picks the layer).""")
D("D02", ref="11.13")
note("N13 [B]", r"""
**The normal-mode trial** (C01's mode with $\mathbf K=(k,0,0)$): $\phi_j=A_j(z)\exp\{\mathrm ik(x-ct)\}$ (11.14) and
$\zeta=\zeta_o\exp\{\mathrm ik(x-ct)\}$; inserted into Laplace's equation it gives R07's $A''=k^2A$.""")
P("P256", "np.roots and np.lib.scimath.sqrt", r"""
`np.roots(coeffs)` returns all roots of a polynomial given its coefficients, highest power first (complex if needed).
`np.sqrt(-1.0)` gives `nan` with a warning, but `np.lib.scimath.sqrt(-1.0)` returns the principal complex root `1j` — what we
want when the quantity under a square root can turn negative.""", code=r"""
import numpy as np                              # arrays
print(np.roots([1, -3, 4]))                     # c^2 - 3c + 4 = 0  ->  1.5 + 1.3229j and 1.5 - 1.3229j
print(np.lib.scimath.sqrt(-1.75))               # 1.3229j: the principal square root of a negative number
""")
D("D03", ref="11.18", check_src=r"""
out = ch11.kh_sympy()                               # builds (11.16)-(11.17) with sympy, eliminates A-, A+ and solves for c
print("the quadratic for c:", out["quadratic"])     # D03 step 10: (rho1 + rho2) c^2 - 2 (rho1 U1 + rho2 U2) c + ... = 0
print("its roots minus (11.18):", out["residual_11_18"])   # [0, 0]: the two roots ARE (11.18)
print("discriminant identity (step 12):", out["identity"]) # 0: (rho1 U1 + rho2 U2)^2 - (rho1 + rho2)(rho1 U1^2 + rho2 U2^2) = -rho1 rho2 (U1 - U2)^2
print("limits (step 14): static", out["static"], "| vortex sheet", out["vortex_sheet"])   # [0, 0] each: (11.19) and (11.20)
""")
note("N16 [B]", r"""
**The stability boundary.** Unstable exactly when the square root's argument is negative:
$\frac{\rho_2-\rho_1}{\rho_2+\rho_1}\frac gk<\frac{\rho_2\rho_1}{(\rho_2+\rho_1)^2}(U_2-U_1)^2$, i.e.
$g(\rho_2^2-\rho_1^2)<k\rho_1\rho_2(U_2-U_1)^2$ — unstable if the velocity jump is large, the density jump small, or the wave
short: $k>k_c=g(\rho_2^2-\rho_1^2)/(\rho_1\rho_2\,\Delta U^2)$ with $\Delta U=U_2-U_1$.""")
code(r"""
k_c = ch11.kh_critical_k(5.0, 0.0, 1.2, 1000.0)            # air (1.2 kg/m^3) at 5 m/s over still water (1000 kg/m^3): k_c [1/m]
print(f"k_c = {k_c:.0f} 1/m  ->  every wave shorter than {100*2*np.pi/k_c:.2f} cm grows (no surface tension)")
""")
note("N17 [B]", r"""
**Growing and decaying twins.** The quadratic for $c$ has real coefficients, so its two roots are complex conjugates whenever
they are not real: $c_+=c_-^*$. For every growing wave there is a decaying one; a nonzero $c_i$ therefore always means instability
(stated; the explainer below shows the pair on the $c$-plane).""")
note("N18 [B]", r"""
**Short waves always lose.** Without surface tension $k_c$ is finite for any $\Delta U\neq0$, so all waves shorter than
$2\pi/k_c$ grow — the flow is always unstable to short waves when $U_1\neq U_2$ — and the growth rate $kc_i$ rises without bound
as $k$ grows (slider figure F1). Surface tension (N120) cures it.""")
nb.worked_example("a two-layer flow with easy numbers", r"""
$\rho_1=1$, $\rho_2=3$ kg/m³ (light over heavy), $g=10$ m/s², $k=1$ m⁻¹, $U_2=0$.

**Case $U_1=4$ m/s.**
1. Mean speed $(\rho_2U_2+\rho_1U_1)/(\rho_1+\rho_2)=4/4=1$ m/s.
2. Gravity term $\frac{3-1}{4}\cdot\frac{10}{1}=5$ m²/s².
3. Shear term $\frac{3\cdot1}{16}\cdot16=3$ m²/s².
4. Under the root $5-3=2>0$ ⇒ $c=1\pm\sqrt2=2.414$ or $-0.414$ m/s: two real speeds, neutral waves.

**Case $U_1=6$ m/s.**
5. Mean $6/4=1.5$ m/s.
6. Shear term $\frac3{16}\cdot36=6.75$.
7. Under the root $5-6.75=-1.75<0$ ⇒ $c=1.5\pm1.3229\,\mathrm i$ m/s.
8. Growth rate $kc_i=1.3229$ s⁻¹: the amplitude grows by $e$ every 0.76 s.
9. Check with the inequality: $g(\rho_2^2-\rho_1^2)=10\times8=80$; $k\rho_1\rho_2\Delta U^2=3\times16=48$ (stable) and
   $3\times36=108$ (unstable) ✓; $k_c=80/108=0.741$ m⁻¹ for $\Delta U=6$ m/s.""")
code(r"""
cp, cm = ch11.kh_phase_speed(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0)        # (11.18) with k = 1, U1 = 6, U2 = 0, rho1 = 1, rho2 = 3
print("U1 = 6 m/s: c =", np.round(cp, 4), "and", np.round(cm, 4), "m/s")   # 1.5 + 1.3229i and 1.5 - 1.3229i
print(ch11.kh_discriminant_terms(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0))   # the pieces under the root [m^2/s^2] and the mean speed [m/s]
print("U1 = 4 m/s: c =", np.round(ch11.kh_phase_speed(1.0, 4.0, 0.0, 1.0, 3.0, g=10.0), 4), "m/s (both real: neutral)")
k = np.logspace(1, 4, 400)                                           # wavenumbers from 10 to 10^4 1/m (logarithmically spaced)
growth = ch11.kh_growth_rate(k, 5.0, 0.0, 1.2, 1000.0)               # k c_i [1/s] for air at 5 m/s over still water
print(f"growth is zero up to k = {k[np.nonzero(growth)[0][0]]:.0f} 1/m (k_c = {k_c:.0f}), then rises to {growth[-1]:.0f} 1/s at k = 1e4")
""", explain=r"""
1. `ch11.kh_phase_speed(k, U1, U2, rho1, rho2, g)` returns the two roots of (11.18), the one with $c_i\ge0$ first.
2. `ch11.kh_discriminant_terms` returns the pieces under the square root — `gravity` 5, `tension` 0, `shear` −6.75, their `total`
   −1.75 — and the `mean` speed 1.5: the worked example, term by term (these are the bars of the explainer below).
3. The neutral case $U_1=4$: two real speeds, 2.4142 and −0.4142.
4. `ch11.kh_growth_rate` gives $kc_i$ for a whole array of $k$: zero below $k_c$, then growing without bound (N18).""")
scratch(r"""
mine = np.roots([4.0, -2*6.0, 36.0 - 10.0*2.0])                    # (rho1+rho2) c^2 - 2 (rho1 U1 + rho2 U2) c + rho1 U1^2 + rho2 U2^2 - (g/k)(rho2-rho1)
assert np.allclose(sorted(mine, key=np.imag), sorted([cp, cm], key=np.imag))   # same two roots as (11.18)
c_static = np.real(ch11.kh_phase_speed(2.0, 0.0, 0.0, 1.0, 3.0, g=9.81)[0])    # U1 = U2 = 0: the static interface (11.19)
assert np.isclose(c_static, WV.interface_omega(2.0, 1.0, 3.0, g=9.81)/2.0)     # Ch. 7's omega/k for the same two layers (recap R09)
res = ch11.kh_residuals(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0)           # the computed mode put back into (11.3), (11.9), (11.13)
assert max(res.values()) < 1e-12                                   # all three residuals vanish
print("quadratic by np.roots:", np.round(mine, 4), "| residuals:", {n: f"{v:.1e}" for n, v in res.items()})
""", r"""
1. D03 step 10's quadratic typed with our numbers ($\rho_1+\rho_2=4$, $\rho_1U_1+\rho_2U_2=6$, $\rho_1U_1^2+\rho_2U_2^2=36$,
   $(g/k)(\rho_2-\rho_1)=20$) and solved with `np.roots` (P256): the same roots.
2. With both streams at rest the function reproduces Ch. 7's interface wave, $c=\omega/k$.
3. `kh_residuals` substitutes the mode into Laplace's equation and the two linearised interface conditions at random points.""")
fig(r"""
k = np.logspace(1.3, 3.5, 300)                                     # wavenumbers 20 ... 3000 1/m
dU0 = ch11.kh_stability_boundary(k, 1.2, 1000.0)                   # smallest unstable shear without surface tension [m/s]
dUs = ch11.kh_stability_boundary(k, 1.2, 1000.0, surface_tension=0.074)   # ... with surface tension 0.074 N/m
ms = ch11.kh_min_shear(1.2, 1000.0)                                # the minimum of that curve (D04): dU_min, k_star, wavelength
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.6))
a.fill_between(k, dUs, 14, color=C_GROW, alpha=0.15)
a.semilogx(k, dU0, "--", color=C_SHEAR, label="no surface tension")
a.semilogx(k, dUs, color=C_NEUT, label=r"with $\sigma_s$ = 0.074 N/m")
a.plot(ms["k_star"], ms["dU_min"], "D", color=C_NEUT)              # the diamond: the weakest wind that can raise any wave (D04)
a.text(ms["k_star"]*1.15, ms["dU_min"] - 1.1, f"{ms['dU_min']:.2f} m/s\n$\\lambda$ = {100*ms['wavelength']:.2f} cm", fontsize=8)
for dU in (5.0, 8.0):                                              # our two winds
    a.axhline(dU, color=C_BASE, lw=1, ls=":")
    a.text(22, dU + 0.15, f"$\\Delta U$ = {dU:.0f} m/s", fontsize=8, color=C_BASE)
a.text(900, 12, "unstable", color=C_GROW)
a.set(xlabel="wavenumber $k$ [m$^{-1}$]", ylabel=r"velocity jump $\Delta U$ [m/s]", ylim=(0, 14))
a.set_title("(a) air over water: which waves grow", fontsize=10)
a.legend(fontsize=8, loc="lower left")
k1, k2 = ch11.kh_unstable_band(8.0, 1.2, 1000.0, surface_tension=0.074)   # the unstable wavenumbers at dU = 8 m/s [1/m]
b.axvspan(k1, k2, color=C_GROW, alpha=0.12)
b.semilogx(k, ch11.kh_growth_rate(k, 8.0, 0.0, 1.2, 1000.0), "--", color=C_SHEAR, label="no surface tension")
b.semilogx(k, ch11.kh_growth_rate(k, 8.0, 0.0, 1.2, 1000.0, surface_tension=0.074), color=C_GROW, label="with surface tension")
b.set(xlabel="wavenumber $k$ [m$^{-1}$]", ylabel="growth rate $kc_i$ [s$^{-1}$]", ylim=(0, 90))
b.set_title(r"(b) growth rate at $\Delta U$ = 8 m/s", fontsize=10)
b.legend(fontsize=8)
savefig(fig, "ch11", "c02_kh_boundary")
print(f"band at 8 m/s: k from {k1:.0f} to {k2:.0f} 1/m = wavelengths {100*2*np.pi/k2:.2f} to {100*2*np.pi/k1:.2f} cm")
print("thermocline (densities 1025 and 1027.05 kg/m^3, dU = 0.1 m/s): k_c =", round(ch11.kh_critical_k(0.1, 0.0, 1025.0, 1027.05), 2), "1/m")
""",
    see="(a) The smallest velocity jump that makes a wave of wavenumber $k$ grow: a falling line without surface tension (orange "
        "dashed), a V-shaped curve with it (purple) whose lowest point ◆ is at 6.70 m/s; the rose region is unstable. (b) The "
        "growth rate at 8 m/s: without surface tension it rises forever, with it only a band of wavelengths grows (shaded).",
    read="A horizontal wind line that enters the rose region means: waves of those lengths grow. At 5 m/s the line stays below "
         "the purple curve — with surface tension nothing grows. At 8 m/s it cuts the curve twice: the band of panel (b).",
    change="…the densities differ by only 0.2 % (an ocean thermocline, no surface tension) and the jump is 0.1 m/s: the last "
           "printed line gives $k_c\\approx3.9$ m⁻¹, so every wave shorter than about 1.6 m grows.",
    explain=r"""
1. `ch11.kh_stability_boundary(k, rho1, rho2, surface_tension=…)` solves "discriminant = 0" for the velocity jump.
2. `ch11.kh_min_shear` returns the bottom of the V (derived in D04 below); `ch11.kh_unstable_band` the two edges of the unstable
   band at a given jump.
3. `semilogx` uses a logarithmic $k$ axis (Ch. 1, P13), because wavelengths from 30 cm to 2 mm matter here.""")
nb.md(r"""
**Slider figure F1 — the growth curve as the wind rises.** Every slider position is computed in advance (`slider_figure`,
Ch. 1 P17), so it also works on the web page.""")
nb.plotly(r"""
kF = np.logspace(1.3, 3.5, 160)                                    # wavenumbers [1/m]

def f1(dU):                                                        # the two growth curves for one velocity jump dU [m/s]
    return {"no surface tension": (kF, ch11.kh_growth_rate(kF, dU, 0.0, 1.2, 1000.0)),                       # k c_i(k) from (11.18)
            "with surface tension 0.074 N/m": (kF, ch11.kh_growth_rate(kF, dU, 0.0, 1.2, 1000.0, surface_tension=0.074))}   # ... with the tension term

figF1 = slider_figure(f1, "ΔU", np.linspace(1.0, 12.0, 12 if FAST else 23), unit="m/s", xlabel="wavenumber k [1/m]",   # slider: the velocity jump
                      ylabel="growth rate k c_i [1/s]", title="Air over water: surface tension holds the short waves until 6.7 m/s",   # labels
                      yrange=(0, 80))                              # a fixed vertical range, so the curve is seen to rise
figF1.update_xaxes(type="log")                                     # logarithmic k axis
figF1.show()
""", explain=r"""
1. `f1(dU)` returns the two curves for one slider value; `slider_figure` calls it for every value now.
2. Drag ΔU: below 6.7 m/s the surface-tension curve is flat zero; above it a hump appears near $k\approx360$ m⁻¹ and widens.""")
D("D04", ref="11.18")
note("N120 [B]", r"""
**Depth and surface tension (Exercise 11.1).** With a lower layer of depth $h$ and surface tension $\sigma_s$ [N/m]:
$c=\frac{\rho_1U_1+\rho_2U_2\coth kh}{\rho_1+\rho_2\coth kh}\pm\Big[\frac{(g/k)(\rho_2-\rho_1)+\sigma_sk}{\rho_1+\rho_2\coth kh}-
\frac{\rho_1\rho_2(U_1-U_2)^2\coth kh}{(\rho_1+\rho_2\coth kh)^2}\Big]^{1/2}$; $h\to\infty$ gives the deep-water relation plus the
$\sigma_sk$ term. Minimum wind for deep water: $\Delta U_{\min}^2=2\sqrt{g\Delta\rho\,\sigma_s}(\rho_1+\rho_2)/(\rho_1\rho_2)$ at
$k_*=\sqrt{g\Delta\rho/\sigma_s}$ (D04). Real wind waves start at lower winds: the turbulent wind's pressure fluctuations (the
mechanisms of Miles and Phillips, named only) do the job Kelvin–Helmholtz cannot.""")
code(r"""
print(ch11.kh_min_shear(1.2, 1000.0))                               # air over water: dU_min [m/s], k_star [1/m], wavelength [m]
print(ch11.kh_unstable_band(8.0, 1.2, 1000.0, surface_tension=0.074))   # unstable k-interval at 8 m/s [1/m]
print(ch11.kh_unstable_band(5.0, 1.2, 1000.0, surface_tension=0.074))   # at 5 m/s: (nan, nan) = no unstable wave
c_deep = ch11.kh_phase_speed(300.0, 8.0, 0.0, 1.2, 1000.0, surface_tension=0.074)[0]           # deep water, k = 300 1/m
c_shal = ch11.kh_phase_speed(300.0, 8.0, 0.0, 1.2, 1000.0, surface_tension=0.074, h=0.002)[0]  # the same over water only 2 mm deep
print("growth rate k c_i: deep", round(300*c_deep.imag, 2), "1/s | 2 mm of water", round(300*c_shal.imag, 2), "1/s")
""", explain=r"""
1. The minimum wind 6.70 m/s at $k_*=364$ m⁻¹ (wavelength 1.73 cm) — D04's result with numbers.
2. The band at 8 m/s is about 149–887 m⁻¹ (0.71–4.21 cm); at 5 m/s there is none (`nan` = "not a number" marks "no band").
3. The keyword `h` switches on the finite depth of the lower layer (the $\coth kh$ factors).""")
note("N121 [B]", r"""
**Rayleigh–Taylor with surface tension (Exercise 11.2, our answer).** Heavy fluid above (paint on a ceiling): with $\rho_1>\rho_2$
every wave grows unless surface tension holds the short ones; the longest wave that stays neutral is
$\lambda_c=2\pi\sqrt{\sigma_s/((\rho_h-\rho_l)g)}$ — $2\pi$ times the *capillary length* $\sqrt{\sigma_s/((\rho_h-\rho_l)g)}$, the length at which surface tension and gravity are equally strong (about 2.7 mm for water in air) ($\rho_h$, $\rho_l$: heavy and light
densities). This is the cut-off of C01's figure.""")
code(r"""
lam_c = ch11.rayleigh_taylor_cutoff(0.074, 1000.0, 1.2)             # water above air, surface tension 0.074 N/m: cut-off wavelength [m]
print(f"lambda_c = {100*lam_c:.2f} cm: a water film on a ceiling holds only over patches smaller than this")
""")
note("N19 [B]", r"""
**The moving frame.** Seen from a frame moving at $(U_1+U_2)/2$ (Ch. 3's change of frame), the two streams are $\pm\Delta U/2$ —
a symmetric picture with no preferred direction, so the unstable wave cannot travel: for $\rho_1=\rho_2$ the real part of $c$ is
exactly the mean speed, $c_r=(U_1+U_2)/2$. With different densities the density-weighted mean
$(\rho_2U_2+\rho_1U_1)/(\rho_2+\rho_1)$ takes its place.""")
fig(r"""
U1, U2 = 3.0, 1.0                                                  # upper and lower stream speeds in the laboratory frame [m/s]
u_lab = lambda x, t: np.array([U1 if x[1] > 0 else U2, 0.0])      # the lab-frame velocity field: U1 above z = 0, U2 below
u_mov = KIN.galilean_transform(u_lab, np.array([(U1 + U2)/2, 0.0]))   # Ch. 3: the same flow seen from a frame moving at the mean speed
zz = np.linspace(-1, 1, 9)                                         # heights at which we draw arrows [m]
fig, axs = plt.subplots(1, 2, figsize=(7.6, 2.8), sharey=True)
for ax, fn, ttl in ((axs[0], u_lab, "laboratory frame"), (axs[1], u_mov, "frame moving at $(U_1+U_2)/2$")):
    uu = np.array([fn(np.array([0.0, z]), 0.0)[0] for z in zz if z != 0])   # horizontal speed at each height
    ax.quiver(np.zeros(len(uu)), [z for z in zz if z != 0], uu, 0*uu, angles="xy", scale_units="xy", scale=1, color=C_SHEAR)
    ax.axhline(0, color=COLORS["ink"])
    ax.axvline(0, color=C_BASE, lw=0.8)
    ax.set(xlim=(-1.5, 3.5), xlabel="speed [m/s]")
    ax.set_title(ttl, fontsize=10)
axs[0].set_ylabel("$z$ [m]")
plt.show()
print("c for equal densities:", ch11.vortex_sheet_c(U1, U2))       # (11.20): mean 2 m/s, plus/minus i (U2 - U1)/2
""",
    see="Left: both streams move to the right, 3 m/s above and 1 m/s below. Right: the same flow seen by an observer moving at "
        "2 m/s — +1 m/s above, −1 m/s below.",
    read="In the moving frame the picture is mirror-symmetric, so a growing wave has no reason to drift either way: $c_r=0$ "
         "there, i.e. $c_r=2$ m/s in the laboratory — exactly the real part printed.",
    change="…the lower layer is three times denser: the wave is carried more by the heavy layer, and $c_r$ becomes the "
           "density-weighted mean $(3\\times1+1\\times3)/4=1.5$ m/s.",
    explain=r"""
1. `KIN.galilean_transform(u, U)` (Ch. 3, P96) returns the velocity field seen from a frame moving at the constant velocity `U`.
2. `ch11.vortex_sheet_c(U1, U2)` is the equal-density limit (11.20).""")
note("N20 [C]", r"""
**Shear against stratification in nature.** Laboratory tilting-tube experiments (Thorpe 1971), billow clouds and dye in ocean
thermoclines (Woods 1969) all show this instability; it is a main source of internal waves and of mixing across density
interfaces. With *continuous* stratification the question becomes the Richardson number — C09. Climate hook: Kelvin–Helmholtz
billows in the thermocline and in the stable night-time atmosphere are why mixing parameterisations switch on when the Richardson
number drops below about ¼ (Ch. 12, Ch. 13).""")
nb.md(r"""
**Animation A2 — linear growth against the real roll-up.** Left: linear theory for a vortex sheet ($\rho_1=\rho_2$, jump
$\Delta U=1$, wavelength 1, so $k=2\pi$ and the growth rate is $k\Delta U/2=\pi$). We start from a pure displacement, which is an
equal mix of the growing mode and its decaying twin (N17), so the linear amplitude is $0.01\cosh(\pi t)$. Right: Ch. 5's point-vortex
sheet started from the same displacement.""")
nb.animation(r"""
nf = 14 if FAST else 18                                            # number of frames
tt = np.linspace(0.0, 2.0, nf)                                     # times [s] (wavelength 1 m, jump 1 m/s)
roll = ch05.sheet_rollup(N=80 if FAST else 100, gamma=1.0, amplitude=0.01, delta=0.05, t_eval=tt)   # the nonlinear sheet (Ch. 5)
a_lin = 0.01*np.cosh(np.pi*tt)                                     # linear amplitude: growing + decaying mode, rate k dU/2 = pi [1/s]
a_num = np.max(np.abs(roll["y"]), axis=1)                          # the largest displacement of the nonlinear sheet [m]
xs = np.linspace(0, 1, 200)                                        # one wavelength [m]
fig, (a, b) = plt.subplots(1, 2, figsize=(8.4, 3.2), sharey=True)
(lin,) = a.plot(xs, 0.01*np.sin(2*np.pi*xs), color=C_GROW)         # the linear interface
(pts,) = b.plot(roll["x"][0] % 1.0, roll["y"][0], ".", ms=2.5, color=C_SHEAR)   # the point vortices
(env,) = b.plot(xs, 0.01*np.sin(2*np.pi*xs), "--", color=C_BASE, lw=1)          # the linear interface again, as a ghost
a.set(xlim=(0, 1), ylim=(-0.3, 0.3), xlabel="$x$ [m]", ylabel="$z$ [m]")
a.set_title("linear theory", fontsize=10)
b.set(xlim=(0, 1), xlabel="$x$ [m]")
ttl = b.set_title("point-vortex sheet (Ch. 5)", fontsize=10)

def update(i):                                                     # frame i: time tt[i]
    y_lin = np.clip(a_lin[i]*np.sin(2*np.pi*xs), -0.3, 0.3)        # linear interface, cut at the panel edge
    lin.set_ydata(y_lin)
    env.set_ydata(y_lin)
    pts.set_data(roll["x"][i] % 1.0, roll["y"][i])                 # move every vortex (positions folded into one wavelength)
    ttl.set_text(f"t = {tt[i]:.2f} s: linear {a_lin[i]:.3f} m, sheet {a_num[i]:.3f} m")
    return lin, env, pts, ttl

show_animation(animate(update, frames=nf, fig=fig, interval=110), player="video", dpi=60)
i10 = int(np.argmax(a_num > 0.05))                         # first frame with a sheet amplitude above 5 % of the wavelength
print(f"at t = {tt[i10]:.2f} s: linear {a_lin[i10]:.3f} m vs sheet {a_num[i10]:.3f} m; at the end: {a_lin[-1]:.2f} m vs {a_num[-1]:.3f} m")
""", explain=r"""
1. `ch05.sheet_rollup(N, gamma, amplitude, delta, t_eval)` moves `N` point vortices of a periodic sheet (strength `gamma` = the
   velocity jump) from a sine displacement; `delta` smooths the singular kernel.
2. The linear amplitude is $0.01\cosh(\pi t)$: it grows exponentially forever.
3. The printed line compares the two when the sheet's amplitude passes 5 % of the wavelength, and at the end.""")
nb.figure_notes(
    see="Two growing waves that agree at first; then the left one shoots off the panel while the right one folds over and winds "
        "into a spiral (a cat's eye).",
    read="Linear theory is right until the amplitude is about a tenth of the wavelength; after that the slope is no longer small, "
         "the dropped quadratic terms matter, and the sheet rolls up instead of growing forever. (The sheet lags the ideal law "
         "slightly even early on, because the smoothing δ slows short-wave growth.)",
    change="…we halve δ (less smoothing): the roll-up core tightens and forms earlier; the linear stage is the same.")
note("N21 [C] · N22 [B]", r"""
**Where the energy comes from.** The billows feed on the kinetic energy of the two streams. A mixing example: a stream $U_1$ over
still fluid, both of thickness $h$ and the same density $\rho$, mixed into a linear profile $U(z)=U_1(\tfrac12+\tfrac z{2h})$,
$-h\le z\le h$, keeps its momentum $\int\rho U\,dz=\rho U_1h$ but its kinetic energy drops from $E_i=\tfrac\rho2U_1^2h$ to
$E_f=\tfrac\rho2\int_{-h}^hU^2dz=\tfrac\rho3U_1^2h$ — a third is released, available to lift heavy fluid (mixing a stably
stratified interface raises the potential energy).""")
code(r"""
mix = ch11.kh_mixing_energy(1.0, 1.0, 1000.0)                       # U1 = 1 m/s, h = 1 m, rho = 1000 kg/m^3
E_f_quad = 0.5*1000.0*quad(lambda z: (1.0*(0.5 + z/2.0))**2, -1.0, 1.0)[0]   # the integral of U^2 done numerically (quad, Ch. 3 P87)
print({n: round(v, 3) for n, v in mix.items()})                     # energies [J/m^2], momenta [kg/(m s)], ratio E_f/E_i
assert np.isclose(E_f_quad, mix["E_f"])                             # quadrature = closed form rho U1^2 h / 3
fig, ax = plt.subplots(figsize=(4.6, 2.8))
ax.bar(["before $E_i$", "after $E_f$", "released"], [mix["E_i"], mix["E_f"], mix["E_i"] - mix["E_f"]], color=[C_SHEAR, C_BASE, C_GROW])
ax.set(ylabel="kinetic energy per unit area [J/m$^2$]")
ax.set_title("mixing a shear layer releases a third", fontsize=10)
plt.show()
""", explain=r"""
1. `ch11.kh_mixing_energy(U1, h, rho)` returns the kinetic energy before (`E_i`) and after (`E_f`) mixing into the linear profile,
   their `ratio` (2/3), and the momentum before and after (`M_i` = `M_f`: mixing conserves momentum).
2. `scipy.integrate.quad` integrates $U^2$ numerically as a cross-check.
3. The bar chart: one third of the kinetic energy is released.""")
note("N23 [C]", r"""
**A general rule** (named): for a fixed momentum $\int U\,dz$, smoothing the velocity gradients always lowers $\int U^2dz$
(Cauchy–Schwarz: $\int_{-h}^hU^2dz\ge(\int_{-h}^hU\,dz)^2/(2h)$, with equality only for a uniform flow), so mixing a shear layer
always releases kinetic energy. C14 shows how a growing wave taps that energy (the Reynolds stress).""")
explainer("kelvin_helmholtz_boundary", "Why does wind raise waves only above a threshold?", r"""
Two sliders (velocity jump, density ratio) and two toggles (surface tension, finite depth) move both the stability boundary and
the current point; on the $c$-plane the two real wave speeds slide together, collide and split into a growing/decaying pair
exactly when the point crosses the boundary — the collision is the instability, and only motion shows it.""", [
    "Preset 'air over water 5 m/s': stable. Drag ΔU up and stop when the status turns rose — read ΔU_min ≈ 6.7 m/s and λ ≈ 1.7 cm.",
    "Watch the c-plane while you cross the boundary: two dots on the real axis meet and leave it vertically.",
    "Turn surface tension off: the boundary becomes a falling line — every short wave grows (N18).",
])
nb.pointer(r"""**S01** — Exercises 11.3–11.5 (a porous surface, a compliant surface, membrane flutter) extend the interface
conditions of C02 to other boundaries; they are not solved here.""")
whatif(r"""
…instead of a velocity jump the fluid were at rest and heated from below? Then there is no shear to tap; the energy source is
buoyancy, and viscosity and heat diffusion fight it. The answer is a threshold number, the Rayleigh number — C03.""")
# =====================================================================================================================
# A.4 §11.4 Thermal Instability: The Bénard Problem — R12 R13, C03 · C04 · C05
# =====================================================================================================================
nb.section("11.4", "Thermal Instability: The Bénard Problem", intro=r"""
**What is this section about?** Heat a layer of fluid from below: the bottom fluid expands, becomes lighter, and wants to rise —
yet a thin layer of honey on a warm plate stays still. Viscosity and heat diffusion resist; only when the heating beats them,
measured by one number, the Rayleigh number, does the layer overturn into convection cells. We set up the linear problem (C03),
find the threshold for rigid walls, Ra ≈ 1708 (C04), and solve the stress-free case by hand (C05). The same physics drives shallow
atmospheric convection on a sunny day, cloud streets, and convection in the ocean's mixed layer.""")
nb.recap("R12", "The Boussinesq equations", r"""
Density changes are kept only where they meet gravity (Ch. 4 §4.9): $\nabla\cdot\tilde{\mathbf u}=0$,
$\frac{\partial\tilde{\mathbf u}}{\partial t}+(\tilde{\mathbf u}\cdot\nabla)\tilde{\mathbf u}=-\frac1{\rho_0}\nabla\tilde p-g[1-\alpha(\tilde
T-T_0)]\mathbf e_z+\nu\nabla^2\tilde{\mathbf u}$, $\frac{\partial\tilde T}{\partial t}+(\tilde{\mathbf u}\cdot\nabla)\tilde T=\kappa\nabla^2\tilde T$
— Ch. 4's continuity, momentum and heat equations (4.10, 4.86, 4.89) written for the total fields (a tilde marks "total"), with
$\rho=\rho_0[1-\alpha(\tilde T-T_0)]$, $\alpha$ the thermal expansion coefficient [1/K], $\nu$ the kinematic viscosity and $\kappa$ the
thermal diffusivity [m²/s]. They are the starting line of D05.""", where="Ch. 4 §4.9")
nb.recap("R13", "The conduction state", r"""
With no motion the momentum equation is hydrostatic and the heat equation says the temperature profile has no curvature:
$0=-\frac1{\rho_0}\nabla P-g[1-\alpha(\bar T-T_0)]\mathbf e_z$ and $0=\kappa\frac{\partial^2\bar T}{\partial z^2}$ (11.23) — Ch. 1's
hydrostatics and Ch. 8's steady conduction: a straight temperature line between the plates.""", where="Ch. 1 §1.7, Ch. 8")

core("C03", r"The Bénard amplitude problem $\big(\sigma+K^2-\frac{d^2}{dz^2}\big)\hat T=W$ and $\big(\frac\sigma\Pr+K^2-\frac{d^2}{dz^2}\big)\big(\frac{d^2}{dz^2}-K^2\big)W=-\mathrm{Ra}K^2\hat T$ (11.36, 11.37)",
     "A layer heated from below: which disturbances grow, and how fast?")
problem(r"""
A pan of water on a stove, a 5 mm layer of oil on a hot plate, the lowest kilometre of the atmosphere on a sunny morning: warm
fluid underneath, cool on top. A blob nudged upward finds itself warmer and lighter than its new surroundings and keeps rising —
unless viscosity slows it and heat leaks out of it first. We want to know, for a given layer, whether that tug-of-war is won, and
if so how fast the convection starts.""")
idea(r"""
T₀ − ΔT  ─────────────── cold plate
          ↑ warm blob: buoyancy  gαT′           pushes it up      (∝ ΔT, d)
          │ viscosity ν∇²w                       holds it back
          │ heat leaks out: κ∇²T′                erases T′
T₀       ─────────────── hot plate
one number compares them:  Ra = buoyancy / (viscous × diffusive) = gαΓd⁴/(κν)
small disturbance + normal modes in x, y  ⇒  two ODEs in z  ⇒  σ is an EIGENVALUE that depends on (K, Ra, Pr)
""")
remind("C03")
note("N24 [C]", r"""
**History.** Bénard's 1900 hexagons in millimetre-thin layers with a free surface were mostly driven by the variation of surface
tension with temperature (Marangoni convection), not by buoyancy (Drazin & Reid 1981); Rayleigh (1916) solved the buoyancy
problem with free surfaces, Jeffreys (1928) with rigid ones. The patterns themselves are N51.""")
note("N25 [B]", r"""
**The Rayleigh number** $\mathrm{Ra}=g\alpha\Gamma d^4/\kappa\nu$ (11.21), with $g$ [m/s²], $\alpha$ thermal expansion [1/K],
$\Gamma=-d\bar T/dz$ [K/m] (positive when heated from below), $d$ depth [m], $\kappa$ thermal diffusivity [m²/s], $\nu$ kinematic
viscosity [m²/s] — dimensionless. Since $\Gamma=\Delta T/d$ (N28), $\mathrm{Ra}=g\alpha\,\Delta T\,d^3/(\kappa\nu)$.""")
confusion(r"""
**which Γ?** (slip #10) The Rayleigh number $\mathrm{Ra}=g\alpha\Gamma d^4/\kappa\nu$ (11.21) uses $\Gamma\equiv-d\bar T/dz$ — the
*opposite* sign to Ch. 1's Kundu lapse rate $\Gamma\equiv dT/dz$; here Γ happens to agree with the meteorological $-dT/dz$.

| Convention | Definition | Our water layer ($d$ = 5 mm, bottom 2 K warmer) | Unstable side |
|---|---|---|---|
| this section | $\Gamma=-d\bar T/dz$ | $+400$ K/m | $\Gamma>0$ (heated from below) |
| Ch. 1 (Kundu) | $\Gamma=dT/dz$ | $-400$ K/m | $dT/dz<0$ |
| meteorology | $\Gamma_{met}=-dT/dz$ | $+400$ K/m | $\Gamma_{met}>0$ |

The code avoids the trap: it takes `dT` = T_bottom − T_top.""")
code(r"""
Ra_w = ch11.rayleigh_number(alpha=2.1e-4, dT=2.0, d=0.005, kappa=1.4e-7, nu=1.0e-6, g=9.81)    # 5 mm of water, bottom 2 K warmer: (11.21)
print(f"Ra = {Ra_w:.0f}")                                           # about 3679 (dimensionless)
print(ch11.gamma_conventions(2.0, 0.005))                           # the same gradient in the three conventions [K/m]
Ra_wrong = ch11.rayleigh_number(alpha=2.1e-4, dT=-2.0, d=0.005, kappa=1.4e-7, nu=1.0e-6, g=9.81)   # heated from ABOVE (top 2 K warmer)
print(f"heated from above: Ra = {Ra_wrong:.0f} (negative: stable to every disturbance)")
""", explain=r"""
1. `ch11.rayleigh_number(alpha, dT, d, kappa, nu, g)` evaluates (11.21) with `dT` = T_bottom − T_top [K].
2. `ch11.gamma_conventions(dT, d)` prints the one gradient in the three sign conventions: +400, −400, +400 K/m.
3. Reversing the heating flips the sign of Ra — a mistyped sign of Γ would turn a convecting layer into a stable one.""")
note("N26 [B]", r"""
**Geometry**: a layer of depth $d$ between two plates, $z$ measured from the middle ($-d/2\le z\le d/2$), bottom at temperature
$T_0$, top at $T_0-\Delta T$.""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.4, 3.0))
benard_sketch(ax)                                                  # drawing only: the layer and its conduction temperature line
plt.show()
""",
    see="A fluid layer between a hot plate (below) and a cold plate (above), with the conduction temperature profile drawn beside "
        "it: a straight line, hottest at the bottom.",
    read="The slope of that line, taken positive downward, is Γ = ΔT/d: the larger the temperature difference or the thinner the "
         "layer, the steeper the line.",
    change="…we heat from above: the line tilts the other way, Ra < 0, and nothing can grow (the light fluid is already on top).")
note("N27 [B] · N28 [B]", r"""
**Base state plus disturbance:** $\tilde{\mathbf u}=0+\mathbf u(\mathbf x,t)$, $\tilde T=\bar T(z)+T'(\mathbf x,t)$,
$\tilde p=P(z)+p(\mathbf x,t)$ (11.22) — $\mathbf u=(u,v,w)$, $T'$, $p$ are the small disturbances. Integrating the conduction
equation twice between the plates gives the conduction profile $\bar T(z)=T_0-\tfrac12\Delta T-\Gamma z$, $\Gamma\equiv\Delta T/d$
(11.24).""")
code(r"""
base = ch11.benard_base_state(np.array([-0.0025, 0.0, 0.0025]), T0=300.0, dT=2.0, d=0.005)   # bottom, middle, top of a 5 mm layer [m]
print("conduction temperatures [K]:", base["T"])                    # 300, 299, 298 K: (11.24) with T0 = 300 K at the bottom
""")
D("D05", ref="11.27")
note("N31 [B]", r"""
**Where Ra comes from (scaling).** In $\frac{\partial T'}{\partial t}-w\Gamma=\kappa\nabla^2T'$ (11.27) balance
$w\Gamma\sim\kappa\,\Delta T/d^2$ with $T'\sim\Delta T$ and $\nabla\sim1/d$: $w\sim\kappa/d$. Then in the momentum equation the
ratio of buoyancy to viscous force is $\sim\frac{g\alpha T'}{\nu w/d^2}\sim\frac{g\alpha\Gamma d^4}{\nu\kappa}=\mathrm{Ra}$
(order-of-magnitude scaling, Ch. 4 P130).""")
code(r"""
print({n: float(f"{v:.4g}") for n, v in ch11.benard_scales(alpha=2.1e-4, dT=2.0, d=0.005, kappa=1.4e-7, nu=1e-6, g=9.81).items()})   # water layer
Ra_air = ch11.rayleigh_number(alpha=1/293, dT=1.0, d=0.01, kappa=2.1e-5, nu=1.5e-5, g=9.81)   # 1 cm of air, 1 K across it
print(f"air gap of 1 cm: Ra = {Ra_air:.0f} per kelvin -> about {1708/Ra_air:.0f} K needed to reach Ra = 1708")
print(f"water layer of 5 mm: Ra = {Ra_w/2:.0f} per kelvin -> about {1708/(Ra_w/2):.2f} K needed")
Ra_w_1cm = ch11.rayleigh_number(alpha=2.1e-4, dT=1.0, d=0.01, kappa=1.4e-7, nu=1.0e-6, g=9.81)   # water at the SAME depth as the air gap
print(f"at equal depth (1 cm) and equal temperature difference: Ra(water)/Ra(air) = {Ra_w_1cm/Ra_air:.0f}")
""", explain=r"""
1. `ch11.benard_scales` returns Ra, the Prandtl number $\Pr=\nu/\kappa=7.14$ for water, Γ, the velocity scale
   $\kappa/d=2.8\times10^{-5}$ m/s and the time scale $d^2/\kappa=179$ s.
2. Air is much harder to set convecting than water: about 16 K across 1 cm of air against less than 1 K across 5 mm of water
   (the threshold 1708 is found in C04). At equal depth and equal temperature difference the Rayleigh number of water is about
   140 times that of air (the last printed line) — air's large $\kappa\nu$ outweighs its larger expansion coefficient.""")
D("D06", ref="11.29")
code(r"""
chk = ch11.benard_perturbation_sympy()                              # D05-D07 redone with sympy on generic fields (cached after the first call)
print("residual of (11.29) after eliminating the pressure:", chk["eq_11_29_residual"])     # 0
print("planted mistake (full Laplacian instead of the horizontal one) leaves:", chk["planted_residual"])   # not 0
print("amplitude equations (11.36), (11.37) residuals:", chk["eq_11_36"], chk["eq_11_37"])   # 0 0 (used by D07 below)
""", explain=r"""
1. `ch11.benard_perturbation_sympy()` repeats the pressure elimination symbolically; the residual 0 means D06's result is exactly
   what the moves give.
2. The planted variant — forgetting to differentiate the pressure equation in $z$, which leaves $\nabla^2T'$ instead of
   $\nabla_H^2T'$ — leaves a nonzero term, $g\alpha\,\partial^2T'/\partial z^2$: the check can tell right from wrong.""")
note("N35 [B]", r"""
**Rigid isothermal walls:** $w=\partial w/\partial z=T'=0$ on $z=\pm d/2$ (11.30) — no slip gives $u=v=w=0$ on the whole plate, so
$\partial u/\partial x=\partial v/\partial y=0$ there, and continuity forces $\partial w/\partial z=0$.""")
nb.md(r"""
> ⚠️ **Four ways this chapter makes things dimensionless — we are here: §11.4.** Lengths by $d$, time by $d^2/\kappa$, and the
> vertical velocity $w$ is left dimensional until the very last step of D07, where $W\equiv(\Gamma d^2/\kappa)\hat w$ makes Ra
> appear. (§11.6 will use the gap and $d^2/\nu$; §11.8–11.11 a length $L$ and $L/U_0$; §11.14 rescales once more — the table in
> the conventions cell at the top.)""")
D("D07", ref="11.37")
P("P257", "eigenvalue problem for a differential equation", r"""
In Ch. 2 an eigenvalue problem was a matrix equation $A\mathbf v=\lambda\mathbf v$ with a nonzero $\mathbf v$. Here it is an ODE
with boundary conditions, like $-u''=\lambda u$ with $u(0)=u(\pi)=0$: for most $\lambda$ the only solution is $u=0$; only for
special $\lambda$ (here $\lambda=1,4,9,\dots$ with $u=\sin nz$) does a nonzero solution meet every boundary condition. Those
$\lambda$ are the eigenvalues, and the nonzero solution that goes with each one is its *eigenfunction* (its *mode shape*). In the amplitude problem above the growth rate $\sigma$ plays $\lambda$: only special $\sigma$
allow a nonzero $(W,\hat T)$.""", code=r"""
import numpy as np                                    # arrays
z = np.linspace(0, np.pi, 5)                          # five points from 0 to pi
for n in (1, 2):                                      # u = sin(n z) meets u(0) = u(pi) = 0 ...
    print(n**2, np.allclose(np.sin(n*z)[[0, -1]], 0)) # ... and -u'' = n^2 u: eigenvalues 1 and 4, both True
""")
note("N39 [B]", r"""
**Boundary conditions for the amplitudes:** $W=\frac{dW}{dz}=\hat T=0$ on $z=\pm\tfrac12$ (11.38). With the two amplitude
equations this is a sixth-order problem in $z$ whose eigenvalue is $\sigma$.""")
P("P258", "Chebyshev–Gauss–Lobatto points and the differentiation matrix", r"""
To turn an ODE into a matrix problem we sample the unknown at $N+1$ points $x_j=\cos(j\pi/N)$ — crowded near the ends — and
replace $d/dx$ by a matrix $D$ that differentiates the polynomial through those values exactly. Because the points are Chebyshev
points, the error falls *exponentially* with $N$ for smooth solutions (spectral accuracy) — far faster than Ch. 10's second-order
stencils. $D^2$ is `D @ D`; for $z\in[-\tfrac12,\tfrac12]$ the matrix is just rescaled by 2. Note the order: `x[0]` is the right
(top) end.""", code=r"""
import numpy as np                                    # arrays
from fluidpy.core import stability as ST              # the eigen-solver toolkit
D, x = ST.cheb(4)                                     # 5 points: 1, 0.707, 0, -0.707, -1 and the 5 x 5 matrix D
print(np.round(x, 3))                                 # the Chebyshev points, from +1 down to -1
print(np.max(np.abs(D @ x**3 - 3*x**2)))              # ~1e-15: D differentiates polynomials up to degree N exactly
""")
P("P259", "boundary-row replacement and the generalised non-symmetric eigenproblem scipy.linalg.eig(A, B)", r"""
Collocation turns the ODEs into $A\mathbf v=\sigma B\mathbf v$. Boundary conditions replace the rows at the end points: in $A$ the
row becomes the condition (for example $[1,0,\dots,0]$ for $W=0$), in $B$ the row becomes zero — so those rows say "0 = condition"
for every $\sigma$. $B$ is then singular and `scipy.linalg.eig(A, B)` returns some infinite eigenvalues (drop them) and, at larger
$N$, a few huge spurious ones caused by round-off (keep only eigenvalues that do not move when $N$ grows). Unlike Ch. 10's `eigh`
(P247: symmetric matrices, real eigenvalues), $A$ and $B$ are not symmetric here, so eigenvalues can be complex.""", code=r"""
import numpy as np                                    # arrays
from scipy.linalg import eig                          # the generalised eigen-solver
A = np.array([[2.0, 1.0], [1.0, 1.0]])                # a 2 x 2 example
B = np.array([[1.0, 0.0], [0.0, 0.0]])                # its second row is a 'boundary row': zero in B
print(eig(A, B, right=False))                         # [1, inf]: the boundary row gives an infinite eigenvalue
""")
nb.worked_example("a 5 mm water layer on a hot plate", r"""
$d=5$ mm, $\Delta T=2$ K, $\alpha=2.1\times10^{-4}$ K⁻¹, $\kappa=1.4\times10^{-7}$ m²/s, $\nu=1.0\times10^{-6}$ m²/s, $g=9.81$ m/s².

1. $\Gamma=\Delta T/d=400$ K/m.
2. $\mathrm{Ra}=g\alpha\,\Delta T\,d^3/(\kappa\nu)=9.81\times2.1\times10^{-4}\times2\times1.25\times10^{-7}/(1.4\times10^{-7}\times10^{-6})
   =5.15\times10^{-10}/1.4\times10^{-13}\approx3679$.
3. $\Pr=\nu/\kappa=7.14$.
4. Time unit $d^2/\kappa=25\times10^{-6}/1.4\times10^{-7}=179$ s.
5. The code below finds the largest growth rate $\sigma\approx19.9$ (in units of $\kappa/d^2$) at $K=3.12$: dimensional
   $\sigma=19.9/179\ \text{s}\approx0.111$ s⁻¹, so a disturbance grows by $e$ every ≈ 9 s.
6. Halve $\Delta T$ to 1 K: $\mathrm{Ra}=1839$, just above 1708 — the growth almost stops (C04).""")
code(r"""
Ra = ch11.rayleigh_number(alpha=2.1e-4, dT=2.0, d=0.005, kappa=1.4e-7, nu=1e-6, g=9.81)   # the worked example: Ra = 3679 (11.21)
Pr = 1e-6/1.4e-7                                                    # Prandtl number nu/kappa = 7.14 for water
sig = ch11.benard_growth_rate(3.1163, Ra, Pr, all=True)             # every converged growth rate at K = 3.1163 [units kappa/d^2]
print("three largest growth rates:", np.round(sig[:3], 4))
for N in (16, 24, 32, 40):                                          # how many digits does the leading one gain with N?
    print(f"  N = {N}: sigma_1 = {ch11.benard_growth_rate(3.1163, Ra, Pr, N=N):.9f}")
sig_c = ch11.benard_growth_rate(3.1163, Ra, Pr, all=True, return_complex=True)   # the same spectrum as complex numbers
print("largest imaginary part among the five leading ones:", np.max(np.abs(np.imag(sig_c[:5]))))   # 0: all real (D08 proves it)
print(f"e-folding time = (d^2/kappa)/sigma_1 = {0.005**2/1.4e-7/sig[0]:.2f} s;  at Ra = 1000: sigma_1 = {ch11.benard_growth_rate(3.1163, 1000.0, Pr):.2f}")
""", explain=r"""
1. Ra from (11.21) and Pr for the water layer.
2. `ch11.benard_growth_rate(K, Ra, Pr, all=True)` builds the two amplitude equations with the rigid conditions as a generalised
   eigenproblem on Chebyshev points (P258, P259) and returns the growth rates, largest first: one positive, the rest negative.
3. The leading eigenvalue is the same to nine digits from $N=16$ to $N=40$ — spectral accuracy.
4. `return_complex=True` shows the imaginary parts: zero. Convection starts in place, it does not oscillate (D08).
5. In seconds: the disturbance grows by a factor $e$ about every 9 s; at Ra = 1000 the leading growth rate is negative — every
   disturbance decays.""")
scratch(r"""
def cheb_mine(N):                                                   # Trefethen's Chebyshev differentiation matrix on [-1, 1]
    j = np.arange(N + 1)                                            # point indices 0 ... N
    xc = np.cos(np.pi*j/N)                                          # the Chebyshev points x_j = cos(j pi/N)
    cj = np.where((j == 0) | (j == N), 2.0, 1.0)*(-1.0)**j          # weights c_j (-1)^j, with c = 2 at the two ends
    dX = xc[:, None] - xc[None, :]                                  # all differences x_i - x_j
    Dm = np.outer(cj, 1.0/cj)/(dX + np.eye(N + 1))                  # off-diagonal entries c_i (-1)^(i+j) / (c_j (x_i - x_j))
    return Dm - np.diag(Dm.sum(axis=1)), xc                         # diagonal = minus the sum of the row (a constant has zero derivative)

D8, x8 = cheb_mine(8)                                               # 9 points
assert np.allclose(D8, ST.cheb(8)[0]) and np.allclose(D8 @ x8**3, 3*x8**2)   # same matrix as the library; exact for x^3

def benard_spectrum_mine(K, Ra, Pr, N):                             # the amplitude problem (11.36)-(11.38) as A v = sigma B v, v = [W; T]
    Dz = 2.0*cheb_mine(N)[0]                                        # d/dz on z in [-1/2, 1/2]: the interval is half as long, so D doubles
    I = np.eye(N + 1)                                               # identity
    L = Dz @ Dz - K**2*I                                            # the operator d^2/dz^2 - K^2
    A = np.block([[L @ L, -Ra*K**2*I], [I, L]])                     # (11.37): L^2 W - Ra K^2 T = (sigma/Pr) L W;  (11.36): W + L T = sigma T
    B = np.block([[L/Pr, 0*I], [0*I, I]])                           # the terms multiplied by sigma
    for r, row in ((0, I[0]), (N, I[N]), (1, Dz[0]), (N - 1, Dz[N])):   # W = 0 and dW/dz = 0 at both walls (11.38)
        A[r, :] = 0.0; A[r, :N + 1] = row; B[r, :] = 0.0            # row replacement (P259): condition in A, zeros in B
    for r in (N + 1, 2*N + 1):                                      # T = 0 at both walls
        A[r, :] = 0.0; A[r, r] = 1.0; B[r, :] = 0.0
    lam = linalg.eig(A, B, right=False)                             # generalised eigenvalues (some infinite, from the boundary rows)
    return lam[np.isfinite(lam)]                                    # drop the infinite ones

lam24 = benard_spectrum_mine(3.1163, Ra, Pr, 24)                    # our own spectrum with N = 24
mine = np.max(lam24[np.abs(lam24) < 1e4].real)                      # leading growth rate among the eigenvalues of reasonable size
assert np.isclose(mine, ch11.benard_growth_rate(3.1163, Ra, Pr, N=24), rtol=1e-8)   # same number -> the library does exactly this
print(f"from scratch (N = 24): sigma_1 = {mine:.8f}")
""", r"""
1. `cheb_mine` is the Chebyshev matrix in six lines; the assertion shows it equals `ST.cheb` and differentiates $x^3$ exactly.
2. `benard_spectrum_mine` assembles the two amplitude equations in block form for the unknown vector $[W;\hat T]$, replaces six
   rows by the six boundary conditions (P259) and calls `scipy.linalg.eig(A, B)`.
3. The leading eigenvalue agrees with `ch11.benard_growth_rate` to eight digits. (The library eliminates the boundary unknowns
   instead of replacing rows, which removes the spurious eigenvalues you will see in the next figure.)""")
P("P260", "real and imaginary parts of a complex identity", r"""
One complex equation is two real equations: its real parts agree and its imaginary parts agree. Three facts do most of the work in
this chapter: (i) $\int\lvert f\rvert^2dz>0$ unless $f\equiv0$ (a sum of non-negative numbers); (ii)
$\frac1{U-c}=\frac{U-c^*}{\lvert U-c\rvert^2}$, so $\mathrm{Im}\frac1{U-c}=\frac{c_i}{\lvert U-c\rvert^2}$ for real $U$ (multiply
top and bottom by the conjugate, Ch. 2 P81); (iii) a real quantity times $c_i$ is zero only if $c_i=0$ or the quantity vanishes.""",
  code=r"""
import numpy as np                                    # arrays
U, c = 0.3, 0.1 + 0.2j                                # a real velocity and a complex wave speed
print(np.imag(1/(U - c)), c.imag/abs(U - c)**2)       # both 2.5: Im{1/(U - c)} = c_i/|U - c|^2
""")
D("D08", ref="11.37", check_src=r"""
ex = ch11.exchange_of_stabilities_sympy()             # the two energy relations of D08 built with sympy on trial functions meeting (11.38)
print("relation of step 3 (T-equation times conj(T), integrated):", ex["relation_T"])   # 0: the relation holds identically
print("relation of step 8 (W-equation times -conj(W), integrated):", ex["relation_W"])  # 0
print("imaginary part of step 10:", ex["imag_part"])                                    # sigma_i times a bracket of positive integrals
print("its only solution for Ra > 0: sigma_i in", ex["sigma_i_solutions"])              # [0]: sigma is real
""")
note("N40 [B] · N122 [B]", r"""
**Exchange of stabilities.** D08 is Exercise 11.6 written out: with the positive integrals
$I_1=\int\lvert\hat T\rvert^2dz$, $I_2=\int(\lvert\hat T'\rvert^2+K^2\lvert\hat T\rvert^2)dz$,
$J_1=\int(\lvert W'\rvert^2+K^2\lvert W\rvert^2)dz$ and $J_2=\int\lvert(D^2-K^2)W\rvert^2dz$, the imaginary part of the combined
relation is $\sigma_i\big(\frac{J_1}\Pr+\mathrm{Ra}K^2I_1\big)=0$, so $\sigma$ is real when $\mathrm{Ra}>0$. A real $\sigma$ can
only change sign by passing *through* zero: the marginal state is $\sigma=0$, and what appears is a stationary pattern of cells.""")
fig(r"""
lam36 = benard_spectrum_mine(3.1163, Ra, Pr, 36)                    # the from-scratch spectrum again with more points
mode = ch11.benard_marginal_Ra(3.1163, return_mode=True)[1]         # the marginal eigenfunctions W(z), T(z) at K = 3.1163 (C04)
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.5))
a.plot(lam24.real, np.full(len(lam24), 24), "o", color=C_BUOY, ms=4, label="$N$ = 24")
a.plot(lam36.real, np.full(len(lam36), 36), "x", color=C_GROW, ms=5, label="$N$ = 36")
a.set_xscale("symlog", linthresh=50)                                # logarithmic for large |sigma|, linear near zero
a.set_xticks([-1e15, -1e10, -1e5, 0, 1e5, 1e10, 1e15])               # a few readable tick marks
a.axvline(0, color=C_BASE, lw=0.8)
a.legend(fontsize=8, loc="upper left")
a.set(xlabel=r"growth rate $\sigma$ [$\kappa/d^2$]", ylabel="number of points $N$", yticks=[24, 36], ylim=(18, 42))
a.set_title("(a) which eigenvalues to trust", fontsize=10)
b.plot(mode["W"]/np.max(np.abs(mode["W"])), mode["z"], color=C_BUOY, label="$W(z)$ (velocity)")
b.plot(mode["T"]/np.max(np.abs(mode["T"])), mode["z"], color=C_GROW, label=r"$\hat T(z)$ (temperature)")
b.set(xlabel="amplitude (scaled to 1)", ylabel="$z$ [$d$]")
b.set_title("(b) the leading eigenfunctions", fontsize=10)
b.legend(fontsize=8)
plt.show()
big24, big36 = lam24[np.abs(lam24) > 1e6], lam36[np.abs(lam36) > 1e6]
print(f"eigenvalues larger than 1e6 in size: {len(big24)} at N = 24, {len(big36)} at N = 36 (spurious: they move with N)")
print("leading three, N = 24:", np.round(np.sort(lam24[np.abs(lam24) < 1e4].real)[::-1][:3], 3), "| N = 36:", np.round(np.sort(lam36[np.abs(lam36) < 1e4].real)[::-1][:3], 3))
""",
    see="(a) Every finite eigenvalue of the from-scratch matrices for two resolutions, on an axis that is logarithmic for large "
        "values: near zero the circles and crosses sit on top of each other; far out to the left and right they do not. (b) The "
        "shapes of the leading mode: the velocity amplitude W is flat near the walls (no slip), the temperature amplitude is a "
        "single arch.",
    read="Trust what does not move with N. The top eigenvalue ≈ 19.9 > 0 grows; the next ones (≈ −41, −95, …) decay. The "
         "huge values (10⁶ and more, some positive) are artefacts of the boundary rows and round-off — a different set at every "
         "N. That is why every solver in this chapter keeps only eigenvalues that reappear at a larger N.",
    change="…Ra = 1000: the top eigenvalue drops to ≈ −7.8 — every disturbance decays.",
    explain=r"""
1. `benard_spectrum_mine` (the from-scratch cell) at $N=36$, next to the $N=24$ spectrum.
2. `ch11.benard_marginal_Ra(K, return_mode=True)` returns the marginal Rayleigh number and the profiles `W`, `T` on the nodes `z`.
3. `set_xscale("symlog")` gives an axis that is linear near zero and logarithmic beyond ±50, so that eigenvalues of size 20 and
   of size 10¹⁵ fit in one picture.""")
whatif(r"""
…we only want to know *when* convection starts, not how fast? Since $\sigma$ is real (D08), it passes through zero, not around
it: setting $\sigma=0$ removes Pr and leaves one equation for Ra(K) — C04.""")
core("C04", r"The rigid–rigid neutral curve Ra(K) and the critical point $\mathrm{Ra}_c=1707.76$, $K_c=3.117$ (Chandrasekhar 1961)",
     "At what heating does convection start, and how wide are the first cells?")
problem(r"""
Turn up the stove slowly. For a while nothing moves; then, at a sharp temperature difference, rolls appear — and they are all
about as wide as the layer is deep. Why that width, and why that threshold? Each possible cell width is a separate normal mode,
each with its own threshold; convection starts at the cheapest one.""")
idea(r"""
 Ra
  ▲  \                                /      narrow cells (large K): viscosity and diffusion act over a short distance
  │   \     marginal Ra(K)          /        wide cells (small K): buoyancy drives them only weakly
  │    \___                    ___/
  │        \___  ●  ______/        ← the bottom of the valley: Ra_c ≈ 1708 at K_c ≈ 3.12  (λ_c = 2π/K_c ≈ 2d)
  └────────────────────────────────► K
 heat until Ra reaches the valley bottom: the first K to go unstable is K_c
""")
remind("C04")
D("D09", ref="11.40")
P("P261", "even and odd functions; parity under z → −z", r"""
$f$ is *even* if $f(-z)=f(z)$ (cos, cosh, $z^2$) and *odd* if $f(-z)=-f(z)$ (sin, sinh, $z$). A problem that looks the same after
flipping $z\to-z$ (both walls alike) has eigenfunctions that are either even or odd, so we can look for them separately — even
$W$ gives one row of cells, odd $W$ two.""", code=r"""
import numpy as np                                    # arrays
z = np.linspace(-0.5, 0.5, 5)                         # points placed symmetrically about z = 0
print(np.allclose(np.cosh(2*z), np.cosh(-2*z)), np.allclose(np.sinh(2*z), -np.sinh(-2*z)))   # True True: cosh even, sinh odd
""")
note("N44 [B]", r"""
**Pr drops out; even and odd modes.** At $\sigma=0$ the Prandtl number has disappeared from the marginal problem, so the onset does
not depend on the fluid's Pr — only Ra matters. Both walls are alike, so the eigenfunctions are even ($W$ symmetric about $z=0$:
one row of cells) or odd (two rows); the even one goes unstable first.""")
fig(r"""
ev = ch11.benard_eigenfunction(3.1163, mode="even")                 # marginal even mode at K = 3.1163: fields on an (x, z) grid
od = ch11.benard_eigenfunction(5.3647, mode="odd")                  # marginal odd mode at its own critical K = 5.3647
fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.0))
for ax, e, ttl in ((axs[0], ev, "even mode: one row of rolls"), (axs[1], od, "odd mode: two rows")):
    ax.pcolormesh(e["x"], e["z"], e["T"], cmap="RdBu_r", shading="auto")   # temperature disturbance: red warm, blue cold
    ax.contour(e["x"], e["z"], e["psi"], 8, colors=C_BASE, linewidths=0.8) # streamlines of the rolls
    ax.set(xlabel="$x$ [$d$]", ylabel="$z$ [$d$]")
    ax.set_title(f"{ttl} (Ra = {e['Ra']:.0f})", fontsize=10)
plt.show()
""",
    see="Left: one row of counter-rotating rolls (grey streamlines) with warm fluid (red) rising between them and cold fluid (blue) "
        "sinking. Right: two rows of smaller rolls stacked on top of each other.",
    read="Warm columns sit where the streamlines point upward: that is the feedback of D05 — rising fluid is warmer, and warmer "
         "fluid rises. The odd mode needs smaller cells, so viscosity and diffusion cost it more: its marginal Ra is ten times "
         "larger (N123).",
    change="…free walls (C05): the even mode's vertical shape becomes exactly a cosine.",
    explain=r"""
1. `ch11.benard_eigenfunction(K, mode=…)` solves the marginal problem at wavenumber $K$ and returns the roll fields on a grid:
   `x`, `z`, the stream function `psi`, the temperature disturbance `T`, and the marginal `Ra`.
2. `pcolormesh` paints the temperature, `contour` draws streamlines (Ch. 2, P78).""")
P("P262", "cube roots of a negative number", r"""
Every nonzero number has three cube roots. Write $-1=e^{\mathrm i\pi}$: its cube roots are $e^{\mathrm i\pi/3}$,
$e^{\mathrm i\pi}=-1$ and $e^{\mathrm i5\pi/3}$, i.e. $-1$ and $\tfrac12(1\pm\mathrm i\sqrt3)$. So $s^3=-a$ ($a>0$) has
$s=-a^{1/3}$ and $s=\tfrac12a^{1/3}(1\pm\mathrm i\sqrt3)$ — one real, two complex conjugates (Euler's formula, Ch. 1 P45).""",
  code=r"""
import numpy as np                                    # arrays
print(np.roots([1, 0, 0, 1]))                         # s^3 + 1 = 0  ->  -1 and 0.5 + 0.866j, 0.5 - 0.866j
print((0.5*(1 + 1j*np.sqrt(3)))**3)                   # (-1+0j): the complex root cubed is -1 again
""")
P("P263", "neutral curve as the zero contour of the growth rate; critical point as its minimum", r"""
For each wavenumber $K$ find the parameter value where the leading growth rate crosses zero (`brentq`, Ch. 3 P108: a root finder
that needs two values with opposite signs): that is the neutral (marginal) curve Ra(K). Then minimise Ra(K) over $K$
(`minimize_scalar`, Ch. 7 P170): the minimum is the critical point (Ra_c, K_c). Inside the curve $\sigma>0$, outside $\sigma<0$.
The same two nested searches give the critical Taylor number (C07) and the critical Reynolds number (C13).""", code=r"""
import numpy as np                                    # arrays
from scipy.optimize import minimize_scalar            # minimum of a function of one variable
Ra_of_K = lambda K: (np.pi**2 + K**2)**3/K**2         # a neutral curve given in closed form (C05 derives it)
r = minimize_scalar(Ra_of_K, bounds=(0.5, 6), method="bounded")   # search K between 0.5 and 6
print(round(r.x, 3), round(r.fun, 1))                 # 2.221 657.5: the critical K and the critical Ra
""")
D("D10", ref="11.42", after=r"""
> ⚠️ **slip #1 — the book prints** the third term of $\big(\frac{d^2}{dz^2}-K^2\big)^2W$ as $B(q^{*2}-K^2)^2\cosh q^*z$ **; the
> correct form is** $C(q^{*2}-K^2)^2\cosh q^*z$ — the coefficient of the third solution is $C$ (step 9). The 3 × 3 matrix printed
> under it is right.""")
nb.worked_example("the characteristic roots at K = 2, Ra = 1000", r"""
1. $\mathrm{Ra}/K^4=1000/16=62.5$; $s=(62.5)^{1/3}=3.9685$.
2. First root of the characteristic equation: $q^2=-K^2(s-1)=-4\times2.9685=-11.874$, so $q=\pm\mathrm iq_0$ with
   $q_0=\sqrt{11.874}=3.4459$.
3. The other two: $q^2=K^2[1+\tfrac12s(1\pm\mathrm i\sqrt3)]=4[1+1.9843\pm3.4369\,\mathrm i]=11.937\pm13.748\,\mathrm i$, so
   $q=3.8822+1.7705\,\mathrm i$ and its conjugate.
4. Check: $(q^2-K^2)^3$ for the first root $=(-11.874-4)^3=(-15.874)^3=-4000=-\mathrm{Ra}K^2$ ✓.
5. Below $\mathrm{Ra}=K^4=16$ the cube root $s<1$ and $q_0$ turns imaginary — why the determinant search starts above $K^4$.""")
code(r"""
print("q0, q, q* at Ra = 1000, K = 2:", np.round(ch11.benard_char_roots(1000.0, 2.0), 4))   # the roots (11.42) of the worked example
for K in (2.0, 3.1163, 5.0):                                        # two independent routes to the marginal Ra at three wavenumbers
    print(f"K = {K}: Chebyshev eigenproblem {ch11.benard_marginal_Ra(K):.4f} | book's determinant {ch11.benard_marginal_Ra_det(K):.4f}")
crit = ch11.benard_critical()                                       # minimum of Ra(K) over K (P263): rigid-rigid walls
ref = BENCH["benard_rigid_rigid"]                                   # published value (Chandrasekhar 1961), from reference/ch11
print(f"critical point: Ra_c = {crit['Ra_c']:.3f} at K_c = {crit['K_c']:.4f}  (published {ref['Ra_c']}, {ref['K_c']})")
print(f"relative difference from the published Ra_c: {abs(crit['Ra_c']/ref['Ra_c'] - 1):.1e};  cell pair width 2 pi/K_c = {2*np.pi/crit['K_c']:.3f} d")
""", explain=r"""
1. `ch11.benard_char_roots(Ra, K)` returns $q_0$, $q$, $q^*$ of D10 steps 4–6.
2. `ch11.benard_marginal_Ra(K)` solves the marginal pair as a generalised eigenproblem *for Ra* on Chebyshev points;
   `ch11.benard_marginal_Ra_det(K)` finds the zero of the book's 3 × 3 determinant (D10 steps 11–13). Two routes, one answer.
3. `ch11.benard_critical()` is P263's recipe: minimise over $K$. The result agrees with the published benchmark to better than one
   part in 10⁵; the book rounds it to 1708 and 3.12. Our converged wavenumber is $K_c=3.1163$; the 1961 table (and this block's
   heading) rounds it to 3.117.
4. One wavelength $2\pi/K_c\approx2.02\,d$ holds a *pair* of counter-rotating rolls, so each roll is about as wide as the layer is
   deep.""")
scratch(r"""
def det_mine(Ra, K):                                                # the 3 x 3 determinant of D10 step 11 (even mode, rows at z = +1/2)
    s = (Ra/K**4)**(1/3)                                            # the cube-root factor of (11.42)
    q0 = K*np.sqrt(s - 1)                                           # first root: q = +- i q0 (real q0 for Ra > K^4)
    q = np.sqrt(K**2*(1 + 0.5*s*(1 + 1j*np.sqrt(3))))               # second root (complex; np.sqrt of a complex number = principal root)
    qs = np.conj(q)                                                 # third root: its conjugate
    M = np.array([[np.cos(q0/2), np.cosh(q/2), np.cosh(qs/2)],                                   # W = 0
                  [-q0*np.sin(q0/2), q*np.sinh(q/2), qs*np.sinh(qs/2)],                          # dW/dz = 0
                  [(q0**2 + K**2)**2*np.cos(q0/2), (q**2 - K**2)**2*np.cosh(q/2), (qs**2 - K**2)**2*np.cosh(qs/2)]])   # (D^2 - K^2)^2 W = 0
    return np.linalg.det(M)                                         # zero determinant <=> a nonzero (A, B, C) exists (Ch. 1 P53, P56)

K = 3.1163                                                          # close to the critical wavenumber
Ra_scan = np.linspace(1.0001*K**4 + 1, 4000.0, 400)                 # start just above K^4, where q0 becomes real
im = np.array([det_mine(R, K).imag for R in Ra_scan])               # the determinant is purely imaginary (D10 step 12)
i0 = np.nonzero(np.diff(np.sign(im)))[0][0]                         # first sign change of its imaginary part (Ch. 7 P180)
Ra_det = brentq(lambda R: det_mine(R, K).imag, Ra_scan[i0], Ra_scan[i0 + 1], xtol=1e-10)   # refine the root (P108)
assert np.isclose(Ra_det, ch11.benard_marginal_Ra(K), rtol=1e-8)    # the determinant route = the Chebyshev route
d_ok, d_bad = ch11.benard_determinant(1707.762, K), ch11.benard_determinant(1707.762, K, printed=True)   # correct, and with slip #1 built in
print(f"my determinant root: Ra = {Ra_det:.4f}")
print(f"real part / size of det: correct {abs(d_ok.real)/abs(d_ok):.1e}, as printed (slip #1) {abs(d_bad.real)/abs(d_bad):.2f}")
""", r"""
1. `det_mine` types the matrix of D10 step 11 with `np.cos` and `np.cosh` of complex arguments (they accept complex numbers:
   $\cosh(a+\mathrm ib)=\cosh a\cos b+\mathrm i\sinh a\sin b$).
2. A scan above $\mathrm{Ra}=K^4$ finds the first sign change of the imaginary part; `brentq` refines it. The assertion shows the
   root equals the Chebyshev eigenvalue to eight digits.
3. The library's `printed=True` variant builds slip #1 in: the determinant is then no longer purely imaginary (its real part is
   most of it), so the printed line cannot be the one that was actually solved.""")
fig(r"""
t_ = np.genfromtxt(ROOT / "reference" / "ch11" / "benard_neutral_curves.csv", delimiter=",", skip_header=1, names=True)   # our stored table
tab = {"K": t_["K"], "rigid": t_["Ra_rigid_rigid"], "free": t_["Ra_free_free"], "rigid_free": t_["Ra_rigid_free"], "odd": t_["Ra_odd_rigid_rigid"]}
assert np.isclose(tab["rigid"][15], ch11.benard_marginal_Ra(tab["K"][15]), rtol=6e-4)   # one table entry recomputed now (the table is rounded)
cr = {"rigid": crit, "free": ch11.benard_critical(bc=("free", "free")),
      "rigid_free": ch11.benard_critical(bc=("rigid", "free"), mode="any"), "odd": ch11.benard_critical(mode="odd")}   # the four minima
K1 = brentq(lambda K: ch11.benard_marginal_Ra(K) - 2500.0, 1.0, crit["K_c"])   # lower edge of the unstable band at Ra = 2500 (rigid)
K2 = brentq(lambda K: ch11.benard_marginal_Ra(K) - 2500.0, crit["K_c"], 7.0)   # upper edge
fig, ax = plt.subplots(figsize=(7.2, 4.0))
ax.fill_betweenx(tab["K"], tab["rigid"], 3e4, color=C_GROW, alpha=0.12)
ax.semilogx(tab["rigid"], tab["K"], color=C_NEUT, lw=2.5, label="rigid–rigid")
ax.semilogx(tab["rigid_free"], tab["K"], "--", color=C_NEUT, label="rigid–free")
ax.semilogx(tab["free"], tab["K"], color=C_BASE, label="free–free")
ax.semilogx(tab["odd"], tab["K"], ":", color=C_NEUT, label="odd mode (rigid–rigid)")
for name in ("rigid", "rigid_free", "free"):
    ax.plot(cr[name]["Ra_c"], cr[name]["K_c"], "D", color=COLORS["ink"], ms=5)
    ax.text(cr[name]["Ra_c"]*0.97, cr[name]["K_c"] - 0.55, f"{cr[name]['Ra_c']:.2f}\n$K$ = {cr[name]['K_c']:.3f}", fontsize=7, ha="right")
ax.plot([2500, 2500], [K1, K2], color=C_GROW, lw=3)
ax.text(8000, 3.0, r"$\sigma>0$", color=C_GROW)
ax.text(450, 6.5, r"$\sigma<0$", color=C_DECAY)
ax.set(xlabel="Rayleigh number Ra [–]", ylabel="wavenumber $K$ [$1/d$]", xlim=(400, 3e4), ylim=(0.5, 8))
ax.set_title("Onset is the bottom of a valley", fontsize=11)
ax.legend(fontsize=8, loc="upper right")
savefig(fig, "ch11", "c04_neutral_curves")
plt.show()
print(f"at Ra = 2500 (rigid walls) the wavenumbers from K = {K1:.2f} to {K2:.2f} grow")
print({n: (round(v["Ra_c"], 2), round(v["K_c"], 4)) for n, v in cr.items()})
""",
    see="Neutral curves in the (Ra, K) plane, drawn as the book does (K vertical, Ra horizontal on a logarithmic axis): rigid–rigid "
        "(bold purple), rigid–free (dashed), free–free (grey), and the odd mode (dotted, far to the right). Diamonds mark the "
        "leftmost point of each curve; the rose region is where the rigid–rigid layer is unstable.",
    read="For a given Ra, the K inside the curve grow: the thick rose bar marks the unstable band at Ra = 2500 for rigid walls. "
         "The leftmost point of a curve is the critical point — below that Ra nothing grows. Rigid walls hold the fluid most, so "
         "their valley starts last, at $\\mathrm{Ra}_c=1707.76$; free walls first, at 657.51.",
    change="…we double the depth d: Ra grows 8× (it scales with d³), so a layer twice as deep convects at one eighth of the "
           "temperature difference.",
    explain=r"""
1. `np.genfromtxt` reads our small public table `reference/ch11/benard_neutral_curves.csv` (the four neutral curves on one grid
   of 96 wavenumbers, written by `ch11.benard_neutral_table`); one entry is recomputed live and compared (the stored numbers are
   rounded, hence the relative tolerance).
2. The four minima are computed live with `ch11.benard_critical` (P263); `mode="any"` is needed for rigid–free walls, which have no
   even/odd symmetry.
3. `brentq` finds the two wavenumbers where the rigid curve passes through Ra = 2500.""")
nb.md(r"""
**Slider figure F2 — the unstable band opens at the valley bottom.** Drag Ra; the thick segments show, for each pair of walls,
which wavenumbers grow.""")
nb.plotly(r"""
def band(Ra_curve, Ra_now):                                        # the K-interval where the neutral curve lies below Ra_now
    inside = tab["K"][Ra_curve < Ra_now]                           # wavenumbers of the table that are unstable at this Ra
    return (np.array([inside[0], inside[-1]]), np.array([Ra_now, Ra_now])) if len(inside) else (np.array([np.nan]), np.array([np.nan]))

def f2(Ra_now):                                                    # all traces for one slider value
    out = {"rigid–rigid": (tab["K"], tab["rigid"]), "rigid–free": (tab["K"], tab["rigid_free"]), "free–free": (tab["K"], tab["free"])}
    for name, key in (("unstable K (rigid–rigid)", "rigid"), ("unstable K (rigid–free)", "rigid_free"), ("unstable K (free–free)", "free")):
        out[name] = band(tab[key], Ra_now)                         # a horizontal segment at height Ra_now
    return out

figF2 = slider_figure(f2, "Ra", np.round(np.geomspace(500.0, 2.0e4, 15 if FAST else 30)), xlabel="wavenumber K [1/d]",
                      ylabel="Rayleigh number Ra", title="Which cell widths grow at this Rayleigh number?",
                      xrange=(0.5, 8), yrange=(400, 3e4))
figF2.update_yaxes(type="log", range=[np.log10(400), np.log10(3e4)])   # logarithmic Ra axis
figF2.show()
""", explain=r"""
1. `band` reads the unstable wavenumbers off the table for one Ra; `f2` returns the three static curves and the three bands.
2. Below 657.5 there is no band at all; the free–free band opens first, the rigid–rigid one last, each at its valley bottom.""")
nb.live(r"""
def benard_live(K=3.1, logRa=3.4, walls="rigid–rigid", Pr=7.0):    # free K, Ra (as log10), boundaries and Prandtl number
    bc = {"rigid–rigid": ("rigid", "rigid"), "rigid–free": ("rigid", "free"), "free–free": ("free", "free")}[walls]
    Ra_now = 10.0**logRa                                           # the Rayleigh number from the slider
    try:                                                           # the solver raises ValueError if the leading eigenvalue is not confirmed at a larger N
        sigma = ch11.benard_growth_rate(K, Ra_now, Pr, bc=bc)      # leading growth rate [kappa/d^2]
    except ValueError:
        sigma = np.nan                                             # show 'nan' for that slider position instead of stopping
    Ra_m, m = ch11.benard_marginal_Ra(K, bc=bc, mode="any", return_mode=True)   # the marginal Ra at this K and its mode shape
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    ax.plot(m["W"]/np.max(np.abs(m["W"])), m["z"], color=C_GROW if sigma > 0 else C_DECAY)
    ax.set(xlabel="$W(z)$ (scaled)", ylabel="$z$ [$d$]")
    ax.set_title(f"Ra = {Ra_now:.0f}, marginal {Ra_m:.0f}: sigma = {sigma:.2f}", fontsize=9)
    plt.show()

live(benard_live, K=(1.0, 7.0, 0.1), logRa=(2.6, 4.4, 0.05), walls=["rigid–rigid", "rigid–free", "free–free"], Pr=(0.1, 10.0, 0.1))
""", explain=r"""
1. `benard_live` computes the growth rate at your $(K,\mathrm{Ra},\Pr)$ and the marginal Ra at that $K$, and draws the mode shape
   $W(z)$ — rose if it grows, teal if it decays.
2. `live(...)` (Ch. 1, P47) turns the arguments into sliders and a drop-down. Changing Pr changes $\sigma$ but never its sign (N44).
   The slider figure F2 above carries the same idea on the web page.""")
nb.md(r"""
**Animation A3 — the same rolls just above and just below onset.** Rolls of the critical width at $\mathrm{Ra}=1.3\,\mathrm{Ra}_c$
(left) and $0.7\,\mathrm{Ra}_c$ (right) for water ($\Pr=7.14$); the amplitude follows $e^{\sigma t}$ with $\sigma$ from
`ch11.benard_growth_rate`; time is in units of $d^2/\kappa$.""")
nb.animation(r"""
nf = 14 if FAST else 18                                            # number of frames
Kc, Rac = crit["K_c"], crit["Ra_c"]                                # the critical point (rigid walls)
s_up = ch11.benard_growth_rate(Kc, 1.3*Rac, Pr)                    # growth rate 30 % above onset [kappa/d^2]
s_dn = ch11.benard_growth_rate(Kc, 0.7*Rac, Pr)                    # ... and 30 % below
e = ch11.benard_eigenfunction(Kc)                                  # the roll shape at K_c (fields on an (x, z) grid)
tA = np.linspace(0.0, 0.5, nf)                                     # times [d^2/kappa]
fig, axs = plt.subplots(1, 2, figsize=(8.4, 2.9))
meshes = []                                                        # the two temperature images
for ax, a0 in zip(axs, (np.exp(-s_up*tA[-1]), 1.0)):               # start small on the left (so that it ends at 1), at 1 on the right
    meshes.append(ax.pcolormesh(e["x"], e["z"], a0*e["T"]/np.max(np.abs(e["T"])), cmap="RdBu_r", vmin=-1, vmax=1, shading="auto"))
    ax.contour(e["x"], e["z"], e["psi"], 6, colors=C_BASE, linewidths=0.7)
    ax.set(xlabel="$x$ [$d$]", ylabel="$z$ [$d$]")
titles = [axs[0].set_title("", fontsize=9, color=C_GROW), axs[1].set_title("", fontsize=9, color=C_DECAY)]
Tn = e["T"]/np.max(np.abs(e["T"]))                                 # temperature pattern scaled to 1

def update(i):                                                     # frame i: time tA[i]
    amp_up = np.exp(s_up*(tA[i] - tA[-1]))                         # growing amplitude (reaches 1 at the last frame)
    amp_dn = np.exp(s_dn*tA[i])                                    # decaying amplitude (starts at 1)
    meshes[0].set_array((amp_up*Tn).ravel())
    meshes[1].set_array((amp_dn*Tn).ravel())
    titles[0].set_text(f"Ra = 1.3 Ra_c: grows, sigma = {s_up:+.2f}, t = {tA[i]:.2f}")
    titles[1].set_text(f"Ra = 0.7 Ra_c: decays, sigma = {s_dn:+.2f}, t = {tA[i]:.2f}")
    return meshes + titles

show_animation(animate(update, frames=nf, fig=fig, interval=130), player="video", dpi=60)
print(f"growth rates at K_c: {s_up:+.3f} (1.3 Ra_c), {s_dn:+.3f} (0.7 Ra_c); with Pr = 0.7 (air): "
      f"{ch11.benard_growth_rate(Kc, 1.3*Rac, 0.7):+.3f}, {ch11.benard_growth_rate(Kc, 0.7*Rac, 0.7):+.3f}")
""", explain=r"""
1. Two growth rates at the same wavenumber $K_c$: positive 30 % above the critical Rayleigh number, negative 30 % below.
2. The roll pattern (temperature in colour, streamlines in grey) is the marginal mode; only its amplitude $e^{\sigma t}$ changes.
3. The last line repeats the growth rates for air ($\Pr=0.7$): different numbers, same signs.""")
nb.figure_notes(
    see="The same roll pattern brightening on the left and fading on the right.",
    read="Onset is a sign change of the growth rate at fixed K: nothing about the *shape* of the rolls changes as Ra passes the "
         "critical value, only whether they grow or die.",
    change="…Pr = 0.7 (air): the growth rates change (the printed line), the threshold does not (N44).")
note("N50 [B]", r"""
**Rigid bottom, free top** (a liquid layer open to air): our computation gives $\mathrm{Ra}_c=1100.65$ at $K_c=2.682$
(Chandrasekhar 1961 tabulates this case); one free wall lowers the threshold because a stress-free surface resists less.""")
code(r"""
rf = ch11.benard_critical(bc=("rigid", "free"), mode="any")         # rigid bottom, stress-free top
print(f"rigid-free: Ra_c = {rf['Ra_c']:.3f} at K_c = {rf['K_c']:.4f}  (published {BENCH['benard_rigid_free']['Ra_c']}, {BENCH['benard_rigid_free']['K_c']})")
""")
note("N123 [B]", r"""
**The gravest odd mode** (Exercise 11.7): two rows of cells need $\mathrm{Ra}=17\,610.39$ at $K=5.365$ (Chandrasekhar 1961) — ten
times the even threshold, so it is never the first to appear.""")
code(r"""
odd = ch11.benard_critical(mode="odd")                              # the odd (two-row) mode between rigid walls
print(f"odd mode: Ra_c = {odd['Ra_c']:.2f} at K_c = {odd['K_c']:.4f} - {odd['Ra_c']/crit['Ra_c']:.1f} times the even threshold")
""")
P("P264", "Helmholtz equation in the plane; planforms as sums of cosines", r"""
Linear theory fixes only $\lvert\mathbf K\rvert$. Any horizontal pattern $f(x,y)$ with $\nabla_H^2f=-K^2f$ (the Helmholtz
equation) has the same growth rate: rolls $f=\cos Kx$, squares $\cos Kx+\cos Ky$, hexagons
$\sum_{j=1}^3\cos(\mathbf k_j\cdot\mathbf x)$ with three wave vectors of length $K$ at 120° to each other.""", code=r"""
import numpy as np                                    # arrays
from fluidpy import ch11_instability as ch11          # the chapter module
x = y = np.linspace(0, 2, 41)                         # a 2 x 2 patch of the layer (in depths d)
X, Y = np.meshgrid(x, y)                              # grid of points (Ch. 2 P76)
f = ch11.planform(X, Y, 3.0, "hexagons")              # a hexagonal pattern with |K| = 3
print(f.shape, round(f.max(), 3))                     # (41, 41) and 3.0: the maximum sits at the cell centres
""")
note("N51 [B]", r"""
**Which pattern?** Linear theory cannot choose between rolls, squares and hexagons (P264); experiments near onset often show
hexagons first, rolls as Ra rises, then time-dependent and finally turbulent convection at much larger Ra. In liquids the cell
centres usually rise (viscosity falls with temperature), in gases they sink. Climate hook: cloud streets are rolls aligned with the
wind; open and closed hexagonal cells cover the subtropical oceans in satellite pictures.""")
fig(r"""
xg = yg = np.linspace(0, 6, 121)                                   # a 6 d by 6 d patch seen from above
Xg, Yg = np.meshgrid(xg, yg)                                       # grid of points
fig, axs = plt.subplots(1, 3, figsize=(9.0, 3.0))
for ax, kind in zip(axs, ("rolls", "squares", "hexagons")):
    ax.pcolormesh(Xg, Yg, ch11.planform(Xg, Yg, crit["K_c"], kind), cmap="RdBu_r", shading="auto")   # vertical velocity: red up, blue down
    ax.set(xlabel="$x$ [$d$]", ylabel="$y$ [$d$]")
    ax.set_title(kind, fontsize=10)
    ax.set_aspect("equal")
plt.show()
""",
    see="Three patterns seen from above (red = rising, blue = sinking): parallel rolls, a chequerboard of squares, a honeycomb of "
        "hexagons.",
    read="All three are built from waves of the same wavenumber K_c, so linear theory gives all three the same growth rate.",
    change="…we go beyond linear theory: the nonlinear terms decide which pattern survives (named, N52).")
note("N52 [C]", r"""
**After onset** (named): the cells grow until the nonlinear terms balance them — kinetic energy produced by buoyancy (warm fluid
rising) equals viscous dissipation, fed by the heating. The energy bookkeeping is C14's; the simplest nonlinear model of it is
Lorenz's (C15).""")
explainer("benard_neutral_curve", "Why does convection start near Ra ≈ 1708?", r"""
Dragging a point in the (K, Ra) plane shows that every point is a whole eigenproblem — inside the curve the rolls grow, outside
they decay, and a third view shows the determinant's imaginary part crossing zero as Ra passes the curve; switching boundaries
moves the whole valley: from 657.5 through 1100.65 to 1707.76.""", [
    "Preset 'rigid onset': read Ra_c = 1707.76 at K = 3.117; nudge Ra up by 5 % and watch the rolls start growing.",
    "Keep Ra = 2500 and drag K from 1 to 6: the status flips twice — the unstable band.",
    "Switch the walls to free–free: the valley drops to 657.5 and moves to K = π/√2 ≈ 2.22.",
])
whatif(r"""
…both walls were stress-free? Every even derivative of $W$ vanishes at the walls, sine modes solve the problem exactly, and the
neutral curve becomes a formula you can minimise by hand — C05.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C05", r"Free–free boundaries: $\mathrm{Ra}=(n^2\pi^2+K^2)^3/K^2$ (11.44), $K_c^2=\pi^2/2$, $\mathrm{Ra}_c=\tfrac{27}{4}\pi^4$",
     "Can we find the onset of convection with pencil and paper?")
problem(r"""
Rigid plates are easy to build but hard to solve. A layer whose boundaries cannot hold any shear stress — oil floating on mercury,
or a slice of the atmosphere with no lid — is the case Rayleigh solved in 1916. It gives an exact formula, and its minimum,
$27\pi^4/4\approx657.5$, is the number double diffusion (C06) and Lorenz's model (C15) are built on.""")
idea("", r"""
With $W=W''=W''''=0$ at the walls, a sine wave $\sin n\pi(z+\tfrac12)$ vanishes at both walls together with all its even
derivatives, and $\big(\frac{d^2}{dz^2}-K^2\big)$ acting on it just multiplies it by $-(n^2\pi^2+K^2)$. Every operator becomes a
number, and the sixth-order marginal equation becomes algebra.""")
remind("C05")
P("P265", "quotient rule, and differentiating with respect to K² as the variable", r"""
$\frac{d}{dx}\frac{f}{g}=\frac{f'g-fg'}{g^2}$. When a formula depends on $K$ only through $K^2$, call $x=K^2$ and differentiate in
$x$: the minimum over $K$ and over $x$ is the same point (for $K>0$).""", code=r"""
import sympy as sp                                    # symbolic algebra
x = sp.symbols("x", positive=True)                    # x stands for K^2
Ra_x = (sp.pi**2 + x)**3/x                            # the free-free neutral curve as a function of x = K^2
print(sp.solve(sp.diff(Ra_x, x), x))                  # [pi**2/2]: the minimum is at K^2 = pi^2/2
""")
D("D11", ref="11.44", check_src=r"""
ff = ch11.benard_free_free_sympy()                    # D11 redone with sympy
print("even derivatives of sin n pi (z + 1/2) at the two walls (step 7):", ff["bc_W4"])    # all 0
print("residual of the sixth-order equation with Ra = (n^2 pi^2 + K^2)^3/K^2 (step 8):", ff["residual"])   # 0
print("dRa/dK^2 (step 9):", ff["dRa_dK2"], "| root (step 10):", ff["root"], "| Ra_c (step 11):", ff["Ra_c"])
print("slip #3, as printed:", ff["printed_dRa_dK2"], "-> roots:", ff["printed_root"])       # []: the printed derivative has no root
print("slip #2: sin(n pi z) with n = 1 at z = -1/2, +1/2:", ff["printed_bc_n1"])           # [1, -1]: it does NOT vanish at the walls
""", after=r"""
> ⚠️ **slip #2 — the book prints** $W=A\sin(n\pi z)$ **; the correct form is** $W=A\sin n\pi(z+\tfrac12)$: on
> $z\in[-\tfrac12,\tfrac12]$ the printed sine vanishes at both walls only for even $n$; the family that meets the conditions is
> $\sin n\pi(z+\tfrac12)$, whose $n=1$ member is $\cos\pi z$ (step 7).

> ⚠️ **slip #3 — the book prints** $\frac{d\mathrm{Ra}}{dK^2}=\frac{3(\pi^2+K^2)^2}{K^2}-\frac{3(\pi^2+K^2)^3}{K^4}$ **; the correct
> form is** $\frac{d\mathrm{Ra}}{dK^2}=\frac{3(\pi^2+K^2)^2}{K^2}-\frac{(\pi^2+K^2)^3}{K^4}$ — no 3 on the second term (step 9).
> With the printed 3 the root equation reads $K^2=\pi^2+K^2$, which has no solution.""")
D("D12", ref="11.44")
nb.worked_example("the free–free neutral curve with easy K", r"""
1. $K=2$: $\mathrm{Ra}=(\pi^2+4)^3/4=(13.870)^3/4=667.0$.
2. $K=1$: $(10.870)^3/1=1284.2$.
3. $K=3$: $(18.870)^3/9=746.5$.
4. The minimum is between 1 and 3: $d\mathrm{Ra}/dK^2=0$ at $K^2=\pi^2/2=4.935$, $K=2.2214$.
5. There $\mathrm{Ra}=(3\pi^2/2)^3/(\pi^2/2)=(27\pi^6/8)\cdot(2/\pi^2)=27\pi^4/4=657.51$.
6. Growth rate at $K=\pi/\sqrt2$, $\mathrm{Ra}=2000$, $\Pr=1$ (D12): $a^2=3\pi^2/2=14.804$;
   $(\sigma+a^2)^2=\mathrm{Ra}K^2/a^2=2000\times4.935/14.804=666.7$; $\sigma+a^2=25.82$; $\sigma=11.02$ (in units $\kappa/d^2$).""")
code(r"""
print("Ra(K) for K = 1, 2, 3, 4:", [round(ch11.benard_free_free_Ra(K), 2) for K in (1.0, 2.0, 3.0, 4.0)])   # (11.44) with n = 1
print(ch11.benard_free_free_critical())                             # the minimum: Ra_c = 27 pi^4/4, K_c = pi/sqrt(2), K_c^2 = pi^2/2
Kff = np.pi/np.sqrt(2)                                              # the critical wavenumber 2.2214
print("sigma_+, sigma_- at Ra = 2000, Pr = 1 (D12):", np.round(ch11.benard_free_free_sigma(Kff, 2000.0, 1.0), 4))   # D12: the growing and the decaying root [kappa/d^2]
print("same from the Chebyshev eigen-solver with free walls:", round(ch11.benard_growth_rate(Kff, 2000.0, 1.0, bc=("free", "free")), 4))   # independent check of the growing root
print("critical point from the eigen-solver:", ch11.benard_critical(bc=("free", "free")))   # should equal 27 pi^4/4 at pi/sqrt(2)
print("with Pr = 7:", np.round(ch11.benard_free_free_sigma(Kff, 2000.0, 7.0), 3))   # water instead of Pr = 1: faster growth, same sign
""", explain=r"""
1. `ch11.benard_free_free_Ra(K)` is (11.44): 1284.23, 667.01, 746.53, 1082.06 — the worked example.
2. `ch11.benard_free_free_critical()` returns the minimum in closed form: 657.511 at 2.2214.
3. D12's two growth rates, 11.0155 and −40.6243, and the same leading value from the general Chebyshev solver with free walls —
   the closed form and the eigen-solver agree.
4. With $\Pr=7$ the positive root is larger (22.26): Pr changes how fast, not whether.""")
scratch(r"""
r = minimize_scalar(lambda K: (np.pi**2 + K**2)**3/K**2, bounds=(0.5, 6), method="bounded")   # (11.44) typed out and minimised over K
assert np.isclose(r.fun, 27*np.pi**4/4) and np.isclose(r.x, np.pi/np.sqrt(2), rtol=1e-5)      # 657.51 at K = pi/sqrt(2)
a2 = np.pi**2 + Kff**2                                              # a^2 = pi^2 + K^2 at K_c
s_mine = 0.5*(-2*a2 + np.sqrt(4*2000.0*Kff**2/a2))                  # D12's last step for Pr = 1: sigma_+ = -a^2 + sqrt(Ra K^2/a^2)
assert np.isclose(s_mine, ch11.benard_free_free_sigma(Kff, 2000.0, 1.0)[0])                   # same number as the library
assert not np.allclose(ch11.benard_free_free_mode(np.array([-0.5, 0.5]), 1, printed=True), 0) # slip #2: sin(pi z) is not zero at the walls
assert np.allclose(ch11.benard_free_free_mode(np.array([-0.5, 0.5]), 1), 0)                   # the corrected mode is
print(f"minimum by minimize_scalar: Ra = {r.fun:.3f} at K = {r.x:.4f};  sigma_+ by hand = {s_mine:.4f}")
""", r"""
1. `minimize_scalar` (P263) on (11.44): the same minimum as the closed form.
2. D12's growth rate for $\Pr=1$ typed in one line.
3. The printed mode shape of slip #2 does not vanish at the walls; the corrected one does.""")
fig(r"""
Kp = np.linspace(0.6, 6.5, 300)                                    # wavenumbers [1/d]
fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, 3.6))
a.plot(Kp, ch11.benard_free_free_Ra(Kp), color=C_NEUT, label="free–free, $n$ = 1 (11.44)")
a.plot(Kp, ch11.benard_free_free_Ra(Kp, n=2), "--", color=C_NEUT, label="$n$ = 2")
a.plot(tab["K"], tab["rigid"], ":", color=C_BASE, label="rigid–rigid (C04)")
a.plot(Kff, 27*np.pi**4/4, "D", color=COLORS["ink"])
a.text(Kff + 0.2, 520, "657.51 at $K$ = 2.221", fontsize=8)
a.set_yscale("log")                                                # logarithmic Ra axis: the n = 2 curve (minimum 108 pi^4 = 10520) fits too
a.set(xlabel="wavenumber $K$ [$1/d$]", ylabel="marginal Ra [–]", ylim=(400, 1e5), xlim=(0.5, 6.5))
a.set_title("(a) a valley you can draw by hand", fontsize=10)
a.legend(fontsize=8)
for Ra_now, col in ((500.0, C_DECAY), (27*np.pi**4/4, C_NEUT), (2000.0, C_GROW)):
    b.plot(Kp, [ch11.benard_free_free_sigma(K, Ra_now, 1.0)[0] for K in Kp], color=col, label=f"Ra = {Ra_now:.1f}")   # sigma_+(K), Pr = 1
b.axhline(0, color=C_BASE, lw=0.8)
b.set(xlabel="wavenumber $K$ [$1/d$]", ylabel=r"growth rate $\sigma_+$ [$\kappa/d^2$]", ylim=(-60, 15))
b.set_title("(b) the growth curve lifts through zero", fontsize=10)
b.legend(fontsize=8)
plt.show()
edges = [brentq(lambda K: ch11.benard_free_free_Ra(K) - 2000.0, lo, hi) for lo, hi in ((0.3, Kff), (Kff, 7.0))]   # where Ra(K) = 2000
print("at Ra = 2000 the growing wavenumbers are", np.round(edges, 3), "| largest sigma there:",
      round(max(ch11.benard_free_free_sigma(K, 2000.0, 1.0)[0] for K in Kp), 2))
""",
    see="(a) On a logarithmic Ra axis: the free–free neutral curve for one half-wave across the layer (n = 1, purple), the one "
        "for two half-waves (dashed, sixteen times higher at its minimum), and the rigid curve as a grey dotted ghost between "
        "them; the diamond is the minimum. (b) The growth rate against K for three Rayleigh numbers: "
        "all negative (teal), touching zero at one K (purple), positive over a band (rose).",
    read="This is C01's picture with a real formula: as Ra rises the whole growth curve lifts, touches zero first at "
         "K = 2.22 when Ra = 657.5, and then a band of growing wavenumbers opens (its edges are printed).",
    change="…Pr = 7: the band edges stay where they are (σ = 0 does not involve Pr), the growth rates inside grow: $\\sigma=22.26$ instead "
           "of 11.02 at $K_c$ for Ra = 2000.",
    explain=r"""
1. `ch11.benard_free_free_Ra(K, n)` accepts arrays; `ch11.benard_free_free_sigma(K, Ra, Pr)` takes one $K$ at a time (hence the
   list comprehension) and returns $(\sigma_+,\sigma_-)$.
2. `brentq` finds the two wavenumbers where the neutral curve passes through Ra = 2000 — the edges of the growing band.""")
whatif(r"""
…the density depended on a second, slowly diffusing ingredient — salt? The same marginal problem reappears with Ra replaced by a
difference of two Rayleigh numbers, and diffusion can now *destabilise* a column that is lighter on top — C06.""")
# =====================================================================================================================
# A.5 §11.5 Double-Diffusive Instability — C06
# =====================================================================================================================
nb.section("11.5", "Double-Diffusive Instability", intro=r"""
**What is this section about?** Seawater's density depends on temperature and salt, and heat diffuses about a hundred times faster
than salt. That difference alone can make a column that is lighter on top overturn — in thin "salt fingers" or in oscillating
layers. We reuse the free–free Bénard solution to find the criterion, and meet the subtropical ocean's thermohaline staircases.""")
core("C06", r"The salt-finger criterion $\frac{gd^4}{\nu}\Big[\frac\beta{\kappa_s}\frac{dS}{dz}-\frac\alpha\kappa\frac{d\bar T}{dz}\Big]=657$ (11.46)",
     "How can a column that is lighter on top still overturn?")
problem(r"""
In the subtropical Atlantic, warm salty water from the surface sits on top of cooler, fresher water. The warmth makes the top
lighter, the salt makes it heavier, and the warmth wins: the column is stable by its density. Yet profilers find the water
arranged in "staircases" — uniform layers metres thick separated by sharp steps — the fingerprint of salt fingers. How can
diffusion, which usually smooths everything, start an overturn?""")
idea(r"""
warm, salty  ↑   a parcel pushed DOWN:   heat leaks in/out fast (κ)  → it takes the cold temperature of its new level
─────────────┼                           salt stays (κ_s ≈ κ/100)     → it keeps its extra salt
cold, fresh  ↓                           ⇒ colder AND saltier than its neighbours ⇒ heavier ⇒ keeps sinking: a FINGER
condition: the salt's destabilising effect (÷ its slow diffusivity κ_s) beats heat's stabilising one (÷ κ) by 27π⁴/4
""")
note("N53 [B]", r"""
**Two ingredients of density.** $\tilde\rho=\rho_0[1-\alpha(\tilde T-T_0)+\beta(\tilde s-s_0)]$ (the linear equation of state of
Ch. 1), with $\alpha$ the thermal expansion [1/K], $\beta$ the haline contraction [per g/kg], both positive, $s$ the salinity
[g/kg]; the salt diffusivity is $\kappa_s\approx1.5\times10^{-9}$ m²/s against $\kappa\approx1.4\times10^{-7}$ m²/s for heat:
$\tau\equiv\kappa_s/\kappa\approx0.011$.""")
code(r"""
from fluidpy import ch01_introduction as ch01                       # Ch. 1's linear seawater equation of state
rho = ch11.linear_eos(np.array([20.0, 10.0]), np.array([36.5, 35.0]))   # warm salty (20 C, 36.5 g/kg) and cold fresh (10 C, 35 g/kg) [kg/m^3]
print("densities [kg/m^3]:", np.round(rho, 3), "-> the warm salty water is lighter by", round(rho[1] - rho[0], 3))
rho_ch01 = ch01.seawater_density_linear(293.15, 36.5, alpha_T=2e-4, beta_S=7.6e-4)   # the same water through Ch. 1's function (T in kelvin)
assert np.isclose(rho_ch01, rho[0])                                 # one equation of state, two chapters
print("tau = kappa_s/kappa =", round(1.5e-9/1.4e-7, 4))
""", explain=r"""
1. `ch11.linear_eos(T, S)` evaluates the linear equation of state with $\alpha=2\times10^{-4}$ K⁻¹, $\beta=7.6\times10^{-4}$ per
   g/kg, reference 10 °C and 35 g/kg.
2. Warm salty water on top is *lighter* than the cold fresh water below it by about 0.9 kg/m³ — statically stable.
3. The same number from Ch. 1's function (which takes kelvin): the assertion passes.""")
note("N54 [B]", r"""
**Two regimes.**

| regime | arrangement | what a displaced parcel does | onset |
|---|---|---|---|
| **fingers** | hot salty over cold fresh ($d\bar T/dz>0$, $dS/dz>0$) | loses its temperature difference, keeps its salt, keeps going | stationary ($\sigma$ real): long thin fingers |
| **diffusive** | cold fresh over hot salty ($d\bar T/dz<0$, $dS/dz<0$) | keeps its salt but exchanges heat, overshoots and oscillates with growing amplitude | oscillatory ($\sigma$ complex, N06): convecting layers separated by sharp interfaces — Arctic water under sea ice |

The book gives no dispersion relation; ours (free–free walls, D13's equations with $\sigma$ kept) is the cubic
$(\sigma/\Pr+a^2)a^2(\sigma+a^2)(\sigma+\tau a^2)=K^2[-\mathrm{Ra}(\sigma+\tau a^2)+\tau\mathrm{Rs}(\sigma+a^2)]$ with
$a^2=\pi^2+K^2$ — Ra and Rs are defined in D13 (with §11.5's sign: Ra > 0 means warm on top).""")
code(r"""
s_f = ch11.double_diffusive_sigma(np.pi**2/2, 1000.0, 2000.0, 7.0, 0.0107)      # K^2 = pi^2/2, Ra = 1000, Rs = 2000, Pr = 7, tau = 0.0107
s_d = ch11.double_diffusive_sigma(np.pi**2/2, -2.0e4, -1.9e6, 7.0, 0.0107)      # cold fresh over hot salty (both gradients negative)
print("finger case, the three roots sigma:", np.round(s_f, 4))      # one real positive root: stationary growth
print("diffusive case, the three roots:   ", np.round(s_d, 2))      # a complex pair with positive real part: growing oscillations
print("onset type of the growing root:", "stationary" if abs(s_f[0].imag) < 1e-12 else "oscillatory", "|",
      "stationary" if abs(s_d[0].imag) < 1e-12 else "oscillatory")
""", explain=r"""
1. `ch11.double_diffusive_sigma(K2, Ra, Rs, Pr, tau)` returns the three roots of our cubic, largest real part first.
2. Fingers: one small positive *real* root — slow, stationary growth. Diffusive regime: a *complex pair* with positive real part —
   oscillations of growing amplitude (the "oscillatory" onset of N06), in a column that is statically stable.""")
P("P266", "uniqueness of a linear boundary-value problem", r"""
If two unknowns obey the same linear ODE with the same right-hand side and the same boundary conditions, and the homogeneous
problem has only the zero solution, they are equal: their difference solves the homogeneous problem, so it is zero. Here
$\big(\frac{d^2}{dz^2}-K^2\big)f=0$ with $f=0$ at both walls has only $f=0$ (its solutions are $e^{\pm Kz}$, which cannot vanish
twice).""", code=r"""
import numpy as np                                    # arrays
K = 2.0                                               # a wavenumber
M = np.array([[np.exp(-K/2), np.exp(K/2)], [np.exp(K/2), np.exp(-K/2)]])   # f = A e^{Kz} + B e^{-Kz} at z = -1/2 and z = +1/2
print(np.linalg.det(M))                               # not zero -> only A = B = 0 works: the solution is unique
""")
D("D13", ref="11.46", after=r"""
> ⚠️ **slip #10, second half — the book prints** in §11.5 $\mathrm{Ra}\equiv g\alpha d^4(d\bar T/dz)/(\nu\kappa)$ **; the correct
> form is** the same expression — it is not an error, but it is the *opposite sign* to §11.4's $\mathrm{Ra}=g\alpha\Gamma
> d^4/\kappa\nu$ with $\Gamma=-d\bar T/dz$: here Ra is negative when heated from below. That is why the third marginal equation
> carries $-\mathrm{Ra}$. In code: `ch11.rayleigh_number` (§11.4 sign) and `ch11.thermal_rayleigh_signed` (§11.5 sign).""")
nb.worked_example("a finger layer with easy Rayleigh numbers", r"""
§11.5 signs throughout.

1. Take $\mathrm{Ra}=1000$ (warm on top: stabilising) and $\mathrm{Rs}=2000$ (salty on top: destabilising).
2. Margin $\mathrm{Rs}-\mathrm{Ra}-27\pi^4/4=1000-657.5=342.5>0$ ⇒ finger-unstable.
3. Is the column statically stable? The density gradient is $\propto-\alpha T_z+\beta S_z$, and
   $\alpha T_z/(\beta S_z)=(\mathrm{Ra}\,\kappa)/(\mathrm{Rs}\,\kappa_s)=0.5/0.0107=46.7>1$: the heat term wins by far — lighter on top.
4. Still it overturns: the salt gradient counts $\kappa/\kappa_s\approx93$ times more in the criterion than in the density.""")
code(r"""
r = ch11.salt_finger_unstable(dTdz=0.01, dSdz=0.002, d=0.05, g=9.81)            # our thermocline: +0.01 K/m, +0.002 (g/kg)/m, a 5 cm layer
for name, val in r.items():                                         # every entry of the result
    print(f"  {name:15s} {val if isinstance(val, (bool, np.bool_)) else round(float(val), 3)}")
d_min = 0.05*(27*np.pi**4/4/r["lhs"])**0.25                         # the left side of (11.46) grows like d^4: thickness where it equals 27 pi^4/4
print(f"thinnest finger-unstable layer: d = {100*d_min:.2f} cm")
cases = {"hot salty over cold fresh": (0.01, 0.002, 0.05), "cold fresh over hot salty": (-0.01, -0.0027, 0.11),
         "warm on top, no salt": (0.01, 0.0, 0.05), "cold on top, weakly": (-0.005, 0.0, 0.05), "cold on top, strongly": (-0.01, 0.0, 0.05)}
for label, (Tz, Sz, d) in cases.items():                            # (dT/dz [K/m], dS/dz [(g/kg)/m], thickness [m])
    out = ch11.salt_finger_regime(Tz, Sz, d)                        # regime from the static stability and the cubic
    print(f"{label:28s} -> {out['regime']:11s} sigma_max = {np.round(out['sigma_max'], 2)}, lighter on top: {bool(out['density_stable'])}")
""", explain=r"""
1. `ch11.salt_finger_unstable(dTdz, dSdz, d)` evaluates the criterion with SI gradients ($z$ up): `lhs` is its left side (about
   6 × 10⁴, far above 657), `margin` = lhs − 27π⁴/4, `density_stable` says the column is lighter on top, `R_rho` $=\alpha T_z/(\beta
   S_z)=1.32$ is the *density ratio* oceanographers use (above 1: heat wins the density), and `Ra`, `Rs` are the two Rayleigh numbers.
2. The left side scales like $d^4$, so a layer thicker than about 1.6 cm is already finger-unstable.
3. `ch11.salt_finger_regime` labels the regime: `fingers`, `diffusive` (complex $\sigma$: growing oscillations), `stable`,
   `overturning`. Note the fourth row: cold on top is top-heavy, yet a thin, weakly heated layer is still **stable** — being
   top-heavy is not enough, the threshold $27\pi^4/4$ has to be passed (the fifth row passes it).""")
scratch(r"""
lhs = 9.81*0.05**4/1e-6*(7.6e-4*0.002/1.5e-9 - 2e-4*0.01/1.4e-7)    # (11.46) typed out: (g d^4/nu) [beta S_z/kappa_s - alpha T_z/kappa]
assert np.isclose(lhs, r["lhs"])                                    # same number as the library
s_one = np.max(np.real(ch11.double_diffusive_sigma(np.pi**2/2, -2000.0, 0.0, 1.0, 1.0)))   # no salt, tau = 1, heated from below (Ra_11.5 = -2000)
assert np.isclose(s_one, ch11.benard_free_free_sigma(np.pi/np.sqrt(2), 2000.0, 1.0)[0])    # = C05's free-free growth rate at Ra = 2000
print(f"left side of (11.46) by hand: {lhs:.0f};  single-component limit of the cubic: sigma = {s_one:.4f}")
""", r"""
1. The criterion in one line: the salt term divided by the slow diffusivity minus the heat term divided by the fast one.
2. Without salt and with equal diffusivities, the cubic must reduce to Bénard convection: its growth rate equals C05's 11.0155
   when §11.5's Ra = −2000 is §11.4's Ra = +2000.""")
fig(r"""
al, be, dd = 2e-4, 7.6e-4, 0.25                                    # alpha [1/K], beta [per g/kg], layer thickness [m] (a 25 cm layer)
n = 41 if FAST else 61                                             # grid points per axis
ax_T = np.linspace(-4e-6, 4e-6, n)                                 # alpha dT/dz [1/m]
ax_S = np.linspace(-4e-6, 4e-6, n)                                 # beta dS/dz [1/m]
code_of = {"stable": 0, "fingers": 1, "diffusive": 2, "overturning": 3}
reg = np.array([[code_of[ch11.salt_finger_regime(x/al, y/be, dd)["regime"]] for x in ax_T] for y in ax_S])   # regime at every grid point
from matplotlib.colors import ListedColormap                       # a colour map with one colour per regime
cmap = ListedColormap([COLORS["teal"], COLORS["amber"], COLORS["rose"], COLORS["muted"]])
fig, ax = plt.subplots(figsize=(6.4, 4.6))
ax.pcolormesh(ax_T*1e6, ax_S*1e6, reg, cmap=cmap, vmin=-0.5, vmax=3.5, alpha=0.45, shading="auto")
C_map = 27*np.pi**4/4*1e-6/(G0*dd**4)                              # the constant of (11.46), 27 pi^4 nu/(4 g d^4), for this thickness
ax.plot(ax_T*1e6, ax_T*1e6, color=C_BUOY, label="density line: $\\beta S_z=\\alpha T_z$")
ax.plot(ax_T*1e6, 1.5e-9*(ax_T/1.4e-7 + C_map)*1e6, color=C_SALT, lw=2, label="finger line (11.46)")
ax.plot(al*0.01*1e6, be*0.002*1e6, "o", color=COLORS["ink"])
ax.annotate("our thermocline", (al*0.01*1e6, be*0.002*1e6), (2.1, 0.55), fontsize=8, arrowprops=dict(arrowstyle="-", color=COLORS["ink"]))
for (x, y, txt) in ((2.0, -2.5, "stable"), (-2.2, 2.0, "overturning"), (2.9, 1.3, "fingers")):
    ax.text(x, y, txt, fontsize=9, ha="center")
ax.annotate("diffusive (oscillatory):\nthe thin band under the line", (-2.5, -2.72), (-1.3, -3.6), fontsize=8,
            arrowprops=dict(arrowstyle="-", color=COLORS["ink"]))
ax.set(xlabel=r"$\alpha\,d\bar T/dz$ [$10^{-6}$ m$^{-1}$]", ylabel=r"$\beta\,dS/dz$ [$10^{-6}$ m$^{-1}$]", xlim=(-4, 4), ylim=(-4, 4))
ax.set_title("Stable by density, unstable by diffusion", fontsize=11)
ax.legend(fontsize=8, loc="upper left")
savefig(fig, "ch11", "c06_regimes")
plt.show()
""",
    see="A map of regimes for a 25 cm layer of seawater: horizontal axis the temperature gradient (times α), vertical axis the salt "
        "gradient (times β). Teal: stable. Amber: salt fingers. Rose: the diffusive (oscillatory) regime. Grey: ordinary "
        "overturning. The blue diagonal is where the density does not change with height; below it the column is lighter on top.",
    read="The amber wedge between the nearly horizontal finger line and the blue density line is the surprise: lighter on top, yet "
         "unstable. It is wide because the finger line is almost the horizontal axis — the salt gradient is divided by the tiny "
         "salt diffusivity, so it counts about 93 times more than in the density. The thin rose band just under the blue line in the lower left is its "
         "mirror: cold fresh over hot salty, statically stable, growing oscillations. (For a much thinner layer every border "
         "moves away from the lines, because the threshold 27π⁴/4 then matters: the left side of the criterion scales like d⁴.)",
    change="…the two diffusivities were equal: the finger line would turn parallel to the density line and lie above it — the wedge "
           "closes. No double diffusion without two diffusivities (slider figure F3).",
    explain=r"""
1. Each grid point is one call of `ch11.salt_finger_regime` (static stability + the cubic of N54); the names are turned into
   integers 0–3 so that `pcolormesh` can colour them with a four-colour `ListedColormap`.
2. **Why 25 cm here, when the tiny example and slider F3 use 5 cm?** The borders are shifted away from the two lines by an
   offset proportional to $27\pi^4\nu/(4gd^4)$. At 25 cm that offset is invisible on these axes and the thin diffusive band can be
   resolved, so the map shows the regimes in their clean form; at 5 cm the offset is as large as the plotted gradients, which is
   exactly what F3 needs to show the wedge opening and closing.
3. The amber line is the criterion solved for the salt gradient: $\beta S_z=\frac{\kappa_s}{\kappa}\alpha T_z+\kappa_s\frac{27\pi^4\nu}{4gd^4}$.""")
nb.md(r"""
**Slider figure F3 — two diffusivities open the wedge.** The ratio $\tau=\kappa_s/\kappa$ is the slider; the wedge is the region
between the finger line and the density line.""")
nb.plotly(r"""
xT = np.linspace(-1e-6, 4e-6, 60)                                  # alpha dT/dz [1/m]
kap = 1.4e-7                                                       # thermal diffusivity [m^2/s]
C_f = 27*np.pi**4/4*1e-6/(G0*0.05**4)                              # the constant of (11.46) for a thin 5 cm layer, where it matters

def f3(tau):                                                       # the two lines and the wedge for one diffusivity ratio
    finger = tau*(xT + kap*C_f)                                    # finger line (11.46): beta S_z = tau (alpha T_z + kappa C)
    x0 = tau*kap*C_f/(1 - tau) if tau < 1 else np.inf              # where it meets the density line beta S_z = alpha T_z
    xs = np.linspace(x0, 4e-6, 40) if x0 < 4e-6 else np.array([np.nan])   # the wedge exists only to the right of that point
    wedge_x = np.concatenate([xs, xs[::-1]])                       # outline of the wedge: along the finger line, back along the density line
    wedge_y = np.concatenate([tau*(xs + kap*C_f), xs[::-1]])
    return {"density line (lighter on top below it)": (xT*1e6, xT*1e6), "finger line (11.46)": (xT*1e6, finger*1e6),
            "wedge: lighter on top, finger-unstable": (wedge_x*1e6, wedge_y*1e6)}

figF3 = slider_figure(f3, "κ_s/κ", np.round(np.geomspace(0.005, 1.0, 12 if FAST else 20), 4), xlabel="α dT/dz [1e-6 1/m]",
                      ylabel="β dS/dz [1e-6 1/m]", title="A 5 cm layer: the finger wedge closes as the diffusivities become equal",
                      xrange=(-1, 4), yrange=(-1, 4))
figF3.show()
""", explain=r"""
1. `f3(tau)` returns the density line, the finger line for that ratio, and the outline of the wedge between them.
2. At seawater's $\tau\approx0.011$ the wedge fills almost the whole lower-right triangle; at $\tau=1$ the finger line is parallel
   to the density line and above it — no wedge.""")
note("N57 [C]", r"""
**Fingers, staircases and layers** (named): at onset the cells are as wide as the layer, but well above it the fingers are long and
thin; a thick finger layer breaks into a staircase of convecting layers with fingers only at the steps (observed in the
subtropical North Atlantic and the Tyrrhenian Sea); heating a salt gradient from below builds a stack of diffusive layers (the
Arctic under sea ice). The two requirements: two diffusivities, and opposite contributions to the density gradient. Climate hook:
salt fingers mix heat and salt at different rates — a term ocean models parameterise separately (Ch. 12).""")
explainer("salt_fingers", "How can a stably stratified column overturn?", r"""
Dragging the two gradients across the regime map crosses the density-stable line and the finger line separately, so you see a
region that is stable by density yet finger-unstable; the parcel animation shows the parcel's temperature relaxing while its salt
stays; the diffusivity-ratio slider opens and closes the wedge.""", [
    "Preset 'subtropical thermocline': read 'statically stable but finger-unstable'.",
    "Drag κ_s/κ to 1: the wedge closes and the same gradients are stable.",
    "Preset 'Arctic': the σ-roots leave the real axis as a complex pair with positive real part — growing oscillations.",
])
whatif(r"""
…the "stratification" were not of heat or salt but of angular momentum — fluid spinning between two cylinders? The centrifugal
force plays gravity, and the same viscous threshold 1708 comes back — C07.""")

# =====================================================================================================================
# A.6 §11.6 Centrifugal Instability: Taylor Problem — R14 R15 R16, C07
# =====================================================================================================================
nb.section("11.6", "Centrifugal Instability: Taylor Problem", intro=r"""
**What is this section about?** Fluid between two coaxial cylinders, the inner one spinning: above a certain speed the smooth
circular flow of Ch. 8 breaks into a stack of doughnut-shaped vortices. Rotation acts like gravity — the centrifugal force pulls
fast-spinning fluid outward — and the threshold has the same 1708 as Bénard convection. In the atmosphere and ocean the same
mechanism is "inertial instability" (Ch. 13).""")
nb.recap("R14", "Rayleigh's circulation criterion", r"""
Ch. 8 flagged circular Couette flow as inviscidly stable when $\Omega_2/\Omega_1>(R_1/R_2)^2$. In general (Rayleigh 1916): an
inviscid swirling flow is unstable if the square of the circulation, $\Gamma^2=(2\pi rU_\theta)^2$, decreases outward somewhere —
$d\Gamma^2/dr<0$ — just as a fluid is statically unstable when $d\bar\rho/dz>0$. C07 adds the reason (the ring interchange,
N59). Here $\Gamma$ is a circulation [m²/s], not a temperature gradient.""", where="Ch. 8 §8.2")
nb.recap("R15", "Axisymmetric Navier–Stokes in cylindrical coordinates", r"""
With $\partial/\partial\varphi=0$: $\frac{D\tilde u_R}{Dt}-\frac{\tilde u_\varphi^2}R=-\frac1\rho\frac{\partial\tilde p}{\partial R}+
\nu\Big(\nabla^2\tilde u_R-\frac{\tilde u_R}{R^2}\Big)$, $\frac{D\tilde u_\varphi}{Dt}+\frac{\tilde u_R\tilde u_\varphi}R=\nu\Big(\nabla^2
\tilde u_\varphi-\frac{\tilde u_\varphi}{R^2}\Big)$, $\frac{D\tilde u_z}{Dt}=-\frac1\rho\frac{\partial\tilde p}{\partial z}+\nu\nabla^2\tilde u_z$,
$\frac1R\frac{\partial}{\partial R}(R\tilde u_R)+\frac{\partial\tilde u_z}{\partial z}=0$ (11.47), with
$\frac D{Dt}=\partial_t+\tilde u_R\partial_R+\tilde u_z\partial_z$ and $\nabla^2=\partial_R^2+\frac1R\partial_R+\partial_z^2$ — Ch. 4's
cylindrical operators ($R$ radius, $\varphi$ angle, $z$ along the axis).""", where="Ch. 4 Appendix B, Ch. 8 P186")
slip(12, r"""continuity (in the axisymmetric set and again in the linearised one) as $\frac{\partial}{\partial R}(R\tilde u_R)+
\frac{\partial\tilde u_z}{\partial z}=0$, whose two terms do not have the same units (m/s against 1/s)""",
     r"""$\frac1R\frac{\partial}{\partial R}(R\tilde u_R)+\frac{\partial\tilde u_z}{\partial z}=0$.""")
nb.recap("R16", "Circular Couette flow", r"""
Between cylinders $R_1<R<R_2$ turning at $\Omega_1$, $\Omega_2$: $U_R=U_z=0$, $U_\varphi=AR+B/R$,
$\frac1\rho\frac{dP}{dR}=\frac{U_\varphi^2}R$ (11.49) with $A\equiv\frac{\Omega_2R_2^2-\Omega_1R_1^2}{R_2^2-R_1^2}$,
$B\equiv\frac{(\Omega_1-\Omega_2)R_1^2R_2^2}{R_2^2-R_1^2}$ — Ch. 8's solution $u_\varphi(R)=AR+\frac BR$ (8.9), computed by
`LAM.circular_couette(…, return_coeffs=True)`.""", where="Ch. 8 §8.2")

core("C07", r"Taylor–Couette onset: $\mathrm{Ta}_{cr}=\frac{1708}{\tfrac12(1+\Omega_2/\Omega_1)}$ (11.54) from the narrow-gap equations and the Taylor number $\mathrm{Ta}=4\big(\frac{\Omega_1R_1^2-\Omega_2R_2^2}{R_2^2-R_1^2}\big)\frac{\Omega_1d^4}{\nu^2}$ (11.52)",
     "Why does spinning the inner cylinder stack up vortices, while spinning the outer one does not?")
problem(r"""
G. I. Taylor (1923) filled the gap between two glass cylinders with water and turned the inner one. At low speed the water simply
swirls in circles. Past a sharp speed, the flow organises into a stack of counter-rotating ring vortices, each about as tall as
the gap is wide. Turn the outer cylinder instead and nothing happens. The same physics governs a journal bearing, a Couette
viscometer — and, with Earth's rotation, inertial instability of jets in the atmosphere and ocean (Ch. 13).""")
idea(r"""
swap two thin fluid rings (equal mass) at r₁ < r₂, each keeping its circulation Γ = 2πrU_θ (Kelvin):
    E = Γ²/(8π²r²) per unit mass   ⇒   ΔE = (Γ₂² − Γ₁²)(1/r₁² − 1/r₂²)/(8π²)
    Γ² falls outward (Γ₂² < Γ₁²)  ⇒ ΔE < 0: energy RELEASED ⇒ unstable   (inner cylinder spinning)
    Γ² rises outward              ⇒ ΔE > 0: you must PAY         ⇒ stable     (outer cylinder spinning)
centrifugal force ↔ gravity,  Γ² ↔ −ρ̄  — and viscosity adds a threshold, Ta_c ≈ 1708 / ((1 + Ω₂/Ω₁)/2)
""")
remind("C07")
note("N58 [C]", r"""
**The set-up**: coaxial cylinders, radii $R_1<R_2$ [m], angular speeds $\Omega_1$, $\Omega_2$ [rad/s] (counter-clockwise
positive), gap $d=R_2-R_1$, ratio $\mu\equiv\Omega_2/\Omega_1$; experiments show the first instability is axisymmetric
($\partial/\partial\varphi=0$), so we only look at disturbances that are rings. Wavy vortices come later (N69).""")
note("N59 [B]", r"""
**The ring interchange.** A ring of fluid at radius $r$ with circulation $\Gamma=2\pi rU_\theta$ has kinetic energy per unit mass
$E=U_\theta^2/2=\Gamma^2/(8\pi^2r^2)$. Swap rings 1 and 2 (equal masses), each keeping its own circulation (Kelvin's theorem):
$E_{\rm final}=\frac1{8\pi^2}\big[\frac{\Gamma_2^2}{r_1^2}+\frac{\Gamma_1^2}{r_2^2}\big]$,
$E_{\rm initial}=\frac1{8\pi^2}\big[\frac{\Gamma_1^2}{r_1^2}+\frac{\Gamma_2^2}{r_2^2}\big]$, so
$\Delta E=\frac1{8\pi^2}(\Gamma_2^2-\Gamma_1^2)\Big(\frac1{r_1^2}-\frac1{r_2^2}\Big)$ — since $r_2>r_1$ the second bracket is
positive and the sign of $\Delta E$ is the sign of $\Gamma_2^2-\Gamma_1^2$. Energy released ⇒ the swap happens on its own ⇒
unstable. (Notation: this argument and recap R14 write $r$, $U_\theta$; the derivation D14 below writes $R$, $U_\varphi$ — the same radius
and the same swirl velocity, in the symbols the book uses in each passage.)""")
code(r"""
print(ch11.ring_interchange_energy(4.0, 2.0, 1.0, 2.0))             # circulations 4 and 2 m^2/s at radii 1 and 2 m: energies [m^2/s^2] per unit mass
Rr = np.linspace(0.1, 0.105, 11)                                    # radii across a 5 mm gap [m]
inner = ch11.rayleigh_circulation_criterion(Rr, LAM.circular_couette(Rr, 0.1, 0.105, 0.5, 0.0))   # only the inner cylinder turns (0.5 rad/s)
outer = ch11.rayleigh_circulation_criterion(Rr, LAM.circular_couette(Rr, 0.1, 0.105, 0.4, 0.5))   # the outer cylinder turns faster
print("inner cylinder spinning: Rayleigh-unstable?", inner["unstable"], "| outer faster:", outer["unstable"])
print("Rayleigh line for this gap: mu = (R1/R2)^2 =", round(ch11.couette_rayleigh_line(0.1, 0.105), 4))
""", explain=r"""
1. `ch11.ring_interchange_energy(Gamma1, Gamma2, r1, r2)` returns the energies before and after the swap and their difference:
   $\Delta E=-0.114$ m²/s² — released, because the circulation falls outward.
2. `ch11.rayleigh_circulation_criterion(r, U_theta)` checks the sign of $d\Gamma^2/dr$ on a sampled profile; `LAM.circular_couette`
   (Ch. 8) supplies the profile.
3. With only the inner cylinder turning the flow is Rayleigh-unstable everywhere; with the outer one faster it is stable. The
   border is $\mu=(R_1/R_2)^2=0.907$ for our gap.""")
note("N61 [B]", r"""
**Geometry**: a cut through the gap along the axis — the vortices are stacked rings; neighbouring rings turn opposite ways.""")
fig(r"""
fig, ax = plt.subplots(figsize=(5.2, 3.6))
taylor_sketch(ax)                                                  # drawing only: the gap and a stack of counter-rotating vortices
plt.show()
""",
    see="A meridional cut through the gap between the cylinders (the axis is to the left): a stack of counter-rotating cells.",
    read="Neighbouring cells turn opposite ways; a pair spans one axial wavelength 2π/k, about twice the gap width.",
    change="…the outer cylinder turns faster than the Rayleigh line allows: no cells at any speed.")
P("P267", "narrow-gap (small-curvature) approximation: expand in d/R and keep the leading order", r"""
When the gap $d$ is much smaller than the radius $R_1$, curvature terms such as $\frac1R\frac{d}{dR}$ or $\frac1{R^2}$ are smaller
than $\frac{d^2}{dR^2}$ by factors $d/R_1$, $(d/R_1)^2$; we keep the leading order and drop the rest — the same idea as Ch. 8's
lubrication scaling (P188). With $x=(R-R_1)/d\in[0,1]$, $\frac d{dR}=\frac1d\frac d{dx}$.""", code=r"""
R1, d = 0.10, 0.005                                   # 10 cm inner radius, 5 mm gap
print(d/R1, (d/R1)**2)                                # 0.05 and 0.0025: the sizes of the dropped terms relative to the kept one
""")
D("D14", ref="11.51")
code(r"""
tp = ch11.taylor_perturbation_sympy()                               # D14 redone by the library's sympy engine (cached)
print("linearised equations (residuals R, phi, z):", tp["lin_R"], tp["lin_phi"], tp["lin_z"])   # 0 0 0: step 3-4
print("dU/dR + U/R - 2A:", tp["term_2A"], "| operator identity of step 10:", tp["operator_identity"])   # 0 0
print("Ta coefficient (step 15):", tp["Ta_coefficient"])           # -4 A Omega1 d^4/nu^2 with d = R2 - R1
print("continuity as printed has consistent units?", tp["printed_continuity_units_ok"], "| corrected:", tp["correct_continuity_units_ok"])
""", explain=r"""
1. `ch11.taylor_perturbation_sympy()` linearises the axisymmetric equations about Couette flow and follows D14's eliminations;
   every residual is 0.
2. The last line is slip #12 as a units check: the printed continuity equation fails it, the corrected one passes.""")
D("D15", ref="11.52")
note("N66 [B] · N67 [B]", r"""
**Conditions and the marginal state.** No slip on both cylinders: $\hat u_R=\frac{d\hat u_R}{dR}=\hat u_\varphi=0$ at $x=0,1$
(11.53). Taylor assumed the onset is stationary ($\sigma=0$); this is proved for cylinders turning the same way ($\mu\ge0$) and
checked numerically here, but not in general — for counter-rotation ($\mu<0$) the $\sigma=0$ curve is only what the narrow-gap
equations give *if* $\sigma=0$ is assumed.""")
code(r"""
for mu in (0.0, 0.5):                                               # outer cylinder at rest, and co-rotating at half the inner speed
    s = ch11.taylor_growth_rate(3.12, 4000.0, mu)                   # leading growth rate at k = 3.12, Ta = 4000 [units nu/d^2]
    print(f"mu = {mu}: sigma = {s.real:.4f} + {abs(s.imag):.1e} i  ->  {'unstable' if s.real > 0 else 'stable'}, stationary (no imaginary part)")
""")
nb.worked_example("water between two cylinders", r"""
$R_1=10$ cm, $R_2=10.5$ cm (so $d=5$ mm and $R_2/R_1=1.05$, our choice), water $\nu=10^{-6}$ m²/s, outer cylinder at rest, inner
at $\Omega_1=0.5$ rad/s.

1. Taylor number: $\mathrm{Ta}=4\frac{\Omega_1R_1^2-0}{R_2^2-R_1^2}\cdot\frac{\Omega_1d^4}{\nu^2}=4\times\frac{0.5\times0.01}{0.001025}\times
   \frac{0.5\times6.25\times10^{-10}}{10^{-12}}=19.51\times312.5=6098$.
2. Narrow-gap shortcut $2(\Omega_1R_1d/\nu)^2(d/R_1)=2\times250^2\times0.05=6250$ — 2.5 % higher (the order $d/R_1$ we dropped).
3. Critical value for $\mu=0$: the approximate formula gives $1708/0.5=3416$; the exact narrow-gap value is 3389.9 (0.77 % lower).
4. $6098>3390$ ⇒ vortices.
5. Onset speed: $\mathrm{Ta}\propto\Omega_1^2$, so $\Omega_{1,c}=0.5\sqrt{3389.9/6098}=0.373$ rad/s (about 3.6 revolutions per
   minute).""")
code(r"""
Ta = ch11.taylor_number(0.5, 0.0, 0.1, 0.105, 1e-6)                 # (11.52): Omega1 = 0.5 rad/s, outer at rest, R1 = 0.1 m, R2 = 0.105 m, water
print(Ta, "| narrow-gap shortcut:", round(ch11.taylor_number_narrow_inner(0.5, 0.1, 0.005, 1e-6), 1))
for mu in (1.0, 0.5, 0.0, -0.5):                                    # co-rotating ... outer at rest ... counter-rotating
    c = ch11.taylor_critical(mu)                                    # exact narrow-gap minimum over k (P263's recipe)
    a = ch11.taylor_critical_approx(mu)                             # the approximate formula (11.54)
    print(f"mu = {mu:5.2f}: Ta_c = {c['Ta_c']:8.2f} at k_c = {c['k_c']:.4f} | (11.54) gives {a:7.1f} ({100*(a/c['Ta_c'] - 1):+.2f} %)")
print(f"onset speed for our cylinders: Omega_1 = {0.5*np.sqrt(ch11.taylor_critical(0.0)['Ta_c']/Ta['Ta']):.3f} rad/s")
""", explain=r"""
1. `ch11.taylor_number` returns a dictionary: `Ta` (6098 here), `rayleigh_stable` (True on and beyond the Rayleigh line, where
   Ta ≤ 0 and nothing drives the flow) and `mu`.
2. `ch11.taylor_critical(mu)` minimises the marginal Taylor number over $k$; `ch11.taylor_critical_approx(mu)` is the approximate
   formula. At $\mu=1$ the exact value is Bénard's 1707.76 (D15); the approximation is within 1 % for co-rotation and 6.5 % too
   high at $\mu=-0.5$ — it is a co-rotation formula.""")
scratch(r"""
O1, R1, R2, nu = 0.5, 0.1, 0.105, 1e-6                             # the worked example's numbers (SI)
Ta_mine = 4*(O1*R1**2 - 0.0*R2**2)/(R2**2 - R1**2)*O1*(R2 - R1)**4/nu**2      # (11.52) typed out
assert np.isclose(Ta_mine, ch11.taylor_number(O1, 0.0, R1, R2, nu)["Ta"])     # same number as the library
assert np.isclose(1708/(0.5*(1 + 0.5)), ch11.taylor_critical_approx(0.5))     # (11.54) typed out for mu = 0.5
dE_mine = (2.0**2 - 4.0**2)*(1/1.0**2 - 1/2.0**2)/(8*np.pi**2)                # N59's energy change for circulations 4, 2 at radii 1, 2
assert np.isclose(dE_mine, ch11.ring_interchange_energy(4.0, 2.0, 1.0, 2.0)["dE"])
assert np.isclose(ch11.taylor_critical(1.0)["Ta_c"], ch11.benard_critical()["Ra_c"], rtol=1e-6)   # mu -> 1 is Benard (D15)
print(f"Ta by hand = {Ta_mine:.1f}; ring swap dE = {dE_mine:.5f} m^2/s^2; Ta_c(mu = 1) = Ra_c of Benard")
""", r"""
1. The Taylor number, the approximate critical value and the ring-swap energy typed out and compared with the library.
2. The last assertion is D15's claim with numbers: the co-rotating limit of the Taylor problem has the same eigenvalue as rigid
   Bénard convection.""")
fig(r"""
tt = ch11.taylor_critical_table()                                  # exact narrow-gap Ta_c(mu) for mu = -0.5 ... 1 (read from reference/ch11)
sb = ch11.taylor_stability_boundary(mus=np.linspace(-1.0, 0.85, 6 if FAST else 8))   # the marginal curve for our radius ratio 1.05
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.7))
a.plot(tt["mu"], tt["Ta_c"], "o", ms=3.5, color=C_NEUT, label="exact narrow gap")
a.plot(tt["mu"], tt["approx"], "--", color=C_BASE, label="approximate formula (11.54)")
a.plot(1.0, BENCH["benard_rigid_rigid"]["Ra_c"], "D", color=COLORS["ink"], ms=5)
a.text(0.52, 1150, "Bénard's 1707.76", fontsize=8)
a.set(xlabel=r"$\mu=\Omega_2/\Omega_1$ [–]", ylabel="critical Taylor number [–]", ylim=(0, 7200))
a2 = a.twinx()                                                     # a second vertical axis for the error
a2.plot(tt["mu"], 100*tt["rel_error"], color=C_GROW, lw=1)
a2.set_ylabel("error of (11.54) [%]", color=C_GROW)
a2.set_ylim(0, 8)
a.set_title("(a) threshold: exact and approximate", fontsize=10)
a.legend(fontsize=8, loc="upper right")
b.fill_between(sb["x_outer"], sb["y_inner"], 9000, color=C_GROW, alpha=0.12)
b.plot(sb["x_outer"], sb["y_inner"], color=C_NEUT, label="viscous threshold")
b.plot(sb["rayleigh_x"], sb["rayleigh_y"], "--", color=C_BUOY, label="Rayleigh line")
b.text(-5200, 7600, "vortices", color=C_GROW)
b.text(-4800, 2200, "stable", color=C_DECAY)
b.set(xlabel=r"outer cylinder $\Omega_2R_2^2/\nu$ [–]", ylabel=r"inner cylinder $\Omega_1R_2^2/\nu$ [–]", xlim=(-7000, 7000), ylim=(0, 9000))
b.set_title("(b) who must spin how fast ($R_2/R_1$ = 1.05)", fontsize=10)
b.legend(fontsize=8, loc="lower left")
savefig(fig, "ch11", "c07_taylor")
plt.show()
live_one = ch11.taylor_critical(0.25)                              # one value recomputed now, to compare with the table
i25 = int(np.argmin(np.abs(tt["mu"] - 0.25)))                      # the table row for mu = 0.25
print(f"table: Ta_c(0.25) = {tt['Ta_c'][i25]:.2f}; recomputed now: {live_one['Ta_c']:.2f}")
assert np.isclose(tt["Ta_c"][i25], live_one["Ta_c"], rtol=6e-4)    # the stored table has 4-6 significant figures
""",
    see="(a) The critical Taylor number against the rotation ratio: purple dots exact, grey dashed the approximate formula, the thin "
        "rose line (right axis) its error. (b) The plane of the two cylinder speeds for our radius ratio: above the purple curve "
        "vortices appear (rose); the blue dashed line is Rayleigh's inviscid border.",
    read="(a) At μ = 1 the threshold is Bénard's 1707.76; it rises as the outer cylinder slows and rises steeply for "
         "counter-rotation, where the approximate formula drifts off: 6.5 % too high at μ = −0.5. (b) With the outer cylinder at rest "
         "(the vertical axis) the inner one must exceed about 4100 in these units. To the right (co-rotation) the purple curve (a polygon through 8 computed points) rises "
         "with the Rayleigh line and stays above it — by roughly 1500 units at the right edge; to the left (counter-rotation) the flow is Rayleigh-unstable for *any* inner speed, yet "
         "viscosity holds it until the curve is crossed.",
    change="…a wider gap (say a radius ratio 1.5): the narrow-gap theory no longer applies and the true curve moves (named).",
    explain=r"""
1. `ch11.taylor_critical_table()` reads our stored table of exact narrow-gap minima (31 values of μ); the last three lines
   recompute one of them live and compare (`rtol` ≥ 6e-4 because the table is stored with few digits).
2. `ch11.taylor_stability_boundary(mus=…)` converts $\mathrm{Ta}_c(\mu)$ into the two cylinder speeds for the radius ratio 1.05,
   and returns the Rayleigh line $\Omega_1/\Omega_2=R_2^2/R_1^2$.
3. `a.twinx()` adds a second vertical axis that shares the horizontal one.""")
fig(r"""
kc0 = ch11.taylor_critical(0.0)["k_c"]                             # critical axial wavenumber for mu = 0 [1/d]
tv = ch11.taylor_eigenfunction(kc0, 0.0, x=np.linspace(0, 1, 41), z=np.linspace(0, 2*2*np.pi/kc0, 121))   # fields over two wavelengths
fig, ax = plt.subplots(figsize=(3.6, 5.2))
pc = ax.pcolormesh(tv["x"], tv["z"], tv["u_phi"]/np.max(np.abs(tv["u_phi"])), cmap="PuOr_r", shading="auto")   # swirl disturbance (scaled)
ax.contour(tv["x"], tv["z"], tv["psi"], 8, colors=COLORS["ink"], linewidths=0.7)   # streamlines of the meridional flow
fig.colorbar(pc, ax=ax, label=r"swirl disturbance $u_\varphi$ (scaled)")
ax.set(xlabel="$x=(R-R_1)/d$", ylabel="$z$ [$d$]")
ax.set_title("Taylor vortices at onset", fontsize=10)
plt.show()
print(f"axial wavelength 2 pi/k_c = {2*np.pi/kc0:.3f} d  ->  each vortex is {np.pi/kc0:.3f} d tall")
""",
    see="Streamlines (black) of the flow in the plane through the axis, over two wavelengths; colour shows where the swirl is "
        "faster (orange) or slower (purple) than Couette flow. The inner cylinder is on the left.",
    read="Each cell is about as tall as the gap is wide. Where the cells carry fluid outward (away from the inner cylinder) they "
         "bring fast-swirling fluid with them — the orange jets; inflow brings slow fluid — the purple ones. That is the "
         "centrifugal analogue of warm plumes rising.",
    change="…μ → 1: the cells become exactly Bénard's rolls (D15), with the swirl disturbance in the role of temperature.",
    explain=r"""
1. `ch11.taylor_eigenfunction(k, mu, x, z)` returns the marginal mode on an $(x,z)$ grid: radial velocity `u_R`, swirl disturbance
   `u_phi`, axial velocity `u_z` and the stream function `psi` of the meridional flow ($u_R=\partial\psi/\partial z$).""")
nb.md(r"""
**Slider figure F4 — the marginal curve for each rotation ratio.** For every $\mu$: the exact narrow-gap curve $\mathrm{Ta}(k)$,
its minimum, and the constant value of the approximate formula.""")
nb.plotly(r"""
kT = np.linspace(2.0, 4.6, 12 if FAST else 17)                     # axial wavenumbers [1/d]

def f4(mu):                                                        # traces for one rotation ratio
    Ta_k = np.array([ch11.taylor_marginal_Ta(k, mu) for k in kT])  # exact narrow-gap marginal Taylor number at each k
    j = int(np.argmin(Ta_k))                                       # the lowest point of the curve
    appr = ch11.taylor_critical_approx(mu) if mu > -1 else np.nan  # the approximate formula (11.54)
    return {"exact narrow gap Ta(k)": (kT, Ta_k), "minimum": (np.array([kT[j]]), np.array([Ta_k[j]])),    # the curve and its lowest point
            "approximate formula (11.54)": (kT, np.full(len(kT), appr))}                                # a horizontal line at 1708/((1 + mu)/2)

figF4 = slider_figure(f4, "μ = Ω₂/Ω₁", np.round(np.linspace(-0.5, 1.0, 7 if FAST else 11), 3), xlabel="axial wavenumber k [1/d]",   # slider: mu
                      ylabel="marginal Taylor number", title="Counter-rotation raises the threshold; co-rotation approaches 1708",   # labels
                      yrange=(1500, 9000), modes={"minimum": "markers"})   # fixed range; the minimum drawn as a marker
figF4.show()
""", explain=r"""
1. `ch11.taylor_marginal_Ta(k, mu)` solves the narrow-gap equations at $\sigma=0$ for the Taylor number at one $k$.
2. At $\mu=1$ the minimum is 1707.76; drag left and the valley rises and its floor moves to larger $k$.""")
note("N68 [B]", r"""
**Theory against Taylor's experiment**: Taylor's measured onset speeds lay on his narrow-gap curve — one of the first quantitative
triumphs of linear stability theory. Our panel (b) above redraws the theory for our radius ratio (the book's ratio is not used);
the agreement shown in the book is with Taylor's own data.""")
note("N124 [B]", r"""
**A second route (Exercises 11.8–11.9).** At $\sigma=0$ the narrow-gap pair reads
$\big(\frac{d^2}{dR^2}-k^2\big)^2\hat u_R=(1+\alpha x)\hat u_\varphi$ (11.92) and
$\big(\frac{d^2}{dR^2}-k^2\big)\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R$ (11.93, corrected); expanding
$\hat u_\varphi=\sum_{m=1}^\infty C_m\sin(m\pi x)$ (11.94) and projecting on each sine (a Galerkin method — Ch. 10's idea with
sines instead of hats) gives the critical value with four modes to 0.01 %.""")
slip(7, r"""the second equation of Exercise 11.9 as $\big(\frac{d^2}{dR^2}-k^2\big)^2\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R$""",
     r"""$\big(\frac{d^2}{dR^2}-k^2\big)\hat u_\varphi=-\mathrm{Ta}\,k^2\hat u_R$ — the first power, as the narrow-gap equations give
at $\sigma=0$; with the square the problem would be of eighth order.""")
code(r"""
Ta_gal = ch11.taylor_galerkin_Ta(kc0, 0.0, n_modes=4)               # four sine modes (Exercise 11.9's route)
Ta_ex = ch11.taylor_marginal_Ta(kc0, 0.0)                           # the Chebyshev eigenvalue
Ta_pr = ch11.taylor_galerkin_Ta(kc0, 0.0, n_modes=4, printed=True)  # the same projection with the printed (squared) operator of slip #7
print(f"Galerkin, 4 sines: {Ta_gal:.2f} | Chebyshev: {Ta_ex:.2f} | difference {100*abs(Ta_gal/Ta_ex - 1):.3f} %")
print("with the printed operator:", Ta_pr, "(nan = no positive real Taylor number exists for that form)")
""", explain=r"""
1. Two independent discretisations of the same marginal problem agree to about 0.01 %.
2. With the printed square on the operator, the projected problem has no positive real eigenvalue at all (the function returns
   `nan`, "not a number") — the printed form cannot be the equation that was meant.""")
note("N69 [C]", r"""
**After onset** (named): as $\Omega_1$ rises the vortices become wavy ($\partial/\partial\varphi\neq0$), then modulated, then
turbulent (Coles 1965); the steady equations have many solutions for the same speeds (non-uniqueness). Relatives: Dean vortices in
curved channels, Görtler vortices on concave walls. Pointer: the route to turbulence, N111 and Ch. 12.""")
explainer("taylor_couette_onset", "Why do stacked vortices appear between spinning cylinders?", r"""
Dragging the rotation ratio and the Taylor number moves the point across the Rayleigh line and the viscous boundary separately;
the profile view shows where the squared circulation falls outward and the ring-pair inspector gives the energy change for any two
radii; vortices appear only past the viscous curve.""", [
    "Preset 'inner only': raise Ta until the vortices appear — compare the onset with the approximate formula's 3416 and the exact 3390.",
    "Preset 'counter-rotation': the status says Rayleigh-unstable but viscously stable — explain why with the ring inspector near each wall.",
    "Drag μ to 1: the threshold becomes Bénard's 1708.",
])
whatif(r"""
…the "restoring" stratification and the "driving" shear lived in the same fluid at the same place — a current over a stably
stratified thermocline? Then gravity and shear compete continuously across the layer, and one ODE (Taylor–Goldstein) decides — C08.""")
# =====================================================================================================================
# A.7 §11.7 Instability of Continuously Stratified Parallel Flows — R17 R18 R19, C08 · R20, C09 · C10
# =====================================================================================================================
nb.section("11.7", "Instability of Continuously Stratified Parallel Flows", intro=r"""
**What is this section about?** C02 had a sharp interface; real oceans and atmospheres have smooth profiles of velocity $U(z)$ and
density $\bar\rho(z)$. We derive one equation for small waves on such a flow (Taylor–Goldstein, C08), then prove two theorems that
need no solving: if the gradient Richardson number exceeds ¼ everywhere the flow is stable (C09), and any unstable wave speed
lies inside a semicircle set by the slowest and fastest fluid (C10).""")
nb.recap("R17", "The stratified parallel flow", r"""
Basic state plus disturbance: $\tilde{\mathbf u}=U(z)\mathbf e_x+\mathbf u$, $\tilde p=P+p$, $\tilde\rho=\bar\rho(z)+\rho$, with the
inviscid Boussinesq momentum equation $\frac{\partial\tilde{\mathbf u}}{\partial t}+(\tilde{\mathbf u}\cdot\nabla)\tilde{\mathbf
u}=-\frac1{\rho_0}\nabla\tilde p-g\frac{\bar\rho+\rho}{\rho_0}\mathbf e_z$ and the hydrostatic balance of the basic state,
$0=-\frac1{\rho_0}\frac{\partial P}{\partial z}-g\frac{\bar\rho}{\rho_0}$ — Ch. 4's Boussinesq set and Ch. 7 §7.8's base state. Only
two-dimensional disturbances are considered (Squire's theorem, proved in C11 for unstratified flow, is *assumed* here).""",
         where="Ch. 4 §4.9, Ch. 7 §7.8")
nb.recap("R18", "The buoyancy frequency", r"""
$N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ (7.127) [1/s²] (Ch. 7 §7.8; Ch. 1's $N^2=-\frac g\rho\big(\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\big)$
(1.29) for an incompressible fluid); with it the linearised density equation becomes $\frac{\partial\rho}{\partial t}+U\frac{\partial
\rho}{\partial x}-\frac{\rho_0N^2w}g=0$ (11.56). ⚠️ slip #13: where §11.7 recalls $N^2$ the book prints the cross-reference
"(7.128)"; in Ch. 7 the definition $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ carries the number (7.127), and (7.128) is the
momentum equation after it — a slip of the cross-reference only, the formula is right.""", where="Ch. 1 §1.10, Ch. 7 §7.8")
nb.code(r"""
N2_thermo = STRAT.brunt_vaisala_sq(1025.0, -0.002*1025.0/10.0, 0.0)   # a thermocline: density falls by 0.2 % over 10 m upward [1/s^2]
print(f"N^2 = {N2_thermo:.2e} 1/s^2,  N = {np.sqrt(N2_thermo):.4f} rad/s,  period 2 pi/N = {2*np.pi/np.sqrt(N2_thermo)/60:.1f} min")
""", explain=r"""
1. `STRAT.brunt_vaisala_sq(rho0, drho_dz, drho_a_dz)` is Ch. 1's function: reference density, density gradient [kg/m⁴] and the
   adiabatic gradient (zero for water here).
2. A 0.2 % density change over 10 m gives $N^2\approx2\times10^{-3}$ s⁻² — a parcel would bob with a period of about two minutes.""")
nb.recap("R19", "The stream function, §11.7 version", r"""
$u=\frac{\partial\psi}{\partial z}$, $w=-\frac{\partial\psi}{\partial x}$ (11.57) — Ch. 7's sign. ⚠️ slip #9: §11.8 uses
$u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$ ($y$ is the cross-stream coordinate there) and §11.14 uses
$u=-\partial\psi/\partial z$, $w=\partial\psi/\partial x$; each block repeats its own.""", where="Ch. 4 §4.3, Ch. 7 §7.8")

core("C08", r"The Taylor–Goldstein equation $(U-c)\big(\frac{d^2}{dz^2}-k^2\big)\hat\psi-\frac{d^2U}{dz^2}\hat\psi+\frac{N^2}{U-c}\hat\psi=0$ (11.61)",
     "What single equation decides whether a stratified current is stable?")
problem(r"""
A tidal current flows over a stably stratified thermocline; a jet stream blows over a cold inversion. The velocity changes smoothly
with height, and so does the density. Shear wants to roll the layers up; the density stratification resists lifting heavy water.
We want one equation that contains both, for one wave at a time.""")
idea("", r"""
Exactly C02's recipe with smooth profiles: linearise, use a stream function so continuity is automatic, insert
$e^{\mathrm ik(x-ct)}$, eliminate pressure and density. One second-order ODE in $z$ remains, with $c$ as the eigenvalue: $U''$ is
the shear profile's curvature (C12 will make it the hero), $N^2/(U-c)$ the stratification. With $N^2=0$ it is Rayleigh's equation
of C12.""")
note("N70 [C]", r"""
**History.** Taylor (1915) conjectured from examples that a gradient Richardson number below ¼ is needed for instability; Prandtl,
Goldstein, Richardson, Synge and Chandrasekhar suggested other values; Miles (1961) proved Taylor's conjecture and Howard (1961)
gave the short proof we follow (C09).""")
D("D16", ref="11.57", after=r"""
> ⚠️ **slip #4 — the book prints** the pressure term of the $w$-equation as $-\frac1{\rho_0}\frac{\partial p}{\partial x}$ **; the
> correct form is** $-\frac1{\rho_0}\frac{\partial p}{\partial z}$ — the vertical momentum equation needs the vertical pressure
> gradient, and the next equation in the book indeed uses it (step 5).""")
P("P269", "singular point of an ODE", r"""
Where the coefficient of the highest derivative vanishes, an ODE is *singular*: in the Taylor–Goldstein and Rayleigh equations that
happens where $U(z)=c$. For a growing mode ($c_i\neq0$) it never happens; for a neutral mode (real $c$ inside the velocity range)
the solution may have a kink or a logarithm there — the critical layer of C12. Viscosity (C11) removes the singularity by restoring
the fourth derivative.""", code=r"""
import numpy as np                                    # arrays
z = np.linspace(-1, 1, 5)                             # five heights
U, c = z, 0.5                                         # Couette flow U = z and a real wave speed inside its range
print(U - c)                                          # changes sign: the coefficient (U - c) vanishes at z = 0.5
""")
D("D17", ref="11.61")
note("N73 [B] · N74 [B]", r"""
**Pairs again.** The Taylor–Goldstein equation contains no $\mathrm i$, so if $(\hat\psi,c)$ solves it so does $(\hat\psi^*,c^*)$:
every growing mode has a decaying twin, and any $c_i\neq0$ means instability (as N17). **Rigid lids:** $w=0$ at $z=0$ and $z=d$ means
$\partial\psi/\partial x=\mathrm ik\hat\psi e^{\mathrm ik(x-ct)}=0$ there, so $\hat\psi(0)=\hat\psi(d)=0$ (11.62); for an unbounded
layer, $\hat\psi\to0$ far away.""")
P("P268", "quadratic eigenvalue problem c²M₂ + cM₁ + M₀ and its companion linearisation", r"""
Multiplying the Taylor–Goldstein equation by $(U-c)$ leaves $c^2$ in it, so collocation gives $(c^2M_2+cM_1+M_0)\mathbf v=0$.
Introduce $\mathbf w=c\mathbf v$; then
$\begin{pmatrix}0&I\\-M_0&-M_1\end{pmatrix}\begin{pmatrix}\mathbf v\\\mathbf w\end{pmatrix}=c\begin{pmatrix}I&0\\0&M_2\end{pmatrix}
\begin{pmatrix}\mathbf v\\\mathbf w\end{pmatrix}$ — an ordinary generalised eigenproblem twice as large (P259). For 1 × 1
"matrices" it is the quadratic formula in disguise.""", code=r"""
import numpy as np                                    # arrays
from scipy.linalg import eig                          # generalised eigen-solver (P259)
M2, M1, M0 = 1.0, -3.0, 4.0                           # c^2 - 3c + 4 = 0 (the quadratic of C02's tiny example, divided by 4)
A = np.array([[0, 1], [-M0, -M1]])                    # companion matrix: second row holds -M0, -M1
B = np.array([[1, 0], [0, M2]])                       # second row holds M2
print(eig(A, B, right=False))                         # 1.5 + 1.3229j and 1.5 - 1.3229j: the same two roots
""")
P("P270", "mapping an infinite domain to [−1, 1]", r"""
A tanh shear layer or a jet extends to $z=\pm\infty$, but Chebyshev points live on $[-1,1]$. A map stretches them: our solvers use
$z=s\tan(\theta\xi)$ with $\theta=\arctan(z_{\max}/s)$, which sends $\xi\in[-1,1]$ to $[-z_{\max},z_{\max}]$, crowds the points
near the layer ($\lvert z\rvert\lesssim s$) and still reaches far away, where the disturbance (which decays like
$e^{-k\lvert z\rvert}$) is set to zero. Derivatives follow from the chain rule (Ch. 1 P49):
$\frac{d}{dz}=\frac{d\xi}{dz}\frac{d}{d\xi}$. Always check that the answer does not change when $s$, $z_{\max}$ or $N$ changes — long
waves need a bigger box (`ST.decay_box(k)` = max(30, 12/k)).""", code=r"""
import numpy as np                                    # arrays
from fluidpy.core import stability as ST              # the eigen-solver toolkit
g = ST.cheb_grid(8, map="tan", y_max=30.0, s=0.5)     # nine mapped Chebyshev points on [-30, 30] with scale s = 0.5
print(np.round(g.y, 2))                               # dense near 0 (the layer), out to +-30
""")
nb.worked_example("two small cases by hand", r"""
**(a) Couette flow, no stratification:** $U=z$ on $[0,1]$, $U''=0$, $N^2=0$. Then the equation is $(U-c)(\hat\psi''-k^2\hat\psi)=0$.
If $c$ is not between 0 and 1, $U-c\neq0$, so $\hat\psi''=k^2\hat\psi$ with $\hat\psi(0)=\hat\psi(1)=0$, whose only solution is
$\hat\psi=0$ (P266): no eigenvalue outside the velocity range, no instability.

**(b) The companion trick:** for the quadratic $c^2-3c+4=0$ (C02's tiny example with $U_1=6$), the 2 × 2 companion matrix
$\begin{pmatrix}0&1\\-4&3\end{pmatrix}$ has eigenvalues $1.5\pm1.3229\,\mathrm i$ — the same roots (P268).""")
code(r"""
prof = ch11.richardson_profiles("tanh", J=0.1)                      # U = tanh z, N^2 = J sech^2 z with J = 0.1 (functions U, Up, Upp, N2, Ri)
c_tg = ch11.taylor_goldstein_eigs(0.4, prof["U"], prof["Upp"], prof["N2"], bc="decay", N=100)   # converged unstable eigenvalues at k = 0.4
raw = ch11.taylor_goldstein_eigs(0.4, prof["U"], prof["Upp"], prof["N2"], bc="decay", N=100, unstable_only=False, filter=False)   # everything
print("unstable eigenvalue c =", np.round(c_tg, 4), "-> growth rate k c_i =", round(0.4*c_tg[0].imag, 4))
print("distance from its conjugate to the nearest raw eigenvalue:", f"{np.min(np.abs(raw - np.conj(c_tg[0]))):.1e}")   # the decaying twin (N73)
d2U = lambda z: -2*np.tanh(z)/np.cosh(z)**2                         # U'' for U = tanh z
c_ray = ch11.rayleigh_eigs(0.4, np.tanh, d2U, bc="decay", N=100)    # Rayleigh's equation (no stratification) for the same layer
c_tg0 = ch11.taylor_goldstein_eigs(0.4, np.tanh, d2U, lambda z: 0*z, bc="decay", N=100)   # Taylor-Goldstein with N^2 = 0
print("N^2 = 0: Taylor-Goldstein", np.round(c_tg0[0], 5), "| Rayleigh", np.round(c_ray[0], 5))
""", explain=r"""
1. `ch11.richardson_profiles("tanh", J)` returns the profile functions of a tanh shear layer with stratification
   $N^2=J\,\mathrm{sech}^2z$ (C09's family).
2. `ch11.taylor_goldstein_eigs(k, U, Upp, N2, bc="decay")` solves the Taylor–Goldstein problem on the mapped infinite domain
   (P270) as a quadratic eigenproblem (P268). By default it returns only unstable eigenvalues that survive an $N$-convergence test;
   `unstable_only=False, filter=False` returns the raw spectrum.
3. The unstable $c\approx0.311\,\mathrm i$ has a conjugate twin in the raw spectrum.
4. With $N^2=0$ the solver reproduces the Rayleigh solver: $c\approx0.470\,\mathrm i$.""")
scratch(r"""
ss = ch11.stratified_shear_sympy()                                  # D16-D17 with sympy: (11.55)-(11.57) -> (11.58)-(11.60) -> (11.61)
print("normal-mode equations (11.58)-(11.60), residuals:", ss["normal_modes"])     # [0, 0, 0]
print("Taylor-Goldstein residual after eliminating pressure and density:", ss["tg_residual"])   # 0
assert ss["tg_residual"] == 0                                       # the elimination of D17 is exact
bad = ch11.stratified_shear_sympy(printed=True)                     # the same elimination starting from the w-equation AS PRINTED (slip #4)
assert bad["printed_slip4_reaches_tg"] is False                     # with dp/dx in the w-equation the pressure does not cancel
print("with slip #4 built in, Taylor-Goldstein is reached?", bad["printed_slip4_reaches_tg"])
""", r"""
1. Here the independent check is symbolic rather than a second numerical version: `ch11.stratified_shear_sympy()` repeats D16–D17
   on generic profiles; the residuals are exactly zero.
2. With the printed pressure term of slip #4, eliminating the pressure does *not* give the Taylor–Goldstein equation.""")
fig(r"""
raw2 = ch11.taylor_goldstein_eigs(0.4, prof["U"], prof["Upp"], prof["N2"], bc="decay", N=120 if FAST else 140, unstable_only=False, filter=False)
arc_r, arc_i = ST.howard_semicircle(-1.0, 1.0)                     # Howard's semicircle for U between -1 and 1 (C10)
fig, ax = plt.subplots(figsize=(7.0, 3.6))
ax.plot(arc_r, arc_i, ":", color=C_BASE, label="Howard's semicircle (C10)")   # the arc that must contain every unstable eigenvalue
ax.plot(arc_r, -arc_i, ":", color=C_BASE)   # its mirror image, for the decaying twins
ax.plot(raw.real, raw.imag, "o", ms=3, color=C_BASE, alpha=0.6, label="raw spectrum, $N$ = 100")   # every eigenvalue the matrices give at N = 100
ax.plot(raw2.real, raw2.imag, "x", ms=4, color=C_BUOY, alpha=0.6, label="raw spectrum, larger $N$")   # ... and at the larger N: what coincides is physics
ax.plot(c_tg[0].real, c_tg[0].imag, "o", ms=9, mfc="none", color=C_GROW, label="unstable mode (converged)")   # circle the growing mode
ax.plot(c_tg[0].real, -c_tg[0].imag, "o", ms=9, mfc="none", color=C_DECAY, label="its decaying twin")   # circle its complex conjugate (N73)
ax.set(xlabel="$c_r$ [$U_0$]", ylabel="$c_i$ [$U_0$]", xlim=(-1.3, 1.3), ylim=(-0.45, 0.45))
ax.set_title("A Taylor–Goldstein spectrum ($k$ = 0.4, $J$ = 0.1)", fontsize=10)
ax.legend(fontsize=7, loc="upper right")
plt.show()
p03 = ch11.richardson_profiles("tanh", J=0.3)                      # more stratification
print("J = 0.3: converged unstable eigenvalues:", ch11.taylor_goldstein_eigs(0.4, p03["U"], p03["Upp"], p03["N2"], bc="decay", N=100))
""",
    see="The complex $c$-plane. Two mirror dots off the real axis (circled rose and teal), a dense line of dots on the real axis "
        "between −1 and 1, and a cloud of dots just above and below that line. Grey circles and blue crosses are the same "
        "calculation at two resolutions. The dotted arcs at the left and right edges are parts of Howard's semicircle (C10), which "
        "must contain every unstable eigenvalue.",
    read="The circled pair is the growing mode and its decaying twin: both resolutions give the same point. The line on the axis "
         "is the *continuous spectrum* — neutral disturbances that each travel with the fluid at one height, where $U(z)=c$ "
         "(P269); they are not instabilities. The cloud just off the axis moves when $N$ changes: numerical junk from that "
         "singular line, which the solver's convergence filter removes.",
    change="…$J=0.3$: no converged unstable eigenvalue is left (the last line prints an empty array) — C09 proves that it must be so.",
    explain=r"""
1. The raw spectra at two resolutions are overlaid; `ST.howard_semicircle(Umin, Umax)` returns the arc of C10 as a ghost.
2. Only what coincides at both resolutions is physics.""")
whatif(r"""
…we wanted a *guarantee* of stability without computing any spectrum? Multiply the equation by the conjugate, integrate, and look
at the imaginary part — C09 turns the Taylor–Goldstein equation into the Richardson-number theorem.""")

# ---------------------------------------------------------------------------------------------------------------------
nb.recap("R20", "The gradient Richardson number", r"""
$\mathrm{Ri}(z)\equiv N^2/(dU/dz)^2$ (11.66): stratification's restoring strength over the shear's — Ch. 4's gradient Richardson
number and Ch. 1's $N^2$ together. Large Ri: buoyancy wins; small Ri: shear wins.""", where="Ch. 4 §4.11, Ch. 1 §1.10")
core("C09", r"The Miles–Howard criterion (11.67): $\mathrm{Ri}=\frac{N^2}{(dU/dz)^2}>\tfrac14$ everywhere ⇒ stable",
     "How much stratification guarantees that a shear layer cannot break into billows?")
problem(r"""
Every ocean and climate model has to decide, grid box by grid box, whether the shear is strong enough to mix the water or air. The
rule they use is a number: if the gradient Richardson number is above about ¼, no mixing from shear instability. Where does ¼ come
from — and is it a guarantee, a trigger, or both?""")
idea(r"""
TG equation  ──substitute φ = ψ̂/(U − c)^{1/2}──►  self-adjoint form (11.64)
             ──× φ*, integrate, by parts──────►  (11.65): ∫ (N² − ¼U′²)/(U − c) |φ|² = ∫(U − c)(|φ′|² + k²|φ|²) + ∫ ½U″|φ|²
             ──imaginary part──────────────────►  c_i · ∫ (N² − ¼U′²)/|U − c|² |φ|²  =  − c_i · ∫ (|φ′|² + k²|φ|²)
if N² > ¼U′² everywhere: c_i × (positive) = c_i × (negative)  ⇒  c_i = 0   ⇒  Ri > ¼ everywhere guarantees stability
""", "*Self-adjoint* means an equation of the shape $(p\\phi')'-q\\phi=0$: multiplied by $\\phi^*$ and integrated, its derivative term turns by parts into $-\\int p\\lvert\\phi'\\rvert^2dz$ with no boundary term when $\\phi$ vanishes at both ends — the only property D18 uses. The two numbered stations of the sketch, both derived in D18 below: the self-adjoint form " + E("11.64") + " and the integral identity " + E("11.65") + ".")
remind("C09")
P("P271", "integrating on a Chebyshev grid (Clenshaw–Curtis weights)", r"""
On Chebyshev points, the integral of the polynomial through the samples is a weighted sum $\int f\,dx\approx\sum_jw_jf(x_j)$; the
weights are exact for polynomials up to degree $N$, so smooth integrands converge as fast as the collocation itself. We use them
to check the integral identities of this section and the energy budget of C14 on computed modes (Simpson's rule, Ch. 9 P203, would
be far less accurate on these clustered points).""", code=r"""
import numpy as np                                    # arrays
from fluidpy.core import stability as ST              # the eigen-solver toolkit
D, x = ST.cheb(8)                                     # 9 Chebyshev points on [-1, 1]
w = ST.clenshaw_curtis_weights(8)                     # their integration weights
print(w @ x**2, w @ np.cos(x))                        # 0.666667 (= 2/3) and 1.682942 (= 2 sin 1): both exact to round-off
""")
confusion(r"""
**$\phi$ is not a potential here.** In D18 the symbol $\phi\equiv\hat\psi/(U-c)^{1/2}$ (11.63) is just a rescaled stream-function
amplitude — a different object from C02's velocity potentials $\phi_1,\phi_2$ and from the $\phi(y)$ of C11–C13.""")
D("D18", ref="11.67")
nb.worked_example("Ri for a thermocline and for the tanh layer", r"""
**Thermocline:** a current changing by $\Delta U=0.1$ m/s over $h=2$ m, so $dU/dz\approx0.05$ s⁻¹, with $N^2=2\times10^{-4}$ s⁻².
$\mathrm{Ri}=2\times10^{-4}/0.05^2=2\times10^{-4}/0.0025=0.08<\tfrac14$: instability is *allowed* (not guaranteed).

**tanh layer:** $U=\tanh z$, $N^2=J\,\mathrm{sech}^2z$. $U'=\mathrm{sech}^2z$, so
$\mathrm{Ri}(z)=J\,\mathrm{sech}^2z/\mathrm{sech}^4z=J\cosh^2z$: smallest at the centre, $\mathrm{Ri}_{\min}=J$. $J=0.1$:
$\mathrm{Ri}_{\min}=0.1<\tfrac14$ (the computed growth is positive, C08); $J=0.3$: $\mathrm{Ri}>\tfrac14$ everywhere ⇒ stable by the
theorem, whatever $k$.""")
code(r"""
z = np.linspace(-4, 4, 801)                                         # heights across the layer [L]
for J in (0.1, 0.24, 0.3):                                          # three strengths of stratification
    p = ch11.richardson_profiles("tanh", J=J)                       # U = tanh z, N^2 = J sech^2 z
    mh = ch11.miles_howard_stable(z, U=p["U"], N2=p["N2"])          # (11.66)-(11.67) on this profile
    print(f"J = {J}: Ri_min = {mh['Ri_min']:.3f} at z = {mh['z_min']:.1f} -> guaranteed stable: {mh['guaranteed_stable']}")
print("growth rate k c_i at k = 0.4 for J = 0, 0.1, 0.3:", [round(ch11.tg_growth(0.4, J), 4) for J in (0.0, 0.1, 0.3)])
p = ch11.richardson_profiles("tanh", J=0.1)                         # back to J = 0.1 for the identity check
md = ch11.taylor_goldstein_eigs(0.4, p["U"], p["Upp"], p["N2"], bc="decay", N=100, map_scale=0.5, return_vectors=True)   # mode + its grid
idn = ch11.richardson_identity_check(0.4, md["c"][0], md["psi"][:, 0], md["y"], p["U"], p["Up"], p["Upp"], p["N2"], grid=md["grid"])
print(f"(11.65): left {idn['lhs']:.6f}, right {idn['rhs']:.6f}, relative residual {idn['residual']:.1e}")
print(f"imaginary parts: left {idn['imag_lhs']:.6e}, right {idn['imag_rhs']:.6e}")
""", explain=r"""
1. `ch11.miles_howard_stable(z, U, N2)` evaluates $\mathrm{Ri}(z)$ on the profile, finds its minimum and applies the theorem.
   $J=0.1$ and $0.24$: not guaranteed; $J=0.3$: guaranteed stable.
2. `ch11.tg_growth(k, J)` returns the leading growth rate $kc_i$ (0 if no converged unstable mode): about 0.188, 0.124 and 0 —
   the computation agrees with the theorem.
3. `return_vectors=True` returns the mode $\hat\psi$ on its mapped grid; `ch11.richardson_identity_check` evaluates both sides of
   D18's integral identity with Clenshaw–Curtis weights (P271): they agree to about 10⁻¹¹, imaginary parts included.""")
scratch(r"""
J = 0.1                                                             # the stratification strength
Ri_mine = (J/np.cosh(z)**2)/(1/np.cosh(z)**2)**2                    # (11.66) typed out: N^2/(U')^2 with N^2 = J sech^2 z, U' = sech^2 z
i = int(np.argmin(Ri_mine))                                         # where it is smallest
mh = ch11.miles_howard_stable(z, U=p["U"], N2=p["N2"])              # the library on the same profile
assert np.isclose(Ri_mine[i], mh["Ri_min"]) and np.isclose(z[i], 0.0)   # same minimum, at the centre of the layer
Ri_thermo = 2e-4/0.05**2                                            # the thermocline of the worked example
print(f"tanh layer: Ri_min = {Ri_mine[i]:.3f} at z = {z[i]:.1f};  thermocline: Ri = {Ri_thermo:.2f}")
""", r"""
1. The Richardson number of the tanh layer typed out, $J\cosh^2z$; its minimum equals the library's.
2. The thermocline of the worked example: 0.08.""")
fig(r"""
gm = ch11.tg_growth_map()                                          # k c_i on a (k, J) grid, read from reference/ch11/tg_growth_map.csv (ours)
zz = np.linspace(-2.5, 2.5, 300)                                   # heights [L]
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.8))
a.plot(np.tanh(zz), zz, color=C_SHEAR, label="$U=\\tanh z$")
for J, ls in ((0.1, "-"), (0.3, "--")):
    pj = ch11.richardson_profiles("tanh", J=J)                     # the profile for this J
    a.plot(pj["N2"](zz), zz, ls, color=C_BUOY, label=f"$N^2$, $J$ = {J}")
    Ri = pj["Ri"](zz)                                              # the Richardson number J cosh^2 z
    a.plot(np.where(Ri <= 1.2, Ri, np.nan), zz, ls, color=C_NEUT, label=f"Ri, $J$ = {J}")   # only the part inside the panel
a.axvline(0.25, color=C_BASE, ls=":")
a.text(0.27, -2.3, "¼", color=C_BASE)
a.set(xlabel="$U$ [$U_0$], $N^2$ [$U_0^2/L^2$], Ri [–]", ylabel="$z$ [$L$]", xlim=(-1.1, 1.25))
a.set_title("(a) profiles and their Richardson number", fontsize=10)
a.legend(fontsize=7, loc="upper left")
cf = b.contourf(gm["k"], gm["J"], gm["kci"], levels=np.linspace(0.0, 0.2, 11), cmap="RdPu")
fig.colorbar(cf, ax=b, label="growth rate $kc_i$ [$U_0/L$]")
kk = np.linspace(0.0, 1.0, 200)                                    # wavenumbers for the exact neutral curve
b.plot(kk, ch11.tg_tanh_neutral_J(kk), color=C_NEUT, label="exact neutral curve $J=k(1-k)$")
b.axhline(0.25, color=COLORS["ink"], ls="--", lw=1, label="$J$ = ¼")
b.set(xlabel="wavenumber $k$ [$1/L$]", ylabel="$J=\\mathrm{Ri}_{\\min}$ [–]", xlim=(0.05, 1.0), ylim=(0, 0.3))
b.set_title("(b) growth exists only below ¼", fontsize=10)
b.legend(fontsize=7, loc="upper right")
savefig(fig, "ch11", "c09_richardson")
plt.show()
live_pt = ch11.tg_growth(0.4, 0.1)                                 # one point of the map recomputed now
assert np.isclose(live_pt, gm["kci"][10, 7], rtol=6e-4)            # the stored table has 4 significant figures
print(f"largest growth in the map: {gm['kci'].max():.4f} at J = 0;  map point (k = 0.4, J = 0.1) = {gm['kci'][10, 7]}, recomputed {live_pt:.5f}")
""",
    see="(a) The shear layer $U=\\tanh z$ (orange), the stratification $N^2$ (blue) and the Richardson number (purple) for a weak "
        "(solid) and a stronger (dashed) stratification; the dotted line is ¼. (b) The growth rate $kc_i$ of the leading mode over "
        "the plane of wavenumber and stratification: a tongue that closes on the purple curve, whose top just touches the dashed "
        "line $J=\\tfrac14$ at $k=\\tfrac12$.",
    read="For $J=0.1$ the purple curve in (a) dips below ¼ in the middle of the layer — growth is allowed, and (b) shows it happens. "
         "For $J=0.3$ the curve stays right of ¼ everywhere — nothing can grow, the guarantee. For this family the neutral curve is "
         "known exactly, $J=k(1-k)$. One caution about the map: a zero within about 0.006 (in $J$) of that curve means "
         "'$kc_i<0.004$, too weak for the solver's convergence test', not 'stable'.",
    change="…the density layer were thinner than the shear layer: the minimum of Ri would sit off-centre and the tongue would change "
           "shape (the explainer's profile options).",
    explain=r"""
1. `ch11.tg_growth_map()` reads our stored table (20 wavenumbers × 31 values of $J$; recomputing it takes about two minutes);
   the last lines recompute one point live and compare.
2. `ch11.tg_tanh_neutral_J(k)` is the exact neutral curve of this family, $J=k(1-k)$.
3. `contourf` fills the growth rate; the largest value, 0.19 at $J=0$ near $k=0.45$, is the unstratified shear layer of C12.""")
note("N78 [B]", r"""
**Necessary, not sufficient.** $\mathrm{Ri}<\tfrac14$ somewhere is needed for instability but does not force it; there is no
universal critical Ri (walls and profile shapes lower it). For common shear layers (linear, tanh, erf) instability does start when
the minimum drops below ¼, and as it goes to zero the fastest wave has $kL\approx0.445$ for $U=U_0\tanh(z/L)$, i.e. a wavelength
of about $14L$ (Michalke 1964; our solver gives 0.4449 with growth rate 0.1897). Laboratory (Scotti & Corcos 1972) and ocean
(Eriksen 1978) data support "minimum Ri below ¼" as a useful guide.""")
nb.md(r"""
**Slider figure F5 — the growth curve shrinks as the stratification rises.** One row of the stored map per slider position.""")
nb.plotly(r"""
def f5(J):                                                         # traces for one stratification strength
    row = gm["kci"][int(np.argmin(np.abs(gm["J"] - J)))]           # the stored growth rates at this J
    disc = 1 - 4*J                                                 # J = k(1 - k) has real roots k only if J <= 1/4
    kn = 0.5 + 0.5*np.array([-1, 1])*np.sqrt(disc) if disc >= 0 else np.array([np.nan, np.nan])   # the two neutral wavenumbers
    return {"growth rate k c_i (this J)": (gm["k"], row), "no stratification (J = 0)": (gm["k"], gm["kci"][0]),   # this J, and J = 0 for reference
            "neutral wavenumbers, J = k(1 − k)": (kn, 0*kn)}                                             # two markers on the axis

figF5 = slider_figure(f5, "J = Ri_min", gm["J"][::2] if FAST else gm["J"][:27], xlabel="wavenumber k [1/L]", ylabel="growth rate k c_i [U0/L]",   # slider: J
                      title="Stratification closes the band of growing waves at J = 1/4", yrange=(-0.01, 0.2),   # title and a fixed vertical range
                      modes={"neutral wavenumbers, J = k(1 − k)": "markers"})   # the neutral points drawn as markers
figF5.show()
""", explain=r"""
1. `f5(J)` picks the map row for that $J$ and marks the two exact neutral wavenumbers.
2. Drag $J$ up: the hump sinks and narrows, the two markers move together and meet at $k=\tfrac12$ when $J=\tfrac14$ — beyond
   that nothing grows.""")
explainer("richardson_shear_instability", "Why does Ri = ¼ decide whether a shear layer billows?", r"""
Sliding the stratification shows the Richardson-number profile dip below the ¼ line at the same moment the growth tongue of the
(k, J) map is reached and the billow animation starts growing — the necessary-not-sufficient logic is visible.""", [
    "Preset 'J = 0.24' vs 'J = 0.26': read the status before and after ¼.",
    "Drag k at J = 0.1: growth only in a band of wavenumbers — Ri < ¼ allows, the wavelength decides.",
    "Click the Ri(z) profile: the inspector computes N²/U′² at that height.",
])
whatif(r"""
…we asked not *whether* a wave can grow but *how fast and how fast it travels*? A second substitution, $F=\hat\psi/(U-c)$, bounds
both at once — C10.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C10", r"Howard's semicircle $\big[c_r-\tfrac12(U_{\max}+U_{\min})\big]^2+c_i^2\le\big[\tfrac12(U_{\max}-U_{\min})\big]^2$ and the growth bound $kc_i\le\tfrac k2(U_{\max}-U_{\min})$",
     "Where in the complex plane can the wave speed of an unstable mode be?")
problem(r"""
Before computing a spectrum it helps to know where to look — and how fast anything could possibly grow. In a current whose speed
varies between 0.2 and 1 m/s, can an unstable wave travel at 3 m/s? Can it grow in a microsecond? Howard's answer: no — every
unstable wave speed lies in a half-disc spanned by the slowest and fastest fluid.""")
idea(r"""
 c_i ▲            ___
     │        .-'     '-.          every unstable c lies under this arc
     │      /             \        radius = ½(U_max − U_min),  centre = ½(U_max + U_min)
     │     |       ●       |       ⇒ U_min < c_r < U_max   and   k c_i ≤ (k/2)(U_max − U_min)
     └─────┴───────────────┴────► c_r
         U_min            U_max
""", r"""
(The semicircle itself carries no number in the book; it follows from the numbered inequality
$\int[U^2-c_r^2-c_i^2]Q\,dz>0$ (11.72), step 11 of D19.)""")
remind("C10")
P("P272", "completing the square into a circle (x − a)² + y² ≤ R²", r"""
An inequality $x^2+y^2-2ax+b\le0$ is a disc: add and subtract $a^2$, $(x-a)^2+y^2\le a^2-b$ — centre $(a,0)$, radius
$\sqrt{a^2-b}$ (Ch. 4 completed squares for tensors, P129). Howard's last step is exactly this with $x=c_r$, $y=c_i$.""", code=r"""
import numpy as np                                    # arrays
a, b = 1.0, 0.0                                       # x^2 + y^2 - 2x <= 0
print(np.sqrt(a**2 - b))                              # radius 1 around (1, 0): the half-disc on [0, 2] for y >= 0
""")
D("D19", ref="11.72", after=r"""
> ⚠️ **slip #11 — the book prints** the divergence form as $\frac d{dz}\big[(U-c)^2F'\big]-k^2(U-c)F+N^2F=0$ **; the correct form
> is** $\frac d{dz}\big[(U-c)^2F'\big]-k^2(U-c)^2F+N^2F=0$ — expanding the line above it gives $(U-c)^2$ (step 5), and the integral
> that follows (with $(U-c)^2\lvert F\rvert^2$) confirms it.

> ⚠️ **slip #8 — the book prints** the velocity-range inequality without a weight, $\int[U_{\min}-U][U_{\max}-U]\,dz\le0$, **; the
> correct form is** $\int(U-U_{\min})(U-U_{\max})\,Q\,dz\le0$ — the positive weight $Q=\lvert F'\rvert^2+k^2\lvert F\rvert^2$ must be
> inside the integral from the start, otherwise it cannot be combined with the two identities (step 12).""")
nb.worked_example("a semicircle with easy numbers", r"""
$U$ between $U_{\min}=0$ and $U_{\max}=2$ m/s.

1. Centre $\tfrac12(0+2)=1$, radius $\tfrac12(2-0)=1$.
2. Is $c=1.2+0.5\,\mathrm i$ possible? $(1.2-1)^2+0.5^2=0.04+0.25=0.29\le1$ ✓.
3. $c=2.5+0.1\,\mathrm i$? $(1.5)^2+0.01=2.26>1$ ✗ — no unstable wave can outrun the fastest fluid.
4. Growth bound for $k=2$ m⁻¹: $kc_i\le(2/2)(2-0)=2$ s⁻¹, so the e-folding time is at least 0.5 s.
5. tanh layer ($U$ from −1 to 1): its fastest mode $c\approx0.426\,\mathrm i$ is well inside the unit half-disc.""")
code(r"""
modes = []                                                          # (k, J, c) of every unstable mode found
for J in (0.0, 0.1):                                                # no stratification, and J = 0.1
    pj = ch11.richardson_profiles("tanh", J=J)                      # the tanh layer
    for k in np.arange(0.1, 0.95, 0.2 if FAST else 0.1):            # a set of wavenumbers
        e = ch11.taylor_goldstein_eigs(k, pj["U"], pj["Upp"], pj["N2"], bc="decay", N=100, y_max=ST.decay_box(k), map_scale=ST.decay_map_scale(k))
        if len(e):                                                  # an unstable, converged eigenvalue exists
            modes.append((k, J, e[0]))
cs = np.array([m[2] for m in modes])                                # the unstable wave speeds
print(len(cs), "unstable modes; all inside Howard's semicircle on [-1, 1]:", all(bool(ST.in_howard_semicircle(c, -1.0, 1.0)) for c in cs))
hw = ch11.howard_identity_check(0.4, md["c"][0], md["psi"][:, 0], md["y"], p["U"], p["N2"], grid=md["grid"], Up=p["Up"])   # C09's mode
print({n: (f"{v:.1e}" if isinstance(v, float) else v) for n, v in hw.items()})
""", explain=r"""
1. For each $(k,J)$ the unstable Taylor–Goldstein eigenvalue (if any), computed with the box and node rule that long and
   near-neutral waves need (`ST.decay_box`, `ST.decay_map_scale`).
2. `ST.in_howard_semicircle(c, Umin, Umax)` tests the semicircle inequality; every computed mode passes.
3. `ch11.howard_identity_check` evaluates D19's two integral identities for the mode of C09: residuals `res_69`, `res_70` of
   order 10⁻¹¹; `cr_mean` $=\int UQ/\int Q$ equals $c_r$ (zero for this symmetric layer); `int_11_72` is positive.""")
scratch(r"""
inside = (cs.real - 0.0)**2 + cs.imag**2 <= 1.0 + 1e-9             # the semicircle typed out: centre 0, radius 1 for U from -1 to 1
assert np.all(inside) and np.all(inside == np.array([bool(ST.in_howard_semicircle(c, -1.0, 1.0)) for c in cs]))
ks_m = np.array([m[0] for m in modes])                              # their wavenumbers
assert np.all(ks_m*cs.imag <= ks_m*(1.0 - (-1.0))/2)                # the growth bound k c_i <= (k/2)(U_max - U_min)
print("all", len(cs), "modes satisfy the inequality typed by hand; largest k c_i / bound =", round(np.max(cs.imag), 3))
""", r"""
1. The semicircle inequality written out for this layer and compared with the library function.
2. The growth bound: every $kc_i$ is below $k$; the ratio (printed) is at most about 0.84 — the bound is safe but not tight.""")
fig(r"""
fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, 3.6))
a.plot(arc_r, arc_i, color=C_NEUT, label="Howard's semicircle")
a.plot([-1, 1], [0, 0], color=C_SHEAR, lw=4, alpha=0.6, label="velocity range")
for J, mk in ((0.0, "o"), (0.1, "s")):
    sel = np.array([m[1] == J for m in modes])                     # the modes with this J
    a.plot(cs[sel].real, cs[sel].imag, mk, color=C_GROW, ms=5, label=f"unstable modes, $J$ = {J}")
    b.plot(ks_m[sel], ks_m[sel]*cs[sel].imag, mk + "-", color=C_GROW, ms=4, label=f"$kc_i$, $J$ = {J}")
a.set(xlabel="$c_r$ [$U_0$]", ylabel="$c_i$ [$U_0$]", xlim=(-1.15, 1.15), ylim=(0, 1.1))
a.set_aspect("equal")
a.set_title("(a) every unstable $c$ is under the arc", fontsize=10)
a.legend(fontsize=7, loc="upper right")
b.plot(kk, kk*(1.0 - (-1.0))/2, "--", color=C_BASE, label="bound $\\frac{k}{2}(U_{\\max}-U_{\\min})$")
b.set(xlabel="wavenumber $k$ [$1/L$]", ylabel="growth rate $kc_i$ [$U_0/L$]", xlim=(0, 1), ylim=(0, 0.6))
b.set_title("(b) growth stays below the bound", fontsize=10)
b.legend(fontsize=7)
plt.show()
""",
    see="(a) The half-disc on the velocity range [−1, 1] (purple arc, orange bar) with the computed unstable wave speeds as rose "
        "markers on the imaginary axis. (b) Their growth rates against k, far below the dashed bound.",
    read="The symmetric layer's waves do not travel ($c_r=0$: no preferred direction, as in N19), and they sit well inside the arc. "
         "The arc is far from tight — a sanity window that tells a solver where to look and how fast anything can grow, not a "
         "prediction.",
    change="…an asymmetric layer $U=1+\\tanh z$: the arc and the markers shift right together by 1 (a change of frame).")
nb.md(r"""
**Slider figure F6 — a bounded shear flow: $U=\sin y$ between walls at $y=\pm b$.** For each half-width $b$: Rayleigh's equation
(C12) solved for four wavenumbers, with Howard's semicircle for that velocity range. The unstable eigenvalues appear only when
$2b>\pi$ (N95 in C12 explains why this profile matters).""")
nb.plotly(r"""
def f6(b):                                                         # traces for one half-width b
    ps = ch11.inviscid_profile("sin", b=b)                         # U = sin y on |y| <= b, with walls
    Umax = 1.0 if b >= np.pi/2 else np.sin(b)                      # the fastest fluid (the slowest is -Umax)
    ar, ai = ST.howard_semicircle(-Umax, Umax, n=60)               # the semicircle for this velocity range
    cu = []                                                        # unstable eigenvalues found
    for k in (0.2, 0.4, 0.6, 0.8):
        e = ch11.rayleigh_eigs(k, ps["U"], ps["Upp"], domain=ps["domain"], N=80, tol=1e-4, unstable_only=True)   # (11.81)-(11.82)
        cu += [e[0]] if len(e) else []
    cu = np.array(cu) if cu else np.array([np.nan + 0j])           # nothing unstable: an empty marker
    return {"Howard's semicircle": (ar, ai), "unstable c (k = 0.2, 0.4, 0.6, 0.8)": (cu.real, cu.imag)}

figF6 = slider_figure(f6, "b", np.round(np.linspace(1.5, 3.0, 6 if FAST else 11), 2), xlabel="c_r", ylabel="c_i",
                      title="U = sin y between walls at ±b: unstable only when 2b > π", xrange=(-1.1, 1.1), yrange=(-0.02, 1.05),
                      modes={"unstable c (k = 0.2, 0.4, 0.6, 0.8)": "markers"})
figF6.show()
""", explain=r"""
1. `ch11.inviscid_profile("sin", b=b)` returns the profile functions and the domain; `ch11.rayleigh_eigs(..., unstable_only=True)`
   returns the converged unstable eigenvalues (none for $b=1.5$, where $2b<\pi$).
2. As $b$ grows past $\pi/2\approx1.57$, markers leave the real axis and climb — always under the arc.""")
whatif(r"""
…viscosity came back? The equation gains a fourth derivative, the critical-layer singularity disappears, conjugate pairs break —
and, surprisingly, a profile with no inflection point can become unstable. That is the Orr–Sommerfeld story, C11–C14.""")
# =====================================================================================================================
# A.8 §11.8 Squire's Theorem and the Orr–Sommerfeld Equation — R21 R22, C11
# =====================================================================================================================
nb.section("11.8", "Squire's Theorem and the Orr-Sommerfeld Equation", intro=r"""
**What is this section about?** Viscous parallel flows — water in a channel, air in a boundary layer. Two simplifications make
them tractable: Squire's theorem says two-dimensional disturbances go unstable first, and a stream function then turns the
linearised Navier–Stokes equations into one fourth-order equation, the Orr–Sommerfeld equation, whose eigenvalue $c(k,\mathrm{Re})$
decides stability.""")
nb.recap("R21", "Dimensionless Navier–Stokes", r"""
Lengths by $L$, velocities by $U_0$, time by $L/U_0$, pressure by $\rho U_0^2$, $\mathrm{Re}=U_0L/\nu$ (Ch. 4 §4.11): the
perturbed $x$-momentum equation is $\frac{\partial u}{\partial t}+(U+u)\frac{\partial}{\partial x}(U+u)+v\frac{\partial}{\partial y}(U+u)=
-\frac{\partial}{\partial x}(P+p)+\frac1{\mathrm{Re}}\nabla^2(U+u)$ (11.73) — the starting line of D20; $U(y)$ is the basic flow,
$(u,v,w)$ and $p$ the disturbance. ⚠️ This is the third of the chapter's four scalings (the table in the conventions cell).""",
         where="Ch. 4 §4.11")
nb.recap("R22", "The stream function, §11.8 version", r"""
Here $y$ is across the flow: $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$ (Ch. 4 §4.3). With
$\psi=\phi(y)e^{\mathrm ik(x-ct)}$: $\hat u=\phi'$, $\hat v=-\mathrm ik\phi$ — and $\phi$ is now the stream-function amplitude,
**not** a potential (slip #9).""", where="Ch. 4 §4.3")

core("C11", r"The Orr–Sommerfeld equation $(U-c)\big(\frac{d^2\phi}{dy^2}-k^2\phi\big)-\frac{d^2U}{dy^2}\phi=\frac1{\mathrm ik\mathrm{Re}}\big(\frac{d^4\phi}{dy^4}-2k^2\frac{d^2\phi}{dy^2}+k^4\phi\big)$ (11.79)",
     "What equation decides whether a viscous channel or boundary-layer flow is stable?")
problem(r"""
Reynolds showed in 1883 that water in a pipe stays smooth at low speed and turns turbulent at high speed. The wing of an airliner
keeps a laminar boundary layer only near its nose. To predict where smooth flow breaks down we need the viscous version of C08's
equation: small waves on $U(y)$, now with friction. Friction should only damp — but C13 will show that it can do the opposite.""")
idea(r"""
3-D wave e^{i(kx + mz − kct)} at Re   ──Squire──►  2-D wave with k̄ = √(k² + m²) at the LOWER Re̅ = k Re / k̄
     ⇒ the first instability (lowest Re) is two-dimensional
2-D: u = ∂ψ/∂y, v = −∂ψ/∂x, ψ = φ(y)e^{ik(x−ct)}   ──eliminate p──►   one 4th-order ODE for φ:  Orr–Sommerfeld
     left side: Rayleigh's inviscid operator      right side: viscosity, ∝ 1/(ikRe)
""")
remind("C11")
note("N83 [C]", r"""
**Viscosity can destabilise.** In Bénard and Taylor flows viscosity only raised thresholds; in channels and boundary layers it can
be the very cause of instability (C13), for a reason C14 explains.""")
D("D20", ref="11.77")
note("N88 [B]", r"""
**Squire's theorem.** The transformation $\bar k=\sqrt{k^2+m^2}$, $\bar c=c$, $\bar k\bar u=k\hat u+m\hat w$, $\bar v=\hat v$,
$\bar p/\bar k=\hat p/k$, $\bar k\,\overline{\mathrm{Re}}=k\,\mathrm{Re}$ (11.78) turns the three-dimensional normal-mode equations
into the two-dimensional ones with $m=\hat w=0$ (multiply the $x$-equation by $k$, the $z$-equation by $m$, add, divide by
$\bar k$ — stated, not derived). Because $\bar k\ge k$, the equivalent 2-D problem has the lower Reynolds number
$\overline{\mathrm{Re}}=k\,\mathrm{Re}/\bar k\le\mathrm{Re}$: the critical Re is found with 2-D waves. ⚠️ Not valid with rotation or
stratification (Ch. 13); §11.7 assumed it.""")
code(r"""
U = lambda y: 1 - y**2                                              # plane Poiseuille flow in channel units (half-width, centreline speed)
Upp = lambda y: -2 + 0*y                                            # its second derivative U'' = -2 (as an array of the same shape)
sq = ST.squire_transform(1.0, np.tan(np.pi/6), 1e4)                 # an oblique wave: k = 1, m = tan 30 deg, Re = 10^4  ->  (11.78)
print({n: (round(v, 4) if isinstance(v, float) else v) for n, v in sq.items()})
c3 = ch11.os_3d_eigs(1.0, np.tan(np.pi/6), 1e4, U, N=60)[0]         # the 3-D normal-mode equations (11.77) solved directly
c2 = ch11.orr_sommerfeld_eigs(sq["kbar"], sq["Rebar"], U, Upp)[0]   # the equivalent 2-D problem at (kbar, Rebar)
print("leading c, 3-D problem:", np.round(c3, 8), "| equivalent 2-D problem:", np.round(c2, 8), "| difference", f"{abs(c3 - c2):.1e}")
""", explain=r"""
1. `ST.squire_transform(k, m, Re)` applies (11.78): $\bar k=1.1547$, $\overline{\mathrm{Re}}=8660$ (lower by $\cos30°$), angle 30°.
2. `ch11.os_3d_eigs` solves the three-dimensional equations in the primitive variables $\hat u,\hat v,\hat w,\hat p$;
   `ch11.orr_sommerfeld_eigs` the two-dimensional Orr–Sommerfeld problem derived below (D21).
3. The two leading eigenvalues agree to about 10⁻¹⁰ — Squire's theorem checked with numbers. Here $c_i<0$: this oblique wave decays,
   although the 2-D wave with $k=1$ at the same $\mathrm{Re}=10^4$ grows (the main code cell below).""")
note("N89 [B]", r"""
**What Squire means physically.** An oblique wave sees only the component of the basic flow along the normal to its crests — a
slower flow, a lower effective Reynolds number — and the growth rate of the equivalent 2-D wave, $\bar k\bar c_i$, exceeds $kc_i$.
So two-dimensional disturbances are the most dangerous (for parallel, non-rotating, unstratified flows).""")
D("D21", ref="11.79", check_src=r"""
os_ok = ch11.os_derivation_sympy()                    # D21 with sympy: stream function in, pressure eliminated, compared with (11.79)
print("residual against the Orr-Sommerfeld equation:", os_ok["residual"])            # 0
os_bad = ch11.os_derivation_sympy(printed_v_sign=True)   # the same with the WRONG sign v-hat = +ik phi planted
print("with v-hat = +ik phi the residual is zero?", os_bad["residual"] == 0)         # False: the planted sign error is caught
for i, s in enumerate(os_ok["steps"], 1):             # the engine's own list of moves
    print(f"  {i}. {s}")
""")
note("N90 [B]", r"""
**No slip:** $\phi=\frac{d\phi}{dy}=0$ at $y=y_1$ and $y_2$ (11.80) — $\hat v=-\mathrm ik\phi=0$ and $\hat u=\phi'=0$ on each wall.
Four conditions for a fourth-order equation; the eigenvalue is $c$.""")
nb.worked_example("a channel and an oblique wave", r"""
1. Water ($\nu=10^{-6}$ m²/s) in a channel of half-width $L=1$ cm with centreline speed $U_0=0.58$ m/s:
   $\mathrm{Re}=U_0L/\nu=0.58\times0.01/10^{-6}=5800$ — right at the edge of instability (C13: 5772).
2. An oblique wave at 30° to the flow, $k=1$, $m=\tan30°=0.577$ (in units of $1/L$): $\bar k=\sqrt{1+0.333}=1.1547$.
3. Its 2-D twin sees $\overline{\mathrm{Re}}=\mathrm{Re}\cdot k/\bar k=5800\times0.866=5023<5772$: the oblique wave is stable even
   though the 2-D wave at the same Re is (just) not.""")
code(r"""
c_os = ch11.orr_sommerfeld_eigs(1.0, 1e4, U, Upp, N=100)            # Orr-Sommerfeld eigenvalues for Poiseuille flow, k = 1, Re = 10^4
print("leading eigenvalue c =", f"{c_os[0]:.8f}", "-> growth rate k c_i =", f"{c_os[0].imag:.5f} (in units U0/L)")   # c_i > 0: this wave grows
for N in (60, 80, 100, 120):                                        # the same eigenvalue at four resolutions
    print(f"  N = {N}: c = {ch11.orr_sommerfeld_eigs(1.0, 1e4, U, Upp, N=N)[0]:.8f}")   # should not change with N
ref = BENCH["plane_poiseuille_Re1e4_k1"]                            # the published benchmark, read from reference/ch11 (Orszag 1971)
print(f"published: {ref['c_r']} + {ref['c_i']}i  -> difference {abs(c_os[0] - complex(ref['c_r'], ref['c_i'])):.1e}")   # our value minus Orszag's
print("at Re = 5000:", f"{ch11.orr_sommerfeld_eigs(1.0, 5000.0, U, Upp)[0]:.5f}")   # below the critical Re: c_i < 0, it decays
""", explain=r"""
1. The plane-Poiseuille base flow in channel units: $U=1-y^2$ on $-1\le y\le1$.
2. `ch11.orr_sommerfeld_eigs(k, Re, U, Upp)` solves the Orr–Sommerfeld problem by Chebyshev collocation with the four no-slip rows
   and returns the eigenvalues sorted by $c_i$, largest first. The leading one has $c_i>0$: this wave **grows**.
3. Four resolutions agree to eight digits.
4. The published benchmark (read from a file, Ch. 10 P251) is reproduced — this is *verification* (Ch. 10 P253: are we solving
   the equation right?). At $\mathrm{Re}=5000$ the same mode has $c_i<0$: it decays.""")
scratch(r"""
k_, Re_, N_ = 1.0, 1e4, 100                                         # wavenumber, Reynolds number, number of intervals
Dc, yc = ST.cheb(N_)                                                # Chebyshev matrix and points on [-1, 1] (P258)
I = np.eye(N_ + 1)                                                  # identity
D2 = Dc @ Dc                                                        # second derivative
D4 = D2 @ D2                                                        # fourth derivative
A = np.diag(U(yc)) @ (D2 - k_**2*I) - np.diag(Upp(yc)) - (D4 - 2*k_**2*D2 + k_**4*I)/(1j*k_*Re_)   # (11.79) sorted as A phi = c B phi
B = D2 - k_**2*I                                                    # the operator that multiplies c
A = A.astype(complex); B = B.astype(complex)                        # complex matrices (the viscous term carries an i)
for r, row in ((0, I[0]), (N_, I[N_]), (1, Dc[0]), (N_ - 1, Dc[N_])):   # (11.80): phi = 0 and phi' = 0 at both walls
    A[r, :] = row; B[r, :] = 0.0                                    # row replacement (P259)
lam = linalg.eig(A, B, right=False)                                 # generalised eigenvalues
lam = lam[np.isfinite(lam) & (np.abs(lam) < 2)]                     # drop the infinite ones and the huge spurious ones
c_mine = lam[np.argmax(lam.imag)]                                   # the least stable mode
assert np.isclose(c_mine, c_os[0], atol=1e-8) and np.isclose(c_mine, complex(ref["c_r"], ref["c_i"]), atol=1e-7)
print(f"from scratch: c = {c_mine:.8f}")
""", r"""
1. The Orr–Sommerfeld equation rearranged as $A\phi=cB\phi$ with
   $A=U(D^2-k^2)-U''-\frac1{\mathrm ik\mathrm{Re}}(D^4-2k^2D^2+k^4)$ and $B=D^2-k^2$.
2. Four rows are replaced by the no-slip conditions; `scipy.linalg.eig(A, B)` does the rest.
3. The least stable eigenvalue equals the library's and the published value.""")
fig(r"""
rawA = ch11.orr_sommerfeld_eigs(1.0, 1e4, U, Upp, N=100, filter=False)    # every finite eigenvalue, N = 100
rawB = ch11.orr_sommerfeld_eigs(1.0, 1e4, U, Upp, N=140, filter=False)    # ... and N = 140
fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.plot(rawA.real, rawA.imag, "o", ms=4, color=C_BUOY, label="$N$ = 100")
ax.plot(rawB.real, rawB.imag, "x", ms=5, color=C_BASE, label="$N$ = 140")
ax.plot(c_os[0].real, c_os[0].imag, "o", ms=11, mfc="none", color=C_GROW, label="the one growing mode")
ax.axhline(0, color=COLORS["ink"], lw=0.8)
for (xt, yt, txt) in ((0.12, -0.25, "A (wall modes)"), (0.72, -0.2, "P (centre modes)"), (0.5, -0.85, "S")):
    ax.text(xt, yt, txt, fontsize=8, color=C_BASE)
ax.set(xlabel="$c_r$ [$U_0$]", ylabel="$c_i$ [$U_0$]", xlim=(0, 1.05), ylim=(-1.0, 0.1))
ax.set_title("The Orr–Sommerfeld spectrum is a Y (Re = 10$^4$, $k$ = 1)", fontsize=10)
ax.legend(fontsize=8, loc="lower left")
savefig(fig, "ch11", "c11_os_spectrum")
plt.show()
""",
    see="Eigenvalues in the complex $c$-plane arranged like the letter Y: two arms at the top (left: wave speeds around 0.2–0.4, "
        "right: around 0.9–1) and a stem going straight down near $c_r\\approx0.67$. Circles and crosses are two resolutions. One "
        "circled eigenvalue at the top of the left arm sits just above the axis.",
    read="That weak mode is the Tollmien–Schlichting wave: $c_i=0.0037$, the only growing one; everything else decays. The left arm "
         "holds modes that live near the walls (slow), the right arm modes near the centreline (fast), the stem modes that span "
         "the channel. Where circles and crosses disagree (the lower stem) the eigenvalues are not converged — the stem is "
         "notoriously sensitive.",
    change="…$\\mathrm{Re}=5000$: the circled mode drops below the axis ($c_i=-0.0018$) — C13 finds where it crosses.",
    explain=r"""
1. `filter=False` returns the raw spectrum; two resolutions are overlaid, as for the Taylor–Goldstein spectrum of C08.
2. The labels A, P, S are the customary names of the three branches.""")
whatif(r"""
…$\mathrm{Re}\to\infty$? The right side of the Orr–Sommerfeld equation disappears, the order drops from four to two, two boundary
conditions must be abandoned, and the equation is Rayleigh's — whose theorems need no computer at all (C12).""")

# =====================================================================================================================
# A.9 §11.9 Inviscid Stability of Parallel Flows — C12
# =====================================================================================================================
nb.section("11.9", "Inviscid Stability of Parallel Flows", intro=r"""
**What is this section about?** Drop viscosity from the disturbances (the basic profile may still be a viscous one): the
Orr–Sommerfeld equation becomes Rayleigh's equation. Two theorems then classify profiles by their shape alone — a growing wave
needs an inflection point (Rayleigh) at which the vorticity peaks (Fjørtoft) — and neutral waves come with a "critical layer" and
Kelvin's cat's-eye streamlines.""")
core("C12", r"Rayleigh's inflection-point theorem $c_i\int\frac{1}{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=0$ (11.84): an unstable inviscid parallel flow needs $U''$ to change sign",
     "Can we tell from the shape of a velocity profile whether it can be unstable without viscosity?")
problem(r"""
Jets, wakes and mixing layers break into eddies within a few widths downstream; channel and boundary-layer flows hold out much
longer. A forecaster sees jet streams meander into weather systems; an engineer sees a wake shed vortices. Is there something
about the *shape* of a profile that makes it fragile?""")
idea("", r"""
$U''$ is (minus) the slope of the vorticity $-U'$ across the flow. An inflection point ($U''=0$ with a sign change) is where the
vorticity has an extremum — the smooth version of C02's vortex sheet. Rayleigh's identity says: no such point, no inviscid growth.
Fjørtoft adds that the extremum must be a *maximum* of the vorticity magnitude.""")
remind("C12")
note("N91 [B] · N92 [B]", r"""
**Rayleigh's equation** ($\mathrm{Re}\to\infty$ in the Orr–Sommerfeld equation):
$(U-c)\big(\frac{d^2\phi}{dy^2}-k^2\phi\big)-\frac{d^2U}{dy^2}\phi=0$ (11.81) — C08's Taylor–Goldstein equation with $N^2=0$. The
limit is *singular*: the order drops from 4 to 2 and only two boundary conditions can be kept, $\phi=0$ at $y=y_1,y_2$ (11.82) (no
flow through the walls; slip is now allowed). No $\mathrm i$ appears, so eigenvalues come in conjugate pairs — a property the
viscous term (with its $\mathrm i$) destroys.""")
fig(r"""
d2tanh = lambda y: -2*np.tanh(y)/np.cosh(y)**2                     # U'' for the shear layer U = tanh y
sp_ray = ch11.rayleigh_eigs(0.4, np.tanh, d2tanh, bc="decay", N=100, filter=False)                 # inviscid raw spectrum, k = 0.4
sp_os = {Re: ch11.orr_sommerfeld_eigs(0.4, Re, np.tanh, d2tanh, bc="decay", N=100, filter=False) for Re in (100.0, 1000.0)}   # viscous
fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, 3.5), sharey=True)
a.plot(sp_ray.real, sp_ray.imag, "o", ms=4, color=C_BASE)
a.plot([0, 0], [sp_ray[0].imag, -sp_ray[0].imag], "o", ms=10, mfc="none", color=C_GROW)
a.set(xlabel="$c_r$ [$U_0$]", ylabel="$c_i$ [$U_0$]", xlim=(-1.2, 1.2), ylim=(-0.6, 0.6))
a.set_title("(a) Rayleigh: mirror pairs", fontsize=10)
for Re, mk, col in ((100.0, "s", C_SHEAR), (1000.0, "^", C_BUOY)):
    b.plot(sp_os[Re].real, sp_os[Re].imag, mk, ms=4, color=col, label=f"Re = {Re:.0f}")
b.set(xlabel="$c_r$ [$U_0$]", xlim=(-1.2, 1.2))
b.set_title("(b) Orr–Sommerfeld: no pairs", fontsize=10)
b.legend(fontsize=8)
plt.show()
print("inviscid pair: c =", np.round(sp_ray[0], 4), "and", np.round(sp_ray[-1], 4))
print("viscous leading c: Re = 100:", np.round(sp_os[100.0][0], 4), "| Re = 1000:", np.round(sp_os[1000.0][0], 4))
""",
    see="(a) Inviscid eigenvalues of the tanh shear layer: a line on the real axis and one pair of mirror points above and below "
        "it (circled). (b) The viscous spectra at two Reynolds numbers: one point above the axis, and everything else in the "
        "lower half-plane — no mirror image.",
    read="Without viscosity every growing mode has a decaying twin (the equation has real coefficients). Viscosity breaks the "
         "symmetry: the decaying twin is gone, and the continuous line of real eigenvalues is pushed into the damped half-plane.",
    change="…$\\mathrm{Re}\\to\\infty$: the viscous unstable eigenvalue approaches the inviscid one (0.437 at Re = 100, 0.467 at "
           "Re = 1000, 0.470 inviscid).",
    explain=r"""
1. `ch11.rayleigh_eigs(k, U, Upp, bc="decay")` solves Rayleigh's equation on the mapped infinite domain (P270);
   `ch11.orr_sommerfeld_eigs(..., bc="decay")` the viscous one for the same layer.
2. Raw spectra (`filter=False`) are shown so that the real-axis line is visible.""")
D("D22", ref="11.84")
code(r"""
mR = ch11.rayleigh_eigs(0.4, np.tanh, d2tanh, bc="decay", N=100, return_vectors=True)   # the growing mode of the tanh layer, with its shape phi(y)
idR = ch11.rayleigh_identity_check(0.4, mR["c"][0], mR["phi"][:, 0], mR["y"], np.tanh, d2tanh, grid=mR["grid"])   # both parts of D22's identity
print("c =", np.round(mR["c"][0], 4), "| residuals of the identity of step 5: real part", f"{idR['res_real']:.1e}", "imaginary part", f"{idR['res_imag']:.1e}")
yP = np.linspace(-1, 1, 201)                                        # plane Poiseuille flow: U'' = -2 everywhere, no sign change
cP = ch11.rayleigh_eigs(1.0, lambda y: 1 - y**2, lambda y: -2 + 0*y, N=100, unstable_only=True)   # converged unstable inviscid modes
print("plane Poiseuille, inviscid: unstable eigenvalues found:", len(cP))
""", explain=r"""
1. `return_vectors=True` returns the mode shape; `ch11.rayleigh_identity_check` evaluates the integral identity of D22 step 5
   (real and imaginary parts) with Clenshaw–Curtis weights: both residuals vanish to round-off for the computed growing mode.
2. For plane Poiseuille flow, whose curvature never changes sign, the inviscid solver finds no unstable mode — as the theorem says.""")
note("N93 [B]", r"""
**Fjørtoft's theorem.** The real part of Rayleigh's identity is
$\int\frac{U-c_r}{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=-\int\big(\lvert\phi'\rvert^2+k^2\lvert\phi\rvert^2\big)dy<0$
(11.85); because the integral of D22's step 8 is zero for a growing mode, any constant times it can be added:
$(c_r-U_I)\int\frac1{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy=0$ (11.86); the sum gives
$\int\frac{U-U_I}{\lvert U-c\rvert^2}\frac{d^2U}{dy^2}\lvert\phi\rvert^2dy<0$, so $(U-U_I)U''<0$ somewhere — the vorticity magnitude
has a **maximum** at the inflection point. ⚠️ The book writes $U_1$ for $U$ at the inflection point; we write $U_I$ ($U_1$ was the
upper stream in C02).""")
nb.worked_example("three profiles by hand", r"""
1. **Poiseuille** $U=1-y^2$: $U''=-2$ everywhere — no sign change, no inflection ⇒ inviscidly stable (Rayleigh).
2. **tanh** $U=\tanh y$: $U''=-2\tanh y\,\mathrm{sech}^2y$ — negative for $y>0$, positive for $y<0$ ⇒ inflection at $y=0$,
   $U_I=0$; $(U-U_I)U''=-2\tanh^2y\,\mathrm{sech}^2y\le0$ ⇒ Fjørtoft satisfied: instability possible (and real: $kc_i=0.19$).
3. **$\sin y$ on $\lvert y\rvert\le b$**: $U''=-\sin y$ changes sign at 0, and $(U-0)U''=-\sin^2y<0$: both criteria pass — yet for
   $2b<\pi$ the flow is stable (N95). The criteria are necessary, not sufficient.""")
code(r"""
print("panel | profile | inflection point (Rayleigh) | vorticity maximum (Fjortoft) | y_I")
for v in ch11.fig_11_21_verdicts():                                 # our stand-ins for six classic profiles, with the two verdicts
    print(f"  ({v['panel']}) {v['label']:34s} {str(v['rayleigh']):6s} {str(v['fjortoft']):6s} {np.round(v['y_I'], 3)}")
yT = np.linspace(-3.0, 3.0, 601)                                    # sample points across the tanh layer
print(ch11.rayleigh_criterion(yT, np.tanh, d2tanh))                 # (11.84): is there an inflection point, and where?
print({n: v for n, v in ch11.fjortoft_criterion(yT, np.tanh, d2tanh).items()})   # (11.86): is (U - U_I) U'' < 0 somewhere?
tm = ch11.tanh_max_growth()                                         # the fastest-growing inviscid wave of the tanh layer
bm = BENCH["tanh_shear_layer_inviscid"]                             # published (Michalke 1964)
print(f"fastest tanh wave: k = {tm['k']:.4f}, k c_i = {tm['kci']:.4f}, c = {tm['c'].imag:.4f} i   (published {bm['k_max']}, {bm['kci_max']})")
""", explain=r"""
1. `ch11.fig_11_21_verdicts()` lists six profiles — (a) Blasius, (b) plane Couette, (c) plane Poiseuille, (d) a profile whose
   vorticity is largest at the walls, (e) $\sin y$, (f) a shear layer between walls — with the two verdicts: (a)–(c) have no
   inflection point; (d) has one but fails Fjørtoft; (e) and (f) pass both.
2. `ch11.rayleigh_criterion(y, U, Upp)` and `ch11.fjortoft_criterion(y, U, Upp)` apply the two necessary conditions to one profile.
3. `ch11.tanh_max_growth()` maximises $kc_i$ over $k$ for the tanh layer: $k=0.445$, $kc_i=0.190$, matching the published value.""")
scratch(r"""
yS = np.linspace(-3.0, 3.0, 600)                                    # an even number of points, so that y = 0 is not a sample
s = np.sign(d2tanh(yS))                                             # the sign of U'' at each point
idx = np.nonzero(np.diff(s))[0]                                     # where the sign changes between neighbours (Ch. 7 P180)
y_I = 0.5*(yS[idx] + yS[idx + 1])                                   # the inflection point: midway between the two samples
assert np.allclose(y_I, ST.inflection_points(yS, Upp=d2tanh(yS)), atol=yS[1] - yS[0])   # same point as the library
fj = np.min((np.tanh(yS) - np.tanh(y_I[0]))*d2tanh(yS))             # Fjortoft's product (U - U_I) U'': is it negative somewhere?
print(f"inflection at y = {y_I[0]:.3f}; min of (U - U_I) U'' = {fj:.3f} < 0 -> Fjortoft satisfied")
""", r"""
1. `np.sign` and `np.diff` find where $U''$ changes sign — the inflection point.
2. Fjørtoft's product $(U-U_I)U''$ is negative: the vorticity magnitude peaks at the inflection point.""")
fig(r"""
fig, axs = plt.subplots(2, 3, figsize=(8.6, 4.8))
six_profiles(axs)                                                  # drawing helper: the six profiles with inflection points and verdicts
fig.supxlabel("velocity $U$ (each panel: its own scale; the thin vertical line is $U=0$)", fontsize=9)   # one label for all panels
fig.supylabel("cross-stream coordinate $y$", fontsize=9)
plt.show()
""",
    see="Six velocity profiles ($U$ horizontal, $y$ vertical): (a) Blasius boundary layer, (b) plane Couette, (c) plane Poiseuille, "
        "(d) a profile that is steepest at the walls, (e) a sine profile, (f) a shear layer between walls. Rose dots mark "
        "inflection points; each title gives the two verdicts.",
    read="Only (e) and (f) pass both tests: they have an inflection point *and* the vorticity magnitude is largest there. (d) has "
         "an inflection point, but its vorticity is largest at the walls, not inside — Fjørtoft rules it out.",
    change="…an adverse pressure gradient bends a boundary-layer profile toward the shape of (f) (Ch. 9: the curvature at the "
           "wall changes sign): it becomes inviscidly unstable (R24, C13).")
note("N94 [B]", r"""
**Which flows are suspect?** Jets, wakes, mixing layers and adverse-pressure boundary layers have an inflection point with a
vorticity maximum: potentially unstable at any Re. Couette, Poiseuille and favourable or zero-gradient boundary layers have none:
inviscidly stable — any instability they have must come from viscosity (C13).""")
note("N95 [B]", r"""
**Not sufficient.** $U=\sin y$ between walls at $y=\pm b$ has an inflection point that satisfies both criteria, yet it is stable
for $2b<\pi$ (Tollmien).""")
code(r"""
for b in ((1.4, 2.0) if FAST else (1.4, 1.7, 2.0, 3.0)):            # half-widths: 2b < pi for 1.4, 2b > pi for the others
    print(f"b = {b}: 2b {'<' if 2*b < np.pi else '>'} pi, largest growth rate k c_i = {ch11.sin_profile_max_growth(b):.4f}")
""", explain=r"""
`ch11.sin_profile_max_growth(b)` scans $k$ and returns the largest inviscid growth rate: zero for $b=1.4$, positive beyond
$b=\pi/2\approx1.571$, and growing with $b$.""")
note("N96 [B]", r"""
**Critical layers.** A neutral mode has a real $c$, and Howard's result $U_{\min}<c_r<U_{\max}$ (11.71) puts it inside the velocity
range, so $U(y_c)=c$ at some $y_c$: the *critical layer*, a singular point of Rayleigh's equation (P269) where $\phi$ may be
discontinuous. The full Orr–Sommerfeld equation has no such singularity; a thin viscous layer forms there instead, thinning as Re
grows.""")
code(r"""
print("critical layer of U = tanh y for c = 0.3: y_c =", np.round(ch11.critical_layer(np.linspace(-2, 2, 401), np.tanh, 0.3), 4),
      "| artanh(0.3) =", round(np.arctanh(0.3), 4))
""")
note("N97 [B]", r"""
**Kelvin's cat's eye.** Seen by an observer moving at $c$, the basic flow is $U-c$, with stream function $\int(U-c)\,dy$; adding
the neutral wave, $\hat\psi=\int(U-c)\,dy+A\phi(y)\exp\{\mathrm ikx\}$ (11.87). Near $y_c$ a Taylor expansion ($U(y_c)=c$ kills the
linear term) gives $\hat\psi\cong\frac{(y-y_c)^2}2\big[\frac{dU}{dy}\big]_{y=y_c}+A\phi(y_c)\cos(kx)$: closed "eyes" of half-height
$2\sqrt{A\phi_c/U'_c}$ (our formula) around the critical layer.""")
fig(r"""
Xc, Yc = np.meshgrid(np.linspace(0, 4*np.pi, 241), np.linspace(-1.2, 1.2, 121))   # two wavelengths, around the critical layer
psi_c = ch11.cats_eye_streamfunction(Xc, Yc, y_c=0.0, A=0.1, phi_c=1.0, k=1.0, Uy_c=1.0)   # (11.87) expanded near y_c
half = ch11.cats_eye_width(0.1, 1.0, 1.0)                           # half-height of the eye
fig, ax = plt.subplots(figsize=(7.4, 2.9))
ax.contour(Xc, Yc, psi_c, levels=np.linspace(-0.08, 0.7, 14), colors=C_BASE, linewidths=0.7)
ax.contour(Xc, Yc, psi_c, levels=[0.1], colors=C_NEUT, linewidths=2.2)      # the separatrix: the streamline through the saddle points
ax.plot([np.pi, np.pi], [-half, half], color=C_GROW, lw=2)
ax.text(np.pi + 0.15, 0.3, f"±{half:.3f}", color=C_GROW, fontsize=8)
ax.set(xlabel="$kx$", ylabel="$y-y_c$")
ax.set_title("Kelvin's cat's eyes around a critical layer", fontsize=10)
plt.show()
print("half-height of the eye:", round(half, 4), "= 2 sqrt(A phi_c / U'_c) with A = 0.1")
""",
    see="Streamlines seen by an observer riding on the wave: wavy lines above and below, and between them a chain of closed loops — "
        "the eyes — bounded by the bold purple streamline. The rose bar marks the height of one eye.",
    read="Fluid inside an eye circulates with the wave; fluid outside passes by, to the right above and to the left below. This is "
         "the billow of C02 in its linear infancy.",
    change="…we double the wave amplitude A: the eyes grow taller by √2.",
    explain=r"""
1. `ch11.cats_eye_streamfunction(x, y, y_c, A, phi_c, k, Uy_c)` evaluates the expanded stream function near the critical layer;
   `ch11.cats_eye_width` its half-height $2\sqrt{A\phi_c/U'_c}=0.632$.
2. The level $\hat\psi=A\phi_c$ is the separatrix (it passes through the saddle points between the eyes).""")
note("N126 [B]", r"""
**A shear layer you can solve by hand (Exercise 11.11).** A piecewise-linear layer (speed $U_1$ above, $U_3$ below, linear in
between, thickness $h$) gives $c_o^2=\big(\frac{U_1-U_3}{2kh}\big)^2\{(kh-1)^2-e^{-2kh}\}$ for the wave speed $c_o$ relative to the
mean: unstable when the brace is negative, neutral at $kh=1.2785$ (our root of $(kh-1)^2=e^{-2kh}$).""")
code(r"""
kh_n = ch11.piecewise_neutral_kh()                                  # root of (kh - 1)^2 = exp(-2 kh)
c_pw = ch11.piecewise_shear_layer_c(0.797)                          # wave speed at kh = 0.797 for a unit velocity difference
print(f"neutral at kh = {kh_n:.4f} (reference {BENCH['piecewise_shear_layer_neutral_kh']['kh']}); at kh = 0.797: c = {c_pw:.4f}, growth k c_i h/dU = {0.797*c_pw.imag:.4f}")
""")
note("N127 [B]", r"""
**Symmetric and antisymmetric profiles (Exercise 11.12).** Written for $\hat v=-\mathrm ik\phi$, Rayleigh's equation is
$(U-c)\big(\frac{d^2\hat v}{dy^2}-k^2\hat v\big)-\frac{d^2U}{dy^2}\hat v=0$ (11.95). If $U$ is odd (a shear layer) and $c$ is an
eigenvalue so is $-c^*$ (waves come in mirror pairs); if $U$ is even (a jet) the modes split into *sinuous* ($\hat v$ even: the
jet wiggles) and *varicose* ($\hat v$ odd: it pulses).""")
code(r"""
bk = ch11.inviscid_profile("bickley")                               # the plane jet U = sech^2 y (R23 in the next section)
c_sin = ch11.rayleigh_eigs(1.0, bk["U"], bk["Upp"], bc="decay", parity="even", N=100)[0]   # sinuous mode at k = 1
c_var = ch11.rayleigh_eigs(0.5, bk["U"], bk["Upp"], bc="decay", parity="odd", N=100)[0]    # varicose mode at k = 0.5
print("jet: sinuous c =", np.round(c_sin, 4), "(k = 1) | varicose c =", np.round(c_var, 4), "(k = 0.5): both grow, the sinuous one faster")
""")
explainer("inviscid_shear_criteria", "Which profiles can be unstable without viscosity?", r"""
Switching profiles moves the inflection marker, flips the verdict badges and moves the eigenvalues on the $c$-plane — always
inside the semicircle — and clicking a neutral mode opens the cat's eye in the frame moving with $c$: each theorem is checked
against an example instead of being read.""", [
    "Preset 'wall-vorticity-max': Rayleigh ✓ but Fjørtoft ✗ — find where (U − U_I)U″ is positive.",
    "Preset 'sin y' and drag b across π/2: the eigenvalue leaves the real axis only when 2b > π.",
    "Open the cat's eye and double A: the eye grows by √2.",
])
whatif(r"""
…the profile has no inflection point at all, like plane Poiseuille flow? Rayleigh says no inviscid instability — but C11's
spectrum already showed a growing mode at $\mathrm{Re}=10^4$. Viscosity itself must be doing it: C13.""")
# =====================================================================================================================
# A.10 §11.10 Results for Parallel and Nearly Parallel Viscous Flows — R23 R24, C13 · C14
# =====================================================================================================================
nb.section("11.10", "Results for Parallel and Nearly Parallel Viscous Flows", intro=r"""
**What is this section about?** With the Orr–Sommerfeld solver of C11 we can trace, for any parallel flow, the curve in the
(Re, k) plane where waves neither grow nor decay. Inflectional flows (jets, mixing layers) go unstable at tiny Re; plane
Poiseuille flow, which Rayleigh calls stable, goes unstable at Re = 5772 because of viscosity; Couette and pipe flow never do,
linearly. C14 explains how viscosity can destabilise, through the energy budget of the disturbance.""")
nb.recap("R23", "The Bickley jet", r"""
The plane laminar jet of Ch. 9, $u=u_0\,\mathrm{sech}^2\big(\frac\eta{\sqrt6}\big)$ (9.71), i.e. $U=U_0\,\mathrm{sech}^2(y/L)$
(`JET.free_jet_profile`). **New here:** it is inflectional and very unstable: its sinuous mode goes unstable at
$\mathrm{Re}=U_0L/\nu\approx4$ at $kL\approx0.2$ (Tatsumi & Kakutani 1958, approximate; our value 4.017 at 0.173); its inviscid neutral modes are at $k=2$
(sinuous) and $k=1$ (varicose) with $c=2/3$ (Drazin & Reid 1981).""", where="Ch. 9, free jets")
nb.recap("R24", "Boundary layers with a pressure gradient", r"""
Ch. 9's rule: a favourable pressure gradient keeps the profile free of inflection points,
$\big(\frac{\partial^2u}{\partial y^2}\big)_{wall}<0$ for $dp/dx<0$ (9.51); an adverse one creates one,
$\big(\frac{\partial^2u}{\partial y^2}\big)_{wall}>0$ for $dp/dx>0$ (9.52). **New here:** the neutral curves of Falkner–Skan
profiles — favourable and zero gradient close into a loop as $\mathrm{Re}\to\infty$ (inviscidly stable), adverse ones keep a flat
upper branch (inflectional, inviscidly unstable) and go unstable at lower Re (`ch11.falkner_skan_neutral_curve`).""",
         where="Ch. 9 §9.7")

core("C13", r"Plane Poiseuille flow: viscously unstable above $\mathrm{Re}_c=5772.22$ at $k_c=1.02056$ (Orszag 1971), and the family of Table 11.1",
     "How can a flow that Rayleigh's theorem calls stable become unstable when viscosity is added?")
problem(r"""
Water flowing between two plates has the parabolic profile of Ch. 8 — no inflection point, so by C12 no inviscid instability. Yet
the Orr–Sommerfeld spectrum of C11 had a growing wave at $\mathrm{Re}=10^4$. Where exactly does it start, which wavelength comes
first, and how do other flows of the same family compare?""")
idea(r"""
 k ▲      ____                    neutral curve c_i(k, Re) = 0: a 'thumb'
   │   .-'    '----.___           inside: Tollmien–Schlichting (TS) waves grow
   │  (   unstable      '-----    outside: everything decays
   │   '-.___ ◆_______________    ◆ = the critical point (Re_c, k_c): first instability as Re rises
   └──────────────────────────► log Re
 inviscid profiles (tanh, jets): the thumb reaches Re → 0;  Poiseuille: starts at 5772;  Couette, pipe: no thumb at all
""")
remind("C13")
note("N98 [C]", r"""
**Viscosity cuts both ways** (named): at very large Re it damps the shortest waves; near the walls it shifts the phase of the
motion so that the wave can extract energy (C14). The asymptotic theory (Heisenberg, Lin, Shen, Yih — wall and critical layers as
singular perturbations) is notoriously hard; we compute.""")
nb.worked_example("what the critical Reynolds number means for water", r"""
Channel half-width $L=1$ cm, water $\nu=10^{-6}$ m²/s.

1. $\mathrm{Re}=U_0L/\nu$ with $U_0$ the centreline speed: $\mathrm{Re}_c=5772$ ⇒ $U_0=5772\times10^{-6}/0.01=0.577$ m/s.
2. The mean speed of a parabola is $\tfrac23U_0=0.385$ m/s; a Reynolds number built on the mean speed and the full width $2L$
   would read $\tfrac23\times2\times5772=7696$ — always say which length and speed (Ch. 9, P200).
3. At $\mathrm{Re}=10^4$, $k=1$ (C11): $kc_i=0.00374$ in units $U_0/L$; with $L=1$ cm and $U_0=1$ m/s that is 0.374 s⁻¹, an
   e-folding time of 2.7 s, for a wave $2\pi L=6.3$ cm long travelling at $c_rU_0=0.24$ m/s.""")
code(r"""
crit_P = ch11.poiseuille_critical()                                 # the critical point (stored result of the two nested searches of P263)
ref_P = BENCH["plane_poiseuille_critical"]                          # published (Orszag 1971)
print(f"Re_c = {crit_P['Re_c']:.2f}, k_c = {crit_P['k_c']:.5f}, c_r = {crit_P['c_r']:.5f}   (published {ref_P['Re_c']}, {ref_P['k_c']}, {ref_P['c_r']})")
nc = ch11.poiseuille_neutral_curve()                                # the neutral curve: lower and upper unstable wavenumber for each Re
lo = np.interp(np.log(1e4), np.log(nc["Re"]), nc["k_lower"])        # read the curve at Re = 10^4 (interpolating in log Re)
hi = np.interp(np.log(1e4), np.log(nc["Re"]), nc["k_upper"])
print(f"at Re = 10^4 the wavenumbers from k = {lo:.2f} to {hi:.2f} grow")
print("plane Couette flow, leading c_i at k = 1: Re = 10^3:", round(ch11.couette_max_growth(1.0, 1e3), 4), "| Re = 10^4:", round(ch11.couette_max_growth(1.0, 1e4), 4))
""", explain=r"""
1. `ch11.poiseuille_critical()` returns the critical point — P263's two nested searches (zero of $c_i$ in Re, minimum over $k$)
   on the Orr–Sommerfeld eigenvalue; it is a stored result (Ch. 10 P252: expensive runs are cached), and the next cell recomputes
   it independently. It matches the published benchmark to all its digits.
2. `ch11.poiseuille_neutral_curve()` gives, for each Re, the lower and upper edge of the unstable band — the "thumb".
3. `ch11.couette_max_growth(k, Re)`: for plane Couette flow the leading $c_i$ is negative — it decays.""")
scratch(r"""
Re_mine = brentq(lambda Re: np.imag(ch11.orr_sommerfeld_eigs(1.02056, Re, U, Upp, N=80)[0]), 5000.0, 7000.0, xtol=1e-3)   # where c_i = 0 at k = k_c
assert np.isclose(Re_mine, crit_P["Re_c"], rtol=1e-4)               # one root search reproduces the stored critical point
print(f"c_i = 0 at k = 1.02056 for Re = {Re_mine:.2f}")
""", r"""
The stored critical point is not a black box: one `brentq` on the sign of $c_i(\mathrm{Re})$ at $k=k_c$, with a fresh
Orr–Sommerfeld solve at every trial Re, lands on 5772.22.""")
fig(r"""
def read_curve(name):                                              # one of our small public tables of neutral curves in reference/ch11
    return np.genfromtxt(ROOT / "reference" / "ch11" / f"os_neutral_{name}.csv", delimiter=",", skip_header=2, names=True)

bj, bl = read_curve("bickley"), read_curve("blasius")              # Bickley jet and Blasius boundary layer (Orr-Sommerfeld, ours)
tl = ch11.tanh_shear_layer_neutral_curve()                         # tanh shear layer: the upper neutral wavenumber for each Re
fig, (a, b) = plt.subplots(1, 2, figsize=(9.8, 3.8))
Re_th = np.concatenate([[crit_P["Re_c"]], nc["Re"]])                  # the thumb, closed at its tip by the critical point
k_lo, k_hi = np.concatenate([[crit_P["k_c"]], nc["k_lower"]]), np.concatenate([[crit_P["k_c"]], nc["k_upper"]])
a.fill_between(Re_th, k_lo, k_hi, color=C_GROW, alpha=0.15)
a.semilogx(Re_th, k_lo, color=C_NEUT)
a.semilogx(Re_th, k_hi, color=C_NEUT)
a.plot(crit_P["Re_c"], crit_P["k_c"], "D", color=COLORS["ink"])
a.plot([1e4, 1e4], [lo, hi], color=C_GROW, lw=3)
a.text(6200, 1.19, "critical point\n5772.22, $k$ = 1.02056", fontsize=8)
a.text(4e4, 0.8, "TS waves grow", color=C_GROW, fontsize=9)
a.set(xlabel="Re = $U_0L/\\nu$ [–]", ylabel="wavenumber $k$ [$1/L$]", xlim=(3e3, 1e6), ylim=(0.2, 1.4))
a.set_title("(a) plane Poiseuille: a thumb", fontsize=10)
b.fill_between(tl["Re"], 0*tl["Re"], tl["k_upper"], color=C_SHEAR, alpha=0.12)
b.semilogx(tl["Re"], tl["k_upper"], color=C_SHEAR, label="tanh shear layer")
b.semilogx(bj["Re"], bj["k_upper"], color=C_BUOY, label="Bickley jet (sinuous)")
b.semilogx(bj["Re"], bj["k_lower"], color=C_BUOY)                   # its lower edge where the table resolves it (gaps in the line: edge below k = 0.02)
b.semilogx(bl["Re"], bl["k_upper"], color=C_NEUT, label="Blasius (in $\\delta^*$ units)")
b.semilogx(bl["Re"], bl["k_lower"], color=C_NEUT)
b.semilogx(Re_th, k_hi, ":", color=C_BASE, label="plane Poiseuille")
b.semilogx(Re_th, k_lo, ":", color=C_BASE)
b.set(xlabel="Re [–] (each flow's own $L$ and $U_0$)", ylabel="wavenumber $k$ [$1/L$]", xlim=(1, 1e6), ylim=(0, 2.1))
b.set_title("(b) the family on one axis", fontsize=10)
b.legend(fontsize=7, loc="upper right")
savefig(fig, "ch11", "c13_neutral_curves")
plt.show()
print(f"first Re in each table: Bickley {bj['Re'][0]:.1f}, Blasius {bl['Re'][0]:.0f}, tanh {tl['Re'][0]:.0f} (unstable at every Re computed)")
bcj = ch11.bickley_critical()                                      # the jet's critical point (stored)
print(f"Bickley jet: critical point Re_c = {bcj['Re_c']:.3f} at k = {bcj['k_c']:.3f} (its bands are drawn in detail in the next figure)")
""",
    see="(a) The neutral curve of plane Poiseuille flow on a logarithmic Re axis: a thumb pointing left, tip (◆) at 5772.22; the "
        "rose bar is the unstable band at Re = 10⁴. (b) The same curve (grey dotted, far right) next to the neutral curves of the "
        "tanh shear layer (orange), the Bickley jet (blue) and the Blasius boundary layer (purple). The jet has, besides its long upper "
        "curve, two short blue pieces near the bottom: its lower edge, drawn only where the table resolves it (elsewhere it lies "
        "below the smallest wavenumber computed, k = 0.02). The next figure shows the jet alone, with its bands shaded.",
    read="Inside a curve waves grow. The inflectional flows — shear layer and jet — are unstable down to Re of order 1 and their "
         "curves stay open to the right: inviscid instability, which viscosity merely trims. The flows without an inflection "
         "point — Blasius and Poiseuille — have thumbs that start at hundreds or thousands and narrow as Re grows: "
         "instability only through viscosity.",
    change="…plane Couette flow: no curve at all (N100).",
    explain=r"""
1. `read_curve` reads one of our small public tables (`reference/ch11/os_neutral_*.csv`, computed by
   `scripts/ch11_neutral_curves.py`) with `np.genfromtxt`; the first two lines of each file are comments.
2. `ch11.tanh_shear_layer_neutral_curve()` returns the stored upper neutral wavenumber of the tanh layer (it has no lower branch:
   every longer wave grows).
3. Each flow uses its own length and velocity scale (Table 11.1 below), so only the orders of magnitude should be compared.""")
fig(r"""
bjl = read_curve("bickley_longwave")                               # the long-wave band of the jet (second public table, same Re values)
K0 = 0.02                                                          # the smallest wavenumber of the computation behind the tables [1/L]
Rej = bj["Re"]                                                     # Reynolds numbers of the table
unres = np.isnan(bj["k_lower"])                                    # NaN = the band reaches below K0: UNSTABLE there, lower edge not resolved
has_long = ~np.isnan(bjl["k_long_upper"])                          # where the upper edge of the long-wave band is resolved
gap = ~unres & (Rej > Rej[unres][0])                               # Re values with a stable gap under the main band (after the band has split)
fig, ax = plt.subplots(figsize=(7.6, 4.0))
i_first = int(np.argmax(gap))                                      # first tabulated Re at which the gap exists
Re_f = np.insert(Rej, i_first, 0.999*Rej[i_first])                 # one extra point just left of it, so that the shading does not
lo_f = np.insert(np.where(unres, K0, bj["k_lower"]), i_first, K0)  # ... draw a slanted (invented) edge between two table rows:
up_f = np.insert(bj["k_upper"], i_first, bj["k_upper"][i_first])   # ... up to that row the band is shaded down to K0, as the table says
ax.fill_between(Re_f, lo_f, up_f, color=C_GROW, alpha=0.15, label="main unstable band")
ax.fill_between(Rej[gap], np.where(has_long, bjl["k_long_upper"], K0)[gap], bj["k_lower"][gap], color=C_DECAY, alpha=0.35, label="stable gap")
ax.fill_between(Rej[has_long], K0, bjl["k_long_upper"][has_long], color=C_SHEAR, alpha=0.30, label="long-wave unstable band")
ax.loglog(Rej, bj["k_upper"], color=C_BUOY)                        # upper neutral wavenumber of the main band
ax.loglog(Rej, bj["k_lower"], color=C_BUOY)                        # its lower edge where resolved (before the split, and as the top of the gap after it)
ax.loglog(Rej[has_long], bjl["k_long_upper"][has_long], color=C_SHEAR)   # upper edge of the long-wave band = bottom of the gap
ax.plot(Rej[unres], np.full(unres.sum(), K0), "|", color=C_GROW, ms=9, label="unstable down to $k$ = 0.02 (edge not resolved)")
ax.plot(bcj["Re_c"], bcj["k_c"], "D", color=COLORS["ink"], ms=5)   # the critical point
ax.text(bcj["Re_c"]*1.12, bcj["k_c"]*0.72, f"critical point\nRe = {bcj['Re_c']:.3f}, $k$ = {bcj['k_c']:.3f}", fontsize=8)
ax.axhline(K0, color=C_BASE, lw=0.8, ls=":")
ax.set(xlabel="Re = $U_0L/\\nu$ [–]", ylabel="wavenumber $k$ [$1/L$]", xlim=(3.5, 1000), ylim=(0.014, 2.6))
ax.set_title("The Bickley jet alone: two unstable bands and a gap", fontsize=10)
ax.legend(fontsize=7, loc="center right")
savefig(fig, "ch11", "c13_bickley_bands")
plt.show()
print(f"one band up to Re = {Rej[i_first - 1]:.2f}; at Re = {Rej[i_first]:.2f} the gap spans k = {bjl['k_long_upper'][i_first]:.4f} ... {bj['k_lower'][i_first]:.4f}")
print(f"upper neutral wavenumber at Re = {Rej[-1]:.0f}: {bj['k_upper'][-1]:.3f} (inviscid limit: 2)")
""",
    see="The jet's stability map on logarithmic axes. Rose: the main unstable band, bounded above by the blue curve that rises "
        "toward k = 2. Orange: a second, long-wave unstable band at the bottom. Teal: a stable gap between the two. Small rose "
        "ticks along the dotted line k = 0.02: Reynolds numbers at which the instability reaches down to the smallest "
        "wavenumber computed — the true lower edge is somewhere below and is not resolved. The diamond is the critical point.",
    read="Up to Re ≈ 17 there is one unstable band, and its lower edge drops out of the computed range almost at once: very "
         "long waves are unstable. Our reading of why (not a result of the book): a sinuous wave much longer than the jet is wide "
         "shifts the jet sideways almost as a whole, so it adds little strain, and the viscous damping of a wave, of order k²/Re, "
         "is small. A caution: for waves hundreds of jet widths long the parallel-flow model is itself questionable, because a "
         "real jet spreads over such a distance (note N109 makes the same point for the boundary layer). Then, between the two tabulated Reynolds "
         "numbers printed below, a stable gap opens *inside* the band and splits it in two: wavenumbers in the teal region decay, "
         "while longer and shorter waves both grow (from Re ≈ 19 the two bands are carried by two different modes). The gap "
         "slides to longer waves as Re rises and leaves the computed range near Re ≈ 80. The orange region ends near Re ≈ 27 only "
         "because its upper edge drops below k = 0.02 there: for Re ≳ 29.5 the long-wave band continues below the computed range — "
         "that is the edge of the computation, not the end of the band. Nowhere does a missing line mean "
         "'stable': only the teal region and the white region outside the blue curve are.",
    change="…we could compute to wavenumbers below 0.02 (it needs ever larger boxes, twelve wavelengths-over-2π wide): the rose "
           "ticks would be replaced by the true lower edge, and the gap could be followed beyond Re ≈ 80.",
    explain=r"""
1. Two public tables of ours: `os_neutral_bickley.csv` (the band that contains the fastest-growing wave: `k_lower`, `k_upper`)
   and `os_neutral_bickley_longwave.csv` (the long-wave band: `k_long_upper`). They were computed in a box that grows with the
   wavelength, max(40, 12/k) (`ch11.bickley_neutral_curve`), because a long wave decays slowly away from the jet.
2. `NaN` in `k_lower` means "the unstable band reaches below k = 0.02" — unstable, edge unknown; the cell fills such places
   down to the dotted line and marks them with ticks.
3. `fill_between` shades the three regions; the printed lines give where the gap first appears in the table and how close the
   upper branch has come to the inviscid neutral wavenumber 2 (R23).""")
nb.md(r"""
**Slider figure F7 — the one eigenvalue that crosses.** The Orr–Sommerfeld spectrum of plane Poiseuille flow at $k=1$ for a
range of Reynolds numbers.""")
nb.plotly(r"""
def f7(Re):                                                        # the spectrum at k = 1 for one Reynolds number
    c = ch11.poiseuille_spectrum(1.0, Re, N=100)                   # all finite eigenvalues, sorted by c_i
    c = c[c.imag > -1.0]                                           # the part of the plane we show
    return {"spectrum": (c.real, c.imag), "least stable mode (TS wave)": (np.array([c[0].real]), np.array([c[0].imag])),   # all modes; the leading one
            "c_i = 0": (np.array([0.0, 1.0]), np.array([0.0, 0.0]))}                                     # the real axis: above it, growth

figF7 = slider_figure(f7, "Re", np.round(np.geomspace(1e3, 1e5, 8 if FAST else 16)), xlabel="c_r [U0]", ylabel="c_i [U0]",   # slider: Re (log-spaced)
                      title="Plane Poiseuille, k = 1: the TS eigenvalue crosses c_i = 0 between Re = 5000 and 6000",       # the message
                      xrange=(0, 1.05), yrange=(-1.0, 0.1), modes={"spectrum": "markers", "least stable mode (TS wave)": "markers"})   # fixed window
figF7.show()
""", explain=r"""
1. `ch11.poiseuille_spectrum(k, Re)` returns the whole (unfiltered) spectrum; the first entry is the least stable mode.
2. Drag Re up from 1000: the Y sharpens and the marked mode rises through the axis near Re ≈ 5800, then (far beyond 10⁴) sinks
   back — the two sides of the thumb at $k=1$.""")
note("N99 [B]", r"""
**The two-stream shear layer** $U=U_0\tanh(y/L)$ (orange in panel (b)): unstable for $0<kL<k_u(\mathrm{Re})$, with $k_u\to1$ as
$\mathrm{Re}\to\infty$ and a critical Reynolds number of zero — viscosity is only a small correction to an inviscid instability.
The inviscid neutral mode at $kL=1$ has $c=0$ and $\phi=\mathrm{sech}(y/L)$.""")
code(r"""
ys = sp.symbols("y", real=True)                                     # cross-stream coordinate (in units of L)
Us, phis = sp.tanh(ys), sp.sech(ys)                                 # the tanh layer and the claimed neutral mode
print("Rayleigh residual for phi = sech y, k = 1, c = 0:", sp.simplify((Us - 0)*(sp.diff(phis, ys, 2) - phis) - sp.diff(Us, ys, 2)*phis))
for Re in (10.0, 100.0, 400.0):                                     # the stored upper neutral wavenumber at three Reynolds numbers
    print(f"  Re = {Re:5.0f}: unstable for k < {np.interp(np.log(Re), np.log(tl['Re']), tl['k_upper']):.3f}")
""", explain=r"""
1. sympy confirms that $\phi=\mathrm{sech}\,y$ with $c=0$, $k=1$ satisfies Rayleigh's equation exactly (residual 0).
2. The stored viscous neutral curve: the upper neutral wavenumber rises toward 1 as Re grows.""")
note("N100 [B]", r"""
**Plane Couette flow** is linearly stable at every Re (every Orr–Sommerfeld eigenvalue has $c_i<0$ on the grid below); in
experiments it still becomes turbulent at Reynolds numbers of a few hundred (based on half the wall-speed difference and half the
gap) through finite disturbances — named.""")
code(r"""
ks_c = (0.1, 0.5, 1.0, 2.0, 3.0)                                    # wavenumbers
Res_c = (1e2, 1e3, 1e4, 1e5, 1e6) if not FAST else (1e2, 1e4, 1e6)  # Reynolds numbers
ci = np.array([[ch11.couette_max_growth(k, Re) for Re in Res_c] for k in ks_c])   # leading c_i of plane Couette flow
print(np.round(ci, 4))
print("largest c_i on the grid:", round(ci.max(), 4), "-> every mode decays")
""")
note("N101 [C]", r"""
**Pipe flow** (named): linearly stable too; it becomes turbulent through finite disturbances whose threshold depends on how quiet
the inlet is; recent work describes the transition as a "chaotic saddle" (Eckhardt et al. 2007). Pointer: Ch. 12.""")
note("N102 [B]", r"""
**Table 11.1, recomputed.** Our own critical Reynolds numbers next to published benchmarks, each with its length scale. The
book's rounded column is not reproduced.""")
code(r"""
import pandas as pd                                                 # tables: pd.DataFrame(list_of_dicts) shows rows and columns
pd.set_option("display.width", 200)                                 # let the table use the full width
tab11 = pd.DataFrame(ch11.table_11_1())                             # one row per flow (stored critical points, sources cited)
print(tab11[["flow", "U", "length_scale", "Re_c_ours", "k_c_ours", "benchmark", "source"]].to_string(index=False))
""", explain=r"""
`ch11.table_11_1()` returns one dictionary per flow — jet, shear layer, Blasius, plane Poiseuille, pipe, plane Couette — with our
critical Reynolds number and wavenumber, the published benchmark and its source; `pd.DataFrame` turns the list into a table.
Infinite means "linearly stable at every Re".""")
note("N103 [C]", r"""
**A caveat on a zero critical Reynolds number** (named): the mixing layer's zero assumes a parallel flow; a spreading mixing layer
has a finite one (Bhattacharya et al. 2006).""")
whatif(r"""
…we want to know *why* viscosity can feed the wave? Follow the wave's kinetic energy: what produces it and what dissipates it —
C14.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C14", r"The disturbance kinetic-energy equation $\frac{d}{dt}\int\tfrac12u_i^2\,dV=-\int u_iu_j\frac{\partial U_i}{\partial x_j}dV-\Lambda$ (11.88)",
     "Where does a growing wave get its energy, and how can viscosity help it?")
problem(r"""
A growing wave gains kinetic energy; something must supply it. The only reservoir is the mean flow's shear. And something always
takes energy away: viscous dissipation. Bookkeeping the two tells us when a wave grows — and reveals the trick by which viscosity,
a pure damper, can make the supply larger than the loss.""")
idea(r"""
 d/dt ∫ ½(u² + v²) dV   =   −∫ uv dU/dy dV      −   Λ
      wave's energy          PRODUCTION            DISSIPATION ν∫(∂u_i/∂x_j)² ≥ 0
                             Reynolds stress −uv working on the shear U′
 u and v exactly out of phase (u = sin, v = cos): mean of uv = 0 → NO production (inviscid, no inflection)
 viscosity near the wall tilts the phase of v against u → mean uv < 0 where U′ > 0 → production > 0
""")
remind("C14")
P("P273", "the integral of an x-derivative of a periodic function over one period is zero", r"""
If $f(x)$ repeats every $\lambda$, then $\int_0^\lambda\frac{\partial f}{\partial x}dx=f(\lambda)-f(0)=0$ (the fundamental theorem
of calculus, Ch. 2 P84). Averaged over a whole number of wavelengths, every term of the form $\partial(\dots)/\partial x$
disappears — the two ends of the control volume cancel.""", code=r"""
import numpy as np                                    # arrays
x = np.linspace(0, 2*np.pi, 2001)                     # one full period
f = np.sin(3*x)**2 + np.cos(x)                        # a function that repeats every 2 pi
print(np.trapezoid(np.gradient(f, x), x))             # ~0: the integral of df/dx over a period (to discretisation error)
""")
note("N128 [B]", r"""
**The starting line (Exercise 11.13):** the disturbed Navier–Stokes equation in index form (Ch. 2: a repeated index is summed),
$\frac{\partial}{\partial t}(U_i+u_i)+(U_j+u_j)\frac{\partial}{\partial x_j}(U_i+u_i)=-\frac1\rho\frac{\partial}{\partial
x_i}(P+p)+\nu\frac{\partial^2}{\partial x_j\partial x_j}(U_i+u_i)$ (11.96), with $U_i$, $P$ the basic flow and $u_i$, $p$ the
disturbance.""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.0, 2.8))
cv_sketch(ax)                                                      # drawing only: the control volume between the walls
plt.show()
""",
    see="The control volume used in the derivation: a box between the two walls, a whole number of wavelengths long.",
    read="Nothing crosses the walls (no slip), and whatever leaves through one end enters through the other (periodicity, P273): "
         "only what happens *inside* the box can change the wave's energy.",
    change="…the box were not a whole number of wavelengths long: the two ends would no longer cancel and the budget would not close.")
D("D23", ref="11.88", check_src=r"""
en = ch11.energy_equation_sympy()                     # D23 with sympy on generic 2-D fields: multiply by u_i, sort into the three boxes
print("residual of the energy identity:", en["residual"])   # 0: steps 2-7 are exact
for t in en["terms"]:                                 # the engine's own names for the terms
    print("  -", t)
""")
note("N104 [B]", r"""
**The 2-D form and what it says:** $\frac d{dt}\int\tfrac12(u^2+v^2)\,dV=-\int uv\frac{\partial U}{\partial y}dV-\Lambda$. The
product $uv$ averaged over a wavelength (the *Reynolds shear stress* of Ch. 12, up to a factor $-\rho$) vanishes when $u$ and $v$
are 90° out of phase (Ch. 7 P178). For an inviscid, inflection-free profile the phases are such that the wave cannot take energy
from the shear; viscosity changes the phase relation near the walls and the critical layer so that $-\langle uv\rangle U'$
integrates to a positive number larger than the dissipation — the destabilising mechanism.""")
nb.worked_example("why the phase matters", r"""
Average over one period (Ch. 7 P178: the mean of a product of two sines is half the cosine of their phase difference).

1. $u=\sin t$, $v=\cos t$ (90° apart): mean of $uv=\tfrac12\sin2t$ is 0 — no production.
2. $u=\sin t$, $v=-\sin(t-\pi/6)$ (tilted by 30° from anti-phase): mean $uv=-\tfrac12\cos(\pi/6)=-0.433$; with $U'=+1$ the
   production $-\langle uv\rangle U'=+0.433>0$.
3. For the TS wave at $\mathrm{Re}=10^4$, $k=1$ the code below finds production/dissipation ≈ 1.62: production beats dissipation
   by 62 %, so $dE/dt=2kc_iE>0$.
4. At $\mathrm{Re}=5000$ the same mode has a ratio of about 0.77: it decays.""")
code(r"""
budgets = {}                                                        # the budget at two Reynolds numbers
for Re in (1e4, 5e3):
    m = ch11.os_mode(1.0, Re, U, Upp)                               # leading Orr-Sommerfeld mode at k = 1: c, nodes y, phi, u_hat, v_hat
    b = ch11.disturbance_energy_budget(1.0, m["c"], m["phi"], m["y"], -2*m["y"], Re)   # (11.88) averaged over a wavelength; U' = -2y
    budgets[Re] = (m, b)
    print(f"Re = {Re:.0f}: c = {m['c']:.5f}")
    print(f"   E = {b['E']:.5f}, dE/dt = 2 k c_i E = {b['dEdt']:+.6f}, production P = {b['production']:.6f}, dissipation = {b['dissipation']:.6f}")
    print(f"   P - dissipation - dE/dt = {b['residual']:.1e},  P/dissipation = {b['ratio']:.3f}")
""", explain=r"""
1. `ch11.os_mode(k, Re, U, Upp)` returns the leading mode: `c`, the nodes `y`, `phi`, and the velocity amplitudes
   $\hat u=\phi'$ (`u_hat`) and $\hat v=-\mathrm ik\phi$ (`v_hat`) — §11.8's convention.
2. `ch11.disturbance_energy_budget(k, c, phi, y, Up, Re)` evaluates the wavelength-averaged budget with Clenshaw–Curtis weights
   (P271): $E=\tfrac14\int(\lvert\hat u\rvert^2+\lvert\hat v\rvert^2)dy$, production $P=-\tfrac12\int\mathrm{Re}(\hat u\hat v^*)U'dy$,
   dissipation $\Lambda$, and $dE/dt=2kc_iE$.
3. The closing check — $dE/dt=P-\Lambda$ to about 10⁻¹² — is independent of the eigen-solve: the budget really balances.
4. $P/\Lambda=1.62$ at $\mathrm{Re}=10^4$ (grows), 0.77 at 5000 (decays).""")
scratch(r"""
m, b = budgets[1e4]                                                 # the growing mode
Cheb = np.polynomial.chebyshev.Chebyshev                            # numpy's Chebyshev series (re-interpolation of the nodal values)
yf = np.linspace(-1, 1, 4001)                                       # a fine uniform grid across the channel
deg = len(m["y"]) - 1                                               # the polynomial degree that passes through all nodes
uh = Cheb.fit(m["y"], m["u_hat"].real, deg)(yf) + 1j*Cheb.fit(m["y"], m["u_hat"].imag, deg)(yf)   # u-hat on the fine grid
vh = Cheb.fit(m["y"], m["v_hat"].real, deg)(yf) + 1j*Cheb.fit(m["y"], m["v_hat"].imag, deg)(yf)   # v-hat on the fine grid
P_mine = -0.5*np.trapezoid(np.real(uh*np.conj(vh))*(-2*yf), yf)    # P = -(1/2) int Re(u-hat conj(v-hat)) U' dy with U' = -2y
duh, dvh = np.gradient(uh, yf), np.gradient(vh, yf)                 # their y-derivatives by finite differences
L_mine = np.trapezoid(abs(uh)**2 + abs(vh)**2 + abs(duh)**2 + abs(dvh)**2, yf)/(2*1e4)   # dissipation with k = 1
assert np.isclose(P_mine, b["production"], rtol=1e-4) and np.isclose(L_mine, b["dissipation"], rtol=1e-3)
print(f"production by trapezoid: {P_mine:.6f} (library {b['production']:.6f}); dissipation: {L_mine:.6f} (library {b['dissipation']:.6f})")
""", r"""
1. `np.polynomial.chebyshev.Chebyshev.fit(y, values, deg)` builds the polynomial through the nodal values, so the mode can be
   evaluated on a fine uniform grid (real and imaginary parts separately).
2. Production and dissipation with the plain trapezoid rule (Ch. 9 P203) and `np.gradient` agree with the library's
   Clenshaw–Curtis values.""")
fig(r"""
fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.6))
m, b = budgets[1e4]
y_c = np.sqrt(1 - m["c"].real)                                     # critical layers: U(y) = 1 - y^2 = c_r
axs[0].plot(-b["uv"], b["y"], color=C_SHEAR, label=r"$-\langle uv\rangle$")
axs[0].plot(b["production_density"], b["y"], color=C_GROW, label=r"production $-\langle uv\rangle U'$")
for s in (-1, 1):
    axs[0].axhline(s*y_c, color=C_BASE, ls=":", lw=1)                # the two critical layers, where the flow moves at the wave speed
axs[0].set(xlabel="per unit height [–]", ylabel="$y$ [$L$]")
axs[0].set_title("(a) where the energy is produced", fontsize=10)
axs[0].legend(fontsize=7, loc="center right")
ph = np.degrees(np.angle(np.exp(1j*b["phase_uv"])))                # phase of v relative to u, wrapped to (-180, 180] degrees
ok = np.abs(b["y"][1:-1]) > 0.03                                   # leave out the centreline, where the amplitudes vanish and the phase is undefined
axs[1].plot(ph[1:-1][ok], b["y"][1:-1][ok], ".", ms=3, color=C_NEUT)   # the phase between v and u at each height
for x0 in (-90, 90):
    axs[1].axvline(x0, color=C_BASE, ls=":", lw=1)
axs[1].set(xlabel="phase of $v$ relative to $u$ [deg]", xlim=(-180, 180), xticks=[-180, -90, 0, 90, 180])
axs[1].set_title("(b) the phase is not exactly ±90°", fontsize=10)
for j, Re in enumerate((1e4, 5e3)):
    bb = budgets[Re][1]
    vals = [bb["production"], bb["dissipation"], bb["dEdt"]]
    axs[2].bar(np.arange(3) + 0.38*j - 0.19, vals, 0.36, color=[C_SHEAR, C_BUOY, C_GROW if bb["dEdt"] > 0 else C_DECAY], alpha=1.0 if j == 0 else 0.55)
axs[2].set_xticks(range(3), ["$P$", r"$\Lambda$", "$dE/dt$"])
axs[2].axhline(0, color=COLORS["ink"], lw=0.8)
axs[2].set_title("(c) budget: Re = 10$^4$ (solid), 5000 (pale)", fontsize=10)
savefig(fig, "ch11", "c14_budget")
plt.show()
near_wall = (np.abs(b["y"]) > 0.6) & (np.abs(b["y"]) < 0.999)        # the two bands between the walls and the middle of the channel
print(f"critical layers at y = +-{y_c:.3f}; largest departure of the phase from 90 degrees near the walls: {np.max(np.abs(np.abs(ph[near_wall]) - 90)):.0f} degrees")
""",
    see="(a) Across the channel: the Reynolds stress (orange) and the local production (rose) of the growing wave; dotted lines "
        "mark the critical layers, where the flow speed equals the wave speed. (b) The phase between $v$ and $u$ at each height. "
        "(c) Bars of production, dissipation and their difference for two Reynolds numbers.",
    read="Production is concentrated in two thin bands around the critical layers, near each wall. There the phase departs from "
         "±90° (by the number of degrees printed) — and where it does, the stress is nonzero. In the middle of the channel the "
         "phase is ±90° and nothing is produced. At Re = 10⁴ production exceeds "
         "dissipation and $dE/dt>0$ (rose bar); at 5000 the pale bars show production below dissipation and $dE/dt<0$ (teal).",
    change="…the phase were exactly ±90° everywhere (the inviscid, inflection-free case): the orange curve would vanish and the "
           "wave could only decay.",
    explain=r"""
1. The budget dictionary also holds profiles: `uv` ($\langle uv\rangle$ at each $y$), `production_density`, and `phase_uv` (the
   phase of $\hat v$ relative to $\hat u$ in radians).
2. `np.angle(np.exp(1j*phase))` wraps a phase into (−180°, 180°].""")
nb.md(r"""
**Animation A4 — the Tollmien–Schlichting wave.** The growing mode at $\mathrm{Re}=10^4$, $k=1$: perturbation stream function
(colour) travelling at $c_r\approx0.24$ under the parabolic flow, with the Reynolds-stress profile beside it.""")
nb.animation(r"""
nf = 14 if FAST else 18                                            # number of frames
m, b = budgets[1e4]                                                # the growing mode and its budget
xa = np.linspace(0, 4*np.pi, 90)                                   # two wavelengths in x [L]
period = 2*np.pi/(1.0*m["c"].real)                                 # time for the wave to travel one wavelength [L/U0]
tA4 = np.linspace(0, period, nf)                                   # one period
f0 = ch11.ts_wave_fields(xa, m["y"], 0.0, 1.0, m["c"], m["phi"], m["u_hat"], amp=0.05)   # real fields at t = 0
vmax = 1.15*np.max(np.abs(f0["psi"]))                              # fixed colour range
fig, (a, a2) = plt.subplots(1, 2, figsize=(8.6, 3.0), gridspec_kw={"width_ratios": [3, 1]}, sharey=True)
mesh = a.pcolormesh(xa, m["y"], f0["psi"], cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
a.set(xlabel="$x$ [$L$]", ylabel="$y$ [$L$]")
ttl = a.set_title("", fontsize=9)
a2.plot(-b["uv"], b["y"], color=C_SHEAR)
a2.set(xlabel=r"$-\langle uv\rangle$")
a2.set_title(f"P/Λ = {b['ratio']:.2f}", fontsize=9)

def update(i):                                                     # frame i: time tA4[i]
    f = ch11.ts_wave_fields(xa, m["y"], tA4[i], 1.0, m["c"], m["phi"], m["u_hat"], amp=0.05)   # the wave moved and grew a little
    mesh.set_array(f["psi"].ravel())
    ttl.set_text(f"TS wave, t = {tA4[i]:5.1f} L/U0: amplitude x {np.exp(m['c'].imag*tA4[i]):.3f}")
    return mesh, ttl

show_animation(animate(update, frames=nf, fig=fig, interval=130), player="video", dpi=60)
print(f"one period = {period:.1f} L/U0; growth over one period: factor {np.exp(1.0*m['c'].imag*period):.3f}")
""", explain=r"""
1. `ch11.ts_wave_fields(x, y, t, k, c, phi, phi_y, amp)` returns the real disturbance fields of the mode at time $t$: `psi`, `u`,
   `v` on the $(y,x)$ grid.
2. Over one period the pattern slides one wavelength downstream and its amplitude grows by about 10 %.""")
nb.figure_notes(
    see="Cells of the disturbance stream function sliding downstream; near each wall the cells are tilted — they lean against the "
        "shear. Beside them, the Reynolds-stress profile.",
    read="The tilt is the visible sign of a nonzero Reynolds stress: fluid moving away from the wall carries a deficit of "
         "$x$-velocity. Where the cells lean against the mean shear, the wave takes energy from the mean flow.",
    change="…$k=0.6$ (outside the thumb): the tilt weakens and the wave decays.")
note("N125 [C]", r"""
**With stratification (Exercise 11.10, named):** the budget gains the available potential energy,
$\tfrac12\frac d{dt}\int\big(u^2+w^2+\frac{g^2\rho^2}{\rho_0^2N^2}\big)dV=-\int uw\frac{\partial U}{\partial z}dV$ — the shear's
production now has to lift fluid too (C09's energy view; `ch11.stratified_energy_budget`). Pointer: available potential energy in
Ch. 13.""")
explainer("orr_sommerfeld_neutral_curve", "How can viscosity make a channel flow unstable?", r"""
Clicking a point (Re, k) on the neutral-curve map shows that Tollmien–Schlichting mode travelling over the profile, its
production/dissipation bars and the phase between $u$ and $v$; crossing the curve flips the balance. Switching flows recomputes
the Table 11.1 family row by row, and the oblique slider applies Squire's map live.""", [
    "Preset 'Orszag (Re = 10⁴, k = 1)': read c ≈ 0.2375 + 0.0037i and P/Λ ≈ 1.6.",
    "Drag Re down through 5772 at k = 1.02: the bars flip when the dot leaves the thumb.",
    "Set the oblique angle to 30°: the effective Re drops by cos 30° and the 3-D wave falls outside the curve (Squire).",
])
whatif(r"""
…the wave keeps growing? Its own Reynolds stress starts to change the mean flow it feeds on — the linear theory's end,
§11.12–11.13.""")

# =====================================================================================================================
# A.11 §11.11 — C13 (continued): N105–N109
# =====================================================================================================================
nb.section("11.11", "Experimental Verification of Boundary-Layer Instability", intro=r"""
**What is this section about?** The Blasius boundary layer on a flat plate (Ch. 9) is nearly parallel, so the Orr–Sommerfeld
machinery applies to it as an approximation. Its neutral curve was computed long before anyone saw a Tollmien–Schlichting wave; a
famous wind-tunnel experiment then confirmed it. (These notes continue `C13`.)""")
note("N105 [B]", r"""
**C13 (continued). Tollmien's stand-in for Blasius**: a piecewise profile in $\eta=y/\delta$ — linear near the wall, $a\eta$; a
parabola reaching 1 with zero slope at $\eta=1$, $1-b(1-\eta)^2$; then 1 — with zero curvature at the wall like Blasius. We fix
the breakpoint $\eta_1=0.2$ ourselves and get the coefficients from continuity of value and slope: $a=2/(1+\eta_1)$,
$b=1/(1-\eta_1^2)$.""")
slip(6, r"""the middle branch as $1-b[1-(y/\delta)^2]$ — the square inside the bracket; that branch jumps at the breakpoint""",
     r"""$1-b(1-y/\delta)^2$ — with the square outside, value and slope match the linear branch.""")
code(r"""
co = ch11.tollmien_coefficients()                                   # a and b from continuity at our breakpoint eta_1 = 0.2
eta = np.linspace(0, 1.2, 241)                                      # y/delta
eB = np.linspace(0, 1.2*4.91, 241)                                  # the same heights in Blasius' variable (delta = 99 % thickness = 4.91)
fB = BL.blasius_profile(eB)[1]                                      # Blasius u/U = f'(eta) (Ch. 9)
jump = ch11.tollmien_profile(np.array([0.2 - 1e-6, 0.2 + 1e-6]), printed=True)   # the printed middle branch on both sides of the breakpoint
print({n: round(v, 4) for n, v in co.items()}, "| jump of the printed form at the breakpoint:", round(abs(jump[1] - jump[0]), 3))
fig, ax = plt.subplots(figsize=(5.4, 3.2))
ax.plot(fB, eta, color=C_BASE, lw=3, alpha=0.5, label="Blasius")
ax.plot(ch11.tollmien_profile(eta), eta, color=C_NEUT, label="Tollmien's piecewise profile")
ax.plot(ch11.tollmien_profile(eta, printed=True), eta, ":", color=C_GROW, label="as printed (slip #6)")
ax.set(xlabel="$U/U_\\infty$", ylabel="$y/\\delta$", xlim=(-0.1, 1.1))
ax.legend(fontsize=8)
plt.show()
""", explain=r"""
1. `ch11.tollmien_coefficients()` gives $a=1.6667$, $b=1.0417$ for $\eta_1=0.2$; `ch11.tollmien_profile(eta)` the profile.
2. The purple curve follows Blasius (grey) closely; the dotted curve is the printed form, which jumps by about 0.33 at the
   breakpoint — it cannot be a velocity profile.""")
note("N106 [C]", r"""
**C13 (continued). Seen at last** (named): TS waves stayed undetected until Schubauer & Skramstad (1947) excited them with a
vibrating ribbon in a very quiet wind tunnel and measured them with hot wires; their growth and decay matched the computed neutral
curve.""")
note("N107 [B]", r"""
**C13 (continued). The Blasius neutral curve in frequency**: experiments fix the frequency $\omega$ of the ribbon, so the curve
is drawn as the dimensionless frequency $F=\omega\nu/U_\infty^2=kc_r/\mathrm{Re}_{\delta^*}$ against
$\mathrm{Re}_{\delta^*}=U_\infty\delta^*/\nu$ ($\delta^*$ the displacement thickness, Ch. 9); a wave of fixed $F$ travelling
downstream ($\mathrm{Re}_{\delta^*}$ rising like $\sqrt x$) enters the curve, grows, and leaves it.""")
fig(r"""
bF = ch11.blasius_neutral_curve(in_frequency=True)                  # stored Blasius neutral curve: F_lower, F_upper for each Re
fig, ax = plt.subplots(figsize=(6.2, 3.4))
ax.fill_between(bF["Re"], 1e6*bF["F_lower"], 1e6*bF["F_upper"], color=C_GROW, alpha=0.15)
ax.plot(bF["Re"], 1e6*bF["F_lower"], color=C_NEUT)
ax.plot(bF["Re"], 1e6*bF["F_upper"], color=C_NEUT)
ax.axhline(100, color=C_SHEAR, ls="--")
ax.text(3300, 108, "one ribbon frequency, $F$ = 100 × 10$^{-6}$", fontsize=8, color=C_SHEAR)
ax.set(xlabel="$\\mathrm{Re}_{\\delta^*}=U_\\infty\\delta^*/\\nu$ [–]", ylabel="$F=\\omega\\nu/U_\\infty^2$ [10$^{-6}$]", xlim=(400, 6000), ylim=(0, 260))
ax.set_title("Blasius boundary layer: where a fixed frequency grows", fontsize=10)
plt.show()
enter = np.interp(100, (1e6*bF["F_lower"])[::-1], bF["Re"][::-1])   # where the F = 100e-6 line meets the lower branch
leave = np.interp(100, (1e6*bF["F_upper"])[::-1], bF["Re"][::-1])   # ... and the upper branch
print(f"a wave with F = 100e-6 grows from Re = {enter:.0f} to Re = {leave:.0f}")
""",
    see="A tongue in the plane of Reynolds number (based on the displacement thickness) and dimensionless frequency; inside it "
        "(rose) waves grow. The dashed orange line is one fixed frequency.",
    read="Follow the dashed line from left to right: that is one ribbon frequency carried downstream. It enters the tongue, the "
         "wave grows, and further downstream it leaves the tongue and decays again (the two Reynolds numbers are printed).",
    change="…a higher frequency (say $F=250\\times10^{-6}$): the line misses the tongue — those waves never grow.")
note("N108 [B]", r"""
**C13 (continued). The Blasius critical point** with the parallel-flow assumption: published
$\mathrm{Re}_{\delta^*,c}=519.2$ at $k\delta^*=0.303$ (Thomas; via Gallagher, Griffiths & Stephen 2016; Jordinson 1970: 520). Our value, 519.06 (converged in the size of the computational box; a fixed box 20 displacement thicknesses high gave 519.08), is 0.03 % below it. Our
base flow uses $U=f'$ and $U''\propto-\tfrac12ff''$ from the Blasius equation (Ch. 9), not a numerical second derivative.""")
code(r"""
bc_ = ch11.blasius_critical()                                       # stored critical point of the Blasius layer (delta* units)
rb = BENCH["blasius_critical"]                                      # published
print(f"ours: Re_c = {bc_['Re_c']:.2f}, k_c = {bc_['k_c']:.5f}, c_r = {bc_['c_r']:.5f}, omega_c = {bc_['omega_c']:.5f}   (published {rb['Re_delta_star']}, {rb['alpha_delta_star']}, {rb['omega']})")
assert abs(bc_["Re_c"]/rb["Re_delta_star"] - 1) < 0.01 and abs(bc_["k_c"]/rb["alpha_delta_star"] - 1) < 0.01   # within 1 %
""")
note("N109 [C]", r"""
**C13 (continued). What the parallel model misses** (named): the boundary layer grows downstream (non-parallel corrections,
Nayfeh & Saric 1975); disturbances can grow for a while even where every eigenmode decays (transient growth from the non-normal
Orr–Sommerfeld operator, Reshotko 2001); suction thins the layer and delays transition. Pointer: Ch. 12.""")

# =====================================================================================================================
# A.12 §11.12, A.13 §11.13 — C14 (continued): N110, N111
# =====================================================================================================================
nb.section("11.12", "Comments on Nonlinear Effects", intro=r"""
**What is this section about?** Linear theory says a disturbance grows like $e^{\sigma t}$ forever; real disturbances stop
growing. The reason is that a finite wave changes the flow it lives on. (This note continues `C14`.)""")
note("N110 [B]", r"""
**C14 (continued). The wave changes its own supply.** The Reynolds stress $-\langle uv\rangle$ of C14 (and the heat flux
$\langle wT'\rangle$ in convection) is a *rectified* flux — its average is not zero: it transports momentum and heat across the
layer and so modifies $U(y)$ or $\bar T(z)$ — usually in the direction that weakens the instability, and the growth stops at a
finite amplitude (N52's equilibrium cells). Some flows go the other way: finite disturbances destabilise flows that are linearly
stable (pipe, N101). A rotating annulus heated at the rim shows the whole sequence — steady waves, then irregular ones — the
laboratory model of the atmosphere's baroclinic waves (pointer: Ch. 13 §13.17).""")
nb.section("11.13", "Transition", intro=r"""
**What is this section about?** How a laminar flow becomes turbulent: a sequence of instabilities, each feeding on the state the
previous one left behind. (This note continues `C14`.)""")
note("N111 [B]", r"""
**C14 (continued). A sequence of instabilities** (Landau's idea), drawn as two pipelines:

```
boundary layer / channel:  TS waves (2-D, C13) ──► saturated, then a secondary 3-D instability (Klebanoff's
                           peak–valley pattern) ──► Λ-shaped vortices ──► turbulent spots ──► turbulence
free shear layer:          KH roll-up (C02) ──► vortex pairing ──► 3-D 'braid' vortices ──► turbulence
                           (shear layers transition at Re of order 10; boundary layers need order 10³)
```

Pointer: Ch. 12 (turbulence) and the routes to chaos of C15.""")
# =====================================================================================================================
# A.14 §11.14 Deterministic Chaos — C15, S02, summary
# =====================================================================================================================
nb.section("11.14", "Deterministic Chaos", intro=r"""
**What is this section about?** Far above onset the steady cells and waves of this chapter give way to irregular motion.
Remarkably, irregularity does not need many degrees of freedom: three ordinary differential equations from a three-mode model of
convection (Lorenz 1963) never settle down and amplify any tiny difference in the start. This is why weather has a predictability
limit of about two weeks, and why forecasters run ensembles.""")
core("C15", r"The Lorenz system $\dot X=\Pr(Y-X)$, $\dot Y=-XZ+rX-Y$, $\dot Z=XY-bZ$ (11.91), $r=\mathrm{Ra}/\mathrm{Ra}_{cr}$, $b=4\pi^2/(\pi^2+k^2)$",
     "How can three exact equations be unpredictable?")
problem(r"""
In 1961 Edward Lorenz restarted a small weather model from numbers he had typed with three decimals instead of six. Within a
couple of simulated months the two runs had nothing in common. He then boiled the model down to convection in a layer — the Bénard
problem of C05 — keeping only three modes. The three equations still behaved the same way: completely determined, and yet
unpredictable beyond a horizon.""")
idea(r"""
free–free convection (C05)  ──keep 1 roll + 2 temperature modes──►  3 ODEs for X (flow), Y, Z (temperature)
   r < 1            : conduction (origin attracts)
   1 < r < r_H      : steady rolls C± (two senses of rotation)
   r > r_H ≈ 24.74  : every fixed point unstable ⇒ the state wanders forever between the two rolls: a STRANGE ATTRACTOR
two starts 10⁻⁸ apart: separation ∝ e^{λt} (λ ≈ 0.9) ⇒ after ~20–30 time units they are unrelated
""")
remind("C15")
note("N112 [B]", r"""
**Deterministic chaos** = extreme sensitivity to the initial state + aperiodic motion + a broadband spectrum (not a few sharp
frequencies). In a *linear* dissipative system steady forcing gives a steady response; with nonlinearity it can give a periodic or
an aperiodic one (the heated rotating annulus; the cylinder wake that starts to oscillate near $\mathrm{Re}\approx40$, Ch. 9).""")
fig(r"""
sol28 = ch11.lorenz_integrate((1.0, 1.0, 1.0), 40.0)               # Lorenz's own parameters: Pr = 10, b = 8/3, r = 28 (X, Y, Z on 4001 times)
sol10 = ch11.lorenz_integrate((1.0, 1.0, 1.0), 40.0, r=10.0)       # a gentler heating: r = 10
fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, 3.4))
a.plot(sol28["t"], sol28["X"], color=C_GROW, lw=0.9, label="$r$ = 28")
a.plot(sol10["t"], sol10["X"], color=C_DECAY, label="$r$ = 10")
a.set(xlabel="time $t$ [–]", ylabel="$X(t)$ (roll speed) [–]")
a.set_title("(a) one settles, one never does", fontsize=10)
a.legend(fontsize=8, loc="lower left")
dt = sol28["t"][1] - sol28["t"][0]                                 # sampling interval
freq = np.fft.rfftfreq(len(sol28["t"]) - 500, dt)                  # frequencies of the discrete Fourier transform (Ch. 7 P181)
for s, col, lab in ((sol28, C_GROW, "$r$ = 28"), (sol10, C_DECAY, "$r$ = 10")):
    sig = s["X"][500:] - np.mean(s["X"][500:])                     # drop the first 5 time units and the mean
    b.semilogy(freq, np.abs(np.fft.rfft(sig))/len(sig), color=col, label=lab)   # amplitude spectrum
b.set(xlabel="frequency [cycles per time unit]", ylabel="amplitude [–]", xlim=(0, 6), ylim=(1e-6, 10))
b.set_title("(b) broadband against one low bump", fontsize=10)
b.legend(fontsize=8)
plt.show()
f10 = ch10.dominant_frequency(sol10["t"][:1500], sol10["X"][:1500] - sol10["X"][-1])   # frequency of the decaying oscillation at r = 10
print(f"r = 10: X settles to {sol10['X'][-1]:.4f} after ringing at {f10:.3f} cycles per time unit "
      f"(eigenvalue of the steady roll: {abs(ch11.lorenz_eigs(10.0)[0].imag)/(2*np.pi):.3f})")
""",
    see="(a) The roll speed $X(t)$: for $r=10$ (teal) a ringing that dies away to a constant; for $r=28$ (rose) irregular swings "
        "that switch sign at unpredictable moments. (b) Their spectra: for $r=10$ one low bump near 1 cycle per time unit; for $r=28$ a level a hundred times higher at every frequency.",
    read="Chaos means energy at all frequencies — there is no period to find. The $r=10$ run has one frequency, the one of the "
         "linearised equations near the steady roll (printed); because that oscillation dies away, its line is a broad, low bump.",
    change="…$r=24$ (just below the threshold derived in D25): long irregular transients, then a steady state.",
    explain=r"""
1. `ch11.lorenz_integrate(s0, t_end, r=…)` integrates the three equations with `solve_ivp` (method DOP853, relative tolerance
   10⁻¹⁰) and returns `t`, `X`, `Y`, `Z`.
2. `np.fft.rfft` gives the amplitude at each frequency; `ch10.dominant_frequency` (Ch. 10) the main frequency of the ringing.""")
note("N113 [B]", r"""
**Phase space** with the pendulum: $\ddot X+(g/l)\sin X=0$ becomes the first-order system $\dot X=Y$, $\dot Y=-(g/l)\sin X$
(11.89); the plane $(X,Y)$ is its *phase space*, a solution is a *trajectory* there, and the number of independent initial values
(2) is its number of *degrees of freedom*. ($X$ is the angle, $Y$ the angular velocity, $g/l$ gravity over length.)""")
fig(r"""
starts = [(0.5, 0.0), (1.5, 0.0), (2.5, 0.0), (3.0, 0.0), (-3.14, 1.0), (3.14, -1.0)]   # (angle, angular velocity) at t = 0
tr = ch11.phase_portrait(ch11.pendulum_rhs, starts, t_end=20.0)    # one trajectory per start (g/l = 1, no damping)
Xg, Yg = np.meshgrid(np.linspace(-2*np.pi, 2*np.pi, 200), np.linspace(-3, 3, 120))   # a grid in the phase plane
fig, ax = plt.subplots(figsize=(7.2, 3.4))
ax.contour(Xg, Yg, 0.5*Yg**2 - np.cos(Xg), 12, colors=C_BASE, linewidths=0.5)   # lines of constant energy Y^2/2 - cos X
for p_ in tr:
    ax.plot(p_["X"], p_["Y"], color=C_SHEAR if abs(p_["X"]).max() > np.pi else C_BUOY, lw=1.6)   # orange: over the top; blue: swinging
ax.set(xlabel="angle $X$ [rad]", ylabel="angular velocity $Y$ [1/s]", xlim=(-2*np.pi, 2*np.pi), ylim=(-3, 3))
ax.set_title("Phase space of the pendulum", fontsize=10)
plt.show()
drift = max(np.max(np.abs(ch11.pendulum_energy(np.array([p_["X"], p_["Y"]])) - ch11.pendulum_energy(np.array([p_["X"][0], p_["Y"][0]])))) for p_ in tr)
print(f"largest change of the energy along any trajectory: {drift:.1e}")
""",
    see="The phase plane of a pendulum: closed loops around the origin (blue: swinging back and forth) and wavy lines above and "
        "below (orange: rotating over the top). Thin grey lines are curves of constant energy.",
    read="Each trajectory stays on one energy curve — the printed drift is of order 10⁻⁹. Without friction nothing attracts: the "
         "motion remembers its starting energy forever.",
    change="…we add damping (`damping=0.2`): every trajectory spirals into the point (0, 0) — an *attractor*.",
    explain=r"""
1. `ch11.phase_portrait(rhs, starts, t_end)` integrates a two-variable system from several starting points and returns one
   dictionary `t`, `X`, `Y` per start; `ch11.pendulum_rhs` is the right side of the pendulum system.
2. `ch11.pendulum_energy` evaluates $\tfrac12Y^2-(g/l)\cos X$, which should not change along a trajectory.""")
P("P274", "Jacobian matrix of a nonlinear ODE system and the stability of its fixed points", r"""
For $\dot{\mathbf s}=\mathbf f(\mathbf s)$ a *fixed point* has $\mathbf f(\mathbf s^*)=0$. Near it write
$\mathbf s=\mathbf s^*+\boldsymbol\delta$: to first order $\dot{\boldsymbol\delta}=J\boldsymbol\delta$ with the Jacobian
$J_{ij}=\partial f_i/\partial s_j$ at $\mathbf s^*$ (a first-order Taylor expansion in several variables, Ch. 3 P98). The fixed
point is stable if every eigenvalue of $J$ has negative real part (Ch. 9 P214); a complex pair crossing the imaginary axis starts
an oscillation (a *Hopf bifurcation*).""", code=r"""
import numpy as np                                    # arrays
Pr, r, b = 10.0, 28.0, 8/3                            # Lorenz's parameters
J0 = np.array([[-Pr, Pr, 0], [r, -1, 0], [0, 0, -b]]) # the Lorenz Jacobian at the origin X = Y = Z = 0
print(np.sort(np.linalg.eigvals(J0)))                 # -22.83, -2.67, 11.83: one positive -> the origin is unstable
""")
note("N114 [B]", r"""
**Attractors and bifurcations**: in a dissipative system trajectories crowd onto *attractors* — a fixed point (a steady flow) or
a *limit cycle* (a steady oscillation); as a parameter $R$ passes a critical value a fixed point can turn into a repeller and a
limit cycle grow around it whose size increases with the distance from the critical value (a *bifurcation*). Our model of that
picture is the Hopf normal form $\dot z=(\mu+\mathrm i\omega)z-\lvert z\rvert^2z$ for a complex $z=X+\mathrm iY$ (ours): stable
fixed point for $\mu<0$, limit cycle of radius $\sqrt\mu$ for $\mu>0$.""")
fig(r"""
fig, axs = plt.subplots(1, 3, figsize=(9.8, 3.2))
for ax, mu, st in ((axs[0], -0.5, [(0.9, 0.0)]), (axs[1], 0.5, [(0.1, 0.0), (1.2, 0.0)])):
    for h in ch11.phase_portrait(ch11.hopf_normal_form, st, t_end=30.0, n=1500, mu=mu):   # trajectories of the normal form
        ax.plot(h["X"], h["Y"], color=C_DECAY if mu < 0 else C_SHEAR, lw=1)
    ax.set(xlabel="$X$", ylabel="$Y$", xlim=(-1.3, 1.3), ylim=(-1.3, 1.3))
    ax.set_aspect("equal")
    ax.set_title(f"$\\mu$ = {mu}: " + ("fixed point" if mu < 0 else "limit cycle"), fontsize=10)
mus = np.linspace(-1, 1, 201)                                      # the control parameter
axs[2].plot(mus, ch11.limit_cycle_amplitude(mus), color=C_NEUT)
axs[2].plot(mus[mus > 0], 0*mus[mus > 0], "--", color=C_BASE)
axs[2].set(xlabel="$\\mu$", ylabel="amplitude of the attractor")
axs[2].set_title("bifurcation diagram", fontsize=10)
plt.show()
""",
    see="Left: for negative μ a trajectory spirals into the origin. Middle: for positive μ two trajectories — one starting inside, "
        "one outside — both wind onto the same circle. Right: the amplitude of the attractor against μ: zero up to μ = 0, then a "
        "square-root branch (the dashed line is the fixed point, now unstable).",
    read="Linear theory gives the growth at μ > 0; the nonlinear term fixes the final amplitude √μ. This is the gentle "
         "(supercritical) way a steady state gives way to an oscillation.",
    change="…the sign of the cubic term were reversed: no small stable cycle would exist near the threshold (a *subcritical* "
           "bifurcation) — that is Lorenz's case (D25).",
    explain=r"""
`ch11.hopf_normal_form` is the right side of the normal form in the two real variables $(X,Y)$; `ch11.limit_cycle_amplitude(mu)`
returns $\sqrt\mu$ for $\mu>0$ and 0 otherwise.""")
note("N115 [B]", r"""
**Lorenz's truncation**: stress-free walls, rolls along $y$, $u=-\partial\psi/\partial z$, $w=\partial\psi/\partial x$ (⚠️ slip #9:
the opposite sign to §11.7), and $\psi\propto X(t)\cos(\pi z)\sin(kx)$, $T'\propto Y(t)\cos(\pi z)\cos(kx)+Z(t)\sin(2\pi z)$
(11.90), $z\in[-\tfrac12,\tfrac12]$: $X$ measures the roll speed, $Y$ the temperature difference between rising and sinking
fluid, $Z$ the distortion of the mean temperature profile. $\cos\pi z$ is the free–free mode $\sin\pi(z+\tfrac12)$ of C05 (slip
#2's correct family).""")
P("P275", "Galerkin truncation: keep a few modes and project with orthogonality of sines", r"""
Write the unknown as a sum of a few fixed shapes with time-dependent amplitudes, insert it, and multiply the equation by each shape
and integrate: orthogonality ($\int\sin m\pi z\,\sin n\pi z\,dz=0$ for $m\neq n$, Ch. 6 P151) keeps one equation per amplitude, and
everything that does not fit the chosen shapes is thrown away. Ch. 10's finite elements used hats as shapes; here we use sines and
cosines.""", code=r"""
import numpy as np                                    # arrays
z = np.linspace(-0.5, 0.5, 20001)                     # across the layer
print(np.trapezoid(np.cos(np.pi*z)*np.cos(3*np.pi*z), z), np.trapezoid(np.cos(np.pi*z)**2, z))   # ~0 and 0.5: orthogonal shapes
""")
D("D24", ref="11.91", check_src=PF["D24"]["check_src"] + "\n" + D24_EXTRA.strip("\n"))
P("P276", "purely imaginary roots of a cubic λ³ + a₂λ² + a₁λ + a₀: exactly when a₂a₁ = a₀", r"""
Put $\lambda=\mathrm i\omega$: the real part gives $-a_2\omega^2+a_0=0$ and the imaginary part $-\omega^3+a_1\omega=0$, so
$\omega^2=a_1$ and $a_2a_1=a_0$. For positive coefficients, all three roots have negative real part when $a_2a_1>a_0$ (the
Routh–Hurwitz condition); at $a_2a_1=a_0$ a pair sits on the imaginary axis — the Hopf point.""", code=r"""
import numpy as np                                    # arrays
a2, a1, a0 = 13.6667, 101.3333, 1440.0                # the cubic of the Lorenz convection state at r = 28 (D25 step 9)
print(a2*a1, a0)                                      # 1384.9 < 1440: a pair has crossed the axis (unstable)
print(np.roots([1, a2, a1, a0]))                      # -13.85 and 0.094 + 10.19j, 0.094 - 10.19j
""")
D("D25", ref="11.91")
nb.worked_example("the Lorenz numbers at r = 28", r"""
$\Pr=10$, $b=8/3$, $r=28$ (Lorenz 1963).

1. Convection states: $\bar X=\bar Y=\pm\sqrt{b(r-1)}=\pm\sqrt{\tfrac83\times27}=\pm\sqrt{72}=\pm8.485$, $\bar Z=r-1=27$.
2. Hopf point: $r_H=\Pr(\Pr+b+3)/(\Pr-b-1)=10\times15.667/6.333=24.74<28$ ⇒ the convection states are unstable (eigenvalues
   $0.094\pm10.19\,\mathrm i$).
3. Origin: eigenvalues $-8/3$ and the roots of $\lambda^2+11\lambda-270=0$: $11.83$ and $-22.83$ ⇒ unstable.
4. Volume: the divergence of the right-hand sides, $\frac{\partial\dot X}{\partial X}+\frac{\partial\dot Y}{\partial Y}+\frac{\partial\dot Z}{\partial Z}=-\Pr-1-b$
   (the trace of the Jacobian, P274), is the rate at which a small volume of states changes: every blob shrinks at
   $-(\Pr+1+b)=-13.67$ per unit time, so the attractor has zero volume — yet it
   is not a point or a loop: it is *strange*.
5. $b=8/3$ comes from $k^2=\pi^2/2$, the free–free critical wavenumber of C05: $4\pi^2/(\pi^2+\pi^2/2)=8/3$.""")
code(r"""
print("fixed points at r = 28:", [tuple(round(v, 4) for v in p_) for p_ in ch11.lorenz_fixed_points(28.0)])     # origin and C+, C-
print("eigenvalues at C+:", np.round(ch11.lorenz_eigs(28.0, which="C"), 4))            # a complex pair with positive real part
print("eigenvalues at the origin:", np.round(ch11.lorenz_eigs(28.0, which="origin"), 4))
print("eigenvalues at C+ for r = 20 (below r_H):", np.round(ch11.lorenz_eigs(20.0, which="C"), 3))   # every real part negative: steady rolls are stable
print("largest real part at C+ exactly at r_H:", f"{np.max(ch11.lorenz_eigs(ch11.lorenz_hopf_r(), which='C').real):.1e}")   # ~0: the pair sits on the axis
print("Hopf point r_H =", round(ch11.lorenz_hopf_r(), 4), "| b at k = pi/sqrt(2):", round(ch11.lorenz_b(np.pi/np.sqrt(2)), 4),
      "| volume contraction rate:", round(ch11.lorenz_divergence(), 3))
sep = ch11.lorenz_separation(delta0=1e-8)                           # two runs that start 1e-8 apart in X
print(f"separation grows like e^(lambda t) with fitted slope lambda = {sep['slope']:.2f} over t in {sep['window']}")
print("time at which the two runs differ by 1:", ch11.lorenz_predictability_time(1e-8, 1.0))
""", explain=r"""
1. `ch11.lorenz_fixed_points(r)` returns the steady states of D25; `ch11.lorenz_eigs(r, which=…)` the eigenvalues of the Jacobian
   (P274) there: at $r=28$ both the origin and $C_\pm$ are unstable; at $r=20$ the rolls are stable ($-13.36$, $-0.155\pm8.71\,\mathrm i$);
   exactly at $r_H$ the real part of the pair is zero — the three numbers D25's *Check it* quotes.
2. `ch11.lorenz_hopf_r()` is D25's threshold 24.74; `ch11.lorenz_divergence()` the rate at which volumes of states shrink.
3. `ch11.lorenz_separation(delta0)` integrates two runs started `delta0` apart and fits the slope of $\ln$(separation) — about
   0.9 here (the literature value of the largest Lyapunov exponent is about 0.9; our single-run fit depends on the window and is
   qualitative). The two runs differ by 1 after about 29 time units.""")
scratch(r"""
def lorenz_f(s, Pr=10.0, r=28.0, b=8/3):                            # the right side of (11.91)
    X, Y, Z = s
    return np.array([Pr*(Y - X), -X*Z + r*X - Y, X*Y - b*Z])

s, dt_ = np.array([1.0, 1.0, 1.0]), 0.01                            # start and time step
traj = [s]                                                          # the trajectory, step by step
for n in range(4000):                                               # classical fourth-order Runge-Kutta (Ch. 3 P95)
    k1 = lorenz_f(s); k2 = lorenz_f(s + 0.5*dt_*k1); k3 = lorenz_f(s + 0.5*dt_*k2); k4 = lorenz_f(s + dt_*k3)   # four slopes
    s = s + dt_*(k1 + 2*k2 + 2*k3 + k4)/6                           # weighted average of the slopes
    traj.append(s)
traj = np.array(traj)                                               # shape (4001, 3): the same times as sol28["t"]
diff = np.abs(traj[:, 0] - sol28["X"])                              # difference in X between our RK4 and the library's adaptive solver
assert np.allclose(traj[:501, 0], sol28["X"][:501], atol=1e-3)      # they agree up to t = 5
assert np.allclose(ch11.lorenz_rk4((1, 1, 1), 0.01, 1)[1], traj[1])  # the library's fixed-step RK4 takes the same first step
print(f"largest difference up to t = 5: {diff[:501].max():.1e}; first time the two differ by more than 1: t = {sol28['t'][np.argmax(diff > 1)]:.1f}")
""", r"""
1. A hand-written RK4 loop for the three equations.
2. It agrees with the library's high-accuracy solver to 4 × 10⁻⁴ up to $t=5$ — and then drifts away: by $t\approx18$ the two
   differ by more than 1. Both are correct numerical solutions; their tiny truncation differences grow like any other
   perturbation. **The disagreement is the lesson.**""")
P("P277", "Lyapunov exponent: the slope of log(separation) against time", r"""
If two nearby trajectories separate like $\lvert\delta(t)\rvert\approx\lvert\delta_0\rvert e^{\lambda t}$, then
$\ln\lvert\delta\rvert$ grows linearly with slope $\lambda$, the (largest) *Lyapunov exponent*; $\lambda>0$ means chaos. The time to
grow from $\delta_0$ to an error of size $\Delta$ is about $\ln(\Delta/\delta_0)/\lambda$ — halving the initial error only adds
$\ln2/\lambda$ to the forecast horizon. Fit the slope only while the separation is small (it saturates at the attractor's size).
`np.polyfit(t, y, 1)` fits a straight line by least squares and returns (slope, intercept).""", code=r"""
import numpy as np                                    # arrays
t = np.linspace(0, 10, 11)                            # times
sep_ideal = 1e-8*np.exp(0.9*t)                        # an ideal exponential separation with lambda = 0.9
print(np.polyfit(t, np.log(sep_ideal), 1)[0])         # 0.9: the slope of ln(separation) is lambda
""")
P("P278", "matplotlib 3-D line plots", r"""
`ax = fig.add_subplot(projection="3d")` makes 3-D axes; `ax.plot(X, Y, Z)` draws a curve through space, `ax.view_init(elev, azim)`
sets the viewpoint. (Ch. 1–2's 3-D figures were plotly, P41/P64; for a static page figure matplotlib is lighter.)""", code=r"""
import numpy as np, matplotlib.pyplot as plt           # arrays and figures
fig = plt.figure(figsize=(3.6, 2.8))                   # an empty figure
ax = fig.add_subplot(projection="3d")                  # 3-D axes
t = np.linspace(0, 6*np.pi, 300)                       # a parameter along the curve
ax.plot(np.cos(t), np.sin(t), t/10)                    # a helix
ax.view_init(20, 35)                                   # elevation 20 degrees, azimuth 35 degrees
plt.show()
""")
fig(r"""
fig = plt.figure(figsize=(10.2, 3.6))
ax3 = fig.add_subplot(1, 3, 1, projection="3d")                    # panel (a): the orbit in (X, Y, Z)
ax3.plot(sol28["X"], sol28["Y"], sol28["Z"], color=C_BASE, lw=0.4)
ax3.plot(sol28["X"][-500:], sol28["Y"][-500:], sol28["Z"][-500:], color=C_SHEAR, lw=1.4)   # the last 5 time units
for p_ in ch11.lorenz_fixed_points(28.0):
    ax3.plot([p_[0]], [p_[1]], [p_[2]], "o", color=COLORS["ink"], ms=4)
ax3.view_init(18, -60)
ax3.set(xlabel="$X$", ylabel="$Y$")                                 # the vertical axis is Z (named in the title: its label would collide with panel b)
ax3.set_title("(a) the attractor ($Z$ upward)", fontsize=10)
a = fig.add_subplot(1, 3, 2)
a.plot(sep["t"], sep["a"][0], color=C_SHEAR, lw=0.9, label="run 1")
a.plot(sep["t"], sep["b"][0], color=C_BUOY, lw=0.9, label="run 2 (started $10^{-8}$ away)")
a.set(xlabel="time $t$", ylabel="$X(t)$")
a.set_title("(b) two runs that part", fontsize=10)
a.legend(fontsize=7, loc="lower left")
b = fig.add_subplot(1, 3, 3)
b.semilogy(sep["t"][1:], sep["sep"][1:], color=C_BASE)
w0, w1 = sep["window"]                                             # the window used for the fit
tw = np.array([w0, w1])
i0 = int(np.argmin(np.abs(sep["t"] - w0)))                         # index of the window's start
b.semilogy(tw, sep["sep"][i0]*np.exp(sep["slope"]*(tw - w0)), color=C_GROW, lw=2.5, label=f"slope {sep['slope']:.2f}")
b.set(xlabel="time $t$", ylabel="separation $|\\delta|$")
b.set_title("(c) exponential separation", fontsize=10)
b.legend(fontsize=8)
savefig(fig, "ch11", "c15_lorenz")
plt.show()
""",
    see="(a) The orbit in the three-dimensional phase space: two lobes, each around one convection state (dots); the last five time "
        "units are drawn bold. (b) $X(t)$ for two runs that start 10⁻⁸ apart: on top of each other until about $t=28$, then "
        "unrelated (they differ by 1 at $t=29.1$, printed above). (c) Their separation on a logarithmic axis: flat for the first twelve "
        "time units, then a straight rising line that flattens once it is as large as the "
        "attractor.",
    read="$X$ changes sign irregularly — the roll reverses its sense of rotation — each time the orbit jumps to the other lobe. "
         "The slope in (c) is the Lyapunov exponent (P277), about 0.9 here. The flat start is not a contradiction: for the first "
         "twelve time units both runs are still spiralling regularly around one convection state, where neighbouring states "
         "do not separate on average; the exponential growth begins once the orbit starts switching lobes.",
    change="…the first difference were 10⁻¹⁶ (the round-off of double precision): the two runs would agree only about 20 time "
           "units longer — $\\ln(10^8)/\\lambda$. A hundred million times more precision buys two thirds more time.",
    explain=r"""
1. `sep["a"]` and `sep["b"]` are the two runs (rows $X$, $Y$, $Z$); `sep["sep"]` their distance; `sep["window"]` the time window
   of the fit and `sep["slope"]` the fitted exponent.
2. `fig.add_subplot(1, 3, 1, projection="3d")` (P278) puts a 3-D panel next to two ordinary ones.""")
nb.md(r"""
**Animation A5 — two runs, one roll.** The two runs of panel (b) on the $X$–$Z$ plane (with a short trail), the temperature field
of the convection roll driven by run 1 ($T'$ of Lorenz's truncation), and the separation.""")
nb.animation(r"""
nf = 14 if FAST else 18                                            # number of frames
idx = np.linspace(2200, 3700, nf).astype(int)                      # sample times t = 22 ... 37, where the runs part
xr, zr = np.meshgrid(np.linspace(0, 2*np.sqrt(2), 40), np.linspace(-0.5, 0.5, 20))   # one wavelength 2 pi/k = 2 sqrt(2) of the roll
fld = lambda j: ch11.lorenz_fields(xr, zr, sep["a"][0, j], sep["a"][1, j], sep["a"][2, j])   # roll fields for the state of run 1 at index j
fig, axs = plt.subplots(1, 3, figsize=(8.4, 2.6))
axs[0].plot(sol28["X"][::3], sol28["Z"][::3], color=C_BASE, lw=0.3)  # the attractor as a faint background (every third point)
(tr1,) = axs[0].plot([], [], color=C_SHEAR, lw=1.5)                # trail of run 1
(tr2,) = axs[0].plot([], [], color=C_BUOY, lw=1.5)                 # trail of run 2
(d1,) = axs[0].plot([], [], "o", color=C_SHEAR)
(d2,) = axs[0].plot([], [], "o", color=C_BUOY)
axs[0].set(xlabel="$X$", ylabel="$Z$")
mesh = axs[1].pcolormesh(xr, zr, fld(idx[0])["T"], cmap="RdBu_r", vmin=-45, vmax=45, shading="auto")   # temperature disturbance of the roll
axs[1].set(xlabel="$x$ [$d$]", ylabel="$z$ [$d$]")
ttl = axs[1].set_title("", fontsize=9)
axs[2].semilogy(sep["t"][1::4], sep["sep"][1::4], color=C_BASE)
(mark,) = axs[2].semilogy([sep["t"][idx[0]]], [sep["sep"][idx[0]]], "o", color=C_GROW)
axs[2].set(xlabel="time $t$", ylabel="separation")

def update(i):                                                     # frame i: time index idx[i]
    j = idx[i]
    tr1.set_data(sep["a"][0, j - 150:j + 1], sep["a"][2, j - 150:j + 1])   # the last 1.5 time units of run 1
    tr2.set_data(sep["b"][0, j - 150:j + 1], sep["b"][2, j - 150:j + 1])   # ... and of run 2
    d1.set_data([sep["a"][0, j]], [sep["a"][2, j]])
    d2.set_data([sep["b"][0, j]], [sep["b"][2, j]])
    mesh.set_array(fld(j)["T"].ravel())
    ttl.set_text(f"t = {sep['t'][j]:.1f}: roll turns {'one way' if sep['a'][0, j] > 0 else 'the other way'}")
    mark.set_data([sep["t"][j]], [sep["sep"][j]])
    return tr1, tr2, d1, d2, mesh, ttl, mark

show_animation(animate(update, frames=nf, fig=fig, interval=170), player="video", dpi=60)
""", explain=r"""
1. `ch11.lorenz_fields(x, z, X, Y, Z)` turns a state $(X,Y,Z)$ back into the fields of the truncation: stream function, temperature
   disturbance and velocities of the roll.
2. Each frame moves the two dots (and their trails), repaints the roll's temperature for run 1 and moves the marker on the
   separation curve.""")
nb.figure_notes(
    see="Two dots glued together on the attractor, then flying apart onto different lobes; in the middle the roll's temperature "
        "pattern, which flips left–right whenever the orange dot changes lobe; on the right the separation climbing to order 10.",
    read="A lobe switch *is* a reversal of the roll. Once the separation has saturated the two runs are two different weathers — "
         "both perfectly valid solutions of the same three equations.",
    change="…$r=20$: both dots spiral into the same convection state — no chaos, only a transient.")
nb.md(r"""
**Slider figure F8 — from conduction to chaos.** The $X$–$Z$ trajectory after a transient for a range of $r$: it collapses onto the
origin for $r<1$, onto a convection state $C_\pm$ for $1<r<r_H\approx24.74$ (after long irregular transients just below $r_H$), and
never settles beyond.""")
nb.plotly(r"""
sw = ch11.lorenz_r_sweep() if FAST else ch11.lorenz_r_sweep(np.array([0.5, 2.0, 5.0, 10.0, 15.0, 20.0, 22.0, 24.0, 25.0, 26.0, 28.0, 30.0]))   # X(t), Z(t) for each r
late = sw["t"] > 12.0                                              # the part of each run after the first 12 time units

def f8(r):                                                         # traces for one value of r
    i = int(np.argmin(np.abs(sw["r"] - r)))                        # the stored run for this r
    cpm = np.sqrt(8/3*(r - 1)) if r > 1 else np.nan                # the convection states C+ and C- (D25 step 4)
    return {"trajectory after t = 12": (sw["X"][i][late], sw["Z"][i][late]),                             # where the run ends up
            "start of the run": (sw["X"][i][:150], sw["Z"][i][:150]),                                    # its first 3 time units
            "steady states C+ and C−": (np.array([-cpm, cpm]), np.array([r - 1, r - 1]))}                # the two convection states

figF8 = slider_figure(f8, "r", sw["r"], xlabel="X (roll speed)", ylabel="Z (profile distortion)",       # slider: the heating r
                      title="Lorenz system: conduction (r < 1), steady rolls (1 < r < 24.74), chaos beyond",   # the message
                      xrange=(-22, 22), yrange=(-2, 50), modes={"steady states C+ and C−": "markers"})   # fixed window; fixed points as markers
figF8.show()
""", explain=r"""
1. `ch11.lorenz_r_sweep(r_values)` integrates the system from (1, 1, 1) for each $r$ and returns `X`, `Z` as arrays (one row per
   $r$).
2. Drag $r$: for $r\le24$ the late trajectory is a dot on one of the markers; from 25 on it is a butterfly that never closes.""")
nb.live(r"""
def lorenz_live(r=28.0, Pr=10.0, b=2.667, log10_delta0=-8.0):      # free parameters and the size of the first difference
    sp_ = ch11.lorenz_separation(delta0=10.0**log10_delta0, t_end=40.0, Pr=Pr, r=r, b=b)   # two runs, delta0 apart
    fig, (a, c) = plt.subplots(1, 2, figsize=(8.4, 3.0))
    a.plot(sp_["a"][0], sp_["a"][2], color=C_SHEAR, lw=0.5)        # run 1 on the X-Z plane
    a.plot(sp_["b"][0], sp_["b"][2], color=C_BUOY, lw=0.5, alpha=0.6)   # run 2
    a.set(xlabel="$X$", ylabel="$Z$")
    c.semilogy(sp_["t"][1:], sp_["sep"][1:] + 1e-300, color=C_BASE)   # their separation
    c.set(xlabel="time $t$", ylabel="separation")
    c.set_title(f"r_H = {ch11.lorenz_hopf_r(Pr, b):.2f}" if Pr > b + 1 else "no Hopf point (Pr < b + 1)", fontsize=9)
    plt.show()

live(lorenz_live, r=(0.5, 40.0, 0.5), Pr=(1.0, 20.0, 0.5), b=(0.5, 4.0, 0.1), log10_delta0=(-14.0, -2.0, 1.0))
""", explain=r"""
`lorenz_live` reruns the pair of trajectories for your $r$, $\Pr$, $b$ and first difference, and prints the Hopf threshold for
that $\Pr$ and $b$. Slider figure F8 above carries the $r$-dependence on the web page.""")
P("P279", "iterated maps: fixed point, stability ∣f′(x*)∣ < 1, cobweb diagram", r"""
A map $x_{n+1}=f(x_n)$ is a rule applied again and again. A fixed point has $x^*=f(x^*)$; a small error $\delta$ becomes
$f'(x^*)\delta$ at the next step, so $x^*$ attracts if $\lvert f'(x^*)\rvert<1$ (compare Ch. 10's amplification factor
$\lvert G\rvert\le1$). A *cobweb* draws the iteration: up to the curve $y=f(x)$, across to the diagonal $y=x$, repeat.""", code=r"""
A, x = 2.8, 0.2                                       # the logistic map x -> A x (1 - x), started at 0.2
for n in range(30):                                   # iterate 30 times
    x = A*x*(1 - x)
print(round(x, 4), round(1 - 1/A, 4), abs(2 - A))     # x has almost reached the fixed point 1 - 1/A = 0.6429; |f'| = 0.8 < 1: stable
""")
note("N117 [B] · N129 [B]", r"""
**One route to chaos: period doubling.** The logistic map $x_{n+1}=Ax_n(1-x_n)$ (Exercise 11.14, a toy model of transition; its
steady state is $x^*=1-1/A$) has a stable fixed point for $1<A<3$ ($\lvert f'(x^*)\rvert=\lvert2-A\rvert<1$), a stable 2-cycle from
$A_1=3$, a 4-cycle from $A_2=1+\sqrt6=3.4495$, and so on; the gaps between doublings shrink by a universal ratio,
$\frac{R_n-R_{n-1}}{R_{n+1}-R_n}\to4.6692$ (Feigenbaum 1978; the book writes $R_n$ for the bifurcation values, here $A_n$), the
same in any one-hump map and in convection experiments.""")
code(r"""
pd_ = ch11.period_doubling_points(6)                                # where the 2-, 4-, 8-, ... cycles are born (A_n), and superstable values (S_n)
print("period doublings at A =", np.round(pd_["A_n"], 6))
print("estimates of Feigenbaum's ratio:", np.round(ch11.feigenbaum_estimate(6), 4), "(published", BENCH["feigenbaum"]["delta"], ")")
""", explain=r"""
`ch11.period_doubling_points` returns the bifurcation values $A_n$ = 3, 3.449490, 3.544090, 3.564407, … and the *superstable*
parameters (where $x=\tfrac12$ is on the cycle), from which `ch11.feigenbaum_estimate` forms successive ratios: 4.709, 4.681,
4.663, 4.668, 4.669 → 4.6692.""")
fig(r"""
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.6))
xl = np.linspace(0, 1, 200)                                        # x from 0 to 1
a.plot(xl, xl, color=C_BASE, lw=0.8)                               # the diagonal y = x
for A, col in ((2.8, C_DECAY), (3.2, C_NEUT), (3.5, C_GROW)):
    a.plot(xl, A*xl*(1 - xl), color=col, lw=0.8)                   # the map y = A x (1 - x)
    cw = ch11.cobweb(A, 0.2, 60)                                   # the cobweb path of 60 iterations from x = 0.2
    if A < 3:                                                      # a stable fixed point: the late cobweb has shrunk to a single point
        a.plot(1 - 1/A, 1 - 1/A, "o", color=col, ms=7)             # mark it: x* = 1 - 1/A on the diagonal
    a.plot(cw["xs"][60:], cw["ys"][60:], color=col, lw=1.0, label=f"$A$ = {A}")   # only the late part: the attractor
a.set(xlabel="$x_n$", ylabel="$x_{n+1}$", xlim=(0, 1), ylim=(0, 1))
a.set_title("(a) cobwebs: point, 2-cycle, 4-cycle", fontsize=10)
a.legend(fontsize=8)
bd = ch11.bifurcation_diagram(np.linspace(2.8, 4.0, 500 if FAST else 1000))   # the values visited after a transient, for each A
b.plot(bd["A"], bd["x"], ",", color=COLORS["ink"], alpha=0.25)
for An in pd_["A_n"][:3]:
    b.axvline(An, color=C_NEUT, lw=0.8, ls=":")
b.set(xlabel="$A$", ylabel="$x$ on the attractor", xlim=(2.8, 4.0))
b.set_title("(b) the period-doubling tree", fontsize=10)
plt.show()
""",
    see="(a) The logistic map for three values of $A$ with the late part of the cobweb: a single point (teal), a square that "
        "visits two values (purple), a figure that visits four (rose). (b) All values visited after a transient, against $A$: one "
        "line splitting into 2, 4, 8 … and then a cloud; dotted lines mark the first three doublings.",
    read="Each split comes sooner than the last by the factor 4.669: the intervals between doublings shrink geometrically, so "
         "infinitely many fit before $A\\approx3.57$, where chaos starts.",
    change="…a different one-hump map (for example $x\\to A\\sin\\pi x$): different numbers $A_n$, the same ratio.",
    explain=r"""
1. `ch11.cobweb(A, x0, n)` returns the corner points of the cobweb path; only its late part is drawn, so the transient is hidden.
2. `ch11.bifurcation_diagram(A_values)` iterates the map for each $A$, discards a transient and returns the values visited.""")
note("N118 [C]", r"""
**Other routes** (named): quasi-periodicity — after two incommensurate frequencies the motion becomes chaotic (Ruelle & Takens
1971), unlike Landau's picture of an endless sequence of new frequencies; the Lorenz system follows neither — its chaos appears
after a subcritical Hopf point (D25). Pointer: the turbulence-onset debate in Ch. 12.""")
note("N119 [C]", r"""
**What chaos means for prediction** (named): a deterministic system can still be unpredictable in practice — not because of
quantum uncertainty but because no measurement is exact (Poincaré saw this in 1908). Climate hook: weather forecasts lose skill
after about two weeks; centres therefore run *ensembles* of forecasts from slightly different starts and forecast probabilities;
climate projections ask about statistics of the attractor, not one trajectory (Ch. 13). Whether turbulence itself is "chaos in a
few modes" remains open.""")
explainer("lorenz_attractor", "How can three exact equations be unpredictable?", r"""
The 3-D orbit, the $X(t)$ traces and the log-separation panel run on one clock from two starts 10⁻⁸ apart; the $r$ slider moves
the fixed points and changes the verdict while you watch; a static attractor picture hides both the divergence and the
$r$-dependence.""", [
    "Preset 'r = 10': both runs spiral into a convection state; raise r past 24.74 and watch the status turn rose.",
    "At r = 28 read the time at which the runs differ by 50 % (end card); make δ₀ 100 times smaller and see how few time units you gain.",
    "Click the orbit: the inspector shows the Jacobian's eigenvalues there.",
])
nb.pointer(r"""**S02** — Literature: the benchmarks used here are cited where they appear (Chandrasekhar 1961; Orszag 1971; Thomas,
via Gallagher, Griffiths & Stephen 2016; Jordinson 1970; Tatsumi & Kakutani 1958; Michalke 1964; Drazin & Reid 1981; Lorenz 1963;
Feigenbaum 1978); the book's bibliography is not reproduced.""")

nb.summary(
    clicked=[
        r"**C01** A small disturbance is a sum of waves $e^{\mathrm ikx+\sigma t}$; a flow is stable only if $\sigma_r\le0$ for every $k$, and onset happens where $\sigma_r(k)$ first touches zero.",
        "**C02** Across a velocity jump, gravity's restoring term fights shear's kinetic term under a square root; when shear wins the two wave speeds become a complex pair and one grows — short waves first, unless surface tension holds them (6.7 m/s for wind over water).",
        "**C03** Heating from below gives two ODEs in $z$ with one parameter Ra; $\\sigma$ is a real eigenvalue, computed by Chebyshev collocation.",
        "**C04** Each cell width has its own marginal Ra; the valley's bottom, $\\mathrm{Ra}_c=1707.76$ at $K_c=3.117$, is the onset, with cells about as wide as the layer is deep.",
        "**C05** With stress-free walls sine modes make everything algebra: $\\mathrm{Ra}=(\\pi^2+K^2)^3/K^2$, minimum $27\\pi^4/4$ at $K^2=\\pi^2/2$.",
        "**C06** Because heat diffuses about 100× faster than salt, a column lighter on top can still overturn into salt fingers when $\\mathrm{Rs}-\\mathrm{Ra}>27\\pi^4/4$.",
        "**C07** Swapping rings releases energy when the squared circulation falls outward; the centrifugal force plays gravity and viscosity adds the 1708-type threshold $\\mathrm{Ta}_c\\approx1708/(\\tfrac12(1+\\Omega_2/\\Omega_1))$.",
        "**C08** Stratified shear flow obeys one ODE, Taylor–Goldstein, with $c$ as the eigenvalue; its unstable modes come in conjugate pairs.",
        "**C09** $\\mathrm{Ri}>\\tfrac14$ everywhere guarantees stability (Miles–Howard); $\\mathrm{Ri}<\\tfrac14$ only allows instability.",
        "**C10** Every unstable wave speed lies in the half-disc on $[U_{\\min},U_{\\max}]$; growth is at most $k(U_{\\max}-U_{\\min})/2$.",
        "**C11** Squire makes 2-D waves the first to go unstable; eliminating pressure leaves the fourth-order Orr–Sommerfeld equation.",
        "**C12** Inviscid instability needs an inflection point with a vorticity maximum; neutral waves have a critical layer and cat's eyes.",
        "**C13** Plane Poiseuille flow, inviscidly stable, is viscously unstable above $\\mathrm{Re}_c=5772.22$ ($k_c=1.02056$) — Tollmien–Schlichting waves.",
        "**C14** A disturbance grows when its Reynolds-stress production $-\\int uvU'$ beats viscous dissipation; viscosity can supply the phase shift that makes production positive.",
        "**C15** Three modes of convection give the Lorenz system; past $r_H\\approx24.74$ it is chaotic — deterministic yet unpredictable beyond about 20–30 time units.",
    ],
    feeds_forward=[
        "**Ch. 12:** Reynolds stresses $-\\langle uv\\rangle$ and the production term of the turbulent kinetic energy (C14), transition routes (N111), Richardson-number mixing (C09).",
        "**Ch. 13:** baroclinic instability (C01's recipe with rotation; the semicircle and inflection theorems reappear with $U''-\\beta$, Rayleigh–Kuo), inertial instability (C07), stratified shear (C08–C09), predictability (C15).",
        "**Ch. 14:** boundary-layer transition on airfoils (C13).",
    ],
    left_out=[
        "Exercises 11.3–11.5 (S01).",
        "The stability analysis of pipe flow (named, N101).",
        "Spatial instability and absolute/convective instability (named, N03).",
        "Weakly nonlinear theory (named, N110).",
    ])
# Comment pass: plotting boilerplate lines that carry no comment get a plain-words one (the science lines are commented by hand).
# ---------------------------------------------------------------------------------------------------------------------
_ANNOT_RULES = [
    (r"plt\.subplots\(", "the figure and its panels"),
    (r"plt\.show\(\)", "display the figure"),
    (r"savefig\(", "keep a copy in outputs/ch11 and display"),
    (r"\bfig\.suptitle\(", "the figure's message as its title"),
    (r"set_title", "the panel's title: what this panel shows"),
    (r"set_xlabel|set_ylabel", "axis labels, with units where the quantity has them"),
    (r"set_xlim|set_ylim|set_xticks|set_yscale", "axis range / scale chosen so the feature discussed below is visible"),
    (r"set_aspect", "equal scales in x and y, so shapes are not distorted"),
    (r"\.legend\(", "legend: one entry per curve"),
    (r"\.axhline\(|\.axvline\(|\.axvspan\(", "a reference line or band"),
    (r"\.text\(|\.annotate\(", "a label written on the plot"),
    (r"fill_between|fill_betweenx", "shade the area between two curves"),
    (r"\.contourf\(|\.tricontourf\(", "filled colour contours of the field"),
    (r"\.contour\(", "contour lines of the field"),
    (r"\.tripcolor\(|\.pcolormesh\(|\.imshow\(", "a colour map of the field"),
    (r"\.quiver\(|\.streamplot\(", "arrows / streamlines of the velocity"),
    (r"\.loglog\(|\.semilogy\(|\.semilogx\(", "curve(s) on logarithmic axes"),
    (r"\.bar\(", "bars"),
    (r"\.plot\(|\.triplot\(|\.spy\(", "draw the curve(s)"),
    (r"fig\.show\(\)", "display the interactive figure"),
    (r"^\s*for .* in (?:range|enumerate\(range)", "repeat the indented lines for each index"),
    (r"^\s*for .* in ", "repeat the indented lines for each item listed"),
    (r"slider_figure\(", "one figure whose curves follow the slider; every position precomputed (P17)"),
    (r"^\s*import numpy", "arrays and maths"),
    (r"^\s*import sympy", "symbolic algebra"),
    (r"^\s*from |^\s*import ", "a tool used below"),
    (r"^\s*def ", "a small helper used below"),
    (r"^\s*return ", "the values for this call"),
    (r"sp\.symbols|sp\.Function", "symbols for sympy"),
    (r"np\.linspace|np\.geomspace|np\.logspace|np\.arange", "sample points"),
    (r"\.append\(", "store the value"),
    (r"^\s*assert ", "check: stops with an error if it fails"),
    (r"^\s*print\(", "show the numbers computed above"),
]


def annotate_boilerplate() -> int:
    changed = 0
    for c in nb.cells:
        if c.cell_type != "code" or "setup" in c.metadata.get("tags", []):
            continue
        lines = c.source.split("\n")
        new = []
        depth_str = False
        for ln in lines:
            if ln.count('"""') % 2 == 1:
                depth_str = not depth_str
            if depth_str or "#" in ln or not ln.strip() or ln.rstrip().endswith("\\") or ln.rstrip().endswith("(") \
                    or ln.rstrip().endswith(","):
                new.append(ln)
                continue
            for pat, msg in _ANNOT_RULES:
                if re.search(pat, ln):
                    ln = ln.rstrip() + "   # " + msg
                    changed += 1
                    break
            new.append(ln)
        c.source = "\n".join(new)
    return changed



# ---------------------------------------------------------------------------------------------------------------------
# Hand-written comments for the plot and print lines (lesson review): each says what the curve, marker or printed number
# IS in the physics.  Format: start of the code line ||| comment.  Applied to lines that carry no comment yet; every
# entry must match at least one line (the build fails otherwise, so the table cannot go stale silently).
# ---------------------------------------------------------------------------------------------------------------------
HAND_COMMENTS = r'''
print(len([n for n in dir(ch11) ||| how many tested functions the chapter module offers
(dots[(s, x0)],) = ax.plot( ||| the ball at its starting position on the landscape
print(f"k = {k:.0f}: amplitude at t = 10 s ||| each mode alone has grown or decayed by exactly e^(sigma t)
print(f"sum of all three: max ||| the sum is dominated by the one growing mode
print("from-scratch normal mode ||| all three ways of writing the wave give the same numbers
grow = np.array([ch11.normal_mode_growth("interface", ki, rho1=1000.0, rho2=1.2, surface_tension=0.074).real ||| growth rate of each wavenumber: water over air, held by surface tension at short waves
a.plot(k, grow0, "--", color=C_BASE ||| dashed: without surface tension every shorter wave grows faster
a.plot(k, grow, color=C_GROW ||| with surface tension: only waves longer than the cut-off grow
print(f"fastest of the three ||| the mode at the top of the growth curve wins
b.plot(xs, zeta, color=COLORS["ink"]) ||| the interface itself
print(f"wave speed of this mode c ||| a real wave speed: this wave neither grows nor decays
print(f"k_c = {k_c:.0f} 1/m ||| the wavelength below which shear beats gravity
print("U1 = 4 m/s: c =" ||| a weaker wind: two real wave speeds, no growth
print(f"growth is zero up to k ||| below k_c gravity holds the interface; above it growth rises without bound
print("quadratic by np.roots:" ||| our own roots, and how well the mode satisfies the governing equations
a.semilogx(k, dU0, "--", color=C_SHEAR ||| dashed: smallest wind that makes wavenumber k grow if only gravity resists
a.semilogx(k, dUs, color=C_NEUT ||| with surface tension: short waves are held too, so the curve turns up again
b.semilogx(k, ch11.kh_growth_rate(k, 8.0, 0.0, 1.2, 1000.0), "--" ||| dashed: growth rate at 8 m/s without surface tension - rises forever
b.semilogx(k, ch11.kh_growth_rate(k, 8.0, 0.0, 1.2, 1000.0, surface_tension=0.074) ||| with surface tension: only a band of wavelengths grows
print(f"band at 8 m/s: k from ||| the two edges of the growing band, as wavenumbers and as wavelengths
print("thermocline (densities 1025 ||| the same criterion for a weak density step in the ocean
print("growth rate k c_i: deep" ||| a shallow lower layer slows the growth of this wave
print(f"lambda_c = {100*lam_c:.2f} cm ||| the longest ceiling wave that surface tension can hold
print(f"at t = {tt[i10]:.2f} s: linear ||| where the real sheet starts to lag the linear law, and how far apart they end
print(f"heated from above: Ra ||| the same layer heated from the wrong side: negative Rayleigh number
print(f"air gap of 1 cm: Ra ||| how many kelvin an air gap needs to reach the threshold
print(f"water layer of 5 mm: Ra ||| ... and how few a thin water layer needs
print(f"at equal depth (1 cm) ||| like for like: water convects about 140 times more readily than air
print("three largest growth rates:" ||| one positive growth rate (it grows), all others negative (they decay)
print(f"  N = {N}: sigma_1 ||| the leading growth rate does not change with the number of points
print(f"e-folding time = ||| the growth rate turned into seconds; and a layer below the threshold for contrast
A[r, :] = 0.0; A[r, r] = 1.0; B[r, :] = 0.0 ||| this row now says 'temperature amplitude = 0 at the wall' for every sigma
print(f"from scratch (N = 24) ||| our own matrices give the library's growth rate
a.plot(lam24.real, np.full(len(lam24), 24) ||| every eigenvalue of the coarse discretisation
a.plot(lam36.real, np.full(len(lam36), 36) ||| ... and of the finer one: the physical ones sit in the same place, the spurious ones do not
b.plot(mode["W"]/np.max(np.abs(mode["W"])) ||| vertical-velocity shape of the first cell: flat at the walls (no slip)
b.plot(mode["T"]/np.max(np.abs(mode["T"])) ||| temperature shape of the first cell: one arch, zero at the walls
big24, big36 = ||| the huge eigenvalues that row replacement and round-off create
print(f"eigenvalues larger than 1e6 ||| how many spurious eigenvalues each resolution has
print("leading three, N = 24:" ||| the physical eigenvalues agree between the two resolutions
print(f"K = {K}: Chebyshev eigenproblem ||| two independent routes to the marginal Rayleigh number at this cell width
print(f"critical point: Ra_c = ||| the bottom of the valley, next to the published value
print(f"relative difference from the published ||| agreement with the benchmark, and the width of a pair of rolls
print(f"my determinant root ||| the Rayleigh number at which our hand-typed determinant vanishes
tab = {"K": t_["K"] ||| the four neutral curves under short names
cr = {"rigid": crit ||| the critical point of each pair of walls (first two here, the others on the next line)
ax.semilogx(tab["rigid"], tab["K"] ||| marginal curve for rigid plates: to its right, cells of this width grow
ax.semilogx(tab["rigid_free"], tab["K"] ||| one free surface: convection starts earlier
ax.semilogx(tab["free"], tab["K"] ||| two free surfaces: earliest of all (C05)
ax.semilogx(tab["odd"], tab["K"] ||| two rows of cells: needs ten times the heating
ax.plot(cr[name]["Ra_c"], cr[name]["K_c"], "D" ||| diamond: the critical point - the first cell width to convect
ax.plot([2500, 2500], [K1, K2] ||| rose bar: the cell widths that grow at Ra = 2500 between rigid plates
print(f"at Ra = 2500 (rigid walls) ||| the edges of that bar
print({n: (round(v["Ra_c"], 2) ||| critical Rayleigh number and wavenumber for each pair of walls
out = {"rigid–rigid": (tab["K"], tab["rigid"]) ||| the three marginal curves (the same for every slider position)
figF2 = slider_figure(f2, "Ra" ||| slider: the Rayleigh number, on a logarithmic scale
ylabel="Rayleigh number Ra", title="Which cell widths ||| axis label and the question the figure answers
xrange=(0.5, 8), yrange=(400, 3e4)) ||| fixed axes, so only the bands move
bc = {"rigid–rigid": ("rigid", "rigid") ||| the drop-down text turned into the solver's (bottom, top) wall types
except ValueError: ||| the solver could not confirm the leading mode at this point
ax.plot(m["W"]/np.max(np.abs(m["W"])), m["z"] ||| the cell's vertical-velocity shape: rose if it grows, teal if it decays
live(benard_live, ||| sliders for cell width, heating and Prandtl number, a menu for the walls
f"{ch11.benard_growth_rate(Kc, 1.3*Rac, 0.7) ||| the same two growth rates for air: other values, same signs
print(f"rigid-free: Ra_c ||| a free top lowers the threshold; published value alongside
print(f"odd mode: Ra_c ||| two rows of cells need ten times more heating than one
print("dRa/dK^2 (step 9):" ||| the slope of the neutral curve, where it vanishes, and the value there
print(f"minimum by minimize_scalar ||| numerical minimum = 27 pi^4/4, and the growth rate typed by hand
a.plot(Kp, ch11.benard_free_free_Ra(Kp), color=C_NEUT ||| free-free marginal curve, one half-wave across the layer: below it nothing grows
a.plot(Kp, ch11.benard_free_free_Ra(Kp, n=2) ||| dashed: two half-waves across the layer - sixteen times harder to excite
a.plot(tab["K"], tab["rigid"], ":" ||| dotted ghost: rigid plates, for comparison
a.plot(Kff, 27*np.pi**4/4, "D" ||| diamond: the free-free critical point 27 pi^4/4 at K = pi/sqrt(2)
print("at Ra = 2000 the growing wavenumbers are" ||| the band of growing cell widths at Ra = 2000 ...
round(max(ch11.benard_free_free_sigma(K, 2000.0, 1.0)[0] for K in Kp), 2)) ||| ... and the fastest growth inside it
print("densities [kg/m^3]:" ||| warm salty water on top is lighter: the column is statically stable
print("tau = kappa_s/kappa =" ||| salt diffuses about a hundred times more slowly than heat
print("onset type of the growing root:" ||| fingers grow in place; the diffusive regime oscillates while it grows
"stationary" if abs(s_d[0].imag) ||| (second half of the same statement: the diffusive case)
print(f"  {name:15s} ||| one named result per line
print(f"thinnest finger-unstable layer ||| any layer thicker than this already forms fingers
cases = {"hot salty over cold fresh" ||| five layers: the two double-diffusive arrangements ...
"warm on top, no salt" ||| ... and three with temperature only (stable, weakly top-heavy, strongly top-heavy)
print(f"{label:28s} -> ||| the regime of each layer, its fastest growth rate, and whether it is lighter on top
print(f"left side of (11.46) by hand ||| the criterion typed by hand, and the check that no salt gives Benard convection back
code_of = {"stable": 0 ||| one integer per regime, so the map can be coloured
cmap = ListedColormap( ||| teal stable, amber fingers, rose diffusive, grey overturning
ax.plot(ax_T*1e6, ax_T*1e6, color=C_BUOY ||| blue: salt and heat cancel in the density - below this line the column is lighter on top
ax.plot(ax_T*1e6, 1.5e-9*(ax_T/1.4e-7 + C_map)*1e6 ||| amber: the finger criterion - above this line salt fingers grow
ax.plot(al*0.01*1e6, be*0.002*1e6, "o" ||| our thermocline: inside the wedge - lighter on top, yet fingering
ax.annotate("diffusive (oscillatory) ||| points at the thin band of growing oscillations
arrowprops=dict(arrowstyle="-" ||| a plain line as the pointer
wedge_y = np.concatenate( ||| ... and the heights of that outline
return {"density line (lighter on top below it)" ||| the fixed density line and the finger line for this ratio of diffusivities
"wedge: lighter on top, finger-unstable" ||| the region between them: stable by density, unstable by diffusion
figF3 = slider_figure(f3, ||| slider: salt diffusivity over heat diffusivity
ylabel="β dS/dz [1e-6 1/m]", title="A 5 cm layer ||| axis label and the message
xrange=(-1, 4), yrange=(-1, 4)) ||| fixed axes, so only the wedge moves
print("inner cylinder spinning: Rayleigh-unstable?" ||| inner cylinder driving: unstable to ring swaps; outer faster: stable
print("Rayleigh line for this gap ||| the rotation ratio at which the squared circulation stops falling outward
print("continuity as printed has consistent units?" ||| slip #12 as a units test
print(f"mu = {mu}: sigma = ||| a real, positive growth rate: vortices grow in place
print(Ta, "| narrow-gap shortcut:" ||| the Taylor number of our cylinders, and its thin-gap estimate
print(f"mu = {mu:5.2f}: Ta_c = ||| exact threshold and preferred wavenumber, against the approximate formula
print(f"onset speed for our cylinders ||| the inner-cylinder speed at which vortices first appear
print(f"Ta by hand = ||| the three formulas typed out agree with the library
a.plot(tt["mu"], tt["Ta_c"], "o" ||| dots: the exact narrow-gap threshold for each rotation ratio
a.plot(tt["mu"], tt["approx"], "--" ||| dashed: the book's approximate formula - good for co-rotation, poor for counter-rotation
a.plot(1.0, BENCH["benard_rigid_rigid"]["Ra_c"], "D" ||| diamond: at equal rotation rates the threshold is Benard's 1707.76
a2.plot(tt["mu"], 100*tt["rel_error"] ||| how far the approximate formula is off, in per cent
b.plot(sb["x_outer"], sb["y_inner"], color=C_NEUT ||| the viscous threshold: above this curve Taylor vortices appear
b.plot(sb["rayleigh_x"], sb["rayleigh_y"], "--" ||| dashed: Rayleigh's inviscid border - right of it nothing can drive vortices
print(f"table: Ta_c(0.25) ||| one stored threshold recomputed now
print(f"axial wavelength 2 pi/k_c ||| each vortex is about as tall as the gap is wide
print(f"Galerkin, 4 sines: ||| two different discretisations give the same threshold
print("with the printed operator:" ||| the printed (squared) operator has no threshold at all
print(f"N^2 = {N2_thermo:.2e} ||| stratification of a thermocline: a displaced parcel bobs every couple of minutes
print("unstable eigenvalue c =" ||| a purely imaginary wave speed: the billow grows without travelling
print("N^2 = 0: Taylor-Goldstein" ||| without stratification the two solvers must agree, and do
raw2 = ch11.taylor_goldstein_eigs(0.4, prof["U"] ||| the same raw spectrum with more points, to see what moves
print("J = 0.3: converged unstable eigenvalues:" ||| stronger stratification: no growing mode is left
SA = sp.diff((U - c)*sp.diff(phi, z), z) ||| the self-adjoint form (11.64): its derivative term ...
print(f"J = {J}: Ri_min = ||| the smallest Richardson number of the profile and the theorem's verdict
print("growth rate k c_i at k = 0.4 for J ||| computed growth: positive below 1/4, zero above
idn = ch11.richardson_identity_check( ||| both sides of D18's integral identity for this computed mode
print(f"(11.65): left ||| the identity holds to round-off
print(f"imaginary parts: left ||| ... and so does its imaginary part, the one that gives the theorem
print(f"tanh layer: Ri_min ||| the minimum sits at the centre of the layer; the thermocline example is below 1/4
a.plot(np.tanh(zz), zz, color=C_SHEAR ||| the velocity profile of the shear layer
a.plot(pj["N2"](zz), zz, ls, color=C_BUOY ||| the stratification, concentrated where the shear is
b.plot(kk, ch11.tg_tanh_neutral_J(kk) ||| exact neutral curve of this family: growth only below it
print(f"largest growth in the map ||| fastest growth is without stratification; one stored point recomputed
e = ch11.taylor_goldstein_eigs(k, pj["U"] ||| unstable eigenvalue(s) at this k, in a box and on nodes suited to its wavelength
print(len(cs), "unstable modes; all inside ||| every computed unstable wave speed obeys Howard's bound
print({n: (f"{v:.1e}" if isinstance(v, float) else v) for n, v in hw.items()}) ||| residuals of D19's two integral identities for a computed mode
print("all", len(cs), "modes satisfy ||| the inequality typed by hand agrees; the growth bound is safe but not tight
a.plot(arc_r, arc_i, color=C_NEUT ||| Howard's arc: no unstable wave speed can lie outside it
a.plot([-1, 1], [0, 0], color=C_SHEAR ||| the range of flow speeds, which is the arc's diameter
a.plot(cs[sel].real, cs[sel].imag, mk ||| computed unstable wave speeds: all on the imaginary axis, well inside
b.plot(ks_m[sel], ks_m[sel]*cs[sel].imag ||| their growth rates
b.plot(kk, kk*(1.0 - (-1.0))/2, "--" ||| dashed: the largest growth rate Howard's theorem allows
cu += [e[0]] if len(e) else [] ||| keep the growing eigenvalue if there is one
figF6 = slider_figure(f6, "b" ||| slider: the half-width of the channel
title="U = sin y between walls ||| the message
modes={"unstable c (k = 0.2 ||| eigenvalues drawn as markers
print({n: (round(v, 4) if isinstance(v, float) else v) for n, v in sq.items()}) ||| the equivalent two-dimensional wave: longer wavenumber, lower Reynolds number
print("leading c, 3-D problem:" ||| Squire's theorem in numbers: the oblique wave and its 2-D twin share one eigenvalue
print(f"  {i}. {s}") ||| one line per move of the symbolic derivation
print(f"from scratch: c = ||| our own matrices reproduce the growing Tollmien-Schlichting eigenvalue
ax.plot(rawA.real, rawA.imag, "o" ||| every eigenvalue at N = 100: the three branches of the Y
ax.plot(rawB.real, rawB.imag, "x" ||| ... and at N = 140: where the two differ the mode is not converged
ax.plot(c_os[0].real, c_os[0].imag, "o", ms=11 ||| circled: the single eigenvalue above the axis - the growing wave
a.plot(sp_ray.real, sp_ray.imag, "o" ||| inviscid eigenvalues: symmetric about the real axis
a.plot([0, 0], [sp_ray[0].imag, -sp_ray[0].imag] ||| circled: the growing mode and its decaying mirror image
b.plot(sp_os[Re].real, sp_os[Re].imag, mk ||| viscous eigenvalues: one grows, nothing mirrors it
print("inviscid pair: c =" ||| the growing mode and its conjugate twin
print("viscous leading c: Re = 100:" ||| with viscosity the growing mode approaches the inviscid one as Re rises
print("c =", np.round(mR["c"][0], 4) ||| D22's identity holds for the computed growing mode
print("plane Poiseuille, inviscid: ||| no inflection point, no inviscid instability
print("panel | profile | inflection point ||| header of the verdict table
print(f"  ({v['panel']}) {v['label']:34s} ||| one profile per row: does it pass Rayleigh's and Fjortoft's tests?
print(f"fastest tanh wave: k = ||| the most dangerous wavelength of a shear layer, against the published value
print(f"inflection at y = ||| the vorticity is largest at the inflection point: both criteria met
fig.supylabel("cross-stream coordinate ||| one vertical label for all six panels
print(f"b = {b}: 2b ||| the inflected sine profile grows only when the channel is wider than pi
print("critical layer of U = tanh y ||| the height at which the flow moves exactly at the wave speed ...
"| artanh(0.3) =" ||| ... which for a tanh profile is artanh(c)
ax.plot([np.pi, np.pi], [-half, half] ||| rose bar: the height of one eye
print("half-height of the eye:" ||| the eye grows like the square root of the wave amplitude
print(f"neutral at kh = {kh_n:.4f} ||| the piecewise-linear layer: where it turns neutral, and its growth at one wavelength
print("jet: sinuous c =" ||| a jet has two kinds of growing wave: flapping and pulsing
print(f"Re_c = {crit_P['Re_c']:.2f} ||| the Reynolds number at which a channel flow first becomes unstable
hi = np.interp(np.log(1e4), np.log(nc["Re"]), nc["k_upper"]) ||| ... and its upper edge
print(f"at Re = 10^4 the wavenumbers ||| the band of growing waves well above the critical point
print("plane Couette flow, leading c_i ||| negative at both Reynolds numbers: Couette flow damps this wave
print(f"c_i = 0 at k = 1.02056 ||| the root search lands on the stored critical Reynolds number
k_lo, k_hi = np.concatenate( ||| the lower and upper neutral wavenumbers, both starting at the tip
a.semilogx(Re_th, k_lo, color=C_NEUT ||| lower side of the thumb: longer waves decay
a.semilogx(Re_th, k_hi, color=C_NEUT ||| upper side of the thumb: shorter waves decay
a.plot(crit_P["Re_c"], crit_P["k_c"], "D" ||| diamond: the tip - the first wave to grow as Re rises
a.plot([1e4, 1e4], [lo, hi] ||| rose bar: the growing wavenumbers at Re = 10^4
b.semilogx(tl["Re"], tl["k_upper"], color=C_SHEAR ||| shear layer: every wave longer than this grows, at any Re
b.semilogx(bj["Re"], bj["k_upper"], color=C_BUOY ||| jet: upper neutral wavenumber, rising toward the inviscid value 2
b.semilogx(bl["Re"], bl["k_upper"], color=C_NEUT ||| boundary layer: upper side of its thumb
b.semilogx(bl["Re"], bl["k_lower"], color=C_NEUT ||| ... and lower side
b.semilogx(Re_th, k_hi, ":" ||| dotted: the channel-flow thumb of panel (a), far to the right
b.semilogx(Re_th, k_lo, ":" ||| ... its lower side
print(f"first Re in each table ||| how early each flow becomes unstable
print(f"Bickley jet: critical point ||| the jet's first unstable wave
ax.plot(Rej[unres], np.full(unres.sum(), K0), "|" ||| ticks: the instability reaches the smallest wavenumber computed
print(f"one band up to Re = ||| where the stable gap first appears in the table, and its two edges
print(f"upper neutral wavenumber at Re = ||| viscosity matters less and less: the upper edge approaches the inviscid limit
print("Rayleigh residual for phi = sech y ||| zero: sech y is an exact neutral mode of the tanh layer
print(f"  Re = {Re:5.0f}: unstable for k < ||| viscosity trims the unstable band only slightly
print(np.round(ci, 4)) ||| leading growth indicator for each (k, Re): rows are k, columns Re
print("largest c_i on the grid:" ||| negative everywhere: no growing wave in plane Couette flow
print(tab11[["flow" ||| our recomputed table of critical Reynolds numbers, with sources
print("  -", t) ||| one term of the energy equation per line
budgets[Re] = (m, b) ||| keep the mode and its budget for the figure and the animation
print(f"Re = {Re:.0f}: c = ||| the wave speed: positive imaginary part at 10^4, negative at 5000
print(f"   E = {b['E']:.5f} ||| the wave's energy, its rate of change, what the shear supplies and what viscosity removes
print(f"   P - dissipation - dE/dt ||| the budget closes; a ratio above 1 means the wave grows
print(f"production by trapezoid ||| plain trapezoid integration gives the same budget
m, b = budgets[1e4] ||| the growing wave and its budget
axs[0].plot(-b["uv"], b["y"], color=C_SHEAR ||| the Reynolds stress: nonzero only near the walls
axs[0].plot(b["production_density"], b["y"] ||| where the wave takes energy from the mean shear
bb = budgets[Re][1] ||| the budget at this Reynolds number
vals = [bb["production"], bb["dissipation"], bb["dEdt"]] ||| supply, loss and their difference
print(f"critical layers at y = ||| where the flow moves at the wave speed, and how far the phase is from a quarter cycle
a2.plot(-b["uv"], b["y"], color=C_SHEAR) ||| the Reynolds stress of this wave across the channel
print(f"one period = ||| in one period the wave moves one wavelength and grows by ten per cent
print({n: round(v, 4) for n, v in co.items()} ||| the two coefficients from continuity, and the jump of the printed form
ax.plot(fB, eta, color=C_BASE ||| grey: the exact Blasius profile
ax.plot(ch11.tollmien_profile(eta), eta ||| Tollmien's piecewise stand-in: straight, then a parabola, then uniform
ax.plot(bF["Re"], 1e6*bF["F_lower"] ||| lower edge of the tongue of growing frequencies
ax.plot(bF["Re"], 1e6*bF["F_upper"] ||| upper edge of the tongue
print(f"a wave with F = 100e-6 grows ||| the stretch of plate over which this ribbon frequency is amplified
print(f"ours: Re_c = {bc_['Re_c']:.2f} ||| the boundary layer's critical point, next to the published one
a.plot(sol28["t"], sol28["X"] ||| strong heating: the roll keeps reversing irregularly
a.plot(sol10["t"], sol10["X"] ||| gentle heating: the roll settles to a steady speed
f"(eigenvalue of the steady roll ||| the ringing frequency predicted by linearising about the steady roll
drift = max( ||| how much the energy changes along each trajectory (it should not)
print(f"largest change of the energy ||| conserved to round-off: no friction, no attractor
ax.plot(h["X"], h["Y"], color=C_DECAY ||| a trajectory: into the fixed point (teal) or onto the limit cycle (orange)
axs[2].plot(mus, ch11.limit_cycle_amplitude(mus) ||| size of the attractor: zero, then growing like the square root
axs[2].plot(mus[mus > 0], 0*mus[mus > 0], "--" ||| dashed: the fixed point is still there, but now repels
a2 = sp.pi**2 + k**2 ||| the combination pi^2 + k^2 that every Laplacian of the roll produces
print("eigenvalues at the origin:" ||| one positive eigenvalue: the state of rest is unstable
print("Hopf point r_H =" ||| the heating beyond which the steady rolls are unstable, and b for the critical roll
"| volume contraction rate:" ||| volumes of states shrink at this rate
print(f"separation grows like ||| the fitted exponent of the error growth
print("time at which the two runs differ by 1:" ||| the forecast horizon for a starting error of 1e-8
X, Y, Z = s ||| roll speed, temperature contrast, profile distortion
print(f"largest difference up to t = 5 ||| two correct solvers agree at first and disagree completely later
fig = plt.figure(figsize=(10.2, 3.6)) ||| an empty figure for three panels (the first one three-dimensional)
ax3.plot(sol28["X"], sol28["Y"], sol28["Z"], color=C_BASE ||| the whole orbit: two lobes, never closing
ax3.plot([p_[0]], [p_[1]], [p_[2]], "o" ||| the three steady states, all unstable at r = 28
ax3.view_init(18, -60) ||| a viewpoint from which both lobes are visible
a = fig.add_subplot(1, 3, 2) ||| panel (b)
a.plot(sep["t"], sep["a"][0], color=C_SHEAR ||| the roll speed of the first run
a.plot(sep["t"], sep["b"][0], color=C_BUOY ||| the second run, started a hundred-millionth away
b = fig.add_subplot(1, 3, 3) ||| panel (c)
b.semilogy(sep["t"][1:], sep["sep"][1:] ||| distance between the two runs: a straight line here means exponential growth
tw = np.array([w0, w1]) ||| the two ends of the fitted window
b.semilogy(tw, sep["sep"][i0]*np.exp(sep["slope"]*(tw - w0)) ||| the fitted exponential: its slope is the Lyapunov exponent
(d1,) = axs[0].plot([], [], "o", color=C_SHEAR) ||| the current state of run 1
(d2,) = axs[0].plot([], [], "o", color=C_BUOY) ||| the current state of run 2
axs[2].semilogy(sep["t"][1::4], sep["sep"][1::4] ||| how far apart the two runs are
(mark,) = axs[2].semilogy( ||| a dot that follows the present time along that curve
j = idx[i] ||| the time index shown in this frame
live(lorenz_live, ||| sliders for the heating, the Prandtl number, b and the size of the first error
x = A*x*(1 - x) ||| one application of the logistic map
print("period doublings at A =" ||| the parameter values at which the cycle doubles its period
print("estimates of Feigenbaum's ratio:" ||| the ratios of successive gaps approach a universal number
b.plot(bd["A"], bd["x"], "," ||| every value the map keeps visiting at this A: one branch, two, four, ... chaos
'''


def apply_hand_comments() -> tuple[int, list[str]]:
    table = [tuple(x.strip() for x in ln.split("|||")) for ln in HAND_COMMENTS.strip("\n").splitlines() if "|||" in ln]
    used, changed = set(), 0
    for c in nb.cells:
        if c.cell_type != "code" or "setup" in c.metadata.get("tags", []):
            continue
        new = []
        for ln in c.source.split("\n"):
            body = ln.strip()
            if body and "#" not in ln:
                for j, (pre, com) in enumerate(table):
                    if body.startswith(pre):
                        ln = ln.rstrip() + "   # " + com
                        used.add(j)
                        changed += 1
                        break
            new.append(ln)
        c.source = "\n".join(new)
    return changed, [table[j][0] for j in range(len(table)) if j not in used]


_n_hand, _unused = apply_hand_comments()
print(f"hand-written comments attached: {_n_hand}")
if _unused:
    raise SystemExit("HAND_COMMENTS entries that match no code line:\n  " + "\n  ".join(_unused))


# ---------------------------------------------------------------------------------------------------------------------
# Specific comments first (lesson review): a print line says WHAT it prints, a plot line WHICH curve it draws; only the
# lines this pass cannot describe fall through to the generic boilerplate pass above.
# ---------------------------------------------------------------------------------------------------------------------
def _top_args(s: str) -> list[str]:
    """Top-level comma-separated arguments of a call's argument text."""
    out, depth, cur, quote = [], 0, "", ""
    for ch in s:
        if quote:
            cur += ch
            if ch == quote:
                quote = ""
            continue
        if ch in "\"'":
            quote = ch
        elif ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def _plain(text: str, n: int = 72) -> str:
    text = re.sub(r"\{[^{}]*\}", "…", text).replace("$", "").replace("\\", "")
    text = re.sub(r"\s+", " ", text).strip(" :,;-|")
    return text[:n].rstrip() + ("…" if len(text) > n else "")


def annotate_specific() -> int:
    changed = 0
    for c in nb.cells:
        if c.cell_type != "code" or "setup" in c.metadata.get("tags", []):
            continue
        new, in_str = [], False
        for ln in c.source.split("\n"):
            if ln.count('"""') % 2 == 1:
                in_str = not in_str
                new.append(ln)
                continue
            body = ln.rstrip()
            if in_str or "#" in ln or not body or body.count("(") != body.count(")") or not body.endswith(")"):
                new.append(ln)
                continue
            msg = ""
            m = re.match(r"\s*print\((.*)\)$", body)
            if m:
                msg = "@keep"                                    # a print line: no machine-made comment (the list under the cell explains the output)
            if not msg and re.search(r"\.(plot|semilogx|semilogy|loglog)\(", body):
                msg = "@keep"                                    # a plot line without a hand comment: the notes under the figure say what the curves are
            m = re.search(r"""set_title\(f?r?["']([^"']+)["']""", body)
            if not msg and m:
                msg = "panel title — the message of this panel"
            if not msg:
                for pat, text in ((r"\.set\((?:title|xlabel|ylabel|xlim|ylim)", "axis labels (with units) and fixed ranges for this panel"),
                                  (r"^\s*fig\w*\.show\(\)$", "display the interactive slider figure"),
                                  (r"\.colorbar\(", "the colour scale and what it measures"),
                                  (r"\.set_text\(", "update the title text for this frame"),
                                  (r"\.set_(?:array|data|ydata)\(", "move the drawn data to this frame's values"),
                                  (r"^\s*show_animation\(", "render the frames and show the animation")):
                    if re.search(pat, body):
                        msg = text
                        break
            if msg == "@keep":
                ln = body + "   #@keep"
            elif msg:
                ln = body + "   # " + msg
                changed += 1
            new.append(ln)
        c.source = "\n".join(new)
    return changed


print(f"code lines given a specific comment: {annotate_specific()}")


# ---------------------------------------------------------------------------------------------------------------------
# final pass: labels tidied, every equation named by number written out, self-checks, save (nbkit's coverage checks)
# ---------------------------------------------------------------------------------------------------------------------
def relabel() -> int:
    """Label hygiene: the design's working codes for explainers (E1…E9), animations (A1…A5) and slider figures (F1…F8)
    are replaced by reader-facing words in every markdown cell (Part F texts use the codes)."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        s = re.sub(r"\bE([1-9])(?:'s)?\b(?![_^])", lambda m: f"explainer {m.group(1)}" + ("'s" if m.group(0).endswith("'s") else ""), c.source)
        s = re.sub(r"\(A([1-5])\)", r"(animation A\1)", s)
        s = re.sub(r"(?<![A-Za-z(])(?<!Animation )(?<!animation )\bA([1-5])\b(?![_^,\d])(?= \(| below| above|\)|;|\.)", r"animation A\1", s)
        s = re.sub(r"\bB1's chips\b", "another one-hump map", s)
        if s != c.source:
            c.source, changed = s, changed + 1
    return changed


nb.cells[0].source = nb.cells[0].source.replace("equations are cited by their numbers so you can follow along",
                                                "equations are shown in full with their book numbers so you can follow along")
n_annot = annotate_boilerplate()
for _c in nb.cells:                                         # print lines keep no generic comment either
    if _c.cell_type == "code":
        _c.source = _c.source.replace("   #@keep", "")
print(f"plotting boilerplate lines commented: {n_annot}")
print(f"labels tidied in {relabel()} cells")
n_changed = finalize_equations()
n_tex = tidy_raw_tex()
print(f"plain-text exponents turned into maths in {n_tex} cells; equations written out in {n_changed} cells")
bad = self_check_ctrl() + self_check_numbers() + self_check_prose() + self_check_near()
if "--partial" in sys.argv:                                  # development aid: write what exists so far, unchecked
    import nbformat as _nbf
    for _m in bad:
        print("CHECK", _m)
    _nb2 = _nbf.v4.new_notebook(cells=list(nb.cells))
    _nb2.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    _out = ROOT / "notebooks" / "ch11_instability.ipynb"
    _out.parent.mkdir(parents=True, exist_ok=True)
    for _i, _c in enumerate(_nb2.cells):
        _c["id"] = f"ch11-{_i:03d}"
    _nbf.write(_nb2, _out)
    print("wrote (partial, unchecked)", _out, len(nb.cells), "cells")
    sys.exit(0)
if bad:
    raise SystemExit("markdown cells cite equation numbers without maths (or TeX outside maths):\n  " + "\n  ".join(bad))
miss = self_check_ledger()
if miss:
    raise SystemExit("Part E rows explained by a primer but no primer/reminder in the notebook:\n  " + "\n  ".join(miss))
out = nb.save()
print(f"wrote {out.relative_to(ROOT)} ({len(nb.cells)} cells: "
      f"{sum(c.cell_type == 'markdown' for c in nb.cells)} markdown, {sum(c.cell_type == 'code' for c in nb.cells)} code)")
