"""Build the Chapter 12 teaching notebook: ``notebooks/ch12_turbulence.ipynb``.

Source of truth: ``analysis/ch12_design.md`` Part A (storyboard, one nbkit call per row, with the A.15a table that attaches
further contract functions to their notes), Part C (the function contract — ``fluidpy.ch12_turbulence`` imported as
``ch12`` re-exports ``core.turbstats`` = ``TS`` and ``core.wall_turbulence`` = ``WT``), Part D (runtime budget and FAST
sizes), Part E (prerequisite ledger → primers P280–P306 and one-line reminders of earlier primers), Part F (the 28
derivations D01–D28, one move per step) and ``analysis/ch12_curation.md`` (IDs, depths, section coverage), with the rulings
of ``reports/ch12_verification.md`` written into the text:

* the log-law fit (C11) uses the public DNS subset ``reference/ch12/lee_moser_2015_channel_mean.csv`` (Lee & Moser 2015,
  J. Fluid Mech. 774, 395) with ``window=(350.0, 0.15)``; a small table shows how the fitted pair depends on the window;
* ``WT.total_stress`` returns a dict; ``WT.LOG_LAW_CONSTANTS`` entries are never splatted (``kappa=``, ``B=`` passed explicitly);
* stability is read from the ``verdict`` keys of ``ch12.surface_layer_regime``; the dispersion-regime boundaries are inclusive;
* ``jet_tke_budget`` and the free-shear amplitude constants are labelled illustrative; ``spalding_uplus``,
  ``stress_partition`` and ``channel_mixing_length`` are labelled models (approximate against DNS);
* the lapse-rate sign: Kundu's Γ ≡ dT/dz and the meteorological Γ_met ≡ −dT/dz are both shown wherever a lapse rate appears;
* the printed slips listed by ``ch12.book_slips()`` (18) are taught in corrected form with a short box (``slip #k``).

**Derivations are read from Part F at build time** (``part_f()`` below): goal, start, plan, tools, assumptions, every step's
*did / tex / why / plain*, result, check, meaning and traps are copied word for word; the ★★★ D07, D08, D10, D14 carry sympy
check cells.

Book numbers never printed (rule 9): our own inputs and labelled illustrative constants everywhere (design convention 10);
public benchmark values are read from ``reference/ch12/`` and cited.

Run:  .venv/Scripts/python.exe notebooks/build_ch12.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch12.py --dump     (prints the parsed Part F derivations only)
      .venv/Scripts/python.exe notebooks/build_ch12.py --partial  (writes what exists, unchecked — development aid)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

nb = ChapterNotebook("ch12")
DESIGN = (ROOT / "analysis" / "ch12_design.md").read_text(encoding="utf-8")
_CTRL = sorted({hex(ord(_c)) for _c in DESIGN if ord(_c) < 32 and _c != "\n"})
if _CTRL:                                # a form feed, tab or backspace in the design = a LaTeX backslash eaten by a non-raw string
    _k = next(j for j, _c in enumerate(DESIGN) if ord(_c) < 32 and _c != "\n")
    raise SystemExit(f"analysis/ch12_design.md contains control characters {_CTRL} (line {DESIGN.count(chr(10), 0, _k) + 1}): "
                     f"…{DESIGN[max(0, _k - 40):_k + 20]!r} — repair the eaten backslash (\\f, \\t, \\b, \\r, \\a, \\v) before building")
REMIND: dict[str, str] = {}          # filled below: one-sentence reminders of earlier primers, written for this chapter
NOTES_IN: dict[str, str] = {}        # filled below: notes of the curation taught inside a derivation

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.  The LaTeX
# is the analyst's page-image transcription as carried by the design (42 of the chapter's 80 pages re-read there; the
# pages re-rendered for this builder are listed in the builder's report), in the CORRECTED form wherever the book prints
# a slip (#1 … #15).
# ---------------------------------------------------------------------------------------------------------------------

EQ: dict[str, str] = {
    "12.1": r"\langle u^m(\mathbf x,t)\rangle=\lim_{N\to\infty}\frac1N\sum_{n=1}^N\big(u(\mathbf x,t{:}n)\big)^m",
    "12.2": r"\overline{u^m(\mathbf x)}=\frac1{\Delta t}\int_{t-\Delta t/2}^{t+\Delta t/2}u^m(\mathbf x,t)\,dt",
    "12.3": r"\overline{u^m(t)}=\frac1V\int_Vu^m(\mathbf x,t)\,dV",
    "12.4": r"\overline{u^m+v^m}=\overline{u^m}+\overline{v^m}",
    "12.5": r"\overline{Au^m}=A\overline{u^m}",
    "12.6": r"\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t",
    "12.7": r"\overline{\int u\,dt}=\int\bar u\,dt",
    "12.8": r"\overline{\partial u^m/\partial x_j}=\partial\overline{u^m}/\partial x_j",
    "12.9": r"\overline{\int u^m\,d\mathbf x}=\int\overline{u^m}\,d\mathbf x",
    "12.10": r"\overline{u(\mathbf x,t)}\equiv\frac1N\sum_{n=1}^Nu(\mathbf x,t{:}n)",
    "12.11": r"\overline{(u-\langle u\rangle)^m}\equiv\frac1N\sum_n\big(u(\mathbf x,t{:}n)-\langle u(\mathbf x,t)\rangle\big)^m",
    "12.12": r"R_{ij}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_i(\mathbf x_1,t_1)u_j(\mathbf x_2,t_2)}",
    "12.13": r"R_{11}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_1(\mathbf x_1,t_1)u_1(\mathbf x_2,t_2)}",
    "12.14": r"r_{12}\equiv R_{12}/\sqrt{R_{11}R_{22}}=\overline{u_1u_2}\big/\big(\sqrt{\overline{u_1^2}}\sqrt{\overline{u_2^2}}\big)",
    "12.15": r"r_{11}\equiv R_{11}(\mathbf x_1,t_1,\mathbf x_2,t_2)/\sqrt{R_{11}(1,1)R_{11}(2,2)}",
    "12.16": r"\lvert\overline{uv}\rvert\le\sqrt{\overline{u^2}}\sqrt{\overline{v^2}}",
    "12.17": r"R_{11}(\\tau)=\\overline{u_1(t)u_1(t+\\tau)}=R_{11}(-\\tau)",
    "12.18": r"\Lambda_t\equiv\int_0^\infty r_{11}(\tau)\,d\tau=\frac1{R_{11}(0)}\int_0^\infty R_{11}(\tau)\,d\tau",
    "12.19": r"\lambda_t^2\equiv-2\big/\big[d^2r_{11}/d\tau^2\big]_{\tau=0}",
    "12.20": r"S_e(\omega)\equiv\frac1{2\pi}\displaystyle\int_{-\infty}^{+\infty}R_{11}(\tau)\,e^{-i\omega\tau}\,d\tau",
    "12.21": r"R_{11}(\tau)\equiv\int_{-\infty}^{+\infty}S_e(\omega)e^{+i\omega\tau}d\omega",
    "12.22": r"\\overline{u_1^2}=\\int_{-\\infty}^{+\\infty}S_e(\\omega)\\,d\\omega",
    "12.23": r"R_{ij}(\mathbf r)\equiv\overline{u_i(\mathbf x)u_j(\mathbf x+\mathbf r)}",
    "12.24": r"\tilde u_i=U_i+u_i",
    "12.25": r"\overline{\tilde T}=\bar T",
    "12.26": r"\overline{u_i}=0",
    "12.27": r"\partial U_i/\partial x_i=0",
    "12.28": r"\partial u_i/\partial x_i=0",
    "12.29": r"\frac{\partial(U_i+u_i)}{\partial t}+\frac{\partial}{\partial x_j}\big((U_j+u_j)(U_i+u_i)\big)=-\frac1{\rho_0}\frac{\partial(P+p)}{\partial x_i}-g[1-\alpha(\bar T+T'-T_0)]\delta_{i3}+\nu\frac{\partial^2(U_i+u_i)}{\partial x_j^2}",
    "12.30": r"\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\frac1{\rho_0}\frac{\partial\bar\tau_{ij}}{\partial x_j}",
    "12.31": r"\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}+\frac{\partial}{\partial x_j}(\overline{u_jT'})=\kappa\frac{\partial^2\bar T}{\partial x_j^2}",
    "12.32": r"Q_j=-k\frac{\partial\bar T}{\partial x_j}+\rho_0C_p\overline{u_jT'}",
    "12.33": r"\frac{\partial}{\partial t}(\rho_m\tilde Y)+\frac{\partial}{\partial x_j}(\rho_m\tilde Y\tilde u_j)=\frac{\partial}{\partial x_j}\big(\rho_m\kappa_m\frac{\partial\tilde Y}{\partial x_j}\big)",
    "12.34": r"\frac{\partial\bar Y}{\partial t}+U_j\frac{\partial\bar Y}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\kappa_m\frac{\partial\bar Y}{\partial x_j}-\overline{u_jY'}\big)",
    "12.35": r"\frac{\partial\overline{u_iu_j}}{\partial t}+U_k\frac{\partial\overline{u_iu_j}}{\partial x_k}+\frac{\partial\overline{u_iu_ju_k}}{\partial x_k}=-\overline{u_iu_k}\frac{\partial U_j}{\partial x_k}-\overline{u_ju_k}\frac{\partial U_i}{\partial x_k}-\frac1\rho\Big(\overline{u_i\frac{\partial p}{\partial x_j}}+\overline{u_j\frac{\partial p}{\partial x_i}}\Big)-2\nu\overline{\frac{\partial u_i}{\partial x_k}\frac{\partial u_j}{\partial x_k}}+\nu\frac{\partial^2\overline{u_iu_j}}{\partial x_k^2}+g\alpha\big(\overline{u_jT'}\delta_{i3}+\overline{u_iT'}\delta_{j3}\big)",
    "12.36": r"\overline{(\partial u_1/\partial x_1)^n}=\overline{(\partial u_2/\partial x_2)^n}=\overline{(\partial u_3/\partial x_3)^n}",
    "12.37": r"\overline{(\partial u_1/\partial x_2)^n}=\overline{(\partial u_1/\partial x_3)^n}=\dots=\overline{(\partial u_3/\partial x_2)^n}",
    "12.38": r"f(r)\equiv\overline{u_\parallel(\mathbf x+\mathbf r)u_\parallel(\mathbf x)}/\overline{u_\parallel^2}",
    "12.39": r"\lambda_g^2\equiv-2/[d^2g/dr^2]_{r=0}",
    "12.40": r"R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}",
    "12.41": r"R_{ij}=\overline{u^2}\{f\delta_{ij}+\frac r2f'(\delta_{ij}-r_ir_j/r^2)\}",
    "12.42": r"\bar\varepsilon=\frac\nu2\overline{(\partial u_i/\partial x_j+\partial u_j/\partial x_i)^2}",
    "12.43": r"\bar\varepsilon=15\nu\overline{u^2}/\lambda_g^2",
    "12.44": r"R_\lambda\equiv\lambda_{(g\text{ or }f)}\sqrt{\overline{u^2}}/\nu",
    "12.45": r"S_{11}(k_1)=\frac1{2\pi}\int_{-\infty}^{+\infty}R_{11}(r_1)e^{-ik_1r_1}dr_1=\frac{\overline{u_1^2}}{2\pi}\int_{-\infty}^{+\infty}f(r_1)e^{-ik_1r_1}dr_1",
    "12.46": r"\frac{\partial\bar E}{\partial t}+U_j\frac{\partial\bar E}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac{U_jP}{\rho_0}+2\nu U_i\bar S_{ij}-\overline{u_iu_j}U_i\big)-2\nu\bar S_{ij}\bar S_{ij}+\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}-\frac g{\rho_0}\bar\rho U_3",
    "12.47": r"\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}",
    "12.48": r"\dot W\sim\overline{u_iu_j}(\partial U_i/\partial x_j)\sim(\Delta U)^2[\Delta U/L]=(\Delta U)^3/L",
    "12.49": r"\bar\varepsilon\sim(\Delta U)^3/L",
    "12.50": r"u_K=(\nu\bar\varepsilon)^{1/4}",
    "12.51": r"\eta/L\sim\mathrm{Re}_L^{-3/4}",
    "12.52": r"\frac{\lambda_T}L\propto\mathrm{Re}_L^{-1/2}",
    "12.53": r"\frac{S_{11}(k_1)}{u_K^2\eta}=\Phi(k_1\eta)",
    "12.55": r"\int_{-\infty}^{+\infty}S_{11}(k_1)\,dk_1=\overline{u_1^2}=\int_0^{+\infty}2S_{11}(k_1)\,dk_1",
    "12.56": r"U=U_{CL}(x)F(y/\delta(x))",
    "12.57": r"-\overline{uv}=\Psi(x)G(y/\delta(x))",
    "12.58": r"\partial U/\partial x+\partial V/\partial y=0",
    "12.59": r"U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}=-\frac1\rho\frac{\partial P}{\partial x}+\nu\Big(\frac{\partial^2U}{\partial x^2}+\frac{\partial^2U}{\partial y^2}\Big)-\frac{\partial\overline{u^2}}{\partial x}-\frac{\partial\overline{uv}}{\partial y}",
    "12.60": r"U\frac{\partial V}{\partial x}+V\frac{\partial V}{\partial y}=-\frac1\rho\frac{\partial P}{\partial y}+\nu\Big(\frac{\partial^2V}{\partial x^2}+\frac{\partial^2V}{\partial y^2}\Big)-\frac{\partial\overline{uv}}{\partial x}-\frac{\partial\overline{v^2}}{\partial y}",
    "12.61": r"U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y}",
    "12.62": r"J_s=\rho\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}",
    "12.63": r"\big\{\frac{\delta U'_{CL}}{U_{CL}}\big\}F^2-\big\{\frac{\delta U'_{CL}}{U_{CL}}+\delta'\big\}F'\int_0^\xi F\,d\xi=\big\{\frac\Psi{U_{CL}^2}\big\}G'",
    "12.64": r"\frac\Psi{U_{CL}^2}=C_3",
    "12.65": r"J_s=\rho\int_{-\infty}^{+\infty}U^2dy=\rho U_{CL}^2\delta\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi=\rho C_4^2x^{2\gamma+1}\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi",
    "12.66": r"U=C_5(J_s/\rho)^{1/2}x^{-1/2}F(y/x)",
    "12.67": r"-\overline{uv}=C_3U_{CL}^2G(y/x)=C_3C_5^2(J_s/\rho)x^{-1}G(y/x)",
    "12.68": r"\dot V\propto x^{1/2}",
    "12.69": r"\bar Y(x,y)=Y_{CL}(x)H(y/x)",
    "12.70": r"\dot M_s\cong\rho\int\bar YU\,dy",
    "12.71": r"\bar Y=C_6\big(\dot M_s/\sqrt{\rho J_s}\big)x^{-1/2}H(y/x)",
    "12.72": r"U=C_5U_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}F(y/x)",
    "12.73": r"\bar Y=C_7Y_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}H(y/x)",
    "12.74": r"\frac{\delta U'_{CL}}{U_{CL}}=C_8\big(\frac{\delta U'_{CL}}{U_{CL}}+\delta'\big)=C_9\frac\Psi{U_{CL}^2}",
    "12.75": r"0=-U\frac{\partial\bar e}{\partial x}-V\frac{\partial\bar e}{\partial y}-\overline{uv}\frac{\partial U}{\partial y}-\frac{\partial}{\partial y}\Big(\frac1{\rho_0}\overline{pv}+\frac12\overline{ev}\Big)-\bar\varepsilon",   # from analysis/ch12.md §2
    "12.76": r"\bar\tau=\mu\frac{\partial U}{\partial y}-\rho_0\overline{uv}",
    "12.77": r"\frac{\partial}{\partial x}P(x,y)-\frac{d}{dx}P(x,0)=-\rho\frac{\partial}{\partial x}\overline{v^2}=0",
    "12.78": r"U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}=-\frac1\rho\frac{\partial P}{\partial x}+\frac1\rho\frac{\partial\bar\tau}{\partial y}",
    "12.79": r"U=U(\rho,\tau_0,\nu,y)",
    "12.80": r"U^+=f(y^+)",
    "12.81": r"u_*^2\equiv\tau_0/\rho",
    "12.82": r"U^+=y^+",
    "12.83": r"U=U(\rho,\tau_0,\delta,y)",
    "12.84": r"\frac{U_\infty-U}{u_*}=F(y/\delta)=F(\xi)",
    "12.85": r"\frac{dU}{dy}=\frac{u_*^2}\nu\frac{df}{dy^+}",
    "12.86": r"-\frac{dU}{dy}=\frac{u_*}\delta\frac{dF}{d\xi}",
    "12.87": r"-\xi\frac{dF}{d\xi}=y^+\frac{df}{dy^+}",
    "12.88": r"U^+=\tfrac1\kappa\ln y^++B",
    "12.89": r"F(\xi)=-\frac1\kappa\ln(\xi)+A",
    "12.90": r"dP/dx=-2\tau_0/h",
    "12.91": r"dP/dx=-4\\tau_o/d",
    "12.92": r"\kappa B=1.6[\exp(0.1663B)-1]",
    "12.93": r"U^+=\tfrac1\kappa\ln(y/y_0)",
    "12.94": r"\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)",
    "12.95": r"\overline{u_iT'}=-\kappa_T\,\partial\bar T/\partial x_i",
    "12.96": r"\overline{u_iY'}=-\kappa_{mT}\,\partial\bar Y/\partial x_i",
    "12.97": r"\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-\frac1\rho\frac{\partial P}{\partial x_i}+\frac{\partial}{\partial x_j}\Big([\nu+\nu_T]\Big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\Big)-\frac23\bar e\delta_{ij}\Big)",
    "12.98": r"\nu_T\sim l_Tu_T",
    "12.99": r"0=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big(\nu\frac{\partial U}{\partial y}-\overline{uv}\big)=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big([\nu+\nu_T]\frac{\partial U}{\partial y}\big)",
    "12.100": r"0=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big(\nu\frac{dU}{dy}+\kappa^2y^2\big(\frac{dU}{dy}\big)^2\big)",
    "12.101": r"\frac U{u_*}\cong\frac1\kappa\ln y+\mathit{const.}",
    "12.102": r"Dw/Dt\sim g\alpha T'\sim g\Delta T/T",
    "12.103": r"\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\frac{\nu_T}{\sigma_e}\frac{\partial\bar e}{\partial x_j}\big)-\bar\varepsilon-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}",
    "12.104": r"\nu_T=C_\mu\bar e^2/\bar\varepsilon",
    "12.105": r"\frac{\partial\bar\varepsilon}{\partial t}+U_j\frac{\partial\bar\varepsilon}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\frac{\nu_T}{\sigma_\varepsilon}\frac{\partial\bar\varepsilon}{\partial x_j}\big)-C_{\varepsilon1}\big(\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}\big)\frac{\bar\varepsilon}{\bar e}-C_{\varepsilon2}\frac{\bar\varepsilon^2}{\bar e}",
    "12.107": r"\mathrm{Rf}=\frac{-g\alpha\overline{wT'}}{-\overline{uw}(dU/dz)}",
    "12.108": r"\\mathrm{Ri}=\\alpha g(d\\bar T/dz)/(dU/dz)^2",
    "12.109": r"\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}",
    "12.110": r"L_M\equiv-u_*^3/(\kappa\alpha g\overline{wT'})",
    "12.111": r"\mathrm{Rf}=z/L_M",
    "12.112": r"\frac{\partial}{\partial t}\Big(\frac12\overline{T'^2}\Big)+U\frac{\partial}{\partial x}\Big(\frac12\overline{T'^2}\Big)=-\overline{wT'}\frac{d\bar T}{dz}-\frac{\partial}{\partial z}\Big(\frac12\overline{T'^2w}-\kappa\frac{\partial\overline{T'^2}}{\partial z}\Big)-\bar\varepsilon_T$, $\bar\varepsilon_T=\kappa\overline{(\partial T'/\partial x_j)^2}",   # from analysis/ch12.md §2
    "12.113": r"S_T\propto\bar\varepsilon_T\bar\varepsilon^{-1/3}K^{-5/3}",
    "12.114": r"S_T\propto K^{-1}$ for $2\pi/\eta\ll K\ll2\pi/\eta_T",   # from analysis/ch12.md §2
    "12.115": r"\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha\frac{dX_\alpha}{dt}}",
    "12.116": r"\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha u_\alpha}=2\int_0^t\overline{u_\alpha(t')u_\alpha(t)}\,dt'",
    "12.117": r"\frac{d}{dt}\overline{X_\alpha^2}=2\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau",
    "12.118": r"\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\int_0^tdt'\int_0^{t'}r_\alpha(\tau)\,d\tau",
    "12.119": r"\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)r_\alpha\,d\tau",
    "12.120": r"\overline{X_\alpha^2}\simeq\overline{u_\alpha^2}t^2",
    "12.121": r"(X_\alpha)_{rms}=(u_\alpha)_{rms}t",
    "12.122": r"\overline{X_\alpha^2}=2\overline{u_\alpha^2}\Lambda_tt",
    "12.123": r"(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}",
    "12.124": r"\overline{R_n^2}=\overline{R_{n-1}^2}+L^2+2\overline{\mathbf R_{n-1}\cdot\mathbf L}",
    "12.125": r"(R_n)_{rms}=L\sqrt n",
    "12.126": r"\\nu=\\frac12\\frac{d\\sigma^2}{dt}",
    "12.127": r"D_T\equiv\frac12\frac{d}{dt}(\overline{X_\alpha^2})=\overline{u_\alpha^2}\int_0^tr_\alpha(\tau)\,d\tau",
    "12.128": r"D_T\cong\overline{u_\alpha^2}\,t$ for $t\ll\Lambda_t",   # from analysis/ch12.md §2
    "12.129": r"D_T\cong\overline{u_\alpha^2}\Lambda_t",
}

EQ.update({
    "12.7": r"\overline{\int_a^bu^m\,dt}=\int_a^b\overline{u^m}\,dt",
    "12.17": r"R_{11}(\tau)=\overline{u_1(t)u_1(t+\tau)}=R_{11}(-\tau)",
    "12.22": r"\overline{u_1^2}=\int_{-\infty}^{+\infty}S_e(\omega)\,d\omega",
    "12.39": r"\Lambda_f\equiv\int_0^\infty f(r)\,dr,\ \Lambda_g\equiv\int_0^\infty g(r)\,dr,\ \lambda_f^2\equiv-2/[d^2f/dr^2]_{r=0},\ \lambda_g^2\equiv-2/[d^2g/dr^2]_{r=0}",
    "12.41": r"R_{ij}=\overline{u^2}\Big\{f(r)\delta_{ij}+\frac r2\frac{df}{dr}\Big(\delta_{ij}-\frac{r_ir_j}{r^2}\Big)\Big\}",
    "12.42": r"\bar\varepsilon=2\nu\overline{S'_{ij}S'_{ij}}=\frac\nu2\overline{\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)^2}",
    "12.43": r"\bar\varepsilon=6\nu\Big\{\overline{\Big(\frac{\partial u_1}{\partial x_1}\Big)^2}+\overline{\Big(\frac{\partial u_1}{\partial x_2}\Big)^2}+\overline{\Big(\frac{\partial u_1}{\partial x_2}\Big)\Big(\frac{\partial u_2}{\partial x_1}\Big)}\Big\}=-15\nu\overline{u^2}\Big[\frac{d^2f}{dr^2}\Big]_{r=0}=30\nu\frac{\overline{u^2}}{\lambda_f^2}=15\nu\frac{\overline{u^2}}{\lambda_g^2}",   # p590 re-read
    "12.50": r"\eta=(\nu^3/\bar\varepsilon)^{1/4},\ u_K=(\nu\bar\varepsilon)^{1/4}",   # p595 re-read: two members only
    "12.61": r"U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y},\ 0\cong-\frac1\rho\frac{\partial}{\partial y}\big(P+\rho\overline{v^2}\big)",   # p602 re-read: both members
    "12.76": r"0=-\frac{\partial P}{\partial x}+\frac{\partial\bar\tau}{\partial y},\ 0=-\frac{\partial}{\partial y}\big(P+\rho\overline{v^2}\big)",   # p610 re-read: the definition of the total stress is the unnumbered line after it
    "12.54": r"S_{11}(k_1)=C_1\bar\varepsilon^{2/3}k_1^{-5/3}",
    "12.64": r"\frac{\delta U'_{CL}}{U_{CL}}=C_1,\ \frac{\delta U'_{CL}}{U_{CL}}+\delta'=C_2,\ \frac\Psi{U_{CL}^2}=C_3",
    "12.75": r"0=-U\frac{\partial\bar e}{\partial x}-V\frac{\partial\bar e}{\partial y}-\overline{uv}\frac{\partial U}{\partial y}-\frac{\partial}{\partial y}\Big(\frac1{\rho_0}\overline{pv}+\frac12\overline{u_i^2v}\Big)-\bar\varepsilon",
    "12.88": r"U^+=\frac1\kappa\ln(y^+)+B",
    "12.91": r"dP/dx=-4\tau_0/d",
    "12.93": r"U^+=\frac1\kappa\ln(y/y_0)",
    "12.106": r"\frac{\partial\bar e}{\partial t}+U\frac{\partial\bar e}{\partial x}=-\frac{\partial}{\partial z}\Big(\frac1{\rho_0}\overline{pw}+\frac12\overline{u_i^2w}\Big)-\overline{uw}\frac{\partial U}{\partial z}+g\alpha\overline{wT'}-\bar\varepsilon",
    "12.108": r"\mathrm{Ri}\equiv\frac{N^2}{(dU/dz)^2}=\frac{\alpha g(d\bar T/dz)}{(dU/dz)^2}",
    "12.112": r"\frac{\partial}{\partial t}\Big(\frac12\overline{T'^2}\Big)+U\frac{\partial}{\partial x}\Big(\frac12\overline{T'^2}\Big)=-\overline{wT'}\frac{d\bar T}{dz}-\frac{\partial}{\partial z}\Big(\frac12\overline{T'^2w}-\kappa\frac{\partial}{\partial z}\Big(\frac12\overline{T'^2}\Big)\Big)-\bar\varepsilon_T",
    "12.114": r"S_T\propto K^{-1}\ \text{for}\ 2\pi/\eta\ll K\ll2\pi/\eta_T",
    "12.126": r"\nu=\frac12\frac{d\sigma^2}{dt}",
    "12.128": r"D_T\cong\overline{u_\alpha^2}\,t\ \text{for}\ t\ll\Lambda_t",
    "12.129": r"D_T\cong\overline{u_\alpha^2}\Lambda_t\ \text{for}\ t\gg\Lambda_t",
    # earlier chapters (as their notebooks show them)
    "4.58": r"\varepsilon=2\nu S_{ij}S_{ij}",
    "5.13": r"\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega",
    "8.8": r"\tau_0=\frac a2\frac{dp}{dz}",
})
for _k, _v in EQ.items():
    assert "\\\\" not in _v and "$" not in _v, (_k, _v)
CHAPTER_PREFIX = r"12\."

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




# ---------------------------------------------------------------------------------------------------------------------
# Part F of the design, parsed at build time: every derivation copied word for word (the ch07–ch11 pattern; the ch12
# design writes the fields inline — "**Goal:** … · **Start:** … · **Plan:** • … • …" — and the steps as
# "n. **did** … · **tex** … · **why** … · **plain** …").
# ---------------------------------------------------------------------------------------------------------------------
_FIELD = re.compile(r"(?:^- |\s·\s+)\*\*(Goal|Start|Plan|Tools|Assumptions|Steps|Result|Check|sympy check intent|What it means|Traps)"
                    r"(?::\*\*|\*\*\s*(\([^)]*\))?:)\s*", re.M)


def _start_split(s: str) -> tuple[str, str]:
    """'maths … — "in words".' → (maths, words); without the quoted tail the whole text is the maths line."""
    m = re.search(r"\s—\s[\"“](.+?)[\"”]\.?\s*$", s)
    if m:
        return s[:m.start()].strip(), m.group(1).strip()
    return s.strip(), ""


def part_f() -> dict[str, dict]:
    """Parse Part F into {D01: dict(title, stars, goal, start, plan, tools, assumptions, steps, result, check, meaning,
    traps, check_src, intent)} — every field word for word."""
    part = DESIGN.split("## Part F", 1)[1].split("\n## Closing check", 1)[0]
    chunks = re.split(r"^### (D\d\d) · ", part, flags=re.M)
    out: dict[str, dict] = {}
    for key, body in zip(chunks[1::2], chunks[2::2]):
        code = re.search(r"```python\n(.*?)```", body, flags=re.S)
        check_src = textwrap.dedent(code.group(1)).strip("\n") if code else None
        body = re.sub(r"\s*`check_src` sketch:\s*```python\n.*?```", "", body, flags=re.S)
        body = re.sub(r"```python\n.*?```", "", body, flags=re.S)
        lines = body.splitlines()
        title = _abs_plain(re.split(r" — ★", lines[0])[0].strip())
        stars = lines[0].count("★")
        # the step list is cut out first (its lines start with "  n. **did**"), the rest are fields
        steps_raw: list[list[str]] = []
        rest: list[str] = []
        in_steps = False
        for ln in lines[1:]:
            if ln.startswith("---") or ln.startswith("## "):
                break
            if re.match(r"^- \*\*Steps:\*\*", ln):
                in_steps = True
                rest.append(ln)
                continue
            if in_steps and re.match(r"^\s{2}\d+\. \*\*did\*\*", ln):
                steps_raw.append([re.sub(r"^\s{2}\d+\. ", "", ln)])
                continue
            if in_steps and ln.startswith("- **"):
                in_steps = False
            if in_steps and steps_raw and ln.strip():
                steps_raw[-1].append(ln)
                continue
            rest.append(ln)
        text = "\n".join(rest)
        pieces = _FIELD.split("\n" + text.replace("\n- **", "\n- **"))
        # _FIELD has two groups: name, optional bracket; walk the split result
        fields: dict[str, str] = {}
        flat = _FIELD.split(text)
        k = 1
        while k < len(flat):
            name, bracket, val = flat[k], flat[k + 1], flat[k + 2]
            val = _abs_plain(_join(val.splitlines()))
            if bracket:
                val = f"{bracket} {val}"
            fields[name] = val.strip()
            k += 3
        parsed = []
        for st in steps_raw:
            s = _abs_plain(_join(st))
            bits = re.split(r"(?:^|\s·\s)\*\*(did|tex|why|plain)\*\*\s+", s)
            d = {name: val.strip().rstrip("·").strip() for name, val in zip(bits[1::2], bits[2::2])}
            assert set(d) == {"did", "tex", "why", "plain"}, (key, s[:120])
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=d["why"], plain=d["plain"]))
        start_tex, start_plain = _start_split(fields["Start"].rstrip(" ·"))
        res_tex, res_plain = _start_split(fields["Result"].rstrip(" ·"))
        plan = [p.strip().rstrip(".·").strip() + "." for p in fields["Plan"].split("•") if p.strip(" ·.")]
        tools = [t.strip().rstrip(".·").strip() for t in fields["Tools"].split(";") if t.strip(" ·.")]
        out[key] = dict(title=title, stars=stars, goal=fields["Goal"].rstrip(" ·"), start=(display_tex(start_tex), start_plain),
                        plan=plan, tools=tools, assumptions=fields.get("Assumptions", "").rstrip(" ·"), steps=parsed,
                        result=(display_tex(res_tex), res_plain), check=fields.get("Check", ""),
                        meaning=fields.get("What it means", ""), traps=fields.get("Traps", ""), check_src=check_src,
                        intent=fields.get("sympy check intent", ""), raw_start=start_tex, raw_result=res_tex)
    return out


PF = part_f()
assert len(PF) == 28, len(PF)
N_STEPS = sum(len(d["steps"]) for d in PF.values())
assert N_STEPS == 266, N_STEPS                              # the 266 steps of Part F, none lost by the parser

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
                    if not re.match(r"12\.", m.group(1)):
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


if "--dump" in sys.argv:
    for k, d in PF.items():
        print(f"=== {k} {d['title']}  ({len(d['steps'])} steps, {d['stars']} stars, check_src={bool(d['check_src'])})")
        print("  GOAL", d["goal"][:160])
        print("  START", d["start"][0][:300], "|", d["start"][1][:80])
        for i, s in enumerate(d["steps"], 1):
            print(f"  {i}. {s['did']} :: {s['tex'][:70]}  || why={len(s['why'].split())}w plain={len(s['plain'].split())}w")
        print("  RESULT", d["result"][0][:300], "|", d["result"][1][:80])
        print("  PLAN", len(d["plan"]), "TOOLS", len(d["tools"]), "ASSUME", d["assumptions"][:80])
        print("  CHECK", d["check"][:100], "| MEAN", d["meaning"][:60], "| TRAPS", d["traps"][:60], "| INTENT", d["intent"][:60])
    sys.exit(0)


# ---------------------------------------------------------------------------------------------------------------------
# Steps whose Part F *why* names only the rule: add why we make the move (nbkit needs both halves; Part F text is kept).
# ---------------------------------------------------------------------------------------------------------------------
pf_why_add("D05", 3, "The integral of a sum is the sum of the integrals, so the cosine part and the sine part can be looked at one at a time.")
pf_why_add("D05", 9, "Dividing the correlation by the variance turns the integral into the memory time, which is the quantity we want to read off the spectrum.")
pf_why_add("D06", 2, "Continuity holds for the total field in every realization; writing it with mean and fluctuation prepares it for averaging.")
pf_why_add("D14", 5, "We need this derivative for the second advection term, V times the cross-stream gradient of U.")
pf_why_add("D14", 6, "This is the first of the two advection terms of the thin-layer equation, now written with the similarity functions only.")
pf_why_add("D18", 2, "Five quantities (U, ρ, τ₀, ν, y) in three dimensions leave two independent dimensionless groups, so the law is one function of one variable.")
pf_why_add("D19", 8, "Dividing by a positive number keeps the equation true and isolates the slope, ready to be integrated.")
pf_why_add("D19", 12, "This is the definition of the friction velocity, so the ratio of wall stress to dynamic pressure becomes a ratio of velocities squared.")
pf_why_add("D21", 9, "For a very small wall distance the term with κ²y⁺² under the root is negligible, so only the viscous stress is left — the sublayer must come out.")
pf_why_add("D24", 3, "We only rename the vertical component and the dissipation, so that the budget is written with the symbols used in a stratified layer.")
pf_why_add("D25", 2, "We want Rf as a function of height alone, so every flux is replaced by its surface-layer value.")
pf_why_add("D26", 2, "The derivative of a square is twice the quantity times its derivative; we differentiate because the rate of spreading is easier to find than the spread itself.")
pf_why_add("D26", 4, "The rate of change of a particle's position is its velocity, which brings the velocity statistics into the equation.")


# ---------------------------------------------------------------------------------------------------------------------
# One-sentence reminders written for the way each earlier tool is used in THIS chapter (a concept missing here falls back
# to the gist stored in knowledge/primers.md), and a remind() that names the earlier chapter even when Part E gives
# several primer numbers or none.
# ---------------------------------------------------------------------------------------------------------------------
REMIND.update({
    "Gaussian distribution": r"the bell-shaped density $\propto e^{-x^2/2\sigma^2}$ with mean 0 and standard deviation $\sigma$; `rng.normal(0, sigma, n)` draws samples from it.",
    "sympy (symbols, integrate, simplify)": "algebra with symbols: `sp.symbols`, `sp.diff`, `sp.integrate`, `sp.simplify`; an expression that simplifies to 0 is an identity.",
    "quadratic formula": r"$a\lambda^2+b\lambda+c=0$ has the roots $\big(-b\pm\sqrt{b^2-4ac}\big)/(2a)$; they are real only if the discriminant $b^2-4ac\ge0$.",
    "substitution / change of variable in an integral or average": r"rename the variable ($t=t'+\tau$) and convert limits and differentials with it; a stationary average does not care what the start time is called.",
    "Taylor series to second order": r"near a point, $f(\tau)\approx f(0)+f'(0)\tau+\tfrac12f''(0)\tau^2$: value, slope and curvature.",
    "product rule": r"$(fg)'=f'g+fg'$.",
    "np.trapezoid; Simpson": "`np.trapezoid(y, x)` integrates sampled data with straight pieces (error ∝ step²); Simpson's rule uses parabolas (error ∝ step⁴).",
    "Euler's formula e^{iθ} = cos θ + i sin θ": r"$e^{i\theta}=\cos\theta+i\sin\theta$ with $i^2=-1$: a complex exponential is a cosine (even part) plus $i$ times a sine (odd part).",
    "np.fft.rfft, Fourier modes": r"any record is a sum of waves $e^{i\omega t}$; `np.fft.rfft(u)` returns their complex amplitudes at the frequencies `np.fft.rfftfreq(n, dt)`.",
    "np.isclose (the one-number twin of np.allclose)": "`np.isclose(a, b, rtol=…)` is True when two numbers agree within a tolerance; inside an `assert` it stops the cell if they do not.",
    "isotropic tensors (built from δ_ij and a vector)": r"a tensor with no preferred direction can only be built from $\delta_{ij}$ — and, if one vector $\mathbf r$ is given, from $r_ir_j$ as well, with coefficients that depend on the length $r$ alone.",
    "pandas DataFrame for small tables": "`pd.DataFrame(rows)` is a table with named columns; `.to_string(index=False)` prints it without the row numbers.",
    "separation argument: a function of x equal to a function of ξ is a constant": r"if a function of $x$ alone equals a function of $\xi$ alone for all $x$ and $\xi$, both must be the same constant.",
    "chain rule": r"$\frac{d}{dy}f(g(y))=f'(g)\,g'(y)$ — each change of variable brings the derivative of the inner function as a factor.",
    "scaled (dimensionless) variables": r"writing $y=l\,y^+$ turns $\partial/\partial y$ into $(1/l)\,\partial/\partial y^+$; each derivative brings one factor of the scale.",
    "iterated integrals": r"a double integral done one variable at a time, the inner one first with the outer variable held fixed.",
    "sinh": r"$\sinh x=(e^x-e^{-x})/2\approx x+x^3/6$ for small $x$ — needed once in this block, in the window average of N23.",
    "even and odd functions": r"even: $f(-x)=f(x)$ (cosine); odd: $f(-x)=-f(x)$ (sine); an odd function integrates to zero over a symmetric range.",
    "integral as a limit of sums; fundamental theorem of calculus": r"an integral is the limit of (value × width) sums, so it is linear; and $\int_a^bf'(x)\,dx=f(b)-f(a)$.",
})


def remind(cid: str, lead: str = "") -> None:   # noqa: F811  (replaces the ch11 helper: every line names its earlier chapter)
    """One 🔁 cell reminding the tools of earlier chapters that Part E lists as first used in ``cid`` — one sentence each,
    the concept text verbatim so that the ledger check matches."""
    items = [(r[0], r[2]) for r in LEDGER if "primers.md" in r[2] and r[1] == cid]
    if not items:
        return
    lines = []
    for concept, by in items:
        refs = re.findall(r"ch(\d\d) (P\d+a?)", by)
        where = ", ".join(f"Ch. {int(ch)}, {p}" for ch, p in refs) or "earlier chapters"
        if concept in REMIND:
            text = REMIND[concept]
        else:
            text = " ".join(_gist(p) for _ch, p in refs if _gist(p)) or "used here exactly as in the chapter named."
        shown = re.sub(r"\^\{([^{}]*)\}", r"^(\1)", concept)        # no raw TeX braces in the bold concept name
        lines.append(f"- **{shown}** — {text} *({where})*")
        nb.primers.append(f"{concept} (reminder)")
    head = "> 🔁 **Tools from earlier chapters used in this block** (one line each; the full primers are in the chapters named)"
    nb.md(f"{head}{(' — ' + lead) if lead else ''}\n\n" + "\n".join(lines), tags=["primer"])


_recap_plain = nb.recap


def _recap_with_id(rid: str, title: str, text: str, where: str):
    """A recap whose heading shows its id (Part F points to recaps as "(R01)")."""
    return _recap_plain(rid, f"{title} `{rid}`", text, where)


nb.recap = _recap_with_id
_D_ch11 = D


def D(key: str, ref: str = "", check_src: str | None = None, extra_check: str = "", after: str = "") -> None:   # noqa: F811
    """Part F writes Goal, Plan, What it means and Check as lower-case clauses; a sentence in the notebook starts with a
    capital. Nothing else is changed."""
    d = PF[key]
    def cap(s: str) -> str:
        """Capital first letter only when the text starts with an ordinary lower-case word of two or more letters
        (never a symbol such as λ_t, κ, σ², t or a maths span)."""
        return s[:1].upper() + s[1:] if re.match(r"[a-z]{2,}(?=[\s,;:.—–-])", s) else s
    for fld in ("goal", "assumptions", "check", "meaning", "traps"):
        d[fld] = cap(d[fld])
    d["plan"] = [cap(x) for x in d["plan"]]
    d["start"] = (d["start"][0], cap(d["start"][1]))
    d["result"] = (d["result"][0], cap(d["result"][1]))
    _D_ch11(key, ref=ref, check_src=check_src, extra_check=extra_check, after=after)


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=r"""
Stir milk into coffee and in two seconds the cup is uniform; without stirring, molecular diffusion would need hours. The
difference is turbulence — a flow so irregular that no one can predict where a given drop will be, yet so regular *on average*
that engineers design aircraft with it and climate models run on it. This chapter is about that "on average". We split every
field into a mean and a fluctuation, average the equations, and find that the fluctuations push back on the mean through one new
term — a correlation, the Reynolds stress — for which there is no equation. Everything else is a way of living with that gap:
measuring the fluctuations (correlations, spectra), following their energy from the big eddies that take it from the mean flow
down to millimetre eddies that turn it into heat (the cascade and Kolmogorov's −5/3 law), using symmetry and conservation where
they are enough (jets, the logarithmic law near a wall), guessing the missing term where they are not (eddy viscosity, mixing
length, k–ε), adding buoyancy (Richardson numbers, the Monin–Obukhov length — the vocabulary of every surface-flux scheme in a
climate model), and finally asking how turbulence spreads things (Taylor's dispersion).""",
    roadmap=[
        "**Averages and the one rule that fails** — an average passes through everything linear, not through a product (`C01`).",
        "**The autocorrelation and the memory time** — how long a signal remembers itself (`C02`).",
        "**The spectrum** — the same information, sorted by frequency (`C03`).",
        "**The Reynolds-averaged equations and the Reynolds stress** — the term the averaging leaves behind (`C04`).",
        "**Isotropic turbulence** — one function, one gradient (`C05`).",
        "**The two energy budgets** — and the term they share with opposite signs (`C06`).",
        "**Kolmogorov scales and the cascade** — viscosity sets where, not how much (`C07`).",
        "**The −5/3 law** — a spectrum from dimensional analysis alone (`C08`).",
        "**The self-similar jet** — one invariant fixes both exponents (`C09`).",
        "**Wall units and the law of the wall** (`C10`).",
        "**The logarithmic law** — what is left when neither length may matter (`C11`).",
        "**Eddy viscosity and the mixing length** — the simplest guess for the missing term (`C12`).",
        "**The k–ε model** — two transport equations for a velocity and a length (`C13`).",
        "**Flux and gradient Richardson numbers** — buoyancy against shear (`C14`).",
        "**The Monin–Obukhov length** — the height where buoyancy takes over (`C15`).",
        "**Taylor's dispersion** — why a plume grows like $t$ first and $\\sqrt t$ later (`C16`).",
    ],
    prerequisites=[
        r"index notation, $\delta_{ij}$, trace, symmetric contraction (Ch. 2); Lagrangian vs Eulerian descriptions (Ch. 3)",
        r"the Boussinesq equations and the dissipation $\varepsilon=2\nu S_{ij}S_{ij}$ (4.58) (Ch. 4); vortex stretching (Ch. 5)",
        "similarity solutions and the laminar jet (Ch. 8–9); thin-layer scaling (Ch. 9)",
        "the Π theorem, potential temperature and the two lapse-rate conventions (Ch. 1)",
        "the disturbance-energy budget and the Richardson number (Ch. 11)",
    ])
nb.explainer_index([
    ("reynolds_averaging_window", "How long must you average to get 'the mean'?",
     "a mean is only as good as its number of independent samples: N members, or Δt/Λ_t memory times"),
    ("correlation_and_spectrum", "Why are a correlation and a spectrum the same information?",
     "stretch the memory and the spectrum squeezes; its height at zero is the integral scale, its area the variance"),
    ("reynolds_stress_parcels", "How can fluctuations that average to zero push the mean flow?",
     "going up correlates with being slow: ⟨uv⟩ < 0 in a positive shear"),
    ("energy_cascade_spectrum", "What does a higher Reynolds number change?",
     "the supply ΔU³/L is fixed by the big eddies; viscosity only sets how far down the ladder the energy must go"),
    ("turbulent_energy_budget", "Where does turbulent energy come from and go?",
     "one term, two signs: the mean flow's loss is the turbulence's income"),
    ("turbulent_jet_similarity", "How does a jet spread without a turbulence model?",
     "one invariant + shape-preservation fix both exponents"),
    ("law_of_the_wall", "Where does the logarithm come from?",
     "it is what is left when the answer may depend on neither the viscous nor the outer length"),
    ("mixing_length_closure", "What does a closure constant do?",
     "κ sets the slope, the wall damping sets the intercept"),
    ("stratified_surface_layer", "When does stratification kill turbulence?",
     "Rf = z/L_M: below |L_M| shear rules, above it buoyancy"),
    ("taylor_dispersion", "Why t first and √t later?",
     "a random walk whose step is the particle's memory"),
])
nb.setup()
nb.code(r"""
import json, itertools, logging                          # small reference files; index combinations; Python's message system
import numpy as np                                       # arrays and maths (Ch. 1 primer P03)
import sympy as sp                                       # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                          # static figures (Ch. 1 primer P01)
import pandas as pd                                      # tables (DataFrame) for printed comparisons
from fractions import Fraction                           # exact fractions such as -5/3 (Ch. 1 primer P60)
from scipy import signal                                 # signal.welch: an averaged spectrum estimate (primer P292 below)
from scipy.integrate import quad, solve_ivp, cumulative_trapezoid   # 1-D integral (P87), ODE integrator (P31), running integral (P190)
from scipy.optimize import brentq                        # root of a function between two brackets (P108)
from fluidpy import ch12_turbulence as ch12              # the tested chapter-12 module (re-exports TS and WT)
from fluidpy.core import turbstats as TS, wall_turbulence as WT    # the two new core modules: statistics/spectra and wall layers
from fluidpy.core import dimensional as DIM, stratification as STRAT, laminar as LAM   # Ch. 1 Π theorem and N^2, Ch. 8 laminar flows
from fluidpy.core import boundary_layer as BL, jets as JET, diffusion as DIFF          # Ch. 9 boundary layers and jets, Ch. 8 diffusion
from fluidpy import ch11_instability as ch11             # Ch. 11: the Lorenz system and the Richardson-number criterion
from fluidpy.core.interact import slider_figure, animate_figure, live   # plotly sliders (P17) and live widgets (P47)
from fluidpy.core.anim import animate                    # matplotlib animations (P16); show_animation came with setup
from fluidpy.core.style import COLORS, savefig           # the house palette and a helper that saves PNGs to outputs/ch12
from fluidpy.core.thermo import G0                       # standard gravity 9.80665 m/s^2 (the default g of every function)
sys.path.insert(0, str(ROOT / "scripts"))                # the chapter's drawing helpers live in scripts/
from ch12_drawings import (scales_sketch, parcel_sketch, stress_element_sketch, fg_geometry_sketch,   # drawing only:
                           free_shear_sketch, wall_layers_sketch, surface_layer_sketch, plume_sketch)   # no physics lives there
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless "font not found" notes
plt.rcParams["axes.titlesize"] = 10                      # this chapter's figure titles are whole sentences: a smaller title font keeps them inside the figure
C_MEAN, C_FLUC, C_RS, C_VISC = COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["rose"]   # mean, fluctuation, Reynolds stress, viscous
C_BUOY, C_SCAL, C_REF = COLORS["blue"], COLORS["amber"], COLORS["muted"]                             # buoyancy/temperature, scalar, reference
KAPPA, B_LOG = 0.41, 5.0                                 # a common textbook pair for the log law (illustrative; C11 fits its own)
REF = ROOT / "reference" / "ch12"                        # public benchmark data with their sources (SOURCES.md)
print(len([n for n in dir(ch12) if not n.startswith("_")]), "public names in fluidpy.ch12_turbulence")
print("ch12.log_law is WT.log_law:", ch12.log_law is WT.log_law)   # one function, two names
""", explain=r"""
1. Numerical, symbolic and plotting libraries (primed in Ch. 1); `quad`, `solve_ivp`, `brentq`, `cumulative_trapezoid` are
   reminded where they are first used.
2. `ch12` is the chapter module. `TS` (`core/turbstats.py`: averages, correlations, spectra, synthetic signals) and `WT`
   (`core/wall_turbulence.py`: wall units, log law, composite profiles) are the two new core modules; `ch12` **re-exports both**,
   so `ch12.log_law` and `WT.log_law` are the same function (the last line prints `True`). Every function cites its § and equation
   and is tested in `tests/test_ch12.py`.
3. Earlier chapters' functions are reused, not rewritten: the Π theorem (`DIM`), the buoyancy frequency (`STRAT`), laminar
   channel and pipe flow (`LAM`), boundary layers (`BL`), the laminar jet (`JET`), diffusion (`DIFF`), the Lorenz system (`ch11`).
4. `scripts/ch12_drawings.py` only draws (eddies of many sizes, a parcel in a shear, a stress element, the layers near a wall,
   a plume) — no physics lives there.
5. `C_MEAN, C_FLUC, …` fix **one colour per meaning** for the whole chapter: mean purple, fluctuation and turbulent energy teal,
   Reynolds stress and shear production orange, viscous stress and dissipation rose, buoyancy and temperature blue, scalar amber,
   laminar or reference curves grey.
6. `KAPPA, B_LOG = 0.41, 5.0` is a common textbook pair for the logarithmic law, used as a **labelled illustrative** default;
   `REF` points to the public benchmark files.""")
nb.md(r"""
## ⚠️ Conventions in this chapter (read once, come back when a symbol surprises you)

**$u$ is a fluctuation here.** A tilde marks the total field, a capital letter or an over-bar the mean, a lower-case letter or a
prime the fluctuation: $\tilde u_i=U_i+u_i$ (12.24). So "$u$" is *not* "the velocity" in this chapter.

**One letter, several meanings.** The book reuses symbols between sections; this is what each means and what we write.

| Book symbol | Meanings in this chapter (and earlier) | What we write | Code name |
|---|---|---|---|
| $\kappa$ | von Kármán constant (log law, Monin–Obukhov length); thermal diffusivity (mean heat equation); $\kappa_m$, $\kappa_T$, $\kappa_{mT}$ nearby | $\kappa$ = von Kármán; $\kappa_{th}$ = molecular thermal diffusivity in our own lines; $\kappa_T$ eddy diffusivity | `kappa`, `kappa_th`, `kappa_T` |
| $k$, $K$ | thermal conductivity; wavenumber $k_1$; 3-D wavenumber $K$; the "k" of k–ε | $k_{th}$ conductivity; $k_1$, $K$ wavenumbers; $\bar e$ turbulent kinetic energy | `k_th`, `k1`, `K`, `e` |
| $e$, $E$, $\varepsilon$ | $\bar e=\tfrac12\overline{u_i^2}$ turbulent kinetic energy (Ch. 1: internal energy); $\bar E=\tfrac12U_i^2$ mean-flow energy; $\bar\varepsilon$ dissipation, $\bar\varepsilon_T$ its thermal twin | the same, always with the bar | `e`, `E_mean`, `eps`, `eps_T` |
| $\lambda$, $\Lambda$ | Taylor microscales $\lambda_t,\lambda_f,\lambda_g$ (earlier: wavelength, eigenvalue); integral scales $\Lambda_t,\Lambda_f,\Lambda_g$ | always with the subscript | `lambda_t`, `Lambda_t`, … |
| $\eta$ | Kolmogorov length (earlier: similarity variable, surface elevation) | $\eta$; $\eta_T$ Batchelor scale | `eta`, `eta_T` |
| $f$, $g$, $F$, $G$ | correlation coefficients $f(r)$, $g(r)$; the wall function $U^+=f(y^+)$; gravity $g$; jet profiles $F$, $G$; defect function $F$ | $f(r)$, $g(r)$ in §12.6; $f(y^+)$ for the wall function in §12.9, as in the book (each block says which); $F$, $G$ jet; $F_d$ defect; $f_D$ Darcy friction factor | `f`, `g_corr`, `F`, `G`, `fD` |
| $R_{ij}$, $r$ | correlation tensor (Ch. 2–3: rotation tensor); separation $r$; correlation coefficients $r_{11}$, $r_\alpha$ | the same, said in words where both meet | `R`, `r` |
| $\tau$ | time lag; stress $\bar\tau$, $\tau_0$ | $\tau$ lag; $\tau_c$ memory time of our test signals; $\tau_d$ a decay time; $\bar\tau$, $\tau_0$ stress | `lag`, `tau_c`, `tau0` |
| $\delta$, $\theta$, $\alpha$ | jet width, boundary-layer thickness, $\delta_{ij}$; potential temperature; thermal expansion and the no-sum index of §12.12 | each named at first use; "component α (no sum)" said in words | `delta`, `theta`, `alpha` |
| $S$, $N$, $L$, $\sigma$ | spectra $S_e,S_{11}$ vs strain $\bar S_{ij},S'_{ij}$; $N$ realizations vs buoyancy frequency; outer scale $L$, $L_M$; Gaussian width $\sigma$ | as in the book, each named at first use | `n_members`, `N2`, `L`, `L_M` |
| $h$, $d$ | $h$ = **full** channel height; $d$ slot width, nozzle or pipe diameter | $h$ (full), $\delta=h/2$ | `h`, `delta`, `d` |

**Normalisations.**

| Quantity | Convention here |
|---|---|
| spectra | two-sided in angular frequency [rad/s] or wavenumber [rad/m], the factor $1/2\pi$ in the forward transform, area = variance; "one-sided" doubles the density |
| $\overline{u^2}$ in §12.6 | **one** velocity component, so $\bar e=\tfrac32\overline{u^2}$ |
| skewness, kurtosis | the book's are raw central moments (units $u^3$, $u^4$); `normalized=True` gives the usual ones (Gaussian 0 and 3) |
| Reynolds numbers | always named: $\mathrm{Re}_L=\Delta UL/\nu$, $R_\lambda$ (with which $\lambda$), $\mathrm{Re}_x$, $\mathrm{Re}_\tau=\delta^+$ |
| coordinates | $y$ from the wall; jet variable $\xi=y/x$ (width $\delta=x$ by convention, virtual origin dropped) |

**Three assumptions that switch silently.** (1) §12.5–12.7 use the Boussinesq equations with buoyancy; §12.8–12.10 use constant
density ($\rho_0\to\rho$, $P$ = deviation from hydrostatic); buoyancy returns in §12.11. (2) The theory is written with
ensemble averages, measurements are time averages (they agree for a stationary signal — *ergodicity*, primed in C01). (3) Every
sampled number is an estimate: cells fix their random seed and compare within 5 standard errors.

**Which Γ? Two conventions for one lapse rate.** An atmosphere cooling 6.5 K per km of height:

| Convention | Definition | Our example | Adiabatic value | "Stable" means |
|---|---|---|---|---|
| Kundu (the code computes with this) | $\Gamma\equiv dT/dz$ | $-6.5$ K/km | $\Gamma_a\approx-9.8$ K/km | $dT/dz>\Gamma_a$: $-6.5>-9.8$ ✓ |
| meteorology | $\Gamma_{met}\equiv-dT/dz$ | $+6.5$ K/km | $\Gamma_d\approx+9.8$ K/km | $\Gamma_{met}<\Gamma_d$: $6.5<9.8$ ✓ |

To convert: negate the number and flip the inequality (Ch. 1, P48). Both say **stable**. In every buoyancy term of this chapter
$\bar T$ and $T'$ are **potential** temperature; with the thermometer temperature
$N^2=g\alpha(dT/dz-\Gamma_a)=g\alpha(\Gamma_d-\Gamma_{met})$.""")
nb.code(r"""
pd.set_option("display.max_colwidth", 70, "display.width", 200)   # let long table cells show
display(ch12.conventions())                              # the symbol table above, as the module itself records it
slips = pd.DataFrame(ch12.book_slips())[["id", "where", "printed", "corrected"]]   # the printed slips we teach in corrected form
print(len(ch12.book_slips()), "printed slips are taught in corrected form in this notebook")   # the count comes from the module, not from this text
display(slips)                                           # slip number · where in the book · what is printed · the correct form
""", explain=r"""
1. `ch12.conventions()` returns the module's own record of the symbol table (a pandas DataFrame — a table with named columns).
2. `ch12.book_slips()` lists **the printed slips** of the chapter — the cell prints how many (a wrong exponent, a wrong index, a swapped condition, a wrong sign…). We teach
   each in its corrected form, with a short ⚠️ box (`slip #k`) at the place where it is used. The most visible one: the book prints
   Kolmogorov's law with the exponent $+5/3$; the correct form is $S_{11}(k_1)=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ (12.54).""")
nb.md(r"""
> 🔁 **Tools from earlier chapters used in this one.** Each block reminds the ones it uses in one line where they first appear.
> In short: seeded random generators, the Gaussian and the $1/\sqrt N$ scatter (Ch. 1) carry every statistic; partial derivatives,
> the product and chain rules, Taylor series and integration by parts (Ch. 1, 9) carry the derivations; Fourier modes and the FFT
> (Ch. 7, 10) the spectra; the Π theorem (Ch. 1) and scaled variables (Ch. 8) the wall layer; the summation convention,
> $\delta_{ij}$ and isotropic tensors (Ch. 2) the averaged equations; `quad`, `brentq`, `solve_ivp`, `cumulative_trapezoid`
> (Ch. 1–9) compute; `slider_figure`, `animate` and `show_viz` (Ch. 1) draw.""")

# =====================================================================================================================
# A.1 §12.1 Introduction — N01 N02 N03 N04
# =====================================================================================================================
nb.section("12.1", "Introduction", intro=r"""
**What is this section about?** What "turbulent" means: not just "irregular", but five properties that come together. Then the
plan — since no one can predict a single turbulent flow in detail, we predict its statistics.""")
note("N01 [B]", r"""
**Five marks of turbulence.**

| Mark | What it means | Where we put a number on it |
|---|---|---|
| fluctuations | irregular in space and time, even with steady boundary conditions | C01–C03 |
| nonlinearity | appears above a critical Reynolds or Rayleigh number or below a critical Richardson number (Ch. 11); the nonlinear term makes the Reynolds stress | C04 |
| vorticity | three-dimensional fluctuating vorticity in eddies of many sizes, kept alive by vortex stretching, the first term on the right of $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13) (there $\mathbf u$ is the full velocity) | C07 |
| dissipation | it dies without a continuous supply of energy | C06 |
| diffusivity | rapid mixing of momentum, heat and matter | C12, C16 |

Random waves on a pond are irregular, but they neither dissipate much nor mix: they are not turbulence.""")
P("P286", "random-phase synthetic fields (np.fft.ifftn)", r"""
A field that *looks* turbulent can be built by giving every Fourier mode a chosen amplitude and a random phase, then
transforming back with the inverse FFT (`np.fft.ifft2` in two dimensions, `np.fft.ifftn` in any number). Building the velocity
from a stream function makes it divergence-free. Such a field has the spectrum we chose but **no dynamics — no cascade** — so
every caption says "kinematic".""", code=r"""
rng = np.random.default_rng(0)                 # seeded generator: the same field every run
phase = np.exp(2j*np.pi*rng.random((8, 8)))    # unit amplitude, a random phase for each of the 8 x 8 modes
field = np.fft.ifft2(phase).real               # back to physical space; keep the real part
print(field.shape, round(float(field.std()), 3))   # an 8 x 8 field and its standard deviation
""")
note("N02 [B]", r"""
**A turbulent velocity field still obeys mass conservation**, $\nabla\cdot\mathbf u=0$; white noise does not. The figure
compares the two: `TS.synthetic_solenoidal_field` (random phases, built from a stream function) and `TS.white_noise_field`
(an independent random number at every grid point). The first is built from a prescribed *energy spectrum* $E(K)$: a function that
says how much kinetic energy the eddies of size about $2\pi/K$ hold, $K$ being the wavenumber [rad/m] (spectra are taught in C03
and C08).""")
fig(r"""
nfield = 96 if FAST else 128                                   # grid points per side
E_raw = lambda K: K**4*np.exp(-(K/20.0)**2)                    # shape of the energy spectrum E(K) we prescribe (peak near K = 28 rad/m)
E_norm = quad(E_raw, 0, np.inf)[0]                             # its integral, used to normalise the shape
spec = lambda K: E_raw(K)/E_norm                               # now ∫E dK = 1 m²/s² = ½ mean(u² + v²): each component has variance ≈ 1
us, vs = TS.synthetic_solenoidal_field(nfield, 1.0, spectrum=spec, seed=1)   # divergence-free field in a 1 m box (kinematic)
uw, vw = TS.white_noise_field(nfield, seed=1)                  # white noise: every point independent
dxf = 1.0/nfield                                               # grid spacing [m]
div_s = ch12.divergence_rms(us, vs, dxf)/(np.sqrt(np.mean(us**2 + vs**2))/dxf)   # rms divergence in units of (speed / dx)
div_w = ch12.divergence_rms(uw, vw, dxf)/(np.sqrt(np.mean(uw**2 + vw**2))/dxf)   # the same for the noise
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4.0))
xg = (np.arange(nfield) + 0.5)*dxf                             # cell positions [m]
sk = nfield//24                                                # draw every sk-th arrow
for ax, (uu_, vv_), name, dv in ((a1, (us, vs), "synthetic eddies (kinematic)", div_s), (a2, (uw, vw), "white noise", div_w)):
    ax.pcolormesh(xg, xg, np.hypot(uu_, vv_), cmap="Purples", shading="auto")       # colour = speed
    ax.quiver(xg[::sk], xg[::sk], uu_[::sk, ::sk], vv_[::sk, ::sk], color=COLORS["ink"], width=0.004)   # arrows = velocity
    ax.set(title=f"{name}\nrms divergence = {dv:.1e} × (speed/dx)", xlabel="x [m]", ylabel="y [m]", aspect="equal")
print(f"relative rms divergence: synthetic {div_s:.1e}, white noise {div_w:.2f}")
""",
    see="Left: swirls of several sizes, with arrows that follow each other round. Right: salt-and-pepper, neighbouring arrows "
        "unrelated. The titles give the rms divergence in units of (rms speed)/(grid spacing): round-off on the left, of order one "
        "on the right.",
    read="Eddies are what mass conservation looks like in a random field: fluid that enters a region must leave it, so the "
         "velocity has to curl round. The left field is *kinematic* — it has a prescribed spectrum but no dynamics, no cascade.",
    change="…we move the spectrum's peak to higher $K$ (replace 20 by 40 in `spec`): the eddies shrink, and the field is still "
           "divergence-free to round-off.")
note("N03 [C]", r"""
**The largest eddies are as big as the flow.** In a turbulent boundary layer the largest eddies have the size of the layer,
$l\sim\delta$, and the edge between turbulent and irrotational fluid is ragged, not smooth. The outer scale $L$ returns in C07,
entrainment across that ragged edge in C09, the thickness $\delta$ in C10.""")
note("N04 [C]", r"""
**Scope.** Incompressible flow, no Coriolis force, three-dimensional fluctuations. Rotating stratified flows have nearly
two-dimensional (geostrophic) turbulence, which sends energy to *larger* scales — the opposite of the cascade of C07 — and
waits for Ch. 13.""")

# =====================================================================================================================
# A.2 §12.2 Historical Notes — N05, N06 named
# =====================================================================================================================
nb.section("12.2", "Historical Notes", intro=r"""
**What is this section about?** Who found what — read as a map of the chapter.""")
note("N05 [C]", r"""
**A timeline, by block.**

| Years | Who | Idea | Block |
|---|---|---|---|
| 1883–1895 | Reynolds | the pipe experiment; averaging the equations | C04 |
| 1915–1921 | Taylor | eddy transport; dispersion by continuous movements | C16 |
| 1922–1926 | Richardson | the cascade ("big whirls have little whirls"); the four-thirds law | C07, C16 |
| 1925 | Prandtl | the mixing length | C12 |
| 1930 | von Kármán | the logarithmic velocity profile | C11 |
| 1935 | Taylor | isotropic turbulence, correlations and spectra | C05 |
| 1941 | Kolmogorov, Obukhov | universal small scales, the −5/3 spectrum | C07, C08 |
| 1954 | Monin, Obukhov | the stratified surface layer | C15 |""")
nb.md(r"""
**N06** (Richardson's four-thirds law, an eddy diffusivity that grows with the size $l$ of the cloud,
$K\sim\bar\varepsilon^{1/3}l^{4/3}$) is only named here; it is stated with a number at the end of C16, where both of its
ingredients — the dissipation rate $\bar\varepsilon$ and the eddy diffusivity — exist.""")


# =====================================================================================================================
# A.3 §12.3 Nomenclature and Statistics for Turbulent Flow — C01
# =====================================================================================================================
nb.section("12.3", "Nomenclature and Statistics for Turbulent Flow", intro=r"""
**What is this section about?** What "the mean" of a turbulent signal is, how to compute it from many runs or from one long
record, which operations an average passes through — and the one it does not, which is where the rest of the chapter comes from.""")
core("C01", r"Ensemble averages $\langle u^m\rangle=\lim_{N\to\infty}\frac1N\sum_{n=1}^N(u(\mathbf x,t{:}n))^m$ (12.1) and the rules of averaging",
     "What is 'the mean' of a signal that never repeats, and what may an average be moved through?")
problem(r"""
A wind gauge on a mast reads 4.1, 6.3, 3.8, 5.5 m/s within four seconds. What is "the wind"? A climate normal is a 30-year
average for the same reason: the weather never repeats, so we describe it by its mean and by how it scatters about the mean. In
a laboratory you can repeat an experiment a hundred times and average the hundred runs at the same instant; outdoors you have
one run and can only average over time. When do the two agree?""")
idea(r"""
realization 1:  ~~~/\~~\/~~~      average ACROSS runs at one time   → ensemble average (theory uses this)
realization 2:  ~\/~~~/\~~~~      average ALONG one run over Δt     → time average (measurements use this)
realization N:  ~~~\/~/\~~~~      they agree when: stationary, and Δt ≫ memory, Δt ≪ drift
an average is a SUM (÷N):  it passes through +, ×constant, ∂/∂t, ∂/∂x, ∫   — all linear
                           it does NOT pass through a product:  mean(ũṽ) = Ū V̄ + mean(uv)
""", r"""
Symbols: $u(\mathbf x,t{:}n)$ is the value at position $\mathbf x$ and time $t$ **in realization $n$** (read "$t{:}n$" as "time
$t$ in run $n$"); $N$ is the number of runs; an over-bar is the average over the $N$ runs and $\langle\ \rangle$ its limit
$N\to\infty$; $m$ is a power (1 for the mean, 2 for the mean square).""")
remind("C01")
P("P280", "random variable, probability density and histogram", r"""
A *random variable* takes a different value in every realization. Its *probability density* says how often each value occurs
(area under the density between two values = the fraction of runs that land there). A *histogram* (`np.histogram`, `ax.hist`)
estimates the density by counting samples in bins. A Gaussian (bell-shaped) density has skewness 0 and normalised kurtosis 3
(both defined in N22 below).""", code=r"""
rng = np.random.default_rng(0)                    # seeded generator (Ch. 1, P10)
x = rng.normal(0.0, 1.0, 10_000)                  # 10 000 samples of a Gaussian random variable, mean 0, std 1
print(round(float(x.mean()), 3), round(float(x.std()), 3))   # sample mean ≈ 0 and standard deviation ≈ 1
print(np.histogram(x, bins=5, range=(-2.5, 2.5))[0])         # counts in 5 bins: most samples near the middle
""")
P("P283", "Ornstein–Uhlenbeck signal (a Langevin equation)", r"""
Our test signal throughout the chapter: a random signal with a chosen standard deviation $\sigma$ and a chosen *memory time*
$\tau_c$. Each step keeps a fraction $e^{-\Delta t/\tau_c}$ of the old value and adds fresh Gaussian noise $\xi_n$:
$u_{n+1}=u_ne^{-\Delta t/\tau_c}+\sigma\sqrt{1-e^{-2\Delta t/\tau_c}}\,\xi_n$. The square root is chosen so that the variance
stays $\sigma^2$; the update is exact for any step $\Delta t$. (It is the discrete form of a *Langevin equation*, a relaxation
toward zero driven by random kicks.) `TS.make_ensemble` builds many such runs around a mean of your choice.""", code=r"""
rng = np.random.default_rng(1)                    # seeded generator
sig, tau_c, dt_ = 1.0, 0.5, 0.05                  # standard deviation [m/s], memory time [s], time step [s]
keep = np.exp(-dt_/tau_c)                         # the fraction of the old value kept in one step
u_ou = np.empty(20_000); u_ou[0] = sig*rng.normal()          # start from a random value with the right variance
for n in range(u_ou.size - 1):                    # march in time
    u_ou[n + 1] = u_ou[n]*keep + sig*np.sqrt(1 - keep**2)*rng.normal()   # keep part of the past, add fresh noise
print(round(float(u_ou.std()), 2))                # ≈ 1: the variance is preserved
""")
P("P281", "standard error of a mean", r"""
The mean of $N$ **independent** samples scatters about the true mean with standard deviation $\sigma/\sqrt N$ — the *standard
error*. "Within 5 standard errors" is our test for every sampled number in this chapter (a correct code fails it about once in
two million tries). Samples taken close together in time are **not** independent. For a record of length $T$ whose memory time
is $\Lambda_t$ (the area under its autocorrelation, defined in C02) the variance of the time mean is $2\sigma^2\Lambda_t/T$ when
$T\gg\Lambda_t$ — the same triangle-weighted integral that C16 derives for dispersion — so the record is worth
$N=T/(2\Lambda_t)$ independent samples: one for every two memory times.""", code=r"""
rng = np.random.default_rng(2)                    # seeded generator
batch_means = rng.normal(0.0, 1.0, (2000, 64)).mean(axis=1)   # the means of 2000 batches of 64 Gaussian samples
print(round(float(batch_means.std()), 3), 1/np.sqrt(64))      # their scatter ≈ 0.125 = 1/sqrt(64)
""")
note("N07 [B]", r"""
**Realization, ensemble, ensemble average.** One run of an experiment is a *realization*; the set of all runs under the same
conditions is the *ensemble*; the over-bar is the average over $N$ of them, and $\langle\ \rangle$ its limit $N\to\infty$ (the
*expected value*). Chapter 11's Lorenz system showed why we must think this way: two runs that start $10^{-8}$ apart end up
completely different, so only statistics are reproducible.""",
     equation=r"\langle u^m(\mathbf x,t)\rangle=\lim_{N\to\infty}\frac1N\sum_{n=1}^N\big(u(\mathbf x,t{:}n)\big)^m", ref="12.1")
note("N19 [B]", r"""
**The mean** is the $m=1$ case with a finite number of runs — what the over-bar means from here on.""",
     equation=r"\overline{u(\mathbf x,t)}\equiv\frac1N\sum_{n=1}^Nu(\mathbf x,t{:}n)", ref="12.10")
code(r"""
t = np.arange(0.0, 20.0, 0.01)                    # time axis [s]: 2000 samples, 0.01 s apart
mean_fn = lambda t: np.exp(-t/10.0)               # the true (decaying) mean we hide inside the noise [m/s]
ens_u = TS.make_ensemble(64, t, mean_fn, sigma=0.3, tau_c=0.5, seed=0)   # 64 runs: mean + OU noise (std 0.3 m/s, memory 0.5 s)
ens_v = TS.make_ensemble(64, t, mean_fn, sigma=0.3, tau_c=0.5, seed=1) + 0.5*(ens_u - mean_fn(t))   # a second variable, partly correlated with the first
ubar = TS.ensemble_average(ens_u)                 # Eq. (12.10): average over the member axis (axis 0), one value per time
print(ens_u.shape, ubar.shape)                    # (64 members, 2000 times) -> (2000 times,)
print(f"rms error of the 64-member mean: {np.sqrt(np.mean((ubar - mean_fn(t))**2)):.4f} m/s  (sigma/sqrt(64) = {0.3/8:.4f})")
""", explain=r"""
1. `TS.make_ensemble(n_members, t, mean_fn, sigma, tau_c, seed)` returns an array with **the member index on axis 0**: row $n$ is
   realization $n$.
2. `ens_v` is a second signal that shares half of `ens_u`'s fluctuation, so the two are correlated (we need that for the product
   rule below).
3. `TS.ensemble_average` is the mean over axis 0, $\overline{u(\mathbf x,t)}\equiv\frac1N\sum_{n=1}^Nu(\mathbf x,t{:}n)$ (12.10).
   With 64 members its rms distance from the true mean is close to $\sigma/\sqrt{64}=0.0375$ m/s.""")
note("N08 [B]", r"""
**Stationary** means the statistics do not depend on the time origin. Then a *time* average over a window $\Delta t$ can replace
the ensemble average.""",
     equation=r"\overline{u^m(\mathbf x)}=\frac1{\Delta t}\int_{t-\Delta t/2}^{t+\Delta t/2}u^m(\mathbf x,t)\,dt", ref="12.2")
P("P284", "sliding-window average with np.cumsum", r"""
The sum over a window equals the difference of two cumulative sums, so a sliding average costs one pass over the data.
`np.cumsum(x)` returns the running sum $x_0,\ x_0+x_1,\ \dots$. Within half a window of each end there is no full window:
`TS.time_average` returns NaN there ("not a number", a marker for "no value").""", code=r"""
x = np.arange(10.0)                               # 0, 1, ..., 9
c = np.concatenate([[0.0], np.cumsum(x)])         # running sum with a leading zero: c[k] = x[0] + ... + x[k-1]
print((c[3:] - c[:-3])/3)                         # mean of every 3 consecutive values: 1, 2, ..., 8
""")
note("N09 [B]", r"""
**Homogeneous** is "stationary in space": the statistics do not depend on where you are. Then a *volume* average replaces the
ensemble average. One number: the mean square speed of the synthetic field of §12.1.""",
     equation=r"\overline{u^m(t)}=\frac1V\int_Vu^m(\mathbf x,t)\,dV", ref="12.3")
code(r"""
print(f"volume average of u^2 over the box: {TS.volume_average(us, m=2):.4g} (m/s)^2")   # Eq. (12.3) with m = 2 on the N02 field
print(f"the same with np.mean:              {np.mean(us**2):.4g} (m/s)^2")               # on a uniform grid it is the plain mean
""")
note("N10 [B]", r"""
**Stationary and non-stationary records side by side.** The figure shows two members of an ensemble with a constant mean and two
members of our decaying ensemble; the ensemble mean is dashed purple, the members teal.""")
fig(r"""
ens_c = TS.make_ensemble(64, t, lambda t: 0.6 + 0*t, sigma=0.3, tau_c=0.5, seed=2)   # a stationary ensemble: constant mean 0.6 m/s
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
for ax, ens, name in ((a1, ens_c, "stationary: the mean does not change"), (a2, ens_u, "non-stationary: the mean decays")):
    ax.plot(t, ens[0], color=C_FLUC, lw=0.7, label="member 1")                 # one realization
    ax.plot(t, ens[1], color=C_FLUC, lw=0.7, alpha=0.45, label="member 2")     # another one: different in detail
    ax.plot(t, TS.ensemble_average(ens), "--", color=C_MEAN, lw=1.8, label="ensemble mean (64 members)")   # Eq. (12.10)
    ax.set(title=name, xlabel="time t [s]")
a1.set_ylabel("u [m/s]")
a1.legend(fontsize=8)
""",
    see="Left: two wiggly records around a level dashed line. Right: two records around a dashed line that decays from 1 to "
        "about 0.14 m/s. In both panels the two members differ in every detail but share the dashed mean.",
    read="The ensemble mean is defined at every instant, whether or not the flow is stationary. A time average only makes sense "
         "on the left — or on the right if the window is short against the 10 s decay and long against the 0.5 s memory.",
    change="…the decay time were shorter than the memory time: no window could separate the mean from the fluctuation, and only "
           "an ensemble of repeated runs would do.")
P("P282", "ergodicity", r"""
For a stationary signal, a long time average of **one** run equals the ensemble average over **many** runs. This is called
*ergodicity*. The book assumes it silently every time a measurement (a time average) is compared with theory (ensemble
averages).""", code=r"""
t_long = np.arange(0.0, 2000.0, 0.05)                                 # one long record: 2000 s
one_run = TS.make_ensemble(1, t_long, lambda t: 0*t, 1.0, 0.5, seed=5)[0]   # a single OU run, true mean 0, std 1, memory 0.5 s
many = TS.make_ensemble(2000, np.arange(0.0, 5.0, 0.05), lambda t: 0*t, 1.0, 0.5, seed=6)   # 2000 short runs
print(f"time mean of one long run:        {one_run.mean():+.3f}  (standard error ≈ sqrt(2 × 0.5/2000) = 0.022)")   # ≈ 0 within a few standard errors
print(f"ensemble mean of 2000 runs at t = 2.5 s: {many[:, 50].mean():+.3f}  (standard error 1/sqrt(2000) = 0.022)")   # ≈ 0 as well
""")
note("N11 [B]", r"""
**Choosing a window.** Outdoors we average over time (or space): the window must be long against the memory of the signal and
short against the drift of its mean, $t_c\ll\Delta t\ll$ drift time. **Climate hook:** a 30-year normal is such a window; an
"eddy flux" $\overline{v'T'}$ in the general circulation is a covariance of deviations from a time or zonal mean (Ch. 13).""")
NOTES_IN["D01"] = (r"**N12 [B]** the mean of a sum, $\overline{u^m+v^m}=\overline{u^m}+\overline{v^m}$ (12.4) (step 2) · **N13 [B]** a "
                   r"constant factors out, $\overline{Au^m}=A\overline{u^m}$ (12.5) (step 3) · **N14 [B]** averaging and $\partial/\partial t$ "
                   r"swap, $\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$ (12.6) — ⚠️ exact for ensemble averages; for a "
                   r"finite time window only when the window is short against the drift; used again in C16 (step 4) · **N15 [B]** "
                   r"$\overline{\int_a^bu^m\,dt}=\int_a^b\overline{u^m}\,dt$ (12.7) (step 5) · **N16 [B]** "
                   r"$\overline{\partial u^m/\partial x_j}=\partial\overline{u^m}/\partial x_j$ (12.8) (step 4) · **N17 [B]** "
                   r"$\overline{\int u^m\,d\mathbf x}=\int\overline{u^m}\,d\mathbf x$ (12.9) (step 5) · **N18 [B]** the product does not "
                   r"commute: $\overline{u^m}\neq\bar u^m$ and $\overline{uv}\neq\bar u\,\bar v$ in general (step 7)")
D("D01", ref="12.4")
code(r"""
res = TS.check_averaging_rules(ens_u, ens_v, t, A=2.0)        # apply each operation before and after averaging; return the largest difference
for name, err in res.items():                                 # one line per rule
    print(f"{name:20s} largest |average(op u) - op(average u)| = {err:.1e}")
mean_prod, prod_means, cov = TS.product_average_split(ens_u, ens_v)   # mean(u v), mean(u) mean(v), mean(u' v') at every time
print("mean of product = product of means + covariance:", np.allclose(mean_prod, prod_means + cov))
print(f"time-mean covariance of the fluctuations: {cov.mean():.4f} (m/s)^2  (built in: 0.5 x 0.3^2 = {0.5*0.3**2:.4f})")
""", explain=r"""
1. `TS.check_averaging_rules` applies the **same discrete operator** (a sum, a factor 2, a finite-difference derivative, a
   trapezoid integral) to each member and then averages, and to the average directly; it returns the largest difference for each
   rule. Every entry is at round-off ($10^{-13}$ or smaller) — except `product`.
2. `TS.product_average_split` returns the three pieces of step 7 of the derivation. The leftover covariance is not noise: we
   built `ens_v` to share half of `ens_u`'s fluctuation, so it should be $0.5\times0.3^2=0.045$, and it is.""")
nb.worked_example("two realizations are enough to break the product rule", r"""
Run 1: $\tilde u=1$, $\tilde v=2$. Run 2: $\tilde u=3$, $\tilde v=6$.

1. Means: $\bar u=2$, $\bar v=4$, so $\bar u\,\bar v=8$.
2. Mean of the product: $(1\cdot2+3\cdot6)/2=10$.
3. Fluctuations: $u=-1,+1$; $v=-2,+2$; mean of $uv$ = $(2+2)/2=2$.
4. Check $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$: $10=8+2$ ✓.

The 2 that is left over is the kind of term the Reynolds stress is made of.""")
scratch(r"""
# From scratch: mean and variance by an explicit loop over members, and the product split by hand
sum_u = np.zeros(t.size); sum_u2 = np.zeros(t.size)           # running sums of u and u^2 at every time
for member in ens_u:                                          # one realization at a time
    sum_u += member                                           # add this run's value
    sum_u2 += member**2                                       # add this run's square
mean_mine = sum_u/len(ens_u)                                  # Eq. (12.10): divide the sum by N
var_mine = sum_u2/len(ens_u) - mean_mine**2                   # variance = mean of the square minus square of the mean
up, vp = ens_u - ens_u.mean(axis=0), ens_v - ens_v.mean(axis=0)   # fluctuations: subtract the ensemble mean at each time
cov_mine = (up*vp).mean(axis=0)                               # covariance of the fluctuations
assert np.allclose(mean_mine, TS.ensemble_average(ens_u))     # same numbers -> the library does exactly this
assert np.allclose(var_mine, TS.central_moment(ens_u, 2))     # second central moment = variance
assert np.allclose(cov_mine, TS.product_average_split(ens_u, ens_v)[2])   # the leftover of the product rule
tiny_u, tiny_v = np.array([[1.0], [3.0]]), np.array([[2.0], [6.0]])       # the two runs of the tiny example
print([float(a[0]) for a in TS.product_average_split(tiny_u, tiny_v)])    # [10.0, 8.0, 2.0]
""", r"""
1. The loop adds each run's value and square; dividing by the number of runs gives the mean and the mean square.
2. The three `assert` lines stop with an error if the hand-made numbers differ from `TS.ensemble_average`, `TS.central_moment`
   and `TS.product_average_split` — they do not.
3. The last line runs the tiny example through the library: 10, 8, 2.""")
note("N20 [B]", r"""
**How fast does an average converge?** The scatter of an $N$-member mean falls like $N^{-1/2}$: four times the runs for half
the error. Animation A1 adds members one at a time; the figure after it measures the slope.""")
nb.animation(r"""
ens_big = TS.make_ensemble(64, t, mean_fn, sigma=0.3, tau_c=0.5, seed=7)   # the decaying ensemble again, 64 members
stops = [1, 2, 4, 8, 16, 24, 32, 48, 64]                      # members shown per frame (9 frames)
tavg = TS.time_average(t, ens_big[0], window=2.0)             # a 2 s sliding time average of member 1 alone, Eq. (12.2)
sk4 = slice(None, None, 8)                                    # draw every 8th sample (0.08 s apart): the same picture in a much lighter image
fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.plot(t, mean_fn(t), color=C_REF, lw=3, alpha=0.5, label="true mean $e^{-t/10}$")   # the ghost we are trying to recover
(ln_new,) = ax.plot(t[sk4], ens_big[0][sk4], color="0.75", lw=0.5, label="newest member")      # the member just added (grey)
(ln_mean,) = ax.plot(t[sk4], ens_big[0][sk4], color=C_MEAN, lw=1.4, label="ensemble mean of N members")   # the running ensemble mean
ax.plot(t, tavg, color=C_FLUC, lw=1.4, label="2 s time average of member 1")          # what one run alone can give
ax.set(xlabel="time t [s]", ylabel="u [m/s]", xlim=(0, 20), ylim=(-0.6, 1.7))
ax.legend(fontsize=7, loc="upper right")
ttl = ax.set_title("")

def update(i):                                                # frame i: show the mean of the first N = stops[i] members
    N = stops[i]                                              # number of members in the average
    ln_new.set_ydata(ens_big[N - 1][sk4])                     # the newest member
    ln_mean.set_ydata(ens_big[:N].mean(axis=0)[sk4])          # Eq. (12.10) with N members
    ttl.set_text(f"N = {N}: scatter of the mean ≈ σ/√N = {0.3/np.sqrt(N):.3f} m/s")
    return ln_new, ln_mean, ttl

show_animation(animate(update, frames=len(stops), fig=fig, interval=400), player="frames", dpi=45)   # step with the buttons
""", explain=r"""
1. `stops` lists how many members are in the average at each frame; `update(i)` redraws the newest member (grey) and the mean of
   the first $N$ (purple).
2. The teal curve never changes: it is a 2-second sliding average of member 1 only — the best a single record can do.
3. `player="frames"` gives ◀ ▶ buttons: stop at $N=2,4,8,64$ and compare the purple wiggle with $\sigma/\sqrt N$ in the title.""")
nb.figure_notes(
    see="A thick grey ghost (the true decaying mean), a thin grey member, a purple curve that starts as wiggly as one member and "
        "calms down as $N$ grows, and a teal curve (time average of one member) that is smooth but wanders about the ghost.",
    read="The purple curve's scatter is $\\sigma/\\sqrt N$: 0.21 m/s at $N=2$, 0.0375 m/s at $N=64$. The teal curve averages "
         "about $2\\text{ s}/(2\\times0.5\\text{ s})=2$ independent samples, so it is no better than $N=2$ — a single record needs a long "
         "window, and a long window blurs a changing mean.",
    change="…the members were correlated with each other (say, all started from the same gust): the scatter would fall more "
           "slowly than $N^{-1/2}$, because the runs would not be independent samples.")
fig(r"""
ens_256 = TS.make_ensemble(256, t, lambda t: 0*t, sigma=0.3, tau_c=0.5, seed=8)   # 256 members around a zero mean
Ns = np.array([2, 4, 8, 16, 32, 64, 128, 256])                # ensemble sizes to test
scatter = np.array([np.sqrt(np.mean(ens_256[:N].mean(axis=0)**2)) for N in Ns])   # measured rms error of the N-member mean
est = np.array([np.mean(TS.standard_error(ens_256[:N])) for N in Ns])             # TS.standard_error: sample std / sqrt(N)
slope = np.polyfit(np.log(Ns), np.log(scatter), 1)[0]         # slope of the log-log line (P13)
fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.loglog(Ns, scatter, "o", color=C_MEAN, label=f"measured scatter (slope {slope:.2f})")   # what the ensemble does
ax.loglog(Ns, est, "s", color=C_FLUC, ms=4, label="TS.standard_error (from the sample)")     # what one can estimate from the data
ax.loglog(Ns, TS.standard_error_of_mean(0.3, Ns), "--", color=C_REF, label=r"$\sigma/\sqrt{N}$")   # the −1/2 law
ax.set(xlabel="number of members N [–]", ylabel="scatter of the mean [m/s]", title="An average converges like $N^{-1/2}$")
ax.legend(fontsize=8)
print(f"fitted slope {slope:.3f} (expect -0.50 ± 0.05); sigma/sqrt(64) = {TS.standard_error_of_mean(0.3, 64):.4f} m/s")
assert abs(slope + 0.5) < 0.05                                # the −1/2 law holds for our seeded ensemble
""",
    see="Three sets of points falling along one straight line on logarithmic axes: the measured scatter of the $N$-member mean "
        "(purple), the estimate `TS.standard_error` computes from the sample itself (teal), and $\\sigma/\\sqrt N$ (dashed).",
    read="A slope of −½ on log–log axes is the $N^{-1/2}$ law. To gain one decimal digit in a mean you need 100 times the data — "
         "the reason turbulence statistics are expensive.",
    change="…we averaged over time instead: $N$ becomes the number of *independent* samples in the record, about "
           "$\\Delta t/(2\\Lambda_t)$, where $\\Lambda_t$ is the memory time defined in C02.")
note("N21 [B]", r"""
**Central moments** measure the scatter about the mean.""",
     equation=r"\overline{(u-\langle u\rangle)^m}\equiv\frac1N\sum_n\big(u(\mathbf x,t{:}n)-\langle u(\mathbf x,t)\rangle\big)^m", ref="12.11")
note("N22 [B]", r"""
**Variance, standard deviation, rms, skewness, kurtosis.** The second central moment is the *variance*
$\overline{(u-\bar u)^2}=\overline{u^2}-\bar u^2$; its square root is the *standard deviation* $\sigma$ (equal to the
root-mean-square, *rms*, when the mean is zero). The third and fourth moments measure lopsidedness and the weight of the tails.

| | third moment | fourth moment | Gaussian |
|---|---|---|---|
| the book (raw central moments) | *skewness* $\overline{(u-\bar u)^3}$, units $u^3$ | *kurtosis* $\overline{(u-\bar u)^4}$, units $u^4$ | $0$ and $3\sigma^4$ |
| the usual, `normalized=True` | divided by $\sigma^3$, a pure number | divided by $\sigma^4$, a pure number | $0$ and $3$ |""")
code(r"""
xg_ = np.random.default_rng(3).normal(0.0, 2.0, 100_000)      # 100 000 Gaussian samples with standard deviation 2
raw = TS.statistics(xg_, normalized=False)                    # the book's convention: raw central moments
nrm = TS.statistics(xg_, normalized=True)                     # the usual convention: divided by sigma^3 and sigma^4
print({k: round(float(v), 3) for k, v in raw.items()})        # kurtosis ≈ 3 sigma^4 = 48
print({k: round(float(v), 3) for k, v in nrm.items()})        # skewness ≈ 0, kurtosis ≈ 3
assert abs(nrm["skewness"]) < 0.04 and abs(nrm["kurtosis"] - 3) < 0.08    # Gaussian values within sampling error
assert np.isclose(raw["variance"], np.mean(xg_**2) - np.mean(xg_)**2)    # variance = mean square − squared mean
""", explain=r"""
1. `TS.statistics` returns a dictionary: mean, variance, skewness, kurtosis, std, rms.
2. With `normalized=False` the kurtosis is $3\sigma^4=3\times16=48$; with `normalized=True` it is 3. Same data, two conventions —
   always say which.""")
P("P285", "the sinc function and np.sinc", r"""
$\sin(x)/x$ equals 1 at $x=0$ and is zero at $x=\pi,2\pi,\dots$. numpy's `np.sinc(x)` is $\sin(\pi x)/(\pi x)$, so
$\sin(a)/a$ = `np.sinc(a/np.pi)`. It appears whenever a wave is averaged over a window.""", code=r"""
a = np.array([0.0, 0.4*np.pi, np.pi])             # three arguments: 0, 0.4 pi, pi
print(np.round(np.sinc(a/np.pi), 4))              # sin(a)/a = 1, 0.7568, 0
""")
note("N23 [B]", r"""
**What a time window does to a decaying mean plus a wave** (our numbers). For the signal
$u=Ae^{-t/\tau_d}+B\cos\omega t$ ($\tau_d$ a decay time [s], $\omega$ an angular frequency [rad/s]) the window average is

$$\overline{u(t)}=\Big[\frac{\sinh(\Delta t/2\tau_d)}{\Delta t/2\tau_d}\Big]Ae^{-t/\tau_d}+\Big[\frac{\sin(\omega\Delta t/2)}{\omega\Delta t/2}\Big]B\cos\omega t .$$

Each part keeps its form and is multiplied by a factor: $\sinh(x)/x$ (always $\ge1$; $\sinh x=(e^x-e^{-x})/2$, Ch. 7 P168) for the
mean, the sinc function for the wave. A good window makes the first factor ≈ 1 and the second ≈ 0: $1\ll\omega\Delta t\ll\omega\tau_d$.""")
code(r"""
ts, Dt, A_, B_, td, w_ = sp.symbols("t Delta_t A B tau_d omega", positive=True)   # symbols: time, window, amplitudes, decay time, frequency
s_ = sp.symbols("s")                                          # the integration variable (time inside the window)
u_sym = A_*sp.exp(-s_/td) + B_*sp.cos(w_*s_)                  # the signal: decaying mean + wave
avg = sp.integrate(u_sym, (s_, ts - Dt/2, ts + Dt/2))/Dt      # Eq. (12.2): the window average, done symbolically
target = sp.sinh(Dt/(2*td))/(Dt/(2*td))*A_*sp.exp(-ts/td) + sp.sin(w_*Dt/2)/(w_*Dt/2)*B_*sp.cos(w_*ts)   # the claimed result
print("sympy window integral − formula =", sp.simplify((avg - target).rewrite(sp.exp)))   # 0: the formula is right
val, f_mean, f_wave = ch12.time_average_exp_cos(5.0, 1.0, 0.5, 10.0, 2*np.pi, 0.4)   # t = 5 s, A = 1, B = 0.5, tau_d = 10 s, 1 s period, window 0.4 s
print(f"mean factor {f_mean:.7f}, wave factor {f_wave:.4f}, window average at t = 5 s: {val:.5f} m/s")
""", explain=r"""
1. sympy integrates the signal over the window and subtracts the formula: the difference simplifies to 0.
2. `ch12.time_average_exp_cos(t, A, B, tau, omega, window)` returns the average and the two factors: 1.0000667 (the mean is
   untouched), 0.7568 (three quarters of the wave leaks through a 0.4-period window), average 0.98498 m/s.""")
nb.worked_example("a window of 0.4 periods", r"""
$A=1$, $B=0.5$ m/s, $\tau_d=10$ s, period 1 s ($\omega=2\pi$ rad/s), $\Delta t=0.4$ s.

1. Mean factor: $x=\Delta t/2\tau_d=0.02$, $\sinh(x)/x\approx1+x^2/6=1.00007$ — the mean passes untouched.
2. Wave factor: $x=\omega\Delta t/2=0.4\pi=1.2566$, $\sin(x)/x=0.9511/1.2566=0.757$ — 76 % of the wave leaks into the "mean".
3. With $\Delta t=1$ s (one whole period) the wave factor is $\sin(\pi)/\pi=0$: the wave is removed exactly.
4. With $\Delta t=20$ s the mean factor is $\sinh(1)/1=1.175$: now the window is too long and distorts the decaying mean by 17.5 %.""")
nb.plotly(r"""
tt = np.linspace(0.0, 30.0, 160)                            # time axis for the slider figure [s]
sig_ = np.exp(-tt/10.0) + 0.5*np.cos(2*np.pi*tt)              # the signal: decaying mean (10 s) + a wave of period 1 s

def f1(window):                                               # the four curves for one window length
    sliding = TS.time_average(tt, sig_, window)               # Eq. (12.2) on the sampled signal (NaN near the ends)
    formula = ch12.time_average_exp_cos(tt, 1.0, 0.5, 10.0, 2*np.pi, window)[0]   # the closed form of N23
    return {"signal": (tt, sig_), "sliding average (TS.time_average)": (tt, sliding),
            "true mean exp(-t/10)": (tt, np.exp(-tt/10.0)), "formula (time_average_exp_cos)": (tt, formula)}

windows = np.round(np.geomspace(0.1, 20.0, 12 if FAST else 16), 2)   # window lengths, log-spaced from 0.1 to 20 s
figF1 = slider_figure(f1, "Δt", windows, unit="s", xlabel="time t [s]", ylabel="u [m/s]",
                      title="A good window: longer than the wave, shorter than the decay", yrange=(-0.6, 1.6))
figF1.show()
""", explain=r"""
1. `f1(window)` returns the four curves for one window; `slider_figure` computes them for every slider position now, so the
   slider also works on the web page.
2. The sliding average and the formula lie on top of each other wherever the window fits inside the record.""")
nb.figure_notes(
    see="Drag Δt. Below about 0.3 s the average still carries most of the wave. Near 1, 2, 3 s (whole periods) the wave vanishes "
        "exactly. Beyond about 10 s the average rides visibly *above* the true mean.",
    read="Two factors, two failures: the sinc factor says how much wave survives, the $\\sinh(x)/x$ factor how much a long window "
         "inflates a decaying mean. Between them is the usable range $1\\ll\\omega\\Delta t\\ll\\omega\\tau_d$.",
    change="…the wave were slower (period 5 s): the usable range between 'longer than the wave' and 'shorter than the decay' "
           "would nearly close — the two time scales must be well separated.")
explainer("reynolds_averaging_window", "How long must you average to get 'the mean'?",
          "dragging the window and the number of members shows the two estimates of the mean failing in different ways — scatter "
          "falling like N^(−1/2), and the two window factors — and the product bars show the mean of a product minus the product of "
          "the means refusing to vanish.",
          ["Set the window to exactly one wave period: the teal curve loses the wave completely.",
           "Halve the window: 76 % of the wave is back.",
           "Raise N from 8 to 64: the scatter of the purple curve falls by √8 ≈ 2.8.",
           "Open the product bars and set the correlation to zero: the leftover vanishes."])
whatif(r"""…we multiplied two signals *before* averaging — as the term $\tilde u_j\tilde u_i$ in the momentum equation does? The
leftover covariance appears, and C04 shows that it acts on the mean flow like a stress. First we need to measure how a signal is
correlated with itself (C02).""")


# =====================================================================================================================
# A.4 §12.4 Correlations and Spectra — C02, C03
# =====================================================================================================================
nb.section("12.4", "Correlations and Spectra", intro=r"""
**What is this section about?** Two ways to describe how a fluctuating signal is organised: how long it remembers itself (the
correlation) and which frequencies carry its variance (the spectrum). They are the same information in two forms.""")
core("C02", r"The autocorrelation $R_{11}(\tau)=\overline{u_1(t)u_1(t+\tau)}=R_{11}(-\tau)$ (12.17) and the integral time scale "
            r"$\Lambda_t\equiv\int_0^\infty r_{11}(\tau)\,d\tau$ (12.18)",
     "How long does a turbulent signal remember itself?")
problem(r"""
If the wind is strong now, it is probably still strong half a second from now — and tells you nothing about the wind ten
minutes from now. Somewhere in between the signal "forgets". That forgetting time decides how long a record you need (C01), how
wide an eddy is, and — in C16 — how smoke spreads.""")
idea(r"""
signal u(t)    :  ~~/\~~~\/~~/\~~~      slide a copy by a lag τ, multiply point by point, average
copy u(t + τ)  :    ~~/\~~~\/~~/\~~~
lag τ = 0       → every product is a square        → average = variance      → r = 1
lag τ large     → the signs are unrelated          → average = 0             → r = 0
area under the normalised curve r(τ)  =  the memory time Λ_t
""", r"""
Symbols: $\tau$ is the time lag [s]; $R_{11}(\tau)$ the autocorrelation [m²/s²]; $r_{11}=R_{11}(\tau)/R_{11}(0)$ the
autocorrelation *coefficient*, a pure number between −1 and 1; $\Lambda_t$ the integral time scale [s].""")
remind("C02")
note("N24 [B]", r"""
**The correlation tensor** multiplies two velocity fluctuations — any two components, at any two points and times — and
averages. ⚠️ $R_{ij}$ was the rotation tensor in Ch. 2–3; here it is a correlation. `TS.correlation(a, b)` is the mean of the
product over the member axis.""",
     equation=r"R_{ij}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_i(\mathbf x_1,t_1)u_j(\mathbf x_2,t_2)}", ref="12.12")
note("N26 [B]", r"""
**The autocorrelation** is the case $i=j=1$: one component with itself at two points or two times.""",
     equation=r"R_{11}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_1(\mathbf x_1,t_1)u_1(\mathbf x_2,t_2)}", ref="12.13")
P("P287", "covariance, the covariance matrix and its ellipse", r"""
The *covariance* of two zero-mean variables is the mean of their product. The 2 × 2 *covariance matrix* has the two variances
on its diagonal and the covariance off it; its eigenvectors (Ch. 2's principal axes) point along and across the scatter cloud,
so the cloud is an ellipse tilted by the covariance. `TS.correlated_pair(n, r)` draws such a pair with correlation coefficient
$r$; `np.cov` computes the matrix.""", code=r"""
a_, b_ = TS.correlated_pair(100_000, -0.6, seed=0)        # two unit-variance signals with correlation coefficient −0.6
print(np.round(np.cov(a_, b_), 3))                         # [[1, −0.6], [−0.6, 1]]: variances on the diagonal, covariance off it
print(np.round(np.linalg.eigh(np.cov(a_, b_))[0], 3))      # eigenvalues 0.4 and 1.6: the squared half-axes of the scatter ellipse
""")
note("N27 [B]", r"""
**The correlation coefficient** divides the covariance by the two standard deviations, so that it is a pure number (the
arguments are those of $R_{ij}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_i(\mathbf x_1,t_1)u_j(\mathbf x_2,t_2)}$ (12.12)).
`TS.correlation_coefficient` computes it; numpy's `np.corrcoef` is the same number.""",
     equation=r"r_{12}\equiv R_{12}/\sqrt{R_{11}R_{22}}=\overline{u_1u_2}\big/\big(\sqrt{\overline{u_1^2}}\sqrt{\overline{u_2^2}}\big)", ref="12.14")
note("N28 [B]", r"""
**The autocorrelation coefficient** is the same normalisation for one component; it equals 1 at zero separation.""",
     equation=r"r_{11}\equiv R_{11}(\mathbf x_1,t_1,\mathbf x_2,t_2)/\sqrt{R_{11}(1,1)R_{11}(2,2)}", ref="12.15")
note("N25 [B]", r"""
**What a correlation coefficient looks like.** Four scatter clouds, from uncorrelated to strongly anticorrelated.""")
fig(r"""
fig, axs = plt.subplots(1, 4, figsize=(9.6, 2.7), sharex=True, sharey=True)
for ax, r_set, name in zip(axs, (0.0, 0.3, 0.9, -0.9), ("uncorrelated", "weakly correlated", "strongly correlated", "anticorrelated")):
    a_, b_ = TS.correlated_pair(1500, r_set, seed=1)                     # a pair with the chosen coefficient
    r_meas = TS.correlation_coefficient(a_, b_)                          # Eq. (12.14) on the samples
    assert abs(r_meas - np.corrcoef(a_, b_)[0, 1]) < 0.01                # numpy's built-in gives the same number (it removes the tiny sample means first)
    ax.plot(a_, b_, ".", ms=2, color=C_FLUC)                             # one dot per sample
    ax.set(title=f"{name}\nr = {r_meas:+.2f}", xlim=(-4, 4), ylim=(-4, 4), aspect="equal", xlabel="$u_1$")
axs[0].set_ylabel("$u_2$")
""",
    see="Four clouds of dots: round ($r\\approx0$), slightly tilted ($r\\approx0.3$), a thin ellipse rising to the right "
        "($r\\approx0.9$), a thin ellipse falling to the right ($r\\approx-0.9$).",
    read="The tilt and thinness of the cloud *are* the correlation: $r=\\pm1$ would be a straight line, $r=0$ a round cloud. "
         "Keep the last panel in mind — in C04 the two axes become $u$ and $v$ in a shear flow.",
    change="…we doubled $u_2$: the cloud would stretch vertically and the covariance would double, but $r$ would not change — "
           "that is what dividing by the standard deviations buys.")
P("P288", "a quadratic that is never negative has discriminant ≤ 0", r"""
If $a\lambda^2+b\lambda+c\ge0$ for **every** real $\lambda$ (with $a>0$), the parabola never dips below the axis, so it has at
most one real root. By the quadratic formula (Ch. 6, P159) that means the *discriminant* $b^2-4ac$ is not positive:
$b^2-4ac\le0$.""", code=r"""
lam = np.linspace(-5, 5, 1001)                    # many values of lambda
for a_, b_, c_ in ((1, 2, 2), (1, 2, 0.5)):       # two parabolas a lam^2 + b lam + c
    print(f"discriminant {b_**2 - 4*a_*c_:+.1f}: minimum of the parabola = {np.min(a_*lam**2 + b_*lam + c_):+.2f}")   # −4 → stays positive; +2 → dips below zero
""")
NOTES_IN["D02"] = (r"**N29 [B]** the Schwartz inequality $\lvert\overline{uv}\rvert\le\sqrt{\overline{u^2}}\sqrt{\overline{v^2}}$ (12.16) "
                   r"(steps 4–5); ⚠️ a trap, not a slip: the page prints it without the absolute value — true as printed, but then it says nothing about a negative correlation")
D("D02", ref="12.16")
D("D03", ref="12.17")
P("P289", "correlation by FFT with zero padding", r"""
The direct sum over lags costs (number of samples) × (number of lags). The FFT gives all lags at once: transform the record,
multiply by its complex conjugate, transform back. The record must first be padded with zeros to at least twice its length,
otherwise the end wraps round onto the beginning. "Unbiased" divides lag $m$ by $N-m$, the number of products actually
available; far lags have few products and are noisy.""", code=r"""
x = np.random.default_rng(0).normal(size=16)      # 16 random numbers
direct = np.correlate(x, x, mode="full")[15:]     # direct sums for lags 0 … 15
F = np.fft.rfft(x, 32)                            # FFT of the record zero-padded to 32
via_fft = np.fft.irfft(F*np.conj(F))[:16]         # multiply by the conjugate, transform back, keep lags 0 … 15
print(np.allclose(direct, via_fft))               # True: the same sums
""")
note("N31 [B]", r"""
**The integral time scale** is the area under the autocorrelation coefficient: the width of the rectangle of height 1 with the
same area. It is the *memory* of the signal. `TS.integral_scale(lag, r)` integrates to the first zero of $r$ (a measured tail
only adds noise).""",
     equation=r"\Lambda_t\equiv\int_0^\infty r_{11}(\tau)\,d\tau=\frac1{R_{11}(0)}\int_0^\infty R_{11}(\tau)\,d\tau", ref="12.18")
note("N32 [B]", r"""
**Correlation time and independent samples.** The correlation time $t_c$ is the first zero of $r_{11}$. The book counts the
independent samples of a record of length $\Delta t$ as $N\approx\Delta t/t_c$; we use the statistician's
$N=\Delta t/2\Lambda_t$ of primer P281 (`TS.effective_samples(record length, 2Λ_t)`), because $\Lambda_t$ is an area and far less
noisy than a zero crossing. Both are order-of-magnitude counts — this is the $N$ of C01's $N^{-1/2}$.""")
nb.worked_example("memory of an exponential correlation", r"""
Take $r(\tau)=e^{-\tau/\tau_c}$ with $\tau_c=0.5$ s.

1. $\Lambda_t=\int_0^\infty e^{-\tau/0.5}d\tau=0.5$ s.
2. A 100 s record holds about $100/0.5=200$ memory times.
3. So the standard error of its time mean is about $\sigma/\sqrt{200/2}=0.1\sigma$ (two memory times per independent sample
   for this shape).
4. For a Gaussian shape $r=e^{-\tau^2/\tau_c^2}$: $\Lambda_t=\tfrac{\sqrt\pi}2\tau_c=0.443$ s.""")
code(r"""
dt = 0.01                                                     # sampling interval [s]
t_rec = np.arange(0.0, 500.0 if FAST else 2000.0, dt)         # one long record: 2000 s (500 s in FAST mode)
u_rec = TS.make_ensemble(1, t_rec, lambda t: 0*t, 1.0, 0.5, seed=3)[0]   # one OU run: std 1 m/s, memory 0.5 s
lag, R = TS.autocorrelation(u_rec, dt, max_lag=500)           # Eq. (12.17): R(tau) for lags 0 … 5 s (FFT, unbiased)
r = R/R[0]                                                    # the coefficient: divide by the variance R(0)
Lam = TS.integral_scale(lag, r)                               # Eq. (12.18): area under r up to its first zero [s]
t_c = TS.correlation_time(lag, r)                             # first zero crossing of r [s]
N_eff = TS.effective_samples(t_rec[-1], 2*Lam)                # independent samples in the record, N ≈ record length / (2 Λ_t) (primer P281)
print(f"variance R(0) = {R[0]:.3f} m²/s²;  Lambda_t = {Lam:.3f} s (exact 0.5);  t_c = {t_c:.2f} s;  N ≈ {N_eff:.0f}")
print(f"exact scales of the two shapes: {TS.correlation_spectrum_pair('exponential', 1.0, 0.5)['Lambda_t']:.3f} s, "
      f"{TS.correlation_spectrum_pair('gaussian', 1.0, 0.5)['Lambda_t']:.3f} s")
assert abs(Lam - 0.5) < 0.08                                  # the seeded estimate is within sampling error of the exact memory
""", explain=r"""
1. `TS.autocorrelation(u, dt, max_lag)` returns the lags [s] and $R_{11}(\tau)$ [m²/s²]; it uses the zero-padded FFT of the primer
   and the unbiased divisor.
2. `TS.integral_scale` gives $\Lambda_t\approx0.48$ s for this seeded record (exact value 0.5 s): a statistical estimate, not a
   round-off-exact number.
3. `TS.correlation_time` finds where $r$ first reaches zero — for an exponential it never truly does, so $t_c$ is set by noise in the
   tail; that is why $\Lambda_t$ (an area) is the more robust measure of memory, and why $N$ is counted with $2\Lambda_t$: about 2000
   independent samples in this record.
4. The last line prints the exact scales of the tiny example from `TS.correlation_spectrum_pair`: 0.5 s and 0.443 s.""")
scratch(r"""
# From scratch: the autocorrelation as a direct sum over lags, and the integral scale by the trapezoid rule
um = u_rec - u_rec.mean()                                     # remove the record mean
n_lags = 200                                                  # lags 0 … 1.99 s
R_loop = np.array([np.sum(um[:um.size - m]*um[m:])/(um.size - m) for m in range(n_lags)])   # mean of u(t) u(t + m dt), unbiased
assert np.allclose(R_loop, R[:n_lags])                        # the FFT route gives exactly these sums
r_loop = R_loop/R_loop[0]                                     # coefficient
first_zero = np.argmax(r < 0) if np.any(r < 0) else r.size    # index of the first negative value of r
Lam_mine = np.trapezoid(r[:first_zero], lag[:first_zero])     # area under r up to (just before) its first zero
print(f"Lambda_t by hand {Lam_mine:.4f} s, TS.integral_scale {Lam:.4f} s")
assert np.isclose(Lam_mine, Lam, rtol=0.02)                   # same area (the library interpolates the last partial step)
""", r"""
1. `R_loop[m]` is literally "shift by $m$ samples, multiply, average" — the idea box in code.
2. The first `assert` shows the FFT route of the primer returns the same numbers.
3. The trapezoid area up to the first zero reproduces `TS.integral_scale` (the library also adds the sliver up to the interpolated
   zero, hence the 2 % tolerance).""")
note("N34 [B]", r"""
**The picture of a memory.** The figure shows the measured coefficient, the exact exponential, and the rectangle of height 1
whose width is $\Lambda_t$.""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(lag, r, color=C_FLUC, label="measured $r_{11}(\\tau)$")                       # the estimate from one record
ax.plot(lag, np.exp(-lag/0.5), "--", color=C_REF, label="exact $e^{-\\tau/\\tau_c}$")    # what an infinite record would give
ax.fill_between([0, Lam], [1, 1], color=C_MEAN, alpha=0.18, label=f"equal-area rectangle, width $\\Lambda_t$ = {Lam:.2f} s")   # the memory
ax.plot([t_c], [0], "o", color=C_RS, label=f"first zero $t_c$ = {t_c:.2f} s")           # the correlation time
ax.set(xlabel="lag τ [s]", ylabel="autocorrelation coefficient r [–]", xlim=(0, 4), ylim=(-0.2, 1.05),
       title="The memory of a signal is the area under its autocorrelation")
ax.legend(fontsize=8)
""",
    see="A teal curve falling from 1 toward 0 along the dashed exponential, then wandering slightly about zero; a purple "
        "rectangle of height 1 and width ≈ 0.5 s; an orange dot where the curve first crosses zero.",
    read="The rectangle has the same area as the curve: its width is the memory $\\Lambda_t$. Samples closer than this are "
         "nearly copies of each other; samples several $\\Lambda_t$ apart are independent.",
    change="…the record were 20 s instead of 2000 s: the tail would wander far from zero and $\\Lambda_t$ would depend on where "
           "we stop integrating — the reason `TS.integral_scale` stops at the first zero.")
note("N30 [B]", r"""
**Cross-correlation finds a delay.** Correlating two *different* signals as a function of lag,
$R_{uv}(\tau)=\overline{u(t)v(t+\tau)}$, peaks at the lag that lines them up — the way two anemometers a known distance apart
measure how fast eddies travel. Unlike an autocorrelation it is **not even**: $R_{uv}(\tau)=R_{vu}(-\tau)$ (D03 step 5).
Animation A2 slides a delayed copy past the original.""")
nb.animation(r"""
n_show = 600                                                  # samples shown: 6 s of the record
shift0 = 80                                                   # the built-in delay t_o = 80 samples = 0.8 s
u_a = u_rec[1000:1000 + n_show]                               # the original signal u(t)
u_b = u_rec[1000 - shift0:1000 - shift0 + n_show]             # v(t) = u(t − 0.8 s): a delayed copy
lag_x, R_x = TS.cross_correlation(u_rec[shift0:20000], u_rec[:20000 - shift0], dt, max_lag=200)   # R_uv for lags −2 … 2 s (v delayed by 0.8 s)
r_x = R_x/u_rec[:20000].var()                                 # normalise by the variance
tt_ = np.arange(n_show)*dt                                    # time axis of the window [s]
nfr = 24 if FAST else 45                                      # frames
lags_fr = np.linspace(-1.0, 2.0, nfr)                         # the trial lag swept by the animation [s]
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.8, 4.6))
a1.plot(tt_, u_a, color=C_REF, lw=0.9, label="u(t)")          # the fixed signal
(ln_v,) = a1.plot(tt_, u_b, color=C_FLUC, lw=0.9, label="v(t + τ): the delayed copy, shifted back by τ")   # the sliding copy
a1.set(xlabel="time t [s]", ylabel="signal [m/s]", xlim=(0, 6), ylim=(-3.5, 3.5))
a1.legend(fontsize=7, loc="upper right")
a2.plot(lag_x, r_x, color="0.85", lw=3)                       # the whole cross-correlation as a ghost
(ln_r,) = a2.plot([], [], color=C_RS, lw=1.8)                 # the part traced so far
(dot,) = a2.plot([], [], "o", color=C_RS)                     # the current lag
a2.axvline(0.8, color=C_REF, ls=":")
a2.set(xlabel="lag τ [s]", ylabel="$R_{uv}(τ)/\\sigma^2$ [–]", xlim=(-1, 2), ylim=(-0.2, 1.1))

def update(i):                                                # frame i: shift the copy by lags_fr[i] and mark the correlation there
    m = int(round(lags_fr[i]/dt))                             # the lag in samples
    seg = u_rec[1000 - shift0 + m:1000 - shift0 + m + n_show] # v(t + τ) on the window
    ln_v.set_ydata(seg)
    keep = lag_x <= lags_fr[i]                                # lags already visited
    ln_r.set_data(lag_x[keep], r_x[keep])
    dot.set_data([lags_fr[i]], [np.interp(lags_fr[i], lag_x, r_x)])
    a1.set_title(f"τ = {lags_fr[i]:+.2f} s: mean of the product u·v = {np.mean(u_a*seg):+.2f} m²/s²", fontsize=10)
    return ln_v, ln_r, dot

show_animation(animate(update, frames=nfr, fig=fig, interval=150), player="video", dpi=60)   # a smooth video
""", explain=r"""
1. `u_b` is the same record read 0.8 s earlier — a delayed copy, like a second probe downstream.
2. `TS.cross_correlation(u, v, dt, max_lag)` returns lags of **both signs** and $R_{uv}(\tau)$; if $v(t)=u(t-t_0)$ the peak is at
   $\tau=+t_0$.
3. `update(i)` slides the teal copy by the trial lag, prints the mean of the product in the title and extends the orange curve.""")
nb.figure_notes(
    see="Top: a grey signal and a teal copy that slides sideways. Bottom: an orange curve being traced over a grey ghost; it "
        "climbs to a peak of 1 exactly when the two curves in the top panel coincide, at τ = 0.8 s (dotted line).",
    read="The cross-correlation is largest at the lag that undoes the delay. It is not symmetric about τ = 0: its mirror image is "
         "the correlation with the two signals swapped.",
    change="…the copy were also *distorted* (eddies change while they travel): the peak would still sit at the travel time but "
           "would be lower than 1 — the frozen-turbulence hypothesis of C03 is the statement that it stays near 1.")
NOTES_IN["D04"] = (r"**N33 [B]** the Taylor microscale $\lambda_t^2\equiv-2\big/\big[d^2r_{11}/d\tau^2\big]_{\tau=0}$ (12.19): the lag at "
                   r"which the parabola that hugs the top of $r_{11}$ reaches zero (steps 3–4)")
D("D04", ref="12.19")
code(r"""
gauss = TS.correlation_spectrum_pair("gaussian", 1.0, 0.5)    # the exact pair r = exp(−τ²/τ_c²) with τ_c = 0.5 s
u_s = TS.smooth_signal(2**16, 0.01, spectrum=gauss["S"], seed=4)   # a random record with that (smooth) spectrum
lag_s, R_s = TS.autocorrelation(u_s, 0.01, max_lag=300)       # its autocorrelation for lags 0 … 3 s
r_s = R_s/R_s[0]                                              # coefficient
lam_t = TS.taylor_microscale(lag_s, r_s)                      # Eq. (12.19): from the curvature of r at zero lag [s]
print(f"Taylor microscale {lam_t:.3f} s (exact: tau_c = {gauss['lambda_t']:.3f} s);  integral scale {TS.integral_scale(lag_s, r_s):.3f} s (exact {gauss['Lambda_t']:.3f} s)")
print("an OU signal has no microscale:", TS.correlation_spectrum_pair("exponential", 1.0, 0.5)["lambda_t"])   # nan
""", explain=r"""
1. `TS.smooth_signal(n, dt, spectrum)` builds a random record with a prescribed smooth spectrum (random phases again, P286); we
   give it the spectrum of a Gaussian correlation.
2. `TS.taylor_microscale` fits the parabola $1-\tau^2/\lambda_t^2$ to the first few lags: $\lambda_t\approx0.5$ s $=\tau_c$, as the
   exact pair says.
3. ⚠️ An Ornstein–Uhlenbeck signal has a *cusp* at zero lag, $r\approx1-\lvert\tau\rvert/\tau_c$: no parabola fits and its microscale
   does not exist (`nan`). Real turbulence is smooth at the smallest scales, so its microscale exists.""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.plot(lag_s, r_s, color=C_FLUC, label="measured r(τ) of the smooth signal")                    # the estimate
ax.plot(lag_s, 1 - lag_s**2/lam_t**2, color=C_VISC, label=f"osculating parabola $1-\\tau^2/\\lambda_t^2$, $\\lambda_t$ = {lam_t:.2f} s")   # hugs the top
ax.plot(lag_s, np.exp(-lag_s/0.5), ":", color=C_REF, label="OU signal: a cusp, no parabola")      # for contrast
ax.set(xlabel="lag τ [s]", ylabel="r [–]", xlim=(0, 1.5), ylim=(-0.1, 1.05), title="The Taylor microscale is the foot of the parabola at the top")
ax.legend(fontsize=8)
""",
    see="A teal bell-shaped correlation with a flat top; a rose parabola that follows it near τ = 0 and reaches zero at "
        "$\\lambda_t\\approx0.5$ s; a dotted exponential that starts with a sharp corner instead of a flat top.",
    read="$\\Lambda_t$ (an area) measures how long the memory lasts; $\\lambda_t$ (a curvature) measures how fast the fastest "
         "wiggles are. The cusp of the dotted curve means infinitely fast wiggles — no microscale.",
    change="…the signal contained faster wiggles (a wider spectrum): the top of $r$ would be more sharply curved and "
           "$\\lambda_t$ smaller, with $\\Lambda_t$ hardly changed.")
whatif(r"""…we asked which *frequencies* carry the variance instead of how long the memory is? That is the Fourier transform of
this same curve (C03).""")

# ---------------------------------------------------------------------------------------------------------------------
core("C03", r"The energy spectrum $S_e(\omega)\equiv\frac1{2\pi}\int_{-\infty}^{+\infty}R_{11}(\tau)e^{-i\omega\tau}d\tau$ (12.20) and the variance "
            r"it distributes, $\overline{u_1^2}=\int_{-\infty}^{+\infty}S_e(\omega)\,d\omega$ (12.22)",
     "Which frequencies carry the variance, and why is that the same information as the memory?")
problem(r"""
A graphic equaliser shows how loud each pitch is. A spectrum does the same for a wind record: how much of the variance comes
from slow gusts, how much from fast flutter. Ocean and atmosphere records are almost always shown this way, and the −5/3 law of
C08 is a statement about a spectrum.""")
idea("", r"""
| Memory | The signal | Its spectrum |
|---|---|---|
| long ($\Lambda_t$ large) | slow wiggles | squeezed near $\omega=0$, tall there |
| short ($\Lambda_t$ small) | fast wiggles | flat out to $\omega\approx1/\Lambda_t$, low |

Two facts to remember: **the area under $S_e$ is the variance, always; its height at zero frequency is variance × $\Lambda_t/\pi$.**
Symbols: $\omega$ angular frequency [rad/s]; $S_e(\omega)$ the energy spectrum [m²/s² per rad/s = m²/s], defined for both signs
of $\omega$ ("two-sided").""")
remind("C03")
P("P290", "Fourier-transform pair: where the 2π sits", r"""
A Fourier-transform pair (forward and back) needs one factor $1/2\pi$ in total. The book puts it in the **forward** transform
and uses angular frequency $\omega$ [rad/s]. For a real, even function the complex exponential
$e^{-i\omega\tau}=\cos\omega\tau-i\sin\omega\tau$ (Euler's formula, Ch. 1 P45) reduces to a cosine, because the odd sine part
integrates to zero: $\frac1{2\pi}\int_{-\infty}^{\infty}R\,e^{-i\omega\tau}d\tau=\frac1\pi\int_0^\infty R\cos\omega\tau\,d\tau$.""", code=r"""
tau_c_ = 0.5                                                  # memory time [s]
val = quad(lambda tau: np.exp(-tau/tau_c_)*np.cos(2.0*tau), 0, np.inf)[0]/np.pi   # (1/pi) ∫ R cos(ω τ) dτ at ω = 2 rad/s
print(round(val, 4), round(tau_c_/(np.pi*(1 + (2.0*tau_c_)**2)), 4))   # 0.0796 both ways: the closed form τ_c/[π(1 + ω²τ_c²)]
""")
NOTES_IN["D05"] = (r"**N35 [B]** the inverse transform $R_{11}(\tau)\equiv\int_{-\infty}^{+\infty}S_e(\omega)e^{+i\omega\tau}d\omega$ (12.21) "
                   r"(step 6) · **N36 [B]** the area is the variance, $\overline{u_1^2}\equiv\int_{-\infty}^{+\infty}S_e(\omega)\,d\omega$ (12.22) "
                   r"(step 7) · **N37 [B]** the height at zero is the memory, $S_e(0)=\overline{u_1^2}\Lambda_t/\pi$ (steps 8–9)")
D("D05", ref="12.20")
nb.worked_example("the pair for an exponential memory", r"""
$\sigma^2=1$ m²/s², $\tau_c=0.5$ s.

1. $S_e(\omega)=\sigma^2\tau_c/[\pi(1+\omega^2\tau_c^2)]$.
2. At $\omega=0$: $0.5/\pi=0.159$ m²/s — and $\overline{u^2}\Lambda_t/\pi=1\times0.5/\pi$ ✓.
3. At $\omega=2$ rad/s ($=1/\tau_c$): $0.5/(2\pi)=0.0796$, half the peak: the spectrum is flat to about $1/\tau_c$.
4. Area: $\int_{-\infty}^\infty S_e\,d\omega=\frac{\sigma^2}\pi[\arctan\omega\tau_c]_{-\infty}^{\infty}=\sigma^2$ ✓.""")
code(r"""
pair = TS.correlation_spectrum_pair("exponential", 1.0, 0.5)  # the exact pair: r(τ), R(τ), S(ω) and its scales
lag_e = np.linspace(0.0, 15.0, 3001)                          # lags out to 30 memory times [s]
omega = np.linspace(0.0, 60.0, 601)                           # angular frequencies [rad/s]
S_num = TS.spectrum_from_correlation(lag_e, pair["R"](lag_e), omega)   # Eq. (12.20): cosine transform of R, done numerically
var_num = TS.spectrum_variance(omega, S_num)                  # Eq. (12.22): area under the two-sided spectrum
Lam_S = TS.integral_scale_from_spectrum(S_num[0], pair["R"](0.0))      # Λ_t = π S(0)/variance
R_back = TS.correlation_from_spectrum(omega, S_num, lag_e[:200])       # Eq. (12.21): back to the correlation
print(f"S(0) = {S_num[0]:.4f} m²/s (exact {pair['S0']:.4f});  S(2) = {np.interp(2.0, omega, S_num):.4f}")
print(f"area = {var_num:.3f} m²/s² (the variance is 1; the 2 % missing lies beyond ω = 60 rad/s);  Lambda_t = {Lam_S:.3f} s")
print(f"round trip correlation → spectrum → correlation: largest error {np.abs(R_back - pair['R'](lag_e[:200])).max():.3f}")
""", explain=r"""
1. `TS.correlation_spectrum_pair(kind, sigma, tau_c)` returns the exact functions and scales of a correlation–spectrum pair (a
   dictionary: `r`, `R`, `S`, `Lambda_t`, `S0`, …).
2. `TS.spectrum_from_correlation` does the cosine integral of the primer at every $\omega$; $S(0)=0.159$ and $S(2)=0.080$ are the
   numbers of the tiny example.
3. `TS.spectrum_variance` integrates over both signs of $\omega$: 0.98 of the variance lies below 60 rad/s, the rest in the slowly
   decaying tail.
4. `TS.correlation_from_spectrum` goes back; the small error near zero lag is the missing tail again.""")
P("P291", "the discrete Fourier transform as a Riemann sum; Parseval", r"""
`np.fft.rfft(u)*dt` approximates $\int u\,e^{-i\omega t}dt$ at the frequencies $\omega_k=2\pi k/(N\,dt)$ (a Riemann sum: a sum of
samples times the step). In the book's normalisation the two-sided spectral density of a record of length $T$ is
$\lvert\hat u\rvert^2/(2\pi T)$. *Parseval's theorem*: its sum times $\Delta\omega$, over both signs of $\omega$, is exactly the
variance. A *one-sided* spectrum (positive frequencies only) is twice the two-sided one.""", code=r"""
x = np.random.default_rng(1).normal(size=1024); x -= x.mean()   # a zero-mean record of 1024 samples
dts = 0.01; T_ = x.size*dts                                     # sampling step [s] and record length [s]
uhat = np.fft.rfft(x)*dts                                       # the Riemann sum for the Fourier integral
S2 = np.abs(uhat)**2/(2*np.pi*T_)                               # two-sided density |û|²/(2πT)
dw = 2*np.pi/T_                                                 # frequency step [rad/s]
total = dw*(S2[0] + 2*np.sum(S2[1:-1]) + S2[-1])                # both signs of ω: every interior frequency counts twice
print(round(total/x.var(), 6))                                  # 1.0: Parseval
""")
P("P292", "segment averaging (Welch) and leakage", r"""
One raw periodogram (the $\lvert\hat u\rvert^2$ of the previous primer) is as noisy as its own mean, however long the record:
more data gives more frequencies, not better ones. Cutting the record into $K$ segments and averaging their periodograms reduces
the scatter by $\sqrt K$ at the price of frequency resolution. `scipy.signal.welch` does this with overlapping, tapered segments
(the taper reduces *leakage*, the smearing of power between neighbouring frequencies). It returns a one-sided density per hertz:
divide by $4\pi$ for the book's two-sided density per rad/s.""", code=r"""
x = TS.make_ensemble(1, np.arange(0, 1000, 0.01), lambda t: 0*t, 1.0, 0.5, seed=2)[0]   # a 1000 s OU record
f_hz, P_hz = signal.welch(x, fs=100.0, nperseg=4096)          # Welch estimate: frequency [Hz], one-sided density [m²/s² per Hz]
print(round(float(P_hz[1:4].mean()/(4*np.pi)), 3), "≈ S_e near ω = 0 (exact 0.159; about 48 segments, so ±15 % scatter)")   # convert to two-sided per rad/s
""")
note("N38 [B]", r"""
**The same spectrum straight from a record** (the Wiener–Khinchin theorem, stated): the spectrum is also the squared magnitude
of the Fourier transform of the record itself,
$S_e(\omega)=\lim_{T\to\infty}\frac1{2\pi T}\big\lvert\int_{-T/2}^{+T/2}u(t)e^{-i\omega t}dt\big\rvert^2$. The proof expands the
square as a double integral over two times; changing variables to their difference $\tau$ leaves a triangular weight
$(1-\lvert\tau\rvert/T)$ that tends to 1 (the same triangle returns in C16). The figure computes one spectrum three ways.""")
fig(r"""
om_fit = np.geomspace(0.05, 200.0, 120)                               # frequencies for route 1 [rad/s]
lag_l, R_l = TS.autocorrelation(u_rec, dt, max_lag=3000)              # the measured correlation out to 30 s
taper = 0.5*(1 + np.cos(np.pi*lag_l/lag_l[-1]))                       # a smooth taper to zero at the largest lag (reduces leakage)
S_route1 = TS.spectrum_from_correlation(lag_l, R_l*taper, om_fit)     # route 1: transform of the measured correlation, Eq. (12.20)
om2, S_route2 = TS.periodogram(u_rec, dt, segments=32)                # route 2: averaged periodogram (Exercise 12.8 form)
f_hz, P_hz = signal.welch(u_rec, fs=1/dt, nperseg=8192)               # route 3: scipy's Welch estimate
fig, ax = plt.subplots(figsize=(6.6, 3.9))
ax.loglog(om2[1:], S_route2[1:], ".", ms=2, color="0.7", label="2: TS.periodogram, 32 segments")
ax.loglog(2*np.pi*f_hz[1:], P_hz[1:]/(4*np.pi), color=C_RS, lw=0.9, label="3: scipy.signal.welch (÷ 4π)")
ax.loglog(om_fit, S_route1, color=C_FLUC, lw=2.2, label="1: transform of the measured correlation")
ax.loglog(om_fit, pair["S"](om_fit), "--", color=C_REF, label="exact $\\sigma^2\\tau_c/[\\pi(1+\\omega^2\\tau_c^2)]$")
ax.set(xlabel="angular frequency ω [rad/s]", ylabel="$S_e$ [m²/s]", xlim=(0.05, 200), ylim=(1e-6, 1),
       title="One spectrum, three routes")
ax.legend(fontsize=8)
""",
    see="Grey dots, an orange line and a bold teal line all following the same dashed curve on logarithmic axes: flat at "
        "about 0.16 m²/s up to ω ≈ 2 rad/s, then falling with slope −2.",
    read="The corner at $\\omega\\approx1/\\tau_c=2$ rad/s is the memory seen from the frequency side; the flat level is "
         "$\\sigma^2\\Lambda_t/\\pi$. All three routes are estimates of the same function — they differ only in scatter.",
    change="…we used one segment instead of 32: the grey dots would scatter by 100 % of their value at every frequency, however "
           "long the record (primer P292).")
scratch(r"""
# From scratch: the cosine transform by the trapezoid rule, and a periodogram scaled by hand
om50 = np.linspace(0.0, 20.0, 50)                             # 50 frequencies [rad/s]
S_mine = np.array([np.trapezoid(pair["R"](lag_e)*np.cos(w*lag_e), lag_e)/np.pi for w in om50])   # (1/π) ∫ R cos(ωτ) dτ
S_lib = TS.spectrum_from_correlation(lag_e, pair["R"](lag_e), om50, rule="trapezoid")            # the library with the same rule
assert np.allclose(S_mine, S_lib, rtol=1e-10)                 # identical sums
x = u_rec[:2**14] - u_rec[:2**14].mean()                      # a zero-mean stretch of the record
T_ = x.size*dt                                                # its length [s]
P_mine = np.abs(np.fft.rfft(x)*dt)**2/(2*np.pi*T_)            # two-sided density |û|²/(2πT) (primer P291)
om_p, P_lib = TS.periodogram(x, dt)                           # the library: one boxcar segment
assert np.allclose(P_mine, P_lib)                             # same numbers
dw = om_p[1] - om_p[0]                                        # frequency step
assert np.isclose(dw*(P_mine[0] + 2*np.sum(P_mine[1:-1]) + P_mine[-1]), x.var(), rtol=1e-12)   # Parseval: the area is the variance
print("cosine transform and periodogram reproduced by hand; Parseval holds to round-off")
""", r"""
1. The cosine integral with `np.trapezoid` equals `TS.spectrum_from_correlation(..., rule="trapezoid")` to round-off.
2. The hand-scaled FFT equals `TS.periodogram`; summing it over both signs of $\omega$ returns the variance of the samples exactly
   (Parseval).""")
note("N39 [B]", r"""
**In space instead of time.** For a homogeneous field the correlation depends only on the separation vector $\mathbf r$.
`TS.spatial_correlation(a, b, dx)` computes it for a periodic field by averaging over the whole box.""",
     equation=r"R_{ij}(\mathbf r)\equiv\overline{u_i(\mathbf x)u_j(\mathbf x+\mathbf r)}", ref="12.23")
code(r"""
r_sp, R_sp = TS.spatial_correlation(us, us, dxf, axis=-1)     # R_11(r) of the §12.1 field, separation along x
print(f"R_11(0) = {R_sp[0]:.3f} = mean(u²) = {np.mean(us**2):.3f} m²/s²;  integral length scale {TS.integral_scale(r_sp, R_sp/R_sp[0]):.4f} m")
""", explain=r"""
`TS.spatial_correlation` returns separations $r$ [m] and $R_{11}(r)$; at $r=0$ it is the mean square, and the area under the
coefficient is an integral *length* scale — the size of the eddies of the synthetic field (centimetres to a decimetre in a 1 m box).""")
P("P293", "change of variable in a density; units of a spectral density", r"""
A density is "variance per unit of the axis". Changing the axis from $\omega$ to $k=\omega/U_0$ must keep the area:
$S(k)\,dk=S(\omega)\,d\omega$, so $S(k)=U_0\,S(\omega)$. Units: $S_e$ is [m²/s² per rad/s] = m²/s; a wavenumber spectrum
$S_{11}$ is [m²/s² per rad/m] = m³/s².""", code=r"""
U0 = 10.0                                                     # sweeping speed [m/s]
k_ax = omega/U0                                               # wavenumber axis k = ω/U0 [rad/m]
S_k = U0*pair["S"](omega)                                     # the density per unit k
print(round(np.trapezoid(pair["S"](omega), omega), 4), round(np.trapezoid(S_k, k_ax), 4))   # the two areas agree
""")
note("N40 [B]", r"""
**Taylor's frozen-turbulence hypothesis.** A probe moving at speed $U_0$ through the fluid — or a mean wind sweeping the eddies
past a mast — sees a pattern that hardly changes while it passes. Time then stands for distance: $x=U_0t$,
$\partial u_1/\partial x_1\approx-(1/U_0)\,\partial u_1/\partial t$, so $k_1=\omega/U_0$ and $S_{11}(k_1)=U_0S_e(\omega)$. It is
accurate when the fluctuations are small against the sweeping speed, $u_{rms}/U_0\ll1$. **Climate hook:** tower and aircraft
spectra are all converted to wavenumber this way.""")
code(r"""
k_conv, S_conv = TS.frequency_to_wavenumber_spectrum(omega, pair["S"](omega), 10.0)   # k_1 = ω/U0, S_11 = U0 S_e
x_fr, u_fr = TS.taylor_frozen(t_rec[:5], u_rec[:5], 10.0)     # the first five samples placed in space: x = U0 t
print("x = U0 t [m]:", np.round(x_fr, 2), "| area kept:", round(np.trapezoid(S_conv, k_conv), 4))
pattern = us[0]                                               # one row of the §12.1 field: a periodic pattern F(x)
for ratio in (0.05, 0.2, 0.5):                                # fluctuation level u_rms/U0
    errs = [ch12.frozen_field_probe(pattern, 10.0, ratio*10.0, dt=dxf/10.0, dx=dxf, n_samples=n_s, seed=0)["reconstruction_error"]
            for n_s in (nfield//8, nfield)]                   # a short record (1/8 of the 1 m line, about one eddy) and the whole line
    print(f"u_rms/U0 = {ratio:4.2f}: reconstruction error {errs[0]:.2f} over 0.125 m, {errs[1]:.2f} over 1 m (in units of the pattern's rms)")
""", explain=r"""
1. `TS.frequency_to_wavenumber_spectrum` rescales both axes; the area (the variance) is unchanged.
2. `ch12.frozen_field_probe` sweeps a known spatial pattern past a probe at $U_0+u'$, where $u'$ is a random large-eddy velocity,
   and compares what the frozen hypothesis reconstructs with the true pattern (a simple "random sweeping" model of ours).
3. The error grows with $u_{rms}/U_0$ **and with the length of the record**: a wrong sweeping speed displaces the pattern by
   $(u'/U_0)\times$ distance, so a point-by-point reconstruction fails once that displacement exceeds an eddy size — at
   $u_{rms}/U_0=5$ % the error is 0.07 over one eddy (0.125 m) and 0.35 over the whole 1 m line. A *spectrum* is more forgiving: it is only stretched
   along $k_1$ by the factor $1+u'/U_0$. The hypothesis is for weak turbulence on a strong mean flow.""")
nb.plotly(r"""
om_ax = np.geomspace(0.02, 100.0, 150)                       # frequency axis [rad/s]

def f2(tau_c):                                                # the spectrum and its two landmarks for one memory time
    p = TS.correlation_spectrum_pair("exponential", 1.0, tau_c)   # exact pair with variance 1 m²/s²
    return {"S_e(ω) for r = exp(−τ/τ_c)": (om_ax, p["S"](om_ax)),
            "level at zero frequency: σ²Λ_t/π": ([om_ax[0], 1.0/tau_c], [p["S0"], p["S0"]]),
            "corner frequency ω = 1/τ_c": ([1.0/tau_c, 1.0/tau_c], [1e-5, p["S0"]])}

figF2 = slider_figure(f2, "τ_c", np.round(np.geomspace(0.1, 5.0, 11 if FAST else 17), 3), unit="s", xlabel="angular frequency ω [rad/s]",
                      ylabel="S_e [m²/s]", title="Longer memory: a narrower, taller spectrum")
figF2.update_xaxes(type="log", range=[np.log10(0.02), 2.0])   # logarithmic frequency axis
figF2.update_yaxes(type="log", range=[-5.0, 0.5])             # logarithmic spectrum axis, fixed so the change is visible
figF2.show()
""", explain=r"""
1. `f2(tau_c)` returns the exact spectrum of an exponential memory and two landmarks: the horizontal level $\sigma^2\Lambda_t/\pi$
   and the vertical line at the corner frequency $1/\tau_c$.
2. The correlation that belongs to each slider position is the curve of C02's figure, $r=e^{-\tau/\tau_c}$, with a rectangle of
   width $\tau_c$.""")
nb.plotly(r"""
tau_lin = np.linspace(0.0, 15.0, 150)                        # lags [s], linear axis

def f2b(tau_c):                                               # the correlation that belongs to the same slider value
    p = TS.correlation_spectrum_pair("exponential", 1.0, tau_c)
    return {"r(τ) = exp(−τ/τ_c)": (tau_lin, p["r"](tau_lin)),
            "equal-area rectangle, width Λ_t = τ_c": ([0.0, tau_c, tau_c], [1.0, 1.0, 0.0])}

figF2b = slider_figure(f2b, "τ_c", np.round(np.geomspace(0.1, 5.0, 11 if FAST else 17), 3), unit="s", xlabel="lag τ [s]", ylabel="r [–]",
                       title="…while the correlation and its rectangle widen", xrange=(0, 15), yrange=(0, 1.05))
figF2b.show()
om_area = np.geomspace(1e-4, 1e5, 4000)                       # a very wide frequency range for the area check
for tc_ in (0.1, 1.0, 5.0):                                   # the area under the spectrum at three memory times
    p_ = TS.correlation_spectrum_pair("exponential", 1.0, tc_)
    print(f"τ_c = {tc_:3.1f} s: S_e(0) = {p_['S0']:.4f} m²/s, corner at {1/tc_:5.2f} rad/s, area 2∫S dω = {2*np.trapezoid(p_['S'](om_area)*om_area, np.log(om_area)):.3f} m²/s²")
""", explain=r"""
1. The second slider figure shows the other half of the pair on **linear** axes: the correlation $r(\tau)$ and the rectangle of
   width $\Lambda_t=\tau_c$ with the same area. Use the same slider value in both figures.
2. The printed lines put numbers on what a log–log plot cannot show: the level $S_e(0)$ grows in proportion to $\tau_c$, the
   corner moves as $1/\tau_c$, and the area under the spectrum (both signs of $\omega$) stays at the variance, 1 m²/s².""")
nb.figure_notes(
    see="First figure (log–log): drag $\\tau_c$ up and the spectrum's corner moves to the left while its flat level rises in "
        "proportion. Second figure (linear): at the same slider value the correlation curve and its rectangle widen to the right.",
    read="Width in lag × width in frequency stays constant — corner frequency × $\\Lambda_t$ = 1 at every slider value: the corner sits at $\\omega\\approx1/\\tau_c$ and the level at "
         "$\\sigma^2\\Lambda_t/\\pi$. The area under the spectrum is the variance for every slider value — the cell prints 1.000 "
         "three times; on logarithmic axes you cannot see that by eye.",
    change="…the correlation oscillated (a damped cosine): a peak would leave zero frequency and $\\Lambda_t$ would collapse, "
           "because positive and negative lobes cancel in the area. The explainer below has that preset.")
explainer("correlation_and_spectrum", "Why are a correlation and a spectrum the same information?",
          "stretching the memory narrows the spectrum in front of you while its area stays the variance and its height at zero "
          "follows the integral scale; sliding the lag shows the product being averaged.",
          ["Double the memory time: S(0) doubles, the corner frequency halves, the shaded area does not move.",
           "Switch to 'damped cosine' and raise ω₀: a peak leaves zero frequency and Λ_t collapses — an oscillating correlation has little net area.",
           "Click the spectrum: the inspector shows the cosine-transform arithmetic at that ω.",
           "Switch the axis to wavenumber with U₀ = 10 m/s."])
whatif(r"""…we now put the decomposition mean + fluctuation into the equations of motion and averaged? Every term passes through
the average except the product — C04.""")


# =====================================================================================================================
# A.5 §12.5 Averaged Equations of Motion — R01 R02, C04
# =====================================================================================================================
nb.section("12.5", "Averaged Equations of Motion", intro=r"""
**What is this section about?** We average the equations of motion. The mean flow obeys almost the same equations as before —
plus one new term made of fluctuations, which acts like a stress and for which the averaging gives no equation.""")
nb.recap("R01", "The Boussinesq equations in flux form", r"""
From Ch. 4: continuity $\partial\tilde u_i/\partial x_i=0$ (4.10), momentum
$\frac{\partial\tilde u_i}{\partial t}+\tilde u_j\frac{\partial\tilde u_i}{\partial x_j}=-\frac1{\rho_0}\frac{\partial\tilde p}{\partial x_i}-g[1-\alpha(\tilde T-T_0)]\delta_{i3}+\nu\frac{\partial^2\tilde u_i}{\partial x_j^2}$ (4.86)
and heat $\frac{\partial\tilde T}{\partial t}+\tilde u_j\frac{\partial\tilde T}{\partial x_j}=\kappa\frac{\partial^2\tilde T}{\partial x_j^2}$ (4.89).
Here $\rho_0$ is a constant reference density, $\alpha$ the thermal expansion coefficient [1/K], $T_0$ a reference temperature,
$\delta_{i3}$ picks the vertical component, and repeated indices are summed (Ch. 2). Adding
$\tilde u_i\,\partial\tilde u_j/\partial x_j=0$ turns $\tilde u_j\,\partial\tilde u_i/\partial x_j$ into
$\partial(\tilde u_j\tilde u_i)/\partial x_j$ — the *flux form*:
$\frac{\partial\tilde u_i}{\partial t}+\frac{\partial}{\partial x_j}(\tilde u_j\tilde u_i)=-\frac1{\rho_0}\frac{\partial\tilde p}{\partial x_i}-g[1-\alpha(\tilde T-T_0)]\delta_{i3}+\nu\frac{\partial^2\tilde u_i}{\partial x_j^2}$.
This is step 1 of the derivation D06 below. (Ch. 4 also writes the momentum equation with deviations from a hydrostatic state,
$D\mathbf u/Dt=-\rho_0^{-1}\nabla p'+(\rho'/\rho_0)\mathbf g+\nu\nabla^2\mathbf u$: put $\tilde p=p_0(z)+p'$ with $dp_0/dz=-\rho_0g$ and
$\rho'=-\rho_0\alpha(\tilde T-T_0)$ into the form above and the two are the same equation; D10 quotes that second form.)""", where="Ch. 4 §4.9")
core("C04", r"The Reynolds-averaged momentum equation $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\frac1{\rho_0}\frac{\partial\bar\tau_{ij}}{\partial x_j}$, "
            r"$\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$ (12.30)",
     "How can fluctuations that average to zero push on the mean flow?")
problem(r"""
An aircraft designer wants the mean drag, a climate modeller the mean wind; neither wants every gust. (**N41 [C]**: the flight
lasts hours, the gusts milliseconds, and the range of scales is far too wide to compute — C07 puts a number on it.) So we ask
what equation the *mean* obeys. The surprise: the gusts do not average out of it.""")
idea(r"""
fast  ────────────►   y + ℓ     a parcel coming DOWN  (v < 0) arrives too fast   (u > 0)   → uv < 0
slow  ─────►          y         a parcel going  UP    (v > 0) arrives too slow   (u < 0)   → uv < 0
every exchange carries x-momentum downward: mean(uv) < 0  ⇒  a stress −ρ₀·mean(uv) > 0 on the mean flow
""", r"""
Each of $u$ and $v$ averages to zero, but their *product* does not: in a shear, going up goes with being slow. That
correlation is a transport of momentum, and a transport of momentum is a stress.""")
remind("C04")
note("N42 [B]", r"""
**The Reynolds decomposition**: every field = mean + fluctuation.""",
     equation=r"\tilde u_i=U_i+u_i,\quad\tilde p=P+p,\quad\tilde\rho=\bar\rho+\rho',\quad\tilde T=\bar T+T'", ref="12.24")
note("N43 [B]", r"""
**The mean of a total field is the mean part** (by definition of the split).""",
     equation=r"\overline{\tilde u_i}=U_i,\quad\overline{\tilde p}=P,\quad\overline{\tilde\rho}=\bar\rho,\quad\overline{\tilde T}=\bar T", ref="12.25")
note("N44 [B]", r"""
**A fluctuation has zero mean** — subtract the two lines above.""",
     equation=r"\overline{u_i}=0,\quad\bar p=0,\quad\overline{\rho'}=0,\quad\overline{T'}=0", ref="12.26")
code(r"""
mean_u, fluct_u = TS.reynolds_decompose(ens_u)                # Eq. (12.24): mean (one value per time) and fluctuation (per member and time)
print(mean_u.shape, fluct_u.shape)                            # (2000,) and (64, 2000)
print(f"largest |mean of the fluctuation|: {np.abs(fluct_u.mean(axis=0)).max():.1e}")   # Eq. (12.26): zero to round-off
""", explain=r"""
`TS.reynolds_decompose(samples)` subtracts the ensemble mean from every member. The mean of what is left is zero to round-off —
$\overline{u_i}=0$ (12.26) is an identity, not an approximation.""")
P("P294", "a sympy averaging operator", r"""
We can teach sympy the rules of D01: the average is linear, leaves a mean alone, kills a single fluctuation, and keeps a
product of two fluctuations as a new symbol. `ch12.rans_sympy()` uses exactly such an operator on the full equations; here is
the idea on one product.""", code=r"""
U_, V_, u_, v_, uv_bar = sp.symbols("U V u v uv_bar")         # means U, V; fluctuations u, v; the symbol for mean(u v)
product = sp.expand((U_ + u_)*(V_ + v_))                      # (U + u)(V + v) = UV + Uv + uV + uv
averaged = product.subs(u_*v_, uv_bar).subs({u_: 0, v_: 0})   # mean(uv) → a new symbol; a lone fluctuation averages to 0
print(product, " → ", averaged)                               # U*V + U*v + V*u + u*v  →  U*V + uv_bar
""")
NOTES_IN["D06"] = (r"**N45 [B]** mean continuity $\partial U_i/\partial x_i=0$ (12.27) (step 3) · **N46 [B]** the fluctuation is "
                   r"divergence-free too, $\partial u_i/\partial x_i=0$ (12.28) (step 4) · **N47 [B]** the decomposed momentum equation "
                   r"$\frac{\partial(U_i+u_i)}{\partial t}+\frac{\partial}{\partial x_j}\big((U_j+u_j)(U_i+u_i)\big)=-\frac1{\rho_0}\frac{\partial(P+p)}{\partial x_i}-g[1-\alpha(\bar T+T'-T_0)]\delta_{i3}+\nu\frac{\partial^2(U_i+u_i)}{\partial x_j^2}$ (12.29) "
                   r"(step 5) · **N48 [B]** the term-by-term averages (steps 7–9)")
D("D06", ref="12.30")
code(r"""
eqs = ch12.rans_sympy()                                       # the averaging of D06 done by sympy (cached after the first run)
for label, _expr in eqs["steps"]:                             # the steps it took, by name
    print("•", label)
print("checks against the book's equations:", eqs["checks"])  # each True: the averaged equation equals the book's form
ny_ = 32                                                      # a small grid for an ensemble of fields
members = [TS.synthetic_solenoidal_field(ny_, 1.0, spectrum=lambda K: K**4*np.exp(-(K/12.0)**2), seed=s) for s in range(16)]   # 16 divergence-free fields
ens_fu = np.array([m_[0] for m_ in members]); ens_fv = np.array([m_[1] for m_ in members])   # members on axis 0
div_mean = ch12.mean_divergence(ens_fu, ens_fv, 1.0/ny_, 1.0/ny_)   # divergence of the ensemble-MEAN field
scale = np.sqrt(np.mean(ens_fu**2))*ny_                       # a typical (speed / grid spacing) to compare with
print(f"largest divergence of the mean field / (u_rms/dx): {np.abs(div_mean).max()/scale:.1e}")   # round-off: Eq. (12.27)
""", explain=r"""
1. `ch12.rans_sympy()` inserts mean + fluctuation into the Boussinesq equations and averages with the operator of primer P294. Its
   `checks` say that the results equal $\partial U_i/\partial x_i=0$ (12.27), $\partial u_i/\partial x_i=0$ (12.28), the mean momentum
   equation of this block's title, and the mean temperature and scalar equations stated below. The collected momentum equation has
   exactly one term with no counterpart in the instantaneous equation: $\partial\overline{u_iu_j}/\partial x_j$.
2. `ch12.mean_divergence` checks **N45** on numbers: every member of an ensemble of divergence-free fields is divergence-free, so
   the mean field is too (to round-off), and so is each fluctuation.""")
P("P295", "counting the independent components of a symmetric tensor", r"""
A symmetric 3 × 3 tensor has 6 independent components (3 on the diagonal + 3 above it). A fully symmetric triple product
$\overline{u_iu_ju_k}$ has 10. `itertools.product(range(3), repeat=3)` lists all 27 index triples; sorting each and collecting
them in a set counts the distinct ones.""", code=r"""
pairs = {tuple(sorted(c)) for c in itertools.product(range(3), repeat=2)}     # distinct index pairs (i, j) up to order
triples = {tuple(sorted(c)) for c in itertools.product(range(3), repeat=3)}   # distinct index triples up to order
print(len(pairs), len(triples))                               # 6 and 10
""")
note("N49 [B]", r"""
**The Reynolds stress tensor** $-\rho_0\overline{u_iu_j}$ [Pa] is symmetric, so it has six components. The diagonal ones are
normal stresses that add to the mean pressure; the off-diagonal ones are shear stresses. It is a covariance matrix (primer
P287) times $-\rho_0$. Except within a hair of a wall it is far larger than the viscous stress $2\mu\bar S_{ij}$, where
$\bar S_{ij}=\tfrac12(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ is the mean strain rate. `TS.velocity_covariance`
gives $\overline{u_iu_j}$ from samples, `TS.reynolds_stress` multiplies by $-\rho_0$, and `TS.anisotropy_tensor` subtracts
the isotropic part ($\overline{u_iu_j}/2\bar e-\delta_{ij}/3$; zero for turbulence with no preferred direction).""")
nb.recap("R02", "Sign convention of a stress on a face", r"""
As for the stress cube of Ch. 2: $\tau_{ij}$ is the force per area in direction $i$ on the face whose outward normal points
along $+j$. So a positive $-\rho_0\overline{uv}$ on the face whose outward normal is $+y$ pushes the fluid below it along $+x$ —
drawn in the left panel of the figure below.""", where="Ch. 2 §2.4")
nb.current_core = "C04"
note("N50 [B]", r"""
**The parcel argument.** A parcel displaced upward by $\ell$ keeps its old mean speed for a while, so where it arrives it is
too slow: $u\approx-\ell\,dU/dy$. Multiply by $v$ and average: $\overline{uv}\approx-\overline{v\ell}\,dU/dy<0$ when $dU/dy>0$,
because $v$ and $\ell$ have the same sign (a parcel that is going up has gone up).""")
note("N51 [B]", r"""
**Reynolds stress = momentum flux.** $\rho_0\overline{(U+u)v}=\rho_0U\bar v+\rho_0\overline{uv}=\rho_0\overline{uv}$ is the mean
flux of $x$-momentum in the $y$-direction carried by the fluctuations (the first term vanishes because $\bar v=0$). A downward
flux of $x$-momentum is what a shear stress *is*.

> ⚠️ **slip #16 — the book prints** the mean momentum flux as
> $\rho_0\overline{(U+u)v}=\rho_0U\bar u+\rho_0\overline{uv}=\rho_0\overline{uv}$ **; the correct form is**
> $\rho_0\overline{(U+u)v}=\rho_0U\bar v+\rho_0\overline{uv}=\rho_0\overline{uv}$: expand $(U+u)v=Uv+uv$ and average — $U$ is a
> mean, so $\overline{Uv}=U\bar v=0$. The middle term that vanishes is $\rho_0U\bar v$, not $\rho_0U\bar u$; the result is unchanged.""")
nb.worked_example("the stress carried by parcels", r"""
$dU/dy=2$ s⁻¹; parcels move $\ell_{rms}=0.1$ m with $v_{rms}=0.5$ m/s and keep all their momentum.

1. $u\approx-\ell\,dU/dy$ has rms $0.2$ m/s.
2. $\overline{uv}=-v_{rms}\ell_{rms}\,dU/dy=-0.5\times0.1\times2=-0.1$ m²/s².
3. In air ($\rho_0=1.2$ kg/m³): Reynolds shear stress $-\rho_0\overline{uv}=0.12$ Pa.
4. Viscous stress $\mu\,dU/dy=1.8\times10^{-5}\times2=3.6\times10^{-5}$ Pa — 3000 times smaller.
5. The viscosity that would do the same job: $\nu_T=-\overline{uv}/(dU/dy)=0.05$ m²/s, against $\nu=1.5\times10^{-5}$ m²/s.""")
code(r"""
out = ch12.displaced_parcel_uv(2.0, 0.1, 0.5, n=100_000, seed=0)   # 100 000 parcels: dU/dy = 2 1/s, l_rms = 0.1 m, v_rms = 0.5 m/s
expected = ch12.parcel_uv_expected(2.0, 0.1, 0.5)             # the tiny example's −v_rms l_rms dU/dy [m²/s²]
print(f"mean(uv) = {out['uv']:.4f} ± {out['stderr']:.4f} m²/s²  (expected {expected:.3f});  r_uv = {out['r_uv']:.2f};  nu_T = {out['nu_T']:.3f} m²/s")
assert abs(out["uv"] - expected) < 5*out["stderr"]            # within 5 standard errors (primer P281)
half = ch12.displaced_parcel_uv(2.0, 0.1, 0.5, n=100_000, seed=0, correlation=0.5)   # parcels that only half remember where they came from
print(f"with correlation 0.5 between v and l: mean(uv) = {half['uv']:.4f} m²/s², r_uv = {half['r_uv']:.2f}")
print(f"Reynolds stress in air: {-1.2*out['uv']:.3f} Pa;  viscous stress: {1.8e-5*2.0:.1e} Pa")
""", explain=r"""
1. `ch12.displaced_parcel_uv(dUdy, l_rms, v_rms, n, seed)` draws $n$ parcels: a displacement $\ell$, a vertical velocity $v$ with
   the same sign, and $u=-\ell\,dU/dy$. It returns a dictionary with the samples (`u`, `v`, `l`), the mean product `uv`, its standard
   error, the correlation coefficient `r_uv` and the equivalent viscosity `nu_T`.
2. The mean product is $-0.100$ m²/s² within its standard error — steps 2 and 5 of the tiny example.
3. With `correlation=0.5` (a parcel's $v$ is only half related to where it came from) the stress halves and $r_{uv}=-0.5$: the
   stress is as large as the *correlation* between going up and being slow.""")
fig(r"""
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(10.5, 3.6), gridspec_kw=dict(width_ratios=[1.1, 0.8, 1.2]))
parcel_sketch(a1, shear=1.0)                                  # drawing: two parcels swapping levels in a shear
stress_element_sketch(a2)                                     # drawing: the stress on the faces of a fluid element (recap R02)
a1.set_title("a parcel swaps levels in a shear", fontsize=9); a2.set_title("the stress on an element", fontsize=9)   # short panel titles
uu_, vv_ = half["u"][:4000], half["v"][:4000]                 # 4000 of the half-correlated parcels (a fatter cloud shows the tilt better)
a3.plot(uu_, vv_, ".", ms=1.5, color=C_FLUC)                  # one dot per parcel in the (u, v) plane
a3.fill_between([-0.8, 0], 0, 2.0, color=C_RS, alpha=0.12)    # quadrant 2: u < 0, v > 0 (going up, too slow)
a3.fill_between([0, 0.8], -2.0, 0, color=C_RS, alpha=0.12)    # quadrant 4: u > 0, v < 0 (coming down, too fast)
cov2 = np.cov(uu_, vv_)                                       # the 2 × 2 covariance matrix of the cloud
evals, evecs = np.linalg.eigh(cov2)                           # its principal axes (primer P287)
th = np.linspace(0, 2*np.pi, 100)                             # angle round the ellipse
ell = evecs @ (2*np.sqrt(evals)[:, None]*np.array([np.cos(th), np.sin(th)]))   # the 2-standard-deviation ellipse
a3.plot(ell[0], ell[1], color=C_RS, lw=1.8)                   # the covariance ellipse
a3.set(xlabel="u [m/s]", ylabel="v [m/s]", xlim=(-0.8, 0.8), ylim=(-2.0, 2.0), title=f"mean(uv) = {half['uv']:.3f} m²/s²")
""",
    see="Left: a shear profile with two parcels swapping levels. Middle: a fluid element with the shear stress drawn on its faces. "
        "Right: a cloud of (u, v) dots, tilted so that most of it lies in the two shaded quadrants (u < 0 with v > 0, u > 0 with "
        "v < 0), with its covariance ellipse.",
    read="The tilt of the cloud *is* the Reynolds stress: the mean of $uv$ is negative because the shaded quadrants hold more "
         "and larger products. Each of $u$ and $v$ alone still averages to zero.",
    change="…$dU/dy<0$: the cloud tilts the other way and $\\overline{uv}>0$. With no shear the cloud is round and there is no "
           "shear stress — the isotropic case of C05.")
nb.animation(r"""
nfr = 30 if FAST else 60                                      # frames
per = 25                                                      # parcels added per frame
pu, pv = half["u"][:nfr*per], half["v"][:nfr*per]             # the parcels in the order they are released
run = np.cumsum(pu*pv)/np.arange(1, pu.size + 1)              # running mean of u v after each parcel
se5 = 5*np.std(pu*pv)/np.sqrt(np.arange(1, pu.size + 1))      # ±5 standard errors of that running mean
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.3))
(dots,) = a1.plot([], [], ".", ms=3, color=C_FLUC)            # the (u, v) cloud so far
(new,) = a1.plot([], [], "o", ms=5, color=C_RS)               # the parcels of this frame
a1.axhline(0, color=C_REF, lw=0.6); a1.axvline(0, color=C_REF, lw=0.6)
a1.set(xlabel="u [m/s]", ylabel="v [m/s]", xlim=(-0.8, 0.8), ylim=(-2, 2))
n_ax = np.arange(1, pu.size + 1)                              # number of parcels
a2.fill_between(n_ax, half["expected"] - se5, half["expected"] + se5, color=C_RS, alpha=0.15, label="expected ± 5 standard errors")
a2.axhline(half["expected"], color=C_REF, ls="--")
(ln,) = a2.plot([], [], color=C_RS, label="running mean of uv")
a2.set(xlabel="parcels exchanged [–]", ylabel="mean(uv) [m²/s²]", xlim=(1, pu.size), ylim=(-0.15, 0.05), xscale="log")
a2.legend(fontsize=7, loc="upper right")

def update(i):                                                # frame i: release 25 more parcels
    k = (i + 1)*per                                           # parcels so far
    dots.set_data(pu[:k], pv[:k])
    new.set_data(pu[k - per:k], pv[k - per:k])
    ln.set_data(n_ax[:k], run[:k])
    a1.set_title(f"{k} parcels: mean(uv) = {run[k - 1]:+.3f} m²/s²", fontsize=10)
    return dots, new, ln

show_animation(animate(update, frames=nfr, fig=fig, interval=100), player="video", dpi=60)   # a smooth video
""", explain=r"""
1. Each frame releases 25 more parcels; every parcel lands as a dot at its $(u,v)$.
2. `run` is the running mean of the products $uv$; the orange band is the expected value $\pm5$ standard errors, which narrows
   like $1/\sqrt{\text{parcels}}$ (C01).""")
nb.figure_notes(
    see="Left: dots accumulating into a tilted cloud. Right: an orange line that jumps about at first and then settles inside a "
        "narrowing band around −0.05 m²/s² (these parcels have correlation 0.5, so half the −0.1 of the tiny example).",
    read="No single parcel 'has' a Reynolds stress; the stress is the average of many exchanges, and it converges like any "
         "average — as $N^{-1/2}$.",
    change="…the shear were switched off: the same parcels would land in a round cloud and the orange line would settle on zero.")
scratch(r"""
# From scratch: the 2 × 2 covariance with explicit sums, and the Reynolds stress from it
u_p, v_p, l_p = out["u"], out["v"], out["l"]                  # the parcel samples of the first run (correlation 1)
du, dv = u_p - u_p.mean(), v_p - v_p.mean()                   # fluctuations about the sample means
mine = np.array([[np.sum(du*du), np.sum(du*dv)], [np.sum(dv*du), np.sum(dv*dv)]])/u_p.size   # mean(u_i u_j) by explicit sums
lib = TS.velocity_covariance(np.column_stack([u_p, v_p]))     # the library: samples as rows, components as columns
assert np.allclose(mine, lib)                                 # same numbers -> the library does exactly this
assert np.allclose(-1.2*mine, TS.reynolds_stress(np.column_stack([u_p, v_p]), rho0=1.2))   # Reynolds stress = −rho0 × covariance
assert abs(out["uv"] - (-np.mean(v_p*l_p)*2.0)) < 1e-12       # mean(uv) = −mean(v l) dU/dy, exactly (N50)
print(np.round(mine, 4), "\nanisotropy:", np.round(TS.anisotropy_tensor(np.diag([0.5, 0.2, 0.3])), 3).tolist())
""", r"""
1. The four explicit sums reproduce `TS.velocity_covariance`; multiplying by $-\rho_0$ reproduces `TS.reynolds_stress`.
2. The third `assert` is the parcel argument itself, exact for these samples.
3. The last line shows `TS.anisotropy_tensor` for normal stresses 0.5, 0.2, 0.3 m²/s²: the diagonal departs from zero by
   $+0.167,-0.133,-0.033$ — the streamwise component holds more than its third of the energy.""")
note("N52 [B]", r"""
**The mean temperature equation.** Averaging the heat equation the same way leaves one new term, the divergence of the
velocity–temperature correlation. ⚠️ This $\kappa$ is the *thermal diffusivity* [m²/s] (`kappa_th` in code), not the von Kármán
constant.""",
     equation=r"\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}+\frac{\partial}{\partial x_j}(\overline{u_jT'})=\kappa\frac{\partial^2\bar T}{\partial x_j^2}", ref="12.31")
note("N53 [B]", r"""
**The mean heat flux** has a molecular part (Fourier's law $\mathbf q=-k\nabla T$ (1.2) of Ch. 1, with $k$ the thermal
conductivity, `k_th`) and a turbulent part: $\rho_0C_p\big(\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}\big)=-\frac{\partial Q_j}{\partial x_j}$ with""",
     equation=r"Q_j=-k\frac{\partial\bar T}{\partial x_j}+\rho_0C_p\overline{u_jT'}", ref="12.32")
note("N54 [B]", r"""
**The turbulent heat flux** $\rho_0C_p\overline{u_jT'}$ [W/m²] is upward over a heated surface: rising parcels are warm
($w>0$ with $T'>0$).""")
P("P296", "kinematic vs dynamic fluxes", r"""
Observers report a heat flux $H$ [W/m²] and a wall stress $\tau_0$ [Pa]; the equations carry the *kinematic* fluxes
$\overline{wT'}$ [K m/s] and $u_*^2$ [m²/s²]. The link is the density (and heat capacity): $H=\rho c_p\overline{wT'}$,
$\tau_0=\rho u_*^2$. For a perfect gas the expansion coefficient is $\alpha=1/T$.""", code=r"""
rho_air, cp_air = 1.2, 1005.0                     # air density [kg/m³] and specific heat [J/(kg K)]
wT_kin = 0.1                                      # kinematic heat flux mean(w T') [K m/s]
print(rho_air*cp_air*wT_kin, "W/m²")              # H = rho c_p mean(wT') = 120.6 W/m²: a sunny afternoon over land
""")
code(r"""
Q_mol_turb = ch12.mean_heat_flux(-0.01, 0.1, 0.026, 1.2, 1005.0)   # Eq. (12.32): gradient −0.01 K/m, mean(wT') = 0.1 K m/s, k_th of air
H_parts = ch12.turbulent_heat_flux(0.5, 0.4, 0.5, 1.2, 1005.0)     # rho c_p r w_rms T_rms with w_rms = 0.5 m/s, T_rms = 0.4 K, r = 0.5
print(f"mean heat flux Q = {Q_mol_turb:.2f} W/m² (molecular part {0.026*0.01:.1e} W/m²);  from rms values: {H_parts:.1f} W/m²")
""", explain=r"""
1. `ch12.mean_heat_flux(gradT, uT, k_th, rho0, cp)` adds the two parts of $Q_j$: the molecular part is $2.6\times10^{-4}$ W/m², the
   turbulent part 120.6 W/m² — molecular conduction is irrelevant away from the surface.
2. `ch12.turbulent_heat_flux(w_rms, T_rms, r_wT, rho, cp)` builds the same 120.6 W/m² from rms fluctuations and their correlation
   coefficient: $1.2\times1005\times0.5\times0.5\times0.4$.""")
note("N55 [C]", r"""
**Passive scalar and mixture density.** A passive scalar (dye, a trace gas) is carried by the flow without changing it. If a
volume fraction $\tilde\upsilon$ of source fluid (density $\rho_s$) is mixed into ambient fluid (density $\rho$), the mixture
density is $\rho_m=\tilde\upsilon\rho_s+(1-\tilde\upsilon)\rho$, and the *mass* fraction is $\tilde Y=\tilde\upsilon\rho_s/\rho_m$.
Used for jet dilution in C09.""")
note("N56 [B]", r"""
**Conservation of a scalar** with mass fraction $\tilde Y$ and molecular diffusivity $\kappa_m$ [m²/s].""",
     equation=r"\frac{\partial}{\partial t}(\rho_m\tilde Y)+\frac{\partial}{\partial x_j}(\rho_m\tilde Y\tilde u_j)=\frac{\partial}{\partial x_j}\Big(\rho_m\kappa_m\frac{\partial\tilde Y}{\partial x_j}\Big)", ref="12.33")
note("N57 [B]", r"""
**The mean scalar equation** — the D06 moves a third time, at constant $\rho_m$. The new term is the turbulent scalar flux
$\overline{u_jY'}$.""",
     equation=r"\frac{\partial\bar Y}{\partial t}+U_j\frac{\partial\bar Y}{\partial x_j}=\frac{\partial}{\partial x_j}\Big(\kappa_m\frac{\partial\bar Y}{\partial x_j}-\overline{u_jY'}\Big)", ref="12.34")
code(r"""
print(f"mixture density {ch12.mixture_density(0.1, 1.8, 1.2):.2f} kg/m³;  mass fraction {ch12.mass_fraction_from_volume_fraction(0.1, 1.8, 1.2):.4f}")   # 10 % by volume of a gas of 1.8 kg/m³ in air
flux = ch12.mean_scalar_flux(-0.02, 0.001, 2e-5)              # −kappa_m dY/dx + mean(uY'): gradient −0.02 1/m, mean(uY') = 0.001 m/s
print(f"total scalar flux {flux:.4e} m/s: molecular {2e-5*0.02:.0e}, turbulent 1e-03 — a ratio of {0.001/(2e-5*0.02):.0f}")
""", explain=r"""
1. A 10 % volume fraction of a heavier gas is a 14.3 % mass fraction: the two differ whenever the densities do.
2. `ch12.mean_scalar_flux(gradY, uY, kappa_m)` adds molecular and turbulent fluxes; here the turbulent one is 2500 times larger.""")
note("N58 [B]", r"""
**RANS** (Reynolds-averaged Navier–Stokes) is the name for the set: mean continuity $\partial U_i/\partial x_i=0$ (12.27), the
mean momentum equation
$\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\frac1{\rho_0}\frac{\partial\bar\tau_{ij}}{\partial x_j}$ (12.30),
the mean heat equation $\rho_0C_p\big(\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}\big)=-\frac{\partial Q_j}{\partial x_j}$ (12.32)
and the mean scalar equation
$\frac{\partial\bar Y}{\partial t}+U_j\frac{\partial\bar Y}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\kappa_m\frac{\partial\bar Y}{\partial x_j}-\overline{u_jY'}\big)$ (12.34).

> ⚠️ **An assumption that switches silently.** From §12.8 to §12.10 the book drops gravity and the subscript 0:
> $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-\frac1\rho\frac{\partial P}{\partial x_i}+\frac{\partial}{\partial x_j}\big(2\nu\bar S_{ij}-\overline{u_iu_j}\big)$,
> with $P$ now the deviation from hydrostatic pressure. Buoyancy returns in §12.11.""")
note("N59 [B]", r"""
**A transport equation for the Reynolds stress** can be derived (multiply the fluctuation momentum equation for $u_i$ by $u_j$,
add the same with $i$ and $j$ swapped, average). We state it term by term:

$$\frac{\partial\overline{u_iu_j}}{\partial t}+U_k\frac{\partial\overline{u_iu_j}}{\partial x_k}+\frac{\partial\overline{u_iu_ju_k}}{\partial x_k}=-\overline{u_iu_k}\frac{\partial U_j}{\partial x_k}-\overline{u_ju_k}\frac{\partial U_i}{\partial x_k}-\frac1\rho\Big(\overline{u_i\frac{\partial p}{\partial x_j}}+\overline{u_j\frac{\partial p}{\partial x_i}}\Big)-2\nu\overline{\frac{\partial u_i}{\partial x_k}\frac{\partial u_j}{\partial x_k}}+\nu\frac{\partial^2\overline{u_iu_j}}{\partial x_k^2}+g\alpha\big(\overline{u_jT'}\delta_{i3}+\overline{u_iT'}\delta_{j3}\big)\qquad\text{(12.35)}$$

| Term | Name | Colour in this chapter |
|---|---|---|
| $\partial_t\overline{u_iu_j}+U_k\partial_k\overline{u_iu_j}$ | change following the mean flow | teal |
| $\partial_k\overline{u_iu_ju_k}$ | turbulent transport (a **triple** correlation: new unknown) | teal |
| $-\overline{u_iu_k}\partial_kU_j-\overline{u_ju_k}\partial_kU_i$ | production by mean shear | orange |
| $-\frac1\rho(\overline{u_i\partial_jp}+\overline{u_j\partial_ip})$ | pressure–velocity correlation (redistributes energy among components) | grey |
| $-2\nu\overline{\partial_ku_i\,\partial_ku_j}$ | viscous dissipation | rose |
| $\nu\partial_k^2\overline{u_iu_j}$ | viscous diffusion | rose |
| $g\alpha(\overline{u_jT'}\delta_{i3}+\overline{u_iT'}\delta_{j3})$ | buoyant production | blue |

Half of its trace ($i=j$, summed) is the turbulent kinetic-energy budget of C06,
$\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47);
D10 derives that half-trace in full.""")
code(r"""
rsb = ch12.reynolds_stress_budget_sympy()                     # sympy derives the budget of each of the six components (cached)
print("each of the six component budgets equals Eq. (12.35):", rsb["checks"], "| all:", rsb["all_ok"])
uu_ex = np.array([[0.5, -0.1, 0.0], [-0.1, 0.2, 0.0], [0.0, 0.0, 0.3]])   # an example mean(u_i u_j) [m²/s²]: streamwise 0.5, uv = −0.1
gradU = np.zeros((3, 3)); gradU[0, 1] = 2.0                   # simple shear: dU_1/dx_2 = 2 1/s, every other gradient zero
P_ij = ch12.reynolds_stress_production(uu_ex, gradU)          # the production term of (12.35) for every component [m²/s³]
print("production tensor P_ij:\n", P_ij)
print(f"half its trace = {0.5*np.trace(P_ij):.2f} = shear production of kinetic energy {ch12.shear_production(uu_ex, gradU):.2f} m²/s³")
""", explain=r"""
1. `ch12.reynolds_stress_budget_sympy()` repeats the derivation symbolically for the six independent components and compares each
   with the stated equation: all `True`.
2. `ch12.reynolds_stress_production(uu, gradU)` evaluates the orange term. In a simple shear only $P_{11}=0.4$ and $P_{12}=-0.4$
   m²/s³ are non-zero: **shear feeds only the streamwise normal stress**; the other two components get their energy from the pressure
   term.
3. Half the trace, 0.2 m²/s³, is the production of turbulent kinetic energy `ch12.shear_production` — the term C06 is about.""")
note("N60 [B]", r"""
**The closure problem.** Count equations and unknowns. The mean equations are 4 (continuity + 3 momentum) for 4 + 6 = 10
unknowns ($U_i$, $P$ and six $\overline{u_iu_j}$). The Reynolds-stress equation of N59 above gives 6 more equations but brings 10
triple correlations plus pressure and dissipation correlations. Writing an equation for the triple correlations brings the
quadruple ones. **Averaging always loses.** Three responses: model the unknown correlations (RANS models: C12, C13); compute
every scale (direct numerical simulation, DNS: its cost is in C07); or compute the large eddies and model the small ones
(large-eddy simulation, LES — named in Ch. 10).""")
code(r"""
print(pd.DataFrame(ch12.closure_count(3))[["level", "equations", "unknowns", "new_moment", "closed"]].to_string(index=False))   # equations vs unknowns at three levels
""", explain=r"""
`ch12.closure_count(3)` counts equations and unknowns at three levels of averaging: 4 vs 10, 10 vs 20, 20 vs 35 (before even
counting the pressure and dissipation correlations). The gap widens at every level; `closed` is never `True`.""")
note("N61 [C]", r"""
**No Bernoulli integral for the mean flow.** With Reynolds stresses present the mean momentum equation has no Bernoulli
integral, even for steady "inviscid" mean flow: the turbulence takes energy out of the mean motion (C06 says where it goes).""")
explainer("reynolds_stress_parcels", "How can fluctuations that average to zero push the mean flow?",
          "you release parcels in a shear and watch each land in the (u, v) scatter; turning the shear down, to zero and negative "
          "tilts the cloud and the stress bar follows — the sign argument becomes a count.",
          ["Set the shear to zero: the cloud is round and the orange bar vanishes.",
           "Reverse the shear: the cloud tilts the other way.",
           "Lower 'keeps its momentum' to 0.3: the cloud fattens, r_uv and the stress fall together.",
           "Click a parcel: its ℓ, u, v and uv."])
whatif(r"""…the turbulence had no preferred direction at all? Then the cloud is round, the shear stresses vanish, and the whole
two-point tensor collapses to one function — C05.""")


# =====================================================================================================================
# A.6 §12.6 Homogeneous Isotropic Turbulence — C05
# =====================================================================================================================
nb.section("12.6", "Homogeneous Isotropic Turbulence", intro=r"""
**What is this section about?** The simplest turbulence: statistically the same at every point and in every direction. Symmetry
then does most of the work — nine correlation functions become one, and the dissipation rate becomes one measurable gradient.""")
core("C05", r"Isotropic dissipation $\bar\varepsilon=-15\nu\overline{u^2}[d^2f/dr^2]_{r=0}=30\nu\overline{u^2}/\lambda_f^2=15\nu\overline{u^2}/\lambda_g^2$ (12.43) "
            r"from the correlation coefficients $f$ and $g$",
     "How much can symmetry alone tell us about turbulence?")
problem(r"""
Behind a grid in a wind tunnel the turbulence has no mean shear to remember and soon looks the same in every direction.
(**N62 [B]**: *homogeneous* = no preferred place; *isotropic* = no preferred direction; grid turbulence is the laboratory
approximation.) The small eddies of every flow look like this too (*local isotropy*), which is why the results of this section
are used far beyond wind tunnels — for instance to measure the dissipation rate in the ocean from one velocity gradient.""")
idea("", r"""
Take two velocity vectors a distance $r$ apart. Only two arrangements are really different: both components **along** the line
joining the points (the longitudinal correlation $f$), or both **across** it (the transverse correlation $g$).
Incompressibility ties $g$ to $f$. So one curve $f(r)$ holds everything. (**N71 [B]**: the sketch below shows the two
arrangements.)""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.0, 2.8))
fg_geometry_sketch(ax)                                        # drawing: the longitudinal (f) and transverse (g) arrangements
""",
    see="Two pairs of points a distance $r$ apart. In one pair the velocity components point along the separation "
        "(longitudinal, $f$); in the other they point across it (transverse, $g$).",
    read="In isotropic turbulence every two-point velocity correlation is a combination of these two functions of the distance "
         "$r$ alone — direction does not matter.",
    change="…the turbulence were sheared (C04): the correlation would also depend on the *direction* of the separation, and two "
           "functions would not be enough.")
remind("C05")
note("N63 [B]", r"""
**Isotropy seen in a scatter plot.** A round cloud of $(u,v)$ means $\overline{uv}=0$ — no shear stress; a cloud stretched
along $v=-u$ means $\overline{uv}<0$ (the fourth panel of the figure in C02, and the parcels of C04).""")
note("N64 [B]", r"""
**What homogeneity and isotropy say about single-point moments**: no moment changes from place to place,
$\frac{\partial}{\partial x_i}\overline{u_j^n}=0$; the three normal stresses are equal,
$\overline{u_1^2}=\overline{u_2^2}=\overline{u_3^2}$; and so are the moments of the three *longitudinal* gradients.""",
     equation=r"\overline{(\partial u_1/\partial x_1)^n}=\overline{(\partial u_2/\partial x_2)^n}=\overline{(\partial u_3/\partial x_3)^n}", ref="12.36")
note("N65 [B]", r"""
**The six cross (transverse) gradient moments are equal as well.**""",
     equation=r"\overline{(\partial u_1/\partial x_2)^n}=\overline{(\partial u_1/\partial x_3)^n}=\dots=\overline{(\partial u_3/\partial x_2)^n}", ref="12.37")
code(r"""
n3 = 48 if FAST else 64                                       # grid points per side of a periodic cube
L3 = 2*np.pi                                                  # box size [m] (so the wavenumbers are the integers 1, 2, …)
E3_raw = lambda K: K**4*np.exp(-2.0*(K/6.0)**2)               # shape of the prescribed energy spectrum (peak at K = 6 rad/m)
E3_norm = quad(E3_raw, 0, np.inf)[0]                          # its integral
E3 = lambda K: 1.5*E3_raw(K)/E3_norm                          # ∫E dK = 1.5 m²/s² = turbulent kinetic energy ē → each component has variance 1
u3, v3, w3 = TS.synthetic_solenoidal_field(n3, L3, spectrum=E3, seed=0, dim=3)   # a divergence-free 3-D field (kinematic — no cascade)
dx3 = L3/n3                                                   # grid spacing [m]
rep = ch12.isotropy_report(u3, v3, w3, dx3)                   # single-point checks of Eqs. (12.36)–(12.37)
print("normal-stress ratios u_i²/mean:", np.round(rep["normal_stress_ratios"], 3))
print("shear correlation coefficients:", {k: round(float(v), 4) for k, v in rep["shear_coefficients"].items()})
print("mean square of the 3 longitudinal gradients:", np.round(rep["longitudinal"], 2), "\n   and of the 6 transverse ones:", np.round(rep["transverse"], 2))
K_sh, E_sh = TS.shell_spectrum((u3, v3, w3), L3)              # energy per unit wavenumber, measured from the field
print(f"sum of the measured spectrum {np.sum(E_sh)*(K_sh[1] - K_sh[0]):.4f} = ½ mean(u_i u_i) = {0.5*np.mean(u3**2 + v3**2 + w3**2):.4f} m²/s²")
smp = np.random.default_rng(0).normal(size=(100_000, 3))     # 100 000 samples of three independent unit-variance components
print(f"turbulent kinetic energy of unit-variance samples: {TS.turbulent_kinetic_energy(smp):.3f} = 3/2 × 1")
""", explain=r"""
1. `TS.synthetic_solenoidal_field(..., dim=3)` builds a random divergence-free field with the spectrum we prescribe (primer P286) —
   a stand-in for isotropic turbulence that has the right *kinematics* and no dynamics.
2. `ch12.isotropy_report` returns the three normal stresses over their mean (all ≈ 1), the three shear correlation coefficients
   (all ≈ 0), and the mean squares of the three longitudinal and six transverse gradients: equal within each group, as
   $\overline{(\partial u_1/\partial x_1)^n}=\overline{(\partial u_2/\partial x_2)^n}=\overline{(\partial u_3/\partial x_3)^n}$ (12.36) and its transverse
   twin say — and the transverse ones are **twice** the longitudinal ones, which D08 will derive.
3. `TS.shell_spectrum` measures the spectrum back from the field; its sum is half the mean square velocity to round-off: the field
   has the energy we prescribed.
4. `TS.turbulent_kinetic_energy` of three unit-variance components is $1.5=\tfrac32\overline{u^2}$ (**N72** below).""")
note("N66 [B]", r"""
**The two correlation coefficients.** $u_\parallel$ is the velocity component along the separation $\mathbf r$, $u_\perp$ a
component perpendicular to it; $\overline{u_\parallel^2}=\overline{u_\perp^2}=\overline{u^2}$ and $f(0)=g(0)=1$:
$g(r)\equiv\overline{u_\perp(\mathbf x+\mathbf r)u_\perp(\mathbf x)}\big/\overline{u_\perp^2}$ and""",
     equation=r"f(r)\equiv\overline{u_\parallel(\mathbf x+\mathbf r)u_\parallel(\mathbf x)}\big/\overline{u_\parallel^2}", ref="12.38")
note("N67 [B]", r"""
**Their integral scales and Taylor microscales** are C02's definitions with the distance $r$ in place of the lag $\tau$.""",
     equation=r"\Lambda_f\equiv\int_0^\infty f\,dr,\quad\Lambda_g\equiv\int_0^\infty g\,dr,\quad\lambda_f^2\equiv-2/[d^2f/dr^2]_{r=0},\quad\lambda_g^2\equiv-2/[d^2g/dr^2]_{r=0}", ref="12.39")
note("N72 [B]", r"""
**Energy and the one-component convention.** At zero separation the trace of the correlation tensor is twice the turbulent
kinetic energy, $R_{ii}(0)=\overline{u_iu_i}=2\bar e$. ⚠️ In this section $\overline{u^2}$ is **one** component, so
$\bar e=\tfrac32\overline{u^2}$; and $\bar e$ is not Ch. 1's internal energy $e$.""")
P("P297", "derivatives of functions of r = ∣r∣", r"""
The distance is $r=\sqrt{r_kr_k}$ (sum over $k$). Differentiating $r^2=r_kr_k$ gives $\partial r/\partial r_j=r_j/r$. Also
$\partial r_i/\partial r_j=\delta_{ij}$, $\delta_{jj}=3$ and $r_jr_j=r^2$. So for any function of the distance alone, the chain
rule gives $\partial F(r)/\partial r_j=F'(r)\,r_j/r$.""", code=r"""
r1, r2, r3 = sp.symbols("r1 r2 r3", positive=True)            # the three components of the separation
r_len = sp.sqrt(r1**2 + r2**2 + r3**2)                        # its length r
print(sp.simplify(sp.diff(r_len, r1) - r1/r_len))             # 0: ∂r/∂r_1 = r_1/r
print(sum(sp.diff(ri, ri) for ri in (r1, r2, r3)))            # 3: the divergence of the vector r_i is δ_jj = 3
""")
D07_CHECK = r"""
import sympy as sp                                             # symbolic algebra
r1, r2, r3, u2, s = sp.symbols("r1 r2 r3 u2 s", positive=True) # separation components, one-component variance, and the distance s = |r|
rv = [r1, r2, r3]                                              # the separation vector
f = sp.Function("f")(s)                                        # the longitudinal correlation, left completely general
g = f + s/2*sp.diff(f, s)                                      # step 11: g = f + (r/2) f'
R = sp.Matrix(3, 3, lambda i, j: u2*((f - g)/s**2*rv[i]*rv[j] + g*sp.KroneckerDelta(i, j)))   # steps 2 and 6: F r_i r_j + G δ_ij, with F, G from f, g
# step 7: the divergence ∂R_ij/∂r_j. R depends on r_j directly AND through the distance s, with ∂s/∂r_j = r_j/s (primer P297)
div = [sum(sp.diff(R[i, j], rv[j]) + sp.diff(R[i, j], s)*rv[j]/s for j in range(3)) for i in range(3)]
div = [sp.simplify(d.subs(r1**2 + r2**2 + r3**2, s**2).subs(r3, sp.sqrt(s**2 - r1**2 - r2**2))) for d in div]   # use s² = r_k r_k
print("divergence of the tensor for ANY f:", div)              # [0, 0, 0]
assert div == [0, 0, 0]                                        # the incompressible tensor of step 12 is divergence-free
gfree = sp.Function("g")(s)                                    # the converse (steps 8–11): leave g FREE and ask what zero divergence demands
Rfree = sp.Matrix(3, 3, lambda i, j: u2*((f - gfree)/s**2*rv[i]*rv[j] + gfree*sp.KroneckerDelta(i, j)))   # the tensor of step 6 with an unknown g
div1 = sum(sp.diff(Rfree[0, j], rv[j]) + sp.diff(Rfree[0, j], s)*rv[j]/s for j in range(3))               # its divergence, first component
bracket = sp.simplify((div1.subs(r3, sp.sqrt(s**2 - r1**2 - r2**2))/(u2*r1)).expand())                   # divide out the common factor u² r_1 (step 10)
print("divergence / (u² r_1) =", bracket)                      # f'/s + 2(f − g)/s²
g_solved = sp.solve(sp.Eq(bracket, 0), gfree)[0]               # step 11: the g that makes it vanish
assert sp.simplify(g_solved - g) == 0                          # exactly g = f + (r/2) f'
print("zero divergence forces g =", g_solved)
L = sp.symbols("L", positive=True); fG = sp.exp(-s**2/L**2)    # a Gaussian f to test steps 13 and 14
gG = fG + s/2*sp.diff(fG, s)                                   # its transverse partner
assert sp.simplify(sp.integrate(gG, (s, 0, sp.oo)) - sp.integrate(fG, (s, 0, sp.oo))/2) == 0   # step 13: Λ_g = Λ_f/2
assert sp.simplify(sp.diff(gG, s, 2).subs(s, 0) - 2*sp.diff(fG, s, 2).subs(s, 0)) == 0         # step 14: g''(0) = 2 f''(0), so λ_g = λ_f/√2
print("Gaussian f: Λ_g = Λ_f/2 and g''(0) = 2 f''(0) confirmed")
"""
NOTES_IN["D07"] = (r"**N68 [B]** the most general isotropic two-point tensor, $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40) (steps 1–2) · "
                   r"**N69 [B]** its incompressible form $R_{ij}=\overline{u^2}\big\{f\delta_{ij}+\frac r2\frac{df}{dr}\big(\delta_{ij}-\frac{r_ir_j}{r^2}\big)\big\}$ (12.41) "
                   r"(step 12) · **N70 [B]** $g=f+\frac r2f'$, $\Lambda_g=\Lambda_f/2$, $\lambda_g=\lambda_f/\sqrt2$ — ⚠️ three-dimensional results "
                   r"(steps 11, 13, 14)")
D("D07", ref="12.41", check_src=D07_CHECK, after=r"""
> ⚠️ **slip #14 — the book prints** (in its Exercise 12.18a) a reference to the scale definitions
> $\Lambda_f\equiv\int_0^\infty f\,dr,\ \dots,\ \lambda_g^2\equiv-2/[d^2g/dr^2]_{r=0}$ (12.39) as the starting tensor **; the correct
> form is** the general isotropic tensor $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40), from which this derivation starts.""")
code(r"""
Rt = ch12.isotropic_correlation_tensor(np.array([0.5, 0.0, 0.0]), lambda r: np.exp(-r**2))   # Eq. (12.41) for f = exp(−r²), separation 0.5 m along x
print(np.round(Rt, 4))                                        # diag(f, g, g) = diag(0.7788, 0.5841, 0.5841), off-diagonals 0
dv = ch12.isotropic_tensor_divergence_sympy()                 # sympy: the divergence ∂R_ij/∂r_j of Eq. (12.41) for a general f
print("divergence:", dv["divergence"], "| check:", dv["check"])
""", explain=r"""
1. `ch12.isotropic_correlation_tensor(rvec, f)` evaluates the incompressible isotropic tensor: for a separation along $x$,
   $R_{11}=f(0.5)=e^{-0.25}=0.7788$ and $R_{22}=R_{33}=g(0.5)=(1-r^2)e^{-r^2}=0.5841$.
2. `ch12.isotropic_tensor_divergence_sympy()` confirms that the tensor is divergence-free for **any** $f$ — the content of step 7 of
   the derivation.""")
fig(r"""
r_fg, f_meas, g_meas = TS.longitudinal_transverse_correlation(u3, v3, dx3, w3)   # f(r) and g(r) measured from the 3-D field, Eq. (12.38)
g_pred = ch12.transverse_from_longitudinal(r_fg, f_meas)      # g = f + (r/2) f' computed from the measured f (D07 step 11)
Lam_f = TS.integral_scale(r_fg, f_meas, upto="all")           # Λ_f = ∫ f dr over the whole range
Lam_g = TS.integral_scale(r_fg, g_meas, upto="all")           # Λ_g = ∫ g dr over the WHOLE range: g has a negative lobe that must be counted
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(r_fg, f_meas, color=C_MEAN, label=f"f(r) longitudinal, $\\Lambda_f$ = {Lam_f:.3f} m")
ax.plot(r_fg, g_meas, color=C_FLUC, label=f"g(r) transverse, $\\Lambda_g$ = {Lam_g:.3f} m")
ax.plot(r_fg, g_pred, "--", color=C_RS, label="$f+\\frac{r}{2}f'$ from the measured f")
ax.axhline(0, color=C_REF, lw=0.6)
ax.set(xlabel="separation r [m]", ylabel="correlation coefficient [–]", xlim=(0, 1.5), title="Incompressibility ties g to f: g dips below zero")
ax.legend(fontsize=8)
print(f"Lambda_g/Lambda_f = {Lam_g/Lam_f:.3f} (theory 0.5);  stopping at the first zero of g would give {TS.integral_scale(r_fg, g_meas)/Lam_f:.3f}")
""",
    see="A purple curve $f$ that falls from 1 to 0 without going negative; a teal curve $g$ that falls faster, dips below zero "
        "and returns; an orange dashed curve — $g$ predicted from $f$ alone — lying on the teal one through the dip (the "
        "prediction differentiates the measured $f$ on its 0.1 m grid, so it is a little rough, and at larger separations, where both "
        "are small, it wanders by a few hundredths).",
    read="Fluid that crosses the line joining the two points must come back somewhere (mass conservation), so the transverse "
         "correlation needs a negative lobe — large enough that its integral is exactly half of $f$'s. That is why we integrate "
         "$g$ over the **whole** range (`upto=\"all\"`): stopping at the first zero would leave the negative lobe out and "
         "overestimate $\\Lambda_g$.",
    change="…the field were two-dimensional: the relation would be $g=f+rf'$ instead (no factor ½), and $\\Lambda_g$ would be "
           "zero. The factor ½ and the results $\\Lambda_g=\\Lambda_f/2$, $\\lambda_g=\\lambda_f/\\sqrt2$ are three-dimensional.")
note("N73 [B]", r"""
**The dissipation rate of turbulent kinetic energy** is the average of Ch. 4's $\varepsilon=2\nu S_{ij}S_{ij}$ (4.58) for the
fluctuating strain rate $S'_{ij}=\tfrac12(\partial u_i/\partial x_j+\partial u_j/\partial x_i)$: a sum of nine squared gradients
(`ch12.dissipation_rate`).""",
     equation=r"\bar\varepsilon=2\nu\overline{S'_{ij}S'_{ij}}=\frac\nu2\overline{\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)^2}", ref="12.42")
NOTES_IN["D08"] = (r"**N74 [B]** gradient correlations are curvatures of the two-point tensor at zero separation, "
                   r"$\overline{\frac{\partial u_i}{\partial x_k}\frac{\partial u_j}{\partial x_l}}=-\big(\frac{\partial^2R_{ij}}{\partial r_k\partial r_l}\big)_{r=0}$, "
                   r"(steps 6–7), and the three kinds of moment stand in the ratio 2 : 4 : −1 (step 11)")
D("D08", ref="12.43", check_src=PF["D08"]["check_src"] + r"""
print("D08 check — moments × λ_f²/u² (longitudinal, transverse, cross):", [sp.simplify(m_*lam**2/u2) for m_ in (a, b, c)])   # 2, 4, −1
print("D08 check — ε̄ λ_f²/(ν u²) from all 18 terms:", sp.simplify(eps*lam**2/(nu*u2)), "| cross sum:", sp.simplify(sum(M(i, j, j, i) for i in range(3) for j in range(3))))   # 30 and 0
""")
nb.worked_example("dissipation from a Taylor microscale", r"""
Air, $\nu=1.5\times10^{-5}$ m²/s; $\overline{u^2}=1$ m²/s² (one component); $\lambda_f=1$ cm.

1. $\bar\varepsilon=30\nu\overline{u^2}/\lambda_f^2=30\times1.5\times10^{-5}\times1/10^{-4}=4.5$ m²/s³.
2. $\lambda_g=\lambda_f/\sqrt2=7.07$ mm; $15\times1.5\times10^{-5}/(5\times10^{-5})=4.5$ ✓.
3. One gradient: $\overline{(\partial u_1/\partial x_1)^2}=2\overline{u^2}/\lambda_f^2=2\times10^4$ s⁻²; $15\nu\times2\times10^4=4.5$ ✓.
4. $R_\lambda=\lambda_g\sqrt{\overline{u^2}}/\nu=471$.
5. Kolmogorov scale (C07): $\eta=(\nu^3/\bar\varepsilon)^{1/4}=0.165$ mm — 43 times smaller than $\lambda_g$.""")
code(r"""
nu_air = 1.5e-5                                               # kinematic viscosity of air [m²/s]
e1 = ch12.dissipation_isotropic(nu_air, 1.0, lambda_f=0.01)   # Eq. (12.43): 30 ν u²/λ_f² with u² = 1 m²/s², λ_f = 1 cm
e2 = ch12.dissipation_isotropic(nu_air, 1.0, lambda_g=0.01/np.sqrt(2))   # 15 ν u²/λ_g²
e3 = ch12.dissipation_isotropic(nu_air, dudx_sq=2.0e4)        # 15 ν mean((∂u_1/∂x_1)²): one measurable gradient
print(f"dissipation three ways: {e1:.3f}, {e2:.3f}, {e3:.3f} m²/s³;  R_lambda = {ch12.taylor_reynolds_number(1.0, 0.01/np.sqrt(2), nu_air):.1f}")
mom = ch12.gradient_moments_isotropic_sympy()                 # sympy repeats D08 for a general f
print("moments (longitudinal : transverse : cross) =", mom["m11"], ":", mom["m12"], ":", mom["cross"], "| ε̄ λ_f²/(ν u²) =", mom["eps_factor"])
grads = ch12.velocity_gradient_samples(u3, v3, w3, dx3)       # all nine gradients at every grid node of the synthetic field
eps_full = ch12.dissipation_rate(grads, 1e-3)                 # Eq. (12.42): the full nine-gradient sum (ν = 10⁻³ m²/s)
eps_one = 15*1e-3*np.mean(grads[:, 0, 0]**2)                  # Eq. (12.43): 15 ν mean((∂u_1/∂x_1)²) from ONE gradient
print(f"synthetic field: full sum {eps_full:.4f}, one-gradient formula {eps_one:.4f} m²/s³ (ratio {eps_one/eps_full:.3f})")
""", explain=r"""
1. `ch12.dissipation_isotropic(nu, u2, lambda_f=…, lambda_g=…, dudx_sq=…)` accepts any one of the three forms; all give 4.5 m²/s³.
2. `ch12.gradient_moments_isotropic_sympy()` returns the moment ratios 2 : 4 : −1 (in units of $\overline{u^2}/\lambda_f^2$) and the
   factor 30.
3. On the synthetic field the nine-gradient sum and the one-gradient formula agree within about a per cent — the residual is the
   sampling scatter of one finite field. This is how dissipation is measured at sea: one shear probe, one gradient, the factor 15.""")
scratch(r"""
# From scratch: a Gaussian f on a grid — f''(0) by a centred difference, g by np.gradient, and the dissipation three ways
lam_f = 0.01                                                  # the Taylor microscale we build in [m]
rg = np.linspace(0.0, 0.06, 6001)                             # separations from 0 to 6 λ_f [m]
fG = np.exp(-rg**2/lam_f**2)                                  # a Gaussian longitudinal correlation: f''(0) = −2/λ_f²
h = rg[1]                                                     # grid step
f2 = 2*(fG[1] - fG[0])/h**2                                   # f''(0) by a centred difference (f is even: f(−h) = f(h))
gG = fG + 0.5*rg*np.gradient(fG, rg)                          # g = f + (r/2) f'
g2 = 2*(gG[1] - gG[0])/h**2                                   # g''(0) the same way
lam_f_mine, lam_g_mine = np.sqrt(-2/f2), np.sqrt(-2/g2)       # Eq. (12.39): microscales from the curvatures
eps_a = -15*nu_air*1.0*f2                                     # −15 ν u² f''(0)
eps_b = 30*nu_air*1.0/lam_f_mine**2                           # 30 ν u²/λ_f²
eps_c = 15*nu_air*1.0/lam_g_mine**2                           # 15 ν u²/λ_g²
sc = ch12.isotropic_scales(rg, fG)                            # the library: scales of f and of the g it implies
assert np.allclose([eps_a, eps_b, eps_c], ch12.dissipation_isotropic(nu_air, 1.0, lambda_f=lam_f), rtol=1e-3)   # all three = 4.5
assert np.isclose(sc["Lambda_ratio"], 0.5, atol=1e-3) and np.isclose(sc["lambda_ratio"], 1/np.sqrt(2), atol=1e-3)   # Λ_g/Λ_f, λ_g/λ_f
assert np.isclose(np.trapezoid(gG, rg)/np.trapezoid(fG, rg), 0.5, atol=1e-3)   # the area under g is half the area under f
print(f"λ_g/λ_f = {lam_g_mine/lam_f_mine:.4f}, Λ_g/Λ_f = {sc['Lambda_ratio']:.4f}, ε̄ = {eps_a:.3f}, {eps_b:.3f}, {eps_c:.3f} m²/s³")
""", r"""
1. The curvature at zero separation is found with a centred difference, using that $f$ is even.
2. `np.gradient` gives $f'$ on the grid, and with it $g=f+\frac r2f'$.
3. The asserts confirm the three forms of the dissipation, $\Lambda_g/\Lambda_f=0.5$ and $\lambda_g/\lambda_f=0.7071$ against
   `ch12.dissipation_isotropic` and `ch12.isotropic_scales`.""")
note("N75 [B]", r"""
**The Taylor-scale Reynolds number** — say which $\lambda$: the two choices differ by a factor $\sqrt2$.""",
     equation=r"R_\lambda\equiv\lambda_{(g\text{ or }f)}\sqrt{\overline{u^2}}/\nu", ref="12.44")
note("N76 [B]", r"""
**The one-dimensional wavenumber spectrum** is C03's pair with $(r_1,k_1)$ in place of $(\tau,\omega)$; its units are m³/s²
(primer P293).""",
     equation=r"S_{11}(k_1)=\frac1{2\pi}\int_{-\infty}^{+\infty}R_{11}(r_1)e^{-ik_1r_1}dr_1=\frac{\overline{u_1^2}}{2\pi}\int_{-\infty}^{+\infty}f(r_1)e^{-ik_1r_1}dr_1", ref="12.45")
note("N77 [B]", r"""
**The measurement recipe**, as a pipeline: a time record from one probe → Taylor's hypothesis turns it into a space record →
autocorrelation → cosine transform → $S_{11}(k_1)$. The cell does it for lines through our synthetic (periodic) field and
compares with the direct periodogram.""")
fig(r"""
line = u3[0, 0, :]                                            # u along one x-line of the field: what a towed probe would record
t_line = np.arange(line.size)*dx3/10.0                        # the probe moves at U0 = 10 m/s: sample times [s]
x_line, u_line = TS.taylor_frozen(t_line, line, 10.0)         # step 1: time → distance, x = U0 t
u_all = u3.reshape(-1, n3)                                    # every x-line of the field (to average the estimate)
rows3 = u_all[::7] - u_all[::7].mean(axis=1, keepdims=True)   # every 7th line, each with its own mean removed
# step 2: the PERIODIC (circular) correlation of each line, then the average over lines. Wiener–Khinchin (C03):
# the inverse FFT of |FFT|² is the correlation of a periodic record, exactly — no zero padding, because the field IS periodic
R_circ = np.mean([np.fft.ifft(np.abs(np.fft.fft(row))**2).real/n3 for row in rows3], axis=0)   # R_11(r_1) for r_1 = 0 … L
r_circ = np.arange(n3)*dx3                                    # separations over the whole period 0 … L [m]
k1 = np.arange(1, n3//2)*2*np.pi/L3                           # wavenumbers resolved by the box [rad/m]
S_corr = np.array([np.sum(R_circ*np.cos(k_*r_circ))*dx3/(2*np.pi) for k_ in k1])   # step 3: Eq. (12.45) as a sum over one period (exact for a periodic function)
S_per = np.mean([TS.periodogram(row, dx3)[1] for row in rows3], axis=0)   # the direct route: averaged periodograms
k_per = TS.periodogram(rows3[0], dx3)[0]                      # their wavenumbers
R_bias = np.mean([TS.autocorrelation(row, dx3, max_lag=n3//2, unbiased=False)[1] for row in rows3], axis=0)   # what a NON-periodic record would give: sums of n − m products divided by n
r_bias = np.arange(R_bias.size)*dx3                           # its separations 0 … L/2 [m]
S_bias = np.abs(TS.spectrum_from_correlation(r_bias, R_bias, k1))        # the cosine transform of that "biased" estimate
tail = R_circ[0]/(np.pi*L3*k1**2)                             # the tail the kink predicts: R(0)/(π L k²) (see the reading notes)
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.loglog(k_per[1:], S_per[1:], "o", ms=4, color=C_REF, label="periodogram of the lines")
ax.loglog(k1, np.abs(S_corr), color=C_FLUC, lw=2, label="periodic correlation → cosine transform")
ax.loglog(k1, S_bias, ":", color=C_RS, label="biased estimate of a non-periodic record")
ax.loglog(k1[k1 > 9], tail[k1 > 9], "--", color=C_REF, lw=0.8, label="$R_{11}(0)/(\\pi Lk_1^2)$")
ax.set(xlabel="$k_1$ [rad/m]", ylabel="$S_{11}$ [m³/s²]", ylim=(1e-12, 1), title="The two routes to $S_{11}(k_1)$ give the same spectrum")
shown = np.interp(k1, k_per, S_per) > 1e-12                   # the plotted range (above the bottom of the axis)
ratio = np.abs(S_corr[shown])/np.interp(k1[shown], k_per, S_per)
print(f"periodic correlation route / periodogram over the plotted range: between {ratio.min():.4f} and {ratio.max():.4f}")
print(f"biased estimate at k1 = 10 and 30 rad/m: {np.interp(10.0, k1, S_bias):.1e}, {np.interp(30.0, k1, S_bias):.1e};  R(0)/(π L k²) there: {np.interp(10.0, k1, tail):.1e}, {np.interp(30.0, k1, tail):.1e} m³/s²")
far = k1 >= 14.0                                              # the range where only the tail is left
print(f"biased estimate / guide for k1 ≥ 14 rad/m: between {(S_bias[far]/tail[far]).min():.2f} and {(S_bias[far]/tail[far]).max():.2f}")
ax.legend(fontsize=8)
""", explain=r"""
1. `u_all` holds every $x$-line of the field; `rows3` is every seventh line with its mean removed.
2. `R_circ` is the correlation route: FFT each line, take $\lvert\cdot\rvert^2$, transform back, divide by the number of points. This is
   the FFT correlation of primer P289 **without** zero padding — right here because the field is periodic, so a shift that wraps
   round the box is a genuine shift.
3. `S_corr` is the cosine sum of that correlation over one period — the discrete form of
   $S_{11}(k_1)=\frac1{2\pi}\int R_{11}(r_1)e^{-ik_1r_1}dr_1$; `S_per` is the direct route, `TS.periodogram` averaged over the lines.
4. `R_bias`, `S_bias` repeat the recipe with the estimator a non-periodic record needs (`TS.autocorrelation(..., unbiased=False)`:
   the sum of the $n-m$ available products divided by $n$); `tail` is the guide line $R_{11}(0)/(\pi Lk_1^2)$.
5. The three printed lines: the ratio of the two routes over the plotted range; the biased estimate and the guide at two
   wavenumbers; and their ratio over the tail.""",
    see="Grey circles (periodogram) with a bold teal line (periodic correlation, then cosine transform) through them over the "
        "whole range: a plateau at small $k_1$, then a fall of ten decades beyond $k_1\\approx6$ rad/m, where our prescribed spectrum "
        "ends. An orange dotted curve — the same recipe with the estimator one would use for a non-periodic record — follows them "
        "down to about $10^{-3}$ and then bends onto a line of slope −2, close to the thin dashed guide $R_{11}(0)/(\\pi Lk_1^2)$ "
        "(the cell prints how close for $k_1\\ge14$ rad/m).",
    read="A one-dimensional spectrum is flat at small $k_1$ even though the three-dimensional spectrum $E(K)$ vanishes there: a "
         "straight cut through a field of eddies sees large-$K$ eddies obliquely, as long waves (aliasing by geometry). This is "
         "why measured $S_{11}$ never shows the energy-containing peak sharply. The two routes are the same information (the Wiener–Khinchin theorem of C03): "
         "our field is periodic, so its circular correlation is exact and the two spectra agree to round-off (the cell prints the "
         "ratio). **A measured time series is not periodic.** Its usual 'biased' correlation estimate divides the sum of the "
         "$n-m$ available products by $n$, which multiplies the true correlation by the triangle $1-r/L$. A triangle has a kink "
         "at $r=0$, and a kink transforms into a tail $R_{11}(0)/(\\pi Lk_1^2)$ that falls only like $k_1^{-2}$ — the orange dotted "
         "curve. At $k_1=30$ rad/m it is still about $6\\times10^{-5}$ m³/s² (the guide gives $5\\times10^{-5}$; the cell prints both), "
         "where the true spectrum has fallen below $10^{-12}$. So a spectrum from a measured "
         "correlation cannot be trusted below that floor; a longer record ($L$ larger) lowers it.",
    change="…the probe speed were comparable with the fluctuations: step 1 (Taylor's hypothesis) would fail and the spectrum "
           "would be smeared along $k_1$ (C03, N40).")
whatif(r"""…there *is* a mean shear? Then the cloud of C04 tilts, the turbulence can take energy from the mean flow — and we need
a budget for that energy (C06).""")

# =====================================================================================================================
# A.7 §12.7 Turbulent Energy Cascade and Spectrum — C06, C07, C08
# =====================================================================================================================
nb.section("12.7", "Turbulent Energy Cascade and Spectrum", intro=r"""
**What is this section about?** Energy. The mean flow loses it to the large eddies, they hand it to smaller ones, and the
smallest turn it into heat. We write the two budgets, find the size of the smallest eddies, and derive the most famous formula
of the subject, the −5/3 spectrum.""")
core("C06", r"The turbulent kinetic-energy budget $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47)",
     "Where does the energy of the turbulence come from, and where does it go?")
problem(r"""
Switch off the fan and the turbulence in a room dies within seconds: it needs feeding. In a river the feed is the mean current
rubbing against the bed; over a hot road it is buoyancy. An ocean-mixing or boundary-layer scheme in a climate model is, at
heart, a bookkeeping of this energy.""")
idea(r"""
mean flow  Ē = ½U_i²   ──[ shear production  −mean(u_i u_j) ∂U_i/∂x_j ]──►  turbulence  ē = ½ mean(u_i²)  ──[ ε̄ ]──► heat
     │  (tiny direct viscous loss 2ν S̄S̄ ~ 1/Re of production)                 ▲  buoyancy gα·mean(wT′): + heated below, − stable
the SAME term appears in both budgets with OPPOSITE signs
""", r"""
Symbols: $\bar E=\tfrac12U_i^2$ is the kinetic energy of the mean flow and $\bar e=\tfrac12\overline{u_i^2}$ that of the
turbulence, both per unit mass [m²/s²]; $\bar\varepsilon$ is the dissipation rate [m²/s³ = W/kg]. Colours: shear production
orange, dissipation rose, buoyancy blue, mean purple.""")
remind("C06")
NOTES_IN["D09"] = (r"**N78 [B]** the kinetic-energy budget of the mean flow, the Result of this derivation (steps 1–9)")
D("D09", ref="12.46")
note("N79 [B]", r"""
**Reading the mean-flow budget.** Divergence terms only move energy from place to place (Gauss's theorem, Ch. 2: their volume
integral is a surface flux). The exchange term sees only the symmetric part of the mean gradient,
$\overline{u_iu_j}\,\partial U_i/\partial x_j=\overline{u_iu_j}\bar S_{ij}$ (a symmetric tensor contracted with an antisymmetric
one gives zero, Ch. 2). For a shear flow $U(y)$ it is $\overline{uv}\,dU/dy$ — negative, a **loss** for the mean flow, because
$\overline{uv}<0$ when $dU/dy>0$ (C04).""")
note("N80 [B]", r"""
**The mean flow's direct viscous loss is negligible** next to its loss to turbulence:
$\frac{2\nu\bar S_{ij}\bar S_{ij}}{\overline{u_iu_j}(\partial U_i/\partial x_j)}\sim\frac{\nu(U/L)^2}{u_{rms}^2(U/L)}\sim\frac\nu{UL}=\frac1{\mathrm{Re}}\ll1$
(with $u_{rms}\sim U$). At $\mathrm{Re}=10^6$ the ratio is $10^{-6}$ (`ch12.mean_to_turbulent_dissipation_ratio(1e6)`).

⚠️ **Where this holds.** The estimate takes the mean gradient to be $\Delta U/L$ — true in free shear flows (jets, wakes, mixing
layers) and in the outer part of a wall flow. There the mean flow does not dissipate its energy; it hands it to the turbulence.
**Next to a wall it is different**: in the viscous wall layer the mean gradient is $u_*^2/\nu$ (C10), far larger than
$\Delta U/L$, and a channel flow dissipates a third to a half of its pressure work there *directly* — a share that falls only
slowly (like $1/\ln\mathrm{Re}_\tau$) as the Reynolds number rises. The channel cell below prints it.""")
D("D10", ref="12.47", check_src=PF["D10"]["check_src"] + r"""
print("D10 check — step 9, triple term minus its divergence form:", sp.simplify(lhs9 - rhs9))                 # 0
print("D10 check — step 10, pressure term minus its divergence form:", sp.simplify(sum(u[i]*sp.diff(p, X[i]) for i in rng) - sum(sp.diff(p*u[i], X[i]) for i in rng)))   # 0
print("D10 check — steps 12–14, viscous term: left − right simplifies to", sp.simplify(lhs - rhs))           # 0
""", after=r"""
> ⚠️ **slip #2 — the book prints** the label "change of $\bar E$" under the left side of this budget **; the correct form is**
> "change of the turbulent energy $\bar e$" ($\bar E$ is the mean-flow energy of D09).

> ⚠️ **slip #10 — the book prints** the triple correlation of this budget, in its later reduced forms, as $\tfrac12\overline{ev}$
> and as $\overline{ew}$ **; the correct form is** $\tfrac12\overline{u_i^2u_j}$ throughout (C09's N128 and C14's D24 use it).""")
note("N81 [B]", r"""
**The three terms that are not transport.**

| Term | Sign | Name | Colour |
|---|---|---|---|
| $-\overline{u_iu_j}\,\partial U_i/\partial x_j$ | $+$ here, $-$ in the mean-flow budget of D09 above | shear production | orange |
| $\bar\varepsilon=2\nu\overline{S'_{ij}S'_{ij}}$ | always a loss ($>0$, and **not** small) | viscous dissipation | rose |
| $g\alpha\overline{u_3T'}$ | $>0$ for an upward heat flux, $<0$ for a downward one | buoyant production or destruction | blue |

This is the turbulent twin of Ch. 11's disturbance-energy equation
$\frac{d}{dt}\int\frac12u_i^2\,dV=-\int u_iu_j\frac{\partial U_i}{\partial x_j}\,dV-\Lambda$ (11.88), where the same
$-\overline{uv}\,dU/dy$ fed a growing wave. `ch12.shear_production` and `ch12.buoyant_production` evaluate the orange and blue
terms.""")
note("N83 [B]", r"""
**Isotropic turbulence has no shear production**: $\overline{u_iu_j}=\overline{u_1^2}\delta_{ij}$ contracts with
$\partial U_i/\partial x_j$ to $\overline{u_1^2}\,\partial U_i/\partial x_i=0$ by mean continuity. Turbulence must be
*anisotropic* — the tilted cloud of C04 — to feed on the mean flow (checked in the cell below).""")
nb.worked_example("production in a log layer", r"""
$u_*=0.3$ m/s (the friction velocity, defined in C10: $-\overline{uw}=u_*^2$ near the ground), height $z=10$ m, $\kappa=0.41$:
$-\overline{uw}=0.09$ m²/s², $dU/dz=u_*/(\kappa z)=0.0732$ s⁻¹.

1. Production $=0.09\times0.0732=6.59\times10^{-3}$ m²/s³ ($=u_*^3/\kappa z$).
2. Night, heat flux $\overline{wT'}=-0.02$ K m/s, $\alpha=1/300$ K⁻¹, $g\approx9.81$ m/s²: buoyancy term
   $g\alpha\overline{wT'}=-6.5\times10^{-4}$ m²/s³ — it removes about 10 % of the production.
3. Steady, no transport: $\bar\varepsilon=6.59\times10^{-3}-0.65\times10^{-3}=5.9\times10^{-3}$ m²/s³.""")
code(r"""
uu_log = np.array([[0.4, 0.0, -0.09], [0.0, 0.25, 0.0], [-0.09, 0.0, 0.15]])   # mean(u_i u_j) near the ground [m²/s²]: uw = −u*² = −0.09
gU_log = np.zeros((3, 3)); gU_log[0, 2] = 0.3/(KAPPA*10.0)    # the only mean gradient: dU/dz = u*/(κ z) = 0.0732 1/s
P_shear = ch12.shear_production(uu_log, gU_log)               # −mean(u_i u_j) ∂U_i/∂x_j [m²/s³]
P_buoy = ch12.buoyant_production(-0.02, 1/300.0, 9.81)        # g α mean(wT') for a downward (night-time) heat flux [m²/s³]
print(f"shear production {P_shear:.2e}, buoyant term {P_buoy:.2e}, dissipation if steady {P_shear + P_buoy:.2e} m²/s³")
assert abs(ch12.shear_production(0.3*np.eye(3), gU_log)) < 1e-15   # N83: isotropic stresses extract nothing from a shear
print("direct/turbulent loss of the mean flow at Re = 1e6:", ch12.mean_to_turbulent_dissipation_ratio(1e6))
ch = ch12.channel_mixing_length(1000, kappa=KAPPA, A_plus=26.0)   # a model channel at Re_tau = 1000 (mixing length with wall damping, C12)
bud = ch12.channel_energy_budget(1000, KAPPA, 26.0)           # both energy budgets on that channel, in wall units u*⁴/ν
yp = bud["yplus"]                                             # distance from the wall in wall units y⁺ = y u*/ν (C10)
i_pk = np.argmax(bud["production"])                           # where production is largest
print(f"production peaks at y+ = {yp[i_pk]:.1f} with P+ = {bud['production'][i_pk]:.3f};  −mean(uv)+ there = {ch['minus_uv_plus'][i_pk]:.2f}")
print(f"integrals over the half-channel: {({k: round(float(v), 4) for k, v in bud['integrals'].items()})}")
print(f"identity (pressure work = mean dissipation + production) residual: {bud['identity_residual']:.1e};  direct share {bud['direct_fraction']:.2f}")
shares = {Re_: ch12.channel_energy_budget(Re_, KAPPA, 26.0, n=200 if FAST else 400)["direct_fraction"] for Re_ in (180, 1000, 5200, 20000)}   # the same budget at four Reynolds numbers
print("direct viscous share of the pressure work at Re_τ =", {k: round(float(v), 2) for k, v in shares.items()})
""", explain=r"""
1. The first three lines are the tiny example: $6.59\times10^{-3}$, $-6.5\times10^{-4}$ and $5.9\times10^{-3}$ m²/s³.
2. The `assert` is **N83**: an isotropic stress tensor contracted with a shear gives zero production.
3. `ch12.channel_mixing_length(Re_tau, kappa, A_plus)` (with $\mathrm{Re}_\tau=\delta u_*/\nu$, the half-height in wall units) is a **model** channel flow (the mixing-length closure of C12 with wall
   damping — approximate against DNS, not exact). It returns profiles in *wall units* (C10): `yplus`, `Uplus`, `dUdy_plus`, and
   `minus_uv_plus` $=-\overline{uv}/u_*^2$, **minus** the correlation — positive, because $\overline{uv}<0$ where $dU/dy>0$.
4. `ch12.channel_energy_budget` evaluates both budgets on it. Its flat arrays (`pressure_work`, `mean_dissipation`, `production`)
   are magnitudes; the signed terms are in `bud["mean"]` and `bud["turbulence"]`.
5. Production peaks at $y^+\approx10.5$ with $P^+\approx0.245$ — where viscous and Reynolds stress are equal. The reason:
   $P^+=(-\overline{uv}^+)(dU^+/dy^+)$ is a product of two numbers whose sum is the total stress $1-y^+/\mathrm{Re}_\tau\approx0.99$, and
   such a product is largest when the two are equal: $(0.99/2)^2=0.245$, just under ¼. (DNS puts the peak near $y^+\approx12$.)
6. Integrated over the half-channel, the work of the pressure gradient equals direct dissipation by the mean flow plus production
   of turbulence. The direct share is **not** small: the last line prints it at four Reynolds numbers — about 0.56, 0.45, 0.37 and
   0.33 from $\mathrm{Re}_\tau=180$ to 20 000. It falls only like $1/\ln\mathrm{Re}_\tau$, because the direct dissipation sits in the
   viscous wall layer, whose contribution in wall units hardly changes, while the pressure work grows only as the bulk velocity
   $U_{bulk}^+$ (N80's ⚠️).""")
scratch(r"""
# From scratch: production = −mean(uv) dU/dy in wall units, with np.gradient and np.trapezoid
slope = np.gradient(ch["Uplus"], ch["yplus"])                 # dU⁺/dy⁺ by finite differences
P_mine = ch["minus_uv_plus"]*slope                            # P⁺ = (−mean(uv)⁺)(dU⁺/dy⁺)
assert np.allclose(P_mine[2:-2], bud["production"][2:-2], rtol=2e-2, atol=2e-4)   # the library uses the exact slope; ours is a finite difference
I_mine = np.trapezoid(P_mine, ch["yplus"])                    # production integrated across the half-channel
print(f"integrated production by hand {I_mine:.3f}, library {bud['integrals']['production']:.3f}")
assert np.isclose(I_mine, bud["integrals"]["production"], rtol=1e-2)
""", r"""
1. `np.gradient` differentiates the model velocity profile; multiplying by $-\overline{uv}^+$ gives the production.
2. It matches the library's `production` (which uses the model's exact slope) to the accuracy of the finite difference, and so does
   its integral.""")
code(r"""
delta_c, us_c, nu_c = 0.05, 0.02, 1.0e-6                      # a water channel: half-height 5 cm, u* = 2 cm/s, ν = 1e-6 m²/s → Re_τ = 1000
y_c = ch["y_over_delta"]*delta_c                              # heights [m]
U_c = ch["Uplus"]*us_c                                        # mean velocity [m/s]
uv_c = -ch["minus_uv_plus"]*us_c**2                           # mean(uv) [m²/s²]: NEGATIVE (minus the tabulated −mean(uv)⁺)
mb = ch12.mean_energy_budget(y_c, U_c, uv_c, nu_c)            # the mean-flow budget of D09 for U(y), term by term [m²/s³]
tb = ch12.tke_budget(y_c, U_c, uv_c, eps=-uv_c*np.gradient(U_c, y_c))   # the turbulence budget with ε̄ set equal to production (local equilibrium)
print("mean-flow budget, integrated over the half-channel [m³/s³]:", {k: float(f"{v:.3g}") for k, v in mb["integrals"].items()})
print(f"turbulence budget: production integral {np.trapezoid(tb['production'], y_c):.3g} m³/s³ = −(loss of the mean flow) {-mb['integrals']['loss_to_turbulence']:.3g}")
print("the two symbolic derivations agree with the stated budgets:", ch12.mean_energy_budget_sympy()["check"], ch12.tke_budget_sympy()["check"])
""", explain=r"""
1. The same model channel in SI units. `ch12.mean_energy_budget(y, U, uv, nu)` returns the five terms of the mean-flow budget for
   $U(y)$ (pressure work, viscous and turbulent transport, direct dissipation, loss to turbulence) and their integrals; the two
   transport terms integrate to (nearly) zero — they only move energy.
2. `ch12.tke_budget(y, U, uv, eps)` does the same for the turbulence. Its production integral equals the mean flow's loss to
   turbulence with the opposite sign: one term, two budgets.
3. `ch12.mean_energy_budget_sympy()` and `ch12.tke_budget_sympy()` are the sympy engines behind D09 and D10: both `True`.""")
fig(r"""
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 5.4), sharex=True)
m_, tb_ = bud["mean"], bud["turbulence"]                      # the signed terms of the two budgets (wall units)
a1.semilogx(yp, m_["pressure_work"], color=C_MEAN, label="pressure work (gain)")
a1.semilogx(yp, m_["transport"], color=C_REF, label="transport (moves energy toward the wall)")
a1.semilogx(yp, m_["viscous_dissipation"], color=C_VISC, label="direct viscous dissipation (loss)")
a1.semilogx(yp, m_["loss_to_turbulence"], color=C_RS, label="loss to turbulence $+\\overline{uv}\\,dU/dy$")
a1.set(ylabel="mean-flow budget [$u_*^4/\\nu$]", title="One term, two budgets: the mean flow's loss is the turbulence's income")
a1.legend(fontsize=7, loc="upper right")
a2.semilogx(yp, tb_["production"], color=C_RS, label="shear production $-\\overline{uv}\\,dU/dy$")
a2.semilogx(yp, tb_["dissipation_plus_transport"], color=C_VISC, label="dissipation + transport (the residual that closes the budget)")
a2.axvline(yp[i_pk], color=C_REF, ls=":")
a2.set(xlabel="distance from the wall $y^+$ [–]", ylabel="turbulence budget [$u_*^4/\\nu$]", xlim=(0.5, 1000))
a2.legend(fontsize=7)
""",
    see="Top (mean flow): a purple gain everywhere (the pressure work, of size $U^+/\\mathrm{Re}_\\tau\\approx0.02$ — too small to "
        "see on this ±1 scale, yet it supplies everything when summed over the whole channel), a rose loss concentrated at the wall, an orange loss peaking near "
        "$y^+\\approx10$, and a grey transport curve that carries energy from the core toward the wall. Bottom (turbulence): an "
        "orange gain that is the mirror image of the orange loss above, and a rose sink that balances it.",
    read="The two orange curves are the same function with opposite signs — the exchange term. Nearly all of it happens within "
         "$y^+<50$ of the wall. ⚠️ The rose curve in the bottom panel is the **residual** that closes the budget (dissipation "
         "plus transport together): a mixing-length model carries no information on how that sink splits.",
    change="…$\\mathrm{Re}_\\tau=5200$: in wall units the peak stays at $y^+\\approx10$ with a height just under ¼, and the direct "
           "viscous share of the total falls only from 0.45 to about 0.37 — a wall flow never hands *all* of its energy to the turbulence.")
note("N82 [B]", r"""
**Buoyant production is potential energy released.** A layer that is unstably stratified (heavy over light) lowers its centre
of mass when it mixes; the lost potential energy appears as turbulent kinetic energy. Mixing a stable layer *costs* energy.
⚠️ The temperature here is **potential** temperature. In terms of the thermometer temperature our unstable example cools by
15 K per km: Kundu's $dT/dz=-15<\Gamma_a=-9.8$ K/km ⇔ meteorology's $\Gamma_{met}=15>\Gamma_d=9.8$ K/km — unstable in both
conventions.""")
fig(r"""
z_mix = np.linspace(0.0, 1000.0, 201)                         # a 1 km deep layer [m]
Gam_a = STRAT.adiabatic_lapse_rate()                          # adiabatic gradient, Kundu sign [K/m] ≈ −9.8e-3
cases = {"unstable: dT/dz = −15 < Γa = −9.8 K/km  ⇔  Γ_met = 15 > Γd = 9.8 K/km": -15e-3,
         "stable: dT/dz = −6.5 > Γa = −9.8 K/km  ⇔  Γ_met = 6.5 < Γd = 9.8 K/km": -6.5e-3}   # two in-situ temperature gradients [K/m]
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for (name, dTdz), col in zip(cases.items(), (C_VISC, C_BUOY)):
    theta0 = 300.0 + (dTdz - Gam_a)*z_mix                     # potential temperature: its gradient is dT/dz − Γa
    res = ch12.mixing_potential_energy_change(z_mix, theta0, 1/300.0, 1.2)   # mix the column to its mean; α = 1/T, ρ0 = 1.2 kg/m³
    ax.plot(theta0, z_mix, color=col, label=f"{name}\n   mixing changes the potential energy by {res['dPE']:+.0f} J/m²")
    ax.plot(res["T_final"] + 0*z_mix, z_mix, ":", color=col)
    print(f"dT/dz = {dTdz*1e3:+.1f} K/km (Γ_met = {-dTdz*1e3:.1f}): dθ/dz = {(dTdz - Gam_a)*1e3:+.2f} K/km, dPE = {res['dPE']:+.0f} J/m², stable initially: {res['stable_initially']}")
ax.set(xlabel="potential temperature θ [K]", ylabel="height z [m]", title="Mixing an unstable layer releases potential energy")
ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5))
""",
    see="Two initial profiles of potential temperature (solid) and the uniform mixed states they end in (dotted). The rose "
        "profile decreases with height (unstable); the blue one increases (stable). The legend gives the change of potential "
        "energy: negative for the rose case, positive for the blue.",
    read="A negative change is energy handed to the turbulence — the blue term $g\\alpha\\overline{wT'}>0$ of the budget, with an "
         "upward heat flux. A positive change must be paid for by shear production: stable stratification is a sink (C14).",
    change="…the in-situ gradient were exactly adiabatic, $dT/dz=\\Gamma_a\\approx-9.8$ K/km ($\\Gamma_{met}=\\Gamma_d$): the "
           "potential temperature would be uniform and mixing would neither release nor cost energy.")
explainer("turbulent_energy_budget", "Where does turbulent energy come from and go?",
          "clicking a height shows the two budgets' bars side by side with the exchange term mirrored; sliding the friction "
          "Reynolds number shrinks the mean flow's direct viscous loss — 'same term, opposite sign' is seen, not read.",
          ["Click y⁺ ≈ 10: production is at its maximum, ¼ in wall units.",
           "Go to the centreline: production vanishes with the shear.",
           "Raise Re_τ from 180 to 5200 and watch the integrated direct dissipation's share fall.",
           "Turn the wall damping off and see where the peak moves."])
whatif(r"""…we ask how *much* energy flows, and how small the eddies that finally dissipate it are? C07.""")


# ---------------------------------------------------------------------------------------------------------------------
core("C07", r"Kolmogorov scales $\eta=(\nu^3/\bar\varepsilon)^{1/4},\ u_K=(\nu\bar\varepsilon)^{1/4}$ (12.50) with $\bar\varepsilon\sim(\Delta U)^3/L$ (12.49) "
            r"and $\eta/L\sim\mathrm{Re}_L^{-3/4}$ (12.51)",
     "If viscosity does the dissipating, why does the dissipation rate not depend on viscosity?")
problem(r"""
Stir a cup of water and a cup of honey with the same spoon at the same speed until both are turbulent (you would need a much
bigger spoon for honey). The power you put in is set by the spoon and the speed, not by the fluid. All of it ends as heat. So
the fluid must arrange for viscosity to remove exactly that power — by making eddies small enough. **Climate hook:** in the
atmosphere those eddies are about a millimetre across while a model grid box is kilometres wide: six to seven decades that can
only be parameterised.""")
idea(r"""
L (big eddies, speed ΔU) → L/2 → L/4 → …  each tier strained by the one above, no viscosity yet  … → η (viscous: Re = 1)
energy passes down at the SAME rate ε̄ ~ ΔU³/L through every tier:   u′(l′) ~ (ε̄ l′)^(1/3),   Re(l′) = u′l′/ν ↓
""", r"""
**N87 [B] — the cascade.** Big eddies are unstable and break into smaller ones, which break again. No tier dissipates
appreciably until the eddy Reynolds number $u'l'/\nu$ falls to about 1. So the rate at which energy is handed down is fixed at
the top, by the big eddies, and viscosity only decides how far down the ladder goes. Symbols: $L$ outer length [m], $\Delta U$
velocity difference across it [m/s], $\eta$ the Kolmogorov length [m], $u_K$ and $\tau_\eta$ the Kolmogorov velocity and time.""")
remind("C07")
note("N84 [B]", r"""
**The outer scale.** $L$ is the cross-stream extent over which the mean velocity changes by $\Delta U$ (the width of a jet, the
thickness of a boundary layer), so mean gradients are of order $\Delta U/L$ and the largest eddies have size $L$ and speed
$\Delta U$. The sketch shows eddies of many sizes inside one flow.""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.2, 3.2))
scales_sketch(ax, n_tiers=5, seed=0)                          # drawing: eddies of five sizes, from the outer scale L down
""",
    see="A few large swirls of size $L$, more of half that size inside them, and so on down to many tiny ones.",
    read="Turbulence is not eddies of one size but a whole population, each tier living inside (and strained by) the tier above.",
    change="…the Reynolds number were raised at fixed $L$ and $\\Delta U$: more tiers would be added at the small end — the "
           "large eddies would not change. The animation below shows it.")
NOTES_IN["D11"] = (r"**N85 [B]** the work done on the turbulence, $\dot W\sim\overline{u_iu_j}(\partial U_i/\partial x_j)\sim(\Delta U)^2[\Delta U/L]=(\Delta U)^3/L$ (12.48) "
                   r"(steps 2–3) · **N86 [B]** in a steady state all of it is dissipated, $\dot W=\bar\varepsilon$, so $\bar\varepsilon\sim(\Delta U)^3/L$ (12.49) (steps 1 and 3) · "
                   r"**N88 [B]** the Reynolds number built on the Kolmogorov scales is one, $\eta u_K/\nu=1$ (step 8) · **N89 [B]** the scale "
                   r"separation $\eta/L\sim\mathrm{Re}_L^{-3/4}$ with $\mathrm{Re}_L=\Delta UL/\nu$ (12.51) (step 9)")
D("D11", ref="12.50")
nb.worked_example("a stirred tank of water", r"""
$\Delta U=1$ m/s, $L=1$ m, $\nu=10^{-6}$ m²/s.

1. $\mathrm{Re}_L=\Delta UL/\nu=10^6$.
2. $\bar\varepsilon\sim\Delta U^3/L=1$ m²/s³ (1 W per kg).
3. $\eta=(\nu^3/\bar\varepsilon)^{1/4}=(10^{-18})^{1/4}=3.2\times10^{-5}$ m = 32 µm.
4. $u_K=(\nu\bar\varepsilon)^{1/4}=3.2$ cm/s; $\tau_\eta=(\nu/\bar\varepsilon)^{1/2}=1$ ms.
5. Check: $\eta u_K/\nu=3.2\times10^{-5}\times3.2\times10^{-2}/10^{-6}=1$ ✓; $\eta/L=3.2\times10^{-5}=(10^6)^{-3/4}$ ✓ — 4.5 decades.
6. A simulation resolving it all needs about $\mathrm{Re}^{9/4}=3\times10^{13}$ grid points.""")
code(r"""
eps_tank = ch12.dissipation_outer_scaling(1.0, 1.0)           # Eq. (12.49): ε̄ ~ ΔU³/L with ΔU = 1 m/s, L = 1 m [m²/s³]
eta, uK, tau_eta = ch12.kolmogorov_scales(1e-6, eps_tank)     # Eq. (12.50): Kolmogorov length [m], velocity [m/s], time [s] in water
print(f"ε̄ = {eps_tank} m²/s³;  η = {eta:.3e} m, u_K = {uK:.4f} m/s, τ_η = {tau_eta:.1e} s;  η u_K/ν = {eta*uK/1e-6:.3f}")
sep = ch12.scale_separation(1e6)                              # Eq. (12.51) and its relatives at Re_L = 10⁶
print(f"η/L = {sep['eta_over_L']:.3e} = Re^(-3/4);  DNS grid points ≈ {ch12.dns_grid_points(1e6):.2e} = Re^(9/4)")
tiers = ch12.cascade_tiers(1.0, 1.0, 1e-6)                    # the ladder: sizes halving from L down to η
print(f"{tiers['n_tiers']} tiers; eddy Reynolds numbers from {tiers['Re'][0]:.0e} down to {tiers['Re'][-1]:.1f}")
table = pd.DataFrame([ch12.scale_table(c) for c in ch12.SCALE_CASES])[["name", "L", "dU", "nu", "Re_L", "eps", "eta", "tau_eta", "decades"]]   # our four cases
print(table.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
""", explain=r"""
1. `ch12.dissipation_outer_scaling(dU, L)` is $\bar\varepsilon\sim(\Delta U)^3/L$ (12.49) with a constant of order one (1 by
   default); `ch12.kolmogorov_scales(nu, eps)` returns the three scales of the tiny example: 32 µm, 3.2 cm/s, 1 ms.
2. `ch12.cascade_tiers` lists the ladder: each tier half the size of the one above, speed $(\bar\varepsilon l')^{1/3}$, until the
   eddy Reynolds number reaches about 1.
3. **N90 [B] — our own table of cases** (`ch12.SCALE_CASES`): a kitchen mixer, a small wind tunnel, the atmospheric boundary layer
   ($\Delta U=5$ m/s, $L=1$ km: $\bar\varepsilon=0.125$ m²/s³, $\eta\approx0.4$ mm, 6.4 decades between $L$ and $\eta$) and the ocean
   thermocline ($\bar\varepsilon=10^{-4}$ m²/s³, $\eta\approx0.3$ mm). The Kolmogorov scale is a fraction of a millimetre almost
   everywhere in nature — only the *range* of scales changes.""")
scratch(r"""
# From scratch: the exponents of η = ν^a ε̄^b from dimensions alone (a 2 × 2 linear system)
# units: ν [m² s⁻¹], ε̄ [m² s⁻³], η [m¹ s⁰] →  metres: 2a + 2b = 1,  seconds: −a − 3b = 0
a_, b_ = np.linalg.solve([[2.0, 2.0], [-1.0, -3.0]], [1.0, 0.0])   # solve the two exponent equations
print(f"η = ν^{a_:.2f} ε̄^{b_:.2f}")                           # 0.75 and −0.25: η = (ν³/ε̄)^(1/4)
assert np.allclose([a_, b_], [float(x) for x in ch12.kolmogorov_exponents()["eta"]])   # the library's exact fractions 3/4, −1/4
grp = DIM.solve_exponents("eta", ["nu", "eps"], {"eta": {"L": 1}, "nu": {"L": 2, "T": -1}, "eps": {"L": 2, "T": -3}})   # Ch. 1's Π-theorem solver
print("Ch. 1 solver: the group", grp, "is dimensionless")      # η ν^(−3/4) ε̄^(1/4)
abl = ch12.scale_table("atmospheric_boundary_layer")          # the atmospheric case of the table
eta_mine = (abl["nu"]**3/abl["eps"])**0.25                    # Eq. (12.50) by hand
assert np.isclose(eta_mine, ch12.kolmogorov_scales(abl["nu"], abl["eps"])[0]) and np.isclose(eta_mine, abl["eta"])
print(f"atmospheric boundary layer: η = {eta_mine*1e3:.2f} mm")
""", r"""
1. Matching metres and seconds gives two equations for the two exponents; `np.linalg.solve` returns $a=\tfrac34$, $b=-\tfrac14$.
2. Chapter 1's `DIM.solve_exponents` finds the same group as a dimensionless product (so its exponents are the negatives).
3. The last lines evaluate $\eta=(\nu^3/\bar\varepsilon)^{1/4}$ by hand for the atmospheric boundary layer and compare with the
   library: 0.4 mm.""")
nb.animation(r"""
Re_list = 10.0**np.arange(3.0, 8.01, 1.0 if FAST else 0.5)    # Re_L from 10³ to 10⁸, one frame per (half-)decade
K_ax = np.logspace(0, 8, 300)                                 # wavenumbers K L [–] (L = 1 m)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 3.5))
a1.set(xlim=(2e-7, 3), ylim=(0, 1), xscale="log", yticks=[], xlabel="eddy size l / L [–]")
a1.invert_xaxis()                                             # big eddies on the left, small on the right
(tier_dots,) = a1.plot([], [], "o", color=C_FLUC, ms=5)       # the tiers of the cascade
vl_L = a1.axvline(1.0, color=C_MEAN, lw=2); vl_eta = a1.axvline(1.0, color=C_VISC, lw=2); vl_lam = a1.axvline(1.0, color=C_RS, lw=1.2, ls="--")
a1.text(1.0, 0.92, " L", color=C_MEAN); t_eta = a1.text(1.0, 0.92, " η", color=C_VISC); t_lam = a1.text(1.0, 0.8, " λ", color=C_RS)
(sp_line,) = a2.loglog([], [], color=C_FLUC)                  # the model spectrum
a2.loglog(K_ax, 1.5*K_ax**(-5/3), ":", color=C_REF, label="$C\\bar\\varepsilon^{2/3}K^{-5/3}$")
a2.set(xlim=(1, 1e8), ylim=(1e-14, 1), xlabel="wavenumber K L [–]", ylabel="E(K) / (ΔU² L) [–]")
a2.legend(fontsize=8)

def update(i):                                                # frame i: the ladder and the spectrum at Re_list[i]
    Re = Re_list[i]
    tr = ch12.cascade_tiers(1.0, 1.0, 1.0/Re)                 # L = 1 m, ΔU = 1 m/s, ν = 1/Re
    so = ch12.scale_ordering(Re)                              # η/L and λ/L at this Reynolds number
    tier_dots.set_data(tr["size"], 0.5 + 0*tr["size"])
    vl_eta.set_xdata([so["eta_over_L"]]*2); t_eta.set_x(so["eta_over_L"])
    vl_lam.set_xdata([so["lambda_over_L"]]*2); t_lam.set_x(so["lambda_over_L"])
    sp_line.set_data(K_ax, ch12.model_spectrum(K_ax, 1.0, 1.0/Re, L=1.0, kind="pope"))
    a1.set_title(f"Re_L = 10^{np.log10(Re):.1f}: {tr['n_tiers']} tiers, η/L = {so['eta_over_L']:.1e}", fontsize=10)
    a2.set_title(f"−5/3 range ≈ {ch12.inertial_range_decades(Re):.1f} decades (L to η)", fontsize=10)
    return tier_dots, sp_line

show_animation(animate(update, frames=len(Re_list), fig=fig, interval=600), player="frames", dpi=40)   # step through the Reynolds numbers
""", explain=r"""
1. Each frame is one Reynolds number; $L$ and $\Delta U$ are fixed (so $\bar\varepsilon=1$ m²/s³ in every frame) and only $\nu$
   changes.
2. Left: the tiers from `ch12.cascade_tiers` on a logarithmic size axis, with $L$ (purple), the Taylor microscale $\lambda$
   (orange dashed) and $\eta$ (rose) from `ch12.scale_ordering`.
3. Right: a model spectrum for the whole range (`ch12.model_spectrum`, explained in C08) against the −5/3 line.""")
nb.figure_notes(
    see="Step through: on the left the purple line at $L$ never moves, the rose line at $\\eta$ slides to the right (smaller "
        "sizes) and teal tiers are added at the small end. On the right the spectrum's left end stays put while its straight "
        "−5/3 part grows longer toward large wavenumbers.",
    read="The left end never moves: a higher Reynolds number makes the cascade **longer, not stronger**. The energy supply "
         "$\\Delta U^3/L$ is the same in every frame; viscosity only sets where the ladder ends.",
    change="…we doubled $\\Delta U$ instead of lowering $\\nu$: the supply would rise eightfold *and* $\\eta$ would shrink by "
           "$2^{-3/4}=0.59$ — both ends would move.")
note("N91 [B]", r"""
**Where the Taylor microscale sits.** Eliminate $\bar\varepsilon$ between the isotropic relation
$\bar\varepsilon=15\nu\overline{u^2}/\lambda_g^2$ (12.43) and the outer estimate $\bar\varepsilon\sim(\Delta U)^3/L$ (12.49):
$\frac{(\Delta U)^3}L\propto\frac{\nu\overline{u^2}}{\lambda_T^2}\Rightarrow\frac{\lambda_T^2}{L^2}\propto\frac{\overline{u^2}}{(\Delta U)^2}\Big(\frac\nu{\Delta UL}\Big)\propto\frac1{\mathrm{Re}_L}$, or

$$\frac{\lambda_T}L\propto\mathrm{Re}_L^{-1/2}\qquad\text{(12.52)}$$

⚠️ $\lambda_T$ is **not** the size of the dissipating eddies: the relation $\bar\varepsilon=15\nu\overline{u^2}/\lambda_g^2$ pairs
it with the *large-eddy* velocity, not with $u_K$. It lies between the two ends, $\lambda_T/\eta\sim\mathrm{Re}_L^{1/4}$.""")
note("N92 [B]", r"""
**The ordering of the scales.** $R_\lambda\sim\mathrm{Re}_L^{1/2}$, and at high Reynolds number
$\eta<\lambda_T<\Lambda_{(f\text{ or }g)}<L$.""")
code(r"""
for Re in (1e2, 1e4, 1e6, 1e8):                               # four outer Reynolds numbers
    so = ch12.scale_ordering(Re)                              # η/L, λ_g/L, R_λ and whether η < λ_g < L
    print(f"Re_L = {Re:.0e}: η/L = {so['eta_over_L']:.1e}, λ/L = {so['lambda_over_L']:.1e}, λ/η = {so['lambda_over_L']/so['eta_over_L']:5.1f}, R_λ = {so['R_lambda']:.0f}, ordered: {so['ordered']}")
""", explain=r"""
`ch12.scale_ordering(Re_L)` evaluates the three ratios. Each factor 100 in $\mathrm{Re}_L$ divides $\eta/L$ by $100^{3/4}=31.6$,
$\lambda/L$ by 10, and multiplies $\lambda/\eta$ by $100^{1/4}=3.2$ and $R_\lambda$ by 10.""")
whatif(r"""…we ask how the energy is shared among the tiers in between? That is a spectrum, and dimensional analysis gives its
shape (C08).""")

# ---------------------------------------------------------------------------------------------------------------------
core("C08", r"Kolmogorov's inertial-subrange law $S_{11}(k_1)=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ for $2\pi/L\ll k_1\ll2\pi/\eta$ (12.54, slip #1 corrected)",
     "Why do spectra from a tidal channel, a jet and the atmosphere all fall with the same slope?")
problem(r"""
Plot the spectrum of the wind on a mast, of the current in a tidal channel, of the flow behind a grid: on log–log paper, over
the middle range of scales, all three are straight lines with the same slope. Something universal is going on. Eddies in that
range are too small to know about the geometry and too big to feel viscosity; all they know is the rate at which energy is
handed through them.""")
idea("", r"""
Only $\bar\varepsilon$ [m²/s³] and $k_1$ [1/m] are available; $S_{11}$ has units m³/s². There is exactly one way to build m³/s²
from them.

| Range | Wavenumbers | What sets the spectrum |
|---|---|---|
| energy-containing | $k_1\sim2\pi/L$ | the particular flow — not universal |
| inertial | $2\pi/L\ll k_1\ll2\pi/\eta$ | $\bar\varepsilon$ alone: slope −5/3, universal |
| dissipation | $k_1\eta\sim1$ | $\bar\varepsilon$ and $\nu$: universal in Kolmogorov units |""")
remind("C08")
NOTES_IN["D12"] = (r"**N93 [B]** the universal form at high wavenumber, $\frac{S_{11}(k_1)}{\nu^{5/4}\bar\varepsilon^{1/4}}=\Phi\big(\frac{k_1\nu^{3/4}}{\bar\varepsilon^{1/4}}\big)$, or "
                   r"$\frac{S_{11}(k_1)}{u_K^2\eta}=\Phi(k_1\eta)$ for $k_1\gg2\pi/L$ (12.53) (steps 3–4) · **N94 [B]** the normalisation "
                   r"$\int_{-\infty}^{+\infty}S_{11}(k_1)\,dk_1=\overline{u_1^2}=\int_0^{+\infty}2S_{11}(k_1)\,dk_1$ (12.55) (step 9)")
D("D12", ref="12.54")
slip(1, r"$S_{11}(k_1)=\mathit{const}\cdot\bar\varepsilon^{2/3}\cdot k_1^{5/3}$",
     r"$S_{11}(k_1)=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ (12.54). Three checks: units (D12 step 5 and the cell below); the sentence after "
     r"it in the book calls it the $k^{-5/3}$ law; and a spectrum that *rose* with $k_1$ would have infinite variance.")
code(r"""
k_two = np.array([10.0, 100.0])                               # two wavenumbers a decade apart [rad/m]
ok = ch12.inertial_spectrum_1d(k_two, 1.0)                    # Eq. (12.54), corrected: ε̄ = 1 m²/s³, two-sided constant
bad = ch12.inertial_spectrum_1d(k_two, 1.0, printed=True)     # the form as printed (+5/3), kept only to show it fails
print(f"corrected: S_11 = {ok} m³/s², ratio over one decade {ok[1]/ok[0]:.4f} = 10^(-5/3)")
print(f"printed:   ratio over one decade {bad[1]/bad[0]:.1f} — the spectrum would RISE 46-fold per decade")
dims = {"eps": np.array([Fraction(2), Fraction(-3)]), "k": np.array([Fraction(-1), Fraction(0)])}   # exponents of (metre, second) in ε̄ and k_1
for p in (Fraction(-5, 3), Fraction(5, 3)):                   # the corrected and the printed exponent of k_1
    m_exp, s_exp = Fraction(2, 3)*dims["eps"] + p*dims["k"]   # units of ε̄^(2/3) k_1^p
    print(f"ε̄^(2/3) k^({p}) has units m^({m_exp}) s^({s_exp})" + ("  = m³/s², the units of S_11 ✓" if (m_exp, s_exp) == (3, -2) else "  ✗ not m³/s²"))
""", explain=r"""
1. `ch12.inertial_spectrum_1d(k1, eps)` is the corrected law; `printed=True` reproduces the page so that we can watch it fail.
2. The loop does the unit arithmetic with exact fractions (`Fraction`, Ch. 1 P60): $(\mathrm m^2\mathrm s^{-3})^{2/3}(\mathrm m^{-1})^{-5/3}=\mathrm m^{3}\mathrm s^{-2}$ ✓,
   while the printed exponent gives $\mathrm m^{-1/3}\mathrm s^{-2}$ ✗.""")
note("N95 [B]", r"""
**The three-dimensional spectrum.** Summing the energy over all directions of the wavevector gives the energy spectrum $E(K)$
we have used since §12.1 (the book's letter for it is $S(K)$), with $\bar e=\int_0^\infty E(K)\,dK$ and, in the inertial range,
$E(K)=C\bar\varepsilon^{2/3}K^{-5/3}$ with the
*Kolmogorov constant* $C\approx1.5$ (a measured, public value: Sreenivasan 1995, *Phys. Fluids* 7, 2778). Isotropy links it to
the one-dimensional spectrum: the one-sided constant is $C_1=\tfrac{18}{55}C$, the book's two-sided one $\tfrac{9}{55}C$.""")
code(r"""
print(ch12.kolmogorov_constants(1.5))                         # C = 1.5 → one-sided C_1 = 18C/55 = 0.491, two-sided 9C/55 = 0.245
ratio = ch12.one_dimensional_from_3d(10.0, lambda K: 1.5*K**(-5/3))/(1.5*10.0**(-5/3))   # integrate the 3-D law down to 1-D at k_1 = 10
print(f"1-D / 3-D at the same wavenumber: {ratio:.5f} = 18/55 = {18/55:.5f};  S(K = 10) = {ch12.inertial_spectrum_3d(10.0, 1.0):.5f} m³/s²")
""", explain=r"""
`ch12.one_dimensional_from_3d(k1, E_fn)` does the isotropic integral that turns a three-dimensional spectrum into the
one-dimensional one a probe measures; for a −5/3 law the ratio is exactly $18/55=0.327$.""")
note("N98 [C]", r"""
**First confirmation.** The −5/3 law was first confirmed convincingly in a tidal channel, where the Reynolds number is
enormous and the inertial range spans decades; oceanic and atmospheric spectra have confirmed it ever since (the temperature
spectra of C15 are its scalar twin).""")
note("N99 [C]", r"""
**Why models can work at all.** Because the small scales are universal, a turbulence model (C12, C13) or a large-eddy
simulation (Ch. 10) may replace them by a formula that knows only $\bar\varepsilon$ and $\nu$.""")
P("P298", "quadrature on a logarithmic grid", r"""
A spectrum spans many decades, so we integrate it on a logarithmic grid. With $K=e^s$, $dK=K\,ds$, so
$\int E\,dK=\int E\,K\,ds$ with $s=\ln K$ evenly spaced.""", code=r"""
K_ = np.logspace(0, 6, 601)                       # wavenumbers from 1 to 10⁶, evenly spaced in ln K
approx = np.trapezoid(K_**(-5/3)*K_, np.log(K_))  # ∫ K^(−5/3) dK done as ∫ K^(−5/3) K d(ln K)
print(round(approx, 4), 1.5*(1 - 1e-4))           # exact: 1.5 (1 − 10⁻⁴)
""")
note("N97 [B]", r"""
**A model spectrum for the whole range.** Pao's (1965) form multiplies the inertial law by a viscous roll-off,
$E(K)=C\bar\varepsilon^{2/3}K^{-5/3}\exp\{-\tfrac32C(K\eta)^{4/3}\}$; with an outer scale $L$ we add an energy-range factor of
the kind used by Pope (2000) (`ch12.model_spectrum`). These formulas are from the literature, not from the book. Two integral
checks any spectrum must pass: $\int_0^\infty E\,dK=\bar e$ and $2\nu\int_0^\infty K^2E\,dK=\bar\varepsilon$.""")
nb.worked_example("one point on the −5/3 line", r"""
$\bar\varepsilon=1$ m²/s³, $k_1=100$ rad/m (a 6 cm eddy), two-sided constant $C_1=(9/55)\times1.5=0.245$.

1. $\bar\varepsilon^{2/3}=1$.
2. $k_1^{-5/3}=100^{-5/3}=4.64\times10^{-4}$.
3. $S_{11}=0.245\times4.64\times10^{-4}=1.14\times10^{-4}$ m³/s².
4. Ten times the wavenumber: $10^{-5/3}=1/46.4$ — the spectrum falls 46-fold per decade.
5. Units: $(\mathrm m^2/\mathrm s^3)^{2/3}(1/\mathrm m)^{-5/3}=\mathrm m^{4/3+5/3}\mathrm s^{-2}=\mathrm m^3/\mathrm s^2$ ✓.""")
code(r"""
K = np.logspace(0, 6, 400)                                    # wavenumbers [rad/m]
E = ch12.model_spectrum(K, 1.0, 1e-6, L=1.0, kind="pope")     # model E(K): ε̄ = 1 m²/s³, water, outer scale 1 m (the stirred tank of C07)
slope, const = ch12.fit_inertial_range(K, E, band=(50, 2000)) # least-squares straight line on log–log axes in the band 50 … 2000 rad/m
print(f"fitted slope in the inertial band: {slope:.3f} (−5/3 = {-5/3:.3f});  S_11(100) = {ch12.inertial_spectrum_1d(100.0, 1.0):.3e} m³/s²")
E_pao = ch12.model_spectrum(K, 1.0, 1e-6, kind="pao")         # Pao's form without an energy range
print(f"dissipation integral 2ν∫K²E dK = {2e-6*np.trapezoid(K**3*E_pao, np.log(K)):.4f} (should be ε̄ = 1);  energy ∫E dK = {np.trapezoid(E*K, np.log(K)):.3f} m²/s²")
k_eta, S_norm = ch12.kolmogorov_normalize_spectrum(K, E, 1e-6, 1.0)   # Eq. (12.53): (K η, E/(u_K² η))
print(f"in Kolmogorov units the spectrum ends near Kη = {k_eta[np.argmax(S_norm < 1e-3*S_norm[200])]:.2f}")
""", explain=r"""
1. `ch12.model_spectrum(K, eps, nu, L, kind)` returns $E(K)$; `ch12.fit_inertial_range` fits a power law in a band and returns
   (slope, prefactor): the slope is −1.68, i.e. −5/3 within 0.02 (the band still feels a little of both ends).
2. The dissipation integral of Pao's form returns $\bar\varepsilon$ to 4 digits — on the logarithmic grid of primer P298.
3. `ch12.kolmogorov_normalize_spectrum` puts the spectrum in the universal variables of $\frac{S_{11}(k_1)}{u_K^2\eta}=\Phi(k_1\eta)$ (12.53).""")
scratch(r"""
# From scratch: the slope of the inertial range by a straight-line fit on log–log axes
band = (K >= 50) & (K <= 2000)                                # the wavenumbers inside the band
slope_mine = np.polyfit(np.log(K[band]), np.log(E[band]), 1)[0]   # slope of ln E against ln K (P13)
assert abs(slope_mine + 5/3) < 0.02                           # −5/3 within 0.02
assert np.isclose(slope_mine, ch12.fit_inertial_range(K, E, (50, 2000))[0])   # same number as the library
print(f"slope by np.polyfit: {slope_mine:.4f}")
""", r"""
A power law is a straight line on log–log axes, so `np.polyfit` of $\ln E$ against $\ln K$ gives the exponent; the library does
exactly this.""")
note("N96 [B]", r"""
**Spectra in Kolmogorov scaling collapse at high wavenumber.** Plot $E/(u_K^2\eta)$ against $K\eta$: the dissipation end is one
curve for every flow and every Reynolds number, and only the low-wavenumber end (the energy-containing eddies) moves. The
slider figure shows our model spectra; measured spectra from very different flows do the same.""")
nb.plotly(r"""
Keta = np.logspace(-6, 0.5, 150)                             # the universal axis K η [–]

def f3(logRe):                                                # curves for one outer Reynolds number
    Re = 10.0**logRe
    nu_ = 1.0/Re                                              # L = 1 m, ΔU = 1 m/s → ε̄ = 1 m²/s³, ν = 1/Re
    eta_, uK_, _ = ch12.kolmogorov_scales(nu_, 1.0)           # Kolmogorov scales at this Reynolds number
    E_ = ch12.model_spectrum(Keta/eta_, 1.0, nu_, L=1.0, kind="pope")   # the model spectrum on the universal axis
    lo, hi = 2*np.pi*eta_*6.0, 0.1                            # a rough inertial range: well above 2π/L, well below 1/η
    return {"model spectrum E/(u_K² η)": (Keta, E_/(uK_**2*eta_)),
            "−5/3 law: C (Kη)^(−5/3)": (Keta, 1.5*Keta**(-5/3)),
            "inertial range (2π/L ≪ K ≪ 2π/η)": ([lo, hi] if lo < hi else [hi, hi], [1e-2, 1e-2])}

figF3 = slider_figure(f3, "log10 Re_L", np.linspace(3.0, 8.0, 11), xlabel="K η [–]", ylabel="E / (u_K² η) [–]",
                      title="In Kolmogorov units the small-scale end is one curve")
figF3.update_xaxes(type="log", range=[-6, 0.5])               # logarithmic K η axis
figF3.update_yaxes(type="log", range=[-3, 10])                # logarithmic spectrum axis
figF3.show()
""", explain=r"""
1. `f3(logRe)` builds the model spectrum for one Reynolds number and expresses it in Kolmogorov units.
2. The horizontal bar at the bottom marks a rough inertial range; it does not exist at $\mathrm{Re}_L=10^3$ and spans about four
   decades at $10^8$.""")
nb.figure_notes(
    see="Drag the Reynolds number up: the right-hand end of the teal curve (the roll-off near $K\\eta\\approx0.3$) does not move "
        "at all; the hump on the left slides further left along the dotted −5/3 line, and the bar marking the inertial range grows.",
    read="The right-hand end is the same curve for every Reynolds number — that is the universality of "
         "$\\frac{S_{11}(k_1)}{u_K^2\\eta}=\\Phi(k_1\\eta)$ (12.53). The inertial range is the straight stretch between the two ends; "
         "it needs $\\mathrm{Re}_L\\gtrsim10^5$ to span even one clean decade.",
    change="…we plotted against $KL$ instead (outer scaling): the *left* end would stay fixed and the roll-off would move — the "
           "view of animation A4. Both are the same spectra; the scaling decides which end looks universal.")
nb.live(r"""
def cascade_live(dU=1.0, L=1.0, log_nu=-6.0):                 # free velocity difference [m/s], outer scale [m] and log10 of viscosity [m²/s]
    nu_ = 10.0**log_nu
    eps_ = ch12.dissipation_outer_scaling(dU, L)              # ε̄ ~ ΔU³/L
    eta_, uK_, tau_ = ch12.kolmogorov_scales(nu_, eps_)       # Kolmogorov scales
    Kk = np.logspace(np.log10(0.3/L), np.log10(3/eta_), 300)  # from the outer scale to beyond η
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.loglog(Kk, ch12.model_spectrum(Kk, eps_, nu_, L=L, kind="pope"), color=C_FLUC)   # the model spectrum, dimensional
    ax.loglog(Kk, 1.5*eps_**(2/3)*Kk**(-5/3), ":", color=C_REF)                          # the −5/3 law
    ax.axvline(2*np.pi/L, color=C_MEAN); ax.axvline(1/eta_, color=C_VISC)                # the two ends
    ax.set(xlabel="K [rad/m]", ylabel="E(K) [m³/s²]", title=f"Re_L = {dU*L/nu_:.1e}, ε̄ = {eps_:.2g} m²/s³, η = {eta_*1e3:.3g} mm")
    plt.show()

live(cascade_live, dU=(0.1, 10.0, 0.1), L=(0.01, 100.0, 0.01), log_nu=(-6.0, -4.0, 0.1))   # sliders (kernel only)
""", explain=r"""
The same spectrum with three free sliders (it needs a running kernel; the slider figure above is its web-page twin). Changing
$\nu$ moves only the rose line at $1/\eta$; changing $\Delta U$ or $L$ moves the level and both ends.""")
explainer("energy_cascade_spectrum", "What does a higher Reynolds number change?",
          "dragging the Reynolds number slides η away from L, adds tiers and widens the −5/3 range while its level stays put; "
          "changing the viscosity at fixed ΔU and L moves the right end only.",
          ["Switch air ↔ water at the same ΔU and L: ε̄ does not move, η does.",
           "Load 'atmospheric boundary layer': 6.4 decades between L and η.",
           "Toggle the printed +5/3 ghost and read the units inspector.",
           "Double ΔU: ε̄ × 8, η × 2^(−3/4) = 0.59."])
whatif(r"""…the turbulence lives in a jet or a wake, where the big eddies and the mean flow change downstream? Conservation and
symmetry still give the answer (C09).""")


# =====================================================================================================================
# A.8 §12.8 Free Turbulent Shear Flows — C09
# =====================================================================================================================
nb.section("12.8", "Free Turbulent Shear Flows", intro=r"""
**What is this section about?** Jets, wakes, plumes and mixing layers — turbulence with no wall nearby. Far downstream they
forget their source and keep one shape that only stretches. One conserved quantity plus that shape-keeping gives how fast they
widen and slow down, without any turbulence model.

> ⚠️ **An assumption switches here.** From this section to §12.10 the density is constant ($\rho_0\to\rho$) and $P$ is the
> deviation from hydrostatic pressure; buoyancy returns in §12.11.""")
core("C09", r"The far field of the plane turbulent jet $U(x,y)=C_5(J_s/\rho)^{1/2}x^{-1/2}F(y/x)$ (12.66)",
     "How can anyone predict how a turbulent jet spreads without solving for the turbulence?")
problem(r"""
Smoke from a chimney, the exhaust of an engine, a river entering the sea, the wake of an island in the trade winds: each widens
downstream by swallowing the still fluid around it (**N101 [B]** *entrainment*), and each looks the same at every distance once
you rescale it (*self-preservation*). **N100 [B]**: a *jet* is fast fluid leaving a nozzle, a *wake* the slow fluid behind a
body, a *mixing (shear) layer* the zone between two streams, a *plume* a jet driven by buoyancy. A *free* shear flow has no
wall; §12.9 adds one.""")
fig(r"""
ILLUS = {"C_U": 2.5, "C_Y": 4.0, "xi_half_U": 0.10, "xi_half_Y": 0.15, "dxi80_coeff": 0.1}   # LABELLED ILLUSTRATIVE constants (ours — not the tabulated, measured values)
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.2))
for ax, flow in zip(axs, ("plane_jet", "plane_wake", "shear_layer")):
    free_shear_sketch(ax, flow=flow, constants=ILLUS)         # drawing: mean profiles at four stations, from ch12.free_shear_flow
    ax.set_title(flow.replace("_", " "), fontsize=10)
""",
    see="Three flows drawn from left to right: a jet (a bump of fast fluid that widens and weakens), a wake (a dip of slow "
        "fluid that widens and fills in), a mixing layer (a step between two streams that smears out). Each panel shows the "
        "mean profile at four stations and the growing edge of the turbulent region.",
    read="In each flow the profiles at different stations are the same shape, stretched in width and scaled in height. The "
         "amplitudes here use **labelled illustrative constants** — the picture shows the growth laws, not measured numbers.",
    change="…the surroundings were not still (a jet in a co-flow): the excess velocity would behave like a wake far downstream, "
           "with the wake's exponents.")
idea(r"""
conservation:   momentum flux  J_s = ρ ∫U² dy   is the same at every x        →  U_CL² · δ = const
shape-keeping:  U = U_CL(x) · F(y/δ(x))  must satisfy the equations at every x →  dδ/dx = const
together:       δ ∝ x,   U_CL ∝ x^(−1/2),   volume flux ∝ U_CL·δ ∝ x^(+1/2)  (it grows: entrainment)
""", r"""
Symbols: $x$ downstream distance from the slot [m], $y$ across the jet [m]; $U_{CL}(x)$ the centreline mean velocity [m/s];
$\delta(x)$ the width [m]; $\xi=y/\delta$ the similarity variable; $J_s$ the momentum flux per unit span [N/m = kg/s²].""")
remind("C09")
note("N102 [B]", r"""
**The similarity form of the mean velocity**: an amplitude that depends on $x$ times a shape $F$ that depends only on
$\xi=y/\delta(x)$; $F$ is even with $F(0)=1$.""",
     equation=r"U(x,y)=U_{CL}(x)F(y/\delta(x))", ref="12.56")
note("N103 [B]", r"""
**The same for the Reynolds shear stress**, with its own amplitude $\Psi(x)$ [m²/s²] and an odd shape $G$. ⚠️ These $F$, $G$
are jet profiles, not the tensor functions of $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40).""",
     equation=r"-\overline{uv}=\Psi(x)G(y/\delta(x))", ref="12.57")
note("N104 [B]", r"""**Mean continuity in two dimensions.**""", equation=r"\partial U/\partial x+\partial V/\partial y=0", ref="12.58")
note("N105 [B]", r"""**The mean $x$-momentum equation** (steady, constant density).""",
     equation=r"U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}=-\frac1\rho\frac{\partial P}{\partial x}+\nu\Big(\frac{\partial^2U}{\partial x^2}+\frac{\partial^2U}{\partial y^2}\Big)-\frac{\partial\overline{u^2}}{\partial x}-\frac{\partial\overline{uv}}{\partial y}", ref="12.59")
note("N106 [B]", r"""**The mean $y$-momentum equation.**""",
     equation=r"U\frac{\partial V}{\partial x}+V\frac{\partial V}{\partial y}=-\frac1\rho\frac{\partial P}{\partial y}+\nu\Big(\frac{\partial^2V}{\partial x^2}+\frac{\partial^2V}{\partial y^2}\Big)-\frac{\partial\overline{uv}}{\partial x}-\frac{\partial\overline{v^2}}{\partial y}", ref="12.60")
code(r"""
JET_KW = dict(C5="from_invariant", xi_half=0.10)              # amplitude fixed by the momentum invariant; ILLUSTRATIVE half-width ξ½ = 0.10
print("−mean(uv) at (x, y) = (1, 0.1):", round(ch12.plane_jet_reynolds_stress(1.0, 0.1, 1.0, 1.0, **JET_KW), 4), "m²/s²;",   # the Reynolds shear stress at three points
      "at (4, 0.4):", round(ch12.plane_jet_reynolds_stress(4.0, 0.4, 1.0, 1.0, **JET_KW), 4), "; on the axis:", ch12.plane_jet_reynolds_stress(1.0, 0.0, 1.0, 1.0, **JET_KW))   # four times farther at the same ξ, and on the axis
print("V at (1, 0.1):", round(ch12.plane_jet_cross_velocity(1.0, 0.1, 1.0, 1.0, **JET_KW), 4), "m/s; at the edge (1, 1.0):",   # the cross-stream mean velocity inside the jet …
      round(ch12.plane_jet_cross_velocity(1.0, 1.0, 1.0, 1.0, **JET_KW), 4), "m/s = −entrainment velocity", round(ch12.plane_jet_entrainment_velocity(1.0, 1.0, 1.0, **JET_KW), 4))   # … and at its edge: inflow
terms = ch12.thin_shear_layer_terms(1.0, 0.1, 1.0, 1.0, 1.5e-5, **JET_KW)   # sizes of the terms of (12.59) at x = 1 m, y = 0.1 m (air)
print({k: float(f"{float(v):.3g}") for k, v in terms.items()})
""", explain=r"""
1. Everything in this block uses $J_s/\rho=1$ m³/s² and an **illustrative** half-width $\xi_{1/2}=0.10$ (our choice, not a
   measured value); `C5="from_invariant"` lets the momentum invariant fix the amplitude.
2. `ch12.plane_jet_reynolds_stress` returns $-\overline{uv}$ (minus the correlation, the quantity of N103). Above the axis, where
   $dU/dy<0$, it is **negative** — so $\overline{uv}>0$ there: the stress has the sign of the mean shear, as an eddy viscosity
   would give ($-\overline{uv}=\nu_T\,dU/dy$). It is one quarter as large at four times the distance and the same $\xi$ (it decays
   like $x^{-1}$) and zero on the axis.
3. `ch12.plane_jet_cross_velocity`: $V$ is slightly outward inside the jet and **inward** at the edge — the surroundings flow in.
4. `ch12.thin_shear_layer_terms` compares the sizes of the terms of the $x$-momentum equation: the two advection terms balance the
   stress gradient exactly (`residual` 0), the viscous terms are a thousand times smaller, and $V/U\approx0.02$. These are the facts
   the next derivation uses.""")
NOTES_IN["D13"] = (r"**N107 [B]** the thin-layer equations $U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y}$, "
                   r"$0\cong-\frac1\rho\frac{\partial}{\partial y}(P+\rho\overline{v^2})$ (12.61) (steps 2 and 5) · **N108 [B]** the conservative "
                   r"(flux) form of the momentum equation (steps 6–7) · **N109 [B]** the momentum-flux invariant "
                   r"$J_s\equiv\rho_s\int_{-\infty}^{+\infty}[U^2]_{x=0}dy\cong\rho\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}$ (12.62) (steps 8–10)")
D("D13", ref="12.62", after=r"""
> ⚠️ **slip #18 — the book prints** that the boundary terms $[VU+\overline{uv}]_{y=-\infty}^{y=+\infty}$ vanish because $U$, $V$
> and $\overline{uv}$ all go to zero as $y\to\pm\infty$ **; the correct form is** that they vanish because $U\to0$ and
> $\overline{uv}\to0$. $V$ does **not** go to zero: by continuity $\partial U/\partial x+\partial V/\partial y=0$ (12.58) it tends to
> $V(\pm\infty)=\mp v_e$, the entrainment velocity with which the surroundings flow in (step 9). The product $VU$ still vanishes,
> so the invariant $J_s=\rho\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}$ (12.62) is unchanged. The cell above the derivation
> shows it in numbers: far outside the jet `ch12.plane_jet_cross_velocity` returns $-0.137$ m/s, exactly minus
> `ch12.plane_jet_entrainment_velocity`, not 0.""")
NOTES_IN["D14"] = (r"**N110 [B]** the cross-stream velocity eliminated with continuity (steps 1–2 and 7–9) · **N111 [B]** the similarity equation, "
                   r"the Result of this derivation (step 14)")
D("D14", ref="12.63", check_src=PF["D14"]["check_src"] + r"""
print("D14 check — advection computed directly minus the similarity form:", sp.simplify(lhs.subs(y, xi*delta) - rhs))   # 0 for every x and ξ
print("D14 check — the first two brackets for δ ~ x, U_CL ~ x^(−1/2):", sp.simplify(A1), sp.simplify(A2))               # −1/2 and +1/2: constants
""")
code(r"""
sim = ch12.plane_jet_similarity_sympy()                       # sympy repeats D14 with general δ(x), U_CL(x), Ψ(x) (cached)
print("check:", sim["check"])
print("the three brackets:", sim["coefficients"])
print("power family δ ~ x^m, U_CL ~ x^n:", sim["power_family"], "| momentum-flux exponent:", sim["momentum_flux_exponent"])
print("exponential family δ ~ e^(ax), U_CL ~ e^(−ax):", sim["exponential_family"])
print("G(ξ = 0.1) for the Gaussian F with ξ½ = 0.10:", round(float(ch12.plane_jet_stress_profile(0.1, xi_half=0.10)), 5))
""", explain=r"""
1. `ch12.plane_jet_similarity_sympy()` substitutes the similarity forms into the thin-layer equation with **general** functions
   $\delta(x)$, $U_{CL}(x)$, $\Psi(x)$ and returns the three brackets of D14's Result.
2. For powers $\delta\sim x^m$, $U_{CL}\sim x^n$ the brackets are $(n,\ m+n,\ 1)\,x^{m-1}$: constant only if $m=1$; and the momentum
   flux varies like $x^{m+2n}$: constant only if $n=-\tfrac12$.
3. `ch12.plane_jet_stress_profile(xi, xi_half)` integrates the similarity equation once for the stress shape,
   $C_3G=-\tfrac12F\int_0^\xi F\,d\xi$ (D15 step 8): $-0.02025$ at $\xi=0.1$.""")
NOTES_IN["D15"] = (r"**N112 [B]** the three brackets must be constants, $\frac{\delta U'_{CL}}{U_{CL}}=C_1$, $\frac{\delta U'_{CL}}{U_{CL}}+\delta'=C_2$, "
                   r"$\frac\Psi{U_{CL}^2}=C_3$ (12.64) (step 1) · **N113 [B]** linear growth from a virtual origin, $\delta=(C_2-C_1)(x-x_o)$, with "
                   r"$x_o$ of the order of the slot width (steps 2–3) · **N114 [B]** the invariant fixes the decay, "
                   r"$J_s=\rho U_{CL}^2\delta\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi=\rho C_4^2x^{2\gamma+1}\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi$ (12.65) "
                   r"⇒ $2\gamma+1=0$ (steps 5–6) · **N115 [B]** the stress, $-\overline{uv}=C_3U_{CL}^2G(y/x)=C_3C_5^2(J_s/\rho)x^{-1}G(y/x)$ (12.67) "
                   r"(step 8) · **N116 [B]** the volume flux grows, $\dot V(x)=\int_{-\infty}^{+\infty}U\,dy=C_5(J_s/\rho)^{1/2}x^{+1/2}\int_{-\infty}^{+\infty}F(\xi)\,d\xi$ (12.68) "
                   r"(step 9; its derivative, the entrainment velocity, in step 10)")
D("D15", ref="12.66", after=r"""
> ⚠️ **slip #4 — the book prints**, after the stress law, that the constants are "$C_3$ and $C_4$" and that "$C_5$, $C_6$" are
> tabulated **; the correct form is** $C_3$ and $C_5$ (with $C_5=C_4(\rho/J_s)^{1/2}$), and the tabulated pair is the one that
> multiplies the nozzle forms $U=C_5U_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}F(y/x)$ (12.72) and
> $\bar Y=C_7Y_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}H(y/x)$ (12.73).""")
note("N122 [B]", r"""
**A Gaussian is a good fit to the measured shape.** With the half-width $\xi_{1/2}$ (where $F=\tfrac12$),
$F(y/x)=\exp\{-\ln(2)(y/x)^2/(\xi_{1/2})_U^2\}$. The Gaussian integral (Ch. 8, P194) gives
$\int F\,d\xi=\xi_{1/2}\sqrt{\pi/\ln2}$ and $\int F^2d\xi=\xi_{1/2}\sqrt{\pi/(2\ln2)}$; with such an $F$ the invariant fixes the
amplitude, $C_5=\big(\int F^2d\xi\big)^{-1/2}$ (`ch12.gaussian_profile`, `ch12.profile_integrals`).""")
nb.worked_example("a plane jet with an illustrative half-width", r"""
Take $J_s/\rho=1$ m³/s² and an **illustrative** half-width $\xi_{1/2}=0.10$ (our choice, not the book's value).

1. $\int F^2d\xi=0.10\times\sqrt{\pi/(2\ln2)}=0.1505$, so $C_5=0.1505^{-1/2}=2.577$.
2. Centreline speed at $x=1$ m: $U_{CL}=2.577$ m/s; at $x=4$ m: $2.577/\sqrt4=1.289$ m/s.
3. Half-width: 0.10 m at $x=1$ m, 0.40 m at $x=4$ m.
4. Momentum flux at both stations: $\rho U_{CL}^2x\int F^2d\xi=\rho\times1.000$ ✓ (the same).
5. Volume flux: $\int F\,d\xi=0.2129$; at $x=1$ m $\dot V=2.577\times0.2129=0.549$ m²/s; at $x=4$ m it has doubled to
   1.097 m²/s — the extra fluid was entrained.
6. Entrainment velocity at $x=1$ m: $\tfrac12\,d\dot V/dx=\dot V/(4x)=0.137$ m/s on each side.""")
code(r"""
y_j = np.linspace(-3.0, 3.0, 6001)                            # cross-stream grid [m], wide enough for the profile to vanish
stations = (1.0, 2.0, 4.0)                                    # three downstream stations [m]
ints = ch12.profile_integrals(0.10, 0.15)                     # ∫F dξ, ∫F² dξ and ∫HF dξ for half-widths 0.10 (velocity) and 0.15 (scalar)
print("profile integrals:", {k: round(float(v), 4) for k, v in ints.items()}, "→ C5 =", round(ints["I2"]**-0.5, 3))
for x in stations:
    U = ch12.plane_jet_mean_velocity(x, y_j, 1.0, 1.0, **JET_KW)          # Eq. (12.66): U(x, y) for J_s = 1 N/m, ρ = 1 kg/m³
    J = ch12.jet_momentum_flux_per_span(y_j, U, 1.0)                      # Eq. (12.62): ρ ∫U² dy
    Vdot = ch12.plane_jet_volume_flux(x, 1.0, 1.0, **JET_KW)              # Eq. (12.68): ∫U dy
    ve = ch12.plane_jet_entrainment_velocity(x, 1.0, 1.0, **JET_KW)       # inflow speed at each edge
    print(f"x = {x:.0f} m: U_CL = {U.max():.3f} m/s, J = {J:.4f} N/m, volume flux = {Vdot:.3f} m²/s, entrainment velocity = {ve:.3f} m/s")
""", explain=r"""
1. `ch12.plane_jet_mean_velocity(x, y, Js, rho, C5, xi_half)` is the far-field law of this block's title with a Gaussian $F$.
2. `ch12.jet_momentum_flux_per_span` integrates $\rho U^2$ across the jet: 1.000 N/m at every station — the invariant.
3. The volume flux grows like $\sqrt x$ (0.549, 0.776, 1.097 m²/s); the entrainment velocity is what feeds it.""")
scratch(r"""
# From scratch: the two fluxes by the trapezoid rule at the three stations
J_mine, V_mine = [], []
for x in stations:
    U = 2.5774*x**-0.5*np.exp(-np.log(2)*(y_j/x)**2/0.10**2)  # U = C5 x^(−1/2) F(y/x) typed out (C5 from step 1 of the tiny example)
    J_mine.append(np.trapezoid(1.0*U**2, y_j))                # ρ ∫U² dy
    V_mine.append(np.trapezoid(U, y_j))                       # ∫U dy
assert np.allclose(J_mine, 1.0, rtol=1e-4)                    # the momentum flux is the same at every station
assert np.isclose(V_mine[2]/V_mine[0], 2.0, rtol=1e-6)        # the volume flux doubles from x = 1 m to x = 4 m
assert np.allclose(V_mine, [ch12.plane_jet_volume_flux(x, 1.0, 1.0, **JET_KW) for x in stations], rtol=1e-4)   # same as the library
print("momentum flux:", np.round(J_mine, 4), "volume flux:", np.round(V_mine, 4))
""", r"""
The profile is typed out from the formula, the two integrals done with `np.trapezoid`: the momentum flux is flat, the volume
flux doubles over a factor 4 in distance, and both agree with the library.""")
nb.plotly(r"""
y_f = np.linspace(-1.6, 1.6, 161)                            # cross-stream positions [m]
nu_lam = 0.02                                                 # a (large, illustrative) viscosity for the laminar comparison jet [m²/s]

def f4(x):                                                    # raw profiles at station x: turbulent and laminar jet of the same momentum flux
    Ut = ch12.plane_jet_mean_velocity(x, y_f, 1.0, 1.0, **JET_KW)         # turbulent: δ ∝ x, U_CL ∝ x^(−1/2)
    Ul = JET.free_jet(x, y_f, 1.0, 1.0, nu_lam)["u"]                      # laminar (Ch. 9): δ ∝ x^(2/3), U_CL ∝ x^(−1/3)
    return {"turbulent jet U(x, y)": (y_f, Ut), "laminar jet of Ch. 9 (same J_s)": (y_f, Ul),
            "turbulent, rescaled to the station x = 1 m: U × √x against y/x": (y_f/x, Ut*np.sqrt(x))}

figF4 = slider_figure(f4, "x", np.round(np.linspace(0.5, 8.0, 16), 2), unit="m", xlabel="y [m]   (rescaled curve: y/x × 1 m)",
                      ylabel="U [m/s]", title="Raw profiles spread and decay; the rescaled one stays", xrange=(-1.6, 1.6), yrange=(0, 4))
figF4.show()
""", explain=r"""
1. `f4(x)` returns the raw turbulent profile, the laminar jet of Ch. 9 (`JET.free_jet`, Bickley's $\mathrm{sech}^2$ solution) with
   the same momentum flux, and the turbulent profile **rescaled**: its height multiplied by $\sqrt x$ (undoing the decay
   $U_{CL}\propto x^{-1/2}$) and its width divided by $x$ (undoing the spreading). Rescaled like this, the profile of any station
   becomes the profile at $x=1$ m.
2. The rescaled curve is the same at every slider position — and coincides with the raw profile when the slider is at 1 m: that
   is self-preservation.""")
nb.figure_notes(
    see="Drag $x$ downstream: the purple turbulent profile widens in proportion to $x$ and its peak falls like $x^{-1/2}$; the "
        "laminar profile widens more slowly and decays more slowly; the third (rescaled) curve does not move at all.",
    read="One curve $F(y/x)$ describes the jet at every station. The area under $U^2$ is the same at every $x$ (the invariant); "
         "the area under $U$ grows (entrainment).",
    change="…the jet were laminar: it widens more slowly at first and its exponents are 2/3 and −1/3, because viscosity — a "
           "constant — sets the spreading, not eddies that grow with the jet.")
note("N117 [B]", r"""
**The scalar field carried by the jet** (the mass fraction $\bar Y$ of nozzle fluid) is self-similar too, with its own shape
$H$ (wider than $F$: scalars spread a little faster than momentum).""",
     equation=r"\bar Y(x,y)=Y_{CL}(x)H(y/x)", ref="12.69")
note("N118 [B]", r"""
**Its invariant is the flux of nozzle fluid.**

$$\dot M_s=\rho_s\int_{-\infty}^{+\infty}[U]_{x=0}dy\cong\rho\int_{-\infty}^{+\infty}\bar YU\,dy=\rho Y_{CL}C_5(J_s/\rho)^{1/2}x^{+1/2}\int_{-\infty}^{+\infty}H(\xi)F(\xi)\,d\xi\qquad\text{(12.70)}$$

> ⚠️ **slip #3 — the book prints** the source integral as $\rho_s\int_{-\infty}^{+\infty}[U]_{y=0}\,dy$ **; the correct form is**
> $\rho_s\int_{-\infty}^{+\infty}[U]_{x=0}\,dy$ — taken across the slot, at $x=0$, exactly as in
> $J_s\equiv\rho_s\int_{-\infty}^{+\infty}[U^2]_{x=0}dy$ (12.62).""")
note("N119 [B]", r"""**So the centreline concentration falls like $x^{-1/2}$**, the same exponent as the velocity.""",
     equation=r"\bar Y=C_6\big(\dot M_s/\sqrt{\rho J_s}\big)x^{-1/2}H(y/x)", ref="12.71")
note("N120 [B]", r"""**In nozzle variables** (slot width $d$, exit speed $U_0$, exit density $\rho_s$).""",
     equation=r"U=C_5U_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}F(y/x)", ref="12.72")
note("N121 [B]", r"""**And the scalar in nozzle variables** ($Y_0$ the exit mass fraction).""",
     equation=r"\bar Y=C_7Y_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}H(y/x)", ref="12.73")
code(r"""
Y1 = ch12.plane_jet_mass_fraction(1.0, 0.0, 1.0, 1.0, 1.0, C6=1.0, xi_half_Y=0.15)   # Eq. (12.71) on the axis at x = 1 m; C6 = 1 ILLUSTRATIVE
Y4 = ch12.plane_jet_mass_fraction(4.0, 0.0, 1.0, 1.0, 1.0, C6=1.0, xi_half_Y=0.15)   # the same at x = 4 m
Yinv = ch12.plane_jet_mass_fraction(1.0, 0.0, 1.0, 1.0, 1.0, C6="from_invariant", xi_half_Y=0.15, **JET_KW)   # C6 fixed by the flux invariant (12.70)
print(f"centreline mass fraction: {Y1:.3f} at 1 m, {Y4:.3f} at 4 m (x^(−1/2));  with C6 from the invariant: {Yinv:.3f}")
for x in stations:                                            # the scalar flux at three stations
    Yp = ch12.plane_jet_mass_fraction(x, y_j, 1.0, 1.0, 1.0, C6="from_invariant", xi_half_Y=0.15, **JET_KW)   # Ȳ(x, y) for Ṁ_s = 1 kg/(m s)
    Up = ch12.plane_jet_mean_velocity(x, y_j, 1.0, 1.0, **JET_KW)
    print(f"x = {x:.0f} m: ρ∫ȲU dy = {ch12.scalar_flux_per_span(y_j, Up, Yp, 1.0):.6f} kg/(m s)")
print(f"slot of d = 1 cm, U0 = 10 m/s, air: J_s = {ch12.slot_momentum_flux(1.2, 10.0, 0.01):.2f} N/m, Ṁ_s = {ch12.slot_mass_flux(1.2, 10.0, 0.01):.2f} kg/(m s)")
cl = [ch12.free_shear_centerline("plane_jet", x, constants=ILLUS, d=0.01, U0=10.0, rho_s=1.2, rho=1.2)["U_CL"] for x in (0.4, 1.6)]   # Eq. (12.72), illustrative C_U
print(f"centreline speed at 0.4 m / at 1.6 m = {cl[0]/cl[1]:.3f} (x^(−1/2) → 2, whatever the illustrative constant)")
""", explain=r"""
1. `ch12.plane_jet_mass_fraction` is the scalar law: with an illustrative $C_6=1$ it gives 1.0 and 0.5 — a factor 2 over a factor 4
   in distance. With `C6="from_invariant"` the constant follows from the flux invariant and our two illustrative half-widths.
2. `ch12.scalar_flux_per_span` integrates $\rho\bar YU$ across the jet: the same at every station (1 to six digits).
3. `ch12.slot_momentum_flux` and `ch12.slot_mass_flux` give the two invariants of a real slot ($\rho_sU_0^2d$ and $\rho_sU_0d$);
   `ch12.free_shear_centerline` evaluates the nozzle form with constants **we** supply (`ILLUS`): the ratio of two centreline speeds
   depends only on the exponent.""")
NOTES_IN["D16"] = (r"**N124 [B]** the exponents of every free shear flow follow from one invariant and one growth law; the local Reynolds "
                   r"number $U_s\delta/\nu$ varies like $x^{n+m}$ (steps 1–12)")
D("D16", ref="")
note("N123 [B]", r"""
**Our own table of exponents.** The exponents are mathematics and public; the book's amplitude constants and half-widths are
measured data and are **not reproduced here** — `ch12.FREE_SHEAR_CONSTANTS` is an empty dictionary, and every constant in this
notebook is passed explicitly and labelled illustrative.""")
code(r"""
rows = [ch12.free_shear_exponents(f) for f in ch12.FREE_SHEAR_FLOWS]      # exact fractions for the seven flows
tab = pd.DataFrame(rows)[["flow", "width", "velocity", "scalar", "reynolds", "invariant"]]   # width m, velocity n, scalar, local Reynolds number n + m
print(tab.to_string(index=False))
assert [r["reynolds"] for r in rows] == [ch12.local_reynolds_number_exponent(f) for f in ch12.FREE_SHEAR_FLOWS]   # n + m, flow by flow
print("public table of constants adopted:", ch12.FREE_SHEAR_CONSTANTS)   # {} — none
""", explain=r"""
1. `ch12.free_shear_exponents(flow)` solves the two linear equations of D16 in exact fractions: width $\propto x^m$, velocity
   scale $\propto x^n$, scalar, and the local Reynolds number $\propto x^{n+m}$.
2. Read the last column of numbers: a plane jet's Reynolds number *grows* downstream ($+\tfrac12$), a round jet's stays constant
   (0), a round wake's *falls* ($-\tfrac13$) — far enough downstream a round wake stops being turbulent.""")
note("N125 [B]", r"""
**A more general similarity.** The three brackets need not be constants; it is enough that they share one $x$-dependence:

$$\frac{\delta U'_{CL}}{U_{CL}}=C_8\Big(\frac{\delta U'_{CL}}{U_{CL}}+\delta'\Big)=C_9\frac\Psi{U_{CL}^2}\qquad\text{(12.74)}$$

for instance $\delta\sim x^m$, $U_{CL}\sim x^n$, $\Psi\sim x^{2n+m-1}$. So the constants of a self-similar flow need not be
universal — they may remember the nozzle.

> ⚠️ **slip #5 — the book prints**, as a second example of this general similarity, the exponential family
> $\delta\sim e^{ax}$, $U_{CL}\sim e^{-ax}$, $\Psi\sim e^{-ax}$ **; the correct form is** that with these exponents the middle
> bracket $\delta U'_{CL}/U_{CL}+\delta'$ is identically zero, so the three brackets are *not* proportional to one another and the
> condition above is not met; and $U_{CL}^2\delta\sim e^{-ax}$ is not constant, which violates the jet's invariant
> $J_s=\rho\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}$ (12.62). An exponential family that *does* satisfy both is
> $\delta\sim e^{ax}$, $U_{CL}\sim e^{-ax/2}$, $\Psi=$ const. (The book's next sentence itself remarks that the invariant rules out
> some of the possibilities.) The cell below checks all three families with sympy.""")
code(r"""
pw = ch12.general_similarity_check(lambda x: 0.1*x, lambda x: 2.577*x**-0.5, lambda x: 0.3/x, 2.0)   # power family m = 1, n = −1/2 at x = 2 m
ex = ch12.general_similarity_check(lambda x: 0.1*np.exp(0.5*x), lambda x: 2.0*np.exp(-0.5*x), lambda x: 0.3*np.exp(-0.5*x), 2.0)   # exponential family, a = 0.5
for name, r_ in (("power", pw), ("exponential", ex)):
    print(f"{name:12s}: c1 = {r_['c1']:+.4f}, c2 = {r_['c2']:+.4f}, c3 = {r_['c3']:+.4f}, d ln(U_CL² δ)/d ln x = {r_['momentum_flux_exponent']:+.3f}")
xs_, a_s = sp.symbols("x a", positive=True)                   # downstream distance and the growth rate of the exponential families
for label, d_s, U_s, Psi_s in (("printed: δ ~ e^(ax), U_CL ~ e^(−ax), Ψ ~ e^(−ax)", sp.exp(a_s*xs_), sp.exp(-a_s*xs_), sp.exp(-a_s*xs_)),
                               ("ours:    δ ~ e^(ax), U_CL ~ e^(−ax/2), Ψ = const", sp.exp(a_s*xs_), sp.exp(-a_s*xs_/2), sp.Integer(1))):
    c1_s = sp.simplify(d_s*sp.diff(U_s, xs_)/U_s)             # first bracket  δ U_CL'/U_CL
    c2_s = sp.simplify(c1_s + sp.diff(d_s, xs_))              # second bracket δ U_CL'/U_CL + δ'
    c3_s = sp.simplify(Psi_s/U_s**2)                          # third bracket  Ψ/U_CL²
    print(f"{label}: brackets ({c1_s}, {c2_s}, {c3_s});  U_CL² δ = {sp.simplify(U_s**2*d_s)}")
slope_vo, x0_vo = ch12.virtual_origin_fit(np.array([1.0, 2.0, 4.0]), 0.10*(np.array([1.0, 2.0, 4.0]) - 0.05))   # half-widths that start from x_o = 0.05 m
print(f"virtual-origin fit: spreading rate {slope_vo:.3f}, x_o = {x0_vo:.3f} m")
""", explain=r"""
1. `ch12.general_similarity_check(delta_fn, UCL_fn, Psi_fn, x)` evaluates the three brackets for trial functions by finite
   differences. Power family: $(n,\ m+n,\ 1)\times x^{m-1}\times$ amplitudes, momentum-flux exponent 0 ✓. Exponential family with
   the printed exponents: the middle bracket is 0 and the momentum flux falls with $x$ ✗ (slip #5).
2. The sympy lines do the same symbolically. Printed family: brackets $(-ae^{ax},\ 0,\ e^{ax})$ and $U_{CL}^2\delta=e^{-ax}$ — not
   allowed. Our family: $(-\tfrac a2e^{ax},\ +\tfrac a2e^{ax},\ e^{ax})$, all proportional to $e^{ax}$, and $U_{CL}^2\delta=1$ — allowed.
   So exponentials are not excluded as such; the printed exponents are.
3. `ch12.virtual_origin_fit` recovers the spreading rate and the virtual origin $x_o$ of **N113** from half-widths at three stations.""")
note("N126 [B]", r"""
**How far until a fuel jet is lean enough to burn?** (our own gas, nozzle and speed). A round jet of hydrogen leaves a nozzle of
diameter $d=2$ mm at about 100 m/s into air. Hydrogen burns completely when its *mass fraction* in the mixture is the
stoichiometric value — the reaction $\mathrm H_2+\tfrac12\mathrm O_2\to\mathrm H_2\mathrm O$ needs half a mole of O₂ per mole of
fuel, and air is 21 % O₂ by mole. (A *mole fraction* counts molecules, a *mass fraction* kilograms; they differ because the molecular
masses do.) On the axis of a round jet the mass fraction falls like $1/x$ (the table above), so there is one distance where it
reaches that value.""")
code(r"""
st = ch12.stoichiometric_mass_fraction(2.016, 28.97, 0.5, 0.21)            # H2 (2.016 g/mol) in air (28.97 g/mol): ½ O2 per fuel, air is 21 % O2
x_st = ch12.round_jet_distance_for_mass_fraction(st["Y_fuel"], 0.002, 0.083, 1.2, 1.0, C_Y=4.0)   # d = 2 mm, ρ_s = 0.083 (H2), ρ = 1.2 kg/m³ (air); C_Y = 4 ILLUSTRATIVE (not a tabulated value)
print(f"stoichiometric: mole fraction {st['v_fuel']:.3f}, mass fraction {st['Y_fuel']:.4f};  reached {x_st*100:.1f} cm = {x_st/0.002:.0f} diameters downstream (illustrative C_Y)")
""", explain=r"""
`ch12.stoichiometric_mass_fraction(M_fuel, M_air, moles of O₂ per mole of fuel, O₂ mole fraction of air)`: one mole of fuel needs
$0.5/0.21=2.38$ moles of air, so the fuel is $1/3.38=29.6$ % of the mixture by mole — but, hydrogen being so light, under 3 % by
mass. `ch12.round_jet_distance_for_mass_fraction` solves
$Y_{CL}=C_YY_0(\rho_s/\rho)^{1/2}(x/d)^{-1}$ for $x$. The amplitude constant `C_Y=4.0` is **illustrative** (not a tabulated value), so the distance
the cell prints (a few tens of nozzle diameters) shows the method, not a design value.""")
note("N127 [C]", r"""
**Turbulence quantities across the jet** (qualitative). The turbulent energy $\bar e$ is spread over the whole width with a
shallow maximum off the axis; the shear stress $\overline{uv}$ is zero on the axis by symmetry and largest near the steepest
mean shear, $\xi\approx\pm\xi_{1/2}$.""")
note("N128 [B]", r"""
**The turbulent-energy budget of the jet** (the budget of C06 for a thin steady layer, written with the triple correlation
$\tfrac12\overline{u_i^2u_j}$ of the exact budget):

$$0=-U\frac{\partial\bar e}{\partial x}-V\frac{\partial\bar e}{\partial y}-\overline{uv}\frac{\partial U}{\partial y}-\frac{\partial}{\partial y}\Big(\frac1{\rho_0}\overline{pv}+\frac12\overline{u_i^2v}\Big)-\bar\varepsilon\qquad\text{(12.75)}$$

advection by the mean flow, production, turbulent transport, dissipation. (⚠️ slip #10 again: the page prints the triple
correlation as $\tfrac12\overline{ev}$.)""")
note("N129 [C]", r"""
**Three balances across the jet.** On the axis there is no production (no shear): advection from upstream balances
dissipation. In mid-layer production ≈ dissipation. At the edge neither is large: turbulent transport carries energy outward
and advection (the entrained flow) carries it back.""")
fig(r"""
xi_b = np.linspace(-0.3, 0.3, 241)                            # similarity variable ξ = y/x across the jet
bj = ch12.jet_tke_budget(np.abs(xi_b), 0.10)                  # a MODEL budget (label: qualitative), terms in units of U_CL³/x
G_prof = np.array([float(ch12.plane_jet_stress_profile(abs(s), xi_half=0.10))*np.sign(s) for s in xi_b])   # the stress shape G(ξ), odd in ξ
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 3.5))
a1.plot(xi_b, ch12.gaussian_profile(xi_b, 0.10), color=C_MEAN, label="mean velocity $F=U/U_{CL}$")
a1.plot(xi_b, bj["e"]/bj["e"].max(), color=C_FLUC, label="turbulent energy $\\bar e$ (model, scaled)")
a1.plot(xi_b, G_prof/np.abs(G_prof).max(), color=C_RS, label="shear-stress shape $G$ (scaled)")
a1.set(xlabel="ξ = y/x [–]", ylabel="scaled profiles [–]", title="Across the jet")
a1.legend(fontsize=7)
for key, col in (("production", C_RS), ("dissipation", C_VISC), ("advection", C_MEAN), ("transport", C_REF)):
    a2.plot(xi_b, bj[key], color=col, label=key)
a2.set(xlabel="ξ = y/x [–]", ylabel="budget terms [$U_{CL}^3/x$]", title=f"Energy budget — {bj['label']}")
a2.legend(fontsize=7)
""",
    see="Left: the bell-shaped mean velocity, a broader turbulent-energy curve with a shallow off-axis maximum, and an "
        "S-shaped stress curve that is zero on the axis and extreme near $\\xi\\approx\\pm0.1$. Right: four budget curves — "
        "production (two humps where the shear is), dissipation (negative everywhere), advection and transport.",
    read="On the axis production is zero and advection feeds the dissipation; under the humps production ≈ dissipation; at the "
         "edges transport and advection trade places. ⚠️ Only the production curve follows from the similarity solution. The "
         "energy and dissipation shapes are **assumed** (`ch12.jet_tke_budget` is a qualitative model) and transport is the "
         "residual — read the shapes, not the numbers.",
    change="…this were a wake instead of a jet: production would again vanish on the axis, but advection would be relatively "
           "larger everywhere, because a wake's turbulence decays faster than it is produced.")
note("N130 [B]", r"""
**The simplest closure (ours, a forward pointer to C12).** Assume a constant eddy viscosity across the jet,
$\nu_T\propto U_{CL}\delta$. The similarity equation then becomes the laminar jet's with $\nu\to\nu_T$, so
$F=\mathrm{sech}^2(a\xi)$ — the laminar *shape* with the turbulent *exponents*.""")
code(r"""
ev = ch12.plane_jet_eddy_viscosity_profile(xi_b, 0.10)        # F = sech²(a ξ) with the same half-width, and the ν_T it implies
diff = np.abs(ev["F"] - ch12.gaussian_profile(xi_b, 0.10))    # distance from the Gaussian fit
print(f"a = {ev['a']:.2f}, ν_T/(U_CL x) = {ev['nu_hat']:.4f};  largest |sech² − Gaussian| inside the half-width: {diff[np.abs(xi_b) <= 0.10].max():.3f} of U_CL, anywhere: {diff.max():.3f}")
""", explain=r"""
`ch12.plane_jet_eddy_viscosity_profile` returns the $\mathrm{sech}^2$ profile with our half-width and the eddy viscosity that goes
with it, $\nu_T\approx0.003\,U_{CL}x$ — hundreds of times the molecular value for a laboratory jet. Inside the half-width it
differs from the Gaussian fit by a few per cent of $U_{CL}$; its tails are fatter.""")
explainer("turbulent_jet_similarity", "How does a jet spread without a turbulence model?",
          "toggling raw ↔ rescaled makes five profiles collapse onto one F(ξ); a wrong-exponent slider breaks either the collapse "
          "or the invariant bar; switching the flow changes the invariant and the exponents follow.",
          ["Toggle 'rescaled': five curves become one.",
           "Set the trial decay exponent to −0.4: the momentum bar grows as x^0.2 — not allowed.",
           "Switch to 'round jet': the velocity now falls as 1/x and the local Reynolds number stays constant.",
           "Switch to 'plane wake' and watch the width grow as √x."])
whatif(r"""…there is a wall? Then a second length — a tiny viscous one — enters, the flow is never Reynolds-number independent,
and the argument must be made twice (C10, C11).""")


# =====================================================================================================================
# A.9 §12.9 Wall-Bounded Turbulent Shear Flows — R03 R04, C10, C11
# =====================================================================================================================
nb.section("12.9", "Wall-Bounded Turbulent Shear Flows", intro=r"""
**What is this section about?** Turbulence next to a wall: pipes, channels, the boundary layer on a wing, the wind over the
ground. Two lengths matter — a viscous one of microns to millimetres at the wall and the thickness of the whole layer — and
where the two descriptions overlap the mean velocity is a logarithm.""")
core("C10", r"The law of the wall $U^+\equiv U/u_*=f(yu_*/\nu)=f(y^+)$ (12.80) with $u_*^2\equiv\tau_0/\rho$ (12.81)",
     "Why do profiles from different flows and Reynolds numbers fall on one curve near a wall?")
problem(r"""
Measure the wind 1 mm above a smooth plate in a wind tunnel, the water speed 0.1 mm from a pipe wall, the flow just above a
ship's hull. Plotted in m/s against mm they have nothing in common. Divide each speed by a "friction velocity" and each
distance by a "viscous length", both built from the wall stress alone, and they fall on one curve. **N131 [B]**: a wall brings
two length scales, a viscous one $l_\nu$ and the layer thickness $\delta$ with $l_\nu\ll\delta$, and — unlike a jet — the flow
never becomes independent of the Reynolds number on a smooth wall.""")
idea("", r"""
| Layer | What matters | Where | Velocity law |
|---|---|---|---|
| viscous sublayer | $\nu$, $\tau_0$ | $y^+<5$ | $U^+=y^+$ |
| buffer layer | viscous and Reynolds stress both | $5<y^+<30$ | neither simple law |
| logarithmic (overlap) layer | neither length | $y^+>30$ and $y/\delta<0.15$ | a logarithm (C11) |
| wake (outer) layer | $\delta$ | $y/\delta>0.15$ | defect law (C11) |

Symbols: $y$ distance from the wall [m]; $\tau_0$ wall shear stress [Pa]; $u_*=\sqrt{\tau_0/\rho}$ the friction velocity [m/s];
$l_\nu=\nu/u_*$ the viscous length [m]; a superscript $+$ means "in wall units": $y^+=y/l_\nu$, $U^+=U/u_*$.""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.6, 3.4))
wall_layers_sketch(ax, Re_tau=2000.0, Pi=0.2)                 # drawing: the four layers on a semi-log profile (illustrative κ, B, Π)
""",
    see="A mean-velocity profile in wall units on a logarithmic distance axis, with four bands: sublayer, buffer layer, "
        "logarithmic layer, wake.",
    read="The first three bands are the *inner layer*, where the wall stress and viscosity set the scales; the last two the "
         "*outer layer*, where the thickness $\\delta$ does. The logarithmic band belongs to both — the overlap that C11 exploits.",
    change="…the Reynolds number $\\delta^+=\\delta u_*/\\nu$ were smaller: the logarithmic band would shrink and vanish near "
           "$\\delta^+\\approx200$, where the buffer layer meets the wake.")
remind("C10")
note("N132 [B]", r"""
**A turbulent profile is blunt.** The figure compares the mean velocity of a turbulent channel flow (a composite model profile,
`WT.composite_profile`, with the illustrative pair $\kappa=0.41$, $B=5.0$ and no wake) with the laminar parabola of Ch. 8
(`LAM.channel_flow`) carrying the same flow rate. Three numbers appear in the model profiles of this block before C11 derives what
they mean: $\kappa$ (the von Kármán constant, the inverse slope of the logarithmic part of the profile) and $B$ (its intercept),
for which we use the labelled illustrative pair 0.41 and 5.0; and $\Pi$, the strength of the extra "wake" velocity near the outer
edge (0 = none). A *composite* profile is an inner formula plus that outer correction (primer P301 in C11).""")
fig(r"""
h_ch, nu_w = 0.1, 1.0e-6                                      # full channel height [m] and viscosity of water [m²/s]
u_st = 0.05                                                   # friction velocity [m/s] → δ⁺ = (h/2) u*/ν = 2500
y_ch = np.linspace(0.0, h_ch/2, 400)                          # from the wall to the centreline [m]
U_turb = WT.composite_profile(y_ch, h_ch/2, u_st, nu_w, kappa=KAPPA, B=B_LOG, Pi=0.0)   # composite inner law, no wake (Π = 0, illustrative)
U_bulk = np.trapezoid(U_turb, y_ch)/(h_ch/2)                  # bulk (cross-section mean) speed [m/s]
G_lam = 12*1e-3*U_bulk/h_ch**2                                # the pressure gradient −dp/dx that gives a laminar flow the same bulk speed (μ = 1e-3 Pa s)
U_lam = LAM.channel_flow(y_ch, h_ch, U=0.0, G=G_lam, mu=1e-3) # Ch. 8's parabola between fixed plates
if isinstance(U_lam, dict): U_lam = U_lam["u"]                # (the Ch. 8 function may return a dict of fields)
fig, ax = plt.subplots(figsize=(5.8, 3.6))
ax.plot(U_turb/U_bulk, y_ch/h_ch, color=C_MEAN, label="turbulent mean (composite model)")
ax.plot(np.asarray(U_lam)/U_bulk, y_ch/h_ch, "--", color=C_REF, label="laminar parabola, same flow rate")
ax.set(xlabel="U / U_bulk [–]", ylabel="y / h [–]", title="Blunter in the core, far steeper at the wall")
ax.legend(fontsize=8)
print(f"centreline/bulk: turbulent {U_turb[-1]/U_bulk:.2f}, laminar {np.asarray(U_lam)[-1]/U_bulk:.2f};  wall slope ratio turbulent/laminar ≈ {(U_turb[1]/y_ch[1])/(np.asarray(U_lam)[1]/y_ch[1]):.0f}")
""",
    see="Two profiles over the lower half of a channel. The dashed laminar parabola peaks at 1.5 times the bulk speed; the "
        "purple turbulent profile is nearly flat across the core (about 1.1 times the bulk speed) and drops to zero in a very "
        "thin layer at the wall.",
    read="Turbulent mixing evens out the momentum across the core, so the whole velocity change is squeezed against the wall: "
         "the wall shear — and the drag — is tens of times the laminar value at the same flow rate.",
    change="…the Reynolds number were higher: the core gets flatter and the wall layer thinner still.")
NOTES_IN["D17"] = (r"**N133 [B]** the mean momentum balance of a fully developed channel, $0=-\frac{\partial P}{\partial x}+\frac{\partial\bar\tau}{\partial y}$, "
                   r"$0=-\frac{\partial}{\partial y}(P+\rho\overline{v^2})$ with the total stress $\bar\tau=\mu\frac{\partial U}{\partial y}-\rho_0\overline{uv}$ (12.76) "
                   r"(steps 2–3) · **N134 [B]** $P(x,y)-P(x,0)=-\rho\overline{v^2}$ (step 4) · **N135 [B]** the pressure gradient is the same at "
                   r"every height, $\frac{\partial}{\partial x}P(x,y)-\frac{d}{dx}P(x,0)=-\rho\frac{\partial}{\partial x}\overline{v^2}=0$ (12.77) (step 5) · "
                   r"**N136 [B]** the total stress is a straight line, $\bar\tau(y)=\tau_0(1-2y/h)$ with $h$ the **full** height (steps 6–8)")
PF["D17"]["title"] = r"The linear total stress $\bar\tau=\tau_0(1-2y/h)$ and the pressure gradient it fixes"
D("D17", ref="12.76")
nb.recap("R03", "Channel: pressure gradient fixed by the wall stress", r"""
Ch. 8's force balance on a slug of fluid holds for the turbulent mean too: the pressure drop pushes on the cross-section, the
two walls hold back, so $dP/dx=-2\tau_0/h$ (12.90) — the last step of D17. (⚠️ slip #7: the book cites its Exercise 12.31 for
the proof; it is Exercise 12.32.)""", where="Ch. 8 §8.2")
nb.recap("R04", "Pipe: the same balance", r"""
For a pipe of diameter $d$: $dP/dx=-4\tau_0/d$ (12.91). This is Ch. 8's $\tau_0=\frac a2\frac{dp}{dz}$ (8.8) with radius
$a=d/2$ and $z\to x$. (There $\tau_0$ carries the sign of the pressure gradient; here $\tau_0>0$ is its magnitude, hence the
minus sign.)""", where="Ch. 8 §8.2")
nb.current_core = "C10"
code(r"""
print(f"channel, τ0 = 0.3 Pa, full height 0.1 m: dP/dx = {WT.channel_pressure_gradient(0.3, 0.1):.1f} Pa/m")   # Eq. (12.90): −2 τ0/h
print(f"pipe, τ0 = 0.3 Pa, diameter 0.1 m:       dP/dx = {WT.pipe_pressure_gradient(0.3, 0.1):.1f} Pa/m")     # Eq. (12.91): −4 τ0/d
tau_back = WT.wall_stress_from_pressure_gradient(-12.0, 0.1)  # the inverse: τ0 from the pipe's pressure gradient
assert np.isclose(tau_back, abs(LAM.pipe_wall_stress(0.05, dpdz=-12.0)))   # Ch. 8's (8.8) with a = d/2 gives the same magnitude
print(f"back to the wall stress: {tau_back:.2f} Pa;  total stress a quarter of the way up the channel: {WT.channel_total_stress(0.025, 0.1, 0.3):.2f} Pa")
st_ch = WT.total_stress(ch["y_over_delta"], ch["Uplus"], -ch["minus_uv_plus"], 1.0/ch["Re_tau"], 1.0)   # model channel of C06 in outer units: μ = 1/Re_τ, ρ = 1
print("viscous + Reynolds = total:", np.allclose(st_ch["viscous"] + st_ch["reynolds"], st_ch["total"]),
      "| total stress at y/δ = 0, 0.5, 1:", np.round(np.interp([0.0, 0.5, 1.0], ch["y_over_delta"], st_ch["total"]), 3))
""", explain=r"""
1. `WT.channel_pressure_gradient(tau0, h)` and `WT.pipe_pressure_gradient(tau0, d)` are the two recaps: −6 and −12 Pa/m.
2. `WT.wall_stress_from_pressure_gradient` inverts the pipe relation and agrees with Ch. 8's `LAM.pipe_wall_stress`.
3. `WT.channel_total_stress(y, h, tau0)` is **N136**: a quarter of the way up a channel of full height 0.1 m the total stress is
   half the wall stress, 0.15 Pa.
4. `WT.total_stress(y, U, uv, mu, rho)` returns a dictionary with the `viscous` part $\mu\,dU/dy$, the `reynolds` part
   $-\rho\overline{uv}$ and their sum `total`. On the model channel of C06 (lengths in $\delta$, velocities in $u_*$, so $\mu=1/\mathrm{Re}_\tau$)
   the total falls linearly from 1 at the wall to 0 on the centreline. Note the sign: we pass $\overline{uv}=-(\texttt{minus\_uv\_plus})$,
   which is negative.""")
note("N137 [B]", r"""
**Who carries the stress?** The total is a straight line, but its two parts trade places within a hair of the wall. The figure
uses `WT.stress_partition`, a **model** partition (Spalding's inner profile for the viscous part — approximate against DNS).""")
note("N138 [B]", r"""
**In a boundary layer** the stress is not exactly linear, because the mean flow accelerates or decelerates:
$U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}=-\frac1\rho\frac{\partial P}{\partial x}+\frac1\rho\frac{\partial\bar\tau}{\partial y}$ (12.78)
— the turbulent form of Ch. 9's laminar boundary-layer equation
$u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ (9.9) (there $u$, $v$ are the full velocities).
With no pressure gradient the left side vanishes at the wall, so the total stress is nearly **constant** across the inner layer.
The cell below puts a number on "nearly": `WT.boundary_layer_stress_from_profile` integrates this equation from the wall for
sampled model profiles of a boundary layer in air.""")
code(r"""
x_bl = np.linspace(2.9, 3.1, 5)                               # five stations around x = 3 m [m] (we need ∂U/∂x)
blx = WT.zpg_boundary_layer(x_bl, 60.0, 1.5e-5, kappa=KAPPA)  # thickness and friction velocity at each station (air, U_inf = 60 m/s; fits of N150 in C11)
us_mid, d_mid = float(blx["u_star"][2]), float(blx["delta99"][2])   # friction velocity [m/s] and thickness [m] at the middle station
y_bl = np.concatenate([[0.0], np.geomspace(0.5*1.5e-5/us_mid, 1.3*d_mid, 160)])   # heights from the wall (first point at y⁺ = 0.5) to beyond the edge [m]
U_bl = np.stack([WT.composite_profile(y_bl, float(blx["delta99"][i]), float(blx["u_star"][i]), 1.5e-5, kappa=KAPPA, B=B_LOG, Pi=0.2)
                 for i in range(5)], axis=1)                  # model profile at each station as one column: U[j, i] = U(y_j, x_i); Π = 0.2 illustrative
dUdx = np.gradient(U_bl, x_bl, axis=1, edge_order=2)          # ∂U/∂x by finite differences along the stations (np.gradient, Ch. 1 P22)
V_bl = np.zeros_like(U_bl)                                    # the wall-normal mean velocity, zero at the wall
V_bl[1:] = -np.cumsum(0.5*(dUdx[1:] + dUdx[:-1])*np.diff(y_bl)[:, None], axis=0)   # continuity: V = −∫ ∂U/∂x dy (running trapezoid sum with np.cumsum)
excess = WT.boundary_layer_stress_from_profile(x_bl, y_bl, U_bl, V_bl, 0.0, 1.2)   # τ̄(x, y) − τ0(x) from Eq. (12.78) with dP/dx = 0, ρ = 1.2 kg/m³ [Pa]
ratio_bl = 1.0 + excess[:, 2]/(1.2*us_mid**2)                 # total stress / wall stress at the middle station (τ0 = ρ u*²)
yp_bl = y_bl*us_mid/1.5e-5                                    # the heights in wall units
print(f"δ+ = {d_mid*us_mid/1.5e-5:.0f};  total stress / τ0 at y+ = 5: {np.interp(5.0, yp_bl, ratio_bl):.5f}, at y+ = 30: {np.interp(30.0, yp_bl, ratio_bl):.5f}, at y/δ = 0.1: {np.interp(0.1*d_mid, y_bl, ratio_bl):.4f}")
""", explain=r"""
1. `WT.zpg_boundary_layer` (the fitted formulas stated in C11) gives $\delta_{99}$ and $u_*$ at five stations; `WT.composite_profile`
   builds a model mean profile at each (inner law plus a wake; $\kappa=0.41$, $B=5.0$ "a common textbook pair", $\Pi=0.2$
   illustrative). `np.stack(..., axis=1)` puts the five profiles side by side as columns of one array.
2. `np.gradient(..., axis=1)` differentiates along the stations to get $\partial U/\partial x$; continuity then gives $V$ as a running
   integral in $y$ (`np.cumsum` of trapezoid pieces).
3. `WT.boundary_layer_stress_from_profile(x, y, U, V, dPdx, rho)` integrates the momentum equation from the wall and returns
   $\bar\tau-\tau_0$. In the inner layer the total stress stays **within 0.05 % of $\tau_0$ at $y^+=30$ and within about 2.4 % at
   $y/\delta=0.1$** — a *constant-stress layer*.
4. ⚠️ We print inner-layer values only: this model profile is not consistent with the outer fits it is combined with, so its
   stress in the outer layer and at the edge would mean nothing.""")
fig(r"""
yp_s = np.geomspace(0.3, 1000.0, 300)                         # y⁺ from the wall to the centreline of a channel with Re_τ = 1000
part = WT.stress_partition(yp_s, 1000.0, kappa=KAPPA, B=B_LOG)   # total, viscous and Reynolds parts in units of τ0 (a model)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.semilogx(yp_s, part["total"], color=COLORS["ink"], label="total $\\bar\\tau/\\tau_0=1-y^+/\\mathrm{Re}_\\tau$")
ax.semilogx(yp_s, part["viscous"], color=C_VISC, label="viscous part $\\mu\\,dU/dy$")
ax.semilogx(yp_s, part["reynolds"], color=C_RS, label="Reynolds part $-\\rho\\overline{uv}$")
ax.axvspan(5, 30, color=C_REF, alpha=0.12)
ax.set(xlabel="$y^+$ [–]", ylabel="stress / $\\tau_0$ [–]", title="Viscous stress carries the load only within y⁺ ≈ 30 of the wall")
ax.legend(fontsize=8)
i_eq = np.argmin(np.abs(part["viscous"] - part["reynolds"])[yp_s < 100])  # where the two parts are equal (searched inside y⁺ < 100)
print(f"viscous = Reynolds at y+ ≈ {yp_s[i_eq]:.1f};  at y+ = 30 the viscous share is {np.interp(30, yp_s, part['viscous'])/np.interp(30, yp_s, part['total']):.2f}")
""",
    see="A black total-stress curve that stays near 1 until $y^+\\approx100$ and then falls to zero at the centreline "
        "(a straight line in $y$, bent by the logarithmic axis); a rose viscous curve that starts at 1 and dies by "
        "$y^+\\approx30$; an orange Reynolds curve that takes over. The grey band is the buffer layer.",
    read="At the wall the fluctuations vanish, so viscosity carries everything; by $y^+\\approx10$ the two parts are equal — "
         "the place where C06 found the peak of production — and beyond $y^+\\approx30$ the Reynolds stress carries almost all.",
    change="…this were a zero-pressure-gradient boundary layer: the total would stay within a few per cent of $\\tau_0$ across "
           "the whole inner layer (a *constant-stress layer*), which is the assumption behind the mixing-length profile of C12.")
P("P300", "inner, outer and overlap: two descriptions that must agree where both hold", r"""
Close to the wall one set of variables works (*inner*), far away another (*outer*). If there is a region where both are valid,
the two formulas must give the same answer there — *matching*. It is the separation argument of Ch. 7 (a function of one
variable can equal a function of another only if both are constant) applied to two approximations.""", code=r"""
eps_ = 0.01                                       # a small ratio of the two lengths (inner length / outer length)
Y = 0.1                                           # a point in the middle: outer coordinate Y, inner coordinate y = Y/eps
exact = Y/(Y + eps_) + Y                          # a function with an inner part y/(1 + y) and an outer part Y
inner = 1.0                                       # inner description for large y:  y/(1 + y) → 1
outer = 1.0 + Y                                   # outer description for small Y: 1 + Y
print(round(exact, 3), inner, outer)              # in the overlap all three agree to within the small parameters
""")
note("N139 [B]", r"""
**Inner layer, outer layer, overlap.** In the *inner layer* viscosity matters and the right length is the viscous length
$l_\nu$; in the *outer layer* it does not and the right length is $\delta$; where both descriptions hold is the *overlap*
(`WT.layer_name(yplus, y_over_delta)` names the layer of a point with the nominal limits 5, 30 and 0.15).""")
note("N140 [B]", r"""
**The inner list.** Near the wall the mean velocity can depend only on the density, the wall stress, the viscosity and the
distance — the thickness $\delta$ and the free-stream speed are deliberately left out.""",
     equation=r"U=U(\rho,\tau_0,\nu,y)", ref="12.79")
NOTES_IN["D18"] = (r"**N141 [B]** the friction velocity $u_*^2\equiv\tau_0/\rho$ (12.81), the viscous length $l_\nu=\nu/u_*$ and the friction "
                   r"Reynolds number $\delta^+=\delta u_*/\nu$ (steps 3–4) · **N142 [B]** the viscous sublayer: $\mu(dU/dy)=\tau_0\Rightarrow U=\tau_0y/\mu$, "
                   r"or $U^+=y^+$ (12.82) (steps 5–7)")
D("D18", ref="12.80")
nb.worked_example("wall units for air over a plate", r"""
$\tau_0=0.3$ Pa, $\rho=1.2$ kg/m³, $\nu=1.5\times10^{-5}$ m²/s.

1. $u_*=\sqrt{0.3/1.2}=0.5$ m/s.
2. $l_\nu=\nu/u_*=3\times10^{-5}$ m = 30 µm.
3. The sublayer ($y^+<5$) is 0.15 mm thick.
4. At $y=0.09$ mm: $y^+=3$, $U^+=3$, so $U=1.5$ m/s.
5. A 3 cm boundary layer has $\delta^+=0.03/(3\times10^{-5})=1000$.""")
code(r"""
yp1, Up1, u_star, l_nu = WT.wall_units(9.0e-5, 1.5, 0.3, 1.2, 1.5e-5)   # y = 0.09 mm, U = 1.5 m/s, τ0 = 0.3 Pa in air
print(f"u* = {u_star:.2f} m/s, l_ν = {l_nu:.1e} m, y+ = {yp1:.1f}, U+ = {Up1:.1f};  sublayer law U+(3) = {WT.viscous_sublayer(3.0):.1f}")
print(f"δ+ of a 3 cm layer: {WT.friction_reynolds_number(0.5, 0.03, 1.5e-5):.0f};  back to SI: {WT.from_wall_units(3.0, 3.0, 0.5, 1.5e-5)}")
print("dimensionless groups of U = U(ρ, τ0, ν, y):", WT.law_of_the_wall_groups())   # exponents of each variable in the two groups
print([WT.layer_name(yp_, yd_) for yp_, yd_ in ((3.0, 0.003), (15.0, 0.015), (100.0, 0.1), (600.0, 0.6))])
""", explain=r"""
1. `WT.wall_units(y, U, tau0, rho, nu)` returns $(y^+, U^+, u_*, l_\nu)$ — the numbers of the tiny example; `WT.from_wall_units`
   goes back to metres and m/s.
2. `WT.law_of_the_wall_groups()` runs Ch. 1's Π theorem on the inner list: two groups, $U\rho^{1/2}\tau_0^{-1/2}=U/u_*$ and
   $y\rho^{-1/2}\tau_0^{1/2}\nu^{-1}=yu_*/\nu$ — exactly $U^+$ and $y^+$.
3. `WT.layer_name` names the four layers for four points of a $\delta^+=1000$ layer.""")
scratch(r"""
# From scratch: the four wall-unit numbers by hand
u_star_mine = np.sqrt(0.3/1.2)                    # Eq. (12.81): friction velocity [m/s]
l_nu_mine = 1.5e-5/u_star_mine                    # viscous length ν/u* [m]
yp_mine, Up_mine = 9.0e-5/l_nu_mine, 1.5/u_star_mine   # Eq. (12.80): y⁺ = y/l_ν and U⁺ = U/u*
assert np.allclose([yp_mine, Up_mine, u_star_mine, l_nu_mine], WT.wall_units(9.0e-5, 1.5, 0.3, 1.2, 1.5e-5))   # same numbers as the library
print(u_star_mine, l_nu_mine, yp_mine, Up_mine)
""", r"""
Four one-line formulas reproduce `WT.wall_units`: wall units are nothing more than dividing by $u_*$ and $l_\nu$.""")
P("P299", "semi-log axes: a logarithm is a straight line", r"""
On `ax.semilogx` (logarithmic horizontal axis, linear vertical axis) the law $U^+=\frac1\kappa\ln y^++B$ is a straight line; it
rises by $\ln(10)/\kappa=2.303/\kappa$ per decade (5.6 for $\kappa=0.41$). The linear law $U^+=y^+$ looks *curved* there.""", code=r"""
vals = np.log(np.array([10.0, 100.0, 1000.0]))/0.41 + 5.0   # U⁺ at y⁺ = 10, 100, 1000 for κ = 0.41, B = 5.0 (illustrative pair)
print(np.round(vals, 2), np.round(np.diff(vals), 2))        # each decade adds 2.303/0.41 = 5.62
""")
fig(r"""
yp_w = np.geomspace(0.3, 2000.0, 400)                         # y⁺ over almost four decades
U_sp = WT.spalding_uplus(yp_w, kappa=KAPPA, B=B_LOG)          # Spalding's single formula for the whole inner layer (a model; κ, B illustrative)
fig, ax = plt.subplots(figsize=(6.6, 3.8))
ax.semilogx(yp_w, U_sp, color=C_MEAN, lw=2.2, label="law of the wall $f(y^+)$ (Spalding's formula)")
ax.semilogx(yp_w[yp_w < 20], WT.viscous_sublayer(yp_w[yp_w < 20]), "--", color=C_VISC, label="sublayer $U^+=y^+$")
ax.semilogx(yp_w[yp_w > 5], WT.log_law(yp_w[yp_w > 5], kappa=KAPPA, B=B_LOG), "--", color=C_FLUC, label="logarithmic law (C11)")
for lo, hi, name in ((0.3, 5, "sublayer"), (5, 30, "buffer"), (30, 2000, "logarithmic")):
    ax.axvspan(lo, hi, color=C_REF, alpha=0.06 if name != "buffer" else 0.15)
    ax.text(np.sqrt(lo*hi), 23.5, name, ha="center", fontsize=8, color=C_REF)
ax.set(xlabel="$y^+=yu_*/\\nu$ [–]", ylabel="$U^+=U/u_*$ [–]", xlim=(0.3, 2000), ylim=(0, 25), title="One curve near every smooth wall")
ax.legend(fontsize=8, loc="center left")
err5 = abs(WT.spalding_uplus(5.0, kappa=KAPPA, B=B_LOG) - 5.0)/5.0   # how far the sublayer line is from the curve at y⁺ = 5
print(f"sublayer law at y+ = 5 is {100*err5:.1f} % above the curve; at y+ = 7: {100*abs(WT.spalding_uplus(7.0, kappa=KAPPA, B=B_LOG) - 7)/7:.1f} %")
""",
    see="A bold purple curve that follows the rose dashed line $U^+=y^+$ (curved on these axes) up to $y^+\\approx5$, bends "
        "through the buffer band, and joins the teal dashed straight line beyond $y^+\\approx30$.",
    read="The sublayer line is accurate to a few per cent up to $y^+\\approx5$ (7 % off by $y^+=7$); the straight line takes over after the buffer "
         "layer. In between neither simple law holds. (Spalding's formula is a smooth interpolation — a model that is within "
         "about half a unit of $U^+$ of channel DNS.)",
    change="…the wall were rough: the roughness elements poke through the sublayer, which disappears, and viscosity with it "
           "(C11's rough-wall law).")
whatif(r"""…we go far from the wall, where viscosity cannot matter? A second law holds there, and making the two agree gives
the logarithm (C11).""")

# ---------------------------------------------------------------------------------------------------------------------
core("C11", r"The logarithmic law $U^+=\frac1\kappa\ln(y^+)+B$ (12.88); rough wall $U^+=\frac1\kappa\ln(y/y_0)$ (12.93)",
     "Where does the logarithm come from?")
problem(r"""
A wind profile over open country, plotted against the logarithm of height, is a straight line — the fact behind every
wind-turbine siting study, every "wind at 10 m" in a weather report and every bulk formula for surface drag in a climate model.
No turbulence model is needed to get it: only the statement that two descriptions of the same flow must agree.""")
idea("", r"""
Near the wall the velocity depends on $y$ and the viscous length, not on $\delta$. Far away it depends on $y$ and $\delta$, not
on viscosity. In between, the *gradient* can depend on neither length. The only thing left with the right units is
$dU/dy=u_*/(\kappa y)$ — and the integral of $1/y$ is a logarithm. ⚠️ $\kappa$ here is the **von Kármán constant** (a pure
number near 0.4), not the thermal diffusivity of the mean heat equation in C04 (our $\kappa_{th}$).""")
remind("C11")
NOTES_IN["D19"] = (r"**N143 [B]** the outer list $U=U(\rho,\tau_0,\delta,y)$ (12.83) (step 1) · **N144 [B]** the velocity-defect law "
                   r"$\frac{U_\infty-U}{u_*}=F(y/\delta)=F(\xi)$ (12.84) (step 2) · **N145 [B]** the inner gradient $\frac{dU}{dy}=\frac{u_*^2}\nu\frac{df}{dy^+}$ (12.85) "
                   r"(step 3) · **N146 [B]** the outer gradient $-\frac{dU}{dy}=\frac{u_*}\delta\frac{dF}{d\xi}$ (12.86) (step 4) · **N147 [B]** the matching "
                   r"condition $-\xi\frac{dF}{d\xi}=y^+\frac{df}{dy^+}$ (12.87) (steps 5–6) · **N148 [B]** the outer logarithm "
                   r"$F(\xi)=-\frac1\kappa\ln(\xi)+A$ (12.89) (step 10)")
D("D19", ref="12.88")
nb.worked_example("reading a log law", r"""
Illustrative pair $\kappa=0.41$, $B=5.0$.

1. At $y^+=100$: $U^+=\ln(100)/0.41+5.0=11.23+5.0=16.23$.
2. At $y^+=1000$: 21.85 — one decade adds $2.303/0.41=5.62$.
3. Where do $U^+=y^+$ and the log law cross? Solve $y^+=\ln(y^+)/0.41+5.0$: $y^+\approx10.8$ — the middle of the buffer layer.
4. With $u_*=0.5$ m/s, $l_\nu=30$ µm: at $y=3$ mm ($y^+=100$) the mean speed is $16.23\times0.5=8.1$ m/s.""")
code(r"""
print(f"U+(100) = {WT.log_law(100.0, kappa=KAPPA, B=B_LOG):.3f};  sublayer and log law cross at y+ = {WT.log_law_crossing(kappa=KAPPA, B=B_LOG):.2f}")
print(f"friction law: U_inf+ at δ+ = 1000 with A = 1.0 (illustrative): {WT.friction_law_from_overlap(1000, kappa=KAPPA, A=1.0, B=B_LOG):.2f}")
print("outer groups of U = U(ρ, τ0, δ, y):", WT.defect_law_groups())
xi_d, defect = WT.velocity_defect(np.array([0.01, 0.05]), np.array([8.0, 9.0]), 10.0, 0.5, 0.1)   # two points of a 10 cm layer: U = 8, 9 m/s; U_inf = 10; u* = 0.5
print(f"defect variables: ξ = {xi_d}, (U_inf − U)/u* = {defect};  outer log law F(0.1) = {WT.log_law_defect(0.1, kappa=KAPPA, A=1.0):.3f}")
om = ch12.overlap_matching_sympy()                            # sympy repeats the matching of D19 (cached)
print("inner solution f =", om["f"], "| outer solution F =", om["F_outer"], "| checks:", om["checks"])
""", explain=r"""
1. `WT.log_law(yplus, kappa=, B=)` — both constants are **required keywords** everywhere in the module (no silent defaults);
   `WT.log_law_crossing` solves step 3 of the tiny example: 10.80.
2. `WT.friction_law_from_overlap` adds the inner and outer logarithms (D19 step 11): $U_\infty^+=\frac1\kappa\ln\delta^++A+B=22.85$;
   the outer constant $A=1.0$ is illustrative.
3. `WT.defect_law_groups()` is the Π theorem for the outer list: $U/u_*$ and $y/\delta$. `WT.velocity_defect` converts a profile to
   those variables; `WT.log_law_defect` is the outer logarithm $F(\xi)=-\frac1\kappa\ln(\xi)+A$ (12.89): 6.616 at $\xi=0.1$.
4. `ch12.overlap_matching_sympy()` solves the matching condition symbolically and confirms both logarithms and the friction law.""")
P("P301", "composite profile: inner law + outer correction", r"""
Add to an inner formula a correction that vanishes near the wall and grows to its full size at the outer edge: the sum is
usable across the whole layer. `WT.coles_wake(xi)` is such a correction shape: 0 at the wall, 1 at the edge.""", code=r"""
print(WT.coles_wake(np.array([0.0, 0.5, 1.0])))   # the cubic wake shape W = 3ξ² − 2ξ³ at the wall, mid-layer and edge: 0, 0.5, 1
""")
note("N153 [B]", r"""
**One formula for the whole inner layer** (Spalding 1961): it gives $y^+$ as a function of $U^+$,
$y^+=U^++e^{-\kappa B}\big[e^{\kappa U^+}-1-\kappa U^+-\tfrac12(\kappa U^+)^2-\tfrac16(\kappa U^+)^3\big]$, and tends to
$U^+=y^+$ at the wall and to the log law far from it. `WT.spalding_uplus` inverts it with `brentq` (Ch. 3, P108). It is a
**model** (approximate against DNS, within about 0.7 in $U^+$).""")
note("N154 [B]", r"""
**The outer profile** adds a *wake* to the logarithm (Coles 1956):
$U^+=\frac1\kappa\ln(y^+)+B+\frac{2\Pi}\kappa W(y/\delta)$ with $W=3\xi^2-2\xi^3$ or $\sin^2(\pi\xi/2)$. The wake strength $\Pi$
depends on the flow (it is larger in a boundary layer than in a pipe); `WT.composite_profile` takes it as a required keyword,
and our examples pass labelled illustrative values.""")
note("N155 [C]", r"""
**Open questions.** Whether the overlap is exactly logarithmic or a weak power law, and how layers with strong stress
gradients behave, are still debated. The sublayer is universal; the wake is not — hence slightly different constants are
reported for pipes, channels and boundary layers.""")
note("N156 [B]", r"""
**Which constants?** We print only cited public values (`WT.LOG_LAW_CONSTANTS`, each entry with its citation) and pass
$\kappa$ and $B$ explicitly everywhere.""")
note("N157 [B]", r"""
**An empirical link between the two constants** across flows and pressure gradients (Nagib & Chauhan 2008):""",
     equation=r"\kappa B=1.6[\exp(0.1663B)-1]", ref="12.92")
code(r"""
print("Spalding U+ at y+ = 1, 5, 12, 30, 100:", np.round(WT.spalding_uplus(np.array([1.0, 5.0, 12.0, 30.0, 100.0]), kappa=KAPPA, B=B_LOG), 3))
print(f"its inverse y+(U+ = 10) = {WT.spalding_yplus(10.0, kappa=KAPPA, B=B_LOG):.2f};  slope dU+/dy+ at U+ = 10: {WT.spalding_slope(10.0, kappa=KAPPA, B=B_LOG):.3f}")
for name, entry in WT.LOG_LAW_CONSTANTS.items():              # the cited public values (never unpacked with **: each entry also holds its citation)
    B_txt = entry['B'] if entry['B'] is not None else "not quoted by the source (we fit it below)"   # the DNS entry gives κ only
    print(f"{name}: kappa = {entry['kappa']}, B = {B_txt}  — {entry['citation'][:70]}…")
# PUBLIC-REPO RULE: keep κ = 0.40 here — at κ = 0.384 the printed pair coincides with a row of the book's table (tools/check_public.py fails).
B_nc = WT.nagib_chauhan_B(0.40)                               # Eq. (12.92) solved for B at an illustrative κ = 0.40
print(f"Eq. (12.92): κ = 0.400 → B = {B_nc:.2f} → back to κ = {WT.nagib_chauhan_kappa(B_nc):.3f}")
""", explain=r"""
1. Spalding's curve at five heights: 1.000 (sublayer), 4.87, 9.16 (buffer), 12.63, 16.08 (already close to the log law's 16.23).
2. `WT.LOG_LAW_CONSTANTS` holds the pairs we may show: a common textbook pair, and $\kappa=0.384$ from channel DNS; each entry
   carries a `citation` string, which is why we read the keys one by one.
3. `WT.nagib_chauhan_B` and `WT.nagib_chauhan_kappa` are the two directions of $\kappa B=1.6[\exp(0.1663B)-1]$ (12.92): a round
   trip returns the $\kappa$ we started from. The printed $B$ is **our evaluation** of that formula at an illustrative $\kappa=0.40$
   (our choice, between the textbook 0.41 and the DNS 0.384), not a number taken from a table.""")
note("N149 [B]", r"""
**The profile at many Reynolds numbers.** In wall units the inner part never moves; a higher $\delta^+$ only lengthens the
logarithmic stretch before the wake peels off. The slider figure shows a composite model profile; the dots are public channel
DNS (Lee & Moser 2015, *J. Fluid Mech.* 774, 395; file `reference/ch12/lee_moser_2015_channel_mean.csv`).""")
nb.plotly(r"""
dns = pd.read_csv(REF / "lee_moser_2015_channel_mean.csv", comment="#")   # public DNS mean profiles (case, Re_tau, y/δ, y⁺, U⁺, dU⁺/dy⁺)
dns_pts = {c: dns[dns["case"] == c].iloc[::8] for c in (180, 1000, 5200)}   # three cases, every 8th point

def f5(logRe):                                                # curves for one friction Reynolds number δ⁺ = 10^logRe
    Re_tau = 10.0**logRe
    yp_ = np.geomspace(0.5, Re_tau, 120)                     # from inside the sublayer to the edge of the layer
    out = {"composite model profile (κ = 0.41, B = 5.0, Π = 0.1; illustrative)": (yp_, WT.composite_profile_plus(yp_, Re_tau, kappa=KAPPA, B=B_LOG, Pi=0.1)),
           "log law": (yp_[yp_ > 8], WT.log_law(yp_[yp_ > 8], kappa=KAPPA, B=B_LOG)),
           "wake begins: y+ = 0.15 δ+": ([0.15*Re_tau]*2, [0.0, 32.0])}
    for c, d_ in dns_pts.items():                             # the DNS dots are the same at every slider position
        out[f"channel DNS, Re_τ ≈ {c} (Lee & Moser 2015)"] = (d_["yplus"].values[1:], d_["Uplus"].values[1:])
    return out

figF5 = slider_figure(f5, "log10 δ+", np.round(np.linspace(np.log10(180), 5.0, 10 if FAST else 14), 3), xlabel="y+ [–]", ylabel="U+ [–]",
                      title="The inner curve never moves; the logarithm lengthens",
                      modes={f"channel DNS, Re_τ ≈ {c} (Lee & Moser 2015)": "markers" for c in (180, 1000, 5200)})
figF5.update_xaxes(type="log", range=[np.log10(0.5), 5.0])    # logarithmic y⁺ axis
figF5.update_yaxes(range=[0, 32])                             # fixed U⁺ range
figF5.show()
for c in (1000, 5200):                                        # how far are the DNS dots from the model line at the SAME Reynolds number?
    d_ = dns[dns["case"] == c]
    Re_c = float(d_["Re_tau"].iloc[0])                        # the exact Re_τ of this DNS case
    at = np.array([30.0, 100.0, 300.0, 1000.0, 3000.0]); at = at[at < 0.6*Re_c]   # heights inside the layer
    gap = np.interp(at, d_["yplus"], d_["Uplus"]) - WT.composite_profile_plus(at, Re_c, kappa=KAPPA, B=B_LOG, Pi=0.1)   # DNS minus model
    print(f"Re_τ ≈ {c}: DNS − model at y+ = {at.astype(int).tolist()}: {np.round(gap, 2).tolist()}")
""", explain=r"""
1. `pd.read_csv(..., comment="#")` reads the public DNS table (lines starting with `#` hold the citation).
2. `f5(logRe)` returns the composite profile (`WT.composite_profile_plus`: Spalding's inner law plus a small illustrative wake), the
   log law, a vertical marker where the wake begins, and the DNS dots.
3. `modes=` draws the DNS as markers instead of lines.
4. The printed lines compare DNS and model at the same Reynolds number: the DNS is **above** the model by a few tenths of a wall
   unit everywhere in the log region.""")
nb.figure_notes(
    see="Drag $\\delta^+$ from 180 to 10⁵: the left part of the curve stays exactly where it is; its right end travels along "
        "the straight log line, and the vertical marker ($y^+=0.15\\,\\delta^+$) travels with it. The three sets of DNS dots lie "
        "along the same inner curve and leave it one after the other.",
    read="This is the law of the wall as data: three simulations at very different Reynolds numbers share one inner curve. "
         "The DNS dots sit slightly *above* the model line through the log region (by 0.3 to 0.8 in $U^+$, the numbers the "
         "cell prints). The gap is largest at $y^+=30$, where Spalding's curve is still low, smallest near $y^+\\approx300$, and "
         "grows again further out: the DNS line is **steeper** — a smaller $\\kappa$, about 0.384 rather than the textbook 0.41 — "
         "and the illustrative pair 0.41, 5.0 was not fitted to this simulation (next cell).",
    change="…we plotted in outer variables ($(U_\\infty-U)/u_*$ against $y/\\delta$): the wakes would collapse and the wall "
           "region would fan out — the opposite collapse, which the explainer below lets you switch on.")
scratch(r"""
# From scratch: κ and B from a straight-line fit of U⁺ against ln y⁺ — on public DNS, and on the model profile
d52 = dns[dns["case"] == 5200]                                # channel DNS at Re_τ ≈ 5200 (Lee & Moser 2015, J. Fluid Mech. 774, 395)
yp_d, Up_d, Re_d = d52["yplus"].values, d52["Uplus"].values, float(d52["Re_tau"].iloc[0])   # y⁺, U⁺ and the exact Re_τ of the case
fit = WT.fit_log_law(yp_d, Up_d, window=(350.0, 0.15), Re_tau=Re_d)       # the library: window 350 < y⁺ < 0.15 δ⁺
win = (yp_d >= 350.0) & (yp_d <= 0.15*Re_d)                               # the same window by hand
slope, intercept = np.polyfit(np.log(yp_d[win]), Up_d[win], 1)            # straight line U⁺ = slope · ln y⁺ + intercept (P13)
assert np.allclose((1/slope, intercept), (fit["kappa"], fit["B"]), rtol=1e-8)   # κ = 1/slope, B = intercept: same as the library
print(f"DNS Re_τ ≈ {Re_d:.0f}, window 350 < y+ < 0.15 δ+: κ = {fit['kappa']:.4f}, B = {fit['B']:.2f} ({fit['n_points']} points)")
assert abs(fit["kappa"] - 0.384) < 0.004                                  # the published value for this simulation
yp_m = np.geomspace(1.0, 5000.0, 600)                                     # a composite MODEL profile at δ⁺ = 5000, built with κ = 0.41, B = 5.0
Up_m = WT.composite_profile_plus(yp_m, 5000.0, kappa=KAPPA, B=B_LOG, Pi=0.0)
rows = [("DNS Re_τ ≈ 5200", 350.0, fit), ("DNS Re_τ ≈ 5200", 30.0, WT.fit_log_law(yp_d, Up_d, window=(30.0, 0.15), Re_tau=Re_d)),
        ("model built with 0.41, 5.0", 30.0, WT.fit_log_law(yp_m, Up_m, window=(30.0, 0.15), Re_tau=5000.0)),
        ("model built with 0.41, 5.0", 300.0, WT.fit_log_law(yp_m, Up_m, window=(300.0, 0.15), Re_tau=5000.0))]
print(pd.DataFrame([{"profile": n, "window starts at y+": lo, "fitted κ": round(f_["kappa"], 3), "fitted B": round(f_["B"], 2)} for n, lo, f_ in rows]).to_string(index=False))
""", r"""
**What does the code above do, and what does it show?**

1. `np.polyfit` of $U^+$ against $\ln y^+$ inside the window gives slope $=1/\kappa$ and intercept $=B$; `WT.fit_log_law` does the
   same and the `assert` shows the two agree.
2. On the public DNS at $\mathrm{Re}_\tau\approx5200$, in the window $350<y^+<0.15\,\delta^+$, the fit returns $\kappa\approx0.384$ —
   the value its authors report.
3. **The constants you fit depend on where the window starts** (the table). Start the DNS window at $y^+=30$ and $\kappa$ rises
   to about 0.40. Fit a model profile that was *built* with 0.41 and 5.0 from $y^+=30$ and you get back a smaller pair, because at
   $y^+=30$ the buffer layer has not quite joined the logarithm; start at 300 and the built-in pair nearly returns. This is one
   reason published values differ from one experiment to the next (N155, N156).""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.4, 3.6))
for Re_tau, col in ((180.0, C_REF), (1000.0, C_FLUC), (10000.0, C_MEAN)):
    yp_ = np.geomspace(1.0, Re_tau, 300)
    ind = WT.log_law_indicator(yp_, WT.composite_profile_plus(yp_, Re_tau, kappa=KAPPA, B=B_LOG, Pi=0.1))   # y⁺ dU⁺/dy⁺ of the model profile
    ax.semilogx(yp_, ind, color=col, label=f"model, δ+ = {Re_tau:.0f}")
ax.semilogx(d52["yplus"].values[1:], (d52["yplus"]*d52["dUdy_plus"]).values[1:], ".", ms=3, color=C_RS, label="DNS Re_τ ≈ 5200 (Lee & Moser 2015)")
ax.axhline(1/KAPPA, color=C_REF, ls=":"); ax.axhline(1/0.384, color=C_RS, ls=":")
ax.set(xlabel="y+ [–]", ylabel="$y^+dU^+/dy^+$ [–]", ylim=(0, 6), title="A logarithm shows as a plateau at 1/κ — only when δ+ is large")
ax.legend(fontsize=7)
""",
    see="The quantity $y^+dU^+/dy^+$ against $y^+$: every curve rises to a peak near $y^+\\approx10$ and falls. At "
        "$\\delta^+=180$ it never levels off; at $10^4$ the model shows a long plateau at $1/0.41=2.44$ (grey dotted). The DNS "
        "dots level off near $1/0.384=2.60$ (orange dotted) only for $y^+$ of several hundred.",
    read="This *indicator function* is the matching condition $-\\xi\\frac{dF}{d\\xi}=y^+\\frac{df}{dy^+}$ (12.87) itself: it "
         "equals $1/\\kappa$ exactly where a logarithm holds. It is the honest way to look for a log law — and it shows that even "
         "at $\\mathrm{Re}_\\tau\\approx5200$ the plateau is short.",
    change="…the profile had a stronger wake (a boundary layer in an adverse pressure gradient): the curve would turn upward "
           "sooner at the outer end, shortening the plateau.")
NOTES_IN["D20"] = (r"**N159 [B]** the rough-wall logarithm $U^+=\frac1\kappa\ln(y/y_0)$ (12.93): the roughness length $y_0$ is the height at "
                   r"which the logarithm extrapolates to zero velocity (steps 1–3)")
D("D20", ref="12.93")
note("N158 [B]", r"""
**Hydrodynamically smooth or rough?** If the roughness elements are buried in the viscous sublayer the wall is *smooth* and
the logarithm sits on a sublayer; if they poke through it the wall is *rough*, the sublayer is destroyed, viscosity drops out
and the logarithm extrapolates to $U=0$ at the roughness length $y_0$.""")
fig(r"""
y_r = np.geomspace(1e-4, 1.0, 300)                            # heights from 0.1 mm to 1 m
u_st_r, nu_air = 0.4, 1.5e-5                                  # friction velocity [m/s] and viscosity of air [m²/s]
U_smooth = u_st_r*WT.spalding_uplus(y_r*u_st_r/nu_air, kappa=KAPPA, B=B_LOG)     # smooth wall: inner law in SI units
U_rough = WT.rough_wall_log_law(y_r, u_st_r, 0.003, kappa=KAPPA)                 # rough wall with y0 = 3 mm (NaN below y0)
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.semilogx(y_r, U_smooth, color=C_MEAN, label="smooth wall: sublayer + logarithm")
ax.semilogx(y_r, U_rough, color=C_RS, label="rough wall, $y_0$ = 3 mm: $U=(u_*/\\kappa)\\ln(y/y_0)$")
ax.plot([0.003], [0], "o", color=C_RS)
ax.set(xlabel="height y [m]", ylabel="U [m/s]", ylim=(0, 12), title="Same friction velocity, same slope — a rough wall only shifts the line")
ax.legend(fontsize=8)
""",
    see="Two lines with the same slope on semi-log axes. The purple (smooth) one bends into a sublayer at the far left; the "
        "orange (rough) one is straight all the way down to $U=0$ at $y_0=3$ mm (the dot).",
    read="The slope is $u_*/\\kappa$ in both cases — it measures the stress. Roughness lowers the whole line: at the same height "
         "the wind over a rough surface is slower for the same stress, i.e. the drag coefficient is larger.",
    change="…$y_0$ were ten times larger (crops instead of short grass): the line would shift down by "
           "$(u_*/\\kappa)\\ln10\\approx2.2$ m/s at every height.")
nb.worked_example("wind over grass", r"""
$z_0=0.03$ m (we write $z$ for height in the atmosphere), $u_*=0.4$ m/s, $\kappa=0.41$.

1. $U(10\text{ m})=(0.4/0.41)\ln(10/0.03)=0.976\times5.81=5.67$ m/s.
2. Neutral drag coefficient at 10 m: $C_D=[\kappa/\ln(z/z_0)]^2=(0.41/5.81)^2=4.98\times10^{-3}$.
3. Check: $\tau_0=\rho u_*^2=1.2\times0.16=0.192$ Pa $=\rho C_DU^2=1.2\times4.98\times10^{-3}\times32.1$ ✓.

**Climate hook:** this $C_D$ is the neutral value every bulk surface-flux formula starts from; C15 corrects it for stability.""")
code(r"""
U10 = WT.rough_wall_log_law(10.0, 0.4, 0.03, kappa=KAPPA)     # Eq. (12.93): wind at 10 m over grass (z0 = 3 cm) for u* = 0.4 m/s
CD = WT.drag_coefficient_neutral(10.0, 0.03, kappa=KAPPA)     # C_D = [κ/ln(z/z0)]²
print(f"U(10 m) = {U10:.2f} m/s;  C_D = {CD:.2e};  u* back from the wind: {WT.friction_velocity_from_wind(U10, 10.0, 0.03, kappa=KAPPA):.3f} m/s")
print(f"stress two ways: ρ u*² = {1.2*0.4**2:.3f} Pa, ρ C_D U² = {1.2*CD*U10**2:.3f} Pa")
""", explain=r"""
`WT.rough_wall_log_law`, `WT.drag_coefficient_neutral` and `WT.friction_velocity_from_wind` are the three uses of the rough-wall
law: wind from stress, the bulk coefficient, and stress from one wind measurement — the numbers of the tiny example.""")
note("N150 [B]", r"""
**Zero-pressure-gradient boundary layer: fitted formulas** (Monkewitz, Chauhan & Nagib 2007, for $\mathrm{Re}_x>10^6$, as
quoted in the book): momentum thickness $\theta\approx0.016\,x\,\mathrm{Re}_x^{-0.15}$; displacement thickness
$\delta^*\approx\theta\exp\{7.11\kappa/\ln(\mathrm{Re}_\theta)\}$;
$\delta_{99}=0.2\,\delta^*[\kappa^{-1}\ln(\mathrm{Re}_{\delta^*})+3.30]$; and the skin-friction coefficient
$C_f\cong2.0/[\kappa^{-1}\ln(\mathrm{Re}_{\delta^*})+3.30]^2$ — the friction law of D19 again, since $C_f=2/(U_\infty^+)^2$
(`WT.zpg_boundary_layer`). The thicknesses $\theta$, $\delta^*$, $\delta_{99}$ are those of Ch. 9.""")
note("N151 [B]", r"""
**Two older skin-friction fits**: $C_f\cong0.370(\log_{10}\mathrm{Re}_x)^{-2.584}$ (Schultz-Grunow) and
$C_f\cong0.455/[\ln(0.06\,\mathrm{Re}_x)]^2$ (White). ⚠️ One uses $\log_{10}$, the other the natural logarithm
(`WT.skin_friction_zpg(Re_x, law=…)`).""")
note("N152 [B]", r"""
**A boundary layer with our own numbers**: air at $U_\infty=60$ m/s, 3 m from the leading edge, $\nu=1.5\times10^{-5}$ m²/s, so
$\mathrm{Re}_x=1.2\times10^7$ — compared with what a laminar (Blasius) layer would be at the same place.""")
code(r"""
bl = WT.zpg_boundary_layer(3.0, 60.0, 1.5e-5, kappa=0.384)    # the fits of N150 at x = 3 m, U_inf = 60 m/s in air (κ = 0.384: the cited DNS value)
print(f"Re_x = {bl['Re_x']:.1e}: θ = {bl['theta']*1e3:.2f} mm, δ* = {bl['delta_star']*1e3:.2f} mm, δ99 = {bl['delta99']*1e3:.1f} mm, H = {bl['H']:.2f}, C_f = {bl['Cf']:.2e}, u* = {bl['u_star']:.2f} m/s, δ+ = {bl['Re_tau']:.0f}")
print(f"laminar (Blasius) at the same Re_x: C_f = {BL.blasius_skin_friction(bl['Re_x']):.2e}, δ99 = {5.0*3.0/np.sqrt(bl['Re_x'])*1e3:.1f} mm")
for law in ("schultz_grunow", "white"):                       # the two older fits
    print(f"C_f by {law}: {WT.skin_friction_zpg(bl['Re_x'], law=law):.2e}")
""", explain=r"""
1. `WT.zpg_boundary_layer(x, U_inf, nu, kappa=)` returns a dictionary of thicknesses and friction: $\theta\approx4.2$ mm,
   $\delta_{99}\approx32$ mm, $C_f\approx2.3\times10^{-3}$; $H=\delta^*/\theta\approx1.3$ is the *shape factor* of Ch. 9 (2.59 for the
   laminar Blasius layer — a fuller profile has a smaller $H$).
2. A laminar layer at the same $\mathrm{Re}_x$ would be about 4 mm thick with a twelve times smaller $C_f$ (Blasius,
   $C_f=0.664/\sqrt{\mathrm{Re}_x}$, Ch. 9): turbulence makes the layer much thicker **and** the wall shear much larger.
3. The three skin-friction fits agree within a few per cent.""")
note("N161 [B]", r"""
**Laminar against turbulent friction over eight decades of Reynolds number.**""")
fig(r"""
Re_x = np.geomspace(1e4, 1e10, 200)                           # Reynolds numbers based on distance from the leading edge
fig, ax = plt.subplots(figsize=(6.2, 3.6))
ax.loglog(Re_x, BL.blasius_skin_friction(Re_x), "--", color=C_REF, label="laminar (Blasius) $0.664/\\sqrt{\\mathrm{Re}_x}$")
turb = Re_x > 1e6                                             # the fits are for Re_x > 10⁶
ax.loglog(Re_x[turb], [WT.skin_friction_zpg(r_, law="monkewitz", kappa=0.384) for r_ in Re_x[turb]], color=C_MEAN, label="turbulent: Monkewitz et al. fit")
ax.loglog(Re_x[turb], WT.skin_friction_zpg(Re_x[turb], law="white"), ":", color=C_RS, label="turbulent: White's fit")
ax.set(xlabel="$\\mathrm{Re}_x=U_\\infty x/\\nu$ [–]", ylabel="$C_f$ [–]", title="Turbulent skin friction falls slowly — like 1/ln², not like a power")
ax.legend(fontsize=8)
""",
    see="A dashed laminar line of slope −½ on log–log axes and, from $\\mathrm{Re}_x=10^6$, two nearly identical turbulent "
        "curves that lie well above it and fall much more slowly, bending gently.",
    read="The turbulent $C_f$ depends on the logarithm of the Reynolds number — the friction law that follows from the log law. "
         "It never becomes independent of Reynolds number on a smooth wall, but it changes very slowly.",
    change="…the wall were rough: beyond some Reynolds number the curve would level off completely, because the roughness "
           "length, not the viscous length, would set the inner scale.")
note("N160 [B]", r"""
**Pipe friction from the log law** (stated). Integrating the log law over a pipe's cross-section (one integration by parts,
Ch. 9 P218a) gives the bulk speed $U_{av}\cong u_*[(1/\kappa)\ln(au_*/\nu)+B-3/(2\kappa)]$ ($a$ the radius); with the Darcy
friction factor $f_D=8(u_*/U_{av})^2$ this becomes an implicit friction law of the form of Prandtl's
$f_D^{-1/2}=2.0\log_{10}(\mathrm{Re}_df_D^{1/2})-0.8$.""")
code(r"""
a_p, us_p = 0.05, 0.5                                         # pipe radius [m] and friction velocity [m/s] (air)
Uav = WT.pipe_bulk_velocity_loglaw(a_p, us_p, 1.5e-5, kappa=KAPPA, B=B_LOG)   # the closed form for the bulk speed
yq = np.geomspace(1e-9, a_p, 4000)                            # distance from the wall for a numerical check
Uq = us_p*WT.log_law(np.maximum(yq*us_p/1.5e-5, 1e-12), kappa=KAPPA, B=B_LOG)   # log law at every y (integrable at the wall)
Uav_num = 2/a_p**2*np.trapezoid(Uq*(a_p - yq), yq)            # U_av = (2/a²) ∫ U (a − y) dy
print(f"bulk speed: closed form {Uav:.3f} m/s, quadrature {Uav_num:.3f} m/s;  derived coefficient of log10: {np.log(10)/(KAPPA*np.sqrt(8)):.2f} (Prandtl: 2.0)")
fD = WT.pipe_friction_factor_turbulent(1e5)                   # Prandtl's law at Re_d = 10⁵ (constants as cited in the docstring)
print(f"f_D(Re_d = 1e5): turbulent {fD:.4f}, laminar formula 64/Re would give {LAM.pipe_friction_factor(1e5):.5f} — the turbulent value is {fD/LAM.pipe_friction_factor(1e5):.0f} times larger")
""", explain=r"""
1. `WT.pipe_bulk_velocity_loglaw` is the closed form; the trapezoid integral of the log law over the cross-section reproduces it.
2. Turning the natural logarithm into $\log_{10}$ gives the coefficient $\ln(10)/(\kappa\sqrt8)=1.99$ for $\kappa=0.41$ — Prandtl's 2.0.
3. `WT.pipe_friction_factor_turbulent(1e5)` solves the implicit law: $f_D\approx0.018$, about 28 times what the laminar formula
   $64/\mathrm{Re}$ of Ch. 8 would give at that Reynolds number.""")
explainer("law_of_the_wall", "Where does the logarithm come from?",
          "sliding the friction Reynolds number stretches the log region while the inner curve stays fixed in wall units; switching "
          "to outer units makes the opposite collapse; the indicator shows the plateau 1/κ appear only when δ⁺ is large.",
          ["Raise Re_τ from 180 to 5200: the straight part grows from nothing to more than a decade.",
           "Switch to outer scaling: now the wakes collapse and the wall region fans out.",
           "Click y⁺ = 12: buffer layer — neither law holds.",
           "Switch to a rough wall and change y₀⁺."])
whatif(r"""…we tried to *compute* this profile from the averaged equations? We must guess the Reynolds stress. The simplest guess
that returns the logarithm is next (C12).""")


# =====================================================================================================================
# A.10 §12.10 Turbulence Modeling — C12, C13
# =====================================================================================================================
nb.section("12.10", "Turbulence Modeling", intro=r"""
**What is this section about?** Closing the gap: replacing the unknown Reynolds stress by a formula in terms of the mean flow.
We build the simplest such model, see exactly what its one constant does, then meet the two-equation model used in most
engineering codes.""")
core("C12", r"The turbulent-viscosity hypothesis $\overline{u_iu_j}=\frac23\bar e\delta_{ij}-\nu_T\big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\big)$ (12.94) and the mixing length",
     "What does 'modelling the Reynolds stress' mean in practice, and how can one constant give the whole profile?")
problem(r"""
Every ocean and atmosphere model has a number called the eddy viscosity (Ch. 13's Ekman layer cannot be written down without
it). It is not a property of water or air — honey has a viscosity, a storm does not — but of the *flow*: big energetic eddies
mix momentum fast. **N162 [B]**: the *turbulent-viscosity hypothesis* for momentum and the *gradient-diffusion hypothesis* for
heat and scalars assume, by analogy with Newton's, Fourier's and Fick's laws of Ch. 1 and 4, that a turbulent flux is
proportional to the gradient of the mean quantity.""")
idea("", r"""
A diffusivity is (speed of the carriers) × (distance they travel before mixing). For molecules: thermal speed × mean free
path. For turbulence: eddy speed $u_T$ × eddy size $l_T$. Near a wall an eddy cannot be bigger than its distance from the wall:
$l_T=\kappa y$. One guessed length — and the log law comes out. Symbols: $\nu_T$ eddy (turbulent) viscosity [m²/s]; $l_T$ the
*mixing length* [m]; $u_T$ a turbulent velocity scale [m/s].""")
remind("C12")
code(r"""
gU_s = np.zeros((3, 3)); gU_s[0, 1] = 2.0                     # a simple shear dU/dy = 2 1/s
print(np.round(ch12.eddy_viscosity_stress(gU_s, 0.01, 0.3), 3))   # Eq. (12.94) with ν_T = 0.01 m²/s and ē = 0.3 m²/s²: mean(u_i u_j) [m²/s²]
print(f"heat flux by gradient diffusion: mean(wT') = {ch12.gradient_diffusion_flux(-0.01, 1.0):.3f} K m/s for dT/dz = −0.01 K/m, κ_T = 1 m²/s")   # Eq. (12.95)
print(f"ν_T ~ l_T u_T = {ch12.eddy_diffusivity_estimate(0.1, 0.5):.3f} m²/s for l_T = 0.1 m, u_T = 0.5 m/s")   # Eq. (12.98): the parcels of C04
""", explain=r"""
1. `ch12.eddy_viscosity_stress(gradU, nu_T, e)` evaluates the hypothesis of this block's title: the shear component is
   $\overline{uv}=-\nu_T\,dU/dy=-0.02$ m²/s² (negative for a positive shear, as in C04), and each normal stress is
   $\tfrac23\bar e=0.2$ m²/s² — the model makes the normal stresses equal, which real shear flows are not (C04's anisotropy).
2. `ch12.gradient_diffusion_flux(grad, K)` is "flux = −diffusivity × gradient"; `ch12.eddy_diffusivity_estimate(l_T, u_T)` the
   product of a length and a velocity: 0.05 m²/s for the parcels of C04's tiny example.""")
note("N163 [B]", r"""**Gradient diffusion of heat** with an eddy thermal diffusivity $\kappa_T$ [m²/s].""",
     equation=r"\overline{u_iT'}=-\kappa_T\,\partial\bar T/\partial x_i", ref="12.95")
note("N164 [B]", r"""**Gradient diffusion of a scalar** with an eddy mass diffusivity $\kappa_{mT}$ [m²/s].""",
     equation=r"\overline{u_iY'}=-\kappa_{mT}\,\partial\bar Y/\partial x_i", ref="12.96")
note("N165 [B]", r"""
**The closed mean momentum equation.** Put the hypothesis
$\overline{u_iu_j}=\frac23\bar e\delta_{ij}-\nu_T\big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\big)$ (12.94)
into the constant-density RANS equation:

$$\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-\frac1\rho\frac{\partial P}{\partial x_i}+\frac{\partial}{\partial x_j}\Big([\nu+\nu_T]\Big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\Big)-\frac23\bar e\delta_{ij}\Big)\qquad\text{(12.97)}$$

The eddy viscosity simply adds to the molecular one, and the $\tfrac23\bar e$ term acts as an extra pressure.

> ⚠️ **slip #6 — the book prints** the pressure gradient as $\partial P/\partial x_j$ **; the correct form is**
> $\partial P/\partial x_i$: $i$ is the free index of the equation, and an index that appears once in a term cannot be $j$.""")
code(r"""
nu_m, nuT_m, rho_m = 1.0e-3, 0.05, 1.0                        # molecular viscosity, a constant eddy viscosity [m²/s], density [kg/m³]
U_f = lambda x: np.array([x[1]**2, 0.0])                      # a manufactured mean flow: U = y², V = 0
P_f = lambda x: 2*(nu_m + nuT_m)*rho_m*x[0]                   # the pressure that balances it: dP/dx = 2 (ν + ν_T) ρ
res_ok = ch12.rans_eddy_viscosity_residual(U_f, P_f, lambda x: nuT_m, lambda x: 0.3, np.array([0.4, 0.7]), nu=nu_m, rho=rho_m)   # corrected index
res_bad = ch12.rans_eddy_viscosity_residual(U_f, P_f, lambda x: nuT_m, lambda x: 0.3, np.array([0.4, 0.7]), nu=nu_m, rho=rho_m, printed=True)   # as printed
print("residual of (12.97), corrected:", np.round(res_ok, 8), "| with the printed index:", np.round(res_bad, 4))
""", explain=r"""
`ch12.rans_eddy_viscosity_residual` evaluates left side minus right side of the closed equation by finite differences for fields
given as functions. Our manufactured field $U=y^2$ with $dP/dx=2(\nu+\nu_T)\rho$ satisfies it (residual ≈ 0 in both components).
With `printed=True` — one reading of the page's ill-formed index — a spurious pressure force appears in the $y$-component.""")
note("N166 [B]", r"""
**Why the analogy with molecules is imperfect.** Molecules travel a mean free path far smaller than the flow (small Knudsen
number, Ch. 1), so molecular viscosity is a property of the fluid. Eddies are as big as the shear layer itself,
$l_T/L=O(1)$ — there is no separation of scales, so $\nu_T$ cannot be a universal constant; it must be modelled flow by flow.""")
note("N167 [B]", r"""
**The size of an eddy diffusivity.** C16 will derive the eddy diffusivity $D_T=\overline{u^2}\Lambda_t$ for dispersion — the
same product of a velocity squared and a time, i.e. a velocity and a length.""",
     equation=r"\nu_T,\ \kappa_T,\ \text{or}\ \kappa_{mT}\sim l_Tu_T", ref="12.98")
NOTES_IN["D21"] = (r"**N168 [B]** the mean momentum equation of a unidirectional flow with an eddy viscosity, "
                   r"$0=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big(\nu\frac{\partial U}{\partial y}-\overline{uv}\big)=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big([\nu+\nu_T]\frac{\partial U}{\partial y}\big)$ (12.99) "
                   r"(the Start and step 1) · **N169 [B]** Prandtl's mixing-length stress $-\overline{uv}=l_T^2(dU/dy)^2$ (steps 2–3) · **N170 [B]** the wall model "
                   r"$0=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big(\nu\frac{dU}{dy}+\kappa^2y^2\big(\frac{dU}{dy}\big)^2\big)$ (12.100) (step 4) · "
                   r"**N171 [B]** its first integral and the exact slope (steps 5–8) · **N172 [B]** the logarithmic limit "
                   r"$\frac{dU}{dy}\cong\sqrt{\frac{\tau_0}\rho}\frac1{\kappa y}$, or $\frac U{u_*}\cong\frac1\kappa\ln y+\mathit{const.}$ (12.101) (step 10)")
D("D21", ref="12.101")
nb.worked_example("the model at y⁺ = 10", r"""
$\kappa=0.41$.

1. The first integral in wall units: $s+\kappa^2y^{+2}s^2=1$ with $s=dU^+/dy^+$.
2. $\kappa^2y^{+2}=0.1681\times100=16.81$.
3. Positive root: $s=2/(1+\sqrt{1+4\times16.81})=2/(1+8.261)=0.216$.
4. So the viscous stress carries 21.6 % and the Reynolds stress 78.4 % of $\tau_0$ here.
5. Compare the pure log slope $1/(\kappa y^+)=0.244$: already close.
6. Eddy viscosity $\nu_T/\nu=\kappa^2y^{+2}s=3.63$.""")
code(r"""
plain = ch12.mixing_length_wall_profile(10.0, KAPPA)                          # the model of D21 at y⁺ = 10, l_T = κ y
damped = ch12.mixing_length_wall_profile(10.0, KAPPA, damping="van_driest")   # the same with van Driest's wall damping (A⁺ = 26)
for name, r_ in (("l_T = κy", plain), ("with van Driest damping", damped)):   # one line of numbers per model
    print(f"{name:24s}: slope {r_['slope']:.4f}, U+ = {r_['Uplus']:.2f}, −mean(uv)+ = {r_['minus_uv_plus']:.3f}, ν_T/ν = {r_['nuT_over_nu']:.2f}")
print(f"intercept B of the model's log law: {ch12.mixing_length_intercept(KAPPA):.2f} without damping "
      f"(exact: [ln(4κ) − 1]/κ = {(np.log(4*KAPPA) - 1)/KAPPA:.2f}), {ch12.mixing_length_intercept(KAPPA, A_plus=26.0):.2f} with A+ = 26")   # the intercept with and without damping
""", explain=r"""
1. `ch12.mixing_length_wall_profile(yplus, kappa)` returns a dictionary: the slope $s$ (= the viscous share of the stress), $U^+$,
   `minus_uv_plus` $=-\overline{uv}/u_*^2$ (the Reynolds share, **minus** the correlation, positive) and $\nu_T/\nu$ — the numbers of the
   tiny example.
2. `ch12.mixing_length_intercept(kappa)` is the constant $B$ of the logarithm the model tends to. **The plain model has the right
   slope $1/\kappa$ but its line sits about 6 units too low** ($B=-1.23$ instead of ≈ 5): with $l_T=\kappa y$ the eddies mix right down to
   the wall and the sublayer is far too thin.
3. *Van Driest's damping* is an empirical fix: $l_T=\kappa y[1-e^{-y^+/A^+}]$ switches the eddies off within about $A^+$ wall units
   of the wall. With the customary $A^+=26$ the intercept becomes 5.28. This is what "tuning a closure" means.""")
scratch(r"""
# From scratch: the positive root of s + κ² y⁺² s² = 1 on a log grid, integrated with the running trapezoid rule
yp_g = np.concatenate([[0.0], np.geomspace(1e-3, 1000.0, 4000)])              # y⁺ from the wall outward
s_g = 2.0/(1.0 + np.sqrt(1.0 + 4.0*KAPPA**2*yp_g**2))                         # D21 step 8: the exact slope dU⁺/dy⁺
U_mine = cumulative_trapezoid(s_g, yp_g, initial=0.0)                         # U⁺ = ∫ s dy⁺ (running integral, Ch. 8 P190)
check_at = np.array([1.0, 10.0, 100.0, 1000.0])                               # four heights
U_lib = np.array([ch12.mixing_length_wall_profile(y_, KAPPA)["Uplus"] for y_ in check_at])   # the library (closed form)
assert np.allclose(np.interp(check_at, yp_g, U_mine), U_lib, rtol=1e-4)       # same profile
print("U+ at y+ = 1, 10, 100, 1000:", np.round(U_lib, 3), "| log law with B = −1.23:", np.round(np.log(check_at[2:])/KAPPA - 1.2324, 3))
""", r"""
The quadratic's positive root is evaluated point by point and integrated with `cumulative_trapezoid`; it matches the library's
closed form at four heights. At $y^+=100$ and 1000 the profile is already the logarithm with $B=-1.23$.""")
nb.plotly(r"""
yp_f = np.geomspace(0.5, 1000.0, 50 if FAST else 80)          # y⁺ on a logarithmic grid
kappas = (0.384, 0.40, 0.41)                                  # three values of the von Kármán constant

def f6(A_plus):                                               # model profiles for one damping constant A⁺
    out = {}
    for k_ in kappas:
        prof = ch12.mixing_length_wall_profile(yp_f, k_, damping=("van_driest" if A_plus > 0 else None), A_plus=max(A_plus, 1e-9))
        out[f"mixing-length model, κ = {k_}"] = (yp_f, prof["Uplus"])
    out["reference log line κ = 0.41, B = 5.0 (illustrative pair)"] = (yp_f[yp_f > 10], WT.log_law(yp_f[yp_f > 10], kappa=KAPPA, B=B_LOG))
    return out

A_values = np.linspace(0.0, 40.0, 11)        # the damping constant swept by the slider
figF6 = slider_figure(f6, "A+", A_values, xlabel="y+ [–]", ylabel="U+ [–]", title="κ turns the line, the wall damping A+ slides it up and down", yrange=(0, 28))
figF6.update_xaxes(type="log")                                # logarithmic y⁺ axis
figF6.show()
print(pd.DataFrame([{"A+": A_, **{f"B(κ={k_})": round(ch12.mixing_length_intercept(k_, A_plus=(A_ if A_ > 0 else None)), 2) for k_ in kappas}} for A_ in (0.0, 20.0, 26.0, 30.0)]).to_string(index=False))
""", explain=r"""
1. `f6(A_plus)` returns the model profile for three values of $\kappa$ and the reference log line.
2. The printed table is the intercept $B$ the model implies: at $\kappa=0.41$ it goes from −1.23 (no damping) through 4.01
   ($A^+=20$) and 5.28 ($A^+=26$) to 6.06 ($A^+=30$).""")
nb.figure_notes(
    see="Drag $A^+$ from 0 to 40: the three model curves keep their slopes at large $y^+$ and slide upward together; near "
        "$A^+\\approx25$ the $\\kappa=0.41$ curve lies on the reference line. The three $\\kappa$ curves fan out — a smaller "
        "$\\kappa$ gives a steeper line.",
    read="**$\\kappa$ sets the slope, the wall damping sets the intercept.** Neither constant comes from theory: $\\kappa$ is "
         "read from measured profiles and $A^+$ is tuned until $B$ comes out right.",
    change="…the wall were rough: no damping would be needed at all (there is no viscous sublayer to protect), and the "
           "intercept would be set by the roughness length instead.")
P("P302", "a nonlinear diffusion problem by Picard iteration on the eddy viscosity", r"""
The eddy viscosity depends on the answer (on $dU/dy$). *Picard iteration* (Ch. 8's idea of repeated substitution, applied to a
boundary-value problem): freeze $\nu_T$ from the last guess, solve the now **linear** problem
$\frac{d}{dy}\big([\nu+\nu_T]\frac{dU}{dy}\big)=\frac1\rho\frac{dP}{dx}$, update $\nu_T$, repeat until nothing changes.
`ch12.shear_flow_eddy_viscosity_solve` does this on a grid. The demo below iterates the first integral at one point.""", code=r"""
l2, s = (0.41*10.0)**2, 1.0                       # squared mixing length at y⁺ = 10 and a first guess for the slope s = dU⁺/dy⁺
for sweep in range(6):                            # six Picard sweeps on  (1 + ν_T/ν) s = 1  with  ν_T/ν = l² s
    s_new = 1.0/(1.0 + l2*s)                      # freeze ν_T from the old slope, solve the linear equation for the new slope
    s_new = 0.5*s + 0.5*s_new                     # under-relaxation: move only half-way (keeps the iteration from oscillating)
    print(f"sweep {sweep + 1}: s = {s_new:.4f}, change {abs(s_new - s):.4f}"); s = s_new   # the change shrinks every sweep; s → 0.216
""")
note("N181 [B]", r"""
**The whole channel with one mixing length** (ours). Combine the linear total stress of D17 with the damped mixing length,
$(1-y/\delta)u_*^2=\nu\frac{dU}{dy}+l_T^2\big(\frac{dU}{dy}\big)^2$ with $l_T=\kappa y[1-e^{-y^+/A^+}]$ (capped in the core, where
eddies are limited by the channel, not by the wall). `ch12.channel_mixing_length` solves it — a **model**, approximate against
DNS.""")
fig(r"""
ch1000 = ch12.channel_mixing_length(1000, kappa=KAPPA, A_plus=26.0)           # the model channel at Re_τ = 1000
d1000 = dns[dns["case"] == 1000]                                              # public DNS at Re_τ ≈ 1000 (Lee & Moser 2015)
y_half = np.concatenate([[0.0], np.geomspace(5e-4, 1.0, 120)])                # from one wall to the centreline, points crowded toward the wall
y_lam = np.concatenate([y_half, 2.0 - y_half[-2::-1]])                        # mirrored: the full channel in units of the half-height δ
lam = ch12.shear_flow_eddy_viscosity_solve(y_lam, lambda yf, sf: 0.0*yf, -1.0, 1.0, 1.0/1000.0)   # ν_T ≡ 0: laminar flow at the SAME pressure gradient
mix = ch12.shear_flow_eddy_viscosity_solve(y_lam, lambda yf, sf: (0.41*np.minimum(yf, 2.0 - yf))**2*np.abs(sf), -1.0, 1.0, 1.0/1000.0)   # plain mixing length l_T = κ × wall distance
print(f"laminar centreline U+ = {lam['U'].max():.1f} (= Re_τ/2), model centreline U+ = {ch1000['U_cl_plus']:.1f}: ratio {lam['U'].max()/ch1000['U_cl_plus']:.0f};  Picard converged: {mix['converged']} in {mix['iterations']} sweeps")
fig, ax = plt.subplots(figsize=(6.6, 3.8))
ax.semilogx(ch1000["yplus"], ch1000["Uplus"], color=C_MEAN, lw=2, label="mixing-length model with wall damping")
ax.semilogx(d1000["yplus"].values[1::4], d1000["Uplus"].values[1::4], "o", ms=3, color=C_RS, label="channel DNS, Re_τ ≈ 1000 (Lee & Moser 2015)")
ax.semilogx(y_lam[1:121]*1000, mix["U"][1:121], ":", color=C_FLUC, label="plain $l_T=\\kappa y$ (no damping, Picard solve)")
ax.semilogx(y_lam[1:121]*1000, lam["U"][1:121], "--", color=C_REF, label="laminar parabola, same pressure gradient")
ax.set(xlabel="y+ [–]", ylabel="U+ [–]", xlim=(0.5, 1000), ylim=(0, 30), title="One guessed length reproduces the mean profile — approximately")
ax.legend(fontsize=7)
print(f"largest |model − DNS| in U+: {np.abs(np.interp(d1000['yplus'].values[1:], ch1000['yplus'], ch1000['Uplus']) - d1000['Uplus'].values[1:]).max():.2f}")
""",
    see="The purple model curve passes through the orange DNS dots from the sublayer to the centreline. The dotted curve "
        "(no wall damping) has the same slope in the log region but lies about six units lower. The dashed laminar parabola "
        "shoots off the top of the plot — its centreline value is 500.",
    read="For the same pressure gradient the laminar flow would be about 20 times faster: turbulence is an enormous extra "
         "resistance. The damped model is within about one unit of $U^+$ of the DNS — good for a one-constant-plus-one-fix "
         "model, and **only** because the two numbers were tuned on flows like this one.",
    change="…the flow separated from the wall, or curved strongly: $l_T=\\kappa y$ would have no justification, and the model "
           "fails. That is why models that carry their own scales were invented (C13).")
note("N173 [B]", r"""
**A velocity scale for convection.** A parcel warmer than its surroundings by $T'$ accelerates upward (for a perfect gas
$\alpha=1/T$):""",
     equation=r"Dw/Dt\sim g\alpha T'\sim g\Delta T/T", ref="12.102")
note("N174 [B]", r"""
**Free fall over the layer depth.** If the parcel accelerates over the depth $L$ of the layer, $w^2/L\sim g\Delta T/T$, so
$w\sim\sqrt{gL\Delta T/T}$ and the eddy diffusivity is $\kappa_T\sim wL$. Our own example: a 1 km deep convecting layer with a
2 K temperature contrast at 300 K.""")
code(r"""
w_c = ch12.convective_velocity_scale(1000.0, 2.0, 300.0)      # w ~ sqrt(g L ΔT/T) for L = 1 km, ΔT = 2 K, T = 300 K [m/s]
K_c = ch12.convective_eddy_diffusivity(1000.0, 2.0, 300.0)    # κ_T ~ w L [m²/s]
print(f"w ≈ {w_c:.1f} m/s, κ_T ≈ {K_c:.1e} m²/s — {K_c/2e-5:.0e} times the molecular 2e-5 m²/s")
""", explain=r"""
Thermals of about 8 m/s and an eddy diffusivity near $8\times10^3$ m²/s, some $4\times10^8$ times the molecular value: an
order-of-magnitude estimate (no constant of proportionality is claimed), but it explains why a convective boundary layer is
well mixed within tens of minutes.""")
explainer("mixing_length_closure", "What does a closure constant do?",
          "sliding κ rotates the log line and sliding A⁺ shifts it, with the implied intercept B read off live; switching the "
          "model off returns the laminar parabola at the same pressure gradient.",
          ["Set A⁺ = 0: the line keeps its slope but B drops to −1.2.",
           "Set A⁺ = 26: B = 5.3.",
           "Change κ from 0.41 to 0.38 and watch the slope steepen.",
           "Switch the model off: the laminar parabola is 20 times faster at Re_τ = 1000."])
whatif(r"""…no single length like $\kappa y$ is available — a separated flow, a jet in a cross-wind? Then the model must carry
its own velocity and length scales: C13.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C13", r"The modelled energy equation $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\frac{\nu_T}{\sigma_e}\frac{\partial\bar e}{\partial x_j}\big)-\bar\varepsilon-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}$ (12.103) — the k–ε model",
     "Where do the constants of a turbulence model come from?")
problem(r"""
Most engineering flow software, and the turbulence closures inside ocean and boundary-layer schemes, carry two extra fields:
how much turbulent energy there is and how fast it is being dissipated. From those two they build an eddy viscosity
everywhere. The price is five constants. Are they fudge factors? Two of them are pinned by experiments you can picture.
(⚠️ the "k" of k–ε is the book's $\bar e$; the book's $k$ is a wavenumber or a conductivity.)""")
idea("", r"""
A velocity scale $u_T=\sqrt{\bar e}$. A length scale: how far an eddy of that speed goes in its lifetime $\bar e/\bar\varepsilon$,
so $l_T=\bar e^{3/2}/\bar\varepsilon$. Their product is the eddy viscosity, $\nu_T\propto u_Tl_T=\bar e^2/\bar\varepsilon$.""")
remind("C13")
note("N175 [B]", r"""
**Ingredients of a one-equation model** (one transport equation, for $\bar e$; the length $l_T$ is still prescribed): a velocity
scale $u_T=c\sqrt{\bar e}$; a modelled dissipation $\bar\varepsilon=C_\varepsilon\bar e^{3/2}/l_T$; and the three transport terms of
the exact budget (C06) lumped into one gradient flux,
$-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}=\frac{\nu_T}{\sigma_e}\frac{\partial\bar e}{\partial x_j}$
(`ch12.one_equation_closure`).

> ⚠️ **slip #8 — the book prints** the viscous transport here as $2\nu\overline{u_jS'_{ij}}$ **; the correct form is**
> $2\nu\overline{u_iS'_{ij}}$, as in the exact budget.""")
NOTES_IN["D22"] = (r"**N176 [B]** the eddy viscosity of the k–ε model, $\nu_T=C_\mu[\bar e^{3/2}/\bar\varepsilon]\sqrt{\bar e}=C_\mu\bar e^2/\bar\varepsilon$ (12.104) (steps 4–5)")
D("D22", ref="12.103")
note("N177 [B]", r"""
**The modelled dissipation equation** is built by analogy with the energy equation, not derived:

$$\frac{\partial\bar\varepsilon}{\partial t}+U_j\frac{\partial\bar\varepsilon}{\partial x_j}=\frac{\partial}{\partial x_j}\Big(\frac{\nu_T}{\sigma_\varepsilon}\frac{\partial\bar\varepsilon}{\partial x_j}\Big)-C_{\varepsilon1}\Big(\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}\Big)\frac{\bar\varepsilon}{\bar e}-C_{\varepsilon2}\frac{\bar\varepsilon^2}{\bar e}\qquad\text{(12.105)}$$

Term by term: change following the mean flow = gradient diffusion + (production of $\bar e$) ÷ (the time scale
$\bar e/\bar\varepsilon$) − ($\bar\varepsilon$) ÷ (the same time scale), each with its own constant.""")
note("N178 [B]", r"""
**The standard constants** (Launder & Sharma 1974; public values, `ch12.K_EPSILON_CONSTANTS`): $C_\mu=0.09$,
$C_{\varepsilon1}=1.44$, $C_{\varepsilon2}=1.92$, $\sigma_e=1.0$, $\sigma_\varepsilon=1.3$.

> ⚠️ **slip #9 — the book prints** "five" constants and then lists six entries **; the correct form is** five:
> $C_\mu,C_{\varepsilon1},C_{\varepsilon2},\sigma_e,\sigma_\varepsilon$.""")
note("N179 [C]", r"""
**The closed set.** Mean continuity $\partial U_i/\partial x_i=0$ (12.27); the mean momentum equation
$\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=\frac1\rho\frac{\partial\bar\tau_{ij}}{\partial x_j}$ (12.30, constant density);
the closure $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94); the $\bar e$
equation of D22 with $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ (12.104); and the $\bar\varepsilon$ equation of N177 above. At a wall the log
law of C11 is imposed as a boundary condition (a "wall function"); results are sensitive to the inlet values of $\bar e$ and
$\bar\varepsilon$. *Reynolds-stress closures* model the transport equation for $\overline{u_iu_j}$ (N59 in C04) instead of assuming
an eddy viscosity (named in Ch. 10).""")
code(r"""
print(ch12.K_EPSILON_CONSTANTS)                               # the five standard constants
one = ch12.one_equation_closure(1.0, 0.1, 0.55, 0.17, 1.0)    # ē = 1 m²/s², l_T = 0.1 m; c, C_ε, σ_e ILLUSTRATIVE
print({k: round(float(v), 3) for k, v in one.items()})        # u_T, ν_T = l_T u_T, modelled ε̄, transport diffusivity
print("dē/dt, dε̄/dt for ē = ε̄ = 1 and production 0.5:", ch12.k_epsilon_rhs(1.0, 1.0, 0.5), "| with no production:", ch12.k_epsilon_rhs(1.0, 1.0, 0.0))
""", explain=r"""
1. `ch12.one_equation_closure(e, l_T, c, C_eps, sigma_e)` returns the ingredients of **N175** for illustrative constants.
2. `ch12.k_epsilon_rhs(e, eps, production)` evaluates the source terms of the two model equations for homogeneous turbulence:
   $d\bar e/dt=P-\bar\varepsilon=-0.5$ and $d\bar\varepsilon/dt=(1.44\times0.5-1.92)\,\bar\varepsilon/\bar e=-1.2$; with no production,
   $(-1.0,-1.92)$ — the pair of ODEs the next derivation integrates.""")
P("P303", "dividing two ODEs to eliminate time", r"""
If $dy/dt=f$ and $dx/dt=g$, then $dy/dx=f/g$: time disappears and $y$ is found as a function of $x$. If the result is
$dy/dx=n\,y/x$, the solution is the power law $y\propto x^n$.""", code=r"""
x_s = sp.symbols("x", positive=True); y_s = sp.Function("y")  # the independent variable and the unknown function
print(sp.dsolve(sp.Eq(y_s(x_s).diff(x_s), 2*y_s(x_s)/x_s)))   # dy/dx = 2 y/x  →  y = C1 x²
""")
NOTES_IN["D23"] = (r"**N180 [B]** what fixes two of the constants: the decay exponent of grid turbulence, $n=1/(C_{\varepsilon2}-1)$, and the "
                   r"von Kármán constant in the log layer, $\kappa^2=\sqrt{C_\mu}(C_{\varepsilon2}-C_{\varepsilon1})\sigma_\varepsilon$ (steps 1–12)")
D("D23", ref="12.105")
nb.worked_example("what the standard constants imply", r"""
1. Decay exponent $n=1/(C_{\varepsilon2}-1)=1/0.92=1.087$ — grid-turbulence experiments give $n\approx1.1$–1.3:
   $C_{\varepsilon2}$ was chosen to match.
2. Log layer: $\kappa^2=\sqrt{C_\mu}(C_{\varepsilon2}-C_{\varepsilon1})\sigma_\varepsilon=0.3\times0.48\times1.3=0.187$, $\kappa=0.433$ — a
   little above the measured 0.38–0.41: the constants are a compromise.
3. In the log layer $\bar e/u_*^2=1/\sqrt{C_\mu}=3.33$.
4. Decay from $\bar e_0=1$ m²/s², $\bar\varepsilon_0=1$ m²/s³: $t_0=n\bar e_0/\bar\varepsilon_0=1.087$ s; at $t=10$ s,
   $\bar e=(1+10/1.087)^{-1.087}=0.080$ m²/s².""")
code(r"""
t_d = np.array([0.0, 1.0, 10.0])                              # three times [s]
e_d, eps_d, n_d, t0_d = ch12.k_epsilon_decay(1.0, 1.0, t_d)   # decaying homogeneous turbulence from ē0 = 1 m²/s², ε̄0 = 1 m²/s³
print(f"n = {n_d:.4f}, t0 = {t0_d:.4f} s;  ē = {np.round(e_d, 4)} m²/s²,  ε̄ = {np.round(eps_d, 5)} m²/s³")
print(f"κ implied by the constants: {ch12.k_epsilon_loglayer_kappa(0.09, 1.44, 1.92, 1.3):.4f}")
print(f"ν_T = C_μ ē²/ε̄ at t = 0, 10 s: {ch12.k_epsilon_eddy_viscosity(e_d[0], eps_d[0]):.4f}, {ch12.k_epsilon_eddy_viscosity(e_d[2], eps_d[2]):.4f} m²/s;  "
      f"l_T = ē^1.5/ε̄: {ch12.k_epsilon_length_scale(e_d[0], eps_d[0]):.2f} → {ch12.k_epsilon_length_scale(e_d[2], eps_d[2]):.2f} m")
""", explain=r"""
1. `ch12.k_epsilon_decay(e0, eps0, t)` returns $\bar e(t)$, $\bar\varepsilon(t)$, the exponent $n$ and the virtual origin $t_0$ (the
   closed form of D23): $\bar e$ falls to 0.49 after 1 s and 0.080 after 10 s.
2. `ch12.k_epsilon_loglayer_kappa` is the second relation of D23: 0.433.
3. During the decay the eddy viscosity falls slowly while the length scale $l_T=\bar e^{3/2}/\bar\varepsilon$ **grows**: the small
   eddies die first, leaving the big ones.""")
scratch(r"""
# From scratch: the two decay ODEs with a hand-written fourth-order Runge–Kutta loop (Ch. 3, P95)
def rhs(state):                                               # dē/dt = −ε̄,  dε̄/dt = −C_ε2 ε̄²/ē  (no production, no transport)
    e_, eps_ = state
    return np.array([-eps_, -1.92*eps_**2/e_])

state, dt_rk = np.array([1.0, 1.0]), 0.01                     # initial ē, ε̄ and the time step [s]
for _ in range(1000):                                         # 1000 steps to t = 10 s
    k1 = rhs(state); k2 = rhs(state + 0.5*dt_rk*k1); k3 = rhs(state + 0.5*dt_rk*k2); k4 = rhs(state + dt_rk*k3)   # four slopes
    state = state + dt_rk*(k1 + 2*k2 + 2*k3 + k4)/6           # RK4 update
assert np.allclose(state, [e_d[2], eps_d[2]], rtol=1e-6)      # same ē(10) and ε̄(10) as the closed form
print("RK4 at t = 10 s:", np.round(state, 5))
""", r"""
A plain RK4 loop on the two ODEs lands on the closed form at $t=10$ s to six digits: the closed form of D23 is the solution.""")
fig(r"""
t_ax = np.geomspace(0.05, 100.0, 200)                         # times [s]
e_c, eps_c, n_c, t0_c = ch12.k_epsilon_decay(1.0, 1.0, t_ax)  # closed form
t_iv = np.geomspace(0.05, 100.0, 12)                          # a few times for the numerical cross-check
e_iv, eps_iv, _, _ = ch12.k_epsilon_decay(1.0, 1.0, t_iv, method="ivp")   # the same by solve_ivp
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 3.6))
a1.loglog(t_ax, e_c, color=C_FLUC, label="$\\bar e$ (closed form)"); a1.loglog(t_iv, e_iv, "o", ms=4, color=C_FLUC)
a1.loglog(t_ax, eps_c, color=C_VISC, label="$\\bar\\varepsilon$ (closed form)"); a1.loglog(t_iv, eps_iv, "o", ms=4, color=C_VISC)
a1.loglog(t_ax[t_ax > 8], 9*t_ax[t_ax > 8]**(-n_c), ":", color=C_REF, label=f"slopes −n = −{n_c:.2f} and −(n+1)")
a1.loglog(t_ax[t_ax > 8], 9*t_ax[t_ax > 8]**(-n_c - 1), ":", color=C_REF)
a1.set(xlabel="time t [s]", ylabel="$\\bar e$ [m²/s²], $\\bar\\varepsilon$ [m²/s³]", title="Decaying turbulence: a power law")
a1.legend(fontsize=7)
a2.loglog(t_ax, ch12.k_epsilon_eddy_viscosity(e_c, eps_c), color=C_RS, label="$\\nu_T=C_\\mu\\bar e^2/\\bar\\varepsilon$ [m²/s]")
a2.loglog(t_ax, ch12.k_epsilon_length_scale(e_c, eps_c), color=C_MEAN, label="$l_T=\\bar e^{3/2}/\\bar\\varepsilon$ [m]")
a2.set(xlabel="time t [s]", title="The eddies that are left grow")
a2.legend(fontsize=7)
""",
    see="Left: on log–log axes $\\bar e$ (teal) and $\\bar\\varepsilon$ (rose) bend over near $t\\approx t_0=1.09$ s and then fall "
        "along straight lines of slopes −1.09 and −2.09 (dotted guides); the dots (`solve_ivp`) sit on the lines (closed form). "
        "Right: the eddy viscosity falls slowly while the length scale grows.",
    read="A straight line of slope −1.09: one constant, one measurable exponent. This is how $C_{\\varepsilon2}$ was fixed — by "
         "the decay of turbulence behind a grid in a wind tunnel.",
    change="…$C_{\\varepsilon2}=2.0$: then $n=1.0$ exactly, energy $\\propto1/t$, and the eddy viscosity would stay constant during "
           "the decay.")
whatif(r"""…gravity acts on density fluctuations? The budget gains a buoyancy term that can feed the turbulence or starve it
(C14).""")


# =====================================================================================================================
# A.11 §12.11 Turbulence in a Stratified Medium — R05 R06, C14, C15
# =====================================================================================================================
nb.section("12.11", "Turbulence in a Stratified Medium", intro=r"""
**What is this section about?** Turbulence where density changes with height — the atmosphere and the ocean. Buoyancy can feed
turbulence (heated ground at noon) or drain it (a clear calm night). Two numbers measure the contest, and one length says at
what height buoyancy takes over.

> ⚠️ **Buoyancy is back.** The constant-density assumption of §12.8–12.10 ends here; we return to the Boussinesq equations of
> §12.5–12.7.""")
nb.recap("R05", "Stratification, potential temperature and the two lapse-rate conventions", r"""
Stability is set by how the temperature gradient compares with the adiabatic one (Ch. 1), so "temperature" in this section
means **potential** temperature $\theta$: stable when $d\theta/dz>0$. In terms of the thermometer temperature:

| Convention | Lapse rate | Adiabatic value | Buoyancy frequency | Stable when |
|---|---|---|---|---|
| Kundu (computed with) | $\Gamma\equiv dT/dz$ | $\Gamma_a=-g/C_p\approx-9.8$ K/km | $N^2=g\alpha(dT/dz-\Gamma_a)$ | $dT/dz>\Gamma_a$ |
| meteorology (shown alongside) | $\Gamma_{met}\equiv-dT/dz$ | $\Gamma_d\approx+9.8$ K/km | $N^2=g\alpha(\Gamma_d-\Gamma_{met})$ | $\Gamma_{met}<\Gamma_d$ |

Negate the number, flip the inequality (Ch. 1, P48). Where these come from: Ch. 1's
$N^2=-\frac g{\rho(z_o)}\big(\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\big)$ (1.29) and
$\frac T\theta\frac{d\theta}{dz}=\frac{dT}{dz}+\frac g{C_p}=\Gamma-\Gamma_a$ (1.32); Ch. 7's incompressible form
$N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ (7.127).""", where="Ch. 1 §1.10; Ch. 7; Ch. 11 §11.7")
nb.code(r"""
Gamma_a = STRAT.adiabatic_lapse_rate()                        # the adiabatic gradient with Kundu's sign: −g/c_p ≈ −9.76e-3 K/m
print(f"Γa = {Gamma_a*1e3:.2f} K/km (Kundu)  ⇔  Γd = {-Gamma_a*1e3:.2f} K/km (meteorology)")
for dTdz in (-6.5e-3, Gamma_a, -12e-3):                       # three in-situ gradients: standard atmosphere, adiabatic, super-adiabatic
    kundu = STRAT.lapse_rate_stability(dTdz, Gamma_a=Gamma_a, convention="kundu")          # verdict with the criterion in Kundu's convention
    met = STRAT.lapse_rate_stability(dTdz, Gamma_a=Gamma_a, convention="meteorology")      # the same layer in the meteorological convention
    print(f"dT/dz = {dTdz*1e3:6.2f} K/km | Kundu: {kundu.text} | meteorology: {met.text}")
""", explain=r"""
1. `STRAT.adiabatic_lapse_rate()` (Ch. 1) returns $\Gamma_a=-g/C_p$ with Kundu's sign; we never type −0.0098 by hand.
2. `STRAT.lapse_rate_stability` writes the criterion out with the numbers, in either convention. Each printed text **begins with
   the criterion for stability** ("stable ⇔ …") and then shows the actual comparison: `>` (stable) for the standard atmosphere,
   `=` (neutral) for the adiabatic gradient, `<` (unstable) for −12 K/km. ⚠️ In the meteorological line the library prints
   "Γ < Γa": read it as $\Gamma_{met}<\Gamma_d$ — its "Γ" is our $\Gamma_{met}\equiv-dT/dz$ and its "Γa" the positive dry-adiabatic
   rate $\Gamma_d\approx+9.8$ K/km. The same strings appear in the verdicts printed in C14 and C15.""")
core("C14", r"The flux Richardson number $\mathrm{Rf}=\frac{-g\alpha\overline{wT'}}{-\overline{uw}(dU/dz)}$ (12.107), $\mathrm{Ri}\equiv\frac{N^2}{(dU/dz)^2}$ (12.108) and $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109)",
     "When does stratification switch turbulence off?")
problem(r"""
On a clear night the ground cools, the air near it becomes heavy, and the wind you felt at dusk dies near the surface while it
still blows at tree-top height: the turbulence that carried momentum down has been suppressed. In the ocean the same contest
decides how much heat is mixed below the thermocline. A mixing scheme must decide, level by level, whether the shear can keep
turbulence alive against the stratification.""")
idea("", r"""
Lifting heavy fluid and pushing light fluid down costs energy; the turbulence pays it out of what the shear supplies.
**Rf = cost / income.** If the cost takes more than about a quarter of the income, dissipation takes the rest and more — the
turbulence starves.

| Rf | What it means |
|---|---|
| $\mathrm{Rf}<0$ | heat flux upward: buoyancy *adds* energy (convective) |
| $0<\mathrm{Rf}<\approx\tfrac14$ | shear-driven turbulence, weakened by stratification |
| $\mathrm{Rf}\gtrsim\tfrac14$ | decaying: the turbulence cannot sustain itself |

Symbols: $z$ height [m], $w$ the vertical velocity fluctuation, $\overline{wT'}$ the kinematic heat flux [K m/s] (positive
upward), $\overline{uw}$ the kinematic momentum flux [m²/s²], $\alpha$ the thermal expansion coefficient [1/K].""")
remind("C14")
NOTES_IN["D24"] = (r"**N182 [B]** the reduced budget of a horizontally uniform stratified shear flow, "
                   r"$\frac{\partial\bar e}{\partial t}+U\frac{\partial\bar e}{\partial x}=-\frac{\partial}{\partial z}\big(\frac1{\rho_0}\overline{pw}+\frac12\overline{u_i^2w}\big)-\overline{uw}\frac{\partial U}{\partial z}+g\alpha\overline{wT'}-\bar\varepsilon$ (12.106) "
                   r"— ⚠️ slip #10: the page prints the triple correlation as $\overline{ew}$ (steps 1–3) · **N184 [B]** the link between the "
                   r"two Richardson numbers, $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109) (steps 6–8)")
D("D24", ref="12.107")
nb.recap("R06", "The gradient Richardson number", r"""
Ch. 11's $\mathrm{Ri}\equiv N^2/(dU/dz)^2$ (11.66), with linear stability guaranteed if, everywhere in the flow,
$\mathrm{Ri}>\tfrac14$ (11.67). Here the same number is written with temperature,
$\mathrm{Ri}=\alpha g(d\bar T/dz)/(dU/dz)^2$ (12.108), with $\bar T$ the **potential** temperature. New in this chapter:
`ch12.gradient_richardson_thermal` takes the thermometer (in-situ) gradient and $\Gamma_a$ (Kundu sign) and returns Ri with the
verdict in both conventions.""", where="Ch. 11 §11.7")
nb.current_core = "C14"
confusion(r"""**which temperature gradient goes into Ri?** In $\mathrm{Ri}=N^2/(dU/dz)^2$ (12.108) the stratification is
$N^2=g\alpha(dT/dz-\Gamma_a)$ — the in-situ gradient **minus the adiabatic one** (Kundu: $\Gamma_a\approx-9.8$ K/km; meteorology:
$N^2=g\alpha(\Gamma_d-\Gamma_{met})$ with $\Gamma_d\approx+9.8$ K/km). Forgetting $\Gamma_a$ (using $dT/dz$ alone) would call the
standard atmosphere, $dT/dz=-6.5$ K/km, *unstable* — it is stable: $-6.5>-9.8$ ⇔ $6.5<9.8$. That is why
`ch12.gradient_richardson_thermal` has no default for `Gamma_a`: pass the adiabatic value for a thermometer gradient, and
`Gamma_a=0.0` only when the gradient you pass is already a potential-temperature gradient.""")
note("N183 [B]", r"""
**The critical flux Richardson number.** Turbulence stops being self-supporting near $\mathrm{Rf}_{cr}\approx0.25$ — an
**observation**, not a theorem; a large negative Rf means convection dominates (`ch12.turbulence_regime`).

> ⚠️ **Common confusion:** two different "one quarter" statements. Ch. 11's theorem — $\mathrm{Ri}>\tfrac14$ everywhere ⇒ a
> *laminar* stratified shear flow is linearly stable — is sufficient and exact. The observed $\mathrm{Rf}_{cr}\approx\tfrac14$ is
> about *existing turbulence*; by $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109) it corresponds to $\mathrm{Ri}=\mathrm{Pr}_T/4$.""")
note("N185 [B]", r"""
**The turbulent Prandtl number** $\mathrm{Pr}_T=\nu_T/\kappa_T$ compares how well eddies mix momentum and heat. It is larger
than 1 in stable conditions (internal waves, Ch. 7, carry momentum but not heat), small in unstable ones, and about 1 when
neutral (the *Reynolds analogy*). ⚠️ In the book's $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109) the eddy viscosity $\nu_T$ is
set like an italic "v".""")
nb.worked_example("a stable night", r"""
$\overline{wT'}=-0.02$ K m/s (downward), $\overline{uw}=-0.09$ m²/s² ($u_*=0.3$ m/s), $dU/dz=0.1$ s⁻¹, $\alpha=1/300$ K⁻¹,
$g\approx9.81$ m/s².

1. Buoyant destruction $-g\alpha\overline{wT'}=9.81\times0.02/300=6.54\times10^{-4}$ m²/s³.
2. Shear production $-\overline{uw}\,dU/dz=0.09\times0.1=9.0\times10^{-3}$ m²/s³.
3. $\mathrm{Rf}=0.073$ — shear-driven.
4. Thermometer gradient $dT/dz=+10$ K/km (an inversion). Kundu: $d\theta/dz=dT/dz-\Gamma_a=0.010-(-0.0098)=0.0198$ K/m.
   Meteorology: $\Gamma_{met}=-10$ K/km, $\Gamma_d-\Gamma_{met}=9.8+10=19.8$ K/km — the same number.
5. $N^2=g\alpha\,d\theta/dz=6.47\times10^{-4}$ s⁻²; $\mathrm{Ri}=N^2/(dU/dz)^2=0.065$.
6. $\nu_T=0.09/0.1=0.9$ m²/s; $\kappa_T=0.02/0.0198=1.01$ m²/s; $\mathrm{Pr}_T=0.89$; and $\mathrm{Ri}/\mathrm{Rf}=0.065/0.073=0.89$ ✓.
   (These are our round numbers, chosen for easy arithmetic, not a measured case — in measured stable layers $\mathrm{Pr}_T$ is
   usually somewhat above 1, see N185.)""")
code(r"""
Rf = ch12.flux_richardson(-0.02, -0.09, 0.1, 1/300.0)        # Eq. (12.107): mean(wT') = −0.02 K m/s, mean(uw) = −0.09 m²/s², dU/dz = 0.1 1/s
ri = ch12.gradient_richardson_thermal(0.010, 0.1, 1/300.0, Gamma_a=Gamma_a)   # Eq. (12.108) from the IN-SITU gradient +10 K/km; Γa from the Ch. 1 function
print(f"Rf = {Rf:.4f} → {ch12.turbulence_regime(Rf)};  Ri = {ri['Ri']:.4f}, N² = {ri['N2']:.2e} 1/s², dθ/dz = {ri['dthetadz']*1e3:.2f} K/km")
print("Kundu:      ", ri["verdict_kundu"])                    # the criterion, then the numbers
print("meteorology:", ri["verdict_met"])
Pr_T = ch12.turbulent_prandtl(0.9, 0.02/ri["dthetadz"])       # ν_T/κ_T with ν_T = 0.9 m²/s and κ_T = −mean(wT')/(dθ/dz)
print(f"Pr_T = {Pr_T:.3f};  Ri/Rf = {ri['Ri']/Rf:.3f};  Rf back from Ri: {ch12.flux_from_gradient_richardson(ri['Ri'], Pr_T):.4f}")
bud_s = ch12.stratified_tke_budget(10.0, 3.0, -0.09, -0.02, 8.3e-3, 1/300.0, dUdz=0.1)   # the budget of D24 at one height (ε̄ chosen to close it)
print({k: (round(float(v), 5) if not isinstance(v, str) else v) for k, v in bud_s.items()})
wrong = ch12.gradient_richardson_thermal(-6.5e-3, 0.02, 1/288.0, Gamma_a=0.0)["verdict"]     # the trap: the standard atmosphere WITHOUT the adiabatic correction
right = ch12.gradient_richardson_thermal(-6.5e-3, 0.02, 1/288.0, Gamma_a=Gamma_a)["verdict"]  # with it
print(f"standard atmosphere (dT/dz = −6.5 K/km ⇔ Γ_met = 6.5 K/km): forgetting Γa says '{wrong}', the correct verdict is '{right}'")
""", explain=r"""
1. `ch12.flux_richardson(wT, uw, dUdz, alpha)` is cost/income: 0.073 (the function's default $g$ is the standard 9.80665 m/s², so
   the last digit differs from the hand value with 9.81).
2. `ch12.gradient_richardson_thermal(dTdz, dUdz, alpha, Gamma_a=…)` takes the **thermometer** gradient; `Gamma_a` is a required
   keyword. It returns Ri, $N^2$, $d\theta/dz$ and the verdict in both conventions (each string starts with the criterion "stable ⇔ …";
   the meteorological line's "Γ < Γa" means $\Gamma_{met}<\Gamma_d$: here $-10<9.8$ K/km).
3. `ch12.turbulent_prandtl` and `ch12.flux_from_gradient_richardson` close the loop: $\mathrm{Ri}/\mathrm{Rf}=\mathrm{Pr}_T=0.89$.
4. `ch12.stratified_tke_budget` lists the terms of the reduced budget at one height with their Rf and regime: buoyancy removes
   about 7 % of the shear production.
5. The last two lines are the trap of the ⚠️ box above, run on the standard atmosphere.""")
scratch(r"""
# From scratch: Rf and Ri by hand (Kundu sign), and Ri again from the potential-temperature gradient directly
g_ = G0                                                       # standard gravity [m/s²], the library's default
Rf_mine = (-g_*(1/300.0)*(-0.02))/(-(-0.09)*0.1)              # Eq. (12.107): buoyant destruction / shear production
dth = 0.010 - Gamma_a                                         # dθ/dz = dT/dz − Γa (Kundu)  =  Γd − Γ_met (meteorology) [K/m]
Ri_mine = g_*(1/300.0)*dth/0.1**2                             # Eq. (12.108) with the potential-temperature gradient
Ri_pot = ch12.gradient_richardson_thermal(dth, 0.1, 1/300.0, Gamma_a=0.0)["Ri"]   # the library fed a POTENTIAL-temperature gradient: Γa = 0 on purpose
assert np.isclose(Rf_mine, Rf) and np.isclose(Ri_mine, ri["Ri"]) and np.isclose(Ri_pot, ri["Ri"])   # all routes agree
assert np.isclose(Ri_mine/Rf_mine, Pr_T, rtol=1e-3)           # Eq. (12.109): Ri = Pr_T Rf
print(f"Rf = {Rf_mine:.4f}, Ri = {Ri_mine:.4f}, Ri/Rf = {Ri_mine/Rf_mine:.3f}")
""", r"""
Two one-line formulas reproduce the functions; feeding the library the potential-temperature gradient with `Gamma_a=0.0` (the one
legitimate use of that value) gives the same Ri.""")
fig(r"""
u_st, z_b, alpha_b = 0.3, 10.0, 1/300.0                       # friction velocity [m/s], height [m], expansion coefficient [1/K]
P_sh = u_st**3/(KAPPA*z_b)                                    # shear production in the log layer: u*² × u*/(κ z) [m²/s³]
cases = (("unstable noon, H = +150 W/m²\ndT/dz < Γa = −9.8 K/km ⇔ Γ_met > Γd = 9.8 K/km", 150.0),
         ("neutral, H = 0\ndT/dz = Γa ⇔ Γ_met = Γd", 0.0),
         ("stable night, H = −30 W/m²\ndT/dz > Γa = −9.8 K/km ⇔ Γ_met < Γd = 9.8 K/km", -30.0))   # three surface heat fluxes, both conventions in each label
fig, ax = plt.subplots(figsize=(8.6, 3.8))
for i, (name, H) in enumerate(cases):
    B_b = ch12.buoyant_production(H/(1.2*1005.0), alpha_b)    # g α mean(wT') with mean(wT') = H/(ρ c_p) [m²/s³]
    Rf_b = 0.0 if B_b == 0 else -B_b/P_sh                     # Eq. (12.107) (written so that the neutral case prints +0.00)
    ax.bar(i - 0.25, P_sh, 0.25, color=C_RS, label="shear production" if i == 0 else None)
    ax.bar(i, B_b, 0.25, color=C_BUOY, label="buoyancy $g\\alpha\\overline{wT'}$" if i == 0 else None)
    ax.bar(i + 0.25, -(P_sh + B_b), 0.25, color=C_VISC, label="dissipation (the residual)" if i == 0 else None)
    ax.text(i, 0.0125, f"Rf = {Rf_b:+.2f}\n{ch12.turbulence_regime(Rf_b)}", ha="center", fontsize=9)
ax.axhline(0, color=C_REF, lw=0.6)
ax.set_xticks(range(3)); ax.set_xticklabels([c[0] for c in cases], fontsize=7)
ax.set(ylabel="budget terms at z = 10 m [m²/s³]", ylim=(-0.013, 0.016), title="Buoyancy feeds the turbulence by day and drains it by night")
ax.legend(fontsize=7, loc="lower left")
""",
    see="Three groups of bars at 10 m height. The orange bar (shear production) is the same in all three. The blue bar "
        "(buoyancy) is positive at noon — about 60 % of the shear production —, zero when neutral and slightly negative at "
        "night. The rose bar (dissipation, computed as the residual) balances the other two.",
    read="Rf is minus the ratio of the blue bar to the orange one: −0.6 at noon (convective), 0, and +0.12 at night "
         "(shear-driven but weakened). Each label states the stratification in both lapse-rate conventions.",
    change="…the wind dropped by half at night: production $\\propto u_*^3$ falls 8-fold, Rf rises 8-fold to about 1 — far "
           "beyond ¼ — and the turbulence collapses. This is why calm clear nights are so still.")
whatif(r"""…we ask at what *height* buoyancy catches up with shear? In the log layer the answer is a single length (C15).""")

# ---------------------------------------------------------------------------------------------------------------------
core("C15", r"The Monin–Obukhov length $L_M\equiv-u_*^3/(\kappa\alpha g\overline{wT'})$ (12.110), with $\mathrm{Rf}=z/L_M$ (12.111)",
     "What does the Monin–Obukhov length measure?")
problem(r"""
A weather or climate model knows the wind and temperature at its lowest level, some tens of metres up, and needs the stress and
heat flux at the ground. The bridge is a wind profile that is logarithmic near the ground and bends with stability — and the
one length that says how fast it bends is $L_M$. Every bulk surface-flux formula is built on it.""")
idea("", r"""
In the log layer the shear production falls with height like $u_*^3/(\kappa z)$, while the buoyancy term
$g\alpha\overline{wT'}$ does not change with height in the surface layer (the fluxes are constant there). They are equal in
size at the height $z=\lvert L_M\rvert$. Below it the wind makes the turbulence; above it buoyancy makes it (unstable) or kills
it (stable). Sign: $L_M<0$ unstable (heat flux up), $L_M>0$ stable, $L_M=\pm\infty$ neutral.""")
remind("C15")
note("N187 [B]", r"""
**Forced and free convection.** On an unstable day the layer below $\lvert L_M\rvert$ is *forced convection* (turbulence made by
the wind shear, heat carried along); above it is *free convection*, with thermal plumes that do not care about the wind.""")
fig(r"""
fig, ax = plt.subplots(figsize=(6.6, 3.4))
surface_layer_sketch(ax, L_M=-20.0, top=100.0, seed=1)        # drawing: an unstable surface layer with |L_M| = 20 m
print([ch12.surface_layer_regime(z_, L_) for z_, L_ in ((2.0, -20.0), (100.0, -20.0), (10.0, 100.0), (300.0, 100.0), (10.0, np.inf))])   # the layer a height lies in
""",
    see="A surface layer on an unstable day: small shear-made eddies below a dashed line at $z=\\lvert L_M\\rvert=20$ m, tall "
        "buoyant plumes above it.",
    read="$\\lvert L_M\\rvert$ is the height at which buoyancy takes over from shear. The printed list is "
         "`ch12.surface_layer_regime(z, L_M)` for five (height, $L_M$) pairs: forced convection, free convection, forced "
         "convection, stable, neutral. Note that it says \"forced convection\" for any height well below $\\lvert L_M\\rvert$ "
         "**whatever the sign** of $L_M$ — the layer there is shear-driven; whether the stratification is stable or unstable must "
         "be read from the sign of $L_M$ (or from the `verdict` keys of `ch12.surface_layer_state`).",
    change="…the heat flux were weaker or the wind stronger: $\\lvert L_M\\rvert\\propto u_*^3/\\overline{wT'}$ grows and the "
           "shear-driven layer deepens — on a windy overcast day it fills the whole boundary layer.")
P("P304", "the stability parameter ζ = z/L and integrating a flux–profile relation", r"""
Surface-layer schemes write the dimensionless shear $\phi_m=(\kappa z/u_*)\,dU/dz$ as a function of $\zeta=z/L$ alone: 1 when
neutral, $>1$ stable, $<1$ unstable (*Monin–Obukhov similarity*). Integrating $dU/dz=(u_*/\kappa z)\phi_m$ from the roughness
height $z_0$ gives the wind profile. With $\phi_m=1+5\zeta$ the integral is $\ln(z/z_0)+5(z-z_0)/L$.""", code=r"""
z0_, L_ = 0.03, 83.0                              # roughness length [m] and a (stable) Monin–Obukhov length [m]
num = quad(lambda z: (1 + 5*z/L_)/z, z0_, 10.0)[0]    # ∫ φ_m dz/z from z0 to 10 m, done numerically
print(round(num, 4), round(np.log(10.0/z0_) + 5*(10.0 - z0_)/L_, 4))   # the closed form gives the same number
""")
NOTES_IN["D25"] = (r"**N186 [B]** in the surface layer the flux Richardson number is a height in units of $L_M$, $\mathrm{Rf}=z/L_M$ (12.111) "
                   r"(steps 1–3) · **N188 [B]** the log-linear wind profile $U=\frac{u_*}\kappa\big[\ln\frac z{z_o}+5\frac z{L_M}\big]$ — the book's "
                   r"coefficient 5, the stable-side $\phi_m=1+5z/L_M$ (steps 4–7)")
D("D25", ref="12.111", after=r"""
> ⚠️ **slip #17 — the book prints**, in the sentence that leads to this result, "using $\overline{uw}=u_*^2$" **; the correct form
> is** $-\overline{uw}=u_*^2$ (step 1): the momentum flux is downward, $\overline{uw}<0$. With the printed sign the definition
> $\mathrm{Rf}=\frac{-g\alpha\overline{wT'}}{-\overline{uw}(dU/dz)}$ (12.107) would give $-z/L_M$ instead of
> $\mathrm{Rf}=z/L_M$ (12.111). The cell below shows both.""")
code(r"""
us17, wT17, z17 = 0.3, -0.02, 10.0                            # friction velocity [m/s], a downward heat flux [K m/s], height [m]
L17 = ch12.monin_obukhov_length(us17, wT17, T=300.0, kappa=KAPPA)   # Eq. (12.110) [m]
shear17 = us17/(KAPPA*z17)                                    # the neutral log-layer shear dU/dz = u*/(κ z) [1/s]
Rf_right = ch12.flux_richardson(wT17, -us17**2, shear17, 1/300.0)    # Eq. (12.107) with mean(uw) = −u*²  (the correct sign)
Rf_printed = ch12.flux_richardson(wT17, +us17**2, shear17, 1/300.0)  # the same with mean(uw) = +u*²  (the sign as printed: slip #17)
print(f"z/L_M = {ch12.flux_richardson_surface_layer(z17, L17):+.4f};  Rf with mean(uw) = −u*²: {Rf_right:+.4f};  with the printed +u*²: {Rf_printed:+.4f}")
assert np.isclose(Rf_right, z17/L17) and np.isclose(Rf_printed, -z17/L17)   # the correct sign gives z/L_M, the printed one −z/L_M
""", explain=r"""
`ch12.flux_richardson(wT, uw, dUdz, alpha)` evaluated with the log-layer shear: with $\overline{uw}=-u_*^2$ it returns exactly
$z/L_M$ (the same number as `ch12.flux_richardson_surface_layer`); with the printed $+u_*^2$ it returns $-z/L_M$ — a stable night
would be reported as unstable.""")
nb.worked_example("a clear night over grass", r"""
$u_*=0.3$ m/s, $H=-30$ W/m² (downward), $\rho=1.2$ kg/m³, $c_p=1005$ J/(kg K), $T=300$ K, $\kappa=0.4$ (the round value customary
in micrometeorology; 0.41 elsewhere in this chapter — a 2.5 % difference in every wind speed here), $z_0=0.03$ m.

1. $\overline{wT'}=H/(\rho c_p)=-30/1206=-0.0249$ K m/s (primer P296).
2. $\alpha=1/T=3.33\times10^{-3}$ K⁻¹.
3. $L_M=-u_*^3/(\kappa\alpha g\overline{wT'})=-0.027/(0.4\times3.33\times10^{-3}\times9.81\times(-0.0249))=+83$ m (positive: stable).
4. Rf at 10 m $=z/L_M=0.12$.
5. Rf reaches ¼ at $z=L_M/4=21$ m: above that the turbulence cannot sustain itself.
6. Wind at 10 m: $(0.3/0.4)[\ln(333)+5\times10/83]=0.75\times(5.81+0.60)=4.81$ m/s, against 4.36 m/s for a neutral profile —
   more shear for the same stress.""")
code(r"""
L_night = ch12.monin_obukhov_from_fluxes(tau=1.2*0.3**2, H=-30.0, rho=1.2, cp=1005.0, T=300.0, kappa=0.4)   # Eq. (12.110) from the stress [Pa] and heat flux [W/m²]
print(f"L_M = {L_night:.1f} m;  Rf(10 m) = {ch12.flux_richardson_surface_layer(10.0, L_night):.3f};  "
      f"U(10 m) = {WT.surface_layer_wind(10.0, 0.3, 0.03, L_night, kappa=0.4):.2f} m/s (neutral: {WT.surface_layer_wind(10.0, 0.3, 0.03, kappa=0.4):.2f})")
st_n = ch12.surface_layer_state(0.3, -30.0, 300.0, 0.03, 10.0, 1.2, 1005.0, 0.4)   # everything at once: u*, H, T, z0, z, ρ, c_p, κ
print(f"regime: {st_n['regime']};  Rf = ¼ at z = {st_n['z_crit']:.1f} m;  verdict: {st_n['verdict']}")
print("Kundu:      ", st_n["verdict_kundu"]); print("meteorology:", st_n["verdict_met"])   # stability is read from the verdict keys
L_noon = ch12.monin_obukhov_from_fluxes(tau=1.2*0.3**2, H=150.0, rho=1.2, cp=1005.0, T=300.0, kappa=0.4)    # a sunny noon: H = +150 W/m²
print(f"noon: L_M = {L_noon:.1f} m (negative: unstable), Rf(10 m) = {ch12.flux_richardson_surface_layer(10.0, L_noon):.2f}")
print("from kinematic fluxes (κ = 0.41):", [round(float(ch12.monin_obukhov_length(0.3, wT_, T=300.0, kappa=KAPPA)), 1) for wT_ in (0.1, -0.02, 0.0)], "m")
""", explain=r"""
1. `ch12.monin_obukhov_from_fluxes(tau, H, rho, cp, T, kappa=)` converts the observer's fluxes (Pa, W/m²) to kinematic ones and
   returns $L_M=+83$ m; `ch12.flux_richardson_surface_layer` is $\mathrm{Rf}=z/L_M$ (12.111); `WT.surface_layer_wind` the log-linear
   profile: 4.81 m/s against 4.36 m/s neutral.
2. `ch12.surface_layer_state` returns all of it in one dictionary, including the height where Rf reaches ¼ (20.8 m) and the
   stability verdict in both conventions (the temperature gradient it implies is a strong inversion; in the meteorological line
   read the printed "Γ < Γa" as $\Gamma_{met}<\Gamma_d$).
3. At noon $L_M\approx-17$ m: unstable, and at 10 m buoyancy already supplies more than half as much energy as the shear.
4. `ch12.monin_obukhov_length(u_star, wT, T=, kappa=)` works from kinematic fluxes: −20 m (unstable), +101 m (stable), `inf`
   (neutral, no heat flux).""")
scratch(r"""
# From scratch: L_M in one line and the log-linear wind at 10 m
wT_n = -30.0/(1.2*1005.0)                                     # kinematic heat flux mean(wT') = H/(ρ c_p) [K m/s]
L_mine = -0.3**3/(0.4*(1/300.0)*G0*wT_n)                      # Eq. (12.110) with α = 1/T
U_mine = 0.3/0.4*(np.log(10.0/0.03) + 5*10.0/L_mine)          # the log-linear profile at z = 10 m
assert np.isclose(L_mine, L_night) and np.isclose(U_mine, WT.surface_layer_wind(10.0, 0.3, 0.03, L_night, kappa=0.4))   # same as the library
print(f"L_M = {L_mine:.2f} m, U(10 m) = {U_mine:.3f} m/s")
""", r"""
The definition and the profile, typed out, agree with `ch12.monin_obukhov_from_fluxes` and `WT.surface_layer_wind` (83.0 m with
the standard $g$; 82.98 m with the hand value 9.81).""")
note("N189 [B]", r"""
**Wind profiles for every stability.** On a logarithmic height axis the neutral profile is a straight line; a stable layer
bends toward larger wind (more shear), an unstable one toward smaller wind (the air is well mixed).""")
nb.plotly(r"""
z_ax = np.geomspace(0.05, 100.0, 80)                        # heights from just above z0 = 0.03 m to 100 m
U_neut = WT.surface_layer_wind(z_ax, 0.3, 0.03, kappa=0.4)    # the neutral logarithm for u* = 0.3 m/s

def f7(invL):                                                 # profiles for one value of 1/L_M [1/m]
    L_ = np.inf if abs(invL) < 1e-9 else 1.0/invL
    nan = np.full_like(z_ax, np.nan)
    stable = WT.surface_layer_wind(z_ax, 0.3, 0.03, L_, kappa=0.4) if invL > 0 else nan                              # log-linear (book's form)
    unstable = WT.surface_layer_wind(z_ax, 0.3, 0.03, L_, kappa=0.4, unstable="businger_dyer") if invL < 0 else nan   # Businger–Dyer form
    book_uns = WT.surface_layer_wind(z_ax, 0.3, 0.03, L_, kappa=0.4) if invL < 0 else nan                            # the book's log-linear form used outside its range
    z_q = [0.25*L_]*2 if invL > 0 and 0.25*L_ < 100 else [np.nan]*2                                                   # the height where Rf = z/L_M = ¼
    return {"neutral: dT/dz = Γa ⇔ Γ_met = Γd (logarithm)": (U_neut, z_ax),
            "stable, L_M > 0: dT/dz > Γa = −9.8 K/km ⇔ Γ_met < Γd = 9.8 K/km": (stable, z_ax),
            "unstable, L_M < 0: dT/dz < Γa = −9.8 K/km ⇔ Γ_met > Γd = 9.8 K/km (Businger–Dyer)": (unstable, z_ax),
            "unstable: the book's log-linear form (ends where 1 + 5 z/L_M = 0, z = |L_M|/5)": (book_uns, z_ax),
            "height where Rf = z/L_M = ¼ (stable side)": ([0.0, 12.0], z_q)}

figF7 = slider_figure(f7, "1/L_M", np.round(np.linspace(-0.1, 0.05, 16), 4), unit="1/m", xlabel="wind speed U [m/s]",
                      ylabel="height z [m]", title="Stable: more shear than neutral; unstable: less",
                      xrange=(0, 12))
figF7.update_yaxes(type="log", range=[np.log10(0.05), 2.0])   # logarithmic height axis
figF7.show()
""", explain=r"""
1. `f7(invL)` returns the neutral logarithm and, depending on the sign of $1/L_M$, the stable (log-linear) or unstable profile.
   For the unstable side we use a commonly used form (often called Businger–Dyer; `unstable="businger_dyer"`) and show the book's log-linear form with it. The log-linear shear
   $1+5z/L_M$ reaches zero at $z=\lvert L_M\rvert/5$; above that height it has no meaning and the function returns no value (NaN).
2. Each trace name carries the stratification in **both** lapse-rate conventions; the horizontal marker is the height
   $z=L_M/4$ where $\mathrm{Rf}=\tfrac14$.""")
nb.figure_notes(
    see="Drag $1/L_M$ from negative (unstable) through zero to positive (stable). At zero only the straight neutral line is "
        "there. To the right the profile bends toward larger wind aloft and the ¼-marker comes down from the top; to the left "
        "it bends toward smaller wind, and the book's log-linear form simply stops at $z=\\lvert L_M\\rvert/5$.",
    read="Stable: the same surface stress needs **more** shear, so the 10 m wind is stronger for a given $u_*$ — or, the way a "
         "model uses it, a given 10 m wind produces **less** stress and less heat flux. Unstable: the reverse.",
    change="…the surface were rougher (larger $z_0$): every curve shifts to smaller wind by the same amount at all heights; "
           "the *bending* depends only on $z/L_M$.")
nb.live(r"""
def surface_live(u_star=0.3, H=-30.0, log_z0=-1.5):           # free friction velocity [m/s], surface heat flux [W/m²], log10 of roughness length [m]
    z0_ = 10.0**log_z0
    L_ = ch12.monin_obukhov_from_fluxes(tau=1.2*u_star**2, H=H, rho=1.2, cp=1005.0, T=300.0, kappa=0.4)   # Eq. (12.110)
    zz = np.geomspace(1.5*z0_, 100.0, 100)                    # heights [m]
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    ax.semilogy(WT.surface_layer_wind(zz, u_star, z0_, kappa=0.4), zz, "--", color=C_REF, label="neutral")
    ax.semilogy(WT.surface_layer_wind(zz, u_star, z0_, L_, kappa=0.4, unstable="businger_dyer"), zz, color=C_BUOY, label=f"L_M = {L_:.0f} m")
    name = "stable: dT/dz > Γa ⇔ Γ_met < Γd" if L_ > 0 else "unstable: dT/dz < Γa ⇔ Γ_met > Γd"   # both conventions in the title
    ax.set(xlabel="U [m/s]", ylabel="z [m]", title=f"{name};  Rf(10 m) = {10.0/L_:+.2f}")
    ax.legend(fontsize=8)
    plt.show()

live(surface_live, u_star=(0.05, 1.0, 0.05), H=(-80.0, 400.0, 10.0), log_z0=(-4.0, 0.0, 0.25))   # sliders (kernel only)
""", explain=r"""
The same profile with three free sliders (it needs a running kernel; the slider figure above is its web-page twin). Lower $u_*$
on a night with $H<0$ and watch $\mathrm{Rf}(10\text{ m})$ climb past ¼.""")
note("N190 [B]", r"""
**A budget for temperature fluctuations.** The recipe of D10 with $T'$ in place of $u_i$ gives the budget of the temperature
variance:

$$\frac{\partial}{\partial t}\Big(\frac12\overline{T'^2}\Big)+U\frac{\partial}{\partial x}\Big(\frac12\overline{T'^2}\Big)=-\overline{wT'}\frac{d\bar T}{dz}-\frac{\partial}{\partial z}\Big(\frac12\overline{T'^2w}-\kappa\frac{\partial}{\partial z}\Big(\frac12\overline{T'^2}\Big)\Big)-\bar\varepsilon_T,\qquad\bar\varepsilon_T=\kappa\overline{(\partial T'/\partial x_j)^2}\qquad\text{(12.112)}$$

production by the heat flux acting on the mean gradient, transport, and the "thermal dissipation" $\bar\varepsilon_T$ [K²/s]
(here $\kappa$ is the thermal diffusivity, our $\kappa_{th}$).

> ⚠️ **slip #15 — the book prints** the molecular transport as $\kappa\,\partial\overline{T'^2}/\partial z$ **; the correct form is**
> $\kappa\,\partial(\tfrac12\overline{T'^2})/\partial z$ — the factor ½ belongs to every term.""")
note("N191 [B]", r"""**The temperature spectrum** $S_T(K)$ distributes the temperature variance over wavenumbers:
$\overline{T'^2}\equiv\int_0^\infty S_T(K)\,dK$.""")
note("N192 [B]", r"""
**The Obukhov–Corrsin law.** D12's dimensional argument with temperature as an extra dimension: in the inertial range the
temperature spectrum can depend only on $\bar\varepsilon_T$, $\bar\varepsilon$ and $K$.""",
     equation=r"S_T\propto\bar\varepsilon_T\bar\varepsilon^{-1/3}K^{-5/3}\quad\text{for}\quad2\pi/L\ll K\ll2\pi/\eta", ref="12.113")
note("N193 [B]", r"""
**The Batchelor scale.** When heat diffuses more slowly than momentum ($\nu/\kappa_{th}\gg1$: sea water has a Prandtl number
near 7) temperature keeps finer structure than velocity, down to $\eta_T=\eta(\kappa_{th}/\nu)^{1/2}$.""")
note("N194 [B]", r"""**Between the two scales** the smallest eddies only stretch the temperature field, and its spectrum flattens.""",
     equation=r"S_T\propto K^{-1}\quad\text{for}\quad2\pi/\eta\ll K\ll2\pi/\eta_T", ref="12.114")
code(r"""
tv = ch12.temperature_variance_sympy()                        # sympy derives the temperature-variance budget (cached)
print("derived budget = corrected form:", tv["check"], "| = the form as printed:", tv["printed_check"])   # True, False (slip #15)
z_t = np.linspace(1.0, 50.0, 50)                              # heights [m]
tb = ch12.temperature_variance_budget(z_t, 300.0 + 0.02*z_t, np.full(50, -0.02), np.full(50, 4.0e-4))   # dθ/dz = 0.02 K/m, mean(wT') = −0.02 K m/s, ε̄_T = 4e-4 K²/s
print(f"production −mean(wT') dθ/dz = {tb['production'][0]:.1e} K²/s (> 0: the flux runs down the gradient);  dissipation {tb['dissipation'][0]:.1e} K²/s")
print("exponents of S_T ∝ ε̄_T^a ε̄^b K^c:", ch12.scalar_spectrum_exponents())   # 1, −1/3, −5/3
print(f"sea water (ν = 1e-6, κ_th = 1.4e-7 m²/s, ε̄ = 1e-6 m²/s³): η = {ch12.kolmogorov_scales(1e-6, 1e-6)[0]*1e3:.2f} mm, η_T = {ch12.batchelor_scale(1e-6, 1.4e-7, 1e-6)*1e3:.3f} mm")
""", explain=r"""
1. `ch12.temperature_variance_sympy()` confirms the budget with the ½ and shows that the printed form (without it) does not follow.
2. `ch12.temperature_variance_budget(z, T_mean, wT, eps_T)`: production is positive whenever the heat flux runs down the mean
   gradient, as here on a stable night; dissipation is returned as $-\bar\varepsilon_T$.
3. `ch12.scalar_spectrum_exponents()` solves the dimensional problem of **N192** in exact fractions.
4. `ch12.batchelor_scale(nu, kappa_th, eps)`: with $\eta=1$ mm in sea water the temperature field has structure down to 0.37 mm.""")
note("N195 [B]", r"""
**Velocity and temperature spectra side by side** (our model curves for sea water).""")
fig(r"""
nu_sw, kth_sw, eps_sw = 1.0e-6, 1.4e-7, 1.0e-6                # sea water: viscosity, thermal diffusivity [m²/s]; dissipation rate [m²/s³]
eta_sw = ch12.kolmogorov_scales(nu_sw, eps_sw)[0]             # Kolmogorov length [m]
K_sw = np.logspace(0, 4.3, 300)                               # wavenumbers [rad/m]
E_sw = ch12.model_spectrum(K_sw, eps_sw, nu_sw, kind="pao")   # velocity spectrum (model of C08)
S_T = ch12.scalar_spectrum(K_sw, eps_sw, 1.0e-8, nu_sw, kth_sw, cutoff=True)   # temperature spectrum with ε̄_T = 1e-8 K²/s (illustrative)
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.loglog(K_sw*eta_sw, E_sw/E_sw[0], color=C_FLUC, label="velocity spectrum (scaled)")
ax.loglog(K_sw*eta_sw, S_T/S_T[0], color=C_BUOY, label="temperature spectrum (scaled)")
ax.loglog(K_sw*eta_sw, (K_sw/K_sw[0])**(-5/3), ":", color=C_REF, label="$K^{-5/3}$")
ax.axvline(1.0, color=C_VISC, ls="--"); ax.axvline(eta_sw/ch12.batchelor_scale(nu_sw, kth_sw, eps_sw), color=C_BUOY, ls="--")
ax.set(xlabel="K η [–]", ylabel="spectrum / its value at the left edge [–]", ylim=(1e-9, 3), title="Temperature keeps finer structure than velocity when ν/κ ≫ 1")
ax.legend(fontsize=8)
""",
    see="Two curves that follow the dotted −5/3 line together up to $K\\eta\\approx1$ (rose dashed line). There the teal "
        "velocity spectrum falls off steeply, while the blue temperature spectrum turns to the flatter slope −1 and continues to "
        "the Batchelor scale (blue dashed line) before it too is cut off.",
    read="Velocity fluctuations are smoothed out at $\\eta$; temperature fluctuations survive to $\\eta_T<\\eta$ because heat "
         "diffuses about seven times more slowly than momentum in sea water. Ocean microstructure probes exploit this range to "
         "measure mixing.",
    change="…this were air ($\\nu/\\kappa_{th}\\approx0.7$): there would be no −1 range — both spectra would end near $\\eta$.")
explainer("stratified_surface_layer", "When does stratification kill turbulence?",
          "dragging the surface heat flux through zero flips the sign of the Monin–Obukhov length, bends the wind profile to the "
          "other side of the neutral logarithm and moves the Rf = ¼ height up and down, while the badge states the regime in both "
          "lapse-rate conventions.",
          ["Load 'clear calm night' and lower u_*: the Rf = ¼ height drops toward the ground.",
           "Drag H through zero: L_M jumps from +∞ to −∞ and the profile crosses the neutral line.",
           "Toggle the convention: the same layer, the other inequality.",
           "Switch z₀ from sea to forest at the same 10 m wind."])
whatif(r"""…instead of momentum and heat we follow a puff of smoke released into the turbulence? C16.""")


# =====================================================================================================================
# A.12 §12.12 Taylor's Theory of Turbulent Dispersion — R07 R08, C16
# =====================================================================================================================
nb.section("12.12", "Taylor's Theory of Turbulent Dispersion", intro=r"""
**What is this section about?** How far turbulence carries a marked particle. One exact formula shows the cloud spreads in
proportion to time at first and to the square root of time later — and that an "eddy diffusivity" is not a constant.""")
nb.recap("R07", "The molecular yardstick", r"""
A line vortex switched on at $t=0$ spreads by viscosity as $u_\theta=(\Gamma/2\pi r)\exp(-r^2/4\nu t)$: a Gaussian whose
standard deviation per coordinate is $\sigma=\sqrt{2\nu t}$, i.e. $\sigma^2=2\nu t$. (⚠️ Ch. 3's vortex *core radius* used
$4\nu t$ — a different definition of width.) Molecular diffusion always spreads like $\sqrt t$.""", where="Ch. 8 §8.4; Ch. 3 and Ch. 5 (Lamb–Oseen vortex)")
nb.recap("R08", "A diffusivity from the growth of a variance", r"""
Ch. 1's $\sigma^2=2Dt$ read backwards gives a way to *measure* a diffusivity: $\nu=\frac12\frac{d\sigma^2}{dt}$ (12.126)
(`ch12.diffusivity_from_variance`). This is the starting line of the derivation D28 below.""", where="Ch. 1 §1.5")
core("C16", r"Taylor's formula $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\,t\int_0^t\big(1-\frac\tau t\big)r_\alpha(\tau)\,d\tau$ (12.119)",
     "Why does a puff spread like t at first and like √t later?")
problem(r"""
Smoke leaves a chimney as a narrow cone that opens like a wedge, then, far downwind, widens ever more slowly. An oil slick,
volcanic ash, a cloud of drifting buoys in the ocean behave the same way. The switch happens when each smoke particle has
forgotten the velocity it started with.""")
idea(r"""
t ≪ Λ_t : the particle still has its first velocity  → flies straight:  X ≈ u t          → X_rms = u_rms · t
t ≫ Λ_t : many independent "steps" of duration ~2Λ_t  → random walk:    X_rms = step·√n  → X_rms = u_rms √(2Λ_t t)
""", r"""
Symbols: $X_\alpha(t)$ is one coordinate of a particle's displacement from its release point [m]; $u_\alpha(t)$ the same
component of **that particle's** velocity; $r_\alpha(\tau)$ its autocorrelation coefficient; $\Lambda_t$ its integral (memory)
time. ⚠️ $\alpha$ is a component label here (no sum over it), not the expansion coefficient.""")
remind("C16")
note("N196 [B]", r"""
**The set-up.** We follow particles (the Lagrangian description of Ch. 3): $\mathbf X(\mathbf a,t)$ is the position at time $t$ of
the particle that started at $\mathbf a$. The turbulence is stationary and homogeneous with zero mean velocity.
`ch12.langevin_particles` makes such paths by integrating the Ornstein–Uhlenbeck velocity of C01 (primer P283) in time.""")
note("N199 [B]", r"""
**The Lagrangian autocorrelation** is C02's coefficient for the velocity of one particle:
$r_\alpha(\tau)\equiv\overline{u_\alpha(t)u_\alpha(t+\tau)}/\overline{u_\alpha^2}$; its integral $\Lambda_t$ is the **Lagrangian**
integral time scale (not the one a fixed probe measures).""")
fig(r"""
t_p = np.linspace(0.0, 100.0, 1001)                           # times [s]
Xp, up = ch12.langevin_particles(12, t_p, 1.0, 10.0, seed=3)  # 12 particles: u_rms = 1 m/s, memory Λ_t = 10 s → positions [m], velocities [m/s]
fig, ax = plt.subplots(figsize=(6.6, 3.5))
ax.plot(t_p, Xp.T, lw=0.8, color=C_FLUC, alpha=0.7)           # each particle's path X(t)
ax.plot(t_p, np.sqrt(ch12.taylor_dispersion_exponential(t_p, 1.0, 10.0)), color=C_MEAN, lw=2, label="$+X_{rms}(t)$ (Taylor's formula)")
ax.plot(t_p, -np.sqrt(ch12.taylor_dispersion_exponential(t_p, 1.0, 10.0)), color=C_MEAN, lw=2)
ax.axvline(10.0, color=C_REF, ls=":")
ax.set(xlabel="time t [s]", ylabel="displacement X [m]", title="Paths from a point source: straight at first, wandering later")
ax.legend(fontsize=8)
""",
    see="Twelve teal paths leaving the origin. For the first few seconds each is nearly a straight line with its own slope; after "
        "the memory time (dotted line at 10 s) they wander. The purple envelope $\\pm X_{rms}(t)$ opens like a wedge and then "
        "like a parabola on its side.",
    read="Early on a particle's displacement is (its initial velocity) × (time); later it is a sum of many unrelated pieces. The "
         "envelope is what this block derives.",
    change="…the memory were longer: the straight part of every path, and the wedge-shaped part of the envelope, would last "
           "longer.")
P("P305", "a double integral over a triangle", r"""
$\int_0^tdt'\int_0^{t'}g(\tau)\,d\tau$ covers the triangle $0<\tau<t'<t$. Each value of $\tau$ is counted for every $t'$ between
$\tau$ and $t$ — a strip of length $t-\tau$ — so the double integral equals the single integral $\int_0^t(t-\tau)g(\tau)\,d\tau$.""", code=r"""
from scipy.integrate import dblquad                          # a double integral over a region with variable inner limits
g = lambda tau: np.exp(-tau)                                  # any function of the inner variable
double = dblquad(lambda tau, tp: g(tau), 0, 2.0, lambda tp: 0.0, lambda tp: tp)[0]   # ∫_0^2 dt' ∫_0^t' g dτ
single = quad(lambda tau: (2.0 - tau)*g(tau), 0, 2.0)[0]      # ∫_0^2 (2 − τ) g dτ
print(round(double, 4), round(single, 4))                     # 1.1353 both ways
""")
NOTES_IN["D26"] = (r"**N197 [B]** $\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha\frac{dX_\alpha}{dt}}$ (12.115) (steps 2–3) · **N198 [B]** "
                   r"$\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha u_\alpha}=2\int_0^t\overline{u_\alpha(t')u_\alpha(t)}\,dt'$ (12.116) (steps 4–5) · "
                   r"**N200 [B]** the rate of spreading is the integral of the autocorrelation, $\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{u_\alpha^2}\int_0^tr_\alpha(\tau)\,d\tau$ (12.117) "
                   r"(steps 6–7) · **N201 [B]** the double integral $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\int_0^tdt'\int_0^{t'}r_\alpha(\tau)\,d\tau$ (12.118) (step 8) · "
                   r"**N202 [B]** the integration by parts that turns it into the single integral of this block's title (steps 9–11)")
PF["D26"]["title"] = "Taylor's formula: from the particle's velocity to the double integral and the single integral"
D("D26", ref="12.119")
NOTES_IN["D27"] = (r"**N204 [B]** short times, $\overline{X_\alpha^2}\simeq\overline{u_\alpha^2}t^2$ (12.120) (steps 1–2) · **N205 [B]** "
                   r"$(X_\alpha)_{rms}=(u_\alpha)_{rms}t$ for $t\ll\Lambda_t$ (12.121) (step 2) · **N206 [B]** long times, $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\Lambda_tt$ (12.122) "
                   r"(steps 4–5) · **N207 [B]** $(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ for $t\gg\Lambda_t$ (12.123) (step 5) · **N208 [B]** the closed "
                   r"form for an exponential memory $r=e^{-\tau/\Lambda_t}$ (steps 6–8); for a Gaussian memory $r=e^{-\tau^2/t_c^2}$ "
                   r"($\Lambda_t=\tfrac{\sqrt\pi}2t_c$) the closed form is `ch12.taylor_dispersion_gaussian`")
D("D27", ref="12.121", after=r"""
> ⚠️ **slip #11 — the book prints** a reference to "(11.119)" in this argument **; the correct form is** this chapter's
> $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)r_\alpha\,d\tau$ (12.119).""")
note("N203 [B]", r"""
**"Small $t$" and "large $t$" on the correlation curve.** The two limits are two ways of looking at the same integral: for
$t\ll\Lambda_t$ only the top of the curve, $r\approx1$, is inside the range of integration; for $t\gg\Lambda_t$ the whole curve
is, and its area is $\Lambda_t$.""")
fig(r"""
tau_c16 = np.linspace(0.0, 60.0, 400)                         # lags [s]
r_c16 = np.exp(-tau_c16/10.0)                                 # the Lagrangian autocorrelation with Λ_t = 10 s
fig, ax = plt.subplots(figsize=(6.2, 3.2))
ax.plot(tau_c16, r_c16, color=C_FLUC, label="$r_\\alpha(\\tau)=e^{-\\tau/\\Lambda_t}$")
ax.fill_between(tau_c16[tau_c16 <= 2], r_c16[tau_c16 <= 2], color=C_RS, alpha=0.4, label="small t = 2 s: r ≈ 1 over the whole range")
ax.fill_between(tau_c16[tau_c16 <= 50], r_c16[tau_c16 <= 50], color=C_MEAN, alpha=0.12, label="large t = 50 s: the area is ≈ Λ_t")
ax.set(xlabel="lag τ [s]", ylabel="r [–]", title="The two limits of Taylor's formula")
ax.legend(fontsize=8)
""",
    see="The exponential correlation with two shaded ranges of integration: a narrow orange strip under the flat top (small $t$) "
        "and a wide purple area that contains nearly the whole curve (large $t$).",
    read="Small $t$: $\\int_0^tr\\,d\\tau\\approx t$, so the spreading *rate* grows with time — ballistic. Large $t$: "
         "$\\int_0^tr\\,d\\tau\\approx\\Lambda_t$, a constant rate — diffusive.",
    change="…the correlation had a negative lobe (particles that tend to come back, as in a wave field): the area, and with it "
           "the long-time spreading, would be reduced.")
nb.worked_example("a puff in turbulence with a 10 s memory", r"""
$u_{rms}=1$ m/s, $\Lambda_t=10$ s, exponential correlation: $\overline{X^2}=2\overline{u^2}\Lambda_t^2[t/\Lambda_t-1+e^{-t/\Lambda_t}]$.

1. $t=1$ s: $200\times(0.1-1+0.90484)=0.967$ m² — the ballistic value $u^2t^2=1.000$ is 3 % high: $X_{rms}\approx0.98$ m.
2. $t=10$ s: $200\times0.3679=73.6$ m² ($X_{rms}=8.6$ m) — neither limit works (ballistic 100, diffusive 200).
3. $t=100$ s: $200\times9.000=1800$ m² ($X_{rms}=42$ m); the diffusive formula $2u^2\Lambda_tt=2000$ is 11 % high; with the
   constant offset $-2u^2\Lambda_t^2=-200$ it is exact.
4. Eddy diffusivity at $t=10$ s: $D_T=u^2\Lambda_t(1-e^{-1})=6.3$ m²/s; final value 10 m²/s.""")
code(r"""
n_part = 2000 if FAST else 10000                              # number of particles
t_l = np.linspace(0.0, 100.0, 401 if FAST else 1001)          # times [s]
X, u = ch12.langevin_particles(n_part, t_l, 1.0, 10.0, seed=0)   # positions [m] and velocities [m/s], particles on axis 0
X2 = (X**2).mean(axis=0)                                      # mean square displacement over the particles
se = (X**2).std(axis=0)/np.sqrt(n_part)                       # its standard error
exact = ch12.taylor_dispersion_exponential(t_l, 1.0, 10.0)    # the closed form for an exponential memory
for t_ in (1.0, 10.0, 100.0):
    i = np.argmin(np.abs(t_l - t_))
    print(f"t = {t_:5.0f} s ({ch12.dispersion_regime(t_, 10.0):10s}): closed form {exact[i]:8.2f} m², particles {X2[i]:8.2f} ± {se[i]:.2f} m²")
assert np.all(np.abs(X2 - exact)[1:] < 5*se[1:] + 1e-9)       # particles within 5 standard errors at every time
tt3 = np.array([1.0, 10.0, 100.0])                            # three times for the two integral forms
single = ch12.taylor_dispersion(tt3, lambda s: np.exp(-s/10.0), 1.0, form="single")   # Eq. (12.119)
double = ch12.taylor_dispersion(tt3, lambda s: np.exp(-s/10.0), 1.0, form="double")   # Eq. (12.118)
print("single-integral form:", np.round(single, 3), "| double-integral form equal:", np.allclose(single, double, rtol=1e-8))
lhs, rhs = ch12.dispersion_rate_from_particles(t_l, X, u)     # both sides of d mean(X²)/dt = 2 mean(X u), Eq. (12.116)
print(f"rate of spreading at t = 50 s: d mean(X²)/dt = {np.interp(50, t_l, lhs):.1f}, 2 mean(Xu) = {np.interp(50, t_l, rhs):.1f}, theory {ch12.taylor_dispersion_rate(50.0, lambda s: np.exp(-s/10.0), 1.0):.1f} m²/s")
print("rate at t = 1, 10, 100 s:", np.round(ch12.taylor_dispersion_rate(tt3, lambda s: np.exp(-s/10.0), 1.0), 2), "m²/s")
""", explain=r"""
1. `ch12.langevin_particles(n, t, u_rms, Lambda_t)` returns positions and velocities; averaging $X^2$ over particles estimates
   $\overline{X_\alpha^2}$, with a standard error that the `assert` uses (5 standard errors, primer P281).
2. `ch12.taylor_dispersion_exponential` is the closed form of the tiny example: 0.967, 73.6, 1800 m². `ch12.dispersion_regime`
   names the regime (ballistic for $t\le0.3\Lambda_t$, diffusive for $t\ge3\Lambda_t$, transition between).
3. `ch12.taylor_dispersion(t, r_fn, u2, form=…)` evaluates Taylor's formula for **any** correlation function, as the single integral
   or as the double integral of D26: the same numbers.
4. `ch12.dispersion_rate_from_particles` checks the first step of the derivation on the particles: the rate of growth of
   $\overline{X^2}$ equals $2\overline{Xu}$; `ch12.taylor_dispersion_rate` is $2\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$ (12.117): it rises
   from $2u^2t$ (1.9 m²/s at 1 s) to $2u^2\Lambda_t=20$ m²/s.""")
scratch(r"""
# From scratch: Taylor's integral by the trapezoid rule, and a random walk by a cumulative sum
mine = []
for t_ in (1.0, 10.0, 100.0):
    tau = np.linspace(0.0, t_, 20001)                         # lags from 0 to t
    mine.append(2*1.0*t_*np.trapezoid((1 - tau/t_)*np.exp(-tau/10.0), tau))   # 2 u² t ∫ (1 − τ/t) r dτ
assert np.allclose(mine, ch12.taylor_dispersion_exponential(np.array([1.0, 10.0, 100.0]), 1.0, 10.0), rtol=1e-5)   # same as the closed form
rng = np.random.default_rng(5)                                # seeded generator
theta = rng.uniform(0, 2*np.pi, (4000, 100))                  # 4000 walkers × 100 steps: a random direction for every step
R = np.cumsum(np.stack([np.cos(theta), np.sin(theta)], axis=-1), axis=1)   # positions after each unit step (cumulative sum of steps)
R2 = (R[:, -1]**2).sum(axis=1)                                # squared distance from the start after 100 steps
lib = ch12.random_walk(100, 4000, seed=0)                     # the library's walkers: array (walker, step, component)
R2_lib = (lib[:, -1]**2).sum(axis=1)
se_rw = R2.std()/np.sqrt(4000)                                # standard error of the mean of R²
assert abs(R2.mean() - 100.0) < 5*se_rw and abs(R2_lib.mean() - 100.0) < 5*se_rw   # mean(R_n²) = n L² within 5 standard errors
print("Taylor integral:", np.round(mine, 3), f"| random walk: rms distance after 100 steps = {np.sqrt(R2.mean()):.2f} (mine), {np.sqrt(R2_lib.mean()):.2f} (library); L√n = 10")
""", r"""
1. The single integral with `np.trapezoid` reproduces the closed form at the three times.
2. A random walk is a cumulative sum (`np.cumsum`) of unit steps in random directions; after 100 steps the rms distance is
   $L\sqrt n=10$ step lengths, for our walkers and for `ch12.random_walk` alike.""")
nb.animation(r"""
nfr = 30 if FAST else 60                                      # frames
Xa, ua = ch12.langevin_particles(400, t_l, 1.0, 10.0, seed=1) # 400 particles for the movie
idx = np.unique(np.round(np.geomspace(1, t_l.size - 1, nfr)).astype(int))   # frame times, log-spaced so that the early stage is visible
rng_y = np.random.default_rng(2).uniform(-1, 1, 400)          # a random vertical position for each dot (display only)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 3.5))
dots = a1.scatter(Xa[:, 0], rng_y, s=4, color=C_FLUC)         # the particles
edge_l = a1.axvline(0.0, color=C_MEAN, lw=1.5); edge_r = a1.axvline(0.0, color=C_MEAN, lw=1.5)   # the two edges at ± X_rms
a1.set(xlim=(-130, 130), ylim=(-1, 1), yticks=[], xlabel="displacement X [m]")
a2.loglog(t_l[1:], exact[1:], color=C_MEAN, label="Taylor (closed form)")
a2.loglog(t_l[1:], t_l[1:]**2, ":", color=C_RS, label="ballistic $u^2t^2$ (slope 2)")
a2.loglog(t_l[1:], 20*t_l[1:], "--", color=C_REF, label="diffusive $2u^2\\Lambda_tt$ (slope 1)")
(pt,) = a2.plot([], [], "o", color=C_FLUC, label="400 particles")
a2.set(xlabel="time t [s]", ylabel="mean(X²) [m²]", xlim=(0.1, 100), ylim=(0.01, 3000))
a2.legend(fontsize=7)

def update(k):                                                # frame k: the cloud at time t_l[idx[k]]
    i = idx[min(k, idx.size - 1)]
    dots.set_offsets(np.c_[Xa[:, i], rng_y])
    xr = np.sqrt(np.mean(Xa[:, i]**2))                        # the cloud's rms half-width now
    edge_l.set_xdata([-xr, -xr]); edge_r.set_xdata([xr, xr]) # move the two edge lines
    pt.set_data(t_l[idx[:k + 1]], np.mean(Xa[:, idx[:k + 1]]**2, axis=0))
    a1.set_title(f"t = {t_l[i]:.1f} s ({ch12.dispersion_regime(t_l[i], 10.0)}): X_rms = {xr:.1f} m", fontsize=10)
    return dots, pt

show_animation(animate(update, frames=idx.size, fig=fig, interval=120), player="video", dpi=60)   # a smooth video
""", explain=r"""
1. Left: 400 particles released at one point (the vertical position is only for display), with two purple lines at $\pm X_{rms}$.
2. Right: their mean square displacement, dot by dot, against the closed form and the two limits on log–log axes (slope 2, then 1).
3. The frame times are log-spaced, so the ballistic stage is not over in the first frame.""")
nb.figure_notes(
    see="The cloud opens quickly at first — its width grows in proportion to time — and then ever more slowly. On the right the "
        "dots climb along the dotted slope-2 line, bend near $t\\approx\\Lambda_t=10$ s, and settle parallel to the dashed slope-1 line "
        "(a little below it: the constant offset of D27).",
    read="One formula covers both regimes; the bend is at the particles' memory time. The title names the regime "
         "(`ch12.dispersion_regime`).",
    change="…$\\Lambda_t=100$ s: the cloud would still be a wedge at the end of the movie — and ten times wider than a "
           "$\\sqrt t$ law fitted to late times would predict for early times.")
P("P306", "proof by induction", r"""
To prove a statement for every whole number $n$: show it for $n=1$; show that if it holds for $n-1$ it holds for $n$; then it
holds for every $n$ (each case hands over to the next, like dominoes).""", code=r"""
for n in range(1, 6):                             # the statement "1 + 2 + … + n = n(n+1)/2", checked for n = 1 … 5
    print(n, sum(range(1, n + 1)), n*(n + 1)//2)  # the induction step adds n to both sides: (n−1)n/2 + n = n(n+1)/2
""")
note("N209 [B]", r"""
**The random walk.** A walker takes steps of length $L$ in random directions: $\mathbf R_n=\mathbf R_{n-1}+\mathbf L$. Squaring
and averaging, the cross term vanishes when the new step is uncorrelated with the past — the product rule of C01 again.""",
     equation=r"\overline{R_n^2}=\overline{R_{n-1}^2}+L^2+2\overline{\mathbf R_{n-1}\cdot\mathbf L}", ref="12.124")
note("N210 [B]", r"""**By induction** $\overline{R_n^2}=\overline{R_{n-1}^2}+L^2=\dots=nL^2$: every step adds $L^2$ to the mean
square distance (the figure shows three walkers and the ensemble rms).""")
note("N211 [B]", r"""
**The rms distance grows like $\sqrt n$.** With $n=t/\Delta t$ steps of length $L=u_{rms}\Delta t$ and $\Delta t=2\Lambda_t$ this
is exactly $(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ (12.123): a dispersing particle makes one independent step every
two memory times.""", equation=r"(R_n)_{rms}=L\sqrt n", ref="12.125")
fig(r"""
walk = ch12.random_walk(100, 4000, seed=0)                    # 4000 walkers, 100 unit steps in two dimensions
rms_n = np.sqrt((walk**2).sum(axis=2).mean(axis=0))           # ensemble rms distance after each step
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 3.4))
for k, col in zip(range(3), (C_FLUC, C_RS, C_BUOY)):
    a1.plot(walk[k, :, 0], walk[k, :, 1], lw=0.9, color=col)  # three individual paths
a1.add_patch(plt.Circle((0, 0), 10.0, fill=False, color=C_MEAN, ls="--"))   # the circle of radius L√n = 10
a1.set(aspect="equal", xlim=(-22, 22), ylim=(-22, 22), xlabel="x / L", ylabel="y / L", title="three walkers, 100 steps")
a2.loglog(np.arange(1, 101), rms_n[1:], "o", ms=3, color=C_MEAN, label="rms distance of 4000 walkers")
a2.loglog(np.arange(1, 101), np.sqrt(np.arange(1, 101)), color=C_REF, label="$L\\sqrt{n}$")
a2.set(xlabel="number of steps n [–]", ylabel="$(R_n)_{rms}/L$ [–]", title=f"after 100 steps: {rms_n[-1]:.2f} L")
a2.legend(fontsize=8)
""",
    see="Left: three tangled paths of 100 steps each and a dashed circle of radius 10. Right: the rms distance of 4000 walkers "
        "on a straight line of slope ½ on log–log axes.",
    read="No single walker follows the $\\sqrt n$ law — some end inside the circle, some outside; the law is for the "
         "**ensemble**. One hundred steps take you only ten step lengths away: random steps mostly cancel.",
    change="…successive steps were correlated (`persistence=0.5`): the walk would still spread like $\\sqrt n$, but as if each "
           "step were longer — which is what a long velocity memory does to dispersion.")
note("N212 [B]", r"""
**A smoke plume in a wind** is the same picture with time read as distance, $t=x/U$: the time-averaged plume widens linearly
near the source and like $\sqrt x$ far away.

> ⚠️ **slip #13 — the book's figure caption prints** width $\propto x^{1/2}$ near the source and $\propto x$ far away **; the
> correct statement is** $\propto x$ near (from $(X_\alpha)_{rms}=(u_\alpha)_{rms}t$ (12.121)) and $\propto x^{1/2}$ far (from
> $(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ (12.123)).""")
fig(r"""
x_pl = np.linspace(1.0, 5000.0, 400)                          # distance downwind [m]
z_pl = np.linspace(-150.0, 150.0, 201)                        # cross-wind coordinate [m]
Zrms = ch12.smoke_plume_width(x_pl, 5.0, 0.5, 20.0)           # rms half-width for wind 5 m/s, w_rms = 0.5 m/s, Λ_t = 20 s
conc = ch12.plume_concentration(x_pl[None, :], z_pl[:, None], 1.0, 5.0, Zrms[None, :])   # Gaussian cross-section, source strength 1 kg/(m s)
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.6), gridspec_kw=dict(width_ratios=[1.5, 1]))
x_log = np.geomspace(1.0, 5000.0, 200)                        # the same distances on a logarithmic grid
ax2.loglog(x_log, ch12.smoke_plume_width(x_log, 5.0, 0.5, 20.0), color=C_RS, label="$Z_{rms}(x)$")
ax2.loglog(x_log[x_log < 60], 0.1*x_log[x_log < 60], ":", color=C_REF, label="∝ x (wedge)")
ax2.loglog(x_log[x_log > 300], 0.5*np.sqrt(8.0*x_log[x_log > 300]), "--", color=C_REF, label="∝ √x (parabola)")
ax2.axvline(100.0, color=C_REF, lw=0.6)                       # x = U Λ_t = 100 m, where the particles have forgotten their first velocity
ax2.set(xlabel="distance downwind x [m]", ylabel="half-width [m]", title="the same on log–log axes")
ax2.legend(fontsize=7)
pc = ax.pcolormesh(x_pl, z_pl, np.log10(np.maximum(conc, 1e-6)), cmap="Purples", shading="auto", vmin=-4.5, vmax=-1.5)
ax.plot(x_pl, Zrms, color=C_RS, label="$\\pm Z_{rms}(x)$"); ax.plot(x_pl, -Zrms, color=C_RS)
ax.plot(x_pl[x_pl < 800], 0.5*x_pl[x_pl < 800]/5.0, ":", color=C_REF, label="near: ∝ x")
ax.plot(x_pl, 0.5*np.sqrt(2*20.0*x_pl/5.0), "--", color=C_REF, label="far: ∝ √x")
ax.axvline(5.0*20.0, color=C_REF, lw=0.6)
fig.colorbar(pc, label="log10 concentration [kg/m³]")
ax.set(xlabel="distance downwind x [m]", ylabel="cross-wind z [m]", ylim=(-150, 150), title="A time-averaged plume: a wedge near the source, a parabola far away")
ax.legend(fontsize=7, loc="upper left")
print("Z_rms at x = 10, 100, 1000, 5000 m:", np.round(ch12.smoke_plume_width(np.array([10.0, 100.0, 1000.0, 5000.0]), 5.0, 0.5, 20.0), 2), "m")
""",
    see="A plume (colour = concentration) spreading to the right, with its rms half-width in orange. Up to about "
        "$x=U\\Lambda_t=100$ m (thin vertical line) the edge follows the dotted straight line; far downstream it follows the dashed "
        "$\\sqrt x$ curve. Right: the half-width on log–log axes, where the near-source wedge (slope 1, the first 100 m — only 2 % of "
        "the left panel) and the far-field parabola (slope ½) are both plain to see.",
    read="Read $x$ as time: $t=x/U$. The half-width is about 1 m at 10 m, 8.6 m at 100 m, 42 m at 1 km and 99 m at 5 km — "
         "tenfold in distance gives nearly tenfold in width near the source and only $\\sqrt{10}\\approx3$-fold far away.",
    change="…the atmosphere were stable (a night-time inversion): $w_{rms}$ and the vertical memory time fall, and the plume "
           "stays a thin ribbon for kilometres.")
NOTES_IN["D28"] = (r"**N213 [B]** the eddy diffusivity $D_T\equiv\frac12\frac{d}{dt}(\overline{X_\alpha^2})=\overline{u_\alpha^2}\int_0^tr_\alpha(\tau)\,d\tau$ (12.127) "
                   r"(steps 3–4) · **N214 [B]** $D_T\cong\overline{u_\alpha^2}t$ for $t\ll\Lambda_t$ (12.128) (step 5) · **N215 [B]** "
                   r"$D_T\cong\overline{u_\alpha^2}\Lambda_t$ for $t\gg\Lambda_t$ (12.129) (step 6)")
PF["D28"]["title"] = r"The eddy diffusivity: $\nu=\tfrac12\,d\sigma^2/dt$, $D_T=\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$ and its two limits (the printed condition corrected)"
D("D28", ref="12.127", after=r"""
> ⚠️ **slip #12 — the book prints** the condition $t\ll\Lambda_t$ on the constant value $D_T\cong\overline{u_\alpha^2}\Lambda_t$ (12.129)
> **; the correct form is** $t\gg\Lambda_t$: the constant is the long-time limit (at $t=0.01\Lambda_t$ the constant would be 100
> times the true diffusivity).""")
code(r"""
D_num = ch12.eddy_diffusivity_taylor(tt3, lambda s: np.exp(-s/10.0), 1.0)     # Eq. (12.127) for the exponential memory at t = 1, 10, 100 s
print("D_T at t = 1, 10, 100 s:", np.round(D_num, 3), "m²/s = closed form", np.round(ch12.eddy_diffusivity_exponential(tt3, 1.0, 10.0), 3))
print("short-time form u²t at t = 1 s:", ch12.eddy_diffusivity_asymptote(1.0, 1.0, 10.0, which="short"), "| long-time form u²Λ_t at t = 100 s:", ch12.eddy_diffusivity_asymptote(100.0, 1.0, 10.0, which="long"))
print("the constant under the PRINTED condition, at t = 0.1 s:", ch12.eddy_diffusivity_asymptote(0.1, 1.0, 10.0, which="long", printed=True),
      "m²/s — the true value there is", round(float(ch12.eddy_diffusivity_exponential(0.1, 1.0, 10.0)), 4))
D_part = ch12.diffusivity_from_variance(t_l, X2)              # Eq. (12.126) applied to the particles: ½ d mean(X²)/dt
print(f"from the particles at t = 100 s: D_T ≈ {D_part[-5:].mean():.1f} m²/s (theory {ch12.eddy_diffusivity_exponential(100.0, 1.0, 10.0):.1f})")
""", explain=r"""
1. `ch12.eddy_diffusivity_taylor(t, r_fn, u2)` is half the spreading rate: 0.95, 6.3 and 10.0 m²/s — it **grows** from $u^2t$
   toward $u^2\Lambda_t$ (step 4 of the tiny example).
2. `ch12.eddy_diffusivity_asymptote(..., which=…)` returns each limiting form where its condition holds; `printed=True` applies
   the constant under the printed (wrong) condition: 10 m²/s at $t=0.1$ s, a hundred times the true 0.1 m²/s (slip #12).
3. `ch12.diffusivity_from_variance` (recap R08) applied to the particle cloud returns about 10 m²/s at late times.""")
nb.plotly(r"""
t_f8 = np.geomspace(0.1, 1000.0, 100)                        # times [s]

def f8(Lam):                                                  # curves for one memory time Λ_t [s] (u_rms = 1 m/s)
    Xr = np.sqrt(ch12.taylor_dispersion_exponential(t_f8, 1.0, Lam))   # exact X_rms(t)
    return {"X_rms (Taylor)": (t_f8, Xr), "ballistic u_rms t": (t_f8, 1.0*t_f8),
            "diffusive u_rms √(2 Λ_t t)  (= a constant diffusivity u²Λ_t from the start)": (t_f8, np.sqrt(2*Lam*t_f8))}

figF8 = slider_figure(f8, "Λ_t", np.round(np.geomspace(1.0, 100.0, 13), 2), unit="s", xlabel="time t [s]",
                      ylabel="X_rms [m]", title="The bend from t to √t sits at the memory time")
figF8.update_xaxes(type="log", range=[-1, 3])                 # logarithmic time axis
figF8.update_yaxes(type="log", range=[-1.5, 3.2])             # logarithmic vertical axis
figF8.show()

def f8d(Lam):                                                 # the eddy diffusivity for the same slider values, on its own axes
    return {"eddy diffusivity D_T(t)": (t_f8, ch12.eddy_diffusivity_exponential(t_f8, 1.0, Lam)),
            "short-time form u² t": (t_f8[t_f8 < Lam], 1.0*t_f8[t_f8 < Lam]), "long-time value u² Λ_t": ([Lam, 1000.0], [Lam, Lam])}

figF8d = slider_figure(f8d, "Λ_t", np.round(np.geomspace(1.0, 100.0, 13), 2), unit="s", xlabel="time t [s]",
                       ylabel="D_T [m²/s]", title="The eddy diffusivity grows, then levels off at u²Λ_t")
figF8d.update_xaxes(type="log", range=[-1, 3])                # logarithmic time axis
figF8d.update_yaxes(type="log", range=[-1.2, 2.3])            # logarithmic diffusivity axis
figF8d.show()
""", explain=r"""
1. `f8(Lam)` returns the exact rms displacement and its two limits; the diffusive line is also what a **constant** diffusivity
   $D=u^2\Lambda_t$ would give from the start — the "Fickian ghost".
2. `f8d(Lam)` puts the eddy diffusivity $D_T(t)$ on its own axes (its unit is m²/s, not m), with its short-time form and its
   long-time value.""")
nb.figure_notes(
    see="Drag $\\Lambda_t$: the purple curve follows the ballistic line up to $t\\approx\\Lambda_t$ and the diffusive line after; "
        "the bend moves with the slider. In the second figure the diffusivity rises linearly and levels off at $u^2\\Lambda_t$ — ten "
        "times higher for ten times the memory.",
    read="Near the source the constant-diffusivity (diffusive) line lies far **above** the true width: Fickian diffusion spreads "
         "a fresh puff much too fast, because it assumes the particles have already forgotten their velocities.",
    change="…the turbulence were twice as intense at the same memory time: every curve shifts up by a factor 2 ($X_{rms}$) or 4 "
           "($D_T$); the bend stays at $t\\approx\\Lambda_t$.")
note("N216 [B]", r"""
**Turbulent diffusion is not molecular diffusion with a bigger constant.** Near the source a constant-$D$ model spreads a puff
like $\sqrt t$ — far too wide a start — and a patch that keeps growing is stirred by ever larger eddies, so its effective
diffusivity keeps growing too (Ch. 10 drew the same line between a numerical and an eddy diffusivity).""")
note("N06 [B]", r"""
**Richardson's four-thirds law**, stated here where both ingredients exist. In the inertial range (C08) the effective
diffusivity $K$ [m²/s] for a patch of size $l$ can depend only on $l$ and $\bar\varepsilon$; dimensions leave one combination,
$K\sim\bar\varepsilon^{1/3}l^{4/3}$. **Climate hook:** this is why tracer diffusivities in ocean models are chosen to depend on
the grid size.""")
code(r"""
print("dimensionless group:", DIM.solve_exponents("K", ["eps", "l"], {"K": {"L": 2, "T": -1}, "eps": {"L": 2, "T": -3}, "l": {"L": 1}}))   # K ε̄^(−1/3) l^(−4/3)
print(f"1 km patch in the ocean (ε̄ = 1e-8 m²/s³): K ≈ {ch12.richardson_diffusivity(1000.0, 1e-8):.0f} m²/s;  in the atmosphere (ε̄ = 1e-3): K ≈ {ch12.richardson_diffusivity(1000.0, 1e-3):.0f} m²/s")
""", explain=r"""
Chapter 1's Π-theorem solver finds the one dimensionless group, i.e. $K\propto\bar\varepsilon^{1/3}l^{4/3}$;
`ch12.richardson_diffusivity(l, eps)` evaluates it with a constant of order one: about 22 m²/s for a 1 km patch in the ocean and
1000 m²/s in the atmosphere.""")
explainer("taylor_dispersion", "Why t first and √t later?",
          "the particles spread on one clock with the mean-square-displacement curve changing slope from 2 to 1 as t passes the "
          "memory time; dragging the memory time moves the bend; a constant-diffusivity ghost shows how wrong Fickian diffusion is "
          "near the source; the plume mode maps t to x/U.",
          ["Press ▶ and watch the local slope readout fall from 2 to 1.",
           "Drag Λ_t up tenfold: the bend moves tenfold later and the final D_T is ten times larger.",
           "Switch to 'plume in a wind': the same curve drawn in space.",
           "Compare with the constant-D ghost at t = 0.1 Λ_t."])
whatif(r"""…the planet rotates and the fluid is stratified on the large scale? The eddies become flat, nearly two-dimensional,
and their energy moves to *larger* scales instead of smaller ones — C07's cascade turned round. The tools of this chapter
(averages, Reynolds stresses, eddy viscosity, Ri and $L_M$, dispersion) all carry over; the next chapter adds the Coriolis
force to them.""")

# =====================================================================================================================
# A.13 §12.13 Concluding Remarks — N217; S01, S02; summary
# =====================================================================================================================
nb.section("12.13", "Concluding Remarks", intro=r"""
**What is this section about?** Where the subject goes from here.""")
note("N217 [C]", r"""
**Turbulence remains an open research field.** What this chapter left out comes next: Ch. 13 adds rotation and large-scale
stratification — eddy-viscosity Ekman layers (C12), surface layers (C11, C15), and geostrophic turbulence, whose energy flows
toward *larger* scales, with a steeper spectrum at small scales than C08's −5/3. Direct and large-eddy simulation are named in
Ch. 10.""")
nb.pointer(r"""**S01** — The chapter's exercises are not reproduced. Results this lesson derived or stated that appear among
them: the moment identities (N22), the real, even spectrum (D05), the periodogram (N38), the Taylor microscale (D04), frozen
turbulence (N40), the mean scalar equation (N57), the mean-flow energy budget (D09), the Reynolds-stress equation (N59), the
isotropic tensor and dissipation (D07, D08), the exponents of the free shear flows (D16), laminar vs turbulent skin friction
(N161), the channel force balance (D17), the pipe friction law (N160), dispersion for a Gaussian correlation (N208).""")
nb.pointer(r"""**S02** — Literature: public sources are cited where used (Sreenivasan 1995; Lee & Moser 2015; Launder & Sharma
1974; McKeon, Zagarola & Smits 2005 for Prandtl's law; Monkewitz, Chauhan & Nagib 2007; Nagib & Chauhan 2008; Pao 1965; Spalding
1961; Coles 1956; Taylor 1921). For more: Pope, *Turbulent Flows*; Tennekes & Lumley, *A First Course in Turbulence*.""")
nb.summary(
    clicked=[
        r"**C01** An average passes through every linear operation and fails on a product: $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$.",
        r"**C02** The area under the autocorrelation is the memory time $\Lambda_t$; a record is worth about $\Delta t/\Lambda_t$ independent samples.",
        r"**C03** The spectrum is the Fourier transform of the autocorrelation: long memory, narrow spectrum; its area is the variance.",
        r"**C04** Averaging the momentum equation leaves one new term, the Reynolds stress $-\rho_0\overline{u_iu_j}$ — a covariance, with $\overline{uv}<0$ in a positive shear, and no equation of its own.",
        r"**C05** With no preferred direction, one function $f(r)$ fixes every two-point correlation, and $\bar\varepsilon=15\nu\overline{u^2}/\lambda_g^2$.",
        r"**C06** Shear production leaves the mean-flow budget and enters the turbulent one with the opposite sign; dissipation ends it.",
        r"**C07** The big eddies set $\bar\varepsilon\sim\Delta U^3/L$; viscosity only sets $\eta=(\nu^3/\bar\varepsilon)^{1/4}$, and $\eta/L\sim\mathrm{Re}_L^{-3/4}$.",
        r"**C08** Between $L$ and $\eta$ only $\bar\varepsilon$ and $k_1$ matter, so $S_{11}\propto\bar\varepsilon^{2/3}k_1^{-5/3}$.",
        r"**C09** A conserved momentum flux plus shape-keeping give $\delta\propto x$ and $U_{CL}\propto x^{-1/2}$ for the plane jet.",
        r"**C10** Near a smooth wall $\tau_0$ and $\nu$ alone set the scales: $U^+=f(y^+)$, and $U^+=y^+$ in the sublayer.",
        r"**C11** Where the answer may depend on neither the viscous nor the outer length, $y\,dU/dy=u_*/\kappa$: a logarithm.",
        r"**C12** An eddy viscosity is a length times a velocity of the flow; $l_T=\kappa y$ returns the log law, and wall damping sets its intercept.",
        r"**C13** The k–ε model builds $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ from two transport equations; its constants are tied to the decay exponent of grid turbulence and to $\kappa$.",
        r"**C14** Rf = buoyant destruction / shear production; turbulence starves near $\mathrm{Rf}\approx\tfrac14$; $\mathrm{Ri}=\mathrm{Pr}_T\,\mathrm{Rf}$.",
        r"**C15** $L_M$ is the height where buoyancy matches shear: $\mathrm{Rf}=z/L_M$.",
        r"**C16** A cloud spreads like $t$ while particles remember their velocity and like $\sqrt t$ afterwards; the eddy diffusivity grows to $\overline{u^2}\Lambda_t$.",
    ],
    feeds_forward=[
        "Ch. 13 (geophysical fluid dynamics): eddy viscosity in Ekman layers; the surface layer and bulk drag; Ri-dependent mixing; geostrophic turbulence spectra; `TS.periodogram` for wave and eddy spectra.",
        "Ch. 14 (aerodynamics): turbulent skin friction on airfoils.",
        "Ch. 15 (compressible flow): friction in ducts.",
    ],
    left_out=[
        "Direct and large-eddy simulation (named in Ch. 10; Pope, *Turbulent Flows*).",
        "Intermittency and structure functions (Frisch, *Turbulence*).",
        "Reynolds-stress closures (Pope, ch. 11).",
        "Two-dimensional and geostrophic turbulence (Ch. 13).",
    ])



# ---------------------------------------------------------------------------------------------------------------------
# Hand-written comments for lines that were left bare (lesson review, round 1): "start of the line ||| what it does".
# ---------------------------------------------------------------------------------------------------------------------
HAND_COMMENTS = r'''
for x in stations: ||| one line of results per station
J_mine, V_mine = [], [] ||| lists for the two fluxes at the three stations
for name, r_ in (("power", pw), ("exponential", ex)): ||| the two trial families, one line each
for label, d_s, U_s, Psi_s in (("printed: ||| the printed exponential family and one that works, treated symbolically
for key, col in (("production", C_RS) ||| the four terms of the model budget, each in its colour
Re_tau = 10.0**logRe ||| the friction Reynolds number δ⁺ of this slider position
d_ = dns[dns["case"] == c] ||| the rows of this DNS case
Up_m = WT.composite_profile_plus(yp_m, 5000.0 ||| the model profile (no wake) in wall units
for Re_tau, col in ((180.0, C_REF) ||| three friction Reynolds numbers, one colour each
yp_ = np.geomspace(1.0, Re_tau, 300) ||| y⁺ from the sublayer to the edge of this layer
out = {} ||| the curves of this slider position, by name
for k_ in kappas: ||| one model profile per value of κ
prof = ch12.mixing_length_wall_profile(yp_f, k_ ||| D21's profile; van Driest damping only when A⁺ > 0
out[f"mixing-length model, κ = {k_}"] ||| store the curve under a name that carries κ
out["reference log line κ = 0.41 ||| the log law with the illustrative pair, for comparison
e_, eps_ = state ||| the two unknowns: turbulent energy ē and dissipation rate ε̄
for i, (name, H) in enumerate(cases): ||| one group of three bars per case
L_ = np.inf if abs(invL) < 1e-9 else 1.0/invL ||| the Monin–Obukhov length of this slider position (infinite when neutral)
nan = np.full_like(z_ax, np.nan) ||| an all-NaN curve: plotly draws nothing for it
for t_ in (1.0, 10.0, 100.0): ||| three times: early, at the memory time, late
i = np.argmin(np.abs(t_l - t_)) ||| index of the sample nearest to this time
mine = [] ||| Taylor's integral at the three times, done by hand
R2_lib = (lib[:, -1]**2).sum(axis=1) ||| squared end-to-end distance of the library's walkers
for k, col in zip(range(3), (C_FLUC, C_RS, C_BUOY)): ||| three individual walkers, one colour each
print("check:", sim["check"]) ||| True: the symbolic similarity equation equals D14's Result
print("the three brackets:", sim["coefficients"]) ||| the three x-dependent coefficients of the similarity equation
print("power family δ ~ x^m ||| what the brackets become for powers of x, and the exponent of the momentum flux
print("exponential family δ ~ e^(ax), U_CL ~ e^(−ax):" ||| the same for the exponential family as printed (middle bracket 0: slip #5)
print("G(ξ = 0.1) for the Gaussian F ||| the stress shape G at ξ = 0.1 from the once-integrated similarity equation
print(f"u* = {u_star:.2f} m/s ||| friction velocity, viscous length and the point in wall units; the sublayer law at y⁺ = 3
print(f"δ+ of a 3 cm layer ||| the friction Reynolds number δ⁺ = δ u*/ν, and the conversion back to metres and m/s
print([WT.layer_name(yp_, yd_) ||| the layer each of four points lies in
ax.semilogx(yp_w, U_sp, ||| the whole inner profile U⁺ = f(y⁺)
ax.semilogx(yp_w[yp_w < 20], WT.viscous_sublayer ||| Eq. (12.82): U⁺ = y⁺, drawn up to y⁺ = 20
ax.semilogx(yp_w[yp_w > 5], WT.log_law ||| Eq. (12.88): the logarithmic law, drawn from y⁺ = 5
for lo, hi, name in ((0.3, 5, "sublayer") ||| shade and name the three layers of the inner region
ax.axvspan(lo, hi, color=C_REF, alpha=0.06 ||| a light band for this layer (darker for the buffer layer)
ax.text(np.sqrt(lo*hi), 23.5, name ||| its name, centred on the logarithmic axis
print(f"sublayer law at y+ = 5 is ||| how far U⁺ = y⁺ is from the curve at y⁺ = 5 and 7 [%]
print(f"U+(100) = ||| Eq. (12.88) at y⁺ = 100, and where it crosses U⁺ = y⁺
print(f"friction law: U_inf+ ||| D19 step 11: the free-stream speed in wall units at δ⁺ = 1000
print("outer groups of U = U(ρ, τ0, δ, y):" ||| Π theorem on the outer list: U/u* and y/δ
print(f"defect variables: ||| a profile in outer variables, and the outer logarithm (12.89) at ξ = 0.1
print("inner solution f =" ||| what sympy finds for the two overlap laws, and its checks
print("Spalding U+ at y+ = 1, 5, 12, 30, 100:" ||| the single-formula inner profile at five heights
print(f"its inverse y+(U+ = 10) ||| Spalding's formula gives y⁺ directly from U⁺; its slope there
print(f"{name}: kappa = {entry['kappa']} ||| one cited pair per line, with the start of its citation
print(f"Eq. (12.92): κ = 0.400 |||the round trip κ → B → κ through the Nagib–Chauhan relation
print(f"Re_x = {bl['Re_x']:.1e}: ||| thicknesses, shape factor, skin friction, u* and δ⁺ of the turbulent layer
print(f"laminar (Blasius) at the same Re_x: ||| what a laminar layer would have at the same place (Ch. 9)
print(f"C_f by {law}: ||| the older skin-friction fits at the same Re_x
print(f"n = {n_d:.4f}, t0 = ||| decay exponent, virtual origin, and ē, ε̄ at t = 0, 1, 10 s
print(f"κ implied by the constants: ||| D23's log-layer relation κ² = √C_μ (C_ε2 − C_ε1) σ_ε
print(f"ν_T = C_μ ē²/ε̄ at t = 0, 10 s: ||| Eq. (12.104) at the start and after 10 s …
f"l_T = ē^1.5/ε̄: ||| … and the length scale ē^(3/2)/ε̄ at the same two times [m]
print(f"L_M = {L_night:.1f} m; ||| the night-time Monin–Obukhov length, Rf at 10 m …
f"U(10 m) = {WT.surface_layer_wind ||| … and the wind at 10 m with and without the stability correction
print(f"regime: {st_n['regime']}; ||| regime, the height where Rf = ¼, and the stability verdict
print(f"noon: L_M = ||| the same for a sunny noon: L_M negative
print("from kinematic fluxes (κ = 0.41):" ||| L_M for an upward, a downward and a zero heat flux
print("D_T at t = 1, 10, 100 s:" ||| Eq. (12.127) by quadrature and in closed form
print("short-time form u²t at t = 1 s:" ||| Eqs. (12.128) and (12.129), each inside its own range
print("the constant under the PRINTED condition ||| slip #12: the constant used at short times (as printed) …
"m²/s — the true value there is" ||| … against the true diffusivity at t = 0.1 s
print(f"from the particles at t = 100 s: ||| ½ d mean(X²)/dt of the particle cloud at late times
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



# =====================================================================================================================
# Final passes: plotting boilerplate commented, labels tidied, equations written next to their numbers, self-checks, save
# =====================================================================================================================
_BOILER = (
    (r"^\s*fig\w*, .*= plt\.subplots\(", "a new figure and its panel(s); figsize in inches"),
    (r"^\s*fig\w* = plt\.figure\(", "a new empty figure; figsize in inches"),
    (r"\.set\((?:title|xlabel|ylabel|xlim|ylim|xscale|yscale|aspect|xticks|yticks)", "axis labels (with units), title and fixed ranges of this panel"),
    (r"\.set_(?:xlabel|ylabel|title|xlim|ylim|xscale|yscale|aspect|xticks|yticks|xticklabels|yticklabels)\(", "label / range / scale of this axis"),
    (r"\.legend\(", "the legend names every curve"),
    (r"\.suptitle\(", "the message of the whole figure"),
    (r"\.tight_layout\(", "tidy the spacing so labels do not overlap"),
    (r"\.colorbar\(", "the colour scale and what it measures"),
    (r"\.axhline\(|\.axvline\(", "a thin reference line"),
    (r"\.grid\(", "light grid lines to read values"),
    (r"^\s*plt\.show\(\)$", "display the figure"),
    (r"^\s*fig\w*\.show\(\)$", "display the interactive figure (the slider works on the web page too)"),
    (r"^\s*display\(", "show the table"),
    (r"\.set_(?:data|ydata|xdata|offsets|array|text)\(", "move the drawn objects to this frame's values"),
    (r"^\s*return\b", "hand the result back to the caller"),
)


def annotate_boilerplate() -> int:
    """A plotting-boilerplate line left without a comment gets a short generic one (hand comments are never touched)."""
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
            if in_str or "#" in ln or not body or body.count("(") != body.count(")") or body.count("[") != body.count("]"):
                new.append(ln)
                continue
            m_pr = re.match(r"""\s*print\(\s*f?r?["']([^"'{]{4,})""", body)
            if m_pr:                                             # a print line keeps no machine-made comment (the list under the cell explains the output)
                new.append(ln)
                continue
            for pat, text in _BOILER:
                if re.search(pat, body):
                    ln = body + "   # " + text
                    changed += 1
                    break
            new.append(ln)
        c.source = "\n".join(new)
    return changed


def uncommented_report() -> list[str]:
    """Code cells the coverage tool would warn about (≥ 4 lines, fewer than half commented) and single lines without a
    comment — printed so they can be fixed by hand."""
    out = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code" or "setup" in c.metadata.get("tags", []):
            continue
        lines = [ln for ln in c.source.splitlines() if ln.strip()]
        if len(lines) >= 4 and sum("#" in ln for ln in lines) < len(lines) / 2:
            out.append(f"cell {i}: {len(lines)} lines, {sum('#' in ln for ln in lines)} commented")
    return out


def relabel() -> int:
    """Label hygiene: the design's working codes for explainers (E1…E10) are replaced by reader-facing words in every
    markdown cell (Part F texts use the codes)."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        s = re.sub(r"\bE(10|[1-9])(?:'s)?\b(?![_^(\[])", lambda m: f"explainer {m.group(1)}" + ("'s" if m.group(0).endswith("'s") else ""), c.source)
        if s != c.source:
            c.source, changed = s, changed + 1
    return changed


_EATEN = re.compile(r"(?<![\\A-Za-z])(rac\{|dfrac\{|tfrac\{|ext\{|heta\b|abla\b|ightarrow|imes\b|ilde\{|oldsymbol|egin\{|overline\{|partial\b|mathrm\{|"
                    r"varepsilon|lambda\b|Lambda\b|kappa\b|cong\b|approx\b|quad\b|sqrt\{|infty\b|displaystyle|bar\{)")


def self_check_eaten() -> list[str]:
    """TeX commands that lost their backslash (and first letter): inside $…$ a bare 'rac{', 'ext{', 'partial' … is the
    remnant of '\\frac', '\\text', … after a form feed or tab was stripped. Also any control character in a cell."""
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "markdown":
            continue
        for mm in re.finditer(r"\$\$?(.+?)\$\$?", c.source, flags=re.S):
            for m in _EATEN.finditer(mm.group(1)):
                seg = mm.group(1)
                bad.append(f"cell {i}: TeX command without its backslash: …{seg[max(0, m.start() - 30):m.end() + 20]!r}")
    return bad


def id_coverage() -> list[str]:
    """Every id of the curation (C, N, R, S and D rows) must be named in the notebook."""
    from nbkit import curation_items
    text = "\n".join(c.source for c in nb.cells if c.cell_type == "markdown")
    ids = set(re.findall(r"\b([CNRSD]\d{2,3})\b", text)) | set(nb.cores) | set(nb.recaps) | set(nb.derivations)
    return sorted((k for k in curation_items("ch12") if k not in ids), key=lambda s: (s[0], int(s[1:])))


nb.cells[0].source = nb.cells[0].source.replace("equations are cited by their numbers so you can follow along",
                                                "equations are shown in full with their book numbers so you can follow along")
def no_stray_reprs() -> int:
    """A figure cell whose last statement is a bare call such as ``ax.legend(...)`` would print that object's repr under
    the figure: end such cells with ``plt.show()`` instead."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "code" or "plt." not in c.source and "ax" not in c.source:
            continue
        if _bare_last_expression(c.source):
            c.source = c.source.rstrip("\n") + "\nplt.show()"
            changed += 1
    return changed


_QUIET_CALLS = {"print", "display", "show_animation", "show_viz", "live"}


def _bare_last_expression(src: str) -> bool:
    """Is the cell's last top-level statement a bare expression whose value Jupyter would echo (an Axes returned by a
    sketch call, a Legend, a list of lines …)? print / display / show_* calls and anything ending in ``.show()`` are not."""
    import ast
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return False
    if not tree.body or not isinstance(tree.body[-1], ast.Expr):
        return False
    val = tree.body[-1].value
    if isinstance(val, ast.Call):
        fn = val.func
        name = fn.id if isinstance(fn, ast.Name) else (fn.attr if isinstance(fn, ast.Attribute) else "")
        if name in _QUIET_CALLS or name == "show":
            return False
    return not isinstance(val, ast.Constant)


def self_check_reprs() -> list[str]:
    """No code cell may end in a bare expression (its repr would be shown), and every figure cell ends in plt.show()."""
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code":
            continue
        if _bare_last_expression(c.source):
            bad.append(f"cell {i}: the last statement is a bare expression — its repr would be printed: {c.source.splitlines()[-1][:80]!r}")
    return bad


def self_check_silent() -> list[str]:
    """Every code cell shows something: a print, a display, a figure, an animation, a widget or an explainer."""
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type == "code" and not re.search(r"print\(|display\(|\.show\(\)|show_animation\(|show_viz\(|^live\(|plt\.", c.source, flags=re.M):
            bad.append(f"cell {i}: no visible output: {c.source.splitlines()[0][:80]!r}")
    return bad


print(f"figure cells closed with plt.show(): {no_stray_reprs()}")
print(f"plotting boilerplate lines commented: {annotate_boilerplate()}")
print(f"labels tidied in {relabel()} cells")
n_changed = finalize_equations()
n_tex = tidy_raw_tex()
print(f"plain-text exponents turned into maths in {n_tex} cells; equations written out in {n_changed} cells")
bad = self_check_reprs() + self_check_silent() + self_check_ctrl() + self_check_eaten() + self_check_numbers() + self_check_prose() + self_check_near()
missing_ids = id_coverage()
print(f"curation ids not named in the notebook: {len(missing_ids)} {missing_ids}")
for _m in uncommented_report():
    print("UNCOMMENTED", _m)
if "--partial" in sys.argv:                                  # development aid: write what exists so far, unchecked
    import nbformat as _nbf
    for _m in bad:
        print("CHECK", _m)
    print("LEDGER", self_check_ledger())
    _nb2 = _nbf.v4.new_notebook(cells=list(nb.cells))
    _nb2.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    _nb2.metadata["fluidpy"] = {"chapter": "ch12", "explainers": nb.explainers, "cores": sorted(nb.cores),
                                "recaps": sorted(nb.recaps), "primers": nb.primers, "derivations": nb.derivations}
    _out = ROOT / "notebooks" / "ch12_turbulence.ipynb"
    for _i, _c in enumerate(_nb2.cells):
        _c["id"] = f"ch12-{_i:03d}"
    _nbf.write(_nb2, _out)
    print("wrote (partial, unchecked)", _out, len(nb.cells), "cells")
    sys.exit(0)
if bad:
    raise SystemExit("markdown cells cite equation numbers without maths (or TeX outside maths):\n  " + "\n  ".join(bad))
if missing_ids:
    raise SystemExit("curation ids with no home in the notebook: " + ", ".join(missing_ids))
miss = self_check_ledger()
if miss:
    raise SystemExit("Part E rows explained by a primer but no primer/reminder in the notebook:\n  " + "\n  ".join(miss))
out = nb.save()
_kinds = {k: sum(k in c.metadata.get("tags", []) for c in nb.cells) for k in ("figure", "animation", "plotly", "live-only", "explainer", "primer", "derivation", "derivation-check", "from-scratch", "recap")}
print(f"wrote {out.relative_to(ROOT)} ({len(nb.cells)} cells: "
      f"{sum(c.cell_type == 'markdown' for c in nb.cells)} markdown, {sum(c.cell_type == 'code' for c in nb.cells)} code) {_kinds}")
