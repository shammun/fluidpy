"""Build the Chapter 16 teaching notebook: ``notebooks/ch16_biofluid_mechanics.ipynb``.

Source of truth: ``analysis/ch16_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract —
``fluidpy.ch16_biofluid_mechanics`` imported as ``ch16`` re-exports ``core.womersley``, ``core.rheology`` and the
Poiseuille, Stokes-layer and acoustics helpers of Ch. 8 and Ch. 15), Part D (runtime budget and FAST sizes), Part E
(prerequisite ledger → primers P374–P394 and one-line reminders of earlier primers), Part F (the 34 derivations D01–D34,
one move per step) and ``analysis/ch16_curation.md`` (IDs, depths, section coverage, the 34 slips and 27 traps), with
these rulings written into the text:

* the twelve physics-changing book slips (S3, S4, S5, S6, S8, S10, S15, S23, S26, S27, S28 and S34, the sign of the
  advection term of (16.174) found in verification) are taught CORRECTED, each with a "slip Sk — the book prints … ;
  the correct form is …" box; the 22 loose slips are named in one clause where they occur;
* traps (true as printed, easy to misread, ``ch16.traps()``) are "Common confusion" callouts; slips and traps kept apart;
* Dean's K⁴ flux coefficient: the book's 0.0120 is shown with the Navier–Stokes-consistent 0.0119312
  (``ch16.DEAN_C4_DERIVED``) beside it (orchestrator note S34); the first axial correction w₁ has the PLUS sign;
* no book number is reproduced (rule 9): worked numbers use the design's own input dicts, written literally in the
  imports cell; ``ch16.illustrative_inputs()`` (the module's own set, also ours) is shown once beside them.

**Derivations are read from Part F at build time** (``part_f()`` below): goal, start, plan, tools, assumptions, every
step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word (step titles exactly as in
Part F, because the explainers copy them); the eight ★★★ derivations (D13, D19, D20, D21, D24, D27, D30, D33) carry the
sympy check cells listed in Part F.  Every bare equation number left in a derivation's text is written out with its
equation (``_expand_refs``), so no equation is cited by number alone.

**Equations.** ``EQ`` below is harvested from the design (the analyst's page-image transcription, re-read on nine pages by
the designer and on eleven more by this builder: 16.5–16.7, 16.30–16.36, 16.153–16.160, 16.172–16.182, 16.187–16.193,
16.200–16.206, 16.221), in the CORRECTED form where the book prints a slip; the printed forms carry a ``p`` suffix.

Run:  .venv/Scripts/python.exe notebooks/build_ch16.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch16.py --dump     (prints the parsed Part F derivations and EQ only)
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

if __name__ == "__main__":                 # unknown flags (e.g. --help by mistake) must not rebuild the notebook
    import argparse
    _ap = argparse.ArgumentParser(description="Build notebooks/ch16_biofluid_mechanics.ipynb from analysis/ch16_design.md.")
    _ap.add_argument("--dump", action="store_true", help="print the parsed Part F derivations and stop (no notebook written)")
    _ap.parse_args()                       # --help prints and exits; an unknown flag is an error

nb = ChapterNotebook("ch16")
DESIGN = (ROOT / "analysis" / "ch16_design.md").read_text(encoding="utf-8")
_CTRL = sorted({hex(ord(_c)) for _c in DESIGN if ord(_c) < 32 and _c != "\n"})
if _CTRL:                                # a form feed, tab or backspace in the design = a LaTeX backslash eaten by a non-raw string
    raise SystemExit(f"analysis/ch16_design.md contains control characters {_CTRL} — repair the eaten backslash before building")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# «16.62» inside any text passed through X() becomes "$c=…$ (16.62)".  EQ is harvested from the design text (every
# "$…$ (16.N)" pair), then corrected / completed by hand below.
# ---------------------------------------------------------------------------------------------------------------------
_PAIR = re.compile(r"\$([^$]+)\$\s*\((16\.\d+)((?:, [^)]*)?)\)")
EQ: dict[str, str] = {}
for _m in _PAIR.finditer(DESIGN.split("## Part B", 1)[0] + DESIGN.split("## Part F", 1)[1]):
    _tex, _num, _lab = _m.group(1).strip(), _m.group(2), _m.group(3)
    if "printed" in _lab:
        continue
    if len(_tex) < 9 or _tex.endswith(("=", "+", "-", r"\neq")) or r"\neq" in _tex or "≠" in _tex:
        continue                                     # fragments and slip-table "≠" entries are never an equation
    if _num not in EQ or "corrected" in _lab:
        EQ[_num] = _tex

EQ_FIX: dict[str, str] = {           # hand-set forms (page-checked) where the harvest caught a fragment or a shorthand
    "16.1": r"SV=EDV-ESV",
    "16.2": r"EF=\dfrac{SV}{EDV}\times100\%",
    "16.3": r"CO=SV\times HR",
    "16.4": r"CI=\dfrac{CO}{BSA}",
    "16.5": r"\sqrt{\tau_y}=(H-0.1)(C_F+0.5)",
    "16.6": r"\sqrt{\dfrac{\tau}{\mu_p}}=k_c\sqrt{\dot\gamma}+\sqrt{\dfrac{\tau_y}{\mu_p}}",
    "16.7": r"\sqrt{\dfrac{\tau}{\mu_p}}=1.53\sqrt{\dot\gamma}+2.0",
    "16.8": r"Q_FH_F=Q_cH_c,\ \ Q_c+Q_p=Q_F,\ \ H_Ta^2=H_c(a-\delta)^2",
    "16.9": r"R=\dfrac{\Delta\bar p}{Q}",
    "16.10": r"\Delta\bar p=\bar p_A-\bar p_V",
    "16.11": r"\bar p_A=\tfrac13p_S+\tfrac23p_D",
    "16.12": r"\alpha=a\sqrt{\omega/\nu}",
    "16.13": r"u(r)=\dfrac{r^2-a^2}{4\mu}\left(\dfrac{dp}{dx}\right)",
    "16.14": r"\dfrac{dp}{dx}=-\dfrac{\Delta p}{L}",
    "16.16": r"u_{max}=\dfrac{\Delta p\,a^2}{4\mu L}",
    "16.17": r"Q=\dfrac{\pi a^4}{8\mu}\dfrac{\Delta p}{L}",
    "16.18": r"V=\dfrac{Q}{\pi a^2}=\dfrac{u_{max}}{2}",
    "16.19": r"\tau=-\mu\left(\dfrac{du}{dr}\right)\Big\rvert_{r=a}=+\dfrac a2\dfrac{\Delta p}{L}",
    "16.20": r"R=\dfrac{8\mu L}{\pi a^4}",
    "16.21": r"\tau_{xr}=-\mu_c\dfrac{du^c}{dr}=+\dfrac r2\dfrac{\Delta p}{L}",
    "16.22": r"\dfrac{du^c}{dr}=0\ \text{at}\ r=0",
    "16.23": r"\tau_{xr}\rvert_c=\tau_{xr}\rvert_p\ \text{at}\ r=a-\delta",
    "16.24": r"\tau_{xr}=-\mu_p\dfrac{du^p}{dr}=+\dfrac r2\dfrac{\Delta p}{L}",
    "16.25": r"u^c=u^p\ \text{at}\ r=a-\delta",
    "16.26": r"u^p=0\ \text{at}\ r=a",
    "16.27": r"u^p=\dfrac{a^2}{4\mu_p}\dfrac{\Delta p}{L}\left[1-\left(\dfrac ra\right)^2\right]",
    "16.30": r"Q_c=\dfrac{\pi(a-\delta)^2\Delta p}{4\mu_pL}\left[a^2-\left(1-\dfrac{\mu_p}{2\mu_c}\right)(a-\delta)^2\right]",
    "16.30p": r"Q_c=\dfrac{\pi a^2\Delta p}{4\mu_pL}\left[a^2-\left(1-\dfrac{\mu_p}{2\mu_c}\right)\dfrac{(a-\delta)^4}{a^2}\right]",
    "16.31": r"Q=Q_p+Q_c=\dfrac{\pi a^4\Delta p}{8\mu_pL}\left[1-\left(1-\dfrac\delta a\right)^4\left(1-\dfrac{\mu_p}{\mu_c}\right)\right]",
    "16.32": r"Q=\dfrac{\pi a^4\Delta p}{8\mu_{app}L}",
    "16.33": r"\mu_{app}=\mu_p\left[1-\left(1-\dfrac\delta a\right)^4\left(1-\dfrac{\mu_p}{\mu_c}\right)\right]^{-1}",
    "16.34": r"\mu_{app}\approx\mu_c\left[1+4\dfrac\delta a\left(\dfrac{\mu_c}{\mu_p}-1\right)\right]^{-1}",
    "16.35": r"\mu_c=\mu_p\dfrac{1}{1-\alpha_cH_c}",
    "16.36": r"\alpha_c=0.070\exp\left[2.49H_c+\dfrac{1107}{T}\exp(-1.69H_c)\right]",
    "16.37": r"\dfrac\ell d=0.06\,Re\ \text{(laminar)},\ \ \dfrac\ell d=0.693\,Re^{1/4}\ \text{(turbulent)}",
    "16.38": r"2\sigma_{\theta\theta}h=\displaystyle\int_0^\pi(p-p_e)\,a\sin\theta\,d\theta",
    "16.39": r"\sigma_{\theta\theta}=\dfrac{(p-p_e)a}{h}",
    "16.40": r"e_{\theta\theta}=\dfrac{\sigma_{\theta\theta}}{E}=\dfrac{a-a_0}{a_0}",
    "16.41": r"a=a_0\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-1}",
    "16.42": r"Q=-\dfrac{\pi}{8\mu}\left(\dfrac{dp}{dx}\right)a(x)^4",
    "16.43": r"\dfrac{dp}{dx}=-\dfrac{8\mu Q}{\pi a(x)^4}",
    "16.44": r"\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-4}dp=-\dfrac{8\mu}{\pi a_0^4}Q\,dx",
    "16.44p": r"\left[1-\dfrac{a_0}{Eh}(p-p_e)\right]^{-4}dp=-\dfrac{8\mu}{\pi a_0}Q\,dx",
    "16.45": (r"\dfrac{Eh}{3a_0}\left\{\left[1-\dfrac{a_0}{Eh}(p_2-p_e)\right]^{-3}-\left[1-\dfrac{a_0}{Eh}(p_1-p_e)\right]^{-3}"
              r"\right\}=-\dfrac{8\mu}{\pi a_0^4}LQ"),
    "16.45p": (r"\dfrac{Eh}{3a_0}\left\{\left[1-\dfrac{a_0}{Eh}(p_2-p_e)\right]^{-3}-\left[1-\dfrac{a_0}{Eh}(p_1-p_e)\right]^{-3}"
               r"\right\}=-\dfrac{8\mu}{\pi a_0}LQ"),
    "16.46": r"Q=\dfrac{\pi a_0^3Eh}{24\mu L}\left\{\left[1-\dfrac{a_0}{Eh}(p_1-p_e)\right]^{-3}-\left[1-\dfrac{a_0}{Eh}(p_2-p_e)\right]^{-3}\right\}",
    "16.47": r"Q(t)=\dfrac{dv}{dt}+Q_{out}(t)",
    "16.48": r"\dfrac{dv}{dt}=\left(\dfrac{dv}{dp}\right)\left(\dfrac{dp}{dt}\right)=\mathrm K\dfrac{dp}{dt}",
    "16.49": r"Q(t)=\mathrm K\dfrac{dp}{dt}+\dfrac{p(t)}{R}",
    "16.50": r"Q=\dfrac{dy}{dx}+Py",
    "16.51": r"ye^{\int P\,dx}=A+\displaystyle\int Qe^{\int P\,dx}dx",
    "16.52": r"p(t)=\dfrac1{\mathrm K}e^{-t/R\mathrm K}\displaystyle\int_0^tQ(\tau)e^{\tau/R\mathrm K}d\tau+p_0e^{-t/R\mathrm K}",
    "16.53": r"\dfrac{\partial A}{\partial t}+\dfrac{\partial(Au)}{\partial x}=0",
    "16.54": r"\rho A\left(\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}\right)=-A\dfrac{\partial p}{\partial x}",
    "16.54p": r"\rho A\left(\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}\right)=-\dfrac{\partial((p-p_e)A)}{\partial x}",
    "16.55": r"p-p_e=\dfrac{Eh}{a_0}\left[1-\left(\dfrac{A_0}{A}\right)^{1/2}\right]",
    "16.56": r"\dfrac{\partial A}{\partial t}+A_0\dfrac{\partial u}{\partial x}=0",
    "16.57": r"\rho\dfrac{\partial u}{\partial t}=-\dfrac{\partial p}{\partial x}",
    "16.58": r"p-p_e=\dfrac{Eh}{2a_0A_0}(A-A_0)",
    "16.61": r"\dfrac{\partial^2p}{\partial t^2}=c^2\dfrac{\partial^2p}{\partial x^2},\quad c^2=\dfrac{A_0}{\rho}\dfrac{\partial p}{\partial A}",
    "16.62": r"c=\sqrt{\dfrac{Eh}{2\rho a_0}}",
    "16.63": r"c=\sqrt{\dfrac{Eh}{2\rho(a_0+h/2)}}",
    "16.65": r"p=f_1\left(t-\dfrac xc\right)+f_2\left(t+\dfrac xc\right)",
    "16.68": r"\hat p_1=\rho c\hat u_1,\ \ \hat p_2=-\rho c\hat u_2",
    "16.71": r"\rho\dfrac{\partial u}{\partial t}=-\dfrac{\partial p}{\partial x}+\mu\left[\dfrac{\partial^2u}{\partial r^2}+\dfrac1r\dfrac{\partial u}{\partial r}\right]",
    "16.72": r"\mu\left[\dfrac{\partial^2u}{\partial r^2}+\dfrac1r\dfrac{\partial u}{\partial r}\right]-\rho\dfrac{\partial u}{\partial t}=\dfrac{\partial p}{\partial x}",
    "16.73": r"\dfrac{\partial p}{\partial x}=-Ge^{i\omega t}",
    "16.75": r"u(r,t)=U(r)e^{i\omega t}",
    "16.76": r"\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=-\dfrac G\mu",
    "16.76p": r"\dfrac{d^2U}{dr^2}+\dfrac1r\dfrac{dU}{dr}-\dfrac{i\omega\rho}{\mu}U=+\dfrac G\mu",
    "16.81": r"F_1(\alpha)=\dfrac{J_0(i^{3/2}\alpha\,r/a)}{J_0(i^{3/2}\alpha)}",
    "16.82": r"u(r,t)=\dfrac{Ga^2}{i\mu\alpha^2}\left(1-\dfrac{J_0(i^{3/2}\alpha r/a)}{J_0(i^{3/2}\alpha)}\right)e^{i\omega t}",
    "16.83": r"\dfrac{\partial p}{\partial x}=-\sum_{n=0}^NG_ne^{in\omega t}",
    "16.85": r"u_0(r)=\dfrac{G_0a^2}{4\mu}\left(1-\dfrac{r^2}{a^2}\right)",
    "16.86": r"u_n(r)=\dfrac{G_na^2}{i\mu\alpha_n^2}(1-F_1(\alpha_n))",
    "16.87": r"u_n(r)\approx\dfrac{G_na^2}{4\mu}\left(1-\dfrac{r^2}{a^2}\right)",
    "16.88": r"u\approx\dfrac{Ga^2}{i\mu\alpha^2}\left\{1-\exp\left[-\sqrt{\dfrac{\omega}{2\nu}}(1+i)(a-r)\right]\right\}e^{i\omega t}",
    "16.90": r"Q(t)=\dfrac{\pi a^4}{i\mu\alpha^2}G(1-F_2(\alpha))e^{i\omega t}",
    "16.91": r"F_2(\alpha)=\dfrac{2J_1(i^{3/2}\alpha)}{i^{3/2}\alpha J_0(i^{3/2}\alpha)}",
    "16.92": r"\tau_w=-\mu\dfrac{\partial u}{\partial r}\Big\rvert_a=\dfrac{G_0a}2+\dfrac a2\sum_1^NG_nF_2(\alpha_n)e^{in\omega t}",
    "16.93": r"F_2(\alpha)\approx1-\dfrac{i\alpha^2}8",
    "16.94": r"Q=\dfrac{\pi Ga^4}{8\mu}e^{i\omega t}",
    "16.94p": r"Q=\dfrac{\pi Ga^4}{8\mu}e^{in\omega t}",
    "16.95": r"|Q_0|=\dfrac{\pi Ga^4}{8\mu}",
    "16.96": r"F_2(\alpha)\approx\dfrac{2}{i^{1/2}\alpha}\left(1+\dfrac{i^{3/2}}{2\alpha}\right)",
    "16.96p": r"F_2(\alpha)\approx\dfrac{2}{i^{1/2}\alpha}\left(1+\dfrac1{2\alpha}\right)",
    "16.97": r"1-F_2(\alpha)=Z(\alpha)=X(\alpha)+iY(\alpha)",
    "16.98": r"Q=\dfrac{\pi Ga^4}{\mu\alpha^2}\left\{[Y\cos\omega t+X\sin\omega t]-i[X\cos\omega t-Y\sin\omega t]\right\}",
    "16.99": r"|Q|=\dfrac{\pi Ga^4}{\mu\alpha^2}\sqrt{X^2+Y^2}",
    "16.100": r"\tan\phi=\dfrac XY",
    "16.101": r"\dfrac{\partial u}{\partial x}+\dfrac1r\dfrac{\partial(rv)}{\partial r}=0",
    "16.102": r"\rho\dfrac{\partial u}{\partial t}=-\dfrac{\partial p}{\partial x}+\mu\left(\dfrac{\partial^2u}{\partial r^2}+\dfrac1r\dfrac{\partial u}{\partial r}+\dfrac{\partial^2u}{\partial x^2}\right)",
    "16.103": r"\rho\dfrac{\partial v}{\partial t}=-\dfrac{\partial p}{\partial r}+\mu\left(\dfrac{\partial^2v}{\partial r^2}+\dfrac1r\dfrac{\partial v}{\partial r}-\dfrac{v}{r^2}+\dfrac{\partial^2v}{\partial x^2}\right)",
    "16.104": r"T_\theta=\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac\eta{a_0}+\hat\nu\dfrac{\partial\xi}{\partial x}\right)",
    "16.106": r"\rho_wh\dfrac{\partial^2\eta}{\partial t^2}=\sigma_{rr}\rvert_a-\dfrac{T_\theta}{a_0}",
    "16.108": r"\rho_wh\dfrac{\partial^2\eta}{\partial t^2}=\sigma_{rr}\rvert_a-\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac{\eta}{a_0^2}+\dfrac{\hat\nu}{a_0}\dfrac{\partial\xi}{\partial x}\right)",
    "16.110": r"\sigma_{rr}=+p-2\mu\dfrac{\partial v}{\partial r}",
    "16.111": r"\sigma_{rx}=\mu\left(\dfrac{\partial u}{\partial r}+\dfrac{\partial v}{\partial x}\right)",
    "16.112": (r"\rho_wh\dfrac{\partial^2\eta}{\partial t^2}=+p\rvert_{r=a}-2\mu\dfrac{\partial v}{\partial r}\Big\rvert_{r=a}"
               r"-\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac\eta{a_0^2}+\dfrac{\hat\nu}{a_0}\dfrac{\partial\xi}{\partial x}\right)"),
    "16.113": (r"\rho_wh\dfrac{\partial^2\xi}{\partial t^2}=-\mu\left(\dfrac{\partial u}{\partial r}+\dfrac{\partial v}{\partial x}"
               r"\right)\Big\rvert_{r=a}+\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac{\partial^2\xi}{\partial x^2}+\dfrac{\hat\nu}{a_0}"
               r"\dfrac{\partial\eta}{\partial x}\right)"),
    "16.114": r"u\rvert_{a_0}=\dfrac{\partial\xi}{\partial t}",
    "16.115": r"v\rvert_{a_0}=\dfrac{\partial\eta}{\partial t}",
    "16.116": r"(u,v,p,\xi,\eta)=(\hat u(r),\hat v(r),\hat p,\hat\xi,\hat\eta)\,e^{i(kx-\omega t)}",
    "16.117": r"\left|\dfrac{v}{u}\right|=O(|ak|)",
    "16.119": r"\rho\dfrac{\partial u}{\partial t}=-\dfrac{\partial p}{\partial x}+\mu\left(\dfrac{\partial^2u}{\partial r^2}+\dfrac1r\dfrac{\partial u}{\partial r}\right)",
    "16.120": r"\dfrac{\partial p}{\partial r}=0",
    "16.122": r"\dfrac{d(r\hat v)}{dr}=-ikr\hat u",
    "16.123": r"\dfrac{d^2\hat u}{dr^2}+\dfrac1r\dfrac{d\hat u}{dr}+\dfrac{i\omega\rho}{\mu}\hat u=\dfrac{ik\hat p}{\mu}",
    "16.124": r"\hat u(a_0)=-i\omega\hat\xi",
    "16.125": r"\hat v(a_0)=-i\omega\hat\eta",
    "16.126": r"-\rho_wh\omega^2\hat\eta=\hat p-2\mu\dfrac{d\hat v}{dr}\Big\rvert_{a_0}-\dfrac{Eh}{1-\hat\nu^2}\left(\dfrac{\hat\eta}{a_0^2}+\dfrac{i\hat\nu k}{a_0}\hat\xi\right)",
    "16.127": r"-\rho_wh\omega^2\hat\xi=-\mu\left(\dfrac{d\hat u}{dr}+ik\hat v\right)\Big\rvert_{a_0}+\dfrac{Eh}{1-\hat\nu^2}\left(-k^2\hat\xi+\dfrac{i\hat\nu k}{a_0}\hat\eta\right)",
    "16.128": r"\hat p=\dfrac{Eh}{a_0^2}\hat\eta-\dfrac{i\hat\nu}{a_0k}\mu\dfrac{d\hat u}{dr}\Big\rvert_{a_0}",
    "16.129": r"\hat\xi=\dfrac{i\hat\nu}{ka_0}\hat\eta-\dfrac{1-\hat\nu^2}{Ehk^2}\mu\dfrac{d\hat u}{dr}\Big\rvert_{a_0}",
    "16.132": r"\hat v=-\dfrac{ik}{r}\displaystyle\int_0^rs\,\hat u(s)\,ds",
    "16.136": r"X^2\left[2\hat\beta\dfrac{J_0(\hat\beta)}{J_1(\hat\beta)}-4\right]+X\left[4\hat\nu-1-2\hat\beta\dfrac{J_0(\hat\beta)}{J_1(\hat\beta)}\right]+(1-\hat\nu^2)=0",
    "16.137": r"E^*=E-i\omega E',\ \ \hat\nu^*=\hat\nu-i\omega\hat\nu'",
    "16.138": r"\hat E_1=\dfrac{8\mu L}{\pi a^4}Q^2+K\pi a^2L",
    "16.140": r"a=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q^{1/3}",
    "16.140p": r"a=\left[\dfrac{16\mu}{\pi^2}K\right]^{1/6}Q^{1/3}",
    "16.141": r"\hat E_{1,\min}=\dfrac{3\pi}2KLa^2",
    "16.142": r"\hat E_{2,\min}=\dfrac{3\pi}2K\left(a_0^2L_0+a_1^2L_1+a_2^2L_2\right)",
    "16.143": r"a_i=\left[\dfrac{16\mu}{\pi^2K}\right]^{1/6}Q_i^{1/3}",
    "16.144": r"a_0^3=a_1^3+a_2^3",
    "16.145": r"\delta\hat E_{2,\min}=\dfrac{3\pi}2K\left(\delta L_0\,a_0^2+\delta L_1\,a_1^2+\delta L_2\,a_2^2\right)",
    "16.146": r"\cos\theta=\dfrac{a_0^4+a_1^4-a_2^4}{2a_0^2a_1^2},\ \cos\phi=\dfrac{a_0^4-a_1^4+a_2^4}{2a_0^2a_2^2},\ \cos(\theta+\phi)=\dfrac{a_0^4-a_1^4-a_2^4}{2a_1^2a_2^2}",
    "16.147": r"\cos\theta=\dfrac{a_0^4+a_1^4-(a_0^3-a_1^3)^{4/3}}{2a_0^2a_1^2}",
    "16.148": r"p_I+p_R=p_{T1}=p_{T2}",
    "16.149": r"Q_I-Q_R=Q_{T1}+Q_{T2}",
    "16.151": r"Q=Au=\pm Yp=\pm\dfrac pZ",
    "16.152": r"Y_1(p_I-p_R)=\sum_jY_jp_{Tj}",
    "16.153": r"\mathrm R=\dfrac{p_R}{p_I}=\dfrac{Y_1-\sum Y_j}{Y_1+\sum Y_j},\ \ \mathrm T=\dfrac{p_{Tj}}{p_I}=\dfrac{2Y_1}{Y_1+\sum Y_j}",
    "16.154": r"\dfrac{p_R}{p_I}=\dfrac{Z_1^{-1}-\sum Z_j^{-1}}{Z_1^{-1}+\sum Z_j^{-1}}",
    "16.159": r"\dfrac p{P_I}=f\left(t-\dfrac x{c_1}\right)+\dfrac{P_R}{P_I}f\left(t+\dfrac x{c_1}\right)",
    "16.160": r"Q=Y_1P_I\left[f\left(t-\dfrac x{c_1}\right)-\dfrac{P_R}{P_I}f\left(t+\dfrac x{c_1}\right)\right]",
    "16.161": r"\text{incident power}=\dfrac{p_I^2}{Z_1}",
    "16.162": r"\text{reflected power}=\mathrm R^2\dfrac{p_I^2}{Z_1}",
    "16.163": r"\dfrac{p_{T2}^2/Z_2+p_{T3}^2/Z_3}{p_I^2/Z_1}=\dfrac{Z_2^{-1}+Z_3^{-1}}{Z_1^{-1}}\mathrm T^2",
    "16.164": r"\delta=\dfrac aR",
    "16.165": r"p=p_0(s)+\delta p_1(r,\alpha,s)+\dots",
    "16.166": r"\dfrac{\partial u}{\partial r}+\dfrac ur+\dfrac1r\dfrac{\partial v}{\partial\alpha}=0",
    "16.167": (r"u\dfrac{\partial u}{\partial r}+\dfrac vr\dfrac{\partial u}{\partial\alpha}-\dfrac{v^2}r-w^2\cos\alpha="
               r"-\dfrac{\partial p_1}{\partial r}-\dfrac2\kappa\dfrac1r\dfrac{\partial}{\partial\alpha}\left(\dfrac{\partial v}"
               r"{\partial r}+\dfrac vr-\dfrac1r\dfrac{\partial u}{\partial\alpha}\right)"),
    "16.168": (r"u\dfrac{\partial v}{\partial r}+\dfrac vr\dfrac{\partial v}{\partial\alpha}+\dfrac{uv}r+w^2\sin\alpha="
               r"-\dfrac1r\dfrac{\partial p_1}{\partial\alpha}+\dfrac2\kappa\dfrac{\partial}{\partial r}\left(\dfrac{\partial v}"
               r"{\partial r}+\dfrac vr-\dfrac1r\dfrac{\partial u}{\partial\alpha}\right)"),
    "16.169": (r"u\dfrac{\partial w}{\partial r}+\dfrac vr\dfrac{\partial w}{\partial\alpha}=-\dfrac{\partial p_0}{\partial s}"
               r"+\dfrac2\kappa\left(\dfrac{\partial^2w}{\partial r^2}+\dfrac1r\dfrac{\partial w}{\partial r}+\dfrac1{r^2}"
               r"\dfrac{\partial^2w}{\partial\alpha^2}\right)"),
    "16.170": r"u=v=w=0\ \text{at}\ r=1,\quad\text{no singularity at}\ r=0",
    "16.171": r"\kappa=\sqrt\delta\,\dfrac{2a\bar W_0}\nu",
    "16.172": r"u=\dfrac1r\dfrac{\partial\psi}{\partial\alpha},\ \ v=-\dfrac{\partial\psi}{\partial r}",
    "16.173": r"\nabla_1^2w-\dfrac\kappa2\dfrac{\partial p_0}{\partial s}=\dfrac\kappa{2r}\left(\dfrac{\partial\psi}{\partial\alpha}\dfrac{\partial w}{\partial r}-\dfrac{\partial\psi}{\partial r}\dfrac{\partial w}{\partial\alpha}\right)",
    "16.174": r"\dfrac2\kappa\nabla_1^4\psi+\dfrac1r\left(\dfrac{\partial\psi}{\partial r}\dfrac{\partial}{\partial\alpha}-\dfrac{\partial\psi}{\partial\alpha}\dfrac{\partial}{\partial r}\right)\nabla_1^2\psi=-2w\left(\sin\alpha\dfrac{\partial w}{\partial r}+\dfrac{\cos\alpha}r\dfrac{\partial w}{\partial\alpha}\right)",
    "16.174p": r"\dfrac2\kappa\nabla_1^4\psi-\dfrac1r\left(\dfrac{\partial\psi}{\partial r}\dfrac{\partial}{\partial\alpha}-\dfrac{\partial\psi}{\partial\alpha}\dfrac{\partial}{\partial r}\right)\nabla_1^2\psi=-2w\left(\sin\alpha\dfrac{\partial w}{\partial r}+\dfrac{\cos\alpha}r\dfrac{\partial w}{\partial\alpha}\right)",
    "16.175": r"\nabla_1^2=\dfrac{\partial^2}{\partial r^2}+\dfrac1r\dfrac{\partial}{\partial r}+\dfrac1{r^2}\dfrac{\partial^2}{\partial\alpha^2}",
    "16.176": r"\psi=\dfrac{\partial\psi}{\partial r}=w=0\ \text{at}\ r=1",
    "16.177": r"w=\sum\kappa^{2n}w_n,\ \ \psi=\kappa\sum\kappa^{2n}\psi_n",
    "16.178": r"\dfrac{Q}{\pi a^2\bar W}=1-0.0306\left(\dfrac K{576}\right)^2+0.0120\left(\dfrac K{576}\right)^4+O(K^6)",
    "16.179": r"D=\sqrt{2\delta}\,\dfrac{\hat Ga^2}{\mu}\dfrac a\nu",
    "16.180": r"\hat G=-\dfrac{8\mu\bar W}{a^2}",
    "16.181": r"\dfrac Q{\pi a^2\bar W}=1-0.0306\left(\dfrac D{96}\right)^4+0.0120\left(\dfrac D{96}\right)^8+O(D^{12})",
    "16.182": r"\lambda=\dfrac{\lambda_c}{\lambda_s}=1+0.0306\left(\dfrac K{576}\right)^2-0.0110\left(\dfrac K{576}\right)^4+\dots",
    "16.183": r"p-p_e=\hat P(A)",
    "16.185": r"c^2=\dfrac A\rho\dfrac{d(p-p_e)}{dA}",
    "16.186": r"\mathrm P=\dfrac{p-p_e}{K_p},\ \ \alpha=\dfrac A{A_0}",
    "16.187": r"K_p\propto EI",
    "16.188": r"-\mathrm P\approx\alpha^{-n}-1,\ \ n=\tfrac32",
    "16.189": r"c^2=\dfrac{nK_p\alpha^{-n}}{\rho}",
    "16.190": r"\dfrac{dA}A+\dfrac{du}u=0",
    "16.191": r"-A\,dp-\tau_ws\,dx-\rho gA\,dz=\rho Au\,du",
    "16.192": r"\dfrac{\tau_ws\,dx}{A}=\dfrac12\rho u^2\dfrac{4f_T\,dx}{d_e}",
    "16.193": r"\dfrac{\tau_ws\,dx}{A}=\dfrac{\mu u}{d_0}\dfrac1\alpha\dfrac{4f_L'\,dx}{d_0}",
    "16.194": r"d(p+\rho gz)+\dfrac{\tau_ws\,dx}A-\rho u^2\dfrac{dA}A=0",
    "16.195": r"S=\dfrac uc,\ \ \dfrac{dS^2}{S^2}=2\dfrac{du}u-2\dfrac{dc}c",
    "16.196": r"\dfrac{d\alpha}\alpha=\dfrac{dA}A-\dfrac{dA_0}{A_0}",
    "16.197": r"dp=dp_e+K_p\,d\mathrm P+\mathrm P\,dK_p",
    "16.198": r"c^2(\alpha,x)=\dfrac\alpha\rho K_p\dfrac{\partial\mathrm P}{\partial\alpha}\Big\rvert_x",
    "16.199": r"d\mathrm P=\dfrac{\partial\mathrm P}{\partial\alpha}d\alpha+\dfrac{\partial\mathrm P}{\partial x}dx",
    "16.200": r"dp=dp_e+\rho c^2\dfrac{d\alpha}\alpha+K_p\dfrac{\partial\mathrm P}{\partial x}dx+\mathrm P\,dK_p",
    "16.201": r"2\dfrac{dc}c=\left(1+\dfrac{\alpha\,\partial^2\mathrm P/\partial\alpha^2}{\partial\mathrm P/\partial\alpha}\right)\dfrac{d\alpha}\alpha+\dfrac{dK_p}{K_p}+\dfrac{\alpha K_p}{\rho c^2}\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial\alpha}\right)dx",
    "16.201p": r"\dots+\dfrac{\alpha K_p}{\rho c^2}\dfrac{\partial}{\partial x}\left(\dfrac{\partial\mathrm P}{\partial x}\right)dx",
    "16.202": r"\dfrac{dS^2}{S^2}=-2\dfrac{d\alpha}\alpha-2\dfrac{dA_0}{A_0}-2\dfrac{dc}c",
    "16.203": r"p-p_e(x)=K_p(x)\mathrm P(\alpha)",
    "16.204": r"(1-S^2)\dfrac1\alpha\dfrac{d\alpha}{dx}=\dfrac{S^2}{A_0}\dfrac{dA_0}{dx}-\dfrac1{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ+\mathrm P\dfrac{dK_p}{dx}\right]",
    "16.205": r"(1-S^2)\dfrac1{S^2}\dfrac{dS^2}{dx}=\dfrac1{A_0}\dfrac{dA_0}{dx}[-2+(2-\mathrm M)S^2]+\dfrac{\mathrm M}{\rho c^2}\left[\dfrac{dp_e}{dx}+\rho g\dfrac{dz}{dx}+\mathrm RQ\right]+\dfrac1{\rho c^2}\dfrac{dK_p}{dx}\left[\mathrm M\mathrm P-(1-S^2)\alpha\dfrac{d\mathrm P}{d\alpha}\right]",
    "16.206": r"\mathrm M=3+\dfrac{\alpha\,\partial^2\mathrm P/\partial\alpha^2}{\partial\mathrm P/\partial\alpha}",
    "16.207": r"\tau^{1/2}=\tau_y^{1/2}+K_c\dot\gamma^{1/2}\ (\tau\ge\tau_y),\ \ \dot\gamma=0\ (\tau<\tau_y)",
    "16.208": r"\tau_{rx}\,2\pi r\,\Delta x=-\pi r^2\left(p\rvert_{x+\Delta x}-p\rvert_x\right)",
    "16.209": r"\tau_{rx}(r)=-\dfrac r2\dfrac{dp}{dx}=\dfrac{(p_1-p_2)r}{2L}",
    "16.210": r"a_c=\dfrac{2L\tau_y}{p_1-p_2}=a\left(\dfrac{\tau_y}{\tau_w}\right)",
    "16.211": r"\dot\gamma=-\dfrac{du}{dr}=\dfrac{\tau+\tau_y-2\sqrt{\tau\tau_y}}{K_c^2}",
    "16.215": r"u=\dfrac1{4K_c^2}\dfrac{\Delta p}L\left[(a^2-r^2)-\dfrac83\sqrt{a_c}\left(\sqrt{a^3}-\sqrt{r^3}\right)+2a_c(a-r)\right]",
    "16.216": r"u=\dfrac{a\tau_w}{2K_c^2}\left\{\left[1-\left(\dfrac ra\right)^2\right]-\dfrac83\sqrt{\dfrac{\tau_y}{\tau_w}}\left[1-\left(\dfrac ra\right)^{3/2}\right]+2\dfrac{\tau_y}{\tau_w}\left(1-\dfrac ra\right)\right\}",
    "16.217": r"\dfrac ra=\dfrac{a_c}a=\dfrac{\tau_y}{\tau_w}",
    "16.218": r"u_c=\dfrac1{4K_c^2}\dfrac{\Delta p}L\left(\sqrt a-\sqrt{a_c}\right)^3\left(\sqrt a+\tfrac13\sqrt{a_c}\right)",
    "16.219": r"Q=\pi a_c^2u_c+\displaystyle\int_{a_c}^a2\pi ru\,dr",
    "16.220": r"Q=\dfrac\pi8\dfrac1{K_c^2}\dfrac{\Delta p}La^4\left[1-\dfrac{16}7\left(\dfrac{a_c}a\right)^{1/2}+\dfrac43\left(\dfrac{a_c}a\right)-\dfrac1{21}\left(\dfrac{a_c}a\right)^4\right]",
    "16.221": r"\psi=p-RTc+\rho gz",
    "16.222": r"\Delta p=\dfrac{8\mu Q}{\pi a^4}\left[L+\dfrac\ell N\left(\dfrac ar\right)^4\right]+2\,\Delta p'",
    # earlier chapters (written-out forms, as in their notebooks)
    "1.3": r"\tau=\mu\dfrac{du}{dy}",
    "1.8": r"\dfrac{dp}{dz}=-\rho g",
    "8.6": r"u_z=\dfrac{R^2-a^2}{4\mu}\dfrac{dp}{dz}",
    "8.8": r"\tau_0=\dfrac a2\dfrac{dp}{dz}",
    "8.38": r"u(y,t)=U\exp\{-y\sqrt{\omega/2\nu}\}\cos(\omega t-y\sqrt{\omega/2\nu})",
    "15.15": r"p'(x,t)=f(x-ct)+g(x+ct)",
    "15.16": r"u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)",
}
EQ.update(EQ_FIX)
LABEL = {
    "16.19": "16.19, last member's sign corrected", "16.21": "16.21, sign corrected", "16.24": "16.24, sign corrected",
    "16.30": "16.30, corrected", "16.30p": "16.30, as printed", "16.44": "16.44, corrected", "16.44p": "16.44, as printed",
    "16.54": "16.54, corrected", "16.54p": "16.54, as printed", "16.76": "16.76, right side corrected",
    "16.76p": "16.76, as printed", "16.92": "16.92, corrected", "16.96": "16.96, corrected", "16.96p": "16.96, as printed",
    "16.140": "16.140, corrected", "16.140p": "16.140, as printed", "16.143": "16.143, with the corrected bracket",
    "16.174": "16.174, advection sign corrected — S34", "16.174p": "16.174, as printed",
    "16.201": "16.201, corrected", "16.201p": "16.201, last term as printed",
    "16.209": "16.209, first member's sign corrected",
    "16.45": "16.45, corrected", "16.45p": "16.45, as printed", "16.94": "16.94, single-mode exponent",
    "16.94p": "16.94, as printed", "16.180": "16.180, as printed",
}
for _k, _v in EQ.items():
    assert "$" not in _v, (_k, _v)
    assert "≠" not in _v and r"\neq" not in _v and not _v.rstrip().endswith("="), (_k, _v)   # never a slip-table fragment


def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ ({LABEL.get(num, num)})"


def X(s: str) -> str:
    """Replace every «key» in a text by the written-out equation and its number (E)."""
    return re.sub(r"«([^»]+)»", lambda m: E(m.group(1)), textwrap.dedent(s).strip())


def _join(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(x.strip() for x in lines)).strip()


def _abs_plain(s: str) -> str:
    r"""Part F writes absolute values as \lvert … \rvert in places; outside $…$ they become plain |…|."""
    parts = re.split(r"(\$[^$]+\$)", s)
    return "".join(q if q.startswith("$") else re.sub(r"\\rvert", "|", re.sub(r"\\lvert\s?", "|", q)) for q in parts)


def _tidy(s: str) -> str:
    """Plain-text exponents of Part F outside maths and code become readable markdown."""
    parts = re.split(r"(\$[^$]+\$|`[^`\n]*`)", s)
    out = []
    for q in parts:
        if q.startswith("$") or q.startswith("`"):
            out.append(q)
            continue
        q = re.sub(r"\)\^\{([^{}]*)\}", r")^(\1)", q)
        q = re.sub(r"([A-Za-zα-ωζ0-9]+)\^\{([^{}]*)\}", r"$\1^{\2}$", q)
        out.append(q)
    return "".join(out)


def _tex_escape(s: str) -> str:
    return (s.replace("\\", r"\backslash ").replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")
            .replace("_", r"\_").replace("{", r"\{").replace("}", r"\}"))


def display_tex(s: str, width: int = 78) -> str:
    """A Part F line (maths in `$…$`, possibly with words around it) as the body of one display-math block."""
    s = s.strip()
    parts = [p for p in re.split(r"(\$[^$]+\$)", s) if p.strip()]
    if len(parts) == 1 and parts[0].startswith("$"):
        return parts[0][1:-1]
    tokens: list[tuple[str, str]] = []
    for p in parts:
        if p.startswith("$"):
            tokens.append(("m", p[1:-1]))
        else:
            tokens += [("w", w.replace("**", "").replace("*", "").replace("`", "")) for w in p.split()]
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


# Bare equation numbers inside derivation prose: "(16.43)" → "$…$ (16.43)" — the first mention in each text field, unless
# the number is already written out next to its equation there or in the step's own line.  Ranges "(16.21)–(16.26)" name
# groups that the derivation's start writes out and stay as they are.
_BARE = re.compile(r"(?<![\$\w])(?<!\$ )\((16\.\d+)((?:, [^)]*)?)\)(?!\s*[–-]\s*\(16\.)")


def _expand_refs(text: str, skip: set[str] | None = None) -> str:
    skip = set(skip or ())
    seen: set[str] = set()
    parts = re.split(r"(\$[^$]+\$|`[^`\n]*`)", text)
    shown = {m.group(1) for m in re.finditer(r"\$\s*\((16\.\d+)", text)}
    out = []
    prev_end_dollar = False
    for q in parts:
        if q.startswith("$") or q.startswith("`"):
            out.append(q)
            prev_end_dollar = q.startswith("$")
            continue

        def rep(m, _first=[prev_end_dollar]):
            num = m.group(1)
            start = m.start()
            # "$…$ (16.N)": the equation is already beside its number
            if start == 0 and _first[0] or q[:start].strip() == "" and _first[0]:
                return m.group(0)
            if num in skip or num in seen or num in shown or num not in EQ:
                return m.group(0)
            if re.match(r"\s*[–-]\s*\(16\.", q[m.end():]) or re.search(r"\(16\.\d+\)\s*[–-]\s*$", q[:start]):
                return m.group(0)
            seen.add(num)
            lab = m.group(2) or ""
            key = num
            if re.search(r"print(?:s|ed)\s*(?:the\s+)?$|print(?:s|ed)\s+\S+\s*$", q[:start]) and num + "p" in EQ:
                key, lab = num + "p", ", as printed"              # "the printed (16.54)" names the printed form
            return f"${EQ[key]}$ ({num}{lab})"
        out.append(_BARE.sub(rep, q))
        prev_end_dollar = False
    return "".join(out)


# ---------------------------------------------------------------------------------------------------------------------
# Part F of the design, parsed at build time: every derivation copied word for word (the ch07–ch15 pattern; the design
# writes the fields inline — "**Goal:** … · **Start:** … · **Plan:** • … • …" — and the steps as
# "n. **did** … · **tex** … · **why** … · **plain** …"; a ★★★ block carries "- **sympy check** (…):", a listing and a
# "printed:" line).
# ---------------------------------------------------------------------------------------------------------------------
_FIELD = re.compile(r"(?:^- |\s·\s+)\*\*(Goal|Start|Plan|Tools|Assumptions|Steps|Result|Check|What it means|Traps):\*\*\s*", re.M)


def _start_split(s: str) -> tuple[str, str]:
    """'maths … — "in words"; tail.' → (maths, words + tail); without the quoted part the whole text is the maths line."""
    m = re.search(r"\s—\s[\"“](.+?)[\"”]", s)
    if m:
        tail = s[m.end():].strip(" .;")
        words = m.group(1).strip() + ("; " + tail if tail else "")
        return s[:m.start()].strip(), words
    return s.strip(), ""


def part_f() -> dict[str, dict]:
    """Parse Part F into {D01: dict(title, stars, goal, start, plan, tools, assumptions, steps, result, check, meaning,
    traps, check_src)} — every field word for word."""
    part = DESIGN.split("## Part F", 1)[1].split("\n## Closing check", 1)[0]
    chunks = re.split(r"^### (D\d\d) · ", part, flags=re.M)
    out: dict[str, dict] = {}
    for key, body in zip(chunks[1::2], chunks[2::2]):
        code = re.search(r"```python\n(.*?)```", body, flags=re.S)
        check_src = textwrap.dedent(code.group(1)).strip("\n") if code else None
        # the listing, its "printed:" lines and the orchestrator's S34 note are handled separately (D27 below)
        body = re.sub(r"^- \*\*sympy check\*\*.*?```python\n.*?```[ \t]*\n(?:  (?!- ).*\n)*", "", body, flags=re.S | re.M)
        body = re.sub(r"^> \*\*Orchestrator note.*?\n\n", "", body, flags=re.S | re.M)
        lines = body.splitlines()
        title = _abs_plain(re.split(r" — ★", lines[0])[0].strip())
        stars = lines[0].split(" — ", 1)[1].count("★")
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
        flat = _FIELD.split("\n".join(rest))
        fields: dict[str, str] = {}
        for k in range(1, len(flat), 2):
            fields[flat[k]] = _abs_plain(_join(flat[k + 1].splitlines())).strip()
        parsed = []
        for st in steps_raw:
            s = _abs_plain(_join(st))
            bits = re.split(r"(?:^|\s·\s)\*\*(did|tex|why|plain)\*\*\s+", s)
            d = {name: val.strip().rstrip("·").strip() for name, val in zip(bits[1::2], bits[2::2])}
            assert set(d) == {"did", "tex", "why", "plain"}, (key, s[:120])
            nums = set(re.findall(r"\((16\.\d+)", d["tex"]))
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=_tidy(_expand_refs(d["why"], nums)),
                               plain=_tidy(_expand_refs(d["plain"], nums)), _nums=nums))
        start_tex, start_plain = _start_split(fields["Start"].rstrip(" ·"))
        res_tex, res_plain = _start_split(fields["Result"].rstrip(" ·"))
        plan = [p.strip().rstrip(".·").strip() + "." for p in fields["Plan"].split("•") if p.strip(" ·.")]
        tools = [_expand_refs(t.strip().rstrip(".·").strip()) for t in fields["Tools"].split(";") if t.strip(" ·.")]
        out[key] = dict(title=title, stars=stars, goal=_tidy(_expand_refs(fields["Goal"].rstrip(" ·"))),
                        start=(display_tex(start_tex), _tidy(start_plain)), plan=[_expand_refs(p) for p in plan],
                        tools=tools, assumptions=_tidy(_expand_refs(fields.get("Assumptions", "").rstrip(" ·"))),
                        steps=parsed, result=(display_tex(res_tex), _tidy(res_plain)),
                        check=_tidy(_expand_refs(fields.get("Check", ""))),
                        meaning=_tidy(_expand_refs(fields.get("What it means", ""))),
                        traps=_tidy(_expand_refs(fields.get("Traps", ""))), check_src=check_src)
    return out


PF = part_f()
assert len(PF) == 34, len(PF)
N_STEPS = sum(len(d["steps"]) for d in PF.values())
assert N_STEPS == 266, N_STEPS                              # Part F lists 266 steps (its closing line says 265: the per-derivation counts add to 266)
assert [k for k, d in PF.items() if d["check_src"]] == ["D13", "D19", "D20", "D21", "D24", "D27", "D30", "D33"]
assert [k for k, d in PF.items() if d["stars"] == 3] == ["D13", "D19", "D20", "D21", "D24", "D27", "D30", "D33"]
assert [len(PF[k]["steps"]) for k in sorted(PF)] == [6, 5, 4, 5, 9, 8, 6, 6, 4, 9, 7, 7, 13, 6, 8, 7, 7, 9, 10, 10, 13,
                                                      8, 3, 11, 6, 6, 14, 4, 9, 15, 5, 11, 10, 5]

if "--dump" in sys.argv:
    for _k, _d in PF.items():
        print(_k, _d["stars"], len(_d["steps"]), _d["title"])
        print("   start:", _d["start"])
        print("   result:", _d["result"])
        for _i, _s in enumerate(_d["steps"], 1):
            print(f"   {_i}. {_s['did']} | {_s['tex'][:90]}")
    print(sorted(EQ, key=lambda k: [float(x) if x.replace('.', '').isdigit() else 0 for x in re.findall(r'[\d.]+', k)][:1]))
    sys.exit(0)

for _k in ("16.15", "16.28", "16.29", "16.31", "16.131", "16.133"):
    assert _k in EQ, _k                                     # harvested from the derivation lines of Part F


def EQS(*nums: str) -> str:
    """Several numbered equations as one display block, one per line, each with its number."""
    rows = [EQ[n] + r"\qquad\text{(" + LABEL.get(n, n).replace("—", "-") + ")}" for n in nums]
    return rows[0] if len(rows) == 1 else r"\begin{array}{l}" + r" \\[6pt] ".join(rows) + r"\end{array}"


# Starts and results that Part F gives as bare numbers or in words: the same statement as display maths (the words of
# Part F stay in "In words").
START_TEX = {
    "D01": r"\text{force on a wall patch}=p\,dA_w,\qquad\text{patch moves out by }ds\ \Rightarrow\ dV=dA_w\,ds",
    "D05": EQS("16.21", "16.22", "16.23", "16.24", "16.25", "16.26"),
    "D06": EQS("16.27", "16.28") + r"\qquad b\equiv a-\delta",
    "D07": EQS("16.31", "16.32"),
    "D09": EQS("16.47"),
    "D11": EQS("16.53", "16.54", "16.55"),
    "D12": EQS("16.56", "16.57", "16.58"),
    "D14": EQS("16.82") + r"\qquad J_0(z)=\sum_{m\ge0}\dfrac{(-1)^m(z/2)^{2m}}{(m!)^2}",
    "D18": EQS("16.101", "16.102", "16.103", "16.114", "16.115"),
    "D20": EQS("16.123", "16.128", "16.122"),
    "D21": EQS("16.131", "16.133", "16.124", "16.125", "16.129"),
    "D22": r"P_{pump}=\Delta p\,Q,\qquad P_{up}=K\times(\text{volume of blood held})",
    "D26": r"\text{power of a travelling wave}=p\,Q,\qquad Q=Yp",
    "D27": EQS("16.173", "16.174", "16.176", "16.177"),
    "D31": r"\textstyle\sum F_x=0:\ \text{pressure on the two end faces of a cylinder }(r,\Delta x)+\text{shear on its side}=0",
    "D33": EQS("16.219", "16.215", "16.218"),
}
START_WORDS = {
    "D05": "Three conditions in the core and three in the plasma layer, with the signs of the two stress laws corrected (slip S2).",
    "D06": "The two layer profiles found in D05, with the shorthand b for the radius of the interface.",
    "D07": "The total flux of the two-layer model and the definition of an apparent viscosity by Poiseuille's formula.",
    "D09": "Inflow into the aorta = rate of change of the chamber's volume + outflow into the periphery.",
    "D11": "Area conservation, the corrected momentum law and the Hookean tube law for a blood-filled elastic tube.",
    "D12": "The three linear laws of D11: area, momentum and the linear tube law.",
    "D14": "Womersley's profile and the power series of J₀ (primer P374).",
    "D18": "Linearised continuity and the two momentum equations in (r, x), and the no-slip conditions at the moving wall.",
    "D20": "The mode equation for the axial velocity, the wall's pressure law of D19 and continuity in mode form.",
    "D21": "The mode solution of D20, the two wall conditions and the axial wall displacement of D19.",
    "D22": "The two costs of a segment — pumping blood through it and keeping its blood alive.",
    "D26": "A travelling wave's power is pressure × flux, with Q = Yp.",
    "D27": ("The axial and vorticity equations of the curved pipe, written with the advection sign of (16.174) corrected "
            "(slip S34, found when we re-derived it: eliminating the pressure from the two momentum equations gives +), the "
            "wall conditions and Dean's series in κ."),
    "D31": "A fluid cylinder of radius r and length Δx in steady, fully developed flow.",
    "D33": "Plug plus ring, with the ring profile and the plug speed of D32.",
}
RESULT_TEX = {
    "D05": EQS("16.27", "16.28"),
    "D06": EQS("16.29", "16.30", "16.31"),
    "D07": EQS("16.33", "16.34"),
    "D10": (r"\begin{array}{l}" + EQ["16.52"] + r"\qquad\text{(16.52)}\\[4pt]"
            r"p=p(T_s)\,e^{-(t-T_s)/R\mathrm K}\ \ (T_s\le t<T)\qquad\text{(diastole)}\\[4pt]"
            r"p_0=\dfrac{\frac1{\mathrm K}e^{-T/R\mathrm K}\int_0^TQ\,e^{\tau/R\mathrm K}d\tau}{1-e^{-T/R\mathrm K}},\qquad"
            r"\bar p=R\bar Q\qquad\text{(ours)}\end{array}"),
    "D11": EQS("16.56", "16.57", "16.58"),
    "D15": EQS("16.88"),
    "D16": EQS("16.90", "16.91"),
    "D17": EQS("16.99", "16.100") + r"\qquad\dfrac{|Q|}{|Q_0|}=\dfrac{8\sqrt{X^2+Y^2}}{\alpha^2}\ \text{(ours)}",
    "D18": EQS("16.120", "16.122", "16.123", "16.124", "16.125"),
    "D19": EQS("16.128", "16.129"),
    "D20": EQS("16.131", "16.133") + r"\qquad\beta=\sqrt{i\omega/\nu}",
    "D21": EQS("16.136") + r"\qquad\dfrac{c}{c_0}=\dfrac1{\mathrm{Re}\sqrt X},\ \ \dfrac{k_2}{k_1}=\dfrac{\mathrm{Im}\sqrt X}{\mathrm{Re}\sqrt X}",
    "D23": EQS("16.144"),
    "D24": EQS("16.146", "16.147") + r"\qquad\text{symmetric: }\theta=\arccos2^{-1/3}=37.47^\circ\ \text{(ours)}",
    "D25": EQS("16.153", "16.154"),
    "D27": (r"\begin{array}{l}\psi_0=\dfrac{r}{144}(1-r^2)^2(4-r^2)\sin\alpha,\qquad w_1=+\dfrac{r(1-r^2)(19-21r^2+9r^4-r^6)}{5760}"
            r"\cos\alpha\qquad\text{(ours)}\\[6pt]\dfrac{Q}{\pi a^2\bar W}=1-\underbrace{0.0306}_{=1541/50400}\left(\dfrac K{576}"
            r"\right)^2+\underbrace{0.0120}_{\text{printed; }0.0119312\text{ with S34}}\left(\dfrac K{576}\right)^4+O(K^6)"
            r"\qquad\text{(16.178)}\end{array}"),
    "D29": EQS("16.196", "16.197", "16.200", "16.201", "16.202"),
    "D30": EQS("16.204", "16.205", "16.206") + r"\qquad\text{power law: }\mathrm M=2-n=\tfrac12",
    "D31": EQS("16.209", "16.210") + r"\qquad\Delta p_{\min}=\dfrac{2L\tau_y}{a}\ \text{(ours)}",
    "D32": EQS("16.215", "16.216", "16.218"),
    "D33": EQS("16.220"),
}
RESULT_WORDS = {
    "D05": "Speeds match at the interface, slopes do not.",
    "D06": "The plasma layer's flux, the core's flux (corrected form) and their sum.",
    "D07": "The two-layer apparent viscosity and its thin-layer approximation.",
    "D10": "Today's pressure is the faded memory of every past stroke; between beats it decays with time constant RK; a "
           "steady heartbeat settles to one starting pressure; the mean is set by R alone.",
    "D11": "Linear laws for a small pulse: area, momentum and a linear spring for the wall.",
    "D15": "A flat plug lagging by a quarter cycle, joined to the wall by a Stokes layer.",
    "D16": "Womersley's oscillating flow rate, with F₂ the area average of the Bessel ratio.",
    "D17": "How much flow and how late, as functions of α only.",
    "D18": "Ordinary differential equations in r for one travelling mode, with the pressure uniform over each section.",
    "D19": "Hoop stiffness holds the pressure; the axial wall motion follows from Poisson coupling and the blood's drag.",
    "D20": "The axial and radial velocity of the mode.",
    "D21": "The characteristic quadratic; the pressure root gives the speed and damping of the viscous pulse.",
    "D23": "Murray's cube law.",
    "D24": "Branch angles fixed by the radii alone.",
    "D25": "Reflection and transmission set by the admittance mismatch.",
    "D27": ("The vortex pair, the outward-skewed axial correction and the flux coefficient 1541/50400 = 0.030575 of "
            "(16.178); the K⁴ coefficient is printed 0.0120, and 0.0119312 once the advection sign of (16.174) is "
            "corrected (slip S34)."),
    "D29": "Small changes of area ratio, pressure, wave speed and speed index along the tube.",
    "D30": "Shapiro's equations: the factor 1 − S² decides every response, and at S = 1 the drivers must balance.",
    "D31": "Linear stress, a plug of radius a_c, and a smallest pressure drop below which nothing flows.",
    "D32": "The velocity in the sheared ring and the speed of the plug.",
    "D33": "The Casson flux: Poiseuille's with μ = K_c², times a bracket below 1.",
}
for _k, _v in START_TEX.items():
    PF[_k]["start"] = (_v, START_WORDS.get(_k, PF[_k]["start"][1]))
for _k, _v in RESULT_TEX.items():
    PF[_k]["result"] = (_v, RESULT_WORDS.get(_k, PF[_k]["result"][1]))
PF["D27"]["result"] = (PF["D27"]["result"][0], _expand_refs(PF["D27"]["result"][1]))

# Part F's *why* of these steps is a bare formula or a word; each gets the sentence that says why the move is allowed and
# what it is for (nbkit requires both).  The text of Part F stays, word for word, in front.
WHY_ADD = {
    ("D06", 4): "Integrating a constant over the core disc is the ring sum of D04 step 1; we do it because the plug part moves as one piece.",
    ("D07", 3): "Taking reciprocals of two equal positive numbers keeps them equal, and it gives μ_app itself.",
    ("D07", 5): "Multiplying out a product is ordinary algebra (no approximation beyond step 4); it separates the core's share from the layer's.",
    ("D07", 6): "Factoring a common term out of a bracket changes nothing; it shows the core's viscosity as the reference value.",
    ("D08", 3): "Dividing both sides by the positive number 2h keeps the balance true and leaves the hoop stress alone.",
    ("D11", 6): "Replacing the square root by its expansion of step 5 keeps the law valid to first order; the result is linear in A′.",
    ("D13", 10): "Blood sticks to a rigid wall that does not move (no slip), so the profile must vanish at r = a; this fixes C₁.",
    ("D13", 11): "Renaming a group by its definition changes nothing; dividing by J₀(i^{3/2}α), which is never zero, solves for C₁.",
    ("D14", 2): "Squaring a product squares each factor; we need z² because the series of step 1 is in powers of z².",
    ("D14", 5): "Subtracting both sides from 1 keeps the equation true; 1 − F₁ is the factor that multiplies the profile.",
    ("D16", 2): "The power rule gives ∫₀ᵃ r dr = a²/2; it is the contribution of the constant 1 in 1 − F₁.",
    ("D16", 7): "Multiplying the constants together is ordinary algebra; it gives the compact form of the book.",
    ("D17", 4): "Multiplying out two complex numbers and grouping terms is exact; the real part is the physical flow.",
    ("D19", 9): "Adding η̂/a₀² to both sides of step 8 keeps it true; it builds the bracket of step 7 so that it can be replaced.",
    ("D20", 1): "Substituting one equation into another is allowed; the factors of i cancel and only a constant is left on the right.",
    ("D20", 10): "Integrating each constant and multiplying the whole by the same factor keeps the equality; it completes the integral of step 8.",
    ("D21", 8): "Rewriting a coefficient by an identity changes nothing; it puts the radial condition in the variables of step 7.",
    ("D24", 10): "Replacing a quantity by an equal one keeps the equation true; the angle then depends on a₀ and a₁ only.",
    ("D24", 11): "Setting a₁ = a₂ is a special case of the general result; the fourth powers then cancel exactly.",
    ("D25", 4): "Adding the same terms to both sides keeps the equation true; it gathers the unknown p_R on one side.",
    ("D25", 5): "Dividing both sides by the sum of admittances is allowed because the sum is positive; it leaves the ratio p_R/p_I.",
    ("D26", 2): "Applying the power formula of step 1 to the incident wave in the parent tube (impedance Z₁) is direct substitution.",
    ("D26", 6): "Expanding the square shows the numerator equals the denominator, and a fraction with equal top and bottom is 1.",
    ("D27", 8): "Factoring a polynomial changes nothing; the factored form shows at once that ψ₀ and its slope vanish at r = 1.",
    ("D28", 4): "Substituting the slope of step 3 into step 2 is direct; powers of the same base multiply by adding exponents.",
    ("D29", 2): "Differentiating a product is exact; it splits the pressure change into three causes.",
    ("D30", 1): "Replacing du by its equal from continuity is exact; it removes the speed change from the momentum law.",
    ("D30", 9): "Adding the two dα terms is exact; it shows that the speed index responds to the shape through M alone.",
    ("D30", 12): "Grouping like terms after multiplying out is exact; it isolates how the area driver acts on S.",
    ("D30", 14): "Adding three exact lines and dividing by a small length dx gives rates per metre, the book's form.",
    ("D31", 3): "Setting r = a in the linear stress of step 2 gives its value at the wall, where it is largest.",
    ("D32", 4): "Replacing τ and τ_y by their expressions in r (step 3) is direct substitution; the root of a product is the product of the roots.",
    ("D32", 6): "Blood at the rigid wall does not move (no slip), so u(a) = 0; this fixes the one constant.",
    ("D33", 2): "Constants come out of the integral, so only the bracket that depends on r stays inside.",
    ("D33", 4): "The power rule holds for fractional powers when r > 0; it produces the 2/7 behind the 16/7 of the result.",
    ("D33", 10): "Undoing the substitution is allowed because t and a_c/a correspond one to one; it returns to the book's variable.",
}
WHY_SET = {   # Part F's one-word "Elementary." becomes the rule it stands for
    ("D33", 5): r"The power rule, $\int r\,dr=r^2/2$ and $\int r^2dr=r^3/3$, evaluated between $a_c$ and $a$; it gives the "
                r"part of the yield correction that is linear in $a_c$.",
}
for (_k, _i), _txt in WHY_ADD.items():
    _st = PF[_k]["steps"][_i - 1]
    assert len(_st["why"].split()) < 6, (_k, _i, _st["why"])        # only the steps that need it
    _st["why"] = _st["why"].rstrip(" .") + ". " + _txt
for (_k, _i), _txt in WHY_SET.items():
    PF[_k]["steps"][_i - 1]["why"] = _txt
for _k, _d in PF.items():
    for _i, _st in enumerate(_d["steps"], 1):
        assert len(_st["why"].split()) >= 6, (_k, _i, _st["why"])
        _st.pop("_nums", None)

# Learner-facing corrections to Part F's text (the design keeps its pipeline shorthand; the notebook says it in plain words,
# and two numbers the lesson review recomputed are corrected here).  Every (old, new) pair must match exactly once or more.
PF_FIX: list[tuple[str, str, str]] = [
    ("D01", "Our schematic loop (header convention 10)", "Our schematic loop (the `heart` inputs above)"),
    ("D01", "line integral (knowledge/primers.md P35, reminder)", "line integral along a path (Ch. 1, reminder)"),
    ("D01", "loop area (P138, reminder)", "loop area (Ch. 5, reminder)"),
    ("D02", "square roots (P393 is later", r"squaring a sum, $(x+y)^2=x^2+2xy+y^2$"),
    ("D03", "R02's link to separation is not needed", "the Stokes-layer depth of Ch. 8 (recap above)"),
    ("D06", r"$=\pi\Delta p\,a^2\delta(\delta-2a)/(4\mu_pL)\neq0$ (analysis S3)",
     r"$=+\pi\Delta p\,a^2\delta(2a-\delta)/(4\mu_pL)>0$ — the printed form overstates the core's flux (slip S3)"),
    ("D07", "binomial approximation (knowledge/primers.md P26 and P342, reminder)",
     r"binomial approximation $(1+\varepsilon)^n\approx1+n\varepsilon$ (Taylor series, Ch. 1, reminder)"),
    ("D10", "the closed form of C.12 equals the RK4 march", "the closed periodic form (`ch16.windkessel_summary`) equals the RK4 march"),
    ("D15", "**expect** at α = 20, r = 0.9a the asymptote is 5 % below the exact profile (the curvature correction ~δ/r)",
     "**expect** at α = 20, r = 0.9a the asymptote's amplitude is about 0.1 % from the exact one (0.14 %; 1.3 % as a complex "
     "number, which also counts the phase), and the gap grows toward the wall (at 0.97a: 2.1 % complex, 1.1 % in amplitude), "
     "because the thin-layer form drops the wall's curvature"),
    ("D27", "**plus sign** (header correction)", "**plus sign** (the fastest fluid moves to the outer side)"),
    ("D27", "(the series diverges; the reliability criterion of E8 is ours)",
     "(the series diverges; our reliability test `ch16.dean_series_reliability` asks that the K⁴ term stay below a quarter of the K² term)"),
    ("D30", "(ours, Part C prm)", "(ours: the resistance per length that `ch16.collapsible_march` uses)"),
    ("D30", "(trap T10, the ch15 lesson)", "(trap T10, the same care as in Ch. 15's duct-flow derivations)"),
    ("D34", "A negative absolute-gauge pressure is a tension, held by the cohesion of water (trap T22).",
     "Here p is a gauge pressure; the absolute pressure p + p_atm turns negative — a true tension, held by the cohesion of "
     "water — once p < −0.1 MPa, which for ψ_soil = −0.1 MPa happens from about 13 cm above the soil (trap T22)."),
    ("D34", "(cavitation, embolism)", "(cavitation: a vapour bubble forms and the conduit fills with gas — an *embolism* — "
                                      "after which it no longer carries water)"),
]


def _fix_all(obj, old: str, new: str, hits: list[int]):
    """Replace ``old`` by ``new`` in every string inside a parsed derivation (dicts, lists, tuples)."""
    if isinstance(obj, str):
        if old in obj:
            hits[0] += 1
        return obj.replace(old, new)
    if isinstance(obj, dict):
        return {k: (v if k == "check_src" else _fix_all(v, old, new, hits)) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return type(obj)(_fix_all(v, old, new, hits) for v in obj)
    return obj


for _k, _old, _new in PF_FIX:
    _hits = [0]
    PF[_k] = _fix_all(PF[_k], _old, _new, _hits)
    assert _hits[0], f"PF_FIX for {_k} did not match: {_old!r}"
PF["D02"]["tools"] = [t for t in PF["D02"]["tools"] if not t.startswith("here only")]   # the tail of the split D02 tool

# The ★★★ Check fields carry the sympy listing inline (it renders as one run-on line); the listing runs, readable, in the
# next cell, so the Check text says what it does and what it prints instead.
_INLINE_SYMPY = re.compile(r"\s*-\s*\*\*sympy check \(`check_src`, executed[^)]*\):\*\*\s*```python.*?```\s*printed:\s*", re.S)
for _k, _d in PF.items():
    if _d["check_src"]:
        _new, _n = _INLINE_SYMPY.subn(" The sympy cell below re-runs the derivation's own steps (every line is commented); it prints ",
                                      _d["check"])
        assert _n == 1 and "```" not in _new, (_k, _d["check"][-300:])
        _d["check"] = _new

# The D21 listing names its quadratic `quad`, which would hide scipy's `quad` for the rest of the session.
PF["D21"]["check_src"] = re.sub(r"\bquad\b", "quadratic", PF["D21"]["check_src"])
# Every line of a check listing carries a comment: Part F leaves a few lines bare (matched by their start).
_BARE_COMMENT = {
    "res = sp.diff(u, r, 2)": "residual of (16.130): left side minus right side",
    "v133 = (-sp.I*k*A/beta": "(16.133) written out, to compare",
    "r, w, rho, k, E, h, a0, nu, nuk, A, eta, X, bh": "the symbols of D18–D21",
    "s = sp.symbols('s', positive=True)": "integration variable",
    "    ratio = complex(": "determinant over (16.136), scaled: should be the same constant",
    "a0, a1, a2, th = sp.symbols": "radii and the angle θ",
    "B = minimize(cost": "the cheapest branch point, found numerically (P387)",
    "eqs = sp.Poly(sp.expand((2*lap": "coefficients of each power of r must vanish (step 7)",
    "eqs = sp.Poly(sp.simplify": "the same matching for w1 (step 11)",
    "K = sp.symbols('K', positive=True)": "Dean's K",
    "dal, dA0, dpe, dz, dx, dKp": "differentials and coefficients of D29–D30",
    "rho = sp.symbols('rho', positive=True)": "density",
    "t204 = al*(S2*dA0/A0": "dα from (16.204), to compare",
    "rhs = (dA0/A0*(-2": "right side of (16.205), continued on the next line",
    "u = G/(4*Kc**2)*((a**2": "(16.215), continued on the next line",
    "xi = ac/a": "ξ = a_c/a",
}
for _key in [k for k, d in PF.items() if d["check_src"]]:
    _lines = PF[_key]["check_src"].splitlines()
    for _i, _ln in enumerate(_lines):
        if _ln.strip() and "#" not in _ln:
            _c = next(v for k_, v in _BARE_COMMENT.items() if _ln.startswith(k_))
            _lines[_i] = _ln + "   # " + _c
    PF[_key]["check_src"] = "\n".join(_lines)

STARS = {1: "★", 2: "★★", 3: "★★★"}


def _cap(s: str) -> str:
    """Capital first letter when the text starts with an ordinary lower-case word (never a symbol or a maths span)."""
    return s[:1].upper() + s[1:] if re.match(r"[a-z]{2,}(?=[\s,;:.—–'-])", s) else s


# ---------------------------------------------------------------------------------------------------------------------
# Small wrappers so that every block has the same shape (every text passes through X: «16.62» → the equation)
# ---------------------------------------------------------------------------------------------------------------------
def D(key: str, title: str, ref: str = "", notes: str = "", after: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then the notes it carries and its traps as a ⚠️ callout.
    ``title`` is our heading; ``ref`` = "16.62" shows the key equation in the heading."""
    d = PF[key]
    if ref in EQ:
        ref = f"{LABEL.get(ref, ref)}: ${EQ[ref]}$"
    goal = _cap(d["goal"]) + (f"\n\n**Assumptions.** {_cap(d['assumptions'])}" if d["assumptions"] else "")
    nb.derivation(key, f"{X(title)} `{key}` {STARS[d['stars']]}", ref=ref, goal=goal,
                  start=(d["start"][0], _cap(d["start"][1])), plan=[_cap(x) for x in d["plan"]], uses=d["tools"],
                  steps=d["steps"], result=(d["result"][0], _cap(d["result"][1])), interpret=_cap(d["meaning"]),
                  check=_cap(d["check"]), check_src=d["check_src"])
    if notes:
        nb.md(X(f"📝 **Notes taught inside `{key}`:** {notes}"))
    if after:
        nb.md(X(after))
    if d["traps"]:
        nb.md(f"> ⚠️ **Common confusion (traps in `{key}`):** {_cap(d['traps'])}")


def core(cid: str, title: str, question: str) -> None:
    nb.core(cid, f"{X(title)} `{cid}`", question=X(question))


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` starts with the Part E concept text (the ledger check matches it)."""
    nb.primer(term, f"`{pid}` · {X(text)}", code=textwrap.dedent(code).strip("\n") if code else None)


def recap(rid: str, title: str, text: str, where: str, core_id: str | None = None) -> None:
    """A RECAP row. ``nb.recap`` closes the current CORE block; ``core_id`` re-enters it (the design's re-entry rule)."""
    nb.recap(rid, f"{title} `{rid}`", X(text), where)
    if core_id:
        nb.current_core = core_id


def recap_md(title: str, where: str, text: str) -> None:
    """An id-less reminder of an earlier chapter's idea in the 🔁 recap format; it does not close the block."""
    nb.md(f"> 🔁 **Recap — {title}** ({where}). {X(text)}", tags=["recap"])


def note(nid: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    nb.note(f"**{nid}** {X(text)}", equation=equation, ref=ref)


def md(text: str) -> None:
    nb.md(X(text))


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + X(text))


def idea(sketch: str = "", words: str = "") -> None:
    sk = textwrap.dedent(sketch).strip("\n")
    body = ("```\n" + sk + "\n```") if sk else ""
    nb.md("#### The idea\n\n" + body + (("\n\n" if body else "") + X(words) if words else ""))


def whatif(text: str) -> None:
    nb.md(f"**What would change if…** {X(text)}")


def confusion(text: str) -> None:
    nb.md("> ⚠️ **Common confusion:** " + X(text))


def slip(k: str, printed: str, correct: str, why: str = "") -> None:
    """A printed-vs-correct box (book slip Sk, taught in corrected form)."""
    nb.md(f"> ⚠️ **slip {k} — the book prints** {X(printed)}; **the correct form is** {X(correct)}"
          + (f"\n>\n> {X(why)}" if why else ""))


def code(src: str, explain: str | None = None, tags: list[str] | None = None) -> None:
    nb.code(src, explain=X(explain) if explain else None, tags=tags)


def plotly(src: str, explain: str | None = None) -> None:
    nb.plotly(src, explain=X(explain) if explain else None)


def animation(src: str, explain: str | None = None) -> None:
    nb.animation(src, explain=X(explain) if explain else None)


def live_cell(src: str, explain: str | None = None) -> None:
    nb.live(src, explain=X(explain) if explain else None)


def scratch(src: str, explain: str) -> None:
    nb.code(src, explain=X(explain), tags=["from-scratch"])


def figure(src: str, see: str, read: str, change: str, explain: str | None = None) -> None:
    nb.figure(src, see=X(see), read=X(read), change=X(change), explain=X(explain) if explain else None)


def notes3(see: str, read: str, change: str) -> None:
    nb.figure_notes(see=X(see), read=X(read), change=X(change))


def worked(title: str, steps_md: str) -> None:
    nb.worked_example(title, X(steps_md))


EXPLAINER_NOTE = {   # one conventions sentence per explainer, about what that explainer actually shows
    "casson_blood_flow": "Stresses are positive magnitudes; τ_y and K_c² are illustrative (Whitmore's constants with our plasma viscosity); the plug radius a_c is a fraction of the pipe radius a.",
    "windkessel_heartbeat": "Pressures are transmural and gauge (venous pressure 0), shown in mmHg; K is the compliance (not Murray's or Dean's K); every input is illustrative.",
    "pulse_wave_artery": "The wall bulge is drawn exaggerated (the true radius change is stated in the readout); c is the wave speed, not the blood's speed.",
    "womersley_profiles": "Time convention e^{+iωt} (the rigid tube); profiles are in units of Poiseuille's centre speed Ga²/4μ for the same gradient amplitude.",
    "womersley_flow_phase": "|Q|/|Q₀| is the flow amplitude over Poiseuille's at the same gradient amplitude; φ is the lag of Q behind −∂p/∂x = G cos ωt.",
    "murray_branching": "Here K is the upkeep cost per volume of blood (W/m³), not the compliance; angles are measured from the extension of the parent.",
    "junction_reflection": "R and T are pressure ratios; the energy bars use R² and (ΣY_j/Y₁)T², which always add to 1.",
    "dean_curved_pipe": "α is the angle around the cross-section (α = 0 at the outer wall); the Dean-number chips switch between κ, K = 2κ² and D = 4√2κ; the flux shows both the printed 0.0120 and the corrected 0.0119312.",
    "collapsible_tube_choking": "α is the area ratio A/A₀ (not the Womersley number); P is negative when the tube is collapsed; the tube-law fit is for collapse only (α ≤ 1); the Starling resistor is our model.",
    "xylem_water_potential": "Pressures are gauge, in MPa; ψ < 0 means water is held under tension; the cavitation limit is illustrative.",
}


def explainer(slug: str, heading: str, why_static: str, tries: list[str]) -> None:
    why = (f"**Why interactive rather than a static figure:** {X(why_static)}\n\n"
           "Start with the **Walkthrough**; the **Explain** tab works every number out with your settings, and the "
           "**Derivation** tab (where there is one) steps through the same derivation as above. " + X(EXPLAINER_NOTE[slug]))
    nb.explainer(slug, heading=heading, why=why, tries=[X(t) for t in tries])


# ---------------------------------------------------------------------------------------------------------------------
# One-line reminders of tools primed in earlier chapters (Part E rows that point to knowledge/primers.md)
# ---------------------------------------------------------------------------------------------------------------------
def _ledger_rows() -> list[list[str]]:
    part = DESIGN.split("## Part E", 1)[1].split("## Part F", 1)[0]
    rows = [[x.strip() for x in ln.strip().strip("|").split("|")] for ln in part.splitlines() if ln.strip().startswith("|")]
    return [r for r in rows if len(r) >= 3 and not set("".join(r)) <= set("-: ") and r[0].lower() != "concept"]


LEDGER = _ledger_rows()

REMIND: dict[str, str] = {        # our own sentence about how the earlier tool is used HERE
    "line integral along a path": "∫ F·ds adds up a quantity along a path, piece by piece; here the path is the heart's pressure–volume loop and the quantity is p dV.",
    "loop area ½∮x × dx": "the area enclosed by a closed curve is a line integral around it; here it turns the PV loop into the work of one beat.",
    "pint quantities `UN.Q_`, `.to()`, `UN.check_dimensions`": "`UN.Q_(120, 'mmHg').to('Pa')` attaches a unit to a number and converts it; `UN.check_dimensions(q, '[pressure]')` tests what a formula's result measures — here it catches two printed slips.",
    "matplotlib figures, `np.linspace`, f-strings": "`np.linspace(a, b, n)` makes n evenly spaced values, `plt.subplots` makes a figure, `f\"{x:.3g}\"` prints a number with 3 significant figures — the three tools of every cell below.",
    "`animate` / `show_animation`, `slider_figure`, `show_viz`, `live`": "`slider_figure(fn, name, values)` precomputes a plotly figure for every slider position (it works on the web page); `animate` + `show_animation` play a matplotlib animation; `show_viz` embeds an explainer; `live` makes kernel-only sliders.",
    "`assert np.allclose` / `np.isclose`": "`assert np.allclose(a, b)` stops with an error if two results differ beyond a tiny tolerance; here it proves each from-scratch cell matches the tested `fluidpy` function.",
    "ordinary derivative as a slope": "dy/dx is the slope of the tangent to the curve y(x); here it is the stress needed per extra unit of shear rate.",
    "power laws and log–log axes": "on log–log axes y = Cxⁿ is a straight line of slope n; here the resistance a⁻⁴ shows as slope −4 and the apparent viscosity bends toward its asymptote.",
    "`np.where` for a piecewise law": "`np.where(cond, A, B)` picks A where the condition holds and B elsewhere, element by element; here it gives zero shear rate below the yield stress.",
    "definite integral; a ring of area 2πr dr": "∫ₐᵇ f dr adds f dr over a range; a thin ring of radius r and width dr has area 2πr dr, so a flux is ∫u 2πr dr.",
    "substitution in an integral": "replace the integration variable by a new one, converting the differential and the limits too; here s = a² − r² turns ring integrals of the two-layer profile into simple powers.",
    "`scipy.integrate.quad`": "`quad(f, a, b)` integrates a function numerically and returns (value, error estimate); here it re-adds a profile's rings as a check of the closed-form flux.",
    "sympy symbols, `diff`, `integrate`, `simplify`, `solve`, `expand`": "sympy does algebra with symbols: `sp.symbols` names them, `sp.diff`/`sp.integrate` differentiate and integrate, `sp.solve` solves, `sp.simplify(a - b) == 0` proves two expressions equal; every ★★★ derivation below is re-run with it.",
    "first-order Taylor expansion and the binomial approximation": "f(x + ε) ≈ f(x) + f′(x)ε, and (1 + ε)ⁿ ≈ 1 + nε for small ε; here a thin plasma layer turns (1 − δ/a)⁴ into 1 − 4δ/a.",
    "`scipy.optimize.brentq`": "`brentq(f, a, b)` finds the root of f between a and b where f changes sign; here it closes the two-layer model by finding the core's hematocrit.",
    "chain rule": "if v depends on t only through p(t), then dv/dt = (dv/dp)(dp/dt); here it turns the aorta's volume change into a pressure change.",
    "first-order linear ODE and the integrating factor": "y′ + Py = Q becomes d(y e^{∫P})/dx = Q e^{∫P} after multiplying by e^{∫P}, and one integration solves it; here y is the aortic pressure.",
    "product rule": "d(uv) = u dv + v du; here it shows that the integrating factor makes the left side an exact derivative, and it splits dp into three causes in a collapsible tube.",
    "`scipy.integrate.solve_ivp`": "`solve_ivp(rhs, (t0, t1), y0)` integrates dy/dt = rhs(t, y) numerically; here it marches the Windkessel and the collapsible tube.",
    "ipywidgets `live` and `_ = live(...)`": "`live(fn, a=(lo, hi, step), …)` draws sliders that re-run fn while a kernel runs; writing `_ = live(...)` keeps the cell from printing a stray object.",
    "limits and orders of smallness": "when a quantity ε is small, ε² is smaller still and is dropped beside ε; here products of two small disturbances are dropped to linearise the tube equations.",
    "Schwarz's theorem (mixed partials commute)": "for a smooth function ∂²u/∂t∂x = ∂²u/∂x∂t; here it lets the area law and the momentum law share one term so the velocity can be eliminated.",
    "the second-order wave equation and its two characteristic coordinates": "p_tt = c²p_xx is solved by any f(x − ct) + g(x + ct), two shapes sliding apart at ±c; the pulse in an artery is exactly this.",
    "Laplacian in cylindrical coordinates": "for a field that depends on r only, ∇²u = u″ + u′/r; it is the viscous term of every pipe-flow equation in this chapter.",
    "complex amplitudes; Euler's formula": "e^{iθ} = cos θ + i sin θ, so a complex number times e^{iωt} carries an amplitude and a phase; the physical quantity is the real part.",
    "complex square roots and the quadratic formula": "aX² + bX + c = 0 has the roots (−b ± √(b² − 4ac))/2a also for complex coefficients; the principal square root has a positive real part.",
    "binomial approximation with a complex small quantity": "1/(1 + ε) ≈ 1 − ε holds for a complex ε of small modulus too; here ε = iα²/4 at small Womersley number.",
    "linear second-order ODE: trying e^{λy}": "a constant-coefficient linear ODE is solved by e^{λy} with λ a root of the characteristic equation; here it gives the Stokes layer at the wall.",
    "Fourier modes": "a periodic signal is a sum of harmonics G_n e^{inωt}; for a linear equation each harmonic can be solved alone and the answers added.",
    "Crank–Nicolson check; `ch16.womersley_pde_solve`; caching": "Crank–Nicolson marches a diffusion equation in time with second-order accuracy; a slow result can be cached (`functools.lru_cache`) so it runs once per session — here the PDE route is an independent check of the Bessel profile.",
    "the complex plane in numpy; `np.abs`, `np.angle`": "`np.abs(z)` is the modulus and `np.angle(z)` the argument of a complex number; here they read amplitude and phase off Womersley's complex flow.",
    "`np.arctan2`": "`np.arctan2(y, x)` returns the angle of the point (x, y) in the correct quadrant; here it gives the phase lag from X and Y.",
    "phase of a wave; lag in time": "a cosine shifted by φ, cos(ωt − φ), peaks φ/ω seconds later; here the flow peaks after the pressure gradient.",
    "non-trivial solutions of a homogeneous 2 × 2 system": "two linear equations with zero right side have a solution other than zero only if their determinant vanishes; here that condition picks the allowed wavenumbers.",
    "the minimum of a falling plus a rising power": "a/xⁿ + bxᵐ has one minimum, where the derivative vanishes; at that point the two terms stand in the fixed ratio m : n — here pumping : upkeep = 1 : 2.",
    "level sets and the directional derivative": "the rate of change of f along a unit direction m is ∇f·m; here the change of a length when the branch point moves.",
    "the cosine rule": "in a triangle with sides p, q, s, s² = p² + q² − 2pq cos(angle between p and q); here three weighted pulls close a triangle and fix the branch angles.",
    "logarithmic differentiation": "d(ln f) = df/f turns a product of powers into a sum of relative changes; here c² = αK_pP′/ρ becomes 2dc/c = dα/α + dK_p/K_p + dP′/P′.",
    "differentials of algebraic laws as a small linear system": "differentiate each law, divide by it, and solve the small linear system for the changes; in Ch. 15 the common factor was 1 − M², here it is 1 − S².",
}


def _norm(c: str) -> str:
    return re.sub(r"[`*$\\\"“”]", "", c).lower().strip()


REMIND_N = {_norm(k): v for k, v in REMIND.items()}
_NO_ROW = sorted(k for k in REMIND_N if k not in {_norm(r[0]) for r in LEDGER})
assert not _NO_ROW, f"REMIND keys that match no Part E row: {_NO_ROW}"


def remind(cid: str, lead: str = "") -> None:
    """One 🔁 cell reminding the tools primed earlier that Part E lists as first used in ``cid``."""
    items = [(r[0], r[2]) for r in LEDGER if r[1] == cid and "primer" in r[2].lower() and not r[2].lower().startswith("primer (in")]
    if not items:
        return
    lines = []
    for concept, by in items:
        chs = re.findall(r"ch(\d\d)", by)
        pids = re.findall(r"P\d+a?", by)
        where = ", ".join(([f"Ch. {int(chs[0])}"] if chs else ["this chapter"]) + pids)
        text = REMIND_N.get(_norm(concept))
        assert text, f"no reminder sentence for ledger row {concept!r}"
        lines.append(f"- **{concept}** — {text} *({where})*")
        nb.primers.append(re.sub(r"[`*$\\]", "", concept) + " (reminder)")
    head = "> 🔁 **Tools from earlier chapters used in this block** (one line each; the full primers are in the chapters named)"
    nb.md(f"{head}{(' — ' + lead) if lead else ''}\n\n" + "\n".join(lines), tags=["primer"])


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=X(r"""
The last chapter turns the whole book on your own body. Blood is a slightly odd fluid, pushed in pulses through elastic,
branching, curved and sometimes collapsing tubes — and almost every answer is either Poiseuille's resistance «16.20» or a
wave whose speed the wall sets, «16.62», adjusted by one number, the Womersley number «16.12», that compares a vessel with
its own Stokes layer. At the end the same tools lift water up a tree."""),
    roadmap=[
        "**The heart as a pump** — its work per beat is the area of its pressure–volume loop (`C01`).",
        "**Blood thins when sheared** — Casson's law, a straight line in square roots (`C02`).",
        "**One number for every vessel** — the Womersley number α (`C03`).",
        "**Resistance ∝ a⁻⁴** — why small arteries control blood pressure (`C04`).",
        "**A plasma layer makes narrow tubes slippery** — the Fåhræus–Lindqvist effect (`C05`).",
        "**The tube law** — an artery is a spring (`C06`).",
        "**The Windkessel** — the aorta stores each beat and drains it (`C07`).",
        "**The pulse is a wave of wall stretch** — the Moens–Korteweg speed (`C08`).",
        "**Womersley's oscillating profile** — a breathing parabola or a lagging plug (`C09`).",
        "**How much flow, how late** — amplitude and phase lag of pulsatile flow (`C10`).",
        "**Viscous pulse waves** — slower and fading, by amounts ∝ 1/α (`C11`).",
        "**Murray's law** — Q ∝ a³ and fixed branch angles (`C12`).",
        "**Reflections at forks** — an admittance mismatch sends the pulse back (`C13`).",
        "**Dean's vortices in a bend** — two cells and a core flung outward (`C14`).",
        "**Collapsible tubes choke at S = 1** — flow limitation (`C15`).",
        "**The Casson plug** — a yield stress freezes the centre of a pipe flow (`C16`).",
        "**Water potential** — sap is pulled up a tree under tension (`C17`).",
    ],
    prerequisites=[
        X(r"Newton's viscosity law «1.3» and hydrostatics «1.8» (Ch. 1)"),
        "continuity and Navier–Stokes in cylindrical coordinates (Ch. 4)",
        X(r"Poiseuille flow «8.6» and Stokes' second problem «8.38» (Ch. 8)"),
        "boundary-layer separation (Ch. 9)",
        X(r"d'Alembert's solution «15.15» and the $1-M^2$ structure of duct flow (Ch. 15)"),
    ],
)
nb.explainer_index([
    ("casson_blood_flow", "Does a yield stress stop the blood?", "the centre freezes into a plug; the flow stops only when the wall stress drops below τ_y"),
    ("windkessel_heartbeat", "Why doesn't blood pressure fall to zero between beats?", "a stretchy aorta stores each stroke and drains it with time constant RK"),
    ("pulse_wave_artery", "Is the pulse at your wrist the blood arriving?", "a wave of wall stretch at metres per second; the blood barely moves"),
    ("womersley_profiles", "Is pulsing flow just Poiseuille flow sloshing?", "α turns a breathing parabola into a lagging plug"),
    ("womersley_flow_phase", "Does the flow follow the pressure?", "small α a damper, large α a mass: 1/α² and 90° late"),
    ("murray_branching", "Why isn't a wider vessel always better?", "pumping ∝ a⁻⁴ against upkeep ∝ a²: Q ∝ a³ and fixed angles"),
    ("junction_reflection", "When does a fork send the pulse back?", "an admittance mismatch, not size alone; energy always balances"),
    ("dean_curved_pipe", "What does a bend do to the flow?", "two vortices, a core flung outward, a flux loss of order K²"),
    ("collapsible_tube_choking", "Why can't you suck harder through a floppy straw?", "at S = 1 the flow stops hearing the outlet"),
    ("xylem_water_potential", "How does water climb a 30-metre tree?", "it is pulled: ψ falls from soil to leaf and p goes below zero"),
])
nb.setup()

md(r"""
## 📦 Chapter imports and the numbers we will use

Every worked number in this notebook uses **our own illustrative inputs**, written out in the cell below as small
dictionaries: a resting adult heart (`heart`), blood (`blood`), Whitmore's Casson constants (`casson`), a medium artery
(`artery`), an aorta (`aorta`), a Windkessel (`windkessel`), a 40-µm tube (`micro`), an arteriole, Murray's costs
(`murray`), a bend (`bend`), a collapsible vein (`vein`), a Casson pipe, a tree and a phloem sieve tube. None of them is a
number from the book. The module carries its own illustrative set too; it is printed once at the end of the cell — **both
sets are ours**; where they differ, the notebook uses the dictionaries below.""")
code(r"""
import numpy as np                                          # arrays and maths
import sympy as sp                                          # symbolic algebra (derivation checks)
import matplotlib.pyplot as plt                             # static figures and animations
import pandas as pd                                         # small tables
import math, functools                                      # factorials; lru_cache (remember a slow result)
from scipy.integrate import quad, solve_ivp                 # numerical integral; ODE integrator
from scipy.optimize import brentq, minimize                 # root in a bracket; minimum of a function of several variables
from scipy.special import jv, jve                           # Bessel functions J_n (and the scaled version for large arguments)
from fluidpy import ch16_biofluid_mechanics as ch16         # this chapter's tested functions (re-exports Ch. 8 and Ch. 15 helpers)
from fluidpy.core import laminar as LAM, acoustics as AC, similarity as SIM, thermo as TH, units as UN   # earlier chapters' tools
from fluidpy import ch01_introduction as ch01               # Ch. 1: Newton's law of viscosity, capillary rise
from fluidpy.core.interact import slider_figure, animate_figure, live   # plotly slider figures and kernel-only widgets
from fluidpy.core.anim import animate                       # builds a matplotlib animation from an update function
from fluidpy.core.style import COLORS, savefig              # house colours; save a figure to outputs/ch16/
import logging, warnings; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # silence harmless font messages

MMHG = 133.322387415                                    # Pa per mmHg (exact definition)
heart = dict(EDV=130e-6, ESV=55e-6, HR=66.0, BSA=1.9,   # m^3, m^3, beats/min, m^2 (illustrative adult, ours)
             pS=125*MMHG, pD=75*MMHG, systole_share=0.35,                  # systolic, diastolic [Pa]; share of a beat in systole
             p_min=5*MMHG, p_ED=10*MMHG, p_open=80*MMHG, p_peak=125*MMHG, p_ES=100*MMHG)   # loop corners [Pa]
blood = dict(rho=1060.0, mu=3.5e-3, mu_p=1.2e-3, T=310.15)       # = ch16.BLOOD (cited typical values): kg/m^3, Pa s, Pa s, K
casson = dict(k_c=1.53, s0=2.0)                          # Whitmore's constants in (16.7): k_c [-], sqrt(tau_y/mu_p) [s^-1/2]
artery = dict(a=2e-3)                                    # a medium artery radius [m]
aorta = dict(a0=12e-3, h=1.5e-3, E=0.4e6, L=0.3)          # radius [m], wall [m], Young's modulus [Pa], length [m]
windkessel = dict(R=1.45e8, K=1.2e-8)                    # peripheral resistance [Pa s/m^3], compliance [m^3/Pa]
micro = dict(a=20e-6, delta=2e-6, mu_c=4.0e-3)            # a 40-µm tube with a 2-µm plasma layer; core viscosity [Pa s]
arteriole = dict(a=25e-6, L=2e-3)                        # radius, length [m]
murray = dict(mu=3.5e-3, b_metabolic=1000.0, L=0.1)        # Pa s, W/m^3 (the book's K), m
bend = dict(delta=0.25)                                  # curvature ratio a/R_curv
vein = dict(a0=4e-3, Kp=200.0, n=1.5, L=0.1)              # m, Pa, -, m (a collapsible vein segment)
casson_pipe = dict(a=0.5e-3, L=0.5)                      # L/a = 1000
tree = dict(H=30.0, T=293.15, psi_soil=-0.1e6, c_leaf=300.0, a_xylem=25e-6, u_xylem=1e-3, mu_w=1.0e-3, nu_w=1.0e-6)   # m, K, Pa, mol/m^3, m, m/s, Pa s, m^2/s
phloem = dict(a=10e-6, r_pore=1e-6, L=300e-6, ell=1e-6, N=50, u=0.3e-3, mu=1.5e-3)   # sieve tube: m, m, m, m, -, m/s, Pa s
assert np.isclose(MMHG, ch16.MMHG) and blood == ch16.BLOOD  # our constants are the module's
mod = ch16.illustrative_inputs()                            # the module's own illustrative set (its scripts and tests use it; also ours)
print("module's own set (ours):", ", ".join(f"{k}={v:.4g}" for k, v in mod.items() if isinstance(v, float))[:600], "…")
""", explain=r"""
1. The first block imports numpy, sympy, matplotlib, pandas and four scipy tools, then the chapter module `ch16` and the
   earlier tools it reuses: `LAM` (Poiseuille flow and Stokes layers, Ch. 8), `AC` (acoustics, Ch. 15), `SIM`
   (dimensionless groups), `TH` (constants), `UN` (units with pint) and `ch01`. `ch16` re-exports the Ch. 8 Poiseuille and
   Stokes-layer functions and the Ch. 15 wave functions it uses, so every call below reads `ch16.<name>`.
2. The dictionaries are the inputs of every worked number. Pressures are stored in pascals (`125*MMHG`) and quoted in mmHg.
3. The assert proves that our blood properties and the mmHg factor are the module's; the last line prints the module's own
   (different) illustrative set once, so you know both exist.""")

md(r"""
## ⚠️ Conventions in this chapter — read once, return when a symbol surprises you

**(a) Here against elsewhere.** The book's conventions (and this notebook's), with the rival you may meet in other sources:

| # | What | Here | Elsewhere | Trap |
|---|---|---|---|---|
| i | time convention | rigid tube «16.73», «16.75»; elastic tube «16.116» | either sign | T1 |
| ii | pressure gradient | $G=-dp/dx>0$ drives flow in +x; «16.14» | signed $dp/dx$ (Ch. 8 code) | T4 |
| iii | wall stress | $\tau=-\mu\,du/dr>0$, a magnitude | signed $\tau=\mu\,du/dR$ (Ch. 8) | T5 |
| iv | Reynolds, Strouhal | $Re=2au/\nu$, $St=2a\omega/u$ (diameter) ⇒ $\alpha^2=Re\,St/4$ | $\alpha^2=2\pi\,Re\,St$ | T6 |
| v | Dean number | «16.171»; $K=2\kappa^2$; «16.179» | $De=Re\sqrt\delta$ (= κ) | T8 |
| vi | apparent viscosity | the secant $\tau/\dot\gamma$ | sometimes the tangent $d\tau/d\dot\gamma$ | T27 |
| vii | Casson constants | $k_c$ dimensionless in «16.6»; $K_c=k_c\sqrt{\mu_p}$ in (Pa s)$^{1/2}$ in «16.207» | $\tau^{1/2}=a\dot\gamma^{1/2}+b^{1/2}$ | T14 |
| viii | Shapiro's P | $\mathrm P=(p-p_e)/K_p<0$ when collapsed ($p_e$: the pressure outside the vessel, defined in C06); «16.188» | some plot $-\mathrm P$ | T9 |

No lapse rate appears in this chapter; the plants block uses only the hydrostatic term $\rho gz$.

**(b) Overloaded letters (trap T2).**

| Letter | Meanings in this chapter | Code name |
|---|---|---|
| α | Womersley number «16.12»; Shapiro's area ratio $A/A_0$ «16.186»; the angle round a curved pipe (Dean); $\alpha_c$ of «16.36» | `alpha`, `area_ratio`, `alpha_t`, `alpha_c` |
| K | compliance «16.48»; upkeep cost per volume «16.138»; Dean's $K=2\kappa^2$; $K_p$ bending stiffness; $K_c$, $k_c$ Casson | `K`, `b_metabolic`, `K_dean`, `Kp`, `Kc`, `k_c` |
| R | peripheral resistance «16.9»; radius of curvature; reflection coefficient «16.153»; resistance per length in «16.204»; gas constant «16.221» | `R`, `R_curv`, `refl`, `R_per_length`, `R_U` |
| δ | cell-free layer (C05); curvature ratio «16.164» | `delta`, `curvature` |
| ψ | secondary stream function «16.172»; water potential «16.221» | `psi_sec`, `psi` |
| c | wave speed; solute concentration in «16.221» | `c`, `c_solute` |
| T | period; temperature; tensions $T_\theta$, $T_x$; transmission coefficient «16.153» | `T_period`, `T`, `trans` |
| M, S, P | Shapiro's «16.206» (not a Mach number); speed index «16.195»; Shapiro's pressure «16.186» | `M_shapiro`, `S`, `P` |
| X, Y | parts of «16.97»; $X=(k/k_0)^2$ in «16.136» | `X`, `Y`, `Xk` |
| E, H, h, η, ν̂ | Young's modulus; hematocrit; wall thickness; wall displacement; Poisson's ratio (ν is kinematic viscosity) | `E`, `H`, `h`, `eta`, `poisson` |

**(c) Slips and traps.** *Slips* are false as printed; *traps* are true but easy to misread. **Twelve slips are taught corrected**
(eleven change a result or a number; the twelfth, S6, is a sign the book's own solution already gets right) — for example «16.140» and «16.76» — each in a box "the book prints … ; the
correct form is …" where it is used; the other 22 are named in one clause where they occur.

**(d) Which inputs?** Every worked number uses our own illustrative inputs (the dictionaries above), never the book's.

**(e) Which units?** SI in every function; mmHg, poise, dyn/cm², bar and the PRU only through converters (primer P394 in
C01).

**(f) Two signs in Dean's curved-pipe flow (C14).** The first axial correction $w_1$ of Dean's series has a **plus** sign, so
the fastest fluid moves to the outer side of the bend (D27); and the advection term of «16.174p» has the wrong sign as printed
(slip S34) — the corrected form is «16.174».""")
code(r"""
pd.set_option("display.max_colwidth", 110)                  # show long table cells without cutting them
display(ch16.conventions_table())                           # 18 rows: symbol, meaning, code name, where
slips = ch16.book_slips()                                   # 34 rows: id, kind, where, printed, correct, code, changes_physics
print(slips["kind"].value_counts().to_dict(), "| change the physics:", list(slips.loc[slips["changes_physics"], "id"]))
display(slips[["id", "kind", "where", "correct", "changes_physics"]])   # the corrected forms (the printed column is in the module)
display(ch16.traps())                                       # 27 rows: id, trap, guard (what the code does about it)
""", explain=r"""
1. `ch16.conventions_table()` is a pandas table (`pd.DataFrame`: rows and named columns) of the symbols and their code names.
2. `ch16.book_slips()` lists all 34 slips; the twelve of kind "slip" are taught in corrected form (`changes_physics` is True for the eleven that change a result) — the code keeps
   the printed form only as an option (`printed_…=True`) that a test must reject.
3. `ch16.traps()` lists the 27 traps with what the code does about each; each reappears as a "Common confusion" callout where
   it first bites.""")

# =====================================================================================================================
# A.1 §16.1 Introduction — N01
# =====================================================================================================================
nb.section("16.1", "Introduction", intro="**What is this section about?** Biofluid mechanics is the whole book applied to "
           "living things: blood in vessels from the aorta to capillaries, and water rising in plants. Nothing new is "
           "assumed beyond Chapters 1–15; what is new is the material (blood, elastic walls) and the questions.")
note("N01 [C]", r"""The scales run from a red cell (micrometres) to the aorta (centimetres) and the whole circulation
(metres). Every block below picks one scale: the heart and its loop (C01), the blood itself (C02), single vessels (C03–C11,
C16), branching and bending trees (C12–C15), and finally plants (C17).""")

# =====================================================================================================================
# A.2 §16.2 — C01, C02
# =====================================================================================================================
nb.section("16.2", "The Circulatory System in the Human Body", intro="**What is this section about?** The pump (the "
           "heart), the pipes (vessels) and the fluid (blood): the vocabulary of physiology turned into numbers we can "
           "compute with — stroke volume, cardiac output, peripheral resistance, and blood's odd viscosity.")
core("C01", r"The heart as a pump: «16.1» and the work $W=-\oint p\,dV$ (ours)",
     "How much blood does the heart move each beat, and how much work does that cost it?")
problem(r"""
Squeeze a water balloon through a narrow neck, let it refill, squeeze again. That is the left ventricle, about once a
second: it fills, closes its valves, squeezes until it can open the outlet valve, pushes blood into the aorta, closes,
relaxes. Draw its pressure against its volume over one beat and you get a closed loop. **The loop's width is the blood
moved, its area is the work done** — the same p dV work as a piston in Chapter 1.""")
idea(r"""
 p ▲
   │   D ◄──────────────── C      C→D ejection: the volume falls at high pressure
   │   │                   ▲
   │   ▼   area = work     │      B→C and D→A: both valves shut, the volume is fixed
   │   A ────────────────► B      A→B filling: the volume rises at low pressure
   └───┼───────────────────┼──► V
      ESV                 EDV     width = the blood moved per beat
""", r"""**Read two numbers off one loop.** The loop runs anticlockwise. Its width, from the emptiest volume ESV to the fullest
volume EDV, is the stroke volume «16.1»; its area is the work per beat, $W=-\oint p\,dV$ (ours) — positive, because the heart
pushes out at high pressure (top edge) the blood it took in at low pressure (bottom edge).""")
note("N02 [C]", r"""Our two-loop sketch of the circulation:

```
left heart ──► aorta ──► body (systemic) ──► veins ──► right heart ──► lungs (pulmonary) ──► left heart
```

A closed system: the same blood goes round both loops. The pulmonary loop returns in C06 (N40–N42), branching in C12 and
C13.""")
note("N04 [B]", r"""The four phases of one beat (A, B, C, D are the corners of the loop below):

| Phase | Valves | Volume | Pressure |
|---|---|---|---|
| filling A→B | mitral open, aortic shut | rises | low |
| isovolumic contraction B→C | both shut | fixed (EDV) | rises fast |
| ejection C→D | aortic open | falls | near the aortic pressure |
| isovolumic relaxation D→A | both shut | fixed (ESV) | falls fast |

**Systole** is the working part of the beat (contraction and ejection, B→C→D) and takes about a third of it; **diastole** is
the resting part (relaxation and filling, D→A→B), when the ventricle refills and the aorta is fed only by what it stored.
"End diastole" (corner B) is the moment the filling stops; "end systole" (corner D) the moment the ejection stops. The
half-sine ejection used in C07 lasts as long as systole (our `systole_share` = 0.35).""")
note("N03 [C]", r"""A resting adult heart beats about once a second; the cardiac cycle is one beat. The book quotes resting
values; our worked numbers use our own illustrative adult (the `heart` dictionary).""")
remind("C01")
P("P378", "the shoelace formula", r"""The area of a polygon from its corners $(x_i,y_i)$:
$A=\tfrac12\sum_i(x_iy_{i+1}-x_{i+1}y_i)$, positive when the corners run counterclockwise. It is the discrete form of the loop
area of Ch. 5 (P138) and of the line integral ∮y dx (P35). `np.roll(x, -1)` gives the next corner of every corner at once.""",
  code=r"""
x = np.array([0, 2, 2, 0]); y = np.array([0, 0, 1, 1])     # a 2 × 1 rectangle, corners counterclockwise
print(0.5*np.sum(x*np.roll(y, -1) - np.roll(x, -1)*y))      # 2.0 (clockwise order would give −2.0)
""")
P("P394", "units of physiology", r"""Clinicians measure pressure in mmHg (1 mmHg = 133.322387415 Pa, the weight of a 1 mm
mercury column), viscosity in poise or centipoise (1 cP = 10⁻³ Pa s), stresses in dyn/cm² (0.1 Pa), plant pressures in bar
(10⁵ Pa) and resistance in PRU (1 mmHg per cm³/s = 1.333 × 10⁸ Pa s m⁻³). Every function here works in SI; converters do
the rest (trap T3).""", code=r"""
print(ch16.mmHg_to_Pa(120), ch16.cP_to_Pa_s(3.5), ch16.dyn_cm2_to_Pa(1.0), f"{ch16.PRU:.4g}")   # 15999 Pa, 0.0035 Pa s, 0.1 Pa, 1.333e+08
print(UN.Q_(120, "mmHg").to("Pa"))                          # pint agrees (to its own definition of mmHg)
""")
D("D01", "Work per beat as the loop area", ref="")
note("N06 [B]", r"""Stroke volume «16.1»: the loop's width. Ours: 130 − 55 = 75 mL.""")
note("N07 [B]", r"""«16.2», «16.3», «16.4» — the unit chain written out: 75 mL/beat × 66 beats/min = 4950 mL/min =
4.95 L/min = 8.25 × 10⁻⁵ m³/s; per square metre of body surface (1.9 m²) 2.61 L/(min m²). (The book once writes "mL/m"
for mL/min — slip S32.)""")
code(r"""
SV = ch16.stroke_volume(heart["EDV"], heart["ESV"])        # (16.1): SV = EDV − ESV [m³]
EF = ch16.ejection_fraction(SV, heart["EDV"])               # (16.2) as a FRACTION (×100 for per cent)
CO = ch16.cardiac_output(SV, heart["HR"])                   # (16.3): SV × HR, HR passed in beats/min → m³/s
CI = ch16.cardiac_index(CO, heart["BSA"])                   # (16.4): CO per m² of body surface [m³/(s m²)]
print(f"SV = {SV*1e6:.1f} mL, EF = {100*EF:.1f} %, CO = {CO:.3e} m³/s = {CO*6e4:.2f} L/min, CI = {CI*6e4:.2f} L/(min m²)")
""", explain=r"""
1. `stroke_volume` subtracts the two volumes of the loop; `ejection_fraction` returns 0.577 (a fraction — 57.7 %).
2. `cardiac_output` multiplies by the heart rate in beats per **minute** and divides by 60 inside, so the result is in m³/s;
   × 6 × 10⁴ turns m³/s into L/min (10³ L per m³ × 60 s per min).
3. **What does this show?** 75 mL per beat at 66 beats/min is about 5 L of blood a minute — the whole blood volume of an
   adult, once a minute.""")
worked("a rectangular loop, by hand", r"""
1. Take a loop 100 mmHg tall and 75 mL wide.
2. Area = 100 × 75 = 7500 mmHg·mL.
3. One mmHg·mL = 133.32 Pa × 10⁻⁶ m³ = 1.333 × 10⁻⁴ J (P394).
4. W = 7500 × 1.333 × 10⁻⁴ = **1.000 J**.
5. At 66 beats per minute (1.1 per second) that is **1.1 W** of hydraulic power — about as much as a small torch bulb.""")
code(r"""
corners = {k: heart[k] for k in ("p_min", "p_ED", "p_open", "p_peak", "p_ES")}   # the five pressures that shape the loop [Pa]
loop = ch16.pv_loop_schematic(heart["EDV"], heart["ESV"], **corners)   # our schematic loop (ours): closed polygon A→B→C→D
w = ch16.pv_loop_work(loop["p"], loop["V"])                 # W = −∮p dV by the trapezoid rule along the polygon [J]
W_exact = ch16.pv_loop_schematic_work(heart["EDV"], heart["ESV"], **corners)   # the same work from the legs' closed forms [J]
P_card = ch16.cardiac_power(w["work"], heart["HR"])         # work per beat × beats per second [W]
print(f"W = {w['work']:.4f} J (exact {W_exact:.5f} J), orientation: {w['orientation']}, power = {P_card:.2f} W")
""", explain=r"""
1. `pv_loop_schematic` builds our loop (ours): filling A→B a gentle parabola from $p_{\min}$ to $p_{ED}$, the isovolumic
   rise B→C at EDV, ejection C→D a parabola through $p_{open}$, $p_{peak}$ and $p_{ES}$, and the isovolumic fall D→A at ESV.
2. `pv_loop_work` adds $-p\,dV$ around the polygon (trapezoids) and reports the orientation.
3. `pv_loop_schematic_work` integrates the parabolas exactly; the two agree to 4 digits.
4. `cardiac_power` multiplies by beats per second. **What does this show?** About one joule per beat, a little over a watt.""")
scratch(r"""
V, p = loop["V"], loop["p"]                                   # the loop's corners in order [m³], [Pa]
area = 0.5*np.sum(V*np.roll(p, -1) - np.roll(V, -1)*p)        # shoelace (P378): positive for a counterclockwise loop in (V, p)
assert np.isclose(area, w["work"], rtol=1e-9)                 # the area IS the work
rect_V = np.array([55e-6, 130e-6, 130e-6, 55e-6])             # a rectangle 75 mL wide …
rect_p = np.array([0.0, 0.0, 100*MMHG, 100*MMHG])             # … and 100 mmHg tall, counterclockwise in (V, p)
W_rect = ch16.pv_loop_work(rect_p, rect_V)["work"]            # its work by the library
assert np.isclose(W_rect, 75e-6*100*MMHG, rtol=1e-12)         # = Δp ΔV = 1.000 J (the tiny example)
assert np.isclose(ch16.pv_loop_work(rect_p[::-1], rect_V[::-1])["work"], -W_rect)   # run it backwards: −W
print(f"shoelace area = {area:.6f} J; rectangle = {W_rect:.4f} J")
""", explain=r"""
1. The shoelace sum of the loop's corners equals the library's work to 10⁻⁹: the work per beat **is** the loop's area.
2. The rectangle of the tiny example gives 1.000 J, and the same rectangle traced the other way gives −1.000 J (trap T11:
   the sign belongs to the direction).""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.0), gridspec_kw=dict(width_ratios=[1.6, 1]))   # loop and the rectangle
Vml, pmm = loop["V"]*1e6, loop["p"]/MMHG                         # mL and mmHg for reading
ax.fill(Vml, pmm, color=COLORS["amber"], alpha=0.25, label=f"area = W = {w['work']:.2f} J")   # the work, shaded amber
ax.plot(np.r_[Vml, Vml[0]], np.r_[pmm, pmm[0]], color=COLORS["accent"], lw=2)   # the loop itself
names = {"A": "A: valves switch\n(filling starts)", "B": "B: end diastole", "C": "C: aortic valve opens", "D": "D: end systole"}   # corner labels (diastole/systole: N04 above)
offs = {"A": (-8, 0), "B": (6, 6), "C": (6, 6), "D": (-6, 6)}   # label offsets in points: A sits left of its corner, clear of the loop
for k, (Vc, pc) in loop["corners"].items():                      # label the four corners
    ax.plot(Vc*1e6, pc/MMHG, "o", color=COLORS["ink"])           # corner dot
    ax.annotate(names[k], (Vc*1e6, pc/MMHG), xytext=offs[k], textcoords="offset points",   # the corner's name beside its dot (text only, no arrow)
                ha="left" if k in "BC" else "right", va="center" if k == "A" else "baseline", fontsize=8)
for (x0, y0, x1, y1) in [(110, 7, 120, 8.5), (130, 40, 130, 60), (95, 124, 85, 122), (55, 60, 55, 40)]:   # direction arrows
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=COLORS["accent"]))   # an arrow with no text: the direction round the loop
ax.annotate("", xy=(130, -8), xytext=(55, -8), arrowprops=dict(arrowstyle="<->", color=COLORS["blue"]))   # width = SV
ax.text(92.5, -6, "SV = 75 mL", ha="center", va="bottom", color=COLORS["blue"], fontsize=9)
ax.set(xlabel="left-ventricle volume V [mL]", ylabel="pressure p [mmHg]", ylim=(-14, 140), xlim=(25, 150),
       title="One beat traces a counterclockwise loop (ours)")
ax.legend(loc="center", fontsize=8)
bx.fill([55, 130, 130, 55], [0, 0, 100, 100], color=COLORS["amber"], alpha=0.25)   # the tiny example's rectangle
bx.plot([55, 130, 130, 55, 55], [0, 0, 100, 100, 0], color=COLORS["accent"])   # its outline: 75 mL wide, 100 mmHg tall
bx.text(92.5, 50, "75 mL × 100 mmHg\n= 1.00 J", ha="center", fontsize=9)
bx.set(xlabel="V [mL]", ylabel="p [mmHg]", xlim=(40, 150), ylim=(-14, 140), title="The tiny example")
savefig(fig, "ch16", "c01_pv_loop"); plt.show()
""", see=r"""Left: our schematic pressure–volume loop of the left ventricle, traced counterclockwise A → B → C → D, its
inside shaded amber; the blue arrow under it is the stroke volume. Right: the rectangle of the tiny example.""",
    read=r"""The width (EDV − ESV) is the blood moved; the height is the pressure the ventricle must reach to open the aortic
valve; the shaded area is the work per beat, 1.07 J here. The two vertical sides (isovolumic phases) do no work although
the pressure changes most there — nothing moves.""",
    change=r"""…the aorta were stiffer (C07)? The pressure during ejection would rise, lifting the top of the loop and its
area: the heart would work harder to push the same stroke volume.""")
animation(r"""
nfr = 24 if FAST else 40                                      # frames for one beat
idx = np.linspace(0, len(loop["V"]) - 1, nfr).astype(int)     # points of the loop shown in turn
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.5, 3.4), gridspec_kw=dict(width_ratios=[1.6, 1]))
ax.plot(Vml, pmm, color=COLORS["grid"], lw=2)                 # the whole loop, faint
(trace,) = ax.plot([], [], color=COLORS["accent"], lw=2.5)    # the part traced so far
(dot,) = ax.plot([], [], "o", color=COLORS["orange"], ms=8)   # where the beat is now
ax.set(xlabel="V [mL]", ylabel="p [mmHg]", xlim=(40, 150), ylim=(-5, 140), title="The loop, traced over one beat (ours)")
bx.set(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2), aspect="equal"); bx.axis("off")   # the ventricle as a breathing chamber
chamber = plt.Circle((0, 0), 1.0, color=COLORS["rose"], alpha=0.4); bx.add_patch(chamber)
label = bx.text(0, -1.15, "", ha="center", fontsize=9)
def phase(k):                                                 # which valves are open at loop point k
    V_, p_ = loop["V"][k], loop["p"][k]                       # volume [m³] and pressure [Pa] at that point
    if np.isclose(V_, heart["EDV"]): return "both valves shut: contraction"   # right vertical leg (volume fixed at EDV)
    if np.isclose(V_, heart["ESV"]): return "both valves shut: relaxation"    # left vertical leg (volume fixed at ESV)
    return "aortic valve open: ejection" if p_ > 40*MMHG else "mitral valve open: filling"
def update(i):                                                # draw frame i
    k = idx[i]                                                # which loop point this frame shows
    trace.set_data(Vml[:k + 1], pmm[:k + 1]); dot.set_data([Vml[k]], [pmm[k]])
    chamber.set_radius((loop["V"][k]/heart["EDV"])**(1/3))    # chamber size ∝ V^(1/3)
    label.set_text(phase(k))
    return trace, dot, chamber, label
fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it (faster frames)
show_animation(animate(update, frames=nfr, fig=fig, interval=80), player="video")
""", explain=r"""**What you see.** A dot runs round the loop over one beat while the chamber on the right swells (filling) and
shrinks (ejection); the line under the dot says which valves are open. **How to read it.** The chamber only changes size
on the bottom and top legs — on the vertical legs it stays the same size while the pressure jumps. **What would change
if…** a bigger stroke volume? The loop widens, the chamber's swing grows, and the area — the work — grows with it.""")
note("N08 [B]", r"""Pressures are quoted in mmHg; our illustrative 125/75 mmHg (systolic/diastolic, ours) = 16.7/10.0 kPa.
**N09 [B]** The pulse pressure $p_S-p_D$ = 50 mmHg is the size of the pressure wave; C07 shows how compliance sets it, C08 how
it travels.""")
code(r"""
print(ch16.mmHg_to_Pa(np.array([125, 75]))/1e3, "kPa")       # systolic and diastolic in kPa: 16.67, 10.0
print(ch16.pulse_pressure(heart["pS"], heart["pD"])/MMHG, "mmHg pulse pressure")   # p_S − p_D = 50 mmHg
""")
P("P381", "the electrical analogy", r"""A pressure difference drives a flow through a resistance as a voltage drives a
current: Δp ↔ V, Q ↔ I, R = Δp/Q ↔ Ohm's law. A vessel that stores volume when pressurised is a capacitor (compliance
K = dV/dp, C07); a tube's wave admittance Y = Q/p is the inverse of an impedance (C13). Resistors in series add; in parallel
their inverses add.""", code=r"""
R = np.array([1.0, 2.0, 4.0])                  # three resistances (any units)
print(R.sum(), 1/np.sum(1/R))                  # series 7.0, parallel 0.571
""")
note("N10 [B]", r"""Mean pressure falls continuously from the aorta to the veins: each vessel class is a resistance in
series, $\Delta p_i=QR_i$ (ours). Below, our illustrative shares of the total resistance (labelled ours).""")
figure(r"""
classes = ["aorta + large arteries", "small arteries", "arterioles", "capillaries", "venules + veins"]   # five classes in series
share = np.array([0.10, 0.15, 0.45, 0.20, 0.10])             # our illustrative shares of the total resistance
R_tot = 92*MMHG/CO                                            # total resistance that drops a 92 mmHg mean to 0 at our CO [Pa s/m³]
p_after = ch16.series_pressure_drop(list(share*R_tot), CO, 92*MMHG)   # pressure after each class [Pa]: p_in − Q ΣR_i
fig, ax = plt.subplots(figsize=(7.5, 3.4))
ax.bar(range(6), np.r_[92, p_after/MMHG], color=COLORS["orange"], alpha=0.8)   # pressure entering, then after each class
ax.set_xticks(range(6), ["entering"] + [c.replace(" + ", "\n+ ") for c in classes], fontsize=8)
ax.set(ylabel="mean pressure [mmHg]", title="The arterioles take the biggest drop (illustrative shares, ours)")
for i, v in enumerate(np.r_[92, p_after/MMHG]): ax.text(i, v + 1.5, f"{v:.0f}", ha="center", fontsize=8)
savefig(fig, "ch16", "c01_series_pressure"); plt.show()
""", see=r"""Orange bars: the mean pressure entering the arteries and after each vessel class, for our illustrative shares
of the total resistance.""",
    read=r"""Each drop is the flow times that class's resistance, $\Delta p_i=QR_i$; the tallest step is across the
arterioles, which carry 45 % of the resistance in our split.""",
    change=r"""…the arterioles narrowed by 10 %? Their resistance would grow by $1/0.9^4=1.52$ (C04): that is where the
body controls its blood pressure — C04 says why.""")
note("N11 [B]", r"""Continuity across the tree: the mean speed is $\bar u=Q/A_{tot}$ (ours). The capillaries are slow because
there are so many of them that their total area is huge.""")
recap_md("Mass conservation", "Ch. 4", r"""What goes in comes out: for an incompressible fluid the flow rate Q is the same
through every level of the tree, so speed × total area is the same everywhere.""")
figure(r"""
lev = ["aorta", "large arteries", "arterioles", "capillaries", "venules", "veins", "venae cavae"]   # the levels of the tree, in flow order
A_tot = np.array([4.5, 20, 40, 3000, 250, 80, 18])*1e-4      # total cross-section of each level [m²] (illustrative orders of magnitude, ours)
u_bar = ch16.mean_speed(CO, A_tot)                            # ū = Q/A_tot [m/s]
fig, ax = plt.subplots(figsize=(7.5, 3.4)); bx = ax.twinx()  # area on the left axis, speed on the right
ax.bar(range(7), A_tot*1e4, color=COLORS["teal"], alpha=0.5, label="total area [cm²]")
bx.plot(range(7), u_bar*100, "o-", color=COLORS["blue"], label="mean speed [cm/s]")
ax.set_yscale("log"); bx.set_yscale("log")
ax.set_xticks(range(7), lev, fontsize=8, rotation=20)
ax.set(ylabel="total area [cm²] (teal)", title="Speed × area is the same at every level")
bx.set_ylabel("mean speed ū [cm/s] (blue)")
savefig(fig, "ch16", "c01_continuity"); plt.show()
print("ū [cm/s]:", np.round(u_bar*100, 3))                    # 18 cm/s in the aorta, 0.03 cm/s in the capillaries
""", see=r"""Teal bars: total cross-section of each level of the tree (log scale, our orders of magnitude); blue dots: the mean
speed $\bar u=Q/A_{tot}$ at our cardiac output.""",
    read=r"""The two curves are mirror images on log axes because their product is the same Q everywhere: the capillaries,
with by far the largest total area, are about 700 times slower than the aorta — time for exchange with the tissue.""",
    change=r"""…the cardiac output doubled in exercise? Every blue dot would move up by a factor 2; the teal bars would
not move.""")
note("N12 [B]", r"""Total peripheral resistance «16.9» — the Ohm's law of the whole body; unit PRU = mmHg per cm³/s.
**N13 [B]** «16.10» with $\bar p_V$ taken as 0 (trap T12). **N14 [B]** «16.11»: an estimate of the mean for a typical
waveform (trap T15).""")
code(r"""
pA = ch16.mean_arterial_pressure(heart["pS"], heart["pD"])  # (16.11): the ⅓-rule mean [Pa]
R_body = ch16.peripheral_resistance(pA, CO)                 # (16.9) with (16.10), venous pressure 0 [Pa s/m³]
print(f"mean ≈ {pA/MMHG:.1f} mmHg; R = {R_body:.3e} Pa s/m³ = {R_body/ch16.PRU:.2f} PRU (ours)")
""")
confusion(r"""(trap T15) «16.11» is a rule of thumb for a typical waveform. C07 computes the true time mean of a model
waveform — the ⅓ rule is about 5 % low there.""")
note("N05 [C]", r"""Only a small fraction of the heart's metabolic power becomes pressure–volume work; most becomes heat. The
work per beat above (≈ 1 J) times the heart rate is the useful part.""")
whatif(r"""…the aorta were a rigid pipe? Then the pressure would follow the heart's ejection exactly and fall to zero between
beats. It does not — because the aorta stretches. C07 (the Windkessel) makes that precise; first, the fluid.""")

# ---------------------------------------------------------------------------------------------------------------------
# C02 — Casson
# ---------------------------------------------------------------------------------------------------------------------
core("C02", r"Blood thins when sheared: the Casson law «16.6»", "Is blood just thick water?")
problem(r"""
Honey is thick, water thin, and each has one viscosity whatever you do to it. Blood is different: stir it slowly and it
seems very thick, stir it fast and it thins. At rest its red cells stack into coin-like piles (rouleaux) that resist;
shearing breaks them up. And a tiny stress is needed before it flows at all — a **yield stress**.""")
idea(r"""
 √τ ↑        blood: a straight line that starts at √τ_y
    |       ╱
    |     ╱    ╱  Newtonian: a line through 0
    |   ╱    ╱
 √τ_y ╱   ╱
    | ╱
    +──────────────► √γ̇
""", r"""**A straight line in square roots, with an intercept** — that is Casson's law.""")
note("N15 [B]", r"""Blood in numbers (cited typical values, `ch16.BLOOD`): ρ = 1060 kg/m³, high-shear viscosity 3.5 mPa s,
plasma 1.2 mPa s, body temperature 310.15 K. **N16 [B]** A red cell is a flexible disc with far more surface than a sphere of
the same volume — a 90 µm³ cell would need only 97 µm² as a sphere (ours) — so it can fold through capillaries narrower than
itself (C05). **N17 [B]** Hematocrit H: the volume fraction of red cells, about 0.45 (45 %). Fraction in code, per cent in
prose (trap T13).""")
code(r"""
print(pd.Series(ch16.BLOOD).to_dict())                      # rho [kg/m³], mu [Pa s], mu_p [Pa s], T [K]
print(f"{ch16.equal_volume_sphere_area(90e-18)*1e12:.0f} µm²")   # area of a sphere holding 90 µm³: (36π V²)^(1/3) = 97 µm²
try:                                                        # a hematocrit passed as a per cent by mistake …
    ch16._check_hematocrit(45)                              # 45 instead of 0.45: the module's guard raises an error
except ValueError as e:                                     # … is refused with a message (trap T13)
    print("refused:", e)
""", explain=r"""
1. `pd.Series(...).to_dict()` prints the blood constants as a plain dictionary.
2. `equal_volume_sphere_area` gives the area of a sphere holding a red cell's volume: 97 µm², far less than the cell's.
3. `try` / `except ValueError` (Python's error handling, primed in earlier chapters): Python runs the `try` block; if a call
   inside it raises a `ValueError`, the notebook does not stop — it jumps to the `except` block, where `e` is the error and
   its message. Here the guard refuses a hematocrit of 45 (a per cent) because the code expects the fraction 0.45.""")
note("N18 [B]", r"""Non-Newtonian: the ratio of stress to shear rate is not a constant. The units box:
$1\,\mathrm{cP}=10^{-2}\,\mathrm P=10^{-2}\,\mathrm{dyn\,s/cm^2}=10^{-3}\,\mathrm{Pa\,s}$.""")
slip("S27", r"the chain 1 cP = 0.01 P = 1 dyn s/cm² = 0.1 N s/m²",
     r"$1\,\mathrm{cP}=10^{-2}\,\mathrm P=10^{-2}\,\mathrm{dyn\,s/cm^2}=10^{-3}\,\mathrm{Pa\,s}$",
     "The printed chain is 100 times too large from its third member on (1 P = 1 dyn s/cm² = 0.1 Pa s).")
code(r"""
print(ch16.cP_to_Pa_s(1.0), UN.Q_(1, "cP").to("Pa*s"))    # 0.001 Pa s twice: the corrected chain
""")
recap("R01", "Newton's law of viscosity", r"""«1.3»: stress proportional to shear rate, one constant μ. **New here:** for blood
τ/γ̇ is not constant — the Newtonian line becomes the ghost every figure of this block compares with.""",
      where="Ch. 1 §1.5 (C12)", core_id="C02")
code(r"""
print(ch01.newton_shear_stress(3.5e-3, 100.0), ch16.newtonian_stress(100.0, 3.5e-3))   # (1.3): 0.35 Pa twice (μ = 3.5 mPa s, γ̇ = 100 s⁻¹)
""")
note("N19 [B]", r"""The yield stress: below τ_y blood does not shear at all. An empirical correlation «16.5» — with **hidden
units said aloud**: τ_y in dyn/cm², H a fraction above 0.1, fibrinogen $C_F$ in g/100 mL. The book's quoted range at H = 0.45
does not quite match what «16.5» gives over its stated $C_F$ range (slip S29).""")
code(r"""
ty_corr = ch16.yield_stress_correlation(0.45, 0.3)          # (16.5) at H = 0.45, C_F = 0.3 g/100 mL, converted to Pa
print(f"tau_y = {ty_corr:.3e} Pa = {ty_corr/0.1:.4f} dyn/cm²")   # 7.84e-3 Pa = 0.0784 dyn/cm²
""")
note("N20 [C]", r"""Why it thins: rouleaux break up as the shear rate rises; above a threshold the cells move singly. Drawn as two
shaded regions on the figure below; the asymptote is $K_c^2$ (D02).""")
remind("C02")
P("P383", "secant slope and tangent slope", r"""On a curve y(x) the secant from the origin to a point has slope y/x; the
tangent at that point has slope dy/dx (P19). For a straight line through 0 they agree; for a bending curve they do not.
Blood's "apparent viscosity" τ/γ̇ is the secant (trap T27).""", code=r"""
x = 4.0; y = lambda x: (0.07 + 0.053*np.sqrt(x))**2          # a Casson-like curve (lambda: a one-line function, Ch. 1 P29)
h = 1e-6                                                     # a small step for the centred difference
print(y(x)/x, (y(x + h) - y(x - h))/(2*h))                  # secant 0.00774 > tangent 0.00466
""")
note("N21 [B]", r"""Apparent viscosity $\mu_{app}=\tau/\dot\gamma$ — the secant (`ch16.apparent_viscosity`); the tangent
$d\tau/d\dot\gamma$ is `ch16.differential_viscosity`.""")
D("D02", r"The Casson apparent viscosity and its limits", ref="16.6")
note("N22 [B]", r"""Whitmore's fit «16.7»: $k_c$ = 1.53 (a pure number) and $\sqrt{\tau_y/\mu_p}$ = 2.0 s$^{-1/2}$, so
$\tau_y=4\mu_p$ (in Pa when μ_p is in Pa s and the 2.0 in s$^{-1/2}$) and $K_c^2=1.53^2\mu_p=2.34\mu_p$. Two ways to write the
constants (trap T14):

| Constant | Units | In |
|---|---|---|
| $k_c$ | none | «16.6» |
| $K_c=k_c\sqrt{\mu_p}$ | (Pa s)$^{1/2}$ | «16.207» |""")
code(r"""
wc = ch16.whitmore_casson(np.array([4.0, 400.0]), blood["mu_p"], k_c=casson["k_c"], s0=casson["s0"])   # (16.7) at γ̇ = 4 and 400 s⁻¹
ty, Kc = wc["tau_y"], wc["Kc"]                               # yield stress [Pa] and K_c [(Pa s)^½]
print(f"tau_y = {ty*1e3:.2f} mPa, Kc = {Kc:.5f} (Pa s)^½, Kc² = {Kc**2*1e3:.3f} mPa s")
print("mu_app [mPa s]:", np.round(wc["mu_app"]*1e3, 2), "| tangent [mPa s]:", np.round(wc["mu_tangent"]*1e3, 2))
""")
worked("apparent viscosity at two shear rates", r"""
With Whitmore's constants and μ_p = 1.2 mPa s:
1. √τ_y = √(4.8 × 10⁻³) = 0.0693; K_c = 1.53 × √0.0012 = 0.0530.
2. γ̇ = 4 s⁻¹: √τ = 0.0693 + 0.0530 × 2 = 0.1753, τ = 30.7 mPa, μ_app = 30.7/4 = **7.68 mPa s** (6.4 μ_p).
3. γ̇ = 400 s⁻¹: √τ = 0.0693 + 0.0530 × 20 = 1.129, τ = 1.275 Pa, μ_app = **3.19 mPa s** (2.66 μ_p).
4. γ̇ → ∞: μ_app → K_c² = **2.81 mPa s**.

Slow shear makes blood look 2.7 times as viscous as fast shear.""")
code(r"""
gd = np.logspace(-1, 3, 200)                                # shear rates 0.1 … 1000 s⁻¹
tau = ch16.casson_stress(gd, ty, Kc)                        # (16.207): (√τ_y + K_c √γ̇)² [Pa]
mu_app = ch16.casson_apparent_viscosity(gd, ty, Kc)         # secant τ/γ̇ = (√(τ_y/γ̇) + K_c)² [Pa s] (D02)
mu_tan = ch16.differential_viscosity(gd, ty, Kc)            # tangent dτ/dγ̇ [Pa s] (trap T27)
gd_back = ch16.casson_shear_rate(tau, ty, Kc)               # invert the law: γ̇ from τ
print("round trip max error:", np.max(np.abs(gd_back - gd)/gd))   # ≈ 1e-15
print("below yield:", ch16.casson_shear_rate(np.array([0.5*ty, 0.99*ty]), ty, Kc))   # zero shear rate under τ_y
""")
scratch(r"""
tau_mine = (np.sqrt(ty) + Kc*np.sqrt(gd))**2                 # (16.207) by hand
gd_mine = np.where(tau_mine >= ty, ((np.sqrt(tau_mine) - np.sqrt(ty))/Kc)**2, 0.0)   # inverse, zero below yield (np.where, P46)
assert np.allclose(tau_mine, ch16.casson_stress(gd, ty, Kc))     # same stress
assert np.allclose(gd_mine, ch16.casson_shear_rate(tau_mine, ty, Kc))   # same inverse
print("from-scratch Casson law agrees")
""", explain=r"""The law written by hand in one line, its inverse with `np.where` (zero shear rate below the yield stress), and two
asserts: the library computes exactly this.""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.8))
sg = np.linspace(0, 32, 200)                                   # √γ̇ axis [s^-1/2]
ax.plot(sg, np.sqrt(ty) + Kc*sg, color=COLORS["blue"], lw=2, label="Casson (16.207)")
ax.plot(sg, Kc*sg, "--", color=COLORS["muted"], label="Newtonian, μ = K_c²")
ax.plot(0, np.sqrt(ty), "o", color=COLORS["rose"]); ax.annotate("intercept √τ_y", (0, np.sqrt(ty)), xytext=(6, 0.1), fontsize=8,
       arrowprops=dict(arrowstyle="->", color=COLORS["rose"]))   # rose dot and arrow: the yield-stress intercept √τ_y
ax.set(xlabel="√γ̇ [s$^{-1/2}$]", ylabel="√τ [Pa$^{1/2}$]", title="A straight line in square roots", xlim=(0, 32))
ax.legend(fontsize=8)
bx.semilogx(gd, mu_app*1e3, color=COLORS["blue"], lw=2, label="μ_app = τ/γ̇ (secant)")
bx.axhline(Kc**2*1e3, ls="--", color=COLORS["muted"], label=f"asymptote K_c² = {Kc**2*1e3:.2f} mPa s")
bx.plot(4, wc["mu_app"][0]*1e3, "o", color=COLORS["rose"], label="secant at 4 s⁻¹")
bx.plot(4, wc["mu_tangent"][0]*1e3, "s", color=COLORS["accent"], label="tangent at 4 s⁻¹")
bx.axvspan(0.1, 50, color=COLORS["grid"], alpha=0.5); bx.text(1, 0.8, "rouleaux", fontsize=8, color=COLORS["muted"])
bx.axvspan(200, 1000, color=COLORS["teal"], alpha=0.08); bx.text(220, 0.8, "single cells", fontsize=8, color=COLORS["muted"])
bx.set(xlabel="shear rate γ̇ [s$^{-1}$]", ylabel="viscosity [mPa s]", ylim=(0, 15),
       title="Slowly sheared blood looks thickest (ours)")
bx.legend(fontsize=7, loc="upper right")
savefig(fig, "ch16", "c02_casson"); plt.show()
""", see=r"""Left: Casson's law as a straight blue line in √τ against √γ̇, cutting the axis at √τ_y; the Newtonian line of the
same high-shear viscosity passes through the origin (grey dashed). Right: the apparent viscosity (secant, blue) falling from
large values at low shear to the asymptote $K_c^2$; the rose dot and the purple square are the secant and tangent at
4 s⁻¹; the shaded bands mark (illustrative) the rouleaux and single-cell regimes.""",
    read=r"""The intercept is the yield stress: you need a finite √τ before any shear starts. At 4 s⁻¹ the secant is 7.7 mPa s
but the tangent only 4.6 mPa s — they are different numbers, and "apparent viscosity" means the secant.""",
    change=r"""…we doubled τ_y? Only the low-shear end of the right panel rises; the asymptote $K_c^2$ does not move, because
at high shear the yield term $\sqrt{\tau_y/\dot\gamma}$ vanishes.""")
plotly(r"""
ty_vals = list(np.round(np.linspace(0.0, 0.02, 25), 5))      # yield stresses 0 … 20 mPa (25 positions)
ty_vals[6] = round(ty, 5)                                     # include Whitmore's 4.8 mPa
def curves(t_y):                                              # three traces for one yield stress
    return {"μ_app (secant)": (gd, ch16.casson_apparent_viscosity(gd, t_y, Kc)*1e3),
            "K_c² (asymptote)": (gd, Kc**2*1e3 + 0*gd),
            "tangent dτ/dγ̇": (gd, ch16.differential_viscosity(gd, t_y, Kc)*1e3)}
fig = slider_figure(curves, "τ_y", ty_vals, unit="Pa", xlabel="shear rate γ̇ [1/s] (log)", ylabel="viscosity [mPa s]",
                    title="The yield stress lifts only the slow-shear end", xrange=[-1, 3], yrange=[0, 30], active=6)
fig.update_xaxes(type="log")                                  # log shear-rate axis (its range is in powers of ten: 0.1 … 1000)
fig.show()
for t_y in (ty, 0.02):                                        # Whitmore's τ_y and the slider's largest, 20 mPa
    over = ch16.casson_apparent_viscosity(np.array([200.0, 1000.0]), t_y, Kc)/Kc**2 - 1   # how far above K_c² at 200 and 1000 s⁻¹
    print(f"τ_y = {t_y*1e3:.1f} mPa: μ_app is {100*over[0]:.0f} % above K_c² at 200 s⁻¹, {100*over[1]:.0f} % at 1000 s⁻¹")
""", explain=r"""**Figure F1.** Slider: the yield stress τ_y from 0 to 20 mPa (opens on Whitmore's 4.8 mPa); every slider position
was computed in Python when the cell ran (P17), so the figure works on the web page. **What you see.** The secant (apparent)
and tangent viscosities against shear rate, with the asymptote $K_c^2$. **How to read it.** At τ_y = 0 the three lines are
one: a Newtonian fluid. **What would change if…** you slide τ_y up: the curves lift most at the left and approach $K_c^2$
only slowly at high shear — with Whitmore's τ_y the secant is still 19 % above $K_c^2$ at 200 s⁻¹ and 8 % above at
1000 s⁻¹, and with τ_y = 20 mPa 41 % and 18 % (printed under the figure), because the yield term
$\sqrt{\tau_y/\dot\gamma}$ fades only like $\dot\gamma^{-1/2}$.""")
note("N23 [C]", r"""Casson's law treats blood as a homogeneous fluid. In tubes below about a third of a millimetre the cells
leave a plasma layer at the wall — C05 models exactly that.""")
whatif(r"""…we pushed this blood through a pipe? The stress is zero on the axis and largest at the wall, so near the axis it is
below τ_y: the centre moves as a solid plug. C16 works that out; the explainer embedded there lets you watch it.""")

# =====================================================================================================================
# A.3 §16.3 — C03 … C16
# =====================================================================================================================
nb.section("16.3", "Modeling of Flow in Blood Vessels", intro="**What is this section about?** The long modelling section, "
           "read in four arcs: steady flow (C04, C05); storage and waves in elastic tubes (C06–C08); oscillating flow "
           "(C09–C11); geometry and odd materials (C12–C16). One number, α, decides which arc's tools a vessel needs.")
core("C03", r"One number for every vessel: the Womersley number «16.12»",
     "When can the pulsing flow in a vessel be treated as steady, and when not?")
problem(r"""
Shake a glass of water gently and the whole glass of water moves together; shake a bathtub and only the water near the
walls notices the walls. Viscosity spreads a disturbance only so far in one cycle — about $\sqrt{\nu/\omega}$. Compare that
distance with the vessel's radius: **that comparison is α**. A capillary is far thinner than this distance (α ≪ 1); the
aorta far wider (α ≫ 1).""")
idea(r"""
|←──────── a ────────→|      the vessel radius
|←δ_e→|                      the Stokes layer √(2ν/ω): how far viscosity reaches in one cycle

small α:  the layer fills the tube       large α:  a thin skin at the wall, an inertial core
""", r"""$\alpha=\sqrt2\,a/\delta_e$ (D03 below).""")
recap("R02", "Separation in an adverse pressure gradient", r"""A boundary layer pushed against a rising pressure can stop and
turn back (Ch. 9). **New here:** in a pulse the deceleration phase supplies that rising pressure every beat — reversal near
the wall in C09 is its oscillating cousin.""", where="Ch. 9 §9.7", core_id="C03")
recap_md("Stokes' second problem", "Ch. 8 C11", r"""«8.38»: a wall oscillating in its own plane drags a layer whose
e-folding depth is $\delta_e=\sqrt{2\nu/\omega}$ (`LAM.stokes_layer`).""")
D("D03", r"What the Womersley number compares", ref="16.12")
worked("α for an artery", r"""
1. Radius a = 2 mm, heart rate 66 per minute → ω = 2π × 66/60 = 6.91 rad/s.
2. ν = μ/ρ = 3.5 × 10⁻³/1060 = 3.30 × 10⁻⁶ m²/s.
3. ω/ν = 2.09 × 10⁶ m⁻², its square root 1447 m⁻¹.
4. α = 0.002 × 1447 = **2.9**.
5. δ_e = √(2ν/ω) = **0.98 mm**: the viscous layer reaches halfway to the axis — neither limit.""")
code(r"""
omega = 2*np.pi*heart["HR"]/60                              # angular frequency of the heartbeat [rad/s]
nu = blood["mu"]/blood["rho"]                               # kinematic viscosity of blood [m²/s]
alpha_art = ch16.womersley_number(artery["a"], omega, nu)   # (16.12): α = a √(ω/ν) for our 2-mm artery
d_e = LAM.stokes_layer(nu, omega)["delta_e"]                # Stokes-layer depth √(2ν/ω) [m] (Ch. 8)
print(f"omega = {omega:.3f} rad/s, nu = {nu:.4e} m²/s, alpha = {alpha_art:.4f}, delta_e = {d_e*1e3:.3f} mm")
print("sqrt(2) a/delta_e =", np.sqrt(2)*artery["a"]/d_e)    # D03: the same α another way
print("from Re, St:", ch16.womersley_from_re_st(500.0, 0.2), "| regime:", ch16.womersley_regime(alpha_art))   # √(Re·St/4) = 5
print("our aorta (a = 12 mm):", round(ch16.womersley_number(aorta["a0"], omega, nu), 1))
tab = ch16.vessel_scales_table(HR=66.0)                     # cited typical vessel sizes (ours): Re, α, h/a per class
display(tab.round(4))                                       # the table, rounded to 4 decimals
""", explain=r"""
1. `womersley_number` evaluates «16.12» for the artery: 2.89. `LAM.stokes_layer` (Ch. 8) gives the depth 0.98 mm, and
   $\sqrt2\,a/\delta_e$ gives the same α — D03 in numbers.
2. `womersley_from_re_st` uses the book's diameter-based Re and St: with Re = 500 and St = 0.2, α = √(100/4) = 5.
3. `vessel_scales_table` lists seven vessel classes with cited typical sizes and speeds (labelled cited/ours) and their Re
   and α at 66 beats per minute. Our 12-mm aorta gives 17.4; the table's 11-mm aorta 15.9.""")
note("N24 [B]", r"""(From §16.2, out of book order.) Vessel walls have three layers (intima, media, adventitia); sizes run from
centimetre aortas to micrometre capillaries; veins hold most of the blood. Our scales table uses cited typical sizes
(labelled cited/ours) and gives Re, α and h/a for each class — the table the figure below plots.""")
figure(r"""
fig, ax = plt.subplots(figsize=(8, 3.8))
xs = np.arange(len(tab))                                      # one position per vessel class
ax.semilogy(xs, tab["alpha"], "o-", color=COLORS["blue"], label="Womersley number α")
ax.semilogy(xs, tab["Re"], "s--", color=COLORS["orange"], label="Reynolds number Re (diameter)")
ax.axhspan(1e-3, 1, color=COLORS["teal"], alpha=0.07); ax.axhspan(10, 1e4, color=COLORS["rose"], alpha=0.05)
ax.text(len(tab) - 0.6, 0.15, "quasi-steady α < 1", ha="right", fontsize=8, color=COLORS["teal"])
ax.text(len(tab) - 0.6, 40, "inertial α > 10", ha="right", fontsize=8, color=COLORS["rose"])
ax.set_xticks(xs, tab["vessel"], rotation=20, fontsize=8)
ax.set(ylabel="dimensionless number [–]", ylim=(1e-3, 1e4), title="Big vessels are inertial, small ones quasi-steady (cited sizes)")
ax.legend(fontsize=8, loc="upper center")
savefig(fig, "ch16", "c03_alpha_re"); plt.show()
cap, ao = tab.set_index("vessel").loc["capillary"], tab.set_index("vessel").loc["aorta"]   # the two ends of the tree
print(f"aorta → capillary: Re falls {np.log10(ao['Re']/cap['Re']):.1f} decades, α falls {np.log10(ao['alpha']/cap['alpha']):.1f} decades")
""", see=r"""For seven vessel classes (cited typical sizes), the Womersley number (blue) and the diameter Reynolds number
(orange) on a log axis, with the quasi-steady band α < 1 and the inertial band α > 10 shaded.""",
    read=r"""From the aorta to the capillaries the Reynolds number falls by about six decades (1666 → 0.0012) and α by about
3.4 decades (15.9 → 0.0058) — printed under the figure; Re falls faster because it carries the speed as well as the size
(Re ∝ ua, α ∝ a). The arterioles, capillaries and
venules sit deep in the quasi-steady band (Poiseuille at every instant); the aorta and venae cavae are inertial; medium
arteries and veins sit in between, where neither limit is good.""",
    change=r"""…the heart rate doubled? Every α grows by √2 = 1.41 (α ∝ √ω): the bands do not move, the blue dots do; Re does
not depend on the heart rate at all.""")
confusion(r"""(trap T6) **Which Re and St?**

| Source | Definitions | Link to α |
|---|---|---|
| this chapter | $Re=2au/\nu$, $St=2a\omega/u$ (diameter) | $\alpha^2=Re\,St/4$ |
| other sources | Re and St on other length scales | $\alpha^2=2\pi\,Re\,St$ |
| Exercise 16.2 | $St=a\omega/U$ | yet another factor |""")
note("N25 [C]", r"""Three ideas named here are taught later: compliance (C07, N45), the pulse as a wave faster than the flow
(C08) and reflections at branches (C13).""")
whatif(r"""…α were large in every vessel? Steady Poiseuille flow would never apply. It applies in the small vessels, which is where
the body's resistance sits — C04.""")

# ---------------------------------------------------------------------------------------------------------------------
# C04 — resistance
# ---------------------------------------------------------------------------------------------------------------------
core("C04", r"Resistance grows as $a^{-4}$: «16.20»", "Why do tiny changes in the radius of small arteries change blood pressure so much?")
problem(r"""
Pinch a garden hose to half its width and the flow drops far more than by half. Poiseuille's law says by how much: the flow
through a tube grows as the **fourth power** of its radius. Write it as "pressure drop = resistance × flow" and the
resistance carries $a^{-4}$. The body's arterioles use this lever.""")
idea(r"""
pressure drop = resistance × flow rate        Δp = R Q
resistance of a tube ∝ 1 / radius⁴
radius × 0.8  ──►  R × 0.8⁻⁴ ≈ 2.4  ──►  2.4 × the pressure drop for the same flow
radius × 0.5  ──►  R × 16
""", r"""**A tube's resistance is a steep power of its radius**: the tube's version of "voltage = resistance × current", with
the resistance «16.20». A 20 % narrowing of the arterioles more than doubles the pressure needed to push the same flow.""")
recap("R03", "When Poiseuille's law holds", r"""Steady, fully developed, laminar flow of a Newtonian fluid in a straight rigid tube
(Ch. 8), with $Re=\rho\bar ud/\mu$ (`SIM.reynolds_number`). **New here:** a checklist of how wrong each assumption is in an
artery and which block relaxes it:

| Assumption | In an artery | Relaxed in |
|---|---|---|
| steady | pulsatile | C09–C11 |
| rigid | elastic | C06–C08 |
| straight | curved | C14 |
| unbranched | branches every few cm | C12, C13 |
| Newtonian | shear-thinning, yield stress | C02, C16 |
| fully developed | entrance lengths comparable with segments | R08, N27 |""", where="Ch. 8 §8.2 (C03); Ch. 4 (C15)", core_id="C04")
recap("R04", "The Poiseuille profile", r"""«16.13» — the same as «8.6». **New here:** written with the overall drop.""",
      where="Ch. 8 §8.2", core_id="C04")
recap("R05", "Gradient from the overall drop", r"""«16.14», so «16.15». $G=-dp/dx=\Delta p/L>0$ (trap T4).""", where="Ch. 8",
      core_id="C04")
recap("R06", "Centreline speed, flux and mean speed", r"""«16.16», «16.17», «16.18». The integral is D04 step 2.""",
      where="Ch. 8 (D06)", core_id="C04")
recap("R07", "Wall shear stress", r"""«16.19» (slip S1: the book's last member has the wrong sign): this chapter's τ is a
positive magnitude (trap T5), the Ch. 8 code uses signed τ, «8.8».""", where="Ch. 8 (8.8)", core_id="C04")
code(r"""
r = np.linspace(0, artery["a"], 5)                          # five radii across the 2-mm artery [m]
u = ch16.poiseuille_velocity(r, artery["a"], 100.0, 0.1, blood["mu"])   # (16.15) with Δp = 100 Pa over L = 0.1 m [m/s]
u8 = LAM.pipe_poiseuille(r, artery["a"], mu=blood["mu"], G=100.0/0.1)   # Ch. 8's (8.6) with G = Δp/L
Q = ch16.poiseuille_flow_rate(artery["a"], 100.0, 0.1, blood["mu"])     # (16.17) [m³/s]
tau_w = ch16.poiseuille_wall_stress(artery["a"], 100.0, 0.1)            # (16.19), a positive magnitude [Pa]
print(np.allclose(u, u8), f"Q = {Q:.4e} m³/s = Ch. 8: {LAM.pipe_flow_rate(artery['a'], mu=blood['mu'], G=1000.0)[0]:.4e}; tau_w = {tau_w:.3f} Pa")
print("Re of the mean flow:", SIM.reynolds_number(Q/(np.pi*artery["a"]**2), 2*artery["a"], nu=nu))   # ρūd/μ
""", explain=r"""
1. `poiseuille_velocity` is «16.15» and agrees with Ch. 8's `pipe_poiseuille` given $G=\Delta p/L$ (trap T4: G is the
   positive drop per length).
2. `poiseuille_flow_rate` and Ch. 8's `pipe_flow_rate` return the same Q; the wall stress is +1.0 Pa, a magnitude (trap T5).""")
recap("R08", "Entrance length", r"""A flat inlet profile needs a length to become parabolic; wall shear is higher there (Ch. 8
sketch, Ch. 9 boundary-layer growth). **New here:** some vessels are too short ever to be fully developed.""",
      where="Ch. 8 §8.2; Ch. 9 (C02)", core_id="C04")
remind("C04")
D("D04", r"Vascular resistance from Poiseuille's flux", ref="16.20")
worked("halve the radius", r"""
1. μ = 4 mPa s, L = 0.1 m, a = 1 mm: R = 8 × 0.004 × 0.1/(π × 10⁻¹²) = **1.02 × 10⁹ Pa s m⁻³**.
2. a = 0.5 mm: R × 2⁴ = **1.63 × 10¹⁰**.
3. Ten thousand such 0.5 mm tubes in parallel: R/10⁴ = 1.63 × 10⁶ — parallel paths divide the resistance (P381).""")
code(r"""
R1 = ch16.vascular_resistance(4e-3, 0.1, 1e-3)              # (16.20): a = 1 mm [Pa s/m³]
R2 = ch16.vascular_resistance(4e-3, 0.1, 0.5e-3)            # half the radius: 16 times the resistance
Rpar = ch16.network_resistance([R2]*10000, arrangement="parallel")   # 10 000 such tubes side by side
R_art = ch16.vascular_resistance(blood["mu"], arteriole["L"], arteriole["a"])   # one 25-µm, 2-mm arteriole
print(f"{R1:.4g}  {R2:.4g}  ratio {R2/R1:.1f}  parallel {Rpar:.4g}  arteriole {R_art:.3g} Pa s/m³")
""")
scratch(r"""
a, mu_, L_, dp_ = 1e-3, 4e-3, 0.1, 100.0                     # tube radius [m], viscosity [Pa s], length [m], drop [Pa]
r = np.linspace(0, a, 2001)                                  # radii from the axis to the wall
u = a**2*dp_/(4*mu_*L_)*(1 - r**2/a**2)                      # (16.15) by hand [m/s]
Q_mine = np.trapezoid(2*np.pi*r*u, r)                        # flux = Σ u × ring area 2πr dr (D04 step 1)
assert np.isclose(dp_/Q_mine, ch16.vascular_resistance(mu_, L_, a), rtol=1e-6)   # R = Δp/Q matches (16.20)
assert np.isclose(ch16.vascular_resistance(mu_, L_, a/2)/ch16.vascular_resistance(mu_, L_, a), 16.0)   # the a⁻⁴ lever
print("R from the ring sum:", f"{dp_/Q_mine:.6g}")
""", explain=r"""The rings are added with `np.trapezoid` (the trapezoid rule, a sum of thin strips — P27): the resulting R = Δp/Q equals
«16.20» to 10⁻⁶, and halving the radius multiplies R by exactly 16.""")
figure(r"""
a_ = np.logspace(np.log10(5e-6), np.log10(15e-3), 200)       # radii 5 µm … 15 mm
Rr = ch16.vascular_resistance(blood["mu"], 0.01, a_)          # resistance of a 1-cm segment [Pa s/m³]
fig, ax = plt.subplots(figsize=(7, 3.8))
ax.loglog(a_*1e3, Rr, color=COLORS["accent"], lw=2, label="R = 8μL/(πa⁴), L = 1 cm")
for name, row in tab.set_index("vessel").iterrows():          # vessel classes of the scales table
    Rv = ch16.vascular_resistance(blood["mu"], 0.01, row["a"])
    ax.plot(row["a"]*1e3, Rv, "o", color=COLORS["teal"]); ax.annotate(name, (row["a"]*1e3, Rv), fontsize=7, xytext=(4, 2),   # a teal dot and its vessel name beside it (no arrow)
                                                                     textcoords="offset points")   # a teal dot and its vessel name
ax.annotate("slope −4", (0.05, ch16.vascular_resistance(blood["mu"], 0.01, 5e-5)), xytext=(0.3, 1e17), fontsize=9,
            arrowprops=dict(arrowstyle="->"))
ax.set(xlabel="radius a [mm]", ylabel="resistance of 1 cm [Pa s/m³]", title="Halve the radius, sixteen times the resistance")
ax.legend(fontsize=8)
savefig(fig, "ch16", "c04_resistance"); plt.show()
v_ = tab.set_index("vessel")                                   # the scales table indexed by vessel name
print(f"capillary / aorta resistance per length: (a_aorta/a_cap)⁴ = {(v_.loc['aorta', 'a']/v_.loc['capillary', 'a'])**4:.2e}")
""", see=r"""The resistance of a 1-cm segment against its radius on log–log axes (a straight line of slope −4), with the
vessel classes of the scales table as teal dots.""",
    read=r"""Each decade of radius is four decades of resistance: a 4-µm capillary segment resists about 10¹⁴ times
(6 × 10¹³, printed under the figure) more than an 11-mm aortic one of the same length. On log–log axes a power law is a
straight line whose slope is the exponent (P13).""",
    change=r"""…an arteriole narrowed by 19 %? $1/0.81^4=2.3$ — its resistance would roughly double; that is how a small
change of tone controls the pressure.""")
note("N27 [B]", r"""Entrance length «16.37» (laminar for Re > 50). **Which correlation?** (slip S31)

| Source | laminar | turbulent |
|---|---|---|
| the book | $\ell/d=0.06\,Re$ | $\ell/d=0.693\,Re^{1/4}$ |
| Çengel (cited via Wikipedia) | $\ell/d\approx0.0575\,Re$ | $\ell/d=1.359\,Re^{1/4}$ |""")
code(r"""
print(ch16.entrance_length(500.0, 4e-3, regime="laminar"), "m")   # (16.37): 0.06 × 500 × 4 mm = 0.12 m
""", explain=r"""**What does this show?** At Re = 500 in a 4-mm artery the flow needs 12 cm to become parabolic — longer than many
arterial segments between branches, so Poiseuille's profile is an idealisation there.""")
note("N26 [C]", r"""Poiseuille's law fits muscular arteries best; heat transfer in perfused tissue is named in the book and not
treated here. **N28 [C]** In pulsatile flow the entrance length depends on Re and α; the wall shear settles sooner than the
centreline speed (α is C03, the oscillating profile C09).""")
whatif(r"""…the blood were not one fluid? In tubes below about 0.3 mm the cells crowd to the centre and leave plasma at the wall —
then a single μ no longer describes the flow: C05.""")

# ---------------------------------------------------------------------------------------------------------------------
# C05 — two-layer model
# ---------------------------------------------------------------------------------------------------------------------
core("C05", r"A plasma layer lubricates narrow tubes: «16.33»", "Why does blood seem less viscous in narrower tubes?")
problem(r"""
Measure blood's viscosity in glass tubes of shrinking diameter and, below about a third of a millimetre, it seems to get
*thinner*. The red cells drift away from the wall, leaving a thin layer of plain plasma. **The thick core slides on that
runny skin** — and the narrower the tube, the larger the skin's share.""")
idea(r"""
 wall ─────────────────────────────  u = 0
       plasma layer, μ_p, thickness δ    (runny skin: shears a lot)
 ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─   r = a − δ: speed and stress continuous
       core, μ_c > μ_p                   (thick: nearly a plug)
 axis ─────────────────────────────
""")
note("N29 [B]", r"""(From §16.2.) The **Fåhræus–Lindqvist effect**: μ_app falls with diameter over a range of small tubes,
rising again in the very smallest. The **Fåhræus effect**: the tube hematocrit $H_T$ is below the discharge hematocrit
$H_D$ (= the feed $H_F$). The curve below uses the cited in-vitro law of Pries et al. via Secomb & Pries (2013),
`ch16.pries_relative_viscosity` (labelled cited).""")
note("N30 [B]", r"""(From §16.2.) Red cells and volume are conserved: «16.8». Eliminating $Q_F$ gives $H_c/H_F=1+Q_p/Q_c$
(result given): the core is richer in cells than the feed, and the tube holds fewer cells than flow through it, because the
cell-rich core moves faster than the mean.""")
code(r"""
cc = ch16.cell_conservation(1.0e-12, 0.2e-12, 0.5, 20e-6, 2e-6)   # (16.8): Q_c, Q_p [m³/s], H_c, a, δ [m]
print({k: float(f"{v:.4g}") for k, v in cc.items()})               # H_F = 0.417 < H_c = 0.5; H_T/H_D < 1 (Fåhræus effect)
""")
note("N31 [B]", r"""The model: a core of viscosity $\mu_c$ inside a plasma layer of viscosity $\mu_p$ and thickness δ; stress
continuous at the interface; regular on the axis: «16.21», «16.22», «16.23» — the printed «16.21» has the wrong sign on the
right (slip S2). **N32 [B]** «16.24», «16.25», «16.26».""")
remind("C05")
D("D05", r"The two-layer profiles", ref="16.27", notes=r"**N33 [B]** the two profiles «16.27» and «16.28».")
D("D06", r"The layer fluxes", ref="16.31", notes=r"**N34 [B]** the plasma flux «16.29», the core flux «16.30» and the total «16.31».")
slip("S3", r"«16.30p»", r"«16.30»", r"""Only the corrected core flux adds to the plasma flux «16.29» to give the total «16.31».
The printed core flux is too large: printed minus corrected is
$\dfrac{\pi\Delta p}{4\mu_pL}(a^4-a^2b^2)=+\dfrac{\pi\Delta p\,a^2\delta(2a-\delta)}{4\mu_pL}>0$ (with $b=a-\delta$), so the
printed form overstates the core's flux — by 75 % in our 40-µm tube (the cell below).""")
code(r"""
args = (micro["a"], micro["delta"], 1.0, 1e-3, blood["mu_p"], micro["mu_c"])   # a, δ [m], Δp [Pa], L [m], μ_p, μ_c [Pa s]
good = ch16.two_layer_flow_rates(*args)                     # (16.29), (16.30, corrected), (16.31)
bad = ch16.two_layer_flow_rates(*args, printed_16_30=True)  # the printed (16.30)
print(f"Qp + Qc - Q: corrected {good['Qp'] + good['Qc'] - good['Q']:.2e}, printed {bad['Qp'] + bad['Qc'] - good['Q']:.2e} m³/s")
print(f"relative error of the printed core flux: {bad['Qc']/good['Qc'] - 1:+.3f}")
""")
D("D07", r"The apparent viscosity and its thin-layer limit", ref="16.33")
worked("a 40 µm tube", r"""
Relative units μ_p = 1, μ_c = 4, δ/a = 0.1:
1. (1 − 0.1)⁴ = 0.6561.
2. 1 − μ_p/μ_c = 0.75.
3. Bracket 1 − 0.6561 × 0.75 = 0.508.
4. μ_app = 1/0.508 = **1.97** μ_p — half the core's 4.
5. The thin-layer form «16.34»: 4/(1 + 0.4 × 3) = 1.82 (8 % low).

With our SI inputs (a = 20 µm, δ = 2 µm, μ_p = 1.2, μ_c = 4.0 mPa s): **2.22 mPa s**; in a 500 µm tube with the same δ,
3.73 mPa s.""")
code(r"""
r = np.linspace(0, micro["a"], 6)                           # six radii from the axis to the wall [m]
u2 = ch16.two_layer_velocity(r, micro["a"], micro["delta"], 1.0, 1e-3, blood["mu_p"], micro["mu_c"])   # (16.27)/(16.28) [m/s]
mu_two = ch16.apparent_viscosity_two_layer(micro["a"], micro["delta"], blood["mu_p"], micro["mu_c"])   # (16.33)
mu_thin = ch16.apparent_viscosity_thin_layer(micro["a"], micro["delta"], blood["mu_p"], micro["mu_c"]) # (16.34)
mu_500 = ch16.apparent_viscosity_two_layer(250e-6, micro["delta"], blood["mu_p"], micro["mu_c"])       # a 500-µm tube
print("u [mm/s]:", np.round(u2*1e3, 4))
print(f"mu_app = {mu_two*1e3:.3f} mPa s (thin-layer {mu_thin*1e3:.3f}); in a 500-µm tube {mu_500*1e3:.2f} mPa s")
""")
scratch(r"""
a, d = micro["a"], micro["delta"]                             # radius, layer [m]
Q_quad = quad(lambda r: 2*np.pi*r*ch16.two_layer_velocity(r, a, d, 1.0, 1e-3, blood["mu_p"], micro["mu_c"]), 0, a,   # ∫ 2πr u dr from axis to wall [m³/s] …
              points=[a - d])[0]                              # … told about the kink at the interface r = a − δ; [0] keeps the value
assert np.isclose(Q_quad, good["Q"], rtol=1e-10)              # the closed-form total (16.31)
mu_hand = np.pi*a**4*1.0/(8*1e-3*Q_quad)                      # (16.32) solved for μ_app
assert np.isclose(mu_hand, mu_two, rtol=1e-10)                # equals (16.33)
print(f"quad flux {Q_quad:.6e} m³/s; mu_app by hand {mu_hand*1e3:.4f} mPa s")
""", explain=r"""`quad` adds the profile's rings numerically (told where the kink is); the flux equals «16.31» to 10⁻¹⁰, and the
apparent viscosity solved from «16.32» by hand equals «16.33».""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.8))
rr = np.linspace(-micro["a"], micro["a"], 401)                # across the whole diameter [m]
ux = ch16.two_layer_velocity(np.abs(rr), micro["a"], micro["delta"], 1.0, 1e-3, blood["mu_p"], micro["mu_c"])
ug = ch16.poiseuille_velocity(np.abs(rr), micro["a"], 1.0, 1e-3, micro["mu_c"])   # one fluid of the core's viscosity
ax.plot(rr*1e6, ux*1e3, color=COLORS["blue"], lw=2, label="two-layer (model)")
ax.plot(rr*1e6, ug*1e3, "--", color=COLORS["muted"], label="Poiseuille with μ_c")
for s in (-1, 1):                                             # the plasma layers at both walls
    ax.axvspan(s*(micro["a"] - micro["delta"])*1e6, s*micro["a"]*1e6, color=COLORS["accent"], alpha=0.15)
ax.text(0, 0.02, "core", ha="center", fontsize=8); ax.text(micro["a"]*1e6*0.95, 0.02, "plasma", ha="center", fontsize=7, rotation=90)
ax.set(xlabel="r [µm]", ylabel="u [mm/s]", title="The core slides on a plasma skin (Δp/L = 1 kPa/m)")
ax.legend(fontsize=8, loc="upper left")
D_um = np.logspace(1, 3, 40)                                  # diameters 10 … 1000 µm
def fl_rel(D, delta=micro["delta"]):                          # our closure's μ_app/μ_p (NaN if out of range)
    try:                                                      # the closure can fail when δ is too large for the tube …
        return ch16.fahraeus_lindqvist_model(D*0.5e-6, delta, 0.45, blood["mu_p"], blood["T"])["mu_app"]/blood["mu_p"]
    except ValueError:                                        # … then no point is drawn (NaN leaves a gap in the curve)
        return np.nan
bx.semilogx(D_um, [fl_rel(D) for D in D_um], color=COLORS["blue"], lw=2, label="two-layer, δ = 2 µm (model, ours)")
bx.semilogx(D_um, ch16.pries_relative_viscosity(D_um, 0.45), "--", color=COLORS["muted"], label="Pries in vitro (cited)")
bx.set(xlabel="tube diameter D [µm]", ylabel="μ_app / μ_p [–]", title="Narrow tubes look less viscous (H = 0.45)")
bx.legend(fontsize=8)
savefig(fig, "ch16", "c05_two_layer"); plt.show()
""", see=r"""Left: the two-layer profile across a 40-µm tube (blue), with the plasma layers shaded purple, and the Poiseuille
profile of a single fluid with the core's viscosity (grey dashed). Right: the relative apparent viscosity against diameter
for our model closure (blue, δ = 2 µm) and the cited in-vitro law (grey dashed).""",
    read=r"""The blue profile is flatter in the middle and steeper in the layer: speeds match at the interface but slopes do
not (the kink). On the right both curves fall as the tube narrows: the model gets the direction and rough size of the
Fåhræus–Lindqvist effect.""",
    change=r"""…we went below about 10 µm? The cited curve rises again (cells deform and plug the tube); our model has no
deforming cells and cannot show that.""")
note("N35 [B]", r"""The core's viscosity from its hematocrit: «16.35» with «16.36», T in kelvin, valid to $H_c$ = 0.6. $\alpha_c$
is not the Womersley number (trap T2).""")
code(r"""
print(f"{ch16.charm_kurland_core_viscosity(0.5, 310.15, 1.2e-3)*1e3:.3f} mPa s")   # (16.35)-(16.36) at H_c = 0.5, 37 °C
""")
note("N36 [B]", r"""Closing the model (ours): given a, δ, $H_F$, $\mu_p$, T, solve $H_c=H_F[1+Q_p/Q_c]$ with $\mu_c(H_c)$ from
«16.35» — a monotone root find (`brentq`, P108).""")
code(r"""
fl = ch16.fahraeus_lindqvist_model(micro["a"], micro["delta"], 0.45, blood["mu_p"], blood["T"])   # closure (model, ours)
print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in fl.items()})   # H_c, μ_c, μ_app, H_T, residual
""")
confusion(r"""(trap T19) The curve depends on our choice of δ, which the book leaves free — it is a model, not a prediction.""")
plotly(r"""
d_vals = list(np.round(np.linspace(0.5e-6, 5e-6, 19), 8))     # plasma-layer thickness 0.5 … 5 µm (19 positions)
D_f = np.logspace(1.2, 3, 30) if FAST else np.logspace(1.2, 3, 40)   # diameters [µm]
pries = ch16.pries_relative_viscosity(D_f, 0.45)              # the cited curve does not depend on δ
def fl_curve(delta):                                          # model curve for one layer thickness
    return {"two-layer model (ours)": (D_f, [fl_rel(D, delta) for D in D_f]), "Pries in vitro (cited)": (D_f, pries)}
fig = slider_figure(fl_curve, "δ", d_vals, unit="m", xlabel="tube diameter D [µm] (log)", ylabel="μ_app/μ_p [–]",
                    title="A thicker plasma layer lowers the apparent viscosity", xrange=[1.2, 3], yrange=[1, 3.6],
                    active=d_vals.index(2e-6))
fig.update_xaxes(type="log"); fig.show()
""", explain=r"""**Figure F2.** Slider: the plasma-layer thickness δ from 0.5 to 5 µm (opens on our 2 µm). **What you see.** Our model's
relative apparent viscosity against diameter beside the fixed cited curve. **How to read it.** A thicker layer pulls the
model curve down at every diameter, most in the narrow tubes. **What would change if…** δ shrank to zero: the model would
give the core viscosity at every diameter — no effect at all.""")
whatif(r"""…the wall were not rigid? A pressure inside an elastic tube stretches it, the radius grows and, with R ∝ a⁻⁴, the
resistance falls — C06.""")

# ---------------------------------------------------------------------------------------------------------------------
# C06 — tube law
# ---------------------------------------------------------------------------------------------------------------------
core("C06", r"An artery is a spring: the tube law «16.41»", "How much does an artery widen when the pressure inside rises?")
problem(r"""
Blow into a long balloon: the higher the pressure, the wider it gets, and the rubber is stretched around its circumference.
Arteries are like that, with a stiffer rubber. **The stretch around the circumference — the hoop stress — balances the
pressure**; Hooke's law turns stretch into a wider radius.""")
md(r"""**Three symbols for the rest of the chapter.** $p_e$ is the pressure **outside** the vessel (the surrounding tissue, or the
chamber of a lab rig; we take $p_e=0$ in most examples). The difference $p-p_e$, inside minus outside, is the **transmural
pressure** ("across the wall"): it is what stretches the wall — raise the inside and the outside by the same amount and
nothing happens to the vessel. $a_0$ is the vessel's radius when $p=p_e$, the unstretched radius from which every bulge is
measured.""")
idea(r"""
   p_e (outside)          wall: thickness h, Young's modulus E
  ───────────────────────────────────────────
   p (inside)  ──►  p − p_e > 0 pushes the wall out:  a₀ → a > a₀
  ───────────────────────────────────────────
""", r"""**The wall is a spring loaded by the transmural pressure**: the higher $p-p_e$, the larger the radius.""")
P("P379", "Young's modulus and Poisson's ratio", r"""A solid bar stretched by a stress σ (force per area) stretches by the
strain ε = σ/E; E is Young's modulus (Pa) — rubber ~1 MPa, an artery wall ~0.1–1 MPa, steel 200 GPa. Stretching one way
thins the material the other way: Poisson's ratio ν̂ (0.5 for an incompressible wall, used in C11).""", code=r"""
E = 0.4e6; sigma = 1.0e5                       # artery wall [Pa], a 100 kPa stress
print(sigma/E)                                 # strain 0.25: a 25 % stretch
""")
P("P380", "a stretched curved membrane", r"""A curved membrane under tension T (force per length) pushes inward with a pressure
T/a, a the radius of curvature — Laplace's law, the same law as the surface-tension jump of Ch. 1 (Δp = σ/R for a cylinder).
For a wall of thickness h, T = σh.""", code=r"""
sigma, h, a = 1.0e5, 1.5e-3, 12e-3             # hoop stress [Pa], wall [m], radius [m]
print(sigma*h/a)                               # 12500 Pa ≈ 94 mmHg held
""")
D("D08", r"The hoop stress and the tube law", ref="16.41", notes=r"**N37 [B]** the half-tube balance «16.38» and the hoop stress «16.39» (steps 1–3).")
worked("an aorta at 100 mmHg", r"""
a₀ = 12 mm, h = 1.5 mm, E = 0.4 MPa:
1. Eh/a₀ = 0.4 × 10⁶ × 1.5/12 = 50 kPa (375 mmHg).
2. p − p_e = 100 mmHg = 13.3 kPa; a₀(p − p_e)/Eh = 13.3/50 = 0.267.
3. a = 12/(1 − 0.267) = **16.4 mm** (the linear estimate 12 × 1.267 = 15.2 mm).
4. Hoop stress (p − p_e)a/h = 13.3 × 16.4/1.5 = **145 kPa**.""")
code(r"""
a100 = ch16.tube_radius_hookean(100*MMHG, 0.0, aorta["a0"], aorta["E"], aorta["h"])   # (16.41) at 100 mmHg [m]
s100 = ch16.hoop_stress(100*MMHG, 0.0, a100, aorta["h"])    # (16.39) with the deformed radius [Pa]
print(f"a = {a100*1e3:.3f} mm, hoop stress = {s100/1e3:.1f} kPa")
try:                                                        # beyond Eh/a0 = 375 mmHg the law has no solution …
    ch16.tube_radius_hookean(400*MMHG, 0.0, aorta["a0"], aorta["E"], aorta["h"])   # 400 mmHg: past the singular pressure
except ValueError as e:                                     # … and the function says so (trap T16)
    print("refused:", e)
""", explain=r"""
1. `tube_radius_hookean` is «16.41»: 16.36 mm at 100 mmHg; `hoop_stress` uses that deformed radius.
2. The `try` / `except ValueError` pattern of C02 again: at 400 mmHg the bracket $1-a_0(p-p_e)/Eh$ is negative, the law has
   no solution, and the function raises an error that the `except` block catches and prints (trap T16).""")
figure(r"""
pp = np.linspace(0, 350, 300)*MMHG                            # transmural pressures 0 … 350 mmHg [Pa]
aa = ch16.tube_radius_hookean(pp, 0.0, aorta["a0"], aorta["E"], aorta["h"])   # (16.41) [m]
lin = aorta["a0"]*(1 + aorta["a0"]*pp/(aorta["E"]*aorta["h"]))  # linear (small-pressure) form [m]
fig, ax = plt.subplots(figsize=(7.5, 3.8)); bx = ax.twinx()
ax.plot(pp/MMHG, aa/aorta["a0"], color=COLORS["teal"], lw=2, label="tube law (16.41)")
ax.plot(pp/MMHG, lin/aorta["a0"], "--", color=COLORS["muted"], label="linear form")
ps = aorta["E"]*aorta["h"]/aorta["a0"]/MMHG                   # singular pressure Eh/a0 [mmHg]
ax.axvline(ps, color=COLORS["rose"], ls=":"); ax.text(ps - 3, 1.2, f"law breaks at\nEh/a₀ = {ps:.0f} mmHg", ha="right", fontsize=8, color=COLORS["rose"])
bx.plot(pp/MMHG, ch16.hoop_stress(pp, 0.0, aa, aorta["h"])/1e3, color=COLORS["rose"], lw=1.2, alpha=0.7, label="hoop stress")
ax.set(xlabel="transmural pressure p − p_e [mmHg]", ylabel="a / a₀ [–] (teal)", xlim=(0, 390), ylim=(1, 8),
       title="Stretch grows faster than pressure (our aorta)")
bx.set_ylabel("hoop stress [kPa] (rose)"); ax.legend(fontsize=8, loc="upper left")
savefig(fig, "ch16", "c06_tube_law"); plt.show()
""", see=r"""The radius ratio a/a₀ of our aorta against transmural pressure (teal), its linear form (grey dashed), and the hoop
stress (rose, right axis); the dotted line marks the pressure Eh/a₀ where the law breaks.""",
    read=r"""At small pressures the two radius curves agree; by 100 mmHg the full law is 8 % above the linear one, and it runs
away toward the singular pressure because the stress uses the deformed radius, which grows with the stress.""",
    change=r"""…the wall were thicker (larger h)? Eh/a₀ grows, the curves flatten and the singular pressure moves right; real
arteries avoid the singularity differently — collagen stiffens them at large strain.""")
note("N38 [B]", r"""In a slowly widening tube each slice obeys Poiseuille's law with its own radius: «16.42», «16.43» — the
lubrication idea of Ch. 8 (C05).""")
recap_md("Lubrication: slowly varying channels", "Ch. 8 C05", r"""When a channel changes width slowly compared with its
length, each cross-section carries the parallel-flow profile of its own width; only the pressure gradient varies.""")
code(r"""
x = np.linspace(0, 0.05, 6)                                 # positions along a 5-cm soft tube [m]
Qd = ch16.elastic_tube_flow_rate(20*MMHG, 10*MMHG, 0.0, 1e-3, 0.2e6, 0.1e-3, blood["mu"], 0.05)   # (16.46) [m³/s]
px = ch16.elastic_tube_pressure(x, Qd, 20*MMHG, 0.0, 1e-3, 0.2e6, 0.1e-3, blood["mu"])          # (16.44) integrated [Pa]
print("p(x) [mmHg]:", np.round(px/MMHG, 3))                 # falls from 20 to 10 mmHg, faster where the tube is narrower
""")
note("N39 [B]", r"""Put the tube law into «16.43» and separate: «16.44»; integrating gives «16.46» (result given).""")
slip("S4", r"«16.44p» and the same $\pi a_0$ in its integral «16.45p»", r"«16.44» and «16.45»",
     r"""A units check catches it: $8\mu Q\,dx/(\pi a_0^4)$ is a pressure, $8\mu Q\,dx/(\pi a_0)$ is a pressure times m³. Solving
the corrected «16.45» for Q gives «16.46», which the book prints correctly.""")
code(r"""
mu_q, Q_q, dx_q, a_q = UN.Q_(3.5e-3, "Pa*s"), UN.Q_(1e-6, "m**3/s"), UN.Q_(1e-3, "m"), UN.Q_(1e-3, "m")   # quantities with units
print("corrected is a pressure:", UN.check_dimensions(8*mu_q*Q_q*dx_q/(np.pi*a_q**4), "[pressure]"))   # True
print("printed is a pressure:  ", UN.check_dimensions(8*mu_q*Q_q*dx_q/(np.pi*a_q), "[pressure]"))      # False
Qr = ch16.poiseuille_flow_rate(1e-3, 10*MMHG, 0.05, blood["mu"])   # a rigid tube of the unstretched radius
print(f"distensible / rigid flux: {Qd/Qr:.4f}")               # 1.53: the stretched tube passes more
Q_stiff = ch16.elastic_tube_flow_rate(20*MMHG, 10*MMHG, 0.0, 1e-3, 1e12, 0.1e-3, blood["mu"], 0.05)   # E → ∞
print(f"rigid limit: {Q_stiff/Qr - 1:.1e}")                   # back to Poiseuille (16.17)
""")
plotly(r"""
E_vals = list(np.round(np.logspace(np.log10(0.05e6), np.log10(2e6), 20), -2))   # Young's modulus 0.05 … 2 MPa (log steps)
dps = np.linspace(0, 30, 61)*MMHG                             # pressure drops, outlet at p_e = 0 [Pa]
def flux(E):                                                  # distensible against rigid flux for one stiffness
    Qe = [ch16.elastic_tube_flow_rate(d, 0.0, 0.0, 1e-3, E, 0.1e-3, blood["mu"], 0.05) for d in dps]   # (16.46) for each drop [m³/s]
    Qr_ = ch16.poiseuille_flow_rate(1e-3, dps, 0.05, blood["mu"])   # (16.17), the rigid tube [m³/s]
    return {"distensible (16.46)": (dps/MMHG, np.array(Qe)*1e6), "rigid (16.17)": (dps/MMHG, Qr_*1e6)}
fig = slider_figure(flux, "E", E_vals, unit="Pa", xlabel="Δp = p₁ − p₂ [mmHg] (p₂ = p_e)", ylabel="Q [mL/s]",
                    title="Softer walls pass more flow at high pressure", yrange=[0, 6], active=9)
fig.show()
""", explain=r"""**Figure F3.** Slider: Young's modulus of a 1-mm tube with a 0.1-mm wall, 0.05 to 2 MPa. **What you see.** The flux of
the distensible tube «16.46» against the pressure drop, beside the rigid line «16.17». **How to read it.** The rigid flux is
a straight line; the distensible one curves upward because the inlet end inflates. **What would change if…** you slide E
up: the curve straightens onto the rigid line — a stiff tube is a rigid tube.""")
note("N40 [B]", r"""(From the end of §16.3.) The lungs' vessels are soft: their resistance falls as the pulmonary pressure
rises — «16.46» at work.""")
figure(r"""
pm = np.linspace(5, 25, 41)*MMHG                              # mean transmural pressure of a soft vessel [Pa]
Reff = [ch16.elastic_tube_resistance(p + 2*MMHG, p - 2*MMHG, 0.0, 1e-3, 0.2e6, 0.1e-3, blood["mu"], 0.05) for p in pm]   # (p₁ − p₂)/Q for a 4 mmHg drop [Pa s/m³]
fig, ax = plt.subplots(figsize=(6.5, 3.2))
ax.plot(pm/MMHG, np.array(Reff)/1e9, color=COLORS["accent"], lw=2, label="soft tube (16.46)")   # the purple curve: effective resistance
ax.axhline(ch16.vascular_resistance(blood["mu"], 0.05, 1e-3)/1e9, ls="--", color=COLORS["muted"], label="rigid tube")
ax.set(xlabel="mean transmural pressure [mmHg]", ylabel="(p₁ − p₂)/Q [10⁹ Pa s/m³]",
       title="A soft vessel resists less when inflated (ours, qualitative)"); ax.legend(fontsize=8)
savefig(fig, "ch16", "c06_soft_resistance"); plt.show()
""", see=r"""The effective resistance (p₁ − p₂)/Q of a soft 1-mm tube with a 4 mmHg drop against its mean pressure, and the
rigid tube's resistance (grey dashed).""",
    read=r"""As the mean pressure rises the tube inflates and, through R ∝ a⁻⁴, its resistance falls well below the rigid
value — the reason the low-pressure pulmonary circulation copes with large flow changes.""",
    change=r"""…E were doubled? The curve would flatten toward the rigid line: less inflation, less change.""")
note("N41 [C]", r"""The pulmonary circulation runs from the right ventricle through the lungs to the left atrium: larger,
softer, lower-resistance vessels (the loop sketch of N02). **N42 [C]** Its pressure pulses look like the systemic ones at a
lower pressure scale (C07).""")
whatif(r"""…the pressure changed in time? A stretchy vessel then stores blood when the pressure rises and gives it back when it
falls — the Windkessel, C07.""")

# ---------------------------------------------------------------------------------------------------------------------
# C07 — Windkessel
# ---------------------------------------------------------------------------------------------------------------------
core("C07", r"Store and drain: the Windkessel «16.49»",
     "The heart pushes for a third of each beat — why doesn't the pressure drop to zero in between?")
problem(r"""
Fill a balloon from a tap that you open for one second in every three, while the balloon leaks through a pinhole: the jet
from the pinhole never stops, it only weakens a little between fills. The heart pushes blood out for only about a third of
each beat, yet your arterial pressure swings between about 80 and 120 mmHg instead of falling to zero. **The stretchy aorta
is the balloon, and the many small vessels downstream are the pinhole.**""")
note("N43 [C]", r"""The opening picture: old fire engines had an air chamber (German *Windkessel*) between the pump and the hose;
the pump's strokes compressed the air, which kept water flowing between strokes. Otto Frank (1899) saw the aorta doing the
same.

```
heart ──► [aorta: stretchy store K] ──► [small vessels: resistor R] ──► veins (p ≈ 0)
```""")
idea("", r"""**A capacitor in parallel with a resistor** (P381): inflow = charging + leaking. During systole the inflow
charges the aorta; during diastole the stored blood leaks out through R, and the pressure decays slowly.""")
note("N44 [B]", r"""«16.47» in words: inflow into the aorta = rate of change of its volume + outflow into the periphery.
**N45 [B]** Compliance K = dv/dp; «16.48». For a Hookean segment (ours): $dV/dp=2\pi a^3L/(Eh)$.""")
code(r"""
K_aorta = ch16.compliance_from_tube(100*MMHG, 0.0, aorta["a0"], aorta["E"], aorta["h"], aorta["L"])   # dV/dp of our aorta at 100 mmHg [m³/Pa]
print(f"K = {K_aorta:.3e} m³/Pa = {K_aorta*1e6*MMHG:.2f} mL/mmHg (our Windkessel uses {windkessel['K']:.1e} m³/Pa)")
""")
remind("C07")
D("D09", r"The two-element Windkessel from a volume balance", ref="16.49")
P("P382", "a sum of decaying kicks", r"""For K dp/dt + p/R = Q(t) the solution adds up every past inflow Q(τ)dτ, each one
raising p by Q dτ/K and then fading as $e^{-(t-\tau)/RK}$ — a convolution. It extends the integrating factor of Ch. 9
(P210).""", code=r"""
RK, K = 1.74, 1.2e-8                           # time constant [s], compliance [m³/Pa]
kick = 1e-6                                    # 1 mL injected at t = 0
for t in (0, 1.74, 3.48):
    print(t, kick/K*np.exp(-t/RK))             # 83 Pa, 31 Pa, 11 Pa: the kick fades by e per RK
""")
D("D10", r"The Windkessel pressure, the diastolic decay and the periodic state", ref="16.52",
  notes=r"**N46 [B]** the solution «16.52» through «16.50» and «16.51» (steps 1–5).")
worked("diastolic decay", r"""
R = 1.45 × 10⁸ Pa s m⁻³, K = 1.2 × 10⁻⁸ m³/Pa:
1. RK = 1.74 s.
2. The beat lasts 60/66 = 0.909 s; systole 0.35 of it = 0.318 s; diastole 0.591 s.
3. Over diastole p falls by $e^{-0.591/1.74}$ = **0.712**.
4. From 105 mmHg at the valve's closing to 105 × 0.712 = **75 mmHg** before the next beat — the diastolic pressure.""")
code(r"""
T = 60/heart["HR"]; Ts = heart["systole_share"]*T           # beat and systole durations [s]
Qin = lambda t: ch16.ejection_half_sine(t, SV, T, Ts)       # half-sine ejection of one stroke volume per beat [m³/s]
t = np.linspace(0, 8*T, 4001)                               # eight beats
p_ode = ch16.windkessel_ode(t, Qin, windkessel["R"], windkessel["K"], 60*MMHG)   # (16.49) integrated numerically, from 60 mmHg [Pa]
tq = t[::40]                                                # a coarser grid for the quadrature form
p_q = ch16.windkessel_pressure(tq, Qin, windkessel["R"], windkessel["K"], 60*MMHG)   # (16.52) by quadrature [Pa]
print("ODE vs (16.52):", np.max(np.abs(p_q - p_ode[::40]))/MMHG, "mmHg")
s = ch16.windkessel_summary(SV=SV, HR=heart["HR"], systole_share=heart["systole_share"], **windkessel)   # periodic state, closed form
print({k: round(s[k]/MMHG, 2) for k in ("p0", "pS", "pD", "p_mean", "map_third")}, "mmHg; tau =", round(s["tau"], 3), "s")
print(f"R × mean inflow = {windkessel['R']*SV/T/MMHG:.2f} mmHg = the mean (D10 step 9)")
""", explain=r"""
1. `ejection_half_sine` is our inflow (ours): a half sine lasting 35 % of the beat that delivers one stroke volume.
2. `windkessel_ode` integrates «16.49» with `solve_ivp`; `windkessel_pressure` evaluates the solution «16.52» by quadrature.
   They agree to about 10⁻⁸ mmHg.
3. `windkessel_summary` uses the closed periodic form of D10: the pressure settles to 106.0/74.5 mmHg, its true time mean
   is 89.7 mmHg = R × mean inflow, and the ⅓ rule «16.11» would say 85.0 mmHg — 5 % low (trap T15).""")
scratch(r"""
dt = T/400; n = 8*400                                         # RK4 with 400 steps per beat, 8 beats
p_rk = np.empty(n + 1); p_rk[0] = 60*MMHG                     # start from 60 mmHg
f = lambda tt, pp: (Qin(tt) - pp/windkessel["R"])/windkessel["K"]   # dp/dt from (16.49)
for i in range(n):                                            # classical Runge–Kutta (Ch. 10)
    ti, pi_ = i*dt, p_rk[i]                                   # time [s] and pressure [Pa] at the start of step i
    k1 = f(ti, pi_); k2 = f(ti + dt/2, pi_ + dt*k1/2); k3 = f(ti + dt/2, pi_ + dt*k2/2); k4 = f(ti + dt, pi_ + dt*k3)   # four slopes dp/dt [Pa/s]: start, two midpoints, end
    p_rk[i + 1] = pi_ + dt*(k1 + 2*k2 + 2*k3 + k4)/6          # step with their weighted mean (weights 1, 2, 2, 1)
t_rk = np.arange(n + 1)*dt                                    # its times
assert np.allclose(p_rk, ch16.windkessel_ode(t_rk, Qin, windkessel["R"], windkessel["K"], 60*MMHG), rtol=1e-6)   # same as the library's solve_ivp run
tl = np.linspace(0, T, 2001)                                  # one beat of the periodic state
pl = ch16.windkessel_periodic_state(tl, SV=SV, HR=heart["HR"], systole_share=heart["systole_share"], **windkessel)
assert abs(pl[-1] - pl[0]) < 1e-6*pl[0]                       # periodic: p(T) = p(0)
assert np.isclose(ch16.time_mean(pl, tl), windkessel["R"]*SV/T, rtol=1e-6)   # mean = R Q̄
print("RK4 march, periodicity and mean all agree")
""", explain=r"""A 20-line Runge–Kutta loop (the classical fourth-order scheme of Ch. 10: four slope samples per step) marches
«16.49» by hand and matches the library to 10⁻⁶; the closed periodic state repeats after one beat and its mean is R times
the mean inflow.""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 4.4), gridspec_kw=dict(width_ratios=[1.6, 1]))
t3 = np.linspace(0, 3*T, 1500)                                # three periodic beats
p3 = ch16.windkessel_periodic_state(t3, SV=SV, HR=heart["HR"], systole_share=heart["systole_share"], **windkessel)
Q3 = Qin(t3)                                                  # inflow [m³/s]
ax.plot(t3, p3/MMHG, color=COLORS["orange"], lw=2, label="aortic pressure (Windkessel)")
ax.plot(t3, Q3*windkessel["R"]/MMHG, color=COLORS["muted"], ls="--", label="K → 0 (rigid): p = QR")
ax.axhline(s["p_mean"]/MMHG, color=COLORS["amber"], label="true mean")
ax.axhline(s["map_third"]/MMHG, color=COLORS["amber"], ls="-.", label="⅓ rule (16.11)")
ax2 = ax.twinx(); ax2.plot(t3, Q3*1e6, color=COLORS["blue"], lw=1, alpha=0.6); ax2.set_ylabel("inflow Q [mL/s] (blue)")
ax.set(xlabel="t [s]", ylabel="p [mmHg]", ylim=(0, 300), title="The aorta keeps the pressure up between beats")
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)   # legend under the panel, clear of every curve
Ks = np.linspace(0.3, 4, 40)/(1e6*MMHG)                      # compliances 0.3 … 4 mL/mmHg [m³/Pa]
ss = [ch16.windkessel_summary(SV=SV, HR=heart["HR"], systole_share=heart["systole_share"], R=windkessel["R"], K=k) for k in Ks]   # periodic state for each K
bx.plot(Ks*1e6*MMHG, [x["pS"]/MMHG for x in ss], color=COLORS["orange"], label="systolic p_S")
bx.plot(Ks*1e6*MMHG, [x["pD"]/MMHG for x in ss], color=COLORS["orange"], ls="--", label="diastolic p_D")
bx.plot(Ks*1e6*MMHG, [x["p_mean"]/MMHG for x in ss], color=COLORS["amber"], label="mean")
bx.axvline(windkessel["K"]*1e6*MMHG, color=COLORS["muted"], ls=":")
bx.set(xlabel="compliance K [mL/mmHg]", ylabel="p [mmHg]", title="Stiffer aorta, wider swing"); bx.legend(fontsize=7)
savefig(fig, "ch16", "c07_windkessel"); plt.show()
""", see=r"""Left: three periodic beats of our Windkessel — the inflow (blue, right axis), the aortic pressure (orange), what a
rigid aorta would give (grey dashed: p = QR, zero in diastole), the true mean (amber) and the ⅓-rule estimate (dash-dot).
Right: systolic, diastolic and mean pressure against the compliance at fixed R; the dotted line is our K.""",
    read=r"""The pressure never falls to zero: during diastole it decays by $e^{-(T-T_s)/RK}=0.71$. On the right the mean is
flat — it is R times the mean inflow, whatever K — while the swing (the pulse pressure) widens as K falls.""",
    change=r"""…the aorta stiffened with age (K halved)? The systolic pressure would rise and the diastolic fall at the same
mean — the isolated systolic hypertension of older people.""")
note("N48 [B]", r"""(From §16.2.) Along the circulation the large swings of the ventricle are damped: the ventricle's pressure
falls near zero each beat while the aorta stays high — our curves (a schematic ventricle and the Windkessel aorta), labelled
ours.""")
figure(r"""
p_ao = ch16.windkessel_periodic_state(t3, SV=SV, HR=heart["HR"], systole_share=heart["systole_share"], **windkessel)   # aorta [Pa]
phase_ = (t3 % T)                                             # time within the beat [s]
p_lv = np.where(phase_ < Ts, np.maximum(p_ao, 5*MMHG + (p_ao - 5*MMHG)*np.abs(np.sin(np.pi*phase_/Ts))**0.3),   # ejection: ≈ aortic
                5*MMHG + 5*MMHG*(phase_ - Ts)/(T - Ts))      # filling: 5 → 10 mmHg (schematic, ours)
fig, ax = plt.subplots(figsize=(7.5, 3.2))
ax.plot(t3, p_lv/MMHG, color=COLORS["rose"], label="left ventricle (schematic, ours)")
ax.plot(t3, p_ao/MMHG, color=COLORS["orange"], lw=2, label="aorta (Windkessel, ours)")
ax.set(xlabel="t [s]", ylabel="p [mmHg]", title="The ventricle swings to zero; the aorta does not"); ax.legend(fontsize=8)
savefig(fig, "ch16", "c07_ventricle_aorta"); plt.show()
""", see=r"""A schematic left-ventricle pressure (rose) and the Windkessel aortic pressure (orange) over three beats (both ours).""",
    read=r"""While the aortic valve is open the two pressures nearly coincide; when it shuts the ventricle drops to a few
mmHg to refill, while the aorta, holding its stored blood, decays only slowly.""",
    change=r"""…the aortic valve leaked? Blood would flow back into the ventricle in diastole and the aortic pressure would fall
faster — a wider pulse pressure.""")
P("P395", "`np.searchsorted`", r"""`np.searchsorted(t, x)` takes a **sorted** array `t` and returns the index at which `x` would
have to be inserted to keep it sorted — the index of the first element ≥ x. The animations below use it to find how much of
a fine time grid has elapsed at each (coarser) frame time, and draw the trace up to there.""", code=r"""
t_demo = np.array([0.0, 0.5, 1.0, 1.5])        # a sorted time grid [s]
print(np.searchsorted(t_demo, 0.7), np.searchsorted(t_demo, 1.0))   # 2 2: 0.7 s and 1.0 s both fall at index 2
""")
animation(r"""
nfr = 32 if FAST else 48                                      # frames over eight beats
tA = np.linspace(0, 8*T, nfr)                                 # frame times [s]
pA_t = ch16.windkessel_ode(t, Qin, windkessel["R"], windkessel["K"], 60*MMHG)   # the run from 60 mmHg (as above)
terms = ch16.windkessel_terms(tA % T, SV=SV, HR=heart["HR"], systole_share=heart["systole_share"], **windkessel)   # periodic split
fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(8, 3.0), gridspec_kw=dict(width_ratios=[2.2, 1, 1]))
ax.plot(t, pA_t/MMHG, color=COLORS["grid"]); (tr,) = ax.plot([], [], color=COLORS["orange"], lw=2)
ax.set(xlabel="t [s]", ylabel="p [mmHg]", xlim=(0, 8*T), ylim=(40, 120), title="Settling into the beat")
bx.set(xlim=(-1.3, 1.3), ylim=(-1.3, 1.3), aspect="equal"); bx.axis("off"); bx.set_title("aortic balloon", fontsize=9)
ball = plt.Circle((0, 0), 1.0, color=COLORS["teal"], alpha=0.5); bx.add_patch(ball)
bars = cx.bar([0, 1], [0, 0], color=[COLORS["teal"], COLORS["blue"]]); cx.set_xticks([0, 1], ["storage\nK dp/dt", "run-off\np/R"], fontsize=7)
cx.set(ylim=(-150, 350), ylabel="mL/s"); cx.set_title("inflow split (periodic)", fontsize=9)
def update(i):                                                # frame i
    k = np.searchsorted(t, tA[i])                             # index of the first t ≥ this frame's time (primer P395 above)
    tr.set_data(t[:k + 1], pA_t[:k + 1]/MMHG)
    ball.set_radius(0.6 + 0.5*(pA_t[k]/MMHG - 40)/80)         # bigger when the pressure (stored volume) is higher
    for b, v in zip(bars, (terms["storage"][i], terms["runoff"][i])): b.set_height(v*1e6)
    return (tr, ball, *bars)
fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it (faster frames)
show_animation(animate(update, frames=nfr, fig=fig, interval=90), player="video")
""", explain=r"""**Figure A1.** **What you see.** Eight beats starting from 60 mmHg: the pressure trace climbs and settles into the
periodic beat within about three time constants; the balloon breathes with the pressure; the bars split the inflow into
storage (K dp/dt, negative in diastole: the balloon gives blood back) and run-off (p/R, always positive). **How to read
it.** In diastole the inflow is zero, so storage = −run-off: the aorta feeds the periphery. **What would change if…** the
aorta were stiffer: RK shorter, the transient settles faster and the swing is larger.""")
live_cell(r"""
def wk_live(R=1.45, K=1.6, HR=66.0, SV_mL=75.0):              # R in 1e8 Pa s/m³, K in mL/mmHg, HR in beats/min, SV in mL
    Kp = K/(1e6*MMHG)                                         # compliance in m³/Pa
    sm = ch16.windkessel_summary(SV=SV_mL*1e-6, HR=HR, systole_share=0.35, R=R*1e8, K=Kp)   # periodic state
    Tb = 60/HR; tt = np.linspace(0, 3*Tb, 900)                # three beats
    pp = ch16.windkessel_periodic_state(tt, SV=SV_mL*1e-6, HR=HR, systole_share=0.35, R=R*1e8, K=Kp)
    fig, ax = plt.subplots(figsize=(6.5, 3)); ax.plot(tt, pp/MMHG, color=COLORS["orange"])
    ax.set(xlabel="t [s]", ylabel="p [mmHg]", ylim=(0, 250),
           title=f"p_S/p_D = {sm['pS']/MMHG:.0f}/{sm['pD']/MMHG:.0f} mmHg, mean {sm['p_mean']/MMHG:.0f}, ⅓ rule error {100*(sm['map_third']/sm['p_mean'] - 1):+.1f} %")
    plt.show()                                                # draw this setting's figure
import ipywidgets as W                                        # slider widgets (kernel only); each opens on our resting values
sl = lambda lo, hi, st, v: W.FloatSlider(min=lo, max=hi, step=st, value=v, continuous_update=False)   # one slider: range lo…hi, step st, opening value v
_ = live(wk_live, R=sl(0.5, 3.0, 0.05, 1.45), K=sl(0.3, 4.0, 0.1, 1.6), HR=sl(40.0, 180.0, 2.0, 66.0), SV_mL=sl(40.0, 120.0, 1.0, 75.0))   # four sliders re-run wk_live; `_ =` keeps the cell from printing an object
""", explain=r"""Four sliders (kernel only): R, K, heart rate and stroke volume; the title reports the periodic systolic/diastolic
pressure, the mean and the error of the ⅓ rule. The same idea runs on the web page in figure A1 above and in the explainer
below.""")
note("N47 [C]", r"""The model assumes the pressure is the same everywhere in the arteries at each instant — no travel time, no
reflection. C08 gives the pulse a finite speed; C13 adds reflections.""")
explainer("windkessel_heartbeat", "Why doesn't blood pressure fall to zero between beats?",
          r"""one clock drives the ventricle, the aortic balloon and the pressure trace; you change K and R and watch the beat
settle into a new rhythm, with the inflow split into storage and run-off at every instant — a static figure shows one
setting, the explainer shows the settling.""",
          ["Press play and watch the pressure settle within a few beats.", "Halve K (stiff aorta): the pulse widens, the mean stays put.",
           "Raise the heart rate to 120: diastole shortens and p_D climbs."])
whatif(r"""…the aorta were long enough for the pressure to take time to travel along it? Then the pressure would be a wave, not a
single number — C08.""")

# ---------------------------------------------------------------------------------------------------------------------
# C08 — Moens–Korteweg
# ---------------------------------------------------------------------------------------------------------------------
core("C08", r"The pulse is a wave of wall stretch: «16.62»", "Is the pulse you feel at your wrist the blood arriving?")
problem(r"""
Press your fingers on your wrist: the beat arrives about a tenth of a second after the heart's. The blood itself takes many
seconds to get there. **What travels is a bulge**: the heart's push stretches the aorta, the stretched wall squeezes the
blood ahead, which stretches the next piece of wall. A stiffer wall passes the push on faster.""")
idea(r"""
wall stiffness (a spring)  +  blood inertia (a mass)  ──►  a wave,   c² = stiffness / inertia
""")
note("N49 [B]", r"""Model: inviscid, incompressible blood in a long thin elastic tube; long waves; the pressure uniform over each
section. Area continuity «16.53» — the area is now an unknown (Ch. 15's duct laws in time).""")
recap_md("The material derivative in one dimension", "Ch. 3", r"""$Du/Dt=\partial u/\partial t+u\,\partial u/\partial x$: the
acceleration of a fluid particle, its local change plus the change it meets by moving.""")
note("N50 [B]", r"""Momentum per slice: «16.54».""")
slip("S5", r"«16.54p»", r"«16.54»", r"""The printed form misses the push of the sloping wall. Linearised, both give «16.57»,
so the wave speed is unaffected; the nonlinear run below uses the corrected form.""")
note("N51 [B]", r"""The tube law in A: «16.55» — the inverse of «16.41» with A = πa².""")
code(r"""
A0a = np.pi*aorta["a0"]**2                                  # undisturbed area of our aorta [m²]
a_test = ch16.tube_radius_hookean(25*MMHG, 0.0, aorta["a0"], aorta["E"], aorta["h"])   # (16.41) at 25 mmHg
print(ch16.tube_law_hookean(np.pi*a_test**2, A0a, aorta["E"], aorta["h"], aorta["a0"])/MMHG, "mmHg")   # (16.55) gives 25 back
""")
remind("C08")
D("D11", r"The linearised tube equations", ref="16.58", notes=r"**N52 [B]** what was dropped and how small it is (steps 3–5).")
code(r"""
print(ch16.linearisation_ratios(25*MMHG, 100*MMHG, 0.65, 4.86, Eh_over_a0=50e3))   # pulse/mean, u/c, strain: 0.25, 0.13, 0.067
""")
D("D12", r"The wave equation and the Moens–Korteweg speed", ref="16.62")
worked("the pulse speed in an aorta", r"""
1. E = 0.72 MPa, h/a₀ = 0.1, ρ = 1000 kg/m³ (easy numbers): c² = 0.72 × 10⁶ × 0.1/2000 = 36, c = **6 m/s**.
2. Our aorta (E = 0.4 MPa, h/a₀ = 0.125, ρ = 1060): c = √(0.4 × 10⁶ × 0.125/2120) = **4.86 m/s**.
3. Over 0.5 m to the wrist: 0.10 s.
4. Under a 25 mmHg crest the blood moves at p/(ρc) = 3333/(1060 × 4.86) = **0.65 m/s**, only while the crest passes.""")
code(r"""
c0 = ch16.moens_korteweg_speed(aorta["E"], aorta["h"], blood["rho"], aorta["a0"])   # (16.62) [m/s]
c_mid = ch16.moens_korteweg_speed(aorta["E"], aorta["h"], blood["rho"], aorta["a0"], midwall=True)   # (16.63) thick wall
dpdA = aorta["E"]*aorta["h"]/(2*aorta["a0"]*A0a)             # slope of the linear tube law (16.58) [Pa/m²]
c_bh = ch16.bramwell_hill_speed(A0a, dpdA, blood["rho"])     # c² = (A/ρ) dp/dA, the general form
E_back = ch16.youngs_modulus_from_pwv(c0, aorta["h"], blood["rho"], aorta["a0"])   # invert (16.62) from a measured speed
u_crest = ch16.velocity_from_pulse(25*MMHG, blood["rho"], c0)   # (16.68): u = p/(ρc) under a 25 mmHg crest
print(f"c0 = {c0:.4f} m/s, thick-wall {c_mid:.4f} m/s, Bramwell–Hill {c_bh:.4f} m/s, E back = {E_back:.4g} Pa, u = {u_crest:.3f} m/s")
""")
note("N53 [B]", r"""Thick walls: evaluate the strain at mid-wall, «16.63» (result given) — 3 % slower for h/a₀ = 0.125.
Inverting «16.62» gives E from a measured transit time (`youngs_modulus_from_pwv`).""")
recap("R09", "d'Alembert's solution", r"""«15.15»: any two shapes sliding apart at ±c. **New here:** the book writes shapes in
time, «16.65» — a pressure trace recorded at one place.""", where="Ch. 15 §15.2 (C02, D03)", core_id="C08")
code(r"""
F = lambda s: np.exp(-(s/0.1)**2)                           # a space shape [–] of width 0.1 m
f1 = lambda tau: F(-c0*tau)                                 # the same pulse as a time shape: f(x − ct) = f1(t − x/c)
xg = np.linspace(-0.5, 1.0, 7)                              # positions [m]
print(np.allclose(ch16.pulse_wave(xg, 0.1, c0, f1), AC.dalembert(xg, 0.1, c0, f=F)))   # (16.65) = (15.15): True
""")
recap("R10", "Pressure and particle velocity", r"""«15.16». **New here:** «16.68» with c the tube's wave speed — the seed of the
admittance of C13.""", where="Ch. 15 §15.2 (C02, D04)", core_id="C08")
scratch(r"""
h_ = 1e-4                                                     # step for finite differences [m] and [s]
X_, T_ = np.meshgrid(np.linspace(-0.3, 0.9, 61), np.linspace(0.02, 0.12, 21))   # a grid of (x, t)
p = lambda x, t: F(x - c0*t)                                  # a right-going pulse
p_tt = (p(X_, T_ + h_) - 2*p(X_, T_) + p(X_, T_ - h_))/h_**2  # second time difference
p_xx = (p(X_ + h_, T_) - 2*p(X_, T_) + p(X_ - h_, T_))/h_**2  # second space difference
assert np.max(np.abs(p_tt - c0**2*p_xx)) < 1e-4*np.max(np.abs(p_tt))   # the wave equation (16.61) holds
assert np.allclose(ch16.pulse_wave(X_, T_, c0, f1), p(X_, T_))          # and the library's (16.65) is this pulse
print("residual of p_tt − c²p_xx:", np.max(np.abs(p_tt - c0**2*p_xx))/np.max(np.abs(p_tt)))
""", explain=r"""Second differences of a sliding Gaussian on a grid (`np.meshgrid` builds every (x, t) pair) show
$p_{tt}=c^2p_{xx}$ to round-off; the library's «16.65» returns the same pulse.""")
figure(r"""
Es = np.logspace(5, np.log10(3e6), 200)                       # Young's modulus 0.1 … 3 MPa
fig, ax = plt.subplots(figsize=(7, 3.8))
for ratio, col in ((0.05, "#93c5fd"), (0.125, COLORS["blue"]), (0.25, "#1e3a8a")):   # wall-thickness ratios h/a0
    ax.semilogx(Es/1e6, ch16.moens_korteweg_speed(Es, ratio*aorta["a0"], blood["rho"], aorta["a0"]), color=col, label=f"h/a₀ = {ratio}")
    ax.semilogx(Es/1e6, ch16.moens_korteweg_speed(Es, ratio*aorta["a0"], blood["rho"], aorta["a0"], midwall=True), "--", color=col, lw=0.8)   # dashed: the thick-wall form (16.63)
ax.axvspan(0.2, 0.6, color=COLORS["teal"], alpha=0.07); ax.text(0.22, 13, "young … older aorta\n(illustrative)", fontsize=8)
ax.plot(aorta["E"]/1e6, c0, "o", color=COLORS["orange"]); ax.annotate(f"our aorta {c0:.2f} m/s", (aorta["E"]/1e6, c0), xytext=(5, -12), textcoords="offset points", fontsize=8)   # orange dot: our aorta, labelled with its speed
ax.set(xlabel="Young's modulus E [MPa] (log)", ylabel="pulse speed c [m/s]", title="Stiffer, thicker walls carry faster pulses (dashed: thick wall)")
ax.legend(fontsize=8)
savefig(fig, "ch16", "c08_pulse_speed"); plt.show()
""", see=r"""The Moens–Korteweg speed against Young's modulus for three wall-thickness ratios (solid), the thick-wall form
«16.63» dashed, an illustrative band of aortic stiffnesses, and our aorta (orange dot).""",
    read=r"""c ∝ √E: each factor 4 in stiffness doubles the speed; a thicker wall (larger h/a₀) is stiffer for the same E.
Pulse speeds of a few to ten metres per second — far faster than the blood itself.""",
    change=r"""…E quadrupled (a stiff, old aorta)? The pulse speed doubles, the reflected wave returns earlier (C13) and the
systolic peak rises — pulse-wave velocity is a clinical measure of arterial stiffness.""")
animation(r"""
nfr = 24 if FAST else 40                                      # frames
xs_ = np.linspace(0, 1.0, 300)                                # 1 m of artery [m]
tt_ = np.linspace(0, 0.16, nfr)                               # 0.16 s: the crest travels about 0.8 m
amp = 25*MMHG                                                 # crest pressure [Pa]
st0 = ch16.pulse_state(xs_, 0.0, c=c0, amplitude=amp, width=0.08, rho=blood["rho"], a0=aorta["a0"], E=aorta["E"], h=aorta["h"], x0=0.15)
xp = np.linspace(0.05, 0.95, 19)                              # tracer particles in the blood [m]
fig, ax = plt.subplots(figsize=(7.5, 3.0))
ex = 4.0                                                      # wall bulge drawn ×4
(top,) = ax.plot([], [], color=COLORS["teal"], lw=2); (bot,) = ax.plot([], [], color=COLORS["teal"], lw=2)
dots = ax.scatter(xp, np.zeros_like(xp), s=14, color=COLORS["blue"])
ax.set(xlim=(0, 1), ylim=(-2.2, 2.2), xlabel="x [m]", ylabel="r / a₀", title=f"The bulge runs at c = {c0:.2f} m/s; the blood barely moves (bulge ×{ex:.0f})")
xp_now = xp.copy(); prev_t = 0.0                              # tracer positions now [m] and the time of the last frame [s]
def update(i):                                                # frame i
    global xp_now, prev_t
    st = ch16.pulse_state(xs_, tt_[i], c=c0, amplitude=amp, width=0.08, rho=blood["rho"], a0=aorta["a0"], E=aorta["E"], h=aorta["h"], x0=0.15)
    y = 1 + ex*(st["a"]/aorta["a0"] - 1)                      # exaggerated wall radius
    top.set_data(xs_, y); bot.set_data(xs_, -y)
    u_p = np.interp(xp_now, xs_, st["u"])                     # blood speed at each tracer [m/s]
    xp_now = xp_now + u_p*(tt_[i] - prev_t); prev_t = tt_[i]  # move the tracers (Euler step)
    dots.set_offsets(np.c_[xp_now, np.zeros_like(xp_now)])
    return top, bot, dots
fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it (faster frames)
show_animation(animate(update, frames=nfr, fig=fig, interval=60), player="video")
print(f"true radius change under the crest: {(st0['a'].max()/aorta['a0'] - 1)*100:.1f} % (drawn ×{ex:.0f})")
shift = st0["u"].max()*0.08*np.sqrt(np.pi)/c0                 # tracer shift = ∫u dt = (u_max/c)∫e^{−(x/w)²}dx = u_max w√π/c [m]
print(f"blood speed under the crest {st0['u'].max():.3f} m/s; each tracer moves on by {shift*100:.1f} cm as the pulse passes")
""", explain=r"""**Figure A2.** **What you see.** A 25 mmHg pressure pulse travelling along 1 m of our aorta: the teal walls bulge
(drawn four times larger than the true 7 % — the printed line gives the true size) and the bulge runs at c = 4.86 m/s; the
blue tracers in the blood shift forward by about 2 cm (printed: $u_{max}w\sqrt\pi/c$ = 0.647 × 0.08 × 1.77/4.86 m) and
stop. **How to read it.** The wave moves, the blood hardly does: what reaches your wrist first is the bulge. **What would
change if…** the wall were four times stiffer: the bulge would run twice as fast ($c\propto\sqrt E$) and bulge only a
quarter as much for the same pressure (the radius change $a_0^2(p-p_e)/(Eh)$ is ∝ 1/E).""")
note("N54 [B]", r"""Measured pulses travel at a few metres per second in the aorta and faster in stiffer peripheral arteries;
waves in collapsible veins are slower still (C15); real pulses peak by reflection (C13), steepen (nonlinearity) and damp
(C11). Below, a large pulse in our nonlinear 1-D model (ours — **DEVIATION**: MacCormack's scheme of Ch. 10 with a small
artificial viscosity; the book stops at the linear wave equation).""")
recap_md("Riemann invariants", "Ch. 15 P358; characteristics Ch. 7 P174", r"""In Ch. 15 a gas carried $u\pm\dfrac{2c}{\gamma-1}$
unchanged along its characteristics. The same idea holds for our elastic tube: the quantity $u\pm\displaystyle\int\frac{c\,dA}{A}$
keeps its value along the characteristics $\dfrac{dx}{dt}=u\pm c$. **Where the 4 comes from:** the Hookean tube law «16.55»
gives $c=c_0\left(\dfrac{A_0}{A}\right)^{1/4}$, so $dc/c=-\tfrac14\,dA/A$, i.e. $\dfrac{c\,dA}{A}=-4\,dc$ and
$\displaystyle\int\frac{c\,dA}{A}=-4c$ (+ a constant). The invariant on the left-going characteristics $dx/dt=u-c$ is therefore
$u-(-4c)=u+4c$. If we start with $u+4c$ the same everywhere, every left-going characteristic carries the same value, so
**no left-going wave exists** — the pulse below is purely right-going. Check with easy numbers: at the foot $c=c_0$ and
$u=0$, so $u+4c=4c_0$; under the crest c is smaller, so $u=4(c_0-c)>0$ — the blood moves forward under the bulge.""")
figure(r"""
nx_ = 200 if FAST else 400                                    # grid points
x1 = np.linspace(0, 2.0, nx_)                                 # 2 m of tube [m]
A_in = A0a*(1 + 0.6*np.exp(-((x1 - 0.4)/0.08)**2))            # a large bulge: area up 60 % (ours)
p_in = ch16.tube_law_hookean(A_in, A0a, aorta["E"], aorta["h"], aorta["a0"])   # its pressure (16.55) [Pa]
c_in = c0*(A0a/A_in)**0.25                                    # local wave speed of the Hookean tube law: c = c₀(A₀/A)^{1/4} [m/s]
u_in = 4*(c0 - c_in)                                          # Riemann invariant u + 4c = 4c₀ everywhere (u + ∫c dA/A with ∫c dA/A = −4c), carried by the left-going characteristics u − c: uniform → no left-going wave [m/s]
run = ch16.tube1d_run(x1, A_in, u_in, A0a, aorta["E"], aorta["h"], aorta["a0"], blood["rho"], 0.22, snapshots=(0.0, 0.11, 0.22))
fig, ax = plt.subplots(figsize=(7.5, 3.2))
for k, col in zip(range(3), (COLORS["muted"], COLORS["orange"], COLORS["rose"])):   # the three snapshot times, grey → orange → rose
    ax.plot(x1, run["snapshots"]["p"][k]/MMHG, color=col, label=f"t = {run['snapshots']['t'][k]:.2f} s")
ax.set(xlabel="x [m]", ylabel="p − p_e [mmHg]", title="A large pulse steepens at its front (nonlinear 1-D model, ours)")
ax.legend(fontsize=8)
savefig(fig, "ch16", "c08_steepening"); plt.show()
lam = u_in + c_in                                              # characteristic speed u + c of each point of the start [m/s]
t_break = -1/np.min(np.gradient(lam, x1))                      # when the faster crest first catches the slower foot [s]
k_c = np.argmax(p_in)                                          # index of the crest
print(f"crest: c = {c_in[k_c]:.2f} m/s (c/c₀ = {c_in[k_c]/c0:.2f}), u = {u_in[k_c]:.2f} m/s, u + c = {lam[k_c]:.2f} m/s; foot: u + c = c₀ = {c0:.2f} m/s")
print(f"the front breaks at t ≈ {t_break:.3f} s; peak p − p_e at t = 0, 0.11, 0.22 s:", np.round(run["snapshots"]["p"].max(axis=1)/MMHG, 1), "mmHg")
print(f"mass kept: {run['mass']/run['mass0'] - 1:.1e} (relative)")
""", see=r"""Pressure along 2 m of our aorta at three times for a large right-going pulse in the nonlinear 1-D model (ours). The
start is a pure right-going wave. 🔁 Riemann invariant (Ch. 15 P358, which covered a gas, $u\pm2c/(\gamma-1)$):
$u\pm\int c\,dA/A$ keeps its value along $dx/dt=u\pm c$ (Ch. 7 P174); here $c\propto A^{-1/4}$ gives $\int c\,dA/A=-4c$, so
$u+4c$ rides the left-going characteristics — make it uniform and no left-going wave exists. So the blood speed is set so
that $u+4c$ is the same everywhere, $u=4(c_0-c)$ (zero at the foot, where $c=c_0$), and nothing travels left.""",
    read=r"""The crest moves faster than the foot, so the front steepens and the back stretches. The reason is not that the
wave itself is faster in a wider tube: with the tube law «16.55» the wave speed is $c^2=\dfrac{Eh}{2\rho a_0}\left(\dfrac{A_0}{A}\right)^{1/2}$,
which **falls** as the tube widens (c/c₀ = 0.89 under the crest, where A/A₀ = 1.6). What a point of the pulse travels at is
$u+c$, the wave speed plus the blood's own speed, and u rises faster than c falls: about 6.5 m/s at the crest against
4.86 m/s at the foot (printed). The front breaks at about t = 0.06 s; after that it is a jump the scheme cannot resolve, and
the peak printed above its starting 78.5 mmHg (84–86 mmHg at t = 0.11 s, depending on the grid), the small wiggles just
behind it and the slight dip below zero just ahead of it are MacCormack's numerical ringing at that jump, not physics — the
exact solution never raises the crest. The later fall of the peak is real in this
inviscid model: the jump eats into the crest. Total volume is conserved to round-off.""",
    change=r"""…the pulse were small (a few mmHg)? u and the change of c would both be tiny, every point would travel at c₀, and
the pulse would keep its shape — the linear limit of D11–D12.""")
explainer("pulse_wave_artery", "Is the pulse at your wrist the blood arriving?",
          r"""you release a pulse and watch the wall bulge race ahead while the blood barely moves; stiffness, thickness and radius
change the speed live — the separation between the wave and the blood is a motion, which a still image cannot show.""",
          ["Click under the crest: the blood moves at 0.6 m/s, the bulge at 4.9 m/s.", "Quadruple E: the speed doubles.",
           "Release a bump at rest: two half-height pulses run apart."])
whatif(r"""…the blood were viscous and the tube rigid? The "wave" speed becomes infinite (the whole column moves together) and the
question becomes the profile across the tube — C09.""")

# ---------------------------------------------------------------------------------------------------------------------
# C09 — Womersley profile
# ---------------------------------------------------------------------------------------------------------------------
core("C09", r"Womersley's oscillating profile: «16.82»", "What does the velocity profile look like when the pressure gradient swings back and forth?")
problem(r"""
Push a syringe in and out slowly and the liquid in the needle flows back and forth with a parabolic profile, in step with
your hand. Do it fast in a wide tube and something else happens: **the middle moves like a solid slug that cannot keep up,
and only a thin layer near the wall feels the wall**. Which picture applies is decided by α (C03).""")
idea(r"""
α small (slow beat, narrow tube):  viscosity reaches the axis within one cycle
      ──► a Poiseuille parabola that rises and falls in step with the gradient
α large (fast beat, wide tube):    viscosity reaches only a wall layer about a/α thick
      ──► a flat core that sloshes like a slug, lagging the gradient by almost a quarter cycle
""", r"""**α compares the radius with the distance viscosity spreads in one cycle**: «16.12» is $a$ divided by
$\sqrt{\nu/\omega}$. Womersley's exact profile, derived below, joins the two pictures smoothly.""")
recap("R11", "Unsteady flow in a rigid tube", r"""Rigid, axisymmetric, long-wave flow with the convective term dropped: «16.71»
(Ch. 4 cylindrical Navier–Stokes; Ch. 8 parallel flows). **New here:** u depends on r and t and p on x and t, so ∂p/∂x can
depend on t only (result given). ⚠️ (trap T23) In a rigid tube this is **not a wave**: the whole column moves at once — waves
need the elastic wall (C08, C11).""", where="Ch. 4 §4.5; Ch. 8 (C03, C11)", core_id="C09")
remind("C09")
confusion(r"""(trap T1) **Two time conventions** in this chapter — never mix them in one formula:

| Problem | Convention | Bessel argument |
|---|---|---|
| rigid tube (this block, C10) | «16.73», «16.75»: $e^{+i\omega t}$ | $i^{3/2}\alpha$ |
| elastic tube (C11) | «16.116»: $e^{i(kx-\omega t)}$ | $\hat\beta=e^{i\pi/4}\alpha$ |

The two are complex conjugates of each other.""")
note("N55 [B]", r"""One harmonic of the gradient: «16.73» — the real part $-G\cos\omega t$ is meant (complex amplitudes, P176);
the trial «16.75».""")
P("P374", "Bessel's equation and J₀, Y₀, J₁", r"""The equation $y''+y'/z+y=0$ is what the Laplacian in a circle turns into for
oscillating problems. Its two solutions are $J_0(z)$, smooth at z = 0 (= 1 there), and $Y_0(z)$, infinite at z = 0. For real
z, J₀ wiggles and decays like a damped cosine; for a complex z such as $i^{3/2}\alpha$ it grows exponentially. Its power
series is $J_0(z)=\sum_m(-1)^m(z/2)^{2m}/(m!)^2$; $J_1$ is its partner, $J_0'=-J_1$.""", code=r"""
from scipy.special import jv, yv               # Bessel functions of the first and second kind
print(jv(0, [0, 2.405, 5.0]), yv(0, 1e-3))     # 1, ≈0 (first zero), −0.178; Y0 is −4.47 near 0 (infinite at 0)
""")
P("P377", "powers of i", r"""$i^{1/2}=e^{i\pi/4}=(1+i)/\sqrt2$ and $i^{3/2}=e^{3i\pi/4}=(-1+i)/\sqrt2$; also $i\cdot
i^{1/2}=i^{3/2}$ and $(i^{3/2})^2=i^3=-i$. numpy's `1j**1.5` takes the principal branch and gives the same number (P159, P45).""",
  code=r"""
print(1j**0.5, 1j**1.5, np.exp(3j*np.pi/4))     # (0.707+0.707j) (−0.707+0.707j) twice
print((1j**1.5)**2)                             # −1j
""")
P("P375", "scipy.special.jv with a complex argument", r"""`jv(n, z)` evaluates $J_n$ for complex z. For a complex z such as
$i^{3/2}\alpha$ the values grow fast ($|J_0(i^{3/2}\alpha)|$ grows like $e^{\alpha/\sqrt2}$): plain `jv` stays finite and accurate
up to α ≈ 950 ($|J_0|\approx7\times10^{289}$) and returns `inf` only near α ≈ 1000, where it passes the largest double,
about $1.8\times10^{308}$. **Our own conservative habit**: once α passes about 30 we use `jve(n, z)` =
$J_n(z)e^{-|\mathrm{Im}\,z|}$ and take ratios, so the huge factors cancel before they can overflow (trap T25).""", code=r"""
z = np.exp(3j*np.pi/4)*40.0                     # i^{3/2} α at α = 40
print(abs(jv(0, z)), abs(jve(0, z)))            # 1.21e11 and 0.0632: the same number, scaled by e^{−|Im z|} (e^{−28.3})
""")
D("D13", r"Womersley's profile", ref="16.82", notes=r"**N56 [B]** the Bessel equation «16.76» (step 4) and its solution «16.77»–«16.82».")
slip("S6", r"«16.76p»", r"«16.76»", r"""The book's own solution, with the particular part $G/(\omega\rho i)$ of «16.77», needs the
minus sign on the right.""")
code(r"""
res = ch16.womersley_ode_residual(0.5e-3, 1e-3, 2*np.pi, 1060.0, 3.5e-3, 1.0)   # left side of (16.76) for the coded profile, G = 1 Pa/m
print(res, "should be −G/μ =", -1.0/3.5e-3)                  # −285.7: the corrected right side
""")
worked("the centreline at α = 3", r"""
1. z = i^{3/2} × 3: z² = −9i.
2. J₀(z) from the series: 1 + 9i/4 − 81/64 − 729i/2304 + … = −0.2214 + 1.9376i (the series converges fast).
3. On the axis J₀(0) = 1, so F₁ = 1/J₀(z) = −0.0582 − 0.5095i.
4. 1 − F₁ = 1.0582 + 0.5095i.
5. u/(Ga²/4μ) = (4/(9i))(1.0582 + 0.5095i) = 0.2264 − 0.4703i: amplitude **0.522** of Poiseuille's centre speed, lagging by
   arctan(0.4703/0.2264) = **64.3°**.""")
code(r"""
F1_axis = ch16.womersley_F1(0.0, 1.0, 3.0)                  # (16.81) on the axis at α = 3
u_c = ch16.womersley_velocity_nd(0.0, 0.0, 3.0)             # Re[u e^{iωt}]/(Ga²/4μ) on the axis at ωt = 0
u_c90 = ch16.womersley_velocity_nd(0.0, np.pi/2, 3.0)       # … a quarter period later
amp, lag = np.hypot(u_c, u_c90), np.degrees(np.arctan2(u_c90, u_c))   # amplitude √(x²+y²) and lag; np.arctan2(y, x) = angle of (x, y) in the right quadrant (Ch. 2 P70)
print(f"F1 = {F1_axis:.4f}; centre: {u_c:.4f} at ωt = 0, amplitude {amp:.3f}, lag {lag:.1f}°")
r_ = np.linspace(0, artery["a"], 5)                          # radii of our artery [m]
print("u [mm/s] at ωt = 0 for G = 100 Pa/m:", np.round(ch16.womersley_velocity(r_, 0.0, artery["a"], omega, blood["rho"], blood["mu"], 100.0)*1e3, 3))
""", explain=r"""
1. `womersley_F1` is the Bessel ratio «16.81»; on the axis it is $1/J_0(i^{3/2}\alpha)$.
2. `womersley_velocity_nd` returns the real profile in units of Poiseuille's centre speed $Ga^2/4\mu$; evaluating it at ωt = 0
   and a quarter period later gives the real and (minus the) imaginary part, hence amplitude 0.522 and lag 64.3° — the tiny
   example.
3. `womersley_velocity` is the dimensional profile for our artery with a gradient amplitude of 100 Pa/m.""")
scratch(r"""
i32 = np.exp(3j*np.pi/4)                                       # i^{3/2}
J0 = lambda z: sum((-1)**m*(z/2)**(2*m)/math.factorial(m)**2 for m in range(60))   # power series of J0 (P374), complex z
s_ = np.linspace(0, 1, 11)                                     # r/a
for al in (0.5, 3.0, 8.0):                                     # three Womersley numbers
    u_mine = (4/(1j*al**2))*(1 - J0(i32*al*s_)/J0(i32*al))     # (16.82) in units of Ga²/4μ, complex
    assert np.allclose(u_mine.real, ch16.womersley_velocity_nd(s_, 0.0, al), rtol=1e-9, atol=1e-12)        # ωt = 0
    assert np.allclose(-u_mine.imag, ch16.womersley_velocity_nd(s_, np.pi/2, al), rtol=1e-9, atol=1e-12)   # ωt = π/2
    assert np.allclose(J0(i32*al), jv(0, i32*al), rtol=1e-12)  # the series is scipy's J0
@functools.lru_cache(maxsize=None)                             # remember the PDE run (P252): it runs once per session
def pde(nr):
    return ch16.womersley_pde_solve(1e-3, 2*np.pi, 1060.0, 3.5e-3, 100.0, nr=nr, periods=4 if FAST else 8)
sol = pde(80 if FAST else 200)                                 # Crank–Nicolson in r over several periods (ours)
RR, TT = np.meshgrid(sol["r"], sol["t"])                       # the same (r, t) points
u_bes = ch16.womersley_velocity(RR, TT, 1e-3, 2*np.pi, 1060.0, 3.5e-3, 100.0)   # the Bessel profile there
print(f"series = library at α = 0.5, 3, 8; PDE vs Bessel: max error {np.max(np.abs(sol['u'] - u_bes))/np.max(np.abs(u_bes)):.1e} of the peak")
""", explain=r"""
1. A 60-term power series of $J_0$ written by hand reproduces `scipy.special.jv` to 10⁻¹² and the library's profile at
   three Womersley numbers, at ωt = 0 and a quarter period later.
2. An independent route: march the PDE «16.72» itself in time with a Crank–Nicolson scheme (`womersley_pde_solve`, ours,
   cached with `functools.lru_cache` so a re-run is free) until it forgets its start; its last period matches the Bessel
   profile to its grid error (printed).""")
D("D14", r"The small-α (quasi-steady) limit", ref="16.87")
recap_md("Stokes' second problem", "Ch. 8 C11", r"""«8.38»: a wall oscillating in its own plane drags a layer of fluid that
decays and turns its phase over the depth $\sqrt{2\nu/\omega}$ (`LAM.stokes_second_problem`).""")
D("D15", r"The large-α limit by a thin wall layer", ref="16.88", notes=r"""**N58 [B]** the large-α profile «16.88»; for harmonic
n the layer scale is $\sqrt{n\omega/2\nu}$ and the harmonic sum stops at N (the book keeps ω and sums to ∞ — slip S7).""")
code(r"""
for al in (0.5, 20.0):                                       # one small and one large Womersley number
    w_ = al**2*blood["mu"]/(blood["rho"]*artery["a"]**2)     # the ω that gives this α in our artery [rad/s]
    rr = np.array([0.0, 0.9])*artery["a"]                    # axis and 0.9 a
    ex_ = ch16.womersley_velocity(rr, 0.0, artery["a"], w_, blood["rho"], blood["mu"], 100.0)   # exact (16.82)
    sm_ = ch16.womersley_small_alpha(rr, 0.0, artery["a"], w_, blood["rho"], blood["mu"], 100.0) # (16.87)
    lg_ = ch16.womersley_large_alpha(rr, 0.0, artery["a"], w_, blood["rho"], blood["mu"], 100.0) # (16.88)
    amp = lambda f: np.hypot(f(rr, 0.0), f(rr, np.pi/(2*w_)))   # amplitude from two instants a quarter period apart [m/s]
    ex_, sm_, lg_ = (amp(lambda r, t, F=F: F(r, t, artery["a"], w_, blood["rho"], blood["mu"], 100.0))   # amplitudes [m/s] of …
                     for F in (ch16.womersley_velocity, ch16.womersley_small_alpha, ch16.womersley_large_alpha))   # … exact (16.82), (16.87), (16.88)
    print(f"α = {al}: amplitude relative to exact at r = 0, 0.9a — small-α {np.round(sm_/ex_, 3)}, large-α {np.round(lg_/ex_, 3)}")
""", explain=r"""**What does this show?** Amplitudes relative to the exact profile: at α = 0.5 the quasi-steady parabola «16.87» is within
0.1 % and the wall-layer form «16.88» is useless; at α = 20 it is the other way round — the wall-layer form is within about 0.1 %
at these two radii (0.14 % at 0.9a) (it drops the wall's curvature, which matters most just inside the layer), the parabola is 20–100 times
too large.""")
note("N57 [B]", r"""A real pressure gradient is a sum of harmonics: «16.83»; because the equation is linear the velocity is the
sum of their profiles, «16.86» with $\alpha_n=\alpha\sqrt n$, and the mean (n = 0) is Poiseuille's «16.85» (result given).""")
code(r"""
G_n = [400.0, 900.0*np.exp(-0.3j), 450.0*np.exp(-1.6j), 200.0*np.exp(-2.6j)]   # a synthetic gradient (ours): mean + 3 harmonics [Pa/m]
rr = np.linspace(0, artery["a"], 4)                          # radii [m]
for wt in (0.0, np.pi/2, np.pi, 3*np.pi/2):                  # four instants of the beat
    print(f"ωt = {wt:.2f}: u [cm/s] =", np.round(ch16.womersley_harmonics(rr, wt/omega, artery["a"], omega, blood["rho"], blood["mu"], G_n)*100, 2))
""")
figure(r"""
s_ = np.linspace(-1, 1, 201)                                  # across the whole diameter r/a (trap T20: the book's figure runs wall–centre)
fig, axs = plt.subplots(1, 4, figsize=(11, 3.2), sharey=False)
for ax, al in zip(axs, (1.0, 4.0, 8.0, 14.0)):                 # our four Womersley numbers
    for k, wt in enumerate(np.radians(np.arange(0, 180, 15))):  # every 15° of phase over half a cycle
        ax.plot(ch16.womersley_velocity_nd(np.abs(s_), wt, al), s_, color=plt.cm.viridis(k/12), lw=1.2)   # one profile; colour runs dark → light with phase
    ax.axvline(0, color=COLORS["muted"], lw=0.6)
    ax.set(title=f"α = {al:g}", xlabel="u / (Ga²/4μ)"); axs[0].set_ylabel("r / a")
fig.suptitle("Small α: parabolas breathing in step; large α: flat cores lagging (half a cycle, every 15°)", fontsize=10)
savefig(fig, "ch16", "c09_profiles"); plt.show()
""", see=r"""**N59 [B]** Velocity profiles across the full diameter at four Womersley numbers (our own values), drawn every 15° of
phase over half a cycle (colour from dark to light), in units of Poiseuille's centre speed for the same gradient
amplitude.""",
    read=r"""At α = 1 the profiles are parabolas that grow and shrink in step with the gradient. At α = 8 and 14 the core is
flat, its speed is much smaller (the scale shrinks like 1/α²), and near the wall the profile reverses before the core does
— the wall layer turns first. Our figure shows the whole diameter; the book's runs wall–centre–wall mirrored (trap T20).""",
    change=r"""…the heart rate doubled? Every α grows by √2: α = 4 would look like α = 5.7 — flatter and later.""")
note("N60 [C]", r"""Read twice in the book: reversal near the wall, a flat core, an almost solid core at the largest α.""")
animation(r"""
nfr = 24 if FAST else 36                                      # frames for one cycle
wts = np.linspace(0, 2*np.pi, nfr, endpoint=False)            # phases of the gradient
fig, axs = plt.subplots(1, 3, figsize=(8, 3.0))
lines = []
for ax, al in zip(axs, (1.0, 4.0, 12.0)):                      # three Womersley numbers side by side
    umax = 1.05*max(np.max(np.abs(ch16.womersley_velocity_nd(np.abs(s_), w, al))) for w in wts[::6])   # fixed axis per panel
    (ln,) = ax.plot([], [], color=COLORS["blue"], lw=2); lines.append((ln, al))
    ax.set(xlim=(-umax, umax), ylim=(-1, 1), title=f"α = {al:g}", xlabel="u/(Ga²/4μ)"); ax.axvline(0, color=COLORS["muted"], lw=0.6)
axs[0].set_ylabel("r / a")
clock = fig.suptitle("", fontsize=10)                          # the gradient's phase, written out
def update(i):                                                # frame i
    for ln, al in lines:                                      # each panel's line and its Womersley number
        ln.set_data(ch16.womersley_velocity_nd(np.abs(s_), wts[i], al), s_)
    clock.set_text(f"gradient −∂p/∂x = G cos ωt at ωt = {np.degrees(wts[i]):.0f}° (largest at 0°)")
    return (*[l for l, _ in lines], clock)
fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it (faster frames)
show_animation(animate(update, frames=nfr, fig=fig, interval=60), player="video")
""", explain=r"""**Figure A3.** **What you see.** One full cycle of the profile at α = 1, 4 and 12 (each panel with its own fixed axis),
with the gradient's phase written above. **How to read it.** At α = 1 the parabola peaks at ωt = 0° (in
phase); at α = 12 the flat core peaks almost a quarter cycle later and the wall layer reverses first. **What would change
if…** α were 30: the core would lag by nearly 90° and the wall layer would be a thin skin.""")
explainer("womersley_profiles", "Is pulsing flow just Poiseuille flow sloshing?",
          r"""dragging α morphs the profile from a breathing parabola to a lagging plug while the phase runs; the two limits snap on
as ghosts and the vessel axis shows where each vessel sits — the morph between the limits is what a set of still curves
hides.""",
          ["Drag α from 0.5 to 17 and watch the profile flatten.", "Watch the purple near-wall trace turn before the centre.",
           "Click a radius to see the Bessel arithmetic."])
whatif(r"""…we asked for the total flow instead of the profile? Integrate «16.82» over the section: C10.""")

# ---------------------------------------------------------------------------------------------------------------------
# C10 — amplitude and phase
# ---------------------------------------------------------------------------------------------------------------------
core("C10", r"How much flow, how late: «16.99», «16.100»", "Does the flow follow the pressure gradient?")
problem(r"""
Push a child on a swing gently and often: the swing moves with your hand. Push a heavy cart back and forth quickly: it
barely moves, and its motion trails your push by a quarter cycle — your force goes into accelerating it. **Blood in a large
artery is the cart.**""")
idea(r"""
flow rate = pressure-gradient amplitude × (one complex number that depends only on α)
   its size   ──► how much flows
   its angle  ──► how late the flow peaks
α small:  Poiseuille's flux, in step with the gradient      (friction wins: the swing)
α large:  8/α² × Poiseuille's flux, a quarter cycle late    (inertia wins: the cart)
""", r"""**Add the oscillating profile over the cross-section and you get one complex number**: its size gives the amplitude
«16.99» and its angle the lag «16.100». With the module: the lag is 2° at α = 0.5, 79° at α = 8 and 86° at α = 20.""")
remind("C10")
P("P376", "two Bessel identities", r"""$\dfrac{d}{dz}J_0(z)=-J_1(z)$ and $\displaystyle\int_0^zsJ_0(s)\,ds=zJ_1(z)$. The first
gives the wall shear, the second the flow rate. Both hold for complex z.""", code=r"""
from scipy.integrate import quad               # numerical integral
z = 2.0; h = 1e-6
print((jv(0, z + h) - jv(0, z - h))/(2*h), -jv(1, z))      # −0.5767 twice
print(quad(lambda s: s*jv(0, s), 0, z)[0], z*jv(1, z))     # 1.1534 twice
""")
D("D16", r"The flow rate and F₂", ref="16.90", notes=r"**N61 [B]** the flux «16.90» with «16.91».")
D("D17", r"Amplitude and phase lag", ref="16.100")
worked("α = 8", r"""
1. F₂(8) = 0.1772 − 0.1608i (`womersley_F2`).
2. Z = 1 − F₂ = 0.8228 + 0.1608i: X = 0.8228, Y = 0.1608.
3. √(X² + Y²) = 0.8384.
4. |Q|/|Q₀| = 8 × 0.8384/64 = **0.105** — a tenth of the steady flow for the same gradient amplitude.
5. φ = arctan2(0.8228, 0.1608) = **78.9°**; at 66 beats per minute (ω = 6.91 rad/s) the flow peaks 1.378/6.91 = **0.199 s**
   after the gradient.

(The book's own worked example does not quite match «16.91» — slip S11; the comparison stays private; ours is computed.)""")
code(r"""
F2 = ch16.womersley_F2(8.0)                                  # (16.91) at α = 8
X8, Y8 = ch16.womersley_XY(8.0)                              # (16.97): 1 − F2 = X + iY
ratio8 = ch16.womersley_flow_ratio(8.0)                      # |Q|/|Q0| = 8√(X²+Y²)/α² (ours, from (16.99) and (16.95))
phi8 = ch16.womersley_phase_lag(8.0)                         # (16.100): lag of Q behind −∂p/∂x [rad]
print(f"F2 = {F2:.4f}, X = {X8:.4f}, Y = {Y8:.4f}, |Q|/|Q0| = {ratio8:.4f}, phi = {np.degrees(phi8):.2f}°, delay = {phi8/omega:.3f} s")
tq_ = np.linspace(0, 2*np.pi/omega, 5)                       # five instants of one beat [s]
print("Q(t) [mL/s] in our artery, G = 100 Pa/m:", np.round(ch16.womersley_flow_rate(tq_, artery["a"], omega, blood["rho"], blood["mu"], 100.0)*1e6, 4))
""")
scratch(r"""
rr = np.linspace(0, artery["a"], 2001)                        # radii [m]
for t0 in (0.0, 0.1, 0.37):                                   # three instants [s]
    u_ = ch16.womersley_velocity(rr, t0, artery["a"], omega, blood["rho"], blood["mu"], 100.0)   # the profile (16.82)
    Q_mine = np.trapezoid(2*np.pi*rr*u_, rr)                  # add the rings
    assert np.isclose(Q_mine, ch16.womersley_flow_rate(t0, artery["a"], omega, blood["rho"], blood["mu"], 100.0), rtol=1e-6)   # = (16.90)
print("ring sums equal the flux (16.90)")
""", explain=r"""At three instants the trapezoid sum of $2\pi r u$ over the profile equals the closed-form flux «16.90» (to the
trapezoid rule's error).""")
note("N62 [B]", r"""Wall shear: «16.92» (result given, from $dJ_0/dz=-J_1$).""")
slip("S8", r"this as a shear *rate* $\partial u/\partial r$ with an unnamed F", r"«16.92» — the wall shear *stress* $-\mu\,\partial u/\partial r$ with $F_2$")
code(r"""
a_, hh = artery["a"], 1e-7                                   # wall radius, small step [m]
for t0 in (0.0, 0.2):                                        # two instants [s]
    uu = ch16.womersley_velocity(np.array([a_, a_ - hh, a_ - 2*hh]), t0, a_, omega, blood["rho"], blood["mu"], 100.0)
    tau_fd = -blood["mu"]*(3*uu[0] - 4*uu[1] + uu[2])/(2*hh)  # −μ ∂u/∂r at the wall, one-sided second-order difference
    tau_lib = ch16.womersley_wall_shear(t0, a_, omega, blood["rho"], blood["mu"], [0, 100.0])   # (16.92, corrected)
    print(f"t = {t0}: FD {tau_fd:.6f} Pa, (16.92) {tau_lib:.6f} Pa")
""", explain=r"""**What does this show?** The finite-difference wall stress of the profile equals the corrected «16.92» — a stress
in Pa, not a shear rate. The wall shear leads the flow (it responds to the thin wall layer, which turns first).""")
note("N63 [B]", r"""Small α: «16.93», so the flow tends to «16.94» with magnitude «16.95» — Poiseuille's flux for the gradient
amplitude G. (The book writes the exponent as $e^{in\omega t}$, «16.94p», for this single mode — slip S9.)""")
code(r"""
als = np.array([0.25, 0.5, 1.0])                             # small Womersley numbers
err = np.abs(ch16.womersley_F2(als) - ch16.womersley_F2_series(als))   # error of 1 − iα²/8
print("errors:", err, "| slopes:", np.round(np.diff(np.log(err))/np.diff(np.log(als)), 2))   # ≈ 4: the next term is α⁴
""")
note("N64 [B]", r"""Large α: «16.96».""")
slip("S10", r"«16.96p»", r"«16.96»", r"""The printed real correction is worse than no correction at all (table below).""")
code(r"""
rows = []                                                    # one row of errors per α
for al in (10.0, 30.0, 100.0):                               # large Womersley numbers
    ex = ch16.womersley_F2(al)                                # exact (16.91)
    lead = 2/(1j**0.5*al)                                     # the leading term only
    rows.append(dict(alpha=al, leading=abs(lead - ex)/abs(ex), corrected=abs(ch16.womersley_F2_asymptote(al) - ex)/abs(ex),
                     printed=abs(ch16.womersley_F2_asymptote(al, printed_16_96=True) - ex)/abs(ex)))
display(pd.DataFrame(rows).style.format({"leading": "{:.1e}", "corrected": "{:.1e}", "printed": "{:.1e}"}))   # relative errors
""", explain=r"""**What does this show?** Relative error of three approximations of F₂: the leading term (∝ 1/α), the corrected
«16.96» (∝ 1/α², much better) and the printed form (larger than the leading term's: its "correction" points the wrong way).""")
figure(r"""
al_ = np.logspace(-1, np.log10(50), 300)                      # Womersley numbers 0.1 … 50
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.6))
ax.loglog(al_, ch16.womersley_flow_ratio(al_), color=COLORS["blue"], lw=2, label="|Q|/|Q₀|")
ax.loglog(al_, np.minimum(1, 8/al_**2), "--", color=COLORS["muted"], label="asymptotes 1 and 8/α²")
bx.semilogx(al_, np.degrees(ch16.womersley_phase_lag(al_)), color=COLORS["accent"], lw=2, label="lag φ")
bx.axhline(90, ls="--", color=COLORS["muted"]); bx.axhline(0, ls="--", color=COLORS["muted"])
for name, al0 in (("our artery", alpha_art), ("our aorta", ch16.womersley_number(aorta["a0"], omega, nu))):   # mark our two vessels
    ax.axvline(al0, color=COLORS["orange"], ls=":"); bx.axvline(al0, color=COLORS["orange"], ls=":")
    ax.text(al0*1.05, 2e-3, name, rotation=90, fontsize=8, color=COLORS["orange"])
ax.set(xlabel="Womersley number α (log)", ylabel="flow amplitude / Poiseuille's [–]", ylim=(1e-3, 2), title="A damper below α ≈ 1, a mass above α ≈ 10")
bx.set(xlabel="Womersley number α (log)", ylabel="lag of Q behind −∂p/∂x [°]", ylim=(-5, 95), title="…and up to a quarter cycle late")
ax.legend(fontsize=8); bx.legend(fontsize=8)
savefig(fig, "ch16", "c10_bode"); plt.show()
""", see=r"""Left: the flow amplitude relative to Poiseuille's for the same gradient amplitude, against α (log–log), with the
asymptotes 1 and 8/α². Right: the phase lag of the flow behind the gradient, from 0° to 90°. Dotted: our artery and our
aorta.""",
    read=r"""Below α ≈ 1 the flow is Poiseuille's and in step (a driven damper: viscosity limits it). Above α ≈ 10 the flow
falls as 8/α² and lags by almost 90° (a driven mass: inertia limits it). Our 2-mm artery sits in the transition; the aorta
is well inside the inertial regime.""",
    change=r"""…the vessel were twice as wide? α doubles: in the inertial regime the flow amplitude relative to Poiseuille's
falls four times — though Poiseuille's own flux grows sixteen times.""")
explainer("womersley_flow_phase", "Does the flow follow the pressure?",
          r"""dragging α moves a cursor on the amplitude and lag curves while the time traces slide apart and the complex number Z
turns — the link between the curve, the traces and the complex plane is the point.""",
          ["Set α = 1: flow and gradient peak together.", "Set α = 20: the flow is 2 % of Poiseuille's and 86° late.",
           "Switch the asymptote to the printed «16.96p» at α = 10 and compare the errors."])
whatif(r"""…the tube were elastic as well? Then the pressure travels (C08) and the viscous profile of C09 rides on it: C11.""")

# ---------------------------------------------------------------------------------------------------------------------
# C11 — viscous waves in an elastic tube
# ---------------------------------------------------------------------------------------------------------------------
core("C11", r"Where waves meet viscosity: the characteristic quadratic «16.136»",
     "How much do viscosity and the vessel wall slow the pulse and make it fade?")
problem(r"""
C08 gave the pulse a speed from the wall's stiffness, ignoring viscosity. C09 gave the oscillating profile, ignoring the
wall's give. A real artery has both: **the pulse drags a viscous layer along the wall, which slows it and makes it fade
with distance**. How much depends — once more — on α.""")
idea(r"""
wall stiffness  ──► sets the speed c₀ (C08)
viscous layer   ──► a complex wavenumber k = k₁ + i k₂: slower (c = ω/k₁ < c₀) and fading (e^{−k₂x})
α large: thin layer         ──► c ≈ c₀, little fading per wavelength     (α = 20: c/c₀ = 0.98)
α small: layer fills tube   ──► much slower, strongly damped              (α = 1:  c/c₀ = 0.53)
""", r"""**The wall sets the speed, the viscous layer trims it and damps it.** Both come out of one quadratic for the
unknown $X$ (which carries $k^2$), «16.136»; α decides how thick the layer is, so α decides the trim (numbers from
`ch16.elastic_tube_wave`).""")
confusion(r"""(trap T1, repeated) This block uses «16.116», the opposite sign to C09's «16.75»: here $\partial/\partial t\to
-i\omega$, so +iωρ/μ appears where C09 had −iωρ/μ. The two-row table of C09 again: rigid tube $e^{+i\omega t}$ with
$i^{3/2}\alpha$; elastic tube $e^{i(kx-\omega t)}$ with $\hat\beta=e^{i\pi/4}\alpha$.""")
recap("R12", "Continuity and Navier–Stokes in (r, x)", r"""«16.101» and «16.102», linearised. **New here:** the radial velocity v
is kept, because the wall moves.""", where="Ch. 4 §4.5; Appendix B", core_id="C11")
note("N65 [B]", r"""The wall is a membrane with displacements η (radial) and ξ (axial); its tensions follow Hooke's law with
Poisson's ratio: «16.104»; its radial motion «16.106» — the hoop tension's inward push T/a (P380) (stated).""")
note("N66 [B]", r"""With the fluid's stresses «16.110» and «16.111» put in, the wall equations become «16.112» and «16.113»; in
mode form they are D19's start, «16.126» and «16.127». The book's «16.108» prints 1 − ν̂ for 1 − ν̂² (slip S12). This chapter counts a
normal stress positive in compression here. **N67 [B]** No slip at the moving wall, applied at the undisturbed radius:
«16.114», «16.115» — Ch. 7's linearised boundary condition.""")
P("P385", "a complex wavenumber", r"""With $k=k_1+ik_2$, the wave $e^{i(kx-\omega t)}=e^{-k_2x}e^{i(k_1x-\omega t)}$ travels at
$c=\omega/k_1$ and its amplitude falls by $e^{-2\pi k_2/k_1}$ over one wavelength — the spatial twin of the complex
frequency of Ch. 11 (growth in time).""", code=r"""
k1, k2 = 1.0, 0.0622                           # per metre (α = 8 below)
print(np.exp(-2*np.pi*k2/k1))                  # 0.677: a third lost per wavelength
""")
D("D18", r"The long-wave normal-mode reduction", ref="16.123", notes=r"""**N68 [B]** the travelling-wave modes «16.116» with $|k|=\sqrt{k_1^2+k_2^2}$ (the
book prints $\sqrt{k_1+k_2}$ — slip S13); **N69 [B]** the long-wave ordering «16.117», «16.120».""")
D("D19", r"The reduced wall equations", ref="16.128", notes=r"**N70 [B]** «16.128», «16.129» (the shear that survives is $\mu\,d\hat u/dr$ — slip S14).")
D("D20", r"The mode solution", ref="16.131", notes=r"**N71 [B]** «16.131» and «16.133».")
P("P386", "choosing a root by continuation", r"""A quadratic has two roots; which one is "the" wave? Start where the answer is
known (here α → ∞: X → 1, the inviscid wave of C08), then lower α in small steps, each time picking the root nearest the
previous one. It extends the complex quadratic formula of Ch. 6 (P159).""", code=r"""
def roots(b, c): d = np.sqrt(b*b - 4*c + 0j); return (-b + d)/2, (-b - d)/2   # roots of X² + bX + c
prev = 1.0                                     # the known root at eps = 0
for eps in (0.0, 0.1, 0.2):                    # a family X² − (1 + eps(1+1j))X + 0 = 0
    r = roots(-(1 + eps*(1 + 1j)), 0.0); prev = min(r, key=lambda z: abs(z - prev)); print(prev)   # follow the root near 1
""")
remind("C11")
D("D21", r"The characteristic quadratic", ref="16.136")
worked("α = 8, ν̂ = 0.5", r"""
1. β̂ = e^{iπ/4} × 8 = 5.657 + 5.657i; g = β̂J₀(β̂)/J₁(β̂) = 6.19 − 5.62i (`jv`).
2. «16.136» has two roots; the one near 1 is X = 1.094 + 0.137i.
3. √X = 1.048 + 0.065i: c/c₀ = 1/1.048 = **0.954** — 5 % slower than Moens–Korteweg.
4. k₂/k₁ = 0.065/1.048 = **0.062**: the amplitude falls to e^{−2π × 0.062} = **0.68** in one wavelength.
5. In the aorta-like α = 17 the pulse would lose only about a fifth per wavelength; at α = 3 it loses four fifths.""")
code(r"""
print(ch16.elastic_tube_roots(8.0, 0.5))                    # both roots of (16.136), labelled by continuation (trap T7)
for al in (3.0, 8.0, 20.0):                                  # three Womersley numbers
    wv = ch16.elastic_tube_wave(al, c0=1.0, poisson=0.5)     # pressure root → c/c0, k2/k1, damping per wavelength
    print(f"α = {al:4.1f}: c/c0 = {wv['c_over_c0']:.4f}, k2/k1 = {wv['k2_over_k1']:.4f}, kept per wavelength = {wv['damping_per_wavelength']:.3f}")
wo = ch16.elastic_tube_wave(8.0, c0=1.0, poisson=0.5, branch="other")   # the other root (excluded)
print(f"other root at α = 8: c/c0 = {wo['c_over_c0']:.2f}, kept per wavelength = {wo['damping_per_wavelength']:.2e}")
""", explain=r"""
1. `elastic_tube_roots` solves the quadratic and labels the root that tends to X = 1 as α → ∞ ("pressure") by following it
   down in α (P386); the other root is a fast, heavily damped wall mode (trap T7).
2. `elastic_tube_wave` turns X into the speed ratio c/c₀ = 1/Re√X and the damping k₂/k₁ = Im√X/Re√X.""")
scratch(r"""
Xp = ch16.elastic_tube_roots(8.0, 0.5)["pressure"]          # the closed-form root
d_root = ch16.elastic_tube_determinant_numeric(Xp, 8.0, 0.5)    # mode ODEs integrated numerically, no Bessel closed form (ours)
d_one = ch16.elastic_tube_determinant_numeric(1.0, 8.0, 0.5)    # the same determinant at X = 1 (the inviscid guess)
assert abs(d_root) < 1e-8*abs(d_one)                          # vanishes at the root, not at X = 1
print(f"|det| at the root {abs(d_root):.1e}, at X = 1 {abs(d_one):.1e}")
""", explain=r"""An independent route: the mode equations of D18–D20 are integrated numerically from the axis (no Bessel
functions) and the 2 × 2 wall determinant is evaluated; it vanishes at the closed-form root of «16.136» and not at X = 1.""")
figure(r"""
al_ = np.logspace(0, np.log10(50), 200)                       # Womersley numbers 1 … 50
pw = ch16.elastic_tube_wave_curve(al_, poisson=0.5)           # pressure root along the whole curve (continuation)
ow = ch16.elastic_tube_wave_curve(al_, poisson=0.5, branch="other")   # the other root
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.6))
ax.semilogx(al_, pw["c_over_c0"], color=COLORS["blue"], lw=2.5, label="pressure root")
ax.semilogx(al_, ow["c_over_c0"], color=COLORS["muted"], lw=1, label="excluded: fast, heavily damped")
ax.axhline(1, ls="--", color=COLORS["muted"]); ax.set(ylim=(0, 6), xlabel="α (log)", ylabel="c / c₀ [–]", title="Viscosity slows the pulse")
bx.loglog(al_, pw["k2_over_k1"], color=COLORS["accent"], lw=2.5, label="k₂/k₁ (pressure root)")
bx.loglog(al_, 0.5/al_, "--", color=COLORS["muted"], label="∝ 1/α (guide)")
bx.set(xlabel="α (log)", ylabel="k₂ / k₁ [–]", title="…and damps it, both ∝ 1/α at large α")
ax.legend(fontsize=8); bx.legend(fontsize=8)
savefig(fig, "ch16", "c11_viscous_waves"); plt.show()
""", see=r"""Left: the wave speed relative to Moens–Korteweg's for the pressure root (bold blue) and the excluded root (faint);
right: the damping ratio k₂/k₁ of the pressure root with a 1/α guide.""",
    read=r"""At large α the pressure wave approaches c₀ and its damping falls like 1/α — the aorta carries nearly the inviscid
wave. Toward α ≈ 1–3 the wave slows markedly and loses most of its amplitude within a wavelength. The excluded root is far
faster and dies within a fraction of a wavelength — not the pulse.""",
    change=r"""…the wall were viscoelastic (N72)? Both curves shift toward more damping at every α.""")
plotly(r"""
nus = [0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]   # Poisson's ratio (9 positions)
al_p = np.logspace(0, np.log10(50), 120)                      # α grid
def roots_curves(nu_):                                        # both roots for one Poisson ratio
    p_ = ch16.elastic_tube_wave_curve(al_p, poisson=nu_); o_ = ch16.elastic_tube_wave_curve(al_p, poisson=nu_, branch="other")
    return {"c/c₀ pressure root": (al_p, p_["c_over_c0"]), "k₂/k₁ pressure root": (al_p, p_["k2_over_k1"]),
            "c/c₀ other root": (al_p, o_["c_over_c0"]), "k₂/k₁ other root": (al_p, o_["k2_over_k1"])}
fig = slider_figure(roots_curves, "ν̂", nus, xlabel="Womersley number α (log)", ylabel="ratio [–]",
                    title="Poisson's ratio barely moves the pulse", xrange=[0, np.log10(50)], yrange=[0, 2.5], active=8)
fig.update_xaxes(type="log"); fig.show()
""", explain=r"""**Figure F4.** Slider: Poisson's ratio of the wall, 0.3 to 0.5 (opens on 0.5, an incompressible wall). **What you
see.** Speed ratio and damping of both roots against α. **How to read it.** The pressure root's curves hardly move with
ν̂; the other root (cut at 2.5) changes more. **What would change if…** ν̂ = 0: the Poisson coupling disappears and the
wall acts only through its hoop stiffness.""")
animation(r"""
nfr = 24 if FAST else 40                                      # frames for one period
xw = np.linspace(0, 3*2*np.pi, 400)                           # three inviscid wavelengths (k0 = 1, so x in units of 1/k0)
wts = np.linspace(0, 2*np.pi, nfr, endpoint=False)            # phases ωt
waves = {}
for al in (3.0, 15.0):                                        # a small artery and an aorta-like case
    sX = np.sqrt(ch16.elastic_tube_roots(al, 0.5)["pressure"])   # k/k0 = √X (Re > 0)
    waves[al] = sX                                            # keep the complex wavenumber ratio for this α
fig, axs = plt.subplots(2, 1, figsize=(7.5, 3.6), sharex=True)
lines = []
for ax, al in zip(axs, (3.0, 15.0)):                          # top panel α = 3, bottom α = 15
    ax.plot(xw, np.exp(-waves[al].imag*xw), ":", color=COLORS["rose"], lw=0.8)   # the decaying envelope e^{−k2 x}
    (lv,) = ax.plot([], [], color=COLORS["blue"], lw=2, label=f"viscous, α = {al:g}")
    (li,) = ax.plot([], [], "--", color=COLORS["muted"], label="inviscid (C08)")
    lines.append((lv, li, al)); ax.set(ylim=(-1.1, 1.1), ylabel="p / p_in"); ax.legend(fontsize=7, loc="upper right")
axs[1].set_xlabel("distance k₀x [rad]")
def update(i):                                                # frame i
    out = []
    for lv, li, al in lines:                                  # viscous line, inviscid line and α of each panel
        lv.set_data(xw, np.real(np.exp(1j*(waves[al]*xw - wts[i]))))   # Re e^{i(kx − ωt)} (16.116)
        li.set_data(xw, np.cos(xw - wts[i])); out += [lv, li]
    return out
fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it (faster frames)
show_animation(animate(update, frames=nfr, fig=fig, interval=70), player="video")
""", explain=r"""**Figure A5.** **What you see.** One period of a pressure wave along three inviscid wavelengths for α = 3 (top) and
α = 15 (bottom), against the inviscid wave of the same frequency (grey dashed); the dotted line is the envelope
$e^{-k_2x}$. **How to read it.** At α = 3 the wave falls behind the inviscid one and fades within a wavelength; at α = 15 it
almost keeps pace and keeps most of its amplitude. **What would change if…** α = 50: the two curves would nearly coincide
over the three wavelengths.""")
note("N72 [B]", r"""A viscoelastic wall: replace E and ν̂ by «16.137». This is the **Kelvin–Voigt** model of a material: a spring
and a damper side by side, so the wall's stress is $E\varepsilon+E'\,d\varepsilon/dt$ — the elastic part plus a part
proportional to the rate of stretching. For a mode $e^{-i\omega t}$, $d/dt\to-i\omega$, which gives the complex modulus
$E-i\omega E'$ (the sign of the imaginary part belongs to the convention, trap T1).""")
code(r"""
for Er in (1.0, 1 - 0.1j):                                   # elastic, and a wall with 10 % viscous modulus
    wv = ch16.elastic_tube_wave(8.0, E_ratio=Er)             # (16.136) with E* = E·E_ratio (16.137)
    print(f"E*/E = {Er}: c/c0 = {wv['c_over_c0']:.4f}, kept per wavelength = {wv['damping_per_wavelength']:.3f}")   # more damping
""")
note("N73 [C]", r"""A steady stream adds to the wave speed downstream and subtracts upstream (Ch. 15's U ± c, `AC.wave_speeds`);
small in arteries.""")
whatif(r"""…the tube branched? A pulse meeting a fork is partly reflected — C13. But first: why are the branches the sizes they
are? C12.""")

# ---------------------------------------------------------------------------------------------------------------------
# C12 — Murray
# ---------------------------------------------------------------------------------------------------------------------
core("C12", r"Murray's law: «16.140» and «16.144»", "Why are blood vessels the width they are?")
problem(r"""
A wider pipe is easier to push water through — so why not make every vessel wide? Because blood is expensive: the body must
make it, warm it and carry it around. Murray asked which radius makes **the total cost — pumping plus upkeep — smallest**.
The answer is a rule you can test with a ruler on any cast of a vascular tree.""")
idea(r"""
 cost ↑  pumping ∝ a⁻⁴ ╲                ╱ upkeep ∝ a²
         ╲             ╱
           ╲_minimum__╱
         ──────────────────────────► radius a
""", r"""A falling power plus a rising power has one minimum (P349, Ch. 14).""")
note("N74 [B]", r"""The cost of a segment: pumping power Δp·Q = RQ² with R from «16.20», plus an upkeep power K per unit volume of
blood (W/m³ — the book's K, not the compliance; trap T2): «16.138».""")
code(r"""
print(ch16.murray_cost(5e-3, CO, murray["L"], murray["mu"], murray["b_metabolic"]), "W")   # (16.138) for a 5-mm, 10-cm segment carrying our CO
""")
remind("C12")
D("D22", r"Murray's optimal radius and the minimum cost", ref="16.140")
slip("S15", r"«16.140p» and the same bracket in (16.143)", r"«16.140»", r"""A units check catches it: μ/K has units s², so
$(\mu Q^2/K)^{1/6}$ is a length; the printed $(\mu KQ^2)^{1/6}$ is not.""")
code(r"""
mu_q, Q_q, K_q = UN.Q_(3.5e-3, "Pa*s"), UN.Q_(CO, "m**3/s"), UN.Q_(1000.0, "W/m**3")   # viscosity, flow, upkeep with units
print("corrected is a length:", UN.check_dimensions(((16*mu_q*Q_q**2/(np.pi**2*K_q)).to_base_units())**(1/6), "[length]"))   # True: μQ²/K is m⁶
print("printed is a length:  ", UN.check_dimensions(((16*mu_q*K_q*Q_q**2/np.pi**2).to_base_units())**(1/6), "[length]"))      # False
print("units of the printed a:", ((16*mu_q*K_q*Q_q**2/np.pi**2).to_base_units())**(1/6))   # not metres
""")
worked("an aorta-sized optimum", r"""
Q = 82.5 mL/s, μ = 3.5 mPa s, K = 1000 W/m³:
1. 16μ/(π²K) = 0.056/9870 = 5.67 × 10⁻⁶ s² (μ/K is (Pa s)/(Pa/s)); its sixth root 0.1335 s^{1/3}.
2. Q^{1/3} = (8.25 × 10⁻⁵)^{1/3} = 0.04355 m s^{−1/3}.
3. a = 0.1335 × 0.04355 = **5.81 mm** (ours).
4. Over 10 cm: pumping 5.31 mW, upkeep 10.62 mW — exactly 1 : 2.
5. The printed formula would give 58.1 mm.""")
code(r"""
a_star = ch16.murray_radius(CO, murray["mu"], murray["b_metabolic"])                 # (16.140, corrected) [m]
terms = ch16.murray_terms(a_star, CO, murray["L"], murray["mu"], murray["b_metabolic"])   # pumping and upkeep at the optimum [W]
E_min = ch16.murray_min_cost(CO, murray["L"], murray["mu"], murray["b_metabolic"])   # (16.141) [W]
a_bad = ch16.murray_radius(CO, murray["mu"], murray["b_metabolic"], printed_16_140=True)   # the printed (16.140)
print(f"a* = {a_star*1e3:.4f} mm; pumping {terms['pumping']*1e3:.3f} mW, upkeep {terms['upkeep']*1e3:.3f} mW (ratio {terms['pumping']/terms['upkeep']:.3f}); "
      f"min cost {E_min*1e3:.3f} mW; printed formula {a_bad*1e3:.2f} mm")
""")
scratch(r"""
a = np.linspace(1e-3, 20e-3, 20001)                            # candidate radii [m]
cost = 8*murray["mu"]*murray["L"]*CO**2/(np.pi*a**4) + murray["b_metabolic"]*np.pi*a**2*murray["L"]   # (16.138) by hand [W]
assert abs(a[np.argmin(cost)] - a_star) < 2e-6                 # the grid minimum is Murray's radius
assert ch16.murray_cost(a_bad, CO, murray["L"], murray["mu"], murray["b_metabolic"]) > 10*E_min   # the printed radius costs > 10×
print(f"grid minimum at {a[np.argmin(cost)]*1e3:.4f} mm")
""", explain=r"""A brute-force search over 20 001 radii finds the minimum of «16.138» exactly where «16.140» puts it; the printed
radius costs more than ten times the minimum.""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw=dict(width_ratios=[1.5, 1]))
aa = np.logspace(np.log10(2e-3), np.log10(15e-3), 200)         # radii [m]
tt_ = ch16.murray_terms(aa, CO, murray["L"], murray["mu"], murray["b_metabolic"])   # the two costs [W]
ax.loglog(aa*1e3, tt_["pumping"]*1e3, color=COLORS["orange"], label="pumping ∝ a⁻⁴")
ax.loglog(aa*1e3, tt_["upkeep"]*1e3, color=COLORS["teal"], label="upkeep ∝ a²")
ax.loglog(aa*1e3, tt_["total"]*1e3, color=COLORS["amber"], lw=2.5, label="total")
ax.axvline(a_star*1e3, color=COLORS["muted"], ls=":"); ax.text(a_star*1e3*1.03, 300, f"a* = {a_star*1e3:.2f} mm", fontsize=8)
ax.set(xlabel="radius a [mm] (log)", ylabel="power [mW] per 10 cm", ylim=(1, 1e3), title="The cheapest vessel: pumping = half the upkeep")
ax.legend(fontsize=8)
th = ch16.murray_symmetric_angle()                             # arccos 2^{-1/3} [rad]
for s in (1, -1):                                              # a symmetric Murray fork drawn to scale
    bx.plot([0, 2*np.cos(th)], [0, s*2*np.sin(th)], color=COLORS["teal"], lw=8*2**(-1/3), solid_capstyle="round")   # a daughter, line width ∝ radius
bx.plot([-2, 0], [0, 0], color=COLORS["teal"], lw=8, solid_capstyle="round")   # the parent vessel
bx.plot([0, 1.2], [0, 0], ":", color=COLORS["muted"])          # the parent's extension, from which the angles are measured
bx.text(0.9, 0.18, f"{np.degrees(th):.2f}°", fontsize=9); bx.set(xlim=(-2.2, 2.2), ylim=(-1.6, 1.6), aspect="equal", title="Symmetric fork: a₁ = a₂ = 2^(−1/3) a₀")
bx.axis("off")
savefig(fig, "ch16", "c12_murray"); plt.show()
for x_ in (0.7, 1.3):                                          # radii 30 % below and above the optimum
    extra = ch16.murray_cost(x_*a_star, CO, murray["L"], murray["mu"], murray["b_metabolic"])/E_min - 1   # (16.138) over (16.141)
    print(f"a = {x_:.1f} a*: {100*extra:+.0f} % above the minimum cost (formula ⅓x⁻⁴ + ⅔x² − 1 = {100*(x_**-4/3 + 2*x_**2/3 - 1):+.0f} %)")
""", see=r"""Left: pumping (orange), upkeep (teal) and total (amber) power of a 10-cm segment carrying our cardiac output, against
its radius, with the optimum marked. Right: a symmetric Murray fork drawn to scale, each daughter leaving at 37.47° from the
parent's extension.""",
    read=r"""The total cost has one minimum where the pumping curve sits at half the upkeep curve (step 8 of D22). It is
lopsided: relative to the optimum a segment of radius $x\,a^*$ costs $\tfrac13x^{-4}+\tfrac23x^2$ times the minimum, so a
vessel 30 % too wide costs 24 % more but one 30 % too narrow costs 71 % more (printed) — the steep $a^{-4}$ pumping side
punishes narrowing hardest. The printed formula's tenfold radius costs more than ten times the minimum.""",
    change=r"""…the flow were eight times larger? The optimum would move to twice the radius ($Q\propto a^3$) — slide it in F5.""")
plotly(r"""
Qs_ = list(np.round(np.logspace(np.log10(1e-6), np.log10(200e-6), 25), 9))   # flows 1 … 200 mL/s
aa_p = np.logspace(np.log10(0.5e-3), np.log10(20e-3), 200)    # radii [m]
def costs(Q_):                                                 # three costs for one flow
    tt = ch16.murray_terms(aa_p, Q_, murray["L"], murray["mu"], murray["b_metabolic"])
    return {"pumping": (aa_p*1e3, tt["pumping"]*1e3), "upkeep": (aa_p*1e3, tt["upkeep"]*1e3), "total": (aa_p*1e3, tt["total"]*1e3)}
fig = slider_figure(costs, "Q", Qs_, unit="m³/s", xlabel="radius a [mm] (log)", ylabel="power [mW] per 10 cm (log)",
                    title="The optimum slides as Q^(1/3)", xrange=[np.log10(0.5), np.log10(20)], yrange=[-3, 3],
                    active=int(np.argmin(np.abs(np.array(Qs_) - CO))))
fig.update_xaxes(type="log"); fig.update_yaxes(type="log"); fig.show()
""", explain=r"""**Figure F5.** Slider: the flow carried, 1 to 200 mL/s (opens near our 82.5 mL/s). **What you see.** The three costs
against radius. **How to read it.** The minimum of the total moves right by a factor 2 for every factor 8 in Q. **What would
change if…** blood were cheaper (smaller K): the upkeep line drops and every optimum widens as $K^{-1/6}$.""")
D("D23", r"The cube law", ref="16.144", notes=r"**N75 [B]** the cube law «16.144».")
code(r"""
a1, a2 = ch16.murray_daughters(a_star, 0.5)                  # an even split: each daughter 2^{-1/3} a0
print(f"{a1/a_star:.4f}, {a2/a_star:.4f}; check a0³ − a1³ − a2³ = {a_star**3 - a1**3 - a2**3:.1e}")   # 0.7937; cube law holds
""")
P("P387", "scipy.optimize.minimize for a point in the plane", r"""`minimize(f, x0, method='Nelder-Mead')` finds a minimum of a
function of several variables by moving a small triangle downhill; no derivatives needed. Here it finds the cheapest
branch point B. `minimize_scalar` (P170) did the one-variable case.""", code=r"""
from scipy.optimize import minimize            # several-variable minimiser
f = lambda B: (B[0] - 1)**2 + (B[1] + 2)**2    # a bowl centred at (1, −2)
print(minimize(f, [0, 0], method="Nelder-Mead").x)   # ≈ [1, −2]
""")
P("P388", "the weighted Fermat point", r"""Move a point B by a small step ds in a direction m: its distance to a fixed point P
changes by −(unit vector from B to P)·m ds (a directional derivative, P75). So the weighted sum Σwᵢ|Pᵢ − B| is least where
the weighted unit vectors wᵢeᵢ add to zero — they close a triangle of "forces", whose angles the cosine rule (P340)
gives.""", code=r"""
w = np.array([1.0, 0.63, 0.63])                 # weights a_i² of a symmetric Murray fork
print(np.degrees(np.arccos((w[0]**2 + w[1]**2 - w[2]**2)/(2*w[0]*w[1]))))   # 37.47°
""")
D("D24", r"The branching angles", ref="16.147", notes=r"**N76 [B]** the angles «16.146» and «16.147».")
code(r"""
a0_, a1_, a2_ = 1.0, 0.7**(1/3), 0.3**(1/3)                   # a 70/30 split by the cube law (radii relative to the parent)
ang = ch16.murray_angles(a0_, a1_, a2_)                       # (16.146) [rad]
bp = ch16.bifurcation_point(np.array([0.0, 0.0]), np.array([3.0, 1.2]), np.array([2.6, -1.5]), a0_, a1_, a2_)   # numerical optimum (ours)
print(f"symmetric: {np.degrees(ch16.murray_symmetric_angle()):.3f}°; 70/30: θ = {np.degrees(ang['theta']):.2f}°, φ = {np.degrees(ang['phi']):.2f}°; "
      f"numeric θ = {np.degrees(bp['theta']):.2f}°, (16.147) θ = {np.degrees(ch16.murray_theta(a0_, a1_)):.2f}°")
""")
explainer("murray_branching", "Why isn't a wider vessel always better?",
          r"""you drag the radius past the minimum and watch the two cost bars trade, then drag the branch point and see the cost
bottom out exactly where three weighted pulls close a triangle — the balance is found by moving, which a picture cannot
do.""",
          ["Put the radius cursor at the minimum: pumping is half the upkeep.", "Pick the 70/30 split: the thin daughter turns more.",
           "Switch to the printed «16.140p» and raise K."])
whatif(r"""…a pressure pulse arrived at this fork? Its daughters, Murray-sized, are almost matched to the parent — C13 says what
"matched" means.""")

# ---------------------------------------------------------------------------------------------------------------------
# C13 — junctions
# ---------------------------------------------------------------------------------------------------------------------
core("C13", r"Echoes at a fork: «16.153»", "When a pulse reaches a branch, how much comes back?")
problem(r"""
Send a wave along a rope tied to a thinner rope: part of it bounces back. Tie it to a thicker one and the echo comes back
upside down. At an arterial fork the "thickness" that matters is not the radius alone but **how easily each tube accepts
flow for a given pressure — its admittance**.""")
idea(r"""
incoming pulse ──►  ┌ junction ┐ ──► transmitted into each daughter
            ◄── echo└──────────┘
daughters together accept less flow than the parent  ──► like a closed end ──► upright echo   (R > 0)
daughters accept exactly as much                       ──► matched          ──► no echo        (R = 0)
daughters accept more                                  ──► like an open end ──► inverted echo  (R < 0)
""", r"""**Only the mismatch of admittances matters.** A tube's admittance $Y=A/(\rho c)$ (cross-section area over density
times wave speed) is the flow it accepts per unit pressure in a travelling wave; the echo's share is «16.153».""")
note("N77 [B]", r"""At the junction the pressure is single-valued and flux is conserved: «16.148», «16.149» (the book's minus with
a negative $Q_R$ counts the sign twice — slip S16; we keep signs once, D25). **N78 [B]** Multiply «16.68» by A/(ρc):
«16.151», admittance Y = A/(ρc), impedance Z = ρc/A — the electrical analogy (P381).""")
code(r"""
A_art = np.pi*artery["a"]**2                                  # parent area [m²]
Y1 = ch16.tube_admittance(A_art, blood["rho"], c0)            # Y = A/(ρc) [m⁴ s/kg]
print(f"Y = {Y1:.3e}, Z = {ch16.tube_impedance(A_art, blood['rho'], c0):.3e}, Y·Z = {Y1*ch16.tube_impedance(A_art, blood['rho'], c0):.1f}")
""")
D("D25", r"Reflection and transmission at a junction", ref="16.153")
D("D26", r"Energy at a junction", ref="16.163", notes=r"**N80 [B]** the energy shares «16.161», «16.162», «16.163».")
worked("a symmetric fork", r"""
Daughters of radius 0.8a, same wave speed:
1. Y_j/Y₁ = (A_j/A₁)(c₁/c_j) = 0.64 each; ΣY/Y₁ = 1.28.
2. R = (1 − 1.28)/(1 + 1.28) = **−0.123**: a small inverted echo.
3. T = 1 + R = 0.877.
4. Energy: R² = 0.015 reflected, 1.28 × 0.877² = 0.985 transmitted; sum 1.
5. A Murray fork (daughters $2^{-1/3}a$): ΣY/Y₁ = 1.26, R = **−0.115** — nearly matched.""")
code(r"""
Rj, Tj = ch16.junction_coefficients(1.0, [0.64, 0.64])        # (16.153) with admittances relative to the parent
en = ch16.junction_energy(1.0, [0.64, 0.64])                  # reflected R², transmitted (ΣY/Y1)T², total
geo = ch16.junction_from_geometry(artery["a"], [0.8*artery["a"]]*2, c0, [c0, c0], blood["rho"])   # the same from radii and speeds
print(f"R = {Rj:.4f}, T = {Tj:.4f}; energy {en['reflected']:.4f} + {en['transmitted']:.4f} = {en['total']:.4f}; from geometry R = {geo['R']:.4f}")
print("Murray fork:", round(ch16.junction_coefficients(1.0, [2**(-2/3)]*2)[0], 4))   # Y ∝ a² at equal c: 2 × 2^{-2/3} = 1.26
""")
scratch(r"""
Y1_, Yd = 1.0, [0.64, 0.64]                                    # parent and daughter admittances
sol = np.linalg.solve([[1.0, -1.0], [Y1_, sum(Yd)]], [-1.0, Y1_])   # p_I + p_R = p_T and Y1(p_I − p_R) = ΣY p_T, per unit p_I
assert np.allclose(sol, [Rj, Tj])                              # (p_R, p_T) = (R, T)
rng = np.random.default_rng(0)                                 # random admittances (reproducible)
for _ in range(1000):                                          # 1000 random forks with 1–4 daughters
    y1, ys = rng.uniform(0.1, 5), list(rng.uniform(0.01, 5, rng.integers(1, 5)))   # a random parent admittance and 1–4 random daughter admittances
    r_, t_ = ch16.junction_coefficients(y1, ys)
    assert abs(r_**2 + sum(ys)/y1*t_**2 - 1) < 1e-13          # the energy identity of D26
print("linear solve = (16.153); energy identity holds for 1000 random forks")
""", explain=r"""`np.linalg.solve` solves the two junction conditions as a 2 × 2 linear system and returns R and T; then 1000 random
forks (`np.random.default_rng(0)`: a reproducible random generator) satisfy $\mathrm R^2+(\sum Y_j/Y_1)\mathrm T^2=1$ to
10⁻¹³.""")
figure(r"""
s_ = np.linspace(0, 4, 400)                                    # ΣY_j/Y1
fig, ax = plt.subplots(figsize=(7, 3.6))
ax.plot(s_, (1 - s_)/(1 + s_), color=COLORS["rose"], lw=2, label="R (reflected pressure)")
ax.plot(s_, 2/(1 + s_), color=COLORS["blue"], lw=2, label="T (transmitted pressure)")
ax.axhline(0, color=COLORS["muted"], lw=0.6); ax.axvline(1, color=COLORS["muted"], ls=":")
ax.axvspan(0, 1, color=COLORS["orange"], alpha=0.06); ax.axvspan(1, 4, color=COLORS["teal"], alpha=0.06)
ax.text(0.5, 1.75, "closed-end-like\nR > 0", ha="center", fontsize=8); ax.text(2.5, 1.3, "open-end-like, R < 0", ha="center", fontsize=8)
ax.text(1.03, -0.85, "matched", fontsize=8)
ax.plot(2**(1/3), (1 - 2**(1/3))/(1 + 2**(1/3)), "o", color=COLORS["accent"]); ax.annotate("Murray fork", (2**(1/3), -0.115), xytext=(1.6, -0.5), arrowprops=dict(arrowstyle="->"), fontsize=8)   # purple dot: a Murray fork, ΣY/Y₁ = 2^{1/3}
ax.set(xlabel="ΣY_j / Y₁ [–]", ylabel="ratio to the incident pressure [–]", ylim=(-1.05, 2.05), title="Reflection is set by the admittance mismatch")
ax.legend(fontsize=8, loc="upper right")
savefig(fig, "ch16", "c13_reflection"); plt.show()
""", see=r"""R (rose) and T (blue) of «16.153» against the daughters' total admittance relative to the parent's, with the
closed-end-like and open-end-like regimes shaded and our Murray fork marked.""",
    read=r"""At ΣY_j = Y₁ there is no echo. Daughters that accept less flow than the parent (left) reflect an upright echo and
transmit a pressure up to twice the incident one; daughters that accept more (right) reflect an inverted echo. The Murray
fork sits close to matched.""",
    change=r"""…the daughters were stiffer than the parent (larger c_j)? Their admittance A/(ρc) falls, the point moves left and
the echo becomes upright — one reason pressure pulses grow toward the periphery.""")
note("N79 [B]", r"""In the parent the pressure and flux are «16.159» and «16.160»: once the echo returns, p and Q no longer share a
shape. (The book writes the reflected shape as f; it equals f here only because R does not depend on frequency — slip
S17.)""")
code(r"""
fpulse = lambda s: np.exp(-(s/0.05)**2)                       # time shape of the incident pulse (width 0.05 s)
ts = np.linspace(-0.1, 0.4, 6)                                # instants [s]
wf = [ch16.junction_waveforms(-0.5, tt, fpulse, 5.0, Rj) for tt in ts]   # (16.159)-(16.160) at x = −0.5 m in the parent
print("p/P_I:", np.round([w_["p"] for w_ in wf], 3)); print("Q/(Y1 P_I):", np.round([w_["Q"] for w_ in wf], 3))
""", explain=r"""**What does this show?** At the probe the incident pulse passes first (p and Q equal), then the inverted echo
returns: p dips by R while Q rises by −R — the two traces part ways.""")
note("N81 [B]", r"""Three regimes (trap T18):

| Daughters | R | Name |
|---|---|---|
| ΣY_j < Y₁ | R > 0 | closed-end-like — the pressure rises at the fork (a closed end, ΣY = 0, doubles it: R = +1) |
| ΣY_j = Y₁ | 0 | matched |
| ΣY_j > Y₁ | R < 0 | open-end-like |

Reflections from the leg arteries help make the aortic pressure peak; repeated forks make a tree.""")
animation(r"""
nfr = 20 if FAST else 30                                      # frames (the reader steps through them)
net = ch16.network_pulse_run(1.0, [1.0, 1.0], 5.0, [5.0, 5.0], 1.0, [0.64, 0.64], width=0.03, n_frames=nfr)   # grid run (ours)
xp_, xd_ = net["x"]["parent"], net["x"]["daughters"][0]       # parent x < 0, a daughter x > 0 [m]
fig, (ax, bx) = plt.subplots(1, 2, figsize=(8.5, 3.0), gridspec_kw=dict(width_ratios=[1.6, 1]))
(lp,) = ax.plot([], [], color=COLORS["orange"], lw=2, label="parent"); (ld,) = ax.plot([], [], color=COLORS["teal"], lw=2, label="each daughter")
ax.axvline(0, color=COLORS["muted"], ls=":"); ax.axvline(net["x_probe"], color=COLORS["accent"], ls="--", lw=0.8)
ax.set(xlim=(-1, 1), ylim=(-0.3, 1.1), xlabel="x [m] (fork at 0)", ylabel="p / P_I", title="A pulse meets a fork (grid run, ours)"); ax.legend(fontsize=7)
(pp,) = bx.plot([], [], color=COLORS["orange"], label="p at the probe"); (qq,) = bx.plot([], [], color=COLORS["blue"], label="Q/Y₁ at the probe")
bx.set(xlim=(0, net["probe_t"][-1]), ylim=(-0.3, 1.3), xlabel="t [s]", title="probe (dashed line)"); bx.legend(fontsize=7)
def update(i):                                                # frame i
    lp.set_data(xp_, net["p_parent"][i]); ld.set_data(xd_, net["p_daughters"][0][i])
    k = np.searchsorted(net["probe_t"], net["t"][i])          # how many probe samples have been recorded by this frame (P395)
    pp.set_data(net["probe_t"][:k], net["probe_p"][:k]); qq.set_data(net["probe_t"][:k], net["probe_Q"][:k]/1.0)
    return lp, ld, pp, qq
fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it (faster frames)
show_animation(animate(update, frames=nfr, fig=fig, interval=150), player="frames")
print(f"grid R = {net['R_measured']:.4f}, theory {net['R_theory']:.4f}")
""", explain=r"""**Figure A4** (frame player — step through it). **What you see.** A pressure pulse runs down the parent (orange) into two
equal daughters (teal); at the fork a small inverted echo returns while 88 % of the pressure goes on in each daughter; the
right panel records pressure and flow at the probe. **How to read it.** Stop on the frame where the echo reaches the probe:
the pressure dips below zero while the flow rises above zero — the echo carries flow backward with a negative pressure.
**What would change if…** the daughters were narrower (ΣY < Y₁): the echo would come back upright.""")
explainer("junction_reflection", "When does a fork send the pulse back?",
          r"""changing the daughters' radius or stiffness you watch the echo shrink, vanish and flip while the energy bars always add
to the incident pulse — the sign change of the echo is a motion through zero, best seen live.""",
          ["Find the radius ratio with no echo.", "Make the daughters narrow: the echo comes back upright.",
           "Watch the probe: p and Q part ways when the echo arrives."])
whatif(r"""…the vessel curved instead of branching? A bend does not reflect a long pulse much, but it reshapes the steady flow inside
it — C14.""")

# ---------------------------------------------------------------------------------------------------------------------
# C14 — Dean
# ---------------------------------------------------------------------------------------------------------------------
core("C14", r"A bend makes two vortices: Dean's expansion «16.177»", "What does a bend do to the flow inside a pipe?")
note("N82 [C]", r"""The scene: the aortic arch turns the blood through half a circle (and out of plane). Curvature distorts the
profile, sets up secondary motions, raises dissipation and delays turbulence; Ch. 11 (N69) named Dean flow among the
centrifugal instabilities.""")
problem(r"""
Stir tea and the leaves gather in the middle of the bottom — a secondary flow. In a bent pipe the fast fluid at the centre
is flung toward the outside of the bend harder than the slow fluid near the walls; it piles up at the outer wall and returns
along the walls: **two counter-rotating vortices, and the fastest axial flow sits off-centre, toward the outer wall**.""")
idea(r"""
centre of the pipe:  fast fluid ──► strong centrifugal push ──► flung toward the outer wall
near the walls:      slow fluid ──► weak push              ──► squeezed back inward along the walls
result: two counter-rotating cells (Dean vortices), fastest axial flow shifted toward the outer wall
""", r"""**A bend sorts the fluid by speed**: the centrifugal push grows with the square of the speed, so the fast core wins
and the slow wall fluid gives way. How strong the secondary flow is depends on the speed and on how tight the bend is,
combined in one number defined below; the book solves the flow as a series in that number, «16.177».""")
P("P389", "toroidal coordinates of a curved pipe", r"""Inside a pipe bent into a circle of radius R: r and α locate a point in the
cross-section (α = 0 toward the outside of the bend here), s runs along the centreline. The axial speed w is carried round
the bend, so it feels a centrifugal acceleration w²/R pointing outward — that is the term $-w^2\cos\alpha$ of the radial
momentum equation of the curved pipe, written out in full in N86 below (trap T2: this α is an angle).""", code=r"""
w, R = 0.3, 0.05                               # axial speed [m/s], bend radius [m]
print(w**2/R)                                  # 1.8 m/s²: a fifth of g
""")
note("N83 [B]", r"""The secondary flow forms Dean vortices; for a developed flow the axial maximum skews to the outer wall (for a
potential-vortex entrance flow, to the inner wall); heart–lung machines use coiled tubes for mixing.

```
            inner wall (α = π)
          ╭───────────────╮
         │  ↖    ⟲    ↙   │      upper cell turns one way,
         │ ←  ←  ←  ←  ←  │      fluid returns along the walls
  axis ──│ →  →  →  → ●→  │──    the core is flung outward (●: fastest axial flow)
         │ ←  ←  ←  ←  ←  │
         │  ↙    ⟳    ↖   │      lower cell the other way
          ╰───────────────╯
            outer wall (α = 0) is on the right of the bend
```

**N84 [B]** Curvature ratio «16.164» — not the cell-free layer δ of C05 (trap T2).""")
code(r"""
print(ch16.curvature_ratio(12.5e-3, 50e-3))                  # (16.164): δ = a/R for our illustrative arch: 0.25
""")
note("N85 [B]", r"""The book rescales so that centrifugal, inertial and viscous terms are the same size (√δ scaling) and splits the
pressure «16.165» (stated). **N86 [B]** The resulting equations, stated (u, v the secondary velocities along r and α, w the
axial one, $p_1$ the cross-section pressure, all rescaled) — continuity, radial and angular momentum, axial momentum and
the wall conditions:

$$ """ + EQS("16.166", "16.167", "16.168", "16.169", "16.170") + r""" $$

The one term that drives everything is the centrifugal push $-w^2\cos\alpha$ in the radial equation (and its partner
$w^2\sin\alpha$ in the angular one): without the axial flow w there is no secondary flow. **N88 [B]** A stream function for the secondary flow:
«16.172» (Ch. 4's stream function), and the axial equation «16.173» with $p_0=-Gs$. **N89 [B]** The vorticity equation
«16.174» with «16.175» and «16.176». Read its right side in words: the axial speed's own gradient drives the secondary
vorticity.""")
slip("S34", r"«16.174p»", r"«16.174»", r"""Found when we checked the algebra: eliminating the pressure $p_1$ from the radial and
angular momentum equations of N86 — the ones carrying $-w^2\cos\alpha$ (16.167) and $+w^2\sin\alpha$ (16.168) — with
«16.172» gives **+** in front of the advection term. ψ₀, w₁ and the K² coefficient 1541/50400
are unchanged (the term first enters at order κ²); the K⁴ coefficient of «16.178» becomes
6471982981/542442700800 = 0.0119312 instead of the 0.0121130 the printed sign gives. The book's own 0.0120 sits between
them; we show 0.0119312 beside it (`ch16.DEAN_C4_DERIVED`).""")
note("N87 [B]", r"""The Dean number «16.171». **Which Dean number?** (trap T8)

| Name | Definition | Relation |
|---|---|---|
| κ | «16.171» | — |
| K | used in «16.178» | $K=2\kappa^2$ |
| D | «16.179» | $K=D^2/16$ |
| De (literature) | $Re\sqrt\delta$ with the diameter Re | = κ |

With the chapter's diameter Re, (16.171) reads $\kappa=\sqrt\delta\,Re$ (the book writes $\sqrt\delta\,2Re$ with a radius
Re — slip S18).""")
code(r"""
kap = ch16.dean_number(0.25, 5e-3, 12.5e-3, 3.3e-6, convention="kappa")   # (16.171) at W̄ = 5 mm/s in our arch
print(f"kappa = {kap:.3f}, K = {ch16.dean_number(0.25, 5e-3, 12.5e-3, 3.3e-6, convention='K'):.1f}, D = {ch16.dean_number(0.25, 5e-3, 12.5e-3, 3.3e-6, convention='D'):.2f}")
print("K = 576 is D =", ch16.dean_convert(576.0, from_="K", to="D"))   # 96
""")
P("P391", "a regular perturbation series", r"""When a small parameter ε appears in an equation, write the answer as
$f=f_0+\varepsilon f_1+\varepsilon^2f_2+\dots$, put it in, and collect equal powers of ε: each power gives a simpler equation,
solved in turn. Ch. 2's orders of smallness (P68) only dropped terms; here we keep and solve them.""", code=r"""
eps, x = sp.symbols('epsilon x')                # a small parameter and the unknown
root = sp.solve(x**2 + eps*x - 1, x)[1]          # the root that is x = 1 at eps = 0 (exact, with a square root)
print(sp.series(root, eps, 0, 3))                 # 1 − eps/2 + eps²/8 + …: the terms a perturbation series finds one order at a time
""")
P("P390", "the biharmonic operator in polar coordinates", r"""$\nabla^4=\nabla^2\nabla^2$. On $r^n\sin\alpha$ the polar
Laplacian gives $(n^2-1)r^{n-2}\sin\alpha$, so ∇⁴ lowers the power by four; a polynomial forcing is matched term by term
(undetermined coefficients), and the terms ∇⁴ kills (r, r³) carry the boundary conditions. Ch. 8's Stokes operator E² (P199)
is a different operator.""", code=r"""
r, a = sp.symbols('r alpha', positive=True)     # radius and angle
lap = lambda F: sp.diff(F, r, 2) + sp.diff(F, r)/r + sp.diff(F, a, 2)/r**2   # polar Laplacian
print(sp.simplify(lap(lap(r**5*sp.sin(a)))))    # 192*r*sin(alpha)
""")
md(r"""**Integration by parts** (used in step 13 below; Ch. 9): $\int_0^1 f\,g'\,dr=[fg]_0^1-\int_0^1 f'g\,dr$ — it moves a
derivative from one factor to the other.""")
D("D27", r"Dean's leading orders and the flux coefficient", ref="16.178",
  notes=r"""**N90 [B]** ψ₀ and w₁ (ours; **w₁ with the plus sign**, so the fast core sits toward the outer wall) and **N91 [B]** the flux
«16.178» with its K² coefficient reproduced exactly as 1541/50400; its K⁴ coefficient is printed 0.0120, and 0.0119312 with
the advection sign of «16.174» corrected (slip S34).""")
worked("the flux loss at K = 576", r"""
1. K/576 = 1.
2. Q/(πa²W̄) = 1 − 0.0306 + 0.0120 = **0.9814** (with the exact 1541/50400 and the printed 0.0120: 0.98142; with 0.0119312
   from slip S34: 0.98136).
3. The resistance ratio «16.182»: 1 + 0.0306 − 0.0110 = 1.0196.
4. In D: K = D²/16 gives D = 96, and «16.181» gives the same numbers.
5. For our arch (a = 12.5 mm, δ = 0.25, ν = 3.3 × 10⁻⁶ m²/s) K = 576 needs a mean speed of only W̄ = κν/(2a√δ) = 16.97 ×
   3.3 × 10⁻⁶/(2 × 0.0125 × 0.5) ≈ **4.5 mm/s** — the series is a slow-flow theory; real arch flows are far beyond it.""")
code(r"""
r_, al_ = np.array([0.5]), np.array([0.0])                   # a point halfway out on the outer side
print(ch16.dean_secondary_streamfunction(r_, np.pi/2), ch16.dean_axial_correction(r_, al_))   # ψ0 at α = π/2; w1 = +0.00093 (plus sign)
print("w at r = 0.5, α = 0 and π, K = 288:", ch16.dean_axial_velocity(np.array([0.5, 0.5]), np.array([0.0, np.pi]), 288.0))
Ks_ = np.array([100.0, 288.0, 576.0])                        # three Dean numbers K
print("flux ratio (book 0.0120):  ", [f"{v:.6f}" for v in ch16.dean_flux_ratio(Ks_, exact=True)])   # 1541/50400 and the printed K⁴ coefficient
print("flux ratio (S34 0.0119312):", [f"{v:.6f}" for v in ch16.dean_flux_ratio(Ks_, exact=True, exact_c4=True)], float(ch16.DEAN_C4_DERIVED))
print("lambda(576) =", ch16.curved_resistance_ratio(576.0), "| reliability (ours):", ch16.dean_series_reliability(576.0, threshold=0.25))
""", explain=r"""
1. `dean_secondary_streamfunction` and `dean_axial_correction` are ψ₀ and w₁ of D27 (ours); w₁ is positive on the outer side
   (α = 0) — the plus sign.
2. `dean_axial_velocity` = $w_0+\kappa^2w_1$ with κ² = K/2: at K = 288 the outer point is faster than the inner one.
3. `dean_flux_ratio` evaluates «16.178» with the exact 1541/50400; `exact_c4=True` uses the corrected K⁴ coefficient.
4. `dean_series_reliability` is our own criterion (the K⁴ term over the K² term, threshold ¼ passed explicitly): K = 576 is
   past it.""")
scratch(r"""
h_ = 1e-3                                                      # finite-difference step
rr, aa = 0.4, 0.7                                              # a test point (r, α)
lapf = lambda F, r, a: ((F(r + h_, a) - 2*F(r, a) + F(r - h_, a))/h_**2 + (F(r + h_, a) - F(r - h_, a))/(2*h_*r)   # ∂²F/∂r² + (1/r)∂F/∂r by centred differences …
                        + (F(r, a + h_) - 2*F(r, a) + F(r, a - h_))/(h_**2*r**2))   # … + (1/r²)∂²F/∂α²: the polar Laplacian (16.175)
psi0 = lambda r, a: ch16.dean_secondary_streamfunction(r, a)   # ψ0(r, α) from the module (D27 step 8)
lap_psi = lambda r, a: lapf(psi0, r, a)                        # ∇²ψ0
bih = lapf(lap_psi, rr, aa)                                    # ∇⁴ψ0 by nested differences
assert abs(bih - 8*rr*(1 - rr**2)*np.sin(aa)) < 1e-3           # D27 step 4: ∇⁴ψ0 = 8r(1 − r²) sin α
f_ = rr*(1 - rr**2)**2*(4 - rr**2)/144                         # f(r) of step 8
assert abs(lapf(ch16.dean_axial_correction, rr, aa) + 2*f_*np.cos(aa)) < 1e-6   # step 10: ∇²w1 = −2f cos α
R_, A_ = np.meshgrid(np.linspace(0, 1, 401), np.linspace(0, 2*np.pi, 401))    # a polar grid over the cross-section
w2d = ch16.dean_axial_velocity(R_, A_, 288.0)                  # w0 + κ² w1 at K = 288
flux = np.trapezoid(np.trapezoid(w2d*R_, R_[0], axis=1), A_[:, 0])/np.pi   # mean axial speed over the disc
assert abs(flux - 1) < 1e-4                                    # w1 carries no flux (cos α averages to 0): the loss is O(K²) and comes from w̄2 (step 12)
print(f"∇⁴ψ0 = {bih:.6f} vs {8*rr*(1 - rr**2)*np.sin(aa):.6f}; mean of w0 + κ²w1 = {flux:.6f}")
""", explain=r"""Finite differences of the polar Laplacian «16.175» confirm the two equations D27 solves (∇⁴ψ₀ = 8r(1 − r²) sin α and
∇²w₁ = −2f cos α), and a 2-D trapezoid over the disc shows that w₀ + κ²w₁ still carries exactly the straight-pipe flux: the
flux loss is a K⁴-order effect of the mean correction w̄₂ (steps 12–14).""")
figure(r"""
nr_, na_ = (40, 80) if FAST else (60, 120)                      # polar grid sizes
Rg, Ag = np.meshgrid(np.linspace(0, 1, nr_), np.linspace(0, 2*np.pi, na_))   # every (r, α) pair of the polar grid (np.meshgrid, C08)
Xg, Yg = Rg*np.cos(Ag), Rg*np.sin(Ag)                          # outer wall (α = 0) on the right
Wg = ch16.dean_axial_velocity(Rg, Ag, 576.0)                   # axial speed w0 + κ²w1 at K = 576 (units of W̄)
Pg = ch16.dean_secondary_streamfunction(Rg, Ag)                 # ψ0
fig, ax = plt.subplots(figsize=(6.0, 4.8))
cf = ax.contourf(Xg, Yg, Wg, levels=20, cmap="Blues"); fig.colorbar(cf, label="axial speed w / W̄ [–]")
ax.contour(Xg, Yg, Pg, levels=np.linspace(-0.0075, 0.0075, 12), colors=COLORS["accent"], linewidths=0.9)
k_ = np.unravel_index(np.argmax(Wg), Wg.shape)                  # (row, column) of the largest w: np.argmax counts in the flattened array, np.unravel_index turns that count back into 2-D indices
ax.plot(Xg[k_], Yg[k_], "o", color=COLORS["orange"])           # orange dot: the fastest point of the cross-section
ax.annotate(f"max at {Xg[k_]:.2f} a", (Xg[k_], Yg[k_]), xytext=(0.1, -0.45), fontsize=8, arrowprops=dict(arrowstyle="->"))
ax.add_patch(plt.Circle((0, 0), 1, fill=False, color=COLORS["ink"])); ax.text(1.02, 0, "outer\nwall", fontsize=8, va="center")
ax.set(aspect="equal", xlim=(-1.15, 1.35), ylim=(-1.1, 1.1), xlabel="x / a (toward the outside of the bend →)", ylabel="y / a",
       title="Two vortices push the core outward\n(K = 576)")
savefig(fig, "ch16", "c14_dean"); plt.show()
""", see=r"""The pipe's cross-section with the outer wall on the right: the axial speed $w_0+\kappa^2w_1$ at K = 576 (blue shading)
and the streamlines of the secondary flow ψ₀ (purple), two counter-rotating cells; the orange dot is the fastest point. To
find that point the code uses `np.unravel_index`: `np.argmax` of a 2-D array returns one count into the array read row by
row, and `np.unravel_index(count, shape)` converts it back into a (row, column) pair — e.g. count 5 in a 2 × 3 array is
(1, 2).""",
    read=r"""The secondary flow crosses the middle toward the outer wall and returns along the top and bottom walls. It carries
fast fluid outward, so the axial maximum shifts toward the outer wall (here by about 0.19 a) — the plus sign of w₁.""",
    change=r"""…K doubled? The shift doubles (it is proportional to κ² = K/2) while the vortex pattern keeps its shape — try F6.""")
note("N92 [B]", r"""With D: «16.181» — a change of variable, K = D²/16. In D's definition «16.179» the book takes «16.180»;
with that sign D would be negative, so Ĝ must be read as the magnitude $8\mu\bar W/a^2$ of the gradient (slip S19). **N93 [B]** The
resistance ratio «16.182» is the reciprocal of «16.178» (0.0306² − 0.0120 = −0.0111 — checked below; the text swaps
"curved" and "straight" — slip S20). At large Dean numbers numerical solutions show thin wall layers, a wall shear growing
with D, and even four-vortex states — named only.""")
code(r"""
for K_ in (100.0, 300.0, 576.0):                              # three Dean numbers K
    print(K_, f"{1/ch16.dean_flux_ratio(K_) - ch16.curved_resistance_ratio(K_):.1e}")   # small: (16.182) is the reciprocal of (16.178) to the order kept
print(f"0.0306² − 0.0120 = {0.0306**2 - 0.0120:.4f}")         # the K⁴ coefficient of the reciprocal
""")
plotly(r"""
K_vals = list(np.linspace(0, 600, 25))                         # Dean numbers K (25 positions)
xs_ = np.linspace(-1, 1, 200)                                  # along the plane of symmetry (200 points: the axis itself is skipped): x = r cos α, α = 0 (outer) or π (inner)
def cuts(K_):                                                  # axial speed along two diameters, and the secondary flow
    r_ = np.abs(xs_); al_ = np.where(xs_ >= 0, 0.0, np.pi)   # x > 0 is the outer half (α = 0), x < 0 the inner half (α = π)
    w_sym = ch16.dean_axial_velocity(r_, al_, K_)              # across the bend (outer wall at +1)
    w_perp = ch16.dean_axial_velocity(r_, np.where(xs_ >= 0, np.pi/2, 3*np.pi/2), K_)   # perpendicular diameter
    u_, v_ = ch16.dean_secondary_velocity(r_, al_, K_)         # secondary velocity components on the symmetry line
    return {"w across the bend (outer → right)": (xs_, w_sym), "w on the perpendicular diameter": (xs_, w_perp),
            "secondary speed ×10 on the symmetry line": (xs_, 10*u_*np.cos(al_))}   # x-component: toward the outer wall
fig = slider_figure(cuts, "K", K_vals, xlabel="position across the pipe x / a [–]", ylabel="speed / W̄ [–]",
                    title="The skew grows with K; the pattern keeps its shape", yrange=[-1.5, 2.6], active=24)
fig.show()
""", explain=r"""**Figure F6.** Slider: Dean's K from 0 to 600 (opens on 600). Instead of contour maps, three cuts that carry the same
information: the axial speed across the bend (outer wall at +1), the axial speed on the perpendicular diameter, and ten
times the secondary speed along the plane of symmetry. **How to read it.** At K = 0 both axial cuts are the same parabola
and there is no secondary flow; as K grows the cross-bend cut tilts toward the outer wall, linearly in K, while the
secondary flow grows as √K. **What would change if…** the series were pushed beyond about K = 460: our reliability
criterion says the next term is no longer small.""")
explainer("dean_curved_pipe", "What does a bend do to the flow?",
          r"""the secondary flow is three-dimensional: orbit the bend to see the tracers corkscrew, drag K to grow the skew, and watch
the flux cursor approach the edge of the series — a helix is hard to read from any single picture.""",
          ["Set K = 0 then 288: the bullseye slides outward.", "Orbit the 3-D view: helices.",
           "Switch the chip to D: 96 at K = 576."])
whatif(r"""…the vessel's wall could not hold its shape at all — a vein above the heart, an airway in a cough? Then the tube's area
depends on the flow: C15.""")

# ---------------------------------------------------------------------------------------------------------------------
# C15 — collapsible tubes
# ---------------------------------------------------------------------------------------------------------------------
core("C15", r"Floppy tubes choke at S = 1: «16.204»", "Why can't you suck harder through a soft straw?")
problem(r"""
Suck hard on a soft drinking straw: past a point, sucking harder brings no more drink — the straw pinches shut near your
mouth. Veins above the heart, airways in a forced cough and an artery under a blood-pressure cuff do the same. **A tube whose
wall gives way can cap its own flow, however hard you pull.**""")
idea(r"""
lower the outlet pressure ──► the tube narrows there ──► the blood speeds up (same flow through less area)
the speed u reaches the local wave speed c (S = u/c = 1) ──► no signal can travel back upstream
──► the upstream flow no longer hears the outlet: the flow is choked
""", r"""**The speed index S = u/c plays the role of the Mach number of Ch. 15** (it is «16.195», written out in N101 below):
the tube's response to every push carries the factor $1-S^2$ of «16.204», which vanishes at S = 1, as $1-M^2$ did at a
nozzle throat.""")
note("N94 [B]", r"""The scene: a tube whose inside pressure falls below the outside collapses — circle → ellipse → two lobes that
touch. It happens to veins above the heart, coronary vessels squeezed by the heart muscle, an artery under a blood-pressure
cuff, airways during a forced cough. Collapsible tubes can also flutter and make sounds.

```
   ◯            ⬭              ∞
 open       elliptic       two lobes touching
 α = 1      α ≈ 0.6        α ≈ 0.25          (α = A/A₀, the area ratio)
```""")
note("N95 [C]", r"""Korotkoff sounds: a cuff pressed between systolic and diastolic pressure makes the artery a collapsible tube
that opens and closes each beat; the sounds a stethoscope hears mark the two pressures. The flow-limitation curve below is
the reason such tubes behave strangely.""")
note("N96 [B]", r"""The Starling resistor — a collapsible tube in a pressure chamber (inlet $p_1$, outlet $p_2$, chamber $p_e$):

| Case | What happens |
|---|---|
| 1 flow limitation | lowering $p_2$ stops raising Q |
| 2 pressure-drop limitation | the drop stops growing with the flow |
| 3 rising resistance | the tube narrows as the flow is forced |
| 4 flutter | with a downstream resistor the tube oscillates |
| 5 self-excited oscillations | the tube oscillates on its own |

Case 1 is computed below (ours, model).""")
note("N97 [B]", r"""Shapiro's 1-D theory: a general tube law «16.183»; small waves travel at «16.185» — C08's speed again
(Bramwell–Hill). The analogy: speed index S ↔ Mach number M · tube law ↔ equation of state · throat ↔ Laval throat (Ch. 15).""")
code(r"""
A0v = np.pi*vein["a0"]**2                                    # neutral area of our vein [m²]
print(ch16.bramwell_hill_speed(A0v, 1.5*vein["Kp"]/A0v, blood["rho"]), "m/s")   # (16.185) at α = 1: dp/dA = nK_p/A0 → 0.53 m/s
""")
P("P392", "bending stiffness", r"""A thin wall resists bending with a stiffness proportional to EI, where the second moment of area
of a wall of thickness h grows as h³: double the thickness, eight times the bending stiffness. A collapsed tube holds its
shape by bending, not by hoop tension.""", code=r"""
E, h1, h2 = 1e6, 1e-3, 2e-3                    # modulus [Pa], two thicknesses [m]
print((E*h2**3)/(E*h1**3))                     # 8.0
""")
note("N98 [B]", r"""Dimensionless pressure and area: «16.186», with «16.187» — α here is an **area ratio** (trap T2); the constant
in front of EI is not given, so the code requires it as a keyword. **N99 [B]** The collapsible tube law «16.188» with n = 3/2
— **a fit for collapse only (α ≤ 1)**; it is not meant for a distended tube.""")
code(r"""
print(ch16.bending_stiffness_parameter(1e6, 0.5e-3, 4e-3, 0.5, coeff=1.0), "Pa")   # K_p = coeff × E (h/a0)³/(1 − ν̂²)
""")
confusion(r"""(trap T9) P is **negative** when the tube is collapsed (α < 1), so we plot −P, which is positive there.""")
figure(r"""
al_c = np.linspace(0.12, 1.0, 200)                            # collapse branch: the fit's range
al_x = np.linspace(1.0, 1.6, 60)                              # beyond α = 1: outside the fit (formula drawn dashed, not the library)
fig, ax = plt.subplots(figsize=(7, 3.6))
ax.plot(al_c, -ch16.collapsible_tube_law(al_c, n=1.5), color=COLORS["teal"], lw=2, label="−P = α⁻ⁿ − 1 (16.188)")
ax.plot(al_c, -ch16.collapsible_tube_law(al_c, n=1.5, form="pure"), "--", color=COLORS["muted"], label="−P = α⁻ⁿ")
ax.plot(al_x, al_x**-1.5 - 1, ":", color=COLORS["rose"], label="outside the fit (α > 1)")
ax.axhline(0, color=COLORS["muted"], lw=0.6); ax.axvline(1, color=COLORS["muted"], ls=":", lw=0.8)
for a_, s_ in ((0.25, "two lobes"), (0.6, "ellipse"), (0.95, "circle")):          # cross-section shapes along the axis
    ax.text(a_ + 0.08*(a_ < 0.3), 10, s_, ha="center", fontsize=9)               # "two lobes" nudged right, clear of the steep curves
ax.set(xlabel="area ratio α = A/A₀ [–]", ylabel="−P = −(p − p_e)/K_p [–]", ylim=(-1, 15), title="A collapsing tube gets stiffer (fit for α ≤ 1)")
ax.legend(fontsize=8)
savefig(fig, "ch16", "c15_tube_law"); plt.show()
""", see=r"""−P against the area ratio for the fit «16.188» (teal) and the pure power law (grey dashed), on the collapse branch
α ≤ 1; the rose dotted piece beyond α = 1 is the same formula outside the range it was fitted for.""",
    read=r"""As the tube flattens (α falls) it takes ever larger negative transmural pressures to squeeze it further: the slope
dP/dα grows, so the wave speed «16.189» grows too (D28).""",
    change=r"""…n were 2? The curve would rise more steeply at small α — a tube that resists the last stage of collapse harder.""")
D("D28", r"The wave speed of a collapsible tube", ref="16.189", notes=r"**N100 [B]** «16.189».")
code(r"""
print(ch16.collapsible_wave_speed(np.array([1.0, 0.5, 0.25]), vein["Kp"], blood["rho"]), "m/s")   # (16.189): 0.53, 0.89, 1.50 m/s
""")
recap("R13", "Duct laws with friction", r"""Steady 1-D continuity «16.190» and momentum with friction and gravity «16.191», the wall
stress for turbulent «16.192» or laminar «16.193» flow — Ch. 15's duct laws with a wall force. **New here:** the area is set
by the tube law, not by the wall's shape. «16.193» uses an area $A_e$ the book does not define (slip S21; coded only as a
keyword).""", where="Ch. 15 §15.3, §15.8 (C03, C11)", core_id="C15")
note("N101 [B]", r"""The speed index «16.195» — the Mach number of a floppy tube; in Ch. 15 every response carried the factor
1 − M² (C06), here 1 − S².""")
code(r"""
print(ch16.speed_index(0.6, ch16.collapsible_wave_speed(0.5, vein["Kp"], blood["rho"])))   # S = u/c at α = 0.5 for u = 0.6 m/s
""")
remind("C15")
D("D29", r"The differential relations", ref="16.201", notes=r"**N102 [B]** «16.196»–«16.202» (the book writes $dA_e$ for $dA_0$ just before (16.196) — slip S22).")
slip("S23", r"the last term of «16.201p»", r"«16.201»", r"""The inner derivative is ∂P/∂α (the slope of the tube law), not ∂P/∂x.""")
code(r"""
print(ch16.collapsible_differentials_sympy())               # every residual 0, M = 2 − n for the power law
print(ch16.collapsible_differentials_sympy(printed_16_201=True)["dc_difference"])   # the printed term differs when P depends on x
""")
D("D30", r"Shapiro's equations and the condition at S = 1", ref="16.204")
confusion(r"""(trap T9, slip S24) The condition at S = 1 needs $\mathrm P\,dK_p/dx<0$; a collapsed tube has P < 0, so the stiffness must
**increase** along the flow — the text's "decrease" reads the sign of P wrongly.""")
worked("a vein at its limit", r"""
1. A collapsed segment at α = 0.25 with K_p = 250 Pa (easy numbers): c² = 1.5 × 250 × 8/1000 = 3, c = **1.73 m/s**.
2. Its area A = 0.25A₀; with A₀ = 50 mm² the choking flux is Q* = αA₀c = 0.25 × 50 × 10⁻⁶ × 1.73 = **21.7 mL/s**.
3. Our vein (A₀ = 50.3 mm², K_p = 200 Pa, L = 0.1 m, open inlet): friction narrows it until S = 1 at the outlet at
   **Q = 20.0 mL/s**, α = 0.31, p₂ − p_e = −946 Pa = −7.1 mmHg (ours, model).
4. Lower p₂ further: Q stays at 20.0 mL/s.""")
code(r"""
prm = dict(L=vein["L"], A0=A0v, Kp=vein["Kp"], n=vein["n"], rho=blood["rho"], mu=blood["mu"])   # our vein (model, ours)
print(f"Q* = {ch16.choking_flux(0.25, 50e-6, 250.0, 1000.0)*1e6:.2f} mL/s")   # the tiny example's choking flux
mr = ch16.collapsible_march(15e-6, 1.0, (0.0, prm["L"]), prm)  # (16.204) marched from an open inlet at 15 mL/s
print(f"15 mL/s: outlet α = {mr['alpha'][-1]:.4f}, S = {mr['S'][-1]:.3f}, choked: {mr['choked']}")
fl = ch16.flow_limitation_curve(np.linspace(-15, 0, 15 if FAST else 31)*MMHG, p1=0.0, pe=0.0, prm=prm)   # Starling Case 1
print(f"Q_choke = {fl['Q_choke']*1e6:.3f} mL/s, reached for p2 − pe ≤ {fl['p2_choke']:.1f} Pa = {fl['p2_choke']/MMHG:.2f} mmHg")
print("M =", ch16.shapiro_M(0.5, 1.5))                       # (16.206) for the power law: 2 − n = 0.5
""", explain=r"""
1. `collapsible_march` integrates «16.204» along the tube (friction only: $dA_0=dp_e=dz=dK_p=0$) and stops if S reaches 1.
2. `flow_limitation_curve` finds, for each outlet pressure, the flow whose outlet state matches it; below a critical
   outlet pressure the outlet is at S = 1 and the flow cannot rise any more.
3. `shapiro_M` evaluates «16.206»: ½ for n = 3/2.""")
scratch(r"""
def rhs(x, al, Q):                                             # dα/dx from (16.204) with friction only (laminar R′ = 8πμ/(αA0)²)
    c2 = 1.5*vein["Kp"]*al**-1.5/blood["rho"]                  # (16.189)
    S2 = (Q/(al*A0v))**2/c2                                    # speed index squared
    Rp = 8*np.pi*blood["mu"]/(al*A0v)**2                       # laminar resistance per length [Pa s/m⁴]
    return -al*Rp*Q/(blood["rho"]*c2*(1 - S2))                 # α·driver/(1 − S²)
n_e = 20000; dx = prm["L"]/n_e; al = 1.0                        # explicit Euler, 20 000 steps
for i in range(n_e):                                            # march down the tube, one small step dx at a time
    al += dx*rhs(i*dx, al, 15e-6)                               # Euler step: α ← α + dx·dα/dx at Q = 15 mL/s
assert np.isclose(al, mr["alpha"][-1], rtol=1e-3)               # the library's march
# a frictionless throat (ours; the orchestrator's case): S = 1 exactly where A0(x) is smallest
A0t, Kpt, rhot, al_t = np.pi*6e-3**2, 200.0, 1000.0, 0.35       # neutral area, stiffness, density, area ratio at the throat
Qs = ch16.choking_flux(al_t, A0t, Kpt, rhot, 1.5)              # critical flux: S = 1 at the throat
H = Kpt*ch16.collapsible_tube_law(al_t, n=1.5) + 0.5*rhot*(Qs/(al_t*A0t))**2   # Bernoulli constant K_pP + ½ρu²
A0f = lambda x: A0t*(1 + (x/0.05)**2)                           # a smooth throat at x = 0
a_in = ch16.collapsible_frictionless_alpha(Qs, H, A0f(-0.12), Kpt, rhot, 1.5, branch="sub")   # inlet state (subcritical)
mt = ch16.collapsible_march(Qs, a_in, (-0.12, 0.12), dict(A0=A0f, Kp=Kpt, n=1.5, rho=rhot, mu=1e-3, friction=False), rtol=1e-11, atol=1e-13)
assert a_in <= 1 and np.all(mt["alpha"] <= 1)                   # upstream stays on the collapse branch (α_t ≤ 0.25^(2/3) = 0.397)
assert abs(mt["x_choke"]) < 2e-3                                # S reaches 1 at the throat
print(f"Euler {al:.5f} vs library {mr['alpha'][-1]:.5f}; throat: Q* = {Qs*1e6:.2f} mL/s, inlet α = {a_in:.4f}, S = 1 at x = {mt['x_choke']*1e3:.2f} mm")
""", explain=r"""
1. An explicit Euler march of «16.204» with laminar friction (20 000 steps) reproduces the library's outlet area ratio.
2. A frictionless tube with a smooth throat (ours; ρ = 1000 kg/m³, $K_p$ = 200 Pa, a₀ = 6 mm, throat area ratio 0.35, small
   enough that the upstream tube stays collapsed, α ≤ 1): at the critical flux the march reaches S = 1 exactly at the throat —
   the floppy-tube twin of a Laval nozzle (Ch. 15).""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.8))
ax.plot(fl["p2"]/MMHG, fl["Q"]*1e6, "o-", color=COLORS["blue"], ms=3)   # blue: flow against outlet pressure (Starling Case 1)
ax.axhline(fl["Q_choke"]*1e6, ls="--", color=COLORS["muted"]); ax.axvline(fl["p2_choke"]/MMHG, ls=":", color=COLORS["rose"])
ax.text(fl["p2_choke"]/MMHG - 0.3, 5, "S = 1 at the outlet\n(flow limited)", ha="right", fontsize=8, color=COLORS["rose"])
ax.set(xlabel="outlet pressure p₂ − p_e [mmHg]", ylabel="flow Q [mL/s]", ylim=(0, 25), title="Starling Case 1: the flow stops rising (model, ours)")
for p2_, col in ((-2.0, COLORS["teal"]), (-5.0, COLORS["blue"]), (fl["p2_choke"]/MMHG, COLORS["rose"])):   # three outlet pressures
    st = ch16.starling_state(0.0, p2_*MMHG, 0.0, prm)          # the flow for that outlet pressure
    m_ = ch16.collapsible_march(min(st["Q"], fl["Q_choke"]*(1 - 1e-6)), 1.0, (0.0, prm["L"]), prm)
    bx.plot(m_["x"]*100, m_["alpha"], color=col, label=f"α, p₂ = {p2_:.1f} mmHg"); bx.plot(m_["x"]*100, m_["S"], "--", color=col)
bx.axhline(1, color=COLORS["muted"], lw=0.6)
bx.set(xlabel="x [cm]", ylabel="α (solid), S (dashed) [–]", ylim=(0, 1.1), title="The tube narrows until S = 1"); bx.legend(fontsize=7)
savefig(fig, "ch16", "c15_starling"); plt.show()
""", see=r"""Left: the flow through our collapsible vein against the outlet pressure (open inlet at the chamber pressure); right: the
area ratio (solid) and speed index (dashed) along the tube for three outlet pressures, the last one at the critical value.""",
    read=r"""Lowering the outlet pressure first raises the flow steeply, as in any pipe; below about −4 mmHg the flow is within a
few per cent of its limit, and from −7.1 mmHg on it is exactly flat: the outlet has reached S = 1 and news of a lower outlet pressure can no longer travel upstream against the flow. On the right,
friction narrows the tube along its length and raises S; at the critical case S reaches 1 exactly at the outlet.""",
    change=r"""…the blood were more viscous? Friction would narrow the tube faster and the flow would limit at a lower value.""")
explainer("collapsible_tube_choking", "Why can't you suck harder through a floppy straw?",
          r"""lowering the outlet pressure you watch Q rise and stall as the outlet reaches S = 1; in the throat mode the drivers of
«16.204» change sign exactly where 1 − S² does; the giraffe mode balances gravity against friction — three regimes of one
equation, switched live.""",
          ["Drag p₂ from 0 to −15 mmHg: the flow plateaus at about −7 mmHg.", "Throat mode: watch the term bars cross zero at S = 1.",
           "Giraffe mode: the vein stays partly collapsed."])
md(r"""**The giraffe** (part of C15). A giraffe's jugular vein runs downhill about 2 m from head to heart, dz/dx = −1: gravity
pulls the blood faster than friction alone allows, and the vein runs partly collapsed, its pressure near the outside
pressure — the gravity driver of «16.204» balancing friction (ours, a model in the explainer above).""")
whatif(r"""…the fluid could not flow at all below a certain stress? Then the blood in the middle of a pipe stops shearing — C16, the
payoff of C02.""")

# ---------------------------------------------------------------------------------------------------------------------
# C16 — Casson pipe flow
# ---------------------------------------------------------------------------------------------------------------------
core("C16", r"A plug riding on a ring: «16.210» and the Casson flux «16.220»", "Does a yield stress stop blood flowing in a vessel?")
note("N103 [C]", r"""Below about 200 s⁻¹ blood is clearly non-Newtonian (C02); Casson's law is the usual model.""")
problem(r"""
Squeeze toothpaste from a tube: it comes out as a solid cylinder that slides on a thin sheared layer at the nozzle wall. A
fluid with a yield stress does that in any pipe: **near the centre the stress is too small to shear it, so the centre moves
as a plug**. Blood's yield stress is tiny, but the same picture holds — and it costs flow.""")
idea(r"""
stress τ = (Δp/2L) r :  0 on the axis ─────────► τ_w at the wall
   r < a_c :  τ < τ_y  →  no shear: a solid plug
   r > a_c :  τ > τ_y  →  a sheared ring
""")
note("N104 [B]", r"""Casson's law with its yield branch: «16.207»; $K_c$ is dimensional (trap T14).""")
code(r"""
print(ch16.casson_shear_rate(np.array([0.5, 1.0, 4.0])*ty, ty, Kc))   # zero below τ_y, positive above (16.207)
""")
D("D31", r"The stress line, the plug radius and the minimum pressure drop", ref="16.210",
  notes=r"**N105 [B]** the linear stress «16.209» (its first member's sign printed without the minus — slip S25) and the plug radius «16.210».")
slip("S26", r"its worked minimum pressure drop in dyn/cm² 1000 times too small", r"$\Delta p_{\min}=2L\tau_y/a$, computed below with our numbers in Pa and mmHg")
note("N106 [B]", r"""(trap T21) If $a_c\ge a$ the whole cross-section is below the yield stress: no shear anywhere, so **no flow** —
the functions return zero flux, not a negative plug.""")
P("P393", "algebra with square roots", r"""Three moves the Casson algebra uses in every step: $(\sqrt x-\sqrt y)^2=x+y-2\sqrt{xy}$;
$\int r^{1/2}dr=\tfrac23r^{3/2}$; and a polynomial in $t=\sqrt\xi$ is still a polynomial (t² = ξ, t⁸ = ξ⁴).""", code=r"""
x, y = 9.0, 4.0
print((np.sqrt(x) - np.sqrt(y))**2, x + y - 2*np.sqrt(x*y))     # 1.0 1.0
print(quad(np.sqrt, 0, 4)[0], 2/3*4**1.5)                        # 5.333 twice
""")
D("D32", r"The Casson profile and the plug speed", ref="16.218", notes=r"**N107 [B]** the profile «16.215»/«16.216» and the plug speed «16.218».")
D("D33", r"The Casson flux", ref="16.220", notes=r"**N108 [B]** «16.220».")
worked("a plug a quarter of the radius", r"""
a = 0.5 mm, τ_y = 4.8 mPa, K_c² = 2.81 mPa s, Δp/L = 76.8 Pa/m:
1. τ_w = (Δp/L)a/2 = 19.2 mPa.
2. a_c/a = τ_y/τ_w = 0.25 → a_c = 0.125 mm.
3. Flux bracket: 1 − (16/7)(0.5) + (4/3)(0.25) − (1/21)(0.25)⁴ = **0.190**.
4. Poiseuille with μ = K_c²: πa⁴(Δp/L)/(8K_c²) = 6.71 × 10⁻¹⁰ m³/s; Casson: **1.28 × 10⁻¹⁰ m³/s**.
5. The minimum drop for L/a = 1000: 2Lτ_y/a = 2 × 1000 × 0.0048 = **9.6 Pa = 0.072 mmHg**.""")
code(r"""
Kc = wc["Kc"]                                                 # K_c again (the sympy cell above used the name for a symbol)
aP, LP = casson_pipe["a"], casson_pipe["L"]                  # pipe radius and length [m]
dpP = 76.8*LP                                                 # pressure drop [Pa] (gradient 76.8 Pa/m)
ac = ch16.casson_plug_radius(dpP, LP, ty)                     # (16.210) [m]
uc = ch16.casson_core_velocity(aP, dpP, LP, ty, Kc)           # (16.218) [m/s]
QC = ch16.casson_pipe_flow_rate(aP, dpP, LP, ty, Kc)          # (16.220) [m³/s]
dpmin = ch16.casson_min_pressure_drop(LP, aP, ty)             # 2Lτ_y/a [Pa]
print(f"a_c = {ac*1e3:.4f} mm ({ac/aP:.3f} a), u_c = {uc*1e3:.4f} mm/s, Q = {QC:.4e} m³/s, Δp_min = {dpmin:.2f} Pa = {dpmin/MMHG:.3f} mmHg")
print("bracket:", ch16.casson_flux_ratio(0.25), "| Poiseuille (μ = K_c²):", ch16.poiseuille_flow_rate(aP, dpP, LP, Kc**2))
""")
scratch(r"""
r_ring = np.linspace(ac, aP, 4001)                            # the sheared ring
u_ring = ch16.casson_pipe_velocity(r_ring, aP, dpP, LP, ty, Kc)   # (16.215) [m/s]
Q_mine = np.pi*ac**2*uc + np.trapezoid(2*np.pi*r_ring*u_ring, r_ring)   # plug + ring (16.219)
assert np.isclose(Q_mine, QC, rtol=1e-6)                      # equals (16.220)
assert np.isclose(ch16.casson_pipe_flow_rate(aP, dpP, LP, 0.0, Kc), ch16.poiseuille_flow_rate(aP, dpP, LP, Kc**2), rtol=1e-12)   # τ_y → 0: Poiseuille
assert ch16.casson_pipe_flow_rate(aP, 0.9*dpmin, LP, ty, Kc) == 0.0   # below Δp_min: no flow (trap T21)
print(f"plug + ring = {Q_mine:.6e} m³/s")
""", explain=r"""Plug flux plus a trapezoid sum over the ring equals «16.220»; a zero yield stress gives Poiseuille's flux with
μ = K_c²; and below the minimum pressure drop the flux is exactly zero.""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.8))
rr = np.linspace(-aP, aP, 801)                                 # across the diameter [m]
uC = ch16.casson_pipe_velocity(np.abs(rr), aP, dpP, LP, ty, Kc)   # Casson profile (plug value inside a_c)
uP = ch16.poiseuille_velocity(np.abs(rr), aP, dpP, LP, Kc**2)  # Poiseuille with μ = K_c²
ax.plot(rr*1e3, uC*1e3, color=COLORS["blue"], lw=2, label="Casson (16.215)")
ax.plot(rr*1e3, uP*1e3, "--", color=COLORS["muted"], label="Poiseuille, μ = K_c²")
ax.axvspan(-ac*1e3, ac*1e3, color=COLORS["teal"], alpha=0.15); ax.text(0, 0.05, "plug", ha="center", fontsize=8)
ax2 = ax.twinx(); ax2.plot(rr*1e3, dpP/(2*LP)*np.abs(rr)*1e3, color=COLORS["rose"], lw=1); ax2.axhline(ty*1e3, color=COLORS["rose"], ls=":")
ax2.set_ylabel("stress τ [mPa] (rose; dotted τ_y)")
ax.set(xlabel="r [mm]", ylabel="u [mm/s]", title="A plug rides on a sheared ring (a_c = a/4)"); ax.legend(fontsize=7, loc="upper left")
xi = np.linspace(0, 1, 300)                                    # a_c/a
bx.plot(xi, ch16.casson_flux_ratio(xi), color=COLORS["amber"], lw=2, label="Q/Q_Pois (16.220)")
bx.plot(xi, np.maximum(0, 1 - 16/7*np.sqrt(xi)), "--", color=COLORS["muted"], label="first-order 1 − (16/7)√ξ")
bx.plot(0.25, ch16.casson_flux_ratio(0.25), "o", color=COLORS["blue"])   # blue dot: our example, a_c = a/4
bx.set(xlabel="plug size a_c/a [–]", ylabel="flux / Poiseuille's [–]", title="The first few per cent of plug cost the most"); bx.legend(fontsize=8)
savefig(fig, "ch16", "c16_casson_pipe"); plt.show()
""", see=r"""Left: the Casson velocity profile across the pipe (blue) with the plug shaded teal, Poiseuille's profile for
μ = K_c² at the same pressure drop (grey dashed), and the linear stress (rose, right axis) crossing τ_y at the plug's edge.
Right: the flux relative to Poiseuille's against the plug size, with the first-order √ξ approximation.""",
    read=r"""Inside the plug the profile is flat; outside it bends to zero at the wall. The flux falls steeply at small plug
sizes because of the $\sqrt{a_c/a}$ term: a plug 1 % of the radius already costs a fifth of the flow.""",
    change=r"""…the pressure gradient were halved? $a_c/a$ would double to 0.5 and the bracket of «16.220» would fall to
`ch16.casson_flux_ratio(0.5)` = 0.047: the flux would drop to 4.7 % of Poiseuille's (printed below), and Poiseuille's own
flux halves too.""")
code(r"""
print(f"bracket of (16.220) at a_c/a = 0.25: {ch16.casson_flux_ratio(0.25):.4f}; at 0.5: {ch16.casson_flux_ratio(0.5):.4f}")   # 0.190 and 0.047
""")
plotly(r"""
xis = list(np.round(np.linspace(0.0, 0.95, 20), 3))           # τ_y/τ_w = a_c/a (20 positions)
sP = np.linspace(-1, 1, 301)                                   # r/a across the diameter
def casson_prof(xi_):                                          # profiles in units of the Poiseuille centre speed a τ_w/(2K_c²)
    tw = 1.0; ty_ = xi_*tw; s = np.abs(sP)                     # wall stress 1, yield stress xi
    u_ = np.where(s <= xi_, (1 - xi_**2) - 8/3*np.sqrt(xi_)*(1 - xi_**1.5) + 2*xi_*(1 - xi_),   # inside the plug (|r|/a ≤ a_c/a): the flat plug speed
                  (1 - s**2) - 8/3*np.sqrt(xi_)*(1 - s**1.5) + 2*xi_*(1 - s))   # (16.216) divided by aτ_w/(2K_c²)
    return {"Casson (16.216)": (sP, u_), "Poiseuille, μ = K_c²": (sP, 1 - sP**2)}
fig = slider_figure(casson_prof, "τ_y/τ_w", xis, xlabel="r / a [–]", ylabel="u / (aτ_w/2K_c²) [–]",
                    title="The plug widens as the yield stress approaches the wall stress", yrange=[0, 1.05], active=5)
fig.show()
""", explain=r"""**Figure F7.** Slider: the ratio τ_y/τ_w = a_c/a from 0 to 0.95 (opens on 0.25, our example). **What you see.** The Casson
profile «16.216» in units of the Poiseuille centre speed, beside Poiseuille's parabola. **How to read it.** At 0 the two are
one; as the ratio grows the flat plug widens and the whole profile sinks. **What would change if…** the ratio reached 1: the
plug would fill the pipe and the flow would stop (trap T21).""")
explainer("casson_blood_flow", "Does a yield stress stop the blood?",
          r"""dragging the pressure gradient down you watch the plug widen until it fills the pipe and the flow stops; dragging a
shear-rate cursor on the law shows the secant and tangent separate — two linked views of one material.""",
          ["Lower Δp/L toward 2τ_y/a: the plug fills the pipe.", "Set τ_y = 0: Poiseuille returns.",
           "Click inside the plug: the shear rate is zero there."])
whatif(r"""…the fluid were sap in a tree instead of blood in an artery? The flow is Poiseuille's again, but now gravity and dissolved
sugars join the pressure — C17.""")

# =====================================================================================================================
# A.4 §16.4 — C17
# =====================================================================================================================
nb.section("16.4", "Introduction to the Fluid Mechanics of Plants", intro="**What is this section about?** How water rises from "
           "soil to leaf (xylem) and how sugar travels from leaf to root (phloem): one 'water potential' organises both.")
core("C17", r"Water potential «16.221»", "How does water reach the top of a 30-metre tree without a pump?")
note("N109 [B]", r"""The scene: plants make up most of the Earth's living mass. Their vascular tissue has two pipe systems: xylem
(water up) and phloem (sugar sap wherever it is needed).

```
   air  ~~~  water vapour leaves the leaf (transpiration)
   leaf  ◄── lowest ψ
    │  xylem: water pulled up under tension        phloem: sugar sap pushed from leaf to where it is needed
   stem
    │
   root ◄── water enters from the soil (highest ψ)
   soil
```""")
problem(r"""
Water in a tree climbs tens of metres with no heart. It is not pushed up; **it is pulled**. Water evaporates from the leaves;
the water column in the narrow xylem pipes is pulled up after it, hanging from the leaves under tension — a negative
pressure. To keep track of where water wants to go we need one number that combines pressure, dissolved salt and height.""")
P("P384", "osmosis and van 't Hoff's law", r"""Water moves through a membrane toward the side with more dissolved solute. The
pressure that stops it is the osmotic pressure; for a dilute solution $\Pi=cR_uT$ (c in mol/m³) — the solute lowers the
water's "potential" by Π. (Chemical potential is the deeper idea; named only.)""", code=r"""
Ru, T = 8.314462618, 293.15                    # gas constant [J/(mol K)], temperature [K]
for c in (10, 100, 300):                       # concentrations [mol/m³]
    print(c, c*Ru*T/1e6)                       # 0.024, 0.244, 0.731 MPa
""")
idea(r"""
ψ  =   p     −   R_u T c    +   ρ g z
     pressure    dissolved      height
     pushes      salt pulls     costs
water moves from high ψ to low ψ:   soil (highest) ──► root ──► stem ──► leaf ──► air (lowest)
""", r"""**One number, the water potential «16.221», says which way water moves** (c is the solute concentration of the
osmosis primer above, R_u the gas constant). A tree needs no pump: evaporation keeps the leaf's ψ lowest, and the column is
pulled up toward it.""")
recap_md("Hydrostatics", "Ch. 1 D05", r"""«1.8»: in a fluid at rest the pressure falls by ρg per metre of height.""")
D("D34", r"The three terms of the water potential and the tension that lifts sap", ref="16.221")
worked("a 30 m tree", r"""
1. Soil water potential ψ_soil = −0.1 MPa (moist soil, ours).
2. Still column: p(30 m) = −0.1 − 1000 × 9.81 × 30/10⁶ = **−0.394 MPa** (about −3.9 bar).
3. Moving sap, 1 mm/s in a 25 µm conduit: Poiseuille gradient 8μu/a² = 8 × 10⁻³ × 10⁻³/(25 × 10⁻⁶)² = 12.8 kPa/m; over 30 m
   0.384 MPa more: **−0.778 MPa** at the top.
4. A leaf cell with 0.3 mol/L of solutes: −R_uTc = −0.731 MPa — its ψ is low enough to draw water from the xylem.""")
code(r"""
psi_top = ch16.water_potential(-0.394e6, 0.0, tree["T"], 30.0)   # (16.221) at the top of a still column [Pa]
terms30 = ch16.water_potential_terms(-0.394e6, 0.0, tree["T"], 30.0)   # pressure, osmotic, gravity parts [Pa]
print(f"psi(30 m) = {psi_top/1e6:.4f} MPa (≈ ψ_soil)", {k: round(v/1e6, 4) for k, v in terms30.items()})
print(f"osmotic pressure of 300 mol/m³: {ch16.osmotic_pressure(300.0, tree['T'])/1e6:.4f} MPa (R_u = TH.R_U/1000 per mole)")
col = ch16.xylem_column(np.linspace(0, tree["H"], 61), psi_soil=tree["psi_soil"], u=tree["u_xylem"], a=tree["a_xylem"], mu=tree["mu_w"])
print(f"moving sap: p(top) = {col['p'][-1]/1e6:.4f} MPa; margin to an illustrative −2 MPa cavitation limit: {ch16.cavitation_margin(col['p'][-1], p_cav=-2.0e6)/1e6:.3f} MPa")
""", explain=r"""
1. `water_potential` is «16.221» in SI (pressure gauge in Pa, concentration in mol/m³, height in m); at the top of a still
   column the pressure term and the height term cancel to give back ψ_soil.
2. `osmotic_pressure` is van 't Hoff's law with the universal gas constant per mole (`TH.R_U` is per kmol, so /1000).
3. `xylem_column` adds the Poiseuille drop of moving sap (ours): −0.778 MPa at the top, still 1.2 MPa away from our
   illustrative cavitation limit.""")
scratch(r"""
Ru = TH.R_U/1000                                               # J/(mol K)
p_, c_, z_ = -0.394e6, 50.0, 30.0                              # pressure [Pa], solutes [mol/m³], height [m]
psi_mine = p_ - Ru*tree["T"]*c_ + 1000.0*TH.G0*z_              # (16.221) term by term, SI
assert np.isclose(psi_mine, ch16.water_potential(p_, c_, tree["T"], z_))   # the library computes exactly this sum
psi_bar = ch16.water_potential_book_units(p_/1e5, c_*1e-6, tree["T"], z_*100)   # the book's units: bar, mol/cm³, cm, g = 980 cm/s²
assert np.isclose(ch16.bar_to_Pa(psi_bar), ch16.water_potential(p_, c_, tree["T"], z_, g=9.80), rtol=1e-4)   # ρgz ÷ 10⁶ to be in bar (slip S30)
print(f"SI: {psi_mine/1e6:.5f} MPa; book units: {psi_bar:.4f} bar")
""", explain=r"""The water potential written term by term in SI equals the library's; the book's mixed units (bar, mol/cm³, cm, with
ρgz in dyn/cm² divided by 10⁶ to be in bar — slip S30) give the same number once converted.""")
figure(r"""
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.0), gridspec_kw=dict(width_ratios=[1.3, 1]))
still = ch16.xylem_column(col["z"], psi_soil=tree["psi_soil"], u=0.0, a=tree["a_xylem"], mu=tree["mu_w"])   # no flow
ax.plot(still["p"]/1e6, col["z"], color=COLORS["orange"], ls="--", label="p, still sap")
ax.plot(col["p"]/1e6, col["z"], color=COLORS["orange"], lw=2, label="p, sap at 1 mm/s")
ax.plot(still["psi"]/1e6, col["z"], color=COLORS["accent"], ls="--", label="ψ, still sap (constant)")
ax.plot(col["psi"]/1e6, col["z"], color=COLORS["accent"], lw=2, label="ψ, moving sap (falls upward)")
hcap = ch01.capillary_rise(0.0728, np.pi/2, 1000.0, tree["a_xylem"])   # Ch. 1: capillary rise in a 25-µm conduit [m]
ax.axhline(hcap, color=COLORS["teal"], ls=":"); ax.text(-0.62, hcap + 0.6, f"capillarity alone: {hcap:.2f} m", fontsize=8, color=COLORS["teal"])
ax.axvline(0, color=COLORS["muted"], lw=0.6)
ax.set(xlabel="pressure or water potential [MPa]", ylabel="height z [m]", xlim=(-1.0, 0.45), title="Sap hangs under tension (ours)")
ax.legend(fontsize=7, loc="upper right")
places = ["soil", "root", "stem top", "leaf cell"]               # four stations along the water's path
pr = np.array([-0.1, -0.12, col["p"][-1]/1e6, -0.10])           # pressure term [MPa] (illustrative, ours)
os_ = np.array([0.0, 0.0, 0.0, -ch16.osmotic_pressure(300.0, tree["T"])/1e6])   # solute term [MPa]
gr = np.array([0.0, 0.0, 1000*TH.G0*30/1e6, 1000*TH.G0*30/1e6])   # height term [MPa]
xs4 = np.arange(4)                                              # their bar positions 0, 1, 2, 3
bx.bar(xs4, pr, color=COLORS["orange"], label="pressure p"); bx.bar(xs4, os_, bottom=np.minimum(pr, 0), color=COLORS["accent"], label="solutes −R_uTc")
bx.bar(xs4, gr, color=COLORS["teal"], alpha=0.7, label="height ρgz")
bx.plot(xs4, pr + os_ + gr, "k_", ms=25, mew=2, label="ψ (total)")
bx.set_xticks(xs4, places, fontsize=8); bx.axhline(0, color=COLORS["muted"], lw=0.6)
bx.set(ylabel="[MPa]", title="ψ falls from soil to leaf (illustrative)"); bx.legend(fontsize=7, loc="lower left")
savefig(fig, "ch16", "c17_xylem"); plt.show()
""", see=r"""Left: pressure (orange) and water potential (purple) up a 30 m tree, for still sap (dashed) and sap moving at 1 mm/s
in 25-µm conduits (solid); the dotted line is how high capillarity alone could lift water. Right: the three terms of ψ at
four places along the path (illustrative values, ours) and their total.""",
    read=r"""With still sap ψ is the same everywhere — no flow without a gradient of ψ — and the pressure falls by ρg per metre.
The pressure plotted is gauge: it is already −0.1 MPa (= ψ_soil) at the ground, and the absolute pressure p + p_atm turns
negative — a true tension — about 13 cm above the soil, where p passes −0.1013 MPa. Moving sap needs ψ to fall upward, which pushes the pressure further negative (−0.78 MPa at the
top). Water moves from high ψ (soil) to low ψ (leaf), where the solute term pulls hardest.""",
    change=r"""…the tree were twice as tall? The tension at the top would roughly double — trees stop growing roughly where their
xylem can no longer hold the column.""")
note("N110 [B]", r"""Xylem: dead conduits (tracheids, vessels); evaporation at the leaves (transpiration) pulls a cohesive column
under tension; too much tension and the column breaks (cavitation); at night root pressure can push water out at leaf edges
(guttation). Capillarity alone: Ch. 1's $h=2\sigma\sin\alpha/(\rho gR)$ (Example 1.1) gives 0.59 m for a 25 µm conduit — far
too little (ours).""")
code(r"""
print(f"{ch16.capillary_rise(0.0728, np.pi/2, 1000.0, 25e-6):.3f} m")   # Ch. 1: σ = 0.0728 N/m, full wetting, R = 25 µm
""")
note("N111 [B]", r"""Xylem flow is creeping flow: Re = u·2a/ν = 10⁻³ × 50 × 10⁻⁶/10⁻⁶ = **0.05**; the Poiseuille drop per metre
(12.8 kPa) is comparable with ρg (9.8 kPa).""")
slip("S28", r"a kinematic viscosity ten times water's beside its xylem Reynolds number", r"water's $\nu=10^{-6}$ m²/s, which the quoted Reynolds number needs (the book's numbers stay private)")
code(r"""
print(ch16.xylem_reynolds(tree["u_xylem"], tree["a_xylem"], tree["nu_w"]))   # Re = u·2a/ν = 0.05
""")
note("N112 [B]", r"""Phloem: living sieve tubes joined by perforated sieve plates; sugar loaded at a source raises the solute term,
water follows by osmosis, the raised pressure pushes the sap to a sink where sugar is unloaded (the Münch pressure-flow idea)
— ψ again.

```
 source (leaf): sugar in → ψ falls → water enters → p rises ══► sap flows ══► sink (root, fruit): sugar out → p falls
```""")
note("N113 [B]", r"""The pressure drop along a sieve tube with its plates: «16.222» (Δp′ the entry loss at a plate, given in the
book).""")
code(r"""
Qph = phloem["u"]*np.pi*phloem["a"]**2                        # sap flux in one sieve tube [m³/s]
st_ = ch16.sieve_tube_pressure_drop(Qph, phloem["mu"], phloem["a"], phloem["r_pore"], phloem["L"], phloem["ell"], phloem["N"])   # (16.222)
print({k: round(float(v), 3) for k, v in st_.items()}, "Pa")  # the plates carry most of the drop (illustrative sizes, ours)
lim = ch16.sieve_tube_pressure_drop(Qph, phloem["mu"], phloem["a"], phloem["a"], phloem["L"], 0.0, 1)["total"]   # no plates (r = a, ℓ = 0)
print(f"Poiseuille limit {lim:.3f} Pa vs 8μQL/(πa⁴) = {8*phloem['mu']*Qph*phloem['L']/(np.pi*phloem['a']**4):.3f} Pa")
""")
plotly(r"""
us = list(np.round(np.linspace(0, 3e-3, 13), 5))              # sap speeds 0 … 3 mm/s
zz = np.linspace(0, tree["H"], 61)                            # heights [m]
def column(u_):                                               # ψ terms and p along the trunk for one sap speed
    c_ = ch16.xylem_column(zz, psi_soil=tree["psi_soil"], u=u_, a=tree["a_xylem"], mu=tree["mu_w"])
    return {"p": (c_["p"]/1e6, zz), "ψ": (c_["psi"]/1e6, zz), "height term ρgz": (c_["gravity"]/1e6, zz),
            "viscous loss": (-c_["viscous"]/1e6, zz)}
fig = slider_figure(column, "u", us, unit="m/s", xlabel="[MPa]", ylabel="height z [m]",
                    title="Faster sap, deeper tension", xrange=[-1.4, 0.35], yrange=[0, 30], active=4)
fig.show()
""", explain=r"""**Figure F8.** Slider: the sap speed in a 25-µm conduit, 0 to 3 mm/s (opens on our 1 mm/s). **What you see.** Pressure,
water potential, the height term and (minus) the accumulated viscous loss against height. **How to read it.** At u = 0 the ψ
line is vertical (no flow without a gradient); every extra mm/s tilts ψ and pushes the top pressure 0.38 MPa further
negative. **What would change if…** the conduits were twice as wide: the viscous loss would fall fourfold
(Poiseuille, 8μu/a²).""")
explainer("xylem_water_potential", "How does water climb a 30-metre tree?",
          r"""dragging height, sap speed and soil dryness you watch the three terms of ψ trade along the trunk while the pressure
goes below zero and a warning appears near the illustrative cavitation limit — three inputs at once, which one static
figure cannot show.""",
          ["Set the sap speed to 0: ψ is the same at every height.", "Pick the 100 m giant.", "Drought: lower ψ_soil to −1.5 MPa."])
whatif(r"""…the tree were twice as tall? The tension doubles — trees stop growing roughly where their xylem can no longer hold the
column.""")

# =====================================================================================================================
# A.5 exercises, literature, summary
# =====================================================================================================================
md(r"""---

## Exercises, literature and what comes next""")
nb.pointer(X(r"""**S01** Exercises 16.1–16.6 (not reproduced): 16.1 extends C06 (an elastic tube whose radius varies along it); 16.2 is
the dimensional analysis of C03 (careful: its Strouhal number is aω/U — trap T6); 16.3 a stenosis; 16.4 Shapiro's pure
pressure–gravity flows (C15; the citation's year is misprinted — slip S33); 16.5–16.6 generalise C16 to power-law and
Herschel–Bulkley fluids (`ch16.power_law_pipe_flow_rate`, `ch16.herschel_bulkley_pipe_flow_rate`, ours)."""))
nb.pointer(r"""**S02** Acknowledgment, literature cited and supplemental reading: not reproduced. Sources we read first-hand are listed
in the next cell.""")
md(r"""**Sources read first-hand (public).** Secomb & Pries (2013), *C. R. Physique* 14, 470 (the in-vitro viscosity law of C05);
Wikipedia articles "Womersley number", "Pulsatile flow", "Moens–Korteweg equation", "Windkessel effect", "Murray's law",
"Hemorheology", "Entrance length (fluid dynamics)"; Lubashevsky & Gafiychuk, arXiv:adap-org/9909003 (Murray's cost).
Open items: Dean (1928) and Womersley (1955) were not read first-hand; the Navier–Stokes check of Dean's K⁴ coefficient
(0.0119312, slip S34) is our own calculation, not a published value.""")

nb.summary(
    clicked=[X(s) for s in (
        r"**C01** The heart's work per beat is the area of its pressure–volume loop, about a joule: $W=-\oint p\,dV$ (ours).",
        r"**C02** Blood thins with shear: Casson's law «16.6» is a straight line in square roots, with the yield stress as intercept.",
        r"**C03** «16.12» compares a vessel with its own Stokes layer: $\alpha=\sqrt2\,a/\delta_e$.",
        r"**C04** «16.20»: halve the radius, sixteen times the resistance — the arterioles' lever on blood pressure.",
        r"**C05** A thin plasma layer makes blood seem thinner in narrow tubes: «16.33».",
        r"**C06** An artery is a spring: the tube law «16.41», with a hoop stress «16.39».",
        r"**C07** The aorta stores each beat and drains it with time constant RK: «16.49»; the mean is $R\bar Q$, compliance sets the swing.",
        r"**C08** The pulse is a wave of wall stretch, «16.62», metres per second while the blood barely moves.",
        r"**C09** Womersley's profile «16.82»: a breathing parabola at small α, a lagging plug at large α.",
        r"**C10** The flow amplitude falls as 8/α² and lags up to 90°: «16.99», «16.100».",
        r"**C11** Viscosity slows and damps the pulse by amounts ∝ 1/α: the pressure root of «16.136».",
        r"**C12** Murray: at the best radius pumping is half the upkeep, so «16.140» gives Q ∝ a³ and «16.144».",
        r"**C13** A fork reflects by admittance mismatch, «16.153», and energy balances: $\mathrm R^2+(\sum Y_j/Y_1)\mathrm T^2=1$.",
        r"**C14** A bend makes two vortices and pushes the fast core outward; the flux falls as «16.178» (K⁴ coefficient 0.0119312 with slip S34 corrected).",
        r"**C15** A floppy tube chokes at S = 1 — flow limitation: «16.204».",
        r"**C16** A yield stress makes a plug and costs flow through $\sqrt{a_c/a}$: «16.220».",
        r"**C17** Sap is pulled up under tension: «16.221».",
    )],
    feeds_forward=["This is the book's last chapter: the tools reused here (Poiseuille, Stokes layers, d'Alembert, "
                   "differential relations with a 1 − S² factor, perturbation series) are the toolkit for any applied fluid problem — including the "
                   "climate-system flows of Ch. 13 (Ekman layers are Stokes layers in a rotating frame; wave reflection and "
                   "choking reappear in ocean straits)."],
    left_out=["swimming and flying (not in this chapter; Ch. 14's thrust-from-lift explainer is the closest)",
              "turbulence in the aorta (Ch. 12)", "large-Dean-number numerics (named in N93)",
              "collapsible-tube oscillations (named in N96)"],
)

# A "> ⚠️" callout written in the middle of a line prints a literal ">": start every such callout on its own line.
_MID = re.compile(r"(?m)^(?P<pre>[^\n]*?\S)[ \t]+> ⚠️")


def _own_line(text: str) -> str:
    def fix(m):
        pre = m.group("pre")
        if re.fullmatch(r"\s*(\d+\.|[-*])", pre):
            return m.group(0)
        return pre + "\n\n> ⚠️"
    prev = None
    while prev != text:
        prev, text = text, _MID.sub(fix, text)
    return text


for _c in nb.cells:
    if _c.cell_type == "markdown":
        _c.source = _own_line(_c.source)
        # A plain number in brackets, "(0.1 Pa)", reads to the coverage check like an equation number cited without its
        # equation; in a cell with no maths such a number is set as maths, "($0.1$ Pa)" — it is a value, not a citation.
        if "$" not in _c.source:
            _c.source = re.sub(r"\((\d+\.\d+)(?=[\s,)–-])", r"($\1$", _c.source)

# Plotting lines left without a comment get the comment that says what they draw (novice-grade: every line commented).
_AUTO = [
    (r"^fig, .*plt\.subplots\(", "a new figure and its axes"), (r"twinx\(\)", "a second y-axis sharing the x-axis"),
    (r"\.(plot|semilogx|semilogy|loglog)\(\[\], \[\]", "an empty line, filled in frame by frame by update()"),
    (r"\.(plot|semilogx|semilogy|loglog)\(.*label=f?\"([^\"]+)\"", None),
    (r"\.scatter\(", "markers (the tracers)"),
    (r"\.bar\(", "bars"), (r"\.(axhline|axvline)\(", "a reference line"), (r"\.(axhspan|axvspan)\(", "a shaded band"),
    (r"\.annotate\(", "a label with an arrow"), (r"\.text\(", "a text label on the plot"),
    (r"\.contourf\(", "filled contours of the axial speed"), (r"\.contour\(", "streamlines (contours of ψ0)"),
    (r"\.set\(", "axis labels, limits and title"), (r"set_(x|y)label|set_title|suptitle", "axis label or title"),
    (r"\.legend\(", "the legend"), (r"set_xticks", "tick positions and their labels"),
    (r"set_(x|y)scale", "logarithmic axis"), (r"savefig\(", "save to outputs/ch16/ and show the figure"),
    (r"\.axis\(\"off\"\)", "no axes: a drawing, not a graph"), (r"add_patch|plt\.Circle", "a circle"),
    (r"^return\b", "what this function hands back"), (r"slider_figure\(", "build the slider figure (every position computed now)"),
    (r"update_xaxes|update_yaxes", "logarithmic axis on the web figure"), (r"fig\.show\(\)", "display the figure"),
    (r"show_animation\(", "render the frames and embed the player"), (r"^global\b", "the tracers' state is kept between frames"),
    (r"\.set_(data|offsets|height|radius|text)\(", "move the artist to this frame"), (r"^print\(", "print the result"),
    (r"^(lines|out|waves) = (\[\]|\{\})", "collect the artists / results"),
    (r"^(\w+(?:, \w+)*) = ch16\.(\w+)\(", None),
    (r"^[\]\)}]|[,(\[{]$", "(continued)"),
]


_KWARG = [   # keyword-argument lines of a call split over several lines
    (r"title=", "the figure's title: its message"),
    (r"(xrange|yrange|active)=", "fixed axis ranges and the slider position it opens on"),
    (r"arrowprops=", "the arrow's style"), (r"textcoords=", "the label offset is measured in points"),
    (r"ha=", "text alignment and size"),
]


def _bracket_balance(line: str) -> int:
    """Opened minus closed brackets on a line, ignoring those inside string literals and comments."""
    bare = re.sub(r"(\"[^\"]*\"|'[^']*')", "", line.split("  #")[0])
    return sum(bare.count(c) for c in "([{") - sum(bare.count(c) for c in ")]}")


def _auto_comment(src: str) -> str:
    out = []
    depth = 0                                # open brackets carried over from the lines above
    for ln in src.split("\n"):
        st = ln.strip()
        inside, depth = depth > 0, depth + _bracket_balance(ln)
        if not st or "#" in ln or st.startswith(("@", "def ")):
            out.append(ln)
            continue
        txt = None
        for pat, msg in _KWARG:
            if re.match(pat, st):
                txt = msg
                break
        if txt is None and inside and not re.match(r"(for|if)\b", st):
            txt = "(continued from the line above)"
        for pat, msg in ([] if txt else _AUTO):
            m = re.search(pat, st)
            if not m:
                continue
            if msg is None and pat.startswith(r"\.(plot"):
                txt = f"curve: {m.group(2)}"
            elif msg is None:
                txt = f"`ch16.{m.group(2)}` (the chapter module)"
            else:
                txt = msg
            break
        if txt is None:                      # no generic placeholder: a line we cannot name precisely is commented by hand
            _UNCOMMENTED.append(st)
        out.append(ln + ("   # " + txt if txt else ""))
    return "\n".join(out)


_UNCOMMENTED: list[str] = []


for _c in nb.cells:
    if _c.cell_type == "code" and not ({"setup", "derivation-check", "primer"} & set(_c.metadata.get("tags", []))):
        _c.source = _auto_comment(_c.source)
if _UNCOMMENTED:
    raise SystemExit("code lines without a comment (write one by hand that says what the line means):\n  "
                     + "\n  ".join(_UNCOMMENTED))

path = nb.save()
print(f"wrote {path.relative_to(ROOT)}: {len(nb.cells)} cells "
      f"({sum(c.cell_type == 'markdown' for c in nb.cells)} markdown, {sum(c.cell_type == 'code' for c in nb.cells)} code), "
      f"{len(nb.cores)} CORE blocks, {len(nb.derivations)} derivations ({sum(d['steps'] for d in nb.derivations.values())} steps), "
      f"{len(nb.recaps)} recaps, {len(nb.explainers)} explainers, {len([p for p in nb.primers if not p.endswith('(reminder)')])} primers")

