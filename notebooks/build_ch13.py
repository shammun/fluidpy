"""Build the Chapter 13 teaching notebook: ``notebooks/ch13_geophysical_fluid_dynamics.ipynb``.

Source of truth: ``analysis/ch13_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract —
``fluidpy.ch13_geophysical_fluid_dynamics`` imported as ``ch13`` re-exports ``core.gfd`` = ``GFD``, ``core.vertical_modes``
= ``VM`` and ``core.shallow_water`` = ``SW``), Part D (runtime budget and FAST sizes), Part E (prerequisite ledger →
primers P307–P335 and one-line reminders of earlier primers), Part F (the 29 derivations D01–D29, one move per step) and
``analysis/ch13_curation.md`` (IDs, depths, section coverage), with the rulings of ``reports/ch13_verification.md`` and of
the orchestrator written into the text:

* slips (false as printed, ``ch13.book_slips()``; rows of kind "loose" are taught as loose statements) are taught in corrected form with a ``slip #k`` box; traps (true
  as printed, easy to misread, ``ch13.traps()``) are ``Common confusion`` callouts; the two are kept apart;
* §13.13, westward flow over a step: taught as slip #14, with the caveat that the description is right for a finite ridge;
* vertical modes are orthogonal with weight 1 for BOTH lids; the surface term belongs to the energy relation only;
  a rigid-lid ``Modes`` starts at the first baroclinic mode (index 0);
* the bottom Ekman layer's largest u is 1.067 U at z = 3πδ/4;
* ``shallow_water_omega`` returns NaN where ``shallow_water_discriminant`` < 0 — taught as a caveat, never printed bare;
* the cached turbulence runs are 64²: qualitative statements only, no −3 slope read off them;
* geostrophic adjustment: the inertial-period MEAN of the run is compared with the analytic end state (a snapshot of a
  sharp step carries grid-scale ripples);
* the lapse-rate sign: Kundu's Γ ≡ dT/dz and the meteorological Γ_met ≡ −dT/dz are both shown wherever a lapse rate appears;
* items that are ours and not in the book are labelled so; no paper is cited for them, except the two Eady numbers 0.3098
  and 1.606, which the verifier read first-hand (reference/ch13/SOURCES.md).

**Derivations are read from Part F at build time** (``part_f()`` below): goal, start, plan, tools, assumptions, every step's
*did / tex / why / plain*, result, check, meaning and traps are copied word for word; the nine ★★★ derivations carry the
sympy check cells listed in Part F.

**Equations.** ``EQ`` below is the analyst's page-image transcription (``analysis/ch13.md`` §2; 21 pages re-read by the
designer), in the CORRECTED form wherever the book prints a slip; a number stands only beside the equation printed under it
(the sin²/cos² form of the inertia–gravity relation is unnumbered and never written under (13.112)).

Book numbers are never printed (rule 9): every worked number uses ``ch13.illustrative_inputs()`` or easy round values.

Run:  .venv/Scripts/python.exe notebooks/build_ch13.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch13.py --dump     (prints the parsed Part F derivations only)
      .venv/Scripts/python.exe notebooks/build_ch13.py --partial  (writes what exists, unchecked — development aid)
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

nb = ChapterNotebook("ch13")
DESIGN = (ROOT / "analysis" / "ch13_design.md").read_text(encoding="utf-8")
_CTRL = sorted({hex(ord(_c)) for _c in DESIGN if ord(_c) < 32 and _c != "\n"})
if _CTRL:                                # a form feed, tab or backspace in the design = a LaTeX backslash eaten by a non-raw string
    _k = next(j for j, _c in enumerate(DESIGN) if ord(_c) < 32 and _c != "\n")
    raise SystemExit(f"analysis/ch13_design.md contains control characters {_CTRL} (line {DESIGN.count(chr(10), 0, _k) + 1}): "
                     f"…{DESIGN[max(0, _k - 40):_k + 20]!r} — repair the eaten backslash (\\f, \\t, \\b, \\r, \\a, \\v) before building")
REMIND: dict[str, str] = {}          # filled below: one-sentence reminders of earlier primers, written for this chapter
NOTES_IN: dict[str, str] = {}        # filled below: notes of the curation taught inside a derivation

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# ---------------------------------------------------------------------------------------------------------------------

EQ: dict[str, str] = {
    "13.1": r"\frac{d\rho_\theta}{dz}=\frac{d\rho}{dz}+\frac{g\rho}{c^2}",
    "13.2": r"\nabla\cdot\mathbf u=0,\ \ \frac{D\mathbf u}{Dt}+2\boldsymbol\Omega\times\mathbf u=-\frac{1}{\rho_0}\nabla p-\frac{g\rho}{\rho_0}\mathbf e_z+\mathbf F,\ \ \frac{D\rho}{Dt}=0",
    "13.3": r"\frac{d\bar p}{dz}=-\bar\rho g",
    "13.4": r"\rho(\mathbf x,t)=\bar\rho(z)+\rho'(\mathbf x,t),\ \ p(\mathbf x,t)=\bar p(z)+p'(\mathbf x,t)",
    "13.5": r"\tau_{xz}=\tau_{zx}=\rho\nu_v\frac{\partial u}{\partial z}+\rho\nu_H\frac{\partial w}{\partial x},\ \ \tau_{yz}=\tau_{zy}=\rho\nu_v\frac{\partial v}{\partial z}+\rho\nu_H\frac{\partial w}{\partial y},\ \ \tau_{xy}=\tau_{yx}=\rho\nu_H\Big(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\Big),\ \ \tau_{xx}=2\rho\nu_H\frac{\partial u}{\partial x},\ \ \tau_{yy}=2\rho\nu_H\frac{\partial v}{\partial y},\ \ \tau_{zz}=2\rho\nu_v\frac{\partial w}{\partial z}",
    "13.6": r"F_x=\nu_H\Big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\Big)+\nu_v\frac{\partial^2u}{\partial z^2},\ \ F_y=\nu_H\Big(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\Big)+\nu_v\frac{\partial^2v}{\partial z^2},\ \ F_z=\nu_H\Big(\frac{\partial^2w}{\partial x^2}+\frac{\partial^2w}{\partial y^2}\Big)+\nu_v\frac{\partial^2w}{\partial z^2}",
    "13.7": r"2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)",
    "13.8": r"f=2\Omega\sin\theta",
    "13.9": r"\frac{Du}{Dt}-fv=-\frac{1}{\rho_0}\frac{\partial p}{\partial x}+\nu_H\Big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\Big)+\nu_v\frac{\partial^2u}{\partial z^2},\ \ \frac{Dv}{Dt}+fu=-\frac{1}{\rho_0}\frac{\partial p}{\partial y}+\nu_H\Big(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\Big)+\nu_v\frac{\partial^2v}{\partial z^2},\ \ \frac{Dw}{Dt}=-\frac{1}{\rho_0}\frac{\partial p}{\partial z}-\frac{g\rho}{\rho_0}+\nu_H\Big(\frac{\partial^2w}{\partial x^2}+\frac{\partial^2w}{\partial y^2}\Big)+\nu_v\frac{\partial^2w}{\partial z^2}",
    "13.10": r"f=f_0+\beta y,\ \ \ \beta\equiv\Big(\frac{df}{dy}\Big)_{\theta_0}=\Big(\frac{df}{d\theta}\frac{d\theta}{dy}\Big)_{\theta_0}=\frac{2\Omega\cos\theta_0}{R}",
    "13.11": r"-fv=-\frac{1}{\rho_0}\frac{\partial p}{\partial x}",
    "13.12": r"fu=-\frac{1}{\rho_0}\frac{\partial p}{\partial y}",
    "13.13": r"\mathrm{Ro}=\frac{U^2/L}{fU}=\frac{U}{fL}",
    "13.14": r"0=-\frac{\partial p}{\partial z}-g\rho",
    "13.15": r"\frac{\partial v}{\partial z}=-\frac{g}{\rho_0f}\frac{\partial\rho}{\partial x},\ \ \frac{\partial u}{\partial z}=\frac{g}{\rho_0f}\frac{\partial\rho}{\partial y}",
    "13.16": r"-2\Omega v=-\frac{1}{\rho}\frac{\partial p}{\partial x}",
    "13.17": r"2\Omega u=-\frac{1}{\rho}\frac{\partial p}{\partial y}",
    "13.18": r"E=\frac{\rho\nu U/L^2}{\rho fU}=\frac{\nu}{fL^2}",
    "13.19": r"\frac{\partial w}{\partial z}=0",
    "13.20": r"\frac{\partial v}{\partial z}=\frac{\partial u}{\partial z}=0",
    "13.21": r"\partial\mathbf u/\partial z=0",
    "13.22": r"-fv=\nu_v\frac{d^2u}{dz^2}",
    "13.23": r"fu=\nu_v\frac{d^2v}{dz^2}",
    "13.24": r"\rho\nu_v\frac{du}{dz}=\tau\ \text{at}\ z=0",
    "13.25": r"\frac{dv}{dz}=0\ \text{at}\ z=0",
    "13.26": r"u,v\to0\ \text{as}\ z\to-\infty",
    "13.27": r"\frac{d^2V}{dz^2}=\frac{if}{\nu_v}V,\ \ V\equiv u+iv",
    "13.28": r"V=A\,e^{(1+i)z/\delta}+B\,e^{-(1+i)z/\delta}",
    "13.29": r"\delta=\sqrt{2\nu_v/f}",
    "13.30": r"\int_{-\infty}^{0}u\,dz=0,\ \ \int_{-\infty}^{0}v\,dz=-\frac{\tau}{\rho f}",
    "13.31": r"-f\frac{dv}{dz}=\nu_v\frac{d^2\omega_y}{dz^2},\ \ -f\frac{du}{dz}=\nu_v\frac{d^2\omega_x}{dz^2}",
    "13.32": r"fU=-\frac{1}{\rho}\frac{dp}{dy}",
    "13.33": r"-fv=\nu_v\frac{d^2u}{dz^2}",
    "13.34": r"fu=\nu_v\frac{d^2v}{dz^2}+fU",
    "13.35": r"u=U,\ v=0\ \text{as}\ z\to\infty",
    "13.36": r"u=0,\ v=0\ \text{at}\ z=0",
    "13.37": r"\frac{d^2V}{dz^2}=\frac{if}{\nu_v}(V-U)",
    "13.38": r"V=U\ \text{as}\ z\to\infty",
    "13.39": r"V=0\ \text{at}\ z=0",
    "13.40": r"V=A\,e^{-(1+i)z/\delta}+B\,e^{(1+i)z/\delta}+U",
    "13.41": r"u=U\big[1-e^{-z/\delta}\cos(z/\delta)\big],\ \ v=Ue^{-z/\delta}\sin(z/\delta)",
    "13.42": r"\frac{\partial p}{\partial x}=\rho g\frac{\partial\eta}{\partial x},\ \ \frac{\partial p}{\partial y}=\rho g\frac{\partial\eta}{\partial y}",
    "13.43": r"(H+\eta)\frac{\partial u}{\partial x}+(H+\eta)\frac{\partial v}{\partial y}+w(\eta)-w(0)=0",
    "13.44": r"\frac{\partial\eta}{\partial t}+\frac{\partial}{\partial x}\big[u(H+\eta)\big]+\frac{\partial}{\partial y}\big[v(H+\eta)\big]=0",
    "13.45": r"\frac{\partial\eta}{\partial t}+H\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)=0,\ \ \frac{\partial u}{\partial t}-fv=-g\frac{\partial\eta}{\partial x},\ \ \frac{\partial v}{\partial t}+fu=-g\frac{\partial\eta}{\partial y}",
    "13.46": r"c^2=gH_e",
    "13.47": r"\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}=0",
    "13.48": r"\frac{\partial u}{\partial t}-fv=-\frac{1}{\rho_0}\frac{\partial p}{\partial x}",
    "13.49": r"\ \frac{\partial v}{\partial t}+fu=-\frac{1}{\rho_0}\frac{\partial p}{\partial y}",
    "13.50": r"0=-\frac{\partial p}{\partial z}-g\rho",
    "13.51": r"\ \frac{\partial\rho}{\partial t}-\frac{\rho_0N^2}{g}w=0",
    "13.52": r"[u,v,p/\rho_0]=\sum_{n=0}^{\infty}[u_n,v_n,p_n]\,\psi_n(z)",
    "13.53": r"w=\sum_{n=0}^{\infty}w_n\int_{-H}^{z}\psi_n(z)\,dz",
    "13.54": r"\rho=\sum_{n=0}^{\infty}\rho_n\frac{d\psi_n}{dz}",
    "13.55": r"\frac{d\psi_n/dz}{N^2\int_{-H}^{z}\psi_n\,dz}=\frac{\rho_0}{g}\frac{w_n}{\partial\rho_n/\partial t}\equiv-\frac{1}{c_n^2}",
    "13.56": r"\frac{d}{dz}\Big(\frac{1}{N^2}\frac{d\psi_n}{dz}\Big)+\frac{1}{c_n^2}\psi_n=0",
    "13.57": r"\frac{\partial u_n}{\partial x}+\frac{\partial v_n}{\partial y}+\frac{1}{c_n^2}\frac{\partial p_n}{\partial t}=0",
    "13.58": r"\frac{\partial u_n}{\partial t}-fv_n=-\frac{\partial p_n}{\partial x}",
    "13.59": r"\ \frac{\partial v_n}{\partial t}+fu_n=-\frac{\partial p_n}{\partial y}",
    "13.60": r"p_n=-\frac{g}{\rho_0}\rho_n",
    "13.61": r"\ w_n=\frac{1}{c_n^2}\frac{\partial p_n}{\partial t}",
    "13.62": r"c_n^2\equiv gH_e",
    "13.63": r"w=\frac{g(\partial\rho/\partial t)}{\rho_0N^2}=-\frac{1}{\rho_0N^2}\frac{\partial^2p}{\partial z\,\partial t}=-\frac{1}{N^2}\sum_{n=0}^{\infty}\frac{\partial p_n}{\partial t}\frac{d\psi_n}{dz}",
    "13.64": r"\frac{d\psi_n}{dz}=0\ \text{at}\ z=-H",
    "13.65": r"\frac{d\psi_n}{dz}+\frac{N^2}{g}\psi_n=0\ \text{at}\ z=0",
    "13.66": r"\frac{d^2\psi_n}{dz^2}+\frac{N^2}{c_n^2}\psi_n=0",
    "13.67": r"\psi_n=A_n\cos\frac{Nz}{c_n}+B_n\sin\frac{Nz}{c_n}",
    "13.68": r"B_n=-\frac{c_nN}{g}A_n",
    "13.69": r"\tan\frac{NH}{c_n}=\frac{c_nN}{g}",
    "13.70": r"c_0=\sqrt{gH}",
    "13.71": r"c_n=\frac{NH}{n\pi},\ n=1,2,3,\dots",
    "13.72": r"\frac{\partial^2u}{\partial t^2}-f\frac{\partial v}{\partial t}=gH\frac{\partial}{\partial x}\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)",
    "13.73": r"\ \frac{\partial^2v}{\partial t^2}+f\frac{\partial u}{\partial t}=gH\frac{\partial}{\partial y}\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)",
    "13.74": r"\frac{\partial^3v}{\partial t^3}+f\Big[f\frac{\partial v}{\partial t}+gH\frac{\partial}{\partial x}\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)\Big]=gH\frac{\partial^2}{\partial y\,\partial t}\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)",
    "13.75": r"\frac{\partial^3v}{\partial t^3}-gH\frac{\partial}{\partial t}\nabla_H^2v+f_0^2\frac{\partial v}{\partial t}-gH\beta\frac{\partial v}{\partial x}=0",
    "13.76": r"\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0",
    "13.77": r"-i\omega\hat u-f\hat v=-ikg\hat\eta",
    "13.78": r"\ -i\omega\hat v+f\hat u=-ilg\hat\eta",
    "13.79": r"\ -i\omega\hat\eta+iH(k\hat u+l\hat v)=0",
    "13.80": r"\hat u=\frac{g\hat\eta}{\omega^2-f^2}(\omega k+ifl),\ \ \hat v=\frac{g\hat\eta}{\omega^2-f^2}(-ifk+\omega l)",
    "13.81": r"\omega^2-f^2=gH(k^2+l^2)",
    "13.82": r"\omega^2=f^2+gHK^2,\ \ \ K=\sqrt{k^2+l^2}",
    "13.83": r"u=\frac{\omega\hat\eta}{kH}\cos(kx-\omega t),\ \ v=\frac{f\hat\eta}{kH}\sin(kx-\omega t)",
    "13.84": r"\frac{\partial\eta}{\partial t}+H\frac{\partial u}{\partial x}=0,\ \ \frac{\partial u}{\partial t}=-g\frac{\partial\eta}{\partial x},\ \ fu=-g\frac{\partial\eta}{\partial y}",
    "13.85": r"-i\omega\hat\eta+iHk\hat u=0,\ \ -i\omega\hat u=-igk\hat\eta,\ \ f\hat u=-g\frac{d\hat\eta}{dy}",
    "13.86": r"c=\sqrt{gH}",
    "13.87": r"\eta=\eta_0e^{-fy/c}\cos k(x-ct),\ \ u=\eta_0\sqrt{\frac{g}{H}}\,e^{-fy/c}\cos k(x-ct)",
    "13.88": r"\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}-fv=-g\frac{\partial\eta}{\partial x}",
    "13.89": r"\frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}+fu=-g\frac{\partial\eta}{\partial y}",
    "13.90": r"\frac{\partial h}{\partial t}+\frac{\partial}{\partial x}(uh)+\frac{\partial}{\partial y}(vh)=0",
    "13.91": r"\frac{\partial}{\partial t}\Big(\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}\Big)+\frac{\partial}{\partial x}\Big[u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}\Big]-\frac{\partial}{\partial y}\Big[u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}\Big]+f_0\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)+\beta v=0",
    "13.92": r"\frac{D\zeta}{Dt}+(\zeta+f_0)\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)+\beta v=0",
    "13.93": r"\frac{D(\zeta+f)}{Dt}=\frac{\zeta+f_0}{h}\frac{Dh}{Dt}",
    "13.94": r"\frac{D}{Dt}\Big(\frac{\zeta+f}{h}\Big)=0",
    "13.95": r"\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}=0,\ \ \frac{\partial u}{\partial t}-fv=-\frac{1}{\rho_0}\frac{\partial p}{\partial x},\ \ \frac{\partial v}{\partial t}+fu=-\frac{1}{\rho_0}\frac{\partial p}{\partial y},\ \ \frac{\partial w}{\partial t}=-\frac{1}{\rho_0}\frac{\partial p}{\partial z}-\frac{\rho g}{\rho_0},\ \ \frac{\partial\rho}{\partial t}-\frac{\rho_0N^2}{g}w=0",
    "13.96": r"\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w+f^2\frac{\partial^2w}{\partial z^2}=0",
    "13.97": r"[u,v,w]=[\hat u(z),\hat v(z),\hat w(z)]\,e^{i(kx+ly-\omega t)}",
    "13.98": r"\frac{d^2\hat w}{dz^2}+\frac{(N^2-\omega^2)(k^2+l^2)}{\omega^2-f^2}\hat w=0",
    "13.99": r"m^2(z)\equiv\frac{(k^2+l^2)\,[N^2(z)-\omega^2]}{\omega^2-f^2}",
    "13.100": r"\frac{d^2\hat w}{dz^2}+m^2\hat w=0",
    "13.101": r"\frac{d^2A}{dz^2}+A\Big[m^2-\Big(\frac{d\phi}{dz}\Big)^2\Big]=0",
    "13.102": r"\ 2\frac{dA}{dz}\frac{d\phi}{dz}+A\frac{d^2\phi}{dz^2}=0",
    "13.103": r"\frac{d\phi}{dz}=\pm m,\ \ \ \phi=\pm\int^{z}m\,dz",
    "13.104": r"\hat w=\frac{A_0}{\sqrt m}\,e^{\pm i\int^{z}m\,dz}",
    "13.105": r"ik\hat u+\frac{d\hat w}{dz}=0",
    "13.106": r"\hat u=\mp\frac{A_0\sqrt m}{k}\,e^{\pm i\int^{z}m\,dz}",
    "13.107": r"\hat v=\pm\frac{if}{\omega}\frac{A_0\sqrt m}{k}\,e^{\pm i\int^{z}m\,dz}",
    "13.108": r"u=\mp\frac{A_0\sqrt m}{k}\cos\Big(kx\pm\int^{z}m\,dz-\omega t\Big),\ \ v=\mp\frac{A_0f\sqrt m}{\omega k}\sin\Big(kx\pm\int^{z}m\,dz-\omega t\Big),\ \ w=\frac{A_0}{\sqrt m}\cos\Big(kx\pm\int^{z}m\,dz-\omega t\Big)",
    "13.109": r"m^2=\frac{k^2(N^2-\omega^2)}{\omega^2-f^2}",
    "13.110": r"u=\mp\cos\omega t,\ \ v=\pm\frac{f}{\omega}\sin\omega t",
    "13.111": r"\frac{u}{w}=\mp\frac{m}{k}=\mp\tan\theta",
    "13.112": r"\omega^2-f^2=\frac{k^2}{m^2}(N^2-\omega^2)",
    "13.113": r"\omega^2=\frac{N^2k^2}{m^2+k^2}",
    "13.114": r"(H+\eta)\Big(\frac{\partial\zeta}{\partial t}+u\frac{\partial\zeta}{\partial x}+v\frac{\partial\zeta}{\partial y}+\beta v\Big)-(\zeta+f_0)\Big(\frac{\partial\eta}{\partial t}+u\frac{\partial\eta}{\partial x}+v\frac{\partial\eta}{\partial y}\Big)=0",
    "13.115": r"H\frac{\partial\zeta}{\partial t}+H\beta v-f_0\frac{\partial\eta}{\partial t}=0",
    "13.116": r"u\simeq-\frac{g}{f_0}\frac{\partial\eta}{\partial y},\ \ v\simeq\frac{g}{f_0}\frac{\partial\eta}{\partial x}",
    "13.117": r"\frac{\partial}{\partial t}\Big(\frac{\partial^2\eta}{\partial x^2}+\frac{\partial^2\eta}{\partial y^2}-\frac{f_0^2}{c^2}\eta\Big)+\beta\frac{\partial\eta}{\partial x}=0",
    "13.118": r"\omega=-\frac{\beta k}{k^2+l^2+f_0^2/c^2}",
    "13.119": r"c_x=\frac{\omega}{k}=-\frac{\beta}{k^2+l^2+f_0^2/c^2}",
    "13.120": r"c_x=U-\frac{\beta}{k^2+l^2+f_0^2/c^2}",
    "13.121": r"\omega^3-c^2\omega(k^2+l^2)-f_0^2\omega-c^2\beta k=0",
    "13.122": r"\Big(\frac{\partial}{\partial t}+\mathbf u\cdot\nabla\Big)(\zeta+f)=0",
    "13.123": r"\frac{\partial}{\partial t}(\nabla^2\psi)+U\frac{\partial}{\partial x}(\nabla^2\psi)+\Big(\beta-\frac{d^2U}{dy^2}\Big)\frac{\partial\psi}{\partial x}=0",
    "13.124": r"\frac{d}{dy}(\bar\zeta+f)=\beta-\frac{d^2U}{dy^2}",
    "13.125": r"\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}-fv=-\frac{1}{\rho_0}\frac{\partial p}{\partial x},\ \ \frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}+fu=-\frac{1}{\rho_0}\frac{\partial p}{\partial y},\ \ 0=-\frac{\partial p}{\partial z}-\rho g,\ \ \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}=0,\ \ \frac{\partial\rho}{\partial t}+u\frac{\partial\rho}{\partial x}+v\frac{\partial\rho}{\partial y}+w\frac{\partial\rho}{\partial z}=0",
    "13.126": r"u=U(z)+u'(x,y,z),\ \ v=v'(x,y,z),\ \ w=w'(x,y,z),\ \ \rho=\bar\rho(y,z)+\rho'(x,y,z),\ \ p=\bar p(y,z)+p'(x,y,z)",
    "13.127": r"fU=-\frac{1}{\rho_0}\frac{\partial\bar p}{\partial y},\ \ 0=-\frac{\partial\bar p}{\partial z}-\bar\rho g",
    "13.128": r"\frac{dU}{dz}=\frac{g}{f\rho_0}\frac{\partial\bar\rho}{\partial y}",
    "13.129": r"\frac{\partial\zeta}{\partial t}+u\frac{\partial\zeta}{\partial x}+v\frac{\partial\zeta}{\partial y}-(\zeta+f)\frac{\partial w}{\partial z}=0",
    "13.130": r"\frac{\partial\zeta'}{\partial t}+U\frac{\partial\zeta'}{\partial x}-f\frac{\partial w'}{\partial z}=0",
    "13.131": r"u'\simeq-\frac{1}{\rho_0f}\frac{\partial p'}{\partial y},\ \ v'\simeq\frac{1}{\rho_0f}\frac{\partial p'}{\partial x}",
    "13.132": r"\zeta'=\frac{1}{\rho_0f}\nabla_H^2p'",
    "13.133": r"\frac{\partial\rho'}{\partial t}+U\frac{\partial\rho'}{\partial x}+v'\frac{\partial\bar\rho}{\partial y}-\frac{\rho_0N^2w'}{g}=0",
    "13.134": r"0=-\frac{\partial p'}{\partial z}-\rho'g",
    "13.135": r"w'=-\frac{1}{\rho_0N^2}\Big[\Big(\frac{\partial}{\partial t}+U\frac{\partial}{\partial x}\Big)\frac{\partial p'}{\partial z}-\frac{dU}{dz}\frac{\partial p'}{\partial x}\Big]",
    "13.136": r"\Big(\frac{\partial}{\partial t}+U\frac{\partial}{\partial x}\Big)\Big[\nabla_H^2p'+\frac{f^2}{N^2}\frac{\partial^2p'}{\partial z^2}\Big]=0",
    "13.137": r"p'=\hat p(z)\,e^{i(kx+ly-\omega t)}",
    "13.138": r"\frac{d^2\hat p}{dz^2}-\alpha^2\hat p=0",
    "13.139": r"\alpha^2\equiv\frac{N^2}{f^2}(k^2+l^2)",
    "13.140": r"\hat p=A\cosh\alpha\Big(z-\frac H2\Big)+B\sinh\alpha\Big(z-\frac H2\Big)",
    "13.141": r"c=\frac{U_0}{2}\pm\frac{U_0}{\alpha H}\sqrt{\Big(\frac{\alpha H}{2}-\tanh\frac{\alpha H}{2}\Big)\Big(\frac{\alpha H}{2}-\coth\frac{\alpha H}{2}\Big)}",
    "13.142": r"\frac{HN}{f}<\frac{\alpha_cH}{k}",
    "13.143": r"\frac{d}{dt}\int_0^\infty S(K)\,dK=0",
    "13.144": r"\ \frac{d}{dt}\int_0^\infty K^2S(K)\,dK=0",
    "13.145": r"\frac{S_1}{S_2}=\frac{K_2-K_0}{K_0-K_1}\,\frac{K_2+K_0}{K_1+K_0},\ \ \frac{K_1^2S_1}{K_2^2S_2}=\frac{K_1^2}{K_2^2}\,\frac{K_2^2-K_0^2}{K_0^2-K_1^2}",
}

# A multi-part equation is stored with its parts separated by ",\ \ " (eq_rows() stacks them in a display).

for _k, _v in EQ.items():
    assert "\\\\" not in _v and "$" not in _v, (_k, _v)
CHAPTER_PREFIX = r"13\."

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
# Part F of the design, parsed at build time: every derivation copied word for word (the ch07–ch11 pattern; the ch13
# design writes the fields inline — "**Goal:** … · **Start:** … · **Plan:** • … • …" — and the steps as
# "n. **did** … · **tex** … · **why** … · **plain** …").
# ---------------------------------------------------------------------------------------------------------------------
_FIELD = re.compile(r"(?:^- |\s·\s+)\*\*(Goal|Start|Plan|Tools|Assumptions|Steps|Result|Check|sympy check intent|sympy check|What it means|Traps)"
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
        body = re.sub(r"^- \*\*sympy check\*\*.*$", "", body, flags=re.M)   # the sympy check line
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
assert len(PF) == 29, len(PF)
N_STEPS = sum(len(d["steps"]) for d in PF.values())
assert N_STEPS == 253, N_STEPS                              # the 253 steps of Part F, none lost by the parser

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
        nb.primers.append(re.sub(r"[`*$]", "", concept) + " (reminder)")
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
                    if not re.match(r"13\.", m.group(1)):
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
# One-sentence reminders written for the way each earlier tool is used in THIS chapter (a concept missing here falls back
# to the gist stored in knowledge/primers.md). Keys are the Part E concept texts, verbatim.
# ---------------------------------------------------------------------------------------------------------------------
REMIND.update({
    "negate the number, flip the inequality": r"multiplying both sides of an inequality by $-1$ reverses it: $a>b\Leftrightarrow-a<-b$ — the whole difference between the two lapse-rate conventions.",
    "second partial derivatives; mixed partials commute": r"for a smooth field the order of two partial derivatives does not matter, $\frac{\partial^2p}{\partial x\,\partial z}=\frac{\partial^2p}{\partial z\,\partial x}$ (Schwarz's theorem).",
    "partial derivative ∂/∂x": r"$\partial F/\partial x$ is the rate of change of $F$ along $x$ with every other variable held fixed.",
    "latitude θ and the local vertical": r"latitude $\theta$ is the angle between the local vertical and the equatorial plane; the earth's axis makes the angle $\theta$ with the local horizontal.",
    "plotly 3-D figures with a dropdown": "a `go.Figure` holds several precomputed cases; a dropdown button only switches which traces are visible, so it works on the web page without Python.",
    "scaling with two length scales": r"in a thin layer a horizontal derivative is of size $1/L$ and a vertical one $1/H$, with $H\ll L$; keep the two apart when estimating terms.",
    "cross product as a determinant": r"$\mathbf a\times\mathbf b$ is the determinant with rows $(\mathbf e_x,\mathbf e_y,\mathbf e_z)$, $\mathbf a$, $\mathbf b$; it is perpendicular to both vectors.",
    "logarithmic bar chart": "on a logarithmic axis equal steps are equal factors, so bars that differ by powers of ten can be compared in one picture.",
    "first-order Taylor expansion": r"near a point, $F(y)\approx F(0)+F'(0)\,y$: the value plus the slope times the distance; the error grows like $y^2$.",
    "`slider_figure` (plotly)": "`slider_figure(fn, name, values, …)` precomputes one set of curves per slider position, so the slider also works on the web page.",
    "`assert np.allclose` / `np.isclose`": "`np.allclose(a, b)` is True when two arrays agree within a tolerance; inside `assert` it stops the cell if they do not.",
    "sympy engines returning a residual and \"ok\"": "sympy does algebra with symbols; our engines build an identity, return what is left over (`residual`) and `ok = True` when that is zero.",
    "order-of-magnitude scaling": r"replace each derivative by (typical change)/(typical distance) — $\partial u/\partial x\sim U/L$ — to compare the sizes of terms without solving anything.",
    "`np.meshgrid`, the [j, i] grid layout, array slicing for centred differences": "`X, Y = np.meshgrid(x, y)` gives arrays indexed `[j, i]` = (row y, column x); `(p[:, 2:] - p[:, :-2])/(2*dx)` is a centred x-difference on the interior columns.",
    "contour, quiver": "`ax.contour` draws lines of equal value of a field; `ax.quiver` draws an arrow for a vector at each grid point.",
    "inertial loops of a released parcel; `GFD.parcel_adjust`; complex numbers in numpy": "numpy handles complex numbers natively: `1j` is $i$, `z.real` and `z.imag` are the two parts, `abs(z)` the length — one complex number carries a horizontal vector.",
    "`show_viz` explainers": "`show_viz(chapter, slug)` embeds an interactive explainer in the cell output (full window on the web page).",
    "cross-differentiation to eliminate a variable": r"differentiate one equation in $z$ and another in $x$ so that both contain the same mixed derivative, then equate (or subtract) to remove the pressure.",
    "`scipy.integrate.cumulative_trapezoid`": "`cumulative_trapezoid(y, x, initial=0)` returns the running integral of sampled data — here a shear integrated upward into a wind profile.",
    "hypotheses of a theorem (necessary conditions)": "a theorem holds only under its stated hypotheses; drop one (here: uniform density) and the conclusion may fail.",
    "trying e^{rz} in a linear ODE with constant coefficients": r"a linear equation with constant coefficients is solved by trying $e^{rz}$: the derivatives become powers of $r$, and an algebraic equation for $r$ is left.",
    "square root of i; conjugates; polar form; Euler's formula": r"$i=e^{i\pi/2}$, so $\sqrt i=e^{i\pi/4}=(1+i)/\sqrt2$; and $e^{is}=\cos s+i\sin s$ turns a complex exponential into a cosine and a sine.",
    "integral of an exponential over a half-line": r"$\int_{-\infty}^0e^{az}\,dz=1/a$ whenever the real part of $a$ is positive, because the integrand dies away at the far end.",
    "`np.trapezoid`": "`np.trapezoid(y, x)` integrates sampled data by joining the samples with straight lines.",
    "Ekman layer with K_v(z); finite-volume tridiagonal solve; `GFD.ekman_solve`": "a second-order boundary-value problem on a grid couples each node to its two neighbours only, so it is one tridiagonal (banded) linear solve.",
    "particular solution of a forced linear ODE": "the general solution of a forced linear equation is any one solution of the full equation plus the general solution of the unforced one.",
    "product rule": r"$(ab)'=a'b+ab'$.",
    "staggered arrays (η at centres, u at faces)": "with `n` cells there are `n` centre values and `n + 1` face values; `np.diff(u)/dx` of the face array lands exactly on the centres.",
    "separation of variables": r"try a product such as $F(x,y,t)\,\psi(z)$; if the equation splits into a part depending only on $z$ and a part that does not, each must equal the same constant.",
    "eigenvalue problem for a differential operator": "a differential equation with a free parameter and conditions at both ends has non-zero solutions only for special values of the parameter — the eigenvalues.",
    "integration by parts": r"$\int_a^bF\,G'\,dz=[FG]_a^b-\int_a^bF'G\,dz$: move a derivative from one factor to the other at the price of a boundary term.",
    "Dirichlet and Neumann conditions": "a Dirichlet condition fixes the value at a boundary, a Neumann condition fixes the slope.",
    "`scipy.optimize.brentq`": "`brentq(F, a, b)` finds the root of `F` between `a` and `b` when `F(a)` and `F(b)` have opposite signs.",
    "eigenvalues and eigenvectors; generalised symmetric eigenproblem": r"$A\mathbf x=\lambda\mathbf x$: a vector the matrix only stretches; a symmetric matrix has real eigenvalues and orthogonal eigenvectors.",
    "plane wave e^{i(kx+ly−ωt)}; wavenumbers k, l; K": r"a field $\hat v\,e^{i(kx+ly-\omega t)}$ (real part understood) is a wave with crests $2\pi/K$ apart, $K=\sqrt{k^2+l^2}$; $\partial/\partial t\to-i\omega$, $\partial/\partial x\to ik$, $\partial/\partial y\to il$.",
    "`np.roots`": "`np.roots([a, b, c, d])` returns the (possibly complex) roots of $a\\omega^3+b\\omega^2+c\\omega+d$.",
    "dominant balance; three frequency regimes": "when the terms of an equation differ greatly in size, the largest two must balance each other and the smallest may be dropped — then check the dropped one really was small.",
    "phase speed ω/k": r"a crest of $\cos(kx-\omega t)$ moves at $c=\omega/k$.",
    "`np.linalg.solve`": "`np.linalg.solve(M, b)` solves the linear system `M x = b` (complex entries allowed).",
    "chain rule (implicit differentiation of ω²)": r"differentiating $\omega^2=F(k)$ with respect to $k$ gives $2\omega\,\partial\omega/\partial k=F'(k)$ — no square root needed.",
    "a linear map of a circle is an ellipse": r"$(a\cos s,\ b\sin s)$ is a circle stretched by $a$ along one axis and $b$ along the other: an ellipse with those semi-axes.",
    "`animate` and `show_animation`": "`animate(update, frames, fig)` calls `update(i)` once per frame; `show_animation(anim, player=…)` shows it as a video or as a frame player with step buttons.",
    "first-order linear ODE": r"$d\hat\eta/dy=a\,\hat\eta$ has the solution $\hat\eta=\hat\eta(0)\,e^{ay}$.",
    "cached model runs; `ch13.load_reference_run`": "an expensive run is made once, saved as a compressed `.npz` file with its parameters, and loaded where it is shown.",
    "reduced gravity; internal Kelvin waves; internal Rossby radius; three radii": r"across an interface between two layers gravity acts only through the density difference, $g'=g\,\Delta\rho/\rho$ — a hundred to a thousand times weaker than $g$.",
    "frames player for stop-and-look animations": "`show_animation(anim, player=\"frames\")` gives ◀ ▮▮ ▶ and single-step buttons, for processes worth stopping at.",
    "quotient rule": r"$(a/b)'=(a'b-ab')/b^2$.",
    "`scipy.integrate.solve_ivp`": "`solve_ivp(rhs, (z0, z1), y0, t_eval=…)` integrates a system of first-order ordinary differential equations from given starting values.",
    "completing the square into a circle": r"$k^2+ak+l^2=b$ is $(k+a/2)^2+l^2=b+a^2/4$: a circle centred at $(-a/2,0)$.",
    "Fourier modes": r"any periodic field is a sum of waves $e^{ikx}$; a linear constant-coefficient equation evolves each one on its own.",
    "Rayleigh equation with β; normal modes with complex c": r"perturb a steady flow, drop products of small quantities, try $e^{ik(x-ct)}$: the flow is unstable if some mode has $c_i=\mathrm{Im}\,c>0$, because it then grows like $e^{kc_it}$.",
    "real and imaginary parts of a complex equation": "a complex equation is two real ones: its real part and its imaginary part must each hold.",
    "cosh, sinh, tanh": r"$\cosh x=(e^x+e^{-x})/2$, $\sinh x=(e^x-e^{-x})/2$, $\tanh x=\sinh x/\cosh x$; $(\cosh)'=\sinh$, $(\sinh)'=\cosh$.",
    "`scipy.optimize.minimize_scalar`": "`minimize_scalar(F, bounds=(a, b), method=\"bounded\")` finds the minimum of a function of one variable on an interval (maximise by minimising `-F`).",
    "Chebyshev eigen-solve as an independent route": "on Chebyshev points a derivative becomes a dense matrix that is accurate to many digits with few points, so a differential eigenproblem becomes a small matrix one.",
    "mean of a product of two wave fields": r"over a wavelength, the mean of $\mathrm{Re}(\hat a\,e^{i\phi})\,\mathrm{Re}(\hat b\,e^{i\phi})$ is $\tfrac12\mathrm{Re}(\hat a\,\hat b^*)$: only the in-phase parts contribute.",
})


_tnn = _title_no_numbers


def _title_no_numbers(title: str) -> str:   # noqa: F811  (tidies the punctuation left where an equation number was removed)
    s = _tnn(title)
    s = re.sub(r":\s*,", ":", s)
    s = re.sub(r",\s*,", ",", s)
    s = re.sub(r"\s+and\s*;", ";", s)
    s = re.sub(r":\s*and\b", ":", s)
    return re.sub(r"\s{2,}", " ", s).strip(" ,:;")


def _pf_all(key: str, old: str, new: str) -> None:
    """Replace a phrase in every text field of one Part F block (fails loudly if it is not there)."""
    d, n = PF[key], 0
    for fld in ("goal", "assumptions", "check", "meaning", "traps"):
        n += d[fld].count(old)
        d[fld] = d[fld].replace(old, new)
    d["tools"] = [t.replace(old, new) for t in d["tools"]]
    for st in d["steps"]:
        for part in ("did", "tex", "why", "plain"):
            n += st[part].count(old)
            st[part] = st[part].replace(old, new)
    for fld in ("start", "result"):
        n += d[fld][0].count(old) + d[fld][1].count(old)
        d[fld] = (d[fld][0].replace(old, new), d[fld][1].replace(old, new))
    assert n or any(new in t for t in d["tools"]), (key, old)


# reader-facing wording of Part F: no pipeline language, no quoted book phrases, a few precise words
_pf_all("D02", r"\ll1", r"\approx0.1") if any(r"\mathrm{m\,s^{-2}}}\ll1" in st["tex"] for st in PF["D02"]["steps"]) else None
PF["D02"]["steps"][5]["why"] = PF["D02"]["steps"][5]["why"].replace("So it is dropped.", "A tenth is small rather than negligible (against the full weight g the term is 10⁻⁴). So it is dropped.")
_w9 = PF["D09"]["steps"][5]["why"]
PF["D09"]["steps"][5]["why"] = (r"From D08, with the complex geostrophic velocity written $W_g\equiv U_g+iV_g$, "
                                r"$\int(V-W_g)\,dz=-W_g\delta/(1+i)=\tfrac\delta2(-1+i)W_g$" + _w9[_w9.index("; only this"):])
PF["D24"]["check"] = PF["D24"]["check"].replace(r"\mathrm{sech}^2(y/L)$", r"\mathrm{sech}^2(y/L)$ (where $\mathrm{sech}\,y\equiv1/\cosh y$, a smooth bump)", 1)
assert "equiv1/\\cosh" in PF["D24"]["check"]
PF["D11"]["tools"] = PF["D11"]["tools"] + ["the bottom and free-surface conditions on the modes (notes N66 and N67 above) and the Robin condition (P315)"]
_pf_all("D18", "matching at a junction (P323)", "matching two solutions at a junction (the value and the slope must be continuous there)")
_pf_all("D18", " (reported to the verifier)", "")
PF["D11"]["check"] = re.sub(r"; verified by the designer's own finite-volume solve of the\s+thermocline profile with a free surface: off-diagonal overlaps below 10⁻⁸\)",
                            "; the code cell further down measures it for free-surface modes)", PF["D11"]["check"])
assert "designer" not in PF["D11"]["check"], PF["D11"]["check"][:300]
PF["D16"]["check"] = re.sub(r"Model — our forward–backward march averaged over the fifth inertial period matches.*?(?=$)",
                            "Model — our forward–backward march, averaged over one inertial period, matches the closed form to a fraction of a per cent (the animation cell below prints the measured numbers).",
                            PF["D16"]["check"], flags=re.S)
PF["D26"]["goal"] = PF["D26"]["goal"].replace('The book says "after some straightforward algebra".', "The book leaves the algebra to the reader.")
PF["D27"]["goal"] = PF["D27"]["goal"].replace('The book says "it can be shown" and stops at the wavelength', "The book states the result without the steps and stops at the wavelength")
for _k in ("D26", "D27"):
    assert "The book says" not in PF[_k]["goal"], PF[_k]["goal"]
PF["D29"]["traps"] = re.sub(r"^α is the enstrophy flux", "The letter α is the enstrophy flux", PF["D29"]["traps"])
PF["D18"]["start"] = (PF["D18"]["start"][0], "Each column keeps its potential vorticity; upstream the flow is uniform")
PF["D28"]["start"] = (PF["D28"]["start"][0], "Energy and enstrophy are both conserved while energy moves between three wavenumbers")
PF["D01"]["result"] = (PF["D01"]["result"][0], PF["D01"]["result"][1] or "Friction is diffusion of momentum: fast sideways, slow vertically")

D06_EXTRA = r"""
# step 14: the southern hemisphere, f = -|f| (the same symbol f now stands for the magnitude |f|)
lam_S = (1 - sp.I) / delta                                    # the root that decays downward when f < 0
A_S = tau * delta * (1 + sp.I) / (2 * rho * nu)               # the constant from the same stress condition
V_S = A_S * sp.exp(lam_S * z)                                 # the southern solution
assert sp.simplify(sp.diff(V_S, z, 2) + sp.I * f / nu * V_S) == 0     # it solves V'' = (i f/nu) V with f -> -|f|
assert sp.simplify(rho * nu * sp.diff(V_S, z).subs(z, 0) - tau) == 0  # the same surface stress, along x
assert sp.simplify(sp.im(V_S.subs(z, 0))) > 0                 # surface current to the LEFT of the wind (v > 0 for an eastward wind)
print("D06 step 14 passed: for f < 0 the spiral is the mirror image (surface current 45 degrees to the left)")
"""
D11_EXTRA = r"""
# steps 11-12 for ANY two modes and ANY N(z): the identity behind orthogonality
pm, pn, N2f = sp.Function("psi_m")(z), sp.Function("psi_n")(z), sp.Function("N2")(z)   # two modes and N^2(z), left general
flux = (pm * sp.diff(pn, z) - pn * sp.diff(pm, z)) / N2f      # the bracket that is evaluated at the two boundaries
lhs = sp.diff(flux, z)                                        # its derivative …
rhs = pm * sp.diff(sp.diff(pn, z) / N2f, z) - pn * sp.diff(sp.diff(pm, z) / N2f, z)    # … equals psi_m (psi_n'/N^2)' - psi_n (psi_m'/N^2)'
assert sp.simplify(lhs - rhs) == 0                            # the cross terms psi_m' psi_n'/N^2 cancel
print("D11 steps 11-12 passed for general modes; the free-surface case is measured numerically in the orthogonality cell below")
"""
PF["D06"]["check_src"] = PF["D06"]["check_src"].rstrip("\n") + "\n" + D06_EXTRA.strip("\n")
PF["D11"]["check_src"] = PF["D11"]["check_src"].rstrip("\n") + "\n" + D11_EXTRA.strip("\n")


def trap(tid: str, text: str) -> None:
    """A trap of ``ch13.traps()`` (true as printed, easy to misread) as a callout."""
    nb.md(f"> ⚠️ **Common confusion (trap {tid}):** {textwrap.dedent(text).strip()}")


def choice(text: str) -> None:
    """A numerical or modelling choice that is ours, not the book's."""
    nb.md(f"> 🔧 **Our choice** — {textwrap.dedent(text).strip()}")


def loose(k: int, quoted: str, formula: str) -> None:
    """A row of the slip table of kind "loose": an inconsistent order-of-magnitude statement, not an error."""
    nb.md(f"> ⚠️ **slip table, row #{k} (a loose order-of-magnitude statement, not an error) — the book quotes** "
          f"{textwrap.dedent(quoted).strip()} **; what the formula gives is** {textwrap.dedent(formula).strip()}")


for _st in PF["D17"]["steps"]:                    # the word "partial" inside \text{…} trips the eaten-backslash check
    _st["tex"] = _st["tex"].replace("subscripts = partial derivatives", "subscripts mean derivatives")


def ours(text: str) -> None:
    """A result that is not in the book."""
    nb.md(f"> 🧭 **Not in the book — ours.** {textwrap.dedent(text).strip()}")


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
        nb.primers.append(re.sub(r"[`*$]", "", concept) + " (reminder)")
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
Look at a weather map: the wind does not blow from high pressure to low, it blows *around* the highs and lows. Look at a map of
ocean currents: the water piles up in the middle of each basin and circles it. Two things make the atmosphere and the ocean
behave so unlike water in a sink — the earth turns, and the fluid is layered, light over heavy. This chapter adds those two
ingredients to everything the book has built. Rotation first: on a thin shell only the local vertical part of the earth's spin
matters, and slow flows stay close to a stand-off between the Coriolis force and the pressure gradient (geostrophy), with a
vertical shear wherever temperature changes horizontally (thermal wind) and thin friction layers at the sea surface and the
ground (Ekman layers) that move water at right angles to the wind. Then the waves that a rotating layer supports — fast ones
that cannot go slower than the earth's own turning rate, a coast-hugging one, and one slow planetary wave that exists only
because the Coriolis parameter changes with latitude. One conserved quantity, potential vorticity, organises all the slow
motion. Finally, why weather exists at all (baroclinic instability) and why its eddies grow into jets instead of breaking down
(two-dimensional turbulence).""",
    roadmap=[
        "**The equations on a thin rotating shell** — only the vertical part of the earth's spin turns the wind (`C01`).",
        "**Geostrophic balance** — Coriolis force against pressure gradient (`C02`).",
        "**Thermal wind** — a horizontal temperature gradient sets the change of wind with height (`C03`).",
        "**The surface Ekman spiral** — friction against rotation makes a layer of fixed thickness (`C04`).",
        "**Ekman transport and pumping** — the water goes at right angles to the wind (`C05`).",
        "**The bottom Ekman layer** — why air spirals into a low (`C06`).",
        "**The shallow-water equations** — three equations for a thin layer (`C07`).",
        "**Vertical normal modes** — a stratified ocean as a stack of shallow layers (`C08`).",
        "**One cubic for all the waves** — two fast roots, one slow (`C09`).",
        "**Poincaré waves** — gravity waves with a floor under their frequency (`C10`).",
        "**The Kelvin wave** — the wave that leans on a coast (`C11`).",
        "**The Rossby radius and geostrophic adjustment** — what rotation holds up (`C12`).",
        "**Potential vorticity** — what every column keeps (`C13`).",
        "**Inertia–gravity waves** — internal waves between $f$ and $N$ (`C14`).",
        "**Rossby waves (and barotropic instability)** — crests west, energy either way (`C15`).",
        "**The Eady problem** — where storms come from (`C16`).",
        "**Two conserved quantities and the inverse cascade** — why eddies merge (`C17`).",
    ],
    prerequisites=[
        r"the rotating-frame momentum equation $\rho\big(\frac{D'\mathbf u'}{Dt}\big)_{O'1'2'3'}=-\nabla'p+\rho\big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\big]+\mu\nabla'^2\mathbf u'$ (4.45) — for the steadily turning earth only the Coriolis and centrifugal terms of the bracket survive — and the Boussinesq approximation (Ch. 4)",
        "static stability, potential temperature and the two lapse-rate conventions (Ch. 1)",
        "vorticity, stretching and the column argument (Ch. 5)",
        "surface and internal gravity waves, group velocity (Ch. 7); the Stokes layer (Ch. 8); boundary-layer scaling (Ch. 9)",
        "Rayleigh's criterion and normal modes with a complex phase speed (Ch. 11); eddy viscosity and the energy cascade (Ch. 12)",
    ],
)
nb.explainer_index([
    ("geostrophic_balance", "Why doesn't air flow straight from high to low pressure?", "the Coriolis force keeps turning the moving air, so that on average it runs along the isobars; only with drag does it settle"),
    ("thermal_wind", "Why is there a jet stream, and why is it westerly in both hemispheres?", "a horizontal temperature gradient sets the change of wind with height"),
    ("ekman_spiral", "The wind blows east — why does the water go south?", "the eddy viscosity sets how deep the spiral reaches, not how much water it carries"),
    ("ekman_force_balance", "Why does air spiral into a low?", "friction weakens the Coriolis force near the ground and the pressure force wins"),
    ("vertical_modes", "How can one layer stand for a stratified ocean?", "each vertical mode is a shallow-water system of its own, much shallower, depth"),
    ("shallow_water_dispersion", "Are Poincaré, Kelvin and Rossby waves separate theories?", "three roots of one cubic, decades apart in frequency"),
    ("kelvin_wave", "Why does this wave run only one way along a coast?", "only one direction gives a slope that decays offshore"),
    ("geostrophic_adjustment", "Why doesn't a pile of water flatten out on a rotating planet?", "rotation holds up whatever lies beyond a Rossby radius"),
    ("rossby_waves", "The crests go west — so how does the energy go east?", "phase and group velocity part company at the Rossby radius"),
    ("eady_instability", "Where do storms come from?", "two boundary waves that hold each other in place and grow, feeding on the tilt of the density surfaces"),
])
nb.setup()

nb.md(r"""
## 📦 Imports and our own inputs

One cell loads everything the chapter uses. The physics lives in three new core modules — `GFD` (rotation, balances, Ekman
layers, waves, instability), `VM` (vertical modes) and `SW` (our numerical shallow-water and vorticity models) — and the
chapter module `ch13` re-exports all of them, so `ch13.rossby_omega` and `GFD.rossby_omega` are the same function.""")
code(r"""
import logging                                                 # to silence one harmless library message
import numpy as np                                             # arrays and maths
import sympy as sp                                             # symbolic algebra (derivation checks)
import matplotlib.pyplot as plt                                # static figures
import pandas as pd                                            # small tables
import plotly.graph_objects as go                              # interactive 3-D figures
from scipy.integrate import solve_ivp, cumulative_trapezoid    # ODE integration and running integrals
from scipy.optimize import brentq, minimize_scalar             # root finding and one-variable minimisation
from scipy.linalg import eigh_tridiagonal                      # eigenvalues of a symmetric tridiagonal matrix
from fluidpy import ch13_geophysical_fluid_dynamics as ch13    # the tested chapter module (re-exports GFD, VM, SW)
from fluidpy.core import gfd as GFD, vertical_modes as VM, shallow_water as SW   # the three new core modules
from fluidpy.core import rotating as ROT, stratification as STRAT, similarity as SIM   # Ch. 4 and Ch. 1 primitives
from fluidpy.core import waves as WAV, stability as ST, laminar as LAM, turbstats as TS, dimensional as DIM   # Ch. 7, 11, 8, 12, 1
from fluidpy import ch05_vorticity_dynamics as ch05, ch11_instability as ch11, ch12_turbulence as ch12   # earlier chapter modules
from fluidpy.core.interact import slider_figure, animate_figure, live   # plotly sliders (work on the page) and live widgets
from fluidpy.core.anim import animate                          # builds a matplotlib animation
from fluidpy.core.style import COLORS, savefig                 # the shared palette and a figure saver
from fluidpy.core.thermo import G0                             # standard gravity 9.80665 m/s^2
sys.path.insert(0, "scripts")                                  # the folder with our sketch helpers
from ch13_drawings import *                                    # draw_tangent_plane, draw_taylor_column, ... (pure matplotlib)

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # no font-substitution messages under figures
plt.rcParams["figure.dpi"] = min(plt.rcParams["figure.dpi"], 80)       # modest figure resolution keeps the notebook and the page light
inp = ch13.illustrative_inputs()                               # OUR illustrative inputs (never the book's numbers)
lat = inp["lat"]                                               # 35 degrees north, in radians
f35 = GFD.coriolis_parameter(lat)                              # Coriolis parameter f = 2 Omega sin(latitude) at 35 N [1/s]
f60 = GFD.coriolis_parameter(inp["lat_ekman"])                 # the same at 60 N (used for the Ekman layers)
f12 = GFD.coriolis_parameter(inp["lat_rossby"])                # the same at 12 N (used for Rossby waves)
beta35 = GFD.beta_parameter(lat)                               # beta = df/dy at 35 N [1/(m s)]
beta12 = GFD.beta_parameter(inp["lat_rossby"])                 # beta at 12 N
Gamma_a = STRAT.adiabatic_lapse_rate()                         # adiabatic temperature gradient dT/dz = -g/c_p [K/m] (Kundu sign)
print(f"f(35N) = {f35:.4e} 1/s   f(60N) = {f60:.4e} 1/s   f(12N) = {f12:.4e} 1/s")   # the three Coriolis parameters
print(f"beta(35N) = {beta35:.4e} 1/(m s)   beta(12N) = {beta12:.4e} 1/(m s)")        # and the two beta values
print(f"Gamma_a = {Gamma_a*1e3:.2f} K/km (Kundu)  <=>  Gamma_d = {-Gamma_a*1e3:.2f} K/km (meteorology)")   # both conventions
""", explain=r"""
1. The first block imports numerical, symbolic and plotting tools; each is recalled where it is first used.
2. `GFD`, `VM`, `SW` are the new modules of this chapter; `ROT`, `STRAT`, `SIM`, `WAV`, `ST`, `LAM`, `TS`, `DIM` are primitives
   written for earlier chapters and reused here unchanged.
3. `inp = ch13.illustrative_inputs()` holds **our own** illustrative numbers (35° N and 35° S for most examples, 60° N for the
   Ekman layers, 12° N for Rossby waves; an ocean 4200 m deep, an atmosphere 9 km deep, and so on).
4. `f35`, `f60`, `f12` are the Coriolis parameter $f=2\Omega\sin\theta$ *(13.8)* at three latitudes; `beta35`, `beta12` its
   northward rate of change; `Gamma_a` the adiabatic temperature gradient in Kundu's sign.""")

remind("C01", "§13.1–§13.4")

nb.md(r"""
## ⚠️ Conventions in this chapter (read once, come back often)

**(a) Axes.** The coordinate $x$ points east, $y$ north, $z$ up; the velocity is $(u, v, w)$.

**(b) Where is $z=0$?** It moves from section to section (trap **T4**):

| Section | $z=0$ is | The fluid is in |
|---|---|---|
| §13.6 (surface Ekman layer), §13.9 (vertical modes), §13.14 (internal waves) | the sea surface | $z<0$ |
| §13.7 (bottom Ekman layer) | the solid surface | $z>0$ |
| §13.8 (shallow water) | the flat bottom | $0<z<H+\eta$ |
| §13.17 (Eady problem) | the lower lid | $0<z<H$ |

**(c) Hemisphere.** Every "to the right", "clockwise" and "coast on the right" in the book is the northern-hemisphere case
$f>0$. In this notebook each such sentence says so and adds *mirror for $f<0$* (trap **T5**); the code takes the sign of $f$.

**(d) Primes are dropped.** From the thin-shell equations of `C01` on, $p$ and $\rho$ are *perturbations* from a state of
rest (trap **T3**).

**(e) Overloaded letters** (trap **T13**) — the table printed by the next cell: $f$, $\Omega/\omega$, $\beta$, $N$,
$H/h/\eta$, $c/c_n$, $k/l/m/K$, $\zeta$, $\theta$, $\delta$, $\Lambda/\lambda$, $\psi$, $E/\mathrm{Ro}/R$, $\alpha$, $U/V/\tau$.

**(f) Four more traps, one line each.** Three Rossby radii share one name — external $\sqrt{gH}/f$, internal $NH/(n\pi f)$,
and the Eady radius $NH/f$ without the $\pi$ (trap **T11**; `C12`, `C16`). Two Rossby numbers — $U/(fL)$ here, $U/(2\Omega L)$
in Chapter 4 (trap **T17**; `C02`). The stream function here has $u=-\partial\psi/\partial y$, $v=\partial\psi/\partial x$, the
opposite sign to Chapters 4 and 11 (trap **T14**; `C02`, `C15`). The spectrum of §13.18 is one-sided with no factor ½
(trap **T16**; `C17`).

**(g) Cross-references.**

> ⚠️ **slip #9 — the book prints** cross-references to Section 4.18 (for the Boussinesq set), Section 8.7 (for the
> impulsively started plate), Section 5.7 (for vortex stretching) and to equation number 7.128 **; the correct form is**
> §4.9, §8.4, §5.6 and the definition
> $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ (7.127). We cite our own section numbers and write the equation out.

**(h) Slips and traps are kept apart.** A *slip* is false as printed (each is taught in corrected form in a `slip #k` box;
two rows of the table are only *loose* order-of-magnitude statements and are labelled so); a *trap* is true as printed but
easy to misread (each is a ⚠️ *Common confusion* callout). The next cell counts and lists them.

**(i) Which inputs?** All worked numbers use our own illustrative inputs (`ch13.illustrative_inputs()`) or easy round
values, never the book's.""")
code(r"""
pd.set_option("display.max_colwidth", 90)                      # let long table cells show
display(ch13.conventions_table())                              # the overloaded letters: book meanings, our symbol, code name
""", explain=r"""
**What does this show?** One row per overloaded letter: what it means in this chapter (and meant earlier), the symbol this
notebook uses to keep the meanings apart, and the name in the code. Each block repeats the row it needs.""")
code(r"""
slips = pd.DataFrame(ch13.book_slips()).T                      # the slip table, one row each
display(slips[["kind", "where", "taught_in", "how_to_tell"]])  # slip or loose statement, where it is, which block treats it, how to tell
traps_tab = pd.DataFrame(ch13.traps())                         # the traps (true as printed, easy to misread)
display(traps_tab[["id", "what", "where"]])                    # what is easy to misread and where the callout sits
print(f"{len(slips)} rows in the slip table ({(slips['kind'] == 'loose').sum()} of them loose statements), {len(traps_tab)} traps — kept in two separate tables on purpose")   # the counts
""", explain=r"""
**What does this show?** The first table lists the places where the printed text is wrong (kind "slip") or only loosely
consistent (kind "loose"), and the block that treats each; the second the places where it is right but invites a misreading. Slips are false as printed;
traps are true but misleading — they are never mixed.""")

# =====================================================================================================================
# A.1  §13.1 Introduction
# =====================================================================================================================
nb.section("13.1", "Introduction", intro=r"""
**What is this section about?** What makes the atmosphere and the ocean a subject of their own: the earth's rotation and the
layering of density. And the words — easterly, eastward, cyclonic — that the rest of the chapter uses.""")
note("N01 [C]", r"""
**The picture.** Geophysical fluid dynamics is the dynamics of the atmosphere and the ocean. Two features set it apart:
rotation and vertical density stratification. Their most visible signature is that large-scale flow runs *along* lines of
constant pressure, not across them. Rotation is the subject of `C01`–`C06`; stratification of the §13.2 recap, `C08`, `C14`
and `C16`.""")
fig(r"""
x = np.linspace(-2.0e6, 2.0e6, 41)                             # east-west distance from the centre [m]
y = np.linspace(-2.0e6, 2.0e6, 41)                             # north-south distance [m]
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
for ax, kind in zip(axes, ("low", "high")):                    # one panel for a low, one for a high
    pc = ch13.pressure_centre(x, y, dp=400.0, R_c=6.0e5, lat_rad=lat, kind=kind)   # a 400 Pa = 4 hPa pressure anomaly (1 hPa = 100 Pa = 1 millibar) and its balanced wind at 35 N
    ax.contour(x/1e3, y/1e3, pc["p"]/100, 7, colors=COLORS["orange"], linewidths=1.2)   # isobars [hPa anomaly], orange = pressure
    ax.quiver(x[::4]/1e3, y[::4]/1e3, pc["u"][::4, ::4], pc["v"][::4, ::4], color=COLORS["teal"], scale=70)   # the wind, every fourth point (an arrow as long as the panel is wide would be 70 m/s)
    ax.text(0, 0, "L" if kind == "low" else "H", ha="center", va="center", fontsize=16, weight="bold")   # mark the centre
    ax.set(title=f"35° N, a {kind}: wind circles the centre", xlabel="x east [km]", ylabel="y north [km]", aspect="equal")
plt.show()
""",
    see="Isobars (orange; pressure in hectopascals, 1 hPa = 100 Pa) of a low and of a high at 35° N, and the wind (teal arrows) that balances each of them. The arrows "
        "circle the centres: counter-clockwise round the low, clockwise round the high (northern hemisphere, $f>0$; mirror for $f<0$).",
    read="No arrow points at the low. Each one runs along an isobar, with low pressure on its left. Why that is the only "
         "steady possibility is the subject of `C02`.",
    change="…the same two pressure patterns sat at 35° S? Both circulations would reverse (`C02` shows all four cases).")
P("P307", "geographic vocabulary: zonal / meridional, easterly wind vs eastward current, poleward / equatorward, cyclonic / anticyclonic", r"""
Zonal = along a latitude circle (east–west); meridional = north–south. A wind is named by where it comes FROM (a westerly
blows toward the east); a current by where it goes TO (an eastward current). Poleward and equatorward avoid saying north or
south. Cyclonic = turning in the same sense as the earth below: counter-clockwise in the northern hemisphere, clockwise in the
southern.""", code=r"""
print(ch13.wind_from_to(10.0, 0.0))        # u = +10 m/s: a 'westerly' wind, an 'eastward' current
print(GFD.hemisphere(np.deg2rad(-35.0)))   # 35 S: sign -1, turns 'left', cyclonic = 'clockwise'
""")
note("N02 [B]", r"""
**Conventions in words.** Height $z$ is positive upward. Winds are named by origin, currents by destination. All the book's
figures are drawn for the northern hemisphere; the sign of every rotational effect follows the sign of $f$.""")
code(r"""
rows = [dict(latitude=f"{np.rad2deg(L):+.0f}°", **GFD.hemisphere(L)) for L in (lat, -lat)]   # 35 N and 35 S side by side
display(pd.DataFrame(rows).rename(columns={"sign": "sign of f", "name": "hemisphere", "turns": "moving fluid is turned to the", "cyclonic": "sense round a low"}))   # the two-row table
""", explain=r"""
**What does this show?** The same four facts for the two hemispheres: the sign of $f$, the side to which moving fluid is
deflected, and the sense of the circulation round a low. Everything "right-handed" at 35° N is "left-handed" at 35° S.""")

# =====================================================================================================================
# A.2  §13.2 Vertical variation of density
# =====================================================================================================================
nb.section("13.2", "Vertical Variation of Density in the Atmosphere and Ocean", intro=r"""
**What is this section about?** A reminder of static stability from Chapter 1, in the form this chapter needs: which density
gradient decides stability, how the lapse rate says the same thing in two sign conventions, and the buoyancy frequency $N(z)$
that later sets the speed of internal modes (`C08`), the frequency range of internal waves (`C14`) and the growth of storms
(`C16`).""")
nb.recap("R01", "Static stability is decided by the potential density", r"""
A parcel moved up expands and its density falls even if nothing else changes. What decides whether it sinks back is the density
it would have at a common pressure — the potential density $\rho_\theta$:
$\dfrac{d\rho_\theta}{dz}=\dfrac{d\rho}{dz}+\dfrac{g\rho}{c^2}$ (13.1), where $c$ is the speed of sound ($c_s$ in our own
lines). Stable when $d\rho_\theta/dz<0$.""",
         where=r"Ch. 1 §1.10, where the same statement is $\frac{d\rho_\theta}{dz}=\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\cong\frac{d\rho}{dz}+\frac{\rho g}{c^2}$ (1.35)")
code(r"""
rho_w, c_sound, drho_dz = 1027.0, 1500.0, -5.24e-3             # sea-water density [kg/m^3], sound speed [m/s], in-situ gradient [kg/m^4] (our inputs)
adiab = STRAT.isentropic_density_gradient(rho_w, c_sound)      # adiabatic gradient -g rho / c^2: what compression alone does
pot = STRAT.ocean_potential_density_gradient(drho_dz, rho_w, c_sound)   # Eq. (13.1): potential-density gradient
print(f"adiabatic gradient  = {adiab:.3e} kg m^-4")            # the part due to compression
print(f"potential gradient  = {pot:.3e} kg m^-4  (negative -> stable)")   # the part that decides stability
""", explain=r"""
1. `isentropic_density_gradient` returns $-g\rho/c_s^2$: how fast density rises with depth just because deeper water is squeezed.
2. `ocean_potential_density_gradient` adds $g\rho/c_s^2$ to the in-situ gradient, as $\frac{d\rho_\theta}{dz}=\frac{d\rho}{dz}+\frac{g\rho}{c^2}$ *(13.1)* says; a negative result means stable.""")
note("N03 [B]", r"""
**Most of the increase of in-situ density with depth is compression.** The budget for our inputs is computed below: only a
small part of the measured gradient is real stratification. So $N^2$ must be built from the *potential* density — the in-situ
gradient would overstate it several-fold.""")
code(r"""
b = ch13.ocean_density_gradient_budget(drho_dz=-5.24e-3, rho=1027.0, c=1500.0)   # split the in-situ gradient into its two parts
print(f"in situ {b['in_situ']:.3e} = adiabatic {b['adiabatic']:.3e} + potential {b['potential']:.3e}  [kg m^-4]")   # the budget
print(f"compression share = {100*b['compression_share']:.1f} %;  in-situ / potential = {b['in_situ']/b['potential']:.1f}")   # how misleading the raw gradient is
""", explain=r"""
**What does this show?** About 85 % of the in-situ gradient is compression; using it for $N^2$ would overstate the
stratification almost sevenfold (the printed ratio).""")
nb.recap("R02", "The adiabatic gradient, in two sign conventions", r"""
The adiabatic density gradient is $-g\rho/c^2$. For a gas the same idea is the adiabatic temperature gradient.
**Kundu:** $\Gamma\equiv dT/dz$, adiabatic value $\Gamma_a=-g/C_p\approx-9.8$ K/km, stable when $dT/dz>\Gamma_a$.
**Meteorology:** $\Gamma_{met}\equiv-dT/dz$, dry adiabatic value $\Gamma_d=+g/C_p\approx+9.8$ K/km, stable when
$\Gamma_{met}<\Gamma_d$. Same physics: negate the number, flip the inequality. This section of the book writes no symbol
$\Gamma$; it quotes rates of decrease in words, which are the meteorological magnitudes.""",
         where="Ch. 1 §1.10 (C52–C55) and Ch. 12 (C14)")
code(r"""
tab = ch13.lapse_rate_table(-6.5e-3, Gamma_a=Gamma_a)          # the standard troposphere, dT/dz = -6.5 K/km, in both conventions
display(tab)                                                   # two rows: Kundu and meteorology, same verdict
""", explain=r"""
**What does this show?** One atmosphere, two rows. Kundu's row compares $dT/dz=-6.5$ with $\Gamma_a=-9.8$ K/km (stable
because $-6.5>-9.8$); the meteorological row compares $\Gamma_{met}=6.5$ with $\Gamma_d=9.8$ K/km (stable because $6.5<9.8$).
The number is negated and the inequality flipped; the verdict is the same.""")
trap("T2, the library string", r"""
The Chapter 1 function `STRAT.lapse_rate_stability(..., convention='meteorology').text` prints "Γ < Γa". Read it as
$\Gamma_{met}<\Gamma_d$: in that convention the library's Γa is the positive number +9.8 K/km. We do not show that raw string
again; every label below is built by this notebook and carries both forms.""")
nb.worked_example("is a −8 K/km layer stable?", r"""
1. **Kundu:** $dT/dz=-8$ K/km; is $-8>-9.8$? Yes → stable.
2. **Meteorology:** negate the number, $\Gamma_{met}=+8$ K/km; flip the inequality: is $8<9.8$? Yes → stable.
3. **Margin:** 1.8 K/km in both.
4. **A −11 K/km layer:** $-11>-9.8$ is false; $11<9.8$ is false → unstable in both.""")
note("N04 [B]", r"""
**The standard atmosphere** (our version of the book's temperature sketch, from the 1976 standard atmosphere): the
*troposphere*, where temperature falls with height, ends at the *tropopause*; above it the *stratosphere* is first isothermal
and then warms up to the *stratopause*. The right panel is the buoyancy frequency squared, $N^2$.""")
fig(r"""
z_atm = np.linspace(0.0, 5.0e4, 501)                           # height from the ground to 50 km [m]
atm = ch13.atmosphere_layers(z_atm)                            # temperature T [K], dT/dz [K/m] (Kundu sign), N^2 [1/s^2], layer names
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
a.plot(atm["T"], z_atm/1e3, color=COLORS["blue"], label="T(z); troposphere: dT/dz = −6.5 > −9.8 K/km ⇔ Γ_met = 6.5 < 9.8 K/km")   # temperature profile, both conventions in the label
for name, zc in (("troposphere", 5.5), ("tropopause", 11.0), ("stratosphere", 30.0), ("stratopause", 47.0)):   # name the layers
    a.text(atm["T"].max() - 2, zc, name, ha="right", fontsize=9, color=COLORS["muted"])   # label at its height
a.set(xlabel="temperature T [K]", ylabel="height z [km]", title="Temperature falls, then rises")
a.legend(fontsize=7, loc="upper center")
b.plot(atm["N2"]*1e4, z_atm/1e3, color=COLORS["accent"])       # N^2 in units of 1e-4 1/s^2
b.set(xlabel="$N^2$ [10$^{-4}$ s$^{-2}$]", title="$N^2$ jumps at the tropopause")
plt.show()
N2_ground = atm["N2"][0]                                       # N^2 at the ground [1/s^2]
i20 = np.argmin(abs(z_atm - 2.0e4))                            # index of the level z = 20 km
print(f"ground: N^2 = {N2_ground:.3e} 1/s^2, N = {np.sqrt(N2_ground):.3e} 1/s, period 2 pi/N = {2*np.pi/np.sqrt(N2_ground)/60:.1f} min")   # tropospheric values
print(f"20 km (lower stratosphere): N^2 = {atm['N2'][i20]:.3e} 1/s^2 = {atm['N2'][i20]/N2_ground:.1f} x the ground value")   # the jump
""",
    see="Left: temperature against height, falling at 6.5 K/km through the troposphere, constant in the lower stratosphere, "
        "rising above. Right: $N^2$, a few times larger in the stratosphere than in the troposphere.",
    read="The troposphere is stable in both conventions — $dT/dz=-6.5>-9.8$ K/km, i.e. $\\Gamma_{met}=6.5<9.8$ K/km — and "
         "$N^2$ jumps at the tropopause: the stratosphere acts as a lid on the weather layer. The printed lines give the numbers.",
    change="…the troposphere followed the dry adiabat exactly? Then $N^2=0$ there, parcels would feel no restoring force, and "
           "the storm growth rate of `C16` (which is proportional to $1/N$) would be unbounded.")
note("N05 [B]", r"""
⚠️ **Common confusion (trap T2).** "The troposphere is close to neutral" compares the observed lapse rate with the **moist**
adiabat (named only: saturated air cools more slowly as it rises, because condensation releases heat). To *dry* displacements
the standard troposphere is statically stable in both conventions: $\Gamma_a<dT/dz<0$ (Kundu), i.e.
$0<\Gamma_{met}<\Gamma_d$ (meteorology); $N^2>0$. The gradient Richardson number of Chapter 12, $N^2$ divided by the squared
shear, says the same for a typical shear:""")
code(r"""
ri = ch12.gradient_richardson_thermal(-6.5e-3, 3.0e-3, 1/288.15, Gamma_a=Gamma_a)   # dT/dz = -6.5 K/km, shear 3 m/s per km, alpha = 1/T
badge = f"stable · dT/dz = −6.5 > −{-Gamma_a*1e3:.1f} K/km · Γ_met = 6.5 < {-Gamma_a*1e3:.1f} K/km"   # OUR label, both conventions
print(f"Ri = N^2/(dU/dz)^2 = {ri['Ri']:.1f}   ->  {badge}")    # far above the critical value 1/4 of Ch. 11
""", explain=r"""
**What does this show?** With a shear of 3 m/s per kilometre the Richardson number is about 12 — some fifty times the
critical value ¼ of Chapter 11. The label is built here so that it carries both sign conventions.""")
nb.md(r"""
#### 🎚️ The same verdict in two conventions

Drag the environmental temperature gradient. A parcel is displaced by $\zeta$ metres and the figure shows the buoyancy
acceleration it then feels, $-N^2\zeta$, with $N^2=\frac{g}{T}\big(\frac{dT}{dz}-\Gamma_a\big)$.""")
nb.plotly(r"""
zeta_ax = np.linspace(-200.0, 200.0, 9)                        # vertical displacement of the parcel [m]
nan_ax = np.full_like(zeta_ax, np.nan)                         # an empty curve (hides a trace)

def f2(dTdz_km):                                               # curves for one environmental gradient dT/dz [K/km] (Kundu sign)
    N2 = STRAT.brunt_vaisala_sq_from_lapse(288.15, dTdz_km*1e-3)   # N^2 = (g/T)(dT/dz - Gamma_a) [1/s^2]
    acc = -N2*zeta_ax                                          # buoyancy acceleration of the displaced parcel [m/s^2]
    return {"stable: dT/dz > Γa = −9.8 K/km ⇔ Γ_met < Γd = +9.8 K/km (pushed back)": (zeta_ax, acc if N2 > 0 else nan_ax),
            "unstable: dT/dz < Γa = −9.8 K/km ⇔ Γ_met > Γd = +9.8 K/km (pushed further)": (zeta_ax, acc if N2 < 0 else nan_ax),
            "neutral reference: dT/dz = Γa ⇔ Γ_met = Γd (no force)": (zeta_ax, 0*zeta_ax)}

figF2 = slider_figure(f2, "dT/dz (Kundu; Γ_met is minus this)", np.arange(-12.0, 2.01, 0.5), unit="K/km",
                      xlabel="displacement ζ of the parcel [m]", ylabel="buoyancy acceleration −N²ζ [m/s²]",
                      title="The verdict flips at the adiabatic value, whichever sign you write it with")   # 29 slider positions
figF2.show()
for g_km in (-6.5, -9.76, -11.0):                              # three sample gradients [K/km]
    t2 = ch13.lapse_rate_table(g_km*1e-3, Gamma_a=Gamma_a)     # both conventions for this gradient
    print(f"dT/dz = {g_km:+.2f} K/km:  {t2.loc['Kundu', 'criterion']}  ⇔  {t2.loc['meteorology', 'criterion']}  →  {t2.loc['Kundu', 'verdict']}")   # one line each
""", explain=r"""
1. `f2` computes $N^2$ for one gradient and returns the acceleration $-N^2\zeta$ on the *stable* trace or on the *unstable*
   one (the other is left empty), so the legend always names the regime in **both** conventions.
2. The three printed lines are `ch13.lapse_rate_table` at a stable, a neutral and an unstable gradient.""")
nb.figure_notes(
    see="A straight line through the origin whose slope is $-N^2$. For gradients to the right of −9.8 K/km on the slider the "
        "slope is negative (a displaced parcel is pushed back); to the left it is positive (pushed further).",
    read="The verdict flips where the line passes through the horizontal — at the adiabatic value, $dT/dz=\\Gamma_a=-9.8$ K/km "
         "in Kundu's sign, $\\Gamma_{met}=\\Gamma_d=+9.8$ K/km in the meteorological one. The slider shows Kundu's number; the "
         "meteorological number is its negative.",
    change="…the air were an inversion, $dT/dz=+2$ K/km ($\\Gamma_{met}=-2$ K/km)? The line is at its steepest: the strongest "
           "restoring force on the slider, and the largest $N^2$.")
note("N06 [B]", r"""
**The ocean** (our idealised version of the book's sketch; `ch13.ocean_profile_idealized`): a warm *mixed layer* stirred by the
wind, the *thermocline* where temperature falls quickly, and the cold, nearly uniform *abyss*. Salinity is secondary and left
out of our profile. The same shape of $N(z)$ feeds `C08`.""")
fig(r"""
z_oc = np.linspace(-4200.0, 0.0, 841)                          # depth from the bottom (-4200 m) to the surface [m]
prof = ch13.ocean_profile_idealized(z_oc)                      # T [deg C], potential density [kg/m^3], N^2 and N (analytic)
fig, axes = plt.subplots(1, 3, figsize=(9.6, 4.0), sharey=True)
axes[0].plot(prof["T"], z_oc, color=COLORS["blue"])            # temperature profile
axes[0].set(xlabel="temperature T [°C]", ylabel="height z [m] (surface at 0)", title="Warm on top")
axes[1].plot(prof["rho_theta"], z_oc, color=COLORS["blue"])    # potential density profile
axes[1].set(xlabel=r"potential density $\rho_\theta$ [kg m$^{-3}$]", title="Light over heavy")
axes[2].plot(prof["N"]*1e3, z_oc, color=COLORS["accent"])      # buoyancy frequency in units of 1e-3 1/s
axes[2].set(xlabel="$N$ [10$^{-3}$ s$^{-1}$]", title="$N$ peaks in the thermocline")
for ax in axes:                                                # name the three layers in every panel
    ax.axhspan(-50, 0, color=COLORS["grid"], alpha=0.7)        # the mixed layer (top 50 m)
    ax.axhspan(-450, -50, color=COLORS["amber"], alpha=0.12)   # the thermocline
axes[0].text(6, -800, "thermocline\n(amber band, above)", fontsize=8)   # label, placed below the band so that it can be read
axes[0].text(4, -2500, "abyss", fontsize=9)                    # label
plt.show()
""",
    see="Temperature, potential density and buoyancy frequency against depth. The change is fastest just under the thin mixed "
        "layer and fades with depth; the amber band marks the upper 450 m, and below about 1500 m the ocean is cold and nearly "
        "uniform.",
    read="The buoyancy frequency follows the slope of the density curve: zero in the mixed layer (the thin band at the very "
         "top), largest at the top of the thermocline just beneath it, and decaying smoothly into the abyss.",
    change="…the thermocline were deeper? The peak of $N$ moves down, and the first internal mode of `C08` travels faster.")
nb.recap("R03", "The buoyancy frequency", r"""
Its definition is $N^2\equiv-\dfrac{g}{\rho_0}\dfrac{d\rho}{dz}$ with $\rho$ the **potential** density and $\rho_0$ a constant reference value:
the squared frequency at which a displaced parcel bobs. Where it enters this chapter: `C08` (mode speeds
$c_n\approx NH/n\pi$), `C14` (internal waves need $\omega<N$), `C16` (the Eady radius $NH/f$ and the growth rate).""",
         where=r"Ch. 1, $N^2=-\frac g\rho\big(\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\big)$ (1.29); Ch. 7, $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ (7.127)")
code(r"""
N2_num = GFD.buoyancy_frequency_sq(z_oc, prof["rho_theta"], 1027.0)   # N^2 = -(g/rho0) d(rho_theta)/dz by finite differences
inner = (z_oc < -60.0) & (z_oc > -4190.0)                      # below the base of the mixed layer (a kink at -50 m) and off the bottom
assert np.allclose(N2_num[inner], prof["N2"][inner], rtol=2e-2, atol=1e-12)   # the numerical derivative matches the analytic profile
print(f"✓ N^2 from differences and the analytic N^2 agree to 2 % (largest N^2 = {prof['N2'].max():.2e} 1/s^2)")   # visible confirmation
""", explain=r"""
**What does this show?** Differentiating the potential-density profile numerically reproduces the profile's own analytic
$N^2$: the function does exactly what the definition says.""")
note("N07 [C]", r"""
**The fastest bobbing in our thermocline.** The largest $N$ of the profile gives the shortest buoyancy period, computed below.
It is the upper frequency limit of the internal waves of `C14` and the profile maximum of `C08`.""")
code(r"""
N_max = np.sqrt(prof["N2"].max())                              # largest buoyancy frequency of the profile [1/s]
print(f"N_max = {N_max:.2e} 1/s  ->  buoyancy period 2 pi/N = {STRAT.stability_timescale(prof['N2'].max())[1]/60:.0f} min")   # the shortest period
""", explain=r"""
**What does this show?** A parcel in the sharpest part of our thermocline bobs with a period of some tens of minutes; nothing
internal can oscillate faster.""")

# =====================================================================================================================
# A.3  §13.3 Equations of motion
# =====================================================================================================================
nb.section("13.3", "Equations of Motion", intro=r"""
**What is this section about?** The starting equations — Chapter 4's Boussinesq set written in a rotating frame — and the one
new modelling choice: friction by eddies, with a large horizontal and a small vertical eddy viscosity.""")
nb.recap("R04", "The rotating Boussinesq equations", r"""
Continuity $\nabla\cdot\mathbf u=0$; momentum
$\dfrac{D\mathbf u}{Dt}+2\boldsymbol\Omega\times\mathbf u=-\dfrac1{\rho_0}\nabla p-\dfrac{g\rho}{\rho_0}\mathbf e_z+\mathbf F$;
density $\dfrac{D\rho}{Dt}=0$ (13.2). Term by term (momentum equation): acceleration seen from the turning earth, the Coriolis
acceleration, the pressure-gradient force, weight, friction per unit mass.""",
         where=r"Ch. 4: the rotating-frame equation $\rho\big(\frac{D'\mathbf u'}{Dt}\big)_{O'1'2'3'}=-\nabla'p+\rho\big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\big]+\mu\nabla'^2\mathbf u'$ (4.45) and the Boussinesq momentum equation $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ (4.86)")
slip(1, r"the pressure term of the momentum equation as $+\frac1{\rho_0}\nabla p$",
     r"$-\frac1{\rho_0}\nabla p$ (a fluid at rest must satisfy $0=-\nabla p-g\rho\,\mathbf e_z$).")
code(r"""
ok_corr = ch13.boussinesq_rotating_sympy()                     # the corrected momentum equation, tested on the state of rest
ok_prnt = ch13.boussinesq_rotating_sympy(printed=True)         # the same test with the sign as printed
print("corrected sign: rest state satisfies it ->", ok_corr["ok"])                 # True
print("printed sign:   rest state satisfies it ->", ok_prnt["ok"], "| left over:", ok_prnt["residual"])   # False, with a residual
""", explain=r"""
**What does this show?** A fluid at rest in hydrostatic balance satisfies the corrected equation exactly. With the printed
sign a term twice the weight is left over — the printed form cannot even describe a motionless ocean.""")
nb.recap("R05", "When Boussinesq holds", r"""
The layer must be thin against the scale height $c^2/g$ (written $c_s^2/g$ in our lines), the height over which density
changes by its own size through compression.""", where="Ch. 4 §4.9")
code(r"""
print(f"scale height, sea water (c_s = 1500 m/s): {GFD.scale_height(1500.0)/1e3:.0f} km")   # far deeper than any ocean
print(f"scale height, air       (c_s =  340 m/s): {GFD.scale_height(340.0)/1e3:.1f} km")    # comparable with the troposphere
""", explain=r"""
**What does this show?** The condition holds for any ocean depth, but only roughly for the troposphere — which is why
meteorology often uses pressure as the vertical coordinate (not in the book).""")
nb.recap("R06", "Density is carried with the fluid", r"""
The statement $D\rho/Dt=0$ follows from $DT/Dt=0$ (or $DS/Dt=0$) and a linear equation of state, $\delta\rho/\rho_0=-\alpha\,\delta T$ and
$\delta\rho/\rho_0=\beta_S\,\delta S$. ⚠️ The book writes the haline coefficient as $\beta$; from §13.4 on $\beta$ is $df/dy$,
so we write $\beta_S$ here. The expansion coefficient $\alpha$ comes back in `C03`.""",
         where=r"Ch. 4, the Boussinesq heat equation $DT/Dt=\kappa\nabla^2T$ (4.89) with its diffusion term dropped; Ch. 1 linear equation of state")
nb.recap("R07", "The state of rest", r"""
It obeys $\dfrac{d\bar p}{dz}=-\bar\rho g$ (13.3): the hydrostatic balance of the motionless reference state.""",
         where=r"Ch. 1, $dp/dz=-\rho g$ (1.8)")
nb.recap("R08", "Rest state plus perturbation", r"""
Write $\rho(\mathbf x,t)=\bar\rho(z)+\rho'(\mathbf x,t)$ and $p(\mathbf x,t)=\bar p(z)+p'(\mathbf x,t)$ (13.4).
⚠️ **Trap T3:** from the thin-shell equations of `C01` on, the book drops the primes: $p$ and $\rho$ there are *perturbations*.""",
         where=r"Ch. 4, where $p'=p-p_s$, $\rho'=\rho-\rho_s$ give $\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u$ (4.84); Ch. 7, $\rho=\bar\rho(z)+\rho'$ (7.124)")
nb.recap("R09", "The pressure and gravity terms keep their form", r"""
Subtract the rest state: $-\dfrac1{\rho_0}\nabla p-\dfrac{g\rho}{\rho_0}\mathbf e_z=-\dfrac1{\rho_0}\nabla p'-\dfrac{g\rho'}{\rho_0}\mathbf e_z$,
because $\nabla\bar p=(d\bar p/dz)\,\mathbf e_z=-\bar\rho g\,\mathbf e_z$ cancels the weight of the rest state.""",
         where="Ch. 4, derivation D27")
code(r"""
print("perturbation form of pressure + gravity terms holds ->", ch13.perturbation_form_sympy()["ok"])   # the engine agrees: True
""", explain=r"""
**What does this show?** A sympy engine subtracts the rest state symbolically and finds nothing left over.""")
nb.recap("R10", "Friction per unit mass, and the eddy-viscosity idea", r"""
The friction force per unit mass is the divergence of the stress divided by density,
$F_i=\dfrac1\rho\dfrac{\partial\tau_{ij}}{\partial x_j}$. For large-scale flow the stress is carried by turbulent eddies and is
modelled as an eddy viscosity times a velocity gradient.""",
         where=r"Ch. 4 §4.4 and Ch. 12, the eddy-viscosity hypothesis $\overline{u_iu_j}=\frac23\bar e\delta_{ij}-\nu_T\big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\big)$ (12.94)")
slip(2, r"$F_i=\partial\tau_{ij}/\partial x_j$ for a force per unit mass",
     r"$F_i=\frac1\rho\,\partial\tau_{ij}/\partial x_j$ (N/m³ divided by kg/m³ gives m/s²).")
note("N08 [B]", r"""
**The six eddy stresses:**
$\tau_{xz}=\tau_{zx}=\rho\nu_v\frac{\partial u}{\partial z}+\rho\nu_H\frac{\partial w}{\partial x}$,
$\tau_{yz}=\tau_{zy}=\rho\nu_v\frac{\partial v}{\partial z}+\rho\nu_H\frac{\partial w}{\partial y}$,
$\tau_{xy}=\tau_{yx}=\rho\nu_H\big(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\big)$,
$\tau_{xx}=2\rho\nu_H\frac{\partial u}{\partial x}$, $\tau_{yy}=2\rho\nu_H\frac{\partial v}{\partial y}$,
$\tau_{zz}=2\rho\nu_v\frac{\partial w}{\partial z}$ (13.5). Two coefficients because the eddies are wide and flat: horizontal exchange
is strong, vertical exchange weak, $\nu_H\gg\nu_v$.""")
code(r"""
rng = np.random.default_rng(13)                                # a fixed random seed, so the numbers repeat
Gm = rng.normal(size=(3, 3))                                   # a velocity-gradient matrix G[i, j] = d u_i / d x_j [1/s]
Gm -= np.trace(Gm)/3*np.eye(3)                                 # make it divergence-free (incompressible)
tau_iso = ch13.anisotropic_eddy_stress(Gm, 5.0, 5.0, 1027.0)   # Eq. (13.5) with nu_H = nu_v = 5 m^2/s
S = 0.5*(Gm + Gm.T)                                            # strain-rate tensor S_ij of Ch. 4
assert np.allclose(tau_iso, 2*1027.0*5.0*S)                    # equal coefficients -> the Newtonian form 2 rho nu S_ij
print("✓ with nu_H = nu_v the eddy stress equals 2 rho nu S_ij (max difference", f"{abs(tau_iso - 2*1027.0*5.0*S).max():.1e} Pa)")   # visible check
""", explain=r"""
**What does this show?** With the two coefficients equal, the six stresses collapse to the Newtonian stress
$2\rho\nu S_{ij}$ of Chapter 4. The anisotropy is the only new ingredient.""")
note("N09 [C]", r"""
**A defect the book accepts.** This stress is not frame-indifferent: a rigid rotation in a vertical plane, $u=\omega_rz$,
$w=-\omega_rx$, has no deformation at all, yet it gives $\tau_{xz}=\rho\,\omega_r(\nu_v-\nu_H)\neq0$ (the cell below). The
book accepts the defect for simplicity; so do we. Pointer: Ch. 4 §4.5.""")
code(r"""
G_rot = np.zeros((3, 3)); G_rot[0, 2], G_rot[2, 0] = 1.0, -1.0 # rigid rotation at 1 rad/s: du/dz = +1, dw/dx = -1, no strain
tau_rot = ch13.anisotropic_eddy_stress(G_rot, 300.0, 0.03, 1027.0)   # nu_H = 300, nu_v = 0.03 m^2/s (our illustrative values)
print(f"tau_xz for a rigid rotation = {tau_rot[0, 2]:.3e} Pa  (rho*(nu_v - nu_H) = {1027.0*(0.03 - 300.0):.3e})")   # non-zero: the defect
""", explain=r"""
**What does this show?** A motion without any deformation produces a stress — impossible for a true material law, and
harmless here because large-scale flows never rotate rigidly in a vertical plane.""")
note("N10 [B]", r"""
**The friction force** per unit mass, component by component:
$F_x=\nu_H\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\big)+\nu_v\frac{\partial^2u}{\partial z^2}$,
$F_y=\nu_H\big(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\big)+\nu_v\frac{\partial^2v}{\partial z^2}$,
$F_z=\nu_H\big(\frac{\partial^2w}{\partial x^2}+\frac{\partial^2w}{\partial y^2}\big)+\nu_v\frac{\partial^2w}{\partial z^2}$ (13.6).
The book states this without the steps. The six differentiations and the use of continuity are written out as derivation
`D01` inside the `C01` block below.""")
code(r"""
Fx = GFD.eddy_friction(lap_h=-1.0e-11, d2z=-2.0e-6, nu_H=3.0e4, nu_v=7.0)   # one component: nu_H * (horizontal Laplacian) + nu_v * (d2/dz2)
print(f"F_x = {Fx:.2e} m/s^2  (horizontal part {3.0e4*-1.0e-11:.1e}, vertical part {7.0*-2.0e-6:.1e})")   # the vertical part dominates
""", explain=r"""
**What does this show?** For a wind that varies over 1000 km horizontally and 1 km vertically, the vertical mixing term is the
larger one even though $\nu_v\ll\nu_H$ — because the vertical scale is so much shorter.""")
note("N11 [C]", r"""
**How big are the eddy coefficients?** A small table of **our own illustrative** values (not the book's) against molecular
ones. They are used where $\nu_v$ sets the Ekman thickness (`C04`, `C06`).

| Fluid | $\nu_v$ (eddy, illustrative) [m² s⁻¹] | $\nu_H$ (eddy, illustrative) [m² s⁻¹] | molecular $\nu$ [m² s⁻¹] |
|---|---|---|---|
| lower atmosphere | 7 | 3 × 10⁴ | 1.5 × 10⁻⁵ |
| upper ocean | 0.03 | 300 | 1.0 × 10⁻⁶ |""")


# =====================================================================================================================
# A.4  §13.4 — C01 the thin-shell equations
# =====================================================================================================================
nb.section("13.4", "Approximate Equations for a Thin Layer on a Rotating Sphere", intro=r"""
**What is this section about?** Standing at one latitude, we replace the sphere by its tangent plane and ask which parts of the
earth's rotation a thin layer of fluid can feel. The answer is one number per latitude, $f$, and a set of equations that every
later section simplifies further.""")
core("C01", r"The equations of motion on a thin rotating shell: $\dfrac{Du}{Dt}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}+F_x$, "
            r"$\dfrac{Dv}{Dt}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}+F_y$, "
            r"$\dfrac{Dw}{Dt}=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial z}-\dfrac{g\rho}{\rho_0}+F_z$ (13.9) — $(F_x,F_y,F_z)$ being the eddy friction of note N10 — with $f=2\Omega\sin\theta$ (13.8)",
     "Which part of the earth's rotation does a thin layer of air or water actually feel?")
problem(r"""
The atmosphere is about 10 km deep and its weather systems are 1000 km and more across: a sheet thinner, in proportion, than the
paper of a wall map. In such a sheet the fluid can hardly move up or down. Of the earth's spin, the part that matters is the
part about the *local vertical* — a turntable under your feet that turns once a day at the pole and not at all on the equator.
Climate models call dropping the rest the "traditional approximation"; every dynamical core makes it.""")
idea(r"""
            Ω (earth's axis)
            ↑        at latitude θ the axis is tilted θ above the northward horizontal
     Ω sin θ ↑  ↗ Ω
   (vertical)| /
             |/____→  Ω cos θ  (northward, horizontal)

   Ω sin θ : spins the ground under you            → turns horizontal motion
   Ω cos θ : tips the ground about a north–south line → only couples to vertical motion
""", r"""
| Place | Vertical part $\Omega\sin\theta$ | Horizontal part $\Omega\cos\theta$ | Coriolis parameter |
|---|---|---|---|
| pole ($\theta=90°$) | all of it | none | $f=2\Omega$ |
| equator ($\theta=0$) | none | all of it | $f=0$ |""")
nb.recap("R11", "How fast the earth turns", r"""
One turn against the stars takes 86 164 s (a sidereal day), so $\Omega=2\pi/86\,164\ \mathrm{s}=7.292\times10^{-5}$ rad/s
(`ROT.OMEGA_EARTH`). ⚠️ **Trap T1:** the book takes $\Omega=2\pi$ rad/day with a 24-hour solar day, 0.27 % smaller. Every
number of ours uses the sidereal value.""", where="Ch. 4 §4.7")
nb.current_core = "C01"
P("P308", "sidereal day vs solar day", r"""
In one year the earth turns 366.25 times against the stars but the sun crosses the sky only 365.25 times, because one turn is
used up going round the sun. So one true rotation takes 24 h × 365.25/366.25 = 23 h 56 min 4 s. Dynamics cares about the true
rotation.""", code=r"""
T_sid = 86400.0 * 365.25 / 366.25                # one rotation against the stars [s]
print(T_sid, 2*np.pi/T_sid)                      # 86164.1 s and 7.2921e-05 rad/s
print(GFD.OMEGA_EARTH / (2*np.pi/86400.0) - 1)   # 0.00274: the solar-day value is 0.27 % smaller
""")
note("N13 [B]", r"""
**The rotation vector on the tangent plane.** With $x$ east, $y$ north and $z$ up at latitude $\theta$,
$\boldsymbol\Omega=(0,\ \Omega\cos\theta,\ \Omega\sin\theta)$: no eastward part, a northward horizontal part and a vertical
part. The figure lets you move the tangent plane to four latitudes.""")
nb.plotly(r"""
uu, vv = np.meshgrid(np.linspace(0, 2*np.pi, 40), np.linspace(-np.pi/2, np.pi/2, 20))   # longitude and latitude grids for a unit sphere
figS = go.Figure(go.Surface(x=np.cos(vv)*np.cos(uu), y=np.cos(vv)*np.sin(uu), z=np.sin(vv), opacity=0.25, showscale=False,
                            colorscale=[[0, COLORS["grid"]], [1, COLORS["grid"]]], hoverinfo="skip"))   # the earth, pale grey
figS.add_trace(go.Scatter3d(x=[0, 0], y=[0, 0], z=[-1.4, 1.4], mode="lines", line=dict(color=COLORS["ink"], width=5), name="rotation axis"))   # the axis
lats_deg = [0.0, 35.0, 60.0, 90.0]                              # the four dropdown latitudes
for j, ld in enumerate(lats_deg):                               # three arrows per latitude
    th = np.deg2rad(ld)                                         # latitude in radians
    comp = GFD.earth_rotation_local(th)/GFD.OMEGA_EARTH         # (0, cos(theta), sin(theta)): local components of Omega, in units of Omega
    P0 = np.array([np.cos(th), 0.0, np.sin(th)])                # the point on the sphere (longitude 0)
    up = P0                                                     # local vertical unit vector
    north = np.array([-np.sin(th), 0.0, np.cos(th)])            # local northward unit vector
    for vec, col, nm in ((comp[2]*up, COLORS["teal"], f"vertical part Ω sin θ = {comp[2]:.2f} Ω (= f/2)"),
                         (comp[1]*north, COLORS["orange"], f"northward part Ω cos θ = {comp[1]:.2f} Ω"),
                         (comp[1]*north + comp[2]*up, COLORS["accent"], "Ω itself (parallel to the axis)")):
        Q = P0 + 0.6*vec                                        # tip of the arrow (scaled to fit the picture)
        figS.add_trace(go.Scatter3d(x=[P0[0], Q[0]], y=[P0[1], Q[1]], z=[P0[2], Q[2]], mode="lines+markers", name=nm,
                                    marker=dict(size=[0, 5], color=col), line=dict(color=col, width=7), visible=(j == 1)))   # one arrow
buttons = [dict(label=f"latitude {ld:.0f}°", method="update",
                args=[{"visible": [True, True] + [k//3 == j for k in range(3*len(lats_deg))]}]) for j, ld in enumerate(lats_deg)]   # show one latitude's arrows
figS.update_layout(height=520, title="The local vertical part of the earth's spin shrinks toward the equator",
                   updatemenus=[dict(buttons=buttons, active=1, x=0.0, y=1.08, xanchor="left")],
                   scene=dict(aspectmode="data", xaxis_visible=False, yaxis_visible=False, zaxis_visible=False), margin=dict(l=0, r=0, t=60, b=0))   # layout and dropdown
figS.show()
""", explain=r"""
1. The pale sphere and the black axis are drawn once.
2. For each of four latitudes, `GFD.earth_rotation_local` gives the components $(0,\Omega\cos\theta,\Omega\sin\theta)$; three
   arrows are drawn at the point: the vertical part (teal), the northward part (orange) and their sum (purple), which is
   parallel to the axis.
3. The dropdown only switches which three arrows are visible.""")
nb.figure_notes(
    see="A point on the globe with three arrows: the earth's rotation vector (purple, always parallel to the axis) and its two "
        "local parts — along the local vertical (teal) and toward the north (orange).",
    read="The teal arrow is $f/2$. At 35° it is 57 % of $\\Omega$; at the pole it is all of $\\Omega$; on the equator it vanishes "
         "and the whole rotation is 'horizontal'.",
    change="…you stood at 35° S? The teal arrow would point *into* the ground: $f<0$, and everything that turns right in the "
           "north turns left.")
note("N12 [B]", r"""
**How small is the vertical velocity?** Continuity $\partial u/\partial x+\partial v/\partial y+\partial w/\partial z=0$ with
horizontal scale $L$, depth $H$ and velocity scales $U$, $W$ gives $U/L\sim W/H$, so $W/U\sim H/L$: the aspect ratio. It is
used in step 3 of `D02`.""")
code(r"""
W = GFD.vertical_velocity_scale(14.0, 9000.0, 1.4e6)           # W = U H / L for U = 14 m/s, H = 9 km, L = 1400 km [m/s]
print(f"W = {W:.3f} m/s,  aspect ratio H/L = {9000.0/1.4e6:.1e}")   # centimetres per second against metres per second
""", explain=r"""
**What does this show?** A 14 m/s wind in a weather system comes with vertical motion of only about 9 cm/s: the flow is
almost horizontal.""")
D("D01", ref="13.6")
code(r"""
print("friction force (13.6) follows from the stresses (13.5) ->", ch13.eddy_friction_force_sympy()["ok"])   # sympy repeats D01: True
""", explain=r"""
**What does this show?** The sympy engine differentiates the six stresses, uses continuity and finds the three friction
components of `D01` with nothing left over.""")
D("D02", ref="13.9")
note("N14 [B]", r"""
**The approximation has a name.** Before any approximation the $x$-component of the Coriolis acceleration is
$2\Omega(w\cos\theta-v\sin\theta)$. Dropping $w\cos\theta$ against $v\sin\theta$ (steps 3–4 of `D02`) and the vertical
component against gravity (step 6) leaves $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7). This
is the "traditional approximation" (the book does not name it).""")
trap("T18", r"""
The text calls $-2\Omega u\cos\theta$ the vertical component of the Coriolis *force*. It is the third component of the Coriolis
*acceleration* $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7), a term on the left-hand side of the
vertical equation; the force per unit mass is its negative, $+2\Omega u\cos\theta$ (upward for an eastward wind).
`GFD.coriolis_acceleration_local` returns the acceleration.""")
code(r"""
print("thin-shell equations (13.9) follow from the rotating Boussinesq set ->", ch13.thin_layer_sympy()["ok"])   # the engine agrees: True
""", explain=r"""
**What does this show?** The same reduction done symbolically: resolving $\boldsymbol\Omega$, taking the cross product and
dropping the two small terms gives the three thin-shell equations.""")
nb.recap("R12", "The Coriolis parameter", r"""
It is $f=2\Omega\sin\theta$ (13.8): twice the local vertical component of the earth's rotation, also called the planetary vorticity.
Positive in the northern hemisphere, negative in the southern, zero on the equator, odd in latitude.""",
         where="Ch. 4 §4.7, Ch. 5 §5.6")
nb.current_core = "C01"
code(r"""
for name, L in (("35 N", lat), ("35 S", -lat), ("pole", np.pi/2)):   # three latitudes [rad]
    print(f"f at {name}: {GFD.coriolis_parameter(L):+.4e} 1/s")       # Eq. (13.8): f = 2 Omega sin(latitude)
""", explain=r"""
**What does this show?** The parameter is odd in latitude (equal and opposite at 35° N and 35° S) and largest at the pole,
where it equals $2\Omega$.""")
note("N15 [B]", r"""
**The inertial period**, $T_i=2\pi/f$, is the time the Coriolis force takes to turn a moving parcel through a full circle
(`C10` shows the circle). It is half a "pendulum day" and grows toward the equator.""")
code(r"""
for name, f_ in (("35 N", f35), ("60 N", f60), ("12 N", f12), ("pole", GFD.coriolis_parameter(np.pi/2))):   # four latitudes
    print(f"inertial period at {name}: {GFD.inertial_period(f_)/3600:.2f} h")   # T_i = 2 pi / |f|, in hours
""", explain=r"""
**What does this show?** About 21 hours at 35°, 14 hours at 60°, 12 hours at the pole (half a sidereal day), and more than two
days at 12°.""")
slip(5, r"after $T_i=2\pi/f$, a remark that the subscript of $T_i$ stands for a vector component",
     r'that it does **not**: in $T_i=2\pi/f$ the letter i is only a label, for "inertial".')
nb.worked_example("f, β and the size of what we dropped at 30° N", r"""
Take $\Omega\approx7.3\times10^{-5}$ s⁻¹, $R\approx6400$ km, $\sin30°=0.5$, $\cos30°\approx0.87$.

1. $f=2\Omega\sin\theta=2\times7.3\times10^{-5}\times0.5=7.3\times10^{-5}$ s⁻¹.
2. Inertial period $2\pi/f\approx86\,000$ s ≈ 24 h (at 30° the inertial period is one day — half the pendulum day).
3. $\beta=2\Omega\cos\theta/R=1.46\times10^{-4}\times0.87/6.4\times10^{6}\approx2.0\times10^{-11}$ m⁻¹ s⁻¹.
4. A storm with $U=20$ m/s, $L=2000$ km, $H=10$ km: $W\approx UH/L=0.1$ m/s.
5. Dropped term over kept term: $(W\cos\theta)/(U\sin\theta)=(0.1\times0.87)/(20\times0.5)\approx0.009$ — under one per cent.""")
code(r"""
terms = GFD.thin_layer_terms(14.0, 1.4e6, 9000.0, lat, nu_H=3.0e4, nu_v=7.0)   # sizes of every term for U = 14 m/s, L = 1400 km, H = 9 km at 35 N
display(ch13.term_table_thin_layer(14.0, 1.4e6, 9000.0, lat, nu_H=3.0e4, nu_v=7.0))   # the same as a table [m/s^2] with ratios to the Coriolis term
full = GFD.coriolis_acceleration_local(14.0, 0.0, 0.09, lat, thin=False)       # exact 2 Omega x u for (u, v, w) = (14, 0, 0.09) m/s
thin = GFD.coriolis_acceleration_local(14.0, 0.0, 0.09, lat, thin=True)        # the thin-layer form (-f v, f u, -2 Omega u cos(theta))
print(f"Ro = {terms['Ro']:.3f};  dropped x-term 2 Omega cos(theta) W = {terms['x']['coriolis_w']:.2e} m/s^2 = {100*terms['x']['coriolis_w']/terms['x']['coriolis']:.1f} % of f U")   # Rossby number, size of the dropped term
print(f"vertical Coriolis 2 Omega U cos(theta) = {terms['z']['coriolis']:.2e} m/s^2 = {terms['z']['coriolis']/G0:.1e} of g")   # tiny against gravity
print("exact [m/s^2]:", [f"{a_:+.2e}" for a_ in full], " thin:", [f"{a_:+.2e}" for a_ in thin])   # the two forms differ only in the x-component
""", explain=r"""
1. `GFD.thin_layer_terms` estimates every term of the three thin-shell equations from the scales $U$, $L$, $H$; the table shows
   them in m/s² and as fractions of the Coriolis term. An entry NaN ("not a number") means only that this term does not
   occur in that equation (there is no buoyancy in the horizontal equation, for instance).
2. The Rossby number $U/(fL)$ is about 0.12: the acceleration is an eighth of the Coriolis term.
3. The dropped piece $2\Omega\cos\theta\,W$ is under 1 % of $fU$ — step 3 of `D02` in numbers.
4. In the vertical, the Coriolis term is some 10⁻⁴ of $g$: rotation is negligible next to weight and pressure.""")
scratch(r"""
# From scratch: f, beta and the full Coriolis acceleration typed out
f_mine = 2*GFD.OMEGA_EARTH*np.sin(lat)                         # Eq. (13.8): f = 2 Omega sin(latitude)
beta_mine = 2*GFD.OMEGA_EARTH*np.cos(lat)/GFD.EARTH_RADIUS_MEAN   # Eq. (13.10): beta = 2 Omega cos(latitude) / R
cor_mine = 2*np.cross(GFD.earth_rotation_local(lat), [14.0, 0.0, 0.09])   # 2 Omega x u with Omega = (0, Omega cos, Omega sin)
assert np.isclose(f_mine, f35) and np.isclose(beta_mine, beta35)           # same numbers as the library
assert np.allclose(cor_mine, GFD.coriolis_acceleration_local(14.0, 0.0, 0.09, lat, thin=False))   # same exact acceleration
gap = np.linalg.norm(cor_mine - np.array(thin))/np.linalg.norm(cor_mine)    # relative size of what the thin-layer form drops
print(f"✓ f, beta and 2 Omega x u agree with the library; the thin-layer form differs by {100*gap:.2f} % of the acceleration")   # visible confirmation
""", r"""
`np.cross(a, b)` is the cross product of two 3-vectors. The definitions typed out reproduce `GFD.coriolis_parameter`,
`GFD.beta_parameter` and `GFD.coriolis_acceleration_local`; the last line measures how little the traditional approximation
changes.""")
fig(r"""
x_terms = [("acceleration", terms["x"]["acceleration"], COLORS["accent"], None), ("Coriolis f U", terms["x"]["coriolis"], COLORS["teal"], None),
           ("pressure", terms["x"]["pressure"], COLORS["orange"], None), ("vertical friction", terms["x"]["friction_v"], COLORS["rose"], None),
           ("horizontal friction", terms["x"]["friction_H"], COLORS["rose"], None), ("dropped: 2Ω cos θ W", terms["x"]["coriolis_w"], COLORS["muted"], "//")]   # the x-equation
z_terms = [("acceleration", terms["z"]["acceleration"], COLORS["accent"], None), ("pressure (balances weight)", terms["z"]["buoyancy"], COLORS["orange"], None),
           ("weight g ρ'/ρ0", terms["z"]["buoyancy"], COLORS["blue"], None), ("vertical friction", terms["z"]["friction_v"], COLORS["rose"], None),
           ("horizontal friction", terms["z"]["friction_H"], COLORS["rose"], None), ("dropped: 2Ω cos θ U", terms["z"]["coriolis"], COLORS["muted"], "//")]   # the z-equation
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), sharey=True)
for ax, rows, ttl in zip(axes, (x_terms, z_terms), ("horizontal (x) equation", "vertical (z) equation")):   # one panel per equation
    ax.bar(range(len(rows)), [r[1] for r in rows], color=[r[2] for r in rows], hatch=[r[3] for r in rows])   # one bar per term, coloured by kind
    ax.set_xticks(range(len(rows)), [r[0] for r in rows], rotation=35, ha="right", fontsize=8)   # the names of the terms
    ax.set(yscale="log", title=ttl)
axes[0].set_ylabel("size of the term [m s$^{-2}$]")
fig.suptitle("Two tall pairs: Coriolis–pressure in the horizontal, weight–pressure in the vertical")
plt.show()
""",
    see="The estimated size of every term of the horizontal and vertical thin-shell equations for a mid-latitude weather system, "
        "on a logarithmic axis. Hatched grey bars are the pieces of the Coriolis acceleration that were dropped.",
    read="In the horizontal, Coriolis (teal) and pressure (orange) tower over everything else; in the vertical, pressure and "
         "weight do. These two tall pairs are the two balances of the next blocks: geostrophic (`C02`) and hydrostatic.",
    change="…$L$ shrinks to 10 km (a thunderstorm)? Then Ro ≈ 17: the acceleration bar overtakes Coriolis and rotation stops "
           "mattering.")
note("N16 [C]", r"""
**The f-plane:** hold $f$ at its value at the central latitude, $f=f_0=2\Omega\sin\theta_0$ (`GFD.f_plane(lat)`). Used in
`C04`–`C12`, `C14`, `C16`.""")
note("N17 [B]", r"""
**The β-plane:** keep the first term of a Taylor expansion in the northward distance $y$:
$f=f_0+\beta y$ with $\beta\equiv(df/dy)_{\theta_0}=(df/d\theta\cdot d\theta/dy)_{\theta_0}=2\Omega\cos\theta_0/R$ (13.10),
since $dy=R\,d\theta$. The result is stated, not derived; its error is measured in the next figure.""")
code(r"""
print(f"f-plane value at 35 N: f0 = {GFD.f_plane(lat):.4e} 1/s")                       # f held constant
print(f"beta at 35 N = {GFD.beta_parameter(lat):.4e},  at 12 N = {GFD.beta_parameter(inp['lat_rossby']):.4e}  1/(m s)")   # Eq. (13.10)
for y_km in (500.0, 1000.0, 2000.0, -2000.0):                                          # northward distances from 35 N [km]
    print(f"  y = {y_km:+6.0f} km: beta-plane error = {100*GFD.beta_plane_error(y_km*1e3, lat):+.2f} %")   # (f_beta - f_exact)/f_exact
""", explain=r"""
**What does this show?** The straight-line approximation of $f$ is good to a fraction of a per cent within 500 km of the
central latitude and to a few per cent within 2000 km; the error grows like $y^2$ and is larger on the equatorward side.""")
nb.plotly(r"""
y_bp = np.linspace(-2.0e6, 2.0e6, 33)                          # northward distance from the central latitude [m]

def f1(lat0_deg):                                              # curves for one central latitude [degrees]
    lat0 = np.deg2rad(lat0_deg)                                # in radians
    f_exact = 2*GFD.OMEGA_EARTH*np.sin(lat0 + y_bp/GFD.EARTH_RADIUS_MEAN)   # the true f at each y [1/s]
    f_beta = GFD.beta_plane(y_bp, lat0)                        # Eq. (13.10): f0 + beta y
    return {"exact f = 2Ω sin(θ0 + y/R)": (y_bp/1e3, f_exact*1e4),
            "β-plane f0 + βy": (y_bp/1e3, f_beta*1e4),
            "10 × (β-plane − exact)": (y_bp/1e3, 10*(f_beta - f_exact)*1e4)}

figF1 = slider_figure(f1, "central latitude θ0", np.arange(5.0, 75.1, 5.0), unit="°", xlabel="northward distance y [km]",
                      ylabel="f [10⁻⁴ s⁻¹]", title="The β-plane is a straight line through the true curve", active=6)   # 15 slider positions, starting at 35 degrees
figF1.show()
""", explain=r"""
1. `f1` computes the exact $f$ along a north–south line and the β-plane straight line through the same central point.
2. The third trace is their difference magnified ten times, so that it can be seen on the same axis.""")
nb.figure_notes(
    see="The true Coriolis parameter (a piece of a sine curve), the β-plane straight line through its central point, and ten "
        "times their difference.",
    read="The difference is a parabola: the error grows like $y^2$. At low central latitudes the line is steep ($\\beta$ large) "
         "while $f_0$ itself is small, so the *relative* error is largest there.",
    change="…the central latitude were 75°? The curve flattens (β → 0 toward the pole) and bends more: the β-plane is least "
           "useful where β matters least.")
nb.md(r"""
#### Reading the three equations

$$ \underbrace{\frac{Du}{Dt}}_{\text{acceleration}}\underbrace{-fv}_{\text{Coriolis}}=\underbrace{-\frac1{\rho_0}\frac{\partial p}{\partial x}}_{\text{pressure}}+\underbrace{F_x}_{\text{friction}},\qquad
\frac{Dv}{Dt}+fu=-\frac1{\rho_0}\frac{\partial p}{\partial y}+F_y,\qquad
\frac{Dw}{Dt}=-\frac1{\rho_0}\frac{\partial p}{\partial z}-\frac{g\rho}{\rho_0}+F_z \qquad\text{(13.9)} $$

Rotation enters through $f$ alone, and only in the horizontal; the two Coriolis terms have opposite signs, which is what makes
things turn. Colour code used from here on: acceleration purple, Coriolis teal, pressure orange, friction rose, density blue.

**What would change if…** the flow is slow and wide, so that the acceleration and friction bars are negligible? Only the two
tall bars are left in each horizontal equation — geostrophic balance (`C02`).""")


# =====================================================================================================================
# A.5  §13.5 — C02 geostrophic balance, C03 thermal wind
# =====================================================================================================================
nb.section("13.5", "Geostrophic Flow", intro=r"""
**What is this section about?** The balance that rules every slow, large-scale flow: Coriolis force against pressure gradient.
Then its two consequences — a wind that changes with height wherever temperature changes horizontally, and, when it does not,
fluid that moves in rigid columns.""")
core("C02", r"Geostrophic balance: $-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$ and $fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$ (13.11)–(13.12)",
     "Why doesn't the air simply flow from high to low pressure?")
problem(r"""
Open any weather map. The wind arrows run *along* the isobars, with low pressure on their left in the northern hemisphere. A
ball on a hill rolls downhill; air on a pressure "hill" goes round it. The reason is that on a turning earth anything that
moves is pushed sideways, and the push grows with speed. The air speeds up toward low pressure, is turned, overshoots and is turned
back: released from rest, a frictionless parcel traces a row of arches along the isobar. Averaged over an arch it moves along
the isobars, at the one speed for which the sideways push cancels the pressure force. Friction, or the radiation of waves
(`C12`), removes the arches and leaves that steady drift. Oceanographers use
the same balance backwards: an altimeter measures the slope of the sea surface, and the slope gives the current.""")
idea(r"""
  (1) at rest             (2) moving, being turned           (3) on average: along the isobar
   L                        L                                  L
   ↑ pressure force         ↑ pressure force                   ↑ pressure force
   ●                        ●→  velocity                       ●━━━▶ velocity (along the isobar)
                            ↘ Coriolis (right of motion,       ↓ Coriolis = − pressure force
                              grows with speed)
""", r"""
In the northern hemisphere ($f>0$); mirror for $f<0$. The Coriolis force is always at right angles to the motion, so it can
cancel the pressure force only when the motion is at right angles to the pressure force — along the isobars. A parcel
released from rest does not stop there: it overshoots (at its fastest it moves at twice the balanced speed), and its
*velocity* circles the balanced value once per inertial period for ever, unless something damps it. Its *path* is then a row
of cycloid arches with cusps — no closed loops — drifting along the isobar. Panel (3) is the *average* over an arch.""")
remind("C02")
P("P309", "reading a pressure map: isobars, the pressure-gradient force points from high to low, tight spacing = strong force", r"""
An isobar joins points of equal pressure, like a height contour on a hiking map. The pressure-gradient force per unit mass is
$-\nabla p/\rho_0$: it points straight across the isobars from high to low, and it is large where the lines are close.""", code=r"""
dp, dn, rho0 = 400.0, 3.0e5, 1.2      # 4 hPa between isobars 300 km apart; air density [kg/m^3]
print(dp / dn / rho0)                 # force per unit mass: 1.1e-3 m/s^2 (compare g = 9.8)
""")
nb.recap("R13", "The Rossby number", r"""
It is $\mathrm{Ro}=\dfrac{U^2/L}{fU}=\dfrac{U}{fL}$ (13.13): the acceleration of the flow divided by the Coriolis acceleration.
Small Ro = rotation rules. ⚠️ **Trap T17:** Chapter 4's `SIM.rossby_number(U, Omega, l)` is $U/(2\Omega l)$; the two differ
by $\sin\theta$.""", where="Ch. 4 §4.7")
nb.current_core = "C02"
code(r"""
Ro_here = GFD.rossby_number(14.0, f35, 1.4e6)                  # Eq. (13.13): U/(f L) for U = 14 m/s, L = 1400 km at 35 N
Ro_ch4 = SIM.rossby_number(14.0, GFD.OMEGA_EARTH, 1.4e6)       # Chapter 4's version U/(2 Omega L)
print(f"Ro = U/(fL) = {Ro_here:.3f};  Ch. 4's U/(2 Omega L) = {Ro_ch4:.4f};  ratio = {Ro_ch4/Ro_here:.3f} = sin(35 deg)")   # the two definitions
""", explain=r"""
**What does this show?** Both numbers are small, so rotation dominates a weather system; they differ by the factor
$\sin 35°=0.574$, because this chapter uses $f$ where Chapter 4 used $2\Omega$.""")
note("N22 [B]", r"""
**The Ekman number** $E=\dfrac{\rho\nu U/L^2}{\rho fU}=\dfrac{\nu}{fL^2}$ (13.18) is friction over Coriolis. Geostrophy needs
both $\mathrm{Ro}\ll1$ and $E\ll1$.""")
code(r"""
E_vert = GFD.ekman_number(7.0, f35, 1000.0)                    # vertical eddy friction over a depth of 1 km (nu_v = 7 m^2/s)
E_horiz = GFD.ekman_number(3.0e4, f35, 1.4e6)                  # horizontal eddy friction over 1400 km (nu_H = 3e4 m^2/s)
print(f"E (vertical, lowest km) = {E_vert:.3f};  E (horizontal) = {E_horiz:.1e}")   # the first is not small
""", explain=r"""
**What does this show?** Horizontal friction is utterly negligible, but within the lowest kilometre vertical friction is
about a tenth of the Coriolis force — not small. That is the Ekman layer of `C06`. In `C04` the same number appears as
$(\delta/L)^2/2$.""")
D("D03", ref="13.11")
note("N18 [B]", r"""
**The wind blows along the isobars, and pressure is its stream function** (steps 5–6 of `D03`): $\mathbf u\cdot\nabla p=0$, and
for constant $f$ the function $\psi=p/(f\rho_0)$ gives $u=-\partial\psi/\partial y$, $v=\partial\psi/\partial x$
(`GFD.geostrophic_streamfunction`).""")
trap("T14", r"""
Chapters 4 and 11 used $u=+\partial\psi/\partial y$. Here the sign is the other way, so that $\psi$ is high where $p$ is high
(for $f>0$).""")
nb.worked_example("the wind from two isobars", r"""
Isobars 4 hPa (400 Pa) apart lie 400 km apart; $\rho_0\approx1.25$ kg/m³; $f\approx10^{-4}$ s⁻¹ (about 43° N).

1. Pressure-gradient force per unit mass: $400/(4\times10^{5}\times1.25)=8\times10^{-4}$ m/s².
2. Balance $f\,U=\frac1{\rho_0}\lvert\nabla p\rvert$: $U=8\times10^{-4}/10^{-4}=8$ m/s.
3. Direction: along the isobars, low pressure on the left (northern hemisphere, $f>0$; on the right for $f<0$).
4. Halve the spacing → 16 m/s.
5. The same isobars at 20° N ($f\approx5\times10^{-5}$ s⁻¹): 16 m/s — and at the equator the formula gives infinity: geostrophy
   fails there.""")
code(r"""
u_n, v_n = GFD.geostrophic_velocity(0.0, 400.0/3.0e5, f35, 1.2)     # pressure rising northward by 4 hPa in 300 km, at 35 N
u_s, v_s = GFD.geostrophic_velocity(0.0, 400.0/3.0e5, -f35, 1.2)    # the same pressure field at 35 S (f < 0)
u_oc, v_oc = GFD.geostrophic_from_height(np.outer([0.0, 0.1, 0.2], np.ones(3)), 1.0e5, 1.0e5, f35)   # a sea surface rising northward by 0.1 m per 100 km
print(f"35 N: u = {u_n:+.2f} m/s ({ch13.wind_from_to(u_n, v_n)['wind_name']} wind; high pressure on its right)")   # Eq. (13.12)
print(f"35 S: u = {u_s:+.2f} m/s ({ch13.wind_from_to(u_s, v_s)['wind_name']} wind; high pressure on its left)")    # the mirror
print(f"ocean: surface slope 0.1 m per 100 km -> current u = {u_oc[1, 1]:+.3f} m/s")   # Eq. (13.116) used as a current meter
""", explain=r"""
1. `GFD.geostrophic_velocity(dpdx, dpdy, f, rho0)` solves the two balance equations for $(u, v)$.
2. With high pressure to the north the wind at 35° N is an easterly of about 13 m/s (low pressure on the left of the motion);
   at 35° S the same pressure field drives a westerly of the same speed.
3. `GFD.geostrophic_from_height` is the oceanographer's version: a surface slope of 0.1 m in 100 km means a current of about
   0.12 m/s.""")
scratch(r"""
# From scratch: the geostrophic wind of a low, by centred differences written with slices
xg = np.linspace(-2.0e6, 2.0e6, 81); yg = np.linspace(-2.0e6, 2.0e6, 81)   # grid [m]; arrays are indexed [j, i] = (y, x)
dxg, dyg, rho0 = xg[1] - xg[0], yg[1] - yg[0], 1.2             # grid spacings [m] and air density [kg/m^3]
pc = ch13.pressure_centre(xg, yg, dp=400.0, R_c=6.0e5, lat_rad=lat, kind="low")   # a 4 hPa low at 35 N
p = pc["p"]                                                    # pressure anomaly [Pa]
dpdx = np.zeros_like(p); dpdy = np.zeros_like(p)               # arrays for the two pressure gradients
dpdx[:, 1:-1] = (p[:, 2:] - p[:, :-2])/(2*dxg)                 # centred difference along x (columns)
dpdy[1:-1, :] = (p[2:, :] - p[:-2, :])/(2*dyg)                 # centred difference along y (rows)
u_mine = -dpdy/(rho0*f35)                                      # Eq. (13.12): f u = -(1/rho0) dp/dy
v_mine = dpdx/(rho0*f35)                                       # Eq. (13.11): -f v = -(1/rho0) dp/dx
u_lib, v_lib = GFD.geostrophic_from_field(p, dxg, dyg, f35, rho0)   # the tested function
assert np.allclose(u_mine[1:-1, 1:-1], u_lib[1:-1, 1:-1]) and np.allclose(v_mine[1:-1, 1:-1], v_lib[1:-1, 1:-1])   # same interior values
cross = np.abs(u_mine*dpdx + v_mine*dpdy)[1:-1, 1:-1].max()    # u . grad p: zero if the wind is along the isobars
print(f"✓ hand-written differences agree with GFD.geostrophic_from_field; max |u·∇p| = {cross:.1e} (round-off: wind along isobars)")   # visible confirmation
""", r"""
The two balance equations, with the derivatives written as array slices, reproduce the library; and $\mathbf u\cdot\nabla p$ is
zero to round-off — the wind crosses no isobar.""")
nb.recap("R14", "Highs and lows", r"""
Round a low the geostrophic wind turns counter-clockwise in the northern hemisphere ($f>0$) and clockwise in the southern
($f<0$); round a high the other way. The pressure force $-\nabla p/\rho_0$ points inward for a low; the Coriolis force
$-f\,\mathbf e_z\times\mathbf u$ must point outward.""", where="Ch. 4 §4.7")
nb.current_core = "C02"
fig(r"""
fig, axes = plt.subplots(2, 2, figsize=(8.6, 8.0), sharex=True, sharey=True)
for row, (L, hemi) in enumerate(((lat, "35° N (f > 0)"), (-lat, "35° S (f < 0)"))):   # one row per hemisphere
    for col, kind in enumerate(("low", "high")):               # one column per kind of centre
        pc4 = ch13.pressure_centre(x, y, dp=400.0, R_c=6.0e5, lat_rad=L, kind=kind)   # pressure anomaly and its geostrophic wind
        ax = axes[row, col]                                    # this panel
        ax.contour(x/1e3, y/1e3, pc4["p"]/100, 7, colors=COLORS["orange"], linewidths=1.1)   # isobars (orange = pressure)
        ax.quiver(x[::4]/1e3, y[::4]/1e3, pc4["u"][::4, ::4], pc4["v"][::4, ::4], color=COLORS["teal"], scale=70)   # geostrophic wind, every fourth point
        sense = GFD.hemisphere(L)["cyclonic"] if kind == "low" else ("clockwise" if L > 0 else "counter-clockwise")   # sense of rotation
        ax.set(title=f"{hemi}, {kind}: {sense}", aspect="equal")
for ax in axes[1]:
    ax.set_xlabel("x east [km]")
for ax in axes[:, 0]:
    ax.set_ylabel("y north [km]")
plt.show()
""",
    see="Four panels: a low and a high, at 35° N (top) and 35° S (bottom). Isobars are orange, the geostrophic wind teal; "
        "each title names the sense of rotation.",
    read="In the top row low pressure is on the left of every arrow; in the bottom row it is on the right. Lows are cyclonic "
         "in both hemispheres — but cyclonic means counter-clockwise in the north and clockwise in the south.",
    change="…the pressure difference doubles? Every arrow doubles in length; no arrow turns.")
nb.md(r"""
#### Where geostrophy fails, and how it is set up

(The verbal part of note **N19**; its quantitative part is the adjustment problem of `C12`.) Geostrophy fails near the equator
($f\to0$), where friction matters ($E$ not small), and when the flow changes in less than about a day (the acceleration is not
small). And how is the balance reached? Release a parcel from rest in a uniform pressure gradient and it does not go to the
low — it swings from side to side of the isobar in arches. `GFD.parcel_adjust` is the exact solution of $dV/dt+(r+if)V=-G$ for the complex velocity $V=u+iv$, with $G$
the pressure-gradient acceleration and $r$ a linear drag (ours).""")
fig(r"""
G = (0.0 + 1j*400.0/3.0e5)/1.2                                 # complex pressure-gradient acceleration (1/rho0)(dp/dx + i dp/dy) [m/s^2]
T_i = GFD.inertial_period(f35)                                 # inertial period at 35 N [s]
t_p = np.linspace(0.0, 10*T_i, 4001)                           # ten inertial periods [s]
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.8))
for ax, f_, ttl in ((a, f35, "35° N"), (b, -f35, "35° S")):    # one panel per hemisphere
    for r_, col, lab in ((0.0, COLORS["accent"], "no drag"), (0.3*abs(f_), COLORS["rose"], "drag r = 0.3 |f|")):   # two parcels
        V = GFD.parcel_adjust(t_p, G, f_, r=r_)                # complex velocity u + iv of the released parcel [m/s]
        Z = cumulative_trapezoid(V, t_p, initial=0)            # its position x + iy by integrating the velocity [m]
        ax.plot(Z.real/1e3, Z.imag/1e3, color=col, label=lab)  # the path
    ax.set(title=f"{ttl}: arches that drift along the isobars", xlabel="x east [km]", ylabel="y north [km]")
    ax.legend(fontsize=8)
fig.suptitle("pressure rises toward +y (north), so the pressure force points toward −y", fontsize=9)
plt.show()
V0 = GFD.parcel_adjust(t_p, G, f35, r=0.0)                     # the frictionless parcel at 35 N again
drift = np.trapezoid(V0, t_p)/t_p[-1]                          # its mean velocity over ten inertial periods
Vg = complex(*GFD.geostrophic_velocity(0.0, 400.0/3.0e5, f35, 1.2))   # the geostrophic velocity for the same gradient
print(f"mean drift of the frictionless parcel = {drift.real:+.2f} {drift.imag:+.2f}i m/s; geostrophic velocity = {Vg.real:+.2f} {Vg.imag:+.2f}i m/s")   # they agree
""",
    see="Paths of a parcel released from rest in a uniform pressure gradient (pressure increasing northward, so the pressure "
        "force points south). Without drag (purple) it traces a row of arches, coming to rest for an instant at each cusp; with "
        "drag (rose) the wiggles die out within a few hundred kilometres and the path becomes a straight line.",
    read="Each arch is one inertial circle added to a steady drift *along* the isobars (a cycloid) — and the printed line shows "
         "that the mean drift is exactly the geostrophic velocity. Geostrophy is the average of the motion; the arches are what "
         "is left of the start. With drag, the parcel ends up moving steadily at an angle to the isobars, toward low pressure "
         "(a preview of `C06`).",
    change="…$f<0$ (right panel)? The parcel is turned to the left instead: the arches and the drift go east instead of west.")
explainer("geostrophic_balance", "Why doesn't air flow straight from high to low pressure?",
          "releasing a parcel and watching the Coriolis arrow grow and turn with the velocity shows how the balance is approached — "
          "the parcel overshoots and swings about it in arches, and only with drag do the bars settle; a static map shows only the end state.",
          ["Press play with the default low and watch the teal Coriolis arrow swing round and overshoot; with drag switched on it settles opposite the orange pressure arrow.",
           "Switch to S: the circulation reverses at once.",
           "Drag the latitude toward 5°: Ro climbs past 1 and the parcel cuts across the isobars.",
           "Turn on drag: the wind settles at an angle to the isobars, toward low pressure — a preview of `C06`."])
whatif(r"""
…the density, and so the pressure pattern, differs from one height to the next? Then the geostrophic wind differs too: the
thermal wind (`C03`).""")

# ---------------------------------------------------------------------------------------------------------------------
core("C03", r"Thermal wind: $\dfrac{\partial v}{\partial z}=-\dfrac{g}{\rho_0f}\dfrac{\partial\rho}{\partial x}$ and $\dfrac{\partial u}{\partial z}=\dfrac{g}{\rho_0f}\dfrac{\partial\rho}{\partial y}$ (13.15)",
     "Why is there a jet stream high up, and why is it westerly in both hemispheres?")
problem(r"""
Airliners flying east across the Atlantic ride a river of air moving at 150 km/h or more; the wind at the ground below is far
gentler. The jet sits above the place where the temperature changes fastest from south to north. That is no coincidence. Cold
air is dense, so pressure falls faster with height in the cold column than in the warm one; the pressure difference between
the columns therefore grows with height, and with it the geostrophic wind. The temperature gradient does not set the wind — it
sets how the wind *changes with height*.

*Two words.* A state in which surfaces of constant pressure and constant density coincide is **barotropic** (no thermal wind);
one in which they cross is **baroclinic** (the wind changes with height).""")
idea(r"""
   pole (cold)                         equator (warm)
   ───────────────── p = 300 hPa ──────────────╱      ← the upper pressure surface slopes down toward the pole
     short, dense column          tall, light column
   ───────────────── p = 1000 hPa ─────────────       ← flat at the ground
        ⊙ weak wind at the ground,  ⊙⊙⊙ strong wind aloft (out of the page = eastward)
""", r"""
| Hemisphere | Cold air lies to the | Sign of $f$ | Shear $\partial u/\partial z$ |
|---|---|---|---|
| northern | north | $f>0$ | westerly (eastward wind grows with height) |
| southern | south | $f<0$ | westerly again — both signs flip |""")
nb.recap("R15", "The perturbation is hydrostatic", r"""
The balance is $0=-\dfrac{\partial p}{\partial z}-g\rho$ (13.14). ⚠️ **Trap T3 again:** the symbols $p$ and $\rho$ are the perturbations from the state
of rest (primes dropped).""", where=r"Ch. 1 hydrostatics $dp/dz=-\rho g$ (1.8); Ch. 7 $p'=\rho g\eta$ (7.52)")
nb.current_core = "C03"
remind("C03", "the only two tools `D04` needs are the first two")
D("D04", ref="13.15")
note("N20 [B]", r"""
**Not in the book — ours** (the book states it in words): the last two steps of `D04` put the thermal wind in temperature form.
With $\rho'/\rho_0=-\alpha T'$, $\dfrac{\partial u}{\partial z}=-\dfrac{g\alpha}{f}\dfrac{\partial T}{\partial y}$ and
$\dfrac{\partial v}{\partial z}=\dfrac{g\alpha}{f}\dfrac{\partial T}{\partial x}$, where $\alpha=1/T_0$ for a perfect gas and
$T$ is the potential temperature in a deep layer.""")
nb.worked_example("a jet from a temperature gradient", r"""
Take $g\approx10$ m/s², $T_0=250$ K so $g\alpha=g/T_0=0.04$ m s⁻² K⁻¹, $f=10^{-4}$ s⁻¹, and let temperature fall northward by
1 K per 100 km: $\partial T/\partial y=-10^{-5}$ K/m.

1. Shear: $\frac{\partial u}{\partial z}=-\frac{g\alpha}{f}\frac{\partial T}{\partial y}=-(0.04/10^{-4})\times(-10^{-5})=+4\times10^{-3}$ s⁻¹ = 4 m/s per km.
2. With calm air at the ground, the wind at 10 km is 40 m/s from the west.
3. Southern hemisphere: $\partial T/\partial y=+10^{-5}$ K/m (cold toward the south pole) and $f=-10^{-4}$ s⁻¹: the two sign
   changes cancel, again +4 m/s per km.
4. No gradient → no shear.""")
code(r"""
dudz, dvdz = GFD.thermal_wind_from_temperature(0.0, -7.0e-6, f35, alpha=1/280.0)    # temperature falling northward by 0.7 K per 100 km at 35 N
dudz_s, _ = GFD.thermal_wind_from_temperature(0.0, +7.0e-6, -f35, alpha=1/280.0)    # 35 S: cold toward the south pole, f < 0
du_oc, dv_oc = GFD.thermal_wind_shear(0.0, 1.0e-5, f35, 1027.0)                     # Eq. (13.15) for an ocean front: density rising northward by 1 kg/m^3 per 100 km
print(f"35 N: du/dz = {dudz:.3e} 1/s  ->  wind at 9 km above calm ground = {dudz*9000.0:.1f} m/s (westerly)")   # the jet
print(f"35 S: du/dz = {dudz_s:.3e} 1/s  (the same: both signs flipped)")             # westerly in both hemispheres
print(f"ocean front: du/dz = {du_oc:.2e} 1/s = {du_oc*1000:.2f} m/s per 1000 m")     # the density form
""", explain=r"""
1. `GFD.thermal_wind_from_temperature(dTdx, dTdy, f, alpha=…)` is our temperature form; `alpha` must be passed **by name**,
   because the letter α means three different things in this chapter and a positional argument would invite a mix-up.
2. A modest gradient of 0.7 K per 100 km gives a 26 m/s westerly at 9 km.
3. The southern twin gives exactly the same shear: the temperature gradient and $f$ both change sign.
4. `GFD.thermal_wind_shear` is the book's density form $\frac{\partial u}{\partial z}=\frac{g}{\rho_0f}\frac{\partial\rho}{\partial y}$ *(13.15)*; an ocean front gives about 1 m/s per kilometre of depth.""")
scratch(r"""
# From scratch: the shear in one line, then the wind profile by a hand-written running integral
shear_mine = -G0*(1/280.0)*(-7.0e-6)/f35                       # du/dz = -(g alpha/f) dT/dy [1/s]
assert np.isclose(shear_mine, dudz)                            # same as GFD.thermal_wind_from_temperature
z_tw = np.linspace(0.0, 9000.0, 91)                            # heights from the ground to 9 km [m]
drho_dy = 1.2*(1/280.0)*7.0e-6*np.exp(-z_tw/8000.0)            # a density gradient that weakens with height [kg/m^4]
U_loop = np.zeros_like(z_tw)                                   # wind profile, calm at the ground
for k in range(1, len(z_tw)):                                  # trapezoid rule, one layer at a time
    s0 = G0*drho_dy[k-1]/(1.2*f35)                             # shear at the bottom of the layer, Eq. (13.15)
    s1 = G0*drho_dy[k]/(1.2*f35)                               # shear at the top of the layer
    U_loop[k] = U_loop[k-1] + 0.5*(s0 + s1)*(z_tw[k] - z_tw[k-1])   # add the layer's contribution
assert np.allclose(U_loop, GFD.thermal_wind_integrate(z_tw, drho_dy, f35, 1.2))   # same as the library's integral
print(f"✓ shear {shear_mine:.3e} 1/s and the integrated profile agree with the library (U at 9 km = {U_loop[-1]:.1f} m/s)")   # visible confirmation
""", r"""
The thermal-wind relation is one multiplication; turning a shear into a wind is a running integral from a level where the wind
is known. Both hand-written versions match `GFD.thermal_wind_from_temperature` and `GFD.thermal_wind_integrate`.""")
fig(r"""
y_j = np.linspace(-4.0e6, 4.0e6, 121)                          # north-south distance from the front [m]
z_j = np.linspace(0.0, 1.2e4, 61)                              # height [m]
sec = ch13.jet_section(y_j, z_j, dT=28.0, width=2.0e6, lat_rad=lat, alpha=1/280.0)   # temperature and thermal-wind-balanced zonal wind (ours)
fig, ax = plt.subplots(figsize=(7.6, 4.2))
cf = ax.contourf(y_j/1e3, z_j/1e3, sec["U"], 12, cmap="Purples")   # zonal wind [m/s], filled
cs = ax.contour(y_j/1e3, z_j/1e3, sec["T"], 12, colors=COLORS["blue"], linewidths=0.9)   # isotherms (blue = temperature)
ax.clabel(cs, fmt="%.0f K", fontsize=7)                        # label the isotherms
jz, jy = np.unravel_index(np.argmax(sec["U"]), sec["U"].shape) # where the wind is strongest
ax.plot(y_j[jy]/1e3, z_j[jz]/1e3, "o", color=COLORS["orange"]) # mark the jet core
fig.colorbar(cf, label="eastward wind u [m/s]")
ax.set(xlabel="y north [km]  (cold air to the north)", ylabel="height z [km]", title="The wind grows with height above the tightest isotherms")
plt.show()
i9 = np.argmin(abs(z_j - 9000.0)); i0 = np.argmin(abs(y_j))    # indices of z = 9 km and y = 0
print(f"wind at 9 km above the front: {sec['U'][i9, i0]:.1f} m/s;  dT/dy there: {sec['dTdy'][i0]:.1e} K/m")   # compare with the code cell above
""",
    see="A north–south section through a front at 35° N: isotherms (blue) dipping toward the cold side, and the eastward wind "
        "(purple shading) that is in thermal-wind balance with them. The orange dot is the strongest wind.",
    read="Follow one isotherm: where it slopes most steeply, the wind contours are closest in the vertical. The wind maximum "
         "sits above the tightest isotherms, at the top of the section.",
    change="…the surface air were already blowing at 5 m/s? Every contour would shift by 5 m/s — the thermal wind fixes the "
           "shear, not the wind.")
nb.plotly(r"""
z_f3 = np.linspace(0.0, 1.2e4, 13)                             # height [m]

def f3(dT):                                                    # wind profiles above the front for one temperature contrast dT [K]
    out = {}
    for sgn, nm in ((+1, "35° N: cold air to the north, f > 0"), (-1, "35° S: cold air to the south, f < 0")):   # both hemispheres
        dTdy = -sgn*dT/(2*2.0e6)                               # temperature gradient at the centre of a front 2000 km wide [K/m]
        drho = -1.2*(1/280.0)*dTdy*np.ones_like(z_f3)          # density gradient from rho' = -rho0 alpha T' [kg/m^4]
        out[nm] = (GFD.thermal_wind_integrate(z_f3, drho, sgn*f35, 1.2), z_f3/1e3)   # Eq. (13.15) integrated up from calm ground
    return out

figF3 = slider_figure(f3, "temperature contrast across the front", np.arange(0.0, 40.1, 2.0), unit="K", xlabel="eastward wind u [m/s]",
                      ylabel="height z [km]", title="More contrast, more shear — and westerly in both hemispheres",
                      modes={"35° S: cold air to the south, f < 0": "markers"})   # 21 slider positions
figF3.show()
""", explain=r"""
1. `f3` builds the wind profile above the centre of a front for one temperature contrast, in both hemispheres.
2. The southern profile is drawn with markers so that you can see it lies exactly on the northern line.""")
nb.figure_notes(
    see="The eastward wind against height above the centre of the front, for the northern (line) and southern (markers) "
        "hemisphere.",
    read="The profile is a straight line whose slope is the shear; it steepens in proportion to the temperature contrast. "
         "Markers and line coincide: the jet is westerly in both hemispheres.",
    change="…the contrast is zero? The profile is vertical: no change of wind with height at all (the Taylor–Proudman limit "
           "below).")
note("N21 [B]", r"""
**A rapidly rotating tank of uniform density.** Geostrophy with $f=2\Omega$ in a homogeneous fluid:
$-2\Omega v=-\dfrac1\rho\dfrac{\partial p}{\partial x}$ (13.16) and $2\Omega u=-\dfrac1\rho\dfrac{\partial p}{\partial y}$ (13.17).""")
slip(3, r"$-2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$",
     r"$2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$ (it is the second geostrophic equation with $f=2\Omega$; with the printed sign the next step of the book cannot be reached).")
code(r"""
print("corrected pair leads to dw/dz = 0 ->", ch13.taylor_proudman_sympy()["ok"])                # True
print("printed pair leads to dw/dz = 0   ->", ch13.taylor_proudman_sympy(printed=True)["ok"])    # False: the printed sign fails
""", explain=r"""
**What does this show?** With the corrected sign, cross-differentiation and continuity give "no vertical stretching"; with the
printed sign they do not — the engine reaches a different (wrong) statement.""")
D("D05", ref="13.21")
note("N23 [B]", r"""
**No stretching** (step 4 of `D05`): $\dfrac{\partial w}{\partial z}=0$ (13.19).""")
note("N24 [B]", r"""
**No shear** (step 6 of `D05`): $\dfrac{\partial v}{\partial z}=\dfrac{\partial u}{\partial z}=0$ (13.20).""")
note("N25 [B]", r"""
**The Taylor–Proudman theorem:** $\partial\mathbf u/\partial z=0$ (13.21), under four hypotheses — steady, $\mathrm{Ro}\ll1$,
$E\ll1$, uniform density. It is the thermal wind with nothing to drive a shear.""")
note("N26 [B]", r"""
**Taylor's experiment, as our sketch.** A short obstacle on the floor of a rapidly rotating tank acts as if it reached the
surface: dye divides ahead of the column above it and goes round. (The tank's dimensions are not reproduced.)""")
fig(r"""
figT = draw_taylor_column()                                    # our sketch: the same two-dimensional flow at every height, column shaded
plt.show()
res = GFD.taylor_proudman_residual(lambda xx, yy, zz: (-yy, xx, 0.0), (0.3, 0.2, 0.5))   # a rigid rotation about the vertical axis
print("d(u, v, w)/dz of a columnar flow:", {k: round(float(v_), 12) for k, v_ in res.items()})   # all three are zero
""",
    see="A side view of a rapidly rotating tank (rotation axis $\\Omega$ vertical): a short obstacle on the floor (grey), the "
        "column of fluid above it (shaded, 'moves with the obstacle'), and the oncoming flow drawn as equal arrows at three "
        "heights.",
    read="The arrows are the same at every height — nothing varies along the rotation axis — so the flow must go round the "
         "whole shaded column, not just round the obstacle: its 'shadow' extends to the surface. The printed line checks a "
         "columnar flow numerically — all three vertical derivatives vanish.",
    change="…the fluid were stratified? Density gradients would allow a vertical shear (the thermal wind) and the column would "
           "fade with height.")
explainer("thermal_wind", "Why is there a jet stream, and why is it westerly in both hemispheres?",
          "dragging the temperature contrast tilts the isotherms and grows the shear at once; flipping the hemisphere flips both "
          "$f$ and the gradient and the jet stays westerly — three coupled changes a single section cannot show.",
          ["Drag the contrast to zero: every wind profile collapses to a vertical line (the Taylor–Proudman limit).",
           "Switch N ↔ S and watch which two signs change.",
           "Click a point in the section to see its local gradient and shear worked out.",
           "Switch to the ocean mode: the same relation in density, with a front instead of a jet."])
whatif(r"""
…friction is not negligible? Near the sea surface and near the ground it is not; the balance gains a third force and the flow
turns — the Ekman layers (`C04`–`C06`). (The sheared flow of this block is also the starting state of the Eady problem, `C16`.)""")


# =====================================================================================================================
# A.6  §13.6 — C04 the surface Ekman spiral, C05 Ekman transport
# =====================================================================================================================
nb.section("13.6", "Ekman Layer at a Free Surface", intro=r"""
**What is this section about?** What a steady wind does to the top of the ocean. The answer is a thin layer in which the
current turns and weakens with depth — and whose total transport is at right angles to the wind. ⚠️ In this section $z=0$ is
the sea surface and the ocean is $z<0$.""")
core("C04", r"The surface Ekman spiral: $\dfrac{d^2V}{dz^2}=\dfrac{if}{\nu_v}V$, $V\equiv u+iv$ (13.27), with $V=A\,e^{(1+i)z/\delta}+B\,e^{-(1+i)z/\delta}$, $\delta=\sqrt{2\nu_v/f}$ (13.28)–(13.29)",
     "A steady wind blows over the sea. Which way does the water go, and how deep does the motion reach?")
problem(r"""
Nansen noticed that Arctic ice drifts not downwind but well to the right of the wind; Ekman explained it in 1905.
The wind drags the top layer of water; that layer is turned by the Coriolis force and drags the layer below, which is turned
a little more, and so on downward. The result is a staircase of currents that rotates and fades with depth: a spiral.""")
idea(r"""
   wind stress τ  ─────────▶
  ═══════════════════════════  z = 0 (sea surface)
   slab 1   pushed by the wind above, held back by slab 2, turned by Coriolis
   slab 2   pushed by slab 1,         held back by slab 3, turned by Coriolis
   slab 3   …                          each one slower and turned further
""", r"""
In the steady state the net drag on each slab (push from above minus hold from below) balances its Coriolis force. Friction
needs *something* to balance it — and there are three candidates (next recap).""")
nb.recap("R16", "Three ways to balance friction", r"""
Friction can be balanced by unsteadiness (the layer thickens as $\delta\sim\sqrt{\nu t}$), by advection (it thickens
downstream, $\delta\sim\sqrt{\nu x/U}$), or — new here — by the Coriolis force, and then the thickness
$\delta=\sqrt{2\nu_v/f}$ does **not** grow at all.""", where="Ch. 8 §8.4 (impulsively started plate), Ch. 9 §9.3 (Blasius layer)")
nb.current_core = "C04"
code(r"""
th = ch13.viscous_layer_thicknesses(0.03, t=86400.0, x=1.0e5, U=0.1, f=f60)   # nu = 0.03 m^2/s; after one day; after 100 km at 0.1 m/s; at 60 N
display(pd.DataFrame({"balanced by": ["unsteadiness (one day)", "advection (100 km at 0.1 m/s)", "Coriolis (Ekman, 60° N)"],
                      "thickness [m]": [round(th["diffusive"], 1), round(th["boundary_layer"], 1), round(th["ekman"], 1)],
                      "keeps growing?": ["yes, like √t", "yes, like √x", "no"]}))   # the three-row comparison
""", explain=r"""
**What does this show?** The same eddy viscosity gives a layer that would thicken for ever if only time or distance balanced
friction — and a layer about 22 m thick that stays put when the Coriolis force does the balancing.""")
note("N27 [B]", r"""
**The balance.** Steady, horizontally uniform, no pressure gradient: only Coriolis (teal) and vertical friction (rose) remain,
$-fv=\nu_v\dfrac{d^2u}{dz^2}$ and $fu=\nu_v\dfrac{d^2v}{dz^2}$ (13.22)–(13.23).""")
note("N28 [B]", r"""
**The conditions.** The wind stress $\tau$ acts along $x$ at the surface, $\rho\nu_v\dfrac{du}{dz}=\tau$ at $z=0$ (13.24) and
$\dfrac{dv}{dz}=0$ at $z=0$ (13.25); and the current dies out at depth, $u,v\to0$ as $z\to-\infty$ (13.26, in its corrected
form).""")
slip(4, r"the condition (13.26) as $u,v\to0$ for $z\to\infty$",
     r"$u,v\to0$ for $z\to-\infty$: the ocean lies below the surface $z=0$.")
P("P310", "hodograph: the tip of the velocity vector traced as depth (or time) changes", r"""
Plot $v$ against $u$ for every depth and join the points. The curve is the hodograph: each point is the tip of the velocity
arrow at one depth. A straight hodograph means the current keeps its direction; a spiral means it turns.""", code=r"""
s = np.linspace(0, 3, 4)                         # four depths, in units of delta
V = np.exp(-s) * np.exp(-1j*s)                   # speed decays, direction turns clockwise
print(np.round(V.real, 2), np.round(V.imag, 2))  # the points of a spiral
""")
P("P311", "complex velocity V = u + iv: multiplying by i turns a vector 90° to the left; e^{(1+i)s} decays and turns at the same rate", r"""
Pack the two velocity components into one complex number $V=u+iv$. Then $iV$ is the same arrow turned 90° counter-clockwise —
exactly what the Coriolis *acceleration* term $ifV$ on the left-hand side does (the Coriolis *force* is $-ifV$: to the right of
the motion for $f>0$) — so two coupled real equations become one complex one. A factor $e^{(1+i)s}=e^{s}e^{is}$
changes the length by $e^s$ and the direction by $s$ radians together.""", code=r"""
V = 3 + 4j                                       # u = 3, v = 4
print(1j*V)                                      # (-4+3j): turned 90 degrees to the left
print(abs(np.exp((1+1j)*(-1.0))), np.angle(np.exp((1+1j)*(-1.0))))   # 0.368 and -1.0 rad: shorter and turned
print(np.sqrt(1j))                               # (0.707+0.707j) = (1+i)/sqrt(2)
""")
remind("C04")
D("D06", ref="13.27")
note("N29 [B]", r"""
**One complex equation** (step 2 of `D06`): $\dfrac{d^2V}{dz^2}=\dfrac{if}{\nu_v}V$, $V\equiv u+iv$ (13.27).""")
note("N30 [B]", r"""
**The general solution and the thickness** (steps 3–6 of `D06`): $V=A\,e^{(1+i)z/\delta}+B\,e^{-(1+i)z/\delta}$ with
$\delta=\sqrt{2\nu_v/f}$ (13.28)–(13.29).""")
trap("T6", r"""
The thickness $\delta$ is the depth over which the current falls by a factor $e$ (and turns by one radian). Oceanographers'
"Ekman depth" is $\pi\delta$, where the current first points opposite to the surface current
(`GFD.ekman_depth(nu_v, f, convention='efold' or 'pi')`). The laminar cousin is Chapter 8's Stokes layer,
$u=Ue^{-y\sqrt{\omega/2\nu}}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38), with thickness $\sqrt{2\nu/\omega}$: replace the
oscillation frequency by $f$.""")
nb.worked_example("thickness and surface speed", r"""
Take $\nu_v=0.05$ m²/s, $f=10^{-4}$ s⁻¹, $\tau=0.1$ N/m², $\rho=1000$ kg/m³.

1. $\delta=\sqrt{2\nu_v/f}=\sqrt{0.1/10^{-4}}=\sqrt{1000}\approx31.6$ m.
2. Ekman depth $\pi\delta\approx99$ m.
3. Surface speed $V_0=\dfrac{\tau/\rho}{\sqrt{f\nu_v}}=10^{-4}/\sqrt{5\times10^{-6}}=10^{-4}/2.24\times10^{-3}\approx0.045$ m/s.
4. Direction: 45° to the right of the wind (northern hemisphere, $f>0$; to the left for $f<0$), so each component has the
   magnitude $0.045/\sqrt2\approx0.032$ m/s.
5. At $z=-\delta$: speed × $e^{-1}$ = 0.016 m/s, turned a further 57°.""")
code(r"""
tau, nu_o, rho_o = inp["tau"], inp["nu_v_ocean"], inp["rho_ocean"]   # wind stress 0.07 N/m^2 (eastward), nu_v = 0.03 m^2/s, rho = 1027 kg/m^3
delta = GFD.ekman_depth(nu_o, f60)                             # Eq. (13.29): e-folding thickness sqrt(2 nu_v/f) at 60 N [m]
z_e = np.linspace(-5*delta, 0.0, 201)                          # depths from 5 delta below the surface up to the surface [m]
u_e, v_e = GFD.ekman_surface(z_e, tau, 0.0, rho_o, nu_o, f60)  # the spiral at 60 N [m/s]
u_s, v_s = GFD.ekman_surface(z_e, tau, 0.0, rho_o, nu_o, -f60) # the spiral at 60 S
u_p, v_p = ch13.ekman_surface_printed(z_e, tau, rho_o, nu_o, f60)   # the book's cos / sin form (valid for f > 0 only)
i_pi = np.argmin(abs(z_e + np.pi*delta))                       # index of the depth z = -pi delta
print(f"delta = {delta:.2f} m;  Ekman depth pi*delta = {GFD.ekman_depth(nu_o, f60, convention='pi'):.2f} m")   # the two thicknesses (trap T6)
print(f"60 N surface current (u, v) = ({u_e[-1]:+.5f}, {v_e[-1]:+.5f}) m/s, speed {np.hypot(u_e[-1], v_e[-1]):.5f} m/s: 45 deg to the RIGHT of the wind")   # f > 0
print(f"60 S surface current (u, v) = ({u_s[-1]:+.5f}, {v_s[-1]:+.5f}) m/s: 45 deg to the LEFT of the wind")   # f < 0
print(f"speed at z = -pi*delta = {100*np.hypot(u_e[i_pi], v_e[i_pi])/np.hypot(u_e[-1], v_e[-1]):.1f} % of the surface speed, pointing the opposite way")   # e^-pi
print(f"printed cos/sin form vs ours (f > 0): max difference {max(abs(u_p - u_e).max(), abs(v_p - v_e).max()):.1e} m/s")   # a trap, not a slip
""", explain=r"""
1. `GFD.ekman_depth` is $\delta=\sqrt{2\nu_v/f}$ *(13.29)*; with `convention="pi"` it returns the oceanographic depth $\pi\delta$.
2. `GFD.ekman_surface(z, tau_x, tau_y, rho, nu_v, f)` returns the two velocity components; it works for either sign of $f$.
3. At 60° N the surface water moves at about 3.5 cm/s, 45° to the right of the wind; at 60° S, 45° to the left.
4. At the Ekman depth the current is 4 % of its surface value and points backwards.
5. `ch13.ekman_surface_printed` is the book's cosine–sine form; for $f>0$ it is identical to ours (it is a *trap*, not a slip:
   true as printed, but only for the northern hemisphere).""")
scratch(r"""
# From scratch: the constant A of D06 and the complex exponential, typed out
A = tau*delta*(1 - 1j)/(2*rho_o*nu_o)                          # D06 step 10: A = tau delta (1 - i)/(2 rho nu_v) [m/s]
V_mine = A*np.exp((1 + 1j)*z_e/delta)                          # D06 step 7: V = A exp((1+i) z/delta)
assert np.allclose(V_mine.real, u_e) and np.allclose(V_mine.imag, v_e)   # same numbers as GFD.ekman_surface
r_x, r_y = GFD.ekman_residual(z_e, u_e, v_e, nu_o, f60)        # what is left of -f v - nu u'' and f u - nu v'' on the grid [m/s^2]
print(f"✓ A e^((1+i)z/δ) equals GFD.ekman_surface; residual of (13.22)–(13.23): {max(abs(r_x).max(), abs(r_y).max())/(f60*abs(A)):.1e} of f|V| (finite-difference error)")   # visible confirmation
""", r"""
Two lines reproduce the library: the constant from the stress condition and one complex exponential. `GFD.ekman_residual`
then puts the profile back into $-fv=\nu_v\frac{d^2u}{dz^2}$ and $fu=\nu_v\frac{d^2v}{dz^2}$ *(13.22)–(13.23)*; what is left is
the truncation error of the second differences.""")
note("N31 [B]", r"""
**The spiral** (our version of the book's figure, `ch13.fig_ekman_surface`): the hodograph with depth marks, and the two
components against depth.""")
fig(r"""
figE = ch13.fig_ekman_surface()                                # hodograph with depth marks + u(z), v(z), for our inputs at 60 N
plt.show()
""",
    see="Left: the tip of the current vector at each depth (the hodograph), starting 45° to the right of the eastward wind and "
        "winding clockwise into the origin. Right: the two velocity components against depth in units of $\\delta$.",
    read="The surface arrow is 45° to the right of the wind (northern hemisphere, $f>0$; mirror for $f<0$); one $\\delta$ down "
         "it has turned one more radian and shrunk by $e$.",
    change="…$\\nu_v$ were four times larger? The spiral would reach twice as deep and the surface current would halve (the "
           "slider figure below).")
nb.plotly(r"""
figSp = go.Figure()                                            # the spiral in three dimensions
z16 = np.linspace(0.0, -4*delta, 17)                           # 17 depths from the surface down to 4 delta [m]
for j, (f_, nm) in enumerate(((f60, "60° N"), (-f60, "60° S"))):   # two hemispheres
    uu_, vv_ = GFD.ekman_surface(z16, tau, 0.0, rho_o, nu_o, f_)   # current at the 17 depths [m/s]
    xs, ys, zs = [], [], []                                    # coordinates of the arrows (as line segments)
    for a_, b_, c_ in zip(uu_, vv_, z16):                      # one arrow per depth, from the axis to the tip
        xs += [0.0, a_*100, None]; ys += [0.0, b_*100, None]; zs += [c_, c_, None]   # in cm/s; None breaks the line
    figSp.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=COLORS["teal"], width=5), name=f"current vectors, {nm}", visible=(j == 0)))   # the arrows
    figSp.add_trace(go.Scatter3d(x=uu_*100, y=vv_*100, z=z16, mode="lines", line=dict(color=COLORS["accent"], width=3), name="spiral through the tips", visible=(j == 0)))   # the curve
figSp.add_trace(go.Scatter3d(x=[0, 3.5], y=[0, 0], z=[0, 0], mode="lines+text", text=["", "wind"], line=dict(color=COLORS["orange"], width=8), name="wind stress (eastward)"))   # the wind
figSp.update_layout(height=520, title="The current turns and fades with depth",
                    updatemenus=[dict(x=0.0, y=1.08, xanchor="left", buttons=[
                        dict(label="60° N (f > 0)", method="update", args=[{"visible": [True, True, False, False, True]}]),
                        dict(label="60° S (f < 0)", method="update", args=[{"visible": [False, False, True, True, True]}])])],
                    scene=dict(xaxis_title="u east [cm/s]", yaxis_title="v north [cm/s]", zaxis_title="z [m]"), margin=dict(l=0, r=0, t=60, b=0))   # axes and the hemisphere dropdown
figSp.show()
""", explain=r"""
1. For each hemisphere the current is computed at 17 depths and drawn as horizontal arrows from the vertical axis.
2. The purple curve joins the arrow tips: the Ekman spiral. The dropdown switches hemisphere.""")
nb.figure_notes(
    see="Horizontal current arrows at 17 depths under an eastward wind (orange), and the spiral through their tips.",
    read="Going down, the arrow shortens and swings clockwise in the northern hemisphere ($f>0$); choose 60° S to see it swing "
         "counter-clockwise.",
    change="…you rotate the picture to look straight down? You see the hodograph of the figure above.")
nb.plotly(r"""
z_f4 = np.linspace(-150.0, 0.0, 76)                            # depths [m]
M_fixed = GFD.ekman_transport(tau, 0.0, rho_o, f60)            # Eq. (13.30): the transport (M_x, M_y) [m^2/s] — no nu_v in it

def f4(log_nu):                                                # hodograph for one eddy viscosity, nu_v = 10**log_nu [m^2/s]
    uu_, vv_ = GFD.ekman_surface(z_f4, tau, 0.0, rho_o, 10**log_nu, f60)   # the spiral for this nu_v
    return {"hodograph (tip of the current at each depth)": (uu_*100, vv_*100),
            "surface current": ([0.0, uu_[-1]*100], [0.0, vv_[-1]*100]),
            f"direction of the transport, τ/(ρf) = {abs(M_fixed[1]):.3f} m²/s — does not change": ([0.0, 0.0], [0.0, -6.0])}

figF4 = slider_figure(f4, "log10 of ν_v [m²/s]", np.round(np.linspace(np.log10(0.005), np.log10(0.3), 25), 3), xlabel="u east [cm/s]",
                      ylabel="v north [cm/s]", title="More mixing: a deeper, slower spiral — the same transport",
                      xrange=(-2, 9), yrange=(-9, 2))          # 25 slider positions, logarithmically spaced
figF4.show()
for nu_ in (0.0075, 0.03, 0.12):                               # three eddy viscosities, each four times the last [m^2/s]
    us_, vs_ = GFD.ekman_surface(0.0, tau, 0.0, rho_o, nu_, f60)   # surface current
    print(f"nu_v = {nu_:6.4f} m^2/s: delta = {GFD.ekman_depth(nu_, f60):5.2f} m, surface speed = {np.hypot(us_, vs_):.5f} m/s, transport = {M_fixed[1]:+.4f} m^2/s")   # delta doubles, speed halves
""", explain=r"""
1. `f4` returns the hodograph for one value of $\nu_v$ (the slider is logarithmic), the surface current as a line from the
   origin, and a fixed arrow in the direction of the total transport.
2. The printed lines quadruple $\nu_v$ twice: $\delta$ doubles each time, the surface speed halves, the transport does not move.""")
nb.figure_notes(
    see="The hodograph for the chosen eddy viscosity. As $\\nu_v$ grows the spiral shrinks toward the origin (slower currents); "
        "the transport arrow, pointing 90° to the right of the wind, stays where it is.",
    read="The eddy viscosity sets how deep the spiral reaches ($\\delta\\propto\\sqrt{\\nu_v}$) and how fast the surface moves "
         "($\\propto1/\\sqrt{\\nu_v}$) — their product, the transport, is independent of it.",
    change="…the latitude were lower? Smaller $f$: a thicker layer and a larger transport, $\\tau/(\\rho f)$.")
note("N33 [B]", r"""
**A pressure gradient only shifts the spiral.** The problem is linear, so a horizontal pressure gradient merely adds a
depth-independent geostrophic velocity to the solution.""")
fig(r"""
u_g0, v_g0 = GFD.ekman_surface(z_e, tau, 0.0, rho_o, nu_o, f60)                # no interior flow
u_g1, v_g1 = GFD.ekman_surface(z_e, tau, 0.0, rho_o, nu_o, f60, U_g=0.05)      # with a geostrophic interior current of 0.05 m/s eastward
fig, ax = plt.subplots(figsize=(5.4, 3.8))
ax.plot(u_g0*100, v_g0*100, color=COLORS["accent"], label="wind only")         # hodograph without interior flow
ax.plot(u_g1*100, v_g1*100, color=COLORS["teal"], label="wind + geostrophic current 5 cm/s")   # the same curve, shifted
ax.plot([0, 5], [0, 0], "o", color=COLORS["muted"])                            # the two deep-water end points
ax.set(xlabel="u east [cm/s]", ylabel="v north [cm/s]", title="The same spiral, shifted by the interior current", aspect="equal")
ax.legend(fontsize=8)
plt.show()
""",
    see="Two hodographs: the wind-driven spiral, and the same spiral with a uniform geostrophic current added.",
    read="The second curve is the first one moved 5 cm/s to the east; it ends at the interior velocity instead of at zero.",
    change="…the interior current were northward? The spiral would shift north instead.")
note("N35 [B]", r"""
**Why the layer does not grow.** With the horizontal vorticity components $\omega_x=-\dfrac{dv}{dz}$, $\omega_y=\dfrac{du}{dz}$,
the balance reads $-f\dfrac{dv}{dz}=\nu_v\dfrac{d^2\omega_y}{dz^2}$ and $-f\dfrac{du}{dz}=\nu_v\dfrac{d^2\omega_x}{dz^2}$ (13.31):
tilting of the planetary vorticity — the $(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u$ term of Chapter 5's
$\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ (5.30)
— balances the downward diffusion of horizontal vorticity. The result is stated; the cell draws the two sides.""")
fig(r"""
vb = GFD.ekman_vorticity_balance(z_e, tau, rho_o, nu_o, f60)   # the two sides of Eq. (13.31) from the closed-form spiral
fig, ax = plt.subplots(figsize=(5.6, 3.8))
ax.plot(vb["tilt_x"], z_e/delta, color=COLORS["teal"], lw=3, label=r"tilting $-f\,dv/dz$")   # tilting of planetary vorticity (Coriolis colour)
ax.plot(vb["diff_x"], z_e/delta, "--", color=COLORS["rose"], label=r"diffusion $\nu_v\,d^2\omega_y/dz^2$")   # diffusion of vorticity (friction colour)
ax.set(xlabel="term of the vorticity balance [s$^{-2}$]", ylabel=r"$z/\delta$", title="Tilting supplies what diffusion removes")
ax.legend(fontsize=8)
plt.show()
print(f"largest difference between the two curves: {abs(vb['tilt_x'] - vb['diff_x']).max():.1e} 1/s^2")   # round-off
""",
    see="The two sides of the first vorticity balance against depth: tilting of planetary vorticity (teal) and diffusion of "
        "horizontal vorticity (rose, dashed).",
    read="The curves coincide at every depth. Without rotation there is no tilting term, nothing can balance the diffusion, and "
         "the layer keeps thickening as in Chapter 8.",
    change="…$f\\to0$? The tilting term vanishes and the steady solution disappears with it ($\\delta\\to\\infty$).")
whatif(r"""
…we do not care about the shape of the spiral, only about how much water moves? Then we should add the arrows up — and the
eddy viscosity drops out (`C05`).""")

# ---------------------------------------------------------------------------------------------------------------------
core("C05", r"Ekman transport: $\displaystyle\int_{-\infty}^{0}u\,dz=0$ and $\displaystyle\int_{-\infty}^{0}v\,dz=-\dfrac{\tau}{\rho f}$ (13.30)",
     "Summed over the whole layer, where does the wind-driven water go — and what happens at a coast?")
problem(r"""
Off Peru and California the wind blows toward the equator along the coast, and the water at the shore is cold and full of
fish. The wind is not pulling the cold water up; it is moving the surface water *offshore*, at right angles to itself, and
deeper water has to rise to replace it. The same sideways transport, driven by the trade winds and the westerlies, piles water
up in the middle of each ocean basin and sets the great gyres turning.""")
idea(words=r"""
Add up all the arrows of the spiral. The along-wind parts cancel; what is left points 90° to the right of the wind (northern
hemisphere, $f>0$; to the left for $f<0$). Reason in one line: over the whole layer the only horizontal forces are the wind
stress on top and the Coriolis force on the total transport — they must cancel, and the Coriolis force is at right angles to
the transport.""")
P("P313", "depth-integrated quantities: transport per unit width [m²/s] and its divergence as a vertical velocity", r"""
Integrating a velocity over depth gives a transport per unit width, in m²/s: the volume crossing a one-metre-wide gate each
second. If more leaves a column sideways than enters, the difference must come through its top or bottom: the divergence of the
transport [m/s] is a vertical velocity.""", code=r"""
z = np.linspace(-100, 0, 1001)            # depth [m]
u = 0.1*np.exp(z/20)                      # a current that decays with depth [m/s]
print("transport =", round(np.trapezoid(u, z), 3), "m^2/s")   # 1.987 (exactly 0.1*20*(1 - e^-5))
""")
P("P312", "first integral: integrating an ODE once across a layer to get a transport without solving it", r"""
If an equation says "something = d(flux)/dz", integrating across the layer gives "total something = flux at the top − flux at
the bottom". You need the flux only at the two ends, not the solution in between.""", code=r"""
tau_top, tau_bottom, rho, f = 0.07, 0.0, 1027.0, 1.0e-4   # stress at the surface and at depth [N/m^2], density, Coriolis parameter
print("y-transport =", round(-(tau_top - tau_bottom)/(rho*f), 3), "m^2/s")   # with no profile in sight
""")
remind("C05")
D("D07", ref="13.30")
note("N32 [B]", r"""
**The transport does not depend on the turbulence model.** The second route of `D07` (steps 6–8) integrates
$-\rho fv=\dfrac{d\tau}{dz}$ once; only "the stress vanishes at depth" is used, so the result holds for any eddy viscosity,
constant or not.""")
nb.worked_example("how much water, and which way", r"""
A stress $\tau=0.1$ N/m² toward the east, $\rho=1000$ kg/m³, $f=10^{-4}$ s⁻¹.

1. Transport $\int v\,dz=-\tau/(\rho f)=-0.1/(1000\times10^{-4})=-1$ m²/s: one cubic metre per second through every metre of
   an east–west line, toward the **south** (to the right of the wind).
2. Along a 1000 km line: $10^{6}$ m³/s = 1 Sverdrup.
3. Southern hemisphere ($f=-10^{-4}$ s⁻¹): +1 m²/s, toward the north — to the left of the wind.
4. Double the eddy viscosity: same answer.""")
code(r"""
Mx, My = GFD.ekman_transport(tau, 0.0, rho_o, f60)             # Eq. (13.30): total transport under an eastward stress at 60 N [m^2/s]
Md = complex(*GFD.ekman_transport_partial(-delta, tau, 0.0, rho_o, nu_o, f60))        # transport carried between z = -delta and the surface
Mp = complex(*GFD.ekman_transport_partial(-np.pi*delta, tau, 0.0, rho_o, nu_o, f60))  # the same down to the Ekman depth pi*delta
print(f"total transport (M_x, M_y) = ({Mx:+.4f}, {My:+.4f}) m^2/s: at right angles to the wind, to its right (60 N)")   # nothing along the wind
print(f"above z = -delta:    {100*abs(Md)/abs(My):.0f} % of the magnitude, direction {np.degrees(np.angle(Md)):.0f} deg (final: -90 deg)")   # not yet in the final direction
print(f"above z = -pi delta: {100*abs(Mp)/abs(My):.0f} % of the magnitude (a small overshoot, removed by the opposed currents below)")   # 104 %
""", explain=r"""
1. `GFD.ekman_transport(tau_x, tau_y, rho, f)` is $\int u\,dz=0$, $\int v\,dz=-\tau/(\rho f)$ *(13.30)* for any wind direction.
2. `GFD.ekman_transport_partial` is the running integral from a depth up to the surface (ours, a closed form): the top
   $\delta$ carries most of the magnitude but not yet in the final direction.
3. Down to $\pi\delta$ the running transport slightly exceeds the total; the weak, backward currents below bring it back.""")
scratch(r"""
# From scratch: add the arrows up numerically, for two very different eddy viscosities
for nu_ in (0.03, 0.12):                                       # nu_v and four times nu_v [m^2/s]
    d_ = GFD.ekman_depth(nu_, f60)                             # the layer thickness for this nu_v [m]
    zz = np.linspace(-40*d_, 0.0, 40001)                       # a grid that reaches far below the layer
    uu_, vv_ = GFD.ekman_surface(zz, tau, 0.0, rho_o, nu_, f60)   # the spiral
    Mx_num, My_num = np.trapezoid(uu_, zz), np.trapezoid(vv_, zz) # the two transports by the trapezoid rule [m^2/s]
    assert np.isclose(My_num, My, rtol=1e-4) and abs(Mx_num) < 1e-6   # same as Eq. (13.30), whatever nu_v is
    print(f"✓ nu_v = {nu_:.2f}: delta = {d_:.2f} m, surface speed = {np.hypot(uu_[-1], vv_[-1]):.5f} m/s, integrated transport = ({Mx_num:+.1e}, {My_num:+.4f}) m^2/s")   # visible confirmation
""", r"""
Numerical integration of the spiral gives the closed-form transport for both eddy viscosities: $\delta$ doubles, the surface
speed halves, the transport stays at $-\tau/(\rho f)$.""")
note("N37 [B]", r"""
**Not in the book — ours** (the book remarks that the constant-viscosity spiral is rarely observed): the layer with a
depth-dependent eddy viscosity, $\dfrac{d}{dz}\Big(K_v\dfrac{dV}{dz}\Big)=if\,(V-V_g)$, solved numerically by `GFD.ekman_solve`.""")
choice(r"""
a second-order finite-volume scheme on a stretched grid with $K_v$ at the cell faces, one complex tridiagonal solve (the
pattern of Chapter 12's eddy-viscosity solver).""")
fig(r"""
z_k = -600.0*np.linspace(1.0, 0.0, 401)**2                     # 401 nodes from -600 m to the surface, crowded near the surface [m]
profiles = {"constant 0.03": lambda zf: 0.03 + 0*zf,           # three eddy-viscosity profiles K_v(z) [m^2/s], given at the cell faces
            "growing with depth 0.03(1 + |z|/20)": lambda zf: 0.03*(1 + abs(zf)/20.0),
            "surface-intensified 0.003 + 0.03 e^(z/30)": lambda zf: 0.003 + 0.03*np.exp(zf/30.0)}
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.9))
M_num = {}                                                     # the three numerical transports
for (nm, Kf), col in zip(profiles.items(), (COLORS["accent"], COLORS["teal"], COLORS["rose"])):   # one curve per profile
    Vk = GFD.ekman_solve(z_k, Kf, f60, tau=tau + 0j, rho=rho_o)   # complex velocity u + iv at the nodes [m/s]
    M_num[nm] = np.trapezoid(Vk, z_k)                          # its depth integral: the transport M_x + i M_y [m^2/s]
    a.plot(Vk.real*100, Vk.imag*100, color=col, label=f"K_v {nm}")   # hodograph
a.set(xlabel="u east [cm/s]", ylabel="v north [cm/s]", title="Three turbulence models, three spirals")
a.legend(fontsize=7)
b.bar(range(3), [m_.imag for m_ in M_num.values()], color=[COLORS["accent"], COLORS["teal"], COLORS["rose"]])   # the three transports
b.axhline(My, color=COLORS["ink"], ls="--")                    # the closed-form value -tau/(rho f)
b.set_xticks(range(3), ["constant", "growing", "surface-\nintensified"])   # the names of the profiles
b.set(ylabel="transport $M_y$ [m$^2$ s$^{-1}$]", title="…one transport")
plt.show()
for nm, m_ in M_num.items():                                   # compare each numerical transport with -tau/(rho f)
    print(f"K_v {nm}: M_y = {m_.imag:+.4f} m^2/s ({100*(m_.imag/My - 1):+.2f} % from -tau/(rho f)), M_x = {m_.real:+.4f}")   # the measured agreement
V_const = GFD.ekman_solve(z_k, profiles["constant 0.03"], f60, tau=tau + 0j, rho=rho_o)   # the constant case once more
u_ex, v_ex = GFD.ekman_surface(z_k, tau, 0.0, rho_o, 0.03, f60)   # its closed form
print(f"constant K_v: numerical vs closed-form spiral, max difference = {abs(V_const - (u_ex + 1j*v_ex)).max()/np.hypot(u_ex[-1], v_ex[-1]):.1e} of the surface speed")   # scheme accuracy
""",
    see="Left: hodographs for three eddy-viscosity profiles — very different shapes and surface angles. Right: their depth-"
        "integrated transports, with the dashed line at $-\\tau/(\\rho f)$.",
    read="The shape of the spiral belongs to the turbulence model; the transport belongs to Newton's law. The printed lines "
         "give the measured agreement of each numerical transport with the closed form.",
    change="…the water were only a few $\\delta$ deep? Then the stress would not vanish at the bottom and the transport would "
           "change (note N34 in `C06`).")
D("D09")
note("N36 [B]", r"""
**Not in the book — ours** (the book gives only the consequence, rising air in a low): Ekman pumping
$w_E=\dfrac1\rho\,\mathbf e_z\cdot\nabla\times\Big(\dfrac{\boldsymbol\tau}{f}\Big)$ under a surface layer, and
$w=\dfrac\delta2\zeta_g$ above a bottom layer (`D09`). One line, labelled beyond the book: **Sverdrup balance** — in the
interior below the surface layer, stretching by this pumping is balanced by moving columns north or south,
$\beta V=f\,w_E$ per unit width (depth-integrated, with no vertical velocity at depth): the wind's curl sets the north–south
transport of a whole gyre.""")
fig(r"""
xg2 = np.linspace(0.0, 4.0e6, 81); yg2 = np.linspace(-1.0e6, 1.0e6, 81)   # a band 4000 km long and 2000 km wide [m]
X2, Y2 = np.meshgrid(xg2, yg2)                                 # grids indexed [j, i] = (y, x)
tau_x = -0.07*np.cos(np.pi*(Y2 + 1.0e6)/2.0e6)                 # easterlies in the south, westerlies in the north [N/m^2]
tau_y = np.zeros_like(tau_x)                                   # no north-south stress
wE = GFD.ekman_pumping(tau_x, tau_y, xg2[1] - xg2[0], yg2[1] - yg2[0], rho_o, f35)   # w_E = (1/rho) curl(tau/f) at 35 N [m/s]
fig, ax = plt.subplots(figsize=(7.6, 3.8))
pm = ax.pcolormesh(xg2/1e3, yg2/1e3, wE*86400*365, cmap="RdBu", shading="auto", vmin=-45, vmax=45)   # pumping in metres per year
ax.quiver(xg2[::8]/1e3, yg2[::8]/1e3, tau_x[::8, ::8], tau_y[::8, ::8], color=COLORS["orange"])   # the wind stress
fig.colorbar(pm, label="Ekman pumping $w_E$ [m per year]")
ax.set(xlabel="x east [km]", ylabel="y north [km]", title="Easterlies south, westerlies north: water is pumped down (35° N)")
plt.show()
cu = ch13.coastal_upwelling(-0.07, "east", lat)                # an equatorward wind along an eastern boundary at 35 N
print(f"open ocean: strongest pumping = {wE[1:-1, 1:-1].min():.2e} m/s = {wE[1:-1, 1:-1].min()*86400*365:.0f} m per year (downward)")   # compare pi tau0/(Ly rho f)
print(f"coast: offshore transport = {cu['transport_offshore']:.3f} m^2/s, upwelling: {cu['upwelling']}; fed over a 40 km strip -> w = {cu['transport_offshore']/4.0e4:.1e} m/s = {cu['transport_offshore']/4.0e4*86400:.1f} m per day")   # coastal upwelling
""",
    see="A band of wind stress (orange arrows) that blows westward in the south and eastward in the north — the pattern over a "
        "subtropical gyre — and the vertical velocity it forces at the base of the Ekman layer (colour).",
    read="The Ekman transports under the two wind belts point toward each other (each to the right of its wind), the water "
         "converges and is pumped **down** everywhere in the band — tens of metres per year. At a coast the same mechanism is "
         "far stronger: metres per day (second printed line).",
    change="…the same wind band lay in the southern hemisphere? The transports would point apart and the water would be pumped "
           "**up**.")
note("N50 [C]", r"""
**Named, not presented.** The wind-driven circulation of the gyres and Stommel's explanation of strong western boundary
currents are mentioned by the book only in passing. Pointer: the Sverdrup line above; any dynamical-oceanography text for the
gyre.""")
nb.live(r"""
def ekman_live(tau=0.07, nu_v=0.03, lat_deg=60.0):             # free wind stress [N/m^2], eddy viscosity [m^2/s], latitude [degrees, either sign]
    if lat_deg == 0:                                           # on the equator f = 0 and there is no steady Ekman layer
        print("f = 0 on the equator: no steady Ekman layer"); return
    f_ = GFD.coriolis_parameter(np.deg2rad(lat_deg))           # Coriolis parameter at this latitude [1/s]
    d_ = GFD.ekman_depth(nu_v, f_)                             # layer thickness [m]
    zz = np.linspace(-5*d_, 0.0, 200)                          # depths [m]
    uu_, vv_ = GFD.ekman_surface(zz, tau, 0.0, 1027.0, nu_v, f_)   # the spiral
    Mx_, My_ = GFD.ekman_transport(tau, 0.0, 1027.0, f_)       # Eq. (13.30)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.plot(uu_*100, vv_*100, color=COLORS["accent"])          # hodograph
    ax.annotate("", xy=(0, np.sign(My_)*abs(uu_[-1])*100), xytext=(0, 0), arrowprops=dict(color=COLORS["teal"], width=2))   # direction of the transport
    ax.set(xlabel="u east [cm/s]", ylabel="v north [cm/s]", aspect="equal", title=f"δ = {d_:.1f} m, transport = {My_:+.3f} m²/s")
    plt.show()

_ = live(ekman_live, tau=(0.01, 0.3, 0.01), nu_v=(0.005, 0.3, 0.005), lat_deg=(-80, 80, 5))   # three sliders
""", explain=r"""
`live(fn, name=(min, max, step), …)` makes one slider per argument and redraws when you release it. Here: the hodograph and the
direction of the transport (teal arrow) for a free wind stress, eddy viscosity and latitude of either sign. On the web page the
slider figure of `C04` carries the same idea.""")
explainer("ekman_spiral", "The wind blows east — why does the water go south?",
          "orbiting the spiral in 3-D, sliding a depth cursor down the hodograph and watching the running transport swing round "
          "to 90° connects three pictures of one solution; changing $\\nu_v$ changes $\\delta$ and the surface speed and leaves "
          "the transport untouched.",
          ["Press play: the depth cursor descends and the running transport arrow swings to 90° from the wind.",
           "Use the preset 'four times ν_v': twice as deep, half as fast at the surface, same transport.",
           "Switch to S: surface current and transport flip to the left of the wind.",
           "Choose 'wind along a coast' and read the status: upwelling or downwelling?"])
whatif(r"""
…the boundary is solid and the fluid above it is already moving? The same equation with a different boundary condition: the
bottom Ekman layer (`C06`).""")


# =====================================================================================================================
# A.7  §13.7 — C06 the bottom Ekman layer
# =====================================================================================================================
nb.section("13.7", "Ekman Layer on a Rigid Surface", intro=r"""
**What is this section about?** The friction layer under a geostrophic flow — the lowest kilometre of the atmosphere, or the
bottom of the ocean. ⚠️ Here $z=0$ is the solid surface and the fluid is $z>0$.""")
core("C06", r"The bottom Ekman layer: $u=U\big[1-e^{-z/\delta}\cos(z/\delta)\big]$, $v=Ue^{-z/\delta}\sin(z/\delta)$ (13.41)",
     "If the wind follows the isobars, why does air spiral into a low?")
problem(r"""
Satellite pictures show clouds spiralling *into* a depression, not circling it. Aloft the wind does follow the isobars. But
near the ground friction slows the air; slower air feels a weaker Coriolis force; the pressure force, which does not care about
speed, is no longer cancelled, and the air drifts toward low pressure. Air converging on a low has nowhere to go but up —
clouds and rain. Stir a cup of tea and the leaves collect in the middle for the same reason (Ch. 9).""")
idea(r"""
   aloft (no friction)                 near the ground (friction)
        LOW                                  LOW
         ↑ pressure force                     ↑ pressure force (unchanged)
         ●━━━▶ wind along the isobar           ●━━▶ slower wind, turned toward LOW
         ↓ Coriolis (equal and opposite)      ↘ Coriolis (shorter, turned with the wind)
                                              ← friction (roughly against the wind) closes the triangle
""", r"""
Three forces must add to zero. Aloft two suffice; near the ground the weakened Coriolis force cannot cancel the pressure force
alone, and friction supplies the rest — which is only possible if the wind has a component toward low pressure.""")
note("N38 [C]", r"""
**The interior is geostrophic:** $fU=-\dfrac1\rho\dfrac{dp}{dy}$ (13.32). It enters the layer problem as the constant $U$.""")
note("N39 [B]", r"""
**Three forces in the layer:** $-fv=\nu_v\dfrac{d^2u}{dz^2}$ and $fu=\nu_v\dfrac{d^2v}{dz^2}+fU$ (13.33)–(13.34) — the pressure
gradient is hidden in the term $fU$.""")
note("N40 [C]", r"""
**Far-field condition:** $u=U,\ v=0$ as $z\to\infty$ (13.35).""")
note("N41 [C]", r"""
**No slip:** $u=0,\ v=0$ at $z=0$ (13.36).""")
remind("C06")
D("D08", ref="13.41")
note("N42 [B]", r"""
**Complex form** (step 2 of `D08`): $\dfrac{d^2V}{dz^2}=\dfrac{if}{\nu_v}(V-U)$ (13.37).""")
note("N43 [C]", r"""
**Its far-field condition:** $V=U$ as $z\to\infty$ (13.38).""")
note("N44 [C]", r"""
**Its no-slip condition:** $V=0$ at $z=0$ (13.39).""")
note("N45 [B]", r"""
**General solution** (step 3 of `D08`): $V=A\,e^{-(1+i)z/\delta}+B\,e^{(1+i)z/\delta}+U$ (13.40) — the homogeneous part plus the
particular solution $V=U$.""")
note("N47 [B]", r"""
**Transport toward low pressure** (steps 8–9 of `D08`):
$\displaystyle\int_0^\infty v\,dz=U\Big[\dfrac{\nu_v}{2f}\Big]^{1/2}=\dfrac12U\delta$.""")
nb.worked_example("the layer under a 10 m/s wind", r"""
Take $U=10$ m/s, $\nu_v=5$ m²/s, $f=10^{-4}$ s⁻¹.

1. $\delta=\sqrt{2\nu_v/f}=\sqrt{10/10^{-4}}=\sqrt{10^{5}}\approx316$ m.
2. At $z=\delta$: $e^{-1}=0.368$, $\cos1=0.540$, $\sin1=0.841$ → $u=10(1-0.368\times0.540)=8.0$ m/s,
   $v=10\times0.368\times0.841=3.1$ m/s: the wind points 21° to the left of the isobars' direction (toward low pressure;
   northern hemisphere).
3. Very near the ground $u\approx v$: 45°.
4. Transport toward low pressure: $\tfrac12U\delta=0.5\times10\times316\approx1580$ m²/s.""")
code(r"""
U_g, nu_a = inp["U_g"], inp["nu_v_atm"]                        # geostrophic wind 12 m/s aloft and nu_v = 7 m^2/s (our inputs), at 60 N
delta_b = GFD.ekman_depth(nu_a, f60)                           # Eq. (13.29) again: layer thickness [m]
for s_, nm in ((0.5, "0.5 δ"), (1.0, "δ"), (np.pi/2, "πδ/2"), (np.pi, "πδ")):   # four heights in units of delta
    ub, vb_ = GFD.ekman_bottom(s_*delta_b, U_g, 0.0, nu_a, f60)   # Eq. (13.41): (u, v) at this height [m/s]
    print(f"z = {nm:>5}: (u, v) = ({ub:5.2f}, {vb_:4.2f}) m/s, {np.degrees(np.arctan2(vb_, ub)):4.1f} deg from the isobars")   # the wind turns toward the isobars with height
Mbx, Mby = GFD.ekman_bottom_transport(U_g, 0.0, nu_a, f60)     # transport of the departure from U, integrated over the layer [m^2/s]
fb = GFD.ekman_force_balance(0.5*delta_b, U_g, nu_a, f60)      # the three forces per unit mass at z = 0.5 delta [m/s^2]
print(f"delta = {delta_b:.1f} m (pi*delta = {np.pi*delta_b:.0f} m); transport (M_x, M_y) = ({Mbx:+.0f}, {Mby:+.0f}) m^2/s; (1/2) U delta = {0.5*U_g*delta_b:.0f}")   # thickness and transport
print(f"sum of the three forces at 0.5 delta: ({fb['sum'][0]:.1e}, {fb['sum'][1]:.1e}) m/s^2 (round-off)")   # the triangle closes
""", explain=r"""
1. `GFD.ekman_bottom(z, U_g, V_g, nu_v, f)` is $u=U[1-e^{-z/\delta}\cos(z/\delta)]$, $v=Ue^{-z/\delta}\sin(z/\delta)$ *(13.41)*.
2. Low down the wind blows about 30° across the isobars toward low pressure; the angle shrinks with height and is zero at
   $z=\pi\delta$.
3. `GFD.ekman_bottom_transport` integrates the departure from $U$: $+\tfrac12U\delta$ toward low pressure (north here, for
   $f>0$) and $-\tfrac12U\delta$ along the isobars (the layer carries less than the interior flow would).
4. `GFD.ekman_force_balance` returns Coriolis, pressure and friction as vectors; their sum is zero to round-off.""")
scratch(r"""
# From scratch: the three forces at z = 0.5 delta, with the friction from second differences
h_ = 1.0e-3*delta_b                                            # a small height step [m]
zz3 = 0.5*delta_b + np.array([-h_, 0.0, h_])                   # three neighbouring heights
u3, v3 = GFD.ekman_bottom(zz3, U_g, 0.0, nu_a, f60)            # the profile there
d2u = (u3[0] - 2*u3[1] + u3[2])/h_**2                          # second derivative of u by differences [1/(m s)]
d2v = (v3[0] - 2*v3[1] + v3[2])/h_**2                          # second derivative of v
cor = f60*np.array([v3[1], -u3[1]])                            # Coriolis force per unit mass: f (v, -u)
pgf = f60*np.array([0.0, U_g])                                 # pressure force per unit mass: (0, f U), toward low pressure
fric = nu_a*np.array([d2u, d2v])                               # friction per unit mass: nu_v (u'', v'')
assert np.allclose(cor + pgf + fric, 0, atol=1e-7)             # Eqs. (13.33)-(13.34): the three add to zero
assert np.allclose(cor, fb["coriolis"]) and np.allclose(pgf, fb["pressure"])   # the same vectors as the library returns
print(f"✓ Coriolis {np.round(cor*1e4, 2)}, pressure {np.round(pgf*1e4, 2)}, friction {np.round(fric*1e4, 2)} [1e-4 m/s^2] add to zero")   # visible confirmation
""", r"""
The balance $-fv=\nu_v\frac{d^2u}{dz^2}$, $fu=\nu_v\frac{d^2v}{dz^2}+fU$ *(13.33)–(13.34)* rewritten as "Coriolis + pressure +
friction = 0", with the friction taken from second differences of the profile. The three vectors close to within the
finite-difference error.""")
note("N46 [B]", r"""
**The profile** (our version of the book's figure, `ch13.fig_ekman_bottom`). **Stated precisely:** the component $u$ overshoots
$U$. Its largest value is at $z=3\pi\delta/4$, where $u/U=1+e^{-3\pi/4}/\sqrt2=1.067$ (6.7 % above $U$); at $z=\pi\delta$, where
$v$ first returns to zero, $u/U=1+e^{-\pi}=1.043$.""")
fig(r"""
figB = ch13.fig_ekman_bottom()                                 # hodograph with height marks + u(z), v(z) for our inputs at 60 N
plt.show()
zz_b = np.linspace(0.0, 6*delta_b, 6001)                       # a fine grid of heights [m]
ub_all, _ = GFD.ekman_bottom(zz_b, U_g, 0.0, nu_a, f60)        # the along-isobar component
k_max = np.argmax(ub_all)                                      # where u is largest
print(f"largest u = {ub_all[k_max]:.2f} m/s = {ub_all[k_max]/U_g:.3f} U at z/delta = {zz_b[k_max]/delta_b:.3f} (3 pi/4 = {3*np.pi/4:.3f})")   # the overshoot
print(f"at z = pi*delta: u = {GFD.ekman_bottom(np.pi*delta_b, U_g, 0.0, nu_a, f60)[0]:.2f} m/s = {1 + np.exp(-np.pi):.3f} U")   # where v returns to zero
""",
    see="Left: the hodograph from the origin (no slip at the ground) to the point $(U, 0)$ (the geostrophic wind aloft), with "
        "marks at $z/\\delta=\\pi/4$, $\\pi/2$, $\\pi$. Right: $u$ and $v$ against height.",
    read="The hodograph leaves the origin at 45° to the left of $U$ (northern hemisphere, $f>0$; mirror for $f<0$) and winds "
         "into the point $(U,0)$, passing slightly beyond it: $u$ peaks at 1.067 $U$ at $z=3\\pi\\delta/4$ (printed line).",
    change="…the surface were rougher (a larger eddy viscosity)? The thickness $\\delta$ grows; in units of $\\delta$ the picture is unchanged.")
note("N48 [B]", r"""
**Eddies, not molecules, carry the stress.** With the molecular viscosity of air the layer would be half a metre thick; the
observed layer is of the order of a kilometre. Turning the formula round gives the eddy viscosity that such a depth implies.""")
code(r"""
print(f"laminar air (nu = 1.5e-5 m^2/s) at 60 N: delta = {GFD.ekman_depth(1.5e-5, f60):.2f} m")       # absurdly thin
print(f"eddy viscosity implied by delta = 1000 m at 60 N: {GFD.eddy_viscosity_from_depth(1000.0, f60):.0f} m^2/s")   # nu_v = |f| delta^2 / 2
""", explain=r"""
**What does this show?** A kilometre-deep friction layer needs an eddy viscosity some four million times the molecular one.""")
loose(10, r"a laminar thickness for air about a quarter below what the formula gives for the inputs stated beside it",
      r"the value of $\delta=\sqrt{2\nu_v/f}$ (13.29) — our numbers above are computed, not quoted.")
note("N49 [B]", r"""
**The force triangle at three heights:**
$-f\,\mathbf e_z\times\mathbf u-\dfrac1\rho\nabla p+\nu_v\dfrac{d^2\mathbf u}{dz^2}=0$. Inflow and rising air in a low, outflow and
sinking in a high — in both hemispheres. The cell also prints a *spin-down time* (ours, one line, not in the book): the pumping
$w=\tfrac\delta2\zeta_g$ out of the layer squashes the column of height $H$ above it, and by the stretching term of `C13` that
changes its vorticity at the rate $d\zeta_g/dt=-f\,w/H=-\big(f\delta/2H\big)\zeta_g$: the vortex decays by a factor $e$ in the time
$2H/(f\delta)$.""")
fig(r"""
fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.9), sharex=True, sharey=True)
for ax, s_ in zip(axes, (0.1, 1.0, 3.0)):                      # three heights in units of delta
    fb3 = GFD.ekman_force_balance(s_*delta_b, U_g, nu_a, f60)  # the three force vectors there [m/s^2]
    start = np.zeros(2)                                        # draw the vectors head to tail, starting at the origin
    for key, col in (("pressure", COLORS["orange"]), ("coriolis", COLORS["teal"]), ("friction", COLORS["rose"])):   # house colours
        vec = np.array(fb3[key])*1e3                           # this force in units of 1e-3 m/s^2
        ax.annotate("", xy=start + vec, xytext=start, arrowprops=dict(color=col, width=2.5, headwidth=8))   # one arrow
        start = start + vec                                    # the next arrow starts where this one ends
    ax.annotate("", xy=(fb3["u"]/U_g, fb3["v"]/U_g), xytext=(0, 0), arrowprops=dict(color=COLORS["muted"], width=1, headwidth=6))   # the wind, as a fraction of U (grey)
    ax.set(title=f"z = {s_} δ: wind {np.degrees(fb3['angle_to_isobars']):.0f}° off the isobars", xlim=(-1.3, 1.8), ylim=(-0.4, 1.9), aspect="equal", xlabel="x (along the isobars)")
axes[0].set_ylabel("y (toward LOW pressure)")
for ax in axes:
    ax.title.set_fontsize(9)                                   # small titles, so that the three do not collide
fig.suptitle("orange pressure + teal Coriolis + rose friction = 0 at every height (grey: the wind, in units of U)", fontsize=10)
plt.show()
w_n = GFD.ekman_pumping_bottom(1.0e-5, nu_a, f60)              # w = (delta/2) zeta_g above a northern low with zeta_g = 1e-5 1/s
w_s = GFD.ekman_pumping_bottom(-1.0e-5, nu_a, -f60)            # a southern low: f < 0 and zeta_g < 0
print(f"pumping out of the layer under a low: north {w_n*1e3:+.2f} mm/s, south {w_s*1e3:+.2f} mm/s (upward in both)")   # rising air in lows
print(f"spin-down time of a 9 km deep vortex, 2H/(f delta) = {2*9000.0/(f60*delta_b)/86400:.2f} days")   # how fast friction kills a low
""",
    see="At three heights the pressure force (orange), the Coriolis force (teal) and friction (rose) drawn head to tail: each "
        "triangle closes. The grey arrow is the wind.",
    read="Near the ground the pressure force is balanced mostly by friction and the wind blows 40° across the isobars; at "
         "$3\\delta$ friction has almost vanished and the triangle has collapsed to the geostrophic pair. The inflow toward "
         "the low must rise: the printed pumping is upward for a low in *either* hemisphere.",
    change="…the centre were a high? All horizontal vectors reverse; the outflow near the ground is fed by sinking air — clear "
           "skies.")
note("N34 [B]", r"""
**Shallow water: two layers in one.** Observed current profiles in shallow seas look like a surface Ekman layer on top of a
bottom one. The observation is not reproduced; **our** exact finite-depth solution (`GFD.ekman_finite_depth`: wind stress at
the surface, no slip at the bottom, a geostrophic interior current) shows the same structure.""")
fig(r"""
H_sh = 100.0                                                   # water depth [m], about 4.6 delta for the ocean inputs
z_sh = np.linspace(-H_sh, 0.0, 301)                            # from the sea bed to the surface
V_sh = GFD.ekman_finite_depth(z_sh, tau, 0.0, 0.05, 0.0, H_sh, nu_o, rho_o, f60)   # stress 0.07 N/m^2 on top, interior current 0.05 m/s, no slip below
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.8))
a.plot(V_sh.real*100, z_sh, color=COLORS["accent"], label="u (east)")     # eastward component against depth
a.plot(V_sh.imag*100, z_sh, color=COLORS["teal"], label="v (north)")      # northward component
a.set(xlabel="velocity [cm/s]", ylabel="z [m] (surface at 0, bed at −100)", title="A surface layer and a bottom layer")
a.legend(fontsize=8)
b.plot(V_sh.real*100, V_sh.imag*100, color=COLORS["accent"])   # the hodograph
b.plot([5], [0], "o", color=COLORS["muted"])                   # the interior (geostrophic) velocity
b.set(xlabel="u [cm/s]", ylabel="v [cm/s]", title="Hodograph: bed → interior → surface", aspect="equal")
plt.show()
""",
    see="Left: the two velocity components through 100 m of water. Right: the hodograph, from zero at the bed, through the "
        "interior velocity (grey dot), to the wind-driven surface current.",
    read="Near the bed the profile is the bottom spiral of this block; near the surface it is the wind-driven spiral of `C04`; "
         "in between the current is the geostrophic interior value.",
    change="…the depth were less than about $2\\delta$? The two layers would merge and the current would run nearly downwind.")
nb.plotly(r"""
s_f5 = np.linspace(0.0, 6.0, 61)                               # height in units of delta

def f5(U_):                                                    # hodograph for one geostrophic wind speed U [m/s] (nu_v = 7 m^2/s, 60 N)
    uu_, vv_ = GFD.ekman_bottom(s_f5*delta_b, U_, 0.0, nu_a, f60)   # Eq. (13.41)
    u1, v1 = GFD.ekman_bottom(delta_b, U_, 0.0, nu_a, f60)     # the wind at z = delta
    return {"hodograph: ground (origin) → aloft (U, 0)": (uu_, vv_),
            "wind at z = δ (always 21° from the isobars)": ([0.0, u1], [0.0, v1]),
            "geostrophic wind aloft": ([0.0, U_], [0.0, 0.0])}

figF5 = slider_figure(f5, "geostrophic wind U", np.arange(2.0, 30.1, 2.0), unit="m/s", xlabel="u along the isobars [m/s]",
                      ylabel="v toward low pressure [m/s]", title="A stronger wind stretches the spiral; no angle changes",
                      xrange=(0, 33), yrange=(-2, 11))         # 15 slider positions
figF5.show()
for nu_ in (2.0, 7.0, 20.0):                                   # three eddy viscosities [m^2/s]
    d_ = GFD.ekman_depth(nu_, f60)                             # thickness for each
    print(f"nu_v = {nu_:4.1f} m^2/s: delta = {d_:5.1f} m, cross-isobar transport (1/2) U delta = {0.5*U_g*d_:6.0f} m^2/s for U = {U_g:.0f} m/s")   # thickness and transport grow with nu_v
""", explain=r"""
1. `f5` returns the hodograph for one wind speed, the wind vector at $z=\delta$ and the geostrophic wind.
2. The printed lines show what the eddy viscosity changes: the thickness and the transport toward low pressure — not the shape.""")
nb.figure_notes(
    see="The hodograph of the bottom layer for the chosen geostrophic wind, with the wind at $z=\\delta$ and the wind aloft.",
    read="The whole figure scales with $U$: the cross-isobar angle at a given $z/\\delta$ does not depend on the wind speed at "
         "all (nor on $\\nu_v$).",
    change="…$\\nu_v$ doubled? In this picture nothing — only the height at which each point of the spiral sits, and the "
           "transport (printed lines).")
explainer("ekman_force_balance", "If the wind follows the isobars, why does air spiral into a low?",
          "sliding the height cursor from the ground to the top of the layer turns the force triangle continuously from "
          "'pressure against friction' to 'pressure against Coriolis'; the triangle always closes.",
          ["Slide the height to 0: the wind is 45° from the isobars.",
           "Jump to the preset z = πδ: v is zero and u is 4 % above U.",
           "Switch the centre from low to high and read the pumping sign.",
           "Switch to S and check that air still flows *into* the low."])
whatif(r"""
…we forget friction again and ask how such a layer of fluid *moves in time* — its waves? For that we need equations for a
thin layer with a free surface (`C07`).""")

# =====================================================================================================================
# A.8  §13.8 — C07 the shallow-water equations
# =====================================================================================================================
nb.section("13.8", "Shallow-Water Equations", intro=r"""
**What is this section about?** The simplest model that still has gravity waves and rotation: one thin layer of uniform
density with a free surface. Three equations for the surface height and the two horizontal velocities. Sections 13.10–13.15
are all about its solutions. ⚠️ Here $z=0$ is the flat bottom and the surface is at $z=H+\eta$.""")
core("C07", r"The linear shallow-water equations: $\dfrac{\partial\eta}{\partial t}+H\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$, $\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x}$, $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.45)",
     "What is the least we must keep to describe a tide or a tsunami crossing an ocean on a turning earth?")
problem(r"""
A tsunami in mid-ocean is a few hundred kilometres long in water four kilometres deep. To such a wave the ocean is a puddle:
the water moves almost horizontally, the same at every depth, and the pressure at any point is just the weight of water above
it. Then the whole three-dimensional problem collapses to a map: how high is the surface, and which way is the column moving.
Tides, storm surges, the adjustment of the ocean to a change of wind — and, with one substitution, each internal mode of a
stratified ocean — all obey these three equations.""")
idea(r"""
   columns converge → the surface rises → the slope pushes the columns apart again → …
        →  ←                 ▲                        ←  →
""", r"""
A column of water is the unit. Its height changes when its neighbours crowd in; the resulting slope accelerates the columns.""")
note("N54 [C]", r"""
**The geometry** (our sketch): mean depth $H$, surface displacement $\eta$, height $z$ measured from the flat bottom. Which
$z=0$ each section uses is in the table of the conventions block (trap T4).""")
fig(r"""
figL = draw_shallow_layer()                                    # our sketch of the layer: H, eta, and z measured from the bottom
plt.show()
""",
    see="A layer of mean depth $H$ whose surface is displaced by $\\eta$; the horizontal velocity is the same at every height "
        "of a column.",
    read="Everything in this section is a function of $(x, y, t)$ only: the vertical coordinate has been integrated away.",
    change="…the bottom were uneven? The depth $H$ becomes $H(x,y)$ and the total depth $h=H+\\eta$ takes its place (`C13`).")
nb.recap("R17", "Hydrostatic pressure under a displaced surface", r"""
With zero pressure at the surface $z=H+\eta$, the pressure at height $z$ is the weight of the water above:
$p=\rho g(H+\eta-z)$.""", where=r"Ch. 7, long waves are hydrostatic: $p'=\rho g\eta$ (7.52)")
nb.current_core = "C07"
nb.recap("R18", "The pressure gradient is the surface slope", r"""
Differentiating, $\dfrac{\partial p}{\partial x}=\rho g\dfrac{\partial\eta}{\partial x}$ and $\dfrac{\partial p}{\partial y}=\rho g\dfrac{\partial\eta}{\partial y}$ (13.42):
neither depends on $z$, so every level of a column is pushed alike and a flow that starts depth-independent stays so.""",
         where="Ch. 7 §7.2")
nb.current_core = "C07"
remind("C07")
nb.md(r"""
Two more reminders: the transport per unit width and its divergence were primed in `C05` (P313); and Chapter 7's kinematic
surface condition, $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}\cong\frac{\partial\eta}{\partial t}$ (7.17), is used here in
its exact form $w(\eta)=D\eta/Dt$: a parcel on the surface stays on the surface.""")
D("D10", ref="13.45")
note("N51 [B]", r"""
**Continuity integrated over a column** (step 3 of `D10`):
$(H+\eta)\dfrac{\partial u}{\partial x}+(H+\eta)\dfrac{\partial v}{\partial y}+w(\eta)-w(0)=0$ (13.43).""")
note("N52 [B]", r"""
**The nonlinear continuity equation** (steps 4–6 of `D10`):
$\dfrac{\partial\eta}{\partial t}+\dfrac{\partial}{\partial x}\big[u(H+\eta)\big]+\dfrac{\partial}{\partial y}\big[v(H+\eta)\big]=0$ (13.44)
— the divergence of the transport lowers the surface.""")
nb.worked_example("how fast, and how strong a current", r"""
Take $H=4000$ m, $g\approx10$ m/s², no rotation.

1. Try $\eta=F(x-ct)$: the two one-dimensional equations give $c^2=gH$, so $c=\sqrt{40\,000}=200$ m/s (720 km/h — a jet
   aircraft).
2. An ocean 6000 km wide is crossed in 30 000 s ≈ 8 h.
3. A crest 1 m high carries a current $u=c\eta/H=200\times1/4000=0.05$ m/s: the water barely moves; the *shape* travels.
4. With $H=1.3$ m instead (the equivalent depth of an internal mode, `C08`): $c=\sqrt{13}\approx3.6$ m/s.""")
P("P335", "a C-grid shallow-water step and its CFL limit with c = √(gH) (our choice of scheme)", r"""
To step the equations in time we store $\eta$ at the centres of grid cells and the velocities on the cell faces (a "C-grid"),
so every difference is taken across exactly one cell. Forward–backward stepping: update $\eta$ with the old velocities, then
the velocities with the new $\eta$. It is stable only if a wave cannot cross more than one cell per step:
$c\,\Delta t/\Delta x\le1$ with $c=\sqrt{gH}$ (the CFL limit of Ch. 10). The ratio $c\,\Delta t/\Delta x$ is called the *Courant number*. This scheme is OUR choice — the book has no numerical
model.""", code=r"""
dx, H = 1.0e4, 1.33                    # 10 km cells; a layer 1.33 m deep
c = np.sqrt(G0*H)                      # long-wave speed: 3.61 m/s
print(f"c = {c:.2f} m/s; time step at Courant number 0.5: {0.5*dx/c:.0f} s")   # 3.61 m/s and 1385 s
""")
choice(r"""
forward–backward time stepping on a staggered grid in one dimension here (`SW.linear_1d_run`); an energy-conserving C-grid
scheme with third-order Runge–Kutta steps for the two-dimensional model `SW.ShallowWater`, used only to make the cached runs of
`C11` and `C13`.""")
code(r"""
He1 = 1.3286                                                   # layer depth [m]: the equivalent depth of our first internal mode (C08)
dx1, dt1, n1 = 1.0e4, 1385.0, 400                              # 10 km cells, 1385 s steps (Courant number 0.5), 400 cells
x1 = (np.arange(n1) + 0.5)*dx1 - 2.0e6                         # cell centres from -2000 km to +2000 km [m]
eta_bump = np.exp(-(x1/1.0e5)**2)                              # a bump 1 m high and 100 km wide
out0 = SW.linear_1d_run(eta_bump, dx=dx1, dt=dt1, n_steps=144, H=He1, f=0.0)   # march 144 steps (2e5 s) without rotation
eta_end = out0["eta"][-1]                                      # the surface at the last step
right = x1 > 0                                                 # the right half of the line
x_peak = x1[right][np.argmax(eta_end[right])]                  # where the right-going pulse is
print(f"c = sqrt(gH) = {GFD.long_wave_speed(He1):.3f} m/s for H = {He1} m;  {GFD.long_wave_speed(4200.0):.2f} m/s for H = 4200 m")   # Eq. (13.86)
print(f"after {out0['t'][-1]:.0f} s: two pulses of height {eta_end[right].max():.3f} m at x = ±{x_peak/1e3:.0f} km -> speed {x_peak/out0['t'][-1]:.2f} m/s")   # the bump has split
""", explain=r"""
1. `SW.linear_1d_run` marches the linear equations with no variation in $y$ from an initial surface shape at rest.
2. Without rotation the bump splits into two pulses of half the height, running apart at the long-wave speed — Chapter 7's
   result, recovered by our scheme (the small speed deficit is the dispersion of a 10 km grid).
3. `GFD.long_wave_speed(H)` is $c=\sqrt{gH}$ *(13.86)*: 3.6 m/s for a layer 1.33 m deep, about 200 m/s for the real ocean depth.
4. The time step obeys the CFL limit of the primer above.""")
scratch(r"""
# From scratch: the forward-backward march, every line of the scheme written out
def march(eta0, f_, n_steps):                                  # returns the surface after n_steps
    eta = eta0.copy()                                          # surface height at the cell centres [m]
    u = np.zeros(len(eta0) + 1)                                # x-velocity at the cell faces [m/s]; the two end faces are walls
    v = np.zeros(len(eta0))                                    # y-velocity at the cell centres [m/s]
    for _ in range(n_steps):                                   # one time step per pass
        u_c = 0.5*(u[1:] + u[:-1])                             # old u averaged to the centres
        eta = eta - dt1*He1*np.diff(u)/dx1                     # Eq. (13.45), first member: d(eta)/dt = -H du/dx (old u)
        v = v - dt1*f_*u_c                                     # third member with d/dy = 0: dv/dt = -f u (old u)
        v_f = 0.5*(v[1:] + v[:-1])                             # new v averaged to the interior faces
        u[1:-1] = u[1:-1] + dt1*(f_*v_f - G0*np.diff(eta)/dx1) # second member: du/dt = f v - g d(eta)/dx (new eta, new v)
    return eta

assert np.allclose(march(eta_bump, 0.0, 144), out0["eta"][-1], atol=1e-12)   # same scheme -> same numbers, without rotation
out_f = SW.linear_1d_run(eta_bump, dx=dx1, dt=dt1, n_steps=144, H=He1, f=f35)   # the library run with rotation (35 N)
assert np.allclose(march(eta_bump, f35, 144), out_f["eta"][-1], atol=1e-12)     # and with rotation
print("✓ the hand-written march reproduces SW.linear_1d_run with and without rotation (to round-off)")   # visible confirmation
""", r"""
Six lines inside the loop are the whole model: continuity with the old velocity, the Coriolis turn of $v$, then the momentum
equation with the new surface. Being the same scheme, it gives the same numbers as the library to round-off.""")
fig(r"""
n_long = 380                                                   # steps for the picture (5.3e5 s, about 6 days: just before the pulses reach the end walls)
runs = [SW.linear_1d_run(eta_bump, dx=dx1, dt=dt1, n_steps=n_long, H=He1, f=f_, save_every=4) for f_ in (0.0, f35)]   # without and with rotation
c1 = GFD.long_wave_speed(He1)                                  # long-wave speed [m/s]
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), sharey=True)
for ax, r_, ttl in zip(axes, runs, ("no rotation: both pulses leave", "35° N: part of the bump stays")):   # one panel each
    pm = ax.pcolormesh(x1/1e3, r_["t"]/86400, r_["eta"], cmap="RdBu_r", vmin=-0.5, vmax=0.5, shading="auto")   # x-t diagram of the surface height
    ax.plot(c1*r_["t"]/1e3, r_["t"]/86400, "--", color=COLORS["ink"], lw=0.8)    # the ray x = +ct
    ax.plot(-c1*r_["t"]/1e3, r_["t"]/86400, "--", color=COLORS["ink"], lw=0.8)   # the ray x = -ct
    ax.set(xlabel="x [km]", title=ttl, xlim=(-2000, 2000))
axes[0].set_ylabel("time [days]")
fig.colorbar(pm, ax=axes, label="surface height η [m]")
plt.show()
print(f"what is left at x = 0 after {runs[1]['t'][-1]/86400:.1f} days: {runs[0]['eta'][-1][n1//2]:+.3f} m without rotation, {runs[1]['eta'][-1][n1//2]:+.3f} m with rotation")   # the part that stays
""",
    see="Two *x–t* (Hovmöller) diagrams of the surface height: distance across, time upward. Left, without rotation: two "
        "straight rays leave the origin. Right, with rotation: the rays spread into wave trains, and a part of the bump stays "
        "at $x=0$.",
    read="In an x–t diagram the slope of a ray is a speed: the dashed lines are $x=\\pm ct$. Whatever still sits on the line "
         "$x=0$ at late time never left — the printed line gives its height.",
    change="…$H$ were quadrupled? The speed doubles, so — with time on the vertical axis — the rays would be half as steep: the "
           "pulses would reach the edge in half the time. (The line has walls at its ends; the picture stops just before the "
           "pulses get there.) Why something stays at all is the subject of `C12`.")
note("N53 [B]", r"""
**Equivalent depth.** A stratified mode with long-wave speed $c$ behaves like a homogeneous layer of depth $H_e$ defined by
$c^2=gH_e$ (13.46). Made precise in `C08`.""")
code(r"""
print(f"a mode with c = 3.6096 m/s has the equivalent depth H_e = c^2/g = {GFD.equivalent_depth(3.6096):.3f} m")   # Eq. (13.46)
""", explain=r"""
**What does this show?** An internal wave that travels at 3.6 m/s behaves, horizontally, exactly like a surface wave on a
layer of water only 1.3 m deep — the depth used in the runs above.""")
nb.code(r"""
if not FAST:                                                   # a slow cell: skipped in FAST runs
    m2d = SW.ShallowWater(64, 64, 2.0e6, 2.0e6, He1, f35, bc="closed")       # a closed 2000 km square basin on a 64 x 64 C-grid
    st0 = SW.geostrophic_state(m2d, np.zeros((64, 64)))                      # a state dictionary {"eta", "u", "v"} at rest
    st0["eta"] = SW.gaussian_bump(m2d, 0.05, 1.0e6, 1.0e6, 1.5e5)            # put a 5 cm bump in the middle
    hist = SW.run(m2d, st0, t_end=20*SW.dt_limit(m2d))                       # twenty time steps of the two-dimensional model
    vol = hist["eta"].sum(axis=(1, 2))                                       # total volume anomaly at each saved step
    print(f"2-D model: {len(hist['t'])} saved steps; volume changes by {abs(vol - vol[0]).max()/abs(vol[0]):.1e} of itself (flux form conserves it)")   # conservation check
else:
    print("FAST run: the two-dimensional model demo is skipped (its cached output is used in C11 and C13)")   # say what was skipped
""", explain=r"""
**What does the code above do?** It shows the interface of our two-dimensional model in six lines: make a grid, put a bump on
it, run twenty steps and check that the volume is conserved to round-off (the continuity equation is in flux form). Everything
else two-dimensional in this notebook comes from cached runs.""", tags=["slow"])
whatif(r"""
…the fluid is not one layer but continuously stratified? Surprisingly little: it splits into a stack of such layers, each with
its own depth (`C08`).""")


# =====================================================================================================================
# A.9  §13.9 — C08 vertical normal modes
# =====================================================================================================================
nb.section("13.9", "Normal Modes in a Continuously Stratified Layer", intro=r"""
**What is this section about?** How a continuously stratified ocean or atmosphere can be replaced by a handful of
shallow-water systems: the vertical structure separates out as an eigenvalue problem set by $N(z)$. ⚠️ Here $z=0$ is the
surface and the flat bottom is at $z=-H$.""")
core("C08", r"Vertical normal modes: $\dfrac{d}{dz}\Big(\dfrac{1}{N^2}\dfrac{d\psi_n}{dz}\Big)+\dfrac{1}{c_n^2}\psi_n=0$ (13.56), each a shallow-water system with $c_n^2\equiv gH_e$ (13.62)",
     "How can a one-layer model say anything about an ocean whose density changes all the way down?")
problem(r"""
When El Niño begins, a bulge in the thermocline crosses the Pacific in about two months — far slower than a tsunami (hours),
yet it is the same kind of wave. The ocean has many ways to wobble: all together, top to bottom (fast), or with the upper ocean
moving against the deep ocean (slow), or in three alternating layers (slower still). Each of these *modes* has a fixed vertical
shape and moves horizontally exactly like a single shallow layer — of a different, much smaller, depth. For the first internal
mode that depth is about a metre.""")
idea(r"""
   mode 0 (barotropic)     mode 1 (first baroclinic)     mode 2
      →  →  →  →               →  →  →                     →  →
      →  →  →  →               →  →                        ←  ←
      →  →  →  →               ←  ←                        ←  ←
      →  →  →  →               ←  ←  ←                     →  →
   no sign change           one sign change            two sign changes
""", r"""
Mode $n$ has $n$ levels where the horizontal velocity changes sign. Each mode has its own speed $c_n$ and its own equivalent
depth $H_e=c_n^2/g$ — the table is filled in by the code below.""")
nb.recap("R19", "Continuity", r"""
Mass conservation: $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$ (13.47).""", where="Ch. 4")
nb.current_core = "C08"
note("N55 [B]", r"""
**Horizontal momentum, linear, with rotation:**
$\dfrac{\partial u}{\partial t}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$ and $\dfrac{\partial v}{\partial t}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$ (13.48)–(13.49).""")
nb.recap("R20", "Hydrostatic balance and the density equation", r"""
The two are $0=-\dfrac{\partial p}{\partial z}-g\rho$ (13.50) and $\dfrac{\partial\rho}{\partial t}-\dfrac{\rho_0N^2}{g}w=0$ (13.51): density at
a point changes only because the background gradient is carried up or down, $w\,d\bar\rho/dz=-\rho_0N^2w/g$.""",
         where=r"Ch. 7, $\frac{d\bar\rho}{dz}=-\frac{\rho_0N^2}{g}\Rightarrow\frac{\partial\rho'}{\partial t}-\frac{N^2\rho_0}gw=0$ (7.131), there without rotation")
nb.current_core = "C08"
code(r"""
print("the five linear hydrostatic equations are consistent (sympy engine) ->", ch13.hydrostatic_linear_set_sympy()["ok"])   # True
""", explain=r"""
**What does this show?** The five equations above — continuity, two momentum equations, hydrostatic balance, density — are
the Start of derivation `D11`; the engine confirms the set as coded.""")
P("P314", "Sturm–Liouville problem: a ladder of eigenvalues, eigenfunctions with n zero crossings, orthogonality", r"""
A string fixed at both ends can only vibrate in certain shapes — one hump, two, three… — each with its own frequency. Any
equation of the form $(p\,\psi')'+\lambda w\,\psi=0$ with conditions at two ends behaves the same way: only special values
$\lambda_0<\lambda_1<\dots$ allow a solution, the $n$-th solution crosses zero $n$ times, and two different solutions are
orthogonal (their weighted product integrates to zero). Here $\lambda=1/c^2$ and $p=1/N^2$.""", code=r"""
z = np.linspace(-1, 0, 2001)                                    # a layer of unit depth
p1, p2 = np.cos(np.pi*z), np.cos(2*np.pi*z)                     # the first two shapes for uniform N and a rigid lid
print("integral of p1*p2 =", abs(round(np.trapezoid(p1*p2, z), 12)), "; integral of p1*p1 =", round(np.trapezoid(p1*p1, z), 3))   # 0.0 and 0.5: orthogonal, not normalised
""")
remind("C08")
nb.md(r"""
**Two boundary conditions first** — derivation `D11` uses them in its last step, so here is where they come from.""")
note("N65 [B]", r"""
**The link the boundary conditions need:**
$w=\dfrac{g(\partial\rho/\partial t)}{\rho_0N^2}=-\dfrac1{\rho_0N^2}\dfrac{\partial^2p}{\partial z\,\partial t}=-\dfrac1{N^2}\sum_{n=0}^\infty\dfrac{\partial p_n}{\partial t}\dfrac{d\psi_n}{dz}$ (13.63).""")
note("N66 [B]", r"""
**Flat bottom**, $w=0$: $\dfrac{d\psi_n}{dz}=0$ at $z=-H$ (13.64).""")
nb.recap("R21", "The linearised free surface", r"""
At the surface $w=\partial\eta/\partial t$ and $p=\rho_0g\eta$ at $z=0$ (the book labels this pair (13.65′)), which combine to
$\partial p/\partial t=\rho_0gw$ at $z=0$.""",
         where=r"Ch. 7: the kinematic condition $\big(\frac{\partial\phi}{\partial z}\big)_{z=0}\cong\frac{\partial\eta}{\partial t}$ (7.18) and the dynamic condition $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}\cong-g\eta$ (7.21)")
nb.current_core = "C08"
note("N67 [B]", r"""
**Hence the surface condition on the modes:** $\dfrac{d\psi_n}{dz}+\dfrac{N^2}{g}\psi_n=0$ at $z=0$ (13.65).""")
P("P315", "Robin (mixed) boundary condition, between Dirichlet and Neumann", r"""
A Dirichlet condition fixes the value ($\psi=0$), a Neumann condition fixes the slope ($\psi'=0$), a Robin condition ties them
together ($\psi'+a\psi=0$). With $a=0$ it is Neumann. Here $a=N^2/g$ is tiny — about 10⁻⁶ per metre in the ocean — so for
internal modes the free surface is almost a rigid lid.""", code=r"""
N, g = 2.7e-3, 9.80665                 # our ocean's buoyancy frequency [1/s] and gravity [m/s^2]
print(f"a = N^2/g = {N**2/g:.1e} per metre; times the depth 4200 m: {N**2/g*4200:.4f}")   # 7.4e-07 and 0.0031
""")
D("D11", ref="13.56")
note("N56 [B]", r"""
**Separation of variables** (step 1 of `D11`): $[u,v,p/\rho_0]=\sum_{n=0}^\infty[u_n,v_n,p_n]\,\psi_n(z)$ (13.52).""")
note("N57 [B]", r"""
**The vertical velocity** (step 2): $w=\sum_{n=0}^\infty w_n\int_{-H}^z\psi_n\,dz$ (13.53).""")
note("N58 [B]", r"""
**The density** (step 3): $\rho=\sum_{n=0}^\infty\rho_n\dfrac{d\psi_n}{dz}$ (13.54).""")
note("N59 [B]", r"""
**The separation constant** (step 5):
$\dfrac{d\psi_n/dz}{N^2\int_{-H}^z\psi_n\,dz}=\dfrac{\rho_0}{g}\dfrac{w_n}{\partial\rho_n/\partial t}\equiv-\dfrac1{c_n^2}$ (13.55).""")
note("N61 [B]", r"""
**Modal continuity** (step 8): $\dfrac{\partial u_n}{\partial x}+\dfrac{\partial v_n}{\partial y}+\dfrac1{c_n^2}\dfrac{\partial p_n}{\partial t}=0$ (13.57).""")
note("N62 [B]", r"""
**Modal momentum** (step 9): $\dfrac{\partial u_n}{\partial t}-fv_n=-\dfrac{\partial p_n}{\partial x}$ and
$\dfrac{\partial v_n}{\partial t}+fu_n=-\dfrac{\partial p_n}{\partial y}$ (13.58)–(13.59).""")
note("N63 [B]", r"""
**The remaining amplitudes** (steps 3 and 7): $p_n=-\dfrac g{\rho_0}\rho_n$ and $w_n=\dfrac1{c_n^2}\dfrac{\partial p_n}{\partial t}$ (13.60)–(13.61).""")
note("N64 [B]", r"""
**Each mode is a shallow-water system** (step 10): identify $p_n\leftrightarrow g\eta$ and $c_n^2\leftrightarrow gH$; the
definition $c_n^2\equiv gH_e$ (13.62) names the equivalent depth.""")
note("N60 [B]", r"""
**Orthogonality** (steps 11–13): $\int_{-H}^0\psi_m\psi_n\,dz=0$ for $m\neq n$, with weight 1, **for the rigid lid and for the
free surface alike** (the book only states it). The surface term $\psi_m(0)\psi_n(0)/g$ belongs to a second, "energy"
relation — $\int\psi_m'\psi_n'/N^2\,dz+\psi_m(0)\psi_n(0)/g=0$ — not to this one.""")
trap("T12", r"""
The modal amplitudes do not share units: $\psi_n$ is a pure number, so $u_n$ and $v_n$ are velocities [m/s] and $p_n$ is
$p/\rho_0$ [m²/s²]; but $w_n$ multiplies an integral of $\psi_n$ over $z$ [m], so $w_n$ is in s⁻¹, and $\rho_n$ multiplies
$d\psi_n/dz$ [1/m], so $\rho_n$ is in kg/m².""")
code(r"""
amp = VM.modal_amplitudes(3.61, 1.0e-5, p_n=0.1, rho0=1027.0)  # mode speed 3.61 m/s; dp_n/dt = 1e-5 m^2/s^3; p_n = 0.1 m^2/s^2
print(f"w_n = (1/c_n^2) dp_n/dt = {amp['w_n']:.2e} 1/s;   rho_n = -(rho0/g) p_n = {amp['rho_n']:.2f} kg/m^2")   # Eqs. (13.61), (13.60)
print("without p_n the density amplitude is not computed:", VM.modal_amplitudes(3.61, 1.0e-5)["rho_n"])    # None
""", explain=r"""
**What does this show?** The units of the trap, in numbers: $w_n$ comes out in s⁻¹ and $\rho_n$ in kg/m². If `p_n` is not
passed, the entry `"rho_n"` is `None` (not zero).""")
note("N68 [C]", r"""
**Uniform $N$, worked case** (stated, as the book writes it out): the equation becomes
$\dfrac{d^2\psi_n}{dz^2}+\dfrac{N^2}{c_n^2}\psi_n=0$ (13.66).""")
note("N69 [C]", r"""
**Its solution:** $\psi_n=A_n\cos\dfrac{Nz}{c_n}+B_n\sin\dfrac{Nz}{c_n}$ (13.67).""")
note("N70 [C]", r"""
**The surface condition** fixes one constant: $B_n=-\dfrac{c_nN}{g}A_n$ (13.68).""")
note("N71 [B]", r"""
**The bottom condition** then gives the eigenvalue condition $\tan\dfrac{NH}{c_n}=\dfrac{c_nN}{g}$ (13.69).""")
P("P316", "graphical roots of a transcendental equation such as tan x = εx", r"""
Some equations have the unknown both inside and outside a function and cannot be solved by algebra. Plot both sides against
$x$: the roots are where the curves cross. Reading the picture also tells you roughly where they are — here, just above 0,
$\pi$, $2\pi$, … — which is all a root-finder such as `brentq` needs.""", code=r"""
eps = 0.0031                                   # N^2 H / g for our ocean
F = lambda X: np.tan(X) - eps/X                # zero where tan X = eps/X
print(f"first root X = {brentq(F, 1e-6, 1.0):.4f}; second root = pi + {brentq(F, np.pi + 1e-9, np.pi + 1.0) - np.pi:.5f}")   # 0.0558 and just above pi
""")
note("N72 [B]", r"""
**The roots, graphically** (our version of the book's figure, `ch13.fig_mode_roots`): with $X=NH/c_n$ the condition
$\tan\frac{NH}{c_n}=\frac{c_nN}{g}$ *(13.69)* reads $\tan X=(N^2H/g)/X$.""")
fig(r"""
figR = ch13.fig_mode_roots()                                   # tan X and (N^2 H/g)/X against X = N H / c_n, roots marked (our ocean)
plt.show()
""",
    see="The branches of $\\tan X$ and the hyperbola $(N^2H/g)/X$ (drawn exaggerated, as the figure says, so that it can be "
        "seen), with the crossings marked.",
    read="The hyperbola is so low that it meets each tangent branch almost at its foot: the roots sit just above "
         "$X=0,\\pi,2\\pi,\\dots$ — one tiny root (the barotropic mode) and then $X\\approx n\\pi$.",
    change="…$g$ were a thousand times smaller? The roots would climb the branches and the lid approximation below would fail.")
nb.recap("R22", "The barotropic mode", r"""
The first root has $NH/c_0\ll1$; there $\tan X\approx X$ gives $c_0=\sqrt{gH}$ (13.70), the long-wave speed of a homogeneous
ocean, with a nearly uniform structure $\psi_0\simeq1-N^2z/g\simeq1$.""", where="Ch. 7, long waves")
nb.current_core = "C08"
slip(6, r'that the first root of $\tan\frac{NH}{c_n}=\frac{c_nN}{g}$ (13.69) occurs "for $NH/c_n=1$"',
     r"$NH/c_0\ll1$ (ours, computed in the cell below: 0.0558 for our ocean).")
note("N73 [B]", r"""
**The baroclinic roots.** Since $c_nN/g\ll1$, the condition is $\tan\dfrac{NH}{c_n}=0$ to a very good approximation, so
$c_n=\dfrac{NH}{n\pi}$, $n=1,2,3,\dots$ (13.71) (`GFD.baroclinic_mode_speed`).""")
note("N74 [B]", r"""
**Numbers with our inputs** ($H=4200$ m, $N=2.7\times10^{-3}$ s⁻¹): computed in the code cell below — a factor of about 56 in
speed, and 3200 in equivalent depth, between the barotropic and the first baroclinic mode.""")
nb.worked_example("the first baroclinic mode by hand", r"""
Take $N=\pi\times10^{-3}\approx3.14\times10^{-3}$ s⁻¹ and $H=4000$ m.

1. $c_1=NH/\pi=(\pi\times10^{-3}\times4000)/\pi=4$ m/s.
2. Equivalent depth $H_e=c_1^2/g\approx16/10=1.6$ m.
3. Mode 2 travels at half that speed, $c_2=c_1/2$; mode 3 at a third, $c_1/3\approx1.33$ m/s.
4. Barotropic: $\sqrt{gH}=\sqrt{40\,000}=200$ m/s.
5. Size of the surface term: $N^2H/g=10^{-5}\times4000/10=0.004\ll1$ — the surface is almost a lid for these modes.""")
P("P317", "scipy.linalg.eigh_tridiagonal for a discretised eigenproblem", r"""
Replace $\psi$ by its values on a grid and the second derivative by differences: the differential eigenproblem becomes a
matrix one, with only the diagonal and its two neighbours filled. `eigh_tridiagonal(d, e)` returns all eigenvalues (sorted)
and eigenvectors of a symmetric tridiagonal matrix from its diagonal `d` and off-diagonal `e`.""", code=r"""
n = 200; h = 1.0/n                                                     # interior points of a string of unit length
lam, vec = eigh_tridiagonal(2*np.ones(n-1)/h**2, -np.ones(n-2)/h**2)   # the matrix of minus the second derivative
print("sqrt(eigenvalue)/pi of the first three modes:", np.round(np.sqrt(lam[:3])/np.pi, 4))   # 1, 2, 3: the modes of a string
""")
code(r"""
N_o, H_o = inp["ocean_N"], inp["ocean_H"]                      # our ocean: N = 2.7e-3 1/s, H = 4200 m
z_m = np.linspace(-H_o, 0.0, 401)                              # 401 nodes from the bottom to the surface [m]
N2_u = np.full_like(z_m, N_o**2)                               # uniform N^2 at the nodes [1/s^2]
m_free = VM.vertical_modes(z_m, N2_u, n_modes=4, lid="free")   # numerical modes of Eq. (13.56) with a free surface: index 0 = barotropic
m_exact = VM.modes_uniform_N(N_o, H_o, n_modes=4)              # exact speeds for uniform N (roots of Eq. (13.69) by brentq)
print("numerical c_n [m/s]:", np.round(m_free.c, 4), "  exact:", np.round(m_exact.c, 4))   # the ladder of speeds
print("equivalent depths H_e = c_n^2/g [m]:", np.round(m_free.He, 3))                      # Eq. (13.62)
print(f"barotropic: exact c_0 is {100*(m_exact.c[0]/GFD.long_wave_speed(H_o) - 1):.3f} % above sqrt(gH); first root X_0 = N H/c_0 = {N_o*H_o/m_exact.c[0]:.5f}")   # slip #6 in numbers
print(f"baroclinic: N H/(n pi) = {[round(GFD.baroclinic_mode_speed(N_o, H_o, n_), 4) for n_ in (1, 2, 3)]} m/s; ratio c_0/c_1 = {m_exact.c[0]/m_exact.c[1]:.0f}, H_e ratio = {m_exact.He[0]/m_exact.He[1]:.0f}")   # Eq. (13.71)
print("rigid-lid error (c_free - c_rigid)/c_rigid for n = 1, 2, 3:", [f"{VM.rigid_lid_error(N_o, H_o, n=n_):.1e}" for n_ in (1, 2, 3)])   # how good the lid is
print(f"Rossby radius of mode 1 at 35 N: c_1/f = {GFD.rossby_radius(m_free.c[1], f35)/1e3:.1f} km")   # used in C12
""", explain=r"""
1. `VM.vertical_modes(z, N2, lid=…)` solves the eigenproblem on a grid and returns a `Modes` named tuple with fields `c`
   (speeds, fastest first), `psi` (shapes, normalised to $\psi(0)=1$), `He`, `z`, `weights`.
2. `VM.modes_uniform_N` gives the exact speeds for uniform $N$ from the roots of $\tan\frac{NH}{c_n}=\frac{c_nN}{g}$ *(13.69)*; the two agree to
   four or five digits.
3. The barotropic speed is within a twentieth of a per cent of $\sqrt{gH}$, and its root is $X_0\approx0.056$ — far below 1 (slip #6).
4. The baroclinic speeds are $NH/(n\pi)$ to a few parts in 10⁴; the first is some 56 times slower than the barotropic mode.
5. `VM.rigid_lid_error` measures what replacing the free surface by a lid costs: 3 parts in 10⁴ for $n=1$, less for higher modes.""")
scratch(r"""
# From scratch: the rigid-lid modes for uniform N as a small symmetric matrix eigenproblem
n_c = 400                                                      # number of grid intervals
h_c = H_o/n_c                                                  # grid spacing [m]
K_st = (np.diag(2*np.ones(n_c + 1)) - np.diag(np.ones(n_c), 1) - np.diag(np.ones(n_c), -1))/h_c   # stiffness matrix of -d2/dz2 …
K_st[0, 0] = K_st[-1, -1] = 1.0/h_c                            # … with zero slope at both ends (half cells there)
w_c = np.full(n_c + 1, h_c); w_c[[0, -1]] = h_c/2              # the weights of the trapezoid rule
A_sym = K_st/np.sqrt(np.outer(w_c, w_c))                       # symmetric form W^(-1/2) K W^(-1/2)
lam_c = np.linalg.eigh(A_sym)[0]                               # eigenvalues N^2/c^2 in increasing order; the first is 0 (a constant)
c_mine = N_o/np.sqrt(lam_c[1])                                 # speed of the first baroclinic mode [m/s]
m_rigid = VM.vertical_modes(z_m, N2_u, lid="rigid")            # the library's rigid-lid modes: index 0 is the FIRST BAROCLINIC mode
assert np.isclose(c_mine, GFD.baroclinic_mode_speed(N_o, H_o), rtol=1e-4)   # same as N H / pi
assert np.isclose(c_mine, m_rigid.c[0], rtol=1e-4)             # same as the library (note the index 0)
print(f"✓ c_1 by hand = {c_mine:.4f} m/s;  N H/pi = {GFD.baroclinic_mode_speed(N_o, H_o):.4f};  VM rigid-lid c[0] = {m_rigid.c[0]:.4f} m/s")   # visible confirmation
""", r"""
The mode equation for uniform $N$ with a rigid lid is "minus the second derivative of $\psi$ = $(N^2/c^2)\,\psi$ with zero
slope at both ends" — a small symmetric matrix problem. **Mind the index:** a rigid-lid `Modes` has no barotropic entry, so
its first element (`c[0]`) is the first baroclinic mode; free-surface modes keep the barotropic mode at index 0 and the first
baroclinic one at index 1.""")
note("N75 [B]", r"""
**The rigid-lid approximation.** Replace the surface condition by $w=0$ at $z=0$, that is, the slope of $\psi_n$ vanishes
there. Then $\psi_n=\cos\dfrac{n\pi z}{H}$, $n=0,1,2,\dots$: the baroclinic speeds change by 3 parts in 10⁴ or less (cell above)
and the barotropic mode disappears ($c_0\to\infty$). What the lid does **not** mean: the surface pressure still varies under
it — the lid pushes back.""")
P("P318", "projecting a profile on modes: the inner product of two functions", r"""
Two vectors are orthogonal when their dot product is zero; two functions when the integral of their product is zero. Because
the modes are orthogonal, the amount of mode $n$ in any profile $F(z)$ is found by one integral,
$a_n=\int F\psi_n\,dz\big/\int\psi_n^2\,dz$ — like reading off a component along an axis.""", code=r"""
z = np.linspace(-1, 0, 2001); F = 1 + z                        # a linear profile on a layer of unit depth
p1 = np.cos(np.pi*z)                                           # the first rigid-lid mode
print("mode-1 coefficient a_1 =", round(np.trapezoid(F*p1, z)/np.trapezoid(p1*p1, z), 4))   # 0.4053 = 4/pi^2
""")
code(r"""
O_psi = VM.orthogonality_matrix(m_free)                        # normalised integrals of psi_m psi_n, weight 1, no surface term
O_en = VM.orthogonality_matrix(m_free, kind="energy", N2=N2_u) # the second ("energy") relation, which carries the surface term
off = lambda M: abs(M - np.diag(np.diag(M))).max()             # largest off-diagonal entry of a matrix
print(f"free-surface modes: weight-1 orthogonality off-diagonal max = {off(O_psi):.1e}; energy relation = {off(O_en):.1e}")   # both vanish
profile = np.exp(z_m/500.0)                                    # a surface-trapped current profile, e-folding depth 500 m
a_n = VM.project(m_free, profile)                              # its modal coefficients a_n
print("coefficients a_0 … a_3:", np.round(a_n, 3), f"; rms misfit of the 4-mode sum = {np.sqrt(np.mean((VM.reconstruct(m_free, a_n) - profile)**2)):.3f}")   # how much four modes capture
print(f"w-structure of mode 1 (integral of psi_1) peaks at z = {z_m[np.argmax(abs(VM.w_structure(m_free, 1)))]:.0f} m; rho-structure dpsi_1/dz peaks at z = {z_m[np.argmax(abs(VM.rho_structure(m_free, 1)))]:.0f} m")   # mid-depth
print(f"rigid-lid modes projected on a constant profile: {abs(VM.project(m_rigid, np.ones_like(z_m))).max():.1e} (they cannot hold a depth mean)")   # add profile.mean() back
""", explain=r"""
1. `VM.orthogonality_matrix(m)` is the weight-1 integral $\int\psi_m\psi_n\,dz$, normalised; it is the identity matrix for the
   free-surface modes too — no surface term is needed. With `kind="energy"` it is the second relation, which does carry the
   term $\psi_m(0)\psi_n(0)/g$.
2. `VM.project` gives the coefficients of a profile on the modes; `VM.reconstruct` sums them back.
3. `VM.w_structure(m, 1)` and `VM.rho_structure(m, 1)` are the integral and the derivative of $\psi_1$: where mode 1 moves
   the water up and down the most, and where it changes the density the most — at mid-depth, where its horizontal velocity
   changes sign.
4. Rigid-lid modes have no barotropic member, so a constant projects to zero: add the depth mean of the profile back by hand.""")
note("N76 [B]", r"""
**The shapes** (our version of the book's figure): the first three modes for uniform $N$ (cosines) and for a thermocline
profile, `ch13.thermocline_N2(z)` — ours, the kind of $N(z)$ shown in §13.2. The cell also prints a *WKB estimate* of the
first speed. WKB is the approximation for a slowly varying medium (its primer is in `C14`): if $N(z)$ changed little over one
vertical wavelength, the uniform-$N$ result $c_n=NH/(n\pi)$ would hold with $NH$ replaced by its integral,
$c_n\approx\frac1{n\pi}\int_{-H}^0N\,dz$ (`VM.wkb_mode_speed`).""")
fig(r"""
N2_th = ch13.thermocline_N2(z_m)                               # N^2(z) with a thermocline at -300 m (defaults; ours)
m_th = VM.vertical_modes(z_m, N2_th, n_modes=4, lid="free")    # its modes
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=True)
for n_, col in zip((1, 2, 3), (COLORS["accent"], COLORS["teal"], COLORS["rose"])):   # the first three baroclinic modes
    a.plot(m_free.psi[n_], z_m, color=col, label=f"mode {n_}: c = {m_free.c[n_]:.2f} m/s")   # uniform N: cosines
    b.plot(m_th.psi[n_], z_m, color=col, label=f"mode {n_}: c = {m_th.c[n_]:.2f} m/s")       # thermocline profile
b.fill_betweenx(z_m, 0, np.sqrt(N2_th)/np.sqrt(N2_th).max(), color=COLORS["blue"], alpha=0.15, label="N(z), scaled")   # the stratification behind
for ax, ttl in ((a, "uniform N: cosines"), (b, "thermocline: crossings crowd into it")):
    ax.axvline(0, color=COLORS["muted"], lw=0.6)
    ax.set(xlabel=r"mode shape $\psi_n$ (normalised to 1 at the surface)", title=ttl)
    ax.legend(fontsize=7, loc="lower right")
a.set_ylabel("z [m]")
plt.show()
zero_x = z_m[np.where(np.diff(np.sign(m_th.psi[1])) != 0)[0][0]]   # the depth where mode 1 changes sign
print(f"thermocline modes: c = {np.round(m_th.c[1:4], 3)} m/s; H_e of mode 1 = {m_th.He[1]:.3f} m; Rossby radius at 35 N = {GFD.rossby_radius(m_th.c[1], f35)/1e3:.1f} km")   # speeds and scales
print(f"mode 1 changes sign at z = {zero_x:.0f} m; WKB estimate of c_1 = {VM.wkb_mode_speed(z_m, np.sqrt(N2_th), n=1):.2f} m/s ({100*(VM.wkb_mode_speed(z_m, np.sqrt(N2_th), n=1)/m_th.c[1] - 1):+.0f} %)")   # how good the slowly-varying estimate is
""",
    see="The first three internal modes for uniform stratification (left) and for a profile with a sharp thermocline (right, "
        "with $N(z)$ shaded in blue).",
    read="Mode $n$ crosses zero $n$ times. With a thermocline the crossings crowd into it, the deep ocean moves as one slab, "
         "and the modes are slower. The printed WKB estimate $\\int N\\,dz/(n\\pi)$ misses the first speed by a quarter — a "
         "sharp thermocline is not a slowly varying medium.",
    change="…the thermocline were deeper? Next figure.")
nb.plotly(r"""
z_f6 = np.linspace(-H_o, 0.0, 141)                             # a coarser grid (141 nodes) for fifteen eigen-solves

def f6(depth):                                                 # curves for one thermocline depth [m]
    N2_ = ch13.thermocline_N2(z_f6, z_t=-depth)                # N^2(z) with the thermocline centred at -depth
    md = VM.vertical_modes(z_f6, N2_, n_modes=4, lid="free")   # its modes
    top = z_f6 >= -2000.0                                      # only the top 2000 m are drawn
    return {"N(z), scaled to 1": ((np.sqrt(N2_)/np.sqrt(N2_).max())[top], z_f6[top]),
            "mode 1": (md.psi[1][top], z_f6[top]), "mode 2": (md.psi[2][top], z_f6[top]), "mode 3": (md.psi[3][top], z_f6[top])}

figF6 = slider_figure(f6, "thermocline depth", np.arange(100.0, 801.0, 50.0), unit="m", xlabel="ψ_n (and N scaled)",
                      ylabel="z [m]", title="A deeper thermocline is a faster first mode", xrange=(-4, 4), yrange=(-2000, 0))   # 15 slider positions
figF6.show()
for d_ in (150.0, 300.0, 600.0):                               # three thermocline depths [m]
    md_ = VM.vertical_modes(z_f6, ch13.thermocline_N2(z_f6, z_t=-d_), n_modes=4, lid="free")   # modes for each
    print(f"thermocline at {d_:.0f} m: c_1 = {md_.c[1]:.3f} m/s, H_e = {md_.He[1]:.3f} m, c_2 = {md_.c[2]:.3f}, c_3 = {md_.c[3]:.3f} m/s")   # the speeds the slider cannot show in its legend
""", explain=r"""
1. `f6` rebuilds $N^2(z)$ with the thermocline at the chosen depth, solves the eigenproblem and returns the first three
   baroclinic modes (only the top 2000 m are shown).
2. The printed lines give the speeds and the equivalent depth for three positions of the slider.""")
nb.figure_notes(
    see="The stratification (scaled) and the first three internal modes in the top 2000 m, for the chosen thermocline depth.",
    read="The zero crossing of mode 1 follows the thermocline down, and its speed grows (printed lines). A deeper thermocline "
         "is a faster first mode — the El Niño thermocline signal in one slider.",
    change="…the thermocline were sharper? The modes would bend more abruptly there; the WKB estimate would be worse still.")
note("N77 [B]", r"""
**The fine print.** The decomposition needs a flat bottom and no sheared mean current; it is hydrostatic, so it holds for
$\omega\ll N$ and the shapes do not depend on frequency; it is valid with or without $f$ and $\beta$.""")
explainer("vertical_modes", "How can one layer stand for a stratified ocean?",
          "reshaping $N(z)$ and watching the shapes, speeds and equivalent depths respond shows what the eigenproblem does; "
          "stepping $n$ shows the zero crossings appear one by one; the lid toggle removes the barotropic mode and barely moves "
          "the others.",
          ["Step n from 0 to 3 and count the zero crossings.",
           "Deepen the thermocline and watch c₁ and the Rossby radius grow.",
           "Toggle the rigid lid: the top rung of the ladder disappears, the others do not move.",
           "Click a depth to see the two terms of the mode equation cancel there."])
whatif(r"""
…we now ask what waves one of these shallow-water systems supports when the earth turns? One cubic equation answers for all
of them (`C09`).""")


# =====================================================================================================================
# A.10  §13.10 — C09 the complete dispersion relation
# =====================================================================================================================
nb.section("13.10", "High- and Low-Frequency Regimes in Shallow-Water Equations", intro=r"""
**What is this section about?** One equation that contains every wave of a rotating shallow layer, and how to tell, from the
frequency alone, which terms matter.""")
core("C09", r"The complete dispersion relation of rotating shallow water: $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$, with $K^2=k^2+l^2$, $c=\sqrt{gH}$ (13.76)",
     "Poincaré, Kelvin, Rossby — how many different kinds of wave does one layer of water really have?")
problem(r"""
A weather forecast must not let fast gravity waves swamp the slow weather; a tide model wants exactly those fast waves. Both
run the same equations. The waves separate by frequency: two fast ones (period hours) that are gravity waves bent by rotation,
and one slow one (period days to years) that exists only because the Coriolis parameter changes with latitude. Knowing which
is which tells a modeller which terms may be dropped.""")
idea(words=r"""
The cubic has four terms, each with its own physics:

| Term | What it stands for | Fast waves ($\omega>f_0$) | Slow wave ($\omega\ll f_0$) |
|---|---|---|---|
| $\omega^3$ | pure change in time | large — balances the next two | negligible |
| $-c^2K^2\omega$ | gravity (the surface slope) | large | large |
| $-f_0^2\omega$ | rotation | large | large |
| $-c^2\beta k$ | the change of $f$ with latitude | negligible | balances the two above |""")
remind("C09")
D("D12", ref="13.75")
note("N78 [C]", r"""
**Steps 2–3 of `D12`:** $\dfrac{\partial^2u}{\partial t^2}-f\dfrac{\partial v}{\partial t}=gH\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ and
$\dfrac{\partial^2v}{\partial t^2}+f\dfrac{\partial u}{\partial t}=gH\dfrac{\partial}{\partial y}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ (13.72)–(13.73).""")
note("N79 [C]", r"""
**Step 4:** $\dfrac{\partial^3v}{\partial t^3}+f\Big[f\dfrac{\partial v}{\partial t}+gH\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)\Big]=gH\dfrac{\partial^2}{\partial y\,\partial t}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ (13.74).""")
note("N80 [B]", r"""
**The linear vorticity equation on the β-plane** (step 6):
$\dfrac{\partial}{\partial t}\Big(\dfrac{\partial u}{\partial y}-\dfrac{\partial v}{\partial x}\Big)-f_0\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)-\beta v=0$
— the first appearance of the idea that `C13` turns into potential vorticity.""")
note("N81 [B]", r"""
**One equation for $v$** (the result of `D12`):
$\dfrac{\partial^3v}{\partial t^3}-gH\dfrac{\partial}{\partial t}\nabla_H^2v+f_0^2\dfrac{\partial v}{\partial t}-gH\beta\dfrac{\partial v}{\partial x}=0$ (13.75).""")
trap("T10", r"""
The parameter $f$ is replaced by the constant $f_0$ everywhere except where it is differentiated (that derivative is $\beta$).""")
code(r"""
print("the single equation for v (13.75) follows from the linear set on a beta-plane ->", ch13.v_equation_sympy()["ok"])   # True
""", explain=r"""
**What does this show?** The engine carries out the eliminations of `D12` on symbolic fields and is left with nothing.""")
P("P319", "three real roots of a cubic: the discriminant and the trigonometric form", r"""
A cubic $\omega^3+p\omega+q=0$ with real $p$, $q$ has three real roots exactly when its discriminant $-4p^3-27q^2$ is positive
(this needs $p<0$). They are then given without complex arithmetic by
$\omega_j=2\sqrt{-p/3}\,\cos\big[\tfrac13\arccos\big(\tfrac{3q}{2p}\sqrt{-3/p}\big)-\tfrac{2\pi j}3\big]$, $j=0,1,2$. When one root is
thousands of times smaller than the others this form keeps its digits; a general polynomial solver may not.""", code=r"""
p, q = -7.0, 6.0                                               # w^3 - 7w + 6 = (w - 1)(w - 2)(w + 3)
print("discriminant =", -4*p**3 - 27*q**2)                     # 400 > 0: three real roots
m = 2*np.sqrt(-p/3); th = np.arccos(3*q/(p*m))/3               # the amplitude and the angle of the trigonometric form
print("roots:", np.sort(m*np.cos(th - 2*np.pi*np.arange(3)/3)))   # [-3.  1.  2.]
""")
D("D13", ref="13.76")
note("N82 [B]", r"""
**All three roots are real** (steps 4–5 of `D13`): the discriminant $4(c^2K^2+f_0^2)^3-27(c^2\beta k)^2$ is positive under
β-plane scaling (the book says only that this "can be shown").""")
trap("T9", r"""
"Two superinertial, one subinertial" is true under β-plane scaling; the two fast roots have opposite signs and only satisfy
$\lvert\omega\rvert>\lvert f_0\rvert$, not $\omega\gg f$.""")
note("N83 [B]", r"""
**Three regimes** (steps 6–7 of `D13`), from the two scale ratios $\dfrac{c^2\beta k}{c^2\omega K^2}\sim\dfrac\beta{\omega K}$ and
$\dfrac{\omega^3}{c^2\beta k}\ll1$: high frequency (drop $\beta$ and, further up, $f_0$), near-inertial (drop $\beta$), low frequency
(drop $\omega^3$).""")
slip(7, r"$\omega\gg f$ as the range in which the first term of $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76) is negligible",
     r"$\omega\ll f$.")
P("P320", "reading a dispersion diagram with several branches: signed ω for signed k, a logarithmic frequency axis", r"""
A dispersion diagram plots frequency against wavenumber; each curve ("branch") is one kind of wave. The slope of the line from
the origin to a point is the phase speed, the slope of the curve itself the group velocity. We let $k$ carry the direction
($k>0$ eastward) and keep $\omega\ge0$ for the plot; when branches differ by factors of a thousand, the frequency axis must be
logarithmic or the slow branch hides on the axis.""", code=r"""
k = np.array([1e-7, 1e-6, 1e-5])                               # three wavenumbers [rad/m]
print("fast branch [1e-4 1/s]:", np.round(1e4*np.sqrt(1e-8 + 4e4*k**2), 3))   # flattens to f = 1e-4 as k -> 0
print("slow branch [1e-6 1/s]:", np.round(1e6*2e-11*k/(k**2 + 2.5e-13), 2))   # rises, peaks, falls
""")
nb.worked_example("fast and slow root by hand", r"""
Our 35° N inputs, rounded: $f_0=8.4\times10^{-5}$ s⁻¹, $\beta=1.9\times10^{-11}$ m⁻¹ s⁻¹, external speed $c=203$ m/s
($c^2=4.12\times10^{4}$ m² s⁻²), $k=2\times10^{-6}$ m⁻¹ (wavelength 3140 km), $l=0$.

1. **Fast:** drop the β term → $\omega^2=f_0^2+c^2k^2=0.07\times10^{-7}+1.65\times10^{-7}=1.72\times10^{-7}$ → $\omega=\pm4.15\times10^{-4}$ s⁻¹ (period 4.2 h).
2. **Slow:** drop $\omega^3$ → $\omega=-c^2\beta k/(c^2k^2+f_0^2)=-(4.12\times10^{4}\times1.9\times10^{-11}\times2\times10^{-6})/(1.72\times10^{-7})=-9.1\times10^{-6}$ s⁻¹
   (period 8 days; minus = westward for $k>0$).
3. Was dropping $\omega^3$ fair? Compare $(9.1\times10^{-6})^3=7.5\times10^{-16}$ against $c^2\beta k=1.57\times10^{-12}$: one part in two thousand.
4. The sum of the three roots must be 0 (there is no $\omega^2$ term): so the fast pair is not exactly symmetric; the slow root
   is the difference.""")
code(r"""
c_ext = GFD.long_wave_speed(H_o)                               # external (barotropic) long-wave speed of a 4200 m ocean [m/s]
k31 = 2*np.pi/3.1e6                                            # eastward wavenumber of a 3100 km wave [rad/m]
w3 = GFD.shallow_water_omega(k31, 0.0, c_ext, f35, beta35)     # Eq. (13.76): the three real roots, in ascending order [rad/s]
print("external mode, three roots [1/s]:", [f"{w_:+.4e}" for w_ in w3])                      # fast westward, slow westward, fast eastward
print(f"periods: fast {2*np.pi/w3[2]/3600:.2f} h, slow {2*np.pi/abs(w3[1])/86400:.2f} days; slow/fast = {abs(w3[1])/w3[2]:.3f}")   # decades apart
print("discriminant > 0:", GFD.shallow_water_discriminant(k31, 0.0, c_ext, f35, beta35) > 0, "| regime of the slow root:", GFD.shallow_water_regime(w3[1], f35))   # three real roots; which term to drop
sz = GFD.dispersion_term_sizes(k31, 0.0, c_ext, f35, beta35, w3[1])   # the four terms of the cubic at the slow root
print("terms at the slow root:", {k_: f"{v_:+.2e}" for k_, v_ in sz.items()})                # omega^3 is by far the smallest
w3b = GFD.shallow_water_omega(k31, 0.0, 3.6096, f35, beta35)   # the same wave on the first baroclinic mode (c = 3.61 m/s)
print("first baroclinic mode [1/s]:", [f"{w_:+.4e}" for w_ in w3b], f"-> {2*np.pi/w3b[2]/3600:.1f} h and {2*np.pi/abs(w3b[1])/86400/365.25:.1f} years")   # even further apart
br = GFD.shallow_water_branches(k31, 0.0, c_ext, f35, beta35)  # the same roots by name, plus the Kelvin line c k
print("branches:", {k_: f"{v_:+.3e}" for k_, v_ in br.items()})                              # Poincare (two), Rossby, Kelvin
k20 = 2*np.pi/2.0e7                                            # a 20 000 km wave at 12 N: outside the range where the cubic has three real roots
disc20 = GFD.shallow_water_discriminant(k20, 0.0, c_ext, f12, beta12)   # the flag to test first
w20 = GFD.shallow_water_omega(k20, 0.0, c_ext, f12, beta12)    # there the function returns "not a number" for all three roots
print(f"12 N, 20 000 km, external mode: discriminant = {disc20:.2e} (negative) -> three real roots exist: {not np.isnan(w20).any()}")   # test the discriminant first
""", explain=r"""
1. `GFD.shallow_water_omega(k, l, c, f0, beta)` returns the three real roots of $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ *(13.76)* by the trigonometric form
   of the primer.
2. For a 3100 km wave on the external mode the fast roots have periods of about 4 hours and the slow one about 8 days.
3. `GFD.dispersion_term_sizes` evaluates the four terms at a root: at the slow root $\omega^3$ is the smallest by orders of
   magnitude; `GFD.shallow_water_regime` names the term that may be dropped.
4. On the first baroclinic mode ($c=3.6$ m/s) the same wave has periods of about 21 hours and nearly 3 years.
5. **The caveat:** where `GFD.shallow_water_discriminant` is negative the cubic has only one real root, and
   `GFD.shallow_water_omega` returns NaN ("not a number") for all three. Test the discriminant first; NaN is the function's
   way of saying "no three real roots here". The figure after next explains where that happens.""")
scratch(r"""
# From scratch: the same cubic handed to a general polynomial solver
roots_mine = np.sort(np.roots([1.0, 0.0, -(c_ext**2*k31**2 + f35**2), -c_ext**2*beta35*k31]).real)   # coefficients of w^3, w^2, w, 1
assert np.allclose(roots_mine, w3, rtol=1e-9)                  # same three roots as the trigonometric form
assert abs(sum(w3))/max(abs(np.array(w3))) < 1e-12             # no w^2 term -> the roots sum to zero
print("✓ np.roots gives", [f"{r_:+.4e}" for r_ in roots_mine], f"1/s; the three roots sum to {sum(w3):.1e}")   # visible confirmation
""", r"""
`np.roots` finds the same three roots, and their sum vanishes as it must for a cubic with no $\omega^2$ term — which is why
the two fast roots are not exactly equal and opposite.""")
fig(r"""
ks = np.geomspace(2*np.pi/4.0e7, 2*np.pi/2.0e5, 600)           # wavenumbers from a 40 000 km to a 200 km wave [rad/m]
roots_all = np.array([GFD.shallow_water_omega(k_, 0.0, c_ext, f35, beta35) for k_ in ks])   # the three roots for each (35 N, external mode)
cases = []                                                     # (label, wavenumber, root) for three frequencies
for target, col_ in ((3.0, 2), (1.02, 2), (0.1, 1)):           # 3f and about f on the fast branch; 0.1 f on the slow branch
    j_ = np.argmin(abs(abs(roots_all[:, col_])/f35 - target))  # the wavenumber whose root is closest to the target
    cases.append((f"ω = {abs(roots_all[j_, col_])/f35:.2f} f", ks[j_], roots_all[j_, col_]))
fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.6), sharey=True)
for ax, (lab, k_, w_) in zip(axes, cases):                     # one panel per frequency
    t4 = GFD.dispersion_term_sizes(k_, 0.0, c_ext, f35, beta35, w_)   # the four terms of Eq. (13.76) at this root
    vals = [abs(t4["omega3"]), abs(t4["gravity"]), abs(t4["rotation"]), abs(t4["beta"])]   # their magnitudes
    ax.bar(range(4), vals, color=[COLORS["accent"], COLORS["orange"], COLORS["teal"], COLORS["amber"]])   # purple, orange, teal, amber
    ax.set_xticks(range(4), ["ω³", "c²K²ω", "f₀²ω", "c²βk"])   # the names of the four terms
    ax.set(yscale="log", title=f"{lab}  (λ = {2*np.pi/k_/1e3:.0f} km)")
axes[0].set_ylabel("size of the term [s$^{-3}$]")
fig.suptitle("Each regime is 'drop the shortest bar'")
plt.show()
for nm, c_, f_, b_ in (("35 N, external", c_ext, f35, beta35), ("35 N, first baroclinic", 3.6096, f35, beta35), ("12 N, external", c_ext, f12, beta12)):   # three settings
    print(f"{nm}: beta c / f0^2 = {b_*c_/f_**2:.4f}  ->  three real roots for every wavelength: {b_*c_/f_**2 < 1}")   # the condition of the caveat
""",
    see="The magnitudes of the four terms of the cubic at three frequencies: three times $f$ and just above $f$ on the fast "
        "branch, and a tenth of $f$ on the slow branch.",
    read="At $3f$ the $\\omega^3$ bar is among the tallest and the β bar is the shortest; at $0.1f$ the $\\omega^3$ bar is the "
         "smallest. Each regime of the book is \"drop the shortest bar\" — and slip #7 is settled by the right-hand panel.",
    change="…$\\beta=0$? The amber bar vanishes and the slow root becomes $\\omega=0$ — a steady geostrophic flow, the state "
           "that the adjustment of `C12` ends in.")
nb.md(r"""
> ⚠️ **An honest caveat (ours, computed).** All three roots are real for every wavelength exactly when $\beta c<f_0^2$, i.e.
> when the Rossby radius $c/f_0$ is smaller than $R\tan\theta_0$. The lines printed above give $\beta c/f_0^2$ for three
> settings: it is below 1 at 35° N for both modes, but far above 1 for the external mode at 12° N. There the discriminant is
> negative for a band of very long waves and the cubic has no three real roots — an artefact of freezing $f$ at $f_0$, not a
> real instability. `GFD.shallow_water_discriminant` is the flag; `GFD.shallow_water_omega` returns NaN there, and the slider
> figure below simply leaves that band blank.""")
nb.plotly(r"""
k_f7 = np.geomspace(2*np.pi/4.0e7, 2*np.pi/1.0e5, 50)          # wavenumbers from a 40 000 km to a 100 km wave [rad/m]

def f7(lat_deg):                                               # branches for one latitude [degrees]
    f_ = GFD.coriolis_parameter(np.deg2rad(lat_deg)); b_ = GFD.beta_parameter(np.deg2rad(lat_deg))   # f0 and beta there
    out = {}
    for c_, nm in ((c_ext, "external mode"), (3.6096, "first baroclinic mode")):   # two long-wave speeds
        w_ = np.array([GFD.shallow_water_omega(k_, 0.0, c_, f_, b_) for k_ in k_f7])   # three roots per wavenumber (NaN where there are not three)
        out[f"Poincaré (fast), {nm}"] = (k_f7, w_[:, 2])       # the eastward fast root
        out[f"Rossby (slow), {nm}"] = (k_f7, abs(w_[:, 1]))    # magnitude of the slow root
        out[f"Kelvin line ω = ck, {nm}"] = (k_f7, GFD.kelvin_omega(k_f7, c_))   # explained in C11
    out["ω = f₀"] = (k_f7, f_ + 0*k_f7)                        # the inertial frequency
    return out

figF7 = slider_figure(f7, "latitude", np.arange(5.0, 75.1, 5.0), unit="°", xlabel="wavenumber k [rad/m]", ylabel="|ω| [rad/s]",
                      title="Fast and slow branches are decades apart", active=6)   # 15 slider positions, starting at 35 degrees
figF7.update_xaxes(type="log", range=[np.log10(k_f7[0]), np.log10(k_f7[-1])])       # logarithmic wavenumber axis
figF7.update_yaxes(type="log", range=[-9.5, -1.5])                                  # logarithmic frequency axis
figF7.show()
""", explain=r"""
1. `f7` computes, for one latitude and for two long-wave speeds, the fast root, the magnitude of the slow root and the
   straight Kelvin line $\omega=ck$ (its meaning is the subject of `C11`).
2. Both axes are logarithmic. Where the cubic has no three real roots the curves are left blank (the caveat above).""")
nb.figure_notes(
    see="For each mode: a fast branch that flattens onto $\\omega=f_0$ for long waves, a slow branch three or more decades "
        "below it that rises, peaks and falls, and the straight Kelvin line. At low latitudes the external-mode curves have a "
        "gap at the longest waves.",
    read="The empty band of frequencies between the fast and slow branches is the reason filtered (\"quasi-geostrophic\") "
         "models work: one can drop the fast waves without touching the slow one.",
    change="…you slide toward the pole? The parameter $f_0$ rises (the floor of the fast branch moves up) and β falls (the slow branch sinks).")
whatif(r"""
…we look at the fast roots alone, on an f-plane? They are the gravity waves of Chapter 7 with a floor under their frequency
(`C10`).""")

# =====================================================================================================================
# A.11  §13.11 — C10 Poincaré waves
# =====================================================================================================================
nb.section("13.11", "Gravity Waves with Rotation", intro=r"""
**What is this section about?** What rotation does to a long gravity wave: it cannot oscillate more slowly than $f$, it becomes
dispersive, and the water under it moves in ellipses.""")
core("C10", r"Poincaré waves: $\omega^2=f^2+gHK^2$, $K=\sqrt{k^2+l^2}$ (13.82)",
     "What does the earth's rotation do to a long gravity wave?")
problem(r"""
After a storm passes, current meters in the open ocean show the water going round in circles a few kilometres across, once
every 17 to 20 hours at mid-latitudes, for days. Nothing is pushing it; it is coasting, and the Coriolis force keeps bending
its path. That circle is the longest possible gravity wave. Shorter waves — the tides among them — are mixtures: partly a
gravity wave sloshing back and forth, partly this turning.""")
idea(words=r"""
Two restoring forces, two frequencies. Gravity alone gives $\omega=cK$; rotation alone gives $\omega=f$; together
$\omega^2=f^2+c^2K^2$ — add the squares, like the sides of a right triangle.""")
note("N84 [B]", r"""
**Plane-wave amplitudes on the f-plane:** $-i\omega\hat u-f\hat v=-ikg\hat\eta$, $-i\omega\hat v+f\hat u=-ilg\hat\eta$,
$-i\omega\hat\eta+iH(k\hat u+l\hat v)=0$ (13.77)–(13.79). These are the Start of `D14`.""")
P("P321", "polarisation relations: solving a 2 × 2 complex system for the velocity amplitudes", r"""
For a plane wave every field is a complex amplitude times the same $e^{i(kx+ly-\omega t)}$. The momentum equations become two
linear equations for the two velocity amplitudes in terms of the height amplitude. Solving them tells how the
current is oriented and timed relative to the surface — the wave's "polarisation". A factor $i$ means a quarter-period shift. *Cramer's rule* for two equations
$a\,x+b\,y=e$, $c\,x+d\,y=g$: $x=(ed-bg)/(ad-bc)$ and $y=(ag-ec)/(ad-bc)$ — each unknown is a ratio of two determinants.""", code=r"""
w, f, k, g = 2.0, 1.0, 1.0, 1.0                 # easy numbers, l = 0
M = np.array([[-1j*w, -f], [f, -1j*w]])         # the unknowns are (u_hat, v_hat)
print("(u_hat, v_hat) =", np.round(np.linalg.solve(M, [-1j*k*g, 0.0]), 3))   # [0.667, -0.333j]: v lags u by a quarter period
""")
remind("C10")
D("D14", ref="13.82")
note("N85 [B]", r"""
**The velocity amplitudes** (step 3 of `D14`): $\hat u=\dfrac{g\hat\eta}{\omega^2-f^2}(\omega k+ifl)$ and
$\hat v=\dfrac{g\hat\eta}{\omega^2-f^2}(-ifk+\omega l)$ (13.80).""")
note("N86 [C]", r"""
**The same relation before $K$ is introduced** (step 5): $\omega^2-f^2=gH(k^2+l^2)$ (13.81). Steps 7–8 (ours) add the group
velocity $\mathbf c_g=c^2\mathbf K/\omega$ and the product $c_pc_g=c^2$.""")
nb.worked_example("a long wave at mid-latitude", r"""
Our 35° N inputs, rounded: $f=8.4\times10^{-5}$ s⁻¹, $c=203$ m/s, wavelength 3140 km ($K=2\times10^{-6}$ m⁻¹).

1. Frequency: $\omega^2=(8.4\times10^{-5})^2+(203\times2\times10^{-6})^2=0.07\times10^{-7}+1.65\times10^{-7}$ → $\omega=4.15\times10^{-4}$ s⁻¹ = 4.9 $f$ (the code below, with unrounded inputs and a 3100 km wave, gives 5.02 $f$).
2. Phase speed $\omega/K=207$ m/s — faster than $c$.
3. Group velocity $c^2K/\omega=4.12\times10^{4}\times2\times10^{-6}/4.15\times10^{-4}=199$ m/s — slower than $c$; product
   $207\times199=4.12\times10^{4}=c^2$.
4. Current ellipse: axis ratio $\omega/f=4.9$, turning clockwise (northern hemisphere, $f>0$; counter-clockwise for $f<0$).
5. A wave ten times longer ($K=2\times10^{-7}$ m⁻¹): $\omega^2=0.71\times10^{-8}+0.16\times10^{-8}$ → $\omega=1.1\,f$ — almost a pure
   inertial circle.""")
code(r"""
w_p = GFD.poincare_omega(K=k31, f=f35, c=c_ext)                     # Eq. (13.82): frequency of a 3100 km wave at 35 N [rad/s]
cgx, cgy = GFD.poincare_group_velocity(k31, 0.0, f35, c_ext)   # group velocity c^2 K / omega [m/s] (ours, D14 step 7)
uh, vh = GFD.poincare_amplitudes(k31, 0.0, w_p, f35, G0, 0.5)  # Eq. (13.80): complex velocity amplitudes for a 0.5 m wave
t_orb = np.linspace(0.0, 2*np.pi/w_p, 181)                     # one wave period [s]
orb = GFD.poincare_orbit(t_orb, k31, 0.5, H_o, f35)            # velocity and particle path at a fixed place
print(f"omega = {w_p:.4e} 1/s = {w_p/f35:.3f} f;  c_p = {w_p/k31:.2f} m/s, c_g = {cgx:.2f} m/s, c_p*c_g/c^2 = {w_p/k31*cgx/c_ext**2:.6f}")   # faster and slower than c
print(f"velocity amplitudes: |u| = {abs(uh):.4f} m/s, |v| = {abs(vh):.4f} m/s, ratio = {abs(uh)/abs(vh):.3f} = omega/f")   # the ellipse
print("sense of rotation at 35 N:", orb["sense"], "| at 35 S:", GFD.poincare_orbit(t_orb, k31, 0.5, H_o, -f35)["sense"])   # mirror in the south
""", explain=r"""
1. `GFD.poincare_omega(K, f, c)` is $\omega^2=f^2+gHK^2$ *(13.82)*: a 3100 km wave on a 4200 m ocean at 35° N oscillates at five times $f$.
2. The phase speed is above $c$ and the group velocity below it; their product is exactly $c^2$.
3. `GFD.poincare_amplitudes` gives the polarisation: the along-wave current is $\omega/f$ times the cross-wave one.
4. `GFD.poincare_orbit` returns the velocity and the particle path over time; the current vector turns clockwise in the
   northern hemisphere and counter-clockwise in the southern.""")
scratch(r"""
# From scratch: the frequency and the real velocity components at x = 0, typed out
w_mine = np.sqrt(f35**2 + G0*H_o*k31**2)                       # Eq. (13.82)
u_mine = w_mine*0.5/(k31*H_o)*np.cos(-w_mine*t_orb)            # Eq. (13.83), first member at x = 0: (omega eta/(kH)) cos(kx - omega t)
v_mine = f35*0.5/(k31*H_o)*np.sin(-w_mine*t_orb)               # second member: (f eta/(kH)) sin(kx - omega t)
assert np.isclose(w_mine, w_p) and np.allclose(u_mine, orb["u"]) and np.allclose(v_mine, orb["v"])   # same as the library
print(f"✓ omega = {w_mine:.4e} 1/s and the velocity ellipse agree with GFD.poincare_orbit (axis ratio {orb['axis_ratio']:.3f})")   # visible confirmation
""", r"""
The dispersion relation and the two real velocity components $u=\frac{\omega\hat\eta}{kH}\cos(kx-\omega t)$, $v=\frac{f\hat\eta}{kH}\sin(kx-\omega t)$ *(13.83)*,
written out, reproduce the library.""")
note("N87 [B]", r"""
**The dispersion diagram** (our version of the book's figure, `ch13.fig_poincare_kelvin_dispersion`), drawn on logarithmic
axes in natural units — frequency in units of $\lvert f\rvert$, wavenumber in units of $1/\Lambda$ with $\Lambda=c/\lvert f\rvert$:
the Poincaré curve $\omega^2=f^2+gHK^2$ *(13.82)*, the Kelvin line $\omega=cK$ (dashed; `C11`), and — **ours**, added for
comparison — the slow Rossby branch.""")
fig(r"""
figP = ch13.fig_poincare_kelvin_dispersion()                   # one log-log panel: Poincare curve, dashed Kelvin line, Rossby branch (35 N, external mode)
plt.show()
""",
    see="Three curves on logarithmic axes: the Poincaré curve (purple), which is flat at $\\omega=\\lvert f\\rvert$ on the left and "
        "rises on the right; the straight dashed Kelvin line $\\omega=cK$; and the Rossby branch (orange), a low arch that "
        "peaks at $K\\Lambda=1$. The thin horizontal line is $\\omega=\\lvert f\\rvert$.",
    read="For long waves ($K\\Lambda\\ll1$) the Poincaré curve flattens onto $\\omega=\\lvert f\\rvert$: no such wave is slower than "
         "$f$. For short waves ($K\\Lambda\\gg1$) it joins the Kelvin line, which is its asymptote $\\omega=cK$ — rotation no "
         "longer matters. The Rossby branch lies below $\\lvert f\\rvert$ everywhere. On logarithmic axes a straight line of "
         "slope 1 means *non-dispersive* (frequency proportional to wavenumber): the Kelvin line everywhere, the Poincaré "
         "curve only on the right.",
    change="…$f\\to0$? Then $\\Lambda\\to\\infty$: every wave is 'short' in these units and the Poincaré curve lies on the "
           "Kelvin line — Chapter 7's non-dispersive long waves.")
note("N88 [B]", r"""
**The real fields** for a wave along $x$ ($l=0$), with $\eta=\hat\eta\cos(kx-\omega t)$:
$u=\dfrac{\omega\hat\eta}{kH}\cos(kx-\omega t)$ and $v=\dfrac{f\hat\eta}{kH}\sin(kx-\omega t)$ (13.83) — stated (they are the real parts
of the amplitudes of `D14`, using $\omega^2-f^2=gHk^2$).""")
code(r"""
eta_f, u_f, v_f = GFD.poincare_fields(np.array([0.0, 7.75e5]), 0.0, 0.0, k31, 0.0, 0.5, H_o, f35)   # fields at x = 0 and a quarter wavelength on, at t = 0
print(f"at the crest:        eta = {eta_f[0]:+.3f} m, u = {u_f[0]:+.4f} m/s, v = {v_f[0]:+.4f} m/s")   # u in phase with eta
print(f"a quarter wave on:   eta = {eta_f[1]:+.3f} m, u = {u_f[1]:+.4f} m/s, v = {v_f[1]:+.4f} m/s")   # v a quarter period out of phase
""", explain=r"""
**What does this show?** Under the crest the water moves with the wave ($u$ in phase with $\eta$); a quarter wavelength away
$u$ and $\eta$ vanish and the cross-wave current $v$ is at its largest — the quarter-period shift that the factor $i$ meant.""")
P("P322", "sense of rotation of (a cos ωt, b sin ωt) from the sign of the swept area", r"""
The point $(a\cos\omega t,\ b\sin\omega t)$ runs round an ellipse. Which way? Its "angular momentum" $x\,dy/dt-y\,dx/dt=ab\omega$ is
positive for counter-clockwise motion. So $(a\cos\omega t, b\sin\omega t)$ with $a,b,\omega>0$ is counter-clockwise — and
$(a\cos(-\omega t), b\sin(-\omega t))$, which is what a wave gives at a fixed place, is clockwise.""", code=r"""
t = np.linspace(0, 1, 5); a, b, w = 2.0, 1.0, -2*np.pi          # the wave's phase at a fixed place is -wt
x, y = a*np.cos(w*t), b*np.sin(w*t)                             # the point at five times
print("sign of the swept area:", np.sign((x[:-1]*np.diff(y) - y[:-1]*np.diff(x)).sum()))   # -1.0: clockwise
""")
note("N89 [B]", r"""
**The velocity ellipse** (our version of the book's figure): the tip of the current vector at a fixed place, with axes
$2\omega\hat\eta/(kH)$ and $2f\hat\eta/(kH)$.""")
fig(r"""
fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.plot(orb["u"]*100, orb["v"]*100, color=COLORS["accent"])    # the velocity hodograph over one period [cm/s]
for ph, lab in ((0.0, "ωt = 0"), (np.pi/2, "ωt = π/2"), (np.pi, "ωt = π")):   # three marked phases
    j_ = np.argmin(abs(w_p*t_orb - ph))                        # the time index of this phase
    ax.plot(orb["u"][j_]*100, orb["v"][j_]*100, "o", color=COLORS["orange"])   # the mark
    ax.annotate(lab, (orb["u"][j_]*100, orb["v"][j_]*100), textcoords="offset points", xytext=(6, 6), fontsize=8)   # its label
ax.set(xlabel="u along the wave [cm/s]", ylabel="v across the wave [cm/s]", aspect="equal", title=f"35° N: {orb['sense']}, axis ratio ω/f = {orb['axis_ratio']:.2f}")
plt.show()
""",
    see="The tip of the current vector over one wave period, with the phases $\\omega t=0$, $\\pi/2$, $\\pi$ marked.",
    read="From $\\omega t=0$ to $\\pi/2$ the vector moves from the positive $u$ axis toward negative $v$: clockwise in the "
         "northern hemisphere ($f>0$; mirror for $f<0$). The long axis lies along the direction of propagation.",
    change="…the wave were longer? Then $\\omega\\to f$ and the ellipse fattens into a circle (animation below).")
trap("T7", r"""
The book's "orbit" figure is a velocity hodograph. The particle's path has the same shape and sense, with both axes divided by
$\omega$. And §13.11 quotes the axis ratio as $\omega/f$ (long over short), §13.14 as $f/\omega$ (short over long) — the same
ellipse.""")
nb.animation(r"""
nfr = 30 if FAST else 40                                       # number of frames
lams = np.geomspace(3.0e5, 3.0e7, nfr)                         # wavelengths from 300 km to 30 000 km [m]
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9))
lines, dots = [], []                                           # the drawn objects of the two panels
for ax, ttl in zip(axes, ("35° N (f > 0)", "35° S (f < 0)")):  # one panel per hemisphere
    (ln,) = ax.plot([], [], color=COLORS["accent"])            # the particle path
    (dt_,) = ax.plot([], [], "o", color=COLORS["orange"])      # the particle now
    ax.set(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2), aspect="equal", xlabel="x / (long semi-axis)", ylabel="y / (long semi-axis)")
    ax.set_title(ttl, fontsize=9)
    lines.append(ln); dots.append(dt_)

fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it: frames render much faster

def update(i):                                                 # frame i: the path for wavelength lams[i]
    k_ = 2*np.pi/lams[i]                                       # wavenumber [rad/m]
    for ax, ln, dt_, f_, hemi in zip(axes, lines, dots, (f35, -f35), ("35° N", "35° S")):   # both hemispheres
        w_ = GFD.poincare_omega(K=k_, f=f_, c=c_ext)                # Eq. (13.82)
        tt = np.linspace(0.0, 2*np.pi/w_, 121)                 # one period
        o_ = GFD.poincare_orbit(tt, k_, 0.5, H_o, f_)          # particle displacement about its mean position
        scale = abs(o_["x"]).max()                             # the long semi-axis, used to normalise the picture
        ln.set_data(o_["x"]/scale, o_["y"]/scale)              # the ellipse
        j_ = (3*i) % 121                                       # the particle advances a little each frame
        dt_.set_data([o_["x"][j_]/scale], [o_["y"][j_]/scale]) # its position now
        ax.set_title(f"{hemi}: λ = {lams[i]/1e3:.0f} km, ω/f = {w_/abs(f_):.2f}, {o_['sense']}", fontsize=8)   # the numbers of this frame
    return lines + dots

show_animation(animate(update, frames=nfr, fig=fig, interval=120), player="video", dpi=50)   # a smooth video
plt.close(fig)                                                 # do not show the last frame a second time
""", explain=r"""
1. The wavelength grows from 300 km to 30 000 km over the frames; for each, `GFD.poincare_orbit` gives the path of a water
   parcel, drawn normalised by its long semi-axis so that only the *shape* changes.
2. The left panel is 35° N, the right panel 35° S; the orange dot moves along the path.""")
nb.figure_notes(
    see="The path of one water parcel under a Poincaré wave as the wavelength grows, in both hemispheres. The title gives the "
        "wavelength, the frequency ratio $\\omega/f$ and the sense of rotation.",
    read="For short waves $\\omega/f$ is large and the path is a thin ellipse — almost the back-and-forth sloshing of a "
         "non-rotating wave. As the wave lengthens, $\\omega/f$ falls toward 1 and the ellipse fattens into the inertial "
         "circle. The dot goes clockwise in the north, counter-clockwise in the south.",
    change="…the layer were the first baroclinic mode ($c=3.6$ m/s)? The same sequence would happen at wavelengths 56 times "
           "shorter.")
note("N90 [B]", r"""
**Inertial motion, the $K\to0$ limit:** $\partial u/\partial t-fv=0$, $\partial v/\partial t+fu=0$, solved by $u=q\cos ft$,
$v=-q\sin ft$: a circle of radius $r=q/f$ traced in one inertial period. Seen from space this is Chapter 4's particle moving in
a straight line while the earth turns under it.""")
code(r"""
q_in = inp["q_inertial"]                                       # speed of the inertial current: 0.23 m/s (our input)
t_in = np.linspace(0.0, GFD.inertial_period(f35), 5)           # five times through one inertial period [s]
u_in, v_in, x_in, y_in = GFD.inertial_oscillation(t_in, q_in, 0.0, f35)   # velocity and position on the inertial circle
print(f"radius = q/f = {GFD.inertial_radius(q_in, f35)/1e3:.2f} km; period = {GFD.inertial_period(f35)/3600:.2f} h at 35 N")   # size and time of the circle
print("speed stays", np.round(np.hypot(u_in, v_in), 3), "m/s; back at the start after one period:", bool(np.hypot(x_in[-1] - x_in[0], y_in[-1] - y_in[0]) < 1.0))   # a closed circle
""", explain=r"""
**What does this show?** A current of 0.23 m/s left to itself at 35° N goes round a circle of radius about 2.8 km once every
21 hours, at constant speed.""")
explainer("shallow_water_dispersion", "Are Poincaré, Kelvin and Rossby waves separate theories?",
          "moving a wavenumber cursor along the diagram while the four term bars of the cubic rebalance shows which terms each "
          "wave is made of; sliding $f$ to zero or $\\beta$ to zero removes a whole branch.",
          ["Select the slow root and read the status: which term is negligible?",
           "Use the preset 'f → 0': the hyperbola collapses onto ω = cK.",
           "Switch β off: the slow root drops to ω = 0.",
           "Toggle 'printed vs corrected' to see slip #7 in the regime text."])
whatif(r"""
…a coast stands in the way, so that the water cannot move across it? A wave with no cross-shore velocity at all becomes
possible — and it can have any frequency, even below $f$ (`C11`).""")


# =====================================================================================================================
# A.12  §13.12 — C11 the Kelvin wave, C12 the Rossby radius and geostrophic adjustment
# =====================================================================================================================
nb.section("13.12", "Kelvin Wave", intro=r"""
**What is this section about?** The wave that leans on a coast: it travels at the ordinary long-wave speed, only in one
direction, and dies away offshore over a distance that turns out to be the most important length scale of rotating fluids.""")
core("C11", r"The Kelvin wave: $\eta=\eta_0e^{-fy/c}\cos k(x-ct)$, $u=\eta_0\sqrt{g/H}\,e^{-fy/c}\cos k(x-ct)$ (13.87)",
     "How can a wave exist below the frequency f, and why does it run only one way along a coast?")
problem(r"""
The tide in the North Sea does not slosh in and out; it travels round the basin, counter-clockwise, as a wave whose range is
largest at the shore. In the Pacific, a relaxation of the trade winds sends a signal along the equator and then up and down
the coast of the Americas. Both are Kelvin waves: gravity waves that use a boundary to dodge the Coriolis force.""")
note("N91 [B]", r"""
**The idea.** At a wall the water cannot move across the shore, so $v=0$ there — try $v=0$ everywhere. Then the Coriolis force
on the along-shore current has nothing to turn; it must be *held* by a pressure force: the surface tilts across the shore,
$fu=-g\dfrac{\partial\eta}{\partial y}$. Under a crest the current flows with the wave, so the surface must be high at the wall
and fall offshore — and that fixes the direction of travel: coast on the right (northern hemisphere, $f>0$; on the left for
$f<0$).

```
   cross-shore section through a crest, looking ALONG the direction of travel, f > 0 (coast on the RIGHT)

                              .~~~~~~|▓ coast     the surface is highest at the wall
   ______.....~~~~~~~                |▓           and falls offshore (to the left) over the distance c/f
        ⊗ current under the crest    |▓           (into the page: with the wave)
        Coriolis force →             |▓           pushes the water toward the coast (to the right of the current);
        ← pressure force of the slope             the surface slope pushes it back

   for f < 0 mirror the picture: coast on the LEFT, Coriolis force to the left of the current
```""")
note("N94 [B]", r"""
**The equations with $v\equiv0$:** $\dfrac{\partial\eta}{\partial t}+H\dfrac{\partial u}{\partial x}=0$, $\dfrac{\partial u}{\partial t}=-g\dfrac{\partial\eta}{\partial x}$,
$fu=-g\dfrac{\partial\eta}{\partial y}$ (13.84).""")
note("N93 [C]", r"""
**Compare the two waves** through the cross-shore equation $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$:
for a Poincaré wave the Coriolis term is partly balanced by the acceleration $\partial v/\partial t$; for a Kelvin wave $v=0$ and
it is balanced entirely by the slope — geostrophically.""")
P("P323", "trapped solutions: keeping the exponential that decays away from a boundary", r"""
A first- or second-order equation in the offshore coordinate often has two exponential solutions, one growing and one decaying
away from the wall. In a half-space only the decaying one is physical (finite energy). Which one decays depends on signs in
the problem — here on the sign of $f$ and on the direction of travel.""", code=r"""
y = np.array([0.0, 1.0, 2.0])                  # distance from the wall in units of c/|f|
for s in (+1, -1):                             # s = sign(f) * direction of travel
    print("s =", s, "-> profile offshore:", np.round(np.exp(-s*y), 3))   # +1 decays offshore (kept), -1 grows (rejected)
""")
remind("C11")
D("D15", ref="13.87")
note("N95 [B]", r"""
**Steps 2–4 of `D15`:** with $[u,\eta]=[\hat u(y),\hat\eta(y)]e^{i(kx-\omega t)}$, the three equations become
$-i\omega\hat\eta+iHk\hat u=0$, $-i\omega\hat u=-igk\hat\eta$, $f\hat u=-g\dfrac{d\hat\eta}{dy}$ (13.85), hence
$\hat\eta\,[\omega^2-gHk^2]=0$.""")
nb.recap("R23", "The Kelvin wave is not dispersive", r"""
Its speed is $c=\sqrt{gH}$ (13.86), i.e. $\omega=\pm k\sqrt{gH}$: the non-rotating long-wave speed, for every frequency — including below
$f$, where no Poincaré wave exists.""", where="Ch. 7, long-wave speed")
nb.current_core = "C11"
nb.worked_example("a Kelvin wave on a deep ocean", r"""
Our inputs, rounded: $H=4200$ m, $g\approx9.8$ m/s², $f=8.4\times10^{-5}$ s⁻¹ (35° N), crest height $\eta_0=1$ m at the coast.

1. Speed: $c=\sqrt{gH}=\sqrt{41\,160}\approx203$ m/s.
2. Trapping width $\Lambda=c/f=203/8.4\times10^{-5}=2.42\times10^{6}$ m ≈ 2400 km.
3. Current under the crest at the coast: $u=\eta_0\sqrt{g/H}=1\times\sqrt{9.8/4200}=0.048$ m/s.
4. At $y=\Lambda$ offshore both are $e^{-1}$ = 37 % of that.
5. Check the cross-shore balance at the coast: $fu=8.4\times10^{-5}\times0.048=4.1\times10^{-6}$ m/s²;
   $-g\,\partial\eta/\partial y=g\eta_0/\Lambda=9.8/2.42\times10^{6}=4.1\times10^{-6}$ m/s² ✓.""")
code(r"""
eta_k, u_k = GFD.kelvin_wave(0.0, 0.0, 0.0, 0.5, k31, H_o, f35)   # Eq. (13.87) at the coast under the crest: 0.5 m wave, 3100 km long, 35 N
Lam_ext = GFD.rossby_radius(c_ext, f35)                        # the trapping width c/|f| [m]
print(f"trapping width c/f = {Lam_ext/1e3:.0f} km; current under the crest at the coast u = {u_k:.4f} m/s; period = {2*np.pi/GFD.kelvin_omega(k31, c_ext)/3600:.2f} h")   # scales of the wave
print("35 N, wave toward +x (fluid in y > 0):", GFD.kelvin_decay_side(f35, +1))    # trapped, coast on the right
print("35 S, wave toward +x (fluid in y > 0):", GFD.kelvin_decay_side(-f35, +1))   # not trapped: it would have to run the other way
try:                                                           # ask for the wave that runs the wrong way at 35 N
    GFD.kelvin_wave(0.0, 0.0, 0.0, 0.5, k31, H_o, f35, direction=-1)
except ValueError as err:                                      # the function refuses: that solution grows offshore
    print("wrong direction ->", err)
""", explain=r"""
1. `GFD.kelvin_wave(x, y, t, eta0, k, H, f)` returns the surface height and the along-shore current, with the fluid in $y\ge0$
   and the coast at $y=0$; by default it travels in the trapped direction.
2. For the external mode at 35° N the trapping width is about 2400 km and a 0.5 m wave carries a current of 2.4 cm/s.
3. `GFD.kelvin_decay_side(f, direction)` says whether a wave in that direction is trapped and on which side the coast must be.
4. `try: … except ValueError:` runs a call that may fail and catches the failure: asking for the wave that runs with the coast
   on its left at 35° N raises an error, because that solution grows offshore without bound.""")
scratch(r"""
# From scratch: put the formula into the three equations (13.84) by finite differences
x0, y0, t0, hx, ht = 2.0e5, 5.0e5, 3.0e3, 10.0, 1.0            # a point offshore, a time, and small steps [m, m, s, m, s]
Kw = lambda xx, yy, tt: GFD.kelvin_wave(xx, yy, tt, 0.5, k31, H_o, f35)   # (eta, u) at any point and time
deta_dt = (Kw(x0, y0, t0 + ht)[0] - Kw(x0, y0, t0 - ht)[0])/(2*ht)       # d(eta)/dt by a centred difference
deta_dx = (Kw(x0 + hx, y0, t0)[0] - Kw(x0 - hx, y0, t0)[0])/(2*hx)       # d(eta)/dx
deta_dy = (Kw(x0, y0 + hx, t0)[0] - Kw(x0, y0 - hx, t0)[0])/(2*hx)       # d(eta)/dy
du_dt = (Kw(x0, y0, t0 + ht)[1] - Kw(x0, y0, t0 - ht)[1])/(2*ht)         # du/dt
du_dx = (Kw(x0 + hx, y0, t0)[1] - Kw(x0 - hx, y0, t0)[1])/(2*hx)         # du/dx
u0 = Kw(x0, y0, t0)[1]                                                    # u itself
res = [abs(deta_dt + H_o*du_dx)/abs(deta_dt), abs(du_dt + G0*deta_dx)/abs(du_dt), abs(f35*u0 + G0*deta_dy)/abs(f35*u0)]   # each divided by its largest term
assert max(res) < 1e-5                                         # all three equations are satisfied
lib = GFD.kelvin_residuals(x0, y0, t0, 0.5, k31, H_o, f35)     # the library's own check
print(f"✓ residuals by hand: {[f'{r_:.0e}' for r_ in res]}; library: {[f'{abs(v_):.0e}' for v_ in lib.values()]} (continuity, x-momentum, cross-shore geostrophy)")   # visible confirmation
""", r"""
Differentiating the closed form numerically and putting the result into the three equations with $v\equiv0$ leaves relative
residuals at the level of the finite-difference error; `GFD.kelvin_residuals` does the same check inside the library.""")
note("N92 [B]", r"""
**Cross-shore sections and a plan view** (our version of the book's figures, `ch13.fig_kelvin_sections`): the surface across the
shore through a crest and through a trough, and a map of one coast with the surface height in colour and the current as
arrows.""")
fig(r"""
figK = ch13.fig_kelvin_sections()                              # cross-shore sections of eta under a crest and a trough, and a plan view of one coast (our inputs, 35 N)
plt.show()
""",
    see="One panel: the sea surface across the shore under a crest and under a trough — highest (or lowest) at the wall, "
        "relaxing exponentially offshore. Other panel: a map of the coast, colours for the surface height and arrows for the "
        "current, both fading within the distance $\\Lambda=c/f$ of the wall.",
    read="The slope across the shore reverses between crest and trough — and so does the along-shore current, so the Coriolis "
         "force is always held by the slope. In the map the current runs with the wave under crests and against it under "
         "troughs; the wave travels with the coast on its right (northern hemisphere, $f>0$; mirror for $f<0$).",
    change="…there were a second coast less than a few $\\Lambda$ away (a channel)? That wall would carry its own Kelvin wave, "
           "travelling the opposite way; in a channel much narrower than $c/f$ the two overlap and rotation hardly matters.")
nb.plotly(r"""
xs_k = np.linspace(0.0, 6.2e6, 60); ys_k = np.linspace(0.0, 3*Lam_ext, 40)   # two wavelengths alongshore, three trapping widths offshore [m]
Xk, Yk = np.meshgrid(xs_k, ys_k)                               # grids indexed [j, i] = (y, x)
figKs = go.Figure()
for j, (f_, nm) in enumerate(((f35, "35° N"), (-f35, "35° S"))):   # the trapped wave of each hemisphere
    eta_s, _ = GFD.kelvin_wave(Xk, Yk, 0.0, 0.5, k31, H_o, f_) # Eq. (13.87); the default direction is the trapped one
    figKs.add_trace(go.Surface(x=xs_k/1e3, y=ys_k/1e3, z=eta_s, colorscale="RdBu_r", cmin=-0.5, cmax=0.5, visible=(j == 0),
                               colorbar=dict(title="η [m]"), name=nm))   # the sea surface
figKs.update_layout(height=520, title="The Kelvin wave hugs the coast (the wall is the edge y = 0)",
                    updatemenus=[dict(x=0.0, y=1.08, xanchor="left", buttons=[
                        dict(label="35° N: travels toward +x (coast on its right)", method="update", args=[{"visible": [True, False]}]),
                        dict(label="35° S: travels toward −x (coast on its left)", method="update", args=[{"visible": [False, True]}])])],
                    scene=dict(xaxis_title="x alongshore [km]", yaxis_title="y offshore [km]", zaxis_title="η [m]",
                               aspectratio=dict(x=1.6, y=1.0, z=0.4)), margin=dict(l=0, r=0, t=60, b=0))   # axes and the hemisphere dropdown
figKs.show()
""", explain=r"""
The surface $\eta(x,y)$ of the trapped wave at one instant, for either hemisphere. In both the fluid is on the side $y>0$ of
the wall; what changes is the direction of travel.""")
nb.figure_notes(
    see="A wave whose crests and troughs are tallest along the edge $y=0$ (the coast) and fade offshore.",
    read="The shape is the same in both hemispheres; the direction of travel is not. With the fluid on the same side of the "
         "wall, the northern wave runs toward $+x$ (coast on its right) and the southern one toward $-x$ (coast on its left).",
    change="…the wave were on the first baroclinic mode? The offshore scale would shrink from 2400 km to 43 km.")
choice(r"""
the animation below is **our numerical model, cached**: a linear C-grid run (`SW.ShallowWater`, 64 × 64, closed basin) whose
equivalent depth was chosen so that the trapping width is a fifth of the basin. It is loaded with
`ch13.load_reference_run("kelvin_basin")`; if the file were missing, the straight-coast closed form would be animated instead.""")
nb.animation(r"""
kb = ch13.load_reference_run("kelvin_basin")                   # our cached run: eta(t, y, x) in a closed square basin, or None
fig, ax = plt.subplots(figsize=(4.8, 4.5))
if kb is not None:                                             # the cached model output
    eta_fr = np.asarray(kb["eta"], dtype=float)                # 24 frames of the surface height [m]
    ext = [kb["x"][0]/1e3, kb["x"][-1]/1e3, kb["y"][0]/1e3, kb["y"][-1]/1e3]   # basin extent [km]
    src = f"cached C-grid run, Λ = {float(kb['Lambda'])/1e3:.0f} km"            # caption
    t_fr = np.asarray(kb["t"], dtype=float)                    # times of the frames [s]
else:                                                          # fallback: the closed form along a straight coast
    t_fr = np.linspace(0.0, 2*np.pi/GFD.kelvin_omega(k31, c_ext), 24)           # one period
    eta_fr = np.array([GFD.kelvin_wave(Xk, Yk, t_, 0.5, k31, H_o, f35)[0] for t_ in t_fr])   # 24 frames
    ext = [0, xs_k[-1]/1e3, 0, ys_k[-1]/1e3]                   # extent [km]
    src = "closed form, straight coast (cache absent)"         # caption
vmax = abs(eta_fr).max()                                       # symmetric colour limits
im = ax.imshow(eta_fr[0], origin="lower", extent=ext, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")   # the first frame
fig.colorbar(im, label="surface height η [m]")
ax.set(xlabel="x [km]", ylabel="y [km]")

fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it: frames render much faster

def update(i):                                                 # frame i: the surface at time t_fr[i]
    im.set_data(eta_fr[i])                                     # new field
    ax.set_title(f"{src}\nt = {t_fr[i]/86400:.2f} days", fontsize=9)   # the time of this frame
    return (im,)

show_animation(animate(update, frames=len(eta_fr), fig=fig, interval=250), player="video", dpi=50)   # a smooth video
plt.close(fig)                                                 # do not show the last frame a second time
""", explain=r"""
1. `ch13.load_reference_run("kelvin_basin")` returns the saved arrays of our model run (or `None` if the file is absent, in
   which case the closed form along a straight coast is used and the caption says so).
2. `imshow` draws the first frame; `update(i)` swaps in frame `i`.""")
nb.figure_notes(
    see="A bulge of the sea surface travelling round a closed square basin, pressed against the wall.",
    read="The bulge goes round counter-clockwise ($f>0$) and turns the corners: at every wall the coast is on the right of the "
         "direction of travel. Away from the walls the surface hardly moves.",
    change="…$f<0$? The bulge would go round clockwise, with the coast on its left.")
note("N96 [B]", r"""
**Upwelling sends a Kelvin wave** (joining `C05` and `C11`). An equatorward wind along an eastern ocean boundary drives Ekman
transport offshore; deeper water rises; the lifted thermocline is a disturbance that leaves **poleward** along the coast as an
internal Kelvin wave (coast on its right in the northern hemisphere). Its speed and width, with our two-layer inputs, are
computed in `C12`.""")
code(r"""
cu2 = ch13.coastal_upwelling(-0.07, "east", lat)               # equatorward wind stress of 0.07 N/m^2 along an eastern boundary at 35 N
print(f"offshore Ekman transport {cu2['transport_offshore']:.3f} m^2/s -> upwelling: {cu2['upwelling']}; the thermocline signal leaves {cu2['kelvin_direction']} as an internal Kelvin wave")   # the chain of events
""", explain=r"""
**What does this show?** One function strings the two blocks together: the wind moves the surface water offshore, the coast
upwells, and the disturbance propagates poleward.""")
explainer("kelvin_wave", "Why does this wave run only one way along a coast?",
          "reversing the direction of travel (or the hemisphere) makes the offshore profile grow instead of decay and the status "
          "badge reject it; widening the channel separates the two coastal waves; switching to the internal mode shrinks the "
          "trapping width from thousands of kilometres to tens.",
          ["Reverse the direction: the status says 'would grow offshore — not a solution'.",
           "Switch to S and see which way the trapped wave now runs.",
           "Choose the channel 2Λ wide and watch the two walls carry opposite waves (a case the static figure above does not show).",
           "Switch to the internal mode and read the new Λ."])
whatif(r"""
…there is no coast? The width $c/f$ still means something: it is how far a gravity wave gets before rotation takes over
(`C12`).""")

# ---------------------------------------------------------------------------------------------------------------------
core("C12", r"The Rossby radius of deformation $\Lambda\equiv c/f$",
     "If I pile water up and let go, why doesn't it flatten out as it would without rotation — and how wide is what stays?")
problem(r"""
Pour a bucket of water into a bath and the bump spreads until the surface is flat. Do the same on a planet-sized turntable and
it does not: the water starts to spread, the Coriolis force turns the outflow into a current running *along* the edge of the
bump, and that current's own Coriolis force holds the rest of the bump up. The ocean is full of such fronts and eddies standing
in slopes that gravity cannot flatten. Their width is the Rossby radius — tens of kilometres in the ocean, about a thousand in
the atmosphere — and it is why an ocean model needs a grid of a few kilometres to "resolve eddies".""")
idea(words=r"""
In a time $1/f$ a gravity wave travels a distance $c/f$. Disturbances narrower than that spread before rotation notices them;
wider ones are caught by rotation before gravity can flatten them.

| Size of the disturbance | What happens |
|---|---|
| much smaller than $\Lambda$ | behaves as if there were no rotation: it radiates away as gravity waves |
| much larger than $\Lambda$ | stays, in geostrophic balance |""")
nb.recap("R24", "Internal Kelvin waves and the internal radius", r"""
On the interface between a thin upper layer and a deep lower one the long-wave speed is $c=\sqrt{g'H}$ with the reduced
gravity $g'=g(\rho_2-\rho_1)/\rho_2$; in a continuously stratified layer the speeds are $c=NH/n\pi$ (`C08`). The corresponding
radius $\Lambda=NH/(\pi f)$ for $n=1$ is far smaller than the external one; the interface moves far more than the surface, in
the opposite sense.""", where=r"Ch. 7, $g'=g(\rho_2-\rho_1)/\rho_2$ (7.117)")
nb.current_core = "C12"
trap("T11", r"""
Three Rossby radii share one name: external $\sqrt{gH}/f$; internal $\sqrt{g'H}/f$ or $NH/(n\pi f)$; and in `C16` the Eady
radius $NH/f$ **without** the $\pi$.""")
loose(11, r"a typical internal radius about three times what its own typical $N$, $H$ and $f$ give",
      r"$\Lambda=NH/(\pi f)$ evaluated with the inputs one states — the mismatch is between the quoted inputs and the quoted radius, not in the formula.")
remind("C12")
nb.worked_example("three radii by hand", r"""
Our 35° N value, rounded: $f=8.4\times10^{-5}$ s⁻¹.

1. **External**, $H=4200$ m: $c=\sqrt{9.8\times4200}\approx203$ m/s → $\Lambda\approx2400$ km.
2. **Two-layer**, upper layer 120 m, $\Delta\rho/\rho=0.003$: $g'=0.029$ m/s², $c=\sqrt{0.029\times120}\approx1.9$ m/s → $\Lambda\approx22$ km.
3. **Uniform** $N=2.7\times10^{-3}$ s⁻¹ over 4200 m: $c_1=NH/\pi=11.3/\pi\approx3.6$ m/s → $\Lambda\approx43$ km.
4. Ratio external / internal: 56 to 110 — one ocean, two utterly different scales.""")
code(r"""
rho2 = inp["rho_ocean"]; rho1 = rho2 - inp["drho"]             # lower-layer and upper-layer density [kg/m^3] (difference 3.1)
L_ext = GFD.rossby_radius(GFD.long_wave_speed(H_o), f35)                       # external: sqrt(gH)/f
L_int = GFD.rossby_radius_internal(N_o, H_o, f35)                              # internal, n = 1: N H/(pi f)
L_eady = GFD.rossby_radius_internal(N_o, H_o, f35, with_pi=False)              # N H/f, no pi (the Eady radius of this ocean)
L_two = GFD.rossby_radius_two_layer(inp["H1"], rho1, rho2, f35)                # two-layer: sqrt(g' H1)/f for a 120 m upper layer
c_two = L_two*f35                                                               # the internal long-wave speed sqrt(g' H1) [m/s]
print(f"35 N: external {L_ext/1e3:.0f} km | internal N H/(pi f) {L_int/1e3:.2f} km | without pi {L_eady/1e3:.1f} km | two-layer {L_two/1e3:.1f} km")   # four numbers, one name
print(f"two-layer: g' = {G0*inp['drho']/rho2:.4f} m/s^2, c = {c_two:.3f} m/s = {c_two*86.4:.0f} km per day")   # the internal Kelvin wave of note N96
print(f"12 N: internal N H/(pi f) = {GFD.rossby_radius_internal(N_o, H_o, f12)/1e3:.1f} km")                     # larger toward the equator
""", explain=r"""
1. `GFD.rossby_radius(c, f)` is $\Lambda=c/\lvert f\rvert$ for any wave speed.
2. `GFD.rossby_radius_internal` is $NH/(n\pi f)$; with `with_pi=False` it returns $NH/f$, the radius used in the Eady problem.
3. `GFD.rossby_radius_two_layer` uses the reduced gravity of a two-layer ocean.
4. External and internal radii differ by a factor of about 56 at the same latitude; and every radius grows toward the equator,
   where $f$ is small.""")
P("P324", "an adjustment problem: what a steady end state can remember (a conserved quantity pins it)", r"""
Release an unbalanced state and wait. Waves carry away whatever can travel. What is left must be steady — but there are
infinitely many steady states. The one nature picks is singled out by a quantity that cannot change at any point while the
waves pass. Find that quantity, demand it has its initial value, and the end state follows without solving for the transient
at all.""", code=r"""
# a bank analogy: transfers move money between branches (the waves); the total cannot change
start = np.array([10.0, 0.0, 0.0]); end = np.full(3, start.sum()/3)   # all in one branch; then spread evenly
print(end, end.sum() == start.sum())                                  # the end state is fixed by what is conserved
""")
D("D16")
note("N19 [B]", r"""
**Not in the book — ours** (the book gives only a verbal account of how geostrophy is set up, restated in `C02`): the
adjustment of a step $\eta_0\,\mathrm{sgn}(x)$ to the steady state $\eta=\eta_0\,\mathrm{sgn}(x)\,(1-e^{-\lvert x\rvert/\Lambda})$ with a
jet $v=(g\eta_0/c)\,e^{-\lvert x\rvert/\Lambda}$ along the step (`D16`). Label: analytic (ours); no source is cited for it.""")
code(r"""
x_a = (np.arange(800) + 0.5)*5.0e3 - 2.0e6                     # 800 cell centres, 5 km apart, from -2000 km to +2000 km [m]
eta_end_a, v_end_a = GFD.geostrophic_adjustment_1d(x_a, 0.05, He1, f35)   # closed-form end state of a +/-5 cm step (first baroclinic mode as a layer)
Lam1 = GFD.rossby_radius(GFD.long_wave_speed(He1), f35)        # its Rossby radius [m]
en = GFD.adjustment_energy(0.05, He1, f35)                     # energy budget per unit length of the step [J/m]
_, v_south = GFD.geostrophic_adjustment_1d(x_a, 0.05, He1, -f35)   # the same step at 35 S
print(f"Rossby radius = {Lam1/1e3:.2f} km; jet maximum g eta0/c = {G0*0.05/GFD.long_wave_speed(He1):.3f} m/s (along +y at 35 N, {'-y' if v_south[400] < 0 else '+y'} at 35 S)")   # width and strength of what stays
print(f"energy: released PE {en['pe_released']:.2e}, jet KE {en['ke_jet']:.2e}, radiated {en['radiated']:.2e} J/m; KE/PE = {en['ratio']:.4f}")   # one third stays
""", explain=r"""
1. `GFD.geostrophic_adjustment_1d(x, eta0, H, f)` is the closed-form end state of `D16`: a front of width $\Lambda$ and a jet
   along it.
2. For the first baroclinic mode at 35° N the front is 43 km wide and the jet reaches about 14 cm/s; it flows the other way in
   the southern hemisphere.
3. `GFD.adjustment_energy` gives the budget: of the potential energy released, exactly one third stays as kinetic energy of
   the jet; two thirds leave with the waves.""")
scratch(r"""
# From scratch: the three radii, the adjusted step, and the energy ratio by numerical integration
assert np.isclose(np.sqrt(G0*H_o)/f35, L_ext)                  # external radius sqrt(gH)/f
assert np.isclose(np.sqrt(G0*(rho2 - rho1)/rho2*inp["H1"])/f35, L_two)   # two-layer radius sqrt(g' H1)/f
assert np.isclose(N_o*H_o/(np.pi*f35), L_int)                  # uniform-N radius N H/(pi f)
eta_mine = 0.05*np.sign(x_a)*(1 - np.exp(-abs(x_a)/Lam1))      # D16: the adjusted surface
v_mine = (G0*0.05/np.sqrt(G0*He1))*np.exp(-abs(x_a)/Lam1)      # D16: the jet (f > 0)
assert np.allclose(eta_mine, eta_end_a) and np.allclose(v_mine, v_end_a)   # same as the library
xx = np.linspace(-20*Lam1, 20*Lam1, 40001)                     # a fine grid over +/- 20 Rossby radii
e_x = 0.05*np.sign(xx)*(1 - np.exp(-abs(xx)/Lam1)); v_x = (G0*0.05/np.sqrt(G0*He1))*np.exp(-abs(xx)/Lam1)   # end state on it
KE = 0.5*1000.0*He1*np.trapezoid(v_x**2, xx)                   # kinetic energy of the jet per unit length [J/m]
PE = 0.5*1000.0*G0*np.trapezoid(0.05**2 - e_x**2, xx)          # potential energy released: initial minus final [J/m]
assert np.isclose(KE/PE, en["ratio"], rtol=1e-3)               # one third
print(f"✓ three radii, the end state and the energy ratio KE/PE = {KE/PE:.4f} agree with the library (1/3)")   # visible confirmation
""", r"""
The three radii are one-line formulas; the end state is two exponentials; and integrating $\tfrac12\rho Hv^2$ and
$\tfrac12\rho g(\eta_0^2-\eta^2)$ numerically gives the ratio one third.""")
nb.animation(r"""
dx_a, dt_a = 5.0e3, 600.0                                      # 5 km cells, 600 s steps (Courant number 0.43)
T_in = GFD.inertial_period(f35)                                # inertial period at 35 N [s]
n_per = int(round(T_in/dt_a))                                  # time steps per inertial period
step0 = 0.05*np.sign(x_a)                                      # the initial step: +5 cm on the right, -5 cm on the left
runs_a = [SW.linear_1d_run(step0, dx=dx_a, dt=dt_a, n_steps=6*n_per, H=He1, f=f_) for f_ in (0.0, f35)]   # six inertial periods, without and with rotation
mean_eta = runs_a[1]["eta"][5*n_per:6*n_per + 1].mean(axis=0)  # the MEAN surface over the sixth inertial period (rotating run)
mean_v = runs_a[1]["v"][5*n_per:6*n_per + 1].mean(axis=0)      # and the mean along-step velocity
stops = [0, 1, 3, 5]                                           # the snapshots: after 0, 1, 3, 5 inertial periods
fig, axes = plt.subplots(2, 1, figsize=(6.6, 5.2), sharex=True)
ln0, = axes[0].plot(x_a/1e3, step0*100, color=COLORS["accent"])   # top row: no rotation
ln1, = axes[1].plot(x_a/1e3, step0*100, color=COLORS["accent"], label="our model")   # bottom row: with rotation
axes[1].plot(x_a/1e3, eta_end_a*100, "--", color=COLORS["ink"], lw=1, label="closed-form end state (D16)")   # the analytic ghost
for s_ in (-1, 1):
    axes[1].axvline(s_*Lam1/1e3, color=COLORS["muted"], lw=0.7)   # mark x = +/- Lambda
axes[0].set(ylim=(-8, 8), ylabel="η [cm]"); axes[1].set(ylim=(-8, 8), xlim=(-600, 600), xlabel="x [km]", ylabel="η [cm]")
axes[1].legend(fontsize=7, loc="upper left")

fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it: frames render much faster

def update(i):                                                 # frames 0-3: snapshots; frame 4: the period mean
    if i < len(stops):
        ln0.set_ydata(runs_a[0]["eta"][stops[i]*n_per]*100)    # snapshot without rotation
        ln1.set_ydata(runs_a[1]["eta"][stops[i]*n_per]*100)    # snapshot with rotation
        axes[0].set_title(f"no rotation — snapshot after {stops[i]} inertial periods: the fronts run off, the middle goes flat", fontsize=8)
        axes[1].set_title("f at 35° N — snapshot (waves and grid-scale ripples still present)", fontsize=8)
    else:
        ln0.set_ydata(runs_a[0]["eta"][-1]*100)                # without rotation: the last snapshot
        ln1.set_ydata(mean_eta*100)                            # with rotation: the mean over the sixth inertial period
        axes[0].set_title("no rotation — after 6 inertial periods", fontsize=8)
        axes[1].set_title("f at 35° N — MEAN over the sixth inertial period: a front of width Λ stays", fontsize=8)
    return ln0, ln1

show_animation(animate(update, frames=len(stops) + 1, fig=fig, interval=1200), player="frames", dpi=42)   # step with the buttons
plt.close(fig)                                                 # do not show the last frame a second time
near = abs(x_a) < 3*Lam1                                       # the region within three Rossby radii of the step
print(f"mean over the sixth inertial period vs closed form, within 3 Lambda: eta differs by {100*abs(mean_eta - eta_end_a)[near].max()/0.05:.1f} % of eta0, v by {100*abs(mean_v - v_end_a)[near].max()/abs(v_end_a).max():.1f} % of the jet maximum")   # measured agreement
for s_ in (1, 3, 5):                                           # how far single snapshots still are from the end state
    print(f"snapshot after {s_} periods: eta differs by up to {100*abs(runs_a[1]['eta'][s_*n_per] - eta_end_a)[near].max()/0.05:.0f} % of eta0")   # waves and ripples
""", explain=r"""
1. `SW.linear_1d_run` marches the step for six inertial periods, once with $f=0$ and once with $f$ at 35° N.
2. Frames 0–3 are **snapshots** after 0, 1, 3 and 5 inertial periods. A snapshot of the rotating run still carries
   near-inertial waves, which die away only slowly, and — because the step is sharp and our scheme has no damping — grid-scale
   ripples.
3. The last frame therefore shows the **mean over the sixth inertial period**, which is what should be compared with the
   steady end state (dashed). The printed lines give the measured agreement of the mean, and how far single snapshots are.
4. `player="frames"` gives step buttons, so that you can stop at each stage.""")
nb.figure_notes(
    see="Top: without rotation the step splits into two fronts that run off at $\\pm c$; between them the surface is flat. "
        "Bottom: with rotation, waves leave too, but a front of width $\\Lambda$ (grey lines) remains.",
    read="Step to the last frame: the period-mean of the model lies on the dashed closed form of `D16`. Rotation has held up "
         "everything beyond about one Rossby radius from the step.",
    change="…the step were replaced by a bump much narrower than $\\Lambda$? Almost all of it would radiate away (try it in the "
           "explainer below).")
fig(r"""
x_face = 0.5*(x_a[1:] + x_a[:-1])                              # the interior cell faces, where the potential vorticity lives [m]
fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, 3.6))
a.plot(x_face/1e3, runs_a[1]["pv"][0]*1e6, color=COLORS["amber"], lw=3, label="at the start")       # linear PV at the first step (amber = vorticity)
a.plot(x_face/1e3, runs_a[1]["pv"][1]*1e6, "--", color=COLORS["ink"], lw=1, label="after 6 inertial periods")   # and at the last
a.set(xlim=(-300, 300), xlabel="x [km]", ylabel=r"$\partial v/\partial x-f\eta/H$ [10$^{-6}$ s$^{-1}$]", title="The quantity that pinned the end state")
a.legend(fontsize=8)
b.bar(["PE released", "KE of the jet", "radiated"], [en["pe_released"], en["ke_jet"], en["radiated"]], color=[COLORS["blue"], COLORS["accent"], COLORS["muted"]])   # the energy budget [J/m]
b.set(ylabel="energy per unit length [J m$^{-1}$]", title="One third stays, two thirds leave")
plt.show()
print(f"largest change of the linear potential vorticity during the run: {abs(runs_a[1]['pv'][1] - runs_a[1]['pv'][0]).max():.1e} 1/s")   # round-off
""",
    see="Left: the linear potential vorticity $\\partial v/\\partial x-f\\eta/H$ against $x$ at the first and at the last time "
        "step of the rotating run. Right: the energy budget of the adjustment.",
    read="The two curves on the left coincide: while the waves rearranged the surface and the velocity, this combination did "
         "not change anywhere. It is what the end state \"remembers\". `C13` makes it exact and nonlinear.",
    change="…$f$ were halved? The radius $\\Lambda$ doubles: twice as much potential energy is released and the jet is twice as wide; the "
           "ratio stays one third.")
explainer("geostrophic_adjustment", "Why doesn't a pile of water flatten out on a rotating planet?",
          "pressing play with $f=0$ (everything radiates, nothing stays) and then with rotation (waves leave, a front of width "
          "$\\Lambda$ stays) is the whole idea; the linked bars show how much potential energy was released and how much stayed.",
          ["Play the preset 'no rotation', then the mid-latitude preset, and compare the end cards.",
           "Halve f: the front is twice as wide.",
           "Switch from the step to the narrow bump: what fraction survives?",
           "Watch the potential-vorticity curve: it does not move."])
whatif(r"""
…the motion is not small and $f$ varies with latitude? The quantity that pinned the end state here survives all of that:
potential vorticity (`C13`).""")


# =====================================================================================================================
# A.13  §13.13 — C13 potential vorticity
# =====================================================================================================================
nb.section("13.13", "Potential Vorticity Conservation in Shallow-Water Theory", intro=r"""
**What is this section about?** The one conservation law behind all slow, large-scale motion: every column of fluid carries
the ratio of its absolute vorticity to its depth unchanged. Stretch it, squash it or move it to another latitude, and its
spin must respond.""")
core("C13", r"Potential vorticity: $\dfrac{D(\zeta+f)}{Dt}=\dfrac{\zeta+f_0}{h}\dfrac{Dh}{Dt}$ (13.93), hence $\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=0$ with $f=f_0+\beta y$ (13.94)",
     "What does a column of fluid keep as it moves across an ocean of changing depth and latitude?")
problem(r"""
A skater who pulls her arms in spins faster. A column of ocean that is stretched taller gets thinner (its volume is fixed) and
spins faster too — except that the column already carries a spin it never chose: the earth's, $f$, which depends on latitude.
So there are two ways to change a column's own spin $\zeta$: change its height, or carry it north or south. Westerlies
crossing the Rockies are squashed, turn toward the equator, overshoot and meander downstream for thousands of kilometres.
Dynamicists treat potential vorticity as the dye that marks large-scale air and water masses.""")
idea(words=r"""
| What happens to the column | Northern hemisphere ($f>0$) | Southern hemisphere ($f<0$) |
|---|---|---|
| stretched ($h$ grows) | $\zeta+f$ grows: more counter-clockwise spin | $\zeta+f$ becomes more negative: more clockwise spin |
| carried poleward ($\lvert f\rvert$ grows) at fixed $h$ | $\zeta$ falls (clockwise relative spin) | $\zeta$ rises (counter-clockwise relative spin) |

In both hemispheres stretching makes the absolute vorticity larger in magnitude, and a poleward trip gives anticyclonic
relative spin.""")
note("N104 [C]", r"""
**The geometry** (our sketch): a column of total depth $h$ (over an uneven bottom, $h$ changes as the column moves) drawn
three times with the same volume — as it starts, stretched tall and thin, and squashed short and wide.""")
fig(r"""
figC = draw_pv_column()                                        # our sketch: the same column, tall and thin or short and fat
plt.show()
""",
    see="Three columns of the same volume: the reference column ($h_0$, $\\zeta=0$), the same column stretched (taller, "
        "thinner, with a counter-clockwise arrow: $\\zeta>0$) and squashed (shorter, wider, with a clockwise arrow: $\\zeta<0$). "
        "The caption says $(\\zeta+f)/h$ is the same for all three.",
    read="The sketch is drawn for the northern hemisphere ($f>0$). When a column moves into deeper water it becomes taller and "
         "thinner and, like the skater, gains spin; over shallower water it loses spin. For $f<0$ the two arrows reverse.",
    change="…the bottom were flat and the surface nearly so? Then $h$ hardly changes and only a change of latitude can alter "
           "$\\zeta$ — the Rossby wave of `C15`.")
note("N97 [B]", r"""
**Nonlinear $x$-momentum:** $\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}-fv=-g\dfrac{\partial\eta}{\partial x}$ (13.88).""")
note("N98 [B]", r"""
**Nonlinear $y$-momentum:** $\dfrac{\partial v}{\partial t}+u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.89).""")
note("N99 [B]", r"""
**Continuity for the total depth:** $\dfrac{\partial h}{\partial t}+\dfrac{\partial}{\partial x}(uh)+\dfrac{\partial}{\partial y}(vh)=0$ (13.90),
with $h$ the total depth over an uneven bottom and $f=f_0+\beta y$. These three are the Start of `D17`; unlike `C07` they are
nonlinear.""")
nb.recap("R25", "Relative and absolute vorticity", r"""
The quantity $\zeta\equiv\dfrac{\partial v}{\partial x}-\dfrac{\partial u}{\partial y}$ is the vertical vorticity measured by someone turning with
the earth (relative); adding the planet's own, $\zeta+f$, gives the absolute vorticity.""", where="Ch. 3 and Ch. 5 §5.6")
nb.current_core = "C13"
P("P325", "materially conserved: Dq/Dt = 0 labels a parcel; it does not mean q is steady at a point", r"""
The operator $D/Dt=\frac\partial{\partial t}+u\frac\partial{\partial x}+v\frac\partial{\partial y}$ (Chapter 3's material derivative)
follows one parcel. The statement $Dq/Dt=0$ says: each parcel keeps its own value of $q$ for ever, like a dye. At a fixed place $q$ can still
change, because parcels with different values pass by. Steady means $\partial q/\partial t=0$ — a different statement.""", code=r"""
x = np.linspace(0, 10, 11); q0 = lambda x: x**2      # each parcel's label, set by where it started
U, t = 2.0, 1.5                                      # a uniform flow carries the parcels 3 units to the right
print(q0(x - U*t)[5], q0(x)[5])                      # at the fixed point x = 5: 4.0 now, 25.0 before — not steady, yet every parcel kept its q
""")
remind("C13")
D("D17", ref="13.93")
note("N100 [B]", r"""
**Cross-differentiated momentum** (steps 1–3 of `D17`):
$\dfrac{\partial}{\partial t}\Big(\dfrac{\partial v}{\partial x}-\dfrac{\partial u}{\partial y}\Big)+\dfrac{\partial}{\partial x}\Big[u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}\Big]-\dfrac{\partial}{\partial y}\Big[u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}\Big]+f_0\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)+\beta v=0$ (13.91).""")
note("N101 [B]", r"""
**The vorticity equation** (step 6): $\dfrac{D\zeta}{Dt}+(\zeta+f_0)\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)+\beta v=0$ (13.92).""")
trap("T10", r"""
The book replaces $f$ by $f_0$ in the coefficient and then restores it. The exact shallow-water result needs no such step — the
check cell above keeps $f=f_0+\beta y$ throughout.""")
nb.recap("R26", "Potential vorticity is conserved following the motion", r"""
$\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=0$, $f=f_0+\beta y$ (13.94) — the boxed result of `D17`. Chapter 5 obtained the
f-plane version from Kelvin's theorem for a column; here it comes from the full nonlinear equations, with $f$ varying.""",
         where=r"Ch. 5 §5.6, derivation D20: $(\omega_z+2\Omega)/h$ constant for a column")
nb.current_core = "C13"
code(r"""
print("D/Dt[(zeta + f)/h] = 0 follows from the nonlinear shallow-water set ->", ch13.pv_conservation_sympy()["ok"])   # the engine agrees: True
""", explain=r"""
**What does this show?** The engine forms $(\zeta+f)/h$ from symbolic fields, takes its material derivative and substitutes
the three equations: nothing is left.""")
nb.recap("R27", "How to read it", r"""
Stretching creates relative vorticity even from $\zeta=0$, because the column carries $f$. There is no tilting term: a shallow
layer moves in vertical columns. An increase of $h$ makes $\zeta+f$ more positive in the northern hemisphere ($f>0$) and more
negative in the southern ($f<0$) — larger in magnitude in both.""", where="Ch. 5 §5.6 (C11)")
nb.current_core = "C13"
nb.worked_example("a column crossing a ridge", r"""
Take $f=10^{-4}$ s⁻¹; a column with no relative vorticity ($\zeta=0$) and depth $h_0=4000$ m moves onto a plateau where
$h_1=3600$ m, without changing latitude.

1. Before: $q=(0+f)/h_0=10^{-4}/4000=2.5\times10^{-8}$ m⁻¹ s⁻¹.
2. After: $(\zeta+f)/h_1=q$ → $\zeta+f=2.5\times10^{-8}\times3600=0.9\times10^{-4}$ s⁻¹.
3. So $\zeta=-0.1\times10^{-4}$ s⁻¹ $=-0.1f$: clockwise (anticyclonic) spin in the northern hemisphere.
4. The same in one step: $\zeta=f(h_1-h_0)/h_0=10^{-4}\times(-400/4000)$.
5. Instead keep the depth and move the column 500 km north, where $f$ is larger by $\beta\Delta y=2\times10^{-11}\times5\times10^{5}=10^{-5}$ s⁻¹:
   $\zeta=-10^{-5}$ s⁻¹ again. A 10 % squash equals a 500 km trip north.""")
code(r"""
q_0 = GFD.potential_vorticity(0.0, f35, 9000.0)                # Eq. (13.94): (zeta + f)/h for a column at rest in a 9000 m layer at 35 N
z_1 = GFD.step_vorticity(f35, 9000.0, 8550.0)                  # its relative vorticity after being squashed by 5 % [1/s]
z_5 = ch05.column_relative_vorticity(8550.0, 9000.0, 0.0, f35) # Chapter 5's column function: the same number
z_s = GFD.step_vorticity(-f35, 9000.0, 8550.0)                 # the southern twin (f < 0)
print(f"q = {q_0:.3e} 1/(m s); after a 5 % squash at 35 N: zeta = {z_1:+.3e} 1/s = {z_1/f35:+.2f} f (clockwise); Ch. 5 function: {z_5:+.3e}")   # anticyclonic in the north
print(f"at 35 S: zeta = {z_s:+.3e} 1/s (counter-clockwise, which is anticyclonic there too)")   # anticyclonic in the south
""", explain=r"""
1. `GFD.potential_vorticity(zeta, f, h)` is $(\zeta+f)/h$; `GFD.step_vorticity(f, h0, h1)` solves "same $q$ before and after"
   for the new $\zeta$.
2. A 5 % squash produces a relative vorticity of 5 % of $f$, anticyclonic in both hemispheres.
3. Chapter 5's `column_relative_vorticity` — written from Kelvin's theorem — returns the same number.""")
scratch(r"""
# From scratch: keep (zeta + f)/h and solve for zeta
q_before = (0.0 + f35)/9000.0                                  # potential vorticity before the step
zeta_after = q_before*8550.0 - f35                             # the relative vorticity that keeps it over the shallower depth
assert np.isclose(zeta_after, GFD.step_vorticity(f35, 9000.0, 8550.0))               # same as the library
assert np.isclose(GFD.potential_vorticity(zeta_after, f35, 8550.0), q_before)        # and q is indeed unchanged
print(f"✓ zeta after the step = {zeta_after:+.3e} 1/s; q before = q after = {q_before:.3e} 1/(m s)")   # visible confirmation
""", r"""
Two lines of arithmetic are the whole content of the conservation law for a column that changes depth at a fixed latitude.""")
D("D18")
note("N102 [B]", r"""
**Eastward flow over a step.** Just downstream of a step from depth $h_0$ to $h_1<h_0$, $\dfrac{f}{h_0}=\dfrac{\zeta+f}{h_1}$ gives
$\zeta=\dfrac{f(h_1-h_0)}{h_0}<0$; the column turns equatorward, where $f$ is smaller, overshoots, and a standing meander of
wavelength $\lambda=2\pi\sqrt{U/\beta}$ follows.""")
note("N103 [B]", r"""
**Westward flow over the same step:** no oscillation — the displacement is exponential on both sides and begins *before* the
step. **Ours** (the book argues in words): the streamline equation $Y''+(\beta/U)Y=\text{const}$ and its solutions (`D18`).""")
fig(r"""
x_s = np.linspace(-4.0e6, 1.2e7, 801)                          # distance along the flow; the step is at x = 0 [m]
east = ch13.flow_over_step(x_s, 17.0, beta35, f35, 9000.0, 8550.0)    # eastward flow at 17 m/s over a 5 % step (ours, D18)
west = ch13.flow_over_step(x_s, -17.0, beta35, f35, 9000.0, 8550.0)   # westward flow over the same step
fig, axes = plt.subplots(2, 2, figsize=(9.6, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
for col, (r_, ttl, arrow) in enumerate(((east, "eastward flow: a wake of standing waves", "→ flow"), (west, "westward flow: felt upstream, no wake", "← flow"))):   # two columns
    for y0 in np.linspace(-600.0, 600.0, 7):                   # seven streamlines, starting at different latitudes [km]
        axes[0, col].plot(x_s/1e6, y0 + r_["Y"]/1e3, color=COLORS["accent"], lw=1)   # each displaced by Y(x)
    lo_hi = (0, x_s[-1]/1e6) if col == 0 else (x_s[0]/1e6, 0)  # the shallower side lies DOWNSTREAM of the step: x > 0 for eastward, x < 0 for westward flow
    axes[0, col].axvspan(*lo_hi, color=COLORS["grid"], alpha=0.6)   # shade it
    axes[0, col].set(title=ttl, ylim=(-1200, 800))
    axes[0, col].text(-3.8, 680, arrow, fontsize=9)            # the direction of the flow
    axes[1, col].plot(x_s/1e6, r_["zeta"]/f35, color=COLORS["amber"])   # relative vorticity in units of f (amber = vorticity)
    axes[1, col].set(xlabel="x [1000 km] (shaded above: shallower by 5 %)")
axes[0, 0].set_ylabel("y north [km]"); axes[1, 0].set_ylabel("ζ / f")
plt.show()
Y_p = f35*(8550.0 - 9000.0)/(beta35*9000.0)                    # the permanent shift f0 (h1 - h0)/(beta h0) [m]
print(f"eastward: Y swings between {east['Y'].max()/1e3:.0f} and {east['Y'].min()/1e3:.0f} km, wavelength {east['wavelength']/1e3:.0f} km")   # the standing meander
print(f"westward: Y = {west['Y'][np.argmin(abs(x_s))]/1e3:.0f} km at the step, tends to {west['Y'][0]/1e3:.0f} km far downstream; e-folding length {west['decay_length']/1e3:.0f} km; Y_p = {Y_p/1e3:.0f} km")   # the permanent shift
""",
    see="Top: streamlines of a uniform current meeting a step to shallower water (shaded: the shallower, downstream side), "
        "eastward flow on the left, westward flow (coming from the right) on the right. Bottom: the relative vorticity along a "
        "streamline.",
    read="Eastward: nothing happens upstream; behind the step the stream swings equatorward and meanders about a line shifted "
         "toward the equator — a wake of stationary Rossby waves. Westward: the deflection begins *before* the step, there is "
         "no wake, and the stream settles at the shifted latitude.",
    change="…$U$ were four times larger? The meander would be twice as long ($\\lambda\\propto\\sqrt{U}$), and the upstream "
           "influence of the westward case would reach twice as far.")
slip(14, r"(as a description and a sketch) that for westward flow over a step the stream returns to its original latitude far downstream",
     r"""a **permanent** shift. Conservation of potential vorticity, $\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=0$ (13.94), with the
flow uniform again far downstream ($\zeta\to0$) requires $f_0+\beta Y_p=f_0h_1/h_0$ there, i.e.
$Y_p=f_0(h_1-h_0)/(\beta h_0)$ — about 220 km toward the equator for our 5 % step (printed above).""")
nb.md(r"""
**The caveat to slip #14, stated honestly.** The printed description *is* right for a **finite ridge**: if the depth returns to
$h_0$ behind the obstacle, the column returns to its latitude. For a *step* — a permanent change of depth, which is what the
text describes — it cannot. No equation of the book is affected; only the description of this one case.""")
choice(r"""
the next figure uses **our numerical model, cached**: a nonlinear C-grid run on a 48 × 48 grid with 24 marked particles
(`ch13.load_reference_run("pv_particles")`).""")
fig(r"""
pvr = ch13.load_reference_run("pv_particles")                  # our cached nonlinear run: particle tracks and their potential vorticity, or None
if pvr is None:                                                # no cache: say so and stop
    print("the cached run is absent — the closed-form figure above stands on its own")
else:
    q_p = np.asarray(pvr["q"], dtype=float)                    # potential vorticity of each particle at 13 times [1/(m s)]
    t_d = np.asarray(pvr["t"], dtype=float)/86400              # the 13 saved times [days]
    ix = np.abs(pvr["x"][None, None, :] - np.asarray(pvr["xp"], float)[:, :, None]).argmin(axis=2)   # nearest grid column of each particle
    iy = np.abs(pvr["y"][None, None, :] - np.asarray(pvr["yp"], float)[:, :, None]).argmin(axis=2)   # nearest grid row
    h_p = float(pvr["H"]) + np.asarray(pvr["eta"], float)[np.arange(len(t_d))[:, None], iy, ix]      # local depth h = H + eta at each particle [m]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, 3.6), sharey=True)
    a.plot(t_d, q_p/q_p[0] - 1, color=COLORS["amber"], lw=0.8) # relative change of (zeta + f)/h along each track
    a.set(xlabel="time [days]", ylabel="relative change since t = 0", title="Potential vorticity of 24 particles")
    b.plot(t_d, h_p/h_p[0] - 1, color=COLORS["blue"], lw=0.8)  # relative change of the local depth along each track
    b.set(xlabel="time [days]", title="Their depth h (same scale)")
    plt.show()
    rms = lambda A_: np.sqrt(np.mean((A_/A_[0] - 1)**2))       # root-mean-square relative change over all particles and times
    print(f"rms relative change: potential vorticity {100*rms(q_p):.2f} %, depth {100*rms(h_p):.2f} %; largest: {100*abs(q_p/q_p[0] - 1).max():.2f} % and {100*abs(h_p/h_p[0] - 1).max():.2f} %")   # q is conserved, h is not
""",
    see="Left: the potential vorticity of each of 24 particles in a nonlinear model run, relative to its starting value. "
        "Right: the depth of the layer at the same particles.",
    read="The depth at a particle changes as waves and eddies pass; its potential vorticity changes about half as much (the "
         "printed line gives both; what is left for the potential vorticity is the interpolation error of the diagnostic on a "
         "coarse 48 × 48 grid). The combination $(\\zeta+f)/h$ is what each column keeps.",
    change="…the run had friction? Then potential vorticity would slowly change along each track — it is conserved only for "
           "inviscid flow.")
whatif(r"""
…the stratification is continuous and we ask for waves of any slope, fast or slow? First the fast ones (`C14`); the slow
ones, which are this conservation law in motion, are Rossby waves (`C15`).""")

# =====================================================================================================================
# A.14  §13.14 — C14 inertia–gravity waves
# =====================================================================================================================
nb.section("13.14", "Internal Waves", intro=r"""
**What is this section about?** Chapter 7's internal gravity waves with the earth's rotation added. The frequency is still set
by the direction of the wavevector alone, but now it is caught between $f$ and $N$. ⚠️ In this section the book's $\theta$ is
the angle of the wavevector with the horizontal, not latitude, and its $H$ is the scale over which $N$ changes, not a depth.""")
core("C14", r"Inertia–gravity waves: $\omega^2-f^2=\dfrac{k^2}{m^2}(N^2-\omega^2)$ (13.112); equivalently, in a form the book prints without a number, $\omega^2=f^2\sin^2\theta+N^2\cos^2\theta$",
     "What does rotation change about the internal waves of Chapter 7?")
problem(r"""
Lower a current meter into the thermocline and the record is full of oscillations with periods from about ten minutes to most
of a day. The short end is the buoyancy period $2\pi/N$; the long end is the inertial period $2\pi/f$. In between, every period
belongs to a wave whose crests have a definite slope: steep crests, fast wave; nearly flat crests, slow wave, almost an
inertial circle. These waves carry the energy of wind and tide into the deep ocean, where their breaking does much of the
mixing that ocean models have to parameterise.""")
idea(words=r"""
Two springs again: buoyancy (stiffness $N^2$) resists vertical motion, rotation (stiffness $f^2$) resists horizontal motion.
Fluid moves along the crests. With $\theta_K$ the angle of the wavevector above the horizontal (the book's $\theta$ here):

| Wavevector | Crests | Motion | Frequency |
|---|---|---|---|
| horizontal, $\theta_K=0°$ | vertical | vertical: only buoyancy acts | $\omega=N$ |
| $\theta_K=45°$ | tilted 45° | both | $\omega^2=(f^2+N^2)/2$ |
| vertical, $\theta_K=90°$ | horizontal | horizontal: only rotation acts | $\omega=f$ |""")
nb.recap("R28", "Internal waves without rotation", r"""
Anisotropic; $\omega\le N$; frequency set by the angle of the wavevector, $\omega=N\cos\theta$; phase and group velocity
perpendicular, with opposite vertical components.""", where="Ch. 7 §7.8 (C15–C16)")
nb.current_core = "C14"
nb.recap("R29", "The linear rotating Boussinesq set", r"""
Continuity $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$;
$\dfrac{\partial u}{\partial t}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$; $\dfrac{\partial v}{\partial t}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$;
$\dfrac{\partial w}{\partial t}=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial z}-\dfrac{\rho g}{\rho_0}$;
$\dfrac{\partial\rho}{\partial t}-\dfrac{\rho_0N^2}{g}w=0$ (13.95). New against Chapter 7: the two Coriolis
terms. Not hydrostatic.""", where=r"Ch. 7, the same set with $f=0$, ending in $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ (7.134)")
nb.current_core = "C14"
code(r"""
print("the five-equation set (13.95) as coded is consistent (sympy engine) ->", ch13.rotating_internal_wave_set_sympy()["ok"])   # True
""", explain=r"""
**What does this show?** The starting set of `D19`, checked by its engine before we eliminate four of its five unknowns.""")
remind("C14")
D("D19", ref="13.96")
note("N105 [B]", r"""
**The single equation for $w$:** $\dfrac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w+f^2\dfrac{\partial^2w}{\partial z^2}=0$ (13.96) —
written out in `D19` because the book refers to its §7.8 and never shows the steps with rotation.""")
slip(8, r'that $N$ is taken "depth independent" in §13.14',
     r"depth dependent, $N(z)$, as its own definition $m^2(z)\equiv\frac{(k^2+l^2)[N^2(z)-\omega^2]}{\omega^2-f^2}$ (13.99) requires.")
code(r"""
k_w, m_w = 2*np.pi/1.0e4, 2*np.pi/200.0                        # a wave 10 km long and 200 m tall [rad/m]
w_ig = GFD.inertia_gravity_omega(k_w, m_w, N_o, f35)           # its frequency from Eq. (13.112) [rad/s]
w_fn = lambda xx, yy, zz, tt: np.cos(k_w*xx + m_w*zz - w_ig*tt)   # a plane wave of vertical velocity
res_w = ch13.w_equation_rotating_residual(w_fn, 0.0, 0.0, 0.0, 0.0, N_o, f35)   # what is left of Eq. (13.96) by finite differences
print("the w-equation (13.96) follows from the set (13.95) ->", ch13.w_equation_rotating_sympy()["ok"])   # the engine agrees: True
print(f"a plane wave with omega from (13.112) put into (13.96): residual = {abs(res_w)/(N_o**2*k_w**2):.1e} of the N^2 k^2 term")   # the dispersion relation solves it
""", explain=r"""
**What does this show?** Two checks of `D19`: the symbolic elimination gives the printed equation, and a plane wave whose
frequency comes from the dispersion relation satisfies it to finite-difference accuracy.""")
D("D20", ref="13.112")
note("N106 [C]", r"""
**The trial solution** (step 1 of `D20`): $[u,v,w]=[\hat u(z),\hat v(z),\hat w(z)]\,e^{i(kx+ly-\omega t)}$ (13.97).""")
note("N107 [B]", r"""
**The vertical-structure equation** (step 2): $\dfrac{d^2\hat w}{dz^2}+\dfrac{(N^2-\omega^2)(k^2+l^2)}{\omega^2-f^2}\hat w=0$ (13.98).""")
note("N108 [B]", r"""
**The local vertical wavenumber** (step 3): $m^2(z)\equiv\dfrac{(k^2+l^2)[N^2(z)-\omega^2]}{\omega^2-f^2}$ (13.99)
(`GFD.inertia_gravity_m2`).""")
note("N109 [B]", r"""
**The oscillator form:** $\dfrac{d^2\hat w}{dz^2}+m^2\hat w=0$ (13.100).""")
note("N119 [C]", r"""
**For $l=0$** (step 4): $m^2=\dfrac{k^2(N^2-\omega^2)}{\omega^2-f^2}$ (13.109).""")
note("N110 [B]", r"""
**The band.** Internal waves exist only for $f<\omega<N$: there $m^2>0$ and the structure oscillates in $z$; outside the band
$m^2<0$ and it decays (it is *evanescent*).""")
code(r"""
for w_, nm in ((0.5*f35, "0.5 f"), (5*f35, "5 f"), (1.2*N_o, "1.2 N")):   # three frequencies: below, inside and above the band
    bd = GFD.inertia_gravity_band(w_, N_o, f35)                # where this frequency lies
    m2 = GFD.inertia_gravity_m2(k_w, 0.0, w_, N_o, f35)        # Eq. (13.99): m^2 for a 10 km wave [1/m^2]
    print(f"omega = {nm:>5}: {bd['where']:>8} | m^2 = {m2:+.2e} 1/m^2 | {'propagates' if bd['propagating'] else 'evanescent'}")   # our own label, built from the dictionary
""", explain=r"""
**What does this show?** Below $f$ and above $N$ the squared vertical wavenumber is negative — no wave; in between it is
positive. (`N_o` is our ocean's $N=2.7\times10^{-3}$ s⁻¹; $N/f=32$ at 35° N.)""")
nb.worked_example("frequency from the slope of the crests", r"""
Take $f=10^{-4}$ s⁻¹, $N=10^{-2}$ s⁻¹ ($N/f=100$).

1. Wavevector 60° from the horizontal: $\omega^2=f^2\sin^260°+N^2\cos^260°=0.75\times10^{-8}+0.25\times10^{-4}\approx0.25\times10^{-4}$ →
   $\omega=5\times10^{-3}$ s⁻¹ $=N/2$ (period 21 min): rotation is invisible.
2. Wavevector almost vertical, $m/k=100$ ($\tan\theta=100$): $\cos^2\theta\approx10^{-4}$, $\sin^2\theta\approx1$ →
   $\omega^2=10^{-8}+10^{-4}\times10^{-4}=2\times10^{-8}$ → $\omega=1.41f$ (period 12.3 h): half rotation, half buoyancy.
3. With $m/k=1000$: $\omega=1.005f$ — an inertial oscillation.""")
code(r"""
for th_deg in (5, 45, 85, 89):                                 # angle of the wavevector above the horizontal [degrees]
    w_ = GFD.inertia_gravity_omega(np.cos(np.deg2rad(th_deg)), np.sin(np.deg2rad(th_deg)), N_o, f35)   # only the direction of (k, m) matters
    print(f"angle {th_deg:2d} deg: omega = {w_/N_o:.3f} N = {w_/f35:5.2f} f")   # from N (flat wavevector) down to f (vertical wavevector)
cg = GFD.inertia_gravity_group_velocity(k_w, m_w, N_o, f35)    # group velocity (c_gx, c_gz) of the 10 km x 200 m wave [m/s]
cp = w_ig/(k_w**2 + m_w**2)*np.array([k_w, m_w])               # phase velocity vector: (omega/K^2) K [m/s]
print(f"10 km x 200 m wave: omega = {w_ig:.4e} 1/s = {w_ig/f35:.3f} f (period {2*np.pi/w_ig/3600:.1f} h)")   # near-inertial
print(f"group velocity = ({cg[0]:+.4f}, {cg[1]:+.2e}) m/s; phase speed = {np.hypot(*cp)*1e3:.2f} mm/s; c_p · c_g = {np.dot(cp, cg):.1e}")   # perpendicular
print(f"f = 0 check against Chapter 7: {GFD.inertia_gravity_omega(k_w, m_w, N_o, 0.0):.6e} vs {WAV.internal_wave_omega(k_w, m_w, N_o):.6e} 1/s")   # the non-rotating limit
""", explain=r"""
1. `GFD.inertia_gravity_omega(k, m, N, f)` is the dispersion relation solved for $\omega$; only the direction of $(k, m)$ matters.
2. As the wavevector turns from horizontal to vertical the frequency falls from $N$ to $f$.
3. For a wave 10 km long and 200 m tall the frequency is 1.19 $f$; its energy travels almost horizontally and slightly
   *downward* while its crests move upward; phase and group velocity are perpendicular.
4. With $f=0$ the function reduces to Chapter 7's `WAV.internal_wave_omega`.""")
P("P326", "group velocity as the gradient of ω in wavenumber space, read off a contour plot", r"""
In one dimension the energy of a wave packet travels at $d\omega/dk$. In two or three, at the vector
$(\partial\omega/\partial k,\ \partial\omega/\partial l,\ \partial\omega/\partial m)$: the gradient of $\omega$ in wavenumber space. On a
contour plot of $\omega$ over the $(k, m)$ plane it points straight uphill, at right angles to the contours — and its length is
how crowded the contours are.""", code=r"""
w = lambda k, m: np.sqrt((k**2 + 0.01*m**2)/(k**2 + m**2))    # N = 1, f = 0.1
h = 1e-6; k0, m0 = 1.0, 2.0                                    # a wavevector pointing forward and up
print("group velocity (c_gx, c_gz) =", round((w(k0+h, m0)-w(k0-h, m0))/(2*h), 3), round((w(k0, m0+h)-w(k0, m0-h))/(2*h), 3))   # (0.354, -0.177): energy goes forward and DOWN while crests go forward and up
""")
D("D21")
note("N124 [B]", r"""
**The group velocity** (the book leaves it to an exercise; `D21`):
$[c_{gx},c_{gz}]=\dfrac{(N^2-f^2)\,km}{(m^2+k^2)^{3/2}(m^2f^2+k^2N^2)^{1/2}}\,[m,\,-k]$ — phase velocity along the wavevector,
group velocity at right angles to it with the opposite vertical sign, fluid motion along the group velocity.""")
scratch(r"""
# From scratch: the unnumbered sin^2/cos^2 form, and the group velocity by finite differences
th_K = np.arctan2(m_w, k_w)                                    # angle of the wavevector above the horizontal [rad]
w_mine = np.sqrt(f35**2*np.sin(th_K)**2 + N_o**2*np.cos(th_K)**2)   # omega^2 = f^2 sin^2 + N^2 cos^2
hk, hm = 1e-6*k_w, 1e-6*m_w                                    # small steps in k and m
cg_fd = ((GFD.inertia_gravity_omega(k_w + hk, m_w, N_o, f35) - GFD.inertia_gravity_omega(k_w - hk, m_w, N_o, f35))/(2*hk),
         (GFD.inertia_gravity_omega(k_w, m_w + hm, N_o, f35) - GFD.inertia_gravity_omega(k_w, m_w - hm, N_o, f35))/(2*hm))   # (d omega/dk, d omega/dm)
assert np.isclose(w_mine, w_ig) and np.allclose(cg_fd, cg, rtol=1e-6)   # same frequency, same group velocity
print(f"✓ omega = {w_mine:.4e} 1/s from the angle form; group velocity by differences ({cg_fd[0]:+.4f}, {cg_fd[1]:+.2e}) m/s matches D21")   # visible confirmation
""", r"""
The frequency from the angle of the wavevector equals the library's, and differencing the dispersion relation in $k$ and $m$
reproduces the closed-form group velocity of `D21`.""")
note("N123 [B]", r"""
**Three regimes** (our version of the book's figure): high frequency, $m^2\simeq\dfrac{k^2(N^2-\omega^2)}{\omega^2}$ (non-rotating,
$\omega=N\cos\theta$); low frequency, $\omega^2\simeq f^2+\dfrac{k^2N^2}{m^2}$ (hydrostatic); mid frequency,
$m^2\simeq\dfrac{k^2N^2}{\omega^2}$ (both approximations at once).""")
fig(r"""
k_ax = np.geomspace(1e-6, 1e-1, 300)                           # horizontal wavenumbers [rad/m]
fig, ax = plt.subplots(figsize=(6.8, 3.9))
for m_, col in zip((2*np.pi/2000.0, 2*np.pi/200.0, 2*np.pi/20.0), (COLORS["accent"], COLORS["teal"], COLORS["rose"])):   # three vertical wavelengths
    ax.loglog(k_ax, GFD.inertia_gravity_omega(k_ax, m_, N_o, f35), color=col, label=f"vertical wavelength {2*np.pi/m_:.0f} m")   # Eq. (13.112)
ax.axhline(f35, color=COLORS["muted"], ls="--"); ax.axhline(N_o, color=COLORS["muted"], ls="--")   # the two limits
ax.text(1.2e-6, f35*1.15, "ω = f", fontsize=8); ax.text(1.2e-6, N_o*1.15, "ω = N", fontsize=8)     # their labels
ax.text(2e-6, 1.6e-4, "near-inertial:\nhydrostatic form", fontsize=7); ax.text(2e-2, 1.2e-3, "near N:\nnon-rotating form", fontsize=7, ha="center")   # the regimes
ax.set(xlabel="horizontal wavenumber k [rad/m]", ylabel="ω [rad/s]", title="Every curve climbs from f to N", ylim=(5e-5, 6e-3))
ax.legend(fontsize=7, loc="center right")
plt.show()
rows = [dict(omega=nm, **{k_: (v_ if isinstance(v_, str) else round(float(v_), 4)) for k_, v_ in GFD.inertia_gravity_regime(w_, N_o, f35).items()})
        for nm, w_ in (("2 f", 2*f35), ("√(f N)", np.sqrt(f35*N_o)), ("N/2", N_o/2))]   # the three approximations at three frequencies
display(pd.DataFrame(rows))                                    # relative error in m^2 of each approximation
""",
    see="Frequency against horizontal wavenumber for three vertical wavelengths, on logarithmic axes; the table gives the "
        "relative error in $m^2$ of each approximate relation at three frequencies.",
    read="Long, flat waves sit just above $f$ (rotation matters, buoyancy enters hydrostatically); short, steep ones approach "
         "$N$ (rotation is invisible). In the middle of the band, far from both limits, both approximations hold at once.",
    change="…the latitude were lower? Then $f$ falls, the band widens downward and the near-inertial regime moves to longer periods.")
P("P327", "slowly varying medium (WKB): amplitude and phase ansatz, valid when the medium changes little in one wavelength", r"""
If the "stiffness" of an oscillator equation $\hat w''+m^2(z)\hat w=0$ changes slowly, the solution still looks locally like a
wave, with a local wavenumber $m(z)$ and a slowly changing amplitude. Write $\hat w=A(z)e^{i\phi(z)}$, and demand that $A$
changes little over one wavelength. It fails where $m\to0$ (a turning point: the wave reflects).""", code=r"""
z = np.linspace(0, 50, 5001); m = 1 + 0.02*z                  # the wavenumber doubles over 50 units
phase = np.cumsum(m)*(z[1]-z[0])                              # the phase is the running integral of m
w_wkb = np.cos(phase)/np.sqrt(m)                              # the amplitude falls like m^(-1/2)
print(f"amplitude at the start {w_wkb[0]:.2f}, near the end {np.abs(w_wkb[-200:]).max():.2f}")   # 1.00 and 0.71
""")
note("N111 [B]", r"""
**The WKB condition:** $H_Nm\gg1$ (the book writes $Hm\gg1$, with $H$ the scale over which $N$ varies): many vertical
wavelengths fit into the distance over which the stratification changes.""")
note("N112 [C]", r"""
**Substituting** $\hat w=A(z)e^{i\phi(z)}$ gives $\dfrac{d^2A}{dz^2}+A\Big[m^2-\Big(\dfrac{d\phi}{dz}\Big)^2\Big]=0$ and
$2\dfrac{dA}{dz}\dfrac{d\phi}{dz}+A\dfrac{d^2\phi}{dz^2}=0$ (13.101)–(13.102).""")
note("N113 [B]", r"""
**The eikonal.** Neglecting $A''$ gives $\dfrac{d\phi}{dz}=\pm m$, $\phi=\pm\displaystyle\int^zm\,dz$ (13.103): the phase gradient is
the local wavenumber.""")
note("N114 [B]", r"""
**The WKB solution.** The second equation then integrates to $\hat w=\dfrac{A_0}{\sqrt m}\,e^{\pm i\int^zm\,dz}$ (13.104): the
amplitude grows where the wave is long, because the vertical energy flux must be the same at every level. (Stated; the book
prints these steps.)""")
fig(r"""
k_wk, w_wk = 1.0e-3, 4.0e-4                                    # horizontal wavenumber [rad/m] and frequency [rad/s] of the test wave
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.8))
Hm_list, err_list = [], []                                     # collected for the right panel
for D_ in (500.0, 1000.0, 2000.0, 4000.0):                     # four thicknesses of the layer over which N doubles [m]
    N_fn = lambda zz, D_=D_: N_o*(1 + 0.5*zz/D_)               # N rises from N/2 at z = -D to N at z = 0
    z_w = np.linspace(-D_, 0.0, 801)                           # the layer
    we = ch13.wkb_error(z_w, N_fn, k_wk, w_wk, f35)            # largest relative error of the WKB solution, and the value of H_N m
    Hm_list.append(we["Hm"]); err_list.append(we["max_rel_error"])
    if D_ == 1000.0:                                           # show the two solutions for one case
        m_z = np.sqrt(GFD.inertia_gravity_m2(k_wk, 0.0, w_wk, N_fn(z_w), f35))   # local vertical wavenumber, Eq. (13.99)
        w_num = ch13.vertical_structure_solve(z_w, N_fn, k_wk, w_wk, f35)        # numerical solution of Eq. (13.100): w = 0, dw/dz = 1 at the bottom
        w_wkb = GFD.wkb_vertical_structure(z_w, m_z, A0=1/np.sqrt(m_z[0])).imag  # Eq. (13.104), the combination with the same starting values
        a.plot(np.real(w_num), z_w, color=COLORS["accent"], lw=3, label="numerical")
        a.plot(w_wkb, z_w, "--", color=COLORS["ink"], lw=1, label="WKB, Eq. (13.104)")
a.set(xlabel=r"$\hat w$ (arbitrary units)", ylabel="z [m]", title="N doubles over 1000 m: WKB follows the wave")
a.legend(fontsize=8)
b.loglog(Hm_list, err_list, "o-", color=COLORS["accent"], label="measured")   # error against H_N m
b.loglog(Hm_list, err_list[0]*Hm_list[0]/np.array(Hm_list), "--", color=COLORS["muted"], label="slope −1")   # a guide proportional to 1/(H_N m)
b.set(xlabel=r"$H_Nm$", ylabel="largest relative error", title="The error falls like $1/(H_Nm)$")
b.legend(fontsize=8)
plt.show()
print("H_N m:", np.round(Hm_list, 2), "| error:", np.round(err_list, 4), "| error x H_N m:", np.round(np.array(err_list)*np.array(Hm_list), 3))   # measured
""",
    see="Left: the vertical structure of a wave in a layer whose $N$ doubles from bottom to top — numerical solution (purple) "
        "and WKB formula (dashed). Right: the largest relative error of the WKB formula against $H_Nm$.",
    read="The wave is shorter and smaller where $N$ is larger, exactly as $A_0/\\sqrt m$ says. The error is inversely "
         "proportional to $H_Nm$ (the points follow the slope −1 guide; the printed product is roughly constant) — not to its "
         "square.",
    change="…the frequency approached the local $N$ somewhere? Then $m\\to0$ there, the WKB formula blows up and the wave "
           "reflects (a turning point).")
note("N115 [C]", r"""
**The velocity field** (stated): continuity for $l=0$, $ik\hat u+\dfrac{d\hat w}{dz}=0$ (13.105).""")
note("N116 [C]", r"""
$\hat u=\mp\dfrac{A_0\sqrt m}{k}e^{\pm i\int^zm\,dz}$ (13.106). ⚠️ Trap T15: $\sqrt m$ is treated as constant when $\hat w$ is
differentiated.""")
note("N117 [C]", r"""
$\hat v=\pm\dfrac{if}{\omega}\dfrac{A_0\sqrt m}{k}e^{\pm i\int^zm\,dz}$ (13.107), from $\hat u/\hat v=i\omega/f$.""")
note("N118 [B]", r"""
**The real fields:** $u=\mp\dfrac{A_0\sqrt m}{k}\cos\Big(kx\pm\displaystyle\int^zm\,dz-\omega t\Big)$,
$v=\mp\dfrac{A_0f\sqrt m}{\omega k}\sin\Big(kx\pm\displaystyle\int^zm\,dz-\omega t\Big)$,
$w=\dfrac{A_0}{\sqrt m}\cos\Big(kx\pm\displaystyle\int^zm\,dz-\omega t\Big)$ (13.108), upper signs for upward phase propagation
(`GFD.inertia_gravity_fields`).""")
note("N120 [B]", r"""
**The horizontal hodograph at a point:** $u=\mp\cos\omega t$, $v=\pm\dfrac f\omega\sin\omega t$ (13.110): a clockwise ellipse in the
northern hemisphere ($f>0$; counter-clockwise for $f<0$), long axis along the direction of propagation, axis ratio $f/\omega$
(trap T7: the same ellipse as the $\omega/f$ of `C10`, quoted the other way up).""")
note("N121 [B]", r"""
**The motion lies along the phase lines:** $\dfrac uw=\mp\dfrac mk=\mp\tan\theta$ (13.111), with $\theta=\tan^{-1}(m/k)$ the angle of
the wavevector with the horizontal (**not** latitude — trap T13).""")
code(r"""
N_uni = lambda zz: N_o + 0*zz                                  # uniform stratification, as a function of z
m_uni = np.sqrt(GFD.inertia_gravity_m2(k_w, 0.0, w_ig, N_o, f35))   # the vertical wavenumber that belongs to (k_w, w_ig): 2 pi/200 m
x_g, z_g = np.linspace(0.0, 1.0e4, 5), np.linspace(-200.0, 0.0, 5)  # a small grid: one horizontal and one vertical wavelength
u_g, v_g, w_g = GFD.inertia_gravity_fields(x_g, z_g, 0.0, k_w, w_ig, N_uni, f35)   # Eq. (13.108), upward phase propagation
u_h, v_h = GFD.inertia_gravity_hodograph(np.linspace(0, 2*np.pi/w_ig, 5), w_ig, f35)   # Eq. (13.110) over one period
print(f"k u + m w (motion along the crests) = {abs(k_w*u_g + m_uni*w_g).max()/abs(k_w*u_g).max():.1e} of k u")   # Eq. (13.111): zero
print(f"amplitude ratio |v|/|u| = {abs(v_g).max()/abs(u_g).max():.3f} = f/omega = {f35/w_ig:.3f}; hodograph points (u, v): {np.round(u_h, 2)}, {np.round(v_h, 2)}")   # the ellipse
""", explain=r"""
**What does this show?** On a small grid in a uniformly stratified layer: the velocity is perpendicular to the wavevector
($ku+mw=0$, i.e. the water moves along the crests), and the cross-wave current is $f/\omega$ times the along-wave one.""")
note("N122 [B]", r"""
**The hodograph in both hemispheres** (our version of the book's sketches, `ch13.fig_inertia_gravity_orbit`): the horizontal
velocity at a fixed point over one wave period, for $f>0$ and for $f<0$. (That the water also moves along the sloping crests,
and that phase and group velocity are at right angles, is the content of note N121 and of derivation `D21` above.)""")
fig(r"""
figO = ch13.fig_inertia_gravity_orbit()                        # two horizontal velocity hodographs: f > 0 and f < 0, for omega = 2.5 |f|
plt.show()
""",
    see="Two ellipses: the tip of the horizontal velocity vector $(u, v)$, in units of the amplitude of $u$, over one period — "
        "left for $f>0$, right for $f<0$. Dots mark the phases $\\omega t=0$, $\\pi/2$ and $\\pi$.",
    read="Follow the dots from $\\omega t=0$ to $\\pi/2$ to $\\pi$: on the left the vector passes through positive $v$ on its "
         "way from $-u$ to $+u$, which is clockwise; on the right it passes through negative $v$, counter-clockwise. The long "
         "axis is along the direction of propagation and the axis ratio is $f/\\omega=0.40$ here ($\\omega=2.5\\lvert f\\rvert$).",
    change="…$\\omega$ approached $\\lvert f\\rvert$? The axis ratio tends to 1 and the ellipse becomes the inertial circle of `C10`.")
nb.plotly(r"""
z_hx = np.linspace(-400.0, 0.0, 161)                           # two vertical wavelengths of the 200 m wave [m]
figH = go.Figure()
for j, (sg, nm) in enumerate(((+1, "upward phase (energy downward)"), (-1, "downward phase (energy upward)"))):   # the two signs of Eq. (13.108)
    uh_, vh_, _ = GFD.inertia_gravity_fields(np.array([0.0]), z_hx, 0.0, k_w, w_ig, N_uni, f35, sign=sg)   # velocity at x = 0, t = 0, against depth
    scale = abs(uh_).max()                                     # normalise by the largest u
    figH.add_trace(go.Scatter3d(x=uh_[:, 0]/scale, y=vh_[:, 0]/scale, z=z_hx, mode="lines", line=dict(color=COLORS["accent"], width=6), name=nm, visible=(j == 0)))   # the helix
figH.update_layout(height=520, title="The current vector turns with depth: a helix (35° N)",
                   updatemenus=[dict(x=0.0, y=1.08, xanchor="left", buttons=[
                       dict(label="upward phase propagation", method="update", args=[{"visible": [True, False]}]),
                       dict(label="downward phase propagation", method="update", args=[{"visible": [False, True]}])])],
                   scene=dict(xaxis_title="u (normalised)", yaxis_title="v (normalised)", zaxis_title="z [m]"), margin=dict(l=0, r=0, t=60, b=0))   # axes and dropdown
figH.show()
""", explain=r"""
The horizontal velocity at one place and time, plotted against depth: its tip traces a helix. The dropdown switches between
the two signs of the vertical wavenumber.""")
nb.figure_notes(
    see="The tip of the horizontal current vector at successive depths: a flattened helix (an ellipse of axis ratio $f/\\omega$ "
        "seen from above).",
    read="Looking down from above in the northern hemisphere, the vector turns clockwise with increasing depth when the phase "
         "moves upward (energy going down) and counter-clockwise when it moves downward — the test oceanographers use to tell "
         "which way near-inertial energy is travelling.",
    change="…$f<0$? Both senses reverse.")
note("N125 [B]", r"""
**Lee waves.** With $f$ negligible, $\omega^2=\dfrac{N^2k^2}{m^2+k^2}$ (13.113). A wave that stands still over the ground has its
frequency Doppler-shifted to zero, $\omega_0=\omega+\mathbf K\cdot\mathbf U=0$ (Chapter 7's $\omega_0=\omega+\mathbf U\cdot\mathbf K$ (7.9)), so
$\omega=kU$ and $U=\dfrac{N}{\sqrt{k^2+m^2}}$: only waves with $k<N/U$ exist. (Trap T15: $k$ here is a magnitude; the wavevector
points upstream.)""")
code(r"""
U_m, N_a = inp["U_mean"], inp["atm_N"]                         # a 17 m/s wind over a mountain; atmospheric N = 1.1e-2 1/s (our inputs)
m_lee = GFD.lee_wave_m(U_m, N_a, 2*np.pi/2.0e4)                # vertical wavenumber of a stationary wave 20 km long [1/m]
print(f"shortest stationary wavelength 2 pi U/N = {2*np.pi*U_m/N_a/1e3:.2f} km; a 20 km wave has m = {m_lee:.2e} 1/m (vertical wavelength {2*np.pi/m_lee/1e3:.1f} km)")   # the cut-off and one wave
try:                                                           # a 5 km wave is shorter than the cut-off
    GFD.lee_wave_m(U_m, N_a, 2*np.pi/5.0e3)
except ValueError as err:                                      # the function refuses
    print("a 5 km wave ->", err)
""", explain=r"""
**What does this show?** In a 17 m/s wind with $N=1.1\times10^{-2}$ s⁻¹ no stationary wave shorter than about 10 km exists; a
20 km wave has a vertical wavelength of about 11 km; asking for a 5 km wave raises an error ("evanescent").""")
note("N126 [B]", r"""
**The lee-wave pattern** (our version of the book's sketch, `ch13.lee_wave_field`): streamlines over sinusoidal terrain 300 m
high.""")
fig(r"""
x_l = np.linspace(0.0, 6.0e4, 201); z_l = np.linspace(0.0, 1.2e4, 101)   # three wavelengths downstream, 12 km up [m]
lw = ch13.lee_wave_field(x_l, z_l, U_m, N_a, 2*np.pi/2.0e4, 300.0)        # stream function of wind + stationary wave (ours)
fig, ax = plt.subplots(figsize=(7.6, 3.8))
ax.contour(x_l/1e3, z_l/1e3, lw["psi"], 16, colors=COLORS["accent"], linewidths=0.9)   # streamlines
ax.plot([0, 2*np.pi/lw["m"]/1e3*0 + 5.0], [0, 0], alpha=0)     # (keeps the axes starting at the ground)
ax.set(xlabel="x downstream [km]  (wind →)", ylabel="height z [km]", title=f"Lee waves: the crests tilt {lw['tilt']} with height")
plt.show()
""",
    see="Streamlines of a 17 m/s wind over wavy terrain: each streamline is a copy of the terrain, shifted in phase more and "
        "more with height.",
    read="Join the crests of successive streamlines: the line leans upstream (to the left) with height. That tilt is what a "
         "wave with downward phase velocity and **upward** energy flux looks like — the mountain is the source.",
    change="…the wind were stronger, so that $N/U<k$? The disturbance would decay with height instead of tilting: no wave, no "
           "upward energy flux.")
whatif(r"""
…the frequency is far below $f$? Then none of these waves is possible — but the slow root of the cubic of `C09` is. It needs
$\beta$ (`C15`).""")


# =====================================================================================================================
# A.15  §13.15 — C15 Rossby waves;  A.16  §13.16 — barotropic instability (C15 continued)
# =====================================================================================================================
nb.section("13.15", "Rossby Wave", intro=r"""
**What is this section about?** The slow wave that owes its existence to the change of the Coriolis parameter with latitude.
Its crests always drift west; its energy can go either way.""")
core("C15", r"Rossby waves: $\omega=-\dfrac{\beta k}{k^2+l^2+f_0^2/c^2}$ (13.118)",
     "The crests go west — so how can a storm track's energy go east, and how long does the ocean take to hear about a change in the wind?")
problem(r"""
The jet stream meanders in four or five great waves round the hemisphere; when one of them stalls, a region gets weeks of the
same weather. In the ocean, a change of wind in the east Pacific is felt in the west a year later at low latitudes, a decade
later at mid-latitudes. Both are Rossby waves. They are not held up by gravity like the waves of `C10` but by the conservation
of potential vorticity on a planet where $f$ grows toward the pole.""")
idea(r"""
   north ↑        ↻ (pushed north: f larger, so ζ must fall → clockwise spin)
   ─ ─ ─ ─ ─ ● ─ ─ ─ ─ ─ ─ ● ─ ─ ─ ─ ─   a line of columns along a latitude circle
                            ↺ (pushed south: f smaller, so ζ must rise → counter-clockwise spin)

   between them the two spins push fluid NORTH on the western side of the northern bulge
   and SOUTH on its eastern side  →  the whole pattern shifts WEST
""", r"""
This is the conservation law of `C13` in motion: a column displaced north or south must change its relative spin, and the
induced flow moves the displacement pattern westward (northern hemisphere shown; the direction is westward in both).""")
note("N127 [B]", r"""
**Planetary waves.** Rossby waves exist only because of $\beta$, at $\omega\ll f$. The motion is *quasi-geostrophic*:
geostrophic to lowest order, with its slow evolution set by the small departures from geostrophy. (The book's observed height
map is replaced here by a synthetic wavenumber-5 field with its geostrophic wind.)""")
fig(r"""
xr = np.linspace(0.0, 2.8e7, 141); yr = np.linspace(-2.5e6, 2.5e6, 41)   # once round the 35 N latitude circle (about 28 000 km), 5000 km in y [m]
Xr, Yr = np.meshgrid(xr, yr)                                   # grids indexed [j, i] = (y, x)
eta_r = -300.0*np.tanh(Yr/1.2e6) + 80.0*np.cos(2*np.pi*5*Xr/2.8e7)*np.exp(-(Yr/1.5e6)**2)   # a height field: low toward the pole plus five waves [m]
u_r, v_r = GFD.geostrophic_from_height(eta_r, xr[1] - xr[0], yr[1] - yr[0], f35)   # its geostrophic wind, Eq. (13.116)
fig, ax = plt.subplots(figsize=(9.2, 3.0))
ax.contour(xr/1e6, yr/1e6, eta_r, 12, colors=COLORS["orange"], linewidths=0.9)   # height contours (orange = pressure)
ax.quiver(xr[::5]/1e6, yr[::4]/1e6, u_r[::4, ::5], v_r[::4, ::5], color=COLORS["teal"])   # the wind
ax.set(xlabel="x east [1000 km]", ylabel="y north [1000 km]", title="A synthetic wavenumber-5 pattern: the wind follows the meandering height contours")
plt.show()
""",
    see="A synthetic height field (orange contours; dashed where the height is below average, toward the pole) with five waves "
        "round a latitude circle, and its geostrophic wind (teal arrows, strongest where the contours are closest).",
    read="The westerly jet meanders: poleward on one side of each trough, equatorward on the other. A column riding this flow is "
         "carried north and south — and must change its spin as it goes. That is the restoring mechanism of the wave.",
    change="…there were only two or three waves round the circle? They would be longer, and — as the dispersion relation shows — "
           "drift westward faster relative to the air.")
P("P328", "ordering in a small parameter: lowest order gives the balance, next order gives the evolution", r"""
When a small number $\varepsilon$ (here the Rossby number, or $\omega/f$) multiplies some terms, sort every term by how many
factors of $\varepsilon$ it carries. The largest terms must balance among themselves — that gives a *diagnostic* relation with
no time derivative (geostrophy). To learn how things change you must go to the next size of terms. So the small terms cannot
all be thrown away: the right ones decide the future.""", code=r"""
eps = 0.1                              # a small parameter, e.g. the Rossby number
print(1.0, eps, eps**2)                # sizes of: Coriolis and pressure | acceleration and the divergent wind | products of wave amplitudes
""")
remind("C15")
D("D22", ref="13.117")
note("N128 [B]", r"""
**Potential-vorticity conservation expanded** (step 1 of `D22`):
$(H+\eta)\Big(\dfrac{\partial\zeta}{\partial t}+u\dfrac{\partial\zeta}{\partial x}+v\dfrac{\partial\zeta}{\partial y}+\beta v\Big)-(\zeta+f_0)\Big(\dfrac{\partial\eta}{\partial t}+u\dfrac{\partial\eta}{\partial x}+v\dfrac{\partial\eta}{\partial y}\Big)=0$ (13.114).""")
note("N129 [B]", r"""
**Linearised** (step 4): $H\dfrac{\partial\zeta}{\partial t}+H\beta v-f_0\dfrac{\partial\eta}{\partial t}=0$ (13.115).""")
note("N130 [B]", r"""
**Geostrophic velocities** (steps 5–6): $u\simeq-\dfrac g{f_0}\dfrac{\partial\eta}{\partial y}$, $v\simeq\dfrac g{f_0}\dfrac{\partial\eta}{\partial x}$ (13.116),
hence $\zeta=\dfrac g{f_0}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}\Big)$.""")
note("N131 [B]", r"""
**The quasi-geostrophic vorticity equation** (the result):
$\dfrac{\partial}{\partial t}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}-\dfrac{f_0^2}{c^2}\eta\Big)+\beta\dfrac{\partial\eta}{\partial x}=0$, $c=\sqrt{gH}$ (13.117).""")
code(r"""
print("the quasi-geostrophic vorticity equation (13.117) follows from (13.94) ->", ch13.qg_vorticity_sympy()["ok"])   # the engine agrees: True
""", explain=r"""
**What does this show?** The engine repeats `D22` symbolically: expand, linearise, insert the geostrophic velocities.""")
D("D23", ref="13.118")
note("N134 [B]", r"""
**The zonal phase speed** (step 3 of `D23`): $c_x=\dfrac\omega k=-\dfrac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.119) — westward for every wave.""")
note("N132 [B]", r"""
**Circles, group velocity, maximum frequency** (steps 4–8): curves of constant $\omega$ are the circles
$\Big(k+\dfrac\beta{2\omega}\Big)^2+l^2=\Big(\dfrac\beta{2\omega}\Big)^2-\dfrac{f_0^2}{c^2}$; the group velocity is
$\mathbf c_g=\mathbf e_x\dfrac{\partial\omega}{\partial k}+\mathbf e_y\dfrac{\partial\omega}{\partial l}$ with (ours)
$c_{gx}=\dfrac{\beta(k^2-l^2-f_0^2/c^2)}{(k^2+l^2+f_0^2/c^2)^2}$, $c_{gy}=\dfrac{2\beta kl}{(k^2+l^2+f_0^2/c^2)^2}$; and the largest
frequency is $\omega_{max}=\beta c/(2f_0)$, at $kc/f_0=-1$.""")
note("N136 [B]", r"""
**On a mean current $U$** (step 9): $c_x=U-\dfrac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.120); the wave stands still for
$\lambda=2\pi\sqrt{U/\beta}$ (barotropic case).""")
trap("T8", r"""
With $\omega$ taken positive the zonal wavenumber is negative. The "maximum" phase speed is a maximum of magnitude. The group
velocity is westward for long waves and eastward for short ones; drawn on the constant-$\omega$ circles its arrows point
inward.""")
nb.worked_example("two Rossby waves by hand", r"""
Our 35° N value, rounded: $\beta=1.9\times10^{-11}$ m⁻¹ s⁻¹.

**(a) A barotropic atmospheric wave**, wavelength 6300 km ($k=-10^{-6}$ m⁻¹), $l=0$, $f_0^2/c^2$ negligible:

1. $\omega=-\beta k/k^2=\beta/\lvert k\rvert=1.9\times10^{-5}$ s⁻¹ → period 3.8 days.
2. $c_x=-\beta/k^2=-19$ m/s: westward at 19 m/s relative to the air.
3. In a westerly wind of 19 m/s it stands still — a stationary wave ($U=\beta/k^2$).

**(b) A long baroclinic ocean wave** with Rossby radius $\Lambda=43$ km:

4. $c_x\simeq-\beta\Lambda^2=-1.9\times10^{-11}\times1.85\times10^{9}=-0.035$ m/s ≈ −3 km per day.
5. An ocean 10 000 km wide is crossed in $10^{7}/0.035=2.9\times10^{8}$ s ≈ 9 years.""")
code(r"""
k_r = -k31                                                     # a 3100 km wave with k < 0, so that omega > 0 [rad/m]
c_bc = 3.6096                                                  # first-baroclinic speed N H/pi of our ocean [m/s]
for nm, c_ in (("external mode", c_ext), ("first baroclinic mode", c_bc)):   # two long-wave speeds
    w_ = GFD.rossby_omega(k_r, 0.0, beta35, f35, c_)           # Eq. (13.118) [rad/s]
    cx_ = GFD.rossby_phase_speed(k_r, 0.0, beta35, f35, c_)    # Eq. (13.119) [m/s]
    cgx_, _ = GFD.rossby_group_velocity(k_r, 0.0, beta35, f35, c_)   # group velocity (ours, D23) [m/s]
    print(f"{nm}: omega = {w_:.3e} 1/s (period {2*np.pi/w_/86400:.1f} days), c_x = {cx_:+.4f} m/s, c_gx = {cgx_:+.4f} m/s; 2 pi Lambda = {2*np.pi*c_/f35/1e3:.0f} km")   # crests west; energy east or west
wm = GFD.rossby_max_frequency(beta35, f35, c_bc)               # the largest frequency of baroclinic Rossby waves at 35 N
print(f"baroclinic: omega_max = {wm['omega_max']:.2e} 1/s -> shortest period {2*np.pi/wm['omega_max']/86400:.0f} days")   # nothing faster exists
print(f"long-wave speed -beta Lambda^2: {100*GFD.rossby_long_wave_speed(beta35, f35, c_bc):.2f} cm/s at 35 N, {100*GFD.rossby_long_wave_speed(beta12, f12, c_bc):.1f} cm/s at 12 N")   # non-dispersive and westward
print(f"a 10 000 km basin is crossed in {ch13.basin_crossing_time(1.0e7, inp['lat_rossby'], c_bc)/86400/365.25:.2f} years at 12 N and {ch13.basin_crossing_time(1.0e7, lat, c_bc)/86400/365.25:.2f} years at 35 N")   # the ocean's memory
print(f"stationary wavelength in a {inp['U_mean']:.0f} m/s westerly: {GFD.stationary_rossby_wavelength(inp['U_mean'], beta35)/1e3:.0f} km")   # Eq. (13.120)
print("circle of constant omega = 1e-6 1/s (external):", {k_: (v_ if isinstance(v_, bool) else f"{v_:.3e}") for k_, v_ in GFD.rossby_omega_circle(1.0e-6, beta35, f35, c_ext).items()})   # centre and radius
""", explain=r"""
1. `GFD.rossby_omega`, `GFD.rossby_phase_speed`, `GFD.rossby_group_velocity` are the dispersion relation, the phase speed and
   the group velocity; they return signed values for signed $k$.
2. The same 3100 km wave sends its energy **east** on the external mode (it is shorter than $2\pi\Lambda$) and **west** on the
   first baroclinic mode (it is longer than $2\pi\Lambda$); its crests go west in both.
3. `GFD.rossby_max_frequency`: baroclinic Rossby waves at 35° N have periods of half a year or more.
4. `GFD.rossby_long_wave_speed` is the long-wave limit $-\beta\Lambda^2$: centimetres per second, nine times faster at 12° than
   at 35°.
5. `ch13.basin_crossing_time`: a year at 12° N, nine years at 35° N — how long the ocean takes to hear about a change of wind.
6. `GFD.stationary_rossby_wavelength` is the wave that a westerly wind holds still, $2\pi\sqrt{U/\beta}$: about 6000 km in a
   17 m/s wind. This is the rigid-lid (barotropic) value; with a finite Rossby radius the condition is $U=\beta/(k^2+f_0^2/c^2)$
   and the stationary wave is slightly longer — the explainer below shows both.""")
note("N135 [B]", r"""
**Long waves:** $c_x\simeq-\dfrac{\beta c^2}{f_0^2}=-\beta\Lambda^2$ — non-dispersive and westward (the pair of numbers in the cell
above).""")
slip(12, r"an exercise answer for this long-wave speed that follows only from a round value of $\beta$",
     r"$c_x\simeq-\beta c^2/f_0^2$ with $\beta$ evaluated at the stated latitude, which is about a fifth smaller (a remark on rounding; the exercise is not reproduced).")
scratch(r"""
# From scratch: the dispersion relation typed out, and the group velocity by the complex-step trick
om = lambda k_, l_: -beta35*k_/(k_**2 + l_**2 + f35**2/c_bc**2)   # Eq. (13.118) for the first baroclinic mode
w_mine = om(k_r, 0.0)                                          # the frequency
h_cs = 1e-20                                                   # a tiny imaginary step
cg_cs = (np.imag(om(k_r + 1j*h_cs, 0.0))/h_cs, np.imag(om(k_r, 0.0 + 1j*h_cs))/h_cs)   # d omega/dk and d omega/dl
assert np.isclose(w_mine, GFD.rossby_omega(k_r, 0.0, beta35, f35, c_bc))                        # same frequency
assert np.allclose(cg_cs, GFD.rossby_group_velocity(k_r, 0.0, beta35, f35, c_bc), rtol=1e-12)   # same group velocity
print(f"✓ omega = {w_mine:.4e} 1/s; complex-step group velocity = ({cg_cs[0]:+.5f}, {cg_cs[1]:+.1e}) m/s agrees with the library")   # visible confirmation
""", r"""
The *complex-step derivative*: perturb the argument by a tiny imaginary amount $ih$, take the imaginary part of the result and
divide by $h$. There is no subtraction, so no round-off — the derivative comes out to machine precision. It reproduces the
closed-form group velocity of `D23`.""")
note("N133 [B]", r"""
**The dispersion diagram** (our version of the book's figure, `ch13.fig_rossby_dispersion`): on the left $\omega/\omega_{max}$
against $k\Lambda$ for $l=0$, with $\omega_{max}=\beta c/(2\lvert f_0\rvert)$ and $\Lambda=c/\lvert f_0\rvert$; on the right circles of
constant $\omega$ in the $(k\Lambda,\ l\Lambda)$ plane with group-velocity arrows.""")
fig(r"""
figRo = ch13.fig_rossby_dispersion()                           # omega/omega_max against k Lambda (l = 0), and circles of constant omega (35 N)
plt.show()
""",
    see="Left: the frequency (as a fraction of the largest possible one) against the zonal wavenumber — a curve that rises "
        "from zero, peaks at 1 where $k\\Lambda=-1$, and falls again. Right: circles of constant frequency, with the group "
        "velocity drawn as arrows.",
    read="Between the origin and the maximum (long waves) the curve is nearly a straight line: non-dispersive, energy westward. "
         "Beyond the maximum (short waves) the slope has the other sign: energy eastward. On the circles the arrows point inward.",
    change="…the Rossby radius were infinite (barotropic, rigid lid)? There would be no maximum: every wave would send its energy "
           "east.")
nb.plotly(r"""
lam_ax = np.geomspace(5.0e4, 5.0e7, 80)                        # wavelengths from 50 km to 50 000 km [m]

def f8(logLam):                                                # curves for one Rossby radius, Lambda = 10**logLam km
    Lam_ = 10**logLam*1e3                                      # the radius [m]
    c_ = Lam_*f35                                              # the long-wave speed that gives this radius at 35 N [m/s]
    cx_ = GFD.rossby_phase_speed(-2*np.pi/lam_ax, 0.0, beta35, f35, c_)        # Eq. (13.119)
    cg_ = GFD.rossby_group_velocity(-2*np.pi/lam_ax, 0.0, beta35, f35, c_)[0]  # c_gx
    unit = beta35*Lam_**2                                      # the long-wave speed beta Lambda^2 used as the unit
    return {"phase speed c_x / (βΛ²): always westward": (lam_ax/1e3, cx_/unit),
            "group velocity c_gx / (βΛ²): changes sign at λ = 2πΛ": (lam_ax/1e3, cg_/unit),
            "zero": (lam_ax/1e3, 0*lam_ax)}

figF8 = slider_figure(f8, "log10 of the Rossby radius Λ [km]", np.round(np.linspace(np.log10(20.0), np.log10(3000.0), 20), 3),
                      xlabel="wavelength λ [km]", ylabel="speed in units of βΛ²", title="Short waves send energy east, long waves west",
                      yrange=(-1.1, 0.2))                      # 20 slider positions, logarithmically spaced
figF8.update_xaxes(type="log", range=[np.log10(50.0), np.log10(5.0e4)])   # logarithmic wavelength axis
figF8.show()
""", explain=r"""
`f8` returns the phase speed and the zonal group velocity against wavelength for one Rossby radius, both in units of the
long-wave speed $\beta\Lambda^2$. The slider moves the radius from 20 km (an ocean's internal mode) to 3000 km (an external
mode).""")
nb.figure_notes(
    see="Phase speed (always negative: westward) and group velocity against wavelength. The group velocity is positive "
        "(eastward) for short waves, crosses zero at $\\lambda=2\\pi\\Lambda$, and joins the phase speed at $-1$ for long waves.",
    read="Sliding the radius shifts the crossing: for an internal mode almost every wave longer than a few hundred kilometres "
         "carries energy west; for the external mode even a 10 000 km wave carries it east.",
    change="…a mean eastward current were added? Both curves would shift up by $U$ (in these units), and the wave with $c_x=0$ "
           "would stand still.")
note("N137 [B]", r"""
**The cubic of `C09` again:** $\omega^3-c^2\omega(k^2+l^2)-f_0^2\omega-c^2\beta k=0$ (13.121). For $\omega\ll f$ its first term drops and
the relation of this block follows.""")
fig(r"""
lam_sw = np.geomspace(3.0e5, 2.0e7, 60)                        # wavelengths [m] (35 N: three real roots for all of them)
fig, ax = plt.subplots(figsize=(6.0, 3.6))
for nm, c_, col in (("external mode", c_ext, COLORS["accent"]), ("first baroclinic mode", c_bc, COLORS["teal"])):   # two modes
    slow = np.array([GFD.shallow_water_omega(2*np.pi/L_, 0.0, c_, f35, beta35)[1] for L_ in lam_sw])   # the slow root of the cubic
    ross = GFD.rossby_omega(2*np.pi/lam_sw, 0.0, beta35, f35, c_)   # Eq. (13.118) for the same waves
    ax.loglog(abs(ross)/f35, abs(slow/ross - 1), "o", ms=3, color=col, label=nm)   # relative difference against omega/f
ax.loglog([1e-4, 0.3], [1e-8, 0.09], "--", color=COLORS["muted"], label="(ω/f)² (upper bound)")   # a guide proportional to (omega/f)^2
ax.set(xlabel="ω / f", ylabel="relative difference", title="The Rossby formula is the cubic's slow root, to (ω/f)² or better")
ax.legend(fontsize=8)
plt.show()
""",
    see="The relative difference between the slow root of the full cubic and the Rossby formula, for a sweep of wavelengths on "
        "two modes, against $\\omega/f$.",
    read="The guide $(\\omega/f)^2$ is an upper bound. The baroclinic points sit on it. The external-mode points lie below it, "
         "because the error of dropping $\\omega^3$ is really $\\omega^2/(c^2K^2+f_0^2)$ and for these waves $c^2K^2\\gg f_0^2$; "
         "it reaches a few per cent only at the longest external waves.",
    change="…the latitude were low and the mode external? There the cubic may not even have three real roots (the caveat of `C09`).")
P("P329", "two-dimensional FFT wavenumber grids (np.fft.fftfreq, np.fft.fft2)", r"""
A field on a periodic grid is a sum of plane waves. `np.fft.fft2` returns their complex amplitudes; `2π·np.fft.fftfreq(n, d)`
gives the wavenumber that belongs to each array index (positive first, then negative). A linear equation with constant
coefficients moves each plane wave independently — multiply its amplitude by $e^{-i\omega t}$ and transform back: an exact
solution with no time stepping.""", code=r"""
k = 2*np.pi*np.fft.fftfreq(8, d=1.0)       # wavenumbers of an 8-point periodic grid
print(np.round(k, 2))                      # [0, 0.79, 1.57, 2.36, -3.14, -2.36, -1.57, -0.79]
""")
choice(r"""
the animation below is an exact spectral evolution on a periodic line (`SW.qg_linear_evolve_1d`): each Fourier mode of a packet
(`GFD.rossby_packet_spectrum`) is advanced with its own frequency from $\omega=-\frac{\beta k}{k^2+l^2+f_0^2/c^2}$ *(13.118)*. No time stepping, no
numerical error beyond round-off.""")
nb.animation(r"""
nfr = 30 if FAST else 36                                       # number of frames
t_pk = np.linspace(0.0, 7.0e7, nfr)                            # 2.2 years of evolution [s]
packs = []                                                     # one wave packet per row: (x grid, wavenumbers, amplitudes, speeds, label)
for lam0, nm in ((1.5e5, "short wave, 150 km"), (1.0e6, "long wave, 1000 km")):   # shorter and longer than 2 pi Lambda = 271 km
    k0 = 2*np.pi/lam0                                          # central wavenumber [rad/m]
    kk, a0 = GFD.rossby_packet_spectrum(k0, k0/10, 48)         # 48 modes with Gaussian amplitudes around k0
    xx = np.linspace(-9*lam0, 9*lam0, 600)                     # eighteen wavelengths of the line [m]
    cpx = GFD.rossby_phase_speed(k0, 0.0, beta35, f35, c_bc)   # speed of the crests [m/s]
    cgx = GFD.rossby_group_velocity(k0, 0.0, beta35, f35, c_bc)[0]   # speed of the envelope [m/s]
    packs.append((xx, kk, a0, cpx, cgx, nm))
fig, axes = plt.subplots(2, 1, figsize=(6.8, 4.8))
lines, marks = [], []
for ax, (xx, kk, a0, cpx, cgx, nm) in zip(axes, packs):        # set up the two rows
    (ln,) = ax.plot(xx/1e3, 0*xx, color=COLORS["accent"])      # the wave
    mk_g = ax.axvline(0, color=COLORS["orange"], lw=2)         # marker moving at the group velocity (the envelope)
    mk_p = ax.axvline(0, color=COLORS["muted"], ls="--")       # marker moving at the phase speed (a crest)
    ax.set(ylim=(-1.1, 1.1), xlim=(xx[0]/1e3, xx[-1]/1e3), ylabel="η (normalised)")
    ax.set_title(f"{nm}: crests {100*cpx:+.2f} cm/s (dashed), envelope {100*cgx:+.2f} cm/s (orange)", fontsize=8)
    lines.append(ln); marks.append((mk_g, mk_p))
axes[1].set_xlabel("x east [km]")

fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it: frames render much faster

def update(i):                                                 # frame i: both packets at time t_pk[i]
    for ln, (mk_g, mk_p), (xx, kk, a0, cpx, cgx, nm) in zip(lines, marks, packs):
        a_t = SW.qg_linear_evolve_1d(a0, kk, 0.0, t_pk[i], beta35, f35, c_bc)   # each mode's amplitude times exp(-i omega t)
        ln.set_ydata(np.real(np.exp(1j*np.outer(xx, kk)) @ a_t))                # the field: real part of the sum of modes
        mk_g.set_xdata([cgx*t_pk[i]/1e3]*2)                    # where the envelope should be
        mk_p.set_xdata([cpx*t_pk[i]/1e3]*2)                    # where a crest that started at x = 0 should be
    return lines

show_animation(animate(update, frames=nfr, fig=fig, interval=100), player="video", dpi=50)   # a smooth video
plt.close(fig)                                                 # do not show the last frame a second time
""", explain=r"""
1. Two packets on the first baroclinic mode at 35° N: one of 150 km waves (shorter than $2\pi\Lambda=271$ km), one of 1000 km
   waves (longer).
2. `SW.qg_linear_evolve_1d` multiplies each mode by $e^{-i\omega t}$; the field is the real part of the sum.
3. The dashed line moves at the phase speed, the orange line at the group velocity; the titles give both.""")
nb.figure_notes(
    see="Two wave packets over about two years. In both, the crests (follow the dashed line) march west. The envelope (orange "
        "line) of the short-wave packet creeps **east**; that of the long-wave packet moves **west**, almost with its crests.",
    read="Crests and energy part company: watch crests appear at the eastern edge of the short packet, travel through it and "
         "vanish at the western edge. At $\\lambda=2\\pi\\Lambda$ the envelope would stand still.",
    change="…the packets rode an eastward current of a few centimetres per second? Every speed shifts by $U$; the crests of one "
           "particular wavelength would stand still.")
fig(r"""
xe = np.linspace(-1.5e6, 1.5e6, 128, endpoint=False); ye = np.linspace(-7.5e5, 7.5e5, 64, endpoint=False)   # a periodic box 3000 km x 1500 km [m]
Xe, Ye = np.meshgrid(xe, ye)                                   # grids indexed [j, i]
eta_e0 = np.exp(-(Xe**2 + Ye**2)/(1.0e5)**2)                   # a Gaussian eddy 100 km in radius (height normalised to 1)
fig, axes = plt.subplots(1, 3, figsize=(9.6, 2.9), sharey=True)
for ax, yrs in zip(axes, (0.0, 0.5, 1.0)):                     # three times [years]
    eta_t = SW.qg_linear_evolve(eta_e0, xe, ye, yrs*3.156e7, beta35, f35, c_bc)   # exact evolution of Eq. (13.117) on the periodic box
    ax.pcolormesh(xe/1e3, ye/1e3, eta_t, cmap="RdBu_r", vmin=-1.0, vmax=1.0, shading="auto")   # the height field (red high, blue low)
    ax.axvline(0, color=COLORS["muted"], lw=0.6)
    ax.set(title=f"after {yrs:.1f} year", xlabel="x east [km]")
axes[0].set_ylabel("y north [km]")
plt.show()
""",
    see="A Gaussian eddy on the first baroclinic mode, evolved exactly by the linear quasi-geostrophic equation for one year "
        "on a box that is periodic in both directions (whatever leaves through one edge re-enters through the opposite one).",
    read="The eddy drifts west at a few kilometres per day (the long-wave speed), weakens, and leaves a train of shorter waves "
         "to its **east** — the short waves, whose energy travels eastward. The run is stopped before the eddy reaches the "
         "western edge of the periodic box.",
    change="…$\\beta=0$? The eddy would simply sit there: a steady geostrophic vortex.")
note("N138 [C]", r"""
**Where it fails.** This relation fails within a few degrees of the equator, where geostrophy itself breaks down. Equatorial
waves — the equatorial Kelvin wave and the equatorial radius $\sqrt{c/\beta}$ — are beyond the book; pointer: any text on
equatorial dynamics.""")
explainer("rossby_waves", "The crests go west — so how does the energy go east?",
          "a packet whose crests march west while its envelope moves east cannot be drawn still; marked columns show each one's "
          "vorticity change as it is displaced; a mean flow slows the crests to a standstill.",
          ["Play the short-wave preset and follow one crest and the envelope separately.",
           "Drag the wavelength through 2πΛ and watch the status flip from 'energy east' to 'energy west'.",
           "Add a mean flow until the crests stand still; read the wavelength.",
           "Click a column to see its displacement, βy and the vorticity it must have."])

# ---------------------------------------------------------------------------------------------------------------------
nb.section("13.16", "Barotropic Instability", intro=r"""
**What is this section about?** Chapter 11's question — when is a shear flow unstable? — asked again on the β-plane. The
answer is Rayleigh's criterion with one change: what must change sign is the gradient of *absolute* vorticity. (This section
continues block `C15`.)""")
nb.current_core = "C15"
note("N139 [B]", r"""
**Constant depth** turns potential-vorticity conservation into conservation of absolute vorticity:
$\Big(\dfrac{\partial}{\partial t}+\mathbf u\cdot\nabla\Big)(\zeta+f)=0$ (13.122). `SW.barotropic_vorticity_rhs` is its spectral
right-hand side (used in `C17`).""")
note("N140 [B]", r"""
**Linearised about a zonal current $U(y)$**, with $u'=-\partial\psi/\partial y$, $v'=\partial\psi/\partial x$, $\bar\zeta=-dU/dy$:
$\dfrac{\partial}{\partial t}(\nabla^2\psi)+U\dfrac{\partial}{\partial x}(\nabla^2\psi)+\Big(\beta-\dfrac{d^2U}{dy^2}\Big)\dfrac{\partial\psi}{\partial x}=0$ (13.123)
(steps 1–3 of `D24`).""")
nb.recap("R30", "Rayleigh's equation, now with β", r"""
For normal modes $\psi=\hat\psi(y)e^{ik(x-ct)}$:
$(U-c)\Big[\dfrac{d^2}{dy^2}-k^2\Big]\hat\psi+\Big[\beta-\dfrac{d^2U}{dy^2}\Big]\hat\psi=0$ — Chapter 11's Rayleigh equation with $-U''$
replaced by $\beta-U''$. ⚠️ Trap T14: the stream function has the opposite sign to Chapter 11's; the eigenvalue $c$ does not
care.""", where=r"Ch. 11, $(U-c)\big(\frac{d^2\phi}{dy^2}-k^2\phi\big)-\frac{d^2U}{dy^2}\phi=0$ (11.81)")
nb.current_core = "C15"
D("D24", ref="13.124")
note("N141 [B]", r"""
**The Rayleigh–Kuo criterion:** a necessary (not sufficient) condition for instability is that
$\dfrac{d}{dy}(\bar\zeta+f)=\beta-\dfrac{d^2U}{dy^2}$ (13.124) changes sign somewhere in the flow (the book says the analysis of
Chapter 11 carries over and does not show it; `D24` does).""")
code(r"""
y_j = np.linspace(-3.0e6, 3.0e6, 1201)                         # north-south distance from the jet axis [m]
for L_j in (4.0e5, 1.2e6):                                     # two jet half-widths: 400 km and 1200 km
    U_j = -14.0/np.cosh(y_j/L_j)**2                            # an easterly jet of 14 m/s at 12 N (ours)
    rk = GFD.rayleigh_kuo_criterion(y_j, U_j, beta12)          # does beta - U'' change sign?
    Upp = np.gradient(np.gradient(U_j, y_j), y_j)              # U'' by differences, for the hand-written test
    mine = bool(np.any(np.diff(np.sign(beta12 - Upp[2:-2])) != 0))   # the same test in one line
    assert mine == rk["changes_sign"]                          # same verdict as the library
    print(f"L = {L_j/1e3:.0f} km: largest U'' = 2 U0/L^2 = {2*14.0/L_j**2:.2e} against beta = {beta12:.2e} -> gradient changes sign: {rk['changes_sign']} ({'may be unstable' if rk['changes_sign'] else 'stable by the criterion'})")   # the verdict
print(f"widest easterly jet that can be unstable: L = sqrt(2 U0/beta) = {np.sqrt(2*14.0/beta12)/1e3:.0f} km; the same jet as a westerly needs L < {np.sqrt(2*14.0/(3*beta12))/1e3:.0f} km")   # easterly jets are easier to destabilise
print("beta = 0 (Chapter 11's criterion): inflection points at y =", np.round(ch11.rayleigh_criterion(y_j, U=-14.0/np.cosh(y_j/4.0e5)**2)["y_I"]/1e3), "km")   # Rayleigh's own test
Us = lambda yy: 1/np.cosh(yy)**2                               # non-dimensional westerly jet U = sech^2 y
Usp = lambda yy: -2*np.tanh(yy)/np.cosh(yy)**2                 # its first derivative
Uspp = lambda yy: (4*np.tanh(yy)**2 - 2/np.cosh(yy)**2)/np.cosh(yy)**2   # its second derivative
k_list = (1.4,) if FAST else (1.0, 1.4, 1.8)                   # wavenumbers in units of 1/L
b_list = (0.0, 0.3, 0.6, 0.7)                                  # beta in units of U0/L^2
growth = {k_: [ch13.rayleigh_kuo_eigs(k_, Us, Usp, Uspp, b_, bc="decay", y_max=16, parity="even")["growth_rate"] for b_ in b_list] for k_ in k_list}   # k c_i for each pair
print("growth rate k c_i [U0/L] of the sech^2 westerly jet at k = 1.4, beta = 0, 0.3, 0.6, 0.7:", np.round(growth[1.4], 4))   # beta stabilises
""", explain=r"""
1. `GFD.rayleigh_kuo_criterion(y, U, beta)` tests whether $\beta-U''$ changes sign; the one-line hand test gives the same
   verdict (the `assert`).
2. A 400 km wide easterly jet of 14 m/s at 12° N passes the test for possible instability; a 1200 km wide one does not.
3. Because $\beta$ adds to $-U''$, an easterly jet is easier to destabilise than a westerly one of the same shape.
4. `ch13.rayleigh_kuo_eigs` solves the Rayleigh equation with β numerically (Chapter 11's contour solver, with β passed as its
   own argument) for the non-dimensional westerly jet $U=\mathrm{sech}^2y$, where $\mathrm{sech}\,y\equiv1/\cosh y$ is a smooth bump
   of height 1 and width about 1 (with $\cosh y=(e^y+e^{-y})/2$): the growth rate falls as β rises and no growing mode
   is found at β = 0.7, beyond the largest $U''$ of this jet (2/3) — as the criterion requires.""")
note("N142 [B]", r"""
**The profile and its vorticity** (our version of the book's figure): $U(y)$, $\bar\zeta$, $f$ and $\bar\zeta+f$ for the 400 km
jet, with the growth rates computed above.""")
fig(r"""
U_4 = -14.0/np.cosh(y_j/4.0e5)**2                              # the 400 km easterly jet again [m/s]
zeta_b = -np.gradient(U_4, y_j)                                # its relative vorticity, zeta_bar = -dU/dy [1/s]
f_y = f12 + beta12*y_j                                         # planetary vorticity on the beta-plane [1/s]
fig, (a, b) = plt.subplots(1, 2, figsize=(9.4, 3.8))
a.plot(U_4/14.0, y_j/1e3, color=COLORS["accent"], label="U / 14 m/s")                       # the jet
a.plot(f_y*1e5, y_j/1e3, color=COLORS["teal"], label="f [10⁻⁵ s⁻¹]")                        # planetary vorticity
a.plot((zeta_b + f_y)*1e5, y_j/1e3, "--", color=COLORS["amber"], label="ζ̄ + f [10⁻⁵ s⁻¹]")  # absolute vorticity (amber)
a.set(ylim=(-1500, 1500), xlim=(-2, 8), xlabel="(see legend)", ylabel="y north [km]", title="Where ζ̄ + f turns back, its gradient changes sign")
a.legend(fontsize=7)
for k_, col in zip(growth, (COLORS["accent"], COLORS["teal"], COLORS["rose"])):   # growth rate against beta
    b.plot(b_list, growth[k_], "o-", color=col, label=f"k L = {k_}")
b.set(xlabel="β L² / U₀", ylabel="growth rate k c_i [U₀/L]", title="β stabilises the westerly sech² jet")
b.legend(fontsize=8)
plt.show()
print(f"for scale: our 400 km, 14 m/s jet at 12 N has beta L^2/U0 = {beta12*(4.0e5)**2/14.0:.2f}")   # where the real jet sits on the right panel
""",
    see="Left: the easterly jet, the planetary vorticity $f$ (a straight line) and the absolute vorticity $\\bar\\zeta+f$ "
        "(dashed). Right: computed growth rates of the $\\mathrm{sech}^2$ westerly jet against β.",
    read="Where the dashed curve turns back, its gradient $\\beta-U''$ changes sign: there, fluid columns can exchange places "
         "without violating conservation of absolute vorticity, and a wave can feed on the jet. On the right, growth weakens "
         "with β and stops near $\\beta L^2/U_0=2/3$.",
    change="…the jet were wider? The curvature $U''$ shrinks like $1/L^2$; once it is everywhere below β the dashed curve is monotonic and the "
           "flow is stable.")
whatif(r"""
…the flow has no horizontal shear at all, only a vertical one — the thermal wind of `C03`? There is no inflection point and
Ri is far above ¼, yet it is unstable: baroclinic instability (`C16`).""")


# =====================================================================================================================
# A.17  §13.17 — C16 the Eady problem
# =====================================================================================================================
nb.section("13.17", "Baroclinic Instability", intro=r"""
**What is this section about?** Why mid-latitude weather exists. A wind that increases with height over sloping density
surfaces is unstable to waves a few thousand kilometres long; they grow by carrying warm air poleward and upward and cold air
equatorward and downward. ⚠️ Here $z=0$ is the lower lid and $z=H$ the upper one; $\alpha$ is a scaled wavenumber, not the
expansion coefficient; $\Lambda=NH/f$ has no $\pi$.""")
core("C16", r"The Eady problem: $c=\dfrac{U_0}{2}\pm\dfrac{U_0}{\alpha H}\sqrt{\Big(\dfrac{\alpha H}{2}-\tanh\dfrac{\alpha H}{2}\Big)\Big(\dfrac{\alpha H}{2}-\coth\dfrac{\alpha H}{2}\Big)}$ (13.141)",
     "The jet stream has no inflection point and a Richardson number far above ¼ — so where do storms come from?")
problem(r"""
The tropics are warm, the poles cold, and in between the surfaces of constant density tilt — a little, about one part in five
hundred. That tilt is stored energy: if cold air could slide down the slope under warm air sliding up it, the centre of mass
would fall. Rotation forbids the direct slide (`C03`: the tilted state is in balance). But a wave of the right size can do it
sideways, and once it starts, it feeds itself. The waves are the highs and lows of the weather map; their size, about 4600 km for our
inputs, and their growth time, a day or two, come out of one formula. *The problem solved here is known as the
Eady problem.*""")
nb.md(r"""
> ⚠️ **slip #13 — the book credits** the theory of this instability to one named author, whose first name is misspelt
> **; only the spelling is asserted here.** Who should be credited is not asserted (we have not read the original papers
> first-hand). We say only that the problem solved in this section, with the phase speed $c=\dfrac{U_0}{2}\pm\dots$ of the
> heading, is known as the Eady problem.""")
idea(words=r"""
Two lids, and on each lid a wave that exists because temperature varies along it (an "edge wave"). The lower one, left alone,
drifts east slowly; the upper one sits in a fast eastward wind but itself propagates west relative to it. If they can feel each
other through the layer, they lock: each one's flow strengthens the other. They can feel each other only if the wavelength is
long compared with the layer's Rossby radius.""")
note("N143 [B]", r"""
**The basic state** (`ch13.eady_basic_state`): density surfaces sloping up toward the pole; by the thermal wind (`C03`) the
eastward flow increases with height.""")
code(r"""
N_a, H_a, U0_a = inp["atm_N"], inp["atm_H"], inp["atm_U0"]     # our atmosphere: N = 1.1e-2 1/s, H = 9 km, U0 = 27 m/s at the upper lid
bs = ch13.eady_basic_state(np.linspace(-1.0e6, 1.0e6, 5), np.linspace(0.0, H_a, 5), N=N_a, f=f35, H=H_a, U0=U0_a, rho0=1.2)   # density rho(z, y), wind U(z), slope
print(f"slope of the density surfaces f U0/(N^2 H) = {bs['slope']:.2e} (about 1 in {1/bs['slope']:.0f}); Richardson number N^2 H^2/U0^2 = {(N_a*H_a/U0_a)**2:.1f}")   # gentle slope, large Ri
print("wind at the five levels [m/s]:", np.round(bs["U"], 2))   # U = U0 z/H
""", explain=r"""
**What does this show?** The density surfaces slope by about 1 in 480, the wind rises linearly from 0 to 27 m/s, and the
Richardson number is about 13 — far above the ¼ that Chapter 11's shear instability needs. Whatever makes storms, it is not
that.""")
P("P332", "available potential energy and the wedge of sloping convection", r"""
Only part of a fluid's potential energy can ever be released: the part that would be freed by flattening the density surfaces.
A parcel exchange releases energy only if the heavier parcel ends up lower *and* was the denser one at the same level: its path
must lie inside the wedge between the horizontal and the sloping density surface. Steeper than the surface: it costs energy
(ordinary static stability). Flatter than horizontal: impossible.""", code=r"""
slope_rho = 2.07e-3                               # slope of our density surfaces (computed above)
for s in (0.5, 1.0, 1.5):                         # slope of the exchange path / slope of the density surface
    print(s, "releases energy" if 0 < s < 1 else "does not")
""")
fig(r"""
figW = draw_eady_wedge()                                       # our sketch: the wedge between the horizontal and a density surface
plt.show()
""",
    see="Sloping density surfaces (blue lines; light fluid above, dense below), the shaded wedge between the horizontal and one "
        "of those surfaces, and two arrows along a path inside the wedge: a light parcel rising toward the dense side, a dense "
        "parcel sinking toward the light side.",
    read="Both arrows are flatter than the density surfaces but not level: such an exchange puts light fluid higher and dense "
         "fluid lower, so the centre of mass falls. Growing baroclinic waves move air along such paths: poleward and up, "
         "equatorward and down.",
    change="…the density surfaces were level? The wedge closes: nothing to release, no instability.")
note("N144 [C]", r"""
**The starting set**, for the total flow on an f-plane, hydrostatic and inviscid:
$\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$,
$\dfrac{\partial v}{\partial t}+u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$,
$0=-\dfrac{\partial p}{\partial z}-\rho g$, $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$,
$\dfrac{\partial\rho}{\partial t}+u\dfrac{\partial\rho}{\partial x}+v\dfrac{\partial\rho}{\partial y}+w\dfrac{\partial\rho}{\partial z}=0$ (13.125)
(⚠️ trap T15: no $w\,\partial u/\partial z$ in the momentum equations).""")
note("N145 [C]", r"""
**Basic state plus perturbation** (13.126): $u=U(z)+u'$, $v=v'$, $w=w'$, $\rho=\bar\rho(y,z)+\rho'$, $p=\bar p(y,z)+p'$.""")
note("N146 [C]", r"""
**The basic state is geostrophic and hydrostatic** (13.127): $fU=-\dfrac1{\rho_0}\dfrac{\partial\bar p}{\partial y}$,
$0=-\dfrac{\partial\bar p}{\partial z}-\bar\rho g$.""")
P("P330", "non-trivial solutions of a homogeneous 2 × 2 system: the determinant must vanish", r"""
Two equations $a\,A+b\,B=0$ and $c\,A+d\,B=0$ always have the solution $A=B=0$. They have another one only if the two equations
say the same thing, i.e. the determinant $ad-bc$ is zero. When the coefficients contain an unknown (a frequency, a phase
speed), "determinant = 0" is the equation that fixes it.""", code=r"""
c = sp.symbols("c")                               # an unknown, like the phase speed
M = sp.Matrix([[c, -1], [-4, c]])                 # a toy system whose coefficients contain it
print("values of c with a non-trivial solution:", sp.solve(M.det(), c))   # [-2, 2]
""")
P("P331", "hyperbolic half-angle identities and coth", r"""
The function $\coth x=\cosh x/\sinh x=1/\tanh x$ is large near 0 and tends to 1 for large $x$. Two identities turn products at
the half argument into functions of the full argument: $2\sinh x\cosh x=\sinh2x$ and $\cosh^2x+\sinh^2x=\cosh2x$; hence
$\tanh x+\coth x=2\coth2x$.""", code=r"""
x = 0.8                                                    # any positive number
print("tanh x + coth x =", round(np.tanh(x) + 1/np.tanh(x), 4), "; 2 coth 2x =", round(2/np.tanh(2*x), 4))   # equal
print("x = coth x at x =", round(brentq(lambda x: x - 1/np.tanh(x), 0.5, 3), 5))   # 1.19968
""")
remind("C16")
D("D25", ref="13.136")
note("N147 [B]", r"""
**The thermal wind of the basic state** (step 3 of `D25`): $\dfrac{dU}{dz}=\dfrac g{f\rho_0}\dfrac{\partial\bar\rho}{\partial y}$ (13.128), so
$U=\dfrac{U_0z}H$.""")
note("N148 [C]", r"""
**Vorticity with stretching** (step 4): $\dfrac{\partial\zeta}{\partial t}+u\dfrac{\partial\zeta}{\partial x}+v\dfrac{\partial\zeta}{\partial y}-(\zeta+f)\dfrac{\partial w}{\partial z}=0$ (13.129).""")
note("N149 [B]", r"""
**Its linear form** (step 5): $\dfrac{\partial\zeta'}{\partial t}+U\dfrac{\partial\zeta'}{\partial x}-f\dfrac{\partial w'}{\partial z}=0$ (13.130).""")
note("N150 [C]", r"""
**Geostrophic perturbations** (step 6): $u'\simeq-\dfrac1{\rho_0f}\dfrac{\partial p'}{\partial y}$, $v'\simeq\dfrac1{\rho_0f}\dfrac{\partial p'}{\partial x}$ (13.131).""")
note("N151 [B]", r"""
**Vorticity from pressure** (step 7): $\zeta'=\dfrac1{\rho_0f}\nabla_H^2p'$ (13.132).""")
note("N152 [B]", r"""
**The linear density equation** (step 8): $\dfrac{\partial\rho'}{\partial t}+U\dfrac{\partial\rho'}{\partial x}+v'\dfrac{\partial\bar\rho}{\partial y}-\dfrac{\rho_0N^2w'}{g}=0$ (13.133).""")
note("N153 [C]", r"""
**Hydrostatic perturbations** (step 9): $0=-\dfrac{\partial p'}{\partial z}-\rho'g$ (13.134).""")
note("N154 [B]", r"""
**The vertical velocity from the pressure** (step 10):
$w'=-\dfrac1{\rho_0N^2}\Big[\Big(\dfrac\partial{\partial t}+U\dfrac\partial{\partial x}\Big)\dfrac{\partial p'}{\partial z}-\dfrac{dU}{dz}\dfrac{\partial p'}{\partial x}\Big]$ (13.135)
(`GFD.eady_vertical_velocity`).""")
note("N155 [B]", r"""
**The perturbation equation** (the result of `D25`):
$\Big(\dfrac\partial{\partial t}+U\dfrac\partial{\partial x}\Big)\Big[\nabla_H^2p'+\dfrac{f^2}{N^2}\dfrac{\partial^2p'}{\partial z^2}\Big]=0$ (13.136). The bracket is the
perturbation potential vorticity; it is carried by the basic flow and is zero in the interior for a normal mode — so
everything happens at the two lids.""")
code(r"""
print("the Eady perturbation equation (13.136) follows from the set (13.125) ->", ch13.eady_qg_sympy()["ok"])   # the engine agrees: True
""", explain=r"""
**What does this show?** The engine eliminates $w'$ between the vorticity and density equations and is left with the
potential-vorticity equation, as in `D25`.""")
D("D26", ref="13.141")
note("N156 [C]", r"""
**The wave solution** (step 1 of `D26`): $p'=\hat p(z)e^{i(kx+ly-\omega t)}$ (13.137).""")
note("N157 [B]", r"""
**Vertical structure** (step 2): $\dfrac{d^2\hat p}{dz^2}-\alpha^2\hat p=0$ (13.138).""")
note("N158 [B]", r"""
**The scaled wavenumber** (step 3): $\alpha^2\equiv\dfrac{N^2}{f^2}(k^2+l^2)$ (13.139) (`GFD.eady_alpha`).""")
note("N159 [C]", r"""
**Solution about mid-depth** (step 4): $\hat p=A\cosh\alpha\Big(z-\dfrac H2\Big)+B\sinh\alpha\Big(z-\dfrac H2\Big)$ (13.140).""")
note("N160 [B]", r"""
**The lids** (steps 5–8). No flow through them, $w'=0$, gives
$\dfrac{\partial^2p'}{\partial t\,\partial z}-\dfrac{U_0}H\dfrac{\partial p'}{\partial x}=0$ at $z=0$ and
$\dfrac{\partial^2p'}{\partial t\,\partial z}-\dfrac{U_0}H\dfrac{\partial p'}{\partial x}+U_0\dfrac{\partial^2p'}{\partial x\,\partial z}=0$ at $z=H$, and with the
solution above the 2 × 2 system
$A\Big[\alpha c\sinh\dfrac{\alpha H}2-\dfrac{U_0}H\cosh\dfrac{\alpha H}2\Big]+B\Big[-\alpha c\cosh\dfrac{\alpha H}2+\dfrac{U_0}H\sinh\dfrac{\alpha H}2\Big]=0$,
$A\Big[\alpha(U_0-c)\sinh\dfrac{\alpha H}2-\dfrac{U_0}H\cosh\dfrac{\alpha H}2\Big]+B\Big[\alpha(U_0-c)\cosh\dfrac{\alpha H}2-\dfrac{U_0}H\sinh\dfrac{\alpha H}2\Big]=0$
(`ch13.eady_matrix(c, alphaH, U0, H)`).""")
nb.worked_example("storm size and growth time by hand", r"""
Our inputs, rounded: $f=8.4\times10^{-5}$ s⁻¹ (35° N), $N=1.1\times10^{-2}$ s⁻¹, $H=9$ km, $U_0=27$ m/s.

1. Eady radius $\Lambda_E=NH/f=1.1\times10^{-2}\times9000/8.4\times10^{-5}=1.18\times10^{6}$ m ≈ 1180 km.
2. Fastest wave (ours, computed): $\alpha H=1.6061$, so $k=1.6061/\Lambda_E$ and the wavelength is $2\pi/k=3.9120\times1180$ km ≈ 4600 km.
3. Growth rate $\sigma=0.30982\times fU_0/(NH)=0.30982\times8.4\times10^{-5}\times27/99=7.1\times10^{-6}$ s⁻¹.
4. e-folding time $1/\sigma=1.4\times10^{5}$ s ≈ 1.6 days.
5. Shortest unstable wave: $2.6187\times1180$ km ≈ 3100 km; anything shorter is two neutral edge waves.""")
D("D27", ref="13.142")
note("N161 [B]", r"""
**The cut-off.** The marginal condition is $\dfrac{\alpha_cH}2=\coth\Big(\dfrac{\alpha_cH}2\Big)$; the flow is unstable for
$\alpha H<\alpha_cH$, i.e. $\dfrac{HN}f<\dfrac{\alpha_cH}{\sqrt{k^2+l^2}}$, and for $l=0$ $\dfrac{HN}f<\dfrac{\alpha_cH}k$ (13.142) (the book prints
$\alpha_cH$ as a rounded number). With $\Lambda\equiv\dfrac{HN}f$ the unstable wavelengths are $\lambda>(2\pi/\alpha_cH)\,\Lambda$.""")
trap("T11", r"""
This $\Lambda=NH/f$ has **no** $\pi$: it is about three times the internal radius $NH/(\pi f)$ of `C12`.""")
note("N162 [B]", r"""
**Not in the book — ours** (the book stops at the wavelength): the maximum growth rate
$\sigma_{max}=0.30982\,\dfrac fN\dfrac{dU}{dz}$ at $\alpha H=1.6061$ ($l=0$), and its e-folding time. The five-digit constants are
ours, computed by `GFD.eady_fastest()`; the values 0.3098 and 1.606 were checked first-hand against
K. Emanuel's MIT OpenCourseWare notes for course 12.803, Lecture 19 (`reference/ch13/SOURCES.md`).""")
code(r"""
fast = GFD.eady_fastest()                                      # the fastest-growing wave (l = 0), by a bounded maximisation
fac = GFD.eady_factors(fast["alphaH"])                         # the two brackets under the root of Eq. (13.141)
print(f"cut-off alpha_c H = {GFD.eady_critical():.5f}; fastest alpha H = {fast['alphaH']:.4f}, growth rate = {fast['sigma_nd']:.5f} f U0/(N H), c_i = {fast['ci_over_U0']:.4f} U0")   # ours, computed
print(f"factors at the fastest wave: (x - tanh x) = {fac['tanh_factor']:+.4f}, (x - coth x) = {fac['coth_factor']:+.4f} -> negative product, complex c")   # why it grows
c1 = GFD.eady_phase_speed(1.0, U0_a); c3 = GFD.eady_phase_speed(3.0, U0_a)   # Eq. (13.141) at alpha H = 1 (unstable) and 3 (neutral)
print(f"alpha H = 1: c = {c1.real:.2f} + {c1.imag:.2f} i m/s (growing);  alpha H = 3: c = {c3.real/U0_a:.4f} U0, real (neutral)")   # complex vs real
Lam_E = GFD.rossby_radius_internal(N_a, H_a, f35, with_pi=False)             # the Eady radius N H/f (no pi) [m]
sig = GFD.eady_max_growth_rate(f=f35, N=N_a, dUdz=U0_a/H_a)                  # sigma_max = 0.30982 (f/N) dU/dz [1/s]
t_e = GFD.eady_time_scale(f=f35, N=N_a, dUdz=U0_a/H_a)                       # its e-folding time [s]
print(f"atmosphere, 35 N: Lambda_E = {Lam_E/1e3:.0f} km, sigma_max = {sig:.3e} 1/s, e-folding {t_e/86400:.3f} days, fastest wavelength {2*np.pi/fast['alphaH']*Lam_E/1e3:.0f} km, cut-off {2*np.pi/GFD.eady_critical()*Lam_E/1e3:.0f} km")   # weather systems
print(f"check: growth rate from k = 1.6061/Lambda_E: {GFD.eady_growth_rate(fast['alphaH']/Lam_E, 0.0, N_a, f35, H_a, U0_a):.3e} 1/s")   # the same number
Lam_O = GFD.rossby_radius_internal(5.0e-3, 1000.0, f35, with_pi=False)       # an ocean twin: N = 5e-3 1/s over the top 1000 m
print(f"ocean twin (U0 = 0.1 m/s): Lambda_E = {Lam_O/1e3:.1f} km, e-folding {GFD.eady_time_scale(f=f35, N=5.0e-3, dUdz=0.1/1000.0)/86400:.1f} days, fastest wavelength {2*np.pi/fast['alphaH']*Lam_O/1e3:.0f} km")   # mesoscale eddies
""", explain=r"""
1. `GFD.eady_critical()` solves $x=\coth x$; `GFD.eady_fastest()` maximises the growth rate (ours, computed).
2. `GFD.eady_factors` shows why there is growth: below the cut-off the second bracket is negative, the root is imaginary and
   $c$ is complex.
3. `GFD.eady_phase_speed(alphaH, U0)` is $c$ from the formula of this block: complex at $\alpha H=1$, real at $\alpha H=3$.
4. `GFD.eady_max_growth_rate` and `GFD.eady_time_scale` are called **by keyword** (`f=`, `N=`, `dUdz=`): sibling functions
   take $N$ before $f$, and swapping two positional numbers would give a wrong answer without any error.
5. For our atmosphere: an Eady radius of about 1200 km, storms about 4600 km long that grow by a factor $e$ in 1.6 days.
6. For an ocean thermocline: a radius of 60 km, eddies about 230 km across that take three weeks — the ocean's "weather".""")
scratch(r"""
# From scratch: the determinant of the 2 x 2 system is a quadratic in c — three evaluations give its coefficients
detM = lambda c_: np.linalg.det(ch13.eady_matrix(c_, 1.4, U0_a, H_a))   # the determinant for a trial phase speed c [m/s] at alpha H = 1.4
d0, d1, d2 = detM(0.0), detM(1.0), detM(2.0)                   # its value at three trial speeds
qa = (d2 - 2*d1 + d0)/2; qb = d1 - d0 - qa; qc = d0            # the quadratic q_a c^2 + q_b c + q_c through the three values
roots = np.roots([qa, qb, qc])                                 # its two zeros: the two phase speeds
c_grow = roots[np.argmax(roots.imag)]                          # the one with positive imaginary part (the growing wave)
assert np.isclose(c_grow, GFD.eady_phase_speed(1.4, U0_a))     # same as Eq. (13.141)
k14 = 1.4/Lam_E                                                # the wavenumber with alpha H = 1.4 [rad/m]
c_cheb = ch13.eady_numeric_eigs(k14, 0.0, N_a, f35, H_a, U0_a)["c"]   # an independent route: Chebyshev eigen-solve of Eq. (13.136) with the lid conditions
print(f"✓ determinant route: c = ({c_grow.real/U0_a:.4f} + {c_grow.imag/U0_a:.4f} i) U0; formula: same; Chebyshev route differs by {abs(c_cheb - c_grow)/U0_a:.1e} U0")   # visible confirmation
""", r"""
"The determinant must vanish" (primer above), done numerically: the determinant of `ch13.eady_matrix` is a quadratic in $c$,
so three evaluations fix it and `np.roots` gives its zeros. The growing root equals the closed form; a completely independent
numerical solution of the differential equation (`ch13.eady_numeric_eigs`, on Chebyshev points) agrees as well.""")
fig(r"""
aH = np.linspace(0.05, 4.0, 400)                               # the scaled wavenumber alpha H
x_h = aH/2                                                     # half of it
growth_nd = np.array([aH_*GFD.eady_phase_speed(aH_, 1.0).imag for aH_ in aH])   # k c_i in units of f U0/(N H): alpha H times c_i/U0
fig, (a, b) = plt.subplots(2, 1, figsize=(6.6, 5.4), sharex=True)
a.plot(aH, 1/np.tanh(x_h), color=COLORS["rose"], label="coth(αH/2)")       # coth x
a.plot(aH, np.tanh(x_h), color=COLORS["teal"], label="tanh(αH/2)")         # tanh x
a.plot(aH, x_h, color=COLORS["ink"], label="αH/2")                         # the line x
a.axvline(GFD.eady_critical(), color=COLORS["muted"], ls="--")             # the cut-off, where x = coth x
a.set(ylim=(0, 3), ylabel="value", title=f"Unstable where αH/2 lies between tanh and coth: αH < {GFD.eady_critical():.4f}")
a.legend(fontsize=8)
b.plot(aH, growth_nd, color=COLORS["accent"])                              # the growth rate
b.plot(fast["alphaH"], fast["sigma_nd"], "o", color=COLORS["orange"])      # its maximum
b.axvline(GFD.eady_critical(), color=COLORS["muted"], ls="--")
b.set(xlabel="scaled wavenumber αH", ylabel="growth rate $kc_i$ [f U₀/(N H)]", title=f"Fastest growth {fast['sigma_nd']:.5f} at αH = {fast['alphaH']:.4f}")
plt.show()
""",
    see="Upper: $\\coth x$, $\\tanh x$ and the line $x=\\alpha H/2$. Lower: the growth rate against the scaled wavenumber, with "
        "its maximum marked; the dashed line is the cut-off.",
    read="The radicand of the phase-speed formula is negative exactly where the straight line lies between the two curves — to "
         "the left of the crossing $x=\\coth x$. There the growth rate rises from zero, peaks, and falls back to zero at the "
         "cut-off.",
    change="…$l\\neq0$? The same curve read at $\\alpha=(N/f)\\sqrt{k^2+l^2}$; the growth rate carries the factor $k/K$ and is "
           "smaller.")
confusion(r"""
The growth rate is $kc_i$, not $c_i$. The imaginary part $c_i$ is largest as $k\to0$, where nothing grows.""")
nb.plotly(r"""
lam_e = np.linspace(1.0e6, 1.2e7, 70)                          # wavelengths from 1000 km to 12 000 km [m]

def f9(U0_):                                                   # growth rate against wavelength for one wind difference U0 [m/s]
    sig_ = np.array([GFD.eady_growth_rate(2*np.pi/L_, 0.0, N_a, f35, H_a, U0_) for L_ in lam_e])*86400   # [1/day]
    return {"growth rate k c_i (atmosphere: N = 1.1e-2 1/s, H = 9 km, 35° N)": (lam_e/1e3, sig_),
            "cut-off wavelength 2.62 Λ_E (shorter waves are neutral)": ([2*np.pi/GFD.eady_critical()*Lam_E/1e3]*2, [0.0, 1.2])}

figF9 = slider_figure(f9, "wind difference U₀ between the lids", np.arange(5.0, 50.1, 2.5), unit="m/s", xlabel="wavelength [km]",
                      ylabel="growth rate [1/day]", title="More shear: faster storms of the same size", yrange=(0, 1.25))   # 19 slider positions
figF9.show()
for U0_ in (10.0, 27.0, 50.0):                                 # three shears
    print(f"U0 = {U0_:4.0f} m/s: e-folding time {GFD.eady_time_scale(f=f35, N=N_a, dUdz=U0_/H_a)/86400:.2f} days")   # inversely proportional to the shear
""", explain=r"""
`f9` evaluates `GFD.eady_growth_rate` over a range of wavelengths for one value of the wind difference between the lids; the
vertical line is the cut-off. The printed lines give the e-folding time for three shears.""")
nb.figure_notes(
    see="The growth rate (per day) against wavelength. Nothing grows below the cut-off near 3100 km; the maximum sits near "
        "4600 km for every shear.",
    read="The shear sets how fast storms grow (the curve scales with $U_0$); the stratification, depth and latitude set how big "
         "they are (through $\\Lambda_E=NH/f$).",
    change="…$N$ doubled? The radius $\\Lambda_E$ doubles — storms twice as long — and the growth rate halves.")
nb.live(r"""
def eady_live(N=1.1e-2, H=9000.0, U0=27.0, lat_deg=35.0):      # free buoyancy frequency [1/s], depth [m], wind difference [m/s], latitude [degrees]
    f_ = GFD.coriolis_parameter(np.deg2rad(lat_deg))           # Coriolis parameter there
    Lam_ = GFD.rossby_radius_internal(N, H, f_, with_pi=False) # Eady radius N H/f [m]
    lam_ = np.linspace(0.5, 8.0, 200)*Lam_                     # wavelengths from 0.5 to 8 Eady radii
    sig_ = np.array([GFD.eady_growth_rate(2*np.pi/L_, 0.0, N, f_, H, U0) for L_ in lam_])*86400   # growth rate [1/day]
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    ax.plot(lam_/1e3, sig_, color=COLORS["accent"])            # the growth curve
    ax.set(xlabel="wavelength [km]", ylabel="growth rate [1/day]", title=f"Λ_E = {Lam_/1e3:.0f} km, e-folding {GFD.eady_time_scale(f=f_, N=N, dUdz=U0/H)/86400:.2f} days")
    plt.show()

_ = live(eady_live, N=(0.004, 0.02, 0.001), H=(1000.0, 12000.0, 500.0), U0=(1.0, 50.0, 1.0), lat_deg=(15.0, 75.0, 5.0))   # four sliders
""", explain=r"""
Four free sliders — stratification, depth, wind difference, latitude — and the growth curve with the Eady radius and the
e-folding time in its title. On the web page the slider figure above carries the same idea.""")
note("N163 [B]", r"""
**Energetics** (stated; the integrations are Chapter 11's energy-equation moves, whose result there was
$\frac d{dt}\int\frac12u_i^2\,dV=-\int u_iu_j\frac{\partial U_i}{\partial x_j}\,dV-\Lambda$ (11.88) with $\Lambda$ the dissipation). Here the perturbation
kinetic energy $K_E\equiv\dfrac{\rho_0}2\displaystyle\int(u'^2+v'^2)\,dx\,dy\,dz$ grows through the buoyancy flux alone,
$\dfrac{dK_E}{dt}=-g\displaystyle\int w'\rho'\,dx\,dy\,dz$ — not through the Reynolds stress against the shear.""")
fig(r"""
z_e = np.linspace(0.0, H_a, 121)                               # heights between the lids [m]
fl = GFD.eady_fluxes(z_e, fast["alphaH"], U0_a, H_a, N_a, f35, 1.2)   # wave-averaged fluxes of the fastest mode (unit pressure amplitude)
md = GFD.eady_mode(z_e, fast["alphaH"], U0_a, H_a)             # its vertical structure: amplitude and phase of p-hat
fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.6), sharey=True)
axes[0].plot(fl["w_rho"]/abs(fl["w_rho"]).max(), z_e/1e3, color=COLORS["blue"])   # mean of w' rho' (normalised)
axes[0].set(xlabel=r"$\overline{w'\rho'}$ (normalised)", ylabel="z [km]", title="Light air rises")
axes[1].plot(fl["v_rho"]/abs(fl["v_rho"]).max(), z_e/1e3, color=COLORS["blue"])   # mean of v' rho' (normalised)
axes[1].set(xlabel=r"$\overline{v'\rho'}$ (normalised)", title="Warm air goes poleward", xlim=(-1.2, 0.2))
axes[2].plot(np.degrees(md["phase"]), z_e/1e3, color=COLORS["accent"])            # phase of the pressure wave against height
axes[2].set(xlabel="phase of $\\hat p$ [degrees]", title="Troughs tilt west with height")
plt.show()
print(f"w'rho' <= 0 at every height: {bool((fl['w_rho'] <= 1e-12*abs(fl['w_rho']).max()).all())}; sign of the heat flux v'T': {fl['v_T_sign']:+.0f} (poleward); phase difference between the lids: {np.degrees(fl['phase_tilt']):.1f} degrees")   # the three signatures
""",
    see="For the fastest-growing mode: the vertical density flux, the northward density flux, and the phase of the pressure "
        "wave, each against height.",
    read="The vertical flux is negative everywhere inside the layer (light fluid up, heavy fluid down — potential energy is "
         "released); the northward density flux is negative (light, warm air poleward); and the phase increases with height: "
         "the troughs lean westward by a quarter wavelength between the lids.",
    change="…the wave were shorter than the cut-off? No tilt, no fluxes, no growth: two edge waves that do not feel each other.")
nb.animation(r"""
nfr = 30 if FAST else 36                                       # number of frames
k_f = fast["alphaH"]/Lam_E                                     # wavenumber of the fastest mode [rad/m]
c_f = GFD.eady_phase_speed(fast["alphaH"], U0_a)               # its complex phase speed [m/s]
x_e = np.linspace(0.0, 2*2*np.pi/k_f, 120)                     # two wavelengths in x [m]
t_an = np.linspace(0.0, 3*t_e, nfr)                            # three e-folding times [s]
p_hat = md["p_hat"][:, None]                                   # complex vertical structure p-hat(z)
fig, ax = plt.subplots(figsize=(6.8, 3.4))
im = ax.imshow(np.zeros((len(z_e), len(x_e))), origin="lower", aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1,
               extent=[0, x_e[-1]/1e3, 0, H_a/1e3])            # the pressure perturbation in a longitude-height section
ax.set(xlabel="x east [km]", ylabel="z [km]")
fig.colorbar(im, label="p' (renormalised each e-folding)")

fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it: frames render much faster

def update(i):                                                 # frame i: the growing, travelling wave at time t_an[i]
    grow = np.exp((k_f*c_f.imag*t_an[i]) % 1.0)/np.e           # amplitude e^(sigma t), reset after every e-folding
    field = np.real(p_hat*np.exp(1j*k_f*(x_e[None, :] - c_f.real*t_an[i])))*grow   # Re[p-hat exp(ik(x - c_r t))] e^(sigma t)
    im.set_data(field)                                         # new picture
    ax.set_title(f"fastest Eady mode: t = {t_an[i]/86400:.2f} days = {t_an[i]/t_e:.2f} e-folding times (basic wind → 0 to {U0_a:.0f} m/s, bottom to top)", fontsize=8)
    return (im,)

show_animation(animate(update, frames=nfr, fig=fig, interval=100), player="video", dpi=50)   # a smooth video
plt.close(fig)                                                 # do not show the last frame a second time
""", explain=r"""
The pressure perturbation of the fastest-growing mode, $\mathrm{Re}\,[\hat p(z)e^{ik(x-ct)}]$, in a longitude–height section over
three e-folding times. The amplitude is reset after each e-folding (otherwise it would leave the colour scale); the title is
the clock.""")
nb.figure_notes(
    see="Highs (red) and lows (blue) that lean westward with height, travel east at half the upper-lid wind, and brighten "
        "(grow) until the amplitude is reset.",
    read="The westward tilt is the signature of a growing baroclinic wave: it puts poleward flow where the air is warm. The "
         "pattern moves at $U_0/2$ — the wind at mid-level, the \"steering level\".",
    change="…$\\alpha H=3$ (a wave shorter than the cut-off)? No tilt and no growth: two separate edge waves passing each other.")
explainer("eady_instability", "Where do storms come from?",
          "sliding the wavenumber through the cut-off turns a growing, westward-tilting mode into two separate neutral edge "
          "waves; changing the shear, the stratification or the latitude rescales the growth curve and prints the e-folding time "
          "in days and the preferred wavelength in kilometres.",
          ["Use the preset 'fastest-growing wave' and read the e-folding time.",
           "Drag the wavelength shorter until the status changes to 'two neutral edge waves'.",
           "Double N: what happens to the growth rate and to the preferred wavelength?",
           "Switch to the ocean mode: mesoscale eddies instead of weather systems."])
whatif(r"""
…the growing eddies become strong enough to interact with each other? Then we are in turbulence — but of an unusual, nearly
two-dimensional kind (`C17`).""")


# =====================================================================================================================
# A.18  §13.18 — C17 Fjørtoft's argument; pointers; summary
# =====================================================================================================================
nb.section("13.18", "Geostrophic Turbulence", intro=r"""
**What is this section about?** What turbulence does when rotation and stratification keep it nearly two-dimensional. Chapter
12's cascade runs backwards: energy moves to *larger* scales. ⚠️ In this section $\alpha$ is the enstrophy flux, $\eta$ (once)
the Kolmogorov length, $l$ an eddy size.""")
core("C17", r"Fjørtoft's argument: with $S_0=S_1+S_2$ and $K_0^2S_0=K_1^2S_1+K_2^2S_2$, $\dfrac{S_1}{S_2}=\dfrac{K_2-K_0}{K_0-K_1}\,\dfrac{K_2+K_0}{K_1+K_0}$ and $\dfrac{K_1^2S_1}{K_2^2S_2}=\dfrac{K_1^2}{K_2^2}\,\dfrac{K_2^2-K_0^2}{K_0^2-K_1^2}$ (13.145)",
     "Why do the eddies of the atmosphere and ocean merge into bigger eddies and jets instead of breaking down into smaller ones?")
problem(r"""
Stir a cup of coffee and the swirls break into smaller swirls until viscosity erases them (Ch. 12). Jupiter's atmosphere does
the opposite: small storms merge into large ones and into bands that have lasted centuries. So do ocean eddies. The difference
is one extra conservation law. In flat, two-dimensional flow a vortex cannot be stretched, so the total squared vorticity —
the enstrophy — is conserved along with the energy. Two conserved quantities are too many for a simple downhill cascade.""")
idea(words=r"""
A see-saw: energy and enstrophy sit on the same wavenumbers, but enstrophy weighs them by $K^2$. Move some energy to higher
$K$ and the enstrophy goes up; to keep it fixed you must move *more* energy to lower $K$.""")
fig(r"""
figA = draw_cascade_arrows()                                   # our sketch: energy to small wavenumbers, enstrophy to large ones
plt.show()
""",
    see="A sketch of the energy spectrum against wavenumber (both logarithmic) with the stirring scale $K_0$ in the middle: "
        "the energy arrow points to small wavenumbers (large eddies), the enstrophy arrow to large wavenumbers (small eddies). "
        "The two slopes written on the branches are explained in `D29` below.",
    read="The two conserved quantities leave the injection scale in opposite directions — the opposite of Chapter 12's "
         "three-dimensional cascade, where energy goes to small scales.",
    change="…the flow could stretch vortices (three dimensions)? Enstrophy would no longer be conserved and the energy would be "
           "free to go downscale.")
nb.recap("R31", "Where the energy finally goes", r"""
Energy must still reach the Kolmogorov scale to be dissipated; three-dimensional turbulence does that through the cascade with
spectrum $S_{11}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ (12.54, in its corrected form). Internal waves are a suggested route from the
large eddies to it.""", where="Ch. 12 §12.7 (C07–C08)")
nb.current_core = "C17"
P("P333", "enstrophy: mean-square vorticity, and its spectrum K²S(K)", r"""
Enstrophy is to vorticity what kinetic energy is to velocity: its mean square. A Fourier mode of wavenumber $K$ with velocity
amplitude $\hat u$ has vorticity amplitude $K\hat u$ (a derivative multiplies by $K$), so its share of the enstrophy is $K^2$
times its share of the energy.""", code=r"""
K = np.array([1.0, 2.0, 4.0]); S = np.array([4.0, 2.0, 1.0])   # energy in three modes
print("energy =", S.sum(), "; enstrophy =", (K**2*S).sum(), "; enstrophy per mode:", K**2*S)   # 7; 28, most of it in the SMALLEST scale
""")
note("N164 [B]", r"""
**Geostrophic turbulence** is nearly two-dimensional (rotation, stratification and thinness suppress $w$) and has no vortex
stretching. With isotropic spectra, $\overline{u^2}=\displaystyle\int_0^\infty S(K)\,dK$ and
$\overline{\zeta^2}=\displaystyle\int_0^\infty K^2S(K)\,dK$ (`GFD.enstrophy_spectrum(K, S)` is $K^2S$).""")
trap("T16", r"""
This $S$ is one-sided in $K$ with no factor ½; Chapter 12's spectra were two-sided, with $\int S$ equal to the variance of one
velocity component.""")
note("N165 [B]", r"""
**Two conservation statements** (inviscid, two-dimensional): $\dfrac d{dt}\displaystyle\int_0^\infty S(K)\,dK=0$ and
$\dfrac d{dt}\displaystyle\int_0^\infty K^2S(K)\,dK=0$ (13.143)–(13.144) (`SW.barotropic_invariants(zeta, L)` returns both integrals
for a vorticity field).""")
D("D28", ref="13.145")
nb.worked_example("where the energy of wavenumber 6 goes", r"""
Energy $S_0=1$ at $K_0=6$ is moved to $K_1=2$ and $K_2=9$ ($K_1=K_0/3$, $K_2=3K_0/2$).

1. Energy: $S_1+S_2=1$.
2. Enstrophy: $4S_1+81S_2=36$.
3. Subtract 4 × (1) from (2): $77S_2=32$ → $S_2=0.416$, $S_1=0.584$.
4. Energy ratio $S_1/S_2=1.41$: more energy went to the **larger** scale.
5. Enstrophy shares: $4\times0.584=2.34$ against $81\times0.416=33.66$, ratio 0.069: fourteen times more enstrophy went to the
   **smaller** scale.
6. Formula check: $\frac{K_2-K_0}{K_0-K_1}\frac{K_2+K_0}{K_1+K_0}=(3/4)(15/8)=1.406$ ✓.""")
code(r"""
r_f = GFD.fjortoft_transfer(1.0, 1/3, 1.5)                     # Eq. (13.145): energy at K0 = 1 shared between K1 = 1/3 and K2 = 3/2
print(f"S1 = {r_f['S1']:.4f}, S2 = {r_f['S2']:.4f}; energy ratio S1/S2 = {r_f['energy_ratio']:.5f}; enstrophy ratio = {r_f['enstrophy_ratio']:.5f} (= 1/{1/r_f['enstrophy_ratio']:.1f})")   # the triad of the worked example
for K2 in (1.2, 1.5, 2.0, 4.0):                                # move the small-scale partner further out
    rr = GFD.fjortoft_transfer(1.0, 1/3, K2)                   # same K1
    print(f"K2/K0 = {K2}: energy to large scale {100*rr['S1']:.0f} %, enstrophy to small scale {100/(1 + rr['enstrophy_ratio']):.0f} %")   # the two shares
print("enstrophy spectrum K^2 S of the primer's three modes:", GFD.enstrophy_spectrum(np.array([1.0, 2.0, 4.0]), np.array([4.0, 2.0, 1.0])))   # K^2 S
""", explain=r"""
1. `GFD.fjortoft_transfer(K0, K1, K2)` solves the two conservation equations for the shares $S_1$, $S_2$.
2. For the triad of the worked example, 58 % of the energy goes to the larger scale and 94 % of the enstrophy to the smaller.
3. The further out the small-scale partner, the larger the share of the energy that goes up-scale (from a third to more than
   nine tenths in the sweep); the enstrophy goes mostly down-scale in every case (90 % or more).""")
scratch(r"""
# From scratch: two equations, two unknowns
K0, K1, K2, S0 = 1.0, 1/3, 1.5, 1.0                            # the triad and the energy to be shared
A_f = np.array([[1.0, 1.0], [K1**2, K2**2]])                   # rows: energy conservation, enstrophy conservation
S1, S2 = np.linalg.solve(A_f, [S0, K0**2*S0])                  # the two shares
assert np.allclose([S1, S2], [r_f["S1"], r_f["S2"]])           # same as the library
print(f"✓ S1 = {S1:.4f}, S2 = {S2:.4f} from a 2 x 2 linear solve — the same as GFD.fjortoft_transfer")   # visible confirmation
""", r"""
Fjørtoft's argument is nothing more than two linear equations in two unknowns.""")
fig(r"""
k1 = np.linspace(0.05, 0.95, 200); k2 = np.linspace(1.05, 4.0, 200)   # K1/K0 below 1 and K2/K0 above 1
K1g, K2g = np.meshgrid(k1, k2)                                 # all triads
ratio = (K2g - 1)/(1 - K1g)*(K2g + 1)/(K1g + 1)                # Eq. (13.145), first member: S1/S2 with K0 = 1
fig, ax = plt.subplots(figsize=(6.0, 3.8))
pm = ax.pcolormesh(k1, k2, np.log10(ratio), cmap="PuOr", vmin=-2, vmax=2, shading="auto")   # log10 of the energy ratio
ax.contour(k1, k2, ratio, [1.0], colors=COLORS["ink"])         # the line S1 = S2
fig.colorbar(pm, label="log₁₀(S₁/S₂)")
ax.set(xlabel="K₁/K₀ (the larger scale)", ylabel="K₂/K₀ (the smaller scale)", title="More energy goes to the larger scale, except for lopsided triads")
plt.show()
print(f"fraction of this triad plane with S1 > S2: {100*(ratio > 1).mean():.0f} %")   # how common up-scale transfer is
""",
    see="The energy ratio $S_1/S_2$ for every triad: purple where more energy goes to the larger scale, orange where more goes "
        "to the smaller; the black line is $S_1=S_2$.",
    read="Over most of the plane the ratio exceeds 1. It falls below 1 only when the large-scale partner is very far from "
         "$K_0$ and the small-scale partner very near — and even then the *enstrophy* still goes mostly to the small scale.",
    change="…a triad symmetric in wavenumber, $K_0-K_1=K_2-K_0$? The ratio is $(K_2+K_0)/(K_1+K_0)>1$: always up-scale.")
D("D29")
note("N166 [B]", r"""
**Two inertial ranges** (ours, computed by the dimensional argument of `D29`; no source is cited): $S(K)\propto\varepsilon^{2/3}K^{-5/3}$
where an energy flux $\varepsilon$ passes toward small $K$, and $S(K)\propto\alpha^{2/3}K^{-3}$ where an enstrophy flux $\alpha$ passes
toward large $K$. The exponents follow from units alone, as in Chapter 1's Π theorem.""")
note("N167 [B]", r"""
**Where the up-scale transfer stops** (ours, by the scale argument of `D29`): at the length $l\sim\sqrt{u/\beta}$ where an eddy's
turnover rate matches the Rossby-wave frequency; beyond it the motion is wave-like and elongates into zonal jets
(`GFD.rhines_length`).""")
fig(r"""
K_ax = np.geomspace(1.0, 200.0, 200)                           # wavenumbers (arbitrary units)
S_ax = GFD.two_d_cascade_spectrum(K_ax, 15.0, 1.0, 15.0**2)    # the two power laws, joined at the injection wavenumber K0 = 15 (shape only)
fig, ax = plt.subplots(figsize=(6.0, 3.6))
ax.loglog(K_ax, S_ax, color=COLORS["accent"])                  # the spectrum
ax.axvline(15.0, color=COLORS["muted"], ls="--")               # the injection scale
ax.text(2.0, S_ax[0]*0.3, "energy range\nS ∝ ε^(2/3) K^(−5/3)\n← energy", fontsize=8)        # left range
ax.text(30.0, S_ax[0]*0.02, "enstrophy range\nS ∝ α^(2/3) K^(−3)\nenstrophy →", fontsize=8)  # right range
ax.set(xlabel="wavenumber K", ylabel="energy spectrum S(K)", title="Two ranges on either side of the injection scale (dimensional argument)")
plt.show()
u_atm, u_oc = inp["u_rms"]                                     # rms eddy speeds: 13 m/s (atmosphere) and 0.08 m/s (ocean), our inputs
print(f"stopping length sqrt(u/beta) at 35 N: atmosphere {GFD.rhines_length(u_atm, beta35)/1e3:.0f} km, ocean {GFD.rhines_length(u_oc, beta35)/1e3:.0f} km")   # where eddies turn into waves and jets
ex = np.linalg.solve(np.array([[2.0, -1.0], [-3.0, 0.0]]), [3.0, -2.0])   # exponents (a, b) of eps^a K^b: metres 2a - b = 3, seconds -3a = -2
print(f"dimensional check: eps^a K^b has the units of S (m^3 s^-2) for a = {ex[0]:.4f}, b = {ex[1]:.4f}")   # 2/3 and -5/3
""",
    see="The shape of the energy spectrum of two-dimensional turbulence forced at one scale (dashed line): a gentler slope on "
        "the large-scale side, a steeper one on the small-scale side.",
    read="Energy leaks toward the left, enstrophy toward the right. The exponents come from units alone (last printed line: "
         "$\\varepsilon$ in m² s⁻³ and $K$ in m⁻¹ combine to m³ s⁻² only as $\\varepsilon^{2/3}K^{-5/3}$). The stopping length is several hundred "
         "kilometres in the atmosphere and tens in the ocean.",
    change="…$\\beta$ were zero? Nothing would stop the growth of the eddies short of the size of the domain.")
P("P334", "the Jacobian J(ψ, ζ) and a pseudo-spectral step with 2/3 de-aliasing", r"""
The advection of vorticity by the flow it induces is $\mathbf u\cdot\nabla\zeta=\psi_x\zeta_y-\psi_y\zeta_x\equiv J(\psi,\zeta)$ (subscripts are
derivatives). A pseudo-spectral model takes derivatives in Fourier space (multiply by $ik$ — exact), multiplies fields on the
grid (cheap), and zeroes the top third of the wavenumbers before each product so that the product's new short waves do not
fold back onto long ones (aliasing). This model is OUR choice; the book describes such calculations only in words.""", code=r"""
n = 8; k = np.fft.fftfreq(n, 1/n)                  # integer wavenumbers 0..3, -4..-1
print("kept by the 2/3 rule:", np.abs(k) <= n/3)   # only |k| <= 2 survive
""")
note("N168 [B]", r"""
**Not in the book — ours; a qualitative demonstration.** Decaying two-dimensional turbulence,
$\dfrac{\partial\zeta}{\partial t}+J(\psi,\zeta)+\beta\dfrac{\partial\psi}{\partial x}=\nu\nabla^2\zeta$ with $\zeta=\nabla^2\psi$, without and with $\beta$.""")
choice(r"""
**our numerical model, cached**: two pseudo-spectral runs on a 64 × 64 grid (`SW.barotropic_run`, fourth-order Runge–Kutta,
2/3 rule), ten saved frames each, loaded with `ch13.load_reference_run`. At this resolution the steep small-scale range is
**not** resolved: everything said about these runs is qualitative, and no spectral slope is read off them.""")
nb.animation(r"""
runs_t = [ch13.load_reference_run("turbulence_f"), ch13.load_reference_run("turbulence_beta")]   # our two cached runs (beta = 0, beta > 0), or None
if any(r_ is None for r_ in runs_t):                           # a cache is missing
    print("a cached turbulence run is absent — regenerate it with scripts/ch13_make_caches.py; the sketch and the spectrum figure above stand on their own")
else:
    zet = [np.asarray(r_["zeta"], dtype=float) for r_ in runs_t]   # vorticity fields zeta(t, y, x)
    shown = [0, 1, 2, 4, 6, 9]                                 # six of the ten saved times (keeps the page light)
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.8))
    ims = [axes[j].imshow(zet[j][0], origin="lower", cmap="RdBu_r", vmin=-abs(zet[j][0]).max(), vmax=abs(zet[j][0]).max()) for j in (0, 1)]   # the two fields
    for ax, ttl in zip(axes[:2], ("no β", "with β")):
        ax.set_xticks([]); ax.set_yticks([])                   # no tick marks on the pictures
        ax.set_title(ttl, fontsize=9)
    for r_, ls, nm in ((runs_t[0], "-", "no β"), (runs_t[1], "--", "with β")):   # the two invariants against time
        axes[2].plot(r_["t"], r_["energy"]/r_["energy"][0], ls, color=COLORS["accent"], label=f"energy, {nm}")
        axes[2].plot(r_["t"], r_["enstrophy"]/r_["enstrophy"][0], ls, color=COLORS["amber"], label=f"enstrophy, {nm}")
    clock = axes[2].axvline(0, color=COLORS["muted"])          # a moving time marker
    axes[2].set(xlabel="time (model units)", ylabel="fraction of the initial value", title="Enstrophy decays, energy much less")
    axes[2].legend(fontsize=6)

    fig.canvas.draw(); fig.set_layout_engine("none")              # lay the figure out once, then freeze it: frames render much faster

    def update(n_):                                            # frame n_: the two vorticity fields at the saved time shown[n_]
        i = shown[n_]                                          # index of that saved time
        for j in (0, 1):
            ims[j].set_data(zet[j][i])                         # new field
            ims[j].set_clim(-abs(zet[j][i]).max(), abs(zet[j][i]).max())   # rescale the colours to this frame
        clock.set_xdata([runs_t[0]["t"][i]]*2)                 # move the time marker
        return ims

    show_animation(animate(update, frames=len(shown), fig=fig, interval=700), player="frames", dpi=36)   # step with the buttons
    plt.close(fig)                                             # do not show the last frame a second time
    for r_, z_, nm in zip(runs_t, zet, ("no beta", "with beta")):   # the numbers behind the pictures
        E_hat = abs(np.fft.fft2(z_[-1]))**2                    # squared vorticity amplitudes of the last frame …
        kk = 2*np.pi*np.fft.fftfreq(z_.shape[-1], d=float(r_["L"])/z_.shape[-1])   # … and their wavenumbers
        K2 = kk[None, :]**2 + kk[:, None]**2; K2[0, 0] = 1.0   # K^2 on the grid (the mean mode is skipped)
        E_k = E_hat/K2; E_k[0, 0] = 0.0                        # energy per mode = |zeta-hat|^2 / K^2
        print(f"{nm}: energy falls to {r_['energy'][-1]/r_['energy'][0]:.2f}, enstrophy to {r_['enstrophy'][-1]/r_['enstrophy'][0]:.2f} of the start; energy centroid K_E {r_['K_E'][0]:.1f} -> {r_['K_E'][-1]:.1f}; zonal (k_x = 0) share of the energy at the end: {E_k[:, 0].sum()/E_k.sum():.2f}")   # qualitative summary
""", explain=r"""
1. The two cached runs are loaded (if a file were missing the cell would say so and stop).
2. The left and middle panels show the vorticity field without and with β; the right panel the energy and the enstrophy as
   fractions of their initial values, with a marker at the time of the frame.
3. `player="frames"` lets you step through six of the ten saved times.
4. The printed lines summarise each run: how much energy and enstrophy are left, where the energy-weighted mean wavenumber
   $K_E$ has moved, and what share of the energy sits in the zonal ($k_x=0$) modes at the end.""")
nb.figure_notes(
    see="Many small vortices merging into a few large ones (left); with β (middle) the pattern stretches east–west into bands. "
        "Right: enstrophy drops to about a fifth of its initial value while energy keeps most of its own.",
    read="Qualitatively — and only qualitatively, at 64 × 64 — this is Fjørtoft's argument at work: enstrophy is drained at "
         "small scales by the viscosity while the energy survives and moves to larger scales (the centroid $K_E$ falls, printed "
         "lines). With β the zonal share of the energy is markedly larger: the beginning of jets.",
    change="…the grid were much finer? The small-scale range would be resolved and its slope could be compared with the "
           "dimensional argument of `D29`; we make no such claim from these runs.")
nb.pointer(r"""**S01** Exercises 13.1–13.8 are not reproduced. The one result the text relies on — the group velocity of
inertia–gravity waves — is derivation `D21` in `C14`. · **S02** Literature and supplemental reading: see the book's list; this
notebook cites only what was read first-hand.""")
whatif(r"""
…the fluid is a gas moving so fast that its own compressibility matters, or a wing rather than a planet? Those are Chapters
14 and 15. For climate dynamics, the next steps beyond this book are the equatorial waves named in `C15` and the wind-driven
gyre named in `C05`.""")

nb.summary(
    clicked=[
        r"`C01` — On a thin shell only the vertical part of the earth's spin, $f=2\Omega\sin\theta$ *(13.8)*, turns the wind.",
        "`C02` — Slow flow runs along the isobars because that is the only direction in which Coriolis can cancel the pressure force.",
        "`C03` — A horizontal temperature gradient sets how the wind changes with height.",
        "`C04` — Friction against rotation makes a layer of fixed thickness in which the current spirals.",
        "`C05` — The wind-driven transport is at right angles to the wind and does not depend on the eddy viscosity; its divergence pumps water up or down.",
        "`C06` — Near the ground friction lets the wind cross isobars toward low pressure: air rises in lows.",
        "`C07` — Three equations for surface height and two velocities hold all the dynamics of a thin layer.",
        "`C08` — A stratified ocean is a stack of such layers; the first internal one is about a metre deep.",
        "`C09` — One cubic: two fast gravity waves and one slow wave that needs β.",
        "`C10` — Rotation puts a floor $f$ under gravity-wave frequencies and turns currents into ellipses.",
        "`C11` — A coast replaces the missing cross-shore flow with a slope; the wave runs with the coast on its right ($f>0$).",
        "`C12` — Beyond a Rossby radius $c/f$ rotation holds a slope up against gravity.",
        r"`C13` — Each column keeps $(\zeta+f)/h$.",
        "`C14` — Internal waves live between $f$ and $N$; the slope of the crests sets the frequency.",
        "`C15` — Rossby crests drift west because displaced columns must change their spin; energy goes east for short waves.",
        r"`C16` — Storms grow by sliding fluid along the wedge between level and density surfaces, at a rate $0.30982\,(f/N)\,dU/dz$ (ours, computed).",
        "`C17` — Two conserved quantities send energy to large scales and enstrophy to small.",
    ],
    feeds_forward=[
        r'Ch. 15: the shallow-water ↔ gas-dynamics analogy uses $\sqrt{gH}$ as the "sound speed".',
        "Climate-dynamics work beyond the book: `fluidpy.core.gfd`, `vertical_modes` and `shallow_water` are importable on their own.",
    ],
    left_out=[
        "equatorial waves and the equatorial radius of deformation (any text on equatorial dynamics)",
        "the wind-driven gyre and western boundary currents (dynamical oceanography texts)",
        "the moist adiabat and pressure coordinates (dynamic meteorology texts)",
        "quasi-geostrophic theory in a continuously stratified fluid beyond the Eady problem, and the Charney problem",
        "observed data figures (replaced here by our own model output)",
    ],
)


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
    return sorted((k for k in curation_items("ch13") if k not in ids), key=lambda s: (s[0], int(s[1:])))


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
    _nb2.metadata["fluidpy"] = {"chapter": "ch13", "explainers": nb.explainers, "cores": sorted(nb.cores),
                                "recaps": sorted(nb.recaps), "primers": nb.primers, "derivations": nb.derivations}
    _out = ROOT / "notebooks" / "ch13_geophysical_fluid_dynamics.ipynb"
    for _i, _c in enumerate(_nb2.cells):
        _c["id"] = f"ch13-{_i:03d}"
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
