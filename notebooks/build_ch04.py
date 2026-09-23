"""Build the Chapter 4 teaching notebook: ``notebooks/ch04_conservation_laws.ipynb``.

Source of truth: ``analysis/ch04_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract
— every call goes to ``fluidpy.ch04_conservation_laws``, imported as ``ch04``, which re-exports the nine new core
modules), Part E (prerequisite ledger → 23 primers P111–P133 and one-line reminders of earlier primers), Part F (the 30
derivations D01–D30, one move per step) and ``analysis/ch04_curation.md`` (IDs, depths, section coverage, §6 animations
and plotly figures, §7 from-scratch moments). Physics lives in ``fluidpy``; cells only call it.

**Derivations are read from Part F at build time** (``part_f()`` below, the ch03 parser): goal, start, plan, tools,
assumptions, every step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word, so the
notebook and the design cannot drift apart. Explainer-only fields (*live*, *set*, *watch*) are dropped. Where a Part F
Start/Result names an equation by number only, the equation itself is written (``EQ``). The four sympy checks (D06,
D09 ★★★, D15 ★★★, D30) are written here, every line commented.

Conventions (design header, binding): stress first index = face normal, (∇·τ)_j = ∂τ_ij/∂x_i; z up, g = 9.81 m/s²;
χ = −z ⇒ ρu = ∂ψ/∂y, ρv = −∂ψ/∂x; +2Ω×u′ is the Coriolis *term* of (4.43), −2Ω×u′ the Coriolis *force* per unit mass
of (4.45); four meanings of the prime; the shear rate is written γ̇, the extension rate s; Kundu's lapse rate
Γ = dT/dz with the meteorological −dT/dz alongside; the book's slips are taught in corrected form.

Run:  .venv/Scripts/python.exe notebooks/build_ch04.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch04.py --dump     (prints the parsed Part F derivations only)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch04")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# (Transcribed from the rendered pages p124–p173 via the design; earlier chapters' equations from their notebooks.)
# ---------------------------------------------------------------------------------------------------------------------
EQ = {
    "1.2": r"\mathbf q=-k\nabla T",
    "1.3": r"\tau=\mu\,du/dy",
    "1.5": r"\Delta p=\sigma\Big(\frac1{R_1}+\frac1{R_2}\Big)",
    "1.8": r"dp/dz=-\rho g",
    "1.10": r"de=\delta q+\delta w",
    "1.18": r"T\,ds=de+p\,dv",
    "1.20": r"\alpha=-\frac1\rho\Big(\frac{\partial\rho}{\partial T}\Big)_p",
    "2.15": r"f_j=n_i\tau_{ij}",
    "2.19": r"\varepsilon_{ijk}\varepsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}",
    "2.30": r"\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA",
    "3.5": r"\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F",
    "3.11": r"\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}",
    "3.12": r"S_{ij}=\tfrac12\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)",
    "3.13": r"R_{ij}=\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}",
    "3.14": r"\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}",
    "3.17": r"\boldsymbol\omega=\nabla\times\mathbf u=0",
    "3.35": r"\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA",
    "4.1": r"\frac{d}{dt}\int_{V(t)}\rho(\mathbf x,t)\,dV=0",
    "4.2": r"\int_{V(t)}\frac{\partial\rho}{\partial t}dV+\int_{A(t)}\rho\,\mathbf u\cdot\mathbf n\,dA=0",
    "4.3": r"\frac{d}{dt}\int_{V^*(t)}\rho\,dV-\int_{V^*(t)}\frac{\partial\rho}{\partial t}dV-\int_{A^*(t)}\rho\,\mathbf b\cdot\mathbf n\,dA=0",
    "4.4": r"\int_{V^*}\frac{\partial\rho}{\partial t}dV=\int_{V}\frac{\partial\rho}{\partial t}dV=-\int_{A}\rho\mathbf u\cdot\mathbf n\,dA=-\int_{A^*}\rho\mathbf u\cdot\mathbf n\,dA",
    "4.5": r"\frac{d}{dt}\int_{V^*(t)}\rho\,dV+\int_{A^*(t)}\rho\,(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0",
    "4.6": r"\int_{V}\Big\{\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)\Big\}dV=0",
    "4.7": r"\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0",
    "4.8": r"\frac{1}{\rho}\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0",
    "4.9": r"\frac{D\rho}{Dt}\equiv\frac{\partial\rho}{\partial t}+\mathbf u\cdot\nabla\rho=0",
    "4.10": r"\nabla\cdot\mathbf u=0",
    "4.11": r"\nabla\cdot(\rho\mathbf u)=0",
    "4.12": r"\rho\mathbf u=\nabla\times\boldsymbol\Psi=\nabla\chi\times\nabla\psi",
    "4.13": r"\frac{d}{dt}\int_{V(t)}\rho\mathbf u\,dV=\int_{V(t)}\rho\mathbf g\,dV+\int_{A(t)}\mathbf f\,dA",
    "4.14": r"\int_{V(t)}\frac{\partial}{\partial t}(\rho\mathbf u)dV+\int_{A(t)}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{V(t)}\rho\mathbf g\,dV+\int_{A(t)}\mathbf f\,dA",
    "4.15": r"\int_{V^*}\frac{\partial}{\partial t}(\rho\mathbf u)dV=\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV-\int_{A^*}\rho\mathbf u\,(\mathbf b\cdot\mathbf n)\,dA",
    "4.16a": r"\int_{V}\frac{\partial(\rho\mathbf u)}{\partial t}dV=\int_{V^*}\frac{\partial(\rho\mathbf u)}{\partial t}dV",
    "4.16b": r"\int_{A}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{A^*}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA",
    "4.16c": r"\int_{V}\rho\mathbf g\,dV=\int_{V^*}\rho\mathbf g\,dV",
    "4.16d": r"\int_{A}\mathbf f\,dA=\int_{A^*}\mathbf f\,dA",
    "4.17": r"\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA",
    "4.18": r"\mathbf g=-\nabla\Phi",
    "4.19": r"\tfrac12U^2+gz+p/\rho=\text{const along a streamline}",
    "4.20a": r"\int_{A}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{V}\frac{\partial}{\partial x_i}(\rho u_iu_j)dV",
    "4.20b": r"\int_{A}n_i\tau_{ij}dA=\int_{V}\frac{\partial\tau_{ij}}{\partial x_i}dV",
    "4.21": r"\int_{V}\Big\{\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)-\rho g_j-\frac{\partial\tau_{ij}}{\partial x_i}\Big\}dV=0",
    "4.22": r"\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}",
    "4.23": r"\frac{\partial(\rho u_j)}{\partial t}+\frac{\partial(\rho u_iu_j)}{\partial x_i}=\rho\frac{Du_j}{Dt}",
    "4.24": r"\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}",
    "4.25": r"\tau_{ij}=\tau_{ji}",
    "4.26": r"\tau_{ij}=-p\,\delta_{ij}",
    "4.27": r"\tau_{ij}=-p\,\delta_{ij}+\sigma_{ij}",
    "4.28": r"\sigma_{ij}=K_{ijmn}S_{mn}",
    "4.29": r"K_{ijmn}=\lambda\delta_{ij}\delta_{mn}+\mu\delta_{im}\delta_{jn}+\gamma\delta_{in}\delta_{jm}",
    "4.30": r"\gamma=\mu",
    "4.31": r"\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}",
    "4.32": r"p=-\tfrac13\tau_{ii}+\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u",
    "4.33": r"\bar p\equiv-\tfrac13\tau_{ii}",
    "4.34": r"p-\bar p=\big(\tfrac23\mu+\lambda\big)\nabla\cdot\mathbf u",
    "4.35": r"\tau_{ij}=-p\,\delta_{ij}+2\mu S_{ij}",
    "4.36": r"\lambda+\tfrac23\mu=0",
    "4.37": r"\tau_{ij}=-p\,\delta_{ij}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)+\mu_vS_{mm}\delta_{ij}",
    "4.38": r"\rho\frac{Du_j}{Dt}=-\frac{\partial p}{\partial x_j}+\rho g_j+\frac{\partial}{\partial x_i}\Big[\mu\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\big(\mu_v-\tfrac23\mu\big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big]",
    "4.39a": r"\rho\frac{Du_j}{Dt}=-\frac{\partial p}{\partial x_j}+\rho g_j+\mu\frac{\partial^2u_j}{\partial x_i^2}+\big(\mu_v+\tfrac13\mu\big)\frac{\partial}{\partial x_j}\frac{\partial u_m}{\partial x_m}",
    "4.39b": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u",
    "4.39": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u",
    "4.40": r"\mu\nabla^2\mathbf u=2\mu\frac{\partial S_{ij}}{\partial x_i}=-\mu\nabla\times\boldsymbol\omega",
    "4.41": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g",
    "4.42": r"\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x'",
    "4.43": r"\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')",
    "4.44": r"\frac{D\mathbf u}{Dt}=\frac{D'\mathbf u'}{Dt}+\frac{d\mathbf U}{dt}+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')",
    "4.45": r"\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho\Big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\Big]+\mu\nabla'^2\mathbf u'",
    "4.46": r"\frac{d}{dt}\int_{V(t)}\rho\big(e+\tfrac12\lvert\mathbf u\rvert^2\big)dV=\int_{V(t)}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A(t)}\mathbf f\cdot\mathbf u\,dA-\int_{A(t)}\mathbf q\cdot\mathbf n\,dA",
    "4.47": r"\int_{V}\frac{\partial}{\partial t}\big(\rho e+\tfrac\rho2\lvert\mathbf u\rvert^2\big)dV+\int_{A}\big(\rho e+\tfrac\rho2\lvert\mathbf u\rvert^2\big)(\mathbf u\cdot\mathbf n)dA=\int_{V}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A}\mathbf f\cdot\mathbf u\,dA-\int_{A}\mathbf q\cdot\mathbf n\,dA",
    "4.48": r"\frac{d}{dt}\int_{V^*}\rho\big(e+\tfrac12\lvert\mathbf u\rvert^2\big)dV+\int_{A^*}\big(\rho e+\tfrac\rho2\lvert\mathbf u\rvert^2\big)(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A^*}\mathbf f\cdot\mathbf u\,dA-\int_{A^*}\mathbf q\cdot\mathbf n\,dA",
    "4.49": r"\int_{A}\rho\big[e+\tfrac12u_j^2\big]u_in_i\,dA=\int_{V}\frac{\partial}{\partial x_i}\Big(\rho\big[e+\tfrac12u_j^2\big]u_i\Big)dV",
    "4.50": r"\int_{A}n_i\tau_{ij}u_j\,dA=\int_{V}\frac{\partial}{\partial x_i}\big(\tau_{ij}u_j\big)dV",
    "4.51": r"\int_A q_in_i\,dA=\int_V\frac{\partial q_i}{\partial x_i}dV",
    "4.52": r"\int_{V}\Big\{\frac{\partial(\rho E)}{\partial t}+\frac{\partial(\rho Eu_i)}{\partial x_i}-\rho g_iu_i-\frac{\partial(\tau_{ij}u_j)}{\partial x_i}+\frac{\partial q_i}{\partial x_i}\Big\}dV=0",
    "4.53": r"\frac{\partial}{\partial t}\big(\rho\big[e+\tfrac12u_j^2\big]\big)+\frac{\partial}{\partial x_i}\big(\rho\big[e+\tfrac12u_j^2\big]u_i\big)=\rho g_iu_i+\frac{\partial}{\partial x_i}(\tau_{ij}u_j)-\frac{\partial q_i}{\partial x_i}",
    "4.54": r"\frac{\partial(\tau_{ij}u_j)}{\partial x_i}=\Big(-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}\Big)+\Big(-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}\Big)",
    "4.55": r"\rho\frac{D}{Dt}\big(e+\tfrac12u_j^2\big)=\rho g_iu_i+\Big(-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}\Big)+\Big(-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}\Big)-\frac{\partial q_i}{\partial x_i}",
    "4.56": r"\rho\frac{D}{Dt}\big(\tfrac12u_j^2\big)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}",
    "4.57": r"\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}",
    "4.58": r"\varepsilon\equiv\frac1\rho\sigma_{ij}S_{ij}=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}{\rho}S_{mm}^2\ge0",
    "4.59": r"\sigma_{ij}=\mu\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)+\big(\mu_v-\tfrac23\mu\big)\frac{\partial u_m}{\partial x_m}\delta_{ij}",
    "4.60": r"\rho\frac{De}{Dt}=-p\frac{\partial u_m}{\partial x_m}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\mu_vS_{mm}^2+\frac{\partial}{\partial x_i}\Big(k\frac{\partial T}{\partial x_i}\Big)",
    "4.61": r"\frac{De}{Dt}=T\frac{Ds}{Dt}-p\frac{D(1/\rho)}{Dt}",
    "4.62": r"\frac{Ds}{Dt}=-\frac1\rho\frac{\partial}{\partial x_i}\Big(\frac{q_i}{T}\Big)-\frac{q_i}{\rho T^2}\frac{\partial T}{\partial x_i}+\frac{\varepsilon}{T}",
    "4.63": r"\frac{k}{\rho T^2}\lvert\nabla T\rvert^2+\frac\varepsilon T\ge0",
    "4.64": r"d\mathbf H/dt=\mathbf M",
    "4.65": r"\frac{d}{dt}\int_{V_o}(\mathbf r\times\rho\mathbf u)dV+\int_{A_o}(\mathbf r\times\rho\mathbf u)(\mathbf u\cdot\mathbf n)dA=\int_{V_o}(\mathbf r\times\rho\mathbf g)dV+\int_{A_o}(\mathbf r\times\mathbf f)dA",
    "4.66": r"\frac{\partial u_j}{\partial t}+u_i\frac{\partial u_j}{\partial x_i}=-\frac1\rho\frac{\partial p}{\partial x_j}-\frac{\partial\Phi}{\partial x_j}",
    "4.67": r"\frac1\rho\frac{\partial p}{\partial x_j}=\frac{\partial}{\partial x_j}\int_{p_o}^{p}\frac{dp'}{\rho(p')}",
    "4.68": r"u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\big(\tfrac12u_i^2\big)",
    "4.69": r"\frac{\partial u_j}{\partial t}+\frac{\partial}{\partial x_j}\Big[\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\Big]=(\mathbf u\times\boldsymbol\omega)_j",
    "4.70": r"\nabla B=\mathbf u\times\boldsymbol\omega",
    "4.71": r"\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const along streamlines and vortex lines}",
    "4.72": r"\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const everywhere}",
    "4.73": r"\mathbf u\equiv\nabla\phi",
    "4.74": r"\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=B(t)",
    "4.75": r"\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const}",
    "4.76": r"\rho u_i\frac{\partial}{\partial x_i}\big(e+\tfrac12u_j^2\big)=\rho u_ig_i-\frac{\partial}{\partial x_j}\big(\rho u_j\tfrac p\rho\big)",
    "4.77": r"\rho u_i\frac{\partial}{\partial x_i}\big(e+\tfrac p\rho+\tfrac12u_j^2+gz\big)=0",
    "4.78": r"h+\tfrac12\lvert\mathbf u\rvert^2+gz=\text{const on streamlines}",
    "4.79": r"\rho\frac{\partial\mathbf u}{\partial t}+\rho\nabla\big(\tfrac12\lvert\mathbf u\rvert^2\big)-\rho\mathbf u\times\boldsymbol\omega=-\nabla p+\rho\mathbf g-\mu\nabla\times\boldsymbol\omega",
    "4.80": r"\rho\frac{\partial\mathbf u}{\partial t}+\nabla\big(\tfrac12\rho\lvert\mathbf u\rvert^2+\rho gz+p\big)=0",
    "4.81": r"\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s+\int_1^2\frac{\partial}{\partial s}\big(\tfrac12\lvert\mathbf u\rvert^2+gz+\tfrac p\rho\big)ds=0",
    "4.82": r"\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s+\big(\tfrac12\lvert\mathbf u\rvert^2+gz+\tfrac p\rho\big)_2=\big(\tfrac12\lvert\mathbf u\rvert^2+gz+\tfrac p\rho\big)_1",
    "4.83": r"\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+\frac p\rho=\text{const}",
    "4.84": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u",
    "4.85": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\mu\nabla^2\mathbf u",
    "4.86": r"\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u",
    "4.87": r"\rho\frac{De}{Dt}=-p\nabla\cdot\mathbf u+\rho\varepsilon-\nabla\cdot\mathbf q",
    "4.88": r"\rho C_p\frac{DT}{Dt}=\rho\varepsilon-\nabla\cdot\mathbf q",
    "4.89": r"\frac{DT}{Dt}=\kappa\nabla^2T",
    "4.90": r"\frac{\partial\eta}{\partial t}+(\mathbf u_s\cdot\nabla)\eta=0\ \text{on}\ \eta=0",
    "4.91": r"\frac{\partial\eta}{\partial t}+(\mathbf u\cdot\nabla)\eta\equiv\frac{D\eta}{Dt}=0\ \text{on}\ \eta=0",
    "4.92": r"(u_{rel})_n=\mathbf u\cdot\mathbf n-\mathbf u_s\cdot\mathbf n=\frac{1}{\lvert\nabla\eta\rvert}\frac{D\eta}{Dt}",
    "4.93": r"\frac{\rho}{\lvert\nabla\eta\rvert}\frac{D\eta}{Dt}\ \text{on}\ \eta=0",
    "4.94": r"f=e-Ts",
    "4.95": r"df=de-T\,ds-s\,dT",
    "4.96": r"F=\rho_1V_1f_1+\rho_2V_2f_2+A\sigma",
    "4.97": r"(F_p)_z=-\pi\Delta p\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}",
    "4.98": r"(F_{st})_z=\pi\sigma\sqrt{2R_1\zeta}\sqrt{2R_2\zeta}\Big(\frac1{R_1}+\frac1{R_2}\Big)",
    "4.99": r"\frac{F_D}{\rho U^2d^2}=\Psi\Big(\frac{\rho Ud}{\mu}\Big)",
    "4.100": r"x_i^*=\frac{x_i}{l},\ t^*=\Omega t,\ u_j^*=\frac{u_j}{U},\ p^*=\frac{p-p_\infty}{\rho U^2},\ g_j^*=\frac{g_j}{g}",
    "4.101": r"\Big[\frac{\Omega l}{U}\Big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\Big[\frac{gl}{U^2}\Big]\mathbf g^*+\Big[\frac{\mu}{\rho Ul}\Big]\nabla^{*2}\mathbf u^*",
    "4.102": r"\mathrm{St}=\frac{\Omega l}{U}",
    "4.103": r"\mathrm{Re}=\frac{\rho Ul}{\mu}",
    "4.104": r"\mathrm{Fr}=\frac{U}{\sqrt{gl}}",
    "4.105": r"\mathrm{Fr}'=\frac{U}{\sqrt{g'l}},\ \mathrm{Ri}=\frac1{\mathrm{Fr}'^2}",
    "4.106": r"C_p\equiv\frac{p-p_\infty}{\tfrac12\rho U^2}",
    "4.107": r"C_D\equiv\frac{F_D}{\tfrac12\rho U^2A}",
    "4.108": r"C_L\equiv\frac{F_L}{\tfrac12\rho U^2A}",
    "4.109": r"t^*=Ut/l,\ p^*=(p-p_\infty)/\rho_oU^2,\ \rho^*=\rho/\rho_o",
    "4.110": r"\nabla^*\cdot\mathbf u^*=-\Big[\frac{U^2}{c^2}\Big]\frac1{\rho^*}\frac{Dp^*}{Dt^*}",
    "4.111": r"M\equiv U/c",
    "4.112": r"\rho\frac{Dh}{Dt}=\frac{Dp}{Dt}+\rho\varepsilon+\frac{\partial}{\partial x_i}\Big(k\frac{\partial T}{\partial x_i}\Big)",
    "4.113": r"\varepsilon^*=\frac{\rho_ol^2}{\mu_oU^2}\varepsilon,\ T^*=\frac{T-T_o}{T_w-T_o}",
    "4.114": r"\rho^*\frac{DT^*}{Dt^*}=\mathrm{Ec}\frac{Dp^*}{Dt^*}+\frac{\mathrm{Ec}}{\mathrm{Re}}\varepsilon^*+\frac{1}{\mathrm{Pr}\,\mathrm{Re}}\nabla^*\cdot(k^*\nabla^*T^*)",
    "4.115": r"\mathrm{Ec}\equiv\frac{U^2}{C_p(T_w-T_o)}",
    "4.116": r"\mathrm{Pr}\equiv\frac{\nu}{\kappa}",
    "4.117": r"\mathrm{We}\equiv\frac{\rho U^2l}{\sigma}",
    "4.118": r"\mathrm{Bo}\equiv\frac{\rho l^2g}{\sigma}",
    "4.119": r"\mathrm{Ca}\equiv\frac{\mu U}{\sigma}",
    "B3.6": r"\nabla\cdot(\rho\mathbf u)=\mathbf u\cdot\nabla\rho+\rho\nabla\cdot\mathbf u",
}
EQ["4.16a–d"] = EQ["4.16a"] + r",\ \ldots"
EQ["4.16"] = EQ["4.16a"] + r",\ \ldots"


def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ *({num})*"


_EQ_REF = re.compile(r"\(((?:\d|B\d)\.\d+[a-d]?)\)")


def show_eqs(text: str) -> str:
    """Write the equation next to the first mention of a book equation number in a text that names it without
    writing it — the house rule "show the equation, not just its number"."""
    done: set[str] = set()

    def rep(m):
        n = m.group(1)
        nxt = text[m.end():m.end() + 16]
        if (n in done or n not in EQ or EQ[n] in text or nxt.lstrip(", :").startswith("$")
                or nxt.startswith(" — the line above") or nxt.startswith(" — the result")
                or nxt[:1] in "/–-)" or text[max(0, m.start() - 1):m.start()] == "("):
            return m.group(0)
        done.add(n)
        return f"({n}), ${EQ[n]}$"
    return _EQ_REF.sub(rep, text)


# ---------------------------------------------------------------------------------------------------------------------
# Part F reader: the derivations, word for word (the ch03 parser)
# ---------------------------------------------------------------------------------------------------------------------
def _join(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(x.strip() for x in lines)).strip()


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
            tokens += [("w", w.replace("**", "").replace("`", "")) for w in p.split()]
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


def _numbers_as_maths(s: str) -> str:
    """A Start/Result that names equations by number only ("(4.55) and (4.56)") → the equations themselves."""
    if "$" in s:
        return s
    nums = re.findall(r"\(((?:\d)\.\d+[a-d]?)\)", s)
    if not nums:
        return s
    return " ".join(f"$({n})\\quad {EQ[n]}$" for n in nums if n in EQ) or s


_SELF = re.compile(r"(This is (?:the book's )?|This chain is the book's |is the book's |it is the book's |This is exactly )\(((?:\d)\.\d+[a-d]?)\)")


def _self_ref(text: str) -> str:
    """"This is (4.49)." inside a step's *why* points at the line just displayed: say so instead of repeating it."""
    return _SELF.sub(lambda m: f"{m.group(1)}({m.group(2)}) — the line above", text)


def part_f() -> dict[str, dict]:
    """Parse Part F of the design into {D01: dict(title, goal, start, plan, tools, assumptions, steps, result, check,
    meaning, traps)}."""
    text = (ROOT / "analysis" / "ch04_design.md").read_text(encoding="utf-8")
    part = text.split("## Part F", 1)[1]
    chunks = re.split(r"^### (D\d\d) · ", part, flags=re.M)
    out: dict[str, dict] = {}
    for key, body in zip(chunks[1::2], chunks[2::2]):
        lines = body.splitlines()
        title = re.split(r" — ★", lines[0])[0].strip()
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
        f = {k: _join(v) for k, v in fields.items()}
        parsed = []
        for st in steps:
            s = _join(st)
            bits = re.split(r"(?:^|\s·\s)\*(did|tex|why|plain|live|set|watch):\*\s*", s)
            d = {name: val.strip() for name, val in zip(bits[1::2], bits[2::2])}
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=_self_ref(d["why"]),
                               plain=_self_ref(d["plain"])))
        start_tex, start_plain = _split_words(f["Start"])
        res_tex, res_plain = _split_words(f["Result"])
        res_tex = re.sub(r"\s*\((?:\d)\.\d+[a-d]?\)\s*$", "", res_tex) if res_tex.count("$") >= 2 else res_tex
        plan = [p.strip() for p in re.split(r"\(\d+\)\s*", f["Plan"]) if p.strip()]
        tools = [t.strip() for t in f["Tools"].split(" · ") if t.strip()]
        out[key] = dict(title=title, goal=f["Goal"], start=(display_tex(_numbers_as_maths(start_tex)), start_plain),
                        plan=plan, tools=tools, assumptions=f.get("Assumptions", ""), steps=parsed,
                        result=(display_tex(_numbers_as_maths(res_tex)), res_plain), check=f.get("Check", ""),
                        meaning=f.get("What it means", ""), traps=f.get("Traps", ""))
    return out


PF = part_f()


def pf_sub(key: str, field: str, old: str, new: str) -> None:
    """Edit one Part F field after parsing (lesson-review fixes); fails loudly if the design text changed.
    ``field`` = goal/check/meaning/traps/assumptions, ``tools``, or ``stepN.did|why|plain``."""
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


# lesson review round 1: every check names a computation that the notebook really runs
pf_sub("D13", "check", "the three forms agree to stencil accuracy (`viscous_force_forms`)",
       "the three forms agree to stencil accuracy (`viscous_force_forms`; the viscous-paradox figure cell of this block asserts it)")
pf_sub("D21", "check", "`ch04.kinetic_energy_budget` residual ≈ 0 on a test field (C10).",
       "`ch04.kinetic_energy_budget` residual ≈ 0 on the Taylor–Green test field (the cell right after this derivation) ✓.")
pf_sub("D24", "check", "Rankine core (r = 0.5 m, Γ = 2π, σ = 1)",
       "Rankine core (r = 0.5 m, Γ = 2π, σ = 1; the vortex is worked with numbers further down this block)")
pf_sub("D24", "check", "`ch04.lamb_identity_sym` returns zero vectors ✓.",
       "`ch04.lamb_identity_sym` returns zero vectors for a general 3-D field (the cell right after this derivation) ✓.")
pf_sub("D28", "check", "sympy (C13 note N110): e = C_vT",
       "sympy (the cell after note N110 below): e = C_vT")
pf_sub("D12", "check", "The toy field u = (ax², 0, 0) has ∇·u = 2ax and a nonzero (μ_v + ⅓μ)·2a term ✓.",
       "The toy field u = (ax², 0, 0) has ∇·u = 2ax and a nonzero (μ_v + ⅓μ)·2a term (both computed in the sympy cell right after this derivation) ✓.")
pf_sub("D11", "check", "Poiseuille: 0 = 100 + 0 − 100 ✓.",
       "Poiseuille (the plane channel flow worked with numbers later in this block, G = 100 Pa/m): 0 = 100 + 0 − 100 ✓.")
pf_sub("D05", "check", "Bore with b = U: storage 0 ✓ (E1).",
       "Bore with b = U: storage 0 ✓ (the `control_volume_budgets` explainer's bore mode).")
pf_sub("D10", "check", "(E3 parity row)", "(the `newtonian_stress_lab` explainer checks the same numbers)")
pf_sub("D28", "traps", "(high Eckert number: E7's oil bearing)",
       "(high Eckert number: the oil-bearing preset of the `viscous_dissipation_heating` explainer)")
pf_sub("D20", "check", "on a random field, C10", "on the Taylor–Green test field, the code cell after note N74")
pf_sub("D06", "check", "and prints `U*dU + g*dz + dp/rho` = 0.",
       "and prints `U*U_s + g*sin(theta) + p_s/rho` (U_s = ∂U/∂s, p_s = ∂p/∂s; with sin θ = dz/ds this is d(½U² + gz + p/ρ)/ds) together with a residual 0 against (4.19).")
pf_sub("D04", "tools", "curl of a product (gloss with the sympy check in N14's cell)",
       "curl of a product (gloss in step 2's why)")
pf_sub("D04", "step4.why", "the minus sign is the usual convention.",
       "we pick the minus sign because it makes ρu = +∂ψ/∂y (step 6), the book's convention.")
pf_sub("D14", "meaning", "why vorticity seen from the Earth is ζ, not ζ + f",
       "why the vorticity seen from the Earth, ζ (the relative vorticity, Ch. 13), differs from the absolute vorticity ζ + f, "
       "f = 2Ω sin φ being the Coriolis parameter primed later in this block (P126)")
pf_sub("D18", "step5.why", "Sum of two conservative forces is conservative;",
       "Here g_n = −∇Φ_n is the Newtonian gravitation alone (the attraction of the Earth's mass, potential Φ_n; N64 below gives numbers). The sum of two conservative forces is conservative;")
pf_sub("D26", "step3.did", "Insert into the steady-less (4.69)", "Insert into (4.69) with ω = 0")
pf_sub("D26", "step3.why", "ω = 0 removes u × ω;", "In (4.69) the right side u × ω vanishes because ω = 0;")

if "--dump" in sys.argv:
    for k, d in PF.items():
        print(f"=== {k} {d['title']}  ({len(d['steps'])} steps)")
        print("  START", d["start"][0][:150], "|", d["start"][1][:80])
        for i, s in enumerate(d["steps"], 1):
            print(f"  {i}. {s['did']} :: {s['tex'][:140]}")
        print("  RESULT", d["result"][0][:150], "|", d["result"][1][:80])
        print("  CHECK", d["check"][:300])
    sys.exit(0)

# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch03.py)
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
    """One 🔁 cell reminding several tools primed in earlier chapters (knowledge/primers.md) — one sentence each.
    ``items`` = [(Part E concept text, reminder sentence)]; the concepts are registered so the ledger check finds them."""
    body = "\n".join(f"- **{c}** — {textwrap.dedent(t).strip()}" for c, t in items)
    head = "> 🔁 **Tools from earlier chapters used here** (see `knowledge/primers.md`)"
    nb.md(f"{head}{(' — ' + lead) if lead else ''}\n\n{body}", tags=["primer"])
    for c, _ in items:
        nb.primers.append(f"{c} (reminder)")


def note(nid: str, title: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note with its curation id."""
    body = textwrap.dedent(text).strip()
    sep = " " if body.startswith("—") or not body else " — "
    nb.note(f"**{title}** `{nid}`{sep}{body}", equation=equation, ref=ref)


def _tidy_check(s: str) -> str:
    s = re.sub(r"\*\*sympy check \(`check_src`, every line commented\):\*\*", "**sympy check (the cell below, every line commented):**", s)
    s = re.sub(r"`check_src` \(optional, ★★\):", "sympy check (the cell below):", s)
    s = re.sub(r"sympy \(`check_src` optional\):", "sympy check (the cell below):", s)
    return s


def D(key: str, ref: str = "", check_src: str | None = None, extra_check: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then its traps as a ⚠️ callout."""
    d = PF[key]
    S = show_eqs
    check = S(_tidy_check(d["check"])) + (f" {extra_check}" if extra_check else "")
    goal = S(d["goal"]) + (f"\n\n**Assumptions.** {S(d['assumptions'])}" if d["assumptions"] else "")
    steps = [dict(st, why=S(st["why"]), plain=S(st["plain"])) for st in d["steps"]]
    nb.derivation(key, f"{d['title']} `{key}`", ref=ref, goal=goal, start=(d["start"][0], S(d["start"][1])),
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


def gloss(concepts: list[str]) -> None:
    """Register Part E concepts that the primer just written explains in passing (its text or code comment)."""
    for c in concepts:
        nb.primers.append(f"{c} (glossed in the primer above)")


_MATH_SPLIT = re.compile(r"(```.*?```|\$\$.*?\$\$|\$[^$\n]+\$)", re.S)
_NUM = re.compile(r"\*?\(((?:\d|B\d)\.\d+[a-d]?)\)\*?")


def _show_all_eqs(cell_src: str) -> str:
    """Final pass over a markdown cell: a book equation named by number in plain text gets the equation written next to
    it (first mention only), unless the equation is already shown in the cell (the same LaTeX, or a displayed line
    tagged with that number) or sits right next to the mention."""
    parts = _MATH_SPLIT.split(cell_src)
    done: set[str] = set()
    out = []
    for k, seg in enumerate(parts):
        if k % 2 == 1:                                   # math or code: leave alone
            out.append(seg)
            continue

        def rep_(m, seg=seg, k=k):
            n = m.group(1)
            if n not in EQ or n in done:
                return m.group(0)
            if EQ[n] in cell_src or ("\\text{(" + n + ")}") in cell_src:
                return m.group(0)                        # the equation is displayed in this cell already
            before, after = seg[:m.start()], seg[m.end():]
            line = before[before.rfind("\n") + 1:]
            if line.startswith("**Step ") or line.startswith("#"):
                return m.group(0)                        # headings stay short; the equation goes to a later mention
            if after.startswith(" — the line above") or after.startswith(" — the result"):
                return m.group(0)                        # a derivation step pointing at the line it just displayed
            if k > 0 and re.fullmatch(r"\s*(const(ant)?)?[\s*(]*(Eq\.\s*)?", before):
                return m.group(0)                        # the equation is written just before the number
            if k + 1 < len(parts) and re.fullmatch(r"[\s*,:)]*", after):
                return m.group(0)                        # the equation follows right after the number
            if after[:1] in "–-/" or before[-1:] == "(":
                return m.group(0)                        # a range like (4.49)–(4.52)
            done.add(n)
            star = "*" if (m.group(0).startswith("*") and not m.group(0).endswith("*")) else ""
            end = "*" if (m.group(0).endswith("*") and not m.group(0).startswith("*")) else ""
            return f"{star}({n}), ${EQ[n]}${end}"      # no italics round the maths
        out.append(_NUM.sub(rep_, seg))
    return "".join(out)


def finalize_equations() -> int:
    """Apply ``_show_all_eqs`` to every markdown cell; returns the number of cells changed."""
    changed = 0
    for c in nb.cells:
        if c.cell_type == "markdown":
            new = _show_all_eqs(c.source)
            if new != c.source:
                c.source = new
                changed += 1
    return changed


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + textwrap.dedent(text).strip())


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=r"""
Every flow obeys three bookkeeping rules: mass is neither made nor lost, momentum changes only when a force acts
(Newton), and energy is conserved (the first law). Written for a fixed set of fluid particles and carried to any control
volume by the Reynolds transport theorem, they become integral budgets that give forces from fluxes; shrunk to a point,
they become the continuity equation $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(Eq. 4.7)*, Cauchy's
equation $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(Eq. 4.24)* and the energy equation. A
Newtonian stress law turns Cauchy's equation into Navier–Stokes,
$\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(Eq. 4.39b)* — the equation of the rest of the
book. The chapter then shows what a rotating observer adds (Coriolis), how viscosity turns motion into heat, when
Bernoulli's equation holds, why the ocean and atmosphere can ignore density except where it meets gravity (Boussinesq),
what happens at the edges of a flow, and which dimensionless numbers decide whether two flows are the same flow.
""",
    roadmap=[
        "§4.1 integral vs differential laws; counting equations and unknowns",
        r"§4.2 mass: a moving control volume, $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ *(4.5)* (C01), and the continuity equation $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(4.7)* (C02)",
        "§4.3 stream functions: one scalar carries a 2-D flow (C03)",
        r"§4.4 momentum: forces from fluxes *(4.17)* (C04), Bernoulli along a streamline $\tfrac12U^2+gz+p/\rho=$ const *(4.19)* (C05), Cauchy's equation *(4.24)* (C06)",
        r"§4.5 the Newtonian stress law $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ *(4.31)* (C07)",
        r"§4.6 Navier–Stokes $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(4.39b)* (C08)",
        "§4.7 a rotating, accelerating observer: Coriolis and centrifugal terms *(4.45)* (C09)",
        r"§4.8 energy: dissipation $\varepsilon\ge0$ turns kinetic energy into heat *(4.57)*, *(4.58)* (C10)",
        "§4.9 Bernoulli's forms (C11, C12) and the Boussinesq approximation (C13)",
        r"§4.10 boundary and interface conditions, the kinematic condition $D\eta/Dt=0$ *(4.91)* (C14)",
        "§4.11 dimensionless Navier–Stokes *(4.101)* and dynamic similarity (C15)",
    ],
    prerequisites=[
        r"the Reynolds transport theorem and material volumes (Ch. 3 §3.6)",
        r"the material derivative $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ (Ch. 3 §3.2)",
        "strain rate S, vorticity ω and ∇·u as the volume growth rate (Ch. 3 §3.4)",
        r"stress tensor, Cauchy traction $f_j=n_i\tau_{ij}$, Gauss' theorem, the ε–δ identity (Ch. 2)",
        "Newton's viscosity law, Fourier's law, hydrostatics, first law and Gibbs relation, N², Π groups (Ch. 1)",
    ],
)
nb.explainer_index([
    ("control_volume_budgets", "Weigh a force by counting what flows through a box",
     "C01 C04: the moving-CV mass and momentum budgets — wake drag, bore, jet, rocket"),
    ("stream_function_spacing", "Can one number field hold a whole 2-D flow?",
     "C03: ψ contours are streamlines, their spacing is the speed, Δψ is the flux"),
    ("newtonian_stress_lab", "How does a fluid decide its stress?",
     "C07: G → S → τ, traction on a plane, bulk viscosity, the spinning cube"),
    ("navier_stokes_term_balance", "Which terms of Navier–Stokes are awake here?",
     "C08: the five terms of incompressible Navier–Stokes at a probe for six exact solutions"),
    ("rotating_frame_coriolis", "Why does a straight throw curve on a merry-go-round?",
     "C09: the rotating-frame terms, Coriolis deflection Ωut², effective gravity"),
    ("which_bernoulli", "Bernoulli is constant along what, exactly?",
     "C05 C11 C12: the four Bernoulli equations and their hypotheses"),
    ("viscous_dissipation_heating", "Where does the energy go when viscosity stops a flow?",
     "C10: dissipation ε ≥ 0, heating of a sheared channel, entropy production"),
    ("boussinesq_buoyancy", "Density hardly changes — so why does it drive the flow?",
     "C13: the Boussinesq momentum and heat equations and when they hold"),
    ("dynamic_similarity_models", "When does a model behave like the real thing?",
     "C15: dimensionless Navier–Stokes, Re, Fr, St and model testing"),
])
nb.setup()
nb.code(r"""
import numpy as np                                      # arrays (Ch. 1 primer P03)
import sympy as sp                                      # symbolic algebra (Ch. 1 primer P40)
import pandas as pd                                     # labelled tables (a DataFrame is a table with named rows/columns)
import matplotlib.pyplot as plt                         # static figures (Ch. 1 primer P01)
import plotly.graph_objects as go                       # rotatable 3-D figures (Ch. 1 P41, Ch. 2 P64)
from fluidpy import ch04_conservation_laws as ch04      # the tested chapter-4 module: every function cites its § and Eq.
from fluidpy import ch01_introduction as ch01           # Ch. 1 functions reused here (fluid properties, Newton's law, …)
from fluidpy import ch03_kinematics as ch03             # Ch. 3 functions reused as test fields (cylinder, vortices, paths)
from fluidpy.core.interact import slider_figure         # plotly figure with a slider that works on the web page (P17)
from fluidpy.core.anim import animate                   # matplotlib animations (P16); show_animation came with the setup
from fluidpy.core.style import COLORS, savefig          # the house palette and a helper that saves PNGs to outputs/ch04
from fluidpy.core import transport, grids, operators, tensors, kinematics, statics   # Ch. 2–3 machinery, called by module name
from fluidpy.core import integral_theorems as itg       # Ch. 2's Gauss and Stokes helpers
from tools.convergence import observed_order            # slope of log(error) vs log(step): the observed order (P13)
import logging                                          # standard library: controls library log messages
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless font-substitution notes


def recolor(fig, colors, dashes=None):                  # give plotly slider traces the notebook's colours by trace name
    for tr in fig.data:                                 # every trace of every slider step
        if tr.name in colors:                           # a name we assigned a colour to
            tr.line.color = colors[tr.name]             # same colour meaning as in the matplotlib figures
            tr.marker.color = colors[tr.name]           # markers too (for "markers" traces)
        if dashes and tr.name in dashes:                # optional dash pattern ("dash", "dot")
            tr.line.dash = dashes[tr.name]
    return fig                                          # the same figure, restyled


print(len([n for n in dir(ch04) if not n.startswith("_")]), "public names in fluidpy.ch04_conservation_laws")   # the toolbox
""", explain="""
1. Numerical, symbolic, table and plotting libraries (all primed in Ch. 1–2; `pandas` is used only for one small
   labelled table, the equation ledger below).
2. `ch04` is the chapter module; it re-exports the nine new core modules of this chapter (conservation budgets, stream
   functions, constitutive laws, Navier–Stokes residuals, rotating frames, curvilinear operators, Bernoulli,
   interfaces, similarity), so one name covers the whole chapter. Every function is tested in `tests/test_ch04.py`.
3. Ch. 1 and Ch. 3 functions are reused as test fields; the core modules of Ch. 2–3 are called by module name.
4. `recolor` only restyles plotly traces so that a colour always means the same thing.
""")
nb.md(r"""
**Notation, conventions and colours used in this notebook.**

| Symbol | Meaning | Unit |
|---|---|---|
| $\rho$ | density | kg/m³ |
| $\mathbf u=(u,v,w)$ | velocity | m/s |
| $\mathbf b$ | velocity of a control surface (the walls of our box) | m/s |
| $V(t)$, $A(t)$ | a material volume (always the same particles) and its surface | m³, m² |
| $V^*(t)$, $A^*(t)$ | a control volume (a box we choose) and its surface | m³, m² |
| $\mathbf n$ | outward unit normal | – |
| $\mathbf g$ | body force per unit mass; $z$ up, $\mathbf g=-g\mathbf e_z$, $\Phi=gz$ | m/s² |
| $\tau_{ij}$ | stress; **first index = the face's normal**, second = force direction | Pa |
| $p$, $\sigma_{ij}$ | pressure; viscous stress | Pa |
| $S_{ij}$, $\boldsymbol\omega$ | strain rate; vorticity | 1/s |
| $\mu$, $\nu=\mu/\rho$, $\mu_v$ | dynamic, kinematic and bulk viscosity | Pa s, m²/s, Pa s |
| $e$, $h$ | internal energy, enthalpy per unit mass | J/kg |
| $\mathbf q$, $k$ | heat flux; thermal conductivity | W/m², W/(m K) |
| $\varepsilon$ | viscous dissipation rate | W/kg |
| $\boldsymbol\Omega$ | angular velocity of a rotating frame | rad/s |
| $\dot\gamma$, $s$ | shear rate, extension rate | 1/s |

> ⚠️ **Four primes.** $\mathbf x'$, $\mathbf u'$ in §4.7 are measured in the **rotating** frame; $p'$ inside
> $\int_{p_o}^{p}dp'/\rho(p')$ in §4.9 is only an **integration variable**; $p'$, $\rho'$ later in §4.9 are **departures
> from the hydrostatic state**; and Ch. 3's $\mathbf x'$ was a **translating** frame. The code never reuses one name:
> `u_rot`, `x_rot` (rotating frame), `p_pert`, `rho_pert` (perturbations).

> ⚠️ **Symbols reused.** σ is the viscous stress (§4.5), the surface tension (§4.10) and a vortex core radius (§4.9);
> ε is the dissipation (§4.8) and the alternating tensor $\varepsilon_{ijk}$; γ is a material coefficient in the stress
> law (§4.5), the ratio $C_p/C_v$ and a scaled height in Ex. 4.7 — so we write the shear rate as $\dot\gamma$ and an
> extension rate as $s$. Φ, Ψ, φ, Ω, h, H, M, R and f also carry two or more meanings; each second meaning gets a ⚠️
> where it appears.

**Colours** (one meaning each, also in the explainers): mass or momentum **inflow blue**, **outflow orange**,
**storage purple**, **force on the fluid rose** · in Navier–Stokes: local $\partial\mathbf u/\partial t$ **blue**,
advective $(\mathbf u\cdot\nabla)\mathbf u$ **teal**, pressure **orange**, viscous **rose**, gravity **grey** ·
Coriolis **amber**, centrifugal **purple** · kinetic energy **teal**, heat and dissipation **rose** · in Bernoulli
stacks ½u² **teal**, p/ρ **orange**, gz **blue**, ∂φ/∂t **amber** · references and ghosts **grey dashed**.
""")
nb.md(r"""
**Where this chapter is used later**

| Result here | Used in |
|---|---|
| continuity $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ (4.7) | every flow of Ch. 5–16 |
| moving-box budgets (4.5) and $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ (4.17) | momentum integral (Ch. 9), layer budgets (Ch. 13), shocks (Ch. 15) |
| $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ (4.39b) | Ch. 5–13 |
| the Coriolis force $-2\boldsymbol\Omega\times\mathbf u'$ of (4.45) | geostrophy, Ekman layers, Rossby waves (Ch. 13) |
| dissipation $\varepsilon$ (4.58) | turbulence (Ch. 12) |
| $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=$ const (4.75) | water waves (Ch. 7) |
| Boussinesq $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ (4.86), $\frac{DT}{Dt}=\kappa\nabla^2T$ (4.89) | internal waves (Ch. 7), convection (Ch. 11), ocean and atmosphere (Ch. 13) |
| $D\eta/Dt=0$ on a surface (4.91) | every free surface |
| Re, Fr, Ri, Ro | the regime maps of every later chapter |

*Climate hook:* the equations of an ocean or atmosphere model are continuity (4.7), Navier–Stokes in a rotating frame
(4.45) with the Boussinesq simplification (4.86), and the heat equation (4.89) — this chapter writes each one down and
says when it holds.
""")

# =====================================================================================================================
# A.1 §4.1 Introduction — no A item
# =====================================================================================================================
nb.section("4.1", "Introduction", intro="""
**What is this section about?** The chapter's plan in two sentences: the same three principles (mass, momentum,
energy) are written twice — once as a budget for a finite volume (integral form) and once at a single point
(differential form) — and we keep a running count of equations against unknowns until the count closes.
""")
note("N01", "Integral vs differential form", r"""
An *integral* law is a budget for a finite region — what is inside, what flows through its surface, what forces act on
it. A *differential* law is the same statement at one point, a field equation. Gauss' theorem plus "true for every
volume" turns the first into the second (C02 does it for mass). Integral forms return in Ch. 9 (momentum integral),
Ch. 13 (layer budgets) and Ch. 15 (shock control volumes). *(N01 belongs to C01.)*
""")
nb.md(r"""
**The equation ledger (`N02`, filled in at C06, C08 and C10).** A set of equations can determine a flow only if there
are as many independent equations as unknown fields. We count as we go, starting from an empty table.
""")
nb.code(r"""
stages = ["cauchy", "navier_stokes", "barotropic", "full"]            # the four stages at which we will count
ledger = pd.DataFrame(0, index=stages, columns=["equations", "unknowns"])   # a labelled table of zeros
empty = ch04.closure_count("empty")                                   # the skeleton: nothing derived yet
print(empty["equations"], empty["unknowns"])                          # 0 0
print(ledger)                                                         # one row per stage, all zeros for now
""", explain="""
1. `pd.DataFrame(0, index=…, columns=…)` builds a small table of zeros with named rows (the stages) and columns.
2. `ch04.closure_count("empty")` is the skeleton; the same call with a stage name fills a row after C06, C08 and C10 —
   the table closes (7 equations = 7 unknowns) at C10.
""")
# =====================================================================================================================
# A.2 §4.2 Conservation of Mass — R01, C01, R02, C02
# =====================================================================================================================
nb.section("4.2", "Conservation of Mass", intro="""
**What is this section about?** Mass can be neither created nor destroyed. Written for a sealed bag of fluid that moves
with the flow, this is one line; the work is to rewrite it for any box we choose (moving, deforming or fixed) and then
for a single point.
""")
remind([
    ("functions as arguments and lambda", "`lambda x, t: …` is a one-line function; fluidpy takes fields as functions of (x, t) (Ch. 1 P29)."),
    ("f-strings", "`f\"{x:+.3f}\"` prints a number with a sign and 3 decimals (Ch. 1 P04)."),
    ("assert np.allclose", "`assert np.allclose(a, b)` stops the notebook if two results differ beyond rounding (Ch. 1 P15)."),
    ("path lines by solve_ivp", "`ch03.pathline(u, r0, t0, t_eval)` integrates $d\\mathbf r/dt=\\mathbf u(\\mathbf r,t)$ with `solve_ivp` (Ch. 3 P94)."),
    ("definite integral and np.trapezoid", "`np.trapezoid(y, x)` adds trapezoids under sampled values — a definite integral (Ch. 1 P27, P37)."),
    ("central finite differences in time", "$(F(t+\\Delta t)-F(t-\\Delta t))/2\\Delta t$ approximates $dF/dt$ with error ∝ Δt² (Ch. 1 P21)."),
    ("volume and surface integrals as midpoint sums", "a volume or surface integral is a sum of value × small volume (or area) over cells (Ch. 2 P83)."),
    ("animate and show_animation", "`animate(update, frames, fig)` builds a movie from an update function; `show_animation` plays it (Ch. 1 P16)."),
    ("matplotlib figures", "`fig, ax = plt.subplots()`, `ax.plot`, labels with units (Ch. 1 P01)."),
    ("numpy arrays", "vectors and tables of numbers that do arithmetic element by element (Ch. 1 P03)."),
])
nb.recap("R01", "Material volume and material surface", r"""
A *material volume* V(t) always contains the same fluid particles; its surface A(t) moves with the local fluid velocity,
so in the Reynolds transport theorem its surface velocity is b = u. Picture a perfectly sealed, infinitely flexible
balloon: it drifts, stretches and twists, but nothing crosses its skin. The theorem we use,
$\frac{d}{dt}\int_{V^*(t)}F\,dV=\int_{V^*(t)}\frac{\partial F}{\partial t}dV+\int_{A^*(t)}F\,\mathbf b\cdot\mathbf n\,dA$
*(Eq. 3.35)*, holds for any control volume; for the balloon, b = u.
""", where="Ch. 3 §3.6")
nb.code(r"""
u = lambda x, t: ch04.expanding_flow(x, t, a=1.0)[0]                 # 1-D test flow u = ax/(1 + at), a = 1 1/s [m/s]
ends = [ch03.pathline(u, [X], 0.0, [0.0, 1.0])[0, -1] for X in (1.0, 2.0)]   # follow the particles starting at 1 m and 2 m
print(np.round(ends, 8))                                            # [2. 4.]: after 1 s they sit at 2 m and 4 m
""", explain="""
Two particles that start at x = 1 m and x = 2 m end at 2 m and 4 m after 1 s: the material interval between them
(our 1-D balloon) stretches by the factor 1 + at.
""")

# ---- C01 ---------------------------------------------------------------------------------------------------------
core("C01", "Mass in a box that moves the way you choose (4.5)", """
A sealed balloon of air drifts and swells, so its mass never changes. A fixed box sits where the balloon passes. What
does the box record — and what would a box that slides along at its own speed record?
""", eqs=("4.5",))
problem("""
Engineers rarely follow a particular lump of fluid: they draw a box round a pump, a wing or a river reach and count what
enters and leaves. Meteorologists do the same with a grid cell, and oceanographers with a layer between two density
surfaces that moves up and down. We need one mass-conservation statement that works for *any* box — fixed, sliding,
stretching or riding with the fluid — and that reduces to "the balloon keeps its mass" when the box *is* the balloon.
""")
idea("""
material volume V(t)            control volume V*(t)                what the box sees
moves with the fluid (b = u)    moves as you choose (b = anything)   d/dt(mass inside) = −(net mass flux through
mass fixed: d/dt ∫ρ dV = 0      mass can change                                           its walls, RELATIVE to them)
                   the trick:  make V* and V coincide at one instant  →  (4.5)
""", r"""
The flux through a wall depends on how fast the fluid moves **relative to the wall**, $\mathbf u-\mathbf b$: a wall
that moves with the fluid lets nothing through, however fast both go.
""")
P("P111", "dataclasses and named results", r"""
A `dataclass` is a small Python class whose only job is to hold named fields; fluidpy returns every budget this way, so
you write `budget.residual` instead of remembering that the residual is the fourth number. `@dataclass(frozen=True)`
makes the result read-only (a result should not be edited by accident).
""", code="""
from dataclasses import dataclass                     # the decorator that builds the class
@dataclass(frozen=True)
class Budget:                                         # three named numbers
    storage: float; outflux: float; residual: float
b = Budget(storage=-1.0, outflux=1.0, residual=0.0)   # build one
print(b.residual, b)                                  # 0.0 Budget(storage=-1.0, outflux=1.0, residual=0.0)
""")
note("N03", "Mass of a material volume is constant", r"""
— our starting line. Number with our test flow $u=ax/(1+at)$, $\rho=\rho_0/(1+at)$ ($a$ = 1 s⁻¹, $\rho_0$ = 1 kg/m³, a
1-D flow, mass per unit cross-section area): the particles that occupy [1, 2] m at t = 0 occupy [2, 4] m at t = 1 s,
where ρ = 0.5 kg/m³ — mass 1 × 1 = 0.5 × 2 = 1 kg/m² both times.
""", equation=EQ["4.1"], ref="4.1")
remind([
    ("signed b·n and the swept volume of a moving wall", r"a wall moving at $\mathbf b$ sweeps volume at the rate $\mathbf b\cdot\mathbf n\,dA$ — positive when it moves outward (Ch. 3 P110)."),
    ("differentiation under the integral sign", r"$\frac{d}{dt}\int_a^bF\,dx=\int_a^b\frac{\partial F}{\partial t}dx$ when the limits are fixed; moving limits add boundary terms (Ch. 3 P109)."),
], lead="needed in the derivation below")
D("D01", ref="4.5")
note("N04", "The book's (4.2)", r"""
Step 2 of the derivation is the transport theorem with F = ρ and b = u. Correct but awkward: V(t) is carried by the
flow, which is usually what we do not know.
""", equation=EQ["4.2"], ref="4.2")
note("N05", "The book's (4.3)", r"Step 3 is the same theorem for a box of our choosing:", equation=EQ["4.3"], ref="4.3")
note("N06", "The instantaneously coincident control volume", r"""
Steps 4–6 are the book's (4.4): at the one instant when V* and V fill the same region, integrals of the same field over
the same region are equal — but the rates $\frac{d}{dt}\int\rho\,dV$ are **not** equal: a moment later the two volumes
have moved apart. That difference is exactly the flux term.
""", equation=EQ["4.4"], ref="4.4")
nb.md(r"""
> ⚠️ **Common confusion:** "$\frac{d}{dt}\int_{V^*}\rho\,dV$ and $\int_{V^*}\frac{\partial\rho}{\partial t}dV$ are the same
> thing." Only for a box that does not move or deform. The first is how fast the mass *in the box* changes; the second
> adds up how fast the density changes *at each fixed point*. For a moving box they differ by
> $\int_{A^*}\rho\,\mathbf b\cdot\mathbf n\,dA$ — the mass the moving walls sweep in or out (Ch. 3's signed swept volume).
""")
nb.worked_example("the expanding flow and three boxes, all at t = 0", r"""
Flow $u=ax/(1+at)$, $\rho=\rho_0/(1+at)$ with $a$ = 1 s⁻¹, $\rho_0$ = 1 kg/m³; box [1, 2] m.

1. **Fixed box** (b = 0): mass flux in at x = 1: ρu = 1 × 1 = 1 kg/(m² s); out at x = 2: 1 × 2 = 2 kg/(m² s); net
   outflux 2 − 1 = +1. By $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ *(4.5)*
   the mass inside changes at −1 kg/(m² s). Check: $\partial\rho/\partial t=-a\rho_0/(1+at)^2=-1$ kg/(m³ s) everywhere,
   times the length 1 m = −1 ✓.
2. **Material box** (ends move at the local u: $\dot x_0$ = 1, $\dot x_1$ = 2 m/s): u − b = 0 at both ends, no flux, mass
   constant — that is $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$ *(4.1)*.
3. **Box sliding at b = 1 m/s**: at x = 1, ρ(u − b) = 0; at x = 2, ρ(u − b) = 1 × (2 − 1) = 1; net outflux +1, so the
   mass inside falls at 1 kg/(m² s) — but for a new reason: its right wall moves slower than the fluid there.
""")
nb.code(r"""
print(ch04.material_interval(1.0, 2.0, 1.0))                      # where the particles of [1, 2] m are after 1 s: (2.0, 4.0)
print(ch04.material_mass(1.0, 2.0, 0.0), ch04.material_mass(1.0, 2.0, 1.0))   # their mass per area at t = 0 and 1 s [kg/m²]
boxes = {"fixed": (0.0, 0.0), "sliding b=1": (1.0, 1.0), "material": (1.0, 2.0)}   # wall speeds (left, right) [m/s]
for name, (d0, d1) in boxes.items():                                # the same flow, the same instant, three boxes
    bud = ch04.interval_mass_budget(1.0, 2.0, 0.0, dx0dt=d0, dx1dt=d1, flow="expanding", a=1.0, rho0=1.0)   # (4.5) in 1-D
    print(f"{name:12s} storage {bud['storage']:+.3f}  in {bud['flux_left']:+.3f}  out {bud['flux_right']:+.3f}"  # show the numbers
          f"  residual {bud['residual']:+.1e}")                     # storage + (out − in) = 0 every time [kg/(m² s)]
""", explain=r"""
1. `material_interval`: where the particles of [1, 2] m are after 1 s.
2. `material_mass`: their mass per unit area at t = 0 and t = 1 s — unchanged, $\frac{d}{dt}\int_{V(t)}\rho\,dV=0$ *(4.1)*.
3. `interval_mass_budget`: the budget (4.5), $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$,
   for three boxes on the same flow at the same instant — storage is measured by differencing the mass in the box in
   time, the fluxes by ρ(u − b) at the two ends; in every case storage + net outflux = 0.
""")
nb.md("**From scratch — two routes, one budget:** our own trapezoid sums of the mass in each box (differenced in time), "
      "then the library's 3-D quadrature on a box borrowed from Ch. 3.")
nb.check_agree(r"""
M = lambda t, x0, x1: np.trapezoid(ch04.expanding_flow(np.linspace(x0, x1, 2001), t)[1], np.linspace(x0, x1, 2001))   # mass per area in [x0, x1] at time t
dt = 1e-5                                                            # time step for the central difference [s]
fixed = (M(dt, 1, 2) - M(-dt, 1, 2)) / (2*dt)                       # fixed box: d/dt of the mass inside
sliding = (M(dt, 1 + dt, 2 + dt) - M(-dt, 1 - dt, 2 - dt)) / (2*dt)   # box whose walls move at b = 1 m/s
material = (M(dt, *ch04.material_interval(1, 2, dt)) - M(-dt, *ch04.material_interval(1, 2, -dt))) / (2*dt)   # walls ride on particles
print(round(fixed, 6), round(sliding, 6), round(material, 6))        # −1, −1, 0 kg/(m² s)
assert np.allclose([fixed, sliding, material], [-1.0, -1.0, 0.0], atol=1e-6)   # same numbers as interval_mass_budget above
cv = transport.MovingBox((1.0, 1.0, 1.0), origin=(1.0, 1.0, 1.0))    # Ch. 3's box: the unit cube [1, 2]³ m, not moving
rho3, u3 = ch04.expanding_flow_fields(dim=3)                         # the same flow in 3-D: u = a x/(1 + at), ρ = ρ0/(1 + at)³
B = ch04.mass_budget(rho3, u3, cv, 0.0, n=24 if not FAST else 16)   # (4.5) by midpoint quadrature over the six faces
print(B)                                                             # storage −3, outflux +3, residual ≈ 0 [kg/s]; local = nan: ∫∂ρ/∂t dV is only computed if you pass drho_dt
assert abs(B.residual) < 1e-6                                        # the library's 3-D budget closes too
assert abs(ch04.mass_budget(rho3, u3, cv, 0.0, material=True, n=16).outflux) < 1e-12   # b = u on the walls: nothing crosses
""")
nb.animation(r"""
a, frames = 1.0, (24 if not FAST else 12)                            # stretching rate [1/s]; number of frames
times = np.linspace(0.0, 1.2, frames)                                # the clock [s]
xg = np.linspace(0, 5, 400)                                          # the x-axis [m]
fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.8, 4.2), gridspec_kw=dict(height_ratios=[1.3, 1]), dpi=80)  # the figure and its panels
strip = ax.imshow(ch04.expanding_flow(xg, 0.0)[1][None, :], extent=(0, 5, 0, 1), aspect="auto",     # draw once; update() moves it
                  cmap="Blues", vmin=0.0, vmax=1.2)                   # density as a shaded strip (darker = denser)
ax.add_patch(plt.Rectangle((1, 0.05), 1, 0.9, fill=False, ec=COLORS["muted"], lw=2))   # the fixed box [1, 2] m
(bracket,) = ax.plot([1, 2], [0.5, 0.5], color=COLORS["accent"], lw=5, solid_capstyle="butt")   # the material interval
arrows = []                                                          # inflow/outflow arrows, redrawn every frame
txt = ax.text(0.02, 1.08, "", transform=ax.transAxes, fontsize=9)    # the running readout
ax.set_xlim(0, 5); ax.set_yticks([]); ax.set_xlabel("$x$ [m]")                                      # axis labels with units
bars = bx.bar(["storage", "net outflux", "sum"], [0, 0, 0], color=[COLORS["accent"], COLORS["orange"], COLORS["ink"]])  # three budget bars, heights set in update()
bx.set_ylim(-1.2, 1.2); bx.axhline(0, color=COLORS["muted"], lw=0.8); bx.set_ylabel("fixed box [kg/(m² s)]")  # axis labels with units


def update(i):                                                       # draw frame i
    t = times[i]                                                     # the time of this frame [s]
    strip.set_data(ch04.expanding_flow(xg, t, a=a)[1][None, :])      # density fades as 1/(1 + at)
    bracket.set_xdata(ch04.material_interval(1.0, 2.0, t, a=a))      # the balloon's ends ride on particles
    bud = ch04.interval_mass_budget(1.0, 2.0, t, flow="expanding", a=a)   # the fixed box's budget at time t
    for ar in arrows:                                                # remove last frame's arrows
        ar.remove()                                        # redraw the arrows each frame
    arrows.clear()                                        # redraw the arrows each frame
    arrows.append(ax.arrow(0.4, 0.25, 0.5*bud["flux_left"], 0, width=0.04, color=COLORS["blue"]))    # inflow ∝ ρu at x = 1
    arrows.append(ax.arrow(2.1, 0.25, 0.5*bud["flux_right"], 0, width=0.04, color=COLORS["orange"]))  # outflow ∝ ρu at x = 2
    for bar, v in zip(bars, [bud["storage"], bud["net_outflux"], bud["storage"] + bud["net_outflux"]]):  # storage, net outflux and their sum
        bar.set_height(v)                                            # the three budget bars
    txt.set_text(f"t = {t:.2f} s   mass in the material interval = {ch04.material_mass(1.0, 2.0, t, a=a):.3f} kg/m²")  # move this artist to the new frame
    return strip, bracket                                        # the artists that changed in this frame


show_animation(animate(update, frames=frames, fig=fig, interval=250), player="frames")   # step through with ◀ ▶
""")
see_read_change(
    see="""The purple bracket stretches to the right while the blue strip pales; the grey box keeps its place, blue
    arrow in, larger orange arrow out, and its bars keep adding to zero (black bar).""",
    read="""The balloon (bracket) keeps its mass because it grows exactly as fast as the density falls; the fixed box loses
    mass because more leaves through its right face than enters on the left. Both obey (4.5),
    $\\frac{d}{dt}\\int_{V^*}\\rho\\,dV+\\int_{A^*}\\rho(\\mathbf u-\\mathbf b)\\cdot\\mathbf n\\,dA=0$.""",
    change="""…a = −0.5 s⁻¹ (a converging flow): the bracket shrinks, the strip darkens, and the fixed box *gains* mass —
    its bars flip sign together.""")
whatif(r"""
…the box shrank to a point? The storage term would become ∂ρ/∂t times a tiny volume, the flux term a divergence times
the same volume — dividing by the volume leaves one equation at every point. That is C02.
""")

# ---- R02 + C02 ----------------------------------------------------------------------------------------------------
nb.recap("R02", "Gauss' divergence theorem", r"""
For a smooth vector field Q in a volume V with outward normal n on its surface A,
$\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA$ *(Eq. 2.30)*: the net outflow through the skin equals
the sum of the sources inside. Here Q = ρu, the mass flux per unit area.
""", where="Ch. 2 §2.12")
nb.code(r"""
Q = lambda X, Y, Z: (ch04.expanding_flow(np.array([X, Y, Z]), 0.0, dim=3)[1]
                     * ch04.expanding_flow(np.array([X, Y, Z]), 0.0, dim=3)[0])   # ρu of the 3-D expanding flow at t = 0 [kg/(m² s)]
print(itg.divergence_theorem_box(Q, [(1, 2), (1, 2), (1, 2)], n=24))   # (∫∇·Q dV, ∮Q·n dA) over the unit cube: both 3.0 kg/s
""", explain=r"""
Both sides of $\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA$ *(2.30)* for the mass flux in a unit
cube: the divergence of ρu is $3a\rho_0$ = 3 kg/(m³ s) at t = 0 and the cube's volume is 1 m³.
""")
core("C02", "The continuity equation: mass conservation at a point (4.7)", """
A crowd pours out of a square faster than people arrive. The square empties. What does "mass is conserved" say about one
single point of a flow?
""", eqs=("4.7",))
problem("""
A numerical weather model, a pipe-flow solver and a tidal model all store ρ and u at points of a grid. They need a rule
that ties the density's rate of change at a point to how the flow spreads out or converges there — a field equation,
not a box budget. Every exact solution in later chapters is checked against it.
""")
idea("""
box budget (4.5), fixed box       →  Gauss: surface flux = ∫ ∇·(ρu) dV     →  ∫ [∂ρ/∂t + ∇·(ρu)] dV = 0
true for EVERY box, however small  →  the integrand itself must vanish      →  ∂ρ/∂t + ∇·(ρu) = 0   (4.7)
""", r"""
**Flux divergence $\nabla\cdot(\rho\mathbf u)$** = the net outflow of mass per unit volume from a point: where it is
positive, the density there must fall.
""")
P("P112", "continuity of a function and the small-ball argument", r"""
A function f is *continuous* at x₀ if its values near x₀ are close to f(x₀): pick any margin, and a small enough ball
around x₀ keeps f within that margin. So if f(x₀) = 3, some small ball has f > 1.5 everywhere inside it — and the
integral of f over that ball is then at least 1.5 × (the ball's volume) > 0. That is the whole trick of the
"localisation lemma" below. (`ch04.ball_integral` integrates over a ball with `scipy.integrate.tplquad`, the triple
version of `quad`/`dblquad`.)
""", code="""
f = lambda X: 3.0 - 50*((X[0] - 0.5)**2 + X[1]**2 + X[2]**2)   # f(0.5, 0, 0) = 3, continuous; X = (x, y, z) [m]
for r in (0.2, 0.1, 0.05):                                     # shrink the ball [m]
    I = ch04.ball_integral(f, (0.5, 0.0, 0.0), r)              # ∫ f dV over the ball (scipy tplquad)
    print(r, round(I, 6), round(I/(4/3*np.pi*r**3), 4))        # the mean over the ball → f(x0) = 3 (it is 3 − 30 r²)
""")
gloss(["scipy.integrate.tplquad"])
remind([
    ("scipy.integrate.quad and dblquad", "adaptive 1-D and 2-D quadrature; `tplquad` adds a third, innermost variable (Ch. 3 P87)."),
    ("mean-value theorem for integrals", r"$\int_Vf\,dV=f(\mathbf x^*)\,V$ for some point $\mathbf x^*$ inside V when f is continuous (Ch. 2 P85)."),
    ("sympy symbols, Function, diff, simplify", "`sp.symbols`, `sp.Function('f')(x)` (an unknown function), `sp.diff`, `sp.simplify` (Ch. 1 P40)."),
    ("partial derivative", r"$\partial F/\partial x$ = the rate of change of F with x while the other variables are held fixed (Ch. 1 P25)."),
    ("tuple unpacking", "`rho, u = f(...)` names the parts of a returned pair (Ch. 1 P14)."),
])
note("N07", "Localisation lemma", r"""
If $\int_Vf\,dV=0$ for **every** volume V and f is continuous, then f = 0 at every point. Proof in the derivation below
(steps 5–7): a point where f ≠ 0 would have a small ball where f keeps its sign, and that ball's integral could not be
zero. We use it again for momentum (C06) and energy (C10).
""", equation=r"\int_Vf\,dV=0\ \ \forall V\ \Rightarrow\ f\equiv0")
D("D02", ref="4.7")
note("N08", "Flux divergence (transport) term", r"""
$\nabla\cdot(\rho\mathbf u)$ is the net mass outflow per unit volume at a point. Over a whole domain with no flux through
its boundary it integrates to zero (Gauss): transport moves mass around, it cannot change the total.
""")
remind([
    ("np.meshgrid and the project grid layout", "`grids.grid2d` returns `X, Y` arrays with x along axis 0 and spacing `h` (Ch. 2 P76)."),
    ("periodic grids and numpy stencils", "differences of shifted arrays (`np.roll` on a periodic grid) act on the whole grid at once — numpy broadcasting (Ch. 2 P77)."),
])
nb.code(r"""
n = 64 if not FAST else 48                                           # grid points per side
g = grids.grid2d(bounds=((0, 2*np.pi), (0, 2*np.pi)), n=(n, n), periodic=True)   # a periodic 2-D box [m]
rho = 1 + 0.2*np.sin(g.X)*np.cos(g.Y)                                # a lumpy density field [kg/m³]
U = np.array([np.sin(g.Y), np.cos(g.X)])                             # a smooth divergence-free velocity (u, v) [m/s]
flux_div = operators.divergence(rho*U, g.h, bc="periodic")           # ∇·(ρu) on the grid [kg/(m³ s)]
print(f"domain total {flux_div.sum()*g.h[0]*g.h[1]:.1e}   min {flux_div.min():+.3f}   max {flux_div.max():+.3f}")  # show the numbers
""", explain=r"""
$\partial\rho/\partial t=-\nabla\cdot(\rho\mathbf u)$ raises the density where the flux converges (negative divergence)
and lowers it elsewhere — local changes of ±0.15 kg/(m³ s) (the printed min and max) — while the domain total (the first number, ≈ 0 to rounding)
never changes: transport only moves mass around.
""")
P("P113", "product rule for a divergence", r"""
The divergence of a scalar times a vector splits like a product rule:
$\nabla\cdot(\rho\mathbf u)=\mathbf u\cdot\nabla\rho+\rho\,\nabla\cdot\mathbf u$ — in index form
$\partial(\rho u_i)/\partial x_i=u_i\,\partial\rho/\partial x_i+\rho\,\partial u_i/\partial x_i$ (the book's (B3.6)).
Here: the flux of mass changes because the density varies along the flow *or* because the flow spreads out.
""", code="""
x, y, z = sp.symbols('x y z'); rho = sp.Function('rho')(x, y, z)          # an unknown density field
u = [sp.Function(f'u{i}')(x, y, z) for i in (1, 2, 3)]; X = (x, y, z)     # an unknown velocity field
lhs = sum(sp.diff(rho*u[i], X[i]) for i in range(3))                      # ∇·(ρu)
rhs = sum(u[i]*sp.diff(rho, X[i]) + rho*sp.diff(u[i], X[i]) for i in range(3))   # u·∇ρ + ρ∇·u
print(sp.simplify(lhs - rhs))                                             # 0
""")
D("D03", ref="4.10")
note("N09", "Continuity with the particle's own rate", r"""
D03's step 3 is the book's (4.8). Read it as: the fractional rate at which a particle's density rises equals minus the
rate at which its volume grows (Ch. 3's $\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\nabla\cdot\mathbf u$ *(3.14)*).
Number: in the expanding flow at t = 0, $(1/\rho)D\rho/Dt=-a/(1+at)=-1$ s⁻¹ and ∇·u = +1 s⁻¹.
""", equation=EQ["4.8"], ref="4.8")
note("N10", "Incompressible flow", r"""
means each particle keeps its density. Different particles may have different densities: a stratified lake flowing
horizontally, $\mathbf u=(U(z),0,0)$, $\rho=\rho(z)$, has $D\rho/Dt=0+U\,\partial\rho/\partial x=0$ although ρ varies
with depth.

> ⚠️ §4.11 re-prints (4.9) as $\nabla\cdot\mathbf u=-\frac1\rho\frac{D\rho}{Dt}=-\frac{1}{\rho c^2}\frac{Dp}{Dt}$ — that is
> a general identity for any flow whose density changes are isentropic ($dp=c^2d\rho$), not the incompressibility
> condition (C15 uses it).
""", equation=EQ["4.9"], ref="4.9")
note("N11", "Continuity of an incompressible flow", r"""
With (4.9) in (4.8) it is simply the line below — zero volumetric strain rate (Ch. 3 §3.4): particles change shape but
not volume. Number: the Ch. 3 cylinder flow has ∇·u = 0 to round-off at every point off the body (code below).
""", equation=EQ["4.10"], ref="4.10")
note("N12", "Incompressible flow vs incompressible fluid", r"""
ρ = const everywhere is a special case of $D\rho/Dt=0$ *(4.9)* (every particle has the *same* density). A fluid whose
density does not respond to pressure is an incompressible *fluid*; a gas is not, but its *flow* is nearly incompressible
when the speed is small compared with the speed of sound c (Ch. 1 §1.9).

| case | $D\rho/Dt=0$? | ρ = const? |
|---|---|---|
| lake with a thermocline | yes | no |
| tap water in a pipe | yes | yes |
| air at 100 m/s | nearly (M = 0.29) | no |
| air at 300 m/s | no (M ≈ 0.88) | no |

Number: air at 288 K has c = 340.3 m/s, so U = 100 m/s gives M = U/c = 0.294 < 0.3 (the usual rule of thumb; C15
shows the departure from ∇·u = 0 scales with M²).
""")
nb.md(r"""
> ⚠️ **Common confusion:** "incompressible" means "constant density". It means *each particle* keeps its density,
> $D\rho/Dt=0$ *(4.9)*; the ocean's density varies by a few kg/m³ and its flow is still incompressible to excellent
> accuracy.
""")
nb.worked_example("continuity at one point of the expanding flow (t = 0, x = 1.5 m)", r"""
Flow $u=ax/(1+at)$, $\rho=\rho_0/(1+at)$, $a$ = 1 s⁻¹, $\rho_0$ = 1 kg/m³.

1. Local term: $\partial\rho/\partial t=-a\rho_0/(1+at)^2=-1$ kg/(m³ s).
2. Flux: $\rho u=\rho_0ax/(1+at)^2=1.5$ kg/(m² s); its x-derivative $\partial(\rho u)/\partial x=a\rho_0/(1+at)^2=+1$ kg/(m³ s).
3. Sum: −1 + 1 = 0 ✓ — $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(4.7)*.
4. Particle form: $u\,\partial\rho/\partial x=0$ (ρ does not vary in x), so $D\rho/Dt=-1$, $(1/\rho)D\rho/Dt=-1$ s⁻¹, and
   $\nabla\cdot\mathbf u=a/(1+at)=+1$ s⁻¹: $\frac1\rho\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0$ *(4.8)* holds. The flow is
   compressible ($D\rho/Dt\ne0$) — it is our test of the compressible form.
""")
nb.code(r"""
rho, u = ch04.expanding_flow_fields(a=1.0, rho0=1.0)                # ρ(x, t) and u(x, t) of the 1-D expanding flow
x0 = np.array([1.5])                                                 # the point x = 1.5 m
print(ch04.continuity_terms(rho, u, x0, 0.0, ht=1e-4))              # ∂ρ/∂t ≈ −1, ∂(ρu)/∂x ≈ +1, residual ≈ 0 [kg/(m³ s)]
print(ch04.continuity_material_terms(rho, u, x0, 0.0, ht=1e-4))     # (1/ρ)Dρ/Dt ≈ −1 and ∇·u ≈ +1 [1/s]: (4.8)
r3, u3 = ch04.stratified_shear_fields()                              # a stratified lake: u = U0 + shear·z, ρ = ρ0 + (dρ/dz) z
print(ch04.density_material_rate(r3, u3, np.array([0.3, 0.0, -2.0]), 0.0))   # Dρ/Dt = 0: incompressible (N10)
print(ch04.divergence_free_check(ch03.cylinder_velocity_field(1.0, 1.0), np.array([1.5, 0.7]), 0.0))   # (True, ~1e-9): N11
xs, ts, a = sp.symbols('x t a', positive=True)                       # symbols for the exact check
print(sp.simplify(ch04.continuity_residual_sym(1/(1 + a*ts), [a*xs/(1 + a*ts)], [xs], ts)))   # 0 exactly
print(ch04.is_incompressible_regime(100.0), round(ch04.mach_number(100.0, T=288.15), 4))   # True 0.2939 (N12)
""", explain=r"""
1. The two terms of $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(4.7)* by second-order stencils in x
   and t — they cancel to about 1e-8 (the stencil error).
2. The particle form $\frac1\rho\frac{D\rho}{Dt}+\nabla\cdot\mathbf u=0$ *(4.8)*: fractional density rate and divergence
   cancel.
3. A flow with varying density that is still incompressible (N10).
4. ∇·u = 0 for the Ch. 3 cylinder flow (N11).
5. sympy proves the residual is identically 0 for the test flow.
6. The Mach-number rule of thumb (N12): 100 m/s in air at 288 K is just inside it.
""")
nb.check_agree(r"""
hx, ht = 1e-4, 1e-4                                                  # space and time steps [m], [s]
drho_dt = (rho(x0, ht) - rho(x0, -ht)) / (2*ht)                     # ∂ρ/∂t by a central difference in time
flux = lambda x, t: rho(x, t)*u(x, t)[0]                             # mass flux ρu [kg/(m² s)]
dflux_dx = (flux(x0 + hx, 0.0) - flux(x0 - hx, 0.0)) / (2*hx)       # ∂(ρu)/∂x by a central difference in space
assert np.allclose(drho_dt + dflux_dx, 0.0, atol=1e-7)              # (4.7) holds at this point
T = ch04.continuity_terms(rho, u, x0, 0.0, ht=ht)                    # the library's two terms
assert np.allclose([drho_dt, dflux_dx], [T.local, T.flux_divergence], rtol=1e-7)   # same numbers → same stencils
wrong = lambda x, t: 1.0                                             # a WRONG density: constant, with the same u
print(ch04.continuity_residual(wrong, u, x0, 0.0))                   # +1: the equation catches an impossible field
""")
nb.md('**What does the code above do?** It rebuilds ∂ρ/∂t and ∂(ρu)/∂x at x = 1.5 m with our own central differences, checks they cancel and match `continuity_terms`, then feeds a wrong (constant) density: the residual +1 kg/(m³ s) shows the equation catching an impossible field.')
remind([("slider_figure", "`slider_figure(fn, name, values)` precomputes the curves for every slider position, so the figure works on the web page without Python (Ch. 1 P17).")])
nb.plotly(r"""
xg = np.linspace(0, 5, 41)                                           # 41 points along the line [m]
def terms(t):                                                        # the two terms of (4.7) at time t
    T = ch04.continuity_terms(rho, u, xg[None, :], t, ht=1e-4)      # stencils at all 41 points at once
    return {"∂ρ/∂t": (xg, T.local), "∂(ρu)/∂x": (xg, T.flux_divergence), "sum": (xg, T.residual)}   # the curves for this slider value
fig = slider_figure(terms, "t", np.linspace(0, 2, 21 if not FAST else 11), unit="s", xlabel="x [m]",  # precompute every slider position (works on the web page)
                    ylabel="[kg/(m³ s)]", title="Two mirror images that always add to zero", yrange=(-1.1, 1.1))  # labels, title and a fixed range
recolor(fig, {"∂ρ/∂t": COLORS["blue"], "∂(ρu)/∂x": COLORS["orange"], "sum": COLORS["ink"]}, {"sum": "dash"}).show()  # house colours for each trace, then draw
""")
see_read_change(
    see="A flat blue line below zero and a flat orange line above it; both shrink toward 0 as the slider time grows; the black dashed sum sits on zero.",
    read=r"""The density falls everywhere at the same rate because the flux diverges everywhere at the same rate; the sum
    is $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(4.7)* and stays at zero.""",
    change="…ρ₀ doubled: both lines double, the sum stays at zero — (4.7) is linear in ρ for a given u.")
nb.figure(r"""
z = np.linspace(-10, 0, 101)                                         # depth below the surface [m]
uz, rz = ch04.stratified_shear_flow(z, U0=0.2, shear=0.02, rho0=1000.0, drho_dz=-0.5)   # u(z) [m/s] and ρ(z) [kg/m³]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw=dict(width_ratios=[1, 1.6]))        # the figure and its panels
a.plot(rz, z, color=COLORS["blue"], lw=2.5, label="ρ(z)")          # denser water deeper
a.set_xlabel("ρ [kg/m³]"); a.set_ylabel("z [m]")                                        # axis labels with units
a2 = a.twiny(); a2.plot(uz, z, color=COLORS["teal"], lw=2.5); a2.set_xlabel("u [m/s] (teal)")   # sheared current
a.set_title("a lake: density and current vary with depth", fontsize=10)                             # the panel's message
Z0 = np.repeat(np.linspace(-9.5, -0.5, 8), 6)                        # 48 particles on 8 levels
X0 = np.tile(np.linspace(0, 5, 6), 8)                                # 6 per level [m]
u0, r0 = ch04.stratified_shear_flow(Z0, U0=0.2, shear=0.02, rho0=1000.0, drho_dz=-0.5)   # each particle's speed and density
for tt, al in ((0.0, 0.35), (60.0, 1.0)):                            # start and after 60 s
    b.scatter(X0 + u0*tt, Z0, c=r0, cmap="Blues", vmin=999, vmax=1005.5, s=40, alpha=al, edgecolors="k", linewidths=0.3)  # draw the data
b.set_xlabel("x [m]"); b.set_ylabel("z [m]")                                        # axis labels with units
b.set_title("after 60 s every particle keeps its colour (Dρ/Dt = 0)", fontsize=10)                  # the panel's message
savefig(fig, "ch04", "stratified_flow"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Horizontal layers sliding over each other (faded: start; solid: after 60 s); every dot keeps its colour.",
    read=r"""Each particle keeps its density ($D\rho/Dt=0$, incompressible) although density varies from layer to layer
    (ρ ≠ const): the surface layer moved farthest because the current is strongest there.""",
    change="""…a vertical velocity w pushed dense water up: ρ at a fixed height would rise (∂ρ/∂t > 0), but each particle's
    density would still be its own — Dρ/Dt = 0 is about particles, not places.""")
whatif(r"""
…the flow were steady and two-dimensional? Then (4.7) reduces to $\partial(\rho u)/\partial x+\partial(\rho v)/\partial y=0$
— two unknown fields tied by one equation. One scalar can then satisfy it automatically and carry the whole flow
picture: the stream function (C03).
""")
# =====================================================================================================================
# A.3 §4.3 Stream Functions — C03
# =====================================================================================================================
nb.section("4.3", "Stream Functions", intro="""
**What is this section about?** Steady continuity is a constraint on the velocity field. Writing the mass flux as a
cross product of two gradients satisfies it identically; in two dimensions a single function ψ remains, its contours are
the streamlines, and the flux between two of them is the difference of their labels.
""")
core("C03", "One number field for a whole 2-D flow: the stream function", """
On a weather map the streamlines crowd together where the wind is strong. Is that a drawing convention — or can one
scalar field carry both the direction and the speed of a 2-D flow?
""")
problem("""
A 2-D flow has two velocity components tied by continuity. If we could describe the whole flow by *one* function whose
contours are the streamlines, every picture would be a contour plot, and the continuity equation could never be violated
by a numerical slip. Potential-flow theory (Ch. 6), the Blasius boundary layer (Ch. 9) and quasi-geostrophic ocean and
atmosphere models (Ch. 13) are all written this way.
""")
idea("""
ψ = const   along every streamline            u = ∂ψ/∂y ,  v = −∂ψ/∂x    (ρ = 1)
ψ₂ − ψ₁   = volume flux between two streamlines (per metre of depth)
equal steps Δψ between contours  →  contours close together = fast flow   (u ≈ Δψ/Δn)
""")
note("N13", "Steady continuity", r"Start from continuity $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(4.7)* with ∂ρ/∂t = 0 — steady flow:",
     equation=EQ["4.11"], ref="4.11")
note("N14", "A mass flux written as a curl", r"""
A divergence of a curl is always zero ($\nabla\cdot\nabla\times\mathbf A=0$, Ch. 2 §2.13), so any mass flux written as a
curl satisfies (4.11) automatically. Choosing the vector potential as $\chi\nabla\psi$ and using $\nabla\times\nabla\psi=0$
gives the book's (4.12). (That every steady solenoidal ρu can be written like this — locally — is asserted by the book,
not proved; we only use the converse.) The sympy line below checks $\nabla\cdot(\nabla\chi\times\nabla\psi)=0$ for
arbitrary χ, ψ.
""", equation=r"\rho\mathbf u=\nabla\times\boldsymbol\Psi,\qquad\boldsymbol\Psi=\chi\nabla\psi\ \Rightarrow\ \rho\mathbf u=\nabla\chi\times\nabla\psi", ref="4.12")
nb.code(r"""
x, y, z = sp.symbols('x y z')                                        # Cartesian coordinates
chi = sp.Function('chi')(x, y, z); psi = sp.Function('psi')(x, y, z)  # two arbitrary smooth functions
m = ch04.mass_flux_from_stream_functions(chi, psi, (x, y, z))        # ρu = ∇χ × ∇ψ, three components
print(sp.simplify(sum(sp.diff(m[i], v) for i, v in enumerate((x, y, z)))))   # ∇·(ρu) = 0 identically
""")
remind([
    ("cross product with e_z; right-hand rule", r"$\mathbf e_z\times\mathbf e_x=\mathbf e_y$, $\mathbf e_z\times\mathbf e_y=-\mathbf e_x$ (Ch. 2 P74)."),
    ("chain rule along a curve dψ = ∇ψ·ds", r"along a curve $\mathbf x(s)$, $d\psi/ds=\nabla\psi\cdot d\mathbf x/ds$ (Ch. 3 P91)."),
    ("line integral along a path", r"$\int_1^2\mathbf F\cdot d\mathbf s$ adds the tangential part of F along a path (Ch. 2 P86)."),
], lead="needed in the derivation below")
D("D04")
nb.md(r"""
> ⚠️ **Sign convention.** With χ = −z the book (and this notebook) uses $u=\partial\psi/\partial y$,
> $v=-\partial\psi/\partial x$; many meteorology and oceanography texts use $u=-\partial\psi/\partial y$,
> $v=\partial\psi/\partial x$ (the flux is then ψ₁ − ψ₂). Same contours, opposite labels — check the convention before
> comparing numbers.
""")
nb.worked_example("the stagnation-point flow ψ = kxy with k = 1 s⁻¹", r"""
1. Velocity: $u=\partial\psi/\partial y=kx$, $v=-\partial\psi/\partial x=-ky$; at (1, 2) m: u = 1 m/s, v = −2 m/s.
2. Streamlines: xy = const — hyperbolas.
3. Flux between ψ = 1 and ψ = 2 m²/s: 2 − 1 = 1 m²/s per metre of depth, whichever path joins them.
4. Along the diagonal y = x the two contours cross at x = 1 m and $x=\sqrt2=1.414$ m, 0.586 m apart along the
   diagonal, where the speed is about $kx\sqrt2\approx1.7$ m/s — flux ≈ speed × spacing = 1.7 × 0.586 ≈ 1.0 ✓.
5. Where the contours are 0.1 m apart for Δψ = 1 m²/s, the speed is ≈ 10 m/s.
""")
nb.code(r"""
print(ch04.velocity_from_streamfunction_2d("stagnation", 1.0, 2.0, k=1.0))   # (u, v) = (1, −2) m/s from ψ = kxy
print(ch04.flux_between_streamlines("stagnation", (1.0, 1.0), (np.sqrt(2), np.sqrt(2)), k=1.0))   # straight gate: 1.0 m²/s
gate = np.array([[1.0, 1.2, 1.3, 1.414213562], [1.0, 0.9, 1.25, 1.414213562]])   # a bent gate with the same two ends
print(ch04.flux_along_path("stagnation", gate, k=1.0))              # 1.0 m²/s again: only the end labels matter
psi_cyl = lambda x, y: ch03.cylinder_streamfunction(x, y, 1.0, 1.0)  # Ch. 3's cylinder ψ (U = 1 m/s, a = 1 m)
print(ch04.velocity_from_streamfunction_2d(psi_cyl, 2.0, 0.5), ch03.cylinder_flow(2.0, 0.5, 1.0, 1.0))   # the same pair twice
X, Y, k = sp.symbols('x y k')                                        # symbols for the exact version
print(ch04.velocity_from_streamfunction_2d_sym(k*X*Y, X, Y))         # (k*x, -k*y)
""", explain=r"""
1. Velocity by central differences of ψ: $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$.
2. Flux through a straight gate = ψ₂ − ψ₁.
3. The same through a bent gate — the flux depends only on the end labels.
4. The cylinder's ψ from Ch. 3 gives back the Ch. 3 velocity (≈ (0.7924, −0.1107) m/s).
5. The symbolic version.
""")
nb.check_agree(r"""
h = 1e-5                                                             # finite-difference step [m]
psi = lambda x, y: ch04.streamfunction_preset("stagnation", x, y, k=1.0)   # ψ = kxy [m²/s]
u_fd = (psi(1, 2 + h) - psi(1, 2 - h)) / (2*h)                      # u = ∂ψ/∂y by a central difference
v_fd = -(psi(1 + h, 2) - psi(1 - h, 2)) / (2*h)                     # v = −∂ψ/∂x
assert np.allclose([u_fd, v_fd], [1.0, -2.0], rtol=1e-8)             # (1, −2) m/s, as the library said
p1, p2 = np.array([1.0, 1.0]), np.array([np.sqrt(2), np.sqrt(2)])     # the gate's ends (ψ = 1 and ψ = 2)
s = np.linspace(0, 1, 401)                                           # parameter along the gate
P = np.outer(p1, 1 - s) + np.outer(p2, s)                            # 401 points on the straight segment
t = (p2 - p1) / np.linalg.norm(p2 - p1); n = np.array([t[1], -t[0]])   # unit tangent and the normal to its right
un = ch04.velocity_preset("stagnation", *P, k=1.0)                   # exact (u, v) at those points
flux = np.trapezoid(un[0]*n[0] + un[1]*n[1], s) * np.linalg.norm(p2 - p1)   # ∫ u·n ds by trapezoids
print(flux)                                                          # 1.0 m²/s
assert np.isclose(flux, 1.0, rtol=1e-6)                              # = ψ₂ − ψ₁
""")
remind([
    ("plt.contour and streamplot", "`ax.contour(X, Y, F, levels)` draws level sets; `contourf` fills between them (Ch. 2 P78)."),
])
nb.figure(r"""
n = 241 if not FAST else 121                                         # grid points along x
xg, yg = np.linspace(-3, 3, n), np.linspace(-2, 2, (2*n)//3)          # the window [m]
Xg, Yg = np.meshgrid(xg, yg)                                         # 2-D grids (rows = y)
PSI = ch03.cylinder_streamfunction(Xg, Yg, 1.0, 1.0)                 # ψ of the cylinder flow, U = 1 m/s, a = 1 m [m²/s]
Uc, Vc = ch03.cylinder_flow(Xg, Yg, 1.0, 1.0)                        # its velocity [m/s]
inside = np.hypot(Xg, Yg) < 1.0                                      # points inside the body
speed = np.ma.masked_where(inside, np.hypot(Uc, Vc))                 # speed |u|/U, hidden inside the cylinder
fig, ax = plt.subplots(figsize=(8.5, 5.2))                                        # the figure and its panels
im = ax.pcolormesh(Xg, Yg, speed, cmap="viridis", vmin=0, vmax=2, shading="auto")   # speed heat map
ax.contour(Xg, Yg, np.ma.masked_where(inside, PSI), levels=np.arange(-2, 2.01, 0.2), colors="k", linewidths=0.7,
           negative_linestyles="solid")                               # equal Δψ = 0.2 (solid for ψ < 0 too)
ax.contourf(Xg, Yg, np.ma.masked_where(inside, PSI), levels=[0.5, 1.0], colors=[COLORS["teal"]], alpha=0.35)   # the band ψ ∈ [0.5, 1]
ax.contour(Xg, Yg, np.ma.masked_where(inside, PSI), levels=[0.5, 1.0], colors=[COLORS["teal"]], linewidths=2.5)   # its edges
y1, y2 = [(c + np.sqrt(c**2 + 4))/2 for c in (0.5, 1.0)]             # where ψ = 0.5 and 1 cross x = 0: y − 1/y = c
F = ch04.flux_between_streamlines(lambda x, y: ch03.cylinder_streamfunction(x, y, 1.0, 1.0), (0.0, y1), (0.0, y2))   # gate flux
ax.plot([0, 0], [y1, y2], color="k", lw=3)                           # the gate across the band
ax.annotate(f"flux = {F:.3f} m²/s", (0, y2), (0.4, 1.85), fontsize=10, arrowprops=dict(arrowstyle="->"))  # a labelled arrow
ax.add_patch(plt.Circle((0, 0), 1.02, color=COLORS["muted"]))       # the cylinder (drawn 2 % larger to hide the mask edge)
ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")                              # axis labels with units
fig.colorbar(im, ax=ax, label="speed |u|/U [–]")
ax.set_title("Equal steps of ψ: crowded contours = fast flow")                                      # the panel's message
savefig(fig, "ch04", "cylinder_psi"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="""Contours bunch together over the top and bottom of the cylinder and spread out in front and behind; the heat
map is brightest exactly where they bunch. The teal band between ψ = 0.5 and 1.0 m²/s narrows over the cylinder.""",
    read="""Equal steps in ψ carry equal flux, so a narrow gap means fast flow: at the top of the cylinder the speed is 2U
    and the gap is about half its far-field value. The black gate across the teal band carries exactly
    ψ₂ − ψ₁ = 0.5 m²/s.""",
    change="""…U doubled: every contour label doubles (ψ ∝ U) — with the same levels the drawing would show twice as many
    contours, each gap carrying the same 0.2 m²/s at twice the speed.""")
nb.plotly(r"""
import contourpy                                                     # matplotlib's contour engine: level sets as point lists
xs_, ys_ = np.linspace(-3, 5, 161), np.linspace(-2.5, 2.5, 101)      # the window [m]
Xs, Ys = np.meshgrid(xs_, ys_)                                       # 2-D grids
def level_lines(Z, levels):                                          # all contour pieces of several levels, NaN-separated
    gen = contourpy.contour_generator(xs_, ys_, Z)                   # a contour generator for this field
    px, py = [], []  # x and y of all contour pieces
    for c in levels:  # one contour level at a time
        for seg in gen.lines(c):                                     # each piece is an (n, 2) array of points
            px += list(seg[:, 0]) + [None]; py += list(seg[:, 1]) + [None]
    return px, py  # all pieces, ready for one plotly trace
def source_stream(m, U=1.0):                                         # source of strength m [m²/s] in a stream U [m/s]
    PS = ch04.streamfunction_preset("source_stream", Xs, Ys, U=U, m=m)   # ψ = Uy + (m/2π)θ
    th = np.linspace(1e-3, np.pi - 1e-3, 300)                        # polar angle along the body
    r = (m/2)*(1 - th/np.pi) / (U*np.sin(th)) if m > 0 else 0*th     # the dividing streamline ψ = ±m/2 in polar form
    bx = list(r*np.cos(th)) + [None] + list(r*np.cos(th)); by = list(r*np.sin(th)) + [None] + list(-r*np.sin(th))
    out = {"dividing streamline ψ = ±m/2": (bx, by) if m > 0 else ([-3, 5], [0, 0])}
    for k, dpsi in enumerate((0.3, 0.9)):                            # two families of streamlines outside the body
        out[f"ψ = ±(m/2 + {dpsi})"] = level_lines(PS, [m/2 + dpsi, -(m/2 + dpsi)])
    out["stagnation point"] = ([-m/(2*np.pi*U)], [0.0])              # where u = 0 on the axis
    return out  # the traces for this source strength
fig = slider_figure(source_stream, "m", np.linspace(0, 4, 17 if not FAST else 9), unit="m²/s", xlabel="x [m]",  # precompute every slider position (works on the web page)
                    ylabel="y [m]", title="A source in a stream: the half-body grows as m/(2U)", xrange=(-3, 5),  # labels, title and a fixed range
                    yrange=(-2.5, 2.5), modes={"stagnation point": "markers"})                      # labels, title and a fixed range
recolor(fig, {"dividing streamline ψ = ±m/2": COLORS["rose"], "ψ = ±(m/2 + 0.3)": COLORS["teal"],   # house colours for each trace, then draw
              "ψ = ±(m/2 + 0.9)": COLORS["blue"], "stagnation point": COLORS["ink"]}).show()        # colour of each named trace
""")
see_read_change(
    see="A rose half-body opening to the right, with teal and blue streamlines bending round it; the black stagnation point moves upstream as m grows.",
    read=r"""The dividing streamline (rose) encloses the source's fluid; its width far downstream is m/U — the flux m
    spread over the stream speed — and the stagnation point sits at $x=-m/(2\pi U)$ where the source's outflow just
    cancels the stream.""",
    change="…U doubled at fixed m: the body halves in width (m/U) and the stagnation point moves halfway back toward the source.")
note("N15", "Two stream functions", r"""
Since $\rho\mathbf u=\nabla\chi\times\nabla\psi$ is perpendicular to both gradients, the flow runs *along* every surface
χ = const and every surface ψ = const, so a 3-D streamline is where one of each meets (our version of Fig. 4.1, below).
""", equation=r"\rho\mathbf u\cdot\nabla\chi=\rho\mathbf u\cdot\nabla\psi=0")
remind([
    ("plotly 3-D surfaces and lines", "`go.Surface(x, y, z)` draws a sheet, `go.Scatter3d(mode='lines')` a curve; drag to rotate (Ch. 2 P64)."),
])
nb.plotly(r"""
chi_e, psi_e, (x, y, z) = ch04.stream_function_pair("parabolic")      # χ = y, ψ = z − x²
m_e = ch04.mass_flux_from_stream_functions(chi_e, psi_e, (x, y, z))   # ρu = ∇χ × ∇ψ = (1, 0, 2x)
print("ρu =", m_e)                                        # show the numbers
u3 = sp.lambdify((x, y, z), m_e, "numpy")                            # a numeric velocity from the sympy one
uf = lambda X, t: np.array([np.ones_like(X[0]), 0*X[0], 2*X[0]], dtype=float)   # the same field as fluidpy expects: u(x, t)
nn = 40 if not FAST else 24                                          # surface resolution
xs_ = np.linspace(-1, 1, nn); zz = np.linspace(-0.5, 2.5, nn)          # plotting ranges [m]
a_, b_, c_, d_ = 0.0, 1.0, 0.0, 1.0                                  # the four labels χ = a, b and ψ = c, d
fig = go.Figure()
for yv in (a_, b_):                                                  # the planes χ = y = a and y = b (blue)
    XX, ZZ = np.meshgrid(xs_, zz)                                        # 2-D grids of the coordinates
    fig.add_trace(go.Surface(x=XX, y=0*XX + yv, z=ZZ, opacity=0.25, colorscale=[[0, COLORS["blue"]], [1, COLORS["blue"]]], showscale=False, name=f"χ = {yv}"))  # add a trace to the plotly figure
for cv in (c_, d_):                                                  # the parabolic sheets ψ = z − x² = c and d (orange)
    XX, YY = np.meshgrid(xs_, np.linspace(-0.3, 1.3, nn))                                        # 2-D grids of the coordinates
    fig.add_trace(go.Surface(x=XX, y=YY, z=XX**2 + cv, opacity=0.35, colorscale=[[0, COLORS["orange"]], [1, COLORS["orange"]]], showscale=False, name=f"ψ = {cv}"))  # add a trace to the plotly figure
for yv in (a_, b_):                                                  # four streamlines traced numerically from the velocity
    for cv in (c_, d_):  # the two sheets ψ = c and ψ = d
        sl = ch03.streamline(uf, np.array([-1.0, yv, 1.0 + cv]), s_max=3.5, both=False)   # start on the intersection at x = −1
        fig.add_trace(go.Scatter3d(x=sl[0], y=sl[1], z=sl[2], mode="lines", line=dict(color=COLORS["ink"], width=6), showlegend=False))  # add a trace to the plotly figure
        assert np.allclose(sl[2] - sl[0]**2, cv, atol=1e-6) and np.allclose(sl[1], yv)   # it stays on both surfaces
yy, zp = np.meshgrid(np.linspace(a_, b_, 2), np.linspace(c_, d_, 2))  # the patch in the plane x = 0 bounded by the four sheets
fig.add_trace(go.Surface(x=0*yy, y=yy, z=zp, opacity=0.8, colorscale=[[0, COLORS["muted"]], [1, COLORS["muted"]]], showscale=False, name="patch"))  # add a trace to the plotly figure
fig.update_layout(height=520, title="Streamlines where χ = const meets ψ = const (flux through the grey patch = (b − a)(d − c))",  # axis titles and layout
                  scene=dict(xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]"))
fig.show()                                        # draw the interactive figure
""")
see_read_change(
    see="Blue planes and orange parabolic sheets crossing along parabolas; the black streamlines lie exactly on the crossings; a grey patch spans the tube between them.",
    read=r"""The flow $\mathbf u=(1,0,2x)$ climbs as it goes right; every streamline stays in its plane y = const and on
    its sheet $z-x^2$ = const (the cell asserts it).""",
    change="…χ = y replaced by χ = y + x: the planes tilt, the streamlines leave the planes y = const, and the flux through the patch stays (b − a)(d − c).")
note("N16", "Flux through a stream tube", r"""
The tube bounded by χ = a, b and ψ = c, d carries $\dot m=\oint_C\chi\,d\psi=b(d-c)+a(c-d)=(b-a)(d-c)$ — on the
ψ = const legs dψ = 0, on the χ = b leg the integral is b(d − c), on χ = a it is a(c − d) (Stokes' theorem, Ch. 2 §2.13).
""")
nb.code(r"""
print(ch04.stream_tube_mass_flux("parabolic", "parabolic", 0.0, 1.0, 0.0, 2.0))   # (numeric patch flux, (b − a)(d − c)) = (2.0, 2.0) kg/s
""")
remind([
    ("cylindrical unit vectors", r"$\mathbf e_R$ points away from the z-axis, $\mathbf e_z$ along it; $u_R$, $u_z$ are projections on them (Ch. 3 P88)."),
])
note("N17", "Axisymmetric flows", r"""
(no swirl, streamlines in planes through the z-axis) use χ = −φ. Test: ψ = ½UR² gives $u_z=U$, $u_R=0$ — a uniform
stream along the axis (code: `ch04.velocity_from_streamfunction_axisym("axisym_uniform", 0.5, 0.3, U=2.0)` → (0.0, 2.0),
printed below). Used for the sphere (Ch. 6) and pipe flow (Ch. 8); the flux between two stream surfaces is 2π(ψ₂ − ψ₁).
""", equation=r"\rho u_R=-\frac1R\frac{\partial\psi}{\partial z},\qquad\rho u_z=\frac1R\frac{\partial\psi}{\partial R}")
nb.code(r"""
print(ch04.velocity_from_streamfunction_axisym("axisym_uniform", 0.5, 0.3, U=2.0))   # (u_R, u_z) = (0, 2) m/s at R = 0.5 m
""")
note("N18", "Constant density", r"""
With constant density the same construction holds for u itself, and ψ differences are volume fluxes (m²/s per metre of
depth in 2-D) — every function here takes `rho=1.0` by default. Ch. 6 uses it throughout.
""")
remind([("show_viz", "embeds an interactive explainer; on the web page it fills the window (Ch. 1 P18).")])
nb.explainer("stream_function_spacing", heading="Can one number field hold a whole 2-D flow?", why=r"""
**Why interactive rather than a static figure:** a static contour plot shows one flow and one gate; here you move the gate, bend it and switch flows, and see the flux stay equal to ψ₂ − ψ₁ while the spacing–speed link updates under your cursor.

Seven classic flows drawn as ψ contours at equal steps over a speed map, with tracers moving along them. Drag a gate
across the flow: its flux stays equal to ψ₂ − ψ₁ however you tilt or bend it; click a point to read
$u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$ from the contour slopes.
""", tries=[
    "Choose 'stagnation' and drag the gate between the ψ = 1 and ψ = 2 contours along two different paths — the flux readout stays 1 m²/s.",
    "Choose 'cylinder': where do the contours crowd? Read the speed there in Explain.",
    "Choose 'line vortex': the contours are circles packed tighter toward the centre — the speed grows as 1/r.",
    "Open the Derivation tab and step through D04 with the uniform stream.",
])
whatif(r"""
…the flow were not steady? In 2-D incompressible flow ψ still exists at every instant (∇·u = 0 is enough), but its
contours are then the *instantaneous* streamlines, not paths (Ch. 3 §3.3). Knowing the flow pattern, the next question is
what forces it takes to maintain it — momentum (C04).
""")
# =====================================================================================================================
# A.4 §4.4 Conservation of Momentum — R03, R04, C04, C05, C06
# =====================================================================================================================
nb.section("4.4", "Conservation of Momentum", intro="""
**What is this section about?** Newton's second law for a material volume, carried to any control volume (forces from
fluxes: wakes, bores, rockets, sprinklers), Bernoulli's equation from a thin stream tube, and — shrunk to a point —
Cauchy's equation of motion, Newton's law for every continuum.
""")
nb.recap("R03", "Surface force from the stress tensor", r"""
On a surface element with unit normal n the force per unit area is $f_j=n_i\tau_{ij}$ *(Eq. 2.15)*: the **first** index
of τ is contracted with the normal. Its normal part is $\mathbf n\cdot\mathbf f=n_if_i$ and its tangential part
$f_k-(n_if_i)n_k$. For a fluid at rest $\tau_{ij}=-p\delta_{ij}$ and $\mathbf f=-p\mathbf n$, a push along the inward
normal.
""", where="Ch. 2 §2.6")
nb.code(r"""
tau = np.array([[-1e5, 3.0, 0], [3.0, -1e5, 0], [0, 0, -1e5]])        # pressure 1e5 Pa plus a 3 Pa shear pair [Pa]
n = np.array([1, 1, 0]) / np.sqrt(2)                                  # the unit normal of a 45° plane
print(tensors.traction(tau, n))                                       # f_j = n_i τ_ij ≈ (−70708.6, −70708.6, 0) Pa
print(tensors.normal_shear_stress(tau, n)[:2])                        # normal −99997.0 Pa, shear 0.0 Pa
""", explain="""
Traction and its normal and shear parts for pressure plus a small shear stress: on the 45° plane the shear pair adds
+3 Pa to the normal part (−100000 + 3) and leaves no shear — that plane is a principal plane of this τ.
""")
nb.recap("R04", "Surface tension acts through boundary conditions", r"""
Surface tension is a force per unit length acting along lines in an interface (Ch. 1 §1.6), not a force on the bulk
fluid, so it never appears in the field equations of this chapter; it enters through the conditions at a free surface or
an interface — derived from a force balance on a curved cap in C14.
""", where="Ch. 1 §1.6")

# ---- C04 ---------------------------------------------------------------------------------------------------------
core("C04", "Forces from fluxes: momentum in a moving box (4.17)", """
A bar is held in a wind tunnel on a thin sting. You may not touch the bar or the sting — you may only measure the air
speed behind it. Can you tell how hard the air pushes on the bar?
""", eqs=("4.17",))
problem("""
Measuring the force on a body directly is often impossible (a ship's hull, a turbine blade in a machine, the Earth's
surface under a wind). But the fluid leaving a region carries momentum away, and Newton's law says the momentum a box
loses per second must be supplied by forces. Count the momentum flowing in and out of a box drawn round the body and you
have weighed the force. The same bookkeeping gives the speed of a tidal bore, the thrust of a rocket and the torque of a
lawn sprinkler.
""")
idea("""
rate of change of momentum inside  +  momentum carried OUT through the walls (relative to them)
     d/dt ∫ ρu dV                   +  ∮ ρu (u − b)·n dA
=   body forces on the fluid       +  surface forces on the fluid
     ∫ ρg dV                        +  ∮ f dA                                   (4.17)
steady flow, fixed box:  (momentum out − momentum in) = forces on the fluid  → force on the body = −(that)
""")
note("N26", "Body forces vs surface forces", r"""

| body forces | surface forces |
|---|---|
| act without contact, ∝ mass | act by contact, ∝ area |
| given per unit mass (an acceleration, m/s²) | given per unit area (a stress, Pa) |
| gravity, electromagnetic; in accelerating frames the fictitious forces of §4.7 | pressure and viscous stress |

In $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$
*(4.17)* the first is $\int\rho\mathbf g\,dV$, the second $\int\mathbf f\,dA$ with $f_j=n_i\tau_{ij}$ (R03).
""")
P("P114", "momentum flux through a surface", r"""
Fluid crossing a surface carries its momentum with it. Through a patch dA the *volume* crossing per second is
$(\mathbf u-\mathbf b)\cdot\mathbf n\,dA$ (velocity relative to the patch, component along the normal); its mass is ρ
times that, and its momentum is u times that mass: $\rho\mathbf u\,[(\mathbf u-\mathbf b)\cdot\mathbf n]\,dA$ — a vector
(the direction of u) times a scalar (the rate of crossing). A jet of water 1 cm² in area at 10 m/s through a fixed
nozzle carries ρU²A = 1000 × 100 × 10⁻⁴ = 10 N of momentum per second.
""", code="""
rho, U, A = 1000.0, 10.0, 1e-4          # water [kg/m³], speed [m/s], area [m²]
mdot = rho*U*A                          # mass flux [kg/s]
print(mdot, mdot*U)                     # 1.0 kg/s carries 10.0 N of x-momentum per second
""")
remind([("Newton's second law and momentum", "force = rate of change of momentum, $\\mathbf F=d(m\\mathbf u)/dt$ (Ch. 1 P09).")])
note("N19", "Newton's second law for a material volume", r"""
— momentum per volume ρu, body force per mass g, surface force per area f, outward n:
""", equation=EQ["4.13"], ref="4.13")
D("D05", ref="4.17")
note("N20", "The book's (4.14)", r"""
Step 2 of D05 is the transport theorem applied to each component of ρu with b = u. It is also the starting line of
Cauchy's equation (C06).
""", equation=EQ["4.14"], ref="4.14")
note("N21", "The book's (4.15)", r"""
Step 3 is the same theorem for our box, rearranged.

> ⚠️ **Book slip:** the printed (4.15) ends with an extra "= 0". That would claim ∂(ρu)/∂t integrates to zero; (4.15)
> is an identity between three terms (compare $\frac{d}{dt}\int_{V^*(t)}\rho\,dV-\int_{V^*(t)}\frac{\partial\rho}{\partial t}dV-\int_{A^*(t)}\rho\,\mathbf b\cdot\mathbf n\,dA=0$
> *(4.3)*, where "= 0" is right because the mass statement was already substituted).
""", equation=EQ["4.15"], ref="4.15")
nb.note(r"""**The coincidence equalities** `N22` `N23` `N24` `N25` — step 4 is the book's (4.16a–d): at the instant
V* = V, the same integrand over the same region gives the same number — exactly the move of D01, now for a vector.""",
        equation=r"\begin{aligned}" + r"\\ ".join(EQ[k] + rf"\quad&\text{{({k})}}" for k in ("4.16a", "4.16b", "4.16c", "4.16d")) + r"\end{aligned}")
P("P115", "conservative force and its potential", r"""
A force field is *conservative* if the work it does on a particle moving from A to B does not depend on the route. Then
it is minus the gradient of a potential energy per unit mass Φ, and the work round any closed loop is zero. Gravity near
the ground: Φ = gz, force per mass −∇Φ = (0, 0, −g).
""", code="""
loop = itg.planar_loop(center=(0.0, 0.0, 0.0), normal=(1.0, 0.0, 0.0), radius=2.0)   # a vertical circle of radius 2 m
g = lambda X, Y, Z: np.array([0*X, 0*Y, -9.81 + 0*Z])        # uniform gravity field per unit mass [m/s²]
print(itg.circulation(g, loop))                               # work of g round the loop per unit mass: ≈ 0 J/kg
""")
note("N27", "Conservative body force", r"""
with z up. For gravity Φ = gz gives $\mathbf g=-g\mathbf e_z$; `ch04.body_force_from_potential("gravity", (0, 0, 5))` →
(0, 0, −9.81) m/s² (printed below). C11 uses Φ to fold gravity into the Bernoulli function.
""", equation=r"\mathbf g=-\nabla\Phi\quad\text{or}\quad g_j=-\partial\Phi/\partial x_j", ref="4.18")
nb.code(r"""
print(ch04.body_force_from_potential("gravity", np.array([0.0, 0.0, 5.0])))   # −∇(gz) at z = 5 m: (0, 0, −9.81) m/s²
""")
nb.md(r"""
> ⚠️ **Common confusion — whose force?** `N28` Drag $F_D$ is the force *on the body*, positive downstream. The forces in
> (4.17), $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$,
> are forces *on the fluid in the box*. By Newton's third law the body pushes the fluid with $-F_D\mathbf e_x$. Forget the
> minus and the drag comes out negative.
""")
nb.worked_example("drag of a bar from its wake (Ex. 4.1, our synthetic wake)", r"""
Far downstream the speed is $U(y)=U_\infty-\Delta\,e^{-y^2/b^2}$ with $U_\infty$ = 10 m/s, Δ = 2 m/s, b = 0.1 m; air
ρ = 1.2 kg/m³; the box is H = 2 m tall, its inlet far upstream (uniform $U_\infty$) and its outlet in the wake; pressure
is $p_\infty$ on all faces (so the pressure forces cancel — Gauss on a constant).

1. **Mass:** less leaves through the outlet than enters the inlet, by $\int(U_\infty-U)dy=\Delta b\sqrt\pi=2\times0.1\times1.7725=0.354$
   m²/s per metre of span — that must leave through the top and bottom (the Gaussian integral
   $\int_{-\infty}^{\infty}e^{-y^2/b^2}dy=b\sqrt\pi$).
2. That side outflow carries x-momentum $U_\infty$ per unit mass: $\rho U_\infty\times0.354$.
3. **Momentum x:** out − in $=\rho\int U^2dy+\rho U_\infty(0.354)-\rho\int U_\infty^2dy$, which must equal the force on the
   fluid, $-F_D/l$.
4. Collect: $F_D/l=\rho\int U(U_\infty-U)\,dy$ (N29).
5. Numbers: $\int U(U_\infty-U)dy=U_\infty\Delta b\sqrt\pi-\Delta^2b\sqrt{\pi/2}=3.5449-0.5013=3.0436$ m³/s²; × 1.2 =
   **3.65 N/m** of span.
""")
note("N29", "Ex. 4.1 result", r"""
(given; the derivation is the worked example above). It needs the top and bottom faces far enough out that shear and
pressure deviations there vanish; the integral of $U(1-U/U_\infty)/U_\infty$ is the *momentum thickness* of Ch. 9.
""", equation=r"F_D/l=\rho\int_{-H/2}^{+H/2}U(y)\big(U_\infty-U(y)\big)dy")
nb.code(r"""
y = np.linspace(-1, 1, 5)                                            # five heights across the box [m]
print(ch04.gaussian_wake(y, 10.0, 2.0, 0.1))                          # U(y): 10, 10, 8, 10, 10 m/s
print(ch04.wake_side_outflow("gaussian", 10.0, 2.0))                  # mass leaking through the sides: 0.3545 m²/s (U∞ = 10, H = 2 m)
print(ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0, return_error=True))   # (F_D/l, quadrature error) = (3.6523 N/m, ~1e-12)
sc = ch04.cv_scenario("wake", U_inf=10.0, deficit=2.0, width=0.1, rho=1.2, H=2.0)   # the same budget face by face
for f in sc["faces"]:                                                  # inlet, outlet, top, bottom
    print(f"{f['name']:6s} mass {f['mass_flux']:+9.4f} kg/(m s)   x-mom {f['momentum_flux_x']:+9.4f} N/m")   # out +, in −
print(sc["result"], sc["residual_mass"], sc["residual_momentum"])    # 3.6523 N/m, residuals ≈ 0
""", explain=r"""
1. Our synthetic wake profile.
2. The mass leaking through the sides, $\Delta b\sqrt\pi$.
3. The drag formula of N29 by adaptive quadrature with its error estimate.
4. The same budget face by face (signs: out positive, in negative) — per metre of span the inlet brings −240 N/m of
   x-momentum, the outlet takes +232.09 N/m and the sides +2.13 N/m each; the sum −3.652 N/m is the force on the fluid,
   $-F_D/l$, so the bar feels +3.652 N/m.
""")
nb.check_agree(r"""
y = np.linspace(-1, 1, 4001)                                         # fine grid across the box [m]
U = ch04.gaussian_wake(y, 10, 2, 0.1)                                # outlet profile [m/s]
m_in = -1.2*np.trapezoid(10*np.ones_like(y), y)                      # mass flux in through the inlet (negative) [kg/(m s)]
m_out = 1.2*np.trapezoid(U, y)                                        # mass flux out through the outlet
m_side = -(m_in + m_out)                                              # continuity: what is missing leaves through the sides
p_in = -1.2*np.trapezoid(100*np.ones_like(y), y)                     # x-momentum in: −ρU∞² H [N/m]
p_out = 1.2*np.trapezoid(U**2, y)                                     # x-momentum out through the outlet
p_side = m_side*10.0                                                  # side outflow carries U∞ per unit mass
FD = -(p_in + p_out + p_side)                                         # force on the bar = −(net momentum outflow)
print(FD)                                                             # 3.6523 N/m
assert np.isclose(FD, ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0), rtol=1e-8)   # same as the library
""")
nb.figure(r"""
sc = ch04.cv_scenario("wake", U_inf=10.0, deficit=2.0, width=0.1, rho=1.2, H=2.0)   # the four faces of the box
F = {f["name"]: f for f in sc["faces"]}                               # faces by name
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[1.5, 1]))        # the figure and its panels
a.add_patch(plt.Rectangle((0, -1), 4, 2, fill=False, ls="--", ec=COLORS["ink"]))   # the control volume (inlet x = 0, outlet x = 4)
a.add_patch(plt.Rectangle((0.9, -0.08), 0.2, 0.16, color=COLORS["muted"]))   # the bar
ys = np.linspace(-0.9, 0.9, 13)                                       # arrow rows
a.quiver(np.zeros_like(ys) - 0.9, ys, 10*np.ones_like(ys), 0*ys, color=COLORS["blue"], scale=80, width=0.004)   # uniform inflow U∞
yw = np.linspace(-0.9, 0.9, 37)                                       # finer rows in the wake
a.quiver(np.full_like(yw, 4.05), yw, ch04.gaussian_wake(yw, 10, 2, 0.1), 0*yw, color=COLORS["orange"], scale=80, width=0.004)   # outflow U(y)
for yy in (-1, 1):                                                    # small outflow through top and bottom
    a.annotate("", (2.0, yy + 0.35*np.sign(yy)), (2.0, yy), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=1.5))  # a labelled arrow
a.annotate("−F_D/l on the fluid", (0.3, 0.45), (1.6, 0.45), color=COLORS["rose"], fontsize=10, va="center",       # a labelled arrow
           arrowprops=dict(arrowstyle="-|>", color=COLORS["rose"], lw=2.5))   # the force on the fluid points upstream
a.set_xlim(-1, 5.4); a.set_ylim(-1.5, 1.5); a.set_xlabel("x [m]"); a.set_ylabel("y [m]")            # axis labels with units
a.set_title("(a) the box: uniform in, a narrow deficit out, a small side leak", fontsize=10)        # the panel's message
vals = [F["outlet"]["momentum_flux_x"] + F["inlet"]["momentum_flux_x"],        # outlet − inlet (net through the ends)
        F["top"]["momentum_flux_x"] + F["bottom"]["momentum_flux_x"],          # top + bottom
        -sc["result"]]                                                        # the sum = −F_D/l
b.bar(["outlet − inlet", "top + bottom", "sum = −F_D/l"], vals, color=[COLORS["orange"], COLORS["orange"], COLORS["rose"]])  # draw the data
for i, v in enumerate(vals):                                          # label each bar with its value
    b.text(i, v + 0.2*np.sign(v), f"{v:+.3f}", ha="center", fontsize=9)                             # a label on the plot
b.axhline(0, color=COLORS["muted"], lw=0.8); b.set_ylabel("x-momentum outflow [N/m]")               # axis labels with units
b.set_title(f"(b) inlet {F['inlet']['momentum_flux_x']:.0f}, outlet +{F['outlet']['momentum_flux_x']:.2f} N/m", fontsize=10)  # the panel's message
savefig(fig, "ch04", "wake_budget"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="A thin deficit in the orange outlet profile and a small but nonzero leak through the top and bottom; on the right, a negative net end flux, a positive side flux and their rose sum.",
    read="""The fluid leaves with less momentum than it brought; the difference, including what leaks out sideways at speed
    U∞, is the force on the fluid, −F_D/l = −3.652 N/m — so the bar is pushed downstream with 3.652 N per metre of span.""",
    change="…the box made taller: the side leak and every face flux change, the drag does not (it depends only on the deficit).")
nb.plotly(r"""
yy = np.linspace(-1.0, 1.0, 201)                                     # across the wake [m]
Uy = ch04.gaussian_wake(yy, 10.0, 2.0, 0.1)                           # U(y) [m/s]
def box(H):                                                           # what a box of height H sees
    inside = np.abs(yy) <= H/2                                        # the part of the wake inside the box
    return {"U(y) [m/s]": (yy, Uy),                                        # the curves for this slider value
            "integrand ρU(U∞ − U) inside the box [N/m²]": (yy[inside], (1.2*Uy*(10.0 - Uy))[inside]),  # another named curve (x, y)
            "faces y = ±H/2": ([-H/2, -H/2, None, H/2, H/2], [0, 11, None, 0, 11])}                 # another named curve (x, y)
fig = slider_figure(box, "H", np.linspace(0.1, 2.0, 20 if not FAST else 10), unit="m", xlabel="y [m]", ylabel="U [m/s] · integrand [N/m²]",  # precompute every slider position (works on the web page)
                    title="A box must enclose the whole wake (F_D/l → 3.652 N/m)", xrange=(-1, 1), yrange=(0, 11))  # labels, title and a fixed range
recolor(fig, {"U(y) [m/s]": COLORS["teal"], "integrand ρU(U∞ − U) inside the box [N/m²]": COLORS["rose"],  # house colours for each trace, then draw
              "faces y = ±H/2": COLORS["muted"]}, {"faces y = ±H/2": "dash"}).show()                # colour of each named trace
for H in (0.1, 0.2, 0.3, 0.5, 2.0):                                   # the drag a box of height H would report
    print(f"H = {H:.1f} m: F_D/l = {ch04.wake_drag_per_span('gaussian', 10.0, 1.2, H):.4f} N/m")    # show the numbers
""")
see_read_change(
    see="The teal wake profile and the rose integrand; the grey faces close in as H shrinks and cut off the rose curve's tails.",
    read="""Below H ≈ 3b = 0.3 m the box cuts the wake and misses part of the deficit (the printed F_D/l falls below 3.652
    N/m); above it F_D/l stops changing — and the formula's assumptions (uniform pressure, no shear on the sides) also
    need the faces outside the wake.""",
    change="…a wider wake (b = 0.2 m): the box must be twice as tall before the drag stops changing.")
note("N31", "Ex. 4.3 — a small bore (surge) with a moving control volume", r"""
Dimensional analysis first (the Π theorem of Ch. 1, R12): with g, h and the speed U there is one group, $U^2/(gh)$, so
$U\propto\sqrt{gh}$. Riding with the wave (b = U e_x) the flow in the box is steady; hydrostatic pressure on the vertical
faces, $p_o$ on the free surface (its net push cancels), mass and momentum give the result below. Number: still water
$h_{in}$ = 1 m ahead of the wave, $h_{out}$ = 1.1 m behind it → U = √(9.81 × 1.1 × 2.1/2) = 3.37 m/s, close to
√(gh) = 3.13 m/s (the shallow-water wave speed of Ch. 7 and Ch. 13).
""", equation=r"U^2\frac{h_{in}}{h_{out}}=\frac g2(h_{in}+h_{out}),\qquad U=\sqrt{\frac{gh_{out}}{2h_{in}}(h_{in}+h_{out})}\approx\sqrt{gh}")
remind([("Python dictionaries", "`{'name': value}` maps names to values and `d['name']` reads one — fluidpy returns many results this way (Ch. 1 P23)."),
        ("np.random.default_rng", "`np.random.default_rng(seed)` makes reproducible random numbers (Ch. 1 P10)."),
        ("hydrostatic pressure on a vertical face ½ρgh²", r"$\int_0^h\rho g(h-y)\,dy=\tfrac12\rho gh^2$ per metre of width (Ch. 1 P28).")])
nb.code(r"""
print(ch04.bore_speed(1.0, 1.1), np.sqrt(9.81*1.0))                  # 3.3661 m/s vs √(gh) = 3.1321 m/s
for p_o in (1.0e5, 0.0):                                              # two atmospheric pressures [Pa]
    print(ch04.bore_pressure_force(1.0, 1.1, p_o=p_o))                # end faces and top [N per m width]; `net` does not change
""", explain=r"""
1. The bore speed from the momentum budget and the shallow-water estimate √(gh).
2. The horizontal pressure forces on the moving box: each end face carries $p_o h$ plus the hydrostatic
   $\tfrac12\rho gh^2$ (4905 N and 5935 N per metre width), the sloping top carries $p_o(h_{out}-h_{in})$; the net,
   −1030 N/m, does not depend on $p_o$.
""")
nb.animation(r"""
h_in, h_out, g0 = 1.0, 1.1, 9.81                                      # depths ahead of / behind the front [m]
U = ch04.bore_speed(h_in, h_out)                                     # front speed [m/s], moving to the left
V = U - ch04.bore_outlet_velocity(h_in, h_out)                       # lab speed of the water behind the front [m/s]
frames = 60 if not FAST else 30                                      # number of frames
T = 3.0                                                              # duration [s]
x = np.linspace(0, 16, 400)                                          # the channel [m]
x_front = lambda t: 12.0 - U*t                                       # front position [m]
eta = lambda t: h_in + (h_out - h_in)*0.5*(1 + np.tanh((x - x_front(t))/0.3))   # surface: 1 m ahead (left), 1.1 m behind
rng = np.random.default_rng(0)  # reproducible random parcel positions (Ch. 1 P10)
xp0 = rng.uniform(0.5, 15.5, 60); yp = rng.uniform(0.1, 0.9, 60)      # parcels [m]
fig, ax = plt.subplots(figsize=(7, 2.8), dpi=80)                                        # the figure and its panels
(surf,) = ax.plot(x, eta(0), color=COLORS["blue"], lw=2)             # the free surface
dots = ax.scatter(xp0, yp, s=10, color=COLORS["teal"])                # water parcels
(box,) = ax.plot([], [], "k--", lw=1.2)                              # the control volume riding with the front
ax.set_xlim(0, 16); ax.set_ylim(0, 1.35); ax.set_xlabel("x [m]"); ax.set_ylabel("height [m]")       # axis labels with units
ax.set_title(f"A bore moving left at U = {U:.2f} m/s; in the box's frame the flow is steady", fontsize=10)  # the panel's message


def update(i):                                                       # draw frame i
    t = i*T/frames                                                   # time [s]
    surf.set_ydata(eta(t))                                           # move the front
    t_pass = (12.0 - xp0)/U                                          # when the front reached each parcel
    xp = np.where(t > t_pass, xp0 - V*(t - t_pass), xp0)             # at rest ahead; moving left at V behind
    dots.set_offsets(np.c_[xp, yp])                                  # move the parcels
    xf = x_front(t)
    box.set_data([xf - 1.5, xf + 1.5, xf + 1.5, xf - 1.5, xf - 1.5], [0.02, 0.02, 1.25, 1.25, 0.02])   # the riding box
    return surf, dots, box                                        # the artists that changed in this frame


show_animation(animate(update, frames=frames, fig=fig, interval=50))                                # build the movie from update() and play it
""")
see_read_change(
    see="The dashed box slides left with the front; parcels ahead of it are at rest, parcels behind it creep slowly left.",
    read=f"""Seen from the box, water enters fast (at U = 3.37 m/s) and 1 m deep and leaves slower and 1.1 m deep — a
    steady budget whose solution is the front's speed. Behind the front the water moves left at only
    U(1 − h_in/h_out) ≈ 0.31 m/s.""",
    change="…h_out → h_in: the front becomes a small wave and its speed tends to √(gh).")
note("N32", "Ex. 4.4 — a rocket as an accelerating control volume", r"""
(b = b(t) e_z). Exhaust leaves the nozzle at $V_e$ relative to the rocket; with $F_S$ the sum of drag and pressure
thrust, mass and vertical momentum give the pair below. With g = $F_S$ = 0 it integrates to Tsiolkovsky's
$\Delta b=V_e\ln(M_0/M_1)$ (our extension; ln is the natural logarithm): $M_0$ = 2 kg, $M_1$ = 1 kg, $V_e$ = 800 m/s
→ 554.5 m/s.
""", equation=r"\frac{dM}{dt}+\rho_eV_eA_e=0,\qquad M\frac{d^2z_R}{dt^2}=-V_e\frac{dM}{dt}-Mg+F_S")
remind([
    ("natural logarithm", r"$\ln$ undoes $e^x$; $\int dM/M=\ln M$ (Ch. 1 P36)."),
    ("separation of variables (rocket, tank, meniscus)", r"move all of one variable to one side and integrate each side (Ch. 1 P42)."),
])
nb.figure(r"""
r = ch04.rocket_trajectory(2.0, 0.1, 800.0, 10.0, t_eval=np.linspace(0, 10, 201))   # M0 = 2 kg, 0.1 kg/s, Ve = 800 m/s, 10 s burn
ghost = 800.0*np.log(2.0/r["M"]) - 9.81*r["t"]                        # Tsiolkovsky minus the gravity loss gt [m/s]
fig, ax = plt.subplots(figsize=(6.5, 3.6))                                        # the figure and its panels
ax.plot(r["t"], r["b"], color=COLORS["accent"], lw=2.5, label="b(t) from solve_ivp")       # the rocket's speed
ax.plot(r["t"], ghost, "--", color=COLORS["muted"], label=r"$V_e\ln(M_0/M)-gt$")          # closed form
ax.set_xlabel("t [s]"); ax.set_ylabel("speed b [m/s]"); ax.legend()                                 # axis labels with units
ax.set_title("A rocket speeds up faster as it gets lighter")                                        # the panel's message
print(f"burn-out speed {r['b'][-1]:.1f} m/s; ΔV without gravity {ch04.rocket_delta_v(2.0, 1.0, 800.0):.1f} m/s")  # show the numbers
savefig(fig, "ch04", "rocket"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="A purple curve bending upward, with the grey dashed closed form on top of it.",
    read="The thrust V_e|dM/dt| = 80 N is constant but it pushes less and less mass; the numerical and closed-form speeds agree: 456.4 m/s at burn-out, 554.5 m/s without gravity.",
    change="…V_e doubled: every speed doubles minus the gravity loss gt, which does not change.")
P("P116", "torque, moment arm and moment of inertia", r"""
A force F applied at position r (from an axis point) turns a body with torque M = r × F; only the part of F perpendicular
to r, times the distance, counts (the moment arm). The angular analogue of mass is the moment of inertia: for a cube of
side h and density ρ spinning about an axis through its centre, I = ρh⁵/6. Rate of change of angular momentum = net
torque. (`np.cross` is the cross product of Ch. 2 on arrays.)
""", code="""
r = np.array([0.2, 0.0, 0.0])                  # lever 0.2 m along x
F = np.array([0.0, 5.0, 0.0])                  # 5 N pushing along y
print(np.cross(r, F))                          # torque (0, 0, 1.0) N m about z
rho, h = 1000.0, 0.01; print(rho*h**5/6)       # I of a 1 cm water cube: 1.67e-08 kg m²
""")
gloss(["np.cross"])
note("N82", "Angular momentum of a rigid body", r"For a rigid body the rate of change of angular momentum equals the applied torque. Its fluid version for a fixed control volume follows.",
     equation=EQ["4.64"], ref="4.64")
note("N83", "Angular momentum for a stationary control volume", r"""
(the transport theorem with F = r × ρu, the same move as D05). Internal torques cancel in pairs — because the stress is
symmetric, which C07 proves.
""", equation=EQ["4.65"], ref="4.65")
note("N84", "Ex. 4.6 — a lawn sprinkler held still", r"""
Two jets of area A at radius a leave at speed U, tilted by α from the tangential direction; only the tangential part
U cos α has a moment arm, so the torque needed to hold the arm is as below. Number: a = 0.2 m, water, A = 1 cm², U = 5 m/s,
α = 30° → M = 2 × 0.2 × 1000 × 10⁻⁴ × 25 × 0.866 = 0.866 N m. Let go of it and (with no friction) it spins until the
jets leave with no tangential speed, at Ω = U cos α/a = 21.7 rad/s (our extension).
""", equation=r"M=2a\rho AU^2\cos\alpha")
nb.code(r"""
print(ch04.sprinkler_torque(0.2, 1000.0, 1e-4, 5.0, np.pi/6), ch04.sprinkler_free_spin_rate(0.2, 5.0, np.pi/6))   # 0.866 N m, 21.65 rad/s
""")
nb.explainer("control_volume_budgets", heading="Weigh a force by counting what flows through a box", why=r"""
**Why interactive rather than a static figure:** the figure above shows one box; the lesson of a control volume is that you may choose the box — its size, its speed, which fluid it holds — and only by moving it yourself do you see which face fluxes change while the force stays put.

One scene, one control volume, every face's flux as a bar: storage + outflow − inflow on the left of
$\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ (4.17),
forces on the right, and a residual that stays at zero. Four scenarios (wake, bore, jet, rocket) and three box types
(fixed, riding, material). (For a box of finite height the Gaussian wake's integrals use the error function
erf, the finite-limit version of $\int e^{-y^2/b^2}dy=b\sqrt\pi$.)
""", tries=[
    "In 'wake', shrink the box height below three wake widths — the status turns amber and the drag readout drops: the box cuts the wake.",
    "In 'bore', drag the box speed b until the storage bar vanishes: that speed is the bore's speed (Explain shows the formula).",
    "Switch the budget to mass and choose the material box (b = u): every flux bar vanishes — that is (4.1), $\\frac{d}{dt}\\int_{V(t)}\\rho\\,dV=0$.",
    "Open the Derivation tab and step through D01: the coincidence step freezes the balloon on top of the box.",
])
whatif(r"""
…the box were a thin stream tube in a steady, frictionless flow? Its side walls carry no flux, only pressure; the
momentum budget along it becomes a relation between speed, pressure and height — Bernoulli's equation (C05).
""")

# ---- C05 ---------------------------------------------------------------------------------------------------------
core("C05", "Bernoulli along a streamline: fast flow, low pressure (4.19)", """
An airliner measures its airspeed with a small tube pointing into the wind and a hole in its side — no moving parts. How
does a pressure difference become a speed?
""", eqs=("4.19",))
problem("""
Pitot tubes on aircraft and in rivers, the speed of water leaving a hole in a tank, the suction on the top of a wing: all
trade pressure for speed along the path of the fluid. The trade follows from momentum conservation on a very thin stream
tube in a steady, frictionless, constant-density flow.
""")
idea("""
along ONE streamline, per unit mass:   ½U²  +  g z  +  p/ρ  = constant           (4.19)
                                       kinetic  potential  "pressure energy"
speed up (U ↑) → pressure down (p ↓) at the same height;  stop the flow (U = 0) → pressure up by ½ρU²
""")
P("P117", "sympy expand, series, removeO, collect and subs", r"""
Four sympy moves for approximations: `expand` multiplies brackets out, `series(expr, ds, 0, 2)` keeps powers of ds up to
ds¹ and adds an O(ds²) marker that `removeO()` drops, `collect(expr, ds)` groups terms by powers of ds, and
`subs({a: b})` substitutes. Here: keep only first-order terms of a small element of length ds.
""", code="""
ds, U, dU = sp.symbols('ds U dU')                               # a small length and two symbols
e = sp.expand((U + dU*ds)**2)                                   # U**2 + 2*U*dU*ds + dU**2*ds**2
print(e.series(ds, 0, 2).removeO())                             # U**2 + 2*U*dU*ds (ds² dropped)
print(sp.collect(e, ds), e.subs({U: 2, dU: 1, ds: 0.1}))        # grouped by powers of ds; 4.41
""")
remind([
    ("first-order Taylor expansion", r"$f(s+ds)\approx f(s)+f'(s)\,ds$ for a small step ds (Ch. 1 P26)."),
    ("limits and orders of smallness", r"terms in $ds^2$ vanish faster than ds as ds → 0 and drop out after dividing by ds (Ch. 2 P68)."),
])
note("N30", "Ex. 4.2 builds (4.19) from a stream-tube element", r"""
of length ds: mass with first-order changes of U and A, the streamwise momentum with gravity (sin θ ds = dz) and the
extra pressure force on the slowly widening side, then drops the (ds)² terms to get the line below. The derivation fills
in the moves the book skips (where the side pressure force comes from, why the mean area multiplies gravity).
""", equation=r"U\frac{\partial U}{\partial s}ds=-g\,dz-\frac1\rho\frac{\partial p}{\partial s}ds")
D("D06", ref="4.19", check_src=r"""
ds, U, A, p, rho, g = sp.symbols('ds U A p rho g', positive=True)   # length, speed, area, pressure, density, gravity
Us, As, ps, th = sp.symbols('U_s A_s p_s theta')                 # U_s = ∂U/∂s, A_s = ∂A/∂s, p_s = ∂p/∂s (any sign), tilt θ
mass = -rho*U*A + rho*(U + Us*ds)*(A + As*ds)                                   # step 3: mass out − mass in
mom = -rho*U**2*A + rho*(U + Us*ds)**2*(A + As*ds)                              # step 4: momentum out − momentum in
forces = (-rho*g*sp.sin(th)*(A + As*ds/2)*ds + p*A                             # step 5: weight on the mean area, inlet push,
          + (p + ps*ds/2)*As*ds - (p + ps*ds)*(A + As*ds))                      #         side push, outlet push
mass1 = sp.expand(mass).coeff(ds, 1)                                            # first-order part of the mass balance
As_sol = sp.solve(sp.Eq(mass1, 0), As)[0]                                       # it fixes A_s = −A U_s/U
balance = sp.series(sp.expand(mom - forces), ds, 0, 2).removeO()               # steps 6–8: keep terms up to ds¹
balance = sp.simplify(balance.subs(As, As_sol) / (rho*A*ds))                   # step 9: use mass, divide by ρA ds
print(balance)                                                                   # U*U_s + g*sin(theta) + p_s/rho = 0 is (4.19) per unit length
assert sp.simplify(balance - (U*Us + g*sp.sin(th) + ps/rho)) == 0               # step 10: d(½U²)/ds + g dz/ds + (1/ρ)dp/ds
lib = ch04.stream_tube_element_balance_sym()                                    # the library's own re-run of Ex. 4.2
print(lib["result"], lib["residual_vs_4_19"])                                   # the same expression, residual 0
""")
nb.md(r"""
> ⚠️ **Book slip:** the statement of Ex. 4.2 writes the result as "½ρU² + gz + p/ρ = constant"; the ρ in the first term is
> a misprint (J/m³ added to J/kg). The derived $\tfrac12U^2+gz+p/\rho=$ constant along a streamline *(4.19)* is per unit
> mass; multiplied by ρ it reads $p+\tfrac12\rho U^2+\rho gz=$ constant (pressures).
""")
note("N104", "Stagnation and dynamic pressure", r"""
Where a streamline meets a body head-on the fluid stops (a *stagnation point*). Between a far point (p, U) and the
stagnation point (p₀, 0) at the same height, $\tfrac12U^2+gz+p/\rho=$ const *(4.19)* gives the line below: p₀ is the
*stagnation* (total) pressure and ½ρU² the *dynamic* pressure.
""", equation=r"p_0=p+\tfrac12\rho\lvert\mathbf u\rvert^2")
note("N103", "Pitot tube", r"""
(numbers in the worked example and code below): the tube's mouth is a stagnation point, the side hole reads the static
pressure; a manometer between them gives the speed. Number: Δp = 500 Pa in air (ρ = 1.2 kg/m³) → 28.9 m/s (104 km/h).
Only two points on one streamline are used — irrotationality is not needed (the book's wording says "irrotational", a
stronger assumption than required).
""", equation=r"\lvert\mathbf u\rvert_1=\sqrt{2(p_2-p_1)/\rho}=\sqrt{2g(h_2-h_1)}")
note("N105", "Orifice (Torricelli)", r"""
Water leaving a hole a depth h below the free surface of a large tank: at the surface U ≈ 0, p = p_atm; in the jet
p = p_atm (straight, parallel streamlines have no pressure difference across them). $A_c$ is the area of the jet at its
narrowest (the *vena contracta*, where the jet stops contracting) — about 0.61 of a sharp-edged hole's area (a measured
value), close to 1 for a rounded one. Number: h = 1 m → 4.43 m/s.
""", equation=r"u=\sqrt{2gh},\qquad\dot m=\rho A_c\sqrt{2gh}")
nb.worked_example("pitot tube and tank", r"""
1. Pitot: Δp = 500 Pa, ρ = 1.2 kg/m³: $U=\sqrt{2\times500/1.2}=\sqrt{833.3}=28.87$ m/s; dynamic pressure
   ½ × 1.2 × 28.87² = 500 Pa ✓.
2. Tank: h = 1 m, g = 9.81 m/s²: $u=\sqrt{2\times9.81\times1}=4.43$ m/s; a 1 cm² sharp hole ($C_c$ = 0.61) passes
   $\rho A_cu$ = 1000 × 0.61 × 10⁻⁴ × 4.43 = 0.270 kg/s.
3. Both from $\tfrac12U^2+gz+p/\rho=$ constant *(4.19)* between two points of one streamline.
""")
nb.code(r"""
print(ch04.bernoulli_head(28.87, 0.0, 1.0e5, 1.2), ch04.bernoulli_head(0.0, 0.0, 1.0e5 + 500, 1.2))   # ½U² + gz + p/ρ at two points [m²/s²]
print(ch04.pitot_speed(1.0e5 + 500.0, 1.0e5, 1.2), ch04.dynamic_pressure(28.8675, 1.2),             # show the numbers
      ch04.stagnation_pressure(1.0e5, 28.8675, 1.2))                   # 28.8675 m/s, 500 Pa, 100500 Pa
print(ch04.torricelli_speed(1.0), ch04.orifice_mass_flow(1.0, 1e-4, Cc=0.61))   # 4.4294 m/s, 0.2702 kg/s
print(ch04.bernoulli_solve({"U": 0.0, "z": 1.0, "p": 1.0e5}, {"z": 0.0, "p": 1.0e5}, "U"))   # solve (4.19) for the jet speed
""", explain=r"""
1. The Bernoulli constant $\tfrac12U^2+gz+p/\rho$ at the free stream and at the stagnation point agree (≈ 8.375×10⁴ m²/s²).
2. The pitot inversion and the two pressures.
3. Torricelli and the jet's mass flow.
4. The general solver: solve $\tfrac12U^2+gz+p/\rho=$ const *(4.19)* for the one unknown (the jet speed 1 m below the
   surface).
""")
remind([("explicit Euler stepping (tank draining)", r"$h_{n+1}=h_n+\Delta t\,f(h_n)$ marches an ODE forward in small steps (Ch. 1 P30).")])
nb.check_agree(r"""
At, Ao, g0, dt = 1.0, 1e-3, 9.81, 0.1                                # tank and hole areas [m²], gravity, time step [s]
h, t = 1.0, 0.0                                                      # start full (1 m) at t = 0
while h > 0:                                                         # explicit Euler on dh/dt = −(Ao/At)√(2gh)
    h -= dt*(Ao/At)*np.sqrt(2*g0*max(h, 0.0))
    t += dt
t_closed = ch04.tank_drain(1.0, At, Ao)["t_empty"]                    # closed form (At/Ao)√(2h0/g) = 451.5 s
print(round(t, 1), round(t_closed, 1))                                # the two emptying times [s]
assert abs(t - t_closed) < 1.0                                       # agree to within a few time steps
""")
nb.figure(r"""
th = np.linspace(0, np.pi, 181)                                     # angle round the cylinder from the front [rad]
us, vs = ch03.cylinder_flow(-(1 + 1e-9)*np.cos(th), (1 + 1e-9)*np.sin(th), 1.0, 1.0)   # velocity on the surface r = a (front at θ = 0) [m/s]
p_s = 0.0 + 0.5*1.0*(1.0**2 - (us**2 + vs**2))                      # Bernoulli (4.19) from the free stream, ρ = 1 [Pa]
Cp = ch04.pressure_coefficient(p_s, 0.0, 1.0, 1.0)                   # C_p = (p − p∞)/(½ρU²)
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6))                                        # the figure and its panels
a.plot(np.rad2deg(th), Cp, color=COLORS["teal"], lw=2.5)              # the surface pressure
a.plot([0, 180], [1, 1], "o", color=COLORS["ink"]); a.plot([90], [-3], "o", color=COLORS["rose"])   # draw the data
a.annotate("stagnation points: C_p = 1", (0, 1), (30, 0.2), fontsize=9)                             # a labelled arrow
a.annotate("shoulder: C_p = −3", (90, -3), (100, -2.4), fontsize=9)                                 # a labelled arrow
a.set_xlabel("angle from the front θ [°]"); a.set_ylabel("$C_p$ [–]")                               # axis labels with units
a.set_title("(a) cylinder: $C_p=1-4\\sin^2\\theta$", fontsize=10)                                   # the panel's message
for Cc, ls in ((1.0, "--"), (0.611, "-")):                          # rounded vs sharp-edged hole
    d = ch04.tank_drain(1.0, 1.0, 1e-3, Cc=Cc, t_eval=np.linspace(0, 740, 300))   # h(t) for a 1 m² tank, 10 cm² hole
    b.plot(d["t"], d["h"], ls, color=COLORS["blue"], label=f"$C_c$ = {Cc}")                         # draw the data
b.set_xlabel("t [s]"); b.set_ylabel("level h [m]"); b.legend(); b.set_title("(b) a draining tank", fontsize=10)  # axis labels with units
savefig(fig, "ch04", "bernoulli_uses"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="(a) The pressure is highest where the flow stops and lowest where it is fastest; (b) a parabola in t that reaches zero sooner without the vena contracta.",
    read=r"""(a) $C_p=1-(U/U_\infty)^2$ is $\tfrac12U^2+gz+p/\rho=$ const *(4.19)* in dimensionless form; (b) the outflow
    slows as √h, so the level falls fast at first and slowly at the end (451.5 s to empty with $C_c$ = 1, 739 s with 0.611).""",
    change="…viscosity added: behind a real cylinder the flow separates and the pressure never recovers to C_p = 1 (Ch. 9) — Bernoulli fails where friction acts.")
whatif(r"""
…we asked for momentum at every *point* rather than along one streamline in frictionless flow? We would need the stresses
on a tiny cube — pressure and friction together — and Newton's law for it: Cauchy's equation (C06). Bernoulli's *general*
forms return in C11 and C12, where the `which_bernoulli` explainer compares all four.
""")

# ---- C06 ---------------------------------------------------------------------------------------------------------
core("C06", "Newton's law for a fluid particle: Cauchy's equation (4.24)", """
A small cube of fluid is pushed by its neighbours on six faces and pulled by gravity. Before we know anything about what
kind of fluid it is, what is F = ma for it?
""", eqs=("4.24",))
problem(r"""
The integral law (4.17) needs a box; to predict a flow everywhere we need the law at every point — a differential
equation for u. We can get it without any knowledge of how stress depends on motion: that part comes later (C07). The
result holds for water, air, honey, even a steel beam.
""")
idea("""
(4.14) for a material volume  →  Gauss on BOTH surface integrals  →  one volume integral = 0 for every volume
                              →  localise  →  flux form (4.22)  →  subtract u_j × continuity  →  ρ Du_j/Dt = ρ g_j + ∂τ_ij/∂x_i
mass × acceleration          =   gravity   +   net push of the neighbours (divergence of stress, FIRST index)
""")
remind([("stress tensor τ_ij, first index = face", r"$\tau_{ij}$ is the j-component of the force per area on a face whose normal points along $x_i$ (Ch. 1 P05, R03).")])
P("P118", "tensor divergence over the first index", r"""
The net surface force per unit volume on a small cube is the divergence of the stress. Because the traction on a face
with normal n is $f_j=n_i\tau_{ij}$ (first index = face), Gauss' theorem turns $\oint n_i\tau_{ij}dA$ into
$\int\partial\tau_{ij}/\partial x_i\,dV$ — the derivative acts on the **first** index. For a symmetric τ it does not
matter; for a non-symmetric one it does, and C07 has not yet proved τ symmetric. (Ch. 2's `tensor_divergence`
differentiated the second index by default; here we pass `index=0`.)
""", code="""
tau_ns = lambda x, t: np.array([[0.0, x[0], 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])   # τ12 = x1, not symmetric [Pa]
first, second = ch04.divergence_first_index_demo(tau_ns, np.array([0.3, 0.2, 0.1]))   # ∂τ_ij/∂x_i and ∂τ_ij/∂x_j
print(first, second)                                                                   # [0. 1. 0.] vs [0. 0. 0.]: different
""")
D("D07", ref="4.24")
note("N33", "The book's (4.20a)", "D07's step 2: Gauss applied to each component j of the momentum flux, the normal's index contracted.",
     equation=r"\int_{A}\rho\mathbf u(\mathbf u\cdot\mathbf n)dA=\int_{V}\nabla\cdot(\rho\mathbf u\mathbf u)dV=\int_{V}\frac{\partial}{\partial x_i}(\rho u_iu_j)dV", ref="4.20a")
note("N34", "The book's (4.20b)", "Step 3: the same for the surface force — the first index again.",
     equation=r"\int_{A}\mathbf f\,dA=\int_{A}n_i\tau_{ij}dA=\int_{V}\frac{\partial\tau_{ij}}{\partial x_i}dV", ref="4.20b")
note("N35", "The book's (4.21)", "Step 4: everything under one integral.", equation=EQ["4.21"], ref="4.21")
note("N36", "The conservative (flux) form", "Step 5 (localisation) gives the form finite-volume codes discretise (Ch. 10):",
     equation=EQ["4.22"], ref="4.22")
note("N37", "The bracket is continuity", r"""
Steps 6–8 expand the left side; the bracket is the continuity equation $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$
*(4.7)* and vanishes. The sympy cell below shows the flux form minus $\rho\,Du_j/Dt$ is exactly $u_j$ × (continuity).
""", equation=r"\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho\frac{\partial u_j}{\partial t}+u_j\Big[\frac{\partial\rho}{\partial t}+\frac{\partial}{\partial x_i}(\rho u_i)\Big]+\rho u_i\frac{\partial u_j}{\partial x_i}=\rho\frac{Du_j}{Dt}", ref="4.23")
nb.code(r"""
x, y, z, t = sp.symbols('x y z t')                                   # coordinates and time
r = sp.Function('rho')(x, y, z, t)                                   # an arbitrary density field
Uf = [sp.Function(f'u{i}')(x, y, z, t) for i in (1, 2, 3)]           # an arbitrary velocity field
diff_ = ch04.conservative_to_advective_sym(r, Uf, (x, y, z), t)      # flux form − ρ Du_j/Dt, one component at a time
cont = sp.diff(r, t) + sum(sp.diff(r*Uf[i], v) for i, v in enumerate((x, y, z)))   # ∂ρ/∂t + ∂(ρu_i)/∂x_i
print([sp.simplify(d - Uf[j]*cont) for j, d in enumerate(diff_)])    # [0, 0, 0]: the difference is u_j × continuity
""")
nb.md(r"""
> ⚠️ **Book slip (index):** the sentence after (4.24) calls the net surface force ∂τ_ij/∂x_j. The equation itself,
> $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(4.24)*, contracts the first index. The two agree only
> because τ is symmetric — proved in C07, *after* (4.24).
""")
note("N02", "The ledger, first row", r"""
Continuity (1 equation), Cauchy (3) and two thermodynamic equations of state give 6 equations; the unknowns are ρ (1),
u_j (3) and τ_ij (9): 13. Not solvable yet — a *constitutive law* must tie τ to the motion (C07).
""")
nb.code(r"""
c = ch04.closure_count("cauchy")                                     # equations and unknowns after Cauchy's equation
ledger.loc["cauchy"] = [c["equations"], c["unknowns"]]                # fill the first row of the §4.1 table
print(c["equation_names"], c["unknown_names"])                       # what is being counted
print(ledger)                                                        # cauchy: 6 equations, 13 unknowns
""")
remind([("rigid (solid-body) rotation and its pressure", r"$\mathbf u=\boldsymbol\Omega\times\mathbf x$ turns every element without deforming it (Ch. 3 P101).")])
nb.worked_example("water spinning like a solid body (Ω = 1 rad/s) at r = 0.1 m", r"""
Velocity $u_\varphi=\Omega r=0.1$ m/s; every particle moves on a circle, so its acceleration is centripetal:
$\mathbf a=-\Omega^2r\,\mathbf e_R=-0.1$ m/s², ρa = −100 N/m³. The pressure that holds it on the circle is
$p=\rho\Omega^2R^2/2-\rho gz$ (plus a constant): $-\partial p/\partial R=-\rho\Omega^2R=-100$ N/m³ radially, and
$-\partial p/\partial z=+\rho g$ vertically balancing gravity $\rho g_z=-\rho g$. No viscous stress (a rigid motion has
S = 0, Ch. 3). So $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(4.24)* with $\tau_{ij}=-p\delta_{ij}$
reads −100 = 0 + (−100) radially and 0 = −9810 + 9810 vertically ✓.
""")
nb.code(r"""
u_sb, p_sb = ch04.exact_field("solid_body", Omega=1.0, rho=1000.0)   # u = Ω × x and p = ρΩ²R²/2 − ρgz
tau_sb = lambda x, t: -p_sb(x, t)*np.eye(3)                          # stress of a fluid with no viscous part [Pa]
x0 = np.array([0.1, 0.0, 0.0])                                       # a point 0.1 m from the axis
gvec = np.array([0.0, 0.0, -9.81])                                   # gravity [m/s²]
T = ch04.cauchy_terms(1000.0, u_sb, tau_sb, gvec, x0, 0.0, h=1e-4)    # the three terms of (4.24) per unit volume [N/m³]
print(T.inertia, T.body, T.stress_divergence, T.residual)            # (−100, 0, 0), (0, 0, −9810), (−100, 0, 9810), ≈ 0
print(ch04.momentum_conservative_residual(1000.0, u_sb, tau_sb, gvec, x0, 0.0))   # the flux form (4.22) balances too
""", explain=r"""
1. The exact solid-body field and its pressure.
2. The stress of a fluid with no viscous part, $\tau_{ij}=-p\delta_{ij}$.
3. The three terms of $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(4.24)* at a point — mass ×
   acceleration equals gravity plus the net push.
4. The flux form $\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$
   *(4.22)* gives the same balance (continuity holds).
""")
nb.check_agree(r"""
div = np.zeros(3)                                                    # ∂τ_ij/∂x_i, built by hand
for i in range(3):                                                   # differentiate row i along x_i
    e = np.eye(3)[i]*1e-4                                            # a small step along x_i [m]
    div += (tau_sb(x0 + e, 0)[i, :] - tau_sb(x0 - e, 0)[i, :]) / 2e-4
print(div)                                                           # (−100, 0, 9810) N/m³
assert np.allclose(div, T.stress_divergence, rtol=1e-6)              # the library contracts the FIRST index too
first, second = ch04.divergence_first_index_demo(tau_ns, np.array([0.3, 0.2, 0.1]))   # the primer's non-symmetric τ
assert not np.allclose(first, second)                                # the second-index version would be a different force
""")
nb.md('**What does the code above do?** It differentiates each row i of τ along x_i by hand (the first-index divergence), matches `cauchy_terms`, and confirms that for the non-symmetric τ of the primer the second-index divergence would give a different force.')
nb.figure(r"""
Rs = np.array([0.05, 0.10, 0.15, 0.20])                               # distances from the axis [m]
terms = [ch04.cauchy_terms(1000.0, u_sb, tau_sb, gvec, np.array([R, 0, 0]), 0.0) for R in Rs]   # (4.24) at each point
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw=dict(width_ratios=[2.2, 1]))        # the figure and its panels
w = 0.012                                                             # bar width [m]
a.bar(Rs - w, [T_.inertia[0] for T_ in terms], w, color=COLORS["teal"], label="ρDu/Dt")          # mass × acceleration
a.bar(Rs, [T_.body[0] for T_ in terms], w, color=COLORS["muted"], label="ρg (radial: 0)")         # gravity has no radial part
a.bar(Rs + w, [T_.stress_divergence[0] for T_ in terms], w, color=COLORS["orange"], label="∂τ_ij/∂x_i")   # net push
a.plot(Rs, [T_.inertia[0] - T_.body[0] - T_.stress_divergence[0] for T_ in terms], "ko", label="residual")  # draw the data
a.set_xlabel("R [m]"); a.set_ylabel("radial component [N/m³]"); a.legend(fontsize=8)                # axis labels with units
a.set_title("radial: the pressure push supplies the centripetal acceleration", fontsize=10)         # the panel's message
T1 = terms[1]                                                         # the point R = 0.1 m
b.bar(["ρDw/Dt", "ρg_z", "∂τ_i3/∂x_i"], [T1.inertia[2], T1.body[2], T1.stress_divergence[2]],       # draw the data
      color=[COLORS["teal"], COLORS["muted"], COLORS["orange"]])  # one colour per term
b.set_ylabel("vertical component [N/m³]"); b.set_title("vertical at R = 0.1 m", fontsize=10)        # axis labels with units
savefig(fig, "ch04", "cauchy_bars"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Teal and orange bars of equal length growing with R, nothing grey radially, residual dots on zero; vertically a grey and an orange bar of equal size and opposite sign.",
    read="The inward push of the pressure (orange) supplies exactly the centripetal acceleration (teal); vertically the pressure holds the water up against gravity.",
    change="…Ω doubled: both radial bars grow four times (Ω²); the vertical balance is unchanged.")
whatif(r"""
…we knew how τ depends on the motion? Then the 13 unknowns collapse to 5 (ρ, p, u_j) and Cauchy's equation becomes a
closed PDE for the velocity. The simplest possible law — stress linear in the strain rate, the same in every direction —
is the Newtonian fluid (C07).
""")
# =====================================================================================================================
# A.5 §4.5 Constitutive Equation for a Newtonian Fluid — R05, R06, C07
# =====================================================================================================================
nb.section("4.5", "Constitutive Equation for a Newtonian Fluid", intro="""
**What is this section about?** Cauchy's equation needs a rule for the stress. We show first that the stress tensor is
symmetric, then build the simplest law allowed by physics — stress linear in the rate of strain and the same in every
direction — and find it has only two material constants, μ and the bulk viscosity μ_v.
""")
nb.recap("R05", "Only the strain rate can create viscous stress", r"""
Viscous stress must be the same for every observer moving at constant velocity (Galilean invariance, Ch. 3 §3.3), so it
cannot depend on u itself — only on the velocity gradient G. And of $\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}$
*(Eq. 3.11)* only the strain-rate tensor S deforms an element; the rotation part R describes a rigid spin, which deforms
nothing (Ch. 3 §3.4: a rigid motion U + Ω × x has S = 0). So the viscous stress is a function of S alone and vanishes
when S = 0.
""", where="Ch. 3 §3.3–3.4")
nb.code(r"""
G_rigid = kinematics.velocity_gradient_preset("solid_body_rotation", Gamma=2.0, dim=3)   # G of a rigid spin (antisymmetric) [1/s]
print(tensors.strain_rate_tensor(G_rigid))                           # S = 0: no deformation
print(ch04.viscous_stress(G_rigid, mu=1e-3))                         # σ = 0: no viscous stress [Pa]
""", explain="A rigid rotation has no strain rate and so no viscous stress: two 3 × 3 zero matrices.")
nb.recap("R06", "Non-Newtonian fluids exist", r"""
Paint, ketchup, blood, polymer melts and mud do not obey a linear stress–strain-rate law: they shear-thin, have a yield
stress or remember their past (Ch. 1 §1.3 named Bingham and Maxwell materials). This chapter treats the Newtonian fluid
only; Ch. 16 returns to the others.
""", where="Ch. 1 §1.3")
core("C07", "How a fluid decides its stress: the Newtonian law (4.31)", """
Honey between two plates: slide the top plate and it resists; squeeze the plates together and it resists differently;
spin the whole sandwich and it does not resist at all. What rule turns "how the fluid deforms" into "what stress it
carries"?
""", eqs=("4.31",))
problem(r"""
Newton's law $\tau=\mu\,du/dy$ *(1.3)* covers one flow: parallel layers sliding. A real flow shears, stretches and swells
in all directions at once, and Cauchy's equation needs all nine stress components. We want the general rule, as simple
as physics allows, with as few material constants as possible — the one used for water and air in every later chapter,
in CFD codes and in the eddy-viscosity models of turbulence.
""")
idea("""
motion  G = ∂u_i/∂x_j  →  keep only S (R05)  →  linear: σ_ij = K_ijmn S_mn (81 numbers)
                       →  isotropic: K = λδδ + μδδ + γδδ (3 numbers)  →  symmetric: 2 numbers
τ_ij = −p δ_ij  +  2μ S_ij  +  λ S_mm δ_ij                                      (4.31)
       static     shear/stretch  swelling
""")
note("N39", "Constitutive equation", "A material's rule linking stress to deformation; the *Newtonian* fluid has the simplest linear one. Non-Newtonian rules: Ch. 16.")
note("N38", "The stress tensor is symmetric", r"""
so it has six independent components. The book leaves the proof to an exercise; the derivation below writes it out: an
unequal pair τ₁₂ ≠ τ₂₁ would spin a small cube faster and faster as it shrinks. The only exception is a fluid with body
*couples* (torques per unit mass, e.g. polarised molecules in an electric field).
""", equation=EQ["4.25"], ref="4.25")
remind([("exponent rules", r"$h^5/h^3=h^2$, so a torque ∝ h³ over an inertia ∝ h⁵ grows as 1/h² (Ch. 1 P43); torque and moment of inertia were primed in C04 above (P116).")])
D("D08", ref="4.25")
remind([
    ("log–log slope and observed_order", "a power law $y\\propto h^k$ is a straight line of slope k on log–log axes; `observed_order` fits k (Ch. 1 P13)."),
    ("np.logspace", "`np.logspace(-4, -1, 7)` = 7 values evenly spaced in the exponent, 10⁻⁴ … 10⁻¹ (Ch. 1 P06)."),
])
nb.code(r"""
h = np.logspace(-4, -1, 7)                                           # cube sides from 0.1 mm to 10 cm [m]
alpha = ch04.cube_spin_acceleration(1.0, 0.0, 1000.0, h)             # 6(τ12 − τ21)/(ρh²) for a 1 Pa imbalance in water [rad/s²]
print(np.c_[h, alpha])                                                # 60 rad/s² at 1 cm, 6×10⁵ at 0.1 mm
print(round(observed_order(h, alpha), 3))                            # slope −2: it grows as 1/h²
""", explain="""
1. The angular acceleration a 1 Pa stress imbalance would give a water cube of side h.
2. Its log–log slope: it grows as 1/h² — no finite imbalance can survive at a point.
""")
note("N40", "At rest", r"""
a fluid's stress is the same on every plane (isotropic, Ch. 1 §1.6), and the only isotropic second-order tensor is
$\delta_{ij}$ (Ch. 2 §2.5). p is the thermodynamic pressure (e.g. p = ρRT); the minus sign because tension counts
positive. Traction on any plane: −pn (`tensors.traction(ch04.static_stress(1e5), n)` = −10⁵ n for any unit n).
""", equation=EQ["4.26"], ref="4.26")
note("N41", "Moving fluid", r"""
add a viscous part σ_ij that vanishes at rest. p stays the thermodynamic pressure because fluid particles are in local
equilibrium (Ch. 1 §1.8).

> ⚠️ The book calls σ the "deviatoric" stress; it is traceless (a true deviator) only when μ_v = 0 or ∇·u = 0, since its
> trace is 3μ_v∇·u (shown after D10).
""", equation=EQ["4.27"], ref="4.27")
remind([("np.einsum index strings", "`np.einsum('ijmn,mn->ij', K, S)` sums over the repeated letters m, n (Ch. 2 P62).")])
note("N42", "Most general linear law with σ = 0 when S = 0", r"""
K has 3⁴ = 81 components (each of the nine σ's may depend on each of the nine S's); code: `ch04.linear_stress(K, S)` =
`np.einsum('ijmn,mn->ij', K, S)`.
""", equation=EQ["4.28"], ref="4.28")
remind([("random rotations and transform_tensor", "`tensors.random_rotation(rng)` makes a random orthogonal matrix C (via `np.linalg.qr`); `transform_tensor(T, C)` rotates a tensor of any order (Ch. 2 P67).")])
P("P119", "isotropic fourth-order tensor", r"""
A tensor is *isotropic* if its components are the same in every rotated frame. For second order only multiples of
$\delta_{ij}$ qualify. For fourth order the only possibilities are combinations of products of two δ's — three of them:
$\delta_{ij}\delta_{mn}$, $\delta_{im}\delta_{jn}$, $\delta_{in}\delta_{jm}$ (a classical result we cite, not prove). In a
fluid with no preferred direction K must be isotropic, which cuts 81 constants to 3. We check it numerically: rotate K
with 50 random rotations and nothing changes.
""", code="""
K = ch04.isotropic_fourth_order(2.0, 1.0, 0.5)                   # λ, μ, γ
C = tensors.random_rotation(np.random.default_rng(1))            # a random rotation matrix
print(np.abs(tensors.transform_tensor(K, C) - K).max())          # ~1e-15: unchanged
print(tensors.is_isotropic(K, n_rotations=50))                   # (True, largest change ~1e-15)
""")
note("N43", "Isotropy leaves three scalars", r"(⚠️ γ here is a material coefficient — not $C_p/C_v$ and not a shear rate.)", equation=EQ["4.29"], ref="4.29")
note("N44", "Symmetry of σ in i, j", r"""
The subtlety (derivation steps 6–9): on a *symmetric* S the μ and γ terms both give $S_{ij}$, so only the sum μ + γ ever
acts, and σ is automatically symmetric. "γ = μ" is the convention that calls the sum 2μ — which is exactly what makes μ
the viscosity of Newton's law $\tau=\mu\,du/dy$ *(1.3)*. (If one also demands that K itself be symmetric in i, j, γ = μ
is forced — the book's argument.)
""", equation=EQ["4.30"], ref="4.30")
D("D09", ref="4.31", check_src=r"""
lam, mu, gam, th = sp.symbols('lambda mu gamma theta', real=True)             # the three coefficients and a rotation angle
d = lambda i, j: sp.KroneckerDelta(i, j)                                       # δ_ij
K = {(i, j, m, n): lam*d(i, j)*d(m, n) + mu*d(i, m)*d(j, n) + gam*d(i, n)*d(j, m)   # step 4: the isotropic K_ijmn
     for i in range(3) for j in range(3) for m in range(3) for n in range(3)}  # all 81 components
s11, s22, s33, s12, s13, s23 = sp.symbols('S11 S22 S33 S12 S13 S23')          # six independent strain rates
S = sp.Matrix([[s11, s12, s13], [s12, s22, s23], [s13, s23, s33]])             # a general SYMMETRIC S
sig = sp.Matrix(3, 3, lambda i, j: sum(K[i, j, m, n]*S[m, n] for m in range(3) for n in range(3)))   # steps 5–7: σ_ij = K_ijmn S_mn
Smm = S.trace()                                                                 # S_mm = ∇·u
step8 = sig - (lam*Smm*sp.eye(3) + (mu + gam)*S)                               # step 8: σ − [λ S_mm δ + (μ + γ) S]
print("step 8:", sp.simplify(step8) == sp.zeros(3, 3))                         # True: only the sum μ + γ acts
print("σ symmetric:", sp.simplify(sig - sig.T) == sp.zeros(3, 3))              # True for any λ, μ, γ
step10 = sig.subs(gam, mu) - (2*mu*S + lam*Smm*sp.eye(3))                      # steps 9–10: γ = μ gives 2μS + λ S_mm δ
print("steps 9-10:", sp.simplify(step10) == sp.zeros(3, 3))                    # True
Cz = sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])   # rotation about z by θ
rng = np.random.default_rng(4)                                                  # 20 random index sets (keeps the cell fast)
ok = True
for i, j, m, n in rng.integers(0, 3, size=(20, 4)):                             # step 3: K'_ijmn = C_ia C_jb C_mc C_nd K_abcd
    Kp = sum(Cz[i, a]*Cz[j, b]*Cz[m, c]*Cz[n, e]*K[a, b, c, e]
             for a in range(3) for b in range(3) for c in range(3) for e in range(3))  # sum over the four rotated indices
    ok &= sp.simplify(Kp - K[int(i), int(j), int(m), int(n)]) == 0              # the rotated component equals the original
print("isotropic under rotation (20 components):", ok)                         # True
assert ok and sp.simplify(step8) == sp.zeros(3, 3) and sp.simplify(step10) == sp.zeros(3, 3)
""")
P("P120", "deviatoric (traceless) part of a tensor", r"""
Any second-order tensor splits into an isotropic part (its average diagonal times δ) and a traceless remainder:
$A_{ij}=\tfrac13A_{mm}\delta_{ij}+\big(A_{ij}-\tfrac13A_{mm}\delta_{ij}\big)$. For the strain rate the first part is pure
swelling (volume change at rate $S_{mm}=\nabla\cdot\mathbf u$) and the second is pure shape change at constant volume. Its
trace is zero because $\delta_{ii}=3$.
""", code="""
S = np.array([[1.0, 0.5, 0], [0.5, 2.0, 0], [0, 0, 3.0]])     # a strain rate [1/s]
D = ch04.deviatoric_part(S)                                     # S − (tr S/3) δ
print(np.trace(D), np.trace(S)/3)                               # 0.0 and 2.0: the average stretching
""")
D("D10", ref="4.37")
note("N45", "Taking the trace", r"""
of $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ *(4.31)* (set i = j and sum; $\delta_{ii}=3$) gives
$\tau_{ii}=-3p+(2\mu+3\lambda)S_{mm}$, so the line below. ($S_{mm}=\nabla\cdot\mathbf u$ is the volumetric strain rate of
Ch. 3 §3.4, $\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}$ *(3.14)* — the book's pointer
"Section 3.6" should read §3.4.)
""", equation=EQ["4.32"], ref="4.32")
note("N46", "The mean (mechanical) pressure", "is minus the average normal stress:", equation=EQ["4.33"], ref="4.33")
note("N47", "They differ only in an expanding or compressing flow", r"""
by the amount in the line below.

> ⚠️ **Mean vs thermodynamic pressure.** p comes from the equation of state; p̄ is what the normal stresses average to.
> For an incompressible fluid there is no equation of state for p: only a mechanical pressure exists, and only its
> *gradient* matters — adding a constant to p changes nothing in $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$
> *(4.39b)* (demonstrated in C08's code).
""", equation=EQ["4.34"], ref="4.34")
note("N48", "Incompressible Newtonian stress", r"""
($S_{mm}=0$ removes λ). `ch04.newtonian_stress(G, p, mu, incompressible=True)` raises an error if tr G ≠ 0 — a guard
against using (4.35) on a compressible field.
""", equation=EQ["4.35"], ref="4.35")
note("N49", "Bulk viscosity", r"""
$\mu_v=\lambda+\tfrac23\mu$: the resistance to pure expansion. It matters for sound absorption and the inside of shock
waves (Ch. 15) and is nonzero in polyatomic gases (slow exchange of energy with molecular rotation).
""")
note("N50", "Stokes' assumption", r"""
— the code's default (`mu_v=0.0`). Accurate whenever μ_v or the expansion rate is small — nearly always outside acoustics
and shocks.
""", equation=EQ["4.36"], ref="4.36")
note("N51", "The bulk-viscosity form", r"""
(D10's result) separates shape change from volume change. For a parallel flow u = (u(y), 0, 0) it gives
$\tau_{12}=\mu\,du/dy$ — Newton's law $\tau=\mu\,du/dy$ *(1.3)*. Number: shear rate $\dot\gamma=10$ s⁻¹ in water
(μ = 1.0×10⁻³ Pa s) → τ₁₂ = 0.010 Pa, the same as `ch01.newton_shear_stress(1e-3, 10)`.
""", equation=EQ["4.37"], ref="4.37")
nb.worked_example("three motions of water (μ = 1.0×10⁻³ Pa s, p = 0 gauge)", r"""
1. **Shear** $\mathbf u=(\dot\gamma y,0,0)$, $\dot\gamma=10$ s⁻¹: S₁₂ = S₂₁ = 5 s⁻¹, $S_{mm}=0$ → τ₁₂ = 2μS₁₂ =
   2 × 10⁻³ × 5 = 0.010 Pa, all normal stresses 0.
2. **Extension** u = (sx, −sy, 0), s = 1 s⁻¹: S₁₁ = 1, S₂₂ = −1, $S_{mm}=0$ → τ₁₁ + p = +2 mPa (a pull along x),
   τ₂₂ + p = −2 mPa (a push along y); p̄ = p.
3. **Pure expansion** u = (ex, ey, ez), e = 1 s⁻¹: S = diag(1, 1, 1), $S_{mm}=3$ → the shear part
   $S-\tfrac13S_{mm}\delta$ is zero; with Stokes (μ_v = 0) σ = 0 and p̄ = p; with μ_v = 1.0×10⁻³ Pa s, every normal stress
   gains $\mu_vS_{mm}$ = 3 mPa and $p-\bar p=\mu_v\nabla\cdot\mathbf u$ = 3 mPa.
4. **Rigid spin**: S = 0, σ = 0 — no resistance at all.
""")
nb.code(r"""
G_sh = np.array([[0, 10.0, 0], [0, 0, 0], [0, 0, 0]])                # shear, γ̇ = du/dy = 10 1/s (G[i, j] = ∂u_i/∂x_j)
G_ext = np.diag([1.0, -1.0, 0.0])                                    # extension s = 1 1/s
G_exp = np.eye(3)                                                    # pure expansion e = 1 1/s
for name, G in {"shear": G_sh, "extension": G_ext, "expansion": G_exp}.items():  # the three motions of the worked example
    tau = ch04.newtonian_stress(G, p=0.0, mu=1e-3, mu_v=1e-3)          # (4.37) with μ_v = 1e-3 Pa s, p = 0 gauge [Pa]
    print(name, np.round(tau*1e3, 3).tolist(), "mPa;  p − p̄ =", round((0.0 - ch04.mean_pressure(tau))*1e3, 3), "mPa")  # show the numbers
print(ch04.stress_on_plane(G_sh, 0.0, 1e-3, 0.0, np.pi/4))           # (σ_n, τ_s) on the 45° plane: (0.010, 0.0) Pa
print(ch04.bulk_viscosity(-2e-3/3, 1e-3), ch04.lam_from_bulk(0.0, 1e-3))   # Stokes: μ_v = 0 ⇔ λ = −⅔μ
print(ch01.newton_shear_stress(1e-3, 10.0))                          # Ch. 1's τ = μ du/dy: 0.01 Pa
""", explain=r"""
1. The three motions of the worked example through $\tau_{ij}=-p\delta_{ij}+2\mu(S_{ij}-\tfrac13S_{mm}\delta_{ij})+\mu_vS_{mm}\delta_{ij}$
   *(4.37)* with μ_v = 10⁻³ Pa s.
2. Pressure minus mean pressure is nonzero only for expansion: $p-\bar p=\mu_v\nabla\cdot\mathbf u$ = +3 mPa (p̄ = −⅓τ_ii
   = −3 mPa).
3. The traction on a 45° plane in shear flow is a pure pull — that plane is a principal direction.
4. μ_v and λ convert into each other (Stokes: λ = −⅔μ).
5. Ch. 1's Newton law agrees.
""")
nb.check_agree(r"""
rng = np.random.default_rng(7)                                       # a random velocity gradient
G = rng.normal(size=(3, 3)); p, mu, mu_v = 2.0, 1e-3, 5e-4           # pressure [Pa], viscosities [Pa s]
lam = mu_v - 2*mu/3                                                  # λ from the bulk viscosity
S = 0.5*(G + G.T); Smm = np.trace(S)                                 # strain rate and its trace
tau = np.zeros((3, 3))                                               # (4.31) by explicit loops
for i in range(3):  # row i of τ
    for j in range(3):  # column j of τ
        tau[i, j] = -p*(i == j) + 2*mu*S[i, j] + lam*Smm*(i == j)
assert np.allclose(tau, ch04.newtonian_stress(G, p, mu, lam=lam))    # same as the library
assert np.allclose(ch04.linear_stress(ch04.isotropic_fourth_order(lam, mu, mu), S), tau + p*np.eye(3))   # the 81-coefficient route
print(np.round(tau, 5))                                              # the stress [Pa]
""")
remind([
    ("quadratic form n·τ·n; principal directions", r"$\sigma_n=\mathbf n\cdot\boldsymbol\tau\cdot\mathbf n$ is largest and smallest along the eigenvectors of τ, where the shear vanishes (Ch. 2 P82)."),
    ("eigenvalues of S (principal directions)", "the eigenvectors of S are the directions of pure stretching (Ch. 2 P80)."),
])
nb.figure(r"""
th = np.linspace(0, np.pi, 181)                                      # angle of the plane's normal from x [rad]
sh = np.array([ch04.stress_on_plane(G_sh, 0.0, 1e-3, 0.0, t) for t in th])    # shear flow γ̇ = 10 1/s: (σ_n, τ_s) [Pa]
ex = np.array([ch04.stress_on_plane(G_ext, 0.0, 1e-3, 0.0, t) for t in th])   # extension s = 1 1/s
fig, (a, b, c) = plt.subplots(1, 3, figsize=(12, 3.6), gridspec_kw=dict(width_ratios=[1.3, 1.3, 1]))  # the figure and its panels
for ax, dat, lab in ((a, sh, "(a) shear flow γ̇ = 10 s⁻¹"), (b, ex, "(b) extension s = 1 s⁻¹")):  # the two flows, one panel each
    ax.plot(np.rad2deg(th), dat[:, 0]*1e3, color=COLORS["orange"], lw=2, label="normal σ_n")   # normal part [mPa]
    ax.plot(np.rad2deg(th), dat[:, 1]*1e3, color=COLORS["rose"], lw=2, label="shear τ_s")      # signed shear part [mPa]
    ax.axhline(0, color=COLORS["muted"], lw=0.6); ax.set_xlabel("normal angle θ [°]"); ax.set_title(lab, fontsize=10)  # axis labels with units
for v in (45, 135):  # the principal directions of the shear flow
    a.axvline(v, ls="--", color=COLORS["muted"])                     # principal directions of the shear flow
a.set_ylabel("stress [mPa]"); a.legend(fontsize=8)                                        # axis labels with units
c.plot(sh[:, 0]*1e3, sh[:, 1]*1e3, color=COLORS["accent"], label="shear flow")    # Mohr's circle of each flow
c.plot(ex[:, 0]*1e3, ex[:, 1]*1e3, color=COLORS["teal"], label="extension")                         # draw the data
c.set_aspect("equal"); c.set_xlabel("σ_n [mPa]"); c.set_ylabel("τ_s [mPa]"); c.legend(fontsize=8)   # axis labels with units
c.set_title("(σ_n, τ_s): Mohr's circles", fontsize=10)                                        # the panel's message
savefig(fig, "ch04", "stress_on_planes"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Two sinusoids of period 180°, a quarter-period apart; the shear-flow circle is five times the extension-flow circle.",
    read=r"""Shear flow *is* an extension along 45° and a compression along 135°, seen on tilted planes: the shape of the
    curves is the same, only rotated by 45° — the stress follows the principal axes of S. The parametric curve
    (σ_n, τ_s) is a circle (Mohr's circle, a picture of how the traction turns with the plane, Ch. 2 §2.11).""",
    change="…the fluid spun rigidly on top: G changes (it gains an antisymmetric part), the curves do not — only S matters (R05).")
nb.plotly(r"""
th = np.linspace(0, np.pi, 91)                                       # plane normal angle [rad]
def bulk(mu_v):                                                      # pure expansion e = 1 1/s with bulk viscosity μ_v
    sn = np.array([ch04.stress_on_plane(G_exp, 0.0, 1e-3, mu_v, t)[0] for t in th])   # normal stress on each plane, p = 0 [Pa]
    tau = ch04.newtonian_stress(G_exp, p=0.0, mu=1e-3, mu_v=mu_v)     # the full stress tensor
    return {"normal stress + p [mPa]": (np.rad2deg(th), sn*1e3),                                    # the curves for this slider value
            "p − p̄ [mPa]": (np.rad2deg(th), 0*th - ch04.mean_pressure(tau)*1e3)}                   # another named curve (x, y)
fig = slider_figure(bulk, "μ_v", np.linspace(0, 5e-3, 11), unit="Pa s", xlabel="plane normal angle θ [°]",  # precompute every slider position (works on the web page)
                    ylabel="[mPa]", title="Bulk viscosity acts only on expansion", yrange=(-1, 16))  # labels, title and a fixed range
recolor(fig, {"normal stress + p [mPa]": COLORS["orange"], "p − p̄ [mPa]": COLORS["muted"]}, {"p − p̄ [mPa]": "dash"}).show()  # house colours for each trace, then draw
""")
see_read_change(
    see="Two flat lines that rise together as μ_v grows (they overlap).",
    read=r"""With μ_v = 0 the orange line sits at 0 (p̄ = p, Stokes); raising μ_v lifts it uniformly: every plane feels the
    same extra pull $\mu_vS_{mm}=\mu_v\times3\ \text{s}^{-1}$ (in Pa), and $p-\bar p=\mu_v\nabla\cdot\mathbf u$ grows with it.""",
    change="…the flow switched to shear: μ_v would do nothing ($S_{mm}=0$).")
nb.figure(r"""
h = np.logspace(-5, -1, 50)                                          # cube side [m]
fig, ax = plt.subplots(figsize=(6.5, 3.8))                                        # the figure and its panels
for dtau, col in ((0.01, COLORS["teal"]), (0.1, COLORS["blue"]), (1.0, COLORS["rose"])):   # stress imbalances [Pa]
    ax.loglog(h, ch04.cube_spin_acceleration(dtau, 0.0, 1000.0, h), color=col, lw=2, label=f"τ12 − τ21 = {dtau} Pa")  # draw the data
ax.set_xlabel("cube side h [m]"); ax.set_ylabel("angular acceleration [rad/s²]"); ax.legend(fontsize=8)  # axis labels with units
ax.set_title("An unbalanced stress pair would spin a shrinking cube without limit")                 # the panel's message
savefig(fig, "ch04", "cube_spin"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Straight lines of slope −2 climbing to absurd values as h shrinks.",
    read=r"""A stress imbalance at a point would spin the fluid there infinitely fast — impossible, so $\tau_{ij}=\tau_{ji}$ *(4.25)*.""",
    change="…body couples present: a couple (torque) per unit *mass* times the cube's mass ρh³ scales like h³, like (τ₁₂ − τ₂₁)h³, and could balance the imbalance at every size (D08 step 8) — then τ need not be symmetric.")
nb.explainer("newtonian_stress_lab", heading="How does a fluid decide its stress?", why=r"""
**Why interactive rather than a static figure:** the stress law is a chain G → S → τ with nine numbers at each link; turning the plane and the presets yourself shows at once which entries of τ respond to which motion, something no single static figure can hold.

Set a velocity gradient with presets (shear, extension, rotation, expansion) or sliders and watch the chain G → S (and
R) → τ of $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ (4.31) respond together: the element deforms, the
traction on a plane you rotate swings through its normal and shear values, and the matrices light up by source. A second
mode shrinks the spinning cube of D08.
""", tries=[
    "Preset 'pure rotation': G is not zero but every stress entry is — only S matters.",
    "Preset 'expansion' and raise μ_v: p̄ separates from p; switch to 'shear' and μ_v does nothing.",
    "Rotate the plane in shear flow until the shear readout is zero: that is 45°, a principal direction.",
    "Open the Derivation tab at D09 and watch the δ-substitution steps light the matrix entries they produce.",
])
whatif(r"""
…we put (4.37) into Cauchy's equation $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(4.24)*? The
divergence of 2μS becomes μ∇²u plus a term in ∇(∇·u), and we get the Navier–Stokes equation — 4 equations for 5 unknowns
(C08).
""")
# =====================================================================================================================
# A.6 §4.6 Navier–Stokes Momentum Equation — C08
# =====================================================================================================================
nb.section("4.6", "Navier–Stokes Momentum Equation", intro="""
**What is this section about?** Substitute the Newtonian stress into Cauchy's equation: the Navier–Stokes equation. With
constant viscosity and ∇·u = 0 it takes its famous compact form; its viscous force can be written three ways, and with
the viscosity switched off it is Euler's equation.
""")
core("C08", "The Navier–Stokes equation: five accelerations in balance (4.39b)", """
Honey oozing down a pipe and air racing over a wing obey the same equation. Which of its terms matter where?
""", eqs=("4.39b",))
problem("""
Weather models, ocean models, aircraft design codes and blood-flow simulations all solve one equation. Every exact
solution of Ch. 8, every boundary layer of Ch. 9 and every instability of Ch. 11 is a flow in which some of its terms
balance and the others are asleep. Learning to read the terms — which is big, which is zero, and why — is the skill this
block builds.
""")
idea("""
ρ ( ∂u/∂t  +  (u·∇)u )  =  −∇p    +   ρg     +   μ∇²u                        (4.39b)
    local    advective     pressure   gravity    viscous
    blue     teal          orange     grey       rose
Poiseuille: orange ↔ rose     Stokes' first problem: blue ↔ rose     potential flow: teal ↔ orange (rose = 0)
""", r"""
The three example flows are worked below: plane Poiseuille flow (a channel driven by a pressure gradient), Stokes' first
problem (a plate set suddenly in motion under still fluid) and the potential flow round a cylinder (Ch. 3).
""")
D("D11", ref="4.38")
note("N52", "Navier–Stokes with variable viscosity", r"""
(D11's result). **Ledger:** continuity (1) + Navier–Stokes (3) = 4 equations; unknowns ρ, p, u_j = 5. Closed when ρ is
constant or a known function of p alone (a *barotropic* flow, met properly in C11): 5 = 5.
""", equation=EQ["4.38"], ref="4.38")
nb.code(r"""
for s in ("navier_stokes", "barotropic"):                            # two more rows of the §4.1 ledger
    c = ch04.closure_count(s)                                        # call the tested chapter function
    ledger.loc[s] = [c["equations"], c["unknowns"]]
print(ledger)                                                        # navier_stokes 4/5, barotropic 5/5
""")
P("P121", "Schwarz's theorem", r"""
For a smooth function the order of partial derivatives does not matter:
$\frac{\partial}{\partial x}\frac{\partial f}{\partial y}=\frac{\partial}{\partial y}\frac{\partial f}{\partial x}$. Here it
lets us move $\partial/\partial x_j$ outside $\partial u_i/\partial x_i$, turning a term into the gradient of the
divergence ∇(∇·u). It needs continuous second derivatives (true for every viscous flow we meet).
""", code="""
x, y = sp.symbols('x y')                                        # two coordinates
f = sp.sin(x*y) + x**3*sp.exp(y)                                # any smooth function
print(sp.simplify(sp.diff(f, x, y) - sp.diff(f, y, x)))         # 0
""")
D("D12", ref="4.39b")
nb.code(r"""
x, y, z, a, mu, mu_v = sp.symbols('x y z a mu mu_v')                  # coordinates, a toy constant, the two viscosities
X = (x, y, z)
for name, u in {"Poiseuille-like u = (y(1 − y), 0, 0)": [y*(1 - y), 0, 0], "toy u = (a x², 0, 0)": [a*x**2, 0, 0]}.items():
    div = sum(sp.diff(u[m], X[m]) for m in range(3))                 # ∇·u = ∂u_m/∂x_m
    extra = [sp.simplify((mu_v + mu/3)*sp.diff(div, X[j])) for j in range(3)]   # the compressible term of (4.39a)
    print(f"{name}:  ∇·u = {div},  (μ_v + μ/3) ∂(∇·u)/∂x_j = {extra}")   # 0 for the first, (2a(μ_v + μ/3), 0, 0) for the toy
""", explain=r"""
The last term of the constant-viscosity form (4.39a) for two fields: it vanishes for the divergence-free Poiseuille-like
profile, and for the toy field $\mathbf u=(ax^2,0,0)$ (∇·u = 2ax) it is $2a(\mu_v+\tfrac13\mu)$ along x — the D12 check.
""")
note("N53", "The constant-viscosity compressible form", r"""
D12's middle line. Its last term is the only trace of compressibility; ∇·u = 0 removes it and leaves (4.39b),
$\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$.
""", equation=EQ["4.39a"], ref="4.39a")
P("P122", "curl of a curl identity", r"""
For any smooth vector field, $\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u$. It
follows from the ε–δ identity $\varepsilon_{ijk}\varepsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}$ *(2.19)*
of Ch. 2. With ω = ∇×u and ∇·u = 0 it says $\nabla^2\mathbf u=-\nabla\times\boldsymbol\omega$: the Laplacian of an
incompressible velocity is minus the curl of its vorticity.
""", code="""
x, y, z = sp.symbols('x y z'); X = sp.Matrix([x, y, z])          # coordinates
u = sp.Matrix([y**2*z, sp.sin(x)*z, x*y**3])                      # any smooth field
curl = lambda F: sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)])
div = sum(sp.diff(u[i], X[i]) for i in range(3))                  # ∇·u
lap = sp.Matrix([sum(sp.diff(u[k], v, 2) for v in X) for k in range(3)])   # ∇²u component by component
print(sp.simplify(curl(curl(u)) - (sp.Matrix([sp.diff(div, v) for v in X]) - lap)))   # zero vector
""")
gloss(["ε–δ identity (2.19)"])
D("D13", ref="4.40")
note("N54", "The viscous force three ways (incompressible flow)", r"""
**The paradox resolved:** the viscous force contains vorticity although rotation was excluded from the stress law — but
it contains its *derivative*. Solid-body rotation has uniform ω, so ∇×ω = 0 and no viscous force; and then S = 0
everywhere too. Viscosity acts only where vorticity (or strain) varies in space — the vorticity diffusion of Ch. 5.
""", equation=r"(\mu\nabla^2\mathbf u)_j=\mu\frac{\partial^2u_j}{\partial x_i^2}=2\mu\frac{\partial S_{ij}}{\partial x_i}=\mu\frac{\partial}{\partial x_i}\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)=-\mu\varepsilon_{jik}\frac{\partial\omega_k}{\partial x_i}=-\mu(\nabla\times\boldsymbol\omega)_j", ref="4.40")
note("N55", "Euler's equation", r"""
— viscosity negligible (far from walls). Check below: the Ch. 3 cylinder potential flow with its Bernoulli pressure
satisfies it to round-off (and its viscous term is exactly zero anyway, C12).
""", equation=EQ["4.41"], ref="4.41")
P("P123", "complementary error function erfc", r"""
$\operatorname{erfc}(\eta)=1-\operatorname{erf}(\eta)=\frac{2}{\sqrt\pi}\int_\eta^\infty e^{-s^2}ds$ falls smoothly from 1 at
η = 0 to 0 as η grows (0.48 at η = 0.5, 0.16 at 1, 0.005 at 2). It is the shape of a quantity diffusing into a region from
a boundary held at a fixed value: here momentum diffusing up from a plate set suddenly in motion.
`scipy.special.erfc` evaluates it.
""", code="""
from scipy.special import erfc                   # the complementary error function
print(erfc(np.array([0.0, 0.5, 1.0, 2.0])))      # [1. 0.4795 0.1573 0.0047]
""")
nb.worked_example("plane Poiseuille flow of water (gap h = 1 mm, G = −dp/dx = 100 Pa/m)", r"""
Walls at y = 0 and y = h, μ = 1.0×10⁻³ Pa s.

1. Try $u=\frac{G}{2\mu}y(h-y)$, v = w = 0 (steady, parallel).
2. Local term 0 (steady); advective term $u\,\partial u/\partial x=0$ (u does not change along x).
3. Pressure force per volume $-\partial p/\partial x=+G=+100$ N/m³ everywhere.
4. Viscous force $\mu\,\partial^2u/\partial y^2=\mu\times(-G/\mu)=-100$ N/m³ everywhere.
5. Sum 0 ✓ — $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(4.39b)* holds with only two awake
   terms.
6. Peak speed at y = h/2: $u_{max}=Gh^2/(8\mu)=100\times10^{-6}/(8\times10^{-3})=0.0125$ m/s = 12.5 mm/s.
""")
nb.code(r"""
u_po, p_po = ch04.exact_field("poiseuille", G=100.0, h=1e-3, mu=1e-3)   # the Poiseuille velocity and pressure fields
for y in (0.0, 2.5e-4, 5e-4):                                        # three heights in the 1 mm gap [m]
    T = ch04.ns_incompressible_terms(u_po, p_po, np.array([0.0, y, 0.0]), 0.0, rho=1000.0, mu=1e-3,  # call the tested chapter function
                                     g=np.zeros(3), h=1e-5, per="volume")   # the five terms of (4.39b) [N/m³]
    print(y, "local", T.local[0], "adv", T.advective[0], "press", T.pressure[0], "visc", round(T.viscous[0], 6), "res", f"{T.residual[0]:.1e}")  # show the numbers
print(ch04.ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0, per="volume", G=100.0, h=1e-3, mu=1e-3))   # the one-number wrapper the navier_stokes_term_balance explainer mirrors
p_shift = lambda x, t: p_po(x, t) + 1e5                              # the same pressure plus a constant 1e5 Pa
print(ch04.ns_incompressible_terms(u_po, p_shift, np.array([0, 2.5e-4, 0]), 0.0, 1000.0, 1e-3, g=np.zeros(3),  # show the numbers
                                   per="volume").residual)            # still ≈ 0: only ∇p enters (N47)
X_, y_, Z_, T_ = sp.symbols('x y z t'); mu0, b = sp.symbols('mu0 b', positive=True)   # symbols for a variable viscosity
print(ch04.navier_stokes_sym(1000, [y_*(1 - y_), 0, 0], -2*X_, (X_, y_, Z_), T_, mu=mu0*(1 + b*y_))[0])   # x-residual of (4.38)
print(ch04.ns_terms_preset("cylinder", 1.5, 0.7, component=0, U=1.0, a=1.0, rho=1.0))   # an Euler flow: viscous term 0
""", explain=r"""
1. The Poiseuille field.
2. The five terms at three heights: pressure +100 and viscous −100 N/m³ everywhere, the rest zero.
3. The parity wrapper used by the `navier_stokes_term_balance` explainer.
4. A constant pressure offset changes nothing (only ∇p enters).
5. With μ = μ₀(1 + by) the viscous term keeps a dμ/dy part: the residual $-b\mu_0(1-2y)+2\mu_0(1+by)-2$ vanishes only for
   b = 0 and μ₀ = 1 — the variable-viscosity form (4.38), not (4.39b), is needed when μ varies.
6. An Euler flow (the Ch. 3 cylinder with its Bernoulli pressure): the advective acceleration (−0.0321 m/s²) is balanced
   by the pressure term, and the viscous force is exactly zero.
""")
nb.check_agree(r"""
dy, y0 = 1e-6, 2.5e-4                                                # step and height [m]
u = lambda y: 100/(2e-3)*y*(1e-3 - y)                                # Poiseuille profile by hand [m/s]
visc = 1e-3*(u(y0 + dy) - 2*u(y0) + u(y0 - dy)) / dy**2              # μ ∂²u/∂y² by a second difference [N/m³]
press = 100.0                                                        # −∂p/∂x = G [N/m³]
print(press, visc)                                                   # 100 and −100
assert abs(press + visc) < 1e-3                                      # (4.39b): they cancel
assert np.isclose(visc, ch04.ns_terms_preset("poiseuille", 0.0, y0, per="volume", G=100.0, h=1e-3, mu=1e-3)["viscous"], rtol=1e-5)
""")
nb.md("**What does the code above do?** It writes the Poiseuille profile by hand, takes μ ∂²u/∂y² by a second difference (−100 N/m³), checks it cancels the pressure push +100 N/m³ and equals the library's viscous term.")
nb.figure(r"""
yy = np.linspace(0, 1e-3, 101)                                       # across the gap [m]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw=dict(width_ratios=[1, 1.5]))        # the figure and its panels
a.plot(ch04.plane_poiseuille(yy, G=100.0, h=1e-3, mu=1e-3)*1e3, yy*1e3, color=COLORS["teal"], lw=2.5)   # u(y) [mm/s]
a.set_xlabel("u [mm/s]"); a.set_ylabel("y [mm]"); a.set_title("(a) the Poiseuille parabola", fontsize=10)  # axis labels with units
ys = np.linspace(1e-4, 9e-4, 5)                                      # five heights [m]
res = [ch04.ns_terms_preset("poiseuille", 0.0, y, per="volume", G=100.0, h=1e-3, mu=1e-3) for y in ys]  # call the tested chapter function at each point
w = 0.03                                                              # bar width [mm]
for k, (name, col) in enumerate((("local", COLORS["blue"]), ("advective", COLORS["teal"]), ("pressure", COLORS["orange"]), ("viscous", COLORS["rose"]))):  # one bar per term, side by side
    b.bar(ys*1e3 + (k - 1.5)*w, [r_[name] for r_ in res], w, color=col, label=name)   # each term at each height
b.plot(ys*1e3, [r_["residual"] for r_ in res], "ko", label="residual")                              # draw the data
b.axhline(0, color=COLORS["muted"], lw=0.6); b.set_xlabel("y [mm]"); b.set_ylabel("[N/m³]"); b.legend(fontsize=8, ncol=3)  # axis labels with units
b.set_title("(b) at every height: pressure push = viscous drag", fontsize=10)                       # the panel's message
savefig(fig, "ch04", "poiseuille_terms"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="A parabola; at every height the same two equal and opposite bars (orange up, rose down) and zero-height blue and teal bars.",
    read="In fully developed pipe/channel flow nothing accelerates: the pressure push is used up by friction at every point.",
    change="…G doubled: both bars double and the parabola doubles (u ∝ G) — the balance is linear.")
nb.plotly(r"""
yy = np.linspace(0, 10e-3, 150)                                      # height above the plate [m]
nu = 1e-6                                                            # water [m²/s]
def stokes(t):                                                       # Stokes' first problem at time t [s]
    u = ch04.stokes_first_problem(yy, t, U=1.0, nu=nu)               # u(y, t)/U = erfc(y/2√(νt))
    dt = 1e-3*t                                                      # a small time step for ∂u/∂t
    local = (ch04.stokes_first_problem(yy, t + dt, 1.0, nu) - ch04.stokes_first_problem(yy, t - dt, 1.0, nu)) / (2*dt)
    visc = nu*np.gradient(np.gradient(u, yy), yy)                   # ν ∂²u/∂y² on the grid
    scale = 1/np.max(np.abs(local))                                  # show both on the profile's scale
    return {"u(y, t)/U": (yy*1e3, u), "local ∂u/∂t (scaled)": (yy*1e3, local*scale),                # the curves for this slider value
            "viscous ν∂²u/∂y² (scaled)": (yy*1e3, visc*scale)}                                      # another named curve (x, y)
fig = slider_figure(stokes, "t", np.logspace(-2, 1, 20 if not FAST else 10), unit="s", xlabel="y [mm]", ylabel="[–]",  # precompute every slider position (works on the web page)
                    title="Momentum diffuses from the wall: the local term is paid by viscosity", yrange=(-0.05, 1.05))  # labels, title and a fixed range
recolor(fig, {"u(y, t)/U": COLORS["teal"], "local ∂u/∂t (scaled)": COLORS["blue"],                  # house colours for each trace, then draw
              "viscous ν∂²u/∂y² (scaled)": COLORS["rose"]}, {"viscous ν∂²u/∂y² (scaled)": "dash"}).show()  # colour of each named trace
""")
see_read_change(
    see="The teal profile thickens as t grows; the blue curve and the dashed rose curve lie on top of each other at every t.",
    read=r"""They balance: $\partial u/\partial t=\nu\,\partial^2u/\partial y^2$ is (4.39b) for this flow. The profile's
    thickness grows like $2\sqrt{\nu t}$: 2 mm after 1 s, 6.3 mm after 10 s.""",
    change="…oil (ν 100 × larger): the same shapes appear 100 × sooner — only the combination y/√(νt) matters (C15).")
nb.figure(r"""
n = 21 if not FAST else 15                                           # arrows per side
xs = np.linspace(-2, 2, n); Xg, Yg = np.meshgrid(xs, xs)             # a grid in the plane [m]
P = np.array([Xg.ravel(), Yg.ravel()])                               # points as a (2, N) array
Pz = np.vstack([P, 0*P[0]])                                          # the same points in 3-D (z = 0)
u_sb2, _ = ch04.exact_field("solid_body", Omega=1.0)                 # rigid rotation Ω = 1 rad/s
u_lo, _ = ch04.exact_field("lamb_oseen", Gamma=2*np.pi, nu=0.01, t0=25.0)   # Lamb–Oseen vortex, σ² = 4ν(t + t0) = 1 m²
fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))                                        # the figure and its panels
for ax, uf, pts, lab in ((axs[0], u_sb2, Pz, "solid-body rotation"), (axs[1], u_lo, P, "Lamb–Oseen vortex")):  # the two vortices, one panel each
    w = kinematics.vorticity(uf, pts, 0.0)[2].reshape(Xg.shape)      # ω_z on the grid [1/s]
    lap, f2, f3 = ch04.viscous_force_forms(uf, pts, 0.0, mu=0.01)    # μ∇²u, 2μ ∂S_ij/∂x_i, −μ∇×ω (μ = ρν, ρ = 1) [N/m³]
    assert np.allclose(lap, f2, atol=1e-6) and np.allclose(lap, f3, atol=1e-6)   # D13: the three forms of (4.40) agree
    ax.pcolormesh(Xg, Yg, w, cmap="Purples", vmin=0, vmax=2.2, shading="auto")   # vorticity heat map
    ax.quiver(Xg, Yg, lap[0].reshape(Xg.shape), lap[1].reshape(Xg.shape), color=COLORS["rose"], scale=0.25)   # viscous force
    ax.set_aspect("equal"); ax.set_title(f"{lab}: max |μ∇²u| = {np.abs(lap).max():.1e}", fontsize=10); ax.set_xlabel("x [m]")  # axis labels with units
axs[0].set_ylabel("y [m]")                                        # axis labels with units
savefig(fig, "ch04", "viscous_paradox"); plt.show()                                        # save the PNG to outputs/ch04, then draw
r_ = np.array([0.2, 0.5, 1.0, 1.5, 2.0])                               # radii along the x-axis [m]
f_th = [ch04.viscous_force_forms(u_lo, np.array([r, 0.0]), 0.0, mu=0.01)[0][1] for r in r_]   # θ-component of μ∇²u at (r, 0)
print("viscous force along e_θ [N/m³]:", np.round(f_th, 5))            # all negative: against the counter-clockwise swirl
""", see="A uniform vorticity disc with no arrows; a peaked vorticity blob whose arrows all run clockwise, against the counter-clockwise swirl, largest near r ≈ σ = 1 m (the printed θ-components are all negative).",
    read="""Every arrow points against the swirl: viscosity slows the vortex at every radius, most strongly near r ≈ σ,
    while its vorticity spreads outward (ω rises outside the core as σ² = 4ν(t + t₀) grows — Ch. 5). A rigid rotation,
    whose vorticity is uniform, feels nothing: $\\mu\\nabla^2\\mathbf u=-\\mu\\nabla\\times\\boldsymbol\\omega$ = 0 there.""",
    change="…ν doubled: the arrows double and the core spreads twice as fast (σ² = 4νt).")
nb.explainer("navier_stokes_term_balance", heading="Which terms of Navier–Stokes are awake here?", why=r"""
**Why interactive rather than a static figure:** each exact solution balances different terms at different points and times; clicking points and scrubbing time lets you find the regions where each term wakes up, instead of reading one pre-chosen probe.

Pick an exact solution — Couette, Poiseuille, Stokes' first problem, Taylor–Green (a decaying array of vortices),
Lamb–Oseen or the potential flow round a cylinder — click any point and read the five terms of
$\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ (4.39b) as bars that add to zero; scrub time for the
unsteady ones and watch the balance move. Its Explain tab also computes a local Reynolds number, speed × length/ν
(C15 treats it properly).
""", tries=[
    "Poiseuille: click anywhere — only orange and rose, equal and opposite.",
    "Stokes' first problem: press ▶ and watch the blue (local) bar trade with the rose (viscous) bar.",
    "Cylinder (Euler): the viscous bar is exactly zero although μ ≠ 0 — why? (Explain, section 5).",
    "Open the third view: the viscous force computed three ways of (4.40), $\\mu\\nabla^2\\mathbf u=2\\mu\\partial S_{ij}/\\partial x_i=-\\mu\\nabla\\times\\boldsymbol\\omega$, gives the same arrow.",
])
whatif(r"""
…the observer were on a turntable, or on the rotating Earth? Newton's law holds only in an inertial frame; in a rotating
one extra "apparent" forces appear — Coriolis and centrifugal (C09). For the atmosphere and ocean they are not small.
""")
# =====================================================================================================================
# A.7 §4.7 Noninertial Frame of Reference — C09
# =====================================================================================================================
nb.section("4.7", "Noninertial Frame of Reference", intro="""
**What is this section about?** The equations so far hold for an observer who is not accelerating. We rewrite
Navier–Stokes for an observer who accelerates and rotates — the frame every geophysical flow is described in — and meet
the Coriolis and centrifugal terms.
""")
core("C09", "Seen from a rotating frame: Coriolis and centrifugal forces (4.45)", """
Roll a ball straight across a spinning merry-go-round. Someone standing on the ground sees a straight line; the rider
sees it curve. Who is right — and what force does the rider need to invent?
""", eqs=("4.45",))
problem("""
We live on a rotating planet and measure winds and currents relative to the ground. For a thrown ball the rotation hardly
matters; for a hurricane or an ocean gyre that lasts days it dominates. Turbomachines and centrifuges are also analysed in
frames that turn with them. We need Navier–Stokes as seen by an observer who translates with acceleration dU/dt and
rotates at Ω(t).
""")
idea("""
inertial observer:  a = F/m                       rotating observer: a' = F/m + (apparent forces)
the rotating axes e'_i turn:  de'_i/dt = Ω × e'_i   →  velocity gains Ω × x'      (4.42)
differentiate again           →  a = dU/dt + a' + 2Ω×u' + dΩ/dt×x' + Ω×(Ω×x')      (4.43)
                                    frame   rel.  Coriolis angular  centripetal
move them to the force side   →  −dU/dt, −2Ω×u' (amber), −dΩ/dt×x', −Ω×(Ω×x') (purple)   in (4.45)
""")
note("N56", "An inertial frame", r"""
is one in which Newton's law holds as written — fixed to the distant stars, or moving at constant velocity relative to
them. A laboratory on the Earth is nearly inertial over seconds and metres; over hours and hundreds of kilometres it is
not (Ch. 13).
""")
nb.md(r"""
> ⚠️ **The prime changes meaning again.** In this section $\mathbf x'$, $\mathbf u'$, $\mathbf a'$ and $D'/Dt$ are
> measured in the **rotating** frame O′1′2′3′ (Ch. 3's primes were a translating frame; §4.9's p′, ρ′ will be
> perturbations).
""")
remind([
    ("frames of reference and relative velocity", r"velocities add between frames: $\mathbf u=\mathbf U+\mathbf u'$ for a frame moving at $\mathbf U$ (Ch. 3 P96)."),
    ("rotating frame of reference (first look)", r"a rotating observer measures the vorticity $\omega'=\omega-2\Omega$ (Ch. 3 P103)."),
    ("product rule", r"$d(fg)=f\,dg+g\,df$ (Ch. 1 P38)."),
])
P("P124", "derivative of a rotating unit vector", r"""
A unit vector e′ fixed to a frame turning at angular velocity Ω keeps its length but changes direction: its tip moves on a
circle round the Ω axis. In time dt it moves a distance (sin α)|Ω| dt (α the angle between e′ and Ω), perpendicular to
both — exactly the cross product: $d\mathbf e'/dt=\boldsymbol\Omega\times\mathbf e'$. Here it is why a rotating
observer's axes add terms to every time derivative. (`ch04.rotating_basis` turns the three axes with Rodrigues' rotation
formula, the rotation of Ch. 2 §2.2 about an arbitrary axis.)
""", code="""
Om = np.array([0.0, 0.0, 0.5])                                # 0.5 rad/s about z
E = ch04.rotating_basis(Om, 1.0); dE = ch04.basis_rate(Om, 1.0)   # rows e'_i(t) and d e'_i/dt (central difference)
print(np.allclose(dE, np.cross(Om, E)))                       # True: de'/dt = Ω × e'
""")
gloss(["Rodrigues rotation of a basis"])
P("P125", "product rule for a cross product", r"""
$\frac{d}{dt}(\mathbf a\times\mathbf b)=\dot{\mathbf a}\times\mathbf b+\mathbf a\times\dot{\mathbf b}$ — like the ordinary
product rule, but the order of the factors must be kept (a × b = −b × a).
""", code="""
t = sp.symbols('t')                                              # time
a = sp.Matrix([sp.cos(t), t, 1]); b = sp.Matrix([t**2, 0, sp.sin(t)])   # two moving vectors
print(sp.simplify(a.cross(b).diff(t) - (a.diff(t).cross(b) + a.cross(b.diff(t)))))   # zero vector
""")
D("D14", ref="4.42")
note("N57", "Velocity seen from the inertial frame", r"(D14's result). u′ is how fast the *components* x′_i change — what the rotating observer measures.",
     equation=r"\mathbf u=\mathbf U+\mathbf u'+\boldsymbol\Omega\times\mathbf x',\qquad\frac{d\mathbf e'_i}{dt}=\boldsymbol\Omega\times\mathbf e'_i", ref="4.42")
D("D15", ref="4.43", check_src=r"""
t = sp.symbols('t', real=True)                                     # time
th = sp.Function('theta')(t)                                       # rotation angle about z: Ω = θ̇ e_z, dΩ/dt = θ̈ e_z
R = sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])   # columns = the rotating axes e'_i(t)
X = sp.Matrix([sp.Function(f'X{i}')(t) for i in (1, 2, 3)])        # position of the moving origin O' (inertial components)
xp = sp.Matrix([sp.Function(f'xp{i}')(t) for i in (1, 2, 3)])      # the particle's components x'_i in the rotating frame
x = X + R*xp                                                        # D14 step 1–2: x = X + x'_i e'_i
a_inertial = x.diff(t, 2)                                          # the true acceleration, differentiated twice
a_in_primed = sp.simplify(R.T*a_inertial)                          # the same vector written in rotating components
Om = sp.Matrix([0, 0, th.diff(t)]); Omdot = Om.diff(t)             # Ω and dΩ/dt (the axis z is shared by both frames)
up, ap = xp.diff(t), xp.diff(t, 2)                                 # u' = dx'_i/dt, a' = d²x'_i/dt² (what the rider measures)
rhs = R.T*X.diff(t, 2) + ap + 2*Om.cross(up) + Omdot.cross(xp) + Om.cross(Om.cross(xp))   # the five terms of (4.43)
print([sp.simplify(e) for e in (a_in_primed - rhs)])               # [0, 0, 0]: (4.43) holds exactly
one = rhs - Om.cross(up)                                           # the same with 1 × Ω × u' instead of 2 ×
print([sp.simplify(e) for e in (a_in_primed - one)])               # = Ω × u' ≠ 0: the factor 2 is necessary
assert all(sp.simplify(e) == 0 for e in (a_in_primed - rhs))
""")
note("N58", "Acceleration in the inertial frame, term by term", r"""
Names: frame acceleration, acceleration seen in the rotating frame, **Coriolis** term, angular-acceleration term,
**centripetal** term. The factor 2: one Ω × u′ comes from the turning axes (differentiating u′'s basis), one from the
moving position (differentiating Ω × x′).
""", equation=EQ["4.43"], ref="4.43")
note("N59", "For a fluid particle", r"a = Du/Dt in the inertial frame and a′ = D′u′/Dt in the rotating one:",
     equation=r"\Big(\frac{D\mathbf u}{Dt}\Big)_{O123}=\Big(\frac{D'\mathbf u'}{Dt}\Big)_{O'1'2'3'}+\frac{d\mathbf U}{dt}+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')", ref="4.44")
D("D16", ref="4.45")
nb.md(r"""
> ⚠️ **Coriolis *term* vs Coriolis *force*.** On the acceleration side,
> $\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\dots$ *(4.43)* has **+2Ω × u′**; moved to
> the force side, (4.45) has **−2Ω × u′** per unit mass. The book calls both "Coriolis acceleration". We say *Coriolis
> term* for the first and *Coriolis force per unit mass* for the second; the second is the one that deflects a moving
> parcel to the **right** in the Northern Hemisphere (Ω up). Likewise Ω × (Ω × x′) is *centripetal* (points to the axis)
> and −Ω × (Ω × x′) *centrifugal* (away from it). In code: `frame_acceleration_terms` returns the + terms,
> `apparent_body_forces` the − terms.
""")
note("N60", "Frame acceleration −dU/dt", r"""
the push back into your seat when a car accelerates. An aircraft on a parabolic arc has dU/dt = g, so g − dU/dt = 0:
weightlessness (the cell prints the total apparent body force: zero).
""")
nb.code(r"""
print(ch04.apparent_body_forces(np.zeros(3), np.zeros(3), np.zeros(3), dU_dt=np.array([0, 0, -9.81]))["total"])   # [0 0 0]: weightless
""")
P("P126", "latitude, Earth's rotation rate and the local vertical", r"""
The Earth turns once per sidereal day (86 164 s), Ω = 2π/86164 = 7.292×10⁻⁵ rad/s, about the axis through the poles. At
latitude φ the local vertical makes an angle with that axis, so Ω splits into a vertical part Ω sin φ (full at the
poles, zero at the equator) and a horizontal part Ω cos φ along the meridian. For horizontal winds the vertical part does
the steering — twice it, f = 2Ω sin φ, is the *Coriolis parameter* of Ch. 13 (named here only).
""", code="""
Om = 2*np.pi/86164.0                                           # Earth's rotation rate [rad/s]
for lat in (0, 30, 45, 90):                                    # latitude [°]
    print(lat, round(Om*np.sin(np.deg2rad(lat)), 9), round(ch04.coriolis_parameter(np.deg2rad(lat)), 9))   # Ω sin φ, f = 2Ω sin φ [1/s]
""")
gloss(["Coriolis parameter f = 2Ω sin φ (named)"])
remind([("small-angle approximation", r"$\sin x\approx x$ for |x| ≪ 1 rad (Ch. 3 P100).")])
D("D17")
note("N61", "The Coriolis force per unit mass −2Ω × u′", r"""
depends on the velocity, not the position; it is perpendicular to u′, so it changes the direction but never the speed (it
does no work). A projectile fired horizontally from the North Pole at speed u is deflected sideways by Ωut² after time t
(D17), an angle Ωt — exactly the Earth's rotation in that time: seen from space the path is straight and the ground
turned under it. Number: u = 10 m/s, t = 1 h: Coriolis acceleration 2Ωu = 1.46×10⁻³ m/s², forward 36 km, deflection
Ωut² = 9.45 km (the exact offset ut sin Ωt = 9.34 km), angle Ωt = 15.0°.
""")
note("N62", "Highs and lows", r"""
In cylindrical coordinates about a pressure centre, with $\Omega_z>0$ (NH), flow leaving a high ($u_R>0$) feels
$-2\boldsymbol\Omega\times\mathbf u=-2\Omega_zu_R\,\mathbf e_\varphi$: clockwise. Flow into a low ($u_R<0$) is turned
counter-clockwise. In the Southern Hemisphere $\Omega_z<0$ and both reverse. (The steady wind round a high or low —
geostrophic and gradient wind — is Ch. 13.)
""", equation=r"-2\boldsymbol\Omega\times\mathbf u=-2\Omega_zu_R\,\mathbf e_\varphi")
note("N63", "Angular-acceleration term −(dΩ/dt) × x′", "only when the rotation rate or axis changes (a spinning-up centrifuge); zero for the Earth. One bar in the explainer below.")
D("D18")
note("N64", "Centrifugal term and effective gravity", r"""
For steady rotation about z, $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')=\Omega^2R\,\mathbf e_R$ —
outward, growing with the distance R from the axis — and it is minus the gradient of $-\tfrac12\Omega^2R^2$ (D18). It
adds to gravitation: $\mathbf g_e=\mathbf g_n+\Omega^2R\,\mathbf e_R$ is the "gravity" a plumb line shows, and the sea
surface is an equipotential of $\Phi_e=gz-\tfrac12\Omega^2R^2$ — the Earth bulges at the equator by 2(a − b) = 42.77 km
(WGS-84; the book's "42 km" is rounded; its Fig. 4.9 caption says "budge" for bulge). Number: at the equator
Ω²a = (7.292×10⁻⁵)² × 6.378×10⁶ = 0.0339 m/s², 0.35 % of g. Meteorologists fold it into g and the geopotential (Ch. 13).
""")
nb.code(r"""
print(round(ch04.earth_oblateness_diameter()/1e3, 2), "km")          # 2(a − b) for the WGS-84 ellipsoid: 42.77 km
""")
note("N65", "Ex. 4.5 — a von Kármán viscous pump", r"""
a disc rotating under fluid, written in the frame that turns with the disc, in cylindrical components (the Appendix-B
operators: curvilinear coordinates add terms such as $-u_\varphi^2/R$ because the unit vectors turn). The rotation terms
are $\rho[2\Omega_zu_\varphi+\Omega_z^2R]$ in the radial equation and $\rho[-2\Omega_zu_R]$ in the azimuthal one; the
axial equation has none. The cell prints the three equations from `ch04.rotating_pump_equations()` and the rotation terms
alone. Rotating flows between discs and cylinders: Ch. 8.
""")
nb.code(r"""
eqs = ch04.rotating_pump_equations()                                # steady, axisymmetric NS in the disc's frame (sympy)
for name, e in zip(("R", "φ", "z"), eqs):                            # radial, azimuthal and axial momentum
    print(f"{name}:", e)                                        # show the numbers
print("rotation terms:", ch04.rotating_pump_terms()["rotation_terms"])   # [ρ(Ω_z²R + 2Ω_z u_φ), −2ρΩ_z u_R, 0]
""")
nb.worked_example("a projectile at the North Pole", r"""
Ω = 7.292×10⁻⁵ rad/s (up), u = 10 m/s along +x, no friction.

1. Coriolis force per mass $-2\boldsymbol\Omega\times\mathbf u'$: Ω = Ωe_z, u′ = u e_x, e_z × e_x = e_y, so it is
   −2Ωu e_y = −1.458×10⁻³ e_y m/s² — to the right of the motion.
2. In t = 3600 s the ball goes ut = 36 000 m forward.
3. Sideways, with a nearly constant acceleration: ½(2Ωu)t² = Ωut² = 7.292×10⁻⁵ × 10 × 3600² = 9451 m.
4. Angle 9451/36000 = 0.2625 rad = 15.0° = Ωt.
5. Exact (straight inertial line seen from the turning frame): ut sin Ωt = 9342 m — the small-angle estimate is 1.2 % high
   after an hour.
6. Centrifugal at the equator for scale: Ω²a = 0.0339 m/s².
""")
nb.code(r"""
Om = np.array([0, 0, ch04.OMEGA_EARTH])                              # Earth's rotation at the pole [rad/s]
terms = ch04.frame_acceleration_terms(np.zeros(3), np.array([10.0, 0, 0]), np.array([1000.0, 0, 0]), Om)   # the + terms of (4.43)
print({k: np.round(v, 8).tolist() for k, v in terms.items()})        # coriolis (0, +1.458e-3, 0), centripetal (−5.3e-6, 0, 0)
forces = ch04.apparent_body_forces(np.array([10.0, 0, 0]), np.array([1000.0, 0, 0]), Om)   # the − terms of (4.45)
print([f"{v:+.3e}" for v in forces["coriolis"]], [f"{v:+.3e}" for v in forces["centrifugal"]])   # (0, −1.458e-3, 0) and (+5.3e-6, 0, 0) [m/s²]
pr = ch04.coriolis_projectile(10.0, ch04.OMEGA_EARTH, np.array([600.0, 1800.0, 3600.0]))   # after 10, 30, 60 min
print(pr["forward"], np.round(pr["deflection_small"], 1), np.round(pr["deflection"], 1), np.round(np.rad2deg(pr["angle"]), 2))  # show the numbers
print(ch04.centrifugal_acceleration(Om, np.array([ch04.EARTH_RADIUS, 0, 0])), ch04.effective_gravity(np.deg2rad(45.0)))  # show the numbers
""", explain=r"""
1. The five acceleration terms of $\mathbf a=\frac{d\mathbf U}{dt}+\mathbf a'+2\boldsymbol\Omega\times\mathbf u'+\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$
   *(4.43)* for a particle moving at 10 m/s, 1 km from the axis.
2. The same as forces per mass in (4.45) — every sign flipped.
3. The projectile after 10, 30, 60 min: forward distance, the small-angle Ωut² and the exact offset [m], and the angle [°].
4. The centrifugal acceleration at the equator (0.0339 m/s²) and the effective gravity at 45° (9.783 m/s², tilted
   1.73×10⁻³ rad from the radial direction, with g_n = 9.80 m/s² — a toy spherical Earth with uniform gravitation; the real standard gravity 9.80665 m/s² at 45° also includes the Earth's flattening).
""")
nb.check_agree(r"""
u_rot, x_rot = np.array([10.0, 0, 0]), np.array([1000.0, 0, 0])      # velocity and position seen by the rider
cor = 2*np.cross(Om, u_rot)                                          # +2Ω × u' by hand
cen = np.cross(Om, np.cross(Om, x_rot))                              # Ω × (Ω × x') by hand
assert np.allclose(cor, terms["coriolis"]) and np.allclose(cen, terms["centripetal"])   # same as the library
Omz = 1e-3                                                           # a faster rotation so the curvature is large [rad/s]
tt = np.linspace(0, 600, 601); dt = tt[1] - tt[0]                    # 10 minutes in 1 s steps
x_in = np.array([10.0*tt, 0*tt, 0*tt])                               # a straight inertial path, a = 0 [m]
x_r = np.array([ch04.rotating_basis(np.array([0, 0, Omz]), t_) @ x_in[:, k] for k, t_ in enumerate(tt)]).T   # the same path in the turning frame
u_r = np.gradient(x_r, dt, axis=1); a_r = np.gradient(u_r, dt, axis=1)   # velocity and acceleration the rider measures
k = 300                                                              # a point in the middle
pred = -2*np.cross([0, 0, Omz], u_r[:, k]) - np.cross([0, 0, Omz], np.cross([0, 0, Omz], x_r[:, k]))   # (4.43) with a = 0
print(a_r[:, k], pred)                                               # the measured and predicted a'
assert np.allclose(a_r[:, k], pred, rtol=1e-4, atol=1e-8)            # the rider's acceleration is all apparent force
""")
nb.animation(r"""
u0, boost = 10.0, 10.0                                               # speed [m/s]; Ω exaggerated ×10 so the curve shows
Omx = boost*ch04.OMEGA_EARTH                                         # the exaggerated rotation rate [rad/s]
frames = 60 if not FAST else 30                                        # number of frames (halved when FAST)
tt = np.linspace(0, 3600.0, frames)                                  # one hour [s]
inert, rot = ch04.projectile_paths(u0, Omx, tt)                      # (2, n) paths in the two frames [m]
L = u0*tt[-1]/1e3                                                    # disc radius [km]
fig, (a, b) = plt.subplots(1, 2, figsize=(8.5, 4.2), dpi=80)                                        # the figure and its panels
for ax, ttl in ((a, "inertial view: straight path, turning floor"), (b, "rotating view: the path bends right")):  # same disc in both panels
    ax.add_patch(plt.Circle((0, 0), L, color=COLORS["grid"], alpha=0.5)); ax.set_aspect("equal")
    ax.set_xlim(-1.1*L, 1.1*L); ax.set_ylim(-1.1*L, 1.1*L); ax.set_title(ttl, fontsize=9); ax.set_xlabel("km")  # axis labels with units
(mer,) = a.plot([0, L], [0, 0], color=COLORS["muted"], lw=2)          # a meridian painted on the turning floor
(pa,) = a.plot([], [], color=COLORS["rose"], lw=2.5)                 # inertial path
(pb,) = b.plot([], [], color=COLORS["rose"], lw=2.5)                 # rotating-frame path
b.plot(u0*tt/1e3, -Omx*u0*tt**2/1e3, "--", color=COLORS["muted"], label="Ωut² (small angle)")   # ghost parabola
b.plot([0, L], [0, 0], color=COLORS["muted"], lw=2); b.legend(fontsize=7, loc="lower left")         # legend
arrow = [None]                                                       # the Coriolis arrow, redrawn each frame


def update(i):                                                       # draw frame i
    th = Omx*tt[i]                                                   # angle turned so far [rad]
    mer.set_data([0, L*np.cos(th)], [0, L*np.sin(th)])               # the floor turns counter-clockwise
    pa.set_data(inert[0, :i+1]/1e3, inert[1, :i+1]/1e3)                                        # move this artist to the new frame
    pb.set_data(rot[0, :i+1]/1e3, rot[1, :i+1]/1e3)                                        # move this artist to the new frame
    if arrow[0] is not None:  # remove last frame's Coriolis arrow
        arrow[0].remove()
    j = max(i, 1)
    vel = (rot[:, j] - rot[:, j-1]) / (tt[j] - tt[j-1])              # the rider's velocity u' [m/s]
    F = ch04.coriolis_acceleration(np.array([0, 0, Omx]), np.array([vel[0], vel[1], 0.0]))   # −2Ω × u' [m/s²]
    s = 0.25*L/np.linalg.norm(F[:2])                                 # arrow length on screen [km per m/s²]
    arrow[0] = b.arrow(rot[0, i]/1e3, rot[1, i]/1e3, s*F[0], s*F[1], width=0.4, color=COLORS["amber"])
    return pa, pb, mer                                        # the artists that changed in this frame


show_animation(animate(update, frames=frames, fig=fig, interval=60))                                # build the movie from update() and play it
""")
see_read_change(
    see="Left: a straight rose line over a turning floor; right: the same throw curving to the right, the amber arrow always at right angles to it.",
    read=r"""The two pictures are one motion; the rotating observer needs the amber arrow $-2\boldsymbol\Omega\times\mathbf u'$ to
    explain the curve (plus a small centrifugal push). The grey parabola Ωut² is the small-angle estimate of D17 — good for
    the first part of the hour, too big later.""",
    change="…the Southern Hemisphere (Ω pointing down relative to the local vertical): the right panel bends left.")
remind([("plt.quiver", "`ax.quiver(X, Y, U, V)` draws an arrow field (Ch. 2 P78).")])
nb.figure(r"""
ang = np.linspace(0, 2*np.pi, 12, endpoint=False)                   # 12 directions round the centre [rad]
Xg = np.concatenate([R*np.cos(ang) for R in (8e4, 1.6e5)])           # points on two rings, 80 and 160 km out [m]
Yg = np.concatenate([R*np.sin(ang) for R in (8e4, 1.6e5)])
fig, axs = plt.subplots(1, 2, figsize=(10, 4.8))                                        # the figure and its panels
Om_NH = np.array([0, 0, ch04.OMEGA_EARTH])                           # Northern-Hemisphere pole: Ω up
for ax, sense in zip(axs, ("high", "low")):  # a high and a low
    U, V = ch04.high_low_flow(Xg, Yg, U_R=1.0, sense=sense)          # radial outflow (high) or inflow (low) [m/s]
    F = np.array([ch04.coriolis_acceleration(Om_NH, np.array([u_, v_, 0.0])) for u_, v_ in zip(U, V)])   # −2Ω × u at each point [m/s²]
    ax.quiver(Xg/1e3, Yg/1e3, U, V, color=COLORS["muted"], scale=10, width=0.006, label="flow u")   # the flow (grey)
    ax.quiver(Xg/1e3, Yg/1e3, F[:, 0], F[:, 1], color=COLORS["amber"], scale=1.5e-3, width=0.006,
              label="−2Ω × u (Coriolis force per mass)")             # Coriolis force (amber), at right angles to the flow
    ax.plot(0, 0, "o", color=COLORS["ink"]); ax.text(8, 8, "H" if sense == "high" else "L", fontsize=14)   # the pressure centre
    ax.set_xlim(-230, 230); ax.set_ylim(-230, 230); ax.legend(fontsize=7, loc="lower right")               # fixed limits, legend
    ax.set_aspect("equal"); ax.set_xlabel("x [km]"); ax.set_title(f"flow {'out of a high' if sense == 'high' else 'into a low'} (NH): amber turns it "  # axis labels with units
                                                               + ("clockwise" if sense == "high" else "counter-clockwise"), fontsize=9)
axs[0].set_ylabel("y [km]")                                        # axis labels with units
savefig(fig, "ch04", "highs_lows"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Grey arrows pointing out (high) or in (low); amber arrows at right angles to them.",
    read="The amber arrows turn outflow clockwise and inflow counter-clockwise in the Northern Hemisphere — the sense of the winds round highs and lows on a weather map.",
    change="…Southern Hemisphere: every amber arrow reverses.")
nb.plotly(r"""
from plotly.subplots import make_subplots                            # a 3-D scene next to a 2-D plot
boost = 18.0                                                         # Ω exaggerated so that Ω²a ≈ g/3 is visible
fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "xy"}]], column_widths=[0.62, 0.38],  # a 3-D scene next to a 2-D plot
                    subplot_titles=("exaggerated Ω: g_n (blue), Ω²R (purple), g_e (black)", "true tilt of the plumb line"))
nn = 40 if not FAST else 24
th_, ph_ = np.meshgrid(np.linspace(0, np.pi, nn), np.linspace(0, 2*np.pi, nn))   # sphere parameters
fig.add_trace(go.Surface(x=np.sin(th_)*np.cos(ph_), y=np.sin(th_)*np.sin(ph_), z=0.85*np.cos(th_), opacity=0.25,  # add a trace to the plotly figure
                         colorscale=[[0, COLORS["blue"]], [1, COLORS["blue"]]], showscale=False, name="flattened equipotential"), 1, 1)  # a translucent flattened Earth
gn, cen_scale = 0.5, 0.5*(boost*ch04.OMEGA_EARTH)**2*ch04.EARTH_RADIUS/9.80   # arrow lengths (g_n drawn as 0.5)
for lat in (0, 30, 60, 90):                                          # four latitudes on one meridian
    la = np.deg2rad(lat); P = np.array([np.cos(la), 0, np.sin(la)])  # the point on the unit sphere
    g_arrow = -gn*P                                                   # gravitation toward the centre
    c_arrow = np.array([cen_scale*np.cos(la)/1.0, 0, 0])              # centrifugal: away from the axis, ∝ R = cos(lat)
    for vec, col in ((g_arrow, COLORS["blue"]), (c_arrow, COLORS["accent"]), (g_arrow + c_arrow, COLORS["ink"])):  # gravitation, centrifugal and their sum
        fig.add_trace(go.Scatter3d(x=[P[0], P[0] + vec[0]], y=[P[1], P[1] + vec[1]], z=[P[2], P[2] + vec[2]], mode="lines",  # add a trace to the plotly figure
                                   line=dict(color=col, width=6), showlegend=False), 1, 1)
lats = np.linspace(0, 90, 91)                                        # latitude [°]
ge, dev = ch04.effective_gravity(np.deg2rad(lats))                   # true |g_e| and tilt from the radial direction
fig.add_trace(go.Scatter(x=lats, y=np.rad2deg(dev), mode="lines", line=dict(color=COLORS["ink"]), name="tilt [°]"), 1, 2)  # add a trace to the plotly figure
fig.update_xaxes(title_text="latitude [°]", row=1, col=2); fig.update_yaxes(title_text="tilt of g_e [°]", row=1, col=2)  # axis titles and layout
fig.update_layout(height=520, showlegend=False, scene=dict(aspectmode="data"))                      # axis titles and layout
fig.show()                                        # draw the interactive figure
""")
see_read_change(
    see="On the exaggerated Earth the black effective-gravity arrows lean away from the centre except at the pole and the equator; on the right the true tilt peaks near 45°.",
    read="Plumb lines do not point at the Earth's centre except at the poles and the equator; the real tilt is small (≈ 0.1° at 45°), but it is why the Earth's surface is a flattened spheroid. (The right panel is a toy: a spherical Earth with uniform gravitation g_n = 9.80 m/s²; real gravity also depends on the flattening, so standard gravity at 45° is 9.80665 m/s².)",
    change="…a planet spinning twice as fast: the centrifugal arrows grow four times (Ω²) and so does the tilt.")
nb.explainer("rotating_frame_coriolis", heading="Why does a straight throw curve on a merry-go-round?", why=r"""
**Why interactive rather than a static figure:** the whole idea is that one motion looks different to two observers at the same instant; running both views on one clock, and changing Ω and the hemisphere, makes the apparent forces visible as they act.

One ball, one clock, two observers: the inertial view (straight path, table turning) beside the rotating view (curved
path) with the apparent forces of (4.45),
$\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho\big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\big]+\mu\nabla'^2\mathbf u'$,
drawn as arrows and bars. Modes: turntable, Earth's pole, flow out of a high, effective gravity.
""", tries=[
    "Preset 'Earth pole 1 h': compare the deflection readout with Ωut² = 9.45 km.",
    "Set Ω = 0: the rotating view becomes the inertial one — no amber arrow.",
    "Mode 'flow out of a high': switch to the Southern Hemisphere and watch the turning reverse.",
    "Derivation tab, D15 steps 6 and 9: each adds one half of the Coriolis arrow.",
])
whatif(r"""
…we asked not for forces but for energy? Every force above does work at some rate — except Coriolis, which is always
perpendicular to the motion. The next block follows the energy of a fluid particle and finds where viscosity sends it
(C10).
""")
# =====================================================================================================================
# A.8 §4.8 Conservation of Energy — R07, R08, C10
# =====================================================================================================================
nb.section("4.8", "Conservation of Energy", intro="""
**What is this section about?** The first law of thermodynamics for a moving fluid particle. We derive the total-energy
equation, split off its mechanical part (Cauchy's equation dotted with u) and are left with an equation for the internal
energy whose viscous term — the dissipation ε — can only be positive: viscosity turns motion into heat and never the
reverse.
""")
nb.recap("R07", "A symmetric tensor contracted with an antisymmetric one gives zero", r"""
For symmetric σ and antisymmetric A, $\sigma_{ij}A_{ij}=0$ (Ch. 2 §2.10: swap the dummy indices, use σ_ji = σ_ij and
A_ji = −A_ij, so the sum equals minus itself). Chapter 2 flagged that this is exactly what makes the viscous work depend
only on the strain rate — used in step 7 of D22 below.
""", where="Ch. 2 §2.10")
nb.code(r"""
rng = np.random.default_rng(0)                                       # reproducible random numbers (Ch. 1 P10)
B = rng.normal(size=(3, 3)); sig = B + B.T; A = B - B.T              # a symmetric and an antisymmetric tensor
print(tensors.symmetric_double_contraction(sig, A))                  # σ_ij A_ij = 0 (to rounding)
""")
nb.recap("R08", "The Gibbs relation", r"""
For a substance in equilibrium $T\,ds=de+p\,dv$ *(Eq. 1.18)*, with v = 1/ρ the specific volume. Applied along the path of
a fluid particle (each particle is in local equilibrium) it reads $\frac{De}{Dt}=T\frac{Ds}{Dt}-p\frac{D(1/\rho)}{Dt}$
*(Eq. 4.61)* — the bridge from energy to entropy.
""", where="Ch. 1 §1.8")
core("C10", "Where the energy goes: internal energy and viscous dissipation (4.57)", """
Stir a cup of coffee and let go. The swirl dies away. Its kinetic energy has to go somewhere — where, and can it ever
come back?
""", eqs=("4.57",))
problem("""
Energy is conserved, but not every kind of energy is equally useful. The swirl's kinetic energy ends up as a (tiny)
warming of the coffee; the reverse — coffee spontaneously starting to swirl while cooling — never happens. In the
atmosphere the same one-way conversion ends the kinetic-energy cycle: winds are driven by heating differences and die by
friction, which returns the energy as heat. Lubricated bearings get hot for the same reason, and in turbulence (Ch. 12)
the rate of this conversion, ε, is the most important single number.
""")
idea("""
total energy e + ½u²:  rate = work of gravity + work of surface stresses − heat conducted out        (4.53)
kinetic ½u² alone:     Cauchy (4.24) · u  = work of gravity + work of the NET pressure and viscous forces (4.56)
subtract:  internal e:  De/Dt = −p Dv/Dt  +  ε  −  (1/ρ)∇·q                                          (4.57)
                                compression  dissipation ≥ 0  conduction
""", r"""
Stress does two kinds of work: *force work* (the net force times the velocity) speeds particles up or slows them down;
*deformation work* (stress times the rate of deformation) changes their shape — and the viscous part of that is lost to
heat.
""")
remind([("internal and kinetic energy per unit mass", r"e is the molecular (thermal) energy per kg, ½|u|² the bulk kinetic energy per kg (Ch. 1 P32).")])
P("P127", "power of a force and heat flux through a surface", r"""
A force F acting on something moving at velocity u does work at the rate F·u (watts). Per unit area of a surface the
stress does work at f·u; per unit volume gravity does work at ρg·u. Heat crossing a surface with outward normal n at flux
q (W/m²) leaves at the rate q·n per unit area — so a *minus* sign appears when we count heat *gained*.
""", code="""
f = np.array([2.0, 0.0, 0.0]); u = np.array([3.0, 4.0, 0.0])     # traction [Pa], velocity [m/s]
print(f @ u)                                                      # 6.0 W/m²: only the part of u along f counts
q = np.array([0.0, 0.0, 50.0]); n = np.array([0, 0, 1.0])         # heat flux [W/m²] and a lid's outward normal
print(-q @ n)                                                     # −50 W/m²: heat leaving through a lid
""")
note("N66", "First law for a material volume", r"""
(Ch. 1's first law $de=\delta q+\delta w$ *(1.10)* per unit time: heat in + work done = rise of internal + kinetic
energy):
""", equation=EQ["4.46"], ref="4.46")
D("D19", ref="4.53")
note("N67", "The book's (4.47)", r"D19's step 2: the transport theorem with F = ρ(e + ½|u|²), b = u.",
     equation=r"\int_{V}\frac{\partial}{\partial t}\Big(\rho e+\frac\rho2\lvert\mathbf u\rvert^2\Big)dV+\int_{A}\Big(\rho e+\frac\rho2\lvert\mathbf u\rvert^2\Big)(\mathbf u\cdot\mathbf n)dA=\int_{V}\rho\mathbf g\cdot\mathbf u\,dV+\int_{A}\mathbf f\cdot\mathbf u\,dA-\int_{A}\mathbf q\cdot\mathbf n\,dA", ref="4.47")
note("N68", "The energy budget of any control volume", r"""
The same coincidence move as D01 and D05 gives it (stated, not re-derived). `ch04.energy_budget` evaluates it (used for
nozzles and shocks in Ch. 15).
""", equation=EQ["4.48"], ref="4.48")
nb.note(r"""**Gauss on each surface term** `N69` `N70` `N71` `N72` — steps 3–7 of D19 are the book's (4.49)–(4.52): Gauss on
the energy flux, on the stress work and on the heat flux, then everything under one integral. The heat-flux line is the
one below.

> ⚠️ **Book slip:** (4.51) prints dA inside its two volume integrals; Gauss' theorem gives dV (as (4.52) then correctly
> uses).""", equation=r"\int_A\mathbf q\cdot\mathbf n\,dA=\int_Aq_in_i\,dA=\int_V\nabla\cdot\mathbf q\,dV=\int_V\frac{\partial q_i}{\partial x_i}dV", ref="4.51")
note("N73", "Total-energy equation", r"(D19's result), the conservative form used by compressible codes (Ch. 10, 15). $u_j^2=u_1^2+u_2^2+u_3^2$ (summed).",
     equation=EQ["4.53"], ref="4.53")
D("D20", ref="4.55")
note("N74", "Stress work = deformation work + force work", r"""
(D20's steps 5–6). The first bracket deforms particles (and heats them), the second changes their kinetic energy. The
same split, applied to fluctuations, is the turbulent kinetic-energy budget of Ch. 12. Code: `ch04.stress_work_split`
adds the four parts back to $\partial(\tau_{ij}u_j)/\partial x_i$ on a test field (below).
""", equation=r"\frac{\partial}{\partial x_i}(\tau_{ij}u_j)=\tau_{ij}\frac{\partial u_j}{\partial x_i}+u_j\frac{\partial\tau_{ij}}{\partial x_i}=\Big(-p\frac{\partial u_j}{\partial x_j}+\sigma_{ij}\frac{\partial u_j}{\partial x_i}\Big)+\Big(-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}\Big)", ref="4.54")
nb.code(r"""
u_tg, p_tg = ch04.exact_field("taylor_green", U0=1.0, k=1.0, nu=1e-3, rho=1000.0)   # a smooth test field and its pressure
sig_tg = ch04.newtonian_viscous_stress_field(u_tg, mu=1.0)           # σ(x, t) of that field for μ = 1 Pa s
sw = ch04.stress_work_split(p_tg, sig_tg, u_tg, np.array([0.3, 0.7]), 0.0)   # the four pieces of (4.54) at a point of this 2-D field [W/m³]
print({k: round(float(v), 8) for k, v in sw.items()})                # they add up to `total`: residual ≈ 0
""")
note("N75", "Total energy following a particle", "(D20's result; the book leaves it to an exercise).", equation=EQ["4.55"], ref="4.55")
remind([("chain rule", r"$\frac{d}{dt}f(g(t))=f'(g)\,\dot g$ — here $\frac{D}{Dt}(\tfrac12u^2)=u\frac{Du}{Dt}$ (Ch. 1 P49).")])
P("P128", "chain rule for the kinetic energy", r"""
D/Dt obeys the ordinary chain rule, so $\frac{D}{Dt}\big(\tfrac12u_j^2\big)=u_j\frac{Du_j}{Dt}$ (sum over j): the rate of
change of kinetic energy per unit mass is the velocity dotted with the acceleration. (A one-line reminder of Ch. 1's
chain rule P49, kept as a primer because D21 hinges on it.)
""", code="""
t = sp.symbols('t'); u = [sp.Function(f'u{j}')(t) for j in (1, 2, 3)]            # a velocity following a particle
print(sp.simplify(sp.diff(sum(v**2 for v in u)/2, t) - sum(v*sp.diff(v, t) for v in u)))   # 0
""")
D("D21", ref="4.56")
nb.code(r"""
u_tg, p_tg = ch04.exact_field("taylor_green", U0=1.0, k=1.0, nu=1e-3, rho=1000.0)   # an exact Navier–Stokes solution (ρ = 1000, μ = 1)
sig_tg = ch04.newtonian_viscous_stress_field(u_tg, mu=1.0)           # its viscous stress σ(x, t) [Pa]
ke = ch04.kinetic_energy_budget(1000.0, u_tg, p_tg, sig_tg, g=(0.0, 0.0), x=np.array([0.3, 0.7]), t=0.2)   # the terms of (4.56) [W/m³]
print({k: round(float(v), 6) for k, v in ke.items()})               # ρD(½u²)/Dt = pressure work + viscous-force work; residual ≈ 0
""", explain=r"""
The mechanical-energy equation $\rho\frac{D}{Dt}(\tfrac12u_j^2)=\rho g_ju_j-u_j\frac{\partial p}{\partial x_j}+u_j\frac{\partial\sigma_{ij}}{\partial x_i}$
*(4.56)* term by term on the decaying Taylor–Green vortices: the kinetic-energy rate (−240.0 W/m³) is the pressure work
(−239.1) plus the work of the net viscous force (−0.86); the residual is ≈ 10⁻⁹ — the D21 check.
""")
note("N76", "Mechanical-energy equation", r"""
— Cauchy's equation dotted with u.

> ⚠️ The book says "multiply (4.22) by u_j"; (4.22),
> $\frac{\partial}{\partial t}(\rho u_j)+\frac{\partial}{\partial x_i}(\rho u_iu_j)=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$,
> is the flux form, from which you must also subtract u_j × continuity. Starting from Cauchy's (4.24),
> $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$, is one step shorter (D21 does this). Mean and
> turbulent kinetic energy: Ch. 12.
""", equation=EQ["4.56"], ref="4.56")
D("D22", ref="4.57")
note("C10", "Internal energy — the thermal half of the energy budget", r"""
C10's equation, restated. It is the first law of Ch. 1 for a particle, per unit time: compression work −p dv, viscous
heating, heat conduction. Only the first and last can change sign.
""", equation=EQ["4.57"], ref="4.57")
P("P129", "completing the square for tensors", r"""
The double sum $A_{ij}A_{ij}=\sum_{i,j}A_{ij}^2$ is a sum of squares, so it is ≥ 0 and zero only when every component is.
To show an expression is never negative, rewrite it as such sums: here
$S_{ij}S_{ij}-\tfrac13S_{mm}^2=\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)$ (expand the
right side and use $\delta_{ij}\delta_{ij}=3$).
""", code="""
S = np.random.default_rng(3).normal(size=(3, 3)); S = (S + S.T)/2          # a random symmetric S
D = S - np.trace(S)/3*np.eye(3)                                             # its deviatoric part
print(np.isclose((S*S).sum() - np.trace(S)**2/3, (D*D).sum()), (D*D).sum() >= 0)   # True True
""")
D("D23", ref="4.58")
note("N77", "The dissipation rate", r"""
— kinetic energy turned into heat per unit mass and time. Squares with positive coefficients: ε ≥ 0 whenever μ ≥ 0 and
μ_v ≥ 0. Incompressible: $\varepsilon=2\nu S_{ij}S_{ij}$. It grows with the square of the velocity gradients — hot
bearings, glowing re-entry shields, and turbulence (Ch. 12).
""", equation=r"\varepsilon\equiv\frac1\rho\sigma_{ij}S_{ij}=2\nu\Big(S_{ij}-\tfrac13\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big)^2+\frac{\mu_v}{\rho}\Big(\frac{\partial u_m}{\partial x_m}\Big)^2\ge0", ref="4.58")
note("N78", "The Newtonian viscous stress", "used there (the σ part of (4.37)):", equation=EQ["4.59"], ref="4.59")
note("N79", "With Fourier's law", r"""
$\mathbf q=-k\nabla T$ *(1.2)* (Ch. 1 §1.5) the internal energy equation becomes the line below. **Ledger closes:**
continuity (1) + Navier–Stokes (3) + energy (1) + two equations of state = 7; unknowns ρ, e, p, T, u_j = 7.
""", equation=EQ["4.60"], ref="4.60")
nb.code(r"""
c = ch04.closure_count("full")                                       # the complete set
ledger.loc["full"] = [c["equations"], c["unknowns"]]                  # the last row of the §4.1 table
print(ledger)                                                        # cauchy 6/13, navier_stokes 4/5, barotropic 5/5, full 7/7
""")
note("N80", "The entropy equation", r"""
Combining $\frac{De}{Dt}=-p\frac{Dv}{Dt}+\frac1\rho\sigma_{ij}S_{ij}-\frac1\rho\frac{\partial q_i}{\partial x_i}$ *(4.57)*
with the Gibbs relation (R08) and the quotient rule gives it (stated; the moves are the book's):
""", equation=EQ["4.62"], ref="4.62")
note("N81", "Entropy production", r"""
The last two terms of (4.62); with $\mathbf q=-k\nabla T$ they are the line below. The second law demands it for every
flow, so **μ ≥ 0, μ_v ≥ 0, k ≥ 0**.

> ⚠️ The book writes "μ, κ, k > 0": κ there is a leftover symbol for the bulk viscosity μ_v (κ is the thermal diffusivity
> elsewhere), and the law requires ≥, not >. An inviscid, non-conducting flow has Ds/Dt = 0: every particle keeps its
> entropy (isentropic) — used for (4.78) and in Ch. 15.
""", equation=EQ["4.63"], ref="4.63")
nb.worked_example("water sheared between plates (plane Couette, U = 1 m/s, gap h = 1 mm)", r"""
μ = 1.0×10⁻³ Pa s, ρ = 1000 kg/m³, k = 0.6 W/(m K), both walls held at T₀.

1. Velocity u = Uy/h, shear rate du/dy = 1000 s⁻¹.
2. S₁₂ = 500 s⁻¹, so ε = 2ν(S₁₂² + S₂₁²) = 2 × 10⁻⁶ × 2 × 500² = 1.0 W/kg, ρε = μ(U/h)² = 1000 W/m³.
3. Work done by the moving wall per unit area: τU = μ(U/h)U = 1.0 × 1 = 1.0 W/m²; heat generated in the gap:
   ρε × h = 1000 × 10⁻³ = 1.0 W/m² — the same (all the wall's work becomes heat).
4. Steady temperature: k T″ = −ρε ⇒ a parabola, peak rise at mid-gap
   $\Delta T_{max}=\rho\varepsilon h^2/(8k)=\mu U^2/(8k)=10^{-3}/4.8=2.1\times10^{-4}$ K; each wall conducts away 0.5 W/m².
5. Entropy production at mid-gap: ε/T = 1/293 = 3.4×10⁻³ W/(kg K) > 0.
""")
nb.code(r"""
G = np.array([[0, 1000.0, 0], [0, 0, 0], [0, 0, 0]])                # velocity gradient of the Couette flow [1/s]
print(ch04.dissipation_rate(G, 1000.0, 1e-3, form="both"))           # ε by contraction and by sum of squares: (1.0, 1.0) W/kg
c = ch04.couette_heating(np.linspace(0, 1e-3, 5), U=1.0, h=1e-3, mu=1e-3, k=0.6, T0=293.15)   # steady heated Couette flow
print(c["eps"], np.round((c["T"] - 293.15)*1e4, 3), "×1e-4 K")      # ρε = 1000 W/m³ everywhere; T − T0 at five heights
print(c["work_in"], c["heat_out"], c["q_bottom"], c["q_top"], c["dT_max"])   # 1 W/m² in = 1 W/m² out (0.5 through each wall)
print(ch04.entropy_production(np.zeros(3), 293.15, 0.6, 1000.0, 1.0))   # at mid-gap ∇T = 0: ε/T = 3.41e-3 W/(kg K)
print(ch04.dissipation_rate(G, 1000.0, -1e-3))                       # a "negative viscosity" gives ε < 0 — forbidden
""", explain=r"""
1. ε by contraction $\sigma_{ij}S_{ij}/\rho$ and by the sum of squares of (4.58) — identical.
2. The steady Couette solution with viscous heating and its energy bookkeeping: the wall's work (1 W/m²) all leaves as
   heat, half through each wall; the peak rise is 2.08×10⁻⁴ K.
3. The entropy production at mid-gap (no temperature gradient there, so all of it from ε).
4. With μ < 0 dissipation would be negative — the second law rules it out.
""")
nb.check_agree(r"""
sig = ch04.viscous_stress(G, 1e-3); S = 0.5*(G + G.T)                # viscous stress [Pa] and strain rate [1/s]
eps = sum(sig[i, j]*S[i, j] for i in range(3) for j in range(3)) / 1000.0   # ε = σ_ij S_ij/ρ by the explicit double sum
assert np.isclose(eps, 1.0)                                          # 1.0 W/kg
rng = np.random.default_rng(11)                                      # 1000 random compressible velocity gradients
Gs = rng.normal(size=(1000, 3, 3))
both = np.array([ch04.dissipation_rate(g_, 1000.0, 1e-3, mu_v=5e-4, form="both") for g_ in Gs])   # (contraction, sum of squares)
print(np.abs(both[:, 0] - both[:, 1]).max(), both.min())             # agreement ~1e-19 and the smallest ε ≥ 0
assert np.allclose(both[:, 0], both[:, 1], atol=1e-12) and np.all(both >= 0)
""")
nb.figure(r"""
yy = np.linspace(0, 1e-3, 201)                                       # across the gap [m]
fig, axs = plt.subplots(1, 4, figsize=(13, 3.6), gridspec_kw=dict(width_ratios=[1, 1, 1, 0.9]))    # the figure and its panels
for dpdx, ls, lab in ((0.0, "-", "Couette"), (-2000.0, "--", "Couette + dp/dx = −2000 Pa/m")):  # pure Couette and with a pressure gradient
    c = ch04.couette_heating(yy, U=1.0, h=1e-3, mu=1e-3, k=0.6, dpdx=dpdx)   # the steady heated flow
    axs[0].plot(c["u"], yy*1e3, ls, color=COLORS["teal"], label=lab)          # velocity
    axs[1].plot(c["eps"], yy*1e3, ls, color=COLORS["rose"])                    # ρε [W/m³]
    axs[2].plot((c["T"] - 293.15)*1e3, yy*1e3, ls, color=COLORS["rose"])       # T − T0 [mK]
axs[0].set_xlabel("u [m/s]"); axs[0].set_ylabel("y [mm]"); axs[0].legend(fontsize=7)                # axis labels with units
axs[1].set_xlabel("ρε [W/m³]"); axs[2].set_xlabel("T − T₀ [mK]")                                    # axis labels with units
c = ch04.couette_heating(yy, U=1.0, h=1e-3, mu=1e-3, k=0.6)            # the pure Couette budget
ins = axs[3]                                                          # panel (d): the energy budget of pure Couette flow
ins.bar([0, 1, 2, 3], [c["work_in"], -c["q_bottom"], -c["q_top"], c["work_in"] - c["heat_out"]],
        color=[COLORS["teal"], COLORS["rose"], COLORS["rose"], COLORS["ink"]])   # work in, heat out (bottom, top), sum
ins.set_xticks([0, 1, 2, 3], ["wall work\nin", "heat out\nbottom", "heat out\ntop", "sum"])   # label the four bars
ins.tick_params(axis="x", labelsize=8); ins.set_ylabel("[W/m²]"); ins.axhline(0, color=COLORS["muted"], lw=0.6)   # axis label, zero line
axs[0].set_title("(a) velocity", fontsize=10); axs[1].set_title("(b) dissipation", fontsize=10)     # the panel's message
axs[2].set_title("(c) temperature rise", fontsize=10); ins.set_title("(d) budget (Couette)", fontsize=10)   # the panels' messages
savefig(fig, "ch04", "couette_heating"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="A straight velocity profile, a flat dissipation strip, a parabolic temperature bump of a fifth of a millikelvin (solid); with the pressure gradient (dashed) the dissipation piles up near the lower wall and the bump shifts toward it.",
    read="""The moving wall pumps in 1 W per m²; all of it becomes heat inside and leaves by conduction through the two walls
    (the budget bars add to zero); the bump is where the heat has farthest to go.""",
    change="…oil instead of water (μ 100 × larger, k 4 × smaller): the bump is 400 × larger — 0.08 K — which is why bearings are cooled.")
nb.plotly(r"""
yy = np.linspace(0, 1e-3, 101)                                       # across the gap [m]
def heating(U):                                                      # wall speed U [m/s]
    w = ch04.couette_heating(yy, U=U, h=1e-3, mu=1e-3, k=0.6)          # water
    o = ch04.couette_heating(yy, U=U, h=1e-3, mu=0.1, k=0.15)          # a light oil: μ ×100, k ÷4
    return {"oil T − T₀ [K]": (yy*1e3, o["T"] - 293.15), "water T − T₀ × 400 [K]": (yy*1e3, 400*(w["T"] - 293.15))}  # the curves for this slider value
fig = slider_figure(heating, "U", np.logspace(-1, 1, 20 if not FAST else 10), unit="m/s", xlabel="y [mm]",  # precompute every slider position (works on the web page)
                    ylabel="temperature rise [K]", title="Viscous heating grows as U² (water curve ×400 lies on the oil curve)")  # labels, title and a fixed range
recolor(fig, {"oil T − T₀ [K]": COLORS["rose"], "water T − T₀ × 400 [K]": COLORS["teal"]}, {"water T − T₀ × 400 [K]": "dash"}).show()  # house colours for each trace, then draw
""")
see_read_change(
    see="Two identical parabolas (water scaled by 400) that grow by a factor 100 when the slider moves from 1 to 10 m/s.",
    read=r"""The peak rise is $\mu U^2/(8k)$: at 10 m/s water warms by only 0.02 K, oil by 8 K; the shape never changes, only
    its size, as U².""",
    change="…the gap halved at the same U: the shear rate doubles but the path for the heat halves — the peak rise $\mu U^2/8k$ does not change.")
nb.explainer("viscous_dissipation_heating", heading="Where does the energy go when viscosity stops a flow?", why=r"""
**Why interactive rather than a static figure:** the heating is a process in time — the temperature relaxes to its parabola while the energy bars fill — and changing fluid and wall speed shows the U² and μ/k scalings as you go.

A sheared channel on one clock: the velocity profile, the dissipation ε(y) of
$\varepsilon=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}\rho S_{mm}^2$ (4.58), and the temperature relaxing to
its steady parabola, with energy-budget bars (wall work in = heat out) and the entropy production. The start-up of the
temperature is a sine series — the same method as Ch. 1's Couette start-up. A "negative μ" preset breaks the second law
on purpose.
""", tries=[
    "Water Couette at 1 m/s: read the peak temperature rise in Explain (2.1×10⁻⁴ K).",
    "Oil bearing preset: the same flow heats 400 times more.",
    "Press ▶ and watch T(y) grow to its steady ghost; the end card says shear work in = heat out.",
    "Preset 'negative μ': ε and the entropy production turn negative and the status flags the second law.",
])
whatif(r"""
…the flow were frictionless and heat did not flow? Then ε = 0, q = 0, each particle keeps its entropy, and the energy and
momentum equations can be integrated along streamlines — the Bernoulli family (C11, C12).
""")
# =====================================================================================================================
# A.9 §4.9 Special Forms of the Equations — C11, R09, C12, R10, C13
# =====================================================================================================================
nb.section("4.9", "Special Forms of the Equations", intro="""
**What is this section about?** Useful consequences of the equations under special conditions: the angular-momentum
principle, four Bernoulli equations (each with its own hypotheses), gravity absorbed into the pressure, and the
Boussinesq approximation that ocean and atmosphere models use.
""")
nb.pointer(r"""The angular-momentum principle, $d\mathbf H/dt=\mathbf M$ (4.64) and its control-volume form (4.65), with the
lawn sprinkler (Ex. 4.6), was taught with the momentum budgets in C04 above (notes N82–N84).""")
nb.pointer(r"""The pitot tube, $\lvert\mathbf u\rvert=\sqrt{2(p_2-p_1)/\rho}$, stagnation and dynamic pressure
$p+\tfrac12\rho\lvert\mathbf u\rvert^2$, and the orifice $u=\sqrt{2gh}$ were taught with Bernoulli's equation
$\tfrac12U^2+gz+p/\rho=$ const (4.19) in C05 above (notes N103–N105).""")

# ---- C11 ---------------------------------------------------------------------------------------------------------
core("C11", "Bernoulli is constant along what, exactly? (4.71)", """
In a bathtub whirlpool the water surface dips in the middle, although every particle goes round at a steady speed. If
Bernoulli held across the circles the surface would be level. Along what, exactly, is Bernoulli's sum constant?
""", eqs=("4.71",))
problem("""
C05 derived Bernoulli along one streamline of a steady, frictionless, constant-density flow. People use it far more
widely — across streamlines, in gases, in unsteady flows — sometimes correctly, sometimes not. We derive it again from
Euler's equation, in a way that shows precisely what it is constant along, and why irrotational flows are special.
""")
idea("""
Euler + Lamb identity:   ∂u/∂t + ∇B = u × ω         B = ½u² + ∫dp/ρ + gz          (4.69)
steady:                  ∇B = u × ω   ⊥ u and ⊥ ω   →  B constant along streamlines AND vortex lines   (4.71)
irrotational (ω = 0):    ∇B = 0                     →  B constant everywhere                          (4.72)
""")
note("N85", "Bernoulli equations are not new laws", r"""
each is a consequence of Navier–Stokes (4.38) or the energy equation (4.60) under stated conditions — so each comes with
its own list of hypotheses (the table at the end of this block).
""")
note("N86", "Start from Euler's equation", r"""
$\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g$ *(4.41)*, with gravity from its potential Φ = gz,
$\mathbf g=-\nabla\Phi$ *(4.18)*:
""", equation=EQ["4.66"], ref="4.66")
remind([("fundamental theorem of calculus with a variable limit", r"$\frac{d}{dp}\int_{p_o}^{p}f(p')\,dp'=f(p)$ (Ch. 2 P84).")])
note("N87", "Barotropic flow", r"""
— density a function of pressure alone, ρ = ρ(p) (constant density, isothermal gas, isentropic gas). Then dp/ρ is an exact
differential and the line below holds (the fundamental theorem of calculus with a variable upper limit).

> ⚠️ The prime on p′ here is only the integration variable. Closed forms (`ch04.pressure_function`): constant ρ →
> (p − p_o)/ρ; isothermal → RT ln(p/p_o); isentropic → the enthalpy difference h − h_o. A *baroclinic* fluid (ρ depends on
> T too, like the real atmosphere) has no such single function — the root of Ch. 13's thermal wind.
""", equation=EQ["4.67"], ref="4.67")
remind([("permutations and the alternating tensor ε_ijk", r"$\varepsilon_{ijk}$ = +1 for cyclic (123), −1 for anticyclic, 0 with a repeat; $(\mathbf a\times\mathbf b)_i=\varepsilon_{ijk}a_jb_k$ (Ch. 2 P72).")])
D("D24", ref="4.69")
nb.code(r"""
x, y, z = sp.symbols('x y z')                                        # coordinates
u_any = [y**2*z, sp.sin(x)*z, x*y**3 + sp.exp(z)]                    # an arbitrary smooth 3-D velocity field
print(ch04.lamb_identity_sym(u_any, (x, y, z)))                      # (u·∇)u − [−u × ω + ∇(½u²)] = [0, 0, 0]
""", explain=r"""
The Lamb identity $u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}(\tfrac12u_i^2)$
*(4.68)* checked symbolically for a general field (D24 steps 3–7); the numerical check at a point of the Rankine vortex
follows below.
""")
note("N88", "The Lamb identity", r"""
(D24 steps 3–7; the book leaves it to an exercise). The advective acceleration = gradient of the kinetic energy minus the
"Lamb vector" u × ω. Ch. 5 starts the vorticity equation from here.
""", equation=EQ["4.68"], ref="4.68")
note("N89", "Euler in Bernoulli-function form", r"(D24's result). The bracket is the **Bernoulli function** B (`ch04.bernoulli_function`).",
     equation=EQ["4.69"], ref="4.69")
D("D25", ref="4.71")
note("N90", "Steady flow", r"""
∇B is normal to the surfaces B = const, and u × ω is perpendicular to both u and ω, so each B-surface (a *Lamb surface*)
is woven from streamlines and vortex lines. Picture: a Rankine vortex with an axial current — its Lamb surfaces are
cylinders carrying helical streamlines and axial vortex lines.
""", equation=EQ["4.70"], ref="4.70")
note("C11", "Bernoulli 1 (steady, inviscid, barotropic, conservative body force)", "C11's result, restated.",
     equation=r"\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{constant along streamlines and vortex lines}", ref="4.71")
note("N91", "Irrotational as well", r"""
(ω = 0 everywhere in a connected region): ∇B = 0 and B is one constant. Number (code below): B at about 900 random points of
the Ch. 3 cylinder flow agrees to 1e-15; inside a Rankine core B varies from circle to circle.
""", equation=r"\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{constant everywhere}", ref="4.72")
note("N92", "Persistence of irrotationality", r"""
An inviscid, barotropic flow in an inertial frame that starts irrotational stays irrotational (proved with Kelvin's
circulation theorem in Ch. 5). That is why the flow outside a thin boundary layer round a body is irrotational and (4.72)
applies there (Ch. 6, Ch. 9).
""")
nb.note(r"""**Energy form** `N94` `N95` — start instead from the total-energy equation (4.55) with σ = q = 0 and steady flow;
with $\mathbf g=-\nabla(gz)$ and steady continuity $\partial(\rho u_i)/\partial x_i=0$ it becomes the second line.""",
        equation=r"\begin{aligned}" + EQ["4.76"] + r"\quad&\text{(4.76)}\\ " + EQ["4.77"] + r"\quad&\text{(4.77)}\end{aligned}")
note("N96", "Bernoulli 3 (steady, inviscid, non-conducting; h = e + p/ρ)", r"""
It is an *energy* statement: it holds in a gas where kinetic and thermal energy trade places. For isentropic flow
dp/ρ = dh (Gibbs), so (4.71) and (4.78) agree. Number: air at 300 K moving at 100 m/s brought to rest heats to
$T_0=T+U^2/2C_p=300+10^4/2009=304.98$ K (Ch. 15).
""", equation=EQ["4.78"], ref="4.78")
note("N102", "Which Bernoulli? — a decision table", r"""
(`ch04.BERNOULLI_FORMS`, `ch04.which_bernoulli`):

| form | steady? | inviscid? | irrotational? | barotropic / ρ const? | constant along |
|---|---|---|---|---|---|
| (4.19) $\tfrac12U^2+gz+p/\rho$ | yes | yes | not needed | ρ const | one streamline |
| (4.71) $\tfrac12u^2+\int dp/\rho+gz$ | yes | yes | not needed | barotropic | streamlines and vortex lines |
| (4.72) same B | yes | yes | yes | barotropic | everywhere |
| (4.75) $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int dp/\rho+gz$ | **unsteady allowed** | yes | yes | barotropic | everywhere (a function of t absorbed into φ) |
| (4.78) $h+\tfrac12\lvert\mathbf u\rvert^2+gz$ | yes | inviscid and non-conducting | not needed | any | streamlines |
| (4.82)/(4.83) | unsteady | **viscous allowed** (ρ, μ const) | yes | ρ const | along a streamline at one instant / everywhere |
""")
remind([("polar coordinates as a moving basis", r"$\mathbf e_r$, $\mathbf e_\theta$ turn with the point; for circular motion the radial balance is $dp/dr=\rho u_\theta^2/r$ (Ch. 3 P105).")])
nb.worked_example("a Rankine vortex in water (Γ = 2π m²/s, core radius σ = 1 m, p∞ = 0)", r"""
1. Inside the core (r < σ) the fluid turns like a solid body: $u_\theta=\Gamma r/(2\pi\sigma^2)=r$ (s⁻¹ × m), vorticity
   ω = Γ/(πσ²) = 2 s⁻¹. Outside it is irrotational: $u_\theta=\Gamma/(2\pi r)=1/r$.
2. Pressure from the radial balance $dp/dr=\rho u_\theta^2/r$: outside p/ρ = −1/(2r²); inside p/ρ = −1 + r²/2 (matched
   at r = 1).
3. B = ½u_θ² + p/ρ: outside ½/r² − ½/r² = **0** everywhere; inside r²/2 − 1 + r²/2 = **r² − 1**: −1 at the centre, −0.75
   at r = 0.5, 0 at the edge.
4. Along a circle (a streamline) B is constant — (4.71) holds; across circles it changes inside the rotational core and
   not outside — $\tfrac12u_i^2+\int dp/\rho+gz$ = const everywhere *(4.72)* holds only where ω = 0.
5. Check D25: dB/dr = 2r = 1 at r = 0.5, and |u × ω| = u_θω = 0.5 × 2 = 1 ✓.
""")
nb.code(r"""
r = np.array([0.0, 0.5, 1.0, 2.0, 4.0])                              # radii [m]
print(ch04.rankine_bernoulli(r, 2*np.pi, 1.0)["B"])                  # B = −1, −0.75, 0, 0, 0 m²/s²
u_rk = ch03.vortex_velocity_field("rankine", Gamma=2*np.pi, sigma=1.0)   # the Rankine vortex as a 2-D field u(x, t)
x0 = np.array([0.5, 0.0])                                            # a point inside the core [m]
print(ch04.lamb_vector(u_rk, x0)[:2], ch04.lamb_identity_terms(u_rk, x0)["residual"])   # u × ω = (1, 0) outward; residual ≈ 0
u_cyl = ch03.cylinder_velocity_field(1.0, 1.0)                        # Ch. 3's cylinder flow (U = 1 m/s, a = 1 m)
pts = np.random.default_rng(0).uniform(-3, 3, size=(2, 1000)); pts = pts[:, np.hypot(*pts) > 1.05]   # points outside the body
B = ch04.bernoulli_along_line(u_cyl, lambda x, t: ch04.exact_solution("cylinder", x, U=1.0, a=1.0, rho=1.0)[1],  # call the tested chapter function
                              pts, rho=1.0, g=0.0)                   # B at every point (pressure from the exact solution)
print(f"spread of B over {pts.shape[1]} points: {B.max() - B.min():.1e}")   # ≈ 1e-16: one constant (4.72)
print(ch04.which_bernoulli(steady=True, viscous=False, irrotational=False, barotropic=True))   # ['4.71']
print(round(ch04.pressure_function(2e5, 1e5, "isothermal", T=300.0), 1))   # RT ln 2 = 59692 J/kg for an isothermal gas
print(round(ch04.stagnation_temperature(300.0, 100.0), 2))            # 304.98 K (N96)
""", explain=r"""
1. B across the vortex — varies inside, constant outside.
2. The Lamb identity at a point inside the core: u × ω points outward, equal to dB/dr = 1 m/s².
3. B at about 900 points of an irrotational flow — one constant.
4. The decision function for a steady, inviscid, barotropic rotational flow: only (4.71).
5. An isothermal gas's pressure function $\int_{p_o}^{p}dp'/\rho(p')=RT\ln(p/p_o)$.
6. The energy form's stagnation temperature.
""")
nb.check_agree(r"""
w = kinematics.vorticity(u_rk, x0, 0.0)                              # ω by central differences: (0, 0, 2) 1/s
u3 = np.array([*u_rk(x0, 0.0), 0.0])                                 # u at the point, as a 3-vector
h = 1e-5                                                             # step for the kinetic-energy gradient [m]
ke = lambda x: 0.5*np.sum(u_rk(x, 0.0)**2)                           # ½|u|²
grad_ke = np.array([(ke(x0 + e) - ke(x0 - e))/(2*h) for e in np.eye(2)*h])   # ∇(½|u|²)
adv = kinematics.acceleration(u_rk, x0, 0.0).advective               # (u·∇)u
print(adv, -np.cross(u3, w)[:2] + grad_ke)                           # the two sides of the Lamb identity (4.68)
assert np.allclose(adv, -np.cross(u3, w)[:2] + grad_ke, atol=1e-6)
assert np.allclose(np.cross(u3, w)[:2], ch04.lamb_vector(u_rk, x0)[:2], atol=1e-6)   # the library's u × ω
""")
nb.md('**What does the code above do?** It builds both sides of the Lamb identity at a point inside the Rankine core — (u·∇)u from the acceleration stencil, −u × ω + ∇(½|u|²) from `np.cross` and a central difference — and checks they agree with each other and with `ch04.lamb_vector`.')
nb.figure(r"""
rr = np.linspace(0, 4, 401)                                          # radius [m]
rb = ch04.rankine_bernoulli(rr, 2*np.pi, 1.0)                        # u_θ, p, B across the circles
fig, (a, b, c) = plt.subplots(1, 3, figsize=(13, 3.6), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1.4]))  # the figure and its panels
a.axvspan(0, 1, color=COLORS["accent"], alpha=0.12, label="rotational core")   # the core
a.plot(rr, rb["B"], color=COLORS["rose"], lw=2.5, label="B across circles")                         # draw the data
a.set_xlabel("r [m]"); a.set_ylabel("B [m²/s²]"); a.legend(fontsize=8); a.set_title("(a) B changes across streamlines only in the core", fontsize=9)  # axis labels with units
phi = np.linspace(0, 2*np.pi, 181)                                   # angle round a circle
for rad, lab in ((0.5, "r = 0.5 m (inside)"), (2.0, "r = 2 m (outside)")):  # one circle inside the core, one outside
    P = np.array([rad*np.cos(phi), rad*np.sin(phi)])                # points on the circle
    Bc = ch04.bernoulli_along_line(u_rk, lambda x, t: ch04.rankine_vortex_pressure(np.hypot(*x), 2*np.pi, 1.0, rho=1.0), P, rho=1.0, g=0.0)  # call the tested chapter function
    b.plot(np.rad2deg(phi), Bc, lw=2.5, label=lab)                  # B along one streamline
b.set_ylim(-1.1, 0.2); b.set_xlabel("angle [°]"); b.set_ylabel("B [m²/s²]"); b.legend(fontsize=8)   # axis labels with units
b.set_title("(b) B along a streamline: flat", fontsize=9)                                        # the panel's message
c.axis("off")                                                        # the decision table, (4.71) and (4.72) lit
rows = [["(4.19)", "1 streamline", "steady, inviscid, ρ const"], ["(4.71)", "streamlines + vortex lines", "steady, inviscid, barotropic"],  # the table's rows
        ["(4.72)", "everywhere", "… and irrotational"], ["(4.75)", "everywhere", "unsteady, irrotational"],  # table row
        ["(4.78)", "streamlines", "steady, adiabatic, energy form"]]                                # table row
tb = c.table(cellText=rows, colLabels=["form", "constant along", "needs"], loc="center", cellLoc="left",
             colWidths=[0.16, 0.42, 0.52])                           # a small table drawn in the panel
tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1, 1.4)   # small font so every entry fits
for k in (2, 3):                                                     # light the rows this vortex illustrates
    for j in range(3):  # column j of τ
        tb[k, j].set_facecolor("#fde2ec")
c.set_title("(c) which Bernoulli? (lit: this vortex)", fontsize=9)                                  # the panel's message
savefig(fig, "ch04", "rankine_bernoulli"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="(a) A parabola inside the core joining a flat line outside; (b) two flat lines; (c) the decision table with two rows lit.",
    read="""Along any streamline B is constant (4.71); across streamlines it changes only where there is vorticity — outside
    the core the flow is irrotational and one B holds everywhere (4.72). The whirlpool's dip is exactly the core's B deficit.""",
    change="…a uniform flow instead: B would be one constant everywhere; the whirlpool needs vorticity for its dip.")
nb.plotly(r"""
rr = np.linspace(0, 4, 300)                                          # radius [m]
def rankine(sig):                                                    # core radius σ [m], Γ = 2π m²/s
    rb = ch04.rankine_bernoulli(rr, 2*np.pi, sig)                                        # call the tested chapter function
    return {"B(r) [m²/s²]": (rr, rb["B"]), "u_θ(r) [m/s]": (rr, rb["u_theta"])}                     # the curves for this slider value
fig = slider_figure(rankine, "σ", np.linspace(0.2, 2.0, 19 if not FAST else 10), unit="m", xlabel="r [m]",  # precompute every slider position (works on the web page)
                    ylabel="B [m²/s²] · u_θ [m/s]", title="Rotational core: B varies across streamlines only where ω ≠ 0",  # labels, title and a fixed range
                    yrange=(-26, 6))                                        # labels, title and a fixed range
recolor(fig, {"B(r) [m²/s²]": COLORS["rose"], "u_θ(r) [m/s]": COLORS["teal"]}).show()               # house colours for each trace, then draw
""")
see_read_change(
    see="A rose dip inside the core and a teal speed profile peaking at the core's edge.",
    read="Shrinking σ deepens the dip, $B(0)=-\\Gamma^2/(4\\pi^2\\sigma^2)$, and narrows it; outside the core B is zero whatever σ is.",
    change="…Γ doubled: the speed doubles and the dip deepens four times (B ∝ Γ²).")
whatif(r"""
…the flow were irrotational but unsteady — water sloshing in a U-tube, a balloon pushed suddenly through water? Then the
∂u/∂t term of (4.69) stays, but it is a gradient too, and a new Bernoulli equation holds everywhere at each instant (C12).
""")

# ---- R09 + C12 ----------------------------------------------------------------------------------------------------
nb.recap("R09", "Velocity potential", r"""
An irrotational flow (ω = ∇×u = 0) in a simply connected region can be written as the gradient of a scalar,
$\mathbf u\equiv\nabla\phi$ *(Eq. 4.73, first met as $\boldsymbol\omega=0$ (3.17) in Ch. 3)*, because the curl of a
gradient is zero. φ is defined up to an added constant — or an added function of time only, which changes no velocity.
""", where="Ch. 3 §3.4")
core("C12", "Unsteady Bernoulli: the pressure of an accelerating flow (4.75)", """
Water sloshes back and forth in a U-tube. At the instant the column is momentarily still, the pressures at the two ends
still differ. Why — and by how much?
""", eqs=("4.75",))
problem("""
Steady Bernoulli says pressure is high where flow is slow. In an accelerating flow that is not enough: pushing fluid to
speed it up takes a pressure difference even at zero speed. Water waves (Ch. 7), a body starting to move (its "added
mass", Ch. 6), a bubble collapsing and sound (Ch. 15) all need the unsteady form — valid everywhere, but only for
irrotational flow.
""")
idea("""
u = ∇φ  →  Lamb term vanishes, ∂u/∂t = ∇(∂φ/∂t)  →  ∇[∂φ/∂t + ½|∇φ|² + ∫dp/ρ + gz] = 0
       →  the bracket is the same everywhere: a function of time B(t) only  →  absorb it into φ  →  (4.75)
U-tube:  p₁ − p₂ = ρ L dU/dt   (the ∂φ/∂t term, amber, carries the difference when U = 0)
""")
D("D26", ref="4.75")
note("N93", "The book's (4.74)", r"""
D26 steps 4–5.

> ⚠️ **Book slip (sign):** to absorb B(t) the potential must be redefined as $\phi\to\phi-\int_{t_o}^{t}B(t')dt'$; the
> printed "+" makes the bracket 2B(t) instead of 0 (step 6 and the code below). The result (4.75) is unaffected; the
> velocity ∇φ does not change either way.
""", equation=r"\nabla\Big[\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz\Big]=0,\quad\text{or}\quad\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=B(t)", ref="4.74")
note("C12", "Bernoulli 2 (unsteady, inviscid, irrotational, barotropic)", "C12's result, restated. Ch. 7 applies it at the free surface of every water wave.",
     equation=EQ["4.75"], ref="4.75")
note("N97", "Viscosity too?", r"""
Put the Lamb identity (4.68) and the viscous force $-\mu\nabla\times\boldsymbol\omega$ of
$\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega$ *(4.40)* into (4.39b):
""", equation=EQ["4.79"], ref="4.79")
note("N98", "With ω = 0 and constant ρ", r"""
**both** the Lamb term and the viscous term vanish — viscosity is present but exerts no net force on an irrotational flow.
Code: `ch04.viscous_irrotational_residual(u_cyl, (1.5, 0.7))` → 0 to round-off (printed below). (Viscosity still
*dissipates* energy in such a flow — ε ≠ 0 — which damps water waves slowly, Ch. 7.)
""", equation=EQ["4.80"], ref="4.80")
remind([("directional derivative along a streamline", r"$\partial F/\partial s=\hat{\mathbf s}\cdot\nabla F$, the rate of change of F along the unit tangent (Ch. 2 P75).")])
note("N99", "Along a streamline", "Dot (4.80) with a line element ds along a streamline and integrate from point 1 to point 2:",
     equation=EQ["4.81"], ref="4.81")
note("N100", "Bernoulli 4 (constant ρ and μ, irrotational, along a streamline at one instant)", r"""
Number: a water column L = 1 m long accelerating at dU/dt = 1 m/s² along a straight tube needs p₁ − p₂ = ρL dU/dt =
1000 Pa, even at the instant U = 0.
""", equation=EQ["4.82"], ref="4.82")
note("N101", "In potential form", r"(u = ∇φ, swap ∂/∂t and ∇) — the constant-ρ (4.75), now shown to hold with constant viscosity as well.",
     equation=EQ["4.83"], ref="4.83")
nb.worked_example("the U-tube and an accelerating sphere", r"""
1. **U-tube:** a water column of length L = 1 m with its two surfaces displaced ±h oscillates at ω = √(2g/L) = 4.43 rad/s
   (period 1.42 s). At the instant the column is still (U = 0) its acceleration is largest: with h₀ = 5 cm,
   dU/dt = ω²h₀ = 19.62 × 0.05 = 0.98 m/s², and (4.82),
   $\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s+(\tfrac12\lvert\mathbf u\rvert^2+gz+\tfrac p\rho)_2=(\dots)_1$, gives
   p₁ − p₂ = ρL dU/dt = 1000 × 1 × 0.98 = 981 Pa — exactly the weight of the 2h₀ = 10 cm height difference,
   ρg(2h₀) = 981 Pa ✓.
2. **Sphere** of radius a = 0.1 m accelerating from rest at 1 m/s² through still water: the unsteady term gives a surface
   pressure ½ρa cos θ dU/dt (high in front, low behind) whose net force is −½ρV dU/dt, V = 4πa³/3 = 4.19×10⁻³ m³: the
   sphere must push as if it carried half its volume of water with it — "added mass" 2.09 kg, force 2.09 N (Ch. 6).
""")
nb.code(r"""
c = ch04.u_tube_column(np.array([0.0, 0.3547, 0.7093]), L=1.0, h0=0.05)   # t = 0, a quarter and a half period [s]
print(np.round(c["h"], 4), np.round(c["U"], 4), np.round(c["dUdt"], 3), np.round(c["dp"], 1))   # level, speed, acceleration, p1 − p2
print(ch04.unsteady_streamline_bernoulli(np.full(11, 1.0), np.linspace(0, 1, 11), {"U": 0.0, "z": 0.0, "p": 1000.0},  # evenly spaced values
                                         {"U": 0.0, "z": 0.0, "p": 0.0}))   # residual of (4.82): 0 — 1000 Pa is what the acceleration needs
s = ch04.accelerating_sphere_pressure(np.array([0.0, np.pi/2, np.pi]), a=0.1, dUdt=1.0)   # front, side, back of the sphere
print(np.round(s["p"], 3), round(s["force"], 3), round(s["added_mass"], 3))   # +50, 0, −50 Pa; −2.094 N; 2.094 kg
print(ch04.gauge_absorbed_bracket(3.0, sign=-1.0), ch04.gauge_absorbed_bracket(3.0, sign=+1.0))   # 0.0 and 6.0: the printed + doubles B
print(ch04.viscous_irrotational_residual(u_cyl, np.array([1.5, 0.7]), mu=1e-3))   # −μ∇×ω ≈ 0: no net viscous force (N98)
""", explain=r"""
1. The column's level, speed, acceleration and the pressure difference that drives it (signs as in the function's
   docstring: h = h₀ cos ωt, U = dh/dt): the pressure difference is largest (±981 Pa) when the column stops.
2. The streamline form (4.82) balanced by ρL dU/dt.
3. The unsteady surface pressure on the sphere and the added-mass force.
4. The gauge sign test of (4.74): the minus sign absorbs B(t) = 3; the printed plus leaves 2B = 6.
5. Viscosity exerts no net force on an irrotational flow.
""")
nb.check_agree(r"""
t = sp.symbols('t')                                                  # time
phi = sp.Function('phi')(t); Bt = sp.Function('B')(t)                # the potential (its time part) and the bracket's value B(t)
rest = Bt - sp.diff(phi, t)                                          # everything in the bracket except ∂φ/∂t, so ∂φ/∂t + rest = B(t)
new_minus = sp.diff(phi - sp.Integral(Bt, t), t) + rest              # bracket after φ → φ − ∫B dt (our sign)
new_plus = sp.diff(phi + sp.Integral(Bt, t), t) + rest               # bracket after φ → φ + ∫B dt (the printed sign); sp.Integral is an unevaluated ∫
print(sp.simplify(new_minus), sp.simplify(new_plus))                 # 0 and 2*B(t)
assert sp.simplify(new_minus) == 0 and sp.simplify(new_plus - 2*Bt) == 0
assert np.isclose(ch04.gauge_absorbed_bracket(3.0, -1.0), 0.0)       # the library agrees
""")
nb.animation(r"""
L_, h0 = 1.0, 0.05                                                   # column length [m], initial displacement [m]
frames = 60 if not FAST else 30                                        # number of frames (halved when FAST)
tt = np.linspace(0, 2*np.pi/np.sqrt(2*9.81/L_), frames)             # one period [s]
col = ch04.u_tube_column(tt, L=L_, h0=h0)                            # h(t), U(t), dU/dt, p1 − p2
fig, (a, b) = plt.subplots(1, 2, figsize=(7.5, 3.6), gridspec_kw=dict(width_ratios=[1, 1.1]), dpi=80)  # the figure and its panels
a.plot([0, 0, 1, 1], [0.4, 0, 0, 0.4], color=COLORS["ink"], lw=3); a.plot([0.25, 0.25, 0.75, 0.75], [0.4, 0.25, 0.25, 0.4], color=COLORS["ink"], lw=3)  # draw the data
(left,) = a.plot([], [], color=COLORS["blue"], lw=16, solid_capstyle="butt")   # water in the left leg
(right,) = a.plot([], [], color=COLORS["blue"], lw=16, solid_capstyle="butt")  # water in the right leg
a.fill_between([0, 1], 0, 0.25, color=COLORS["blue"], alpha=0.6)    # the bottom of the U (always full)
a.set_xlim(-0.2, 1.2); a.set_ylim(0, 0.45); a.axis("off"); a.set_title("water in a U-tube", fontsize=9)  # the panel's message
bar = b.bar([0], [0], color=COLORS["amber"], width=0.5)[0]           # p1 − p2 = ρL dU/dt (the ∂φ/∂t part)
b.axhline(0, ls="--", color=COLORS["muted"], label="steady Bernoulli at U = 0: p1 − p2 = 0")        # a reference line
b.set_xlim(-1, 1); b.set_ylim(-1100, 1100); b.set_xticks([]); b.set_ylabel("p1 − p2 [Pa]"); b.legend(fontsize=7, loc="lower right")  # axis labels with units
txt = b.text(-0.95, 950, "", fontsize=8)                                        # a label on the plot


def update(i):                                                       # draw frame i
    h = col["h"][i]                                                  # left level up by h, right down by h [m]
    left.set_data([0.125, 0.125], [0.25, 0.32 + 0.6*h]); right.set_data([0.875, 0.875], [0.25, 0.32 - 0.6*h])   # levels (drawn ×0.6)
    bar.set_height(col["dp"][i])                                     # pressure difference that accelerates the column
    txt.set_text(f"U = {col['U'][i]:+.3f} m/s   dU/dt = {col['dUdt'][i]:+.3f} m/s²")                # move this artist to the new frame
    return left, right, bar  # the artists that changed in this frame


show_animation(animate(update, frames=frames, fig=fig, interval=50))                                # build the movie from update() and play it
""")
see_read_change(
    see="The amber bar is largest when the column stops at the top of its swing and zero when it moves fastest.",
    read=r"""The pressure difference tracks the acceleration, not the speed — the ∂φ/∂t term of
    $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=$ const *(4.75)*, invisible in any
    steady Bernoulli (whose prediction, the grey line, would be zero whenever U = 0).""",
    change="…a column twice as long: ω falls by √2, and for the same h₀ the bar keeps its size (ρL dU/dt = 2ρgh₀) — the weight of the height difference.")
nb.figure(r"""
th = np.linspace(0, np.pi, 181)                                      # angle from the front of the sphere [rad]
s = ch04.accelerating_sphere_pressure(th, a=0.1, dUdt=1.0)            # unsteady surface pressure, starting from rest [Pa]
fig = plt.figure(figsize=(9, 3.6))                                        # an empty figure
a = fig.add_subplot(1, 2, 1)                                        # a panel
a.plot(np.rad2deg(th), s["p"], color=COLORS["amber"], lw=2.5, label="unsteady part ½ρa cosθ dU/dt")  # draw the data
a.plot(np.rad2deg(th), 0*th, "--", color=COLORS["muted"], label="steady part (U = 0 at the start)")  # draw the data
a.set_xlabel("θ from the front [°]"); a.set_ylabel("p − p∞ [Pa]"); a.legend(fontsize=8)             # axis labels with units
a.set_title(f"net force {s['force']:.2f} N: added mass {s['added_mass']:.2f} kg", fontsize=10)      # the panel's message
b = fig.add_subplot(1, 2, 2, projection="polar")                    # pressure arrows round the sphere
tq = np.linspace(0, 2*np.pi, 24, endpoint=False)                                        # evenly spaced values
pq = ch04.accelerating_sphere_pressure(np.abs(np.angle(np.exp(1j*tq))), a=0.1, dUdt=1.0)["p"]   # same formula round the circle
b.plot(np.linspace(0, 2*np.pi, 200), np.ones(200), color=COLORS["ink"])                             # draw the data
for t_, p_ in zip(tq, pq):                                           # inward arrow where p > 0, outward where p < 0
    b.annotate("", (t_, 1.0), (t_, 1.0 + 0.012*p_), arrowprops=dict(arrowstyle="->", color=COLORS["amber"]))  # a labelled arrow
b.set_ylim(0, 1.8); b.set_yticks([]); b.set_title("pushed in front (θ = 0), pulled behind", fontsize=9)  # the panel's message
savefig(fig, "ch04", "sphere_added_mass"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="A cosine: pushed in front, pulled behind.",
    read="Summed over the surface it resists the acceleration like extra mass: half the displaced water.",
    change="…the sphere moving at steady speed: the amber curve vanishes and the steady Bernoulli pressure 1 − (9/4)sin²θ times ½ρU² remains, symmetric front to back — no net force (d'Alembert's paradox, Ch. 6).")
nb.explainer("which_bernoulli", heading="Bernoulli is constant along what, exactly?", why=r"""
**Why interactive rather than a static figure:** whether a Bernoulli form holds depends on the hypotheses you grant and on where you put the two probes; toggling them yourself shows which constant survives, which a fixed table cannot show.

Six scenarios — pitot tube, tank orifice, Rankine vortex, cylinder, U-tube, hot-gas nozzle — each with its streamlines,
two probes, and B plotted along and across streamlines as stacked bars (½u², p/ρ, gz, ∂φ/∂t). Toggle the hypotheses and
the status picks the valid equation: $\tfrac12U^2+gz+p/\rho$ (4.19), the Bernoulli function (4.71)/(4.72), the unsteady
form $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int dp/\rho+gz$ (4.75) or the energy form
$h+\tfrac12\lvert\mathbf u\rvert^2+gz$ (4.78).
""", tries=[
    "Rankine vortex: B is flat along each circle but not across the core — the status says '(4.71) along streamlines only'.",
    "Cylinder: B is flat everywhere ((4.72)); drag the probes anywhere.",
    "U-tube: stop the column (U = 0) — the amber ∂φ/∂t bar carries the pressure difference ((4.75)).",
    r"Hot nozzle: B of (4.71) with ρ = const fails; the energy form $h+\tfrac12\lvert\mathbf u\rvert^2+gz$ (4.78) stays flat.",
])
whatif(r"""
…the density varied a little from place to place, as in a lake warmed at the top? Then gravity no longer cancels into the
pressure, and a small density difference, multiplied by the large g, drives motion — the Boussinesq approximation (C13).
""")

# ---- R10 + C13 ----------------------------------------------------------------------------------------------------
nb.recap("R10", "Hydrostatic balance", r"""
A fluid at rest in gravity has $0=-\nabla p_s+\rho_s\mathbf g$ — pressure increases downward at the rate ρ_s g (Ch. 1,
$dp/dz=-\rho g$ *(1.8)*). ρ_s(z) and p_s(z) may vary with height (a stratified reference state).
""", where="Ch. 1 §1.7")
nb.code(r"""
z = np.linspace(0, -100, 5)                                          # depths 0 … 100 m in seawater [m]
print(np.round(statics.integrate_hydrostatic(z, lambda z, p: 1025.0, 1.013e5, g=9.81)))   # p(z) with g = 9.81 as elsewhere: +1 atm about every 10 m [Pa]
""")
core("C13", "Boussinesq: density matters only where it meets gravity (4.86)", """
A lake is 2 °C warmer at the top than at the bottom — its density changes by less than one part in a thousand. Yet that
tiny difference stops the wind from mixing the lake. Why does such a small density change matter so much, and when may
we ignore it?
""", eqs=("4.86",))
problem("""
In the ocean and the atmosphere, density varies by a percent or less, but those variations drive the overturning
circulation, sea breezes, convection and internal waves. Solving the full compressible equations for them is wasteful
(and fills the solution with sound waves). Boussinesq's idea keeps the density variation in exactly one place — where it
is multiplied by g — and treats the fluid as incompressible everywhere else. It is the equation set of Ch. 7's internal
waves, Ch. 11's convection and Ch. 13's ocean and atmosphere.
""")
idea("""
ρ = ρ₀ + ρ'   with  ρ'/ρ₀ ~ αδT ~ 10⁻³
inertia   ρ Du/Dt  ≈ ρ₀ Du/Dt           (error 10⁻³: drop ρ')
gravity   ρ g      = ρ_s g + ρ' g         (ρ_s g balanced by the hydrostatic p_s; ρ' g KEPT: g ≫ Du/Dt)
continuity  ∇·u ≈ 0;   heat  DT/Dt = κ∇²T  (pressure work turns C_v into C_p)
""")
P("P130", "order-of-magnitude scaling", r"""
To compare terms of an equation without solving it, replace each quantity by its typical size: a velocity varying by U
over a distance L has ∂u/∂x ~ U/L and ∂²u/∂x² ~ U/L²; a temperature varying by δT has ∂T/∂x ~ δT/L. Ratios of such
estimates tell which terms can be dropped (≪ 1) — the method of C13 and C15.
""", code="""
U, L, nu = 0.1, 1.0, 1e-6                         # a lake current [m/s], its scale [m], water's ν [m²/s]
advective = U**2/L; viscous = nu*U/L**2           # (u·∇)u ~ U²/L,  ν∇²u ~ νU/L²
print(advective, viscous, advective/viscous)      # 0.01, 1e-07, 1e5 = Re
""")
P("P131", "reduced gravity and buoyancy", r"""
A parcel lighter than its surroundings by Δρ feels a net upward force per unit mass g′ = gΔρ/ρ₀ — gravity *reduced* by the
density contrast. The buoyancy b = −gρ′/ρ₀ (upward positive) is the same idea as a field. For a density contrast of
10⁻³, g′ ≈ 0.01 m/s²: small, but it acts on the whole mass and for as long as the contrast lasts. N² = ∂b/∂z is the
buoyancy frequency of Ch. 1.
""", code="""
rho0, d_rho = 1000.0, -0.4                          # a parcel 0.4 kg/m³ lighter than 1000 kg/m³ water
print(ch04.buoyancy(d_rho, rho0), 9.81*0.4/1000)    # +0.003924 m/s² upward, both ways
""")
note("N106", "Subtract the hydrostatic state", r"""
(R10) from $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(4.39b)*; with p′ = p − p_s and
ρ′ = ρ − ρ_s:

> ⚠️ **Primes again**: here p′ and ρ′ are departures from the hydrostatic state — not the rotating frame of C09, not the
> dummy variable of (4.67).
""", equation=EQ["4.84"], ref="4.84")
note("N107", "Constant density", r"""
ρ′ = 0 and gravity disappears from the equation of motion. Gravity returns as soon as there is a free surface, an
interface or any density variation. Demo below: the Poiseuille velocity in a channel with gravity along −y solves
(4.39b) with p including −ρgy, and (4.85) with p′ alone.
""", equation=EQ["4.85"], ref="4.85")
nb.code(r"""
u_g, p_g = ch04.exact_field("poiseuille", G=100.0, h=1e-3, mu=1e-3, g=9.81)   # p = −Gx − ρgy: hydrostatic part included
X0 = np.array([0.0, 2.5e-4, 0.0])                                    # a point in the gap
full = ch04.ns_incompressible_terms(u_g, p_g, X0, 0.0, rho=1000.0, mu=1e-3, g=np.array([0, -9.81, 0]), per="volume")   # (4.39b) with g
p_prime = lambda x, t: p_g(x, t) + 1000.0*9.81*x[1]                  # p' = p − p_s: remove the hydrostatic part ρgy
pert = ch04.ns_incompressible_terms(u_g, p_prime, X0, 0.0, rho=1000.0, mu=1e-3, g=np.zeros(3), per="volume")   # (4.85): no g
print(full.pressure, full.gravity, full.residual)                    # pressure (100, 9810, 0) balances gravity (0, −9810, 0)
print(pert.pressure, pert.gravity, pert.residual)                    # (100, 0, 0), no gravity at all, residual ≈ 0
""")
note("N108", "When Boussinesq holds", r"""
(stated; the scaling chain is the book's): low Mach number and no sound; a vertical scale L ≪ c²/g (the height over which
hydrostatic compression changes the density — ≈ 11.8 km for air at 288 K); and small temperature differences, αδT ≪ 1,
because then the ratio below is small and continuity reduces to ∇·u = 0. The code below computes the three small numbers
for five situations (`ch04.boussinesq_validity`): only the deep atmosphere fails.
""", equation=r"\frac{(1/\rho)(D\rho/Dt)}{\nabla\cdot\mathbf u}\sim\frac{(U/\rho)(\delta\rho/L)}{U/L}=\frac{\delta\rho}{\rho}=\alpha\,\delta T\ll1")
remind([("perfect gas p = ρRT, R = C_p − C_v, α = 1/T", r"for a perfect gas $p=\rho RT$, $R=C_p-C_v$ and the expansion coefficient $\alpha=-\frac1\rho(\partial\rho/\partial T)_p=1/T$ (Ch. 1 P33).")])
D("D27", ref="4.86")
note("C13", "The Boussinesq momentum equation", r"""
C13's momentum equation, restated. ρ₀ is a constant reference density. In buoyant convection the term ρ′g/ρ₀ is as large
as ∂w/∂t or ν∇²w.
""", equation=EQ["4.86"], ref="4.86")
note("N109", "The energy equation (4.60) in vector form", "— the start of D28:", equation=EQ["4.87"], ref="4.87")
D("D28", ref="4.89")
note("N110", "Why C_p and not C_v", r"""
(D28 steps 2–6): although ∇·u ≈ 0 in continuity, the pressure work p∇·u is *not* negligible in the energy equation,
because p is large: for a perfect gas $-p\nabla\cdot\mathbf u\cong-p\alpha\frac{DT}{Dt}=-\rho(C_p-C_v)\frac{DT}{Dt}$. It moves
to the left and turns ρC_v DT/Dt into ρC_p DT/Dt.

> ⚠️ For liquids the argument is different: p∇·u is small, but C_p ≈ C_v anyway, so the ocean uses the same equation with
> water's C_p.
""")
nb.code(r"""
rho, T, R, Cp, DTDt = sp.symbols('rho T R C_p DT_Dt', positive=True)  # density, temperature, gas constant, C_p, the rate DT/Dt
Cv = Cp - R                                                          # perfect gas: R = C_p − C_v
p = rho*R*T                                                          # perfect-gas law p = ρRT
alpha = 1/T                                                          # its expansion coefficient α = −(1/ρ)(∂ρ/∂T)_p = 1/T
div_u = -(1/rho)*(-rho*alpha*DTDt)                                   # ∇·u = −(1/ρ)Dρ/Dt with Dρ/Dt = −ρα DT/Dt (D28 steps 2–3)
lhs = rho*Cv*DTDt + p*div_u                                          # ρ De/Dt + p∇·u with e = C_v T (step 6)
print(sp.simplify(lhs - rho*Cp*DTDt))                                # 0: the pressure work turns ρC_v into ρC_p
""", explain="The D28 check in sympy: with e = C_vT, p = ρRT and α = 1/T, ρC_v DT/Dt + p∇·u equals ρC_p DT/Dt exactly.")
note("N111", "Heat equation with C_p", "", equation=EQ["4.88"], ref="4.88")
note("N112", "Viscous heating is negligible under Boussinesq conditions", r"""
Water, U = 0.1 m/s, L = 1 m, δT = 1 K: 10⁻⁷/4186 ≈ 2×10⁻¹¹ (the book's "typical" value is 10⁻⁷ — tiny either way). This
is the low-Eckert-number condition of C15.
""", equation=r"\frac{\rho\varepsilon}{\rho C_p(DT/Dt)}\sim\frac{2\mu S_{ij}S_{ij}}{\rho C_pu_i(\partial T/\partial x_i)}\sim\frac{\mu U^2/L^2}{\rho C_pU(\delta T/L)}=\frac{\nu U}{C_p\,\delta T\,L}")
note("N113", "The Boussinesq heat equation", r"""
(D28's result), κ = k/ρC_p. Exact test solution (the animation below): a warm Gaussian blob carried by a uniform current U
and spreading, $T'=T_a(\sigma_0/\sigma)^2\exp(-\lvert\mathbf x-\mathbf Ut\rvert^2/2\sigma^2)$ in 2-D with
$\sigma^2=\sigma_0^2+2\kappa t$.
""", equation=EQ["4.89"], ref="4.89")
note("N114", "The Boussinesq set", r"""
(boxed): continuity $\nabla\cdot\mathbf u=0$ *(4.10)*; momentum (4.86) with $\mathbf g=-g\mathbf e_z$; heat (4.89); and a
linear equation of state (below) with ρ′ = −ρ₀α(T − T₀), so the buoyancy is b = −gρ′/ρ₀ = gα(T − T₀) and the
stratification $N^2=\partial b/\partial z=g\alpha\,dT/dz$ (Kundu's Γ ≡ dT/dz; in the meteorological convention
Γ_met = −dT/dz, N² = −gαΓ_met). For seawater, salinity enters too (`ch01.seawater_density_linear`). The Coriolis term is
omitted here — Ch. 13 adds it back.
""", equation=r"\rho=\rho_0\big[1-\alpha(T-T_0)\big]")
nb.worked_example("three Boussinesq numbers", r"""
One set of water numbers throughout: α = 2×10⁻⁴ K⁻¹, a warm layer δT = 10 K above cold water (this is the code's
`"thermocline"` row below; the question's lake, δT = 2 K, has αδT = 4×10⁻⁴, and the code's gentler `"lake"` row uses
α = 1.5×10⁻⁴, δT = 5 K).

1. **Thermocline** (water, α ≈ 2×10⁻⁴ K⁻¹, δT = 10 K): αδT = 2×10⁻³ ≪ 1; reduced gravity g′ = gαδT = 9.81 × 2×10⁻³ =
   0.0196 m/s² (the code prints the same).
2. **Air** at 300 K (perfect gas, α = 1/T): δT = 10 K → αδT = 0.033 — still small.
3. **Depth limit**: c²/g for air at 288 K = 340.3²/9.81 = 11.8 km — a flow 1 km deep is fine, the whole troposphere is not.
4. **Rise speed** of a 1 m warm patch in the thermocline: √(g′L) ≈ 0.14 m/s.
5. **Viscous heating** for U = 0.1 m/s, L = 1 m, δT = 1 K: νU/(C_pδT L) = 10⁻⁶ × 0.1/(4186 × 1 × 1) ≈ 2×10⁻¹¹ (the code's
   scenarios use their own U, L, δT and give even smaller ratios, 10⁻¹³–10⁻¹¹).
""")
remind([("np.where", "`np.where(cond, a, b)` picks a where cond is true, b elsewhere (Ch. 1 P46).")])
nb.code(r"""
for name in ("lake", "thermocline", "lab_tank", "abl", "deep_atmosphere"):   # five real situations
    p = ch04.boussinesq_scenario(name)                                # fluid, α, δT, L, U, ν, c_p, c
    v = ch04.boussinesq_validity(p["alpha"], p["dT"], p["L"], p["U"], c=p["c"], nu=p["nu"], cp=p["cp"])   # the three small numbers
    print(f"{name:16s} αδT={v['alpha_dT']:.1e}  L/(c²/g)={v['L_over_Hc']:.2e}  heating={v['heating_ratio']:.1e}"  # show the numbers
          f"  g'={v['g_prime']:.3g} m/s²  → {v['verdict']}")
print(ch04.boussinesq_density(303.15, 1000.0, 2e-4, 293.15), ch04.buoyancy(-2.0, 1000.0))   # 998.0 kg/m³; +0.0196 m/s² upward
Tb = lambda x, t: ch04.gaussian_blob_advection_diffusion(x, t, U=np.array([0.01, 0.0]), kappa=1e-4, sigma0=0.05)   # a warm blob T'(x, t)
ub = lambda x, t: np.array([0.01, 0.0])                              # the uniform current carrying it [m/s]
print(ch04.temperature_equation_residual(Tb, ub, np.array([0.02, 0.01]), 10.0, kappa=1e-4, h=1e-4, ht=1e-3))   # ≈ 0: it solves (4.89)
""", explain=r"""
1. The three small numbers for five real situations: αδT, the depth over the scale height c²/g, and the viscous-heating
   ratio; the deep atmosphere fails on the first two.
2. The linear equation of state and the buoyancy of a 2 kg/m³ lighter parcel.
3. The blob solves $\frac{DT}{Dt}=\kappa\nabla^2T$ *(4.89)* to stencil accuracy (residual ~3×10⁻⁸ against terms of
   3×10⁻³ K/s; κ = 10⁻⁴ m²/s is an eddy diffusivity, to see it spread in seconds).
""")
nb.check_agree(r"""
x0, t0, h, ht = np.array([0.02, 0.01]), 10.0, 1e-4, 1e-3             # point [m], time [s], steps
dTdt = (Tb(x0, t0 + ht) - Tb(x0, t0 - ht)) / (2*ht)                  # ∂T/∂t
grad = np.array([(Tb(x0 + e, t0) - Tb(x0 - e, t0))/(2*h) for e in np.eye(2)*h])   # ∇T
lap = sum((Tb(x0 + e, t0) - 2*Tb(x0, t0) + Tb(x0 - e, t0))/h**2 for e in np.eye(2)*h)   # ∇²T
adv = ub(x0, t0) @ grad                                              # u·∇T
print(dTdt + adv, 1e-4*lap)                                          # DT/Dt and κ∇²T: equal
assert abs(dTdt + adv - 1e-4*lap) < 1e-4*abs(dTdt)                   # (4.89) holds to stencil accuracy
res = ch04.heat_equation_terms(Tb, ub, x0, t0, rho=1000.0, cp=4186.0, k=1e-4*1000*4186, eps=0.0)["residual"]   # library, (4.88) form
assert abs(res) < 1e-4*abs(dTdt)
""")
nb.md("**What does the code above do?** It evaluates DT/Dt = ∂T/∂t + u·∇T and κ∇²T for the warm blob with our own stencils, checks they are equal (the heat equation holds), and that the library's residual is as small.")
nb.animation(r"""
n = 121 if not FAST else 81                                          # grid points per side
xs_, ys_ = np.linspace(-0.2, 0.6, n), np.linspace(-0.25, 0.25, (5*n)//8)   # a 2-D tank [m]
Xg, Yg = np.meshgrid(xs_, ys_)                                        # 2-D grids of the coordinates
frames = 60 if not FAST else 30                                        # number of frames (halved when FAST)
tt = np.linspace(0, 30, frames)                                      # 30 s [s]
kap, s0, Uc = 1e-4, 0.05, 0.01                                       # κ [m²/s], initial width [m], current [m/s]
fig, ax = plt.subplots(figsize=(6.5, 3.4), dpi=80)                                        # the figure and its panels
im = ax.imshow(ch04.gaussian_blob_advection_diffusion(np.array([Xg, Yg]), 0.0, U=np.array([Uc, 0.0]), kappa=kap, sigma0=s0),  # draw once; update() moves it
               extent=(xs_[0], xs_[-1], ys_[0], ys_[-1]), origin="lower", cmap="Reds", vmin=0, vmax=1)   # T'/T_a
(ring,) = ax.plot([], [], "--", color=COLORS["ink"], lw=1)           # the circle of radius σ(t)
ttl = ax.set_title(""); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")                              # axis labels with units
ph = np.linspace(0, 2*np.pi, 100)                                        # evenly spaced values


def update(i):                                                       # draw frame i
    t = tt[i]                                        # the time of this frame [s]
    im.set_data(ch04.gaussian_blob_advection_diffusion(np.array([Xg, Yg]), t, U=np.array([Uc, 0.0]), kappa=kap, sigma0=s0))  # move this artist to the new frame
    sig = np.sqrt(s0**2 + 2*kap*t)                                   # σ² = σ0² + 2κt
    ring.set_data(Uc*t + sig*np.cos(ph), sig*np.sin(ph))                                        # move this artist to the new frame
    ttl.set_text(f"t = {t:4.1f} s   peak T'/T_a = (σ0/σ)² = {(s0/sig)**2:.2f}")                     # move this artist to the new frame
    return im, ring                                        # the artists that changed in this frame


show_animation(animate(update, frames=frames, fig=fig, interval=60))                                # build the movie from update() and play it
""")
see_read_change(
    see="The blob drifts right and flattens; its dashed σ-circle grows; its total heat (area under the bump) stays constant.",
    read=r"""Advection moves it ($\mathbf u\cdot\nabla T$), diffusion spreads it ($\kappa\nabla^2T$); the peak falls as
    $1/(1+2\kappa t/\sigma_0^2)$ — both halves of $\frac{DT}{Dt}=\kappa\nabla^2T$ *(4.89)* at work.""",
    change="…κ halved: the circle grows half as fast; the drift is unchanged.")
nb.figure(r"""
names = ["lake", "thermocline", "lab_tank", "abl", "deep_atmosphere"]  # the five scenarios
vals = []                                        # an empty list to fill
for name in names:                                                    # the three small numbers of each
    p = ch04.boussinesq_scenario(name)                                        # call the tested chapter function
    v = ch04.boussinesq_validity(p["alpha"], p["dT"], p["L"], p["U"], c=p["c"], nu=p["nu"], cp=p["cp"])  # call the tested chapter function
    vals.append([v["alpha_dT"], v["L_over_Hc"], v["heating_ratio"]])  # keep the three small numbers
vals = np.array(vals)  # 5 scenarios × 3 numbers
fig, ax = plt.subplots(figsize=(8, 3.8))                                        # the figure and its panels
labels = ["αδT", "L/(c²/g)", "νU/(C_p δT L)"]                                        # one colour / label per series
cols = [COLORS["blue"], COLORS["accent"], COLORS["rose"]]                                        # one colour / label per series
for k in range(3):                                                    # three bars per scenario, log scale
    bars = ax.barh(np.arange(5) + (k - 1)*0.25, vals[:, k], 0.25, color=cols[k], label=labels[k], log=True)  # one row of bars per number (log axis)
    for j, b_ in enumerate(bars):  # mark the bars beyond the threshold
        if vals[j, k] > 0.1:  # neglected effect not small
            b_.set_hatch("//"); b_.set_edgecolor(COLORS["ink"])       # crossing the threshold: hatched, colour kept
ax.axvline(0.1, ls="--", color=COLORS["ink"]); ax.text(0.12, 4.4, "threshold 0.1", fontsize=8)      # a label on the plot
ax.set_yticks(np.arange(5), names); ax.set_xlabel("size of the neglected effect [–]"); ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))  # labels; legend outside the bars
ax.set_title("Boussinesq is excellent except for flows as deep as the atmosphere")                  # the panel's message
savefig(fig, "ch04", "boussinesq_validity"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Most bars far to the left of the dashed line; the deep atmosphere's two bars crossing it, hatched.",
    read="""Boussinesq is excellent for lakes, oceans and boundary layers; it fails for flows as deep as the atmosphere's
    scale height and with large temperature contrasts (then use the anelastic or full compressible equations).""",
    change="…a 50 K temperature difference in air: αδT ≈ 0.17 — marginal.")
nb.explainer("boussinesq_buoyancy", heading="Density hardly changes — so why does it drive the flow?", why=r"""
**Why interactive rather than a static figure:** the approximation is a judgement about three small numbers; moving the fluid, the temperature contrast and the depth shows when each crosses its threshold while the blob keeps (or loses) its buoyancy.

Choose a fluid, a temperature contrast and a depth: three small numbers are compared with their thresholds while a warm
blob rises and spreads under exactly the terms Boussinesq keeps,
$\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ (4.86) and
$\frac{DT}{Dt}=\kappa\nabla^2T$ (4.89). The rise uses a one-line parcel model (our extension, labelled in the explainer):
$dw/dt=g'-w/\tau_d$. Switch the buoyancy term off and the blob stops rising, though ρ changed by only 0.2 %.
""", tries=[
    "Lake preset: read αδT and g′ in Explain.",
    "Toggle 'keep buoyancy' off: the blob only diffuses.",
    "Deep-atmosphere preset: the status turns amber — L is comparable to c²/g.",
    "Derivation tab D27, the 'why keep ρ′g' step: the grey dropped bar is a thousand times smaller than the blue buoyancy bar.",
])
whatif(r"""
…we asked what happens at the edges — the walls, the free surface, the interface between two fluids? The field equations
are complete; to solve them we need conditions on the boundaries (C14).
""")
# =====================================================================================================================
# A.10 §4.10 Boundary Conditions — R11, C14
# =====================================================================================================================
nb.section("4.10", "Boundary Conditions", intro="""
**What is this section about?** The field equations need conditions at the edges of the flow. Conservation laws applied
to a thin "pillbox" straddling an interface give continuity of mass flux, traction and heat flux; experience adds no-slip
and no temperature jump; a moving surface made of fluid obeys the kinematic condition Dη/Dt = 0; and surface tension adds
a pressure jump across a curved interface.
""")
nb.recap("R11", "Surface tension", r"""
Molecules at a liquid surface are pulled inward more than outward; the surface behaves like a stretched membrane with
tension σ [N/m], the energy per unit area of surface. Water at 20 °C has σ ≈ 0.0728 N/m (Ch. 1 §1.6,
`ch01.surface_tension_water`). Here we *derive* the pressure jump it causes.
""", where="Ch. 1 §1.6")
core("C14", "The surface is made of fluid: the kinematic condition (4.91) and other edge rules", """
A wave runs across the sea. The water on the surface stays on the surface — it neither sinks into the sea nor flies off.
What does that sentence say as an equation?
""", eqs=("4.91",))
problem("""
Every solution of Navier–Stokes is fixed by what happens at its edges: the wall of a pipe, the sea surface, the boundary
between warm and cold layers of the ocean, the skin of a raindrop. Some conditions follow from the conservation laws
(what crosses an interface must be continuous), some from experiment (fluid sticks to walls), and one from geometry (a
free surface moves with the fluid on it). Ch. 7 builds water waves on the last one; Ch. 13's layered ocean models on the
same condition at each layer interface.
""")
idea("""
pillbox across an interface, height l → 0:  volume and side terms vanish  →  flux in = flux out:
      ρu·n continuous, traction n_iτ_ij continuous, k ∂T/∂n continuous
experiment at a solid wall:  no slip (u·t = 0 relative to the wall), no temperature jump (T₁ = T₂)
a surface η(x, t) = 0 made of fluid:  Dη/Dt = ∂η/∂t + u·∇η = 0   (4.91)  — normal velocities match
curved interface with tension σ:  Δp = σ(1/R₁ + 1/R₂)
""")
remind([("boundary conditions", "the values (or fluxes) a field must take on the edges of its domain; with them a differential equation has one solution (Ch. 1 P20).")])
note("N115", "What must be specified", r"""
the velocity on every bounding surface; for an external flow, the velocity and thermodynamic state on a distant closed
surface. The solvers of Ch. 8–10 show it in action.
""")
note("N116", "Interface conditions from a pillbox", r"""
(stated; the book's Fig. 4.18 argument): apply the budgets (4.5), (4.17) and (4.48) to a flat cylinder straddling the
interface and let its height l → 0 — the volume integrals and the side-wall fluxes (∝ l) vanish, leaving the two end faces
(interface at rest; traction continuous if surface tension is neglected — with it, the Laplace jump below is added).
**Two worked cases** (code): (a) a composite wall, k₁ = 1, k₂ = 0.1 W/(m K), each 0.1 m thick, 30 °C and 20 °C outside:
the heat flux 10/(0.1/1 + 0.1/0.1) = 9.09 W/m² is the same in both layers, the interface is at 29.09 °C and the
temperature slopes differ by the factor k₁/k₂ = 10 — the profile has a kink; (b) two-fluid Couette, water (μ₁ = 10⁻³ Pa s,
1 mm) under oil (μ₂ = 0.1 Pa s, 1 mm), top plate at 1 m/s: the shear stress τ = U/(h₁/μ₁ + h₂/μ₂) = 0.990 Pa is
continuous, the velocity slopes differ by μ₂/μ₁ = 100, the interface moves at 0.990 m/s.
""", equation=r"\rho_1\mathbf u_1\cdot\mathbf n=\rho_2\mathbf u_2\cdot\mathbf n,\qquad n_i\tau^{(1)}_{ij}=n_i\tau^{(2)}_{ij},\qquad k_1\frac{\partial T_1}{\partial n}=k_2\frac{\partial T_2}{\partial n}")
nb.code(r"""
print(ch04.two_layer_conduction(np.array([0.05, 0.15]), 1.0, 0.1, 0.1, 0.1, 303.15, 293.15))   # q = 9.09 W/m² in both layers, interface 302.24 K
print(ch04.two_fluid_couette(np.array([5e-4, 1.5e-3]), 1e-3, 0.1, 1e-3, 1e-3, 1.0))           # τ = 0.990 Pa in both, interface at 0.990 m/s
print(ch04.pillbox_limit(9.09, 9.09, side_rate=2.0, volume_rate=5.0, l=np.array([1e-1, 1e-2, 1e-3]))["residual"])   # ∝ l → 0
""", explain="""
1. The composite wall: one heat flux through both layers, a kinked temperature profile.
2. Two-fluid Couette: one shear stress in both fluids, a kinked velocity profile.
3. The pillbox balance for three heights: the side and volume terms shrink with l, leaving flux in = flux out.
""")
note("N117", "No slip and no temperature jump", r"""
at a solid wall — not consequences of conservation laws, but of experiment. Exceptions: superfluid helium (no viscosity),
super-hydrophobic textured surfaces (apparent slip), and rarefied gases whose mean free path is not small compared with
the flow (Knudsen number Kn ≳ 0.01, Ch. 1). The figure below adds a Couette profile with a Navier slip length b (our
extension, labelled): u = U(y + b)/(h + 2b) — invisible for b ≪ h.
""", equation=r"\mathbf u_1\cdot\mathbf t=0\ \ (\text{relative to the wall}),\qquad T_1=T_2")
remind([("unit normal ∇η/‖∇η‖", r"the gradient of a level-set function is perpendicular to its level surfaces; dividing by its length gives the unit normal (Ch. 2 P75).")])
P("P132", "moving level set and its normal speed", r"""
A surface can be described as the set of points where a function is zero, η(x, t) = 0 (a level set, Ch. 2 P75). Its unit
normal is $\mathbf n=\nabla\eta/\lvert\nabla\eta\rvert$ (pointing toward increasing η). If it moves, a point riding on it
keeps η = 0, and the chain rule gives its speed along n: $-(\partial\eta/\partial t)/\lvert\nabla\eta\rvert$. Only this
normal speed is defined — sliding along the surface does not change the surface.
""", code="""
V = 1.5                                                        # a wall moving at 1.5 m/s along +x
eta = lambda x, t: x[0] - V*t                                  # the plane x = Vt as a level set
print(ch04.surface_normal_speed(eta, np.array([2.0, 0.3, 0.0]), 1.0))   # 1.5: it moves at V along n = +x
""")
note("N118", "An observer riding on the moving surface", r"at velocity $\mathbf u_s$ sees η stay zero:",
     equation=r"d\eta/dt=\partial\eta/\partial t+(\mathbf u_s\cdot\nabla)\eta=0\quad\text{on}\quad\eta=0", ref="4.90")
D("D29", ref="4.91")
note("C14", "The kinematic boundary condition — no fluid crosses the surface", r"""
C14's result, restated. A free surface, an interface between immiscible fluids and a solid wall all obey it.
""", equation=r"\partial\eta/\partial t+(\mathbf u\cdot\nabla)\eta\equiv D\eta/Dt=0\quad\text{on}\quad\eta=0", ref="4.91")
note("N119", "If mass does cross", r"(evaporation, a shock), the normal velocity of the fluid relative to the surface is",
     equation=r"(u_{rel})_n=\mathbf u\cdot\mathbf n-\mathbf u_s\cdot\mathbf n=\big(\mathbf u\cdot\nabla\eta+\partial\eta/\partial t\big)/\lvert\nabla\eta\rvert=(1/\lvert\nabla\eta\rvert)D\eta/Dt", ref="4.92")
note("N120", "and the mass flux per unit area through it", "Zero flux is again Dη/Dt = 0. Moving shocks: Ch. 15. (Code: `ch04.interface_mass_flux`, below.)",
     equation=r"(\rho/\lvert\nabla\eta\rvert)D\eta/Dt\quad\text{on}\quad\eta=0", ref="4.93")
remind([("Gibbs free energy (as the model for f)", r"g = h − Ts traded an unmeasurable variable for a measurable one; f = e − Ts does the same at fixed T and v (Ch. 1 P50).")])
note("N121", "Surface tension as free energy", r"""
At constant temperature the useful quantity is the Helmholtz free energy per unit mass (a Legendre device, like Ch. 1's
Gibbs g = h − Ts):
""", equation=EQ["4.94"], ref="4.94")
note("N122", "Its differential", r"""
At constant T and reversibly ($T\,ds=de+p\,dv$ *(1.18)*): df = −p dv — the work done on the system. (The sympy line
below: for a perfect gas at fixed T, e is constant and s = C_v ln T + R ln v + const, so df/dv = −RT/v = −p.)
""", equation=EQ["4.95"], ref="4.95")
nb.code(r"""
v, T, R, Cv = sp.symbols('v T R C_v', positive=True)                  # specific volume, temperature, gas constants
e = Cv*T                                                             # internal energy of a perfect gas: fixed at fixed T
s = Cv*sp.log(T) + R*sp.log(v)                                       # its entropy (up to a constant)
f = e - T*s                                                          # Helmholtz free energy (4.94)
print(sp.simplify(sp.diff(f, v)))                                    # −R*T/v = −p: df = −p dv at constant T
""")
note("N123", "Free energy of two fluids and their interface", r"""
of area A: σ is free energy per unit area; with σ > 0 (immiscible fluids) the system lowers F by shrinking A — the
interface contracts as far as the constraints allow. Number (code): the surface area of a spheroid of fixed volume
(`ch04.spheroid_area`) is smallest at aspect ratio 1 — a sphere.
""", equation=EQ["4.96"], ref="4.96")
nb.code(r"""
for asp in (0.5, 0.8, 1.0, 1.25, 2.0):                               # oblate … sphere … prolate, volume 1 m³
    print(asp, round(ch04.spheroid_area(1.0, asp), 4))               # smallest area (4.836 m²) at aspect 1
""")
note("N124", "Marangoni flows", r"""
Where σ varies along an interface (with temperature or a surfactant) the surface pulls toward high σ and drags fluid with
it — *Marangoni* flows (tears of wine). Not treated in the book; see Ch. 16 and interfacial-flow texts.
""")
note("N125", "The Laplace jump derived", r"""
(stated with numeric checks; the book writes the construction): cut a small cap $z=x^2/2R_1+y^2/2R_2$ (R₁, R₂ its two
principal radii of curvature, Ch. 1 §1.6) at height ζ; its rim is an ellipse with semi-axes $\sqrt{2R_1\zeta}$,
$\sqrt{2R_2\zeta}$. The pressure difference Δp pushes on it with the force below.
""", equation=EQ["4.97"], ref="4.97")
note("N126", "Surface tension on the rim", r"""
The surface tension pulls along the rim ($\sigma\oint_C\mathbf t\times\mathbf n\,ds$) with a net upward force, for a
small cap, as below.

> ⚠️ **Book slip:** the curve C is printed as ζ = x²/2R₁ − y²/2R₂; it must be **+** (the cap equation); the book's later
> lines use +.
""", equation=EQ["4.98"], ref="4.98")
note("N127", "Both forces are ∝ ζ; they balance when", r"""
— Laplace's law (1.5), which Ch. 1 stated, now derived. Higher pressure on the concave side. Number: a water drop of
radius 1 mm (R₁ = R₂ = 1 mm, σ = 0.0728 N/m): Δp = 146 Pa. Capillary waves: Ch. 7.
""", equation=EQ["1.5"], ref="1.5")
note("N128", "Capillary length", r"""
— where gravity and surface tension are equally strong. Water at 20 °C: √(0.0728/(998 × 9.81)) = 2.73 mm. Smaller than
ℓ_c, surface tension wins (a drop is round); larger, gravity wins (a puddle is flat). Its square ratio is the Bond number
(C15).
""", equation=r"\ell_c=\big(\sigma/\rho g\big)^{1/2}")
note("N129", "Ex. 4.7 — the meniscus at a vertical wall", r"""
(stated; the closed form is checked against an ODE solution instead of re-derived): a liquid meets a wall at contact
angle θ; the surface height ζ(x) obeys the first relation below, and the second is the height at the wall. Water on clean
glass (θ = 0): h = √2 ℓ_c = 3.86 mm.

> ⚠️ **Book slips in Ex. 4.7:** the separated equation drops a minus sign (the slope is negative), and "evaluated when
> η = h/δ" should read γ = h/δ — here γ = ζ/δ, yet another γ; the final closed form is right (verified by
> differentiation in the tests).
""", equation=r"\Big(\frac{\rho g}{2\sigma}\Big)\zeta^2+\big(1+\zeta'^2\big)^{-1/2}=1,\qquad h^2=\frac{2\sigma}{\rho g}(1-\sin\theta)")
nb.worked_example("three edges", r"""
1. **A moving wall** η = x − Vt, V = 1 m/s, with fluid at u = (1, 0.3, 0) m/s on it: ∂η/∂t = −1, ∇η = (1, 0, 0), so
   Dη/Dt = −1 + 1 × 1 = 0 ✓ — the fluid's normal speed (1 m/s) equals the wall's; the tangential 0.3 m/s is not
   constrained by $D\eta/Dt=0$ *(4.91)* (no-slip is a separate condition).
2. **A small water wave** η = a cos(kx − ωt) with ka = 0.1: the linearised condition ∂η/∂t = w at z = 0 holds exactly for
   the linear solution; the full condition (4.91) at the true surface leaves a residual of order (ka)² ≈ 0.01 of the
   phase speed — a quarter of that at ka = 0.05.
3. **A drop**: R = 1 mm, σ = 0.0728 N/m → Δp = 2σ/R = 146 Pa; ℓ_c = 2.73 mm.
""")
nb.code(r"""
u_w = lambda x, t: np.array([1.0, 0.3, 0.0])                         # fluid velocity at the wall [m/s]
print(ch04.kinematic_bc_residual("moving_wall", u_w, np.array([0.5, 0.0, 0.0]), 0.5, V=1.0),        # show the numbers
      ch04.surface_normal_speed("moving_wall", np.array([0.5, 0, 0]), 0.5, V=1.0))   # Dη/Dt = 0; the wall moves at 1 m/s
u_evap = lambda x, t: np.array([1.2, 0.3, 0.0])                      # fluid crossing the wall at 0.2 m/s relative to it
print(ch04.relative_normal_velocity("moving_wall", u_evap, np.array([0.5, 0, 0]), 0.5, V=1.0),      # show the numbers
      ch04.interface_mass_flux("moving_wall", u_evap, 1000.0, np.array([0.5, 0, 0]), 0.5, V=1.0))   # (4.92): 0.2 m/s, (4.93): 200 kg/(m² s)
for ka in (0.2, 0.1, 0.05):                                          # three wave steepnesses
    k = 1.0; a = ka/k; c = np.sqrt(9.81/k)                           # wavenumber [1/m], amplitude [m], phase speed [m/s]
    res = max(abs(ch04.wave_kinematic_residual(x, 0.0, a, k)) for x in np.linspace(0, 2*np.pi, 64))   # full (4.91) at the true surface
    print(ka, round(res/c, 5), ch04.wave_kinematic_residual(0.3, 0.0, a, k, full=False))   # ∝ (ka)²; the linearised one is 0
print(ch04.laplace_jump_from_balance(0.0728, 1e-3, 1e-3), ch01.laplace_pressure_jump(0.0728, 1e-3, 1e-3))   # 145.6 Pa both ways
print(ch04.cap_pressure_force(2*0.0728/1e-3, 1e-3, 1e-3, 1e-6), ch04.cap_surface_tension_force(0.0728, 1e-3, 1e-3, 1e-6))   # ≈ ∓9.2e-7 N with the exact Δp = 2σ/R
print(ch04.capillary_length(0.0728, 998.0), ch04.meniscus_height(0.0, 0.0728, 998.0))   # 2.727 mm and 3.856 mm
""", explain=r"""
1. The wall condition: Dη/Dt = 0 and the wall's normal speed; then a fluid that crosses the wall — (4.92) gives its
   relative normal speed and (4.93) its mass flux.
2. The full condition's residual for the linear wave shrinks like (ka)², while the linearised one is exact — the error made
   by linearising (Ch. 7).
3. Laplace from the cap balance equals Ch. 1's formula $\Delta p=\sigma(1/R_1+1/R_2)$ *(1.5)*.
4. The two cap forces nearly cancel (with the exact Δp = 2σ/R their small difference is the finite-ζ, higher-order part).
5. Capillary length and the meniscus height at θ = 0.
""")
nb.check_agree(r"""
a, k = 0.1, 1.0; w = np.sqrt(9.81*k)                                 # amplitude [m], wavenumber [1/m], frequency [rad/s]
eta = lambda x, t: x[1] - a*np.cos(k*x[0] - w*t)                     # the surface as the level set z − ζ(x, t) = 0
uw = lambda x, t: np.array([ch04.linear_wave_surface(x[0], x[1], t, a, k)["u"],  # the linear wave's velocity (u, w)
                            ch04.linear_wave_surface(x[0], x[1], t, a, k)["w"]])   # the linear wave's velocity (u, w)
p = np.array([0.3, a*np.cos(0.3)])                                   # a point ON the surface at t = 0
d = 1e-6                                                             # finite-difference step
deta_dt = (eta(p, d) - eta(p, -d)) / (2*d)                          # ∂η/∂t
grad = np.array([(eta(p + e, 0) - eta(p - e, 0))/(2*d) for e in np.eye(2)*d])   # ∇η
mine = deta_dt + uw(p, 0.0) @ grad                                   # Dη/Dt = ∂η/∂t + u·∇η by hand
print(mine, ch04.kinematic_bc_residual(eta, uw, p, 0.0))             # the same small O((ka)²) residual
assert np.isclose(mine, ch04.kinematic_bc_residual(eta, uw, p, 0.0), rtol=1e-5)
""")
nb.animation(r"""
a, k = 0.15, 1.0; w = np.sqrt(9.81*k); lam = 2*np.pi/k              # ka = 0.15 (exaggerated), deep water
uw2 = lambda x, t: np.array([ch04.linear_wave_surface(x[0], x[1], t, a, k)["u"],
                             ch04.linear_wave_surface(x[0], x[1], t, a, k)["w"]])   # velocity field of the wave [m/s]
frames = 60 if not FAST else 30                                        # number of frames (halved when FAST)
tt = np.linspace(0, 2*np.pi/w, frames)                               # one period [s]
x0s = np.linspace(0, 2*lam, 13)[:-1]                                 # particles along two wavelengths
top = [ch03.pathline(uw2, [x0, a*np.cos(k*x0)], 0.0, tt) for x0 in x0s]   # start ON the surface; follow them (solve_ivp)
deep = [ch03.pathline(uw2, [x0, -0.3*lam], 0.0, tt) for x0 in x0s]  # a row 0.3 λ below
xg = np.linspace(0, 2*lam, 300)                                        # evenly spaced values
fig, ax = plt.subplots(figsize=(7, 3.2), dpi=80)                                        # the figure and its panels
(surf,) = ax.plot(xg, a*np.cos(k*xg), color=COLORS["blue"], lw=2)                                   # draw once; update() moves it
dt_ = ax.scatter([p[0, 0] for p in top], [p[1, 0] for p in top], color=COLORS["orange"], s=18, zorder=3)  # draw once; update() moves it
dd = ax.scatter([p[0, 0] for p in deep], [p[1, 0] for p in deep], color=COLORS["muted"], s=12)      # draw once; update() moves it
arrows = []                                        # an empty list to fill
ax.set_xlim(0, 2*lam); ax.set_ylim(-0.3*lam - 0.3, 0.6); ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")  # axis labels with units
ax.set_title("ka = 0.15: surface particles stay on the surface", fontsize=10)                       # the panel's message


def update(i):                                                       # draw frame i
    t = tt[i]                                        # the time of this frame [s]
    surf.set_ydata(a*np.cos(k*xg - w*t))                             # the moving surface
    dt_.set_offsets([[p[0, i], p[1, i]] for p in top]); dd.set_offsets([[p[0, i], p[1, i]] for p in deep])  # move this artist to the new frame
    for ar in arrows:  # remove last frame's arrows
        ar.remove()                                        # redraw the arrows each frame
    arrows.clear()                                        # redraw the arrows each frame
    P = np.array([top[3][0, i], top[3][1, i]])                       # one marked surface particle
    eta_f = lambda x, t_: x[1] - a*np.cos(k*x[0] - w*t_)             # level set of the surface
    g_ = np.array([-a*k*np.sin(k*P[0] - w*t), 1.0]); n = g_/np.linalg.norm(g_)   # unit normal ∇η/|∇η|
    un = uw2(P, t) @ n                                                # fluid normal speed u·n (teal)
    us = ch04.surface_normal_speed(eta_f, P, t)                       # surface normal speed (rose)
    arrows.append(ax.arrow(P[0], P[1], 1.5*un*n[0], 1.5*un*n[1], width=0.02, color=COLORS["teal"]))  # redraw the arrows each frame
    arrows.append(ax.arrow(P[0] + 0.15, P[1], 1.5*us*n[0], 1.5*us*n[1], width=0.02, color=COLORS["rose"]))  # redraw the arrows each frame
    return surf, dt_, dd                                        # the artists that changed in this frame


show_animation(animate(update, frames=frames, fig=fig, interval=60))                                # build the movie from update() and play it
""")
see_read_change(
    see="Surface particles (orange) bob and circle but stay on the moving surface; deeper particles (grey) make smaller circles; at the marked particle the teal and rose arrows have the same length.",
    read=r"""That coincidence *is* (4.91), $D\eta/Dt=0$: the surface moves normal to itself exactly as fast as the fluid on
    it (up to the small O((ka)²) error of the linear solution).""",
    change="…an evaporating surface: the rose arrow would lag the teal one by the evaporation speed and (4.92) would be nonzero.")
nb.figure(r"""
fig, (a, b, c) = plt.subplots(1, 3, figsize=(12, 3.6))                                        # the figure and its panels
yy = np.linspace(0, 0.2, 201)                                        # through the composite wall [m]
cw = ch04.two_layer_conduction(yy, 1.0, 0.1, 0.1, 0.1, 303.15, 293.15)   # temperature through two layers
a.plot(cw["T"] - 273.15, yy*100, color=COLORS["rose"], lw=2.5); a.axhline(10, ls="--", color=COLORS["muted"])  # a reference line
for yv in (5, 15):                                                    # equal heat-flux arrows on both sides
    a.annotate("", (24.5, yv + 3), (24.5, yv - 3), arrowprops=dict(arrowstyle="->", color=COLORS["rose"], lw=2))  # heat-flux arrow: upward, from the 30 °C wall (y = 0) to the 20 °C wall
a.set_xlabel("T [°C]"); a.set_ylabel("y [cm]"); a.set_title(f"(a) composite wall: q = {cw['q']:.2f} W/m² in both", fontsize=9)  # axis labels with units
y2 = np.linspace(0, 2e-3, 201)                                       # through water (below) and oil (above) [m]
tf = ch04.two_fluid_couette(y2, 1e-3, 0.1, 1e-3, 1e-3, 1.0)                                        # call the tested chapter function
b.plot(tf["u"], y2*1e3, color=COLORS["teal"], lw=2.5); b.axhline(1.0, ls="--", color=COLORS["muted"])  # a reference line
b.set_xlabel("u [m/s]"); b.set_ylabel("y [mm]"); b.set_title(f"(b) water under oil: τ = {tf['tau']:.3f} Pa in both", fontsize=9)  # axis labels with units
y3 = np.linspace(0, 1e-3, 101)                                       # Couette with a slip length [m]
for bl, ls in ((0.0, "-"), (1e-4, "--"), (5e-4, ":")):  # three slip lengths [m]
    c.plot(ch04.navier_slip_couette(y3, U=1.0, h=1e-3, slip_length=bl), y3*1e3, ls, color=COLORS["accent"], label=f"b = {bl*1e3:.1f} mm")  # draw the data
c.set_xlabel("u [m/s]"); c.set_title("(c) Couette with a slip length (our extension)", fontsize=9); c.legend(fontsize=8)  # axis labels with units
savefig(fig, "ch04", "interfaces"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Kinked profiles in (a) and (b), with equal flux arrows on both sides; in (c) the no-slip line (solid) and two slipping profiles that no longer start at u = 0.",
    read="In (a) heat flows up, from hot (30 °C at y = 0) to cold, at the same 9.09 W/m² in both layers. What is continuous is the *flux* (heat, momentum); the gradients jump in inverse proportion to k or μ. Slip is invisible unless b is a sizeable fraction of the gap.",
    change="…k₂ = k₁: no kink, one straight line.")
nb.figure(r"""
lc = ch04.capillary_length(0.0728, 998.0)                            # capillary length of water [m]
fig, ax = plt.subplots(figsize=(6.5, 3.8))                                        # the figure and its panels
for deg, col in ((0, COLORS["blue"]), (30, COLORS["teal"]), (60, COLORS["orange"]), (90, COLORS["muted"])):  # four contact angles
    th = np.deg2rad(deg)                                             # contact angle [rad]
    h = ch04.meniscus_height(th, 0.0728, 998.0)                      # height at the wall [m]
    if h > 1e-9:  # a meniscus that rises
        zeta = np.linspace(h, 0.02*h, 200)                           # heights from the wall down to nearly flat
        ax.plot(ch04.meniscus_profile_x(zeta, th)/lc, zeta/lc, color=col, lw=2, label=f"θ = {deg}°: h/ℓ_c = {h/lc:.2f}")   # closed form
        if deg in (0, 60) or not FAST:  # ODE markers (fewer when FAST)
            x_o, z_o = ch04.meniscus_profile_ode(th, x_max=4*lc)     # the ODE solution, for comparison
            ax.plot(x_o[::12]/lc, z_o[::12]/lc, "o", ms=3, color=col)                               # draw the data
    else:  # θ = 90°: no rise
        ax.plot([0, 4], [0, 0], color=col, lw=2, label="θ = 90°: flat")                             # draw the data
ax.set_xlim(0, 4); ax.set_xlabel("distance from the wall x/ℓ_c [–]"); ax.set_ylabel("ζ/ℓ_c [–]"); ax.legend(fontsize=8)  # axis labels with units
ax.set_title("A meniscus rises within a few capillary lengths of the wall")                         # the panel's message
savefig(fig, "ch04", "meniscus"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Curves rising to the wall, lower for larger contact angle, flat at 90°; the ODE markers sit on the closed-form curves.",
    read="Surface tension lifts the liquid within a few capillary lengths of the wall; gravity flattens it farther away.",
    change="…mercury on glass (θ ≈ 140°): 1 − sin θ < 1 still, but the meniscus is depressed — the formula's h² gives |h|, the sign comes from cos θ < 0.")
nb.md(r"""
> ⚠️ **Common confusion:** "$D\eta/Dt=0$ *(4.91)* is the no-slip condition." It only matches the **normal** velocities;
> the tangential velocity at a free surface is set by the stress condition (no shear from the air), at a solid wall by
> no-slip.
""")
whatif(r"""
…we measured every length in units of the body's size, every speed in units of the free stream, every time in units of
the forcing period? The equations would lose their dimensions and only a few numbers would remain — and two flows with
the same numbers would be the same flow (C15).
""")
# =====================================================================================================================
# A.11 §4.11 Dimensionless Forms of the Equations and Dynamic Similarity — R12, C15
# =====================================================================================================================
nb.section("4.11", "Dimensionless Forms of the Equations and Dynamic Similarity", intro="""
**What is this section about?** Scale every variable by the problem's own sizes and Navier–Stokes keeps only three
numbers — St, Fr, Re. Two flows with equal numbers and the same shape are the same flow: model testing works, data
collapse, and every later chapter names its regimes by these numbers.
""")
nb.recap("R12", "Dimensional analysis of sphere drag", r"""
The drag $F_D$ on a sphere of diameter d in a fluid (ρ, μ) moving at U depends on 5 variables with 3 dimensions, so on
5 − 3 = 2 groups (Ch. 1's Π theorem): $\frac{F_D}{\rho U^2d^2}=\Psi\Big(\frac{\rho Ud}{\mu}\Big)$ or equally
$\frac{F_D\rho}{\mu^2}=\Phi\Big(\frac{\mu}{\rho Ud}\Big)$ *(Eq. 4.99)* — two choices of repeating variables; the second group
is the first times Re².
""", where="Ch. 1 §1.11")
nb.code(r"""
from fluidpy.core.dimensional import pi_groups, SPHERE_DRAG          # Ch. 1's Π-theorem machinery; F, D, U, ρ, μ
print(pi_groups(SPHERE_DRAG, solution="F", repeating=["rho", "U", "D"]))   # F/(ρU²D²) and μ/(ρUD) = 1/Re
print(pi_groups(SPHERE_DRAG, solution="F", repeating=["rho", "mu", "D"]))  # Fρ/μ² and ρUD/μ = Re
""", explain="Each group is printed as its exponents: {F: 1, ρ: −1, U: −2, D: −2} is F/(ρU²D²). Two choices of repeating variables give two equivalent pairs of groups.")
core("C15", "When is a model the real thing? Dimensionless Navier–Stokes (4.101)", """
A ship designer tows a 4-m model in a tank to predict the drag of a 100-m ship. What must be the same in the tank for the
answer to be right — and can everything be matched at once?
""", eqs=("4.101",))
problem("""
Wind tunnels, towing tanks, rotating tanks and computer runs all study small or slow copies of real flows. The copy is
only useful if the physics is the same. Equations tell us exactly which combinations of sizes must agree. The same
combinations — Reynolds, Froude, Richardson, Rossby, Mach — classify the flows of every later chapter: creeping or
boundary-layer flow (Ch. 8–9), waves (Ch. 7), stratified instability (Ch. 11), geostrophic flow (Ch. 13).
""")
idea("""
x = l x*,  t = t*/Ω,  u = U u*,  p − p∞ = ρU² p*,  g = g g*                                (4.100)
(4.39b) ÷ (ρU²/l)  →  [Ωl/U] ∂u*/∂t* + (u*·∇*)u* = −∇*p* + [gl/U²] g* + [μ/ρUl] ∇*²u*     (4.101)
                          St (blue)                           1/Fr² (grey)   1/Re (rose)
same St, Fr, Re + same shape + same boundary conditions  →  same dimensionless solution
""")
note("N130", "Two routes to the groups", r"""
Dimensional analysis (R12) finds *that* groups exist; scaling the equations tells *which* groups and what each measures (a
ratio of two terms). Dynamic similarity = geometric similarity + equal groups. The figure at the end of the block shows two
dimensional Stokes-layer solutions collapsing onto one dimensionless curve.
""")
note("N132", "The parameters of a general unsteady flow", "a length l, a speed U, a frequency Ω (pulsating flow in a tube, a swimmer's stroke, a turbine's rotation), and the fluid's ρ and μ.")
P("P133", "scaled variables and the chain rule", r"""
If t* = Ωt then ∂/∂t = Ω ∂/∂t* (each second of t is Ω units of t*), and if x* = x/l then ∂/∂x = (1/l)∂/∂x*,
∂²/∂x² = (1/l²)∂²/∂x*². Scaling a derivative brings out its scale factor — the move that turns each term of an equation
into (a size) × (a dimensionless term of order one).
""", code="""
t, Om = sp.symbols('t Omega', positive=True); f = sp.Function('f')     # time, frequency, an unknown function
expr = sp.diff(f(Om*t), t)                                             # d/dt of f(Ωt)
print(sp.simplify(expr))                                               # Omega*Subs(Derivative(f(_xi_1), _xi_1), _xi_1, Omega*t): the factor Ω
""")
note("N133", "Scaled variables", r"""
(Ω an imposed frequency).

> ⚠️ The time and pressure scales change inside this section: t* = Ωt here but Ut/l for steady boundary conditions and in
> (4.109); p scaled by ρU² here, but μU/l (slow viscous flow) or ρgl (hydrostatics) are also used. `ch04.Scales` records
> which convention each use takes.
""", equation=r"x_i^*=x_i/l,\quad t^*=\Omega t,\quad u_j^*=u_j/U,\quad p^*=(p-p_\infty)/\rho U^2,\quad g_j^*=g_j/g", ref="4.100")
D("D30", ref="4.101", check_src=r"""
Om, l, U, g, mu, rho = sp.symbols('Omega l U g mu rho', positive=True)   # the scales of (4.100)
sizes = {"unsteady": rho*U*Om, "advective": rho*U**2/l, "pressure": rho*U**2/l,   # step 4: the factor each term carries
         "gravity": rho*g, "viscous": mu*U/l**2}
mine = {k: sp.simplify(v/(rho*U**2/l)) for k, v in sizes.items()}    # step 5: divide by ρU²/l
lib = ch04.nondimensional_ns_coefficients()                          # the library's brackets
print(mine)                                                          # Ωl/U, 1, 1, gl/U², μ/(ρUl)
print(ch04.nondimensional_ns_sym())                                  # the whole of (4.101)
assert all(sp.simplify(mine[k] - lib[k]) == 0 for k in mine)         # identical
""")
note("N134", "Strouhal number", r"— unsteady over advective acceleration. Vortex shedding behind a cylinder happens at St ≈ 0.2 (Ch. 9).",
     equation=r"\mathrm{St}\equiv\frac{\text{unsteady acceleration}}{\text{advective acceleration}}\propto\frac{\partial u/\partial t}{u(\partial u/\partial x)}\propto\frac{\Omega U}{U^2/l}=\frac{\Omega l}{U}", ref="4.102")
note("N135", "Reynolds number", r"""
— inertia over viscous force. Our numbers (code below): a bacterium (1 µm, 30 µm/s, water) Re ≈ 3×10⁻⁵ · a 1-cm marble
at 1 m/s in water 10⁴ · a car (4 m, 30 m/s, air) 8×10⁶ · an ocean liner (300 m, 10 m/s, water) 3×10⁹ · the Gulf Stream
(100 km, 1 m/s) 10¹¹.
""", equation=r"\mathrm{Re}\equiv\frac{\text{inertia force}}{\text{viscous force}}\propto\frac{\rho u(\partial u/\partial x)}{\mu(\partial^2u/\partial x^2)}\propto\frac{\rho U^2/l}{\mu U/l^2}=\frac{\rho Ul}{\mu}", ref="4.103")
note("N136", "Froude number", r"""
— the square root of inertia over gravity. Gravity matters dynamically only with a free surface or density differences —
a submarine deep down does not feel Fr; a ship does (waves, Ch. 7).
""", equation=r"\mathrm{Fr}\equiv\Big[\frac{\text{inertia force}}{\text{gravity force}}\Big]^{1/2}\propto\Big[\frac{\rho U^2/l}{\rho g}\Big]^{1/2}=\frac{U}{\sqrt{gl}}", ref="4.104")
note("N137", "Stratified flows", r"""
use the reduced gravity g′ = g(ρ₂ − ρ₁)/ρ₁ or the buoyancy frequency N. Numbers: an ocean thermocline Δρ/ρ = 10⁻³
(g′ = 9.81×10⁻³ m/s²), l = 100 m, U = 0.1 m/s → Ri = 98 — stratification dominates. The atmosphere with dT/dz = −6.5 K/km
(Kundu's Γ; meteorology's lapse rate Γ_met = −dT/dz = +6.5 K/km) at 288 K:
$N^2=\frac gT\big(\frac{dT}{dz}+\frac g{C_p}\big)=\frac{9.81}{288}(-6.5+9.76)\times10^{-3}=1.11\times10^{-4}$ s⁻² (stable since
dT/dz = −6.5 > Γ_a = −9.8 K/km ⇔ Γ_met = 6.5 < 9.8 K/km); with a wind shear dU/dz = 0.01 s⁻¹, Ri_g = 1.1 — above the ¼
of Ch. 11's stability criterion.

| convention | lapse rate | adiabatic | stable when |
|---|---|---|---|
| Kundu Γ ≡ dT/dz | −6.5 K/km | Γ_a = −9.8 K/km | Γ > Γ_a |
| meteorology Γ_met ≡ −dT/dz | +6.5 K/km | +9.8 K/km | Γ_met < 9.8 K/km |
""", equation=r"\mathrm{Fr}'=\frac{U}{\sqrt{g'l}}\ \text{or}\ \frac{U}{Nl},\qquad\mathrm{Ri}=\frac1{\mathrm{Fr}'^2},\qquad\mathrm{Ri}_g=\frac{N^2}{(dU/dz)^2}", ref="4.105")
nb.code(r"""
from fluidpy.core import stratification                              # Ch. 1's buoyancy-frequency helpers
N2 = stratification.brunt_vaisala_sq_from_lapse(288.0, -6.5e-3)     # N² for dT/dz = −6.5 K/km (Kundu's Γ) [1/s²]
print(f"N² = {N2:.3e} 1/s²  (dT/dz = −6.5 K/km ⇔ Γ_met = +6.5 K/km)")                               # show the numbers
print(ch04.richardson_number(9.81e-3, 100.0, 0.1), round(ch04.gradient_richardson_number(N2, 0.01), 3))   # Ri = 98.1, Ri_g ≈ 1.11
""")
note("N138", "Every output agrees", "Under dynamic similarity **every** dimensionless output agrees too — local ones (a pressure coefficient at a point) and overall ones (a drag coefficient).")
note("N139", "Pressure coefficient", r"""
Demo (code below): the ideal cylinder's surface C_p = 1 − 4 sin²θ is identical for five (U, a) pairs
(`ch04.pressure_coefficient` on `ch03.cylinder_flow` with the Bernoulli pressure).
""", equation=r"C_p\equiv\frac{p-p_\infty}{\tfrac12\rho U^2}=\Psi\Big(\mathrm{St},\mathrm{Fr},\mathrm{Re};\frac{\mathbf x}{l},\Omega t\Big)", ref="4.106")
nb.code(r"""
th = np.linspace(0, np.pi, 7)                                        # angles round the cylinder [rad]
for U, a in ((1.0, 1.0), (10.0, 0.1), (0.5, 3.0), (2.0, 0.2), (7.0, 5.0)):   # five different cylinders and streams
    us, vs = ch03.cylinder_flow(-a*(1 + 1e-9)*np.cos(th), a*(1 + 1e-9)*np.sin(th), U, a)   # surface velocity (just outside r = a) [m/s]
    p = 0.5*1.2*(U**2 - (us**2 + vs**2))                             # Bernoulli pressure relative to p∞ (air) [Pa]
    print((U, a), np.round(ch04.pressure_coefficient(p, 0.0, 1.2, U), 4))   # identical rows: C_p = 1 − 4 sin²θ
""")
note("N140", "Steady boundary conditions", "With steady boundary conditions the only time scale is l/U: t* = Ut/l and St drops out — but a flow can still become unsteady by itself (vortex shedding), at a frequency that scales with U/l (Ch. 9).")
note("N141", "A purely oscillating body", r"""
with amplitude l at frequency Ω has U = lΩ, so St = 1, Re = Ωl²/ν and Fr = Ω(l/g)^{1/2} (`ch04.Scales.from_oscillation`,
code below) — the Stokes-layer and Womersley parameters of Ch. 8 and Ch. 16.
""")
nb.code(r"""
S = ch04.Scales.from_oscillation(0.1, 2.0, 1000.0, 1e-3)             # l = 0.1 m, Ω = 2 rad/s, water
print(S.St, round(S.Re), round(S.Fr, 4))                             # St = 1, Re = Ωl²/ν = 2×10⁴, Fr = Ω√(l/g) = 0.2019
""")
note("N142", "Drag coefficient", "", equation=EQ["4.107"], ref="4.107")
note("N143", "Lift coefficient", r"""
A is a reference area: the frontal area πd²/4 for a sphere or a car, the planform (span × chord) for a plate or a wing
(`ch04.reference_area`: a 1-cm sphere 7.85×10⁻⁵ m², a 2 m × 0.5 m plate 1 m²). Wings: Ch. 14.
""", equation=EQ["4.108"], ref="4.108")
note("N131", "Sphere drag", r"""
(our version of the sphere-drag curve: the slider figure and the data-collapse figure at the end of this block): $C_D=D/(\tfrac12\rho U^2\cdot\pi d^2/4)$ against Re = ρUd/μ — 24/Re (Stokes, Ch. 8) at
small Re, a plateau near 0.4–0.5 for 10³ ≲ Re ≲ 2×10⁵, then the *drag crisis* (Ch. 9). We use the published correlation
of Morrison (2013) (a cited V5 benchmark). Number: d = 1 cm, U = 1 m/s in water → Re = 10⁴, C_D = 0.39.
""")
note("N144", "Ship drag", "For a ship, C_D = C_D(Fr, Re); far from any free surface and at low Mach number, C_D = C_D(Re) alone — the bridge to Ex. 4.8 below.")
note("N145", "Compressible flow scalings", "(t* = Ut/l now):", equation=EQ["4.109"], ref="4.109")
note("N146", "Departure from ∇·u = 0", r"""
Using (4.9)'s §4.11 re-display, $\nabla\cdot\mathbf u=-\frac1\rho\frac{D\rho}{Dt}=-\frac{1}{\rho c^2}\frac{Dp}{Dt}$ (isentropic
density changes dp = c²dρ — a general identity, not the incompressibility condition), the scalings give the line below:
the departure from ∇·u = 0 is of order M² — 0.09 at M = 0.3 (`ch04.compressibility_parameter`).
""", equation=EQ["4.110"], ref="4.110")
note("N147", "Mach number", r"""
(the square root of inertia over compressibility forces). Incompressible treatment is good below M ≈ 0.3 (C02's N12):
≈ 100 m/s in air. Ch. 14, 15.
""", equation=EQ["4.111"], ref="4.111")
note("N148", "Energy equation in enthalpy form", r"""
(stated; the book calls it "a mild revision" of (4.60) and does not write the steps; a sympy one-liner proves the two
agree, `ch04.energy_forms_sym()` → 0, below):
""", equation=EQ["4.112"], ref="4.112")
nb.code(r"""
print(ch04.energy_forms_sym())                                       # 0: (4.60) and (4.112) are the same equation
""")
note("N149", "Thermal scalings", "", equation=r"\varepsilon^*=\frac{\rho_ol^2}{\mu_oU^2}\varepsilon,\quad\mu^*=\mu/\mu_o,\quad k^*=k/k_o,\quad T^*=\frac{T-T_o}{T_w-T_o}", ref="4.113")
note("N150", "Dimensionless energy equation", r"""
(dh ≅ C_p dT). The brackets are Ec, Ec/Re and 1/(Pr Re) (`ch04.nondimensional_energy_coefficients()`).

> ⚠️ The book introduces it with "those defined in (4.106), (4.107)"; it means the scalings (4.109),
> $t^*=Ut/l,\ p^*=(p-p_\infty)/\rho_oU^2,\ \rho^*=\rho/\rho_o$, and the equation (4.112).
""", equation=r"\rho^*\frac{DT^*}{Dt^*}=\Big[\frac{U^2}{C_p(T_w-T_o)}\Big]\frac{Dp^*}{Dt^*}+\Big[\frac{\mu_oU}{\rho_oC_p(T_w-T_o)l}\Big]\varepsilon^*+\Big[\frac{k_o}{\rho_oC_pUl}\Big]\frac{\partial}{\partial x_i^*}\Big(k^*\frac{\partial T^*}{\partial x_i^*}\Big)", ref="4.114")
note("N151", "Eckert number", r"""
Small Ec ⇒ pressure work and viscous heating drop out and (4.112) becomes the Boussinesq heat equation
$\frac{DT}{Dt}=\kappa\nabla^2T$ *(4.89)* — the hidden assumption of C13.
""", equation=EQ["4.115"], ref="4.115")
note("N152", "Prandtl number", r"""
— momentum diffusivity over heat diffusivity. Our values from `ch01.FLUIDS` (`ch04.prandtl_of`, code below): air ≈ 0.709,
water at 20 °C ≈ 6.95; kinetic theory for a monatomic gas 2/3; Eucken's estimate 4γ/(9γ − 5) = 0.74 for γ = 1.4.
Published property tables (the book's included) give values 1.5–2.2 % different at 20 °C — Pr is a ratio of three
measured properties, and tables disagree at that level. Thermal boundary layers: Ch. 9; convection: Ch. 11.
""", equation=r"\mathrm{Pr}\equiv\frac{\nu}{\kappa}=\frac{\mu_oC_p}{k_o}", ref="4.116")
note("N153", "Weber number", "(a ratio of *forces*, not per volume). Break-up of drops and jets (Ch. 7, 16).", equation=EQ["4.117"], ref="4.117")
note("N154", "Bond number", "(the capillary length of C14).", equation=r"\mathrm{Bo}\equiv\frac{\rho l^2g}{\sigma}=(l/\ell_c)^2", ref="4.118")
note("N155", "Capillary number", "Coating flows and flow in porous media (Ch. 16).", equation=r"\mathrm{Ca}\equiv\frac{\mu U}{\sigma}=\mathrm{We}/\mathrm{Re}", ref="4.119")
note("N156", "Ex. 4.8 — the ship model", r"""
(stated with our own numbers; the book's are in the private test data). A 1:25 model of a 100-m ship designed for
10 m/s. **Froude matching** $U_m=U_p\sqrt{l_m/l_p}=10/5=2$ m/s makes the wave pattern similar — but then
Re_m/Re_p = (l_m/l_p)^{3/2} = 1/125: the model's boundary layers are far too viscous. Froude's way out: split the drag into
friction (depends on Re, from a flat-plate friction line) and wave drag (depends on Fr, scaled by
$(\rho_p/\rho_m)(l_p/l_m)^2(U_p/U_m)^2=(l_p/l_m)^3$ = 25³ for equal densities, i.e. λ⁻³ with λ = l_m/l_p = 1/25). Numbers in the worked example;
`ch04.ship_drag_extrapolation`.

> ⚠️ The book's total 9.14×10⁵ N is 9.15×10⁵ N without the intermediate rounding — a detail, not an error in the method.
""")
nb.worked_example("similarity by hand", r"""
1. **Sphere**: d = 1 cm, U = 1 m/s, water ν = 10⁻⁶ m²/s: Re = Ud/ν = 10⁴; C_D ≈ 0.39 (plateau) → drag
   ½ρU²C_D πd²/4 = 0.5 × 1000 × 1 × 0.39 × 7.85×10⁻⁵ = 0.015 N.
2. **The same Re in air** (ν = 1.5×10⁻⁵): a 1-cm sphere needs U = 15 m/s — then its drag coefficient is the same 0.39.
3. **Ship**: λ = 1/25, U_p = 10 m/s → U_m = 2 m/s (Fr = 0.32 both); Re_p = 10 × 100/10⁻⁶ = 10⁹, Re_m = 2 × 4/10⁻⁶ = 8×10⁶
   — off by 125.
4. **Drag extrapolation** (model wetted area 4 m², prototype 2500 m²; the ITTC 1957 friction line
   C_f = 0.075/(log₁₀Re − 2)², a standard flat-plate correlation): C_f,m = 0.00312, C_f,p = 0.00153. Model measures 40 N
   total; friction ½ρU²SC_f = ½ × 1000 × 4 × 4 × 0.00312 = 25.0 N; wave 15.0 N. Wave drag scales by
   (1025/1000) × 25² × 5² = 16 016 → 240.9 kN. Prototype friction ½ × 1025 × 100 × 2500 × 0.00153 = 196.0 kN. Total ≈ 437
   kN. Scaling the whole 40 N by 16 016 would give 641 kN — 47 % too much, because friction does not scale with Fr.
""")
nb.code(r"""
for name, (U, l, nu) in {"bacterium": (3e-5, 1e-6, 1e-6), "marble": (1.0, 0.01, 1e-6), "car": (30.0, 4.0, 1.5e-5),  # five flows from microbes to ocean currents
                         "liner": (10.0, 300.0, 1e-6), "Gulf Stream": (1.0, 1e5, 1e-6)}.items():   # speed [m/s], size [m], ν [m²/s]
    print(f"{name:12s} Re = {ch04.reynolds_number(U, l, nu):.1e}")  # 16 orders of magnitude
print(ch04.sphere_drag_coefficient(1e4), ch04.sphere_drag_coefficient(0.1, model="stokes"))   # 0.3926 (Morrison), 240 (24/Re)
print(ch04.froude_scaled_speed(10.0, 100.0, 4.0), ch04.froude_number(10.0, 100.0), ch04.froude_number(2.0, 4.0))   # 2 m/s; Fr 0.3193 twice
mp = ch04.model_prototype(100.0, 10.0, 1/25)                          # the 1:25 model, Froude-matched
print(mp["U_m"], round(mp["Re_ratio"], 1), mp["matched"])             # 2.0 m/s, Re off by 125, only Fr matched
print(ch04.ship_drag_extrapolation(100.0, 10.0, 2500.0, 1/25, 40.0, 0.00312, 0.00153))   # Froude's extrapolation [N]
print(ch04.prandtl_of("air"), ch04.prandtl_of("water"), ch04.eucken_prandtl(1.4))   # 0.709, 6.95, 0.737
print(ch04.nondimensional_ns_coefficients())                          # the brackets of (4.101)
""", explain=r"""
1. Reynolds numbers across 16 orders of magnitude.
2. The sphere's drag coefficient in two regimes.
3. Froude scaling keeps Fr equal.
4. The model–prototype summary: Fr matched, Re not (Re_p/Re_m = $(l_p/l_m)^{3/2}$ = 125).
5. Froude's extrapolation: friction and wave drag scaled separately (437 kN instead of the naive 641 kN).
6. The Prandtl numbers and Eucken's estimate.
7. The brackets of (4.101), $\big[\frac{\Omega l}{U}\big]$, $\big[\frac{gl}{U^2}\big]$, $\big[\frac{\mu}{\rho Ul}\big]$.
""")
nb.check_agree(r"""
U_p, l_p, nu = 10.0, 100.0, 1e-6                                     # the ship [m/s], [m], water ν [m²/s]
U_m, l_m = 2.0, 4.0                                                  # the model
Re_hand = (U_p*l_p/nu, U_m*l_m/nu); Fr_hand = (U_p/np.sqrt(9.81*l_p), U_m/np.sqrt(9.81*l_m))   # by hand
assert np.isclose(Fr_hand[0], Fr_hand[1]) and np.isclose(Re_hand[0]/Re_hand[1], 125.0)   # Fr matched, Re off by 125
x_, t_, l, U, Om, rho, mu = sp.symbols('x t l U Omega rho mu', positive=True)   # dimensional symbols
us = sp.Function('u_s')                                              # the dimensionless velocity u*(x*, t*)
u = U*us(x_/l, Om*t_)                                                # (4.100): u = U u*(x/l, Ωt)
dus_dts = sp.diff(us(x_/l, Om*t_), t_) / Om                         # ∂u*/∂t* (the chain rule brought out a factor Ω)
d2us_dxs2 = sp.diff(us(x_/l, Om*t_), x_, 2) * l**2                   # ∂²u*/∂x*² (the chain rule brought out 1/l²)
unsteady = sp.simplify(rho*sp.diff(u, t_) / (rho*U**2/l) / dus_dts)   # ρ ∂u/∂t divided by ρU²/l, per unit ∂u*/∂t*
viscous = sp.simplify(mu*sp.diff(u, x_, 2) / (rho*U**2/l) / d2us_dxs2)   # μ ∂²u/∂x² divided by ρU²/l, per unit ∂²u*/∂x*²
print(unsteady, viscous)                                             # Omega*l/U and mu/(U*l*rho)
lib = ch04.nondimensional_ns_coefficients()                                        # call the tested chapter function
assert sp.simplify(unsteady - lib["unsteady"]) == 0 and sp.simplify(viscous - lib["viscous"]) == 0
""")
nb.figure(r"""
yy = np.linspace(0, 12e-3, 200)                                      # height above the plate [m]
cases = [(1.0, 1e-6, 1.0, COLORS["teal"], "U = 1 m/s, ν = 1e-6 m²/s, t = 1 s"),                     # the cases to compare
         (5.0, 1e-5, 0.1, COLORS["orange"], "U = 5 m/s, ν = 1e-5 m²/s, t = 0.1 s")]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6))                                        # the figure and its panels
for U, nu, t, col, lab in cases:                                     # two dimensional Stokes layers
    u = ch04.stokes_first_problem(yy, t, U=U, nu=nu)                                        # call the tested chapter function
    a.plot(u, yy*1e3, color=col, lw=2.5, label=lab)                  # different curves
    b.plot(u/U, yy/(2*np.sqrt(nu*t)), "o" if U > 1 else "-", color=col, ms=3, label=lab)   # the same curve
a.set_xlabel("u [m/s]"); a.set_ylabel("y [mm]"); a.legend(fontsize=7); a.set_title("dimensional: two flows", fontsize=10)  # axis labels with units
b.set_xlabel("u/U [–]"); b.set_ylabel("y/(2√(νt)) [–]"); b.set_ylim(0, 3); b.legend(fontsize=7)     # axis labels with units
b.set_title("dimensionless: one erfc curve", fontsize=10)                                        # the panel's message
savefig(fig, "ch04", "stokes_collapse"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Two different curves on the left, one curve on the right (dots on the line).",
    read="After scaling, the equation has no parameters left, so every Stokes layer is the same curve u/U = erfc(y/2√(νt)).",
    change="…a plate oscillating instead of started: a new group (Ωl²/ν) appears and the collapse needs it to match too.")
nb.plotly(r"""
Re_axis = np.logspace(-1, 6, 200)                                    # Reynolds numbers
CD_curve = ch04.sphere_drag_coefficient(Re_axis)                     # Morrison (2013) correlation
def point(logRe):                                                    # the slider picks log10(Re)
    Re = 10**logRe  # the slider's Reynolds number
    return {"C_D(Re), Morrison 2013": (np.log10(Re_axis), np.log10(CD_curve)),                      # the curves for this slider value
            "Stokes 24/Re": (np.log10(Re_axis[Re_axis < 10]), np.log10(24/Re_axis[Re_axis < 10])),  # another named curve (x, y)
            "this Re": ([logRe], [np.log10(ch04.sphere_drag_coefficient(Re))])}                     # another named curve (x, y)
fig = slider_figure(point, "log10 Re", np.linspace(-1, 6, 36 if not FAST else 18), xlabel="log10 Re [–]", ylabel="log10 C_D [–]",  # precompute every slider position (works on the web page)
                    title="One curve for every sphere: C_D depends on Re alone", modes={"this Re": "markers"})  # labels, title and a fixed range
recolor(fig, {"C_D(Re), Morrison 2013": COLORS["ink"], "Stokes 24/Re": COLORS["muted"], "this Re": COLORS["rose"]},  # house colours for each trace, then draw
        {"Stokes 24/Re": "dash"}).show()                                        # dashed reference line
""")
see_read_change(
    see="A black curve falling steeply (on the grey Stokes line) at small Re, flattening near C_D ≈ 0.4, then dropping sharply near Re ≈ 3×10⁵; the rose point rides on it.",
    read="At the slider's Re the point shows C_D; any sphere, any fluid, any speed with that Re sits on it.",
    change="…a rough sphere: the drag crisis moves to lower Re (a new group, roughness/d, enters).")
nb.figure(r"""
d = ch04.synthetic_sphere_drag_data(60, seed=0, noise=0.03)          # 60 "experiments": fluid, d, U, ρ, μ, F, Re, C_D (3 % noise)
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6))                                        # the figure and its panels
cols = {f: c for f, c in zip(sorted(set(d["fluid"])), (COLORS["blue"], COLORS["orange"], COLORS["teal"], COLORS["accent"]))}  # one colour / label per series
for f in cols:                                                       # one colour per fluid
    m = np.array(d["fluid"]) == f  # the experiments done in this fluid
    a.loglog(np.asarray(d["U"])[m], np.asarray(d["CD"])[m], "o", color=cols[f], ms=4, label=f)      # draw the data
    b.loglog(np.asarray(d["Re"])[m], np.asarray(d["CD"])[m], "o", color=cols[f], ms=4)              # draw the data
Rg = np.logspace(-1, 6, 300); b.loglog(Rg, ch04.sphere_drag_coefficient(Rg), color=COLORS["ink"], lw=1)  # the Morrison curve for reference
a.set_xlabel("speed U [m/s]"); a.set_ylabel("C_D [–]"); a.legend(fontsize=7); a.set_title("against a dimensional variable: a cloud", fontsize=10)  # axis labels with units
b.set_xlabel("Re [–]"); b.set_title("against the right group: one curve", fontsize=10)              # axis labels with units
savefig(fig, "ch04", "drag_collapse"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="Chaos on the left, order on the right.",
    read="Plotting against the right group removes the dimensional scatter — that is dynamic similarity measured.",
    change="…plotting against d alone: still a cloud; only the combination ρUd/μ collapses the data.")
nb.figure(r"""
sd = ch04.ship_drag_extrapolation(100.0, 10.0, 2500.0, 1/25, 40.0, 0.00312, 0.00153)   # Froude's method, our numbers [N]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6))                                        # the figure and its panels
a.bar(["model total", "friction", "wave"], [40.0, sd["D_m_friction"], sd["D_m_wave"]],              # draw the data
      color=[COLORS["ink"], COLORS["rose"], COLORS["blue"]])         # the model's drag split
a.set_ylabel("model drag [N]"); a.set_title("(a) split the model's drag", fontsize=10)              # axis labels with units
b.bar(["wave × 16 016", "prototype friction", "prototype total", "naive: total × 16 016"],          # draw the data
      np.array([sd["D_p_wave"], sd["D_p_friction"], sd["D_p_total"], sd["D_p_uncorrected"]])/1e3,  # prototype drag parts [kN]
      color=[COLORS["blue"], COLORS["rose"], COLORS["ink"], COLORS["muted"]], hatch=["", "", "", "//"])  # the naive bar hatched
b.set_ylabel("prototype drag [kN]"); b.tick_params(axis="x", labelsize=7); b.set_title("(b) scale each part by its own law", fontsize=10)  # axis labels with units
savefig(fig, "ch04", "ship_drag"); plt.show()                                        # save the PNG to outputs/ch04, then draw
""", see="The corrected prototype total (437 kN) much shorter than the naive grey bar (641 kN).",
    read="Friction scales with Re, waves with Fr; the model cannot match both, so each part is scaled by its own law.",
    change="…a model in a fluid of much lower ν (none practical exists): both Fr and Re could match — the reason towing tanks live with the correction.")
nb.explainer("dynamic_similarity_models", heading="When does a model behave like the real thing?", why=r"""
**Why interactive rather than a static figure:** similarity is about matching groups between two flows; changing scale, speed and fluid lets you watch which badges turn green and which cannot, and the data collapse only when the axis is the right group.

A prototype and a model side by side, the groups of
$\big[\frac{\Omega l}{U}\big]\frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\big[\frac{gl}{U^2}\big]\mathbf g^*+\big[\frac{\mu}{\rho Ul}\big]\nabla^{*2}\mathbf u^*$
(4.101) as paired bars with "matched" badges, and the consequences: sphere data collapsing on C_D(Re), a ship's drag split
and scaled. Modes: sphere, ship, stratified flow (Fr′, Ri), rotating tank (Rossby number U/(2Ωl), the advective term
U²/l over the Coriolis term 2ΩU of (4.45) — a forward pointer to Ch. 13).
""", tries=[
    "Ship preset 1:25: match Fr and read the Re mismatch ×125 in Explain.",
    "Sphere mode: toggle the axes from dimensional to Re — the points collapse.",
    "Try to match Re too by changing the model fluid: which fluid would you need?",
    "Rotating-tank mode: match the Rossby number U/(2Ωl) of an atmospheric system in a 1-m tank.",
])
whatif(r"""
…the flow had rotation and stratification together, as the ocean and atmosphere do? Two more groups join — the Rossby
number U/(2Ωl), from the Coriolis term −2Ω × u′ of (4.45), $\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho[\mathbf g-2\boldsymbol\Omega\times\mathbf u'-\dots]+\mu\nabla'^2\mathbf u'$,
and the Richardson number above — and the regime map of Ch. 13 is drawn in their plane.
""")

# =====================================================================================================================
# A.12 End matter
# =====================================================================================================================
nb.pointer("""
S01 · Exercises: the book's Exercises 4.1–4.63 are for practice. The derivations the text defers to Exercises 4.7, 4.8,
4.30, 4.38, 4.42, 4.43, 4.45, 4.46, 4.47 and 4.50 are written out above in our own words (D04, D08, D13, D15, D18, D20,
D21, D23, D24).
""")
nb.pointer("""
S02 · Literature: the chapter's references — Aris (isotropic tensors, used in C07), Batchelor and Lamb (classical
treatments), Spiegel & Veronis (the Boussinesq approximation, C13) — are listed at the end of the book's chapter.
""")
nb.summary(
    clicked=[
        r"**C01** A box of any shape and motion keeps a mass budget: $\frac{d}{dt}\int_{V^*}\rho\,dV+\int_{A^*}\rho(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0$ *(4.5)* — flux counts relative to the walls.",
        r"**C02** At a point, $\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0$ *(4.7)*; incompressible means each particle keeps its density, so $\nabla\cdot\mathbf u=0$ *(4.10)*.",
        r"**C03** In 2-D one scalar ψ carries the flow: $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$; contours are streamlines and Δψ is the flux between them.",
        r"**C04** Forces follow from fluxes: $\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA$ *(4.17)* weighs the drag on a bar from its wake.",
        r"**C05** Along one streamline of a steady, frictionless, constant-density flow $\tfrac12U^2+gz+p/\rho$ is constant *(4.19)*.",
        r"**C06** Newton for every continuum: $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(4.24)*, stress divergence over the first index.",
        r"**C07** A Newtonian fluid feels only its strain rate: $\tau_{ij}=-p\delta_{ij}+2\mu S_{ij}+\lambda S_{mm}\delta_{ij}$ *(4.31)*, two constants.",
        r"**C08** $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(4.39b)* is a budget of five accelerations; each flow wakes only some of them.",
        r"**C09** A rotating observer adds $-2\boldsymbol\Omega\times\mathbf u'$ (Coriolis) and $-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')$ (centrifugal) to the forces of (4.45); the ball flies straight, the floor turns.",
        r"**C10** Viscosity turns kinetic energy into heat at the rate $\varepsilon=2\nu(S_{ij}-\tfrac13S_{mm}\delta_{ij})^2+\frac{\mu_v}\rho S_{mm}^2\ge0$ *(4.58)*, a sum of squares — never the reverse.",
        r"**C11** Steady inviscid barotropic flow keeps $B=\tfrac12u^2+\int dp/\rho+gz$ constant along streamlines and vortex lines *(4.71)*, everywhere if irrotational *(4.72)*.",
        r"**C12** Irrotational unsteady flow obeys $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=$ const *(4.75)*: accelerating fluid needs pressure even at rest.",
        r"**C13** Boussinesq keeps density only where it meets gravity: $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ *(4.86)* and $\frac{DT}{Dt}=\kappa\nabla^2T$ *(4.89)*.",
        r"**C14** A surface made of fluid moves with it: $D\eta/Dt=0$ on η = 0 *(4.91)*; fluxes are continuous across interfaces; walls impose no-slip.",
        r"**C15** Scaled, Navier–Stokes keeps St, Fr, Re *(4.101)*; equal groups make a model the real thing.",
    ],
    feeds_forward=[
        r"Ch. 5: the vorticity equation from the Lamb identity $u_i\partial u_j/\partial x_i=-(\mathbf u\times\boldsymbol\omega)_j+\partial(\tfrac12u_i^2)/\partial x_j$ *(4.68)* and (4.39b).",
        r"Ch. 6: potential flow with B constant everywhere *(4.72)* and unsteady Bernoulli *(4.75)*.",
        r"Ch. 7: free-surface conditions — $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+p/\rho=$ const *(4.75)* and $D\eta/Dt=0$ *(4.91)*.",
        "Ch. 8: exact solutions checked term by term with `core.navier_stokes`.",
        r"Ch. 9: the momentum integral from the box budget *(4.17)*.",
        r"Ch. 10: conservative forms — the flux form *(4.22)* and the total-energy equation *(4.53)*.",
        "Ch. 11: Boussinesq instabilities and the Richardson number.",
        "Ch. 12: ε and the kinetic-energy budget (the force-work/deformation-work split of (4.54)).",
        r"Ch. 13: Navier–Stokes in a rotating frame *(4.45)* with Boussinesq, Ro and Ri.",
        "Ch. 15: the energy equation and the Mach number.",
    ],
    left_out=[
        "Non-Newtonian fluids (Ch. 16).",
        "The proof of the isotropic-tensor form (cited, Aris).",
        "Kelvin's circulation theorem (Ch. 5).",
        "Marangoni flows (interfacial-flow texts).",
        "Anelastic equations for deep atmospheres (Ch. 13 texts).",
    ],
)

if __name__ == "__main__":
    print("equations written next to their numbers in", finalize_equations(), "markdown cells")
    path = nb.save()
    print(f"wrote {path.relative_to(ROOT)}: {len(nb.cells)} cells, {len(nb.cores)} CORE, {len(nb.recaps)} RECAP, "
          f"{len(nb.derivations)} derivations, {len(nb.primers)} primers, {len(nb.explainers)} explainers")
