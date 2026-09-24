"""Build the Chapter 6 teaching notebook: ``notebooks/ch06_ideal_flow.ipynb``.

Source of truth: ``analysis/ch06_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract
— every call goes to ``fluidpy.ch06_ideal_flow``, imported as ``ch06``, which re-exports the four new core modules
``core.potential`` (pf), ``core.conformal`` (cm), ``core.laplace_solvers`` (ls) and ``core.panels`` (pn)), Part E
(prerequisite ledger → 16 primers P149–P164 and one-line reminders of earlier primers), Part F (the 31 derivations
D01–D31, one move per step), **Part G (errata from verification, which override A, B and F: Example 6.2 converges at
order 4/3 because of the 270° corner; the D21 principal-root check excludes the slit; the axial-method strengths
alternate in sign and only their moments converge; constant-source panels are exact on a circle)** and
``analysis/ch06_curation.md`` (IDs, depths, section coverage). Physics lives in ``fluidpy``; cells only call it.

**Derivations are read from Part F at build time** (``part_f()`` below, the ch03–ch05 parser): goal, start, plan, tools,
assumptions, every step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word, so the
notebook and the design cannot drift apart. Explainer-only fields (*live*, *set*, *watch*) are dropped. The "check"
texts are edited (``pf_sub``) so that every check names a cell of this notebook that really runs it. The five ★★★
sympy checks (D17, D18, D21, D29, D31) are written here, every line commented, and re-run the construction.

Conventions (design header, binding): Γ counterclockwise positive in code (``pf.Vortex``); the book's (6.36)–(6.40),
(6.52), (6.61)–(6.62), (6.68) and Example 6.1 use a clockwise Γ (``Gamma_cw=``); the 2-D dipole points from the sink to
the source; θ from +x; in §6.8 z is the horizontal axis along the stream; in §6.9 w is the z-velocity; D, L are forces
on the body per unit depth; the book's slips are taught in corrected form ((6.61) coefficient, (6.104) bracket, (6.108)
stray dφ, the §6.3 and §6.4 cross-references, Example 6.2's loop index and unit).

Run:  .venv/Scripts/python.exe notebooks/build_ch06.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch06.py --dump     (prints the parsed Part F derivations only)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch06")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# Ch. 6 from analysis/ch06.md §2 (transcribed from the rendered pages p225–p267; (6.61) re-read on p248 here);
# earlier chapters' equations from their notebooks ((3.17) re-read on ch03 p106: it is ω = 0).
# ---------------------------------------------------------------------------------------------------------------------
EQ = {
    "6.1": r"\nabla\cdot\mathbf u=0\ \text{ and }\ \rho\,\frac{D\mathbf u}{Dt}=-\nabla p",
    "6.2": r"\partial u/\partial x+\partial v/\partial y=0",
    "6.3": r"u\equiv\partial\psi/\partial y,\ v\equiv-\partial\psi/\partial x",
    "6.4": r"\omega_z=\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}=-\nabla^2\psi",
    "6.5": r"\nabla^2\psi=0",
    "6.6": r"\nabla^2\psi=-\Gamma\,\delta(x-x')\,\delta(y-y')",
    "6.7": r"\psi=-Vx+Uy",
    "6.8": r"\psi=-\frac{\Gamma}{2\pi}\ln\sqrt{(x-x')^2+(y-y')^2}",
    "6.9": r"\partial v/\partial x-\partial u/\partial y=0",
    "6.10": r"u\equiv\partial\phi/\partial x,\ v\equiv\partial\phi/\partial y",
    "6.11": r"\nabla^2\phi=q(x,y)",
    "6.12": r"\nabla^2\phi=0",
    "6.13": r"\nabla^2\phi=m\,\delta(x-x')\,\delta(y-y')",
    "6.14": r"\phi=Ux+Vy",
    "6.15": r"\phi=\frac{m}{2\pi}\ln\sqrt{(x-x')^2+(y-y')^2}",
    "6.16": r"\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{surface}};\ \partial\phi/\partial n=0\ \text{or}\ \partial\psi/\partial s=0\ \text{(stationary)}",
    "6.17": r"\partial\phi/\partial x=U\ \text{or}\ \partial\psi/\partial y=U\ \text{far away}",
    "6.18": r"p+\tfrac12\rho\lvert\nabla\phi\rvert^2=p+\tfrac12\rho\lvert\nabla\psi\rvert^2=\text{const}",
    "6.19": r"\frac1r\frac{\partial}{\partial r}(ru_r)+\frac1r\frac{\partial u_\theta}{\partial\theta}=0",
    "6.20": r"\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}=0",
    "6.21": r"u_r=\frac{\partial\phi}{\partial r}=\frac1r\frac{\partial\psi}{\partial\theta}",
    "6.22": r"u_\theta=\frac1r\frac{\partial\phi}{\partial\theta}=-\frac{\partial\psi}{\partial r}",
    "6.23a": r"\nabla^2\psi=\frac1r\frac{\partial}{\partial r}\Big(r\frac{\partial\psi}{\partial r}\Big)+\frac1{r^2}\frac{\partial^2\psi}{\partial\theta^2}=0",
    "6.23b": r"\nabla^2\phi=\frac1r\frac{\partial}{\partial r}\Big(r\frac{\partial\phi}{\partial r}\Big)+\frac1{r^2}\frac{\partial^2\phi}{\partial\theta^2}=0",
    "6.24": r"\psi=2Axy",
    "6.25": r"\phi=2Axy",
    "6.26": r"\psi=A(x^2-y^2)",
    "6.27": r"\phi=A(x^2-y^2)",
    "6.28": r"\phi=\frac{m}{2\pi}\ln\sqrt{(x+\varepsilon)^2+y^2}-\frac{m}{2\pi}\ln\sqrt{(x-\varepsilon)^2+y^2}",
    "6.29": r"\phi=\frac{m\varepsilon}{\pi}\frac{x}{r^2}=-\frac{\mathbf d\cdot\mathbf x}{2\pi r^2}=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{\cos\theta}{r}",
    "6.30": r"\phi=Ux+\frac{m}{2\pi}\ln\sqrt{x^2+y^2}=Ur\cos\theta+\frac{m}{2\pi}\ln r",
    "6.31": r"\psi=Uy+\frac{m}{2\pi}\tan^{-1}\frac yx=Ur\sin\theta+\frac{m}{2\pi}\theta",
    "6.32": r"C_p=\frac{p-p_\infty}{\frac12\rho U^2}=1-\frac{\lvert\mathbf u\rvert^2}{U^2}",
    "6.33": r"\phi=U\Big(r+\frac{a^2}{r}\Big)\cos\theta,\ \psi=U\Big(r-\frac{a^2}{r}\Big)\sin\theta",
    "6.34": r"u_r=U\Big(1-\frac{a^2}{r^2}\Big)\cos\theta,\ u_\theta=-U\Big(1+\frac{a^2}{r^2}\Big)\sin\theta",
    "6.35": r"C_p(r=a,\theta)=1-4\sin^2\theta",
    "6.36": r"\psi=U\Big(r-\frac{a^2}{r}\Big)\sin\theta+\frac{\Gamma}{2\pi}\ln\Big(\frac ra\Big)",
    "6.37": r"u_\theta(r=a,\theta)=-2U\sin\theta-\Gamma/2\pi a",
    "6.38": r"\sin\theta=-\Gamma/4\pi aU",
    "6.39": r"p(r=a,\theta)=p_\infty+\tfrac12\rho\Big[U^2-\Big(-2U\sin\theta-\frac{\Gamma}{2\pi a}\Big)^2\Big]",
    "6.40": r"L=\rho U\Gamma",
    "6.41": r"\psi=\frac m{2\pi}\Big[\tan^{-1}\Big(\frac y{x+a}\Big)+\tan^{-1}\Big(\frac y{x-a}\Big)\Big]",
    "6.42": r"w\equiv\phi+i\psi",
    "6.43": r"z\equiv x+iy=re^{i\theta}",
    "6.44": r"\frac{\partial\phi}{\partial x}=\frac{\partial\psi}{\partial y},\ \frac{\partial\phi}{\partial y}=-\frac{\partial\psi}{\partial x}",
    "6.45": r"dw/dz=u-iv",
    "6.46": r"w(z)=Az^n=Ar^n(\cos n\theta+i\sin n\theta)",
    "6.47": r"w=-\frac{i\Gamma}{2\pi}\ln(z-z')",
    "6.48": r"w=\frac m{2\pi}\ln(z-z')",
    "6.49": r"w=\frac{d}{2\pi(z-z')}",
    "6.50": r"w=Uz+\frac m{2\pi}\ln z",
    "6.51": r"w=U\Big(z+\frac{a^2}z\Big)",
    "6.52": r"w=U\Big(z+\frac{a^2}z\Big)+\frac{i\Gamma}{2\pi}\ln(z/a)",
    "6.53": r"w=\frac m{2\pi}\ln\Big(\frac{z^2-a^2}{a^2}\Big)",
    "6.54": r"\int_{A^*}\rho\mathbf u(\mathbf u\cdot\mathbf n)\,dA=-\int_{A^*}p\,\mathbf n\,dA+\mathbf F",
    "6.55": r"D\mathbf e_x+L\mathbf e_y=-\frac1B\int_{A^*}p\,\mathbf n\,dA",
    "6.56": r"D=-\oint_Cp\,dy,\ L=\oint_Cp\,dx",
    "6.57": r"D-iL=-i\oint_Cp\,dz^*",
    "6.58": r"D-iL=-i\oint_C\big[p_\infty+\tfrac12\rho U^2-\tfrac12\rho(u-iv)(u+iv)\big]dz^*",
    "6.59": r"(u+iv)\,dz^*=(u-iv)\,dz=(dw/dz)\,dz",
    "6.60": r"D-iL=\frac{i\rho}{2}\oint_C\Big(\frac{dw}{dz}\Big)^2dz",
    "6.61": r"D-iL=\frac{i\rho}{2}\oint_C\Big(U+\frac{i\Gamma}{2\pi z}-\frac{d}{2\pi z^2}+\dots\Big)^2dz",
    "6.62": r"D-iL=\frac{i\rho}{2}2\pi i\Big(\frac{iU\Gamma}{\pi}\Big)=-i\rho U\Gamma,\ D=0,\ L=\rho U\Gamma",
    "6.63": r"\delta w=\frac{dw}{dz}\delta z",
    "6.64": r"\delta'w=\frac{dw}{dz}\delta'z,\ \alpha=\beta",
    "6.65": r"z=\zeta+\frac{b^2}{\zeta}",
    "6.66": r"z=ae^{i\theta}+\frac{b^2}{a}e^{-i\theta}",
    "6.67": r"\frac{x^2}{(a+b^2/a)^2}+\frac{y^2}{(a-b^2/a)^2}=1",
    "6.68": r"w=U\Big(\zeta+\frac{a^2}{\zeta}\Big)+\frac{i\Gamma}{2\pi}\ln(\zeta/a)",
    "6.69": r"\zeta=\tfrac12z+\tfrac12(z^2-4b^2)^{1/2}",
    "6.70": r"\Big(\frac{\partial^2\psi}{\partial x^2}\Big)_{i,j}\simeq\frac{\psi_{i+1,j}-2\psi_{i,j}+\psi_{i-1,j}}{\Delta x^2}",
    "6.71": r"\Big(\frac{\partial^2\psi}{\partial y^2}\Big)_{i,j}\simeq\frac{\psi_{i,j+1}-2\psi_{i,j}+\psi_{i,j-1}}{\Delta y^2}",
    "6.72": r"\psi_{i,j}=\tfrac14\big[\psi_{i-1,j}+\psi_{i+1,j}+\psi_{i,j-1}+\psi_{i,j+1}\big]",
    "6.73": r"\psi_{2,2}=\tfrac14\big[\psi^B_{1,2}+\psi_{3,2}+\psi^B_{2,1}+\psi_{2,3}\big],\ \dots",
    "6.74": r"\frac1R\frac{\partial}{\partial R}(Ru_R)+\frac{\partial u_z}{\partial z}=0",
    "6.75": r"u_R=-\frac1R\frac{\partial\psi}{\partial z},\ u_z=\frac1R\frac{\partial\psi}{\partial R}",
    "6.76": r"\omega_\varphi=\frac{\partial u_R}{\partial z}-\frac{\partial u_z}{\partial R}",
    "6.77": r"\frac{\partial}{\partial R}\Big(\frac1R\frac{\partial\psi}{\partial R}\Big)+\frac1R\frac{\partial^2\psi}{\partial z^2}=-\omega_\varphi=0",
    "6.78": r"dQ=2\pi R(\mathbf u\cdot\mathbf n)ds=2\pi\,d\psi",
    "6.79": r"u_R=\partial\phi/\partial R,\ u_z=\partial\phi/\partial z",
    "6.80": r"\frac1R\frac{\partial}{\partial R}\Big(R\frac{\partial\phi}{\partial R}\Big)+\frac{\partial^2\phi}{\partial z^2}=0",
    "6.81": r"R=r\sin\theta,\ z=r\cos\theta",
    "6.82": r"\frac1r\frac{\partial}{\partial r}(r^2u_r)+\frac1{\sin\theta}\frac{\partial}{\partial\theta}(u_\theta\sin\theta)=0",
    "6.83": r"u_r=\frac1{r^2\sin\theta}\frac{\partial\psi}{\partial\theta}=\frac{\partial\phi}{\partial r},\ u_\theta=-\frac1{r\sin\theta}\frac{\partial\psi}{\partial r}=\frac1r\frac{\partial\phi}{\partial\theta}",
    "6.84": r"\omega_\varphi=\frac1r\Big[\frac{\partial}{\partial r}(ru_\theta)-\frac{\partial u_r}{\partial\theta}\Big]",
    "6.85": r"\frac1{r^2}\frac{\partial}{\partial r}\Big(r^2\frac{\partial\phi}{\partial r}\Big)+\frac1{r^2\sin\theta}\frac{\partial}{\partial\theta}\Big(\sin\theta\frac{\partial\phi}{\partial\theta}\Big)=0",
    "6.86": r"\phi=Ur\cos\theta,\ \psi=\tfrac12Ur^2\sin^2\theta",
    "6.87": r"\phi=-\frac Q{4\pi r},\ \psi=-\frac Q{4\pi}\cos\theta",
    "6.88": r"\phi=\frac d{4\pi r^2}\cos\theta,\ \psi=-\frac d{4\pi r}\sin^2\theta",
    "6.89": r"\psi=\tfrac12Ur^2\Big(1-\frac{a^3}{r^3}\Big)\sin^2\theta",
    "6.90": r"u_r=U\Big[1-\Big(\frac ar\Big)^3\Big]\cos\theta,\ u_\theta=-U\Big[1+\frac12\Big(\frac ar\Big)^3\Big]\sin\theta",
    "6.91": r"C_p=1-\frac94\sin^2\theta",
    "6.92": r"\phi=\Big(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3}\Big)\cdot\mathbf x",
    "6.93": r"\psi_{\text{sink}}=\frac k{4\pi}\int_0^a\cos\alpha\,d\xi",
    "6.94": r"\psi_{\text{sink}}=\frac{kR}{4\pi}\Big[\frac1{\sin\theta}-\frac1{\sin\alpha_1}\Big]=\frac k{4\pi}(r-r_1)",
    "6.95": r"\psi=-\frac Q{4\pi}\cos\theta+\frac Q{4\pi a}(r-r_1)+\tfrac12Ur^2\sin^2\theta",
    "6.96": r"\phi=-\frac{\mathbf d\cdot(\mathbf x-\mathbf x_s(t))}{4\pi\lvert\mathbf x-\mathbf x_s(t)\rvert^3}",
    "6.97": r"\phi=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi",
    "6.98": r"\mathbf F_s=-\int_{\text{sphere's surface}}(p-p_\infty)\,\mathbf n\,dA",
    "6.99": r"\Big[\frac{\partial\phi}{\partial t}+\frac12\lvert\nabla\phi\rvert^2+\frac p\rho\Big]_{\text{surface}}=\frac{p_\infty}\rho",
    "6.100": r"\frac{p_a-p_\infty}\rho=-\Big(\frac{\partial\phi}{\partial t}\Big)_a-\frac12\lvert\nabla\phi\rvert^2_a",
    "6.101": r"\frac{\partial\phi}{\partial t}=-\mathbf u\cdot\mathbf u_s-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\boldsymbol\xi\cdot\frac{d\mathbf u_s}{dt}",
    "6.102": r"\Big(\frac{\partial\phi}{\partial t}\Big)_a=-\mathbf u_a\cdot\mathbf u_s-\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}",
    "6.103": r"\nabla\phi=-\frac{a^3}2\Big[-\frac{3\boldsymbol\xi}{\lvert\boldsymbol\xi\rvert^5}\mathbf u_s\cdot\boldsymbol\xi+\frac{\mathbf u_s}{\lvert\boldsymbol\xi\rvert^3}\Big]",
    "6.104": r"\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s",
    "6.105": r"\frac{p_a-p_\infty}{\rho}=\frac12\lvert\mathbf u_s\rvert^2\Big(\frac94\frac{(\mathbf u_s\cdot\mathbf e_\xi)^2}{\lvert\mathbf u_s\rvert^2}-\frac54\Big)+\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}",
    "6.106": r"\Big(\frac{p_a-p_\infty}{\frac12\rho\lvert\mathbf u_s\rvert^2}\Big)_{\text{steady}}=1-\frac94\sin^2\theta_s",
    "6.107": r"\mathbf F_s=-\rho\frac a2\Big\lvert\frac{d\mathbf u_s}{dt}\Big\rvert\int_0^\pi\!\!\int_0^{2\pi}\cos\theta\,\mathbf e_\xi\,a^2\sin\theta\,d\varphi\,d\theta",
    "6.108": r"\mathbf F_s=-M\frac{d\mathbf u_s}{dt},\ M=\frac{2\pi a^3\rho}{3}",
    "6.109": r"\mathbf F_E=\Big(m+\frac{2\pi}{3}\rho a^3\Big)\frac{d\mathbf u_s}{dt}",
    # earlier chapters
    "2.30": r"\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA",
    "2.34": r"\int_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA=\oint_C\mathbf u\cdot d\mathbf x",
    "3.7": r"dx/u=dy/v=dz/w",
    "3.16": r"\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}",
    "3.17": r"\boldsymbol\omega=0",
    "3.18": r"\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA",
    "4.7": r"\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0",
    "4.10": r"\nabla\cdot\mathbf u=0",
    "4.12": r"\rho\mathbf u=\nabla\chi\times\nabla\psi",
    "4.17": r"\frac{d}{dt}\int_{V^*}\rho\mathbf u\,dV+\int_{A^*}\rho\mathbf u(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=\int_{V^*}\rho\mathbf g\,dV+\int_{A^*}\mathbf f\,dA",
    "4.38": r"\rho\frac{Du_j}{Dt}=-\frac{\partial p}{\partial x_j}+\rho g_j+\frac{\partial}{\partial x_i}\Big[\mu\Big(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\Big)+\big(\mu_v-\tfrac23\mu\big)\frac{\partial u_m}{\partial x_m}\delta_{ij}\Big]",
    "4.39b": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u",
    "4.40": r"\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega",
    "4.72": r"\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const everywhere}",
    "4.75": r"\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const}",
    "4.83": r"\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+\frac p\rho=\text{const}",
    "4.103": r"\mathrm{Re}=\frac{\rho Ul}{\mu}",
    "4.106": r"C_p\equiv\frac{p-p_\infty}{\tfrac12\rho U^2}",
    "4.111": r"M\equiv U/c",
    "5.2": r"u_\theta=\frac{\Gamma}{2\pi r}",
    "5.8": r"\frac{D\Gamma}{Dt}=0",
    "5.11": r"\frac{D\Gamma}{Dt}=\oint_C\Big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i",
}


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
                or nxt[:1] in "/–-)" or text[max(0, m.start() - 1):m.start()] in ("(", "–", "-")):
            return m.group(0)
        done.add(n)
        return f"({n}), ${EQ[n]}$"
    return _EQ_REF.sub(rep, text)


# ---------------------------------------------------------------------------------------------------------------------
# Part F reader: the derivations, word for word (the ch03–ch05 parser)
# ---------------------------------------------------------------------------------------------------------------------
def _join(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(x.strip() for x in lines)).strip()


def _abs_plain(s: str) -> str:
    r"""Part F writes absolute values as \lvert … \rvert everywhere; outside $…$ they must become plain |…| (in prose a
    TeX command would print as raw text)."""
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
    return s.strip(), ""


def _numbers_as_maths(s: str) -> str:
    """A Start/Result that names equations by number only ("(6.72) at the nodes") → the equations themselves."""
    if "$" in s:
        return s
    nums = re.findall(r"\(((?:\d)\.\d+[a-d]?)\)", s)
    if not nums:
        return s
    return " ".join(f"$({n})\\quad {EQ[n]}$" for n in nums if n in EQ) or s


_SELF = re.compile(r"(This is (?:the book's )?|This chain is the book's |is the book's |it is the book's |This is exactly )\(((?:\d)\.\d+[a-d]?)\)")


def _self_ref(text: str) -> str:
    """"This is (6.29)." inside a step's *why* points at the line just displayed: say so instead of repeating it."""
    return _SELF.sub(lambda m: f"{m.group(1)}({m.group(2)}) — the line above", text)


def part_f() -> dict[str, dict]:
    """Parse Part F of the design into {D01: dict(title, goal, start, plan, tools, assumptions, steps, result, check,
    meaning, traps)}."""
    text = (ROOT / "analysis" / "ch06_design.md").read_text(encoding="utf-8")
    part = text.split("## Part F", 1)[1].split("## Part G", 1)[0]
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
        f = {k: _abs_plain(_join(v)) for k, v in fields.items()}
        parsed = []
        for st in steps:
            s = _abs_plain(_join(st))
            bits = re.split(r"(?:^|\s·\s)\*(did|tex|why|plain|live|set|watch):\*\s*", s)
            d = {name: val.strip() for name, val in zip(bits[1::2], bits[2::2])}
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=_self_ref(d["why"]),
                               plain=_self_ref(d["plain"])))
        start_tex, start_plain = _split_words(f["Start"])
        res_tex, res_plain = _split_words(f["Result"])
        res_tex = re.sub(r"\s*\*?\((?:\d)\.\d+[a-d]?\)\*?\s*$", "", res_tex) if res_tex.count("$") >= 2 else res_tex
        plan = [p.strip() for p in re.split(r"\(\d+\)\s*", f["Plan"]) if p.strip()]
        tools = [t.strip() for t in f["Tools"].split(" · ") if t.strip()]
        out[key] = dict(title=title, goal=f["Goal"], start=(display_tex(_numbers_as_maths(start_tex)), start_plain),
                        plan=plan, tools=tools, assumptions=f.get("Assumptions", ""), steps=parsed,
                        result=(display_tex(_numbers_as_maths(res_tex)), res_plain), check=f.get("Check", ""),
                        meaning=f.get("What it means", ""), traps=f.get("Traps", ""))
    return out


PF = part_f()


def pf_sub(key: str, field: str, old: str, new: str) -> None:
    """Edit one Part F field after parsing (pointers to the cells that run a check, Part G errata); fails loudly if the
    design text changed. ``field`` = goal/check/meaning/traps/assumptions, ``tools``, or ``stepN.did|why|plain``."""
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
            print(f"  {i}. {s['did']} :: {s['tex'][:120]}  || why={len(s['why'].split())}w")
        print("  RESULT", d["result"][0][:150], "|", d["result"][1][:80])
        print("  CHECK", d["check"])
        print("  TOOLS", d["tools"])
    sys.exit(0)


# ---------------------------------------------------------------------------------------------------------------------
# Every "check" names a cell of this notebook that really runs it, and Part G errata override the design text.
# ---------------------------------------------------------------------------------------------------------------------
pf_sub("D01", "check", "has μ∇²u ≈ 1e-10 N/m³ by stencils; plane Poiseuille has 1 N/m³ (C01 code).",
       "has μ∇²u of order 1e-8 N/m³ by stencils (round-off, against an inertia term of about 10³ N/m³); plane Poiseuille "
       "has 1 N/m³ (the code cell after the tiny example of this block prints both).")
pf_sub("D03", "check", "(C02 tiny example). `delta_flux_check` gives", "(C02 tiny example). `delta_flux_check` (the code cell after that example) gives")
pf_sub("D04", "check", "`orthogonality_check` ≲ 1e-10 on the cylinder with circulation.",
       "`orthogonality_check` ≲ 1e-9 on the cylinder with circulation (the code cell after note N09).")
pf_sub("D05", "check", "(C03 code)", "(the code cell after the tiny example of C03)")
pf_sub("D06", "check", "`doublet_limit_error` has observed order 2.00.",
       "`doublet_limit_error` has observed order 2.0 (the code cell after that example).")
pf_sub("D08", "check", "(C05 code prints −0.4053)", "(the code cell after note N29 prints −0.4053)")
pf_sub("D10", "check", "θ = −9.16° and −170.84° ✓;", "θ = −9.16° and −170.84° ✓ (the code cell after the tiny example of C07);")
pf_sub("D11", "check", "(C07 code)", "(the from-scratch cell of C07)")
pf_sub("D12", "check", "(C08 code, ≤ 1e-14)", "(the code cell after note N39, ≤ 1e-14)")
pf_sub("D13", "check", "✓ (C08 code); the numerical route", "✓ (the code cell after the tiny example of C08); the numerical route")
pf_sub("D13", "check", "agrees to 1e-7.", "agrees to about seven digits (same cell).")
pf_sub("D15", "check", "(C09 figure)", "(the log–log figure of C09)")
pf_sub("D16", "check", "Cylinder: 24 N/m by (6.56) (C10 code) ✓.",
       "Cylinder: 24 N/m by (6.56) (`contour_force` in the code cell after the tiny example of C10) ✓.")
pf_sub("D17", "check", "(C10 code)", "(the code cell after the tiny example of C10)")
pf_sub("D18", "check", "(C10 code)", "(the same cell)")
pf_sub("D19", "check", "(C11 code)", "(the code cell after note N65)")
pf_sub("D21", "check", "for 10⁴ random ζ outside the circle in all four quadrants (tests); normal velocity on the ellipse ≲ 1e-12 and Blasius lift ρUΓ (C10, C11 code) ✓.",
       "for 10⁴ random ζ outside the circle in all four quadrants (in `tests/test_ch06.py`); normal velocity on the ellipse ≲ 1e-12 "
       "(the code cell after the tiny example of C11) and Blasius lift ρUΓ (the code cell after the tiny example of C10) ✓.")
pf_sub("D22", "check", "Observed order 2.00 ± 0.05 on sin(πx) sinh(πy) (R21 code) ✓.",
       "Observed order ≈ 2 on sin(πx) sinh(πy) (1.94, then 1.97 as the grid is halved again: the code cell under recap R21) ✓.")
pf_sub("D23", "check", "(C12 code)", "(the code cell after note N73)")
pf_sub("D24", "check", "(C13 code)", "(the code cell after note N75)")
pf_sub("D25", "check", "✓ (C13 code). Spherical Laplace (6.85) of both potentials = 0 (R30 code).",
       "✓ (the code cell after note N84). Spherical Laplace (6.85) of both potentials = 0 (the code cell under recap R30).")
pf_sub("D26", "check", "(C14 code)", "(the code cell after the tiny example of C14)")
# D17 steps 8–9 (lesson review): on a counterclockwise contour the velocity is parallel OR antiparallel to dz.
PF["D17"]["steps"][7]["tex"] = r"u+iv=\pm\lvert q\rvert e^{i\vartheta},\quad dz=\lvert dz\rvert e^{i\vartheta}"
PF["D17"]["steps"][7]["why"] = ("On the body the velocity is tangent to the surface (D05), so it points along dz or exactly "
                                "against it: the same direction angle ϑ up to a sign ± (against dz where the flow runs "
                                "clockwise round the body, e.g. the top of a cylinder in a stream from the left).")
PF["D17"]["steps"][7]["plain"] = "On the body, flow and surface point along the same line — forwards or backwards."
PF["D17"]["steps"][8]["tex"] = r"(u+iv)\,dz^*=\pm\lvert q\rvert\lvert dz\rvert=(u-iv)\,dz"
PF["D17"]["steps"][8]["why"] = ("The phases e^{iϑ} and e^{−iϑ} cancel either way round, and the same sign ± appears on both "
                                "sides, so it cancels too; this is (6.59).")
# D26 trap (lesson review M5): d(cot α) brings a minus sign, but so does −dξ — the limits keep their order.
pf_sub("D26", "traps", "d(cot α) = −dα/sin²α reverses the direction of the limits.",
       "The minus sign of d(cot α) = −dα/sin²α is cancelled by the one in −dξ: do not flip the limits "
       "(ξ = 0 → a becomes α = θ → α₁, in the same order).")
# Titles whose sense depends on a number (headings show the key equation instead of bare numbers).
PF["D10"]["title"] = "The cylinder with circulation: from its stream function to the stagnation points and the free point"
PF["D11"]["title"] = "The lift on the cylinder: from its surface pressure to L = ρUΓ, and D = 0"
# D09's start names (6.49) inside a display: write the equation there too.
assert r"\text{(from Im of (6.49))}" in PF["D09"]["start"][0]
PF["D09"]["start"] = (PF["D09"]["start"][0].replace(r"\text{(from Im of (6.49))}",
                      r"\text{(from the imaginary part of }w=\frac{d}{2\pi(z-z')}\ \text{(6.49))}"), PF["D09"]["start"][1])
# D17 step 12 and check (E5 builder note): cutting the body only matters where a singularity is crossed.
pf_sub("D17", "step12.why", "provided C′ encloses the body and no other singularity.",
       "provided C′ encloses every singularity of (dw/dz)² and crosses none. What counts is the singularities, not the "
       "outline: the cylinder's all sit at its centre, so even a circle that cuts the body gives ρUΓ; for the Zhukhovsky "
       "ellipse such a circle crosses the map's branch cut and gives a wrong number (R = 0.07 m: 8.13 N/m instead of 24).")
pf_sub("D17", "check", "the trapezoid on a circle converges exponentially.",
       "the trapezoid on a circle is exact for the cylinder with 4 or more points (its (dw/dz)² ends at z⁻⁴) and converges "
       "exponentially for bodies with an infinite Laurent series, such as the ellipse.")
# D23 (E7 builder measurement): on fine grids the residual is as small, relative to the true error, as the change.
pf_sub("D23", "step8.did", "Stop on the residual", "Test the residual, with a tight tolerance")
pf_sub("D23", "step8.why", "The change per sweep is about (1 − ρ) × error: when ρ ≈ 1 a small change hides a large error. The residual measures how badly (6.72) is violated.",
       "The residual measures how badly (6.72) is violated, so it is the right quantity to test. But when ρ ≈ 1 (fine "
       "grids) both the residual and the change per sweep are only about (1 − ρ) × the error: on the ×4 contraction grid "
       "after 120 Gauss–Seidel sweeps the change is 1.2e-3 and the residual 6.2e-4 while the true error is 5.6e-2. So the "
       "tolerance must be tight (10⁻¹⁰ here), far below the accuracy you want.")
pf_sub("D23", "step8.plain", 'Ask "is it solved?", not "did it stop moving?".',
       "Ask \"is (6.72) satisfied?\" — and demand it to many more digits than you need, because a slowly converging sweep "
       "looks solved long before it is.")
assert "and the residual says when to stop." in PF["D23"]["result"][1]
PF["D23"]["result"] = (PF["D23"]["result"][0], PF["D23"]["result"][1].replace(
    "and the residual says when to stop.", "and a tight tolerance on the residual says when to stop."))
pf_sub("D23", "traps", "Stopping on the change instead of the residual stops too early on fine grids.",
       "On fine grids a small residual or a small change does not mean a small error (both ≈ (1 − ρ) × error): use a "
       "tight tolerance.")
# D10 step 8 typo in the design ("r²/1").
pf_sub("D10", "step8.why", "multiply by r²/1", "multiply by r²")
# Part G3 (a), (b), (d): the axial strengths alternate in sign; only their moments converge; the matrix gets touchy.
pf_sub("D27", "check",
       "Rankine-oval target: the recovered k_n converge to a point source and sink as N grows; Σk_nΔξ → 0 (C14 code) ✓. "
       "Airship target: k_n reproduce the flat line sink of (6.95). Sphere: small residual but cond ≫ 10⁸ — fenced as qualitative.",
       "Rankine-oval target (the code cell after this derivation): the body error falls (2.5e-3, 1.7e-4, 1.8e-6 for N = 10, "
       "20, 40) and Σk_nΔξ → 0, but the individual k_n **alternate in sign** and grow (about ±5, ±10, ±19): only their "
       "sums and moments (the net source, the dipole) converge, not the bars themselves ✓. The condition number climbs "
       "from ~3×10³ to ~2×10¹², so N = 40 is the sensible ceiling here. Airship target: the fitted bars do not reproduce "
       "the point source and flat line sink of (6.95) — many strength patterns draw almost the same surface. Sphere: small "
       "residual but cond ≈ 10¹⁸ — qualitative only.")
pf_sub("D28", "check", "(C15 code: 1.6 and 0.0 for u_s = 2 e_z)", "(the code cell after note N94: 1.6 and 0.0 for u_s = 2 e_z)")
pf_sub("D29", "check", "Numbers (C15 code):", "Numbers (the code cell after the tiny example of C15):")
pf_sub("D29", "check", "Parity with `ch04.accelerating_sphere_pressure` for u_s ∥ du_s/dt;",
       "Parity with `ch04.accelerating_sphere_pressure` for u_s ∥ du_s/dt (same cell: 600 Pa at the front);")
pf_sub("D30", "check", "(C15 code)", "(the from-scratch cell of C15)")
pf_sub("D31", "check", "(C15 code)", "(the code cell after the tiny example of C15)")


# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch05.py)
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
    return re.sub(r"\s*`check_src` sketch.*$", "", s, flags=re.S)


def _title_no_numbers(title: str) -> str:
    """Drop bare equation numbers from a derivation title (the heading shows the key equation itself instead)."""
    num = r"\(\d\.\d+[a-d]?\)"
    t_ = re.sub(r"\s*(?:—\s*)?" + num + r"(?:\s*(?:→|–|,|and|or)\s*" + num + r")*", "", title)
    t_ = re.sub(r"\s+([,:)])", r"\1", t_)
    t_ = re.sub(r"\s*[—:,]\s*$", "", t_)
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
            if after[:1] in "–-/" or before[-1:] in "(–-":
                return m.group(0)                        # a range like (6.19)–(6.23)
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


_PROSE_SPLIT = re.compile(r"(```.*?```|\$\$.*?\$\$|\$[^$]+\$|`[^`\n]*`)", re.S)
_LABEL = re.compile(r"(?<![\w.§/−-])(\d\.\d+[a-d]?)(?![\w.]?\d)")
_NOT_EQ = re.compile(r"(?:§|Figs?\.|Figures?|Exercises?|Examples?|Ch\.|Chapters?|Sections?|chapter|section|version)\s*"
                     r"(?:\d\.\d+[a-d]?\s*(?:–|-|,|and|or)\s*)*$")


def _shown(n: str, cell: str) -> bool:
    """Is equation n written in this cell (its LaTeX, a display tagged with it, or maths right next to the number)?"""
    if EQ[n] in cell or ("\\text{(" + n + ")}") in cell or ("(" + n + ")}") in cell:
        return True
    e = re.escape(n)
    return bool(re.search(r"\$\s*\*?\((?:Eq\.\s*)?" + e + r"[\s,)*]", cell)
                or re.search(e + r"\)?\*?[,:]?\s*\$", cell))


def _bare_mentions(cell: str) -> list[tuple[int, int, str]]:
    """(start, end, label) of every equation number named in prose — with or without parentheses, in ranges
    "(6.21)–(6.22)", lists "(6.89, 6.92)", "Im of (6.49)", "(6.16 with …)", "from 4.72" — outside maths, code,
    headings and derivation step titles, and not a section, figure, exercise or example number."""
    out, pos = [], 0
    for k, seg in enumerate(_PROSE_SPLIT.split(cell)):
        if k % 2 == 0:
            for m in _LABEL.finditer(seg):
                n = m.group(1)
                if n not in EQ:
                    continue
                before = seg[:m.start()]
                line = before[before.rfind("\n") + 1:]
                if line.lstrip().startswith("#") or line.startswith("**Step ") or _NOT_EQ.search(before):
                    continue
                out.append((pos + m.start(), pos + m.end(), n))
        pos += len(seg)
    return out


def finalize_bare_numbers() -> int:
    """Second pass (lesson review M3): every equation number still named without its equation gets the equation
    written right after the mention (after the closing parenthesis when there is one)."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        units = c.source.split("\n---\n")               # a derivation's steps are separate units: each shows its own
        for u, src in enumerate(units):
            for s, e, n in reversed(_bare_mentions(src)):
                if _shown(n, src):
                    continue
                first = [x for x in _bare_mentions(src) if x[2] == n][0]
                if (s, e) != first[:2]:
                    continue                             # only the first mention in the unit gets the equation
                if src[e:e + 1] == ")":
                    e += 1
                src = src[:e] + f", ${EQ[n]}$" + src[e:]
            units[u] = src
        new = "\n---\n".join(units)
        if new != c.source:
            c.source = new
            changed += 1
    return changed


def self_check_prose() -> list[str]:
    """Lesson review M1/M3: no TeX command outside maths or code, and no equation named without being shown."""
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "markdown":
            continue
        prose = "".join(seg for k, seg in enumerate(_PROSE_SPLIT.split(c.source)) if k % 2 == 0)
        cmds = re.findall(r"\\[A-Za-z]+", prose)
        if cmds:
            bad.append(f"cell {i}: TeX outside maths {sorted(set(cmds))[:5]}: {c.source[:80]!r}")
        miss = sorted({n for u in c.source.split("\n---\n") for _s, _e, n in _bare_mentions(u) if not _shown(n, u)})
        if miss:
            bad.append(f"cell {i}: equations named but not shown {miss}: {c.source[:80]!r}")
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
Air over a wing, water past a bridge pier, the ocean round an island: far from walls and wakes these fluids move almost
without friction and without spin. There the flow obeys $\nabla\cdot\mathbf u=0$ and
$\rho\,D\mathbf u/Dt=-\nabla p$ *(Eq. 6.1)*, and two scalar functions — the stream function ψ and the velocity
potential φ — obey the **linear** Laplace equation $\nabla^2\psi=0$ *(Eq. 6.5)*, $\nabla^2\phi=0$ *(Eq. 6.12)*. So flows
can be *added*: a stream plus a source makes a blunt nose, a stream plus a doublet makes a cylinder, a vortex added to
it makes lift, $L=\rho U\Gamma$ *(Eq. 6.40)*. Complex numbers pack φ and ψ into one function $w=\phi+i\psi$ *(Eq. 6.42)*,
turn force calculations into a residue and let a map carry the circle's flow onto other shapes. A grid of averages
solves the same Laplace equation numerically, sources hidden on an axis draw airships, and a moving sphere shows the one
force an ideal fluid *does* exert on a body moving through it at rest far away: added mass,
$\mathbf F_E=\big(m+\tfrac{2\pi}{3}\rho a^3\big)\,d\mathbf u_s/dt$ *(Eq. 6.109)*.
""",
    roadmap=[
        r"§6.1 the ideal-flow equations $\nabla\cdot\mathbf u=0$, $\rho D\mathbf u/Dt=-\nabla p$ *(6.1)* and where to trust them (C01)",
        r"§6.2 ψ and φ: $\omega_z=-\nabla^2\psi$ *(6.4)*, Laplace, vortices and sources as point sources (C02); superposition — any streamline can be a wall *(6.16)* (C03)",
        r"§6.3 the element kit and the doublet $\phi=\lvert\mathbf d\rvert\cos\theta/2\pi r$ *(6.29)* (C04); the half-body *(6.31)* (C05); the cylinder and d'Alembert's paradox, $C_p=1-4\sin^2\theta$ *(6.35)* (C06); circulation and lift $L=\rho U\Gamma$ *(6.40)* (C07); images and the wall pressure under a passing vortex (C08)",
        r"§6.4 the complex potential $w=\phi+i\psi$ *(6.42)* and $dw/dz=u-iv$ *(6.45)* (C09)",
        r"§6.5 Blasius $D-iL=\frac{i\rho}{2}\oint(dw/dz)^2dz$ *(6.60)* and Kutta–Zhukhovsky lift for any body (C10)",
        r"§6.6 conformal maps and the Zhukhovsky map $z=\zeta+b^2/\zeta$ *(6.65)* (C11)",
        r"§6.7 Laplace on a grid, $\psi_{i,j}=\tfrac14[\dots]$ *(6.72)*, and Gauss–Seidel relaxation (C12)",
        r"§6.8 axisymmetric flow and the sphere *(6.90)* (C13); bodies from axial sources (C14)",
        r"§6.9 the accelerating sphere and added mass $M=2\pi\rho a^3/3$ *(6.108)* (C15)",
    ],
    prerequisites=[
        "continuity and the stream function (Ch. 4 §4.2–4.3)",
        "Bernoulli, steady and unsteady (Ch. 4 §4.9)",
        "vorticity, circulation and Kelvin's theorem (Ch. 3 §3.4, Ch. 5 §5.2)",
        "the point vortex and its image (Ch. 5 §5.7)",
        "the Laplacian, Gauss' theorem, the Dirac delta (Ch. 2)",
        "the cylinder flow in two frames (Ch. 3 §3.3)",
    ],
)
nb.explainer_index([
    ("superposition_sandbox", "Where does the body come from?",
     "C03 C04 C05: add elements; the stagnation streamline becomes a wall; closed only when the sources cancel"),
    ("cylinder_circulation_lift", "How does spin turn into lift?",
     "C06 C07: stagnation points slide, L = ρUΓ, drag stays zero"),
    ("vortex_wall_images", "What does a wall feel as an eddy passes?",
     "C08: image vortex, drift Γ/4πh, suction then over-pressure"),
    ("complex_potential_corners", "Why is the velocity u − iv?",
     "C09: Cauchy–Riemann as two difference quotients; calm vs violent corners"),
    ("blasius_kutta_contour", "Why doesn't lift depend on the shape?",
     "C10: the contour can move; only U × Γ/z survives"),
    ("conformal_joukowski", "How does a map carry a flow?",
     "C11: angles kept, circle → ellipse, the right square root"),
    ("laplace_relaxation", "How does a grid solve Laplace's equation?",
     "C12: average of neighbours, Jacobi vs Gauss–Seidel vs SOR, stop on the residual"),
    ("axial_singularity_bodies", "Given a shape, which sources draw it?",
     "C14: N unknown strengths, ψ = 0 at N points, one linear solve"),
    ("added_mass_sphere", "Why does an ideal fluid resist acceleration?",
     "C15: the speed pressure cancels, the acceleration pressure gives half the displaced mass"),
])
nb.setup()
nb.code(r"""
import numpy as np                                      # arrays (Ch. 1 primer P03)
import sympy as sp                                      # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                         # static figures (Ch. 1 primer P01)
import plotly.graph_objects as go                       # plotly figures that keep working on the web page (Ch. 1 P41)
from fluidpy import ch06_ideal_flow as ch06             # the tested chapter-6 module: every function cites its § and Eq.
from fluidpy.core import potential as pf                # new core module: superposable ideal-flow elements and bodies
from fluidpy.core import conformal as cm                # new core module: conformal maps (Zhukhovsky) and mapped flows
from fluidpy.core import laplace_solvers as ls          # new core module: Laplace on a masked grid, Jacobi / Gauss–Seidel / SOR
from fluidpy.core import panels as pn                   # new core module: axial singularities and source panels
from fluidpy import ch03_kinematics as ch03             # Ch. 3: the cylinder flow in two frames
from fluidpy import ch04_conservation_laws as ch04      # Ch. 4: incompressibility test, the accelerating sphere
from fluidpy import ch05_vorticity_dynamics as ch05     # Ch. 5: Kelvin's hypotheses, point vortices and images
from fluidpy.core.interact import slider_figure, live   # plotly slider figures (P17) and live widgets (P47)
from fluidpy.core.anim import animate, ffmpeg_path      # matplotlib animations (P16); show_animation came with the setup
from fluidpy.core.style import COLORS, savefig          # the house palette and a helper that saves PNGs to outputs/ch06
from tools.convergence import observed_order            # slope of log(error) vs log(step): the observed order (P13)
from scripts.ch06_drawings import flow_net, draw_body, pressure_arrows, separated_cp_band, circle   # drawing only, no physics
import warnings, logging                                # standard library: warning filters and log levels
warnings.filterwarnings("ignore", category=RuntimeWarning)             # 0/0 at a singular point is handled by NaN on purpose
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless font-substitution notes
ffmpeg_path()                                           # find ffmpeg once (sets matplotlib's animation.ffmpeg_path)
plt.rcParams["contour.negative_linestyle"] = "solid"    # negative ψ levels are streamlines too: draw them solid (φ is dashed on purpose)


def recolor(fig, colors, dashes=None):                  # give plotly slider traces the notebook's colours by trace name
    for tr in fig.data:                                 # every trace of every slider step
        if tr.name in colors:                           # a name we assigned a colour to
            tr.line.color = colors[tr.name]             # same colour meaning as in the matplotlib figures
            tr.marker.color = colors[tr.name]           # markers too (for "markers" traces)
        if dashes and tr.name in dashes:                # optional dash pattern ("dash", "dot")
            tr.line.dash = dashes[tr.name]
    return fig                                          # the same figure, restyled


def step_titles(fig, titles):                           # a different title for every slider position (numbers per step)
    for step, text in zip(fig.layout.sliders[0].steps, titles):   # one slider step per title
        step.args = [step.args[0], {"title.text": text}]          # the step now also rewrites the figure title
    fig.layout.title.text = titles[fig.layout.sliders[0].active]  # the title of the starting position
    return fig                                          # the same figure


print(len([n for n in dir(ch06) if not n.startswith("_")]), "public names in fluidpy.ch06_ideal_flow")   # the toolbox
""", explain="""
1. Numerical, symbolic and plotting libraries (all primed in Ch. 1–2).
2. `ch06` is the chapter module; it re-exports the four new core modules of this chapter — `pf` (elements such as
   `Uniform`, `Source`, `Vortex`, `Doublet`, `Corner`, their sum `Flow`, the bodies `half_body`, `cylinder`, `sphere`
   and the force tools), `cm` (conformal maps), `ls` (grid solvers) and `pn` (singularity distributions). We call them
   by these short names; every function is tested in `tests/test_ch06.py`.
3. The chapter modules of Ch. 3–5 are used for recaps and parity checks.
4. `scripts/ch06_drawings.py` only draws: flow nets (ψ solid, φ dashed), bodies, pressure arrows and the *qualitative*
   separated-pressure band of C06.
5. `recolor` and `step_titles` only restyle plotly figures (colours by meaning, a title per slider position).
""")
nb.md(r"""
### ⚠️ Conventions in this chapter (read once; each is repeated where it bites)

| Topic | This notebook and its code | Where the book differs |
|---|---|---|
| circulation Γ [m²/s] | **counterclockwise positive**: a vortex $\psi=-\frac{\Gamma}{2\pi}\ln r$ *(6.8)* turns the fluid at $u_\theta=+\Gamma/2\pi r$; `pf.Vortex(Gamma)` | the cylinder with circulation (6.36)–(6.40), (6.52), (6.61)–(6.62), (6.68) and Example 6.1 use a **clockwise** Γ; body functions accept `Gamma_cw=` (the book's Γ) or `Gamma_ccw=` (ours), $\Gamma_{ccw}=-\Gamma_{cw}$ |
| 2-D doublet | the dipole vector $\mathbf d=\sum_i\mathbf x_im_i$ points **from the sink to the source**; the cylinder needs $\mathbf d=-2\pi Ua^2\mathbf e_x$ (upstream) | (6.49)'s scalar d means the dipole $-d\,\mathbf e_x$ (`pf.Doublet.from_book_scalar`) |
| angles | θ measured from +x (downstream): a body's nose in a stream from the left is at θ = π | the measured cylinder pressure is plotted against the angle from the front, 180° − θ |
| §6.8 axes | z is the **horizontal** symmetry axis, along the stream | — |
| §6.9 letters | w is the z-velocity, not the complex potential; M is the added mass, not the Mach number | — |
| forces | D (along the stream) and L (across it) are forces **on the body** per unit depth [N/m]; B is a span, not Ch. 4's Bernoulli function | (6.54)'s F acts on the fluid |
| units | plane ψ, φ, Γ, m [m²/s]; Stokes ψ and 3-D source Q [m³/s]; 2-D dipole [m³/s], 3-D dipole [m⁴/s] | — |

**Numbers for the Γ convention** (used again in C07): a cylinder of radius a = 0.1 m in air (ρ = 1.2 kg/m³) and a
stream U = 10 m/s with the book's Γ_cw = 2 m²/s — i.e. our Γ_ccw = −2 m²/s — has its stagnation points **below** the
axis at −9.16° and −170.84° and feels the lift L = ρUΓ_cw = +24 N/m. We compute in the project convention and show the
book's alongside wherever they differ.

**Colours** (one meaning each, also in the explainers): uniform stream **blue** · sources **teal** · sinks **rose** ·
vortices **amber** · doublets **purple** · body and dividing streamline **black** · over-pressure **orange**, suction
**blue** · steady pressure part **teal**, acceleration part **orange** · lift **amber**, drag **grey** · ψ contours solid,
φ contours dashed · references and ghosts **grey dashed** · a "wrong variant" (numpy's principal root, a printed
coefficient) **rose dashed**.
""")
nb.md(r"""
**Book slips we correct in this notebook** (each is explained where it is used):

| The book prints | We use |
|---|---|
| (6.61)'s squared series with the $1/z^2$ coefficient $\big(\tfrac{Ud}{\pi}-\tfrac{\Gamma^2}{4\pi^2}\big)$ and an extra outer square | the coefficient is $-\big(\tfrac{Ud}{\pi}+\tfrac{\Gamma^2}{4\pi^2}\big)$; only the $1/z$ term matters (D18) |
| the middle bracket of (6.104), $\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s$: $-\tfrac1{a^3}\mathbf u_s$ | $+\tfrac1{a^3}\mathbf u_s$, which that final form needs (D29) |
| (6.108), $\mathbf F_s=-M\frac{d\mathbf u_s}{dt},\ M=\frac{2\pi a^3\rho}{3}$, keeps $d\varphi$ after the φ-integration | the φ-integral has become its factor 2π (D30) |
| §6.3 "differentiation of (6.8)" for the source velocities | (6.15), $\phi=\frac m{2\pi}\ln r$ (note N24) |
| §6.4 "(6.5) and (6.12), respectively" for φ and ψ; "(6.43) ensures the equality" | the other way round; the Cauchy–Riemann conditions (6.44) (note N43) |
| §6.5 "Section 3"; §6.6 "(see Figure 6.5)" | §6.3; Figure 6.3 (notes N53, N66) |
| §6.7 "first-order central differences" | first-*derivative* differences, second-order accurate (recap R19, D22) |
| Example 6.2's FORTRAN loop over I, and the inlet Δψ in m² | the loop should run over J; Δψ is in m²/s (note N74) |
""")
nb.md(r"""
**Where this chapter is used later**

| Result here | Used in |
|---|---|
| ideal outer flow and $C_p=1-\lvert\mathbf u\rvert^2/U^2$ (6.32) | the pressure that drives boundary layers, wedge flows $w=Az^n$ (6.46), separation (Ch. 9) |
| velocity potentials, images at a wall | water waves over a flat bottom (Ch. 7) |
| Laplace on a grid, relaxation, panels | elliptic solvers and panel methods (Ch. 10) |
| $\omega_z=-\nabla^2\psi$ (6.4), deformation (corner) flow, images, Poisson in a closed basin | stream-function inversion, fronts, coastal walls (Ch. 13) |
| Kutta–Zhukhovsky $L=\rho U\Gamma$ (6.62), conformal maps | airfoils, the Kutta condition, induced drag (Ch. 14) |
| added mass $M=2\pi\rho a^3/3$ (6.108) | bubbles, fish and structures (Ch. 16) |

*Climate hook:* the atmosphere and ocean are rotating and stratified, so they are almost never irrotational — but the
*tools* of this chapter carry over directly: recovering ψ from vorticity (a Poisson problem, C02 and C12), the strain
field of a corner flow that sharpens fronts (C04), walls as mirror images (C08).
""")
nb.md(r"""
> 🔁 **Tools from earlier chapters used in this one** (one line each where first used; the full primers are in
> `knowledge/primers.md`): partial derivatives (P25), Taylor series (P26, P98), definite integrals (P27), the net
> pressure force −∮p n dA (P28), `lambda` (P29), the trapezoid rule (P37), the product rule (P38), sympy (P40),
> i² = −1 and Euler's formula (P45), `np.where` (P46), the chain rule (P49, P91), `np.linalg.solve` (P57), orders of
> smallness (P68), Vieta's formulas (P71), `np.arctan2` (P70), the directional derivative (P75), `np.meshgrid` (P76),
> broadcasting (P77), contour and streamline plots (P78), eigenvalues (P80), the complex conjugate (P81), `quad` (P87),
> polar and spherical unit vectors (P88, P105), parametric curves (P92), substitution in an integral (P106), `brentq`
> (P108), dataclasses (P111), momentum flux (P114), Schwarz's theorem (P121), the curl of a curl (P122), the periodic
> trapezoid rule (P136), Poisson and Green's functions (P139), the gradient of 1/distance (P140), the FFT (P142),
> Gauss–Legendre quadrature (P143), improper integrals (P144), `animate` (P16), `slider_figure` (P17), `show_viz` (P18),
> `solve_ivp` (P31, P94), log–log slopes (P13), `assert np.allclose` (P15). New primers of this chapter: P149–P164.
""")

# =====================================================================================================================
# A.1 §6.1 — R01, C01 (N01–N04, D01)
# =====================================================================================================================
nb.section("6.1", "Relevance of Irrotational Constant-Density Flow Theory", intro=r"""
**What is this section about?** Real fluids are viscous, so why study a theory with no friction? Because at high
Reynolds number viscosity only matters in thin layers next to walls and in wakes; everywhere else the flow keeps no spin
and obeys a much simpler pair of equations. This section says what those equations are, what they drop, and where they
may be trusted.
""")
nb.recap("R01", "Why the fluid stays irrotational — Kelvin", r"""
In a fluid of constant density pushed only by conservative forces, and away from viscous regions, the circulation of
every material loop is constant, $D\Gamma/Dt=0$ *(Eq. 5.8)* — so fluid that starts without vorticity (a uniform stream
far upstream) stays without it until it enters a boundary layer, a wake or a separated region. Constant density needs a
low Mach number, $M\equiv U/c\ll1$ *(Eq. 4.111)*, and no baroclinic torque $\nabla\rho\times\nabla p/\rho^2$ from a
density field (Ch. 5 §5.2 lists the restrictions). Below, `ch05.kelvin_hypotheses_text` says which hypothesis fails and
`ch04.is_incompressible_regime(U)` tests whether a speed U [m/s] is slow enough (U/c below 0.3 in air).
""", where="Ch. 5 §5.2")
nb.code(r"""
print(ch05.kelvin_hypotheses_text(barotropic=False))    # a stratified (baroclinic) fluid: Kelvin's theorem fails
print(ch04.is_incompressible_regime(10.0))              # a car at 10 m/s in air: M = 10/343 ≈ 0.03 → constant density is fine
""")

core("C01", "The ideal-flow equations — and where to trust them", r"""
Water and air are viscous. How can a theory with no friction describe them at all — and where does it fail?
""", eqs=("6.1",))
problem(r"""
Put your hand out of a car window at 10 m/s: the air feels the car only through a skin of slowed air about a millimetre
thick on its surface; a few centimetres away the air moves as if it had no viscosity at all. Engineers design wings,
hulls and cars from that outer flow; oceanographers use it for flow round islands and seamounts. We want the equations
of that outer flow and an honest list of where they fail.
""")
idea("""
far from walls:    no spin (w = 0)  -->  viscous force mu*lap(u) = -mu*curl(w) = 0  -->  Euler + div u = 0  = (6.1)
near the wall:     thin boundary layer, w != 0, friction matters                     -->  Ch. 9
behind a blunt body: separated wake, w != 0 everywhere                               -->  ideal flow fails
""", r"""
(Here w stands for the vorticity ω.) **Irrotational does not mean inviscid — it means the viscous forces cancel.**
""")
remind([
    ("partial derivative ∂/∂x", r"$\partial f/\partial x$ changes x and holds every other variable fixed (Ch. 1 P25)."),
    ("product rule for a divergence", r"$\nabla\cdot(\rho\mathbf u)=\mathbf u\cdot\nabla\rho+\rho\nabla\cdot\mathbf u$ (Ch. 4 P113)."),
    ("curl of a curl identity", r"$\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u$ for any smooth field (Ch. 4 P122)."),
    ("boundary conditions", "what the flow must do on the edges of its region — no flow through a wall, no slip, a given stream far away (Ch. 1 P20)."),
], lead="the vector identities D01 uses; the hydrostatic split p = p_h + p′ is Ch. 4 §4.9")
D("D01", ref="6.1")
note("N02", "Only one boundary condition survives", r"""
$\rho\,D\mathbf u/Dt=-\nabla p$ contains only first derivatives of u, so a wall can impose one condition — no flow
through it:
""", equation=r"\mathbf n\cdot\mathbf u=\mathbf n\cdot\mathbf U_s\quad(\text{kept}),\qquad \mathbf t\cdot\mathbf u=\mathbf t\cdot\mathbf U_s\quad(\text{dropped})")
nb.md(r"""
The fluid slides along the wall. D01 step 8 showed why: the term that needed a second condition, μ∇²u, is zero. The
real fluid obeys no-slip, which is why a boundary layer must exist (Ch. 9).
""")
nb.worked_example("a corner flow that is viscous but feels no viscous force", r"""
Take $\phi=A(x^2-y^2)$ with A = 1 s⁻¹ (a stagnation flow, met again in C04).

1. $u=\partial\phi/\partial x=2x$, $v=\partial\phi/\partial y=-2y$ [m/s].
2. $\nabla\cdot\mathbf u=2-2=0$ ✓ incompressible.
3. $\omega_z=\partial v/\partial x-\partial u/\partial y=0-0=0$ ✓ irrotational.
4. $\nabla^2u=\partial^2(2x)/\partial x^2+\partial^2(2x)/\partial y^2=0$, likewise $\nabla^2v=0$: no net viscous force.
5. But the element is being stretched: the viscous normal stress $2\mu\,\partial u/\partial x=2\times10^{-3}\times2=0.004$ Pa
   in water is **not** zero — it is the same on the two opposite faces of each element, so its net push cancels.
""")
remind([
    ("central differences at a point", r"$f'(x)\approx[f(x+h)-f(x-h)]/2h$ and $f''(x)\approx[f(x+h)-2f(x)+f(x-h)]/h^2$ (Ch. 1 P21)."),
    ("Python dictionaries (fluidpy results)", "fluidpy returns several results at once as a dictionary: `r[\"vorticity\"]` (Ch. 1 P23)."),
    ("f-strings and printing", "`f\"{x:.2e}\"` prints a number in scientific notation inside a sentence (Ch. 1 P04)."),
    ("functions as arguments and lambda", "`lambda X, Y: …` is a one-line function; fluidpy takes fields as functions (Ch. 1 P29)."),
    ("assert np.allclose", "`assert np.allclose(a, b)` stops the notebook if two results differ beyond rounding (Ch. 1 P15)."),
], lead="`{k: f(v) for k, v in d.items()}` is a *dict comprehension*: it builds a new dictionary by applying f to every value of d")
nb.code(r"""
cases = [("corner", (0.2, 0.15)), ("cylinder", (0.2, 0.15)), ("poiseuille", (0.2, 0.005))]   # three flows and a probe point [m]
for name, pt in cases:                                   # the Poiseuille flow is probed inside its 2 cm gap
    r = ch06.ideal_flow_residuals(name, x=pt)            # residuals of ∇·u = 0, of Euler, the viscous force μ∇²u, the vorticity
    print(f"{name:10s}", {k: f"{float(np.max(np.abs(v))):.2e}" for k, v in r.items()})   # size of each term
""", explain=r"""
1. Three velocity fields are built in as presets: the corner flow $u=(2Ax,-2Ay)$ of the tiny example, the cylinder
   flow (U = 1 m/s, a = 0.1 m, C06) with its Bernoulli pressure, and **plane Poiseuille flow** — the viscous flow
   between two plates driven by a pressure gradient G = −dp/dx = 1 Pa/m (Ch. 4 §4.6), here in a 2 cm gap.
2. `ideal_flow_residuals` evaluates, by central differences at the point: the continuity residual ∇·u [1/s], the Euler
   residual ρDu/Dt + ∇p [N/m³], the viscous force μ∇²u [N/m³] (water, μ = 10⁻³ Pa s), the vorticity ω_z [1/s], and for
   scale the size of the inertia term ρDu/Dt.
3. The two potential flows pass: every residual is round-off or stencil error (at most ~10⁻⁴ N/m³ against an inertia
   term of ~10³ N/m³). The Poiseuille flow is rotational (ω_z = Gy/μ = 5 s⁻¹ at y = 5 mm) and carries a net viscous force
   of exactly G = 1 N/m³ — the force that balances its pressure gradient, which is also why the Euler equation alone
   fails for it by the same 1 N/m³.
""")
nb.md("**From scratch — the viscous force by our own five-point sum** (the second differences of Ch. 1 P21 in x and y):")
nb.check_agree(r"""
h, mu = 1e-3, 1e-3                                        # stencil step [m] and water's viscosity [Pa s]
x0, y0 = 0.2, 0.15                                        # the probe point [m]
u = lambda X, Y: pf.cylinder(1.0, 0.1).velocity(X, Y)[0]  # x-velocity of the cylinder flow, U = 1 m/s, a = 0.1 m
lap = (u(x0 + h, y0) + u(x0 - h, y0) + u(x0, y0 + h) + u(x0, y0 - h) - 4*u(x0, y0))/h**2   # five-point ∇²u [1/(m s)]
print(f"our μ∇²u = {mu*lap:.2e} N/m³")                   # ≈ 0 up to the stencil error
assert np.allclose(mu*lap, ch06.ideal_flow_residuals("cylinder", x=(x0, y0))["viscous_force"][0], atol=1e-6)   # both ≈ 0
G, h0 = 1.0, 0.01                                         # Poiseuille: driving gradient [Pa/m], half gap [m]
up = lambda Y: G/(2*mu)*(h0**2 - Y**2)                    # the parabolic Poiseuille profile u(y) [m/s]
lap_p = (up(0.005 + h) - 2*up(0.005) + up(0.005 - h))/h**2   # only ∂²u/∂y² is non-zero [1/(m s)]
print(f"Poiseuille μ∇²u = {mu*lap_p:.4f} N/m³")           # −G: the viscous force balances the pressure push
assert np.isclose(mu*lap_p, -G)                           # exact for a parabola
""")
nb.md("Our own five-point sums give the library's viscous force: zero for the potential flow, −G for Poiseuille.")
note("N01", "Where the ideal-flow equations work", r"""
— with a large Reynolds number $\mathrm{Re}=\rho Ul/\mu$ *(Eq. 4.103)* viscosity and spin stay in thin layers.

| ideal flow predicts | it does not predict |
|---|---|
| the velocity away from walls | skin friction |
| pressure forces, through a thin attached boundary layer | dissipation |
| streamline patterns that minimise form drag | pipe and duct flow |
| unsteady inertia (added mass, C15) | wakes, turbulence, low-Re flow |

**Number:** a car, U = 10 m/s, L = 1 m, air ν = 1.5×10⁻⁵ m²/s: Re = UL/ν ≈ 6.7×10⁵; a laminar boundary layer is about
L/√Re ≈ 1.2 mm thick.
""")
nb.code(r"""
import re                                               # regular expressions, to trim the threshold from each reason
for case in [dict(Re=6.7e5), dict(Re=50), dict(Re=6.7e5, M=0.5), dict(Re=6.7e5, baroclinic=True), dict(Re=6.7e5, region="wake")]:
    res_ = ch06.ideal_flow_applicability(**case)          # the decision of N01 for this case (a dictionary)
    why = "; ".join(re.sub(r"^[^:]*:\s*", "", s) for s in res_["reasons"])   # the reason in words, without the numbers
    print(case, "→", "applies" if res_["ok"] else "does not apply: " + why)
""", explain="""
The five cases: the car's outer flow (applies); a creeping flow at Re = 50 (boundary layers are not thin); a fast
flow at M = 0.5 (the density changes); a stratified fluid (Kelvin's theorem fails, vorticity is created inside); a wake
(the fluid has passed through boundary layers and carries vorticity).
""")
remind([
    ("np.meshgrid and the grid layout", "`X, Y = np.meshgrid(x, y)` gives every (x, y) pair; rows run along y, columns along x (Ch. 2 P76)."),
    ("contour and streamline plots", "`ax.contour(X, Y, psi)` draws level lines — for ψ they are streamlines (Ch. 2 P78)."),
    ("logarithmic colour and axes", "a colour or axis in log₁₀ shows numbers that span many powers of ten; −9 means 10⁻⁹ (Ch. 1 P13)."),
], lead="used by the figure below; `np.linalg.norm(v)` is the length of a vector (Ch. 2 P67); `ax.pcolormesh(x, y, Z)` "
         "colours every grid cell by its value Z (a heat map)")
nb.figure(r"""
nx, ny = (61, 41) if not FAST else (41, 27)               # probe grid (fewer points in FAST mode)
xs, ys = np.linspace(-0.3, 0.3, nx), np.linspace(-0.2, 0.2, ny)   # region round the cylinder [m]
VF = np.full((ny, nx), np.nan)                            # |μ∇²u| at each probe [N/m³]; NaN inside the body
for j, y in enumerate(ys):                                # loop over rows (y) …
    for i, x in enumerate(xs):                            # … and columns (x)
        if np.hypot(x, y) > 0.105:                        # outside the cylinder (a = 0.1 m) plus a stencil margin
            VF[j, i] = np.linalg.norm(ch06.ideal_flow_residuals("cylinder", x=(x, y))["viscous_force"])
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6))     # two panels: potential flow | Poiseuille flow
im = a1.pcolormesh(xs, ys, np.log10(VF + 1e-16), cmap="Greys", vmin=-12, vmax=1, shading="auto")   # log₁₀ of the force
flow_net(a1, pf.cylinder(1.0, 0.1), xlim=(-0.3, 0.3), ylim=(-0.2, 0.2), n=161, n_levels=15,
         body=circle(0.1))                                # streamlines and the body on top
a1.set_title("cylinder flow: μ∇²u at round-off level", fontsize=10)
fig.colorbar(im, ax=a1, label=r"$\log_{10}|\mu\nabla^2\mathbf{u}|$ [N/m$^3$]")
yp = np.linspace(-0.01, 0.01, 41)                         # across the 2 cm Poiseuille gap [m]
up = 1.0/(2e-3)*(0.01**2 - yp**2)                         # u(y) = G(h0² − y²)/2μ [m/s]
vfp = [np.linalg.norm(ch06.ideal_flow_residuals("poiseuille", x=(0.0, y))["viscous_force"]) for y in yp]   # |μ∇²u| [N/m³]
a2.pcolormesh([0, 0.06], yp, np.log10(np.array(vfp))[:, None]*np.ones((1, 2)), cmap="Greys", vmin=-12, vmax=1, shading="auto")
a2.quiver(0*yp[::4], yp[::4], up[::4], 0*yp[::4], color=COLORS["blue"], angles="xy", scale_units="xy", scale=1.0, width=0.006)
a2.plot(up, yp, color=COLORS["blue"], lw=2)               # the parabolic profile
a2.set_xlim(0, 0.06); a2.set_xlabel("u [m/s]"); a2.set_ylabel("y [m]")   # axes with units
a2.set_title("Poiseuille: |μ∇²u| = G = 1 N/m³ everywhere", fontsize=10)
fig.suptitle("No spin, no net viscous force", fontweight="bold")
savefig(fig, "ch06", "c01_viscous_force"); plt.show()    # save the PNG to outputs/ch06, then draw
""", see="Left, streamlines round the cylinder over a grey background: log₁₀ of the viscous force is round-off (about "
         "−12 … −8) far away and ≈ −6 … −4 right next to the body (the stencil's truncation error there, still ~10⁸ times "
         "smaller than the inertia term); right, the parabolic Poiseuille profile over a uniformly dark strip (log₁₀ 1 = 0).",
    read=r"A potential flow has $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega=0$ everywhere, although the fluid "
         "is viscous; Poiseuille flow is all spin, and its viscous force is the same G at every height.",
    change="…μ multiplied by 100: the right panel's force stays G (the profile flattens 100 times so that μ∇²u still "
           "balances G), the left panel stays at round-off — no μ can create a net viscous force where ω = 0.")
note("N03", "Separation breaks the limit", r"""
As Re grows, an attached boundary layer thins toward zero and real flow approaches ideal flow. But where the layer
separates (behind a cylinder, off a sharp corner) a wake of spinning fluid forms whose size does not shrink with Re: the
limit μ → 0 of real flow is then *not* the flow with μ = 0. Ch. 9 explains separation; C06 shows its pressure signature
on a cylinder.
""")
note("N04", "In and out", r"""
Excluded: fluids of varying density, high subsonic and supersonic speeds, boundary layers, wakes, flows inside pipes,
any region where elements spin. Included: flight (Ch. 14), water waves (Ch. 7), flow round vehicles and structures.
""")
confusion(r"""
"ideal flow = inviscid fluid". The fluid can be as viscous as honey; what matters is ω = 0. D01 shows the viscous
*stress* survives (step 5 of the tiny example: 0.004 Pa); only its net force vanishes. What ideal flow gives up is the
*no-slip* condition.
""")
whatif(r"""
…the fluid were stratified (density varying with height)? Then the baroclinic torque $\nabla\rho\times\nabla p/\rho^2$
(Ch. 5 §5.2) creates vorticity in the interior and Kelvin's theorem fails: the outer flow is no longer irrotational. That
is the everyday case in the atmosphere and ocean — which is why Ch. 13 needs potential *vorticity*, not a potential.
""")

# =====================================================================================================================
# A.2 §6.2 — R02 R03 R04 R06 R07 R08 R09, C02 (N05–N15, N18, N19, D02–D04), R05, C03 (N16, N17, D05)
# =====================================================================================================================
nb.section("6.2", "Two-Dimensional Stream Function and Velocity Potential", intro=r"""
**What is this section about?** In a plane flow two scalar functions carry the whole velocity field: the stream
function ψ (its contours are streamlines) and the velocity potential φ (its gradient is the velocity). Irrotational flow
turns both into solutions of Laplace's equation, which is *linear* — and that is the key to the whole chapter.
""")
nb.recap("R02", "Plane continuity", r"""
For constant density in two dimensions, $\partial u/\partial x+\partial v/\partial y=0$ *(Eq. 6.2)* — the plane form of
$\nabla\cdot\mathbf u=0$ *(Eq. 4.10)*.
""", where="Ch. 4 §4.2")
nb.recap("R03", "The stream function", r"""
Define $u\equiv\partial\psi/\partial y$, $v\equiv-\partial\psi/\partial x$ *(Eq. 6.3)*; then
$\partial u/\partial x+\partial v/\partial y=0$ holds for any smooth ψ, because mixed partial derivatives commute. Along a
curve ψ = const the **total differential** vanishes, $0=d\psi=\frac{\partial\psi}{\partial x}dx+\frac{\partial\psi}{\partial y}dy
=-v\,dx+u\,dy$, so $(dy/dx)_\psi=v/u$ — the slope of a streamline. The difference of ψ between two streamlines is the
volume flow between them per unit depth [m²/s]. Ch. 4's `velocity_from_streamfunction_2d` differentiates a ψ function
numerically.
""", where="Ch. 4 §4.3")
nb.code(r"""
from fluidpy.core import streamfunction as sf           # Ch. 4's stream-function tools
print(sf.velocity_from_streamfunction_2d(lambda x, y: 2.0*y, 0.3, 0.7))   # ψ = Uy with U = 2 m/s → (u, v) = (2, 0) m/s
""")
nb.recap("R04", "Plane irrotationality", r"""
A 2-D flow has no spin when $\partial v/\partial x-\partial u/\partial y=0$ *(Eq. 6.9)* — the third component of the
vorticity, $\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}$ *(Eq. 3.16)*, set to zero. When
$\boldsymbol\omega=0$ *(Eq. 3.17)* the velocity is a gradient, $\mathbf u=\nabla\phi$ — and in a **simply connected**
region (one without holes, where every loop can be shrunk to a point inside the fluid) φ is single-valued.
""", where="Ch. 3 §3.4")
nb.recap("R06", "Polar continuity", r"""
$\frac1r\frac{\partial}{\partial r}(ru_r)+\frac1r\frac{\partial u_\theta}{\partial\theta}=0$ *(Eq. 6.19)* — the divergence
in polar coordinates (Appendix B).
""", where="App. B; operators recapped in Ch. 4")
nb.recap("R07", "Polar irrotationality", r"""
$\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}=0$ *(Eq. 6.20)*.
""", where="Ch. 3 §3.4")
nb.recap("R08", "Polar Laplacian of ψ", r"""
$\nabla^2\psi=\frac1r\frac{\partial}{\partial r}\Big(r\frac{\partial\psi}{\partial r}\Big)+\frac1{r^2}\frac{\partial^2\psi}{\partial\theta^2}=0$
*(Eq. 6.23a)* — used in D03.
""", where="App. B; operators recapped in Ch. 4")
nb.recap("R09", "Polar Laplacian of φ", r"""
The same operator on φ, $\nabla^2\phi=\frac1r\frac{\partial}{\partial r}\Big(r\frac{\partial\phi}{\partial r}\Big)+\frac1{r^2}\frac{\partial^2\phi}{\partial\theta^2}=0$
*(Eq. 6.23b)*. Below, sympy applies this polar Laplacian to four functions we will meet as flows.
""", where="App. B; operators recapped in Ch. 4")
remind([
    ("sympy diff, simplify, symbols", "`sp.symbols` makes symbols, `sp.diff` differentiates, `sp.simplify` tidies an expression (Ch. 1 P40)."),
])
nb.code(r"""
r, th = sp.symbols("r theta", positive=True)            # polar coordinates as symbols (r > 0)
lap = lambda f: sp.diff(r*sp.diff(f, r), r)/r + sp.diff(f, th, 2)/r**2   # the polar Laplacian of (6.23a, b)
print([sp.simplify(lap(f)) for f in (r*sp.sin(th), sp.log(r), sp.cos(th)/r, r**2*sp.sin(2*th))])   # [0, 0, 0, 0]
""", explain=r"""
The four functions are the ψ of a uniform stream ($r\sin\theta=y$), the ψ or φ of a vortex or source (ln r), the φ of a
doublet (cos θ/r, C04) and the ψ of a corner flow ($r^2\sin2\theta=2xy$, C04). All four have zero Laplacian for r > 0:
they are **harmonic**.
""")

# ---- C02 ------------------------------------------------------------------------------------------------------------
core("C02", "Vorticity is minus the Laplacian of ψ — and vortices and sources are point sources of Laplace's equation", r"""
Every irrotational plane flow obeys one equation. Which — and what makes a whirlpool or a spring different from the flow
around them?
""", eqs=("6.4",))
problem(r"""
Look down at a bath draining (a whirlpool) and at a garden sprinkler lying flat (a spring). Around both the water moves
without spinning; yet one stirs the water round and the other pushes it out. We want one equation for all such flows,
and a way to say exactly how strong the whirl or the spring is.
""")
nb.md(r"""
#### The idea

| | stream function ψ | potential φ |
|---|---|---|
| defined by | $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$ | $u=\partial\phi/\partial x$, $v=\partial\phi/\partial y$ |
| automatic | continuity | irrotationality |
| the other law gives | $\nabla^2\psi=-\omega_z$ | $\nabla^2\phi=q$ (sources) |
| point singularity | vortex: $-\Gamma\delta$ | source: $m\delta$ |

**Away from a few special points both ψ and φ are harmonic: ∇² = 0. The special points carry all the information.**
""")
remind([
    ("natural logarithm, d(ln r)/dr = 1/r", r"$\ln$ is the inverse of $e^x$, and $\frac{d}{dr}\ln r=\frac1r$ (Ch. 1 P36)."),
], lead="the Laplacian itself is Ch. 2 §2.9; the Dirac delta as a concentrated total is Ch. 2 and Ch. 5 P139")
D("D02", ref="6.4")
note("N05", "Laplace for ψ", "In irrotational flow the vorticity equation becomes Laplace's equation for ψ; its "
     "solutions are called *harmonic*:", equation=EQ["6.5"], ref="6.5")
note("N06", "A point vortex", r"""
A point (ideal) vortex of strength Γ at (x′, y′) is all its vorticity squeezed into one point — a delta function. D03
below shows that $\psi=-\frac{\Gamma}{2\pi}\ln r$ does exactly this, although the book never writes it out:
""", equation=EQ["6.6"], ref="6.6")
P("P149", "2-D divergence theorem and the Dirac delta in the plane", r"""
In the plane, Gauss' theorem says the total of ∇²f over a region equals the outward flux of ∇f through its edge:
$\int_A\nabla^2 f\,dA=\oint_C\nabla f\cdot\mathbf n\,ds$. A 2-D delta $\delta(x)\delta(y)$ is zero everywhere except at
the origin, with total 1 over any region that contains it. So a function whose Laplacian is zero away from the origin
but whose flux through every circle round it is the same number F has $\nabla^2 f=F\,\delta(x)\delta(y)$. Demo: the flux
of ∇ ln r (which is $\mathbf e_r/r$) through a circle of radius R.
""", code=r"""
th = np.linspace(0, 2*np.pi, 400, endpoint=False)       # 400 equally spaced angles on a circle [rad]
R = 3.0                                                  # the circle's radius [m]
integrand = np.full(th.size, 1/R) * R                    # ∇ln r · n = 1/R times the arc element ds = R dθ
flux = np.sum(integrand) * (2*np.pi/400)                 # periodic trapezoid sum over the full turn
print(flux, 2*np.pi)                                     # 6.283185… both: the flux of ∇ln r is 2π for any R
""")
D("D03", ref="6.6")
note("N07", "Uniform flow in any direction", r"""
— U = 2, V = 1 m/s: ψ(1, 1) = −1 + 2 = 1 m²/s. Its twin (note `N13`) is $\phi=Ux+Vy$ *(Eq. 6.14)*.
""", equation=EQ["6.7"], ref="6.7")
note("N08", "The ideal vortex at (x′, y′)", r"""
It spins the fluid round at $u_\theta=\frac{\Gamma}{2\pi r}$ *(Eq. 5.2)* (Ch. 3 §3.5, Ch. 5 §5.1).
> ⚠️ Γ is **counterclockwise positive** here and in all our code; the cylinder with circulation (C07) will use the
> book's clockwise Γ.
""", equation=EQ["6.8"], ref="6.8")
note("N14", "The point source (m > 0) or sink (m < 0)", r"""
m is the volume flow out per unit depth [m²/s]: with m = 2π m²/s, $u_r=m/2\pi r=1/r$ — 1 m/s at 1 m, 0.5 m/s at 2 m.
""", equation=EQ["6.15"], ref="6.15")
remind([
    ("Poisson equation (6.11)", r"$\nabla^2 f=$ a given source term; Ch. 5 P139 solved it in 3-D with a Green's function."),
])
note("N10", "Continuity with φ is Poisson's equation", r"""
— with a source density q [1/s] (the φ-twin of $\omega_z=-\nabla^2\psi$):
""", equation=r"\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=\nabla^2\phi=q(x,y)", ref="6.11")
note("N11", "Laplace for φ", r"""
without sources (note `N12`: with a point source, $\nabla^2\phi=m\,\delta(x-x')\,\delta(y-y')$ *(Eq. 6.13)*, which D03
step 8 derived):
""", equation=EQ["6.12"], ref="6.12")
nb.worked_example("a whirlpool and a spring with round numbers", r"""
1. Vortex, Γ = 2π m²/s at the origin: $\psi=-\ln r$ [m²/s]; $u_\theta=-\partial\psi/\partial r=1/r$ → 1 m/s at r = 1 m.
2. Flux of ∇ψ out of the circle r = 2 m: $\partial\psi/\partial r=-1/2$ m/s, times the length 2π × 2 = 4π m → −2π m²/s
   $=-\Gamma$ ✓ (the delta weight of (6.6)), for **any** radius (try r = 5: −1/5 × 10π = −2π).
3. Source, m = 2π m²/s: $u_r=1/r$; flow out of the circle r = 2: 0.5 × 4π = 2π m²/s = m ✓ (6.13).
4. Away from the centres, $\nabla^2\ln r=\frac1r\frac{d}{dr}(r\cdot\frac1r)=0$ (6.23a).
""")
remind([
    ("dataclasses and named results", "fluidpy's elements are small frozen objects with methods: `pf.Vortex(G).psi(x, y)` (Ch. 4 P111)."),
    ("periodic trapezoid rule on a closed loop", "on a closed loop the plain sum × Δθ is spectrally accurate — the demo in the primer above (Ch. 5 P136)."),
])
nb.code(r"""
els = {"uniform": pf.Uniform(2.0, 1.0), "source": pf.Source(2*np.pi), "vortex": pf.Vortex(2*np.pi)}   # three elements
pts = (np.array([0.7, -1.2, 2.0]), np.array([0.4, 0.9, -1.5]))            # three ordinary points (x, y) [m]
for k, e in els.items():                                                  # the Laplacian of ψ and φ of each element
    print(f"{k:8s} max|∇²ψ| = {np.max(np.abs(pf.laplacian_residual(e.psi, *pts))):.1e}",
          f"max|∇²φ| = {np.max(np.abs(pf.laplacian_residual(e.phi, *pts))):.1e}")   # ≈ 0: harmonic away from centres
print("vortex, flux of ∇ψ on r = 0.01…100 m:", ch06.delta_flux_check("vortex", Gamma=2*np.pi))           # −Γ every time
print("source, flux of ∇φ on r = 0.01…100 m:", ch06.delta_flux_check("source", kind="phi", m=2*np.pi))    # +m every time
print("vortex, circle that misses the centre:", ch06.delta_flux_check("vortex", center=(5.0, 0.0), radii=(1.0,), Gamma=2*np.pi))   # 0
""", explain=r"""
1. Three elements as objects with `psi(x, y)` and `phi(x, y)` methods: a uniform stream (U, V) = (2, 1) m/s, a source
   m = 2π m²/s and a vortex Γ = 2π m²/s, both at the origin.
2. `laplacian_residual` applies a 4th-order 9-point stencil at three ordinary points: every element is harmonic there
   (the ~10⁻⁹ is round-off).
3. `delta_flux_check` sums ∇ψ·n (or ∇φ·n) round circles of radius 0.01 … 100 m by the periodic trapezoid rule: −Γ for
   the vortex and +m for the source on every circle — the delta weights of (6.6) and (6.13).
4. A circle centred at (5, 0) with radius 1 m does not enclose the centre: its flux is 0.
""")
nb.code(r"""
x = np.linspace(-0.3, 0.3, 121)                         # a 121 × 121 grid across a Rankine vortex [m]
X, Y = np.meshgrid(x, x)                                # all grid points (rows = y, columns = x)
om = ch06.vorticity_from_psi("rankine", X, Y, Gamma=1.0, a=0.1)   # ω_z = −∇²ψ of a Rankine vortex, Γ = 1 m²/s, core a = 0.1 m
print(f"centre {om[60, 60]:.3f} 1/s, at r = 0.2 m {om[60, 100]:.1e} 1/s")   # Γ/πa² = 31.831 inside, 0 outside
""", explain=r"""
A **Rankine vortex** (Ch. 3 §3.5) has a core of radius a turning like a solid body and an ideal vortex outside;
`rankine_vortex_psi` gives its ψ, $-\Gamma r^2/(4\pi a^2)$ inside and $-\frac{\Gamma}{2\pi}[\ln(r/a)+\tfrac12]$ outside.
`vorticity_from_psi` takes $-\nabla^2\psi$ by a stencil: the uniform core vorticity $\Gamma/\pi a^2=1/(\pi\cdot0.01)=31.83$
s⁻¹ inside and zero outside — Poisson inside, Laplace outside.
""")
remind([
    ("numpy slicing for stencils", "`psi[1:-1, 2:]` is the whole grid shifted one column right; neighbour sums are shifted slices, no loops (Ch. 2 P77)."),
])
nb.md("**From scratch — our own five-point −∇²ψ and our own loop sum of ∂ψ/∂n:**")
nb.check_agree(r"""
psi = ch06.rankine_vortex_psi(X, Y, Gamma=1.0, a=0.1)    # ψ of the Rankine vortex on the grid [m²/s]
h = x[1] - x[0]                                          # grid spacing [m]
lap = (psi[1:-1, 2:] + psi[1:-1, :-2] + psi[2:, 1:-1] + psi[:-2, 1:-1] - 4*psi[1:-1, 1:-1])/h**2   # five-point ∇²ψ [1/s]
rr = np.hypot(X, Y)[1:-1, 1:-1]                          # distance of each interior node from the centre [m]
ok = np.abs(rr - 0.1) > 2*h                              # skip the core edge, where ψ'' jumps
dev = np.max(np.abs(-lap[ok] - om[1:-1, 1:-1][ok]))    # largest difference from the library's 4th-order vorticity [1/s]
print(f"largest difference {dev:.1e} 1/s (core vorticity 31.83 1/s)")   # the five-point stencil's h² error, near the core
assert dev < 0.01*31.83                                 # −∇²ψ by hand = the library's vorticity, to better than 1 %
th = np.linspace(0, 2*np.pi, 512, endpoint=False)        # 512 points on a circle [rad]
R = 0.7                                                  # its radius [m]
u, v = pf.Vortex(1.0).velocity(R*np.cos(th), R*np.sin(th))   # the vortex's velocity there, Γ = 1 m²/s
dpsidn = -(u*(-np.sin(th)) + v*np.cos(th))               # ∂ψ/∂n = ∂ψ/∂r = −u_θ, from (6.22)
mine = np.sum(dpsidn)*R*2*np.pi/512                      # ∮ ∂ψ/∂n ds with ds = R dθ
print(mine)                                              # −1.0 = −Γ
assert np.isclose(mine, ch06.delta_flux_check("vortex", radii=(0.7,), Gamma=1.0)[0])   # = the library's flux
""")
nb.md("Our own stencil and our own loop sum agree with the library.")
remind([
    ("level sets", r"a curve f = const is a level set; ∇f is perpendicular to it (Ch. 2 P75)."),
], lead="before D04: the total differential dψ is recap R03; the dot product and perpendicular vectors are Ch. 2 §2.2")
D("D04", ref="6.10")
note("N09", "The potential", r"""
(in any irrotational flow; single-valued if the region is simply connected, recap R04). D04 showed the lines of constant
φ cross the streamlines at right angles, and ∇ψ is ∇φ turned by +90°:
""", equation=r"u\equiv\partial\phi/\partial x,\ v\equiv\partial\phi/\partial y;\qquad (dy/dx)_{\phi=\text{const}}=-u/v", ref="6.10")
nb.code(r"""
flow = pf.cylinder(1.0, 1.0, Gamma_cw=2.0)               # the cylinder with circulation of C07 (U = 1 m/s, a = 1 m)
pts = (np.array([1.5, -2, 0.3]), np.array([0.4, 1.1, -1.8]))   # three points in the fluid [m]
print(ch06.orthogonality_check(flow, pts))               # max |∇φ·∇ψ|/|u|² ≈ 0: equipotentials ⟂ streamlines
""")
remind([
    ("chain rule", r"$\frac{d}{dt}f(x(t))=f'(x)\,x'(t)$; in polar coordinates it turns ∂/∂x, ∂/∂y into ∂/∂r, ∂/∂θ (Ch. 1 P49)."),
])
note("N18", "Polar velocities", r"""
(note `N19` is the θ-part; a chain-rule exercise the book gives "for quick reference"):
""", equation=r"u_r=\frac{\partial\phi}{\partial r}=\frac1r\frac{\partial\psi}{\partial\theta}\ \ (6.21),\qquad u_\theta=\frac1r\frac{\partial\phi}{\partial\theta}=-\frac{\partial\psi}{\partial r}\ \ (6.22)")
nb.code(r"""
rs, ts = sp.symbols("r theta", positive=True)           # fresh polar symbols (th was reused for numbers above)
print(pf.polar_velocity_sym(-sp.log(rs), rs, ts, kind="psi"))   # the vortex with Γ = 2π: ψ = −ln r → (u_r, u_θ) = (0, 1/r)
""")
nb.figure(r"""
n = 241 if not FAST else 161                             # grid points per side of each flow net
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))         # three panels
flow_net(axs[0], pf.Uniform(1.0, 0.5), xlim=(-2, 2), ylim=(-2, 2), n=n, n_levels=13, phi=True, title="uniform, U = 1, V = 0.5 m/s")
flow_net(axs[1], pf.Source(2*np.pi, cut_angle=0.0), xlim=(-2, 2), ylim=(-2, 2), n=n, n_levels=17, phi=True,
         cut_ray=True, title="source, m = 2π m²/s")      # the branch cut of θ put along +x and blanked
flow_net(axs[2], pf.Vortex(2*np.pi, cut_angle=0.0), xlim=(-2, 2), ylim=(-2, 2), n=n, n_levels=13, phi=True,
         cut_ray=True, title="vortex, Γ = 2π m²/s")      # ψ circles, φ rays
fig.suptitle("Two families of curves, crossing at right angles", fontweight="bold")
savefig(fig, "ch06", "c02_flow_nets"); plt.show()        # save, then draw
""", see="In each panel two families of curves crossing at right angles: straight lines for the stream; rays (ψ) and "
         "circles (φ) for the source; circles (ψ) and rays (φ) for the vortex.",
    read="Solid lines are streamlines (equal Δψ, so crowded lines mean fast flow); dashed teal lines are equipotentials "
         "(equal Δφ); D04's right angles hold everywhere except at the centres. A source is a vortex with ψ and φ swapped.",
    change="…the vortex strength doubled: the same circles, twice as crowded — the velocity doubles (spacing = speed, "
           "the Ch. 4 `stream_function_spacing` explainer).")
nb.figure(r"""
radii = np.logspace(-2, 2, 41)                           # circle radii from 0.01 to 100 m
fv = ch06.delta_flux_check("vortex", radii=radii, Gamma=2*np.pi)              # ∮∇ψ·n ds of the vortex [m²/s]
fs = ch06.delta_flux_check("source", radii=radii, kind="phi", m=2*np.pi)      # ∮∇φ·n ds of the source [m²/s]
ro = radii[np.abs(radii - 5.0) > 0.5]                   # skip circles that pass right through the vortex
fo = ch06.delta_flux_check("vortex", center=(5.0, 0.0), radii=ro, Gamma=2*np.pi)   # vortex at origin, circles round (5, 0)
fig, ax = plt.subplots(figsize=(7, 3.6))                 # one panel
ax.semilogx(radii, fv, color=COLORS["amber"], lw=2.5, label=r"vortex: $\oint\nabla\psi\cdot\mathbf{n}\,ds$")
ax.semilogx(radii, fs, color=COLORS["teal"], lw=2.5, label=r"source: $\oint\nabla\phi\cdot\mathbf{n}\,ds$")
ax.semilogx(ro, fo, color=COLORS["muted"], lw=2, ls="--", label="vortex, circles centred 5 m away")
ax.axvline(5.0, color=COLORS["muted"], lw=0.8, ls=":")   # where the off-centre circles start to enclose the vortex
ax.set_xlabel("circle radius [m]"); ax.set_ylabel("flux [m$^2$/s]")   # axes with units
ax.set_title("The flux does not care about the radius"); ax.legend(fontsize=8)
savefig(fig, "ch06", "c02_flux_vs_radius"); plt.show()   # save, then draw
""", see="Two flat lines at −2π (amber) and +2π (teal); the grey dashed line is 0 until the radius passes 5 m, then −2π.",
    read="The whole Laplacian of the vortex and the source sits at one point — the delta of (6.6) and (6.13): any circle "
         "round it measures the same strength, a circle that misses it measures nothing.",
    change="…a Rankine core of radius 0.1 m: the amber line would rise from 0 at r = 0 to −Γ at r = 0.1 m and stay flat — "
           "the delta is smeared into a disc.")
nb.figure(r"""
fig, ax = plt.subplots(figsize=(5.2, 4.2))               # one panel
im = ax.pcolormesh(x, x, om, cmap="Purples", shading="auto")   # ω_z = −∇²ψ [1/s]
ax.contour(X, Y, psi, levels=14, colors=COLORS["ink"], linewidths=0.7)   # streamlines (circles)
fig.colorbar(im, ax=ax, label=r"$\omega_z=-\nabla^2\psi$ [1/s]")
ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")   # axes with units
ax.grid(False)                                           # no grid lines over the core
ax.set_title("Rankine vortex: spin only in the core", fontsize=10)
savefig(fig, "ch06", "c02_rankine_vorticity"); plt.show()   # save, then draw
""", see="A uniformly coloured disc of radius 0.1 m (31.8 s⁻¹) on a white background, with circular streamlines inside "
         "and outside.",
    read="Outside the core ψ is harmonic and the flow is irrotational although it goes round; inside, ∇²ψ = −ω_z is a "
         "Poisson equation with a uniform right-hand side.",
    change="…a shrunk with Γ fixed: the disc gets smaller and darker (Γ/πa²) — in the limit the point vortex of (6.6).")
note("N15", "Two descriptions, one flow", r"""
ψ and φ each describe the same plane ideal flow and will fuse into one complex potential $w=\phi+i\psi$ in C09. ψ also
works for rotational plane flow (the Rankine core above); φ works for unsteady and 3-D flow (C13, C15).
""")
whatif(r"""
…the vorticity were spread out (a smooth vortex)? Then $\nabla^2\psi=-\omega_z$ is a Poisson equation with a smooth
right-hand side, solved by the Green's-function sum of Ch. 5 §5.5 or on a grid (`ls.solve_poisson`, C12's method with a
source term). That is exactly how ocean and atmosphere models recover the stream function from the vorticity (Ch. 13's
potential-vorticity inversion).
""")

nb.recap("R05", "Bernoulli everywhere in irrotational flow", r"""
Steady, constant density, irrotational: the Bernoulli constant is the same on every streamline,
$p+\tfrac12\rho\lvert\nabla\phi\rvert^2=p+\tfrac12\rho\lvert\nabla\psi\rvert^2=\text{const}$ *(Eq. 6.18, from 4.72)*;
$\lvert\nabla\psi\rvert=\lvert\nabla\phi\rvert$ by D04. Unsteady flow adds ρ∂φ/∂t:
$\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\int\frac{dp}{\rho}+gz=\text{const}$ *(Eq. 4.75)*. So once ψ or
φ is known, the pressure is free — the nonlinear momentum equation is solved by one algebraic line. Symbols: ρ is the
density [kg/m³], p∞ the pressure far upstream [Pa]. Below, `bernoulli_function` (Ch. 4) returns ½|u|² + p/ρ + gz at
four points on different streamlines of the half-body of C05.
""", where="Ch. 4 §4.9")
nb.code(r"""
from fluidpy.core import bernoulli as be                 # Ch. 4's Bernoulli tools
f = pf.half_body(1.0, 2*np.pi)                           # a stream U = 1 m/s plus a source m = 2π m²/s (C05)
X = np.array([-3.0, 0.5, 2.0, -1.0]); Y = np.array([2.0, 3.5, -4.0, 1.5])   # four points on different streamlines [m]
p = f.pressure(X, Y, rho=1000.0, p_inf=0.0, U=1.0)       # p − p∞ from (6.18) [Pa], water
print(be.bernoulli_function(f.speed(X, Y), p, 0.0, rho=1000.0, g=0.0))   # ½|u|² + p/ρ: the same 0.5 J/kg at all four
""")

# ---- C03 ------------------------------------------------------------------------------------------------------------
core("C03", "Superposition — and any streamline can be a wall", r"""
If every element solves Laplace's equation, where does the body come from — and who imposes its boundary condition?

*In one line:* $\nabla^2(\phi_1+\phi_2)=\nabla^2\phi_1+\nabla^2\phi_2=0$, and on a stationary wall
$\partial\phi/\partial n=0$ or $\partial\psi/\partial s=0$ *(Eq. 6.16)*.
""")
problem(r"""
Dye streaks in a river bend round a stone you cannot see. Could we produce that pattern without ever drawing the stone —
by adding flows we already know and then *declaring* one of the streamlines to be a wall?
""")
idea("""
uniform stream  +  source  =  a flow with a stagnation point S
                                   |
     the streamline through S -----+---> call it a wall: no flow crosses a streamline,
                                          so the no-through-flow condition (6.16) already holds there
""", r"""
**Laplace's equation is linear, so solutions add; the boundary condition is met by choosing which streamline is the
body.** A *stagnation point* is a point where the fluid is at rest, u = v = 0 (Ch. 3 §3.3).
""")
remind([
    ("directional derivative ∂φ/∂n", r"$\partial\phi/\partial n=\mathbf n\cdot\nabla\phi$ is the rate of change of φ along the unit vector n (Ch. 2 P75)."),
    ("unit normal and tangent of a curve", "a curve's unit tangent t points along it; turning t by 90° gives the unit normal n (Ch. 3 P92)."),
    ("complex numbers as points x + iy in code (`-1+0j`)", r"""fluidpy stores a point (x, y) as the complex number
     x + iy, written `-1+0j` for (−1, 0), and returns stagnation points the same way (i² = −1 was met in Ch. 1 P45;
     C09 explains complex numbers properly, primer P153) — here read x + iy simply as the pair of coordinates."""),
])
D("D05", ref="6.16")
note("N16", "No flow through a solid surface", r"""
in its general and its stationary form. A stationary wall is a streamline, and any streamline may be replaced by a wall:
""", equation=r"\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{surface}};\qquad \partial\phi/\partial n=0\ \text{or}\ \partial\psi/\partial s=0\ \text{on a stationary surface}", ref="6.16")
remind([
    ("observed order from a log–log slope", "if an error behaves like C·hᵖ, the slope of log(error) against log(h) is p (Ch. 1 P13)."),
])
note("N17", "The far field", r"""
Far away the flow must return to the given stream (U along x; U = 0 is still fluid). How fast? For a closed body the
disturbance dies like 1/R² (a doublet, C04); for a half-body with a net source like 1/R:
""", equation=r"\partial\phi/\partial x=U,\quad\text{or}\quad\partial\psi/\partial y=U", ref="6.17")
nb.code(r"""
Rs = np.array([10.0, 20.0, 40.0, 80.0])                  # radii of four far circles [m]
e_cyl = [pf.far_field_check(pf.cylinder(1.0, 1.0), R) for R in Rs]        # max |(u − iv) − U| on |z| = R: closed body
e_half = [pf.far_field_check(pf.half_body(1.0, 2*np.pi), R) for R in Rs]  # the same for the open half-body
print(f"closed cylinder: slope {observed_order(Rs, e_cyl):.2f}   half-body: slope {observed_order(Rs, e_half):.2f}")   # −2 and −1
""")
nb.worked_example("a stagnation point from two elements", r"""
Stream U = 1 m/s plus a source m = 2π m²/s at (−1, 0) m. At the point (−2, 0): the source is 1 m away, straight downstream
of the point, so it pushes the fluid there **upstream** at
$u=\frac{m}{2\pi}\frac{x-x'}{(x-x')^2+(y-y')^2}=\frac{2\pi}{2\pi}\cdot\frac{-1}{1}=-1$ m/s.

1. Sum: $u=1+(-1)=0$.
2. On the axis y = 0 the source adds no v, so $v=0$.
3. The fluid stops: a stagnation point, 1 m in front of the source — the nose of a body that nobody drew (C05 finds its
   whole shape).
""")
nb.code(r"""
f = pf.Flow([pf.Uniform(1.0), pf.Source(2*np.pi, -1+0j, cut_angle=2*np.pi)])   # stream U = 1 m/s + source m = 2π m²/s at (−1, 0)
print("velocity at (−2, 0):", f.velocity(-2.0, 0.0))     # (0, 0): a stagnation point
print("stagnation points:", f.stagnation_points())       # found by Newton's method on dw/dz: [-2+0j]
print("ψ there:", f.psi(-2.0, 0.0))                      # m/2 = π m²/s names the dividing streamline
th = np.linspace(0.05, 2*np.pi - 0.05, 4000)             # polar angles along the body (C05's shape) [rad]
xb, yb = ch06.half_body_shape(1.0, 2*np.pi, th)          # the dividing streamline of a source at the origin [m]
body = (xb - 1.0) + 1j*yb                                # shifted to the source at x = −1, as complex points
print("spread of ψ along it:", np.ptp(f.psi(body.real, body.imag)))   # np.ptp = max − min ≈ 0: one streamline
print("max normal velocity:", pf.normal_velocity_on(f, body, closed=False))   # ≈ 0: no flow crosses it
""", explain=r"""
1. A flow is a list of elements; its velocity, ψ and φ are the sums of theirs. The source's ψ contains its polar angle
   θ′, which must jump by 2π somewhere; `cut_angle=2*np.pi` puts that jump along +x behind the source, *inside* the
   body, and counts θ′ from 0 to 2π, so ψ is smooth everywhere in the fluid (C05 says more about choosing the angle's range).
2. `stagnation_points` solves u − iv = 0 by Newton's method from several starting points (the method is primer P152 in
   C07) and checks each root.
3. The streamline through the stagnation point has ψ = m/2 (D07 will show why).
4. Along that curve ψ is constant to round-off, and the velocity component normal to the curve (normals from the 4000
   points of the polyline) is tiny — the small remainder is the error of approximating the curve by straight pieces. So
   no flow crosses the curve: by D05 it can be a wall.
""")
nb.md("**From scratch — the two velocities added by hand at one point:**")
nb.check_agree(r"""
x, y = 0.5, 1.5                                          # a point in the flow [m]
u_hand = 1.0 + (2*np.pi/(2*np.pi))*(x + 1)/((x + 1)**2 + y**2)   # stream + source x-velocity (source at x = −1) [m/s]
v_hand = (2*np.pi/(2*np.pi))*y/((x + 1)**2 + y**2)       # the source's y-velocity [m/s]
print(u_hand, v_hand)                                    # our numbers
assert np.allclose((u_hand, v_hand), f.velocity(x, y))   # the library does nothing more than add
""")
nb.figure(r"""
n = 241 if not FAST else 161                             # grid points per side
fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))         # stream | source | sum
for yl in np.linspace(-3, 3, 13):                        # panel 1: streamlines of the stream are horizontal lines
    axs[0].axhline(yl, color=COLORS["blue"], lw=0.8)     # ψ = Uy = const
for ang in np.linspace(0, 2*np.pi, 16, endpoint=False):  # panel 2: streamlines of the source are rays
    axs[1].plot([-1, -1 + 8*np.cos(ang)], [0, 8*np.sin(ang)], color=COLORS["teal"], lw=0.8)   # ψ = mθ′/2π = const
axs[1].plot([-1], [0], "o", color=COLORS["teal"], ms=6)  # the source itself
for a_, t_ in zip(axs[:2], ("uniform stream, U = 1 m/s", "source m = 2π m²/s at (−1, 0)")):
    a_.set_xlim(-4, 4); a_.set_ylim(-3, 3); a_.set_aspect("equal"); a_.set_title(t_, fontsize=10)   # same frame
    a_.set_xlabel("x [m]"); a_.set_ylabel("y [m]")       # axes with units
xs_, ys_ = np.meshgrid(np.linspace(-4, 4, n), np.linspace(-3, 3, n))   # grid for the sum [m]
P3 = f.psi(xs_, ys_)                                     # ψ of the sum [m²/s]
P3 = np.where((np.abs(ys_) < 0.05) & (xs_ > -1), np.nan, P3)   # blank the source's cut, which lies inside the body
tb = np.linspace(0.02, 2*np.pi - 0.02, 800)              # angles along the body [rad]
xb2, yb2 = ch06.half_body_shape(1.0, 2*np.pi, tb)        # the body ψ = m/2 of a source at the origin [m]
axs[2].fill(xb2 - 1.0, yb2, color="#e8ebf2", zorder=0)   # shade the source's fluid (the body is shifted to x = −1)
axs[2].contour(xs_, ys_, P3, levels=np.arange(-3.0, 9.5, 0.5), colors=COLORS["ink"], linewidths=0.7)   # streamlines, Δψ = 0.5 m²/s
axs[2].plot(xb2 - 1.0, yb2, color="k", lw=2.5)           # the dividing streamline ψ = m/2, bold
axs[2].plot([-2], [0], "o", color="k", ms=6)             # the stagnation point S
axs[2].annotate("S", (-2, 0), (-2.6, 0.5)); axs[2].set_aspect("equal")
axs[2].set_xlim(-4, 4); axs[2].set_ylim(-3, 3); axs[2].set_xlabel("x [m]"); axs[2].set_ylabel("y [m]")
axs[2].set_title("sum: a rounded nose appears", fontsize=10)
fig.suptitle("Add two flows, then choose a streamline as the wall", fontweight="bold")
savefig(fig, "ch06", "c03_stream_plus_source"); plt.show()   # save, then draw
""", see="Straight blue lines, teal rays, and their sum: the streamlines part round a rounded nose; the bold black curve "
         "through S encloses the grey region filled by the source's fluid.",
    read="No element knows about a body; the bold line is simply the streamline through S (ψ = m/2) — by D05 it is a wall, "
         "and everything outside it is the flow past a blunt body.",
    change="…a sink of equal strength added 3 m downstream: the bold curve closes into an oval (the Rankine oval — try it "
           "in the superposition explainer after C05).")
whatif(r"""
…we added a vortex instead of a source? A vortex alone in a uniform stream makes a stagnation point off to one side but
no closed dividing streamline; combined with a closed body it produces lift (C07). The rule stays the same: add, then
find the streamline that can be a wall.
""")

# =====================================================================================================================
# A.3 §6.3 (first half) — R10 R11, C04 (N20–N26, D06), C05 (N27–N29, D07, D08, E1)
# =====================================================================================================================
nb.section("6.3", "Construction of Elementary Flows in Two Dimensions", intro=r"""
**What is this section about?** A kit of elementary flows — stream, source, vortex, corner, doublet — and what you get
by adding them: a blunt half-body, a cylinder with zero drag (d'Alembert's paradox), a spinning cylinder with lift, and
walls made by mirror images. Every body here is a streamline of a sum.
""")
nb.recap("R10", "The stagnation (corner) flow", r"""
The simplest quadratic stream function, $\psi=2Axy$ *(Eq. 6.24)* (A in 1/s), gives u = 2Ax, v = −2Ay: fluid comes down
the y-axis, stops at the origin and leaves along the x-axis; streamlines are the hyperbolae xy = const. It is Ch. 3's
irrotational strain and the Ch. 4 "stagnation" preset. ⚠️ In Ch. 13 the same field is the deformation flow that sharpens
weather fronts. In code it is `pf.Corner(A, n=2)` (the n is explained in C09).
""", where="Ch. 3 §3.4")
nb.code(r"""
c = pf.Corner(1.0, n=2)                                  # the corner flow with A = 1 1/s
print(c.velocity(1.0, 0.5), c.psi(1.0, 0.5))             # (u, v) = (2Ax, −2Ay) = (2, −1) m/s and ψ = 2Axy = 1 m²/s
""")
nb.recap("R11", "The vortex velocity", r"""
Differentiating $\psi=-\frac{\Gamma}{2\pi}\ln\sqrt{(x-x')^2+(y-y')^2}$ *(Eq. 6.8)* gives $u=-\frac{\Gamma}{2\pi}\frac{y}{x^2+y^2}$,
$v=\frac{\Gamma}{2\pi}\frac{x}{x^2+y^2}$, i.e. $u_\theta=\frac{\Gamma}{2\pi r}$ *(Eq. 5.2)*: circles round the centre,
irrotational everywhere else. Ch. 5's Biot–Savart tool gives the same velocity.
""", where="Ch. 5 §5.1")
nb.code(r"""
from fluidpy.core import biot_savart as bs               # Ch. 5's point-vortex tools
print(pf.Vortex(1.0).velocity(0.3, 0.4))                 # our vortex, Γ = 1 m²/s at the origin, seen from (0.3, 0.4) m
print(bs.point_vortex_velocity(np.array([0.3, 0.4]), np.array([[0.0], [0.0]]), np.array([1.0])))   # Ch. 5: the same
""")

# ---- C04 ------------------------------------------------------------------------------------------------------------
core("C04", "The element kit, and the doublet as a merged source–sink pair", r"""
What is left when a source and a sink merge — and why is that the piece that makes closed bodies?
""", eqs=("6.29",))
problem(r"""
Hold two garden hoses nose to nose, one blowing water out, one sucking the same amount in. From far away the water just
loops from one to the other. Push them closer while turning both up so the loops do not vanish: the limit is a new
element, the doublet — the only one that disturbs the flow without adding or removing water or circulation. It is what
turns a stream into flow round a cylinder, a sphere, an airship.
""")
idea("""
source +m at (-eps, 0)   sink -m at (+eps, 0)        eps -> 0 with 2*m*eps = |d| fixed
      *---------------------o             ---->      (.)  doublet: streamlines are circles tangent to the x-axis
dipole vector d = sum of x_i m_i = -2*m*eps e_x  (it points from the sink to the source)
""")
remind([
    ("first-order Taylor expansion, ln(1 + s) ≈ s", r"for small s, $\ln(1+s)=s-s^2/2+\dots\approx s$ (Ch. 1 P26)."),
    ("limits and orders of smallness O(ε²)", "O(ε²) means \"shrinks at least as fast as ε²\" — halving ε quarters it (Ch. 2 P68)."),
], lead=r"log rules used below: $\ln\sqrt s=\tfrac12\ln s$, $\ln(ab)=\ln a+\ln b$, $\frac{\partial}{\partial x}\ln r=x/r^2$")
note("N20", "Harmonic polynomials", r"""
A constant (no flow), the linear stream $\psi=-Vx+Uy$ *(Eq. 6.7)* / $\phi=Ux+Vy$ *(Eq. 6.14)*, and two quadratic families
solve Laplace; in general the real and imaginary parts of $z^n=(x+iy)^n$ do (C09). Higher powers make sharper corners,
fractional powers wider ones.
""")
nb.code(r"""
for d in (2, 3):                                         # degree of the polynomial
    print(d, ch06.harmonic_polynomials(d)[-1])           # (Re zⁿ, Im zⁿ): each is harmonic (sympy checks ∇² = 0)
""", explain=r"""
sympy splits $z^n=(x+iy)^n$ into its real and imaginary parts; each is asserted harmonic inside the function (a preview
of D14 in C09): $x^2-y^2$ and $2xy$ for n = 2, $x^3-3xy^2$ and $3x^2y-y^3$ for n = 3.
""")
note("N21", "The corner family (notes `N22`, `N23` too)", r"""
Swapping which function is which, or turning the axes by 45°, gives the relatives of $\psi=2Axy$ *(Eq. 6.24)*:
$\phi=2Axy$ *(6.25)*, $\psi=A(x^2-y^2)$ *(6.26)*, $\phi=A(x^2-y^2)$ *(6.27)*. (6.25) is the (6.24) flow turned by 45°
(x² − y² in axes turned by 45° becomes 2x′y′); (6.27) *is* (6.24): both are parts of the one function Az².
""", equation=r"\phi=2Axy\ \ (6.25),\qquad \psi=A(x^2-y^2)\ \ (6.26),\qquad \phi=A(x^2-y^2)\ \ (6.27)")
nb.code(r"""
x, y = 0.7, 0.3                                          # a point [m]
print(pf.Corner(1.0, n=2).phi(x, y), x**2 - y**2)        # (6.27) is the φ of the (6.24) flow: 0.4 both
c45 = pf.Corner(1.0, n=2, rotate=np.pi/4)                # the same corner turned by 45°: w = A(z e^{−iπ/4})² = −iAz²
print(c45.phi(x, y), 2*x*y, c45.psi(x, y), -(x**2 - y**2))   # φ = 2Axy is (6.25); ψ = −A(x² − y²) is (6.26) with A → −A
""")
note("N24", "Source velocities", r"""
Differentiating the source potential $\phi=\frac{m}{2\pi}\ln\sqrt{x^2+y^2}$ *(Eq. 6.15)* gives $u=\frac{m}{2\pi}\frac{x}{x^2+y^2}$,
$v=\frac{m}{2\pi}\frac{y}{x^2+y^2}$, so $u_r=m/2\pi r$, $u_\theta=0$. ⚠️ The book says "differentiation of (6.8)" here; it
means (6.15). ∇·u = 0 except at r = 0, where all m m²/s come out (D03).
""")
note("N25", "The pair", r"""
A source +m at (−ε, 0) and a sink −m at (+ε, 0). The dipole vector points from the sink to the source:
""", equation=r"\phi=\frac{m}{2\pi}\ln\sqrt{(x+\varepsilon)^2+y^2}-\frac{m}{2\pi}\ln\sqrt{(x-\varepsilon)^2+y^2},\qquad\mathbf d=\sum_i\mathbf x_im_i=-2m\varepsilon\,\mathbf e_x", ref="6.28")
P("P150", "a limit with a product held fixed", r"""
Sometimes two quantities go to extremes together: here the gap ε → 0 while the strength m → ∞, with the product 2mε
held at a fixed value |d|. Either limit alone gives nothing useful (ε → 0 at fixed m cancels the pair; m → ∞ at fixed ε
blows up). Held together, the leading term survives and the rest shrinks like ε².
""", code=r"""
d = 2.0                                                  # fixed dipole strength |d| = 2mε [m³/s]
for eps in (0.2, 0.02, 0.002):                           # shrink the gap …
    m = d/(2*eps)                                        # … and raise the strength to keep 2mε = d [m²/s]
    phi = m/(2*np.pi)*np.log((1 + eps)/(1 - eps))        # the pair's potential at (1, 0) m [m²/s]
    print(eps, phi, d/(2*np.pi))                         # → 0.31831 (the doublet value); the gap shrinks like ε²
""")
D("D06", ref="6.29")
nb.worked_example("a strong pair, 4 cm apart", r"""
m = 50 m²/s, ε = 0.02 m, so $\lvert\mathbf d\rvert=2m\varepsilon=2$ m³/s. At (1, 0) m:

1. pair, (6.28): $\phi=\frac{50}{2\pi}[\ln1.02-\ln0.98]=7.9577\times0.0400053=0.318352$ m²/s.
2. doublet, (6.29): $\phi=\frac{\lvert\mathbf d\rvert}{2\pi}\frac{\cos0}{1}=0.318310$ m²/s.
3. relative difference 1.33×10⁻⁴ ≈ ε²/3 — the dropped terms of D06 step 9.
""")
nb.code(r"""
eps = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])          # five gaps, each half the last [m]
err = ch06.doublet_limit_error(eps, d=2.0)               # max relative |φ_pair − φ_doublet| at three points, 2mε = 2 m³/s
print(err, f"observed order {observed_order(eps, err):.2f}")   # errors fall ×4 per halving: order 2
D = pf.Doublet((-2.0, 0.0))                              # the doublet with dipole vector d = −2 e_x [m³/s], as the pair
print(D.phi(1.0, 0.0), D.psi(0.0, 1.0))                  # φ = |d|cosθ/2πr = 0.31831; ψ = −|d| y/2πr² = −0.31831 [m²/s]
print(pf.Doublet.from_book_scalar(2.0).phi(1.0, 0.0))    # (6.49)'s scalar d = 2 means the same dipole −2 e_x
""", explain=r"""
1. `doublet_limit_error` compares the pair's φ (6.28) with the doublet's (6.29) at three points for five gaps, holding
   2mε = 2 m³/s.
2. The log–log slope is 2: the error ∝ ε², as D06 step 9 predicted.
3. `pf.Doublet` stores the dipole *vector* d; the one made by our pair points along −x.
4. The book's complex-potential form (6.49) uses a scalar d that means the dipole −d e_x (conventions cell);
   `from_book_scalar` converts it.
""")
nb.md("**From scratch — the pair formula at shrinking ε, by hand, against `pf.Doublet`:**")
nb.check_agree(r"""
for e in eps:                                            # the five gaps [m]
    m = 2.0/(2*e)                                        # strength that keeps 2mε = 2 m³/s [m²/s]
    pair = m/(2*np.pi)*(np.log(np.hypot(0.6 + e, 0.8)) - np.log(np.hypot(0.6 - e, 0.8)))   # (6.28) at (0.6, 0.8) m
    print(f"ε = {e:.4f}: pair {pair:.7f}, doublet {D.phi(0.6, 0.8):.7f}")   # the pair closes in on the doublet
    assert np.isclose(pair, D.phi(0.6, 0.8), rtol=e**2)  # relative error below ε²
""")
nb.md("The pair tends to the doublet with an error that falls like ε², as D06 predicted.")
note("N26", "Where elements are singular", r"""
(6.7), (6.14) and (6.24)–(6.27) are harmonic everywhere; (6.8), (6.15) and (6.29) are harmonic for r > 0 only — their
centres are where the flow's "information" sits (C02). `laplacian_residual` returns NaN within five stencil steps of a
centre on purpose.
""")
nb.figure(r"""
n = 201 if not FAST else 141                             # grid points per side
fig, axs = plt.subplots(1, 5, figsize=(13, 2.9))         # five elements side by side
kit = [("uniform", pf.Uniform(1.0), False), ("source", pf.Source(2*np.pi, cut_angle=0.0), True),
       ("vortex", pf.Vortex(2*np.pi), False), ("corner n = 2", pf.Corner(1.0, n=2), False),
       ("doublet d = −2 e$_x$", pf.Doublet((-2.0, 0.0)), False)]   # (title, element, blank its cut?)
cols = [COLORS["blue"], COLORS["teal"], COLORS["amber"], COLORS["ink"], COLORS["accent"]]   # one colour per element
for ax, (name, el, cut), col in zip(axs, kit, cols):
    X_, Y_ = np.meshgrid(np.linspace(-2, 2, n), np.linspace(-2, 2, n))   # grid [m]
    Pk = np.asarray(el.psi(X_, Y_), float)               # ψ of the element [m²/s]
    if cut:
        Pk = np.where((np.abs(Y_) < 0.03) & (X_ > 0), np.nan, Pk)   # blank the source's branch cut (+x axis)
    lv = np.nanpercentile(Pk, np.linspace(3, 97, 17))    # equally spaced levels between the 3rd and 97th percentiles
    ax.contour(X_, Y_, Pk, levels=np.unique(lv), colors=col, linewidths=0.9)   # streamlines
    ax.set_aspect("equal"); ax.set_title(name, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
axs[3].axhline(0, color="k", lw=2); axs[3].axvline(0, color="k", lw=2)   # the corner's walls (the axes)
fig.suptitle("Five elements, all harmonic (4 m × 4 m windows)", fontweight="bold")
savefig(fig, "ch06", "c04_element_kit"); plt.show()      # save, then draw
""", see="Five streamline patterns: parallel lines, rays, circles, hyperbolae in a right-angle corner, and circles that "
         "all touch the x-axis at the origin.",
    read="Every panel solves Laplace; only the singular centres differ. The doublet's circles are the loops between a "
         "source (left) and a sink (right) squeezed to a point.",
    change="…the doublet turned to point along +y: its circles turn by 90° and touch the y-axis.")
remind([
    ("animate and show_animation", "`animate(update, frames, fig)` calls `update(i)` once per frame; `show_animation(anim, player=\"frames\")` embeds it with step buttons (Ch. 1 P16)."),
])
nb.animation(r"""
eps_list = tuple(np.geomspace(0.5, 0.02, 12 if not FAST else 8))   # gaps from 0.5 m down to 2 cm
frames = ch06.doublet_limit_frames(eps_list, d=2.0, n=161 if not FAST else 121)   # ψ of the pair on a grid, per gap
fig, ax = plt.subplots(figsize=(4.8, 3.7))               # one panel
def update(i):                                           # draw frame i
    fr = frames[i]                                       # the i-th gap
    ax.clear()                                           # contours are redrawn each frame
    ax.contour(fr["X"], fr["Y"], fr["psi"], levels=np.linspace(-1.5, 1.5, 25), colors=COLORS["accent"], linewidths=0.8)
    ax.plot([-fr["eps"]], [0], "o", color=COLORS["teal"], ms=7)   # the source (+m) at (−ε, 0)
    ax.plot([fr["eps"]], [0], "o", color=COLORS["rose"], ms=7)    # the sink (−m) at (+ε, 0)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.2, 1.2); ax.set_aspect("equal")   # fixed frame
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"ε = {fr['eps']:.3f} m, 2mε = 2 m³/s: relative error {fr['err']:.1e}", fontsize=9)
show_animation(animate(update, frames=len(frames), fig=fig), player="frames")   # step through with the buttons
""")
see_read_change(
    see="A teal source and a rose sink moving together; the loops between them tighten into circles tangent to the x-axis.",
    read="The title's error (pair against doublet) falls by about 4 each time ε halves — the ε² of D06 step 9.",
    change="…m held fixed instead of 2mε: the loops would fade away as the pair cancels itself.")
whatif(r"""
…the pair were a vortex pair (+Γ at (0, ε), −Γ at (0, −ε)) squeezed with 2Γε fixed? You get a doublet again, turned by
90° — the same element (an exercise-style check, see the S02 pointer in §6.10).
""")

# ---- C05 ------------------------------------------------------------------------------------------------------------
core("C05", "The half-body: a stream plus a source", r"""
A stream meets a source. Where does the fluid stop, what shape is the body, and how wide does it get?
""", eqs=("6.31",))
problem(r"""
River water parts round the rounded nose of a bridge pier; wind parts round the leading edge of a wing or a cliff top.
The simplest model of such a nose is a uniform stream plus one source: the source's fluid fills a body, the stream's
fluid flows round it. We want the stagnation point, the shape, the width — and the pressure along the surface.
""")
idea("""
stream U  -->              psi = m/2 (dividing streamline, bold)
             S *-- a --* source m     ______________________  h --> h_max = m/2U
stagnation where the source's push m/2*pi*a equals U     inside: all the source's fluid, carried off at speed U
""")
remind([
    ("np.arctan2", "`np.arctan2(y, x)` returns the polar angle of (x, y) in (−π, π], quadrant-correct (Ch. 2 P70)."),
    ("scipy.optimize.brentq", "`brentq(f, a, b)` finds a root of f between a and b where f changes sign (Ch. 3 P108)."),
], lead=r"""**choosing the angle's range** — ψ contains θ itself, which jumps by 2π across a cut. We take θ ∈ [0, 2π)
(`np.mod(np.arctan2(y, x), 2*np.pi)`), so the jump lies along +x behind the source, *inside* the body, and ψ is smooth
everywhere in the fluid; polar coordinates are Ch. 3 P105""")
note("N27", "The same flow in φ", "", equation=EQ["6.30"], ref="6.30")
D("D07", ref="6.31")
nb.worked_example("a pier nose with round numbers", r"""
U = 1 m/s, m = 2π m²/s.

1. $a=m/2\pi U=1$ m: the stagnation point is at (−1, 0) m.
2. The dividing streamline has ψ = m/2 = π m²/s.
3. Right above the source (θ = 90°): $h=m(\pi-\theta)/2\pi U=2\pi\cdot(\pi/2)/(2\pi)=\pi/2\approx1.571$ m.
4. Far downstream $h_{\max}=m/2U=\pi\approx3.142$ m.
5. Mass check: 2 × 3.142 m × 1 m/s = 6.283 m²/s = m ✓.
""")
nb.code(r"""
f = pf.half_body(1.0, 2*np.pi)                           # stream U = 1 m/s + source m = 2π m²/s at the origin (θ ∈ [0, 2π))
print(ch06.half_body_numbers(1.0, 2*np.pi))              # a = m/2πU, the nose x, ψ on the body = m/2, h_max = m/2U
print(f.stagnation_points())                             # Newton finds the nose: [-1+0j]
th = np.array([np.pi/2, np.pi/4, 0.05, 3*np.pi/2])       # four polar angles, two above and one below the axis [rad]
xb, yb = ch06.half_body_shape(1.0, 2*np.pi, th)          # body points from D07 step 6 [m]
print(np.round(yb, 4), f.psi(xb, yb))                    # y = 1.5708, 2.3562, 3.0916, −1.5708; ψ = π at all four
""", explain=r"""
1. `half_body` is the stream plus the source, with θ ∈ [0, 2π) so that ψ = m/2 on both halves of the body.
2. Its three numbers: a = 1 m, ψ_body = π m²/s, h_max = π m.
3. Newton's method finds the stagnation point (−1, 0).
4. Points on the body from D07 step 6 lie on ψ = m/2, above **and below** the axis (with numpy's principal angle range
   the lower half would show −m/2 — the trap of D07).
""")
note("N28", "Pressure coefficient", r"""
Bernoulli $p+\tfrac12\rho\lvert\nabla\phi\rvert^2=\text{const}$ *(Eq. 6.18)* with the constant set far upstream, p∞ + ½ρU².
It is the Euler number of $C_p\equiv\frac{p-p_\infty}{\tfrac12\rho U^2}$ *(Eq. 4.106)* (Ch. 4 §4.11): C_p = 1 at every
stagnation point, 0 where the speed equals U, negative (suction) where the flow is faster than the stream.
""", equation=EQ["6.32"], ref="6.32")
D("D08", ref="6.32")
note("N29", "Along the body", r"""
C_p starts at +1 at the nose, crosses zero where tan θ = −2(π − θ) — our root 113.2° from +x (the curve the book plots
agrees; its printed value stays in the private test file) — has its lowest value −0.587 near θ = 63°, and returns to 0
far downstream. The net pressure force on the whole (infinitely long) body is zero: the tail contributes nothing because
C_p → 0 there.
""")
nb.code(r"""
print(f"C_p = 0 at θ = {np.degrees(ch06.half_body_cp_zero_angle()):.3f}°")   # the root of tan θ = −2(π − θ): 113.218°
print(ch06.half_body_surface_cp(np.radians([179.0, 120.0, 90.0, 63.0, 10.0])))   # D08's C_p at five body angles
for xe in (10.0, 100.0, 1000.0):                         # integrate the body pressure from the nose to x = x_end [m]
    print(xe, ch06.half_body_net_force(1.0, 2*np.pi, xe, rho=1.0))   # drag D → 0 (about ∝ 1/x_end²), lift L = 0 [N/m]
""", explain=r"""
1. `half_body_cp_zero_angle` solves D08 step 8 with `brentq` on (π/2, π).
2. `half_body_surface_cp` evaluates $C_p=-(2k\cos\theta+k^2)$, $k=\sin\theta/(\pi-\theta)$: +0.9997 near the nose,
   −0.4053 at 90°, −0.5865 at 63°.
3. `half_body_net_force` integrates −(p − p∞)n over the body from the nose to x_end: the drag shrinks toward zero as the
   body is lengthened (ρ = 1 kg/m³, U = 1 m/s).
""")
nb.md("**From scratch — our own bracket and our own function for the zero of C_p:**")
nb.check_agree(r"""
from scipy.optimize import brentq                        # root finder with a sign-change bracket (P108)
f_root = lambda t: np.tan(t) + 2*(np.pi - t)             # C_p = 0 on the body (D08 step 8), t in radians
t0 = brentq(f_root, np.pi/2 + 1e-9, np.pi - 1e-9)        # tan changes sign between 90° and 180°
print(np.degrees(t0))                                    # 113.218°
assert np.isclose(t0, ch06.half_body_cp_zero_angle())    # = the library's angle
""")
nb.figure(r"""
n = 241 if not FAST else 161                             # grid points per side
X_, Y_ = np.meshgrid(np.linspace(-3, 8, n), np.linspace(-5, 5, n))   # grid [m]
Ph = f.psi(X_, Y_)                                       # ψ of the half-body [m²/s]
Ph = np.where((np.abs(Y_) < 0.06) & (X_ > 0), np.nan, Ph)   # blank the cut along +x (inside the body)
tb = np.linspace(0.02, 2*np.pi - 0.02, 800)              # angles along the body [rad]
xb_, yb_ = ch06.half_body_shape(1.0, 2*np.pi, tb)        # the body ψ = m/2 [m]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.8))  # shape | pressure
a1.fill(xb_, yb_, color="#e8ebf2", zorder=0)             # the source's fluid
a1.contour(X_, Y_, Ph, levels=np.pi + np.arange(-8, 9)*0.6, colors=COLORS["blue"], linewidths=0.8)   # equal Δψ streamlines
a1.plot(xb_, yb_, "k", lw=2.5)                           # the dividing streamline ψ = m/2
a1.plot([-1], [0], "ko", ms=5); a1.annotate("S", (-1, 0), (-1.9, 0.4))   # the stagnation point
a1.annotate("", (0, 0), (-1, 0), arrowprops=dict(arrowstyle="<->")); a1.text(-0.8, -0.5, "a = 1 m", fontsize=8)
a1.axhline(np.pi, color=COLORS["muted"], ls="--", lw=1); a1.text(4.5, np.pi + 0.2, "$h_{max}$ = π m", fontsize=8)
a1.set_xlim(-3, 8); a1.set_ylim(-5, 5); a1.set_aspect("equal"); a1.set_xlabel("x [m]"); a1.set_ylabel("y [m]")
a1.set_title("U = 1 m/s, m = 2π m²/s", fontsize=10)
tt = np.radians(np.linspace(179.9, 5, 400))              # body angles from the nose (180°) to far downstream [rad]
cp = ch06.half_body_surface_cp(tt)                       # C_p on the body (D08)
a2.plot(np.degrees(tt), cp, color=COLORS["ink"], lw=2)   # the pressure curve
t0d = np.degrees(ch06.half_body_cp_zero_angle())         # zero crossing [deg]
a2.plot([t0d], [0], "o", color=COLORS["orange"]); a2.annotate(f"C_p = 0 at {t0d:.1f}°", (t0d, 0), (t0d + 5, 0.35), fontsize=8)
imin = np.argmin(cp); a2.plot(np.degrees(tt[imin]), cp[imin], "o", color=COLORS["blue"])   # the suction minimum
a2.annotate(f"min {cp[imin]:.3f} at {np.degrees(tt[imin]):.0f}°", (np.degrees(tt[imin]), cp[imin]), (70, -0.35), fontsize=8)
a2.axhline(0, color=COLORS["muted"], lw=0.8); a2.invert_xaxis()   # nose on the left, as in the flow picture
a2.set_xlabel("polar angle θ on the body [deg] (nose = 180°)"); a2.set_ylabel("$C_p$ [–]")
fig.suptitle("Blunt nose: pressure rises at the front, suction on the shoulders", fontweight="bold")
savefig(fig, "ch06", "c05_half_body"); plt.show()        # save, then draw
""", see="Left, a rounded nose opening to a width 2h_max with streamlines crowding over its shoulders; right, C_p falling "
         "from +1 at the nose through zero at 113° to a suction dip of −0.587 and back toward zero downstream.",
    read="High pressure where the stream is stopped (the nose), low where it is fastest (the shoulder), recovery downstream; "
         "streamline spacing is speed.",
    change="…m doubled: every length (a, h, h_max) doubles, the C_p curve against θ is unchanged — only m/U sets the size.")
remind([
    ("slider_figure", "`slider_figure(fn, name, values, …)` computes a plotly figure for every slider position up front, so it works on the web page (Ch. 1 P17)."),
])
nb.plotly(r"""
def half_body_traces(mU):                                # all curves for one value of m/U [m], with U = 1 m/s
    m = mU                                               # source strength [m²/s]
    tb = np.linspace(0.01, 2*np.pi - 0.01, 150)          # body angles [rad]
    xb, yb = ch06.half_body_shape(1.0, m, tb)            # the body ψ = m/2
    xs, ys = [], []                                      # four outer streamlines, joined with NaN gaps
    for c in (m/2 + 1.0, m/2 + 2.5, m/2 - 1.0, m/2 - 2.5):   # ψ values outside the body [m²/s]
        t = np.linspace(0.01, np.pi - 0.01, 90) if c > m/2 else np.linspace(np.pi + 0.01, 2*np.pi - 0.01, 90)
        rr = (c - m*t/(2*np.pi))/np.sin(t)               # r(θ) from Ur sinθ + mθ/2π = c (6.31)
        xs += list(rr*np.cos(t)) + [np.nan]; ys += list(rr*np.sin(t)) + [np.nan]
    a = m/(2*np.pi)                                      # stagnation distance [m]
    return {"body ψ = m/2": (xb, yb), "streamlines": (np.array(xs), np.array(ys)),
            "stagnation point": (np.array([-a]), np.array([0.0]))}
vals = np.linspace(0.5, 5.0, 20 if not FAST else 10)     # slider values of m/U [m]
fig = slider_figure(half_body_traces, "m/U", vals, unit="m", xlabel="x [m]", ylabel="y [m]",
                    title="The half-body scales with m/U", xrange=[-3, 8], yrange=[-6, 6],
                    modes={"stagnation point": "markers"})
recolor(fig, {"body ψ = m/2": COLORS["ink"], "streamlines": COLORS["blue"], "stagnation point": COLORS["ink"]})
step_titles(fig, [f"The half-body scales with m/U: a = {v/(2*np.pi):.2f} m, h_max = {v/2:.2f} m" for v in vals])
fig.show()                                               # draw it (works on the web page too)
""")
see_read_change(
    see="A body that grows as you drag the slider: its nose moves upstream, its far width grows, the outer streamlines "
        "are pushed aside.",
    read="Every length is proportional to m/U: a = m/2πU and h_max = m/2U (the title), so the shape is the same, only scaled.",
    change="…U doubled at fixed m: the same as halving m/U — a body half as big.")
remind([
    ("show_viz", "`show_viz(\"ch06\", slug)` embeds an interactive explainer; on the web page it fills the window (Ch. 1 P18)."),
], lead="explainers have a Walkthrough, Explore with presets, Explain, Derivation, Equations, Code and Check tabs (Ch. 1 §1.4 explainer guide)")
explainer("superposition_sandbox", "Where does the body come from?",
          "a static figure shows a finished sum. Here you build it: toggle elements, slide their strengths and watch the "
          "stagnation points and the dividing streamline appear, move and close; click a point to see each element's "
          "velocity add up.", "", [
              "Preset 'half-body': read a = 1 m and h_max = 3.14 m in the Explain tab, then halve U — both double.",
              "Add the sink (preset 'Rankine oval'): the body closes and the status says so: Σm = 0.",
              "Preset 'source–sink → doublet': press ▶ and watch the pair collapse into the doublet's circles.",
              "Click a point in the flow: the term bars show the stream's and the source's velocities adding to the total.",
          ])
whatif(r"""
…a sink of the same strength sat 2 m downstream of the source? The source's fluid would be swallowed again, the dividing
streamline would close behind the body, and the drag of this closed body would be exactly zero — the **Rankine oval**
(an exercise idea; `ch06.rankine_oval(1.0, 2*np.pi, 1.0)` gives a half-length $\sqrt3\approx1.732$ m and a half-width 1.307 m).
Closing the body is the step from C05 to C06.
""")

# =====================================================================================================================
# A.3 §6.3 (second half) — R12 R13 R14, C06 (N30 N31, D09), C07 (N32–N38, D10 D11, E2), R15 R16 R17, C08 (N39 N40, D12 D13, E3)
# =====================================================================================================================
nb.recap("R12", "The cylinder = stream + doublet", r"""
A uniform stream U plus a doublet $\mathbf d=-2\pi Ua^2\mathbf e_x$ at the origin gives
$\phi=U\big(r+\frac{a^2}{r}\big)\cos\theta,\ \psi=U\big(r-\frac{a^2}{r}\big)\sin\theta$ *(Eq. 6.33)*; ψ = 0 on r = a for every θ,
so the circle is a streamline, and the body closes because the doublet adds no net source. You met this flow in Ch. 3
(Fig. 3.2); D09 builds it. Our `pf.cylinder(U, a)` and Ch. 3's `cylinder_flow` give the same velocity.
""", where="Ch. 3 §3.3")
nb.code(r"""
print(pf.cylinder(2.0, 0.5).velocity(0.3, 0.8), ch03.cylinder_flow(0.3, 0.8, U=2.0, a=0.5))   # U = 2 m/s, a = 0.5 m: equal
""")
nb.recap("R13", "The cylinder's velocity", r"""
$u_r=U\big(1-\frac{a^2}{r^2}\big)\cos\theta,\ u_\theta=-U\big(1+\frac{a^2}{r^2}\big)\sin\theta$ *(Eq. 6.34)*: on r = a the
radial part vanishes and the fluid slides round at $\lvert u_\theta\rvert=2U\lvert\sin\theta\rvert$.
""", where="Ch. 3 §3.3")
nb.recap("R14", "A cylinder moving through still fluid", r"""
Seen from the fluid far away, the moving cylinder's flow is the cylinder flow minus the stream:
$w=U(z+a^2/z)-Uz=Ua^2/z$ — the instantaneous streamlines of a doublet (Ch. 3's `galilean_frames_cylinder` explainer shows
both frames; open it from the Ch. 3 page). (The complex notation w(z) is C09's; here read it as "the cylinder flow minus
the uniform stream".)
""", where="Ch. 3 §3.3")
nb.code(r"""
print(ch03.cylinder_flow(0.3, 0.8, U=2.0, a=0.5, frame="fluid"),   # Ch. 3: the cylinder seen from the fluid far away
      pf.Doublet((-2*np.pi*2.0*0.25, 0.0)).velocity(0.3, 0.8))     # = a doublet with d = −2πUa² e_x alone
""")

# ---- C06 ------------------------------------------------------------------------------------------------------------
core("C06", "The circular cylinder and d'Alembert's paradox", r"""
Wind presses hard on the front of a chimney. Why does ideal flow predict no net push at all?
""", eqs=("6.35",))
problem(r"""
Stand behind a chimney on a windy day and you are sheltered; the chimney itself feels a strong push downwind. Yet the
ideal-flow solution for a circular cylinder — a flow you have already drawn in Ch. 3 — predicts that the push is exactly
zero. This is d'Alembert's paradox (1752), and understanding why the prediction fails is the reason Ch. 9 exists.
""")
idea("""
      suction (-3q)               pressure on the surface, q = rho U^2 / 2:
 +q  (  O  )  +q     --->         front +q, shoulders -3q, back +q
      suction (-3q)               the back pushes forward exactly as hard as the front pushes back
""", r"""
**The pattern is the same front and back, top and bottom — so every push is cancelled.** D is the drag (the force along
the stream) and L the lift (across it), both per metre of cylinder [N/m].
""")
remind([
    ("net pressure force −∮p n dA", "a pressure p pushes on a surface element along −n (n points out of the body); the net force is −∮p n dA (Ch. 1 P28)."),
    ("trapezoid rule, np.trapezoid", "`np.trapezoid(f, x)` integrates sampled values; on a closed loop the plain sum × Δθ is the periodic version (Ch. 1 P37)."),
], lead="the polar velocities are notes N18, N19 of C02")
D("D09", ref="6.35")
note("N30", "Surface pressure", r"""
Stagnation points (C_p = 1) at θ = 0 and π; the fluid is twice as fast as the stream at the shoulders θ = ±90°, where
C_p = −3. Air (ρ = 1.2 kg/m³), U = 10 m/s: ½ρU² = 60 Pa, so +60 Pa front and back, −180 Pa at the shoulders.
""", equation=EQ["6.35"], ref="6.35")
P("P151", "integrals of sines and cosines over a full period", r"""
Over one full turn, 0 ≤ θ ≤ 2π, the positive and negative halves of sin θ cancel, and so do those of sin³θ, cos θ,
sin θ cos θ and sin²θ cos θ: all integrate to zero. Squares do not cancel:
$\int_0^{2\pi}\sin^2\theta\,d\theta=\int_0^{2\pi}\cos^2\theta\,d\theta=\pi$ (their average is ½). Force integrals round a
circle pick out exactly the terms that survive.
""", code=r"""
th = np.linspace(0, 2*np.pi, 2001)                       # one full turn [rad]
for g in (np.sin(th), np.sin(th)**3, np.sin(th)*np.cos(th), np.sin(th)**2):   # four integrands
    print(round(np.trapezoid(g, th), 6))                 # 0, 0, 0, then π = 3.141593
""")
nb.worked_example("the chimney in numbers", r"""
Air, U = 10 m/s, a = 0.1 m, q = ½ρU² = 60 Pa.

1. Front (θ = 180°): C_p = 1 − 4·0 = 1 → p − p∞ = +60 Pa.
2. Shoulder (θ = 90°): C_p = 1 − 4 = −3 → −180 Pa.
3. Back (θ = 0): +60 Pa.
4. Drag = −∮(p − p∞) cos θ a dθ: the points θ and π − θ have the same p but opposite cos θ, so they cancel in pairs → D = 0.
5. Lift: θ and −θ have the same p and opposite sin θ → L = 0.
""")
nb.code(r"""
f = pf.cylinder(10.0, 0.1)                               # air stream U = 10 m/s past a cylinder a = 0.1 m
print(ch06.cylinder_surface_cp(np.radians([180.0, 135.0, 90.0, 45.0, 0.0])))   # (6.35): 1, −1, −3, −1, 1
circle64 = 0.1*np.exp(1j*np.linspace(0, 2*np.pi, 64, endpoint=False))   # 64 surface points as complex numbers [m]
print(ch06.surface_pressure_force(f, circle64, rho=1.2)) # (D, L) = −∮(p − p∞) n dl by the periodic trapezoid: ≈ (0, 0) N/m
print(pf.normal_velocity_on(f, circle64))                # the surface is a streamline: u·n ≈ 0
""", explain=r"""
1. The cylinder object (stream + doublet, R12).
2. $C_p=1-4\sin^2\theta$ at five angles.
3. `surface_pressure_force` integrates the Bernoulli pressure times the outward normal over 64 surface points (air,
   ρ = 1.2 kg/m³): zero drag and lift, to round-off (compare ½ρU²a = 6 N/m).
4. No flow crosses the surface.
""")
nb.md("**From scratch — 64 pressure pushes added by hand:**")
nb.check_agree(r"""
th = np.linspace(0, 2*np.pi, 64, endpoint=False)          # 64 angles round the cylinder [rad]
p = 0.5*1.2*10.0**2*(1 - 4*np.sin(th)**2)                 # p − p∞ from (6.35) [Pa]
D_hand = -np.sum(p*np.cos(th)*0.1)*(2*np.pi/64)           # −∮(p − p∞) n_x dl, dl = a dθ [N/m]
print(D_hand)                                             # ≈ 0
assert np.isclose(D_hand, ch06.surface_pressure_force(f, circle64, rho=1.2)[0], atol=1e-10)   # = the library's drag
""")
nb.md("Adding the 64 pressure pushes by hand gives zero drag too.")
nb.figure(r"""
n = 241 if not FAST else 161                             # grid points per side
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.9))  # flow and arrows | pressure curves
flow_net(a1, f, xlim=(-0.3, 0.3), ylim=(-0.22, 0.22), n=n, n_levels=17, body=circle(0.1),
         mask_body=lambda X, Y: np.hypot(X, Y) < 0.1)    # streamlines round the cylinder
zs = circle(0.1, 36)                                     # 36 surface points [m]
ps = f.pressure(zs.real, zs.imag, rho=1.2, U=10.0)       # p − p∞ there [Pa]
pressure_arrows(a1, zs, ps, zs/0.1, scale=2.5e-4)        # −(p − p∞) n: orange pushes in, blue pulls out
a1.set_title("air, U = 10 m/s: arrows = −(p − p∞) n", fontsize=10)
beta = np.linspace(0, 180, 361)                          # angle from the front stagnation point [deg] = 180° − θ
lo, hi = separated_cp_band(beta, sep_deg=80.0)           # QUALITATIVE band of a separated high-Re flow (not data)
a2.fill_between(beta, lo, hi, color=COLORS["muted"], alpha=0.3, label="qualitative sketch of a measured high-Re curve — not data")
a2.plot(beta, 1 - 4*np.sin(np.radians(beta))**2, color="k", lw=2, label=r"ideal: $1-4\sin^2\theta$ (6.35)")
a2.axhline(0, color=COLORS["muted"], lw=0.8)             # C_p = 0 (p = p∞)
a2.set_xlabel("angle from the front stagnation point = 180° − θ [deg]"); a2.set_ylabel("$C_p$ [–]")
a2.legend(fontsize=7, loc="lower left")
fig.suptitle("Symmetric pressure, no drag — until the flow separates", fontweight="bold")
savefig(fig, "ch06", "c06_cylinder_pressure"); plt.show()   # save, then draw
""", see="Left, arrows mirror-symmetric front/back and top/bottom: orange pushing in at the front and back, long blue arrows "
         "pulling out at the shoulders. Right, the ideal curve from +1 down to −3 at 90° and back to +1; the grey band "
         "follows it on the front and stays low on the back.",
    read="In ideal flow the back half pushes forward as hard as the front pushes back. A real wake keeps the back at "
         "suction (the grey band is a labelled *qualitative* sketch, not measurements: it follows the ideal curve down to "
         "a shallower minimum near −1, stays there until the separation angle, 80° here, then settles on a flat low wake "
         "pressure), so the front wins: drag.",
    change="…the flow separated earlier (drag the separation angle in the plotly figure below): the low-pressure plateau "
           "widens and the imbalance — the drag — grows.")
nb.plotly(r"""
beta = np.linspace(0, 180, 181)                          # angle from the front [deg]
b = np.radians(beta)                                     # the same in radians
def band_traces(sep):                                    # curves for one separation angle [deg]
    lo, hi = separated_cp_band(beta, sep_deg=sep)        # the qualitative band
    return {"ideal 1 − 4 sin²θ": (beta, 1 - 4*np.sin(b)**2), "band, lower edge": (beta, lo), "band, upper edge": (beta, hi)}
seps = np.linspace(60, 120, 13)                          # separation angles [deg]
fig = slider_figure(band_traces, "separation angle", seps, unit="°", xlabel="angle from the front [deg]",
                    ylabel="C_p [–]", title="Where ideal and real flow part company (qualitative)", yrange=[-3.2, 1.2])
recolor(fig, {"ideal 1 − 4 sin²θ": COLORS["ink"], "band, lower edge": COLORS["muted"], "band, upper edge": COLORS["muted"]})
cd = [np.trapezoid(0.5*sum(separated_cp_band(beta, sep_deg=s))*np.cos(b), b) for s in seps]   # sketch drag coefficient
step_titles(fig, [f"Where ideal and real flow part company (qualitative): sketch C_D = {c:.2f}, ideal C_D = 0" for c in cd])
fig.show()                                               # draw it (works on the web page too)
""", explain=r"""
For each separation angle the band's middle curve is integrated like the drag in the tiny example:
$C_D=D/(\tfrac12\rho U^2\cdot2a)=\int_0^\pi C_p\cos\beta\,d\beta$ with β the angle from the front. The ideal curve gives
exactly 0; the sketch gives a positive drag coefficient that grows as separation moves forward. The band is a
*qualitative* sketch (no dataset is cited here), so read the numbers as an illustration of the mechanism, not a prediction.
""")
note("N31", "Ideal vs real", r"""
A measured pressure distribution on a cylinder at high Reynolds number follows 1 − 4 sin²θ on the front, but behind the
widest point the boundary layer separates and the rear pressure stays low (the book compares with a measured curve; no
dataset is cited here, so we draw only a labelled band). ⚠️ That comparison measures the angle from the upstream
stagnation point, 180° − θ in our convention. Ch. 9 explains separation.
""")
confusion(r"""
"no drag means no pressure". The pressures are large (−180 Pa at the shoulders in the example); they just cancel in
total. Each half of the cylinder feels a big force; the halves fight each other to a draw.
""")
whatif(r"""
…we added circulation? The top would speed up and the bottom slow down; the front–back symmetry survives, the top–bottom
one does not — lift without drag (C07).
""")

# ---- C07 ------------------------------------------------------------------------------------------------------------
core("C07", "Circulation gives lift", r"""
A spinning ball curves; a rotor ship sails on spinning towers. How does circulation turn into a sideways force, and how
big is it?
""", eqs=("6.40",))
problem(r"""
A football struck with spin bends in flight; Flettner's rotor ships used tall spinning cylinders instead of sails. In
both, the air is carried round the body on one side and held back on the other. The ideal-flow model adds a vortex to the
cylinder: the stagnation points slide round, the pressure drops on the fast side — and the net force comes out as a
strikingly simple product.
""")
idea("""
stream -->        + clockwise vortex        =  faster on top, slower below
   (  O  )             (circulation -Gamma)     stagnation points slide DOWN: sin(theta) = -Gamma/(4 pi a U)
                                                pressure lower on top  =>  force UP:  L = rho U Gamma, D = 0
""")
nb.md(r"""
> ⚠️ **Convention callout (numbers).** The book now names a **clockwise** circulation Γ (the flow's circulation is −Γ; its
> footnote says the sign is chosen to give the usual $L=\rho U\Gamma$ *(6.40)*). Our code keeps Γ counterclockwise in
> `pf.Vortex` and asks body functions for `Gamma_cw=` (book) or `Gamma_ccw=` (project). Γ_cw = 2 m²/s ⇔ Γ_ccw = −2 m²/s:
> stagnation points at −9.16° and −170.84° (below the axis), L = +24 N/m for U = 10 m/s, a = 0.1 m, air.
""")
note("N32", "Stream + doublet + a clockwise vortex at the centre", r"""
The a inside the logarithm only adds a constant to ψ, so that ψ = 0 on r = a still.
""", equation=EQ["6.36"], ref="6.36")
remind([
    ("Vieta's formulas (product of roots)", r"the two roots of $r^2-Br+C=0$ multiply to C (Ch. 2 P71)."),
], lead="two solutions of sin θ = s: θ = arcsin s and θ = π − arcsin s; `brentq` is Ch. 3 P108")
D("D10", ref="6.38")
note("N33", "Surface speed", r"""
On top (θ = 90°) the fluid moves at 2U + Γ/2πa, below at 2U − Γ/2πa: 23.18 and 16.82 m/s in the example.
""", equation=EQ["6.37"], ref="6.37")
note("N34", "Stagnation points", r"""
Two on the surface for Γ < 4πaU, one (at the bottom) when Γ = 4πaU, one off the body on the negative y-axis at
$r=\frac1{4\pi U}[\Gamma+\sqrt{\Gamma^2-(4\pi aU)^2}]$ for Γ > 4πaU. U = 10 m/s, a = 0.1 m, Γ = 2 m²/s: sin θ = −0.159,
θ = −9.16° and −170.84°; the critical circulation 4πaU = 12.57 m²/s; at Γ = 6πaU the free point sits at r = 0.2618 m
(its partner inside at 0.0382 m, product a² = 0.01 m²).
""", equation=EQ["6.38"], ref="6.38")
P("P152", "Newton's method for complex zeros", r"""
To find where a complex function f(z) vanishes, start at a guess z₀ and repeat $z\leftarrow z-f(z)/f'(z)$: each step
follows the tangent line to zero, and near a simple root the number of correct digits doubles each time. Here
f = dw/dz (a stagnation point is where u − iv = 0) and f′ = d²w/dz². To find a second root, divide f by (z − z₁)
("deflation") so the first is not found again. (Complex arithmetic in numpy: `1j` is i; C09 has the full primer.)
""", code=r"""
f_ = lambda z: z**2 + 1                                  # zeros at ±i
fp = lambda z: 2*z                                       # its derivative
z = 0.5 + 0.5j                                           # a guess in the upper half plane
for _ in range(6):                                       # six Newton steps
    z = z - f_(z)/fp(z)                                  # follow the tangent to zero
print(z)                                                 # 1j (to machine precision)
""")
note("N35", "Surface pressure", r"""
from Bernoulli (6.18) with the constant p∞ + ½ρU² and the surface speed (6.37):
""", equation=EQ["6.39"], ref="6.39")
note("N36", "The lift integral", r"""
With the outward normal n = e_r and the arc element dl = a dθ, the upward force is (no book number; it sits between
(6.39) and (6.40)) the one below. The minus sign: pressure on the upper surface, where n has an upward component, pushes
the cylinder **down**.
""", equation=r"L=-\int_0^{2\pi}p(r=a,\theta)\,\mathbf n\,dl\cdot\mathbf e_y=-\int_0^{2\pi}p(r=a,\theta)\sin\theta\,a\,d\theta")
D("D11", ref="6.40")
nb.worked_example("a spinning cylinder in wind", r"""
U = 10 m/s, a = 0.1 m, Γ = 2 m²/s (clockwise), air ρ = 1.2 kg/m³.

1. Critical Γ = 4πaU = 12.57 m²/s — we are below it.
2. sin θ = −2/12.566 = −0.159 → θ = −9.16° and −170.84°.
3. Top speed 2U + Γ/2πa = 20 + 3.18 = 23.18 m/s; bottom 20 − 3.18 = 16.82 m/s.
4. Lift L = ρUΓ = 1.2 × 10 × 2 = 24 N per metre of cylinder, upward.
5. Drag 0.
""")
nb.code(r"""
f = pf.cylinder(10.0, 0.1, Gamma_cw=2.0)                 # the book's clockwise Γ = 2 m²/s (our Γ_ccw = −2)
print(ch06.cylinder_stagnation_points(10.0, 0.1, Gamma_cw=2.0))   # closed form (6.38): two points on |z| = 0.1 m
print(np.degrees(np.angle(f.stagnation_points())))       # Newton on dw/dz: the same angles, −170.84° and −9.16°
print(ch06.surface_pressure_force(f, circle64, rho=1.2), ch06.lift_per_span(1.2, 10.0, Gamma_cw=2.0))   # (0, 24) and ρUΓ = 24 N/m
print(f.velocity(0.0, 0.1), pf.cylinder(10.0, 0.1, Gamma_ccw=2.0).velocity(0.0, 0.1))   # top of the body, both signs of Γ
""", explain=r"""
1. The cylinder with the book's clockwise Γ.
2. The closed-form stagnation points of (6.38) and Newton's method agree.
3. The surface-pressure integral gives L = ρUΓ = 24 N/m and D = 0 (air).
4. At the top, the clockwise Γ speeds the flow up to 23.18 m/s; a counterclockwise Γ of the same size slows it to 16.82
   m/s — and gives L = −24 N/m, pushing the cylinder down.
""")
nb.md("**From scratch — the lift integral of note N36 with our own 64-point sum:**")
nb.check_agree(r"""
U, a, G, rho = 10.0, 0.1, 2.0, 1.2                       # stream [m/s], radius [m], clockwise Γ [m²/s], air [kg/m³]
th = np.linspace(0, 2*np.pi, 64, endpoint=False)         # 64 surface angles [rad]
p = 0.5*rho*(U**2 - (-2*U*np.sin(th) - G/(2*np.pi*a))**2)   # (6.39) with p∞ = 0 [Pa]
L_hand = -np.sum(p*np.sin(th)*a)*(2*np.pi/64)            # L = −∮ p sinθ a dθ [N/m]
print(L_hand)                                            # 24.000…
assert np.isclose(L_hand, rho*U*G) and np.isclose(L_hand, ch06.lift_per_span(rho, U, Gamma_cw=G))   # ρUΓ, both ways
""")
nb.md("64 pressure pushes added by hand give 24 N/m to machine precision: ρUΓ.")
note("N38", "Which Γ? Uniqueness needs topology", r"""
In a region with no holes an ideal flow is fixed by its boundary conditions. Round a cylinder the region has a hole: every
member of the family (6.36) has u·n = 0 on r = a and U far away, whatever Γ — the boundary conditions cannot choose. Real
flow chooses, through viscosity at a sharp trailing edge: the Kutta condition (Ch. 14). (Simply connected regions: recap
R04.) The circulation of any loop round the body, $\Gamma=\oint_C\mathbf u\cdot d\mathbf s$ *(Eq. 3.18)*, is the same
for every loop.
""")
nb.code(r"""
chk = ch06.circulation_family_check(10.0, 0.1, [0.0, 2*np.pi, 4*np.pi, 8*np.pi])   # Γ_cw = 0, 2πaU, 4πaU, 8πaU (aU = 1)
print("max |u·n| on the body :", chk["max_normal"])      # ≈ 0 for every Γ: all are valid flows round the cylinder
print("far-field error at 1000a:", chk["far_field"])     # the vortex part Γ/2πR decays like 1/R: every Γ returns to U
print("loop circulation      :", chk["circulation"])     # ∮u·ds = −Γ_cw (counterclockwise loops): 0, −6.28, −12.57, −25.13
""", explain=r"""
**What does this show?** Every member of the family is a valid flow round the cylinder: the normal velocity on the body
is zero to round-off for all four Γ. Far away every one returns to the stream; the leftover error on r = 1000a = 100 m is
exactly the vortex's own speed Γ/2πR (0, 0.01, 0.02, 0.04 m/s). The circulation of counterclockwise loops round the body
is −Γ_cw — the book's Γ is clockwise — whatever the loop's radius. Nothing in the boundary conditions picks Γ.
""")
nb.figure(r"""
n = 201 if not FAST else 141                             # grid points per side
fig, axs = plt.subplots(2, 2, figsize=(8.5, 7.2))        # four regimes
for ax, k in zip(axs.ravel(), (0.0, 0.5, 1.0, 1.5)):     # Γ/(4πaU) = 0, ½, 1, 1.5
    G = k*4*np.pi                                        # clockwise Γ for U = 1 m/s, a = 1 m [m²/s]
    fl = pf.cylinder(1.0, 1.0, Gamma_cw=G)               # the flow
    st = ch06.cylinder_stagnation_points(1.0, 1.0, Gamma_cw=G)   # its stagnation points (closed form)
    flow_net(ax, fl, xlim=(-3, 3), ylim=(-3, 3), n=n, n_levels=25, body=circle(1.0), stag=st,
             mask_body=lambda X, Y: np.hypot(X, Y) < 1.0, title=f"Γ = {k*4:.0f}πaU (clockwise)")
    if k > 1:                                            # a free stagnation point: draw the streamline through it
        Xs_, Ys_ = np.meshgrid(np.linspace(-3, 3, n), np.linspace(-3, 3, n))   # grid [m]
        Ps_ = np.where(np.hypot(Xs_, Ys_) < 1.0, np.nan, fl.psi(Xs_, Ys_))   # ψ outside the body
        ax.contour(Xs_, Ys_, Ps_, levels=[float(fl.psi(st[0].real, st[0].imag))], colors=COLORS["orange"], linewidths=1.8)   # the separatrix
fig.suptitle("The stagnation points slide down, meet, and leave the body", fontweight="bold")
savefig(fig, "ch06", "c07_circulation_regimes"); plt.show()   # save, then draw
""", see="Two stagnation points (orange) at the sides; lower, at −30° and −150°; merged at the bottom; and a free point "
         "below the body where two streamlines cross.",
    read=r"$\sin\theta=-\Gamma/4\pi aU$: at Γ = 2πaU they sit at −30° and −150°; at 4πaU at −90°; above that the root "
         r"$r_+=[\Gamma+\sqrt{\Gamma^2-(4\pi aU)^2}]/4\pi U=2.618a$ (for Γ = 6πaU) lies in the fluid, and a ring of fluid "
         "circulates round the body forever.",
    change="…U doubled at fixed Γ: the points move back up (Γ/4πaU halves) and the lift doubles (ρUΓ).")
nb.plotly(r"""
thd = np.linspace(-180, 180, 181)                        # surface angle θ [deg]
U_, a_ = 10.0, 0.1                                       # the air example [m/s, m]
def cp_traces(k):                                        # curves for Γ/(4πaU) = k
    G = k*4*np.pi*a_*U_                                  # clockwise Γ [m²/s]
    cp = ch06.cylinder_surface_cp(np.radians(thd), U_, a_, Gamma_cw=G)   # surface C_p from (6.39)
    st = ch06.cylinder_stagnation_points(U_, a_, Gamma_cw=G)   # stagnation points (complex)
    on = st[np.isclose(np.abs(st), a_)]                  # only those on the body
    return {"surface C_p": (thd, cp), "Γ = 0 ghost": (thd, 1 - 4*np.sin(np.radians(thd))**2),
            "stagnation points": (np.degrees(np.angle(on)) if on.size else np.array([np.nan]), np.ones(max(on.size, 1)))}
ks = np.linspace(0, 2, 21 if not FAST else 11)           # slider values of Γ/(4πaU)
fig = slider_figure(cp_traces, "Γ/4πaU", ks, xlabel="θ on the surface [deg]", ylabel="C_p [–]",
                    title="Circulation shifts the suction to the top", yrange=[-12, 1.5],
                    modes={"stagnation points": "markers"})
recolor(fig, {"surface C_p": COLORS["ink"], "Γ = 0 ghost": COLORS["muted"], "stagnation points": COLORS["orange"]},
        dashes={"Γ = 0 ghost": "dash"})
titles = []                                              # one title per slider step: lift and regime
for k in ks:
    s = ch06.cylinder_circulation_state(U_, a_, Gamma_cw=k*4*np.pi*a_*U_, rho=1.2)   # air
    titles.append(f"Circulation shifts the suction to the top: L = ρUΓ = {s['L_KJ']:.1f} N/m, D = 0 — {s['regime']}")
step_titles(fig, titles)
fig.show()                                               # draw it (works on the web page too)
""")
see_read_change(
    see="The C_p curve against θ: as Γ grows, the dip at +90° (the top) deepens and the one at −90° (the bottom) fills in; "
        "the orange stagnation markers slide toward −90° and vanish from the body above Γ/4πaU = 1.",
    read="Lower pressure on top than below gives the upward lift L = ρUΓ printed in the title; the curve stays symmetric "
         "about ±90° (front–back), which is why the drag stays zero.",
    change="…the Γ sign flipped (counterclockwise): the whole picture mirrors top-to-bottom and the lift points down.")
explainer("cylinder_circulation_lift", "How does spin turn into lift?",
          "four snapshots sample a continuous story: dragging Γ moves the stagnation points round the body, makes them "
          "meet at 4πaU and leave it, while the pressure arrows and the force bars move with them — the lift bar lands "
          "on ρUΓ and the drag bar never leaves zero.", "", [
              "Drag Γ/4πaU from 0 to 1: the two stagnation points slide down and meet.",
              "Keep going to 1.5: a free stagnation point leaves the body — find its radius in Explain.",
              "Switch the convention chips to 'project Γ (ccw)': the same picture, the numbers change sign.",
              "Open Derivation → D11 step 7: the only term of the pressure that survives the lift integral lights up on the C_p curve.",
          ])
note("N37", "History and spin", r"""
Kutta and Zhukhovsky found $L=\rho U\Gamma$ independently just after 1900. A rotating cylinder develops circulation
through its boundary layer; at high spin the real flow resembles the Γ > 4πaU pattern, at low spin it separates. The
Magnus effect on spinning balls is mostly delayed separation (Ch. 9 §9.9); how a sharp edge fixes Γ is Ch. 14.
""")
whatif(r"""
…the body were an ellipse or an airfoil instead of a circle? The lift would still be ρUΓ — the shape only changes *which*
Γ nature picks. C10 proves it for any body; C11 builds the ellipse.
""")

# ---- recaps for C08 ---------------------------------------------------------------------------------------------------
nb.recap("R15", "Images for circles", r"""
A circle can be made a streamline by more than one image (a vortex outside a cylinder needs an image vortex at the
inverse point and one at the centre), and images move when the flow is unsteady; Ch. 5's `circle_image_system` does it.
Milne-Thomson's circle theorem, $w(z)+\overline{w(a^2/\bar z)}$, is the general rule (our addition, `pf.circle_theorem`;
w is C09's complex potential).
""", where="Ch. 5 §5.7")
nb.recap("R16", "A vortex beside a wall", r"""
A vortex of strength −Γ (clockwise) at (h, 0) beside the wall x = 0 is modelled by the vortex plus an image of the
opposite sign at (−h, 0); then u = 0 on the wall.
""", where="Ch. 5 §5.7")
nb.recap("R17", "Its drift", r"""
The vortex moves with the velocity its image induces at its centre (its own is taken as zero): $d\xi_x/dt=0$,
$d\xi_y/dt=\Gamma/4\pi\xi_x$, so $\boldsymbol\xi(t)=(h,\ \Gamma t/4\pi h)$ — along the wall at Γ/4πh (0.0796 m/s for
Γ = 1 m²/s, h = 1 m). Ch. 5's walls are horizontal (y = 0), so we turn the picture by 90° counterclockwise: the wall
x = 0 becomes y = 0, the vortex at (h, 0) moves to (0, h), and the drift along +y becomes a drift along −x.
""", where="Ch. 5 §5.7")
nb.code(r"""
sol = bs.point_vortex_evolve(np.array([[0.0], [1.0]]), np.array([-1.0]), np.array([0.0, 10.0]), boundary="wall")   # Ch. 5
print(sol[-1].ravel())                                   # after 10 s: x = −0.7958 m (= −Γt/4πh), y = 1.0 m
""", explain=r"""
`point_vortex_evolve(positions, strengths, times, boundary="wall")` integrates the motion of point vortices above the
wall y = 0 with their images (Ch. 5 §5.7): a clockwise vortex (Γ = −1 m²/s in our counterclockwise sign) at height 1 m
drifts 0.7958 m in 10 s — the drift Γ/4πh = 0.0796 m/s along the wall.
""")

# ---- C08 ------------------------------------------------------------------------------------------------------------
core("C08", "Images, and what the wall feels as a vortex passes (Example 6.1)", r"""
An eddy sweeps along a wall. What does a pressure sensor on the wall record as it passes?

*In one line:* $\frac{p(0,0,t)-p_\infty}{\rho}=\frac{\Gamma^2}{4\pi^2}\frac{(\Gamma t/4\pi h)^2-h^2}{\big((\Gamma t/4\pi h)^2+h^2\big)^2}$ (Example 6.1).
""")
problem(r"""
An eddy rolls along a seabed, a harbour wall or the side of a building; a pressure sensor mounted in the wall records a
blip. The wall is solid, so the flow must slide along it — and a mirror trick, the image vortex, makes that automatic.
With the image in place, unsteady Bernoulli tells us the whole pressure signal: suction as the eddy arrives, a small
over-pressure after it has passed.
""")
idea("""
 image +G (ccw)  |  (cw) vortex -G       mirror images: vortex -> opposite sign, source -> same sign
   (-h, s)       |    (h, s)             the pair drifts along the wall at G/(4 pi h) (each pushes the other)
            wall x = 0                   p at the wall = -rho dphi/dt (unsteady, +) - rho v^2/2 (speed, -)
""", "(G stands for Γ.)")
nb.md(r"""
> 🔁 **Glosses used below.** *Odd and even functions under a mirror:* $\psi(x,-y)=-\psi(x,y)$ makes ψ = 0 on y = 0; an
> even φ, $\phi(x,-y)=\phi(x,y)$, has ∂φ/∂y = 0 there. *Derivative of an arctangent with a moving argument:*
> $\frac{d}{dt}\tan^{-1}(Y/X)=\frac{X\dot Y}{X^2+Y^2}$ at fixed X (the chain rule, P49). *Unsteady Bernoulli* with
> gravity absorbed: $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+gz+\frac p\rho=\text{const}$ *(Eq. 4.83)*,
> with ∂φ/∂t at a fixed point (recap R05).
""")
D("D12")
note("N39", "The rules", r"""
Wall y = 0: vortex images flip sign, $\psi_2=\psi_1(x,y)-\psi_1(x,-y)$; source images keep it,
$\phi_2=\phi_1(x,y)+\phi_1(x,-y)$ — two different conditions (ψ₂ = 0, ∂φ₂/∂y = 0) that both make the wall impermeable.
""")
nb.code(r"""
els = pf.mirror([pf.Vortex(1.0, 0.5+1.0j), pf.Source(2.0, -1.0+0.7j)], wall="y=0")   # a vortex and a source + their images
fw = pf.Flow(els)                                        # the flow above the wall y = 0
xw = np.linspace(-5, 5, 1001)                            # 1001 points along the wall [m]
print(np.max(np.abs(fw.velocity(xw, 0*xw)[1])))          # max |v| on the wall ≈ 0: nothing crosses it
""")
note("N40", "A source by a wall = two sources", r"""
Equal sources at x = ±a. With the tangent addition formula ($\tan(A+B)=\frac{\tan A+\tan B}{1-\tan A\tan B}$) the
streamlines are $x^2-y^2-2xy\cot(2\pi\psi/m)=a^2$. The axes are streamlines and the origin a stagnation point, so the
same picture is two sources, one source beside a wall (the y-axis), or a slit flow into a right-angled corner:
""", equation=r"\psi=\frac m{2\pi}\Big[\tan^{-1}\Big(\frac y{x+a}\Big)+\tan^{-1}\Big(\frac y{x-a}\Big)\Big],\quad \phi=\frac m{2\pi}\Big[\ln\sqrt{(x+a)^2+y^2}+\ln\sqrt{(x-a)^2+y^2}\Big]", ref="6.41")
nb.code(r"""
f2 = ch06.two_sources(2*np.pi, 1.0)                      # sources m = 2π m²/s at x = ±1 m
xq = np.linspace(1.2, 4, 5)                              # five x positions [m]
yq = ch06.two_source_streamline(np.pi/2, 2*np.pi, 1.0, xq)   # y ≥ 0 on the curve x² − y² − 2xy cot(2πψ/m) = a², ψ = π/2
print(np.mod(f2.psi(xq, yq), 2*np.pi))                   # ψ = π/2 at all five points (mod m, where a branch is crossed)
""")
D("D13")
nb.worked_example("an eddy by a harbour wall", r"""
Water ρ = 1000 kg/m³, Γ = 1 m²/s, h = 1 m.

1. Drift Γ/4πh = 1/12.566 = 0.0796 m/s.
2. At t = 0 (eddy level with the sensor): $\frac{p-p_\infty}{\rho}=\frac{1}{39.48}\cdot\frac{0-1}{1}$ → p − p∞ = −25.33 Pa:
   suction. Split: the unsteady part −ρ∂φ/∂t = +25.33 Pa, the speed part −½ρv² = −50.66 Pa.
3. The signal crosses zero when Γt/4πh = h, i.e. t = 4πh²/Γ = 12.57 s (the eddy is 1 m past the sensor, at 45°).
4. The largest over-pressure comes at Γt/4πh = √3 h, t = 21.77 s: p − p∞ = ρΓ²/(32π²h²) = +3.17 Pa.
5. Then it fades as t⁻².
""")
nb.code(r"""
for t in (0.0, 4*np.pi, 4*np.sqrt(3)*np.pi, 60.0):       # four instants [s]
    r = ch06.example_6_1(t)                              # Example 6.1 in closed form (Γ = 1 m²/s, h = 1 m, water)
    print(f"t = {t:6.2f} s: ξ_y = {r['xi_y']:.3f} m, p − p∞ = {r['p_origin']:+.4f} Pa "
          f"(unsteady {r['unsteady']:+.3f}, speed {r['speed_part']:+.3f})")
print(ch06.example_6_1(5.0, route="numeric")["p_origin"], ch06.example_6_1(5.0)["p_origin"])   # two routes, one number [Pa]
""", explain=r"""
1. The closed form of Example 6.1 and its split into the unsteady part −ρ∂φ/∂t (positive) and the speed part −½ρv²
   (negative): −25.33 Pa at t = 0, zero at 12.57 s, +3.166 Pa at 21.77 s, fading after.
2. The numeric route moves a `Vortex` and its image, takes ∂φ/∂t at the origin by a central difference in time and
   applies Bernoulli — an independent route that agrees to about seven digits.
""")
nb.md("**From scratch — Bernoulli by hand at one instant:**")
nb.check_agree(r"""
Gm, h, rho = 1.0, 1.0, 1000.0                            # Γ [m²/s], distance from the wall [m], water [kg/m³]
phi0 = lambda t: (Gm/(2*np.pi))*(np.arctan2(0 - Gm*t/(4*np.pi*h), 0 + h) - np.arctan2(0 - Gm*t/(4*np.pi*h), 0 - h))   # φ at the origin (D13 step 4)
dt = 1e-4                                                # time step for the central difference [s]
dphidt = (phi0(5 + dt) - phi0(5 - dt))/(2*dt)            # ∂φ/∂t at t = 5 s [m²/s²]
v = Gm*h/(np.pi*(h**2 + (Gm*5/(4*np.pi*h))**2))          # the wall velocity at the origin at t = 5 s [m/s]
p = -rho*(dphidt + 0.5*v**2)                             # unsteady Bernoulli, p∞ = 0 [Pa]
print(p)                                                 # −15.89 Pa
assert np.isclose(p, ch06.example_6_1(5.0)["p_origin"], rtol=1e-6)   # = the library's pressure
""")
nb.md("Bernoulli by hand at one instant gives the library's pressure. (`np.arctan2` keeps the angle continuous here.)")
nb.figure(r"""
tt = np.linspace(-60, 60, 601)                           # time [s]; negative t = the eddy approaching from below
res = [ch06.example_6_1(t) for t in tt]                  # the closed form at every time
P0 = np.array([r["p_origin"] for r in res]); Pu = np.array([r["unsteady"] for r in res]); Ps = np.array([r["speed_part"] for r in res])   # total, unsteady and speed parts [Pa]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.8), gridspec_kw=dict(width_ratios=[1, 1.7]))  # geometry | signal
yv = np.linspace(-6, 6, 50)                              # the vortex path [m]
a1.plot(1 + 0*yv, yv, color=COLORS["amber"], lw=2, label="vortex path (h, Γt/4πh)")   # the vortex moves along x = h
a1.plot(-1 + 0*yv, yv, color=COLORS["amber"], lw=2, alpha=0.35, label="image path")   # its mirror image behind the wall
a1.axvline(0, color="k", lw=3); a1.fill_betweenx([-6, 6], -1.6, 0, color="#e8ebf2")   # the wall and the "other side"
tc = 10.0                                                # instant of the snapshot [s]
fv = pf.Flow(pf.mirror([pf.Vortex(-1.0, 1.0 + 1j*tc/(4*np.pi))], wall="x=0"))   # clockwise vortex + image at t = 10 s
Xg, Yg = np.meshgrid(np.linspace(0.02, 3, 121), np.linspace(-6, 6, 241))   # the fluid side x > 0 [m]
Pv = fv.psi(Xg, Yg)                                      # ψ of the vortex and its image [m²/s]
a1.contour(Xg, Yg, Pv, levels=np.nanpercentile(Pv, np.linspace(3, 97, 17)), colors=COLORS["ink"], linewidths=0.6)   # streamlines at t = 10 s
a1.plot([0], [0], "s", color=COLORS["orange"], ms=8, label="pressure sensor")   # the sensor at the origin
a1.set_xlim(-1.6, 3); a1.set_ylim(-6, 6); a1.set_aspect("equal"); a1.set_xlabel("x [m]"); a1.set_ylabel("y [m]")   # frame and axes with units
a1.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1); a1.set_title("streamlines at t = 10 s", fontsize=10)   # legend and title
a2.plot(tt, P0, color="k", lw=2.2, label="p − p∞ at the sensor")   # the sensor signal
a2.plot(tt, Pu, color=COLORS["orange"], ls="--", label="unsteady part −ρ∂φ/∂t")   # always positive here
a2.plot(tt, Ps, color=COLORS["blue"], ls="--", label="speed part −½ρv²")   # always suction
for t_, lab in ((0.0, "−25.3 Pa"), (4*np.pi, "0"), (4*np.sqrt(3)*np.pi, "+3.17 Pa")):   # marked instants
    a2.plot([t_, -t_], [ch06.example_6_1(t_)["p_origin"]]*2, "o", color=COLORS["ink"], ms=4)   # mark ±t
    a2.annotate(lab, (t_, ch06.example_6_1(t_)["p_origin"]), (t_ + 3, ch06.example_6_1(t_)["p_origin"] - 6), fontsize=8)   # its value
a2.axhline(0, color=COLORS["muted"], lw=0.8); a2.set_xlabel("time t [s]"); a2.set_ylabel("p − p∞ [Pa]")   # p = p∞ line and axes
a2.legend(fontsize=7); a2.set_title("Γ = 1 m²/s, h = 1 m, water", fontsize=10)   # legend and the case
fig.suptitle("Suction as it passes, a small push before and after", fontweight="bold")   # the message
savefig(fig, "ch06", "c08_wall_pressure"); plt.show()    # save, then draw
""", see="Left, the vortex path parallel to the wall, its faint image behind the wall and closed streamlines round the "
         "vortex; right, a deep negative dip at t = 0 flanked by two small positive humps.",
    read="The speed part (blue) is always suction; the unsteady part (orange) is always positive here and wins when the "
         "eddy is far away — it decays like t⁻², the speed part like t⁻⁴. The formula is even in t: approach and departure "
         "look alike.",
    change="…h halved: the dip deepens ×4 (the Γ²/h² scaling) and the whole signal runs 4× faster (time scale h²/Γ).")
nb.figure(r"""
f2 = ch06.two_sources(2*np.pi, 1.0)                      # sources m = 2π m²/s at x = ±1 m
Xg, Yg = np.meshgrid(np.linspace(-3, 3, 241 if not FAST else 161), np.linspace(0.005, 3, 121 if not FAST else 81))   # upper half [m]
fig, ax = plt.subplots(figsize=(7, 3.6))                 # one panel
ax.contour(Xg, Yg, f2.psi(Xg, Yg), levels=np.arange(1, 16)*np.pi/8, colors=COLORS["teal"], linewidths=0.8)   # ψ = kπ/8, includes π/4, π/2, 3π/4
for psi_v in (np.pi/4, np.pi/2, 3*np.pi/4):              # three streamlines of (6.41)'s algebraic form, right half
    xq = np.linspace(0.0, 3, 600)                        # x positions [m] (NaN where the curve has no point)
    ax.plot(xq, ch06.two_source_streamline(psi_v, 2*np.pi, 1.0, xq), color="k", ls="--", lw=1.4)
ax.axvline(0, color="k", lw=3); ax.axhline(0, color="k", lw=1)   # the y-axis is a wall; the x-axis a streamline
ax.plot([-1, 1], [0, 0], "o", color=COLORS["teal"], ms=7)   # the two sources
ax.set_xlim(-3, 3); ax.set_ylim(0, 3); ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
ax.set_title("Two sources = one source beside a wall", fontsize=10)
savefig(fig, "ch06", "c08_two_sources"); plt.show()      # save, then draw
""", see="Teal streamlines leaving both sources; the dashed black curves (the algebraic form of (6.41)) lie exactly on "
         "three of them; nothing crosses the y-axis.",
    read="By symmetry the y-axis is a streamline, so the right half is a source beside a wall — the source's image keeps "
         "its sign. The origin is a stagnation point where the two outflows meet.",
    change="…the left source made a sink: the image rule would be broken — the y-axis would no longer be a streamline "
           "(flow would cross it from the source to the sink).")
nb.animation(r"""
nfr = 60 if not FAST else 30                             # number of frames
times = np.linspace(-40, 40, nfr)                        # time [s]
yw = np.linspace(-6, 6, 201)                             # wall points [m]
Pw = np.array([ch06.example_6_1_wall_pressure(yw, t) for t in times])   # p − p∞ along the wall at every time [Pa]
Pt = np.array([ch06.example_6_1(t)["p_origin"] for t in times])         # the sensor signal [Pa]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6))     # wall | sensor trace
sc = a1.scatter(0*yw, yw, c=Pw[0], cmap="coolwarm_r", vmin=-26, vmax=26, s=18)   # blue = suction, red = over-pressure
vort, = a1.plot([1], [0], "o", color=COLORS["amber"], ms=10)   # the vortex
img, = a1.plot([-1], [0], "o", color=COLORS["amber"], ms=10, alpha=0.3)   # its image
a1.axvline(0, color="k", lw=1); a1.set_xlim(-1.6, 1.6); a1.set_ylim(-6, 6); a1.set_xlabel("x [m]"); a1.set_ylabel("y [m]")
fig.colorbar(sc, ax=a1, label="p − p∞ on the wall [Pa]")
a2.plot(times, Pt, color=COLORS["muted"], lw=1)           # the whole signal, faint
dot, = a2.plot([times[0]], [Pt[0]], "o", color="k")       # where we are now
a2.set_xlabel("t [s]"); a2.set_ylabel("p − p∞ at the sensor [Pa]"); a2.set_ylim(-28, 6)
def update(i):                                           # draw frame i
    s = times[i]/(4*np.pi)                               # vortex height ξ_y = Γt/4πh [m]
    sc.set_array(Pw[i]); vort.set_data([1], [s]); img.set_data([-1], [s])   # recolour the wall, move the pair
    dot.set_data([times[i]], [Pt[i]])                    # move the dot on the trace
    a1.set_title(f"t = {times[i]:+.0f} s", fontsize=10)
show_animation(animate(update, frames=nfr, fig=fig, interval=80), player="video")   # smooth MP4
""")
see_read_change(
    see="The vortex and its image slide up together; a blue (suction) patch on the wall travels with them, with faint red "
        "(over-pressure) patches leading and trailing; the black dot traces the sensor signal.",
    read="The sensor at the origin feels the travelling pattern as the dip-and-humps signal of the static figure.",
    change="…the eddy twice as strong (Γ = 2 m²/s): it drifts twice as fast and every pressure is four times larger.")
explainer("vortex_wall_images", "What does a wall feel as an eddy passes?",
          "the result is a time series tied to a moving vortex: scrubbing time links the vortex position to the sensor "
          "trace and to the pressure along the whole wall, and the source mode shows the other image rule on the same stage.",
          "", [
              "Press ▶ and watch the blue suction patch ride along the wall with the eddy.",
              "Preset 't = 4πh²/Γ': the signal crosses zero — the terms panel shows the unsteady and speed parts cancelling.",
              "Halve h: the dip is four times deeper and four times shorter.",
              "Switch to 'source by a wall': the image keeps its sign, and the wall still has no flow through it.",
          ])
whatif(r"""
…the wall were the sea floor and the eddy a storm-driven vortex in the ocean? The same image model gives the
bottom-pressure signal that ocean-bottom pressure recorders see as eddies pass — a small signal (tens of pascals here)
compared with the tides, which is why it is averaged out or modelled (Ch. 13 adds rotation and stratification).
""")

# =====================================================================================================================
# A.4 §6.4 — C09 (N41–N52, D14 D15, E4)
# =====================================================================================================================
nb.section("6.4", "Complex Potential", intro=r"""
**What is this section about?** φ and ψ are two halves of one complex function w(z) of the complex position z = x + iy.
If w has a derivative that does not depend on the direction you take it in, its real and imaginary parts are
automatically a potential and a stream function of an ideal flow, and the derivative is the velocity — with v flipped:
dw/dz = u − iv. Every element of §6.3 becomes one line.
""")
core("C09", "The complex potential and the complex velocity", r"""
Why does a single complex function carry a whole flow — and why is the velocity u − iv, not u + iv?
""", eqs=("6.42", "6.45"))
problem(r"""
So far every flow needed two functions, one for streamlines and one for equipotentials, and every new body meant new
algebra. Engineers designing wings in the 1900s found a shortcut: write the position as one complex number and the flow
as one complex function. Then any "nice" function of z is a flow, forces become a single contour integral (C10), and
maps between planes carry flows onto new shapes (C11).
""")
idea("""
z = x + iy  --->  w(z) = phi + i psi        Re w: equipotentials (dashed)    Im w: streamlines (solid)
derivative along x  =  derivative along iy   <=>   Cauchy-Riemann   <=>   phi _|_ psi and both harmonic
dw/dz = u - iv   (the velocity, mirrored in the x-axis)
""")
remind([
    ("Euler's formula e^{iθ} = cos θ + i sin θ", r"$e^{i\theta}=\cos\theta+i\sin\theta$, with $i^2=-1$ (Ch. 1 P45)."),
    ("complex conjugate", r"$z^*=\bar z=x-iy$; $zz^*=\lvert z\rvert^2$ (Ch. 2 P81)."),
])
P("P153", "the complex plane in numpy", r"""
A complex number z = x + iy is a point (x, y); its modulus |z| = √(x² + y²) is the distance from 0 and its argument
arg z = atan2(y, x) the angle, so $z=re^{i\theta}$ (Euler's formula). Multiplying two complex numbers multiplies their
moduli and **adds their angles** — a rotation plus a stretch. The conjugate z* = x − iy mirrors in the x-axis. numpy:
`1j`, `np.abs`, `np.angle`, `np.conj`; arrays of complex numbers work like any other array.
""", code=r"""
z = 1 + 1j                                               # the point (1, 1)
print(np.abs(z), np.degrees(np.angle(z)))                # modulus 1.414 and angle 45°
q = z * 2j                                               # times 2i: stretch ×2, turn +90°
print(q, np.abs(q), np.degrees(np.angle(q)), np.conj(z)) # (-2+2j) 2.828 135.0 (1-1j)
""")
note("N41", "The complex position", r"""
in Cartesian and polar form: r = (x² + y²)^½, θ = tan⁻¹(y/x) (numpy: `np.abs`, `np.angle`).
""", equation=EQ["6.43"], ref="6.43")
P("P154", "complex derivative and analytic functions", r"""
The derivative of w(z) is the limit of $[w(z+\delta z)-w(z)]/\delta z$ as the small complex step δz shrinks — but δz can
point in any direction in the plane. If the limit is the same for every direction, w is **analytic** there.
Polynomials, exp, log (away from its cut) and 1/z (away from 0) are analytic; z* is not: along x its quotient is 1,
along iy it is −1.
""", code=r"""
w = lambda z: z**2                                       # an analytic function
z0, h = 1 + 1j, 1e-6                                     # a point and a small step
print((w(z0 + h) - w(z0))/h, (w(z0 + 1j*h) - w(z0))/(1j*h))   # both ≈ 2+2j: one derivative
s = np.conj                                              # z* is not analytic
print((s(z0 + h) - s(z0))/h, (s(z0 + 1j*h) - s(z0))/(1j*h))   # 1 and −1: no derivative
""")
remind([
    ("Schwarz's theorem", r"mixed partial derivatives of a smooth function commute: $\partial^2\psi/\partial x\partial y=\partial^2\psi/\partial y\partial x$ (Ch. 4 P121)."),
])
D("D14", ref="6.44")
note("N42", "Cauchy–Riemann", r"""
They make the φ- and ψ-lines cross at right angles — except where w or dw/dz is zero or infinite (stagnation points and
singular points).
""", equation=EQ["6.44"], ref="6.44")
note("N43", "The complex velocity", r"""
Applying the Cauchy–Riemann conditions to it gives back continuity $\partial u/\partial x+\partial v/\partial y=0$
*(6.2)* and irrotationality $\partial v/\partial x-\partial u/\partial y=0$ *(6.9)*, and Laplace's equation for both φ and
ψ: *any analytic function of z is a plane ideal flow.* ⚠️ The book writes "Laplace equations for φ and ψ, (6.5) and
(6.12), respectively" — it is the other way round ((6.5) is ψ's, (6.12) is φ's); and "(6.43) ensures the equality of the
u, v components" means the Cauchy–Riemann conditions (6.44).
""", equation=EQ["6.45"], ref="6.45")
nb.worked_example("w = z² at z = 1 + i", r"""
1. $dw/dz=2z=2+2i$.
2. So $u-iv=2+2i$: u = 2 m/s and **v = −2 m/s** (the sign flips).
3. Check with ψ = Im z² = 2xy (the corner flow (6.24) with A = 1): $u=\partial\psi/\partial y=2x=2$ ✓,
   $v=-\partial\psi/\partial x=-2y=-2$ ✓.
4. Cauchy–Riemann: φ = x² − y², $\partial\phi/\partial x=2x=2=\partial\psi/\partial y$ ✓,
   $\partial\phi/\partial y=-2y=-2=-\partial\psi/\partial x$ ✓.
""")
P("P155", "complex logarithm, powers and branch cuts", r"""
Because θ is only defined up to multiples of 2π, $\ln z=\ln r+i\theta$ has many values; a program picks θ in a range
(numpy: (−π, π]), so ln z jumps by 2πi across a **branch cut** (numpy's runs along the negative real axis). Powers inherit
it: $z^n=e^{n\ln z}=r^ne^{in\theta}$. In a flow we rotate the cut so it lies inside a body or outside the fluid we draw —
then every quantity we look at is smooth (the `cut_angle` of C03 and C05).
""", code=r"""
z = np.array([-1 + 1e-9j, -1 - 1e-9j])                   # just above and below the negative real axis
print(np.log(z).imag)                                    # +π and −π: a jump of 2π across the cut
print(np.sqrt(z))                                        # ±i: the square root jumps too
print(np.log(z*np.exp(-1j*np.pi/2)).imag + np.pi/2)      # cut turned to the −y axis: no jump at z = −1
""")
D("D15", ref="6.46")
note("N44", "Flow in a corner of angle α = π/n", r"""
Its velocity $dw/dz=nAz^{n-1}=(A\pi/\alpha)z^{(\pi-\alpha)/\alpha}$ is zero at the corner when α < π (a stagnation point)
and infinite when α > π. n = 2 is $\psi=2Axy$ *(6.24)* / $\phi=A(x^2-y^2)$ *(6.27)*; n = ½ is the flow round the edge of
a flat plate. ⚠️ Ch. 9's wedge flows U ∝ x^m start from exactly this outer flow.
""", equation=r"w(z)=Az^n=A(re^{i\theta})^n=Ar^n(\cos n\theta+i\sin n\theta),\quad n\ge\tfrac12", ref="6.46")
nb.code(r"""
for n in (4, 2, 1, 2/3, 0.5):                            # corner exponents
    print(f"n = {n:.3f}: corner angle α = π/n = {np.degrees(np.pi/n):5.0f}°, speed ∝ r^{ch06.corner_speed_exponent(n):+.3f}")
""")
nb.md(r"""
📝 **Note.** **Every element in one line** `N45` `N46` `N47` `N48` `N49` `N50` `N51` — each row is one line of algebra
with $\ln(z-z')=\ln r'+i\theta'$ and $\frac1{z-z'}=\frac{(x-x')-i(y-y')}{r'^2}$:

| element | complex potential | its imaginary part |
|---|---|---|
| vortex (⚠️ counterclockwise Γ; φ = Γθ′/2π jumps by Γ across the cut — the circulation) | $w=-\frac{i\Gamma}{2\pi}\ln(z-z')=\frac{\Gamma}{2\pi}\theta'-i\frac{\Gamma}{2\pi}\ln r'$ *(6.47)* | $\psi=-\frac{\Gamma}{2\pi}\ln r'$, the ψ of (6.8) |
| source | $w=\frac m{2\pi}\ln(z-z')=\frac m{2\pi}\ln r'+i\frac{m\theta'}{2\pi}$ *(6.48)* | $\psi=m\theta'/2\pi$ |
| doublet (⚠️ scalar d = the dipole −d e_x) | $w=\frac{d}{2\pi(z-z')}$ *(6.49)* | Re w is the φ of (6.29) |
| half-body | $w=Uz+\frac m{2\pi}\ln z$ *(6.50)* | $\psi=Uy+\frac m{2\pi}\theta$, the ψ of (6.31) |
| cylinder | $w=U\big(z+\frac{a^2}z\big)$ *(6.51)* | $\psi=U(r-a^2/r)\sin\theta$ of (6.33) |
| cylinder with a **clockwise** Γ | $w=U\big(z+\frac{a^2}z\big)+\frac{i\Gamma}{2\pi}\ln(z/a)$ *(6.52)* | ψ of (6.36) |
| two sources | $w=\frac m{2\pi}\ln\big(\frac{z^2-a^2}{a^2}\big)$ *(6.53)* (ln A + ln B = ln AB up to 2πi) | ψ of (6.41), mod m |
""")
nb.code(r"""
zp = np.array([1+1j, -0.5+2j, 2-0.3j])                   # three points as complex numbers [m]
flows = {"z²": pf.Corner(1.0, n=2), "vortex": pf.Vortex(2*np.pi), "source": pf.Source(2*np.pi),
         "doublet": pf.Doublet.from_book_scalar(2*np.pi), "cylinder": pf.cylinder(1.0, 1.0, Gamma_cw=2.0)}   # five flows
for k, e in flows.items():                               # for each flow …
    cr = max(np.max(np.abs(ch06.cauchy_riemann_residual(e.w, z0))) for z0 in zp)   # (φ_x − ψ_y, φ_y + ψ_x) by differences
    u, v = e.velocity(zp.real, zp.imag)                  # the element's own velocity [m/s]
    print(f"{k:9s} Cauchy–Riemann residual {cr:.1e}   max|dw/dz − (u − iv)| = {np.max(np.abs(e.dwdz(zp) - (u - 1j*v))):.1e}")
print("z* (not a flow):", ch06.cauchy_riemann_residual("conj", 1 + 1j))   # (2, 0): φ_x − ψ_y = 1 − (−1)
""", explain=r"""
1. Five flows written as complex functions (`w`, `dwdz` methods).
2. `cauchy_riemann_residual` measures $\phi_x-\psi_y$ and $\phi_y+\psi_x$ by difference quotients along x and along iy:
   zero to difference-quotient accuracy for every flow.
3. dw/dz equals u − iv computed from the element's own velocity (to round-off).
4. z* fails: its first residual is 2 — it is not analytic, so it is not a flow.
""")
nb.md("**From scratch — the two difference quotients by hand:**")
nb.check_agree(r"""
c = pf.Corner(1.0, n=2)                                  # w = z²
z0, h = 1 + 1j, 1e-6                                     # a point and a small step
qx = (c.w(z0 + h) - c.w(z0))/h                           # quotient along x
qy = (c.w(z0 + 1j*h) - c.w(z0))/(1j*h)                   # quotient along iy
print(qx, qy, c.dwdz(z0))                                # all ≈ 2+2j
assert np.allclose([qx, qy], c.dwdz(z0), atol=1e-5)      # one derivative, whatever the direction
qx_c = (np.conj(z0 + h) - np.conj(z0))/h                 # the same for z*: along x …
qy_c = (np.conj(z0 + 1j*h) - np.conj(z0))/(1j*h)         # … and along iy
assert not np.isclose(qx_c, qy_c)                        # 1 ≠ −1: no derivative
""")
nb.md("Two directions, one derivative — that is what \"analytic\" means; z* gives two different answers.")
nb.figure(r"""
n = 201 if not FAST else 141                             # grid points per side
fig, axs = plt.subplots(1, 4, figsize=(12, 3.2))         # four corners
for ax, nn in zip(axs, (2.0, 1.0, 2/3, 0.5)):            # corner exponents
    el = pf.Corner(1.0, n=nn)                            # w = z^n with its cut outside the wedge
    X_, Y_ = np.meshgrid(np.linspace(-2, 2, n), np.linspace(-2, 2, n))   # grid [m]
    Z = X_ + 1j*Y_                                       # as complex numbers
    th_ = np.mod(np.angle(Z), 2*np.pi)                   # polar angle in [0, 2π)
    out = (th_ > np.pi/nn + 1e-9) | ((np.abs(Y_) < 0.03) & (X_ > 0))   # outside the wedge, or on the lower wall
    Pc = np.where(out, np.nan, el.psi(X_, Y_)); Fc = np.where(out, np.nan, el.phi(X_, Y_))   # ψ and φ in the fluid
    sp_ = np.where(out, np.nan, np.abs(el.dwdz(Z)))      # speed |dw/dz| [m/s]
    ax.pcolormesh(X_, Y_, np.log10(sp_ + 1e-9), cmap="Oranges", vmin=-0.7, vmax=0.5, shading="auto", alpha=0.7)   # dark = fast
    ax.contour(X_, Y_, Pc, levels=np.linspace(0.05, np.nanmax(Pc)*0.9, 10), colors=COLORS["ink"], linewidths=0.8)
    ax.contour(X_, Y_, Fc, levels=12, colors=COLORS["teal"], linewidths=0.6, linestyles="--")
    ax.plot([0, 2], [0, 0], "k", lw=2.5)                 # the wall θ = 0
    ax.plot([0, 2.8*np.cos(np.pi/nn)], [0, 2.8*np.sin(np.pi/nn)], "k", lw=2.5)   # the wall θ = π/n
    ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"α = {np.degrees(np.pi/nn):.0f}° (n = {nn:.2g})", fontsize=10)
fig.suptitle("The corner angle decides: calm or violent (orange = log₁₀ speed)", fontweight="bold")
savefig(fig, "ch06", "c09_corners"); plt.show()          # save, then draw
""", see="A right-angle corner, a flat wall, a re-entrant 270° corner and the edge of a plate, each with streamlines "
         "(solid), equipotentials (dashed teal) and the speed as a faint orange shade.",
    read="Crowded streamlines near the corner mean fast flow: none in the 90° corner (a stagnation point, pale), uniform "
         "for the flat wall, and a dark orange spot at the tip of the 270° and 360° corners (the speed grows without bound).",
    change="…α = 60° (n = 3): the fluid near the corner is even calmer, speed ∝ r².")
remind([
    ("np.polyfit for a slope", "`np.polyfit(np.log(r), np.log(q), 1)[0]` is the slope of the log–log line, i.e. the power p in q ∝ rᵖ (Ch. 1 P13)."),
])
nb.figure(r"""
r = np.logspace(-3, 0, 40)                               # distance from the corner along its bisector [m]
fig, ax = plt.subplots(figsize=(6.5, 3.8))               # one panel
for nn, col in zip((2.0, 1.0, 2/3, 0.5), (COLORS["blue"], COLORS["muted"], COLORS["orange"], COLORS["rose"])):
    zb = r*np.exp(1j*np.pi/(2*nn))                       # points on the bisector θ = π/2n
    sp_ = np.abs(pf.Corner(1.0, n=nn).dwdz(zb))          # speed |dw/dz| [m/s]
    slope = np.polyfit(np.log(r), np.log(sp_), 1)[0]     # fitted power
    ax.loglog(r, sp_, color=col, lw=2, label=f"α = {np.degrees(np.pi/nn):.0f}°: slope {slope:+.3f} (n − 1 = {nn - 1:+.3f})")
ax.set_xlabel("distance from the corner r [m]"); ax.set_ylabel("speed |dw/dz| [m/s]")
ax.set_title("Speed near a corner ∝ r^(n−1)"); ax.legend(fontsize=8)
savefig(fig, "ch06", "c09_corner_speed"); plt.show()     # save, then draw
""", see="Four straight lines on log–log axes: rising (90°), flat (180°), falling (270° and 360°).",
    read=r"The slopes are $n-1=(\pi-\alpha)/\alpha$ (D15): the speed vanishes at a convex corner and blows up at a "
         "re-entrant one. The 270° corner is the step of Example 6.2 (C12), where the grid solution will struggle.",
    change="…a rounded corner instead of a sharp one: the lines would bend over to a finite speed below the rounding "
           "radius — real corners are never perfectly sharp.")
nb.plotly(r"""
def corner_traces(nn):                                   # walls and five streamlines of w = z^n
    al = np.pi/nn                                        # corner angle [rad]
    walls = (np.array([2.5, 0, 2.5*np.cos(al)]), np.array([0, 0, 2.5*np.sin(al)]))   # two rays from the origin
    t = np.linspace(1e-3, al - 1e-3, 90)                 # angles inside the wedge [rad]
    xs, ys = [], []                                      # streamlines ψ = r^n sin(nθ) = c, joined with NaN gaps
    for c in (0.1, 0.3, 0.6, 1.0, 1.5):                  # five ψ values [m²/s]
        rr = (c/np.sin(nn*t))**(1/nn)                    # r(θ) on the streamline [m]
        rr = np.where(rr < 3, rr, np.nan)                # clip far pieces
        xs += list(rr*np.cos(t)) + [np.nan]; ys += list(rr*np.sin(t)) + [np.nan]
    return {"walls": walls, "streamlines": (np.array(xs), np.array(ys))}
ns = np.linspace(0.5, 4.0, 22 if not FAST else 11)       # corner exponents n
fig = slider_figure(corner_traces, "n", ns, xlabel="x [m]", ylabel="y [m]", title="One exponent, every corner",
                    xrange=[-2.6, 2.6], yrange=[-2.6, 2.6], height=520)
recolor(fig, {"walls": COLORS["ink"], "streamlines": COLORS["accent"]})
step_titles(fig, [f"One exponent, every corner: α = π/n = {180/v:.0f}°, speed ∝ r^{v - 1:+.2f}" for v in ns])
fig.show()                                               # draw it (works on the web page too)
""")
see_read_change(
    see="Two black walls opening from a narrow 45° wedge (n = 4) to a full plate (n = ½) as you drag n, with five "
        "streamlines filling the wedge.",
    read="The streamlines hug the walls far out and turn smoothly near a convex corner; near a re-entrant corner (n < 1) "
         "they bunch round the tip, where the speed grows like r^(n−1) (the title).",
    change="…n = 1: the walls become one straight wall and the streamlines straight lines — the uniform stream.")
note("N52", "Next: forces", "Complex variables also turn force calculations into one contour integral (§6.5, C10).")
explainer("complex_potential_corners", "Why is the velocity u − iv?",
          "click a point and watch the two difference quotients — along x and along iy — converge to the same dw/dz, "
          "with the velocity drawn as its mirror image; drag n and the walls bend from a 45° wedge to a flat plate while "
          "the speed at the corner flips from zero to infinite.", "", [
              "Click near the corner with n = 2: dw/dz → 0 (stagnation).",
              "Preset 'n = ⅔ (270°)': the speed readout explodes as you click closer to the tip.",
              "Preset 'z* — not a flow': the two quotients disagree and the Cauchy–Riemann bars do not match.",
              "Watch the purple dw/dz arrow and the teal velocity arrow: mirror images in the horizontal.",
          ])
confusion(r"""
plotting dw/dz as the velocity. It is the velocity *mirrored*: u − iv. Use `np.conj(dwdz)` (or
`(dwdz.real, -dwdz.imag)`) for arrows. A wrong sign sends every arrow the wrong way in y.
""")
whatif(r"""
…w were not analytic somewhere inside the region — a source or a vortex? Then w has a singular point there, the flow has
a source of mass or of circulation, and contour integrals round it pick up exactly that strength (C10's residues).
""")

# =====================================================================================================================
# A.5 §6.5 — R18, C10 (N53–N63, D16 D17 D18, E5)
# =====================================================================================================================
nb.section("6.5", "Forces on a Two-Dimensional Body", intro=r"""
**What is this section about?** The pressure round any 2-D body, written with complex numbers, becomes one contour
integral of (dw/dz)² — Blasius's theorem — and that integral can be moved away from the body to a large circle, where
every body looks alike: a uniform stream, its circulation and a doublet. Only the product of the stream and the
circulation survives: no drag, and lift ρUΓ for every shape.
""")
nb.recap("R18", "Momentum balance on a control volume", r"""
For a stationary control volume in steady flow, the momentum carried out through the surface equals the pressure force
on the surface plus the other forces F:
$\int_{A^*}\rho\mathbf u(\mathbf u\cdot\mathbf n)\,dA=-\int_{A^*}p\,\mathbf n\,dA+\mathbf F$ *(Eq. 6.54, from 4.17)*.
`cv_force_on_body` applies it on a large circle round the cylinder of C07 (an exercise idea).
""", where="Ch. 4 §4.4")
nb.code(r"""
print(ch06.cv_force_on_body(pf.cylinder(10.0, 0.1, Gamma_cw=2.0), R_outer=1.0, rho=1.2))   # (D, L) = (0, 24) N/m, air
""", explain="The momentum-flux route on a far circle already gives L = ρUΓ; we now derive the short route.")

core("C10", "Blasius's theorem and Kutta–Zhukhovsky lift for any body", r"""
A round tube and a flat ellipse with the same circulation feel exactly the same lift. Why does the shape not matter?
""", eqs=("6.60", "6.62"))
problem(r"""
C07 found L = ρUΓ for a circle by integrating its surface pressure. Is that a coincidence of the circle? Aircraft wings,
turbine blades and sails are not circles, and integrating their surface pressure by hand is hopeless. We want a proof
that works for *any* cross-section — and it comes from moving the integral off the body.
""")
idea("""
pressure on the body --complex form--> (i rho/2) closed-integral of (dw/dz)^2 dz on the body   (Blasius, (6.60))
      analytic between body and a big circle --> the same integral on the big circle (Cauchy)
far away every body = U + i G/(2 pi z) - d/(2 pi z^2) + ...  --square-->  only U * (i G/pi z) has a 1/z --> L = rho U G, D = 0
""", "(G stands for Γ.)")
P("P156", "Cauchy's integral theorem", r"""
If f(z) is analytic everywhere inside and on a closed curve C, then ∮_C f dz = 0. Consequence: for two curves round the
same singular region, ∮ over one equals ∮ over the other when f is analytic in the ring between them — a contour can be
squeezed or stretched freely, as long as it does not cross a singularity. On a circle z = Re^{iθ}, dz = iz dθ.
""", code=r"""
f_ = lambda z: (1 + 0.3j/z - 0.2/z**2)**2                # analytic except at z = 0
for R in (0.5, 1.0, 4.0):                                # three different circles round 0
    th = np.linspace(0, 2*np.pi, 256, endpoint=False)    # angles [rad]
    z = R*np.exp(1j*th)                                  # points on the circle
    print(R, np.sum(f_(z)*1j*z)*(2*np.pi/256))           # ∮ f dz with dz = iz dθ: the same number each time
""")
nb.md(r"""
The printed number is $2\pi i\times$ (the coefficient of 1/z) $=2\pi i\times0.6i=-1.2\pi\approx-3.770$ on all three
circles — the residue idea of primer P157 below.
""")
note("N53", "The setting", r"""
A body of span B at rest in a steady stream; D (along x) and L (along y) are the forces per unit depth **on the body**;
the body pushes the fluid with F = −B(De_x + Le_y) (Newton's third law). ⚠️ B here is a length, not Ch. 4's Bernoulli
function; "Section 3" in the book means §6.3.
""")
remind([
    ("momentum flux through a surface", r"$\rho\mathbf u(\mathbf u\cdot\mathbf n)\,dA$ is the momentum carried through dA per second (Ch. 4 P114)."),
], lead=r"""**orientation and the outward normal:** walking a closed curve counterclockwise, the body is on your left; the
tangent (dx, dy)/ds turned clockwise by 90°, (dy, −dx)/ds, points out of the body (check at the rightmost point, where
dy > 0 and the normal is +x); parametric curves are Ch. 3 P92""")
D("D16", ref="6.56")
note("N54", "The force components (note `N55`)", r"""
On the body u·n = 0 kills the momentum flux, leaving the pressure force $D\mathbf e_x+L\mathbf e_y=-\frac1B\int_{A^*}p\,\mathbf n\,dA$
*(6.55)*, and with $\mathbf n=(\mathbf e_xdy-\mathbf e_ydx)/ds$ on a counterclockwise contour:
""", equation=EQ["6.56"], ref="6.56")
D("D17", ref="6.60", check_src=r"""
import sympy as sp                                       # symbolic algebra
th = sp.symbols("theta", real=True)                      # angle round the contour (real)
U, a, G, rho = sp.symbols("U a Gamma rho", positive=True)   # stream, radius, clockwise Γ, density (positive)
z = a*sp.exp(sp.I*th)                                    # a point on the body r = a
dwdz = U*(1 - a**2/z**2) + sp.I*G/(2*sp.pi*z)            # (6.52) differentiated: the cylinder with clockwise Γ
u_minus_iv = dwdz                                        # the complex velocity u − iv (6.45)
u_plus_iv = sp.conjugate(dwdz)                           # its conjugate u + iv
dz = sp.diff(z, th)                                      # the step along the circle (times dθ)
print("step 9 on the body:", sp.simplify(sp.expand_complex(u_plus_iv*sp.conjugate(dz) - u_minus_iv*dz)))   # → 0: (6.59)
print("step 6, ∮dz*:", sp.integrate(sp.conjugate(dz), (th, 0, 2*sp.pi)))   # → 0: a uniform pressure gives no force
I_body = sp.integrate(sp.expand(dwdz**2*dz), (th, 0, 2*sp.pi))   # Blasius's integral ∮(dw/dz)² dz on the body
F = sp.simplify(sp.I*rho/2*I_body)                       # D − iL from (6.60)
print("D − iL on r = a :", F)                            # → −I*Gamma*U*rho: D = 0, L = ρUΓ
z2 = 2*a*sp.exp(sp.I*th)                                 # a bigger contour, r = 2a (step 11)
dwdz2 = U*(1 - a**2/z2**2) + sp.I*G/(2*sp.pi*z2)         # the same dw/dz evaluated there
I_big = sp.integrate(sp.expand(dwdz2**2*sp.diff(z2, th)), (th, 0, 2*sp.pi))   # the integral on r = 2a
print("D − iL on r = 2a:", sp.simplify(sp.I*rho/2*I_big))   # the same value: the contour may move (Cauchy)
""")
note("N56", "From pressure to Blasius (notes `N57`–`N60`)", r"""
The complex force $D-iL=-i\oint_Cp\,dz^*$ *(6.57)*; Bernoulli in complex form
$p_\infty+\tfrac12\rho U^2=p+\tfrac12\rho(u-iv)(u+iv)$ gives
$D-iL=-i\oint_C\big[p_\infty+\tfrac12\rho U^2-\tfrac12\rho(u-iv)(u+iv)\big]dz^*$ *(6.58)*; on the body the velocity is
tangent, $(u+iv)\,dz^*=(u-iv)\,dz=(dw/dz)\,dz$ *(6.59)*; so we get Blasius's theorem below. The contour may be any curve
round the body with no singularity of (dw/dz)² in between. The cell checks (6.59) pointwise on the body.
""", equation=EQ["6.60"], ref="6.60")
nb.code(r"""
thb = np.linspace(0, 2*np.pi, 16, endpoint=False)        # 16 angles on the body [rad]
zb = 0.1*np.exp(1j*thb)                                  # body points of the cylinder a = 0.1 m
fc = pf.cylinder(10.0, 0.1, Gamma_cw=2.0)                # C07's spinning cylinder
u, v = fc.velocity(zb.real, zb.imag)                     # the velocity on the body [m/s]
dz = 1j*zb                                               # the direction of dz along the circle
print(np.max(np.abs((u + 1j*v)*np.conj(dz) - (u - 1j*v)*dz)))   # ≈ 0: (6.59) holds point by point on the body
""")
remind([
    ("FFT of samples on a circle", "the FFT turns N equally spaced samples into the coefficients of e^{ikθ} (Ch. 5 P142); on a circle z = Re^{iθ} these are the c_k R^k of a power series."),
])
P("P157", "Laurent series and residues", r"""
Outside a disc that contains all the singular points, an analytic function is a sum of powers of z including negative
ones: $f(z)=\sum_k c_kz^k$ (a Laurent series). Integrating term by term round a circle, every power gives zero
**except** $z^{-1}$, which gives 2πi — so ∮f dz = 2πi c₋₁, and c₋₁ is called the **residue**. Sampling f on a circle
and taking an FFT reads off all the c_k at once.
""", code=r"""
R, n = 2.0, 64                                           # a circle and the number of samples
z = R*np.exp(2j*np.pi*np.arange(n)/n)                    # n points on it
f_ = 3 + 0.5/z - 2/z**2                                  # c0 = 3, c−1 = 0.5, c−2 = −2
ck = np.fft.fft(f_)/n                                    # FFT coefficient of e^{ikθ} …
print(ck[0], ck[-1]*R, ck[-2]*R**2)                      # … rescaled by R^k: 3, 0.5, −2
""")
note("N63", "The residue theorem on a circle", r"""
With z = Re^{iθ}, dz = iRe^{iθ}dθ (no book number). Only n = 1 survives: the integral round any closed curve is 2πi times
the sum of the residues inside.
""", equation=r"\oint z^{-n}dz=iR^{1-n}\int_0^{2\pi}e^{i(1-n)\theta}d\theta=2\pi i\,\delta_{n1}")
nb.code(r"""
th = np.linspace(0, 2*np.pi, 64, endpoint=False)         # 64 angles [rad]
z = 1.5*np.exp(1j*th)                                    # a circle of radius 1.5
for n in (0, 1, 2, 3):                                   # powers z^(−n)
    print(n, np.round(np.sum(z**(-n)*1j*z)*(2*np.pi/64), 12))   # ∮ z^(−n) dz: 0, 2πi, 0, 0
""")
note("N61", "Every body from far away", r"""
Outside a circle round the body, $w=Uz+\frac{m}{2\pi}\ln z+\frac{i\Gamma}{2\pi}\ln z+\frac{d}{2\pi z}+\dots$ (stream, net
source, clockwise vortex, doublet); a closed body has m = 0. The cell reads the coefficients of dw/dz for C07's
cylinder by FFT.
""")
nb.code(r"""
c = pf.laurent_coefficients(pf.cylinder(10.0, 0.1, Gamma_cw=2.0).dwdz, R=0.3)   # c_k of dw/dz on |z| = 0.3 m
print(np.round(c[0], 10), np.round(c[-1], 10), 1j*2.0/(2*np.pi), np.round(c[-2], 10))   # U, iΓ/2π (twice), −Ua²
""", explain=r"""
$c_0=U=10$ m/s (the stream), $c_{-1}=i\Gamma/2\pi=0.3183i$ m²/s (the clockwise vortex — matched by the third number),
$c_{-2}=-Ua^2=-0.1$ m³/s $=-d/2\pi$ with $d=2\pi Ua^2$ — the cylinder's doublet.
""")
P("P158", "sympy for complex series and residues", r"""
sympy handles the algebra of D18: `sp.I` is i, `sp.expand` multiplies out a squared series, `.coeff(z, -1)` reads a
coefficient, and `sp.residue(f, z, 0)` returns the coefficient of 1/z directly.
""", code=r"""
z, U, G, d = sp.symbols("z U Gamma d")                   # symbols
f_ = sp.expand((U + sp.I*G/(2*sp.pi*z) - d/(2*sp.pi*z**2))**2)   # the squared far field
print(f_.coeff(z, -1), "|", f_.coeff(z, -2))             # I*Gamma*U/pi | −Gamma²/(4π²) − U d/π
print(sp.residue(f_, z, 0))                              # I*Gamma*U/pi
""")
D("D18", ref="6.62", check_src=r"""
import sympy as sp                                       # symbolic algebra
z = sp.symbols("z")                                      # the complex position
U, G, d, rho = sp.symbols("U Gamma d rho", real=True)    # stream, clockwise Γ, doublet strength, density (real)
f = U + sp.I*G/(2*sp.pi*z) - d/(2*sp.pi*z**2)            # step 5: the far-field dw/dz
sq = sp.expand(f**2)                                     # step 6: its square
print("1/z coefficient:", sq.coeff(z, -1))               # → I*Gamma*U/pi
true_c2 = sq.coeff(z, -2)                                # the true 1/z² coefficient
printed_c2 = U*d/sp.pi - G**2/(4*sp.pi**2)               # the coefficient the book prints
print(sp.simplify(true_c2 - (-(U*d/sp.pi + G**2/(4*sp.pi**2)))),   # → 0: the correct form −(Ud/π + Γ²/4π²)
      "|", sp.simplify(true_c2 - printed_c2))            # → −2Ud/π ≠ 0: the printed coefficient is wrong
res = sp.residue(sq, z, 0)                               # step 8: the residue
F = sp.I*rho/2*2*sp.pi*sp.I*res                          # steps 9–10: D − iL = (iρ/2) 2πi × residue
print("D =", sp.simplify(sp.re(F)), "| L =", sp.simplify(-sp.im(F)))   # → 0 and Gamma*U*rho
print(ch06.kutta_zhukhovsky_sym()["L"])                  # the library's sympy route: the same Gamma*U*rho
""")
note("N62", "Blasius with the far field", r"""
⚠️ The book prints the squared series with the 1/z² coefficient $(Ud/\pi-\Gamma^2/4\pi^2)$ and an extra outer square;
expanding gives $-(Ud/\pi+\Gamma^2/4\pi^2)$ (D18 step 6 and its sympy cell). Nothing depends on it: only the 1/z term,
$iU\Gamma/\pi$, has a residue.
""", equation=EQ["6.61"], ref="6.61")
nb.worked_example("the same spinning cylinder, by residues", r"""
U = 10 m/s, Γ = 2 m²/s (clockwise), ρ = 1.2 kg/m³.

1. Residue of (dw/dz)² at 0: $iU\Gamma/\pi=i\times20/\pi=6.366i$ m³/s².
2. $\oint(dw/dz)^2dz=2\pi i\times6.366i=-40$ m⁴/s² (= −2UΓ).
3. $D-iL=\frac{i\rho}{2}\times(-40)=-24i$ N/m.
4. D = 0, L = 24 N/m — the same as C07's pressure integral.
""")
nb.code(r"""
fc = pf.cylinder(10.0, 0.1, Gamma_cw=2.0)                # C07's spinning cylinder (a = 0.1 m)
for R in (0.1, 0.2, 1.0):                                # three contour radii [m]
    print("cylinder, R =", R, pf.blasius_force(fc.dwdz, R=R, rho=1.2))   # (D, L) = (0, 24) N/m on every circle
e = ch06.elliptic_cylinder_flow(10.0, 0.12, 0.1, Gamma_cw=2.0)   # a Zhukhovsky ellipse (C11), same U and Γ
for R in (0.5, 1.0):                                     # circles that enclose the ellipse (semi-major axis 0.2033 m)
    print("ellipse,  R =", R, pf.blasius_force(e.dwdz, R=R, rho=1.2))   # the same lift: the shape does not matter
kz = ch06.kutta_zhukhovsky_sym()                         # the sympy route of D18
print(kz["residue"], "|", kz["coeff_z2"], "vs printed", kz["coeff_z2_printed"], "|", kz["L"])
print(ch06.contour_force(lambda zz: fc.pressure(zz.real, zz.imag, rho=1.2, U=10.0),
                         0.1*np.exp(1j*np.linspace(0, 2*np.pi, 128, endpoint=False))))   # (6.56) on the body: (0, 24)
""", explain=r"""
1. Blasius (6.60) on three circles round the cylinder: the same force (Cauchy lets the contour move).
2. A different body — the Zhukhovsky ellipse of C11 with the same U and Γ — gives the same lift on two contours.
3. sympy: the residue $iU\Gamma/\pi$ and the true 1/z² coefficient next to the printed one (the book slip of N62).
4. The surface-pressure route (6.56) on 128 body points agrees (to the accuracy of its 128-segment polygon).
""")
nb.md("**From scratch — Blasius's integral as our own periodic trapezoid sum:**")
nb.check_agree(r"""
R, nq = 0.3, 128                                         # contour radius [m] and number of points
th = np.linspace(0, 2*np.pi, nq, endpoint=False)         # angles [rad]
z = R*np.exp(1j*th)                                      # contour points
integral = np.sum(fc.dwdz(z)**2*1j*z)*(2*np.pi/nq)       # ∮(dw/dz)² dz with dz = iz dθ
F = 1j*1.2/2*integral                                    # D − iL from (6.60), air
print(F)                                                 # ≈ −24j
assert np.allclose([F.real, -F.imag], pf.blasius_force(fc.dwdz, R=R, rho=1.2))   # = the library's (D, L)
""")
nb.md("The periodic trapezoid round a circle is Blasius's integral (Ch. 5 P136). For the cylinder it is even exact "
      "with 4 or more points, because its (dw/dz)² is a finite sum of powers of z down to z⁻⁴; for a body like the "
      "Zhukhovsky ellipse the series is infinite and the error falls exponentially with the number of points.")
nb.figure(r"""
Ra = np.linspace(1.0, 20.0, 40)                          # contour radius / cylinder radius
cyl = np.array([pf.blasius_force(fc.dwdz, R=0.1*k, rho=1.2) for k in Ra])   # (D, L) for the cylinder [N/m]
Re_ = Ra[Ra > 2.1]                                       # contours that enclose the ellipse (semi-major 0.2033 m)
ell = np.array([pf.blasius_force(e.dwdz, R=0.1*k, rho=1.2) for k in Re_])   # (D, L) for the ellipse [N/m]
Dp, Lp = ch06.contour_force(lambda zz: fc.pressure(zz.real, zz.imag, rho=1.2, U=10.0), 0.1*np.exp(1j*np.linspace(0, 2*np.pi, 256, endpoint=False)))   # (6.56) on the body
lc = ch06.laurent_contributions("cylinder", R=0.2, Gamma_cw=2.0, U=10.0, a=0.1, rho=1.2)   # share of each power z^k
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.7))  # force vs contour | Laurent bars
a1.plot(Ra, cyl[:, 1], color=COLORS["amber"], lw=2.5, label="L, cylinder")   # lift: amber
a1.plot(Ra, cyl[:, 0], color=COLORS["muted"], lw=2.5, label="D, cylinder")   # drag: grey
a1.plot(Re_, ell[:, 1], color=COLORS["amber"], ls="--", lw=1.5, label="L, Zhukhovsky ellipse")   # another body, same Γ
a1.plot(Re_, ell[:, 0], color=COLORS["muted"], ls="--", lw=1.5, label="D, ellipse")   # its drag
a1.plot([1.0, 1.0], [Dp, Lp], "D", color=COLORS["ink"], label="surface-pressure route (6.56)")   # the independent route at R = a
a1.set_xlabel("contour radius R / a [–]"); a1.set_ylabel("force per unit depth [N/m]"); a1.legend(fontsize=7)   # axes with units
a1.set_ylim(-3, 28)   # fixed range: 0 and ρUΓ = 24 N/m in view
a2.bar([str(k) for k in lc["powers"]], lc["contrib_im"], color=[COLORS["amber"] if k == -1 else COLORS["muted"] for k in lc["powers"]])   # one bar per power
a2.set_xlabel("power k of z in (dw/dz)$^2$"); a2.set_ylabel("contribution to L [N/m]")   # axes with units
a2.set_title("only the 1/z term survives the integral", fontsize=10)   # the message of the bars
fig.suptitle("Any contour, any body: only U × Γ/z survives", fontweight="bold")   # the message
savefig(fig, "ch06", "c10_blasius"); plt.show()          # save, then draw
""", see="Left, flat lines: L = 24 N/m and D = 0 for every contour radius and for both bodies, with the surface-pressure "
         "diamond on top at R = a; right, one tall bar at k = −1 and nothing elsewhere.",
    read=r"Cauchy lets the contour move; the residue picks the 1/z term; its coefficient is $iU\Gamma/\pi$ whatever the "
         "body's doublet d — so the lift is ρUΓ for any shape.",
    change="…a source hidden in the body (an open body, m ≠ 0): the 1/z coefficient gains Um/π and the drag becomes "
           "D = −ρUm — a thrust, pushing upstream (the held-source force of the S03 exercises).")
explainer("blasius_kutta_contour", "Why doesn't lift depend on the shape?",
          "a proof cannot show invariance under deformation; here you drag the contour — radius, offset, shape — and "
          "watch the integral stay pinned while the integrand changes completely; the Laurent bars show which single "
          "term survives.", "", [
              "Drag the contour outward: the lift readout stays at ρUΓ.",
              "Switch the body to the ellipse with the same Γ: the same lift.",
              "Preset 'contour cuts the body' with the ellipse: the number jumps (R = 0.07 m gives about 8 N/m) because the "
              "contour crosses the map's branch cut. With the cylinder it does not: all its singularities sit at the centre, "
              "which a cutting circle still encloses — what matters is the singularities, not the outline.",
              "Tilt the ellipse: the force stays ρUΓ in size and perpendicular to the stream (its x and y parts change).",
              "Preset 'few quadrature points': exact for the cylinder already with 4 points; for the ellipse the error "
              "falls exponentially as you add points.",
          ])
whatif(r"""
…the flow were three-dimensional (a finite wing)? The 2-D proof no longer applies: the trailing vortices of a finite wing
induce a downwash and a drag even in ideal flow — induced drag, Ch. 14. In 2-D there is no such escape: D = 0 always.
""")

# =====================================================================================================================
# A.6 §6.6 — C11 (N64–N70, D19 D20 D21, E6)
# =====================================================================================================================
nb.section("6.6", "Conformal Mapping", intro=r"""
**What is this section about?** An analytic function maps one plane to another, turning and stretching every small
figure but keeping its angles. A flow net — streamlines crossing equipotentials at right angles — therefore maps to
another flow net. Solve the easy problem, flow round a circle, and let a map carry it: the Zhukhovsky map turns the circle
into an ellipse (and, in Ch. 14, an airfoil).
""")
core("C11", "Conformal maps keep angles, and the Zhukhovsky map carries the circle's flow onto an ellipse", r"""
We can solve the flow round a circle. How can that one solution give the flow round other shapes?
""", eqs=("6.63", "6.65"))
problem(r"""
Draw a circle on a rubber sheet with a fine square grid, and stretch the sheet smoothly without tearing: the circle
becomes an oval, but where you look closely the little squares are still little squares, just turned and resized. A
conformal map is such a stretch done by an analytic function. Because streamlines and equipotentials form a grid of
little squares, the map carries a flow into a flow.
""")
idea("""
zeta-plane: circle |zeta| = a, flow (6.68) known  --z = zeta + b^2/zeta-->  z-plane: ellipse (6.67), flow = the same w at zeta(z)
small cross at zeta0  --x f'(zeta0) = |f'| e^{i chi}-->  turned by chi, stretched by |f'|, right angle kept
inverse: two roots zeta, b^2/zeta  --  keep the one OUTSIDE the circle
""", "🔁 Reminder (primer P153 in C09): multiplying complex numbers multiplies their moduli and adds their angles.")
note("N64", "A small step δz at z₀ becomes", "", equation=EQ["6.63"], ref="6.63")
D("D19", ref="6.64")
note("N65", "Two elements at one point", r"""
A second step δ′z at the same point becomes $\delta'w=\frac{dw}{dz}\delta'z$ *(Eq. 6.64)*, turned by the same angle
arg(dw/dz), so the angle between the two is kept (α = β). It fails where dw/dz = 0 or ∞, and only *small* figures keep
their shape.
""")
nb.code(r"""
f_, fp = (lambda z: z**2), (lambda z: 2*z)               # the map w = z² and its derivative
print(np.degrees(cm.angle_preservation(f_, fp, 1+0.5j, 1e-3, 1e-3j)))   # (before, after) = (90°, 90°): kept
print(np.degrees(cm.angle_preservation(f_, fp, 0j, 1e-3, 1e-3*np.exp(1j*np.pi/4))))   # (45°, 90°): doubled at z = 0
""")
note("N66", "Flow nets as images of a grid", r"""
A rectangular grid of φ = const and ψ = const lines in the w-plane is uniform flow; its image in the z-plane is the flow
net of w = f(z). Chains of maps work too: w = ln ζ with ζ = sin z gives
$u-iv=\frac{dw}{d\zeta}\frac{d\zeta}{dz}=\frac1\zeta\cos z=\cot z$; w = z² gives the 90° corner (the hyperbolae of
$\psi=2Axy$ *(6.24)* — ⚠️ the book's "(see Figure 6.5)" means Figure 6.3).
""")
nb.code(r"""
r_ = ch06.cot_flow(0.7 + 0.2j)                           # w = ln(sin z) at z = 0.7 + 0.2i
print(r_["dwdz_chain"], 1/np.tan(0.7 + 0.2j))            # (dw/dζ)(dζ/dz) = cot z: equal
""")
D("D20", ref="6.67")
note("N67", "The circle of radius a > b maps to an ellipse (note `N68`)", r"""
$z=ae^{i\theta}+\frac{b^2}{a}e^{-i\theta}$ *(6.66)* is the ellipse below, with foci at ±2b. b = 1 m, a = 1.2 m:
semi-axes 2.033 m and 0.367 m, foci ±2 m.
""", equation=EQ["6.67"], ref="6.67")
nb.code(r"""
print(ch06.joukowski_ellipse(1.2, 1.0))                  # semi-axes A = a + b²/a, B = a − b²/a and the focal distance 2b [m]
""")
note("N69", "The circle flow in the ζ-plane", r"""
The circle of radius a in a stream U with a clockwise circulation Γ — $w=U\big(z+\frac{a^2}z\big)+\frac{i\Gamma}{2\pi}\ln(z/a)$
*(Eq. 6.52)* with z replaced by ζ:
""", equation=EQ["6.68"], ref="6.68")
P("P159", "complex square roots and the quadratic formula", r"""
A quadratic $\zeta^2-z\zeta+b^2=0$ has two roots $\zeta=\tfrac12[z\pm(z^2-4b^2)^{1/2}]$, and in the complex plane "the"
square root is a choice: numpy's `np.sqrt` returns the principal root (real part ≥ 0), whose cut lies where the argument
is a negative real number. The two roots here multiply to b² (Vieta, P71), so one is outside the circle |ζ| = b and one
inside — and the principal root does **not** always pick the outside one.
""", code=r"""
b, z = 1.0, -3 + 0.5j                                    # a point to the left of the slit
r = np.sqrt(z**2 - 4*b**2)                               # numpy's principal root
print(abs(0.5*(z + r)), abs(0.5*(z - r)))                # 0.370 and 2.701: '+' gives the INSIDE root here
print(abs(0.5*(z + np.sqrt(z - 2*b)*np.sqrt(z + 2*b))))  # 2.701: the product form picks the outside root
""")
D("D21", ref="6.69", check_src=r"""
import sympy as sp                                       # symbolic algebra
zeta, zs, bs_ = sp.symbols("zeta z b")                   # the unknown ζ, the given z and the constant b
roots = sp.solve(zeta**2 - zs*zeta + bs_**2, zeta)       # step 2: the two roots of ζ² − zζ + b² = 0
print(roots)                                             # ½[z ± √(z² − 4b²)]
print("product of the roots:", sp.simplify(roots[0]*roots[1]))   # step 3: → b**2 (Vieta)
xg = np.linspace(-4, -0.5, 41)                           # numeric part: a grid in the left half-plane …
yg = np.linspace(-2, 2, 41)                              # … x < 0, all y
zz = (xg[:, None] + 1j*yg[None, :]).ravel()              # 1681 points z = x + iy
zz = zz[~((np.abs(zz.imag) < 1e-12) & (np.abs(zz.real) <= 2.0))]   # drop the slit [−2b, 2b] itself: there both roots have |ζ| = b
zz = zz[np.abs(cm.joukowski_inverse(zz, 1.0)) > 1.0]    # keep points outside the circle's image (all of them here)
bad = np.abs(0.5*(zz + np.sqrt(zz**2 - 4))) < 1          # the principal-root variant lands INSIDE the circle
print("fraction where the principal root fails:", bad.mean())   # → 1.0 at every left-half point off the slit
good = cm.joukowski_inverse(zz, 1.0)                     # the library's outside branch (product form, step 6)
assert np.allclose(cm.joukowski(good, 1.0), zz)          # it really inverts the map z = ζ + b²/ζ …
assert np.all(np.abs(good) >= 1)                         # … and always lands outside the circle |ζ| = b
U_, a_, G_ = sp.symbols("U a Gamma", positive=True)      # stream, circle radius, clockwise Γ
W = U_*(zeta + a_**2/zeta) + sp.I*G_/(2*sp.pi)*sp.log(zeta/a_)   # (6.68): the circle flow
zfun = zeta + bs_**2/zeta                                # (6.65): z as a function of ζ
dWdz_chain = sp.diff(W, zeta)/sp.diff(zfun, zeta)        # steps 9–10: dW/dz = (dW/dζ)/(dz/dζ) (implicit differentiation)
step11 = (U_*(1 - a_**2/zeta**2) + sp.I*G_/(2*sp.pi*zeta))/(1 - bs_**2/zeta**2)   # the velocity u − iv of step 11
print("step 11 − chain rule:", sp.simplify(dWdz_chain - step11))   # → 0
""")
note("N70", "The inverse and the velocity", r"""
with the root that falls outside the circle, and $u-iv=\frac{dw}{dz}=\frac{dw}{d\zeta}\frac{d\zeta}{dz}$. ⚠️ Coded as
$\zeta=\tfrac12[z+\sqrt{z-2b}\sqrt{z+2b}]$: at z = −3 + 0.5i (b = 1) it gives |ζ| = 2.70, while numpy's principal
√(z² − 4b²) gives 0.37 — a point inside the cylinder, and a flow that makes no sense.
""", equation=EQ["6.69"], ref="6.69")
nb.worked_example("an ellipse from a circle", r"""
b = 1 m, a = 1.2 m.

1. The circle $\zeta=1.2e^{i\theta}$ maps to $x=(1.2+1/1.2)\cos\theta=2.033\cos\theta$, $y=(1.2-0.833)\sin\theta=0.367\sin\theta$.
2. Foci: c² = 2.033² − 0.367² = 4.134 − 0.134 = 4 → c = 2 = 2b.
3. At θ = 0 the map's derivative dz/dζ = 1 − b²/ζ² = 1 − 1/1.44 = 0.306: lengths along the circle near its right end
   shrink to 31 % — the ellipse's sharp ends.
4. Inverse at z = −3 + 0.5i: the product form gives ζ = −2.638 + 0.579i (|ζ| = 2.70 > 1.2, outside ✓).
""")
nb.code(r"""
b, a = 1.0, 1.2                                          # Zhukhovsky constant and circle radius [m]
th = np.linspace(0, 2*np.pi, 400, endpoint=False)        # angles on the circle [rad]
zc = cm.joukowski(a*np.exp(1j*th), b)                    # the circle's image z = ζ + b²/ζ (6.65)
A_, B_ = a + b**2/a, a - b**2/a                          # semi-axes [m]
print(np.max(np.abs(zc.real**2/A_**2 + zc.imag**2/B_**2 - 1)))   # ≈ 0: the image satisfies (6.67)
zz = np.array([-3+0.5j, -3-0.5j, 3+0.5j, 0.1+2j])         # four points, left and right of the slit
print(np.abs(cm.joukowski_inverse(zz, b)))               # outside branch: every |ζ| > 1
print(np.abs(cm.joukowski_inverse(zz, b, branch="principal")))   # numpy's principal root: 0.37 on the left — WRONG
e = ch06.elliptic_cylinder_flow(1.0, a, b, Gamma_cw=1.0) # the circle flow (6.68) carried to the ellipse, U = 1 m/s, Γ = 1 m²/s
print(pf.normal_velocity_on(e, zc), pf.far_field_check(e, 200.0))   # no flow through the ellipse; far field → U (error ~Γ/2πR)
""", explain=r"""
1. The circle's image satisfies the ellipse equation (6.67) to round-off.
2. The two inverses compared in all four quadrants: the outside branch always lands outside the circle; the principal
   root sends the two left points inside (|ζ| = 0.37).
3. The mapped flow (`elliptic_cylinder_flow`, the chain rule of D21) has no flow through the ellipse and returns to the
   stream far away (the remaining 8×10⁻⁴ m/s at R = 200 m is the vortex's Γ/2πR).
""")
nb.md("**From scratch — pick the root outside the circle by hand:**")
nb.check_agree(r"""
r1 = 0.5*(zz + np.sqrt(zz**2 - 4*b**2))                  # numpy's principal '+' root
r2 = b**2/r1                                             # the other root: the two multiply to b² (Vieta)
mine = np.where(np.abs(r1) >= b, r1, r2)                 # keep the root outside the circle |ζ| = b
print(np.abs(mine))                                      # all > 1
assert np.allclose(mine, cm.joukowski_inverse(zz, b))    # = the library's branch everywhere
""")
nb.md("Picking the root with |ζ| ≥ b by hand gives the library's branch everywhere.")
nb.figure(r"""
n = 241 if not FAST else 161                             # grid points per side
W, dW = cm.circle_flow_zeta(1.0, 1.2, Gamma_cw=1.0)      # the ζ-plane flow (6.68): U = 1 m/s, a = 1.2 m, Γ = 1 m²/s
Xs, Ys = np.meshgrid(np.linspace(-3, 3, n), np.linspace(-2.2, 2.2, n))   # grid [m]
Zs = Xs + 1j*Ys                                          # as complex numbers
Pz = np.where(np.abs(Zs) < 1.2, np.nan, W(Zs).imag)      # ψ in the ζ-plane, blank inside the circle
Pe = np.asarray(e.psi(Xs, Ys), float)                    # ψ of the mapped flow in the z-plane
Pe = np.where(Xs**2/A_**2 + Ys**2/B_**2 < 1, np.nan, Pe) # blank inside the ellipse
lv = np.linspace(-2.5, 2.5, 21)                          # the SAME ψ levels in both planes [m²/s]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.0))    # ζ-plane | z-plane
a1.contour(Xs, Ys, Pz, levels=lv, colors=COLORS["ink"], linewidths=0.8)   # circle-flow streamlines
draw_body(a1, circle(1.2)); a1.plot([-1, 1], [0, 0], "x", color=COLORS["rose"], ms=9, mew=2.5, zorder=10)   # the circle and the points ±b
a1.set_title("ζ-plane: circle a = 1.2 m (× = ±b)", fontsize=10)
a2.contour(Xs, Ys, Pe, levels=lv, colors=COLORS["ink"], linewidths=0.8)   # the same streamlines, mapped
draw_body(a2, zc); a2.plot([-2, 2], [0, 0], "--", color=COLORS["muted"])   # the ellipse and the slit [−2b, 2b] inside it
a2.set_title("z-plane: the ellipse z = ζ + b²/ζ", fontsize=10)
for ax in (a1, a2):
    ax.set_xlim(-3, 3); ax.set_ylim(-2.2, 2.2); ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
fig.suptitle("Map the circle, carry the flow", fontweight="bold")
savefig(fig, "ch06", "c11_two_planes"); plt.show()       # save, then draw
""", see="The same streamline pattern twice: bent round a circle on the left and round a flat ellipse on the right, "
         "slightly asymmetric top/bottom because of the clockwise circulation.",
    read="ψ is the same number at a point ζ and at its image z, so a streamline maps to a streamline and the circle maps "
         "to the ellipse; angles between streamlines and the body are kept.",
    change="…a → b: the ellipse flattens to the plate [−2, 2]; the speed at its ends blows up unless Γ is chosen to cancel "
           "it (the Kutta condition, Ch. 14).")
nb.figure(r"""
ng = 201 if not FAST else 141                            # grid points per side
maps = [("w = z²", lambda z: z**2, (-2, 2), (-2, 2)), ("w = ln z (a source)", np.log, (-2, 2), (-2, 2)),
        ("w = ln(sin z)", lambda z: np.log(np.sin(z)), (0.05, np.pi - 0.05), (-1.5, 1.5))]   # three flows
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))         # one panel per map
for ax, (name, wf, xl, yl) in zip(axs, maps):
    g = cm.grid_image(wf, xlim=xl, ylim=yl, n=ng, levels=18)   # Re w and Im w on a grid, cuts masked
    if name.startswith("w = ln z"):                      # ln z: blank a strip along its branch cut (the negative x-axis)
        cut = (np.abs(g["Y"]) < 0.05) & (g["X"] < 0)
        g["psi"] = np.where(cut, np.nan, g["psi"]); g["phi"] = np.where(cut, np.nan, g["phi"])
    ax.contour(g["X"], g["Y"], g["psi"], levels=g["psi_levels"], colors=COLORS["ink"], linewidths=0.7)   # ψ = const
    ax.contour(g["X"], g["Y"], g["phi"], levels=g["phi_levels"], colors=COLORS["teal"], linewidths=0.6, linestyles="--")   # φ = const
    ax.set_aspect("equal"); ax.set_title(name, fontsize=10); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
fig.suptitle("The image of a square grid is a flow net", fontweight="bold")
savefig(fig, "ch06", "c11_grid_images"); plt.show()      # save, then draw
""", see="Hyperbolae crossing hyperbolae (z²), rays crossing circles (ln z: the flow net of a source), and the cot-flow net "
         "in a strip (ln sin z).",
    read="Every crossing is a right angle except at critical points where dw/dz = 0 — the origin for z² — where angles "
         "are doubled.",
    change="…w = z³: the pattern repeats three times round the origin, and the angles there are tripled.")
nb.plotly(r"""
bq = 1.0                                                 # Zhukhovsky constant [m]
def joukowski_traces(ab):                                # body and five mapped streamlines for a/b = ab (Γ = 0, U = 1 m/s)
    aq = ab*bq                                           # circle radius [m]
    th = np.linspace(0, 2*np.pi, 120)                    # angles [rad]
    body = cm.joukowski(aq*np.exp(1j*th), bq)            # the ellipse (a plate when a = b)
    t = np.linspace(0.005, np.pi - 0.005, 90)            # angles in the upper half [rad]
    xs, ys = [], []                                      # streamlines, joined with NaN gaps
    for c in (0.25, 0.6, 1.0, 1.6, 2.4):                 # ψ values [m²/s]
        rr = (c + np.sqrt(c**2 + 4*aq**2*np.sin(t)**2))/(2*np.sin(t))   # r(θ) on ψ = (r − a²/r) sinθ = c
        zm = cm.joukowski(rr*np.exp(1j*t), bq)           # carried to the z-plane
        xs += list(zm.real) + [np.nan] + list(zm.real) + [np.nan]; ys += list(zm.imag) + [np.nan] + list(-zm.imag) + [np.nan]
    return {"body": (body.real, body.imag), "streamlines": (np.array(xs), np.array(ys))}
abs_ = np.linspace(1.0, 3.0, 21 if not FAST else 11)     # a/b values
fig = slider_figure(joukowski_traces, "a/b", abs_, xlabel="x [m]", ylabel="y [m]",
                    title="One map, a family of bodies", xrange=[-5, 5], yrange=[-3.2, 3.2], height=460)
recolor(fig, {"body": COLORS["ink"], "streamlines": COLORS["accent"]})
step_titles(fig, [f"One map, a family of bodies: semi-axes {v + 1/v:.2f} m × {v - 1/v:.2f} m, foci ±2 m" for v in abs_])
fig.show()                                               # draw it (works on the web page too)
""")
see_read_change(
    see="At a/b = 1 a flat plate from −2 to 2 m with streamlines sliding past it; dragging a/b up fattens it into an "
        "ellipse and then almost a circle.",
    read="Every body of the family shares the foci ±2b (the plate's ends); the streamlines are the circle's, carried "
         "through the map.",
    change="…a/b → ∞: z ≈ ζ everywhere near the body (b²/ζ becomes negligible) and the body is simply the circle.")
explainer("conformal_joukowski", "How does a map carry a flow?",
          "two linked planes on one state is the natural interactive form: drag a small cross in the ζ-plane and watch "
          "its image stay a right angle (and fail at ζ = ±b), slide a/b from a plate to a circle, and toggle numpy's "
          "square root to see half the streamlines jump inside the body.", "", [
              "Drag the cross toward ζ = +b: the image angle opens toward 180°.",
              "Slide a/b to 1: the ellipse becomes a plate.",
              "Switch the branch chips to 'numpy principal': streamlines on the left jump inside — the status turns ⚠️.",
              "Mode 'z²' at 0: angles doubled.",
          ])
whatif(r"""
…the circle were shifted off-centre (up and a little left) before mapping? The image would be a cambered, round-nosed,
sharp-tailed shape — a Zhukhovsky airfoil. Ch. 14 does exactly this, and fixes Γ so the flow leaves the sharp tail
smoothly.
""")

# =====================================================================================================================
# A.7 §6.7 — R19 R20 R21, C12 (N71–N74, D22 D23, E7), pointer N90
# =====================================================================================================================
nb.section("6.7", "Numerical Solution Techniques in Two Dimensions", intro=r"""
**What is this section about?** Exact solutions exist only for simple shapes. On a grid, Laplace's equation becomes a
rule — every value is the average of its four neighbours — and a sweep that keeps enforcing the rule converges to the
solution. This is the book's first numerical PDE solver, and the ancestor of every elliptic solver in Ch. 10 and of
bounded-domain stream-function inversion in ocean and atmosphere models.
""")
nb.recap("R19", "Grids and half-point differences", r"""
Store ψ at grid points, $\psi_{i,j}=\psi(i\Delta x,j\Delta y)$, and approximate a first derivative by the difference
across a cell centred on the point, $(\partial\psi/\partial x)_{i,j}\simeq(\psi_{i+\frac12,j}-\psi_{i-\frac12,j})/\Delta x$
(Ch. 2 §2.9 did central differences). ⚠️ The book calls these "first-order" differences: they approximate a *first
derivative* and are **second-order accurate** (error ∝ Δx², D22 step 3).
""", where="Ch. 2 §2.9")
nb.recap("R20", "Second difference in x", r"""
Differencing twice: $\big(\frac{\partial^2\psi}{\partial x^2}\big)_{i,j}\simeq\frac{\psi_{i+1,j}-2\psi_{i,j}+\psi_{i-1,j}}{\Delta x^2}$ *(Eq. 6.70)*.
""", where="Ch. 2 §2.9")
nb.recap("R21", "Second difference in y", r"""
$\big(\frac{\partial^2\psi}{\partial y^2}\big)_{i,j}\simeq\frac{\psi_{i,j+1}-2\psi_{i,j}+\psi_{i,j-1}}{\Delta y^2}$ *(Eq. 6.71)*.
Below, `ls.laplacian_5pt` adds the two on three grids for the harmonic function sin(πx) sinh(πy), whose exact
Laplacian is 0 — so whatever the stencil returns is pure error.
""", where="Ch. 2 §2.9")
nb.code(r"""
hs = np.array([0.025, 0.0125, 0.00625])                 # three grid spacings [m]
err = []                                                 # max |discrete ∇²ψ| on each grid
for h in hs:
    xg = np.arange(0.0, 1.0 + h/2, h)                    # grid points on [0, 1] [m]
    Xg, Yg = np.meshgrid(xg, xg)                         # all nodes (rows = y, columns = x)
    psi = np.sin(np.pi*Xg)*np.sinh(np.pi*Yg)             # a harmonic function: its exact ∇² is 0
    lap = ls.laplacian_5pt(psi, h, h)                    # (6.70) + (6.71) at interior nodes
    err.append(np.max(np.abs(lap[1:-1, 1:-1])))          # the stencil's error
print(np.array(err), "pairwise orders", np.log2(np.array(err[:-1])/np.array(err[1:])))   # ÷≈4 per halving: order → 2
""")

core("C12", "On a grid, Laplace says \"be the average of your neighbours\" — and Gauss–Seidel relaxes to it", r"""
How does a grid of numbers solve Laplace's equation — and how do we know when it has?
""", eqs=("6.72",))
problem(r"""
Water flows through a channel that narrows abruptly. No formula gives its streamlines — but ψ must be constant on each
wall and harmonic in between. Stretch a rubber membrane over a wire frame bent to the boundary values: every point
settles at the average height of its neighbours. A computer can do the settling one point at a time.
""")
idea("""
          psi(i, j+1)
psi(i-1, j)  psi(i, j)  psi(i+1, j)       psi(i, j) = 1/4 (sum of the four)          <-- (6.72)
          psi(i, j-1)                     sweep: replace every value by the average of its neighbours, repeat
Jacobi: use last sweep's values  .  Gauss-Seidel: use the newest  .  SOR: overshoot by omega
""")
remind([
    ("np.linalg.solve", "`np.linalg.solve(A, b)` solves the linear system Aψ = b directly (Ch. 1 P57)."),
    ("spectral radius", "the largest |eigenvalue| of a matrix G; repeated multiplication by G shrinks a vector by about that factor each time (eigenvalues: Ch. 2 P80)."),
], lead="Taylor series is Ch. 1 P26; log axes are Ch. 1 P13")
note("N71", "Other routes for complicated shapes", r"""
distributions of sources and sinks (C14, and N90's panels), thin-body perturbation theory (Ch. 14), and numerical
solution of Laplace's equation (here; Ch. 10).
""")
D("D22", ref="6.72")
nb.worked_example("four unknowns with ψ = xy on the boundary", r"""
Put the 4 × 4 grid of Fig. 6.23 on the points x, y ∈ {0, 1, 2, 3} (node (i, j) at x = i − 1, y = j − 1) and give the 12
boundary nodes the values of ψ = xy (a harmonic polynomial).

1. ψ₂₂ = ¼(ψ₁₂ + ψ₃₂ + ψ₂₁ + ψ₂₃) = ¼(0 + 2 + 0 + 2) = 1 = 1·1 ✓.
2. ψ₃₂ = ¼(ψ₂₂ + ψ₄₂ + ψ₃₁ + ψ₃₃) = ¼(1 + 3 + 0 + 4) = 2 ✓.
3. Likewise ψ₂₃ = 2, ψ₃₃ = 4. xy satisfies the average rule **exactly** (its fourth derivatives vanish, D22 step 3).
""")
note("N72", "The 16-point grid", r"""
12 known boundary values ψ^B and four unknowns give four linear equations —
""", equation=r"\begin{aligned}\psi_{2,2}&=\tfrac14[\psi^B_{1,2}+\psi_{3,2}+\psi^B_{2,1}+\psi_{2,3}],\\ \psi_{3,2}&=\tfrac14[\psi_{2,2}+\psi^B_{4,2}+\psi^B_{3,1}+\psi_{3,3}],\\ \psi_{2,3}&=\tfrac14[\psi^B_{1,3}+\psi_{3,3}+\psi_{2,2}+\psi^B_{2,4}],\\ \psi_{3,3}&=\tfrac14[\psi_{2,3}+\psi^B_{4,3}+\psi_{3,2}+\psi^B_{3,4}]\end{aligned}", ref="6.73")
nb.code(r"""
sysm = ch06.four_point_system("xy")                      # (6.73) with ψ^B = xy, unknowns ordered ψ22, ψ32, ψ23, ψ33
print(4*sysm["A"], 4*sysm["b"])                          # stored divided by 4; ×4 gives D23's A (4 on the diagonal) and b = (0, 3, 3, 12)
print(np.linalg.solve(sysm["A"], sysm["b"]))             # [1. 2. 2. 4.]: the values of xy
""")
nb.md("> 🔁 **Reminder:** `A @ x` is the matrix–vector product Ax in numpy (Ch. 2 P63).")
P("P160", "iterative solvers: Jacobi, Gauss–Seidel, SOR", r"""
A big linear system Aψ = b from a stencil is solved by repeated sweeps instead of elimination. **Jacobi** updates every
unknown from the *previous* sweep's values; **Gauss–Seidel** uses each new value as soon as it exists; **SOR**
(successive over-relaxation) takes the Gauss–Seidel change and multiplies it by ω (1 < ω < 2). The error shrinks each
sweep by a factor ρ < 1 (the spectral radius of the sweep's matrix): the closer ρ is to 1, the slower. Stop when the
**residual** b − Aψ (how badly the equations are violated) is small — and make the tolerance tight: when ρ is
close to 1 both the residual and the change sit far below the true error.
""", code=r"""
A = np.array([[4., -1], [-1, 4]]); b = np.array([3., 3])   # a 2 × 2 toy system, answer (1, 1)
x = np.zeros(2)                                          # start from zero
for k in range(4):                                       # four Gauss–Seidel sweeps
    x[0] = (b[0] + x[1])/4; x[1] = (b[1] + x[0])/4       # x[1] already uses the new x[0]
    print(k, x, np.max(np.abs(b - A @ x)))               # residual falls ×16 per sweep
""")
D("D23", ref="6.73")
note("N73", "Gauss–Seidel", r"""
sweep through the nodes, always using the latest value available (no book number). The book stops when values stop
changing; we test the residual of $\psi_{i,j}=\tfrac14[\dots]$ (6.72) with a tight tolerance, because on fine grids both it
and the change are far smaller than the true error
(D23 step 8). ω (the SOR factor) is not a vorticity here.
""", equation=r"\psi^{(k+1)}_{i,j}=\tfrac14\big[\psi^{(k+1)}_{i-1,j}+\psi^{(k)}_{i+1,j}+\psi^{(k+1)}_{i,j-1}+\psi^{(k)}_{i,j+1}\big]")
nb.code(r"""
for method in ("jacobi", "gauss_seidel", "sor"):         # the three sweeps of primer P160
    psi, hist = ls.solve_laplace(*ch06.four_point_system("xy", as_grid=True), method=method, tol=1e-10)   # to a residual of 1e-10
    print(f"{method:12s} {hist['sweeps']:3d} sweeps -> ψ22, ψ32, ψ23, ψ33 =", psi[1:3, 1:3].ravel())
print("ω_opt =", ls.optimal_sor_omega(4, 4))             # the SOR factor used: 2/(1 + sin(π/N)) with N = 3 grid intervals
""", explain=r"""
`four_point_system("xy", as_grid=True)` returns the 4 × 4 grid as a boolean mask of unknown nodes and an array of
boundary values (primer P161 below); `solve_laplace` sweeps until the residual is below 10⁻¹⁰. Jacobi needs 34 sweeps,
Gauss–Seidel 19, SOR 12 — all end at 1, 2, 2, 4.
""")
P("P161", "boolean masks and scipy.sparse", r"""
An L-shaped or stepped domain lives in a rectangular array: a boolean array `mask[j, i]` is True at the unknown nodes, and
the boundary values sit in another array. For big grids we skip iteration and hand the whole system to a sparse direct
solver: `scipy.sparse` stores only the non-zero entries (five per row here) and `scipy.sparse.linalg.spsolve` solves it.
""", code=r"""
import scipy.sparse as sps, scipy.sparse.linalg as spl   # sparse matrices and their direct solver
mask = np.ones((3, 4), bool); mask[0, :2] = False        # a small stepped domain: True = unknown node
print(mask.sum(), "unknowns")                            # 10
A = sps.diags([4.0, -1.0, -1.0], [0, 1, -1], shape=(4, 4), format="csr")   # a 1-D toy stencil, stored sparsely
print(spl.spsolve(A, np.ones(4)))                        # [0.3636 0.4545 0.4545 0.3636]
""")
note("N74", "Example 6.2: a sharp contraction", r"""
A channel narrows in one step. ψ = 0 on the lower wall and on the step, ψ = Q on the upper wall; across the inlet and the
outlet the velocity is uniform, so ψ rises linearly there (u = ∂ψ/∂y constant). The region has no holes, so the solution
is unique. At the step's corner the fluid turns through 270°: by $w(z)=Az^n$ *(6.46)* with α = 3π/2, n = ⅔, the speed
there is infinite (∝ r^{−1/3}) and ψ behaves like $r^{2/3}\sin(2\theta/3)$. That singularity **pollutes the whole grid**:
refining the grid, the solution converges at order **4/3 everywhere**, not the order 2 of a smooth problem (the refinement
cell shows it). We run it with Q = 1 m²/s and show ψ/Q (the book's flow rate and printed grid values stay in the private
test file; our Gauss–Seidel reproduces them). ⚠️ The book's FORTRAN listing loops over I where it should loop over J
(harmless there), and prints the inlet Δψ in m² — it is m²/s.
""")
nb.code(r"""
ex = ch06.example_6_2(Q=1.0)                             # the contraction on the book's grid, Gauss–Seidel to 1e-10
print("grid", ex["grid_shape"], "| sweeps", ex["history"]["sweeps"], "| final residual", f"{ex['history']['residual'][-1]:.1e}")
vals = []                                                # ψ/Q at one fixed interior point, on finer and finer grids
for r in (1, 2, 4, 8, 16):                               # refinement factors
    e2 = ch06.example_6_2(Q=1.0, refine=r, method="direct")   # the same problem, solved directly (spsolve)
    vals.append(e2["psi"][e2["probe"]])                  # `probe` = the same physical point at every refinement
    print(f"refine ×{r:2d}: grid {e2['grid_shape']}, ψ/Q at the probe = {vals[-1]:.6f}")
d = np.abs(np.diff(vals))                                # changes between successive grids
print("changes:", d, "| ratios:", d[:-1]/d[1:], "| pairwise orders:", np.log2(d[:-1]/d[1:]))   # → 2^(4/3) ≈ 2.52, order 4/3
""", explain=r"""
1. `example_6_2` builds the contraction on the book's grid (inlet, step and outlet as in Fig. 6.24, with our Q = 1 m²/s)
   and relaxes it by Gauss–Seidel until the residual is below 10⁻¹⁰.
2. The same problem on grids 2, 4, 8 and 16 times finer is solved directly (primer P161).
3. The value at a fixed interior point settles, but slowly: the change shrinks by a factor falling toward 2^{4/3} ≈ 2.52
   per halving (pairwise orders 1.66, 1.50, 1.43, … → 4/3), not the factor 4 of a second-order scheme — the 270° corner's
   r^{2/3} singularity sets the accuracy of the whole grid (verified to order 1.36 at three points in
   `tests/test_ch06.py`).
""")
nb.md("**From scratch — a double-loop Gauss–Seidel on the 4-point problem:**")
nb.check_agree(r"""
P_ = np.array([[i*j for i in range(4)] for j in range(4)], float)   # ψ = xy on the 4 × 4 grid (row j = y, column i = x)
P_[1:3, 1:3] = 0.0                                       # the four unknowns start at zero
for sweep in range(40):                                  # 40 sweeps
    for j in (1, 2):                                     # rows (y) …
        for i in (1, 2):                                 # … columns (x), i fastest: the book's order
            P_[j, i] = 0.25*(P_[j, i-1] + P_[j, i+1] + P_[j-1, i] + P_[j+1, i])   # (6.72), newest values used at once
print(P_[1:3, 1:3])                                      # [[1, 2], [2, 4]]
assert np.allclose([P_[1, 1], P_[1, 2], P_[2, 1], P_[2, 2]], np.linalg.solve(sysm["A"], sysm["b"]))   # ψ22, ψ32, ψ23, ψ33
""")
nb.md("Four lines of loops do what `solve_laplace` does.")
nb.figure(r"""
big = ch06.example_6_2(Q=1.0, refine=4, method="direct") # a fine grid for the picture
hists = {m: ch06.example_6_2(Q=1.0, method=m)["history"] for m in ("jacobi", "gauss_seidel", "sor")}   # book grid, three sweeps
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.6))  # streamlines | convergence
Pm = np.where(big["mask"] | ~np.isnan(big["psi"]), big["psi"], np.nan)   # ψ/Q where it is defined
a1.contour(big["X"], big["Y"], Pm, levels=np.linspace(0.1, 0.9, 9), colors=COLORS["blue"], linewidths=1)   # ψ/Q = 0.1 … 0.9
a1.contourf(big["X"], big["Y"], np.isnan(Pm).astype(float), levels=[0.5, 1.5], colors=["#c9ceda"])   # outside the channel
a1.set_aspect("equal"); a1.set_xlabel("x [grid units]"); a1.set_ylabel("y [grid units]")
a1.set_title("contraction: ψ/Q = 0.1 … 0.9", fontsize=10)
for m, col in (("jacobi", COLORS["muted"]), ("gauss_seidel", COLORS["accent"]), ("sor", COLORS["teal"])):
    h_ = hists[m]                                        # residual and change per sweep
    a2.semilogy(np.arange(1, len(h_["residual"]) + 1), h_["residual"], color=col, lw=2, label=f"{m}: residual")
    a2.semilogy(np.arange(1, len(h_["change"]) + 1), h_["change"], color=col, lw=1, ls="--", label=f"{m}: change")
a2.set_xlabel("sweep"); a2.set_ylabel("max residual / max change [ψ/Q]"); a2.legend(fontsize=7)
fig.suptitle("Streamlines by averaging; the smarter sweep wins", fontweight="bold")
savefig(fig, "ch06", "c12_contraction"); plt.show()      # save, then draw
""", see="Left, streamlines squeezing from the wide inlet into the narrow outlet round the step; right, straight lines on "
         "the log axis with different slopes — solid residuals, dashed changes.",
    read="A straight line on a log axis means a constant factor per sweep, the spectral radius; Gauss–Seidel's slope is "
         "about twice Jacobi's (ρ_GS ≈ ρ_J²); SOR is steeper still. Residual and change fall together, and on fine grids "
         "both sit far below the true error — hence the tight 10⁻¹⁰ tolerance.",
    change="…the grid refined ×2: Jacobi and Gauss–Seidel need ≈ 4× the sweeps (∝ N²), SOR ≈ 2× (∝ N).")
nb.animation(r"""
ks = list(range(1, 9)) + [12, 16, 24, 32, 48] if not FAST else [1, 2, 3, 4, 6, 8, 12, 20, 32, 48]   # sweep counts to show
snaps = [ch06.example_6_2(Q=1.0, n_iter=k)["psi"] for k in ks]   # ψ/Q after k Gauss–Seidel sweeps from zero
res = ch06.example_6_2(Q=1.0)["history"]["residual"]    # the residual after each sweep
geo = ch06.example_6_2(Q=1.0)                            # grid coordinates
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.5, 3.0))   # grid | residual
im = a1.imshow(snaps[0], origin="lower", cmap="Blues", vmin=0, vmax=1)   # ψ/Q on the nodes
a1.set_xlabel("i"); a1.set_ylabel("j"); fig.colorbar(im, ax=a1, label="ψ/Q")
a2.semilogy(np.arange(1, len(res) + 1), res, color=COLORS["muted"], lw=1)   # the whole residual history, faint
dot, = a2.semilogy([1], [res[0]], "o", color=COLORS["accent"])   # where we are
a2.set_xlabel("sweep"); a2.set_ylabel("max residual")
def update(i):                                           # frame i shows sweep ks[i]
    im.set_data(snaps[i])                                # the relaxed values
    dot.set_data([ks[i]], [res[min(ks[i], len(res)) - 1]])   # move the dot
    a1.set_title(f"after {ks[i]} Gauss–Seidel sweeps", fontsize=10)
show_animation(animate(update, frames=len(ks), fig=fig), player="frames")   # step through with the buttons
""")
see_read_change(
    see="The grid starts almost empty; the boundary values leak inward sweep by sweep until the colours stop changing; the "
        "dot runs down the residual line.",
    read="Each sweep replaces every value by the average of its neighbours (newest first); information spreads about one "
         "node per sweep, so large grids need many sweeps.",
    change="…SOR with ω_opt: the colours settle in about a third of the sweeps.")
nb.pointer("§6.7 also names source panels — a boundary-only numerical method; our constant-strength source panels (note "
           "N90) are taught with the axial singularity method in C14 (§6.8).")
explainer("laplace_relaxation", "How does a grid solve Laplace's equation?",
          "the method is a process: stepping one node at a time, then whole sweeps, with the residual curve growing "
          "beside the grid, makes the average rule and the difference between Jacobi, Gauss–Seidel and SOR visible.", "", [
              "Step node by node with Gauss–Seidel on the 4-point problem: the first values are 0, 0.75, 0.75, 3.375.",
              "Switch to Jacobi and compare the residual slopes.",
              "Preset 'SOR at ω_opt' on the refined contraction: an order of magnitude fewer sweeps.",
              "Preset 'ψ = xy': one solve is exact — the harmonic polynomial satisfies (6.72) exactly.",
          ])
whatif(r"""
…the right-hand side were not zero — the vorticity of a rotating flow, $\omega_z=-\nabla^2\psi$ *(6.4)*? The same sweeps
solve Poisson's equation with ωh²/4 added to each average (`ls.solve_poisson`). That is how Ch. 10's
vorticity–stream-function method, and ocean models that invert potential vorticity in a closed basin, recover ψ every time
step.
""")

# =====================================================================================================================
# A.8 §6.8 — R22–R30, C13 (N75–N85, D24 D25), C14 (N86–N90, D26 D27, E8)
# =====================================================================================================================
nb.section("6.8", "Axisymmetric Ideal Flow", intro=r"""
**What is this section about?** Flow round a sphere, a bullet or an airship is the same in every plane through the axis.
One stream function — Stokes' — still describes it, but its equation is *not* Laplace's, so complex variables no longer
help; elements are built with φ and ψ in spherical coordinates instead. ⚠️ In this section the symmetry axis z is
horizontal, along the stream.
""")
nb.recap("R22", "Two stream functions in 3-D", r"""
A general 3-D incompressible flow needs two stream functions, $\rho\mathbf u=\nabla\chi\times\nabla\psi$ *(Eq. 4.12)*; an
axisymmetric one needs only one.
""", where="Ch. 4 §4.3")
nb.recap("R23", "Axisymmetric continuity", r"""
In cylindrical coordinates (R, φ, z) with no swirl and no φ-dependence:
$\frac1R\frac{\partial}{\partial R}(Ru_R)+\frac{\partial u_z}{\partial z}=0$ *(Eq. 6.74)*. ⚠️ z is now horizontal, and φ
here is the azimuth angle round the axis (the potential is also written φ — the context tells which).
""", where="Ch. 4 §4.3")
nb.recap("R24", "The Stokes stream function", r"""
Choosing χ = −φ (the azimuth angle) in the two-stream-function form gives
$u_R=-\frac1R\frac{\partial\psi}{\partial z},\ u_z=\frac1R\frac{\partial\psi}{\partial R}$ *(Eq. 6.75)*, which satisfies
continuity identically. Same signs as Ch. 4's `velocity_from_streamfunction_axisym`.
""", where="Ch. 4 §4.3")
nb.code(r"""
print(sf.velocity_from_streamfunction_axisym(lambda R, z: 0.5*2.0*R**2, 0.3, 0.1))   # ψ = ½UR², U = 2 m/s: (u_R, u_z) = (0, 2)
""")
nb.recap("R25", "The azimuthal vorticity", r"""
$\omega_\varphi=\frac{\partial u_R}{\partial z}-\frac{\partial u_z}{\partial R}$ *(Eq. 6.76)* — Appendix B's curl for an
axisymmetric field.
""", where="App. B; operators recapped in Ch. 4")
nb.recap("R26", "The axisymmetric potential", r"""
$u_R=\partial\phi/\partial R,\ u_z=\partial\phi/\partial z$ *(Eq. 6.79)* — u = ∇φ, which exists when
$\boldsymbol\omega=0$ *(Eq. 3.17)*, written in cylindrical coordinates.
""", where="Ch. 3 §3.4")
nb.recap("R27", "Cylindrical and spherical coordinates", r"""
(R, φ, z) and (r, θ, φ) with z along the axis: $R=r\sin\theta,\ z=r\cos\theta$ *(Eq. 6.81)*, θ measured from +z
(downstream) — the table the book gives.
""", where="Ch. 3 §3.1")
nb.recap("R28", "Spherical continuity", r"""
$\frac1r\frac{\partial}{\partial r}(r^2u_r)+\frac1{\sin\theta}\frac{\partial}{\partial\theta}(u_\theta\sin\theta)=0$
*(Eq. 6.82)* — r times Appendix B's divergence
$\frac1{r^2}\frac{\partial}{\partial r}(r^2u_r)+\frac1{r\sin\theta}\frac{\partial}{\partial\theta}(u_\theta\sin\theta)$; the
factor r does not change "= 0".
""", where="App. B; operators recapped in Ch. 4")
nb.recap("R29", "Spherical vorticity", r"""
$\omega_\varphi=\frac1r\Big[\frac{\partial}{\partial r}(ru_\theta)-\frac{\partial u_r}{\partial\theta}\Big]$ *(Eq. 6.84)*.
""", where="App. B; operators recapped in Ch. 4")
nb.recap("R30", "Spherical axisymmetric Laplace", r"""
$\frac1{r^2}\frac{\partial}{\partial r}\Big(r^2\frac{\partial\phi}{\partial r}\Big)+\frac1{r^2\sin\theta}\frac{\partial}{\partial\theta}\Big(\sin\theta\frac{\partial\phi}{\partial\theta}\Big)=0$
*(Eq. 6.85)* — used in D25's check. Ch. 4's `core.curvilinear` builds this operator symbolically.
""", where="App. B; operators recapped in Ch. 4")
nb.code(r"""
from fluidpy.core import curvilinear as cu               # Ch. 4's symbolic operators in curvilinear coordinates
r_, t_, ph_ = cu.coordinates("spherical")                # the module's own symbols r, θ, φ
print(cu.laplacian(r_*sp.cos(t_), "spherical"), cu.laplacian(sp.cos(t_)/r_**2, "spherical"))   # 0 0: stream and 3-D doublet φ are harmonic
""")

# ---- C13 ------------------------------------------------------------------------------------------------------------
core("C13", "Axisymmetric potential flow: the Stokes stream function, 3-D elements and the sphere", r"""
Air flows round a football and round a lamp-post. What changes when the flow can also go *over* the body?
""", eqs=("6.90",))
problem(r"""
A lamp-post is effectively 2-D: the air can only go round the sides. A football is 3-D: the air can go round, over and
under, every way at once. That extra freedom should make the disturbance smaller and the fastest speed lower. We want the
numbers — and we find along the way that the tidy complex-variable tools of §6.4–§6.6 do not carry over.
""")
nb.md(r"""
#### The idea

| | cylinder (2-D) | sphere (3-D) |
|---|---|---|
| built from | stream + 2-D doublet | stream + 3-D doublet |
| fastest surface speed | 2U | 1.5U |
| minimum C_p | −3 | −1.25 |
| disturbance decays like | (a/r)² | (a/r)³ |
| stream-function equation | Laplace | (6.77), not Laplace |

**The fluid escapes sideways, so it needs to speed up less.**
""")
remind([
    ("cylindrical and spherical unit vectors", "e_R, e_φ, e_z and e_r, e_θ, e_φ point along increasing coordinates; they change direction from point to point (Ch. 3 P88)."),
], lead=r"""**cylindrical cross products:** $\mathbf e_R\times\mathbf e_\varphi=\mathbf e_z$, $\mathbf e_\varphi\times\mathbf e_z=\mathbf e_R$,
$\mathbf e_z\times\mathbf e_R=\mathbf e_\varphi$ (right-handed, cyclic), so $\mathbf e_\varphi\times\mathbf e_R=-\mathbf e_z$;
**flux of a 3-D point source:** a source of Q m³/s sends Q through every sphere round it, $4\pi r^2u_r=Q$; the 3-D
Taylor series is Ch. 3 P98 and Poisson/Green in 3-D Ch. 5 P139""")
D("D24", ref="6.77")
note("N75", "The field equation of the Stokes stream function", r"""
in irrotational flow. It is not the Laplacian of ψ, so ψ is not the imaginary part of an analytic function: no complex
potential in axisymmetric flow.
""", equation=EQ["6.77"], ref="6.77")
nb.code(r"""
sph = pf.sphere(1.0, 1.0)                                # the sphere in a stream: U = 1 m/s, a = 1 m
Rq = np.array([0.5, 2.0, 1.5]); Zq = np.array([2.0, -1.0, 0.8])   # three points (R, z) in the fluid [m]
print(ch06.stokes_operator_residual(sph.psi, Rq, Zq))    # ≈ 0: the sphere's ψ satisfies (6.77)
print(ch06.stokes_operator_residual(lambda R, z: R*z, Rq, Zq))   # a control ψ = Rz: residual −z/R², clearly not 0
""")
note("N76", "Units", r"""
The Stokes ψ is in m³/s (a volume flow), the plane ψ in m²/s (per unit depth); ψ = const is a surface of revolution.
""")
note("N77", "Flow between two stream surfaces", r"""
⚠️ 2π, not 1: dψ is the flow per radian round the axis. The cell compares a quadrature of 2πR u·n ds along a segment
with 2π(ψ₂ − ψ₁):
""", equation=r"dQ=2\pi R(\mathbf u\cdot\mathbf n)ds=2\pi R(-u_Rdz+u_zdR)=2\pi\,d\psi", ref="6.78")
nb.code(r"""
print(ch06.axisym_flux_between(sph.psi, (1.5, 0.0), (3.0, 0.5)))   # (quadrature, 2πΔψ): two equal numbers [m³/s]
""")
note("N78", "The potential obeys the axisymmetric Laplace equation", "", equation=EQ["6.80"], ref="6.80")
note("N79", "In spherical coordinates", "", equation=EQ["6.83"], ref="6.83")
nb.code(r"""
print(pf.axisym_velocity_spherical_sym(sp.Rational(1, 2)*r_**2*sp.sin(t_)**2, r_, t_))   # ψ = ½r²sin²θ (U = 1): (cos θ, −sin θ)
""")
D("D25", ref="6.91")
note("N80", "The 3-D elements (notes `N81`, `N82`)", r"""
Uniform flow along z: $\phi=Uz,\ \psi=\tfrac12UR^2$ or $\phi=Ur\cos\theta,\ \psi=\tfrac12Ur^2\sin^2\theta$ *(6.86)*; point
source Q [m³/s]: $\phi=-\frac{Q}{4\pi\sqrt{R^2+z^2}},\ \psi=-\frac{Qz}{4\pi\sqrt{R^2+z^2}}$ or
$\phi=-\frac Q{4\pi r},\ \psi=-\frac Q{4\pi}\cos\theta$ *(6.87)*; doublet with dipole −d e_z [m⁴/s]:
$\phi=\frac d{4\pi}\frac{z}{(R^2+z^2)^{3/2}},\ \psi=-\frac d{4\pi}\frac{R^2}{(R^2+z^2)^{3/2}}$ or
$\phi=\frac d{4\pi r^2}\cos\theta,\ \psi=-\frac d{4\pi r}\sin^2\theta$ *(6.88)*. In a plane through the axis their
streamlines look like the 2-D ones.
""")
note("N83", "The sphere", r"""
= stream + an opposing doublet of strength d = 2πa³U. ψ = 0 on the whole axis and on r = a. Outside, it is the flow round
Hill's spherical vortex seen from the vortex (Ch. 5 §5.4).
""", equation=r"\psi=\tfrac12Ur^2\sin^2\theta-\frac{d}{4\pi r}\sin^2\theta=\tfrac12Ur^2\Big(1-\frac{a^3}{r^3}\Big)\sin^2\theta;\quad \phi=Ur\Big(1+\frac{a^3}{2r^3}\Big)\cos\theta", ref="6.89")
nb.code(r"""
print(sph.psi(1.3, 0.4), 0.5*1.3**2*(1 - 1/np.hypot(1.3, 0.4)**3))   # U = a = 1: ψ = ½R²(1 − a³/r³) since r² sin²θ = R²
""")
note("N84", "Surface speed and pressure", r"""
With $u_r=U[1-(a/r)^3]\cos\theta$, $u_\theta=-U[1+\frac12(a/r)^3]\sin\theta$ *(6.90)* the surface speed is (3/2)U sin θ,
and the pressure below follows. Fore–aft symmetric: no drag again (the 3-D d'Alembert paradox, an exercise of S06).
""", equation=r"C_p=\frac{p-p_\infty}{\frac12\rho U^2}=1-\Big(\frac{u_\theta}{U}\Big)^2=1-\frac94\sin^2\theta", ref="6.91")
nb.code(r"""
print(ch06.cylinder_vs_sphere(2.0))                      # at r = 2a: disturbance 0.25 vs 0.125; C_p min −3 vs −1.25; max speed 2U vs 1.5U
""")
note("N85", "Without coordinates", r"""
With x = r e_r, cos θ = e_z·e_r, U = U e_z and d = −d e_z the potential needs no coordinates — valid for a stream in any
direction, the starting point of C15:
""", equation=r"\phi=Ur\cos\theta+\frac{d}{4\pi r^2}\cos\theta=\mathbf U\cdot\mathbf x-\frac{\mathbf d\cdot\mathbf x}{4\pi\lvert\mathbf x\rvert^3}=\Big(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3}\Big)\cdot\mathbf x", ref="6.92")
nb.code(r"""
xv = np.array([0.3, -1.2, 0.9])                          # a point in 3-D [m]
print(pf.sphere_potential_vector(xv, [0, 0, 1.0], [0, 0, -2*np.pi]), sph.phi(np.hypot(0.3, -1.2), 0.9))   # equal (U = 1, a = 1)
""")
nb.worked_example("a football and a lamp-post, both a = 0.11 m, in wind U = 10 m/s", r"""
1. Fastest surface speed: lamp-post 2U = 20 m/s; football 1.5U = 15 m/s.
2. Lowest pressure, air ½ρU² = 60 Pa: lamp-post C_p = −3 → −180 Pa; football −1.25 → −75 Pa.
3. At twice the radius (r = 0.22 m) the stream is disturbed by (a/r)² = 25 % past the lamp-post but only (a/r)³ = 12.5 %
   past the ball.
4. Both have zero drag in ideal flow.
""")
nb.code(r"""
sph2 = pf.sphere(10.0, 0.11)                             # a football-sized sphere in a 10 m/s wind
th = np.radians([0, 45, 90, 135, 180])                   # five surface angles from the downstream axis [rad]
print(sph2.velocity_spherical(0.11, th))                 # (u_r, u_θ): u_r ≈ 0, u_θ = −15 sin θ m/s
print(ch06.sphere_surface_cp(th))                        # (6.91): 1, −0.125, −1.25, −0.125, 1
print(np.max(np.abs(sph2.psi(0.11*np.sin(th), 0.11*np.cos(th)))))   # ≈ 0: the sphere is the stream surface ψ = 0
""", explain=r"""
1. The sphere as an axisymmetric flow (stream + 3-D doublet).
2. On the surface the velocity is tangential, 1.5U = 15 m/s at the equator.
3. $C_p=1-\frac94\sin^2\theta$ at five angles.
4. ψ = 0 on the surface.
""")
nb.md("**From scratch — stream plus doublet, added by hand:**")
nb.check_agree(r"""
U, a = 10.0, 0.11                                        # stream [m/s] and radius [m]
d = 2*np.pi*a**3*U                                       # the doublet that makes r = a a stream surface (D25 step 7) [m⁴/s]
r, t = 0.3, 0.7                                          # a point (r [m], θ [rad])
psi_hand = 0.5*U*r**2*np.sin(t)**2 - d/(4*np.pi*r)*np.sin(t)**2   # (6.86) + (6.88)
print(psi_hand)                                          # [m³/s]
assert np.isclose(psi_hand, sph2.psi(r*np.sin(t), r*np.cos(t)))   # = the library's sphere
""")
nb.md("Stream plus doublet, added by hand, is the library's sphere.")
nb.figure(r"""
n = 201 if not FAST else 141                             # grid points per side
Zg, Rg = np.meshgrid(np.linspace(-3, 3, n), np.linspace(-2.2, 2.2, n))   # meridian plane: z horizontal, R vertical [m]
Ps = sph.psi(np.abs(Rg), Zg)*np.sign(Rg)                 # Stokes ψ, mirrored below the axis (U = a = 1)
Ps = np.where(np.hypot(Rg, Zg) < 1.0, np.nan, Ps)        # blank inside the sphere
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.7))  # streamlines | pressure
a1.contour(Zg, Rg, Ps, levels=np.linspace(-2.2, 2.2, 23), colors=COLORS["ink"], linewidths=0.7)   # equal Δψ
draw_body(a1, circle(1.0)); a1.set_aspect("equal"); a1.set_xlabel("z [m] (along the stream)"); a1.set_ylabel("R [m]")
a1.set_title("sphere, U = 1 m/s, a = 1 m (meridian plane)", fontsize=10)
thd = np.linspace(0, 180, 181)                           # surface angle [deg]
a2.plot(thd, 1 - 4*np.sin(np.radians(thd))**2, "--", color=COLORS["muted"], lw=2, label=r"cylinder: $1-4\sin^2\theta$ (6.35)")
a2.plot(thd, ch06.sphere_surface_cp(np.radians(thd)), color=COLORS["accent"], lw=2.5, label=r"sphere: $1-\frac{9}{4}\sin^2\theta$ (6.91)")
a2.set_xlabel("θ [deg]"); a2.set_ylabel("$C_p$ [–]"); a2.legend(fontsize=8, loc="lower left")
ins = a2.inset_axes([0.38, 0.72, 0.24, 0.21])             # inset: disturbance at the side against r/a
ra = np.logspace(0, 1, 30)                               # r/a from 1 to 10
ins.loglog(ra, 1/ra**2, "--", color=COLORS["muted"]); ins.loglog(ra, 0.5/ra**3, color=COLORS["accent"])   # cylinder, sphere
ins.set_title("|u_θ/U| − 1 at θ = 90°", fontsize=7); ins.tick_params(labelsize=6); ins.set_xlabel("r/a", fontsize=7)
fig.suptitle("In 3-D the flow escapes sideways", fontweight="bold")
savefig(fig, "ch06", "c13_sphere"); plt.show()           # save, then draw
""", see="Left, fore–aft symmetric streamlines round the sphere; right, the sphere's C_p dips only to −1.25 against the "
         "cylinder's −3; the inset shows two straight lines with slopes −2 (cylinder) and −3 (sphere).",
    read="The sphere's disturbance is weaker and dies faster: surface speed 1.5U vs 2U, perturbation (a/r)³ vs (a/r)².",
    change="…reading speed from the spacing of equal-Δψ Stokes streamlines: the flux between two of them is 2πΔψ spread "
           "over a ring of circumference 2πR, so near the axis they crowd without the flow being fast — read speed from "
           "spacing × R, not spacing alone.")
nb.plotly(r"""
thd = np.linspace(0, 180, 181)                           # angle from the downstream axis [deg]
tr = np.radians(thd)                                     # the same in radians
def rel_traces(ra):                                      # |u_θ|/U round a circle / sphere of radius r = ra·a (a = 1)
    return {"cylinder |u_θ|/U": (thd, (1 + 1/ra**2)*np.abs(np.sin(tr))),        # from (6.34)
            "sphere |u_θ|/U": (thd, (1 + 0.5/ra**3)*np.abs(np.sin(tr))),        # from (6.90)
            "far stream": (thd, np.abs(np.sin(tr)))}                            # the undisturbed stream's θ-part
ras = np.linspace(1, 5, 17)                              # r/a values
fig = slider_figure(rel_traces, "r/a", ras, xlabel="θ [deg]", ylabel="|u_θ|/U [–]",
                    title="3-D relief: the sphere's disturbance dies as (a/r)³", yrange=[0, 2.1])
recolor(fig, {"cylinder |u_θ|/U": COLORS["muted"], "sphere |u_θ|/U": COLORS["accent"], "far stream": COLORS["ink"]},
        dashes={"cylinder |u_θ|/U": "dash", "far stream": "dot"})
step_titles(fig, [f"3-D relief: at r = {v:.2f}a the side speed is {1 + 1/v**2:.3f}U (cylinder) vs {1 + 0.5/v**3:.3f}U (sphere)" for v in ras])
fig.show()                                               # draw it (works on the web page too)
""")
see_read_change(
    see="Two arches over the dotted far-stream curve: at r = a the cylinder's reaches 2, the sphere's 1.5; dragging r/a "
        "outward, the sphere's arch collapses onto the far stream much sooner.",
    read="The extra speed is the disturbance of the body: (a/r)² for the cylinder, ½(a/r)³ for the sphere (the title's "
         "numbers).",
    change="…a body of revolution longer than a sphere (an airship, C14): the side disturbance is smaller still.")
whatif(r"""
…the stream were not along the axis but across it? The flow would no longer be axisymmetric about the stream — but the
sphere is round, so $\phi=\big(\mathbf U-\frac{\mathbf d}{4\pi\lvert\mathbf x\rvert^3}\big)\cdot\mathbf x$ *(6.92)* handles
any direction of U at once. For a non-spherical body we would need the full 3-D problem (Ch. 10 methods, or the
singularity distributions of C14).
""")

# ---- C14 ------------------------------------------------------------------------------------------------------------
core("C14", "Bodies of revolution from axial singularities: the airship and the axial singularity method", r"""
An airship hull is drawn on paper. Which sources and sinks, hidden on its axis, produce exactly that shape?

*In one line:* $\psi_m=-\sum_n\frac{k_n}{4\pi}(r^m_{n-1}-r^m_n)+\tfrac12UR_m^2=0$ at N body points (the book's Fig. 6.29
system).
""")
problem(r"""
So far we picked sources and found the body. A designer works the other way: the shape of a hull or a fuselage is given,
and the flow round it is wanted. Hide unknown sources and sinks along the axis, demand that the body be a stream surface
at as many points as there are unknowns, and solve — a linear system. The same "unknown strengths + one condition per
point" idea runs every panel code in aerodynamics.
""")
idea("""
point source Q at the nose + line sink k over length a + stream U   -->  closed airship when Q = ak  (6.95)
inverse: N axial segments of unknown k_n --psi = 0 at N body points--> N x N linear system --> k_n
closure check: sum of k_n * d_xi ~ 0 (sources and sinks cancel, or the body does not close)
""")
note("N86", "Closed bodies need zero net source", r"""
Everything the sources emit must be swallowed again, or fluid would stream out of the body forever (the half-body of C05
is open for that reason). A sink spread over a length gives a gently tapering tail.
""")
remind([
    ("substitution in an integral", r"replace the variable, $d\xi=\frac{d\xi}{d\alpha}d\alpha$, and move the limits with it (Ch. 3 P106)."),
    ("scipy.integrate.quad and dblquad", "`quad(f, a, b)` integrates a function numerically to about 1e-10 (Ch. 3 P87)."),
], lead=r"""**the substitution z − ξ = R cot α** with z, R fixed: d(cot α)/dα = −1/sin²α, and this minus cancels the one from
−dξ, so dξ = R dα/sin²α and the limits keep their order (ξ = 0 → a becomes α = θ → α₁); the point-source ψ is note N81
and `brentq` Ch. 3 P108""")
D("D26", ref="6.95")
note("N87", "The line sink and the airship (notes `N88`, `N89`)", r"""
An element k dξ of a line sink of density k [m²/s] on the axis contributes dψ = (k dξ/4π) cos α (the opposite sign of a
source (6.87)), so $\psi_{\text{sink}}=\frac k{4\pi}\int_0^a\cos\alpha\,d\xi$ *(6.93)*, which the substitution turns into
$\psi_{\text{sink}}=\frac{kR}{4\pi}\Big[\frac1{\sin\theta}-\frac1{\sin\alpha_1}\Big]=\frac k{4\pi}(r-r_1)$ *(6.94)*; with a
point source Q at O, k = Q/a and the stream, the airship is the ψ below. Q = 1 m³/s, a = 1 m → k = 1 m²/s.
""", equation=EQ["6.95"], ref="6.95")
nb.worked_example("a small airship", r"""
Q = 1 m³/s, a = 1 m, U = 1 m/s.

1. Closure: k = Q/a = 1 m²/s, so the sink swallows 1 × 1 = 1 m³/s = Q ✓.
2. Nose: on the axis ahead of O the source's push Q/4πz² nearly balances U; with the line sink's pull the stagnation
   point is at z = −0.252 m (a pure source would give −√(Q/4πU) = −0.282 m).
3. Tail: just behind the sink's end, z = 1.070 m.
4. Length 1.322 m — `brentq` on u_z(z) = 0 along the axis.
""")
nb.code(r"""
Rq = np.array([0.1, 0.3]); Zq = np.array([0.4, -0.5])     # two points (R, z) [m]
print(ch06.line_sink_stream_function(Rq, Zq, 1.0, 1.0), ch06.line_sink_stream_function(Rq, Zq, 1.0, 1.0, method="quad"))   # (6.94) vs quad of (6.93)
A_s = ch06.airship(1.0, 1.0, 1.0)                        # U = 1 m/s, Q = 1 m³/s, a = 1 m
print(f"nose {A_s['stagnation_front']:.4f} m, tail {A_s['stagnation_rear']:.4f} m, length {A_s['length']:.4f} m, "
      f"max radius {A_s['R_max']:.4f} m, closure Q − ak = {A_s['closure']}")
""", explain=r"""
1. The closed form (6.94) against a numerical integral of (6.93): equal.
2. The airship's nose and tail from the axis stagnation points (`brentq` on u_z = 0), its thickness from the ψ = 0
   surface, and the closure Q − ak = 0.
""")
P("P162", "collocation and the condition number", r"""
To pin down N unknown strengths, demand the condition (here ψ = 0) at N chosen points: N equations, N unknowns —
**collocation**. How trustworthy the answer is depends on the matrix: `np.linalg.cond(A)` is the factor by which relative
errors in the data can grow in the solution. Around 10³ is harmless; above 10¹⁰ the digits are noise. `np.vander` below
builds a matrix of powers, famous for being ill-conditioned.
""", code=r"""
for N in (4, 8, 16):                                     # three sizes
    x = np.linspace(0, 1, N)                             # N collocation points
    V = np.vander(x, N)                                  # a notoriously ill-conditioned matrix (powers of x)
    print(N, f"{np.linalg.cond(V):.1e}")                 # grows explosively with N
""")
D("D27")
nb.code(r"""
zb, Rb = ch06.axisym_body_target("rankine_oval", N=40)   # 40 points on a target body (a 3-D Rankine oval) [m]
sol = pn.axial_singularity_solve(zb, Rb, U=1.0, N=40)    # the N × N collocation system of D27, solved
print(f"cond {sol['cond']:.1e}, net strength Σk_nΔξ = {sol['net_strength']:.1e} m³/s, max residual {np.max(np.abs(sol['residual'])):.1e}")
for N in (10, 20, 40):                                   # more segments (N = 40 is the ceiling: beyond it cond(A) > 10¹²)
    s = ch06.axial_state("rankine_oval", N=N)            # body error, net strength, condition number, strength range
    print(f"N = {N:2d}: body error {s['body_error']:.1e}, Σk_nΔξ {s['net_strength']:+.1e}, cond {s['cond']:.1e}, "
          f"k from {s['k_min']:.1f} to {s['k_max']:.1f} m²/s")
print("sphere target:", ch06.axial_state("sphere", N=40))   # a blunt body: tiny residual, but cond ≈ 1e18
""", explain=r"""
1. `axisym_body_target` samples a target body; `axial_singularity_solve` builds the influence matrix of D27 step 4 and
   solves it: the residual at the collocation points is round-off.
2. With more segments the body is drawn more accurately (the body error falls from 2.5×10⁻³ to 1.8×10⁻⁶) and the net
   strength Σk_nΔξ stays ≈ 0 — the body closes. **But** the individual strengths k_n alternate in sign and grow (about
   ±5, ±10, ±19 m²/s): only their sums and moments (net source, dipole) converge, not the bars themselves; the condition
   number climbs to ~10¹², so N = 40 is a sensible ceiling here.
3. The sphere: the residual is still small, but the condition number is ~10¹⁸ and the strengths are huge — a sphere's
   exact "source" is a point doublet, which a smooth line of sources cannot copy. Qualitative only.
""")
nb.md("**From scratch — the influence matrix by our own double loop:**")
nb.check_agree(r"""
N = 20                                                   # number of segments
zb, Rb = ch06.axisym_body_target("rankine_oval", N=N)    # 20 body points [m]
xi = pn.axial_singularity_solve(zb, Rb, U=1.0, N=N)["z_nodes"]   # the library's segment ends (N + 1 nodes) [m]
A_ = np.zeros((N, N))                                    # influence matrix
for m_ in range(N):                                      # body point m
    for n_ in range(N):                                  # segment n
        r0 = np.hypot(zb[m_] - xi[n_], Rb[m_])           # distance to the segment's start [m]
        r1 = np.hypot(zb[m_] - xi[n_ + 1], Rb[m_])       # distance to its end [m]
        A_[m_, n_] = (r0 - r1)/(4*np.pi)                 # unit source segment, minus sign moved to the right (D27 step 4)
k = np.linalg.solve(A_, 0.5*1.0*Rb**2)                   # strengths from ψ = 0 at the 20 points, U = 1 m/s [m²/s]
assert np.allclose(k, pn.axial_singularity_solve(zb, Rb, U=1.0, z_nodes=xi)["k"])   # = the library's strengths
print(np.round(k[:6], 3))                                # the first few: alternating signs — only their sums and moments converge
""")
nb.md("Our own N × N matrix gives the library's strengths.")
note("N90", "The same idea on the body surface (our extension, not in the book)", r"""
Cover a 2-D body with N straight panels carrying unknown uniform source strengths λ_j, demand u·n = 0 at each panel's
midpoint, and solve (Hess & Smith 1967). On a **circle** these constant-strength panels are exact at every N (the midpoint
C_p equals 1 − 4 sin²θ of (6.35) to round-off); on an **ellipse** the midpoint C_p converges with error ∝ 1/N² (verified
order 2.06); Σλ_j s_j = 0 (a closed body). Off the body the panel velocity converges only at first order.
""")
nb.code(r"""
Ns = np.array([8, 16, 32, 64])                           # numbers of panels
e_ell = [ch06.panel_cp_error(N) for N in Ns]             # max C_p error at the panel midpoints on an ellipse (axes 1 and 0.5 m)
e_cir = [ch06.panel_cp_error(N, body="circle") for N in Ns]   # the same on a circle
print("ellipse:", np.array(e_ell), f"order {observed_order(1/Ns, e_ell):.2f}")   # falls ≈ ×4 per doubling: order ≈ 2
print("circle: ", np.array(e_cir))                       # round-off at every N: exact
""")
nb.figure(r"""
fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.4))       # airship | target and strengths | convergence
n = 201 if not FAST else 141                             # grid points per side
Zg, Rg = np.meshgrid(np.linspace(-1.0, 2.0, n), np.linspace(0.005, 1.0, n))   # upper meridian half-plane [m]
Pa = A_s["psi_fn"](Rg, Zg)                               # airship ψ in the body frame [m³/s]
axs[0].contour(Zg, Rg, Pa, levels=np.linspace(-0.3, 0.3, 13), colors=COLORS["ink"], linewidths=0.7)   # streamlines
axs[0].contour(Zg, Rg, Pa, levels=[0.0], colors="k", linewidths=2.5)   # the body ψ = 0
Pf = A_s["psi_fluid_frame"](Rg, Zg)                      # the same body seen from still fluid (stream subtracted)
axs[0].contour(Zg, -Rg, Pf, levels=np.linspace(-0.08, 0.08, 13), colors=COLORS["teal"], linewidths=0.7)   # mirrored below
axs[0].contour(Zg, -Rg, Pa, levels=[0.0], colors="k", linewidths=2.5)
axs[0].plot([0, 1], [0, 0], color=COLORS["rose"], lw=3); axs[0].plot([0], [0], "o", color=COLORS["teal"])   # sink line, source
axs[0].set_aspect("equal"); axs[0].set_xlabel("z [m]"); axs[0].set_ylabel("R [m]")
axs[0].set_title("airship (6.95): body frame above, fluid frame below", fontsize=9)
N = 20                                                   # a clean case for the bars
zb, Rb = ch06.axisym_body_target("rankine_oval", N=N); sol = pn.axial_singularity_solve(zb, Rb, U=1.0, N=N)
axs[1].plot(zb, Rb, "--", color=COLORS["muted"], lw=2, label="target body")   # the target
Zc, Rc = np.meshgrid(np.linspace(zb.min() - 0.3, zb.max() + 0.3, n), np.linspace(0.005, 0.8, n))   # grid for the computed body
axs[1].contour(Zc, Rc, pn.axial_singularity_psi(sol, Rc, Zc), levels=[0.0], colors="k", linewidths=1.5)   # computed ψ = 0
xm = 0.5*(sol["z_nodes"][1:] + sol["z_nodes"][:-1])      # segment midpoints [m]
kk = sol["k"]/np.max(np.abs(sol["k"]))*0.3               # strengths, scaled for display
axs[1].bar(xm, kk, width=np.diff(sol["z_nodes"])*0.9, bottom=-0.4, color=[COLORS["teal"] if v > 0 else COLORS["rose"] for v in kk])
axs[1].set_xlabel("z [m]"); axs[1].set_ylabel("R [m]  (bars: k_n, scaled)"); axs[1].legend(fontsize=7)
axs[1].set_title("Rankine-oval target, N = 20 (bars alternate!)", fontsize=9)
NN = np.array([4, 6, 8, 10, 14, 20, 28, 40])             # numbers of segments
st = [ch06.axial_state("rankine_oval", N=int(k_)) for k_ in NN]   # error and conditioning for each
axs[2].loglog(NN, [s["body_error"] for s in st], "o-", color=COLORS["accent"], label="body error")
ax2b = axs[2].twinx(); ax2b.loglog(NN, [s["cond"] for s in st], "s--", color=COLORS["rose"], label="condition number")
axs[2].set_xlabel("number of segments N"); axs[2].set_ylabel("body error [m]"); ax2b.set_ylabel("cond(A) [–]")
axs[2].set_xticks([4, 10, 20, 40]); axs[2].set_xticklabels(["4", "10", "20", "40"]); axs[2].minorticks_off()   # plain tick labels
axs[2].set_title("more segments: better shape, touchier matrix", fontsize=9)
fig.suptitle("Hidden sources draw the body", fontweight="bold")
savefig(fig, "ch06", "c14_axial"); plt.show()            # save, then draw
""", see="Left, a streamlined hull round a teal source and a rose line sink, and below it the same body seen from still "
         "fluid: loops from nose to tail. Centre, the computed ψ = 0 surface on top of the dashed target, and bars that "
         "alternate teal/rose. Right, the body error falling while the condition number climbs.",
    read="Sources push the flow out at the nose, sinks let it close at the tail; the bars sum to zero. The bars themselves "
         "zig-zag: many strength patterns draw almost the same surface, so trust the shape and the sums, not the "
         "individual bars. More segments fit the shape better but make the system touchier.",
    change="…a blunter target (the sphere): the error plateaus and the condition number explodes — its exact source is a "
           "point doublet, which a smooth line of sources cannot copy.")
nb.plotly(r"""
Nv = np.arange(4, 42, 2) if not FAST else np.arange(4, 42, 4)   # segment numbers, capped at 40 (beyond it the matrix is too ill-conditioned)
cache = {}                                               # solve each N once
def axial_traces(Nf):                                    # target, computed body and scaled strengths for N segments
    N = int(round(Nf))
    zb, Rb = ch06.axisym_body_target("rankine_oval", N=N)   # target points
    sol = pn.axial_singularity_solve(zb, Rb, U=1.0, N=N)    # strengths
    zz = np.linspace(zb.min() - 0.05, zb.max() + 0.05, 200)   # z along the body [m]
    Rgrid = np.linspace(0.005, 0.8, 200)                 # trial radii [m]
    Pz = pn.axial_singularity_psi(sol, Rgrid[:, None], zz[None, :])   # ψ on a (R, z) grid
    Rsurf = np.array([Rgrid[np.argmax(Pz[:, j] > 0)] if np.any(Pz[:, j] > 0) else np.nan for j in range(zz.size)])   # first R with ψ > 0
    xn = sol["z_nodes"]; kstep = np.repeat(sol["k"], 2)/np.max(np.abs(sol["k"]))*0.3 - 0.4   # strengths as a step line
    cache[N] = ch06.axial_state("rankine_oval", N=N)
    return {"target": (zb, Rb), "computed ψ = 0": (zz, Rsurf), "k_n (scaled)": (np.repeat(xn, 2)[1:-1], kstep)}
fig = slider_figure(axial_traces, "N", Nv, xlabel="z [m]", ylabel="R [m]", title="The inverse method converging",
                    yrange=[-0.8, 0.8])
recolor(fig, {"target": COLORS["muted"], "computed ψ = 0": COLORS["ink"], "k_n (scaled)": COLORS["teal"]},
        dashes={"target": "dash"})
step_titles(fig, [f"The inverse method converging: N = {int(v)}, body error {cache[int(v)]['body_error']:.1e} m, "
                  f"cond {cache[int(v)]['cond']:.1e}" for v in Nv])
fig.show()                                               # draw it (works on the web page too)
""")
see_read_change(
    see="The computed surface (black) settles onto the dashed target as N grows, while the teal step line of strengths "
        "zig-zags more and more.",
    read="The shape converges; the individual strengths do not — only their net and their moments do. The title's condition "
         "number says how much the solution amplifies round-off.",
    change="…N pushed beyond 40: the condition number passes 10¹² and the bars become noise (the explainer caps N there).")
explainer("axial_singularity_bodies", "Given a shape, which sources draw it?",
          "convergence is a sequence, not a picture: raise N and watch the computed ψ = 0 surface snap onto the target "
          "while the k_n bars change and the condition number climbs; switch to panels to see the same idea on a 2-D "
          "body's surface.", "", [
              "Preset 'airship (6.95)': compare the fitted bars with the exact point source and flat line sink — the "
              "surfaces agree, the bars need not.",
              "Target 'Rankine oval', N from 4 to 40: the error falls, cond rises, the bars alternate in sign.",
              "Target 'sphere': the ⚠️ status shows why a blunt body is hard.",
              "Mode 'panels': on a circle C_p is exact at every N; on the ellipse it converges like 1/N².",
          ])
whatif(r"""
…we used vortices on the axis or on the surface instead of sources? Sources give thickness but no lift; vortex panels add
circulation — the lifting version is how Ch. 14's panel codes compute an airfoil's lift (with the Kutta condition fixing
the total Γ).
""")

# =====================================================================================================================
# A.9 §6.9 — R31 R32, C15 (N91–N105, D28–D31, E9)
# =====================================================================================================================
nb.section("6.9", "Three-Dimensional Potential Flow and Apparent Mass", intro=r"""
**What is this section about?** Steady ideal flow gives no drag in 3-D either. But when a body *accelerates*, the fluid
round it must be accelerated too, and the pressure field that does this pushes back on the body. For a sphere the push
is exactly that of an extra half of the displaced fluid's mass: the added (or apparent) mass.
""")
nb.recap("R31", "The 3-D potential", r"""
When $\boldsymbol\omega=0$ *(Eq. 3.17)*, u = ∇φ with all three components, so the z-component is w ≡ ∂φ/∂z. ⚠️ In this
section w is a velocity component, not the complex potential of §6.4.
""", where="Ch. 3 §3.4")
nb.recap("R32", "Unsteady Bernoulli between the surface and far away", r"""
For unsteady irrotational flow of constant density with the fluid at rest and at pressure p∞ far away:
$\Big[\frac{\partial\phi}{\partial t}+\frac12\lvert\nabla\phi\rvert^2+\frac p\rho\Big]_{\text{sphere's surface}}=\frac{p_\infty}{\rho}$
*(Eq. 6.99, from 4.75)*. ⚠️ ∂φ/∂t is taken at a *fixed* point in space. Ch. 4's `unsteady_bernoulli_pressure` returns p
from ∂φ/∂t and the speed.
""", where="Ch. 4 §4.9")
nb.code(r"""
print(be.unsteady_bernoulli_pressure(dphi_dt=-0.05, speed=1.0, rho=1000.0, g=0.0, C=0.0))   # p = −ρ(∂φ/∂t + ½|u|²) = −450 Pa
""")

core("C15", "The accelerating sphere: surface pressure and added mass", r"""
An ideal fluid exerts no drag on a steadily moving body. So why is it hard to shake a ball under water — and why does a
released bubble not shoot up infinitely fast?
""", eqs=("6.109",))
problem(r"""
Push a beach ball under water and let go: it accelerates up, but far less violently than its tiny mass would suggest. A
bubble has almost no mass, yet it does not leave with infinite acceleration. A submarine, a fish, a kite and a ship in
waves all feel the same effect: to speed up, a body must also speed up the fluid it pushes aside.
""")
idea("""
pressure on a moving sphere = steady part (~ speed^2, symmetric front/back)     --> no net force (d'Alembert)
                            + acceleration part (~ du_s/dt, + in front, - behind) --> force  -M du_s/dt
M = 1/2 x (mass of the displaced fluid) = (2 pi / 3) rho a^3        -->   F_E = (m + M) du_s/dt
""")
note("N91", "3-D d'Alembert and beyond", r"""
A closed body moving steadily through ideal fluid feels no drag in 3-D too (an exercise of S06); forces appear when the
motion is unsteady or vorticity is present (Ch. 14). The inertia of the fluid is what vehicles, fish and bubbles feel:
apparent or added mass.
""")
note("N92", "The set-up", r"""
A sphere of radius a sits at $\mathbf x_s(t)$ and moves with velocity $\mathbf u_s=d\mathbf x_s/dt$ and a known
acceleration $d\mathbf u_s/dt$; the fluid is at rest far away at pressure p∞; an external force $\mathbf F_E$ (a push, or
buoyancy minus weight) acts on the sphere; we want the fluid's force $\mathbf F_s$ on it. Think of a submarine or a fish
manoeuvring.

```
          fluid at rest, p = p_inf far away
                 .-----.
   F_E  ---->   (  x_s  )  ---> u_s(t),  du_s/dt
                 '-----'         xi = x - x_s,  e_xi = xi/|xi|
```
""")
remind([
    ("multivariable chain rule along a path", r"if f depends on t through x(t), $\frac{df}{dt}=\frac{\partial f}{\partial x_i}\frac{dx_i}{dt}$ summed over i (Ch. 3 P91)."),
], lead=r"""Galilean frames are Ch. 3 §3.3 (recap R14); the coordinate-free sphere potential (6.92) is note N85;
**a function of x − x_s:** ∂/∂x_s = −∇; **∇ of u·ξ/|ξ|³:** product rule with ∇(u·ξ) = u and ∇|ξ|⁻³ = −3ξ/|ξ|⁵ (the
gradient of 1/distance, Ch. 5 P140)""")
D("D28", ref="6.97")
note("N93", "The moving dipole (note `N94`)", r"""
With x → x − x_s and no stream, $\phi=-\frac{1}{4\pi\lvert\mathbf x-\mathbf x_s(t)\rvert^3}\mathbf d\cdot(\mathbf x-\mathbf x_s(t))$
*(6.96)*, and with the dipole following the motion, d(t) = 2πa³u_s:
""", equation=r"\phi(\mathbf x,\mathbf x_s,\mathbf u_s)=-\frac{a^3}{2\lvert\mathbf x-\mathbf x_s\rvert^3}\mathbf u_s\cdot(\mathbf x-\mathbf x_s)=-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\mathbf u_s\cdot\boldsymbol\xi", ref="6.97")
nb.code(r"""
xs, us, a = np.zeros(3), np.array([0.0, 0.0, 2.0]), 0.1  # sphere at the origin moving at 2 m/s along z, radius 0.1 m
e = np.array([[0.6, 0.0, 0.8], [0.0, 1.0, 0.0]]).T       # two unit directions on the surface, shape (3, 2)
u_a = pf.moving_sphere_velocity(xs[:, None] + a*e, xs, us, a)   # the fluid velocity there, ∇φ of (6.97)
print(np.sum(u_a*e, axis=0), us @ e)                     # fluid normal velocity = sphere's normal velocity: (1.6, 0.0) twice
""")
note("N95", "The force", r"""
is the pressure integrated over the surface (the 3-D form of $D\mathbf e_x+L\mathbf e_y=-\frac1B\int_{A^*}p\,\mathbf n\,dA$
*(6.55)*; a uniform p∞ gives nothing on a closed surface):
""", equation=EQ["6.98"], ref="6.98")
D("D29", ref="6.105", check_src=r"""
import sympy as sp                                       # symbolic algebra
t, x, y, z = sp.symbols("t x y z", real=True)            # time and a fixed field point
a_ = sp.symbols("a", positive=True)                      # the sphere's radius (positive)
xs_ = sp.Matrix([sp.Rational(1, 2)*t**2, 0, t])          # a concrete motion x_s(t) whose acceleration turns
us_ = xs_.diff(t)                                        # u_s = dx_s/dt
X = sp.Matrix([x, y, z])                                 # the fixed point x
xi = X - xs_                                             # ξ = x − x_s
r = sp.sqrt(xi.dot(xi))                                  # |ξ|
phi = -a_**3/(2*r**3)*us_.dot(xi)                        # (6.97): the moving sphere's potential
dphidt_direct = sp.diff(phi, t)                          # ∂φ/∂t at fixed x, straight from the formula (no chain rule)
grad = sp.Matrix([sp.diff(phi, v) for v in (x, y, z)])   # ∇φ, (6.103)
chain = -grad.dot(us_) - a_**3/(2*r**3)*xi.dot(us_.diff(t))   # steps 2–5, (6.101): −u·u_s − (a³/2|ξ|³) ξ·du_s/dt
pt = {x: 0.3, y: 0.4, z: 1.2, t: 0.7, a_: 0.5}           # a sample point and instant
print("direct − chain rule:", sp.N((dphidt_direct - chain).subs(pt)))   # → 0: steps 2–5 hold
ev = sp.Matrix([sp.Rational(3, 5), 0, sp.Rational(4, 5)])   # a unit vector e_ξ = (0.6, 0, 0.8)
on_surf = {x: xs_[0] + a_*ev[0], y: xs_[1] + a_*ev[1], z: xs_[2] + a_*ev[2]}   # the surface point x = x_s + a e_ξ
g_s = grad.subs(on_surf)                                 # ∇φ on the surface (step 8's left side)
correct = sp.Rational(3, 2)*us_.dot(ev)*ev - us_/2       # (6.104): (3/2)(u_s·e)e − ½u_s
printed = -(a_**3/2)*(-3*ev*us_.dot(ev)/a_**3 - us_/a_**3)   # the book's middle bracket with −u_s/a³
print("∇φ − (6.104):", sp.simplify(g_s - correct).T)     # → 0: the final form is right
print("printed − (6.104):", sp.simplify(printed - correct).T)   # → u_s ≠ 0: the printed middle term has the wrong sign
p_rho = (-dphidt_direct - grad.dot(grad)/2).subs(on_surf)   # (6.100): (p_a − p∞)/ρ from Bernoulli on the surface
c = us_.dot(ev)                                          # u_s·e_ξ
form = sp.Rational(9, 8)*c**2 - sp.Rational(5, 8)*us_.dot(us_) + a_/2*ev.dot(us_.diff(t))   # (6.105), steps 12–13
print("Bernoulli − (6.105):", sp.simplify(p_rho - form)) # → 0 at every t and a
print(ch06.moving_sphere_dphidt_sym()["difference_chain"])   # the library's own sympy route: 0 as well
""")
note("N96", "The chain of D29 (notes `N97`–`N101`)", r"""
Bernoulli on the surface $\frac{p_a-p_\infty}\rho=-\big(\frac{\partial\phi}{\partial t}\big)_a-\frac12\lvert\nabla\phi\rvert^2_a$
*(6.100)*; the chain rule $\frac{\partial\phi}{\partial t}=-\mathbf u\cdot\mathbf u_s-\frac{a^3}{2\lvert\boldsymbol\xi\rvert^3}\boldsymbol\xi\cdot\frac{d\mathbf u_s}{dt}$
*(6.101)*, on the surface $\big(\frac{\partial\phi}{\partial t}\big)_a=-\mathbf u_a\cdot\mathbf u_s-\frac a2\mathbf e_\xi\cdot\frac{d\mathbf u_s}{dt}$
*(6.102)*; the gradient $\nabla\phi=-\frac{a^3}2\big[-\frac{3\boldsymbol\xi}{\lvert\boldsymbol\xi\rvert^5}\mathbf u_s\cdot\boldsymbol\xi+\frac{\mathbf u_s}{\lvert\boldsymbol\xi\rvert^3}\big]$
*(6.103)*, on the surface $\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s$ *(6.104)*; and
the result below. ⚠️ The book's middle expression of (6.104) prints −u_s/a³ inside the bracket; +u_s/a³ is what (6.103)
gives and what the final form needs (the sympy check above: the printed bracket would give +½u_s instead of −½u_s).
Check with your hands: at the sphere's side (e ⟂ u_s) the fluid moves at −½u_s — backwards, filling in.
""", equation=EQ["6.105"], ref="6.105")
nb.code(r"""
d = ch06.moving_sphere_dphidt_sym()                      # sympy re-run of (6.101)–(6.104) by the library
print(d["correct_bracket"].T, d["printed_bracket"].T)    # difference from (6.104): 0 for the correct bracket, u_s for the printed one
""")
note("N102", "Steady motion", r"""
— $C_p=1-\frac94\sin^2\theta$ *(6.91)* again: a sphere moving steadily through still fluid feels the same pressures as a
still sphere in a stream (Galilean invariance), and no drag. ⚠️ θ_s is measured from the sphere's velocity, i.e. from its
front; in (6.91) θ runs from the downstream side — sin² is the same either way.
""", equation=r"\Big(\frac{p_a-p_\infty}{\frac12\rho\lvert\mathbf u_s\rvert^2}\Big)_{\text{steady}}=\frac94\cos^2\theta_s-\frac54=1-\frac94\sin^2\theta_s", ref="6.106")
P("P163", "surface integrals on a sphere", r"""
On a sphere of radius a use the polar angle θ (from the chosen axis) and the azimuth φ: the area element is
$dA=a^2\sin\theta\,d\theta\,d\varphi$ and the outward normal is
$\mathbf e_r=(\sin\theta\cos\varphi,\ \sin\theta\sin\varphi,\ \cos\theta)$. Integrate φ first: anything with cos φ or sin φ
vanishes. A handy result: $\int_0^\pi\cos^2\theta\sin\theta\,d\theta=\big[-\tfrac13\cos^3\theta\big]_0^\pi=\tfrac23$.
`dblquad(f, a, b, c, d)` integrates f(inner, outer) over two variables.
""", code=r"""
from scipy.integrate import dblquad                      # double integrals (Ch. 3 P87)
area = dblquad(lambda ph, th: np.sin(th), 0, np.pi, 0, 2*np.pi)[0]   # ∫∫ sin θ dφ dθ on the unit sphere
I = dblquad(lambda ph, th: np.cos(th)**2*np.sin(th), 0, np.pi, 0, 2*np.pi)[0]   # ∫∫ cos²θ sin θ dφ dθ
print(area, 4*np.pi, I, 4*np.pi/3)                       # 12.566 = 4π and 4.189 = 2π × 2/3
""")
D("D30", ref="6.108")
note("N103", "The force integral (note `N104`)", r"""
With the acceleration along e_z,
$\mathbf F_s=-\rho\frac a2\big\lvert\frac{d\mathbf u_s}{dt}\big\rvert\int_0^\pi\!\int_0^{2\pi}\cos\theta\,\mathbf e_\xi\,a^2\sin\theta\,d\varphi\,d\theta$
*(6.107)*, and after the φ-integration the result below. ⚠️ The book's (6.108) still shows dφ in the integrand after the
φ-integration has produced its 2π. M is **half the mass of the displaced fluid**: a = 0.1 m in water, M = 2.094 kg vs
4.189 kg displaced. In 2-D a circular cylinder's added mass is ρπa² per unit depth — the *whole* displaced mass (an
exercise of S07; `ch06.cylinder_added_mass(0.1, 1000.0)` = 31.42 kg/m). In general M is a tensor,
$(F_s)_i=-M_{ij}\,d(u_s)_j/dt$ (the apparent-mass tensor).
""", equation=r"\mathbf F_s=-\frac23\pi\rho a^3\frac{d\mathbf u_s}{dt}=-M\frac{d\mathbf u_s}{dt},\qquad M=\frac{2\pi a^3\rho}{3}", ref="6.108")
note("N105", "Why the fluid pushes back", r"""
In front of an accelerating sphere the fluid must be pushed aside faster and faster, behind it the hole must be filled
faster and faster; the pressure that does this is high in front and low behind.
""")
P("P164", "Green's first identity", r"""
For any smooth φ, $\nabla\cdot(\phi\nabla\phi)=\lvert\nabla\phi\rvert^2+\phi\nabla^2\phi$ (the product rule for a divergence,
Ch. 4 P113). If ∇²φ = 0 the last term drops, so the kinetic-energy density ½ρ|∇φ|² is a divergence — and Gauss turns its
volume integral into a surface integral (the divergence theorem, Ch. 2 §2.12). For the fluid outside a body the surface
is the body (with the normal pointing **into** the body, out of the fluid) plus a far sphere.
""", code=r"""
x, y, z = sp.symbols("x y z")                            # Cartesian symbols
phi_ = x*y - z**2 + x**2/2 + y**2/2                      # a harmonic polynomial: ∇²φ = 1 + 1 − 2 = 0
grad_ = [sp.diff(phi_, v) for v in (x, y, z)]            # ∇φ
lhs = sum(sp.diff(phi_*g, v) for g, v in zip(grad_, (x, y, z)))   # ∇·(φ∇φ)
print(sp.simplify(lhs - sum(g**2 for g in grad_)))       # 0: it equals |∇φ|² because ∇²φ = 0
""")
nb.md(r"""
> 🔁 **Gloss — kinetic energy and added mass.** A body of mass m moving at U has kinetic energy ½mU²; the fluid it sets
> moving carries ½MU² more. Pushing both up to speed costs ½(m + M)U² — the same M as the force route (D31 shows it).
""")
D("D31", ref="6.108", check_src=r"""
import sympy as sp                                       # symbolic algebra
r, th, ph, a, U, rho = sp.symbols("r theta phi a U rho", positive=True)   # spherical coordinates and constants
phi = -a**3*U*sp.cos(th)/(2*r**2)                        # (6.97) with u_s = U e_z, sphere at the origin
gr = sp.diff(phi, r)                                     # ∂φ/∂r: radial component of ∇φ
gt = sp.diff(phi, th)/r                                  # (1/r)∂φ/∂θ: polar component of ∇φ
print("∇²φ:", sp.simplify(sp.diff(r**2*sp.diff(phi, r), r)/r**2 + sp.diff(sp.sin(th)*sp.diff(phi, th), th)/(r**2*sp.sin(th))))   # → 0 (step 2)
vol = sp.integrate(sp.integrate((gr**2 + gt**2)*r**2*sp.sin(th), (th, 0, sp.pi)), (r, a, sp.oo))*2*sp.pi   # ∫|∇φ|² dV over r > a
surf = -sp.integrate((phi*gr).subs(r, a)*a**2*sp.sin(th), (th, 0, sp.pi))*2*sp.pi   # −∮ φ ∂φ/∂r dA on r = a (steps 5–7)
print("volume − surface:", sp.simplify(vol - surf), "| volume:", sp.simplify(vol))   # → 0 and 2πU²a³/3
M = sp.simplify(rho*vol/U**2)                            # ½MU² = ½ρ∫|∇φ|² dV  ⇒  M = ρ∫|∇φ|² dV / U²
print("M =", M)                                          # → 2πa³ρ/3, the force route's (6.108)
far = sp.simplify((phi*gr).subs(r, sp.Symbol("R", positive=True))*4*sp.pi*sp.Symbol("R", positive=True)**2)   # ~ far-surface size
print("far-surface term ~", far)                         # ∝ R⁻³ → 0 as R → ∞ (step 4)
""")
nb.worked_example("a ball, a sphere of water and a bubble", r"""
a = 0.1 m, water ρ = 1000 kg/m³.

1. Displaced mass ρ(4/3)πa³ = 4.189 kg; added mass M = 2.094 kg.
2. A steel ball (7800 kg/m³): m = 32.67 kg, behaves as m + M = 34.77 kg; released from rest it accelerates (downward) at
   $(m-\rho V)g/(m+M)=(32.67-4.19)\times9.81/34.77=8.04$ m/s² — without added mass it would be 8.55 m/s².
3. A bubble (m ≈ 0): the **buoyancy** ρVg = 41.1 N (the weight of the displaced water, Ch. 1 §1.7) accelerates only M:
   41.1/2.094 = 19.6 m/s² = **2g**, not infinity.
""")
remind([
    ("Gauss–Legendre quadrature", "N well-chosen nodes and weights integrate polynomials of degree 2N − 1 exactly; `np.polynomial.legendre.leggauss(N)` gives them on [−1, 1] (Ch. 5 P143)."),
    ("scipy.integrate.solve_ivp", "integrates dy/dt = f(t, y) from an initial value; `sphere_motion` uses it for (6.109) (Ch. 1 P31)."),
])
nb.code(r"""
a, rho = 0.1, 1000.0                                     # radius [m], water [kg/m³]
print(pf.added_mass_sphere(a, rho), ch06.added_mass_by_energy(a, rho), ch06.added_mass_by_energy(a, rho, method="volume"))   # M three ways [kg]
us = np.array([0.0, 0.0, 1.0]); dus = np.array([0.0, 0.0, 2.0])   # velocity 1 m/s and acceleration 2 m/s², both along z
F = pf.sphere_force_quadrature(lambda e: pf.moving_sphere_surface_pressure(e, us, dus, a, rho), a)   # −∮(p − p∞) n dA [N]
print(F, -pf.added_mass_sphere(a, rho)*dus)              # (0, 0, −4.1888) N twice: F_s = −M du_s/dt
print(pf.moving_sphere_surface_pressure(np.array([0, 0, 1.0]), us, dus, a, rho, split=True))   # at the front: 500 + 100 = 600 Pa
print(ch04.accelerating_sphere_pressure(0.0, a, 2.0, rho, 1.0)["p"])   # Ch. 4's accelerating sphere at its front: 600 Pa too
print(ch06.sphere_motion(0.0, a, rho, g=9.81, mode="bubble", t_eval=np.array([0.0, 0.01]))["du_dt"][0]/9.81)   # bubble: 2.0 g
""", explain=r"""
1. The added mass by the pressure force (6.108), by the kinetic energy on the surface (D31 steps 5–7) and by a `dblquad`
   volume integral of ½ρ|∇φ|²: 2.0944 kg all three.
2. `sphere_force_quadrature` integrates the surface pressure (6.105) over the sphere (Gauss–Legendre in cos θ ×
   trapezoid in φ): the force is −M du_s/dt; the steady part integrated to zero.
3. At the front point: the steady part ½ρ|u_s|²(9/4 − 5/4) = 500 Pa and the acceleration part ρ(a/2)·2 = 100 Pa.
   Ch. 4's accelerating-sphere tool gives the same 600 Pa (parity).
4. `sphere_motion` integrates $\mathbf F_E=(m+M)\,d\mathbf u_s/dt$ *(6.109)*; a massless bubble starts at 2g.
""")
nb.md("**From scratch — our own surface sum of −(p − p∞) n dA** (🔁 `np.meshgrid(mu, ph, indexing=\"ij\")` makes the "
      "first index run along the first argument, so `Mu[i, k]` pairs node i with azimuth k — Ch. 2 P76):")
nb.check_agree(r"""
mu, w = np.polynomial.legendre.leggauss(24)              # 24 Gauss–Legendre nodes μ = cos θ in [−1, 1] and their weights
ph = np.linspace(0, 2*np.pi, 48, endpoint=False)         # 48 azimuths [rad] (periodic trapezoid)
Mu, Ph = np.meshgrid(mu, ph, indexing="ij")              # all (μ, φ) pairs
st = np.sqrt(1 - Mu**2)                                  # sin θ
E = np.array([st*np.cos(Ph), st*np.sin(Ph), Mu])         # outward unit normals, shape (3, 24, 48)
p = pf.moving_sphere_surface_pressure(E.reshape(3, -1), us, dus, a, rho).reshape(Mu.shape)   # p − p∞ from (6.105) [Pa]
Fz = -np.sum(p*E[2]*w[:, None])*a**2*(2*np.pi/48)        # z-force: dA = a² dμ dφ [N]
print(Fz)                                                # −4.18879 N
assert np.isclose(Fz, -pf.added_mass_sphere(a, rho)*2.0) # = −M du_s/dt
""")
nb.md("Our own surface sum is −M du_s/dt: the steady part cancelled, the acceleration part did not.")
nb.figure(r"""
ths = np.radians(np.linspace(0, 180, 181))               # angle from the direction of motion θ_s [rad]
E_ = np.array([np.sin(ths), 0*ths, np.cos(ths)])         # surface normals in the (x, z) plane
parts = pf.moving_sphere_surface_pressure(E_, us, dus, a, rho, split=True)   # steady, acceleration and total [Pa]
T = np.linspace(0, 1.0, 101)                             # time after release [s]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.8))  # pressure | motion
a1.plot(np.degrees(ths), parts["steady"], color=COLORS["teal"], lw=2, label="steady part (speed)")   # symmetric about 90°
a1.plot(np.degrees(ths), parts["acceleration"], color=COLORS["orange"], lw=2, label="acceleration part")   # + in front, − behind
a1.plot(np.degrees(ths), parts["total"], color="k", lw=2.5, label="total")   # (6.105)
a1.axhline(0, color=COLORS["muted"], lw=0.8); a1.set_xlabel("θ_s, angle from the direction of motion [deg]")   # p = p∞ and the x-axis
a1.set_ylabel("p − p∞ [Pa]"); a1.legend(fontsize=8); a1.set_title("u_s = 1 m/s, du_s/dt = 2 m/s², a = 0.1 m, water", fontsize=9)   # the case
for mass, lab, col in ((0.0, "bubble", COLORS["blue"]), (32.67, "steel ball", COLORS["ink"])):   # two bodies released from rest
    for am, ls_ in ((True, "-"), (False, "--")):         # with and without added mass
        if mass == 0.0 and not am:   # a massless bubble without added mass …
            continue                                     # a massless bubble without added mass would accelerate infinitely
        sm = ch06.sphere_motion(mass, a, rho, g=9.81, mode="bubble" if mass == 0.0 else "ball", t_eval=T, added_mass=am)   # (6.109) in time
        a2.plot(T, np.abs(sm["u"]), ls=ls_, color=col, lw=2, label=f"{lab}{'' if am else ' (no added mass)'}")   # speed [m/s]
a2.set_xlabel("t [s]"); a2.set_ylabel("speed |u| [m/s]"); a2.legend(fontsize=8)   # axes with units
a2.set_title("released from rest", fontsize=9)   # the case
fig.suptitle("Speed pushes symmetrically; acceleration pushes back", fontweight="bold")   # the message
savefig(fig, "ch06", "c15_added_mass"); plt.show()       # save, then draw
""", see="Left, a teal curve symmetric about 90° (500 Pa at both ends, −625 Pa at the side), an orange curve that is +100 Pa "
         "in front and −100 Pa behind, and their black sum. Right, speeds after release growing linearly: the bubble at 2g "
         "= 19.6 m/s², the steel ball at 8.04 m/s² with added mass and a little faster without (dashed).",
    read="Only the orange part survives the force integral: F = −M du_s/dt. The lighter the body, the more its motion is "
         "ruled by M (a bubble's initial acceleration 2g instead of infinity; a steel ball changes by only 6 %).",
    change="…the body a long cylinder moving sideways: its added mass equals the whole displaced mass per unit length — "
           "a massless cylinder would start at g, not 2g.")
nb.animation(r"""
nfr = 60 if not FAST else 30                             # number of frames
om, A0 = 2.0, 0.05                                       # angular frequency [rad/s] and amplitude [m] of x_s = A0 sin ωt
tt = np.linspace(0, 2*np.pi/om, nfr)                     # one period [s]
sm = ch06.sphere_motion(1.0, a, rho, mode="oscillating", amplitude=A0, omega=om, t_eval=tt)   # prescribed motion, F_s(t)
Mv = pf.added_mass_sphere(a, rho)                        # added mass [kg]
thc = np.linspace(0, 2*np.pi, 181)                       # angles round the section [rad]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6))     # sphere section | force trace
sc = a1.scatter(0.1*np.sin(thc), 0.1*np.cos(thc), c=0*thc, cmap="coolwarm", vmin=-1.2, vmax=1.2, s=30)   # surface pressure
arr = a1.quiver([0], [0], [0], [0], color=COLORS["orange"], angles="xy", scale_units="xy", scale=2.0)    # the fluid force
a1.set_xlim(-0.25, 0.25); a1.set_ylim(-0.25, 0.25); a1.set_aspect("equal"); a1.set_xlabel("x [m]"); a1.set_ylabel("z [m]")
fig.colorbar(sc, ax=a1, label="p − p∞ [Pa]")
a2.plot(tt, sm["F_s"], color="k", lw=2, label="F_s from the pressure")
a2.plot(tt, -Mv*sm["du_dt"], "--", color=COLORS["orange"], lw=2, label="−M du_s/dt")
dot, = a2.plot([tt[0]], [sm["F_s"][0]], "o", color=COLORS["accent"])   # the current instant
a2.set_xlabel("t [s]"); a2.set_ylabel("force along z [N]"); a2.legend(fontsize=8)
def update(i):                                           # frame i
    e_ = np.array([np.sin(thc), 0*thc, np.cos(thc)])     # surface normals of the section
    pv = pf.moving_sphere_surface_pressure(e_, np.array([0, 0, sm["u"][i]]), np.array([0, 0, sm["du_dt"][i]]), a, rho)   # p − p∞
    sc.set_offsets(np.c_[0.1*np.sin(thc), sm["x"][i] + 0.1*np.cos(thc)]); sc.set_array(pv)   # move and recolour the sphere
    arr.set_offsets([[0, sm["x"][i]]]); arr.set_UVC([0], [sm["F_s"][i]])   # the force arrow
    dot.set_data([tt[i]], [sm["F_s"][i]])                # the dot on the trace
    a1.set_title(f"t = {tt[i]:.2f} s", fontsize=10)
show_animation(animate(update, frames=nfr, fig=fig, interval=60), player="video")   # smooth MP4
""")
see_read_change(
    see="The sphere oscillating up and down with its surface colour flickering; the orange force arrow always points "
        "against the acceleration; on the right, the pressure force (black) and −M du_s/dt (dashed orange) lie on top of "
        "each other.",
    read="The force is in phase with the acceleration, never with the speed: over a period it does no net work (no drag), "
         "but the sphere is harder to shake than its own mass suggests.",
    change="…ω doubled: the accelerations, the pressures and the force grow four times (∝ ω²) at the same amplitude.")
nb.plotly(r"""
T = np.linspace(0, 0.5, 60)                              # time after release [s]
def release_traces(ratio):                               # a sphere of density ratio ρ_body/ρ released from rest in water
    mass = ratio*rho*4/3*np.pi*a**3                      # its mass [kg]
    with_am = ch06.sphere_motion(mass, a, rho, g=9.81, mode="bubble" if ratio < 1 else "ball", t_eval=T)   # (6.109)
    tr = {"with added mass": (T, np.abs(with_am["u"]))}   # speed [m/s]
    if ratio > 0:                                        # without added mass (undefined for a massless body)
        no_am = ch06.sphere_motion(mass, a, rho, g=9.81, mode="bubble" if ratio < 1 else "ball", t_eval=T, added_mass=False)   # M left out
        tr["without added mass"] = (T, np.abs(no_am["u"]))   # its speed [m/s]
    else:
        tr["without added mass"] = (T, np.full(T.size, np.nan))   # nothing to draw for a massless body
    return tr   # {trace name: (x, y)}
ratios = np.linspace(0, 8, 17 if not FAST else 9)        # ρ_body/ρ from a bubble (0) to steel (≈ 8)
fig = slider_figure(release_traces, "ρ_body/ρ", ratios, xlabel="t [s]", ylabel="speed |u| [m/s]",   # precompute every slider step
                    title="When added mass matters: light bodies", yrange=[0, 10])
recolor(fig, {"with added mass": COLORS["accent"], "without added mass": COLORS["muted"]}, dashes={"without added mass": "dash"})   # colours by meaning
g0 = [abs(ch06.sphere_motion(r_*rho*4/3*np.pi*a**3, a, rho, g=9.81, mode="bubble" if r_ < 1 else "ball",
                             t_eval=np.array([0.0, 0.01]))["du_dt"][0])/9.81 for r_ in ratios]   # initial acceleration / g
step_titles(fig, [f"When added mass matters: initial acceleration {g_:.2f} g (ρ_body/ρ = {r_:.1f})" for g_, r_ in zip(g0, ratios)])   # a title per step
fig.show()                                               # draw it (works on the web page too)
""")
see_read_change(
    see="Two straight speed lines after release; for a bubble only the solid one exists (2g); for light bodies the gap to "
        "the dashed 'no added mass' line is large, for steel it almost vanishes.",
    read=r"The acceleration is $\lvert m-\rho V\rvert g/(m+M)$: M matters most when m is small; at ρ_body = ρ nothing moves "
         "(neutral buoyancy).",
    change="…the fluid were air instead of water: M would be 800 times smaller and irrelevant for everything but a balloon.")
remind([
    ("live widgets", "`live(fn, name=(min, max, step), …)` gives sliders that re-run fn — only while a kernel runs (Ch. 1 P47)."),
])
nb.live(r"""
def push_sphere(m=10.0, a=0.1, F_E=10.0):                # mass [kg], radius [m], constant push [N]
    tq = np.linspace(0, 2, 101)                          # time [s]
    w1 = ch06.sphere_motion(m, a, 1000.0, F_E=lambda t: F_E, t_eval=tq)   # (6.109) with added mass, water
    w0 = ch06.sphere_motion(m, a, 1000.0, F_E=lambda t: F_E, t_eval=tq, added_mass=False)   # without
    fig, ax = plt.subplots(figsize=(6, 3.2))             # one panel
    ax.plot(tq, w1["u"], color=COLORS["accent"], lw=2, label=f"with added mass: m + M = {m + pf.added_mass_sphere(a, 1000.0):.2f} kg")
    ax.plot(tq, w0["u"], "--", color=COLORS["muted"], lw=2, label="without")
    ax.set_xlabel("t [s]"); ax.set_ylabel("u [m/s]"); ax.legend(fontsize=8); plt.show()
live(push_sphere, m=(0.5, 40.0, 0.5), a=(0.02, 0.3, 0.01), F_E=(1.0, 50.0, 1.0))   # three sliders (the plotly figure above is the page version)
""")
explainer("added_mass_sphere", "Why does an ideal fluid resist acceleration?",
          "set speed and acceleration independently, even at an angle, and see the two pressure parts on the sphere with "
          "their separate forces — one always zero, the other always −M du_s/dt; the ball and bubble modes integrate "
          "(6.109) in time.", "", [
              "Preset 'steady motion': the teal pressure is strong, the force bar zero.",
              "Preset 'starting from rest': only the orange part, and the force equals −M du_s/dt.",
              "Turn the acceleration 90° away from the velocity: the force follows the acceleration, not the motion.",
              "Preset 'air bubble': the transport shows it starting at 2g.",
          ])
whatif(r"""
…the body were not a sphere? M becomes a tensor: a flat plate accelerated broadside drags along far more fluid than
edgewise — why a kite or a sail reacts differently to gusts from different directions, and why fish shape their bodies as
they do (Ch. 16).
""")

# =====================================================================================================================
# A.10 §6.10 — N106, S01–S07, summary
# =====================================================================================================================
nb.section("6.10", "Concluding Remarks", intro=r"""
**What is this section about?** Where ideal-flow theory came from, where it fails, and where it lives on.
""")
nb.md(r"""
📝 **Note.** **Two and a half centuries** `N106` — Euler, Bernoulli, d'Alembert, Lagrange, Stokes, Helmholtz, Kirchhoff and
Kelvin built potential-flow theory; the same mathematics describes heat conduction, elasticity and electrostatics. Its
zero-drag prediction contradicted every observation until Prandtl's boundary layer (Ch. 9) explained where viscosity
hides. It remains the working tool for pressure on streamlined bodies (Ch. 9) and lift (Ch. 14, with conformal maps and
the Kutta condition). It closes this chapter's force story: d'Alembert (C06) → lift (C07, C10) → added mass (C15).
""")
nb.pointer("S01 — Exercises 6.1–6.5, 6.8 and 6.33 check the §6.2 statements (delta-source integrals, orthogonality of ∇φ "
           "and ∇ψ, Bernoulli from the ideal-flow equations); their ideas are D03 and D04 above.")
nb.pointer("S02 — Exercises 6.6–6.7, 6.9, 6.14–6.15, 6.17 and 6.22: harmonic polynomials (generated by "
           "`ch06.harmonic_polynomials`, N20), the doublet's stream function (N25, N47), sketched flows and a vortex-pair "
           "doublet — left to the reader.")
nb.pointer("S03 — Exercises 6.10–6.13, 6.18–6.20 and 6.23–6.24: forces on held singularities, the half-body's width by a "
           "mass balance (D07 step 9) and its zero net force (N29), the Rankine oval (the superposition explainer's preset) "
           "— the rest left to the reader.")
nb.pointer("S04 — Exercises 6.21, 6.25–6.27 and 6.29–6.32: Blasius for the cylinder (the Blasius explainer's cylinder "
           "preset) and a control-volume proof of Kutta–Zhukhovsky (`ch06.cv_force_on_body`, R18); corners, images and "
           "vortex pairs left to the reader.")
nb.pointer("S05 — Exercises 6.16 and 6.41 are kitchen experiments (a paper cylinder vs an airfoil in a draught; a vacuum "
           "nozzle over sugar grains — a sink and its image). Try them.")
nb.pointer("S06 — Exercises 6.34–6.40 and 6.42–6.43: the 3-D source and doublet (D25), 3-D d'Alembert (N91), the "
           "airship's length (D26 and the tiny example of C14) — the rest left to the reader.")
nb.pointer("S07 — Exercises 6.44–6.49: cavity collapse, the cylinder's added mass (N104), a bubble's 2g start (the added-"
           "mass explainer's bubble preset), oscillation frequencies, and the kinetic-energy route to added mass (D31).")
nb.summary(
    clicked=[
        r"**C01** ω = 0 makes the viscous force vanish, so outside boundary layers and wakes the flow obeys $\nabla\cdot\mathbf u=0$ and $\rho D\mathbf u/Dt=-\nabla p$ (6.1) — at the price of the no-slip condition.",
        r"**C02** $\omega_z=-\nabla^2\psi$ (6.4): irrotational means Laplace for ψ and φ, with vortices and sources as delta sources of strength −Γ and m.",
        r"**C03** Laplace is linear, so flows add, and any streamline can be a wall: $\partial\psi/\partial s=0$ (6.16).",
        r"**C04** a source and a sink merged at fixed 2mε give the doublet $\phi=\lvert\mathbf d\rvert\cos\theta/2\pi r$ (6.29), the error shrinking like ε².",
        r"**C05** stream + source = half-body: stagnation at m/2πU, body ψ = m/2, width m/U far downstream (6.31).",
        r"**C06** the cylinder's $C_p=1-4\sin^2\theta$ (6.35) is symmetric front–back and top–bottom: zero drag, d'Alembert's paradox, which separation resolves.",
        r"**C07** circulation moves the stagnation points ($\sin\theta=-\Gamma/4\pi aU$, (6.38)) and gives $L=\rho U\Gamma$ (6.40), with D = 0.",
        r"**C08** images turn walls into symmetry lines (vortex images flip sign, sources keep it); a passing eddy is felt as suction, then a small over-pressure.",
        r"**C09** $w=\phi+i\psi$ (6.42) is analytic — Cauchy–Riemann — and $dw/dz=u-iv$ (6.45); corners $w=Az^n$ are calm (α < π) or violent (α > π).",
        r"**C10** Blasius (6.60) plus residues: D = 0 and $L=\rho U\Gamma$ for any body (6.62) — only U × Γ/z survives far away.",
        r"**C11** analytic maps keep angles (6.64); $z=\zeta+b^2/\zeta$ (6.65) carries the circle onto an ellipse — take the square root that lands outside the circle.",
        r"**C12** on a grid Laplace is \"average of your neighbours\" (6.72); Gauss–Seidel relaxes to it, SOR faster; test the residual with a tight tolerance; a 270° corner limits the accuracy to order 4/3.",
        r"**C13** the Stokes ψ obeys (6.77), not Laplace; the sphere has maximum speed 1.5U and $C_p=1-\frac94\sin^2\theta$ (6.91).",
        r"**C14** unknown axial sources + ψ = 0 at N points = one linear solve for a given body; the shape converges, the individual strengths do not.",
        r"**C15** the acceleration pressure gives $\mathbf F_s=-M\,d\mathbf u_s/dt$ with $M=2\pi\rho a^3/3$ (6.108), half the displaced mass — by force and by energy.",
    ],
    feeds_forward=[
        "velocity potentials of water waves and images at the bottom (Ch. 7)",
        "the ideal outer flow U_e(x) that drives boundary layers, wedge flows w = Azⁿ and separation (Ch. 9)",
        "elliptic solvers (relaxation, SOR, sparse direct) and panel methods (Ch. 10)",
        "ψ–vorticity inversion, deformation fields, coastal images and bounded-basin Poisson problems (Ch. 13)",
        "the Kutta condition, Zhukhovsky airfoils and induced drag (Ch. 14)",
        "added mass of bubbles, fish and structures (Ch. 16)",
    ],
    left_out=[
        "the exercise details (S01–S07 above)",
        "the measured cylinder pressure curve of the book's figure — only a labelled qualitative band here (no dataset cited)",
        "the moment (torque) form of Blasius's theorem",
        "added-mass tensors of non-spherical bodies (named in N104)",
    ],
)

# ---------------------------------------------------------------------------------------------------------------------
# final pass: every equation named by number is written out; save (nbkit's coverage checks run here)
# ---------------------------------------------------------------------------------------------------------------------
n_changed = finalize_equations()
n_changed += finalize_bare_numbers()
bad = self_check_numbers() + self_check_prose()
if bad:
    raise SystemExit("markdown cells cite equation numbers without maths:\n  " + "\n  ".join(bad))
out = nb.save()
print(f"wrote {out.relative_to(ROOT)} ({len(nb.cells)} cells; equations written out in {n_changed} cells)")
