"""Build the Chapter 15 teaching notebook: ``notebooks/ch15_compressible_flow.ipynb``.

Source of truth: ``analysis/ch15_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract —
``fluidpy.ch15_compressible_flow`` imported as ``ch15`` re-exports ``core.gasdynamics`` = ``GD`` and ``core.acoustics`` =
``AC``; ``scripts/ch15_drawings.py`` only draws), Part D (runtime budget and FAST sizes), Part E (prerequisite ledger →
primers P354–P373 and one-line reminders of earlier primers), Part F (the 31 derivations D01–D31, one move per step) and
``analysis/ch15_curation.md`` (IDs, depths, section coverage), with these rulings written into the text:

* the three physics-changing book slips — $u\,du=-dp/\rho=c^2d\rho/\rho$ (15.32), the numerator of (15.47) and the sign of the
  second relation of (15.48) — are taught CORRECTED, each with a "slip #k — the book prints … ; the correct form is …" box,
  the derivation step that justifies it and a finite-difference cell; the other slips are named where they occur;
* traps (true as printed, easy to misread, ``ch15.traps()``) are "Common confusion" callouts; slips and traps kept apart;
* the book's lumped f = F/(p₁A) and q = Q/h₁ are primary; Fanno and Rayleigh names and ratios are labelled "not in the
  book — ours"; every formula of §15.7 is labelled ours;
* angles are computed in radians and shown in degrees; σ, δ are the β, θ of other books; a turn of the Prandtl–Meyer
  function is counted positive for an expansion in code (turn = −δ of the book) and both numbers are printed;
* no book number is reproduced (rule 9): worked numbers use the design's own six input dicts, written literally in the
  imports cell; ``ch15.illustrative_inputs()`` (the module's own set, also ours) is shown once beside them.

**Derivations are read from Part F at build time** (``part_f()`` below): goal, start, plan, tools, assumptions, every
step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word; the four ★★★ derivations
(D13, D17, D21, D30) carry the sympy check cells listed in Part F.

**Equations.** ``EQ`` below is the analyst's page-image transcription as carried by the design (the designer re-read the
12 pages that hold 31 of the 59 numbered equations), in the CORRECTED form where the book prints a slip.

Run:  .venv/Scripts/python.exe notebooks/build_ch15.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch15.py --dump     (prints the parsed Part F derivations only)
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

nb = ChapterNotebook("ch15")
DESIGN = (ROOT / "analysis" / "ch15_design.md").read_text(encoding="utf-8")
_CTRL = sorted({hex(ord(_c)) for _c in DESIGN if ord(_c) < 32 and _c != "\n"})
if _CTRL:                                # a form feed, tab or backspace in the design = a LaTeX backslash eaten by a non-raw string
    raise SystemExit(f"analysis/ch15_design.md contains control characters {_CTRL} — repair the eaten backslash before building")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# «15.28» inside any text passed through X() becomes "$T_0/T=…$ (15.28)".
# ---------------------------------------------------------------------------------------------------------------------
EQ: dict[str, str] = {
    "15.1h": r"C_v=\dfrac{R}{\gamma-1},\ C_p=\dfrac{\gamma R}{\gamma-1}",
    "15.1c": r"c=\sqrt{\gamma RT}=\sqrt{\gamma p/\rho}",
    "15.1s": r"S_2-S_1=C_p\ln\dfrac{T_2}{T_1}-R\ln\dfrac{p_2}{p_1}=C_v\ln\dfrac{T_2}{T_1}-R\ln\dfrac{\rho_2}{\rho_1}",
    "15.2": r"\dfrac{p_2}{p_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma},\ \dfrac{T_2}{T_1}=\left(\dfrac{\rho_2}{\rho_1}\right)^{\gamma-1}=\left(\dfrac{p_2}{p_1}\right)^{(\gamma-1)/\gamma}",
    "15.3": r"R=287\ \mathrm{m^2/(s^2K)},\ C_v=717,\ C_p=1004,\ \gamma=1.40",
    "15.4": r"u_i=U_i+u_i',\ p=p_0+p',\ \rho=\rho_0+\rho',\ T=T_0+T'",
    "15.5": r"p'\cong c^2\rho'",
    "15.6": r"c^2=\gamma p_0/\rho_0=\gamma RT_0",
    "15.7": r"p'/\rho_0c^2\ll1",
    "15.8": r"\dfrac{\partial\rho'}{\partial t}+U_j\dfrac{\partial\rho'}{\partial x_j}+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0",
    "15.9": r"\dfrac{\partial u_i'}{\partial t}+U_j\dfrac{\partial u_i'}{\partial x_j}+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0",
    "15.10": r"\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)p'+\rho_0\dfrac{\partial u_j'}{\partial x_j}=0",
    "15.11": r"\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)u_i'+\dfrac1{\rho_0}\dfrac{\partial p'}{\partial x_i}=0",
    "15.12": r"\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0",
    "15.13": r"\dfrac1{c^2}\dfrac{\partial^2p'}{\partial t^2}-\dfrac{\partial^2p'}{\partial x^2}=0",
    "15.14": r"u_1'(x,t)=-\dfrac1{\rho_0}\displaystyle\int\dfrac{\partial p'}{\partial x}\,dt",
    "15.15": r"p'(x,t)=f(x-ct)+g(x+ct)",
    "15.16": r"u_1'(x,t)=\dfrac1{\rho_0c}\left(f(x-ct)-g(x+ct)\right)",
    "15.17": r"p'(x,t)=f(x-(c+U)t)+g(x+(c-U)t)",
    "15.18": r"\rho uA=\dot m=\text{const.}",
    "15.19": r"\dfrac{d\rho}{\rho}+\dfrac{du}{u}+\dfrac{dA}{A}=0",
    "15.20": r"\displaystyle\int_{A^*}\rho\left(e+\tfrac12u_i^2\right)u_jn_j\,dA=\int_{A^*}u_i\tau_{ij}n_j\,dA-\int_{A^*}q_jn_j\,dA",
    "15.21": r"-\left(e+\tfrac12u^2\right)_1\dot m+\left(e+\tfrac12u^2\right)_2\dot m=(upA)_1-(upA)_2+\dot mQ",
    "15.22": r"h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q",
    "15.23": r"\displaystyle\int_{A^*}\rho u_iu_jn_j\,dA=\int_{A^*}\tau_{ij}n_j\,dA",
    "15.24": r"-\dot mu_1+\dot mu_2=(pA)_1-(pA)_2+F",
    "15.25": r"\dot m\dfrac{du}{dx}=-A\dfrac{dp}{dx}-F_f",
    "15.26": r"u\,du+\dfrac{dp}{\rho}=0",
    "15.27": r"C_pT_0=C_pT+\tfrac12u^2",
    "15.28": r"\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2",
    "15.29": r"\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}",
    "15.30": r"\dfrac{\rho_0}{\rho}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{1/(\gamma-1)}",
    "15.31": r"\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}",
    "15.32": r"u\,du=-dp/\rho=-c^2d\rho/\rho",
    "15.32p": r"u\,du=-dp/\rho=c^2d\rho/\rho",
    "15.33": r"d\rho/\rho=-M^2\,du/u",
    "15.34": r"\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}",
    "15.35": r"\rho_1u_1=\rho_2u_2",
    "15.36": r"p_1-p_2=-\rho_1u_1^2+\rho_2u_2^2",
    "15.37": r"h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2",
    "15.38": r"\dfrac{\gamma}{\gamma-1}\dfrac{p_1}{\rho_1}+\dfrac12u_1^2=\dfrac{\gamma}{\gamma-1}\dfrac{p_2}{\rho_2}+\dfrac12u_2^2",
    "15.39": r"\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)",
    "15.40": r"M_2^2=\dfrac{(\gamma-1)M_1^2+2}{2\gamma M_1^2+1-\gamma}",
    "15.41": r"\dfrac{\rho_2}{\rho_1}=\dfrac{u_1}{u_2}=\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}",
    "15.42": r"\dfrac{T_2}{T_1}=1+\dfrac{2(\gamma-1)}{(\gamma+1)^2}\dfrac{\gamma M_1^2+1}{M_1^2}\left(M_1^2-1\right)",
    "15.43": r"\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]",
    "15.44a": r"\dfrac{S_2-S_1}{C_v}\cong\dfrac{2\gamma(\gamma-1)}{3(\gamma+1)^2}\left(M_1^2-1\right)^3",
    "15.44b": r"\dfrac{S_2-S_1}{C_v}\cong\dfrac{\gamma^2-1}{12\gamma^2}\left(\dfrac{p_2-p_1}{p_1}\right)^3",
    "15.45": r"\rho_1u_1=\rho_2u_2,\ \ p_1+\rho_1u_1^2=p_2+\rho_2u_2^2+p_1f,\ \ h_1+\tfrac12u_1^2+h_1q=h_2+\tfrac12u_2^2",
    "15.46": r"M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}",
    "15.47": r"-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df,\ \ \dfrac{du}{u}=-\dfrac{d\rho}{\rho}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}",
    "15.47p": r"-\dfrac{dp}{p_1}=\dfrac{1-(\gamma-1)M^2}{1-M^2}df",
    "15.48": r"\dfrac{dT}{T_1}=\left(\dfrac{1-\gamma M^2}{1-M^2}\right)dq,\ \ \dfrac{u\,du}{h_1}=\dfrac{(\gamma-1)M^2}{1-M^2}\,dq",
    "15.48p": r"\dfrac{u\,du}{h_1}=\dfrac{(1-\gamma)M^2}{1-M^2}dq",
    "15.49": r"\sin\mu=1/M",
    "15.50": r"\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2\sin^2\sigma-1\right)",
    "15.51": r"\dfrac{\rho_2}{\rho_1}=\dfrac{(\gamma+1)M_1^2\sin^2\sigma}{(\gamma-1)M_1^2\sin^2\sigma+2}=\dfrac{\tan\sigma}{\tan(\sigma-\delta)}",
    "15.52": r"M_2^2\sin^2(\sigma-\delta)=\dfrac{(\gamma-1)M_1^2\sin^2\sigma+2}{2\gamma M_1^2\sin^2\sigma+1-\gamma}",
    "15.53": r"\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}",
    "15.54": r"\tan\delta=2\cot\sigma\,\dfrac{\gamma+1}{2\gamma}\left(\dfrac{p_2-p_1}{p_1}\right)\dfrac{1}{M_1^2\left(\gamma+1-2\sin^2\sigma\right)+2}",
    "15.55": r"\dfrac{p_2-p_1}{p_1}=\dfrac{\gamma M_1^2}{\sqrt{M_1^2-1}}\,\delta",
    "15.56": r"\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}",
    "15.56i": r"\nu(M)=\displaystyle\int_1^M\dfrac{\sqrt{M^2-1}}{1+\frac12(\gamma-1)M^2}\dfrac{dM}{M}",
    "15.57": r"\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\,\delta}{\sqrt{M_\infty^2-1}}",
    "15.58": r"L\cong\dfrac{2\alpha\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}},\ \ D\cong\dfrac{2\alpha^2\gamma M_\infty^2p_\infty b}{\sqrt{M_\infty^2-1}}",
    "15.59": r"C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}",
    # unnumbered equations of the chapter
    "u15.4": r"S_2-S_1=-R\ln\dfrac{p_{02}}{p_{01}}+C_p\ln\dfrac{T_{02}}{T_{01}}",
    "u15.6": r"p_1\left(1+\gamma M_1^2\right)=p_2\left(1+\gamma M_2^2\right)",
    "u15.6b": r"\dfrac{U_1\ln(U_1-U)-U_2\ln(U-U_2)}{U_1-U_2}=\dfrac{\gamma+1}{2\gamma}\,\eta",
    "u15.8": r"A\equiv\dfrac{M_1^2\left(1+\frac12(\gamma-1)M_1^2+q\right)}{\left(1+\gamma M_1^2-f\right)^2}=\dfrac{M_2^2\left(1+\frac12(\gamma-1)M_2^2\right)}{\left(1+\gamma M_2^2\right)^2}",
    "u15.9n": r"M_{n1}=M_1\sin\sigma,\ \ M_{n2}=M_2\sin(\sigma-\delta)",
    "u15.9d": r"d\delta=-\dfrac{\sqrt{M^2-1}}{M}\dfrac{dM}{1+\frac12(\gamma-1)M^2}",
    "u15.9t": r"\delta_1+\nu(M_1)=\delta_2+\nu(M_2)",
    "ex15.1": r"f(x)=\dfrac12\left[F(x)-\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right],\ g(x)=\dfrac12\left[F(x)+\dfrac1c\displaystyle\int_0^xG(\hat x)\,d\hat x\right]",
    # results that are ours (not in the book)
    "o_sonic": r"\dfrac{T^*}{T_0}=\dfrac{2}{\gamma+1},\ \dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)},\ \dfrac{\rho^*}{\rho_0}=\left(\dfrac{2}{\gamma+1}\right)^{1/(\gamma-1)}",
    "o_pstar": r"\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}",
    "o_flux": r"\dfrac{\dot m\sqrt{RT_0}}{p_0A}=\sqrt\gamma\,M\left(1+\dfrac{\gamma-1}{2}M^2\right)^{-(\gamma+1)/(2(\gamma-1))}",
    "o_mmax": r"\dot m_{\max}=p_0A^*\sqrt{\dfrac{\gamma}{RT_0}}\left(\dfrac{2}{\gamma+1}\right)^{(\gamma+1)/(2(\gamma-1))}",
    "o_p02": r"\dfrac{p_{02}}{p_{01}}=\left[\dfrac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2}\right]^{\gamma/(\gamma-1)}\left[\dfrac{\gamma+1}{2\gamma M_1^2-(\gamma-1)}\right]^{1/(\gamma-1)}",
    "o_laval": r"\dfrac{p_b}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}(A_e/A_t)},\ \dfrac{p_f}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sup}}(A_e/A_t)},\ \dfrac{p_d}{p_0}=\dfrac{p_f}{p_0}\left.\dfrac{p_2}{p_1}\right\rvert_{M_{\text{sup}}}",
    "o_shockpos": r"\dfrac{p_B}{p_0}=\dfrac{p_{02}}{p_{01}}(M_1)\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}\left(\frac{A_e}{A_t}\frac{p_{02}}{p_{01}}\right)}",
    "o_fmax": r"f_{\max}=1+\gamma M_1^2-M_1\sqrt{2(\gamma+1)\left(1+\tfrac{\gamma-1}{2}M_1^2\right)}",
    "o_qmax": r"q_{\max}=\dfrac{(1+\gamma M_1^2)^2}{2(\gamma+1)M_1^2}-1-\dfrac{\gamma-1}{2}M_1^2",
    "o_hug": r"\dfrac{p_2}{p_1}=\dfrac{(\gamma+1)\rho_2/\rho_1-(\gamma-1)}{(\gamma+1)-(\gamma-1)\rho_2/\rho_1}",
    "o_delay": r"t_{\text{delay}}=\dfrac{h\sqrt{M^2-1}}{Mc}",
    "o_umax": r"u_{\max}=\sqrt{2C_pT_0}",
    "o_pitot": r"p_0=p+\tfrac12\rho u^2\left(1+\tfrac14M^2+\dots\right)",
    "o_thin": r"C_p=\pm\dfrac{2\theta}{\sqrt{M_\infty^2-1}},\ \ C_D=\dfrac{4}{\sqrt{M_\infty^2-1}}\left[\alpha^2+\overline{\theta_c^2}+\overline{\theta_t^2}\right]",
    "o_steep": r"\left.\dfrac{dx}{dt}\right\rvert_{\text{wavelet}}=u+c=c_0+\dfrac{\gamma+1}{2}u,\ \ t_{\text{break}}=\dfrac{2}{(\gamma+1)\max(-\partial u_0/\partial x)}",
    "o_dMM": r"\dfrac{dM}{M}=-\dfrac{1+\frac{\gamma-1}{2}M^2}{1-M^2}\dfrac{dA}{A}",
    # earlier chapters (written-out forms of their designs)
    "1.19": r"c^2\equiv(\partial p/\partial\rho)_s",
    "4.110": r"\nabla\cdot\mathbf u=-M^2\left(\dfrac{\rho_0}{\rho}\right)\dfrac{D}{Dt}\left(\dfrac{p-p_0}{\rho_0U^2}\right)",
    "4.111": r"M\equiv U/c",
    "7.81": r"H_2/H_1=\tfrac12\left(-1+\sqrt{1+8Fr_1^2}\right)",
    "14.12": r"C_L\simeq2\pi(\alpha+\beta)",
}
LABEL = {
    "15.1h": "15.1, specific-heat line", "15.1c": "15.1, sound-speed line", "15.1s": "15.1, entropy line",
    "15.32": "15.32, third member with its minus sign restored", "15.32p": "15.32, as printed",
    "15.47": "15.47, numerator sign corrected", "15.47p": "15.47, first relation as printed",
    "15.48": "15.48, second relation's sign corrected", "15.48p": "15.48, second relation as printed",
    "15.56i": "15.56, integral form", "u15.4": "unnumbered, §15.4", "u15.6": "unnumbered, §15.6",
    "u15.6b": "unnumbered, §15.6", "u15.8": "unnumbered, §15.8", "u15.9n": "unnumbered, §15.9",
    "u15.9d": "unnumbered, §15.9", "u15.9t": "unnumbered, §15.9", "ex15.1": "unnumbered, Exercise 15.1",
    "15.45": "15.45, with $f=F/(p_1A)$, $q=Q/h_1$", "o_sonic": "ours, from (15.28)–(15.30) at $M=1$",
}
for _k, _v in EQ.items():
    assert "$" not in _v, (_k, _v)
    if _k.startswith("o_"):
        LABEL.setdefault(_k, "ours")


def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ ({LABEL.get(num, num)})"


def X(s: str) -> str:
    """Replace every «key» in a text by the written-out equation and its number (E)."""
    return re.sub(r"«([^»]+)»", lambda m: E(m.group(1)), textwrap.dedent(s).strip())


def _join(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(x.strip() for x in lines)).strip()


def _abs_plain(s: str) -> str:
    r"""Part F writes absolute values as \lvert … \rvert everywhere; outside $…$ they become plain |…|."""
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


# ---------------------------------------------------------------------------------------------------------------------
# Part F of the design, parsed at build time: every derivation copied word for word (the ch07–ch14 pattern; the design
# writes the fields inline — "**Goal:** … · **Start:** … · **Plan:** • … • …" — and the steps as
# "n. **did** … · **tex** … · **why** … · **plain** …"; a ★★★ block ends with "- **sympy check** (…):" and a listing).
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
        body = re.sub(r"^- \*\*sympy check\*\*.*?```python\n.*?```", "", body, flags=re.S | re.M | re.I)
        lines = body.splitlines()
        title = _abs_plain(re.split(r" — ★", lines[0])[0].strip())
        stars = lines[0].count("★")
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
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=_tidy(d["why"]), plain=_tidy(d["plain"])))
        start_tex, start_plain = _start_split(fields["Start"].rstrip(" ·"))
        res_tex, res_plain = _start_split(fields["Result"].rstrip(" ·"))
        plan = [p.strip().rstrip(".·").strip() + "." for p in fields["Plan"].split("•") if p.strip(" ·.")]
        tools = [t.strip().rstrip(".·").strip() for t in fields["Tools"].split(";") if t.strip(" ·.")]
        out[key] = dict(title=title, stars=stars, goal=_tidy(fields["Goal"].rstrip(" ·")),
                        start=(display_tex(start_tex), _tidy(start_plain)), plan=plan, tools=tools,
                        assumptions=_tidy(fields.get("Assumptions", "").rstrip(" ·")), steps=parsed,
                        result=(display_tex(res_tex), _tidy(res_plain)), check=_tidy(fields.get("Check", "")),
                        meaning=_tidy(fields.get("What it means", "")), traps=_tidy(fields.get("Traps", "")), check_src=check_src)
    return out


PF = part_f()
assert len(PF) == 31, len(PF)
N_STEPS = sum(len(d["steps"]) for d in PF.values())
assert N_STEPS == 242, N_STEPS                              # the 242 steps of Part F, none lost by the parser
assert [k for k, d in PF.items() if d["check_src"]] == ["D13", "D17", "D21", "D30"]
assert [len(PF[k]["steps"]) for k in sorted(PF)] == [10, 6, 8, 5, 5, 8, 9, 6, 5, 9, 7, 7, 13, 7, 7, 5, 11, 7, 7, 9, 13,
                                                      9, 8, 7, 4, 8, 8, 8, 9, 10, 7]

# Every line of a check listing carries a comment: Part F leaves a few assert lines bare.
for _key, _old, _new in (
        ("D13", "assert sp.simplify(r_shock - ((g-1) + 2*g*a)/(g+1)) == 0",
         "assert sp.simplify(r_shock - ((g-1) + 2*g*a)/(g+1)) == 0   # ... equal to step 10's (gamma-1)/(gamma+1) + 2 gamma a/(gamma+1)"),
        ("D17", "assert sp.simplify(ser.coeff(e, 1)) == 0 and sp.simplify(ser.coeff(e, 2)) == 0",
         "assert sp.simplify(ser.coeff(e, 1)) == 0 and sp.simplify(ser.coeff(e, 2)) == 0   # sympy's own series: no eps, no eps^2 term"),
        ("D17", "assert sp.simplify(ser.coeff(e, 3) - c3) == 0",
         "assert sp.simplify(ser.coeff(e, 3) - c3) == 0             # and its eps^3 coefficient is the one of step 8"),
        ("D21", "assert sp.simplify(qa - ((g-1)/2 - A*g**2)) == 0 and sp.simplify(qb - (1 - 2*A*g)) == 0 and qc == -A",
         "assert sp.simplify(qa - ((g-1)/2 - A*g**2)) == 0 and sp.simplify(qb - (1 - 2*A*g)) == 0 and qc == -A   # the three coefficients of step 9"),
        ("D21", "    assert sp.simplify(eq.subs(X, X_root)) == 0",
         "    assert sp.simplify(eq.subs(X, X_root)) == 0     # ... put back into step 8, it leaves exactly zero"),
        ("D30", "vals = [abs(sp.N(resid.subs({g: gv, M: Mv}))) for gv in (1.4, sp.Rational(5, 3)) for Mv in (1.2, 2, 5)]",
         "vals = [abs(sp.N(resid.subs({g: gv, M: Mv}))) for gv in (1.4, sp.Rational(5, 3)) for Mv in (1.2, 2, 5)]   # the residual at six (gamma, M) pairs"),
        ("D30", "assert sp.simplify(integrand_t - split) == 0",
         "assert sp.simplify(integrand_t - split) == 0       # the two simple fractions add back to the integrand"),
):
    assert PF[_key]["check_src"].count(_old) == 1, (_key, _old)
    PF[_key]["check_src"] = PF[_key]["check_src"].replace(_old, _new)
# The D13 listing names its quadratic `quad`, which would hide scipy's `quad` (used in C15) for the rest of the session.
PF["D13"]["check_src"] = re.sub(r"\bquad\b", "quadratic", PF["D13"]["check_src"])
for _key in ("D13", "D17", "D21", "D30"):
    _bare =[ln for ln in PF[_key]["check_src"].splitlines() if ln.strip() and "#" not in ln]
    assert not _bare, (_key, _bare)

# Starts that Part F gives mostly in words: the same statement as one display line (the words stay in "In words").
START_TEX = {
    "D01": r"p=p(\rho)\ \text{at fixed }S,\qquad\dfrac{\partial\rho}{\partial t}+\dfrac{\partial(\rho u_j)}{\partial x_j}=0,\qquad\dfrac{\partial u_i}{\partial t}+u_j\dfrac{\partial u_i}{\partial x_j}+\dfrac1\rho\dfrac{\partial p}{\partial x_i}=0",
    "D10": r"\rho uA=\rho^*c^*A^*\qquad\text{(mass, }\rho uA=\dot m=\text{const. (15.18), between a station and the sonic station)}",
    "D22": r"\rho u=\text{const},\qquad dp+d(\rho u^2)=-p_1\,df,\qquad h+\tfrac12u^2=\text{const}",
    "D23": r"\rho u=\text{const},\qquad p+\rho u^2=\text{const},\qquad dh+u\,du=h_1\,dq",
    "D25": r"\text{a front emitted at } t_e:\ \ \lvert\mathbf x-\mathbf x_s(t_e)\rvert=c\,(t-t_e),\qquad \mathbf x_s(t)=(Ut,0),\ U>c",
    "D26": r"\text{a plane shock at angle }\sigma\text{ to }V_1\text{, turning the flow through }\delta",
    "D29": r"v_1=v_2\qquad\text{(the velocity component along the wave is the same on both sides — D26 step 3)}",
}
START_WORDS = {
    "D01": "The pressure as a function of density at fixed entropy; continuity and Euler's equation without gravity (Ch. 4); and the split $" + EQ["15.4"] + "$ (15.4) — a steady uniform state plus a small disturbance.",
    "D10": "Mass conservation written between the station of Mach number M and the imagined sonic station, where the speed equals the sonic speed c*.",
    "D22": "The three laws of (15.45) with q = 0, applied between two stations a small distance apart: mass flux constant, momentum changed only by the friction of the slice, total enthalpy constant.",
    "D23": "The three laws of (15.45) with f = 0, between neighbouring stations: mass flux constant, pressure plus momentum flux constant, total enthalpy raised by the heat of the slice.",
    "D25": "A source moving at speed U through still gas of sound speed c, emitting all the time; a front emitted at one instant is a sphere growing at c from where the source was.",
    "D26": "The control-volume laws of C03 for a thin box enclosing a piece of the shock.",
    "D29": "Across a wave at angle σ the speed along the wave is unchanged (D26 step 3); step 1 writes this with the flow speeds and angles.",
}
for _k in START_TEX:
    PF[_k]["start"] = (START_TEX[_k], START_WORDS[_k])
# START_WORDS for D22/D23 name (15.45); keep the equation beside its number.
for _k in ("D22", "D23"):
    PF[_k]["start"] = (PF[_k]["start"][0], PF[_k]["start"][1].replace(
        "(15.45)", "$" + EQ["15.45"] + "$ (15.45)"))

# Verification report: at M1 = 1.1 the exact entropy jump is 4.29e-4 (the design prints 4.2e-4).
_old = "6.0 × 10⁻⁴ against the exact 4.2 × 10⁻⁴"
assert _old in PF["D17"]["check"], PF["D17"]["check"]
PF["D17"]["check"] = PF["D17"]["check"].replace(_old, "6.0 × 10⁻⁴ against the exact 4.3 × 10⁻⁴")

# One small move per step: D24 step 4 ("cross-multiply and take roots") made two moves; it becomes steps 4 and 5, and
# every later "step N" pointer of D24 moves up by one.
_STEP_REF = re.compile(r"\b(steps?) (\d+)(?:(\s?(?:[-–]|and)\s?)(\d+))?")


def _renumber(key: str, new_of) -> None:
    """Rewrite every "step N" / "steps N–M" in derivation ``key`` (all text fields and the check listing) with new_of(N)."""
    def fix(m):
        out = f"{m.group(1)} {new_of(int(m.group(2)))}"
        if m.group(4):
            out += f"{m.group(3)}{new_of(int(m.group(4)))}"
        return out

    d = PF[key]
    for f in ("goal", "assumptions", "check", "meaning", "traps"):
        d[f] = _STEP_REF.sub(fix, d[f])
    d["plan"] = [_STEP_REF.sub(fix, x) for x in d["plan"]]
    d["tools"] = [_STEP_REF.sub(fix, x) for x in d["tools"]]
    d["start"] = (d["start"][0], _STEP_REF.sub(fix, d["start"][1]))
    d["result"] = (d["result"][0], _STEP_REF.sub(fix, d["result"][1]))
    for st in d["steps"]:
        for f in ("did", "why", "plain"):
            st[f] = _STEP_REF.sub(fix, st[f])
    if d["check_src"]:
        d["check_src"] = _STEP_REF.sub(fix, d["check_src"])


assert PF["D24"]["steps"][3]["did"].startswith("Take the square root"), PF["D24"]["steps"][3]["did"]
_renumber("D24", lambda n: n + 1 if n >= 4 else n)
_old4 = PF["D24"]["steps"][3]
PF["D24"]["steps"][3:4] = [
    dict(did="Cross-multiply",
         tex=r"\left(1+\gamma M_1^2-f\right)^2=2(\gamma+1)\,M_1^2\,\Theta_1",
         why="Multiply both sides of step 3 by the two denominators, (1 + γM₁² − f)² and 2(γ + 1); both are non-zero, so "
             "nothing is lost. We do it to free the unknown f from the denominator.",
         plain="The square of the momentum bracket must equal a number fixed by the inlet Mach number."),
    dict(did="Take the square root",
         tex=_old4["tex"],
         why="A square equals a positive number exactly when its base is plus or minus that number's root; both signs "
             "are possible algebraically, and step 6 chooses between them.",
         plain=_old4["plain"]),
]
assert len(PF["D24"]["steps"]) == 8
N_STEPS = sum(len(d["steps"]) for d in PF.values())
assert N_STEPS == 243, N_STEPS                              # 242 of Part F + 1 (D24)

# Part F's *why* of these steps is a bare formula or a pointer; each gets the sentence that says why the move is allowed
# and what it is for (nbkit requires both).  The formula of Part F stays, word for word, in front.
WHY_ADD = {
    ("D05", 5): "Substituting the definition of the drifting coordinate back is always allowed; we do it to read the speeds in the fixed frame.",
    ("D06", 8): "Replacing a combination by its name changes nothing; we do it because h is the quantity tables and gas laws give.",
    ("D08", 1): "Dividing both sides of an equation by the same non-zero number keeps it true.",
    ("D08", 5): "We put the temperature ratio of step 3 into step 4, so that the pressure ratio depends on M alone.",
    ("D09", 2): "Subtracting two true equations side by side gives a true equation; it moves the problem to the stagnation states.",
    ("D14", 4): "Writing 1 as (γ + 1)/(γ + 1) puts everything on one denominator, ready to be inverted in step 5.",
    ("D14", 5): "Substituting one equation into another is allowed; it removes the pressures and leaves only Mach numbers.",
    ("D17", 11): "Substituting the exact relation of step 10 into the cubic changes nothing to the order kept; it writes the result with the measurable pressure jump.",
    ("D21", 2): "This trades the kinetic energy for the Mach number, as P359 does for the momentum flux, so that energy too is written in M.",
    ("D21", 6): "For a perfect gas h = C_pT, so the enthalpy ratio is the temperature ratio; substituting step 2 leaves only Mach numbers, f and q.",
    ("D21", 11): "Expanding the square and collecting like terms is ordinary algebra; the cancellation is what makes the root simple.",
    ("D22", 4): "Taking the logarithm of p = ρRT and differentiating is allowed for positive p, ρ, T; it gives relative changes, the language of the other steps.",
    ("D23", 1): "In a constant-area duct ρu is the same at every station, so its logarithmic differential vanishes.",
    ("D23", 5): "Dividing by the local enthalpy h = C_pT is allowed (h > 0) and turns the energy law into relative changes like steps 1–3.",
    ("D27", 7): "Multiplying top and bottom of a fraction by the same non-zero factor leaves it unchanged; this puts it in the book's shape.",
    ("D27", 8): "These are the identities cos² + sin² = 1 and cos 2σ = cos²σ − sin²σ; they shorten the denominator.",
    ("D29", 8): "Substituting steps 6 and 7 into step 5 and putting the bracket on one denominator is ordinary algebra.",
    ("D29", 9): "We put step 8 into step 4; it ties the turn to the Mach-number change alone, ready to integrate (D30).",
    ("D30", 10): "Undoing the substitution of step 1 is allowed because t and M are in one-to-one correspondence for M ≥ 1.",
}
for (_k, _i), _txt in WHY_ADD.items():
    _st = PF[_k]["steps"][_i - 1]
    assert len(_st["why"].split()) < 6, (_k, _i, _st["why"])        # only the steps that need it
    _st["why"] = _st["why"].rstrip(" .") + ". " + _txt
# Lesson review round 1: where "Part F formula + added sentence" said the same thing twice, the whole *why* is rewritten as
# one statement of the rule and the purpose (the maths of Part F kept).
WHY_SET = {
    ("D14", 3): r"Divide both sides of step 2 by $p_2$ (a pressure, never zero), which keeps the equation true; it leaves the "
                r"downstream bracket alone on the left and the pressure jump $p_2/p_1$ where step 4 can supply it.",
    ("D15", 1): r"Divide both sides of $" + EQ["15.35"] + r"$ (15.35) by $\rho_1u_2$, which is never zero; the mass law then "
                r"reads as a ratio, the form in which all the jumps are written.",
    ("D15", 3): r"Put step 2 into step 1 and take the reciprocal of both sides (allowed: neither side is zero); the density "
                r"jump is then a function of $M_1$ alone.",
    ("D16", 3): r"For a perfect gas $C_p-C_v=R$ and $C_p/C_v=\gamma$, so $C_v+R=C_p=\gamma C_v$; the two terms share the "
                r"same logarithm, so their coefficients simply add.",
    ("D16", 4): r"Two logarithm rules, $\gamma\ln x=\ln x^{\gamma}$ and $\ln a+\ln b=\ln(ab)$, merge the two terms; dividing by "
                r"$C_v>0$ makes the entropy jump dimensionless.",
    ("D14", 6): r"Subtract 1 from both sides of step 5 (the equation stays true) and write the right side over its own "
                r"denominator. The new numerator is $(\gamma+1)(1+\gamma M_1^2)-[2\gamma M_1^2-(\gamma-1)]"
                r"=(\gamma+1)+\gamma(\gamma+1)M_1^2-2\gamma M_1^2+(\gamma-1)=2\gamma+\gamma(\gamma-1)M_1^2$; this leaves "
                r"$\gamma M_2^2$ alone on the left.",
    ("D29", 5): r"Take the logarithm of $V=Mc$, $\ln V=\ln M+\ln c$ (allowed: all three are positive), and differentiate; "
                r"this splits the speed change into a Mach-number part and a sound-speed part, each of which steps 6–7 supply.",
    ("D21", 4): r"Divide both sides of step 3 by $p_2M_1/\sqrt{T_2}$ (none of these is zero), which solves the mass law for "
                r"the Mach-number ratio $M_2/M_1$; the pressure ratio and the temperature ratio left on the right are exactly "
                r"what steps 1 and 2 supply.",
    ("D21", 5): r"Divide step 1 by $p_2\left(1+\gamma M_1^2-f\right)$ (both factors positive); the result is the pressure ratio "
                r"that step 4 needs, supplied by the momentum law.",
    ("D27", 1): r"Multiply both sides by $\tan(\sigma-\delta)$ and divide by $\hat\rho$ (a density ratio, never zero); this "
                r"isolates the only factor that contains δ.",
    ("D28", 2): r"Replace the numerator $M_1^2\sin^2\sigma-1$ of $" + EQ["15.53"] + r"$ (15.53) by its equal from step 1; "
                r"swapping a quantity for an equal one keeps the line exact and brings the pressure jump into the turn formula.",
    ("D28", 3): r"The identity $\cos2\sigma=1-2\sin^2\sigma$ changes nothing (still exact); it writes the denominator in "
                r"$\sin\sigma$, the quantity whose limit step 5 supplies.",
    ("D28", 6): r"In the Mach triangle (primer P372) the side opposite μ is 1 and the hypotenuse M, so by Pythagoras the "
                r"adjacent side is $\sqrt{M^2-1}$ and $\cot\mu_1=\sqrt{M_1^2-1}$; σ → μ₁ carries it over.",
    ("D28", 7): r"By step 5, $M_1^2\sin^2\sigma\to1$, so $M_1^2\cdot2\sin^2\sigma\to2$, which cancels the $+2$; only "
                r"$(\gamma+1)M_1^2$ survives.",
    ("D30", 3): r"Put $M^2=1+t^2$ into the bracket and collect: $1+\tfrac{\gamma-1}{2}(1+t^2)=\tfrac{\gamma+1}{2}+\tfrac{\gamma-1}{2}t^2$; "
                r"factoring out $\tfrac{\gamma+1}{2}$ exposes the form $1+\lambda^2t^2$, whose integral is an arctangent.",
    ("D30", 6): r"$1-\lambda^2=\dfrac{(\gamma+1)-(\gamma-1)}{\gamma+1}=\dfrac{2}{\gamma+1}$, so the two constants are "
                r"reciprocals and their product is exactly 1; the integrand becomes a plain difference of two simple fractions.",
    ("D30", 8): r"The arctangent integral of primer P373, $\int_0^t dt/(1+a^2t^2)=\tan^{-1}(at)/a$, with $a=1$.",
    # D20 step 4: the step's line gives p_d > p_f; the reason p_d < p_b needs its own sentence.
    ("D20", 4): r"Just ahead of the exit the pressure is $p_f$, and the jump $" + EQ["15.39"] + r"$ (15.39) exceeds 1 for "
                r"$M_1>1$, so $p_d>p_f$. And $p_d<p_b$: write $p_d/p_0=(p_{02}/p_{01})\,(p/p_{02})$ at the exit. Behind the "
                r"shock $p_{02}<p_{01}$ (D18), and the sonic area grows to $A_2^*=A_t\,p_{01}/p_{02}$ (step 6), so the exit "
                r"area ratio $A_e/A_2^*$ is smaller than $A_e/A_t$ and the subsonic exit Mach number is larger than "
                r"$M_{\text{sub}}(A_e/A_t)$, which makes $p/p_{02}$ smaller. Both factors are below those of $p_b/p_0$.",
}
for (_k, _i), _txt in WHY_SET.items():
    PF[_k]["steps"][_i - 1]["why"] = _txt
# D31 quotes the book's "found from (15.47)": write the equation that number names (the friction differentials) beside it.
_q = '"found from (15.47)"'
_qq = (r'"found from (15.47)" — a slip: (15.47) is the friction differential '
       r'$-\dfrac{dp}{p_1}=\dfrac{1+(\gamma-1)M^2}{1-M^2}\,df$, nothing to do with a plate —')
_hits = 0
for _st in PF["D31"]["steps"]:
    _hits += _st["why"].count(_q)
    _st["why"] = _st["why"].replace(_q, _qq)
for _f in ("traps", "meaning", "check"):
    _hits += PF["D31"][_f].count(_q)
    PF["D31"][_f] = PF["D31"][_f].replace(_q, _qq)
assert _hits == 2, _hits
for _k, _d in PF.items():
    for _i, _st in enumerate(_d["steps"], 1):
        assert len(_st["why"].split()) >= 6, (_k, _i, _st["why"])

if "--dump" in sys.argv:
    for _k, _d in PF.items():
        print(_k, _d["stars"], len(_d["steps"]), _d["title"])
        print("   start:", _d["start"])
        print("   result:", _d["result"])
        for _i, _s in enumerate(_d["steps"], 1):
            print(f"   {_i}. {_s['did']} | {_s['tex'][:90]}")
    sys.exit(0)

STARS = {1: "★", 2: "★★", 3: "★★★"}


def _cap(s: str) -> str:
    """Capital first letter when the text starts with an ordinary lower-case word (never a symbol or a maths span)."""
    return s[:1].upper() + s[1:] if re.match(r"[a-z]{2,}(?=[\s,;:.—–-])", s) else s


def D(key: str, title: str, ref: str = "", notes: str = "", after: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then the notes it carries and its traps as a ⚠️ callout.
    ``title`` is our heading (Part F headings name bare equation numbers); ``ref`` = "15.39" shows the key equation."""
    d = PF[key]
    if ref in EQ:
        ref = f"{LABEL.get(ref, ref)}: ${EQ[ref]}$"         # the heading shows the key equation, not only its number
    goal = _cap(d["goal"]) + (f"\n\n**Assumptions.** {_cap(d['assumptions'])}" if d["assumptions"] else "")
    nb.derivation(key, f"{title} `{key}` {STARS[d['stars']]}", ref=ref, goal=goal,
                  start=(d["start"][0], _cap(d["start"][1])), plan=[_cap(x) for x in d["plan"]], uses=d["tools"],
                  steps=d["steps"], result=(d["result"][0], _cap(d["result"][1])), interpret=_cap(d["meaning"]),
                  check=_cap(d["check"]), check_src=d["check_src"])
    if notes:
        nb.md(X(f"📝 **Notes taught inside `{key}`:** {notes}"))
    if after:
        nb.md(X(after))
    if d["traps"]:
        nb.md(f"> ⚠️ **Common confusion (traps in `{key}`):** {_cap(d['traps'])}")


# ---------------------------------------------------------------------------------------------------------------------
# Small wrappers so that every block has the same shape (every text passes through X: «15.28» → the equation)
# ---------------------------------------------------------------------------------------------------------------------
def core(cid: str, title: str, question: str) -> None:
    """A CORE block heading with its id (the title carries the block's key equation in LaTeX)."""
    nb.core(cid, f"{X(title)} `{cid}`", question=X(question))


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` is exactly the Part E concept text (the ledger check matches it)."""
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
    """A B/C note whose text opens with its curation id(s) in bold ("**N09 [B]**")."""
    nb.note(f"**{nid}** {X(text)}", equation=equation, ref=ref)


def md(text: str) -> None:
    nb.md(X(text))


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + X(text))


def idea(sketch: str = "", words: str = "") -> None:
    """The idea block: an ASCII sketch (or a table) plus one or two sentences."""
    sk = textwrap.dedent(sketch).strip("\n")
    body = ("```\n" + sk + "\n```") if sk else ""
    nb.md("#### The idea\n\n" + body + (("\n\n" if body else "") + X(words) if words else ""))


def whatif(text: str) -> None:
    nb.md(f"**What would change if…** {X(text)}")


def confusion(text: str) -> None:
    nb.md("> ⚠️ **Common confusion:** " + X(text))


def slip(k: int, printed: str, correct: str, why: str = "") -> None:
    """A printed-vs-correct box (the book's printed slip #k, taught in corrected form)."""
    nb.md(f"> ⚠️ **slip #{k} — the book prints** {X(printed)}; **the correct form is** {X(correct)}"
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
    """A from-scratch cell (hand-written version next to the tested function, then an assert) with its explanation."""
    nb.code(src, explain=X(explain), tags=["from-scratch"])


def figure(src: str, see: str, read: str, change: str, explain: str | None = None) -> None:
    nb.figure(src, see=X(see), read=X(read), change=X(change), explain=X(explain) if explain else None)


def notes3(see: str, read: str, change: str) -> None:
    nb.figure_notes(see=X(see), read=X(read), change=X(change))


def worked(title: str, steps_md: str) -> None:
    nb.worked_example(title, X(steps_md))


EXPLAINER_NOTE = {   # one conventions sentence per explainer, about what that explainer actually shows
    "sound_pulse_dalembert": "The pressure is the disturbance p′ about the undisturbed air, its size set as the dimensionless amplitude p′/ρ₀c²; the probe shows p′ and the air speed u′ at one point.",
    "area_mach_nozzle": "Areas are ratios to the sonic area A*; station 1 is the reservoir side, and every ratio is to the stagnation state (subscript 0).",
    "normal_shock_jump": "Station 1 is upstream (ahead of the shock), station 2 downstream; every jump is written downstream over upstream.",
    "nozzle_back_pressure": "Pressures are fractions of the supply (stagnation) pressure p₀; the back pressure p_B is the pressure in the tank outside the exit.",
    "duct_friction_heating": "f and q are the book's lumped parameters F/(p₁A) and Q/h₁; station 1 is the inlet.",
    "mach_cone": "The Mach angle μ is computed in radians and shown in degrees; the source moves to the right.",
    "oblique_shock_wedge": "Station 1 is upstream; σ and δ (β and θ in most other books) are computed in radians and shown in degrees.",
    "prandtl_meyer_fan": "ν, μ and the turn are computed in radians and shown in degrees; the code counts an expansion turn positive (the book's δ is its negative).",
    "supersonic_flat_plate": "The incidence α is computed in radians and shown in degrees; b is the chord, and pressures are over the free-stream p∞.",
}


def explainer(slug: str, heading: str, why_static: str, tries: list[str]) -> None:
    """An embedded explainer, opening with the "why interactive rather than a static figure" sentence."""
    why = (f"**Why interactive rather than a static figure:** {X(why_static)}\n\n"
           "Start with the **Walkthrough**; the **Explain** tab works every number out with your settings, and the "
           "**Derivation** tab (where there is one) steps through the same derivation as above. " + EXPLAINER_NOTE[slug])
    nb.explainer(slug, heading=heading, why=why, tries=[X(t) for t in tries])


# ---------------------------------------------------------------------------------------------------------------------
# One-line reminders of tools primed in earlier chapters (Part E rows that point to knowledge/primers.md)
# ---------------------------------------------------------------------------------------------------------------------
def _ledger_rows() -> list[list[str]]:
    part = DESIGN.split("## Part E", 1)[1].split("## Part F", 1)[0]
    rows = [[x.strip() for x in ln.strip().strip("|").split("|")] for ln in part.splitlines() if ln.strip().startswith("|")]
    return [r for r in rows if len(r) >= 3 and not set("".join(r)) <= set("-: ") and r[0].lower() != "concept"]


LEDGER = _ledger_rows()
assert len(LEDGER) == 191, len(LEDGER)

REMIND: dict[str, str] = {        # our own sentence about how the earlier tool is used HERE
    "Taylor expansion to first and second order": "close to a point, f(ρ₀ + ρ′) ≈ f(ρ₀) + f′(ρ₀)ρ′ + ½f″(ρ₀)ρ′², with each further term smaller by another factor ρ′; here it expands the pressure as a function of the density about the undisturbed state.",
    "orders of smallness; dropping products of primes": "when every primed quantity is small, a product of two primes is smaller still and is dropped beside a term with one prime; this one rule turns the nonlinear laws of mass and momentum into linear ones.",
    "partial derivative; index notation and the summation convention; the Laplacian $\\partial^2/\\partial x_j\\partial x_j$": "∂/∂t holds x fixed and ∂/∂x holds t fixed; a repeated index is summed over its three values (Ch. 2), so U_j ∂/∂x_j is U·∇ and ∂²/∂x_j∂x_j is the Laplacian ∇², the curvature of a field in all three directions.",
    "product rule": "d(ab) = a db + b da; here it expands the mass flux (ρ₀ + ρ′)(U + u′) and the product pA along a duct.",
    "operator elimination for linear equations": "apply a derivative to a whole linear equation and take the divergence of another so that one unknown appears in both with the same operator, then subtract; here it removes the disturbance velocity and leaves one equation for the pressure.",
    "`sympy` symbols, `diff`, `expand`, `simplify`, `series`, `limit`, `solve`": "sympy does algebra with symbols: `sp.symbols` names them, `sp.diff` differentiates, `sp.expand` multiplies out, `sp.simplify(a - b) == 0` proves two expressions equal, `sp.series` expands in a small quantity, `sp.limit` takes a limit and `sp.solve` solves; every ★★★ derivation below is re-run with it.",
    "`assert np.allclose`": "`assert np.allclose(a, b)` stops with an error if two results differ by more than a tiny tolerance; here it is the proof that each from-scratch cell matches the tested `fluidpy` function.",
    "`slider_figure`, `animate` / `show_animation`, `show_viz`, `live`": "`slider_figure(fn, name, values)` draws a plotly figure with every slider position computed in advance (it keeps working on the web page); `animate` + `show_animation` play a matplotlib animation inline; `show_viz` embeds an explainer; `live` makes sliders that need a running kernel.",
    "dataclasses": "a dataclass is a small named record; `GD.PerfectGas(1.4, 287.05)` is one, holding γ and R and computing $C_p$, $C_v$ and the sound speed from them.",
    "first-order wave equation and characteristics": "∂q/∂t + c ∂q/∂x = 0 keeps q constant for an observer moving at c, so its solution is a shape sliding at c along the lines x − ct = const; the sound equation below has two such families.",
    "chain rule (two new variables)": "if p′ depends on ξ and ζ, and these depend on x and t, then ∂p′/∂x = (∂p′/∂ξ)(∂ξ/∂x) + (∂p′/∂ζ)(∂ζ/∂x), and the same for t; it translates the wave equation into the coordinates that ride with the waves.",
    "frames of reference and the Galilean shift": "an observer moving at constant speed U measures positions x′ = x − Ut and the same time; velocities change by U, but pressure, density and temperature do not — used for sound in a stream, for moving shocks and for oblique shocks.",
    "domain of dependence": "the value at a point (x, t) of a hyperbolic problem depends only on the data inside the wedge cut out by the characteristics through it; here the wedge of a sound wave is bounded by the lines of slope U + c and U − c.",
    "definite integral": "∫ f dt between two limits adds up f dt over that range; here the time integral of the pressure gradient gives the air's velocity under a sound wave.",
    "power of a force; heat flux through a surface": "a force does work at the rate force × velocity, and heat q_j crossing a surface delivers −∫q_j n_j dA; here the pressure on the inlet and outlet faces of a duct does work, and the heat added per kilogram is Q.",
    "momentum flux through a surface": "fluid crossing a patch dA carries its momentum with it at the rate ρu(u·n) dA; here the inlet and outlet faces of a duct carry −ṁu₁ in and +ṁu₂ out.",
    "logarithmic differentiation": "d(ln f) = df/f turns a product of powers into a sum of relative changes; here ρuA = const becomes dρ/ρ + du/u + dA/A = 0, and p = ρRT becomes dp/p = dρ/ρ + dT/T.",
    "binomial approximation": "(1 + x)^n ≈ 1 + nx + ½n(n − 1)x² for small x; here it expands the stagnation pressure ratio for a slow flow and recovers the incompressible pitot formula.",
    "`scipy.optimize.brentq`": "`brentq(f, a, b)` finds the root of f inside a bracket [a, b] where f changes sign; here the bracket also chooses which of two Mach numbers you get.",
    "l'Hôpital's rule (the 0/0 at the throat)": "when the top and the bottom of a fraction both vanish, its limit is the limit of the ratio of their derivatives; at a sonic throat dA = 0 and 1 − M² = 0 together, and the area rule has a finite slope there.",
    "Vieta's formulas": "for a quadratic ax² + bx + c the two roots multiply to c/a and add to −b/a; here knowing one root of the shock quadratic (\"nothing happens\") hands you the other.",
    "caching a slow cell (`functools.lru_cache`)": "a function wrapped in `functools.lru_cache` remembers its answers, so a slow symbolic calculation runs once per session; used here for the two slowest sympy helpers when FAST is on.",
    "`scipy.integrate.solve_ivp`": "`solve_ivp(rhs, (t0, t1), y0)` integrates an ordinary differential equation numerically; here it integrates the shock-structure equation as an independent check of the implicit profile.",
    "partial fractions; inverting an implicit solution": "a ratio of polynomials splits into simple fractions that integrate one by one, and a solution given as η(U) is turned round into U(η) by a root-finder; both are used for the inside of a shock and for the Prandtl–Meyer integral.",
    "`np.expm1`, `np.log1p` and cancellation (inside the library)": "`np.log1p(x)` computes ln(1 + x) accurately when x is tiny, where ln(1 + x) would lose digits; the library uses it for the entropy jump of very weak shocks.",
    "piecewise functions with `np.where` / `if`": "a quantity given by different formulas in different ranges is computed with an `if` (one value) or `np.where` (arrays); here the exit of a convergent nozzle is subsonic above the critical back pressure and sonic below it.",
    "the quadratic formula": "aX² + bX + c = 0 has the roots X = [−b ± √(b² − 4ac)]/(2a); here X is the square of the outlet Mach number.",
    "parabolic and elliptic problems": "an elliptic problem (subsonic flow, Laplace's equation) needs conditions all round because every point influences every other; a parabolic one can be marched in one direction — steady supersonic flow is a third kind, hyperbolic (primer below).",
    "small-angle approximation": "for an angle ε in radians, sin ε ≈ tan ε ≈ ε and cos ε ≈ 1, with errors of order ε³ and ε²; here it turns the oblique-shock formula into the weak-wave rule.",
    "`scipy.optimize.minimize_scalar` (cross-check of the maximum)": "`minimize_scalar(f, bounds=(a, b), method=\"bounded\")` finds the minimum of a function of one variable inside a bracket; here (on −δ) it cross-checks the largest deflection a shock can give.",
    "substitution in an integral": "replace the integration variable by a new one and convert the differential and the limits as well; here t = √(M² − 1) removes the square root from the Prandtl–Meyer integral.",
    "`scipy.integrate.quad`": "`quad(f, a, b)` integrates a function numerically and returns the value and an error estimate; here it evaluates the Prandtl–Meyer integral directly, as a check of its closed form.",
}


def _norm(c: str) -> str:
    return re.sub(r"[`*$\\\"“”]", "", c).lower().strip()


REMIND_N = {_norm(k): v for k, v in REMIND.items()}
_NO_ROW = sorted(k for k in REMIND_N if k not in {_norm(r[0]) for r in LEDGER})
assert not _NO_ROW, f"REMIND keys that match no Part E row: {_NO_ROW}"


def remind(cid: str, lead: str = "") -> None:
    """One 🔁 cell reminding the tools primed earlier that Part E lists as first used in ``cid`` — one sentence each, the
    concept text verbatim so that the ledger check matches."""
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
Until now the fluid could not be squeezed. Let it be squeezed and one new speed enters — the speed at which the fluid
passes on the news that something has changed, the **speed of sound** $c$. Everything in this chapter is three
conservation laws (mass, momentum, energy) and the gas law, written between two stations in terms of one number, the
Mach number «4.111». In every answer the same factor $1-M^2$ appears, and **it changes sign at $M=1$**. That one sign
change is why a widening duct speeds up a supersonic stream, why friction speeds up a subsonic duct flow, why adding heat
can cool a gas, why flows jump through shocks, and why a supersonic aircraft is heard only after it has passed."""),
    roadmap=[
        "**Sound as the linear limit** — weak squeezes obey one wave equation with speeds U ± c (`C01`).",
        "**Two shapes sliding apart** — every plane sound wave is one shape going right plus one going left (`C02`).",
        "**The three duct laws** — mass, energy and momentum between two stations of a duct (`C03`).",
        "**The stagnation ruler** — the state a gas would reach if stopped gently; $T_0$ kept, $p_0$ lost (`C04`).",
        "**One area, two Mach numbers** — the area a flow needs is smallest at Mach 1 (`C05`).",
        "**The area–velocity surprise** — above Mach 1 a nozzle must widen (`C06`).",
        "**The normal-shock jump** — three laws, one extra solution, fixed by $M_1$ (`C07`).",
        "**The second law picks the direction** — shocks only compress; they cost stagnation pressure (`C08`).",
        "**Choking** — once the exit is sonic, the nozzle stops listening to the back pressure (`C09`).",
        "**The Laval nozzle** — the back pressure only chooses where a shock stands (`C10`).",
        "**Friction and heat in one formula** — one number A and a quadratic for the outlet Mach number (`C11`).",
        "**Everything drifts to Mach 1** — friction speeds up a subsonic flow; heat can cool a gas (`C12`).",
        "**The Mach cone** — the edge of all the circles of old news (`C13`).",
        "**The oblique shock** — a normal shock on the normal component; two shocks per turn, a largest turn (`C14`).",
        "**The Prandtl–Meyer fan** — turning away one Mach wave at a time (`C15`).",
        "**Lift and wave drag of a flat plate** — the push is normal to the plate, so lift costs drag (`C16`).",
    ],
    prerequisites=[
        X(r"the perfect gas, entropy and the speed of sound «1.19» (Ch. 1)"),
        X(r"the control-volume laws and the Mach number «4.111» (Ch. 4)"),
        X(r"the hydraulic jump «7.81» (Ch. 7)"),
        X(r"the lift slope of «14.12» (Ch. 14)"),
    ],
)
nb.explainer_index([
    ("sound_pulse_dalembert", "What is a sound pulse made of?", "two shapes sliding apart at ±c — and both swept downstream when the stream beats c"),
    ("area_mach_nozzle", "Why does a supersonic nozzle widen?", "above Mach 1 the gas thins faster than it speeds up"),
    ("normal_shock_jump", "Why can a shock only compress?", "three laws balance both ways; entropy allows one"),
    ("nozzle_back_pressure", "Why does a nozzle stop listening?", "once the throat is sonic the back pressure only moves a shock"),
    ("duct_friction_heating", "Can friction speed a flow up?", "friction and heat push every duct flow toward Mach 1"),
    ("mach_cone", "Why do you hear it only after it has passed?", "the cone is the edge of all the circles of old news"),
    ("oblique_shock_wedge", "One wedge, how many shocks?", "two below the largest turn, none above it"),
    ("prandtl_meyer_fan", "How does a supersonic stream turn a corner?", "one Mach wave at a time; ν(M) keeps the count"),
    ("supersonic_flat_plate", "How can there be drag without friction?", "the push is square to the plate, so lift costs α times itself"),
])
nb.setup()

md(r"""
## 📦 Chapter imports and the numbers we will use

Every worked number in this notebook uses **our own illustrative inputs**, written out in the cell below as six small
dictionaries: the gas (`air`), a nozzle (`noz`), a rocket (`rocket`), a duct (`duct`), a wedge (`wedge`) and a flat plate
(`plate`). None of them is a number from the book. The module carries its own illustrative set too; it is shown once at the
end of the cell — both sets are ours; neither is the book's.""")
code(r"""
import numpy as np                                          # arrays and maths
import sympy as sp                                          # symbolic algebra (derivation checks)
import matplotlib.pyplot as plt                             # static figures and animations
import pandas as pd                                         # small tables
import functools                                            # lru_cache: remember the result of a slow call
from scipy.optimize import brentq                           # root of f(x) = 0 inside a bracket
from scipy.integrate import quad, solve_ivp                 # numerical integral; ODE integrator
from fluidpy import ch15_compressible_flow as ch15          # this chapter's tested functions (re-exports GD and AC)
from fluidpy.core import gasdynamics as GD, acoustics as AC # the two new core modules: gas dynamics, acoustics
from fluidpy.core import thermo as TH, bernoulli as BER, similarity as SIM, statics as STAT, airfoils as AF   # earlier chapters' tools
from fluidpy import ch01_introduction as ch01, ch07_gravity_waves as ch07, ch14_aerodynamics as ch14   # earlier chapter modules
from fluidpy.core.interact import slider_figure, animate_figure, live   # plotly slider figures and kernel-only widgets
from fluidpy.core.anim import animate                       # builds a matplotlib animation from an update function
from fluidpy.core.style import COLORS, savefig              # house colours; save a figure to outputs/ch15/
sys.path.insert(0, "scripts")                               # the folder with this chapter's drawing helpers
from ch15_drawings import *                                 # sketches (draw_…) that only draw what the tested functions return
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # silence harmless font messages

air = GD.PerfectGas(1.4, 287.05)                            # the gas: gamma = 1.4, R = 287.05 J/(kg K) (project air)
noz = dict(p0=5.0e5, T0=400.0, A_t=1.0e-3, Ae_over_At=2.0)  # nozzle: supply pressure [Pa], supply temperature [K], throat area [m²], exit/throat area
rocket = dict(thrust=5.0e4, p0=3.0e6, T0=3000.0, p_exit=101325.0)   # rocket: thrust [N], chamber pressure [Pa], chamber temperature [K], exit pressure [Pa]
duct = dict(M1=0.5, f=0.1, q=0.2)                           # duct: inlet Mach number, lumped friction f = F/(p1 A), lumped heat q = Q/h1
wedge = dict(M1=2.0, delta=np.deg2rad(10.0))                # wedge: upstream Mach number, half-angle 10° (stored in rad)
plate = dict(alpha=np.deg2rad(2.0), M_inf=2.0, b=1.0, p_inf=101325.0)   # plate: incidence 2° (rad), Mach number, chord [m], pressure [Pa]
c0 = air.c(288.15)                                          # sound speed of sea-level air [m/s], c = sqrt(gamma R T)
rho0 = 101325.0 / (air.R * 288.15)                          # sea-level density [kg/m³] from p = rho R T
print(f"c0 = {c0:.1f} m/s, rho0 = {rho0:.4f} kg/m³, cp = {air.cp:.1f} J/(kg K), rho0*c0 = {rho0*c0:.1f} kg/(m² s)")
pd.Series(ch15.illustrative_inputs()).to_frame("the module's own illustrative set (its scripts and tests use it; also ours)").T   # shown once, beside ours
""", explain=r"""
1. The first block imports numpy, sympy, matplotlib, pandas and two scipy tools, then the chapter module `ch15` and the two new
   core modules: `GD` (one-dimensional gas dynamics, shocks, nozzles, ducts, oblique waves) and `AC` (acoustics). `ch15`
   re-exports every public name of both.
2. The second block brings back earlier chapters' tools: `TH` (thermodynamics, Ch. 1), `BER` (Bernoulli, Ch. 4), `SIM` (the
   Mach number, Ch. 4), `STAT` (the standard atmosphere, Ch. 1), `AF` (airfoils, Ch. 14) and three chapter modules.
3. `from ch15_drawings import *` loads the sketch helpers (`draw_…`): pure matplotlib on top of the tested functions.
4. The six dictionaries are the inputs of every worked number below. Angles are **stored in radians** (`np.deg2rad`) and
   **quoted in degrees** in the text.
5. The printed line gives sea-level air: $c_0$ = 340.3 m/s, $\rho_0$ = 1.225 kg/m³, $C_p$ = 1004.7 J/(kg K).
6. The last line shows the module's own illustrative set once (for example a 12° wedge at Mach 2.5). Both sets are ours.""")

md(r"""
## ⚠️ Conventions in this chapter — read once, return when a symbol surprises you

**(a) Here against most other sources.** The table printed by the next cell lists the nine places where this chapter's
symbols differ from the usual ones: the subscript 0 (three meanings — trap T1; most tables write $p_t$, $T_t$ for the
stagnation state), the star (three meanings, and an $A$ that is not an area — trap T2), the stations (1 upstream, 2
downstream — trap T6), the oblique-shock angles (**σ = β elsewhere**, **δ = θ elsewhere** — trap T13), the sign of a turn
(the book counts δ > 0 when the wall turns *into* the flow; our `turn=` counts an expansion positive — trap T14), friction and
heat (the book's lumped $f=F/(p_1A)$ — total wall force F over inlet pressure times area — and $q=Q/h_1$; the usual
"Fanno" friction length $4f_FL/D$, with $f_F=\tau_w/(\tfrac12\rho u^2)$ the Fanning friction factor (wall stress over dynamic
pressure), L the duct length and D its hydraulic diameter, and the usual "Rayleigh" ratio $T_0/T_0^*$ are *not in the book —
ours*; both are set side by side in N86, C12), the wall force (trap T11), entropy (S here, s in Ch. 1).

**(b) Overloaded letters (trap T5).**

| Letter | Meanings in this chapter | Code name |
|---|---|---|
| f, g | the two wave shapes of «15.15» (§15.2 only); **f** also the friction parameter of «15.45» (§15.8) | `f`, `g` (callables) in `AC`; `f` (number) in `GD.duct_*` |
| q, Q | Q heat added per unit mass in «15.22»; $q=Q/h_1$ dimensionless; $q_j$ the heat-flux vector of «15.20» | `Q`, `q` |
| A | an area; **a pure number** in «u15.8» | `A`, `A_over_Astar`, `duct_parameter_A` |
| δ | deflection angle (§15.9–§15.10); shock thickness $\delta_{\text{sh}}$ (§15.6) | `delta`, `thickness` |
| μ | Mach angle (§15.9); viscosity $\mu_{\text{visc}}$ (§15.6, with $\mu''=\tfrac43\mu+\mu_v$) | `mu` (angle), `mu_visc` |
| ν | Prandtl–Meyer function (§15.9); kinematic viscosity $\nu_{\text{visc}}$ (§15.6) | `nu`, `nu_visc` |
| σ | shock angle (Ch. 7: surface tension; Ch. 11: growth rate) | `sigma` |
| U, V, v | U stream speed (§15.2); $U=u/V$ inside a shock (§15.6); V flow speed, v its component along a shock (§15.9) | `U`, `U_scaled`, `V`, `v_t` |
| b | the **chord** in «15.58» (Ch. 14: the map constant; span s) | `b` |
| c, m, e, h, η | c sound speed (Ch. 14: chord); m mass flux inside a shock; e internal energy; h enthalpy (and a height in «o_delay»); η shock coordinate | `c`, `m`, `e`, `h`, `height`, `eta` |

**(c) Slips and traps.** *Slips* are false as printed; *traps* are true but easy to misread. **Three slips change the
physics and are taught corrected**: «15.32», «15.47» and «15.48». The other fourteen are named where they occur.

**(d) Which inputs?** All worked numbers use our own illustrative inputs, never the book's.

**(e) Which names?** The book never says "Fanno" or "Rayleigh"; where we add those names and their usual parameters they
are labelled *not in the book — ours*. The same label marks every formula of §15.7, which has no numbered equation.

**(f) Lapse rate.** No lapse rate enters this chapter's equations. Where the standard atmosphere is used (the rocket of
C06) and its lapse rate is quoted, both conventions are given: Kundu $\Gamma\equiv dT/dz=-6.5$ K/km, meteorology
$-dT/dz=+6.5$ K/km.""")
code(r"""
pd.set_option("display.max_colwidth", 120)                  # show long table cells without cutting them
display(ch15.conventions_table())                           # nine rows: item, here, elsewhere (most other sources), code
slips = ch15.book_slips()                                   # 17 rows: n, kind, where, printed, correct, code, changes_physics
display(slips[["n", "kind", "where", "printed", "correct", "changes_physics"]])   # the three physics-changing slips are marked True
display(ch15.traps())                                       # 22 rows: id, trap, guard (what the code does about it)
""", explain=r"""
1. `ch15.conventions_table()` is a pandas table (`pd.DataFrame`: rows and named columns) of the nine conventions of (a).
2. `ch15.book_slips()` lists the seventeen slips; the column `changes_physics` is True for the three whose corrected form
   changes a result (#4, #11, #12) — the code carries the printed form only as an option a test must reject.
3. `ch15.traps()` lists the twenty-two traps with what the code does about each; each one reappears as a "Common confusion"
   callout where it first bites.""")

# =====================================================================================================================
# A.1 §15.1 Introduction
# =====================================================================================================================
nb.section("15.1", "Introduction", intro="**What is this section about?** What 'compressible' adds, the gas we assume all "
           "chapter long, the one number everything is written in, and three surprises the chapter promises to explain.")
note("N01 [B]", r"""
**Gas dynamics** is the study of flows fast enough for the density to change. Standing assumptions for the whole chapter:
a perfect gas with constant $\gamma$; no gravity; no friction outside thin boundary layers (except in §15.8, where friction is
the subject). **Three promises** the chapter keeps:

| Promise | Kept in | The equation that keeps it |
|---|---|---|
| a widening duct speeds up a supersonic stream | C06 | «15.34» |
| friction speeds up a subsonic duct flow | C12 | «15.47» |
| adding heat can lower the temperature | C12 | «15.48» |""")
recap("R01", "The Mach number", r"""«4.111»: the flow speed over the local speed of sound. **New here:** M is a *local*
number that changes along a duct, and every result of this chapter is written in it.""", where="Ch. 4 §4.11 (C15)")
code(r"""
print(SIM.mach_number(250.0, T=220.0))                      # (4.111): M = U/c with c = sqrt(gamma R T) at 220 K — an airliner at cruise
print(GD.mach_from_speed(250.0, 220.0, air), air.c(220.0))  # the same with this chapter's gas object; c = 297.3 m/s
# the two Mach numbers differ in the 5th digit (0.840775 vs 0.840787) because the Ch. 4 module uses R = 287.058 and this
# chapter's gas uses R = 287.05 J/(kg K) — trap T4, explained with the properties of air below
""")
recap("R02", "The speed of sound", r"""«1.19»: how much pressure it takes to squeeze the gas at constant entropy. **New
here:** it is the coefficient of the first term of the expansion of p(ρ) in C01 — "sound speed" and "stiffness at constant
entropy" are one idea.""", where="Ch. 1 §1.8 (C36)")
code(r"""
# `lambda` (Ch. 1 primer P29): a one-line function without a name — `lambda rho, s: ...` takes rho and s and returns the expression
p_of_rho = lambda rho, s:101325.0 * (rho / rho0) ** 1.4    # an isentrope of sea-level air, p ∝ rho^gamma (s is held fixed)
print(TH.sound_speed_from_eos(p_of_rho, rho0))              # (1.19): sqrt of dp/drho at constant entropy -> 340.3 m/s
""")
recap("R03", "Compressibility scales with M²", r"""«4.110»: the divergence of the velocity — the rate of volume change — is
of order $M^2$, and density changes are about $\tfrac12M^2$. **New here:** this chapter is everything above that line; the
exact ratio is «15.30» in C04.""", where="Ch. 4 §4.11; Ch. 14 (N01)")
code(r"""
print(ch14.compressibility_check(0.3 * c0))                 # Ch. 14: density change ≈ ½ M² at M = 0.3 -> 0.045
print(GD.stagnation_density_ratio(0.3) - 1)                 # (15.30): the exact relative density change -> 0.0456
""")
note("N02 [B]", r"""
Five regimes, named by the Mach number:

| Regime | Range | Everyday example |
|---|---|---|
| incompressible | $M\approx0$ | a car, a breeze |
| subsonic | $0<M<1$ | an airliner at 0.84 |
| transonic | $M\approx1$ | the same airliner's wing upper surface |
| supersonic | $M>1$ | a rifle bullet at 2.5 |
| hypersonic | $M\gg1$ | re-entry at 25 |""")
code(r"""
print([ch15.mach_regime(M, transonic=(0.8, 1.2), hypersonic=5.0) for M in (0.1, 0.5, 1.0, 2.5, 25.0)])   # our edges, passed by keyword
print(ch15.mach_regime(0.0))                                # "incompressible" only at M = 0: a limit, not a band
""")
confusion(r"""(trap T3) The edges of the regimes are conventions, not laws of nature. We pass ours by keyword — 0.8 and 1.2,
and the customary 5; the book places the hypersonic edge lower.""")
note("N03 [C]", r"""The name "Mach number" is Ackeret's, after Ernst Mach, who photographed the waves of supersonic
bullets; "transonic" is von Kármán's. Mach's wave-front picture returns in C13 (§15.9), Ackeret's linear theory in C16
(§15.10).""")
recap("R04", "The perfect gas with constant specific heats", r"""One screen: $e=C_vT$, $h=C_pT$, $p=\rho RT$; «15.1h»;
«15.1c»; «15.1s». **New here:** entropy is written S (s in Ch. 1); its two forms are the tool of C04 and C08.

| Relation | `TH` / `GD` function |
|---|---|
| $p=\rho RT$, $C_p$, $C_v$ | `GD.PerfectGas(gamma, R)` with `.cp`, `.cv` |
| $c=\sqrt{\gamma RT}$ | `.c(T)`, `TH.perfect_gas_sound_speed(T)` |
| entropy change | `TH.perfect_gas_entropy_change_p(T1, p1, T2, p2, cp, R)` |""", where="Ch. 1 §1.8–§1.9 (C28–C47)")
code(r"""
print(air.cp, air.cv, air.c(288.15))                        # (15.1): C_p = gamma R/(gamma-1) = 1004.7, C_v = 717.6 J/(kg K); c = 340.3 m/s
""")
recap("R05", "Isentropic relations", r"""«15.2» — between any two states with the same entropy.""", where="Ch. 1 (C45, C46)")
P("P355", "the gamma exponents", r"""Three exponents appear all chapter. For γ = 1.4: $\gamma/(\gamma-1)=3.5$,
$1/(\gamma-1)=2.5$, $(\gamma+1)/(2(\gamma-1))=3$. For a monatomic gas (γ = 5/3): 2.5, 1.5, 2. The pressure exponent is the big
one; the classic slip is to invert them. In numpy a fractional power is just `x**3.5`.""", code=r"""
g = 1.4                                       # ratio of specific heats of air
print(g/(g-1), 1/(g-1), (g+1)/(2*(g-1)))      # 3.5 2.5 3.0
print(0.5**3.5, 0.5**2.5)                     # halve T at constant entropy: p falls to 8.8 %, rho to 17.7 %
""")
code(r"""
print(TH.isentropic_ratios(0.5))                            # (15.2) at p2/p1 = 0.5: T2/T1 = 0.5**(0.4/1.4) = 0.820, rho2/rho1 = 0.5**(1/1.4) = 0.610
""")
recap("R06", "Properties of air", r"""«15.3» are the chapter's rounded values. > ⚠️ Common confusion (trap T4): 1004/717 =
1.4003, not 1.40; the project's constants (R = 287.05, $C_p$ = 1004.7) differ in the fourth digit. `GD.AIR_BOOK` is passed
explicitly whenever a book number is reproduced (only in private tests).""", where="Ch. 1 constants")
code(r"""
print(1004 / 717, GD.AIR_BOOK.cp, air.cp)                   # 1.4003; 1004.5 (R = 287.0); 1004.7 (R = 287.05): rounding shows in the 4th digit
""")

# =====================================================================================================================
# A.2 §15.2 Acoustics — C01
# =====================================================================================================================
nb.section("15.2", "Acoustics", intro="**What is this section about?** A sound is a disturbance so weak that its products "
           "with itself can be dropped. That turns mass and momentum into one linear equation whose only speeds are U ± c — "
           "and its solution is two shapes sliding apart.")
core("C01", r"Sound is the linear limit: $\dfrac1{c^2}\left(\dfrac{\partial}{\partial t}+U_j\dfrac{\partial}{\partial x_j}\right)^2p'-\dfrac{\partial^2p'}{\partial x_j\partial x_j}=0$ (15.12)",
     "Why does a faint squeeze of the air travel at one fixed speed, whatever its shape?")
problem(r"""
Clap your hands. A thin layer of air is squeezed, pushes on the next layer, which pushes on the next. Nobody across the room
feels a wind; they feel the *news* of the squeeze arriving. How fast does news travel in a gas, and what decides it? The
answer is the stiffness of the gas against a quick squeeze — and the same number $c$ will reappear in every later section,
as the speed the flow is compared with.""")
idea(r"""
full laws (nonlinear) ──► split: state = steady + small ′ ──► drop ′×′ ──► two linear laws ──► eliminate u′ ──► one wave equation
""", r"""**Small** means «15.7»: the pressure disturbance is tiny compared with the stiffness $\rho_0c^2=\gamma p_0$.""")
note("N04 [B]", r"""**Acoustics** is the linear theory of small, isentropic disturbances about a steady reference state.
"Small" means that products of two disturbances are dropped; «15.7» makes it precise.""")
note("N05 [B]", r"""The split «15.4»: every quantity = its steady, uniform reference value + a small disturbance (primed).""")
recap_md("Basic state plus small disturbance", "Ch. 11 C01", r"""Write each field as a known steady state plus a small
disturbance, put it into the full equations, and keep only the terms with one disturbance factor — the recipe of the
stability chapter, used here for sound.""")
confusion(r"""(trap T1) The subscript 0 means three things in this chapter:

| Where | $p_0$, $T_0$ mean | Other sources |
|---|---|---|
| §15.1 | a reference value (a scale) | — |
| §15.2 (here) | the undisturbed state of the sound wave | — |
| from §15.4 on | the **stagnation** state | $p_t$, $T_t$ ("total") |""")
recap("R07", "Sound speed at the reference state", r"""«15.6». **New here:** $T_0$ is the *undisturbed* temperature (trap
T1).

> ⚠️ **slip #1 — the book prints** "monotonic" gases; **the correct form is** monatomic (the argument is about γ, which is
5/3 for a monatomic gas).""", where="Ch. 1 (C47)", core_id="C01")
code(r"""
print(TH.perfect_gas_sound_speed(288.15))                   # (15.6): c = sqrt(gamma R T0) for air at 15 °C -> 340.3 m/s
""")
recap("R08", "A number to remember", r"""Air at 15 °C: $c\approx340$ m/s. It is the $c$ of every worked number in C01, C02 and
C13.""", where="Ch. 1 benchmark", core_id="C01")
remind("C01")
D("D01", "The acoustic relation and the linearised laws", ref="15.10", notes=r"""
**N06 [B]** the expansion of p(ρ) at constant entropy (step 1) · **N07 [B]** the acoustic relation «15.5» with «15.6»
(step 2) · **N11 [B]** linearised continuity «15.8» (steps 3–5) · **N12 [B]** linearised momentum «15.9» (steps 6–8; slip
#2: the book writes $P$ and $\rho_o$ in these lines — we write $p_0$, $\rho_0$ throughout) · **N13 [B]** «15.10» (step 9;
slip #3: the book cites «15.6» where «15.5» does the work) · **N14 [B]** «15.11» (step 10)""")
code(r"""
print(AC.acoustic_density(1.0, c0))                          # (15.5): rho' = p'/c² for p' = 1 Pa -> 8.64e-6 kg/m³
print(AC.acoustic_temperature(1.0, 101325.0, 288.15))        # isentropic T' for p' = 1 Pa -> 8.1e-4 K
print(AC.acoustic_impedance(rho0, c0))                       # rho0 c = 416.9 kg/(m² s): pressure per unit air speed
""", explain=r"""
1. `AC.acoustic_density(p_prime, c)` applies «15.5» backwards: a 1 Pa sound changes the density by 8.64 × 10⁻⁶ kg/m³ —
   **7 parts in a million** of 1.225 kg/m³.
2. `AC.acoustic_temperature` gives the matching temperature change of an isentropic squeeze: under a thousandth of a kelvin.
3. `AC.acoustic_impedance` is $\rho_0c$, the number that will turn a pressure into an air speed in C02.

A 1 Pa sound is loud speech at arm's length.""")
note("N08 [B]", r"""«15.7»: the **condensation** $s=\rho'/\rho_0=p'/\rho_0c^2$ must be small. What was dropped: the relative
error of «15.5» is about $\tfrac12(\gamma-1)s$ (ours).""")
code(r"""
print(AC.condensation(1.0, rho0, c0))                        # s = p'/(rho0 c²) for 1 Pa -> 7.05e-6
for s in (1e-3, 1e-2, 1e-1):                                 # three condensations
    est = AC.linearisation_error(s * rho0 * c0**2, rho0, c0) # ½(gamma-1)s: size of the dropped term relative to the kept one
    rp = s * rho0                                            # the density disturbance rho' [kg/m³] for this condensation
    exact = AC.isentropic_pressure_perturbation(rp, rho0, 101325.0, order="exact")   # p' on the true isentrope
    lin = AC.isentropic_pressure_perturbation(rp, rho0, 101325.0, order=1)           # p' from (15.5)
    print(f"s = {s:g}: estimate {100*est:.2f} %, exact gap {100*(exact - lin)/exact:.2f} %")   # 0.02, 0.2, 2 %
""", explain=r"""
1. The condensation of a 1 Pa sound is 7 × 10⁻⁶: tiny.
2. For three bigger condensations the loop compares the estimate $\tfrac12(\gamma-1)s$ with the true gap between the
   isentrope and the straight line «15.5»: 0.02 %, 0.2 %, 2 % against 0.019 %, 0.20 %, 1.9 %.""")
plotly(r"""
s_vals = [0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]   # 14 condensations, 0.1 included
def curves(s):                                              # the three p'(rho') curves for densities up to s*rho0
    r = np.linspace(-s, s, 120) * rho0                      # rho' from -s rho0 to +s rho0 [kg/m³]
    ex = AC.isentropic_pressure_perturbation(r, rho0, 101325.0, order="exact")   # the true isentrope
    o1 = AC.isentropic_pressure_perturbation(r, rho0, 101325.0, order=1)         # (15.5): the straight line
    o2 = AC.isentropic_pressure_perturbation(r, rho0, 101325.0, order=2)         # with the square term kept
    return {"exact isentrope": (r / rho0, ex / 1e3),        # x: rho'/rho0 [–]; y: p' [kPa]
            "linear (15.5)": (r / rho0, o1 / 1e3),          # the acoustic straight line (same name at every slider position)
            "with the square term": (r / rho0, o2 / 1e3)}
fig = slider_figure(curves, "s", s_vals, xlabel="ρ′/ρ₀ [–]", ylabel="p′ [kPa]",
                    title="Acoustics is the straight line through the origin of the isentrope", yrange=None,
                    active=s_vals.index(0.1))                # open on s = 0.1, the case read below
fig.show()                                                  # the slider moves s; every position was computed now
""")
notes3(see=r"""Three curves of pressure disturbance against relative density disturbance: the true isentrope (purple), the
straight line «15.5» (teal) and the parabola with the square term (orange); hover over a curve to read its values.""",
       read=r"""The slider is the condensation s, from 0.001 to 0.5 (14 hand-picked positions, closer together at small s);
the figure opens at s = 0.1. Drag it down: for s below 0.01 the three curves are one. At s = 0.1 the line is 2 % low: that is a 174 dB sound —
**acoustics is a very good approximation for anything you can listen to**.""",
       change=r"""…γ were 5/3 (a monatomic gas)? The estimate $\tfrac12(\gamma-1)s$ grows from 0.2s to 0.33s: the line
departs from the isentrope faster.""")
P("P356", "decibels and sound pressure level", r"""Loudness is quoted on a logarithmic scale:
$L=20\log_{10}(p_{rms}/20\,\mu\mathrm{Pa})$ dB. Every 20 dB is a factor 10 in pressure. Here it turns "quiet" and "loud"
into pascals.""", code=r"""
p_ref = 20e-6                                 # reference pressure [Pa]: about the faintest audible sound
for L in (60, 94, 140):                       # speech, 1 Pa, a jet engine nearby
    print(L, p_ref*10**(L/20))                # 0.02 Pa, 1.0 Pa, 200 Pa
""")
note("N09 [B]", r"""$\rho_0c^2=\gamma p_0$ = 141.9 kPa for sea-level air. Speech at 60 dB: $s=1.4\times10^{-7}$. A jet engine at
140 dB: $s=1.4\times10^{-3}$ — still small. $p'>0$ is a **compression**, $p'<0$ an **expansion** (rarefaction).""")
code(r"""
print(AC.sound_pressure_level(1.0))                         # 1 Pa rms is 94.0 dB
print(1.4 * 101325.0)                                       # rho0 c² = gamma p0 = 141.9 kPa: the stiffness of air
for L in (60, 140):                                         # speech and a jet engine
    p = 20e-6 * 10**(L / 20)                                # rms pressure [Pa]
    print(L, "dB:", AC.condensation(p, rho0, c0))           # condensations 1.4e-7 and 1.4e-3
""")
note("N10 [B]", r"""Newton assumed the squeeze happens at constant temperature: $c_{\text{Newton}}=\sqrt{RT_0}$, which is
$1/\sqrt\gamma$ of the right value. Laplace's correction: a sound wave passes too quickly for heat to flow, so the squeeze
is isentropic.""")
code(r"""
cN = AC.newton_sound_speed(288.15, air.R)                   # isothermal sound speed sqrt(R T) [m/s]
print(cN, 100 * (1 - cN / c0))                              # 287.6 m/s — 15.5 % low: Laplace was right
""")
P("P354", "an operator in brackets, applied twice", r"""$\left(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\right)$
is an instruction: "take the time derivative, add U times the x-derivative". Squared means *do it twice*, not "square the
result". With constant U two such instructions may be done in either order.""", code=r"""
import sympy as sp                            # symbolic derivatives
x, t, U = sp.symbols('x t U'); p = sp.Function('p')(x, t)   # a pressure field p(x, t)
D = lambda F: sp.diff(F, t) + U*sp.diff(F, x)               # the operator as a function
print(sp.expand(D(D(p))))                     # p_tt + 2 U p_xt + U**2 p_xx
""")
D("D02", "The wave equation for sound in a uniform stream", ref="15.12")
worked("a travelling cosine, by hand", r"""
Take $c$ = 340 m/s, $U$ = 170 m/s ($M$ = 0.5) and try $p'=\cos(x-(U+c)t)$ in «15.12» in one dimension.
1. One application of the bracket: $\partial_t$ gives $+(U+c)\sin$, $U\partial_x$ gives $-U\sin$; sum $c\sin(\cdot)$.
2. Apply it again: $c\,[-(U+c)+U]\cos=-c^2\cos$.
3. Divide by $c^2$: $-\cos$.
4. The second x-derivative is $-\cos$.
5. Left side: $-\cos-(-\cos)=0$ ✓. The shape rides at $U+c$ = 510 m/s. Try the speed 400 m/s instead and the bracket leaves a
   remainder: only $U\pm c$ work.""")
code(r"""
out = ch15.linearised_acoustics_sympy()                     # the expansions of D01 done by sympy
print("dropped from continuity:", out["continuity_dropped"])  # every term with two primes
print("dropped from momentum:  ", out["momentum_dropped"])    # every term with two primes
p_good = lambda x, t: np.exp(-(x - (170.0 + c0) * t)**2)     # a Gaussian riding at U + c [shape in metres]
p_bad = lambda x, t: np.exp(-(x - 400.0 * t)**2)             # the same shape riding at the wrong speed 400 m/s
print(AC.convected_wave_residual(p_good, 0.3, 0.001, c0, U=170.0))   # (15.12) left side: ≈ 0
print(AC.convected_wave_residual(p_bad, 0.3, 0.001, c0, U=170.0))    # of order one: not a solution
""", explain=r"""
1. `ch15.linearised_acoustics_sympy()` puts the split «15.4» into continuity and Euler's equation symbolically and returns,
   among others, the terms that D01 drops: each is a product of two primed quantities.
2. `p_good` is a Gaussian bump that rides at $U+c$ = 510 m/s; `p_bad` rides at 400 m/s.
3. `AC.convected_wave_residual` evaluates the left side of «15.12» by central differences at one point (x = 0.3 m,
   t = 1 ms): almost zero for the right speed, of order one for the wrong one.
4. So the equation picks out the speeds; the shape is free.""")
scratch(r"""
c, U = 1.0, 0.5                                             # units in which c = 1 (so M = 0.5)
p = lambda x, t: np.exp(-(x - (U + c)*t)**2)                # a Gaussian riding at U + c
x0, t0, h = 0.3, 0.2, 1e-3                                  # the point and the step of the differences
ptt = (p(x0, t0 + h) - 2*p(x0, t0) + p(x0, t0 - h)) / h**2  # second difference in t
pxx = (p(x0 + h, t0) - 2*p(x0, t0) + p(x0 - h, t0)) / h**2  # second difference in x
pxt = (p(x0 + h, t0 + h) - p(x0 + h, t0 - h) - p(x0 - h, t0 + h) + p(x0 - h, t0 - h)) / (4*h*h)   # mixed difference
lhs = (ptt + 2*U*pxt + U**2*pxx) / c**2 - pxx              # (15.12) in one dimension, operator expanded as in primer P354
lib = AC.convected_wave_residual(p, x0, t0, c, U=U, h=h)    # the tested function
print(lhs, lib)                                             # both ≈ 0 (discretisation error only)
assert abs(lhs) < 1e-5 and np.isclose(lhs, lib, atol=1e-5)  # same numbers -> the library does exactly this
""", explain=r"""
**From scratch.** Three central differences (the Ch. 10 stencils): $(p_{i+1}-2p_i+p_{i-1})/h^2$ for a second derivative and a
four-point stencil for the mixed one. Put together as in primer P354 — $\partial_{tt}+2U\partial_{xt}+U^2\partial_{xx}$ — they give
the left side of «15.12». It vanishes for a shape riding at $U+c$, and agrees with `AC.convected_wave_residual`.""")
figure(r"""
U, c = 0.5, 1.0                                             # M = 0.5 in units where c = 1
speeds = np.linspace(-2.0, 3.0, 501)                        # trial speeds s from -2c to 3c
t0 = 0.2                                                    # the instant at which we test
res = [AC.convected_wave_residual(lambda x, t, s=s: np.exp(-(x - s*t)**2), s*t0, t0, c, U=U, h=1e-3)
       for s in speeds]                                     # residual at the peak of a bump riding at speed s
fig, ax = plt.subplots(figsize=(6.8, 3.4))                  # one panel
ax.plot(speeds, np.abs(res), color=COLORS["accent"])        # |left side of (15.12)| against the trial speed
for s0, lab in ((U - c, "U − c"), (U + c, "U + c")):        # the two zeros the theory predicts
    ax.axvline(s0, color=COLORS["muted"], ls="--", lw=1)    # dashed line at each
    ax.text(s0 + 0.08, 6.0, lab, ha="left", color=COLORS["muted"],
            bbox=dict(fc="white", ec="none", pad=1))        # its label, just right of the line (clear of curve and line)
ax.set(xlabel="speed of the trial shape s / c [–]", ylabel="|residual of (15.12)| [–]",
       title="Only two speeds solve the equation: U − c and U + c")   # units and message
savefig(fig, "ch15", "c01_residual_vs_speed"); plt.show()  # save to outputs/ch15/ and display
""", see=r"""The size of the left side of «15.12» for a Gaussian bump riding at speed s, against s, for a stream at M = 0.5
(in units where c = 1). The residual itself is a parabola in s that crosses zero at two speeds; we plot its absolute value,
so the curve is a W: two V-shaped dips that touch zero at exactly those two speeds.""",
    read=r"""The zeros sit at $U-c=-0.5$ (the shape travels upstream, slowly) and $U+c=1.5$ (downstream, fast). Any other
speed leaves a residual — the shape would have to change as it moves, and «15.12» does not allow that.""",
    change=r"""…U rose above c? Both zeros would be positive: no shape can travel upstream. C02 turns that into "no news
upstream".""")
whatif(r"""…the amplitude were not small? The dropped products come back, compressions travel faster than expansions,
and the wave steepens into a shock — N22 below, then §15.6.""")

# ---------------------------------------------------------------------------------------------------------------------
# C02 — d'Alembert's solution
# ---------------------------------------------------------------------------------------------------------------------
core("C02", r"Two shapes sliding apart: $p'(x,t)=f(x-ct)+g(x+ct)$ (15.15)",
     "What does the most general plane sound wave look like — and which way does the air itself move?")
problem(r"""
Pluck a long string, or pop a balloon in a long pipe. The bump you made does not sit there and does not spread out: it
splits into two copies that run off in opposite directions, each keeping its shape. That is all a plane sound wave can do.
The air itself only shuffles: forward under a compression that runs forward, backward under one that runs backward.""")
idea(r"""
t = 0:            __/‾‾\__
t > 0:   _/‾\_   ⟵ g        f ⟶   _/‾\_        each half as high
         air moves ⟵                air moves ⟶
""", r"""Each half carries half the bump; under the right-going half the air moves right, under the left-going half it
moves left.""")
note("N15 [B]", r"""«15.13»: «15.12» with $U_j=0$ in one space dimension — the second-order relative of the first-order wave
equation of Ch. 7. `AC.wave_equation_residual` evaluates its left side.""")
remind("C02")
P("P357", "the second-order wave equation and its two characteristic coordinates", r"""A first-order wave equation has
one family of lines along which things are carried (Ch. 7, P174). The second-order one has two: $\xi=x-ct$ and $\zeta=x+ct$.
Rewritten in ξ and ζ (chain rule, twice) it becomes $\partial^2p'/\partial\xi\,\partial\zeta=0$ — "the ξ-slope does not
depend on ζ" — whose solutions are a function of ξ plus a function of ζ.""", code=r"""
import sympy as sp                                       # symbolic algebra
x, t, c = sp.symbols('x t c', positive=True)             # position, time, sound speed
f, g = sp.Function('f'), sp.Function('g')                # two arbitrary shapes
p = f(x - c*t) + g(x + c*t)                              # any two shapes, one moving each way
print(sp.simplify(sp.diff(p, t, 2)/c**2 - sp.diff(p, x, 2)))   # 0: it solves the wave equation (15.13)
""")
D("D03", "d'Alembert's solution", ref="15.15")
note("N16 [B]", r"""«15.14» — the momentum law «15.11» (with $U_j=0$) integrated in time. > ⚠️ Common confusion (trap T19):
the integral is indefinite; its "constant" is a function of x, fixed by "the air was at rest before the wave arrived".""")
code(r"""
p_fn = lambda x, t: 50.0 * np.exp(-(x - c0 * t)**2)          # a right-going 50 Pa bump (width 1 m) riding at c0
u_quad = AC.velocity_from_pressure(p_fn, 0.5, 0.0015, rho0, t0=-0.02)   # (15.14): time integral from before the wave arrived
u_formula = AC.plane_wave_velocity(0.5, 0.0015, c0, rho0, f=lambda s: 50.0 * np.exp(-s**2))   # (15.16) with g = 0
print(u_quad, u_formula)                                     # the same air speed [m/s], two routes
""", explain=r"""
1. `p_fn` is a right-going pressure bump.
2. `AC.velocity_from_pressure` integrates «15.14» numerically in time, starting at t₀ = −20 ms when the bump was 7 m away
   and the air at x = 0.5 m was still at rest (that fixes the "constant").
3. `AC.plane_wave_velocity` uses the closed form «15.16» that D04 derives next. Both give ≈ 0.12 m/s.""")
D("D04", "The particle velocity of a plane sound wave", ref="15.16", notes=r"""**N17 [B]** the particle velocity «15.16»""",
  after=r"""
| Wave | Air velocity | In words |
|---|---|---|
| right-going f | $u'=+p'/\rho_0c$ | a compression pushes the air the way the wave goes |
| left-going g | $u'=-p'/\rho_0c$ | a compression pushes the air the other way — still the way *its* wave goes |

$\rho_0c$ is the **acoustic impedance** (our name; 416.9 kg/(m² s) for sea-level air). And $u'/c=p'/\rho_0c^2$ is the
condensation: **the air moves far slower than the wave**.""")
worked("a bump released at rest", r"""
$c$ = 340 m/s, $\rho_0$ = 1.2 kg/m³. At t = 0 the pressure is a bump F(x) of height 100 Pa and the air is at rest, so
$\partial p'/\partial t=G=0$.
1. The initial-value formulas «ex15.1» with G = 0 give $f=g=\tfrac12F$: two bumps of 50 Pa.
2. After 1 s their centres are at x = +340 m and x = −340 m.
3. Under the right-going one $u'=50/(1.2\times340)$ = +0.12 m/s; under the left-going one −0.12 m/s.
4. At t = 0 the two velocities cancel: the air was at rest ✓.
5. Wave speed 340 m/s, air speed 0.12 m/s: a ratio of 3.6 × 10⁻⁴ — the condensation.

This is **N21 [B]**: the initial-value form, stated; result given, not derived.""")
code(r"""
F = AC.pulse("gaussian", width=1.0, amplitude=100.0)          # the initial bump F(x) = 100 exp(-x²) [Pa]
f, g = AC.dalembert_fg(F, lambda s: 0.0 * s, c0)              # Exercise-15.1 formulas with G = 0: the two halves
x = np.linspace(-8.0, 8.0, 801)                               # positions [m]
for t in (0.0, 0.01):                                         # at release and 10 ms later (the halves 3.4 m from the centre)
    p = AC.dalembert_initial_value(x, t, c0, F)               # (15.15): p' = f(x - ct) + g(x + ct)
    u = AC.plane_wave_velocity(x, t, c0, rho0, f=f, g=g)      # (15.16): u' = (f - g)/(rho0 c)
    print(f"t = {1e3*t:4.1f} ms: max p' = {p.max():6.2f} Pa, max |u'| = {np.abs(u).max():.4f} m/s")
print(np.allclose(AC.dalembert_initial_value(x, 0.0, c0, F), F(x)), np.allclose(AC.plane_wave_velocity(x, 0.0, c0, rho0, f=f, g=g), 0.0))   # True True
""", explain=r"""
1. `AC.pulse` builds the bump as a function (a `lambda`-like callable: give it positions, it returns pressures).
2. `AC.dalembert_fg` applies the initial-value formulas «ex15.1»: with no initial rate of change, f and g are each half of
   F.
3. At t = 0 the pressure is the bump itself (100 Pa) and the air is at rest; 10 ms later two 50 Pa halves have separated
   and the air under each moves at 50/416.9 = 0.120 m/s.
4. The last line checks the two starting conditions: True, True.""")
scratch(r"""
c, t = c0, 0.01                                              # sound speed [m/s] and a time [s]
shape = lambda s: 50.0 * np.exp(-s**2)                       # each half: 50 Pa Gaussian
p_mine = shape(x - c*t) + shape(x + c*t)                     # (15.15) written out
u_mine = (shape(x - c*t) - shape(x + c*t)) / (rho0 * c)      # (15.16) written out
assert np.allclose(p_mine, AC.dalembert(x, t, c, f=shape, g=shape))                 # same pressure as the library
assert np.allclose(u_mine, AC.plane_wave_velocity(x, t, c, rho0, f=shape, g=shape)) # same air speed as the library
print("from scratch = AC.dalembert and AC.plane_wave_velocity")
""", explain=r"""**From scratch.** «15.15» and «15.16» are one line each; the assertions confirm that the library does exactly this.""")
animation(r"""
x = np.linspace(-8.0, 8.0, 400)                              # positions [m]
n_frames = 10 if FAST else 16                                # frames (fewer when FAST; each two-panel frame costs 1–3 s to draw)
times = np.linspace(0.0, 0.018, n_frames)                    # 0 to 18 ms: the halves travel 6 m
right = lambda s: 100.0 * np.exp(-(s + 4.0)**2)              # top panel: a right-going 100 Pa pulse starting at x = -4 m
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.8, 4.6), sharex=True)   # two panels, one clock
(pt,) = a1.plot(x, right(x), color=COLORS["orange"], label="p′ [Pa]")             # pressure (orange)
(ut,) = a1.plot(x, rho0 * c0 * 0 * x, color=COLORS["blue"], ls="--", label="ρ₀c·u′ [Pa]")   # air speed × rho0 c (blue)
a1.set(ylim=(-60, 110), ylabel="[Pa]", title="A right-going pulse: the air moves with the compression")
a1.legend(loc="upper right", fontsize=8)                    # legend
F0 = lambda s: 100.0 * np.exp(-s**2)                         # bottom panel: a bump released at rest
(pb,) = a2.plot(x, F0(x), color=COLORS["orange"])            # pressure
(ub,) = a2.plot(x, 0 * x, color=COLORS["blue"], ls="--")     # rho0 c u'
a2.set(ylim=(-60, 110), xlabel="x [m]", ylabel="[Pa]", title="A bump released at rest splits into two half-height pulses")
def update(i):                                               # draw frame i
    t = times[i]                                             # current time [s]
    pt.set_ydata(right(x - c0 * t))                          # (15.15) with g = 0
    ut.set_ydata(rho0 * c0 * AC.plane_wave_velocity(x, t, c0, rho0, f=right))      # (15.16): equals p' for a right-going wave
    pb.set_ydata(AC.dalembert(x, t, c0, f=lambda s: 0.5 * F0(s), g=lambda s: 0.5 * F0(s)))   # two halves
    ub.set_ydata(rho0 * c0 * AC.plane_wave_velocity(x, t, c0, rho0, f=lambda s: 0.5 * F0(s), g=lambda s: 0.5 * F0(s)))
    return pt, ut, pb, ub
show_animation(animate(update, frames=n_frames, fig=fig, interval=60), player="video")   # smooth motion: video player
""", explain=r"""
**N18 [B]** (animation A1). Top: one right-going pulse; the blue dashed curve is the air speed multiplied by $\rho_0c$, so that
it is in pascals — it lies exactly on the pressure. Bottom: the bump released at rest; the two halves carry air speeds of
opposite sign, which cancel at t = 0.""")
notes3(see=r"""Top: a pressure pulse (orange) sliding right, with $\rho_0c\,u'$ (blue dashed) on top of it. Bottom: a 100 Pa
bump splitting into two 50 Pa pulses; the blue curve is positive under the right-going half and negative under the left-going
one.""",
       read=r"""Neither shape changes as it moves — that is «15.15». In the top panel the two curves coincide: $u'=p'/\rho_0c$.
In the bottom panel the air moves outward from the centre under both halves: each half pushes the air its own way.""",
       change=r"""…the bump had been released with the air already moving right at $u'=F/\rho_0c$? Then $G\neq0$ and only a
right-going pulse would appear (all of it in f).""")
note("N19 [B]", r"""In a stream of speed U the two shapes ride at U + c and U − c: «15.17».""")
D("D05", "Sound in a uniform stream", ref="15.17")
code(r"""
print(AC.wave_speeds(0.5 * c0, c0))                          # M = 0.5: (U + c, U - c) = (510.4, -170.1) m/s
print(AC.wave_speeds(2.0 * c0, c0))                          # M = 2: (1020.9, +340.3) m/s — both downstream
print(AC.upstream_signal_possible(0.5 * c0, c0), AC.upstream_signal_possible(2.0 * c0, c0))   # True, False
""")
confusion(r"""(trap T19) "Downstream" presumes U > 0; with the stream flowing in −x the roles of f and g swap.""")
note("N20 [B]", r"""**Range of influence.** For U < c one shape still travels upstream: a point feels what happens on both
sides. For U > c both travel downstream: a point feels only what is upstream of it. That is why a choked nozzle stops
listening (C09) and why you hear a supersonic aircraft only after it has passed (C13).""")
figure(r"""
t = np.linspace(0.0, 0.01, 50)                              # times [s]
x0 = np.linspace(-6.0, 6.0, 13)                             # starting points of the lines [m]
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=True)   # two panels: M = 0.5 and M = 2
for ax, M in zip(axes, (0.5, 2.0)):                         # one Mach number per panel
    ch = AC.xt_characteristics(M * c0, c0, x0, t)           # the two families x = x0 + (U ± c) t
    for k in range(x0.size):                                # every starting point
        ax.plot(ch["x_plus"][k], 1e3 * t, color=COLORS["orange"], lw=0.8)    # f family, slope U + c
        ax.plot(ch["x_minus"][k], 1e3 * t, color=COLORS["blue"], lw=0.8)     # g family, slope U - c
    one = AC.xt_characteristics(M * c0, c0, 0.0, t)         # the two lines through the origin
    ax.fill_betweenx(1e3 * t, one["x_minus"][0], one["x_plus"][0], color=COLORS["accent"], alpha=0.18)   # its range of influence
    ax.set(xlim=(-6, 12), xlabel="x [m]", title=f"M = {M}: speeds {M*c0 + c0:.0f} and {M*c0 - c0:+.0f} m/s")   # labels
axes[0].set_ylabel("t [ms]")                                # time axis
fig.suptitle("Where the news of a disturbance at x = 0, t = 0 can reach (shaded)", fontsize=11)   # the message
savefig(fig, "ch15", "c02_xt_characteristics"); plt.show()  # save and display
""", see=r"""Two x–t diagrams (time upward). Orange lines: the f family, slope U + c. Blue lines: the g family, slope U − c.
The shaded wedge is everything a disturbance made at the origin can reach.""",
    read=r"""At M = 0.5 the blue lines lean left: news travels upstream at 170 m/s and the wedge opens to both sides. At
M = 2 both families lean right: the wedge lies entirely downstream — a point upstream of the disturbance never hears it
(the **domain of dependence** of Ch. 10, read the other way round).""",
    change=r"""…M were exactly 1? The blue lines would be vertical: the upstream-running wave would stand still at its
starting point — the picture of a sonic throat (C09).""")
explainer("sound_pulse_dalembert", "What is a sound pulse made of?",
          r"""you release the bump, scrub time and watch f and g separate with tracer particles moving opposite ways; then you
drag the stream through U = c and one family of characteristics tips over — the change of regime is something you do, not
something you read.""",
          ["Release a bump at rest and step time: each half is half as high.",
           "Click a point under each half: the air moves opposite ways.",
           "Drag the stream Mach number through 1 and watch the x–t view: the upstream line stands up, then leans downstream."])
note("N22 [B]", r"""**Steepening (the book argues in words; the formulas are ours).** In a compression the gas is hotter, so c
is larger, and the gas itself moves forward: each wavelet travels at «o_steep». Crests catch up with troughs; a sine becomes
a sawtooth; the limit is a shock. Result given, not derived.""")
recap_md("The shallow-water simple wave", "Ch. 7 N81", r"""The same steepening happens to a long water wave: its crest moves
at a speed set by the local depth. It is the γ = 2 twin of this gas wave (`ch07.simple_wave_evolve`).""")
P("P358", "Riemann invariant", r"""A combination of the fluid speed and the sound speed that keeps its value along a moving
wavelet. In a wave going one way through still gas, $u-2c/(\gamma-1)$ has the same value everywhere, so
$c=c_0+\tfrac12(\gamma-1)u$: each bit of the wave knows its own speed.""", code=r"""
g, c0_ = 1.4, 340.0                           # air
u = np.array([-5.0, 0.0, 5.0])                # fluid speed in trough, node, crest [m/s]
c = c0_ + 0.5*(g - 1)*u                       # local sound speed from the invariant
print(u - 2*c/(g - 1), u + c)                 # the invariant is constant; wavelet speeds 334, 340, 346 m/s
""")
code(r"""
L = 150.0                                                    # a loud 1 kHz tone [dB]
p_peak = np.sqrt(2) * 20e-6 * 10**(L / 20)                   # peak pressure of a sine [Pa]: rms × sqrt(2) -> 894 Pa
u_a = p_peak / (rho0 * c0)                                   # (15.16): air-speed amplitude -> 2.15 m/s
lam = c0 / 1000.0                                            # wavelength [m]: c / frequency
x = np.linspace(0.0, lam, 400)                               # one wavelength
u0 = u_a * np.sin(2 * np.pi * x / lam)                       # the initial air speed [m/s]
tb = AC.breaking_time(u0, x)                                 # first time two wavelets meet [s]
print(f"peak {p_peak:.0f} Pa, u_a = {u_a:.2f} m/s, t_break = {1e3*tb:.1f} ms, distance {c0*tb:.2f} m = {c0*tb/lam:.0f} wavelengths")
print(AC.simple_wave_speed(np.array([-u_a, 0.0, u_a]), c0))  # wavelet speeds in trough, node, crest [m/s]
""", explain=r"""
1. 150 dB is a peak pressure of 894 Pa; by «15.16» the air moves at 2.15 m/s.
2. `AC.breaking_time` finds where the speed falls fastest with x (the steepest compression) and returns
   $t_{\text{break}}=2/[(\gamma+1)\max(-\partial u_0/\partial x)]$: 21 ms, after travelling 7.2 m — about 21 wavelengths.
3. `AC.simple_wave_speed` gives the three wavelet speeds: the crest runs 5 m/s faster than the trough.""")
animation(r"""
n_frames = 7 if FAST else 13                                 # frames from 0 to the breaking time (13: stops at exactly ¼, ½ and 1)
times = np.linspace(0.0, tb, n_frames)                       # [s]
fig, ax = plt.subplots(figsize=(6.8, 3.4))                   # one panel
ax.plot(1e3 * x, u0, color=COLORS["muted"], ls="--", lw=1, label="acoustic (shape kept)")   # the ghost
(line,) = ax.plot(1e3 * x, u0, color=COLORS["accent"], lw=2, label="simple wave")            # the steepening wave
ax.set(xlim=(0, 1e3 * lam), ylim=(-1.2 * u_a, 1.2 * u_a), xlabel="x − c₀t [mm] (riding with the sound)",
       ylabel="air speed u [m/s]")                           # fixed axes
ax.legend(loc="lower left", fontsize=8)                      # legend
def update(i):                                               # frame i
    s = AC.simple_wave_evolve(u0, x, times[i], c0)           # each wavelet moved at c0 + (gamma+1)/2 u0
    line.set_xdata(1e3 * (s["x_points"] - c0 * times[i]))    # position seen by an observer riding at c0
    ax.set_title(f"t = {times[i]/tb:.2f} × t_break: compressions catch up")   # time as a fraction of the breaking time
    return (line,)
show_animation(animate(update, frames=n_frames, fig=fig, interval=120), player="frames")   # step through: stop at ¼, ½, 1
""", explain=r"""**Animation A2.** Each point of the wave keeps its air speed and moves at $c_0+\tfrac{\gamma+1}{2}u$; we watch from
a frame riding at $c_0$, so the acoustic wave (grey dashed) stands still. Use the frame buttons to stop at ¼, ½ and 1 × the
breaking time.""")
notes3(see=r"""One wavelength of a 150 dB tone, air speed against position, in a frame that rides with the sound. The grey
dashed sine is what acoustics predicts (no change); the purple curve is the simple wave.""",
       read=r"""Points with u > 0 drift right, points with u < 0 drift left. Where the speed falls with x (the compression,
in the middle of the frame) the curve tilts until it is vertical at $t=t_{\text{break}}$: a shock forms there.""",
       change=r"""…you looked at the expansion half? It flattens instead: only compressions become shocks. And a quieter tone
(say 120 dB, about 32 times smaller in pressure) breaks about 32 times later — after some 660 wavelengths.""")
note("N114 [C]", r"""The full nonlinear (Euler) equations can carry a wave *through* breaking, which the simple-wave
formula cannot; Ch. 10's MacCormack scheme (C10) is the tool. Not built here (it is part of the backup explainer).""")
whatif(r"""…the wave were confined in a duct whose area changes? Then we stop following the wave and write the three
conservation laws between two stations — C03.""")

# =====================================================================================================================
# A.3 §15.3 Basic Equations for One-Dimensional Flow — C03
# =====================================================================================================================
nb.section("15.3", "Basic Equations for One-Dimensional Flow", intro="**What is this section about?** Mass, energy and "
           "momentum written between two stations of a duct. The rest of the chapter is these three lines with different "
           "terms switched on.")
core("C03", r"The three duct laws: $h_2+\tfrac12u_2^2-h_1-\tfrac12u_1^2=Q$ (15.22)",
     "What must be the same at two cross-sections of a duct — and which of the three statements survives friction and shocks?")
problem(r"""
A hair dryer, a jet engine intake, a rocket nozzle: gas goes in one end and out the other. We do not want the flow at every
point — only one speed, pressure, density and temperature per cross-section. Three bookkeeping statements tie the outlet to
the inlet: what goes in comes out (mass); what the gas carries in plus the heat added equals what it carries out (energy);
and the change of momentum equals the push of the pressures and the walls.""")
idea(r"""
station 1 │→ u₁ p₁ ρ₁ T₁   A₁ ⟩═══════════⟨ A₂   u₂ p₂ ρ₂ T₂ →│ station 2
                walls: force F on the fluid · heat Q per kg added
""", r"""
| Law | What crosses the stations | Switched on in |
|---|---|---|
| mass | ṁ = ρuA | C05, C06 |
| energy | $h+\tfrac12u^2$ per kg, plus the heat Q | C04, C07, C11 |
| momentum | ṁu, plus the pressure forces pA and the wall force F | C07, C11 |""")
note("N23 [B]", r"""**Hypotheses of one-dimensional duct flow:** (1) a straight centreline; (2) the area A(x) varies slowly;
(3) one value of u, p, ρ, T per cross-section; (4) steady flow.""")
figure(r"""
fig = draw_duct_control_volume()                             # our drawing: control volume between stations 1 and 2, normals, F and Q
plt.show()                                                   # display it
""", see=r"""A duct whose area changes slowly, with the control volume between stations 1 and 2 (dashed), the outward normals
on the two faces, the wall force F and the heat Q.""",
    read=r"""The outward normal at station 1 points *against* the flow and at station 2 *with* it — that is where the
signs of the inlet and outlet terms in D06 and D07 come from. Nothing crosses the side wall except heat.""",
    change=r"""…the area changed abruptly (a step)? Then hypothesis (3) fails near the step: the flow separates and one value
per cross-section is no longer enough.""")
recap_md("Control volume and the Reynolds transport theorem", "Ch. 3 C15, Ch. 4 C01", r"""A fixed region of space through
whose surface the fluid flows; the rate of change of what is inside plus the net outflow through the surface equals what is
produced inside. In steady flow only the surface fluxes remain.""")
recap("R09", "Mass conservation on a control volume", r"""What flows in flows out: «15.18». **New here:** on this control
volume the integral law is one line.""", where="Ch. 4 §4.2 (C01)", core_id="C03")
code(r"""
print(GD.mass_flow(1.0, 100.0, 0.02))                        # (15.18): rho u A = 1.0 × 100 × 0.02 = 2 kg/s
""")
recap("R10", "The steady energy integral", r"""«15.20»: the flux of internal plus kinetic energy out of the surface equals
the rate of work of the surface stresses plus the heat conducted in. > ⚠️ Common confusion (trap T2): here $A^*$ is the
control surface, not the sonic throat area of «15.31».""", where="Ch. 4 §4.8 (C10)", core_id="C03")
remind("C03")
recap_md("Enthalpy", "Ch. 1", r"""$h=e+p/\rho$; for a perfect gas $h=C_pT$ and $e=C_vT$.""")
D("D06", "The energy equation of duct flow", ref="15.22", notes=r"""**N25 [B]** the duct energy law «15.21» → «15.22»
(steps 1–5). Face by face:

| Face | Outward normal | What crosses |
|---|---|---|
| inlet $A_1$ | against the flow ($u_jn_j=-u_1$) | energy in with ṁ; pressure work done *on* the fluid |
| outlet $A_2$ | with the flow ($u_jn_j=+u_2$) | energy out with ṁ; pressure work done *by* the fluid |
| wall | across the flow ($u_jn_j=0$) | heat only (no slip: no work) |

**N26 [B]** flow work $p/\rho$ (step 6)""")
note("N26 [B], continued", r"""**Flow work.** Pushing one kilogram across a face against the pressure costs $upA/\dot m=p/\rho$ of work.
That is why $h=e+p/\rho$ replaces e in a flowing system.""")
code(r"""
T = 293.15                                                   # room temperature [K]
print(air.R * T, air.cv * T, air.cp * T)                     # p/rho = RT = 84.1, e = 210.4, h = 294.5 kJ/kg (×1000)
print((air.R * T) / (air.cp * T), (1.4 - 1) / 1.4)           # flow work is 29 % of h: (gamma - 1)/gamma
""")
recap("R11", "The steady momentum integral", r"""«15.23»: momentum flux out = force of the surface stresses.""",
      where="Ch. 4 §4.4 (C03)", core_id="C03")
D("D07", "Mass and momentum in a duct", ref="15.25", notes=r"""**N24 [B]** «15.19» (steps 1–2) · **N27 [B]** «15.24»
(step 3) · **N28 [B]** «15.25», the wall force split into the pressure on a sloping wall and friction (steps 4–8) · **R12**
«15.26» (step 9)""")
note("N27 [B], continued", r"""> ⚠️ Common confusion (trap T11): three F's.

| Symbol | Where | Meaning | Sign |
|---|---|---|---|
| F | «15.24» | force of the walls **on the fluid** | positive along +x |
| $F_f$ | «15.25» | friction force per unit length | positive = *retarding* |
| F | $f=F/(p_1A)$ in «15.45» | total friction force | positive = *retarding* |

**N24** in one number: 1 % more area at fixed density gives 1 % less speed.""")
code(r"""
print(GD.influence_coefficients(0.0, dA_over_A=0.01)["du_over_u"])   # (15.19) at M = 0 (fixed density): -0.01
""")
recap("R12", "The one-dimensional Euler equation", r"""With no friction, «15.26» — Euler's equation along the duct.""",
      where="Ch. 4 (C07)", core_id="C03")
code(r"""
eqs = ch15.duct_equations_sympy()                            # sympy re-derives (15.19), (15.25), (15.26) from (15.18), (15.24)
for k in ("eq_15_19", "eq_15_25", "eq_15_26"):               # the three derived lines
    print(k, ":", eqs[k])                                    # each printed as a sympy expression (= 0)
""")
recap("R13", "Euler + Gibbs = the energy Bernoulli equation", r"""With $dS=0$ the Gibbs relation $T\,dS=dh-dp/\rho$ gives
$dh=dp/\rho$, so «15.26» becomes $d(h+\tfrac12u^2)=0$: the Q = 0 case of «15.22». > ⚠️ Common confusion (trap T17):
**adiabatic** = no heat (Q = 0); **reversible** = no friction; **isentropic** = both. «15.22» needs only "adiabatic"; this
route also needs "frictionless". Result given, not derived.""", where="Ch. 1 (D10); Ch. 4 §4.9", core_id="C03")
worked("three laws, easy numbers", r"""
Air ($C_p$ ≈ 1000 J/(kg K)) enters a short insulated duct at $u_1$ = 100 m/s, $T_1$ = 300 K, $\rho_1$ = 1.0 kg/m³ through
$A_1$ = 0.02 m² and leaves at $u_2$ = 300 m/s.
1. Mass «15.18»: ṁ = 1.0 × 100 × 0.02 = 2 kg/s.
2. Energy «15.22» with Q = 0: $C_p(T_2-T_1)=\tfrac12(u_1^2-u_2^2)$ = ½(10⁴ − 9 × 10⁴) = −40 000 J/kg, so $T_2$ = 260 K: the
   gas cooled by 40 K to pay for its speed.
3. Momentum «15.24»: ṁ($u_2-u_1$) = 2 × 200 = 400 N must be supplied by the pressure difference on the two faces plus the
   wall force F.
4. Nothing here asked whether the flow was frictionless — the energy line does not care.""")
code(r"""
p0, T0, A1, A2 = noz["p0"], noz["T0"], 2.0 * noz["A_t"], noz["A_t"]   # supply state [Pa, K]; station 1 (20 cm²) and the throat (10 cm²)
M1 = GD.mach_from_area_ratio(A1 / noz["A_t"], branch="subsonic")       # (15.31) inverted on the subsonic branch -> 0.306
states = []                                                  # the two stations
for M, A in ((M1, A1), (1.0, A2)):                           # station 1, then the sonic throat
    st = GD.isentropic_state(M)                              # (15.28)-(15.30): ratios to the stagnation state
    p, T = p0 * st["p_over_p0"], T0 * st["T_over_T0"]        # static pressure [Pa] and temperature [K]
    rho, u = p / (air.R * T), M * air.c(T)                   # gas law; u = M c
    states.append(dict(p=p, rho=rho, u=u, T=T, A=A))         # one station as a dict
s1, s2 = states                                              # unpack
for s in states:                                             # print both stations
    print(f"p = {s['p']/1e3:6.1f} kPa, T = {s['T']:5.1f} K, u = {s['u']:5.1f} m/s, mdot = {GD.mass_flow(s['rho'], s['u'], s['A']):.3f} kg/s, "
          f"h0 = {GD.total_enthalpy(s['T'], s['u'], air.cp)/1e3:.1f} kJ/kg")
state1 = dict(p=s1["p"], rho=s1["rho"], u=s1["u"]); state2 = dict(p=s2["p"], rho=s2["rho"], u=s2["u"])   # what the residual functions take
print("energy residual (15.22):", GD.energy_residual(state1, state2))                     # ≈ 0 J/kg
print("wall force on the fluid (15.24):", GD.wall_force(state1, state2, A1, A2), "N")     # -426 N
""", explain=r"""
1. Station 1 is where the area is twice the throat area; `GD.mach_from_area_ratio(…, branch="subsonic")` gives M = 0.306 there
   (C05 explains the two branches). The throat is sonic.
2. For each station `GD.isentropic_state` returns the ratios to the stagnation state (C04); times $p_0$, $T_0$ they give
   p and T; the gas law gives ρ; u = Mc.
3. A station is stored as a Python `dict` with keys `p`, `rho`, `u` — the form the residual functions take.
4. The printed lines: 468.6 kPa, 392.7 K, 121.5 m/s at station 1; 264.1 kPa, 333.3 K, 366.0 m/s at the throat; the **same
   mass flow** 1.010 kg/s and the **same total enthalpy** 401.9 kJ/kg at both.
5. `GD.wall_force` solves «15.24» for F: **−426 N**. The converging walls push the gas *backward* while it accelerates — the
   pressure difference does all the pushing.""")
scratch(r"""
cp = air.cp                                                  # [J/(kg K)]
mdot1 = s1["rho"] * s1["u"] * A1; mdot2 = s2["rho"] * s2["u"] * A2         # (15.18) at both stations [kg/s]
dh0 = cp * s2["T"] + 0.5 * s2["u"]**2 - cp * s1["T"] - 0.5 * s1["u"]**2    # (15.22) with Q = 0: should vanish [J/kg]
F_mine = mdot1 * (s2["u"] - s1["u"]) - (s1["p"] * A1 - s2["p"] * A2)       # (15.24) solved for F [N]
assert np.allclose([mdot1, F_mine], [mdot2, GD.wall_force(state1, state2, A1, A2)])   # mass conserved; same F as the library
assert abs(dh0) < 1e-6 * cp * s1["T"]                        # energy conserved to round-off
print(f"mdot {mdot1:.4f} = {mdot2:.4f} kg/s, F = {F_mine:.1f} N, change of h0 = {dh0:.2e} J/kg")
""", explain=r"""**From scratch.** The three laws are three lines of arithmetic. Mass and energy balance by themselves (the
isentropic state was built to satisfy them); momentum then *tells* us the wall force, and it equals the library's.""")
figure(r"""
mdot = s1["rho"] * s1["u"] * A1                              # mass flow [kg/s]
F = GD.wall_force(state1, state2, A1, A2)                    # wall force on the fluid [N]
fig, (a, b, c) = plt.subplots(1, 3, figsize=(9.6, 3.4))      # one panel per law
a.bar(["in", "out"], [mdot, s2["rho"] * s2["u"] * A2], color=COLORS["teal"])   # mass flux in and out
a.set(title="mass [kg/s]")                                   # the units are in the title
h = [air.cp * s["T"] / 1e3 for s in states]; k = [0.5 * s["u"]**2 / 1e3 for s in states]   # h and ½u² [kJ/kg]
b.bar(["in", "out"], h, color=COLORS["rose"], label="h")     # enthalpy (rose)
b.bar(["in", "out"], k, bottom=h, color=COLORS["blue"], label="½u²")   # kinetic energy on top (blue)
b.axhline(h[0] + k[0], color=COLORS["amber"], lw=2, label="h₀ (conserved)")   # the total: flat
b.set(title="energy [kJ/kg]"); b.legend(fontsize=8, loc="lower right")
top_in = mdot * s1["u"] + s1["p"] * A1                       # momentum flux plus pressure force at station 1 [N]
c.bar(["in"], [mdot * s1["u"]], color=COLORS["blue"], label="ṁu")                   # momentum flux in
c.bar(["in"], [s1["p"] * A1], bottom=[mdot * s1["u"]], color=COLORS["orange"], label="pA")   # pressure force on face 1
c.bar(["walls F"], [F], bottom=[top_in], color=COLORS["muted"], label=f"F = {F:.0f} N")   # wall force: a step DOWN (waterfall)
c.bar(["out"], [mdot * s2["u"]], color=COLORS["blue"])                              # momentum flux out
c.bar(["out"], [s2["p"] * A2], bottom=[mdot * s2["u"]], color=COLORS["orange"])     # pressure force on face 2
c.set(title="momentum [N]"); c.legend(fontsize=8, loc="upper right")
fig.suptitle("Three budgets between station 1 and the throat of our nozzle: each balances", fontsize=11)
savefig(fig, "ch15", "c03_budgets"); plt.show()              # save and display
""", see=r"""Bars for station 1 ("in") and the throat ("out"): mass flow; enthalpy (rose) plus kinetic energy (blue) with the
conserved total as an amber line; and, on the right, momentum flux (blue) plus pressure force (orange) at each station with
the wall force (grey) as a step between them.""",
    read=r"""Mass and total enthalpy have equal heights in and out: what enters leaves, and enthalpy is traded for kinetic
energy at a fixed total. On the right, read the bars as a waterfall: start at the top of "in", step *down* by the grey wall
force (it is negative, −426 N) and you land exactly on the top of "out".""",
    change=r"""…you added heat Q? Only the middle pair would change height (the outlet bar rises by Q). Add friction: only
the third — and the middle pair still balances.""")
whatif(r"""…you slowed the gas to rest without heat or friction? The energy line says its temperature would rise to a fixed
value $T_0$ — the ruler of the next block.""")

# =====================================================================================================================
# A.4 §15.4 Reference Properties in Compressible Flow — C04, C05
# =====================================================================================================================
nb.section("15.4", "Reference Properties in Compressible Flow", intro="**What is this section about?** Two imagined states "
           "that every point of a flow carries with it: the state it would reach if brought gently to rest (stagnation), and "
           "the state it would reach if brought gently to Mach 1 (sonic). All tables and nozzle results are ratios to them.")
core("C04", r"The stagnation ruler: $\dfrac{T_0}{T}=1+\dfrac{\gamma-1}{2}M^2$ (15.28)",
     "How hot and how compressed would this gas be if you stopped it gently — and which of those two numbers does a real flow keep?")
problem(r"""
Hold a thermometer out of an aircraft window at Mach 0.85. It does not read the air's temperature: the air that hits it is
brought to rest, its kinetic energy becomes heat, and the reading is about 32 K too high. That higher temperature — and the
higher pressure a pitot tube reads — is a property of *each point* of the flow, defined by an imagined experiment. It is the
ruler of the chapter: $T_0$ is what an adiabatic flow keeps; $p_0$ is what it loses whenever entropy is made.""")
idea(r"""
moving gas (T, p, u)  ──[ stop it: no heat, no friction ]──►  at rest (T₀, p₀, 0)
""", r"""
| $T_0$ needs only "no heat" | $p_0$ also needs "no friction" |
|---|---|
| kept across friction and shocks | falls whenever entropy rises |""")
note("N29 [B]", r"""**Stagnation (total) state:** the state the local gas would reach if brought to rest isentropically;
stagnation enthalpy $h_0\equiv h+\tfrac12u^2$ (`GD.total_enthalpy`). > ⚠️ Common confusion (trap T1): this is the third
meaning of the subscript 0; most tables write $T_t$, $p_t$.""")
recap("R14", "Stagnation temperature", r"""«15.27»: the energy law «15.22» with Q = 0, between the point and its imagined
state of rest.""", where="Ch. 4 §4.9 (N94)", core_id="C04")
code(r"""
print(BER.stagnation_temperature(220.0, 250.0, cp=air.cp))   # (15.27): T0 = T + u²/(2 cp) at 220 K, 250 m/s -> 251.1 K (M = 0.841: +31 K)
""")
D("D08", "The stagnation ratios", ref="15.28", notes=r"""**N30 [B]** «15.29» (steps 4–5) · **N31 [B]** «15.30» (step 6).
The exponents 3.5 and 2.5 are the gamma exponents of primer P355 (§15.1).""")
worked("Mach 2 by hand", r"""
γ = 1.4, M = 2.
1. «15.28»: $T_0/T$ = 1 + 0.2 × 4 = **1.8**.
2. «15.29»: $p_0/p$ = 1.8^3.5 = 1.8³ × √1.8 = 5.832 × 1.342 = **7.82**.
3. «15.30»: $\rho_0/\rho$ = 1.8^2.5 = 1.8² × √1.8 = 3.24 × 1.342 = **4.35**.
4. Check with the gas law: 7.82/4.35 = 1.80 ✓. A Mach-2 stream at 288 K would heat to 519 K if stopped.""")
code(r"""
print(GD.stagnation_temperature_ratio(2.0), GD.stagnation_pressure_ratio(2.0), GD.stagnation_density_ratio(2.0))   # 1.8, 7.8244, 4.3469
print(GD.isentropic_state(2.0))                              # all the ratios at M = 2; p/p0 = 0.1278 (NACA Report 1135 row M = 2)
print(GD.mach_from_pressure_ratio(1 / 7.8244))               # the inverse of (15.29): back to M = 2.000
""", explain=r"""
1. The three functions are «15.28»–«15.30»: 1.8, 7.8244, 4.3469.
2. `GD.isentropic_state(M)` returns all ratios at once, both ways round (`p_over_p0` = 0.1278 and `p0_over_p` = 7.82); its
   values are checked in the tests against the public NACA Report 1135 table.
3. `GD.mach_from_pressure_ratio` inverts «15.29»: what a pitot-static pair measures, turned into a Mach number.""")
scratch(r"""
g = 1.4                                                      # gamma of air
M = np.linspace(0.0, 5.0, 51)                                # Mach numbers 0 … 5
t = 1 + 0.5 * (g - 1) * M**2                                 # (15.28): T0/T
mine = dict(T_over_T0=1 / t, p_over_p0=t**(-g / (g - 1)), rho_over_rho0=t**(-1 / (g - 1)))   # (15.28)-(15.30) inverted
lib = GD.isentropic_state(M, g)                              # the tested function (arrays in, arrays out)
assert all(np.allclose(mine[k], lib[k]) for k in mine)       # same numbers for every Mach number
print("from scratch = GD.isentropic_state at 51 Mach numbers")
""", explain=r"""**From scratch.** One line per ratio, written from «15.28»–«15.30»; they agree with `GD.isentropic_state` everywhere.""")
P("P359", "dynamic pressure in a gas", r"""$\tfrac12\rho U^2=\tfrac12\gamma pM^2$, because $c^2=\gamma p/\rho$ gives
$\rho U^2=\rho c^2M^2=\gamma pM^2$. It lets us trade a speed for a Mach number in any momentum term.""", code=r"""
g, p, M = 1.4, 101325.0, 0.3                  # air at sea level, Mach 0.3
rho, c = 1.225, 340.3                         # density, sound speed
print(0.5*rho*(M*c)**2, 0.5*g*p*M**2)         # 6384 Pa both ways
""")
remind("C04")
note("N32 [B]", r"""For small M, «o_pitot» and $\rho_0/\rho=1+\tfrac12M^2+\dots$ (ours; result given, not derived): the
incompressible pitot formula of Ch. 4 is the first term.""")
code(r"""
g_, M_ = sp.symbols("gamma M", positive=True)                # symbols
print(sp.series((1 + (sp.Rational(7, 5) - 1) / 2 * M_**2)**sp.Rational(7, 2), M_, 0, 6))   # p0/p for gamma = 7/5: 1 + 7M²/10 + 7M⁴/40 + …
print(GD.pitot_compressibility_factor(0.3))                  # (p0 - p)/(½ rho u²) at M = 0.3 -> 1.0227: the incompressible formula reads 2.3 % low
print(GD.mach_for_density_change(0.05))                      # the Mach number at which rho0/rho - 1 = 5 % -> 0.314
""", explain=r"""
1. sympy's `series` expands «15.29» for γ = 7/5: $1+\tfrac{7}{10}M^2+\tfrac{7}{40}M^4+\dots$; the first correction is
   $\tfrac12\gamma M^2\cdot p=\tfrac12\rho u^2$ (primer P359), the second is ¼M² of that — the pitot factor.
2. At M = 0.3 the compressible pitot pressure is 2.3 % above the incompressible one.
3. Density differs from its stagnation value by 5 % at M = 0.314 — the origin of the rule of thumb "incompressible below
   M ≈ 0.3".""")
plotly(r"""
M = np.linspace(0.0, 2.0, 150)                               # Mach numbers
gammas = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 5 / 3]                         # gamma from 1.1 to 5/3
fig = slider_figure(lambda g: {"p0/p, (15.29)": (M, GD.stagnation_pressure_ratio(M, g)),     # the exact ratio
                               "incompressible 1 + ½γM²": (M, 1 + 0.5 * g * M**2),           # Bernoulli's ghost
                               "M = 0.3": ([0.3, 0.3], [1.0, 12.0])},                         # where they part
                    "γ", gammas, xlabel="Mach number M [–]", ylabel="p₀/p (stagnation / static) [–]",
                    title="Below M ≈ 0.3 the compressible and incompressible stagnation pressures agree", yrange=[1, 12],
                    active=gammas.index(1.4))                # open on air, gamma = 1.4
fig.show()                                                   # every slider position computed now
""", explain=r"""**Figure F2.** The slider is γ, from 1.1 to 5/3 in 12 positions; the two curves are «15.29» and its incompressible
first term.""")
notes3(see=r"""The stagnation-to-static pressure ratio against Mach number (purple) and the incompressible estimate
$1+\tfrac12\gamma M^2$ (teal), with a dashed line at M = 0.3.""",
       read=r"""Left of the dashed line the two curves are one: Bernoulli's equation is fine there. At M = 1 the exact ratio
is 1.89 against 1.70; at M = 2 it is 7.8 against 3.8. Drag γ: a stiffer gas (larger γ) gains pressure faster.""",
       change=r"""…the gas were monatomic (γ = 5/3)? The pressure exponent γ/(γ − 1) drops from 3.5 to 2.5, yet the curve
rises faster at small M because the first term is ½γM².""")
recap_md("Entropy is a state function", "Ch. 1 C33–C35", r"""S depends only on the state, so «15.1s» holds between *any* two
states, however the gas got from one to the other; the logarithm rules ln(ab) = ln a + ln b and n ln a = ln aⁿ do the rest.""")
D("D09", "$T_0$ is kept, $p_0$ is lost", ref="u15.4", notes=r"""**N33 [B]** «u15.4»""")
confusion(r"""(trap T17) Adiabatic is not isentropic — $T_0$ is kept, $p_0$ is not.""")
code(r"""
dS = GD.entropy_change_from_stagnation(0.9, gas=air)          # (unnumbered, §15.4): a 10 % loss of p0 at fixed T0 -> 30.2 J/(kg K)
print(dS, GD.stagnation_pressure_loss(dS, gas=air))          # ... and back: p02/p01 = exp(-dS/R) = 0.900
""")
figure(r"""
Ms = [0.3, 0.85, 2.0, 3.0]                                   # four Mach numbers
st = [GD.isentropic_state(M) for M in Ms]                    # their ratios to the stagnation state
xpos = np.arange(len(Ms))                                    # bar groups
fig, ax = plt.subplots(figsize=(7.0, 3.6))                   # one panel
w = 0.26                                                     # bar width
ax.bar(xpos - w, [s["T_over_T0"] for s in st], w, color=COLORS["rose"], label="T/T₀")       # temperature (rose)
ax.bar(xpos, [s["p_over_p0"] for s in st], w, color=COLORS["orange"], label="p/p₀")          # pressure (orange)
ax.bar(xpos + w, [s["rho_over_rho0"] for s in st], w, color=COLORS["teal"], label="ρ/ρ₀")    # density (teal)
ax.axhline(1.0, color=COLORS["amber"], lw=2, label="stagnation value")                         # the ruler
ax.set_xticks(xpos, [f"M = {M}" for M in Ms])               # group labels
ax.set(ylabel="static / stagnation [–]", title="The stagnation ruler: pressure falls fastest (exponent 3.5)")
ax.legend(fontsize=8)                                        # legend
for i, s in enumerate(st):                                   # print the numbers on top
    for dx, k in ((-w, "T_over_T0"), (0, "p_over_p0"), (w, "rho_over_rho0")):
        ax.text(i + dx, s[k] + 0.02, f"{s[k]:.3f}", ha="center", fontsize=7)
savefig(fig, "ch15", "c04_ruler"); plt.show()                # save and display
""", see=r"""For four Mach numbers, the static temperature, pressure and density as fractions of their stagnation values,
with the stagnation value as an amber line at 1.""",
    read=r"""At M = 0.3 everything is within about 6 % of the stagnation value (the pressure, the farthest, is at 0.939). At M = 2 the pressure is 13 % of it, the density
23 % and the temperature 56 %; at M = 3 the pressure is under 3 %. **The pressure always falls fastest**: its exponent
γ/(γ − 1) = 3.5 is the largest.""",
    change=r"""…γ were 5/3? The exponents would be 2.5 and 1.5, so pressure and density would fall more slowly for the same
temperature drop — but the temperature drop itself is larger, $\tfrac13M^2$ instead of $\tfrac15M^2$.""")
note("N34 [C]", r"""A large reservoir at rest holds the stagnation values of everything that flows out of it isentropically;
it is the left end of the nozzle figures of C06, C09 and C10.""")
note("N38 [C]", r"""Let $T\to0$ in «15.28»: the speed tends to «o_umax» — `GD.max_speed(288.15, air.cp)` = 761 m/s for a
reservoir at room temperature; the same limit caps the Prandtl–Meyer turn in C15.""")
code(r"""
print(GD.max_speed(288.15, air.cp))                          # sqrt(2 cp T0): all the enthalpy turned into speed -> 761 m/s
""")
whatif(r"""…instead of stopping the gas you brought it to exactly Mach 1? You get a second reference state — and with it
the area a flow needs. C05.""")

# ---------------------------------------------------------------------------------------------------------------------
# C05 — the sonic state and the area ratio
# ---------------------------------------------------------------------------------------------------------------------
core("C05", r"One area, two Mach numbers: $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31)",
     "How much room does a given flow need at each Mach number — and why is there a least amount?")
problem(r"""
Push a fixed number of kilograms per second down a duct. Where the gas is slow it needs a wide passage. Speed it up and it
needs less — until Mach 1. Beyond that the gas thins out so quickly that it needs *more* room again. So the required area
has a minimum, exactly at Mach 1: that minimum is the throat, and a throat of given size can pass only so much.""")
idea(r"""
A/A* │\                 /
     │ \               /
     │  \             /
   1 ┼───\___________/────
     └──────── M = 1 ─────── M
       one height → two Mach numbers (subsonic, supersonic)
""")
remind("C05")
note("N35 [B]", r"""**Sonic (critical) state** — not in the book as formulas, ours: put M = 1 in «15.28», «15.29», «15.30»:
«o_sonic». For γ = 1.4: 0.833, 0.528, 0.634. > ⚠️ Common confusion (trap T2): the star now means "if brought to M = 1
isentropically".""")
code(r"""
print(GD.sonic_ratios())                                     # T*/T0, p*/p0, rho*/rho0, c*/c0 for gamma = 1.4 -> 0.833, 0.528, 0.634, 0.913
print(GD.critical_pressure_ratio())                          # p*/p0 = 0.5283: the number behind choking (C09)
""")
D("D10", "The area ratio", ref="15.31")
worked("Mach 2 and Mach 0.5 by hand", r"""
γ = 1.4, exponent ½(γ+1)/(γ−1) = 3.
1. M = 2: bracket = (2/2.4) × 1.8 = 1.5; 1.5³ = 3.375; ÷ M: **1.6875**.
2. M = 0.5: bracket = (1/1.2) × 1.05 = 0.875; 0.875³ = 0.670; ÷ 0.5: **1.340**.
3. M = 1: bracket = 1: **1**. Both sides of Mach 1 need more area than the throat.""")
P("P360", "a function with two inverses", r"""A U-shaped function has two x for each y above its minimum. A root-finder
needs to be told which: give it a bracket on one side of the minimum. Here the physics chooses the side.""", code=r"""
from scipy.optimize import brentq                 # bracketed root finder (reminder: P108)
A = lambda M: (1/M)*((2/2.4)*(1 + 0.2*M*M))**3    # A/A* for gamma = 1.4
print(brentq(lambda M: A(M) - 2.0, 1e-8, 1.0))    # 0.306 (subsonic side)
print(brentq(lambda M: A(M) - 2.0, 1.0, 100.0))   # 2.197 (supersonic side)
""")
note("N36 [B]", r"""(ours) The inverse of «15.31»: `GD.mach_from_area_ratio(2.0, branch="subsonic")` = 0.306,
`branch="supersonic"` = 2.197. > ⚠️ Common confusion (trap T8, first half): one area ratio, two Mach numbers — which one
occurs is decided downstream, by the back pressure (C09–C10).""")
code(r"""
print(GD.area_ratio(np.array([0.5, 1.0, 2.0, 3.0])))         # (15.31): 1.3398, 1, 1.6875, 4.2346
for AR in (1.5, 2.0):                                        # two area ratios
    print(AR, GD.mach_from_area_ratio(AR, branch="subsonic"), GD.mach_from_area_ratio(AR, branch="supersonic"))   # (0.4303, 1.8541), (0.3059, 2.1972)
""", explain=r"""
1. `GD.area_ratio` evaluates «15.31»: the worked numbers 1.340 and 1.6875, the minimum 1 at M = 1, and 4.23 at M = 3.
2. `GD.mach_from_area_ratio` inverts it on the side you name with `branch=` — we always name it.""")
scratch(r"""
g, M = 1.4, np.array([0.3, 0.5, 2.0, 3.0])                    # gamma; four Mach numbers
Theta = 1 + 0.5 * (g - 1) * M**2                              # D10 step 4: T0/T
rho_star_over_rho0 = (2 / (g + 1))**(1 / (g - 1))             # D10 step 5
rho0_over_rho = Theta**(1 / (g - 1))                          # (15.30)
T_star_over_T = (2 / (g + 1)) * Theta                         # D10 step 3: T*/T0 × T0/T
A_mine = rho_star_over_rho0 * rho0_over_rho * np.sqrt(T_star_over_T) / M   # D10 step 2: a product of four ratios
assert np.allclose(A_mine, GD.area_ratio(M))                  # same as the library

def bisect(f, a, b, n=80):                                    # a ten-line bisection: halve [a, b] n times
    for _ in range(n):                                        # each pass halves the bracket
        m = 0.5 * (a + b)                                     # midpoint
        if f(a) * f(m) <= 0: b = m                            # root in the left half
        else: a = m                                           # root in the right half
    return 0.5 * (a + b)                                      # the midpoint of the last bracket
F = lambda M: GD.area_ratio(M) - 2.0                          # A/A* - 2
M_sub, M_sup = bisect(F, 1e-6, 1.0), bisect(F, 1.0, 50.0)     # one root on each side of the minimum
assert np.allclose([M_sub, M_sup], [GD.mach_from_area_ratio(2.0, branch="subsonic"),
                                    GD.mach_from_area_ratio(2.0, branch="supersonic")], atol=1e-9)
print(A_mine, M_sub, M_sup)                                   # 2.035 1.340 1.6875 4.235 | 0.3059 2.1972
""", explain=r"""**From scratch.** The area ratio is assembled exactly as D10 does it — four ratios multiplied — and inverted by
plain bisection on each side of Mach 1. Both agree with the library.""")
P("P361", "a generated table as a pandas DataFrame", r"""The classical tool is a printed table of ratios against Mach number.
We generate ours: a DataFrame is a table with named columns; `.loc` picks rows; reading it "backwards" means finding the
row whose column is closest to a value.""", code=r"""
tab = GD.isentropic_table(np.arange(0.5, 3.01, 0.5))     # our own table, gamma = 1.4
print(tab.round(4))                                      # columns 'M', 'p/p0', 'rho/rho0', 'T/T0', 'A/A*'
print(tab.loc[(tab['A/A*'] - 1.69).abs().idxmin(), 'M']) # read backwards: nearest row to A/A* = 1.69
""")
note("N37 [B]", r"""The book's table is not reproduced; ours is generated for γ = 1.4 and 5/3 and checked against three rows
of the public-domain NACA Report 1135 (M = 0.5, 2, 3). (Slip #9 concerns two misprinted entries of the book's table; a
generated table cannot have them.) Note the last line of the demo: reading backwards returns only the nearest row — and only
one of the two branches.""")
D("D11", "The sonic ratios and the choked mass flow", ref="o_mmax", notes=r"""**N35 [B]** the sonic ratios (steps 1–2) ·
**N39 [B]** the choked mass flow «o_mmax» (steps 3–7)""")
code(r"""
print(GD.mass_flux_function(np.array([0.5, 1.0, 2.0])))      # mdot sqrt(R T0)/(p0 A) at three Mach numbers -> 0.5111, 0.6847, 0.4058
print(GD.choked_mass_flow(noz["p0"], noz["T0"], noz["A_t"], air))   # our 10 cm² throat at 5 bar, 400 K -> 1.010 kg/s
""", explain=r"""The flux function «o_flux» is largest at M = 1, where it equals $\sqrt\gamma\,(2/(\gamma+1))^3$ = 0.6847;
with our nozzle's supply that is 1.010 kg/s — the same number C03 found at the throat.""")
plotly(r"""
M = np.geomspace(0.05, 5.0, 200)                             # Mach numbers, log-spaced
gammas = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 5 / 3]                         # gamma slider
def area(g):                                                 # the curves for one gamma
    ms, mp = GD.mach_from_area_ratio(2.0, g, branch="subsonic"), GD.mach_from_area_ratio(2.0, g, branch="supersonic")   # the two roots of A/A* = 2
    return {"A/A*, (15.31)": (M, GD.area_ratio(M, g)), "A/A* = 2: the two roots": ([ms, mp], [2.0, 2.0])}   # curve + two dots
fig_a = slider_figure(area, "γ", gammas, xlabel="Mach number M [–] (log axis)", ylabel="A/A* (area / sonic area) [–]",   # one figure, slider over gamma
                      title="Each area above the throat belongs to two Mach numbers", yrange=[0, 6],   # the message; fixed y-range
                      modes={"A/A* = 2: the two roots": "markers"}, active=gammas.index(1.4))   # roots as dots; open on air
fig_a.update_xaxes(type="log", range=[np.log10(0.05), np.log10(5.0)])   # log axis: subsonic and supersonic get equal room
fig_a.show()   # display
fig_b = slider_figure(lambda g: {"mass flux function (ours)": (M, GD.mass_flux_function(M, g))}, "γ", gammas,   # the flux function for each gamma
                      xlabel="Mach number M [–] (log axis)", ylabel="ṁ√(RT₀)/(p₀A) [–]",   # axis labels
                      title="The mass flux through an area is largest when the flow there is sonic",   # the message
                      active=gammas.index(1.4))               # open on air, gamma = 1.4
fig_b.update_xaxes(type="log", range=[np.log10(0.05), np.log10(5.0)])   # same axis
fig_b.show()   # display
""", explain=r"""**Figures F3.** Two slider figures over γ, from 1.1 to 5/3: «15.31» with the two roots of A/A* = 2 marked, and the
mass flux function «o_flux». The horizontal axes are logarithmic (each tick a factor 10), so slow and fast flow get equal
room.""")
notes3(see=r"""Top: the U-shaped area ratio with its minimum 1 at M = 1 and two dots where it equals 2. Bottom: the mass flux
per unit area with its maximum at M = 1.""",
       read=r"""The two figures are one curve turned upside down: flux × (A/A*) is the same at every Mach number, because a
fixed mass flow needs more area exactly where the flux per unit area is smaller. Drag γ: the peak flux and the spread of the
two roots both change.""",
       change=r"""…the throat area were halved at the same supply? The flux function would not change — but the most the throat
can pass, «o_mmax», would halve.""")
whatif(r"""…you asked not "how much area" but "which way does the speed change when the area changes"? Differentiate — and a
factor $1-M^2$ appears. C06.""")

# =====================================================================================================================
# A.5 §15.5 Area-Velocity Relationship in One-Dimensional Isentropic Flow — C06
# =====================================================================================================================
nb.section("15.5", "Area-Velocity Relationship in One-Dimensional Isentropic Flow", intro="**What is this section about?** "
           "Which way the speed changes when a duct narrows or widens. Below Mach 1 a nozzle converges; above Mach 1 it must "
           "diverge. The first of the chapter's three surprises.")
core("C06", r"A supersonic nozzle widens: $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34)",
     "Why does squeezing a duct speed a slow flow up but slow a supersonic one down?")
problem(r"""
Put your thumb over a garden hose: the water speeds up. Everyone's intuition says "squeeze = faster". Now look at the bell of
a rocket engine: the gas is supersonic and the bell *opens out* — and the gas keeps accelerating. The mass flow ρuA is fixed.
In water ρ cannot change, so u must rise when A falls. In a gas ρ falls as u rises — and above Mach 1 it falls faster than u
rises, so the product ρu *drops* and the gas needs more room.""")
idea("", r"""
| M < 1 | M > 1 |
|---|---|
| speed +1 % costs density −M² % (less than 1 %) | speed +1 % costs density −M² % (more than 1 %) |
| ρu rises → the gas needs **less** area | ρu falls → the gas needs **more** area |""")
note("N40 [B]", r"""A **nozzle** is whatever speeds the flow up and drops its pressure; a **diffuser** does the opposite. They
are defined by what they do, not by their shape — and a diffuser has nothing to do with molecular diffusion.""")
recap("R15", "The constant-density limit", r"""With ρ fixed, «15.19» is $\dfrac{dA}{A}+\dfrac{du}{u}=0$: the garden hose.

> ⚠️ **slip #5 — the book prints** $dA/dx+du/dx=0$; **the correct form is** $\dfrac1A\dfrac{dA}{dx}+\dfrac1u\dfrac{du}{dx}=0$
(a metre cannot be added to a per-second).""", where="Ch. 4 continuity", core_id="C06")
code(r"""
A, u, x = 1e-3, 30.0, 1.0                                    # an area [m²], a speed [m/s], a length [m]
print("printed terms: dA/dx in m²/m = m;  du/dx in (m/s)/m = 1/s  -> they cannot be added")   # units of the printed form
print("logarithmic terms: (1/A) dA/dx and (1/u) du/dx are both in 1/m -> they can")            # units of the corrected form
print(GD.influence_coefficients(0.0, dA_over_A=0.01)["du_over_u"])   # at M = 0: du/u = -dA/A = -0.01
""")
remind("C06")
P("P362", "differentials of algebraic laws as a small linear system", r"""Take the differential of each law and divide by
the law: every law becomes a *linear* equation between the relative changes $d\rho/\rho$, $du/u$, $dp/p$, … Then solve them
like simultaneous equations.""", code=r"""
import sympy as sp
dr, du, dA, M = sp.symbols('dr du dA M')              # d rho/rho, du/u, dA/A, Mach number
sol = sp.solve([dr + du + dA, dr + M**2*du], [dr, du])   # mass (15.19) and Euler (15.33)
print(sol[du])                                        # dA/(M**2 - 1), i.e. (15.34)
""")
D("D12", "The area–velocity relation", ref="15.34", notes=r"""**N41 [B]** Euler's equation for isentropic change,
«15.32» (steps 1–3, the minus sign restored) · **N42 [B]** «15.33»: density against speed (step 4)""")
slip(4, r"«15.32p»", r"«15.32»", why=r"""Why: $u\,du=-dp/\rho$ and $dp=+c^2d\rho$ leave no choice, and the book's own next line
«15.33» needs the minus. The cell below shows what the printed sign would do — and a test of the library must fail on it.""")
code(r"""
good = GD.influence_coefficients(2.0, dA_over_A=0.01)                        # corrected (15.32) at M = 2, 1 % more area
bad = GD.influence_coefficients(2.0, dA_over_A=0.01, printed_15_32=True)     # the sign as printed
print(f"corrected: du/u = {good['du_over_u']:+.5f}, drho/rho = {good['drho_over_rho']:+.5f}")   # +0.00333, -0.01333
print(f"printed:   du/u = {bad['du_over_u']:+.5f}, drho/rho = {bad['drho_over_rho']:+.5f}")    # -0.00200, -0.00800
for M in (2.0, 0.3):                                         # N42 read aloud: density cost of a 1 % speed gain
    print(f"M = {M}: a 1 % gain in speed costs {M**2:.2f} % of the density")   # (15.33): drho/rho = -M² du/u
""", explain=r"""
1. With the corrected sign, at M = 2 a 1 % wider duct speeds the gas up by 0.33 % and thins it by 1.33 % — a nozzle.
2. With the printed sign, density and speed would change in the *same* sense ($d\rho/\rho=+M^2du/u$), the factor $1-M^2$
   would become $1+M^2$, and a widening duct would slow every flow: no supersonic nozzle could exist.
3. **N42** in numbers: at M = 2 a 1 % gain in speed costs 4 % of the density; at M = 0.3 only 0.09 %.""")
note("N43 [B]", r"""The four cases of «15.34»:

| | converging (dA < 0) | diverging (dA > 0) |
|---|---|---|
| **M < 1** | nozzle (faster) | diffuser (slower) |
| **M > 1** | diffuser (slower) | nozzle (faster) |""")
figure(r"""
fig = draw_four_ducts()                                      # our sketch of the four cases
plt.show()                                                   # display it
for M in (0.5, 2.0):                                         # subsonic, supersonic
    print(f"M = {M}: converging -> {GD.duct_role(M, -1)}, diverging -> {GD.duct_role(M, +1)}")   # (15.34) decides
print(GD.duct_role(0.5, 0), GD.duct_role(1.0, -1))           # no area change: "none"; at M = 1: "sonic"
""", see=r"""Four ducts: converging and diverging, each with a subsonic and a supersonic stream; the printed lines give the
role `GD.duct_role` assigns to each.""",
    read=r"""Read the table diagonally: a converging subsonic duct and a diverging supersonic duct are both nozzles; the
other two are diffusers. The shape alone never tells you which.""",
    change=r"""…the flow were exactly sonic? Then «15.34» divides by zero: a sonic flow cannot sit in a duct whose area
changes, only at a throat (N44).""")
note("N44 [B]", r"""At a throat dA = 0, so «15.34» gives either du = 0 or M = 1: **sonic speed can occur only at a throat, but
a throat need not be sonic.** (The 0/0 at M = 1 has a finite limit — l'Hôpital's rule.)""")
figure(r"""
x = np.linspace(0.0, 1.0, 201)                               # 201 stations: the throat x = 0.4 is a grid point
fig, (a, b) = plt.subplots(2, 1, figsize=(6.8, 5.0), sharex=True)   # area above, Mach number below
a.plot(x, GD.laval_area(x), color=COLORS["ink"])             # the duct: two parabolas, throat area 1 at x = 0.4
a.set(ylabel="A / A_throat [–]", title="One duct, three isentropic flows")
for branch, A_star, col, lab in (("subsonic", 0.9, COLORS["blue"], "subsonic throughout (a venturi)"),
                                 ("transonic", None, COLORS["accent"], "transonic: sonic at the throat"),
                                 ("supersonic", 0.9, COLORS["rose"], "supersonic throughout")):
    sol = GD.nozzle_isentropic(x, GD.laval_area, branch=branch, A_star=A_star)   # (15.31) inverted station by station
    b.plot(x, sol["M"], color=col, label=lab)                # Mach number along the duct
b.axhline(1.0, color=COLORS["amber"], ls="--", lw=1)         # M = 1
b.axvline(0.4, color=COLORS["muted"], ls=":", lw=1)          # the throat
b.set(xlabel="x [m]", ylabel="Mach number M [–]"); b.legend(fontsize=8)
savefig(fig, "ch15", "c06_three_solutions"); plt.show()      # save and display
""", see=r"""A convergent–divergent duct (top) and three steady isentropic flows through it (bottom): subsonic everywhere
(blue, with a sonic area 10 % below the throat area), the transonic one that is sonic at the throat (purple), and supersonic
everywhere (rose).""",
    read=r"""The blue flow speeds up to the throat and slows after it (a venturi): the throat is a speed *maximum* with
du = 0. The rose flow does the opposite. Only the purple flow crosses M = 1 — and it does so exactly at the throat.""",
    change=r"""…the sonic area of the blue flow rose to the throat area? Its peak would touch M = 1 at the throat; from there
it can continue on the purple curve.""")
note("N45 [B]", r"""(ours; result given, not derived) the Mach-number twin of «15.34»: «o_dMM».""")
code(r"""
M_, g_ = sp.symbols("M gamma", positive=True)                # symbols
lnA = -sp.log(M_) + (g_ + 1) / (2 * (g_ - 1)) * (sp.log(2 / (g_ + 1)) + sp.log(1 + (g_ - 1) / 2 * M_**2))   # ln of (15.31), written as a sum of logs
resid = sp.diff(lnA, M_) * M_ + (1 - M_**2) / (1 + (g_ - 1) / 2 * M_**2)   # d lnA / d lnM + (1 - M²)/Θ: should vanish
print(sp.simplify(sp.together(resid)))                         # 0: so dM/M = -Θ/(1 - M²) dA/A, the Mach-number twin of (15.34)
""")
worked("1 % more area", r"""
«15.34» with dA/A = +0.01.
1. M = 0.5: du/u = −0.01/(1 − 0.25) = **−1.33 %** (slower; a diffuser).
2. M = 2: du/u = −0.01/(1 − 4) = **+0.33 %** (faster; still a nozzle).
3. Density by «15.33»: at M = 2, dρ/ρ = −4 × 0.33 % = −1.33 %.
4. Check mass «15.19»: −1.33 + 0.33 + 1.00 = 0 ✓.""")
code(r"""
print(GD.area_velocity(np.array([0.3, 0.5, 2.0]), 1.0))       # du/u per unit dA/A: -1/(1 - M²) -> -1.099, -1.333, +0.333
c = GD.influence_coefficients(0.5, dA_over_A=0.01)           # all the relative changes at M = 0.5 for 1 % more area
print({k: round(c[k], 5) for k in ("du_over_u", "drho_over_rho", "dp_over_p", "dT_over_T", "dM_over_M")})   # -1.33 %, +0.33 %, …
""", explain=r"""
1. `GD.area_velocity(M, 1.0)` is the coefficient $-1/(1-M^2)$ of «15.34»: negative (diffuser-like for a widening duct) below
   Mach 1, positive above.
2. `GD.influence_coefficients` solves the whole small linear system of P362 at once: at M = 0.5 a 1 % wider duct gives
   du/u = −1.33 %, dρ/ρ = +0.33 %, and dM/M = −1.40 % (the Mach number falls faster than the speed because the gas warms).""")
scratch(r"""
x = np.linspace(0.0, 1.0, 201)                               # stations with the throat on the grid
sol = GD.nozzle_isentropic(x, GD.laval_area, branch="transonic")   # the transonic solution
u = sol["M"] * np.sqrt(sol["T_over_T0"])                     # u ∝ M sqrt(T): u = M sqrt(gamma R T0 · T/T0)
slope = np.gradient(np.log(u), x) / np.gradient(np.log(sol["A"]), x)   # (du/u)/(dA/A) by finite differences (np.gradient)
mask = np.abs(sol["M"] - 1.0) > 0.2                          # stay away from the throat, where both differences vanish
assert np.allclose(slope[mask], GD.area_velocity(sol["M"][mask], 1.0), rtol=2e-2)   # the measured slope obeys (15.34)
print(f"{mask.sum()} stations checked; worst relative gap {np.max(np.abs(slope[mask]/GD.area_velocity(sol['M'][mask], 1.0) - 1)):.1e}")
""", explain=r"""**From scratch.** Along the computed nozzle flow we *measure* the ratio of the relative speed change to the
relative area change with `np.gradient` (centred differences of an array) and compare it with $-1/(1-M^2)$: «15.34» holds
station by station.""")
plotly(r"""
x = np.linspace(0.0, 1.0, 51)                                # stations (the throat x = 0.4 is a grid point)
ratios = np.linspace(1.2, 4.0, 6 if FAST else 10)            # exit-to-throat area ratios on the slider
def three(AR):                                               # the three isentropic solutions for one nozzle
    A_fn = lambda xx: GD.laval_area(xx, A_exit=AR)           # the duct with this exit area (throat area 1)
    out = {}   # trace name -> (x, M(x))
    for branch, A_star, name in (("subsonic", 0.9, "subsonic throughout"), ("transonic", None, "transonic"),   # the three solutions
                                 ("supersonic", 0.9, "supersonic throughout")):   # (sonic area 10 % below the throat for the two one-sided ones)
        out[name] = (x, GD.nozzle_isentropic(x, A_fn, branch=branch, A_star=A_star)["M"])   # (15.31) inverted at every station
    out["M = 1"] = (x, np.ones_like(x))                      # the sonic line
    return out   # all traces for this area ratio
fig = slider_figure(three, "A_e/A_t", ratios, xlabel="x [m] (throat at 0.4 m)", ylabel="Mach number M [–]",   # slider over the exit area
                    title="A bigger exit area makes the supersonic exit faster", yrange=[0, 3.5])   # the message
fig.show()                                                   # every position computed now
""", explain=r"""**Figure F4.** The slider is the exit-to-throat area ratio $A_e/A_t$, from 1.2 to 4.""")
notes3(see=r"""M(x) for the three isentropic solutions of one duct, with the sonic line; the slider changes the exit area.""",
       read=r"""The converging half never changes — it does not know about the exit. The transonic exit Mach number climbs
from 1.5 at $A_e/A_t$ = 1.2 to 2.9 at 4; the subsonic one falls toward zero.""",
       change=r"""…the gas had γ = 5/3? For the same area ratio the supersonic exit would be faster (2.40 instead of 2.20 for
$A_e/A_t$ = 2): the exponent of «15.31» drops from 3 to 2, so the area ratio grows more slowly with M.""")
explainer("area_mach_nozzle", "Why does a supersonic nozzle widen?",
          r"""you move a cursor along a duct you shaped and three bars — density, speed, area — always sum to zero while the
density bar overtakes the speed bar exactly at Mach 1; a static figure shows the result, the bars show the reason.""",
          ["Play the station through the throat and watch the density bar overtake the speed bar.",
           "Switch the solution to 'subsonic throughout': the same throat, no Mach 1.",
           "Find the ghost dot: the other Mach number with the same area."])
note("N46 [B]", r"""**A rocket nozzle with our own inputs** (`rocket`; air's γ and R, a simplification). The chain: $p_e/p_0$
→ $M_e$ (inverse of «15.29») → $A_e/A^*$ («15.31») → $T_e$ → $u_e$ → ṁ = thrust/$u_e$ → $A_e$ → $A^*$. > ⚠️ Common confusion
(trap T9): the chamber *is* the stagnation state, and the nozzle is designed so that the exit pressure equals the ambient
one.""")
code(r"""
r = ch15.rocket_nozzle(rocket["thrust"], rocket["p0"], rocket["T0"], rocket["p_exit"], air)   # the chain above
print({k: round(float(r[k]), 4) for k in ("M_e", "T_e", "u_e", "mdot", "A_e", "Ae_over_At", "A_throat")})
print(GD.choked_mass_flow(rocket["p0"], rocket["T0"], r["A_throat"], air))   # check: that throat passes the same 25.86 kg/s
""", explain=r"""
1. $p_e/p_0$ = 101 325/3 × 10⁶ = 0.0338 → $M_e$ = 2.857 → $A_e/A^*$ = 3.70; $T_e$ = 1140 K; $u_e$ = 1933 m/s.
2. ṁ = 50 000/1933 = 25.9 kg/s; the exit area is 0.043 m² and the throat 0.0117 m² (a 12 cm diameter).
3. `GD.choked_mass_flow` with that throat returns the same mass flow: the chain closes.""")
note("N47 [B]", r"""Thrust of a pressure-matched nozzle: $\rho_eu_e^2A_e=\dot mu_e=\gamma p_eM_e^2A_e$. Off design a pressure
term appears: thrust = $\dot mu_e+(p_e-p_a)A_e$.""")
code(r"""
T10, p10, rho10 = STAT.standard_atmosphere(10000.0)          # the standard atmosphere at 10 km: 223 K, 26.4 kPa (lapse rate -6.5 K/km in Kundu's sign, +6.5 K/km meteorological)
print(ch15.nozzle_thrust(r["mdot"], r["u_e"], 101325.0, 101325.0, r["A_e"]))   # on design at sea level: 50.0 kN
print(ch15.nozzle_thrust(r["mdot"], r["u_e"], 101325.0, p10, r["A_e"]))        # the same nozzle at 10 km: 53.2 kN
""")
whatif(r"""…the supersonic stream must come back to a high pressure? It cannot do so smoothly in general. It jumps — C07.""")

# =====================================================================================================================
# A.6 §15.6 Normal Shock Waves — C07, C08
# =====================================================================================================================
nb.section("15.6", "Normal Shock Waves", intro="**What is this section about?** A supersonic stream can change to a subsonic "
           "one across a layer a few molecular free paths thick. Across it the three duct laws still hold — and they fix "
           "everything from one number, the upstream Mach number.")
core("C07", r"The jump across a normal shock: $\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left(M_1^2-1\right)$ (15.39)",
     "If a flow changes abruptly with no heat, no wall force and no change of area, what can the state behind it be?")
problem(r"""
Cars on a motorway meet a jam: within a few car lengths they go from fast and far apart to slow and close together. No car is
created or destroyed, and the "front" of the jam is thin. A shock is the same thing in a gas: fast, thin gas arrives; slow,
dense, hot gas leaves. We do not need to know what happens inside the thin layer — only that mass, momentum and energy that go
in must come out.""")
idea(r"""
1: u₁ fast, ρ₁, p₁, T₁  ──►│ shock │──►  2: u₂ slow, ρ₂ > ρ₁, p₂ > p₁, T₂ > T₁
""", r"""Three laws, three unknowns (ρ₂, u₂, p₂) → two solutions: nothing happens, or the shock.""")
note("N48 [B]", r"""A **shock** is a compression of finite strength, a few mean free paths thick (the **mean free path** is
the average distance a molecule travels between collisions — about 7 × 10⁻⁸ m in sea-level air; Ch. 1 C06), treated as a jump. Entropy is
made inside it, so the isentropic ratios of C04–C05 cannot be carried across — only the three laws of C03.""")
note("N49 [B] · N50 [B] · N51 [B]", r"""Mass: «15.35» («15.18» at constant area) · momentum: «15.36» («15.24» with F = 0 and
$A_1=A_2$) · energy: «15.37» («15.22» with Q = 0). > ⚠️ Common confusion (trap T17): the energy line looks like Bernoulli's
equation but is only the energy law; entropy is not constant across the shock.""")
remind("C07")
P("P363", "dividing a known root out of a quadratic", r"""If you already know one solution $r_0$ of a quadratic, the
quadratic is $(r-r_0)$ times something linear — divide and the other root drops out. (Vieta: the product of the roots is the
constant term over the leading coefficient.) Here the known root is "nothing happens".""", code=r"""
import sympy as sp
r = sp.symbols('r')
poly = 2*r**2 - 7*r + 5                 # we know r = 1 is a root: 2 - 7 + 5 = 0
print(sp.cancel(poly/(r - 1)))           # 2*r - 5  ->  the other root is 5/2 (and 1 x 5/2 = 5/2, Vieta)
""")
D("D13", "The pressure jump across a normal shock", ref="15.39", notes=r"""**N52 [B]** the energy law with p and ρ, «15.38»
(step 1) · **N53 [B]** the half-way line
$\dfrac{p_2}{p_1}=1+\dfrac{2\gamma}{\gamma+1}\left[\dfrac{\rho_1u_1^2}{\gamma p_1}-1\right]$ (unnumbered, §15.6) (step 12)""")
D("D14", "The Mach number behind a normal shock", ref="15.40", notes=r"""**N54 [B]** «15.40»""")
confusion(r"""(trap T6) «15.40» is its own inverse — feed it the downstream Mach number and it returns the upstream one
without complaint. `GD.normal_shock` refuses $M_1<1$ unless asked. Other sources number the stations 0 and 1.""")
note("N55 [B]", r"""$M_2$ against $M_1$: a map that is its own mirror image.""")
figure(r"""
M1 = np.linspace(0.4, 5.0, 400)                              # upstream Mach numbers (the formula is defined above 0.378)
M2 = GD.normal_shock_mach(M1)                                # (15.40)
fig, ax = plt.subplots(figsize=(5.6, 4.6))                   # a square-ish panel: the mirror symmetry shows
ax.axvspan(0.4, 1.0, color=COLORS["muted"], alpha=0.15)      # the forbidden half
ax.text(1.08, 4.6, "← shaded: M₁ < 1,\nforbidden (reason: C08)", ha="left", va="top", fontsize=8, color=COLORS["muted"])   # label beside the band, clear of the curve
ax.plot(M1, M2, color=COLORS["blue"], lw=2, label="M₂ from (15.40)")                     # the map
ax.plot([0, 5], [0, 5], color=COLORS["muted"], ls="--", lw=1, label="no change (M₂ = M₁)")   # the diagonal
lim = np.sqrt(0.4 / 2.8)                                     # sqrt((gamma-1)/(2 gamma)) = 0.378
ax.axhline(lim, color=COLORS["amber"], ls=":", lw=1.5, label=f"floor √((γ−1)/2γ) = {lim:.3f}")   # the asymptote
ax.plot([2.0, GD.normal_shock_mach(2.0)], [GD.normal_shock_mach(2.0), 2.0], "o", color=COLORS["orange"])   # the pair 2 <-> 0.577
ax.set(xlim=(0, 5), ylim=(0, 5), xlabel="upstream M₁ [–]", ylabel="downstream M₂ [–]",
       title="The shock map is its own mirror image")       # message
ax.legend(fontsize=8, loc="upper right"); savefig(fig, "ch15", "c07_M2_vs_M1"); plt.show()
""", see=r"""The downstream Mach number against the upstream one (blue), the diagonal "nothing happens" (dashed), the
floor 0.378 that $M_2$ approaches for very strong shocks (amber), and the forbidden half $M_1<1$ shaded. The two orange dots
are the pair 2 ↔ 0.577.""",
    read=r"""Reflect the curve in the diagonal and you get the same curve: apply the map twice and you are back. The two
orange dots are mirror images. Only the part right of $M_1=1$ (supersonic in, subsonic out) is physical.""",
    change=r"""…γ were 5/3? The floor would rise to $\sqrt{1/5}$ = 0.447: in a monatomic gas even the strongest shock leaves a
faster outflow.""")
D("D15", "The density and temperature jumps", ref="15.41", notes=r"""**N56 [B]** «15.41» (steps 1–3) · **N57 [B]** «15.42»
(steps 4–7)""", after=r"""As $M_1\to\infty$ the density ratio tends to $(\gamma+1)/(\gamma-1)$ = 6: **no shock squeezes air more
than sixfold.**""")
worked("a Mach 2 shock by hand", r"""
γ = 1.4, $M_1$ = 2.
1. «15.39»: $p_2/p_1$ = 1 + (2.8/2.4)(4 − 1) = **4.5**.
2. «15.41»: $\rho_2/\rho_1$ = 2.4 × 4/(0.4 × 4 + 2) = 9.6/3.6 = **2.667**.
3. Gas law: $T_2/T_1$ = 4.5/2.667 = **1.6875** («15.42» gives the same).
4. «15.40»: $M_2^2$ = (1.6 + 2)/(11.2 − 0.4) = 3.6/10.8 = 1/3, $M_2$ = **0.577**.
5. In sea-level air: 681 m/s in, 255 m/s out; 288 K in, 486 K out.
6. $T_0$: 288 × 1.8 = 519 K before; 486 × (1 + 0.2/3) = 519 K after ✓.""")
code(r"""
s = GD.normal_shock(2.0)                                     # every jump at M1 = 2 (ratios downstream / upstream)
print({k: round(float(s[k]), 4) for k in ("M2", "p_ratio", "rho_ratio", "T_ratio", "p0_ratio", "T0_ratio", "dS_over_cv")})
state1 = dict(p=1.0, rho=1.4, u=2.0)                         # upstream in units where p1 = 1 and c1 = sqrt(1.4 × 1/1.4) = 1
state2 = dict(p=s["p_ratio"], rho=1.4 * s["rho_ratio"], u=2.0 * s["u_ratio"])   # downstream from the ratios
print(GD.normal_shock_residuals(state1, state2))             # (15.35)-(15.37): relative residuals ~1e-16
print({k: round(float(v), 3) for k, v in GD.normal_shock_fluxes(2.0).items()})   # each flux out/in = 1, and how it is shared
try:                                                         # a "shock" from M1 = 0.7 …
    GD.normal_shock(0.7)
except ValueError as err:                                    # … is refused (trap T6)
    print("ValueError:", err)
print(GD.normal_shock(0.7, allow_expansion=True)["admissible"])   # computed on request, but flagged as not admissible: False
nss = functools.lru_cache(maxsize=None)(ch15.normal_shock_sympy)   # cache the slow symbolic derivation (runs once per session)
print(nss()["p_ratio"])                                      # sympy's own (15.39)
""", explain=r"""
1. `GD.normal_shock(2.0)` returns all the jumps: $M_2$ = 0.5774, $p_2/p_1$ = 4.5, $\rho_2/\rho_1$ = 2.667, $T_2/T_1$ = 1.6875,
   $p_{02}/p_{01}$ = 0.7209, $T_{02}/T_{01}$ = 1 and $(S_2-S_1)/C_v$ = 0.1309.
2. Two state dicts in units where $p_1=1$, $c_1=1$ (so $u_1=M_1$) feed `GD.normal_shock_residuals`: mass, momentum and
   energy balance to round-off.
3. `GD.normal_shock_fluxes` shows how each conserved flux is split — the momentum flux is mostly $\rho u^2$ before the shock
   and mostly $p$ after it.
4. `try … except ValueError` catches the error the function raises for $M_1<1$; with `allow_expansion=True` it computes the
   jump but reports `admissible = False` — C08 says why.
5. `ch15.normal_shock_sympy()` re-derives «15.39» symbolically; `functools.lru_cache` remembers the answer so it runs once.""")
scratch(r"""
g, M1 = 1.4, 2.0                                             # gamma; upstream Mach number
a = 1 / (g * M1**2)                                          # D13 step 6: a = p1/(rho1 u1²) = 1/(gamma M1²)
energy = lambda r: g/(g-1)*a + 0.5 - g/(g-1)*r*(a + 1 - r) - 0.5*r*r   # D13 step 6: energy law in r = u2/u1, left - right
r = brentq(energy, 0.05, 0.999)                              # the root that is not r = 1 (the bracket excludes it)
p_ratio = 1 + (1 - r) / a                                    # D13 step 12: back into momentum
assert np.allclose([1 / r, p_ratio], [s["rho_ratio"], s["p_ratio"]])   # same density and pressure jumps as the library
print(f"u2/u1 = {r:.4f}, rho2/rho1 = {1/r:.4f}, p2/p1 = {p_ratio:.4f}")
""", explain=r"""**From scratch, without the closed form.** The energy law in the one unknown $r=u_2/u_1$ (D13 step 6) is
solved numerically by `brentq` on a bracket that leaves out the trivial root r = 1; momentum then gives the pressure. Same
numbers as `GD.normal_shock`.""")
plotly(r"""
M1 = np.linspace(0.6, 5.0, 150)                              # upstream Mach numbers, including the forbidden half
gammas = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 5 / 3]                         # gamma slider
def jumps(g):                                                # the four jumps for one gamma
    return {"M₂, (15.40)": (M1, GD.normal_shock_mach(M1, g)),   # downstream Mach number
            "p₂/p₁ ÷ 10, (15.39)": (M1, GD.normal_shock_pressure_ratio(M1, g) / 10),   # pressure jump, scaled to share the axis
            "ρ₂/ρ₁, (15.41)": (M1, GD.normal_shock_density_ratio(M1, g)),   # density jump
            "T₂/T₁, (15.42)": (M1, GD.normal_shock_temperature_ratio(M1, g)),   # temperature jump
            "density limit (γ+1)/(γ−1)": ([0.6, 5.0], [(g + 1) / (g - 1)] * 2),   # the strong-shock limit of (15.41)
            "M₁ = 1 (left: forbidden, C08)": ([1.0, 1.0], [0.0, 9.0])}   # a vertical marker at M1 = 1
fig = slider_figure(jumps, "γ", gammas, xlabel="upstream Mach number M₁ [–]", ylabel="downstream / upstream [–]",   # slider over gamma
                    title="All the jumps of a normal shock follow from M₁ alone", yrange=[0, 9],   # the message
                    active=gammas.index(1.4))                # open on air, gamma = 1.4
fig.show()   # display
""", explain=r"""**Figure F5a.** Slider: γ, from 1.1 to 5/3. Ratios are downstream / upstream; the pressure ratio is divided by 10 to
share the axis.""")
notes3(see=r"""The four jumps against $M_1$ for the chosen γ, with the density limit and the line $M_1=1$.""",
       read=r"""At $M_1=1$ every ratio is 1 (no shock). Pressure and temperature keep growing like $M_1^2$; density levels
off at $(\gamma+1)/(\gamma-1)$ — 6 for air, 4 for a monatomic gas. Left of $M_1=1$ the formulas still give numbers (an
"expansion shock") — C08 forbids them.""",
       change=r"""…γ → 1 (a very complex molecule)? The density limit $(\gamma+1)/(\gamma-1)$ runs away: such gases can be
squeezed enormously by a strong shock.""")
note("N64 [B]", r"""**The same shock, seen from the still gas** (result given, not derived): subtract $u_1$ everywhere. The
shock runs into still gas at $u_{\text{shock}}=M_1c_1$ and drags the gas behind it at $u_{\text{gas}}=u_1-u_2$. Tiny
example: a shock that doubles the pressure of sea-level air — «15.39» gives $M_1^2$ = 1 + (2.4/2.8) × 1 = 1.857, $M_1$ =
1.363; shock speed 1.363 × 340.3 = 464 m/s; wind behind it 464 × (1 − 1/1.625) = 178 m/s.""")
figure(r"""
ms = GD.moving_shock(101325.0, 288.15, p2_over_p1=2.0, gas=air)   # a shock that doubles the pressure of sea-level air
print({k: round(float(ms[k]), 2) for k in ("M_s", "u_shock", "u_gas", "T2", "T0_behind_lab", "p0_behind_lab")})
print("T0 in the shock's frame (both sides):", round(288.15 * GD.stagnation_temperature_ratio(ms["M_s"]), 1), "K")   # 395.2 K
fig = draw_shock_two_frames(ms)                              # the same shock at rest and running into still gas (arrows to scale)
plt.show()                                                   # display
""", see=r"""The printed numbers ($M_s$ = 1.363, shock speed 463.7 m/s, wind behind it 178.4 m/s, $T_2$ = 354.6 K) and a
two-panel sketch: left, the shock at rest with gas entering at $u_1$ and leaving at $u_2$; right, the same shock running into
still air (front speed 464 m/s, leftward) with the gas behind following at $u_1-u_2$.""",
    read=r"""Each arrow is drawn to scale. The right panel is the left panel with $u_1$ subtracted from every velocity — a
Galilean shift (Ch. 3).""",
    change=r"""…the pressure ratio were 10 (a blast)? $M_s$ = 2.95, and the wind behind the front would itself be supersonic
relative to the ground.""")
confusion(r"""(trap T7) What changes with the frame and what does not:

| Same in both frames | Different |
|---|---|
| p, ρ, T, c, S | u, M, $T_0$, $p_0$ |

For the pressure-doubling shock above, $T_0$ is 395 K on both sides in the shock's frame (the energy law holds there); in
the frame of the still gas it is 288 K ahead of the front and 370 K behind it (`T0_behind_lab`). Stagnation quantities
belong to a frame.""")
recap("R17", "Strong shocks and the blast wave", r"""Ch. 1 found $R(t)\propto(E/\rho_1)^{1/5}t^{2/5}$ by dimensional
analysis. **New here:** for $M_1^2\gg1$, «15.41» tends to $(\gamma+1)/(\gamma-1)$ and «15.39» to $2\gamma M_1^2/(\gamma+1)$
(results given, not derived); a deliberately **crude** energy argument then gives a constant K in $R=K(E/\rho_1)^{1/5}t^{2/5}$.

> ⚠️ (trap T22) **crude**: all the energy as uniform internal energy behind the front. The true similarity solution has a
non-uniform interior; the crude sphere constant is about 9 % below the similarity value (qualitative only). And the
normalisation differs from Ch. 1: `ch01.blast_radius` writes $R=(Et^2/(K\rho))^{1/5}$, so its K is this chapter's K to the
power −5.""", where="Ch. 1 Example 1.4 (C69)", core_id="C07")
code(r"""
print(GD.strong_shock_limits())                              # rho ratio -> 6, M2 -> 0.378, p ratio -> 2 gamma M1²/(gamma+1)
K_h, K_s = ch15.blast_constant_crude(), ch15.blast_constant_crude(geometry="sphere")   # crude constants: hemisphere 1.07, sphere
print(f"crude K: hemisphere {K_h:.3f}, sphere {K_s:.3f}")
E, t, rho = 4.2e12, 0.01, 1.225                              # 1 kt of TNT [J], 10 ms, sea-level air
print(K_s * (E / rho)**0.2 * t**0.4, ch01.blast_radius(E, t, rho, K=K_s**-5))   # the same radius both ways [m]: Ch. 1's K = K**-5
m_ = sp.symbols("M_1", positive=True)                        # a symbol for M1
print(sp.limit((sp.Rational(12, 5) * m_**2) / (sp.Rational(2, 5) * m_**2 + 2), m_, sp.oo))   # (15.41) for gamma = 7/5 as M1 -> oo: 6
""", explain=r"""
1. `GD.strong_shock_limits()` gives the limits of the jump relations as $M_1\to\infty$.
2. `ch15.blast_constant_crude` is the crude constant (1.07 for a hemisphere); the two radius computations agree once Ch. 1's
   K is set to this K to the power −5.
3. `sp.limit` confirms the sixfold density limit symbolically.""")
whatif(r"""…$M_1$ were below 1? The three laws still balance — step 8 of D13 never asked. What forbids it is the second law.
C08.""")

# ---------------------------------------------------------------------------------------------------------------------
# C08 — the entropy jump and the second law
# ---------------------------------------------------------------------------------------------------------------------
core("C08", r"Why a shock only compresses: $\dfrac{S_2-S_1}{C_v}=\ln\left[\dfrac{p_2}{p_1}\left(\dfrac{\rho_1}{\rho_2}\right)^{\gamma}\right]$ (15.43)",
     "Mass, momentum and energy allow a jump in either direction. Which law chooses — and what exactly is lost?")
problem(r"""
Film a traffic jam and play it backwards: slow dense cars suddenly spread out and speed up, with no reason. Nothing in the car
count forbids it — it just never happens. A shock run backwards (subsonic gas jumping to supersonic) conserves mass, momentum
and energy too. What it would do is *destroy* entropy, and the second law forbids that. The same bookkeeping tells us what a
real shock costs: not energy — that is conserved — but stagnation pressure, the ability to do work.""")
idea("", r"""
| Conserved across a shock (flat, amber) | Not conserved |
|---|---|
| ṁ, momentum flux + p, $h_0$, $T_0$ | S rises (purple), $p_0$ falls, $A^*$ grows |""")
recap_md("Entropy and the Gibbs relation", "Ch. 1 C33–C35", r"""Entropy is a state function; the first law with the Gibbs
relation $T\,dS=dh-dp/\rho$; for a perfect gas the two logarithmic forms «15.1s».""")
remind("C08")
D("D16", "The entropy jump across a shock", ref="15.43")
note("N59 [B]", r"""The entropy jump against $M_1$, and its third-order start.""")
figure(r"""
M1 = np.linspace(0.5, 3.0, 300)                              # upstream Mach numbers on both sides of 1
dS = GD.normal_shock_entropy(M1)                             # (15.43) with the jumps inserted (D16 step 5)
cubic = GD.weak_shock_entropy(M1, form="mach")               # (15.44a): the cubic of D17
fig, ax = plt.subplots(figsize=(6.8, 3.8))                   # main panel
ax.axvspan(0.5, 1.0, color=COLORS["muted"], alpha=0.15)      # entropy would fall here: forbidden
ax.plot(M1, dS, color=COLORS["accent"], lw=2, label="exact (15.43)")      # entropy jump (purple)
ax.plot(M1, cubic, color=COLORS["accent"], ls="--", lw=1.2, label="cubic (15.44a)")   # its weak-shock form
ax.axhline(0, color=COLORS["ink"], lw=0.8)                   # zero
ax.set(xlim=(0.5, 3.0), ylim=(-0.15, 0.6), xlabel="upstream Mach number M₁ [–]", ylabel="(S₂ − S₁)/C_v [–]",
       title="Entropy rises only for M₁ > 1 — and starts at third order")
ax.legend(fontsize=8, loc="upper left")
ins = ax.inset_axes([0.62, 0.12, 0.34, 0.4])                 # a small log–log inset
eps = M1[M1 > 1.001]**2 - 1                                  # epsilon = M1² - 1 > 0
ins.loglog(eps, GD.normal_shock_entropy(np.sqrt(1 + eps)), color=COLORS["accent"])   # exact
ins.loglog(eps, 0.0648 * eps**3, color=COLORS["muted"], ls="--")   # slope 3 reference (coefficient of (15.44a))
ins.set_xlabel("M₁² − 1", fontsize=7); ins.set_ylabel("ΔS/C_v", fontsize=7); ins.tick_params(labelsize=6)
ins.set_title("slope 3 on log axes", fontsize=7)
savefig(fig, "ch15", "c08_entropy"); plt.show()
""", see=r"""The entropy jump of a normal shock against $M_1$ (solid), the cubic «15.44a» (dashed), the forbidden half shaded,
and an inset with both axes logarithmic (each tick a factor 10).""",
    read=r"""The curve crosses zero at $M_1=1$ with zero slope *and* zero curvature — it starts like a cube. Left of 1 it is
negative: an expansion shock would destroy entropy. In the inset the exact curve runs parallel to a line of slope 3: a factor
10 in $M_1^2-1$ gives a factor 1000 in entropy.""",
    change=r"""…we kept only the first two orders of the expansion? We would get zero — the reason D17 must go to third
order.""")
note("N61 [B]", r"""Reading «15.43»: $S_2-S_1>0\iff M_1>1$. So $M_1>1$, $M_2<1$, compression only (for a perfect gas;
expansion shocks need another equation of state — N90 in C12). Weak shocks are nearly isentropic, which is why even loud
sound is still acoustics (C01).""")
code(r"""
print(GD.shock_admissible(2.0), GD.shock_admissible(0.7))    # True, False: the second law's verdict
print(GD.normal_shock_entropy(0.7))                          # -0.0295: entropy would fall
""")
P("P364", "the series of ln(1 + x) to third order", r"""$\ln(1+x)=x-\tfrac12x^2+\tfrac13x^3-\dots$ for small x. When several
logarithms are added, their x and x² terms can cancel *exactly*, and the answer is then of third order — you must keep every
x³.""", code=r"""
import sympy as sp
x = sp.symbols('x')
print(sp.series(sp.log(1 + x), x, 0, 4))                    # x - x**2/2 + x**3/3 + O(x**4)
print(sp.series(sp.log(1 + 2*x) - 2*sp.log(1 + x), x, 0, 4))   # -x**2 + 2*x**3 + ...: first order cancelled
""")
D("D17", "The entropy jump of a weak shock", ref="15.44a", notes=r"""**N60 [B]** «15.44a» and «15.44b»""")
code(r"""
x = 0.1                                                      # a 10 % pressure jump, (p2 - p1)/p1
M1x = np.sqrt(1 + x * 2.4 / 2.8)                             # its Mach number from (15.39)
cub, ex = GD.weak_shock_entropy(x, form="pressure"), GD.normal_shock_entropy(M1x)   # (15.44b) and the exact (15.43)
print(f"10 % jump: cubic {cub:.2e}, exact {ex:.2e} (units of C_v); p0 lost {1 - np.exp(-ex / 0.4):.1e}")   # 4.1e-5, 3.5e-5, 8.8e-5
rows = [(M, GD.normal_shock_entropy(M), GD.weak_shock_entropy(M, form="mach")) for M in (1.05, 1.1, 1.2, 1.5)]
print(pd.DataFrame(rows, columns=["M1", "exact (15.43)", "cubic (15.44a)"]).to_string(index=False, float_format="%.2e"))
wss = functools.lru_cache(maxsize=None)(ch15.weak_shock_series)   # cache the slow series (runs once)
print(wss()["coefficients"])                                 # sympy's coefficients of eps, eps², eps³, eps⁴
""", explain=r"""
1. A 10 % pressure jump makes $3.5\times10^{-5}C_v$ of entropy (the cubic says $4.1\times10^{-5}$) and costs $8.8\times10^{-5}$
   of the stagnation pressure — nine parts in a hundred thousand.
2. The table compares exact and cubic: at $M_1$ = 1.05, $5.9\times10^{-5}$ against $7.0\times10^{-5}$; at 1.1, $4.3\times10^{-4}$
   against $6.0\times10^{-4}$; at 1.5 the cubic is four times too large. **The cubic is a local statement.**
3. `ch15.weak_shock_series()` lets sympy expand «15.43» in ε: the coefficients of ε and ε² are exactly zero.""")
D("D18", "The stagnation-pressure ratio across a shock", ref="o_p02", notes=r"""**N62 [B]** «o_p02» and the growth of
$A^*$ (step 7)""")
confusion(r"""(trap T8, second half) $A^*$ is not conserved across a shock — it grows by $p_{01}/p_{02}$.""")
note("N58 [B]", r"""The relations «15.39»–«15.42» are the **Rankine–Hugoniot relations**. (Ours:) eliminate $M_1$ between
«15.39» and «15.41» to get the Hugoniot curve «o_hug».""")
figure(r"""
R = np.linspace(1.0, 5.9, 300)                               # density ratios rho2/rho1 up to just below the limit 6
fig, ax = plt.subplots(figsize=(6.2, 3.8))                   # one panel
ax.loglog(R, GD.hugoniot_pressure_ratio(R), color=COLORS["orange"], lw=2, label="Hugoniot (shock), ours")   # pressure after a shock
ax.loglog(R, R**1.4, color=COLORS["muted"], ls="--", label="isentrope p ∝ ρ^γ")   # pressure after a smooth squeeze
ax.plot(2.667, GD.hugoniot_pressure_ratio(2.667), "o", color=COLORS["orange"])   # the Mach 2 shock
ax.plot(2.667, 2.667**1.4, "o", color=COLORS["muted"])       # the same squeeze, isentropic
ax.set(xlabel="ρ₂/ρ₁ (downstream / upstream) [–]", ylabel="p₂/p₁ [–]",
       title="A shock needs more pressure than a smooth squeeze to compress the same")
ax.legend(fontsize=8); savefig(fig, "ch15", "c08_hugoniot"); plt.show()
print(GD.hugoniot_pressure_ratio(2.667), 2.667**1.4)         # 4.50 against 3.95
""", see=r"""The pressure ratio needed for a given density ratio: by a shock (orange) and by a smooth isentropic squeeze
(grey dashed), on logarithmic axes (each tick a factor 10). The dots are the Mach 2 shock: 4.50 against 3.95.""",
    read=r"""The two curves touch at (1, 1) and stay together for weak compressions — they share value, slope and curvature
there, so they agree through second order and part only at third order in the jump; that third-order gap is the entropy of
D17 seen from another side. For strong ones the Hugoniot climbs steeply and runs to infinity at ρ₂/ρ₁ = 6.""",
    change=r"""…γ were 5/3? The vertical asymptote would move to ρ₂/ρ₁ = 4.""")
code(r"""
gs, ms, Rs = sp.symbols("gamma m R", positive=True)          # gamma; m = M1²; R = rho2/rho1
m_of_R = sp.solve(sp.Eq(Rs, (gs + 1) * ms / ((gs - 1) * ms + 2)), ms)[0]   # (15.41) solved for M1²
p_of_R = 1 + 2 * gs / (gs + 1) * (m_of_R - 1)                # put into (15.39)
print(sp.simplify(p_of_R - ((gs + 1) * Rs - (gs - 1)) / ((gs + 1) - (gs - 1) * Rs)))   # 0: the Hugoniot curve
""")
note("N63 [B]", r"""The book's shock table is not reproduced; `GD.normal_shock_table` generates ours (a DataFrame, primer
P361 in C05), checked against the public NACA Report 1135 rows. (Slip #10 concerns the last digit of one entry of the book's
shock table.)""")
code(r"""
print(GD.normal_shock_table(np.arange(1.0, 5.01, 0.5)).round(4).to_string(index=False))   # our table; row M1 = 2: 0.5774, 4.5, …, 0.7209
""")
recap("R16", "The hydraulic jump", r"""«7.81» — supercritical ($Fr_1>1$) water jumps to subcritical.

| Hydraulic jump (Ch. 7) | Normal shock |
|---|---|
| Froude number Fr | Mach number M |
| supercritical → subcritical | supersonic → subsonic |
| depth rises | density rises |
| mechanical energy lost | stagnation pressure lost |""", where="Ch. 7 §7.6 (C11–C12)", core_id="C08")
code(r"""
print(ch07.hydraulic_jump(1.0, Fr1=2.0)["ratio"])            # (7.81): H2/H1 = 2.372 for Fr1 = 2
""")
confusion(r"""(trap T21) The analogy is partial — the jump loses mechanical energy to heat; the shock conserves energy and
loses $p_0$.""")
worked("what a Mach 2 shock costs", r"""
From C07: $p_2/p_1$ = 4.5, $\rho_2/\rho_1$ = 2.667.
1. «15.43»: $(S_2-S_1)/C_v=\ln[4.5\times(1/2.667)^{1.4}]=\ln(4.5\times0.2533)=\ln1.140$ = **0.131**.
2. In units of R: 0.131/(γ − 1) = 0.327.
3. Stagnation pressure kept: $e^{-0.327}$ = **0.721** — 28 % lost.
4. Cross-check with «o_p02»: (2.667)^3.5 × (2.4/10.8)^2.5 = 30.97 × 0.02328 = 0.721 ✓.
5. $A^*$ behind the shock is 1/0.721 = 1.387 times larger.""")
code(r"""
print(GD.normal_shock_entropy(2.0), GD.normal_shock_stagnation_pressure_ratio(2.0), np.exp(-0.13091 / 0.4))   # 0.1309, 0.7209, 0.7209
""", explain=r"""Three numbers: the entropy jump, the stagnation-pressure ratio from «o_p02», and $e^{-\Delta S/R}$ — the last
two agree, as D18 step 6 says they must.""")
scratch(r"""
T1, p1 = 288.15, 101325.0                                    # an upstream state [K, Pa]
dS = TH.perfect_gas_entropy_change_p(T1, p1, T1 * s["T_ratio"], p1 * s["p_ratio"], cp=air.cp, R=air.R)   # (15.1, entropy line) [J/(kg K)]
p0_ratio = s["p_ratio"] * GD.stagnation_pressure_ratio(s["M2"]) / GD.stagnation_pressure_ratio(2.0)   # D18 step 1: chain of three ratios
assert np.allclose(dS / air.cv, GD.normal_shock_entropy(2.0))                     # route 1: entropy from the Ch. 1 formula
assert np.allclose([p0_ratio, np.exp(-dS / air.R)], GD.normal_shock_stagnation_pressure_ratio(2.0))   # routes 2 and 3
print(f"dS = {dS:.2f} J/(kg K); p02/p01 = {p0_ratio:.4f} = exp(-dS/R) = {np.exp(-dS/air.R):.4f}")
""", explain=r"""**From scratch — three routes, one number.** The entropy change from the Ch. 1 formula (`TH`), the
stagnation-pressure ratio as a chain of three static/stagnation ratios, and $e^{-\Delta S/R}$: all equal the library's.""")
plotly(r"""
M1 = np.linspace(0.6, 5.0, 150)                              # upstream Mach numbers
gammas = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 5 / 3]                         # gamma slider
fig = slider_figure(lambda g: {"(S₂ − S₁)/C_v, (15.43)": (M1, GD.normal_shock_entropy(M1, g)),   # entropy jump
                               "p₀₂/p₀₁ (ours)": (M1, GD.normal_shock_stagnation_pressure_ratio(M1, g)),   # surviving stagnation pressure
                               "T₀₂/T₀₁ = 1 (conserved)": (M1, np.ones_like(M1)),   # flat: energy is conserved
                               "M₁ = 1 (left: forbidden)": ([1.0, 1.0], [-0.3, 1.6])},   # a vertical marker
                    "γ", gammas, xlabel="upstream Mach number M₁ [–]", ylabel="[–]",   # slider over gamma
                    title="What a shock costs: entropy up, stagnation pressure down, T₀ untouched", yrange=[-0.3, 1.6],   # the message
                    active=gammas.index(1.4))                # open on air, gamma = 1.4
fig.show()   # display
""", explain=r"""**Figure F5b.** Slider: γ. The flat line is the conserved stagnation temperature.""")
notes3(see=r"""The entropy jump, the surviving fraction of stagnation pressure and the (flat) stagnation-temperature ratio
against $M_1$, with the line $M_1=1$.""",
       read=r"""Right of the line the entropy rises and the stagnation pressure falls together: at $M_1$ = 2 for air, 28 %
lost; at 3, 67 %. Left of it the "shock" would *gain* stagnation pressure — forbidden. $T_0$ never moves: energy is
conserved.""",
       change=r"""…you wanted to slow a Mach 3 stream with small losses? Do it in steps (oblique shocks, C14): several weak
shocks lose far less than one strong one, because the loss starts at third order.""")
explainer("normal_shock_jump", "Why can a shock only compress?",
          r"""you drag the upstream Mach number below 1: the three flux bars stay perfectly balanced while the entropy bar goes
negative — you find the second law's veto by experiment rather than by reading it.""",
          ["Drag $M_1$ from 2 down to 0.7: which bar breaks?",
           "Sit at $M_1$ = 1.1: how much stagnation pressure is lost?",
           "Switch to the frame of the still gas: which numbers changed, which did not?"])
recap("R18", "Viscosity sets the thickness, not the loss", r"""In turbulence the dissipation rate is set by the large scales,
not by the viscosity. Likewise the entropy rise of a shock is fixed by the jump conditions; viscosity and conduction only set
how thick the layer is (N65, N69 below).""", where="Ch. 12 §12.7 (C07)", core_id="C08")
note("N65 [B]", r"""**Thickness estimate:** a Reynolds number of order one on the velocity jump,
$(u_1-u_2)\,\delta_{\text{sh}}/\nu_{\text{visc}}\sim1$ (δ is a thickness here, a deflection angle in §15.9 — trap T5).""")
code(r"""
lam = ch01.mean_free_path_air(288.15)                        # mean free path of sea-level air [m]
for M in (2.0, 1.2):                                         # a strong and a weak shock
    u1 = M * c0; u2 = u1 / GD.normal_shock_density_ratio(M)  # speeds before and after [m/s]
    d = ch15.shock_thickness_estimate(u1 - u2, 1.46e-5)      # delta ~ nu/(u1 - u2) [m]
    print(f"M1 = {M}: du = {u1-u2:.0f} m/s, thickness ~ {d:.1e} m = {d/lam:.1f} mean free paths")
""", explain=r"""At $M_1$ = 2 the estimate is 3.4 × 10⁻⁸ m, at 1.2 about 1.4 × 10⁻⁷ m — of the order of the mean free path
(≈ 7 × 10⁻⁸ m). The continuum model is at its limit (the Knudsen number of Ch. 1 C06 is of order one). Our numbers, not the
book's.""")
recap("R19", "The equations inside a shock", r"""Steady, one-dimensional mass, momentum and energy with viscosity and
conduction: $\dfrac{d(\rho u)}{dx}=0$, $\rho u\dfrac{du}{dx}+\dfrac{dp}{dx}=\dfrac{d}{dx}\left(\mu''\dfrac{du}{dx}\right)$,
$\rho u\dfrac{dh}{dx}-u\dfrac{dp}{dx}=\mu''\left(\dfrac{du}{dx}\right)^2+\dfrac{d}{dx}\left(k\dfrac{dT}{dx}\right)$ with
$\mu''=\tfrac43\mu+\mu_v$ (reductions of the Ch. 4 equations).""", where="Ch. 4 (C08, C10; N51 bulk viscosity)", core_id="C08")
note("N66 [B]", r"""Each equation integrates once (result given, not derived): $\rho u=m$,
$mu+p=\mu''\dfrac{du}{dx}+mV$, $m\left(h+\tfrac12u^2\right)=\mu''u\dfrac{du}{dx}+k\dfrac{dT}{dx}+mI$. The three constants
of integration are fixed by the upstream state: **m** $=\rho_1u_1$ is the mass flux [kg/(m² s)],
**V** $=u_1+p_1/m$ is a speed [m/s] (momentum flux over mass flux), and **I** $=C_pT_1+\tfrac12u_1^2$ is the total
enthalpy [J/kg]. Far upstream and downstream the gradients vanish and these *are* «15.35», «15.36», «15.37»: the jump
conditions do not care what happens inside.""")
note("N67 [B]", r"""With $Pr''\equiv\mu''C_p/k=1$ the total enthalpy $C_pT+\tfrac12u^2=I$ at every point inside the shock
(Becker; result given). Trap T20: $Pr''$ is built on $\mu''$ — about 0.95 for air with $\mu_v=0$.""")
note("N68 [B]", r"""One first-order equation remains for $U=u/V$; its two rest points are the end states, (ours)
$U_1=\gamma M_1^2/(\gamma M_1^2+1)$ and the same with $M_2$; $U_1+U_2=2\gamma/(\gamma+1)$. (Slip #7: the book factorises
before dividing by $V^2$.)""")
code(r"""
sss = ch15.shock_structure_sympy()                           # the reductions and their first integrals, by sympy
print(sss["structure_equation"])                             # residual of the reduction to one equation for U = u/V: 0 = exact
sc = ch15.shock_structure_constants(2.0, 288.15, 101325.0, air)   # the integration constants for a Mach 2 shock in sea-level air
print({k: round(float(sc[k]), 4) for k in ("m", "V", "I", "U1", "U2")})   # m [kg/(m² s)], V [m/s], I [J/kg], U1, U2 [–]
""", explain=r"""
1. `ch15.shock_structure_sympy()` redoes the reductions of N66–N68 symbolically and returns each as an expression that must
   vanish. The printed **0** is the residual of the last one — the single first-order equation for $U=u/V$: sympy confirms the
   reduction leaves exactly nothing over.
2. `ch15.shock_structure_constants` evaluates the three constants of N66 for a Mach 2 shock in sea-level air:
   m = 833.7 kg/(m² s), V = 802.1 m/s, I = 521 095 J/kg, and the two end states of N68, $U_1$ = 0.8485 and $U_2$ = 0.3182
   (28/33 and 7/22).""")
code(r"""
U1, U2 = ch15.becker_shock_roots(2.0)                        # the two end states of the scaled speed for M1 = 2
print(U1, U2, U1 + U2, 28 / 33, 7 / 22, 2 * 1.4 / 2.4)       # (28/33, 7/22); sum 2 gamma/(gamma + 1) = 1.1667
""")
note("N69 [B]", r"""The profile «u15.6b», inverted numerically, for three shock strengths. **N70 [B]** Our figure uses our own
Mach numbers; for γ = 1.4, $M_1$ = 2 the plateaus are 28/33 and 7/22 (slip #8 concerns a caption of the book's figure; the
comparison is a private test).""")
figure(r"""
eta = np.linspace(-6.0, 6.0, 200)                            # shock coordinate (scaled distance; origin arbitrary — trap T20)
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.6))         # profiles; thickness
for M, col in ((1.2, COLORS["teal"]), (2.0, COLORS["accent"]), (4.0, COLORS["rose"])):
    U = ch15.becker_shock_profile(eta, M)                     # the implicit profile solved for U at each eta (brentq per point)
    a.plot(eta, U, color=col, label=f"M₁ = {M}")              # scaled speed through the shock
    sol = solve_ivp(lambda e, y: ch15.becker_shock_ode(e, y, M), (eta[0], eta[-1]), [U[0]],
                    t_eval=eta[::20], rtol=1e-9, atol=1e-12)  # the same profile by integrating the ODE (independent check)
    a.plot(sol.t, sol.y[0], "o", ms=3, color=col)             # dots on the curve
a.set(xlabel="η (scaled distance) [–]", ylabel="U = u/V [–]", title="Inside a shock: a smooth step (dots: ODE)")
a.legend(fontsize=8)
Ms = np.linspace(1.05, 5.0, 60)                              # Mach numbers
b.semilogy(Ms, ch15.becker_shock_thickness(Ms), color=COLORS["accent"])   # maximum-slope thickness (ours)
b.set(xlabel="M₁ [–]", ylabel="thickness Δη [–] (log axis)", title="Weak shocks are thick")
savefig(fig, "ch15", "c08_becker"); plt.show()
""", see=r"""Left: the scaled speed $U=u/V$ through a shock for $M_1$ = 1.2, 2 and 4 (curves from the implicit relation «u15.6b»,
dots from integrating the differential equation with `solve_ivp`). Right: the thickness of the step against $M_1$ on a
logarithmic axis.""",
    read=r"""Each profile goes smoothly from the upstream plateau $U_1$ to the downstream plateau $U_2$ (for $M_1$ = 1.2 the
step is wider than the window); for $M_1$ = 2 the plateaus are 28/33 = 0.848 and 7/22 = 0.318. The dots lie on the curves: two routes, one profile. The thickness grows without bound as
$M_1\to1$ — a weak shock is a gentle ramp.""",
    change=r"""…$Pr''$ were not 1? The energy integral would no longer be exact and the profile would need a second equation;
for $Pr''\to0$ a thinner shock appears inside the first (N71).""")
note("N71 [C]", r"""For $Pr''\to0$ a second, thinner shock sits inside the first; see the literature the book cites — not
developed in the book or here.""")
whatif(r"""…you put this jump inside a nozzle? Then the back pressure decides where it stands. First the simpler nozzle —
C09.""")

# =====================================================================================================================
# A.7 §15.7 Operation of Nozzles at Different Back Pressures — C09, C10
# =====================================================================================================================
nb.section("15.7", "Operation of Nozzles at Different Back Pressures", intro="**What is this section about?** A reservoir, "
           "a nozzle and a tank whose pressure you lower. The book describes what happens in words and figures; the section "
           "has no numbered equation, so **every formula here is ours**, built from the earlier sections.")
core("C09", r"Choking: $\dfrac{p^*}{p_0}=\left(\dfrac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}$ (ours)",
     "If you keep lowering the pressure downstream of a nozzle, does more and more gas flow?")
problem(r"""
Let the air out of a bicycle tyre. Open the valve a little into a room at lower pressure and air flows; lower the room
pressure and more flows. Intuition says: lower it to vacuum and the flow keeps growing. It does not. Once the pressure
outside is about half the pressure inside, the air at the exit moves at exactly the speed of sound — and news of any further
drop outside cannot travel upstream against it. The nozzle is **choked**.""")
idea(r"""
reservoir p₀, T₀ ══⟩ exit │ tank p_B
  exit subsonic:  exit pressure = p_B        (news of p_B gets in)
  exit sonic:     exit pressure = p* ≥ p_B   (news cannot get in)
""")
note("N72 [B]", r"""Three pressures: the supply (stagnation) pressure $p_0$, the back pressure $p_B$ in the tank, and the
exit-plane pressure. **The exit pressure equals the back pressure if and only if the exit flow is subsonic** — the range of
influence of C02 (N20) put to work.""")
figure(r"""
fig = draw_reservoir_nozzle_tank()                           # our sketch: reservoir (p0, T0), convergent nozzle, tank at p_B
plt.show()                                                   # display it
""", see=r"""A large reservoir at the stagnation state, a convergent nozzle, and a tank held at the back pressure $p_B$.""",
    read=r"""The reservoir is so large that the gas in it is at rest: it *is* the stagnation state (N34). Everything about
the nozzle flow is a ratio to $p_0$ and $T_0$; the only thing you control is $p_B$.""",
    change=r"""…the reservoir were small? Its pressure would fall as it emptied; each instant would still look like this
picture with a lower $p_0$ (a quasi-steady emptying).""")
D("D19", "Mass flow through a convergent nozzle against back pressure", ref="o_flux")
worked("our nozzle, three back pressures", r"""
$p_0$ = 500 kPa, $T_0$ = 400 K, exit area 10 cm², γ = 1.4; critical ratio «o_pstar» = 0.528.
1. $p_B$ = 400 kPa (ratio 0.8 > 0.528: not choked): invert «15.29»: $M_e^2$ = 5 × (0.8^(−0.2857) − 1) = 5 × 0.0658 = 0.329,
   $M_e$ = 0.574.
2. Flux function «o_flux» at 0.574: 1.183 × 0.574/1.0658³ = 0.561 — that is 82 % of its maximum 0.685.
3. ṁ = 0.561 × 500 000 × 0.001/√(287 × 400) = 0.83 kg/s.
4. $p_B$ = 264 kPa (ratio 0.528): $M_e$ = 1, ṁ = 1.01 kg/s.
5. $p_B$ = 100 kPa: still $M_e$ = 1 and ṁ = 1.01 kg/s; the exit pressure stays at 264 kPa and the jet finishes expanding
   outside.""")
code(r"""
print(GD.critical_pressure_ratio())                          # p*/p0 = 0.5283 for gamma = 1.4
pBs = (4.75e5, 4.0e5, 3.0e5, 2.6414e5, 1.0e5)                # five back pressures [Pa]
rows = [GD.convergent_nozzle(noz["p0"], noz["T0"], pB, noz["A_t"], air) for pB in pBs]   # D19 for each
tab = pd.DataFrame(rows, index=[f"p_B = {pB/1e3:.0f} kPa" for pB in pBs])[["M_exit", "p_exit", "mdot", "choked"]]
print(tab.round(3))                                          # M_exit 0.272 … 1; mdot 0.454 … 1.010 kg/s; choked from 264 kPa down
print(GD.mass_flow_vs_back_pressure(0.8))                    # dimensionless: mdot sqrt(R T0)/(p0 A) = 0.5607, 82 % of the choked value
""", explain=r"""
1. The critical ratio «o_pstar» is 0.5283: for our 5 bar supply, 264 kPa.
2. For each back pressure `GD.convergent_nozzle` returns the exit Mach number, exit pressure, mass flow and whether the
   nozzle is choked. Above 264 kPa the exit pressure *is* the back pressure; below it the exit pressure sticks at 264 kPa.
3. The mass flow rises from 0.45 kg/s (475 kPa) to 1.010 kg/s at 264 kPa — and not a gram more at 100 kPa.
4. `GD.mass_flow_vs_back_pressure` gives the same in dimensionless form: 0.5607 at $p_B/p_0$ = 0.8.""")
remind("C09")
scratch(r"""
g, R, p0, T0, A = 1.4, air.R, noz["p0"], noz["T0"], noz["A_t"]   # gas and nozzle
crit = (2 / (g + 1))**(g / (g - 1))                         # p*/p0 (ours)
pB_list = np.linspace(0.05, 0.99, 50) * p0                  # 50 back pressures [Pa]
mdot_mine = []                                               # filled in the loop
for pB in pB_list:                                           # the piecewise law of D19
    if pB / p0 > crit:                                       # subsonic exit: it feels the tank
        M = np.sqrt(2 / (g - 1) * ((pB / p0)**(-(g - 1) / g) - 1))   # D19 step 4: invert (15.29)
    else:                                                    # choked: the exit is sonic
        M = 1.0                                              # D19 step 6
    mdot_mine.append(np.sqrt(g) * M * (1 + 0.5 * (g - 1) * M * M)**(-(g + 1) / (2 * (g - 1))) * p0 * A / np.sqrt(R * T0))   # D19 step 3
assert np.allclose(mdot_mine, [GD.convergent_nozzle(p0, T0, pB, A, air)["mdot"] for pB in pB_list])   # same as the library
print(f"{len(pB_list)} back pressures: from scratch = GD.convergent_nozzle; plateau {max(mdot_mine):.4f} kg/s")
""", explain=r"""**From scratch.** An `if` decides the regime (subsonic exit or choked); then D19 step 3 gives the mass flow. Fifty
back pressures, the same numbers as the library.""")
note("N73 [B]", r"""The pressure along the nozzle for six back pressures, and the mass flow against the back pressure.""")
figure(r"""
x = np.linspace(0.0, 0.4, 81)                                # the converging half of laval_area: inlet at 0, exit (throat) at 0.4
fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.8))         # pressure along the nozzle; mass flow against back pressure
for r_, col in zip((0.95, 0.8, 0.6, 0.528, 0.3, 0.1), (COLORS["blue"], COLORS["teal"], COLORS["accent"], COLORS["amber"], COLORS["orange"], COLORS["rose"])):
    Me = GD.convergent_nozzle(1.0, 300.0, r_, 1.0, air)["M_exit"]   # exit Mach number for this back pressure (dimensionless supply)
    if Me < 0.9999:                                          # subsonic exit: sonic area below the exit area
        sol = GD.nozzle_isentropic(x, GD.laval_area, branch="subsonic", A_star=1.0 / GD.area_ratio(Me))
    else:                                                    # choked: sonic at the exit
        sol = GD.nozzle_isentropic(x, GD.laval_area, branch="transonic")
    a.plot(x, sol["p_over_p0"], color=col, label=f"p_B/p₀ = {r_}")   # pressure along the nozzle
    a.plot(0.43, r_, ">", color=col)                         # the tank pressure, just outside the exit
a.set(xlabel="x [m] (exit at 0.4)", ylabel="p/p₀ [–]", title="Below 0.528 the curves stop changing")
a.legend(fontsize=7)
rr = np.linspace(0.02, 1.0, 200)                             # back-pressure ratios
b.plot(rr, [GD.mass_flow_vs_back_pressure(v)["mdot_over_choked"] for v in rr], color=COLORS["accent"])   # mdot / mdot_max
b.axvline(GD.critical_pressure_ratio(), color=COLORS["amber"], ls="--", lw=1)   # the critical ratio
b.set(xlabel="back pressure p_B/p₀ [–]", ylabel="ṁ / ṁ_max [–]", title="The mass flow stops rising: choked")
savefig(fig, "ch15", "c09_convergent"); plt.show()
""", see=r"""Left: the pressure along a convergent nozzle for six back-pressure ratios (triangles: the tank pressure just
outside the exit). Right: the mass flow, as a fraction of its maximum, against the back-pressure ratio.""",
    read=r"""For 0.95, 0.8 and 0.6 each curve ends exactly at its triangle: the exit pressure is the back pressure. For 0.528,
0.3 and 0.1 the three curves lie on top of one another and end at 0.528 — above the triangles: the gas leaves at a higher
pressure than the tank's and finishes expanding outside. On the right the curve meets its plateau smoothly (D19 step 7).""",
    change=r"""…the gas were monatomic (γ = 5/3)? The critical ratio would be 0.487: it chokes a little later.""")
note("N74 [B]", r"""**Choking:** with the exit sonic the flow upstream no longer responds to $p_B$, and $\dot m=\dot m_{\max}$ of
«o_mmax». > ⚠️ Common confusion (trap T10): "choked" does not mean the flow cannot be raised — it is the most the throat
passes *for this supply*. Raise $p_0$ or lower $T_0$ and it rises ($\dot m\propto p_0/\sqrt{T_0}$).""")
plotly(r"""
def by_p0(p0):                                               # mass flow against back pressure for one supply pressure
    pB = np.linspace(0.01, 1.0, 60) * p0                     # back pressures up to the supply [Pa]
    m = [GD.convergent_nozzle(p0, noz["T0"], v, noz["A_t"], air)["mdot"] for v in pB]   # D19 at each back pressure
    return {"ṁ [kg/s]": (pB / 1e3, m)}                       # x in kPa
fig_a = slider_figure(by_p0, "p₀", np.linspace(2e5, 1e6, 17), unit="Pa", xlabel="back pressure p_B [kPa]",   # slider over the supply pressure
                      ylabel="mass flow ṁ [kg/s]", title="Raise the supply: the plateau rises (choked ≠ capped)")   # the message
fig_a.show()   # display
def by_T0(T0):                                               # the same for one supply temperature
    pB = np.linspace(0.01, 1.0, 60) * noz["p0"]
    m = [GD.convergent_nozzle(noz["p0"], T0, v, noz["A_t"], air)["mdot"] for v in pB]   # D19 at each back pressure
    return {"ṁ [kg/s]": (pB / 1e3, m)}   # x in kPa
fig_b = slider_figure(by_T0, "T₀", np.linspace(300.0, 1200.0, 19), unit="K", xlabel="back pressure p_B [kPa]",   # slider over the supply temperature
                      ylabel="mass flow ṁ [kg/s]", title="Heat the supply: the plateau falls like 1/√T₀")   # the message
fig_b.show()   # display
""", explain=r"""**Figures F6.** Two slider figures with dimensional axes: the supply pressure (200 … 1000 kPa) and the supply
temperature (300 … 1200 K).""")
notes3(see=r"""The mass flow of our 10 cm² nozzle against the back pressure, for the supply pressure (top) and temperature
(bottom) on the sliders.""",
       read=r"""Each curve rises from zero at $p_B=p_0$ and flattens at $p_B=0.528\,p_0$. Doubling $p_0$ doubles the plateau
and moves the kink to the right; quadrupling $T_0$ halves the plateau. The throat is a flow meter that reads the supply.""",
       change=r"""…the throat area were doubled? Every curve would double: ṁ ∝ $p_0A^*/\sqrt{T_0}$.""")
whatif(r"""…you added a diverging part after the throat? Then below the choking pressure the flow has somewhere to go:
supersonic. And the back pressure gets a new job. C10.""")

# ---------------------------------------------------------------------------------------------------------------------
# C10 — the Laval nozzle
# ---------------------------------------------------------------------------------------------------------------------
core("C10", r"The Laval nozzle: three dividing back pressures $\dfrac{p_b}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sub}}(A_e/A_t)},\ \dfrac{p_f}{p_0}=\left.\dfrac{p}{p_0}\right\rvert_{M_{\text{sup}}(A_e/A_t)},\ \dfrac{p_d}{p_0}=\dfrac{p_f}{p_0}\left.\dfrac{p_2}{p_1}\right\rvert_{M_{\text{sup}}}$ (ours)",
     "Once the throat is sonic, what is left for the back pressure to decide?")
problem(r"""
A rocket nozzle on the test stand at sea level, then the same nozzle at altitude. The chamber is the same, the throat is the
same, the mass flow is the same — but the exhaust looks different, and at start-up a sharp front can be seen standing inside
the bell. With the throat sonic the back pressure cannot change the mass flow. All it can do is place a normal shock
somewhere in the diverging part — or, when it is low enough, leave the inside of the nozzle alone and reshape the jet
outside.""")
idea("", r"""The ladder of seven conditions for our nozzle ($A_e/A_t$ = 2):

| Name | $p_B/p_0$ | What you see |
|---|---|---|
| subsonic throughout | > 0.937 | a venturi; the mass flow still depends on $p_B$ |
| just choked | = 0.937 = $p_b$ | sonic at the throat only |
| shock in the diverging part | 0.513 … 0.937 | a normal shock inside the bell; subsonic exit |
| shock at the exit | = 0.513 = $p_d$ | the shock sits in the exit plane |
| over-expanded | 0.094 … 0.513 | supersonic exit; oblique shocks outside |
| design | = 0.094 = $p_f$ | supersonic exit at the back pressure: a clean jet |
| under-expanded | < 0.094 | supersonic exit; expansion fans outside |""")
P("P365", "nested root finding", r"""Sometimes the function you want a root of itself needs a root-finder to evaluate.
Here: "exit pressure for a shock at Mach $M_1$" needs the subsonic inverse of the area ratio; then we search $M_1$ so that this
exit pressure equals the back pressure. One `brentq` inside another.""", code=r"""
from scipy.optimize import brentq
inner = lambda y: brentq(lambda x: x**3 + x - y, 0.0, 10.0)   # inner: solve x^3 + x = y for x
outer = brentq(lambda y: inner(y) - 1.5, 0.0, 100.0)          # outer: which y makes that x equal 1.5?
print(outer, 1.5**3 + 1.5)                                     # 4.875 both
""")
D("D20", "The Laval nozzle: dividing back pressures and the place of the shock", ref="o_shockpos",
  notes=r"""**N75 [B]** a normal shock inside the diverging part, placed by the back pressure (steps 5–9)""")
note("N78 [B]", r"""A worked shock with our inputs. (Slip #6 concerns the book's own example — lower-case m for a Mach number
and a ratio written as an equality of temperatures; ours uses our notation.)""")
worked("a shock where the area is 1.5 throats (our inputs)", r"""
$A_e/A_t$ = 2, $T_0$ = 400 K, a normal shock at $A/A_t$ = 1.5.
1. Ahead of the shock the flow is isentropic and supersonic with $A^*=A_t$: invert «15.31» → $M_1$ = 1.854.
2. Across the shock «15.40»: $M_2$ = 0.605; «o_p02»: $p_{02}/p_{01}$ = 0.788.
3. **The sonic area grows:** $A_2^*=A_t/0.788=1.268A_t$ (trap T8).
4. Exit: $A_e/A_2^*$ = 2/1.268 = 1.577 → **subsonic** inverse → $M_e$ = 0.404.
5. «15.28»: $T_e/T_0$ = 0.968 → $T_e$ = 387 K ($T_0$ is unchanged by the shock).
6. Exit pressure: $p_e/p_{01}$ = 0.788 × 0.894 = 0.704 — and since the exit is subsonic, that *is* the back pressure this shock
   position belongs to.""")
code(r"""
M1 = GD.mach_from_area_ratio(1.5, branch="supersonic")      # Mach number just ahead of a shock at A = 1.5 A_t -> 1.8541
w = ch15.nozzle_with_shock(M1, noz["A_t"], 2 * noz["A_t"], noz["p0"], noz["T0"], air)   # the whole chain of the example
print(round(M1, 4), {k: float(f"{w[k]:.4g}") for k in ("A_shock", "M2", "p02", "A_star2", "M3", "T3", "p3", "mdot")})
print("p02/p01 =", w["p02"] / noz["p0"], " p_exit/p01 =", w["p3"] / noz["p0"])   # 0.7884 and 0.7045
""", explain=r"""`ch15.nozzle_with_shock` runs the six steps above: $M_2$ = 0.6048, $p_{02}$ = 394.2 kPa, $A_2^*$ = 1.268 ×
10⁻³ m², exit Mach 0.4042, $T_3$ = 387.3 K, exit pressure 352.2 kPa = 0.7045 $p_{01}$, and the mass flow is still the choked
1.010 kg/s.""")
code(r"""
th = GD.laval_critical_pressures(2.0)                        # the three dividing back pressures for A_e/A_t = 2
print({k: round(float(th[k]), 4) for k in ("p_b", "p_d", "p_f", "M_sub", "M_sup")})   # 0.9372, 0.5134, 0.0939; 0.306, 2.197
for pB in (0.97, th["p_b"], 0.7, th["p_d"], 0.3, th["p_f"], 0.05):   # one back pressure per rung of the ladder
    print(f"p_B/p0 = {pB:.4f}: {GD.laval_regime(pB, 2.0)}")  # the three equalities need the exact thresholds
loc = GD.laval_shock_location(0.7, 2.0)                      # D20 steps 5-9 for p_B/p0 = 0.7
print({k: round(float(loc[k]), 4) for k in ("M1", "A_shock_over_At", "p02_over_p01", "M_exit", "p_exit_over_p0")})
print({k: v for k, v in GD.laval_exit_state(0.7, 2.0).items() if k in ("regime", "mdot_over_choked", "choked")})
print(GD.laval_shock_location(0.7045, 2.0)["A_shock_over_At"])   # feeding the worked example's exit pressure back: 1.5
""", explain=r"""
1. `GD.laval_critical_pressures(2.0)` gives the three dividing ratios of D20: $p_b$ = 0.9372, $p_d$ = 0.5134, $p_f$ = 0.0939.
2. `GD.laval_regime` names the regime of each back pressure — the seven rungs of the ladder.
3. For $p_B/p_0$ = 0.7, `GD.laval_shock_location` solves the nested problem: the shock stands where $M_1$ = 1.8627, at
   1.5101 throat areas; the flow keeps 78.45 % of its stagnation pressure and leaves at Mach 0.4067.
4. `GD.laval_exit_state` adds that the throat is choked and the mass flow is the maximum.
5. Feeding back the worked example's exit pressure, 0.7045, returns its shock at 1.5 throat areas.""")
scratch(r"""
def area_inv_sub(AR, n=80):                                  # subsonic inverse of (15.31) by bisection on [1e-6, 1]
    a, b = 1e-6, 1.0
    for _ in range(n):                                       # halve the bracket n times
        m = 0.5 * (a + b)
        if GD.area_ratio(m) > AR: a = m                      # A/A* too large: M too small
        else: b = m
    return 0.5 * (a + b)
def exit_p(M1):                                              # D20 steps 6-8: exit pressure for a shock at Mach M1
    r = GD.normal_shock_stagnation_pressure_ratio(M1)        # p02/p01
    Me = area_inv_sub(2.0 * r)                               # A_e/A2* = (A_e/A_t)(p02/p01)
    return r * (1 + 0.2 * Me * Me)**-3.5                     # p_e/p01 = (p02/p01)(p/p02 at Me)
lo, hi = 1.0 + 1e-9, GD.mach_from_area_ratio(2.0, branch="supersonic")   # M1 between 1 and M_sup
for _ in range(60):                                          # outer bisection: exit_p falls as M1 rises
    mid = 0.5 * (lo + hi)
    if exit_p(mid) > 0.7: lo = mid                           # shock too weak: move it downstream
    else: hi = mid
M1_mine = 0.5 * (lo + hi)                                    # the shock's upstream Mach number
assert np.allclose([M1_mine, GD.area_ratio(M1_mine)], [loc["M1"], loc["A_shock_over_At"]], atol=1e-8)   # same as the library
print(M1_mine, GD.area_ratio(M1_mine))                       # 1.8627, 1.5101
""", explain=r"""**From scratch.** A bisection inside a bisection (primer P365): the inner one inverts the area ratio on the
subsonic branch with the *enlarged* sonic area; the outer one moves the shock until the exit pressure is the back pressure.""")
note("N77 [B]", r"""The pressure along the nozzle for nine back pressures, each labelled with its regime.""")
figure(r"""
x = np.linspace(0.0, 1.0, 101)                               # 101 stations; the throat x = 0.4 is a grid point
fig, ax = plt.subplots(figsize=(7.4, 4.2))                   # one panel
cols = plt.cm.viridis(np.linspace(0.0, 0.9, 9))              # nine colours
for pB, col in zip((0.97, th["p_b"], 0.8, 0.7, 0.6, th["p_d"], 0.3, th["p_f"], 0.05), cols):
    sol = GD.laval_nozzle(x, GD.laval_area, pB)              # the flow for this back pressure (shock included)
    ax.plot(x, sol["p_over_p0"], color=col, lw=1.6, label=f"{pB:.3f}: {sol['regime']}")   # the shock shows as a vertical jump
ax.axvline(0.4, color=COLORS["muted"], ls=":", lw=1)         # the throat
ax.set(xlabel="x [m] (throat at 0.4, exit at 1)", ylabel="p/p₀₁ [–]",
       title="Below p_b the back pressure only moves a shock")
ax.legend(fontsize=7, title="p_B/p₀: regime", title_fontsize=7, loc="lower left")
savefig(fig, "ch15", "c10_laval_family"); plt.show()
""", see=r"""The pressure along our Laval nozzle ($A_e/A_t$ = 2) for nine back pressures, each labelled with its regime;
shocks appear as vertical jumps.""",
    read=r"""Upstream of the throat all choked curves coincide: the mass flow is fixed. Lowering $p_B$ below 0.937 starts a
shock at the throat that walks down the bell (0.8, 0.7, 0.6) and reaches the exit at 0.513. Below that the inside of the
nozzle no longer changes: the curves for 0.3, 0.094 and 0.05 are the same supersonic curve.""",
    change=r"""…the exit-to-throat area ratio were 4? The design pressure would fall to 0.0298 and the shock would need a
lower back pressure to reach the exit (0.295 instead of 0.513).""")
animation(r"""
pBs = np.array([0.99, 0.96, th["p_b"], 0.85, 0.75, 0.65, 0.58, th["p_d"], 0.4, 0.25, 0.15, th["p_f"], 0.05])   # back pressures, with the three exact dividing ones
n_frames = len(pBs)                                          # 13 frames: every one of the seven named conditions appears
x = np.linspace(0.0, 1.0, 41 if FAST else 61)                # stations (the throat x = 0.4 is a grid point)
ghost_sub = GD.nozzle_isentropic(x, GD.laval_area, branch="subsonic", A_star=1.0)   # sonic throat, subsonic after (p_b curve)
ghost_sup = GD.nozzle_isentropic(x, GD.laval_area, branch="transonic")              # supersonic after the throat (p_f curve)
fig, ax = plt.subplots(figsize=(6.8, 3.8))                   # one panel
ax.plot(x, ghost_sub["p_over_p0"], color=COLORS["muted"], ls="--", lw=1)   # ghosts: the two isentropic curves
ax.plot(x, ghost_sup["p_over_p0"], color=COLORS["muted"], ls="--", lw=1)
(line,) = ax.plot(x, ghost_sub["p_over_p0"], color=COLORS["accent"], lw=2.2)   # the current flow
(tank,) = ax.plot([1.04], [1.0], ">", color=COLORS["orange"], ms=8)            # the back pressure, just outside the exit
ax.set(xlim=(0, 1.08), ylim=(0, 1.05), xlabel="x [m] (throat 0.4, exit 1)", ylabel="p/p₀₁ [–]")
def update(i):                                               # frame i: one back pressure
    sol = GD.laval_nozzle(x, GD.laval_area, pBs[i])          # the flow, shock included
    line.set_ydata(sol["p_over_p0"]); tank.set_ydata([pBs[i]])   # move the curve and the tank marker
    ax.set_title(f"p_B/p₀ = {pBs[i]:.3f}: {sol['regime']} — ṁ/ṁ_max = {sol['mdot_over_choked']:.3f}", fontsize=9)
    return line, tank
show_animation(animate(update, frames=n_frames, fig=fig, interval=150), player="frames")   # step to stop at each regime
""", explain=r"""**Animation A3.** The back pressure falls from almost $p_0$ to 0.05 $p_0$; step frame by frame to stop at each
of the seven named conditions (the title names the regime and the mass flow as a fraction of its maximum).""")
notes3(see=r"""The pressure along the nozzle (purple) as the back pressure (orange triangle) falls; the two isentropic curves
are grey ghosts.""",
       read=r"""Above 0.937 the curve and the mass flow both change. Then ṁ/ṁ_max sticks at 1.000, a shock is born at the
throat, walks down the bell and leaves at 0.513; after that the inside follows the lower ghost and only the gap between the
exit pressure and the triangle changes — that gap is settled outside, by oblique shocks or fans.""",
       change=r"""…the walls had friction? The shock would not be a clean normal one: real nozzles show the boundary layer
separating behind it (D20, "fails when").""")
live_cell(r"""
def show(pB_over_p0=0.7, Ae_over_At=2.0):                     # draw one case
    x = np.linspace(0.0, 1.0, 101)                            # stations, throat at 0.4
    A_fn = lambda xx: GD.laval_area(xx, A_exit=Ae_over_At)    # the nozzle with this exit area
    sol = GD.laval_nozzle(x, A_fn, pB_over_p0)                # its flow
    fig, ax = plt.subplots(figsize=(6.4, 3.4))                # one panel
    ax.plot(x, sol["p_over_p0"], color=COLORS["accent"])      # pressure along the nozzle
    ax.plot([1.04], [pB_over_p0], ">", color=COLORS["orange"])   # the tank
    ax.set(ylim=(0, 1.05), xlabel="x [m]", ylabel="p/p₀₁ [–]", title=f"{sol['regime']}; shock at x = {sol['x_shock']:.3f} m")
    plt.show()
live(show, pB_over_p0=(0.05, 1.0, 0.01), Ae_over_At=(1.2, 4.0, 0.1))   # two sliders (need a running kernel)
""")
explainer("nozzle_back_pressure", "Why does a nozzle stop listening?",
          r"""you open the valve on a transport: the mass-flow needle rises and stops dead while the pressure curve keeps
reshaping — a shock is born at the throat, walks down the bell and leaves; one figure per regime cannot show that sequence.""",
          ["Play the valve from closed to wide open and watch the needle stop.",
           "Stop when the shock is half-way down the bell: read the stagnation pressure it costs.",
           "Raise the supply pressure: does the plateau move?"])
note("N76 [C]", r"""Outside the exit the jet meets the back pressure through oblique shocks when over-expanded and through
expansion fans when under-expanded; C14 and C15 (§15.9) compute the first wave angle at the lip of our nozzle.""")
figure(r"""
fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.0))         # two sketches side by side
draw_jet_sketch("overexpanded", ax=a)                        # oblique shocks at the lip
draw_jet_sketch("underexpanded", ax=b)                       # expansion fans at the lip
fig.suptitle("The jet outside the exit — a sketch, not a computation", fontsize=10)
plt.show()
""", see=r"""Two sketches of the jet just outside a supersonic exit: over-expanded (left) and under-expanded (right).""",
    read=r"""Over-expanded: the exit pressure is below the back pressure, so the jet is squeezed by oblique shocks from the
lip. Under-expanded: the exit pressure is above it, so the jet spreads through expansion fans.""",
    change=r"""…the back pressure equalled the design pressure? Neither: the jet would leave parallel and clean.""")
whatif(r"""…the duct had constant area but rough walls, or a flame inside? The same factor $1-M^2$ does the same work without
any throat. C11.""")

# =====================================================================================================================
# A.8 §15.8 Effects of Friction and Heating in Constant-Area Ducts — C11, C12
# =====================================================================================================================
nb.section("15.8", "Effects of Friction and Heating in Constant-Area Ducts", intro="**What is this section about?** A straight "
           "pipe of constant area, with wall friction and with heat added or removed. One closed formula gives the outlet "
           "Mach number; its consequences are the last two of the chapter's three surprises.")
core("C11", r"Friction and heat in one formula: $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46)",
     "A gas enters a straight pipe at a known Mach number, rubs on the wall and is heated. What is its Mach number at the outlet?")
problem(r"""
A long gas pipeline; the tail pipe of an engine; the combustion chamber of a jet. No area change at all — yet the flow speeds
up or slows down, because the wall drags on it and heat goes in. The same three duct laws apply, now with a wall force and a
heat input switched on. They collapse into a single number A that the inlet, the friction and the heat fix, and a quadratic
for the outlet Mach number.""")
idea(r"""
inlet M₁ ──[ friction f · heat q ]──► outlet M₂
(M₁, f, q) → one number A → a quadratic in M₂² → two roots: one continuous with M₁, one through a shock
""")
note("N79 [B]", r"""We start from «15.45»: C03's three laws at constant area with a retarding wall force F and a heat input Q,
made dimensionless with inlet values. > ⚠️ Common confusion (traps T5, T11): f is **not** a friction factor and q is not a
heat flux; F here is the *retarding* force (the opposite sign to F in «15.24»). The book's lumped f is not the friction length $4f_FL/D$ of the
usual "Fanno flow" tables ($f_F$ the Fanning friction factor, L the duct length, D the hydraulic diameter — conventions cell (a))
— both are shown side by side, with the bridge between them, in N86 (C12).""")
remind("C11")
P("P366", "a quadratic in M squared and its discriminant", r"""A quadratic $aX^2+bX+c=0$ (reminder: P159) has two roots when
its discriminant $b^2-4ac$ is positive, one double root when it is zero and none when it is negative. Here $X=M_2^2$: two
roots = two possible outlets; a double root = the two outlets merge (choking); none = no steady flow with this inlet.""",
  code=r"""
disc = lambda A, g=1.4: 1 - 2*A*(g + 1)       # the discriminant of (15.46)
for A in (0.144, 0.168, 1/4.8, 0.22):         # A grows with friction or heat
    print(A, disc(A))                         # 0.309, 0.194, 0.0 (choked), -0.056 (no solution)
""")
D("D21", "The outlet Mach number of a duct with friction and heat", ref="15.46", notes=r"""**N80 [B]** the duct number
«u15.8» and «15.46» (steps 1–7). The duct number A is a pure number, not an area (trap T2). Slip #14: the book says "$M_2$
calculated from (15.45)" — the three duct laws «15.45» — where it means their closed-form solution «15.46».""")
note("N83 [B]", r"""**The check that closes D21.** With f = q = 0 the continuous root is $M_2=M_1$ and the other is «15.40» —
the normal shock of C07 drops out of the duct formula. The root with $M_2>1$ from $M_1<1$ has ΔS < 0 and is excluded.
> ⚠️ Common confusion (trap T12): which sign of ± is the continuous root flips at $M_1=1$ (+ below, − above).""")
code(r"""
print(GD.duct_outlet_mach(2.0), GD.duct_outlet_mach(2.0, branch="shocked"), GD.normal_shock_mach(2.0))   # 2.0; 0.5774 = the normal shock
M2x = GD.duct_outlet_mach(0.5, branch="shocked")             # the other root from a subsonic inlet: 2.6458
print(M2x, GD.duct_entropy_change(0.5, M2x))                 # ... with (S2 - S1)/Cv = -0.326 < 0: forbidden
""")
worked("M₁ = 0.5 with friction f = 0.1, by hand", r"""
γ = 1.4, q = 0.
1. Numerator of «u15.8»: 0.25 × (1 + 0.05) = 0.2625.
2. Denominator: (1 + 0.35 − 0.1)² = 1.25² = 1.5625.
3. A = 0.2625/1.5625 = **0.168**.
4. Discriminant: 1 − 2 × 0.168 × 2.4 = 0.1936; its root is 0.44.
5. «15.46», upper sign: $M_2^2$ = (−(1 − 0.4704) + 0.44)/(0.4 − 0.6586) = (−0.0896)/(−0.2586) = 0.3465 → $M_2$ = **0.589**.
6. Lower sign: (−0.9696)/(−0.2586) = 3.75 → $M_2$ = 1.936 (supersonic from subsonic: forbidden).
7. Friction took the flow from Mach 0.50 to 0.59: it **sped up**.""")
code(r"""
A = GD.duct_parameter_A(duct["M1"], f=duct["f"])              # the duct number for M1 = 0.5, f = 0.1 -> 0.168
M2 = GD.duct_outlet_mach(duct["M1"], f=duct["f"])             # (15.46), continuous root -> 0.5887
rf = GD.duct_state_ratios(duct["M1"], M2, f=duct["f"])        # outlet / inlet ratios
print(round(A, 4), round(M2, 4), {k: round(float(rf[k]), 4) for k in ("p_ratio", "T_ratio", "rho_ratio", "p0_ratio", "T0_ratio")})
M2q = GD.duct_outlet_mach(duct["M1"], q=duct["q"])            # heat q = 0.2 instead of friction
rq = GD.duct_state_ratios(duct["M1"], M2q, q=duct["q"])
print(round(GD.duct_parameter_A(duct["M1"], q=duct["q"]), 4), round(M2q, 4), {k: round(float(rq[k]), 4) for k in ("T_ratio", "p0_ratio", "T0_ratio")})
print(GD.duct_outlet_mach(2.0, f=duct["f"]))                 # a supersonic inlet with the same friction: 2.0 -> 1.882 (it slows)
s1 = dict(rho=1.0, u=1.0, p=1.0 / (1.4 * 0.25))              # inlet in units rho1 = u1 = 1, so M1 = 0.5 means p1 = 1/(gamma M1²)
s2 = dict(rho=rf["rho_ratio"], u=1.0 / rf["rho_ratio"], p=s1["p"] * rf["p_ratio"])   # outlet from the ratios
print(GD.duct_residuals(s1, s2, f=duct["f"], q=0.0))         # (15.45): mass, momentum, energy residuals ~1e-16
""", explain=r"""
1. `GD.duct_parameter_A` is the left side of «u15.8»: 0.1680. `GD.duct_outlet_mach` picks the root that tends to $M_1$ when f
   and q go to zero: 0.5887.
2. `GD.duct_state_ratios` gives the outlet state: pressure 0.842 of the inlet's, temperature 0.982, density 0.857, stagnation
   pressure 0.897 — and $T_{02}/T_{01}$ = 1 exactly (no heat).
3. With heat q = 0.2 instead: A = 0.1715, $M_2$ = 0.6038, the gas is 16.5 % hotter, $T_0$ rose by 19 %, and 3.6 % of
   $p_0$ is lost.
4. A supersonic inlet with the same friction *slows* (2.0 → 1.882).
5. The two state dicts satisfy «15.45» to round-off.""")
scratch(r"""
g, M1, f, q = 1.4, 0.5, 0.1, 0.2                              # both friction and heat at once
def gas_law_gap(r):                                           # r = u2/u1; zero when the gas law closes the three laws
    rho = 1 / r                                               # mass: rho2/rho1
    p = 1 + g * M1**2 * (1 - r) - f                           # momentum: p2/p1
    T = 1 + q + 0.5 * (g - 1) * M1**2 * (1 - r * r)           # energy: T2/T1 (cp T1 + u1²/2 + q cp T1 = cp T2 + u2²/2)
    return p - rho * T                                        # gas law: p2/p1 = (rho2/rho1)(T2/T1)
print([round(gas_law_gap(v), 3) for v in (1.0, 1.5, 2.0, 3.0)])   # signs: -, -, +, - -> two roots; the first one is continuous with r = 1
r = brentq(gas_law_gap, 1.0, 2.0)                             # the root that is continuous with r = 1 (accelerating, subsonic)
T = 1 + q + 0.5 * (g - 1) * M1**2 * (1 - r * r)               # its temperature ratio
M2_mine = M1 * r / np.sqrt(T)                                 # M2 = u2/c2 = M1 (u2/u1)/sqrt(T2/T1)
assert np.allclose(M2_mine, GD.duct_outlet_mach(M1, f=f, q=q))   # same as (15.46)
print(M2_mine)
""", explain=r"""**From scratch, without «15.46».** One unknown, the speed ratio r: mass gives the density ratio, momentum the
pressure ratio, energy the temperature ratio, and the gas law must close them. `brentq` finds the root near r = 1 — the same
outlet Mach number as the closed form.""")
note("N81 [B] · N82 [B]", r"""Our versions of the book's two charts, with our own parameter values: $M_2$ against $M_1$ for a
range of friction (left slider) and of heat (right slider). (Slip #15: some supersonic labels on the book's heating chart lack
a decimal point; the turning points from `GD.duct_choking_limits` show which values were meant.)""")
plotly(r"""
M1 = np.linspace(0.02, 4.0, 160)                              # inlet Mach numbers
loc = GD.duct_choking_locus()                                 # f_max(M1), q_max(M1) on its own grid
def chart(f=0.0, q=0.0, key="f_max", val=0.0):                # one curve set of the chart
    c = GD.duct_curve(M1, f=f, q=q)                           # both roots, NaN where choked
    lim = loc[key] - val                                      # where the choking limit equals this f or q
    idx = np.where(np.diff(np.sign(lim)) != 0)[0]             # sign changes: the choked inlets
    return {"M₂ continuous": (M1, c["M2_continuous"]), "M₂ through a shock": (M1, c["M2_shocked"]),   # the two roots of (15.46)
            "M₂ = M₁": (M1, M1), "choked inlet (M₂ = 1)": (loc["M1"][idx], np.ones(idx.size)),   # diagonal; choking dots
            "forbidden: M₁ < 1 < M₂": ([0, 1, 1, 0, 0], [1, 1, 4, 4, 1])}   # outline of the forbidden box
modes = {"choked inlet (M₂ = 1)": "markers"}                  # dots for the choking points
fig_f = slider_figure(lambda f: chart(f=f, key="f_max", val=f), "f", np.linspace(0.0, 0.6, 13),   # friction slider (q = 0)
                      xlabel="inlet M₁ [–]", ylabel="outlet M₂ [–]", title="Friction pushes every M₂ toward 1",   # labels and message
                      xrange=[0, 4], yrange=[0, 4], modes=modes)   # fixed axes; dots for the choking points
fig_f.show()   # display
fig_q = slider_figure(lambda q: chart(q=q, key="q_max", val=q), "q", np.linspace(-0.4, 1.0, 15),   # heat slider (f = 0)
                      xlabel="inlet M₁ [–]", ylabel="outlet M₂ [–]", title="Heating does the same; cooling does the opposite",   # labels and message
                      xrange=[0, 4], yrange=[0, 4], modes=modes)   # fixed axes; dots for the choking points
fig_q.show()   # display
""", explain=r"""**Figures F7.** Sliders: the lumped friction f from 0 to 0.6 (with q = 0) and the lumped heat q from −0.4 to 1 (with f = 0).
`GD.duct_curve` returns both roots of «15.46» for every inlet; the dots are where the inlet is exactly choked
(`GD.duct_choking_locus`).""")
notes3(see=r"""The outlet Mach number against the inlet one: the continuous root (purple), the root through a shock (teal),
the diagonal "nothing happens", the choking points at $M_2=1$ and the box where the root would be supersonic from a subsonic
inlet.""",
       read=r"""At f = 0 the purple curve is the diagonal and the teal one is the normal shock of C07. Drag f up: the subsonic
part of the purple curve rises toward 1 and the supersonic part falls toward 1; between the two choking dots there is no
curve — those inlets cannot pass this much friction. With q < 0 (cooling) the curve moves *away* from 1 and never chokes.""",
       change=r"""…γ were 5/3? The choking limits change (f_max at $M_1$ = 0.5 drops from 0.228 to about 0.21): the dots move,
the picture does not.""")
whatif(r"""…you asked only for the *direction* of each change, for a little friction or a little heat? Then the quadratic
becomes a set of differentials — with the same $1-M^2$ as the area rule. C12.""")

# ---------------------------------------------------------------------------------------------------------------------
# C12 — friction and heat drive the flow to M = 1
# ---------------------------------------------------------------------------------------------------------------------
core("C12", r"Everything drifts to Mach 1: $\dfrac{du}{u}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ and $\dfrac{du}{u}=\dfrac{h_1}{h}\dfrac{dq}{1-M^2}$ (the velocity members of the friction and heat differentials (15.47), (15.48), sign-corrected)",
     "Can friction make a flow faster, and can heating make a gas colder?")
problem(r"""
Rub your hands: friction slows things and warms them. Light a stove: heat raises the temperature. In a gas flowing through a
pipe both everyday rules fail. Friction makes a subsonic flow *faster* (its pressure and density drop so quickly that the speed
must rise to carry the same mass). And in a narrow band of Mach number just below 1, adding heat makes the gas *colder*: the
heat goes into speed faster than it arrives. Both effects push the Mach number toward 1 from either side — and stop there: the
pipe chokes.""")
idea("", r"""
| | friction, M < 1 | friction, M > 1 | heating, M < 1 | heating, M > 1 |
|---|---|---|---|---|
| M | ↑ | ↓ | ↑ | ↓ |
| u | ↑ | ↓ | ↑ | ↓ |
| p | ↓ | ↑ | ↓ | ↑ |
| T | ↓ | ↑ | ↑ (then ↓ above $M^2=1/\gamma$) | ↑ |
| $p_0$ | ↓ | ↓ | ↓ | ↓ |
| S | ↑ | ↑ | ↑ | ↑ |""")
D("D22", "A little friction: the differentials, corrected", ref="15.47", notes=r"""**N84 [B]** «15.47», as corrected""")
slip(11, r"«15.47p»", r"«15.47»", why=r"""Why: mass gives $d\rho/\rho=-du/u$, energy gives $dT/T=-(\gamma-1)M^2\,du/u$, so the gas
law gives $dp/p=-[1+(\gamma-1)M^2]\,du/u$ (D22 steps 3–6). The second relation is right as printed.""")
D("D23", "A little heat: the differentials, corrected", ref="15.48", notes=r"""**N85 [B]** «15.48», as corrected""")
slip(12, r"«15.48p»", r"«15.48»", why=r"""Why: $du/u=dq/(1-M^2)$ and $u^2/h=(\gamma-1)M^2$ (D23 steps 6–8); the printed sign would
make heating *slow* a subsonic flow, against the book's own text and chart. The first relation is right as printed.""")
scratch(r"""
h = 1e-6                                                      # a little friction (or heat): df = dq = 1e-6
g = 1.4
for M in (0.5, 0.9, 2.0):                                     # subsonic, near-sonic, supersonic
    M2 = GD.duct_outlet_mach(M, f=h)                          # (15.46) with a tiny f
    fd = -(GD.duct_state_ratios(M, M2, f=h)["p_ratio"] - 1) / h   # -dp/(p1 df) by finite differences
    corr, printed = (1 + (g - 1) * M**2) / (1 - M**2), (1 - (g - 1) * M**2) / (1 - M**2)   # corrected and printed numerators of (15.47)
    c = GD.influence_coefficients(M, df=h)                    # the library's (corrected) differentials
    bad = GD.influence_coefficients(M, df=h, printed_15_47=True)   # the printed sign, kept only as a planted wrong variant
    assert np.allclose(fd, -c["dp_over_p"] / h, rtol=1e-4) and np.isclose(fd, corr, rtol=1e-4)   # finite differences = corrected form
    assert not np.allclose(fd, -bad["dp_over_p"] / h, rtol=0.05)                                  # … and ≠ the printed form
    M2q = GD.duct_outlet_mach(M, q=h)                         # (15.46) with a tiny q
    fdq = (GD.duct_state_ratios(M, M2q, q=h)["u_ratio"] - 1) / h * (g - 1) * M**2   # u du/(h1 dq) = (du/u)(u²/h1)
    cq = GD.influence_coefficients(M, dq=h); badq = GD.influence_coefficients(M, dq=h, printed_15_48=True)
    assert np.allclose(fdq, cq["u_du_over_h"] / h, rtol=1e-4) and not np.allclose(fdq, badq["u_du_over_h"] / h, rtol=0.05)
    print(f"M = {M}: friction FD {fd:+.4f} (corrected {corr:+.4f}, printed {printed:+.4f}); heat FD {fdq:+.4f} "
          f"(corrected {(g-1)*M**2/(1-M**2):+.4f}, printed {(1-g)*M**2/(1-M**2):+.4f})")
""", explain=r"""**From scratch — the cell that justifies both corrections.** We ask the *exact* duct formula «15.46» what a
tiny friction or heat does, by finite differences, and compare with both forms of each coefficient:

1. friction: 1.4667, 6.968, −0.8667 at M = 0.5, 0.9, 2 — equal to the corrected numerator; the printed one gives 1.2000, 3.558,
   +0.2000 (18 % and 49 % off, and the wrong sign);
2. heat: 0.1333, 1.7053, −0.5333 — equal to the corrected $(\gamma-1)M^2/(1-M^2)$; the printed sign is the opposite.

The library's printed-sign options exist only so that a test can show they fail.""")
P("P367", "a process line on a T–s diagram (Fanno and Rayleigh lines)", r"""Plot temperature against entropy for all the
states one duct flow can pass through as friction (or heat) accumulates: a curve with a "nose" of maximum entropy at M = 1.
Most books call the friction curve the **Fanno line** and the heating curve the **Rayleigh line** — names *not in our book*.
Their star is the sonic state of the *same* process (trap T2).""", code=r"""
M = np.array([0.5, 0.845, 1.0, 2.0])               # Mach numbers along a heating process
T_ratio = ((1 + 1.4)*M/(1 + 1.4*M**2))**2          # Rayleigh T/T*: (2.4 M / (1 + 1.4 M^2))^2
print(T_ratio)                                     # 0.79, 1.029 (the maximum, at M^2 = 1/gamma), 1.0, 0.529
""")
note("N86 [B]", r"""Heating acts like friction (drives M toward 1); cooling does the opposite and always has a solution. The
book's parameters beside the usual ones (**not in the book — ours**):

| The book | Most other sources (ours) |
|---|---|
| lumped $f=F/(p_1A)$ | friction length $4f_FL/D$ (`GD.fanno_ratios(M)["fL"]`; bridge `GD.fanno_length_to_f(M1, M2)`) |
| $q=Q/h_1$ | $T_0/T_0^*$ (`GD.rayleigh_ratios(M)["T0"]`); $T_{02}/T_{01}=1+q/(1+\tfrac12(\gamma-1)M_1^2)$ |

**Reading the right-hand column.** $f_F=\tau_w/(\tfrac12\rho u^2)$ is the **Fanning friction factor**: the wall shear stress
$\tau_w$ [Pa] over the local dynamic pressure. L is the duct length [m] and D its **hydraulic diameter** [m] (4 × area ÷
wetted perimeter; the diameter for a round pipe — not the drag D of C16, nor the L of lift). **The bridge, slice by slice:**
on a slice of length dx the wall force is $\tau_w$ times the wetted area, $dF=\tau_w\,(4A/D)\,dx$; divide by $p_1A$ and use
$\tfrac12\rho u^2=\tfrac12\gamma pM^2$ (primer P359):

$$ df=\dfrac{4f_F\,dx}{D}\cdot\tfrac12\gamma M^2\,\dfrac{p}{p_1}. $$

So the two parameters grow together but not in proportion — the factor $\tfrac12\gamma M^2p/p_1$ changes along the duct;
`GD.fanno_length_to_f` evaluates both for the same duct. In the Rayleigh column, $T_0^*$ is the stagnation temperature at
which *that same heating process* would reach Mach 1 (trap T2: not the isentropic star of C05).""")
code(r"""
b_ = GD.fanno_length_to_f(0.5, GD.duct_outlet_mach(0.5, f=0.1))   # Mach 0.5 -> 0.5887 by friction: the book's f and the usual 4 f_F L/D
print({k: round(float(v), 4) for k, v in b_.items()})        # f = 0.100, friction length = 0.532
print(GD.fanno_ratios(0.5)["fL"], GD.rayleigh_ratios(0.5)["T0"])   # 4 f L*/D = 1.069 to choke from 0.5; T0/T0* = 0.691
print(rq["T0_ratio"], 1 + duct["q"] / (1 + 0.2 * 0.25))      # T02/T01 = 1 + q/(1 + ½(gamma-1) M1²) = 1.1905 both ways
M_a = GD.duct_outlet_mach(0.5, f=0.1)                         # first duct: f = 0.1
M_b = GD.duct_outlet_mach(M_a, f=0.1)                         # a second, identical duct after it
f_eq = 0.1 + 0.1 * GD.duct_state_ratios(0.5, M_a, f=0.1)["p_ratio"]   # each f is normalised with its own inlet pressure
print(M_a, M_b, GD.duct_outlet_mach(0.5, f=0.2), f_eq, GD.duct_outlet_mach(0.5, f=f_eq))   # 0.589, 0.728; 0.775; 0.184 -> 0.728
""", explain=r"""
1. The bridge: our duct (Mach 0.5 → 0.589) is f = 0.100 in the book's parameter and $4f_FL/D$ = 0.532 in the usual one.
2. The usual tables' ratios at M = 0.5: the friction length to choke is 1.069; the Rayleigh $T_0/T_0^*$ is 0.691.
3. The heat parameter q = 0.2 raises $T_0$ by the factor 1.1905, as the formula in the table says.
4. > ⚠️ (trap T11) **Two ducts in series do not add their f**, because each f is normalised with its own inlet pressure:
   0.1 then 0.1 takes Mach 0.5 to 0.589 to 0.728; a single f = 0.2 gives 0.775. The equivalent single value is
   0.1 + 0.1 × $p_2/p_1$ = 0.184, which indeed gives 0.728.""")
note("N87 [B]", r"""From the first relation of «15.48»: $dT/dq<0$ when $1/\gamma<M^2<1$ — for γ = 1.4, 0.845 < M < 1. Where
does the heat go? Into kinetic energy, faster than it arrives.

> ⚠️ **Common confusion (trap T18):** heating *cools* the gas only for $1/\gamma<M^2<1$. That is a statement about one band of
> Mach number, not about heating in general — below it and above Mach 1 heat warms the gas — and in every band the entropy
> still rises.""")
code(r"""
print(GD.rayleigh_temperature_sign(np.array([0.5, 0.9, 2.0])))   # sign of dT/dq: +1, -1, +1
M2h = GD.duct_outlet_mach(0.9, q=0.005)                       # a little heat at M = 0.9
print(GD.duct_state_ratios(0.9, M2h, q=0.005)["T_ratio"])     # 0.9947: heat in, temperature down 0.5 %
""")
D("D24", "The choking limits of friction and of heat", ref="o_fmax", notes=r"""**N88 [B]** the choking limits «o_fmax» and
«o_qmax» (not in the book — ours)""")
code(r"""
print(GD.duct_choking_limits(0.5))                            # f_max = 0.2275, q_max = 0.4688 for an M = 0.5 inlet
print(GD.duct_choking_limits(2.0))                            # at M1 = 2: the same q_max, by a coincidence of gamma = 1.4
""")
note("N89 [C]", r"""The funnel picture: pour faster than the neck allows and the level rises upstream — once choked, more
friction or heat forces the *inlet* state to change. The same idea as area choking in C09 (§15.7).""")
note("N90 [C]", r"""Dense gases near their critical point can show expansion shocks and "antithroats" (the literature the book
cites); for a perfect gas only compression shocks exist (N61 in C08).""")
worked("a little friction at Mach 0.5 and at Mach 2", r"""
«15.47» with df = 0.01 at the station where $p=p_1$.
1. M = 0.5: du/u = 0.01/(1 − 0.25) = **+1.33 %** (faster); −dp/p = (1 + 0.4 × 0.25)/0.75 × 0.01 = **1.47 %** (pressure falls).
2. M = 2: du/u = 0.01/(1 − 4) = **−0.33 %** (slower); −dp/p = (1 + 1.6)/(−3) × 0.01 = −0.87 % (pressure *rises*).
3. A little heat dq = 0.01 at M = 0.9 («15.48»): dT/T = (1 − 1.4 × 0.81)/(1 − 0.81) × 0.01 = (−0.134/0.19) × 0.01 = **−0.71 %**
   (colder).""")
code(r"""
for M, kw in ((0.5, dict(df=0.01)), (2.0, dict(df=0.01)), (0.9, dict(dq=0.01)), (0.5, dict(dq=0.01))):
    c = GD.influence_coefficients(M, **kw)                    # (15.47) / (15.48), corrected, at the station where p = p1
    print(M, kw, {k: round(c[k], 5) for k in ("du_over_u", "dp_over_p", "dT_over_T", "dp0_over_p0", "ds_over_cp")})
""", explain=r"""Each line is one case of the worked example. In every case with df > 0 or dq > 0 the stagnation pressure falls
(`dp0_over_p0` < 0) and the entropy rises (`ds_over_cp` > 0) — the second law, whichever way the speed goes.""")
animation(r"""
n = 5 if FAST else 7                                          # frames per process (the frames player costs ~3 s per frame here)
fm, qm = GD.duct_choking_limits(0.5)["f_max"], GD.duct_choking_limits(0.5)["q_max"]   # the choking limits for M1 = 0.5
steps = [("f", v) for v in np.linspace(0.0, 0.999 * fm, n)] + [("q", v) for v in np.linspace(0.0, 0.999 * qm, n)]
Mline = np.linspace(0.25, 2.5, 300)                           # Mach numbers along the two process lines
fan, ray = GD.fanno_line_Ts(Mline), GD.rayleigh_line_Ts(Mline)   # T/T* against (s - s*)/cp (ours)
fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.8), gridspec_kw=dict(width_ratios=[1.4, 1]))
(lf,) = a.plot(fan["ds_over_cp"], fan["T_over_Tstar"], color=COLORS["blue"], label="Fanno (friction)")
(lr,) = a.plot(ray["ds_over_cp"], ray["T_over_Tstar"], color=COLORS["rose"], label="Rayleigh (heat)")
(pt,) = a.plot([], [], "o", color=COLORS["ink"], ms=7)        # the outlet state
a.set(xlabel="(s − s*)/C_p [–]", ylabel="T/T* [–]", xlim=(-0.9, 0.05), ylim=(0.3, 1.35),
      title="Fanno / Rayleigh line — usual names, not in the book"); a.legend(fontsize=8, loc="lower left")
bars = b.bar(["M₂", "T₂/T₁", "p₀₂/p₀₁"], [0.5, 1.0, 1.0], color=[COLORS["blue"], COLORS["rose"], COLORS["amber"]])
b.axhline(1.0, color=COLORS["muted"], ls="--", lw=1); b.set(ylim=(0, 1.4), title="outlet gauges")
def update(i):                                                # frame i
    kind, v = steps[i]                                        # friction or heat, and how much
    M2 = GD.duct_outlet_mach(0.5, **{kind: v})                # (15.46)
    r_ = GD.duct_state_ratios(0.5, M2, **{kind: v})           # outlet / inlet
    line = GD.fanno_line_Ts(M2) if kind == "f" else GD.rayleigh_line_Ts(M2)   # where the outlet sits on its line
    pt.set_data([line["ds_over_cp"]], [line["T_over_Tstar"]])
    lf.set_alpha(1.0 if kind == "f" else 0.25); lr.set_alpha(1.0 if kind == "q" else 0.25)   # highlight the active line
    for bar, val in zip(bars, (M2, r_["T_ratio"], r_["p0_ratio"])):
        bar.set_height(val)                                   # the three gauges
    b.set_xlabel(f"{'f' if kind == 'f' else 'q'} = {v:.3f}  (M₂ = {M2:.3f})")
    return pt, lf, lr, *bars
show_animation(animate(update, frames=len(steps), fig=fig, interval=150), player="frames")   # step: stop at the nose
""", explain=r"""**Animation A5.** For an inlet at Mach 0.5: first friction grows from 0 to its choking value (the outlet slides
along the blue line), then heat does (the rose line). Right: gauges of the outlet Mach number, temperature ratio and
stagnation-pressure ratio.""")
notes3(see=r"""Two process lines on a temperature–entropy diagram, both scaled with their own sonic state (the "nose" at
(0, 1)), and the outlet state moving along one of them; on the right three gauges.""",
       read=r"""Both runs end at the nose: M = 1, maximum entropy — the duct is choked. In the friction run the temperature
falls a little and $p_0$ drops. In the heating run the temperature climbs, peaks at M = 0.845 (where $M^2=1/\gamma$) and then
*falls* while heat is still going in — the paradox of N87.""",
       change=r"""…the inlet were supersonic? The point would start on the lower branches of the lines and climb to the same
nose from below: friction and heat slow a supersonic flow toward Mach 1.""")
explainer("duct_friction_heating", "Can friction speed a flow up?",
          r"""you turn up f or q and the outlet point slides along the M₂(M₁) map toward the line M₂ = 1; the badge flips to
"choked" exactly when the discriminant bar reaches zero; flipping the inlet to supersonic reverses every arrow — four regimes
and a limit that a single figure cannot hold at once.""",
          ["Raise f from 0 at $M_1$ = 0.5: which way does $M_2$ go?",
           "Set $M_1$ = 0.9 and add a little heat: watch the temperature bar.",
           "Push f up to the choking limit, then beyond."])
whatif(r"""…we left the duct for open air? "No news upstream" becomes "no news outside a cone". C13.""")

# =====================================================================================================================
# A.9 §15.9 Pressure Waves in Planar Compressible Flow — C13, C14, C15
# =====================================================================================================================
nb.section("15.9", "Pressure Waves in Planar Compressible Flow", intro="**What is this section about?** The same ideas in two "
           "dimensions: where the news of a supersonic body can reach (a cone), how a supersonic stream is turned into itself "
           "(an oblique shock) and away from itself (an expansion fan).")
core("C13", r"The Mach cone: $\sin\mu=1/M$ (15.49)", "Why do you hear a supersonic aircraft only after it has passed overhead?")
problem(r"""
A duck paddling slowly sends ripples ahead of itself. A speedboat outruns its own ripples and drags a V behind it. An aircraft
does the same with sound. Each instant it emits a sphere of "news" that grows at the speed of sound from where the aircraft
*was*. If the aircraft is faster than sound the spheres can never get ahead of it; they pile up along a cone. Outside the
cone: silence. On it: the boom.""")
idea(r"""
                         ● source (now)
                   . ´ ╱ |
             . ´     ╱   |  c·Δt (radius of an old front)
       . ´  μ      ╱     |
  ────────────────●──────┘
                 centre of that front (where the source was)
            |◄──────── U·Δt ────────►|
""", r"""$\sin\mu=(c\,\Delta t)/(U\,\Delta t)=1/M$.""")
remind("C13")
note("N91 [B]", r"""Wave fronts sent out at equal intervals by a source at rest, a subsonic source and a supersonic source.""")
figure(r"""
fronts = [AC.source_wavefronts(U, 1.0, 4.0, 8, 0.5) for U in (0.0, 0.5, 2.0)]   # c = 1; fronts emitted every 0.5 s up to t = 4 s
fig = draw_wavefront_panels(fronts)                          # one panel per source speed
plt.show()
for w in fronts:                                             # the numbers behind the pictures
    print(f"M = {w['M']}: Mach angle {np.degrees(w['mu']):.1f}°" if w["M"] >= 1 else f"M = {w['M']}: fronts reach {w['ahead'][0]:.2f} ahead of the source")
""", see=r"""Fronts emitted every half second by a point source at rest (left), at M = 0.5 (middle) and at M = 2 (right); each
circle is centred where the source was when it emitted it.""",
    read=r"""At rest the circles are concentric. At M = 0.5 they bunch up ahead but still enclose the source: the news gets
ahead of it. At M = 2 every circle lies behind the source and they all touch two straight lines through it — the cone, at
30° to the flight path.""",
    change=r"""…the source moved at exactly M = 1? Every circle would touch the source itself: the fronts pile up into a plane
wall through it (μ = 90°).""")
animation(r"""
from matplotlib.patches import Circle                         # circle artists
n_frames = 10 if FAST else 16                                 # frames (each two-panel frame costs 1–3 s to draw)
T, dt, Nmax = 5.0, 0.5, 10                                    # total time, emission interval, most fronts shown
times = np.linspace(0.01, T, n_frames)                        # [s] (c = 1 m/s)
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4))            # M = 0.5 and M = 2 side by side
panels = []
for ax, M in zip(axes, (0.5, 2.0)):                           # set up each panel once
    circles = [ax.add_patch(Circle((0, 0), 0.0, fill=False, color=COLORS["blue"], lw=1)) for _ in range(Nmax)]
    (src,) = ax.plot([], [], "o", color=COLORS["rose"], ms=6)  # the source
    (up,) = ax.plot([], [], color=COLORS["orange"], lw=1.5)    # the cone (upper generator)
    (dn,) = ax.plot([], [], color=COLORS["orange"], lw=1.5)    # the cone (lower generator)
    ax.set(xlim=(-5.5, 11), ylim=(-5.5, 5.5), aspect="equal", title=f"source at M = {M}")
    ax.set_xticks([]); ax.set_yticks([])
    panels.append((M, circles, src, up, dn))
def update(i):                                                # frame i
    t = times[i]
    out = []
    for M, circles, src, up, dn in panels:
        n = min(int(t / dt), Nmax)                            # fronts emitted so far (the source started at x = 0, t = 0)
        w = AC.source_wavefronts(M, 1.0, t, max(n, 1), dt)    # centres and radii of the fronts
        for k, cir in enumerate(circles):
            cir.set_visible(k < n)                            # only fronts already emitted
            if k < n:
                cir.center = (w["x_centres"][k], 0.0); cir.set_radius(w["radii"][k])
        xs = w["x_source"]; src.set_data([xs], [0.0])         # the source now
        if M > 1:                                             # the envelope: two lines at the Mach angle
            L = M * t; mu = np.arcsin(1 / M)
            up.set_data([xs, xs - L * np.cos(mu)], [0, L * np.sin(mu)]); dn.set_data([xs, xs - L * np.cos(mu)], [0, -L * np.sin(mu)])
        out += circles + [src, up, dn]
    return out   # all artists of this frame (circles, sources, Mach lines)
show_animation(animate(update, frames=n_frames, fig=fig, interval=80), player="video")   # smooth motion: video
""", explain=r"""**Animation A4.** Two sources start at x = 0 and send out a front every half second (c = 1). Left M = 0.5, right
M = 2; on the right the envelope (orange) appears by itself.""")
notes3(see=r"""Growing circles of news behind a moving source; on the right, two orange lines through the source that every
circle touches.""",
       read=r"""Nobody draws the cone: it is just where the circles crowd together. A listener below the right-hand source hears
nothing until an orange line sweeps over them — and then hears the accumulated news all at once: the boom.""",
       change=r"""…the source slowed from 2 to 1.2? The cone would open from 30° to 56°: the boom arrives sooner after the
aircraft passes overhead.""")
P("P368", "the envelope of a family of circles", r"""A curve that touches every member of a family is its envelope. For circles
whose centres lie on a line and whose radii grow in proportion to the distance from a point P, the envelope is the pair of
tangents from P. A tangent meets the radius at its touching point at a right angle — that right angle makes the triangle.""",
  code=r"""
d, r = 2.0, 1.0                               # centre 2 behind the apex, radius 1
mu = np.degrees(np.arcsin(r/d))               # angle of the tangent from the apex
print(mu)                                     # 30 degrees: sin(mu) = r/d
""")
D("D25", "The Mach angle", ref="15.49", after=r"""(μ is the Mach angle here; a viscosity in §15.6 — trap T5.)""")
figure(r"""
fig = draw_mach_triangle(2.0)                                 # the right triangle behind sin(mu) = 1/M, for M = 2
plt.show()
""", see=r"""The right triangle of D25 for M = 2: the source, the centre of an old front, and the point where the cone touches
it.""",
    read=r"""Hypotenuse UΔt, opposite side cΔt, right angle at the touching point: $\sin\mu=1/M$ = 0.5, μ = 30°.""",
    change=r"""…M were 3? The opposite side would be a third of the hypotenuse: μ = 19.5°, a slender cone.""")
note("N92 [B]", r"""Inside the cone: the **zone of action**; outside: the **zone of silence**. The cone's generators are **Mach
lines**. (Ours:) an aircraft at height h passes overhead; the boom arrives «o_delay» later.""")
code(r"""
print(AC.mach_cone_delay(1.0e4, 1.5, 295.1))                  # 10 km up, M = 1.5, c = 295 m/s: the boom arrives 25.3 s after overhead
xg = np.array([-1.0, 1.0, 3.0]); yg = np.array([1.0, 1.0, 1.0])   # three points, x measured DOWNSTREAM of the source (behind the aircraft)
print(AC.in_zone_of_action(xg, yg, 2.0))                      # ahead: False; 1 behind: False (sqrt(3) > 1); 3 behind: True
""")
P("P369", "hyperbolic against elliptic equations", r"""Elliptic: every point influences every other (subsonic flow, Laplace's
equation; reminder P201). Hyperbolic: influence travels only along special lines, the characteristics (waves; reminder P174,
P227). Steady supersonic flow is hyperbolic with the downstream distance playing the part of time, and its characteristics
are the Mach lines.""", code=r"""
for M in (0.5, 2.0):                          # the sign of 1 - M^2 decides the type
    print(M, 'elliptic' if 1 - M*M > 0 else 'hyperbolic')
""")
worked("three Mach numbers", r"""
«15.49»: M = 2 → sin μ = 0.5 → μ = **30°**. M = 1.2 → sin μ = 0.833 → μ = 56.4° (a blunt cone). M = 3 → μ = 19.5° (a slender
one). M = 1 → μ = 90°: a plane front. M < 1: sin μ > 1, no cone.""")
code(r"""
print(np.degrees(GD.mach_angle(np.array([1.0, 1.2, 1.5, 2.0, 3.0]))))   # 90, 56.44, 41.81, 30.00, 19.47 degrees
print(GD.mach_angle(0.8))                                     # nan: no cone below Mach 1
print(GD.mach_from_angle(np.deg2rad(30.0)))                   # the inverse: 2.0
""", explain=r"""`GD.mach_angle` is «15.49» in radians (converted to degrees for printing); it returns NaN below Mach 1, and
`GD.mach_from_angle` inverts it.""")
scratch(r"""
w = AC.source_wavefronts(2.0, 1.0, 1.0, 200, 0.005)           # 200 fronts of a source at M = 2 (c = 1)
xs = np.linspace(w["x_source"] - 0.9, w["x_source"] - 0.1, 50)   # stations behind the source
y_top = [max(np.sqrt(np.clip(r**2 - (x - xc)**2, 0, None)) for xc, r in zip(w["x_centres"], w["radii"])) for x in xs]   # the envelope: highest circle at each x
slope = np.polyfit(xs, y_top, 1)[0]                           # a straight-line fit (np.polyfit, degree 1) of that edge
assert np.allclose(np.arctan(-slope), GD.mach_angle(2.0), rtol=1e-3)   # its angle is the Mach angle
print(np.degrees(np.arctan(-slope)))                          # 30.0 degrees
""", explain=r"""**From scratch.** Draw 200 fronts, take the highest circle above each station behind the source (the envelope),
fit a straight line to it with `np.polyfit(x, y, 1)` (slope and intercept of the best line): its angle is 30.0°, the Mach
angle.""")
explainer("mach_cone", "Why do you hear it only after it has passed?",
          r"""you play time and each circle grows from where the source was; above Mach 1 the circles cross and their envelope
appears by itself; you drag an observer and read "not yet … now … heard" — the delay is something you experience.""",
          ["Play at M = 0.5, then M = 1, then M = 2.", "Click a circle: where was the source when it left?",
           "Drag the observer below the flight path and read the delay."])
whatif(r"""…the disturbance were not infinitely weak — a wedge that really turns the flow? The Mach line steepens into an
oblique shock. C14.""")

# ---------------------------------------------------------------------------------------------------------------------
# C14 — the oblique shock
# ---------------------------------------------------------------------------------------------------------------------
core("C14", r"Turning a supersonic stream into itself: $\tan\delta=2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53)",
     "A supersonic stream meets a wedge. At what angle does the shock stand — and how sharp can the wedge be before no such shock exists?")
problem(r"""
The nose of a supersonic aircraft, the lip of an engine intake, a bullet: each forces the air to turn. The air cannot be warned
(no news upstream), so it turns all at once, through a thin sloping front. Here is the trick that makes it easy: walk along the
front at just the right speed and the flow you see comes straight at the front — a *normal* shock. So an oblique shock is a
normal shock on the component of the velocity normal to it; the component along it passes through untouched.""")
idea(r"""
V₁ = (u₁ ⟂ shock, v ∥ shock)  ──►  V₂ = (u₂ < u₁, v unchanged)
""", r"""The normal part is shocked: replace $M_1$ by $M_1\sin\sigma$ everywhere in C07.""")
remind("C14")
note("N93 [B]", r"""The geometry: shock angle σ from the upstream direction, deflection δ;
$(u_1,v)=V_1(\sin\sigma,\cos\sigma)$, $(u_2,v)=V_2(\sin(\sigma-\delta),\cos(\sigma-\delta))$. > ⚠️ Common confusion (trap
T13): the book's σ and δ are the β and θ of most other sources; both pairs are in the legends.""")
figure(r"""
w = GD.oblique_shock(wedge["M1"], delta=wedge["delta"])       # the weak shock for our 10° wedge at Mach 2
fig = draw_shock_triangles(wedge["M1"], w["sigma"], wedge["delta"])   # velocity triangles ahead of and behind the shock
plt.show()
print(f"σ (= β elsewhere) = {np.degrees(w['sigma']):.2f}°, δ (= θ elsewhere) = {np.degrees(wedge['delta']):.1f}°")
""", see=r"""An oblique shock at σ = 39.3° turning a Mach 2 stream through δ = 10°, with the velocity split into the part across
the shock (u) and the part along it (v), ahead and behind.""",
    read=r"""The along-shock part v is the same on both sides; only the across part shrinks (u₂ < u₁). That alone turns the
velocity toward the shock — by δ.""",
    change=r"""…σ were 90°? There would be no along part: a normal shock, no turn.""")
P("P370", "resolving a velocity along and across a tilted line", r"""Any velocity is the sum of a part along a line and a part
across it: for a line at angle σ to the velocity, across = V sin σ, along = V cos σ. An observer sliding along the line at the
"along" speed sees only the "across" part — and sliding at constant speed changes no pressure, density or temperature
(reminder: P96).""", code=r"""
V, sigma = 680.6, np.deg2rad(39.31)            # Mach 2 at sea level; a 39.31 degree shock
print(V*np.sin(sigma), V*np.cos(sigma))         # across 431.2 m/s (Mach 1.267), along 526.6 m/s
""")
D("D26", "An oblique shock is a normal shock on the normal component", ref="15.50", notes=r"""**N94 [B]** the normal Mach
numbers «u15.9n» (step 4) · **N95 [B]** «15.50» (step 6) · **N96 [B]** «15.51» (step 7) · **N97 [B]** «15.52» (step 8; ⚠️ trap
T13: $M_2$ needs the division by $\sin(\sigma-\delta)$)""")
code(r"""
tm = ch15.tangential_momentum_sympy()                         # sympy: mass and tangential momentum force v1 = v2
print(tm["v2"], "| residual:", tm["residual"])                # v2 equals v1; residual 0
""")
P("P371", "the tangent subtraction formula", r"""$\tan(\sigma-\delta)=\dfrac{\tan\sigma-\tan\delta}{1+\tan\sigma\tan\delta}$,
and the double angle $\cos2\sigma=1-2\sin^2\sigma$ (reminder: P341, ch14).""", code=r"""
s, d = np.deg2rad(39.31), np.deg2rad(10.0)
print(np.tan(s - d), (np.tan(s) - np.tan(d))/(1 + np.tan(s)*np.tan(d)))   # 0.5613 both
print(np.cos(2*s), 1 - 2*np.sin(s)**2)                                    # 0.1972 both
""")
D("D27", "Deflection against shock angle", ref="15.53")
worked("a 10° wedge at Mach 2", r"""
γ = 1.4. Trial: σ = 39.31° (sin σ = 0.6336, sin² σ = 0.4014, cos 2σ = 0.1972, cot σ = 1.2213).
1. «15.53»: tan δ = 2 × 1.2213 × (4 × 0.4014 − 1)/(4 × (1.4 + 0.1972) + 2) = 2.4426 × 0.6056/8.389 = 0.1763 → δ = **10.0°** ✓.
2. Normal Mach number: $M_{n1}$ = 2 × 0.6336 = **1.267**.
3. «15.50»: $p_2/p_1$ = 1 + 1.1667 × (1.6056 − 1) = **1.707**.
4. «15.40» on 1.267: $M_{n2}$ = 0.803; divide by sin(39.31° − 10°) = 0.4896: $M_2$ = **1.64** — still supersonic.
5. The other solution of «15.53» for the same 10°: σ = 83.7°, $p_2/p_1$ = 4.44, $M_2$ = 0.60 — nearly a normal shock.""")
code(r"""
w = GD.oblique_shock(2.0, delta=wedge["delta"])                         # weak root (the default branch)
s = GD.oblique_shock(2.0, delta=wedge["delta"], branch="strong")        # strong root
for name, d in (("weak", w), ("strong", s)):
    print(name, f"σ = {np.degrees(d['sigma']):.2f}°", {k: round(float(d[k]), 4) for k in ("Mn1", "p_ratio", "rho_ratio", "T_ratio", "M2", "p0_ratio")})
print({k: round(float(v), 4) for k, v in GD.oblique_shock_components(2.0, w["sigma"]).items()})   # the split velocities, normal Mach numbers
mx = GD.max_deflection(2.0)                                   # the largest turn at Mach 2
print(f"δ_max = {np.degrees(mx['delta_max']):.2f}° at σ = {np.degrees(mx['sigma_max']):.2f}°")   # 22.97° at 64.67°
try:                                                          # a 25° wedge at Mach 2 …
    GD.oblique_shock(2.0, delta=np.deg2rad(25.0))
except ValueError as err:                                     # … has no attached shock
    print("ValueError:", err)
print(GD.detachment_mach(np.deg2rad(10.0)))                   # the slowest stream that keeps a shock attached to a 10° wedge: 1.421
""", explain=r"""
1. The weak shock: σ = 39.31°, $M_{n1}$ = 1.267, pressure × 1.707, density × 1.458, temperature × 1.170, $M_2$ = 1.641, and only
   1.5 % of the stagnation pressure lost. The strong one: σ = 83.70°, pressure × 4.44, $M_2$ = 0.604.
2. `GD.oblique_shock_components` lists the pieces of D26: the normal and tangential parts and both normal Mach numbers.
3. The largest turn at Mach 2 is 22.97°, at σ = 64.67°. Ask for 25° and the function raises "detached".
4. Turned round: a 10° wedge keeps its shock attached only above Mach 1.421.""")
scratch(r"""
g, M1, delta = 1.4, 2.0, np.deg2rad(10.0)                     # gas, Mach number, wedge half-angle
tan_delta = lambda s: 2/np.tan(s)*(M1**2*np.sin(s)**2 - 1)/(M1**2*(g + np.cos(2*s)) + 2)   # (15.53)
mu = np.arcsin(1 / M1)                                        # the Mach angle: tan(delta) = 0 there
grid = np.linspace(mu + 1e-6, np.pi / 2 - 1e-6, 20001)       # shock angles between the two zeros
s_max = grid[np.argmax(tan_delta(grid))]                      # grid search for the largest deflection
f = lambda s: tan_delta(s) - np.tan(delta)                    # zero where the shock gives 10°
s_weak = brentq(f, mu + 1e-9, s_max)                          # left of the maximum: weak
s_strong = brentq(f, s_max, np.pi / 2 - 1e-9)                 # right of it: strong
assert np.allclose([s_weak, s_strong], [w["sigma"], s["sigma"]], atol=1e-8)   # same roots as the library
print(np.degrees([s_weak, s_max, s_strong]))                  # 39.31, ≈ 64.7, 83.70 degrees
""", explain=r"""**From scratch.** «15.53» in one line; a grid search finds the angle of the largest turn, and one bracketed root on
each side of it gives the weak and the strong shock — the two roots of P360 again.""")
note("N98 [B]", r"""(The closed form of $\delta_{\max}$ is ours; result given, not derived.) For δ below a maximum there are two
shocks: **weak** (smaller σ, usually supersonic behind) and **strong** (larger σ, subsonic behind); δ = 0 at σ = μ (a Mach wave)
and at σ = 90° (a normal shock). > ⚠️ (trap T13) "weak" and "strong" name the two roots, not their size: the "weak" shock for
20° at Mach 2 nearly triples the pressure ($p_2/p_1$ = 2.84).""")
code(r"""
from scipy.optimize import minimize_scalar                    # minimum of a function of one variable in a bracket
res = minimize_scalar(lambda s_: -GD.deflection_angle(2.0, s_), bounds=(np.arcsin(0.5), np.pi / 2), method="bounded")   # maximise delta
print(np.degrees(-res.fun), np.degrees(GD.max_deflection(2.0)["delta_max"]))   # 22.97 both: a numerical cross-check of the closed form
print(GD.oblique_shock(2.0, delta=np.deg2rad(20.0))["p_ratio"])   # the "weak" shock for 20°: p2/p1 = 2.84
""")
plotly(r"""
def curve(M1):                                                # delta against sigma for one upstream Mach number
    c = GD.oblique_shock_curve(M1)                            # (15.53) on 400 shock angles, with its special points
    return {"δ(σ), (15.53)": (np.degrees(c["sigma"]), np.degrees(c["delta"])),   # the curve, in degrees
            "largest turn": ([np.degrees(c["sigma_max"])], [np.degrees(c["delta_max"])]),   # its maximum
            "sonic point (M₂ = 1)": ([np.degrees(c["sigma_sonic"])], [np.degrees(c["delta_sonic"])])}   # where the flow behind is sonic
fig = slider_figure(curve, "M₁", np.linspace(1.2, 5.0, 20), xlabel="shock angle σ (= β elsewhere) [°]",   # slider over M1
                    ylabel="deflection δ (= θ elsewhere) [°]", title="Two shock angles for every turn below the largest",   # both conventions in the labels
                    xrange=[0, 90], yrange=[0, 48], modes={"largest turn": "markers", "sonic point (M₂ = 1)": "markers"})   # fixed axes; dots
fig.show()   # display
""", explain=r"""**Figure F8.** Slider: the upstream Mach number, from 1.2 to 5. Each curve is «15.53»; the dots mark the largest turn and
the point where the flow behind becomes exactly sonic.""")
notes3(see=r"""The deflection against the shock angle for the chosen upstream Mach number, with the largest-turn point and the
sonic point.""",
       read=r"""Draw a horizontal line at your wedge angle: it cuts the curve twice (weak on the left, strong on the right) or
not at all. The largest turn is 22.97° at Mach 2, 34.07° at Mach 3, and creeps toward 45.6° for very fast streams. Between the
sonic point and the maximum the weak shock already leaves a subsonic flow behind.""",
       change=r"""…γ were 5/3? The largest turn at Mach 2 would drop (to about 19°): monatomic gases detach sooner.""")
note("N99 [B]", r"""A wedge takes the weak solution; a blunt nose always has a detached shock.""")
figure(r"""
fig, (a, b, c) = plt.subplots(1, 3, figsize=(11.0, 3.4))      # three wedges
for ax, deg in ((a, 10.0), (b, 22.0)):                        # an easy wedge and a nearly detached one
    wf = ch15.wedge_flow(2.0, np.deg2rad(deg))                # the numbers: sigma, M2, p_ratio, attached, subsonic_behind …
    draw_wedge_streamlines(2.0, np.deg2rad(deg), wf["sigma"], ax=ax)   # the helper only draws them
    ax.set_title(f"{deg:.0f}° wedge at M = 2: σ = {np.degrees(wf['sigma']):.1f}°, M₂ = {wf['M2']:.2f}", fontsize=9)
draw_bow_shock_sketch(ax=c)                                   # the blunt case
c.set_title("detached bow shock — a sketch, not a computation", fontsize=9)
plt.show()
""", see=r"""Left and middle: streamlines through the attached weak shock on a 10° and a 22° wedge at Mach 2. Right: a sketch of
the detached bow shock ahead of a body too blunt for an attached shock.""",
    read=r"""At 10° the shock leans back at 39° and the flow behind is still supersonic ($M_2$ = 1.64). At 22° — just below
the 23° limit — the shock stands up (σ ≈ 58°) and the flow behind is only just supersonic ($M_2$ = 1.08). Beyond the limit the
shock leaves the nose and curves: normal on the axis, strong near it, weak far out.""",
    change=r"""…the Mach number were raised to 3 for the 22° wedge? The shock would lean back again (σ ≈ 40°) with a supersonic
flow behind it.""")
live_cell(r"""
def show(M1=2.0, delta_deg=10.0):                             # one wedge
    d = np.deg2rad(delta_deg)                                 # half-angle [rad]
    try:   # the exact method raises when the lower shock detaches
        wk = GD.oblique_shock(M1, delta=d); st = GD.oblique_shock(M1, delta=d, branch="strong")   # both roots
        print(f"weak: σ = {np.degrees(wk['sigma']):.2f}°, p2/p1 = {wk['p_ratio']:.3f}, M2 = {wk['M2']:.3f} | "
              f"strong: σ = {np.degrees(st['sigma']):.2f}°, p2/p1 = {st['p_ratio']:.3f}, M2 = {st['M2']:.3f}")
        draw_wedge_streamlines(M1, d, wk["sigma"])            # the picture with the weak root
    except ValueError:                                        # beyond the largest turn
        print("detached: no attached shock for this wedge")
        draw_bow_shock_sketch()
    plt.show()
live(show, M1=(1.2, 5.0, 0.05), delta_deg=(1.0, 40.0, 0.5))   # two sliders (need a running kernel)
""")
note("N100 [C]", r"""At a corner the speed is finite in supersonic flow (zero or infinite in subsonic flow — Ch. 6 corner flows,
C09), and a supersonic field can be built marching downstream: the hyperbolic property of N92.""")
P("P372", "the Mach triangle", r"""A right triangle with hypotenuse M and opposite side 1 has the Mach angle μ: $\sin\mu=1/M$, so
$\cos\mu=\sqrt{M^2-1}/M$ and $\cot\mu=\sqrt{M^2-1}$. Every $\sqrt{M^2-1}$ of §15.9–§15.10 is this triangle's third side.""",
  code=r"""
M = 2.0
mu = np.arcsin(1/M)                             # Mach angle: 30 degrees
print(1/np.tan(mu), np.sqrt(M*M - 1))           # cot(mu) = sqrt(M^2 - 1) = 1.732
""")
D("D28", "The weak-wave rule", ref="15.55", notes=r"""**N101 [B]** «15.54» (steps 1–3) · **N102 [B]** «15.55» (steps 4–8)""")
note("N102 [B], continued", r"""> ⚠️ Common confusion (trap T16): how small is "small"? The cell below compares the linear rule «15.55» with
the exact weak shock at $M_1$ = 2. **The rule reads LOW**: about 1 % at 0.5°, 4 % at 2°, 10 % at 5°, 20 % at 10° and 39 % at 20°
(equivalently, at 5° the exact value is 12 % above the rule). And the rule fails as $M\to1$: $\sqrt{M^2-1}\to0$.""")
code(r"""
rows = []                                                     # linear rule against the exact weak shock at M1 = 2
for deg in (0.5, 2.0, 5.0, 10.0, 20.0):
    e = GD.weak_wave_error(2.0, np.deg2rad(deg))              # dict: linear, exact, ratio (exact/linear), rel_error (linear/exact - 1)
    rows.append((deg, e["linear"], e["exact"], -100 * e["rel_error"]))   # rel_error = linear/exact - 1 < 0: the rule is low
print(pd.DataFrame(rows, columns=["δ [°]", "linear (15.55)", "exact shock", "rule low by [%]"]).round(4).to_string(index=False))
print(ch15.weak_wave_sympy()["residual_15_55"])               # sympy: the small-delta limit of (15.54) is (15.55) -> 0
""")
explainer("oblique_shock_wedge", "One wedge, how many shocks?",
          r"""you open the wedge angle and the two solution dots on the δ–σ curve run toward each other, merge at the maximum and
vanish — and the picture switches to a detached bow shock; the merger of two roots is a motion, best seen moving.""",
          ["Open the wedge from 2° to 23° at Mach 2 and watch the two dots.",
           "Click a streamline: which velocity component did the shock change?",
           "Switch to the strong root: how far is it from a normal shock?"])
whatif(r"""…the wall turned *away* from the stream? There is no expansion shock (C08). The stream fans out instead. C15.""")

# ---------------------------------------------------------------------------------------------------------------------
# C15 — the Prandtl–Meyer function
# ---------------------------------------------------------------------------------------------------------------------
core("C15", r"Turning away, one Mach wave at a time: $\nu(M)=\sqrt{\dfrac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\dfrac{\gamma-1}{\gamma+1}\left(M^2-1\right)}-\tan^{-1}\sqrt{M^2-1}$ (15.56)",
     "How does a supersonic stream get round a corner that bends away from it — and how fast is it afterwards?")
problem(r"""
Behind the shoulder of a bullet, at the lip of a rocket nozzle in thin air, over the top of a tilted supersonic wing, the wall
falls away from the stream. A slow flow would simply follow, or separate. A supersonic one turns through a *fan* of infinitely
many weak waves spreading from the corner. Each wave turns the stream a little and speeds it up a little; nothing is lost. Add
the little turns up and you get one function, ν(M), and every such turn becomes a subtraction.""")
idea(r"""
  M₁ ──────────►  ╲  ╲   ╲    ╲     (rays of the fan, all from the corner)
 ════════════════●╲   ╲    ╲    ╲
                   ╲═══════════════►  M₂ > M₁ (wall turned away)
  first ray at μ₁ to the old wall · last ray at μ₂ to the new wall · turn = ν(M₂) − ν(M₁)
""")
remind("C15")
note("N105 [B]", r"""**The fan:** all Mach lines start at the corner; the flow is isentropic; along a streamline the pressure falls
smoothly, along the wall it steps at the corner.""")
figure(r"""
fan = ch15.expansion_fan(2.0, np.deg2rad(10.0))               # Mach 2 turned 10° away (the book's delta = -10°)
fig = draw_fan(fan)                                           # rays and one streamline coloured by Mach number
plt.show()
print(f"M: {fan['M'][0]:.3f} -> {fan['M2']:.3f}; first ray at {np.degrees(fan['mu'][0]):.1f}°, last at {np.degrees(fan['mu'][-1]):.1f}° to its wall; p2/p1 = {fan['p2_over_p1']:.3f}")
""", see=r"""A centred expansion fan at a 10° convex corner in a Mach 2 stream: seven rays from the corner and one streamline
that bends smoothly through them.""",
    read=r"""The first ray stands at μ₁ = 30° to the upstream wall, the last at μ₂ = 24.8° to the downstream wall; between them
the Mach number climbs from 2 to 2.385 and the pressure drops to 55 % — smoothly along the streamline, but in one step along
the wall at the corner.""",
    change=r"""…the turn were 20°? The fan would open wider and the stream would leave at Mach 2.83 with 27 % of its pressure.""")
D("D29", "The Prandtl–Meyer differential relation", ref="u15.9d", notes=r"""**N106 [B]** «u15.9d»""")
code(r"""
pm = ch15.prandtl_meyer_sympy()                               # sympy: the integrand, the closed form and their consistency
print(pm["integrand"], "| d(closed form)/dM - integrand:", pm["derivative_residual"], "| nu(1) =", pm["value_at_1"])
""")
P("P373", "arctangent integrals", r"""$\displaystyle\int\dfrac{dt}{1+a^2t^2}=\dfrac1a\tan^{-1}(at)$: the antiderivative of "one
over one plus a square" is an arctangent, with the scaling a divided out. Reached here after a substitution (reminder P106) and
a split into partial fractions (reminder P219).""", code=r"""
import sympy as sp
t, a = sp.symbols('t a', positive=True)
print(sp.integrate(1/(1 + a**2*t**2), t))            # atan(a*t)/a
print(sp.diff(sp.atan(a*t)/a, t))                    # 1/(a**2*t**2 + 1): differentiate back
""")
D("D30", "The Prandtl–Meyer function in closed form", ref="15.56", after=r"""(ν is the Prandtl–Meyer function here; a kinematic
viscosity in §15.6 — trap T5.)""")
note("N107 [B]", r"""**How it is used:** «u15.9t», so $\nu(M_2)=\delta_1-\delta_2+\nu(M_1)$. > ⚠️ Common confusion (trap T14):
the book counts δ positive when the wall turns **into** the flow; our `turn=` (and most tables) count the **expansion**
positive: turn = −δ. Both numbers are in every legend.""")
worked("Mach 2 turned 10° away", r"""
γ = 1.4; $\sqrt{(\gamma+1)/(\gamma-1)}$ = √6 = 2.449.
1. «15.56» at M = 2: $\sqrt{M^2-1}$ = 1.732; first term 2.449 × atan(1.732/2.449) = 2.449 × atan(0.7071) = 2.449 × 35.26° = 86.38°;
   second atan(1.732) = 60°; **ν(2) = 26.38°**.
2. Turn +10° (the book's δ = −10°): ν₂ = 36.38°.
3. Invert (root-finder): **M₂ = 2.385**.
4. Isentropic, so $p_0$ is unchanged: $p_2/p_1$ = [(1 + 0.2 × 4)/(1 + 0.2 × 2.385²)]^3.5 = (1.8/2.1375)^3.5 = **0.548**.
5. Fan: first ray at $\mu_1$ = 30° to the old wall, last at $\mu_2$ = 24.8° to the new wall.""")
code(r"""
print(np.degrees(GD.prandtl_meyer(np.array([1.0, 1.5, 2.0, 3.0]))))   # (15.56): 0, 11.905, 26.380, 49.757 degrees
print(np.degrees(GD.prandtl_meyer_max()))                     # the largest possible turn from Mach 1: 130.45°
t = GD.prandtl_meyer_turn(2.0, np.deg2rad(10.0))              # turn = +10° (expansion) = the book's delta = -10°
print({k: round(float(v), 4) for k, v in t.items() if k in ("M2", "p_ratio", "T_ratio")}, "mu2 =", round(np.degrees(t["mu2"]), 2),
      "°, book's delta =", round(np.degrees(t["delta_book"]), 1), "°")
print(GD.mach_from_prandtl_meyer(np.deg2rad(36.38)))          # inverse: 2.385
""", explain=r"""
1. `GD.prandtl_meyer` is «15.56» in radians: 26.38° at Mach 2.
2. `GD.prandtl_meyer_max()` is the limit as M → ∞: (√6 − 1) × 90° = 130.45° — a stream starting at Mach 1 can turn at most that
   far, when all its heat has become speed (N38).
3. `GD.prandtl_meyer_turn(2.0, 10°)` gives $M_2$ = 2.3849, $p_2/p_1$ = 0.5480, $T_2/T_1$ = 0.8421, $\mu_2$ = 24.79°, and reports
   the book's δ = −10° beside our turn = +10°.""")
scratch(r"""
nu_quad = quad(lambda M: np.sqrt(M * M - 1) / (1 + 0.2 * M * M) / M, 1.0, 2.0)[0]   # (15.56, integral form) by quad
assert np.allclose(nu_quad, GD.prandtl_meyer(2.0))            # = the closed form
M, n = 2.0, 100                                               # march 100 small turns of 0.1°
dturn = np.deg2rad(10.0) / n                                  # each turn [rad] (expansion: d(delta) = -dturn)
for _ in range(n):                                            # (unnumbered, §15.9): dM = -d(delta) M (1 + 0.2 M²)/sqrt(M² - 1)
    M += dturn * M * (1 + 0.2 * M * M) / np.sqrt(M * M - 1)
assert np.allclose(M, t["M2"], rtol=2e-3)                     # the staircase lands on the exact answer
print(np.degrees(nu_quad), M, [round(ch15.fan_staircase(2.0, np.deg2rad(10.0), N)["M"][-1], 4) for N in (1, 3, 7, 21, 100)])
""", explain=r"""**From scratch.** Two routes: `quad` integrates «15.56i» directly (26.380°), and a loop takes 100 small Mach-wave
steps of «u15.9d» — the integral as a limit of sums. `ch15.fan_staircase` does the same with N waves: 2.363, 2.377, 2.381, 2.384,
2.3846 for N = 1, 3, 7, 21, 100, closing on the exact 2.3849.""")
note("N104 [B]", r"""Mach lines from a curved wall. Along each Mach line the flow direction and Mach number are constant (a simple
wave — the steady twin of N22's steepening).""")
figure(r"""
x = np.linspace(-0.5, 1.5, 41)                                # wall stations [m]
ramp = lambda s: np.deg2rad(15.0) * np.clip(s, 0.0, 1.0)      # wall angle grows from 0 to 15° between x = 0 and 1
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))            # away from the stream / into it
for ax, sign, name in ((axes[0], -1, "wall curving away: lines fan out"), (axes[1], +1, "wall curving into the stream: lines run together")):
    cw = ch15.curved_wall_mach_lines(2.0, lambda s: sign * ramp(s), x)   # theta > 0 into the flow (the book's sign)
    ax.fill_between(cw["x_wall"], cw["y_wall"], cw["y_wall"].min() - 0.3, color=COLORS["grid"])   # the wall
    ax.plot(cw["x_wall"], cw["y_wall"], color=COLORS["ink"])
    for xw, yw, ang in zip(cw["x_wall"][::2], cw["y_wall"][::2], cw["line_angle"][::2]):   # every second Mach line
        ax.plot([xw, xw + 2.5 * np.cos(ang)], [yw, yw + 2.5 * np.sin(ang)], color=COLORS["blue"], lw=0.8)
    if cw["first_crossing"] is not None:                      # compression: where the first two lines meet
        ax.plot(*cw["first_crossing"], "*", color=COLORS["rose"], ms=12, label="first crossing: a shock starts")
        ax.legend(fontsize=8, loc="upper left")
    ax.set(xlim=(-0.5, 3.0), ylim=(-0.6, 1.8), aspect="equal", title=name); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("Mach 2 over a curved wall (15° turn)", fontsize=10); plt.show()
""", see=r"""Mach lines leaving a wall that curves through 15°: away from a Mach 2 stream (left) and into it (right).""",
    read=r"""Turning away, each Mach line leans further downstream than the one before (higher Mach number, smaller μ, and a
wall pointing down): they spread — an expansion fan spread along the wall. Turning into the stream they lean the other way and
cross; at the star the flow would have two Mach numbers at once, so a shock must begin there.""",
    change=r"""…the bend were sharpened into a corner? On the left all lines would start at one point: the centred fan of
N105. On the right the crossing would move onto the wall: an oblique shock from the corner (C14).""")
note("N103 [B]", r"""Weak expansions obey the weak-wave rule too. They share value and slope at δ = 0 because the entropy change is
third order: «15.44b» — so the rule holds for weak expansions too, with the sign of δ reversed.""")
figure(r"""
deg = np.linspace(-15.0, 15.0, 121)                           # turns from -15° (expansion) to +15° (compression), the book's sign
d = np.deg2rad(deg)   # the same turns in radians
shock = [GD.oblique_shock(2.0, delta=v)["p_ratio"] - 1 if v > 0 else np.nan for v in d]          # exact weak shock (delta > 0)
fan_ = [GD.prandtl_meyer_turn(2.0, -v)["p_ratio"] - 1 if v <= 0 else np.nan for v in d]          # exact fan (delta < 0)
isen = [GD.prandtl_meyer_turn(2.0, -v)["p_ratio"] - 1 if v > 0 else np.nan for v in d]           # isentropic compression (for comparison)
fig, ax = plt.subplots(figsize=(6.8, 3.8))   # one panel
ax.plot(deg, GD.weak_wave_pressure(2.0, d), color=COLORS["muted"], ls="--", label="linear rule (15.55)")   # the straight line
ax.plot(deg, shock, color=COLORS["rose"], lw=2, label="exact oblique shock")   # compression by a shock
ax.plot(deg, isen, color=COLORS["teal"], ls=":", lw=1.5, label="isentropic compression")   # compression without a shock
ax.plot(deg, fan_, color=COLORS["blue"], lw=2, label="exact Prandtl–Meyer fan")   # expansion
ax.axhline(0, color=COLORS["ink"], lw=0.6); ax.axvline(0, color=COLORS["ink"], lw=0.6)   # axes through the origin
ax.set(xlabel="turn δ [°] (+ into the flow, the book's sign)", ylabel="(p₂ − p₁)/p₁ [–]",   # labels with units and sign
       title="Shock and fan share value and slope at δ = 0")   # the message
ax.legend(fontsize=8); savefig(fig, "ch15", "c15_weak_wave"); plt.show()   # legend, save, display
print([round(float(v), 4) for v in (GD.oblique_shock(2.0, delta=np.deg2rad(10))["p_ratio"] - 1,   # the four quoted numbers:
       GD.prandtl_meyer_turn(2.0, np.deg2rad(-10))["p_ratio"] - 1, GD.weak_wave_pressure(2.0, np.deg2rad(10)),   # shock, isentropic, line …
       GD.prandtl_meyer_turn(2.0, np.deg2rad(10))["p_ratio"] - 1)])   # +10°: shock 0.7066, isentropic 0.7052, line 0.5643; -10°: fan -0.4520
""", see=r"""The relative pressure change against the turn at Mach 2: the exact shock (rose, δ > 0), the exact fan (blue, δ < 0),
an isentropic compression for comparison (teal dotted) and the straight line «15.55» (grey dashed).""",
    read=r"""All curves leave the origin together: same value, same slope. At +10° the shock (0.707) and the isentropic
compression (0.705) are almost the same — the entropy difference is third order — and both are above the line (0.564). At −10°
the fan (−0.452) lies above the line (−0.564): the rule underestimates compressions and overestimates expansions.""",
    change=r"""…$M_1$ were 1.2? The line's slope $\gamma M^2/\sqrt{M^2-1}$ would be 3.0 instead of 3.2, but the curves would
bend away much sooner — and the shock detaches above 3.9°.""")
plotly(r"""
M = np.geomspace(1.0, 10.0, 300)                              # Mach numbers (log-spaced)
fig = slider_figure(lambda g: {"ν(M), (15.56)": (M, np.degrees(GD.prandtl_meyer(M, g))),   # nu(M) in degrees
                               "limit as M → ∞": (M, np.degrees(GD.prandtl_meyer_max(g)) * np.ones_like(M))},   # its limit, as a flat line
                    "γ", [1.2, 1.3, 1.4, 5 / 3], xlabel="Mach number M [–] (log axis)", ylabel="ν [°]",   # four gases
                    title="The Prandtl–Meyer function and the largest possible turn")   # the message
fig.update_xaxes(type="log", range=[0, 1])                    # log axis from 1 to 10
fig.show()   # display
""", explain=r"""**Figure F9b.** Slider: γ = 1.2, 1.3, 1.4 or 5/3, in «15.56».""")
notes3(see=r"""ν(M) in degrees against the Mach number on a logarithmic axis, with its limit as M → ∞ dashed.""",
       read=r"""Read a turn as a horizontal distance between two ν values: from Mach 2 (26.4°) a 10° turn reaches 36.4°, i.e.
Mach 2.385. The limit is 130.5° for air but 90° for a monatomic gas.""",
       change=r"""…γ → 1? The limit $(\sqrt{(\gamma+1)/(\gamma-1)}-1)\cdot90°$ runs away: a gas with many internal degrees of freedom
can turn very far.""")
explainer("prandtl_meyer_fan", "How does a supersonic stream turn a corner?",
          r"""you build the turn from 1, 3, 7, 21, then infinitely many Mach waves and watch the staircase of Mach number settle
onto ν(M): the integral as a limit of sums, made with your own hands.""",
          ["Set N = 1, 3, 7, 21, ∞ and compare $M_2$.", "Reverse the turn: where do the Mach lines cross?",
           "Start at $M_1$ = 1: how far can the stream turn?"])
whatif(r"""…a thin plate were tilted in the stream — a fan on one side, a weak shock on the other? Two pressures, and their
difference is lift. C16.""")

# =====================================================================================================================
# A.10 §15.10 Thin Airfoil Theory in Supersonic Flow — C16
# =====================================================================================================================
nb.section("15.10", "Thin Airfoil Theory in Supersonic Flow", intro="**What is this section about?** Lift and drag of a thin "
           "plate in a supersonic stream from the weak-wave rule alone: two pressure levels, a lift — and a drag that needs "
           "no viscosity.")
core("C16", r"Lift and wave drag of a flat plate: $C_L\cong\dfrac{4\alpha}{\sqrt{M_\infty^2-1}},\ \ C_D\cong\dfrac{4\alpha^2}{\sqrt{M_\infty^2-1}}$ (15.59)",
     "How does a thin wing make lift above Mach 1 — and why does that lift always cost drag, even without friction?")
recap_md("d'Alembert's paradox and the Kutta condition", "Ch. 6 C07; Ch. 14 C02", r"""**d'Alembert's paradox:** in steady,
inviscid, subsonic flow a body feels no drag at all. **The Kutta condition:** potential theory allows a wing any circulation;
the real one is the circulation that makes the flow leave the sharp trailing edge smoothly, and it gives the lift
«14.12» (slope 2π per radian). Neither carries over above Mach 1 — the plate below has drag without friction, and needs no Kutta condition.""")
problem(r"""
Chapter 14's wing made lift by circulation, needed a sharp trailing edge to choose it, and in two dimensions had no drag at
all. A supersonic wing is simpler and harsher. Each surface feels only its own slope: the top, turned away from the stream,
sits at one lower pressure; the bottom, turned into it, at one higher pressure. The difference pushes at right angles to the
plate. Tilt the plate by α and that push leans back by α: every newton of lift costs α newtons of drag — carried away by the
waves, not by friction.""")
idea(r"""
            p₂ < p∞   (fan above the leading edge)
      ════════════════  plate at α  ════════════════
            p₃ > p∞   (shock below)
  force = (p₃ − p₂) b, normal to the plate  →  L = force × cos α,  D = force × sin α
""")
note("N108 [B]", r"""«15.57» is «15.55» written with free-stream values (Ackeret's linear theory). As a pressure coefficient
(ours): $C_p=2\delta/\sqrt{M_\infty^2-1}$.""")
recap_md("Pressure coefficient", "Ch. 4 C15; P359 in C04", r"""$C_p=(p-p_\infty)/(\tfrac12\rho_\infty U_\infty^2)$, and in a gas
$\tfrac12\rho U^2=\tfrac12\gamma pM^2$ (primer P359), so $C_p=\dfrac{2}{\gamma M_\infty^2}\dfrac{p-p_\infty}{p_\infty}$.""")
code(r"""
print(ch15.ackeret_pressure(plate["alpha"], plate["M_inf"]))  # (15.57) for a 2° slope at Mach 2: (p - p_inf)/p_inf = 0.1129 (the angle comes first)
print(ch15.ackeret_cp(plate["alpha"], plate["M_inf"]))        # as a pressure coefficient: 2 delta/sqrt(M² - 1) = 0.0403
""")
D("D31", "The flat plate in supersonic flow", ref="15.59", notes=r"""**N109 [B]** the two surface pressures from «15.57» (steps
1–3) · **N110 [B]** «15.58», the normal force resolved into lift and drag (steps 4–5; ⚠️ trap T15: b is the **chord** here)""")
slip(13, r"""that the plate pressures are "found from (15.47)" — the friction differentials «15.47», which have nothing to do with a plate""",
     r"from «15.57»")
figure(r"""
pat = ch15.flat_plate_wave_pattern(plate["alpha"], plate["M_inf"])   # the four waves and the states between them
fig = draw_plate_waves(pat)                                   # the plate, its four waves, the two pressure levels; inset: force N -> lift L + drag D
plt.show()
print({k: round(float(pat[k]), 4) for k in ("p2", "p3", "M2", "M3")})   # upper 0.892, lower 1.118 (over p_inf); M 2.073 above, 1.928 below
""", see=r"""A flat plate (black, chord b) at 2° in a Mach 2 stream: an expansion fan (thin blue rays) above and a shock
(thick purple line) below at the leading edge, a shock above and a fan below at the trailing edge, and the two pressure levels
(blue label above the plate, purple below). **Inset, lower left:** the force on the plate — the normal force
$N=(p_3-p_2)\,b$ (purple arrow, at right angles to the plate) split into lift L (teal, vertical) and wave drag D (rose,
horizontal, downstream). In the inset the tilt is drawn **five times too large** (10° for α = 2°); at the true 2° the drag
arrow would be only 3.5 % of the lift arrow and invisible.""",
    read=r"""The upper surface sits at 0.892 $p_\infty$ and the lower at 1.118 $p_\infty$, each uniform along the chord. The
force is normal to the plate, so it leans back by the tilt: in the inset the purple arrow, the teal leg and the rose leg form a
right triangle. The angle α sits at the bottom, between the teal (L) and purple (N) arrows where both start on the plate; the
right angle is at the tip of the teal arrow, where the rose leg begins. So D = N sin α ≈ αL. The trailing-edge waves bring both streams back to
$p_\infty$.""",
    change=r"""…α were 0? No waves of either kind (a plate of zero thickness at zero incidence does not disturb the stream): no
lift, no drag.""")
worked("2° at Mach 2", r"""
α = 2° = 0.0349 rad, $M_\infty$ = 2, $\sqrt{M_\infty^2-1}$ = 1.732.
1. «15.57»: each surface differs from $p_\infty$ by 1.4 × 4 × 0.0349/1.732 = **11.3 %** (up below, down above).
2. «15.59»: $C_L$ = 4 × 0.0349/1.732 = **0.0806**.
3. $C_D$ = 4 × 0.0349²/1.732 = **0.00281**.
4. $C_D/C_L$ = α = 0.0349: L/D = 28.6.
5. Slope: 4/1.732 = 2.31 per radian — against Ch. 14's 2π = 6.28 for the same plate in slow flow.""")
code(r"""
for deg in (2.0, 10.0):                                       # our plate, then a large incidence
    a_ = np.deg2rad(deg)
    lin = ch15.flat_plate_supersonic(a_, plate["M_inf"], method="linear")            # (15.57)-(15.59)
    ex = ch15.flat_plate_supersonic(a_, plate["M_inf"], method="shock_expansion")    # exact: oblique shock below, fan above
    print(f"{deg:4.1f}°: linear CL {lin['CL']:.5f} CD {lin['CD']:.6f} p_up {lin['p_upper']:.4f} p_low {lin['p_lower']:.4f} | "
          f"exact CL {ex['CL']:.5f} CD {ex['CD']:.6f} p_up {ex['p_upper']:.4f} p_low {ex['p_lower']:.4f}")
print(ch15.ackeret_flat_plate(plate["alpha"], plate["M_inf"]))   # CL, CD, L/D and the lift slope 4/sqrt(M² - 1) = 2.309
L_per_span = 0.5 * 1.4 * plate["p_inf"] * plate["M_inf"]**2 * plate["b"] * ch15.ackeret_flat_plate(plate["alpha"], plate["M_inf"])["CL"]
print(f"lift per metre of span for a 1 m chord at sea-level pressure: {L_per_span/1e3:.1f} kN/m")   # (15.58)
""", explain=r"""
1. At 2° the linear theory gives $C_L$ = 0.08061 and $C_D$ = 0.002814; the exact shock–expansion answer is 0.08064 and
   0.002816 — the same to 0.05 %.
2. At 10° each surface's pressure *change* is 20–25 % off (linear ±0.564 against the exact +0.707 below and −0.452 above), yet
   the lift is only 1.1 % off. **The relative error of a single surface is first order in α** (it grows in proportion to α:
   the shock below is under-predicted, the fan above over-predicted), but those leading errors cancel in the difference
   $p_3-p_2$, so **the relative error of $C_L$ and $C_D$ is second order in α** — which is why 2° agrees to 0.05 % and
   10° still to about 1 %.
3. `ch15.ackeret_flat_plate` returns the lift slope 2.309 per radian.
4. «15.58» with our chord and pressure: about 23 kN of lift per metre of span.""")
scratch(r"""
g, M, alpha = 1.4, plate["M_inf"], plate["alpha"]             # gas, Mach number, incidence [rad]
dp = g * M**2 * alpha / np.sqrt(M**2 - 1)                     # (15.57): each surface's relative pressure change
p3, p2 = 1 + dp, 1 - dp                                       # lower and upper pressures over p_inf
qinf = 0.5 * g * M**2                                         # dynamic pressure over p_inf (primer P359)
CN_mine = (p3 - p2) / qinf                                    # normal-force coefficient
CL_mine, CD_mine = CN_mine, CN_mine * alpha                   # cos(alpha) ≈ 1, sin(alpha) ≈ alpha: D31 step 5
lin = ch15.flat_plate_supersonic(alpha, M, method="linear")
assert np.allclose([CL_mine, CD_mine], [lin["CL"], lin["CD"]])   # = the linear method
ack = ch15.ackeret_flat_plate(alpha, M)
assert np.allclose([CL_mine, CD_mine], [ack["CL"], ack["CD"]])   # = (15.59)
ex = ch15.flat_plate_supersonic(alpha, M, method="shock_expansion")
pl = GD.oblique_shock(M, delta=alpha)["p_ratio"]              # exact lower surface: weak oblique shock through alpha
pu = GD.prandtl_meyer_turn(M, alpha)["p_ratio"]               # exact upper surface: Prandtl–Meyer turn of alpha
CN = (pl - pu) / qinf                                         # exact normal-force coefficient
assert np.allclose([CN * np.cos(alpha), CN * np.sin(alpha)], [ex["CL"], ex["CD"]])   # = the shock–expansion method
print(CL_mine, CD_mine, CN * np.cos(alpha), CN * np.sin(alpha))
""", explain=r"""**From scratch.** Six lines of D31 reproduce the linear method and «15.59»; four more — an oblique shock below, a
Prandtl–Meyer turn above — reproduce the exact shock–expansion answer.""")
note("N111 [B]", r"""Wave drag exists in inviscid supersonic flow: d'Alembert's paradox (Ch. 6 C07) does not apply, because waves
carry energy away — like a ship's bow waves. No Kutta condition is needed (Ch. 14 C02), yet the flow has circulation. The theory
fails as $M_\infty\to1$.""")
figure(r"""
M = np.linspace(0.0, 4.0, 401)                                # free-stream Mach numbers
fig, ax = plt.subplots(figsize=(6.8, 3.8))   # one panel
ax.plot(M, ch15.lift_slope_vs_mach(M), color=COLORS["accent"], lw=2, label="2π (14.12) below · 4/√(M²−1) (15.59) above")   # NaN in the gap
Msub = M[M < 0.8]                                              # subsonic part
ax.plot(Msub, 2 * np.pi / np.sqrt(1 - Msub**2), color=COLORS["muted"], ls="--", label="Prandtl–Glauert 2π/√(1−M²) — not in the book")
ax.axvspan(0.8, 1.2, color=COLORS["grid"], alpha=0.6)         # the transonic gap: neither theory
ax.text(1.0, 9.5, "transonic:\nneither", ha="center", fontsize=8, color=COLORS["muted"])
ax.plot(2.0, 4 / np.sqrt(3), "o", color=COLORS["orange"])     # our plate: 2.309
ax.set(xlim=(0, 4), ylim=(0, 12), xlabel="free-stream Mach number M∞ [–]", ylabel="lift slope dC_L/dα [1/rad]",
       title="Supersonic wings lift less per radian")
ax.legend(fontsize=8); savefig(fig, "ch15", "c16_lift_slope"); plt.show()
print(ch15.lift_slope_vs_mach(np.array([0.5, 2.0])), np.sqrt(1 + (4 / (2 * np.pi))**2))   # 2π, 2.309; the slopes would be equal at 1.185
""", see=r"""The lift slope per radian against the Mach number: 2π of «14.12» on the subsonic side (with the Prandtl–Glauert
correction, not in the book, dashed), $4/\sqrt{M_\infty^2-1}$ of «15.59» on the supersonic side, and a gap around Mach 1.""",
    read=r"""At Mach 2 the slope is 2.31 per radian — 37 % of the low-speed 2π. The two theories would give the same slope at
Mach 1.185, but both fail near Mach 1 (the shaded band), where the flow is a mixture of the two kinds.""",
    change=r"""…$M_\infty$ were 1.1? The linear theory would promise a slope of 8.7 per radian — above 2π — because
$\sqrt{M^2-1}\to0$; that is the theory failing, not the wing improving.""")
note("N112 [B]", r"""(Ours; result given, not derived.) A thin section of any shape: «o_thin», where θ is the local surface slope
and the bars are chord averages of the squared camber-line and thickness slopes. **Lift depends on the incidence only;
thickness and camber buy nothing but drag.**""")
code(r"""
yu = lambda x: 0.025 * (1 - np.abs(2 * x - 1))                # a 5 % diamond: half-thickness 0.025 c at mid-chord
yl = lambda x: -yu(x)                                         # symmetric lower surface
dia = ch15.ackeret_airfoil(yu, yl, 0.0, 2.0)                  # linear theory for any thin shape, at zero incidence
print(round(dia["CL"], 6), round(dia["CD"], 6), round(dia["CD_thickness"], 6), 4 * 0.05**2 / np.sqrt(3))   # 0, 0.00577 (about twice the 2° plate's 0.00281: the diamond's slopes are 2.86°), no lift
""")
note("N113 [C]", r"""Supersonic sections are thin with sharp leading edges; a rounded nose brings the detached shock of N99 (C14)
and its drag.""")
plotly(r"""
alphas = np.deg2rad(np.linspace(0.0, 12.0, 25))               # incidences 0 … 12°
def exact_or_nan(a_, M, key):                                 # the exact coefficient, NaN once the lower shock detaches
    try:   # the exact method raises when the lower shock detaches
        return ch15.flat_plate_supersonic(a_, M, method="shock_expansion")[key]   # exact coefficient
    except ValueError:                                        # "detached: …"
        return np.nan   # nothing to plot
def plate_curves(M):                                          # four curves for one Mach number
    lin = [ch15.flat_plate_supersonic(a_, M, method="linear") for a_ in alphas]   # (15.59)
    ex = [(exact_or_nan(a_, M, "CL"), exact_or_nan(a_, M, "CD")) if a_ < GD.max_deflection(M)["delta_max"] else (np.nan, np.nan)   # exact (C_L, C_D)
          for a_ in alphas]                                   # exact pair, skipped beyond the largest turn
    return {"C_L linear (15.59)": (np.degrees(alphas), [d["CL"] for d in lin]),   # linear lift
            "C_L exact (shock–expansion)": (np.degrees(alphas), [e[0] for e in ex]),   # exact lift
            "10 C_D linear (15.59)": (np.degrees(alphas), [10 * d["CD"] for d in lin]),   # linear drag × 10
            "10 C_D exact": (np.degrees(alphas), [10 * e[1] for e in ex])}   # exact drag × 10
Ms = ([1.2, 1.5, 2.0, 2.5, 3.0, 4.0] if FAST else                                # Mach numbers on the slider; 2.0 (our plate) is one
      [1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0])
fig = slider_figure(plate_curves, "M∞", Ms, xlabel="incidence α [°]", ylabel="coefficient [–]",   # slider over the Mach number
                    title="Linear theory against the exact plate (exact curves stop where the shock detaches)", yrange=[0, 2.5],   # the message
                    active=Ms.index(2.0))                     # open on Mach 2, the case read below
fig.show()   # display
""", explain=r"""**Figure F9a.** Slider: the free-stream Mach number, from 1.2 to 4 (12 positions; it opens on Mach 2). Linear theory against the exact shock–expansion
answer; the drag is multiplied by 10 to share the axis.""")
notes3(see=r"""Lift and ten times the drag coefficient against incidence, linear and exact, for the chosen Mach number.""",
       read=r"""At Mach 2 the pairs lie on top of each other up to about 6° and stay within a few per cent to 12°. Slide toward
Mach 1.2: the curves steepen, the exact ones bend away early, and they stop where the lower shock detaches (about 4° at Mach
1.2).""",
       change=r"""…you halved the incidence? Lift halves, drag quarters: $C_D=\alpha C_L$ — fly a supersonic wing at small α.""")
explainer("supersonic_flat_plate", "How can there be drag without friction?",
          r"""you tilt the plate: the fan and the shock swap sides at α = 0 and the force arrow, always normal to the plate, leans
back by exactly α — the drag bar grows as α², the lift bar as α; seeing the two bars grow at different rates is the point.""",
          ["Tilt from 0 to 10°: compare how the two bars grow.", "Slide $M_\infty$ toward 1: which bars run away?",
           "Switch to a diamond at zero incidence: drag without lift."])
whatif(r"""…the "gas" were blood and the "duct" an elastic artery? The pulse is a pressure wave of exactly the kind C02 solved,
«15.15», with a wave speed set by the wall — and a collapsible tube chokes like the nozzle of C09. Chapter 16.""")
nb.pointer(X(r"""**S01** The chapter's exercises are not reproduced. The results they hold are taught above: the initial-value form
(N21), the particle velocity (D04), the small-Mach expansions (N32), the area ratio (D10), the weak shock (D17), the maximum speed
(N38), the strong shock (R17), the two sets of differentials (D22, D23), the choked mass flow (D11), the choking limits (D24), the
thin section (N112) and the oblique-shock components (D26). (Slips #16 and #17 concern the wording of one exercise and a loss
quoted in another; nothing of either is published here.) **S02** Sources: the reference list at the end of the book's chapter;
our public sources (NACA Report 1135 and the pages cited in the analysis) are in `reference/ch15/SOURCES.md`."""))

nb.summary(
    clicked=[X(s) for s in (
        r"**C01** Sound is what the laws of mass and momentum become when the disturbance is weak: one wave equation «15.12» with the speeds U ± c.",
        r"**C02** Any plane sound is one shape going right plus one going left, «15.15»; the air moves with a compression and against a rarefaction, «15.16».",
        r"**C03** A duct flow is three budgets between two stations, and the energy budget «15.22» holds even across friction and shocks.",
        r"**C04** «15.28»: $T_0$ is what an adiabatic flow keeps; $p_0$ is what it loses when entropy is made.",
        r"**C05** The area a flow needs, «15.31», is least at Mach 1, so each area belongs to two Mach numbers and a throat passes only so much.",
        r"**C06** Above Mach 1 the gas thins faster than it speeds up, so a supersonic nozzle widens: «15.34».",
        r"**C07** Across a shock the three laws leave one solution besides 'nothing', fixed by $M_1$ alone: «15.39».",
        r"**C08** Only the supersonic-to-subsonic jump makes entropy, «15.43»; its price is stagnation pressure.",
        r"**C09** When the exit is sonic the nozzle cannot hear the back pressure and the mass flow stops rising at «o_mmax».",
        r"**C10** With a sonic throat the back pressure only chooses where a shock stands.",
        r"**C11** Friction and heat in a straight duct reduce to one number A and a quadratic for $M_2^2$: «15.46».",
        r"**C12** The same $1-M^2$ makes friction and heat push every duct flow toward Mach 1 and choke it: «15.47», «15.48».",
        r"**C13** The Mach cone is the edge of all the circles of old news: «15.49».",
        r"**C14** An oblique shock is a normal shock on the normal component, «15.53»: two shocks per turn, a largest turn.",
        r"**C15** A supersonic stream turns away through a fan, and «15.56» keeps the count.",
        r"**C16** On a thin plate the push is normal to the plate, so lift costs α times itself in wave drag: «15.59».",
    )],
    feeds_forward=[X(r"""Ch. 16 — pulse waves in arteries are «15.15» with a wall-compliance wave speed; flow limitation in
collapsible tubes is choking (C09) with a tube law in place of the gas law.""")],
    left_out=["the 1-D Euler solver through breaking (N114; the backup explainer)",
              "the shock-structure derivations (stated only: N66–N69)", "the Riemann-invariant derivation (N22)",
              "non-perfect gases and expansion shocks (N90)", "the method of characteristics in two dimensions"],
)

# A "> ⚠️" callout written in the middle of a line prints a literal ">": start every such callout on its own line (a list
# item "4. > ⚠️ …" is a quote inside the item and stays as it is), and bold the "Common confusion" label.
_MID = re.compile(r"(?m)^(?P<pre>[^\n]*?\S)[ \t]+> ⚠️")


def _own_line(text: str) -> str:
    def fix(m):
        pre = m.group("pre")
        if re.fullmatch(r"\s*(\d+\.|[-*])", pre):     # a list marker: the quote belongs to the item
            return m.group(0)
        return pre + "\n\n> ⚠️"
    prev = None
    while prev != text:
        prev, text = text, _MID.sub(fix, text)
    return re.sub(r"> ⚠️ Common confusion", "> ⚠️ **Common confusion**", text)


for _c in nb.cells:
    if _c.cell_type == "markdown":
        _c.source = _own_line(_c.source)
assert not any(_MID.search(_c.source) and not re.search(r"(?m)^\s*(\d+\.|[-*]) > ⚠️", _c.source)
               for _c in nb.cells if _c.cell_type == "markdown")

path = nb.save()
print(f"wrote {path.relative_to(ROOT)}: {len(nb.cells)} cells "
      f"({sum(c.cell_type == 'markdown' for c in nb.cells)} markdown, {sum(c.cell_type == 'code' for c in nb.cells)} code), "
      f"{len(nb.cores)} CORE blocks, {len(nb.derivations)} derivations ({sum(d['steps'] for d in nb.derivations.values())} steps), "
      f"{len(nb.recaps)} recaps, {len(nb.explainers)} explainers, {len([p for p in nb.primers if not p.endswith('(reminder)')])} primers")
