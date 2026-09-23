"""Build the Chapter 5 teaching notebook: ``notebooks/ch05_vorticity_dynamics.ipynb``.

Source of truth: ``analysis/ch05_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract
— every call goes to ``fluidpy.ch05_vorticity_dynamics``, imported as ``ch05``, which re-exports ``core.vorticity``,
``core.biot_savart`` and the ch05 additions to ``core.vortices`` and ``core.integral_theorems``), Part E (prerequisite
ledger → 15 primers P134–P148 and one-line reminders of earlier primers), Part F (the 23 derivations D01–D23, one move
per step) and ``analysis/ch05_curation.md`` (IDs, depths, section coverage, §6 animations and plotly figures, §7
from-scratch moments). Physics lives in ``fluidpy``; cells only call it.

**Derivations are read from Part F at build time** (``part_f()`` below, the ch03/ch04 parser): goal, start, plan, tools,
assumptions, every step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word, so the
notebook and the design cannot drift apart. Explainer-only fields (*live*, *set*, *watch*) are dropped. The "check"
texts are edited (``pf_sub``) so that every check names a cell of this notebook that really runs it (ch04 lesson). The
four ★★★ sympy checks (D09, D10, D11, D15) are written here, every line commented, and re-run the construction.

Conventions (design header, binding): ω is the vorticity (the solid-body tank turns at ω/2); z up, g = 9.81 m/s²;
counterclockwise-positive circulation; ``Gamma`` [m²/s] is always a circulation, ``gamma`` [m/s] a sheet strength;
σ the viscous stress (σ_c a core radius); the book's slips are taught in corrected form ((5.14) +1/(4π), the sheet sign,
Fig. 5.2's rotation ω/2, Kelvin needs ρ = ρ(p), Fig. 5.11's G is the centre of vorticity, the exercise numbers).

Run:  .venv/Scripts/python.exe notebooks/build_ch05.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch05.py --dump     (prints the parsed Part F derivations only)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch05")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# (Ch. 5 transcribed from the rendered pages p198–p218 via the design (5.14 and 5.30 re-read here); earlier chapters'
# equations from their notebooks.)
# ---------------------------------------------------------------------------------------------------------------------
EQ = {
    "5.1": r"u_\theta=\omega r/2",
    "5.2": r"u_\theta=\frac{\Gamma}{2\pi r}",
    "5.3": r"dx/\omega_x=dy/\omega_y=dz/\omega_z",
    "5.4": r"\int_V\nabla\cdot\boldsymbol\omega\,dV=\int_A\boldsymbol\omega\cdot\mathbf n\,dA=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0",
    "5.5a": r"-\rho u_\theta^2/r=-\partial p/\partial r",
    "5.5b": r"0=-\partial p/\partial z-\rho g",
    "5.6": r"p(r,z)-p_o=\tfrac18\rho\omega^2r^2-\rho gz",
    "5.7": r"p(r,z)-p_\infty=-\frac{\rho\Gamma^2}{8\pi^2r^2}-\rho gz",
    "5.8": r"\frac{D\Gamma}{Dt}=0",
    "5.9": r"\frac{D\Gamma}{Dt}=\oint_C\frac{Du_i}{Dt}dx_i+\oint_Cu_i\frac{D}{Dt}(dx_i)",
    "5.10": r"\oint_C\frac{Du_i}{Dt}dx_i=-\oint_C\frac1\rho dp-\oint_Cd\Phi+\oint_C\Big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i",
    "5.11": r"\frac{D\Gamma}{Dt}=\oint_C\Big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\Big)dx_i",
    "5.12": r"\nabla\times\Big\{\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\mathbf g+\nu\nabla^2\mathbf u\Big\}",
    "5.13": r"\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega",
    "5.14": r"\mathbf u(\mathbf x,t)=+\frac1{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega(\mathbf x',t)}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'\ \ (\text{sign corrected; the book prints }-\tfrac1{4\pi})",
    "5.15": r"\int_{V'}\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)d^3x'=\int_{A'}\frac{\mathbf n\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^2x'",
    "5.16": r"\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\boldsymbol\omega(\mathbf x',t)\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'",
    "5.17": r"d\mathbf u\cong\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}",
    "5.18": r"\omega_{i,i}=\varepsilon_{inq}u_{q,ni}=0",
    "5.19": r"u_{i,i}=0",
    "5.20": r"\frac{\partial u_i}{\partial t}+u_ju_{i,j}+2\varepsilon_{ijk}\Omega_ju_k=-\frac1\rho p_{,i}+g_i+\nu u_{i,jj}",
    "5.21": r"u_ju_{i,j}=-(\mathbf u\times\boldsymbol\omega)_i+\tfrac12(u_j^2)_{,i}",
    "5.22": r"\varepsilon_{ijk}\omega_k=u_{j,i}-u_{i,j}",
    "5.23": r"\nu u_{i,jj}=-\nu\varepsilon_{ijk}\omega_{k,j}",
    "5.24": r"2\varepsilon_{ijk}\Omega_ju_k=-2\varepsilon_{ijk}\Omega_ku_j",
    "5.25": r"\frac{\partial u_i}{\partial t}+\Big(\tfrac12u_j^2+\Phi\Big)_{,i}-\varepsilon_{ijk}u_j(\omega_k+2\Omega_k)=-\frac1\rho p_{,i}-\nu\varepsilon_{ijk}\omega_{k,j}",
    "5.26": r"\frac{\partial}{\partial t}(\varepsilon_{nqi}u_{i,q})+\varepsilon_{nqi}\big(\tfrac12u_j^2+\Phi\big)_{,iq}-\varepsilon_{nqi}\varepsilon_{ijk}\big[u_j(\omega_k+2\Omega_k)\big]_{,q}=-\varepsilon_{nqi}\big(\tfrac1\rho p_{,i}\big)_{,q}-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}",
    "5.27": r"-\varepsilon_{nqi}\varepsilon_{ijk}\big[u_j(\omega_k+2\Omega_k)\big]_{,q}=-u_{n,j}(\omega_j+2\Omega_j)+u_j\omega_{n,j}",
    "5.28": r"-\varepsilon_{nqi}\Big(\frac1\rho p_{,i}\Big)_{,q}=\frac1{\rho^2}[\nabla\rho\times\nabla p]_n",
    "5.29": r"-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}=\nu\omega_{n,jj}",
    "5.30": r"\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega",
    "5.31": r"(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\frac{\partial\mathbf u}{\partial s}",
    "5.32": r"\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s},\ \frac{D\omega_n}{Dt}=\omega\frac{\partial u_n}{\partial s},\ \frac{D\omega_m}{Dt}=\omega\frac{\partial u_m}{\partial s}",
    "5.33": r"\frac{D\Gamma_a}{Dt}=0,\ \Gamma_a\equiv\int_A(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\mathbf n\,dA=\Gamma+2\int_A\boldsymbol\Omega\cdot\mathbf n\,dA",
    # earlier chapters
    "1.8": r"dp/dz=-\rho g",
    "2.19": r"\varepsilon_{ijk}\varepsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}",
    "2.30": r"\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA",
    "2.34": r"\int_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA=\oint_C\mathbf u\cdot d\mathbf x",
    "3.5": r"\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F",
    "3.7": r"dx/u=dy/v=dz/w",
    "3.15": r"R_{ij}=-\varepsilon_{ijk}\omega_k",
    "3.16": r"\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}",
    "3.18": r"\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA",
    "3.22": r"u_r=0,\ u_\theta=\omega_0r",
    "3.25": r"u_r=0,\ u_\theta=B/r",
    "3.28": r"u_\theta=\frac{\Gamma r}{2\pi\sigma^2}\ (r\le\sigma),\ \frac{\Gamma}{2\pi r}\ (r>\sigma)",
    "3.29": r"u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)",
    "4.10": r"\nabla\cdot\mathbf u=0",
    "4.18": r"\mathbf g=-\nabla\Phi",
    "4.24": r"\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}",
    "4.37": r"\tau_{ij}=-p\,\delta_{ij}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)+\mu_vS_{mm}\delta_{ij}",
    "4.39b": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u",
    "4.40": r"\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega\ \ (\nabla\cdot\mathbf u=0)",
    "4.41": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g",
    "4.45": r"\rho\frac{D'\mathbf u'}{Dt}=-\nabla'p+\rho\big[\mathbf g-2\boldsymbol\Omega\times\mathbf u'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')-\dots\big]+\mu\nabla'^2\mathbf u'",
    "4.58": r"\varepsilon=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}{\rho}S_{mm}^2\ge0",
    "4.59": r"\sigma_{ij}=\mu\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)+\big(\mu_v-\tfrac23\mu\big)\frac{\partial u_m}{\partial x_m}\delta_{ij}",
    "4.67": r"\frac1\rho\frac{\partial p}{\partial x_j}=\frac{\partial}{\partial x_j}\int_{p_o}^{p}\frac{dp'}{\rho(p')}",
    "4.68": r"u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\big(\tfrac12u_i^2\big)",
    "4.69": r"B=\tfrac12u^2+\int\frac{dp}{\rho}+gz,\ \ \frac{\partial\mathbf u}{\partial t}+\nabla B=\mathbf u\times\boldsymbol\omega",
    "4.71": r"\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const along streamlines and vortex lines}",
    "4.72": r"\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const everywhere}",
    "4.86": r"\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u",
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
                or nxt[:1] in "/–-)" or text[max(0, m.start() - 1):m.start()] == "("):
            return m.group(0)
        done.add(n)
        return f"({n}), ${EQ[n]}$"
    return _EQ_REF.sub(rep, text)


# ---------------------------------------------------------------------------------------------------------------------
# Part F reader: the derivations, word for word (the ch03/ch04 parser)
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
            tokens += [("w", w.replace("**", "").replace("`", "").replace("\\|", "|")) for w in p.split()]
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
    """"This is (5.9)." inside a step's *why* points at the line just displayed: say so instead of repeating it."""
    return _SELF.sub(lambda m: f"{m.group(1)}({m.group(2)}) — the line above", text)


def part_f() -> dict[str, dict]:
    """Parse Part F of the design into {D01: dict(title, goal, start, plan, tools, assumptions, steps, result, check,
    meaning, traps)}."""
    text = (ROOT / "analysis" / "ch05_design.md").read_text(encoding="utf-8")
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
    """Edit one Part F field after parsing (pointers to the cells that run a check); fails loudly if the design text
    changed. ``field`` = goal/check/meaning/traps/assumptions, ``tools``, or ``stepN.did|why|plain``."""
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


# Every "check" names a cell of this notebook that runs it (ch04 lesson review), with the numbers that cell prints.
pf_sub("D01", "check", "Numbers (narrowing Gaussian tube,",
       "Numbers (the code cell after the tiny example of this block; narrowing Gaussian tube,")
pf_sub("D02", "check", "puts the vertex 6.37 mm below the rest level (code).",
       "puts the vertex 6.37 mm below the rest level (the code cell after note N05 prints it).")
pf_sub("D03", "check", "(code: `vortex_stress_force(\"gaussian\", …)`)",
       "(the code cell after note N07: `vortex_stress_force(\"gaussian\", …)`)")
pf_sub("D03", "check", "The torque 2πr²σ_rθ = −2μΓ is the same at every radius (N08) ✓.",
       "The torque 2πr²σ_rθ = −2μΓ is the same at every radius (the code cell after note N08) ✓.")
pf_sub("D04", "check", "Numerically (`kelvin_rate_terms` on the cellular flow): contour term ≈ 1e-13, acceleration term = dΓ/dt from central differences of Γ(t) ✓.",
       "Numerically (the code cell after the tiny example of this block: `kelvin_rate_terms` on the cellular flow): contour term ≈ 3e-14, acceleration term ≈ 4e-13, and dΓ/dt from central differences of Γ(t) is ≈ 0 as well ✓.")
pf_sub("D05", "check", "Lamb–Oseen (Γ₀ = 0.01 m²/s,", "Lamb–Oseen (the code cell after the tiny example of this block; Γ₀ = 0.01 m²/s,")
pf_sub("D05", "check", "Baroclinic lock-exchange field: ∮dp/ρ round the square equals −∫(∇ρ×∇p/ρ²)·n dA to 1e-10 ✓.",
       "Baroclinic lock-exchange field (the code cell after note N14): −∮dp/ρ round the 0.2 m square equals ∫(∇ρ×∇p/ρ²)·n dA over it to the digits printed there ✓.")
pf_sub("D06", "check", "the numeric disc (4000 rim points, 400 × 800 area cells) agrees to 1e-9 relative and the error of the torque route falls as R² ✓.",
       "the numeric disc (4000 rim points, 400 × 800 area cells; the from-scratch cell of this block) agrees to the eight digits printed there and the error of the torque route falls as R² (the code cell after the tiny example) ✓.")
pf_sub("D07", "check", "Field route (`baroclinic_term` on the tanh step with hydrostatic p) = 2.4222 s⁻² at x = 0 ✓.",
       "Field route (`baroclinic_term` on the tanh step with hydrostatic p, the code cell right after this derivation) = 2.4222 s⁻² at x = 0 ✓.")
pf_sub("D08", "check", "Numerically (after C06):", "Numerically (the frozen-in cell at the end of C06, §5.4):")
pf_sub("D08", "check", "The small ABC loop on a tube wall keeps", "The small ABC loop on a tube wall (the code cell after the tiny example below) keeps")
pf_sub("D09", "check", "Lamb–Oseen: stretching 0, local = diffusion ✓. Burgers (steady): advective + stretching = diffusion, residual ≈ 1e-6 of the largest term at h = 1e-3 ✓.",
       "Lamb–Oseen: stretching 0, local = diffusion ✓. Burgers (steady): advective = stretching + diffusion, residual ≈ 1e-6 of the largest term (the code cell after the tiny example below) ✓.")
pf_sub("D10", "goal", "The book prints the answer (5.14) with a wrong sign",
       r"The book prints the answer, Eq. 5.14, as $\mathbf u=-\frac1{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'$ — with a wrong sign —")
pf_sub("D10", "check", "the + sign gives u_θ = +0.3183 m/s = Γ/2πr, the printed sign −0.3183 m/s ✓ (quadrature 24³–40³ nodes, error 1e-4).",
       "the + sign gives u_θ = +0.3088 m/s — the finite tube of length 4 m, whose exact value is 0.308838 m/s (the infinite line would give Γ/2πr = 0.318310) — and the printed sign −0.3088 m/s ✓ (the sympy cell below and the sign-test cell after note N25).")
pf_sub("D11", "check", "Rankine and Lamb–Oseen planar vortices recovered inside and outside the core by `biot_savart_2d` (1e-6) ✓. (5.15) on a random polynomial field: volume = surface to 1e-12 ✓.",
       "The Gaussian (Lamb–Oseen) planar vortex is recovered by `biot_savart_2d` (the code cell after the tiny example below: to ~10⁻⁸ m/s outside the core, 5×10⁻⁴ m/s inside on a 5 mm grid) ✓. (5.15) on a polynomial field: volume = surface to 1e-12 (the cell after note N27) ✓.")
pf_sub("D12", "check", "Summed round a ring", "Summed round a ring (the code cell after the tiny example of this block)")
pf_sub("D13", "check", "quad of the step-4 integral agrees to 1e-12 ✓.",
       "`quad` of the step-4 integral agrees with the closed form (the code cell after the tiny example below) ✓.")
pf_sub("D14", "check", "`rotating_lamb_form_terms` residual equals the (5.20) residual `rotating_ns_residual` to 1e-10 on random smooth fields ✓. Relabel check (5.24) with random numbers ✓.",
       "`rotating_lamb_form_terms` residual equals the (5.20) residual `rotating_ns_residual` on the inertial-oscillation field (the cell after note N32) ✓. Relabel check (5.24) with numbers (the cell after note N31) ✓.")
pf_sub("D15", "check", "Numbers: a column", "Numbers (the code cell after the tiny example below): a column")
pf_sub("D16", "check", "`stretching_tilting_split`:", "`stretching_tilting_split` (the code cell after the tiny example of this block):")
pf_sub("D17", "check", "Axial strain: ω_z(t) = ω₀e^{αt} (`uniform_strain_vorticity`)",
       "Axial strain (the code cell after the tiny example of this block): ω_z(t) = ω₀e^{αt} (`uniform_strain_vorticity`)")
pf_sub("D18", "check", "Term bars at R = 1 mm: advective + stretching = diffusion, residual ≈ 10⁻⁶ of the largest ✓ (`vorticity_terms`).",
       "Term bars at R = 1 mm (C06's code cell, §5.4): advective = stretching + diffusion, residual ≈ 10⁻⁶ of the largest ✓ (`vorticity_terms`); the code cell right after this derivation prints the core radius, the peak and Γ.")
pf_sub("D19", "check", "The rotating scenario of E3 (Ω = 0.5 rad/s, α = 0.2 s⁻¹)",
       "The rotating scenario of the Kelvin explainer (the code cell after the tiny example of this block; Ω = 0.5 rad/s, α = 0.2 s⁻¹)")
pf_sub("D19", "check", "Fluid at rest in the inertial frame seen from the rotating frame: Γ = −2ΩA, Γ_a = 0 ✓.",
       "Fluid at rest in the inertial frame seen from the rotating frame: Γ = −2ΩA, Γ_a = 0 (same cell) ✓.")
pf_sub("D20", "check", "Number: +10 % height at f = 10⁻⁴ s⁻¹ → ζ = +1.0×10⁻⁵ s⁻¹ ✓.",
       "Number: +10 % height at f = 10⁻⁴ s⁻¹ → ζ = +1.0×10⁻⁵ s⁻¹ (the code cell after the tiny example below) ✓.")
pf_sub("D21", "check", "the numerical integration (`point_vortex_evolve`) finds the same period",
       "the numerical integration (`point_vortex_evolve`, the code cell after the tiny example below) finds the same period")
pf_sub("D21", "check", "For Γ₂ = 3Γ₁ the fluid at G moves at", "For Γ₂ = 3Γ₁ the fluid at G (same cell) moves at")
pf_sub("D22", "check", "`wall_image_system`: v = 0 at 200 wall points to 1e-12 ✓.",
       "`wall_image_system` (the code cell after the tiny example below): v = 0 at 200 wall points to 1e-12 ✓.")
pf_sub("D23", "check", "the L1 error of u(y) against", "the L1 error of u(y) (the code cell after the tiny example below) against")

# Lesson review round 1: checks describe exactly what the cells run; tools stated before use; conventions written out.
pf_sub("D09", "check", "(u = (y z², x z, −x y)·t-free plus a Taylor–Green piece)",
       "(u = (y z² + sin x cos y, x z − cos x sin y, −x y): a polynomial plus the Taylor–Green (cellular) pattern — the identity checked here holds for any divergence-free field, so no time factor is needed)")
pf_sub("D09", "check", "take the curl of each term of (4.39b) separately (step 1)",
       "take the curl of each term of (4.39b) separately (step 1; `vorticity_equation_sym` returns each curl)")
pf_sub("D15", "check", "then that the extra terms vanish for a divergence-free polynomial u (steps 7, 9); the pressure term equals [∇ρ × ∇p]_n/ρ² (step 11); the viscous term equals νω_{n,jj} for that u (step 13); finally (5.30) − curl(5.25) ≡ 0.",
       "then that ∇·ω ≡ 0 for any u, so u_n ω_{k,k} drops (step 7); the pressure term equals [∇ρ × ∇p]_n/ρ² (steps 10–11); the viscous term equals νω_{n,jj} (steps 12–13); finally curl(5.25) − (5.30) equals exactly u_{j,j}(ω_n + 2Ω_n), which vanishes when u_{i,i} = 0 (step 9). Every line uses a generic u that is NOT assumed divergence-free — a stronger test than a single divergence-free field.")
pf_sub("D13", "check", "infinite 0.159155, semi-infinite 0.079577 m/s;",
       "infinite 0.159155, semi-infinite 0.079577 m/s (both printed in the code cell after the tiny example below);")
pf_sub("D21", "check", "and keeps ΣΓx, ΣΓ\\|x\\|² and H flat to 1e-10 ✓.",
       "and keeps ΣΓx, ΣΓ\\|x\\|² and H flat to 1e-10 (the cell prints all three spreads) ✓.")
pf_sub("D22", "check", "v = 0 at 200 wall points to 1e-12 ✓. Channel of height H: (Γ/4H)cot(πh/H) → Γ/4πh as h → 0 and 0 at h = H/2 ✓.",
       "v = 0 at 201 wall points to 1e-12 ✓. Channel of height H: (Γ/4H)cot(πh/H) → Γ/4πh as h → 0 and 0 at h = H/2 (the same cell prints both h = H/4 and h = H/2) ✓.")
# D12 step 3 (2a/|x − x′|, 20 % at 10 core radii) is now written that way in design Part F itself.
pf_sub("D07", "step7.why", "(convention 3)", "(our convention, set in the notation cell: ω_z > 0 means counterclockwise)")
pf_sub("D21", "assumptions", "(convention 3)", "(ω_z > 0 and Γ > 0 mean counterclockwise, as in the notation cell)")
pf_sub("D23", "step5.why", "(convention 4)", "(Γ in m²/s is always a circulation, γ in m/s a sheet strength — notation cell)")
pf_sub("D06", "step6.why", "(P116)", "(the moment-transfer rule glossed before this derivation)")
pf_sub("D05", "step9.why", "(by Stokes −∮dp/ρ = ∫(∇ρ×∇p/ρ²)·n dA)",
       "(by Stokes (2.34), −∮dp/ρ = −∮(∇p/ρ)·dx = ∫(∇ρ×∇p/ρ²)·n dA, because the curl of −∇p/ρ is ∇ρ×∇p/ρ² by the quotient rule ∇(1/ρ) = −∇ρ/ρ² — written out in C09, §5.6; the code cell after note N14 checks it)")
pf_sub("D19", "tools", "vector area of a closed loop (primer in C11)", "vector area of a closed loop (primer P138 in C05)")
pf_sub("D04", "check", "contour term ≈ 3e-14", "contour term ≈ −3e-14")
PF["D13"]["result"] = (r"u_\varphi=\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b);\qquad\text{infinite line: }u_\varphi=\frac{\Gamma}{2\pi d}\ \text{— the line vortex (5.2)}",
                       PF["D13"]["result"][1])

# Part F steps whose *why* was shorter than nbkit's minimum: the same reason, said in full (listed in the phase reply).
pf_sub("D08", "step5.why", "Stokes again, now on S′.",
       "Stokes' theorem (2.34) again, now applied to the carried patch S′, whose edge is the carried loop of step 4: its zero circulation becomes zero flux.")
pf_sub("D10", "step8.why", "Divide step 7 by −4π.",
       "Divide step 7 by −4π: the equation is linear, so scaling 1/r scales its point source by the same factor, and we want a source of strength exactly 1.")
pf_sub("D13", "step3.why", "Pythagoras: e_d ⟂ e_z.",
       "Pythagoras, because the two pieces d e_d and −l e_z of step 1 are perpendicular (e_d ⟂ e_z).")
pf_sub("D15", "step6.why", "Product rule (P38).",
       "Product rule (P38): the derivative of a product u_n(ω_k + 2Ω_k) is the derivative of each factor times the other.")
pf_sub("D16", "step2.why", "Linearity of the dot product.",
       "Linearity of the dot product: ω·(a + b + c) = ω·a + ω·b + ω·c, applied to the three terms of step 1.")
pf_sub("D17", "step2.why", "(5.31) from D16.",
       "Use (5.31) from D16: at a point where ω ≠ 0 the operator ω·∇ is the magnitude ω times the derivative along the vortex line.")
pf_sub("D18", "step6.why", "Substitute step 4.",
       "Substitute step 4 into step 5: u_R = −αR/2 and ∂u_z/∂z = α, both constant in time, turn it into an ODE in R.")
pf_sub("D20", "step2.why", "(5.33): the loop is material.",
       "(5.33) applies because the loop is material (it moves with the column's fluid) and the fluid is inviscid and barotropic.")

if "--dump" in sys.argv:
    for k, d in PF.items():
        print(f"=== {k} {d['title']}  ({len(d['steps'])} steps)")
        print("  START", d["start"][0][:150], "|", d["start"][1][:80])
        for i, s in enumerate(d["steps"], 1):
            print(f"  {i}. {s['did']} :: {s['tex'][:140]}")
        print("  RESULT", d["result"][0][:150], "|", d["result"][1][:80])
        print("  CHECK", d["check"][:400])
    sys.exit(0)


# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch04.py)
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
    s = re.sub(r"\s*`check_src` sketch.*$", "", s, flags=re.S)
    s = s.replace("**sympy check intent (★★★):**", "**sympy check (the cell right after this derivation, every line commented):**")
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


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + textwrap.dedent(text).strip())


def explainer(slug: str, heading: str, why_static: str, body: str, tries: list[str]) -> None:
    """An embedded explainer, opening with the "why interactive rather than a static figure" sentence (ch04 lesson)."""
    why = f"**Why interactive rather than a static figure:** {textwrap.dedent(why_static).strip()}"
    if body:
        why += "\n\n" + textwrap.dedent(body).strip()
    nb.explainer(slug, heading=heading, why=show_eqs(why), tries=[show_eqs(t) for t in tries])


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
                return m.group(0)                        # a range like (5.19)–(5.23)
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
Vorticity — twice the local spin of the fluid — is the quantity that makes a flow swirl, and this chapter follows it:
where it lives (vortex lines and tubes, which cannot end:
$\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ *(Eq. 5.4)*), what it does
to pressure (the bowl of a spinning bucket, the funnel of a drain), when it is conserved (Kelvin:
$D\Gamma/Dt=0$ *(Eq. 5.8)*; Helmholtz: vortex lines move with the fluid), how it changes (it is carried, stretched,
tilted and diffused: $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$
*(Eq. 5.13)*; on a rotating planet with density contrasts it gains the planet's vorticity and a baroclinic source,
$\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$
*(Eq. 5.30)*), and how it makes velocity everywhere (Biot–Savart,
$\mathbf u=\frac1{4\pi}\int\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$
*(Eq. 5.16)*) — so that vortices push each other, slide along walls and roll up.
""",
    roadmap=[
        r"§5.1 vortex lines and tubes, $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ *(5.4)* (C01); the two basic vortices, their pressure $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$ *(5.6)* and $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$ *(5.7)* and their viscous stress (C02)",
        r"§5.2 Kelvin's circulation theorem $D\Gamma/Dt=0$ *(5.8)* and its proof (C03); barotropic vs baroclinic — the pressure torque (C04)",
        "§5.3 Helmholtz's four vortex theorems: vortex lines move with the fluid (C05)",
        r"§5.4 the vorticity equation $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ *(5.13)* (C06)",
        r"§5.5 velocity from vorticity: Poisson, Green's function, Biot–Savart *(5.16)* (C07) and the filament law $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ *(5.17)* (C08)",
        r"§5.6 the vorticity equation in a rotating frame *(5.30)* (C09); stretching and tilting *(5.32)* (C10); absolute circulation $D\Gamma_a/Dt=0$ *(5.33)* and the fluid column (C11)",
        "§5.7 point vortices move each other (C12); images: a wall is a mirror vortex; rings near walls and leap-frogging (C13)",
        r"§5.8 vortex sheets: strength = jump in tangential velocity, $\gamma=u_2-u_1$ (C14)",
    ],
    prerequisites=[
        r"vorticity $\boldsymbol\omega=\nabla\times\mathbf u$, spin ½ω, circulation $\Gamma=\oint\mathbf u\cdot d\mathbf x$ and the solid-body, line, Rankine and Gaussian vortices (Ch. 3 §3.4–3.5)",
        "Gauss' and Stokes' theorems, the ε–δ identity, the curl of a gradient (Ch. 2)",
        r"Navier–Stokes, Euler, the viscous force $-\mu\nabla\times\boldsymbol\omega$, the Lamb identity, the Bernoulli function, rotating frames and Boussinesq (Ch. 4)",
        "hydrostatics (Ch. 1)",
    ],
)
nb.explainer_index([
    ("vortex_tubes_cannot_end", "Can a vortex just stop in the middle of the water?",
     "C01: vortex lines and tubes; the strength is the same at every section"),
    ("vortex_pressure_funnel", "Why is a spinning bucket a bowl and a drain a funnel?",
     "C02: the pressure of the two basic vortices; stress without net force"),
    ("kelvin_material_loop", "Stretch and tangle a loop of dye — what stays the same?",
     "C03 C05: Kelvin's theorem and its hypotheses; Helmholtz"),
    ("baroclinic_torque", "How can density make fluid spin?",
     "C04: the pressure torque = ∇ρ×∇p/ρ²; the lock exchange"),
    ("vorticity_stretching_tilting", "How can a flow spin fluid faster without a torque?",
     "C06 C10: the vorticity equation, stretching and tilting, Burgers' vortex"),
    ("biot_savart_filament", "How does a vortex push water far away?",
     "C07 C08: Biot–Savart, the filament law, the sign of the Green's-function formula"),
    ("vorticity_equation_rotating", "What does a spinning planet add to the vorticity equation?",
     "C09 C11: the rotating-frame vorticity equation, the fluid column, absolute circulation"),
    ("point_vortex_lab", "A vortex cannot push itself — so how do vortices move?",
     "C12 C13: pairs, the centre of vorticity, images"),
    ("vortex_sheet_rollup", "What is a velocity jump made of?",
     "C14: sheet strength = jump; roll-up"),
])
nb.setup()
nb.code(r"""
import numpy as np                                      # arrays (Ch. 1 primer P03)
import sympy as sp                                      # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                         # static figures (Ch. 1 primer P01)
import plotly.graph_objects as go                       # rotatable 3-D figures (Ch. 1 P41, Ch. 2 P64)
from plotly.subplots import make_subplots               # several plotly panels in one figure
from fluidpy import ch05_vorticity_dynamics as ch05     # the tested chapter-5 module: every function cites its § and Eq.
from fluidpy.core import vortices, kinematics, constitutive, rotating, bernoulli, diffusion, tensors   # Ch. 1–4 machinery
from fluidpy.core import navier_stokes as ns            # Ch. 4's Navier–Stokes helpers (Lamb identity, exact fields)
from fluidpy.core import integral_theorems as itg       # Ch. 2's Gauss and Stokes helpers (loops, discs, boxes)
from fluidpy.core.interact import slider_figure         # plotly figure with a slider that works on the web page (P17)
from fluidpy.core.anim import animate, ffmpeg_path      # matplotlib animations (P16); show_animation came with the setup
from fluidpy.core.style import COLORS, savefig          # the house palette and a helper that saves PNGs to outputs/ch05
from tools.convergence import observed_order            # slope of log(error) vs log(step): the observed order (P13)
from scripts.ch05_drawings import disc_element, column_sketch, vortex_marks, sheet_circuit, biot_savart_geometry, helix_with_frame   # sketches drawn by our code
import logging, warnings                                # standard library: log levels and warning filters
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless font-substitution notes
warnings.filterwarnings("ignore", category=RuntimeWarning)             # 0/0 on the axis of a vortex is handled by NaN on purpose
ffmpeg_path()                                           # find ffmpeg once (sets matplotlib's animation.ffmpeg_path)
ANIM_RC = {"axes.titleweight": "bold", "figure.constrained_layout.use": False, "savefig.dpi": 80}   # movies: light settings, 80 dpi frames (faster, smaller)


def recolor(fig, colors, dashes=None):                  # give plotly slider traces the notebook's colours by trace name
    for tr in fig.data:                                 # every trace of every slider step
        if tr.name in colors:                           # a name we assigned a colour to
            tr.line.color = colors[tr.name]             # same colour meaning as in the matplotlib figures
            tr.marker.color = colors[tr.name]           # markers too (for "markers" traces)
        if dashes and tr.name in dashes:                # optional dash pattern ("dash", "dot")
            tr.line.dash = dashes[tr.name]
    return fig                                          # the same figure, restyled


print(len([n for n in dir(ch05) if not n.startswith("_")]), "public names in fluidpy.ch05_vorticity_dynamics")   # the toolbox
""", explain="""
1. Numerical, symbolic and plotting libraries (all primed in Ch. 1–2).
2. `ch05` is the chapter module; it re-exports the new core modules of this chapter (`core.vorticity` — vortex lines,
   tubes, circulation theorems, vorticity-equation terms — and `core.biot_savart` — velocity from vorticity, point
   vortices, rings, sheets), so one name covers the whole chapter. Every function is tested in `tests/test_ch05.py`.
3. The core modules of Ch. 1–4 are called by module name (vortex profiles, stencils, stresses, rotating frames).
4. `scripts/ch05_drawings.py` holds a few drawing helpers (no physics): the fluid disc of the baroclinic torque, the
   columns over a ridge, vortex marks, the sheet circuit, the Biot–Savart geometry and the helix with its frame.
5. `recolor` only restyles plotly traces so that a colour always means the same thing.
""")
nb.md(r"""
**Notation, conventions and colours used in this notebook.**

| Symbol | Meaning | Unit |
|---|---|---|
| $\boldsymbol\omega=\nabla\times\mathbf u$ | vorticity (**never** a rotation rate in this chapter) | 1/s |
| $u_r,\ u_\theta$ | polar velocity components | m/s |
| $\Gamma$ | circulation (of a loop, a vortex, a tube, a filament) | m²/s |
| $\gamma$ | strength of a vortex sheet (circulation per unit length) | m/s |
| $\omega_z$ | plane vorticity, counterclockwise positive | 1/s |
| $p,\ \rho$ | pressure, density | Pa, kg/m³ |
| $\mu,\ \nu=\mu/\rho$ | dynamic and kinematic viscosity | Pa s, m²/s |
| $\sigma_{ij}$ | viscous stress | Pa |
| $\boldsymbol\Omega$ | angular velocity of a rotating frame | rad/s |
| $2\boldsymbol\Omega$, $\boldsymbol\omega+2\boldsymbol\Omega$ | planetary vorticity, absolute vorticity | 1/s |
| $f=2\Omega\sin\varphi$, $\zeta$ | Coriolis parameter at latitude φ, relative vertical vorticity | 1/s |
| $\mathbf x'$ (with $\nabla'$, $V'$, $A'$) | the **source** point being integrated over in §5.5 | m |
| $(s,n,m)$ | natural coordinates along a vortex line (§5.6) | m |

> ⚠️ **ω is the vorticity.** In $u_\theta=\omega r/2$ the fluid turns at ω/2: a tank spinning at 1 rad/s has
> ω = 2 s⁻¹. Ch. 3's `solid_body_rotation(r, omega0)` takes the turning rate ω₀ = ω/2.

> ⚠️ **Γ and γ.** Γ (`Gamma`, m²/s) is always a circulation; the strength of a vortex sheet is γ (`gamma`, m/s) — the
> book writes Γ for it in §5.8.

> ⚠️ **σ and ε.** σ is the viscous stress here (Ch. 3's Gaussian core radius is written σ_c); ε is the Levi-Civita
> symbol $\varepsilon_{ijk}$ except in note N08, where ε is the dissipation rate.

> ⚠️ **Primes in §5.5** mark the source point x′ being integrated over, not a moving frame.

**Colours** (one meaning each, also in the explainers): relative vorticity **teal** · planetary vorticity **amber** ·
baroclinic **orange** · diffusion **rose** · stretching **purple** · tilting **blue** · advection **grey** ·
counterclockwise (positive) vortices **teal**, clockwise **rose** · isobars and pressure **orange**, isopycnals and
density **blue** · references and ghosts **grey dashed**.
""")
nb.md(r"""
**Where this chapter is used later**

| Result here | Used in |
|---|---|
| tubes cannot end, (5.4), and Helmholtz's theorems | trailing vortices of wings (Ch. 14) |
| Kelvin, $D\Gamma/Dt=0$ (5.8) | irrotational flow stays irrotational (Ch. 6), the starting vortex and lift (Ch. 6, 14) |
| $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13) | vortex decay (Ch. 8), vorticity–stream-function CFD (Ch. 10), the turbulent cascade (Ch. 12) |
| Biot–Savart (5.16) and the filament law (5.17) | vortex methods (Ch. 10), lifting line and induced drag (Ch. 14) |
| (5.30), (5.33), the column $(\zeta+f)/h$ | potential vorticity, Rossby waves, baroclinic instability (Ch. 13) |
| vortex sheets | Kelvin–Helmholtz instability (Ch. 11), wakes (Ch. 14) |

*Climate hook:* Ch. 13's potential-vorticity conservation is the fluid column of C11 plus stratification; the sea breeze
and every front start with the baroclinic torque of C04.
""")
nb.md(r"""
**Book slips we correct in this notebook** (each is explained where it is used):

| The book prints | We use |
|---|---|
| printed (5.14): $\mathbf u=-\frac1{4\pi}\int\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'$ | the factor $+\frac1{4\pi}$ (derivation D10) |
| the integrand rewrite's $+\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\times\boldsymbol\omega$ | $-\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\times\boldsymbol\omega$ (D11; the two slips cancel, so (5.16) is right) |
| Fig. 5.2's rotation label "2ω" | with $u_\theta=\omega r/2$ the tank turns at ω/2 (recap R02) |
| "ρ and p single valued" as the reason $\oint dp/\rho=0$ | $\rho=\rho(p)$ makes $dp/\rho$ exact (recap R09, D05) |
| "C in irrotational fluid ⇒ no viscous term" | only for incompressible, constant-μ flow (note N13) |
| the Lamb step's $\nabla(\mathbf u\cdot\mathbf u)$ | $\nabla(\tfrac12\mathbf u\cdot\mathbf u)$ (recap R11) |
| (5.26)'s Π | Φ (note N33) |
| (5.27) drops $u_{j,j}(\omega_n+2\Omega_n)$ silently | zero by continuity, said aloud (D15) |
| "Exercise 5.8" after (5.14) $\mathbf u=\frac{1}{4\pi}\int\frac{\nabla\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x^\prime\rvert}\,dV^\prime$ (sign corrected), "Exercise 5.11" after (5.33) $D\Gamma_a/Dt=0$ | Exercises 5.9 and 5.10 |
| Fig. 5.16's caption $d\Gamma/ds=u_1-u_2$ | the text's $u_2-u_1$, counterclockwise positive (D23) |
| Fig. 5.11's caption: G is where the induced velocity vanishes | G is the centre of vorticity; the fluid there is at rest only for equal strengths (D21) |
| "hyperboloids of the second degree" (Fig. 5.3) | cubic surfaces $(c-z)r^2=$ const (note N07) |
""")

# =====================================================================================================================
# A.1 §5.1 Introduction — R01, C01, R02–R08, C02
# =====================================================================================================================
nb.section("5.1", "Introduction", intro=r"""
**What is this section about?** Vorticity, $\boldsymbol\omega=\nabla\times\mathbf u$, is twice the spin of a fluid
particle. This section gives the words for where it sits — vortex lines, vortex tubes and their strength — proves that
a tube cannot end inside the fluid, and then looks hard at the two simplest vortices: the rigidly spinning bucket and
the ideal line vortex. Their pressure fields explain a bowl-shaped and a funnel-shaped surface, and their viscous
stresses teach that irrotational does not mean stress-free.
""")
remind([
    ("functions as arguments and lambda", "`lambda x, t: …` is a one-line function; fluidpy takes fields as functions of (x, t) (Ch. 1 P29)."),
    ("numpy arrays and vectorised arithmetic", "vectors and tables of numbers that do arithmetic element by element, no loops needed (Ch. 1 P03)."),
    ("f-strings", "`f\"{x:.4f}\"` prints a number with 4 decimals inside a sentence (Ch. 1 P04)."),
])
nb.recap("R01", "Vorticity and spin", r"""
Vorticity is the curl of the velocity, $\boldsymbol\omega=\nabla\times\mathbf u$, and a small fluid element spins at
½ω (Ch. 3: the average turning rate of any two perpendicular material lines). A concentration of nearly parallel
vorticity is a *vortex*; motion on nearly circular streamlines is *vortex motion*. Below, `kinematics.vorticity` takes
the curl of a velocity function by central differences (Ch. 3 §3.4).
""", where="Ch. 3 §3.4")
nb.code(r"""
u_sb = lambda x, t: np.array([-1.0*x[1], 1.0*x[0], 0*x[0]])        # solid body turning at 1 rad/s: u = (−y, x, 0) [m/s]
u_lv = vortices.vortex_velocity_field("line", Gamma=2*np.pi)        # the ideal line vortex with Γ = 2π m²/s (Ch. 3)
print(kinematics.vorticity(u_sb, np.array([0.3, 0.2, 0.0]))[2])     # ω_z of the solid body at a point: 2.0 1/s
print(kinematics.vorticity(u_lv, np.array([1.0, 0.5]))[-1])          # ω_z of the line vortex off its axis: ≈ 0 (1e-9)
""", explain="""
The turning tank has vorticity twice its rotation rate everywhere (ω = 2 s⁻¹ for 1 rad/s); the line vortex has none
off its axis (the 10⁻⁹ is the round-off of the central differences).
""")

# ---- C01 ----------------------------------------------------------------------------------------------------------
core("C01", "Vortex tubes cannot end", r"""
A smoke ring, a tornado, the swirl down a drain — each is a tube of spinning fluid. Can such a tube simply fade out in
the middle of the air or water?
""", eqs=("5.4",))
note("N01", "The chapter's programme", r"""
Think of vorticity as carried by fluid elements: an element's vorticity can be turned, concentrated or spread by how
the element moves and deforms and by the torques of its neighbours. The rest of the chapter makes each of these precise
— carried (C06), stretched and tilted (C10), twisted by baroclinic torques (C04, C09), spread by viscosity (N22).
""")
problem(r"""
When a weather forecaster tracks a cyclone, an engineer a wing-tip vortex or an oceanographer an eddy, they follow a
tube of concentrated vorticity. What rules does such a tube obey whatever the flow does? The first rule is purely
geometric, and it is strict: the tube's strength is the same all along it — so it can thin, fatten, bend and twist, but
it cannot stop.
""")
idea("""
stream tube (Ch. 3)                      vortex tube (here)
walls made of streamlines                walls made of vortex lines
no flow crosses the wall                 no vorticity crosses the wall
div u = 0  =>  same volume flux Q        div w = 0 ALWAYS  =>  same strength Gamma at every section   (5.4)
thinner  =>  faster flow                 thinner  =>  stronger spin (mean w = Gamma / A)
""", r"""
**A vortex tube carries a fixed amount of spin, like a hose carrying water — squeeze it and the spin speeds up.**
""")
remind([
    ("parametric curves, tangent vector and arc length", r"a curve $\mathbf x(s)$ with arc length s has unit tangent $d\mathbf x/ds$; `solve_ivp` traces it (Ch. 3 P92)."),
    ("parallel vectors and the cross-product test", r"two vectors are parallel when their cross product vanishes — the ratio form $dx/\omega_x=dy/\omega_y=dz/\omega_z$ says the same (Ch. 3 P93)."),
    ("right-hand rule and orientation of a loop", "curl the fingers along the loop's direction; the thumb is the positive normal (Ch. 2 P74)."),
])
note("N02", "Vortex lines", r"""
are drawn tangent to ω everywhere, exactly as streamlines $dx/u=dy/v=dz/w$ *(3.7)* are drawn tangent to u:
""", equation=EQ["5.3"], ref="5.3")
nb.md(r"""
The ratio form breaks when a component is zero; we trace lines with the arc-length form
$d\mathbf x/ds=\boldsymbol\omega/\lvert\boldsymbol\omega\rvert$ (Ch. 3 P92, P93). In solid-body rotation every
vertical line is a vortex line; in the line vortex only the axis is. Irrotational flow has no vortex lines at all.
Our test field below is a **helical swirl** $u_\varphi=aRz$ (a swirl that grows with height and radius), whose
vorticity is $\omega_R=-aR$, $\omega_z=2az$; its vortex lines satisfy $dR/(-aR)=dz/(2az)$, i.e. $zR^2=$ const.
""")
nb.code(r"""
u_h = ch05.helical_swirl_field(a=1.0)                               # u_φ = aRz with a = 1 1/(m s): a test velocity field [m/s]
om = ch05.vorticity_field(u_h)                                      # ω(x, t) of that field by central differences [1/s]
s, X = ch05.vortex_line(om, np.array([1.0, 0.0, 1.0]), s_max=1.5)   # one vortex line through (1, 0, 1) m: dx/ds = ω/|ω|  (5.3)
R = np.hypot(X[0], X[1])                                            # np.hypot(x, y) = √(x² + y²): distance from the z-axis [m]
print(X.shape, f"spread of z·R² along the line = {np.ptp(X[2]*R**2):.1e} m³")   # np.ptp = max − min: ≈ 0, so zR² is constant
""", explain=r"""
1. `helical_swirl_field` builds the test velocity; `vorticity_field` turns any velocity function into its vorticity
   function (central differences).
2. `vortex_line` solves the arc-length equation $d\mathbf x/ds=\boldsymbol\omega/\lvert\boldsymbol\omega\rvert$ with
   `solve_ivp` in both directions from the seed point; `X` holds 400 points of the line.
3. Along the traced line $zR^2$ varies by less than 10⁻⁹ — the analytic answer of $dx/\omega_x=dy/\omega_y=dz/\omega_z$ (5.3).
""")
note("N03", "Tube and strength", r"""
The vortex lines through a closed curve form a *vortex tube*. Its *strength* is the circulation round a loop on the
tube that encircles it once — by Stokes' theorem $\int_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA=\oint_C\mathbf u\cdot d\mathbf x$
*(2.34)*, the vorticity flux through a cross-section:
""", equation=r"\Gamma=\oint_C\mathbf u\cdot d\mathbf x=\int_A\boldsymbol\omega\cdot\mathbf n\,dA")
nb.md(r"""
| stream tube (Ch. 3) | vortex tube (here) |
|---|---|
| volume flux $dQ=\mathbf u\cdot\mathbf n\,dA$ | strength $d\Gamma=\boldsymbol\omega\cdot\mathbf n\,dA$ |
| same Q at every section when $\nabla\cdot\mathbf u=0$ | same Γ at every section, always |

**Number:** the Gaussian (Lamb–Oseen) vortex $u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma_c^2}\big)$ *(3.29)* with
Γ = 1 m²/s and σ_c = 0.1 m, measured on the circle r = σ_c: both routes give Γ(1 − e⁻¹) = 0.632 m²/s.
""")
remind([
    ("named results (NamedTuple TubeStrength, TubeFlux)", "fluidpy returns results with named fields (`ts.circulation`, `ts.flux`) instead of bare tuples (Ch. 4 P111)."),
])
nb.code(r"""
u_g = vortices.vortex_velocity_field("gaussian", Gamma=1.0, sigma=0.1)   # Gaussian vortex (3.29): Γ = 1 m²/s, σ_c = 0.1 m (plane field)
loop = itg.planar_loop((0.0, 0.0), None, radius=0.1, n=512)             # the circle r = σ_c, counterclockwise, 512 points
disc = itg.planar_disc((0.0, 0.0), None, radius=0.1)                    # the disc it bounds (midpoint cells, Ch. 2)
ts = ch05.vortex_tube_strength(u_g, loop, disc)                         # Γ by ∮u·dx and by ∫ω·n dA (ω from stencils of u)
print(f"circulation {ts.circulation:.6f}   flux {ts.flux:.6f}   difference {ts.diff:.1e} m²/s")
print(f"exact Γ(1 − 1/e) = {1.0*(1 - np.exp(-1)):.6f} m²/s")             # the enclosed circulation of (3.29) at r = σ_c
Pc = ch05.circle_loop_points((0.0, 0.0), 0.1, n=512)                    # the same circle as a (2, N) array of points
print(f"clockwise loop: {ch05.loop_circulation(u_g, Pc[:, ::-1]):+.6f} m²/s")   # reversed order = opposite orientation: sign flips
""", explain=r"""
1. The Gaussian vortex as a velocity function; a circle of radius σ_c and the disc inside it.
2. `vortex_tube_strength` measures the strength both ways: the line integral round the loop (the periodic trapezoid rule, primer P136 in C03 —
   spectrally accurate on a circle: 0.632121) and the flux of ω through the disc (the disc's midpoint cells and the stencil curl agree to four
   digits, 0.6322 — about 10⁻⁴ relative, the size of that quadrature's error).
3. Walking the loop the other way (points reversed) flips the sign: strength is measured with the right-hand rule.
""")
D("D01", ref="5.4")
nb.worked_example("a tube squeezed to a quarter of its area", r"""
A thin vortex tube of strength Γ = 0.01 m²/s has a cross-section of 1 cm² = 10⁻⁴ m² at one place.

1. Mean vorticity there: $\bar\omega=\Gamma/A=0.01/10^{-4}=100$ s⁻¹.
2. Further along, the section is 0.25 cm² = 2.5×10⁻⁵ m². By (5.4), $-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$,
   the strength is still 0.01 m²/s.
3. So $\bar\omega=0.01/(2.5\times10^{-5})=400$ s⁻¹ — four times the spin for a quarter of the area.
4. The fluid spins at ½ω̄: from 50 to 200 rad/s, about 8 to 32 turns per second.
""")
remind([
    ("np.meshgrid and the project grid layout", "`np.meshgrid(r, th, indexing=\"ij\")` gives every (r, θ) pair; with indexing='ij' the rows run along the first argument, with the default 'xy' along the second (Ch. 2 P76)."),
    ("volume and surface integrals as midpoint sums (polar grid)", "a surface integral is a sum of value × small area over cells; in polar cells the area is r dr dθ (Ch. 2 P83)."),
    ("assert np.allclose", "`assert np.allclose(a, b)` stops the notebook if two results differ beyond rounding (Ch. 1 P15)."),
])
nb.code(r"""
for z in (0.0, 0.5, 1.0):                                            # three sections along the narrowing tube [m]
    sec = ch05.gaussian_tube_section(z, R0=0.1, Gamma=1.0, a0=0.1, L=1.0)   # the tube through R0 = 0.1 m at z = 0
    print(f"z = {z:.1f} m: radius {sec['radius']:.4f} m, area {sec['area']:.6f} m², flux {sec['flux']:.6f} m²/s, "
          f"mean ω {sec['mean_omega']:.2f} 1/s")                     # flux fixed, area down, mean vorticity up
print(ch05.tube_flux_budget("gaussian_tube", 0.0, 1.0, 0.1, Gamma=1.0, a0=0.1, L=1.0))   # Gauss on the piece 0 ≤ z ≤ 1 m: (5.4)
print(ch05.tube_flux_budget("broken", 0.0, 1.0, 0.1, Gamma=1.0, a0=0.1, L=1.0))         # a field that is NOT a curl
""", explain=r"""
1. The **narrowing Gaussian tube** $\omega_z=(\Gamma/\pi a^2)e^{-R^2/a^2}$ with $a^2=a_0^2e^{-z/L}$ (plus the small inward
   ω_R that keeps it divergence-free): the tube through R₀ = 0.1 m keeps its flux 0.632121 m²/s while its area shrinks
   from 0.031416 to 0.011557 m² and its mean vorticity rises from 20.12 to 54.69 s⁻¹.
2. `tube_flux_budget` applies Gauss to the piece between z = 0 and 1 m: lower lid (outward normal against ω) −0.632121,
   wall ≈ 2×10⁻¹², upper lid +0.632121, total ≈ 0 — that is (5.4).
3. The **broken** field $\omega=(\Gamma/\pi a_0^2)e^{-R^2/a_0^2}e^{-z/L}\mathbf e_z$ has $\nabla\cdot\boldsymbol\omega\ne0$: the
   upper lid carries only 0.232544, the total is −0.399576 — "vorticity" would be created inside, so it cannot be the
   curl of any velocity.
""")
nb.md("**From scratch — our own midpoint sums over the two lids** (polar cells of area r Δr Δθ):")
nb.check_agree(r"""
om_t = ch05.gaussian_tube_field(Gamma=1.0, a0=0.1, L=1.0)          # ω(x, t) of the narrowing Gaussian tube [1/s]
nr, nt = (200, 400) if not FAST else (100, 200)                      # radial and angular cells
flux = []                                                            # ∫ω_z dA over each lid [m²/s]
for z in (0.0, 1.0):                                                 # bottom lid, then top lid [m]
    Rt = 0.1*np.exp(-z/2)                                            # the tube's radius at this height [m]
    r = (np.arange(nr) + 0.5)/nr*Rt                                  # midpoints in r [m]
    th = (np.arange(nt) + 0.5)/nt*2*np.pi                            # midpoints in θ [rad]
    RR, TT = np.meshgrid(r, th, indexing="ij")                       # all (r, θ) cell centres
    Pts = np.array([RR*np.cos(TT), RR*np.sin(TT), z + 0*RR]).reshape(3, -1)   # their (x, y, z) positions, shape (3, N)
    dA = (RR*(Rt/nr)*(2*np.pi/nt)).ravel()                           # cell areas r Δr Δθ [m²]
    flux.append(np.sum(om_t(Pts, 0.0)[2]*dA))                        # Σ ω_z dA
print(f"bottom lid {flux[0]:.6f}, top lid {flux[1]:.6f} m²/s")       # the same number twice
assert np.allclose(flux[0], flux[1], rtol=1e-6)                      # (5.4): same strength at both sections
assert np.allclose(flux[1], ch05.tube_flux_budget("gaussian_tube", 0.0, 1.0, 0.1).upper, rtol=1e-4)   # = the library
""")
nb.md("Our own sums over the two lids agree with each other and, to the accuracy of a midpoint rule, with the library.")
remind([
    ("plotly 3-D lines and surfaces", "`go.Scatter3d(mode=\"lines\")` draws a 3-D curve you can rotate with the mouse (Ch. 2 P64)."),
])
nb.md(r"""
**Stream tube vs vortex tube (Fig. 5.1 analogue, note `N45`).** On the left, a stream tube in **Burgers' vortex** — a vortex fed
by a flow drawn in sideways and stretched upward (derived in C10); on the right, a vortex tube of the narrowing
Gaussian field, with its lids at z = 0 (blue) and z = 1 m (orange).
""")
nb.plotly(r"""
nL = 12 if not FAST else 8                                           # lines per tube
fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                    subplot_titles=("stream tube (Burgers) [mm]", "vortex tube (Gaussian) [m]"))
u_burg = vortices.burgers_vortex_field(1e-3, 1.0, 1e-6)             # Burgers: Γ = 1e-3 m²/s, α = 1 1/s, ν = 1e-6 m²/s
for k in range(nL):                                                  # seeds on a circle of radius 3 mm at z = 1 mm
    th = 2*np.pi*k/nL                                                # seed angle [rad]
    x0 = np.array([3e-3*np.cos(th), 3e-3*np.sin(th), 1e-3])          # seed point [m]
    S = kinematics.streamline(u_burg, x0, 0.0, 0.12, False, 300)     # a streamline: tangent to u, 0.12 m of arc (3.7)
    fig.add_trace(go.Scatter3d(x=S[0]*1e3, y=S[1]*1e3, z=S[2]*1e3, mode="lines", showlegend=False,
                               line=dict(color=COLORS["blue"], width=3)), row=1, col=1)   # in mm
tube = ch05.vortex_tube(ch05.gaussian_tube_field(), (0, 0, 0), (0, 0, 1), 0.1, n_lines=nL, s_max=1.2, n=150, both=False)
for Lk in tube.lines:                                                # each vortex line of the tube, shape (3, n)
    fig.add_trace(go.Scatter3d(x=Lk[0], y=Lk[1], z=Lk[2], mode="lines", showlegend=False,
                               line=dict(color=COLORS["teal"], width=4)), row=1, col=2)
ph = np.linspace(0, 2*np.pi, 80)                                     # angles round a lid
for z, col in ((0.0, COLORS["blue"]), (1.0, COLORS["orange"])):      # lower lid blue, upper lid orange (D01 colours)
    Rl = 0.1*np.exp(-z/2)                                            # the tube's radius there [m]
    fig.add_trace(go.Scatter3d(x=Rl*np.cos(ph), y=Rl*np.sin(ph), z=z + 0*ph, mode="lines", showlegend=False,
                               line=dict(color=col, width=8)), row=1, col=2)
fig.update_layout(height=480, margin=dict(l=0, r=0, t=40, b=0), title="Two tubes: one carries volume, one carries spin")  # size, margins, title
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="Left, blue streamlines that spiral round the axis and crowd together as they rise; right, teal vortex lines "
        "that converge toward the axis between a wide blue lid and a narrow orange one.",
    read="Each tube's wall is made of its own lines, so nothing crosses it: the stream tube carries a fixed volume flux, "
         "the vortex tube a fixed strength — the orange lid carries exactly the blue lid's flux on 37 % (e⁻¹) of the area.",
    change="…the twist of the vortex field turned on (`ch05.gaussian_tube_field(twist=5.0)`): the vortex lines become "
           "helices, the lids' fluxes do not change.")
remind([
    ("log axes", "on a logarithmic axis equal ratios are equal distances, so a quantity that doubles and one that halves look symmetric (Ch. 1 P13)."),
    ("matplotlib figures and labels with units", "`fig, ax = plt.subplots()`, `ax.plot`, axis labels always with units (Ch. 1 P01)."),
])
nb.figure(r"""
zz = np.linspace(0.0, 1.5, 31)                                       # heights along the tube [m]
secs = [ch05.gaussian_tube_section(z) for z in zz]                   # R0 = a0 = 0.1 m, Γ = 1 m²/s, L = 1 m
flux = np.array([s["flux"] for s in secs])                           # tube strength [m²/s]
area = np.array([s["area"] for s in secs])                           # section area [m²]
wbar = np.array([s["mean_omega"] for s in secs])                     # mean vorticity [1/s]
broken = np.array([ch05.tube_flux_budget("broken", 0.0, z, 0.1).upper if z > 0 else flux[0] for z in zz])   # "flux" of the broken field
fig, ax = plt.subplots(figsize=(7, 3.8))                             # one panel
ax.semilogy(zz, flux/flux[0], color=COLORS["teal"], lw=2.5, label="strength Γ (flux)")                   # flat
ax.semilogy(zz, area/area[0], color=COLORS["blue"], lw=2, label="section area A")                        # falls as e^{−z/L}
ax.semilogy(zz, wbar/wbar[0], color=COLORS["accent"], lw=2, label="mean vorticity Γ/A")                   # rises as e^{z/L}
ax.semilogy(zz, broken/broken[0], color=COLORS["rose"], ls="--", lw=2, label="broken field (∇·ω ≠ 0)")   # falls: not a tube
ax.set_xlabel("height along the tube $z$ [m]"); ax.set_ylabel("value / value at $z=0$ [–]")               # axis labels with units
ax.set_title("The strength is pinned; the spin pays for the squeeze"); ax.legend(fontsize=8)
savefig(fig, "ch05", "tube_strength"); plt.show()                    # save the PNG to outputs/ch05, then draw
""", see="A flat teal line, a falling blue curve and a rising purple one that mirror each other; the dashed rose line "
         "falls like the area.",
    read=r"Flux = area × mean vorticity stays constant along the tube, $-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ "
         "(5.4); only a field that is not a curl lets the \"strength\" change along a tube.",
    change="…L = 0.5 m (a faster narrowing): the area and mean-ω curves steepen twice as fast; the flux line does not move.")
remind([
    ("show_viz and embedded explainers", "`show_viz(\"ch05\", slug)` displays an interactive explainer in the output; on the web page it fills the window (Ch. 1 P18)."),
])
explainer("vortex_tubes_cannot_end", "Can a vortex just stop in the middle of the water?",
          "a static figure shows one section; here you slide the section along a narrowing, twisting tube and watch "
          "the flux stay pinned while the area and the spin trade off — then break the field and watch Gauss' budget fail.",
          "", [
              "Drag the section from z = 0 to z = 1 m and watch the lid-flux bars: they stay equal while the mean ω grows by e ≈ 2.7.",
              "Switch to the twisted tube: the lines spiral, the wall bar stays at zero.",
              "Choose the broken field: the upper lid carries less than the lower — the badge says why this cannot be a vorticity field.",
              "Open the Derivation tab and step through D01 with the three coloured pieces (lower lid blue, wall grey, upper lid orange).",
          ])
nb.md(r"""
> ⚠️ **Common confusion:** "the side of the tube carries no flux because the vorticity is zero there." The vorticity on
> the wall can be large; it simply points *along* the wall ($\boldsymbol\omega\cdot\mathbf n=0$).
""")
whatif(r"""
…the fluid were viscous? Nothing in (5.4) changes: it is geometry, true for every flow, because
$\nabla\cdot(\nabla\times\mathbf u)=0$ always. Viscosity will change *which* lines form a tube as time goes on (C05) —
but at every instant each tube has one strength.
""")

# ---- recaps for C02 ----------------------------------------------------------------------------------------------
nb.recap("R02", "Solid-body rotation from uniform vorticity", r"""
A uniform plane-normal vorticity ω makes the fluid turn rigidly: $u_\theta=\omega r/2$ *(Eq. 5.1)* — Ch. 3's
$u_r=0,\ u_\theta=\omega_0r$ *(3.22)* with the rotation rate ω₀ = ω/2. ⚠️ Here ω is the vorticity: a tank turning at
Ω_tank = 5 rad/s has ω = 10 s⁻¹. (Fig. 5.2 labels the tank's rotation "2ω" — with $u_\theta=\omega r/2$ the tank turns
at ω/2.)
""", where="Ch. 3 §3.5")
nb.code(r"""
print(ch05.solid_body_from_vorticity(0.1, 10.0))                    # (u_θ, ω_z) at r = 0.1 m for ω = 10 1/s: (0.5 m/s, 10 1/s)
print(vortices.solid_body_rotation(0.1, 5.0))                        # Ch. 3's function takes the turning rate 5 rad/s: 0.5 m/s
""")
nb.recap("R03", "The ideal line vortex", r"""
A perfect concentration of vorticity on the axis with circulation Γ gives $u_\theta=\Gamma/2\pi r$ *(Eq. 5.2)* — Ch. 3's
$u_r=0,\ u_\theta=B/r$ *(3.25)* with Γ = 2πB; irrotational for r > 0, and $\oint\mathbf u\cdot d\mathbf x=\Gamma$ on every
circle round the axis.
""", where="Ch. 3 §3.5")
nb.code(r"""
print(ch05.line_vortex_gamma(0.1, 1.0))                              # (u_θ, ω_z) at r = 0.1 m for Γ = 1 m²/s: (1.5915 m/s, 0)
print([vortices.circulation_circle(lambda r: ch05.line_vortex_gamma(r, 1.0)[0], r) for r in (0.1, 1.0)])   # Γ on two circles: 1, 1
""")
nb.recap("R04", "Paddle wheels in the two vortices", r"""
Both vortices are steady with circular streamlines, but only the first turns a paddle wheel; in the line vortex a small
wheel travels round the axis without spinning. The Ch. 3 explainer `vortex_paddle_wheels` shows it (open it from the
Ch. 3 page; it is not embedded again).
""", where="Ch. 3 §3.5")
nb.recap("R05", "No deformation, no viscous stress", r"""
Solid-body rotation does not deform elements, $S_{ij}=0$, so the Newtonian stress
$\tau_{ij}=-p\,\delta_{ij}+2\mu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)+\mu_vS_{mm}\delta_{ij}$ *(4.37)* is pure pressure,
$\tau_{ij}=-p\delta_{ij}$, and Cauchy's equation $\rho\frac{Du_j}{Dt}=\rho g_j+\frac{\partial\tau_{ij}}{\partial x_i}$ *(4.24)*
reduces to Euler's $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g$ *(4.41)*.
""", where="Ch. 4 §4.5–4.6")
nb.code(r"""
G_sb = np.array([[0, -5.0, 0], [5.0, 0, 0], [0, 0, 0]])               # velocity gradient of a tank turning at 5 rad/s [1/s]
print(constitutive.newtonian_stress(G_sb, p=1.0e5, mu=1e-3))         # stress [Pa]: diag(−1e5, −1e5, −1e5), no shear — pure pressure
""")
nb.recap("R06", "Hydrostatics", r"""
With no vertical motion the vertical momentum balance is hydrostatic: $0=-\partial p/\partial z-\rho g$ *(Eq. 5.5b)* —
Ch. 1's $dp/dz=-\rho g$ *(1.8)*.
""", where="Ch. 1 §1.7")
nb.recap("R07", "The Bernoulli function across streamlines", r"""
$B=\tfrac12u^2+\int dp/\rho+gz$ (4.69) is constant along streamlines in steady inviscid barotropic flow, and constant
everywhere only if the flow is also irrotational — $\tfrac12u_i^2+\int_{p_o}^{p}\frac{dp'}{\rho(p')}+gz=\text{const everywhere}$
*(4.72)*. Ch. 4's `which_bernoulli` explainer shows the four versions.
""", where="Ch. 4 §4.9")
nb.recap("R08", "The Rankine vortex", r"""
Uniform vorticity inside r ≤ a, a line vortex outside: $u_\theta=\frac{\Gamma r}{2\pi\sigma^2}\ (r\le\sigma),\ \frac{\Gamma}{2\pi r}\ (r>\sigma)$
*(3.28)* with σ = a. A solid cylinder of radius a spinning at ω/2 in viscous fluid drives exactly the outer part,
$u_\theta=\omega a^2/2r$, with Γ = πa²ω (Ch. 8 solves it).
""", where="Ch. 3 §3.5")
nb.code(r"""
print(ch05.rotating_cylinder_flow(np.array([0.05, 0.1, 0.2]), 0.1, 10.0))   # (u_θ, ω_z) at r = 5, 10, 20 cm; a = 0.1 m, ω = 10 1/s
print(f"Γ = πa²ω = {np.pi*0.1**2*10.0:.4f} m²/s")                    # the cylinder's circulation
""")

# ---- C02 ----------------------------------------------------------------------------------------------------------
core("C02", "The bowl and the funnel: pressure in the two basic vortices", r"""
Spin a bucket of water and its surface becomes a bowl; pull the plug of a bath and a narrow funnel appears. Why the
different shapes — and which of the two vortices is the viscous one?
""", eqs=("5.6", "5.7"))
problem(r"""
Every vortex in nature, from a stirred cup to a tornado to a hurricane's eye wall, keeps its fluid on curved paths, and
something must push that fluid inward: pressure, lower at the centre. How much lower depends on how the speed varies
with radius. The same two basic profiles, $u_\theta=\omega r/2$ (5.1) and $u_\theta=\Gamma/2\pi r$ (5.2), also settle a
subtle question about viscosity: the vortex whose fluid does *not* rotate is the one with viscous stress.
""")
idea("""
every parcel on a circle needs an inward push  rho u^2/r  ->  pressure rises outward: dp/dr = rho u^2/r        (5.5a)
solid body  u = w r/2       : push grows with r          ->  p ~ r^2    -> a BOWL (paraboloid)     no viscous stress
line vortex u = Gamma/2pi r : push ~ 1/r^3 near the axis ->  p ~ -1/r^2 -> a deep FUNNEL          stress, zero net force
""")
remind([
    ("cylindrical and spherical unit vectors", r"$\mathbf e_r,\ \mathbf e_\theta,\ \mathbf e_z$ point outward, round and up; they change direction from point to point (Ch. 3 P88)."),
    ("derivative of a rotating unit vector ∂e_θ/∂θ = −e_r", r"moving round a circle, $\mathbf e_\theta$ turns toward the axis: $\partial\mathbf e_\theta/\partial\theta=-\mathbf e_r$ (Ch. 4 P124)."),
    ("sympy symbols, diff, simplify, Function", "`sp.symbols`, `sp.Function('f')(z)` (an unknown function), `sp.diff`, `sp.simplify` (Ch. 1 P40)."),
], lead="needed in the derivation below")
P("P134", "partial integration with an unknown function", r"""
When you integrate a partial derivative ∂p/∂r with respect to r, the "constant" of integration may depend on the other
variable z, because ∂/∂r of any f(z) is zero. So ∂p/∂r = 2r gives p = r² + f(z), and a second equation (for ∂p/∂z)
fixes f. Derivation D02 uses exactly this.
""", code=r"""
r, z = sp.symbols('r z')                     # two variables
f = sp.Function('f')                         # an unknown function of z alone
p = sp.integrate(2*r, r) + f(z)              # ∂p/∂r = 2r integrated in r, plus the "constant" f(z)
print(p, sp.diff(p, r))                      # r**2 + f(z), 2*r: the f(z) is invisible to ∂/∂r
""")
note("N04", "Viscosity re-enters", r"""
Viscosity diffuses vorticity and lets vortex lines *reconnect* — cut and rejoin — which an ideal fluid forbids (C05).
The book only names reconnection; it returns with turbulence (Ch. 12). Here we look at the two basic vortices with
viscosity switched on.
""")
D("D02", ref="5.6")
note("N05", "The radial balance", r"""
D02's step 5 is (5.5a): pressure supplies the centripetal acceleration. At r = 0.1 m in a tank with ω = 10 s⁻¹
($u_\theta$ = 0.5 m/s): ∂p/∂r = 1000 × 0.25/0.1 = 2500 Pa/m.
""", equation=EQ["5.5a"], ref="5.5a")
remind([
    ("scipy.optimize.brentq", "finds the root of a function between two points where it changes sign (Ch. 3 P108) — used for the dry radius of a fast-spun tank."),
])
nb.code(r"""
print(ch05.solid_body_pressure_gradients(0.1, 0.0, 10.0))           # (∂p/∂r, ∂p/∂z) at r = 0.1 m, ω = 10 1/s: (2500, −9810) Pa/m
print(ch05.solid_body_pressure(0.1, 0.0, 10.0), ch05.isobar_height(0.1, 10.0, kind="solid"))   # p − p_o = 125 Pa; isobar through the origin 12.74 mm up
tank = ch05.rotating_tank_free_surface(0.1, 0.1, 5.0)                # tank radius 0.1 m, 0.1 m deep, spun at 5 rad/s
print({k: round(tank[k], 6) for k in ("z_vertex", "z_rim", "rim_rise", "centre_drop", "r_dry")}, tank["spills"])   # heights [m]
R, Z, w, rho, g = sp.symbols('r z omega rho g', positive=True)      # symbols for the sympy check
p = rho*w**2*R**2/8 - rho*g*Z                                        # (5.6) with p_o = 0
print(sp.simplify(sp.diff(p, R) - rho*(w*R/2)**2/R), sp.simplify(sp.diff(p, Z) + rho*g))   # both 0: (5.5a) and (5.5b) hold
""", explain=r"""
1. The two gradients of $-\rho u_\theta^2/r=-\partial p/\partial r$ (5.5a) and $0=-\partial p/\partial z-\rho g$ (5.5b).
2. The pressure of $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$ (5.6) at r = 0.1 m (125 Pa) and the height of the isobar
   through the origin there (12.742 mm = ω²R²/8g).
3. A real tank of radius 10 cm filled 10 cm deep, spun at 5 rad/s: the water keeps its volume, so the centre drops by
   exactly as much as the rim rises, 6.371 mm each (vertex at 0.093629 m, rim at 0.106371 m); nothing dries, nothing spills.
4. sympy confirms (5.6) satisfies both (5.5a) and (5.5b).
""")
remind([
    ("trapezoid rule and cumulative_trapezoid", "`scipy.integrate.cumulative_trapezoid(y, x, initial=0)` returns the running integral ∫y dx at every x (Ch. 1 P37)."),
])
nb.md("**From scratch — add up $\\rho u_\\theta^2/r\\,dr$ ourselves** and compare with the library's bowl and funnel:")
nb.check_agree(r"""
from scipy.integrate import cumulative_trapezoid                     # running trapezoid sum
r = np.linspace(0, 0.1, 2001)                                        # radii from the axis to 10 cm [m]
dpdr = 1000*(10*r/2)**2/np.where(r > 0, r, 1)                        # ρu_θ²/r for the tank, ω = 10 1/s [Pa/m] (0 on the axis)
p_num = cumulative_trapezoid(dpdr, r, initial=0.0)                   # p(r) − p(0) [Pa]
assert np.allclose(p_num, ch05.solid_body_pressure(r, 0.0, 10.0), atol=1e-6)   # = (5.6): the bowl
rr = np.linspace(1.0, 0.1, 40001)                                    # the funnel: from r = 1 m inward to 0.1 m [m]
dp = 1000*(1/(2*np.pi*rr))**2/rr                                     # ρu_θ²/r for Γ = 1 m²/s [Pa/m]
p_f = ch05.line_vortex_pressure(1.0, 0.0, 1.0) + cumulative_trapezoid(dp, rr, initial=0.0)   # integrate inward from p(1 m)
print(f"funnel at r = 0.1 m: ours {p_f[-1]:.3f} Pa, library {ch05.line_vortex_pressure(0.1, 0.0, 1.0):.3f} Pa")
assert np.allclose(p_f[-1], ch05.line_vortex_pressure(0.1, 0.0, 1.0), rtol=1e-6)   # = (5.7): the funnel
""")
nb.md(r"""
**Recap R07 used: B across streamlines.** With $\rho\omega^2r^2/8=\rho u_\theta^2/2$, (5.6) reads
$-\tfrac12u_\theta^2+gz+\frac{p(r,z)}{\rho}=\text{const}$ — so the Bernoulli function $B=u_\theta^2/2+gz+p/\rho$ **grows**
outward, $B-B(0)=\omega^2r^2/4$, as it must in a rotational flow (4.71)–(4.72). For the line vortex B is the same
everywhere (irrotational).
""")
nb.code(r"""
print(ch05.bernoulli_across_vortex("solid", 0.1, param=10.0))        # B(0.1 m) − B(0) in the tank: ω²r²/4 = 0.25 m²/s²
print(ch05.bernoulli_across_vortex("line", 0.1, param=1.0, r_ref=1.0))   # line vortex: B(0.1 m) − B(1 m) = 0 (round-off)
""")
note("N06", "The line vortex is sheared", r"""
Its viscous stress is nonzero everywhere because elements deform (Ch. 3's picture of an element sheared as it goes
round) — yet the net viscous force on every element is zero (the book leaves this to Exercise 5.4; D03 below does it
three ways). ⚠️ σ here is the viscous stress, not Ch. 3's core radius.
""", equation=r"\sigma_{r\theta}=\mu\Big[\frac1r\frac{\partial u_r}{\partial\theta}+r\frac{\partial}{\partial r}\Big(\frac{u_\theta}{r}\Big)\Big]=-\frac{\mu\Gamma}{\pi r^2}")
nb.md(r"""
*Gloss — polar strain rate and polar stress divergence (used in D03):* in polar coordinates the shear strain rate is
$S_{r\theta}=\tfrac12\big[\frac1r\frac{\partial u_r}{\partial\theta}+r\frac{\partial}{\partial r}(\frac{u_\theta}{r})\big]$
and, when only $\sigma_{r\theta}(r)$ is nonzero, the θ-force per unit volume is $\frac1{r^2}\frac{\partial}{\partial r}(r^2\sigma_{r\theta})$
(Appendix B; `ch05.polar_viscous_stress` derives the first with sympy through `core.curvilinear`).

*Gloss — the polar vector Laplacian (used in D03 step 8):* for a purely swirling flow $\mathbf u=u_\theta(r)\,\mathbf e_\theta$,
$(\nabla^2\mathbf u)_\theta=u_\theta''+\frac{u_\theta'}{r}-\frac{u_\theta}{r^2}$. The first two terms are the Laplacian of the
number $u_\theta(r)$; the extra $-u_\theta/r^2$ appears because the unit vector $\mathbf e_\theta$ itself turns from point
to point (Ch. 4 P124). The cell checks it on the two basic vortices: both give zero.
""")
nb.code(r"""
r_, G_, w_ = sp.symbols('r Gamma omega', positive=True)              # radius, circulation, vorticity
lap_theta = lambda ut: sp.simplify(sp.diff(ut, r_, 2) + sp.diff(ut, r_)/r_ - ut/r_**2)   # (∇²u)_θ for u = u_θ(r) e_θ
print(lap_theta(G_/(2*sp.pi*r_)), lap_theta(w_*r_/2))                # line vortex and solid body: 0 0 (no net viscous force)
""")
D("D03", ref="5.7")
note("N07", "The funnel", r"""
The same moves as D02 with $u_\theta=\Gamma/2\pi r$ ($\int\rho\Gamma^2/(4\pi^2r^3)dr=-\rho\Gamma^2/(8\pi^2r^2)$, and
$p\to p_\infty$ far away at z = 0) give (5.7). Isobars are $z=-\frac{\Gamma^2}{8\pi^2r^2g}-\frac{p-p_\infty}{\rho g}$, and
$\tfrac12u_\theta^2+gz+p/\rho$ is the same everywhere (irrotational). ⚠️ The book calls them "hyperboloids of
revolution of the second degree"; the formula says $(c-z)r^2=$ const, a cubic surface — read the shape from the
formula. **Number:** Γ = 1 m²/s, water, r = 0.1 m → deficit 1266.5 Pa and the free surface 12.9 cm below its far
level. The ideal funnel has no bottom; a real core (Rankine) gives it one: a tornado.
""", equation=EQ["5.7"], ref="5.7")
nb.code(r"""
print(ch05.line_vortex_pressure(0.1, 0.0, 1.0), ch05.isobar_height(0.1, 1.0, kind="line"))   # (5.7) at r = 0.1 m [Pa]; surface depth [m]
print(ch05.line_vortex_viscous_stress(0.1, 1.0, 1e-3))              # σ_rθ = −μΓ/πr² at r = 0.1 m, water [Pa]
for kind, prm in [("solid", dict(omega=10.0)), ("line", dict(Gamma=1.0)), ("gaussian", dict(Gamma=1.0, sigma=0.05))]:
    f = ch05.vortex_stress_force(kind, 0.1, mu=1e-3, **prm)          # σ_rθ and the net viscous force three ways, at r = 0.1 m
    print(f"{kind:9s} σ_rθ {f['sigma_rtheta']:+.6f} Pa | force: divergence {f['force_divergence']:+.2e}, "
          f"Laplacian {f['force_laplacian']:+.2e}, −μ∇×ω {f['force_curl']:+.2e} N/m³")
print(ch05.rankine_pressure(np.array([0.0, 0.05, 0.1, 0.2]), 0.0, 1.0, 0.1))   # Rankine (a = 0.1 m): p − p_∞ at r = 0, 5, 10, 20 cm [Pa]
""", explain=r"""
1. $p-p_\infty=-\frac{\rho\Gamma^2}{8\pi^2r^2}-\rho gz$ (5.7) at r = 0.1 m: −1266.515 Pa, and its free surface 0.129104 m
   below the far level.
2. The stress of note N06: −0.031831 Pa.
3. Note N09's table in numbers — the net viscous force three ways (D03's routes): the solid body has no stress and no
   force; the line vortex has stress −0.0318 Pa but a force of only ~10⁻⁸ N/m³ (stencil round-off; exactly 0 in the
   −μ∇×ω route); the Gaussian vortex (σ_c = 0.05 m) has both, the three routes agreeing at −0.2342 N/m³.
4. The Rankine composite: continuous pressure, twice as deep at the centre (−2533 Pa) as at the core edge (−1266.5 Pa).
""")
remind([
    ("torque, moment arm and moment of inertia", "torque = force × lever arm; I = Σ m r² resists spinning (Ch. 4 P116)."),
    ("power of a force", "power = force · velocity; for a torque, torque × angular speed (Ch. 4 P127)."),
    ("scipy.integrate.quad and dblquad", "adaptive 1-D (and 2-D) quadrature to near machine precision (Ch. 3 P87)."),
])
note("N08", "Torque to infinity", r"""
Round a spinning cylinder of radius a, the viscous torque per unit length transmitted across every circle of radius r
is the same (below). The power it carries, |T′| × (u_θ/r) = μΓ²/πr², falls outward; the difference between two radii is
dissipated as heat: $\int_a^R2\pi r\,\rho\varepsilon\,dr=\frac{\mu\Gamma^2}{\pi}\big(\frac1{a^2}-\frac1{R^2}\big)$. As R → ∞ all
the power the cylinder puts in is dissipated — and a container at any radius would feel the full torque. ⚠️ ε here is
the dissipation rate $\varepsilon=2\nu\big(S_{ij}-\tfrac13S_{mm}\delta_{ij}\big)^2+\frac{\mu_v}{\rho}S_{mm}^2\ge0$ *(4.58)*, not
the Levi-Civita symbol. **Number:** a = 0.1 m, R = 1 m, Γ = 1 m²/s, water: in 0.031831 W/m, out 0.000318 W/m,
dissipated 0.031513 W/m.
""", equation=r"T'=2\pi r^2\sigma_{r\theta}=-2\mu\Gamma")
nb.code(r"""
print([ch05.torque_per_length(r, 1.0, 1e-3) for r in (0.1, 1.0, 10.0)])   # torque per metre at three radii [N m/m]: −0.002 each
d = ch05.dissipation_outside_cylinder(0.1, 1.0, 1.0, 1e-3)           # a = 0.1 m, R = 1 m, Γ = 1 m²/s, water
print({k: f"{d[k]:.7f}" for k in ("power_in", "power_out", "dissipated")}, "W/m")   # the energy balance
print(f"in − out − dissipated = {d['power_in'] - d['power_out'] - d['dissipated']:.1e} W/m")   # closes to round-off
""")
note("N09", "The principle", r"""
Irrotationality does not mean *no viscous stress*; it means *no net viscous force*: for incompressible flow
$\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega\ \ (\nabla\cdot\mathbf u=0)$ *(4.40)* vanishes where ω = 0.
Solid-body rotation is the only motion with no viscous stress at all.

| vortex | σ_rθ | net viscous force |
|---|---|---|
| solid body | 0 | 0 |
| line vortex | −μΓ/πr² ≠ 0 | 0 |
| Gaussian (Lamb–Oseen) | ≠ 0 | ≠ 0 (it decays, Ch. 8) |
""")
nb.worked_example("a spinning bucket, a bathtub drain and a tornado", r"""
1. **Bucket**: radius 10 cm, spun at 5 rad/s, so ω = 10 s⁻¹ (⚠️ not 5). Rim–centre height of the free surface:
   $\omega^2R^2/8g=100\times0.01/(8\times9.81)=0.0127$ m = 12.7 mm. Volume is conserved, so the centre drops 6.4 mm and the
   rim rises 6.4 mm.
2. **Drain**: Γ = 1 m²/s (strong). At r = 10 cm: $u_\theta=\Gamma/2\pi r$ = 1.59 m/s; pressure deficit
   $\rho\Gamma^2/8\pi^2r^2=1000/(8\times9.87\times0.01)$ = 1266 Pa; the surface there lies 1266/(1000 × 9.81) = 0.129 m
   below the far level. At r = 5 cm the depth quadruples to 0.52 m — a funnel.
3. **Tornado** (air, ρ = 1.2 kg/m³, Γ = 10⁴ m²/s, core a = 50 m): peak wind Γ/2πa = 31.8 m/s, pressure deficit
   608 Pa at the core edge and 1216 Pa at the centre (Rankine) — about 1 % of the atmospheric pressure.
""")
nb.code(r"""
tor = ch05.vortex_pressure_scenario("rankine", 50.0, rho=1.2, Gamma=1e4, a=50.0)   # the tornado at its core edge, r = a = 50 m
print(f"u_max = {tor['u_theta']:.3f} m/s, p − p_∞ at the edge = {tor['p']:.2f} Pa")   # 31.831 m/s, −607.93 Pa
print(f"centre: {ch05.rankine_pressure(0.0, 0.0, 1e4, 50.0, rho=1.2):.2f} Pa")      # −1215.85 Pa: twice the edge deficit
print(tor["status"])                                                  # at r = a the stress jumps: an edge (line) force, no force per volume
print(f"edge line force {tor['edge_line_force']:.3e} N/m², force per volume {tor['net_viscous_force']}")   # NaN on purpose
""", explain=r"""
`vortex_pressure_scenario` gathers the tornado's numbers at the core edge. At exactly r = a the Rankine vorticity jumps,
so the viscous force is concentrated on the edge itself: the function reports that edge line force (the jump of σ_rθ,
−μΓ/πa²) and returns NaN for a force per unit volume, which does not exist there.
""")
remind([
    ("slider_figure", "`slider_figure(fn, name, values, …)` precomputes one set of curves per slider position, so the slider works on the web page (Ch. 1 P17)."),
])
nb.plotly(r"""
r = np.linspace(-0.15, 0.15, 301)                                    # radius across the tank (both sides of the axis) [m]
ra = np.where(np.abs(r) > 4e-3, np.abs(r), np.nan)                   # |r|, blank within 4 mm of the axis (the ideal funnel is infinite there)
g = 9.81                                                             # gravity [m/s²]
funnel = ch05.isobar_height(ra, 0.05, kind="line")                   # free surface of (5.7), Γ = 0.05 m²/s [m]
rank = ch05.rankine_pressure(np.abs(r), 0.0, 0.05, 0.01)/(1000*g)    # Rankine surface, core a = 1 cm: z = p(r, 0)/ρg [m]


def tank(Om):                                                        # Ω_tank [rad/s] → the curves at that rate
    w = 2*Om                                                         # vorticity ω = 2Ω_tank [1/s]
    surf = w**2*r**2/(8*g)                                           # free surface of (5.6), vertex at z = 0 [m]
    return {"free surface (5.6)": (r, surf),                         # the bowl
            "isobar p = 500 Pa (5.6)": (r, surf - 500/(1000*g)),      # constant-pressure surfaces lie under it, same shape
            "isobar p = 1000 Pa (5.6)": (r, surf - 1000/(1000*g)),
            "funnel (5.7), Γ = 0.05 m²/s": (r, funnel),              # fixed: the line vortex's funnel
            "Rankine surface, a = 1 cm": (r, rank)}                  # fixed: the funnel with a finite bottom


fig = slider_figure(tank, "Ω_tank", np.linspace(0, 10, 21 if not FAST else 11), unit="rad/s",       # precompute every slider position
                    xlabel="r [m]", ylabel="z [m]", title="A spinning bucket is a bowl; a drain is a funnel",  # axis labels with units
                    yrange=[-0.12, 0.06])
recolor(fig, {"free surface (5.6)": COLORS["teal"], "isobar p = 500 Pa (5.6)": COLORS["orange"],    # house colours by trace name
              "isobar p = 1000 Pa (5.6)": COLORS["orange"], "funnel (5.7), Γ = 0.05 m²/s": COLORS["rose"],
              "Rankine surface, a = 1 cm": COLORS["rose"]},
        dashes={"isobar p = 500 Pa (5.6)": "dash", "isobar p = 1000 Pa (5.6)": "dot", "Rankine surface, a = 1 cm": "dash"})  # dashed reference curves
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="(Notes `N46`, `N47`: our Figs. 5.2–5.3 analogues.) The teal bowl (and the orange isobars under it) deepens as the slider moves; the rose funnel plunges near the "
        "axis, and the dashed Rankine curve gives it a rounded bottom.",
    read=r"Both shapes come from $-\rho u_\theta^2/r=-\partial p/\partial r$ (5.5a): the bowl because $u_\theta$ grows with r, the "
         "funnel because it grows toward the axis. Isobars are parallel copies of the free surface in both.",
    change="…the tank rate doubled: the bowl's depth quadruples (ω² in (5.6)); its shape stays a paraboloid.")
nb.figure(r"""
r = np.linspace(0.01, 0.3, 120)                                      # radii [m]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8))                  # speed profiles | stress and force
a.plot(r, ch05.solid_body_from_vorticity(r, 10.0)[0], color=COLORS["teal"], label="solid body, ω = 10 1/s")      # (5.1)
a.plot(r, ch05.line_vortex_gamma(r, 1.0)[0], color=COLORS["teal"], ls="--", label="line vortex, Γ = 1 m²/s")     # (5.2)
a.plot(r, vortices.rankine_vortex(r, 1.0, 0.1)[0], color=COLORS["ink"], ls=":", label="Rankine, Γ = 1 m²/s, a = 0.1 m")   # (3.28)
a.set_ylim(0, 3); a.set_xlabel("r [m]"); a.set_ylabel(r"$u_\theta$ [m/s]"); a.legend(fontsize=8); a.set_title("(a) speed", fontsize=10)
line = [ch05.vortex_stress_force("line", ri, mu=1e-3, Gamma=1.0) for ri in r]            # line vortex, water
gaus = [ch05.vortex_stress_force("gaussian", ri, mu=1e-3, Gamma=1.0, sigma=0.05) for ri in r]   # Gaussian, σ_c = 0.05 m
b.plot(r, [q["sigma_rtheta"] for q in line], color=COLORS["rose"], label=r"$\sigma_{r\theta}$, line vortex")        # stress [Pa]
b.plot(r, [q["sigma_rtheta"] for q in gaus], color=COLORS["rose"], ls="--", label=r"$\sigma_{r\theta}$, Gaussian")   # stress [Pa]
b.set_xlabel("r [m]"); b.set_ylabel(r"viscous stress $\sigma_{r\theta}$ [Pa]"); b.set_ylim(-0.12, 0.01)
b2 = b.twinx()                                                       # a second y-axis for the force per volume
b2.plot(r, [q["force_curl"] for q in line], color=COLORS["ink"], label="net force, line vortex")                     # zero
b2.plot(r, [q["force_curl"] for q in gaus], color=COLORS["ink"], ls="--", label="net force, Gaussian")              # nonzero
b2.set_ylabel(r"net viscous force $f_\theta$ [N/m³]")                # axis labels with units
b.legend(loc="lower right", fontsize=8); b2.legend(loc="center right", fontsize=8); b.set_title("(b) stress vs net force", fontsize=10)
savefig(fig, "ch05", "stress_vs_force"); plt.show()                  # save the PNG to outputs/ch05, then draw
""", see="Left: the solid body's speed grows linearly, the line vortex's falls as 1/r, the Rankine curve follows the "
         "first inside 0.1 m and the second outside. Right: the line vortex's stress (solid rose) is largest near the "
         "axis while its net force (solid black) is zero; the Gaussian (dashed) has both stress and force.",
    read="Stress lives on one face of an element; force is the imbalance between faces. The line vortex is sheared "
         "everywhere but the faces balance; the Gaussian vortex has vorticity, so its faces do not balance and it decays.",
    change="…μ doubled: every stress and force doubles; the zeros stay zeros.")
explainer("vortex_pressure_funnel", "Why is a spinning bucket a bowl and a drain a funnel?",
          "the bowl, the funnel, the Rankine tornado and the rotating cylinder share one radial balance; dragging the "
          "rotation, the circulation and the core size while a probe element shows its pressure, stress and net force "
          "makes \"stress without force\" something you watch rather than read.",
          "", [
              "Preset 'bucket 5 rad/s': read the rim–centre height, 12.7 mm, and the zero stress bar.",
              "Switch to the line vortex: the stress bar fills, the net-force bar stays at zero.",
              "Try 'tornado': the same formula in air, 608 Pa at the core edge.",
              "Open Derivation → D03 step 6 and watch the two faces of the probe element.",
          ])
whatif(r"""
…the fluid were spinning but the frame of reference turned with it? The water would be at rest in that frame, yet the
free surface would still be a paraboloid: in the rotating frame there is no Coriolis force (nothing moves relative to
the frame), and the centrifugal force folds into an effective gravity (Ch. 4 §4.7); the free surface is a surface of
constant effective geopotential — the same reason the Earth's equator bulges. (Geostrophic balance, Ch. 13, is different:
Coriolis against a pressure gradient, which needs motion relative to the frame.)
""")

# =====================================================================================================================
# A.2 §5.2 Kelvin's Circulation Theorem — R09, C03, C04
# =====================================================================================================================
nb.section("5.2", "Kelvin's Circulation Theorem", intro=r"""
**What is this section about?** Circulation round a loop that moves with the fluid is conserved — if the fluid is
inviscid, barotropic, pushed only by conservative forces and watched from a non-rotating frame. The proof shows exactly
which term wakes up when each hypothesis fails: those are the only ways vorticity can be created. The last of them,
baroclinicity, has a simple mechanical picture: a torque.
""")
nb.recap("R09", "Barotropic fluid and the pressure function", r"""
A fluid is barotropic when its density depends on pressure alone, ρ = ρ(p) (constant density, isothermal or isentropic
gas). Then dp/ρ is the exact differential of the pressure function $\mathcal P(p)=\int_{p_o}^{p}dp'/\rho(p')$, as in
$\frac1\rho\frac{\partial p}{\partial x_j}=\frac{\partial}{\partial x_j}\int_{p_o}^{p}\frac{dp'}{\rho(p')}$ *(4.67)* — just as
$\mathbf g\cdot d\mathbf x=-d\Phi$ is exact for a conservative body force $\mathbf g=-\nabla\Phi$ *(4.18)*.
""", where="Ch. 4 §4.9")
nb.code(r"""
print(bernoulli.pressure_function(2.0e5, 1.0e5, kind="isothermal", T=288.15))   # 𝒫 from 1 to 2 bar, isothermal air at 15 °C [J/kg]
print(287.058*288.15*np.log(2))                                       # the same by hand, R = 287.058 J/(kg K): R T ln(p/p_o) = 57 334 J/kg
""")
remind([
    ("conservative body force and potential Φ (4.18)", r"a force per unit mass that is the gradient of a potential, $\mathbf g=-\nabla\Phi$ (Φ = gz for gravity), does no net work round a closed loop (Ch. 4 P115)."),
    ("rotating frame of reference", "a frame turning at Ω adds Coriolis and centrifugal forces to what an observer measures (Ch. 3 P103, Ch. 4 §4.7)."),
])

# ---- C03 ----------------------------------------------------------------------------------------------------------
core("C03", "Kelvin's circulation theorem", r"""
Put a ring of dye into water and let the flow stretch it seven times longer and wind it into a spiral. What has stayed
exactly the same?
""", eqs=("5.8",))
problem(r"""
Circulation measures the net swirl inside a loop. If we could prove it never changes for loops that move with the
fluid, we would know that a fluid starting from rest (zero circulation everywhere) can never start swirling on its own —
the basis of all "irrotational" flow in Ch. 6, of how wings get lift, and, on a rotating planet, of how cyclones spin up
(C11). The proof also tells us precisely how swirl *can* be made.
""")
idea("""
Gamma(t) = loop integral of u.dx round a loop of dye    d/dt:  (acceleration along the loop) + (stretching of the loop)
                                                               = loop integral of forces.dx    = loop integral of d(u^2/2) = 0 always
forces: pressure  -grad(p)/rho  ->  -loop integral of dp/rho  = 0 if rho = rho(p)   (barotropic)
        gravity   -grad(Phi)    ->  -loop integral of dPhi    = 0 if Phi exists     (conservative)
        viscous                 ->   loop integral ...        = 0 if inviscid       (ideal fluid)
        Coriolis                ->   not 0 in a rotating frame                       (inertial frame needed)
""")
remind([
    ("np.arctan2 and np.unwrap", "`np.arctan2(y, x)` is the polar angle of (x, y) in (−π, π]; `np.unwrap` removes its 2π jumps so the angle can be followed continuously (Ch. 2 P70)."),
])
P("P135", "closed-loop integral of an exact differential", r"""
If F is a single-valued function, going once round a closed loop brings F back to its starting value, so ∮dF = 0. It
fails for a multi-valued "function" such as the polar angle θ, which grows by 2π on a loop round the origin. Kelvin's
proof uses this twice: for ½u² and for the pressure function.
""", code=r"""
s = np.linspace(0, 2*np.pi, 2001)                  # a loop: the unit circle, first point = last point
x, y = np.cos(s), np.sin(s)                        # its coordinates [m]
F = x**2*y + 3*x                                   # a single-valued F(x, y)
print(np.sum(np.diff(F)))                          # ∮dF = F(end) − F(start) = 0 (round-off)
theta = np.unwrap(np.arctan2(y, x))                # the polar angle, followed continuously
print(theta[-1] - theta[0])                        # 6.283185 = 2π: θ is not single-valued round the origin
""")
P("P136", "periodic trapezoid rule on a closed loop", r"""
For a smooth closed loop parametrised by s ∈ [0, 1) with N equally spaced points, the plain average Σf(s_k)/N is
extraordinarily accurate — the errors of the trapezoid rule cancel round a periodic curve (the error falls faster than
any power of 1/N; for a trigonometric polynomial it is exact once N exceeds its highest frequency). We compute ∮u·dx this
way, with dx/ds from an FFT derivative (numpy's `np.fft` splits the periodic coordinates into waves $e^{2\pi iks}$;
differentiating a wave just multiplies it by 2πik — C07's primer says more).
""", code=r"""
for N in (4, 8, 16):                                        # three loop resolutions
    s = np.arange(N)/N                                      # equally spaced tags, end point not repeated
    x, y = np.cos(2*np.pi*s), np.sin(2*np.pi*s)              # the unit circle [m]
    dxds, dyds = -2*np.pi*y, 2*np.pi*x                        # the tangent dx/ds [m]
    u, v = -y**3, x                                          # a velocity field [m/s]; Stokes gives ∮u·dx = ∫(1 + 3y²)dA = 7π/4
    print(N, np.mean(u*dxds + v*dyds) - 7*np.pi/4)           # error: 0.785 at N = 4 (aliasing), round-off from N = 8
""")
gloss(["FFT derivative of a periodic sequence (dx/ds)"])
remind([
    ("solve_ivp options (rtol, t_eval, DOP853)", "`solve_ivp(f, (t0, t1), y0, t_eval=…, rtol=…, method=\"DOP853\")` integrates dy/dt = f(t, y) with an adaptive high-order method (Ch. 3 P94)."),
])
P("P137", "advecting many points in one solve_ivp call", r"""
To move a whole loop of N particles we give `solve_ivp` one long state vector (all x's, then all y's) and a right-hand
side that reshapes it to (2, N), evaluates the velocity at every point at once (vectorised) and flattens it back. One
adaptive integrator then moves the whole loop consistently.
""", code=r"""
from scipy.integrate import solve_ivp                                # the adaptive ODE solver
X0 = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, 0.0]])                   # three particles, shape (2, 3) [m]
rhs = lambda t, y: np.array([-y.reshape(2, -1)[1], y.reshape(2, -1)[0]]).ravel()   # solid body u = (−y, x), flattened [m/s]
sol = solve_ivp(rhs, (0, np.pi/2), X0.ravel(), rtol=1e-10)           # a quarter turn (1 rad/s for π/2 s)
print(sol.y[:, -1].reshape(2, -1).round(6))                          # each point turned by 90°: [[0, −1, 0], [1, 0, −1]]
""")
nb.md(r"""
*Gloss — material label.* Tag each particle of the loop once with a number s ∈ [0, 1); the tags travel with the
particles, so the loop is always "the particles with tags 0…1" — fixed limits, even though the loop moves (Ch. 3 P89).
""")
remind([
    ("functions of time with parameters", "a quantity like X(s, t) is a function of time for each fixed tag s (Ch. 3 P89)."),
    ("differentiation under the integral sign", r"$\frac{d}{dt}\int_a^bF\,ds=\int_a^b\frac{\partial F}{\partial t}ds$ when the limits a, b are fixed (Ch. 3 P109)."),
    ("material line element D(δx)/Dt = δu", "a short line of particles stretches and turns at the velocity difference across it (Ch. 3 P99)."),
    ("multivariable first-order Taylor expansion", r"$\mathbf u(\mathbf x+d\mathbf x)\approx\mathbf u(\mathbf x)+(d\mathbf x\cdot\nabla)\mathbf u$ for a short step (Ch. 3 P98)."),
    ("chain rule for the kinetic energy u·du = d(½u²)", r"$u_i\,du_i=d(\tfrac12u_iu_i)$ (Ch. 4 P128)."),
], lead="needed in the derivation below")
D("D04", ref="5.9")
note("N10", "Two causes of change", r"""
D04's step 4 is (5.9): the fluid accelerating along the loop, and the loop's elements changing.
""", equation=EQ["5.9"], ref="5.9")
note("N11", "The element's rate is the velocity difference", r"""
(Fig. 5.4): $\mathbf u+d\mathbf u=\frac{D}{Dt}(\mathbf x+d\mathbf x)$, so $D(d\mathbf x)/Dt=d\mathbf u$, and the second integral of
(5.9) is $\oint_Cu_i\,du_i=\oint_Cd(\tfrac12u_i^2)=0$.
""", equation=r"\frac{D}{Dt}(d\mathbf x)=d\mathbf u")
remind([
    ("multivariable chain rule along a path", r"the change of p along a small step is $dp=\nabla p\cdot d\mathbf x$ (Ch. 3 P91)."),
], lead="needed in the derivation below")
D("D05", ref="5.11")
note("N12", "The forces round the loop", r"D05's step 4 is (5.10):", equation=EQ["5.10"], ref="5.10")
note("N13", "What survives", r"""
is (5.11): Kelvin holds if the fluid is inviscid or the net viscous force vanishes along C. **Number**, the
Lamb–Oseen vortex (Γ₀ = 0.01 m²/s, ν = 10⁻⁶ m²/s) on the circle r = 5 mm at t = 10 s:
$\Gamma=\Gamma_0(1-e^{-r^2/4\nu t})$ = 0.00465 m²/s and falling at 3.35×10⁻⁴ m²/s² — exactly the (5.11) integral.
⚠️ The book adds "this occurs when C lies entirely in irrotational fluid" — true for incompressible, constant-μ flow,
where the viscous force is $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega$ (4.40); a compressible flow keeps a
$\nabla(\nabla\cdot\mathbf u)$ force.
""", equation=EQ["5.11"], ref="5.11")
nb.worked_example("circulation that survives a stretching — and one that leaks away", r"""
1. **Cellular flow** $\psi=\sin x\sin y$ ($u=\sin x\cos y$, $v=-\cos x\sin y$, in m and m/s): a steady inviscid flow whose
   vorticity $\omega=2\sin x\sin y$ is constant on streamlines. A circle of radius 0.5 m round (π/2, 0.9 m) has
   $\Gamma(0)=\iint2\sin x\sin y\,dA$ = 1.155 m²/s. Kelvin predicts Γ(t) = 1.155 m²/s forever.
2. **Lamb–Oseen** (viscous): the circle r = 5 mm with Γ₀ = 0.01 m²/s, ν = 10⁻⁶ m²/s at t = 10 s:
   $r^2/4\nu t=25\times10^{-6}/(4\times10^{-5})$ = 0.625; Γ = 0.01(1 − e^{−0.625}) = 0.01 × 0.4647 = 0.00465 m²/s; at
   t = 20 s: $r^2/4\nu t$ = 0.3125, Γ = 0.00268 m²/s — the vorticity diffuses out through the loop.
""")
nb.code(r"""
sc = ch05.kelvin_scenario("cellular")                                # the cellular flow and its initial circle (2048 labelled points)
t = np.linspace(0, sc["t_end"], 7)                                   # 0 … 37.7 s: six turnover times
G = ch05.material_circulation(sc["u"], sc["pts0"], t)                # Γ(t) of the material loop (one solve_ivp call) [m²/s]
L = [ch05.loop_length(Pk) for Pk in ch05.material_loop(sc["u"], sc["pts0"], t)]   # the loop's length at those times [m]
print("Γ(t) =", [f"{g:.6f}" for g in G], " spread", f"{np.ptp(G):.1e} m²/s")   # 1.155130 every time
print("length(t) =", np.round(L, 3))                                 # 3.142 m → 21.0 m: stretched ×6.7
print(ch05.kelvin_rate_terms(sc["u"], sc["pts0"], 0.0))              # the two terms of (5.9) and their sum [m²/s²]
Gc = ch05.material_circulation(sc["u"], sc["pts0"], [0.0, 0.01, 0.02])   # Γ at three close times
print(f"dΓ/dt by central difference: {(Gc[2] - Gc[0])/0.02:.1e} m²/s²")    # ≈ 0, the same as the acceleration term
print(ch05.lamb_oseen_circulation(0.005, 10.0, 0.01, 1e-6))          # Lamb–Oseen (Γ(r, t), ∂Γ/∂t) at r = 5 mm, t = 10 s
print(ch05.lamb_oseen_viscous_loop_integral(0.005, 10.0, 0.01, 1e-6))   # the (5.11) loop integral of the viscous force [m²/s²]
print(ch05.kelvin_scenario_circulation("rotating", 5.0))              # a loop in a rotating frame: Γ has grown to 1.99 m²/s (C11)
""", explain=r"""
1. The scenario's velocity and its initial loop (2048 labelled points).
2. The loop advected by `solve_ivp` in one call (primer P137) and its circulation by the periodic trapezoid rule (P136):
   Γ = 1.155130 m²/s at every time (spread 3×10⁻¹³) while the length grows from 3.142 m to 20.997 m (×6.7).
3. The two terms of (5.9): the acceleration term ≈ 4×10⁻¹³ and the contour term ≈ −3×10⁻¹⁴ m²/s² — and dΓ/dt from
   central differences of Γ(t) is ≈ 0 as well.
4. Lamb–Oseen: Γ = 0.004647 m²/s falling at −3.34538×10⁻⁴ m²/s², and the viscous loop integral of (5.11) is the same
   number.
5. In a rotating frame the relative circulation grows (to 1.99 m²/s after 5 s here) — C11 explains why.
""")
remind([
    ("RK4 by hand", "four velocity evaluations per step combined with weights 1, 2, 2, 1 — fourth-order accurate (Ch. 3 P95)."),
])
nb.md("**From scratch — our own particles and our own loop sum:** RK4 by hand on the labelled points, then the periodic trapezoid rule with an FFT tangent.")
nb.check_agree(r"""
u = sc["u"]                                                          # the cellular flow u(x, t) [m/s]
P = sc["pts0"][:, ::(1 if not FAST else 2)].copy()                   # the labelled loop points (every 2nd with FAST) [m]
dt = 0.01                                                            # RK4 step [s]
for k in range(600):                                                 # 600 steps → t = 6 s
    k1 = u(P, 0); k2 = u(P + 0.5*dt*k1, 0); k3 = u(P + 0.5*dt*k2, 0); k4 = u(P + dt*k3, 0)   # the four RK4 slopes (steady flow)
    P = P + dt/6*(k1 + 2*k2 + 2*k3 + k4)                             # advance every point together
N = P.shape[1]                                                       # number of points
kw = np.fft.fftfreq(N, 1.0/N)                                        # integer wavenumbers of the periodic tag s ∈ [0, 1)
dPds = np.real(np.fft.ifft(2j*np.pi*kw*np.fft.fft(P, axis=1), axis=1))   # dx/ds by FFT (each wave × 2πik)
Gamma = np.mean(np.sum(u(P, 0)*dPds, axis=0))                        # ∮u·dx = mean of u·dx/ds over the tags [m²/s]
print(f"ours {Gamma:.10f}, library {ch05.material_circulation(sc['u'], sc['pts0'], [0, 6.0])[-1]:.10f} m²/s")
assert np.allclose(Gamma, ch05.material_circulation(sc["u"], sc["pts0"], [0, 6.0])[-1], rtol=1e-7)   # same numbers
assert abs(Gamma - 1.155130) < 1e-6                                  # and the circulation did not move
""")
remind([
    ("animate and show_animation", "`animate(update, frames, fig)` builds a movie from an update function; `show_animation` plays it (Ch. 1 P16)."),
])
nb.md(r"""
**The loop and its circulation in motion (Fig. 5.4 analogue, note `N48`).** Left: the cellular flow (grey shading = its stream function, whose level curves are the streamlines) and the
material loop (teal), with one element dx (black arrow) and the velocities u (teal) and u + du (orange) at its ends.
Right: Γ(t)/Γ(0) of this loop (teal) and of the viscous Lamb–Oseen circle of the tiny example (rose).
""")
nb.animation(r"""
nf = 30 if not FAST else 16                                          # number of frames
tf = np.linspace(0, sc["t_end"], nf)                                 # frame times 0 … 37.7 s
Ploop = ch05.material_loop(sc["u"], sc["pts0"], tf)                  # the loop at every frame, shape (nf, 2, 2048) [m]
Gf = np.array([ch05.loop_circulation(sc["u"], Pk) for Pk in Ploop])  # its circulation [m²/s]
GLO = ch05.lamb_oseen_circulation(0.005, tf + 10.0, 0.01, 1e-6)[0]   # Lamb–Oseen circle r = 5 mm, vortex 10 s old at t = 0
xg = np.linspace(0, np.pi, 120)                                      # grid for the streamlines [m]
XG, YG = np.meshgrid(xg, xg)                                         # (x, y) grid
with plt.rc_context(ANIM_RC):                                         # light settings for fast frames
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 4.0), dpi=80, gridspec_kw=dict(width_ratios=[1.1, 1]))
fig.subplots_adjust(left=0.07, right=0.98, bottom=0.13, top=0.9, wspace=0.25)   # fixed margins (no layout engine)
a.imshow(np.sin(XG)*np.sin(YG), origin="lower", extent=(0, np.pi, 0, np.pi), cmap="Greys", alpha=0.35)   # grey shading = ψ = sin x sin y: streamlines are its level curves
(loopline,) = a.plot([], [], color=COLORS["teal"], lw=1.4)           # the material loop
q = a.quiver([0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], color=[COLORS["ink"], COLORS["teal"], COLORS["orange"]],
             angles="xy", scale_units="xy", scale=1, width=0.006)    # dx, u, u + du arrows
a.set_xlim(0, np.pi); a.set_ylim(0, np.pi); a.set_aspect("equal"); a.set_xlabel("x [m]"); a.set_ylabel("y [m]")
b.plot(tf, GLO/GLO[0], color=COLORS["rose"], lw=2, label="Lamb–Oseen circle (viscous)")        # the viscous leak
b.axhline(1.0, color=COLORS["muted"], ls="--", lw=1, label="Γ(0)")    # ghost of the start value
(gline,) = b.plot([], [], color=COLORS["teal"], lw=2.5, label="cellular loop (inviscid)")      # Kelvin
b.set_xlim(0, tf[-1]); b.set_ylim(0, 1.15); b.set_xlabel("t [s]"); b.set_ylabel("Γ(t)/Γ(0) [–]"); b.legend(fontsize=8, loc="lower left")
txt = a.set_title("")                                                # running readout
k0, k1 = 0, 40                                                       # the element's two particles (tags 0 and 40/2048)


def update(i):                                                       # draw frame i
    Pk = Ploop[i]                                                    # loop points at this time
    loopline.set_data(np.r_[Pk[0], Pk[0, :1]], np.r_[Pk[1], Pk[1, :1]])   # close the loop for drawing
    x0, x1 = Pk[:, k0], Pk[:, k1]                                    # the element's tail and head
    u0, u1 = sc["u"](x0, 0.0), sc["u"](x1, 0.0)                      # velocities there [m/s]
    q.set_offsets(np.array([x0, x0, x1]))                            # arrow tails
    q.set_UVC([x1[0] - x0[0], 0.3*u0[0], 0.3*u1[0]], [x1[1] - x0[1], 0.3*u0[1], 0.3*u1[1]])   # dx, u, u + du (scaled)
    gline.set_data(tf[:i + 1], Gf[:i + 1]/Gf[0])                     # Γ(t)/Γ(0) so far
    txt.set_text(f"t = {tf[i]:.1f} s, loop length {ch05.loop_length(Pk):.1f} m, Γ = {Gf[i]:.6f} m²/s")
    return loopline, q, gline


with plt.rc_context(ANIM_RC):                                         # render the frames with the light settings
    show_animation(animate(update, frames=nf, fig=fig, interval=150), player="video")   # smooth movie
""")
see_read_change(
    see="The loop winds into a long thin spiral round the centre of the cell; its teal curve on the right does not move "
        "off 1, while the rose curve of the viscous vortex falls to about a quarter.",
    read=r"However long the loop gets, the acceleration along it integrates to zero (D05): Kelvin's $D\Gamma/Dt=0$ (5.8). "
         "The rose curve is the viscous leak of (5.11).",
    change=r"…ν switched on in the cellular flow: the teal curve would fall as $e^{-2\nu t}$, whatever the loop's shape — "
           r"for this flow $\nabla^2\mathbf u=-2\mathbf u$, so the viscous term of (5.11) is $\nu\oint\nabla^2\mathbf u\cdot d\mathbf x=-2\nu\Gamma$. "
           "The cell below checks it with two loops of different shape.")
nb.md(r"""
**Checking the what-if.** With viscosity, $\mathbf u=e^{-2\nu t}(\sin x\cos y,\ -\cos x\sin y)$ is an exact Navier–Stokes
solution (the Taylor–Green vortex, the viscous version of our cellular flow: $\nabla^2\mathbf u=-2\mathbf u$ and the
pressure balances the rest). For any material loop, (5.11) then gives $D\Gamma/Dt=\nu\oint\nabla^2\mathbf u\cdot d\mathbf x=-2\nu\Gamma$,
so $\Gamma(t)=\Gamma(0)e^{-2\nu t}$ — the same decay for a circle and for a square.
""")
nb.code(r"""
nu_c = 0.01                                                          # a large viscosity so the decay shows in 37.7 s [m²/s]
u_v = lambda x, t: np.exp(-2*nu_c*t)*ch05.cellular_flow(x)           # the decaying cellular (Taylor–Green) flow [m/s]
tv = np.array([0.0, 10.0, 20.0, sc["t_end"]])                        # four times [s]
loops = {"circle": sc["pts0"], "square": ch05.square_loop_points((np.pi/2, 0.9), 0.8, n=1024)}   # two loop shapes (2, N) [m]
for name, P0 in loops.items():
    Gv = ch05.material_circulation(u_v, P0, tv)                      # Γ(t) of the material loop [m²/s]
    print(f"{name:6s} Γ(t)/Γ(0) = {np.round(Gv/Gv[0], 6)}   e^(−2νt) = {np.round(np.exp(-2*nu_c*tv), 6)}")
""", explain=r"""
Both loops — a circle and a square, advected by the decaying flow — lose circulation at exactly the rate
$e^{-2\nu t}$ (0.8187, 0.6703, 0.4705 at 10, 20 and 37.7 s for ν = 0.01 m²/s): the viscous leak of (5.11) depends on
the flow, not on the loop's shape.
""")
note("N14", "Three ways to make or destroy vorticity", r"""
— one per broken hypothesis.

| source | term of (5.10) | example | number |
|---|---|---|---|
| non-conservative body force | $\oint\mathbf g\cdot d\mathbf x\ne0$ (Coriolis) | a drain vortex in a tank on the rotating Earth | the rotating scenario gains 1.99 m²/s in 5 s |
| non-barotropic (baroclinic) | $-\oint dp/\rho\ne0$ | lock exchange (C04), sea breeze | the lock-exchange square: 0.0467 m²/s² (cell below) |
| net viscous force | $\oint(1/\rho)\partial\sigma_{ij}/\partial x_j\,dx_i\ne0$ | boundary layers at walls (Ch. 9) | −3.35×10⁻⁴ m²/s² (Lamb–Oseen) |

The baroclinic number is checked two ways below: the loop integral $-\oint dp/\rho$ round a 0.2 m square straddling
the lock-exchange interface at the first instant, and — by Stokes — the area integral of the baroclinic source
$\frac1{\rho^2}\nabla\rho\times\nabla p$ over the square.
""")
nb.md(r"""
*Gloss — Gauss–Legendre nodes (used in the next cell).* `itg.gauss_legendre_nodes(lo, hi, n)` returns n points and n
weights on [lo, hi] chosen so that Σ wᵢ f(xᵢ) is exact for every polynomial f of degree up to 2n − 1; for an area we
take all pairs of nodes and multiply their weights. Primer P143 in C07 shows it with a 4-node demo.
""")
nb.code(r"""
bsc = ch05.kelvin_scenario("baroclinic")                             # fluid at rest in the lock-exchange field, square loop 0.2 m
rate = ch05.kelvin_scenario_rate("baroclinic")                       # the terms of (5.9)/(5.10) at t = 0 [m²/s²]
print({k: round(v, 10) for k, v in rate.items()})                    # only "pressure" (= −∮dp/ρ) is nonzero
lx = ch05.lock_exchange_fields(1000.0, 1025.0, 0.1)                  # the same density and hydrostatic pressure as functions
(x0, y0), (x1, y1) = bsc["pts0"].min(axis=1), bsc["pts0"].max(axis=1)   # the square's corners [m]
xs, wx = itg.gauss_legendre_nodes(x0, x1, 40); ys, wy = itg.gauss_legendre_nodes(y0, y1, 40)   # 40 × 40 Gauss nodes
XS, YS = np.meshgrid(xs, ys, indexing="ij")                          # all nodes
gr, gp, rr = lx["grad_rho"](XS, YS), lx["grad_p"](XS, YS), lx["rho"](XS, YS)   # ∇ρ, ∇p, ρ at the nodes
source = (gr[0]*gp[1] - gr[1]*gp[0])/rr**2                           # (∇ρ × ∇p)_z/ρ² [1/s²]
print(f"loop: {rate['pressure']:.12f}   area (Stokes): {np.sum(source*wx[:, None]*wy[None, :]):.12f} m²/s²")
""", explain=r"""
Only the pressure term survives at the first instant (the fluid is at rest: no viscous force, gravity is conservative,
the frame is inertial): $-\oint dp/\rho$ = 0.046708 m²/s², and the area integral of the baroclinic source over the square
gives the same number to 12 digits — the baroclinic term of (5.10) is exactly the flux of $\frac1{\rho^2}\nabla\rho\times\nabla p$.
""")
note("N16", "Four restrictions", r"""
keep irrotational flow irrotational: no net viscous force along C (boundary layers diffuse vorticity in), conservative
body forces (they act through the centre of mass), a barotropic fluid (else baroclinic: temperature, salinity,
composition), an inertial frame (§4.7's extra terms). The decision table below tries all 16 combinations.
""")
nb.code(r"""
import itertools                                                     # itertools.product lists every True/False combination
for combo in itertools.product([True, False], repeat=4):             # (inviscid, barotropic, conservative, inertial)
    r = ch05.kelvin_hypotheses(*combo)                               # which term of (5.10) survives
    print(combo, "Kelvin holds" if r["holds"] else "Kelvin fails:", r["surviving_terms"])
""", explain="""
Only (True, True, True, True) keeps Kelvin's theorem; each False adds exactly its own term ("viscous", "baroclinic",
"body", "coriolis") — the three sources of vorticity plus the frame.
""")
nb.md(r"""
> ⚠️ **Common confusion:** "circulation is conserved round any loop." Only round a **material** loop. A loop fixed in
> space in the same inviscid flow sees Γ change as different fluid (with different vorticity) passes through it — the
> explainer's "fixed loop" mode shows it.
""")
explainer("kelvin_material_loop", "Stretch and tangle a loop of dye — what stays the same?",
          "Kelvin's theorem is about time: watching the loop tangle while Γ(t) stays flat, then breaking one hypothesis "
          "at a time and seeing which term of (5.10) wakes up, teaches both the theorem and its limits.",
          "Its Helmholtz mode shows C05 (next section).", [
              "Press ▶ on the cellular flow: the loop grows ×6.7, Γ stays at 1.155 m²/s.",
              "Switch to 'fixed loop': Γ now changes — the loop is not material.",
              "Choose the baroclinic flow: the orange ∮dp/ρ curve lifts off zero.",
              "Try the rotating frame: Γ grows while Γ_a stays at π — a preview of C11.",
          ])
whatif(r"""
…the fluid were stratified (density depending on temperature as well as pressure)? Then $\oint dp/\rho$ need not vanish:
where surfaces of constant density cross surfaces of constant pressure, circulation — and vorticity — is born. C04
shows the torque that does it.
""")

# ---- C04 ----------------------------------------------------------------------------------------------------------
core("C04", "Baroclinic torque: crossed isobars and isopycnals spin fluid up", r"""
On a sunny afternoon the land warms, the sea stays cool, and a sea breeze starts to blow. Why does a density difference
*sideways* make air turn over, while a density difference up–down leaves it at rest?
""", eqs=("5.28",))
problem(r"""
Warm and cold air side by side, fresh river water meeting salty sea water, a front between two air masses: in each,
density changes horizontally while gravity pulls vertically. Kelvin's theorem says the barotropic hypothesis is broken
there, so circulation must change. But *why*, mechanically? A small blob of fluid gives the answer: pressure pushes
through its middle, its weight hangs off-centre, and it is spun.
""")
idea("""
barotropic element                    baroclinic element
isobars parallel to isopycnals        isobars cross isopycnals
pressure force through centre O       pressure force through centre O   (a circle: every push points at O)
centre of mass G on the force line    G shifted toward the heavy side, off the force line
no lever arm -> no torque             lever arm O-G -> torque ~ grad(rho) x grad(p)  ->  spin-up = grad(rho) x grad(p) / rho^2   (5.28)
""")
nb.md(r"""
*Glosses used below.* **Centre of mass** of a body with varying density: $\mathbf x_G=\int\rho\mathbf x\,dA/\int\rho\,dA$ —
pulled toward the heavy side. **Moment of inertia of a disc** about its centre, $I=\tfrac12MR^2$ (Ch. 4 P116 gave the
cube's). **Smoothed step**: $\rho(x)=\bar\rho-(\Delta\rho/2)\tanh(2x/\delta)$ goes from $\bar\rho+\Delta\rho/2$ to
$\bar\rho-\Delta\rho/2$ across a layer of width ≈ δ, with slope −Δρ/δ at its centre (`np.tanh`). **Over a disc,**
$\int x_ix_j\,dA=(\pi R^4/4)\delta_{ij}$. **Moving the reference point of a torque** (used in D06 step 6): the moment of a
force F about a point G equals its moment about O plus $(\mathbf r_O-\mathbf r_G)\times\mathbf F$ — the lever arm
changes by the vector from G to O.
""")
remind([
    ("net force from pressure and the gradient theorem", r"the pressure force on a body is $-\oint p\,\mathbf n\,ds=-\int\nabla p\,dA$ (Ch. 1 P28)."),
], lead="needed in the derivation below")
D("D06", ref="5.28")
nb.worked_example("a 1 cm blob in water with a sideways density gradient", r"""
Disc radius R = 1 cm; ρ₀ = 1000 kg/m³ increasing to the right at 10 kg/m⁴ (∇ρ = (10, 0) kg/m⁴); hydrostatic
∇p = (0, −ρ₀g) = (0, −9810) Pa/m.

1. Centre-of-mass offset: $R^2\lvert\nabla\rho\rvert/4\rho_0=10^{-4}\times10/4000$ = 2.5×10⁻⁷ m to the right.
2. Pressure force: $\pi R^2\lvert\nabla p\rvert=3.1416\times10^{-4}\times9810$ = 3.08 N/m upward, through the centre.
3. Torque about G: 2.5×10⁻⁷ × 3.08 = 7.70×10⁻⁷ N m/m, clockwise (the upward push acts left of G).
4. $I_G=\rho_0\pi R^4/2$ = 1.571×10⁻⁵ kg m.
5. Spin-up: 2 × 7.70×10⁻⁷/1.571×10⁻⁵ = 0.0981 s⁻², clockwise — the heavy right side sinks.
6. Formula: $(\nabla\rho\times\nabla p)_z/\rho_0^2=10\times(-9810)/10^6$ = −0.0981 s⁻² ✓ (negative = clockwise).
""")
remind([
    ("observed order of convergence (log–log slope)", "if an error falls like hᵖ, the slope of log(error) against log(h) is p (Ch. 1 P13)."),
])
nb.code(r"""
r = ch05.pressure_torque_on_element()                               # the example's disc: R = 1 cm, ∇ρ = (10, 0), ∇p = (0, −9810)
print(f"x_G = ({r['x_G'][0]:.4e}, {r['x_G'][1]:.1e}) m, torque {r['torque']:.4e} N m/m, I_G {r['I_G']:.5e} kg m")
print(f"spin-up by the torque {r['spin_up']:.10f}, formula {r['baroclinic']:.10f} 1/s², ratio {r['ratio']:.10f}")   # both routes: −0.0981
print(ch05.pressure_torque_on_element(grad_rho=(0.0, -10.0))["spin_up"])   # density increasing downward: ∇ρ ∥ ∇p → 0
R = np.array([0.1, 0.03, 0.01])                                      # three disc radii [m]
err = [abs(ch05.pressure_torque_on_element(radius=Ri, grad_rho=(100.0, 0.0))["ratio"] - 1) for Ri in R]   # torque route vs formula
print("|ratio − 1| =", [f"{e:.3e}" for e in err], f" observed order {observed_order(R, err):.2f}")  # the difference falls as R²
""", explain=r"""
1. The disc of the example: the torque route (pressure force about the centre of mass, divided by $I_G$, doubled) and
   the formula $(\nabla\rho\times\nabla p)_z/\rho^2$ both give −0.0981 s⁻² (ratio 1 + 1.2×10⁻⁹).
2. Density increasing downward (stable stratification, ∇ρ ∥ ∇p): no torque at all (≈ 6×10⁻¹⁴, round-off).
3. With a larger density gradient the two routes differ by (1.25×10⁻⁵, 1.1×10⁻⁶, 1.25×10⁻⁷) for R = 10, 3, 1 cm — the
   difference falls as R² (observed order 2.00): in the limit of a point element they are the same thing.
""")
nb.md("**From scratch — sum the pushes ourselves:** 4000 rim points for the pressure force and its torque, polar cells for the mass, the centre of mass and the moment of inertia.")
nb.check_agree(r"""
Rd, n = 0.01, 4000                                                   # disc radius [m], rim points
th = (np.arange(n) + 0.5)/n*2*np.pi                                  # rim angles [rad]
ds = 2*np.pi*Rd/n                                                    # rim element length [m]
xb, yb = Rd*np.cos(th), Rd*np.sin(th)                                # rim points [m]
p = 1e5 + 0*xb - 9810*yb                                             # linear hydrostatic pressure on the rim [Pa]
Fx, Fy = np.sum(-p*np.cos(th))*ds, np.sum(-p*np.sin(th))*ds           # net pressure force −∮p n ds [N/m]
nr, nt = (400, 800) if not FAST else (200, 400)                      # polar cells for the area integrals
rc = (np.arange(nr) + 0.5)/nr*Rd; tc = (np.arange(nt) + 0.5)/nt*2*np.pi   # cell centres
RC, TC = np.meshgrid(rc, tc, indexing="ij")                          # all cells
XC, YC = RC*np.cos(TC), RC*np.sin(TC)                                # their positions [m]
dA = RC*(Rd/nr)*(2*np.pi/nt)                                         # their areas [m²]
rho = 1000.0 + 10.0*XC                                               # ρ = ρ0 + ∇ρ·x [kg/m³]
M = np.sum(rho*dA)                                                   # mass per unit depth [kg/m]
xG, yG = np.sum(rho*XC*dA)/M, np.sum(rho*YC*dA)/M                    # centre of mass [m]
IG = np.sum(rho*((XC - xG)**2 + (YC - yG)**2)*dA)                    # moment of inertia about G [kg m]
tau = np.sum(((xb - xG)*(-p*np.sin(th)) - (yb - yG)*(-p*np.cos(th)))*ds)   # z-torque of the pushes about G [N m/m]
print(f"force ({Fx:.2e}, {Fy:.4f}) N/m, x_G = {xG:.3e} m, spin-up ours {2*tau/IG:.8f}, library {r['spin_up']:.8f} 1/s²")
assert np.allclose(2*tau/IG, r["spin_up"], rtol=1e-5)                 # same spin-up (midpoint cells: ~1e-6 relative)
""")
nb.figure(r"""
fig, (a, b) = plt.subplots(1, 2, figsize=(9, 4.2))                   # barotropic | baroclinic
sa = disc_element(a, (0, 0), 1.0, isobars=0.0, isopycnals=0.0)       # isobars and isopycnals both horizontal
sb = disc_element(b, (0, 0), 1.0, isobars=0.0, isopycnals=np.pi/4)   # isopycnals tilted 45° to the isobars
a.set_title(sa["status"], fontsize=9)                                # the scenario's verdict: no torque
b.set_title(f"baroclinic: spins {sb['sense']}", fontsize=9)          # the sense of the torque
fig.suptitle("orange = isobars, blue dashed = isopycnals, grey = pressure pushes, amber G = centre of mass (offset exaggerated)", fontsize=9)
savefig(fig, "ch05", "baroclinic_discs"); plt.show()                 # save the PNG to outputs/ch05, then draw
""", see="(Note `N50`, Fig. 5.6 analogue.) On the left, isobars and isopycnals are parallel and G sits on the vertical force line through the centre; on "
         "the right, the isopycnals are tilted 45°, G sits off the force line, and a rose arrow shows the element "
         "turning.",
    read="The pressure pushes (grey) all point at the geometric centre, so their sum passes through it; only a density "
         "gradient across the isobars moves G off that line — that lever arm is the torque. (G's true offset is "
         "R²∇ρ/4ρ₀, about 10⁻⁵ of R; it is drawn at 0.35R so you can see it.)",
    change="…the isopycnals turned to vertical (90° to the isobars): the offset across the force line, the lever arm and "
           "the torque are largest (sin 90° = 1).")
remind([
    ("np.deg2rad", "converts degrees to radians, which numpy's trigonometric functions expect (Ch. 2 P66)."),
])
nb.plotly(r"""
ang = np.linspace(0, 180, 181)                                       # angle between isopycnals and isobars [deg]
rates = np.array([ch05.baroclinic_element_scenario(np.deg2rad(a_), 10.0, 0.01)["spin_up"] for a_ in ang])   # spin-up [1/s²]


def tilt(a_deg):                                                     # slider value → curves
    rnow = ch05.baroclinic_element_scenario(np.deg2rad(a_deg), 10.0, 0.01)["spin_up"]   # the rate at this angle
    return {"spin-up rate": (ang, rates), "this angle": ([a_deg], [rnow]),                          # the curves for this slider value
            "formula |∇ρ||∇p| sin θ/ρ²": (ang, -10*9810*np.sin(np.deg2rad(ang))/1000.0**2)}         # the formula of (5.28) as a dashed ghost


fig = slider_figure(tilt, "tilt", np.linspace(0, 180, 19 if not FAST else 10), unit="°",            # precompute every slider position
                    xlabel="angle between isopycnals and isobars [deg]", ylabel="spin-up rate Dω/Dt [1/s²]",  # axis labels with units
                    title="The torque switches on as the isolines cross", modes={"this angle": "markers"})  # the message of the figure
recolor(fig, {"spin-up rate": COLORS["orange"], "this angle": COLORS["ink"], "formula |∇ρ||∇p| sin θ/ρ²": COLORS["muted"]},  # house colours by trace name
        dashes={"formula |∇ρ||∇p| sin θ/ρ²": "dash"})                                               # dashed reference curves
for tr in fig.data:                                                  # bigger marker for the current angle
    if tr.name == "this angle":                                                                     # the current-angle marker
        tr.marker.size = 12                                                                         # made larger
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="An orange curve on the grey dashed formula: zero at 0° and 180° (parallel isolines), largest in size at 90°; "
        "the black dot follows the slider.",
    read=r"The spin-up is $\frac1{\rho^2}\lvert\nabla\rho\rvert\lvert\nabla p\rvert\sin\theta$: the cross product of (5.28). Its "
         "sign here is negative (clockwise) because density increases to the right of the upward pressure force.",
    change="…∇ρ doubled: the whole curve doubles; the zeros at 0° and 180° stay — no density contrast, however large, "
           "spins fluid whose isopycnals are parallel to its isobars.")
remind([
    ("reduced gravity", r"$g'=g\,\Delta\rho/\bar\rho$ is the gravity felt by one fluid floating in another (Ch. 4 P131)."),
])
note("N15", "Lock exchange", r"""
(Fig. 5.5): fresh water (ρ₁) and salt water (ρ₂) side by side, a gate pulled out: the heavy fluid slumps under, the
light rides over, the interface tilts — vorticity has been made. Its first-instant rate (derived below; the book leaves
it to Exercise 5.5) is the formula underneath. **Number:** ρ₁ = 1000, ρ₂ = 1025 kg/m³, interface 10 cm thick:
2.42 s⁻². Sea breezes, gravity currents and oceanic fronts start the same way (Ch. 13).
""", equation=r"\frac{D\omega_z}{Dt}=\frac{2(\rho_2-\rho_1)g}{(\rho_2+\rho_1)\delta}")
D("D07", ref="5.28")
nb.code(r"""
print(ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1))   # the closed form of D07 [1/s²]: 2.422222
F = ch05.lock_exchange_fields(1000.0, 1025.0, 0.1)                   # tanh-step density and hydrostatic pressure (fluid at rest)
print([f"{v:.6f}" for v in ch05.baroclinic_term(F["rho_fn"], F["p_fn"], np.array([0.0, 0.5]))])   # the field route at the interface, mid-depth [1/s²]
print([round(ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, d), 4) for d in (0.2, 0.05, 0.01)])   # thinner → faster
""", explain=r"""
1. The closed form $2(\rho_2-\rho_1)g/((\rho_2+\rho_1)\delta)$ = 2.422222 s⁻².
2. The field route: `baroclinic_term` differentiates the tanh-step density and the hydrostatic pressure by central
   differences at the interface and returns $(0,0,(\nabla\rho\times\nabla p)_z/\rho^2)$ = (0, 0, 2.422219) — the same number
   (the 3×10⁻⁶ gap is ρ² at the interface centre versus ρ̄², D07 step 5).
3. A thinner interface spins faster: 1.2111, 4.8444, 24.2222 s⁻² for δ = 20, 5, 1 cm — in the limit it is a vortex sheet (C14).
""")
nb.figure(r"""
F = ch05.lock_exchange_fields(1000.0, 1025.0, 0.1, nx=201 if not FAST else 101, ny=101 if not FAST else 51)   # fields on a grid
ext = (F["x"][0], F["x"][-1], F["y"][0], F["y"][-1])                 # plot extent [m]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6))                  # before | baroclinic source
im = a.imshow(F["rho_grid"], origin="lower", extent=ext, cmap="Blues", aspect="auto")   # density: salt (dark) left, fresh right
a.axvline(0.0, color=COLORS["ink"], ls="--", lw=1.5)                 # the gate
a.set_title("before: salt (left) | fresh (right)", fontsize=10); a.set_xlabel("x [m]"); a.set_ylabel("y [m]")
fig.colorbar(im, ax=a, label="ρ [kg/m³]")                            # density scale
im2 = b.imshow(F["baroclinic_z"], origin="lower", extent=ext, cmap="Oranges", aspect="auto")   # (∇ρ × ∇p)_z/ρ² [1/s²]
fig.colorbar(im2, ax=b, label=r"$(\nabla\rho\times\nabla p)_z/\rho^2$ [1/s²]")   # colour scale with units
for yc in (0.25, 0.5, 0.75):                                         # counterclockwise tumbling arrows along the interface
    b.annotate("", (0.18, yc + 0.08), (0.18, yc - 0.08), arrowprops=dict(arrowstyle="-|>", color=COLORS["ink"],
               connectionstyle="arc3,rad=0.9"))                       # a curved arrow: counterclockwise
b.set_xlim(-0.5, 0.5); b.set_title("the source at t = 0⁺: a strip on the interface", fontsize=10); b.set_xlabel("x [m]")
savefig(fig, "ch05", "lock_exchange"); plt.show()                    # save the PNG to outputs/ch05, then draw
""", see="(Note `N49`, Fig. 5.5 analogue.) Left, dense salt water (dark) beside fresh water with the gate dashed; right, a bright orange vertical strip "
         "on the interface and curved arrows turning counterclockwise.",
    read="Vorticity is created only where density changes sideways — at the interface — in the counterclockwise sense: "
         "heavy fluid slides under to the right along the bottom, light fluid over to the left along the top.",
    change="…the gate between two fluids of equal density: the strip vanishes; nothing happens.")
explainer("baroclinic_torque", "How can density make fluid spin?",
          "turning the isopycnals against the isobars and watching the force line leave the centre of mass, the torque "
          "grow as the sine of the angle and two independent routes agree, makes (5.28) a mechanism instead of a formula.",
          "The lock-exchange mode shows the same torque on a whole interface.", [
              "Set the tilt to 0°: no torque however large ∇ρ.",
              "Drag the tilt to 90°: the badge reads the spin-up and its sense.",
              "Shrink the disc: the two routes converge (the ratio → 1).",
              "Lock exchange: thin the interface and watch the rate grow as 1/δ.",
          ])
whatif(r"""
…the fluid were also rotating, like the atmosphere? The baroclinic torque would still act, but the new spin would be
turned and stretched by the planet's own vorticity — both appear together in the rotating-frame vorticity equation (C09).
""")

# =====================================================================================================================
# A.3 §5.3 Helmholtz's Vortex Theorems — C05
# =====================================================================================================================
nb.section("5.3", "Helmholtz's Vortex Theorems", intro=r"""
**What is this section about?** Under Kelvin's four restrictions, vorticity is glued to the fluid: vortex lines move
with it, and a tube keeps its strength along its length and in time. Two of the four theorems are (5.4) again; the
other two follow from Kelvin in a few lines.
""")
core("C05", "Helmholtz's theorems: vortex lines move with the fluid", r"""
A smoke ring drifts across a room. Is the ring of spin always made of the same smoke — or does the swirl slide through
the air?

*In one line:* vortex lines are material lines when Kelvin's four restrictions hold — and a tube's strength is the same
along it, $-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ *(5.4)*, and constant in time, $D\Gamma/Dt=0$ *(5.8)*.
""")
problem(r"""
If vorticity stays attached to the same fluid, we can follow a vortex simply by following its fluid: the whole of §5.7
(point vortices, rings) and the potential-vorticity thinking of Ch. 13 rest on it. Helmholtz found four rules; we prove
the one that needs work.
""")
idea("""
(1) vortex lines move with the fluid       <- Kelvin on patches of a tube wall (D08)
(2) a tube has one strength along it       <- (5.4), geometry (D01)
(3) a tube cannot end in the fluid         <- (5.4)
(4) a tube keeps its strength in time      <- Kelvin on a loop round the tube (5.8)
(1), (4) need Kelvin's four restrictions; (2), (3) hold always
""")
note("N51", "Fig. 5.7", r"""
Our 3-D picture further down shows a patch S on the wall of a tube and its material image S′ some time later; we make
no separate computation of it.
""")
remind([
    ("small-patch localisation (\"for every S\")", "if an integral of a continuous function vanishes over every small patch, the function itself vanishes everywhere (Ch. 4 P112)."),
], lead="needed in the derivation below")
D("D08", ref="5.8")
note("N17", "Why \"for every patch\" matters", r"""
One patch with zero circulation proves only that the total flux through it vanishes; because *every* patch on the wall
keeps zero circulation, the flux vanishes locally (D08 step 6) — the carried wall is again a tube wall. A second proof
(the field-equation route: ω and a material element obey the same equation when ν = 0) needs the vorticity equation and
is checked at the end of C06.
""")
remind([
    ("np.cross, np.linalg.norm and @", "`np.cross(a, b)` is the cross product (row by row for arrays of vectors), `np.linalg.norm(v)` a vector's length, and `A @ b` a matrix–vector (or dot) product (Ch. 2 P67, Ch. 4)."),
])
P("P138", "loop vector area ½∮x × dx", r"""
For a closed loop, $\mathbf A_{vec}=\tfrac12\oint\mathbf x\times d\mathbf x$ is a vector whose length is the area enclosed (for
a flat loop) and whose direction is the loop's normal by the right-hand rule; for a curved loop it is ∫n dA over any
surface spanning it. The 2-D version is the surveyor's shoelace formula. We use it to ask whether a small loop still
lies on a tube wall (its normal ⟂ ω) and, in C11, for the planetary circulation.
""", code=r"""
s = 2*np.pi*np.arange(400)/400                                  # tags round a loop [rad]
P = np.array([2*np.cos(s), np.sin(s), 0*s])                     # an ellipse in the (x, y) plane, semi-axes 2 and 1 [m]
dP = np.roll(P, -1, axis=1) - P                                 # the chord to the next point (np.roll shifts the columns by one)
print(0.5*np.sum(np.cross(P.T, dP.T), axis=0), 2*np.pi)         # ≈ (0, 0, 6.28293) vs πab = 6.28319 m²: chords cut the corners
print(ch05.loop_vector_area(P))                                 # the library's spectral version: (0, 0, 6.283185)
""")
gloss(["np.roll (next point round a loop)"])
nb.worked_example("a small loop painted on a tube wall", r"""
The **ABC flow** $\mathbf u=(\sin z+\cos y,\ \sin x+\cos z,\ \sin y+\cos x)$ (in m and m/s) is a steady solution of Euler's
equation whose vorticity equals its velocity, ω = u — Kelvin holds exactly.

1. At the origin ω = (1, 1, 1) s⁻¹, |ω| = 1.732 s⁻¹.
2. Paint a circle of radius 1 cm there in the plane spanned by ω and (1, −1, 0)/√2: its normal (1, 1, −2)/√6 is
   perpendicular to ω, so the patch lies on the local tube wall.
3. Its flux ω·n A is zero to first order; the curvature of the field leaves Γ ≈ 8×10⁻¹⁰ m²/s, against
   |ω|A = 5.4×10⁻⁴ m²/s for a patch facing ω.
4. Let the flow carry it for 2 s: Kelvin keeps Γ at 8×10⁻¹⁰ m²/s, so the carried patch still faces across ω — it is
   still on a tube wall.
""")
nb.code(r"""
hsc = ch05.kelvin_scenario("helmholtz_abc")                          # the ABC flow and the painted loop (256 points, r = 1 cm)
t = np.linspace(0, 2.0, 5)                                           # 0, 0.5, …, 2 s
Pl = ch05.material_loop(hsc["u"], hsc["pts0"], t)                    # the loop carried by the flow, shape (5, 3, 256) [m]
Gh = ch05.material_circulation(hsc["u"], hsc["pts0"], t)             # its circulation [m²/s]
cosang = []                                                          # cos of the angle between the loop's normal and ω
for Pk in Pl:                                                        # at each time
    A = ch05.loop_vector_area(Pk)                                    # vector area of the carried loop [m²]
    w = hsc["u"](Pk.mean(axis=1), 0.0)                               # ω = u of the ABC flow at the loop's centre [1/s]
    cosang.append(A @ w/(np.linalg.norm(A)*np.linalg.norm(w)))       # 0 when the loop lies on a tube wall
print("Γ(t) =", [f"{g:.4e}" for g in Gh], "m²/s")                   # 8.0159e-10 m²/s every time
print(f"max |cos(normal, ω)| = {np.max(np.abs(cosang)):.1e}")        # stays ≈ 0: the normal stays at 90° to ω
print("length(t) =", np.round([ch05.loop_length(Pk) for Pk in Pl], 4), "m")   # the loop is really deformed
""", explain=r"""
1. An inviscid steady flow and a small loop lying on a tube wall.
2. The loop advected (one `solve_ivp` call for all 256 points).
3. Its circulation stays at its tiny starting value 8.0159×10⁻¹⁰ m²/s (Kelvin, $D\Gamma/Dt=0$ (5.8)).
4. The angle between the loop's vector area and ω stays 90° (|cos| ≤ 1.3×10⁻⁵) while the loop's length changes by
   about 14 %: the carried patch is still on a tube wall — Helmholtz 1 seen in numbers.
""")
nb.md("**From scratch — the vector area by the shoelace chords of the primer**, for the last carried loop:")
nb.check_agree(r"""
Pk = Pl[-1]                                                          # the loop at t = 2 s, shape (3, 256) [m]
A_hand = 0.5*np.sum(np.cross(Pk.T, (np.roll(Pk, -1, axis=1) - Pk).T), axis=0)   # ½ Σ x × (chord to the next point) [m²]
A_lib = ch05.loop_vector_area(Pk)                                    # the library's spectral version [m²]
print("chords:", [f"{a:+.5e}" for a in A_hand], " library:", [f"{a:+.5e}" for a in A_lib])   # the same vector to ~1e-4 relative
flux = A_hand @ hsc["u"](Pk.mean(axis=1), 0.0)                       # vorticity flux through the carried patch [m²/s]
print(f"flux through the carried patch {flux:.1e} m²/s vs |ω||A| = {np.linalg.norm(A_hand)*np.sqrt(3):.1e} m²/s")
assert np.allclose(A_hand, A_lib, rtol=1e-3)                         # chords vs spectral: same area vector
assert abs(flux) < 1e-8                                              # no vorticity pierces the carried patch
""")
nb.md("Our own chord sum says the same: no vorticity pierces the carried patch (10⁻⁹ against 3×10⁻⁴ m²/s for a patch facing ω).")
nb.plotly(r"""
zs = np.linspace(0, 1.5, 40 if not FAST else 24)                      # heights along the narrowing tube [m]
ph = np.linspace(0, 2*np.pi, 60 if not FAST else 36)                  # angles round it [rad]
ZZ, PH = np.meshgrid(zs, ph)                                          # surface grid
RT = 0.1*np.exp(-ZZ/2)                                                # tube wall radius R_t(z) = R0 e^{−z/2L} [m]
fig = go.Figure(go.Surface(x=RT*np.cos(PH), y=RT*np.sin(PH), z=ZZ, opacity=0.35, showscale=False,   # the translucent tube wall
                           colorscale=[[0, COLORS["teal"]], [1, COLORS["teal"]]], name="tube wall"))  # one flat teal colour
for p0 in np.linspace(0, 2*np.pi, 12, endpoint=False):               # vortex lines on the wall (untwisted tube: φ = const)
    fig.add_trace(go.Scatter3d(x=0.1*np.exp(-zs/2)*np.cos(p0), y=0.1*np.exp(-zs/2)*np.sin(p0), z=zs, mode="lines",  # a vortex line on the wall
                               line=dict(color=COLORS["teal"], width=3), showlegend=False))         # thin teal line


def patch(z0, z1, p0, p1, col, name):                                # a patch of the wall between two heights and angles
    zp, pp = np.meshgrid(np.linspace(z0, z1, 12), np.linspace(p0, p1, 12))                          # the patch's (z, φ) grid
    rp = 0.1*np.exp(-zp/2)*1.002                                      # just outside the wall so it shows
    return go.Surface(x=rp*np.cos(pp), y=rp*np.sin(pp), z=zp, showscale=False, opacity=0.95,        # the patch as a small surface
                      colorscale=[[0, col], [1, col]], name=name)                                   # one flat colour


fig.add_trace(patch(0.10, 0.30, 0.0, 0.6, COLORS["orange"], "S (t = 0)"))              # the painted patch
fig.add_trace(patch(0.60, 1.20, 0.0, 0.6, COLORS["amber"], "S′ (later, schematic)"))  # the same particles, stretched along the tube
fig.update_layout(height=480, margin=dict(l=0, r=0, t=40, b=0), title="Schematic: a patch on a tube wall stays on the wall",  # size, margins, title
                  scene=dict(aspectmode="data", xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]"))  # equal axes, labels with units
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="A translucent teal tube narrowing upward with its vortex lines drawn on the wall; an orange patch low down and "
        "an amber patch higher up, longer along the tube and narrower round it.",
    read="Helmholtz 1: the patch's edge has zero circulation for ever, so it never tilts off the wall — the particles of "
         "S become S′, still lying on the wall between the same vortex lines. (A schematic: S′ is drawn, not computed; "
         "the ABC cell above is the computation.)",
    change="…viscosity: the patch's circulation would leak and the lines could reconnect — the patch would drift off the "
           "wall on the diffusion time scale.")
nb.md(r"""
🎮 Open **Stretch and tangle a loop of dye — what stays the same?** (embedded in C03 above) and switch on its
*Helmholtz* mode: a small loop on a tube wall keeps zero flux while it is carried.
""")
whatif(r"""
…the vortex were a straight line vortex in a flow that carries it along? Then Helmholtz 1 says the vortex goes wherever
the local flow (made by everything else) takes it — the rule of §5.7.
""")

# =====================================================================================================================
# A.4 §5.4 Vorticity Equation in a Nonrotating Frame — R10, R11, C06
# =====================================================================================================================
nb.section("5.4", "Vorticity Equation in a Nonrotating Frame", intro=r"""
**What is this section about?** Taking the curl of Navier–Stokes removes pressure and gravity and leaves an equation
for the vorticity alone: vorticity is carried by the flow, stretched and tilted by velocity gradients along it, and
diffused by viscosity. This section derives it for a fluid of constant density in an inertial frame; §5.6 adds rotation
and density variations.
""")
nb.recap("R10", "Vorticity has no divergence", r"""
$\nabla\cdot\boldsymbol\omega=\nabla\cdot(\nabla\times\mathbf u)=0$ for every smooth flow — the divergence of a curl
vanishes (Ch. 2). §5.6 proves it again in index form as $\omega_{i,i}=\varepsilon_{inq}u_{q,ni}=0$ *(5.18)*.
""", where="Ch. 2 §2.13")
nb.recap("R11", "The Lamb identity", r"""
$(\mathbf u\cdot\nabla)\mathbf u=\nabla(\tfrac12\mathbf u\cdot\mathbf u)+\boldsymbol\omega\times\mathbf u$, the vector form of
$u_i\frac{\partial u_j}{\partial x_i}=-(\mathbf u\times\boldsymbol\omega)_j+\frac{\partial}{\partial x_j}\big(\tfrac12u_i^2\big)$
*(4.68)*: the advective acceleration is the gradient of the kinetic energy plus the Lamb vector ω × u. ⚠️ The book's
§5.4 line writes $\nabla(\mathbf u\cdot\mathbf u)$ without the ½ — harmless there, because the next move takes a curl,
which kills any gradient.
""", where="Ch. 4 §4.9")
nb.code(r"""
lt = ns.lamb_identity_terms(vortices.vortex_velocity_field("gaussian", Gamma=1.0, sigma=0.1), np.array([0.07, 0.02]))   # Gaussian vortex, one point
print("advective (u·∇)u:", lt["advective"], " residual:", [f"{v:.1e}" for v in lt["residual"]])   # residual ~1e-6 against terms of ~16 m/s²: stencil error
""")
remind([
    ("Schwarz's theorem", r"mixed partial derivatives of a smooth function commute: $\partial^2f/\partial x\partial y=\partial^2f/\partial y\partial x$ (Ch. 4 P121)."),
    ("directional derivative (ω·∇)", r"$(\boldsymbol\omega\cdot\nabla)\mathbf u$ is |ω| times the rate of change of u along ω's direction (Ch. 2 P75)."),
])

# ---- C06 ----------------------------------------------------------------------------------------------------------
core("C06", "The vorticity equation", r"""
Stir a cup of tea and the swirl slowly spreads and fades; a tornado tightens as its column is stretched. What can change
a fluid particle's vorticity — and why can pressure never do it?
""", eqs=("5.13",))
problem(r"""
Navier–Stokes tells how a particle's velocity changes. For swirl, we want the rule for its vorticity. Pressure is usually
the hardest unknown in a flow; an equation without it is a gift — the basis of vorticity-based computer models
(Ch. 10), of vortex decay (Ch. 8), of the turbulence cascade (Ch. 12) and of the quasi-geostrophic equations of Ch. 13.
""")
idea("""
momentum:   Du/Dt = -grad(p)/rho + g + nu lap(u)      take the curl  ---->   vorticity:  Dw/Dt = (w.grad)u + nu lap(w)   (5.13)
                    ^^^^^^^^^^^^   ^                                                  carried    stretched/tilted   diffused
            gradients: their curl is 0 (pressure and gravity push through the centre of mass: no twist)
""")
note("N18", "§5.4's hypotheses", r"""
constant density (so the fluid is barotropic), constant viscosity, conservative body force, inertial frame. §5.6 (C09)
relaxes the first and the last.
""")
note("N19", "The first move", r"""
of the derivation is (5.12), the curl of incompressible Navier–Stokes
$\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(4.39b)*. The curls of −∇p/ρ and of g = −∇Φ vanish
because both are gradients.
""", equation=EQ["5.12"], ref="5.12")
note("N20", "A vector identity we need", r"""
(the book's (B.3.10)), in its general form (below). It follows from the ε–δ identity
$\varepsilon_{ijk}\varepsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}$ *(2.19)* and the product rule (Ch. 2 D09's
moves); the sympy cell of D09 checks it for arbitrary fields. For incompressible flow the last two terms drop.
""", equation=r"\nabla\times(\boldsymbol\omega\times\mathbf u)=(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u+\boldsymbol\omega(\nabla\cdot\mathbf u)-\mathbf u(\nabla\cdot\boldsymbol\omega)")
D("D09", ref="5.13", check_src=r"""
x, y, z, t, nu = sp.symbols('x y z t nu')                            # coordinates [m], time [s], kinematic viscosity [m²/s]
P, Q, R = [sp.Function(n)(x, y, z, t) for n in 'PQR']                # a generic velocity u = (P, Q, R): any smooth field
p = sp.Function('p')(x, y, z, t)                                     # a generic pressure
Phi = sp.Function('Phi')(x, y, z)                                    # a generic body-force potential, g = −∇Φ
terms = ch05.vorticity_equation_sym([P, Q, R], (x, y, z), t, nu, p_expr=p, Phi_expr=Phi)   # step 1: the curl of every term of (4.39b)
print("steps 2–3, curl of the pressure term:", [sp.simplify(c) for c in terms['curl_pressure']])   # [0, 0, 0]: a gradient
print("steps 2–3, curl of the gravity term: ", [sp.simplify(c) for c in terms['curl_gravity']])    # [0, 0, 0]: a gradient
print("step 9, identity (B.3.10), generic u:", [sp.simplify(c) for c in terms['identity_B310']])   # [0, 0, 0] for ANY fields
u0 = [y*z**2 + sp.sin(x)*sp.cos(y), x*z - sp.cos(x)*sp.sin(y), -x*y]   # divergence-free: a polynomial plus the Taylor–Green (cellular) pattern
print("its divergence:", sp.simplify(sum(sp.diff(ui, xi) for ui, xi in zip(u0, (x, y, z)))))      # 0
r513 = ch05.vorticity_equation_sym(u0, (x, y, z), t, nu)['residual_513']   # curl of the NS residual minus the (5.13) residual
print("curl(Navier–Stokes) − (5.13):", [sp.simplify(c) for c in r513])     # [0, 0, 0]: (5.13) is the curl of NS
""")
nb.worked_example("the three terms at one point of Burgers' vortex", r"""
Burgers' vortex (derived in C10): $u_R=-\alpha R/2$, $u_z=\alpha z$, $\omega_z=(\alpha\Gamma/4\pi\nu)e^{-\alpha R^2/4\nu}$ with
α = 1 s⁻¹, Γ = 10⁻³ m²/s, ν = 10⁻⁶ m²/s. At R = 1 mm:

1. $\omega_z=79.58\,e^{-0.25}$ = 61.97 s⁻¹.
2. Its radial slope: $\omega_z'=-(\alpha R/2\nu)\omega_z$ = −500 × 61.97 = −3.10×10⁴ s⁻¹/m.
3. Advection: $u_R\omega_z'$ = (−0.5×10⁻³)(−3.10×10⁴) = +15.49 s⁻².
4. Stretching: $\omega_z\,\partial u_z/\partial z$ = 61.97 × 1 = +61.97 s⁻².
5. Steady, so (5.13) says diffusion = advection − stretching = 15.49 − 61.97 = −46.48 s⁻²: viscosity carries vorticity
   outward, stretching makes more, the inflow brings it back.
""")
nb.md(r"""
*Gloss — the Taylor–Green vortex.* The **Taylor–Green vortex** is the viscous version of C03's cellular flow:
$\mathbf u=U_0(\sin x\cos y,\ -\cos x\sin y)\,e^{-2\nu t}$ (the preset uses the same flow shifted, $(\cos x\sin y,\ -\sin x\cos y)$),
an exact 2-D Navier–Stokes solution in which $\nabla^2\mathbf u=-2\mathbf u$: its vorticity is only diffused, never
stretched, and every part of it decays at the same rate. The **Lamb–Oseen** vortex is the spreading Gaussian vortex
$u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)$ *(3.29)* with $\sigma^2=4\nu(t+t_0)$.
""")
nb.code(r"""
for name in ("burgers", "lamb_oseen", "taylor_green"):               # three classic flows, each at its preset probe point
    b = ch05.vorticity_budget_preset(name, component=2)              # every term of (5.13) (z-component) [1/s²]
    print(f"{name:13s}" + "  ".join(f"{k} {v:+.4f}" for k, v in b.items() if k not in ("planetary", "baroclinic")))
u_b = vortices.burgers_vortex_field(1e-3, 1.0, 1e-6)                 # Burgers' vortex as a velocity function (C10)
T = ch05.vorticity_terms(u_b, np.array([1e-3, 0.0, 0.0]), nu=1e-6, h=1e-5)   # the general call at R = 1 mm, stencil h = 10 μm
print({k: round(float(getattr(T, k)[2]), 4) for k in T._fields})     # z-components: advective 15.49, stretching 61.97, diffusion −46.48
""", explain=r"""
1. `vorticity_budget_preset` (one component at the preset's probe point) for three classic flows: Burgers at R = 1 mm —
   advective +15.4937, stretching +61.9750, diffusion −46.4813, residual 3×10⁻⁵ (printed as +0.0000); Lamb–Oseen at r = 5 mm (t₀ = 10 s) —
   no stretching, local = diffusion = −1.5973; Taylor–Green (2-D) — no stretching, local = diffusion.
2. The general call: vorticity by stencils of u, then each term of (5.13) by stencils of ω (h = 10⁻⁵ m, about R/100).
3. The residual (local + advective − stretching − diffusion) is the check that (5.13) holds: ~10⁻³ here against terms of
   ~60 s⁻² (nested stencils).
""")
nb.md("**From scratch — $(\\boldsymbol\\omega\\cdot\\nabla)\\mathbf u$ is just G times ω, and $\\nabla^2\\boldsymbol\\omega$ is three second differences:**")
nb.check_agree(r"""
x0, h = np.array([1e-3, 0.0, 0.0]), 1e-5                             # the probe point (R = 1 mm) and the stencil step [m]
w = lambda x: kinematics.vorticity(u_b, x, h=h)                      # ω(x) by central differences of u [1/s]
G = kinematics.velocity_gradient_at(u_b, x0, h=h)                    # velocity gradient tensor G_ij = ∂u_i/∂x_j [1/s]
stretch = G @ w(x0)                                                  # (ω·∇)u = G ω [1/s²]
lap = sum((w(x0 + h*e) - 2*w(x0) + w(x0 - h*e))/h**2 for e in np.eye(3))   # ∇²ω: second differences along x, y, z
print(f"stretching {stretch[2]:.4f}, diffusion ν∇²ω {1e-6*lap[2]:.4f} 1/s²")
assert np.allclose(stretch[2], 61.97, rtol=1e-3)                     # the tiny example's number
assert np.allclose(1e-6*lap[2], -46.48, rtol=1e-3)                   # the tiny example's number
assert np.allclose([stretch[2], 1e-6*lap[2]], [T.stretching_tilting[2], T.diffusion[2]], rtol=1e-6)   # = the library
""")
nb.figure(r"""
names = ["lamb_oseen", "taylor_green", "burgers"]                     # three flows, three balances
labels = ["Lamb–Oseen, r = 5 mm", "Taylor–Green (2-D)", "Burgers, R = 1 mm"]
keys = ["local", "advective", "stretching_tilting", "diffusion", "residual"]   # the bars of (5.13)
cols = [COLORS["blue"], COLORS["muted"], COLORS["accent"], COLORS["rose"], COLORS["ink"]]   # house colours
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))                     # one panel per flow
for a, nm, lab in zip(axs, names, labels):
    b = ch05.vorticity_budget_preset(nm, component=2)                 # the terms at the preset probe [1/s²]
    a.bar(range(5), [b[k] for k in keys], color=cols)                # one bar per term
    a.set_yscale("symlog", linthresh=1e-3)                           # symmetric log: signs kept, sizes compressed
    a.axhline(0, color=COLORS["muted"], lw=0.8); a.set_xticks(range(5))
    a.set_xticklabels(["∂ω/∂t", "(u·∇)ω", "(ω·∇)u", "ν∇²ω", "residual"], fontsize=8)
    a.set_title(lab, fontsize=10)
axs[0].set_ylabel("term of (5.13), z-component [1/s²]")              # axis label with units
fig.suptitle("Three flows, three balances", fontsize=11)
savefig(fig, "ch05", "vorticity_terms"); plt.show()                  # save the PNG to outputs/ch05, then draw
""", see="Lamb–Oseen: only a local bar (blue) and an equal diffusion bar (rose); Taylor–Green: the same pair, no "
         "stretching bar; Burgers: no local bar — the grey advective bar (+15) equals the purple stretching bar (+62) "
         "plus the negative rose diffusion bar (−46).",
    read=r"The y-axis is symmetric-log (a symlog axis: linear near zero, logarithmic beyond ±10⁻³), so bars of very "
         r"different sizes and both signs fit. In 2-D flows vorticity is only carried and diffused; stretching needs a velocity "
         r"gradient along ω (C10). The residual bars sit at the stencil-error floor (2.6×10⁻⁵ s⁻² for Burgers): $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ holds.",
    change="…ν halved in Burgers: the core narrows by √2 and every bar at a fixed R changes, but advection still equals "
           "stretching plus diffusion.")
remind([
    ("complementary error function erfc", r"$\mathrm{erfc}(\eta)=\frac2{\sqrt\pi}\int_\eta^\infty e^{-s^2}ds$ (Ch. 4 P123); here we use $\mathrm{erf}=1-\mathrm{erfc}$, which runs from −1 to 1."),
    ("finite differences and FTCS", "replace derivatives by differences on a grid and step forward in time (forward-time, centred-space) (Ch. 1 P21)."),
])
note("N22", "A sheet of vorticity diffusing", r"""
Start with all the vorticity on a plane, $\boldsymbol\omega=\gamma\delta(y)\mathbf e_z$ (a velocity jump of γ, C14). (5.13)
reduces to 1-D diffusion, $\partial\omega_z/\partial t=\nu\partial^2\omega_z/\partial y^2$ (the heat equation of Ch. 1 §1.5), whose
solution (the book leaves it to Exercise 5.6) is below, with velocity
$u=-\frac\gamma2\,\mathrm{erf}\big(\frac{y}{2\sqrt{\nu t}}\big)$ (gloss: erf η = 1 − erfc η). The total ∫ω dy = γ never changes; the
width $2\sqrt{\nu t}$ grows. **Number:** γ = 1 m/s in water, after 1 s the layer is 2 mm thick and the peak vorticity
282 s⁻¹. → Stokes' first problem (Ch. 8) is half of this; shear layers (Ch. 11).
""", equation=r"\omega_z(y,t)=\frac{\gamma}{2\sqrt{\pi\nu t}}\exp\Big\{-\frac{y^2}{4\nu t}\Big\}")
nb.code(r"""
y = np.linspace(-0.05, 0.05, 4001)                                   # across the layer, ±5 cm [m]
for t in (0.1, 1.0, 10.0):                                           # three times [s]
    u, w = ch05.diffusing_vortex_sheet(y, t, 1.0, 1e-6)              # u(y, t) [m/s], ω_z(y, t) [1/s]; γ = 1 m/s, water
    print(f"t = {t:4.1f} s: ∫ω dy = {np.trapezoid(w, y):.6f} m/s, peak ω {w.max():7.1f} 1/s, u(−5 cm) = {u[0]:+.4f}, u(+5 cm) = {u[-1]:+.4f} m/s")
yf = np.linspace(-0.01, 0.01, 2001 if not FAST else 1001); dy = yf[1] - yf[0]   # the FTCS grid, ±1 cm [m]
w0 = ch05.diffusing_vortex_sheet(yf, 0.1, 1.0, 1e-6)[1]              # start from the exact profile at t = 0.1 s
dt = 0.4*dy**2/1e-6; nsteps = int(round(0.9/dt)); dt = 0.9/nsteps    # a stable time step (ν dt/dy² = 0.4), 0.9 s of diffusion
w_num = diffusion.ftcs_diffusion_1d(w0, 1e-6, dy, dt, nsteps, bc=("dirichlet", "dirichlet"), values=(0.0, 0.0), save_every=nsteps)[-1]
w_ex = ch05.diffusing_vortex_sheet(yf, 1.0, 1.0, 1e-6)[1]            # the exact profile at t = 1 s
print(f"FTCS vs exact at t = 1 s: max relative difference {np.max(np.abs(w_num - w_ex))/w_ex.max():.1e}")
""", explain=r"""
1. The closed form for three times: the total ∫ω dy = 1.000000 m/s = γ every time; the peak falls as $1/\sqrt t$ —
   892.1, 282.1 and 89.2 s⁻¹; far away the velocity is ∓γ/2 = ∓0.5 m/s.
2. An independent finite-difference (FTCS) solution of the same diffusion equation, started from the exact profile at
   0.1 s, agrees with the exact profile at 1 s to about 10⁻⁵ of the peak (8×10⁻⁶ on the full grid; 3×10⁻⁵ on the
   coarser grid of a FAST run).
""")
nb.plotly(r"""
ys = np.linspace(-0.02, 0.02, 401)                                   # ±2 cm [m]


def sheet(t):                                                        # time [s] → the two profiles
    u, w = ch05.diffusing_vortex_sheet(ys, t, 1.0, 1e-6)             # γ = 1 m/s, ν = 1e-6 m²/s
    return {"ω_z [10³ 1/s]": (ys*100, w/1e3), "u [m/s]": (ys*100, u)}   # y in cm on the x-axis


fig = slider_figure(sheet, "t", np.geomspace(0.05, 10, 20 if not FAST else 10), unit="s", xlabel="y [cm]",  # precompute every slider position
                    ylabel="ω_z [10³ 1/s] (teal) · u [m/s] (orange)",
                    title="A vortex sheet diffusing: ∫ω dy = γ = 1 m/s at every t (printed above)")  # the message of the figure
recolor(fig, {"ω_z [10³ 1/s]": COLORS["teal"], "u [m/s]": COLORS["orange"]})                        # house colours by trace name
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="A sharp teal spike that widens and sinks as t grows; an orange velocity step from +0.5 to −0.5 m/s that softens.",
    read="Vorticity spreads like heat; its total — the velocity jump, γ — is conserved (the area under the teal curve does "
         "not change). The orange curve's slope at y = 0 is −ω_z there.",
    change="…ν fifteen times larger (air, 1.5×10⁻⁵ m²/s, instead of water): the same shapes at one fifteenth of the time.")
note("N23", "Hill's spherical vortex", r"""
is a second exact solution of (5.13), steady and inviscid (the book leaves it to Exercise 5.11): inside a sphere of
radius a, $\psi=\frac{Aa^4}{10}\frac{R^2}{a^2}\Big(1-\frac{R^2}{a^2}-\frac{z^2}{a^2}\Big)$ (an axisymmetric Stokes stream function,
Ch. 4 §4.3) and $\boldsymbol\omega=AR\,\mathbf e_\varphi$ — ring-shaped vortex lines stretched just enough to keep ω/R constant on
each streamline. The flow outside the sphere is Ch. 6's.
""")
nb.code(r"""
print(ch05.vorticity_budget_preset("hill", x=0.3, z=0.2, component=1))   # ω_φ-component at (0.3, 0, 0.2): advective = stretching, residual ≈ 1e-13
""")
nb.md(r"""
*↪ N21 → C10:* Burgers' vortex (an exercise attached to this section) is taught with stretching in C10, where its
stretching–diffusion balance is the point (D18).
""")
nb.md(r"""
**Helmholtz's first theorem from the field equation (C05, second proof).** With ν = 0, (5.13) reads
$D\boldsymbol\omega/Dt=(\boldsymbol\omega\cdot\nabla)\mathbf u=\mathsf G\boldsymbol\omega$ — the same linear equation as a material
line element, $D(\delta\mathbf x)/Dt=\mathsf G\,\delta\mathbf x$ (Ch. 3 P99). Start δx along ω and they stay parallel, with
|ω|/|δx| constant (Cauchy's frozen-in solution): vortex lines are material. `frozen_in_check` integrates the particle,
δx and ω together along one path.
""")
nb.code(r"""
u_abc = ch05.abc_flow()                                              # a steady Euler flow with ω = u (inviscid)
fr = ch05.frozen_in_check(u_abc, np.zeros(3), 1e-3*u_abc(np.zeros(3), 0.0), (0, 5), np.linspace(0, 5, 11))   # δx starts along ω
print(f"inviscid ABC: max angle(δx, ω) {np.max(fr['angle']):.1e} rad, spread of |ω|/|δx| {np.ptp(fr['ratio']):.1e}")
fv = ch05.frozen_in_check(u_b, np.array([1e-3, 0, 1e-3]), np.array([0, 0, 1e-6]), (0, 2), np.linspace(0, 2, 5), nu=1e-6)   # Burgers, ν ≠ 0
print("viscous Burgers: |ω|/|δx| =", np.round(fv["ratio"], 4))       # falls from 1: diffusion unhooks ω from the fluid
""", explain=r"""
1. In the inviscid ABC flow δx stays parallel to ω (angle ~10⁻¹³ rad) and |ω|/|δx| stays 1 to round-off — vortex lines
   are material.
2. In the viscous Burgers vortex δx and ω both stay along the axis (angle 0 by symmetry), but |ω|/|δx| falls from 1 to
   0.17 in 2 s: the element keeps being stretched while the vorticity at the particle is held by the
   stretching–diffusion balance — viscosity lets vorticity move relative to the fluid.
""")
nb.figure(r"""
fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.4))                  # angle | ratio
a.semilogy(fr["t"], np.maximum(fr["angle"], 1e-16), color=COLORS["teal"], marker="o", label="inviscid ABC")   # ~1e-13
a.set_xlabel("t [s]"); a.set_ylabel("angle between δx and ω [rad]"); a.set_ylim(1e-17, 1e-2); a.legend(fontsize=8)
a.set_title("(a) δx stays along ω (log axis)", fontsize=10)
b.plot(fr["t"], fr["ratio"], color=COLORS["teal"], marker="o", label="inviscid ABC")          # flat at 1
b.plot(fv["t"], fv["ratio"], color=COLORS["rose"], marker="s", label="viscous Burgers")       # drifts
b.set_xlabel("t [s]"); b.set_ylabel("|ω|/|δx|, normalised [–]"); b.legend(fontsize=8); b.set_title("(b) the ratio", fontsize=10)
plt.show()
""", see="Left, the inviscid angle sits at the round-off floor (~10⁻¹³ rad); right, the teal ratio stays exactly 1 while "
         "the rose (viscous) ratio falls to about 0.17.",
    read="Inviscid: vorticity and a material element obey the same equation, so they stay locked together (Helmholtz 1 "
         "and 4). Viscous: diffusion moves vorticity relative to the fluid, so the lock breaks.",
    change="…ν → 10ν in Burgers: the core widens (√(4ν/α)) and the ratio falls away from 1 even faster in the first instants.")
nb.md(r"""
> ⚠️ **Common confusion:** "(5.13) is linear in ω." It only looks linear: u is itself determined by ω (Biot–Savart,
> C07), so the advective term $(\mathbf u\cdot\nabla)\boldsymbol\omega$ and the stretching term $(\boldsymbol\omega\cdot\nabla)\mathbf u$
> are both quadratic in the vorticity.
""")
nb.md(r"""
🎮 **How can a flow spin fluid faster without a torque?** (embedded after C10 in §5.6) steps through D09 in its
Derivation tab.
""")
whatif(r"""
…the density varied or the frame rotated? The pressure curl would no longer vanish (a baroclinic source, C04) and the
Coriolis term would add the planet's vorticity to the stretching term —
$\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$
*(5.30)*, C09.
""")

# =====================================================================================================================
# A.5 §5.5 Velocity Induced by a Vortex Filament: Law of Biot and Savart — C07, C08
# =====================================================================================================================
nb.section("5.5", "Velocity Induced by a Vortex Filament: Law of Biot and Savart", intro=r"""
**What is this section about?** If you know the vorticity everywhere, you know the velocity everywhere (up to a
potential flow). This section derives the formula — Biot–Savart, the same law that gives a magnetic field from an
electric current — and its working form for thin vortex filaments. The book's printed (5.14) carries a sign slip; we
derive the correct sign.
""")

# ---- C07 ----------------------------------------------------------------------------------------------------------
core("C07", "Biot–Savart: velocity from vorticity", r"""
A vortex spins in one corner of a pond. Water on the far side moves too, although nothing touches it there. Given where
the vorticity is, is the velocity everywhere decided?
""", eqs=("5.16",))
problem(r"""
Aerodynamicists replace a wing by the vortices it sheds and compute the downwash they induce; oceanographers invert a
map of vorticity for the currents; vortex methods in CFD move blobs of vorticity with the velocity they induce on each
other. All of them need one formula: velocity from vorticity.
""")
idea("""
w = curl u,  div u = 0   --curl again-->   lap(u) = -curl(w)          (a Poisson equation: like gravity or electrostatics)
point-source solution of lap: -1/(4 pi r)   --add up-->   u = (1/4pi) Int (curl' w)/r d3x'      (5.14), sign corrected
move the curl off w (product rule + Gauss)  -->  u = (1/4pi) Int w x (x-x')/r^3 d3x'           (5.16) = current -> magnetic field
""")
P("P139", "Poisson equation and Green's function", r"""
A Poisson equation ∇²φ = q asks for the field produced by a source density q (gravity from mass, voltage from charge).
Its building block is the response to one point source, the Green's function: in 3-D, $G=-1/(4\pi\lvert\mathbf x-\mathbf x'\rvert)$
solves $\nabla^2G=\delta(\mathbf x-\mathbf x')$, the Dirac delta being a unit source concentrated at one point (Ch. 3 met it
as a gloss). Because the equation is linear, the response to any q is the sum of point responses,
$\phi(\mathbf x)=\int G(\mathbf x,\mathbf x')q(\mathbf x')d^3x'$ (superposition).
""", code=r"""
x, y, z = sp.symbols('x y z')                                   # a point in space; the source sits at the origin
r = sp.sqrt(x**2 + y**2 + z**2)                                 # distance from the source [m]
G = -1/(4*sp.pi*r)                                              # the Green's function of ∇²
print(sp.simplify(sum(sp.diff(G, v, 2) for v in (x, y, z))))    # ∇²G = 0 away from the source
""")
gloss(["Dirac delta as a point source", "superposition for a linear equation"])
P("P140", "gradient of 1/distance with respect to the source point", r"""
With x fixed and the source point x′ moving, $\nabla'(1/\lvert\mathbf x-\mathbf x'\rvert)=+(\mathbf x-\mathbf x')/\lvert\mathbf x-\mathbf x'\rvert^3$
— it points from the source toward the field point. Differentiating with respect to x instead flips the sign; the
book's two sign slips in §5.5 are exactly this trap.
""", code=r"""
X = sp.symbols('x y z'); Xp = sp.symbols('xp yp zp')            # field point x, source point x′
d = [a - b for a, b in zip(X, Xp)]; r = sp.sqrt(sum(di**2 for di in d))   # x − x′ and its length
grad_p = [sp.diff(1/r, v) for v in Xp]                          # ∇′(1/r): derivatives with respect to x′
print([sp.simplify(g - di/r**3) for g, di in zip(grad_p, d)])   # [0, 0, 0]: ∇′(1/r) = +(x − x′)/r³
""")
P("P141", "curl of a product (product rule)", r"""
For a scalar f and a vector A, $\nabla\times(f\mathbf A)=f\nabla\times\mathbf A+\nabla f\times\mathbf A$ — the product rule, with the
order of the cross product kept. D11 uses it to move the curl off ω. (Gloss: the cross product is antisymmetric,
**a** × **b** = −**b** × **a**, so reordering it flips the sign.)
""", code=r"""
x, y, z = sp.symbols('x y z'); f = sp.Function('f')(x, y, z)          # a generic scalar field
A = sp.Matrix([sp.Function(n)(x, y, z) for n in 'PQR'])              # a generic vector field
curl = lambda V: sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], z), sp.diff(V[0], z) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)])   # ∇×
grad = sp.Matrix([sp.diff(f, v) for v in (x, y, z)])                  # ∇f
print(sp.simplify(curl(f*A) - f*curl(A) - grad.cross(A)))             # zero vector: the product rule holds
""")
gloss(["cross-product antisymmetry a × b = −b × a"])
remind([
    ("curl of a curl identity", r"$\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u$ (Ch. 4 P122; the book's (B.3.13))."),
])
note("N24", "Velocity obeys a Poisson equation", r"""
For incompressible flow the curl of the vorticity is minus the Laplacian of the velocity (the curl-of-curl identity,
Ch. 4 P122). On a periodic box the Poisson equation is solved in one line with Fourier modes — the next primer and cell.
→ Vorticity–stream-function CFD (Ch. 10), PV inversion (Ch. 13).
""", equation=r"\nabla\times\boldsymbol\omega=\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u=-\nabla^2\mathbf u")
P("P142", "Fourier modes and the FFT Poisson solver", r"""
On a periodic box every smooth field is a sum of waves $e^{i\mathbf k\cdot\mathbf x}$; a derivative ∂/∂x multiplies a wave by
$ik_x$ and ∇² by $-\lvert\mathbf k\rvert^2$. So ∇²ψ = −ω becomes $\hat\psi=\hat\omega/\lvert\mathbf k\rvert^2$ for each wave (k = 0 is
the mean, set to 0). `numpy.fft.fft2` finds the waves, `ifft2` adds them back.
""", code=r"""
n, L = 32, 2*np.pi; x = np.arange(n)*L/n; X, Y = np.meshgrid(x, x)   # a periodic 32 × 32 grid on [0, 2π)² [m]
w = 2*np.sin(X)*np.sin(Y)                                            # vorticity of ψ = sin x sin y [1/s]
k = np.fft.fftfreq(n, L/n)*2*np.pi; KX, KY = np.meshgrid(k, k); K2 = KX**2 + KY**2; K2[0, 0] = 1   # wavenumbers; avoid 0/0
psi = np.real(np.fft.ifft2(np.fft.fft2(w)/K2))                       # ψ̂ = ω̂/k², back to the grid
print(np.max(abs(psi - np.sin(X)*np.sin(Y))))                        # ~1e-15: the exact stream function
""")
nb.code(r"""
U = ch05.velocity_from_vorticity_fft(w, 2*np.pi)                     # the library: (u, v) on the same grid [m/s]
print(np.max(abs(U[0] - np.sin(X)*np.cos(Y))), np.max(abs(U[1] + np.cos(X)*np.sin(Y))))   # cellular flow recovered to round-off
""")
D("D10", ref="5.14", check_src=r"""
x, y, z = sp.symbols('x y z', real=True)                             # a field point; the source sits at the origin
X = (x, y, z)                                                        # the coordinates as a tuple
A, B, C = [sp.Function(n)(x, y, z) for n in 'ABC']                   # a generic velocity field for step 2
U = sp.Matrix([A, B, C])
curl = lambda V: sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], z), sp.diff(V[0], z) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)])   # ∇×
grad = lambda f: sp.Matrix([sp.diff(f, v) for v in X])               # ∇
div = lambda V: sum(sp.diff(V[i], X[i]) for i in range(3))           # ∇·
lapv = lambda V: sp.Matrix([sum(sp.diff(V[i], v, 2) for v in X) for i in range(3)])   # vector Laplacian (Cartesian)
print("step 2, curl curl u − [∇(∇·u) − ∇²u]:", sp.simplify(curl(curl(U)) - (grad(div(U)) - lapv(U))).T)   # zero vector
r = sp.sqrt(x**2 + y**2 + z**2)                                      # distance from the source [m]
print("step 6, ∇²(1/r) for r > 0:", sp.simplify(sum(sp.diff(1/r, v, 2) for v in X)))     # 0
rr, eps, th, ph = sp.symbols('r epsilon theta phi', positive=True)   # radius, sphere radius, sphere angles
dfdr = sp.diff(1/rr, rr).subs(rr, eps)                               # radial derivative of 1/r on the sphere: −1/ε²
flux = sp.integrate(sp.integrate(dfdr*eps**2*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))   # step 7: ∮∇(1/r)·n dA
print("step 7, flux of ∇(1/r) through a sphere of radius ε:", flux)  # −4π, whatever ε is
print("step 8, flux of ∇G with G = −1/(4πr):", sp.simplify(flux*(-1/(4*sp.pi))))           # 1: a unit source
F = ch05.gaussian_tube_fields(Gamma=1.0, sigma=0.1, L=4.0)           # step 11 with numbers: a smooth tube of length 4 m
nodes, weights = ch05.cylinder_quadrature(F["bounds"]["radius"], F["bounds"]["z0"], F["bounds"]["z1"], nr=24, nphi=32, nz=24)
for s in (+1, -1):                                                   # the corrected sign, then the book's printed sign
    print(f"step 11, sign {s:+d}: u_θ(0.5 m) = {ch05.velocity_from_curl_omega(F['curl_omega'], np.array([0.5, 0, 0]), nodes, weights, sign=s)[1]:+.4f} m/s")
print(f"exact for this finite tube: {F['u_theta_reference'](0.5):.6f} m/s (counterclockwise, like the vortex)")
""")
note("N25", "(5.14) with its sign corrected", r"""
D10's result is below. ⚠️ **The book prints**
$\mathbf u(\mathbf x,t)=-\frac{1}{4\pi}\int_{V'}\frac{1}{\lvert\mathbf x-\mathbf x'\rvert}(\nabla'\times\boldsymbol\omega(\mathbf x',t))d^3x'$
**(5.14) with −1/(4π); the correct factor is +1/(4π).** It also cites "Exercise 5.8" for the proof; the Green's-function
exercise is 5.9. The code below computes both signs on a smooth vortex tube — only one of them spins the right way.
""", equation=r"\mathbf u(\mathbf x,t)=\frac1{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega(\mathbf x',t)}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'", ref="5.14")
P("P143", "Gauss–Legendre quadrature in 3-D and a smoothed kernel", r"""
Gauss–Legendre picks n nodes and weights per axis so that polynomials up to degree 2n − 1 integrate exactly; a box uses
all n³ combinations (weights multiply). `ch05.cylinder_quadrature` does the same on a cylinder: Gauss–Legendre in R
(with the weight R) and z, equal steps in φ. Near a singular kernel like 1/r³ we replace r² by r² + ε² (a tiny smoothing
length) so no node ever divides by zero.
""", code=r"""
from fluidpy.core.integral_theorems import gauss_legendre_nodes     # 1-D nodes and weights on an interval
x, w = gauss_legendre_nodes(0.0, 1.0, 4)                             # 4 nodes on [0, 1]
X, Y, Z = np.meshgrid(x, x, x, indexing='ij'); W = w[:, None, None]*w[None, :, None]*w[None, None, :]   # the 4³ box rule
print(np.sum(W*X**3*Y**2*Z), 1/4*1/3*1/2)                            # exact for a polynomial: 0.041667 twice
""")
nb.md("**The sign test:** a straight Gaussian vortex tube (Γ = 1 m²/s, core σ_c = 0.1 m, length 4 m, flat ends normal to ω), velocity at 0.5 m from its axis in its mid-plane, by (5.14) with both signs and by (5.16).")
nb.code(r"""
F = ch05.gaussian_tube_fields(Gamma=1.0, sigma=0.1, L=4.0)           # ω and ∇×ω of the tube as functions; its bounding cylinder
nq = (32, 48, 40) if not FAST else (24, 32, 24)                       # Gauss nodes in R, φ, z
nodes, weights = ch05.cylinder_quadrature(F["bounds"]["radius"], F["bounds"]["z0"], F["bounds"]["z1"], nr=nq[0], nphi=nq[1], nz=nq[2])
xp = np.array([0.5, 0.0, 0.0])                                       # the field point [m]
for s in (+1, -1):                                                   # corrected (+1/4π) and printed (−1/4π) forms of (5.14)
    print(f"(5.14) with sign {s:+d}: u_y = {ch05.velocity_from_curl_omega(F['curl_omega'], xp, nodes, weights, sign=s)[1]:+.6f} m/s")
print(f"(5.16) Biot–Savart:    u_y = {ch05.biot_savart_volume(F['omega'], xp, nodes, weights)[1]:+.6f} m/s")
print(f"exact, this 4 m tube:  u_θ = {F['u_theta_reference'](0.5):.6f} m/s;  infinite line Γ/2πr = {F['u_theta_infinite'](0.5):.6f} m/s")
""", explain=r"""
1. A smooth vortex tube, its vorticity and the curl of its vorticity as callables, and Gauss nodes on a cylinder round it.
2. (5.14) with both signs: the corrected +1/(4π) gives +0.3088 m/s (counterclockwise, the way the vortex turns); the
   book's −1/(4π) gives −0.3088 m/s — the swirl backwards.
3. (5.16), Biot–Savart, gives the same +0.3088 m/s: D11 turns one into the other.
4. The exact value for this **finite** tube is 0.308838 m/s (the segment law of D13 applied to the core); an infinitely
   long line would give Γ(1 − e⁻²⁵)/2πr = 0.318310 m/s — the 3 % gap is the tube's finite length, not an error.
""")
D("D11", ref="5.16", check_src=r"""
xs = sp.symbols('x y z', real=True); xp = sp.symbols('xp yp zp', real=True)   # field point x and source point x′
d = sp.Matrix(xs) - sp.Matrix(xp); r = sp.sqrt(d.dot(d))             # x − x′ and its length
W = sp.Matrix([sp.Function(f'w{i}')(*xp) for i in range(3)])        # a generic vorticity field ω(x′)
phi = 1/r                                                            # step 1: the kernel φ(x′) = 1/|x − x′|
curlp = lambda A: sp.Matrix([sp.diff(A[2], xp[1]) - sp.diff(A[1], xp[2]), sp.diff(A[0], xp[2]) - sp.diff(A[2], xp[0]), sp.diff(A[1], xp[0]) - sp.diff(A[0], xp[1])])   # ∇′×
gradp = sp.Matrix([sp.diff(phi, v) for v in xp])                     # ∇′φ
print("step 2, ∇′×(φω) − φ∇′×ω − ∇′φ×ω:", sp.simplify(curlp(phi*W) - phi*curlp(W) - gradp.cross(W)).T)   # zero vector
print("step 4, ∇′(1/r) − (x − x′)/r³:  ", sp.simplify(gradp - d/r**3).T)                                  # zero vector
print("step 5, −∇′φ×ω − ω×(x − x′)/r³: ", sp.simplify(-gradp.cross(W) - W.cross(d)/r**3).T)                # zero vector
book = (d/r**3).cross(W)                                             # the book's printed "+((x − x′)/r³) × ω"
print("book's term − correct term − 2(x − x′)×ω/r³:", sp.simplify(book - (-gradp.cross(W)) - 2*d.cross(W)/r**3).T)   # zero: the book's term is off by 2(x − x′)×ω/r³
""")
note("N26", "The rewrite of the integrand", r"""
(D11 steps 2–5), with the correct sign. ⚠️ The book's second line prints
$+\big(\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}\big)\times\boldsymbol\omega$; it should be minus that
(= $+\boldsymbol\omega\times(\mathbf x-\mathbf x')/\lvert\mathbf x-\mathbf x'\rvert^3$). This second slip cancels the first, so (5.16) is right.
""", equation=r"\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}=\nabla'\times\Big(\frac{\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}\Big)+\boldsymbol\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}")
note("N27", "Gauss' theorem in curl form", r"""
used in D11 steps 7–9. Any smooth field F: $\int_V\nabla\times\mathbf F\,dV=\oint_A\mathbf n\times\mathbf F\,dA$ (apply
$\int_V\nabla\cdot\mathbf Q\,dV=\int_A\mathbf Q\cdot\mathbf n\,dA$ *(2.30)* to each component $\varepsilon_{kij}F_j$). The cell checks
it on a polynomial field.
""", equation=EQ["5.15"], ref="5.15")
nb.code(r"""
Fpoly = lambda X, Y, Z: np.array([X*Y**2, Y*Z + X**3, Z*X**2*Y])     # a polynomial vector field
cc6 = ch05.curl_theorem_box(Fpoly, [(0, 1), (0, 2), (-1, 1)], n=6)   # ∫∇×F dV vs ∮n×F dA on a box (6³ Gauss nodes)
print("volume:", cc6.volume.round(10), " surface:", cc6.surface.round(10), f" largest difference {np.max(np.abs(cc6.diff)):.1e}")   # equal to ~1e-13
""")
note("N28", "Why the surface term vanishes", r"""
(Fig. 5.8): V′ is a piece of tube with flat ends perpendicular to ω (n × ω = 0 there) and a side wall outside the vortex
(ω = 0 there). The cell below takes a box round a piece of our Gaussian tube (sides 4 core radii out, ends cutting ω at
right angles) and compares the surface term of (5.15) with the Biot–Savart integral over the same box.
""")
nb.code(r"""
kern = lambda X, Y, Z: F["omega"](np.array([X, Y, Z]), 0.0)/np.sqrt((X - 0.5)**2 + Y**2 + Z**2)   # ω/|x − x′| with x = (0.5, 0, 0)
cc = ch05.curl_theorem_box(kern, [(-0.4, 0.4), (-0.4, 0.4), (-0.5, 0.5)], n=24)   # the two sides of (5.15), 24³ Gauss nodes
print("surface term (5.15):", [f"{v:+.2e}" for v in cc.surface], " volume term:", [f"{v:+.2e}" for v in cc.volume])   # tiny: ~2e-6 in y
nb_, wb_ = ch05.cylinder_quadrature(0.4, -0.5, 0.5, nr=24, nphi=32, nz=24)   # the same piece of tube as a cylinder of radius 4σ_c
print("4π × Biot–Savart velocity from that piece:", [f"{v:+.4f}" for v in 4*np.pi*ch05.biot_savart_volume(F["omega"], np.array([0.5, 0, 0]), nb_, wb_)])   # ~2.85 in y
""", explain=r"""
The surface integral of (5.15) over the box is about 2×10⁻⁶ (from the side walls, where the Gaussian vorticity has
fallen to e⁻¹⁶ of its peak), against about 2.85 for the Biot–Savart integral of the same piece: dropping the surface
term, as D11 steps 10–11 do, costs nothing.
""")
note("N52", "Fig. 5.8 analogue", r"""
one element of a curved filament, the vector from it to the field point, and its push (below).
""")
nb.figure(r"""
phq = np.linspace(0, np.pi/2, 25)                                    # a quarter ring of radius 0.5 m in z = 0
poly = np.array([0.5*np.cos(phq), 0.5*np.sin(phq), 0*phq])           # the filament as a polyline (3, 25) [m]
fig = plt.figure(figsize=(6.2, 4.6)); fig.set_layout_engine("none")  # 3-D axes do not mix with constrained layout
ax = fig.add_subplot(projection="3d")                                # matplotlib 3-D axes (rotatable in a live kernel)
xf = np.array([0.6, 0.6, 0.4])                                       # the field point [m]
contrib = ch05.filament_contributions(xf, poly, 1.0, closed=False)    # (5.17) for each of the 24 segments, (3, 24) [m/s]
k = contrib.shape[1]//2; du_k, u_tot = contrib[:, k], contrib.sum(axis=1)   # the middle element's push and the total
ax.plot(*poly, color=COLORS["teal"], lw=2.5)                          # the filament
ax.plot(*poly[:, k:k + 2], color=COLORS["orange"], lw=6)              # the chosen element, thick orange
mid = 0.5*(poly[:, k] + poly[:, k + 1])                               # its midpoint x′
ax.plot(*np.c_[mid, xf], ls="--", color=COLORS["muted"])              # the vector x − x′
sc_ = 0.5/np.linalg.norm(u_tot)                                       # ONE scale for both arrows [m per m/s]
ax.quiver(*(xf + np.array([0, 0, 0.06])), *(du_k*sc_*5), color=COLORS["orange"], lw=2.5)   # the element's du, ×5 so it can be seen, drawn 6 cm above x (orange)
ax.quiver(*xf, *(u_tot*sc_), color=COLORS["accent"], lw=3.5)          # the total u, from x itself (purple)
ax.view_init(elev=18, azim=-20)                                       # a view in which both arrows are seen side-on
ax.scatter(*xf, color=COLORS["ink"], s=15)                            # the field point
ax.set_title(f"|du| = {np.linalg.norm(du_k):.4f} m/s is {100*np.linalg.norm(du_k)/np.linalg.norm(u_tot):.0f} % of |u| = {np.linalg.norm(u_tot):.3f} m/s", fontsize=9)
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")   # axis labels with units
savefig(fig, "ch05", "biot_savart_geometry"); plt.show()             # save the PNG to outputs/ch05, then draw
""", see="A teal quarter circle (the filament, Γ = 1 m²/s), a black field point above it, a grey dashed line from the "
         "middle element (thick orange) to the point, a purple arrow (the sum of all the elements' pushes) and an orange "
         "arrow just above it (the middle element's push, drawn five times longer than its true share, which the title gives).",
    read=r"Each element pushes perpendicular to itself and to the line joining it to x, with strength ∝ 1/distance² — "
         r"$\mathbf u=\frac1{4\pi}\int\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ (5.16); "
         "the purple arrow is their vector sum; the orange push points almost the same way because every element of this "
         "quarter ring sees the point from a similar side. (In D11 the volume V′ is a thin "
         "tube hugging such a filament.)",
    change="…the field point moved onto the axis of a full ring: all the pushes would line up along the axis.")
nb.worked_example("a vortex's velocity from its vorticity, by hand", r"""
A Gaussian vortex with Γ = 1 m²/s and core σ_c = 0.1 m ($\omega_z=(\Gamma/\pi\sigma_c^2)e^{-r^2/\sigma_c^2}$, peak 31.8 s⁻¹).

1. Outside the core, r = 0.5 m, nearly all the vorticity is inside: the planar Biot–Savart sum must give
   $\Gamma(1-e^{-25})/(2\pi\times0.5)$ = 0.3183 m/s.
2. At the core edge r = 0.1 m: $\Gamma(1-e^{-1})/(2\pi\times0.1)$ = 0.6321/0.6283 = 1.006 m/s.
3. Inside, r = 0.05 m: $\Gamma(1-e^{-0.25})/(2\pi\times0.05)$ = 0.2212/0.3142 = 0.704 m/s.

The vorticity *outside* a circle adds nothing there — the planar version of Newton's shell theorem.
""")
nb.code(r"""
xg = np.linspace(-0.5, 0.5, 201); XG, YG = np.meshgrid(xg, xg)       # a 201 × 201 grid, spacing 5 mm [m]
wz = vortices.gaussian_vortex(np.hypot(XG, YG), 1.0, 0.1)[1]         # ω_z of the Gaussian vortex on the grid [1/s]
for r in (0.05, 0.1, 0.5):                                           # three radii [m]
    ub = ch05.biot_savart_2d(wz, XG, YG, np.array([r, 0.0]), eps=1e-4)[1]   # planar Biot–Savart: Σ ω dA e_z×r/(2π|r|²) [m/s]
    ue = vortices.gaussian_vortex(r, 1.0, 0.1)[0]                    # the exact profile (3.29) [m/s]
    print(f"r = {r:.2f} m: grid sum {ub:.6f}, exact {ue:.6f}, difference {ub - ue:+.1e} m/s")
""", explain=r"""
1. The vorticity sampled on a 201 × 201 grid.
2. The planar Biot–Savart sum (kernel $\mathbf e_z\times\mathbf r/2\pi\lvert\mathbf r\rvert^2$, the 2-D form of (5.16)) at three radii.
3. Agreement with the exact profile $u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)$ *(3.29)*: within about 5×10⁻⁴ m/s inside the
   core (the kernel's 1/r singularity next to the grid nodes: 5×10⁻⁴ m/s at r = 5 and 10 cm) and ~10⁻⁸ m/s outside, where ω ≈ 0.
""")
nb.md("**From scratch — forty thousand little pushes, added by hand** (a double loop over the grid cells, `np.ndindex` walks every (row, column) index):")
nb.check_agree(r"""
dA = (xg[1] - xg[0])**2                                              # cell area [m²]
u = v = 0.0                                                          # velocity at the field point (0.1 m, 0) [m/s]
for j, i in np.ndindex(XG.shape):                                    # every grid cell
    dx, dy = 0.1 - XG[j, i], 0.0 - YG[j, i]                          # x − x′ [m]
    r2 = dx*dx + dy*dy + 1e-8                                        # |x − x′|² + ε² (ε = 1e-4 m) [m²]
    u += -wz[j, i]*dy/r2*dA/(2*np.pi)                                # e_z × (dx, dy) = (−dy, dx)
    v += wz[j, i]*dx/r2*dA/(2*np.pi)
print(f"ours ({u:.2e}, {v:.6f}) m/s")                                # (0, 1.005581)
assert np.allclose(v, ch05.biot_savart_2d(wz, XG, YG, np.array([0.1, 0.0]), eps=1e-4)[1], rtol=1e-6)   # = the library
""")
nb.plotly(r"""
rq = np.linspace(0.005, 0.5, 60)                                     # radii where we evaluate [m]
exact = vortices.gaussian_vortex(rq, 1.0, 0.1)[0]                    # (3.29) [m/s]


def grid(n):                                                         # grid points per side → curves
    g1 = np.linspace(-0.5, 0.5, int(n)); Xn, Yn = np.meshgrid(g1, g1)   # an n × n grid [m]
    w_ = vortices.gaussian_vortex(np.hypot(Xn, Yn), 1.0, 0.1)[1]     # vorticity on it [1/s]
    ub = ch05.biot_savart_2d(w_, Xn, Yn, np.array([rq, 0*rq]), eps=1e-4)[1]   # (5.16) summed over the grid [m/s]
    return {"Biot–Savart on the grid": (rq, ub), "exact (3.29)": (rq, exact), "printed sign −1/(4π)": (rq, -ub)}  # grid sum, exact profile, printed-sign mirror


fig = slider_figure(grid, "n", [21, 31, 41, 61, 81, 101, 121, 151, 181, 201] if not FAST else [21, 41, 61, 81, 101, 121],  # precompute every slider position
                    unit="points per side", xlabel="r [m]", ylabel="u_θ [m/s]",                     # slider label and axis labels
                    title="Biot–Savart rebuilds the vortex — with the right sign", modes={"Biot–Savart on the grid": "markers"})  # the message of the figure
recolor(fig, {"Biot–Savart on the grid": COLORS["teal"], "exact (3.29)": COLORS["ink"], "printed sign −1/(4π)": COLORS["rose"]},  # house colours by trace name
        dashes={"printed sign −1/(4π)": "dash"})                                                    # dashed reference curves
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="Teal dots settle onto the black curve as the grid refines (the coarse grids miss the core); the dashed rose "
        "curve is its upside-down twin.",
    read=r"(5.16) recovers the swirl inside and outside the core; the book's printed −1/(4π) in (5.14) would spin every "
         r"vortex backwards.",
    change="…a Rankine vortex instead: the dots converge more slowly (its vorticity jumps at the core edge), to the kink at r = σ_c.")
whatif(r"""
…the vortex were thin and we stood far from it? Then the details of the core stop mattering: only its strength Γ and
its path count — the filament law (C08).
""")

# ---- C08 ----------------------------------------------------------------------------------------------------------
core("C08", "The filament law", r"""
A smoke ring pushes itself forward; a wing-tip vortex drags air down behind a plane. How much velocity does one short
piece of a thin vortex make?
""", eqs=("5.17",))
problem(r"""
Real vortices are often thin compared with the distances that matter: the trailing vortices of an aircraft, a tornado
seen from a kilometre away, a smoke ring. Then we can replace the tube by a line carrying a circulation Γ and add up its
pieces — the lifting-line theory of wings (Ch. 14) is exactly this.
""")
idea("""
thin tube of strength Gamma, seen from far away          one piece of length dl pushes the fluid at x:
--->=======================>  (all vorticity on a line)   du = (Gamma dl / 4pi) e_w x (x - x')/|x - x'|^3        (5.17)
                                                          sum along a straight line -> Gamma/(4 pi d) (cos ta - cos tb)
                                                          infinite line -> Gamma/(2 pi d): the line vortex (5.2) again
""")
nb.md(r"""
*Glosses used below.* **Far field ("frozen kernel").** If a function barely changes across a small region, pull it out
of the integral: moving across a core of radius a changes the kernel $(\mathbf x-\mathbf x')/\lvert\mathbf x-\mathbf x'\rvert^3$,
whose size is $1/\lvert\mathbf x-\mathbf x'\rvert^2$, by about $2a/\lvert\mathbf x-\mathbf x'\rvert$ — 20 % at ten core radii, 2 % at a hundred. **1 + cot²θ = csc²θ** (divide sin²θ + cos²θ = 1 by sin²θ).
""")
D("D12", ref="5.17")
P("P144", "improper integral as a limit", r"""
An integral to infinity means: integrate to a finite end, then let the end run away. $\int_0^\infty e^{-x}dx$ is the
limit of $1-e^{-b}$ as b → ∞, i.e. 1. D13 lets a straight segment grow into an infinite line this way.
""", code=r"""
from scipy.integrate import quad                            # adaptive quadrature (Ch. 3 P87)
f = lambda l: 1.0/(1.0 + l**2)**1.5                         # d = 1 m: the segment integrand of D13 step 4
for b in (1, 10, 100):                                      # longer and longer segments [m]
    print(b, quad(f, -b, b)[0])                             # → 2 as b → ∞ (so Γ·2/(4π) = Γ/2π)
print(quad(f, -np.inf, np.inf)[0])                          # scipy does the limit for you: 2.0
""")
remind([
    ("substitution in an integral", r"replace the variable l by a function of a new one, l = l(θ), and dl by l′(θ)dθ (Ch. 3 P106)."),
], lead="needed in the derivation below")
D("D13", ref="5.2")
note("N29", "The loop closes", r"""
An infinite straight filament from (5.17) gives back $u_\theta=\Gamma/2\pi r$ (5.2); a finite segment gives
$\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ — the book never writes this. **Numbers** (Γ = 1 m²/s, d = 1 m): infinite
line 0.159 m/s, semi-infinite 0.080 m/s, the segment from −1 to +1 m 0.113 m/s; a square loop of side 1 m at its centre
2√2Γ/πL = 0.900 m/s; a ring of radius 0.5 m at its centre Γ/2R = 1.0 m/s, and on its axis the exact
$u_z=\frac{\Gamma R^2}{2(R^2+z^2)^{3/2}}$.
""")
nb.worked_example("the push of a 2 m segment on a point 1 m away", r"""
Γ = 1 m²/s, segment from z = −1 m to z = +1 m, field point at x = 1 m, z = 0 (d = 1 m).

1. The ends are seen at θ_a = 45° and θ_b = 135° (cos = ±0.707).
2. $u=\frac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ = 0.0796 × 1.414 = 0.1125 m/s.
3. Direction $\mathbf e_\omega\times\mathbf e_d=\mathbf e_z\times\mathbf e_x=\mathbf e_y$.
4. Stretch the segment to ±10 m: cosines ±0.995, u = 0.0796 × 1.990 = 0.1584 m/s — within 0.5 % of the infinite line's
   0.1592 m/s.
""")
nb.code(r"""
from scipy.integrate import quad                                     # for D13's step-4 integral
print(ch05.segment_induced_velocity(np.array([1.0, 0, 0]), np.array([0, 0, -1.0]), np.array([0, 0, 1.0]), 1.0))   # the example [m/s]
print(f"D13 step 4 by quad: {1.0*1.0/(4*np.pi)*quad(lambda l: 1/(1 + l**2)**1.5, -1, 1)[0]:.6f} m/s")   # Γd/4π ∫dl/(d²+l²)^{3/2}
print(f"semi-infinite line (0 to 1000 m): {ch05.segment_induced_velocity(np.array([1.0, 0, 0]), np.zeros(3), np.array([0, 0, 1000.0]), 1.0)[1]:.6f} m/s")   # Γ/4πd = 0.079577
for Lh in (1, 10, 100, 1000):                                        # half-lengths [m]
    print(Lh, round(ch05.segment_induced_velocity(np.array([1.0, 0, 0]), np.array([0, 0, -Lh]), np.array([0, 0, Lh]), 1.0)[1], 6))   # → 1/2π
print(f"square loop, side 1 m, at its centre: {ch05.filament_velocity_preset('square', 0.0, 0.0, 0.0, side=1.0, component=2):.6f} m/s")
for M in (8, 16, 64):                                                # a ring of radius 0.5 m drawn as an M-gon
    print(M, round(ch05.filament_velocity(np.zeros(3), ch05.filament_preset("ring", M=M, R=0.5), 1.0)[2], 6))   # → Γ/2R = 1
print(ch05.ring_axis_velocity(0.0, 0.5, 1.0), ch05.ring_axis_velocity(0.5, 0.5, 1.0))   # exact ring-axis formula at z = 0 and 0.5 m
""", explain=r"""
1. The closed form of D13 for the worked example: (0, 0.112540, 0) m/s — and `quad` of D13's step-4 integral gives the same.
2. Longer segments approach Γ/2πd = 0.159155 m/s (0.112540, 0.158365, 0.159147, 0.159155).
3. A square loop (four segments) at its centre: 0.900316 m/s = 2√2Γ/πL.
4. A ring as an M-gon at its centre: 1.054786, 1.013052, 1.000804 — the error falls as 1/M² (8 → 16 divides it by ≈ 4).
5. The exact ring-axis formula $\Gamma R^2/2(R^2+z^2)^{3/2}$: 1.0 m/s at the centre, 0.353553 m/s at z = 0.5 m.
""")
nb.md("**From scratch — twenty thousand pieces of (5.17), summed:**")
nb.check_agree(r"""
N = 20000                                                            # pieces of the segment
l = (np.arange(N) + 0.5)/N*2 - 1                                     # piece midpoints on [−1, 1] m
rv = np.stack([1.0 - 0*l, 0*l, -l])                                  # x − x′ for each piece [m]
du = np.cross(np.array([0, 0, 1.0]), rv.T)/np.linalg.norm(rv, axis=0)[:, None]**3*(2/N)/(4*np.pi)   # (5.17) with Γ = 1, dl = 2/N
print(du.sum(axis=0))                                                # (0, 0.112540, 0) m/s
assert np.allclose(du.sum(axis=0), ch05.segment_induced_velocity(np.array([1.0, 0, 0]), np.array([0, 0, -1.0]), np.array([0, 0, 1.0]), 1.0), rtol=1e-8, atol=1e-12)
""")
nb.plotly(r"""
zq = np.linspace(-1.5, 1.5, 61)                                      # points on the ring's axis [m]
exact = ch05.ring_axis_velocity(zq, 0.5, 1.0)                        # ΓR²/2(R² + z²)^{3/2} [m/s]


def ring(M):                                                         # number of segments → curves
    uz = ch05.filament_velocity(np.array([0*zq, 0*zq, zq]), ch05.filament_preset("ring", M=int(M), R=0.5), 1.0)[2]   # (5.17) summed
    return {"M-gon, (5.17) summed": (zq, uz), "exact ring": (zq, exact)}                            # polygon sum and exact ring


fig = slider_figure(ring, "M", [4, 6, 8, 12, 16, 24, 32, 48, 64, 128] if not FAST else [4, 8, 16, 32, 64, 128],  # precompute every slider position
                    unit="segments", xlabel="z on the axis [m]", ylabel="u_z [m/s]",                # slider label and axis labels
                    title="(5.17) summed round a ring (R = 0.5 m, Γ = 1 m²/s)", modes={"M-gon, (5.17) summed": "markers"})  # the message of the figure
recolor(fig, {"M-gon, (5.17) summed": COLORS["teal"], "exact ring": COLORS["ink"]})                 # house colours by trace name
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="Teal dots above the black curve for a square (M = 4), closing onto it as M grows.",
    read="A ring of vorticity pushes fluid through its middle, strongest at the centre, Γ/2R = 1 m/s; the polygon's "
         "error falls like 1/M².",
    change="…the ring radius doubled: the centre speed halves (Γ/2R) and the curve widens.")
nb.figure(r"""
ell = np.geomspace(0.1, 100, 200)                                    # half-length of the segment [m]
d = 1.0                                                              # distance of the point from the line [m]
full = np.array([ch05.segment_induced_velocity(np.array([d, 0, 0]), np.array([0, 0, -L_]), np.array([0, 0, L_]), 1.0)[1] for L_ in ell])   # segment −ℓ…ℓ
semi = np.array([ch05.segment_induced_velocity(np.array([d, 0, 0]), np.array([0, 0, 0.0]), np.array([0, 0, L_]), 1.0)[1] for L_ in ell])  # 0…ℓ
fig, ax = plt.subplots(figsize=(7, 3.6))                             # one panel
ax.semilogx(ell, full, color=COLORS["teal"], lw=2, label="segment from −ℓ to ℓ")       # climbs to the line vortex
ax.semilogx(ell, semi, color=COLORS["teal"], ls=":", lw=2, label="segment from 0 to ℓ")  # climbs to half of it
ax.axhline(1/(2*np.pi*d), color=COLORS["muted"], ls="--", label="infinite line Γ/2πd")    # the ghost
ax.axhline(1/(4*np.pi*d), color=COLORS["muted"], ls=":", label="semi-infinite Γ/4πd")
ax.set_xlabel("half-length ℓ [m]"); ax.set_ylabel("speed at d = 1 m [m/s]"); ax.legend(fontsize=8)   # axis labels with units
ax.set_title("A straight filament: far pieces matter less and less")
plt.show()
""", see="The solid teal curve climbs to the grey dashed line once ℓ ≫ d; the dotted curve (a line starting at the foot "
         "of the perpendicular) climbs to exactly half of it.",
    read="Γ/4πd (cos θ_a − cos θ_b): the far pieces are seen at angles near 0 and π and add little. A semi-infinite line "
         "gives exactly half — the reason a wing's trailing vortex makes half the downwash at the wing than far behind it (Ch. 14).",
    change="…d halved: the whole curve doubles and reaches the ghost at half the length.")
explainer("biot_savart_filament", "How does a vortex push water far away?",
          "dragging a field point around straight, bent, square, ring and helical filaments while each segment's "
          "contribution and their vector sum are drawn makes the 1/distance² law and the perpendicular direction tangible.",
          r"The sign toggle shows what the book's printed form $\mathbf u=-\frac1{4\pi}\int_{V'}\frac{\nabla'\times\boldsymbol\omega}{\lvert\mathbf x-\mathbf x'\rvert}d^3x'$ would do; the Derivation tab carries the ★★★ D10 and D11.", [
              "Straight segment: lengthen it and watch the speed approach Γ/2πd.",
              "Ring: put the point on the axis and read Γ/2R at the centre.",
              "Click a segment: the inspector shows its dl, r and du.",
              "Flip the sign to the printed −1/(4π): every arrow reverses — the badge says why that is wrong.",
          ])
whatif(r"""
…the filament were curved and we asked for the velocity *on* it? The sum diverges logarithmically near the point
itself: a curved vortex moves itself only because its core has a finite size — the smoke ring's self-speed
$\frac{\Gamma}{4\pi R}\big[\ln\frac{8R}{a}-\frac14\big]$ (note N42, C13).
""")

# =====================================================================================================================
# A.6 §5.6 Vorticity Equation in a Rotating Frame — R12–R17, C09, C10, R18, C11
# =====================================================================================================================
nb.section("5.6", "Vorticity Equation in a Rotating Frame", intro=r"""
**What is this section about?** The ocean and the atmosphere are watched from a rotating Earth and have density
contrasts. Redoing the vorticity equation in a frame turning at Ω, for a nearly incompressible (Boussinesq) fluid of
variable density, adds two things: the planet's own vorticity 2Ω is stretched and tilted like any other, and density
gradients across pressure gradients create vorticity. Two consequences follow: stretching and tilting of vortex lines
(the engine of 3-D turbulence), and the conservation of absolute circulation — the seed of potential vorticity.
""")
nb.recap("R12", "Hypotheses and comma notation of §5.6", r"""
A frame rotating at constant Ω; a non-barotropic but nearly incompressible (Boussinesq) fluid, so ∇·u ≈ 0; and the
comma notation of Ch. 2 §2.14: $u_{i,j}\equiv\partial u_i/\partial x_j$, $u_{i,jk}\equiv\partial^2u_i/\partial x_j\partial x_k$.
Gloss: at the end of a derivation we may rename the free index (n → i); a repeated (dummy) index may be renamed at any
time (Ch. 2 §2.1).
""", where="Ch. 4 §4.9 (Boussinesq), Ch. 2 §2.14")
nb.recap("R13", "Continuity", r"""
Nearly incompressible flow: $u_{i,i}=0$ *(Eq. 5.19)* — Ch. 4's $\nabla\cdot\mathbf u=0$ *(4.10)*.
""", where="Ch. 4 §4.2")
nb.recap("R14", "Momentum in a steadily rotating frame", r"""
$\frac{\partial u_i}{\partial t}+u_ju_{i,j}+2\varepsilon_{ijk}\Omega_ju_k=-\frac1\rho p_{,i}+g_i+\nu u_{i,jj}$ *(Eq. 5.20)* —
Ch. 4's rotating-frame momentum equation (4.45) with Ω constant; g is the effective gravity (true gravity plus the
centrifugal term). Ch. 4's explainer `rotating_frame_coriolis` shows the frame.
""", where="Ch. 4 §4.7")
nb.code(r"""
u_rest = lambda x, t: np.zeros_like(np.asarray(x, float))            # fluid at rest in a frame turning at 0.5 rad/s [m/s]
p_h = lambda x, t: 1e5 - 1000*9.81*np.asarray(x)[2]                  # hydrostatic pressure with effective gravity [Pa]
print(ch05.rotating_ns_residual(u_rest, p_h, 1000.0, np.array([0.3, 0.1, -0.2]), Omega=(0, 0, 0.5)))   # ≈ 0: a solution of (5.20)
""")
nb.recap("R15", "The Lamb identity in index form", r"""
$u_ju_{i,j}=-(\mathbf u\times\boldsymbol\omega)_i+\tfrac12(u_j^2)_{,i}$ *(Eq. 5.21)* — the index form of (4.68), recapped as R11.
""", where="Ch. 4 §4.9")
remind([
    ("np.random.default_rng", "`np.random.default_rng(seed).normal(size=…)` draws reproducible random numbers (Ch. 1 P10)."),
    ("np.einsum index strings", "`np.einsum('ijk,k->ij', E, w)` sums over the repeated index k exactly as the index notation says (Ch. 2 P62)."),
])
nb.recap("R16", "The antisymmetric gradient is the vorticity", r"""
$\varepsilon_{ijk}\omega_k=u_{j,i}-u_{i,j}$ *(Eq. 5.22)* — Ch. 3's rotation tensor $R_{ij}=-\varepsilon_{ijk}\omega_k$ *(3.15)*,
from the ε–δ identity (2.19).
""", where="Ch. 3 §3.4")
nb.code(r"""
G = np.random.default_rng(5).normal(size=(3, 3))                     # a random velocity gradient, G[i, j] = u_{i,j} [1/s]
w = kinematics.vorticity_from_gradient(G)                            # its vorticity vector [1/s]
E = tensors.levi_civita()                                            # ε_ijk as a 3×3×3 array
print(np.max(abs(np.einsum('ijk,k->ij', E, w) - (G.T - G))))         # ε_ijk ω_k − (u_{j,i} − u_{i,j}): 0 (round-off)
""")
nb.recap("R17", "The viscous force as a curl of the vorticity", r"""
$\nu u_{i,jj}=-\nu\varepsilon_{ijk}\omega_{k,j}$ *(Eq. 5.23)*: ν∇²u = −ν∇×ω for ∇·u = 0 — Ch. 4's (4.40).
""", where="Ch. 4 §4.6")

# ---- C09 ----------------------------------------------------------------------------------------------------------
core("C09", "The full vorticity equation: rotation and baroclinic generation", r"""
The air over a warm coast on a rotating planet: where does its spin come from — the winds, the planet, or the warm land?
""", eqs=("5.30",))
problem(r"""
Weather maps are maps of vorticity: cyclones, troughs, fronts. The equation that governs it for the atmosphere and ocean
must include the planet's rotation and the fact that warm and cold, fresh and salty fluid lie side by side. It is the
same idea as $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13) with two
new sources — and every Ch. 13 theory (quasi-geostrophy, potential vorticity, baroclinic instability) starts from it.
""")
idea("""
(5.13)  Dw/Dt = (w . grad)u                                   + nu lap(w)      constant rho, inertial frame
(5.30)  Dw/Dt = ((w + 2 Omega) . grad)u + grad(rho) x grad(p)/rho^2 + nu lap(w)    Boussinesq, frame turning at Omega
                 ^ the planet's vorticity 2 Omega    ^ baroclinic torque (C04)
                   is stretched and tilted like the fluid's own
""", r"""
Here **2Ω is the planetary vorticity**: fluid at rest in the rotating frame spins at Ω as seen from the stars, and
vorticity is twice the spin (recapped in full as R18 before C11).
""")
note("N30", "∇·ω = 0 in index form", r"""
$\omega_i=\varepsilon_{inq}u_{q,n}$, so (5.18) follows: $\varepsilon_{inq}$ is antisymmetric in i, n and $u_{q,ni}$ symmetric
(Schwarz), and the contraction of a symmetric with an antisymmetric pair vanishes (Ch. 2 §2.7). True for compressible
and unsteady flow too.
""", equation=r"\omega_{i,i}=(\varepsilon_{inq}u_{q,n})_{,i}=\varepsilon_{inq}u_{q,ni}=0\quad\text{or}\quad\nabla\cdot\boldsymbol\omega=0", ref="5.18")
nb.code(r"""
print(ch05.vorticity_divergence(ch05.abc_flow(), np.array([0.3, 1.1, -0.4])))   # ∇·(∇×u) of the ABC flow by nested stencils: ≈ 0
""")
nb.md(r"""
*Gloss — quotient rule for a gradient:* $\nabla(1/\rho)=-\nabla\rho/\rho^2$, so
$\nabla\times(\nabla p/\rho)=\nabla(1/\rho)\times\nabla p=-\nabla\rho\times\nabla p/\rho^2$ (the curl of the gradient ∇p itself is zero).
""")
D("D14", ref="5.25")
note("N31", "The Coriolis term reordered", r"""
D14's step 6 is (5.24), the Coriolis term with its dummy indices swapped (Ω × u = −u × Ω). Check with numbers:
Ω = (0, 0, 1), u = (1, 0, 0): both sides give (0, 2, 0).
""", equation=EQ["5.24"], ref="5.24")
nb.code(r"""
print(2*np.cross([0, 0, 1.0], [1.0, 0, 0]), -2*np.cross([1.0, 0, 0], [0, 0, 1.0]))   # 2Ω×u and −2u×Ω: [0 2 0] twice
""")
note("N32", "Rotating Navier–Stokes in Lamb form", r"""
D14's result is (5.25). The cell checks it on an **inertial oscillation** — fluid sliding in circles at the frequency
f = 2Ω, $\mathbf u=U(\cos ft,-\sin ft,0)$, an exact solution of (5.20) with hydrostatic pressure: the residual of (5.25)
equals the residual of (5.20).
""", equation=EQ["5.25"], ref="5.25")
nb.code(r"""
U0, Om = 0.1, 0.5                                                    # speed [m/s], frame rotation [rad/s]; f = 2Ω = 1 1/s
u_io = lambda x, t: np.stack([U0*np.cos(2*Om*t) + 0*np.asarray(x)[0], -U0*np.sin(2*Om*t) + 0*np.asarray(x)[0], 0*np.asarray(x)[0]])   # inertial oscillation
Phi = lambda x, t: 9.81*np.asarray(x)[2]                             # effective-gravity potential Φ = gz [J/kg]
x0 = np.array([0.2, -0.1, 0.3])                                      # a point [m]
lf = ch05.rotating_lamb_form_terms(u_io, p_h, 1000.0, Phi, x0, t=1.0, Omega=(0, 0, Om))   # the five terms of (5.25) [m/s²]
print("(5.25) residual:", [f"{v:.2e}" for v in lf["residual"]])     # ≈ 1e-10: stencil round-off
print("(5.20) residual:", [f"{v:.2e}" for v in ch05.rotating_ns_residual(u_io, p_h, 1000.0, x0, t=1.0, nu=1e-6, Omega=(0, 0, Om))])   # the same numbers
""")
remind([
    ("permutations, cyclic order and parity of ε", r"$\varepsilon_{ijk}$ is unchanged by a cyclic shift (ijk → jki) and changes sign when two indices are swapped (Ch. 2 P72)."),
    ("product rule for differentials", "d(ab) = a db + b da, also for derivatives with a comma (Ch. 1 P38)."),
], lead="needed in the derivation below")
D("D15", ref="5.30", check_src=r"""
X = sp.symbols('x0:3'); t, nu = sp.symbols('t nu'); Om = sp.symbols('Omega0:3')   # coordinates, time, ν, a constant Ω
u = [sp.Function(f'u{i}')(*X, t) for i in range(3)]                  # a generic velocity (NOT assumed divergence-free)
rho = sp.Function('rho')(*X, t); p = sp.Function('p')(*X, t)         # generic density and pressure
Phi = sp.Function('Phi')(*X)                                         # generic effective-gravity potential
d = lambda f, j: sp.diff(f, X[j])                                    # the comma: f_{,j}
eps = sp.LeviCivita                                                  # ε_ijk
curl = lambda V: [sum(eps(n, q, i)*d(V[i], q) for q in range(3) for i in range(3)) for n in range(3)]   # step 1: ε_nqi( )_{,q}
w = curl(u)                                                          # ω_n = ε_nqi u_{i,q}
gradB = [d(sum(uj**2 for uj in u)/2 + Phi, i) for i in range(3)]    # (½u_j² + Φ)_{,i}
lamb = [-sum(eps(i, j, k)*u[j]*(w[k] + 2*Om[k]) for j in range(3) for k in range(3)) for i in range(3)]   # −ε_ijk u_j(ω_k + 2Ω_k)
pres = [-d(p, i)/rho for i in range(3)]                              # −(1/ρ)p_{,i}
visc = [-nu*sum(eps(i, j, k)*d(w[k], j) for j in range(3) for k in range(3)) for i in range(3)]   # −ν ε_ijk ω_{k,j}
print("step 3, curl of the gradient term:", [sp.simplify(c) for c in curl(gradB)])   # [0, 0, 0]
book = [sum(-d(u[n], j)*(w[j] + 2*Om[j]) + u[j]*d(w[n], j) + d(u[j], j)*(w[n] + 2*Om[n]) - u[n]*d(w[j], j) for j in range(3)) for n in range(3)]   # steps 5–8, nothing dropped
print("steps 5–8, curl of the Lamb term minus the expanded form:", [sp.expand(a - b) for a, b in zip(curl(lamb), book)])   # [0, 0, 0]
print("step 7, ∇·ω for ANY u:", sp.expand(sum(d(w[j], j) for j in range(3))))   # 0: so u_n ω_{k,k} drops
gr = sp.Matrix([d(rho, j) for j in range(3)]); gp = sp.Matrix([d(p, j) for j in range(3)])   # ∇ρ, ∇p
baro = gr.cross(gp)/rho**2                                           # the baroclinic vector (5.28)
print("steps 10–11, curl of the pressure term minus ∇ρ×∇p/ρ²:", [sp.simplify(a - b) for a, b in zip(curl(pres), baro)])   # [0, 0, 0]
lapw = [sum(d(d(w[n], j), j) for j in range(3)) for n in range(3)]   # ω_{n,jj}
print("steps 12–13, curl of the viscous term minus ν∇²ω:", [sp.expand(a - nu*b) for a, b in zip(curl(visc), lapw)])   # [0, 0, 0]
res25 = [sp.diff(u[i], t) + gradB[i] + lamb[i] - pres[i] - visc[i] for i in range(3)]   # residual of (5.25), as LHS − RHS
res30 = [sp.diff(w[n], t) + sum(u[j]*d(w[n], j) for j in range(3)) - sum((w[j] + 2*Om[j])*d(u[n], j) for j in range(3))
         - baro[n] - nu*lapw[n] for n in range(3)]                   # residual of (5.30), as LHS − RHS
divu = sum(d(u[j], j) for j in range(3))                             # u_{j,j}: the term the book drops
print("steps 14–15, curl(5.25) − (5.30) − u_{j,j}(ω_n + 2Ω_n):",
      [sp.simplify(sp.expand(a - b - divu*(w[n] + 2*Om[n]))) for n, (a, b) in enumerate(zip(curl(res25), res30))])   # [0, 0, 0]
""", extra_check=r"The last printed line shows exactly what step 9 claims: the curl of (5.25) differs from (5.30) only by $u_{j,j}(\omega_n+2\Omega_n)$, which vanishes when $u_{i,i}=0$ (5.19).")
note("N33", "The curl of (5.25) in index form", r"""
D15's step 1 is (5.26). ⚠️ The book's sentence after it calls the scalar "Π"; it is Φ.
""", equation=EQ["5.26"], ref="5.26")
note("N34", "The Lamb and Coriolis terms", r"""
D15's steps 5–9 are (5.27). ⚠️ Between its second and third lines the book drops $u_{j,j}(\omega_n+2\Omega_n)$ without
comment; it is zero because $u_{i,i}=0$ (5.19).
""", equation=EQ["5.27"], ref="5.27")
note("N35", "The baroclinic vector", r"""
D15's step 11 is (5.28). It is exactly the torque of C04 (D06). Zero when ρ = ρ(p). The cell asks sympy for it in two
fields: density falling with x under hydrostatic pressure, and a density that is any function R of the pressure.
""", equation=r"-\varepsilon_{nqi}\Big(\frac1\rho p_{,i}\Big)_{,q}=-\frac1\rho\varepsilon_{nqi}p_{,iq}+\frac1{\rho^2}\varepsilon_{nqi}\rho_{,q}p_{,i}=0+\frac1{\rho^2}[\nabla\rho\times\nabla p]_n", ref="5.28")
nb.code(r"""
xs, ys, zs = sp.symbols('x y z')                                     # coordinates [m]
print(ch05.baroclinic_term_sym(1000 - 2*xs, 1e5 - 9810*zs, (xs, ys, zs)).T)   # ρ = 1000 − 2x, hydrostatic p: (0, −19620/(1000 − 2x)², 0)
print(ch05.baroclinic_term_sym(sp.Function('R')(1e5 - 9810*zs), 1e5 - 9810*zs, (xs, ys, zs)).T)   # ρ = R(p): (0, 0, 0)
""", explain=r"""
A sideways density gradient under vertical gravity gives a horizontal vorticity source — sympy writes it as
$-4905/(x-500)^2$, which is $-19620/(1000-2x)^2$: a roll about the y-axis, the sea-breeze circulation. A barotropic
field ρ = R(p) gives none.
""")
note("N36", "The viscous term", r"D15's steps 12–13 are (5.29):",
     equation=r"-\nu\varepsilon_{nqi}\varepsilon_{ijk}\omega_{k,jq}=-\nu(\delta_{nj}\delta_{qk}-\delta_{nk}\delta_{qj})\omega_{k,jq}=-\nu\omega_{k,nk}+\nu\omega_{n,jj}=\nu\omega_{n,jj}", ref="5.29")
nb.worked_example("a stretched column on the Earth and a lock exchange", r"""
1. **Planetary stretching.** At the pole, Ω = 7.292×10⁻⁵ rad/s; a column of air at rest stretched at
   $\partial w/\partial z$ = 10⁻⁵ s⁻¹ (a weak convergence: 1 cm/s of updraft gained per kilometre of height): the planetary term
   $2\Omega\,\partial w/\partial z$ = 2 × 7.292×10⁻⁵ × 10⁻⁵ = 1.46×10⁻⁹ s⁻² — kept up for a day (86 400 s) it gives
   ζ ≈ 1.3×10⁻⁴ s⁻¹, comparable with f itself (1.46×10⁻⁴ s⁻¹ at the pole): cyclones spin up this way.
2. **Baroclinic.** The lock exchange of C04: 2.42 s⁻² at the interface.
3. **Relative stretching** in Burgers' vortex (C06): 62 s⁻².

The same equation spans ten orders of magnitude.
""")
nb.code(r"""
for name in ("burgers", "lock_exchange", "rotating_column"):         # three scenes of (5.30), each at its preset probe
    b = ch05.vorticity_budget_preset(name)                           # z-components of every term [1/s²]
    print(f"{name:16s}" + "  ".join(f"{k} {v:+.4g}" for k, v in b.items()))
print(f"2Ω ∂w/∂z on the Earth: {2*rotating.OMEGA_EARTH*1e-5:.6e} 1/s²")   # the tiny example's planetary term
ts, nus, Oms = sp.symbols('t nu Omega')                               # symbols for the symbolic budget
vb = ch05.vorticity_budget_sym([ys*zs**2, xs*zs, -xs*ys], (xs, ys, zs), ts, nus, (0, 0, Oms), 1000 - 2*xs, 1e5 - 9810*zs)   # a divergence-free field
print("(5.30) with Ω = 0 and ρ = ρ(p) is (5.13):", vb["reduces_to_513"])   # True
""", explain=r"""
1. The budget of (5.30) at a probe of three scenes — one term dominates in each: Burgers (stretching +61.97, advective
   +15.49, diffusion −46.48, planetary and baroclinic 0); the lock exchange at the first instant (baroclinic 2.4222 s⁻²,
   and the local tendency it produces); a column at rest in a frame turning at Ω = 0.5 rad/s, stretched at α = 0.2 s⁻¹
   (planetary 2Ωα = 0.2 s⁻² = local).
2. The Earth number of the tiny example: 1.458423×10⁻⁹ s⁻².
3. sympy: with Ω → 0 and no baroclinic term, the residual of (5.30) becomes that of (5.13).
""")
nb.md("**From scratch — the baroclinic vector by central differences** at the lock-exchange interface:")
nb.check_agree(r"""
F = ch05.lock_exchange_fields(1000.0, 1025.0, 0.1)                   # ρ(x, y) and p(x, y) at t = 0⁺ (fluid at rest)
h, x0, y0 = 1e-4, 0.0, 0.5                                           # step [m]; the probe on the interface, mid-depth
drx = (F["rho"](x0 + h, y0) - F["rho"](x0 - h, y0))/(2*h); dry = (F["rho"](x0, y0 + h) - F["rho"](x0, y0 - h))/(2*h)   # ∇ρ [kg/m⁴]
dpx = (F["p"](x0 + h, y0) - F["p"](x0 - h, y0))/(2*h); dpy = (F["p"](x0, y0 + h) - F["p"](x0, y0 - h))/(2*h)           # ∇p [Pa/m]
b = (drx*dpy - dry*dpx)/F["rho"](x0, y0)**2                          # (∇ρ × ∇p)_z/ρ² [1/s²]
print(f"ours {b:.10f}, library {ch05.baroclinic_term(F['rho_fn'], F['p_fn'], np.array([x0, y0]))[2]:.10f} 1/s²")
assert np.allclose(b, ch05.baroclinic_term(F["rho_fn"], F["p_fn"], np.array([x0, y0]))[2], rtol=1e-8)   # same stencil, same number
""")
nb.figure(r"""
scenes = ["burgers", "lock_exchange", "rotating_column"]              # three regimes of (5.30)
titles = ["Burgers (Ω = 0, ρ const)", "lock exchange, t = 0⁺", "rotating column, Ω = 0.5 rad/s"]
keys = ["local", "advective", "stretching_tilting", "planetary", "baroclinic", "diffusion", "residual"]   # the bars of (5.30)
cols = [COLORS["blue"], COLORS["muted"], COLORS["accent"], COLORS["amber"], COLORS["orange"], COLORS["rose"], COLORS["ink"]]
labs = ["∂ω/∂t", "(u·∇)ω", "(ω·∇)u", "(2Ω·∇)u", "∇ρ×∇p/ρ²", "ν∇²ω", "residual"]
fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=True)      # one panel per scene
for a, nm, tt in zip(axs, scenes, titles):
    b = ch05.vorticity_budget_preset(nm)                              # z-components [1/s²]
    big = max(abs(b[k]) for k in keys)                                # the largest term of this scene
    a.bar(range(7), [b[k]/big for k in keys], color=cols)             # normalised by it
    a.axhline(0, color=COLORS["muted"], lw=0.8); a.set_xticks(range(7)); a.set_xticklabels(labs, rotation=45, fontsize=8)
    a.set_title(f"{tt}\n(largest term {big:.3g} 1/s²)", fontsize=9)
axs[0].set_ylabel("term / largest term [–]")                          # normalised, so the three scenes share an axis
fig.suptitle("Same equation, three regimes", fontsize=11)
savefig(fig, "ch05", "vorticity_budget_530"); plt.show()             # save the PNG to outputs/ch05, then draw
""", see="(Note `N37`: reading the terms of (5.30).) Burgers — grey advective bar balanced by purple stretching and a negative rose diffusion bar; lock exchange — "
         "one orange baroclinic bar and the blue local bar it drives; rotating column — the amber planetary bar equals the "
         "blue local bar. Every residual (black) is invisible.",
    read=r"The left side of (5.30) is the rate following a particle (blue + grey); $\nu\nabla^2\boldsymbol\omega$ is molecular "
         r"diffusion (as $\nu\nabla^2\mathbf u$ diffuses velocity); the baroclinic bar exists only where ∇ρ crosses ∇p; stretching "
         "matters even with Ω = 0. Each panel is divided by its own largest term — the titles give the sizes (62, 2.4 and "
         "0.2 s⁻²).",
    change="…Ω switched off in the rotating column: the amber bar vanishes and nothing spins the column up.")
nb.md(r"""
> ⚠️ **Common confusion:** "Boussinesq means the density can be treated as constant." In (5.30) ρ is kept inside the
> pressure term — that is where the baroclinic source lives. To leading order
> $\frac{1}{\rho^2}\nabla\rho\times\nabla p\approx\frac{1}{\rho_0}\nabla\rho'\times\mathbf g$ with g = −g e_z (cf. Ch. 4's
> Boussinesq momentum equation $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$
> *(4.86)*): only the density *perturbation* crossed with gravity matters.
""")
nb.md(r"""
🎮 **What does a spinning planet add to the vorticity equation?** (embedded after C11 below) walks through D14 and D15
with term bars in its Derivation tab.
""")
whatif(r"""
…we looked at one vortex line in an inviscid, barotropic, inertial flow? Only $(\boldsymbol\omega\cdot\nabla)\mathbf u$ would
remain — and it splits cleanly into stretching and tilting (C10).
""")

# ---- C10 ----------------------------------------------------------------------------------------------------------
core("C10", "Stretching and tilting of vortex lines", r"""
A figure skater spins faster by pulling in her arms; a tornado tightens as its column is stretched upward. How can a
flow make vorticity stronger without any torque at all?
""", eqs=("5.32",))
problem(r"""
In three dimensions a flow can amplify the vorticity it carries simply by stretching vortex lines, and redirect it by
tilting them. This is how turbulence builds its thin, intense vortices (Ch. 12) and why a flat, two-dimensional flow
behaves so differently. It is also, with the planet's vorticity, how a converging column of air becomes a cyclone (C11).
""")
idea("""
(w . grad)u = |w| du/ds        (how u changes ALONG the vortex line)                                     (5.31)
   |- along the line  (du_s/ds > 0): the line is stretched -> thinner -> spins faster   (the skater)
   '- across the line (du_n/ds, du_m/ds): the line is turned -> vorticity appears in new directions (tilting)   (5.32)
2-D flow: w perpendicular to the plane, nothing varies along it -> neither
""")
P("P145", "Frenet frame of a curve", r"""
At each point of a smooth curve: the unit tangent $\mathbf e_s$, the principal normal (toward the centre of curvature)
and the binormal (their cross product). The book's Fig. 5.9 takes $\mathbf e_n$ **away** from the centre of curvature
(the opposite of the usual Frenet normal) and $\mathbf e_m$ as the second normal; for a helix of radius a and pitch 2πc
the curvature is $a/(a^2+c^2)$ and the torsion $c/(a^2+c^2)$.
""", code=r"""
d = ch05.helix_frame(0.0, a=1.0, c=0.3)                  # the frame at s = 0 of a helix, radius 1 m, pitch 2π·0.3 m
print(np.round([d['e_s'], d['e_n'], d['e_m']], 4))       # three orthonormal vectors; e_n = (1, 0, 0) points away from the axis
print(round(d['curvature'], 6), round(d['torsion'], 6), np.dot(d['e_s'], d['e_n']))   # 0.917431 0.275229 0.0
""")
gloss(["e_n away from the centre of curvature (book convention)"])
P("P146", "angular momentum of a spinning cylinder", r"""
A cylinder of mass m and radius R spinning at Ω about its axis has angular momentum L = IΩ with I = ½mR² = mA/2π (A its
cross-section). With no torque L is fixed; stretch the cylinder at fixed mass and volume and A shrinks, so Ω grows like
1/A — like the length.
""", code=r"""
m, A0, Om0 = 1.0, 1e-4, 10.0                             # mass [kg], cross-section [m²], spin [rad/s]
for stretch in (1, 2, 4):                                # length × 1, 2, 4 at fixed volume
    A = A0/stretch; Om = Om0*A0/A                        # L = (m A/2π) Ω fixed
    print(stretch, A, Om)                                # the spin doubles when the length doubles
""")
remind([
    ("orthonormal basis, projection and completeness", "the component of a vector along a unit vector e is its dot product with e; three orthonormal vectors rebuild it (Ch. 2 P65)."),
], lead="needed in the derivations below")
D("D16", ref="5.31")
note("N38", "Natural coordinates", r"""
D16's result is (5.31). The stretching part is $(\mathbf e_s\cdot\mathsf G\,\mathbf e_s)\boldsymbol\omega$ along ω; the tilting
part is the rest of $\mathsf G\boldsymbol\omega$, perpendicular to ω.
""", equation=r"(\boldsymbol\omega\cdot\nabla)\mathbf u=\Big[\boldsymbol\omega\cdot\Big(\mathbf e_s\frac{\partial}{\partial s}+\mathbf e_n\frac{\partial}{\partial n}+\mathbf e_m\frac{\partial}{\partial m}\Big)\Big]\mathbf u=\omega\frac{\partial\mathbf u}{\partial s}", ref="5.31")
D("D17", ref="5.32")
nb.worked_example("stretching and tilting of one vortex line", r"""
ω = (1, 0, 1) s⁻¹ in the axial strain G = diag(−½, −½, 1) s⁻¹ (u = −x/2, v = −y/2, w = z).

1. $(\boldsymbol\omega\cdot\nabla)\mathbf u=\mathsf G\boldsymbol\omega$ = (−0.5, 0, 1) s⁻².
2. $\mathbf e_s=\boldsymbol\omega/\lvert\boldsymbol\omega\rvert$ = (1, 0, 1)/√2; stretching rate $\mathbf e_s\cdot\mathsf G\,\mathbf e_s$ = (−0.5 + 1)/2 = 0.25 s⁻¹.
3. Stretching part 0.25 × ω = (0.25, 0, 0.25) — along ω.
4. Tilting part (−0.5, 0, 1) − (0.25, 0, 0.25) = (−0.75, 0, 0.75) — perpendicular to ω (dot product −0.75 + 0.75 = 0):
   the line is being turned toward the stretching axis z.
""")
remind([
    ("scipy.linalg.expm and e^{Gt}ω₀", r"the matrix exponential solves $d\boldsymbol\omega/dt=\mathsf G\boldsymbol\omega$ for a constant G: $\boldsymbol\omega(t)=e^{\mathsf Gt}\boldsymbol\omega_0$ (Ch. 2 P79)."),
])
nb.code(r"""
S = ch05.stretching_tilting_split(np.array([1.0, 0, 1.0]), ch05.strain_preset("axial_stretch"))   # the tiny example
print(f"e_s {S['e_s'].round(4)}, rate {S['rate']:.4f} 1/s, stretching {S['stretching'].round(4)}, tilting {S['tilting'].round(4)}, total Gω {S['total'].round(4)}")
S2 = ch05.stretching_tilting_split(np.array([0, 0, 2.0]), ch05.strain_preset("axial_stretch"))   # ω along the stretching axis
print(f"ω along z: rate {S2['rate']:.4f} 1/s (= α), tilting {S2['tilting']}")                      # rate α = 1, tilting 0
for name, w0 in (("axial_stretch", [0, 0, 1]), ("shear_tilt", [1, 0, 0]), ("planar", [0, 0, 1])):   # three steady linear flows, rate 1 1/s
    print(name, np.round(ch05.uniform_strain_vorticity(w0, name, 2.0), 4))   # ω after 2 s = e^{2G} ω0 (inviscid, frozen-in)
print(ch05.stretched_tube(2.0, 1.0, 10.0, 1e-4))                      # a tube stretched from 1 m to 2 m: ω 10 → 20 1/s, area halves
""", explain=r"""
1. The split of the worked example: $\mathbf e_s$ = (0.7071, 0, 0.7071), rate 0.25 s⁻¹, stretching (0.25, 0, 0.25),
   tilting (−0.75, 0, 0.75), total Gω = (−0.5, 0, 1).
2. With ω along the stretching axis the whole of Gω is stretching (rate α = 1 s⁻¹) and the tilting is zero.
3. Three steady linear flows for 2 s: axial stretching multiplies ω_z by e² = 7.389; the shear w = x tilts ω_x into
   ω_z = s t ω_x = 2.0; the planar strain leaves ω ⟂ plane alone.
4. The stretched tube of D18: ω doubles, the area halves, Γ = ωA stays 10⁻³ m²/s.
""")
nb.md("**From scratch — the projection by hand:** the part of Gω along ω, and the rest.")
nb.check_agree(r"""
w = np.array([1.0, 0, 1.0]); G = ch05.strain_preset("axial_stretch")   # the tiny example's ω [1/s] and G [1/s]
es = w/np.linalg.norm(w)                                             # unit vector along the vortex line
rate = es @ G @ es                                                   # stretching rate e_s·G e_s [1/s]
st = rate*w                                                          # stretching part (along ω) [1/s²]
ti = G @ w - st                                                      # tilting part: the rest of Gω [1/s²]
print(rate, st, ti, ti @ w)                                          # 0.25, (0.25, 0, 0.25), (−0.75, 0, 0.75), 0
assert np.allclose(st, S["stretching"]) and np.allclose(ti, S["tilting"]) and abs(ti @ w) < 1e-14   # = the library, and ⟂ ω
""")
nb.plotly(r"""
a_h, c_h, s0 = 1.0, 0.3, 1.0                                         # helix radius [m], pitch/2π [m], arc length of the marked point [m]
sl = np.linspace(0, 3*2*np.pi*np.hypot(a_h, c_h), 400)               # three turns of arc length [m]
Xh = ch05.helix(sl, a_h, c_h)                                        # the vortex line (3, 400) [m]
P0 = ch05.helix(s0, a_h, c_h); fr = ch05.helix_frame(s0, a_h, c_h)   # the marked point and its natural frame
fig = go.Figure(go.Scatter3d(x=Xh[0], y=Xh[1], z=Xh[2], mode="lines", line=dict(color=COLORS["teal"], width=5), name="vortex line"))  # the helical vortex line
for key, col in (("e_s", COLORS["ink"]), ("e_n", COLORS["orange"]), ("e_m", COLORS["muted"])):   # the frame (book convention: e_n away from the axis)
    e = 0.6*fr[key]                                                                                 # a frame vector drawn 0.6 m long
    fig.add_trace(go.Scatter3d(x=[P0[0], P0[0] + e[0]], y=[P0[1], P0[1] + e[1]], z=[P0[2], P0[2] + e[2]], mode="lines+text",  # draw it from the marked point
                               text=["", key], line=dict(color=col, width=6), name=key))            # labelled at its tip
presets = ["axial_stretch", "shear_tilt", "planar"]                  # three velocity gradients G (rate 1 1/s)
for j, nm in enumerate(presets):                                     # the stretching and tilting arrows for each G, |ω| = 1 1/s
    sp_ = ch05.stretching_tilting_split(fr["e_s"], ch05.strain_preset(nm))                          # split Gω into stretching and tilting (|ω| = 1 1/s)
    for part, col in (("stretching", COLORS["accent"]), ("tilting", COLORS["blue"])):               # purple along the line, blue across it
        v = 0.8*sp_[part]                                                                           # the arrow, scaled for visibility
        fig.add_trace(go.Scatter3d(x=[P0[0], P0[0] + v[0]], y=[P0[1], P0[1] + v[1]], z=[P0[2], P0[2] + v[2]], mode="lines",  # drawn from the marked point
                                   line=dict(color=col, width=10), name=f"{part} ({nm})", visible=(j == 0)))  # only the first flow's arrows start visible
buttons = []                                                          # a dropdown: show the arrows of one G at a time
for j, nm in enumerate(presets):                                                                    # one dropdown button per flow
    vis = [True]*4 + [k//2 == j for k in range(6)]                   # line + 3 frame vectors always; 2 arrows per preset
    buttons.append(dict(label=nm, method="update", args=[{"visible": vis}]))                        # shows that flow's two arrows
fig.update_layout(updatemenus=[dict(buttons=buttons, x=0.02, y=0.98)], height=500, margin=dict(l=0, r=0, t=40, b=0),  # the dropdown, size and title
                  title=f"Natural frame on a helix (κ = {fr['curvature']:.4f} 1/m): stretching (purple) and tilting (blue)",  # the message of the figure
                  scene=dict(aspectmode="data", xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]"))  # equal axes, labels with units
fig.show()                                                                                          # draw it (works on the web page too)
""")
gloss(["plotly dropdown menus"])
see_read_change(
    see="(Note `N53`, Fig. 5.9 analogue.) A teal helical vortex line with, at one point, its frame (e_s black along the line, e_n orange pointing away from "
        "the axis, e_m grey) and a thick purple arrow along the line plus a blue one across it; the dropdown (a plotly "
        "`updatemenus` menu) switches between three flows.",
    read=r"$(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega\frac{\partial\mathbf u}{\partial s}$ (5.31): only the change of u along the line counts; "
         "its along-component stretches (purple), its across-components tilt (blue).",
    change="…the planar preset: G has no z-row or column, so for a line that is mostly vertical both arrows shrink — and "
           "for a line exactly along z they vanish (2-D flow).")
nb.md(r"""
**Stretching and tilting in motion.** Left, a vortex tube in axial strain α = 1 s⁻¹ seen side-on (it lengthens as
$e^{\alpha t}$ and thins as $e^{-\alpha t/2}$) and end-on with a spin marker; middle, a vortex segment in the shear
w = s x (s = 1 s⁻¹) being tilted from x toward z; right, the growth |ω|(t)/|ω₀| (purple) against the ghost $e^{\alpha t}$
and the tilt angle (blue). Use ◀ ▶ to step.
""")
nb.animation(r"""
nf = 19 if not FAST else 10                                          # frames
tt = np.linspace(0, 3*np.log(2), nf)                                 # 0 … 3 ln 2 s, so that frames land on t = ln 2 and 2 ln 2 [s]
alpha, s_sh, w0 = 1.0, 1.0, 10.0                                     # strain rate, shear rate [1/s], initial vorticity [1/s]
grow = np.array([ch05.uniform_strain_vorticity([0, 0, 1], "axial_stretch", t_)[2] for t_ in tt])   # |ω|/|ω0| in axial strain
tilt = np.array([ch05.uniform_strain_vorticity([1, 0, 0], "shear_tilt", t_) for t_ in tt])         # ω(t) of a segment in shear
with plt.rc_context(ANIM_RC):                                         # light settings for fast frames
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(10.5, 3.6), dpi=80)
fig.subplots_adjust(left=0.04, right=0.98, bottom=0.14, top=0.86, wspace=0.3)   # fixed margins
tube = plt.Rectangle((-0.5, -0.1), 1.0, 0.2, fc=COLORS["teal"], alpha=0.4, ec=COLORS["teal"]); a.add_patch(tube)   # side view
circ = plt.Circle((0, 0.75), 0.1, fc="none", ec=COLORS["teal"], lw=2); a.add_patch(circ)   # end view
(spoke,) = a.plot([0, 0.1], [0.75, 0.75], color=COLORS["ink"], lw=2)   # spin marker
a.set_xlim(-4, 4); a.set_ylim(-0.6, 1.1); a.set_aspect("equal"); a.axis("off"); a.set_title("axial strain: length × e^{αt}")
(seg,) = b.plot([0, 1], [0, 0], color=COLORS["blue"], lw=3)          # the vortex segment in the (x, z) plane
b.set_xlim(-0.2, 1.3); b.set_ylim(-0.2, 2.2); b.set_aspect("equal"); b.set_xlabel("ω_x [1/s]"); b.set_ylabel("ω_z [1/s]")
b.set_title("shear w = s x: ω_x tilted into ω_z")
c.plot(tt, np.exp(alpha*tt), color=COLORS["muted"], ls="--", label="ghost e^{αt}")   # the prediction
(gl,) = c.plot([], [], color=COLORS["accent"], lw=2.5, label="|ω|/|ω₀|, axial strain")
(tl,) = c.plot([], [], color=COLORS["blue"], lw=2, label="tilt angle [rad], shear")
c.set_xlim(0, tt[-1]); c.set_ylim(0, 8.5); c.set_xlabel("t [s]"); c.legend(fontsize=7, loc="upper left"); c.set_title("growth and tilt")
theta = (w0/(2*alpha))*(np.exp(alpha*tt) - 1)                        # angle turned by the spin marker: ∫ ω/2 dt [rad]


def update(i):                                                       # draw frame i
    L = np.exp(alpha*tt[i]); r = 0.1*np.exp(-alpha*tt[i]/2)          # length × e^{αt}, radius × e^{−αt/2}
    tube.set_bounds(-0.5*L, -r, L, 2*r)                              # the stretched, thinned tube
    circ.set_radius(r)                                               # its end view
    spoke.set_data([0, r*np.cos(theta[i])], [0.75, 0.75 + r*np.sin(theta[i])])   # the spin marker turns faster and faster
    seg.set_data([0, tilt[i][0]], [0, tilt[i][2]])                   # ω = (1, 0, s t): tilted toward z
    gl.set_data(tt[:i + 1], grow[:i + 1]); tl.set_data(tt[:i + 1], np.arctan2(tilt[:i + 1, 2], tilt[:i + 1, 0]))
    fig.suptitle(f"t = {tt[i]:.2f} s: length × {L:.2f}, |ω| × {grow[i]:.2f}", fontsize=10)
    return tube, circ, spoke, seg, gl, tl


with plt.rc_context(ANIM_RC):                                         # render the frames with the light settings
    show_animation(animate(update, frames=nf, fig=fig, interval=250), player="frames")   # step through with ◀ ▶
""")
see_read_change(
    see="The tube lengthens and thins while its end-view spoke turns faster; the blue segment swings from the x-axis "
        "toward z; the purple growth curve lies on the grey ghost and the blue tilt angle climbs toward π/2.",
    read="Step to t = ln 2 ≈ 0.69 s: length × 2, vorticity × 2 — D18's ω ∝ L; at 2 ln 2 ≈ 1.39 s both are × 4, at the "
         "last frame × 8. The tilted segment's ω_z = s t ω_x grows linearly; no torque acts in either case (D17 step 4).",
    change="…α negative (compression): the tube shortens, fattens and spins down.")
note("N21", "Burgers' vortex", r"""
(the book leaves it to Exercise 5.12): an axial strain $u_R=-\alpha R/2$, $u_z=\alpha z$ with a swirl, as below, is steady
because stretching (which would thin and intensify the core) is balanced by diffusion (which would spread it). D18
derives the balance. **Number:** α = 1 s⁻¹, water, Γ = 10⁻³ m²/s: core radius $\sqrt{4\nu/\alpha}$ = 2 mm, peak vorticity
79.6 s⁻¹. → the classic model of the finest vortices of turbulence (Ch. 12).
""", equation=r"u_\varphi=\frac{\Gamma}{2\pi R}\Big[1-e^{-\alpha R^2/4\nu}\Big],\qquad\omega_z=\frac{\alpha\Gamma}{4\pi\nu}e^{-\alpha R^2/4\nu}")
remind([
    ("separation of variables", "put everything with ω on one side and everything with R on the other, then integrate both (Ch. 1 P42)."),
    ("natural logarithm and exponential", r"$\int dy/y=\ln y$, so y′/y = k gives $y=y_0e^{kx}$ (Ch. 1 P36)."),
], lead="needed in the derivation below")
nb.md(r"""
*Glosses used in D18.* **Incompressible tube:** a material tube of fixed volume V = A L gets thinner as it gets longer.
**Cylindrical Laplacian of an axisymmetric scalar:** $\nabla^2\omega_z=\frac1R\frac{d}{dR}\big(R\frac{d\omega_z}{dR}\big)$ when ω_z depends on R only.
""")
D("D18", ref="5.13")
nb.code(r"""
R = np.linspace(0, 8e-3, 401)                                        # radius [m]
uR, uphi, uz, wz = ch05.burgers_vortex(R, 0.0, 1e-3, 1.0, 1e-6)       # Burgers: Γ = 1e-3 m²/s, α = 1 1/s, ν = 1e-6 m²/s
print(ch05.burgers_core_radius(1.0, 1e-6), round(wz[0], 4), round(np.trapezoid(wz*2*np.pi*R, R), 7))   # core 2 mm, peak 79.5775, Γ
bb = ch05.vorticity_budget_preset("burgers", x=2e-3)                  # the terms of (5.13) at R = core radius
print({k: round(v, 4) for k, v in bb.items()})                       # advective = stretching = 29.27, diffusion ≈ 0
""", explain=r"""
1. The closed-form Burgers vortex. 2. Core radius $\sqrt{4\nu/\alpha}$ = 0.002 m, peak ω_z = 79.5775 s⁻¹, and the total
circulation $\int\omega_z2\pi R\,dR$ = 0.0010000 m²/s (Γ recovered to about 2×10⁻⁵ on this grid, which stops at 4 core radii).
3. At R = √(4ν/α) the diffusion term changes sign (where the cylindrical Laplacian (1/R)(Rω_z′)′ changes sign): inside it removes
vorticity, outside it adds; there advection (+29.27 s⁻²) equals stretching (+29.27 s⁻²).
""")
nb.figure(r"""
Rm = np.linspace(1e-5, 8e-3, 300)                                    # radius [m] (the terms are regular at R → 0)
bal = ch05.burgers_balance(Rm)                                       # ω_z and the three terms of (5.13) [1/s], [1/s²]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6))                  # profile | balance
a.plot(Rm*1e3, bal["omega_z"], color=COLORS["teal"]); a.axvline(bal["core_radius"]*1e3, color=COLORS["muted"], ls="--")  # the vorticity profile and the core radius
a.set_xlabel("R [mm]"); a.set_ylabel(r"$\omega_z$ [1/s]"); a.set_title("(a) Burgers' vortex: a Gaussian core", fontsize=10)  # labels with units
b.plot(Rm*1e3, bal["advective"], color=COLORS["muted"], label=r"advection $u_R\,\partial\omega_z/\partial R$")  # grey: advection (left side of (5.13))
b.plot(Rm*1e3, bal["stretching"], color=COLORS["accent"], label=r"stretching $\omega_z\,\partial u_z/\partial z$")  # purple: stretching
b.plot(Rm*1e3, bal["diffusion"], color=COLORS["rose"], label=r"diffusion $\nu\nabla^2\omega_z$")    # rose: diffusion
b.plot(Rm*1e3, bal["residual"], color=COLORS["ink"], label="advection − stretching − diffusion")    # black: the residual, ≈ 0
b.axvline(bal["core_radius"]*1e3, color=COLORS["muted"], ls="--"); b.set_xlabel("R [mm]"); b.set_ylabel("term [1/s²]")  # the core radius; labels with units
b.legend(fontsize=8); b.set_title("(b) the steady balance of (5.13)", fontsize=10)                  # legend and title
savefig(fig, "ch05", "burgers_balance"); plt.show()                  # save the PNG to outputs/ch05, then draw
""", see="Left, a Gaussian vorticity profile with the dashed core radius 2 mm. Right, purple stretching largest on the "
         "axis, rose diffusion negative inside the core and positive outside, grey advection zero on the axis and peaking "
         "near the core, and a black residual flat at zero.",
    read="Advection brings vorticity inward, stretching amplifies it, diffusion spreads it outward; advection = "
         "stretching + diffusion everywhere, so the vortex is steady. The dashed line (√(4ν/α)) is where diffusion "
         "changes sign.",
    change="…α four times larger: the core halves (√(4ν/α)) and the peak quadruples (αΓ/4πν).")
explainer("vorticity_stretching_tilting", "How can a flow spin fluid faster without a torque?",
          "a vortex line and a material element ride a chosen flow on one clock while |ω| grows or turns, the stretching "
          "and tilting arrows update and the term bars of (5.13) fight — growth in time is what a static figure cannot show.",
          "Switching on ν brings in diffusion and the Burgers balance; the Derivation tab carries the ★★★ D09 and D16–D18.", [
              "Axial strain, press ▶: |ω| grows as e^{αt}; the material element stays parallel to ω.",
              "Shear tilt: a new component appears — tilting.",
              "Planar flow: nothing happens — 2-D.",
              "Burgers: raise ν and watch the core widen to √(4ν/α).",
          ])
nb.md(r"""
> ⚠️ **Common confusion:** "stretching needs a torque." It needs none: angular momentum is conserved while the moment of
> inertia shrinks — the spin rises on its own, exactly as the skater's does.
""")
whatif(r"""
…the frame rotated? Then a vertical fluid column contains the planet's vorticity 2Ω even when it does not spin relative
to the Earth — and stretching it spins it up (C11).
""")

# ---- R18 + C11 ----------------------------------------------------------------------------------------------------
remind([
    ("latitude, Earth's rotation rate and the local vertical", r"at latitude φ only the local vertical part of the Earth's rotation, Ω sin φ, turns horizontal motion: $f=2\Omega\sin\varphi$ (Ch. 4 P126)."),
])
nb.recap("R18", "Planetary and absolute vorticity", r"""
Seen from a frame rotating at Ω, a fluid at rest in the frame has inertial vorticity 2Ω — the *planetary vorticity*;
the *absolute vorticity* is ω + 2Ω (Ch. 3: ω′ = ω − 2Ω). On the Earth only the local vertical part matters for
horizontal motion: f = 2Ω sin φ, and the absolute vertical vorticity is ζ + f.
""", where="Ch. 3 §3.4, Ch. 4 §4.7")
nb.code(r"""
print([f"{v:.4e}" for v in ch05.absolute_vorticity(np.array([0, 0, 1e-5]), np.array([0, 0, rotating.OMEGA_EARTH]))])   # ω + 2Ω at the pole: (0, 0, 1.5584e-4) [1/s]
print([f"{rotating.coriolis_parameter(np.deg2rad(p)):.4e}" for p in (0, 30, 45, 60, 90)])         # f at 0°, 30°, 45°, 60°, 90° [1/s]
w_rel, Om_v = np.array([1e-5, 2e-5, 3e-5]), np.array([0, 0, rotating.OMEGA_EARTH])               # a relative vorticity, the Earth's Ω
print(np.allclose(kinematics.vorticity_in_rotating_frame(ch05.absolute_vorticity(w_rel, Om_v), Om_v), w_rel))   # round trip: True
""")
core("C11", "Absolute circulation and the fluid column", r"""
Air crosses a mountain range and is squashed into a thinner layer; a ring of air drifts toward the pole. On a spinning
planet, what is conserved — and which way must each turn?
""", eqs=("5.33",))
problem(r"""
Kelvin's theorem failed in a rotating frame because the Coriolis force is not conservative (C03). Something must replace
it — otherwise the atmosphere would have no rule for how cyclones spin up over plains and anticyclones sit over
mountains. The replacement is the circulation of the *absolute* vorticity, and for a column it becomes (ζ + f)/h = const:
potential vorticity, the single most useful conserved quantity of Ch. 13.
""")
idea("""
inertial frame:  Gamma conserved (Kelvin)                                                                   (5.8)
rotating frame:  Gamma changes, because the Coriolis loop integral = -2 Omega . (rate of the loop's vector area)
                 ->  Gamma + 2 Omega . A_vec  =  Int (w + 2 Omega) . n dA  =  Gamma_a  is conserved            (5.33)
thin column:     (zeta + f) A = const  and  A h = const   ->   (zeta + f)/h = const      stretch -> cyclonic, squash -> anticyclonic
""")
nb.md(r"""
*Reminder — loop vector area* (primed in C05, P138): $\mathbf A_{vec}=\tfrac12\oint\mathbf x\times d\mathbf x$, the loop's area
times its normal; for a horizontal loop on the Earth only its vertical part meets the vertical planetary vorticity.
""")
remind([
    ("product rule for a cross product", r"$d(\mathbf a\times\mathbf b)=d\mathbf a\times\mathbf b+\mathbf a\times d\mathbf b$, keeping the order (Ch. 4 P125)."),
], lead="needed in the derivation below")
D("D19", ref="5.33")
nb.md(r"""
*Gloss — column mass.* A vertical column of incompressible fluid keeps its volume: section × height = const, so a
taller column is a thinner one.
""")
D("D20", ref="5.33")
note("N39", "Planetary stretching and tilting", r"""
With Ω = Ω e_z, $2(\boldsymbol\Omega\cdot\nabla)\mathbf u=2\Omega\,\partial\mathbf u/\partial z$, so, keeping only this term of (5.30),
the equations below follow. Stretching vertical *fluid* lines (not vortex lines — they contain 2Ω already) raises ω_z;
tilting them makes horizontal vorticity. By continuity $\partial w/\partial z=-(\partial u/\partial x+\partial v/\partial y)$:
horizontal convergence spins fluid up. This is D20 step 6's small-change form. The cell: a convergence that stretches
columns at ∂w/∂z = 2×10⁻⁵ s⁻¹ gives 2Ω × 2×10⁻⁵ = 2.9168×10⁻⁹ s⁻².
""", equation=r"\frac{D\omega_z}{Dt}=2\Omega\frac{\partial w}{\partial z},\qquad\frac{D\omega_x}{Dt}=2\Omega\frac{\partial u}{\partial z},\qquad\frac{D\omega_y}{Dt}=2\Omega\frac{\partial v}{\partial z}")
nb.code(r"""
G = np.array([[-0.1, 0, 0], [0, -0.1, 0], [0, 0, 0.2]])*1e-4         # horizontal convergence feeding vertical stretching [1/s]
print([f"{v:.4e}" for v in ch05.planetary_vorticity_terms(G, np.array([0, 0, rotating.OMEGA_EARTH]))])   # 2(Ω·∇)u: (0, 0, 2Ω × 2e-5) = 2.9168e-9 1/s²
""")
nb.worked_example("a column over a ridge and a ring moving poleward", r"""
1. **Column** at 45° N (f = 1.03×10⁻⁴ s⁻¹), at rest (ζ₀ = 0), 1000 m deep. It moves over a 200 m high ridge: h = 800 m.
   $\zeta=(0+f)(800/1000)-f=-0.2f$ = −2.06×10⁻⁵ s⁻¹: anticyclonic (clockwise in the north). Past the ridge, h = 1000 m
   again and ζ = 0. Into a trough 1200 m deep: ζ = +0.2f, cyclonic.
2. **Ring of air** of radius 500 km (A = 7.85×10¹¹ m²) at rest at 30° N, carried to 60° N keeping its area:
   $\Gamma_1=2\Omega A(\sin30°-\sin60°)$ = 1.458×10⁻⁴ × 7.854×10¹¹ × (−0.366) = −4.19×10⁷ m²/s; mean ζ = Γ₁/A =
   −5.34×10⁻⁵ s⁻¹ — anticyclonic, about half the local f.
""")
nb.code(r"""
f45 = rotating.coriolis_parameter(np.deg2rad(45.0))                  # f at 45° N [1/s]
print([f"{ch05.column_relative_vorticity(h, 1000.0, 0.0, f45):+.4e}" for h in (800.0, 1000.0, 1200.0)])   # ζ over a ridge, flat, trough [1/s]
print(f"{ch05.column_relative_vorticity(1100.0, 1000.0, 0.0, 1e-4):.2e} 1/s")   # D20's number: +10 % height at f = 1e-4 → +1.0e-5
A = np.pi*5e5**2                                                     # a ring of radius 500 km [m²]
G1 = ch05.relative_circulation_after_move(0.0, A, 30.0, A, 60.0)     # at rest at 30° N, moved to 60° N [m²/s]
print(f"Γ₁ = {G1:.5e} m²/s, mean ζ = {G1/A:.5e} 1/s")                 # −4.19261e7, −5.33820e-5: anticyclonic
rsc = ch05.kelvin_scenario("rotating"); t = np.linspace(0, 5, 6)      # a converging flow seen in a frame turning at 0.5 rad/s
Pr = ch05.material_loop(rsc["u"], rsc["pts0"], t)                    # a material loop (unit circle at t = 0) carried by it
print([tuple(round(g, 6) for g in ch05.absolute_circulation(rsc["u"], Pk, np.array([0, 0, 0.5]), t=tk)) for Pk, tk in zip(Pr, t)])   # (Γ, Γ_a)
u_still = lambda x, t: np.stack([0.5*np.asarray(x)[1], -0.5*np.asarray(x)[0]])   # fluid at rest in the inertial frame, seen rotating: −Ω×x
print(ch05.absolute_circulation(u_still, rsc["pts0"], np.array([0, 0, 0.5])))   # (Γ, Γ_a) = (−2ΩA, 0) = (−π, 0)
""", explain=r"""
1. The column rule of D20: over the ridge ζ = −2.0625×10⁻⁵ s⁻¹ (anticyclonic), flat 0, over the trough +2.0625×10⁻⁵ s⁻¹
   (cyclonic); and D20's own number, +1.0×10⁻⁵ s⁻¹ for a 10 % taller column at f = 10⁻⁴ s⁻¹.
2. The poleward ring of $\frac{D\Gamma_a}{Dt}=0$ (5.33): Γ₁ = −4.19261×10⁷ m²/s, mean ζ = −5.33820×10⁻⁵ s⁻¹.
3. A material loop in a converging rotating flow: Γ grows 0 → 1.985865 m²/s in 5 s while Γ_a stays at π = 3.141593 m²/s.
4. Fluid at rest in the inertial frame, seen from the rotating frame, circulates backwards: Γ = −2ΩA = −π, Γ_a = 0.
""")
nb.md("**From scratch — the vector area by chords, then Γ_a by hand:**")
nb.check_agree(r"""
Pk = Pr[-1]                                                          # the loop at t = 5 s, shape (2, 256) [m]
P3 = np.vstack([Pk, 0*Pk[0]])                                        # as 3-D points in the plane z = 0
A_hand = 0.5*np.sum(np.cross(P3.T, (np.roll(P3, -1, axis=1) - P3).T), axis=0)   # ½ Σ x × (chord) [m²]
A_lib = ch05.loop_vector_area(Pk)                                    # the library (spectral) [m²]
print(A_hand, A_lib)                                                 # (0, 0, A): chords cut the corners by ~1e-4
assert np.allclose(A_hand, A_lib, rtol=5e-4)                         # same vector area
Gam = ch05.loop_circulation(rsc["u"], Pk, 5.0)                       # relative circulation Γ [m²/s]
print(Gam + 2*0.5*A_hand[2], Gam + 2*0.5*A_lib[2], np.pi)            # Γ_a = Γ + 2Ω·A_vec, by chords and spectrally, vs π
assert np.allclose(Gam + 2*0.5*A_hand[2], np.pi, rtol=1e-3) and np.allclose(Gam + 2*0.5*A_lib[2], np.pi, rtol=1e-8)
""")
nb.figure(r"""
xs_ = np.linspace(-3e5, 3e5, 241)                                    # distance along the path [m]
hf = lambda x: 1000 - 200*np.exp(-((x + 1.5e5)/5e4)**2) + 200*np.exp(-((x - 1.5e5)/5e4)**2)   # layer depth: ridge, then trough [m]
cs = ch05.column_over_slope(xs_, hf, 45.0)                           # ζ along the path at 45° N, starting at rest [1/s]
fig, (a, b) = plt.subplots(2, 1, figsize=(8, 5.6), sharex=True)      # the layer | the numbers
column_sketch(a, xs_/1e3, cs["h"], cs["zeta"], n_cols=7)             # columns with spin marks (teal ↺ cyclonic, rose ↻ anticyclonic)
a.set_xlabel(""); a.set_ylabel("height [m]"); a.set_title("A rotating layer over a ridge and a trough (45° N)", fontsize=10)
b.plot(xs_/1e3, cs["zeta"]/cs["f"], color=COLORS["teal"], lw=4, label="ζ/f")             # relative vorticity in units of f (wide)
b.plot(xs_/1e3, cs["h"]/1000.0 - 1, color=COLORS["ink"], ls="--", lw=1.5, label="h/h₀ − 1")  # depth change, dashed on top: identical
b.plot(xs_/1e3, cs["ratio"]/cs["ratio"][0] - 1, color=COLORS["amber"], lw=3, label="(ζ + f)/h, relative change")   # conserved
b.set_xlabel("x [km]"); b.set_ylabel("[–]"); b.legend(fontsize=8)
savefig(fig, "ch05", "column_ridge"); plt.show()                     # save the PNG to outputs/ch05, then draw
""", see="(Note `N54`, Fig. 5.10 analogue.) Top, a layer that thins over a ridge and deepens over a trough, with columns turning clockwise (rose) over the "
         "ridge and counterclockwise (teal) over the trough; bottom, the wide teal ζ/f with the dashed black h/h₀ − 1 lying exactly on it, and the "
         "amber (ζ + f)/h flat at zero change.",
    read="(ζ + f)/h is carried unchanged by each column (D20): squashed columns lose absolute vorticity (anticyclonic), "
         "stretched ones gain it (cyclonic); ζ/f equals h/h₀ − 1 for a column that starts at rest.",
    change="…the same layer at 10° N: f is 4.1 times smaller (sin 10° vs sin 45°), so every ζ is 4.1 times smaller — the "
           "equator has little planetary vorticity to lend.")
nb.plotly(r"""
hh = np.linspace(800, 1200, 41)                                      # column height [m], h₀ = 1000 m
z45 = ch05.column_relative_vorticity(hh, 1000.0, 0.0, rotating.coriolis_parameter(np.deg2rad(45.0)))   # the 45° N reference


def lat(phi):                                                        # latitude [deg] → curves
    f = rotating.coriolis_parameter(np.deg2rad(phi))                 # f = 2Ω sin φ [1/s]
    return {"ζ(h) at this latitude": (hh, ch05.column_relative_vorticity(hh, 1000.0, 0.0, f)*1e5),   # [10⁻⁵ 1/s]
            "45° N reference": (hh, z45*1e5)}


fig = slider_figure(lat, "φ", np.linspace(-90, 90, 19 if not FAST else 10), unit="°", xlabel="column height h [m]",  # precompute every slider position
                    ylabel="relative vorticity ζ [10⁻⁵ 1/s]", title="The planet lends more spin toward the pole")
recolor(fig, {"ζ(h) at this latitude": COLORS["teal"], "45° N reference": COLORS["muted"]}, dashes={"45° N reference": "dash"})  # house colours by trace name
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="A teal line through (1000 m, 0) whose slope grows with |φ|, flat at the equator and reversed in the Southern "
        "Hemisphere; the grey dashed line is 45° N.",
    read="For a column starting at rest, ζ = f(h/h₀ − 1): the slope is f/h₀ = 2Ω sin φ/h₀. Stretching (h > h₀) gives "
         "ζ with the sign of f — cyclonic in either hemisphere.",
    change="…the Southern Hemisphere (φ < 0): f < 0, the line tilts the other way — cyclones turn clockwise there.")
explainer("vorticity_equation_rotating", "What does a spinning planet add to the vorticity equation?",
          "dragging a column over a ridge, moving a ring of air in latitude and switching the scene of a term-bar budget "
          "of (5.30) shows the same conservation law — and the same equation — in three regimes that a single figure "
          "cannot hold at once.",
          "The Derivation tab carries the ★★★ D15 and D14, D19, D20.", [
              "Column mode: drag h below h₀ and watch the column turn anticyclonic while (ζ + f)/h stays put.",
              "Ring mode: move the ring from 30° N to 60° N and read Γ = −4.19×10⁷ m²/s.",
              "Budget mode: switch scenes and watch which bar dominates.",
              "Southern Hemisphere preset: every sign flips.",
          ])
nb.md(r"""
> ⚠️ **Common confusion:** "2ΩA is the planetary part." It is $2\boldsymbol\Omega\cdot\mathbf A_{vec}$ — only the component of
> Ω along the loop's normal counts. On the Earth, for a horizontal loop, that is fA = 2Ω sin φ A.
""")
whatif(r"""
…the layer were also stratified? The column would be a slab between two density surfaces and h its thickness: (ζ + f)/h
becomes Ertel's potential vorticity (Ch. 13) — the same bookkeeping.
""")

# =====================================================================================================================
# A.7 §5.7 Interaction of Vortices — C12, C13
# =====================================================================================================================
nb.section("5.7", "Interaction of Vortices", intro=r"""
**What is this section about?** Idealise each vortex as a line vortex. A straight line vortex cannot move itself, but
Helmholtz says it moves with the fluid — and the fluid at its position is moved by all the *other* vortices
(Biot–Savart). That one rule makes like vortices orbit, opposite ones travel together, a vortex slide along a wall (its
mirror image does the pushing) and smoke rings leap-frog.
""")

# ---- C12 ----------------------------------------------------------------------------------------------------------
core("C12", "Point vortices move each other: pairs and the centre of vorticity", r"""
Two hurricanes that come within a thousand kilometres start to circle each other; two dimples made by a paddle stroke
glide away side by side. A vortex cannot push itself — so how do vortices move?

*In one line:* $\frac{d\mathbf x_k}{dt}=\sum_{j\ne k}\frac{\Gamma_j}{2\pi}\frac{\mathbf e_z\times(\mathbf x_k-\mathbf x_j)}{\lvert\mathbf x_k-\mathbf x_j\rvert^2}$
""")
problem(r"""
Weather forecasters watch binary cyclones dance (the Fujiwhara effect); airport controllers space landing aircraft
because each plane leaves a pair of vortices that sinks behind it; modellers of 2-D turbulence and of ocean eddies use
swarms of point vortices. All of these rest on one simple rule.
""")
idea("""
real vortex (core with vorticity)  ~  line vortex of the same Gamma      (outside the core the flow is irrotational)
each vortex k moves with the velocity induced at x_k by all the OTHERS:  dx_k/dt = sum_{j != k} (Gamma_j/2pi) e_z x (x_k - x_j)/|x_k - x_j|^2
two alike     ->  orbit their centre of vorticity         two opposite  ->  translate together at Gamma/(2 pi h)
""")
nb.md(r"""
*Glosses used below.* **Superposition.** Outside the cores the flow is irrotational and (Ch. 6) obeys the linear
Laplace equation, so the velocities of several vortices add. A **point vortex** is a line vortex seen in the plane:
$u_\theta=\Gamma/2\pi r$ (5.2) around one point. **Lever rule.** The balance point of two weights w₁, w₂ a distance h
apart is $h\,w_2/(w_1+w_2)$ from the first. **Circular motion.** A point turning at rate $\dot\theta$ at distance d from
the centre moves at $d\,\dot\theta$.
""")
P("P147", "systems of ODEs for interacting bodies and their invariants", r"""
N vortices give 2N coupled ODEs (each position's rate depends on all the others). `solve_ivp` integrates them as one
state vector. Because the exact motion conserves some quantities — here the linear impulse ΣΓx, the angular impulse
ΣΓ|x|² and the energy $H=-\frac1{2\pi}\sum_{j<k}\Gamma_j\Gamma_k\ln\lvert\mathbf x_j-\mathbf x_k\rvert$ (Kirchhoff's function) —
watching them stay constant tells us the numerical solution is trustworthy.
""", code=r"""
xv0 = np.array([[-0.5, 0.5], [0.0, 0.0]]); G = np.array([1.0, 1.0])   # two like vortices 1 m apart [m], Γ = 1 m²/s each
t = np.linspace(0, 20, 5)                                            # five output times [s]
X = ch05.point_vortex_evolve(xv0, G, t)                              # positions, shape (T, 2, M) [m]
print([round(ch05.point_vortex_invariants(Xk, G)['I'], 10) for Xk in X])   # ΣΓ|x|² = 0.5 m⁴/s every time
""")
gloss(["Kirchhoff energy H with a logarithm"])
D("D21", ref="5.2")
note("N40", "Two like vortices (Fig. 5.11)", r"""
Γ₁, Γ₂ > 0 a distance h apart push each other sideways at the speeds below, in opposite directions, so the pair turns
about its centre of vorticity G, at $h_1=\Gamma_2h/(\Gamma_1+\Gamma_2)$ from vortex 1, at the rate $(\Gamma_1+\Gamma_2)/2\pi h^2$
(D21; the book leaves G to Exercise 5.18). **Numbers:** Γ₁ = Γ₂ = 1 m²/s, h = 1 m → each moves at 0.159 m/s, period
19.7 s; Γ₂ = 3Γ₁ → G is ¾h from vortex 1. ⚠️ The caption of Fig. 5.11 calls G the point where the induced velocity is
zero — true only for equal strengths.
""", equation=r"V_1=\Gamma_1/2\pi h,\qquad V_2=\Gamma_2/2\pi h")
note("N41", "Equal and opposite pair (Figs. 5.12–5.13)", r"""
each moves the other at Γ/2πh in the same direction, so the pair translates at the speed below relative to the fluid —
made by a paddle stroke or a knife blade drawn briefly through a bucket (Lighthill's picture, our simulation in C13).
The wake of an aircraft is such a pair, sinking behind it (Ch. 14).
""", equation=r"V=\Gamma/(2\pi h)")
nb.worked_example("two like vortices, then two opposite", r"""
1. Γ₁ = Γ₂ = 1 m²/s at (−0.5, 0) and (0.5, 0) m: vortex 2 is pushed up at $1/(2\pi\times1)$ = 0.159 m/s, vortex 1 down
   at 0.159 m/s; G is the midpoint; the joining line turns at (0.159 + 0.159)/1 = 0.318 rad/s — one orbit in
   2π/0.318 = 19.7 s.
2. Γ₂ = 3 m²/s instead: V₁ = 0.159, V₂ = 0.477 m/s; G at 0.75 m from vortex 1 (the lever rule with weights 1 and 3);
   rate $(1+3)/2\pi$ = 0.637 rad/s.
3. Γ₂ = −1 m²/s: both move up at 0.159 m/s — a straight march.
""")
nb.code(r"""
for G2 in (1.0, 3.0, -1.0):                                           # the three cases of the example
    vp = ch05.vortex_pair(1.0, G2, 1.0)                               # D21's closed forms, h = 1 m
    print(G2, {k: round(v, 6) for k, v in vp.items()})
Pq = ch05.point_vortex_preset("equal_pair"); t = np.linspace(0, ch05.vortex_pair(1.0, 1.0, 1.0)["period"], 200)   # one predicted period, 2π/rate
Xe = ch05.point_vortex_evolve(Pq["xv"], Pq["Gamma"], t)               # the full integration [m]
print(f"back at the start after one period: {np.abs(Xe[-1] - Xe[0]).max():.1e} m")
inv = [ch05.point_vortex_invariants(Xk, Pq["Gamma"]) for Xk in Xe]    # invariants along the run
print(f"spread of P = ΣΓx: {np.ptp([v['P'] for v in inv], axis=0).max():.1e}, of I = ΣΓ|x|²: {np.ptp([v['I'] for v in inv]):.1e},  of H: {np.ptp([v['H'] for v in inv]):.1e}")   # all flat
xv_u, G_u = ch05.point_vortex_preset("unequal_pair")["xv"], np.array([1.0, 3.0])   # Γ₂ = 3Γ₁
Gc = ch05.centre_of_vorticity(xv_u, G_u); print("centre of vorticity:", Gc)   # (0.25, 0): ¾ of the way from −0.5 to 0.5
print("fluid velocity AT G:", ch05.point_vortex_velocity(Gc, xv_u, G_u))      # (0, −1.70) m/s: G is not a point of rest
""", explain=r"""
1. The closed forms of D21: equal pair V = 0.159155 m/s, centre 0.5 m, rate 0.318310 rad/s, period 19.739 s; unequal
   pair centre 0.75 m, rate 0.636620 rad/s, period 9.870 s; opposite pair: centre at infinity, no rotation, translation
   0.159155 m/s.
2. The full point-vortex integration returns to its start after the predicted period.
3. The invariants stay flat — the integrator is trustworthy.
4. The centre of vorticity of the unequal pair, and the fluid velocity there: −1.70 m/s — the pair turns about G, yet
   the fluid at G is not at rest (the Fig. 5.11 caption's slip).
""")
nb.md("**From scratch — $d\\mathbf x_k/dt$ by a double loop over pairs** ($\\mathbf e_z\\times(a,b)=(-b,a)$: each pair of vortices, one line):")
nb.check_agree(r"""
P3v = ch05.point_vortex_preset("three_vortices"); xv, Gm = P3v["xv"], P3v["Gamma"]   # three vortices, Γ = (1, 1, −0.5)
M = xv.shape[1]                                                      # number of vortices
v = np.zeros((2, M))                                                 # their velocities [m/s]
for k in range(M):                                                   # the vortex being moved
    for j in range(M):                                               # every other vortex pushes it
        if j != k:
            d = xv[:, k] - xv[:, j]                                  # x_k − x_j [m]
            v[:, k] += Gm[j]/(2*np.pi)*np.array([-d[1], d[0]])/(d @ d)   # (Γ_j/2π) e_z × d/|d|²
print(v)
assert np.allclose(v, ch05.point_vortex_rhs(xv, Gm), atol=1e-14)     # = the library
""")
nb.md(r"""
**Pairs in motion (Figs. 5.11–5.12 analogues, notes `N55`, `N56`):** left, the equal pair orbiting G (black cross); middle, Γ₂ = 3Γ₁
orbiting its off-centre G; right, the opposite pair (teal +, rose −) marching up.
""")
nb.animation(r"""
nf = 30 if not FAST else 16                                          # frames
tf = np.linspace(0, 19.7392, nf)                                     # one period of the equal pair [s]
runs = []                                                            # (positions, strengths, centre) for each panel
for name in ("equal_pair", "unequal_pair", "opposite_pair"):                                        # the three panels
    pr = ch05.point_vortex_preset(name)                              # positions and strengths
    runs.append((ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], tf), pr["Gamma"], ch05.centre_of_vorticity(pr["xv"], pr["Gamma"])))  # positions over one period, strengths, centre of vorticity
with plt.rc_context(ANIM_RC):                                         # light settings for fast frames
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.8), dpi=80)                                      # three panels
fig.subplots_adjust(left=0.05, right=0.98, bottom=0.1, top=0.85, wspace=0.25)   # fixed margins
dots, trails = [], []                                                # artists updated per frame
for a, (X, Gm, Gc), ttl in zip(axs, runs, ["Γ₂ = Γ₁: orbit G", "Γ₂ = 3Γ₁: orbit G", "Γ₂ = −Γ₁: march"]):  # set up each panel
    a.set_xlim(-1.0, 1.0); a.set_ylim(-1.0 if ttl != "Γ₂ = −Γ₁: march" else -0.3, 1.0 if ttl != "Γ₂ = −Γ₁: march" else 3.4)  # fixed limits (the marching pair needs more room upward) [m]
    a.set_aspect("equal"); a.locator_params(nbins=4); a.set_title(ttl, fontsize=10)
    if a is not axs[0]:                                              # ticks on the first panel only (fewer labels to draw per frame)
        a.set_xticks([]); a.set_yticks([])                 # equal axes, few ticks, a title
    if np.all(np.isfinite(Gc)):                                                                     # the opposite pair has no finite centre
        a.plot(*Gc, "x", color=COLORS["ink"], ms=10, mew=2)          # the centre of vorticity (fixed)
    cols = [COLORS["teal"] if g > 0 else COLORS["rose"] for g in Gm]                                # teal counterclockwise, rose clockwise
    trails.append([a.plot([], [], color=c, lw=1, alpha=0.6)[0] for c in cols])                      # empty trail lines, filled in update()
    dots.append(a.scatter(X[0, 0], X[0, 1], c=cols, s=60, zorder=3))                                # the vortices at t = 0


def update(i):                                                       # draw frame i
    for (X, Gm, Gc), dt_, tr in zip(runs, dots, trails):                                            # each panel
        dt_.set_offsets(X[i].T)                                      # vortex positions now
        for m, line in enumerate(tr):                                                               # each vortex's trail
            line.set_data(X[:i + 1, 0, m], X[:i + 1, 1, m])          # their trails
    fig.suptitle(f"t = {tf[i]:.1f} s", fontsize=10)                                                 # the clock
    return dots                                                                                     # the artists that changed


with plt.rc_context(ANIM_RC):                                         # render the frames with the light settings
    show_animation(animate(update, frames=nf, fig=fig, interval=150), player="video")   # smooth movie
""")
see_read_change(
    see="Two circles of equal size round a fixed cross; two circles of different sizes round a cross nearer the stronger "
        "vortex; a teal–rose pair moving straight up.",
    read="Every vortex is carried by the others: like ones orbit their centre of vorticity (which never moves), opposite "
         "ones travel together at Γ/2πh. The equal pair completes exactly one orbit in the 19.74 s shown; the unequal "
         "pair two.",
    change="…h halved: the orbits are four times faster ((Γ₁ + Γ₂)/2πh²) and the march twice as fast (Γ/2πh).")
nb.figure(r"""
tl = np.linspace(0, 20, 401 if not FAST else 201)                     # 20 s [s]
names = ["equal_pair", "unequal_pair", "opposite_pair", "three_vortices"]   # four presets
fig = plt.figure(figsize=(11, 5.6)); gs = fig.add_gridspec(2, 4, height_ratios=[1.3, 1])   # 4 trajectory panels + 1 invariant panel
axi = fig.add_subplot(gs[1, :])                                       # the invariants
for k, nm in enumerate(names):
    pr = ch05.point_vortex_preset(nm); X = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], tl)   # run 20 s
    a = fig.add_subplot(gs[0, k])
    vortex_marks(a, X[-1], pr["Gamma"], size=0.06, trails=X)          # trails and the final positions (teal +, rose −)
    a.set_title(nm.replace("_", " "), fontsize=9); a.locator_params(nbins=4)
    inv = [ch05.point_vortex_invariants(Xk, pr["Gamma"]) for Xk in X]
    I = np.array([v["I"] for v in inv]); H = np.array([v["H"] for v in inv])
    axi.plot(tl, I - I[0], label=f"{nm}: ΔI", lw=1.5)                 # angular impulse change
    axi.plot(tl, H - H[0], ls="--", lw=1.5, label=f"{nm}: ΔH")        # energy change
axi.set_xlabel("t [s]"); axi.set_ylabel("change since t = 0"); axi.set_ylim(-1e-8, 1e-8); axi.legend(fontsize=6, ncol=4)
fig.suptitle("Point-vortex motion and its invariants", fontsize=11)
plt.show()
""", see="Top: the equal pair's single circle, the unequal pair's two circles of different radius, the opposite pair's "
         "straight tracks and the three vortices' tangled paths. Bottom: every invariant's change stays inside ±10⁻⁸ "
         "(most lines lie on zero).",
    read="Like vortices orbit, opposite ones march, three can wander in complicated ways — but ΣΓ|x|² (I) and the energy "
         "H never change, which is how we know the integration is right. (The opposite pair has ΣΓ = 0, so its I is not "
         "tied to a fixed centre, yet it is still conserved.)",
    change="…a fourth vortex added: the motion can become chaotic, but the invariants stay flat.")
whatif(r"""
…one of the two vortices were replaced by a wall? The wall acts like a mirror vortex of opposite sign — the pair becomes
a vortex and its image (C13).
""")

# ---- C13 ----------------------------------------------------------------------------------------------------------
core("C13", "Images: a vortex drifts along a wall", r"""
A vortex near the bottom of a tank slides sideways along it; a smoke ring blown at a wall grows wider as it slows. How
does a wall push a vortex without touching it?

*In one line:* the wall is replaced by a mirror vortex of opposite sign at distance 2h, so the vortex drifts along the
wall at $V_A=\Gamma/4\pi h$.
""")
problem(r"""
A wall forbids flow through it. For vortices this boundary condition has a beautiful solution: pretend the wall is a
mirror and put an opposite vortex behind it. The mirror twin then moves the real vortex — and the same trick handles
circular tanks, channels, the ground under a landing aircraft (Ch. 14) and cylinders in a stream (Ch. 6).
""")
idea("""
     vortex A (+Gamma) o                       wall = the line where A's push and its twin's push have
-------------------------- wall                 equal and opposite normal parts  ->  no flow through it
     image B (-Gamma)  o                       A is moved only by B (distance 2h): V = Gamma/(2 pi 2h) = Gamma/(4 pi h), along the wall
""")
nb.md(r"""
*Gloss — mirror symmetry.* Reflect a counterclockwise vortex in the wall and reverse its spin: at every wall point the
two velocity vectors are mirror images, so their normal parts cancel and their tangential parts add (the method of
images).
""")
D("D22", ref="5.2")
nb.worked_example("a vortex half a metre above the floor", r"""
Γ = 1 m²/s at height h = 0.5 m.

1. Image: −1 m²/s at −0.5 m.
2. At the wall point directly below, both vortices are 0.5 m away; each gives $1/(2\pi\times0.5)$ = 0.318 m/s, both
   horizontal and in the same direction: 0.637 m/s along the floor, 0 through it.
3. At the wall point 0.5 m to the side, distance √0.5 = 0.707 m: each gives 0.225 m/s at 45°; the vertical parts
   cancel, the horizontal ones add to 0.318 m/s.
4. The vortex itself moves at $1/(4\pi\times0.5)$ = 0.159 m/s along the wall.
5. In a channel 1 m high with the vortex at 0.25 m, all the images add up to $(\Gamma/4H)\cot(\pi h/H)$ = 0.25 m/s —
   slower than the single-wall 0.318 m/s, because the ceiling's images push back.
""")
nb.code(r"""
xa, Ga = ch05.wall_image_system(np.array([[0.0], [0.5]]), np.array([1.0]))   # the vortex and its image
xw = np.array([np.linspace(-3, 3, 201), np.zeros(201)])               # 201 points on the wall y = 0 (x = 0 among them) [m]
vel = ch05.point_vortex_velocity(xw, xa, Ga)                          # velocity there [m/s]
print(f"largest |v| through the wall: {np.abs(vel[1]).max():.1e} m/s; u below the vortex: {vel[0][100]:.4f} m/s")   # 0; 0.6366
print(ch05.vortex_near_wall_speed(1.0, 0.5), ch05.point_vortex_rhs(np.array([[0.0], [0.5]]), np.array([1.0]), boundary="wall").ravel())   # drift 0.159155
print(ch05.channel_image_velocity(0.25, 1.0, 1.0), 1/(4*1.0)/np.tan(np.pi*0.25))   # channel: image series vs (Γ/4H)cot(πh/H)
print(f"mid-channel, h = H/2: {ch05.channel_image_velocity(0.5, 1.0, 1.0):.1e} m/s")   # 0: the two walls' images cancel
""", explain=r"""
1. The vortex and its image (opposite strength, mirrored position).
2. The velocity at 201 wall points has no normal part (0 to round-off) — the wall is a streamline; directly below the
   vortex it slides at 0.6366 m/s.
3. The vortex's own drift, from the image only: Γ/4πh = 0.159155 m/s along +x.
4. The channel's infinite image series (summed until the next term is below 10⁻¹²) equals its closed form, 0.25 m/s.
""")
nb.md("**From scratch — the image sum by hand**, at the vortex and at one wall point:")
nb.check_agree(r"""
d = np.array([0.0, 0.5]) - np.array([0.0, -0.5])                     # vortex A minus image B [m]
vA = -1.0/(2*np.pi)*np.array([-d[1], d[0]])/(d @ d)                  # B (Γ = −1) pushes A: (Γ_B/2π) e_z × d/|d|²
print(vA)                                                            # (0.159155, 0): along the wall
assert np.allclose(vA, [1/(4*np.pi*0.5), 0.0])                       # = Γ/4πh
P_w = np.array([1.0, 0.0])                                           # a wall point [m]
dA_, dB_ = P_w - np.array([0.0, 0.5]), P_w - np.array([0.0, -0.5])   # from A and from B to the wall point
vAw = 1.0/(2*np.pi)*np.array([-dA_[1], dA_[0]])/(dA_ @ dA_)          # A's push there
vBw = -1.0/(2*np.pi)*np.array([-dB_[1], dB_[0]])/(dB_ @ dB_)         # B's push there
assert abs(vAw[1] + vBw[1]) < 1e-15                                  # normal parts cancel
""")
nb.figure(r"""
xs9 = np.linspace(-2, 2, 9)                                          # 9 wall points [m]
Pw = np.array([xs9, 0*xs9])                                          # on y = 0
vA_ = ch05.point_vortex_velocity(Pw, np.array([[0.0], [0.5]]), np.array([1.0]))    # A's pushes [m/s]
vB_ = ch05.point_vortex_velocity(Pw, np.array([[0.0], [-0.5]]), np.array([-1.0]))  # the image's pushes [m/s]
fig, ax = plt.subplots(figsize=(8, 3.8))
ax.axhline(0, color=COLORS["ink"], lw=2)                             # the wall
ax.axhspan(-1.0, 0, color=COLORS["grid"], alpha=0.6)                 # behind the wall (the mirror world)
vortex_marks(ax, np.array([[0.0], [0.5]]), [1.0], size=0.08)         # A, teal (counterclockwise)
vortex_marks(ax, np.array([[0.0], [-0.5]]), [-1.0], size=0.08)       # B, rose (clockwise), behind the wall
s_ = 0.6                                                             # arrow scale [m per m/s]
ax.quiver(xs9, 0*xs9, vA_[0], vA_[1], color=COLORS["teal"], angles="xy", scale_units="xy", scale=1/s_, width=0.004)   # from A
ax.quiver(xs9, 0*xs9, vB_[0], vB_[1], color=COLORS["rose"], angles="xy", scale_units="xy", scale=1/s_, width=0.004)   # from B
ax.quiver(xs9, 0*xs9, vA_[0] + vB_[0], vA_[1] + vB_[1], color=COLORS["ink"], angles="xy", scale_units="xy", scale=1/s_, width=0.007)   # sum
ax.set_xlim(-2.4, 2.4); ax.set_ylim(-1.0, 1.0); ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
ax.set_title("The image makes the wall a streamline", fontsize=11)
savefig(fig, "ch05", "wall_image"); plt.show()                       # save the PNG to outputs/ch05, then draw
""", see="(Note `N58`, Fig. 5.14 analogue.) The teal vortex A above the wall, its rose image B below (in the grey mirror world); at nine wall points a thin "
         "teal arrow (A's push), a thin rose arrow (B's push) and a thick black arrow (their sum) that lies flat along "
         "the wall, longest under the vortex.",
    read="The image makes the wall a streamline; the black arrows are the slip velocity an inviscid fluid is allowed "
         "along a wall (a real, viscous fluid brings it to rest in a boundary layer, Ch. 9).",
    change="…an image of the same sign: the vertical parts would add instead of cancel — the black arrows would cross "
           "the wall, so it would no longer be a wall.")
note("N57", "The knife-blade pair in a bucket (Fig. 5.13, after Lighthill)", r"""
A knife drawn briefly through a bucket leaves an opposite pair that marches off; nearing the wall the two vortices
separate and run along it in opposite directions, each pushed by its own image. For a circle of radius a the image of a
vortex at x is at the inverse point $a^2\mathbf x/\lvert\mathbf x\rvert^2$ (measured from the centre, strength −Γ) —
below, our simulation.
""")
nb.animation(r"""
PK = ch05.point_vortex_preset("knife_bucket")                        # Γ = (−1, +1) at y = ±0.15 m; circle a = 1 m about (0.3, 0)
nf = 30 if not FAST else 16                                          # frames
tk = np.linspace(0, PK["t_end"], nf)                                 # 0 … 12 s
XK = ch05.point_vortex_evolve(PK["xv"], PK["Gamma"], tk, boundary=PK["boundary"], **PK["bp"])   # with circle images [m]
cx, cy = PK["bp"]["center"]; ac = PK["bp"]["a"]                      # the bucket's centre and radius [m]
with plt.rc_context(ANIM_RC):                                         # light settings for fast frames
    fig, ax = plt.subplots(figsize=(5.6, 5.0), dpi=80)
fig.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.9)     # fixed margins
th = np.linspace(0, 2*np.pi, 200); ax.plot(cx + ac*np.cos(th), cy + ac*np.sin(th), color=COLORS["ink"], lw=2)   # the bucket wall
cols = [COLORS["teal"] if g > 0 else COLORS["rose"] for g in PK["Gamma"]]
trl = [ax.plot([], [], color=c, lw=1.2)[0] for c in cols]             # trails
vd = ax.scatter(XK[0, 0], XK[0, 1], c=cols, s=60, zorder=3)           # the two vortices
im = ax.scatter([], [], c=[], s=30, alpha=0.35)                       # their images, outside the circle
ax.set_xlim(cx - 2.2, cx + 2.2); ax.set_ylim(-2.2, 2.2); ax.set_aspect("equal"); ax.locator_params(nbins=5)


def update(i):                                                       # draw frame i
    vd.set_offsets(XK[i].T)                                          # the vortices now
    for m, line in enumerate(trl):
        line.set_data(XK[:i + 1, 0, m], XK[:i + 1, 1, m])            # their trails
    xa_, Ga_ = ch05.circle_image_system(XK[i], PK["Gamma"], ac, center=(cx, cy), inside=True)   # vortices + images
    im.set_offsets(xa_[:, 2:].T); im.set_color([COLORS["teal"] if g > 0 else COLORS["rose"] for g in Ga_[2:]])   # images faint
    ax.set_title(f"knife-blade pair in a bucket, t = {tk[i]:.1f} s", fontsize=10)
    return vd, im


with plt.rc_context(ANIM_RC):                                         # render the frames with the light settings
    show_animation(animate(update, frames=nf, fig=fig, interval=150), player="video")   # smooth movie
""")
see_read_change(
    see="The pair crosses the bucket toward the left, splits at the wall and the two vortices run round it in opposite "
        "directions; their faint images move outside the circle, mirroring them.",
    read="Far from the wall the pair translates at Γ/2πh (C12); near it each vortex is dragged by its own image — the "
         "wall drift of D22.",
    change="…a stronger stroke (larger Γ): the same paths, travelled faster.")
note("N42", "A smoke ring meeting a wall (Fig. 5.15)", r"""
A ring moves itself (its own curvature, the thin-core self-speed below — Kelvin's formula, which we cite) and is moved
by its image ring behind the wall. The image pushes it outward and holds it back: the ring **widens** and **slows** as it
approaches, never reaching the wall. We compute the path with the velocity one circular filament induces at another
(elliptic integrals — primer below) and the self-speed. **Number:** Γ = 1 m²/s, R = 1 m, core a = 0.1 m: U = 0.329 m/s.
⚠️ This is a **qualitative model**: thin cores that keep their circular shape and their volume a²R; it stops being
trustworthy once the gap to the wall shrinks to a few core radii — `ring_dynamics` stops it when the gap falls below
three core radii (its `valid` mask and `stop_time`), because a core that close to its image is no longer thin and
circular.
""", equation=r"U=\frac{\Gamma}{4\pi R}\Big[\ln\frac{8R}{a}-\frac14\Big]")
P("P148", "complete elliptic integrals with the parameter m", r"""
The velocity induced by a circular filament involves $K(m)=\int_0^{\pi/2}d\theta/\sqrt{1-m\sin^2\theta}$ and
$E(m)=\int_0^{\pi/2}\sqrt{1-m\sin^2\theta}\,d\theta$. ⚠️ `scipy.special.ellipk(m)` and `ellipe(m)` take the **parameter
m = k²**, not the modulus k — passing k is a classic silent bug (K → ∞ as m → 1: a point near the filament).
""", code=r"""
from scipy.special import ellipk, ellipe                  # complete elliptic integrals K(m), E(m)
from scipy.integrate import quad                          # to check the definition
m = 0.5                                                   # the parameter m = k²
K = quad(lambda th: 1/np.sqrt(1 - m*np.sin(th)**2), 0, np.pi/2)[0]   # K(m) from its definition
print(ellipk(m), K)                                       # 1.854075 1.854075: same definition
print(ellipk(0.5**0.5))                                   # 2.085974 — what you get by passing k instead of m
""")
remind([
    ("np.gradient", "`np.gradient(z, t)` differentiates sampled data by central differences (one-sided at the ends) (Ch. 1 P22)."),
])
nb.code(r"""
print(ch05.ring_self_velocity(1.0, 0.1, 1.0))                        # Kelvin's self-speed: R = 1 m, a = 0.1 m, Γ = 1 m²/s → 0.328816 m/s
print(ch05.ring_ring_velocity(0.0, 0.5, 0.5, 0.0, 1.0), ch05.ring_axis_velocity(0.5, 0.5, 1.0))   # on the axis: elliptic formula = closed form
tw = np.linspace(0, 3.5, 141 if not FAST else 71)                    # 3.5 s of the approach [s]
res = ch05.ring_dynamics([dict(R=0.5, z=0.0, Gamma=1.0, a=0.05)], tw, wall_z=2.0)   # Γ > 0 moves toward +z: a wall at z = 2 m
ok = res["valid"]                                                    # samples where the thin-core model still applies
tv, Rr, zr = tw[ok], res["R"][ok, 0], res["z"][ok, 0]                 # time, radius and height of the ring [s], [m]
a_end = res["a"][ok, 0][-1]                                          # the core has thinned as the ring widened (volume a²R fixed) [m]
print(f"model stops at t = {res['stop_time']:.3f} s ({res['stop_reason']}: gap below 3 core radii, 3a = {3*a_end:.3f} m by then)")
for k in (0, len(tv)//3, 2*len(tv)//3, -1):
    print(f"t = {tv[k]:.2f} s: R = {Rr[k]:.4f} m, z = {zr[k]:.4f} m, gap to the wall {2.0 - zr[k]:.4f} m")
vz = np.gradient(zr, tv)                                             # approach speed [m/s]
print("R grows monotonically:", bool(np.all(np.diff(Rr) > 0)), "| approach speed falls monotonically:", bool(np.all(np.diff(vz) < 0)))
""", explain=r"""
1. Kelvin's self-speed of a thin ring, Γ/(4πR)[ln(8R/a) − ¼] = 0.328816 m/s.
2. The velocity a ring induces on its own axis, from the elliptic-integral formula: u_R = 0, u_z = 0.353553 m/s — the
   closed form ΓR²/2(R² + z²)^{3/2} of C08.
3. A ring (R = 0.5 m, a = 5 cm) launched at a wall 2 m away: it travels at ≈ 0.66 m/s at first, then widens from 0.5 m
   to about 0.66 m while its approach speed falls — Fig. 5.15b in numbers. At t ≈ 3.12 s the gap is down to three core
   radii and the model stops: a core that close to its image is squeezed out of shape, which the thin-core formulas
   ignore (continued blindly, they would let R grow without bound).
""")
nb.figure(r"""
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8))                  # path in the (R, z) plane | speed vs gap
a.plot(Rr, zr, color=COLORS["teal"], lw=2.5, label="ring cross-section")        # the real ring (valid samples only)
a.plot(Rr, 4.0 - zr, color=COLORS["rose"], lw=2, alpha=0.5, label="image ring (behind the wall)")   # mirror at z = 2 m
a.plot(Rr[-1], zr[-1], "o", color=COLORS["ink"], label=f"model stops, t = {res['stop_time']:.2f} s")   # gap = 3 core radii
a.axhline(2.0, color=COLORS["ink"], lw=2); a.set_xlabel("ring radius R [m]"); a.set_ylabel("z [m]"); a.legend(fontsize=8)
a.set_title("(a) the path bends outward along the wall", fontsize=10)
b.plot(2.0 - zr, vz, color=COLORS["teal"], lw=2); b.invert_xaxis()   # approach speed against the remaining gap
b.axvline(3*a_end, color=COLORS["muted"], ls="--"); b.text(3*a_end + 0.01, 0.05, "3 core radii:\nmodel stops", fontsize=8, color=COLORS["muted"])
b.set_xlabel("gap to the wall [m]"); b.set_ylabel("approach speed dz/dt [m/s]"); b.set_title("(b) it slows as it closes in", fontsize=10)
plt.show()
""", see="(Note `N59`, Fig. 5.15 analogue.) Left, the teal path rises almost straight, then bends outward toward the black wall while the faded rose image "
         "path mirrors it, ending at the black dot where the model stops; right, the approach speed stays near 0.66 m/s far "
         "away and falls steeply as the gap closes to the dashed three-core-radius line.",
    read="The image ring (opposite circulation) pushes the real ring outward and back: it widens, and a wider ring is "
         "slower (U ∝ ln(8R/a)/R), so it approaches ever more slowly and never reaches the wall — in this ideal, "
         "qualitative model.",
    change="…a fatter core (a = 0.1 m): a slower ring (ln(8R/a) smaller) that widens the same way.")
note("N43", "Leap-frogging rings", r"""
Two coaxial rings of the same sense (the book leaves it to Exercise 5.15): the front ring is widened and slowed by the
rear one, the rear ring narrowed and sped up; it passes through the front ring and the roles swap — for ever in this
thin-core model (real rings deform and break down after a few passes). The total impulse ΣΓπR² is conserved. The same thin-core model as N42 (qualitative).
""")
nb.code(r"""
tlf = np.linspace(0, 40, 401 if not FAST else 201)                   # 40 s [s]
res2 = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1), dict(R=1.0, z=0.6, Gamma=1.0, a=0.1)], tlf)   # two coaxial rings
passes = int(np.sum(np.diff(np.sign(res2["z"][:, 1] - res2["z"][:, 0])) != 0))   # how often they swap order
print(f"impulse ΣΓπR² changes by {np.ptp(res2['impulse'])/res2['impulse'][0]:.1e} (relative); pass-throughs in 40 s: {passes}")
""")
nb.animation(r"""
nf = 30 if not FAST else 16                                          # frames
idx = np.linspace(0, len(tlf) - 1, nf).astype(int)                   # which outputs to show
Rl, zl = res2["R"], res2["z"]                                        # (T, 2) radius and height of each ring [m]
with plt.rc_context(ANIM_RC):                                         # light settings for fast frames
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8), dpi=80, gridspec_kw=dict(width_ratios=[1.4, 1]))
fig.subplots_adjust(left=0.06, right=0.98, bottom=0.14, top=0.88, wspace=0.25)   # fixed margins
a.set_xlim(zl.min() - 0.5, zl.max() + 0.5); a.set_ylim(-1.6, 1.6); a.set_xlabel("z [m]"); a.set_ylabel("R [m]")
a.axhline(0, color=COLORS["muted"], ls="--", lw=0.8)                 # the common axis
cols = [COLORS["teal"], COLORS["accent"]]                             # ring 1 teal, ring 2 purple
pts = [a.scatter([], [], color=c, s=50) for c in cols]               # the two cross-sections of each ring (±R)
trl = [a.plot([], [], color=c, lw=0.8, alpha=0.6)[0] for c in cols]  # trails of the upper cross-sections
b.set_xlim(0, tlf[-1]); b.set_ylim(Rl.min() - 0.05, Rl.max() + 0.05); b.set_xlabel("t [s]"); b.set_ylabel("ring radius R [m]")
rl = [b.plot([], [], color=c, lw=2)[0] for c in cols]                # R₁(t), R₂(t)


def update(i):                                                       # draw frame i
    k = idx[i]                                                       # output index
    for m in range(2):
        pts[m].set_offsets(np.array([[zl[k, m], Rl[k, m]], [zl[k, m], -Rl[k, m]]]))   # the ring cuts the plane at ±R
        trl[m].set_data(zl[:k + 1, m], Rl[:k + 1, m])
        rl[m].set_data(tlf[:k + 1], Rl[:k + 1, m])
    a.set_title(f"t = {tlf[k]:.1f} s; impulse ΣΓπR² = {res2['impulse'][k]:.6f} m⁴/s", fontsize=9)
    return pts


with plt.rc_context(ANIM_RC):                                         # render the frames with the light settings
    show_animation(animate(update, frames=nf, fig=fig, interval=150), player="video")   # smooth movie
""")
see_read_change(
    see="Left, the meridional plane: two rings (teal and purple dots at ±R) travelling to the right, taking turns to "
        "shrink and slip through the other; right, their radii swapping in a regular rhythm while the impulse in the "
        "title does not change.",
    read="Each ring moves by its own curvature and is pushed by the other: the rear one is narrowed and sped up, the front "
         "one widened and slowed, so the rear one passes through — the thin-core model of N42, now with two rings.",
    change="…unequal strengths: the passes become irregular and may stop.")
explainer("point_vortex_lab", "A vortex cannot push itself — so how do vortices move?",
          "dropping vortices of either sign, dragging them and pressing play turns D21 and D22 into motion — the orbit "
          "about a fixed centre of vorticity, the marching pair, the image that moves in lock-step behind a wall — with "
          "the invariants as a live accuracy check.",
          "", [
              "Equal pair: time one orbit, 19.7 s, and watch G stay put.",
              "Set Γ₂/Γ₁ = 3: G moves to ¾ of the way; the orbit speeds up.",
              "Opposite pair near a wall: they march, then split along it.",
              "Click a vortex: the inspector adds up the others' pushes.",
          ])
whatif(r"""
…instead of a few vortices we had infinitely many, packed side by side along a line? Their pushes would add up to a
jump in velocity across the line — a vortex sheet (C14).
""")

# =====================================================================================================================
# A.8 §5.8 Vortex Sheet — C14
# =====================================================================================================================
nb.section("5.8", "Vortex Sheet", intro=r"""
**What is this section about?** Line up infinitely many line vortices side by side and the velocity jumps across the
line: above it the fluid moves one way, below it the other. The strength of such a vortex sheet — its circulation per
unit length — is exactly that jump. Sheets are the idealisation of every shear layer, wake and wing.
""")
core("C14", "A vortex sheet: strength = jump in tangential velocity", r"""
Where a fast stream runs beside a slow one — a river entering a lake, the jet stream over calmer air, the air leaving a
wing's trailing edge — the velocity changes almost abruptly. What is a velocity jump made of?

*In one line:* $\gamma=d\Gamma/ds=u_2-u_1$ (u₁ just above, u₂ just below, counterclockwise positive).
""")
problem(r"""
A thin layer across which the velocity jumps looks like a boundary, not like vorticity. Yet all the spin of the flow is
concentrated there: the jump *is* a sheet of vortices. Seeing it this way explains why such layers roll up into billows
(Kelvin–Helmholtz, Ch. 11) and how a wing's lift is carried by vortices (Ch. 14).
""")
idea("""
u1 = -gamma/2    <- <- <- <- <- <- <- <-        just above
    o  o  o  o  o  o  o  o  o  o  o  o         a row of counterclockwise filaments, gamma ds each
u2 = +gamma/2    -> -> -> -> -> -> -> ->        just below
circulation of a thin box ds x dn round it:  dGamma = (u2 - u1) ds   ->   strength gamma = u2 - u1
""", r"""
**A velocity jump is a sheet of vortices, and the jump is its strength.**
""")
note("N44", "A sheet of line vortices (Fig. 5.16)", r"""
infinitely many parallel filaments on a surface. The tangential velocity jumps across it; the normal velocity is
continuous (zero if the sheet does not move across itself); a real sheet has a finite thickness that spreads the jump
(N22). **Number:** filaments of total strength γ = 2 m/s per metre → just above u = −1 m/s, just below +1 m/s.
""", equation=r"[u_t]\neq0,\qquad[v_n]=0")
D("D23", ref="3.18")
nb.md(r"""
> ⚠️ **Common confusion — two signs and two Γ's.** The text writes $\Gamma\equiv\frac{d\Gamma}{ds}=u_2-u_1$
> (counterclockwise circulation per length, u₁ above, u₂ below); Fig. 5.16's caption prints $d\Gamma/ds=u_1-u_2$ — the
> clockwise sense, which is positive for its clockwise filaments. And the book reuses Γ for a quantity in m/s; we write
> γ. Our code: `vortex_sheet_strength(u_above, u_below, convention="ccw")` = u₂ − u₁ (switch `"caption"` for u₁ − u₂).
""")
nb.md(r"""
*Gloss — from a sum to a sheet.* N filaments of strength γL/N spaced L/N apart become a continuous sheet as N → ∞ — a
Riemann sum turning into an integral (Ch. 2); far from the sheet (more than a spacing away) the row already looks
continuous. We measure the gap between the two with the **L1 error** $\int\lvert u_N-u_\infty\rvert\,dy$ over a strip
|y| < 0.1 m.
""")
nb.figure(r"""
fig, ax = plt.subplots(figsize=(6, 3.4))                              # the Fig. 5.16 analogue
u1, u2, gm = sheet_circuit(ax, ds=1.0, dn=0.4, n_filaments=9, gamma=2.0)   # counterclockwise filaments, γ = 2 m/s
ax.set_title(f"u₁ = {u1:+.1f} m/s above, u₂ = {u2:+.1f} m/s below: γ = u₂ − u₁ = {gm:+.1f} m/s", fontsize=10)
plt.show()
""", see="(Note `N60`, Fig. 5.16 analogue.) A row of teal (counterclockwise) filaments on a line, a purple rectangle ds × dn round a piece of it with "
         "arrows walking it counterclockwise, and orange velocity arrows: leftward above, rightward below.",
    read="Walking the box counterclockwise: +u₂ ds along the bottom, −u₁ ds along the top, and the short sides "
         "(each just v·dn, since v is continuous) shrink away as the box is flattened, dn → 0 — D23. The title's numbers are the sheet's own velocities ∓γ/2.",
    change="…clockwise filaments (γ < 0): the orange arrows swap, u₂ − u₁ becomes negative — the caption's convention "
           "would call it positive.")
nb.worked_example("a row of 10 filaments, then 100", r"""
γ = 2 m/s on a sheet 1 m long, probe at x = 0.

1. N = 10: each filament has Γ = γ × 0.1 = 0.2 m²/s; at y = 5 cm the sum gives −0.854 m/s (the continuous finite sheet:
   −0.937 m/s).
2. N = 100 (Γ = 0.02 m²/s each): −0.936551 m/s, equal to the continuous value to 2×10⁻⁶.
3. Just above the middle of the continuous sheet: $-\gamma/2$ = −1 m/s; just below +1 m/s; jump u₂ − u₁ = 2 m/s = γ ✓.
4. Far above a finite sheet the velocity fades (the sheet looks like a single vortex of circulation γL = 2 m²/s).
""")
nb.code(r"""
print(ch05.vortex_sheet_strength(-1.0, 1.0), ch05.vortex_sheet_strength(-1.0, 1.0, convention="caption"))   # 2.0 (text), −2.0 (caption)
for N in (10, 100):                                                  # discrete rows of N filaments on [−0.5, 0.5] m
    print(N, ch05.discrete_sheet_u(0.0, 0.05, 2.0, N))                # u at (0, 5 cm) [m/s]
print(ch05.continuous_sheet_velocity(0.0, 0.05, 2.0), ch05.continuous_sheet_velocity(0.0, 1e-6, 2.0), ch05.continuous_sheet_velocity(0.0, -1e-6, 2.0))   # the continuous sheet
c = ch05.discrete_sheet_convergence(2.0)                             # L1 error of u(y) for N = 10, 100, 1000
print([f"{e:.5g}" for e in c["l1_error"]], f"observed order in 1/N: {observed_order(1/np.array(c['N']), c['l1_error']):.2f}")
""", explain=r"""
1. Both sign conventions: the text's u₂ − u₁ = +2 m/s, the caption's u₁ − u₂ = −2 m/s.
2. The discrete row at one height: −0.853907 (N = 10), −0.936551 m/s (N = 100).
3. The continuous finite sheet: −0.936549 m/s at 5 cm; ∓1.000 m/s (= ∓γ/2) just above and below the middle.
4. The L1 error of the discrete profile falls as 1/N: 0.04397, 0.004413, 0.0004413 — order 1.00.
""")
nb.md("**From scratch — a hundred line vortices, one line of numpy:**")
nb.check_agree(r"""
N, L, g = 100, 1.0, 2.0                                              # filaments, sheet length [m], γ [m/s]
xs = (np.arange(N) + 0.5)*L/N - L/2                                  # filament positions (cell midpoints) [m]
y = 0.05                                                             # probe height [m]
u = np.sum(g*L/N/(2*np.pi)*(-y)/(xs**2 + y**2))                      # Σ (Γ_j/2π)(−(y − 0))/|r|²: the x-velocity at (0, y)
print(u)                                                             # −0.936551 m/s
assert np.allclose(u, ch05.vortex_sheet_velocity(np.array([0.0, 0.05]), np.array([xs, 0*xs]), 2.0, L/N)[0], rtol=1e-12)   # = the library
""")
nb.plotly(r"""
yq = np.linspace(-0.2, 0.2, 401)                                     # heights across the sheet [m]
ucont = ch05.continuous_sheet_velocity(0.0, yq, 2.0)[0]              # the continuous finite sheet [m/s]


def row(N):                                                          # number of filaments → curves
    return {"row of N filaments": (yq, ch05.discrete_sheet_u(0.0, yq, 2.0, int(N))), "continuous sheet": (yq, ucont)}  # the discrete row and the continuous sheet


Ns = np.unique(np.geomspace(4, 1000, 12 if not FAST else 6).astype(int))   # log-spaced N
fig = slider_figure(row, "N", Ns, unit="filaments", xlabel="y [m]", ylabel="u(0, y) [m/s]",         # precompute every slider position
                    title="From a row of vortices to a velocity jump (γ = 2 m/s)")                  # the message of the figure
recolor(fig, {"row of N filaments": COLORS["teal"], "continuous sheet": COLORS["ink"]}, dashes={"continuous sheet": "dash"})  # house colours by trace name
fig.show()                                                                                          # draw it (works on the web page too)
""")
see_read_change(
    see="For small N a wiggly teal profile (each filament's own swirl shows); as N grows it sharpens into a clean step "
        "from +1 m/s below to −1 m/s above, on top of the black dashed continuous sheet.",
    read="The jump across the sheet is its strength γ = 2 m/s; within about one filament spacing of the row the "
         "discreteness shows, beyond it the row is a sheet. Far from this 1 m sheet the speed falls off, as for a single "
         "vortex of circulation γL.",
    change="…γ doubled: the step doubles; the convergence with N is unchanged.")
nb.md(r"""
**A sheet rolls up.** A periodic sheet (period 1 m, γ = 1 m/s) of 150 point vortices with a small ripple (1 cm), each
vortex moved by all the others. The kernel is smoothed over a length δ = 0.05 m (Krasny's method: $\lvert\mathbf r\rvert^2\to\lvert\mathbf r\rvert^2+\delta^2$
in the periodic kernel), which stands in for a finite sheet thickness. ⚠️ A **qualitative model**: the smoothing sets the
size of the spiral's inner turns; the ripple's growth into a spiral is the robust result.
""")
nb.animation(r"""
nf = 30 if not FAST else 16                                          # frames
ro = ch05.sheet_rollup(N=150 if not FAST else 100, gamma=1.0, amplitude=0.01, delta=0.05, t_eval=np.linspace(0, 2.0, nf))   # (T, N) positions
with plt.rc_context(ANIM_RC):                                         # light settings for fast frames
    fig, ax = plt.subplots(figsize=(7, 3.6), dpi=80)
fig.subplots_adjust(left=0.1, right=0.97, bottom=0.14, top=0.88)     # fixed margins
(ln,) = ax.plot([], [], color=COLORS["teal"], lw=1.5)                 # the sheet (a line through the vortices)
(pt,) = ax.plot([], [], "o", color=COLORS["teal"], ms=2)              # the vortices
ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.25, 0.25); ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
Gtot = ro["Gamma_each"]*ro["N"]                                      # total circulation per period [m²/s]


def update(i):                                                       # draw frame i
    ln.set_data(ro["x"][i], ro["y"][i]); pt.set_data(ro["x"][i], ro["y"][i])
    ax.set_title(f"t = {ro['t'][i]:.2f} s, total circulation per period {Gtot:.3f} m²/s", fontsize=10)
    return ln, pt


with plt.rc_context(ANIM_RC):                                         # render the frames with the light settings
    show_animation(animate(update, frames=nf, fig=fig, interval=150), player="video")   # smooth movie
""")
gloss(["Krasny smoothing and the periodic sheet kernel"])
see_read_change(
    see="The small ripple steepens, the vortices bunch at one point of each period and the sheet winds into a spiral "
        "(a cat's eye — a closed eddy, shaped like an eye, around each rolled-up lump); the total circulation in the title "
        "does not change.",
    read="The sheet's own induced velocity carries its vortices toward the places where they are already crowded — the "
         "Kelvin–Helmholtz instability (Ch. 11). The smoothing δ only stands in for a finite thickness.",
    change="…a smaller δ: the spiral winds tighter and the computation needs more points.")
explainer("vortex_sheet_rollup", "What is a velocity jump made of?",
          "setting the number of filaments and watching u(y) sharpen into a jump whose size is the sheet's strength, then "
          "rippling the sheet and pressing play to watch it roll up, links the jump, the vortex row and the instability — "
          "which needs motion.",
          "", [
              "N = 10 → 1000: the profile becomes a step of height γ.",
              "Drag the circuit's height: its circulation is γ ds once it spans the sheet.",
              "Roll-up mode, ▶: the ripple winds into cat's eyes.",
              "Switch to the caption's convention: only the sign changes.",
          ])
whatif(r"""
…the sheet were a wall? A no-slip wall with a stream above it is a vortex sheet stuck to the wall — its strength the wall
slip it cancels. Viscosity then spreads it into a boundary layer (Ch. 9).
""")

# =====================================================================================================================
# A.9 End matter — S01, S02, summary
# =====================================================================================================================
nb.pointer(r"""
S01 · 📝 Exercises 5.1–5.20 in the book: practise on them. We have written out, in our own words, the results the text
leaves to Exercises 5.4 (the line vortex's zero net viscous force, D03), 5.5 (the lock exchange, D07), 5.9 (the Green's
function, D10), 5.10 (Kelvin in a rotating frame, D19), 5.12 (Burgers' vortex, D18) and 5.18 (the centre of vorticity,
D21), and stated 5.6 (N22) and 5.11 (N23); 5.7–5.8 (Vazsonyi and Crocco) return with compressible flow in Ch. 15.
""")
nb.pointer(r"""
S02 · 📚 Literature: Lighthill (1986) for the knife-blade pair, Sommerfeld (1964) for Helmholtz's first theorem and
leap-frogging rings; for more, Saffman's *Vortex Dynamics*, and Pedlosky for the geophysical side (Ch. 13).
""")
nb.summary(
    clicked=[
        r"**C01** A vortex tube has one strength along its whole length, $\int_V\nabla\cdot\boldsymbol\omega\,dV=-\Gamma_{\text{lower end}}+\Gamma_{\text{upper end}}=0$ *(5.4)*: it can thin and spin faster but never end in the fluid.",
        r"**C02** Pressure falls toward every vortex's centre to turn its fluid: $p-p_o=\tfrac18\rho\omega^2r^2-\rho gz$ *(5.6)* makes a bowl, $p-p_\infty=-\rho\Gamma^2/8\pi^2r^2-\rho gz$ *(5.7)* a funnel; the irrotational vortex is sheared but feels no net viscous force.",
        r"**C03** Round a material loop in an ideal, barotropic, conservatively forced fluid seen from an inertial frame, $D\Gamma/Dt=0$ *(5.8)*; each broken hypothesis is a source of vorticity.",
        r"**C04** Where isopycnals cross isobars, pressure pushes through the centre while weight hangs off it: the spin-up is $\frac1{\rho^2}\nabla\rho\times\nabla p$ *(5.28)*.",
        r"**C05** Vortex lines move with the fluid, and tubes keep their strength (Helmholtz) — under Kelvin's four restrictions.",
        r"**C06** $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ *(5.13)*: vorticity is carried, stretched and tilted, and diffused — never made by pressure in a uniform fluid.",
        r"**C07** Vorticity decides velocity: $\mathbf u=\frac1{4\pi}\int\frac{\boldsymbol\omega\times(\mathbf x-\mathbf x')}{\lvert\mathbf x-\mathbf x'\rvert^3}d^3x'$ *(5.16)* — with +1/(4π) in the Green's-function form (5.14).",
        r"**C08** A thin filament's pieces push as $d\mathbf u=\frac{\Gamma\,dl}{4\pi}\mathbf e_\omega\times\frac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ *(5.17)*; an infinite line gives Γ/2πd back.",
        r"**C09** On a rotating planet with density contrasts, $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ *(5.30)*.",
        r"**C10** Stretching a vortex line spins it up, tilting turns it, $\frac{D\omega_s}{Dt}=\omega\frac{\partial u_s}{\partial s}$ … *(5.32)* — neither happens in 2-D; Burgers' vortex balances stretching against diffusion.",
        r"**C11** The absolute circulation is conserved, $D\Gamma_a/Dt=0$ *(5.33)*, so a column keeps (ζ + f)/h: squashed columns turn anticyclonic, stretched ones cyclonic.",
        r"**C12** Vortices are moved only by each other: like ones orbit their centre of vorticity at $(\Gamma_1+\Gamma_2)/2\pi h^2$, opposite ones march at $\Gamma/2\pi h$.",
        r"**C13** A wall is a mirror vortex of opposite sign: a vortex slides along it at $\Gamma/4\pi h$; a ring widens and slows near it.",
        r"**C14** A velocity jump is a vortex sheet whose strength is the jump, $\gamma=u_2-u_1$; left alone, it rolls up.",
    ],
    feeds_forward=[
        "Ch. 6: irrotational flow stays irrotational (Kelvin); images for cylinders.",
        "Ch. 8: vortex decay and the flow round a rotating cylinder.",
        "Ch. 10: vorticity–stream-function methods and vortex methods (Biot–Savart).",
        "Ch. 11: Kelvin–Helmholtz roll-up of a vortex sheet.",
        "Ch. 12: vortex stretching and the cascade; Burgers' vortex.",
        r"Ch. 13: the rotating vorticity equation $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac{\nabla\rho\times\nabla p}{\rho^2}+\nu\nabla^2\boldsymbol\omega$ *(5.30)*, absolute circulation $D\Gamma_a/Dt=0$ *(5.33)*, $(\zeta+f)/h$ → potential vorticity, Rossby waves, fronts, baroclinic instability.",
        "Ch. 14: lifting line, trailing vortices, induced drag, ground effect.",
    ],
    left_out=[
        "Vazsonyi and Crocco equations (Exercises 5.7–5.8) → Ch. 15.",
        "Vortex reconnection (named only) → Ch. 12.",
        "The viscous rotating-cylinder solution → Ch. 8.",
    ],
)

if __name__ == "__main__":
    print("equations written next to their numbers in", finalize_equations(), "markdown cells")
    bad = self_check_numbers()
    if bad:
        raise SystemExit("markdown cells citing equation numbers without maths:\n  " + "\n  ".join(bad))
    path = nb.save()
    print(f"wrote {path.relative_to(ROOT)}: {len(nb.cells)} cells, {len(nb.cores)} CORE, {len(nb.recaps)} RECAP, "
          f"{len(nb.derivations)} derivations, {len(nb.primers)} primers, {len(nb.explainers)} explainers")
