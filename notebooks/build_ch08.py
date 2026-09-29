"""Build the Chapter 8 teaching notebook: ``notebooks/ch08_laminar_flow.ipynb``.

Source of truth: ``analysis/ch08_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract
— every call goes to ``fluidpy.ch08_laminar_flow``, imported as ``ch08``, which re-exports the three new core modules
``core.laminar``, ``core.lubrication``, ``core.creeping`` and the diffusion solvers), Part D (runtime budget, FAST
sizes), Part E (prerequisite ledger → 15 primers P185–P199 and one-line reminders of earlier primers), Part F (the 33
derivations D01–D33, one move per step) and ``analysis/ch08_curation.md`` (IDs, depths, section coverage).

Decisions from ``reports/ch08_verification.md`` override the design where they differ:
  * O2 — D01's check no longer quotes a ratio: the code cell computes ν_air/ν_water from the Ch. 1 property functions
    (15.007) and prints it;
  * O3 — D32 sizes the far-field inertia as the stream advecting the disturbance, U ∂u′/∂x (the O(U²a/r²) part sits in
    the radial component of u·∇u, (3U²a/16r²)(8cos²θ − 4sin²θ)); the crossover r ~ a/Re_a is the book's order of
    magnitude, and the code computes the exact prefactor live with ``inertia_viscous_ratio`` — ½ on the axis and on the
    side line, so inertia equals friction at r ≈ 2a/Re_a = 4a/Re (coordinator's correction of the brief: the earlier
    "1/8, 10–20 a/Re_a" was wrong and is not used);
  * O4 — one engine film everywhere (ρ = 870 kg/m³, μ = 0.05 Pa s, ν = 5.75 × 10⁻⁵ m²/s, ε²Re_L = 4.35 × 10⁻³, Λ = 49);
    D27's check curls Stokes' vorticity twice with ``core.curvilinear.curl`` (the ``curlcurl_identity`` key of
    ``stokes_sphere_sympy`` holds the D30 identity), and D28's check builds the curl–curl identity for a generic A(r, θ)
    with the same tool;
  * O5 — the Hele-Shaw grid check is first order (staircase disc), so the error roughly halves per grid doubling.
Two Part F steps did not follow from the step above by their stated move alone; one step is inserted in each (the
notebook says so): D14 (substitute C₁ and C₂ before putting p − p_e over one denominator) and D19 (insert A and B into
(8.29) before the substitution ξ = 2ζ).

**Derivations are read from Part F at build time** (``part_f()`` below, the ch07 parser): goal, start, plan, tools,
assumptions, every step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word. The seven
★★★ sympy checks (D10, D13, D14, D22, D28, D30, D31) are written here, every line commented, and re-run the
construction.

Book numbers never printed (rule 9): our own worked numbers only (engine film, glycerol bead, water plates, cloud
droplet); numbers that coincide with a printed value (3.64, 0.06, 5.54) are computed and printed at 4 significant figures.

Run:  .venv/Scripts/python.exe notebooks/build_ch08.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch08.py --dump     (prints the parsed Part F derivations only)
      .venv/Scripts/python.exe notebooks/build_ch08.py --partial  (development: writes what exists, unchecked)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch08")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# Ch. 8 from analysis/ch08.md §2 (transcribed from the rendered pages p338–p374), in the corrected form where the book
# prints a slip ((8.13b) ∂p/∂y, (8.17a) with ν, (8.19) with U₀(1 − y/h)); earlier chapters from their notebooks.
# ---------------------------------------------------------------------------------------------------------------------
EQ = {
    "8.1": r"\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\nu\nabla^2\mathbf u",
    "8.2": r"\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{on the surface}}",
    "8.3": r"\mathbf t\cdot\mathbf U_s=(\mathbf t\cdot\mathbf u)_{\text{on the surface}}",
    "8.4a": r"0=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{d^2u}{dy^2}",
    "8.4b": r"0=-\frac1\rho\frac{\partial p}{\partial y}",
    "8.5": r"u(y)=\frac Uh\,y-\frac1{2\mu}\frac{dp}{dx}\,y(h-y)",
    "8.6": r"u_z(R)=\frac{R^2-a^2}{4\mu}\frac{dp}{dz}",
    "8.7": r"\tau=\mu\frac{\partial u_z}{\partial R}=\frac R2\frac{dp}{dz}",
    "8.8": r"\tau_0=\frac a2\frac{dp}{dz}",
    "8.9": r"u_\varphi(R)=AR+\frac BR",
    "8.10": r"u_\varphi=\frac{1}{R_2^2-R_1^2}\Big\{\big[\Omega_2R_2^2-\Omega_1R_1^2\big]R-\big[\Omega_2-\Omega_1\big]\frac{R_1^2R_2^2}{R}\Big\}",
    "8.11": r"u_\varphi(R)=\frac{\Omega_1R_1^2}{R}",
    "8.12": r"u_\varphi(R)=\Omega_2R",
    "8.13a": r"\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{\partial p}{\partial x}+\frac\mu\rho\Big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\Big)",
    "8.13b": r"\frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}=-\frac1\rho\frac{\partial p}{\partial y}+\frac\mu\rho\Big(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\Big)",
    "8.14": r"x^*=\frac xL,\ y^*=\frac yh=\frac{y}{\varepsilon L},\ t^*=\frac{Ut}{L},\ u^*=\frac uU,\ v^*=\frac{v}{\varepsilon U},\ p^*=\frac{p}{P_a}",
    "8.15": r"\frac{\partial u^*}{\partial x^*}+\frac{\partial v^*}{\partial y^*}=0",
    "8.16a": r"\varepsilon^2\mathrm{Re}_L\Big(\frac{\partial u^*}{\partial t^*}+u^*\frac{\partial u^*}{\partial x^*}+v^*\frac{\partial u^*}{\partial y^*}\Big)=-\frac1\Lambda\frac{\partial p^*}{\partial x^*}+\varepsilon^2\frac{\partial^2u^*}{\partial x^{*2}}+\frac{\partial^2u^*}{\partial y^{*2}}",
    "8.16b": r"\varepsilon^4\mathrm{Re}_L\Big(\frac{\partial v^*}{\partial t^*}+u^*\frac{\partial v^*}{\partial x^*}+v^*\frac{\partial v^*}{\partial y^*}\Big)=-\frac1\Lambda\frac{\partial p^*}{\partial y^*}+\varepsilon^4\frac{\partial^2v^*}{\partial x^{*2}}+\varepsilon^2\frac{\partial^2v^*}{\partial y^{*2}}",
    "8.17a": r"0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}",
    "8.17b": r"0\cong-\frac1\rho\frac{\partial p}{\partial y}",
    "8.18": r"u(x,y,t)\cong\frac1\mu\frac{\partial p(x,t)}{\partial x}\frac{y^2}{2}+Ay+B",
    "8.19": r"u\cong-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh\Big(1-\frac yh\Big)+U_h\frac yh+U_0\Big(1-\frac yh\Big)",
    "8.20": r"\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}",
    "8.21": r"u(y,t=0)=0",
    "8.22": r"u(y=0,t)=\begin{cases}0 & t<0\\ U & t\ge0\end{cases}",
    "8.23": r"u(y\to\infty,t)=0",
    "8.24": r"u/U=f\big(y/\sqrt{\nu t},\ y/Ut\big)",
    "8.25": r"u/U=F\big(y/\sqrt{\nu t}\big)\equiv F(\eta)",
    "8.26": r"-\frac\eta2\frac{dF}{d\eta}=\frac{d}{d\eta}\Big(\frac{dF}{d\eta}\Big)",
    "8.27": r"F(\eta=0)=1",
    "8.28": r"F(\eta\to\infty)=0",
    "8.29": r"F(\eta)=A\int_0^\eta\exp(-\xi^2/4)\,d\xi+B",
    "8.30": r"\frac{u}{U}=1-\mathrm{erf}\Big(\frac{y}{2\sqrt{\nu t}}\Big)",
    "8.31": r"\delta_{99}\sim3.64\sqrt{\nu t}",
    "8.32a": r"\gamma=At^{-n}F\big(\xi/\delta(t)\big)\equiv At^{-n}F(\eta)",
    "8.32b": r"\gamma=A\xi^{-n}F\big(\xi/\delta(t)\big)\equiv A\xi^{-n}F(\eta)",
    "8.33": r"u(y=0,t)=U\cos(\omega t)",
    "8.34": r"u(y\to\infty,t)=\text{bounded}",
    "8.35": r"u(y,t)=\mathrm{Re}\{e^{i\omega t}f(y)\}",
    "8.36": r"i\omega f=\nu\frac{d^2f}{dy^2}",
    "8.37": r"f(y)=A\exp\{-(i+1)y\sqrt{\omega/2\nu}\}+B\exp\{+(i+1)y\sqrt{\omega/2\nu}\}",
    "8.38": r"u=U\exp\Big\{-y\sqrt{\frac{\omega}{2\nu}}\Big\}\cos\Big(\omega t-y\sqrt{\frac{\omega}{2\nu}}\Big)",
    "8.39": r"\rho\mathbf u\cdot\nabla\mathbf u+\nabla p=\mu\nabla^2\mathbf u",
    "8.40": r"\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*=\frac{1}{\mathrm{Re}}\nabla^{*2}\mathbf u^*",
    "8.41": r"\mathrm{Re}\,(\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*)=\nabla^{*2}\mathbf u^*",
    "8.42": r"\mathrm{Re}\,(\mathbf u^*\cdot\nabla^*\mathbf u^*)=-\nabla^*p^*+\nabla^{*2}\mathbf u^*",
    "8.43": r"\nabla p=\mu\nabla^2\mathbf u",
    "8.44": r"\Big[\frac{\partial^2}{\partial r^2}+\frac{\sin\theta}{r^2}\frac{\partial}{\partial\theta}\Big(\frac1{\sin\theta}\frac{\partial}{\partial\theta}\Big)\Big]^2\psi=0",
    "8.45": r"\psi(r=a,\theta)=0",
    "8.46": r"\partial\psi(r=a,\theta)/\partial r=0",
    "8.47": r"\psi(r\to\infty,\theta)=\tfrac12Ur^2\sin^2\theta",
    "8.48": r"\psi=Ur^2\sin^2\theta\Big(\frac12-\frac{3a}{4r}+\frac{a^3}{4r^3}\Big)",
    "8.49": r"u_r=U\cos\theta\Big(1-\frac{3a}{2r}+\frac{a^3}{2r^3}\Big),\ u_\theta=-U\sin\theta\Big(1-\frac{3a}{4r}-\frac{a^3}{4r^3}\Big)",
    "8.50": r"p-p_\infty=-\frac{3\mu aU\cos\theta}{2r^2}",
    "8.51": r"D=6\pi\mu aU",
    "8.52": r"C_D=\frac{D}{\tfrac12\rho U^2\pi a^2}=\frac{24}{\mathrm{Re}}",
    "8.53": r"\frac{\psi}{Ua^2}=\Big[\frac{r^2}{2a^2}+\frac{a}{4r}\Big]\sin^2\theta-\frac{3}{\mathrm{Re}}(1+\cos\theta)\Big\{1-\exp\Big[-\frac{\mathrm{Re}}{4}\frac ra(1-\cos\theta)\Big]\Big\}",
    # earlier chapters (from their notebooks)
    "3.29": r"u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)",
    "4.5": r"\frac{d}{dt}\int_{V^*(t)}\rho\,dV+\int_{A^*(t)}\rho\,(\mathbf u-\mathbf b)\cdot\mathbf n\,dA=0",
    "4.10": r"\nabla\cdot\mathbf u=0",
    "4.39b": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u",
    "4.85": r"\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\mu\nabla^2\mathbf u",
    "4.89": r"\frac{DT}{Dt}=\kappa\nabla^2T",
    "4.100": r"x_i^*=\frac{x_i}{l},\ u_j^*=\frac{u_j}{U},\ p^*=\frac{p-p_\infty}{\rho U^2}",
    "4.107": r"C_D\equiv\frac{F_D}{\tfrac12\rho U^2A}",
    "5.1": r"u_\theta=\omega r/2\ (=\Omega r\ \text{with vorticity}\ \omega=2\Omega)",
    "5.2": r"u_\theta=\frac{\Gamma}{2\pi r}",
    "5.13": r"\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega",
    "6.2": r"\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0",
    "6.10": r"u\equiv\partial\phi/\partial x,\ v\equiv\partial\phi/\partial y",
    "6.12": r"\nabla^2\phi=0",
    "6.77": r"E^2\psi=\frac{\partial^2\psi}{\partial r^2}+\frac{\sin\theta}{r^2}\frac{\partial}{\partial\theta}\Big(\frac1{\sin\theta}\frac{\partial\psi}{\partial\theta}\Big)",
    "6.83": r"u_r=\frac1{r^2\sin\theta}\frac{\partial\psi}{\partial\theta},\ u_\theta=-\frac1{r\sin\theta}\frac{\partial\psi}{\partial r}",
    "6.86": r"\psi=\tfrac12Ur^2\sin^2\theta",
}
EQ["8.4"] = EQ["8.4a"] + r",\ " + EQ["8.4b"]
EQ["8.13"] = EQ["8.13a"]
EQ["8.16"] = EQ["8.16a"]
EQ["8.17"] = EQ["8.17a"] + r",\ " + EQ["8.17b"]
EQ["8.32"] = EQ["8.32a"]


def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ *({num})*"


def EE(*nums: str) -> str:
    """Several equations, each with its number, joined for prose."""
    return " · ".join(E(n) for n in nums)


_EQ_REF = re.compile(r"\(((?:\d)\.\d+[a-d]?)\)")


def _plain_eq_follows(after: str) -> bool:
    """Is the equation already written in plain symbols right after its number ("(8.25) η = y/√(νt)")? Then inserting
    its LaTeX again would only duplicate it."""
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


def show_eqs(text: str) -> str:
    """Write the equation next to the first mention of a book equation number "(8.25)" in a text that names it without
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
# Part F reader: the derivations, word for word (the ch07 parser)
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
    text = (ROOT / "analysis" / "ch08_design.md").read_text(encoding="utf-8")
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
            parsed.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=d["why"], plain=d["plain"]))
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
    """Insert a step the book (and Part F) skipped, so that every step follows from the one above by one stated move."""
    PF[key]["steps"].insert(before - 1, dict(did=did, tex=tex, why=why, plain=plain))


# ---------------------------------------------------------------------------------------------------------------------
# Verification decisions (reports/ch08_verification.md O2–O4), the two inserted steps, pointers from each "Check" to the
# cell that really runs it, and Starts/Results that Part F writes only as numbers or words now show their equations.
# ---------------------------------------------------------------------------------------------------------------------
# O2: compute the ratio, never quote it
pf_sub("D01", "check", "The property functions give 15.1.",
       "The code cell after the tiny example below computes the ratio from Ch. 1's property functions (Sutherland air, "
       "IAPWS water) at 20 °C and prints it.")
# O3: the far-field inertia is the stream advecting the disturbance; the crossover is an order of magnitude
pf_sub("D32", "step2.did", "Size the velocity gradient", "Size the disturbance's gradient")
PF["D32"]["steps"][1]["tex"] = r"\lvert\nabla\mathbf u'\rvert\sim\frac{Ua}{r^2},\qquad \mathbf u'=\mathbf u-U\mathbf e_x"
pf_sub("D32", "step2.why", "Differentiating a 1/r term brings another 1/r.",
       "Differentiating a 1/r term brings another 1/r; we call the part of the velocity that differs from the stream the "
       "disturbance u′.")
PF["D32"]["steps"][2]["tex"] = (r"\rho\,\mathbf u\cdot\nabla\mathbf u\sim\rho U\frac{\partial\mathbf u'}{\partial x}\sim\rho U\cdot"
                                r"\frac{Ua}{r^2}=\frac{\rho U^2a}{r^2}")
PF["D32"]["steps"][2]["why"] = (
    "Far away the velocity is the stream U plus the small disturbance u′, so the largest part of u·∇u is the stream "
    "carrying the disturbance, U ∂u′/∂x (the term Oseen keeps, N95). In (8.49) this O(U²a/r²) inertia sits in the radial "
    "component, (3U²a/16r²)(8cos²θ − 4sin²θ); the θ-component has no such term.")
PF["D32"]["steps"][2]["plain"] = "the stream sweeps the disturbance past: inertia ~ ρU²a/r²."
pf_sub("D32", "step7.why", "Set the ratio to one: that is where inertia stops being negligible.",
       "Set the ratio to one: that is where inertia stops being negligible. This is an order of magnitude: the exact "
       "prefactor is ½ on the axis and on the side line (so the crossover is at r ≈ 2a/Re_a), and it dips near θ ≈ 55°, "
       "where the leading inertia term changes sign.")
pf_sub("D32", "check", "`inertia_viscous_ratio` has log–log slope 1 ± 0.05 for r/a ∈ [50, 500]; Re_a = 0.01 → r ≈ 100a; the 10 µm "
       "droplet (Re_a = 0.008) → r ≈ 125a = 1.25 mm.",
       "`inertia_viscous_ratio` has log–log slope 1 ± 0.05 for r/a ∈ [50, 500] (the from-scratch cell of C15 measures it). "
       "Re_a = 0.01 → r of order 100a, and the 10 µm droplet (Re_a = 0.008) → r of order 125a = 1.25 mm — orders of "
       "magnitude: the code cell of C15 prints ratio ÷ (Re_a r/a) ≈ 0.5 on the axis and on the side line, so inertia "
       "equals friction at r ≈ 2a/Re_a (200a for Re_a = 0.01).")
# O4: D27 is checked directly (the curlcurl_identity key belongs to D30)
pf_sub("D27", "check", "For Stokes' solution ω_φ = −(3Ua/2r²) sin θ (D30 step 3), `stokes_sphere_sympy()[\"curlcurl_identity\"]` gives 0.",
       "For Stokes' solution ω_φ = −(3Ua/2r²) sin θ (D30 step 3), the check cell below applies the spherical curl of "
       "`core.curvilinear` twice and gets 0 in every component.")
pf_sub("D28", "check", "(check_src in C13 item 28)", "(the check cell below)")
# D14: the substitution of C₁ and C₂ is a move of its own (inserted before Part F's step 12)
pf_insert("D14", 12, "Substitute C₁ and C₂ into step 7",
          r"p-p_e=\frac{6\mu LU}{\alpha h_o^2}\Big[-\frac{1+\alpha}{(2+\alpha)(1+s)^2}+\frac{1}{1+s}-\frac{1}{2+\alpha}\Big],\ \ s=\frac{\alpha x}{L}",
          "Insert C₁ = −(1 + α)Uh₀/(2 + α) and C₂ from step 11 into step 7 and subtract p_e; the factor 6μLU/(αh₀²) comes out of "
          "all three terms. (Part F went straight to the common denominator; this step is added so that each line follows "
          "from the one above.)",
          "three simple fractions in s = αx/L.")
pf_sub("D14", "step13.why", "The numerator −(1 + α) + (2 + α)(1 + s) − (1 + s)² simplifies to s(α − s).",
       "Multiply each fraction up to the common denominator (2 + α)(1 + s)²; the numerator −(1 + α) + (2 + α)(1 + s) − (1 + s)² "
       "expands to −1 − α + 2 + 2s + α + αs − 1 − 2s − s² = αs − s² = s(α − s).")
pf_sub("D14", "check", "Numbers: L = 5 cm", "Numbers (the code cell after the tiny example of C07 prints them): L = 5 cm")
# D19: insert A and B into (8.29) before the substitution (a skipped move)
pf_insert("D19", 10, "Insert A and B into (8.29)",
          r"F(\eta)=1-\frac{1}{\sqrt\pi}\int_0^\eta e^{-\xi^2/4}\,d\xi",
          "Put B = 1 (step 6) and A = −1/√π (step 9) back into (8.29), so that the solution is written with both constants "
          "known. (Part F went straight to the substituted integral; this step is added.)",
          "the solution, still in the variable ξ.")
pf_sub("D19", "step11.why", "Same ξ = 2ζ; the upper limit becomes η/2.",
       "Substitute ξ = 2ζ again (dξ = 2dζ, ξ²/4 = ζ²): the factor 2 turns 1/√π into 2/√π, and the upper limit ξ = η becomes ζ = η/2.")
pf_sub("D02", "check", "The same steps applied in `ch08.parallel_flow_sympy(\"channel\")` return these two equations.",
       "The same steps applied in `ch08.parallel_flow_sympy(\"channel\")` (the code cell after the tiny example of C02) "
       "return these two equations.")
pf_sub("D06", "check", "parity with Ch. 3's `pipe_profile` far downstream and Ch. 4's `exact_solution(\"pipe_poiseuille\")`.",
       "parity with Ch. 3's `pipe_profile` far downstream and Ch. 4's `exact_solution(\"pipe_poiseuille\")` (the code cell after "
       "the tiny example of C03).")
pf_sub("D10", "check", "(`lubrication_nondim_sympy`, check_src in C05)", "(`lubrication_nondim_sympy`, the check cell below)")
pf_sub("D13", "check", "sympy (check_src in C06 item 29)", "sympy (the check cell below)")
pf_sub("D22", "check", "✓ (sympy)", "✓ (sympy, the check cell below)")
pf_sub("D30", "check", "sympy: both components of ∇p = μ∇²u hold (check_src in C14 item 45).",
       "sympy: both components of ∇p = μ∇²u hold (the check cell below).")

# Starts and Results that Part F writes as numbers or words: show the equations (house rule)
PF["D01"]["result"] = (r"\nu=\frac{\mu}{\rho},\qquad t_d\sim\frac{L^2}{\nu},\qquad \frac{\nu_{air}}{\nu_w}\approx15", PF["D01"]["result"][1])
PF["D02"]["start"] = (r"\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0\ \text{((4.10) in 2-D)},\qquad " + EQ["8.1"] + r"\ \text{(8.1)}",
                      PF["D02"]["start"][1])
PF["D02"]["result"] = (eq_display(["8.4a", "8.4b"]), PF["D02"]["result"][1])
PF["D04"]["start"] = (r"\mu\frac{d^2u}{dy^2}=\frac{dp}{dx}=\text{const (D03)},\qquad u(0)=0,\ \ u(h)=U", PF["D04"]["start"][1])
PF["D05"]["result"] = (r"\begin{array}{l}\text{Couette: }u=\frac{Uy}{h},\qquad \text{Poiseuille: }u_{max}=-\frac{h^2}{8\mu}\frac{dp}{dx} \\ "
                       r"\tau=\frac{\mu U}{h}-\Big(\frac h2-y\Big)\frac{dp}{dx},\qquad \text{backflow}\iff\frac{dp}{dx}>\frac{2\mu U}{h^2}\end{array}",
                       PF["D05"]["result"][1])
PF["D06"]["start"] = (r"\mathbf u=(u_R,u_\varphi,u_z)=(0,\,0,\,u_z(R)),\qquad " + EQ["8.1"] + r"\ \text{(8.1) in cylindrical form}",
                      "steady, fully developed, axisymmetric flow along a round pipe")
PF["D06"]["result"] = (eq_display(["8.6"]), PF["D06"]["result"][1])
PF["D07"]["result"] = (r"\begin{array}{l}" + EQ["8.7"] + r"\ \text{(8.7)},\qquad " + EQ["8.8"] + r"\ \text{(8.8)} \\ "
                       r"Q=-\frac{\pi a^4}{8\mu}\frac{dp}{dz},\quad V=-\frac{a^2}{8\mu}\frac{dp}{dz},\quad u_{max}=2V,\quad f=\frac{64}{\mathrm{Re}}\end{array}",
                       PF["D07"]["result"][1])
PF["D08"]["start"] = (r"\mathbf u=(u_R,u_\varphi,u_z)=(0,\,u_\varphi(R),\,0)", "steady swirl between two infinitely long concentric cylinders")
PF["D08"]["result"] = (r"\begin{array}{l}" + EQ["8.10"] + r"\ \text{(8.10)} \\ A=\frac{\Omega_2R_2^2-\Omega_1R_1^2}{R_2^2-R_1^2},\quad "
                       r"B=-\frac{(\Omega_2-\Omega_1)R_1^2R_2^2}{R_2^2-R_1^2}\end{array}", PF["D08"]["result"][1])
PF["D09"]["start"] = (EQ["8.9"] + r"\ \text{(8.9)},\qquad A=\frac{\Omega_2R_2^2-\Omega_1R_1^2}{R_2^2-R_1^2},\ \ B=-\frac{(\Omega_2-\Omega_1)R_1^2R_2^2}{R_2^2-R_1^2}",
                      "the general swirl and its constants from D08")
PF["D09"]["result"] = (eq_display(["8.11", "8.12"]), PF["D09"]["result"][1])
PF["D10"]["result"] = (eq_display(["8.15", "8.16a", "8.16b"]), PF["D10"]["result"][1])
PF["D11"]["start"] = (eq_display(["8.16a", "8.16b"]), "the scaled momentum equations of D10")
PF["D11"]["result"] = (eq_display(["8.17a", "8.17b"]), PF["D11"]["result"][1])
PF["D12"]["result"] = (EQ["8.19"] + r"\ \text{(the consistent (8.19))}", PF["D12"]["result"][1])
PF["D12"]["start"] = (eq_display(["8.17a", "8.17b"]) + r",\qquad u(x,0,t)=U_0(t),\ \ u(x,h,t)=U_h(t)",
                      "the lubrication balance, with the lower wall y = 0 moving at U₀ and the upper wall y = h(x, t) at U_h")
PF["D14"]["result"] = (r"\begin{array}{l}p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}\ \ \text{(exact)} \\ "
                       r"p-p_e\cong\frac{3\alpha\mu LU}{h_o^2}\frac xL\Big(1-\frac xL\Big),\qquad W=\frac{\alpha\mu L^2U}{2h_o^2}\ \ (\alpha\ll1)\end{array}",
                       PF["D14"]["result"][1])
PF["D14"]["start"] = (r"h(x)=h_o\Big(1+\frac{\alpha x}{L}\Big),\ \ u(x,0)=0,\ \ u(x,h)=U,\ \ p(0)=p(L)=p_e,\qquad " + EQ["4.5"] + r"\ \text{(4.5)}",
                      "a sloped pad sliding at U over a floor at rest, the same pressure at both ends, and the control-volume mass balance")
PF["D15"]["start"] = (EQ["4.5"] + r"\ \text{(4.5)},\qquad " + EQ["8.18"] + r"\ \text{(8.18)}",
                      "the control-volume mass balance on a slice of width dx, and the generic lubrication profile")
PF["D16"]["start"] = (r"u(x,y,t)\ \text{for}\ y>0,\ \ u(0,t)=U\ (t\ge0),\qquad " + EQ["4.10"] + r"\ \text{(4.10)},\qquad " + EQ["8.1"] + r"\ \text{(8.1)}",
                      "an infinite plate at y = 0 set moving at U at t = 0 under fluid at rest")
PF["D16"]["result"] = (eq_display(["8.20", "8.21", "8.22", "8.23"]), PF["D16"]["result"][1])
PF["D17"]["start"] = (r"u=u(U,y,t,\nu)\ \text{solving}\ " + EQ["8.20"] + r"\ \text{(8.20) with (8.21)–(8.23)}",
                      "the velocity can depend only on the plate speed, the height, the time and the viscosity")
PF["D17"]["result"] = (eq_display(["8.25"]), PF["D17"]["result"][1])
PF["D18"]["result"] = (eq_display(["8.26", "8.27", "8.28"]), PF["D18"]["result"][1])
PF["D19"]["start"] = (eq_display(["8.26", "8.27", "8.28"]), "the similarity ODE and its two conditions from D18")
PF["D19"]["result"] = (EQ["8.30"] + r",\qquad \mathrm{erf}(\zeta)=\frac2{\sqrt\pi}\int_0^\zeta e^{-\xi^2}d\xi\ \ \text{(8.30)}", PF["D19"]["result"][1])
PF["D20"]["result"] = (r"\delta_{99}=2\,\mathrm{erfc}^{-1}(0.01)\sqrt{\nu t}=3.643\sqrt{\nu t}\ \ \text{(8.31)}", PF["D20"]["result"][1])
PF["D21"]["start"] = (EQ["8.32a"] + r"\ \text{(8.32a) with}\ \gamma=u/U,\ \xi=y,\ \text{in}\ " + EQ["8.20"] + r"\ \text{(8.20)}",
                      "the general similarity form applied to Stokes' first problem")
PF["D21"]["result"] = (r"\delta=\sqrt{2C_1\nu t};\quad C_1=\tfrac12:\ \ " + EQ["8.26"] + r"\ \text{(8.26)}", PF["D21"]["result"][1])
PF["D22"]["start"] = (r"\omega_z=-\frac{\partial u}{\partial y},\quad u(y,0)=\begin{cases}+U & y>0\\ -U & y<0\end{cases},\qquad " + EQ["8.20"] + r"\ \text{(8.20)}",
                      "two streams sliding past each other, and the diffusion equation")
PF["D24"]["result"] = (eq_display(["8.38"]), PF["D24"]["result"][1])
PF["D25"]["start"] = (r"u=Ue^{-y/\delta_e}\cos(\omega t-y/\delta_e),\quad \delta_e=\sqrt{2\nu/\omega}\ \ \text{(8.38)}", "Stokes' second problem written with the e-folding depth")
PF["D25"]["result"] = (r"\delta_e=\sqrt{\frac{2\nu}{\omega}},\quad \frac{\lvert u\rvert_{max}}{U}\Big\rvert_{y=4\sqrt{\nu/\omega}}=e^{-2\sqrt2}=0.0591,\quad "
                       r"\phi=\frac{y}{\delta_e},\quad \frac{dy}{dt}\Big\rvert_{crest}=\sqrt{2\nu\omega}", PF["D25"]["result"][1])
PF["D26"]["result"] = (eq_display(["8.43", "4.10"]), PF["D26"]["result"][1])
PF["D29"]["result"] = (eq_display(["8.48"]), PF["D29"]["result"][1])
PF["D27"]["start"] = (EQ["8.43"] + r"\ \text{(8.43)},\qquad \partial_ip=\mu\,\partial_j\partial_ju_i\ \text{(Cartesian components)}",
                      "the Stokes equations, written component by component")
PF["D28"]["start"] = (EQ["6.83"] + r"\ \text{(6.83)},\qquad \boldsymbol\omega=\omega_\varphi\mathbf e_\varphi,\qquad -\nabla\times\nabla\times\boldsymbol\omega=0\ \text{(D27)}",
                      "Stokes' stream function, the azimuthal vorticity, and the vorticity equation of D27")
PF["D30"]["start"] = (EQ["8.43"] + r"\ \text{(8.43)},\qquad " + EQ["8.48"] + r"\ \text{(8.48)},\qquad \omega_\varphi=-\frac{E^2\psi}{r\sin\theta}\ \text{(D28)}",
                      "the Stokes equations, Stokes' stream function and the vorticity written through it")
PF["D30"]["result"] = (eq_display(["8.50"]), PF["D30"]["result"][1])
PF["D32"]["start"] = (EQ["8.49"] + r"\ \text{(8.49)},\qquad r\gg a", "Stokes' velocities far from the sphere")
PF["D32"]["result"] = (r"\frac{\text{inertia}}{\text{viscous}}\sim\frac{\rho Ua}{\mu}\frac ra=\mathrm{Re}_a\frac ra\quad\text{(order of magnitude)}",
                       PF["D32"]["result"][1])
PF["D33"]["result"] = (r"\psi_{Oseen}=\psi_{Stokes}+O\Big(\mathrm{Re}\frac ra\Big),\qquad " + EQ["8.48"] + r"\ \text{(8.48)}", PF["D33"]["result"][1])

if "--dump" in sys.argv:
    for k, d in PF.items():
        print(f"=== {k} {d['title']}  ({len(d['steps'])} steps)")
        print("  START", d["start"][0][:150], "|", d["start"][1][:80])
        for i, s in enumerate(d["steps"], 1):
            print(f"  {i}. {s['did']} :: {s['tex'][:110]}  || why={len(s['why'].split())}w")
        print("  RESULT", d["result"][0][:150], "|", d["result"][1][:80])
        print("  CHECK", d["check"][:200])
        print("  TOOLS", d["tools"])
    sys.exit(0)


# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch07.py)
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
    """One 🔁 cell reminding several tools primed in earlier chapters (knowledge/primers.md) — one sentence each. The
    first element of each item is the exact Part E concept text (the ledger check matches it)."""
    body = "\n".join(f"- **{c}** — {textwrap.dedent(t).strip()}" for c, t in items)
    head = "> 🔁 **Tools from earlier chapters used here** (see `knowledge/primers.md`)"
    nb.md(f"{head}{(' — ' + lead) if lead else ''}\n\n{body}", tags=["primer"])
    for c, _ in items:
        nb.primers.append(f"{c} (reminder)")


def note(nid: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note whose text opens with its curation id(s) in bold ("**N09 [B]**")."""
    nb.note(f"**{nid}** {textwrap.dedent(text).strip()}", equation=equation, ref=ref)


def _title_no_numbers(title: str) -> str:
    """Drop bare equation numbers from a derivation title (the heading shows the key equation itself instead)."""
    num = r"\(\d\.\d+[a-d]?(?:,\s*(?:\d\.\d+)?[a-d]?)*\)"
    t_ = re.sub(r"\s*(?:—\s*)?" + num + r"(?:\s*(?:→|–|,|and|or)\s*" + num + r")*", "", title)
    t_ = re.sub(r"\s+([,:)])", r"\1", t_)
    t_ = re.sub(r"\s*[—:,]\s*$", "", t_)
    t_ = re.sub(r":\s*→\s*", ": ", t_)
    t_ = re.sub(r"\s*→\s*$", "", t_)
    return re.sub(r"\s{2,}", " ", t_).strip()


def D(key: str, ref: str = "", check_src: str | None = None, extra_check: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then its traps as a ⚠️ callout."""
    d = PF[key]
    if ref in EQ:
        ref = f"{ref}: ${EQ[ref]}$"                     # the heading shows the key equation, not only its number
    S = show_eqs
    check = S(d["check"]) + (f" {extra_check}" if extra_check else "")
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


def idea(sketch: str, words: str = "", table: str = "") -> None:
    """The idea block: an ASCII sketch (and/or a table) plus one or two sentences; the book equations the sketch names
    by number are written out underneath it (the sketch itself is plain text)."""
    sk = textwrap.dedent(sketch).strip("\n")
    nums = [n for n in dict.fromkeys(re.findall(r"\((\d\.\d+[a-d]?)\)", sk)) if n in EQ]
    eqs = ("\n\n*The numbered equations in the sketch:* " + " · ".join(E(n) for n in nums)) if nums else ""
    body = ("```\n" + sk + "\n```" if sk else "") + eqs
    if table:
        body += ("\n\n" if body else "") + textwrap.dedent(table).strip()
    nb.md("#### The idea\n\n" + body + (f"\n\n{textwrap.dedent(words).strip()}" if words else ""))


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + textwrap.dedent(text).strip())


def explainer(slug: str, heading: str, why_static: str, tries: list[str]) -> None:
    """An embedded explainer, opening with the "why interactive rather than a static figure" sentence (ch04 lesson)."""
    why = f"**Why interactive rather than a static figure:** {textwrap.dedent(why_static).strip()}"
    nb.explainer(slug, heading=heading, why=show_eqs(why), tries=[show_eqs(t) for t in tries])


def confusion(text: str) -> None:
    nb.md("> ⚠️ **Common confusion:** " + textwrap.dedent(text).strip())


def gloss(title: str, text: str) -> None:
    """A one-paragraph gloss of a small tool used right here (no primer needed: a sentence says it)."""
    nb.md(f"*{title}.* {textwrap.dedent(text).strip()}")


# ---------------------------------------------------------------------------------------------------------------------
# Final passes: every book equation named by number in prose gets the equation written next to it.
# Only parenthesised groups count as equation mentions — "(8.25)", "(8.16a, 8.16b)", "(8.4a, b)", "(Eq. 8.5)" — so a
# number such as "8.5 m/s" is never mistaken for Eq. (8.5).
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
                new.append(src[last:e] + ", $" + r",\ ".join(EQ[n] for n in miss) + "$")
                last = e
            if new:
                units[u] = "".join(new) + src[last:]
        out = "\n---\n".join(units)
        if out != c.source:
            c.source = out
            changed += 1
    return changed


_UNI_TEX = {"−": "-", "₀": "_0", "₁": "_1", "₂": "_2", "θ": r"\theta ", "λ": r"\lambda ", "ω": r"\omega ", "±": r"\pm ",
            "ρ": r"\rho ", "σ": r"\sigma ", "η": r"\eta ", "ν": r"\nu ", "δ": r"\delta ", "√": r"\sqrt "}


def _exp_to_maths(m: re.Match) -> str:
    """'e^{−y²/4νt}' written in plain text (Part F why/in-words lines) → the maths $e^{-y^2/4\\nu t}$."""
    body = "".join(_UNI_TEX.get(ch, ch) for ch in m.group(0)).replace("²", "^2").replace("³", "^3")
    return "$" + re.sub(r"\s+\}", "}", body) + "$"


def tidy_raw_tex() -> int:
    """Plain-text exponents copied from Part F ('e^{−y²/4νt}', 't^{1/5}') would show their braces on the page: an
    exponential becomes inline maths, any other braced power becomes '^(…)'."""
    changed = 0
    for c in nb.cells:
        if c.cell_type != "markdown":
            continue
        parts = _PROSE_SPLIT.split(c.source)
        for k in range(0, len(parts), 2):
            seg = re.sub(r"(?:(?<![A-Za-z\\])[A-Za-z0-9])?(?<!\\)(?:e\^\{[^{}$]*\})+", _exp_to_maths, parts[k])
            seg = re.sub(r"\^\{([^{}$]*)\}", r"^(\1)", seg)
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


def self_check_ledger() -> list[str]:
    """Every Part E row explained by a primer (new or reminded) has a matching primer entry — the coverage tool's rule
    5, run at build time with the same 18-character match."""
    text = (ROOT / "analysis" / "ch08_design.md").read_text(encoding="utf-8")
    part = text.split("## Part E", 1)[1].split("## Part F", 1)[0]
    rows = [[x.strip() for x in ln.strip().strip("|").split("|")] for ln in part.splitlines() if ln.strip().startswith("|")]
    rows = [r for r in rows if len(r) >= 3 and not set("".join(r)) <= set("-: ") and r[0].lower() != "concept"]
    primers = [p.lower() for p in nb.primers]
    miss = []
    for concept, _first, by, *_ in rows:
        if "primer" in by.lower():
            key = re.sub(r"[`*$\\]", "", concept).lower().split("(")[0].strip()
            if key and not any(key[:18] in p or p[:18] in key for p in primers):
                miss.append(concept)
    return miss


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=r"""
Pour honey, drag a spoon through syrup, watch a cloud droplet drift down at about a centimetre a second: when viscosity
matters everywhere, the Navier–Stokes equations can often be solved exactly. The trick is always the same — a symmetry or
a thin geometry kills the nonlinear term, and what is left is a balance between pressure (or wall motion) and viscous
friction. This chapter solves that balance between plates, in pipes and between rotating cylinders; in the thin film under
a bearing and under a spreading drop; above a plate that starts or shakes, where momentum *diffuses* a distance
$\sqrt{\nu t}$ — the scale behind every boundary layer, every Ekman layer and every spin-up time of the ocean and
atmosphere; and round a sphere so small that inertia vanishes, which gives Stokes' drag $D=6\pi\mu aU$ *(Eq. 8.51)* and
the fall speed of cloud droplets, aerosols and silt.
""",
    roadmap=[
        r"§8.1 laminar vs turbulent flow and ν as a momentum diffusivity, $\mathrm{Re}=Ud/\nu\sim2000$–$3000$ (C01)",
        r"§8.2 plane Couette–Poiseuille flow $u(y)=\frac Uh\,y-\frac1{2\mu}\frac{dp}{dx}\,y(h-y)$ *(8.5)* and backflow (C02); Poiseuille pipe flow, $Q\propto a^4$ (C03); circular Couette flow $u_\varphi=AR+B/R$ *(8.9)* and its two limits (C04)",
        r"§8.3 the lubrication approximation, weighed by $\varepsilon^2\mathrm{Re}_L$ (C05); the lubrication profile and the Reynolds equation (C06); the slider bearing (C07); the thin-film (viscous gravity current) equation (C08)",
        r"§8.4 Stokes' first problem $\frac uU=1-\mathrm{erf}\big(\frac{y}{2\sqrt{\nu t}}\big)$ *(8.30)* and similarity (C09); the similarity ansatz and its exponents (C10)",
        r"§8.5 Stokes' second problem: the oscillating plate and its layer of depth $\sqrt{2\nu/\omega}$ (C11)",
        r"§8.6 creeping flow: the Stokes equations $\nabla p=\mu\nabla^2\mathbf u$ *(8.43)* (C12); Stokes' flow round a sphere (C13); Stokes drag $D=6\pi\mu aU$ *(8.51)* and settling (C14); where Stokes fails and Oseen's fix (C15)",
        r"§8.7 final remarks: where exact solutions stop",
    ],
    prerequisites=[
        "Navier–Stokes and the no-slip condition (Ch. 4 §4.5–4.10)",
        "Reynolds number and scaling (Ch. 4 §4.11)",
        "viscosity, ν = μ/ρ, plane Couette start-up (Ch. 1 §1.5)",
        "dimensional analysis (Ch. 1 §1.11)",
        "vorticity diffusion, the ideal vortex, solid-body rotation (Ch. 5 §5.1–5.4)",
        "the Gaussian vortex (Ch. 3 §3.5)",
        "the Stokes stream function and the ideal-flow sphere (Ch. 6 §6.8)",
        "complex amplitudes (Ch. 7 §7.7)",
    ],
)
for _c in nb.cells:                                     # the title cell's promise: equations are shown, not only cited
    _c.source = _c.source.replace("equations are cited by their numbers so you can follow along in your copy",
                                  "every equation is shown in full together with its number, so you can follow along in your copy")
nb.explainer_index([
    ("couette_poiseuille_backflow", "When does fluid flow backwards in a channel?",
     "C02: a line plus a parabola — wall drag and pressure add, and backflow opens once dp/dx > 2μU/h²"),
    ("lubrication_scaling", "Why can a thin film ignore inertia?",
     "C05: term bars of the scaled momentum equations reorder with ε and Re_L; only ε²Re_L must be small"),
    ("slider_bearing", "How does a film thinner than a hair carry a load?",
     "C06 C07: the same flux through a narrowing gap needs a pressure hump — and it turns into suction if the pad slides backwards"),
    ("viscous_gravity_current", "Why does a honey drop spread like t^(1/5)?",
     "C08: a nonlinear diffusion that chokes itself; every initial shape ends on one profile"),
    ("stokes_first_problem", "How can profiles at all times be one curve?",
     "C09: momentum diffuses √(νt) — plot against y/√(νt) and every profile collapses"),
    ("similarity_exponents", "Where do similarity exponents come from?",
     "C10: guess At^(−n)F(ξ/δ); the equation and one conserved quantity fix n and δ(t)"),
    ("oscillating_plate", "How deep does a shaking plate reach?",
     "C11: a 'wave' made only of diffusion — each layer lags by y/δ and shrinks as e^(−y/δ), δ = √(2ν/ω)"),
    ("stokes_sphere_flow", "What does creeping flow round a sphere look like?",
     "C13 C15: fore–aft symmetric, reaching tens of radii, until inertia returns far away: Oseen's wake"),
    ("stokes_drag_settling", "Where does 6πμaU come from?",
     "C14: one third pressure, two thirds friction; drag ∝ U makes a droplet's fall speed grow as a²"),
])
nb.setup()
nb.code(r"""
import numpy as np                                       # arrays and maths (Ch. 1 primer P03)
import sympy as sp                                       # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                          # static figures (Ch. 1 primer P01)
from scipy import integrate, special, optimize, linalg   # quadrature, erf/erfc, root finding, banded solves
from fluidpy import ch08_laminar_flow as ch08            # the tested chapter-8 module (re-exports the three new toolkits)
from fluidpy.core import laminar as LAM, lubrication as LUB, creeping as CRP, diffusion   # the new toolkits + solvers
from fluidpy import ch01_introduction as ch01            # Ch. 1: property functions, Π groups
from fluidpy import ch03_kinematics as ch03              # Ch. 3: the developing pipe profile
from fluidpy import ch04_conservation_laws as ch04       # Ch. 4: plane Poiseuille, Stokes' first problem (earlier code)
from fluidpy import ch05_vorticity_dynamics as ch05      # Ch. 5: the rotating cylinder, the diffusing vortex sheet
from fluidpy.core import navier_stokes as NS, similarity as SIM, vortices as VX   # Ch. 4–5 tools: NS terms, Re, vortices
from fluidpy.core import potential as PF, curvilinear as CL                       # Ch. 6 potential flows, curvilinear calculus
from fluidpy.core.interact import slider_figure, animate_figure, live   # plotly sliders (P17), time players, widgets (P47)
from fluidpy.core.anim import animate, ffmpeg_path       # matplotlib animations (P16); show_animation came with setup
from fluidpy.core.style import COLORS, savefig           # the house palette and a helper that saves PNGs to outputs/ch08
from tools.convergence import observed_order             # slope of log(error) vs log(step) (P13)
sys.path.insert(0, str(ROOT / "scripts"))                # the chapter's drawing helpers live in scripts/
from ch08_drawings import channel, pad, sphere, profile_arrows   # drawing only: walls, the bearing pad, a sphere, arrows
ffmpeg_path()                                            # find ffmpeg once, for the MP4 animations
MU_W, RHO_W, NU_W = 1.0e-3, 1000.0, 1.0e-6              # water: μ [Pa s], ρ [kg/m³], ν = μ/ρ [m²/s] (our default)
plt.rcParams["contour.negative_linestyle"] = "solid"     # negative ψ levels are streamlines too: draw them solid
import logging                                           # Python's message system (matplotlib reports font fallbacks through it)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless "font weight not found" notes


def recolor(fig, colors, dashes=None):                   # give plotly traces the notebook's colours by trace name
    for tr in fig.data:                                  # every trace of every slider step
        if tr.name in colors:                            # a name we assigned a colour to
            tr.line.color = colors[tr.name]              # same colour meaning as in the matplotlib figures
        if dashes and tr.name in dashes:                 # optional dash pattern ("dash", "dot")
            tr.line.dash = dashes[tr.name]
    return fig                                           # the same figure, restyled


def step_titles(fig, titles):                            # give each slider position of a slider_figure its own title
    old = fig.layout.sliders[0].steps                    # the steps slider_figure made (they switch curves on and off)
    fig.layout.sliders[0].steps = [dict(method=s.method, label=s.label, args=[s.args[0], {"title.text": t}])
                                   for s, t in zip(old, titles)]   # same switch + a new title text per step
    fig.layout.title.text = titles[fig.layout.sliders[0].active]   # the title of the step shown first
    return fig                                           # the same figure, now with numbers that follow the slider


print(len([n for n in dir(ch08) if not n.startswith("_")]), "public names in fluidpy.ch08_laminar_flow")
print("ch08.channel_flow is LAM.channel_flow:", ch08.channel_flow is LAM.channel_flow)   # the re-export
""", explain=r"""
1. Numerical, symbolic and plotting libraries (all primed in Ch. 1); from `scipy` we use `integrate` (`quad`,
   `cumulative_trapezoid`, `solve_bvp`, `solve_ivp`), `special` (`erf`, `erfc` and their inverses), `optimize` (`brentq`,
   `minimize_scalar`) and `linalg` (`solve_banded`) — each is reminded or primed where it is first used.
2. `ch08` is the chapter module. It **re-exports the three new toolkits** `core.laminar` (exact parallel and unsteady
   flows), `core.lubrication` (thin gaps and films) and `core.creeping` (low-Reynolds-number flow), plus the diffusion
   solvers, so `ch08.channel_flow` and `LAM.channel_flow` are the same function (the last line prints `True`). Every
   function cites its § and equation and is tested in `tests/test_ch08.py`.
3. Earlier chapters' modules are used in the recaps and for parity checks (the Ch. 3 pipe profile, the Ch. 4 Poiseuille and
   Stokes-problem functions, the Ch. 5 rotating cylinder and vortex sheet, Ch. 6's potential-flow sphere).
4. `scripts/ch08_drawings.py` only draws (channel walls, the pad, a sphere, profile arrows); it computes no physics.
5. `MU_W, RHO_W, NU_W`: water at room temperature, our default fluid (numbers ours, not the book's).
6. `recolor` only restyles plotly figures.
""")
nb.md(r"""
### ⚠️ Conventions in this chapter (read once; each is repeated where it bites)

| Symbol | Meaning here | Before (and where it bites) |
|---|---|---|
| dp/dx | the book's pressure gradient: **favourable < 0** (pressure falls downstream) | Ch. 4's code used G = −dp/dx; new functions take `dpdx` (and accept `G=`) |
| y = 0, y = h | channel walls: y = 0 **fixed**, y = h moving at U | Ch. 1's Couette start-up moved the upper plate too |
| R, r | **R** = cylindrical radius (§8.2), **r** = spherical radius (§8.6) | Ch. 3–6 used r for both |
| θ | in §8.6 measured from the **downstream** axis: rear stagnation point θ = 0, front θ = π | Ch. 6's sphere also had θ from the axis; check the direction |
| ω | **vorticity** in §8.1, §8.4, §8.6 ($\omega_z$, $\omega_\varphi$) but the plate's **frequency** [rad/s] in §8.5 | as in Ch. 7 (frequency) and Ch. 5 (vorticity) |
| η | the **similarity variable** $\eta=y/\sqrt{\nu t}$ *(8.25)*; the book's figures plot $y/(2\sqrt{\nu t})=\eta/2$ | the surface elevation in Ch. 7 |
| ε | the fineness ratio h/L of a thin gap | a density ratio in Ch. 7, the dissipation rate in Ch. 4 |
| Λ | the bearing number $\mu UL/(P_ah^2)$ | new |
| δ₉₉, δ_e, δ(t) | 99 % thickness of Stokes' first problem · e-folding depth of the Stokes layer · the width in a similarity form | three different thicknesses |
| Re | pipe $Ud/\nu$ (diameter, mean speed) · lubrication $\mathrm{Re}_L=\rho UL/\mu$ (passage length) · generic $\rho UL/\mu$ · sphere $\mathrm{Re}=2aU/\nu$ (**diameter**) | the radius-based $\mathrm{Re}_a=\mathrm{Re}/2$ appears in C15 |
| D | an integration constant (C13) **and** the drag (C14) | say which |
| g | `ch08`'s settling and film functions use the SI standard g = 9.80665 m/s² (`G0`) | Ch. 4, 5, 7 used 9.81; the book prints no g-dependent number here |

We compute with the book's conventions and name the other one wherever it differs.

**Colours** (one meaning each, also in the explainers): pressure and the pressure-driven (Poiseuille) part **orange** ·
viscous friction and the wall-driven (Couette) part **rose** · inertia **teal** · velocity profiles **blue** · vorticity
**purple** · backflow and warnings **amber** · rescaled (similarity) curves purple dashed · an analytic ghost grey dashed ·
numerical dots teal · a printed slip grey dotted, labelled "as printed".
""")
nb.md(r"""
**Book slips we correct in this notebook** (each is explained where it is used):

| The book prints | We use | Where |
|---|---|---|
| §8.1 calls μ "the kinematic viscosity" | μ is the **dynamic** viscosity; ν = μ/ρ is the kinematic one | C01 |
| channel mean velocity written $V=\int_0^hu\,dy$ in the middle expression | $V\equiv Q/h$ (the 1/h is missing; units m²/s ≠ m/s) | C02 (N10) |
| (8.13b) with $-\frac1\rho\frac{\partial p}{\partial x}$ | $-\frac1\rho\frac{\partial p}{\partial y}$ in the y-momentum equation | C05 (N24, D10) |
| (8.17a) as $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\frac{\partial^2u}{\partial y^2}$ | $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ (units) | C05 (D11) |
| "(8.16a) can be integrated twice" | it is the simplified (8.17a) that is integrated | C06 (N33) |
| (8.19) ending in $+U_0$ | $+U_0\big(1-\frac yh\big)$ (otherwise u(h) = U_h + U₀) | C06 (N34, D12) |
| Example 8.1: intermediate integrals with $(1-\alpha x/L)$; final pressure with the denominator to the **first** power | $(1+\alpha x/L)$; the denominator **squared** | C07 (D14) |
| Example 8.2's lubrication equations without ν and walls at "y = 0, h" | with ν; the gap coordinate is z | C06 (N36) |
| "the final equation of Example 8.2" (in Example 8.7) | Example 8.3 | C10 (N63) |
| $\int_0^\infty\omega\,dy=-U$ | $+U$ (u falls upward, so ω > 0) | C09 (N53) |
| Example 8.5's ±2.76 | ±2.772 (its own width 5.54 is right) | C10 (N59) |
| "substitution of (8.33) into (8.20)" | it is (8.35) that is substituted | C11 (N68) |
| rear pressure minimum +3μU/2a | **−**3μU/2a | C14 (N88, D30) |
| "(9.63)" and "(9.68)" in §8.6 | (8.43) and (8.48) | C13 (N93), C15 (N98) |
| Oseen's equation with $+\partial p/\partial x_i$ | $-\partial p/\partial x_i$ | C15 (N95) |
| power into the fluid $(2\pi R_1)\tau_{R\varphi}u_\varphi$ | $-2\pi R_1\sigma_{R\varphi}u_\varphi$ (positive) | C04 (R11) |
| the bearing number Λ "near unity" with p* = p/P_a | Λ ≈ 49 for our film (≈ 10³ for thinner, longer films); the natural scale is μUL/h² | C05 (N31) |
""")
nb.md(r"""
> 🔁 **Tools from earlier chapters used in this one** (one line each where first used; the full primers are in
> `knowledge/primers.md`): partial derivative (P25), first-order Taylor (P26), definite integral (P27), the linear
> second-order ODE (P44), separation of variables (P42), chain rule (P49/P91), exponent rules (P43), power laws and
> log–log plots (P13), limits and orders of smallness (P68), order-of-magnitude scaling (P130), scaled variables and the
> chain rule (P133), product rule (P38), fundamental theorem of calculus (P84), differentiation under the integral sign
> (P109), substitution in an integral (P106), erfc (P123), Euler's formula (P45), complex square roots (P159), complex
> amplitudes (P176), curl of a curl (P122), Schwarz's theorem (P121), surface integrals on a sphere (P163), line integrals
> (P35), cylindrical and spherical unit vectors (P88), polar coordinates as a moving basis (P105), boundary conditions
> (P20), `np.linalg.solve` (P57), sympy (P40) and `expand`/`series` (P117), `quad` (P87), `brentq` (P108), `expm1` (P107),
> `minimize_scalar` (P170), `meshgrid` (P76), contour and streamline plots (P78), finite differences (P21), explicit
> stepping (P30), the trapezoid rule (P37), `np.random.default_rng` (P10), `animate` (P16), `slider_figure` (P17),
> `show_viz` (P18), live widgets (P47), `assert np.allclose` (P15). New primers of this chapter: P185–P199.
""")

# =====================================================================================================================
# A.1 §8.1 — R01, R02, C01 (N01, N02, N103, D01, P185), R03–R06
# =====================================================================================================================
nb.section("8.1", "Introduction", intro=r"""
**What is this section about?** Where this chapter lives: flows in which viscosity acts everywhere, not just in thin
layers. We meet Reynolds's experiment (smooth versus chaotic flow in a pipe), learn to read the kinematic viscosity ν as a
*diffusivity* of momentum and vorticity, and restate the equations and wall conditions every solution in the chapter must
satisfy.
""")
nb.recap("R01", "The Reynolds number as inertia over viscosity", r"""
Scaling the Navier–Stokes equation gives inertia $\sim U^2/L$ and viscous acceleration $\sim\mu U/(\rho L^2)$; their ratio
$\frac{U^2/L}{\mu U/\rho L^2}=\frac{\rho UL}{\mu}=\mathrm{Re}$ is the Reynolds number. Earlier chapters dropped viscosity
where Re ≫ 1 (ideal flow, Ch. 6); in this chapter viscosity is kept everywhere.
""", where="Ch. 4 §4.11")
nb.code(r"""
s = ch08.inertia_viscous_scales(0.1, 0.01, NU_W)         # U = 10 cm/s, L = 1 cm, water: the two accelerations [m/s²]
print({k: round(v, 6) for k, v in s.items()})            # inertia U²/L, viscous νU/L², and their ratio
print("Re =", round(SIM.reynolds_number(0.1, 0.01, nu=NU_W), 3))   # Ch. 4's Re = UL/ν — the same ratio
""")
nb.recap("R02", "Vorticity diffuses with ν", r"""
In a plane flow the vorticity equation $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2
\boldsymbol\omega$ *(5.13)* loses its stretching term, because $\boldsymbol\omega=\omega_z\mathbf e_z$ and nothing depends on
z, so $(\boldsymbol\omega\cdot\nabla)\mathbf u=\omega_z\,\partial\mathbf u/\partial z=0$; what is left is
$D\omega_z/Dt=\nu\nabla^2\omega_z$ — the same form as the heat equation. A fluid with larger ν smooths out a vorticity
pattern faster.
""", where="Ch. 5 §5.4")
core("C01", "Laminar vs turbulent flow, and ν as a momentum diffusivity", r"""
When does a flow stay smooth — and what does the kinematic viscosity ν actually measure?
""")
nb.md(r"*In one line:* $\mathrm{Re}=\frac{Ud}{\nu}\sim2000$–$3000$ separates smooth from chaotic pipe flow, and $\nu=\mu/\rho$ is the "
      r"diffusivity of momentum: motion spreads a distance $\sqrt{\nu t}$ in a time t.")
problem(r"""
Open a tap a little and the stream is glassy; open it fully and it turns frothy. Reynolds put a thread of dye into water
flowing through a glass tube and saw the same thing: at low flow the dye stayed a straight thread, at high flow it broke up
and filled the tube. Every exact solution in this chapter is of the first, *laminar* kind — so we need to know when it
applies, and what sets how fast viscosity spreads the motion.
""")
idea(r"""
laminar  (Re < 2000):   ─────────────── dye stays a thread; layers slide past each other
turbulent (Re > 3000):  ~~~≈≈≈∿∿∿≈≈~~~ dye fills the tube; eddies mix it
Re = U d / ν   (U = mean speed, d = diameter, ν = μ/ρ)
""", table=r"""
| quantity | water (20 °C) | air (20 °C) |
|---|---|---|
| μ [Pa s] | 1.0 × 10⁻³ | 1.8 × 10⁻⁵ |
| ρ [kg/m³] | 1000 | 1.2 |
| ν = μ/ρ [m²/s] | 1.0 × 10⁻⁶ | 1.5 × 10⁻⁵ |
| time to diffuse 1 cm, L²/ν | 100 s | 6.7 s |
""", words=r"""
**ν, not μ, says how fast motion spreads** — air is 55 times less viscous than water but spreads momentum 15 times faster.
(Rounded textbook values in the table; the code below computes them from Ch. 1's property functions.)
""")
confusion(r"""
μ (dynamic viscosity, Pa s) is the friction coefficient in $\tau=\mu\,du/dy$; ν = μ/ρ (kinematic viscosity, m²/s) is the
diffusivity. The book's §8.1 once calls μ "the kinematic viscosity" — it means the dynamic one.
""")
remind([
    ("dynamic vs kinematic viscosity μ, ν = μ/ρ", r"Newton's law $\tau=\mu\,du/dy$ defines μ as the stress per unit shear rate (Ch. 1 P05, stress); dividing by ρ gives ν."),
])
P("P185", "diffusivity and the diffusion time L²/ν", r"""
A diffusivity D (units m²/s) says how fast something spreads by random molecular exchange: in a time t it spreads a
distance of order $\sqrt{Dt}$, so crossing a distance L takes about $L^2/D$. Heat has $\kappa=k/\rho C_p$, dye its molecular
diffusivity, and momentum has ν. Double the distance and the time quadruples.
""", code=r"""
L = [0.001, 0.01, 0.1]                   # distances [m]: 1 mm, 1 cm, 10 cm
nu = 1.0e-6                              # water's kinematic viscosity [m²/s]
print([l**2/nu for l in L])              # diffusion times L²/ν [s]: [1.0, 100.0, 10000.0]
""")
note("N01 [B] · Heat and momentum spread the same way.", r"""
The Boussinesq heat equation below, with the thermal diffusivity $\kappa\equiv k/\rho C_p$, has the same form as
$\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\nu\nabla^2\mathbf u$ *(8.1)*; both κ and ν are in m²/s, so ν is the *momentum
diffusivity*. The analogy is not complete: the pressure gradient in (8.1) also moves momentum around, so velocity is not
simply diffused.
""", equation=EQ["4.89"], ref="4.89")
D("D01")
nb.worked_example("water in a 1 cm tube at 10 cm/s", r"""
1. $\mathrm{Re}=Ud/\nu=0.1\times0.01/10^{-6}=1000$ — below 2000, laminar.
2. Double the speed: Re = 2000, the bottom of the transition band.
3. Diffusion across the tube: $t=d^2/\nu=10^{-4}/10^{-6}=100$ s — much longer than the 0.1 s the water needs to travel
   1 cm.
4. Same tube with air at 10 cm/s: $\mathrm{Re}=0.1\times0.01/1.5\times10^{-5}=67$, $t=6.7$ s.
""")
remind([
    ("Python dictionaries (fluidpy results)", "fluidpy returns several named results at once as a dictionary: `s[\"ratio\"]` (Ch. 1 P23)."),
    ("f-strings", "`f\"{x:.3g}\"` puts a number with 3 significant figures into a sentence (Ch. 1 P04)."),
    ("tuple unpacking", "`Re, label = ch08.pipe_flow_regime(...)` gives the two returned values their own names (Ch. 1 P14)."),
    ("assert np.allclose", "stops the notebook if two results differ beyond a tolerance — our proof that a from-scratch version matches the library (Ch. 1 P15)."),
])
nb.code(r"""
for U in (0.1, 0.25, 0.5):                               # three mean speeds [m/s] in a 1 cm tube of water
    Re, label = ch08.pipe_flow_regime(U, 0.01, NU_W)     # Re = Ud/ν and the 2000–3000 band
    print(f"U = {U} m/s: Re = {Re:.0f} → {label}")
nu_air = ch08.momentum_diffusivity("air", 293.15)        # ν of air at 20 °C [m²/s] (Sutherland μ, ideal-gas ρ)
nu_water = ch08.momentum_diffusivity("water", 293.15)    # ν of water at 20 °C [m²/s] (IAPWS μ and ρ)
print(f"ν_air = {nu_air:.4g} m²/s, ν_water = {nu_water:.4g} m²/s, ratio = {nu_air/nu_water:.4g}")   # D01's ratio
print(f"time to diffuse 1 cm: water {ch08.diffusion_time(0.01, nu_water):.1f} s, air {ch08.diffusion_time(0.01, nu_air):.2f} s")
""", explain=r"""
1. `pipe_flow_regime` computes $\mathrm{Re}=Ud/\nu$ and applies the 2000–3000 band: 1000 laminar, 2500 transitional, 5000
   turbulent.
2. `momentum_diffusivity` takes μ and ρ from Ch. 1's property functions (Sutherland's law for air, IAPWS for water) and
   returns ν = μ/ρ. The ratio printed here — **15.01** — is the number D01's step 4 estimated by hand as 15.
3. `diffusion_time` is $L^2/\nu$: about 100 s for a centimetre of water, 6.6 s for a centimetre of air.
""")
nb.md("**From scratch — the same regime labels and diffusion time by hand**, then a check that the library does exactly this:")
nb.check_agree(r"""
Us = (0.1, 0.25, 0.5)                                    # the same three speeds [m/s]
mine = [U*0.01/NU_W for U in Us]                         # Re = U d / ν with d = 1 cm
labels = ["laminar" if R < 2000 else ("turbulent" if R > 3000 else "transitional") for R in mine]   # the band
print(mine, labels)
assert np.allclose(mine, [ch08.pipe_flow_regime(U, 0.01, NU_W)[0] for U in Us], rtol=1e-12)   # same Re
assert labels == [ch08.pipe_flow_regime(U, 0.01, NU_W)[1] for U in Us]                         # same verdicts
assert np.allclose(0.01**2/nu_air, ch08.diffusion_time(0.01, nu_air), rtol=1e-12)              # L²/ν by hand
""")
remind([
    ("np.random.default_rng (dye sketch)", "`rng = np.random.default_rng(1)` makes reproducible random numbers; we use it to draw a schematic chaotic dye streak (Ch. 1 P10)."),
    ("matplotlib figures", "`fig, ax = plt.subplots()`, `ax.plot`, labels with units; every figure is followed by What you see / How to read it / What would change if (Ch. 1 P01)."),
    ("power laws and log–log plots", "on log–log axes a power law $t\\propto L^2$ is a straight line of slope 2 (Ch. 1 P13)."),
])
nb.figure(r"""
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.4), gridspec_kw=dict(width_ratios=[1.3, 1]))   # two panels
x = np.linspace(0, 1, 400)                               # position along the tube [m]
for yc, U, lab in ((0.6, 0.1, "laminar"), (-0.6, 0.5, "turbulent")):     # two tubes, one above the other [schematic units]
    a1.plot([0, 1], [yc + 0.3, yc + 0.3], color=COLORS["ink"], lw=1.5)   # upper tube wall
    a1.plot([0, 1], [yc - 0.3, yc - 0.3], color=COLORS["ink"], lw=1.5)   # lower tube wall
    Re, verdict = ch08.pipe_flow_regime(U, 0.01, NU_W)   # the regime from the real Re of a 1 cm water pipe
    if verdict == "laminar":                             # a straight dye thread on the axis
        a1.plot(x, yc + 0*x, color=COLORS["blue"], lw=2)
    else:                                                # a seeded random walk that fills the tube (schematic only)
        rng = np.random.default_rng(1)                   # reproducible "chaos"
        yy = np.clip(np.cumsum(rng.normal(0, 0.03, x.size)) * (x > 0.25), -0.28, 0.28)   # starts to wander at x = 0.25 m
        a1.plot(x, yc + yy, color=COLORS["blue"], lw=1)
    a1.text(0.02, yc + 0.34, f"Re = {Re:.0f}: {verdict}", fontsize=9)    # label with the computed Re
a1.set_xlim(0, 1); a1.set_ylim(-1.1, 1.2); a1.set_yticks([])            # schematic: no vertical scale
a1.set_xlabel("distance along the tube [m]"); a1.set_title("Reynolds's dye (schematic — not a simulation)", fontsize=10)
L = np.geomspace(1e-5, 1.0, 50)                          # distances from 10 µm to 1 m [m]
for nu, c, name in ((nu_water, COLORS["blue"], "water"), (nu_air, COLORS["teal"], "air")):
    a2.loglog(L, ch08.diffusion_time(L, nu), color=c, lw=2, label=f"{name}, ν = {nu:.2g} m²/s")   # t = L²/ν
    a2.plot(0.01, ch08.diffusion_time(0.01, nu), "o", color=c)           # the 1 cm point
a2.loglog(L, 1e3*L**2, color=COLORS["muted"], ls="--", lw=1, label="slope 2")   # guide line t ∝ L²
a2.set_xlabel("distance L [m]"); a2.set_ylabel("diffusion time L²/ν [s]"); a2.legend(fontsize=8)
a2.set_title("Momentum spreads as √(νt): time ∝ L²", fontsize=10)
fig.suptitle("Smooth below Re ≈ 2000; momentum spreads as √(νt)", fontsize=11)
savefig(fig, "ch08", "c01_reynolds_diffusion"); plt.show()
""", see=r"""(a) a straight dye thread at Re = 1000 and a mixed tube at Re = 5000 — our synthetic sketch of Reynolds's
apparatus `N103`, labelled schematic; (b) two parallel lines of slope 2, air's 15 times lower than water's.""",
    read=r"""On (b) read the time to diffuse a distance: 1 cm takes about 100 s in water and 6.6 s in air (the dots); the slope 2
says ten times farther takes a hundred times longer.""",
    change=r"""…the fluid were glycerine (ν ≈ 10⁻³ m²/s): its line drops by 1000 below water's — 1 cm in 0.1 s — and Re at
10 cm/s in the 1 cm tube is 1, deeply laminar.""")
note("N02 [C] · Where this leads.", r"""
When viscosity matters only in thin layers next to walls, we get boundary layers (Ch. 9); whether these laminar flows
survive small disturbances is Ch. 11's question; fully turbulent flow is Ch. 12; viscous layers in a rotating fluid — the
Ekman layers of the ocean and atmosphere — are Ch. 13. The $\sqrt{\nu t}$ spreading of this chapter reappears in all four.
""")
nb.recap("R03", "The Navier–Stokes equation for constant ρ and μ", r"""
Every solution below satisfies continuity $\nabla\cdot\mathbf u=0$ *(4.10)* and
$\frac{D\mathbf u}{Dt}=-\frac1\rho\nabla p+\nu\nabla^2\mathbf u$ *(8.1)*, which is Ch. 4's
$\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\mu\nabla^2\mathbf u$ *(4.85)* divided by ρ (gravity absorbed into p, R06). We check
each one with Ch. 4's term calculator: plug the velocity and pressure in and the residual must vanish.
""", where="Ch. 4 §4.6")
remind([
    ("functions as arguments and lambda", "`lambda x, t=0.0: …` is a one-line function; Ch. 4's calculator takes the velocity and pressure *as functions* and differentiates them numerically (Ch. 1 P29)."),
    ("component-first array layout (x[0], x[1])", "a point is an array whose first index is the component: `x[0]` is x, `x[1]` is y, `x[2]` is z (Ch. 2 P76, the project grid layout)."),
])
nb.code(r"""
u_fn = lambda x, t=0.0: np.stack([ch08.channel_flow(x[1], 0.01, U=0.1, dpdx=-4.0), 0*x[0], 0*x[0]])   # C02's channel flow u(y)
p_fn = lambda x, t=0.0: -4.0 * x[0]                      # a linear pressure p = −4x [Pa]: dp/dx = −4 Pa/m
terms = NS.ns_incompressible_terms(u_fn, p_fn, np.array([0.3, 0.004, 0.0]), mu=MU_W, g=(0.0, 0.0, 0.0))   # at one point
print("advective ", terms.advective)                    # u·∇u [m/s²]: exactly zero (C02 explains why)
print("pressure  ", terms.pressure)                     # −(1/ρ)∂p/∂x = +4 × 10⁻³ m/s²
print("viscous   ", terms.viscous)                      # ν ∂²u/∂y² = −4 × 10⁻³ m/s²
print(f"residual   max |sum| = {np.max(np.abs(terms.residual)):.1e} m/s²")   # their sum: ≈ 0 (round-off)
""", explain=r"""
1. `u_fn` is the channel flow of C02 written as a function of position (component-first: `x[1]` is y); `p_fn` is a pressure
   falling 4 Pa per metre.
2. `ns_incompressible_terms` (Ch. 4) evaluates every term of (8.1) at the point (0.3, 0.004, 0) m by finite differences.
3. The pressure term $+4\times10^{-3}$ m/s² and the viscous term $-4\times10^{-3}$ m/s² cancel; the advective term is exactly
   zero; the residual is at round-off (10⁻¹⁵). C02 derives this profile; here we only see that it passes the test.
""")
nb.recap("R04", "No flow through a wall", r"""
At a solid surface moving with velocity $\mathbf U_s$ the normal velocities agree:
$\mathbf n\cdot\mathbf U_s=(\mathbf n\cdot\mathbf u)_{\text{on the surface}}$ *(8.2)*, with **n** the unit normal (length 1,
perpendicular to the surface).
""", where="Ch. 4 §4.10")
remind([
    ("unit normal and tangent vectors", "a unit vector has length 1; the normal **n** points out of the wall, a tangent **t** lies along it, and $\\mathbf n\\cdot\\mathbf u$ picks the component of **u** along **n** (Ch. 2 P65, projection)."),
])
nb.recap("R05", "No slip", r"""
The tangential velocities agree too: $\mathbf t\cdot\mathbf U_s=(\mathbf t\cdot\mathbf u)_{\text{on the surface}}$ *(8.3)* —
viscous fluid sticks to walls. Every constant of integration in this chapter is fixed by (8.2) and (8.3).
""", where="Ch. 4 §4.10")
nb.code(r"""
n = np.array([0.0, 1.0, 0.0])                            # unit normal of a floor y = 0 (pointing into the fluid)
wall = np.array([0.1, 0.0, 0.0])                         # the floor slides at 0.1 m/s along x
print(ch08.wall_bc_residuals(np.array([0.1, 0, 0]), wall, n))      # fluid moves with the wall: (normal, tangential) = (0, 0)
print(ch08.wall_bc_residuals(np.array([0.08, 0.01, 0]), wall, n))  # a leaky, slipping wall: both residuals non-zero [m/s]
""", explain=r"""
`wall_bc_residuals(u_wall, U_s, n)` returns how much (8.2) and (8.3) are violated: the normal mismatch
$\mathbf n\cdot(\mathbf u-\mathbf U_s)$ and the size of the tangential mismatch. The first call obeys both; the second has
fluid leaking through (0.01 m/s) and slipping (0.02 m/s).
""")
nb.recap("R06", "Standing assumptions", r"""
Constant density, an inertial frame, and gravity absorbed into the pressure ($p\to p+\rho gz$) whenever there is no free
surface — only the spreading film (C08) and the settling sphere (C14) bring g back.
""", where="Ch. 4 §4.9")

# =====================================================================================================================
# A.2 §8.2 — R07, R08, C02 (N03–N11, D02–D05, IF1, E1), R09, C03 (N12–N16, D06, D07, P186, IF2), R10–R12,
#            C04 (N17–N22, N104, D08, D09, P187, IF3)
# =====================================================================================================================
nb.section("8.2", "Exact Solutions for Steady Incompressible Viscous Flow", intro=r"""
**What is this section about?** Three flows where symmetry removes the nonlinear term $\mathbf u\cdot\nabla\mathbf u$
exactly, so the Navier–Stokes equation becomes a linear ODE that we can integrate by hand: flow between parallel plates
driven by a moving wall and a pressure gradient, flow in a round pipe, and flow between two rotating cylinders. The same
short argument — continuity kills one velocity, "a function of x equals a function of y" makes the pressure gradient
constant — solves all three.
""")
nb.recap("R07", "Plane Couette flow", r"""
A fluid between a fixed plate and a plate sliding at U, with no pressure gradient, moves with the straight profile
$u=Uy/h$ and a uniform shear stress $\tau=\mu\,du/dy=\mu U/h$ (Ch. 1). C02 shows it is one half of a two-part formula.
""", where="Ch. 1 §1.5")
nb.recap("R08", "Plane Poiseuille flow", r"""
With both plates at rest and a pressure gradient, the profile is the parabola $u(y)=-\frac1{2\mu}\frac{dp}{dx}y(h-y)$, the
stress $\tau=\mu\frac{du}{dy}=-\big(\frac h2-y\big)\frac{dp}{dx}$ is linear and its magnitude at both walls is
$\frac h2\lvert dp/dx\rvert$. Ch. 4's `plane_poiseuille` takes $G=-dp/dx$.
""", where="Ch. 4 §4.6")
nb.code(r"""
y = np.linspace(0, 0.01, 5)                              # five heights across a 1 cm channel [m]
print(ch04.plane_poiseuille(y, G=2.0, h=0.01))           # Ch. 4: G = −dp/dx = +2 Pa/m  → u [m/s]
print(ch08.channel_flow(y, 0.01, U=0.0, dpdx=-2.0))      # Ch. 8: the book's dp/dx = −2 Pa/m, fixed walls → the same u
""")
confusion(r"""
the book's $dp/dx=-2$ Pa/m (pressure falling downstream, a *favourable* gradient) is Ch. 4's $G=+2$ Pa/m. Same flow,
opposite sign of the argument — the two identical rows above pin it.
""")
core("C02", "Plane Couette–Poiseuille flow", r"""
What profile results from a moving wall plus a pressure gradient — and when does fluid near the fixed wall flow backwards?
""", eqs=("8.5",))
problem(r"""
A belt drags oil along a channel while a pump pushes back. Near the belt the oil is carried forward; near the fixed floor
the pump may win and push it backwards. Journal bearings, extrusion dies, the wind-driven surface layer of a lake with a
return current underneath — all are this flow. We want the exact profile and the condition for the reversed current.
""")
idea(r"""
y=h ═══════════► U        y=h ═════════        y=h ═══════════► U
     ─────────►                 ──►                    ─────────►
     ──────►        +           ────►        =         ──────►
     ───►                       ──►                    ─►      (backflow if the parabola
y=0 ═════════             y=0 ═════════        y=0 ═══◄══       is pushed the other way)
  Couette: line            Poiseuille: parabola          sum (8.5)
""", words=r"Because the equation turns out linear, the two causes simply **add**.")
note("N03 [B] · Fully developed flow.", r"""
Near the channel inlet, boundary layers grow on both walls and the profile changes with x — inside this *entrance length*
$\partial u/\partial x\ne0$, so continuity $\partial u/\partial x+\partial v/\partial y=0$ forces $v\ne0$ and the flow is not
parallel. Downstream the layers merge and the profile stops changing: u = u(y) alone. The figure below draws the wall-layer
edge with the diffusion estimate $\delta_{99}\sim3.64\sqrt{\nu t}$ *(8.31)* and t = x/U — **our estimate** (C09 derives the
3.64), not the book's; Ch. 9 treats entry flow properly.
""")
remind([
    ("scipy.optimize.brentq", "`optimize.brentq(f, a, b)` finds the root of f between a and b where f changes sign (Ch. 3 P108); here: where the wall layer reaches mid-gap."),
])
nb.figure(r"""
h, Uin = 0.01, 0.01                                       # gap 1 cm, inflow speed 1 cm/s [m, m/s]
edge = lambda x: ch08.diffusion_thickness(x/Uin, NU_W)    # wall-layer edge δ(x) = 3.643√(ν x/U) [m] (our estimate)
x_e = optimize.brentq(lambda x: edge(x) - h/2, 1e-6, 1.0)  # where the two layers meet at mid-gap [m]
print(f"layers meet at x_e = {x_e:.4f} m (about {x_e/h:.1f} gap widths)")
x = np.linspace(1e-5, 0.05, 400)                          # distance from the inlet [m]
fig, ax = plt.subplots(figsize=(8, 3.2))
channel(ax, h, U=0.0, L=0.05)                             # two fixed walls (drawing only)
xl = x[x <= x_e]                                          # the edges exist only until they meet
ax.plot(xl, edge(xl), color=COLORS["accent"], ls="--", lw=1.5, label="wall-layer edge (our estimate)")   # lower wall's layer
ax.plot(xl, h - edge(xl), color=COLORS["accent"], ls="--", lw=1.5)                                     # the upper wall's layer
ax.axvline(x_e, color=COLORS["muted"], lw=0.8); ax.text(x_e + 5e-4, 0.0105, "x_e", fontsize=8)
yy = np.linspace(0, h, 41)                                # heights for the three profiles [m]
V = Uin                                                   # the mean speed is conserved along the channel [m/s]
for x0 in (0.003, 0.012):                                 # two developing stations: plug core, erf-shaped wall layers
    d = np.sqrt(NU_W*x0/Uin)                              # diffusion length √(ν x/U) there [m]
    shape = special.erf(yy/(2*d)) * special.erf((h - yy)/(2*d))   # schematic developing profile (0 at both walls)
    prof = V * shape / np.trapezoid(shape, yy) * h        # scaled to the same mean speed V
    profile_arrows(ax, yy, prof, x0=x0, scale=0.5, every=4, color=COLORS["blue"])
u_dev = ch08.channel_flow(yy, h, U=0.0, dpdx=-12*MU_W*V/h**2)   # fully developed: Poiseuille with mean V ((8.5), U = 0)
profile_arrows(ax, yy, u_dev, x0=0.035, scale=0.5, every=4, color=COLORS["orange"])
ax.set_xlim(0, 0.05); ax.set_xlabel("distance from the inlet x [m]"); ax.set_ylabel("y [m]")
ax.legend(loc="upper right", fontsize=8)
ax.set_title("Entrance region: wall layers grow like √x and meet; then u = u(y) only", fontsize=10)
savefig(fig, "ch08", "c02_entrance"); plt.show()
""", see=r"""Two dashed edges growing like √x from the walls and meeting at $x_e$ (printed above, about 1.9 cm — about twice the
gap); blue profiles still developing (flat core, thin wall layers), then the orange parabola of the fully developed flow.""",
    read=r"""Left of $x_e$ the core still moves as a plug; right of it the profile no longer changes with x — it is the parabola of
(8.5) with U = 0 — and every result of this section applies.""",
    change=r"""…the inflow were 10 times faster: the edge $\delta\propto\sqrt{\nu x/U}$ grows √10 times more slowly, so the entrance
length grows tenfold ($x_e\propto Uh^2/\nu$).""")
remind([
    ("partial derivative", "∂u/∂x is the rate of change of u with x while y, z, t are held fixed (Ch. 1 P25)."),
])
note("N04 [B] · Continuity kills v.", r"""
With $\partial u/\partial x=0$, continuity gives $\partial v/\partial y=0$; v is then the same at every height and equals
its wall value 0: $\mathbf u=(u(y),0,0)$. This is steps 1–3 of the derivation below.
""")
note("N05 [B], N06 [B] · The reduced equations.", r"""
What is left of the momentum equations is (8.4a) and (8.4b) below: no acceleration at all — the advective term vanishes
*exactly*, which is why the problem is linear.
""", equation=EQ["8.4a"] + r",\qquad " + EQ["8.4b"], ref="8.4a, 8.4b")
D("D02", ref="8.4a")
note("N07 [B] · Separation of variables in one line.", r"""
From $0=-\frac1\rho\frac{\partial p}{\partial y}$ *(8.4b)* p depends on x only; then the balance below says "a function of
x = a function of y", so both equal one constant: the pressure falls **linearly** along the channel. The same argument
returns in the pipe (C03), the rotating cylinders (C04) and the lubrication gap (C06).
""", equation=r"\frac{1}{\mu}\frac{dp}{dx}=\frac{d^2u}{dy^2}=\text{const}")
remind([
    ("separation of variables", "if f(x) = g(y) for all x and y, both sides must equal the same constant (Ch. 1 P42)."),
])
D("D03", ref="8.4a")
remind([
    ("integrating an ODE twice and fixing two constants", "a second-order linear ODE has a two-parameter family of solutions; two conditions pick one (Ch. 1 P44, linear second-order ODE)."),
    ("fundamental theorem of calculus", "integrating a derivative gives back the function plus a constant: $\\int u''\\,dy=u'+c_1$ (Ch. 2 P84)."),
    ("sympy and sympy.dsolve", "`sp.dsolve(sp.Eq(lhs, rhs))` solves an ODE symbolically; `C1`, `C2` in its answer are the free constants (Ch. 1 P40)."),
    ("two linear equations in two unknowns, np.linalg.solve", "`np.linalg.solve(M, b)` returns x with M x = b (Ch. 2 P57); we use it for two wall conditions."),
])
gloss("Integrating an ODE twice", r"""
Integrating $u''=c$ twice brings two constants, $u=cy^2/2+c_1y+c_2$; two wall conditions fix them (a 2 × 2 linear system).
sympy's `dsolve` does it in one line:""")
nb.code(r"""
y_, c_ = sp.symbols("y c")                               # height and the constant curvature (symbols)
u_ = sp.Function("u")                                    # the unknown profile u(y)
print(sp.dsolve(sp.Eq(u_(y_).diff(y_, 2), c_)))          # u'' = c  →  u = C1 + C2·y + c·y²/2
""")
note("N08 [B] · The book's constants.", r"""
The book integrates to the line below with B = 0 and $A=\frac h2\frac{dp}{dx}-\frac{\mu U}{h}$; its A and B are minus our
$c_1$ and $c_2$ (D04 step 3).
""", equation=r"0=-\frac{y^2}{2}\frac{dp}{dx}+\mu u+Ay+B")
D("D04", ref="8.5")
nb.worked_example("a 1 cm channel, belt at 10 cm/s, adverse gradient 4 Pa/m", r"""
Take h = 0.01 m, U = 0.1 m/s, μ = 10⁻³ Pa s, dp/dx = +4 Pa/m (pressure rising downstream), at the quarter height
y = h/4 = 0.0025 m.
1. Couette part: $Uy/h=0.1\times0.25=0.025$ m/s.
2. Poiseuille part: $-\frac{1}{2\mu}\frac{dp}{dx}y(h-y)=-\frac{4}{2\times10^{-3}}\times0.0025\times0.0075=-0.0375$ m/s.
3. Sum: u = −0.0125 m/s — the fluid moves **backwards** at the quarter height.
4. Threshold: $2\mu U/h^2=2\times10^{-3}\times0.1/10^{-4}=2$ Pa/m; we are at twice it.
5. Flow rate $Q=\frac{Uh}2-\frac{h^3}{12\mu}\frac{dp}{dx}=5\times10^{-4}-3.33\times10^{-4}=1.67\times10^{-4}$ m²/s: still
   forward overall.
""")
nb.code(r"""
h, U = 0.01, 0.1                                         # gap [m] and wall speed [m/s] of the tiny example
y = np.linspace(0, h, 9)                                 # nine heights, y = 0, h/8, …, h [m]
u = ch08.channel_flow(y, h, U=U, dpdx=4.0)               # (8.5) with an adverse gradient dp/dx = +4 Pa/m
print("u(y) =", np.round(u, 4), "m/s   (u at y = h/4 is the third entry)")
Q, V = ch08.channel_flow_rate(h, U=U, dpdx=4.0)          # flow rate per unit width [m²/s] and mean speed Q/h [m/s]
print(f"Q = {Q:.4g} m²/s, V = {V:.4g} m/s, threshold 2μU/h² = {ch08.channel_backflow_threshold(U, h):.3g} Pa/m")
s = ch08.couette_poiseuille_state(h, U, 4.0)             # everything the explainer shows, in one dictionary
print({k: (round(s[k], 6) if not isinstance(s[k], bool) else s[k]) for k in
       ("Q_couette", "Q_poiseuille", "tau_bottom", "tau_top", "backflow", "y_reversal", "zero_flow_dpdx")})
ps = ch08.parallel_flow_sympy("channel")                 # D02–D04 done symbolically
print("sympy profile:", ps["profile"], "   (G here is the symbol for dp/dx)")
print("reduced x-momentum (8.4a):", ps["momentum_x"])   # with p = p(x): the y-equation (8.4b) is satisfied identically
""", explain=r"""
1. `channel_flow` evaluates (8.5); at y = h/4 it gives −0.0125 m/s, the tiny example's number.
2. `channel_flow_rate` integrates (8.5) exactly: Q = 1.667 × 10⁻⁴ m²/s and V = Q/h.
3. `couette_poiseuille_state` bundles the flow rate split into its two parts ($Uh/2$ and $-h^3(dp/dx)/12\mu$), the signed
   wall stresses (−0.01 Pa at the floor, +0.03 Pa at the top), the verdict `backflow = True`, where u changes sign —
   y = 5 mm, the reversed layer fills the lower half — and the gradient that stops the net flow, $6\mu U/h^2$ = 6 Pa/m.
4. `parallel_flow_sympy("channel")` repeats D02–D04 symbolically: its profile, written with the symbol G for dp/dx, is
   $\frac{Uy}{h}+\frac{G}{2\mu}(y^2-hy)$ — exactly (8.5).
""")
remind([
    ("finite differences (tridiagonal solve from scratch)", "$u''\\approx(u_{i+1}-2u_i+u_{i-1})/\\Delta y^2$ turns the ODE into one linear equation per grid point (Ch. 1 P21)."),
])
gloss("Two numpy helpers", r"""`np.diag(v, k)` builds a matrix with the vector v on its k-th diagonal (k = 0 main, ±1 just above or
below); `np.r_[a, b, c]` glues numbers and arrays end to end into one array.""")
nb.check_agree(r"""
N = 101; yy = np.linspace(0, h, N); dy = yy[1] - yy[0]   # 101 grid points across the gap [m]
Amat = (np.diag(-2*np.ones(N-2)) + np.diag(np.ones(N-3), 1) + np.diag(np.ones(N-3), -1)) / dy**2   # u'' on interior points
rhs = np.full(N-2, 4.0/MU_W)                             # μ u'' = dp/dx  →  u'' = (dp/dx)/μ at every interior point
rhs[-1] -= U/dy**2                                       # the known top-wall value u(h) = U moves to the right side
u_in = np.linalg.solve(Amat, rhs)                        # interior velocities [m/s]
u_fd = np.r_[0.0, u_in, U]                               # add the wall values u(0) = 0, u(h) = U
print("max |FD − (8.5)| =", np.max(np.abs(u_fd - ch08.channel_flow(yy, h, U=U, dpdx=4.0))), "m/s")
assert np.allclose(u_fd, ch08.channel_flow(yy, h, U=U, dpdx=4.0), rtol=1e-10, atol=1e-12)   # same profile
""")
nb.md("Central differences are exact for a parabola, so the hand-built solver reproduces $u(y)=\\frac Uh\\,y-\\frac1{2\\mu}"
      "\\frac{dp}{dx}\\,y(h-y)$ *(8.5)* to round-off.")
note("N09 [B] · Favourable, adverse, backflow (Fig. 8.4).", r"""
$dp/dx<0$ pushes with the wall: a fuller profile. $dp/dx>0$ pushes against it; the fluid next to the fixed wall reverses as
soon as the slope du/dy at y = 0 turns negative, i.e. when the inequality below holds (ours — the book only draws it).
Ch. 9 meets the same fight between wall shear and an adverse gradient: that is how boundary layers separate.
""", equation=r"\frac{dp}{dx}>\frac{2\mu U}{h^{2}}")
remind([
    ("shear stress τ = μ du/dy, signed", "τ is the force per area that the fluid above a level exerts on the fluid below, in the x-direction (Ch. 1 P05); here we keep its sign."),
    ("derivative as a slope", "du/dy at the floor is the slope of the profile there; its sign says which way the fluid next to the wall moves (Ch. 1 P19)."),
    ("inequalities under division by a positive number", "multiplying or dividing both sides by a positive number keeps the direction of an inequality; a negative one flips it (Ch. 1 P48)."),
    ("first-order Taylor expansion", "near a point, $u(y)\\approx u(0)+u'(0)\\,y$ (Ch. 1 P26)."),
])
D("D05", ref="8.5")
remind([
    ("definite integral, scipy.integrate.quad", "`integrate.quad(f, a, b)` returns $\\int_a^bf\\,dy$ and an error estimate (Ch. 3 P87)."),
])
note("N10 [B] · Flow rate and mean velocity.", r"""
Two polynomial integrals of (8.5) give Q and V below — a favourable gradient raises Q, an adverse one lowers it, and Q = 0 at
$dp/dx=6\mu U/h^2$ (6 Pa/m for the tiny example). ⚠️ The book writes the middle expression of V as $V=\int_0^hu\,dy$ without the
1/h; the units give it away (m²/s is not m/s).
""", equation=r"Q=\int_0^hu\,dy=\frac{Uh}{2}\Big[1-\frac{h^2}{6\mu U}\frac{dp}{dx}\Big],\qquad V\equiv\frac Qh=\frac U2\Big[1-\frac{h^2}{6\mu U}\frac{dp}{dx}\Big]")
nb.code(r"""
Q_num = integrate.quad(lambda yy: ch08.channel_flow(yy, h, U=U, dpdx=4.0), 0, h)[0]   # ∫₀ʰ u dy numerically [m²/s]
print(f"quad: {Q_num:.6g} m²/s   formula: {ch08.channel_flow_rate(h, U=U, dpdx=4.0)[0]:.6g} m²/s")   # the same
""")
nb.figure(r"""
h, U = 0.01, 0.1                                          # the tiny example's channel [m, m/s]
yy = np.linspace(0, h, 201)                               # heights [m]
cases = [-4.0, 0.0, 2.0, 4.0, 8.0]                        # dp/dx [Pa/m]: favourable, none, threshold, 2×, 4×
cols = [COLORS["blue"], COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["amber"]]   # blue → amber
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6))
for G_, c in zip(cases, cols):
    u = ch08.channel_flow(yy, h, U=U, dpdx=G_)            # (8.5)
    Q_ = ch08.channel_flow_rate(h, U=U, dpdx=G_)[0]       # net flow rate [m²/s]: its sign matters beyond 6μU/h²
    tag = ", net flow reversed" if Q_ < 0 else ""
    a1.plot(u/U, yy/h, color=c, lw=2, label=f"dp/dx = {G_:+g} Pa/m (Q = {1e4*Q_:+.2f} cm²/s{tag})")
    if G_ > ch08.channel_backflow_threshold(U, h):       # shade the reversed layer (backflow) amber
        a1.fill_betweenx(yy/h, 0, u/U, where=u < 0, color=COLORS["amber"], alpha=0.35)
    a2.plot(1e3*ch08.channel_shear_stress(yy, h, U=U, dpdx=G_), yy/h, color=c, lw=2)   # signed τ(y) [mPa]
a1.plot(ch08.channel_flow(yy, h, U=0.0, dpdx=-4.0)/U, yy/h, color=COLORS["orange"], ls="--", lw=1.2,
        label="pure Poiseuille (U = 0, dp/dx = −4)")      # the parabola alone, for comparison
a1.axvline(0, color=COLORS["muted"], lw=0.8)
a1.set_xlabel("u/U [–]"); a1.set_ylabel("y/h [–]"); a1.legend(fontsize=7, loc="upper left")
a1.set_title("Profiles: the floor flow reverses past 2μU/h²", fontsize=10)
a2.axvline(0, color=COLORS["muted"], lw=0.8)
a2.set_xlabel("shear stress τ [mPa]"); a2.set_ylabel("y/h [–]"); a2.set_title("Stress: straight lines; τ(0) = 0 at the threshold", fontsize=10)
fig.suptitle("Adverse pressure beats the wall near the floor once dp/dx > 2μU/h²", fontsize=11)
savefig(fig, "ch08", "c02_couette_poiseuille"); plt.show()
""", see=r"""Five profiles through (0, 0) and (1, 1); the two most adverse ones dip below zero near the floor (amber). At 4 Pa/m
the net flow is still forward (Q > 0 in the legend); at 8 Pa/m, beyond the zero-flow gradient $6\mu U/h^2$ (6 Pa/m here), the
net flow itself has reversed (Q < 0). The stresses are straight lines, and the dashed orange parabola is pure Poiseuille flow —
our remake of Fig. 8.4.""",
    read=r"""The slope of a profile at y = 0 is proportional to the floor stress; it crosses zero exactly at the 2 Pa/m curve (vertical
at the floor), matching panel (b), where τ(0) = 0 for that case.""",
    change=r"""…the gap were halved to 5 mm: the threshold $2\mu U/h^2$ quadruples to 8 Pa/m, so the 4 Pa/m case no longer reverses.""")
remind([
    ("slider_figure", "`slider_figure(fn, name, values, …)` precomputes one set of curves per slider position, so the plotly figure works on the web page without Python (Ch. 1 P17)."),
])
nb.plotly(r"""
h, U = 0.01, 0.1                                          # channel [m], wall speed [m/s]
yy = np.linspace(0, h, 81)                                # heights [m]


def frame(G_):                                            # the curves for one slider value dp/dx = G_ [Pa/m]
    u = ch08.channel_flow(yy, h, U=U, dpdx=G_)            # (8.5)
    return {"sum u(y)": (u, yy/h),                                                      # the profile
            "Couette part Uy/h": (U*yy/h, yy/h),                                        # the line
            "Poiseuille part": (ch08.channel_flow(yy, h, U=0.0, dpdx=G_), yy/h),        # the parabola
            "reversed layer (u < 0)": (np.where(u < 0, u, np.nan), yy/h)}               # backflow, if any


Gs = np.linspace(-8, 12, 21 if not FAST else 11)          # slider values dp/dx [Pa/m]
fig = slider_figure(frame, "dp/dx", Gs, unit="Pa/m", xlabel="u [m/s]", ylabel="y/h [–]", title="")
titles = []                                               # one title per slider value, with the verdict
for G_ in Gs:
    s = ch08.couette_poiseuille_state(h, U, G_)           # Q, backflow flag, threshold
    verdict = f"backflow below y = {1e3*s['y_reversal']:.1f} mm" if s["backflow"] else "all forward"
    titles.append(f"A line plus a parabola: Q = {1e4*s['Q']:.2f} cm²/s, {verdict}")
step_titles(fig, titles)
colors = {"sum u(y)": COLORS["blue"], "Couette part Uy/h": COLORS["rose"], "Poiseuille part": COLORS["orange"],
          "reversed layer (u < 0)": COLORS["amber"]}     # colours by meaning
recolor(fig, colors, {"Couette part Uy/h": "dash", "Poiseuille part": "dash"})
fig.show()
""", explain=r"""
The slider re-evaluates (8.5) for 21 pressure gradients (11 in FAST mode). Watch the blue sum cross u = 0 at the floor as
dp/dx passes 2 Pa/m: the amber reversed layer appears and the title switches from "all forward" to "backflow".
""")
remind([
    ("show_viz and the explainer tabs", "`show_viz(\"ch08\", slug)` embeds an interactive explainer; its tabs are Walkthrough, Explore, Explain, Derivation, Equations, Code and Check (Ch. 1 P18)."),
])
explainer("couette_poiseuille_backflow", "When does fluid flow backwards in a channel?", r"""
Drag the pressure gradient through zero and past $2\mu U/h^2$: the straight Couette part and the parabolic Poiseuille part add
up in front of you, the reversed layer opens at the floor and the flow rate splits into its two contributions — a family of
curves that a static figure can only sample.""",
          ["Start at pure Couette (dp/dx = 0) and slowly raise dp/dx — watch the floor-slope readout reach 0 at 2 Pa/m.",
           "Press the 'zero net flow' preset: the Couette and Poiseuille bars cancel exactly at $dp/dx=6\\mu U/h^2$.",
           "Click at mid-height in the profile view and read the arithmetic of u from its two parts.",
           "Halve the gap h and predict the new threshold before you look."])
note("N11 [C] · Linear stress everywhere.", r"""
A constant pressure gradient and a *linear* stress τ(y) are general for any fully developed channel flow — they survive in
turbulent flow for the time-averaged stress (Ch. 12 uses exactly this).
""")
whatif(r"""
…the channel were a round pipe? The same three moves (continuity, "function of x = function of R", integrate twice) work,
but the Laplacian brings a 1/R — and a new rule: the solution must stay finite on the axis (C03).
""")

nb.recap("R09", "Shear stress in cylindrical coordinates", r"""
In cylindrical coordinates (R, φ, z) the shear stress on a surface of constant R in the z-direction is
$\tau_{zR}=\mu\big(\frac{\partial u_R}{\partial z}+\frac{\partial u_z}{\partial R}\big)$ (Ch. 4's `strain_rate(..., "cylindrical")`
gives the same rate of strain).
""", where="Ch. 4 §4.5")
core("C03", "Poiseuille pipe flow and Hagen–Poiseuille", r"""
How much flows through a pipe for a given pressure drop — and why does halving the radius cut the flow sixteen-fold?
""", eqs=("8.6",))
problem(r"""
Blood squeezing through a capillary, water in a garden hose, oil in a pipeline, air in the lungs' smallest airways: a
pressure difference pushes fluid through a tube and viscosity at the wall holds it back. Doctors care because a slightly
narrowed artery needs a much larger pressure; engineers because pumping costs follow the same law.
""")
idea(r"""
R=a  ───────────────  u=0,  |τ| max = τ₀
       ─────────►
         ──────────────►    u_max = 2V on the axis
       ─────────►
R=a  ───────────────
""", words=r"A paraboloid of velocity, fastest on the axis, zero at the wall; the stress grows linearly from the axis to the wall.")
remind([
    ("cylindrical coordinates (R, φ, z) and unit vectors", "R is the distance from the axis, φ the angle round it, z along it; the unit vectors $\\mathbf e_R$, $\\mathbf e_\\varphi$ turn with φ (Ch. 3 P88)."),
])
note("N12 [B] · Set-up in (R, φ, z).", r"""
Take $\mathbf u=(0,0,u_z(R))$: continuity is automatic, the R- and φ-equations give $0=\partial p/\partial\varphi$ and
$0=\partial p/\partial R$, so p = p(z). The z-equation is the balance below.
""", equation=r"0=-\frac{dp}{dz}+\frac\mu R\frac{d}{dR}\Big(R\frac{du_z}{dR}\Big)")
P("P186", "Laplacian in cylindrical coordinates", r"""
For a function of R only, $\nabla^2u=\frac1R\frac{d}{dR}\big(R\frac{du}{dR}\big)$ — the R inside the bracket comes from the
circle's circumference growing with R. For an azimuthal velocity $u_\varphi(R)$ the φ-component of the *vector* Laplacian is
$\frac{d}{dR}\big[\frac1R\frac{d(Ru_\varphi)}{dR}\big]$ (it contains an extra $-u_\varphi/R^2$ because the unit vector
$\mathbf e_\varphi$ turns). Both come from `core.curvilinear`.
""", code=r"""
import sympy as sp
R = sp.symbols('R', positive=True)
u = R**2                                                 # a test profile
print(sp.simplify(sp.diff(R*sp.diff(u, R), R)/R))        # (1/R)(R u')' = 4 for u = R²
""")
remind([
    ("natural logarithm and ln R → −∞", "ln R is the antiderivative of 1/R, and ln R → −∞ as R → 0 (Ch. 1 P36)."),
    ("area element 2πR dR of a disc", "a thin ring of radius R and width dR has area 2πR dR, so a flow rate is $\\int u\\,2\\pi R\\,dR$ (Ch. 2 P83)."),
])
gloss("Regularity at the axis", r"""A solution containing A ln R is infinite on the axis; a physical velocity is finite there, so
A = 0.""")
D("D06", ref="8.6")
note("N13 [B] · Two constants, two conditions.", r"""
The general solution $u_z(R)=\frac{R^2}{4\mu}\frac{dp}{dz}+A\ln R+B$ keeps two constants; the axis removes A, the wall fixes
$B=-\frac{a^2}{4\mu}\frac{dp}{dz}$ (steps 7–9 above).
""")
note("N14 [B], N15 [B], N16 [B] · Stress, wall stress, Hagen–Poiseuille.", r"""
Differentiating (8.6) gives the linear stress (8.7) — zero on the axis, largest at the wall. The **wall stress** (8.8) is
negative for forward flow (dp/dz < 0): the wall pulls the fluid back. A force balance on a slug of fluid gives the same line
without the profile, which is why (8.8) also holds for averaged turbulent pipe flow (D07 steps 4–5; Ch. 12 builds the
friction velocity on it). **Hagen–Poiseuille:** $Q=\int_0^au\,2\pi R\,dR=-\frac{\pi a^4}{8\mu}\frac{dp}{dz}$,
$V=\frac{Q}{\pi a^2}=-\frac{a^2}{8\mu}\frac{dp}{dz}$, $u_{max}=2V$, and the Darcy friction factor f = 64/Re (ours).
""", equation=EQ["8.7"] + r"\ \ \text{(8.7)},\qquad " + EQ["8.8"], ref="8.8")
D("D07", ref="8.8")
nb.worked_example("a 2 mm tube at dp/dz = −1000 Pa/m", r"""
a = 1 mm, μ = 10⁻³ Pa s, water.
1. $V=-\frac{a^2}{8\mu}\frac{dp}{dz}=\frac{10^{-6}}{8\times10^{-3}}\times1000=0.125$ m/s.
2. $u_{max}=2V=0.25$ m/s.
3. $Q=\pi a^2V=\pi\times10^{-6}\times0.125=3.93\times10^{-7}$ m³/s (0.39 mL/s).
4. $\tau_0=\frac a2\frac{dp}{dz}=-0.5$ Pa.
5. $\mathrm{Re}=V\cdot2a/\nu=0.125\times0.002/10^{-6}=250$: laminar ✓.
6. f = 64/250 = 0.256, and the definition $8\lvert\tau_0\rvert/(\rho V^2)=8\times0.5/(1000\times0.0156)=0.256$ ✓.
7. Halve a: V falls 4×, the area 4×, so Q falls 16×.
""")
nb.code(r"""
a, G = 1e-3, -1000.0                                     # radius [m], dp/dz [Pa/m] (favourable: pressure falls along z)
Q, V, umax = ch08.pipe_flow_rate(a, G)                   # Hagen–Poiseuille: Q [m³/s], mean V and axis speed [m/s]
print(f"Q = {Q:.4g} m³/s, V = {V:.4g} m/s, u_max = {umax:.4g} m/s")
print("τ₀ =", ch08.pipe_wall_stress(a, G), "Pa;  τ(R = a/2) =", ch08.pipe_shear_stress(0.5e-3, G), "Pa")   # (8.8), (8.7)
Re = V*2*a/NU_W                                          # pipe Reynolds number (diameter, mean speed)
print(f"Re = {Re:.0f}, f = {ch08.pipe_friction_factor(Re):.4g}")   # f = 64/Re
R = np.linspace(0, a, 5)                                 # five radii from the axis to the wall [m]
print("(8.6):         ", ch08.pipe_poiseuille(R, a, G))  # u_z(R) [m/s]
print("Ch. 3 far down:", ch03.pipe_profile(R, z=1e3, U_mean=V, R=a))   # Ch. 3's developing profile, 1 km downstream
print("Ch. 4 on axis: ", NS.exact_solution("pipe_poiseuille", np.array([0.0, 0.0, 0.0]), G=-G, R=a)[0][2])   # Ch. 4 (G = −dp/dz)
pp = ch08.parallel_flow_sympy("pipe")                    # D06 symbolically
print("sympy:", pp["profile"], "  unbounded term:", pp["unbounded_term"], "  (G here is the symbol for dp/dz)")
""", explain=r"""
1. `pipe_flow_rate` returns Q, V and u_max = 2V: 3.93 × 10⁻⁷ m³/s, 0.125 m/s, 0.25 m/s (the tiny example).
2. `pipe_wall_stress` and `pipe_shear_stress` are (8.8) and (8.7): −0.5 Pa at the wall, half of it at mid-radius.
3. Three independent earlier implementations agree with (8.6): this chapter's formula, Ch. 3's developing profile evaluated
   far downstream, and Ch. 4's exact-solution preset on the axis (0.25 m/s).
4. `parallel_flow_sympy("pipe")` returns the profile and the term $C_1\ln R$ it discarded because it is unbounded on the axis.
""")
remind([
    ("trapezoid rule (from scratch)", "`np.trapezoid(f, x)` adds up trapezoids under the sampled curve (Ch. 1 P37)."),
])
nb.check_agree(r"""
RR = np.linspace(0, a, 2001)                             # 2001 radii [m]
Q_trap = np.trapezoid(ch08.pipe_poiseuille(RR, a, G) * 2*np.pi*RR, RR)   # Q = ∫ u 2πR dR by trapezoids [m³/s]
print(f"trapezoid Q = {Q_trap:.6g}, formula Q = {Q:.6g} m³/s")
assert np.allclose(Q_trap, Q, rtol=1e-6)                 # same flow rate
F_pressure = -np.pi*a**2*G*1.0                           # pressure force on a 1 m slug: πa²[p(0) − p(L)] [N]
F_wall = -2*np.pi*a*1.0*ch08.pipe_wall_stress(a, G)      # wall friction force on its side: −2πaL τ₀ [N]
assert np.isclose(F_pressure, F_wall)                    # D07 step 4: the pressure push equals the wall drag
print(f"force balance on a 1 m slug: pressure {F_pressure:.4g} N = wall friction {F_wall:.4g} N")
""")
nb.figure(r"""
a, G = 1e-3, -1000.0                                      # the tiny example's tube
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6))
R = np.linspace(-a, a, 201)                               # a diameter, from wall to wall [m]
u = ch08.pipe_poiseuille(np.abs(R), a, G)                 # (8.6) along the diameter [m/s]
a1.plot(1e3*R, u, color=COLORS["blue"], lw=2, label="u_z(R), (8.6)")
for Rk in np.linspace(-0.8*a, 0.8*a, 9):                  # nine velocity arrows (blue)
    a1.annotate("", xy=(1e3*Rk, ch08.pipe_poiseuille(abs(Rk), a, G)), xytext=(1e3*Rk, 0),
                arrowprops=dict(arrowstyle="->", color=COLORS["blue"], lw=0.8))
a1.axvline(-1, color=COLORS["muted"], lw=3); a1.axvline(1, color=COLORS["muted"], lw=3)   # the pipe walls
a1.set_xlabel("R [mm] (across a diameter)"); a1.set_ylabel("u_z [m/s]", color=COLORS["blue"])
b1 = a1.twinx()                                           # second axis for the stress (rose)
b1.plot(1e3*R, np.sign(R)*ch08.pipe_shear_stress(np.abs(R), G), color=COLORS["rose"], lw=1.5, ls="--")   # τ = (R/2) dp/dz
b1.set_ylabel("τ [Pa] (on the +R side: (R/2)dp/dz)", color=COLORS["rose"]); b1.axhline(0, color=COLORS["grid"], lw=0.8)
b1.text(0.55, ch08.pipe_wall_stress(a, G), "τ₀ = −0.5 Pa", color=COLORS["rose"], fontsize=8, va="bottom")
a1.set_title("Paraboloid and linear stress", fontsize=10)
aa = np.geomspace(1e-4, 1e-2, 60)                         # radii from 0.1 mm to 1 cm [m]
QV = np.array([ch08.pipe_flow_rate(a_, G)[:2] for a_ in aa])   # (Q, V) for each radius (the function takes one a)
Qs, Res = QV[:, 0], QV[:, 1]*2*aa/NU_W                    # Q(a) [m³/s] and Re(a) = V·2a/ν with the mean speed
lam = Res < 2000                                          # where the laminar law is trustworthy
a2.loglog(aa[lam], Qs[lam], color=COLORS["blue"], lw=2, label="Hagen–Poiseuille (laminar)")
a2.loglog(aa[~lam], Qs[~lam], color=COLORS["blue"], lw=1, ls="--", label="Re > 2000: turbulent beyond")
a2.loglog(aa, Qs[0]*(aa/aa[0])**4, color=COLORS["muted"], lw=1, ls=":", label="slope 4")
a_cross = aa[np.argmax(~lam)]                             # first radius where Re passes 2000
a2.axvline(a_cross, color=COLORS["amber"], lw=0.8); a2.text(a_cross*1.1, Qs[0]*10, f"Re = 2000\na ≈ {1e3*a_cross:.1f} mm", fontsize=8, color=COLORS["amber"])
a2.set_xlabel("radius a [m]"); a2.set_ylabel("Q [m³/s] at dp/dz = −1000 Pa/m"); a2.legend(fontsize=7, loc="lower right")
a2.set_title("Q ∝ a⁴", fontsize=10)
fig.suptitle("Flow rate grows like the fourth power of the radius", fontsize=11)
savefig(fig, "ch08", "c03_pipe"); plt.show()
""", see=r"""(a) a parabola across the diameter with the stress line (rose dashed) growing from zero on the axis to −0.5 Pa at the
wall — our remake of the book's pipe sketch `N104`; (b) a straight log–log line of slope 4.""",
    read=r"""Each decade of radius gives four decades of flow; past Re = 2000 (the amber line, a ≈ 2.5 mm here) the laminar line
is no longer trustworthy, so it is dashed.""",
    change=r"""…the artery narrowed by 20 %: $0.8^4=0.41$ — the same pressure pushes 59 % less blood, so the heart must raise the
pressure 2.4× to keep Q.""")
nb.plotly(r"""
G = -1000.0                                               # dp/dz [Pa/m]
s_ = np.linspace(0, 1, 60)                                # R/a from the axis to the wall [–]


def frame(a_):                                            # curves for one radius a_ [m]
    return {"u_z(R), (8.6)": (s_, ch08.pipe_poiseuille(s_*a_, a_, G))}   # the paraboloid across the radius


radii = np.geomspace(2e-4, 2e-3, 13 if not FAST else 7)   # slider values a [m]
fig = slider_figure(frame, "a", radii, unit="m", xlabel="R/a [–]", ylabel="u_z [m/s]", title="")
titles = []                                               # Q and Re for each radius, in the title
for a_ in radii:
    Q_, V_, _ = ch08.pipe_flow_rate(a_, G)                # Hagen–Poiseuille flow rate and mean speed
    titles.append(f"a = {1e3*a_:.2f} mm: Q = {1e6*Q_:.3g} mL/s, Re = {V_*2*a_/NU_W:.0f} — Q ∝ a⁴")
step_titles(fig, titles)
recolor(fig, {"u_z(R), (8.6)": COLORS["blue"]})           # velocity profiles in blue
fig.show()
""", explain=r"""
Each slider step is one radius between 0.2 and 2 mm at the same dp/dz; the title prints Q and Re. Doubling a multiplies the
centre speed by 4 and Q by 16 (compare the steps near a = 0.5 mm and 1 mm); past about 2.5 mm Re would exceed 2000.
""")
whatif(r"""
…the walls moved instead of a pressure pushing — sideways, round the axis? Then the velocity is azimuthal, $u_\varphi(R)$,
and the same reduction gives circular Couette flow (C04).
""")

nb.recap("R10", "A viscous flow with no net viscous force", r"""
Outside a spinning cylinder the ideal-vortex velocity $u_\varphi=\Gamma/2\pi R$ has shear stress
$\sigma_{R\varphi}=\mu R\frac{d}{dR}\big(\frac{u_\varphi}R\big)=-\frac{2\mu\Omega_1R_1^2}{R^2}$ but no net viscous force on any
fluid element, because $\mu\nabla^2\mathbf u=-\mu\nabla\times\boldsymbol\omega=0$ when $\boldsymbol\omega=0$ (Ch. 5).
""", where="Ch. 5 §5.1")
nb.recap("R11", "Power in = dissipation", r"""
The torque that keeps the cylinder spinning does work at the rate $-2\pi R_1\sigma_{R\varphi}u_\varphi$ per unit length
(positive: $\sigma_{R\varphi}<0$), and all of it is dissipated in the fluid: $4\pi\mu\Omega_1^2R_1^2$ (Ch. 5's
`dissipation_outside_cylinder`; Exercise 8.12). The book writes the power as $(2\pi R_1)\tau_{R\varphi}u_\varphi$, which is
negative with its own sign of σ_Rφ — the power delivered *to* the fluid carries a minus sign.
""", where="Ch. 5 §5.1")
nb.recap("R12", "Solid-body rotation", r"""
A tank of fluid spun long enough rotates like a solid: $u_\theta=\omega r/2$ *(5.1)* with the vorticity ω = 2Ω everywhere,
i.e. $u_\varphi=\Omega R$ — no shear at all.
""", where="Ch. 5 §5.1")
remind([
    ("curl of a curl identity (∇²u = −∇×ω)", "$\\nabla\\times(\\nabla\\times\\mathbf u)=\\nabla(\\nabla\\cdot\\mathbf u)-\\nabla^2\\mathbf u$, so for ∇·u = 0 the viscous term is $\\mu\\nabla^2\\mathbf u=-\\mu\\nabla\\times\\boldsymbol\\omega$ (Ch. 4 P122)."),
])
core("C04", "Circular Couette flow and its two limits", r"""
What swirl is set up between two rotating cylinders, and what happens when one of them is taken away?
""", eqs=("8.10",))
problem(r"""
Stir tea between a spoon handle and the cup wall; spin the inner cylinder of a lab rig (a Taylor–Couette cell, the classic
apparatus of Ch. 11) or a viscometer. The fluid between two concentric rotating walls settles into a steady swirl. We want
that swirl, the pressure that holds it in its circle, and what it becomes when the outer wall is at infinity or the inner
one vanishes — the two vortices of Ch. 5.
""")
idea("", table=r"""
| piece | u_φ | vorticity | shear | appears alone when |
|---|---|---|---|---|
| solid-body rotation | AR | 2A | none | no inner cylinder (R₁ → 0) |
| free (line) vortex | B/R | 0 | yes, but no net force | outer wall at infinity (R₂ → ∞, Ω₂ = 0) |
""", words=r"Only two swirls satisfy the viscous balance, and the walls choose the mix.")
remind([
    ("centripetal acceleration in cylindrical coordinates", "a particle going round a circle of radius R at speed u_φ accelerates toward the axis at $u_\\varphi^2/R$ even at constant speed, because $\\mathbf e_\\varphi$ turns (Ch. 3 P105)."),
])
note("N17 [B] · Two equations left.", r"""
With $\mathbf u=(0,u_\varphi(R),0)$: the R-equation (first below) says pressure rises outward to supply the centripetal
acceleration; the φ-equation (second) is a pure viscous balance (no pressure, no advection).
""", equation=r"-\frac{u_\varphi^2}{R}=-\frac1\rho\frac{dp}{dR},\qquad 0=\mu\frac{d}{dR}\Big[\frac1R\frac{d}{dR}(Ru_\varphi)\Big]")
P("P187", "Euler–Cauchy (equidimensional) ODE", r"""
An ODE in which every term has the form $R^k\,d^ku/dR^k$ (each derivative comes with the same power of R) is solved by trying
$u=R^\lambda$: every term becomes a constant times $R^\lambda$, so λ must solve a polynomial. For $R^2u''+Ru'-u=0$ the
polynomial is $\lambda^2-1=0$, so $u=AR+B/R$.
""", code=r"""
import sympy as sp
R, lam = sp.symbols('R lambda')
u = R**lam                                               # trial power
print(sp.factor(sp.simplify((R**2*sp.diff(u, R, 2) + R*sp.diff(u, R) - u)/u)))   # (lambda - 1)*(lambda + 1)
""")
D("D08", ref="8.10")
note("N18 [B], N19 [B] · The general swirl and its constants.", r"""
The general solution (8.9) is solid-body rotation plus a line vortex; the walls give
$A=\frac{\Omega_2R_2^2-\Omega_1R_1^2}{R_2^2-R_1^2}$ and $B=-\frac{(\Omega_2-\Omega_1)R_1^2R_2^2}{R_2^2-R_1^2}$ (D08 steps
8–10).
""", equation=EQ["8.9"], ref="8.9")
remind([
    ("limits and orders of smallness", "as R₂ → ∞ a ratio like R₁²/R₂² → 0; dividing top and bottom by the growing quantity first makes the limit visible (Ch. 2 P68)."),
])
D("D09", ref="8.11")
note("N21 [B] · A viscous flow that is irrotational.", r"""
With R₂ → ∞ and Ω₂ = 0, (8.11) is the ideal vortex $u_\theta=\frac{\Gamma}{2\pi r}$ *(5.2)* with $\Gamma=2\pi\Omega_1R_1^2$ —
shear stress everywhere, yet zero net viscous force (R10): the only completely irrotational viscous solution in the book. If
gravity acts along the axis, the free surface dips toward the cylinder (Fig. 8.7, remade below). ⚠️ Ch. 5's
`rotating_cylinder_flow(r, a, omega)` takes the cylinder's *vorticity* ω = 2Ω₁. The other limit is R12's solid-body
rotation, $u_\varphi=\Omega_2R$ *(8.12)*.
""")
nb.worked_example("a 1 cm cylinder spinning at 1 rad/s inside a fixed 2 cm cup", r"""
R₁ = 0.01 m, R₂ = 0.02 m, Ω₁ = 1 rad/s, Ω₂ = 0.
1. $A=\frac{0-1\times10^{-4}}{4\times10^{-4}-10^{-4}}=-\frac13$ s⁻¹.
2. $B=-\frac{(0-1)\times10^{-4}\times4\times10^{-4}}{3\times10^{-4}}=1.333\times10^{-4}$ m²/s.
3. At mid-gap R = 0.015 m: $u_\varphi=-\frac13\times0.015+\frac{1.333\times10^{-4}}{0.015}=-0.0050+0.0089=0.0039$ m/s.
4. Checks: at R₁, −0.00333 + 0.01333 = 0.0100 = Ω₁R₁ ✓; at R₂, −0.00667 + 0.00667 = 0 ✓.
""")
nb.code(r"""
R1, R2 = 0.01, 0.02                                      # cylinder radii [m]
u, A, B = ch08.circular_couette(0.015, R1, R2, 1.0, 0.0, return_coeffs=True)   # (8.10) at mid-gap, and A, B of (8.9)
print(f"u_φ(1.5 cm) = {u:.4g} m/s, A = {A:.4g} 1/s, B = {B:.4g} m²/s")
print("R₂ = ∞ (free vortex):", ch08.circular_couette(0.05, R1, np.inf, 1.0, 0.0),
      "  Ch. 5 with ω = 2Ω₁:", ch05.rotating_cylinder_flow(0.05, R1, 2*1.0)[0])   # (8.11) vs Ch. 5 at R = 5 cm
print("R₁ = 0 (solid body): ", ch08.circular_couette(0.01, 0.0, R2, 0.0, 3.0),
      "  Ch. 5 solid body:", VX.solid_body_rotation(0.01, 3.0))               # (8.12) vs Ch. 5 at R = 1 cm, Ω₂ = 3 rad/s
print(f"pressure rise across the gap: {ch08.circular_couette_pressure(R2, R1, R2, 1.0, 0.0):.4g} Pa")
pw = ch08.circular_couette_power(R1, R2, 1.0, 0.0, mu=MU_W)   # torques, power and dissipation per metre of length
print({k: f"{pw[k]:.4g}" for k in ("torque_inner", "torque_outer", "power_in", "power_out", "dissipation")})
print("sympy profile:", ch08.parallel_flow_sympy("circular_couette")["profile"])   # D08 symbolically
""", explain=r"""
1. `circular_couette` evaluates (8.10); with `return_coeffs=True` it also returns A = −1/3 s⁻¹ and B = 1.333 × 10⁻⁴ m²/s,
   the tiny example's numbers.
2. The exact limit branches (R₂ = ∞, R₁ = 0) agree with Ch. 5's vortex (called with ω = 2Ω₁) and solid-body functions.
3. The pressure rises by 0.0217 Pa from the inner to the outer wall (N20 below).
4. The torque is the same in size at both walls (opposite signs: the inner wall drives, the outer wall holds) — a steady
   angular-momentum balance; and every watt the inner cylinder puts in is dissipated (`power_in` = `dissipation`,
   1.676 × 10⁻⁶ W per metre), since the outer wall is at rest (`power_out` = 0).
5. The sympy profile is (8.10) over a common denominator.
""")
nb.check_agree(r"""
M = np.array([[R1, 1/R1], [R2, 1/R2]])                   # no slip: A R + B/R = Ω R at each wall
rhs = np.array([1.0*R1, 0.0*R2])                         # wall speeds Ω₁R₁, Ω₂R₂ [m/s]
A_, B_ = np.linalg.solve(M, rhs)                         # the 2 × 2 system of D08 step 8
print(A_, B_)
assert np.allclose([A_, B_], [A, B], rtol=1e-12)         # same constants as the library
""")
note("N20 [B] · The pressure that holds the swirl (ours; the book says only that it 'can be determined').", r"""
Integrating the R-balance $dp/dR=\rho u_\varphi^2/R$ with $u_\varphi=AR+B/R$ gives the pressure below. In the atmosphere this
centripetal balance is the *cyclostrophic* balance of tornadoes and dust devils (Ch. 13).
""", equation=r"p(R)=p_1+\rho\Big[\frac{A^2}{2}(R^2-R_1^2)+2AB\ln\frac{R}{R_1}-\frac{B^2}{2}\Big(\frac{1}{R^2}-\frac{1}{R_1^2}\Big)\Big]")
remind([
    ("np.gradient", "`np.gradient(f, x)` returns df/dx by central differences on the samples (Ch. 1 P22)."),
])
nb.code(r"""
Rg = np.linspace(R1, R2, 2001)                           # radii across the gap [m]
pg = ch08.circular_couette_pressure(Rg, R1, R2, 1.0, 0.0)   # N20's pressure [Pa]
ug = ch08.circular_couette(Rg, R1, R2, 1.0, 0.0)         # (8.10) [m/s]
dpdR, cent = np.gradient(pg, Rg)[5:-5], (1000*ug**2/Rg)[5:-5]   # finite-difference dp/dR and ρu²/R [Pa/m], away from the ends
print(f"max |dp/dR − ρu²/R| = {np.max(np.abs(dpdR - cent)):.1e} Pa/m (the largest value of ρu²/R is {cent.max():.2f} Pa/m)")
assert np.allclose(dpdR, cent, rtol=1e-5, atol=1e-5)     # they agree to the finite-difference error
""")
note("N22 [B] · Why these three flows are solvable.", r"""
All three are confined between walls, and symmetry kills the advective acceleration: $\mathbf u\cdot\nabla\mathbf u=0$ in the
channel and the pipe, while in circular Couette flow only the centripetal part $-(u_\varphi^2/R)\mathbf e_R$ survives,
balanced by the pressure. Other exact solutions exist (§8.4 and the exercises); Ch. 10 uses such exact solutions to test
flow solvers.
""")
nb.code(r"""
for case in ("channel", "pipe", "circular_couette"):     # u·∇u computed symbolically with core.curvilinear
    print(case, "→ (u·∇)u =", ch08.advective_acceleration_check(case)["advective"])   # v(R) is the symbol for u_φ(R)
""")
nb.figure(r"""
R1, R2 = 0.01, 0.02                                       # the annulus of the tiny example [m] (R₁/R₂ = 0.5)
Rr = np.linspace(R1, R2, 200)                             # radii across the gap [m]
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12, 3.8), gridspec_kw=dict(width_ratios=[1.2, 1, 1.1]))
cases = [((1, 0), "inner only", COLORS["blue"]), ((0, 1), "outer only", COLORS["teal"]),
         ((1, 1), "co-rotating (solid body)", COLORS["accent"]), ((1, -1), "counter-rotating", COLORS["amber"])]
for (O1, O2), lab, c in cases:                            # four (Ω₁, Ω₂) pairs [rad/s]
    a1.plot(1e2*Rr, 1e3*ch08.circular_couette(Rr, R1, R2, O1, O2), color=c, lw=2, label=lab)   # (8.10) [mm/s]
_, A, B = ch08.circular_couette(Rr, R1, R2, 1, 0, return_coeffs=True)   # the two parts of the inner-only case
a1.plot(1e2*Rr, 1e3*A*Rr, color=COLORS["orange"], ls="--", lw=1, label="AR part (inner only)")
a1.plot(1e2*Rr, 1e3*B/Rr, color=COLORS["rose"], ls="--", lw=1, label="B/R part (inner only)")
a1.axhline(0, color=COLORS["muted"], lw=0.8)
a1.set_xlabel("R [cm]"); a1.set_ylabel("u_φ [mm/s]"); a1.legend(fontsize=7); a1.set_title("Every swirl is AR + B/R", fontsize=10)
ph = np.linspace(0, 2*np.pi, 13)[:-1]                     # 12 angles
for Rk in np.linspace(R1*1.1, R2*0.95, 5):                # 5 radii: arrows ∝ u_φ, tangent to circles (inner-only case)
    uk = ch08.circular_couette(Rk, R1, R2, 1, 0)
    a2.quiver(Rk*np.cos(ph), Rk*np.sin(ph), -uk*np.sin(ph), uk*np.cos(ph), color=COLORS["blue"], scale=0.06, width=0.006)
t_ = np.linspace(0, 2*np.pi, 200)
a2.fill(R1*np.cos(t_), R1*np.sin(t_), color=COLORS["muted"]); a2.plot(R2*np.cos(t_), R2*np.sin(t_), color=COLORS["ink"], lw=2)
a2.set_aspect("equal"); a2.axis("off"); a2.set_title("Seen from above (inner cylinder turning)", fontsize=10)
Rf = np.linspace(R1, 6*R1, 200)                           # outside a lone cylinder (R₂ = ∞) [m]
uf = ch08.circular_couette(Rf, R1, np.inf, 1.0, 0.0)      # (8.11): Ω₁R₁²/R
Gam = 2*np.pi*1.0*R1**2                                   # circulation Γ = 2πΩ₁R₁² [m²/s]
zs = -Gam**2/(8*np.pi**2*9.81*Rf**2)                      # free-surface height relative to far away (ours, Ch. 5's funnel) [m]
a3.plot(1e2*Rf, 1e3*uf, color=COLORS["blue"], lw=2, label="u_φ = Ω₁R₁²/R (8.11) [mm/s]")
a3.plot(1e2*Rf, 1e6*zs, color=COLORS["accent"], lw=2, label="free-surface dip z_s [µm] (ours)")
a3.axhline(0, color=COLORS["muted"], lw=0.8)
a3.set_xlabel("R [cm]"); a3.legend(fontsize=7); a3.set_title("R₂ → ∞: the ideal vortex and its dip", fontsize=10)
fig.suptitle("Every swirl between cylinders is AR + B/R", fontsize=11)
savefig(fig, "ch08", "c04_circular_couette"); plt.show()
""", see=r"""(a) four curves joining the two wall speeds; the co-rotating one is a straight line; the counter-rotating one crosses
zero inside the gap; (b) the annulus from above with arrows growing toward the inner cylinder (our remake of Fig. 8.6
`N104`); (c) the free vortex outside a lone cylinder and the free surface dipping toward it (our Fig. 8.7).""",
    read=r"""A straight line is pure AR (B = 0, rigid rotation); a curve bending like 1/R carries a vortex part (dashed rose); where
$u_\varphi$ changes sign the fluid turns the other way. In (c) the dip is tiny here (micrometres) because Ω₁ is only 1 rad/s.""",
    change=r"""…the outer cylinder turned faster than the inner in the same sense: $(Ru_\varphi)^2$ grows outward and the flow is stable
(Rayleigh's criterion, Ch. 11); spin only the inner one and it becomes unstable above a critical speed — the Taylor vortices.""")
nb.plotly(r"""
R1, R2 = 0.01, 0.02                                       # annulus [m]
Rr = np.linspace(R1, R2, 80)                              # radii [m]


def frame(r_):                                            # r_ = Ω₂/Ω₁ with Ω₁ = 1 rad/s
    return {"u_φ(R), (8.10)": (1e2*Rr, 1e3*ch08.circular_couette(Rr, R1, R2, 1.0, r_)),   # the swirl [mm/s]
            "solid body Ω₂R": (1e2*Rr, 1e3*r_*Rr),                                          # rigid rotation reference
            "free vortex Ω₁R₁²/R": (1e2*Rr, 1e3*R1**2/Rr)}                                  # ideal-vortex reference


ratios = np.linspace(-2, 2, 21 if not FAST else 11)       # slider values Ω₂/Ω₁ [–]
fig = slider_figure(frame, "Ω₂/Ω₁", ratios, unit="", xlabel="R [cm]", ylabel="u_φ [mm/s]", title="")
titles = []                                               # A, B and the Rayleigh verdict for each ratio
for r_ in ratios:
    st = ch08.circular_couette_state(R1, R2, 1.0, r_)     # constants of (8.9) and (Ru_φ)² increasing outward?
    tag = "Rayleigh-stable" if st["rayleigh_stable"] else "can become unstable (Ch. 11)"
    titles.append(f"A = {st['A']:.2f} 1/s, B = {1e4*st['B']:.2f} cm²/s — {tag}")
step_titles(fig, titles)
recolor(fig, {"u_φ(R), (8.10)": COLORS["blue"], "solid body Ω₂R": COLORS["orange"], "free vortex Ω₁R₁²/R": COLORS["rose"]},
        {"solid body Ω₂R": "dash", "free vortex Ω₁R₁²/R": "dash"})
fig.show()
""", explain=r"""
The slider sets the outer wall's rotation relative to the inner one (Ω₁ = 1 rad/s). The title prints A, B and whether
Rayleigh's criterion (angular momentum growing outward) holds; the dashed references are the two pure limits.
""")
whatif(r"""
…the walls were not exactly parallel — a thin gap that narrows slowly, as under a bearing pad? The flow is then *nearly*
parallel, the nonlinear term is small instead of zero, and we need a way to decide what may be dropped: the lubrication
approximation (C05).
""")

# =====================================================================================================================
# A.3 §8.3 — R13, C05 (N23, N24, N26–N31, D10, D11, P188, E2), C06 (N25, N32–N34, N36, D12, D13, P189, P190),
#            C07 (N35, N105, D14, IF4, live, E3), C08 (N37, D15, P191, P192, A3, E4)
# =====================================================================================================================
nb.section("8.3", "Elementary Lubrication Theory", intro=r"""
**What is this section about?** Real walls are rarely exactly parallel, but in a thin gap they are *nearly* parallel. We
learn how to decide, by scaling with two different lengths, which terms of the Navier–Stokes equations can be dropped; the
result is a local Couette-plus-Poiseuille profile and one equation for the pressure (the Reynolds equation). It explains why
an oil film under a bearing pad carries tonnes, why viscous flow between glass plates draws the streamlines of ideal flow,
and how a honey drop spreads under its own weight. The same scaling argument becomes the boundary-layer approximation in
Ch. 9.
""")
nb.recap("R13", "The field equations in the gap", r"""
In two dimensions: continuity $\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$ *(6.2)* and the x-component of
(8.1), written out: $\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho
\frac{\partial p}{\partial x}+\frac\mu\rho\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\big)$ *(8.13a)*
— nothing new, just written out.
""", where="Ch. 4 §4.6")
core("C05", "The lubrication approximation", r"""
Why can we throw away inertia in a thin film even when UL/ν is in the thousands?
""", eqs=("8.17a",))
problem(r"""
Under a sliding bearing pad the oil film is 50 µm thick and 5 cm long — a thousand times longer than thick. The oil moves
at metres per second and UL/ν ≈ 4000, which in a pipe would be turbulent. Yet lubrication engineers treat the film as
purely viscous. We want the argument that makes this legitimate — and the number that must be small instead of 1/Re.
""")
note("N23 [B] · The lubrication idea.", r"""
In a narrow passage of gap h and length L ≫ h (Fig. 8.8) the flow is nearly parallel; pressure and viscous forces dominate,
and a gentle curvature of the passage does not matter as long as its radius is much larger than h. Exactly the same
reasoning, with the gap replaced by a thin layer next to a wall, gives the boundary-layer approximation (Ch. 9 §9.1).
""")
idea(r"""
↑ y ~ h = εL (small)          ∂/∂x ~ 1/L   (slow along the gap)
│  ═══════════════════        ∂/∂y ~ 1/h   (fast across it)
│  ──► ──► ──► u ~ U          continuity: U/L ~ v/h  ⇒  v ~ εU
└──────────────────────► x ~ L
""", words=r"Two length scales, and continuity fixes the size of v.")
note("N24 [B] · The y-momentum equation.", r"""
The y-equation is written below in its correct form. ⚠️ The book prints $-\frac1\rho\frac{\partial p}{\partial x}$ here; the
y-equation needs ∂p/∂**y** — the scaled (8.16b) confirms it, and `ch08.lubrication_nondim_sympy()` shows that the printed
form gives a coefficient set that does not match (8.16b) (D10 step 12, check cell below).
""", equation=EQ["8.13b"], ref="8.13b")
P("P188", "anisotropic scaling with two length scales", r"""
When a flow is long and thin, measure x in units of the length L and y in units of the thickness h = εL; derivatives then
carry different sizes, $\partial/\partial x=(1/L)\partial/\partial x^*$ but $\partial/\partial y=(1/\varepsilon L)\partial/
\partial y^*$. The velocity across the layer gets its own scale from continuity: U/L must balance v/h, so v ~ εU. Written
this way every starred quantity is of order one and the size of each term sits in its coefficient.
""", code=r"""
L, eps, U = 0.05, 1e-3, 5.0              # length [m], fineness h/L, speed [m/s]
h = eps*L                                # the gap [m]
print(U/L, (eps*U)/h)                    # both 100 1/s: the two continuity terms balance
""")
remind([
    ("scaled variables and the chain rule", "if $x=Lx^*$ then $\\partial/\\partial x=(1/L)\\,\\partial/\\partial x^*$ (Ch. 4 P133)."),
    ("bookkeeping of a small parameter ε", "multiply the whole equation by the factor that gives the term you expect to dominate the coefficient 1; every other coefficient is then a power of ε times a group, and you read off which terms are small (Ch. 2 P68, orders of smallness)."),
    ("order-of-magnitude scaling", "replace each derivative by (change)/(distance over which it happens), e.g. ∂u/∂y ~ U/h (Ch. 4 P130)."),
    ("sympy symbols and simplify in check cells", "`sp.symbols(..., positive=True)` makes symbols; `sp.simplify(expr) == 0` tests an identity (Ch. 4 P117)."),
])
note("N26 [B] · The scalings (8.14),", r"""
with the fineness ratio ε = h/L and the atmospheric pressure $P_a$ as the book's pressure scale.
""", equation=EQ["8.14"], ref="8.14")
D("D10", ref="8.16a", check_src=r"""
sol = ch08.lubrication_nondim_sympy()                    # D10 steps 1–14 repeated symbolically (cached)
L_, h_, U_, rho_, mu_, Pa_, eps, ReL, Lam = sol["symbols"]   # the symbols it used (ε, Re_L, Λ among them)
print("scaled continuity:", sol["continuity_star"])      # (8.15): no coefficient at all
print("x-coefficients:", sol["x_coeffs"])                # (8.16a): ε²Re_L, 1/Λ, ε², 1
print("y-coefficients:", sol["y_coeffs"])                # (8.16b): ε⁴Re_L, 1/Λ, ε⁴, ε²
assert sp.simplify(sol["x_coeffs"]["inertia"] - eps**2*ReL) == 0      # step 9: inertia weighed by ε²Re_L
assert sp.simplify(sol["x_coeffs"]["pressure"] - 1/Lam) == 0          # step 10: pressure 1/Λ
assert sp.simplify(sol["x_coeffs"]["diff_across"] - 1) == 0           # step 8's choice: cross-gap friction has 1
assert sp.simplify(sol["y_coeffs"]["inertia"] - eps**4*ReL) == 0      # step 14: y-inertia ε⁴Re_L
assert sp.simplify(sol["y_coeffs"]["pressure"] - 1/Lam) == 0          # step 13's choice: the same 1/Λ yardstick
assert sp.simplify(sol["y_coeffs"]["diff_across"] - eps**2) == 0      # step 14: y-friction across the gap ε²
print("with the printed ∂p/∂x in (8.13b):", sol["printed_13b_y_coeffs"])   # pressure ε/Λ: not the set of (8.16b)
assert sp.simplify(sol["printed_13b_y_coeffs"]["pressure"] - 1/Lam) != 0   # the printed slip fails
""")
note("N27 [B], N28 [B], N29 [B] · What D10 delivered.", r"""
**N27** The scaled continuity equation (8.15) has no coefficient: mass is conserved exactly — v is small *because* u varies
slowly. **N28** In (8.16a) the inertia carries $\varepsilon^2\mathrm{Re}_L$ with $\mathrm{Re}_L=\rho UL/\mu$, and the pressure
carries $1/\Lambda$ with the bearing number $\Lambda=\mu UL/(P_ah^2)$. **N29** In (8.16b) only the pressure term lacks a
power of ε — the seed of "pressure is constant across a thin layer" (Ch. 9, Ch. 13).
""", equation=EQ["8.15"], ref="8.15")
D("D11", ref="8.17a")
note("N30 [B] · Pressure is uniform across the gap,", r"""
with an error of relative size ε². The same statement across a thin atmospheric or oceanic layer is the hydrostatic
approximation of Ch. 13.
""", equation=EQ["8.17b"], ref="8.17b")
confusion(r"""
the book prints (8.17a) as $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\frac{\partial^2u}{\partial y^2}$ — without the ν.
The units do not balance (m/s² against 1/(m s)); re-dimensionalising (D11 step 5) restores $\nu\frac{\partial^2u}{\partial y^2}$.
Its own next line, (8.18), carries the 1/μ correctly.
""")
nb.worked_example("an engine-bearing oil film", r"""
h = 50 µm, L = 5 cm, U = 5 m/s, a light oil with μ = 0.05 Pa s and ρ = 870 kg/m³, so ν = μ/ρ = 5.75 × 10⁻⁵ m²/s (numbers
ours).
1. $\varepsilon=h/L=5\times10^{-5}/0.05=10^{-3}$.
2. $\mathrm{Re}_L=UL/\nu=5\times0.05/5.75\times10^{-5}=4350$ — in a pipe this would be turbulent.
3. $\varepsilon^2\mathrm{Re}_L=10^{-6}\times4350=4.35\times10^{-3}$ — small: inertia is under half a percent of friction.
4. The natural pressure scale $\mu UL/h^2=0.05\times5\times0.05/(2.5\times10^{-9})=5\times10^6$ Pa ≈ 49 atm.
5. With the book's atmospheric scale, $\Lambda=\mu UL/(P_ah^2)=5\times10^6/101\,325=49$ — not "near unity".
""")
nb.code(r"""
L, h, U, rho, mu = 0.05, 50e-6, 5.0, 870.0, 0.05          # the engine film: length, gap [m], speed [m/s], ρ, μ
s = ch08.lubrication_scales(L, h, U, rho, mu)             # ε, Re_L, ε²Re_L, Λ (with P_a) and μUL/h²
print({k: f"{v:.4g}" for k, v in s.items()})
keys = ("x_inertia", "x_pressure", "x_diff_along", "x_diff_across", "y_inertia", "y_pressure", "y_diff_along", "y_diff_across")
for scale in ("viscous", "atm"):                          # pressure scaled by μUL/h² (Λ = 1) or by P_a (the book's)
    t = ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale=scale)   # the coefficients of (8.16a), (8.16b)
    print(scale, {k: f"{t[k]:.3g}" for k in keys})
""", explain=r"""
1. `lubrication_scales` returns ε = 10⁻³, Re_L = 4350, ε²Re_L = 4.35 × 10⁻³, Λ = 49.3 and the natural pressure scale
   μUL/h² = 5 × 10⁶ Pa.
2. `lubrication_term_magnitudes` returns exactly the coefficients of (8.16a) and (8.16b). With the viscous pressure scale
   μUL/h² every kept term is O(1): x-pressure 1, cross-gap friction 1, and inertia 4.35 × 10⁻³. With the book's P_a both
   pressure coefficients become 1/Λ = 0.0203 — only the bookkeeping changes, not the physics.
""")
nb.check_agree(r"""
eps = h/L; ReL = rho*U*L/mu; Lam = mu*U*L/(101325*h**2)   # ε, Re_L, Λ by hand (P_a = 101 325 Pa)
mine = [eps**2*ReL, 1/Lam, eps**2, 1.0, eps**4*ReL, 1/Lam, eps**4, eps**2]   # the eight coefficients of D10
lib = ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="atm")
assert np.allclose(mine, [lib[k] for k in keys], rtol=1e-12)   # same numbers, in the order of `keys`
print("hand-made coefficients = library coefficients")
""")
nb.figure(r"""
L, U, rho, mu = 0.05, 5.0, 870.0, 0.05                    # fixed: Re_L = ρUL/μ = 4350
epss = np.geomspace(1e-4, 0.3, 60)                        # fineness ratios ε = h/L
T = [ch08.lubrication_term_magnitudes(L, e*L, U, rho, mu, p_scale="viscous") for e in epss]   # coefficients at each ε
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for ax, pre, title in ((a1, "x", "x-momentum (8.16a)"), (a2, "y", "y-momentum (8.16b)")):
    ax.loglog(epss, [t_[pre + "_inertia"] for t_ in T], color=COLORS["teal"], lw=2, label="inertia")
    ax.loglog(epss, [t_[pre + "_pressure"] for t_ in T], color=COLORS["orange"], lw=2, label="pressure")
    ax.loglog(epss, [t_[pre + "_diff_along"] for t_ in T], color=COLORS["rose"], lw=1.2, ls=":", label="friction along the gap")
    ax.loglog(epss, [t_[pre + "_diff_across"] for t_ in T], color=COLORS["rose"], lw=2, label="friction across the gap")
    e_fail = 1/np.sqrt(4350.0)                            # where ε²Re_L = 1
    ax.axvline(e_fail, color=COLORS["muted"], ls="--", lw=1)
    ax.text(e_fail*1.1, 1e-9, "ε²Re_L = 1:\nlubrication\nfails to the right", fontsize=7)
    ax.axvline(1e-3, color=COLORS["amber"], lw=0.8); ax.text(1.1e-3, 1e2, "engine film", fontsize=7, color=COLORS["amber"])
    ax.set_xlabel("ε = h/L [–]"); ax.set_title(title, fontsize=10); ax.set_ylim(1e-14, 1e4)
a1.set_ylabel("coefficient (size of the term) [–]"); a1.legend(fontsize=7, loc="lower right")
fig.suptitle("In a thin gap inertia is weighed with ε²Re_L, not Re_L (Re_L = 4350, pressure scale μUL/h²)", fontsize=11)
savefig(fig, "ch08", "c05_term_magnitudes"); plt.show()
""", see=r"""Two flat lines on the left (pressure and cross-gap friction) and falling lines for everything else; on the right only the
pressure stays at 1 — every other y-term falls like ε² or faster.""",
    read=r"""Read each term's size at your ε: at ε = 10⁻³ (amber line) inertia is 4.35 × 10⁻³ of the kept terms; only near
ε ≈ 0.015 (dashed line) does it catch up.""",
    change=r"""…Re_L were 100 times larger: the teal line moves up two decades and lubrication fails already at ε ≈ 0.0015.""")
note("N31 [B] · How small is small?", r"""
For our engine film $\varepsilon^2\mathrm{Re}_L=4.35\times10^{-3}$, so the lubrication balance is excellent. (The book's own
example uses a different oil film; its number lives in our private test file.) ⚠️ With $p^*=p/P_a$ the bearing number
$\Lambda=\mu UL/(P_ah^2)$ is 49 here and ≈ 10³ for thinner, longer films — not "near unity" as the text assumes. Nothing in the
ordering changes (Λ only rescales p*), but the physical pressure scale is μUL/h², tens to hundreds of atmospheres.
""")
explainer("lubrication_scaling", "Why can a thin film ignore inertia?", r"""
Drag the gap, the length and the speed and watch every term of the scaled momentum equations (8.16a) and (8.16b) as a bar on
a log scale: inertia falls like ε²Re_L, the pressure and the cross-gap friction stay, and the y-equation leaves only the
pressure — an ordering argument you can see move.""",
          ["Start from the engine-film preset and double U three times — the teal inertia bar climbs but stays far below the pressure bar.",
           "Press 'thick gap': ε = 0.2 and Re_L = 10⁴ put inertia on top and the status turns amber.",
           "Switch the pressure scale to P_a (the book's) and read Λ — then back to μUL/h².",
           "Open the Derivation tab and step through D10: each substitution lights the bar it creates."])
whatif(r"""
…we now use (8.17a) to find the velocity everywhere in a gap of any slowly varying shape? Integrating twice across the gap
gives a local Couette-plus-Poiseuille profile (C06).
""")

core("C06", "The lubrication profile, the gap flux and the Reynolds equation", r"""
Given the shape of the gap, what is the velocity at every point — and what single equation decides the pressure?
""", eqs=("8.19",))
problem(r"""
Under a tilted pad, in a knee joint, in the gap of an injection mould, the gap height changes slowly along the flow. At each
station the flow "thinks" it is between parallel plates — so its profile is the Couette–Poiseuille profile of C02 with the
local gap h(x) and the local pressure gradient. But the pressure gradient is unknown; mass conservation across the whole gap
has to decide it.
""")
idea(r"""
station 1 (wide)       station 2            station 3 (narrow)
══════════ U_h         ══════════ U_h       ══════════ U_h
 ─►                     ──►                   ───►
 ►  + Poiseuille        ─► + smaller          ──►      same flux q must pass every station
═══════════ U₀        ═══════════ U₀         ═══════════ U₀
""", words=r"The channel formula, applied station by station.")
note("N25 [B] · Gap boundary conditions.", r"""
$u=U_0(t)$ on the lower wall y = 0 and $u=U_h(t)$ on the upper wall y = h(x, t); the pressure may depend on t. Time only
enters as a parameter: the lubrication balance (8.17) has no ∂/∂t term, yet h, U₀, U_h and p may all change slowly.
""")
D("D12", ref="8.19")
note("N33 [B] · The intermediate form (8.18)", r"""
has "constants" A, B that may depend on x and t (step 3). ⚠️ The text says (8.16a) "can be integrated twice" — it is the
simplified (8.17a) that is integrated.
""", equation=EQ["8.18"], ref="8.18")
note("N34 [B] · The U₀ term.", r"""
The printed (8.19) adds $U_0$ to the Couette part, so at the upper wall it gives $u(h)=U_h+U_0$ instead of $U_h$. The
consistent profile is $U_h\frac yh+U_0\big(1-\frac yh\big)$ (D12 step 7). Every example of the chapter has U₀ = 0, so no
result changes; our code uses the consistent form and keeps the printed one as `form='book'` for a test that must fail.
""")
nb.code(r"""
h = 50e-6                                                # a 50 µm gap [m]
print(ch08.lubrication_velocity(h, h, 0.0, U_h=2.0, U_0=1.0))               # consistent (8.19) at the top wall: U_h = 2 m/s ✓
print(ch08.lubrication_velocity(h, h, 0.0, U_h=2.0, U_0=1.0, form="book"))  # as printed: U_h + U₀ = 3 m/s ✗
""")
P("P189", "Leibniz rule with a moving upper limit", r"""
Differentiating an integral whose upper limit moves adds a boundary term:
$\frac{\partial}{\partial x}\int_0^{h(x)}u\,dy=\int_0^{h}\frac{\partial u}{\partial x}dy+u(x,h)\frac{\partial h}{\partial x}$.
The second term counts what enters or leaves because the limit itself moved (Ch. 3 met the fixed-limit version).
""", code=r"""
import sympy as sp
x, y = sp.symbols('x y')
h = 1 + x**2                                             # a moving upper limit
u = x*y                                                  # a test integrand
lhs = sp.diff(sp.integrate(u, (y, 0, h)), x)             # d/dx of the integral
rhs = sp.integrate(sp.diff(u, x), (y, 0, h)) + u.subs(y, h)*sp.diff(h, x)   # Leibniz: inside + boundary term
print(sp.simplify(lhs - rhs))                            # 0
""")
remind([
    ("differentiation under the integral sign", "with fixed limits, $\\frac{d}{dx}\\int_a^bf\\,dy=\\int_a^b\\frac{\\partial f}{\\partial x}dy$ (Ch. 3 P109) — the primer above adds the moving-limit term."),
    ("kinematic condition at a moving wall", "a particle on the surface y = h(x, t) stays on it, so D(y − h)/Dt = 0 there: $v=\\partial h/\\partial t+u\\,\\partial h/\\partial x$ at y = h (Ch. 4 P132)."),
])
D("D13", check_src=r"""
r = ch08.reynolds_equation_sympy()                       # D13 built symbolically: profile → flux → integrated continuity
print("q =", r["q"])                                     # −h³p_x/(12μ) + (U₀ + U_h)h/2, over a common denominator
assert sp.simplify(r["leibniz_residual"]) == 0           # step 2: the Leibniz rule holds for the lubrication profile
assert sp.simplify(r["kinematic_cancellation"]) == 0     # steps 5–7: the u(h)·h_x terms cancel
assert sp.simplify(r["reynolds_residual"]) == 0          # step 12: h_t + q_x = 0 is the gap integral of u_x + v_y = 0
print("Reynolds equation:", r["reynolds_equation"])      # the 1-D Reynolds equation, expanded
""")
note("N32 [B] · Completing a lubrication problem.", r"""
The profile alone is not a solution: we still need p(x). Integrating continuity across the gap gives
$\frac{\partial h}{\partial t}+\frac{\partial q}{\partial x}=0$ with the gap flux below — the 1-D **Reynolds equation** (the
book leaves it to Exercises 8.19–8.20). For a steady gap q is the same at every station, which is one ODE for p(x); two end
pressures close it (C07).
""", equation=r"q=\int_0^hu\,dy=-\frac{h^3}{12\mu}\frac{\partial p}{\partial x}+\frac{(U_0+U_h)h}{2}")
P("P190", "scipy.integrate.cumulative_trapezoid", r"""
`cumulative_trapezoid(f, x, initial=0)` returns the running integral $\int_{x_0}^{x_i}f\,dx$ at every sample — the numerical
antiderivative. We use it to turn a pressure gradient into a pressure profile.
""", code=r"""
from scipy.integrate import cumulative_trapezoid
import numpy as np
x = np.linspace(0, 1, 5)
print(cumulative_trapezoid(2*x, x, initial=0))           # ≈ x²: [0, 0.0625, 0.25, 0.5625, 1.0]
""")
nb.worked_example("the flux through a 50 µm gap", r"""
h = 50 µm, upper wall $U_h$ = 5 m/s, lower wall fixed, μ = 0.05 Pa s.
1. With no pressure gradient: $q=U_hh/2=5\times5\times10^{-5}/2=1.25\times10^{-4}$ m²/s.
2. The gradient that halves it: $\frac{h^3}{12\mu}\frac{dp}{dx}=\frac{U_hh}{4}$ ⇒ $\frac{dp}{dx}=\frac{3\mu U_h}{h^2}=\frac{3\times0.05\times5}
   {2.5\times10^{-9}}=3\times10^8$ Pa/m — 3 bar per millimetre.
3. At $dp/dx=6\mu U_h/h^2=6\times10^8$ Pa/m the flux is zero: the pressure-driven backflow exactly cancels the drag.
""")
nb.code(r"""
h, Uh, mu = 50e-6, 5.0, 0.05                             # gap [m], upper-wall speed [m/s], oil viscosity [Pa s]
G_half = 3*mu*Uh/h**2                                    # the gradient that halves the flux [Pa/m]
print(ch08.lubrication_flux(h, 0.0, U_h=Uh, mu=mu), ch08.lubrication_flux(h, G_half, U_h=Uh, mu=mu))   # q [m²/s]
yy = np.linspace(0, h, 5)                                # five heights across the gap [m]
print(ch08.lubrication_velocity(yy, h, G_half, U_h=Uh, mu=mu))   # (8.19) at the halving gradient [m/s]
x = np.linspace(0, 0.05, 401); hx = h*(1 + 0.5*x/0.05)   # a slider gap h₀(1 + αx/L) with α = 0.5, L = 5 cm
p, q = ch08.reynolds_pressure_1d(x, hx, U_0=-Uh, mu=mu)  # pad frame: floor slides at −U under the resting pad; p = 0 at both ends
print(f"p_max = {p.max()/1e6:.3f} MPa at x/L = {x[np.argmax(p)]/0.05:.3f}; flux q = {q:.4g} m²/s")
""", explain=r"""
1. `lubrication_flux` is N32's q: 1.25 × 10⁻⁴ m²/s with no pressure gradient, half of it at the "halving" gradient.
2. The profile at that gradient is [0, −0.156, 0.625, 2.344, 5.0] m/s: a thin reversed layer at the fixed wall, because
   $3\mu U/h^2$ exceeds the backflow threshold $2\mu U/h^2$ of C02.
3. `reynolds_pressure_1d` solves the steady Reynolds equation for a gap given as samples: it finds the constant flux that makes
   the two end pressures equal and integrates $dp/dx=12\mu(\bar q-q)/h^3$ (with $\bar q=(U_0+U_h)h/2$) using
   `cumulative_trapezoid`. For the α = 0.5 pad of C07 in the pad frame it gives a 1.00 MPa hump at x/L = 0.40 and the pad-frame
   flux q = −1.5 × 10⁻⁴ m²/s — the $C_1=-\frac{1+\alpha}{2+\alpha}Uh_o$ of D14 (cross-checked there).
""")
nb.check_agree(r"""
q_num = integrate.quad(lambda y_: ch08.lubrication_velocity(y_, h, 1e8, U_h=Uh, mu=mu), 0, h)[0]   # ∫₀ʰ u dy numerically
print(f"quad: {q_num:.6g} m²/s, formula: {ch08.lubrication_flux(h, 1e8, U_h=Uh, mu=mu):.6g} m²/s")
assert np.allclose(q_num, ch08.lubrication_flux(h, 1e8, U_h=Uh, mu=mu), rtol=1e-10)   # N32's flux is the integral of (8.19)
""")
nb.figure(r"""
h0, al, L, U, mu = 50e-6, 1.5, 0.05, 5.0, 0.05            # the slider gap of C07 with a steep taper α = 1.5
k = 200.0                                                 # vertical exaggeration of the drawing
fig, ax = plt.subplots(figsize=(10, 3.4))
xs = np.linspace(0, L, 200)                               # along the pad [m]
ax.fill_between(1e2*xs, k*h0*(1 + al*xs/L)*1e2, k*h0*(1 + al*xs/L)*1e2 + 1.2, color=COLORS["muted"], alpha=0.6)   # the pad
ax.axhline(0, color=COLORS["ink"], lw=2)                  # the floor
for x0 in np.linspace(0.04*L, 0.96*L, 7):                 # seven stations
    hx = h0*(1 + al*x0/L)                                 # local gap [m]
    yy = np.linspace(0, hx, 25)                           # heights across the gap [m]
    u = ch08.slider_gap_velocity(x0, yy, h0, al, L, U, mu=mu, frame="ground")   # (8.19) with the exact slider dp/dx [m/s]
    profile_arrows(ax, 1e2*k*yy, 1e2*0.12*u/U*L, x0=1e2*x0, every=4, color=COLORS["blue"])   # arrows ∝ u
for x0 in (0.04*L, 0.96*L):                               # the two parts at the narrowest and widest stations
    hx = h0*(1 + al*x0/L); yy = np.linspace(0, hx, 25)
    uc = U*yy/hx                                          # Couette part (the pad drags at U) [m/s]
    up = ch08.slider_gap_velocity(x0, yy, h0, al, L, U, mu=mu, frame="ground") - uc   # Poiseuille part [m/s]
    ax.plot(1e2*(x0 + 0.12*uc/U*L), 1e2*k*yy, color=COLORS["rose"], ls="--", lw=1)
    ax.plot(1e2*(x0 + 0.12*up/U*L), 1e2*k*yy, color=COLORS["orange"], ls="--", lw=1)
ax.set_xlabel("x along the pad [cm]"); ax.set_ylabel(f"height × {k:.0f} [cm]")
ax.set_title("Couette plus Poiseuille, station by station (heights exaggerated 200×; pad moving at U, α = 1.5)", fontsize=10)
ax.plot([], [], color=COLORS["blue"], label="u(y) (ground frame)"); ax.plot([], [], color=COLORS["rose"], ls="--", label="Couette part")
ax.plot([], [], color=COLORS["orange"], ls="--", label="Poiseuille part"); ax.legend(fontsize=8, loc="upper left")
savefig(fig, "ch08", "c06_gap_profiles"); plt.show()
""", see=r"""Profiles in a gap that widens from h₀ at the left (the pad's trailing edge) to 2.5h₀ at the right (its leading
edge), the pad moving right at U over a fixed floor: near the narrow end the profiles sag below the straight Couette line (the
orange Poiseuille part points backward), near the wide end they bulge forward. The dashed rose and orange curves are the two
parts at the end stations — our remake of Fig. 8.8.""",
    read=r"""The pressure hump pushes fluid backward where it rises (from the narrow end up to the peak at x = L/(2 + α) ≈ 1.4 cm) and
forward where it falls. In the ground frame the flux grows with the gap, $q=C_1+Uh$, because the pad-frame flux $C_1$ is the
same at every station (D14).""",
    change=r"""…the pad slid the other way (U < 0): every profile flips, and the pressure hump becomes a suction dip (C07).""")
note("N36 [B] · Hele-Shaw flow (Example 8.2): viscous flow that draws ideal streamlines.", r"""
Between two plates z = 0 and z = h the lubrication balances $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}
{\partial z^2}$, $0\cong-\frac1\rho\frac{\partial p}{\partial y}+\nu\frac{\partial^2v}{\partial z^2}$,
$0\cong-\frac1\rho\frac{\partial p}{\partial z}$ give the velocity below, with the potential $\phi=-\frac{z(h-z)}{2\mu}p$.
Integrating continuity over the gap (the 2-D Reynolds equation with constant h, depth-averaged velocity
$\bar{\mathbf u}=-\frac{h^2}{12\mu}\nabla p$) leaves $\frac{\partial^2p}{\partial x^2}+\frac{\partial^2p}{\partial y^2}=0$ — so φ
obeys $u\equiv\partial\phi/\partial x,\ v\equiv\partial\phi/\partial y$ *(6.10)* and $\nabla^2\phi=0$ *(6.12)*, the equations of
2-D ideal flow (Ch. 6). Dye injected between glass plates therefore traces ideal-flow streamlines, except in layers of
thickness ~h at obstacles, where no slip holds (Exercise 8.34). ⚠️ The book prints these lubrication equations without the ν,
and writes the wall conditions at "y = 0, h" — the gap coordinate here is z.
""", equation=r"u\cong-\frac{1}{2\mu}\frac{\partial p}{\partial x}z(h-z)=\frac{\partial\phi}{\partial x},\qquad \phi=-\frac{z(h-z)}{2\mu}p")
nb.figure(r"""
a, gap, Um = 0.01, 1e-3, 1e-3                             # disc radius [m], plate gap [m], mean speed [m/s]
xg = np.linspace(-0.04, 0.04, 241 if not FAST else 121)   # a 8 cm × 8 cm window [m]
X, Y = np.meshgrid(xg, xg)
hs = ch08.hele_shaw_cylinder(X, Y, None, Um, a, gap)      # depth-averaged Hele-Shaw flow (z = None → gap average)
psi_ideal = PF.stream_function(PF.cylinder(Um, a), X, Y)  # Ch. 6: ideal flow past a circular cylinder
inside = X**2 + Y**2 < a**2                               # the disc
lev = np.linspace(-3.5e-5, 3.5e-5, 15)                    # streamline levels [m²/s]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw=dict(width_ratios=[1.4, 1]))
a1.contour(X*1e2, Y*1e2, np.where(inside, np.nan, hs["psi_mean"]), lev, colors=COLORS["blue"], linewidths=1.6)
a1.contour(X*1e2, Y*1e2, np.where(inside, np.nan, psi_ideal), lev, colors=COLORS["muted"], linewidths=1, linestyles="--")
a1.add_patch(plt.Circle((0, 0), a*1e2, color=COLORS["muted"]))
a1.set_aspect("equal"); a1.set_xlabel("x [cm]"); a1.set_ylabel("y [cm]")
a1.set_title("Hele-Shaw dye lines (blue) on ideal streamlines (dashed)", fontsize=10)
z = np.linspace(0, gap, 50)                               # across the gap [m]
dpdx = -12*MU_W*Um/gap**2                                 # the far-field gradient that drives the mean speed Um [Pa/m]
u_z, _ = ch08.hele_shaw_velocity(z, gap, (dpdx, 0.0))     # the parabola across the gap [m/s]
a2.plot(1e3*u_z, 1e3*z, color=COLORS["blue"], lw=2); a2.axvline(1e3*Um, color=COLORS["muted"], ls="--", lw=1)
a2.text(1e3*Um*1.02, 0.1, "mean", fontsize=8, color=COLORS["muted"])
a2.set_xlabel("u [mm/s]"); a2.set_ylabel("z across the gap [mm]"); a2.set_title("Profile across the gap (far away)", fontsize=10)
savefig(fig, "ch08", "c06_hele_shaw"); plt.show()
""", see=r"""The blue gap-averaged dye lines of Hele-Shaw flow lie exactly on the dashed ideal-flow streamlines round a cylinder,
while across the 1 mm gap the flow is a viscous parabola — our remake of Fig. 8.10.""",
    read=r"""The pressure is harmonic, so the depth-averaged velocity is a potential flow; only a thin layer at the disc (too thin to
see here) differs, where no slip holds.""",
    change=r"""…the gap were as wide as the disc: the ε ≪ 1 assumption fails, inertia and the no-slip layer spread, and the pattern
loses its fore–aft symmetry.""")
nb.md(r"""
**A numerical check (our route, not the book's).** `hele_shaw_streamfunction_grid` solves Laplace's equation for the
gap-averaged stream function on a square grid with the disc cut out (Ch. 6's masked five-point solver) and compares it with
the ideal-flow value away from the disc. The disc boundary is a staircase of grid cells, so the scheme is only **first
order**: doubling the grid roughly halves the error.
""")
nb.code(r"""
ns = (33, 65, 129) if not FAST else (33, 65)             # nodes per side of the square grid
errs, hs_ = [], []                                        # far-field error and grid step for each grid
for n_ in ns:
    g_ = ch08.hele_shaw_streamfunction_grid(n=n_, a=1.0, U_mean=1.0, box=4.0)   # disc radius 1, box ±4 (non-dimensional)
    errs.append(g_["err_far"]); hs_.append(g_["h"])        # max |ψ − ψ_ideal| where r ≥ 2a, and the step
    print(f"n = {n_}: step {g_['h']:.3f}, max error {g_['err_far']:.2e}")
print(f"observed order ≈ {observed_order(hs_, errs):.2f}  (first order: the staircase disc)")
""")
whatif(r"""
…the gap were the space under a real, tilted pad with the same pressure $p_e$ at both ends? The flux must be constant along
the pad, which forces a pressure hump inside — and the hump carries a load (C07).
""")

core("C07", "The slider bearing (Example 8.1)", r"""
How does a thin viscous film carry a load — and why only in one direction?
""")
nb.md(r"*In one line:* $p-p_e=\frac{6\mu LU}{h_o^2}\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha x/L)^2}$ and, for a small taper, "
      r"$W=\frac{\alpha\mu L^2U}{2h_o^2}$ (Example 8.1, denominator corrected to the square).")
problem(r"""
A tilted pad skates on a film of oil over a flat surface, as in the thrust bearing of a ship's propeller shaft or a
hydroelectric turbine carrying hundreds of tonnes. The pad never touches the surface. Where does the upward force come from,
and why does the bearing fail — the pad is sucked down — if it slides the wrong way?
""")
idea(r"""
            W ↓            pad moves → at U, wide end leading   (pad frame: floor moves ← at U)
p_e  ┌────────────────────────┐  p_e
     │ h₀ (x = 0)  →  h₀(1+α) (x = L) │
═════════════════════════════════   the floor drags oil in at the wide end, the narrow end lets less out
p(x):  ____/‾‾‾\____                ⇒ p rises until Poiseuille backflow evens the flux; ∫(p − p_e)dx = W
""", words=r"Same flux through a narrowing gap ⇒ the pressure must rise inside.")
remind([
    ("antiderivatives of (1 + αx/L)⁻ⁿ", "substitute s = 1 + αx/L, ds = (α/L)dx: $\\int(1+\\alpha x/L)^{-n}dx=\\frac L\\alpha\\frac{(1+\\alpha x/L)^{1-n}}{1-n}$ for n ≠ 1 (Ch. 3 P106, substitution in an integral)."),
    ("Taylor expansion in a small parameter α", "for small α, $(1+\\alpha x/L)^{-2}\\approx1-2\\alpha x/L$ and $2+\\alpha\\approx2$: keep the first power of α (Ch. 1 P26, first-order Taylor)."),
])
D("D14", check_src=r"""
sb = ch08.slider_bearing_sympy()                         # D14 steps 5–15 rebuilt symbolically (cached)
x_, L_, h0_, mu_, U_, al_, pe_ = sb["symbols"]           # x, L, h₀, μ, U, α, p_e
print("C1 =", sp.simplify(sb["C1"]))                     # step 10: −(1 + α)U h₀/(2 + α)
print("C2 =", sp.simplify(sb["C2"]))                     # step 11
assert sp.simplify(sb["ode_residual_exact"]) == 0        # the squared denominator satisfies dp/dx = −12μC₁/h³ − 6μU/h² (step 5)
assert sp.simplify(sb["ode_residual_book"]) != 0         # the printed first-power form does not
assert all(sp.simplify(r_) == 0 for r_ in sb["bc_residuals"])   # p(0) = p(L) = p_e
print("W_linear =", sb["W_linear"])                      # step 15: αμL²U/(2h₀²)
print("W_exact  =", sb["W_exact_series"])                # N35's exact load expanded in α: starts with W_linear
""")
confusion(r"""
**two printing slips in Example 8.1.** The intermediate integrals are printed with $(1-\alpha x/L)$, although the gap is
$h_o(1+\alpha x/L)$; and the final exact pressure is printed with $(1+\alpha x/L)$ to the first power in the denominator.
Substituting back into $dp/dx=-12\mu C_1/h^3-6\mu U/h^2$ shows the denominator must be **squared** (D14 step 14, checked
above). The printed C₁, C₂, the linear-α pressure and the load W are right.
""")
nb.worked_example("a 5 cm pad on a 50 µm film", r"""
L = 0.05 m, h₀ = 50 µm, U = 5 m/s, μ = 0.05 Pa s, α = 0.1.
1. Scale: $6\mu LU/h_o^2=6\times0.05\times0.05\times5/(2.5\times10^{-9})=3\times10^7$ Pa.
2. At the middle x = L/2: $\frac{\alpha(x/L)(1-x/L)}{(2+\alpha)(1+\alpha/2)^2}=\frac{0.1\times0.25}{2.1\times1.1025}=0.0108$ ⇒
   p − p_e = 3.24 × 10⁵ Pa ≈ 3.2 atm.
3. Linear formula: $\frac{3\alpha\mu LU}{h_o^2}\times0.25=3.75\times10^5$ Pa.
4. Load $W=\alpha\mu L^2U/(2h_o^2)=0.1\times0.05\times0.0025\times5/(5\times10^{-9})=1.25\times10^4$ N/m — 12.5 kN per metre of
   pad width, about 1.3 tonnes.
5. Reverse the motion (U = −5 m/s): W = −12.5 kN/m, the pad is pulled down.
""")
remind([
    ("scipy.optimize.minimize_scalar", "`optimize.minimize_scalar(f, bounds=(a, b), method=\"bounded\")` finds the minimum of a one-variable function; to maximise W we minimise −W (Ch. 7 P170)."),
    ("np.log1p and cancellation for small α", "`np.log1p(a)` computes ln(1 + a) accurately for tiny a, where 1 + a loses digits — the same cancellation problem as `np.expm1` (Ch. 3 P107); the exact load uses a series for |α| < 10⁻³."),
])
nb.code(r"""
L, h0, U, mu = 0.05, 50e-6, 5.0, 0.05                    # the pad of the tiny example
for model in ("exact", "linear", "book"):                # corrected exact, O(α), and the printed slip (a ghost)
    print(f"{model:6s} p(L/2) − p_e = {ch08.slider_bearing(0.025, h0, 0.1, L, U, mu=mu, model=model):.4g} Pa")
st = ch08.slider_bearing_state(h0, 0.1, L, U, mu=mu)     # everything the explainer shows
print({k: (round(st[k], 7) if not isinstance(st[k], bool) else st[k]) for k in
       ("C1", "p_max", "x_pmax", "W_exact", "W_linear", "err_linear", "p_max_atm", "inlet_backflow")})
print(ch08.slider_optimum_taper())                       # the taper that maximises the exact load
print("reversed U:", ch08.slider_bearing_load(h0, 0.1, L, -U, mu=mu, model="linear"), "N/m")   # suction
""", explain=r"""
1. The three models at mid-pad: exact 3.239 × 10⁵ Pa, linear 3.75 × 10⁵ Pa, and the printed first-power form 3.401 × 10⁵ Pa —
   5 % off, and (D14 check) it does not satisfy the pressure equation.
2. The state: the pad-frame flux C₁ = −1.310 × 10⁻⁴ m²/s, the peak 3.247 × 10⁵ Pa at x = 2.38 cm (x/L = 0.476, just before
   mid-pad), the exact load 10.8 kN/m against the linear 12.5 kN/m — the linear one 15.6 % high — the peak in atmospheres, and no inlet
   recirculation (α ≤ 1).
3. `slider_optimum_taper` maximises the exact load: α = 1.189, inlet/outlet gap ratio 1 + α = 2.189, dimensionless load
   W h₀²/(6μUL²) = 0.02671 (the San Andrés benchmark).
4. Reversing U turns the load negative: the pad is sucked down.
""")
nb.check_agree(r"""
h_ = lambda x_: h0*(1 + 0.1*x_/L)                        # the gap h₀(1 + αx/L), α = 0.1 [m]
dpdx = lambda x_, C1: -12*mu*C1/h_(x_)**3 - 6*mu*U/h_(x_)**2   # D14 step 5 (ground frame) [Pa/m]
C1 = optimize.brentq(lambda C: integrate.quad(lambda x_: dpdx(x_, C), 0, L)[0], -U*h0, 0.0)   # choose C₁ so p(L) = p(0)
x = np.linspace(0, L, 401)                               # stations [m]
p_mine = integrate.cumulative_trapezoid(dpdx(x, C1), x, initial=0)   # p − p_e by integrating the ODE
print(f"C1 by root finding = {C1:.6g} m²/s (library {st['C1']:.6g}); max |Δp| = {np.max(np.abs(p_mine - ch08.slider_bearing(x, h0, 0.1, L, U, mu=mu))):.2f} Pa")
assert np.allclose(C1, st["C1"], rtol=1e-8)              # same flux constant
assert np.allclose(p_mine, ch08.slider_bearing(x, h0, 0.1, L, U, mu=mu), rtol=1e-6, atol=1.0)   # same pressure hump
""")
nb.figure(r"""
L, h0, U, mu = 0.05, 50e-6, 5.0, 0.05
x = np.linspace(0, L, 400)                                # along the pad [m]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
for al, c in ((0.1, COLORS["orange"]), (0.5, COLORS["amber"]), (1.189, COLORS["rose"])):   # three tapers
    a1.plot(x/L, ch08.slider_bearing(x, h0, al, L, U, mu=mu)/1e6, color=c, lw=2, label=f"α = {al:g}, exact")
    a1.plot(x/L, ch08.slider_bearing(x, h0, al, L, U, mu=mu, model="linear")/1e6, color=c, lw=1, ls="--")
a1.plot(x/L, ch08.slider_bearing(x, h0, 0.5, L, U, mu=mu, model="book")/1e6, color=COLORS["muted"], lw=1.5, ls=":",
        label="α = 0.5 as printed (first power)")
a1.plot([], [], color=COLORS["ink"], ls="--", lw=1, label="linear in α")
a1.set_xlabel("x/L [–]"); a1.set_ylabel("p − p_e [MPa]"); a1.legend(fontsize=7); a1.set_title("Pressure humps", fontsize=10)
als = np.linspace(-0.9, 3.0, 200)                         # tapers, including negative (pad sliding the wrong way)
W_ex = np.array([ch08.slider_bearing_load(h0, a_, L, U, mu=mu) for a_ in als])       # N35 exact load [N/m]
W_li = np.array([ch08.slider_bearing_load(h0, a_, L, U, mu=mu, model="linear") for a_ in als])
opt = ch08.slider_optimum_taper()
a2.plot(als, W_ex/1e3, color=COLORS["blue"], lw=2, label="exact (N35)")
a2.plot(als, W_li/1e3, color=COLORS["blue"], lw=1, ls="--", label="linear αμL²U/(2h₀²)")
a2.axvline(opt["alpha_opt"], color=COLORS["muted"], lw=0.8)
a2.text(opt["alpha_opt"]*1.03, 5, f"optimum α = {opt['alpha_opt']:.3f}", fontsize=8)
a2.fill_between(als, W_ex/1e3, 0, where=als < 0, color=COLORS["amber"], alpha=0.3, label="W < 0: pad sucked down")
a2.axhline(0, color=COLORS["muted"], lw=0.8); a2.set_ylim(-60, 100)
a2.set_xlabel("taper α [–]"); a2.set_ylabel("load W [kN/m]"); a2.legend(fontsize=7); a2.set_title("Load against taper", fontsize=10)
fig.suptitle("A pressure hump carries the load — best at inlet/outlet gap ratio 2.19", fontsize=11)
savefig(fig, "ch08", "c07_slider_pressure"); plt.show()
""", see=r"""(a) pressure humps peaking just before mid-pad and growing with α, their linear-α parabolas (dashed) and the printed
first-power form for α = 0.5 (dotted grey) missing the exact hump; (b) the exact load curve bending over and peaking near
α = 1.19, while the linear formula keeps rising — our remake of Fig. 8.9 `N105`.""",
    read=r"""The load is the area under a hump; the linear formula is fine only for α ≲ 0.1 — 15.6 % high at 0.1, about 90 % high at 0.5;
negative α (shaded) means the pad slides with its narrow end leading and is pulled down.""",
    change=r"""…the film were halved to 25 µm: every pressure and the load grow 4× (∝ 1/h₀²) — the bearing is "stiff": extra load
squeezes the film and raises the capacity.""")
note("N35 [B] · The exact load for any taper (ours; the book stops at linear α).", r"""
Integrating the exact pressure gives W below, which tends to $\alpha\mu L^2U/(2h_o^2)$ as α → 0 (the series of ln(1 + α)). It
peaks at an inlet/outlet gap ratio 1 + α = 2.189 — the classic optimum of lubrication texts (San Andrés, benchmark V5). With
our numbers: α = 0.1 → 12.5 kN/m (linear) vs 10.8 kN/m (exact); α = 0.5 → 62.5 vs 32.8 kN/m.
""", equation=r"W=\frac{6\mu UL^2}{h_o^2\alpha^2}\Big[\ln(1+\alpha)-\frac{2\alpha}{2+\alpha}\Big]")
nb.code(r"""
W_formula = ch08.slider_bearing_load(h0, 0.5, L, U, mu=mu)     # N35 for α = 0.5 [N/m]
W_quad = integrate.quad(lambda xx: ch08.slider_bearing(xx, h0, 0.5, L, U, mu=mu), 0, L)[0]   # ∫(p − p_e)dx numerically
print(f"W (formula) = {W_formula:.5g} N/m, W (quad) = {W_quad:.5g} N/m")   # both 3.279e4
""")
nb.plotly(r"""
L, h0, U, mu = 0.05, 50e-6, 5.0, 0.05
xs = np.linspace(0, L, 150)                              # stations [m]


def frame(al):                                            # curves for one taper α
    return {"exact": (xs/L, ch08.slider_bearing(xs, h0, al, L, U, mu=mu)/1e6),
            "linear in α": (xs/L, ch08.slider_bearing(xs, h0, al, L, U, mu=mu, model="linear")/1e6),
            "as printed (first power)": (xs/L, ch08.slider_bearing(xs, h0, al, L, U, mu=mu, model="book")/1e6)}


als = np.linspace(-0.8, 3.0, 20 if not FAST else 10)      # slider values α [–]
fig = slider_figure(frame, "α", als, unit="", xlabel="x/L [–]", ylabel="p − p_e [MPa]", title="", yrange=[-2.0, 2.5])
titles = []                                               # the two loads in each title
for al in als:
    s_ = ch08.slider_bearing_state(h0, al, L, U, mu=mu)
    titles.append(f"α = {al:.2f}: W exact = {s_['W_exact']/1e3:.1f} kN/m, W linear = {s_['W_linear']/1e3:.1f} kN/m")
step_titles(fig, titles)
recolor(fig, {"exact": COLORS["orange"], "linear in α": COLORS["orange"], "as printed (first power)": COLORS["muted"]},
        {"linear in α": "dash", "as printed (first power)": "dot"})
fig.show()
""", explain=r"""
Each slider step is one taper α from −0.8 to 3; the title prints the exact and linear loads. For α < 0 the hump becomes a
suction dip; for large α the linear formula wildly overestimates the load; the grey dotted "as printed" curve meets both end
pressures but is the wrong shape.
""")
remind([
    ("live widgets", "`live(fn, name=(min, max, step), …)` makes sliders that call fn when released — only with a running kernel; the plotly figure above carries the same idea on the web page (Ch. 1 P47)."),
])
nb.live(r"""
def bearing(h0_um=50.0, alpha=0.5, U=5.0, mu=0.05):      # film thickness [µm], taper, sliding speed [m/s], oil viscosity [Pa s]
    s_ = ch08.slider_bearing_state(h0_um*1e-6, alpha, 0.05, U, mu=mu)   # a 5 cm pad
    print(f"p_max = {s_['p_max']/1e5:.3g} bar at x/L = {s_['x_pmax']/0.05:.3f}; W exact = {s_['W_exact']/1e3:.3g} kN/m, "
          f"W linear = {s_['W_linear']/1e3:.3g} kN/m ({s_['err_linear']:.1f} % off); inlet backflow: {s_['inlet_backflow']}")


live(bearing, h0_um=(10, 200, 5), alpha=(-0.8, 3.0, 0.05), U=(-10, 10, 0.5), mu=(0.005, 0.5, 0.005))
""")
explainer("slider_bearing", "How does a film thinner than a hair carry a load?", r"""
Drag the taper and flip the sliding direction: the pressure hump, the station-by-station Couette-plus-Poiseuille profiles and
the load change together, the printed-slip ghost visibly departs from the exact hump, and past α = 1 a small recirculation
appears at the wide inlet — the whole of Example 8.1 in one picture.""",
          ["Set α = 0.1 and compare the exact and linear humps; then α = 1.19 (the optimum) and read the load.",
           "Flip U to negative: the hump turns into suction and the status says the pad is sucked down.",
           "Choose 'compare': the printed curve meets both end pressures but sits above the exact one — it fails the pressure equation.",
           "Push α past 1 and look at the inlet profiles in the pad frame: a backward-moving layer appears under the pad."])
whatif(r"""
…there were no upper wall at all — a free surface — and gravity, not a pressure applied at the ends, pushed the fluid? The
thin-layer balance still holds, with the hydrostatic pressure of the layer itself (C08).
""")

core("C08", "The thin-film equation (Example 8.3)", r"""
What equation governs a thin viscous layer spreading under its own weight?
""")
nb.md(r"*In one line:* $\frac{\partial h}{\partial t}=\frac{\rho g}{3\mu}\frac{\partial}{\partial x}\Big(h^3\frac{\partial h}{\partial x}\Big)$ "
      r"(Example 8.3) — a diffusion equation whose diffusivity $\rho gh^3/3\mu$ depends on the thickness itself.")
problem(r"""
Honey poured on a plate, paint on a wall, a lava lobe, a mud flow, and — with a stiffer, nonlinear rheology — an ice sheet: a
thin layer of very viscous fluid spreads because its own weight presses harder where it is thicker. We want one equation for
the thickness h(x, t) that predicts the spreading.
""")
idea(r"""
     h ↑   ____
          /    \      p = ρg(h − y): higher under the crest
         /      \     flux q = −(ρg/3μ) h³ ∂h/∂x   ← h³: a thin film barely moves
════════╱════════╲════════
""", words=r"Thick places have higher hydrostatic pressure, so fluid is pushed from thick to thin; the flux is choked as the layer thins.")
remind([
    ("hydrostatic pressure in a thin layer, stress-free surface", "at a free surface with negligible air drag the shear stress vanishes, μ∂u/∂y = 0; in a thin, slowly varying layer the vertical balance is hydrostatic, p = p_a + ρg(h − y) — valid because the slope is small (the lubrication ordering of C05), not merely because accelerations are ignored (Ch. 1 P20, boundary conditions)."),
])
D("D15")
P("P191", "nonlinear diffusion in flux form", r"""
$\partial h/\partial t=\partial/\partial x\,(D(h)\,\partial h/\partial x)$ is a diffusion equation whose diffusivity depends
on the unknown — here $D=\rho gh^3/(3\mu)$, large where the layer is thick, tiny where it is thin. Written as
$\partial h/\partial t+\partial q/\partial x=0$ with the flux $q=-D\,\partial h/\partial x$, integrating over x shows the total
volume ∫h dx never changes when no fluid leaves the ends.
""", code=r"""
import numpy as np
rho, g, mu = 1260.0, 9.81, 1.0                           # glycerol
for h in (0.01, 0.005, 0.001):                           # layer thickness [m]
    print(h, rho*g*h**3/(3*mu))                          # effective diffusivity [m²/s] falls as h³
""")
remind([
    ("explicit stepping and its stability limit", "an explicit (forward-Euler) step of a diffusion equation is stable only for Δt ≤ Δx²/(2D) (Ch. 1 P30) — with D ∝ h³ that limit changes as the layer thins."),
])
P("P192", "implicit time stepping with Picard iteration", r"""
Explicit steps of a stiff diffusion equation need a tiny time step. Backward Euler evaluates the right-hand side at the *new*
time, which is stable for any step; when the diffusivity depends on the unknown we guess it from the last iterate, solve the
linear system, update the guess and repeat until it stops changing (Picard iteration; our solver uses the faster Newton
version of the same fixed point). A thin "precursor film" $h_{min}$ keeps the diffusivity positive ahead of the front.
""", code=r"""
import numpy as np
D, dt, dx = 1.0, 0.1, 0.1                                # an explicit step would need dt ≤ dx²/(2D) = 0.005
lam = D*dt/dx**2                                         # = 10: explicit would blow up
print(1/(1 + 4*lam))                                     # backward-Euler damping of the shortest mode: 0.024 < 1, stable
""")
nb.worked_example("a glycerol bead 1 cm high", r"""
ρ = 1260 kg/m³, μ = 1 Pa s, g = 9.81 m/s².
1. Coefficient $\rho g/3\mu=1260\times9.81/3=4.12\times10^3$ m⁻¹s⁻¹.
2. Flux where h = 1 cm and the slope is ∂h/∂x = −0.5: $q=-4.12\times10^3\times10^{-6}\times(-0.5)=2.06\times10^{-3}$ m²/s,
   flowing outward.
3. Where h = 1 mm with the same slope the flux is 1000 times smaller — the edge barely moves.
4. Effective diffusivity at the crest $\rho gh^3/3\mu=4.1\times10^{-3}$ m²/s, falling as h³ while the bead thins.
""")
remind([
    ("np.where", "`np.where(cond, a, b)` picks a where the condition holds and b elsewhere — here a box-shaped initial bead (Ch. 1 P46)."),
])
nb.code(r"""
x = np.linspace(-0.2, 0.2, 400 if not FAST else 200)     # nodes across a 40 cm plate [m]
h0 = np.where(np.abs(x) < 0.01, 0.01, 0.0)               # a box 1 cm high and 2 cm wide: area 2 × 10⁻⁴ m² (half-area 10⁻⁴)
t_out = np.r_[0.0, np.geomspace(0.1, 1000, 24 if not FAST else 12)]   # output times on a log clock [s]
run = ch08.thin_film_spread(h0, x, t_out, rho=1260.0, mu=1.0)   # backward Euler + Newton, conservative flux form
print(f"volume: first {run['volume'][0]:.6e}, last {run['volume'][-1]:.6e} m² (relative change {abs(run['volume'][-1]/run['volume'][0] - 1):.1e})")
for t_ in (10, 100, 1000):                               # compare the front with the similarity law of C10
    i = np.argmin(np.abs(run["t"] - t_))                 # the output time closest to t_
    xN = ch08.viscous_current_similarity(0.0, run["t"][i], 1e-4, rho=1260.0, mu=1.0, return_front=True)[1]
    print(f"t = {run['t'][i]:7.1f} s: numerical front {100*run['x_front'][i]:.2f} cm, similarity x_N {100*xN:.2f} cm")
print("flux at h = 1 cm, slope −0.5:", ch08.thin_film_flux(0.01, -0.5, rho=1260.0, mu=1.0), "m²/s")
""", explain=r"""
1. `thin_film_spread` marches $\partial h/\partial t=\frac{\rho g}{3\mu}\partial_x(h^3\partial_xh)$ in the conservative flux form
   (flux differences between cells), so the volume is conserved to round-off (the printed relative change); the precursor film
   of 1 µm adds a constant 4 × 10⁻⁷ m² to it.
2. The steps are implicit (backward Euler, Newton iterations), so large time steps are stable although the diffusivity is
   stiff.
3. The fronts approach the similarity law of C10, $x_N=\eta_N(\beta A^3t)^{1/5}$: a few percent ahead at 10 s (the box has not
   been forgotten yet), within 2 % at 1000 s.
4. `thin_film_flux` is the tiny example's 2.06 × 10⁻³ m²/s.
""")
nb.check_agree(r"""
xe = np.linspace(-0.2, 0.2, 200); dxe = xe[1] - xe[0]    # 200 cells for the explicit version [m]
he = np.where(np.abs(xe) < 0.01, 0.01, 0.0)              # the same box
beta = 1260.0*ch08.G0/(3*1.0)                            # ρg/3μ [1/(m s)] (ch08 uses g = 9.80665 m/s²)
dte = 0.2*dxe**2/(beta*he.max()**3)                      # explicit stability: well below Δx²/(2D_max) [s]
vol0 = he.sum()*dxe                                      # the starting volume (area per unit width) [m²]
for _ in range(2000):                                    # 2000 explicit finite-volume steps
    hf = 0.5*(he[1:] + he[:-1])                          # thickness at the cell faces
    F = -beta*hf**3*np.diff(he)/dxe                      # face flux q = −β h³ ∂h/∂x [m²/s]
    he[1:-1] -= dte/dxe*np.diff(F)                       # conservative update: flux in minus flux out
t_end = 2000*dte                                         # the time reached [s]
ref = ch08.thin_film_spread(np.where(np.abs(xe) < 0.01, 0.01, 0.0), xe, [t_end], rho=1260.0, mu=1.0)["h"][-1]
print(f"t = {t_end:.3f} s: volume by hand {he.sum()*dxe:.8e} m²; max |explicit − library| = {np.max(np.abs(he - ref)):.1e} m")
assert abs(he.sum()*dxe - vol0) < 1e-12*vol0             # the explicit scheme conserves the volume to round-off
assert np.max(np.abs(he - ref)) < 0.01*0.01              # agree to 1 % of the bead height (the library adds a 1 µm precursor film)
""")
remind([
    ("animate and show_animation", "`animate(update, frames, fig)` calls update(i) for every frame; `show_animation(anim, player=\"frames\")` gives ◀ ▮▮ ▶ buttons to stop at any time (Ch. 1 P16)."),
])
nb.animation(r"""
xN = ch08.viscous_current_similarity(0.0, run["t"][1:], 1e-4, rho=1260.0, mu=1.0, return_front=True)[1]   # similarity fronts [m]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.4), gridspec_kw=dict(width_ratios=[1.6, 1]))
fill = [a1.fill_between(100*x, 1000*run["h"][0], color=COLORS["blue"], alpha=0.5)]   # the bead (a list so update can replace it)
front, = a1.plot([], [], "|", color=COLORS["amber"], ms=18, mew=2)    # the front marker
a1.set_xlim(-15, 15); a1.set_ylim(0, 10.5); a1.set_xlabel("x [cm]"); a1.set_ylabel("h [mm]")
a2.loglog(run["t"][1:], 100*xN, color=COLORS["accent"], ls="--", lw=1.5, label="slope 1/5 (similarity)")
dots, = a2.loglog([], [], "o", color=COLORS["teal"], ms=4, label="numerical front")
a2.set_xlim(0.08, 1300); a2.set_ylim(0.8, 20); a2.set_xlabel("t [s]"); a2.set_ylabel("front x_N [cm]"); a2.legend(fontsize=7)


def update(i):                                            # frame i = output time run["t"][i]
    fill[0].remove()                                      # replace the filled bead
    fill[0] = a1.fill_between(100*x, 1000*run["h"][i], color=COLORS["blue"], alpha=0.5)
    front.set_data([100*run["x_front"][i]], [0.3])        # front marker on the plate
    dots.set_data(run["t"][1:i+1], 100*run["x_front"][1:i+1])   # the front history so far
    a1.set_title(f"t = {run['t'][i]:.3g} s · area ∫h dx = {run['volume'][i]*1e4:.4f} cm² (never changes)", fontsize=9)
    return []


show_animation(animate(update, frames=len(run["t"]), fig=fig, interval=200), player="frames")
""")
see_read_change(r"""A box of glycerol slumps into a dome and keeps spreading ever more slowly; the amber front marker creeps outward; on
the right the numerical front (teal dots) settles onto the purple slope-1/5 line; the volume in the title never changes — the
`A3` animation, our remake of Fig. 8.11.""",
                r"""On the log–log panel ten times longer buys only $10^{1/5}=1.58$ times more spread; use the ◀ ▶ buttons to stop at any
time and compare the dome with the box.""",
                r"""…the fluid were ten times more viscous: every time is ten times longer (t enters only as ρgt/μ), the curves are the
same.""")
note("N37 [C] · A similarity solution", r"""
— a change of variables that turns this PDE into an ODE — exists for the spreading bead; §8.4 develops the idea (C09) and
delivers this one (C10: $h=At^{-1/5}F(x/Dt^{1/5})$).
""")
explainer("viscous_gravity_current", "Why does a honey drop spread like t^(1/5)?", r"""
Play the spreading from any initial shape and watch it forget that shape and lock onto one profile, while the front marches
along a slope-1/5 line on log–log axes and the volume readout never moves — convergence to self-similarity is a process in
time, and C10 derives the exponent you measure here.""",
          ["Start from 'two humps': they merge, then the profile turns into the same dome as the box.",
           "Watch the log–log front view: compare the slope-1/5 and slope-1/2 guides.",
           "Make the fluid 1000× more viscous (lava lobe) and predict how the curves change before you press play.",
           "Switch the rescaled view to the wrong exponent m = 1/2 and see the collapse fail."])
whatif(r"""
…a flow had no length or time scale at all — a plate started suddenly in a fluid that fills all space? Then the only way to
make y dimensionless is to combine it with time, $y/\sqrt{\nu t}$, and the whole PDE collapses to an ODE (C09).
""")

# =====================================================================================================================
# A.4 §8.4 — R14, C09 (N38–N56, D16–D20, P193–P196, A1, IF5, E5), C10 (N57–N64, N106, D21–D23, P197, A5, E6)
# =====================================================================================================================
nb.section("8.4", "Similarity Solutions for Unsteady Incompressible Viscous Flow", intro=r"""
**What is this section about?** A plate starts moving under still fluid: how far does the motion reach after a time t?
Because the problem has no length or time scale of its own, y and t can only appear together as $y/\sqrt{\nu t}$; the PDE
collapses to an ODE whose solution is the error function, and every profile at every time is one curve. We then turn this
into a recipe — guess a power-law form, demand that every term scales the same way with t — and use it for a thickening shear
layer, a decaying vortex and the spreading bead of C08. The $\sqrt{\nu t}$ scale is the thickness of every laminar boundary
layer (Ch. 9) and of the Ekman layers of Ch. 13.
""")
nb.recap("R14", "Couette start-up: a flow with a length scale", r"""
Ch. 1 solved the start-up of plane Couette flow — a plate suddenly moving at U with a second, fixed plate at distance h — as
a Fourier series (`couette_startup_profile`, in which the plate at y = h moves). Impulsively started parallel flows keep
$u\,\partial u/\partial x=0$, so they are exact solutions (Exercise 8.31). The gap h is an imposed length: the profiles at
different times are *not* one curve. Remove the second plate and that changes (C09).
""", where="Ch. 1 §1.5")
core("C09", "Stokes' first problem", r"""
What is u(y, t) above a plate that suddenly starts moving — and why do all the profiles look the same when rescaled?
""", eqs=("8.30",))
problem(r"""
Yank a plate sideways under a still bath (or start the wind over a calm lake): at first only a thin layer moves; the "news"
that the wall moves spreads upward by viscous diffusion. How thick is the moving layer after one second, one minute, one
hour? This is the prototype of every boundary layer.
""")
idea(r"""
t = 1 s        t = 4 s          t = 16 s        rescaled: y/√(νt)
│▏             │▎                │▍               │ one curve F(η)
│▎ thin        │▌ 2× thicker     │█ 4× thicker    │ for every t, U, ν
══►U           ══►U              ══►U
""", words=r"No ruler in the problem except $\sqrt{\nu t}$, so the profile can only stretch.")
note("N38 [B], N39 [B] · Set-up (Fig. 8.12; often called Rayleigh's problem).", r"""
An infinite plate at y = 0 starts moving at U at t = 0 under fluid at rest in y > 0; nothing depends on x, so continuity gives
v = 0, and the pressure is the same everywhere because far from the plate the fluid is still at rest. What remains is the 1-D
diffusion equation (8.20) — the heat equation of Ch. 1 with ν in place of κ.
""", equation=EQ["8.20"], ref="8.20")
D("D16", ref="8.20")
note("N40 [B], N41 [B], N42 [B] · Conditions.", r"""
$u(y,t=0)=0$ *(8.21)*, the plate condition (8.22) below, and $u(y\to\infty,t)=0$ *(8.23)*. The problem is well posed: (8.20)
has one time derivative and two space derivatives, so it needs one condition in t and two in y (the order of a PDE in each
variable counts the conditions it needs in that variable).
""", equation=EQ["8.22"], ref="8.22")
note("N43 [B], N44 [B] · Dimensional analysis, then linearity.", r"""
**N43** u depends on U, y, t and ν; with two dimensions (L, T) there are three groups, (8.24). **N44** Linearity removes the
second group, (8.25) — one variable instead of two.
""", equation=EQ["8.24"] + r"\ \ \text{(8.24)},\qquad " + EQ["8.25"], ref="8.25")
D("D17", ref="8.25")
nb.code(r"""
groups = ch01.pi_groups({"u": "m/s", "U": "m/s", "y": "m", "t": "s", "nu": "m^2/s"})   # Ch. 1's Π-group finder
for g_ in groups:                                        # each group as exponents of the variables
    print(" · ".join(f"{k}^{v}" for k, v in g_.items() if v != 0))
""", explain=r"""
`pi_groups` returns three independent groups, as D17 step 3 counts: $u/U$, $tU/y$ and $\nu/(Uy)$. They are equivalent to the
book's choice: $tU/y=(y/Ut)^{-1}$, and $\frac{\nu}{Uy}\cdot\frac{Ut}{y}=\frac{\nu t}{y^2}=\big(\frac{y}{\sqrt{\nu t}}\big)^{-2}$ —
any set of three independent combinations is as good as another; linearity then removes the one that contains U.
""")
remind([
    ("chain rule", "if u depends on t only through η(y, t), then ∂u/∂t = (du/dη)(∂η/∂t) (Ch. 1 P49)."),
    ("exponent rules", "$\\eta=y\\nu^{-1/2}t^{-1/2}$, so $\\partial\\eta/\\partial t=-\\tfrac12y\\nu^{-1/2}t^{-3/2}=-\\eta/(2t)$ (Ch. 1 P43)."),
])
note("N45 [B], N46 [B], N47 [B], N48 [B] · Into the equation.", r"""
The chain rule gives $\frac{\partial u}{\partial t}=-\frac{U\eta}{2t}\frac{dF}{d\eta}$ and $\frac{\partial^2u}{\partial y^2}=
\frac{U}{\nu t}\frac{d^2F}{d\eta^2}$; they turn (8.20) into the ODE (8.26) with $F(\eta=0)=1$ *(8.27)* and
$F(\eta\to\infty)=0$ *(8.28)*: the initial condition (8.21) and the far condition (8.23) both become η → ∞, so three
conditions collapse into two — the test that the similarity form is right.
""", equation=EQ["8.26"], ref="8.26")
D("D18", ref="8.26")
P("P194", "Gaussian integral", r"""
The bell curve $e^{-\zeta^2}$ has total area √π: $\int_{-\infty}^{\infty}e^{-\zeta^2}d\zeta=\sqrt\pi$, and by symmetry half of it,
√π/2, lies on each side. No elementary antiderivative exists, which is why its running integral gets a name (next primer).
""", code=r"""
import numpy as np
from scipy.integrate import quad
print(quad(lambda z: np.exp(-z**2), 0, np.inf)[0], np.sqrt(np.pi)/2)   # 0.886227 0.886227
""")
P("P195", "error function erf and its inverses", r"""
$\mathrm{erf}(\zeta)=\frac2{\sqrt\pi}\int_0^\zeta e^{-\xi^2}d\xi$ rises from 0 to 1; erfc = 1 − erf falls from 1 to 0 and is
computed directly without cancellation. Their inverses answer "where does the profile reach this level?":
`scipy.special.erfinv`, `erfcinv`.
""", code=r"""
from scipy import special
print(special.erf(1.0), special.erfc(1.0))               # 0.8427 0.1573
print(2*special.erfcinv(0.01))                           # 3.643: where erfc(η/2) = 1 %
""")
remind([
    ("natural logarithm, exponentiating", "∫dG/G = ln G, and $e^{\\ln G}=G$: exponentiating turns a constant of integration into a factor (Ch. 1 P36)."),
    ("substitution ξ = 2ζ in an integral", "with ξ = 2ζ, dξ = 2dζ and the limits change with it (Ch. 3 P106)."),
    ("complementary error function erfc", "erfc = 1 − erf, computed without cancellation for large arguments (Ch. 4 P123)."),
])
note("N49 [B], N50 [B], N51 [B] · Solving the ODE.", r"""
Separating gives $\frac{dF}{d\eta}=A\exp(-\eta^2/4)$, a second integration gives (8.29), and the conditions give B = 1 and, with
ξ = 2ζ and the Gaussian integral, $A=-1/\sqrt\pi$.
""", equation=EQ["8.29"], ref="8.29")
D("D19", ref="8.30")
P("P196", "scipy.integrate.solve_bvp", r"""
`solve_bvp` solves an ODE with conditions at two ends (here F(0) = 1 and F(η_max) = 0 on a truncated domain) by collocation on
a mesh; we give it the system as first-order equations (F′ = G, G′ = −ηG/2) and a guess.
""", code=r"""
import numpy as np
from scipy.integrate import solve_bvp
eta = np.linspace(0, 12, 50)                             # the truncated η-domain and initial mesh
sol = solve_bvp(lambda e, Y: np.vstack([Y[1], -e/2*Y[1]]),   # Y = (F, F′): F′ = G, G′ = −(η/2)G  (8.26)
                lambda a, b: np.array([a[0] - 1, b[0]]),     # F(0) = 1, F(12) = 0
                eta, np.vstack([np.exp(-eta), -np.exp(-eta)]))   # a rough guess
print(sol.sol(2.0)[0])                                   # 0.1573 = erfc(1)
""")
note("N52 [B], N55 [B] · One curve (Figs. 8.12–8.13, remade below).", r"""
**N52** Profiles at several times (raw) and the same against the similarity variable. ⚠️ The book's figure axis is
$y/(2\sqrt{\nu t})=\eta/2$, while (8.25) defines $\eta=y/\sqrt{\nu t}$; our axes say which. **N55** *Self-similarity:*
profiles at any t, for any U and ν, fall on one curve when u/U is plotted against η — we measure the spread below
(< 10⁻¹²).
""")
D("D20", ref="8.31")
note("N54 [B] · The 99 % thickness", r"""
grows like $t^{1/2}$: in water 3.64 cm after 100 s, 21.9 cm after an hour — diffusion is slow over large distances (the code
below prints $2\,\mathrm{erfc}^{-1}(0.01)=3.643$).
""", equation=EQ["8.31"], ref="8.31")
nb.worked_example("a plate in water after 100 s", r"""
ν = 10⁻⁶ m²/s, t = 100 s, U = 0.1 m/s.
1. $\sqrt{\nu t}=\sqrt{10^{-4}}=0.01$ m.
2. At y = 1 cm: η = 1, η/2 = 0.5, u/U = erfc(0.5) = 0.4795 — half the plate speed.
3. $\delta_{99}=3.643\times0.01=3.64$ cm.
4. After 400 s: $\sqrt{\nu t}$ doubles to 2 cm, $\delta_{99}=7.29$ cm, and u(2 cm) = 0.4795U again (same η).
5. Wall stress $\tau_w=\mu U/\sqrt{\pi\nu t}=10^{-3}\times0.1/\sqrt{\pi\times10^{-4}}=5.64\times10^{-3}$ Pa, falling like
   $t^{-1/2}$.
""")
nb.code(r"""
t = 100.0; y = np.array([0.0, 0.005, 0.01, 0.02, 0.0364])   # a time [s] and five heights [m]
print("u =", ch08.stokes_first_problem(y, t, U=0.1, nu=NU_W), "m/s")   # (8.30) via erfc
print("η =", ch08.similarity_variable(y, t, NU_W), " η/2 =", ch08.similarity_variable(y, t, NU_W, half=True))
print(f"δ99 = {ch08.diffusion_thickness(t, NU_W):.5f} m at 100 s, {ch08.diffusion_thickness(3600.0, NU_W):.4f} m after an hour")
print("2 erfcinv(0.01) =", round(2*special.erfcinv(0.01), 4))   # the 3.643 of (8.31)
st = ch08.stokes_first_state(t, U=0.1, nu=NU_W)          # the numbers the explainer shows
print({k: round(float(st[k]), 6) for k in ("sqrt_nut", "eta_edge", "delta", "tau_w", "vorticity_content")})
print("∫ω dy =", ch08.vorticity_content(t, U=0.1, nu=NU_W), "m/s (= U)")   # N53 below, by quad
print("Ch. 4:", ch04.stokes_first_problem(y, t, 0.1, NU_W))   # the same u from Ch. 4's earlier function
eta, F = ch08.similarity_ode_solve("stokes1")             # solve_bvp on (8.26)–(8.28), knowing nothing about erf
print(f"max |F_bvp − erfc(η/2)| = {np.max(np.abs(F - special.erfc(eta/2))):.1e}")
""", explain=r"""
1. `stokes_first_problem` is (8.30) computed with erfc: u = [0.1, 0.0724, 0.0480, 0.0157, 0.0010] m/s.
2. `similarity_variable` returns η = y/√(νt), and with `half=True` the book's figure axis η/2.
3. `diffusion_thickness` is (8.31): 3.64 cm after 100 s, 21.9 cm after an hour.
4. `stokes_first_state` bundles √(νt) = 1 cm, the edge η₉₉ = 3.643, δ₉₉, the wall stress 5.64 × 10⁻³ Pa and the vorticity
   content; `vorticity_content` integrates ω over the whole layer and gets U (N53).
5. Ch. 4's earlier function gives the same numbers; the boundary-value solver, knowing nothing about erf, lands on the same
   curve to about 10⁻¹³.
""")
P("P193", "Crank–Nicolson with scipy.linalg.solve_banded", r"""
Crank–Nicolson averages the diffusion term between the old and the new time level: second-order accurate in time and stable
for any step. Each step solves a tridiagonal system, which `scipy.linalg.solve_banded` does in O(N) from the three diagonals
stored as rows. An impulsive start (the plate jumps to U) makes pure Crank–Nicolson ring, so we take the first two steps
with backward Euler (fully implicit), which damps the ringing.
""", code=r"""
import numpy as np
from scipy.linalg import solve_banded
ab = np.array([[0, -1, -1], [4, 4, 4], [-1, -1, 0]], float)   # upper, main, lower diagonals (as rows)
print(solve_banded((1, 1), ab, np.array([3.0, 2.0, 3.0])))    # [1. 1. 1.]
""")
remind([
    ("observed order of convergence", "if the error behaves like C·Δyᵖ, the slope of log(error) against log(Δy) is p; `observed_order` fits it (Ch. 1 P13, log–log slopes)."),
])
nb.check_agree(r"""
U0, t_end = 0.1, 1000.0                                  # plate speed [m/s] and final time [s]


def cn_by_hand(N, nt):                                   # Crank–Nicolson march of (8.20) on N points, nt steps
    yy = np.linspace(0, 8*np.sqrt(NU_W*t_end), N)        # up to 8√(νt_end), where u ≈ 0 [m]
    dy, dt = yy[1] - yy[0], t_end/nt                     # grid step [m], time step [s]
    r = NU_W*dt/dy**2                                    # the diffusion number ν Δt/Δy²
    u = np.zeros(N)                                      # fluid at rest (8.21)
    for k in range(nt):
        th = 1.0 if k < 2 else 0.5                       # two backward-Euler steps, then Crank–Nicolson (θ = ½)
        ab = np.zeros((3, N - 2))                        # the three diagonals of the implicit part
        ab[0, 1:], ab[1, :], ab[2, :-1] = -th*r, 1 + 2*th*r, -th*r
        rhs = u[1:-1] + (1 - th)*r*(u[2:] - 2*u[1:-1] + u[:-2])   # the explicit part of the average
        rhs[0] += th*r*U0                                # the wall value u(0) = U (8.22) moves to the right side
        u = np.r_[U0, linalg.solve_banded((1, 1), ab, rhs), 0.0]   # new profile with the two boundary values
    return yy, u


N_ = 400 if not FAST else 200                            # grid points
yy, u_cn = cn_by_hand(N_, 200 if not FAST else 100)       # our march
exact = ch08.stokes_first_problem(yy, t_end, U0, NU_W)   # (8.30)
print(f"hand CN vs (8.30): max error {np.max(np.abs(u_cn - exact)):.1e} m/s")
assert np.allclose(u_cn, exact, rtol=1e-3, atol=2e-4)    # the march lands on the erfc curve
u_lib = ch08.crank_nicolson_1d(np.zeros(N_), yy, t_end/(200 if not FAST else 100), 200 if not FAST else 100, NU_W, U0, 0.0)
assert np.allclose(u_lib, u_cn, rtol=0, atol=1e-12)      # the library solver does exactly this
dyf = yy[1] - yy[0]; dtf = diffusion.stable_time_step(NU_W, dyf); nf = int(np.ceil(t_end/dtf))   # Ch. 1's explicit FTCS
u_ftcs = diffusion.ftcs_diffusion_1d(np.r_[U0, np.zeros(N_ - 1)], NU_W, dyf, t_end/nf, nf, values=(U0, 0.0), save_every=nf)[-1]
print(f"FTCS ({nf} small steps) vs (8.30): max error {np.max(np.abs(u_ftcs - exact)):.1e} m/s")
assert np.allclose(u_ftcs, exact, atol=1e-3*U0)          # the explicit scheme agrees too, with ~28× more steps
errs, steps = [], []                                     # a convergence table: refine space and time together
for N2 in (100, 200, 400):
    y2, u2 = cn_by_hand(N2, N2)
    errs.append(np.max(np.abs(u2 - ch08.stokes_first_problem(y2, t_end, U0, NU_W)))); steps.append(y2[1] - y2[0])
    print(f"N = {N2}: max error {errs[-1]:.2e} m/s")
print(f"observed order ≈ {observed_order(steps, errs):.2f}")   # ≈ 2: second order in space and time
""")
nb.figure(r"""
U0 = 0.1                                                  # plate speed [m/s]
ts = [10, 30, 100, 300, 1000]                             # five times [s]
shades = plt.cm.Blues(np.linspace(0.4, 1.0, len(ts)))     # darker blue for later times
y = np.linspace(0, 0.14, 400)                             # heights [m]
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12, 3.8))
for t_, c in zip(ts, shades):
    u = ch08.stokes_first_problem(y, t_, U0, NU_W)        # (8.30)
    a1.plot(u/U0, 100*y, color=c, lw=2, label=f"t = {t_} s")
    d99 = ch08.diffusion_thickness(t_, NU_W)              # (8.31)
    a1.plot([0, 0.08], [100*d99]*2, color=COLORS["amber"], lw=2)   # δ99 tick
    a2.plot(u/U0, ch08.similarity_variable(y, t_, NU_W), color=c, lw=2)   # the same profile against η
a2.plot(u_cn[::8]/U0, ch08.similarity_variable(yy[::8], 1000.0, NU_W), "o", color=COLORS["teal"], ms=3, label="Crank–Nicolson, t = 1000 s")
a1.set_xlabel("u/U [–]"); a1.set_ylabel("y [cm]"); a1.legend(fontsize=7); a1.set_title("Raw profiles (amber: δ99)", fontsize=10)
a2.set_xlabel("u/U [–]"); a2.set_ylabel("η = y/√(νt) [–]"); a2.set_ylim(0, 6); a2.legend(fontsize=7)
a2.axhline(3.643, color=COLORS["amber"], lw=0.8); a2.text(0.5, 3.75, "η99 = 3.643", color=COLORS["amber"], fontsize=8)
sec = a2.secondary_yaxis("right", functions=(lambda e: e/2, lambda e: 2*e)); sec.set_ylabel("η/2 = y/(2√(νt)) (book's axis)")
a2.set_title("Against η: one curve", fontsize=10)
tt = np.geomspace(1, 1e4, 50)                             # times [s]
a3.loglog(tt, 100*ch08.diffusion_thickness(tt, NU_W), color=COLORS["blue"], lw=2, label="water")
a3.loglog(tt, 100*ch08.diffusion_thickness(tt, nu_air), color=COLORS["teal"], lw=2, label="air")
a3.set_xlabel("t [s]"); a3.set_ylabel("δ99 [cm]"); a3.legend(fontsize=8); a3.set_title("δ99 ∝ t^½", fontsize=10)
fig.suptitle("One curve for all times: the profile only stretches like √(νt)", fontsize=11)
savefig(fig, "ch08", "c09_stokes_first"); plt.show()
""", see=r"""(a) five profiles of growing thickness with their 99 % heights (amber ticks); (b) the same five profiles against
η — they coincide exactly, and the Crank–Nicolson dots sit on them; (c) two straight lines of slope ½ — our remake of
Figs. 8.12–8.13 `N52`.""",
    read=r"""(b): the 1 % point sits at η = 3.64 (η/2 = 1.82 on the book's axis, the right-hand scale); (c): 100 times longer gives
10 times thicker.""",
    change=r"""…the fluid were air: every thickness grows √15 ≈ 3.9 times faster, but (b) is unchanged.""")
nb.code(r"""
eta_grid = np.linspace(0, 8, 161)                        # a common η grid
spread = 0.0                                             # the largest difference from erfc(η/2)
for U_, nu_, t_ in ((0.1, 1e-6, 10.0), (2.0, 1.5e-5, 3.0), (0.01, 1e-3, 500.0)):   # water, air, syrup at three times
    y_ = eta_grid*np.sqrt(nu_*t_)                        # heights with those η [m]
    spread = max(spread, np.max(np.abs(ch08.stokes_first_problem(y_, t_, U_, nu_)/U_ - special.erfc(eta_grid/2))))
print(f"collapse: max |u/U − erfc(η/2)| over three (U, ν, t) = {spread:.1e}")
assert spread < 1e-12                                    # N55: one curve for every U, ν and t
""")
note("N53 [B] · Reading it as vorticity.", r"""
At t = 0 the wall creates a vortex sheet (all the vorticity in an infinitely thin layer); afterwards it diffuses out:
$\omega=-\partial u/\partial y=\frac{U}{\sqrt{\pi\nu t}}e^{-y^2/4\nu t}>0$. The total below never changes, so no new vorticity is
made after t = 0. ⚠️ The book prints this integral as −U; since u decreases upward, ω is positive and the integral is +U
(checked above by `vorticity_content`). The heat analogue: a cold solid whose face is suddenly heated.
""", equation=r"\int_0^\infty\omega\,dy=u(0)-u(\infty)=U")
note("N56 [B] · Similarity needs the absence of imposed scales.", r"""
Stop the plate at t = T (Exercise 8.30) and a time scale appears; by superposition (ours) the velocity is
$u=U\big[\mathrm{erfc}\frac{y}{2\sqrt{\nu t}}-\mathrm{erfc}\frac{y}{2\sqrt{\nu(t-T)}}\big]$ for t > T. Put a second plate at
y = h (R14) and a length scale appears. Either way the profiles stop collapsing — the spin-down of a stirred tank (Ch. 13) is
of this kind.
""")
nb.figure(r"""
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.4))
eta_ = np.linspace(0, 6, 200)                             # η grid [–]
a1.plot(special.erfc(eta_/2), eta_, color=COLORS["accent"], ls="--", lw=1.5, label="erfc(η/2): plate never stops")
for t_, c in ((50, COLORS["blue"]), (150, COLORS["teal"]), (400, COLORS["orange"])):   # times [s]; the plate stops at T = 100 s
    y_ = eta_*np.sqrt(NU_W*t_)                            # heights at this time [m]
    a1.plot(ch08.stokes_first_stopped(y_, t_, 100.0, U=1.0, nu=NU_W), eta_, color=c, lw=2, label=f"t = {t_} s")
a1.set_xlabel("u/U [–]"); a1.set_ylabel("η = y/√(νt) [–]"); a1.legend(fontsize=7); a1.set_title("Plate stopped at T = 100 s", fontsize=10)
h_ = 0.02                                                 # a second, fixed plate 2 cm away [m]
for t_, c in ((10, COLORS["blue"]), (100, COLORS["teal"]), (400, COLORS["orange"])):
    yb = np.linspace(0, h_, 200)                          # the gap [m]; Ch. 1's function moves the plate at y = h
    ub = diffusion.couette_startup_profile(yb, t_, 0.1, h_, NU_W)   # Fourier series of the start-up [m/s]
    a2.plot(ub/0.1, (h_ - yb)/np.sqrt(NU_W*t_), color=c, lw=2, label=f"t = {t_} s")   # distance from the moving plate, in η
a2.plot(special.erfc(eta_/2), eta_, color=COLORS["accent"], ls="--", lw=1.5, label="erfc(η/2)")
a2.set_ylim(0, 6); a2.set_xlabel("u/U [–]"); a2.set_ylabel("η from the moving plate [–]"); a2.legend(fontsize=7)
a2.set_title("A second wall 2 cm away (Couette start-up)", fontsize=10)
fig.suptitle("An imposed time or length scale breaks the collapse", fontsize=11)
savefig(fig, "ch08", "c09_no_collapse"); plt.show()
""", see=r"""(a) after the plate stops (T = 100 s) the curves no longer lie on the erfc line and the wall speed drops to zero;
(b) the Couette start-up curves follow erfc while the layer is thin (10 s) and leave it once it feels the far wall (400 s).""",
    read=r"""In (b) η₉₉ = 3.64 corresponds to δ₉₉ = 3.64√(νt): 1.2 cm at 10 s (well inside the 2 cm gap), 7.3 cm at 400 s (far
beyond it) — the far wall's no-slip condition bends the curve back to the straight Couette line.""",
    change=r"""…T → ∞ (the plate never stops): the stopped-plate curves rejoin the erfc curve.""")
nb.animation(r"""
nfr = 60 if not FAST else 30                              # frames
tf = np.geomspace(1, 1000, nfr)                           # a log clock t = 1 … 1000 s
U0 = 0.1                                                  # plate speed [m/s]
yv = np.linspace(0, 0.12, 121)                            # heights [m]
sgrid = np.geomspace(1e-4, 1000, 3000)                    # fine times for the dye displacement [s]
X = integrate.cumulative_trapezoid(ch08.stokes_first_problem(yv[None, :], sgrid[:, None], U0, NU_W), sgrid, axis=0, initial=0)
disp = np.array([X[np.searchsorted(sgrid, t_)] / (U0*t_) for t_ in tf])   # dye displacement ÷ plate travel U t [–]
tcn, Ucn = ch08.crank_nicolson_1d(np.zeros(121), yv, 1000.0/400, 400, NU_W, U0, 0.0, return_all=True)   # numerical twin
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
dye, = a1.plot(disp[0], 100*yv, color=COLORS["blue"], lw=2)   # the dyed line, bent into the time-integrated profile
a1.axhline(0, color=COLORS["ink"], lw=4)                  # the plate
a1.set_xlim(0, 1.05); a1.set_xlabel("dye displacement ÷ plate travel Ut [–]"); a1.set_ylabel("y [cm]")
prof, = a2.plot([], [], color=COLORS["blue"], lw=2, label="(8.30)")
dots, = a2.plot([], [], "o", color=COLORS["teal"], ms=3, label="Crank–Nicolson")
mark = a2.axhline(0, color=COLORS["amber"], lw=1.5, label="δ99 = 3.64√(νt)")
a2.set_xlim(0, 1.05); a2.set_ylim(0, 12); a2.set_xlabel("u/U [–]"); a2.legend(fontsize=7, loc="upper right")


def update(i):                                            # frame i: time tf[i]
    t_ = tf[i]
    dye.set_xdata(disp[i])                                # the dye line
    prof.set_data(ch08.stokes_first_problem(yv, t_, U0, NU_W)/U0, 100*yv)   # the exact profile
    k = np.argmin(np.abs(tcn - t_))                       # the nearest Crank–Nicolson output
    dots.set_data(Ucn[k][::4]/U0, 100*yv[::4])
    mark.set_ydata([100*ch08.diffusion_thickness(t_, NU_W)]*2)   # the rising δ99 marker
    a1.set_title(f"t = {t_:6.1f} s", fontsize=10)
    return []


show_animation(animate(update, frames=nfr, fig=fig, interval=80))
""")
see_read_change(r"""On the left a dyed vertical line above the plate tilts and bends as the fluid is dragged along (drawn as a fraction of
the plate's own travel); on the right the profile grows, the Crank–Nicolson dots ride the erfc curve, and the amber δ₉₉ marker
rises ever more slowly — the `A1` animation.""",
                r"""The marker height at time t is 3.64√(νt): it needs four times as long to double; the clock is logarithmic, so equal
steps of the movie are equal factors of time.""",
                r"""…ν were ten times larger: the same movie, ten times faster.""")
remind([
    ("animate_figure (plotly time slider)", "`animate_figure(frame_fn, times, …)` is the time-player version of `slider_figure`: Play/Pause and a time slider, every frame precomputed, so it works on the web page (Ch. 1 P17)."),
])
nb.plotly(r"""
yv = np.linspace(0, 0.08, 120)                            # heights [m]


def frame(t_):                                            # curves at time t_ [s]
    u = ch08.stokes_first_problem(yv, t_, 0.1, NU_W)/0.1  # u/U
    return {"raw: u/U against y [cm]": (u, 100*yv),
            "rescaled: u/U against η = y/√(νt)": (u, ch08.similarity_variable(yv, t_, NU_W))}


fig = animate_figure(frame, np.geomspace(1, 1000, 25 if not FAST else 12), time_label="t", unit="s",
                     xlabel="u/U [–]", ylabel="y [cm] (blue) or η [–] (purple)", title="Raw profiles spread; rescaled ones do not",
                     yrange=[0, 8], xrange=[0, 1])
recolor(fig, {"raw: u/U against y [cm]": COLORS["blue"], "rescaled: u/U against η = y/√(νt)": COLORS["accent"]},
        {"rescaled: u/U against η = y/√(νt)": "dash"})
fig.show()
""", explain=r"""
Press Play: the blue profile (plotted against y in centimetres) thickens as √t, while the purple one (the same profile plotted
against η) never moves — the collapse, without a kernel.
""")
explainer("stokes_first_problem", "How can profiles at all times be one curve?", r"""
Start the plate, watch the profile spread with a moving δ₉₉ marker, then flip to the rescaled mode: every earlier profile falls
on one curve, and changing ν or U does not break it — while stopping the plate or adding a second wall does.""",
          ["Play with water, then switch to air: the layer grows faster but the rescaled curve does not move.",
           "Click a point in the profile view: read η and the erfc arithmetic of u/U.",
           "Choose 'stop the plate at T' and watch the collapse fail after T.",
           "Open the Derivation tab at D19 (the substitution ξ = 2ζ): the factor 2 in y/(2√(νt)) appears."])
whatif(r"""
…the initial state were a velocity jump inside the fluid (two streams), a line vortex, or the spreading bead of C08? The same
trick works, but the similarity form may need a power of t in front — and finding the exponents becomes the whole game (C10).
""")

core("C10", "The similarity ansatz and exponent matching", r"""
How do you find the exponents that make a problem self-similar?
""", eqs=("8.32a",))
problem(r"""
Two streams slide past each other and the shear layer between them thickens; a vortex left alone decays; a honey bead spreads.
Each looks self-similar, but with different growth laws ($t^{1/2}$, $t^{1/5}$) and amplitudes that fall with time. We want a
procedure that finds those laws instead of guessing them.
""")
idea("", table=r"""
| step | move |
|---|---|
| 1 | guess $\gamma=At^{-n}F(\xi/\delta(t))$ (or $A\xi^{-n}F(\xi/\delta)$ if ξ appears in the initial condition) |
| 2 | substitute into the PDE, divide by one coefficient |
| 3 | demand every bracketed coefficient is the same power of t → one relation between n and δ(t) |
| 4 | a conserved quantity (momentum jump, circulation, volume) → the second relation |
""", words=r"**Two unknown exponents, two conditions: nothing is left to choose.**")
P("P197", "exponent matching for similarity forms", r"""
If $c_1t^aF(\eta)+c_2t^bG(\eta)=0$ must hold for every t and every η, the powers of t must agree (a = b) — otherwise dividing
by $t^a$ leaves a t that no function of η can cancel. Power laws $\delta=Dt^m$ turn every bracket into a power of t, and
matching powers gives linear equations for the exponents.
""", code=r"""
import sympy as sp
n, m = sp.symbols('n m')
print(sp.solve([sp.Eq(-n - 1, -4*n - 2*m), sp.Eq(-n + m, 0)], [n, m]))   # {n: 1/5, m: 1/5}
""")
note("N57 [B] · The second form (8.32b)", r"""
puts a power of the space coordinate in front; use it when ξ appears in the initial or boundary condition (the line vortex
below, Γ/2πr).
""", equation=EQ["8.32b"], ref="8.32b")
remind([
    ("separable first-order ODE for δ(t)", "δ dδ/dt = C: move all δ to one side and all t to the other, then integrate both (Ch. 1 P42)."),
])
D("D21", ref="8.32a")
note("N58 [B] · Example 8.4 in one line.", r"""
It recovers Stokes' first problem from the ansatz: A = 1, n = 0, δδ′ = C₁ν, $\delta=\sqrt{2C_1\nu t}$; C₁ = ½ gives the η of
(8.25).
""")
remind([
    ("product rule read backwards, (ηF)′", "$F+\\eta\\frac{dF}{d\\eta}=\\frac{d}{d\\eta}(\\eta F)$ — the product rule read from right to left (Ch. 1 P38)."),
])
D("D22", check_src=r"""
vs = ch08.similarity_reduce_sympy("vortex_sheet")        # D22 in sympy: the ansatz in ∂ω/∂t = ν ∂²ω/∂y²
print("brackets:", vs["brackets"])                       # n/t, δ'/δ, ν/δ² with δ = √(νt): all ∝ 1/t (step 6)
print("ODE:", vs["ode"])                                 # −nF − ½ηF′ = F″
print("n from the conserved jump:", vs["n"])             # 1/2 (step 9)
assert sp.simplify(vs["residual"]) == 0                  # the final ω solves the reduced problem
y_, t_, U_, nu_ = sp.symbols("y t U nu", positive=True)  # fresh symbols for a direct check
omega = -U_/sp.sqrt(sp.pi*nu_*t_)*sp.exp(-y_**2/(4*nu_*t_))   # step 14's vorticity
assert sp.simplify(sp.diff(omega, t_) - nu_*sp.diff(omega, y_, 2)) == 0   # it solves the diffusion equation
print("−∫ω dy =", sp.simplify(-sp.integrate(omega, (y_, -sp.oo, sp.oo))))   # 2U at every t: the jump is conserved
""")
note("N59 [B] · Example 8.5, the viscous vortex sheet.", r"""
The solution is below. Define the layer's width by u = ±0.95U: η = ±2 erfinv(0.95) = ±2.772 and the width is 5.544√(νt).
⚠️ The book prints ±2.76 but also the width 5.54 — the 2.76 is a rounding slip. ⚠️ Ch. 5's `diffusing_vortex_sheet(y, t, gamma,
nu)` takes γ = u_below − u_above = −2U.
""", equation=r"\omega_z(y,t)=-\frac{U}{\sqrt{\pi\nu t}}\exp\Big\{-\frac{y^2}{4\nu t}\Big\},\qquad u(y,t)=U\,\mathrm{erf}\Big\{\frac{y}{2\sqrt{\nu t}}\Big\}")
nb.code(r"""
u, w = ch08.vortex_sheet_diffusion(0.001, 1.0, U=0.01, nu=NU_W)   # at y = 1 mm after 1 s, U = 1 cm/s
print(f"u = {u:.6f} m/s, ω_z = {w:.4f} 1/s;  Ch. 5 (γ = −2U):", ch05.diffusing_vortex_sheet(0.001, 1.0, -0.02, NU_W))
print(f"width = {ch08.transition_width(1.0, NU_W):.4e} m, 2 erfinv(0.95) = {2*special.erfinv(0.95):.4f}")
""")
remind([
    ("Galilean shift, temporally developing boundary layer, C_f", "adding a constant velocity to everything changes no stress or vorticity (Ch. 3 P96, frames of reference); the skin-friction coefficient is $C_f=\\tau_w/(\\tfrac12\\rho U^2)$."),
])
note("N60 [B] · A temporally developing boundary layer.", r"""
Shift the upper half of Example 8.5 by −U and flip the sign: it is Stokes' first problem. Its wall stress and skin friction
are below; with Ut read as a distance x they become $C_f=1.128\,\mathrm{Re}_x^{-1/2}$ — the $\mathrm{Re}_x^{-1/2}$ law of a
laminar boundary layer on a plate (Ch. 9 finds 0.664 for the real, spatially growing layer: same exponent, different constant).
""", equation=r"\tau_w=\mu\Big(\frac{\partial u}{\partial y}\Big)_{y=0}=\frac{\mu U}{\sqrt{\pi\nu t}},\qquad C_f=\frac{\tau_w}{\frac12\rho U^2}=\frac{2}{\sqrt\pi}\sqrt{\frac{\nu}{U^2t}}")
nb.code(r"""
tb = ch08.temporal_bl_wall_stress(1e4*NU_W/1.0**2, 1.0, NU_W)   # U = 1 m/s at the time where Re_x = U²t/ν = 10⁴
print(f"Re_x = {tb['Re_x']:.0f}, C_f = {tb['Cf']:.5f} (1.128/√Re_x); Blasius (Ch. 9) would give {0.664/np.sqrt(1e4):.5f}")
""")
note("N61 [B], N62 [B] · Example 8.6, a decaying line vortex (stated; the same moves as D22 with the form (8.32b), Ar⁻ⁿ = Γ/2πr).", r"""
**N61** The similarity ODE $(\frac1\eta-\frac\eta2)F'=F''$ gives $F=1-e^{-\eta^2/4}$ and the velocity below — the Gaussian
(Lamb–Oseen) vortex $u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)$ *(3.29)* with σ² = 4νt (Ch. 3, Ch. 5): rigid
rotation inside r ≈ 2.24√(νt), an ideal vortex outside. **N62** Exercise 8.26, a vortex switched on:
$u_\theta=\frac{\Gamma}{2\pi r}\exp\{-\frac{r^2}{4\nu t}\}$ — the difference between the steady ideal vortex and the decaying one.
""", equation=r"u_\theta(r,t)=\frac{\Gamma}{2\pi r}\Big[1-\exp\Big\{-\frac{r^2}{4\nu t}\Big\}\Big]")
nb.code(r"""
r_, t_, G_, nu_ = sp.symbols("r t Gamma nu", positive=True)   # symbols for the two line-vortex solutions
for name, u_th in (("decaying (N61)", G_/(2*sp.pi*r_)*(1 - sp.exp(-r_**2/(4*nu_*t_)))),
                   ("switched on (N62)", G_/(2*sp.pi*r_)*sp.exp(-r_**2/(4*nu_*t_)))):
    res = sp.diff(u_th, t_) - nu_*sp.diff(sp.diff(r_*u_th, r_)/r_, r_)   # ∂u/∂t − ν d/dr[(1/r) d(r u)/dr] (P186's azimuthal Laplacian)
    print(name, "residual:", sp.simplify(res))           # 0: both solve the diffusion of u_θ
print(np.allclose(ch08.line_vortex_decay(0.01, 100.0, 1e-3, NU_W), VX.gaussian_vortex(0.01, 1e-3, 2*np.sqrt(NU_W*100.0))[0]))   # = Ch. 3's Gaussian vortex, σ = 2√(νt)
""")
D("D23")
note("N63 [B] · Example 8.7 solved.", r"""
The spreading bead of C08 is self-similar with n = m = 1/5: $h(x,t)=At^{-1/5}F(x/Dt^{1/5})$. Solving the ODE for F (not done in
the book) gives Huppert's (1982) dome $F\propto(1-\eta^2/\eta_N^2)^{1/3}$ and the front $x_N=\eta_N(\beta A^3t)^{1/5}$,
$\eta_N=1.411$ (benchmark V5) — the `viscous_current_similarity` used in C08. ⚠️ The book's last line says "the final equation of
Example 8.2"; it means Example 8.3.
""")
nb.worked_example("a 1 cm/s shear layer after one second", r"""
U = 0.01 m/s on each side, ν = 10⁻⁶ m²/s, t = 1 s.
1. $\sqrt{\nu t}=10^{-3}$ m.
2. Width 5.544 × 1 mm = 5.5 mm.
3. At y = 1 mm: η = 1, u = U erf(0.5) = 0.0052 m/s.
4. Peak vorticity $-U/\sqrt{\pi\nu t}=-0.01/1.772\times10^{-3}=-5.64$ s⁻¹.
5. After 100 s: width 5.5 cm, peak −0.564 s⁻¹, and $-\int\omega\,dy=2U$ still.
""")
nb.code(r"""
for case in ("stokes1", "vortex_sheet", "line_vortex", "spreading"):   # four similarity reductions in sympy
    s = ch08.similarity_reduce_sympy(case)
    print(f"{case:13s} n = {s['n']},  δ(t) = {s['delta']}")
good = ch08.similarity_collapse_error("spreading", 0.2, 0.2, [10, 30, 100, 300, 1000])    # right exponents n = m = 1/5
bad = ch08.similarity_collapse_error("spreading", 0.25, 0.5, [10, 30, 100, 300, 1000])    # a diffusion-like guess m = 1/2
print(f"collapse spread of the bead: {good:.1e} at (n, m) = (1/5, 1/5), {bad:.2f} at (1/4, 1/2)")
""", explain=r"""
1. `similarity_reduce_sympy` repeats Examples 8.4–8.7 symbolically: n = 0 and δ = √(νt) for Stokes' first problem, n = ½ for
   the sheet, n = 1 (a power of r) for the line vortex, and n = m = 1/5 with δ = Dt^(1/5) for the bead.
2. `similarity_collapse_error` samples the exact bead solution at five times, rescales with trial exponents and measures how
   far the curves are apart: round-off at the right pair, order one at a wrong pair.
""")
nb.check_agree(r"""
tw = np.array([1.0, 10.0, 100.0, 1000.0])                # times [s]
widths = [2*optimize.brentq(lambda y_: ch08.vortex_sheet_diffusion(y_, t_, 1.0, NU_W)[0] - 0.95, 1e-6, 1.0) for t_ in tw]   # u = 0.95U
slope_sheet = np.polyfit(np.log(tw), np.log(widths), 1)[0]   # log–log slope of the sheet width
m_ = run["t"] >= 10                                      # the bead of C08 once it has forgotten the box
slope_bead = np.polyfit(np.log(run["t"][m_]), np.log(run["x_front"][m_]), 1)[0]   # slope of the numerical front
print(f"sheet width slope = {slope_sheet:.4f} (expect 1/2); bead front slope = {slope_bead:.4f} (expect 1/5)")
assert abs(slope_sheet - 0.5) < 5e-3                     # D22's n = 1/2 measured
assert abs(slope_bead - 0.2) < 1e-2                      # D23's m = 1/5 measured (the front is found on a grid: ~1 %)
""")
nb.figure(r"""
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12, 3.8))
yv = np.linspace(-0.012, 0.012, 400)                      # across the sheet [m]
for t_, c in ((1.0, COLORS["accent"]), (4.0, COLORS["blue"])):   # two times [s]
    u, w = ch08.vortex_sheet_diffusion(yv, t_, 0.01, NU_W)   # Example 8.5 with U = 1 cm/s
    a1.plot(w, 1e3*yv, color=c, lw=2, label=f"t = {t_:g} s")
    a2.plot(u/0.01, ch08.similarity_variable(yv, t_, NU_W, half=True), color=c, lw=2, ls="-" if t_ == 1 else "--")
a1.set_xlabel("ω_z [1/s]"); a1.set_ylabel("y [mm]"); a1.legend(fontsize=8); a1.set_title("Sheet vorticity: peak halves, width doubles", fontsize=9)
a2.set_ylim(-3, 3); a2.set_xlabel("u/U [–]"); a2.set_ylabel("y/(2√(νt)) (the book's axis)"); a2.set_title("One curve for both times", fontsize=9)
rr = np.linspace(1e-4, 0.03, 300)                          # radius [m]
G0_ = 1e-4                                                 # circulation Γ [m²/s] (ours)
a3.plot(100*rr, 1e3*G0_/(2*np.pi*rr), color=COLORS["muted"], ls="--", lw=1.5, label="νt = 0: ideal vortex")
for nt_, c in ((1e-5, COLORS["blue"]), (4e-5, COLORS["teal"]), (1.6e-4, COLORS["orange"])):   # νt [m²] (our values)
    a3.plot(100*rr, 1e3*ch08.line_vortex_decay(rr, nt_/NU_W, G0_, NU_W), color=c, lw=2, label=f"νt = {nt_:.1e} m²")
    rc = 2.24*np.sqrt(nt_)                                  # the solid-body core edge ≈ 2.24√(νt)
    a3.plot(100*rc, 1e3*ch08.line_vortex_decay(rc, nt_/NU_W, G0_, NU_W), "o", color=c)
a3.set_ylim(0, 10); a3.set_xlabel("r [cm]"); a3.set_ylabel("u_θ [mm/s]"); a3.legend(fontsize=7)
a3.set_title("Decaying line vortex (dots: core edge)", fontsize=9)
fig.suptitle("Three self-similar diffusions", fontsize=11)
savefig(fig, "ch08", "c10_vortex_diffusion"); plt.show()
""", see=r"""(a) the Gaussian vorticity of the thickening sheet: from 1 s to 4 s the peak halves and the width doubles; (b) both
velocity profiles on one curve against the book's axis y/(2√(νt)); (c) the ideal vortex (dashed) smoothed into a solid-body
core that grows with νt — our remake of Figs. 8.14–8.15 `N106` (our values of νt, not the book's).""",
    read=r"""In (c) the dots mark r ≈ 2.24√(νt), where the speed peaks: inside it the fluid turns like a solid, outside it the speed
still falls as Γ/2πr.""",
    change=r"""…Γ doubled: every curve in (c) doubles; its shape and the core radius 2.24√(νt) do not change.""")
nb.animation(r"""
nfr = 60 if not FAST else 30                              # frames
tf = np.geomspace(0.5, 50, nfr)                           # a log clock [s]
yv = np.linspace(-0.02, 0.02, 300); rr = np.linspace(1e-4, 0.03, 300)   # across the sheet, and radius [m]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4))
l1, = a1.plot([], [], color=COLORS["blue"], lw=2, label="u/U")
l2, = a1.plot([], [], color=COLORS["accent"], lw=2, label="ω/ω_peak(0.5 s)")
l3, = a2.plot([], [], color=COLORS["blue"], lw=2)
a1.legend(fontsize=7, loc="upper left"); a1.set_xlim(-1.1, 1.1)
nres = int(0.75*nfr)                                      # the last quarter of the frames switches to rescaled axes
w0 = abs(ch08.vortex_sheet_diffusion(0.0, tf[0], 0.01, NU_W)[1])   # the first frame's peak vorticity [1/s]


def update(i):                                            # frame i
    t_ = tf[i]
    u, w = ch08.vortex_sheet_diffusion(yv, t_, 0.01, NU_W)   # the sheet at t_
    ut = ch08.line_vortex_decay(rr, t_, 1e-4, NU_W)       # the vortex at t_
    if i < nres:                                          # raw axes
        l1.set_data(u/0.01, 1e3*yv); l2.set_data(w/w0, 1e3*yv); a1.set_ylim(-20, 20); a1.set_ylabel("y [mm]")
        l3.set_data(100*rr, 1e3*ut); a2.set_xlim(0, 3); a2.set_ylim(0, 3); a2.set_xlabel("r [cm]"); a2.set_ylabel("u_θ [mm/s]")
        a1.set_title(f"t = {t_:.2f} s: the sheet thickens", fontsize=9); a2.set_title("the vortex decays", fontsize=9)
    else:                                                 # rescaled axes: every profile falls on one curve
        s_ = np.sqrt(NU_W*t_)
        l1.set_data(u/0.01, yv/s_); l2.set_data(w/abs(w).max(), yv/s_); a1.set_ylim(-6, 6); a1.set_ylabel("y/√(νt)")
        l3.set_data(rr/s_, ut*2*np.pi*rr/1e-4); a2.set_xlim(0, 10); a2.set_ylim(0, 1.05)
        a2.set_xlabel("r/√(νt)"); a2.set_ylabel("2πr u_θ/Γ")
        a1.set_title(f"t = {t_:.2f} s: rescaled — one curve", fontsize=9); a2.set_title("rescaled — one curve", fontsize=9)
    return []


show_animation(animate(update, frames=nfr, fig=fig, interval=90))
""")
see_read_change(r"""Left: the vortex sheet's velocity (blue) and vorticity (purple) spreading; right: the line vortex losing its peak. In
the last quarter of the movie both switch to rescaled axes and stop changing — the `A5` animation.""",
                r"""In the rescaled frames the sheet is plotted against y/√(νt) and the vortex as 2πr u_θ/Γ against r/√(νt): if the
curves still moved, the similarity form would be wrong.""",
                r"""…ν were halved: the raw movie runs √2 times slower in space; the rescaled frames are identical.""")
explainer("similarity_exponents", "Where do similarity exponents come from?", r"""
Set trial exponents n and m with two sliders: the rescaled profiles collapse onto one curve only at the right values, and the
powers of t in the reduced equation turn green when they match — the exponent algebra of Examples 8.4–8.7 as a game you can
win or lose.""",
          ["In 'vortex sheet' mode, start at n = 0.4: the curves fan out; slide to 0.5 and watch them merge.",
           "In 'spreading bead' mode try the diffusion guess m = 1/2 — the bracket chip turns red.",
           "Switch to 'line vortex': why is the prefactor now a power of r, not of t?",
           "Check the conserved-integral readout: it is flat only at the right n."])
note("N64 [C] · Two kinds of length.", r"""
Diffusive lengths grow like $(\nu t)^{1/2}$; the bead's $t^{1/5}$ is not a diffusion length but an advective one — how far the
fluid itself has travelled.
""")
whatif(r"""
…the wall did not start once but oscillated forever? Then the problem has an imposed time 1/ω; no similarity variable exists,
but the diffusion length $\sqrt{\nu/\omega}$ sets the depth (C11).
""")

# =====================================================================================================================
# A.5 §8.5 — C11 (N65–N74, D24, D25, A2, IF6, E7)
# =====================================================================================================================
nb.section("8.5", "Flow Due to an Oscillating Plate", intro=r"""
**What is this section about?** A plate shaking back and forth in its own plane drags a thin layer of fluid with it. The
motion decays and lags with height, looking like a wave travelling away from the wall — but it is diffusion. The depth
$\sqrt{2\nu/\omega}$ and the phase structure $(1+i)y/\delta$ are exactly those of tidal bottom boundary layers and of the
Ekman layer of Ch. 13.
""")
core("C11", "Stokes' second problem: the Stokes layer", r"""
How deep does an oscillation reach into a viscous fluid — and why do the deeper layers lag behind?
""", eqs=("8.38",))
problem(r"""
Shake a plate under water at one cycle per second; the tide rubs the sea floor back and forth twice a day; a sound wave slides
air along a wall a thousand times a second. How far above the wall does the fluid feel the shaking, and how late does each
layer follow?
""")
idea(r"""
y = 3δ   ~~ 5 %,  3 rad late
y = 2δ   ~~~ 14 %, 2 rad late
y = δ    ~~~~ 37 %, 1 rad late
y = 0    ◄══►  U cos ωt           δ = δ_e = √(2ν/ω)
""", table=r"""
| shaking | ν [m²/s] | δ_e (computed below) |
|---|---|---|
| plate at 1 Hz in water | 10⁻⁶ | 0.56 mm |
| M₂ tide (12.42 h), molecular ν | 10⁻⁶ | 0.12 m |
| M₂ tide, eddy ν (turbulent) | 10⁻² | 12 m |
| sound at 1 kHz in air | 1.5 × 10⁻⁵ | 0.069 mm |
""", words=r"Each layer repeats the wall's motion, delayed by y/δ radians and shrunk by $e^{-y/\delta}$.")
note("N65 [B], N66 [B], N67 [B] · Set-up.", r"""
The plate at y = 0 oscillates in its own plane; we want only the periodic state after start-up transients have died, so there
is no initial condition. The equation is (8.20) again, $\frac{\partial u}{\partial t}=\nu\frac{\partial^2u}{\partial y^2}$, but
now the wall imposes a time scale 1/ω, through the wall condition (8.33) below and $u(y\to\infty,t)=\text{bounded}$ *(8.34)*.
""", equation=EQ["8.33"], ref="8.33")
confusion(r"""
ω is now a **frequency** [rad/s] (as in Ch. 7), not a vorticity.
""")
remind([
    ("complex exponential solution, taking the real part", "write the wall motion as $\\mathrm{Re}\\{Ue^{i\\omega t}\\}$; for a linear equation with real coefficients we can solve with the complex form and take the real part at the end (Ch. 7 P176, complex amplitudes)."),
    ("square root of i", "$\\sqrt i=(1+i)/\\sqrt2$, since $((1+i)/\\sqrt2)^2=(1+2i-1)/2=i$ (Ch. 6 P159, complex square roots)."),
    ("Euler's formula, i² = −1", "$e^{i\\theta}=\\cos\\theta+i\\sin\\theta$, so $\\mathrm{Re}\\{e^{i\\theta}\\}=\\cos\\theta$ (Ch. 1 P45)."),
])
nb.code(r"""
print(np.sqrt(1j), (1 + 1j)/np.sqrt(2))                  # both (0.7071+0.7071j): the square root of i
""")
note("N68 [B], N69 [B], N70 [B], N71 [B] · The complex route.", r"""
**N68** The complex trial is (8.35) below. ⚠️ The text then says "substitution of (8.33) into (8.20)"; it is (8.35) that is
substituted. **N69** It turns (8.20) into $i\omega f=\nu\frac{d^2f}{dy^2}$ *(8.36)*. **N70** $f=e^{ky}$ with
$k=(i\omega/\nu)^{1/2}=\pm(i+1)(\omega/2\nu)^{1/2}$. **N71** The general solution is (8.37); bounded ⇒ B = 0, wall ⇒ A = U.
""", equation=EQ["8.35"] + r"\ \ \text{(8.35)},\qquad " + EQ["8.37"] + r"\ \ \text{(8.37)}")
D("D24", ref="8.38")
nb.code(r"""
s2 = ch08.stokes_second_sympy()                          # D24 in sympy
print("roots k:", s2["k_roots"], "  bounded:", s2["bounded_root"])
print("residual of (8.20) for (8.38):", s2["residual"])  # 0
""")
note("N72 [B] · Reading (8.38).", r"""
The cosine moves toward +y — it looks like a damped transverse wave — but there is no restoring force: it is diffusion driven
by a periodic wall. At $y=4(\nu/\omega)^{1/2}$ the amplitude is $Ue^{-4/\sqrt2}$ (0.05911U, printed by the code below), so the
book calls $\delta\sim4(\nu/\omega)^{1/2}$ the layer thickness; the literature's e-folding depth is $\delta_e=\sqrt{2\nu/\omega}$,
2√2 = 2.83 times smaller (ours: crests move at $\sqrt{2\nu\omega}$, "wavelength" $2\pi\delta_e$).
""")
remind([
    ("phase of a cosine and a time lag", "in cos(ωt − φ) the phase lag φ [rad] means the peak comes φ/ω seconds later (Ch. 7 P165, phase of a wave)."),
])
D("D25")
nb.worked_example("a plate shaken once a second in water", r"""
ω = 2π rad/s, ν = 10⁻⁶ m²/s, U = 0.1 m/s.
1. $\delta_e=\sqrt{2\nu/\omega}=\sqrt{2\times10^{-6}/6.283}=5.64\times10^{-4}$ m = 0.56 mm.
2. Book depth $4\sqrt{\nu/\omega}=1.60$ mm.
3. At y = 1 mm: amplitude $e^{-1/0.564}=e^{-1.77}=0.17$ → 1.7 cm/s; phase lag 1.77 rad, i.e. 1.77/(2π) × 1 s = 0.28 s late.
4. Crest speed $\omega\delta_e=3.5$ mm/s.
5. Shake 100× faster: every length shrinks 10×.
""")
nb.code(r"""
y = np.array([0.0, 0.5e-3, 1e-3, 2e-3])                  # four heights [m]
print("u(t = 0) =", ch08.stokes_second_problem(y, 0.0, U=0.1, omega=2*np.pi, nu=NU_W), "m/s")   # (8.38)
sl = ch08.stokes_layer(NU_W, 2*np.pi)                    # depths and speeds of the 1 Hz layer
print({k: f"{sl[k]:.4g}" for k in ("delta_e", "delta_book", "phase_speed", "wavelength", "amplitude_at_delta_book")})
ss = ch08.stokes_layer_state(NU_W, 2*np.pi, 1e-3)        # what the layer does at y = 1 mm
print({k: f"{ss[k]:.4g}" for k in ("amplitude", "phase_lag", "time_lag")})
print(f"e^(−2√2) = {np.exp(-2*np.sqrt(2)):.5f}")          # N72's amplitude at the book's depth
wM2 = 2*np.pi/(12.42*3600)                               # the M₂ tide's frequency [rad/s]
print(f"M2 tide: δ_e = {ch08.stokes_layer(1e-6, wM2)['delta_e']:.3f} m (molecular ν), {ch08.stokes_layer(1e-2, wM2)['delta_e']:.1f} m (eddy ν)")
print(f"sound at 1 kHz in air: δ_e = {1e3*ch08.stokes_layer(1.5e-5, 2*np.pi*1e3)['delta_e']:.3f} mm")
for g_ in ch01.pi_groups({"u": "m/s", "U": "m/s", "y": "m", "t": "s", "nu": "m^2/s", "omega": "1/s"}):   # N73's groups
    print("  group:", " · ".join(f"{k}^{v}" for k, v in g_.items() if v != 0))
""", explain=r"""
1. `stokes_second_problem` is (8.38): at t = 0 the wall moves at U and the layer at 1 mm already moves backward (the lag).
2. `stokes_layer` gives δ_e = 0.564 mm, the book's depth 1.596 mm, crest speed 3.545 mm/s, "wavelength" 2πδ_e = 3.545 mm
   and the amplitude 0.05911 at the book's depth.
3. `stokes_layer_state` at 1 mm: amplitude 0.170, phase lag 1.772 rad, time lag 0.282 s — the tiny example.
4. The tide numbers fill the table above: 12 cm with molecular viscosity, 12 m with an eddy viscosity.
5. **N73 [B]** `pi_groups` finds four groups for (u, U, y, t, ν, ω); linearity removes the one that contains U, leaving three,
   u/U, ωt and $y(\omega/\nu)^{1/2}$: y and t never combine, so there is **no similarity variable**; the extent $(\nu/\omega)^{1/2}$
   is viscosity times the imposed period.
""")
nb.check_agree(r"""
de = np.sqrt(2*NU_W/(2*np.pi)); t_ = 0.3                  # e-folding depth [m] and a time [s]
u_mine = np.real(0.1*np.exp(1j*2*np.pi*t_)*np.exp(-(1 + 1j)*y/de))   # Re{U e^{iωt} e^{−(1+i)y/δ_e}} (D24 steps 9–10)
assert np.allclose(u_mine, ch08.stokes_second_problem(y, t_, U=0.1, omega=2*np.pi, nu=NU_W), rtol=1e-12, atol=1e-15)
T = 1.0; nper = 10 if not FAST else 5                    # period [s] and number of periods to march from rest
yy = np.linspace(0, 20*de, 201)                          # 20 e-folding depths: the far end sees nothing [m]
u_num = ch08.crank_nicolson_1d(np.zeros(201), yy, T/200, 200*nper, NU_W, lambda tt: 0.1*np.cos(2*np.pi*tt), 0.0)   # from rest
err = np.max(np.abs(u_num - ch08.stokes_second_problem(yy, nper*T, U=0.1, omega=2*np.pi, nu=NU_W))[yy <= 6*de])
print(f"Crank–Nicolson from rest after {nper} periods vs (8.38): max error {err/0.1:.1e} U within 6δ_e")
assert err < 2e-3*0.1                                    # the start-up transient has died: the periodic state is (8.38)
""")
nb.figure(r"""
yd = np.linspace(0, 6, 300)                               # y/δ_e = y√(ω/2ν), the book's axis [–]
fig, ax = plt.subplots(figsize=(6.5, 4.2))
for k, wt in enumerate((0, np.pi/2, np.pi, 3*np.pi/2)):  # four phases of the wall
    u = ch08.stokes_second_problem(yd*de, wt/(2*np.pi), U=1.0, omega=2*np.pi, nu=NU_W)   # (8.38) with U = 1
    ax.plot(u, yd, color=plt.cm.Blues(0.45 + 0.18*k), lw=2, label=f"ωt = {['0', 'π/2', 'π', '3π/2'][k]}")
ax.plot(np.exp(-yd), yd, color=COLORS["muted"], ls="--", lw=1, label="envelope ±e^(−y/δ_e)")
ax.plot(-np.exp(-yd), yd, color=COLORS["muted"], ls="--", lw=1)
yb = 2*np.sqrt(2)                                          # the book's depth 4√(ν/ω) in units of δ_e
ax.axhline(yb, color=COLORS["amber"], lw=1.2)
ax.text(0.35, yb + 0.1, f"book depth 4√(ν/ω) = 2√2 δ_e: amplitude {np.exp(-yb):.4f}", color=COLORS["amber"], fontsize=8)
ax.set_xlabel("u/U [–]"); ax.set_ylabel("y/δ_e = y√(ω/2ν) [–]"); ax.legend(fontsize=8, loc="upper right")
ax.set_title("Each layer is the wall's motion, delayed and shrunk", fontsize=10)
savefig(fig, "ch08", "c11_stokes_layer"); plt.show()
""", see=r"""Four S-shaped profiles (one per quarter period) inside a narrowing funnel, the envelope $\pm e^{-y/\delta_e}$; the amber
line is the book's depth 4√(ν/ω), where 0.0591 of the motion is left — our remake of Fig. 8.16 (part of `N107`).""",
    read=r"""At any height the spread of the curves is the local amplitude; the height where a curve crosses zero moves up with ωt —
that is the apparent "wave".""",
    change=r"""…the fluid were ten times more viscous: all lengths grow √10 ≈ 3.2 times; in these rescaled axes nothing moves.""")
nb.animation(r"""
nfr = 60 if not FAST else 30                              # frames (two periods)
tf = np.linspace(0, 2.0, nfr)                             # time [s], 1 Hz plate
de = np.sqrt(2*NU_W/(2*np.pi))                            # e-folding depth [m]
yv = np.linspace(0, 5*de, 150)                            # heights [m]
ytr = de*np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0])         # six tracer heights [m]
xs0 = np.linspace(0, 1, 6)                                # their resting x positions [schematic units]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
trac = a1.scatter(xs0, 1e3*ytr, c=[COLORS["teal"]]*6, s=40, zorder=3)   # tracers
dye, = a1.plot(np.zeros(150), 1e3*yv, color=COLORS["blue"], lw=2)      # a dyed vertical line
plate, = a1.plot([-0.2, 1.2], [0, 0], color=COLORS["ink"], lw=5)       # the plate
a1.set_xlim(-0.5, 1.5); a1.set_ylim(-0.1, 1e3*5*de); a1.set_xlabel("x (displacement × 50, schematic)"); a1.set_ylabel("y [mm]")
prof, = a2.plot([], [], color=COLORS["blue"], lw=2)
a2.plot(np.exp(-yv/de), 1e3*yv, color=COLORS["muted"], ls="--", lw=1); a2.plot(-np.exp(-yv/de), 1e3*yv, color=COLORS["muted"], ls="--", lw=1)
probe, = a2.plot([], [], "o", color=COLORS["orange"], ms=8)            # probe at y = δ_e
a2.set_xlim(-1.1, 1.1); a2.set_xlabel("u/U [–]")
Dp = lambda y_, t_: np.imag(np.exp(1j*2*np.pi*t_)*np.exp(-(1 + 1j)*y_/de))/(2*np.pi)   # displacement ∫u dt for U = 1 [m·s/m]


def update(i):                                            # frame i
    t_ = tf[i]
    u = ch08.stokes_second_problem(yv, t_, U=1.0, omega=2*np.pi, nu=NU_W)   # (8.38) with U = 1
    prof.set_data(u, 1e3*yv)
    probe.set_data([ch08.stokes_second_problem(de, t_, U=1.0, omega=2*np.pi, nu=NU_W)], [1e3*de])
    dye.set_xdata(0.5 + 3*Dp(yv, t_))                     # the dye line bends with the lagging layers (displacement × 3)
    trac.set_offsets(np.c_[xs0 + 3*Dp(ytr, t_), 1e3*ytr]) # tracers oscillate with growing lag
    plate.set_xdata(np.array([-0.2, 1.2]) + 3*Dp(0.0, t_))   # the plate itself
    a1.set_title(f"t = {t_:.2f} s (1 Hz plate)", fontsize=10); a2.set_title("profile and probe at y = δ_e", fontsize=10)
    return []


show_animation(animate(update, frames=nfr, fig=fig, interval=60))
""")
see_read_change(r"""Left: tracers at six heights and a dyed line oscillating with the plate, each layer later and smaller than the one
below; right: the profile swinging inside its envelope, with an orange probe at y = δ_e — the `A2` animation.""",
                r"""The probe moves e⁻¹ = 37 % as far as the plate and peaks 1 rad (0.16 s) later; the dye line is an S that travels
upward — the apparent wave of (8.38).""",
                r"""…the probe moved to 2δ_e: its dot would move e⁻² = 14 % as far and peak 2 rad later.""")
nb.plotly(r"""
de = np.sqrt(2*NU_W/(2*np.pi))                            # e-folding depth [m]
yd = np.linspace(0, 6, 120)                               # y/δ_e [–]


def frame(wt):                                            # curves at phase ωt
    return {"u/U, (8.38)": (ch08.stokes_second_problem(yd*de, wt/(2*np.pi), U=1.0, omega=2*np.pi, nu=NU_W), yd),
            "envelope +": (np.exp(-yd), yd), "envelope −": (-np.exp(-yd), yd)}


fig = animate_figure(frame, np.linspace(0, 2*np.pi, 25 if not FAST else 12), time_label="ωt", unit="rad",
                     xlabel="u/U [–]", ylabel="y/δ_e [–]", title="The Stokes layer through one period", xrange=[-1.05, 1.05], yrange=[0, 6])
recolor(fig, {"u/U, (8.38)": COLORS["blue"], "envelope +": COLORS["muted"], "envelope −": COLORS["muted"]},
        {"envelope +": "dash", "envelope −": "dash"})
fig.show()
""", explain=r"""
The page-surviving version of the animation: Play runs one period of (8.38) inside its envelope.
""")
explainer("oscillating_plate", "How deep does a shaking plate reach?", r"""
The animation shows the phase lag and the envelope together; drag ω and the layer shrinks while the 0.059U marker at
4√(ν/ω) follows, and a probe clock shows a deep layer peaking later — none of which a four-phase static figure makes vivid.""",
          ["Play the 1 Hz preset and click at y = 1 mm: read the amplitude and the time lag.",
           "Switch to the M₂ tide with molecular and then eddy viscosity: 12 cm versus 12 m.",
           "Double ω and predict δ_e before you look (it falls by √2).",
           "Open the 'Right now' notes: which real Stokes layer is closest to your setting?"])
note("N74 [C] · Where the Stokes layer reappears.", r"""
The same layer explains why sound is weakly absorbed at flat walls: the viscous acoustic boundary layer (Ch. 15). For the
climate thread: replace the plate's oscillation by the Coriolis force and the (1 + i)/δ structure of (8.37) becomes the Ekman
spiral (Ch. 13).
""")
whatif(r"""
…the plate stood still and a sphere moved slowly through the fluid instead? Then there is no single wall-normal direction and
the flow is 3-D; for tiny Reynolds numbers inertia vanishes altogether and the equations become linear again (C12).
""")

# =====================================================================================================================
# A.6 §8.6 — R15, R16, C12 (N75–N80, D26, P198, IF7), R17–R19, C13 (N81–N87, N93, N107, D27–D29, P199, A4, E8),
#            R20, C14 (N88–N92, D30, D31, E9), C15 (N94–N101, D32, D33, IF8)
# =====================================================================================================================
nb.section("8.6", "Low Reynolds Number Viscous Flow Past a Sphere", intro=r"""
**What is this section about?** Tiny, slow things — cloud droplets, silt, bacteria, a bead sinking in syrup — live at Reynolds
numbers far below one. We find the right way to drop inertia (pressure must be rescaled, not dropped), solve the resulting
linear Stokes equations for a sphere, integrate the surface stresses to Stokes' drag $D=6\pi\mu aU$, and use it for the fall
speed of droplets and for Millikan's measurement of the electron's charge. Finally we see why the solution fails far from the
sphere and how Oseen repaired it.
""")
nb.recap("R15", "Steady Navier–Stokes", r"""
For steady flow round a body of size L at speed U: $\rho\mathbf u\cdot\nabla\mathbf u+\nabla p=\mu\nabla^2\mathbf u$ *(8.39)* —
Ch. 4's $\rho\frac{D\mathbf u}{Dt}=-\nabla p+\rho\mathbf g+\mu\nabla^2\mathbf u$ *(4.39b)* without the time derivative and with
gravity absorbed into p.
""", where="Ch. 4 §4.6")
nb.recap("R16", "The high-Re scaling", r"""
With $u^*=u/U$, $x^*=x/L$ and $p^*=(p-p_\infty)/\rho U^2$ *(4.100)*, (8.39) becomes
$\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*=\frac{1}{\mathrm{Re}}\nabla^{*2}\mathbf u^*$ *(8.40)*, Re = ρUL/μ; for Re ≫ 1 the
viscous term looks negligible (Euler, Ch. 6).
""", where="Ch. 4 §4.11")
nb.code(r"""
print(SIM.nondimensional_ns_coefficients("dynamic", "advective"))   # pressure scaled by ρU²: viscous coefficient μ/(ρUl) = 1/Re
print(SIM.nondimensional_ns_coefficients("viscous", "advective"))   # pressure scaled by μU/l: pressure coefficient μ/(ρUl) too
""")
core("C12", "Creeping flow: the Stokes equations", r"""
Which terms survive when Re → 0 — and why does pressure not disappear together with inertia?
""", eqs=("8.43",))
problem(r"""
A bacterium swims at a few body lengths per second; a grain of silt sinks through a river; a mist droplet drifts down. For
them Re = UL/ν ≪ 1. Inertia is negligible — but if we simply multiply (8.40) by Re and let it go to zero, the pressure
disappears too, and the resulting equation cannot hold a sphere in a stream. What went wrong, and what is the right limit?
""")
idea("", table=r"""
| regime | pressure balances | pressure scale | dimensionless form |
|---|---|---|---|
| Re ≫ 1 | inertia | ρU² | $\mathbf u^*\cdot\nabla^*\mathbf u^*+\nabla^*p^*=\frac{1}{\mathrm{Re}}\nabla^{*2}\mathbf u^*$ *(8.40)* |
| Re ≪ 1 | viscous stress | μU/L | $\mathrm{Re}\,(\mathbf u^*\cdot\nabla^*\mathbf u^*)=-\nabla^*p^*+\nabla^{*2}\mathbf u^*$ *(8.42)* |
""", words=r"The size of pressure differences is set by whatever they balance.")
note("N75 [C], N76 [C] · Expansions in a parameter.", r"""
**N75** Many problems are solved as expansions in a small or large parameter (here Re); Ch. 9 (1/Re) and Ch. 11 (small
disturbances) do the same. **N76** At high Re the expansion in 1/Re fails near walls, where the inviscid solution cannot satisfy
no slip — the boundary layer of Ch. 9; the two-length scaling (8.14) of C05 was the hint.
""")
P("P198", "dominant balance", r"""
Choose the scale of a quantity from the term it must balance. If pressure differences are pushed by viscous stresses, then
$p-p_\infty\sim\mu U/L$, not ρU²; picking the wrong scale makes a term look negligible when it is not.
""", code=r"""
mu, U, L, rho = 1e-3, 1e-3, 1e-5, 1000.0                 # a 10 µm particle at 1 mm/s in water
print(rho*U**2, mu*U/L)                                  # dynamic 1e-3 Pa vs viscous 0.1 Pa: pressure is viscous here
""")
note("N77 [B], N78 [B], N79 [B] · The two scalings.", r"""
**N77** Multiplying (8.40) by Re gives (8.41): as Re → 0 it leaves 0 = μ∇²u — pressure lost too. **N78** The low-Re pressure
scale is $p^*=(p-p_\infty)L/(\mu U)$. **N79** It gives (8.42), whose Re → 0 limit is (8.43).
""", equation=EQ["8.41"] + r"\ \ \text{(8.41)},\qquad " + EQ["8.42"] + r"\ \ \text{(8.42)}")
D("D26", ref="8.43")
nb.worked_example("a 10 µm cloud droplet", r"""
a = 10 µm, falling at U ≈ 1.2 cm/s in air (ν = 1.5 × 10⁻⁵ m²/s, μ = 1.81 × 10⁻⁵ Pa s; the speed is derived in C14).
1. $\mathrm{Re}=2aU/\nu=2\times10^{-5}\times0.012/1.5\times10^{-5}=0.016$.
2. Dynamic pressure $\rho U^2=1.2\times1.44\times10^{-4}=1.7\times10^{-4}$ Pa.
3. Viscous pressure $\mu U/a=1.81\times10^{-5}\times0.012/10^{-5}=0.022$ Pa — 125 times larger (= 2/Re). The viscous scale is
   the right one.
""")
nb.code(r"""
lr = ch08.low_re_scaling_sympy()                         # D26's coefficients in sympy (in terms of Re)
for k in ("dynamic", "dynamic_times_Re", "viscous_times_Re"):
    print(f"{k:17s}", {t_: lr[k][t_] for t_ in ("inertia", "pressure", "viscous")})
u_fn = lambda x, t=0.0: np.stack(ch08.stokes_sphere_velocity_xyz(x[0], x[1], x[2], U=1e-3, a=1e-5))   # C13's field (8.49)
p_fn = lambda x, t=0.0: ch08.stokes_sphere_pressure(np.sqrt((x**2).sum(0)), np.arccos(x[0]/np.sqrt((x**2).sum(0))),
                                                    U=1e-3, a=1e-5, mu=MU_W)   # C14's pressure (8.50)
xp = np.array([3e-5, 1e-5, 0.5e-5])                      # a point 3.3 radii from the centre [m]
res = ch08.stokes_residual(u_fn, p_fn, xp, mu=MU_W, h=1e-8)   # ∇p − μ∇²u by central differences (step a/1000)
scale = 3*MU_W*1e-3*1e-5/np.sqrt((xp**2).sum())**3       # the size of ∇p there, 3μUa/r³ [Pa/m]
print(f"Stokes residual |∇p − μ∇²u| = {np.max(np.abs(res)):.1e} Pa/m against ∇p ~ {scale:.0f} Pa/m (relative {np.max(np.abs(res))/scale:.0e})")
""", explain=r"""
1. `low_re_scaling_sympy` gives the coefficients of inertia, pressure and viscous terms: with the dynamic pressure scale
   {1, 1, 1/Re}; multiplied by Re {Re, Re, 1} — the pressure dies with inertia; with the viscous pressure scale (times Re)
   {Re, 1, 1} — only inertia dies.
2. `stokes_residual` evaluates ∇p − μ∇²u for the sphere solution of C13–C14 by finite differences: about 10⁻⁷ of the pressure
   gradient, i.e. the field satisfies (8.43) to the accuracy of the differences.
""")
nb.check_agree(r"""
Re_ = 0.016                                              # the droplet's Reynolds number
dyn, visc = [1, 1, 1/Re_], [Re_, 1, 1]                   # by hand: dynamic scaling, and viscous scaling × Re
Re_sym = sp.Symbol("Re", positive=True)                  # the symbol the engine uses
assert np.allclose(dyn, [float(lr["dynamic"][k].subs(Re_sym, Re_)) if hasattr(lr["dynamic"][k], "subs") else float(lr["dynamic"][k])
                         for k in ("inertia", "pressure", "viscous")])
assert np.allclose(visc, [float(lr["viscous_times_Re"][k].subs(Re_sym, Re_)) if hasattr(lr["viscous_times_Re"][k], "subs") else float(lr["viscous_times_Re"][k])
                          for k in ("inertia", "pressure", "viscous")])
print("hand coefficients = engine coefficients at Re = 0.016")
""")
nb.figure(r"""
Re = np.geomspace(1e-3, 1e3, 100)                          # Reynolds numbers
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
a1.loglog(Re, Re, color=COLORS["teal"], lw=2, label="inertia (Re)")
a1.loglog(Re, Re, color=COLORS["orange"], lw=2, ls="--", label="pressure (Re)")
a1.loglog(Re, 1 + 0*Re, color=COLORS["rose"], lw=2, label="viscous (1)")
a1.set_title("(8.41): dynamic pressure scale × Re — pressure dies with inertia", fontsize=9)
a2.loglog(Re, Re, color=COLORS["teal"], lw=2, label="inertia (Re)")
a2.loglog(Re, 1 + 0*Re, color=COLORS["orange"], lw=2, label="pressure (1)")
a2.loglog(Re, 1.05 + 0*Re, color=COLORS["rose"], lw=2, ls="--", label="viscous (1)")
a2.set_title("(8.42): viscous pressure scale — only inertia dies", fontsize=9)
for ax in (a1, a2):
    ax.axvline(1, color=COLORS["muted"], lw=0.8); ax.set_xlabel("Re [–]"); ax.legend(fontsize=7)
a1.set_ylabel("coefficient of the term [–]")
fig.suptitle("Rescale pressure, not just inertia", fontsize=11)
savefig(fig, "ch08", "c12_low_re_scaling"); plt.show()
""", see=r"""Left: with the pressure scaled by ρU² and the equation multiplied by Re, the pressure coefficient falls together with
inertia; right: with the pressure scaled by μU/L, pressure and viscous terms stay at 1 and only inertia falls.""",
    read=r"""Read off which terms survive as Re → 0 (the left end): on the right, pressure and friction balance — the Stokes equations
$\nabla p=\mu\nabla^2\mathbf u$ *(8.43)*.""",
    change=r"""…Re = 1: the two scalings coincide — neither term is negligible and the full equation is needed.""")
nb.plotly(r"""
def frame(lr10):                                          # coefficients at Re = 10^lr10
    R_ = 10**lr10
    return {"dynamic scale × Re: inertia, pressure, viscous": ([1, 2, 3], np.log10([R_, R_, 1.0])),
            "viscous scale × Re: inertia, pressure, viscous": ([1, 2, 3], np.log10([R_, 1.0, 1.0]))}


fig = slider_figure(frame, "log10 Re", np.linspace(-3, 3, 25 if not FAST else 13), unit="",
                    xlabel="term: 1 = inertia, 2 = pressure, 3 = viscous", ylabel="log10(coefficient)",
                    title="Which terms survive as Re → 0?", modes={"dynamic scale × Re: inertia, pressure, viscous": "lines+markers",
                                                                  "viscous scale × Re: inertia, pressure, viscous": "lines+markers"})
recolor(fig, {"dynamic scale × Re: inertia, pressure, viscous": COLORS["muted"], "viscous scale × Re: inertia, pressure, viscous": COLORS["orange"]})
fig.show()
""", explain=r"""
Slide log₁₀ Re down to −3: with the dynamic scale (grey) both inertia and pressure drop below the viscous term; with the
viscous scale (orange) the pressure stays level with the viscous term.
""")
note("N80 [C] · The italic rule.", r"""
The right length and time scales depend on the region of the flow and are found by balancing the terms that matter there.
Ch. 9 (boundary layers) and Ch. 13 (rotating, stratified scales) live by it.
""")
whatif(r"""
…we solve (8.43) round a sphere? The equation is linear, so reversing the stream reverses the whole flow: no wake (C13).
""")

nb.recap("R17", "Vorticity in spherical coordinates", r"""
For an axisymmetric flow without swirl only the azimuthal vorticity survives,
$\omega_\varphi=\frac1r\big[\frac{\partial(ru_\theta)}{\partial r}-\frac{\partial u_r}{\partial\theta}\big]$ (the spherical curl of
`core.curvilinear`).
""", where="Ch. 3 §3.4")
nb.recap("R18", "Stokes' stream function", r"""
Axisymmetric velocities come from one function: $u_r=\frac1{r^2\sin\theta}\frac{\partial\psi}{\partial\theta},\
u_\theta=-\frac1{r\sin\theta}\frac{\partial\psi}{\partial r}$ *(6.83)*.
""", where="Ch. 6 §6.8")
nb.recap("R19", "The uniform stream far away", r"""
A uniform stream U along the axis has $\psi=\tfrac12Ur^2\sin^2\theta$ *(6.86)*; the far-field condition is
$\psi(r\to\infty,\theta)=\tfrac12Ur^2\sin^2\theta$ *(8.47)*.
""", where="Ch. 6 §6.8")
core("C13", "Stokes' solution for the sphere", r"""
What is the flow round a sphere when viscosity dominates everywhere?
""", eqs=("8.48",))
problem(r"""
A tiny bead sinks through syrup. Drawn in the bead's frame, the fluid streams past it; drawn in the syrup's frame, the fluid is
pushed aside in front and closes in behind. We want the whole velocity field — it gives the pressure, the stresses and the
drag (C14), and it is the reference for every settling particle and swimming micro-organism.
""")
idea(r"""
─────────►        ψ = ½Ur² sin²θ  −  (3a/4)Ur sin²θ  +  (a³/4)(U/r) sin²θ
────╮ ● ╭─►          stream          'Stokeslet' ~ a/r      'doublet' ~ a³/r³
─────────►        the a/r part decays so slowly that the sphere is felt tens of radii away
""", words=r"A uniform stream minus a slowly decaying disturbance, fore–aft symmetric.")
remind([
    ("spherical coordinates with θ from the downstream axis", "r from the centre, θ the angle from the +x axis — here the **downstream** direction, so the rear stagnation point is θ = 0 and the front θ = π — and φ round the axis (Ch. 3 P88)."),
    ("curl of a gradient is zero, curl commutes with the Cartesian Laplacian", "∇ × ∇p = 0 because mixed partial derivatives are equal (Schwarz, Ch. 4 P121); in Cartesian components curl and ∇² are both built from constant-coefficient derivatives, so their order can be swapped."),
    ("index notation ε_ijk", "$(\\nabla\\times\\mathbf A)_k=\\varepsilon_{kli}\\partial_lA_i$ with the permutation symbol $\\varepsilon_{kli}$ = +1 for cyclic, −1 for anticyclic orders, 0 otherwise (Ch. 2 §2.7, P72)."),
])
note("N81 [B] · Take the curl.", r"""
The curl of (8.43) kills the pressure (curl of a gradient) and leaves ∇²ω = 0 — in Cartesian components; in spherical
coordinates the operator means −∇×∇×ω (the book's footnote).
""")
D("D27", check_src=r"""
S = ch08.stokes_sphere_sympy()                           # the Stokes-sphere results in sympy (cached; shared by D28–D31)
r, th, a, U, mu = S["symbols"]                           # its symbols: r, θ, a, U, μ
X = (r, th, sp.Symbol("phi", real=True))                 # spherical coordinates for core.curvilinear
w = [0, 0, S["omega_phi"]]                               # Stokes' vorticity ω = ω_φ e_φ, ω_φ = −(3Ua/2r²) sin θ
cc = CL.curl(CL.curl(w, "spherical", X), "spherical", X) # ∇ × (∇ × ω), component by component
print("ω_φ =", S["omega_phi"])
print("∇×∇×ω =", [sp.simplify(c) for c in cc])          # [0, 0, 0]: −∇×∇×ω = 0 holds for Stokes' solution
assert all(sp.simplify(c) == 0 for c in cc)
""")
P("P199", "the Stokes operator E² applied twice", r"""
For axisymmetric flow the operator $E^2=\frac{\partial^2}{\partial r^2}+\frac{\sin\theta}{r^2}\frac{\partial}{\partial\theta}\big(
\frac1{\sin\theta}\frac{\partial}{\partial\theta}\big)$ (Ch. 6 used it once, (6.77)) turns ψ into the vorticity:
$\omega_\varphi=-E^2\psi/(r\sin\theta)$. Applying it twice, $E^2(E^2\psi)$, is **not** the biharmonic ∇⁴ψ — a sympy test in D28
shows the difference.
""", code=r"""
import sympy as sp
r, th = sp.symbols('r theta', positive=True)
E2 = lambda f: sp.diff(f, r, 2) + sp.sin(th)/r**2*sp.diff(sp.diff(f, th)/sp.sin(th), th)   # the Stokes operator
print(sp.simplify(E2(r**2*sp.sin(th)**2)))               # 0: the uniform stream has no vorticity
""")
note("N82 [B], N83 [B] · Vorticity from ψ, and the ψ-equation.", r"""
**N82** $\omega_\varphi=-\frac1r\big[\frac{1}{\sin\theta}\frac{\partial^2\psi}{\partial r^2}+\frac1{r^2}\frac{\partial}{\partial\theta}
\big(\frac1{\sin\theta}\frac{\partial\psi}{\partial\theta}\big)\big]=-\frac{E^2\psi}{r\sin\theta}$. **N83** Combining with ∇²ω_φ = 0
in the footnote's sense gives (8.44). ⚠️ It is the square of E², not the biharmonic (footnote 2).
""", equation=EQ["8.44"], ref="8.44")
D("D28", ref="8.44", check_src=r"""
A_ = sp.Function("A")(r, th)                             # a generic azimuthal field A(r, θ) e_φ
cc = CL.curl(CL.curl([0, 0, A_], "spherical", X, False), "spherical", X, False)   # ∇×∇×(A e_φ), not simplified yet
E2 = lambda f: sp.diff(f, r, 2) + sp.sin(th)/r**2*sp.diff(sp.diff(f, th)/sp.sin(th), th)   # the Stokes operator E²
assert sp.simplify(cc[0]) == 0 and sp.simplify(cc[1]) == 0   # steps 7–9: the curl–curl of an azimuthal field is azimuthal …
assert sp.simplify(cc[2] + E2(r*sp.sin(th)*A_)/(r*sp.sin(th))) == 0   # … and equals −E²(r sinθ A)/(r sinθ) e_φ
print("E²(E²ψ) for (8.48):", sp.simplify(E2(E2(S["psi"]))))   # 0: Stokes' ψ satisfies (8.44)
print("∇⁴ψ for (8.48):   ", S["biharmonic_psi"])          # not 0: the scalar biharmonic is a different operator
assert sp.simplify(S["E4psi"]) == 0 and sp.simplify(S["biharmonic_psi"]) != 0
""")
note("N84 [B], N85 [B], N86 [B] · Conditions and the separated ODE.", r"""
**N84, N85** On the sphere: $\psi(r=a,\theta)=0$ *(8.45)* (no flow through it) and $\partial\psi(r=a,\theta)/\partial r=0$ *(8.46)*
(no slip); far away (8.47). **N86** (8.47) suggests ψ = f(r) sin²θ; (8.44) becomes the ODE below, solved by
$f=Ar^4+Br^2+Cr+D/r$ with A = 0, B = U/2, C = −3Ua/4, D = Ua³/4.
""", equation=r"f^{iv}-\frac{4f''}{r^2}+\frac{8f'}{r^3}-\frac{8f}{r^4}=0")
D("D29", ref="8.48")
note("N87 [B] · Velocities", r"""
(one derivative each through (6.83)):
""", equation=EQ["8.49"], ref="8.49")
nb.code(r"""
print("u_r =", sp.factor(S["u_r"]))                      # (8.49), factored: vanishes at r = a
print("u_θ =", sp.factor(S["u_theta"]))
print("on r = a:", S["u_r"].subs(r, a), S["u_theta"].subs(r, a))   # 0 0: no slip
print("f-roots:", S["roots"], " constants:", S["constants"])        # D29 steps 6 and 8–10
""")
nb.worked_example("a bead in syrup, two radii out on the side", r"""
a = 1 mm, U = 1 mm/s, ν = 10⁻³ m²/s (syrup: Re = 2aU/ν = 0.002). At r = 2a, θ = π/2: u_r = 0,
$u_\theta=-U\big(1-\frac38-\frac1{32}\big)=-0.594U$. Ideal flow (Ch. 6) would give $-U\big(1+\frac1{16}\big)=-1.063U$: the viscous
fluid is *slowed* beside the sphere instead of sped up. At r = 10a the Stokes speed is 0.925U — still a 7.5 % deficit; the
ideal-flow disturbance is already 0.05 %.
""")
nb.code(r"""
for rr_ in (2.0, 10.0):                                  # two distances on the side line, in radii
    print(f"r = {rr_:g}a: Stokes |u| = {ch08.side_line_speed(rr_*1e-3, U=1e-3, a=1e-3, model='stokes')/1e-3:.4f} U, "
          f"ideal |u| = {ch08.side_line_speed(rr_*1e-3, U=1e-3, a=1e-3, model='ideal')/1e-3:.4f} U")
print("(u_r, u_θ) at r = 2a, θ = π/2:", ch08.stokes_sphere_velocity(2e-3, np.pi/2, U=1e-3, a=1e-3), "m/s")
print("fluid-frame ψ at r = 2a, θ = π/4:", ch08.stokes_sphere_streamfunction(2e-3, np.pi/4, U=1e-3, a=1e-3, frame="fluid"), "m³/s")
""", explain=r"""
1. `side_line_speed` evaluates the speed at θ = π/2 for Stokes' (8.49) and for Ch. 6's ideal sphere: 0.594U against 1.0625U at
   2a, 0.925U against 1.0005U at 10a — the tiny example.
2. `stokes_sphere_velocity` gives the two components; `stokes_sphere_streamfunction(..., frame="fluid")` subtracts the stream
   (N93 below).
""")
nb.check_agree(r"""
rng = np.random.default_rng(3)                           # reproducible random points
rr = 1e-3*(1.2 + 8*rng.random(20)); tt = 0.1 + 2.9*rng.random(20)   # 20 points outside the sphere: r in (1.2a, 9.2a), θ in (0.1, 3.0)
hh = 1e-7                                                # finite-difference step [m or rad]
psi = lambda R_, T_: ch08.stokes_sphere_streamfunction(R_, T_, U=1e-3, a=1e-3)   # (8.48)
ur_fd = (psi(rr, tt + hh) - psi(rr, tt - hh))/(2*hh)/(rr**2*np.sin(tt))     # u_r = (1/r² sinθ) ∂ψ/∂θ (6.83)
ut_fd = -(psi(rr + hh, tt) - psi(rr - hh, tt))/(2*hh)/(rr*np.sin(tt))       # u_θ = −(1/r sinθ) ∂ψ/∂r
assert np.allclose((ur_fd, ut_fd), ch08.stokes_sphere_velocity(rr, tt, U=1e-3, a=1e-3), rtol=1e-6, atol=1e-12)
print("central differences of ψ through (6.83) reproduce (8.49) at 20 random points")
""")
note("N93 [B] · The sphere moving through still fluid (fluid frame).", r"""
Subtract the stream: ψ below (Fig. 8.19). The pattern is fore–aft symmetric — no wake — because (8.43) is linear: reversing U
maps u → −u and p − p∞ → −(p − p∞). ⚠️ The text calls the equation "(9.63)"; it means (8.43). Ch. 16 returns to this
reversibility (why a swimming bacterium cannot use a reciprocal stroke).
""", equation=r"\psi=Ur^2\sin^2\theta\Big(-\frac{3a}{4r}+\frac{a^3}{4r^3}\Big)")
remind([
    ("meshgrid, contour and streamplot", "`X, Y = np.meshgrid(x, y)` gives every point; `ax.contour(X, Y, psi, levels)` draws streamlines as lines of constant ψ (Ch. 2 P78)."),
])
nb.figure(r"""
n_ = 301 if not FAST else 151                              # grid points per side
xg = np.linspace(-6, 6, n_); yg = np.linspace(-4, 4, n_)   # in sphere radii (a = 1)
X, Y = np.meshgrid(xg, yg)
Rr = np.hypot(X, Y); Th = np.arccos(np.clip(X/np.maximum(Rr, 1e-12), -1, 1))   # spherical r and θ from +x
out = Rr > 1                                              # outside the sphere
psi_b = np.where(out, ch08.stokes_sphere_streamfunction(Rr, Th, 1.0, 1.0, frame="body"), np.nan)   # (8.48)
psi_f = np.where(out, ch08.stokes_sphere_streamfunction(Rr, Th, 1.0, 1.0, frame="fluid"), np.nan)  # N93
ur, ut = ch08.stokes_sphere_velocity(Rr, Th, 1.0, 1.0)    # (8.49) body frame
speed = np.where(out, np.hypot(ur, ut), np.nan)           # |u|/U
psi_i = np.where(out, PF.sphere(1.0, 1.0).psi(np.abs(Y), X) - 0.5*Y**2, np.nan)   # Ch. 6 ideal sphere, fluid frame
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
im = axs[0].pcolormesh(X, Y, speed, cmap="Blues", vmin=0, vmax=1.1, shading="auto")
axs[0].contour(X, Y, psi_b, np.linspace(-4, 4, 17), colors=COLORS["ink"], linewidths=0.8)
fig.colorbar(im, ax=axs[0], label="|u|/U")
axs[1].contour(X, Y, psi_f, np.linspace(-1.2, 1.2, 25), colors=COLORS["blue"], linewidths=0.9)
axs[2].contour(X, Y, psi_i, np.linspace(-0.5, 0.5, 21), colors=COLORS["muted"], linewidths=0.9)
for ax, title in zip(axs, ("Body frame: stream past the sphere", "Fluid frame: sphere moving left (Stokes)", "Fluid frame: ideal flow (Ch. 6)")):
    sphere(ax, 1.0); ax.set_xlim(-6, 6); ax.set_ylim(-4, 4); ax.set_xlabel("x/a"); ax.set_title(title, fontsize=10)
axs[0].set_ylabel("y/a")
fig.suptitle("Creeping flow: symmetric, and felt far away", fontsize=11)
savefig(fig, "ch08", "c13_stokes_sphere"); plt.show()
""", see=r"""(a) the stream parting round the sphere, slowed over a wide region (pale blue); (b) closed loops that look the same in
front of and behind the moving sphere; (c) Ch. 6's ideal flow in the same frame, a much more compact disturbance — our remakes
of Figs. 8.17 (sphere) and 8.19 `N107`.""",
    read=r"""Count the radii where the loops still bend: in (b) the disturbance is visible beyond 5a (it decays like a/r), in (c) it
is nearly gone by 3a (like a³/r³).""",
    change=r"""…Re were 1: inertia breaks the symmetry and a wake forms behind (C15).""")
remind([
    ("advecting many tracers in one solve_ivp call", "stack all tracer positions into one long state vector and let `solve_ivp` move them together (Ch. 5 P137)."),
])
nb.animation(r"""
nfr = 60 if not FAST else 30                              # frames
T_end = 10.0                                              # the sphere travels 10 radii (a = 1, U = 1, time in a/U)
y0s = np.linspace(0.25, 3.0, 10)                          # ten tracers on the line x = 0, above the axis
z0 = np.r_[np.zeros(10), y0s]                              # state: all x's then all y's


def rhs_stokes(t_, z):                                    # fluid-frame Stokes velocity; the sphere centre is at x = 5 − t
    x_, y_ = z[:10] - (5.0 - t_), z[10:]                  # positions relative to the sphere
    u, v, _ = ch08.stokes_sphere_velocity_xyz(x_, y_, 0*x_, 1.0, 1.0, frame="fluid")   # (8.49) minus the stream
    return np.r_[np.nan_to_num(u), np.nan_to_num(v)]


ideal = PF.sphere(1.0, 1.0)                               # Ch. 6's ideal-flow sphere (stream along its z-axis)


def rhs_ideal(t_, z):                                     # the same for ideal flow
    x_, y_ = z[:10] - (5.0 - t_), z[10:]
    uR, uz = ideal.velocity_cyl(np.abs(y_), x_)           # body-frame (u_R, u_z) with z = our x
    return np.r_[np.nan_to_num(uz - 1.0), np.nan_to_num(uR*np.sign(y_))]   # subtract the stream: fluid frame


tf = np.linspace(0, T_end, nfr)                           # frame times
S_ = integrate.solve_ivp(rhs_stokes, (0, T_end), z0, t_eval=tf, rtol=1e-6, atol=1e-8)   # all Stokes tracers at once
I_ = integrate.solve_ivp(rhs_ideal, (0, T_end), z0, t_eval=tf, rtol=1e-6, atol=1e-8)    # all ideal-flow tracers
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.4), sharey=True)
pts, balls, trails = [], [], []
for ax, sol, c, title in ((a1, S_, COLORS["blue"], "Stokes flow"), (a2, I_, COLORS["muted"], "ideal flow")):
    pts.append(ax.plot(sol.y[:10, 0], sol.y[10:, 0], "o", color=c, ms=5)[0])
    trails.append([ax.plot([], [], color=c, lw=0.8)[0] for _ in range(10)])
    balls.append(plt.Circle((5.0, 0.0), 1.0, color=COLORS["muted"], alpha=0.6)); ax.add_patch(balls[-1])
    ax.set_xlim(-6, 6); ax.set_ylim(-0.5, 3.5); ax.set_aspect("equal"); ax.set_xlabel("x/a"); ax.set_title(title, fontsize=10)
a1.set_ylabel("y/a")


def update(i):                                            # frame i
    for sol, p_, tr, ball in zip((S_, I_), pts, trails, balls):
        p_.set_data(sol.y[:10, i], sol.y[10:, i])         # tracers now
        for k in range(10):
            tr[k].set_data(sol.y[k, :i+1], sol.y[10 + k, :i+1])   # their paths so far
        ball.center = (5.0 - tf[i], 0.0)                  # the sphere moves left at U
    return []


show_animation(animate(update, frames=nfr, fig=fig, interval=70))
""")
see_read_change(r"""Ten tracers on a vertical line as a sphere passes from right to left: in Stokes flow they are pushed aside and
forward, then pulled back and end up only slightly ahead of where they started; in ideal flow they barely move — the `A4`
animation.""",
                r"""Each tracer's path is a loop: forward (left) while the sphere approaches, sideways as it passes, back as it leaves. The
Stokes loops are large even three radii away — the slow a/r decay.""",
                r"""…the sphere stopped halfway: the tracers would stop at once — there is no inertia to coast on (reversibility).""")
explainer("stokes_sphere_flow", "What does creeping flow round a sphere look like?", r"""
Toggle body and fluid frames and Stokes, Oseen and ideal flow, and drag Re: symmetric loops, a disturbance that decays only
like a/r, and then a wake that grows from far downstream once r reaches a/Re — three comparisons that need motion and switching
to be believed. (Its Oseen mode previews C15.)""",
          ["Fluid frame, Stokes: follow one tracer — it goes round a loop and comes back.",
           "Read the side-line speed view: at 10a Stokes still shows a 7.5 % deficit, ideal flow none.",
           "Switch to Oseen and raise Re to 1: where does the circle r = a/Re_a sit, and where does the wake begin?",
           "Open the Derivation tab at D29 step 8: zoom out to see why the r⁴ term must go."])
whatif(r"""
…we asked for the force? The pressure and the shear stress on the surface follow from (8.49); integrating them gives the drag
(C14).
""")

nb.recap("R20", "Dimensional analysis without ρ", r"""
If inertia does not matter, ρ cannot appear: D = f(μ, U, a); four variables and three dimensions leave one group,
D/(μUa) = constant, so D ∝ μUa and $C_D\propto1/\mathrm{Re}$ before any flow calculation (Ch. 1's Π method; Exercise 4.60).
""", where="Ch. 1 §1.11")
nb.code(r"""
print(ch01.pi_groups({"D": "N", "mu": "Pa*s", "U": "m/s", "a": "m"}))   # one group: D μ⁻¹ U⁻¹ a⁻¹
""")
core("C14", "Stokes drag and settling", r"""
What force does the fluid exert on the sphere — and where does it come from, pressure or friction?
""", eqs=("8.51",))
problem(r"""
A cloud droplet falls at its terminal speed when the drag equals its weight minus buoyancy; Millikan used exactly this to weigh
single oil drops and count electrons; sediment settles in rivers and the ocean; aerosols stay aloft for days. All need the drag
of a slowly moving sphere — and the surprise is how it splits.
""")
idea(r"""
front (θ = π): high pressure +3μU/2a  ─►  ● ─►  rear (θ = 0): low pressure −3μU/2a
sides: friction σ_rθ = −(3μU/2a) sin θ drags the surface along
drag = ⅓ pressure (2πμaU) + ⅔ friction (4πμaU) = 6πμaU
""", words=r"The traction on the surface has a pressure part and a friction part; their x-components add up to a *uniform* push "
          r"3μU/2a over the whole sphere.")
remind([
    ("line integral of a gradient", "if ∂p/∂r and (1/r)∂p/∂θ are both known, integrate one and check the other: they must come from one function p (Ch. 1 P35)."),
    ("surface integrals on a sphere, dA = 2πa² sin θ dθ", "for an axisymmetric integrand the sphere is a stack of rings of area 2πa² sin θ dθ (Ch. 6 P163)."),
])
gloss("Traction", r"""The force per unit area that the fluid exerts on the sphere is $\mathbf t=\boldsymbol\sigma\cdot\mathbf e_r$
(Cauchy, Ch. 2 §2.6); its x-component is $\sigma_{rr}\cos\theta-\sigma_{r\theta}\sin\theta$ with θ measured from +x.""")
D("D30", ref="8.50", check_src=r"""
ur, uth, p = S["u_r"], S["u_theta"], S["p"]              # (8.49) and (8.50) from the cached engine
L_ = CL.vector_laplacian([ur, uth, 0], "spherical", X)   # ∇²u in spherical components (core.curvilinear)
print("p − p∞ =", p)                                     # −3μaU cosθ/(2r²)
assert sp.simplify(sp.diff(p, r) - mu*L_[0]) == 0        # radial component of ∇p = μ∇²u
assert sp.simplify(sp.diff(p, th)/r - mu*L_[1]) == 0     # polar component: consistent (D30 step 9)
print("rear (θ = 0):", p.subs({r: a, th: 0}), "  front (θ = π):", p.subs({r: a, th: sp.pi}))   # −3μU/2a and +3μU/2a
""")
note("N88 [B], N89 [B] · Pressure and surface stresses.", r"""
**N88** The pressure (8.50) has its maximum +3μU/2a at the front stagnation point (θ = π) and its minimum **−**3μU/2a at the rear
(θ = 0). ⚠️ The book prints the minimum without the minus sign (and its figure labels both extremes 1.5). **N89** Surface stresses
(ours): on r = a the viscous normal stress vanishes, $\sigma_{rr}=-p$, and $\sigma_{r\theta}=-\frac{3\mu U}{2a}\sin\theta$.
""", equation=EQ["8.50"], ref="8.50")
D("D31", ref="8.51", check_src=r"""
Sm = CL.strain_rate([ur, uth, 0], "spherical", X)        # strain-rate tensor of (8.49) in spherical components
srr = sp.simplify((-p + 2*mu*Sm[0, 0]).subs(r, a))       # σ_rr on the sphere (steps 3–4, 8)
srt = sp.simplify((2*mu*Sm[0, 1]).subs(r, a))            # σ_rθ on the sphere (steps 5–7)
print("σ_rr(a) =", srr, "  σ_rθ(a) =", srt, "  viscous part of σ_rr(a):", S["sigma_rr_viscous_a"])
Dp = sp.integrate(srr*sp.cos(th)*2*sp.pi*a**2*sp.sin(th), (th, 0, sp.pi))     # pressure part of the drag (step 12)
Df = sp.integrate(-srt*sp.sin(th)*2*sp.pi*a**2*sp.sin(th), (th, 0, sp.pi))   # friction part (step 13)
print("pressure:", sp.simplify(Dp), "  friction:", sp.simplify(Df))           # 2πμaU and 4πμaU
assert sp.simplify(Dp + Df - 6*sp.pi*mu*a*U) == 0        # (8.51): D = 6πμaU
""")
nb.worked_example("the drag and fall speed of a 10 µm cloud droplet", r"""
a = 10 µm, water droplet ρ′ = 1000 kg/m³ in air ρ = 1.2 kg/m³, μ = 1.81 × 10⁻⁵ Pa s, g = 9.81 m/s².
1. Terminal balance $\frac43\pi a^3g(\rho'-\rho)=6\pi\mu aU$ ⇒ $U_t=\frac{2(\rho'-\rho)ga^2}{9\mu}=\frac{2\times998.8\times9.81\times
   10^{-10}}{9\times1.81\times10^{-5}}=1.20$ cm/s.
2. $\mathrm{Re}=2aU_t\rho/\mu=2\times10^{-5}\times0.012\times1.2/1.81\times10^{-5}=0.016$ ✓ ≪ 1.
3. $D=6\pi\mu aU_t=4.1\times10^{-11}$ N: 1.4 × 10⁻¹¹ N from pressure, 2.7 × 10⁻¹¹ N from friction.
4. $C_D=24/\mathrm{Re}\approx1.5\times10^3$.
5. A 100 µm drizzle drop: 100× faster by a² — 1.2 m/s — but then Re ≈ 16 and Stokes no longer applies.
""")
remind([
    ("terminal velocity, effective weight, buoyancy", "a body falls at constant speed when drag = weight − buoyancy = (4/3)πa³(ρ′ − ρ)g (Ch. 1 P24, weight)."),
    ("Gauss–Legendre quadrature", "`np.polynomial.legendre.leggauss(n)` gives n nodes and weights that integrate polynomials of degree < 2n exactly (Ch. 5 P143)."),
])
nb.code(r"""
st = ch08.settling_state(10e-6, 1000.0, 1.2, 1.81e-5)    # droplet radius [m], ρ′, ρ [kg/m³], μ of air [Pa s] (g = 9.80665)
print({k: (f"{st[k]:.4g}" if not isinstance(st[k], bool) else st[k]) for k in ("U_t", "Re", "D", "C_D", "valid", "D_pressure", "D_friction")})
print(ch08.stokes_drag(1.81e-5, 10e-6, st["U_t"], parts=True))   # (8.51) split into its two parts [N]
srr_, srt_, tx = ch08.stokes_sphere_surface_stresses(np.array([0, np.pi/2, np.pi]), U=st["U_t"], a=10e-6, mu=1.81e-5)
print("x-traction at the rear, side, front:", tx, "Pa   (3μU/2a =", 1.5*1.81e-5*st["U_t"]/10e-6, ")")
print("drag collected from the rear to the side:", ch08.stokes_drag_running(np.pi/2, 1.81e-5, 10e-6, st["U_t"]))
print("C_D:", ch08.stokes_drag_coefficient(st["Re"]), SIM.sphere_drag_coefficient(st["Re"], "stokes"))   # (8.52), and Ch. 4's
""", explain=r"""
1. `settling_state` solves the terminal balance (N90 below) and checks Re < 0.1: U_t = 1.20 cm/s, Re = 0.016, D = 4.10 × 10⁻¹¹ N,
   C_D = 1505, valid.
2. `stokes_drag(..., parts=True)` splits (8.51): one third pressure, two thirds friction.
3. The surface tractions' x-component is the same everywhere, 3μU/2a = 0.0326 Pa — the pressure part and the friction part
   trade off round the sphere.
4. `stokes_drag_running` collects the drag from the rear stagnation point to θ: by the side (θ = π/2) half of each part is in.
5. (8.52) agrees with Ch. 4's drag-coefficient function in its Stokes branch.
""")
nb.check_agree(r"""
M_ = 2000                                                # midpoint sum over 2000 θ-cells
thm = (np.arange(M_) + 0.5)*np.pi/M_; dth = np.pi/M_      # cell midpoints and width [rad]
a_, U_, mu_ = 10e-6, st["U_t"], 1.81e-5                  # the droplet
srr_, srt_, _ = ch08.stokes_sphere_surface_stresses(thm, U=U_, a=a_, mu=mu_)   # σ_rr = −p and σ_rθ on the sphere
dA = 2*np.pi*a_**2*np.sin(thm)*dth                        # ring areas [m²]
Dp_ = np.sum(srr_*np.cos(thm)*dA); Df_ = np.sum(-srt_*np.sin(thm)*dA)   # x-projections, summed
D_ = ch08.stokes_drag(mu_, a_, U_, parts=True)
print(f"midpoint: pressure {Dp_:.6e} N, friction {Df_:.6e} N")
assert np.allclose([Dp_, Df_], [D_["pressure"], D_["friction"]], rtol=1e-6)
tx_fn = lambda th_: ch08.stokes_sphere_surface_stresses(th_, U=U_, a=a_, mu=mu_)[2]   # x-traction as a function of θ
assert np.isclose(ch08.sphere_drag_quadrature(tx_fn, a_, n=16), D_["total"], rtol=1e-12)   # Gauss–Legendre: exact already
""")
note("N90 [B], N91 [B] · Terminal velocity and Millikan's oil drops.", r"""
**N90** The balance $(4/3)\pi a^3g(\rho'-\rho)=6\pi\mu aU$ gives $U_t=\frac{2(\rho'-\rho)ga^2}{9\mu}$ — ∝ a²: halve the droplet,
quarter the speed. This sets how long aerosols and cloud droplets stay aloft and how fast silt settles (Ch. 13). **N91**
Millikan's experiment (Fig. 8.18): with the plates uncharged the fall speed gives the radius; with the upper plate negative the
drop rises at U_u: $6\pi\mu U_ua+(4/3)\pi a^3g(\rho'-\rho)=neE$, E = −V_b/L; charges from many drops are whole multiples of e.
""")
nb.code(r"""
m = ch08.synthetic_millikan(40, seed=0)                  # 40 synthetic oil drops with 1 % noise in the measured speeds
print(f"estimated e = {m['e_est']:.6e} C, relative error vs CODATA {m['e_est']/1.602176634e-19 - 1:+.1e}")
print("charges in units of e (first 10):", np.round(m["q"][:10]/m["e_est"], 2))   # whole numbers
""")
note("N92 [B] · Drag coefficient.", r"""
With the frontal area πa² in $C_D\equiv\frac{F_D}{\tfrac12\rho U^2A}$ *(4.107)*, Stokes' law gives (8.52), with Re = 2aU/ν based on the
**diameter** — the low-Re branch of the drag curve of Ch. 4.
""", equation=EQ["8.52"], ref="8.52")
nb.figure(r"""
import warnings                                           # to silence the "Re > 0.1" warnings of the dashed extensions
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13, 3.8))
th_ = np.linspace(0, 2*np.pi, 13)[:-1]                    # 12 points round the circle in the drawing plane
for t0 in th_:                                            # traction arrows at each point (units μU/a, U = a = μ = 1)
    c, s = np.cos(t0), np.sin(t0)                         # the point (and the outward normal e_r) on the circle
    th_s = np.arccos(c)                                   # the spherical angle θ from the downstream +x axis
    srr1, srt1, _ = ch08.stokes_sphere_surface_stresses(th_s, U=1.0, a=1.0, mu=1.0)   # σ_rr = −p and σ_rθ there
    e_th = np.array([-abs(s), np.sign(s)*c])              # unit vector of increasing θ in the drawing plane
    t_p = srr1*np.array([c, s])                           # pressure (normal) part of the traction
    t_f = srt1*e_th                                       # friction (tangential) part
    for vec, col, lw in ((t_p, COLORS["orange"], 1.0), (t_f, COLORS["rose"], 1.0), (t_p + t_f, COLORS["ink"], 1.6)):
        a1.annotate("", xy=(1.1*c + 0.5*vec[0], 1.1*s + 0.5*vec[1]), xytext=(1.1*c, 1.1*s),
                    arrowprops=dict(arrowstyle="->", color=col, lw=lw))
sphere(a1, 1.0); a1.set_xlim(-2.4, 2.8); a1.set_ylim(-2.2, 2.2); a1.axis("off")
a1.set_title("Tractions: pressure (orange), friction (rose),\nsum (black) = 3μU/2a in x everywhere", fontsize=9)
thv = np.linspace(0, np.pi, 200)                          # from the rear (0) to the front (π)
a2.plot(thv, ch08.stokes_sphere_pressure(1.0, thv, 1.0, 1.0, 1.0), color=COLORS["orange"], lw=2, label="p − p∞ [μU/a]")
a2.plot(thv, ch08.stokes_sphere_surface_stresses(thv, 1.0, 1.0, 1.0)[1], color=COLORS["rose"], lw=2, label="σ_rθ [μU/a]")
b2 = a2.twinx()                                           # running drag on a second axis
run_ = ch08.stokes_drag_running(thv, 1.0, 1.0, 1.0)       # pressure and friction parts collected up to θ [μaU]
b2.plot(thv, run_["pressure"]/np.pi, color=COLORS["orange"], ls="--", lw=1.2)
b2.plot(thv, run_["friction"]/np.pi, color=COLORS["rose"], ls="--", lw=1.2)
b2.set_ylabel("running drag [π μaU] (dashed)"); b2.set_ylim(0, 4.4)
a2.set_xlabel("θ from the rear [rad]"); a2.legend(fontsize=7, loc="lower left"); a2.set_title("Stresses and the drag they build: 2π + 4π", fontsize=9)
radii = np.geomspace(1e-6, 1e-3, 80)                      # particle radii [m]
fams = [("water droplets in air", 1000.0, 1.2, 1.81e-5, COLORS["blue"]), ("quartz sand in water", 2650.0, 1000.0, 1e-3, COLORS["amber"]),
        ("bacteria in water", 1100.0, 1000.0, 1e-3, COLORS["teal"])]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")                       # beyond Re = 0.1 the law is drawn dashed, not trusted
    for name, rp, rf, mf, c in fams:
        Ut = ch08.terminal_velocity(radii, rp, rf, mf, warn=False)   # N90 [m/s]
        Re_ = 2*radii*Ut*rf/mf                            # Reynolds number with the diameter
        ok = Re_ < 0.1
        a3.loglog(1e6*radii[ok], Ut[ok], color=c, lw=2, label=name)
        a3.loglog(1e6*radii[~ok], Ut[~ok], color=c, lw=1, ls="--")
a3.set_xlabel("radius a [µm]"); a3.set_ylabel("U_t [m/s]"); a3.legend(fontsize=7); a3.set_title("Fall speed ∝ a² (dashed: Re > 0.1)", fontsize=9)
fig.suptitle("A third pressure, two thirds friction; fall speed ∝ a²", fontsize=11)
savefig(fig, "ch08", "c14_stokes_drag"); plt.show()
""", see=r"""(a) arrows round the sphere: pressure (orange) pushing at the front and pulling at the rear, friction (rose) along the
surface, and their sum (black) the same everywhere; (b) the pressure and shear stress against θ and their running integrals
(dashed) reaching 2π and 4π μaU; (c) straight lines of slope 2 for droplets, sand and bacteria, dashed where Re > 0.1 — our
remake of Fig. 8.17's stress plot (`N107`).""",
    read=r"""In (b) the orange dashed curve ends at 2 and the rose one at 4 (in units of πμaU): one third and two thirds of 6πμaU. In
(c) read a fall speed: a 10 µm droplet falls about 1 cm/s, a 10 µm sand grain about 0.2 mm/s.""",
    change=r"""…the droplet were oil (ρ′ = 900 kg/m³): $U_t\propto(\rho'-\rho)$, so it falls about 10 % more slowly.""")
explainer("stokes_drag_settling", "Where does 6πμaU come from?", r"""
Sweep around the sphere and watch the local pressure and shear tractions as arrows and curves while their running integrals
build up to 2πμaU and 4πμaU in the term bars; then drag the particle radius and read U_t and Re as the Stokes law's validity limit
is crossed.""",
          ["Play the sweep from the rear (θ = 0) to the front: which bar grows first?",
           "Pick 'drizzle 100 µm': the status turns amber — why?",
           "Switch the third view to C_D(Re) and find where Oseen and the correlation leave 24/Re.",
           "Click a point on the sphere and read the traction arithmetic."])
whatif(r"""
…we looked far from the sphere? The viscous force there decays like a/r³ but inertia only like a/r²; at some distance inertia
wins even for tiny Re (C15).
""")

core("C15", "Where Stokes fails, and Oseen's fix", r"""
If Re is tiny, how can inertia ever matter?
""")
nb.md(r"*In one line:* $\frac{\text{inertia}}{\text{viscous}}\sim\frac{\rho Ua}{\mu}\frac ra=\mathrm{Re}_a\frac ra$ — small near the "
      r"sphere, order one at a distance of order $a/\mathrm{Re}_a$.")
problem(r"""
Far behind a slowly falling bead the fluid is barely disturbed — so inertia should be even less important there. It is the
opposite: the viscous force falls off faster with distance than inertia does, so beyond some radius Stokes' solution is wrong.
For a cylinder it is so wrong that no Stokes solution exists at all.
""")
idea("", table=r"""
| at distance r | size |
|---|---|
| disturbance velocity | ~ Ua/r |
| inertia ρu·∇u | ~ ρU²a/r² |
| viscous μ∇²u | ~ μUa/r³ |
| ratio | ~ (ρUa/μ)(r/a) = Re_a r/a → order one at r ~ a/Re_a |
""", words=r"The slowly decaying disturbance is swept along by the stream; its viscous force falls one power of r faster.")
remind([
    ("asymptotic size of a term far away", "for r ≫ a keep only the slowest-decaying piece of each quantity; a derivative of a 1/r term brings another 1/r (Ch. 4 P130, order-of-magnitude scaling)."),
    ("series of 1 − e^(−s)", "$1-e^{-s}=s-\\frac{s^2}{2}+\\ldots$ for small s; `sp.series` produces it (Ch. 4 P117)."),
])
D("D32")
note("N94 [C] · A singular perturbation at infinity.", r"""
Treated as the first term of an expansion in Re, Stokes' solution is not uniformly valid: the O(Re) correction grows without
bound relative to it as r → ∞ (Whitehead's paradox); for a cylinder the Stokes equations cannot even meet the uniform stream
(Stokes' paradox, Exercise 8.37). Unlike the boundary layer (Ch. 9), where 1/Re multiplies the highest derivative near the wall,
here the trouble is far away.
""")
note("N95 [B], N96 [B] · Oseen's linearisation (1910).", r"""
**N95** Write u = U + u′ with u′ small far away. The advective term splits as
$u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}+w\frac{\partial u}{\partial z}=U\frac{\partial u'}{\partial x}+
\big[u'\frac{\partial u'}{\partial x}+v'\frac{\partial u'}{\partial y}+w'\frac{\partial u'}{\partial z}\big]$; dropping the bracket gives
the equation below. ⚠️ The book prints $+\partial p/\partial x_i$; the minus of (8.43) is needed. **N96** Conditions: u′, v′, w′ → 0
far away; u′ = −U, v′ = w′ = 0 on the sphere.
""", equation=r"\rho U\frac{\partial u_i'}{\partial x}=-\frac{\partial p}{\partial x_i}+\mu\nabla^2u_i'")
nb.code(r"""
ol_ = ch08.oseen_linearisation_sympy()                   # N95 in sympy
print("dropped (quadratic) part:", ol_["dropped"])
print("Oseen x-equation:", ol_["oseen_x"])               # ρU ∂u'/∂x = μ∇²u' − ∂p/∂x
""")
note("N97 [B], N98 [B] · Oseen's stream function.", r"""
**N97** Oseen's solution is (8.53), with Re = 2aU/ν (stated; no slip holds only to O(Re)). **N98** Near the sphere the
exponential's series returns Stokes' (8.48). ⚠️ The text calls it "(9.68)"; it means (8.48).
""", equation=EQ["8.53"], ref="8.53")
D("D33", ref="8.53")
nb.code(r"""
ol = ch08.oseen_limit_sympy()                            # D33 in sympy: expand 1 − e^(−s) and let Re → 0
print("s =", ol["s"])
print("Oseen − Stokes as Re → 0:", sp.simplify(ol["difference"]))   # 0
print("first-order correction:", ol["first_order"])       # O(Re) and growing with r: the far-field change
""")
nb.worked_example("where inertia catches up", r"""
Re = 2aU/ν = 0.02, so $\mathrm{Re}_a=\mathrm{Re}/2=0.01$.
1. The estimate $\mathrm{Re}_a\,r/a\sim1$ puts the crossover at r of order 100a. With the exact prefactor ½ (the code below computes
   it on the axis and the side line) inertia equals friction at r ≈ 2a/Re_a = 200a.
2. For the 10 µm droplet (Re = 0.016, Re_a = 0.008) the same estimate gives r of order 125a = 1.25 mm, and about 250a with the
   prefactor.
3. Oseen C_D at Re = 0.5: $\frac{24}{\mathrm{Re}}\big(1+\frac{3}{16}\mathrm{Re}\big)=48(1+3/32)=52.5$ vs Stokes 48 — a 9 % correction.
""")
nb.code(r"""
Re_a = 1e-3*1e-3/1e-4                                    # U a/ν with U = 1 mm/s, a = 1 mm, ν = 10⁻⁴ m²/s: 0.01
for th_, name in ((0.0, "axis θ = 0"), (np.pi/4, "θ = π/4"), (np.pi/2, "side line θ = π/2")):
    for rr_ in (10.0, 100.0, 400.0):                     # distances in radii
        ratio = ch08.inertia_viscous_ratio(rr_*1e-3, th_, U=1e-3, a=1e-3, nu=1e-4)   # |u·∇u| / |ν∇²u| on Stokes' field
        print(f"{name:18s} r = {rr_:5.0f}a: ratio = {ratio:.4f}, ratio ÷ (Re_a r/a) = {ratio/(Re_a*rr_):.3f}")
print("C_D at Re = 0.5: Oseen", ch08.oseen_drag_coefficient(0.5), " Stokes", ch08.stokes_drag_coefficient(0.5),
      " Proudman–Pearson", round(ch08.proudman_pearson_drag_coefficient(0.5), 3), " Morrison (data fit)", round(SIM.sphere_drag_coefficient(0.5), 3))
print("ψ Oseen (Re = 1e-6) vs Stokes at r = 2a, θ = π/3:", ch08.oseen_streamfunction(2e-3, np.pi/3, U=1e-3, a=1e-3, Re=1e-6),
      ch08.stokes_sphere_streamfunction(2e-3, np.pi/3, U=1e-3, a=1e-3))
""", explain=r"""
1. `inertia_viscous_ratio` evaluates |u·∇u| / |ν∇²u| on Stokes' field by finite differences. It grows in proportion to r, and
   ratio ÷ (Re_a r/a) tends to **½ on the axis and on the side line** — so inertia equals friction at r ≈ 2a/Re_a (200a here). At
   θ = π/4 the leading inertia term nearly cancels (it changes sign near θ ≈ 55°) and the prefactor is smaller, about 0.16. The
   book's r ~ a/Re is the right order of magnitude.
2. The drag laws at Re = 0.5: Oseen 52.5, Stokes 48, Proudman–Pearson and the Morrison correlation of the measurements close to
   Oseen.
3. At a vanishing Re Oseen's ψ equals Stokes' to about 10⁻⁶ (D33).
""")
nb.check_agree(r"""
hh = 1e-4                                                # relative finite-difference step
pts = [(r_*1e-3, 0.0) for r_ in (50.0, 100.0, 200.0)]   # three points on the axis behind the sphere (x, y) [m]


def vel(x_, y_):                                          # Stokes body-frame velocity (u, v) in the plane z = 0
    u_, v_, _ = ch08.stokes_sphere_velocity_xyz(np.array(x_), np.array(y_), np.array(0.0), 1e-3, 1e-3)
    return np.array([u_, v_], float)


mine = []
for x0, y0 in pts:
    d = hh*x0                                             # step [m]
    ux = (vel(x0 + d, y0) - vel(x0 - d, y0))/(2*d); uy = (vel(x0, y0 + d) - vel(x0, y0 - d))/(2*d)   # ∂u/∂x, ∂u/∂y
    uxx = (vel(x0 + d, y0) - 2*vel(x0, y0) + vel(x0 - d, y0))/d**2
    uyy = (vel(x0, y0 + d) - 2*vel(x0, y0) + vel(x0, y0 - d))/d**2
    lap = uxx + 2*uyy                                     # ∇²u on the axis: the z-direction curves like the y-direction (axisymmetry)
    adv = vel(x0, y0)[0]*ux + vel(x0, y0)[1]*uy           # u·∇u (w = 0 in the plane z = 0)
    mine.append(np.linalg.norm(adv)/np.linalg.norm(1e-4*lap))   # |u·∇u| / |ν∇²u|
lib = [ch08.inertia_viscous_ratio(x0, 0.0, U=1e-3, a=1e-3, nu=1e-4) for x0, _ in pts]
print("by hand:", np.round(mine, 4), " library:", np.round(lib, 4))
assert np.allclose(mine, lib, rtol=1e-3)                 # the same ratio
rs = np.geomspace(50, 500, 8)                            # r/a from 50 to 500
slope = observed_order(rs, [ch08.inertia_viscous_ratio(r_*1e-3, np.pi/2, U=1e-3, a=1e-3, nu=1e-4) for r_ in rs])   # log–log slope
print(f"log–log slope of the ratio for r/a in [50, 500]: {slope:.3f}")
assert abs(slope - 1) < 0.05                             # D32: the ratio grows like r
""")
nb.figure(r"""
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13, 3.8))
rs = np.geomspace(1.5, 3000, 60)                          # r/a
for Rea, c in ((1e-3, COLORS["blue"]), (1e-2, COLORS["teal"]), (1e-1, COLORS["orange"])):   # three Re_a
    nu_ = 1e-3*1e-3/Rea                                    # ν that gives this Re_a with U = 1 mm/s, a = 1 mm
    a1.loglog(rs, ch08.inertia_viscous_ratio(rs*1e-3, np.pi/2, U=1e-3, a=1e-3, nu=nu_), color=c, lw=2, label=f"Re_a = {Rea:g}")
a1.axhline(1, color=COLORS["muted"], ls="--", lw=1); a1.text(2, 1.3, "inertia = friction", fontsize=8)
a1.set_xlabel("r/a"); a1.set_ylabel("|u·∇u| / |ν∇²u| (side line)"); a1.legend(fontsize=7); a1.set_title("Inertia returns far away", fontsize=10)
n_ = 200 if not FAST else 100                             # contour grid
xg = np.linspace(-15, 15, n_); yg = np.linspace(-10, 10, n_)
Xo, Yo = np.meshgrid(xg, yg); Ro = np.hypot(Xo, Yo); To = np.arccos(np.clip(Xo/np.maximum(Ro, 1e-12), -1, 1))
outo = Ro > 1
ps_ = np.where(outo, ch08.stokes_sphere_streamfunction(Ro, To, 1.0, 1.0, frame="fluid"), np.nan)
po_ = np.where(outo, ch08.oseen_streamfunction(Ro, To, 1.0, 1.0, Re=1.0, frame="fluid"), np.nan)
lev = np.linspace(-4, 4, 33)
a2.contour(Xo, Yo, ps_, lev, colors=COLORS["muted"], linewidths=0.7, linestyles="--")
a2.contour(Xo, Yo, po_, lev, colors=COLORS["blue"], linewidths=1.0)
sphere(a2, 1.0); a2.set_xlim(-15, 15); a2.set_ylim(-10, 10); a2.set_xlabel("x/a (sphere moving left)")
a2.set_title("Fluid frame, Re = 1: Oseen (blue) vs Stokes (dashed)", fontsize=9)
Res = np.geomspace(0.01, 10, 100)
a3.semilogx(Res, ch08.stokes_drag_coefficient(Res)*Res/24, color=COLORS["muted"], lw=2, label="Stokes 24/Re")
a3.semilogx(Res, ch08.oseen_drag_coefficient(Res)*Res/24, color=COLORS["blue"], lw=2, label="Oseen (1 + 3Re/16)")
a3.semilogx(Res, np.array([ch08.proudman_pearson_drag_coefficient(r_) for r_ in Res])*Res/24, color=COLORS["accent"], lw=1.5, ls="--", label="Proudman–Pearson")
a3.semilogx(Res, np.array([SIM.sphere_drag_coefficient(r_) for r_ in Res])*Res/24, color=COLORS["amber"], lw=2, label="Morrison (measurements)")
a3.set_ylim(0.9, 3); a3.set_xlabel("Re = 2aU/ν"); a3.set_ylabel("C_D·Re/24 [–]"); a3.legend(fontsize=7); a3.set_title("Drag laws", fontsize=10)
fig.suptitle("Stokes is right near the sphere and wrong far away", fontsize=11)
savefig(fig, "ch08", "c15_oseen"); plt.show()
""", see=r"""(a) the inertia/friction ratio growing in proportion to r for three Re_a, crossing 1 at r ≈ 2a/Re_a; (b) at Re = 1 the
Oseen streamlines (blue) are crowded behind the moving sphere into a wake while Stokes' (dashed) stay symmetric — our remake of
Fig. 8.20 `N100`; (c) C_D·Re/24: Stokes flat at 1, Oseen rising, the Morrison correlation of measurements between them below
Re ≈ 5, and Proudman–Pearson `N101` following Oseen at first.""",
    read=r"""In (c) any curve above 1 means more drag than Stokes' law: at Re = 1 Oseen adds 19 %; the measurements (amber) lie between
Stokes and Oseen up to Re ≈ 5, the range where both formulas are "fairly accurate".""",
    change=r"""…Re doubled: the crossover distance in (a) halves and the wake in (b) widens and moves closer to the sphere.""")
note("N99 [B], N100 [B], N101 [C] · Drag beyond Stokes.", r"""
**N99** Oseen's drag coefficient is below; ⚠️ with the radius-based $\mathrm{Re}_a=\mathrm{Re}/2$ the same coefficient reads 3/8.
**N100** The Oseen streamlines (panel b) are asymmetric, with a wake behind. **N101** Matched asymptotic expansions (Kaplun,
Proudman & Pearson 1957) continue the series: $D=6\pi\mu aU\big(1+\frac38\mathrm{Re}_a+\frac9{40}\mathrm{Re}_a^2\ln\mathrm{Re}_a+
\ldots\big)$ (named only).
""", equation=r"C_D=\frac{24}{\mathrm{Re}}\Big(1+\frac3{16}\mathrm{Re}\Big)")
nb.plotly(r"""
n_ = 200 if not FAST else 100                             # contour grid
xg = np.linspace(-15, 15, n_); yg = np.linspace(-10, 10, n_)
Xo, Yo = np.meshgrid(xg, yg); Ro = np.hypot(Xo, Yo); To = np.arccos(np.clip(Xo/np.maximum(Ro, 1e-12), -1, 1))
lev = np.linspace(-4, 4, 25)                              # streamline levels


def polylines(psi_):                                      # contour lines of ψ as one (x, y) polyline with NaN breaks
    cs = plt.figure().add_subplot().contour(Xo, Yo, psi_, lev)   # matplotlib finds the lines (figure discarded)
    xs, ys = [], []
    for path in cs.get_paths():
        for poly in path.to_polygons(closed_only=False):  # every separate piece of the line
            xs += list(poly[:, 0]) + [np.nan]; ys += list(poly[:, 1]) + [np.nan]
    plt.close("all")
    return np.array(xs), np.array(ys)


stokes_lines = polylines(np.where(Ro > 1, ch08.stokes_sphere_streamfunction(Ro, To, 1.0, 1.0, frame="fluid"), np.nan))
cache = {}                                                # Oseen lines for each Re, computed once
for Re_ in (0.1, 0.2, 0.5, 1.0, 2.0):
    cache[Re_] = polylines(np.where(Ro > 1, ch08.oseen_streamfunction(Ro, To, 1.0, 1.0, Re=Re_, frame="fluid"), np.nan))
circle = (np.cos(np.linspace(0, 2*np.pi, 60)), np.sin(np.linspace(0, 2*np.pi, 60)))   # the sphere's outline


def frame(Re_):                                           # the figure for one Re
    return {"Oseen streamlines": cache[Re_], "Stokes streamlines": stokes_lines, "sphere": circle}


Rs = [0.1, 0.2, 0.5, 1.0, 2.0]
fig = slider_figure(frame, "Re", Rs, unit="", xlabel="x/a (sphere moving left)", ylabel="y/a", title="", xrange=[-15, 15], yrange=[-10, 10])
step_titles(fig, [f"Re = {r_:g}: C_D Oseen {ch08.oseen_drag_coefficient(r_):.1f}, Stokes {ch08.stokes_drag_coefficient(r_):.1f}" for r_ in Rs])
recolor(fig, {"Oseen streamlines": COLORS["blue"], "Stokes streamlines": COLORS["muted"], "sphere": COLORS["ink"]},
        {"Stokes streamlines": "dot"})
fig.update_traces(line_width=1.2)
fig.show()
""", explain=r"""
Precomputed Oseen streamlines (fluid frame) for five Reynolds numbers over Stokes' symmetric ones (dotted); the title gives both
drag coefficients. Slide Re up and watch the wake form behind the sphere.
""")
nb.md(r"""
🎮 The explainer of C13 (`stokes_sphere_flow`, above) has an **Oseen mode**: open it again, choose Oseen and drag Re to see the
wake grow from far downstream.
""")

# =====================================================================================================================
# A.7 §8.7 — N102, S01, summary
# =====================================================================================================================
nb.section("8.7", "Final Remarks", intro=r"""
**What is this section about?** Where the exact solutions stop.
""")
note("N102 [C] · After the exact solutions.", r"""
Most laminar problems that can be solved with pencil and paper have been solved; the field moves on with perturbation methods
(flows close to a known one — Ch. 9 boundary layers, Ch. 11 stability) and with numerical solution of the Navier–Stokes equations
(Ch. 10, which uses this chapter's exact solutions as test cases).
""")
nb.pointer("**S01** Exercises 8.1–8.38 and the literature are not reproduced. Ideas from them used here: the Reynolds equation "
           "— Exercises 8.19–8.20 → C06; the power–dissipation identity — Exercise 8.12 → R11; the switched-on vortex — Exercise "
           "8.26 → C10; the stopped plate — Exercise 8.30 → C09; Couette start-up — Exercise 8.31 → R14; Hele-Shaw near an obstacle "
           "— Exercise 8.34 → C06; Stokes drag from the stresses — Exercise 8.35 → C14; Stokes' paradox for a cylinder — Exercise "
           "8.37 → C15; f = 64/Re — Exercise 8.7 → C03.")
nb.summary(
    clicked=[
        r"**C01** $\nu=\mu/\rho$ is a diffusivity: motion spreads $\sqrt{\nu t}$, and a flow stays laminar below Re ≈ 2000.",
        r"**C02** Between plates the profile is a Couette line plus a Poiseuille parabola, $u(y)=\frac Uh\,y-\frac1{2\mu}\frac{dp}{dx}\,y(h-y)$ *(8.5)*; the floor flow reverses when $dp/dx>2\mu U/h^2$.",
        r"**C03** In a pipe $Q=-\frac{\pi a^4}{8\mu}\frac{dp}{dz}$: halve the radius, lose 94 % of the flow; f = 64/Re.",
        r"**C04** Between rotating cylinders the swirl is $u_\varphi=AR+B/R$ *(8.9)*: rigid rotation plus a free vortex, chosen by the walls.",
        r"**C05** In a thin gap inertia is weighed with $\varepsilon^2\mathrm{Re}_L$, not Re_L; pressure is uniform across the gap.",
        r"**C06** Station by station the gap flow is Couette plus Poiseuille; integrating continuity gives the Reynolds equation $h_t+q_x=0$.",
        r"**C07** A constant flux through a narrowing gap needs a pressure hump — that hump carries the load, but only if the pad slides the right way.",
        r"**C08** A spreading layer obeys $h_t=(\rho g/3\mu)(h^3h_x)_x$, a diffusion that chokes itself as h shrinks.",
        r"**C09** With no imposed scale, y and t combine as $y/\sqrt{\nu t}$: one erfc curve for every time, $\delta_{99}=3.64\sqrt{\nu t}$ *(8.31)*.",
        r"**C10** Guess $At^{-n}F(\xi/\delta(t))$; matching powers of t and one conserved quantity fix n and δ.",
        r"**C11** An oscillating wall drives a decaying, lagging layer of depth $\sqrt{2\nu/\omega}$ — diffusion, not a wave.",
        r"**C12** At low Re rescale pressure by μU/L; the Stokes equations $\nabla p=\mu\nabla^2\mathbf u$ *(8.43)* are linear.",
        r"**C13** Stokes' sphere flow is fore–aft symmetric and decays only like a/r.",
        r"**C14** Drag $D=6\pi\mu aU$ *(8.51)*, ⅓ pressure and ⅔ friction; fall speed ∝ a².",
        r"**C15** Inertia returns at r ~ a/Re_a (about 2a/Re_a with the exact prefactor); Oseen's linearisation adds the wake and $C_D=\frac{24}{\mathrm{Re}}(1+\frac{3}{16}\mathrm{Re})$.",
    ],
    feeds_forward=[
        "Ch. 9: √(νt) ↔ √(νx/U), Blasius vs the temporal layer, the lubrication scaling as the boundary-layer approximation, similarity exponents.",
        "Ch. 10: exact solutions as code tests; Crank–Nicolson.",
        "Ch. 11: Couette, Poiseuille and Taylor–Couette base states; the Rayleigh criterion in A and B.",
        "Ch. 12: linear total stress, τ₀ and the friction velocity, f = 64/Re on the Moody chart.",
        "Ch. 13: Ekman layers with the (1 + i)/δ structure of the Stokes layer, spin-up times from √(νt), cyclostrophic balance, viscous gravity currents, settling of droplets and sediment.",
        "Ch. 16: arteries (Q ∝ a⁴), synovial lubrication, micro-swimmers (reversibility of Stokes flow).",
    ],
    left_out=[
        "The exercises (S01).",
        "Higher-order matched expansions for the sphere (named in N101).",
    ],
)

# ---------------------------------------------------------------------------------------------------------------------
# final pass: every equation named by number is written out; save (nbkit's coverage checks run here)
# ---------------------------------------------------------------------------------------------------------------------
if "--partial" in sys.argv:                                  # development aid: write what exists so far, unchecked
    import nbformat as _nbf
    finalize_equations()
    tidy_raw_tex()
    for _m in self_check_numbers() + self_check_prose():
        print("CHECK", _m)
    _nb2 = _nbf.v4.new_notebook(cells=list(nb.cells))
    _nb2.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    _out = ROOT / "notebooks" / "ch08_laminar_flow.ipynb"
    _nbf.write(_nb2, _out)
    print("wrote (partial, unchecked)", _out, len(nb.cells), "cells")
    sys.exit(0)

n_changed = finalize_equations()
n_tex = tidy_raw_tex()
print(f"plain-text exponents turned into maths in {n_tex} cells")
bad = self_check_numbers() + self_check_prose()
if bad:
    raise SystemExit("markdown cells cite equation numbers without maths:\n  " + "\n  ".join(bad))
miss = self_check_ledger()
if miss:
    raise SystemExit("Part E rows explained by a primer but no primer/reminder in the notebook:\n  " + "\n  ".join(miss))
out = nb.save()
print(f"wrote {out.relative_to(ROOT)} ({len(nb.cells)} cells; equations written out in {n_changed} cells)")
