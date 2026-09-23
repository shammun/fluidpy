"""Build the Chapter 3 teaching notebook: ``notebooks/ch03_kinematics.ipynb``.

Source of truth: ``analysis/ch03_design.md`` Part A (storyboard, one nbkit call per row), Part E (prerequisite ledger →
24 primers P87–P110, ch01/ch02 reminders, glosses), Part F (the 24 derivations D01–D24, one move per step) and
``analysis/ch03_curation.md`` (IDs, depths, section coverage). Physics lives in ``fluidpy.ch03_kinematics`` (imported as
``ch03``; it re-exports ``fluidpy/core/{kinematics,coords,vortices,transport}``); cells only call it. Pure drawing
comes from ``scripts/ch03_*.py``.

**Derivations are read from Part F at build time** (``part_f()`` below): goal, start, plan, tools, assumptions, every
step's *did / tex / why / plain*, result, check, meaning and traps are copied word for word, so the notebook and the
design can never drift apart. Explainer-only fields (*live*, *set*, *watch*) are dropped. Where Part F only sketches an
optional sympy check (D02, D17, D21) the builder adds a small commented one.

Labels shown to the reader: CORE blocks carry their id (`C05`), notes their id (`N38`), primers their number (`P87` …,
continuing ch02's P62–P86), derivations their id (`D06`).

Conventions (design header, binding): G[i, j] = ∂u_i/∂x_j; the book's R = G − Gᵀ (no ½), ω = ∇×u, element spin = ½ω;
γ = du₁/dx₂ = 2S₁₂ (Ch. 2's Γ was S₁₂), Γ = circulation; θ in §3.5 is the plane polar angle; RTT swept volume signed.

Run:  .venv/Scripts/python.exe notebooks/build_ch03.py
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch03")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# ---------------------------------------------------------------------------------------------------------------------
EQ = {
    "3.1": r"\mathbf u=d\mathbf r/dt,\ \mathbf a=d^2\mathbf r/dt^2",
    "3.2": r"F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)\ \text{at}\ \mathbf x=\mathbf r(t;\mathbf r_o,t_o)",
    "3.3": r"\frac{d}{dt}F[\mathbf r,t]=\frac{\partial F}{\partial r_i}\frac{dr_i}{dt}+\frac{\partial F}{\partial t}",
    "3.4": r"\frac{d}{dt}F[\mathbf r,t]=(\nabla F)\cdot\mathbf u+\frac{\partial F}{\partial t}\equiv\frac{DF}{Dt}",
    "3.5": r"\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F",
    "3.6": r"\frac{DF}{Dt}=\frac{\partial F}{\partial t}+|\mathbf u|\frac{\partial F}{\partial s}",
    "3.7": r"dx/u=dy/v=dz/w",
    "3.8": r"d\mathbf r/dt=\mathbf u(\mathbf r,t)",
    "3.9": r"\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'",
    "3.10": r"du_i=(\partial u_i/\partial x_j)\,dx_j",
    "3.11": r"\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}",
    "3.12": r"S_{ij}=\tfrac12\big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\big)",
    "3.13": r"R_{ij}=\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}",
    "3.14": r"\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}",
    "3.15": r"R_{ij}=-\varepsilon_{ijk}\omega_k",
    "3.16": r"\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}",
    "3.17": r"\boldsymbol\omega=0",
    "3.18": r"\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA",
    "3.19": r"du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i",
    "3.20": r"d\bar{\mathbf u}=\bar{\mathbf S}\cdot d\bar{\mathbf x}",
    "3.21": r"d\bar u_\alpha=\bar S_{\alpha\alpha}\,d\bar x_\alpha",
    "3.22": r"u_r=0,\ u_\theta=\omega_0r",
    "3.23": r"\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}",
    "3.24": r"\Gamma=2\pi ru_\theta=2\pi r^2\omega_0",
    "3.25": r"u_r=0,\ u_\theta=B/r",
    "3.26": r"\Gamma=2\pi B",
    "3.27": r"[\omega_z]_{r\to0}=\lim_{r\to0}2B/r^2",
    "3.28": r"u_\theta=\frac{\Gamma r}{2\pi\sigma^2}\ (r\le\sigma),\ \frac{\Gamma}{2\pi r}\ (r>\sigma)",
    "3.29": r"u_\theta=\frac{\Gamma}{2\pi r}\big(1-e^{-r^2/\sigma^2}\big)",
    "3.30": r"\frac{d}{dt}\int_{a}^{b}F\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)",
    "3.31": r"\frac{d}{dt}\int_{V^*}F\,dV=\lim_{\Delta t\to0}\frac{1}{\Delta t}\Big\{\int_{V^*(t+\Delta t)}F(t+\Delta t)dV-\int_{V^*(t)}F\,dV\Big\}",
    "3.32": r"\int_{V^*(t+\Delta t)}F(t+\Delta t)dV\cong\int_{V^*}F\,dV+\int_{V^*}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV+\int_{\Delta V}\Delta t\frac{\partial F}{\partial t}dV",
    "3.33": r"\frac{d}{dt}\int_{V^*}F\,dV=\lim_{\Delta t\to0}\frac{1}{\Delta t}\Big\{\int_{V^*}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV\Big\}",
    "3.34": r"\int_{\Delta V}F\,dV\cong\int_{A^*}F\,(\mathbf b\Delta t\cdot\mathbf n)\,dA",
    "3.35": r"\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA",
}


def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ *({num})*"


_EQ_REF = re.compile(r"\((3\.\d+)\)")


def show_eqs(text: str) -> str:
    """For a text that cites a book equation by number without writing it, write the equation next to its first
    mention — the house rule "show the equation, not just its number"."""
    done: set[str] = set()

    def rep(m):
        n = m.group(1)
        nxt = text[m.end():m.end() + 4]
        if (n in done or n not in EQ or EQ[n] in text or nxt.lstrip(", :").startswith("$")
                or nxt[:1] in "/–-)" or text[max(0, m.start() - 1):m.start()] == "("):
            return m.group(0)
        done.add(n)
        return f"({n}), ${EQ[n]}$"
    return _EQ_REF.sub(rep, text)


# ---------------------------------------------------------------------------------------------------------------------
# Part F reader: the derivations, word for word
# ---------------------------------------------------------------------------------------------------------------------
def _join(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(x.strip() for x in lines)).strip()


def _tex_escape(s: str) -> str:
    return (s.replace("\\", r"\backslash ").replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")
            .replace("_", r"\_").replace("{", r"\{").replace("}", r"\}"))


def display_tex(s: str, width: int = 78) -> str:
    """A Part F line (maths in `$…$`, possibly with words around it) as the body of one display-math block: a single
    formula is used as is; mixed words and maths become \\text{…} runs, wrapped into an array for long lines."""
    s = s.strip()
    parts = re.split(r"(\$[^$]+\$)", s)
    parts = [p for p in parts if p.strip()]
    if len(parts) == 1 and parts[0].startswith("$"):
        return parts[0][1:-1]
    tokens: list[tuple[str, str]] = []                       # ("m", latex) or ("w", word)
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


def part_f() -> dict[str, dict]:
    """Parse Part F of the design into {D01: dict(title, goal, start, plan, tools, assumptions, steps, result, check,
    meaning, traps, code)}."""
    text = (ROOT / "analysis" / "ch03_design.md").read_text(encoding="utf-8")
    part = text.split("## Part F", 1)[1]
    chunks = re.split(r"^### (D\d\d) · ", part, flags=re.M)
    out: dict[str, dict] = {}
    for key, body in zip(chunks[1::2], chunks[2::2]):
        lines = body.splitlines()
        title = re.split(r" — ★", lines[0])[0].strip()
        fields: dict[str, list[str]] = {}
        code: list[str] = []
        steps: list[list[str]] = []
        cur = None
        in_code = False
        for ln in lines[1:]:
            if in_code:
                if ln.strip().startswith("```"):
                    in_code = False
                else:
                    code.append(ln[2:] if ln.startswith("  ") else ln)
                continue
            if ln.strip().startswith("```"):
                in_code = True
                continue
            m = re.match(r"^- \*\*(.+?)\*\*\s*(.*)$", ln)
            if m:
                name = m.group(1).strip().rstrip(".").strip()
                cur = name
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
        parsed_steps = []
        for st in steps:
            s = _join(st)
            bits = re.split(r"(?:^|\s·\s)\*(did|tex|why|plain|live|set|watch):\*\s*", s)
            d = {}
            for name, val in zip(bits[1::2], bits[2::2]):
                d[name] = val.strip()
            parsed_steps.append(dict(did=d["did"], tex=display_tex(d["tex"]), why=d["why"], plain=d["plain"]))
        start_tex, start_plain = _split_words(f["Start"])
        res_tex, res_plain = _split_words(f["Result"])
        plan = [p.strip() for p in re.split(r"\(\d+\)\s*", f["Plan"]) if p.strip()]
        tools = [t.strip() for t in f["Tools"].split(" · ") if t.strip()]
        check = re.sub(r"\s*Optional `check_src`:.*$", "", f.get("Check", "")).strip()
        out[key] = dict(title=title, goal=f["Goal"], start=(display_tex(start_tex), start_plain), plan=plan,
                        tools=tools, assumptions=f.get("Assumptions", ""), steps=parsed_steps,
                        result=(display_tex(res_tex), res_plain), check=check, meaning=f.get("What it means", ""),
                        traps=f.get("Traps", ""), code="\n".join(code).strip("\n"))
    return out


PF = part_f()


def _pf_replace(key: str, field: str, old: str, new: str) -> None:
    """Edit one Part F field after parsing (review fixes); fails loudly if the design text changed."""
    d = PF[key]
    if field.startswith("step"):
        i = int(field[4:]) - 1
        assert old in d["steps"][i]["why"], (key, field, old)
        d["steps"][i]["why"] = d["steps"][i]["why"].replace(old, new)
        return
    if field == "tools":
        hit = [j for j, t in enumerate(d["tools"]) if old in t]
        assert hit, (key, field, old)
        d["tools"][hit[0]] = d["tools"][hit[0]].replace(old, new)
        return
    assert old in d[field], (key, field, old)
    d[field] = d[field].replace(old, new)


# D16 check: the exact semi-axes round to 1.0512 and 0.9512 (singular values 1.051249 and 0.951249)
# D18: write (3.24) and (3.26) out instead of "((3.24)/(3.26))"
_pf_replace("D18", "tools", "(3.24)/(3.26) circulation of a centred circle",
            r"the circulation of a centred circle, (3.24), $\Gamma=2\pi ru_\theta=2\pi r^2\omega_0$ (solid body) and (3.26), "
            r"$\Gamma=2\pi ru_\theta=2\pi B$ (line vortex)")
_pf_replace("D18", "step4", "The centred-circle circulation 2πr u_θ ((3.24)/(3.26)) outside the core",
            r"The centred-circle circulation, (3.24), $\Gamma=2\pi ru_\theta=2\pi r^2\omega_0$ (solid body) and (3.26), "
            r"$\Gamma=2\pi ru_\theta=2\pi B$ (line vortex): outside the core $ru_\theta=\Gamma/2\pi$, so the circle gives Γ")
# D19: (3.24)-type circle integral written out; the check now points to a real sympy cell
_pf_replace("D19", "tools", "(3.24)-type circle integral", r"the circle integral $\Gamma=2\pi ru_\theta$ of (3.24)")
_pf_replace("D19", "check", "(the notebook's sympy line)", "(the sympy cell below)")
# D20: drop the pointer to a private working file
_pf_replace("D20", "check", " A published fit of the Lamb–Oseen peak gives r²/σ² ≈ 1.256 (analysis §8) ✓.", "")
# D21/D23: explainer labels by name
_pf_replace("D21", "check", "E7's default", "the transport explainer's default")
_pf_replace("D23", "check", "E7's ellipse", "The transport explainer's ellipse")
_pf_replace("D23", "check", "`ch03.material_volume_rate` returns equal surface flux and ∫∇·u (1e-10) ✓.",
            "`ch03.material_volume_rate` returns equal surface flux and ∫∇·u (to 1e-10, run in the C15 code cell below) ✓.")
_pf_replace("D23", "tools", "Gauss' divergence theorem (2.30)",
            r"Gauss' divergence theorem (2.30), $\int_V\nabla\cdot\mathbf Q\,dV=\oint_A\mathbf Q\cdot\mathbf n\,dA$")
# D10 step 5: the "axes coincide" assumption
_pf_replace("D10", "step5", "so only R changes.",
            "so only R changes. For a rotating observer the comparison is made at the instant the two sets of axes "
            "coincide (as in D13); later the components of S turn as a tensor, but S itself is unchanged.")
# D08 step 6 → three one-move steps
_d8 = PF["D08"]["steps"]
_d8[5:6] = [
    dict(did="Differentiate the squared length",
         tex=r"\ell\,\dfrac{D\ell}{Dt}=\tfrac12\dfrac{D(\delta\mathbf x\cdot\delta\mathbf x)}{Dt}=\delta\mathbf x\cdot\dfrac{D\delta\mathbf x}{Dt}=\delta\mathbf x\cdot\mathbf G\cdot\delta\mathbf x",
         why=r"For a thread $\delta\mathbf x=\ell\mathbf n$ in any direction, $\ell^2=\delta\mathbf x\cdot\delta\mathbf x$; the product rule gives $2\,\delta\mathbf x\cdot D\delta\mathbf x/Dt$, and a material element changes at $D\delta\mathbf x/Dt=\mathbf G\cdot\delta\mathbf x$ (D07, primer P99). The book skips this.",
         plain="The length changes by the part of the ends' relative velocity that lies along the thread."),
    dict(did="Divide by ℓ²",
         tex=r"\dfrac{1}{\ell}\dfrac{D\ell}{Dt}=\mathbf n\cdot\mathbf G\cdot\mathbf n",
         why=r"$\delta\mathbf x\cdot\mathbf G\cdot\delta\mathbf x=\ell^2\,\mathbf n\cdot\mathbf G\cdot\mathbf n$ and $\ell\neq0$; dividing gives a rate per unit length that does not depend on the thread's length.",
         plain="The stretching rate along n."),
    dict(did="Drop the antisymmetric part",
         tex=r"\dfrac{1}{\ell}\dfrac{D\ell}{Dt}=\mathbf n\cdot\mathbf S\cdot\mathbf n",
         why=r"$\mathbf G=\mathbf S+\tfrac12\mathbf R$ by (3.11), and $\mathbf n\cdot\mathbf R\cdot\mathbf n=0$ for an antisymmetric R ($n_in_jR_{ij}=-n_jn_iR_{ji}$ is its own negative).",
         plain="Stretching along any direction is a quadratic form of S; rotation stretches nothing."),
]
# D24 step 4 → two one-move steps
_d24 = PF["D24"]["steps"]
_d24[3:4] = [
    dict(did="Shrink the right side",
         tex=r"\displaystyle\int_{\delta V}\Big[\dfrac{\partial F}{\partial t}+\mathbf u\cdot\nabla F+F\,\nabla\cdot\mathbf u\Big]dV\to\Big(\dfrac{\partial F}{\partial t}+\mathbf u\cdot\nabla F+F\,\nabla\cdot\mathbf u\Big)\delta V",
         why="Mean-value theorem for integrals (P85): the integral of a continuous integrand is its value at some point x* inside times δV, and x* → x as the lump shrinks round x. We want values at one point.",
         plain="For a tiny lump the budget's right side is the integrand at x times the lump's volume."),
    dict(did="Shrink the left side",
         tex=r"\dfrac{D}{Dt}(F\,\delta V)=\Big(\dfrac{\partial F}{\partial t}+\mathbf u\cdot\nabla F+F\,\nabla\cdot\mathbf u\Big)\delta V",
         why="By the same theorem ∫F dV → FδV; because the volume is material (b = u), its time derivative follows the lump, which is what D/Dt means. Equate with the right side of step 2.",
         plain="The budget of one tiny lump."),
]
# step pointers after the split (old steps 4, 6 are now 5, 7)
assert _d24[7]["did"].startswith("Equate steps 4 and 6"), _d24[7]["did"]
_d24[7]["did"] = _d24[7]["did"].replace("Equate steps 4 and 6", "Equate steps 5 and 7")
_pf_replace("D24", "check", "steps 4 and 6 both give", "steps 5 and 7 both give")
_pf_replace("D24", "check", "the left side of step 4 is 0", "the left side of step 5 is 0")
_pf_replace("D24", "assumptions", "(steps 1, 4)", "(steps 1, 4–5)")
_pf_replace("D24", "assumptions", "(step 4)", "(steps 4–5)")



# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch02.py)
# ---------------------------------------------------------------------------------------------------------------------
def core(cid: str, title: str, question: str, eqs: tuple[str, ...] = ()) -> None:
    """A CORE block heading with its id; ``eqs`` = the book equations named in the title, written out under the
    question ("show the equation, not just its number")."""
    q = textwrap.dedent(question).strip()
    if eqs:
        q += "\n\n*In one line:* " + " · ".join(E(n) for n in eqs)
    nb.core(cid, f"{title} `{cid}`", question=q)


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` is exactly the Part E concept text (the ledger check matches it)."""
    nb.primer(term, f"`{pid}` · {textwrap.dedent(text).strip()}", code=code)


def note(nid: str, title: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note with its curation id."""
    body = textwrap.dedent(text).strip()
    sep = " " if body.startswith("—") else " — "              # no doubled dash when the text itself starts with one
    nb.note(f"**{title}** `{nid}`{sep}{body}", equation=equation, ref=ref)


def D(key: str, ref: str = "", check_src: str | None = None, extra_check: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then its traps as a ⚠️ callout."""
    d = PF[key]
    S = show_eqs                                     # write out any equation Part F names by number only
    src = check_src if check_src is not None else (d["code"] or None)
    check = S(d["check"]) + (f" {extra_check}" if extra_check else "")
    goal = S(d["goal"]) + (f"\n\n**Assumptions.** {S(d['assumptions'])}" if d["assumptions"] else "")
    steps = [dict(st, why=S(st["why"]), plain=S(st["plain"])) for st in d["steps"]]
    nb.derivation(key, f"{d['title']} `{key}`", ref=ref, goal=goal, start=d["start"], plan=[S(x) for x in d["plan"]],
                  uses=[S(x) for x in d["tools"]], steps=steps, result=d["result"], interpret=S(d["meaning"]),
                  check=check, check_src=src)
    if d["traps"]:
        nb.md(f"> ⚠️ **Common confusion (traps in `{key}`):** {show_eqs(d['traps'])}")


def see_read_change(see: str, read: str, change: str) -> None:
    nb.figure_notes(see, read, change)


def whatif(text: str) -> None:
    nb.md(f"**What would change if…** {textwrap.dedent(text).strip()}")


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=r"""
Kinematics describes motion without asking what causes it. A fluid can be described by following each particle
(Lagrangian) or by watching fixed points (Eulerian); the material derivative
$\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ *(Eq. 3.5)* is the bridge. Three kinds of flow line
picture what the flow does, and they differ whenever it is unsteady. Near any point the motion splits into a
deformation (the strain-rate tensor $\mathbf S$: stretching, shearing, swelling) and a rigid spin at half the vorticity.
The chapter ends with the Reynolds transport theorem, the rule for differentiating an integral over a moving volume —
the tool Chapter 4 uses to write every conservation law.
""",
    roadmap=[
        "§3.1 steady vs unsteady, 1-/2-/3-D flows, cylindrical and spherical coordinates",
        "§3.2 Lagrangian and Eulerian descriptions (C01); the material derivative D/Dt (C02)",
        "§3.3 streamlines (C03), path lines and streak lines (C04); the acceleration is the same for every steadily moving observer (C05)",
        "§3.4 relative motion near a point (C06): stretching (C07), shearing (C08), volume change (C09), spin = ½ vorticity (C10), deformation + rotation (C11), principal strain axes (C12)",
        "§3.5 shear flow, solid-body rotation and the irrotational vortex (C13); Rankine and Gaussian vortices (C14)",
        "§3.6 Leibniz's rule and the Reynolds transport theorem (C15)",
    ],
    prerequisites=[
        "partial derivatives, the chain rule and first-order Taylor expansion (Ch. 1 primers P25, P49, P26)",
        "index notation, the velocity gradient G, the split G = S + ½R and the vorticity ω = ∇×u (Ch. 2 §2.1, §2.9–2.10)",
        "eigenvalues and principal axes of a symmetric tensor (Ch. 2 §2.11, P80)",
        "Gauss' and Stokes' theorems, circulation (Ch. 2 §2.12–2.13)",
        "scipy.integrate.solve_ivp and scipy.linalg.expm (Ch. 1 P31, Ch. 2 P79)",
    ],
)
nb.explainer_index([
    ("flow_lines_unsteady", "Three lines through one point — why do they disagree?",
     "C03 C04: streamline, path line and streak line in an unsteady flow (Ex. 3.1)"),
    ("material_derivative_probe", "Why does the station warm while the air does not?",
     "C01 C02: DF/Dt = ∂F/∂t + u·∇F with a fixed probe and a drifting float"),
    ("galilean_frames_cylinder", "Steady or not — does the acceleration care?",
     "C05: the local/advective split moves with the observer, the total does not"),
    ("fluid_element_deformation", "What does each number in S measure?",
     "C06–C09: stretching, shearing and swelling measured on a moving element"),
    ("spin_and_principal_axes", "Can a straight flow make a fluid element spin?",
     "C10–C12: spin = ½ω for any pair of lines, du = S·dx + ½ω×dx, principal axes"),
    ("vortex_paddle_wheels", "Going round in circles ≠ spinning",
     "C13 C14: solid-body, irrotational, Rankine and Gaussian vortices"),
    ("reynolds_transport_cv", "What changes inside a moving box?",
     "C15: Leibniz and the Reynolds transport theorem as a budget"),
])
nb.setup()
nb.code("""
import numpy as np                                      # arrays (Ch. 1 primer P03)
import sympy as sp                                      # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                         # static figures (Ch. 1 primer P01)
import plotly.graph_objects as go                       # rotatable 3-D and slider figures (Ch. 1 primer P41)
from fluidpy import ch03_kinematics as ch03             # the tested chapter-3 module: every function cites its § and Eq.
from fluidpy.core.interact import slider_figure         # plotly figure with a slider that works on the web page (P17)
from fluidpy.core.anim import animate                   # matplotlib animations (P16); show_animation came with the setup
from fluidpy.core.style import COLORS, savefig          # the house palette and a helper that saves PNGs to outputs/ch03
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


print(len([n for n in dir(ch03) if not n.startswith("_")]), "public names in fluidpy.ch03_kinematics")   # the toolbox this notebook calls
""", explain="""
1. Numerical, symbolic and plotting libraries (all primed in Ch. 1).
2. `ch03` is the chapter module; it re-exports the Ch. 3 primitives in `fluidpy/core/` (kinematics, coordinates,
   vortices, transport) and the Ch. 2 tensor tools, so one name covers the whole chapter. Every function is tested in
   `tests/test_ch03.py`.
3. The house helpers for slider figures, animations, colours and the observed order of convergence;
   `recolor` only restyles plotly traces so that a colour always means the same thing (below).
""")
nb.md(r"""
**Notation and colours used in this notebook.**

| Symbol | Meaning | Unit |
|---|---|---|
| $\mathbf x=(x_1,x_2,x_3)=(x,y,z)$ | position | m |
| $t$ | time | s |
| $\mathbf u=(u_1,u_2,u_3)=(u,v,w)$ | velocity | m/s |
| $F$ | any field (a temperature, a velocity component …) | its own |
| $G_{ij}=\partial u_i/\partial x_j$ | velocity gradient (row = velocity component) | 1/s |
| $\mathbf S$, $\mathbf R=\mathbf G-\mathbf G^{\rm T}$ | strain-rate tensor; the book's rotation tensor (no ½) | 1/s |
| $\boldsymbol\omega=\nabla\times\mathbf u$ | vorticity; **the spin of a fluid element is ½ω** | 1/s |
| $\gamma=du_1/dx_2$ | shear rate ($=2S_{12}$; Ch. 2's Γ was $S_{12}$ itself) | 1/s |
| $\Gamma$ | circulation | m²/s |
| $\theta$ | plane polar angle (§3.5) · polar angle from $+z$ (spherical) · cone half-angle (Ex. 3.2) | rad |

**Colours** (one meaning each, also in the explainers): streamline **teal** · path line **orange** · streak line
**rose** · local rate $\partial/\partial t$ **blue** · advective rate **amber** · total $D/Dt$ **purple** · strain part
**teal** · rotation part **orange** · stretching **blue** / compressing **rose** · volume term **blue**, surface term
**orange**.
""")
nb.md(r"""
**Where this chapter is used later**

| Result here | Used in |
|---|---|
| $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ (3.5) | every conservation law of Ch. 4, the vorticity equation (Ch. 5), potential vorticity (Ch. 13) |
| path lines $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ (3.8) | particle orbits under waves (Ch. 7), Lagrangian statistics (Ch. 12), trajectories (Ch. 13) |
| $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ (3.19) | the Newtonian stress law (Ch. 4 §4.5), vortex stretching (Ch. 5) |
| circulation $\Gamma=\oint\mathbf u\cdot d\mathbf s$ (3.18) | Kelvin's theorem (Ch. 5), lift (Ch. 6, 14) |
| $\omega'_z=\omega_z-2\Omega$ | absolute vs relative vorticity $f+\zeta$ (Ch. 4 §4.7, Ch. 13) |
| $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$ (3.35) | mass, momentum and energy equations (Ch. 4 §4.2–4.4), layer budgets (Ch. 13) |

*Climate hook:* the advective term $\mathbf u\cdot\nabla T$ is the "warm/cold advection" shaded on every weather map.
""")
nb.md(r"""
🔁 **Tools from earlier chapters used here without a new primer** (see `knowledge/primers.md`): f-strings `f"{x:.3f}"`
(Ch. 1 P04) · `lambda` and functions passed as arguments (P29) · `assert np.allclose(a, b)` — "the two agree to
rounding" (P15) · tuple unpacking `u, a = f(...)` (P14) · `show_viz` embeds an explainer (P18) · `animate` /
`show_animation` (P16) and `slider_figure` (P17) · the ordinary derivative as a slope (P19) and the definite integral
(P27) · `np.trapezoid` (P37) · matrix–vector products `G @ dx` (Ch. 2 P63) · `np.linalg.norm` (P67) · the quadratic form
$\mathbf n\cdot\mathbf A\cdot\mathbf n$ (P82; it vanishes for an antisymmetric A because $n_in_jA_{ij}=-n_jn_iA_{ji}$) ·
eigenvalues and eigenvectors (P80) · `np.arctan2` (P70) · midpoint sums for volume and surface integrals (P83).

**Small idioms, glossed once:** `zip(a, b)` walks two lists in step · `enumerate(a)` gives (index, item) pairs · `np.hypot(x, y)` $=\sqrt{x^2+y^2}$ · `np.ptp(a)` = max − min ("peak to peak") · `np.r_[a, b]` joins arrays end to end, `np.full(n, v)` is n copies of v · `np.ma.masked_invalid(a)` hides NaN entries from a plot · a result such as `ch03.acceleration(...)` is a *named tuple*: unpack it (`a, loc, adv = …`) or read a field (`res.local`) · `plt.contour(...).allsegs` returns the contour lines as lists of points · `sp.Function('f')` is an unknown function, `sp.Lambda(args, expr)` a function built from an expression · `sp.series(e, x, 0, n)` a Taylor series, `sp.limit(e, x, a)` a limit.
""")

# =====================================================================================================================
# A.1 §3.1 Introduction and Coordinate Systems — no A item
# =====================================================================================================================
nb.section("3.1", "Introduction and Coordinate Systems", intro="""
**What is this section about?** The words used to describe *any* flow before we ask what drives it: steady or
unsteady, one-, two- or three-dimensional, and the coordinate systems whose velocity components later chapters use
(Cartesian, plane polar, cylindrical, spherical).
""")
note("N01", "Kinematics", r"""
describes how a fluid moves — positions, velocities, accelerations, stretching and spinning — without asking which
forces cause it. Forces enter in Ch. 4 (§4.4 onward, Newton's second law for a fluid).
""")
note("N02", "Steady and unsteady", r"""
A flow is **steady** when nothing measured at a fixed point changes with time: $\partial(\cdot)/\partial t=0$ for every
field. Otherwise it is **unsteady**. Careful: steadiness depends on the observer — the air flowing past a parked car
seen from the car is steady; the same kind of flow seen from the pavement as a car drives by is not. C05 turns this into
a theorem.
""", equation=r"\frac{\partial(\cdot)}{\partial t}=0\quad\text{(steady)}")
P("P87", "scipy.integrate.quad and dblquad", r"""
Adaptive numerical integration: `quad(f, a, b)` returns $\int_a^b f\,dx$ and an error estimate; `dblquad(f, a, b,
gfun, hfun)` does a double integral, the **inner variable first** in `f`'s arguments (its limits may depend on the outer
one). We use them for section averages here, for the cone of Ex. 3.2 and for exact checks of the Reynolds transport
theorem in §3.6.
""", code="""
from scipy.integrate import quad, dblquad                     # adaptive 1-D and 2-D quadrature
val, err = quad(lambda r: 2*r, 0.0, 1.0)                      # ∫₀¹ 2r dr = 1
print(val, err)                                               # 1.0 and a tiny error estimate
area, _ = dblquad(lambda r, th: r, 0.0, 2*np.pi, 0.0, 1.0)    # ∫∫ r dr dθ over the unit disc (inner variable r first)
print(area, np.pi)                                            # 3.14159… twice: the disc's area
""")
note("N03", "1-D, 2-D and 3-D flows", r"""
A flow is 3-D when it depends on all three coordinates, plane (2-D) when one coordinate can be dropped, and — as an
engineering approximation — 1-D when we keep only the average over each cross-section of a pipe or channel (an
axisymmetric pipe flow is still 3-D in the book's sense). The 1-D description keeps the average $\bar u$ over the
section area $A$ and throws away the profile. Number: the parabolic (Poiseuille) profile $u=U(1-r^2/R^2)$ averages to
exactly $U/2$.
""", equation=r"\bar u(z)=\frac{1}{A}\int_A u\,dA")
nb.code("""
U, R = 2.0, 0.05                                     # centre speed 2 m/s, pipe radius 5 cm
u_pois = lambda r, z: U*(1 - (r/R)**2)               # parabolic profile u(r, z) [m/s] (the same at every z here)
ubar = ch03.cross_section_average(u_pois, R)         # (1/πR²) ∫₀ᴿ u 2πr dr by quad: rings of area 2πr dr
print(f"section average = {ubar:.4f} m/s (U/2 = {U/2})")
r = np.linspace(0, R, 2001)                          # 2001 radii from the axis to the wall [m]
mine = np.trapezoid(u_pois(r, 0)*2*np.pi*r, r)/(np.pi*R**2)   # the same average by the trapezoid rule (Ch. 1 P37)
assert np.allclose(mine, ubar, rtol=1e-6)            # same number → the library does exactly this
""", explain="""
1. `u_pois` is the profile as a function of radius r and axial position z.
2. `cross_section_average` integrates $u\\cdot2\\pi r\\,dr$ (a thin ring of radius r has area $2\\pi r\\,dr$) with `quad` and
   divides by the area $\\pi R^2$.
3. The trapezoid version on 2001 rings agrees — the from-scratch check of this section.
""")
nb.md(r"""
Our model of a developing pipe flow (not a solution of the equations of motion — Ch. 8 has the real one):
$u=U_m\big(1-(r/R)^{n}\big)$ with $n(z)=2+8e^{-z/L_e}$ ($L_e=0.5$ m: flat at the inlet, a parabola far downstream) and
$U_m=\bar u\,(n+2)/n$, which keeps the section average $\bar u$ the same at every z.
""")
nb.figure(r"""
R, Ub = 1.0, 1.0                                     # pipe radius 1 m, section-average speed 1 m/s
r = np.linspace(-R, R, 201)                          # across the pipe, wall to wall [m]
zs = [0.0, 0.5, 2.0]                                 # three stations downstream of the inlet [m]
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8))  # left: profiles (2-D view); right: averages (1-D view)
for z, al in zip(zs, (0.45, 0.7, 1.0)):              # darker teal = further downstream
    a.plot(ch03.pipe_profile(r, z, Ub, R, n0=10.0, L_e=0.5), r, color=COLORS["teal"], alpha=al, lw=2.2,
           label=f"$z$ = {z} m")                     # our developing profile: flat at the inlet → parabola downstream
a.set_xlabel("$u$ [m/s]"); a.set_ylabel("$r$ [m]"); a.legend(fontsize=8)   # axes with units
a.set_title("2-D view: the profile changes along the pipe")
z = np.linspace(0, 3, 61)                            # stations along the pipe [m]
ubar = [ch03.cross_section_average(lambda rr, zz: ch03.pipe_profile(rr, zz, Ub, R, 10.0, 0.5), R, zz) for zz in z]   # ū(z) by quad
n = 2 + (10 - 2)*np.exp(-z/0.5)                      # the profile's exponent n(z) = 2 + 8 e^{−z/L_e}
b.plot(z, ubar, color=COLORS["accent"], lw=2.5, label=r"$\bar u(z)$ — the whole 1-D description")   # constant: 1 m/s
b.plot(z, Ub*(n + 2)/n, "--", color=COLORS["muted"], label="centre speed $U_m(z)$ (lost in 1-D)")  # what 1-D forgets
b.set_ylim(0, 2.2); b.set_xlabel("$z$ [m]"); b.set_ylabel("speed [m/s]"); b.legend(fontsize=8)     # fixed axes
b.set_title("1-D view: one number per section")
savefig(fig, "ch03", "pipe_profiles"); plt.show()    # save to outputs/ch03 and draw
""", see="Three teal profiles that sharpen downstream (left) and one flat purple line with a dashed curve above it (right).",
    read="""The 1-D model sees only the right panel: the same average at every z while the shape changes completely.
Mass is conserved (the average does not change); the detail — where the fluid is fast — is thrown away. The dashed
centre speed rises from ≈ 1.2 to 2 m/s, a change the 1-D description cannot see.""",
    change="""…the pipe narrowed to half the radius: ū would jump ×4 (the same volume per second through a quarter of the
area); what the profile looks like would still be an open question for the 1-D model.""")
nb.md("""
Fig. 3.2's two views of one cylinder — held fixed in a stream, or towed through still water — are the same flow seen by
two observers; we take them up with the Galilean transformation in C05 (§3.3), where note `N04` treats them.
""")
P("P88", "cylindrical and spherical unit vectors", r"""
At every point P the curvilinear coordinates carry their own right-handed unit vectors, and they turn as P moves.
Cylindrical $(R,\varphi,z)$: $\mathbf e_R=(\cos\varphi,\sin\varphi,0)$, $\mathbf e_\varphi=(-\sin\varphi,\cos\varphi,0)$,
$\mathbf e_z=(0,0,1)$. Spherical $(r,\theta,\varphi)$ with θ measured from $+z$:
$\mathbf e_r=(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta)$,
$\mathbf e_\theta=(\cos\theta\cos\varphi,\cos\theta\sin\varphi,-\sin\theta)$, $\mathbf e_\varphi=(-\sin\varphi,\cos\varphi,0)$.
A velocity component is the dot product of $\mathbf u$ with one of them — a projection on an orthonormal basis (Ch. 2
primer P65).
""", code="""
phi = np.pi/4                                                    # azimuth φ = 45°
e_R = np.array([np.cos(phi), np.sin(phi), 0.0])                  # points away from the z-axis
e_phi = np.array([-np.sin(phi), np.cos(phi), 0.0])               # points round the z-axis
print(e_R @ e_phi, np.cross(e_R, e_phi))                         # 0.0 and (0, 0, 1): perpendicular, e_R × e_φ = e_z
""")
note("N05", "Coordinate systems (Fig. 3.3)", r"""
Plane: $(x,y)=(x_1,x_2)$ or polar $(r,\theta)$ (Ch. 2 Ex. 2.1). Cylindrical $(R,\varphi,z)$ with velocity components
$(u_R,u_\varphi,u_z)$; spherical $(r,\theta,\varphi)$, θ the angle from $+z$, φ the azimuth, components
$(u_r,u_\theta,u_\varphi)$. ⚠️ The same letter θ means the plane polar angle in §3.5 but the angle from the $z$-axis in
spherical coordinates. The conversions (the book states them in its exercises):
""", equation=r"R=\sqrt{x^2+y^2},\quad \varphi=\tan^{-1}(y/x);\qquad r=\sqrt{x^2+y^2+z^2},\quad \theta=\tan^{-1}\!\big(\sqrt{x^2+y^2}/z\big)")
nb.code("""
P = (1.0, 1.0, 1.0)                                              # a point P [m]
print(np.round(ch03.cylindrical_from_cartesian(*P), 4))          # (R, φ, z) = (1.4142, 0.7854, 1.0)
print(np.round(ch03.spherical_from_cartesian(*P), 4))            # (r, θ, φ) = (1.7321, 0.9553, 0.7854); θ = 54.74°
u = np.array([1.0, 0.0, 0.0])                                    # a 1 m/s wind along x at P
uc = ch03.velocity_components(u, np.array(P), "cylindrical")     # (u_R, u_φ, u_z): projections on e_R, e_φ, e_z
us = ch03.velocity_components(u, np.array(P), "spherical")       # (u_r, u_θ, u_φ): projections on e_r, e_θ, e_φ
print(np.round(uc, 4), np.round(us, 4))                          # (0.7071, −0.7071, 0) and (0.5774, 0.4082, −0.7071)
print(np.sum(uc**2), np.sum(us**2))                              # 1.0 twice: a projection on an orthonormal basis keeps |u|²
print(np.round(ch03.unit_vectors_spherical(*ch03.spherical_from_cartesian(*P)[1:]), 4))   # rows e_r, e_θ, e_φ at P
""", explain="""
1. Position P in the two curvilinear systems (angles in radians).
2. The same wind vector projected on the local unit vectors of each system.
3. The squares add to $|\\mathbf u|^2=1$ in every system — the three unit vectors are perpendicular and of length one.
""")
nb.plotly(r"""
P = np.array([1.0, 1.0, 1.0])                                    # the point P [m]
_, th, ph = ch03.spherical_from_cartesian(*P)                    # its spherical angles
triads = {"Cartesian": (np.eye(3), COLORS["muted"]),             # rows e_x, e_y, e_z — the same everywhere
          "cylindrical": (ch03.unit_vectors_cylindrical(ph), COLORS["teal"]),   # rows e_R, e_φ, e_z at P
          "spherical": (ch03.unit_vectors_spherical(th, ph), COLORS["orange"])} # rows e_r, e_θ, e_φ at P
L = 0.6                                                          # drawn arrow length [m]
fig = go.Figure()                                                # an empty 3-D plotly figure
groups = []                                                      # which traces belong to which triad (for the menu)
for name, (E3, col) in triads.items():                           # one triad at a time
    xs, ys, zs = [], [], []                                      # shaft coordinates, NaN-separated
    for e in E3:                                                 # one shaft per unit vector, all in one trace
        xs += [P[0], P[0] + L*e[0], None]; ys += [P[1], P[1] + L*e[1], None]; zs += [P[2], P[2] + L*e[2], None]
    fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=col, width=6), name=name))   # shafts
    tips = P + L*E3                                              # arrow heads at the tips (plotly cones: Ch. 2 P64)
    fig.add_trace(go.Cone(x=tips[:, 0], y=tips[:, 1], z=tips[:, 2], u=E3[:, 0], v=E3[:, 1], w=E3[:, 2],
                          anchor="tip", sizemode="absolute", sizeref=0.15, showscale=False,
                          colorscale=[[0, col], [1, col]], name=name))   # one-colour cones
    groups += [name, name]                                       # shafts and heads belong to this triad
fig.add_trace(go.Scatter3d(x=[0, P[0], P[0], P[0]], y=[0, 0, P[1], P[1]], z=[0, 0, 0, P[2]], mode="lines",
                           line=dict(color=COLORS["grid"], dash="dash"), name="guide to the axes"))   # origin → P
fig.add_trace(go.Scatter3d(x=[P[0]], y=[P[1]], z=[P[2]], mode="markers", marker=dict(size=5, color=COLORS["ink"]), name="P"))
groups += ["always", "always"]                                   # guide and P are always shown
buttons = [dict(label=lab, method="update",
                args=[{"visible": [g == "always" or lab == "all" or g == lab for g in groups]}])
           for lab in ("all", "Cartesian", "cylindrical", "spherical")]   # a dropdown to show one triad at a time
fig.update_layout(height=460, title="Unit vectors at P = (1, 1, 1) m: Cartesian (grey), cylindrical (teal), spherical (orange)",
                  updatemenus=[dict(buttons=buttons, x=0, y=1.08, xanchor="left")],
                  scene=dict(aspectmode="cube", xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]"))   # equal axes
fig.show()                                                       # draw it (rotatable)
""", explain="""
1. `unit_vectors_cylindrical` and `unit_vectors_spherical` return the local triads at P as rows.
2. Each triad is drawn as three shafts (one trace, NaN-separated) plus three cone heads; the dashed grey line walks from
   the origin to P along the axes.
3. The dropdown switches the triads on and off.
""")
see_read_change(
    "Three triads of arrows at one point, one grey, one teal, one orange.",
    """The Cartesian triad is the same everywhere; the teal and orange triads are tied to P: $\\mathbf e_R$ points away
from the $z$-axis, $\\mathbf e_r$ away from the origin, and $\\mathbf e_\\varphi$ is shared by both.""",
    """…P moved to (−1, 1, 1) (edit `P` in the cell): $\\mathbf e_R$ and $\\mathbf e_\\varphi$ swing by 90° while the grey
arrows stay put — this is why derivatives of curvilinear components need extra terms (note `N06`).""")
note("N06", "Curvilinear operators", r"""
Appendix B of the book lists ∇, ∇², ∇·u and (u·∇)u in cylindrical and spherical coordinates; we meet the first of
them in C13 (the vorticity in polar coordinates,
$\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ *(Eq. 3.23)*, derived
there by hand) and the rest with Ch. 4's Navier–Stokes equations.
""")

# =====================================================================================================================
# A.2 §3.2 Particle and Field Descriptions — C01, C02
# =====================================================================================================================
nb.section("3.2", "Particle and Field Descriptions of Fluid Motion", intro="""
**What is this section about?** Two ways to describe one flow — ride along with every particle, or stand still and
watch every point — and the single formula, the material derivative, that turns one into the other.
""")
core("C01", "Following a particle or watching a point — and the bridge between them (3.2)", r"""
A float drifts down a river past a thermometer bolted to a bridge pier. Both report the temperature of the same water.
How do the two sets of numbers fit together?
""", eqs=("3.2",))
nb.md(r"""
#### The problem in plain words

Oceanographers throw Argo floats into the sea and read where each one goes (it *is* a water parcel); weather services
keep thermometers at fixed stations; numerical models store fields on a fixed grid. We need to translate freely between
"what happens to this parcel" and "what happens at this place" — fluid dynamics is written in the second language, but
Newton's law is about the first.
""")
nb.md(r"""
#### The idea

```
Lagrangian (follow the particle)            Eulerian (watch the point)
label r_o at time t_o  →  r(t; r_o, t_o)    field F(x, t): four independent variables x, y, z, t
"what does THIS parcel do?"                 "what happens HERE?"
         the bridge:  F[r(t; r_o, t_o), t] = F(x, t)   when   x = r(t; r_o, t_o)      (3.2)
```

In symbols: $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ when $\mathbf x=\mathbf r(t;\mathbf r_o,t_o)$ *(3.2)*.

**A label is not a variable**: $\mathbf r_o$ only names the particle (where it was at $t_o$), just as a float's serial
number does. The bridge says: at time t, the field at $\mathbf x$ *is* whatever the particle that happens to be at
$\mathbf x$ carries.
""")
P("P89", "functions of time with parameters", r"""
A trajectory $x=Xe^{\alpha t}$ is a function of time t alone once the label X is chosen; X is a *parameter* that picks
one curve out of a family. Differentiating "with the label held fixed" means following one member of the family. Here
X is where the particle was at t = 0.
""", code="""
alpha = 0.5                                          # stretching rate α [1/s]
x_of = lambda t, X: X*np.exp(alpha*t)                # one curve per label X [m]
print([round(float(x_of(1.0, X)), 3) for X in (1.0, 2.0)])   # [1.649, 3.297]: two particles at t = 1 s
""")
note("N07", "Lagrangian description", r"""
Each particle is labelled by its position $\mathbf r_o$ at a reference time $t_o$; its later position is
$\mathbf r(t;\mathbf r_o,t_o)$ (Fig. 3.4). Any property it carries is written $F[\mathbf r(t;\mathbf r_o,t_o),t]$. The
labels are fixed numbers for a given particle, not coordinates.
""", equation=r"\mathbf r=\mathbf r(t;\mathbf r_o,t_o)")
note("N08", "Velocity and acceleration of a particle", r"""
are ordinary time derivatives along its own path — exactly single-particle mechanics. Number (the running example):
$x=Xe^{\alpha t}$ with X = 2 m, α = 0.5 s⁻¹, t = 1 s gives $u=\alpha Xe^{\alpha t}=1.649$ m/s and
$a=\alpha^2Xe^{\alpha t}=0.824$ m/s².
""", equation=r"\mathbf u=d\mathbf r(t;\mathbf r_o,t_o)/dt\quad\text{and}\quad\mathbf a=d^2\mathbf r(t;\mathbf r_o,t_o)/dt^2", ref="3.1")
nb.md(r"""
#### The maths — the Eulerian description and the bridge

Watching fixed points, a property is a field $F(\mathbf x,t)$ of four independent variables. The two descriptions must
agree when the particle position and the field point coincide, in the same coordinates and on a common clock:

$$F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)\quad\text{when}\quad\mathbf x=\mathbf r(t;\mathbf r_o,t_o).\qquad(3.2)$$

Read it as a recipe: to get the Eulerian field at $(\mathbf x,t)$, find the label of the particle that sits at
$\mathbf x$ at time t, then read that particle's value. The derivation below does this for a concrete flow.
""")
P("P90", "inverse functions and sympy solve", r"""
To use $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ at $\mathbf x=\mathbf r$ *(Eq. 3.2)* we must answer "which
particle is at x at time t?" — that is, solve $x=r(t;X)$ for the label X. `sympy.solve(equation, unknown)` does it
symbolically; it returns a list of solutions.
""", code="""
X, x, t, alpha = sp.symbols('X x t alpha', positive=True)   # label, position, time, rate (all > 0)
sp.solve(sp.Eq(x, X*sp.exp(alpha*t)), X)                    # [x*exp(-alpha*t)]: the label of the particle now at x
""")
D("D01")
nb.worked_example("the stretching flow x = X e^{αt} with X = 2 m, α = 0.5 s⁻¹, t = 1 s", r"""
1. Position: $x=2e^{0.5}=2\times1.6487=3.297$ m.
2. Lagrangian velocity from (3.1), $u=dr/dt$ with the label fixed: $u=\alpha Xe^{\alpha t}=0.5\times3.297=1.649$ m/s.
3. Eulerian velocity field from D01: $u(x,t)=\alpha x=0.5\times3.297=1.649$ m/s — the same number, as
   $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ at $\mathbf x=\mathbf r$ *(Eq. 3.2)* demands.
4. Acceleration: $a=\alpha^2x=0.25\times3.297=0.824$ m/s².
5. A fixed probe at x = 3.297 m reads 1.649 m/s for ever (the field is steady) while every particle passing it speeds up.
""")
nb.code("""
X, alpha, t = 2.0, 0.5, 1.0                                        # label [m], stretching rate [1/s], time [s]
x = ch03.lagrangian_map_example(X, t, alpha)                       # position x = X e^{αt} [m]
u, a = ch03.lagrangian_velocity_acceleration(lambda s: ch03.lagrangian_map_example(X, s, alpha), t)   # (3.1) by central differences
print(f"x = {x:.4f} m, u = {u:.4f} m/s, alpha*x = {alpha*x:.4f}, a = {a:.4f} m/s^2")
Xs, xs, ts, al = sp.symbols('X x t alpha', positive=True)          # the same map symbolically
res = ch03.lagrangian_to_eulerian([Xs*sp.exp(al*ts)], [Xs], [xs], ts)   # solve for the label, substitute (3.2)
print(res["label_of_x"], res["u"], res["a_lagrangian"], sp.simplify(res["Du_Dt"][0] - res["a_lagrangian"][0]))
""", explain="""
1. The map and its derivatives at fixed label, $\\mathbf u=d\\mathbf r/dt$ and $\\mathbf a=d^2\\mathbf r/dt^2$ *(Eq. 3.1)*,
   by central differences (Ch. 1 P21).
2. The symbolic route solves for the label, substitutes it into $dr/dt$ and $d^2r/dt^2$ (D01), and
3. checks that the material derivative of the Eulerian $u$ — the formula C02 derives — gives the same acceleration (the
   last printed 0).
""")
nb.check_agree("""
h = 1e-4                                                           # time step for the differences [s]
x_ = lambda s: ch03.lagrangian_map_example(X, s, alpha)            # the path of particle X
u_fd = (x_(t + h) - x_(t - h))/(2*h)                               # central difference: slope of the particle's path
assert np.allclose(u_fd, alpha*x, rtol=1e-7)                       # Lagrangian slope = Eulerian field αx at its position
a_fd = (x_(t + h) - 2*x_(t) + x_(t - h))/h**2                      # second difference: the path's curvature in time
assert np.allclose(a_fd, alpha**2*x, rtol=1e-5)                    # = α²x, the acceleration field of D01
print("from scratch:", round(u_fd, 6), round(a_fd, 5))
""")
nb.figure(r"""
alpha = 0.5                                                        # stretching rate [1/s]
t = np.linspace(0, 3, 200)                                         # time [s]
fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4))              # (a) Lagrangian view, (b) Eulerian view
for X in np.linspace(0.25, 1.5, 6):                                # six labelled particles (their x at t = 0) [m]
    a.plot(t, X*np.exp(alpha*t), color=COLORS["orange"], lw=1.8)   # path x(t) = X e^{αt}
    a.annotate(f"X = {X:.2f}", (0, X), xytext=(3, 2), textcoords="offset points", fontsize=7)   # its label
a.axhline(2.0, color=COLORS["muted"], lw=1.2)                      # a fixed probe at x = 2 m
X1 = 1.25; tc = np.log(2.0/X1)/alpha                               # when particle X = 1.25 m reaches the probe [s]
a.plot([tc - 0.4, tc + 0.4], [2 - 0.4, 2 + 0.4], color=COLORS["accent"], lw=2.5)   # tangent of slope u = αx = 1 m/s
a.annotate("slope u = αx = 1 m/s at the probe", (tc, 2.0), xytext=(tc - 1.6, 2.9), fontsize=8, color=COLORS["accent"],
           arrowprops=dict(arrowstyle="->", color=COLORS["accent"]))   # point at the tangent
a.set_xlabel("$t$ [s]"); a.set_ylabel("$x$ [m]"); a.set_ylim(0, 4)  # axes with units
a.set_title("(a) Lagrangian: one path per particle")
xg = np.linspace(0, 4, 50)                                         # positions along the line [m]
b.plot(xg, alpha*xg, color=COLORS["teal"], lw=2.5, label="$u(x)=\\alpha x$ at $t$ = 0 s")          # the field now
b.plot(xg, alpha*xg, "--", color=COLORS["ink"], lw=1.2, label="$u(x)$ at $t$ = 2 s (the same line)")   # … and later
X0 = np.linspace(0.25, 1.5, 6)                                     # the same six particles
b.plot(X0, alpha*X0, "o", color=COLORS["orange"], label="particles at $t$ = 0")                    # where they are at 0 s
b.plot(X0*np.e, alpha*X0*np.e, "s", color=COLORS["orange"], mfc="none", label="particles at $t$ = 2 s")   # x = X e^{1}
b.set_xlabel("$x$ [m]"); b.set_ylabel("$u$ [m/s]"); b.legend(fontsize=8)   # axes with units
b.set_title("(b) Eulerian: one field, the same at every time")
savefig(fig, "ch03", "lagrangian_eulerian"); plt.show()            # save to outputs/ch03 and draw
""", see="Orange curves fanning out in (a); one teal straight line in (b) that does not move, with the particles sliding along it.",
    read="""Every particle crossing the grey probe line crosses it with the same slope 1 m/s (the field is steady), but each
curve bends upward (each particle accelerates). In (b) the particles move up the line — from the circles to the
squares — so each one's velocity grows although the line itself never changes. The two panels are the two descriptions
of one motion.""",
    change="…α doubled: the fan in (a) opens twice as fast and the line in (b) doubles its slope — and still does not move in time.")
nb.md(r"""
> ⚠️ **Common confusion:** "a steady velocity field means the fluid does not accelerate." Here $u=\alpha x$ does not
> depend on t at any fixed point, yet every particle speeds up at $a=\alpha^2x$. Steadiness is a statement about
> points; acceleration is about particles. C02 shows the missing piece is the advective term.
""")
whatif(r"""
…the map were $x=X+Ut$ (uniform drift)? Then $u=U$ everywhere, a = 0, and both descriptions are trivially the same.
The interesting case is when the particle moves *through* a field that varies in space — then a thermometer on the
float and one on the pier disagree. How fast does the float's reading change? That is C02.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C02", "The material derivative: the rate a moving parcel feels (3.4)–(3.5)", r"""
A weather station reports warming of 0.36 K per hour, yet the air blowing past it has kept exactly its own
temperature. How can both be true?
""", eqs=("3.4", "3.5"))
nb.md(r"""
#### The problem in plain words

A southerly wind carries warm air north over a station. The station's thermometer rises; a thermometer tied to a
balloon drifting with the same air does not change at all. Forecasters call the station's warming *warm advection*.
We want the formula that separates "the field changes here" from "the parcel moves to where the field is different" —
the material derivative.
""")
nb.md(r"""
#### The idea

```
thermometer on the pier      reads  ∂T/∂t          (x fixed)
thermometer on the balloon   reads  DT/Dt          (label fixed)
the difference                      u·∇T           (the balloon moves through a gradient)
           DT/Dt  =  ∂T/∂t  +  u·∇T         → blue + amber = purple in every picture below
```
""")
P("P91", "multivariable chain rule along a path", r"""
If $f(t)=F(x(t),y(t),z(t),t)$, then
$\frac{df}{dt}=\frac{\partial F}{\partial x}\frac{dx}{dt}+\frac{\partial F}{\partial y}\frac{dy}{dt}+\frac{\partial F}{\partial z}\frac{dz}{dt}+\frac{\partial F}{\partial t}$:
each argument that changes contributes (sensitivity to it) × (its rate). Time appears twice — through the moving
position and through the explicit clock. This extends the two-variable chain rule of Ch. 1 (P49).
""", code="""
t = sp.symbols('t')                                  # time
x = sp.cos(t); F = lambda x, t: x**2*t               # a path x(t) and a field F(x, t)
direct = sp.diff(F(x, t), t)                         # differentiate after substituting the path
X, T = sp.symbols('X T')                             # F's own arguments
chain = (sp.diff(F(X, T), X)*sp.diff(x, t) + sp.diff(F(X, T), T)).subs({X: x, T: t})   # (∂F/∂x)(dx/dt) + ∂F/∂t
print(sp.simplify(direct - chain))                   # 0: the chain rule holds
""")
D("D02", ref="3.5", check_src="""
t, T0, G, c, v = sp.symbols('t T_0 G c v')           # time, reference temperature, gradient, front speed, wind speed
x, y = sp.symbols('x y')                             # the field's own coordinates
T = T0 - G*(y - c*t)                                 # a linear front moving north at c: T(x, y, t)
y_path = v*t                                         # the balloon's path: carried north by the wind from y = 0
direct = sp.diff(T.subs(y, y_path), t)               # step 1: the balloon's reading f(t), differentiated directly
DTDt = (sp.diff(T, t) + 0*sp.diff(T, x) + v*sp.diff(T, y)).subs(y, y_path)   # step 8: ∂T/∂t + u·∇T with u = (0, v)
print(sp.simplify(direct - DTDt))                    # 0: (3.5) gives the balloon's own rate
print(sp.simplify(direct))                           # G(c − v): zero when the pattern moves with the wind
""", extra_check="The sympy cell below applies the chain rule along the path $y=vt$ to the front $T=T_0-G(y-ct)$.")
note("N09", "The chain rule written out, (3.3)", r"""
Step 3 of the derivation is the book's (3.3) written out — the chain rule along a particle's path, with the four
arguments of F each contributing:
""", equation=r"\frac{d}{dt}F[\mathbf r(t;\mathbf r_o,t_o),t]=\frac{\partial F}{\partial r_1}\frac{dr_1}{dt}+\frac{\partial F}{\partial r_2}\frac{dr_2}{dt}+\frac{\partial F}{\partial r_3}\frac{dr_3}{dt}+\frac{\partial F}{\partial t}=\frac{d}{dt}F(\mathbf x,t)\quad\text{when}\quad\mathbf x=\mathbf r(t;\mathbf r_o,t_o)", ref="3.3")
note("N09", "…and (3.4)", r"""
Steps 4–7 turn it into the book's (3.4) — velocity components in place of $dr_i/dt$, the field's slopes at
$\mathbf x=\mathbf r$, and a name for the result: the **material**, **substantial** or **particle** derivative.
""", equation=r"\frac{d}{dt}F[\mathbf r(t;\mathbf r_o,t_o),t]=\frac{\partial F}{\partial x_1}u_1+\frac{\partial F}{\partial x_2}u_2+\frac{\partial F}{\partial x_3}u_3+\frac{\partial F}{\partial t}=(\nabla F)\cdot\mathbf u+\frac{\partial F}{\partial t}\equiv\frac{D}{Dt}F(\mathbf x,t)", ref="3.4")
note("N11", "Vector and index forms", r"""
the forms used for the rest of the book. The repeated index i is the dot product $\mathbf u\cdot\nabla F$ (summation
convention, Ch. 2 §2.1); the code below expands it with the Ch. 2 index parser (comma notation: `F,i` means
$\partial F/\partial x_i$).
""", equation=r"\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F,\quad\text{or}\quad\frac{DF}{Dt}\equiv\frac{\partial F}{\partial t}+u_i\frac{\partial F}{\partial x_i}", ref="3.5")
nb.code("""
print(ch03.expand_indices_str("u_i F,i"))            # u_i ∂F/∂x_i with the repeated i summed over 1, 2, 3
""", explain="The repeated index i becomes the three-term sum $u_1\\,\\partial F/\\partial x_1+u_2\\,\\partial F/\\partial x_2+u_3\\,\\partial F/\\partial x_3=\\mathbf u\\cdot\\nabla F$ (Ch. 2 §2.1).")
note("N10", "Local and advective parts", r"""
$\partial F/\partial t$ (**blue**) is the *local* or unsteady rate — what a fixed probe sees; it vanishes when F does
not depend on time. $\mathbf u\cdot\nabla F$ (**amber**) is the *advective* rate — the change a particle meets because
it moves to where F is different; it vanishes when F is uniform, when u = 0, or when $\mathbf u\perp\nabla F$ (wind
blowing along the isotherms). The book says *advection* for transport by the flow and keeps *convection* for heat
carried by fluid motion.
""")
note("N12", "The streamwise form", r"""
Because $\mathbf u\cdot\nabla F=|\mathbf u|\,(\mathbf e_u\cdot\nabla F)$ with the unit vector $\mathbf e_u=\mathbf u/|\mathbf u|$,
and $\mathbf e_u\cdot\nabla F=\partial F/\partial s$ is the directional derivative along the path (Ch. 2 primer P75),
the material derivative can be written with the arc length s along the particle's path. ⚠️ **Book typo:** the printed
(3.6) ends with $|\mathbf u|\,\partial/\partial s$ — the F is missing (both terms must be rates of change of F). The
form is undefined at a stagnation point (u = 0, no direction).
""", equation=r"\frac{DF}{Dt}=\frac{\partial F}{\partial t}+|\mathbf u|\frac{\partial F}{\partial s}", ref="3.6")
nb.md(r"""
> ⚠️ **Common confusion:** "$\partial T/\partial t$ is how fast the air warms." It is how fast a *fixed thermometer*
> warms. The air's own rate is $DT/Dt$. They differ by the advective term $\mathbf u\cdot\nabla T$ — often the largest
> term on a weather map.
""")
nb.worked_example("warm advection over a station", r"""
Temperature falls northward by 1 K per 100 km: $\partial T/\partial y=-1/10^5=-1\times10^{-5}$ K/m. A southerly wind
v = 10 m/s blows north.

1. Advective term: $\mathbf u\cdot\nabla T=v\,\partial T/\partial y=10\times(-10^{-5})=-1\times10^{-4}$ K/s.
2. The air keeps its temperature (no heating): $DT/Dt=0$.
3. From $\frac{DT}{Dt}=\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T$ *(Eq. 3.5)*:
   $\partial T/\partial t=0-(-10^{-4})=+1\times10^{-4}$ K/s.
4. Per hour: $10^{-4}\times3600=0.36$ K/h of warming at the station, although no parcel warmed at all.
""")
nb.code("""
G, v = 1e-5, 10.0                                          # gradient: 1 K per 100 km [K/m]; southerly wind [m/s]
terms = ch03.thermal_front_terms(0.0, 0.0, 0.0, 0.0, v, G, heating_K_per_s=0.0, front_speed=v)   # exact terms of (3.5)
print({k: f"{float(terms[k]):+.2e}" for k in ("local", "advective", "total")})   # K/s: blue, amber, purple
print(f"station warming: {float(terms['local'])*3600:.2f} K/h ({terms['regime']})")
F = lambda x, t: ch03.thermal_front(x[0], x[1], t, G, 0.0, v)   # the same field T(x, t) as a callable [K]
u = lambda x, t: np.array([0.0, v]) + 0*x                     # the uniform southerly wind u(x, t) [m/s]
p = np.array([0.0, 0.0])                                      # the station's position [m]
print(ch03.material_derivative_terms(F, u, p, 0.0, h=10.0))   # general stencils (h = 10 m on a 100 km scale)
print(ch03.streamwise_derivative(F, u, p, 0.0, h=10.0))       # (3.6): |u| ∂T/∂s = the advective term
x_, y_, t_ = sp.symbols('x y t')                              # symbolic twin
print(sp.simplify(ch03.material_derivative_sym(288 - sp.Rational(1, 100000)*(y_ - 10*t_), [0, 10], [x_, y_], t_)))
""", explain="""
1. `thermal_front_terms` gives the exact terms for our moving front (the pattern is carried by the wind, so the parcel's
   rate is zero); the dictionary (Ch. 1 P23) holds the local, advective and total rates.
2. `material_derivative_terms` gets the same numbers from general second-order stencils, which work for any F and u.
3. The streamwise form $\\frac{DF}{Dt}=\\frac{\\partial F}{\\partial t}+|\\mathbf u|\\frac{\\partial F}{\\partial s}$ *(3.6)*
   reproduces the advective term.
4. sympy writes $\\frac{DF}{Dt}=\\frac{\\partial F}{\\partial t}+\\mathbf u\\cdot\\nabla F$ *(3.5)* symbolically and gets 0.
""")
nb.check_agree("""
ht, hx = 1.0, 10.0                                            # time step [s] and space step [m] for the differences
local = (F(p, ht) - F(p, -ht))/(2*ht)                         # ∂T/∂t at the fixed station (central difference)
dFdy = (F(p + np.array([0, hx]), 0) - F(p - np.array([0, hx]), 0))/(2*hx)   # ∂T/∂y at the station
adv = v*dFdy                                                  # u·∇T = v ∂T/∂y
assert np.allclose([local, adv], [terms["local"], terms["advective"]], rtol=1e-8)   # stencils = exact terms
ts = np.array([-60.0, 0.0, 60.0])                             # a minute before and after [s]
path = ch03.pathline(u, [0.0, 0.0], 0.0, ts)                  # the balloon's path line (3.8) through the station
dF_along = (F(path[:, 2], ts[2]) - F(path[:, 0], ts[0]))/(ts[2] - ts[0])   # the balloon's own thermometer rate
assert abs(dF_along) < 1e-9                                   # = DT/Dt = 0: the air keeps its temperature
print("fixed-point stencils and the riding thermometer agree:", local, adv, dF_along)
""")
nb.md("Two routes, one number: stencils at a fixed point, and a thermometer riding the path line.")
nb.md(r"""
**A front of finite width.** $T(y,t)=T_0-G\,w\tanh\big((y-ct)/w\big)$: warm to the south, cold to the north, width
w = 100 km, largest gradient $G=10^{-5}$ K/m at its centre, moving north at c = 5 m/s. Its slope is
$\partial T/\partial y=-G\,\mathrm{sech}^2\big((y-ct)/w\big)$ with $\mathrm{sech}^2=1-\tanh^2$. A float starts 200 km south
and drifts north with the wind, v = 10 m/s; a station sits at y = 100 km.
""")
nb.figure(r"""
G, w, c, v = 1e-5, 1e5, 5.0, 10.0                           # largest gradient [K/m], front width [m], front speed, wind [m/s]
ys, y0 = 1e5, -2e5                                          # station at y = 100 km; float released at y = −200 km [m]
hours = np.linspace(0, 24, 241); t = hours*3600             # a day [h] and [s]
T = lambda y, t: ch03.thermal_front(0.0, y, t, G, 0.0, c, width_m=w)   # the tanh front moving north at c [K]
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4))          # (a) snapshots in space, (b) two thermometers in time
yy = np.linspace(-4e5, 5e5, 400)                            # north–south line [m]
for h, al in zip((0, 6, 12), (0.4, 0.7, 1.0)):              # three snapshots, darker = later
    a.plot(yy/1e3, T(yy, h*3600), color=COLORS["blue"], alpha=al, label=f"$T(y)$ at {h} h")    # the front at h hours
    a.plot((y0 + v*h*3600)/1e3, T(y0 + v*h*3600, h*3600), "o", color=COLORS["accent"], ms=8, alpha=al)   # the float then
a.plot(ys/1e3, T(ys, 0), "s", color=COLORS["blue"], ms=9, label="station ($y$ = 100 km)")   # the fixed station
a.plot([], [], "o", color=COLORS["accent"], label="float at those times")                   # legend entry only
a.set_xlabel("$y$ (north) [km]"); a.set_ylabel("$T$ [K]"); a.legend(fontsize=8)             # axes with units
a.set_title("A front slides north; the float overtakes it")
b.plot(hours, T(ys, t), color=COLORS["blue"], lw=2.2, label="station reads $T(y_s, t)$ → slope ∂T/∂t")        # x fixed
b.plot(hours, T(y0 + v*t, t), color=COLORS["accent"], lw=2.2, label="float reads $T(y_0+vt, t)$ → slope DT/Dt")   # label fixed
b.set_xlabel("$t$ [h]"); b.set_ylabel("$T$ [K]"); b.legend(fontsize=8)                      # axes with units
b.set_title("Two thermometers, one moving temperature pattern")
savefig(fig, "ch03", "front_station_float"); plt.show()     # save to outputs/ch03 and draw
""", see="""Left: a smooth step in temperature (warm south, cold north) sliding north, purple dots for the float catching
up with it. Right: the station's curve rises as the front passes (≈ 5.6 h), the float's curve falls as it crosses into
colder air (≈ 11 h).""",
    read="""The station sees $\\partial T/\\partial t=+Gc\\,\\mathrm{sech}^2(\\cdot)$ (the pattern moving past;
$\\mathrm{sech}^2=1-\\tanh^2$ is the slope of the tanh profile). The float sees $DT/Dt=G(c-v)\\,\\mathrm{sech}^2(\\cdot)<0$
because it outruns the pattern. Where the two curves have the same slope the advective term is zero (far from the
front, where ∇T ≈ 0).""",
    change="""…the front moved exactly with the wind (c = v): the float's curve would be flat ($DT/Dt=0$) and the station's
warming would be pure advection — the worked example above.""")
nb.plotly(r"""
G, w, c, v, y0 = 1e-5, 1e5, 5.0, 10.0, -2e5                 # the same front and float as the figure above
hrs = np.linspace(0, 24, 97)                                # time axis [h]
yf = y0 + v*hrs*3600                                        # the float's position along its path line [m]
tt = ch03.thermal_front_terms(0.0, yf, hrs*3600, 0.0, v, G, 0.0, c, width_m=w)   # exact terms along the path
loc, adv, tot = (np.asarray(tt[k])*3600 for k in ("local", "advective", "total"))  # K/s → K/h

def curves(T):                                              # the three terms up to the slider time T [h]
    k = hrs <= T + 1e-9                                     # the part of the day already travelled
    return {"∂T/∂t (local, at the float's point)": (hrs[k], loc[k]), "u·∇T (advective)": (hrs[k], adv[k]),
            "DT/Dt (the float's rate)": (hrs[k], tot[k]), "DT/Dt, whole day (ghost)": (hrs, tot)}   # name → (x, y)

steps = np.linspace(0, 24, 25 if not FAST else 13)          # slider positions [h]
fig = slider_figure(curves, "t", steps, unit="h", xlabel="time [h]", ylabel="rate [K/h]", xrange=[0, 24],
                    title="The three terms along the float: the balance shifts as it crosses the front")   # precomputed
recolor(fig, {"∂T/∂t (local, at the float's point)": COLORS["blue"], "u·∇T (advective)": COLORS["amber"],   # blue, amber
              "DT/Dt (the float's rate)": COLORS["accent"], "DT/Dt, whole day (ghost)": COLORS["grid"]},     # purple, grey
        dashes={"DT/Dt, whole day (ghost)": "dot"})         # the ghost dotted
fig.show()                                                  # draw it
""", explain="""
1. The float's path line is $y=y_0+vt$; `thermal_front_terms` evaluates the exact local, advective and total rates at
   each of its positions (arrays in, arrays out).
2. `slider_figure` precomputes the curves for every slider time, so the figure works on the web page without Python.
3. `recolor` gives the curves the notebook's colours: blue local, amber advective, purple total.
""")
see_read_change(
    "Three curves drawn up to the slider time; all three peak near 11 h, when the float is at the front's centre.",
    """The purple curve is the sum of the blue and the amber ones at every instant, $\\frac{DT}{Dt}=\\frac{\\partial T}{\\partial t}+\\mathbf u\\cdot\\nabla T$
*(Eq. 3.5)*: blue is positive (the pattern warms any fixed point it passes), amber is more negative (the float runs into
colder air twice as fast as the pattern moves), so purple is negative.""",
    "…you stopped the slider at t ≈ 11 h (the crossing): all three curves are at their largest values there. With c = v the purple curve would lie on zero.")
nb.explainer("material_derivative_probe", "Why does the station warm while the air does not?", r"""
**Why interactive:** one static front shows one wind direction and one speed; only by turning the wind yourself do you see the amber term vanish along the isotherms and the float's rate change sign as it outruns the pattern.
Two observers read one moving temperature pattern on one clock: a fixed probe ($\partial T/\partial t$) and a float
carried by the wind ($DT/Dt$). You set the wind, the gradient and any heating; the term bars show the local and
advective parts adding to the material derivative $\frac{DT}{Dt}=\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T$
*(Eq. 3.5)*, and the float's measured rate lands on the purple bar.
""", tries=[
    "Pick the preset 'pure advection': the float's reading stays flat while the probe warms — read the numbers in Explain.",
    "Turn the wind until it blows along the isotherms: the amber bar vanishes (u ⟂ ∇T).",
    "Switch to the stretching-map mode (x = X e^{αt}): now F is the velocity itself and the float's rate is its acceleration α²x.",
    "Open the Derivation tab and step through D02 with your wind speed.",
])
whatif(r"""
…F were a velocity component instead of a temperature? Then D/Dt of u is the particle's acceleration,
$\frac{D\mathbf u}{Dt}=\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u$ — the left side of Newton's
law in Ch. 4. Section 3.3 asks what these rates look like as curves in the flow, and whether the split into local and
advective parts depends on who is watching (C05).
""")

# =====================================================================================================================
# A.3 §3.3 Flow Lines, Fluid Acceleration, and Galilean Transformation — C03, C04, C05
# =====================================================================================================================
nb.section("3.3", "Flow Lines, Fluid Acceleration, and Galilean Transformation", intro="""
**What is this section about?** Three curves that picture a flow — the instantaneous direction (streamline), the track
of one particle (path line) and the line of dye from a fixed port (streak line) — and the proof that a particle's
acceleration is the same for every observer moving at constant velocity, even though "steady" is not.
""")
core("C03", "Streamlines: the flow's direction at one frozen instant (3.7)", r"""
Freeze the flow for an instant and draw curves that follow the velocity arrows everywhere. What equations do those
curves obey, and why do they change from one instant to the next?
""", eqs=("3.7",))
nb.md(r"""
#### The problem in plain words

A satellite image of cloud streaks, or iron filings along a magnet's field lines, shows *directions at one instant*.
A weather chart's wind streamlines are the same idea. To draw one we need a rule: at every point of the curve, its
tangent points along $\mathbf u$ at that instant.
""")
nb.md(r"""
#### The idea

```
freeze the clock at t   →   field of arrows u(x, t)   →   curves everywhere tangent to the arrows
ds ∥ u   ⇔   dx/u = dy/v = dz/w   ⇔   u × ds = 0            (unsteady flow: a new picture every instant)
```
""")
P("P92", "parametric curves, tangent vector and arc length", r"""
A curve can be written $\mathbf x(s)$ with a parameter s; its tangent is $d\mathbf x/ds$. If s is the **arc length**
(distance measured along the curve), the tangent has length 1. To draw a curve tangent to $\mathbf u$, march along
$d\mathbf x/ds=\mathbf u/|\mathbf u|$.
""", code="""
s = np.linspace(0, 2*np.pi, 5)                       # arc length along the unit circle [m]
x, y = np.cos(s), np.sin(s)                          # the curve x(s)
tx, ty = -np.sin(s), np.cos(s)                       # its tangent dx/ds
print(np.hypot(tx, ty))                              # [1. 1. 1. 1. 1.]: unit length, because s is arc length
""")
P("P93", "parallel vectors and the cross-product test", r"""
Two vectors are parallel when one is a multiple of the other, $\mathbf a=\lambda\mathbf b$ — equivalently when
$\mathbf a\times\mathbf b=0$ (the cross product, Ch. 2 §2.7, measures the area they span; in components
$(\mathbf a\times\mathbf b)_i=\varepsilon_{ijk}a_jb_k$). The ratio test $a_x/b_x=a_y/b_y=a_z/b_z$ says the same but fails
when a component of $\mathbf b$ is zero; the cross product never divides.
""", code="""
a, b = np.array([2.0, 4.0, 0.0]), np.array([1.0, 2.0, 0.0])   # a = 2b
print(np.cross(a, b), a[0]/b[0], a[1]/b[1])                   # [0 0 0] 2.0 2.0: parallel (the z-ratio would be 0/0)
""")
D("D03", ref="3.7")
note("N14", "The streamline equations", r"""
are the result just derived. In a plane flow the first equality gives the slope $dy/dx=v/u$. Integrating from many
starting points, upstream and downstream, fills the picture. Number: where $\mathbf u=(1,2)$ m/s the streamline climbs
with slope 2.
""", equation=r"dx/u=dy/v=dz/w", ref="3.7")
nb.worked_example("slopes and a circle", r"""
1. At a point where u = 1 m/s, v = 2 m/s: $dy/dx=v/u=2$ (63.4° to the x-axis).
2. Solid-body rotation u = −y, v = x (s⁻¹ × m): $dy/dx=v/u=-x/y$, so $y\,dy=-x\,dx$ and $x^2+y^2=\text{const}$ — circles.
3. Ex. 3.1 at the instant t′ = π/4 s (ω = 1 s⁻¹): $v/u=\tan(\omega t')=1$ at every point, so every streamline is a
   straight line at 45°.
""")
P("P94", "solve_ivp options: t_eval, dense_output, events, backward integration", r"""
Ch. 1 used `solve_ivp` with its defaults (P31). Here: `t_eval` asks for the solution at chosen times;
`dense_output=True` returns a continuous solution `sol.sol(t)`; an `events` function stops the integration where it
crosses zero (we stop a streamline at a stagnation point, $|\mathbf u|=0$, where the direction is undefined); a time
span that runs backwards (`t_span=(0, -5)`) integrates into the past — how a streamline is drawn upstream and a streak
particle is traced back to the port.
""", code="""
from scipy.integrate import solve_ivp                           # adaptive ODE solver (Ch. 1 P31)
stop = lambda t, y: y[0] - 0.5; stop.terminal = True           # event: stop when y reaches 0.5
sol = solve_ivp(lambda t, y: -y, (0, 5), [1.0], events=stop, dense_output=True, rtol=1e-8)   # dy/dt = −y, y(0) = 1, tight tolerance
print(sol.t_events[0], sol.sol(0.3))                           # [0.693] (= ln 2) and y(0.3) = e^-0.3 ≈ 0.741
""")
nb.code("""
u31 = ch03.preset_field("ex31", omega=1.0, xi0=1.0)          # Ex. 3.1: u = ωξ_o cos ωt, v = ωξ_o sin ωt as a callable u(x, t)
for tp in (0.0, np.pi/4, np.pi/2):                            # three drawing instants t' [s]
    sl = ch03.streamline(u31, [0.0, 0.0], tp, s_max=2.0)      # dx/ds = u/|u| at frozen t', 2 m each way from the origin
    print(f"t' = {tp:.3f} s: streamline ends at {np.round(sl[:, -1], 3)}")
print(ch03.streamline_slope(u31, np.array([0.3, -0.7]), np.pi/4))   # (3.7) in a plane: dy/dx = v/u at any point
usb = lambda x, t: np.array([-x[1], x[0]])                    # solid-body rotation u = (−y, x) [m/s]
c = ch03.streamline(usb, [1.0, 0.0], 0.0, s_max=2*np.pi)      # the streamline through (1, 0)
print(np.ptp(np.hypot(c[0], c[1])))                           # spread of its radius: ≈ 0 → a circle
""", explain="""
1. `preset_field("ex31")` is the Ex. 3.1 field as a callable $\\mathbf u(\\mathbf x,t)$ (the same at every point).
2. `streamline` integrates $d\\mathbf x/ds=\\mathbf u/|\\mathbf u|$ with the clock frozen at t′, both directions from the
   seed (D03 step 5); the end points (2, 0), (1.414, 1.414), (0, 2) are straight lines at 0°, 45°, 90°.
3. `streamline_slope` is the first equality of $dx/u=dy/v=dz/w$ *(3.7)*: the slope v/u = 1 at t′ = π/4.
4. For solid-body rotation the computed streamline keeps its radius (spread ≈ 10⁻⁹ m): a circle.
""")
nb.figure(r"""
xg = np.linspace(-2, 2, 21); X, Y = np.meshgrid(xg, xg)      # a 21 × 21 grid [m] (Ch. 2 P76)
fig, axs = plt.subplots(1, 4, figsize=(13, 3.5))
for ax, tp in zip(axs[:3], (0.0, np.pi/4, np.pi/2)):          # three frozen instants t' [s]
    U, V = ch03.unsteady_flow_preset("ex31", X, Y, tp, xi0=1.0, omega=1.0)   # the Ex. 3.1 field at t'
    ax.streamplot(X, Y, U + 0*X, V + 0*Y, color=COLORS["teal"], density=0.7)   # streamlines (Ch. 2 P78)
    ax.plot(0, 0, "s", color=COLORS["ink"])                   # the port at the origin
    ax.quiver(0, 0, float(np.mean(U)), float(np.mean(V)), color=COLORS["ink"], scale=3)   # velocity at the port
    ax.set_title(f"Ex. 3.1, $t'$ = {tp:.2f} s"); ax.set_aspect("equal")
Us, Vs = -Y, X                                                # solid-body rotation for comparison
axs[3].streamplot(X, Y, Us, Vs, color=COLORS["teal"], density=0.8)
axs[3].set_title("solid-body rotation"); axs[3].set_aspect("equal")
for ax in axs: ax.set_xlabel("$x$ [m]")
axs[0].set_ylabel("$y$ [m]")
savefig(fig, "ch03", "streamlines_ex31"); plt.show()
""", see="Parallel teal lines that turn by 45° from panel to panel; circles in the last panel.",
    read="""Ex. 3.1's velocity is the same at every point at a given instant, so its streamlines are parallel straight lines;
only their direction changes with time — the pattern turns once per period $2\\pi/\\omega$. Solid-body rotation is steady:
its circles never change.""",
    change="…the flow were steady (freeze ωt′): all three panels would be identical — and so, C04 shows, would path and streak lines.")
nb.plotly(r"""
def lines(tp):                                               # the three Ex. 3.1 curves through the origin at t' [s]
    d = ch03.example_3_1(tp, xi0=1.0, omega=1.0, n=120)      # closed forms (D04, D05)
    return {"streamline": tuple(d["streamline"]), "path line": tuple(d["pathline"]), "streak line": tuple(d["streakline"])}   # (x, y) each

fig = slider_figure(lines, "t'", np.linspace(0, 2*np.pi, 24 if not FAST else 12), unit="s", xlabel="x [m]",
                    ylabel="y [m]", xrange=[-2.5, 2.5], yrange=[-2.5, 2.5], height=520,          # fixed axes [m]
                    title="Streamlines turn with t'; path and streak lines are circles that roll round the port")
recolor(fig, {"streamline": COLORS["teal"], "path line": COLORS["orange"], "streak line": COLORS["rose"]})   # house colours
fig.update_yaxes(scaleanchor="x", scaleratio=1)              # equal scales, so circles look round
fig.show()                                                   # draw it
""", explain="""
1. `example_3_1(t′)` returns the closed-form streamline $y=x\\tan\\omega t'$ and the two circles of D04 and D05 at the
   drawing instant t′.
2. `slider_figure` precomputes them for 24 instants over one period.
""")
see_read_change(
    "A teal line and two circles, orange and rose, all passing through the origin.",
    "At every t′ the teal line touches both circles at the origin — the tangency shown in D05. The circles sit on opposite sides of the port.",
    "…drag t′ by a quarter period, π/2 ≈ 1.57 s: all three turn by 90° together.")
note("N15", "Stream tube", r"""
The streamlines through every point of a closed curve C form a tube (Fig. 3.6). No fluid crosses its wall, because the
wall is everywhere tangent to $\mathbf u$; in a steady flow the same volume per second passes every cross-section.
Number below: the axisymmetric straining flow $\mathbf u=(-x/2,-y/2,z)$ s⁻¹ ($\nabla\cdot\mathbf u=0$) with a tube
seeded on a circle of radius 0.5 m at z = 1 m: the tube narrows as $R(z)=0.5\sqrt{1/z}$ and the flux through every
section is $\pi/4=0.785$ m³/s.
""")
nb.plotly(r"""
u3 = lambda x, t: np.array([-x[0]/2, -x[1]/2, x[2]])          # axisymmetric strain [1/s × m]: ∇·u = −½ − ½ + 1 = 0
nl = 12 if not FAST else 8                                    # number of streamlines in the tube
fig = go.Figure()
for k, ph in enumerate(np.linspace(0, 2*np.pi, nl, endpoint=False)):   # seeds on the circle R = 0.5 m at z = 1 m
    c = ch03.streamline(u3, [0.5*np.cos(ph), 0.5*np.sin(ph), 1.0], 0.0, s_max=1.3)   # (3.7), both directions
    fig.add_trace(go.Scatter3d(x=c[0], y=c[1], z=c[2], mode="lines", line=dict(color=COLORS["teal"], width=4),
                               showlegend=(k == 0), name="streamlines of the tube wall"))   # 3-D lines (Ch. 2 P64)
for z0, R0, col in ((1.0, 0.5, COLORS["blue"]), (2.0, 0.5/np.sqrt(2), COLORS["orange"])):   # two cross-sections
    ph = np.linspace(0, 2*np.pi, 41)
    xs = np.r_[0, R0*np.cos(ph)]; ys = np.r_[0, R0*np.sin(ph)]; zs = np.full(42, z0)   # centre + rim points
    fig.add_trace(go.Mesh3d(x=xs, y=ys, z=zs, i=np.zeros(40, int), j=np.arange(1, 41), k=np.arange(2, 42),
                            color=col, opacity=0.45, name=f"section at z = {z0:g} m"))   # a fan of triangles = a disc
q1 = ch03.flux_through_disc(u3, [0, 0, 1], [0, 0, 1], 0.5)             # ∫ u·n dA through the lower section [m³/s]
q2 = ch03.flux_through_disc(u3, [0, 0, 2], [0, 0, 1], 0.5/np.sqrt(2))  # … through the upper, narrower section
print(f"flux through z = 1 m: {q1:.4f} m^3/s, through z = 2 m: {q2:.4f} m^3/s (pi/4 = {np.pi/4:.4f})")
fig.update_layout(height=480, title="A stream tube narrows where the flow speeds up; the flux stays the same",
                  scene=dict(aspectmode="data", xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]"))
fig.show()
""", explain="""
1. Twelve streamlines started on a circle form the wall of a tube.
2. Two discs (drawn as fans of triangles) cut the tube at z = 1 m and z = 2 m, where its radius is 0.5 and 0.354 m.
3. `flux_through_disc` integrates $\\mathbf u\\cdot\\mathbf n\\,dA$ over each disc (midpoint rule): both 0.7854 m³/s.
""")
see_read_change(
    "Teal lines converging upward into a narrower tube; a blue disc below and a smaller orange disc above.",
    "The tube narrows where the flow speeds up along z; the two discs carry the same flux because the wall lets nothing through.",
    "…the seed circle were wider (radius 1 m): a fatter tube carrying four times the flux (∝ πR²), still the same at every section.")
whatif(r"""
…you followed one particle instead of freezing the clock? In a steady flow it would run along a streamline. In Ex. 3.1
the arrows turn while the particle moves, so its track is not any of the straight teal lines — it is a circle (C04).
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C04", "Path lines and streak lines: where one particle goes, where the dye sits (3.8)", r"""
Dye leaks steadily from a port on the floor of a sloshing tank, and you track one speck of it. Why does the photograph
of the dye not show the path of the speck?
""", eqs=("3.8",))
nb.md(r"""
#### The problem in plain words

Every flow picture you see in a lab is one of three things: a long exposure of one particle (a **path line**), a
snapshot of all the dye that came out of one port (a **streak line**), or a snapshot of directions (a **streamline**).
Near a beach, the water under long swell sloshes back and forth; a drop of dye released there and a floating leaf trace
quite different curves. We need both as computable objects.
""")
nb.md(r"""
#### The idea

| line | what is fixed | what varies along it | how you see it |
|---|---|---|---|
| streamline | the instant t | position (arc length s) | arrows at one instant |
| path line | the particle (label $\mathbf r_o$, $t_o$) | time t | long exposure of one speck |
| streak line | the port $\mathbf x_o$ and the instant t | release time $t_o$ | photo of continuously injected dye |

**Steady flow: all three are the same curve. Unsteady: three different curves.**
""")
note("N16", "Path line", r"""
A **path line** is the trajectory of one particle of fixed identity, the Lagrangian $\mathbf x=\mathbf r(t;\mathbf r_o,t_o)$
of C01 drawn in space. Given only the Eulerian velocity field, we get it by solving the ODE below.
""")
nb.md(r"""
#### The maths — the path-line equation

Each particle moves with the velocity of the field at the place it currently occupies:

$$\frac{d\mathbf r}{dt}=[\mathbf u(\mathbf x,t)]_{\mathbf x=\mathbf r}=\mathbf u(\mathbf r,t),\qquad \mathbf r(t_o)=\mathbf r_o\qquad(3.8)$$

— an initial-value problem: three coupled ODEs, one per component. (A discretised version of (3.8) is how **particle
image velocimetry**, PIV, turns pairs of photographs of seeded particles into velocities.)
""")
note("N17", "Streak line", r"""
A **streak line** through the port $\mathbf x_o$ at time t is the set of particles that passed through $\mathbf x_o$
at earlier times $t_o$: solve $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ *(Eq. 3.8)* once per release time with
$\mathbf r(t_o)=\mathbf x_o$, then plot all positions at the same t. The parameter along it is the release time,
$x_i=r_i(t;\mathbf x_o,t_o)$. Sometimes $t_o$ can be eliminated to give one equation for the curve (D05).
""")
note("N18", "Ex. 3.1", r"""
(a spatially uniform, time-periodic flow — a caricature of the back-and-forth motion under long surface waves, not a
wave solution; real wave orbits shrink with depth, Ch. 7): $u=\omega\xi_o\cos\omega t$, $v=\omega\xi_o\sin\omega t$.
Find the three lines through the origin at the instant t = t′. ⚠️ Four names that look alike: **t′** the drawing
instant, **t_o** a release time, **c_x, c_y** integration constants (the book's x_o, y_o), **ξ_o** the amplitude [m].
The streamline comes straight from $dx/u=dy/v=dz/w$ *(Eq. 3.7)*, integrated through the origin (the constant is 0):
""", equation=r"\frac{dy}{dx}=\frac{v}{u}=\frac{\omega\xi_o\sin(\omega t')}{\omega\xi_o\cos(\omega t')}=\tan(\omega t')\ \Rightarrow\ y=x\tan(\omega t')")
D("D04")
D("D05")
nb.worked_example("Ex. 3.1 with ω = 1 s⁻¹, ξ_o = 1 m, t′ = 0", r"""
1. Streamline: $y=x\tan0=0$ — the x-axis.
2. Path line (D04): $(x+\sin0)^2+(y-\cos0)^2=1$, i.e. $x^2+(y-1)^2=1$: a unit circle centred at (0, 1).
3. Streak line (D05): $(x-\sin0)^2+(y+\cos0)^2=1$, i.e. $x^2+(y+1)^2=1$: a unit circle centred at (0, −1).
4. All three touch the origin with slope $\tan0=0$ (horizontal): the circles lie on opposite sides of the x-axis, as in
   the book's Fig. 3.7.
5. A quarter period later (t′ = π/2 s) everything has turned by 90°: the streamline is the y-axis, the centres are
   (−1, 0) and (1, 0).
""")
P("P95", "RK4 by hand", r"""
The classical fourth-order Runge–Kutta step for $d\mathbf r/dt=\mathbf f(\mathbf r,t)$ samples the slope four times per
step: $k_1=f(r,t)$, $k_2=f(r+\tfrac{h}{2}k_1,t+\tfrac h2)$, $k_3=f(r+\tfrac h2k_2,t+\tfrac h2)$, $k_4=f(r+hk_3,t+h)$, then
$r\leftarrow r+\tfrac h6(k_1+2k_2+2k_3+k_4)$. Its error per unit time falls like $h^4$ (Euler's explicit step, Ch. 1 P30,
like h). The explainers use exactly this step in JavaScript.
""", code="""
f = lambda y, t: -y                                  # dy/dt = −y, exact y = e^−t
y, h = 1.0, 0.1                                      # start value and step [s]
for n in range(10):                                  # ten steps to t = 1 s
    k1 = f(y, n*h); k2 = f(y + h/2*k1, n*h + h/2); k3 = f(y + h/2*k2, n*h + h/2); k4 = f(y + h*k3, n*h + h)
    y += h/6*(k1 + 2*k2 + 2*k3 + k4)                 # the weighted average of the four slopes
print(y, np.exp(-1.0))                               # 0.3678798 vs 0.3678794: 5e-7 after ten steps
""")
nb.code("""
d = ch03.example_3_1(0.0, xi0=1.0, omega=1.0)                # closed forms at t' = 0 (D04, D05)
print(d["path_center"], d["streak_center"], d["slope"])      # centres (0, 1) and (0, −1) m, slope tan 0 = 0
T = 2*np.pi                                                  # the period 2π/ω [s]
ts = np.linspace(-T/2, T/2, 201)                             # half a period back and forward [s]
p = ch03.pathline(u31, [0.0, 0.0], 0.0, ts)                  # (3.8): dr/dt = u(r, t), r(0) = 0
print("path residual:", np.max(np.abs(np.hypot(p[0] - d["path_center"][0], p[1] - d["path_center"][1]) - 1.0)))
tr = np.linspace(-T, 0.0, 181 if not FAST else 91)           # dye released over the last period [s]
s = ch03.streakline(u31, [0.0, 0.0], 0.0, tr)                # positions at t' = 0 of every release
print("streak residual:", np.max(np.abs(np.hypot(s[0] - d["streak_center"][0], s[1] - d["streak_center"][1]) - 1.0)))
""", explain="""
1. The closed forms of D04 and D05.
2. `pathline` solves $d\\mathbf r/dt=\\mathbf u(\\mathbf r,t)$ *(3.8)* for one particle (t_eval may lie before $t_0$);
   every point is 1 m from the centre (0, 1).
3. `streakline` solves (3.8) once per release time (one vector ODE for the whole ensemble) and collects the positions at
   t′ = 0: a circle about (0, −1).
""")
nb.check_agree("""
r = np.array([0.0, 0.0]); h = T/400                          # start at the port; 400 steps per period
f = lambda r, t: u31(r, t)                                   # the right side of (3.8)
for n in range(200):                                         # 200 RK4 steps = half a period (primer P95)
    t_ = n*h
    k1 = f(r, t_); k2 = f(r + h/2*k1, t_ + h/2); k3 = f(r + h/2*k2, t_ + h/2); k4 = f(r + h*k3, t_ + h)
    r = r + h/6*(k1 + 2*k2 + 2*k3 + k4)
lib_end = ch03.pathline(u31, [0.0, 0.0], 0.0, [0.0, T/2])[:, -1]   # the library's answer at t = T/2
assert np.allclose(r, lib_end, atol=1e-8)                    # hand-written RK4 = library
assert np.allclose(r, [np.sin(T/2), 1 - np.cos(T/2)], atol=1e-8)   # = D04's circle (x = sin t, y = 1 − cos t)
print("RK4 by hand:", r, " library:", lib_end)
""")
nb.md("Twenty lines of arithmetic reproduce the library to 10⁻⁸ — no magic inside `pathline`.")
note("N19", "Fig. 3.7, our drawing", r"""
the three lines of Ex. 3.1 at ωt′ = 30°, computed from the closed forms (curves) and by the ODE solvers (dots).
""")
nb.figure(r"""
from scripts.ch03_drawings import flow_lines_figure          # drawing helper (no physics inside)
fig, ax, info = flow_lines_figure(t_prime=np.pi/6, xi0=1.0, omega=1.0)   # ωt' = 30°
ax.plot([-2, 2], [-2*np.tan(np.pi/6), 2*np.tan(np.pi/6)], ":", color=COLORS["muted"], lw=1)   # the common tangent
print({k: f"{v:.1e} m" for k, v in info.items()})            # numerical vs closed form: max distance
savefig(fig, "ch03", "fig3_7_flow_lines"); plt.show()
""", see="One straight teal line and two circles of equal radius on opposite sides, all touching at the origin (dots: the numerical curves).",
    read="""Same point, same instant, three different curves — because the flow is unsteady. The common tangent at the
origin (grey dotted) is the velocity there at t′ (D05's last steps). The path circle is run counterclockwise.""",
    change="…t′ larger by π/(2ω): the whole picture turns by 90° counterclockwise; the circles keep their radius ξ_o.")
nb.animation(r"""
from scripts.ch03_drawings import flow_lines_animation       # dye + one particle + turning streamline over a period
anim = flow_lines_animation(frames=60 if not FAST else 30, xi0=1.0, omega=1.0, n_dye=60)   # one period of Ex. 3.1
show_animation(anim, player="video")                         # smooth MP4 (P16)
""", explain="""
1. For each frame the streamline through the port is redrawn at the current instant (teal), the tagged particle
   released at t = 0 advances along its path line (orange), and all dye released so far is placed with `streakline` (rose).
""")
see_read_change(
    "A teal line turning steadily, an orange particle drawing a circle directly above the port, and rose dye on a circle that swings round the port.",
    """The tagged particle (released at t = 0) keeps to its own circle, centred at (0, 1) m directly above the port. The
dye lies on the mirror image of the path circle of *the particle at the port now* (D05): its centre
$(\\xi_o\\sin\\omega t,\\,-\\xi_o\\cos\\omega t)$ goes once round the port per period, so at t = T/2 the dye circle coincides
with the tagged particle's circle. After one period the particle is back at the port and the pattern repeats.""",
    """…there were a mean current (the explainer's '+ current' mode): the particle would not come back; its path becomes a
looping curve drifting downstream and the dye forms a wave.""")
note("N13", "Steady flow ⇒ the three lines coincide", r"""
If $\mathbf u$ does not depend on t, a particle at a point moves along the streamline through it and every later
particle from the same port follows the same track. Check on the steady hyperbolic flow $\mathbf u=(x,-y)$ s⁻¹ from the
port (1, 1) m, whose streamlines are the hyperbolas $xy=$ const:
""")
nb.code("""
us = ch03.preset_field("rotating_strain", s=1.0, Omega=0.0)        # Ω = 0: the steady strain u = (x, −y) [1/s × m]
x0 = [1.0, 1.0]                                                    # the port [m]
c_stream = ch03.streamline(us, x0, 0.0, s_max=1.0)                 # (3.7) at t = 0
c_path = ch03.pathline(us, x0, 0.0, np.linspace(0, 1, 50))         # (3.8) for the particle leaving at t = 0
c_streak = ch03.streakline(us, x0, 1.0, np.linspace(0, 1, 50))     # dye released 0…1 s, seen at t = 1 s
for name, c in (("streamline", c_stream), ("path line", c_path), ("streak line", c_streak)):
    print(f"{name:12s} max |xy − 1| = {np.max(np.abs(c[0]*c[1] - 1)):.1e}")   # distance from the hyperbola xy = 1
""", explain="""
In a steady flow all three curves lie on the same hyperbola $xy=1$ (residuals ≈ 10⁻⁹). Flow-visualisation atlases (Van
Dyke's *Album of Fluid Motion*) are full of such pictures; in unsteady flow you must know which of the three a
photograph shows.
""")
nb.explainer("flow_lines_unsteady", "Three lines through one point — why do they disagree?", r"""
**Why interactive:** a still picture shows three curves at one instant; the running clock shows *how* they are made — the streamline redrawn, the particle's trail growing, the dye laid down — and lets you switch unsteadiness off.
A clock runs: the streamline pattern is redrawn every instant, one particle leaves its trail and a port leaks dye. You
change the period, the amplitude and a mean current, or switch to a steady flow — and watch the three curves separate
or collapse onto one.
""", tries=[
    "Press ▶ with the Ex. 3.1 preset and stop at ωt = 30°: compare with our Fig. 3.7 above.",
    "Choose 'steady': the three curves become one line — the N13 statement as an experiment.",
    "Add a mean current U₀ = 0.5 m/s: the particle no longer returns to the port; the dye forms a wave.",
    "Click a dye particle (inspector) to see when it left the port and the path it took.",
])
whatif(r"""
…you watched Ex. 3.1 from a boat drifting with a mean current? The streamlines would look different again — whether a
flow is steady depends on the observer. The next block asks whether anything important depends on the observer: the
particle's acceleration does not.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C05", "Galilean invariance: the acceleration does not care who is watching (3.9)", r"""
A cylinder is towed through a still lake; a pier stands in a steady river. One flow is unsteady, the other steady, yet
they are the same flow seen from two boats. Does a water particle's acceleration depend on the boat you watch it from?
""", eqs=("3.9",))
nb.md(r"""
#### The problem in plain words

Wind-tunnel engineers hold a model still and blow air past it; the real car moves through still air. Oceanographers
describe a wave in a frame moving with the crest, where the flow is steady. All of this is legitimate only if Newton's
law — force equals mass times *acceleration* — reads the same in both frames. We check that the acceleration, written as
local + advective, is frame-independent for observers moving at constant velocity.
""")
nb.md(r"""
#### The idea

```
body frame (pier fixed):   ∂u/∂t = 0 (steady)          (u·∇)u ≠ 0     ─┐
fluid frame (lake fixed):  ∂u'/∂t' ≠ 0 (unsteady)      (u'·∇')u'       ├─ same particle, same total Du/Dt
only U (constant) separates them: the split moves, the sum does not     ─┘
```
""")
P("P96", "frames of reference and relative velocity", r"""
A frame of reference is a set of axes plus a clock. Two frames whose axes stay parallel and whose origins separate at a
constant velocity $\mathbf U$ are related by a *Galilean transformation*: positions differ by $\mathbf Ut$ (+ a fixed
offset), velocities by $\mathbf U$, and both use the same clock (no relativity here). A frame that rotates is different:
its points move with $\boldsymbol\Omega\times\mathbf x$, which is not uniform.
""", code="""
U = np.array([2.0, 0.0])                             # frame O' moves at 2 m/s along x
u_lab = np.array([2.5, 0.3])                         # a particle's velocity seen from O [m/s]
print(u_lab - U)                                     # [0.5 0.3]: the same particle seen from O'
""")
note("N04", "Fig. 3.2, one flow in two frames", r"""
Ideal flow past a circular cylinder of radius a in a stream U (a uniform stream plus a *doublet*; derived in Ch. 6
§6.3): $u=U[1-a^2(x^2-y^2)/r^4]$, $v=-2Ua^2xy/r^4$. Seen from the cylinder the streamlines are fixed (steady). Seen from
the still fluid, $\mathbf u'=\mathbf u-\mathbf U$: streamlines start and end on the moving body and change as it passes
(unsteady). The two velocity fields differ by the constant vector $\mathbf U$:
""", equation=r"\mathbf u=\mathbf U+\mathbf u'")
nb.figure(r"""
x = np.linspace(-3, 3, 121); y = np.linspace(-2, 2, 81); X, Y = np.meshgrid(x, y)   # 121 × 81 grid [m]
th = np.linspace(0, 2*np.pi, 100)                            # to draw the cylinder's outline
P_ = (0.0, 1.5)                                              # the tagged point P [m]
fig, axs = plt.subplots(1, 2, figsize=(11, 3.9))             # (a) body frame, (b) fluid frame
for ax, frame, Uf, ttl in zip(axs, ("body", "fluid"), (1.0, 0.0), ("(a) body frame: steady", "(b) fluid frame at t = 0: unsteady")):
    Uu, Vv = ch03.cylinder_flow(X, Y, 1.0, 1.0, frame=frame)   # U = 1 m/s, a = 1 m; NaN inside the body
    ax.streamplot(X, Y, np.ma.masked_invalid(Uu), np.ma.masked_invalid(Vv), color=COLORS["teal"], density=1.1, linewidth=0.9)   # NaN masked
    ax.fill(np.cos(th), np.sin(th), color=COLORS["grid"], zorder=3)   # the cylinder
    r = ch03.frame_acceleration_terms(*P_, 1.0, 1.0, Uf)      # velocity and acceleration at P in this frame
    ax.quiver(*P_, *r["u"], color=COLORS["ink"], scale=4, width=0.006, zorder=4, label="velocity at P")   # black arrow
    ax.quiver(*P_, *r["total"], color=COLORS["accent"], scale=4, width=0.008, zorder=4,
              label=f"acceleration ({r['total'][0]:.3f}, {r['total'][1]:.3f}) m/s²")   # purple arrow, same scale
    ax.set_aspect("equal"); ax.set_title(ttl); ax.set_xlabel("$x$ [m]"); ax.legend(fontsize=7, loc="lower left")   # equal axes
axs[0].set_ylabel("$y$ [m]")                                  # shared y label
savefig(fig, "ch03", "cylinder_two_frames"); plt.show()       # save to outputs/ch03 and draw
""", see="Two completely different streamline pictures and one identical purple arrow at P = (0, 1.5) m.",
    read="""(a) is steady — the arrows never change at a fixed point; (b) is unsteady — the body moves through the water
(towards −x) and the loops leave its front and re-enter its back. The black velocity arrows differ by exactly U; the
purple acceleration arrow, (0, −0.856) m/s², is the same in both.""",
    change="…the observer moved at half the towing speed: a third streamline picture, the same purple arrow (the observer slider below and the Galilean-frames explainer).")
note("N21", "The Galilean transformation (Fig. 3.8)", r"""
Frame O′x′y′z′ moves at constant velocity $\mathbf U$ relative to Oxyz with parallel axes; $\mathbf x'_o$ is the vector
from O to O′ at t = 0. `ch03.galilean_transform(u, U)` returns the primed field
$\mathbf u'(\mathbf x',t')=\mathbf u(\mathbf x'+\mathbf Ut'+\mathbf x'_o,t')-\mathbf U$.
""", equation=r"\mathbf u(\mathbf x,t)=\mathbf U+\mathbf u'(\mathbf x',t'),\qquad t=t',\qquad \mathbf x=\mathbf x'+\mathbf Ut+\mathbf x'_o")
P("P97", "chain rule with a moving frame", r"""
If $g(x,t)=f(x-Ut,t)$, then at fixed x,
$\frac{\partial g}{\partial t}=\frac{\partial f}{\partial t'}-U\frac{\partial f}{\partial x'}$: holding x fixed while t
advances means the moving coordinate $x'=x-Ut$ *decreases*. A time derivative "at fixed position" depends on whose
position is held fixed — the classic trap of this section.
""", code="""
x, t, U = sp.symbols('x t U')                        # lab position, time, frame speed
xp, tp = sp.symbols('xp tp')                         # moving-frame coordinates x', t'
f = sp.sin(xp)*sp.exp(-tp)                           # any f(x', t')
g = f.subs({xp: x - U*t, tp: t})                     # the same field seen in the unprimed frame
rhs = (sp.diff(f, tp) - U*sp.diff(f, xp)).subs({xp: x - U*t, tp: t})   # ∂f/∂t' − U ∂f/∂x'
print(sp.simplify(sp.diff(g, t) - rhs))              # 0
""")
_d06 = PF["D06"]["code"]
for _old, _new in (
        ("# f, g are undefined functions; sympy 1.14 keeps two Subs(Derivative) forms of the same ∂f/∂t′, so test on",
         "# f, g are unknown functions and sympy cannot always simplify their derivatives, so test on"),
        ("cs = sp.symbols('c0:20'); X, Y, T = sp.symbols('X Y T')",
         "cs = sp.symbols('c0:20'); X, Y, T = sp.symbols('X Y T')   # 20 coefficients and three dummy arguments"),
        ("mons = [X**i * Y**j * T**k for i in range(4) for j in range(4) for k in range(4) if i + j + k <= 3]",
         "mons = [X**i * Y**j * T**k for i in range(4) for j in range(4) for k in range(4) if i + j + k <= 3]   # the 20 monomials of degree ≤ 3"),
        ("repl = {f: sp.Lambda((X, Y, T), sum(c * m for c, m in zip(cs, mons))),",
         "repl = {f: sp.Lambda((X, Y, T), sum(c * m for c, m in zip(cs, mons))),   # f → a generic cubic"),
        ("g: sp.Lambda((X, Y, T), sum(c * m for c, m in zip(cs[::-1], mons)))}",
         "g: sp.Lambda((X, Y, T), sum(c * m for c, m in zip(cs[::-1], mons)))}   # g → another one")):
    _d06 = _d06.replace(_old, _new)
D("D06", ref="3.9", check_src=_d06)
note("N20", "Back to Fig. 3.2", r"""
In the body frame the unsteady term is zero but the streamlines bend, so particles still accelerate — all of it
advective. In the fluid frame both terms are nonzero. Because the frames differ by a constant $\mathbf U$,
$\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$
*(Eq. 3.9)* guarantees the same acceleration at the same place relative to the cylinder.
""")
note("N23", "The split is observer-dependent, the sum is not", r"""
— the term bars below trade blue for amber as the observer's speed changes, while the purple total stays put.
""")
note("N24", "Two everyday pictures", r"""
a traffic light that switches from east–west to north–south traffic is an *unsteady* (local) change at a fixed
crossing; a roller coaster on a fixed track gives its riders *advective* acceleration although the track never changes.
""")
note("N22", "Linear and quadratic", r"""
The local term $\partial\mathbf u/\partial t$ is linear in $\mathbf u$; the advective term $(\mathbf u\cdot\nabla)\mathbf u$
is quadratic — doubling every velocity doubles the first and quadruples the second. This nonlinearity is the heart of
fluid mechanics (turbulence, waves steepening). When $\mathbf u$ is tiny it can be dropped (acoustics, Ch. 15); when
$\mathbf u=0$ only statics is left (Ch. 1 §1.7).
""")
nb.code("""
f1 = ch03.cylinder_velocity_field(1.0, 1.0, U_frame=0.0)       # the fluid-frame cylinder flow u(x, t) [m/s]
f2 = lambda x, t: 2*f1(x, t)                                   # every velocity doubled
Pp = np.array([0.0, 1.5])                                      # the point P [m]
a1, a2 = ch03.acceleration(f1, Pp, 0.0), ch03.acceleration(f2, Pp, 0.0)   # (a, local, advective) by stencils
print(round(a2.local[1]/a1.local[1], 6), round(a2.advective[1]/a1.advective[1], 6))   # 2.0 and 4.0
""", explain="The scaling test in two numbers: the local term doubles, the advective term quadruples.")
nb.worked_example("a travelling sine wave, two frames", r"""
In a frame moving with the wave the flow is steady: $u'=0.5\sin x'$ m/s (k = 1 m⁻¹). The lab sees the same pattern
moving at U = 2 m/s: $u=2+0.5\sin(x-2t)$. At the lab point $x-2t=\pi/4$:

1. Wave frame: local 0; advective $u'\,du'/dx'=0.5\sin\frac\pi4\times0.5\cos\frac\pi4=0.3536\times0.3536=0.125$ m/s².
2. Lab frame: local $\partial u/\partial t=-2\times0.5\cos\frac\pi4=-0.707$ m/s²; advective
   $u\,\partial u/\partial x=(2+0.354)\times0.354=0.832$ m/s².
3. Sum $-0.707+0.832=0.125$ m/s² — the same, as
   $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$
   *(Eq. 3.9)* says.
4. The cylinder at P = (0, 1.5) m (U = a = 1): body frame local 0, advective −0.856; fluid frame local −0.593, advective
   −0.263 (y-components, m/s²) — again −0.856 in total.
""")
nb.code("""
U = np.array([2.0])                                            # the wave moves at 2 m/s
up = lambda xp, tp: 0.5*np.sin(xp)                             # steady in the wave frame: u' = 0.5 sin x' [m/s]
u = lambda x, t: U + up(x - U*t, t)                            # the lab field u = U + u'(x − Ut)
print("lab:       ", ch03.acceleration(u, np.array([np.pi/4]), 0.0))               # (a, local, advective)
print("wave frame:", ch03.acceleration(ch03.galilean_transform(u, U), np.array([np.pi/4]), 0.0))
for Uf in (1.0, 0.5, 0.0):                                     # three observers: body, half-way, fluid [m/s]
    r = ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, Uf)  # at P = (0, 1.5) m, U = a = 1
    print(f"U_frame = {Uf}: local {np.round(r['local'], 3)}, advective {np.round(r['advective'], 3)}, total {np.round(r['total'], 3)}")
""", explain="""
1. A field that is steady in the moving frame.
2. `acceleration` returns (a, local, advective) by second-order stencils in the lab: $(0.125,\,-0.707,\,0.832)$ m/s².
3. `galilean_transform` builds the primed field; the same call gives the wave-frame split $(0.125,\,0,\,0.125)$.
4. `frame_acceleration_terms` does it for the cylinder for three observers — the totals agree to 10⁻⁸.
""")
nb.figure(r"""
Ufs = np.array([0.0, 0.25, 0.5, 0.75, 1.0])                    # observer speeds: fluid frame → body frame [m/s]
R_ = [ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, Uf) for Uf in Ufs]   # the (3.9) terms at P for each
loc = np.array([r["local"][1] for r in R_]); adv = np.array([r["advective"][1] for r in R_])   # y-components [m/s²]
fig, ax = plt.subplots(figsize=(7, 3.8))                       # one panel of stacked bars
ax.bar(Ufs, loc, width=0.18, color=COLORS["blue"], label=r"local $\partial v/\partial t$")               # blue: local
ax.bar(Ufs, adv, width=0.18, bottom=loc, color=COLORS["amber"], label=r"advective $(\mathbf{u}\cdot\nabla)v$")   # amber on top
ax.plot(Ufs, loc + adv, "D", color=COLORS["accent"], ms=9, label="total $Dv/Dt$")                          # purple: the sum
ax.axhline(-0.856, ls="--", color=COLORS["accent"], lw=1)      # the common total −0.856 m/s²
ax.set_xlabel("observer speed $U_{frame}$ [m/s]  (0 = still lake, 1 = cylinder)"); ax.set_ylabel("y-acceleration at P [m/s²]")
ax.set_title("The split moves with the observer; the total does not"); ax.legend(fontsize=8, loc="lower left")   # message title
savefig(fig, "ch03", "galilean_term_bars"); plt.show()         # save to outputs/ch03 and draw
""", see="Blue bars shrinking to zero as the observer approaches the body frame, amber bars growing, purple diamonds on one horizontal line.",
    read="Reading left to right is moving the observer from the lake to the cylinder; the acceleration of the particle at P never changes — only how it is booked.",
    change="…P were far from the body (P = (0, 5) m): all bars would shrink toward zero (the flow there is almost uniform).")
nb.md(r"""
📎 *Gloss — a stream function as a drawing tool.* For a plane flow with $\nabla\cdot\mathbf u=0$ there is a function
$\psi(x,y)$ (Ch. 6) with $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$. Then
$\mathbf u\cdot\nabla\psi=\frac{\partial\psi}{\partial y}\frac{\partial\psi}{\partial x}-\frac{\partial\psi}{\partial x}\frac{\partial\psi}{\partial y}=0$:
ψ does not change along a streamline, so its contour lines *are* the streamlines. Adding $c\,y$ to ψ adds $u=c$, $v=0$
everywhere — a uniform flow along x. For the cylinder in the body frame $\psi=Uy(1-a^2/r^2)$
(`ch03.cylinder_streamfunction`); the observer moving at $U_{frame}$ sees the far fluid at $U_{frame}$ instead of $U$,
i.e. $\psi+(U_{frame}-U)\,y$. We use it below only to draw many streamlines cheaply; the cell checks both claims.
""")
nb.code("""
xs_, ys_, U_, a_, c_ = sp.symbols('x y U a c', positive=True)       # position, stream speed, radius, a constant
psi = U_*ys_*(1 - a_**2/(xs_**2 + ys_**2))                           # body-frame ψ of the cylinder
print(sp.limit(sp.diff(psi, ys_), xs_, sp.oo), sp.simplify(sp.diff(psi + c_*ys_, ys_) - sp.diff(psi, ys_)))   # U far away; adding c·y adds u = c
""")
nb.plotly(r"""
x = np.linspace(-3, 3, 151); y = np.linspace(-2.5, 2.5, 126); X, Y = np.meshgrid(x, y)   # drawing grid [m]
psi_body = ch03.cylinder_streamfunction(X, Y, 1.0, 1.0)        # body-frame ψ [m²/s] (NaN inside the cylinder)
levels = np.linspace(-1.5, 1.5, 13)                            # which streamlines to draw
th = np.linspace(0, 2*np.pi, 80)

def frame(Uf):                                                 # streamlines + body + P for observer speed Uf [m/s]
    psi = psi_body + (Uf - 1.0)*Y                              # ψ seen by this observer
    cs = plt.contour(X, Y, psi, levels=levels); segs = cs.allsegs; plt.close()   # contour lines as point lists
    xs, ys = [], []                                            # all contour pieces in one NaN-separated line
    for level in segs:                                         # one list of pieces per level
        for sgm in level:                                      # each piece is an (n, 2) array of points
            xs += list(sgm[:, 0]) + [np.nan]; ys += list(sgm[:, 1]) + [np.nan]
    return {"streamlines": (xs, ys), "cylinder": (np.cos(th), np.sin(th)), "point P": ([0.0], [1.5])}   # name → (x, y)

vals = np.linspace(0, 1, 21 if not FAST else 11)               # observer speeds [m/s]
fig = slider_figure(frame, "U_frame", vals, unit="m/s", xlabel="x [m]", ylabel="y [m]", xrange=[-3, 3],
                    yrange=[-2.5, 2.5], modes={"point P": "markers"}, height=480,
                    title="The observer's speed reshapes the streamlines (P's terms are in its legend entry)")
Ps = [tr for tr in fig.data if tr.name == "point P"]           # one P trace per slider step, in slider order
for tr, Uf in zip(Ps, vals):                                   # write that observer's (3.9) terms into P's name
    r = ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, Uf)  # local, advective, total at P [m/s²]
    tr.name = f"P: local {r['local'][1]:.3f} + advective {r['advective'][1]:.3f} = {r['total'][1]:.3f} m/s²"
    tr.marker.size = 11; tr.marker.color = COLORS["accent"]    # a big purple dot
recolor(fig, {"streamlines": COLORS["teal"], "cylinder": COLORS["ink"]})   # house colours
fig.update_yaxes(scaleanchor="x", scaleratio=1)                # equal scales, so the cylinder is round
fig.show()                                                     # draw it
""", explain="""
1. For each observer speed the stream function is shifted by $(U_{frame}-U)\\,y$ and its contour lines (matplotlib's
   `contour`, Ch. 2 P78, used here only to get the curves) are passed to plotly.
2. The point P carries that observer's local, advective and total accelerations in its legend entry.
""")
see_read_change(
    "Loops leaving the front of the body and re-entering its back at U_frame = 0; streamlines flowing round the body at U_frame = 1 m/s.",
    "The loops of the fluid frame open into the flow-past-a-body picture as $U_{frame}\\to U$; the numbers in P's label trade places while their sum stays −0.856 m/s².",
    "…U_frame went beyond U: the far fluid would stream past faster than U (still towards +x) and the cylinder itself would drift downstream (+x); the total at P would still be −0.856 m/s².")
nb.explainer("galilean_frames_cylinder", "Steady or not — does the acceleration care?", r"""
**Why interactive:** the point is that one number stays fixed while everything else changes; you only believe it when you drag the observer yourself and watch the bars trade places around a purple total that never moves.
One slider moves the observer from the still fluid to the cylinder. The streamline picture morphs from "unsteady, loops
on the body" to "steady, around the body", while the acceleration arrow of a tagged particle stays exactly the same and
its two bars trade places —
$\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$
*(Eq. 3.9)* as an experiment.
""", tries=[
    "Drag the observer speed from 0 to U and watch the blue bar hand its share to the amber bar.",
    "Click any point: the inspector gives both frames' terms with numbers.",
    "Set the preset 'body frame': the status says steady — yet the purple arrow is not zero.",
    "Step through D06 in the Derivation tab; at the cancellation step the two U·∇′ bars light up.",
])
nb.md(r"""
> ⚠️ **Common confusion:** "Galilean invariance means any moving frame is as good as any other." Only frames moving at
> *constant* velocity without rotating. In a rotating frame (the Earth!) the acceleration picks up Coriolis and
> centrifugal terms (Ch. 4 §4.7) and even the vorticity changes by 2Ω (N28 in C10) — the heart of geophysical fluid
> dynamics (Ch. 13).
""")
whatif(r"""
…the observer accelerated (U depends on t)? Step 3 of D06 gains $-d\mathbf U/dt$, which does not cancel: an
accelerating observer sees a fictitious force. Next, §3.4 zooms into a small neighbourhood of a point and asks what the
velocity *differences* do to a fluid element.
""")

# =====================================================================================================================
# A.4 §3.4 Strain and Rotation Rates — R01–R03, C06–C09, R04–R07, C10–C12
# =====================================================================================================================
nb.section("3.4", "Strain and Rotation Rates", intro="""
**What is this section about?** Zoom into a tiny blob of fluid. Its neighbours move at slightly different velocities;
those differences stretch it, shear it, swell it and spin it. Each of these motions is one piece of the
velocity-gradient tensor you met in Ch. 2 — here you *measure* them on a moving element.
""")
nb.recap("R01", "The split of the velocity gradient", r"""
Any velocity gradient splits uniquely into a symmetric and an antisymmetric part:
$\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}$ *(Eq. 3.11)* (Ch. 2 proved the split is unique).
`ch03.strain_rate_tensor(G) + 0.5*ch03.rotation_tensor(G)` returns G exactly.
""", where="Ch. 2 §2.10")
nb.recap("R02", "The strain-rate tensor", r"""
The symmetric part $S_{ij}=\tfrac12\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big)$
*(Eq. 3.12)* — `ch03.strain_rate_tensor(G)`. §3.4 shows what each of its entries *measures*.
""", where="Ch. 2 §2.10")
nb.recap("R03", "The rotation tensor — the book's, without ½", r"""
$R_{ij}=\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}$ *(Eq. 3.13)*, `ch03.rotation_tensor(G)` =
G − Gᵀ. ⚠️ The book's R is **twice** the antisymmetric part A of Ch. 2 (R = 2A); that is why a ½ stands in front of R
in $\frac{\partial u_i}{\partial x_j}=S_{ij}+\tfrac12R_{ij}$ *(Eq. 3.11)*.
""", where="Ch. 2 §2.10")

core("C06", "Relative velocity of a neighbour: du = G·dx (3.10)", r"""
Two corks float 1 cm apart. How fast does one move relative to the other — and what does that single matrix of nine
numbers already tell us about the blob of water between them?
""", eqs=("3.10",))
nb.md(r"""
#### The problem in plain words

A drop of ink in a stream is stretched into a thin filament, a cloud is sheared into a streak, air in a rising thermal
swells. All of these are about *differences* of velocity between neighbouring points. Chapter 4's stress law needs
exactly these differences (the rates of deformation), not the velocity itself.
""")
nb.md(r"""
#### The idea

```
O at x, velocity u          P at x + dx, velocity u + du          (Fig. 3.9)
zoom in far enough and every smooth field looks linear:   du = G · dx,   G_ij = ∂u_i/∂x_j   (nine numbers)
G = S + ½R :   S deforms the blob  (stretch, shear, swell)   ½R spins it
```
""")
P("P98", "multivariable first-order Taylor expansion", r"""
Ch. 1's $f(z+dz)\approx f(z)+f'(z)\,dz$ (P26) in several variables:
$u_i(\mathbf x+d\mathbf x)\approx u_i(\mathbf x)+\frac{\partial u_i}{\partial x_j}dx_j$ — one partial derivative per
direction, summed over j. The neglected terms shrink like $|d\mathbf x|^2$, so after dividing by $|d\mathbf x|$ they
vanish as the neighbour comes closer (Ch. 2 P68, orders of smallness).
""", code="""
f = lambda x, y: np.sin(x)*np.exp(y)                 # a smooth function of two variables
x0, y0, dx, dy = 0.3, 0.1, 0.01, -0.02               # a point and a small step
lin = f(x0, y0) + np.cos(x0)*np.exp(y0)*dx + np.sin(x0)*np.exp(y0)*dy   # value + (∂f/∂x)dx + (∂f/∂y)dy
print(f(x0 + dx, y0 + dy) - lin)                     # ≈ −1e-4: the leftover is second order in the step
""")
D("D07", ref="3.10")
note("N25", "Why rates of deformation?", r"""
A solid resists being *deformed* (strain); a fluid resists being deformed *quickly* (strain rate). Ch. 4 §4.5 relates
the stress to S — the Newtonian law that generalises Ch. 1's $\tau=\mu\,du/dy$.
""")
note("N26", "Names and roles", r"""
S (deformation) is linked to the stress in a moving fluid; R (rotation) is not — a rigidly spinning fluid feels no
viscous stress. C11 shows why only S can enter.
""")
nb.worked_example("G = [[1, 2], [0, −1]] s⁻¹ and a neighbour 2 cm away", r"""
1. $d\mathbf x=(0.01,0.02)$ m.
2. $du_1=G_{11}dx_1+G_{12}dx_2=1\times0.01+2\times0.02=0.05$ m/s.
3. $du_2=G_{21}dx_1+G_{22}dx_2=0+(-1)\times0.02=-0.02$ m/s.
4. Split: $\mathbf S=[[1,1],[1,-1]]$, $\tfrac12\mathbf R=[[0,1],[-1,0]]$ s⁻¹, and $\mathbf S+\tfrac12\mathbf R=[[1,2],[0,-1]]=\mathbf G$ ✓
   (R01).
5. The field $\mathbf u=(x_1+2x_2+x_1^2,\ -x_2+x_1x_2)$ has exactly this G at the origin; its true du at that neighbour is
   (0.0501, −0.0198) m/s — the difference (1, 2) × 10⁻⁴ m/s is the second-order remainder.
""")
nb.code("""
ufield = lambda x, t: np.array([x[0] + 2*x[1] + x[0]**2, -x[1] + x[0]*x[1]])   # a nonlinear plane field [m/s]
G = ch03.velocity_gradient_at(ufield, np.zeros(2), 0.0)   # G_ij = ∂u_i/∂x_j at the origin by central stencils [1/s]
dx = np.array([0.01, 0.02])                                 # the neighbour's offset [m]
print(G)
print(ch03.relative_velocity(G, dx), ufield(dx, 0) - ufield(np.zeros(2), 0))   # (3.10) vs the exact difference [m/s]
S, R = ch03.strain_rate_tensor(G), ch03.rotation_tensor(G)  # (3.12) and the book's (3.13)
assert np.allclose(S + 0.5*R, G)                            # (3.11): G = S + ½R
e = dx/np.linalg.norm(dx); hs = 0.02*0.5**np.arange(6)      # one direction, six shrinking distances [m]
err = [np.linalg.norm(ufield(h*e, 0) - ufield(np.zeros(2), 0) - G @ (h*e)) for h in hs]   # |exact du − G·dx|
print(f"observed order of the remainder: {observed_order(hs, err):.2f}")   # 2: the neglected terms are O(|dx|²)
""", explain="""
1. A nonlinear field and its gradient at the origin, $G=[[1,2],[0,-1]]$ s⁻¹, by second-order stencils.
2. $d\\mathbf u=\\mathbf G\\cdot d\\mathbf x$, i.e. $du_i=(\\partial u_i/\\partial x_j)\\,dx_j$ *(3.10)*, gives (0.05, −0.02) m/s;
   the exact difference is (0.0501, −0.0198) m/s.
3. The R01 split holds exactly.
4. The remainder falls like $|d\\mathbf x|^2$: slope 2 on log–log axes (Ch. 1 P13).
""")
nb.figure(r"""
from scripts.ch03_drawings import ring_arrows                # draws du = G·dx on a ring of neighbours (no physics inside)
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2))          # (a) the ring, (b) the remainder
ring_arrows(a, G, radius=0.1, parts=("total",), n=24, scale=0.25)   # 24 neighbours 10 cm from O; arrows × 0.25
a.plot(0, 0, "o", color=COLORS["ink"]); a.set_aspect("equal"); a.set_xlim(-0.2, 0.2); a.set_ylim(-0.2, 0.2)   # O at the centre
a.set_xlabel("$dx_1$ [m]"); a.set_ylabel("$dx_2$ [m]"); a.set_title("(a) relative velocity du = G·dx of each neighbour")
b.loglog(hs, err, "o-", color=COLORS["accent"], label="|exact du − G·dx|")   # measured remainder
b.loglog(hs, err[0]*(hs/hs[0])**2, "--", color=COLORS["muted"], label="slope 2")   # a guide ∝ |dx|²
b.set_xlabel("distance |dx| [m]"); b.set_ylabel("remainder [m/s]"); b.legend(fontsize=8)   # log–log axes (Ch. 1 P13)
b.set_title("(b) the linear picture becomes exact as the ring shrinks")
savefig(fig, "ch03", "relative_velocity_ring"); plt.show()   # save to outputs/ch03 and draw
""", see="(a) Purple arrows on a ring, pointing out along one direction and in along another, plus a swirl; (b) a straight line of slope 2.",
    read="""(a) is everything the element does in the next instant: each neighbour's velocity relative to O, zero at O and
growing with distance. (b) says the linear picture is exact in the limit — the neglected terms fall like the square of
the distance.""",
    change="…G were antisymmetric (a pure rotation, e.g. [[0, −1], [1, 0]]): every arrow on the ring would be tangential — no stretching at all.")
whatif(r"""
…you watched only two neighbours on the $x_1$-axis? Their separation changes at $\partial u_1/\partial x_1$ per unit
length — the diagonal of S. That is C07.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C07", "Linear strain rate: how fast a material line stretches", r"""
A short thread of dye lies along the flow. At what rate, per unit of its length, is it being stretched — and in which
direction would it stretch fastest?
""")
nb.md(r"""
#### The problem in plain words

Stretching is how ocean eddies draw tracers into long filaments and how vortex tubes intensify (Ch. 5). A rubber band
held at both ends and pulled: the stretching rate *per unit length* is what matters, not the length itself.
""")
nb.md(r"""
#### The idea

```
A ●────────● B     ends move at u₁ and u₁ + (∂u₁/∂x₁)δx₁   (Fig. 3.10)
B outruns A by (∂u₁/∂x₁)δx₁ per second   →   stretch per length per time = ∂u₁/∂x₁ = S₁₁
any direction n:   (1/ℓ) Dℓ/Dt = n·S·n
```
""")
P("P99", "material line element", r"""
A *material* line element $\delta\mathbf x$ joins two nearby fluid particles and is carried with them — it stretches
and turns as they move. Its rate of change following the particles is the difference of their velocities:
$\frac{D(\delta\mathbf x)}{Dt}=\delta\mathbf u=\mathbf G\cdot\delta\mathbf x$ (C06). Every rate of §3.4 is measured on
such elements. For a linear flow the element after time t is $e^{\mathbf Gt}\cdot\delta\mathbf x$ (`linear_flow_map`,
the matrix exponential of Ch. 2 P79).
""", code="""
G2 = np.array([[2.0, 0.0], [0.0, -2.0]])             # u = (2x, −2y) [1/s × m]
dx0 = np.array([0.01, 0.0])                          # a 1 cm element along x [m]
dx1 = ch03.linear_flow_map(G2, 0.01) @ dx0           # carried for 0.01 s: e^{G t}·dx0
print(dx1*100)                                       # in cm: [1.0202 0.]: it grew by 2 %
""")
D("D08")
nb.worked_example("u = (2x, −2y) s⁻¹, a 1 cm thread along x for 0.01 s", r"""
1. $S_{11}=\partial u_1/\partial x_1=2$ s⁻¹.
2. First order: new length $\approx1\times(1+2\times0.01)=1.02$ cm.
3. Exact (the ends move as $e^{2t}$): $e^{0.02}=1.0202$ cm.
4. A thread along y: $S_{22}=-2$ s⁻¹ — it shrinks to $e^{-0.02}=0.9802$ cm.
5. A thread at 45°: $\mathbf n\cdot\mathbf S\cdot\mathbf n=\tfrac12(2-2)=0$ — its length does not change at first order.
""")
nb.code("""
G = np.array([[1.0, 2.0], [0.0, -1.0]])              # C06's velocity gradient [1/s]
for deg in (0, 22.5, 45, 90):                        # probe directions [degrees]
    n = [np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))]   # unit vector n (angles in radians inside)
    print(f"{deg:5.1f}°: n·S·n = {ch03.linear_strain_rate(G, n):+.3f} 1/s")   # stretching rate per unit length
""", explain="""
1. The probe direction n (`np.deg2rad` converts degrees to radians, Ch. 2 P66).
2. `linear_strain_rate` returns n·S·n (the antisymmetric part drops out of any n·G·n, D08 step 8):
   $\\cos2\\theta+\\sin2\\theta$ for this G.
3. The largest stretching, √2 = 1.414 s⁻¹ at 22.5°, is along a principal axis (C12).
""")
nb.check_agree("""
dt = 1e-6                                            # a very short time [s]
for deg in (0, 22.5, 45, 90):
    n = np.array([np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))])
    dx0 = 0.01*n                                     # a 1 cm thread along n [m]
    M = ch03.linear_flow_map(G, dt)                  # carry it for dt
    mine = (np.linalg.norm(M @ dx0) - np.linalg.norm(dx0))/(np.linalg.norm(dx0)*dt)   # (ΔL/L)/Δt, measured
    assert np.allclose(mine, ch03.linear_strain_rate(G, n), rtol=1e-5, atol=1e-5)   # = n·S·n
    assert np.allclose(ch03.measured_strain_rates(G, n, dt)["stretch"], mine, rtol=1e-5, atol=1e-5)   # library's ruler
print("tracked threads stretch at n·S·n in every direction")
""")
nb.figure(r"""
th = np.linspace(0, np.pi, 72, endpoint=False)        # 72 directions from 0° to 180° [rad]
meas = [ch03.measured_strain_rates(G, [np.cos(a), np.sin(a)], 1e-6)["stretch"] for a in th]   # tracked threads
form = [ch03.linear_strain_rate(G, [np.cos(a), np.sin(a)]) for a in th]                       # n·S·n
fig, ax = plt.subplots(figsize=(7.5, 3.8))                    # one panel
ax.plot(np.degrees(th), form, color=COLORS["teal"], lw=2.5, label="formula n·S·n")            # the tensor formula
ax.plot(np.degrees(th), meas, "o", color=COLORS["accent"], ms=3, label="measured (1/ℓ)Δℓ/Δt")   # the ruler on moving threads
ax.axhline(0, color=COLORS["muted"], lw=0.8)                  # no stretching
ax.axvline(22.5, ls="--", color=COLORS["blue"], label="fastest stretching (22.5°)")          # principal axis (C12)
ax.axvline(112.5, ls="--", color=COLORS["rose"], label="fastest squeezing (112.5°)")         # the other one, 90° away
ax.set_xlabel("direction θ of the thread [°]"); ax.set_ylabel("stretching rate [1/s]"); ax.legend(fontsize=8)
ax.set_title("Each direction has its own stretching rate; the extremes are 90° apart")
savefig(fig, "ch03", "stretching_rose"); plt.show()           # save to outputs/ch03 and draw
""", see="Purple dots on a smooth teal wave between −√2 and +√2 s⁻¹, with dashed blue and rose lines at the peak and the trough.",
    read="""The measurement on moving threads and the tensor formula are the same curve. Threads near 22.5° stretch fastest,
threads near 112.5° are squeezed fastest; the zero crossings (at 67.5° and 157.5°) are directions whose length does not
change at first order.""",
    change="…G were antisymmetric (solid-body rotation): the curve would be flat at zero — nothing stretches.")
whatif(r"""
…two threads start perpendicular? Besides stretching they turn toward (or away from) each other; the rate at which the
right angle closes is the off-diagonal of S (C08).
""")

# ---------------------------------------------------------------------------------------------------------------------
nb.recap("R04", "The rotation tensor is a vector in disguise", r"""
R is antisymmetric, $R_{ij}=-R_{ji}$: zero diagonal and three independent entries, which pair up with the three
components of the vorticity $\boldsymbol\omega=\nabla\times\mathbf u$.
""", where="Ch. 2 §2.10")
nb.recap("R05", "R ↔ ω", r"""
$R_{ij}=-\varepsilon_{ijk}(\nabla\times\mathbf u)_k=-\varepsilon_{ijk}\omega_k=\begin{bmatrix}0&-\omega_3&\omega_2\\\omega_3&0&-\omega_1\\-\omega_2&\omega_1&0\end{bmatrix}$
*(Eq. 3.15, the book's (2.26)–(2.27) with R = 2A)* — `ch03.antisymmetric_from_vector(omega)` builds the matrix;
`ch03.vorticity_from_gradient(G)` returns ω = the vector of R.
""", where="Ch. 2 §2.10")
core("C08", "Shear strain rate: how fast a right angle closes", r"""
Two dye threads cross at a right angle. How fast does that angle close, and why does the book call half of that rate
$S_{12}$?
""")
nb.md(r"""
#### The problem in plain words

Push the top of a deck of cards sideways: the corners of each card stop being right angles. Fluids flowing past a wall
are sheared this way all the time (Couette and boundary-layer flows, Ch. 8–9), and the shear *rate* is what the viscous
shear stress depends on.
""")
nb.md(r"""
#### The idea

```
     C ┐ dα                 vertical side tilts clockwise by dα        ∝ ∂u₁/∂x₂
       │ ╲                  horizontal side tilts counterclockwise by dβ ∝ ∂u₂/∂x₁
     B └───── dβ            closing rate α + β  →  S₁₂ = ½ D(α+β)/Dt        (Fig. 3.11)
```
""")
P("P100", "small-angle approximation", r"""
For an angle ε measured in radians, $\tan\varepsilon\approx\varepsilon$ and $\cos\varepsilon\approx1$ with errors of
order ε³ and ε²: a tiny tilt is its own tangent. Every angle below is proportional to dt, so the neglected pieces
vanish after dividing by dt.
""", code="""
eps = np.array([0.1, 0.01, 0.001])                   # three small angles [rad]
print([f"{v:.2e}" for v in np.tan(eps) - eps])     # ['3.35e-04', '3.33e-07', '3.33e-10']: the error falls like ε³
""")
D("D09")
P("P101", "rigid-body velocity Ω × x", r"""
A rigid body turning at angular velocity $\boldsymbol\Omega$ about an axis through the origin moves each of its points
with $\mathbf v=\boldsymbol\Omega\times\mathbf x$ — perpendicular to both the axis and the arm, of size Ω times the
distance from the axis. Adding a translation $\mathbf U$ gives every rigid motion: $\mathbf u=\mathbf U+\boldsymbol\Omega\times\mathbf x$.
""", code="""
Omega = np.array([0.0, 0.0, 1.0])                    # 1 rad/s about z
print(np.cross(Omega, [2.0, 0.0, 0.0]))              # [0. 2. 0.]: 2 m from the axis → 2 m/s, tangential
""")
D("D10")
note("N27", "Rigid motion does not deform", r"""
For $\mathbf u=\mathbf U+\boldsymbol\Omega\times\mathbf x$ with $\mathbf U$ and $\boldsymbol\Omega$ uniform, S = 0 and R
carries all of G (ω = 2Ω). Hence adding a rigid motion — watching from a translating or rotating frame — changes G only
by an antisymmetric part: **S is the same for every such observer** (even if U varies in time); ω is not (N28).
""")
nb.worked_example("simple shear u = (γy, 0), γ = 1 s⁻¹, for dt = 0.01 s", r"""
1. The top of a 1 cm vertical side moves $\gamma\,\delta x_2\,dt=1\times0.01\times0.01=10^{-4}$ m further than its foot:
   $\tan d\alpha=10^{-4}/10^{-2}=0.01$, so $d\alpha\approx0.01$ rad.
2. The horizontal side does not tilt: dβ = 0 ($u_2=0$ everywhere).
3. $S_{12}=\tfrac12(0.01+0)/0.01=0.5$ s⁻¹ $=\gamma/2$.
4. So $\gamma=2S_{12}$. ⚠️ Ch. 2's Ex. 2.4 called $S_{12}$ itself Γ; here γ is twice that.
""")
nb.code("""
G = ch03.velocity_gradient_preset("simple_shear", Gamma=1.0)   # G = [[0, γ], [0, 0]], γ = 1 1/s: Gamma= is the preset's rate
                                                               # (γ here, ω0 for solid body) — NOT the circulation Γ
print(ch03.shear_strain_rate(G, [1, 0], [0, 1]))               # S12 = n1·S·n2 = 0.5 1/s
print(ch03.measured_strain_rates(G, [1, 0], 1e-6)["closing"]/2)   # half the measured closing rate of two threads
U, Om = np.array([0.3, -0.1, 0.2]), np.array([0.5, -1.0, 2.0])   # a translation [m/s] and a rotation [rad/s]
Gr = ch03.velocity_gradient_at(lambda x, t: ch03.rigid_body_velocity(U, Om, x), np.array([0.4, 0.1, -0.7]), 0.0)
print(np.round(ch03.strain_rate_tensor(Gr), 10) + 0.0)          # a 3 × 3 zero matrix: rigid motion does not deform
print(ch03.vorticity_from_gradient(Gr))                        # [ 1. -2.  4.] = 2Ω
""", explain="""
1. $S_{12}$ from the formula $\\mathbf n_1\\cdot\\mathbf S\\cdot\\mathbf n_2$.
2. Half the closing rate of two tracked perpendicular threads — the same 0.5 s⁻¹ (D09).
3. A rigid motion $\\mathbf U+\\boldsymbol\\Omega\\times\\mathbf x$ (D10): its gradient has no symmetric part, and its vorticity
   is twice its angular velocity.
""")
nb.figure(r"""
ts = np.linspace(0, 1, 101)                                    # time [s]
flows = {"simple shear γ = 1 s⁻¹": (ch03.velocity_gradient_preset("simple_shear", 1.0), COLORS["accent"]),
         "pure strain along ±45°": (ch03.velocity_gradient_preset("irrotational_strain", 0.5), COLORS["teal"]),
         "solid-body rotation": (ch03.velocity_gradient_preset("solid_body_rotation", 1.0), COLORS["orange"])}
fig, ax = plt.subplots(figsize=(7.5, 3.8))                      # one panel
for name, (Gf, col) in flows.items():                           # three linear flows
    ang = [ch03.material_line_angle(Gf, t, np.pi/2) - ch03.material_line_angle(Gf, t, 0.0) for t in ts]   # angle between the threads
    s12 = ch03.strain_rate_tensor(Gf)[0, 1]                     # S12 of this flow [1/s]
    ax.plot(ts, np.degrees(ang), color=col, lw=2.3, label=f"{name} ($S_{{12}}$ = {s12:.1f} 1/s)")   # measured angle [°]
    ax.plot(ts[:40], 90 - np.degrees(2*s12*ts[:40]), ":", color=col, lw=1.2)   # initial tangent: slope −2S12
ax.set_xlabel("$t$ [s]"); ax.set_ylabel("angle between the threads [°]"); ax.legend(fontsize=8)   # axes with units
ax.set_title("A right angle closes at 2S₁₂ at first")
savefig(fig, "ch03", "closing_angle"); plt.show()               # save to outputs/ch03 and draw
""", see="The orange curve stays at 90°; the purple and teal curves fall, each starting along its dotted tangent.",
    read="""The initial slope is $-2S_{12}$ in rad/s: shear and pure strain close the angle at the same initial rate (both
have $S_{12}=0.5$ s⁻¹); rotation turns both threads together and closes nothing (S = 0). Later the curves part: the
shear keeps turning its threads, the pure strain pulls both toward its stretching axis.""",
    change="""…the threads were turned by 45°: in the pure-strain flow their angle would stay 90° (they lie on the principal
axes) — the off-diagonal of S depends on the axes (C12).""")
whatif(r"""
…three edges of a small box stretch at once? Then its volume changes — at the sum of the three linear rates (C09).
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C09", "Volumetric strain rate: divergence is how fast a blob swells (3.14)", r"""
A small parcel of air rises and expands. At what rate does its volume grow per unit volume, and how is that read off
the velocity field?
""", eqs=("3.14",))
nb.md(r"""
#### The problem in plain words

Rising air expands, sinking air is compressed; in the ocean water hardly changes volume at all. The fractional growth
rate of a parcel's volume is the kinematic half of mass conservation (Ch. 4: $D\rho/Dt=-\rho\nabla\cdot\mathbf u$). We
find it from the three stretching rates of C07.
""")
nb.md(r"""
#### The idea

```
δV = δx₁ δx₂ δx₃   each edge stretches at S₁₁, S₂₂, S₃₃   →   (1/δV) D(δV)/Dt = S₁₁ + S₂₂ + S₃₃ = ∇·u   (3.14)
shear only tilts the faces (second order);  the trace does not depend on the orientation of the box
```

In symbols: $\frac{1}{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_1}{\partial x_1}+\frac{\partial u_2}{\partial x_2}+\frac{\partial u_3}{\partial x_3}=\frac{\partial u_i}{\partial x_i}$ *(3.14)*.
""")
D("D11", ref="3.14")
nb.worked_example("u = (x, y, z) s⁻¹ and a 1 cm³ box for 0.01 s", r"""
1. $\nabla\cdot\mathbf u=1+1+1=3$ s⁻¹.
2. First order: $\delta V\approx1\times(1+3\times0.01)=1.03$ cm³.
3. Exact for this linear flow (Jacobi's formula, D11 step 7): $e^{0.03}=1.03045$ cm³.
4. Simple shear $\mathbf u=(\gamma y,0)$: $\nabla\cdot\mathbf u=0$ — the sheared box keeps its volume although it changes
   shape.
""")
nb.code("""
G3 = np.eye(3)                                            # u = (x, y, z): uniform expansion [1/s]
print(ch03.volumetric_strain_rate(G3), ch03.material_volume_ratio(G3, 0.01))   # tr G = 3 1/s; det e^{Gt} after 0.01 s
rng = np.random.default_rng(3)                            # reproducible random numbers (Ch. 1 P10)
Grand = rng.normal(size=(3, 3))                           # a random 3-D velocity gradient [1/s]
print(ch03.material_volume_ratio(Grand, 0.7), np.exp(0.7*np.trace(Grand)))   # Jacobi: det e^{Gt} = e^{t tr G}
C = ch03.rotation_matrix_2d(0.4)                          # axes turned by 0.4 rad
G2 = np.array([[1.0, 2.0], [0.0, -1.0]])                  # C06's G
print(np.trace(ch03.transform_tensor(G2, C)), np.trace(G2))   # the trace does not care about the axes
""", explain="""
1. The trace of G is the volumetric rate $\\frac{1}{\\delta V}\\frac{D(\\delta V)}{Dt}=\\frac{\\partial u_i}{\\partial x_i}$ *(3.14)*;
   after 0.01 s the exact volume factor is 1.030455.
2. The exact finite-time factor $\\det e^{\\mathbf Gt}$ equals $e^{t\\,\\mathrm{tr}\\,\\mathbf G}$ (Jacobi) for any G.
3. The trace does not change when the axes turn (Ch. 2 §2.5).
""")
nb.check_agree("""
for t in (0.01, 0.5):                                     # two times [s]
    M = ch03.linear_flow_map(G3, t)                       # carries every tracer: x(t) = e^{Gt} x(0)
    edges = M @ (0.01*np.eye(3))                          # the three 1 cm edge vectors of the cube after time t [m]
    vol = abs(np.linalg.det(edges))                       # volume of the carried cube [m³] (Ch. 1 P56)
    assert np.allclose(vol/1e-6, ch03.material_volume_ratio(G3, t))   # δV(t)/δV(0) = det e^{Gt}
    print(f"t = {t}: tracked cube {vol/1e-6:.6f} × its initial volume")
""")
nb.figure(r"""
ts = np.linspace(0, 1, 11)                                # time [s]
fig, ax = plt.subplots(figsize=(7, 3.8))                  # one panel
tracked = [abs(np.linalg.det(ch03.linear_flow_map(G3, t))) for t in ts]   # tracked cube, expansion G = I
ax.plot(ts, tracked, "o", color=COLORS["blue"], label="tracked cube, $\\mathbf{G}=\\mathbf{I}$")   # dots
tt = np.linspace(0, 1, 101)                               # a fine time axis [s]
ax.plot(tt, np.exp(3*tt), color=COLORS["blue"], lw=2, label="$e^{3t}$ (exact)")                  # Jacobi: e^{t tr G}
ax.plot(tt, 1 + 3*tt, "--", color=COLORS["muted"], label="$1+3t$ (first order, (3.14))")         # the tangent
Gsh = ch03.velocity_gradient_preset("simple_shear", 1.0)  # a shear: tr G = 0
ax.plot(tt, [ch03.material_volume_ratio(Gsh, t) for t in tt], color=COLORS["teal"], lw=2, label="simple shear: 1")
ax.set_xlabel("$t$ [s]"); ax.set_ylabel(r"$\delta V(t)/\delta V(0)$ [–]"); ax.legend(fontsize=8)   # dimensionless ratio
ax.set_title("Divergence sets the initial swelling rate; shear changes shape, not volume")
savefig(fig, "ch03", "volume_ratio"); plt.show()          # save to outputs/ch03 and draw
""", see="An exponential (blue) that leaves its dashed tangent line; a flat teal line for shear.",
    read="""$\\frac{1}{\\delta V}\\frac{D(\\delta V)}{Dt}=\\frac{\\partial u_i}{\\partial x_i}$ *(3.14)* gives the initial slope,
3 s⁻¹ per unit volume; over finite time the growth compounds ($e^{3t}$). Shear alone changes shape, never volume.""",
    change="…G = −I (compression): the volume would decay as $e^{-3t}$ — a sinking, compressed parcel.")
nb.explainer("fluid_element_deformation", "What does each number in S measure?", r"""
**Why interactive:** S has several entries and each measures a different motion; changing G one entry at a time and seeing which measured rate responds is the fastest way to learn what each number means.
Drag the four entries of G and watch a square element and a ring of tracers deform. The measured rates — stretching
per length along a probe direction, the closing rate of a right angle, the area growth — land on the formula rates
n·S·n, $2S_{12}$ and tr G, and the rose curve shows the stretching rate in every direction at once.
""", tries=[
    "Preset 'solid-body rotation': every measured rate is zero — S = 0 (N27).",
    "Preset 'simple shear': rotate the probe to 45° and find the fastest stretching.",
    "Preset 'expansion': the area ratio follows the ghost e^{t tr G}.",
    "Click a tracer: the inspector gives its du = G·dx arithmetic (D07).",
])
whatif(r"""
…you looked not at how the element deforms but at how it turns? The two perpendicular threads of C08 turn by −α and
+β; their *average* turning is the spin, and it is half the vorticity (C10).
""")

nb.recap("R06", "The components of the vorticity", r"""
$\omega_1=\frac{\partial u_3}{\partial x_2}-\frac{\partial u_2}{\partial x_3},\ \omega_2=\frac{\partial u_1}{\partial x_3}-\frac{\partial u_3}{\partial x_1},\ \omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}$
*(Eq. 3.16, the curl of Ch. 2's (2.25))* — `ch03.curl` on a grid, `ch03.vorticity(u, x, t)` at a point.
""", where="Ch. 2 §2.9")
nb.recap("R07", "Circulation", r"""
$\Gamma\equiv\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$ *(Eq. 3.18)* — Stokes' theorem,
$\oint_C\mathbf u\cdot\mathbf t\,ds=\int_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA$ *(2.34)*: the circulation round a loop
equals the flux of vorticity through any surface it bounds, so **vorticity is circulation per unit area**,
$\mathbf n\cdot(\nabla\times\mathbf u)=\lim_{A\to0}\frac1A\oint_C\mathbf u\cdot\mathbf t\,ds$ *(2.35)*. ⚠️ Γ here is a
circulation [m²/s], not Ch. 2's shear rate. `ch03.circulation(u, ch03.planar_loop(...))` computes it. Used again in C13
and C14; Kelvin's theorem (Ch. 5) and lift (Ch. 6, 14) are built on it.
""", where="Ch. 2 §2.13")

core("C10", "Vorticity is twice the element's spin — for any pair of lines, in any frame but its own", r"""
Drop a tiny paddle wheel into a river where the water near the bank is slower. Does it turn, how fast, and does the
answer depend on which two sticks of the wheel you watch — or on whether you watch from the spinning Earth?
""")
nb.md(r"""
#### The problem in plain words

The vorticity $\boldsymbol\omega=\nabla\times\mathbf u$ is the most important derived field in atmosphere and ocean
dynamics (cyclones, eddies, potential vorticity). Its meaning is physical: a small paddle wheel carried by the flow
turns at half of it. We prove that from the two threads of C08 — and find that a flow in perfectly straight lines can
spin every element, while (C13) a flow in circles need not.
""")
nb.md(r"""
#### The idea

```
horizontal thread turns at +dβ/dt, vertical thread at −dα/dt (counterclockwise +)
spin of the element ≡ their average  = ½(−α̇ + β̇) = ½ω₃ = R₂₁/2                (Fig. 3.11)
single threads may turn at different rates (shear!) — the average of ANY perpendicular pair is the same
```
""")
P("P102", "angular velocity of a line", r"""
A line segment at angle θ (counterclockwise from +x) turns at $\dot\theta=d\theta/dt$; counterclockwise is positive. If
its tip moves relative to its tail with velocity $\delta\mathbf u$, only the part of $\delta\mathbf u$ perpendicular to
the segment turns it: $\dot\theta=(\mathbf e_\theta\cdot\delta\mathbf u)/\ell$ with $\mathbf e_\theta=(-\sin\theta,\cos\theta)$.
""", code="""
theta, ell = np.pi/2, 0.01                           # a vertical 1 cm thread
du = np.array([0.01, 0.0])                           # its tip moves 1 cm/s to the right relative to its tail
print(np.array([-np.sin(theta), np.cos(theta)]) @ du/ell)   # -1.0 rad/s: clockwise
""")
D("D12")
nb.md(r"""
**The running example: a parallel shear flow.** Near a river bank (or in the wind near the ground) the velocity is
$\mathbf u=(u_1(x_2),0)$ and locally $u_1\approx\gamma x_2$ with the shear rate $\gamma\equiv du_1/dx_2$. By
$\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}$ *(Eq. 3.16)*, $\omega_3=-\gamma$: straight
streamlines, yet nonzero vorticity — clockwise for γ > 0. How do its threads turn? (§3.5 re-states this flow as note `N34`.)
""")
D("D14")
P("P103", "rotating frame of reference", r"""
An observer on a turntable (or on the Earth) turning at Ω about z sees every point that is fixed on the turntable at
rest, although it moves with $\boldsymbol\Omega\times\mathbf x$ in the lab. At the instant the two sets of axes
coincide, velocities are related by $\mathbf u_{\rm lab}=\boldsymbol\Omega\times\mathbf x+\mathbf u_{\rm rot}$. Unlike a
Galilean frame, this correction varies from point to point — so derivatives (and the vorticity) change. Ch. 4 §4.7
develops the full rotating-frame equations.
""", code="""
Omega = np.array([0, 0, 0.5]); x = np.array([2.0, 0.0, 0.0])   # turntable at 0.5 rad/s; a point 2 m from the axis
u_lab = np.array([0.0, 1.0, 0.0])                    # a point on the disc, seen from the lab [m/s]
print(u_lab - np.cross(Omega, x))                    # [0. 0. 0.]: at rest for the rotating observer
""")
D("D13")
nb.code("""
u_lab = lambda x, t: np.array([-0.8*x[1] + 0.3*x[0]*x[1], 1.2*x[0] + 0.1*x[1]**2])   # a smooth lab field [m/s]
Om = 0.7                                                         # the observer's rotation rate Ω about z [rad/s]
u_rot = lambda x, t: u_lab(x, t) - Om*np.array([-x[1], x[0]])     # step 1: u' = u − Ω × x (in the plane)
p = np.array([0.4, -0.3])                                        # a point [m]
w_lab, w_rot = ch03.vorticity(u_lab, p, 0.0)[2], ch03.vorticity(u_rot, p, 0.0)[2]   # numerical curls ω_z, ω'_z
print(f"omega = {w_lab:.6f}, omega' = {w_rot:.6f}, omega - 2*Omega = {w_lab - 2*Om:.6f}")
assert abs(w_rot - (w_lab - 2*Om)) < 1e-8                        # ω'_z = ω_z − 2Ω (D13 step 5)
""", explain=r"The numerical curl of $\mathbf u-\boldsymbol\Omega\times\mathbf x$ equals the lab vorticity minus $2\Omega$ — the check D13 cites.")
note("N28", "Vorticity depends on the frame", r"""
(unlike S, N27): an observer rotating with the element sees ω′ = 0; one rotating at $\Omega\,\mathbf e_z$ sees
$\omega'_z=\omega_z-2\Omega$. **Climate hook:** on the Earth the vorticity measured relative to the ground (the
*relative* vorticity ζ) and the one seen from space (the *absolute* vorticity) differ by the local normal component of
the planetary vorticity, $2\Omega\sin\varphi=f$ ($1.03\times10^{-4}$ s⁻¹ at 45°N) — larger than ζ of a typical
mid-latitude cyclone (≈ 10⁻⁵–10⁻⁴ s⁻¹). Ch. 4 §4.7 and Ch. 13 build on this.
""", equation=r"\omega'_z=\omega_z-2\Omega")
note("N29", "Irrotational flow", r"""
A flow with no vorticity anywhere. Then $\mathbf u$ can be written as the gradient of a *velocity potential*,
$u_i=\partial\phi/\partial x_i$, because $\nabla\times\nabla\phi=0$ (Ch. 2 §2.13). ⚠️ The converse needs a region
without holes (*simply connected*: every loop can be shrunk to a point without leaving the region): the line vortex of
C13 ($u_\theta=B/r$, B in m²/s) has ω = 0 everywhere except its axis, yet the circulation round the axis is 2πB and the "potential" φ = Bθ jumps by
2πB after one turn. Potential flow is Ch. 6.
""", equation=r"\boldsymbol\omega=0,\quad\text{or equivalently}\quad R_{ij}=\partial u_i/\partial x_j-\partial u_j/\partial x_i=0", ref="3.17")
nb.code("""
xs, ys = sp.symbols('x y')                                   # plane coordinates
print(ch03.potential_velocity(xs**2 - ys**2, [xs, ys]))      # u = ∇φ = (2x, −2y) and its curl: 0
ng = 60 if not FAST else 40                                  # grid nodes per side (even: no node on the axis)
xg = np.linspace(-3, 3, ng)                                  # grid around a line vortex [m]
phi, jump = ch03.velocity_potential_2d(ch03.vortex_velocity_field("line", Gamma=2*np.pi, sigma=1.0), (xg, xg),
                                       x_ref=(-3.0, -3.0))   # φ by line integrals along two routes from a corner
print(f"route dependence of phi: {jump:.3f} m^2/s   vs   2*pi*B = {2*np.pi*1.0:.3f} m^2/s")
""", explain="""
1. A potential always gives an irrotational field: ∇φ of $x^2-y^2$ and its curl, 0 (sympy).
2. Around the line vortex (B = Γ/2π = 1 m²/s, irrotational except on its axis) the line integral of u·ds from a corner
   depends on the route: the two routes differ by the circulation 2πB — no single-valued potential exists in a region
   with a hole.
""")
nb.worked_example("solid-body rotation and shear", r"""
1. Solid body, ω₀ = 1 rad/s: $\mathbf G=[[0,-1],[1,0]]$ s⁻¹; $\omega_3=1-(-1)=2$ s⁻¹; every thread turns at
   $\dot\theta=1$ rad/s; spin = ½ω₃ = 1 rad/s — the element spins as fast as it revolves.
2. Shear γ = 1 s⁻¹: $\dot\theta=-\sin^2\theta$: the horizontal thread 0, the vertical −1 rad/s, a thread at 30° −0.25
   rad/s and its partner at 120° −0.75 rad/s; each pair averages −0.5 = ω₃/2.
3. Watch the solid body from a turntable at Ω = 1 rad/s: $\omega'_z=2-2\times1=0$ — the tank looks at rest.
""")
nb.code("""
Gsb = ch03.velocity_gradient_preset("solid_body_rotation", Gamma=1.0)   # G = [[0, −1], [1, 0]]: ω0 = 1 rad/s
Gsh = ch03.velocity_gradient_preset("simple_shear", Gamma=1.0)          # G = [[0, 1], [0, 0]]: γ = 1 1/s
print(ch03.element_rotation_rate(Gsb), ch03.element_rotation_rate(Gsh))   # spin = ½ω: (0, 0, 1) and (0, 0, −0.5)
th = np.deg2rad([0, 30, 90, 120])                                        # four single threads
print(np.round(ch03.material_line_rotation_rate(Gsh, th), 4))            # −γ sin²θ: 0, −0.25, −1, −0.75 rad/s
pairs = np.linspace(0, np.pi, 7)                                         # seven perpendicular pairs
print(0.5*(ch03.material_line_rotation_rate(Gsh, pairs) + ch03.material_line_rotation_rate(Gsh, pairs + np.pi/2)))
print(ch03.vorticity_in_rotating_frame(2.0, 1.0))                        # ω' = ω − 2Ω = 2 − 2 = 0 (D13)
print(ch03.parallel_shear_kinematics(1.0)["omega3"])                     # the shear flow's ω3 = −γ
""", explain="""
1. Spin = ½ω from G (D12).
2. Single threads in the shear flow turn at $-\\gamma\\sin^2\\theta$ (D14 step 3).
3. Every perpendicular pair averages −γ/2 (D14 step 6).
4. The rotating-observer rule $\\omega'_z=\\omega_z-2\\Omega$ (D13).
5. The shear flow's $\\omega_3=-\\gamma$.
""")
nb.plotly(r"""
thd = np.linspace(0, 180, 181)                                 # thread direction [°]

def rates(g):                                                  # turning rates in the shear flow with rate g [1/s]
    Gg = ch03.velocity_gradient_preset("simple_shear", Gamma=g)          # G = [[0, g], [0, 0]]
    one = ch03.material_line_rotation_rate(Gg, np.deg2rad(thd))            # a single thread
    partner = ch03.material_line_rotation_rate(Gg, np.deg2rad(thd) + np.pi/2)   # its perpendicular partner
    return {"one thread θ̇(θ)": (thd, one), "its partner θ̇(θ + 90°)": (thd, partner),
            "pair average = ω₃/2": (thd, 0.5*(one + partner))}            # name → (x, y)

fig = slider_figure(rates, "γ", np.linspace(-2, 2, 21 if not FAST else 11), unit="1/s", xlabel="thread angle θ [°]",
                    ylabel="turning rate [rad/s]", yrange=[-2.2, 2.2],     # fixed y axis
                    title="Single threads disagree, every pair agrees: spin = ω₃/2")
recolor(fig, {"one thread θ̇(θ)": COLORS["teal"], "its partner θ̇(θ + 90°)": COLORS["teal"],   # threads teal
              "pair average = ω₃/2": COLORS["orange"]}, dashes={"its partner θ̇(θ + 90°)": "dash"})   # rotation orange
fig.show()                                                     # draw it
""", explain="`material_line_rotation_rate` gives $\\dot\\theta$ for every direction at once; the slider changes the shear rate γ.")
see_read_change(
    "Two teal curves (solid and dashed) and a flat orange line between them.",
    "The two teal curves mirror each other about −γ/2 and the orange line never bends; at γ = 0 everything is zero.",
    "…the flow were solid-body rotation (not on this slider): both teal curves would be flat at ω₀ — every thread turns alike.")
whatif(r"""
…you split the relative velocity du = G·dx of C06 into the part from S and the part from R? The R part is exactly a
rigid rotation at ω/2 — C11.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C11", "Deformation plus rigid rotation: du = S·dx + ½ω × dx (3.19)", r"""
The arrow from one fluid particle to its neighbour's velocity can be split into two arrows. What are they, and why
does only one of them matter for friction?
""", eqs=("3.19",))
nb.md(r"""
#### The problem in plain words

A fluid that spins like a rigid body (a stirred cup after it settles) feels no internal friction; a fluid that is
sheared does. To write the friction law of Ch. 4 we must separate, near every point, the motion that deforms from the
motion that merely rotates.
""")
nb.md(r"""
#### The idea

```
du   =   S·dx          +     ½ ω × dx                 (3.19)
total    deformation         rigid rotation at angular velocity ω/2 (same as Ω × x with Ω = ω/2)
purple   teal                orange
```

In symbols: $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ *(3.19)*.
""")
D("D15", ref="3.19")
note("N30", "A rigid rotation at half the vorticity", r"""
The second term has the form of the rigid-body velocity $\mathbf v=\boldsymbol\Omega\times\mathbf x$ (primer P101) with
$\boldsymbol\Omega=\boldsymbol\omega/2$: near any point, the fluid turns rigidly at half its vorticity.
""")
nb.worked_example("shear γ = 1 s⁻¹, neighbour dx = (0, 1) m", r"""
1. $d\mathbf u=\mathbf G\cdot d\mathbf x=(\gamma\cdot1,0)=(1,0)$ m/s.
2. Strain part $\mathbf S\cdot d\mathbf x$ with $\mathbf S=[[0,\tfrac12],[\tfrac12,0]]$: (0.5, 0).
3. Rotation part with $\boldsymbol\omega=(0,0,-1)$:
   $\tfrac12\boldsymbol\omega\times d\mathbf x=\tfrac12(0\cdot0-(-1)\cdot1,\ (-1)\cdot0-0\cdot0,\ 0)=(0.5,0,0)$.
4. Sum (1, 0) ✓: half of the shear is stretching along 45°, half is a clockwise spin.
""")
nb.code("""
du, du_s, du_r = ch03.relative_velocity_split(Gsh, [0.0, 1.0])   # (3.19) for the shear γ = 1, dx = (0, 1) m
print(du, du_s, du_r)                                            # total, strain part, rotation part [m/s]
dx3 = np.array([0.2, -0.1, 0.3])                                 # a 3-D neighbour [m]
du, du_s, du_r = ch03.relative_velocity_split(Grand, dx3)        # the random G of C09
w = ch03.vorticity_from_gradient(Grand)                          # ω = vector of R
print(np.allclose(du_r, 0.5*np.cross(w, dx3)), np.dot(du_r, dx3))   # rotation part = ½ω × dx, perpendicular to dx
""", explain="""
1. The split for the worked example: $(1,0)=(0.5,0)+(0.5,0)$ m/s.
2. In 3-D the rotation part is exactly ½ω × dx and is perpendicular to dx (dot product ≈ 0) — it never changes the
   distance to the neighbour.
""")
nb.check_agree("""
eps = ch03.levi_civita()                                         # the ε_ijk array (Ch. 2 §2.7; three axes, Ch. 2 P73)
rot = np.zeros(3)                                                # the rotation part, built term by term [m/s]
for i in range(3):                                               # the index form −½ ε_ijk ω_k dx_j of (3.19), as loops
    for j in range(3):                                           # j: summed (repeated index)
        for k in range(3):                                       # k: summed (repeated index)
            rot[i] += -0.5*eps[i, j, k]*w[k]*dx3[j]              # one of the 27 terms
assert np.allclose(rot, du_r)                                    # = the library's rotation part
assert np.allclose(ch03.strain_rate_tensor(Grand) @ dx3 + rot, Grand @ dx3)   # S·dx + rotation = G·dx
print("triple loop:", rot)
""")
nb.md(r"""
The loops implement $du_i=\big(S_{ij}-\tfrac12\varepsilon_{ijk}\omega_k\big)dx_j$ — the first form of
$du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ *(Eq. 3.19)* in D15.
""")
nb.figure(r"""
from scripts.ch03_drawings import ring_arrows                    # (no physics inside)
fig, axs = plt.subplots(1, 3, figsize=(12, 4))                   # total, strain part, rotation part
for ax, part, ttl in zip(axs, ("total", "strain", "rotation"),
                         ("total du = G·dx", "strain part S·dx", "rotation part ½ω × dx")):   # one panel per part
    ring_arrows(ax, Gsh, radius=0.1, parts=(part,), n=24, scale=0.6)   # 24 neighbours 10 cm away, shear γ = 1
    ax.plot(0, 0, "o", color=COLORS["ink"]); ax.set_aspect("equal")   # the centre point O
    ax.set_xlim(-0.17, 0.17); ax.set_ylim(-0.17, 0.17); ax.set_title(ttl); ax.set_xlabel("$dx_1$ [m]")   # same axes
axs[0].set_ylabel("$dx_2$ [m]")                                  # shared y label
plt.suptitle("Simple shear = pure strain (teal) + clockwise rigid spin at 0.5 rad/s (orange)")   # the message
savefig(fig, "ch03", "ring_split"); plt.show()                   # save to outputs/ch03 and draw
""", see="Purple arrows that are all horizontal; teal arrows that fan out along 45° and in along 135°; orange arrows that circulate clockwise.",
    read="""The shear's horizontal arrows are, arrow by arrow, the sum of the teal and orange ones: a pure strain that would
turn the ring into an ellipse, plus a rigid clockwise spin at 0.5 rad/s that turns it.""",
    change="…G = [[0, −1], [1, 0]] (solid body): the teal arrows would vanish and the purple arrows would be the orange ones.")
note("N33", "Summary of §3.4", r"""
The relative motion near a point = a rigid rotation of the element (at ω/2) + a deformation (S), which itself = pure
stretching along three perpendicular principal axes (next block).
""")
whatif(r"""
…you turned the axes to where S has no off-diagonal entries? Then the deformation is three pure stretchings — a small
sphere becomes an ellipsoid (C12).
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C12", "Principal strain axes: a small sphere becomes an ellipsoid", r"""
A round drop of dye is released into a flow. What shape is it a moment later, and which way do its long and short axes
point?
""")
nb.md(r"""
#### The problem in plain words

Satellite images of ocean colour show round patches of plankton pulled into ellipses and then filaments; weather fronts
form where the wind's deformation squeezes temperature contours together (*frontogenesis*, Ch. 13). Both are the
principal axes of the strain-rate tensor at work.
""")
nb.md(r"""
#### The idea

```
in the principal frame (overbar) S is diagonal:  dū₁ = S̄₁₁ dx̄₁,  dū₂ = S̄₂₂ dx̄₂,  dū₃ = S̄₃₃ dx̄₃     (3.21)
each principal direction stretches in proportion to its own length  →  circle → ellipse on those axes (Fig. 3.13)
```

In symbols: $d\bar u_\alpha=\bar S_{\alpha\alpha}\,d\bar x_\alpha$ (no sum) *(3.21)*.
""")
note("N31", "The strain part in the principal frame", r"""
In the frame of the principal axes of S (Ch. 2 §2.11: a symmetric tensor has three perpendicular eigenvectors, primer
P80; ⚠️ the book's "Section 2.12" reference here means §2.11), the strain part of
$du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ *(3.19)* is diagonal:
""", equation=r"d\bar{\mathbf u}=\bar{\mathbf S}\cdot d\bar{\mathbf x}=\begin{bmatrix}\bar S_{11}&0&0\\0&\bar S_{22}&0\\0&0&\bar S_{33}\end{bmatrix}\begin{bmatrix}d\bar x_1\\d\bar x_2\\d\bar x_3\end{bmatrix}", ref="3.20")
note("N32", "Its three components", r"""
with $\bar S_{\alpha\alpha}$ the eigenvalues of S (Greek index: no sum):
""", equation=r"d\bar u_1=\bar S_{11}d\bar x_1,\quad d\bar u_2=\bar S_{22}d\bar x_2,\quad d\bar u_3=\bar S_{33}d\bar x_3", ref="3.21")
P("P104", "linear map of a circle is an ellipse", r"""
Multiplying every point of a unit circle by a matrix M gives an ellipse. For a symmetric M its axes are M's
eigenvectors and its semi-axes the eigenvalues; for a general M they are the *singular values* and singular vectors
(`np.linalg.svd`) — a gloss we need only for the finite-time caveat (an ellipse with semi-axes $a_\alpha$ along the
coordinate axes obeys $\sum_\alpha X_\alpha^2/a_\alpha^2=1$).
""", code="""
M = np.array([[1.1, 0.0], [0.0, 0.9]])               # stretch x by 10 %, squeeze y by 10 %
s = np.linspace(0, 2*np.pi, 400); Pc = M @ np.stack([np.cos(s), np.sin(s)])   # the mapped unit circle
print(round(Pc[0].max(), 4), round(Pc[1].max(), 4))  # 1.1 0.9: the semi-axes
""")
D("D16")
nb.worked_example("shear γ = 1 s⁻¹, a 1 mm drop after 0.1 s", r"""
1. Principal rates ±γ/2 = ±0.5 s⁻¹ at +45° and −45° (Ch. 2 Ex. 2.4, now with γ = 2S₁₂).
2. First order (D16): semi-axes $1\times(1\pm0.5\times0.1)=1.05$ and 0.95 mm.
3. The strain acting alone: $e^{\pm0.05}=1.0513$ and 0.9512 mm.
4. Strain + rotation together (the true shear flow): singular values of $e^{\mathbf Gt}=[[1,0.1],[0,1]]$: 1.0512 and
   0.9512 mm, long axis at 43.6° — turned 1.4° clockwise off 45°. The rotation part alone would turn it by 0.5 rad/s × 0.1 s ≈ 2.9°, but
   the strain keeps pulling the long axis back toward 45°, so the net early turning rate is about γ/4, not γ/2.
""")
nb.code("""
lam, axes = ch03.principal_strain_rates(Gsh)                     # eigenvalues (ascending) and unit eigenvectors of S
print(lam, np.round(axes[:, 1], 4))                              # [−0.5, 0.5]; the stretching axis (1, 1)/√2
for m in ("first_order", "strain_only", "exact"):                # three answers to "what shape after 0.1 s?"
    a, d = ch03.strain_ellipse_axes(Gsh, 0.1, method=m)          # semi-axes (descending) and their directions
    print(f"{m:12s} semi-axes {np.round(a, 4)}, long axis at {np.degrees(np.arctan2(d[1, 0], d[0, 0])):.1f}°")
du_bar, dx_bar = ch03.strain_velocity_principal(Gsh, [0.001, 0.0])   # (3.21) in the eigenframe
print(du_bar/dx_bar)                                             # each component: its own eigenvalue
""", explain="""
1. Eigenvalues and eigenvectors of S (Ch. 2 P80), ordered like λ; `np.arctan2` (Ch. 2 P70) turns a direction into an angle.
2. Three answers to "what shape after 0.1 s": the book's first-order statement (1 ± λt), the strain alone ($e^{\\lambda t}$)
   and the exact linear flow (singular values of $e^{\\mathbf Gt}$).
3. In the eigenframe each velocity component is its own eigenvalue × its own coordinate,
   $d\\bar u_\\alpha=\\bar S_{\\alpha\\alpha}\\,d\\bar x_\\alpha$ *(3.21)*.
""")
nb.animation(r"""
nfr = 48 if not FAST else 24                                     # frames
times = np.linspace(0, 3, nfr)                                   # t from 0 to 3 s in the shear γ = 1
fig, ax = plt.subplots(figsize=(5.2, 4.6))                      # one panel, fixed limits (never autoscale)
ax.set_aspect("equal"); ax.set_xlim(-3.2, 3.2); ax.set_ylim(-2.2, 2.2)   # a 1 mm drop, axes in mm
L = 2.0                                                          # half-length of the drawn axes [mm]
ax.plot([-L, L], [-L, L], "--", color=COLORS["blue"], lw=1.2, label="stretching axis of S (45°)")      # eigenvector of +γ/2
ax.plot([-L, L], [L, -L], "--", color=COLORS["rose"], lw=1.2, label="compressing axis of S (−45°)")   # eigenvector of −γ/2
(ring,) = ax.plot([], [], "o", color=COLORS["teal"], ms=3, label="72 tracers")          # updated every frame
(axis_line,) = ax.plot([], [], color=COLORS["accent"], lw=2.2, label="current long axis")   # updated every frame
txt = ax.text(-3.0, 1.9, "", fontsize=9)                         # readout
ax.legend(fontsize=7, loc="lower right"); ax.set_xlabel("$x_1$ [mm]"); ax.set_ylabel("$x_2$ [mm]")
ax.set_title("A circle in simple shear: the ellipse starts on S's axes, then tips")

def update(i):                                                   # draw frame i
    t = times[i]                                                 # time of this frame [s]
    c = ch03.deform_circle(Gsh, t, n=72)                         # tracers of a unit circle carried by e^{Gt}
    ring.set_data(c[0], c[1])                                    # move the dots
    a, d = ch03.strain_ellipse_axes(Gsh, t, method="exact")      # exact semi-axes and directions (SVD)
    axis_line.set_data([-a[0]*d[0, 0], a[0]*d[0, 0]], [-a[0]*d[1, 0], a[0]*d[1, 0]])   # the long axis
    txt.set_text(f"t = {t:.2f} s   a = {a[0]:.2f}, b = {a[1]:.2f} mm   long axis {np.degrees(np.arctan2(d[1, 0], d[0, 0])):.1f}°")
    return ring, axis_line, txt                                  # the artists that changed

show_animation(animate(update, frames=nfr, fig=fig, interval=80), player="video")   # smooth MP4
""", explain="""
1. `deform_circle` carries 72 tracers of a circle with the exact linear map $e^{\\mathbf Gt}$.
2. `strain_ellipse_axes(method="exact")` gives the ellipse's semi-axes and long-axis direction at each time.
""")
see_read_change(
    "A ring of teal dots becoming an ever longer ellipse; a purple long axis that starts on the blue 45° line and then tips toward the horizontal.",
    """At first the ellipse grows along the blue axis, as D16 says; as time goes on the purple axis tips toward the flow
direction (the rotation part keeps turning it) — the principal-axis statement is about the *first instant*.""",
    "…the flow were pure strain (the 'irrotational_strain' preset): the purple and blue axes would stay together for ever.")
nb.explainer("spin_and_principal_axes", "Can a straight flow make a fluid element spin?", r"""
**Why interactive:** the claim is "for every pair of lines" — a figure can show a few pairs, but sweeping the pair angle yourself shows the average never moves while the single lines do.
Two perpendicular threads and a small paddle wheel ride in a flow you choose. Drag the pair's angle: the two threads
turn at different rates but their average stays at ½ω₃. Drag a probe round the ring: its relative velocity splits into
a teal strain arrow and an orange rotation arrow, $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$
*(Eq. 3.19)*. The circle turns into an ellipse on the principal axes.
""", tries=[
    "Preset 'parallel shear γ = 1': sweep the pair angle from 0 to 180° and watch the average line stay flat at −0.5 rad/s.",
    "Preset 'rotate with the element' (Ω = ω/2): the paddle wheel stops — ω′ = 0.",
    "Preset 'pure strain': the orange arrows vanish; the ellipse axes never move.",
    "Derivation tab, D15: at the ε-swap step the orange arrow flips sign on screen.",
])
whatif(r"""
…the streamlines were circles instead of straight lines? Then "going round" and "spinning" come apart completely —
§3.5's two vortices.
""")

# =====================================================================================================================
# A.5 §3.5 Kinematics of Simple Plane Flows — N34, N35 (head), C13, C14
# =====================================================================================================================
nb.section("3.5", "Kinematics of Simple Plane Flows", intro="""
**What is this section about?** Two families of plane flow where one coordinate does all the work: the parallel shear
flow (straight streamlines, spinning elements) and circular flows — solid-body rotation, the irrotational vortex and
the realistic vortices between them.
""")
note("N34", "Parallel shear flow", r"""
$\mathbf u=(u_1(x_2),0)$ with the shear rate $\gamma(x_2)\equiv du_1/dx_2$ — the flow you used in C10. By
$\omega_3=\frac{\partial u_2}{\partial x_1}-\frac{\partial u_1}{\partial x_2}$ *(Eq. 3.16)*, $\omega_3=-\gamma$. The
vertical thread AB turns at −γ, the horizontal BC at 0, their average −γ/2 = ω₃/2 — and D14 (in C10) showed that
*every* perpendicular pair gives −γ/2. ⚠️ ω₃ = −γ is clockwise for γ > 0. Couette flow (Ch. 8) and boundary layers
(Ch. 9) are this flow.
""", equation=r"\gamma(x_2)\equiv du_1/dx_2,\qquad \omega_3=-\gamma")
note("N35", "Two elements in the same shear (Fig. 3.14)", r"""
For the square ABCD with sides along the axes, S has only off-diagonal entries — it shears without stretching its sides.
Turned by 45°, onto the principal axes (⚠️ Ch. 2's Ex. 2.4 had $S_{12}=\Gamma$; here $S_{12}=\gamma/2$), S is diagonal —
the square PQRS stretches along $\bar x_1$ and is compressed along $\bar x_2$ (eigenvalue −γ/2, i.e. compression at rate
γ/2) while its corners stay right angles. Both spin at −γ/2.
""", equation=r"S_{ij}=\begin{bmatrix}0&\gamma/2\\\gamma/2&0\end{bmatrix},\qquad \bar S_{ij}=\begin{bmatrix}\gamma/2&0\\0&-\gamma/2\end{bmatrix}")
nb.animation(r"""
from scripts.ch03_drawings import shear_elements_frames      # corners of ABCD and PQRS carried by u = (γy, 0)
nfr = 16 if not FAST else 10                                  # frames to step through
frames = shear_elements_frames(1.0, np.linspace(0, 0.75, nfr))   # γ = 1 1/s, t = 0 … 0.75 s
fig, axs = plt.subplots(1, 2, figsize=(9, 3.9))              # ABCD left, PQRS right
arts = []                                                     # the artists each frame updates
for ax, key, col in zip(axs, ("ABCD", "PQRS"), (COLORS["accent"], COLORS["teal"])):
    ax.set_aspect("equal"); ax.set_xlim(-1.0, 1.0); ax.set_ylim(-0.8, 0.8); ax.set_xlabel("$x_1$ [m]")   # fixed limits
    ax.plot(*np.hstack([frames[0][key], frames[0][key][:, :1]]), color=COLORS["grid"], lw=1)   # the start shape
    (ln,) = ax.plot([], [], color=col, lw=2.4)                # the element now
    (pad,) = ax.plot([], [], color=COLORS["ink"], lw=1.5)     # a paddle at its centre, turned by −γt/2
    tx = ax.text(-0.95, 0.62, "", fontsize=8)                 # numbers of this frame
    arts.append((key, ln, pad, tx))
axs[0].set_ylabel("$x_2$ [m]")                                # shared y label
axs[0].set_title("ABCD: sides along the axes (shear)"); axs[1].set_title("PQRS: sides at 45° (stretch + squeeze)")

def update(i):                                                # draw frame i
    f = frames[i]                                             # corners, sides and angles at this time
    spin = -0.5*1.0*f["t"]                                    # element spin −γ/2 times t [rad]
    for key, ln, pad, tx in arts:                             # both elements
        P = f[key]; ln.set_data(*np.hstack([P, P[:, :1]]))    # closed outline through the four corners
        c, s = np.cos(spin), np.sin(spin)                     # paddle: two crossing blades of half-length 0.15 m
        pad.set_data([-0.15*c, 0.15*c, np.nan, 0.15*s, -0.15*s], [-0.15*s, 0.15*s, np.nan, -0.15*c, 0.15*c])
        sides = ", ".join(f"{v:.3f}" for v in f[key + "_sides"]); ang = f[key + "_angles_deg"]   # four sides [m], angles [°]
        tx.set_text(f"t = {f['t']:.2f} s\nsides {sides} m\ncorner {ang[1]:.1f}°")   # corner B (or Q)
    return [a for _, ln, pad, tx in arts for a in (ln, pad, tx)]   # the artists that changed

show_animation(animate(update, frames=nfr, fig=fig, interval=500), player="frames")   # step frame by frame
""", explain="""
1. `shear_elements_frames` carries the corners of both squares with the exact linear map of the shear flow and reports
   side lengths and corner angles.
2. Each frame redraws the two elements, a paddle turned by −γt/2, and the numbers.
""")
see_read_change(
    "Left, a square whose top slides right (a parallelogram); right, a diamond that stretches along 45° and shrinks along −45°. Both paddles turn clockwise together.",
    """ABCD's vertical sides lengthen only at second order while its corner angle falls at γ per second; PQRS's corner
stays 90.0° to first order while one pair of sides grows and the other shrinks. Step frame by frame: the numbers are the
linear (C07), shear (C08) and spin (C10) rates made visible.""",
    "…you ran longer (t ≫ 1 s): PQRS's angles would drift from 90° too — the finite-time effect of C12's animation.")

core("C13", "Going round is not spinning: vorticity in polar coordinates (3.23)", r"""
Two tanks of water both go round in circles — one spun up like a merry-go-round, one draining through a plughole. Put a
small paddle wheel in each. Which one turns?
""", eqs=("3.23",))
nb.md(r"""
#### The problem in plain words

Hurricanes, bathtub drains, the flow round a stirred cup: circular streamlines everywhere. The vorticity tells us
whether each little parcel spins as it goes round. We need it in polar coordinates, where these flows are simple.
""")
nb.md(r"""
#### The idea

```
solid body  u_θ = ω₀ r :   Γ grows like r²  → ω_z = 2ω₀ everywhere → the wheel spins as it orbits      (Fig. 3.15)
line vortex u_θ = B/r  :   Γ = 2πB for every circle → ω_z = 0 except on the axis → the wheel keeps its heading (Fig. 3.16)
vorticity = circulation per unit area of a small polar sector  →  ω_z = (1/r)∂(r u_θ)/∂r − (1/r)∂u_r/∂θ   (3.23)
```

In symbols: $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ *(3.23)*.
""")
P("P105", "polar coordinates as a moving basis", r"""
In the plane, $\mathbf e_r=(\cos\theta,\sin\theta)$ and $\mathbf e_\theta=(-\sin\theta,\cos\theta)$ change direction
with θ. A small polar sector has area $r\,dr\,d\theta$ (an arc of length $r\,d\theta$ times a width dr); along a circle
of radius r the line element is $d\mathbf s=r\,d\theta\,\mathbf e_\theta$, along a ray $d\mathbf s=dr\,\mathbf e_r$. Arcs
at r and r + dr have *different* lengths — the source of the extra 1/r terms in polar formulas.
""", code="""
r, dr, dth = 2.0, 0.01, 0.02                          # a sector at r = 2 m, 1 cm deep, 0.02 rad wide
inner, outer = r*dth, (r + dr)*dth                    # arc lengths [m]
print(round(outer - inner, 6), dr*dth)                # 0.0002 0.0002: the arcs differ by dr·dθ
print(r*dr*dth)                                       # 0.0004 m²: the sector's area
""")
D("D17", ref="3.23", check_src=r"""
r, th = sp.symbols('r theta', positive=True)          # polar coordinates
x, y = sp.symbols('x y', real=True)                   # Cartesian coordinates
ur, uth = r**2*sp.sin(th), r*sp.cos(th)               # a non-axisymmetric test field (u_r, u_θ)
polar = ch03.polar_vorticity_z_sym(ur, uth, r, th)    # (3.23): (1/r)∂(r u_θ)/∂r − (1/r)∂u_r/∂θ
ux = ur*sp.cos(th) - uth*sp.sin(th)                   # Cartesian components u = u_r e_r + u_θ e_θ
uy = ur*sp.sin(th) + uth*sp.cos(th)
to_xy = {r: sp.sqrt(x**2 + y**2), th: sp.atan2(y, x)} # write them as functions of (x, y)
curl = sp.diff(uy.subs(to_xy), x) - sp.diff(ux.subs(to_xy), y)   # (3.16): ω_z = ∂v/∂x − ∂u/∂y
back = {x: r*sp.cos(th), y: r*sp.sin(th)}             # return to polar coordinates to compare
print(sp.simplify(sp.expand_trig(curl.subs(back)) - polar))   # 0: the sector derivation agrees with the curl
print(sp.simplify(polar))                             # the vorticity of this test field
""", extra_check="The sympy cell below compares $\\omega_z=\\frac1r\\frac{\\partial}{\\partial r}(ru_\\theta)-\\frac1r\\frac{\\partial u_r}{\\partial\\theta}$ (3.23) with the Cartesian curl for the non-axisymmetric pair $u_r=r^2\\sin\\theta$, $u_\\theta=r\\cos\\theta$.")
note("N36", "Solid-body rotation", r"""
— a tank spun steadily until the fluid turns with it. By
$\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ *(3.23)*,
$\omega_z=\frac1r\frac{d}{dr}(\omega_0r^2)=2\omega_0$ everywhere: each element spins about its own centre at the rate
it revolves round the axis; S = 0, nothing deforms (Fig. 3.15).
""", equation=r"u_r=0\quad\text{and}\quad u_\theta=\omega_0r", ref="3.22")
note("N37", "Its circulation round a centred circle", r"""
= vorticity 2ω₀ × area πr². It holds for *any* circuit, centred or not: with uniform ω, Stokes'
$\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$ *(3.18)* gives Γ = ω × enclosed area.
Number: ω₀ = 1 s⁻¹, r = 1 m → Γ = 2π = 6.283 m²/s; an off-centre circle of radius 1 m centred at (2, 0) m also gives
6.283 m²/s (checked below).
""", equation=r"\Gamma=\oint_C\mathbf u\cdot d\mathbf s=\int_0^{2\pi}u_\theta r\,d\theta=2\pi ru_\theta=2\pi r^2\omega_0", ref="3.24")
note("N38", "The irrotational (line) vortex", r"""
— the ideal limit of a drain or a tornado far from its core. (⚠️ Fig. 3.16 writes C for B.) By (3.23),
$\omega_z=\frac1r\frac{d}{dr}(r\cdot B/r)=\frac1r\frac{dB}{dr}=0$ for every r > 0.
""", equation=r"u_r=0\quad\text{and}\quad u_\theta=B/r", ref="3.25")
note("N39", "Yet the circulation round any centred circle is the same nonzero number", r"""
independent of r ($u_\theta r=B$ is constant).
""", equation=r"\Gamma=\int_0^{2\pi}u_\theta r\,d\theta=2\pi ru_\theta=2\pi B", ref="3.26")
note("N40", "All the vorticity sits on the axis", r"""
Taking vorticity as circulation per unit area in a shrinking disc gives a value that is infinite at r = 0 with a finite
area integral 2πB: a *delta function* (Ch. 2 gloss: an infinitely concentrated amount whose total is finite). The point
vortex of Ch. 6 is exactly this.
""", equation=r"[\omega_z]_{r\to0}=\lim_{r\to0}\frac1A\int_A\omega_z\,dA=\lim_{r\to0}\frac{1}{\pi r^2}\oint_C\mathbf u\cdot d\mathbf s=\lim_{r\to0}\frac{2B}{r^2}", ref="3.27")
note("N41", "A loop that does not enclose the axis has zero circulation", r"""
— the sector ABCD of Fig. 3.16: the radial legs BC and DA contribute nothing ($\mathbf u\perp d\mathbf s$), and the arcs
cancel because $u_\theta r=B$. Elements in this flow deform but do not spin.
""", equation=r"\Gamma_{ABCD}=-[u_\theta r]_r\Delta\theta+[u_\theta r]_{r+\Delta r}\Delta\theta=0")
nb.worked_example("the two tanks with easy numbers", r"""
1. Solid body ω₀ = 1 s⁻¹: $u_\theta=r$; $\omega_z=\frac1r\frac{d(r\cdot r)}{dr}=\frac{2r}{r}=2$ s⁻¹; spin ½ω_z = 1 rad/s =
   orbit rate $u_\theta/r$ = 1 rad/s.
2. Line vortex B = 1 m²/s: $u_\theta=1/r$; $\omega_z=\frac1r\frac{d(r\cdot1/r)}{dr}=\frac1r\frac{d(1)}{dr}=0$; spin 0, orbit
   rate $1/r^2$ rad/s.
3. Circulations: circle r = 1 m: solid body $2\pi(1)^2(1)=6.283$ m²/s; line vortex $2\pi(1)=6.283$ m²/s — equal at
   r = 1, but at r = 2 m: 25.13 vs 6.283 m²/s.
4. Mean vorticity in the disc r = 0.1 m for the line vortex: $2B/r^2=200$ s⁻¹; at r = 0.01 m: 20 000 s⁻¹.
""")
nb.code("""
r, th, w0, B = sp.symbols('r theta omega_0 B', positive=True)          # symbols for (3.23)
print(ch03.polar_vorticity_z_sym(0, w0*r, r, th), ch03.polar_vorticity_z_sym(0, B/r, r, th))   # 2ω0 and 0
print(ch03.polar_vorticity_z(lambda r, t: 0.0, lambda r, t: 1.0*r, 1.3, 0.4),        # stencils in r and θ: 2
      ch03.polar_vorticity_z(lambda r, t: 0.0, lambda r, t: 1.0/r, 1.3, 0.4))        # … and ≈ 0
kw = dict(Gamma=2*np.pi, sigma=1.0)                  # "solid": ω0 = Γ/(2πσ²) = 1 1/s; "line": B = Γ/2π = 1 m²/s
print(ch03.circulation_circle("solid", 1.0, **kw), ch03.circulation_circle("solid", 1.0, center=(2.0, 0.0), **kw))
print([round(ch03.circulation_circle("line", r_, **kw), 4) for r_ in (0.5, 1.0, 2.0)])   # 2πB at every radius
print(ch03.annular_sector_circulation("line", 1.0, 0.1, 0.2, **kw))   # a sector away from the axis: 0
print([ch03.mean_vorticity_in_disc("line", r_, **kw) for r_ in (0.1, 0.01)])   # 2B/r²: grows without bound
""", explain="""
1. sympy applies $\\omega_z=\\frac1r\\frac{\\partial}{\\partial r}(ru_\\theta)-\\frac1r\\frac{\\partial u_r}{\\partial\\theta}$
   *(3.23)* to both profiles: $2\\omega_0$ and 0.
2. Stencils in r and θ agree (2.0 and ≈ 10⁻¹²).
3. Γ = ω × area even off-centre (N37): 6.283 m²/s for both circles.
4. 2πB at every radius (N39).
5. Zero for a sector away from the axis (N41).
6. The delta-function core: the mean vorticity in a disc grows like $1/r^2$ (N40).
""")
nb.check_agree("""
r0, dr, dth = 1.3, 0.01, 0.02                         # a small polar sector inside a Rankine core
kwR = dict(Gamma=2*np.pi, sigma=2.0)                  # Γ = 2π m²/s, core radius σ = 2 m
uth = lambda r: ch03.rankine_vortex(r, **kwR)[0]      # u_θ(r) of the Rankine vortex [m/s]
outer = uth(r0 + dr)*(r0 + dr)*dth                    # outer arc, counterclockwise: u_θ (r + dr) dθ
inner = -uth(r0)*r0*dth                               # inner arc, run backwards
legs = 0.0                                            # radial legs: u_r = 0, so u·ds = 0 there
area = (r0 + dr/2)*dr*dth                             # exact sector area ((r+dr)² − r²)dθ/2 [m²]
mine = (outer + inner + legs)/area                    # circulation / area = ω_z (D17 with numbers)
lib = ch03.polar_vorticity_z(None, "rankine", r0, 0.0, **kwR)     # (3.23) by stencils
cart = ch03.vorticity(ch03.vortex_velocity_field("rankine", **kwR), np.array([r0*np.cos(0.4), r0*np.sin(0.4)]), 0.0)[2]
assert np.allclose(mine, lib, rtol=1e-3) and np.allclose(lib, cart, rtol=1e-6)   # all = Γ/πσ² = 0.5 1/s
print(f"four legs: {mine:.5f}   polar formula: {lib:.5f}   Cartesian curl: {cart:.5f}   (Γ/πσ² = {2*np.pi/(np.pi*4):.5f})")
""")
nb.md("D17 with numbers: four line integrals and one division.")
nb.figure(r"""
rr = np.geomspace(1e-3, 10, 60)                        # disc radii from 1 mm to 10 m
kw = dict(Gamma=2*np.pi, sigma=1.0)
fig, ax = plt.subplots(figsize=(7, 4))                  # log–log axes (Ch. 1 P13)
ax.loglog(rr, ch03.mean_vorticity_in_disc("line", rr, **kw), color=COLORS["orange"], lw=2.3, label="line vortex: 2B/r² (slope −2)")
ax.loglog(rr, ch03.mean_vorticity_in_disc("solid", rr, **kw), color=COLORS["teal"], lw=2.3, label="solid body: 2ω₀")   # flat
ax.loglog(rr, ch03.mean_vorticity_in_disc("rankine", rr, **kw), "--", color=COLORS["accent"], lw=2.3, label="Rankine, σ = 1 m")
ax.axvline(1.0, color=COLORS["muted"], ls=":", lw=1)    # the core radius σ
ax.set_xlabel("disc radius $r$ [m]"); ax.set_ylabel(r"mean vorticity $\Gamma(r)/\pi r^2$ [1/s]"); ax.legend(fontsize=8)
ax.set_title("An irrotational vortex keeps all its vorticity on the axis")
savefig(fig, "ch03", "mean_vorticity_disc"); plt.show() # save to outputs/ch03 and draw
""", see="A straight orange line of slope −2, a flat teal line, and a dashed purple curve that switches from one to the other at r = σ = 1 m.",
    read="""The mean vorticity in a disc is its circulation divided by its area. For the line vortex the circulation is the
same for every disc, so the smaller the disc, the larger the average — all of it sits on the axis. A real vortex has a
finite core, inside which it behaves like a solid body.""",
    change="…σ shrank toward zero: the purple curve's flat part would climb (Γ/πσ²) and the Rankine vortex would approach the line vortex.")
nb.animation(r"""
nfr = 60 if not FAST else 30                           # frames
times = np.linspace(0, 2*np.pi, nfr)                   # one revolution of the solid body [s]
radii = np.array([0.5, 1.0, 1.5, 2.0])                 # four paddle wheels [m]
cases = {"solid body, ω₀ = 1 1/s": (lambda r: 1.0*np.ones_like(r), 2.0),      # (orbit rate u_θ/r [rad/s], ω_z [1/s])
         "line vortex, B = 1 m²/s": (lambda r: 1.0/r**2, 0.0)}
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.6))       # solid body left, line vortex right
arts = []                                              # the artists each frame updates
for ax, (name, (Om, wz)) in zip(axs, cases.items()):
    ax.set_aspect("equal"); ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5); ax.set_title(name)   # fixed limits
    for R_ in radii: ax.add_patch(plt.Circle((0, 0), R_, fill=False, color=COLORS["grid"]))    # the streamlines (circles)
    (wheels,) = ax.plot([], [], color=COLORS["ink"], lw=1.6)   # all four paddle wheels in one line
    (el,) = ax.plot([], [], color=COLORS["teal"], lw=2)        # the polar element's outline
    arts.append((Om, wz, wheels, el))
axs[0].set_xlabel("$x$ [m]"); axs[1].set_xlabel("$x$ [m]"); axs[0].set_ylabel("$y$ [m]")
sq = np.linspace(0, 1, 15)                             # points along each edge of the element

def update(i):                                         # draw frame i
    t = times[i]                                       # time of this frame [s]
    out = []                                           # artists that changed
    for Om, wz, wheels, el in arts:                    # both vortices
        ang = Om(radii)*t                              # each wheel's position angle: it orbits at u_θ/r
        spin = 0.5*wz*t                                # each wheel turns at ½ω_z (D12, D17)
        xs, ys = [], []                                # blade end points, NaN-separated
        for R_, a in zip(radii, ang):                  # one wheel per radius
            cx, cy = R_*np.cos(a), R_*np.sin(a)        # its hub
            for b in (spin, spin + np.pi/2):           # two crossing blades
                xs += [cx - 0.22*np.cos(b), cx + 0.22*np.cos(b), np.nan]; ys += [cy - 0.22*np.sin(b), cy + 0.22*np.sin(b), np.nan]
        wheels.set_data(xs, ys)                        # move the wheels
        rr_ = np.r_[1 + 0.4*sq, np.full(15, 1.4), 1.4 - 0.4*sq, np.full(15, 1.0)]   # element outline: r from 1 to 1.4 m
        tt_ = np.r_[np.zeros(15), 0.35*sq, np.full(15, 0.35), 0.35 - 0.35*sq]      # … and θ from 0 to 0.35 rad
        th_now = tt_ + Om(rr_)*t                       # every boundary point orbits at its own rate
        el.set_data(rr_*np.cos(th_now), rr_*np.sin(th_now))   # back to x, y
        out += [wheels, el]
    return out

show_animation(animate(update, frames=nfr, fig=fig, interval=70), player="video")   # smooth MP4
""", explain="""
1. Every point moves on its circle at the orbit rate $u_\\theta/r$ (1 rad/s for the solid body, $1/r^2$ for the line vortex).
2. Each paddle wheel is turned by $\\tfrac12\\omega_z t$ — the element spin of D12, with $\\omega_z$ from (3.23).
3. The teal polar element is drawn by carrying 60 points of its outline.
""")
see_read_change(
    "Left: four wheels going round together, each turning once per revolution, and a teal element that keeps its shape. Right: wheels orbiting at very different speeds but never turning, and an element sheared into a thin wedge.",
    """On the left, going round and spinning go together (ω_z = 2ω₀). On the right the wheels keep pointing the same way
while they circle — ω_z = 0 — although the fluid clearly deforms: fast inside, slow outside.""",
    "…a wheel were added at r = 0.2 m on the right: it would race round (25 rad/s) — and still not turn.")
whatif(r"""
…a real drain vortex were measured? Near the axis the speed does not blow up like 1/r; the core turns like a solid
body and the outside like a line vortex — the Rankine and Gaussian models of C14.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C14", "Real vortices have a core: the Rankine and Gaussian models (3.28)–(3.29)", r"""
A tornado's wind is calm at its centre, fastest a little way out and weaker far away. Which simple formula reproduces
that profile, and where exactly is the fastest wind?
""", eqs=("3.28", "3.29"))
nb.md(r"""
#### The problem in plain words

Tornado chasers, hurricane forecasters and aircraft-wake engineers all quote two numbers: the peak wind and the radius
where it occurs. Pure solid-body rotation grows for ever; the line vortex blows up at the axis. Joining them gives a
vortex with a spinning core and an irrotational skirt — the model behind every cyclone schematic.
""")
nb.md(r"""
#### The idea

```
              core (r < σ): solid body, ω_z = Γ/πσ²       outside: irrotational, u_θ = Γ/2πr
Rankine:   sharp join at r = σ (speed continuous, vorticity jumps)     → peak at r = σ exactly
Gaussian:  smooth join, ω_z = (Γ/πσ²) e^{−r²/σ²}                        → peak at r ≈ 1.1209 σ
```
""")
note("N42", "Real vortices", r"""
— bathtub drains, wing-tip vortices (Ch. 14), tornadoes, tropical cyclones (Ch. 13) — combine a nearly solid-body core
with a nearly irrotational outer flow and have bounded speeds.
""")
nb.md(r"""
#### The maths — the Rankine vortex

Uniform vorticity inside a core of radius σ, none outside; Γ is the total circulation:

$$\omega_z(r)=\begin{Bmatrix}\Gamma/\pi\sigma^2=\text{const.}&r\le\sigma\\0&r>\sigma\end{Bmatrix}\quad\text{and}\quad u_\theta(r)=\begin{Bmatrix}(\Gamma/2\pi\sigma^2)r&r\le\sigma\\\Gamma/2\pi r&r>\sigma\end{Bmatrix}\qquad(3.28)$$

A piecewise formula is evaluated with `np.where(r <= sigma, inside, outside)` (Ch. 1 P46).
""")
D("D18", ref="3.28")
P("P106", "substitution in an integral", r"""
Replace a messy variable by a simpler one and convert dx too: with $s=r^2/\sigma^2$, $ds=2r\,dr/\sigma^2$, so
$2\pi r\,dr=\pi\sigma^2\,ds$ and the limits r = 0 … R become s = 0 … R²/σ². The value of the integral does not change.
""", code="""
from scipy.integrate import quad                     # adaptive quadrature (P87)
sig, R = 1.0, 1.5                                    # core radius and outer radius [m]
lhs = quad(lambda r: np.exp(-r**2/sig**2)*2*np.pi*r, 0, R)[0]      # in r
rhs = quad(lambda s: np.exp(-s)*np.pi*sig**2, 0, R**2/sig**2)[0]   # in s = r²/σ²
print(round(lhs, 4), round(rhs, 4))                  # 2.8105 2.8105 (= π(1 − e^−2.25))
""")
D("D19", ref="3.29", check_src=r"""
r, sig, Gam = sp.symbols('r sigma Gamma', positive=True)            # radius, core radius, circulation
uth = Gam/(2*sp.pi*r)*(1 - sp.exp(-r**2/sig**2))                   # the Gaussian speed profile (3.29)
wz = sp.diff(r*uth, r)/r                                            # (3.23) with u_r = 0: (1/r) d(r u_θ)/dr
print(sp.simplify(wz - Gam/(sp.pi*sig**2)*sp.exp(-r**2/sig**2)))    # 0: it gives back ω_z of (3.29)
print(sp.limit(2*sp.pi*r*uth, r, sp.oo))                            # Γ: all the circulation, far away (step 3)
""")
note("N43", "The Gaussian vortex", r"""
— the smooth version. It is the Lamb–Oseen vortex of viscous flow at one instant, with $\sigma^2=4\nu t$ growing in time
(Ch. 5, Ch. 8) — our pointer, not the book's.
""", equation=r"\omega_z(r)=\frac{\Gamma}{\pi\sigma^2}\exp\!\big(-r^2/\sigma^2\big)\quad\text{and}\quad u_\theta(r)=\frac{\Gamma}{2\pi r}\Big(1-\exp\!\big(-r^2/\sigma^2\big)\Big)", ref="3.29")
P("P107", "np.expm1 and cancellation near zero", r"""
Near r = 0, $1-e^{-x}$ subtracts two numbers that are almost equal and loses digits (*catastrophic cancellation*).
`np.expm1(y)` computes $e^y-1$ accurately for tiny y, so $1-e^{-x}=$ `-np.expm1(-x)`. `gaussian_vortex` uses it and
returns exactly 0 at r = 0.
""", code="""
x = 1e-12                                            # a tiny argument
print(1 - np.exp(-x), -np.expm1(-x))                 # 9.99977878e-13 (wrong from the 5th digit on) vs 9.999999999995e-13
""")
P("P108", "scipy.optimize.brentq", r"""
Finds a root of f(x) = 0 inside a bracket [a, b] where f changes sign; guaranteed to converge, and fast. Choose the
bracket so it excludes roots you do not want.
""", code="""
from scipy.optimize import brentq                    # bracketing root finder
g = lambda x: 1 + 2*x - np.exp(x)                    # the equation of D20
print(round(g(0.5), 4), round(g(3.0), 4))            # 0.3513 > 0, −13.0855 < 0: a sign change in (0.5, 3)
print(brentq(g, 0.5, 3.0))                           # 1.2564312086…
""")
D("D20")
note("N44", "Where the wind peaks", r"""
Rankine: at r = σ exactly, $u_{\max}=\Gamma/2\pi\sigma$. Gaussian: where the equation below holds, i.e.
$r\approx1.1209\,\sigma$, with $u_{\max}\approx0.6382\,\Gamma/(2\pi\sigma)$ — 36 % below the Rankine peak for the same Γ
and σ.
""", equation=r"1+2\,\frac{r^2}{\sigma^2}=\exp\!\big(r^2/\sigma^2\big)")
nb.worked_example("Γ = 2π m²/s, σ = 1 m", r"""
1. Rankine inside: $u_\theta=(2\pi/2\pi\cdot1)r=r$ m/s; $\omega_z=2\pi/\pi=2$ s⁻¹.
2. Rankine peak at r = 1 m: 1 m/s; at r = 2 m: $2\pi/(2\pi\cdot2)=0.5$ m/s.
3. Gaussian at r = 1 m: $u_\theta=(1/1)(1-e^{-1})=0.632$ m/s.
4. Gaussian peak: r = 1.1209 m, $x=r^2=1.2564$, $u_\theta=(1-e^{-1.2564})/1.1209=0.7153/1.1209=0.638$ m/s.
5. Far away both → 1/r: at 5 m, 0.200 m/s (Gaussian $(1-e^{-25})/5=0.2$ to 11 digits).
""")
nb.code("""
Gam, sig = 2*np.pi, 1.0                                      # circulation [m²/s] and core radius [m]
r = np.array([0.5, 1.0, 1.1209, 2.0, 5.0])                   # five radii [m]
print("Rankine  u, w:", [np.round(v, 4) for v in ch03.rankine_vortex(r, Gam, sig)])    # (3.28): u_θ [m/s], ω_z [1/s]
print("Gaussian u, w:", [np.round(v, 4) for v in ch03.gaussian_vortex(r, Gam, sig)])   # (3.29)
rstar = ch03.gaussian_vortex_max_radius(sig)                 # the root of 1 + 2x = e^x by brentq (D20)
print(rstar, ch03.gaussian_vortex(rstar, Gam, sig)[0])       # 1.1209064 m and the peak speed 0.63817 m/s
print(ch03.gaussian_vortex_max_radius(sig, method="lambertw"))   # the same radius from the Lambert-W formula
print(ch03.circulation_circle("gaussian", 3.0, Gamma=Gam, sigma=sig)/Gam)   # Γ(3σ)/Γ = 1 − e^−9
for name, s_, umax in (("bathtub", 5e-3, 0.3), ("tornado", 100.0, 60.0), ("tropical cyclone", 36e3, 50.0)):
    G_ = 2*np.pi*s_*umax/0.6382                              # Γ from the peak wind of a Gaussian vortex [m²/s]
    print(f"{name:17s} sigma = {s_:8.3g} m  u_max = {umax:4.1f} m/s  ->  Gamma = {G_:.2g} m^2/s, r_max = {1.1209*s_:.3g} m")
""", explain="""
1. Both models at five radii: Rankine (0.5, 1, 0.892, 0.5, 0.2) m/s with a vorticity step 2 → 0; Gaussian (0.442, 0.632,
   0.638, 0.491, 0.200) m/s with a smooth bell.
2. The maximum by `brentq` (D20) and its Lambert-W twin: 1.1209064 m.
3. $\\Gamma(r)/\\Gamma\\to1$ outside the core (D19 step 3): 0.99988 at 3σ.
4. Our order-of-magnitude real vortices span nine orders of magnitude in Γ with the same shape.
""")
nb.figure(r"""
x = np.linspace(1e-3, 4, 400)                                # r/σ [–]
uR, wR = ch03.rankine_vortex(x, 2*np.pi, 1.0)                # Γ = 2π, σ = 1: u is already in units of Γ/2πσ = 1 m/s; ω/2 in units of Γ/πσ² = 2 1/s
uG, wG = ch03.gaussian_vortex(x, 2*np.pi, 1.0)
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 3.9))           # (a) speed, (b) vorticity
a.plot(x, uR, color=COLORS["teal"], lw=2.3, label="Rankine (3.28)")          # kink at r = σ
a.plot(x, uG, "--", color=COLORS["teal"], lw=2.3, label="Gaussian (3.29)")   # smooth
a.plot(x, x, ":", color=COLORS["muted"], label="solid body"); a.plot(x, 1/x, ":", color=COLORS["grid"], label="line vortex")   # ghosts
a.plot([1.0], [1.0], "o", color=COLORS["accent"]); a.plot([1.1209], [0.6382], "o", color=COLORS["accent"], label="peaks: σ and 1.1209σ")
a.set_ylim(0, 1.3); a.set_xlabel(r"$r/\sigma$ [–]"); a.set_ylabel(r"$u_\theta\,/\,(\Gamma/2\pi\sigma)$ [–]"); a.legend(fontsize=8)
a.set_title("Speed: a spinning core and an irrotational skirt")
b.plot(x, wR/2, color=COLORS["orange"], lw=2.3, label="Rankine: a step")     # uniform core, nothing outside
b.plot(x, wG/2, "--", color=COLORS["orange"], lw=2.3, label="Gaussian: a bell")   # e^{−r²/σ²}
b.set_xlabel(r"$r/\sigma$ [–]"); b.set_ylabel(r"$\omega_z\,/\,(\Gamma/\pi\sigma^2)$ [–]"); b.legend(fontsize=8)   # dimensionless
b.set_title("Vorticity: where the core is")
savefig(fig, "ch03", "vortex_profiles"); plt.show()           # save to outputs/ch03 and draw
""", see="Two humps with peaks at 1 and ≈ 1.12 (left); a step and a bell (right).",
    read="""The core is where the vorticity is; outside it both speed profiles follow the line vortex $\\Gamma/2\\pi r$. The
smooth Gaussian vorticity moves the peak outward and lowers it (0.638 instead of 1 in units of Γ/2πσ).""",
    change="…σ were doubled at the same Γ: in physical units the peak would halve and move twice as far out ($u_{\\max}\\propto\\Gamma/\\sigma$).")
nb.plotly(r"""
rr = np.linspace(0.005, 5, 300)                              # radius [m]

def profiles(s):                                             # Γ = 2π m²/s, core radius s [m]
    rs = ch03.gaussian_vortex_max_radius(s)                  # Gaussian peak radius (D20)
    return {"Rankine u_θ": (rr, ch03.rankine_vortex(rr, 2*np.pi, s)[0]),        # (3.28)
            "Gaussian u_θ": (rr, ch03.gaussian_vortex(rr, 2*np.pi, s)[0]),      # (3.29)
            "peaks (σ and 1.1209σ)": ([s, rs], [1.0/s, ch03.gaussian_vortex(rs, 2*np.pi, s)[0]])}   # Rankine peak Γ/2πσ = 1/s

fig = slider_figure(profiles, "σ", np.linspace(0.25, 2, 15 if not FAST else 8), unit="m", xlabel="r [m]",
                    ylabel="u_θ [m/s]", xrange=[0, 5], yrange=[0, 4.2], modes={"peaks (σ and 1.1209σ)": "markers"},
                    title="A smaller core spins faster: u_max ∝ Γ/σ")    # fixed axes
recolor(fig, {"Rankine u_θ": COLORS["teal"], "Gaussian u_θ": COLORS["teal"], "peaks (σ and 1.1209σ)": COLORS["accent"]},
        dashes={"Gaussian u_θ": "dash"})                     # Gaussian dashed, as in the static figure
fig.show()                                                   # draw it
""", explain="The slider changes the core radius σ at fixed Γ = 2π m²/s; the markers sit at the two peaks.")
see_read_change(
    "A teal kinked curve (Rankine), a dashed smooth one (Gaussian) and two purple markers.",
    "The Rankine peak sits on the kink at σ, the Gaussian peak at 1.12σ, both moving outward and down as σ grows; beyond about 2σ the two curves coincide with Γ/2πr.",
    "…Γ doubled: every curve would double in height; the peak radii would not move.")
nb.explainer("vortex_paddle_wheels", "Going round in circles ≠ spinning", r"""
**Why interactive:** "going round" and "spinning" can only be told apart in motion; watching wheels orbit and turn (or not) on one clock, and dragging a loop on and off the axis, makes the difference visible.
Four vortices on one stage: paddle wheels ride with the flow and turn at half the local vorticity, a loop you can drag
measures the circulation, and the profiles $u_\theta(r)$ and $\omega_z(r)$ show where the core is. The same clock runs
the tracers, the wheels and the circulation graph.
""", tries=[
    "Mode 'line vortex': the wheels orbit but always point the same way; drag the loop off the axis — Γ drops to 0.",
    "Mode 'solid body': every wheel turns once per revolution.",
    "Mode 'Gaussian': drag σ and watch the peak marker sit at 1.1209σ.",
    "Preset 'tropical cyclone': read the radius of maximum wind and the core vorticity in Explain.",
])
whatif(r"""
…the fluid were viscous? The Gaussian core would spread as $\sigma^2=4\nu t$ while Γ stays fixed (Ch. 5). And all of
this has been about *points* — §3.6 asks how a quantity inside a whole moving volume changes.
""")

# =====================================================================================================================
# A.6 §3.6 Reynolds Transport Theorem — C15
# =====================================================================================================================
nb.section("3.6", "Reynolds Transport Theorem", intro="""
**What is this section about?** How fast does the amount of something inside a volume change when the volume itself
moves and deforms? In one dimension this is Leibniz's rule; in three it is the Reynolds transport theorem — the bridge
from "what happens to a moving lump of fluid" to equations at fixed points, used for every conservation law in Ch. 4.
""")
core("C15", "The Reynolds transport theorem: change inside a moving volume (3.35)", r"""
A balloon is being inflated in a room that is warming up. How fast does the heat content inside the balloon change —
and which part of that comes from the warming, which from the balloon's growing skin?
""", eqs=("3.35",))
nb.md(r"""
#### The problem in plain words

Conservation laws are about *things* — a lump of fluid keeps its mass, Newton's law acts on it — but the lump moves and
deforms. Engineers draw a *control volume* (a pipe section, a jet engine, a layer of ocean) whose walls may move. In
both cases we need d/dt of an integral whose region moves: the answer has an "inside" part and a "swept by the walls"
part.
""")
nb.md(r"""
#### The idea

```
d/dt ∫_{V*(t)} F dV  =  ∫_{V*} ∂F/∂t dV   +   ∮_{A*} F b·n dA           (3.35)
   change of the total    change in place       what the moving wall sweeps in (b·n > 0) or out (b·n < 0)
1-D:  d/dt ∫_a^b F dx = ∫_a^b ∂F/∂t dx + ḃ F(b) − ȧ F(a)                  (3.30)
```

In symbols: $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$ *(3.35)*, and in 1-D
$\frac{d}{dt}\int_{a}^{b}F\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)$ *(3.30)*.
""")
note("N45", "Why we need it", r"""
every integral conservation law of Ch. 4 (mass, momentum, energy) is a time derivative of an integral over a moving,
deforming volume.
""")
P("P109", "differentiation under the integral sign", r"""
If the limits are fixed, the time derivative may pass inside:
$\frac{d}{dt}\int_a^bF(x,t)\,dx=\int_a^b\frac{\partial F}{\partial t}dx$ (F and ∂F/∂t continuous). It is the special case
of Leibniz's rule with ȧ = ḃ = 0.
""", code="""
t = sp.symbols('t'); x = sp.symbols('x')             # time and position
F = sp.sin(t)*x**2 + t**2*x                          # any smooth F(x, t)
print(sp.simplify(sp.diff(sp.integrate(F, (x, 0, 1)), t) - sp.integrate(sp.diff(F, t), (x, 0, 1))))   # 0
""")
note("N46", "Leibniz's theorem", r"""
— d/dt of an integral with moving limits. The book cites a proof it does not give; here it is (D21).
""", equation=r"\frac{d}{dt}\int_{x=a(t)}^{x=b(t)}F(x,t)\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)", ref="3.30")
D("D21", ref="3.30", check_src="""
x, t = sp.symbols('x t', positive=True)             # position and time
F, a, b = x**2*t, t, t**2                           # the check case: F = x²t between a = t and b = t²
lhs = sp.diff(sp.integrate(F, (x, a, b)), t)        # left side: differentiate the integral itself
rhs = sp.integrate(sp.diff(F, t), (x, a, b)) + sp.diff(b, t)*F.subs(x, b) - sp.diff(a, t)*F.subs(x, a)   # right side of (3.30)
print(sp.simplify(lhs - rhs))                       # 0: Leibniz holds
print(lhs.subs(t, 2))                               # 416/3 = 138.667 at t = 2 s
""", extra_check="The sympy cell below runs the case $F=x^2t$, $a=t$, $b=t^2$.")
note("N47", "Fig. 3.17, the picture of (3.30)", r"""
three thin strips — a band of height $\frac{\partial F}{\partial t}dt$ over [a, b] (the interior change), a strip of
width db and height F(b) gained at the upper end, a strip of width da and height F(a) lost at the lower end; the corner
pieces are dt² small. Below: $F=1+0.5x+0.3t$ between $a(t)=1-0.2t$ and $b(t)=3+0.4t$, so
$\frac{d}{dt}\int_a^bF\,dx=\int_a^b\frac{\partial F}{\partial t}dx+\frac{db}{dt}F(b,t)-\frac{da}{dt}F(a,t)$.
""")
nb.animation(r"""
Fw, dFw = ch03.rtt_field("warming", dim=1)                    # F = 1 + 0.5 x + 0.3 t and ∂F/∂t = 0.3
a_ = lambda t: 1.0 - 0.2*t; b_ = lambda t: 3.0 + 0.4*t         # the moving limits [m]
nfr = 12 if not FAST else 8                                   # frames to step through
times = np.linspace(0, 2.2, nfr); dt = 0.4                    # frame times and a (visible) strip time Δt [s]
xx = np.linspace(0, 4.5, 300)
fig, ax = plt.subplots(figsize=(7.5, 4))                     # one panel

def update(i):                                                # redraw frame i from scratch (few artists)
    ax.clear(); t = times[i]; a, b = a_(t), b_(t)             # this frame's time and limits
    ax.plot(xx, Fw(xx, t), color=COLORS["ink"], lw=2, label="$F(x, t)$")    # the field now
    xi = np.linspace(a, b, 100)                               # points between the limits
    ax.fill_between(xi, 0, Fw(xi, t), color=COLORS["grid"], alpha=0.8)                     # ∫_a^b F dx now
    ax.fill_between(xi, Fw(xi, t), Fw(xi, t + dt), color=COLORS["blue"], alpha=0.5, label="interior: Δt ∂F/∂t")   # F rises in place
    xb = np.linspace(b, b_(t + dt), 20); ax.fill_between(xb, 0, Fw(xb, t), color=COLORS["orange"], alpha=0.6, label="gained at b: ḃΔt F(b)")
    xa = np.linspace(a_(t + dt), a, 20); ax.fill_between(xa, 0, Fw(xa, t), color=COLORS["rose"], alpha=0.6, label="a moves left: −ȧΔt F(a) > 0")
    L = ch03.leibniz_terms(Fw, dFw, a, b, -0.2, 0.4, t)       # the three rates of (3.30) at this t
    ax.set_title(f"t = {t:.2f} s: interior {L.interior:.2f} + upper {L.upper:.2f} − lower ({L.lower:.2f}) = {L.total:.2f} [F·m/s]", fontsize=9)
    ax.set_xlim(0, 4.5); ax.set_ylim(0, 4.0); ax.set_xlabel("$x$ [m]"); ax.set_ylabel("$F$"); ax.legend(fontsize=7, loc="upper left")
    return []                                                 # everything was redrawn

show_animation(animate(update, frames=nfr, fig=fig, interval=600), player="frames")   # step through the strips
""", explain="""
1. `rtt_field("warming", dim=1)` is the shared test field of the notebook and the explainer.
2. Each frame shades the three strips of
   $\\frac{d}{dt}\\int_a^bF\\,dx=\\int_a^b\\frac{\\partial F}{\\partial t}dx+\\frac{db}{dt}F(b,t)-\\frac{da}{dt}F(a,t)$ *(3.30)*
   for a visible step Δt = 0.4 s, and prints the exact rates from
   `leibniz_terms` (the lower term is $\\dot aF(a)$, which is *subtracted*).
""")
see_read_change(
    "A grey area under a rising line, a blue band on top of it, an orange strip on the right and a rose strip on the left.",
    """The orange strip grows on the right as b advances; the rose strip is *added* on the left here because a moves left
(ȧ < 0 makes $-\\dot aF(a)$ positive). At t = 0 the three rates are 0.6, 1.0 and 0.3, total 1.9 (checked in the code
below).""",
    "…a moved right instead (ȧ > 0): the rose strip would be taken away and the total would drop.")
note("N48", "Control volume", r"""
$V^*(t)$ with **control surface** $A^*(t)$, outward unit normal $\mathbf n$ and surface velocity $\mathbf b$ (Fig.
3.18). The surface need not follow the fluid: $\mathbf b=\mathbf u$ for a *material* volume, $\mathbf b=0$ for a volume
fixed in space, anything else for a piston, a balloon or a moving layer. `ch03.GrowingSphere`, `GrowingCylinder`,
`MovingBox`, `GrowingCone` are concrete examples (Python objects that know their volume, surface and b at any t).
""")
P("P110", "signed swept volume of a moving surface", r"""
In a short time Δt a surface patch dA moving with velocity $\mathbf b$ sweeps a thin prism of height
$(\mathbf b\Delta t)\cdot\mathbf n$ — only the normal component counts (sliding along the surface sweeps nothing). With
the outward $\mathbf n$ the volume is **signed**: positive where the wall advances outward (the region gains volume),
negative where it retreats.
""", code="""
n = np.array([1.0, 0.0]); dA, dt = 0.01, 0.1         # outward normal, patch area [m²], time step [s]
for b in ([2.0, 0.0], [0.0, 5.0], [-1.0, 3.0]):      # advancing, sliding, retreating wall velocities [m/s]
    print(b, np.dot(b, n)*dt*dA)                     # 0.002, 0.0 (sliding), -0.001 (retreating) [m³]
""")
note("N49", "The start of the derivation", r"""
the definition of a time derivative (Fig. 3.18: the volume solid at t, dashed at t + Δt).
""", equation=r"\frac{d}{dt}\int_{V^*(t)}F(\mathbf x,t)dV=\lim_{\Delta t\to0}\frac{1}{\Delta t}\Big\{\int_{V^*(t+\Delta t)}F(\mathbf x,t+\Delta t)dV-\int_{V^*(t)}F(\mathbf x,t)dV\Big\}", ref="3.31")
note("N50", "Four terms", r"""
Splitting the new volume, $\Delta V\equiv V^*(t+\Delta t)-V^*(t)$, and Taylor-expanding F in time gives:
""", equation=r"\int_{V^*(t+\Delta t)}F(\mathbf x,t+\Delta t)dV\cong\int_{V^*(t)}F\,dV+\int_{V^*(t)}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV+\int_{\Delta V}\Delta t\frac{\partial F}{\partial t}dV", ref="3.32")
note("N51", "The first cancels, the last is second order", r"""
""", equation=r"\frac{d}{dt}\int_{V^*(t)}F\,dV=\lim_{\Delta t\to0}\frac{1}{\Delta t}\Big\{\int_{V^*(t)}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV\Big\}", ref="3.33")
note("N52", "…and the sliver is a surface integral", r"""
Each of these is a step of D22 below, where the moves between them are filled in.
""", equation=r"\int_{\Delta V}F(\mathbf x,t)dV\cong\int_{A^*(t)}F(\mathbf x,t)(\mathbf b\Delta t\cdot\mathbf n)dA\quad\text{as}\quad\Delta t\to0", ref="3.34")
D("D22", ref="3.35")
note("N53", "Two readings of (3.35)", r"""
of $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$.
(1) With F = 1 it says the volume changes at the rate the surface sweeps: $dV^*/dt=\int_{A^*}\mathbf b\cdot\mathbf n\,dA$;
for a small material volume (b = u) this is $\frac{1}{\delta V}\frac{D}{Dt}(\delta V)=\frac{\partial u_i}{\partial x_i}$
*(3.14)* — D23. (2) For a small material volume it contains the material derivative
$\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ *(3.5)*, plus a term $F\nabla\cdot\mathbf u$ from the
volume change — D24. Where b·n > 0 the surface advances, where b·n < 0 it retreats; one signed term covers both. Only for
a fixed volume (b = 0) may d/dt pass inside the integral.
""")
D("D23", ref="3.14")
D("D24", ref="3.5")
nb.worked_example("a balloon growing in a warming room", r"""
Sphere of radius R = 1 m growing at Ṙ = 0.1 m/s; F = t (a "heat content per volume" rising uniformly, ∂F/∂t = 1 per
second) evaluated at t = 1 s.

1. Volume term: $\int_{V^*}\partial F/\partial t\,dV=1\times\tfrac43\pi R^3=4.189$.
2. Surface term: on the sphere $\mathbf b\cdot\mathbf n=\dot R=0.1$ m/s and F = 1: $\oint F\,\mathbf b\cdot\mathbf n\,dA=1\times0.1\times4\pi R^2=1.257$.
3. Total 5.445.
4. Direct: $\int F\,dV=t\cdot\tfrac43\pi R(t)^3$, whose derivative is $\tfrac43\pi R^3+t\cdot4\pi R^2\dot R=4.189+1.257=5.445$ ✓.
""")
nb.code("""
cv = ch03.GrowingSphere(0.9, 0.1)                     # R(t) = 0.9 + 0.1 t [m], so R = 1.0 m at t = 1 s (b·n = Ṙ on the surface)
F = lambda x, t: t + 0*x[0]                           # F(x, t) = t: uniform, rising at 1 per second
dF = lambda x, t: 1.0 + 0*x[0]                        # ∂F/∂t
print(ch03.reynolds_transport(F, dF, cv, 1.0))        # (3.35): volume term, surface term, total
print(ch03.rtt_check(F, dF, cv, 1.0))                 # finite difference of ∫F dV (left of (3.31)) vs the RTT total
Fw, dFw = ch03.rtt_field("warming", dim=1)            # 1-D Leibniz with the explainer's default numbers
print(ch03.leibniz_terms(Fw, dFw, 1.0, 3.0, -0.2, 0.4, 0.0))   # a = 1, b = 3, ȧ = −0.2, ḃ = 0.4 at t = 0
d = ch03.leibniz_example(2.0)                         # the D21 check case F = x²t, a = t, b = t² at t = 2 s
print({k: round(d[k], 3) for k in ("interior", "upper", "lower", "total", "exact")})
F2, dF2 = ch03.rtt_field("ramp", dim=2)               # F = 1 + 0.5 x1 in the plane
print(ch03.rtt_ellipse_2d(F2, dF2, 2.0, 1.0, 0.2, 0.1, (0.0, 0.0), (0.3, 0.0), 0.0))   # a translating, growing ellipse
print(ch03.material_volume_rate(lambda x, t: 0.1*x, ch03.GrowingSphere(1.0, 0.1), 0.0))   # D23: ∮u·n dA = ∫∇·u dV for u = 0.1 x
""", explain="""
1. `GrowingSphere` is a control volume that knows its nodes, its surface normal and its surface velocity b at any t.
2. `reynolds_transport` evaluates the two terms of
   $\\frac{d}{dt}\\int_{V^*}F\\,dV=\\int_{V^*}\\frac{\\partial F}{\\partial t}dV+\\int_{A^*}F\\,\\mathbf b\\cdot\\mathbf n\\,dA$ *(3.35)*
   by quadrature: (4.1888, 1.2566, 5.4454).
3. `rtt_check` gets the same total from a finite difference of the volume integral itself.
4. Leibniz with the transport explainer's default numbers: 0.6 + 1.0 − (−0.3) = 1.9 (the lower term is $\\dot aF(a)$ and is subtracted).
5. The D21 check case: 18.667 + 128 − 8 = 138.667 = the exact derivative.
6. A translating, growing ellipse in a ramp field (the explainer's 2-D mode): $\\pi(\\dot ab+a\\dot b)(1+0.5c_x)+\\pi ab\\cdot0.5\\dot c_x=1.2566+0.9425$.
7. D23 on a sphere: for $\\mathbf u=0.1\\,\\mathbf x$ the outward flux $\\oint\\mathbf u\\cdot\\mathbf n\\,dA$ and $\\int\\nabla\\cdot\\mathbf u\\,dV$ are both $0.3\\times\\tfrac43\\pi=1.2566$ m³/s.
""")
nb.check_agree("""
nr, nth, nph = (16, 16, 32) if not FAST else (12, 12, 24)    # midpoint grid in r, θ, φ (Ch. 2 P83)

def sphere_nodes(R):                                  # midpoint nodes and weights of a ball of radius R
    r = (np.arange(nr) + 0.5)*R/nr; th = (np.arange(nth) + 0.5)*np.pi/nth; ph = (np.arange(nph) + 0.5)*2*np.pi/nph
    Rr, Th, Ph = np.meshgrid(r, th, ph, indexing="ij")
    w = (Rr**2*np.sin(Th)*(R/nr)*(np.pi/nth)*(2*np.pi/nph)).ravel()          # r² sin θ Δr Δθ Δφ
    X = np.stack([Rr*np.sin(Th)*np.cos(Ph), Rr*np.sin(Th)*np.sin(Ph), Rr*np.cos(Th)]).reshape(3, -1)
    return X, w

R, Rd, t0 = 1.0, 0.1, 1.0                             # radius at t0, growth rate, time
X, w = sphere_nodes(R)
vol = np.sum(dF(X, t0)*w)                             # ∫ ∂F/∂t dV as a midpoint sum
Th, Ph = np.meshgrid((np.arange(nth) + 0.5)*np.pi/nth, (np.arange(nph) + 0.5)*2*np.pi/nph, indexing="ij")
wA = (R**2*np.sin(Th)*(np.pi/nth)*(2*np.pi/nph)).ravel()                       # R² sin θ Δθ Δφ
Xs = R*np.stack([np.sin(Th)*np.cos(Ph), np.sin(Th)*np.sin(Ph), np.cos(Th)]).reshape(3, -1)
surf = np.sum(F(Xs, t0)*Rd*wA)                        # ∮ F b·n dA with b·n = Ṙ
assert np.allclose([vol, surf], ch03.reynolds_transport(F, dF, cv, 1.0)[:2], rtol=5e-3)   # midpoint ≈ library
I = lambda tt: np.sum(F(sphere_nodes(0.9 + 0.1*tt)[0], tt)*sphere_nodes(0.9 + 0.1*tt)[1])  # ∫_{V*(t)} F dV
dtt = 1e-3
fd = (I(t0 + dtt) - I(t0 - dtt))/(2*dtt)              # d/dt of the integral, by a central difference
assert np.allclose(fd, vol + surf, rtol=5e-3)         # = volume term + surface term
print(f"midpoint sums: volume {vol:.4f} + surface {surf:.4f} = {vol + surf:.4f};  d/dt of the sum {fd:.4f}\\n", end="")   # one write: no stray empty output
""")
nb.md("Midpoint sums on a spherical grid converge like $h^2$; with 16 × 16 × 32 cells they agree with the library to about 10⁻³.")
nb.plotly(r"""
F3 = lambda x, t: x[0]**2*t                           # a field that varies in space: F = x1² t
dF3 = lambda x, t: x[0]**2                            # ∂F/∂t
dts = np.logspace(-4, -1, 16 if not FAST else 10)     # time steps Δt [s]
rows = {d: ch03.swept_terms_sphere(1.0, 0.1, F3, dF3, 1.0, d) for d in dts}   # the four terms of (3.32), computed once
fig = slider_figure(lambda d: {"terms T1…T4": ([1, 2, 3, 4], [abs(rows[d][k]) for k in ("T1", "T2", "T3", "T4")])},
                    "Δt", dts, unit="s", xlabel="term of (3.32)", ylabel="size [F·m³]", xrange=[0.5, 4.5],
                    modes={"terms T1…T4": "lines+markers"}, height=420,
                    title="Eq. (3.32) for a growing sphere: T4 falls twice as fast as T2 and T3")   # one trace per step
for tr, d in zip(fig.data, dts):                      # one trace per slider step: put T4/T3 into its name
    tr.name = f"T4/T3 = {rows[d]['T4']/rows[d]['T3']:.1e}"; tr.marker.size = 12; tr.line.color = COLORS["blue"]
fig.update_yaxes(type="log", range=[-11, 1])          # log scale, fixed (10⁻¹¹ … 10)
fig.update_xaxes(tickvals=[1, 2, 3, 4], ticktext=["T1 ∫F", "T2 ∫Δt ∂F/∂t", "T3 ∫_ΔV F", "T4 ∫_ΔV Δt ∂F/∂t"])   # term names
fig.update_layout(showlegend=True)                    # the legend carries T4/T3
fig.show()                                            # draw it
""", explain=r"""
1. `swept_terms_sphere` evaluates the four integrals of
   $\int_{V^*(t+\Delta t)}F(t+\Delta t)\,dV\cong\int_{V^*}F\,dV+\int_{V^*}\Delta t\frac{\partial F}{\partial t}dV+\int_{\Delta V}F\,dV+\int_{\Delta V}\Delta t\frac{\partial F}{\partial t}dV$ *(3.32)* for a sphere of radius 1.1 m at t = 1 s growing at
   0.1 m/s, for 16 time steps between 10⁻⁴ and 10⁻¹ s.
2. The slider moves Δt; the legend shows the ratio T4/T3, which shrinks in proportion to Δt.
""")
nb.figure(r"""
T2 = np.array([rows[d]["T2"] for d in dts]); T3 = np.array([rows[d]["T3"] for d in dts])   # the kept terms [F·m³]
T4 = np.array([rows[d]["T4"] for d in dts]); sl = np.abs([rows[d]["sliver_error"] for d in dts])   # the dropped pieces
fig, ax = plt.subplots(figsize=(7, 4))                # log–log axes
for y, col, lab in ((T2, COLORS["blue"], "T2 (kept)"), (T3, COLORS["orange"], "T3 (kept)"),
                    (T4, COLORS["rose"], "T4 (dropped)"), (sl, COLORS["accent"], "|T3 − surface sliver| (3.34)")):
    ax.loglog(dts, y, "o-", color=col, ms=4, label=f"{lab}: slope {observed_order(dts, y):.2f}")   # slope by least squares
ax.set_xlabel(r"$\Delta t$ [s]"); ax.set_ylabel("size [F·m³]"); ax.legend(fontsize=8)            # axes with units
ax.set_title("Orders of smallness: first-order terms survive division by Δt, second-order ones do not")
savefig(fig, "ch03", "rtt_orders"); plt.show()        # save to outputs/ch03 and draw
""", see="Four straight lines on log–log axes: two of slope 1 (blue, orange), two of slope 2 (rose, purple).",
    read="""After dividing by Δt, the slope-1 terms tend to finite rates (the volume and surface terms), while the slope-2
terms — T4 and the error of replacing the sliver by a surface integral — tend to zero. That is exactly why (3.33),
$\\frac{d}{dt}\\int_{V^*}F\\,dV=\\lim\\frac{1}{\\Delta t}\\{\\int_{V^*}\\Delta t\\,\\partial F/\\partial t\\,dV+\\int_{\\Delta V}F\\,dV\\}$,
keeps two terms and drops the rest.""",
    change="…the sphere grew ten times faster: T3 and T4 would both grow ×10 (both ∝ ΔV), but T4/T3 would still be ∝ Δt.")
note("N54", "Ex. 3.2 — a growing cone", r"""
A right circular cone of fixed height h has base radius r(t) growing at ṙ. Directly, $V=\tfrac13\pi hr^2$ gives
$dV/dt=\tfrac23\pi hr_o\dot r$. By (3.35) with F = 1 and the cone as V*, only the sloping side contributes: a side point
at height z moves outward at $(z/h)\dot r\,\mathbf e_R$, the outward normal is $\mathbf n=\mathbf e_R\cos\theta-\mathbf e_z\sin\theta$
(θ = the cone's half-angle, $\tan\theta=r_o/h$ — a fourth θ in this chapter), and the slanted area element is
$dA=z\tan\theta\,d\varphi\,dz/\cos\theta$. ⚠️ Two print slips taught corrected: the base's points move radially in its own
plane, so it is **b·n = 0** there (not b = 0), which is all the theorem needs; the "[?]" in the printed integrand is just
the factor z.
""", equation=r"\frac{dV}{dt}=\int_{z=0}^{h}\int_{\varphi=0}^{2\pi}\frac zh\dot r\,\mathbf e_R\cdot(\mathbf e_R\cos\theta-\mathbf e_z\sin\theta)\,z\tan\theta\,d\varphi\frac{dz}{\cos\theta}=\frac{2\pi\dot r\tan\theta}{h}\int_0^hz^2dz=\tfrac23\pi h^2\dot r\tan\theta=\tfrac23\pi hr_o\dot r")
nb.code("""
print(ch03.example_3_2(1.0, 0.5, 0.1))               # h = 1 m, r_o = 0.5 m, ṙ = 0.1 m/s: dV/dt three ways [m³/s]
""", explain="The direct derivative, the book's surface integral by `dblquad` (P87), the general quadrature on `GrowingCone` and the closed form all give 0.10472 m³/s.")
note("N55", "Fig. 3.18, our 2-D drawing", r"""
a deforming, translating ellipse over the ramp field F = 1 + 0.5x₁.
""")
nb.figure(r"""
from scripts.ch03_drawings import rtt_blob_figure             # drawing helper; the numbers come from ch03.rtt_ellipse_2d
fig, terms, fd = rtt_blob_figure(a=2.0, b=1.0, adot=0.2, bdot=-0.1, cdot=(0.3, 0.0), F_name="ramp", dt=0.3)   # ellipse 2 × 1 m
fig.axes[0].set_title("swept band: blue where b·n > 0, rose where b·n < 0")   # a shorter title that fits the panel
print(f"volume term {terms[0]:.4f} + surface term {terms[1]:.4f} = {terms[2]:.4f};  measured d/dt of the integral {fd:.4f}")   # [F·m²/s]
savefig(fig, "ch03", "fig3_18_rtt"); plt.show()               # save to outputs/ch03 and draw
""", see="""Left: an ellipse (solid) and its position a moment later (dashed) over a grey ramp, orange b arrows on the
boundary, blue slivers on the right where the boundary advances and rose slivers at the top and bottom where it
retreats. Right: an orange surface bar and a purple total bar, both reaching the black diamond; the blue volume bar
has zero height.""",
    read="""The rose slivers are *negative* volume — the same formula counts them with the right sign. Here F does not
depend on time, so the volume term is zero and the whole change comes from the moving surface; the sum of the bars is
the measured rate.""",
    change="…b = 0 (a fixed box): the slivers and the orange bar would vanish, and d/dt could pass inside the integral.")
nb.explainer("reynolds_transport_cv", "What changes inside a moving box?", r"""
**Why interactive:** the theorem is a budget whose terms change sign with the direction the walls move; moving the walls yourself and watching the slivers turn blue or rose shows why one signed surface term covers all cases.
A moving boundary sweeps a band whose colour shows the sign of b·n; bars for the volume and surface terms of
$\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$ *(Eq. 3.35)*
add up to the measured rate of change of ∫F, and a Δt slider shows the dropped second-order term vanish. Modes: Leibniz
in 1-D, a deforming ellipse in 2-D, and the growing cone of Ex. 3.2.
""", tries=[
    "Preset 'fixed volume (b = 0)': the orange bar disappears.",
    "Preset 'rigid translation in a uniform F': the two terms cancel exactly — a uniform F carried along does not change.",
    "Make one side retreat: its sliver turns rose and subtracts.",
    "Derivation tab, D22: at the 'drop second order' step, scrub Δt on the log–log view.",
])
whatif(r"""
…F were the density ρ and V* a material volume (b = u)? The left side is the rate of change of the lump's mass, which
is zero; D24's steps then give $\frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf u=0$ — the continuity equation that opens
Chapter 4. With F = ρu you get Newton's law for the lump; with F = ρ(e + u²/2), the energy equation.
""")

# =====================================================================================================================
# A.7 End matter
# =====================================================================================================================
nb.pointer("""
S01 · Practise on the book's Exercises 3.1–3.30 (curvilinear operators, flow lines, the Galilean proof, strain and
rotation of simple flows, circulation, transport-theorem checks). The derivations the text leaves to Exercises 3.3,
3.12, 3.17–3.20, 3.23, 3.26, 3.28 and 3.30 are written out above in our own words (D03, D06, D10, D11, D13, D17, N37,
D20, D23, D24).
""")
nb.pointer("""
S02 · Further reading: Van Dyke's *An Album of Fluid Motion* and Samimy et al.'s *A Gallery of Fluid Motion* for
photographs of streak and path lines; Riley, Hobson & Bence for Leibniz's rule; Thompson (1972) for the geometric
transport-theorem argument (all cited in the book's literature list).
""")
nb.summary(
    clicked=[
        r"**C01** A particle is a label; a field is a function of (x, t); $F[\mathbf r(t;\mathbf r_o,t_o),t]=F(\mathbf x,t)$ at $\mathbf x=\mathbf r$ *(3.2)* translates between them.",
        r"**C02** $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ *(3.5)*: what a moving parcel feels = what a fixed probe sees + what it meets by moving.",
        r"**C03** A streamline is tangent to $\mathbf u$ at one frozen instant: $dx/u=dy/v=dz/w$ *(3.7)*.",
        r"**C04** A path line follows one particle, $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ *(3.8)*; a streak line collects everything from one port; in unsteady flow the three lines differ.",
        r"**C05** The acceleration $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u$ is the same for every observer moving at constant velocity *(3.9)*; only its split is not.",
        r"**C06** $du_i=(\partial u_i/\partial x_j)\,dx_j$ *(3.10)*: the neighbourhood of a point moves linearly.",
        r"**C07** The diagonal of S is the stretching rate per length; along any direction $\mathbf n$ it is $\mathbf n\cdot\mathbf S\cdot\mathbf n$.",
        r"**C08** The off-diagonal $S_{12}$ is half the closing rate of a right angle; rigid motion has S = 0.",
        r"**C09** $\frac1{\delta V}\frac{D(\delta V)}{Dt}=\frac{\partial u_i}{\partial x_i}=S_{ii}$ *(3.14)*: divergence is the swelling rate.",
        r"**C10** A fluid element spins at ½ω — the average of any perpendicular pair of threads — and the spin depends on the observer's rotation ($\omega'_z=\omega_z-2\Omega$).",
        r"**C11** $du_i=S_{ij}dx_j+\tfrac12(\boldsymbol\omega\times d\mathbf x)_i$ *(3.19)*: deformation plus rigid rotation.",
        r"**C12** S stretches a small sphere into an ellipsoid on its principal axes, $d\bar u_\alpha=\bar S_{\alpha\alpha}d\bar x_\alpha$ *(3.21)*, at the first instant.",
        r"**C13** $\omega_z=\frac1r\frac{\partial}{\partial r}(ru_\theta)-\frac1r\frac{\partial u_r}{\partial\theta}$ *(3.23)*: solid-body rotation spins every element (2ω₀); the line vortex spins none (except its axis).",
        r"**C14** Real vortices have a solid-body core and an irrotational outside; the peak wind is at σ (Rankine, $u_\theta=\Gamma r/2\pi\sigma^2$ inside, *(3.28)*) or 1.1209σ (Gaussian, $u_\theta=\frac{\Gamma}{2\pi r}(1-e^{-r^2/\sigma^2})$ *(3.29)*).",
        r"**C15** $\frac{d}{dt}\int_{V^*}F\,dV=\int_{V^*}\frac{\partial F}{\partial t}dV+\int_{A^*}F\,\mathbf b\cdot\mathbf n\,dA$ *(3.35)*: change in place + what the moving walls sweep (signed).",
    ],
    feeds_forward=[
        r"Ch. 4: the transport theorem *(3.35)* + $\frac{DF}{Dt}=\frac{\partial F}{\partial t}+\mathbf u\cdot\nabla F$ *(3.5)* → continuity, Cauchy's equation, energy; S → the Newtonian stress law; rotating frames ($\omega'=\omega-2\Omega$) → Coriolis.",
        r"Ch. 5: vorticity, circulation $\Gamma=\oint\mathbf u\cdot d\mathbf s$ *(3.18)*, vortex stretching (S along ω), the Lamb–Oseen vortex.",
        r"Ch. 6: irrotational flow $\boldsymbol\omega=0$ *(3.17)*, the point vortex $u_\theta=B/r$ *(3.25)*, the cylinder of Fig. 3.2.",
        r"Ch. 7: particle orbits via $d\mathbf r/dt=\mathbf u(\mathbf r,t)$ *(3.8)*; wave frames via Galilean invariance, $\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=\frac{\partial\mathbf u'}{\partial t'}+(\mathbf u'\cdot\nabla')\mathbf u'$ *(3.9)*.",
        "Ch. 13: absolute vs relative vorticity, frontogenesis on principal strain axes, cyclone profiles.",
    ],
    left_out=[
        "Curvilinear forms of ∇, ∇², (u·∇)u (Appendix B; used from Ch. 4 on).",
        "Solutions for the velocity potential (Ch. 6).",
        "Viscous spreading of the Gaussian vortex (Ch. 5, 8).",
    ],
)

if __name__ == "__main__":
    path = nb.save()
    print(f"wrote {path.relative_to(ROOT)}: {len(nb.cells)} cells, {len(nb.cores)} CORE, {len(nb.recaps)} RECAP, "
          f"{len(nb.derivations)} derivations, {len(nb.primers)} primers, {len(nb.explainers)} explainers")
