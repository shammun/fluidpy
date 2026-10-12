"""Build the Chapter 14 teaching notebook: ``notebooks/ch14_aerodynamics.ipynb``.

Source of truth: ``analysis/ch14_design.md`` Part A (storyboard, one nbkit call per row), Part C (the function contract —
``fluidpy.ch14_aerodynamics`` imported as ``ch14`` re-exports ``core.airfoils`` = ``AF`` and ``core.lifting_line`` = ``LL``),
Part D (runtime budget and FAST sizes), Part E (prerequisite ledger → primers P336–P350 and one-line reminders of earlier
primers), Part F (the 22 derivations D01–D22, one move per step) and ``analysis/ch14_curation.md`` (IDs, depths, section
coverage), with these rulings written into the text:

* Γ is the CLOCKWISE circulation everywhere in this chapter (the book's sign); every readout says so;
* angles are computed in radians and shown in degrees, both on screen; lift-curve slopes are per radian (``K=`` keyword);
* slips (false as printed, ``ch14.book_slips()``) are taught in corrected form with a ``slip #k`` box; traps (true as printed,
  easy to misread, ``ch14.traps()``) are ``Common confusion`` callouts; the two are kept apart;
* results that are ours and not in the book are labelled "Not in the book — ours"; every figure made with a qualitative
  model says "Qualitative model" in its title;
* the Wagner ghost curve is drawn: ``reports/ch14_verification.md`` lists its source as read first-hand (a published fit);
* no book number is reproduced (rule 9): worked numbers use the design's own input set, written literally in the imports
  cell (``sec``, ``ac``, ``wd``, ``boat``, ``fish``); ``ch14.illustrative_inputs()`` is shown once beside it.

**Derivations are read from Part F at build time** (``part_f()`` below): goal, start, plan, tools, assumptions, every step's
*did / tex / why / plain*, result, check, meaning and traps are copied word for word; the four ★★★ derivations carry the
sympy check cells listed in Part F.

**Equations.** ``EQ`` below is the analyst's page-image transcription (``analysis/ch14.md`` §2; the designer re-read the 12
pages that carry all 27 numbered equations; the builder re-read p707 and p714 = (14.9)–(14.12), (14.18)–(14.24)), in the
CORRECTED form where the book prints a slip.

Run:  .venv/Scripts/python.exe notebooks/build_ch14.py            (writes the notebook)
      .venv/Scripts/python.exe notebooks/build_ch14.py --dump     (prints the parsed Part F derivations only)
"""
from __future__ import annotations

import pathlib
import re
import sys
import textwrap

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from nbkit import ChapterNotebook, curation_items  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

nb = ChapterNotebook("ch14")
DESIGN = (ROOT / "analysis" / "ch14_design.md").read_text(encoding="utf-8")
_CTRL = sorted({hex(ord(_c)) for _c in DESIGN if ord(_c) < 32 and _c != "\n"})
if _CTRL:                                # a form feed, tab or backspace in the design = a LaTeX backslash eaten by a non-raw string
    raise SystemExit(f"analysis/ch14_design.md contains control characters {_CTRL} — repair the eaten backslash before building")

# ---------------------------------------------------------------------------------------------------------------------
# Equations used in prose: every mention of a book equation writes the equation itself next to its number.
# ---------------------------------------------------------------------------------------------------------------------
EQ: dict[str, str] = {
    "14.1": r"\Lambda\equiv s^2/A=s/\bar c",
    "14.2": r"\bar c=A/s",
    "14.3": r"z=\zeta+\dfrac{b^2}{\zeta}",
    "14.4": r"x=(R+b^2/R)\cos\theta,\ \ y=(R-b^2/R)\sin\theta",
    "14.5": r"x^2\sin^2\theta-y^2\cos^2\theta=4b^2\sin^2\theta\cos^2\theta",
    "14.6": r"2b\tan\beta\,\sin\theta=R-b^2/R=y/\sin\theta",
    "14.7": r"a\equiv b(1+e),\ e\ll1",
    "14.8": r"c\simeq4b,\ \ \text{camber}\simeq\tfrac12\beta c,\ \ t_{\max}/c\simeq1.3\,e",
    "14.9": r"\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}",
    "14.10": r"u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}",
    "14.11": r"\Gamma=4\pi Ua\sin(\alpha+\beta)",
    "14.12": r"C_L=\dfrac{L}{(1/2)\rho U^2c}\simeq2\pi(\alpha+\beta)",
    "14.13": r"w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\,\dfrac{dy}{(y_1-y)}",
    "14.14": r"\alpha_e=\alpha-\varepsilon=\alpha-\dfrac wU",
    "14.15": r"L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy,\ \ D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy",
    "14.16": r"C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]",
    "14.17": r"\Gamma(y)=\dfrac K2\,U\,c(y)\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]",
    "14.18": r"\Gamma=\displaystyle\sum_{n=1}^{\infty}\Gamma_n\sin(n\gamma)",
    "14.19": r"w(y_1)=\dfrac{1}{2s}\displaystyle\sum_{n=1}^{\infty}n\,\Gamma_n\,\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}",
    "14.20": r"\dfrac K2\,Uc(\alpha+\beta)=\displaystyle\sum_{n=1}^{\infty}\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)",
    "14.21": r"L=\dfrac{\pi s}{4}\,\rho U\Gamma_1",
    "14.22": r"D_i=\dfrac\pi8\,\rho\displaystyle\sum_{n=1}^{N}n\,\Gamma_n^2",
    "14.23": r"\Gamma=\Gamma_1\sin(\gamma)=\Gamma_1\sqrt{1-(2y/s)^2}",
    "14.24": r"w(y)=\Gamma_1/2s",
    "14.25": r"D_i=\dfrac\pi8\rho\Gamma_1^2=\dfrac{2L^2}{\pi\rho U^2s^2}",
    "14.26": r"C_{D_i}=\dfrac{D_i}{(1/2)\rho U^2A}=\dfrac{C_L^2}{\pi(s^2/A)}=\dfrac{C_L^2}{\pi\Lambda}",
    "14.27": r"dL=\rho U\Gamma\,dy=C_L\,\tfrac12\rho U^2c\,dy",
    # earlier chapters (written-out forms of their designs)
    "4.107": r"C_D\equiv\dfrac{D}{\tfrac12\rho U^2A}",
    "4.108": r"C_L\equiv\dfrac{L}{\tfrac12\rho U^2A}",
    "5.8": r"\dfrac{D\Gamma}{Dt}=0",
    "5.17": r"d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}",
    "6.37": r"u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}",
    "6.38": r"\sin\theta=-\dfrac{\Gamma}{4\pi aU}",
    "6.60": r"D-iL=\dfrac{i\rho}{2}\displaystyle\oint\Big(\dfrac{dw}{dz}\Big)^2dz",
    "6.62": r"L=\rho U\Gamma",
    "6.65": r"z=\zeta+\dfrac{b^2}{\zeta}",
    "6.68": r"w=U\Big(\zeta+\dfrac{a^2}{\zeta}\Big)+\dfrac{i\Gamma}{2\pi}\ln(\zeta/a)",
}
for _k, _v in EQ.items():
    assert "$" not in _v, (_k, _v)


def E(num: str) -> str:
    """Inline "equation (number)" for prose: the equation is always written next to its number."""
    return f"${EQ[num]}$ ({num})"


def _join(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(x.strip() for x in lines)).strip()


def _abs_plain(s: str) -> str:
    r"""Part F writes absolute values as \lvert … \rvert everywhere; outside $…$ they become plain |…|."""
    parts = re.split(r"(\$[^$]+\$)", s)
    return "".join(q if q.startswith("$") else re.sub(r"\\rvert", "|", re.sub(r"\\lvert\s?", "|", q)) for q in parts)


def _tidy(s: str) -> str:
    """Plain-text exponents of Part F ("e^{−iθ}", "(2b − x)^{3/2}") outside maths and code become readable markdown."""
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
# Part F of the design, parsed at build time: every derivation copied word for word (the ch07–ch13 pattern; the design
# writes the fields inline — "**Goal:** … · **Start:** … · **Plan:** • … • …" — and the steps as
# "n. **did** … · **tex** … · **why** … · **plain** …").
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
        body = re.sub(r"```python\n.*?```", "", body, flags=re.S)
        body = re.sub(r"^- \*\*Sympy check\*\*.*$", "", body, flags=re.M | re.I)      # the line that introduces the listing
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
assert len(PF) == 22, len(PF)
N_STEPS = sum(len(d["steps"]) for d in PF.values())
assert N_STEPS == 169, N_STEPS                              # the 169 steps of Part F, none lost by the parser
assert [k for k, d in PF.items() if d["check_src"]] == ["D03", "D04", "D15", "D18"]

# Starts that Part F gives mostly in words: the same statement as one display line (the words stay in "In words").
START_TEX = {
    "D03": r"\zeta=Re^{i\theta}\ \ \text{(a point of the raised circle)}\quad\text{put into}\quad z=\zeta+\dfrac{b^2}{\zeta}\ \text{(14.3)}",
    "D04": r"\lvert\zeta+be\rvert=a,\qquad a\equiv b(1+e)\ \text{(14.7)}",
    "D05": r"z=\zeta+\dfrac{b^2}{\zeta}\ \text{(14.3)}\quad\Longrightarrow\quad\dfrac{dz}{d\zeta}\ \text{(the local magnification)}",
    "D07": r"w(z)=w\big(\zeta(z)\big),\qquad z=\zeta+\dfrac{b^2}{\zeta}\ \text{(14.3)}",
    "D11": r"\Gamma_{\text{tube}}=\text{the same all along a vortex tube (Helmholtz)}",
    "D13": r"\mathbf U_e=\mathbf U+\mathbf w,\qquad\mathbf w\perp\mathbf U",
}
START_WORDS = {
    "D03": "A point of the circle of radius a = b/cos β centred at Q = ib tan β (it passes through ζ = ±b), written by its distance R from the origin and its polar angle θ.",
    "D04": "The circle of radius a centred at ζ = −be — slightly larger than b, slightly to the left — so that it still passes through ζ = b.",
    "D05": "How the map stretches and turns a small step dζ.",
    "D07": "The complex potential has the same value at partner points: same potential, different coordinates.",
    "D11": "The strength of a vortex tube is the same all along it; a tube cannot end in the fluid (recap R11).",
    "D13": "The velocity triangle at a strip: the stream U along the flight path and the downwash w at right angles to it — the air the strip actually meets.",
}
for _k in START_TEX:
    PF[_k]["start"] = (START_TEX[_k], START_WORDS[_k])

if "--dump" in sys.argv:
    for _k, _d in PF.items():
        print(_k, _d["stars"], len(_d["steps"]), _d["title"])
        print("   start:", _d["start"])
        print("   result:", _d["result"])
        for _i, _s in enumerate(_d["steps"], 1):
            print(f"   {_i}. {_s['did']} | {_s['tex'][:90]}")
    sys.exit(0)

NOTES_IN: dict[str, str] = {
    "D01": r"**N03 [B]** the mean chord $\bar c=A/s$ (14.2) (step 1) · **N05 [B]** steady level flight, $L=W$ and thrust = drag (steps 3–4) · the stall speed of **N73** (step 5; used in C15)",
    "D03": r"**N19 [B]** a point of the circle in polar form from the origin (the Start) · **N20 [B]** $x=(R+b^2/R)\cos\theta$, $y=(R-b^2/R)\sin\theta$ (14.4) (step 2) · **N21 [B]** $x^2\sin^2\theta-y^2\cos^2\theta=4b^2\sin^2\theta\cos^2\theta$ (14.5) (steps 3–5) · **N22 [B]** the cosine rule in the triangle OQP (steps 6–7) · **N23 [B]** $2b\tan\beta\,\sin\theta=R-b^2/R=y/\sin\theta$ (14.6) (steps 8–9) · **N24 [B]** the circular arc $x^2+(y+2b\cot2\beta)^2=(2b\csc2\beta)^2$ (steps 10–14) · the camber rule of **N26** (step 15)",
    "D05": r"**N27 [B]** the trailing edge of a Zhukhovsky section is a cusp; an edge of finite angle would carry a stagnation point",
    "D06": r"**N28 [B]** the flow past the off-centre circle at incidence, $w=U\Big[(\zeta-\mu)e^{-i\alpha}+\dfrac{a^2e^{i\alpha}}{\zeta-\mu}\Big]+\dfrac{i\Gamma}{2\pi}\ln\dfrac{\zeta-\mu}{a}$ (not in the book — ours: the book writes only the centred case)",
    "D09": r"**N30 [B]** the Kutta circulation $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) (step 6), clockwise",
    "D10": r"**N09 [B]** $F_y=-\displaystyle\int_0^cp_u\,dx+\int_0^cp_l\,dx$ and $C_y\equiv\dfrac{F_y}{\tfrac12\rho U^2c}=\displaystyle\oint C_p\,d\Big(\dfrac xc\Big)$: the normal-force coefficient is the area of the pressure loop",
    "D13": r"**N42 [B]** the induced angle $\varepsilon=\tan^{-1}(w/U)\simeq w/U$ (step 1) · **N43 [B]** the effective angle $\alpha_e=\alpha-\varepsilon=\alpha-\dfrac wU$ (14.14) (step 2) · **N44 [B]** the strip's forces $dL=\rho U\Gamma\,dy$ and $dD_i=\rho w\Gamma\,dy$ (steps 3–6)",
    "D15": r"**N51 [B]** Glauert's integral $\displaystyle\int_0^\pi\dfrac{\cos n\gamma\,d\gamma}{\cos\gamma-\cos\gamma_1}=\pi\,\dfrac{\sin n\gamma_1}{\sin\gamma_1}$ (not in the book — ours, standard: the book says only that the last equality comes from evaluating the integral)",
    "D16": r"**N46 [B]** the section law on a finite wing $C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.16) (step 1) · **N47 [B]** $\Gamma=\tfrac12UcC_L$ (step 2) · **N48 [B]** the lifting-line equation $\Gamma(y)=\dfrac K2Uc(y)\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.17) (step 3)",
    "D17": r"**N53 [B]** the lift needs only the first coefficient, $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21)",
    "D19": r"**N55 [B]** the elliptic loading $\Gamma=\Gamma_1\sin(\gamma)=\Gamma_1\sqrt{1-(2y/s)^2}$ (14.23) (steps 3–4) · **N56 [B]** the uniform downwash $w(y)=\Gamma_1/2s$ (14.24) (step 5) · **N57 [B]** $D_i=\dfrac\pi8\rho\Gamma_1^2=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25) (steps 6–7)",
    "D20": r"**N59 [B]** the lift slope of a finite wing, $\dfrac{dC_L}{d\alpha}=\dfrac{K}{1+K/(\pi\Lambda)}$, i.e. $C_{L,3D}=\dfrac{C_{L,2D}}{1+2/\Lambda}$ for K = 2π (ours)",
    "D21": r"**N74 [B]** the best glide, $\Big(\dfrac LD\Big)_{\max}=\tfrac12\sqrt{\dfrac{\pi\Lambda e_{\text{span}}}{C_{D0}}}$ at $C_L=\sqrt{\pi\Lambda e_{\text{span}}C_{D0}}$ (ours; the book says only that the induced drag grows as α²)",
    "D22": r"**N81 [B]** the keel as a second wing and the closest course, $\beta_a=\varepsilon_A+\varepsilon_H$ with $\varepsilon=\tan^{-1}(D/L)$ of sail and of hull (steps 7–8; not in the book — ours)",
}
STARS = {1: "★", 2: "★★", 3: "★★★"}


def _cap(s: str) -> str:
    """Capital first letter when the text starts with an ordinary lower-case word (never a symbol or a maths span)."""
    return s[:1].upper() + s[1:] if re.match(r"[a-z]{2,}(?=[\s,;:.—–-])", s) else s


def pf_sub(key: str, field: str, old: str, new: str) -> None:
    """Edit one Part F field after parsing; fails loudly if the design text changed. ``field`` = goal/check/meaning/traps/
    assumptions, or ``stepN.did|tex|why|plain``."""
    d = PF[key]
    if field.startswith("step"):
        i, part = field[4:].split(".")
        st = d["steps"][int(i) - 1]
        assert old in st[part], (key, field, old)
        st[part] = st[part].replace(old, new)
        return
    assert old in d[field], (key, field, old)
    d[field] = d[field].replace(old, new)


def D(key: str, title: str, ref: str = "", after: str = "") -> None:
    """A Part F derivation, copied word for word (see ``part_f``), then the notes it carries and its traps as a ⚠️ callout.
    ``title`` is our heading (Part F headings name bare equation numbers); ``ref`` = "14.11" shows the key equation."""
    d = PF[key]
    if ref in EQ:
        ref = f"{ref}: ${EQ[ref]}$"                     # the heading shows the key equation, not only its number
    goal = _cap(d["goal"]) + (f"\n\n**Assumptions.** {_cap(d['assumptions'])}" if d["assumptions"] else "")
    nb.derivation(key, f"{title} `{key}` {STARS[d['stars']]}", ref=ref, goal=goal,
                  start=(d["start"][0], _cap(d["start"][1])), plan=[_cap(x) for x in d["plan"]], uses=d["tools"],
                  steps=d["steps"], result=(d["result"][0], _cap(d["result"][1])), interpret=_cap(d["meaning"]),
                  check=_cap(d["check"]), check_src=d["check_src"])
    if key in NOTES_IN:
        nb.md(f"📝 **Notes taught inside `{key}`:** {NOTES_IN[key]}.")
    if after:
        nb.md(textwrap.dedent(after).strip())
    if d["traps"]:
        nb.md(f"> ⚠️ **Common confusion (traps in `{key}`):** {_cap(d['traps'])}")


# ---------------------------------------------------------------------------------------------------------------------
# Small wrappers so that every block has the same shape
# ---------------------------------------------------------------------------------------------------------------------
def core(cid: str, title: str, question: str) -> None:
    """A CORE block heading with its id (the title carries the block's key equation in LaTeX)."""
    nb.core(cid, f"{title} `{cid}`", question=textwrap.dedent(question).strip())


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` is exactly the Part E concept text (the ledger check matches it)."""
    nb.primer(term, f"`{pid}` · {textwrap.dedent(text).strip()}", code=textwrap.dedent(code).strip("\n") if code else None)


def recap(rid: str, title: str, text: str, where: str, core_id: str | None = None) -> None:
    """A RECAP row. ``nb.recap`` closes the current CORE block; ``core_id`` re-enters it (the design's re-entry rule)."""
    nb.recap(rid, f"{title} `{rid}`", text, where)
    if core_id:
        nb.current_core = core_id


def recap_md(title: str, where: str, text: str) -> None:
    """An id-less reminder of an earlier chapter's idea in the 🔁 recap format; it does not close the block."""
    nb.md(f"> 🔁 **Recap — {title}** ({where}). {textwrap.dedent(text).strip()}", tags=["recap"])


def note(nid: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note whose text opens with its curation id(s) in bold ("**N09 [B]**")."""
    nb.note(f"**{nid}** {textwrap.dedent(text).strip()}", equation=equation, ref=ref)


def problem(text: str) -> None:
    nb.md("#### The problem in plain words\n\n" + textwrap.dedent(text).strip())


def idea(sketch: str = "", words: str = "") -> None:
    """The idea block: an ASCII sketch (or a table) plus one or two sentences."""
    sk = textwrap.dedent(sketch).strip("\n")
    body = ("```\n" + sk + "\n```") if sk else ""
    nb.md("#### The idea\n\n" + body + (("\n\n" if body else "") + textwrap.dedent(words).strip() if words else ""))


def maths(text: str) -> None:
    nb.md("#### The maths, step by step\n\n" + textwrap.dedent(text).strip())


def whatif(text: str) -> None:
    nb.md(f"**What would change if…** {textwrap.dedent(text).strip()}")


def confusion(text: str) -> None:
    nb.md("> ⚠️ **Common confusion:** " + textwrap.dedent(text).strip())


def slip(k: int, printed: str, correct: str) -> None:
    """A printed-vs-correct box (the book's printed slip #k, taught in corrected form)."""
    nb.md(f"> ⚠️ **slip #{k} — the book prints** {textwrap.dedent(printed).strip()}**; the correct form is** {textwrap.dedent(correct).strip()}")


def scratch(src: str, explain: str) -> None:
    """A from-scratch cell (hand-written version next to the tested function, then an assert) with its explanation."""
    nb.code(src, explain=explain, tags=["from-scratch"])


def explainer(slug: str, heading: str, why_static: str, tries: list[str]) -> None:
    """An embedded explainer, opening with the "why interactive rather than a static figure" sentence."""
    why = (f"**Why interactive rather than a static figure:** {textwrap.dedent(why_static).strip()}\n\n"
           "Start with the **Walkthrough**; the **Explain** tab works every number out with your settings, and the "
           "**Derivation** tab (where there is one) steps through the same derivation as above. Γ is the clockwise "
           "circulation in every readout; angles are shown in degrees and computed in radians.")
    nb.explainer(slug, heading=heading, why=why, tries=tries)


# ---------------------------------------------------------------------------------------------------------------------
# One-line reminders of tools primed in earlier chapters (Part E rows that point to knowledge/primers.md)
# ---------------------------------------------------------------------------------------------------------------------
def _ledger_rows() -> list[list[str]]:
    part = DESIGN.split("## Part E", 1)[1].split("## Part F", 1)[0]
    rows = [[x.strip() for x in ln.strip().strip("|").split("|")] for ln in part.splitlines() if ln.strip().startswith("|")]
    return [r for r in rows if len(r) >= 3 and not set("".join(r)) <= set("-: ") and r[0].lower() != "concept"]


LEDGER = _ledger_rows()
assert len(LEDGER) == 201, len(LEDGER)
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
            g = re.sub(r",?\s*\(?D\d\d\)?", "", g)                                 # pointers to another chapter's derivations
            g = g.replace("$", "").replace("\\", "")
            g = re.sub(r"_\{([^{}]*)\}", r"_\1", g)
            g = re.sub(r"\^\{([^{}]*)\}", r"^(\1)", g)
            cut = re.split(r"(?<=[^0-9])[;.] (?=[A-Za-z⚠(`])", g)                # first clause or two
            out = cut[0]
            for extra in cut[1:]:
                if len(out) + len(extra) > 200:
                    break
                out += "; " + extra
            return out.strip().rstrip(";.") + "."
    return ""


REMIND: dict[str, str] = {        # our own sentence where the way the tool is used HERE needs saying
    "true wind, boat velocity, apparent wind; `ch14.apparent_wind`": "the wind a moving body feels is the air's velocity minus its own, $\\mathbf U_r=\\mathbf U-\\mathbf V$ — primed in C10 of this chapter (P350); the sail lives in it.",
    "`slider_figure`, `animate` / `show_animation`, `show_viz`, `live`": "`slider_figure(fn, name, values)` draws a plotly figure with every slider position precomputed (works on the web page); `animate` + `show_animation` play a matplotlib animation inline; `show_viz` embeds an explainer; `live` makes kernel-only sliders.",
    "sympy (`symbols`, `simplify`, `series`, `limit`, `integrate`)": "sympy does algebra with symbols: `sp.symbols` names them, `sp.simplify(a - b) == 0` proves two expressions equal, `sp.series` expands in a small parameter, `sp.limit` takes a limit, `sp.integrate` integrates.",
    "`lambda`": "`lambda y: expr` is a one-line function without a name, handy for handing a loading Γ(y) or a twist α(y) to a solver.",
}


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
        if concept in REMIND:
            text = REMIND[concept]
        else:
            text = " ".join(g for g in (_gist(p) for p in pids) if g) or "used here exactly as in the chapter named."
        lines.append(f"- **{concept}** — {text} *({where})*")
        nb.primers.append(re.sub(r"[`*$\\]", "", concept) + " (reminder)")
    head = "> 🔁 **Tools from earlier chapters used in this block** (one line each; the full primers are in the chapters named)"
    nb.md(f"{head}{(' — ' + lead) if lead else ''}\n\n" + "\n".join(lines), tags=["primer"])


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=r"""
A wing weighing a few tonnes holds up an aircraft weighing hundreds. How? Chapter 6 gave half the answer: a body in an ideal
stream feels a force at right angles to the stream whenever the flow circulates round it, $L=\rho U\Gamma$ (6.62) — but it left
the circulation $\Gamma$ free. This chapter supplies the missing half in three sentences. **Lift is circulation.** **The sharp
trailing edge chooses it** — viscosity, acting only there, lets exactly one circulation survive, and that gives a lift that
grows by $2\pi$ for every radian of incidence. **The wing tips make you pay for it** — the circulation cannot stop at the tips,
so it trails behind as two vortices whose downwash tilts the lift backward: a drag that exists even without friction, and that
falls as the wing gets long and slender. The last two sections turn the first sentence round: a fish's tail, a bird's wing and
a sail all move a wing *across* the stream so that its lift leans forward and becomes thrust.""",
    roadmap=[
        "**The wing in numbers** — span, area, aspect ratio, and the lift coefficient level flight demands (`C01`).",
        "**The Kutta condition** — ideal flow allows any circulation; the sharp trailing edge allows one (`C02`).",
        "**The starting vortex** — the wing gets its circulation by shedding the opposite one (`C03`).",
        "**The Zhukhovsky airfoil** — a circle moved off-centre and pushed through one map (`C04`).",
        "**The flow round the airfoil** — the circle's flow, rescaled point by point (`C05`).",
        "**The lift of the airfoil** — $2\\pi$ per radian, shifted by the camber (`C06`).",
        "**Surface speed and pressure** — the lift as the area of a pressure loop (`C07`).",
        "**Bound vortex and trailing filaments** — what leaves the wing where the circulation drops (`C08`).",
        "**Downwash** — the sinking air the wake makes at the wing (`C09`).",
        "**Induced drag** — lift tilted back: drag without friction (`C10`).",
        "**The sine series and its downwash** — one singular integral, done once (`C11`).",
        "**The lifting-line equation** — N stations, N sines, one linear solve (`C12`).",
        "**Lift and induced drag from the coefficients** — lift in the first sine, drag in all (`C13`).",
        "**The elliptic loading** — uniform downwash, least drag, and why span matters (`C14`).",
        "**The drag breakdown and the polar** — profile against induced, and the best glide (`C15`).",
        "**Thrust from lift** — a fish tail, a wing beat and a sail (`C16`).",
    ],
    prerequisites=[
        r"the complex potential, the cylinder with circulation, the lift theorem $L=\rho U\Gamma$ (6.62) and the map $z=\zeta+\dfrac{b^2}{\zeta}$ (6.65) (Ch. 6)",
        r"Kelvin's theorem $\dfrac{D\Gamma}{Dt}=0$ (5.8), Helmholtz's theorems and the filament law $d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17) (Ch. 5)",
        r"the coefficients $C_L\equiv\dfrac{L}{\tfrac12\rho U^2A}$ and $C_D\equiv\dfrac{D}{\tfrac12\rho U^2A}$ (4.107), (4.108) (Ch. 4)",
        "boundary layers and separation (Ch. 9)",
    ],
)
nb.explainer_index([
    ("zhukhovsky_airfoil_builder", "Where do camber and thickness come from?", "drag the circle's centre: up bends the plate, back fattens the nose"),
    ("kutta_condition", "Which circulation does a wing choose?", "only one puts the rear stagnation point on the sharp edge"),
    ("starting_vortex", "Where does a wing's circulation come from?", "it sheds an equal and opposite vortex and flies away from it"),
    ("airfoil_pressure_lift", "Is a wing pushed up or pulled up?", "the lift is the area between two pressure curves, mostly suction near the nose"),
    ("lift_curve_stall", "Does more angle always mean more lift?", "camber shifts the line, the angle buys 2π per radian, stall ends it"),
    ("downwash_induced_drag", "How can there be drag without friction?", "the wing flies in its own sinking air, and its lift leans back"),
    ("lifting_line_solver", "What makes one wing better than another?", "lift lives in the first sine; every other one is pure drag"),
    ("elliptic_loading_aspect_ratio", "Why are glider wings so long and thin?", "the price of lift is set by the span you spread it over"),
    ("thrust_from_lift", "How can a boat sail against the wind?", "lift is perpendicular to the wind you feel, not the wind that blows"),
])
nb.setup()

nb.md(r"""
## 📦 Chapter imports and the numbers we will use

Every worked number in this notebook uses **our own illustrative inputs**, written out in the cell below as five small
dictionaries: a wing section (`sec`), a light aircraft (`ac`), a small demonstration wing (`wd`), a sailing boat (`boat`)
and a fish tail (`fish`). None of them is a number from the book.""")
nb.code(r"""
import numpy as np                                          # arrays and maths
import sympy as sp                                          # symbolic algebra (derivation checks)
import matplotlib.pyplot as plt                             # static figures and animations
import pandas as pd                                         # small tables
from scipy.integrate import quad                            # numerical integration of a function
from scipy.optimize import brentq                           # root of f(x) = 0 inside a bracket
from fluidpy import ch14_aerodynamics as ch14               # this chapter's tested functions (re-exports AF and LL)
from fluidpy.core import airfoils as AF, lifting_line as LL  # the two new core modules: sections, finite wings
from fluidpy.core import conformal as CONF, potential as POT, biot_savart as BS   # Ch. 6 maps and potentials, Ch. 5 vortex filaments
from fluidpy.core import similarity as SIM, boundary_layer as BL, vortices as VORT  # Ch. 4 coefficients, Ch. 9 layers, Ch. 3 vortices
from fluidpy import ch06_ideal_flow as ch06                  # the cylinder with circulation and the lift theorem
from fluidpy.core.interact import slider_figure, animate_figure, live   # plotly slider figures and kernel-only widgets
from fluidpy.core.anim import animate                        # builds a matplotlib animation from an update function
from fluidpy.core.style import COLORS, savefig               # house colours; save a figure to outputs/ch14/
sys.path.insert(0, "scripts")                                # the folder with this chapter's drawing helpers
from ch14_drawings import *                                  # sketches (draw_…) and composite figures (fig_…) built on the tested functions
from ch14_common import streamlines                          # contours of a section's stream function with the section filled in
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # silence harmless "font weight not found" messages

sec = dict(b=0.25, e=0.08, beta=np.deg2rad(4.0), alpha=np.deg2rad(5.0), U=30.0, rho=1.225)   # section: map constant b [m] (chord ≈ 4b = 1 m), thickness parameter e, camber angle β = 4°, angle of attack α = 5° (both stored in rad), speed [m/s], air density [kg/m³]
ac = dict(W=9000.0, s=12.0, A=15.0, U_slow=30.0, U_cruise=55.0, CD0=0.025, e_span=0.9, CL_max=1.4, CL_max_flap=2.0, K=2*np.pi)   # light aircraft: weight [N], span [m], planform area [m²], two speeds [m/s], profile drag coefficient, span efficiency, maximum lift coefficients, section lift slope [1/rad]
wd = dict(s=6.0, c=1.0, alpha=np.deg2rad(5.0), U=1.0, rho=1.0, K=2*np.pi)   # demo wing: rectangular, span 6 m, chord 1 m (aspect ratio 6), α = 5° in rad, unit speed and density
boat = dict(U_wind=6.0, V_boat=3.0, theta=np.deg2rad(45.0), A_sail=20.0, CL=1.2, CD=0.2, LD_hull=8.0, CDA_hull=0.05, rho_water=1000.0)   # boat: true wind and boat speed [m/s], course 45° off the wind [rad], sail area [m²], sail coefficients, hull lift/drag, hull drag area [m²], water density [kg/m³]
fish = dict(U=1.0, V=0.75, L=10.0, D=1.0)                    # fish: swimming speed, sideways tail speed [m/s], tail lift and drag [N]
nu_air = 1.46e-5                                             # kinematic viscosity of sea-level air [m²/s]
SEC = {k: sec[k] for k in ("U", "b", "e", "beta", "alpha")}   # the five numbers every section function takes, so that they can be passed by name: f(**SEC)
mu, a = AF.zhukhovsky_circle(b=sec["b"], e=sec["e"], beta=sec["beta"])   # centre μ [m] (complex) and radius a [m] of the circle behind our section (taught in C04)
print(f"circle of our section: centre mu = {mu:.4f} m, radius a = {a:.4f} m")   # −0.0200+0.0189j and 0.2707
pd.Series(ch14.illustrative_inputs()).to_frame("the module's own illustrative set (its scripts and tests use it; also ours)").T   # shown once, beside ours
""", explain=r"""
1. The first block imports numpy, sympy, matplotlib, pandas and two scipy tools, then the chapter module `ch14` and the two new
   core modules: `AF` (airfoil sections) and `LL` (lifting line). `ch14` re-exports every public name of both.
2. The second block brings back earlier chapters' tools: `CONF` and `POT` (complex potentials and maps, Ch. 6), `BS`
   (vortex filaments, Ch. 5), `SIM` (force coefficients, Ch. 4), `BL` (boundary layers, Ch. 9), `VORT` (Ch. 3) and `ch06`.
3. `from ch14_drawings import *` loads sketch helpers (`draw_…`, `fig_…`): pure matplotlib on top of the tested functions.
4. The five dictionaries are the inputs of every worked number below. Angles are **stored in radians** (`np.deg2rad`) and
   **quoted in degrees** in the text.
5. `mu, a` is the circle that will become our wing section: centre −0.0200 + 0.0189i m, radius 0.2707 m.
6. The last line shows the module's own illustrative set once. It overlaps with ours; both are ours, neither is the book's.""")

nb.md(r"""
## ⚠️ Conventions in this chapter — read once, return when a sign surprises you

**(a) Two sets of axes.** In §14.3–§14.5 the wing *section* lies in the $x$–$y$ plane with $x$ along the chord, $y$ up, and the
stream inclined at the angle of attack $\alpha$ to $x$ — so **lift is not $F_y$**. In §14.6 $y$ runs along the *span* and the
downwash $w$ is positive **downward** (trap T5). The sketch below shows both.

**(b) $\Gamma$ is the clockwise circulation.** The book counts circulation clockwise-positive, so that a wing with positive lift
has $\Gamma>0$ in $L=\rho U\Gamma$ (6.62), $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10) and $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11).
Chapters 3 and 5 of this project counted counterclockwise-positive (trap T1).

| Convention | Sign of $\Gamma$ for a lifting wing | Lift theorem | In code |
|---|---|---|---|
| the book, this chapter | $\Gamma>0$ (clockwise) | $L=+\rho U\Gamma$ | `Gamma_cw=` ; returned values are named `Gamma_cw` |
| Ch. 3, Ch. 5 of this project | $\Gamma_{ccw}<0$ | $L=-\rho U\Gamma_{ccw}$ | `Gamma_ccw=` |

**(c) The complex velocity is $dw/dz=u-iv$** — with a minus sign (trap T2). Chapter 13's Ekman velocity was $V=u+iv$.

**(d) Overloaded letters** (traps T3, T4, T13): the table printed by the cell below lists every letter that means more than one
thing — $b$, $a$, $\beta$, $c$, $e$, $\Lambda$, $\gamma$, $\theta$, $w$, $A$, $\Gamma$ — with the meaning *here*.

**(e) Slopes of lift curves are per radian** in every formula and function (the keyword `K=`); "per degree" appears only in
one conversion (trap T12). **Angles are computed in radians and shown in degrees**; where a number matters both are printed.

**(f) No lapse rate here.** The letter $\Gamma$ is circulation throughout this chapter. The one place temperature enters (the
speed of sound in §14.1) uses the sea-level temperature only.

**(g) Slips and traps.** *Slips* are statements that are false as printed in the book; we teach the corrected form in a box
"slip #k — the book prints … ; the correct form is …". *Traps* are true as printed but easy to misread; they are the
"Common confusion" callouts.

**(h) Which inputs?** All worked numbers use our own illustrative inputs (the five dictionaries above), never the book's
(trap T17).""")
nb.code(r"""
fig = draw_axes_conventions()                               # two small sketches: the section plane (x along the chord, y up, stream at α) and the span (y along the span, w down)
plt.show()                                                  # display the sketch
display(ch14.conventions_table())                           # the letter table: symbol, its meanings, the one used here, its code name
""", tags=["figure"], explain=r"""
1. `draw_axes_conventions()` draws the two coordinate systems of the chapter side by side.
2. `ch14.conventions_table()` is a pandas table (`pd.DataFrame`: rows and named columns) of the overloaded letters.""")
nb.figure_notes(
    see="Left: a wing section with $x$ along its chord and the stream arriving at the angle $\\alpha$. Right: a wing from "
        "above with $y$ along the span and the downwash $w$ pointing down.",
    read="The left axes belong to §14.3–§14.5 (one section, forces per metre of span), the right ones to §14.6 (a whole wing). "
         "Because the stream is tilted on the left, the lift — which is across the *stream* — is not the $y$-force.",
    change="…you tilted the drawing instead of the stream? Then lift would be $F_y$ — but the map of §14.4 needs the chord "
           "on the $x$-axis, so the book tilts the stream.")
nb.code(r"""
slips = pd.DataFrame(ch14.book_slips()).T[["where", "printed", "correct"]]   # 10 rows: where it is, what is printed, what is right
traps = pd.DataFrame(ch14.traps())[["id", "what", "where"]]                   # 18 rows: true as printed, easy to misread
pd.set_option("display.max_colwidth", 110)                                    # show long cells without cutting them
display(slips)                                                                # slips: false as printed
display(traps)                                                                # traps: true but easy to misread
""", explain=r"""
1. `ch14.book_slips()` returns the ten suspected misprints (a dictionary keyed "1" … "10"); `.T` turns it into one row per slip.
2. `ch14.traps()` returns the eighteen traps T1–T18.
3. Slips are false as printed; traps are true but easy to misread. Each is taught where it first bites.""")

# =====================================================================================================================
# A.1 §14.1 Introduction
# =====================================================================================================================
nb.section("14.1", "Introduction", intro=r"""
**What is this section about?** What "aerodynamics" covers, the three assumptions the whole chapter leans on, and the one
place where they fail.""")
note("N01 [B]", r"""
Aerodynamics asks for the forces and moments a moving fluid exerts on a body. Three standing assumptions:

| # | Assumption | Taught in | Fails when |
|---|---|---|---|
| 1 | no buoyancy or stratification: the air's density is uniform, so the ideal (inviscid, irrotational) flow of Ch. 6 applies | Ch. 6 | the air is compressed — Mach number above about 0.3 |
| 2 | a uniform stream far ahead of the body | Ch. 6 | gusts, a neighbour's wake |
| 3 | viscosity matters only in thin boundary layers on the surface | Ch. 9 | **stall**: the boundary layer leaves the surface (N13, in C07) |

The air may be treated as incompressible while the Mach number $\mathrm{Ma}=U/c_{\text{sound}}$ (speed over the speed of sound,
Ch. 1 and Ch. 4) stays below about 0.3, because density changes are about $\tfrac12\mathrm{Ma}^2$ (ours).""")
nb.code(r"""
rows = []                                                   # one row per flight condition
for name, U, T in (("our aircraft, slow", 30.0, 288.15), ("our aircraft, cruise", 55.0, 288.15), ("an airliner at cruise", 240.0, 220.0)):   # speed [m/s], air temperature [K]
    cc = ch14.compressibility_check(U=U, T=T)                 # Mach number, density change ½Ma², verdict Ma < 0.3
    rows.append({"case": name, "U [m/s]": U, "Ma": round(cc["Ma"], 3), "density change [%]": round(100*cc["density_change"], 1), "incompressible?": cc["incompressible_ok"]})   # collect
pd.DataFrame(rows)                                          # show the three rows as a table
""", explain=r"""
**What does this show?** Our light aircraft passes at both speeds (Ma = 0.088 and 0.162: the density changes by 0.4 % and
1.3 %). An airliner at cruise (Ma ≈ 0.81, a 33 % change) fails — that is Chapter 15. Everything below is for the first two rows.""")

# =====================================================================================================================
# A.2 §14.2 Aircraft Terminology — C01
# =====================================================================================================================
nb.section("14.2", "Aircraft Terminology", intro=r"""
**What is this section about?** The handful of words and numbers that describe a wing — span, chord, area, aspect ratio, angle
of attack — and the two coefficients every later result is written in.""")
core("C01", r"The wing in numbers: $\Lambda\equiv s^2/A=s/\bar c$ (14.1), $\bar c=A/s$ (14.2)",
     "How do you describe a wing with a few numbers — and how much lift must it deliver to stay up?")
problem(r"""
Stand under an airliner and under a glider. One wing is short and broad, the other long and thin. Before any physics we need
words for that difference, and for the angle at which the wing meets the air. Then one balance — in steady level flight the lift
equals the weight — tells us what the wing has to do, and everything else in the chapter is about *how* it does it and *what it
costs*.""")
idea(r"""
            top view of a wing (looking down)

   tip  ___________________________________________  tip
       |                     |                     |
       |        c(y)  ←chord→|   root              |      A = area seen from above (planform area)
       |_____________________|_____________________|
       |<------------------ s (span) ------------->|

   slender:  Λ = s²/A large  (glider ≈ 20)          stubby:  Λ small  (fighter ≈ 3)
""", r"""
Three lengths describe the outline: the **span** $s$ (tip to tip), the **chord** $c(y)$ (front to back at the spanwise station
$y$) and the **planform area** $A$ (the area seen from above). One ratio, the **aspect ratio** $\Lambda$, says how slender it is.""")
remind("C01")
note("N02 [B]", r"""
**Words for a planform.** The *fuselage* is the body; the wing joins it at the **root** and ends at the **tips**; the **span** $s$
is the tip-to-tip distance; the **chord** $c(y)$ is the wing's length along the flight direction at station $y$; the **planform
area** is $A=\int_{-s/2}^{s/2}c(y)\,dy$ (ours, implied by the book's definition). **Taper** is a chord that shrinks toward the tip
(taper ratio = tip chord / root chord); **sweepback** is a wing whose tips lie behind its root.""")
nb.figure(r"""
fig = draw_planform_labelled()                              # a labelled planform: span s, chord c(y), root, tips, area A
plt.show()                                                  # display the sketch
shapes = {"rectangular": LL.planform(kind="rectangular", s=12.0, c_root=1.25),          # span 12 m, chord 1.25 m everywhere
          "tapered (tip/root = 0.4)": LL.planform(kind="tapered", s=12.0, c_root=1.786, taper=0.4),   # root chord 1.786 m, tip chord 0.4 of it
          "elliptic": LL.planform(kind="elliptic", s=12.0, c_root=1.592)}                # root chord 1.592 m, chord falls like an ellipse
fig, axes = plt.subplots(1, 3, figsize=(10.5, 2.6), sharey=True)           # three outlines side by side
for ax, (name, Pf) in zip(axes, shapes.items()):            # one panel per planform
    ax.fill_between(Pf["y"], 0.25*Pf["c"], -0.75*Pf["c"], color=COLORS["teal"], alpha=0.35)   # outline: leading edge ¼ chord ahead of a straight line, trailing edge ¾ behind
    ax.set_title(f"{name}\nA = {Pf['A']:.1f} m², Λ = {Pf['AR']:.1f}, mean chord = {Pf['c_mean']:.2f} m", fontsize=9)   # area, aspect ratio, mean chord
    ax.set_xlabel("spanwise position $y$ [m]")              # along the span
    ax.set_aspect("equal")                                  # true proportions
axes[0].set_ylabel("flight direction [m]")                  # up the page = forward
plt.show()                                                  # display
""",
    see="A labelled wing outline, then three wings with the **same span (12 m) and the same area (15 m²)**: a rectangle, a "
        "tapered wing and an ellipse.",
    read="All three have the same aspect ratio 9.6 and the same mean chord 1.25 m: those two numbers cannot tell them apart. "
         "What differs is how the chord is distributed along the span — which C12 will show decides how the lift is distributed.",
    change="…the taper went to 0? A pointed tip: the chord vanishes there, and C12 shows the loading it produces is too "
           "peaked in the middle.")
D("D01", "Aspect ratio two ways, and what level flight asks of the wing", ref="14.1")
P("P336", "force per unit span, section coefficient and wing coefficient", r"""
A two-dimensional section has no area, so its forces are **per metre of span** [N/m] and are divided by $\tfrac12\rho U^2c$
(chord $c$). A whole wing's forces are in newtons and are divided by $\tfrac12\rho U^2A$ (planform area $A$).
$\tfrac12\rho U^2$ is the **dynamic pressure**: the pressure rise where the stream is brought to rest.""", r"""
rho, U, c, A = 1.225, 30.0, 1.0, 15.0     # air density [kg/m³], speed [m/s], chord [m], planform area [m²]
q = 0.5*rho*U**2                          # dynamic pressure: 551 Pa
print(q, 587.0/(q*c), 9000.0/(q*A))       # a section's 587 N/m -> 1.06 ; a wing's 9000 N -> 1.09
""")
recap("R01", "Lift and drag coefficients", r"""
$C_D\equiv\dfrac{D}{\tfrac12\rho U^2A}$ and $C_L\equiv\dfrac{L}{\tfrac12\rho U^2A}$ (4.107), (4.108): a force divided by dynamic pressure
times an area. **New here:** $A$ is the *planform* area (Chapter 9's bluff bodies used the frontal area — trap T13), and a
section uses the chord per unit span.

| Body | Area in the coefficient | Symbol |
|---|---|---|
| a wing | planform area (seen from above) | $A$ |
| a bluff body (sphere, car, projectile) | frontal area (seen from ahead) | $A_f$ |
| a two-dimensional section | chord × 1 m of span | $c$ |

> ⚠️ **slip #6 — the book prints** that the two definitions $C_D\equiv\dfrac{D}{\tfrac12\rho U^2A}$, $C_L\equiv\dfrac{L}{\tfrac12\rho U^2A}$ (4.107), (4.108) are "in Section 4.3"**; the correct form is** §4.11.""",
      where="Ch. 4 §4.11", core_id="C01")
nb.code(r"""
print(SIM.lift_coefficient(F=9000.0, rho=1.225, U=30.0, A=15.0))     # a wing: 9000 N over (551 Pa × 15 m²) = 1.088
print(AF.section_coefficient(F_per_span=587.0, rho=1.225, U=30.0, c=1.0))     # a section: 587 N/m over (551 Pa × 1 m) = 1.065
""")
recap("R02", "Friction drag and pressure drag", r"""
Drag is the sum of the surface shear stress (friction drag) and of the surface pressure (pressure or form drag); lift is almost
entirely a pressure force. C15 refines this into friction + form + induced.""",
      where="Ch. 9 (form drag); Ch. 4 (surface forces)", core_id="C01")
nb.worked_example("a wing by hand", r"""
$s$ = 10 m, $A$ = 20 m².

1. Mean chord $\bar c=A/s$ (14.2) = 20/10 = 2 m.
2. Aspect ratio $\Lambda\equiv s^2/A$ (14.1) = 100/20 = 5 — and also $s/\bar c$ = 10/2 = 5 ✓.
3. Level flight with $W$ = 10 000 N at $U$ = 40 m/s in air of ρ = 1.25 kg/m³: dynamic pressure ½ × 1.25 × 40² = 1000 Pa;
   $C_L=W/(\tfrac12\rho U^2A)$ = 10 000/(1000 × 20) = 0.5.
4. Halve the speed: the dynamic pressure falls to 250 Pa and the wing must deliver $C_L$ = 2.0 — four times as much.
5. If the wing can reach at most $C_{L,\max}$ = 1.4, the slowest level flight is
   $U_{\min}=\sqrt{2W/(\rho AC_{L,\max})}$ = √(20 000/35) = 23.9 m/s.""")
nb.code(r"""
AR = LL.aspect_ratio(s=ac["s"], A=ac["A"])                      # Eq. (14.1): Λ = s²/A [–]
cbar = LL.mean_chord(A=ac["A"], s=ac["s"])                      # Eq. (14.2): mean chord = A/s [m]
CL_slow = ch14.level_flight(W=ac["W"], rho=1.225, U=ac["U_slow"], A=ac["A"])["CL"]      # lift coefficient needed at 30 m/s: W/(½ρU²A)
CL_cruise = ch14.level_flight(W=ac["W"], rho=1.225, U=ac["U_cruise"], A=ac["A"])["CL"]  # the same at 55 m/s
U_clean = ch14.stall_speed(W=ac["W"], rho=1.225, A=ac["A"], CL_max=ac["CL_max"])             # slowest level speed, clean wing [m/s]
U_flap = ch14.stall_speed(W=ac["W"], rho=1.225, A=ac["A"], CL_max=ac["CL_max_flap"])         # slowest level speed with the flap down [m/s]
print(f"aspect ratio = {AR:.1f}, mean chord = {cbar:.2f} m")                  # 9.6 and 1.25 m
print(f"C_L needed: {CL_slow:.3f} at 30 m/s, {CL_cruise:.3f} at 55 m/s")      # 1.088 and 0.324
print(f"slowest level flight: {U_clean:.1f} m/s clean, {U_flap:.1f} m/s with flap")   # 26.5 and 22.1 m/s
""", explain=r"""
1. `LL.aspect_ratio` and `LL.mean_chord` are $\Lambda\equiv s^2/A$ (14.1) and $\bar c=A/s$ (14.2): our aircraft has Λ = 9.6 and a
   mean chord of 1.25 m.
2. `ch14.level_flight` solves $L=W$ for the lift coefficient: 1.088 at 30 m/s but only 0.324 at 55 m/s — the same weight,
   3.4 times the dynamic pressure.
3. `ch14.stall_speed` is step 5 of the derivation: 26.5 m/s with $C_{L,\max}$ = 1.4.
4. Flaps raise $C_{L,\max}$ to 2.0, so the aircraft can land 4.3 m/s slower (C15 returns to this).""")
scratch(r"""
y = np.linspace(-6.0, 6.0, 201)                             # 201 stations along a 12 m span [m]
c = 1.786*(1 - 0.6*np.abs(y)/6.0)                           # tapered chord: 1.786 m at the root, 0.4 of it at the tips [m]
A_mine = 0.0                                                # running planform area [m²]
for i in range(len(y) - 1):                                 # one trapezoid per pair of neighbouring stations
    A_mine += 0.5*(c[i] + c[i+1])*(y[i+1] - y[i])           # mean chord of the strip times its width
AR_mine = 12.0**2/A_mine                                    # Eq. (14.1): s²/A
cbar_mine = A_mine/12.0                                     # Eq. (14.2): A/s
Pf = LL.planform(kind="tapered", s=12.0, c_root=1.786, taper=0.4)   # the tested function for the same wing
assert np.allclose([A_mine, AR_mine, cbar_mine], [Pf["A"], Pf["AR"], Pf["c_mean"]], rtol=1e-3)   # same numbers → the library does exactly this
print(f"A = {A_mine:.2f} m², aspect ratio = {AR_mine:.2f}, mean chord = {cbar_mine:.3f} m")       # 15.0 m², 9.6, 1.25 m
""", r"""
1. The chord of a tapered wing is written as a formula and the area is summed strip by strip (the trapezoid rule, by hand).
2. Aspect ratio and mean chord follow from $\Lambda\equiv s^2/A$ (14.1) and $\bar c=A/s$ (14.2).
3. `assert np.allclose` confirms that `LL.planform` returns the same three numbers.""")
note("N07 [B]", r"""
**Words for a section.** Cut the wing with a vertical plane along the flight direction: the cut is the **airfoil section**. It
has a rounded **leading edge** and a sharp **trailing edge**. The **chord line** joins the two; the **camber line** runs midway
between the upper and lower surfaces; the **camber** is the largest height of the camber line above the chord, $\max y_c$; the
**thickness** is the largest distance between the surfaces, $t=\max(y_u-y_l)$; the **angle of attack** α is the angle between
the chord line and the stream.""")
P("P337", "splitting a closed contour into upper and lower surfaces", r"""
A section outline is one closed list of points. To measure it, find the trailing edge (largest x) and the leading edge
(smallest x; `np.argmin` returns its index), cut the list there, and interpolate both halves on the same x stations
(`np.interp(x_new, x, y)` reads a curve between its samples and needs x increasing, so one half is reversed with `[::-1]` —
primed in Ch. 7, P182).""", r"""
x = np.array([1, .5, 0, .5, 1.]); y = np.array([0, .1, 0, -.06, 0])   # a diamond, counterclockwise from the tail
i = np.argmin(x)                                                      # index of the leading edge: 2
yu = np.interp(0.5, x[:i+1][::-1], y[:i+1][::-1])                     # upper surface at x = 0.5 (first half, reversed so that x increases)
yl = np.interp(0.5, x[i:], y[i:])                                     # lower surface at x = 0.5 (second half)
print(yu - yl, 0.5*(yu + yl))                                         # thickness 0.16, camber 0.02
""")
nb.figure(r"""
S0 = AF.zhukhovsky_airfoil(b=sec["b"], e=sec["e"], beta=sec["beta"])  # the outline of our section (how it is made: C04)
geom = AF.section_geometry(x=S0["x"], y=S0["y"])                 # chord, camber and thickness measured on the outline [m]
fig = draw_section_labelled(x=S0["x"], y=S0["y"], geom=geom)          # the section with chord line, camber line, camber, thickness, α
plt.show()                                                   # display
print(f"chord = {geom['chord']:.3f} m, camber/chord = {geom['camber']/geom['chord']:.3f}, thickness/chord = {geom['t_max']/geom['chord']:.3f}")   # 1.006 m, 0.035, 0.096
""",
    see="Our wing section with its chord line (straight, nose to tail), its camber line (midway between the surfaces), the "
        "camber and the maximum thickness marked.",
    read="The printed numbers are measured on the outline: chord 1.006 m, camber 3.5 % of the chord, thickness 9.6 % of the "
         "chord. The camber is small — the camber line is only slightly bowed — yet C06 shows it is worth 4° of angle of attack.",
    change="…you changed α? Nothing in this picture: α is not a property of the shape. Tilt the stream, not the drawing.")
note("N08 [B]", r"""
The force on the section is resolved **relative to the stream**: drag $D$ along it, lift $L$ across it. If the stream makes the
angle α with the $x$-axis, $D=F_x\cos\alpha+F_y\sin\alpha$ and $L=-F_x\sin\alpha+F_y\cos\alpha$ (ours).

> ⚠️ **Common confusion (trap T5):** from §14.4 on the $x$-axis is the chord and the stream is tilted, so **lift is not $F_y$**.""")
nb.code(r"""
L_, D_ = AF.resolve_lift_drag(Fx=-51.1, Fy=584.4, alpha=sec["alpha"])    # force components along and across the chord [N/m] (computed in C07), stream at α = 5°
print(f"L = {L_:.1f} N/m, D = {D_:.2f} N/m")                 # 586.6 N/m and ≈ 0: the pressure force leans forward of the y-axis by exactly α
""")
note("N04 [C] · N06 [C] · N83 [C]", r"""
Three things named here and used later.

- **N04 [C] The three axes of an aircraft.** *Pitch* is nose up or down, *roll* is one wing down, *yaw* is nose left or right
  (sketch below). Roll returns with the ailerons below and with the rolling moment in a wake vortex (N64, in C12).
- **N06 [C] Control surfaces.** *Ailerons* (outer trailing edge, moved in opposite senses to roll), *elevator* (on the tail, to
  pitch), *rudder* (on the fin, to yaw), *flaps* (inner trailing edge, lowered together for more lift), *slots* (a gap near the
  nose that feeds fast air to the upper surface), *spoilers* (plates raised to destroy lift), *trim tabs* (small surfaces that
  hold a control in place). All of them work by changing the local camber; "more camber ⇒ more lift at the same α" is N31 in
  C06, flaps and landing speed N73 in C15.
- **N83 [C] A projectile with drag.** With a constant drag coefficient on a shallow trajectory the range is
  $x\simeq\dfrac{2m}{\rho AC_D}\ln\Big(1+\dfrac{\rho AC_DV_o^2\cos\theta_o\sin\theta_o}{mg}\Big)$ — here $A$ is the *frontal* area (trap T13); the
  next cell evaluates it once with our own numbers.""")
nb.code(r"""
fig = draw_three_axes()                                      # pitch, roll and yaw axes on a plan-view silhouette
plt.show()                                                   # display
x_drag = ch14.projectile_range_drag(m=0.145, rho=1.225, A=0.0042, CD=0.4, V0=30.0, theta0=np.deg2rad(10.0))   # a 145 g ball, frontal area 42 cm², C_D = 0.4, 30 m/s at 10° [m]
x_free = 2*30.0**2*np.cos(np.deg2rad(10.0))*np.sin(np.deg2rad(10.0))/9.80665            # the same throw without drag: 2 V² cos θ sin θ / g [m]
print(f"range with drag = {x_drag:.1f} m, without drag = {x_free:.1f} m")               # drag shortens the throw
""", tags=["figure"])
nb.plotly(r"""
def planform_taper(taper):                                   # outline of a wing of span 12 m and area 15 m² for one taper ratio
    c_root = 2*15.0/(12.0*(1 + taper))                       # root chord that keeps the area at 15 m²: A = s·c_root·(1 + taper)/2
    Pf = LL.planform(kind="tapered", s=12.0, c_root=c_root, taper=taper, n=101)  # stations y and chords c(y)
    yy = np.r_[Pf["y"], Pf["y"][::-1], Pf["y"][0]]           # go out along the leading edge and back along the trailing edge
    xx = np.r_[0.25*Pf["c"], -0.75*Pf["c"][::-1], 0.25*Pf["c"][0]]   # leading edge ¼ chord ahead, trailing edge ¾ chord behind
    return {"planform: s = 12 m, A = 15 m², mean chord 1.25 m, Λ = 9.6 at every taper": (yy, xx)}   # one closed outline
figF1a = slider_figure(planform_taper, "taper ratio (tip chord / root chord)", np.linspace(0.1, 1.0, 19), unit="",
                       xlabel="spanwise position y [m]", ylabel="flight direction [m]", yrange=[-2.6, 1.2],
                       title="Taper moves area inboard; span, area and aspect ratio stay the same")   # 19 slider positions
figF1a.show()                                                # display the slider figure

def planform_AR(AR):                                         # outline of a rectangular wing of area 15 m² for one aspect ratio
    s_ = np.sqrt(AR*15.0)                                    # Eq. (14.1) solved for the span: s = sqrt(Λ A) [m]
    c_ = 15.0/s_                                             # Eq. (14.2): the chord = A/s [m]
    return {"rectangle of area 15 m²": (np.array([-s_/2, s_/2, s_/2, -s_/2, -s_/2]), np.array([0.25*c_, 0.25*c_, -0.75*c_, -0.75*c_, 0.25*c_]))}   # four corners, closed
figF1b = slider_figure(planform_AR, "aspect ratio Λ = s²/A", np.arange(3.0, 20.1, 1.0), unit="",
                       xlabel="spanwise position y [m]", ylabel="flight direction [m]", xrange=[-9, 9], yrange=[-2.6, 1.2],
                       title="Same area, stretched: Λ = s²/A grows with the span squared")   # 18 slider positions
figF1b.show()                                                # display
print(pd.DataFrame({"Λ": [3, 6, 9.6, 20], "span s [m]": [np.sqrt(L_*15.0) for L_ in (3, 6, 9.6, 20)], "chord [m]": [15.0/np.sqrt(L_*15.0) for L_ in (3, 6, 9.6, 20)]}).round(2).to_string(index=False))   # the numbers behind four slider positions
""", explain=r"""
1. `planform_taper` builds the outline of a wing of fixed span and area for one taper ratio; `slider_figure` calls it once
   per slider position (the figure also works on the web page — nothing is computed when you drag).
2. `planform_AR` stretches a rectangle at fixed area: from $\Lambda\equiv s^2/A$ (14.1) the span is $\sqrt{\Lambda A}$.
3. The printed table lists span and chord for four aspect ratios.""")
nb.figure_notes(
    see="Two slider figures. First: a wing whose tip chord shrinks while span and area are held. Second: a rectangle of fixed "
        "area stretched from stubby (Λ = 3) to very slender (Λ = 20).",
    read="In the first, the outline changes but Λ = 9.6 does not: the aspect ratio does not see taper. In the second, going "
         "from Λ = 3 to Λ = 20 multiplies the span by √(20/3) ≈ 2.6 and divides the chord by the same factor.",
    change="…you doubled the span at fixed area? Λ quadruples. C14 shows the price of lift falls with the square of the span.")
whatif(r"""…you only knew the section's shape and the angle α — could you predict $C_L$? Not yet: ideal flow has one free number,
the circulation. C02 is about who sets it.""")

# =====================================================================================================================
# A.3 §14.3 Characteristics of Airfoil Sections — C02, C03
# =====================================================================================================================
nb.section("14.3", "Characteristics of Airfoil Sections", intro=r"""
**What is this section about?** Why a section lifts at all. Ideal flow allows any circulation; the sharp trailing edge picks
one (the Kutta condition); and the circulation it picks is paid for by a vortex left behind at the start. The section's
pressure distribution — a §14.3 picture — is computed after the lift, in C07.""")
core("C02", r"The Kutta condition: $\left.\dfrac{dw}{d\zeta}\right\rvert_{\zeta=b}=0$",
     "Ideal flow round a wing has a solution for every circulation. Which one does a real wing choose?")
problem(r"""
Hold a flat board at a slight angle in a stream. Chapter 6 says the force on it is $L=\rho U\Gamma$ (6.62) — and that the equations
are satisfied for *any* Γ, including zero. A theory that cannot say whether a wing lifts is not much use. The missing fact is
about the sharp rear edge: with the wrong circulation the air would have to whip round that edge from one side to the other,
at a speed that grows without limit as the edge gets sharper. Real air cannot do that. It leaves the edge smoothly — and
exactly one circulation lets it.""")
idea(r"""
   rear end of a plate at incidence, three circulations (flow from the left, plate rising to the left)

   Γ = 0                       Γ too small                  Γ = Γ_K  (Kutta)
   ───────●──╮                 ─────────●─╮                 ───────────●
     stagnation  ↺ flow turns     point closer,               flow leaves the
     point on top  the edge       still turning               edge smoothly
""", r"""
● is the rear **stagnation point** (where the fluid speed is zero). Adding clockwise circulation pushes it toward the edge;
the **Kutta condition** says the flow settles where it sits exactly *on* the edge.""")
remind("C02")
recap("R04", "The Kutta–Zhukhovsky lift theorem", r"""
$L=\rho U\Gamma$ (6.62): a two-dimensional body in a steady ideal stream feels a force per unit span of $\rho U\Gamma$ at right angles
to the stream, whatever its shape. $\Gamma$ is the **circulation**, the loop integral of the velocity round the body,
$\oint\mathbf u\cdot d\mathbf x$ (Ch. 3). **New here:** the theorem leaves Γ free; this chapter supplies it.

| Convention | Circulation of a lifting wing | Lift |
|---|---|---|
| the book (this chapter) | $\Gamma>0$, clockwise | $L=+\rho U\Gamma$ |
| Ch. 3 and Ch. 5 of this project | $\Gamma_{ccw}<0$, counterclockwise-positive | $L=-\rho U\Gamma_{ccw}$ |

> ⚠️ **Common confusion (trap T1):** Γ is **clockwise-positive** in this chapter. A wing flying to the left with the air
> faster over its top has clockwise circulation, Γ > 0, and positive lift.""", where="Ch. 6 §6.5", core_id="C02")
nb.code(r"""
print(ch06.lift_per_span(rho=1.225, U=30.0, Gamma_cw=8.214))      # the book's sign: clockwise 8.214 m²/s -> L = +301.9 N/m
print(ch06.lift_per_span(rho=1.225, U=30.0, Gamma_ccw=-8.214))    # the same flow, counted counterclockwise: Γ_ccw = −8.214 -> the same +301.9 N/m
""")
maths(r"""
1. **A flat plate is a mapped circle.** A flat plate of chord $4b$ is the image of the circle of radius $b$ under the map
   $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3) (recap R06, in C04); the trailing edge is the image of the point $\zeta=b$. So we may look for
   the stagnation points on the circle, where Chapter 6 gave them. (Reminder from Ch. 6: the **complex potential**
   $w=\phi+i\psi$ packs the velocity potential and the stream function into one function; streamlines are the lines of constant
   $\psi$.)
2. **Stagnation points on the circle.** They sit where the surface speed $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (6.37) vanishes:
   $\sin\theta=-\dfrac{\Gamma}{4\pi aU}$ (6.38), with θ measured from the stream direction.
3. **Where is the edge?** With the stream tilted up by α, the point $\zeta=b$ is at $\theta=-\alpha$.
4. **The condition.** The rear stagnation point reaches the edge when the right-hand side of step 2 equals $\sin(-\alpha)$. We
   find that Γ numerically first; C06 turns it into a formula.
5. **In one line.** Stated on the complex velocity, the Kutta condition is $\left.\dfrac{dw}{d\zeta}\right\rvert_{\zeta=b}=0$: the circle's
   flow stagnates at the point that becomes the trailing edge.""")
note("N10 [B]", r"""
Building lift means building circulation. In a real flow the circulation is the net vorticity held in the two boundary layers:
$\Gamma=\oint u_s\,ds=\int_0^c(U_{\text{upper}}-U_{\text{lower}})\,dx$ (ours) — the air is faster over the top. (Worked with numbers in C07
below.)""")
nb.worked_example("three circulations on a unit circle", r"""
$U$ = 1 m/s, $a=b$ = 1 m, α = 30°, so the edge is at θ = −30°. From $\sin\theta=-\dfrac{\Gamma}{4\pi aU}$ (6.38):

1. Γ = 0: sin θ = 0, the rear stagnation point is at θ = 0° — 30° round from the edge, on the upper side.
2. Γ = π m²/s: sin θ = −π/(4π) = −0.25, θ = −14.5° — half-way there.
3. Γ = 2π m²/s: sin θ = −0.5, θ = −30° — on the edge.
4. Γ = 3π m²/s: sin θ = −0.75, θ = −48.6° — past the edge, on the lower side.

Only the third leaves the edge alone.""")
nb.code(r"""
U, b, alpha = sec["U"], sec["b"], sec["alpha"]              # speed [m/s], map constant [m] (plate chord 4b = 1 m), angle of attack [rad] (5°)
GK_plate = AF.kutta_circulation_root(U=U, b=b, e=0.0, beta=0.0, alpha=alpha)  # the clockwise Γ that puts the stagnation point on the edge, found by a root-finder [m²/s]
print(f"Kutta circulation of the plate (clockwise): {GK_plate:.3f} m²/s")   # 8.214
for f in (0.0, 0.5, 1.0, 1.2):                              # four trial circulations as fractions of the Kutta value
    rsp = AF.rear_stagnation_point(U=U, a=b, alpha=alpha, beta=0.0, Gamma_cw=f*GK_plate)   # where the rear stagnation point sits on the circle
    Ff = AF.zhukhovsky_flow(U=U, b=b, e=0.0, beta=0.0, alpha=alpha, Gamma_cw=f*GK_plate)   # the plate's flow with this circulation
    q4, q6 = Ff.speed_near_edge(d=1e-4)/U, Ff.speed_near_edge(d=1e-6)/U                    # fluid speed 0.1 mm and 0.001 mm behind the edge, in units of U
    print(f"Γ/Γ_K = {f:3.1f}: stagnation point {np.degrees(rsp['offset']):+5.2f}° from the edge ({rsp['offset']:+.4f} rad); speed at 0.1 mm = {q4:5.2f} U, at 0.001 mm = {q6:6.2f} U")   # angle in both units; two distances
""", explain=r"""
1. `AF.kutta_circulation_root` knows no formula: it searches for the Γ at which the circle's surface speed at $\zeta=b$ is zero.
   For our plate (chord 1 m, 30 m/s, α = 5°) it returns a clockwise 8.214 m²/s.
2. `AF.rear_stagnation_point` reports how far round the circle the stagnation point is from the edge: 5.0°, 2.5°, 0.0° and
   −1.0° for Γ/Γ_K = 0, 0.5, 1 and 1.2 — it moves linearly at these small angles.
3. `speed_near_edge(d)` is the fluid speed at the distance d behind the edge. Go a hundred times closer (0.1 mm → 0.001 mm)
   and the speed of every row except Γ/Γ_K = 1 grows about tenfold: it rises like $d^{-1/2}$ without limit. Only the Kutta row
   stays at about $U$ however close you look.
4. Past the Kutta value the point is on the *lower* side and the corner flow returns, turning the other way.""")
scratch(r"""
def edge_slip(G):                                            # tangential speed of the circle's flow at the point ζ = b, for a clockwise circulation G
    w_, dw_ = CONF.circle_flow_zeta(U=U, a=b, Gamma_cw=G, alpha=alpha)   # complex potential and complex velocity of the circle flow (Ch. 6)
    return (dw_(b + 0j)*1j).real                             # at ζ = b the tangent points along +i, so u_θ = Re[(u − iv)·i]
lo, hi = 0.0, 4*np.pi*U*b                                    # bracket: no circulation … the largest that keeps stagnation points on the circle [m²/s]
for _ in range(40):                                          # forty halvings of the bracket (bisection)
    mid = 0.5*(lo + hi)                                      # middle of the bracket
    if edge_slip(lo)*edge_slip(mid) <= 0:                    # sign change in the lower half?
        hi = mid                                             # yes: the root is there
    else:                                                    # no:
        lo = mid                                             # the root is in the upper half
G_mine = 0.5*(lo + hi)                                       # the root: the Kutta circulation [m²/s]
assert np.isclose(G_mine, AF.kutta_circulation_root(U=U, b=b, e=0.0, beta=0.0, alpha=alpha), rtol=1e-9)   # same number → the library does exactly this
print(f"bisection: Γ_K = {G_mine:.4f} m²/s (clockwise)")     # 8.2142
""", r"""
1. `edge_slip(G)` evaluates the circle flow of Chapter 6 at the one point that becomes the trailing edge and returns the
   speed along the circle there.
2. Forty halvings of the bracket [0, 4πUb] find the circulation at which that speed is zero.
3. The assert shows the tested function returns the same 8.2142 m²/s.""")
nb.animation(r"""
al15 = np.deg2rad(15.0)                                      # a larger angle of attack (15° = 0.262 rad) so that the motion is easy to see
GK15 = AF.kutta_circulation_root(U=U, b=b, e=0.0, beta=0.0, alpha=al15)   # the Kutta circulation of the plate at 15° [m²/s], clockwise
nx, ny = (120, 80) if not FAST else (80, 54)                 # grid for the stream function (smaller when FAST)
X, Y = np.meshgrid(np.linspace(-0.9, 0.9, nx), np.linspace(-0.6, 0.6, ny))   # points of the x–y plane [m]
levels = np.linspace(-0.55*U, 0.55*U, 23)                    # fixed stream-function levels [m²/s]; ψ = 0 is the plate itself
fractions = np.linspace(0.0, 1.2, 13) if not FAST else np.linspace(0.0, 1.2, 7)   # Γ/Γ_K of each frame: 0, 0.1, … 1.2
fig, ax = plt.subplots(figsize=(5.6, 3.6))                   # one panel
def update(i):                                               # draw frame i
    f = fractions[i]                                         # Γ/Γ_K of this frame
    Fp = AF.flat_plate_flow(U=U, b=b, alpha=al15, Gamma_cw=f*GK15)     # the plate's flow with this circulation
    st = AF.stagnation_points(U=U, b=b, e=0.0, beta=0.0, alpha=al15, Gamma_cw=f*GK15)   # where the two stagnation points are
    ax.clear()                                               # wipe the previous frame
    ax.contour(X, Y, Fp.stream(x=X, y=Y), levels=levels, colors=COLORS["blue"], linewidths=0.8)   # streamlines = lines of constant ψ
    ax.plot([-2*b, 2*b], [0, 0], color="k", lw=3)            # the plate: from x = −2b to x = +2b
    ax.plot(st["z_rear"].real, st["z_rear"].imag, "o", color=COLORS["orange"], ms=9, zorder=5)   # rear stagnation point (amber dot)
    ax.set(xlim=(-0.9, 0.9), ylim=(-0.6, 0.6), xlabel="x [m]", ylabel="y [m]", title=f"Γ/Γ_K = {f:4.2f} (Γ clockwise)")   # fixed axes
    ax.set_aspect("equal")                                   # true proportions
show_animation(animate(update, frames=len(fractions), fig=fig, interval=400), player="frames")   # a frame player: step through, stop at 0, 0.5, 1.0
""", explain=r"""
1. The plate is set at 15° so that the stagnation point starts well away from the edge.
2. `update(i)` builds the flow for one value of Γ/Γ_K (0, 0.1, … 1.2), draws its streamlines (contours of the stream
   function `Fp.stream`), the plate and the rear stagnation point.
3. `show_animation(..., player="frames")` gives ◀ ▮▮ ▶ buttons, so you can stop at Γ/Γ_K = 0, 0.5 and 1.""")
nb.figure_notes(
    see="Streamlines round a flat plate at 15° while the clockwise circulation grows from zero through the Kutta value to "
        "1.2 times it. The amber dot is the rear stagnation point.",
    read="At Γ = 0 the dot sits on the upper surface ahead of the edge and the streamlines wrap round the trailing edge — the "
         "flow turns a knife-sharp corner. As Γ grows the dot slides back. At Γ/Γ_K = 1 it is on the edge and the streamline "
         "from it leaves smoothly, straight off the plate. Watch also the streamlines far upstream: they are deflected "
         "upward more and more — that is the lift appearing.",
    change="…Γ goes past the Kutta value (the last two frames)? The dot passes the edge onto the lower side and the corner "
           "flow returns, the other way round.")
nb.figure(r"""
fr = np.linspace(0.0, 2.0, 201)                              # trial circulations as fractions of the Kutta value
fig, ax = plt.subplots(figsize=(6.4, 3.4))                   # one panel
for d_, col in ((1e-4, COLORS["blue"]), (1e-6, COLORS["accent"])):   # two distances from the edge: 0.1 mm and 0.001 mm
    q_edge = [AF.zhukhovsky_flow(U=U, b=b, e=0.0, beta=0.0, alpha=alpha, Gamma_cw=f*GK_plate).speed_near_edge(d=d_)/U for f in fr]   # speed there, in units of U
    ax.semilogy(fr, q_edge, color=col, label=f"{d_*1e3:g} mm from the edge")   # logarithmic vertical axis: each step up is a factor of ten
ax.axvline(1.0, color=COLORS["orange"], ls="--", lw=1)       # the Kutta value
ax.set(xlabel=r"circulation $\Gamma/\Gamma_K$ (clockwise) [–]", ylabel="speed near the edge, $q/U$ [–]", title="Only one circulation keeps the edge speed finite")   # labels
ax.legend(); plt.show()                                      # legend and display
""",
    see="The fluid speed very close to the trailing edge against the circulation, at two distances from the edge, on a "
        "logarithmic axis.",
    read="On a logarithmic axis each step up is a factor of ten. Each curve is a V with its notch at Γ/Γ_K = 1. Moving a "
         "hundred times closer to the edge lifts both arms of the V about tenfold (the speed grows like the inverse square "
         "root of the distance) — but the notch stays near $q=U$. Exactly one circulation gives a finite speed at the edge.",
    change="…you measured closer still? The arms keep rising without limit; the notch does not move.")
note("N11 [C]", r"""
Three ways to make and control lift: a sharp trailing edge (which lets the other two work), the angle of attack, and the
camber. Angle and camber become the sum α + β in $C_L\simeq2\pi(\alpha+\beta)$ (14.12), C06.""")
confusion(r"""(trap T14) The Kutta condition is **not a theorem of ideal flow**. It is a rule that stands in for what the
boundary layers do, and it fails at stall and for motions faster than the time the air takes to cross the chord.""")
whatif(r"""…the foil starts from rest? Then Γ = 0 at first — where does the circulation come from? C03.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C03", r"The starting vortex: $\Gamma_{ABCD}=0=\Gamma_{\text{foil}}+\Gamma_{\text{starting}}$",
     "Kelvin says circulation cannot be created in an inviscid fluid. So where does the wing's come from?")
problem(r"""
Stir a spoon sideways through coffee and stop: a little whirl is left spinning where you began. A wing does the same. In the
first moments of its motion it throws off a vortex from its trailing edge, and the air round the wing is left circulating the
other way. Nothing was created — the two cancel — but the wing flies away from one half of the pair and keeps the other.""")
idea(r"""
     A ─────────────────────────B───────────────────────── (far away)
     │                          ┊                         │
     │      ↻ foil              ┊        ↺ starting       │      total circulation round the big box  = 0
     │   (moving left)          ┊          vortex         │      = (+Γ round the foil) + (−Γ round the starting vortex)
     │                          ┊                         │
     D ─────────────────────────┴───────────────────────── C
                              cut BD
""", r"""
A big closed curve that surrounds everything had zero circulation before the start, so it has zero circulation for ever.
Cut it in two between the wing and what it left behind: the two halves must be equal and opposite.""")
recap_md("Kelvin's circulation theorem", "Ch. 5 §5.2", r"""
$\dfrac{D\Gamma}{Dt}=0$ (5.8): the circulation round a closed curve that moves with the fluid (a *material circuit*) does not change,
provided the curve itself stays in inviscid fluid of uniform density.""")
recap_md("circulation adds over sub-circuits", "Ch. 3; Stokes' theorem, Ch. 2", r"""
Cut a circuit in two with a line: the two circulations add up to the original, because the cut is walked once each way and
its two contributions cancel, $\Gamma_{\text{whole}}=\Gamma_{\text{left}}+\Gamma_{\text{right}}$.""")
note("N14 [B]", r"""
An impulsively started foil at three times: $\Gamma=0\ \to\ 0<\Gamma<\Gamma_{\text{Kutta}}\ \to\ \Gamma=\Gamma_{\text{Kutta}}$ for $t>c/U$ — the steady
pattern is reached after the foil has moved a chord or two.""")
nb.figure(r"""
fig = fig_three_stages(**SEC)   # streamlines of our section at Γ/Γ_K = 0, ½, 1 with both stagnation points
plt.show()                                                   # display
""",
    see="Three streamline pictures of our section at 5°: no circulation, half the Kutta value, the Kutta value. The dots are "
        "the front and rear stagnation points.",
    read="Left: the rear dot is on the upper surface and the flow turns the sharp tail. Middle: it has moved half-way. Right: "
         "it is at the tail and the flow leaves smoothly. The front dot moves too — under the nose — as the circulation grows.",
    change="…the section were symmetric and at α = 0? All three panels would be the same: no circulation is needed.")
note("N15 [B]", r"""
**The mechanism.** At the very start the flow does turn the sharp corner. But the fluid next to the surface is slow (the
boundary layer) and cannot follow so tight a turn against rising pressure: it separates at the edge, rolls up, and the shed
layer pushes the stagnation point back to the edge. A circulation larger than the Kutta value is allowed by ideal flow but is
not observed (trap T14).""")
nb.code(r"""
fig = draw_kelvin_circuit()                                  # the material circuit ABCD round the foil and its starting vortex, cut by BD
plt.show()                                                   # display the sketch used by the derivation below
""", tags=["figure"])
D("D02", "The starting vortex carries the opposite circulation", ref="5.8")
P("P339", "discrete-vortex time stepping", r"""
A cheap model of a wing that sheds vorticity: one "bound" **point vortex** (a swirl whose speed is Γ/2πr at distance r, Ch. 5)
stands for the wing; at every time step one new vortex is released just behind the trailing edge; its strength is whatever
keeps the total circulation zero while the air does not flow through the plate at one control point. The total is the
invariant to watch.""", r"""
G_bound, G_wake = 0.0, []                 # nothing circulates at the start
G_bound = 5.0; G_wake.append(-5.0)        # step 1: the wing gains 5, sheds -5
G_bound = 6.5; G_wake.append(-1.5)        # step 2: gains 1.5 more, sheds -1.5
print(G_bound + sum(G_wake))              # 0.0 at every step
""")
nb.worked_example("how strong is the starting vortex?", r"""
A flat plate of chord $c$ = 1 m at α = 30° reaches $U$ = 10 m/s.

1. Steady circulation (clockwise): $\pi Uc\sin\alpha$ = π × 10 × 1 × 0.5 = 15.7 m²/s (the plate's value of
   $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11), with $a=c/4$ and β = 0 — derived in C06).
2. Starting vortex: −15.7 m²/s, counterclockwise.
3. Lift: $L=\rho U\Gamma$ (6.62) = 1.2 × 10 × 15.7 = 188 N per metre of span.
4. The speed is then doubled to 20 m/s: the wing needs 31.4 m²/s, so it sheds another −15.7.
5. The wing stops: the bound 31.4 m²/s is shed as a +31.4 "stopping vortex". Total at every stage: 0.""")
nb.code(r"""
run = ch14.impulsive_start_lumped_vortex(U=sec["U"], c=1.0, alpha=sec["alpha"], t_end=20/sec["U"], dt=0.05/sec["U"], rho=1.0)   # plate of chord 1 m started at 30 m/s, α = 5°; 400 steps of 0.05 chords, 20 chords travelled
print(f"steady circulation πUc sin α = {run['Gamma_steady'][-1]:.3f} m²/s (clockwise)")      # 8.214, the Kutta value of C02
print("bound/steady after 1, 3, 5, 10 chords:", np.round(run["ratio"][[19, 59, 99, 199]], 3))   # 0.646, 0.816, 0.879, 0.939
print(f"largest |bound + shed| over the run = {np.abs(run['total']).max():.1e} m²/s")       # zero to round-off: Kelvin's theorem
""", explain=r"""
1. `ch14.impulsive_start_lumped_vortex` is the model of the primer: a bound vortex at the quarter chord, a control point at
   three-quarter chord, one new wake vortex per step, total circulation kept at zero.
2. The steady value it approaches is the Kutta circulation of C02, 8.214 m²/s clockwise.
3. The bound circulation is 65 % of that after one chord of travel, 82 % after three, 88 % after five and 94 % after ten:
   about half the lift arrives almost at once, the rest slowly — the starting vortex is still close and its velocity field
   holds the circulation down.
4. **Not in the book — ours:** this is a lumped-vortex sketch, not a simulation of the boundary layer.""")
scratch(r"""
U_, c_, al_ = sec["U"], 1.0, sec["alpha"]                    # speed [m/s], chord [m], angle of attack [rad]
dt = 0.05*c_/U_                                              # time step: the plate travels 0.05 chords per step [s]
xP = 0.75*c_                                                 # control point at three-quarter chord (the bound vortex is at c/4)
wx, wy, wG = [], [], []                                      # wake vortices: positions [m] and clockwise strengths [m²/s]
G_mine = []                                                  # history of the bound circulation
for k in range(400):                                         # 400 steps = 20 chords
    wx = [x + U_*np.cos(al_)*dt for x in wx]                 # (i) every wake vortex drifts with the stream, x-part
    wy = [y + U_*np.sin(al_)*dt for y in wy]                 #     … and y-part
    xn, yn = c_ + 0.25*U_*dt*np.cos(al_), 0.25*U_*dt*np.sin(al_)   # (ii) the new vortex is placed a quarter-step behind the trailing edge
    k_new = -(xP - xn)/(2*np.pi*((xP - xn)**2 + yn**2))      # y-velocity at the control point per unit clockwise strength of the new vortex
    W = sum(G*(-(xP - x))/(2*np.pi*((xP - x)**2 + y**2)) for G, x, y in zip(wG, wx, wy))   # y-velocity there from the older wake vortices
    S = sum(wG)                                              # circulation already in the wake
    G_new = (-U_*np.sin(al_) - S/(np.pi*c_) - W)/(1/(np.pi*c_) + k_new)   # (iii) no flow through the plate at P, with Γ_bound = −Γ_new − S
    G_b = -G_new - S                                         # total circulation zero: the bound vortex balances the whole wake
    wx.append(xn); wy.append(yn); wG.append(G_new)           # the new vortex joins the wake
    G_mine.append(G_b)                                       # record
assert np.allclose(G_mine, run["Gamma_bound_cw"])            # same numbers → the library does exactly this
print(f"bound circulation after 1 and 10 chords: {G_mine[19]:.3f}, {G_mine[199]:.3f} m²/s")   # 5.309 and 7.711
""", r"""
1. Each step moves the wake with the stream, places one new vortex just behind the edge, and solves two equations for two
   unknowns: *no flow through the plate at the control point* and *total circulation zero*.
2. A clockwise vortex of strength Γ at distance r induces the speed Γ/(2πr); its vertical part at the control point is the
   factor `k_new` (and the sum `W` for the older vortices). The bound vortex, half a chord ahead of the control point,
   contributes −Γ_b/(πc).
3. The assert confirms the tested function is this loop.""")
nb.animation(r"""
hist = ch14.impulsive_start_lumped_vortex(U=sec["U"], c=1.0, alpha=np.deg2rad(10.0), t_end=6/sec["U"], dt=0.1/sec["U"], history=True)   # a shorter run at 10° (0.175 rad) with the wake positions kept at every step
n_frames = 60 if not FAST else 30                            # number of frames
idx = np.linspace(0, len(hist["t"]) - 1, n_frames).astype(int)   # which steps to show
Gs = hist["Gamma_steady"][0]                                 # the steady (Kutta) circulation [m²/s], clockwise
sizes = 4 + 60*np.abs(hist["wake_Gamma_cw"])/np.abs(hist["wake_Gamma_cw"]).max()   # dot size of each wake vortex ∝ its strength
fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 3.4))     # left: the plate and its wake; right: the circulations
axL.plot([0, 1], [0, 0], color="k", lw=3)                    # the plate, chord 1 m, in its own frame: the stream comes from the left
dots = axL.scatter([], [], s=[], color=COLORS["orange"])     # the wake vortices (filled in by update)
ring = plt.Circle((0.25, 0), 0.12, fill=False, color=COLORS["teal"], lw=1)   # the bound vortex at the quarter chord, drawn as a ring
axL.add_patch(ring)                                          # put the ring on the panel
axL.set(xlim=(-0.5, 7.5), ylim=(-1.2, 1.6), xlabel="x [chords]", ylabel="y [chords]")   # fixed axes
axL.set_aspect("equal")                                      # true proportions
(lb,) = axR.plot([], [], color=COLORS["teal"], label="bound Γ (clockwise)")            # the wing's circulation
(lw_,) = axR.plot([], [], color=COLORS["rose"], label="total shed (clockwise-positive, so negative)")   # everything in the wake
(lt,) = axR.plot([], [], color=COLORS["orange"], label="their sum")                    # flat at zero
axR.axhline(Gs, color="0.6", ls="--", lw=1)                  # the steady value
axR.set(xlim=(0, 6), ylim=(-1.15*Gs, 1.15*Gs), xlabel="chords travelled", ylabel="circulation [m²/s]")   # fixed axes
axR.legend(loc="center right", fontsize=7)                   # which curve is which
fig.canvas.draw(); fig.set_layout_engine("none")             # lay the figure out once, then freeze it: frames render much faster
def update(j):                                               # draw frame j
    k = idx[j]                                               # step shown in this frame
    dots.set_offsets(np.c_[hist["wake_x_hist"][k, :k+1], hist["wake_y_hist"][k, :k+1]])   # positions of the vortices shed so far
    dots.set_sizes(sizes[:k+1])                              # their sizes
    ring.set_linewidth(1 + 3*hist["ratio"][k])               # the ring thickens as the bound circulation grows
    axL.set_title(f"{hist['s_chords'][k]:.1f} chords travelled")   # how far the plate has gone
    lb.set_data(hist["s_chords"][:k+1], hist["Gamma_bound_cw"][:k+1])   # bound circulation so far
    lw_.set_data(hist["s_chords"][:k+1], hist["Gamma_wake_cw"][:k+1])   # shed circulation so far
    lt.set_data(hist["s_chords"][:k+1], hist["total"][:k+1])            # their sum so far
show_animation(animate(update, frames=n_frames, fig=fig, interval=80), player="video")   # a smooth video
""", explain=r"""
1. `history=True` keeps the position of every wake vortex at every step (`wake_x_hist[k, :k+1]` are those present at step k).
2. Left panel: the plate, the wake vortices (amber dots, sized by strength) drifting downstream, and the bound vortex as a
   teal ring that thickens as its circulation grows.
3. Right panel: bound circulation (teal), total shed circulation (rose) and their sum (amber).
4. `update(j)` only moves the existing dots and extends the three curves, which keeps the animation quick to build.""")
nb.figure_notes(
    see="A plate started impulsively at 10°. Left: the vortices it sheds drift downstream; the first one is by far the "
        "strongest. Right: the bound circulation climbs toward the dashed steady value while the shed total falls as its "
        "mirror image.",
    read="The amber line never leaves zero: at every instant what the wing has gained equals what the wake holds, with the "
         "opposite sign. The big first dot is the *starting vortex*; the later ones are weak corrections.",
    change="…α were doubled half-way through? A second strong vortex of the same sign as the first would be shed, and the "
           "teal curve would step up to twice its value.")
nb.figure(r"""
s_ch = run["s_chords"]                                       # distance travelled, in chords
lift_ratio = run["lift_circulatory"]/(1.0*sec["U"]*run["Gamma_steady"][0])   # the model's quasi-steady lift ρUΓ_bound over the steady lift ρUΓ_steady (ρ = 1 in this run)
fig, (axG, axL) = plt.subplots(1, 2, figsize=(10.0, 3.6))    # left: circulations; right: lift ratio
axG.plot(s_ch, run["Gamma_bound_cw"], color=COLORS["teal"], label="bound circulation (clockwise) — lumped-vortex model, ours")   # the wing
axG.plot(s_ch, run["Gamma_wake_cw"], color=COLORS["rose"], label="total shed circulation")        # the wake
axG.plot(s_ch, run["total"], color=COLORS["orange"], label="their sum (Kelvin: zero)")            # the invariant
axG.axhline(run["Gamma_steady"][0], color="0.6", ls="--", lw=1, label="steady Kutta value πUc sin α")   # the only reference on this panel
axG.set(xlabel="distance travelled [chords]", ylabel="circulation [m²/s]", title="The wing's circulation is paid for by the wake")   # labels
axG.legend(fontsize=7)                                       # which curve is which
axL.plot(s_ch, lift_ratio, color=COLORS["teal"], label="model: ρUΓ_bound / steady lift (quasi-steady part only)")   # our model's lift ratio
axL.plot(s_ch, ch14.wagner_function(s=2*s_ch), color="0.45", ls=":", lw=2, label="Wagner's function — circulatory lift / steady lift (Jones' fit)")   # the argument is in semi-chords = 2 × chords
axL.axhline(1.0, color="0.6", ls="--", lw=1)                 # the steady lift
axL.set(xlabel="distance travelled [chords]", ylabel="lift / steady lift [–]", ylim=(0, 1.08), title="Lift builds up over several chords")   # labels
axL.legend(fontsize=7, loc="lower right"); plt.show()      # legend and display
s90 = np.interp(0.9, run["ratio"], s_ch)                      # distance at which the bound circulation reaches 90 % of the steady value (the ratio rises steadily, so it can be read backwards)
print(f"first value of Γ_bound/Γ_steady = {run['ratio'][0]:.3f}; 90 % reached after {s90:.2f} chords")   # about 1/3 (0.344 with this step) and about 6.2 chords
print(f"model lift ratio minus Wagner's function at 1, 2, 5, 10 chords: {np.round((lift_ratio - ch14.wagner_function(s=2*s_ch))[[19, 39, 99, 199]], 3)}")   # the gap closes after two chords
""",
    see="The 20-chord run of the code cell. Left: bound circulation, shed circulation, their sum, and the steady Kutta value "
        "(dashed) — the only reference on that panel. Right: the model's lift ratio $\\rho U\\Gamma_{\\text{bound}}$ over the steady "
        "lift, with Wagner's function (dotted).",
    read="Wagner's function is a **lift** ratio — circulatory lift over steady lift, here in a published two-exponential fit "
         "(Jones' fit; its argument is in semi-chords) — so it belongs on the right panel, not on the circulation panel. The "
         "model's first $\\Gamma_{\\text{bound}}$ is 1/3 of the steady value (0.34 with this time step); in the exact thin-airfoil solution the circulation "
         "starts from zero and it is the circulatory lift that starts at one half; this model's $\\rho U\\Gamma_{\\text{bound}}$ "
         "follows Wagner's lift function to about 1 % after two chords and runs above the exact circulation (0.65 against "
         "0.55 after one chord). `lift_circulatory` is the quasi-steady part $\\rho U\\Gamma_{\\text{bound}}$ only — the "
         "$\\rho\\,d\\Gamma/dt$ term is omitted. The model reaches 90 % of the steady value after about 6.2 chords.",
    change="…the plate were ten times longer? Both panels are the same in chords: the transient lasts ten times longer in metres.")
note("N16 [B]", r"""
**Viscosity is responsible for the drag — and for the lift.** Without the boundary layer at the edge no vortex would be shed,
no circulation would arise, and the wing would not lift. Once the circulation is set, the rest of the flow is nearly ideal.""")
note("N17 [C]", r"""
The idea of circulation as the source of lift goes back to Lanchester, Kutta and Zhukhovsky, working independently;
Lanchester and Prandtl return in N69 (C08).""")
slip(5, "a birth year for Lanchester that would make him a child when he founded his motor company", "1868.")
explainer("starting_vortex", "Where does a wing's circulation come from?",
          "a new vortex of the right sign rolls off the edge at every change of speed or incidence while the amber total "
          "stays flat at zero — something only time can show.",
          ["Press play on 'impulsive start' and watch the teal and rose curves mirror each other.",
           "Choose 'sudden stop': which way does the stopping vortex turn?",
           "Shrink the circuit until it cuts between foil and starting vortex and read the circulation it now encloses.",
           "Step α up half-way through a run and predict the sign of the new vortex first."])
whatif(r"""…we want Γ for a real section shape, not a plate? We need a family of shapes whose flow we can write down: C04.""")


# =====================================================================================================================
# A.4 §14.4 Conformal Transformation for Generating Airfoil Shapes — C04
# =====================================================================================================================
nb.section("14.4", "Conformal Transformation for Generating Airfoil Shapes", intro=r"""
**What is this section about?** A recipe for wing sections whose flow we can write down exactly: take a circle, move its
centre, and push it through one simple map. Where the circle sits decides the camber and the thickness.""")
core("C04", r"The Zhukhovsky airfoil: the circle $\lvert\zeta-\mu\rvert=a$ through $\zeta=b$, mapped by $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3)",
     "How can one map turn a circle into a flat plate, a curved plate, a fat symmetric fin or a real-looking wing section?")
problem(r"""
We know the flow round a circle (Ch. 6). If we can find a change of coordinates that turns a circle into a wing section and
leaves the far stream alone, the circle's flow *becomes* the wing's flow, for free. One such change is the Zhukhovsky map.
Chapter 6 used it on a circle centred at the origin and got a flat plate. The new move here is to slide the circle's centre: a
little up, a little back, always keeping the circle through the point $\zeta=b$. That one point is special — the map folds the
plane there — and it becomes the sharp tail.""")
nb.md(r"""
#### The idea

| Where the circle's centre is | What comes out | Controlled by |
|---|---|---|
| at the origin | a flat plate of length $4b$ | — |
| raised (on the imaginary axis) | a circular arc | $\beta$ (camber) |
| shifted back (on the negative real axis) | a symmetric section: round nose, sharp tail | $e$ (thickness) |
| both | a cambered airfoil | $\beta$ and $e$ |

Two knobs, two properties. In this block the letter $b$ is the constant of the map (the chord is about $4b$), $a$ is the
radius of the circle, $e$ is the thickness parameter and $\beta$ the camber angle of the circle.""")
remind("C04")
recap("R05", "The Zhukhovsky transformation", r"""
$z=\zeta+\dfrac{b^2}{\zeta}$ (14.3) — the same map as $z=\zeta+\dfrac{b^2}{\zeta}$ (6.65). It is a *conformal map*: it sends each point ζ of one
complex plane to a point $z$ of another and preserves angles between curves (except where its derivative vanishes). Far away
$z\approx\zeta$, so the stream is untouched; near the origin the $b^2/\zeta$ term dominates. **New here:** circles that are not centred
on the origin.""", where="Ch. 6 §6.6", core_id="C04")
nb.code(r"""
print(CONF.joukowski(zeta=2j, b=1.0))                        # the point 2i maps to 2i + 1/(2i) = 1.5i
print(CONF.joukowski_inverse(z=1.5j, b=1.0))                 # and back: the root outside the unit circle is 2i
""")
recap("R06", "The centred circle becomes a slit", r"""
For $\zeta=be^{i\theta}$: $z=be^{i\theta}+be^{-i\theta}=2b\cos\theta$ — real, running from $2b$ to $-2b$ and back: a flat plate of length $4b$,
traced twice (upper and lower side).""", where="Ch. 6 §6.6", core_id="C04")
note("N18 [C]", r"""
Other routes to the flow round a *given* section: solve Laplace's equation numerically (Ch. 6 panels, Ch. 10), **thin-airfoil
theory** (a vortex sheet on the camber line; reused in Ch. 15), or **Theodorsen's method** (map the section back to a circle).
We do the reverse: choose the map, accept the shape. One check line with thin-airfoil theory, for a parabolic camber line of
3.5 % camber:""")
nb.code(r"""
ta = AF.thin_airfoil(dyc_dx=lambda xc: 4*0.035*(1 - 2*xc), alpha=0.0)   # slope of a parabolic camber line with camber/chord = 0.035, at zero angle of attack
print(f"zero-lift angle = {ta['alpha_L0']:.4f} rad = {np.degrees(ta['alpha_L0']):.2f}°")   # ≈ −0.070 rad = −4.0° = −2 × camber/chord (N31 in C06 gets the same)
""")
note("N25 [B]", r"""
Shift the centre to the left along the real axis by $be$ and enlarge the radius slightly: $a\equiv b(1+e)$, $e\ll1$ (14.7). The circle
still passes through $\zeta=b$. **Not in the book — ours:** one parametrisation for all four cases, centre $\mu=b-ae^{-i\beta}$ and
radius $a=b(1+e)/\cos\beta$. Points of the circle are $\zeta=\mu+ae^{it}$ with $t=-\beta$ at the point that becomes the trailing edge.""")
confusion(r"""(trap T7) the book defines $a=b(1+e)$ only for the symmetric case and $a=b/\cos\beta$ only for the arc;
`AF.zhukhovsky_circle` gives both as limits and also accepts `a=` directly.""")
nb.code(r"""
print(AF.zhukhovsky_circle(b=1.0, e=0.1, beta=0.0))          # symmetric case: centre −0.1, radius 1.1 — the book's a = b(1 + e)
print(AF.zhukhovsky_circle(b=1.0, e=0.0, beta=np.arctan(0.75)))   # raised circle with tan β = 3/4: centre 0.75i, radius 1.25 = b/cos β
print(AF.zhukhovsky_circle(b=sec["b"], e=sec["e"], beta=sec["beta"]))   # our section: centre −0.0200 + 0.0189i m, radius 0.2707 m
""", explain=r"""
**What does this show?** One function covers the book's two special cases and the general one. Each returned pair is
(centre μ, radius a); every circle passes through $\zeta=b$ because $\mu+ae^{-i\beta}=b$.""")
nb.figure(r"""
fig = fig_four_circles()                                     # four circles (ζ-plane, above) and their images under z = ζ + b²/ζ (z-plane, below)
plt.show()                                                   # display
""",
    see="Top row: four circles in the ζ-plane, with the ghost circle of radius $b$ and the points $\\pm b$ marked. Bottom row: "
        "what the map $z=\\zeta+b^2/\\zeta$ makes of each — a flat plate, a circular arc, a symmetric section, a cambered section.",
    read="Every circle passes through $\\zeta=b$, and that point always lands on the tail $z=2b$. Raising the centre bends the "
         "plate; moving it back swallows the point $-b$ and rounds the nose.",
    change="…the circle missed $\\zeta=b$? The tail would be rounded and the Kutta condition would have nothing to hold on to.")
P("P340", "the cosine rule", r"""
In any triangle with sides $p$, $q$, $r$ and the angle $C$ between $p$ and $q$: $r^2=p^2+q^2-2pq\cos C$. With $C$ = 90° it is
Pythagoras.""", r"""
p, q, C = 3.0, 4.0, np.deg2rad(60)            # two sides and the angle between them (60° = 1.047 rad)
print(np.sqrt(p**2 + q**2 - 2*p*q*np.cos(C)))   # third side: 3.606
""")
P("P341", "cot, csc and the double-angle identities", r"""
$\cot x=\cos x/\sin x=1/\tan x$; $\csc x=1/\sin x$. Two facts we need: $1+\cot^2x=\csc^2x$ and $\cot\beta-\tan\beta=2\cot2\beta$.""", r"""
B = 0.3                                       # any angle [rad]
print(1/np.tan(B) - np.tan(B), 2/np.tan(2*B))   # both 2.923: cot β − tan β = 2 cot 2β
print(1 + 1/np.tan(B)**2, 1/np.sin(B)**2)     # both 11.45: 1 + cot² = csc²
""")
nb.code(r"""
fig = draw_triangle_OQP()                                    # the triangle O (origin) – Q (centre of the raised circle) – P (a point on it)
plt.show()                                                   # display the sketch used by steps 6–7 of the derivation below
""", tags=["figure"])
confusion(r"""(trap T4) in the derivation below θ is the polar angle of ζ seen from the **origin** O, as in
$x=(R+b^2/R)\cos\theta$, $y=(R-b^2/R)\sin\theta$ (14.4) — not the angle seen from the circle's centre Q.""")
D("D03", "A raised circle maps to a circular arc, and the camber rule", ref="14.6")
nb.code(r"""
print(ch14.circular_arc_sym()["ok"])                         # the module's own symbolic engine re-does (14.4) → (14.5) → (14.6) → the arc: True
arc = AF.circular_arc(b=sec["b"], beta=sec["beta"])          # the arc for our b = 0.25 m and β = 4°: centre, radius, height, end angle
S_arc = AF.zhukhovsky_airfoil(b=sec["b"], e=0.0, beta=sec["beta"])   # the mapped raised circle (no thickness)
res = np.abs(np.abs(S_arc["z"] - arc["centre"]) - arc["radius"]).max()   # largest distance of a mapped point from that circle [m]
print(f"height = {arc['height']:.4f} m, centre = {arc['centre'].imag:.3f}i m, radius = {arc['radius']:.3f} m, end angle = {np.degrees(arc['end_angle']):.1f}° ({arc['end_angle']:.4f} rad)")   # 0.0350, −3.558i, 3.593, 8°
print(f"largest distance of a mapped point from the arc's circle = {res:.1e} m")   # round-off: every mapped point lies on it
""", explain=r"""
1. `ch14.circular_arc_sym()` is a sympy engine that repeats the derivation's chain; `"ok"` is True when every residual is zero.
2. `AF.circular_arc` returns the circle of the result for our numbers: its centre is far below the chord (−3.558 m), its
   radius 3.593 m, and the arc rises 35 mm above a chord of 1 m. The tangent at each end makes the angle 2β = 8° with the chord.
3. The residual is at round-off level: every point of the mapped circle lies on that arc. The arc is traced twice — the lower
   half of the circle lands on the same arc.""")
P("P342", "binomial approximation with a complex small quantity", r"""
For a small number ε — real or complex — $\dfrac1{1+\epsilon}\approx1-\epsilon$, with an error of about ε². "First order in $e$" means:
keep $e$, drop $e^2$ (orders of smallness, primed in Ch. 2, P68).""", r"""
eps = 0.05*(1 - np.exp(-1j*2.0))              # a small complex number, size 0.08
print(1/(1 + eps), 1 - eps)                   # agree to about 0.007 = abs(eps)**2
""")
D("D04", "A circle shifted back gives a symmetric section: chord, thickest point and thickness rule", ref="14.8")
note("N26 [B]", r"""
The three shape rules: $c\simeq4b$, camber $\simeq\tfrac12\beta c$, $t_{\max}/c\simeq1.3\,e$ (14.8), with the thickest point near the quarter chord,
$x=-b$.""")
slip(4, r"""“camber = ≅ ½βc” with two relation signs, and the coefficient 1.3 in $t_{\max}/c\cong1.3\,e$ (14.8)""",
     r"""camber $\simeq\tfrac12\beta c$, and the coefficient is $3\sqrt3/4=1.299$ (derived above) — a first-order rule in $e$.""")
nb.code(r"""
rows = []                                                    # one row per thickness parameter
for e_ in (0.02, 0.05, 0.10, 0.25):                          # from thin to very fat
    meas = AF.zhukhovsky_shape(b=1.0, e=e_, beta=0.0, exact=True)     # measured on the exact contour
    rule = AF.zhukhovsky_shape(b=1.0, e=e_, beta=0.0, exact=False)    # the first-order rules of (14.8)
    rows.append({"e": e_, "t/c measured": meas["thickness_ratio"], "t/c rule 1.299 e": rule["thickness_ratio"],
                 "rule too high by [%]": 100*(rule["thickness_ratio"]/meas["thickness_ratio"] - 1), "chord/4b": meas["chord"]/4.0})   # collect
print(pd.DataFrame(rows).round(4).to_string(index=False))    # the approximation made measurable
cam = [(np.degrees(B_), AF.zhukhovsky_shape(b=1.0, e=0.0, beta=B_, exact=True)["camber_ratio"], B_/2) for B_ in np.deg2rad([2.0, 4.0, 6.0, 12.0])]   # camber/chord measured against β/2
print(pd.DataFrame(cam, columns=["β [°]", "camber/c measured", "rule β/2"]).round(4).to_string(index=False))   # the camber rule
""", explain=r"""
**What does this show?** The thickness rule is a first-order statement: it overshoots the measured thickness ratio by 2 %,
5 %, 10 % and 26 % at $e$ = 0.02, 0.05, 0.10 and 0.25 (the fourth column; put the other way round, the measured value is
2 %, 5 %, 9 % and 21 % below the rule). The chord is within 1 % of $4b$ up to $e$ = 0.1. The camber rule holds within 1 % up
to β = 6°.""")
nb.worked_example("a raised circle by hand", r"""
$b$ = 1, tan β = 3/4 (β ≈ 36.9°; far from thin, but the numbers are easy).

1. Radius $a=b/\cos\beta$ = 1.25; centre Q at 0.75i.
2. Top of the circle: ζ = 0.75i + 1.25i = 2i. Its image under $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3): 2i + 1/(2i) = 2i − 0.5i = 1.5i.
3. Bottom of the circle: ζ = −0.5i → z = −0.5i + 1/(−0.5i) = −0.5i + 2i = 1.5i — the *same* point: the arc is covered twice.
4. Ends: ζ = ±1 → z = ±2.
5. The arc's circle: centre at −2b cot 2β = −0.583, radius 2b csc 2β = 2.083; check the top: −0.583 + 2.083 = 1.5 ✓; check an
   end: √(2² + 0.583²) = 2.083 ✓.
6. Height 2b tan β = 1.5 ✓.""")
nb.code(r"""
S = AF.zhukhovsky_airfoil(b=sec["b"], e=sec["e"], beta=sec["beta"])            # our section: the circle's points and their images
sh = AF.zhukhovsky_shape(b=sec["b"], e=sec["e"], beta=sec["beta"], exact=True)  # chord, camber and thickness measured on the contour
rule = AF.zhukhovsky_shape(b=sec["b"], e=sec["e"], beta=sec["beta"], exact=False)   # the same from the rules of (14.8)
print(f"chord            {sh['chord']:.4f} m   (rule 4b = {rule['chord']:.4f} m)")          # 1.0055 against 1.0
print(f"camber/chord     {sh['camber_ratio']:.4f}     (rule β/2 = {rule['camber_ratio']:.4f})")   # 0.0346 against 0.0349
print(f"thickness/chord  {sh['thickness_ratio']:.4f}     (rule 1.299 e = {rule['thickness_ratio']:.4f})")   # 0.0964 against 0.1039
print(f"thickest at      {100*sh['x_tmax_over_c']:.1f} % of the chord")                     # 25 %: the quarter chord
""", explain=r"""
1. `AF.zhukhovsky_airfoil` returns the circle's points `zeta`, their images `z` (and `x`, `y`) for our $b$ = 0.25 m, $e$ = 0.08,
   β = 4°.
2. `AF.zhukhovsky_shape(exact=True)` measures the contour; `exact=False` returns the three rules of $c\simeq4b$, camber
   $\simeq\tfrac12\beta c$, $t_{\max}/c\simeq1.3\,e$ (14.8).
3. Chord and camber agree with the rules to better than 1 %; the thickness rule is 8 % high at $e$ = 0.08; the thickest
   point is at the quarter chord, as derived.""")
scratch(r"""
b, e, beta = sec["b"], sec["e"], sec["beta"]                 # map constant [m], thickness parameter, camber angle [rad]
t = -beta + 2*np.pi*np.arange(400)/400                       # 400 angles round the circle from its centre, starting at the tail point (t = −β)
zeta = mu + a*np.exp(1j*t)                                   # the circle: centre μ, radius a (both from the imports cell)
z = zeta + b**2/zeta                                         # Eq. (14.3): the Zhukhovsky map, applied to all 400 points at once
chord_mine = z.real.max() - z.real.min()                     # chord = extent in x [m]
assert np.allclose(z, S["z"])                                # same points → the library does exactly this
assert np.isclose(chord_mine, sh["chord"], rtol=1e-4)        # and the same chord
print(f"chord = {chord_mine:.4f} m; the tail point maps to z = {z[0]:.4f} (= 2b)")   # 1.0055 m; 0.5
""", r"""
1. One line builds the circle from its centre and radius as complex numbers; one line maps it with
   $z=\zeta+\dfrac{b^2}{\zeta}$ (14.3).
2. The asserts confirm that `AF.zhukhovsky_airfoil` and `AF.zhukhovsky_shape` do the same.""")
recap("R07", "Corner flows", r"""
A flow $w=Az^n$ turns a corner; here as the flow *leaving* a trailing edge of finite included angle $2\pi(1-1/n)$:
$\phi=\dfrac{Ud}{n}\Big(\dfrac rd\Big)^{n}\cos n\theta$, with speed $U(r/d)^{n-1}$ → 0 at the edge for $n>1$.""", where="Ch. 6 §6.4", core_id="C04")
nb.code(r"""
wedge = ch14.trailing_edge_wedge(n=1.05)                     # the corner flow leaving an edge with exponent n = 1.05
print(f"included angle of the edge = {wedge['included_angle']:.3f} rad = {np.degrees(wedge['included_angle']):.1f}°")   # 0.299 rad = 17.1°
print(f"speed at r/d = 1e-2, 1e-4, 1e-8: {[round(float(wedge['speed'](r_)), 3) for r_ in (1e-2, 1e-4, 1e-8)]} U")        # falls only slowly: r^0.05
""")
D("D05", "The trailing edge is a cusp; a finite edge angle would need a stagnation point", ref="14.3")
nb.figure(r"""
S_fine = AF.zhukhovsky_airfoil(b=sec["b"], e=sec["e"], beta=sec["beta"], n=4000)   # our section with many points, to zoom on the tail
print(f"included angle at the tail measured on 400 and 4000 points: {np.degrees(AF.trailing_edge_angle(x=S['x'], y=S['y'])):.2f}°, {np.degrees(AF.trailing_edge_angle(x=S_fine['x'], y=S_fine['y'])):.2f}°")   # tends to 0: a cusp
fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.2))          # zoom on the tail; an edge of finite angle; a cusp
axes[0].plot(S_fine["x"], S_fine["y"], color="k")            # the outline
axes[0].set(xlim=(0.48, 0.5005), ylim=(-0.004, 0.004), xlabel="x [m]", ylabel="y [m]", title="Last 2 % of our section:\nthe two surfaces become tangent")   # zoom window (the tail is at x = 2b = 0.5 m)
rr, tt = np.meshgrid(np.linspace(0.01, 1.0, 80), np.linspace(-np.pi, np.pi, 361))   # polar grid round the edge: radius r/d and angle θ
for ax, n_ in zip(axes[1:], (1.05, 1.0)):                    # exponent 1.05 (a 17° edge) and 1 (a cusp)
    psi = np.where(np.abs(tt) < np.pi/n_, rr**n_*np.sin(n_*tt), np.nan)   # stream function r^n sin(nθ) in the fluid; the body fills |θ| > π/n
    ax.contour(rr*np.cos(tt), rr*np.sin(tt), psi, levels=np.linspace(-0.9, 0.9, 13), colors=COLORS["blue"], linewidths=0.8)   # streamlines leaving the edge to the right
    ax.fill([0, -np.cos(np.pi*(1 - 1/n_)), -np.cos(np.pi*(1 - 1/n_))], [0, np.sin(np.pi*(1 - 1/n_)), -np.sin(np.pi*(1 - 1/n_))], color="0.3")   # the solid edge (a wedge; a line for the cusp)
    ax.set(xlabel="x/d", ylabel="y/d", title=("edge of 17°: speed ∝ r^0.05 → 0 at the edge" if n_ > 1 else "cusp: speed = U at the edge"))   # labels
    ax.set_aspect("equal")                                   # true proportions
plt.show()                                                   # display
""",
    see="Left: the last 2 % of our section, strongly magnified — upper and lower surface meet tangent to each other. Middle: "
        "streamlines leaving an edge with a 17° included angle. Right: streamlines leaving a cusp.",
    read="At an edge of finite angle the streamlines spread apart as they approach the edge: the fluid slows to rest there "
         "(a stagnation point). At a cusp they stay evenly spaced: the fluid leaves at a finite speed. The printed angles "
         "show the measured tail angle shrinking toward zero as the outline is resolved more finely.",
    change="…the edge were blunter (larger $n$)? The speed would vanish faster, like $r^{n-1}$.")
nb.animation(r"""
n_frames = 60 if not FAST else 30                            # number of frames
half = n_frames//2                                           # first half: raise the centre; second half: move it back
tt = np.linspace(0, 2*np.pi, 241)                            # angle round the circle
fig, (axZ, axz) = plt.subplots(1, 2, figsize=(9.0, 3.3))     # left: the ζ-plane (circle); right: the z-plane (section)
axZ.plot(sec["b"]*np.cos(tt), sec["b"]*np.sin(tt), color="0.7", ls=":")   # ghost circle of radius b
axZ.plot([sec["b"]], [0], "o", color=COLORS["orange"])       # the point ζ = b, which becomes the tail
(circ,) = axZ.plot([], [], color=COLORS["accent"])           # the moving circle
(cen,) = axZ.plot([], [], "x", color=COLORS["accent"])       # its centre μ
(outline,) = axz.plot([], [], color="k")                     # the section
axZ.set(xlim=(-0.4, 0.4), ylim=(-0.36, 0.4), xlabel="ξ [m]", ylabel="η [m]")   # fixed axes
axZ.set_aspect("equal")                                      # a circle looks like a circle
axz.set(xlim=(-0.6, 0.6), ylim=(-0.2, 0.25), xlabel="x [m]", ylabel="y [m]")   # fixed axes
axz.set_aspect("equal")                                      # true proportions
fig.canvas.draw(); fig.set_layout_engine("none")             # lay the figure out once, then freeze it: frames render much faster
def update(i):                                               # draw frame i
    beta_ = np.deg2rad(8.0)*min(i/(half - 1), 1.0)           # β grows from 0 to 8° during the first half
    e_ = 0.12*max((i - half)/(n_frames - half - 1), 0.0)     # then e grows from 0 to 0.12
    Sx = AF.zhukhovsky_airfoil(b=sec["b"], e=e_, beta=beta_, n=240)   # circle and section for these two parameters
    circ.set_data(Sx["zeta"].real, Sx["zeta"].imag)          # the circle
    cen.set_data([Sx["mu"].real], [Sx["mu"].imag])           # its centre
    outline.set_data(Sx["x"], Sx["y"])                       # its image
    axZ.set_title(f"circle: β = {np.degrees(beta_):.1f}° ({beta_:.3f} rad), e = {e_:.3f}")   # the two knobs, angle in both units
    axz.set_title("flat plate → arc → cambered airfoil")     # what the right panel shows
show_animation(animate(update, frames=n_frames, fig=fig, interval=80), player="video")   # a smooth video
""", explain=r"""
1. The first half of the frames raises the circle's centre (β: 0 → 8°) at zero thickness; the second half moves it back
   ($e$: 0 → 0.12) at fixed β.
2. Each frame calls `AF.zhukhovsky_airfoil` and updates three curves: the circle, its centre, and the mapped section.""")
nb.figure_notes(
    see="Left: the circle and its centre (×) moving while it stays pinned at the amber point $\\zeta=b$. Right: its image.",
    read="While the centre rises, the plate bends into an arc but stays infinitely thin. When the centre then moves back, the "
         "arc fattens from the nose; the tail stays sharp because the circle never lets go of $\\zeta=b$.",
    change="…the centre moved back first and then up? The same final section: the two knobs are independent.")
nb.plotly(r"""
def section_e(e_):                                           # outline of a symmetric section for one thickness parameter (b = 0.25 m)
    Sx = AF.zhukhovsky_airfoil(b=0.25, e=e_, beta=0.0, n=200)   # circle shifted back by b·e
    return {"symmetric section (β = 0)": (Sx["x"], Sx["y"])}   # one closed outline
figF2a = slider_figure(section_e, "thickness parameter e", np.linspace(0.0, 0.25, 26), unit="", xlabel="x [m]", ylabel="y [m]",
                       xrange=[-0.62, 0.62], yrange=[-0.3, 0.3], title="Move the circle back: the nose fattens, the tail stays sharp")   # 26 slider positions
figF2a.show()                                                # display

def section_beta(beta_deg):                                  # outline of a cambered section for one camber angle, and its skeleton
    B_ = np.deg2rad(beta_deg)                                # the slider is in degrees; the functions take radians
    Sx = AF.zhukhovsky_airfoil(b=0.25, e=0.08, beta=B_, n=200)   # the section with e = 0.08
    Sa = AF.zhukhovsky_airfoil(b=0.25, e=0.0, beta=B_, n=200)    # the same circle angle with no thickness: a circular arc
    return {"section (e = 0.08)": (Sx["x"], Sx["y"]), "its skeleton: the circular arc of the same β": (Sa["x"], Sa["y"])}   # two outlines
figF2b = slider_figure(section_beta, "camber angle β", np.linspace(0.0, 12.0, 25), unit="°", xlabel="x [m]", ylabel="y [m]",
                       xrange=[-0.62, 0.62], yrange=[-0.3, 0.3], title="Raise the circle: the section bends round a circular arc")   # 25 slider positions
figF2b.show()                                                # display
""", explain=r"""
1. `section_e` returns the outline for one value of $e$; the slider runs from a flat plate ($e$ = 0) to a very fat fin ($e$ = 0.25).
2. `section_beta` returns the outline for one camber angle (shown in degrees, converted to radians inside) together with
   the circular arc of the same β, which is close to the section's camber line.
3. The measured numbers behind both sliders are the two tables printed under N26 above.""")
nb.figure_notes(
    see="Two slider figures: thickness growing at zero camber, and camber growing at fixed thickness.",
    read="In the first the thickest point stays near the quarter chord while the section grows; in the second the whole "
         "section bends, following the circular arc that the same β would give without thickness.",
    change="…$e$ were raised to 0.25? The rule $t_{\\max}/c\\simeq1.3\\,e$ would be a quarter too high (table above) — the section is "
           "fat, but less fat than the rule says.")
explainer("zhukhovsky_airfoil_builder", "Where do camber and thickness come from?",
          "dragging the circle's centre and watching the section reshape continuously shows that β and e are two independent "
          "knobs, and how far the rules drift as the section gets fat.",
          ["Drag the centre straight up: the plate bends but stays thin.",
           "Drag it left: the nose fattens, the tail stays sharp.",
           "Press 'too fat' and read the gap between the rule and the measured thickness.",
           "Click a point of the circle and follow its arithmetic to the section."])
whatif(r"""…a stream now blows past this shape? We already know the flow past the circle: C05 carries it over.""")

# =====================================================================================================================
# A.5 §14.5 Lift of a Zhukhovsky Airfoil — C05, C06, C07
# =====================================================================================================================
nb.section("14.5", "Lift of a Zhukhovsky Airfoil", intro=r"""
**What is this section about?** Carry the circle's flow to the airfoil, find the one circulation that keeps the trailing-edge
speed finite, and read off the lift — then look at the pressure that produces it.""")
core("C05", r"The flow round the airfoil: $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9)",
     "If we know the flow round the circle, what is the flow round the airfoil?")
problem(r"""
The map is a change of address: every point ζ outside the circle has a partner $z$ outside the airfoil. The complex potential
— the pair (potential, stream function) — keeps its *value* at partner points, so streamlines map to streamlines. The *speed*
does not keep its value: where the map stretches the plane the fluid slows down, where it squeezes the plane the fluid speeds
up. At the tail the map squeezes infinitely hard.""")
idea(r"""
   ζ-plane: circle + flow   ──(same w at partner points)──►   z-plane: airfoil + flow

   speed_z  =  speed_ζ × |dζ/dz|        ← the factor is 1 far away and ∞ at ζ = b (the tail)
""", r"""In this block $w$ is the **complex potential** $w=\phi+i\psi$ (not the downwash of §14.6).""")
remind("C05")
D("D06", "The flow past an off-centre circle at incidence", ref="6.68")
nb.code(r"""
tcirc = np.linspace(0, 2*np.pi, 721)                         # angle round the circle
w_c, dw_c = CONF.circle_flow_zeta(U=30.0, a=a, Gamma_cw=10.0, alpha=sec["alpha"], center=mu)   # potential and complex velocity: 30 m/s at α = 5°, clockwise Γ = 10 m²/s, our circle
print("complex velocity far away:", np.round(dw_c(1e6 + 0j), 3), " expected 30·e^(−iα) =", np.round(30*np.exp(-1j*sec["alpha"]), 3))   # 29.886 − 2.615i: speed 30 at +5°
print(f"spread of the stream function round the circle: {np.ptp(np.imag(w_c(mu + a*np.exp(1j*tcirc)))):.1e} m²/s")   # ≈ 1e-13: the circle is a streamline
""", explain=r"""
**What does this show?** The two checks of the derivation's steps 5 and 6, in numbers: far away the complex velocity is
$Ue^{-i\alpha}$ = 29.886 − 2.615i (remember $u-iv$: so $v$ = +2.615 m/s, the stream points *up*), and the imaginary part of $w$ —
the stream function — is the same all round the circle.""")
recap("R08", "The chain rule for the complex velocity", r"""
$\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\dfrac{d\zeta}{dz}$: the velocity in the $z$-plane is the velocity in the ζ-plane times the map's local rescaling.
The chain rule holds for analytic functions exactly as for real ones.

> ⚠️ **Common confusion (trap T2):** here $dw/dz=u-iv$ — note the minus sign. Chapter 13's Ekman velocity was $V=u+iv$.""",
      where="Ch. 6 §6.6", core_id="C05")
D("D07", "The velocity in the airfoil plane", ref="14.9")
D("D08", "The circulation is the same in both planes", ref="6.62")
nb.code(r"""
F10 = AF.zhukhovsky_flow(**SEC, Gamma_cw=10.0)               # our section's flow with a clockwise circulation of 10 m²/s
tq = np.linspace(0, 2*np.pi, 2001)[:-1]                      # 2000 equally spaced angles round a loop (periodic trapezoid rule)
for R_ in (1.5, 3.0, 10.0):                                  # three loops round the circle, of radius 1.5a, 3a and 10a
    zeta_l = mu + R_*a*np.exp(1j*tq)                         # the loop in the ζ-plane
    dzeta = 1j*R_*a*np.exp(1j*tq)*(tq[1] - tq[0])            # its steps dζ
    circ_zeta = np.sum(F10.dwdzeta(zeta=zeta_l)*dzeta).real  # ∮ (dw/dζ) dζ: counterclockwise circulation in the circle plane
    dz = CONF.joukowski_derivative(zeta=zeta_l, b=sec["b"])*dzeta   # the partner loop's steps dz = (dz/dζ) dζ
    circ_z = np.sum(F10.dwdz_from_zeta(zeta=zeta_l)*dz).real # ∮ (dw/dz) dz: counterclockwise circulation in the airfoil plane
    ratio = abs(F10.dwdz_from_zeta(zeta=zeta_l[0]))/abs(F10.dwdzeta(zeta=zeta_l[0]))   # speed in the z-plane over speed in the ζ-plane at one partner point
    print(f"loop radius {R_:4.1f} a: circulation {circ_zeta:+.4f} (ζ-plane), {circ_z:+.4f} (z-plane) m²/s; speed ratio at one point = {ratio:.3f}")   # −10.000 twice; ratio → 1 only far away
""", explain=r"""
1. Each loop is sampled at 2000 points; the sum of (complex velocity × step) is the periodic trapezoid rule for the loop
   integral, and its real part is the counterclockwise circulation.
2. In both planes and for all three loops the result is −10.000 m²/s: a clockwise circulation of 10, as put in.
3. The last column is the trap (T6): the *speed* at partner points is **not** the same — only its loop integral is. The
   ratio approaches 1 only far from the body.""")
nb.worked_example("three points by hand", r"""
Flat plate ($a=b$ = 1), $U$ = 1, using $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9).

1. α = 0, Γ = 0, ζ = 2: circle velocity $1-a^2/\zeta^2$ = 0.75; map factor $\zeta^2/(\zeta^2-b^2)$ = 4/3; product 1.00 — a plate edge-on
   does not disturb the stream ✓.
2. The same at ζ = 2i: (1 + 0.25) × (−4)/(−5) = 1.25 × 0.8 = 1.00 ✓.
3. Now α = 30° and ζ = 1.1, just off the trailing edge: the circle velocity is $e^{-i\alpha}-e^{i\alpha}/1.21$ = 0.150 − 0.913i
   (size 0.926); the factor is 1.21/0.21 = 5.76; the speed is 5.33 $U$ — and it grows without limit as ζ → 1.
4. Add the Kutta circulation Γ = 2π (C02's hand example): the extra term $i\Gamma/(2\pi\zeta)$ = +0.909i leaves 0.150 − 0.004i, and
   the speed is 0.150 × 5.76 = 0.87 $U$: finite.""")
nb.code(r"""
F = AF.zhukhovsky_flow(**SEC)                                # our section's flow; with no Gamma_cw given it takes the Kutta value
print("velocity above mid-chord (0, 0.3) m:", np.round(F.velocity(x=0.0, y=0.3), 2), "m/s")   # faster than the stream: about (36.3, 0.1)
print("velocity far upstream (−5, 0) m:    ", np.round(F.velocity(x=-5.0, y=0.0), 2), "m/s")   # (29.9, 3.1): the free stream (29.9, 2.6) plus a little upwash ahead of the lifting section
print(f"speed 1 cm behind the tail:          {F.speed(z=2*sec['b'] + 0.01):.2f} m/s")          # finite, a little below U
print("is (0, 0) inside the section?       ", F.inside(x=0.0, y=0.0))                        # True: no fluid there (the velocity is NaN)
""", explain=r"""
1. `AF.zhukhovsky_flow` returns an object that knows the circle, the map and the circulation; called without `Gamma_cw` it
   uses the Kutta value (C06).
2. `F.velocity(x, y)` returns the two velocity components $(u, v)$ at a point of the airfoil plane.
3. Off the surface the function needs ζ for a given $z$: the inverse map has two roots and the one outside the *shifted*
   circle is taken (trap T8) — `np.sqrt` alone would not always choose it.
4. Points inside the section are flagged by `F.inside`.""")
scratch(r"""
U, b, al = sec["U"], sec["b"], sec["alpha"]                  # speed [m/s], map constant [m], angle of attack [rad]
zeta0 = mu + 1.3*a*np.exp(0.7j)                              # a point outside the circle: 1.3 radii from the centre, 0.7 rad round
dw = U*(np.exp(-1j*al) - a**2*np.exp(1j*al)/(zeta0 - mu)**2) + 1j*F.Gamma_cw/(2*np.pi*(zeta0 - mu))   # circle velocity: derivative of the potential of D06
q = dw*zeta0**2/(zeta0**2 - b**2)                            # Eq. (14.9): rescale by the map
assert np.isclose(q, F.dwdz(z=zeta0 + b**2/zeta0))           # same complex velocity at the partner point z → the library does exactly this
print(f"u − iv = {q:.3f} m/s at z = {zeta0 + b**2/zeta0:.3f} m")   # the velocity there
""", r"""
1. `dw` is the derivative of the potential of derivation D06: stream, its image in the circle, and the swirl.
2. `q` applies $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9).
3. The assert compares with `F.dwdz` at the partner point $z=\zeta+b^2/\zeta$.""")
nb.figure(r"""
GK = F.Gamma_kutta                                           # the Kutta circulation of our section [m²/s], clockwise
w_k, dw_k = CONF.circle_flow_zeta(U=sec["U"], a=a, Gamma_cw=GK, alpha=sec["alpha"], center=mu)   # the circle's flow with that circulation
st = AF.stagnation_points(**SEC)                             # front and rear stagnation points in both planes
nx, ny = (160, 100) if not FAST else (100, 64)               # grid size
XZ, YZ = np.meshgrid(np.linspace(-0.9, 0.9, nx), np.linspace(-0.6, 0.6, ny))   # the same window in both planes [m]
psi_zeta = np.imag(w_k(XZ + 1j*YZ))                          # stream function in the circle plane
psi_zeta[np.abs(XZ + 1j*YZ - mu) < a] = np.nan               # no fluid inside the circle
psi_z = F.stream(x=XZ, y=YZ)                                 # stream function in the airfoil plane (NaN inside the section)
levels = np.linspace(-14.0, 14.0, 15)                        # the same 15 streamlines in both panels [m²/s]; 0 is the body
fig, (axa, axb) = plt.subplots(1, 2, figsize=(10.0, 3.6))    # left: ζ-plane; right: z-plane
axa.contour(XZ, YZ, psi_zeta, levels=levels, colors=COLORS["blue"], linewidths=0.8)   # streamlines round the circle
axa.plot((mu + a*np.exp(1j*tcirc)).real, (mu + a*np.exp(1j*tcirc)).imag, color="k")   # the circle
axa.plot([st["zeta_front"].real, st["zeta_rear"].real], [st["zeta_front"].imag, st["zeta_rear"].imag], "o", color=COLORS["orange"])   # its stagnation points
axa.set(xlabel="ξ [m]", ylabel="η [m]", title="ζ-plane: the circle's flow")   # labels
axb.contour(XZ, YZ, psi_z, levels=levels, colors=COLORS["blue"], linewidths=0.8)      # streamlines round the airfoil
axb.fill(S["x"], S["y"], color="0.3")                        # the section
axb.plot([st["z_front"].real, st["z_rear"].real], [st["z_front"].imag, st["z_rear"].imag], "o", color=COLORS["orange"])   # its stagnation points
axb.set(xlabel="x [m]", ylabel="y [m]", title="z-plane: the same streamlines, mapped")   # labels
for ax in (axa, axb): ax.set_aspect("equal")                 # true proportions
plt.show()                                                   # display
""",
    see="The same fifteen streamlines in the two planes: round the circle (left) and round our section (right), with the "
        "front and rear stagnation points as amber dots.",
    read="Partner points carry the same stream-function value, so each streamline on the left is the pre-image of one on the "
         "right. The rear dot sits at the point of the circle that maps to the tail. Far from the body the two pictures are "
         "the same; near it the right one is squeezed flat.",
    change="…Γ were not the Kutta value? The rear dot would sit beside the tail and a streamline would wrap round it: C06.")
nb.plotly(r"""
def edge_speeds(f):                                          # speed along the last 5 % of the upper surface for one circulation Γ = f·Γ_K
    Sf = AF.zhukhovsky_surface(**SEC, n=2000, Gamma_cw=f*GK) # surface speeds with this circulation
    up = Sf["upper"]                                         # the upper surface, x/c increasing
    keep = (up["xc"] > 0.95) & np.isfinite(up["speed"])      # last 5 % of the chord, without the (possibly infinite) edge value
    return {"upper surface, last 5 % of the chord": (up["xc"][keep], np.log10(up["speed"][keep]/sec["U"]))}   # log10 of q/U: +1 means ten times U
figF3 = slider_figure(edge_speeds, "Γ/Γ_K (Γ clockwise)", np.linspace(0.5, 1.5, 21), unit="", xlabel="x/c [–]",
                      ylabel="log10(q/U) [–]", xrange=[0.95, 1.0005], yrange=[-0.3, 1.3], active=10,
                      title="Only at Γ/Γ_K = 1.00 does the speed stay finite at the edge")   # 21 slider positions, starting at 1.00
figF3.show()                                                 # display
""", explain=r"""
1. `AF.zhukhovsky_surface` returns the speed on the contour; `Gamma_cw=f*GK` sets a circulation off the Kutta value.
2. The trace is the base-10 logarithm of speed over $U$ on the last 5 % of the upper surface: 0 means $q=U$, 1 means $10\,U$.""")
nb.figure_notes(
    see="The speed on the upper surface close to the trailing edge, on a logarithmic scale, for circulations from half to "
        "one and a half times the Kutta value.",
    read="At 1.00 the curve runs flat into the edge at a value a little below zero (a speed just under $U$). At every other "
         "slider position it turns upward at the right-hand end — the closer to the edge, the faster — because the map's "
         "factor $\\zeta^2/(\\zeta^2-b^2)$ is infinite there and the circle's flow no longer vanishes to cancel it.",
    change="…we *demand* a finite speed at the tail? That fixes Γ: C06.")
whatif(r"""…we *demand* a finite speed at the tail? That single demand fixes Γ, and with it the lift: C06.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C06", r"The lift of a Zhukhovsky airfoil: $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) ⇒ $C_L\simeq2\pi(\alpha+\beta)$ (14.12)",
     "How much lift does a section give — and what do angle of attack and camber each contribute?")
problem(r"""
A pilot pulls the nose up two degrees and the lift rises by a fifth. A flap comes down and the lift rises without the nose
moving at all. Both facts are in one line of geometry: the lift is proportional to the angle between the stream and one
particular direction of the circle — the line from its centre to the point that becomes the tail. Tilting the wing changes
that angle by α; cambering the wing changes it by β.""")
idea(r"""
                  circle in the ζ-plane, centre Q
                        .  -  ~  -  .
                     .                 .
                    (        Q ─ ─ ─ ─ ─)─ ─ ─  horizontal through Q
                     .         ╲  β    .
                        ` -  . _╲_ . '
     stream ↗ at α                B  (ζ = b: becomes the tail)

     angle between the stream direction and the line QB  =  α + β   ← the "absolute" angle of attack
""")
recap("R09", "Surface speed on a circle with circulation", r"""
$u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (14.10), the same formula as $u_\theta=-2U\sin\theta-\dfrac{\Gamma}{2\pi a}$ (6.37): the first term is the stream
sliding round the circle, the second the clockwise swirl. **New here:** θ is measured from the diameter through the centre Q
that is aligned with the stream.

> ⚠️ **Common confusion (trap T4):** three different θ's in this chapter — from the origin in
> $x=(R+b^2/R)\cos\theta$ (14.4); from the stream-aligned diameter through Q here; and a course angle in §14.9.""",
      where="Ch. 6 §6.3", core_id="C06")
nb.code(r"""
fig = draw_circle_angles()                                   # the off-centre circle: centre Q, the tail point B, the angles α and β, the stream-aligned diameter
plt.show()                                                   # display the sketch the derivation below refers to
GK = AF.kutta_circulation(U=sec["U"], a=a, alpha=sec["alpha"], beta=sec["beta"])   # Eq. (14.11): the Kutta circulation of our section [m²/s], clockwise
print(ch06.cylinder_surface_speed(theta=-(sec["alpha"] + sec["beta"]), U=sec["U"], a=a, Gamma_cw=GK))   # surface speed at B, θ = −(α + β) = −9°: zero to round-off
""", tags=["figure"])
P("P338", "radians and degrees for angles and for slopes", r"""
Formulas want radians: 1° = π/180 = 0.01745 rad. A slope "per radian" is 57.3 times a slope "per degree": 2π per radian is
0.110 per degree.""", r"""
print(np.deg2rad(9.0))                    # 9 degrees = 0.1571 rad
print(2*np.pi, 2*np.pi*np.pi/180)         # lift slope: 6.283 per radian = 0.1097 per degree
""")
D("D09", "The Kutta circulation and the lift coefficient", ref="14.12")
slip(8, r"""“the chord length $c$ (≃ 4a)” beside $c\cong4b$ in (14.8)""",
     r"""both are first-order statements; exactly $C_L=8\pi(a/c)\sin(\alpha+\beta)$ (step 8), and the slope exceeds 2π at finite thickness.""")
note("N34 [B]", r"""
**The flat plate** ($e=\beta=0$, $a=b=c/4$): $\Gamma=4\pi Ub\sin\alpha$ and $C_L=2\pi\sin\alpha$ **exactly** — no thin-section or small-angle
step; this is the flow of C02.""")
nb.code(r"""
print(AF.flat_plate_flow(U=30.0, b=0.25, alpha=sec["alpha"]).Gamma_kutta)   # 8.214 m²/s (clockwise): the number the root-finder of C02 found
""")
P("P343", "l'Hôpital's rule for a 0/0 limit", r"""
If top and bottom of a fraction both go to zero at a point, the limit of the fraction is the limit of (slope of the top)/(slope
of the bottom).""", r"""
x = sp.symbols("x")                                # a symbol
print(sp.limit(sp.sin(3*x)/sp.sin(x), x, 0))       # 3: slopes 3cos(0)/cos(0)
print(sp.limit((x**2 - 1)/(x - 1), x, 1))          # 2: slopes 2x/1 at x = 1
""")
note("N29 [B]", r"""
At the trailing edge $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) is 0/0 when the Kutta condition holds. **Not in the book — ours:** the
limit is finite, $\Big\lvert\dfrac{dw}{dz}\Big\rvert_{B'}=U\,\dfrac ba\cos(\alpha+\beta)$ — stated, not derived; verified two ways in the next cell.""")
nb.code(r"""
zq = sp.symbols("zeta")                                      # a point of the ζ-plane, as a symbol
al30 = sp.pi/6                                               # flat plate (a = b = 1, U = 1) at α = 30°
dw_plate = sp.exp(-sp.I*al30) - sp.exp(sp.I*al30)/zq**2 + sp.I*2*sp.sin(al30)/zq   # circle velocity with the Kutta circulation Γ = 4π sin α, i.e. iΓ/(2πζ) = 2i sin α/ζ
top, bottom = dw_plate*zq**2, zq**2 - 1                      # Eq. (14.9) as a fraction: (dw/dζ)·ζ² over ζ² − b²
assert abs(complex(sp.N(top.subs(zq, 1)))) < 1e-12 and bottom.subs(zq, 1) == 0   # both vanish at the tail point ζ = 1: a 0/0 form
lim = sp.diff(top, zq).subs(zq, 1)/sp.diff(bottom, zq).subs(zq, 1)   # l'Hôpital: slope of the top over slope of the bottom
print(f"flat plate at 30°: |dw/dz| at the edge = {abs(complex(sp.N(lim))):.6f} U; cos 30° = {float(sp.cos(al30)):.6f}")   # 0.866025 both
print(f"our section: U (b/a) cos(α + β) = {AF.trailing_edge_speed(U=sec['U'], a=a, b=sec['b'], alpha=sec['alpha'], beta=sec['beta']):.2f} m/s = {AF.trailing_edge_speed(U=sec['U'], a=a, b=sec['b'], alpha=sec['alpha'], beta=sec['beta'])/sec['U']:.3f} U")   # 27.37 m/s = 0.912 U
tab = pd.DataFrame({f"Γ/Γ_K = {f}": [ch14.edge_speed(**SEC, d=d_, Gamma_cw=f*GK)/sec["U"] for d_ in (1e-2, 1e-3, 1e-4, 1e-5)] for f in (1.0, 0.9, 0.5)}, index=["d = 1e-2 m", "1e-3 m", "1e-4 m", "1e-5 m"])   # speed/U at four distances d from the edge, three circulations
print(tab.round(3))                                          # the Kutta column settles; the others grow like d^(−1/2)
""", explain=r"""
1. The sympy part builds the circle velocity of a flat plate with the Kutta circulation, checks that top and bottom of
   $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9) both vanish at the edge, and applies l'Hôpital's rule (the primer just above): the
   result has the size cos 30°, as the formula of N29 says.
2. `AF.trailing_edge_speed` evaluates $U(b/a)\cos(\alpha+\beta)$ for our section: 27.37 m/s = 0.912 $U$.
3. The table (trap T8): with the Kutta circulation the speed tends to 0.912 $U$ as the distance $d$ shrinks; with 0.9 or 0.5 of
   it the speed grows like $d^{-1/2}$.""")
note("N31 [B]", r"""
α + β is the **absolute angle of attack**; the lift vanishes at $\alpha_{L=0}=-\beta$, and for these sections
$\beta=2\,\dfrac{\text{camber}}{c}$. More camber slides the whole lift line to the left.

> ⚠️ **Common confusion (trap T3):** β is the circle's camber angle here; in §14.6 the same letter is *defined* as minus the
> zero-lift angle of whatever section the wing has.""")
nb.code(r"""
aL0 = AF.zero_lift_angle(camber_ratio=0.035)                 # zero-lift angle of a section with 3.5 % camber [rad]
print(f"zero-lift angle = {aL0:.3f} rad = {np.degrees(aL0):.1f}°")   # −0.070 rad = −4.0° = −β of our section
""")
nb.worked_example("lift by hand", r"""
$U$ = 10 m/s, $a$ = 0.25 m (chord ≈ 1 m), ρ = 1.2 kg/m³, and α + β = 30° so that the sine is easy.

1. $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11) = 4π × 10 × 0.25 × 0.5 = 15.7 m²/s (clockwise).
2. $L=\rho U\Gamma$ (6.62) = 1.2 × 10 × 15.7 = 188 N per metre.
3. $C_L=L/(\tfrac12\rho U^2c)$ = 188/(0.5 × 1.2 × 100 × 1) = 3.14 (= 2π sin 30°).
4. The small-angle form $C_L\simeq2\pi(\alpha+\beta)$ (14.12) gives 2π × 0.524 = 3.29 — 5 % high at 30°, which is far beyond where a
   real section has stalled.
5. A realistic case: α = 5°, β = 4°: α + β = 0.157 rad, $C_L$ ≈ 0.99, and sine and angle differ by 0.4 %.""")
nb.code(r"""
GK = AF.kutta_circulation(U=sec["U"], a=a, alpha=sec["alpha"], beta=sec["beta"])   # Eq. (14.11): Γ = 4πUa sin(α + β) [m²/s], clockwise
L_sec = ch06.lift_per_span(rho=sec["rho"], U=sec["U"], Gamma_cw=GK)                 # Eq. (6.62): L = ρUΓ [N per metre of span]
CL_rule = AF.zhukhovsky_lift_coefficient(alpha=sec["alpha"], beta=sec["beta"])      # Eq. (14.12): 2π(α + β)
CL_exact = AF.zhukhovsky_lift_coefficient(alpha=sec["alpha"], beta=sec["beta"], e=sec["e"], exact=True)   # 8π(a/c) sin(α + β) with the measured chord
G_root = AF.kutta_circulation_root(**SEC)                    # the root-finder of C02 (it knows no formula), for our section
print(f"Γ_K = {GK:.3f} m²/s (clockwise); root-finder: {G_root:.3f} m²/s")           # 15.962 twice
print(f"L = {L_sec:.1f} N/m; C_L = {CL_rule:.3f} (rule), {CL_exact:.3f} (exact), ratio {CL_exact/CL_rule:.3f}")   # 586.6; 0.987; 1.058; 1.072
""", explain=r"""
1. `AF.kutta_circulation` is $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11): 15.962 m²/s clockwise for our section (α = 5°, β = 4°, so
   α + β = 9° = 0.157 rad).
2. `ch06.lift_per_span` is $L=\rho U\Gamma$ (6.62): 586.6 N per metre of span.
3. `AF.zhukhovsky_lift_coefficient` gives the rule $C_L\simeq2\pi(\alpha+\beta)$ (14.12), 0.987, and with `exact=True` the exact
   value 1.058.
4. The exact value is 7.2 % above the rule: thickness makes $a$ larger than $c/4$ (slip #8). The root-finder of C02 agrees
   with the formula to round-off.""")
scratch(r"""
U, alpha, beta, e = sec["U"], sec["alpha"], sec["beta"], sec["e"]   # speed [m/s], angle of attack [rad], camber angle [rad], thickness parameter
G = 4*np.pi*U*a*np.sin(alpha + beta)                         # Eq. (14.11): the Kutta circulation [m²/s], clockwise
c = sh["chord"]                                              # the measured chord [m]
CL = 2*G/(U*c)                                               # lift coefficient: ρUΓ/(½ρU²c) = 2Γ/(Uc)
assert np.isclose(CL, AF.zhukhovsky_lift_coefficient(alpha=alpha, beta=beta, e=e, exact=True), rtol=1e-6)   # same number → the library does exactly this
print(f"C_L = {CL:.4f}")                                     # 1.0583
""", r"""
Three lines: the circulation from $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11), the measured chord, and $C_L=2\Gamma/(Uc)$. The assert
compares with the tested function.""")
nb.plotly(r"""
al_deg = np.linspace(-10.0, 15.0, 101)                       # angle of attack on the horizontal axis [degrees]
def lift_lines(beta_deg):                                    # lift coefficient against α for one camber angle
    B_ = np.deg2rad(beta_deg)                                # degrees on the slider, radians in the formulas
    rule = AF.zhukhovsky_lift_coefficient(alpha=np.deg2rad(al_deg), beta=B_)                       # Eq. (14.12): 2π(α + β)
    k8 = AF.zhukhovsky_lift_coefficient(alpha=0.1, beta=B_, e=0.08, exact=True)/np.sin(0.1 + B_)   # the factor 8π a/c of the exact formula, read off one call (the chord does not depend on α)
    exact = k8*np.sin(np.deg2rad(al_deg) + B_)               # exact value 8π(a/c) sin(α + β) at e = 0.08
    return {"rule 2π(α + β): slope 0.110 per degree": (al_deg, rule), "exact, e = 0.08 (slope about 7 % higher)": (al_deg, exact),
            "zero-lift angle α = −β": (np.array([-beta_deg, -beta_deg]), np.array([-0.4, 0.4]))}   # two lines and a marker
figF4 = slider_figure(lift_lines, "camber angle β", np.linspace(0.0, 8.0, 17), unit="°", xlabel="angle of attack α [°]",
                      ylabel="lift coefficient C_L [–]", yrange=[-1.3, 2.9], title="Camber shifts the lift line; it never changes the slope")   # 17 slider positions
figF4.show()                                                 # display
""", explain=r"""
1. `lift_lines` returns the rule $C_L\simeq2\pi(\alpha+\beta)$ (14.12), the exact curve for $e$ = 0.08, and a short vertical marker at
   the zero-lift angle $\alpha=-\beta$.
2. Angles are in degrees on the axis and the slider, and converted with `np.deg2rad` before each call.""")
nb.figure_notes(
    see="Lift coefficient against angle of attack, with a slider for the camber angle.",
    read="The slider moves the *intercept* — the zero-lift angle slides to −β — and never the slope: 2π per radian, 0.110 per "
         "degree for the rule, about 7 % more for the section with thickness.",
    change="…the section were thinner? The dashed exact line would close on the rule; the flat plate has exactly 2π sin α.")
note("N32 [B]", r"""
**Theory against measurement.** Measured slopes are a little below 2π — the boundary layers' *displacement thickness* (how far
the layer pushes the outer flow away from the wall, Ch. 9) changes the effective shape — and the measured curve ends in
stall. The figure compares $C_L\simeq2\pi(\alpha+\beta)$ (14.12), the exact curve and a **qualitative model** of a measured curve whose
constants are ours.""")
nb.figure(r"""
al = np.deg2rad(np.linspace(-8.0, 20.0, 281))                # angles of attack [rad]
model = ch14.generic_lift_curve(alpha=al, K=5.7, beta=sec["beta"], alpha_stall=np.deg2rad(14.0), CL_post=0.9, width=np.deg2rad(2.0))   # a generic measured curve: slope 5.7/rad, stall at 14°, plateau 0.9 — constants ours
k8 = AF.zhukhovsky_lift_coefficient(alpha=0.1, beta=sec["beta"], e=sec["e"], exact=True)/np.sin(0.1 + sec["beta"])   # the factor 8π a/c of the exact formula (the chord does not depend on α)
exact = k8*np.sin(al + sec["beta"])                          # exact ideal-flow value 8π(a/c) sin(α + β)
fig, ax = plt.subplots(figsize=(6.6, 3.8))                   # one panel
ax.plot(np.degrees(al), AF.zhukhovsky_lift_coefficient(alpha=al, beta=sec["beta"]), color=COLORS["teal"], label=r"rule $2\pi(\alpha+\beta)$ (14.12)")   # the rule
ax.plot(np.degrees(al), exact, color=COLORS["teal"], ls="--", label=r"exact $8\pi(a/c)\sin(\alpha+\beta)$, e = 0.08")   # the exact ideal value
ax.plot(np.degrees(al), model, color=COLORS["rose"], label="qualitative model of a measured curve")   # the model with stall
ax.set(xlabel="angle of attack α [°]", ylabel="lift coefficient $C_L$ [–]", title="Theory (lines) and a Qualitative model of a measured curve — constants ours")   # labels
ax.legend(fontsize=8); plt.show()                            # legend and display
""",
    see="Two ideal-flow lines (rule and exact) and a rose curve that follows them at a slightly lower slope, then bends over "
        "and drops at 14°.",
    read="All three cross zero at α = −4° = −β. Below stall the rose curve is a straight line a little flatter than the "
         "theory. At the stall angle it leaves the line: the Kutta condition no longer describes the flow. The stall angle, "
         "the plateau and the slope 5.7 per radian are our illustrative constants, not data.",
    change="…the Reynolds number were lower? The stall would come earlier (N72 in C15); the straight part would hardly move.")
nb.live(r"""
def show_section(alpha_deg=5.0, beta_deg=4.0, e=0.08):       # redraw for the slider values (angles in degrees)
    A_, B_ = np.deg2rad(alpha_deg), np.deg2rad(beta_deg)     # the functions take radians
    Fl = AF.zhukhovsky_flow(U=30.0, b=0.25, e=e, beta=B_, alpha=A_)   # the section's Kutta flow
    fig, ax = plt.subplots(figsize=(6.4, 3.4))               # one panel
    streamlines(ax=ax, flow=Fl, xlim=(-0.9, 0.9), ylim=(-0.5, 0.5), n=(120, 80), levels=21)   # streamlines and the filled section
    ax.set(xlabel="x [m]", ylabel="y [m]", title=f"α = {alpha_deg:.0f}°, β = {beta_deg:.0f}°, e = {e:.2f}: C_L = {AF.zhukhovsky_lift_coefficient(alpha=A_, beta=B_, e=e, exact=True):.3f} (rule {AF.zhukhovsky_lift_coefficient(alpha=A_, beta=B_):.3f}), Γ_K = {Fl.Gamma_kutta:.2f} m²/s clockwise")   # three numbers
    plt.show()                                               # display
live(show_section, alpha_deg=(-5.0, 15.0, 1.0), beta_deg=(0.0, 10.0, 1.0), e=(0.0, 0.2, 0.02))   # three sliders (kernel only)
""", explain=r"""
**What does the code above do?** Three sliders — angle of attack, camber angle (both in degrees) and thickness parameter —
redraw the section with its streamlines and print the exact lift coefficient, the rule and the Kutta circulation. The slider
figure above carries the same idea on the web page.""")
explainer("kutta_condition", "Which circulation does a wing choose?",
          "sliding Γ and watching the rear stagnation point crawl to the edge while the edge-speed readout falls from "
          "unbounded to a finite number — and rises again past the Kutta value — is the hypothesis as an experiment.",
          ["Drag Γ/Γ_K from 0 to 1 and watch the amber dot reach the tail.",
           "Press 'too much': where is the rear stagnation point now?",
           "Set α = −β: what circulation does the edge ask for?",
           "Open the Derivation tab at step 3 and find the angle α + β in the circle view."])
explainer("lift_curve_stall", "Does more angle always mean more lift?",
          "dragging α along the curve while the section view goes from attached to separated, with camber shifting the line "
          "and the regime chips moving the stall point, lets you find slope, shift and limit by experiment.",
          ["Raise the camber: does the slope change?",
           "Find the zero-lift angle and compare it with −β.",
           "Put the flap down and read the new stall speed.",
           "Go just past stall: how far below the line is the model?"])
note("N33 [C] · N35 [C]", r"""
- **N33 [C]** Zhukhovsky sections are not used on aircraft: a cusp cannot be built, and the camber peaks at mid-chord. The
  Kármán–Trefftz map (a finite edge angle) is the standard relative; this book does not develop it.
- **N35 [C]** Plotting a cambered section with its streamlines, and moving a known lift coefficient to a second incidence by
  $C_L(\alpha_2)=C_L(\alpha_1)+2\pi(\alpha_2-\alpha_1)$: `scripts/ch14_zhukhovsky_airfoil.py` does both with our parameters.""")
whatif(r"""…we ask *where on the surface* this lift is produced? C07 looks at the pressure.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C07", r"Surface speed and pressure: $C_p=1-\dfrac{\lvert dw/dz\rvert^2}{U^2}$ on the contour (ours)",
     "Is a wing pushed up from below or pulled up from above — and where?")
problem(r"""
Put your hand out of a car window, palm tilted up. It feels pushed from below. A real wing is mostly *pulled* from above: the
air speeds up over the curved upper surface, its pressure drops, and that suction — strongest just behind the nose — does most
of the lifting. The same picture shows the wing's weak spot: behind the suction peak the pressure has to climb back (an
*adverse pressure gradient*), and a boundary layer made to climb too steeply gives up. That is stall. **Not in the book —
ours:** the book shows this pressure figure; we compute it.""")
idea(r"""
        ↑ ↑ ↑                                        −C_p
       ↑↑ ↑  ↑   ↑    ↑     suction (arrows point       │╲
     ╭──────────────────────╮   away from the surface)  │ ╲___ upper surface
    (                        ───────                    │     ╲______
     ╰───────────────────────╯                          │  area = normal force
       ↑    ↑      ↑       pressure pushes up, weakly   │_________________  lower surface
                                                        └──────────────────── x/c
""")
recap("R03", "The pressure coefficient", r"""
$C_p=\dfrac{p-p_\infty}{\tfrac12\rho U^2}$: +1 at a stagnation point, 0 where the speed equals $U$, negative where the flow is faster. **New here:**
aerodynamicists plot $-C_p$ upward, so suction is "up", and the two surfaces make a loop.""", where="Ch. 6 §6.3", core_id="C07")
recap_md("Bernoulli", "Ch. 4 §4.9", r"""
Along a streamline of steady ideal flow $p+\tfrac12\rho q^2$ is constant, so $C_p=1-q^2/U^2$ with $q$ the local speed.""")
nb.code(r"""
print(POT.pressure_coefficient(speed=45.0, U=30.0))          # a point where the air moves at 1.5 U: C_p = 1 − 1.5² = −1.25
""")
remind("C07")
maths(r"""
1. Points of the circle: $\zeta=\mu+ae^{it}$, with $t=-\beta$ at the tail point.
2. Circle velocity $dw/d\zeta$ there (derivation D06) with the Kutta circulation $\Gamma=4\pi Ua\sin(\alpha+\beta)$ (14.11).
3. Rescale by the map: $\dfrac{dw}{dz}=\dfrac{dw}{d\zeta}\,\dfrac{\zeta^2}{\zeta^2-b^2}$ (14.9).
4. Speed: $q=\lvert dw/dz\rvert$.
5. Pressure coefficient: $C_p=1-q^2/U^2$ (Bernoulli).
6. At $t=-\beta$ the 0/0 is replaced by its limit $U(b/a)\cos(\alpha+\beta)$ (N29, trap T8).""")
D("D10", "The normal-force coefficient is the area of the pressure loop", ref="14.12")
nb.worked_example("pressure by hand", r"""
$U$ = 30 m/s, ρ = 1.225 kg/m³: dynamic pressure 551 Pa.

1. A point where the speed is 45 m/s: $C_p$ = 1 − (45/30)² = −1.25, i.e. 689 Pa *below* ambient.
2. The stagnation point: $C_p$ = 1, 551 Pa above.
3. A toy loop: $C_p$ = −0.6 all along the top, +0.4 all along the bottom: $C_y$ = ∫(0.4 − (−0.6)) d(x/c) = 1.0 — and 60 % of it
   is suction.
4. Force normal to a 1 m chord: 1.0 × 551 × 1 = 551 N per metre of span.""")
nb.code(r"""
Sf = AF.zhukhovsky_surface(**SEC, n=2000)                    # speed and pressure on the contour of our section (Kutta circulation), 2000 steps round the circle
up, lo = Sf["upper"], Sf["lower"]                            # the two surfaces, each with x/c increasing
Cy = AF.normal_force_coefficient(xc_upper=up["xc"], cp_upper=up["cp"], xc_lower=lo["xc"], cp_lower=lo["cp"])    # area of the pressure loop
CL_ideal = AF.lift_from_cp(xc_upper=up["xc"], cp_upper=up["cp"], xc_lower=lo["xc"], cp_lower=lo["cp"], alpha=sec["alpha"], suction="ideal")   # C_y/cos α: no drag
CL_none = AF.lift_from_cp(xc_upper=up["xc"], cp_upper=up["cp"], xc_lower=lo["xc"], cp_lower=lo["cp"], alpha=sec["alpha"], suction="none")     # C_y cos α: no leading-edge suction
print(f"suction peak C_p = {Sf['cp_min']:.2f} at x/c = {100*Sf['xc_cp_min']:.2f} % (upper surface)")   # −2.24 at 0.45 % of the chord
print(f"loop area C_y = {Cy:.3f}; C_L = {CL_ideal:.3f} (ideal), {CL_none:.3f} (no leading-edge suction)")   # 1.054; 1.058; 1.050
print(f"2Γ/(Uc) = {2*Sf['Gamma_cw']/(sec['U']*Sf['chord_x']):.3f}")                                    # 1.058: the lift coefficient from the circulation
print(f"∮u_s ds round the contour = {AF.surface_circulation(surface=Sf):.2f} m²/s (clockwise)")        # 15.96 = Γ_K: the number promised in N10
""", explain=r"""
1. `AF.zhukhovsky_surface` carries out the six steps above on 2000 points of the contour and splits them into `upper` and
   `lower`.
2. The suction peak is −2.24 and sits at only 0.45 % of the chord: right behind the nose, on the upper surface.
3. `AF.normal_force_coefficient` is the area of the loop, $C_y$ = 1.054.
4. The loop gives the force normal to the *chord*; the lift is normal to the *stream* (trap T5). Dividing by cos α (ideal
   flow, no drag) gives 1.058 — exactly $2\Gamma/(Uc)$ from the circulation.
5. `AF.surface_circulation` adds up the surface speed round the contour: 15.96 m²/s, the Kutta circulation. The faster flow
   over the top *is* the circulation (N10).""")
scratch(r"""
U, b, al, beta = sec["U"], sec["b"], sec["alpha"], sec["beta"]   # speed [m/s], map constant [m], angle of attack and camber angle [rad]
t = -beta + 2*np.pi*np.arange(1, 2000)/2000                  # 1999 angles round the circle, leaving out the tail point itself (0/0 there)
zeta = mu + a*np.exp(1j*t)                                   # points of the circle
dw = U*(np.exp(-1j*al) - a**2*np.exp(1j*al)/(zeta - mu)**2) + 1j*GK/(2*np.pi*(zeta - mu))   # circle velocity with the Kutta circulation
q = np.abs(dw*zeta**2/(zeta**2 - b**2))                      # Eq. (14.9), then the speed on the section [m/s]
cp = 1 - (q/U)**2                                            # Bernoulli: pressure coefficient
assert np.allclose(cp, Sf["cp"][1:], rtol=1e-9)              # same numbers → the library does exactly this
print(f"C_p from {cp.min():.3f} (suction peak) to {cp.max():.3f} (stagnation point)")   # −2.24 … 1.000
""", r"""
The six steps of "the maths" as five lines of complex arithmetic on a whole array of points, compared with
`AF.zhukhovsky_surface`.""")
nb.figure(r"""
fig = fig_pressure_two_ways(S=AF.zhukhovsky_surface(**SEC, n=400))   # left: pressure arrows on the section; right: −C_p against x/c for both surfaces
plt.show()                                                   # display
""",
    see="Left: our section with arrows normal to the surface — pointing outward where the pressure is below ambient "
        "(suction), inward where it is above. Right: $-C_p$ against $x/c$ for both surfaces, with the loop between them shaded.",
    read="The tall thin spike at the nose is the suction peak; behind it the upper curve slopes down all the way to the tail "
         "— the pressure is *rising* there, which the boundary layer has to climb. The lower surface pushes up only weakly. "
         "The shaded area is the normal-force coefficient, and most of it lies above the axis: the wing is pulled up.",
    change="…α were 10°? The peak would pass −7 and the slope behind it would steepen: closer to stall.")
nb.plotly(r"""
def cp_loop(alpha_deg):                                      # the two pressure curves for one angle of attack
    Sx = AF.zhukhovsky_surface(U=sec["U"], b=sec["b"], e=sec["e"], beta=sec["beta"], alpha=np.deg2rad(alpha_deg), n=400)   # surface pressure at this α
    return {"upper surface": (Sx["upper"]["xc"], -Sx["upper"]["cp"]), "lower surface": (Sx["lower"]["xc"], -Sx["lower"]["cp"])}   # −C_p, so suction is up
figF5 = slider_figure(cp_loop, "angle of attack α", np.linspace(-4.0, 12.0, 17), unit="°", xlabel="x/c [–]", ylabel="−C_p [–]",
                      yrange=[-1.1, 6.0], active=9, title="The loop between the two curves is the normal force")   # 17 slider positions, starting at 5°
figF5.show()                                                 # display
rows = []                                                    # the loop area at four angles
for A_deg in (-4.0, 0.0, 5.0, 10.0):                         # angles of attack [degrees]
    Sx = AF.zhukhovsky_surface(U=sec["U"], b=sec["b"], e=sec["e"], beta=sec["beta"], alpha=np.deg2rad(A_deg), n=2000)   # surface pressure
    Cy_ = AF.normal_force_coefficient(xc_upper=Sx["upper"]["xc"], cp_upper=Sx["upper"]["cp"], xc_lower=Sx["lower"]["xc"], cp_lower=Sx["lower"]["cp"])   # loop area
    rows.append({"α [°]": A_deg, "α [rad]": np.deg2rad(A_deg), "loop area C_y": Cy_, "2Γ/(Uc)": 2*Sx["Gamma_cw"]/(sec["U"]*Sx["chord_x"]), "peak −C_p": -Sx["cp_min"]})   # collect
print(pd.DataFrame(rows).round(3).to_string(index=False))    # zero net area at α = −β = −4°
""", explain=r"""
1. `cp_loop` returns $-C_p$ on both surfaces for one angle of attack (degrees on the slider, radians in the call).
2. The table lists the loop area and $2\Gamma/(Uc)$ at four angles: they differ only by the factor cos α.""")
nb.figure_notes(
    see="The pressure loop of our section as the angle of attack goes from −4° to 12°.",
    read="At α = −4° (= −β) the two curves cross and the net area is zero: no lift. As α grows the suction peak near the nose "
         "shoots up (the vertical axis is cut at 6) and the loop opens. The area grows in step with $2\\Gamma/(Uc)$.",
    change="…the section were a flat plate? The peak would be infinite at the leading edge for every α ≠ 0 — a rounded nose "
           "is what keeps it finite (N12 below).")
recap("R10", "The Blasius force", r"""
$D-iL=\dfrac{i\rho}{2}\displaystyle\oint\Big(\dfrac{dw}{dz}\Big)^2dz$ (6.60): the force from a contour integral of the squared complex velocity. **New here:** with
the stream inclined, the force comes out perpendicular to the *stream* — size $\rho U\Gamma$, no drag.""", where="Ch. 6 §6.5", core_id="C07")
nb.code(r"""
fb = AF.zhukhovsky_force(**SEC, rho=1.225, method="blasius")    # force from the Blasius contour integral, Eq. (6.60)
fp = AF.zhukhovsky_force(**SEC, rho=1.225, method="pressure")   # force from adding up the surface pressure
for name, f_ in (("Blasius ", fb), ("pressure", fp)):         # the two routes share no line of code
    print(f"{name}: Fx = {f_['Fx']:7.2f}, Fy = {f_['Fy']:7.2f} N/m → L = {f_['L']:.2f}, D = {f_['D']:.1e} N/m; centre of pressure at {100*f_['xcp_over_c']:.1f} % chord")   # −51.1, 584.4 → 586.6, ≈ 0; 35.7 %
""", explain=r"""
**What does this show?** Two routes that never share a line of code agree: the force has components $F_x$ = −51.1 and $F_y$ =
584.4 N/m along and across the chord, which resolve (N08) into a lift of 586.6 N/m $=\rho U\Gamma$ and no drag. The force leans
*forward* of the chord normal by exactly α. It acts at 35.7 % of the chord.""")
note("N12 [B]", r"""
The pressure on a flat plate is normal to the plate, so why is the force normal to the stream? Because the flow whipping round
the sharp *leading* edge pulls the plate forward: the **leading-edge suction** $S=L\sin\alpha=\pi\rho U^2c\sin^2\alpha$, with the normal
force $N=L\cos\alpha$ (ours; stated, not derived). A rounded nose spreads that suction over a real surface — the −2.24 peak above
— and keeps the boundary layer attached.""")
nb.code(r"""
fpl = AF.flat_plate_forces(rho=1.225, U=30.0, c=1.0, alpha=sec["alpha"])   # flat plate of chord 1 m at 5° in air at 30 m/s
print(f"L = {fpl['L']:.1f}, N = {fpl['N']:.1f}, S = {fpl['S']:.1f} N/m")     # 301.9, 300.7, 26.3 N/m
cv = ch14.edge_flow_suction(rho=1.225, Uo=30.0*np.sin(sec["alpha"]), a=1.0)  # control volume round the edge: near the nose the speed is U sin α · sqrt(c/r)
print(f"edge suction from the control volume = {cv['closed_form']:.1f} N/m") # 26.3 N/m: the same S
""", explain=r"""
**What does this show?** The plate's lift is 301.9 N/m; its normal force is 300.7 N/m and the leading-edge suction 26.3 N/m.
A momentum balance on a tiny circle round the edge (`ch14.edge_flow_suction`) gives the same 26.3 N/m.""")
note("N13 [B]", r"""
**Stall.** At high α the suction-side boundary layer cannot climb the pressure rise behind the peak; it *separates* ahead of
the trailing edge (where the wall shear stress $\tau_0=0$), the Kutta condition fails and the lift collapses.""")
recap_md("Thwaites' method", "Ch. 9 §9.6", r"""
It marches a laminar boundary layer along a given outer speed $U_e(x)$ and reports where the wall shear stress $\tau_0$ reaches
zero — the separation point.""")
nb.figure(r"""
angles_deg = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0]) if not FAST else np.array([0.0, 5.0, 10.0])   # angles of attack to try [degrees]
xsep = [ch14.suction_side_separation(alpha=np.deg2rad(A_), e=sec["e"], beta=sec["beta"], Re_c=1e6, n=120)["x_sep_over_c"] for A_ in angles_deg]   # laminar separation point on the suction side, as a fraction of the chord, at chord Reynolds number 1e6
fig, ax = plt.subplots(figsize=(6.2, 3.4))                   # one panel
ax.plot(angles_deg, 100*np.array(xsep, dtype=float), "o-", color=COLORS["rose"])   # separation position against α
ax.set(xlabel="angle of attack α [°]", ylabel="separation point [% of chord]", title="Qualitative model: laminar separation point on the ideal pressure distribution")   # labels
plt.show()                                                   # display
print(pd.DataFrame({"α [°]": angles_deg, "α [rad]": np.deg2rad(angles_deg), "x_sep/c [%]": 100*np.array(xsep, dtype=float)}).round(3).to_string(index=False))   # the numbers
""",
    see="Where a *laminar* boundary layer, marched by Thwaites' method along the ideal-flow surface speed, would separate on "
        "the upper surface, for several angles of attack.",
    read="As α grows, the separation point moves forward toward the suction peak: the pressure rise behind the peak gets "
         "steeper and the layer gives up earlier. That forward march is the mechanism of stall.",
    change="…the layer were turbulent? It would survive a much steeper rise and separate far later (N72 in C15).")
confusion(r"""(trap T16) this model **cannot predict the maximum lift**: the real layer turns turbulent, and once it separates the
outer flow itself changes. It shows the trend only.""")
explainer("airfoil_pressure_lift", "Is a wing pushed up or pulled up?",
          "turning α and seeing the suction peak sharpen, the stagnation point slide under the nose and the loop area grow in "
          "step with the lift readout connects three pictures of one flow; clicking a point traces Bernoulli's arithmetic there.",
          ["Raise α from 0° to 8° and watch the peak.",
           "Click the stagnation point: what is $C_p$ there?",
           "Compare the two bars 'upper suction' and 'lower pressure'.",
           "Choose 'flat plate': what happens to the peak?"])
whatif(r"""…the wing has ends? Everything so far was per metre of an endless span. §14.6.""")


# =====================================================================================================================
# A.6 §14.6 Elementary Lifting Line Theory for Wings of Finite Span — C08 … C14
# =====================================================================================================================
nb.section("14.6", "Elementary Lifting Line Theory for Wings of Finite Span", intro=r"""
**What is this section about?** A real wing has tips. The circulation must fall to zero there, a vortex cannot simply end, so
vorticity streams off the wing and trails behind. That wake makes the air at the wing sink, the sinking tilts the lift
backward, and the wing has a drag even in ideal flow. Seven steps: the trailing filaments, their downwash, the induced drag, a
sine series, one linear system for any wing, the forces from its coefficients, and the best loading.

> ⚠️ **Axes change here (trap T5).** From here on $y$ runs **along the span** ($-s/2$ at the port tip … $+s/2$ at the starboard
> tip) and the downwash $w$ is positive **downward**. The letter $w$ is now a velocity, no longer the complex potential.""")
core("C08", r"Bound vortex and trailing filaments: $d\Gamma_{\text{trail}}=-\dfrac{d\Gamma}{dy}\,dy$",
     "The circulation round a wing must vanish at the tips. Where does it go?")
problem(r"""
On a damp day you can see two white cords streaming back from an airliner's wing tips or flap ends: tip vortices. Below a
lifting wing the pressure is high, above it low; at the tips the air leaks round from below to above. That leak is the visible
end of something that happens all along the trailing edge: wherever the wing carries a little less circulation than it did a
step further inboard, the difference leaves as a thin vortex trailing downstream.""")
idea(r"""
   top view (flight direction up the page)

            Γ(y):        _.-""""""-._              the bound circulation: largest at mid-span, zero at the tips
                       .'            '.
   lifting line   ════════════════════════
                  ║  ║  ║          ║  ║  ║         many weak horseshoe vortices, one inside the next
                  ║  ║  ╚══════════╝  ║  ║         = one bound vortex of varying strength
                  ║  ╚════════════════╝  ║         + a sheet of trailing filaments
                  ╚══════════════════════╝
                  ↓  ↓  ↓          ↓  ↓  ↓         (the legs trail downstream)
""")
remind("C08")
note("N36 [B]", r"""
Air flows outward under the wing and inward over it. Behind the trailing edge the two streams meet with different spanwise
velocities: a surface across which the tangential velocity jumps — a **vortex sheet** — of opposite sense on the two half-wings.""")
recap_md("vortex sheet", "Ch. 5 §5.8", r"""
A jump in tangential velocity across a surface is a sheet of vorticity whose strength per unit length equals the jump,
$\gamma_{\text{sheet}}=\Delta u_{\text{tangential}}$.""")
nb.figure(r"""
fig = draw_tip_flow_rear_view()                              # rear view of a lifting wing: high pressure below, low above, the leak round each tip
plt.show()                                                   # display
""",
    see="A wing seen from behind: high pressure under it, low pressure over it, and the air moving outward below, inward "
        "above and round each tip.",
    read="The two tip flows turn in opposite senses — the left one clockwise, the right one counterclockwise as drawn — and "
         "between the tips the air they set in motion goes *down*.",
    change="…the wing carried no lift? No pressure difference, no leak, no trailing vortices.")
recap("R11", "A vortex cannot end in the fluid", r"""
Helmholtz: the strength of a vortex tube is the same all along it, $\Gamma_{\text{tube}}=\text{const}$, so a tube must close on itself or end on
a boundary. **New here:** the wing's vortex system is one closed loop — the bound vortex on the wing, two trailing vortices,
and far behind the starting vortex of C03.""", where="Ch. 5 §5.3", core_id="C08")
note("N39 [B]", r"""
**The three assumptions of lifting-line theory** (later blocks cite them by number):

1. the aspect ratio is large, so each spanwise strip sees a two-dimensional flow — rules out stubby and delta wings;
2. the wing is one straight bound vortex of strength $\Gamma(y)$, the **lifting line** — rules out sweep;
3. the trailing filaments are straight, parallel to the stream, and do not act on one another — rules out the roll-up we are
   about to watch.""")
note("N40 [B]", r"""
The **bound vortex** moves with the wing, is strongest at mid-span and vanishes at the tips, $\Gamma(\pm s/2)=0$.""")
nb.code(r"""
print(LL.circulation_series(Gn=np.array([12.0, 0.0, 1.0]), gamma=np.array([0.0, np.pi/2, np.pi])))   # a loading made of two sine shapes (C11), at port tip, mid-span, starboard tip: 0, 11, 0 m²/s
""")
D("D11", "The strength of a shed filament", ref="14.13")
nb.worked_example("a staircase of circulation", r"""
On the starboard half the bound circulation is 10 m²/s at the root, then 8, then 5, then 0 at the tip.

1. Where it drops from 10 to 8 a filament of strength 2 leaves.
2. From 8 to 5: strength 3.
3. From 5 to 0 at the tip: strength 5 — the strongest.
4. Total shed on this half: 2 + 3 + 5 = 10 = the root circulation $\Gamma_0$.
5. The port half sheds the same amounts with the opposite sense of rotation.

The trailing vorticity is concentrated near the tips because that is where Γ falls fastest.""")
nb.code(r"""
g9 = np.linspace(0.05, np.pi - 0.05, 9)                      # nine stations along the span, as span angles γ (0 = port tip, π = starboard tip; C11 explains γ)
y9 = LL.span_position(gamma=g9, s=ac["s"])                   # the same stations in metres: y = −(s/2) cos γ
sheet = LL.trailing_sheet_strength(Gn=np.array([14.17]), gamma=g9, s=ac["s"])   # −dΓ/dy [m/s] for an elliptic loading with mid-span circulation 14.17 m²/s on our 12 m span
print(pd.DataFrame({"y [m]": y9, "−dΓ/dy [m/s]": sheet}).round(2).T.to_string(header=False))   # antisymmetric: zero at mid-span, largest toward the tips
tv = ch14.tip_vortex_estimate(W=ac["W"], rho=1.225, U=ac["U_cruise"], s=ac["s"])   # what our aircraft's wake rolls up into at 55 m/s
print(f"Γ0 = {tv['Gamma0']:.2f} m²/s, spacing of the two tip vortices = {tv['spacing']:.2f} m, they sink at {tv['descent_speed']:.2f} m/s")   # 14.17, 9.42 (= πs/4), 0.24
""", explain=r"""
1. `LL.trailing_sheet_strength` returns $-d\Gamma/dy$, the strength per metre of span of the trailing sheet. It is antisymmetric
   (negative on the port half, positive on the starboard half), zero at mid-span and largest toward the tips.
2. `ch14.tip_vortex_estimate` (ours) gives the mid-span circulation our aircraft needs at cruise, 14.17 m²/s clockwise, and what
   the sheet rolls up into (N37 below): two vortices 9.42 m apart sinking at 0.24 m/s.
3. The span angle γ is a convenient way to place stations; C11 introduces it properly.""")
scratch(r"""
s = ac["s"]                                                  # span [m]
gf = np.linspace(0.05, np.pi - 0.05, 2001)                   # a fine set of span angles
yf = -(s/2)*np.cos(gf)                                       # the stations in metres
Gf = 14.17*np.sin(gf)                                        # the elliptic loading Γ(y) [m²/s]
mine = -(Gf[2:] - Gf[:-2])/(yf[2:] - yf[:-2])                # −dΓ/dy by a centred difference: (next − previous)/(distance between them)
assert np.allclose(mine, LL.trailing_sheet_strength(Gn=np.array([14.17]), gamma=gf, s=s)[1:-1], rtol=1e-5)   # same numbers → the library does exactly this
print(f"−dΓ/dy half-way out on the starboard side (y = {yf[1334]:.2f} m): {mine[1333]:.3f} m/s")   # positive on the starboard half
""", r"""
The slope of the loading is estimated by finite differences on a fine grid and compared with the tested function (which
differentiates the sine series exactly).""")
nb.figure(r"""
fig = draw_horseshoe_staircase(N=7)                          # seven nested horseshoe vortices and the staircase Γ(y) they add up to
plt.show()                                                   # display
yy = np.linspace(-6.0, 6.0, 401)[1:-1]                       # stations along the 12 m span [m], tips excluded (the slope is infinite there)
Gam = LL.elliptic_circulation(y=yy, s=12.0, Gamma1=14.17)    # Eq. (14.23): the elliptic loading [m²/s]
fig, ax = plt.subplots(figsize=(6.4, 3.2))                   # one panel
ax.plot(yy, Gam, color=COLORS["amber"], label="bound circulation Γ(y) [m²/s], clockwise")          # the loading
ax.plot(yy, np.clip(-np.gradient(Gam, yy), -12, 12), color=COLORS["blue"], label="shed strength −dΓ/dy [m/s] (cut at ±12)")   # the trailing sheet
ax.set(xlabel="spanwise position $y$ [m]", ylabel="Γ and −dΓ/dy", title="The wing sheds vorticity where its circulation falls")   # labels
ax.legend(fontsize=8); plt.show()                            # legend and display
""",
    see="First: seven nested horseshoes seen from above, with the staircase of circulation they add up to. Second: the "
        "smooth elliptic loading (amber) and the strength of the sheet it sheds (blue).",
    read="Each horseshoe adds one step to the staircase and two trailing legs. In the smooth limit the legs become a sheet "
         "whose strength is the slope of the loading: nothing leaves at mid-span, almost everything leaves near the tips.",
    change="…the loading were flat-topped? The filaments would crowd even closer to the tips.")
note("N37 [B]", r"""
Assumption (iii) is not what happens: each half of the sheet rolls up into one **tip vortex** of circulation
$\Gamma_{\text{tip}}=\Gamma_0$. **Not in the book — ours:** the two are a distance $b'=\dfrac1{\Gamma_0}\displaystyle\int\Gamma\,dy$ apart ($\tfrac\pi4s$ for the elliptic
loading) and sink at $\Gamma_0/2\pi b'$. This is why airports keep a small aircraft minutes behind a large one, and what a
contrail makes visible.""")
recap_md("point vortices move with the velocity the others induce", "Ch. 5 §5.6", r"""
A vortex is carried by the velocity all the *other* vortices induce at its position, $d\mathbf x_k/dt=\sum_{j\ne k}\mathbf u_j(\mathbf x_k)$;
the total circulation and the impulse $\sum\Gamma_ky_k$ do not change.""")
nb.animation(r"""
roll = ch14.wake_rollup(Gn=np.array([14.17]), s=12.0, n=(200 if not FAST else 80), n_t=41)   # the sheet of our aircraft cut into point vortices and followed in time (time stands for distance behind the wing)
stops = [0, 5, 10, 20, 40]                                   # five stop-frames out of the 41 saved times
fig, ax = plt.subplots(figsize=(6.4, 3.6))                   # rear view
sc = ax.scatter(roll["y"][0], roll["z"][0], s=6, c=np.where(roll["sigma"] > 0, COLORS["blue"], COLORS["rose"]))   # one dot per vortex, coloured by its sense of rotation
(cen,) = ax.plot([], [], "k+", ms=12)                        # the centroid of each half
ax.set(xlim=(-8, 8), ylim=(-4.5, 1.5), xlabel="spanwise position $y$ [m]", ylabel="height $z$ [m]")   # fixed axes
ax.set_aspect("equal")                                       # true proportions
fig.canvas.draw(); fig.set_layout_engine("none")             # lay the figure out once, then freeze it
def update(j):                                               # draw stop-frame j
    k = stops[j]                                             # index of the saved time
    sc.set_offsets(np.c_[roll["y"][k], roll["z"][k]])        # move the vortices
    cen.set_data([-roll["centroid_y"][k], roll["centroid_y"][k]], [roll["centroid_z"][k], roll["centroid_z"][k]])   # the two centroids
    ax.set_title(f"Qualitative past the first turn — t = {roll['t'][k]:.1f} s: centroids {2*roll['centroid_y'][k]:.2f} m apart, {-roll['centroid_z'][k]:.2f} m down")   # readout
show_animation(animate(update, frames=len(stops), fig=fig, interval=900), player="frames")   # a frame player: stop at each of the five times
print(f"theory (ours): spacing b' = {roll['spacing_theory']:.2f} m, descent speed = {roll['descent_speed_theory']:.3f} m/s; impulse drift over the run = {np.ptp(roll['impulse'])/abs(roll['impulse'][0]):.1e}")   # the invariants
""", explain=r"""
1. `ch14.wake_rollup` cuts the trailing sheet of an elliptic loading into point vortices (with a small core so that the
   computation stays smooth) and moves each with the velocity the others induce. Time stands for distance behind the wing.
2. Five of the saved times are shown as stop-frames; blue and rose are the two senses of rotation; the + marks are the
   centroids of the two halves.
3. The printed line gives the theoretical spacing $b'=\tfrac\pi4s$ = 9.42 m and descent speed, and how little the impulse
   $\sum\Gamma_ky_k$ drifted.""")
nb.figure_notes(
    see="The wake seen from behind at five times. A flat row of vortices curls up at both ends into two spirals that sink.",
    read="The centroids (+) stay the same distance apart — 9.42 m, the invariant $b'$ — while they move down. **The turns of "
         "the spiral are qualitative**: their number depends on how many vortices and what core size the computation uses. "
         "What is robust is the spacing, the total circulation of each half and the sinking.",
    change="…the aircraft were heavier or slower? More circulation (it scales as $W/U$): stronger vortices that sink faster.")
note("N69 [C]", r"""
The picture of a bound vortex with trailing vortices is Lanchester's; the usable theory is Prandtl's (see N17 in C03).""")
whatif(r"""…we stand on the wing and ask what all those filaments do to the air there? C09.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C09", r"Downwash: $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int_{-s/2}^{s/2}\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13)",
     "How fast is the air sinking where the wing meets it?")
problem(r"""
A wing holds the aircraft up by throwing air down — Newton's third law. Some of that downward motion is already present *at
the wing itself*, caused by the trailing vortices behind it. Every filament swirls the air round its own axis; between the
tips all of them push down. Add up their pushes at one point of the wing and you have the downwash there.""")
idea(r"""
   top view                         lifting line
   ─────────────●───────────────────────○──────────────      ● station y₁ where we want the downwash
                y₁    |y − y₁| = d      y                    ○ a filament leaves at y …
                                        ║
                                        ║   … and runs downstream for ever
                                        ▼

   an infinite line vortex would give   Γ/(2π d)
   this one starts at the wing, so the wing sees half:   Γ/(4π d)
""")
recap_md("the filament law", "Ch. 5 §5.5", r"""
$d\mathbf u=\dfrac{\Gamma\,d\mathbf l}{4\pi}\times\dfrac{\mathbf x-\mathbf x'}{\lvert\mathbf x-\mathbf x'\rvert^3}$ (5.17): every piece of a vortex filament induces a velocity round itself; added along a
straight piece it gives the speed $\dfrac{\Gamma}{4\pi d}(\cos\theta_a-\cos\theta_b)$ at distance $d$ from the line, with $\theta_a$, $\theta_b$ the angles to its
two ends.""")
recap("R12", "Half an infinite line vortex", r"""
For a filament that starts at the lifting line and runs downstream for ever the two angles are 0 and 90°: speed $\Gamma/(4\pi d)$,
half the value $\Gamma/(2\pi d)$ of a line infinite both ways. With the strength of the shed filament from the derivation in C08:
$dw(y_1)=\dfrac{-(d\Gamma/dy)\,dy}{4\pi(y-y_1)}$.""", where="Ch. 5 §5.5", core_id="C09")
nb.code(r"""
print(BS.segment_speed(d=1.0, theta_a=0.0, theta_b=np.pi/2, Gamma=4*np.pi))   # a half-line of strength 4π at distance 1 m: 1.0 m/s
print(BS.segment_speed(d=1.0, theta_a=0.0, theta_b=np.pi, Gamma=4*np.pi))     # the full infinite line: 2.0 m/s — twice as much
""")
slip(3, r"""that this velocity is induced “at location $y$ (< $y_1$)”""",
     r"""it is induced **at $y_1$**, and it is downward for a filament at $y>y_1$ on the starboard half, as the book's own figure caption says.""")
P("P344", "Cauchy principal value", r"""
An integrand like $1/(x-c)$ is +∞ on one side of $c$ and −∞ on the other. The two infinite areas cancel if we cut out a
*symmetric* gap $(c-\epsilon,\ c+\epsilon)$ and let ε → 0: that limit is the principal value. `scipy.integrate.quad` (primed in Ch. 3,
P87) computes it with `weight="cauchy"`.""", r"""
val = quad(lambda x: 1.0, -1, 1, weight="cauchy", wvar=0.3)[0]   # principal value of the integral of 1/(x - 0.3) from -1 to 1
print(val, np.log(0.7/1.3))                                      # both -0.619
""")
remind("C09")
D("D12", "The downwash at the lifting line", ref="14.13")
nb.worked_example("the simplest wing: one horseshoe", r"""
**N60 [B]** One bound vortex of constant strength Γ = 4π ≈ 12.6 m²/s and span $s$ = 2 m, with one trailing leg from each tip.

1. Each leg is a half-line at distance $s/2$ = 1 m from the centre: it induces Γ/(4π × 1) = 1 m/s there.
2. Both legs push down at the centre: $w(0)$ = 2 m/s.
3. In symbols $w(0)=\dfrac{\Gamma}{\pi s}$ (ours) = 4π/(2π) = 2 ✓.
4. Move the station toward a tip, to $y_1$ = 0.9 m: the near leg is 0.1 m away and gives 10 m/s, the far leg 1.9 m away gives
   0.53 m/s — the downwash grows without limit at the tips.""")
nb.code(r"""
print(f"horseshoe (Γ = 10 m²/s, s = 12 m): downwash at the centre = {LL.horseshoe_downwash(x=0.0, y=0.0, Gamma=10.0, s=12.0):.4f} m/s")   # Γ/(πs) = 0.2653, positive down
print("the same as a velocity vector (x downstream, y along the span, z up):", np.round(LL.horseshoe_velocity(x=0.0, y=0.0, z=0.0, Gamma=10.0, s=12.0), 4))   # (0, 0, −0.2653): with z up, a downwash is negative
hw = LL.horseshoe_wing(Gamma=10.0, s=12.0, rho=1.225, U=30.0)   # lift and 'induced drag' of the horseshoe wing
print(f"horseshoe wing: L = {hw['L']:.0f} N; its drag over the best possible wing's = {hw['ratio']:.2f} ({hw['note']})")   # 0.5 — an artefact
ell = lambda y: np.sqrt(1 - (2*y/12.0)**2)                    # an elliptic loading with mid-span circulation 1 m²/s on a 12 m span
print(f"elliptic loading: w(0) = {LL.downwash(y1=0.0, Gamma_fn=ell, s=12.0):.5f} m/s, w(3 m) = {LL.downwash(y1=3.0, Gamma_fn=ell, s=12.0):.5f} m/s")   # 0.04167 = Γ1/2s at both stations
""", explain=r"""
1. `LL.horseshoe_downwash` returns the downwash (positive down) of one horseshoe: Γ/(πs) = 0.2653 m/s at the centre.
2. `LL.horseshoe_velocity` returns true velocity components with $z$ **up**, so the same downwash appears as −0.2653.
3. `LL.horseshoe_wing` reports a drag ratio of one half — see the callout below.
4. `LL.downwash` evaluates the principal-value integral $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13) for any loading given as a
   function. For the elliptic loading it returns the same number at both stations: 0.04167 = 1/(2 × 12).""")
confusion(r"""(trap T10) the horseshoe's induced drag, $C_{D_i}=C_L^2/(2\pi\Lambda)$ (ours), looks half that of the best possible
wing. It is an artefact: it uses the centre downwash everywhere, while the true downwash of a constant Γ is infinite at the
tips.""")
scratch(r"""
def staircase_downwash_mine(y1, ye, Gs):                     # downwash at y1 of a staircase loading: steps Gs between the edges ye
    w_ = 0.0                                                 # running sum [m/s]
    for k in range(len(ye)):                                 # one trailing filament at every edge
        left = Gs[k-1] if k > 0 else 0.0                     # bound strength just inboard (port side) of this edge; zero outside the wing
        right = Gs[k] if k < len(Gs) else 0.0                # bound strength on the other side
        w_ += (right - left)/(4*np.pi*(y1 - ye[k]))          # half-line formula of R12 with the jump of Γ across the edge
    return w_                                                # downwash, positive down
ye, Gs = LL.staircase_from_loading(Gamma_fn=ell, s=12.0, N=21)   # 21 horseshoes approximating the elliptic loading: 22 edges, 21 step values
assert np.isclose(staircase_downwash_mine(0.0, ye, Gs), LL.staircase_downwash(y1=0.0, y_edges=ye, Gamma_steps=Gs))   # same number → the library does exactly this
for N in (1, 3, 7, 21, 101):                                 # more and more horseshoes
    ye, Gs = LL.staircase_from_loading(Gamma_fn=ell, s=12.0, N=N)   # the staircase with N steps
    print(f"N = {N:3d}: mid-span downwash = {staircase_downwash_mine(0.0, ye, Gs)/(1/24.0):.4f} × Γ1/2s")   # 0.637, 0.961, 0.991, 0.998, 0.9998 → 1
""", r"""
1. `staircase_downwash_mine` adds the half-line contribution of each trailing filament; the filament at an edge has the
   strength of the jump in bound circulation across that edge.
2. The assert compares with `LL.staircase_downwash` for 21 horseshoes.
3. The loop shows the sum converging to the value of the integral: 0.637 with one horseshoe (= 2/π), 0.998 with 21.
   **The downwash integral is the limit of this sum.**""")
nb.plotly(r"""
y_st = np.linspace(-5.93, 5.93, 60)                          # 60 stations along the 12 m span [m]
def staircase(N):                                            # loading and downwash of N nested horseshoes
    ye_, Gs_ = LL.staircase_from_loading(Gamma_fn=ell, s=12.0, N=int(N))   # edges and step values
    w_ = np.array([LL.staircase_downwash(y1=y1, y_edges=ye_, Gamma_steps=Gs_) for y1 in y_st])*24.0   # downwash at the stations, in units of Γ1/2s
    return {"staircase Γ(y)/Γ1": (np.repeat(ye_, 2), np.r_[0.0, np.repeat(Gs_, 2), 0.0]),   # the steps, drawn with vertical risers
            "smooth limit: elliptic Γ(y)/Γ1": (y_st, ell(y_st)),                            # what the staircase tends to
            "downwash of the staircase / (Γ1/2s)": (y_st, np.clip(w_, -1.0, 4.0)),           # cut at −1 and 4: it is unbounded next to a filament
            "smooth limit: uniform downwash": (y_st, np.ones_like(y_st))}                    # the elliptic loading's Γ1/2s
figF6 = slider_figure(staircase, "number of horseshoes N", [1, 3, 5, 7, 11, 21, 41], unit="", xlabel="spanwise position y [m]",
                      ylabel="Γ/Γ1 and w/(Γ1/2s) [–]", yrange=[-1.05, 4.05], active=3,
                      title="A wing built from nested horseshoes: loading and downwash settle onto smooth curves")   # 7 slider positions
figF6.show()                                                 # display
""", explain=r"""
1. `staircase(N)` builds N nested horseshoes that approximate the elliptic loading and evaluates their downwash at 60
   stations with `LL.staircase_downwash`.
2. Four traces: the staircase and its smooth limit; the staircase's downwash and the uniform value the limit has.""")
nb.figure_notes(
    see="A slider for the number of horseshoes. The staircase loading closes on the half-ellipse; its downwash, a jagged "
        "curve with a spike next to every trailing leg, closes on the flat line at 1.",
    read="With one horseshoe the downwash is 0.64 at mid-span and unbounded at the tips. Every added horseshoe adds a pair of "
         "weaker spikes but pulls the level between them toward 1. The principal-value integral is what remains when the "
         "spikes have become infinitely many and infinitely weak.",
    change="…the loading were not elliptic? The limit curve of the downwash would no longer be flat (C11, C12).")
note("N41 [B]", r"""
Along the flight direction: $w(x\to-\infty)=0$, $w(0)=w$, $w(x\to+\infty)=2w$ — far behind the wing every filament looks infinite in both
directions, so the downwash there is twice its value at the wing.""")
nb.figure(r"""
xx = np.linspace(-5.0, 5.0, 401)                             # distance along the flight direction in spans (x > 0 is behind the wing)
wx = LL.downwash_along_x(x=xx, Gn=np.array([1.0]), s=1.0)    # downwash on the centre line for an elliptic loading with Γ1 = s = 1 (so w at the wing = 0.5)
fig, ax = plt.subplots(figsize=(6.2, 3.2))                   # one panel
ax.plot(xx, wx, color=COLORS["blue"])                        # the downwash
for lvl in (0.5, 1.0): ax.axhline(lvl, color="0.7", ls="--", lw=1)   # the value at the wing and twice it
ax.set(xlabel="distance behind the wing $x/s$ [–]", ylabel="downwash $w$ [in units of Γ1/s]", title="Half the final downwash is already there at the wing")   # labels
plt.show()                                                   # display
print("w far ahead, at the wing, far behind:", np.round(LL.downwash_along_x(x=np.array([-50.0, 0.0, 50.0]), Gn=np.array([1.0]), s=1.0), 3))   # 0.000, 0.500, 1.000
""",
    see="The downwash on the centre line from five spans ahead of the wing to five spans behind it.",
    read="Far ahead the air is undisturbed. At the wing the downwash is 0.5 (in these units); far behind it is 1.0 — exactly "
         "twice. The wing meets air that has received half of the downward push it will eventually get.",
    change="…you flew a second wing two spans behind the first? It would sit in nearly the full 2w — a tailplane works in "
           "strongly sinking air.")
note("N65 [B]", r"""
Outboard of the tips the trailing vortices push the air **up**. A bird flying just outboard and behind a neighbour's wing tip
flies in rising air — the reason for formation flight.""")
slip(10, "the “Λ-pattern” of migrating birds", "a V formation — the letter there is a shape, not the aspect ratio $\\Lambda$.")
nb.figure(r"""
Xg, Yg = np.meshgrid(np.linspace(-6.0, 24.0, 121), np.linspace(-14.0, 14.0, 141))   # a patch of the wing's plane: x downstream, y along the span [m]
wmap = LL.horseshoe_downwash(x=Xg, y=Yg, Gamma=10.0, s=12.0) # downwash (positive down) of one horseshoe of span 12 m, Γ = 10 m²/s
fig, ax = plt.subplots(figsize=(6.6, 3.8))                   # one panel
pc = ax.pcolormesh(Xg, Yg, np.clip(wmap, -0.6, 0.6), cmap="RdBu", shading="auto")   # blue = air going down, red = air going up (cut at ±0.6 m/s)
ax.plot([0, 0], [-6, 6], color="k", lw=3)                    # the bound vortex (the wing)
for yt in (-6, 6): ax.plot([0, 24], [yt, yt], color="k", lw=1)   # the two trailing legs
fig.colorbar(pc, ax=ax, label="downwash $w$ [m/s] (positive down)")   # colour scale
ax.set(xlabel="distance behind the wing $x$ [m]", ylabel="spanwise position $y$ [m]", title="Down between the tips, up outboard of them")   # labels
plt.show()                                                   # display
print(f"beside the wing, 2 m outboard of the tip (x = 0, y = 8 m): w = {float(LL.horseshoe_downwash(x=0.0, y=8.0, Gamma=10.0, s=12.0)):.3f} m/s")   # −0.341: an upwash
""",
    see="A map of the vertical air motion in the plane of a horseshoe vortex, seen from above: the wing is the thick line on "
        "the left, the two legs trail to the right.",
    read="Blue (sinking air) fills the lane between the legs; red (rising air) lies outside them. The printed value, −0.341 "
         "m/s two metres outboard of the tip, is an *upwash*: free lift for whoever flies there.",
    change="…the neighbour flew directly behind instead? It would sit in the blue lane, in air sinking at up to twice the "
           "leader's own downwash.")
note("N66 [B]", r"""
**Ground effect.** A wall can be replaced by a mirror vortex of opposite strength (the method of images, Ch. 5 §5.7 and Ch. 6).
An image horseshoe below the ground reduces the downwash at the wing's centre to
$w=\dfrac{\Gamma}{\pi s}\Big[1-\dfrac1{1+(4h/s)^2}\Big]$ (ours), $h$ being the height above the ground.""")
nb.code(r"""
for h in (1.2, 3.0, 6.0):                                    # heights above the ground [m] for a 12 m span: s/10, s/4, s/2
    own = LL.horseshoe_downwash(x=0.0, y=0.0, Gamma=10.0, s=12.0)               # downwash of the wing's own horseshoe at its centre [m/s]
    image = LL.horseshoe_downwash(x=0.0, y=0.0, Gamma=-10.0, s=12.0, z=2*h)     # the mirror horseshoe (opposite Γ) lies 2h below: we are 2h above its plane
    print(f"h = {h:3.1f} m (h/s = {h/12:.2f}): factor = {LL.ground_effect_factor(h=h, s=12.0):.3f}; from the two horseshoes: {(own + image)/own:.3f}")   # 0.138, 0.5, 0.8 — twice each
""", explain=r"""
**What does this show?** The closed form and the sum of two horseshoes agree. At a height of half a span the downwash is 80 %
of its free-air value, at a quarter span 50 %, at a tenth of a span 14 %: less downwash means less induced drag (C10) — a
wing just above the runway *floats*.""")
note("N61 [C]", r"""
A parabolic loading $\Gamma=\Gamma_0\big(1-(2y/s)^2\big)$: the next cell compares its mid-span downwash with that of the elliptic loading of the same
$\Gamma_0$; any $\Gamma(y)$ is handled by the series of C11.""")
nb.code(r"""
w_par = LL.downwash(y1=0.0, Gamma_fn=lambda y: 1 - (2*y/12.0)**2, s=12.0)   # mid-span downwash of a parabolic loading with Γ0 = 1 m²/s [m/s]
print(f"parabolic / elliptic mid-span downwash = {w_par/LL.downwash(y1=0.0, Gamma_fn=ell, s=12.0):.3f}")   # 1.273 = 4/π
""")
whatif(r"""…a wing section flies in air that is sinking at $w$? Its lift tilts. C10.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C10", r"Induced drag: $L=\displaystyle\int_{-s/2}^{s/2}\rho U\Gamma\,dy$, $D_i=\displaystyle\int_{-s/2}^{s/2}\rho w\Gamma\,dy$ (14.15)",
     "How can a wing have drag in a fluid with no friction?")
problem(r"""
Walk up a down escalator: you work just to stay level. A finite wing flies through air that its own wake is pushing downward.
Relative to that air the wing is climbing a gentle slope, and the force on it, which is perpendicular to the air it actually
meets, leans back by the slope angle. The backward part is a drag — with no friction anywhere. It is the price of lift on a
wing with tips.""")
idea(r"""
        U  (flight path)
     ─────────────────────►•
          ε  ╲             │ w  (downwash)
               ╲           ▼
        U_e      ╲►                    the air the strip actually meets is tilted down by ε

        dL ▲   ◄── dL_e is perpendicular to U_e, so it leans BACK by ε
           │  ╱
           │╱____► dD_i               its backward part is the induced drag
""")
P("P350", "relative velocity as a vector subtraction", r"""
The wind a moving body feels is the air's velocity minus the body's own: $\mathbf U_r=\mathbf U-\mathbf V$. Its direction, not the
direction of the true wind, decides where "lift" and "drag" point. (`np.hypot(x, y)` is the length of the vector $(x, y)$;
`np.arctan2(y, x)` its angle in the right quadrant — primed in Ch. 2, P70.)""", r"""
Uv = np.array([50.0, 0.0]); wv = np.array([0.0, -1.0])           # the stream [m/s], and a 1 m/s downwash
Ue = Uv + wv                                                     # the air the strip meets
print(np.hypot(*Ue), np.degrees(np.arctan2(1.0, 50.0)))          # 50.01 m/s, tilted 1.15 degrees
""")
remind("C10")
nb.code(r"""
fig = draw_force_triangle()                                  # the stream U, the downwash w, the tilted stream U_e, and the force normal to it split into dL and dD_i
plt.show()                                                   # display the sketch the derivation below refers to
""", tags=["figure"])
D("D13", "Effective angle, the strip's forces, and lift and induced drag of the wing", ref="14.15")
slip(2, r"""$\varepsilon=\tan\dfrac wU$""", r"""$\varepsilon=\tan^{-1}\dfrac wU$ (the inverse tangent; for $w\ll U$ both are $\simeq w/U$).""")
nb.worked_example("one strip by hand", r"""
$U$ = 50 m/s, downwash $w$ = 1 m/s, a strip of width 1 m carrying Γ = 20 m²/s, ρ = 1.2 kg/m³, geometric angle α = 5°.

1. Tilt: ε = $w/U$ = 0.02 rad = 1.15°.
2. Effective angle $\alpha_e=\alpha-\dfrac wU$ (14.14) = 5° − 1.15° = 3.85°.
3. Lift of the strip: $\rho U\Gamma\,dy$ = 1.2 × 50 × 20 × 1 = 1200 N.
4. Induced drag: $\rho w\Gamma\,dy$ = 1.2 × 1 × 20 × 1 = 24 N.
5. Ratio 24/1200 = 0.02 = ε: the drag is the lift times the tilt.""")
nb.code(r"""
sf_ = LL.strip_forces(Gamma=10.0, w=1.0, U=30.0, rho=1.225, dy=0.1)               # a strip 0.1 m wide carrying 10 m²/s in a 1 m/s downwash: small-angle forms
sx_ = LL.strip_forces(Gamma=10.0, w=1.0, U=30.0, rho=1.225, dy=0.1, exact=True)   # the same without the small-angle step
print(f"strip: dL = {sf_['dL']:.2f} N, dD_i = {sf_['dDi']:.3f} N, ε = {sf_['eps']:.4f} rad (exact: dL = {sx_['dL']:.2f}, dD_i = {sx_['dDi']:.3f})")   # 36.75, 1.225, 0.0333
W, s, rho, U = ac["W"], ac["s"], 1.225, ac["U_slow"]         # our aircraft at 30 m/s: weight [N], span [m], density [kg/m³], speed [m/s]
G1 = 4*W/(np.pi*s*rho*U)                                     # the elliptic loading that carries the weight: Γ1 from L = (πs/4)ρUΓ1 (derived in C13)
w_el = LL.elliptic_downwash(Gamma1=G1, s=s)                  # Eq. (14.24): its uniform downwash Γ1/2s [m/s]
eps_el = LL.induced_angle(w=w_el, U=U)                       # the tilt ε ≈ w/U [rad]
a_eff = LL.effective_angle(alpha=np.deg2rad(5.0), w=w_el, U=U)   # Eq. (14.14): the angle a strip set at 5° really works at [rad]
yy = np.linspace(-s/2, s/2, 4001)                            # stations along the span [m]
L_w, Di_w = LL.lift_and_induced_drag(y=yy, Gamma=LL.elliptic_circulation(y=yy, s=s, Gamma1=G1), w=w_el + 0*yy, rho=rho, U=U)   # Eq. (14.15): both integrals by the trapezoid rule
print(f"Γ1 = {G1:.2f} m²/s, w = {w_el:.3f} m/s, ε = {eps_el:.4f} rad = {np.degrees(eps_el):.2f}°, effective angle = {np.degrees(a_eff):.2f}° ({a_eff:.4f} rad)")   # 25.98, 1.083, 2.07°, 2.93°
print(f"L = {L_w:.0f} N, D_i = {Di_w:.1f} N = {100*Di_w/L_w:.1f} % of the lift")   # 9000 N, 324.8 N, 3.6 %
""", explain=r"""
1. `LL.strip_forces` returns the lift and induced drag of one strip; `exact=True` keeps $U_e$ and ε exactly. The two differ
   in the fourth digit: the small-angle forms are all we need.
2. For our aircraft at 30 m/s an elliptic loading needs Γ₁ = 25.98 m²/s (clockwise); its downwash is 1.083 m/s everywhere,
   a tilt of 2.07°.
3. `LL.effective_angle` is $\alpha_e=\alpha-\dfrac wU$ (14.14): a strip set at 5° works at 2.93°.
4. `LL.lift_and_induced_drag` evaluates $L=\int\rho U\Gamma\,dy$, $D_i=\int\rho w\Gamma\,dy$ (14.15): 9000 N of lift cost 324.8 N of induced drag,
   3.6 % of the lift — with no friction anywhere.""")
scratch(r"""
G_, w_, U_, rho_, dy_ = 10.0, 1.0, 30.0, 1.225, 0.1          # circulation [m²/s], downwash [m/s], speed [m/s], density [kg/m³], strip width [m]
Ue = np.hypot(U_, w_)                                        # speed of the tilted stream [m/s]
eps = np.arctan2(w_, U_)                                     # its tilt [rad], exactly
dLe = rho_*Ue*G_*dy_                                         # lift theorem in the tilted stream [N], perpendicular to U_e
assert np.isclose(dLe*np.sin(eps), rho_*w_*G_*dy_)           # backward part = ρ w Γ dy exactly, because U_e sin ε = w
assert np.isclose(dLe*np.cos(eps), rho_*U_*G_*dy_)           # upward part = ρ U Γ dy exactly, because U_e cos ε = U
assert np.allclose([dLe*np.cos(eps), dLe*np.sin(eps)], [sx_["dL"], sx_["dDi"]])   # and the tested function returns these two numbers
print(f"dL_e = {dLe:.3f} N at ε = {np.degrees(eps):.3f}° → dL = {dLe*np.cos(eps):.3f} N, dD_i = {dLe*np.sin(eps):.4f} N")   # 36.770 N at 1.909° → 36.750, 1.2250
""", r"""
Both resolutions are **exact**: $U_e\sin\varepsilon=w$ and $U_e\cos\varepsilon=U$. Only calling them "lift" and "drag" relative to the flight
path, with a section law at the reduced angle, is approximate.""")
nb.figure(r"""
rect = LL.solve_lifting_line(c=1.25, s=12.0, U=30.0, alpha=0.1, K=2*np.pi, N=32, rho=1.225)   # a rectangular wing (chord 1.25 m, span 12 m) at α = 0.1 rad = 5.7°; the solver is taught in C12
fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.0))          # one force triangle per station
for ax, y1 in zip(axes, (0.0, 3.0, 5.5)):                    # mid-span, half-way out, near the tip [m]
    w1 = np.interp(y1, rect["y"], rect["w"])                 # downwash at this station [m/s]
    e5 = 5*w1/30.0                                           # the tilt ε = w/U, exaggerated five times for the drawing [rad]
    ax.annotate("", xy=(1, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["blue"]))              # the stream U
    ax.annotate("", xy=(1, -np.tan(e5)), xytext=(1, 0), arrowprops=dict(arrowstyle="->", color=COLORS["blue"]))    # the downwash w (×5)
    ax.annotate("", xy=(0.5 + 0.8*np.sin(e5), 0.8*np.cos(e5)), xytext=(0.5, 0), arrowprops=dict(arrowstyle="->", color="k"))   # the force, normal to the tilted stream
    ax.annotate("", xy=(0.5, 0.8*np.cos(e5)), xytext=(0.5, 0), arrowprops=dict(arrowstyle="->", color=COLORS["teal"]))         # its lift part
    ax.annotate("", xy=(0.5 + 0.8*np.sin(e5), 0.8*np.cos(e5)), xytext=(0.5, 0.8*np.cos(e5)), arrowprops=dict(arrowstyle="->", color=COLORS["rose"]))   # its induced-drag part
    ax.set(xlim=(-0.05, 1.15), ylim=(-0.45, 0.95), title=f"y = {y1} m: w = {w1:.2f} m/s, ε = {np.degrees(w1/30.0):.2f}° (drawn ×5)")   # readout
    ax.set_aspect("equal"); ax.axis("off")                   # true proportions, no axes
plt.show()                                                   # display
""",
    see="The force triangle at three stations of a rectangular wing: stream and downwash (blue), the force on the strip "
        "(black), its lift part (teal) and its induced-drag part (rose). The tilt is drawn five times too large.",
    read="On a rectangular wing the downwash is not uniform: it is larger near the tip than at mid-span, so the force leans "
         "back more there. The rose arrow is the induced drag of that strip.",
    change="…the loading were elliptic? All three triangles would be identical: the downwash would be the same everywhere (C14).")
note("N45 [B]", r"""
Read $L=\int\rho U\Gamma\,dy$ and $D_i=\int\rho w\Gamma\,dy$ (14.15) as one rule used twice: $U$ with Γ gives a force normal to $U$; $w$ with Γ gives a force
normal to $w$. With an infinite span there are no trailing vortices, $w=0$ and the drag vanishes. (The book's remark that this
is the largest part of a wing's drag is tested in C15, N75.)""")
note("N38 [B]", r"""
A second reading: the tip vortices get longer every second, so kinetic energy is fed into the wake at the rate
$D_iU=\dfrac{d}{dt}(\text{wake kinetic energy})$ (stated, not derived). Power is force times speed (primed in Ch. 4, P127), so the induced drag
equals the kinetic energy left behind *per metre flown*.""")
nb.code(r"""
Gn_r = rect["Gn"]                                            # the sine coefficients of the rectangular wing's loading [m²/s] (C11–C12 explain them)
print(f"wake kinetic energy per metre = {LL.wake_kinetic_energy_per_length(Gn=Gn_r, s=12.0, rho=1.225):.3f} J/m; induced drag = {LL.induced_drag_from_series(Gn=Gn_r, rho=1.225):.3f} N")   # equal: 1 J/m = 1 N
print(f"one unit coefficient, ρ = 1: {LL.wake_kinetic_energy_per_length(Gn=np.array([1.0]), s=1.0, rho=1.0):.4f} and {LL.induced_drag_from_series(Gn=np.array([1.0]), rho=1.0):.4f} (π/8 = {np.pi/8:.4f})")   # 0.3927 three times
""")
note("N67 [C]", r"""
A tip vortex with a larger core holds less energy per metre, roughly $\propto\Gamma^2\big[\ln(s/\sigma)+\tfrac14\big]$ (ours) for a pair of Rankine
vortices (solid-body core of radius σ, Ch. 3) a distance $s$ apart — one line:""")
nb.code(r"""
print([round(float(ch14.vortex_pair_energy_per_length(Gamma=14.17, s=9.42, sigma=sig, rho=1.225)), 1) for sig in (0.1, 0.3, 1.0)])   # energy per metre [J/m] for core radii 0.1, 0.3, 1 m: it falls as the core grows
""")
explainer("downwash_induced_drag", "How can there be drag without friction?",
          "building the wing from 1, 3, 7, 21 nested horseshoes and watching Γ(y) and w(y) settle onto smooth curves shows "
          "what the integral is; moving the station turns the local force triangle.",
          ["Start with one horseshoe and step N up.",
           "Slide the station to the tip: what happens to the bars of the nearest filaments?",
           "Move outboard of the tip: which way does the air move?",
           "Lower the wing toward the ground and watch the downwash."])
whatif(r"""…Γ(y) is not given but has to be *found*? The downwash depends on Γ, and Γ depends on the downwash. We need a way to
do the integral once and for all: C11.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C11", r"The sine series and its downwash: $\Gamma=\displaystyle\sum_{n=1}^{\infty}\Gamma_n\sin(n\gamma)$ (14.18) ⇒ $w(y_1)=\dfrac{1}{2s}\displaystyle\sum_{n=1}^{\infty}n\,\Gamma_n\,\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19, last member)",
     "The downwash integral is singular and awkward. Can it be done once, for every loading?")
problem(r"""
A violin string is fixed at both ends, and every shape it can take is a sum of sine waves that vanish at both ends. The
circulation along a wing is "fixed" the same way: zero at both tips. Write it as a sum of sines of an angle that runs from 0
at one tip to π at the other, and the awkward integral for the downwash can be done for each sine separately — the answer is
another sine, divided by sin γ.""")
nb.md(r"""
#### The idea

Put a semicircle over the span. A point at the angle γ on it projects straight down to the spanwise position
$y=-(s/2)\cos\gamma$. Equal steps in γ crowd toward the tips — exactly where Γ changes fastest.""")
nb.figure(r"""
fig = draw_span_semicircle()                                 # the span angle γ: a point of the span is the foot of a point on a semicircle of radius s/2
plt.show()                                                   # display
""",
    see="A semicircle of radius $s/2$ over the span; a point at the angle γ on it and its foot on the span.",
    read="γ = 0 is the port tip, γ = 90° mid-span, γ = 180° the starboard tip. Equal steps in γ give short steps in $y$ near "
         "the tips and long ones in the middle.",
    change="…you used $y$ itself? The downwash integral would keep its awkward $1/(y_1-y)$; in γ it becomes a known integral.")
note("N49 [B]", r"""
The span angle: $y=-(s/2)\cos\gamma$, $y_1=-(s/2)\cos\gamma_1$: γ = 0 at the port tip, π/2 at mid-span, π at the starboard tip.

> ⚠️ **Common confusion (trap T3):** this γ is an angle — not Chapter 5's sheet strength, not the ratio of specific heats.""")
nb.code(r"""
print(np.round(LL.span_angle(y=np.array([-6.0, 0.0, 3.0, 6.0]), s=12.0), 3))   # port tip, mid-span, half-way out, starboard tip of a 12 m span: 0, 1.571, 2.094, 3.142 rad (0°, 90°, 120°, 180°)
""")
P("P346", "half-range sine series", r"""
On 0 ≤ γ ≤ π any function that vanishes at both ends can be written as $\sum a_n\sin n\gamma$, and the sines are orthogonal there:
$\int_0^\pi\sin n\gamma\,\sin m\gamma\,d\gamma=\tfrac\pi2$ if $n=m$ and 0 otherwise. (Full-period integrals of sines: Ch. 6, P151; keeping only the
first N terms of such a series: Ch. 11, P275.)""", r"""
g = np.linspace(0, np.pi, 2001)                       # the half range 0…π
print(np.trapezoid(np.sin(2*g)*np.sin(3*g), g))       # 0: different sines
print(np.trapezoid(np.sin(3*g)**2, g), np.pi/2)       # 1.5708: the same sine
""")
remind("C11")
note("N50 [B]", r"""
The loading as a sine series: $\Gamma=\displaystyle\sum_{n=1}^{\infty}\Gamma_n\sin(n\gamma)$ (14.18). Every term vanishes at both tips; the numbers $\Gamma_n$ [m²/s] are
the unknowns of the theory. Odd $n$ are symmetric about mid-span (even functions of $y$), even $n$ antisymmetric (odd functions).""")
nb.figure(r"""
gg = np.linspace(0, np.pi, 401)                              # span angle from port tip to starboard tip
yy1 = LL.span_position(gamma=gg, s=1.0)                      # the same as a position on a span of 1: y/s from −0.5 to 0.5
fig, ax = plt.subplots(figsize=(6.4, 3.4))                   # one panel
for n_, col in zip((1, 2, 3, 4), (COLORS["teal"], COLORS["rose"], COLORS["blue"], COLORS["amber"])):   # the first four shapes
    ax.plot(yy1, np.sin(n_*gg), color=col, label=f"sin({n_}γ)")   # one sine shape against the spanwise position
ax.set(xlabel="spanwise position $y/s$ [–]", ylabel="shape of the loading [–]", title="The building blocks of any loading")   # labels
ax.legend(fontsize=8); plt.show()                            # legend and display
""",
    see="The first four sine shapes plotted against the position along the span.",
    read="$n$ = 1 is a half-ellipse (a circle squashed, since $\\sin\\gamma=\\sqrt{1-(2y/s)^2}$). $n$ = 2 is antisymmetric: more "
         "lift on one wing than the other — a rolling wing. $n$ = 3 is symmetric with a dip in the middle, $n$ = 4 "
         "antisymmetric again. Every one is zero at both tips.",
    change="…a wing is symmetric and flies straight? Only the odd shapes appear (C12 confirms: the even coefficients vanish).")
D("D14", "The sine series inside the downwash integral", ref="14.13")
slip(1, r"""the first two members of its three-member downwash equation (14.19) with the denominator $\cos\gamma_1-\cos\gamma$,
$w(y_1)=\dfrac1{2\pi s}\displaystyle\int_0^\pi\sum n\Gamma_n\dfrac{\cos(n\gamma)\,d\gamma}{\cos\gamma_1-\cos\gamma}$""",
     r"""$\cos\gamma-\cos\gamma_1$ (step 5 above). The last member, $w(y_1)=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$, and everything the book builds on it are right.""")
nb.code(r"""
w_last = LL.downwash_series(Gn=np.array([1.0]), gamma1=np.pi/2, s=1.0)                         # last member of (14.19) for the elliptic loading (Γ1 = s = 1), at mid-span
w_printed = LL.downwash_series(Gn=np.array([1.0]), gamma1=np.pi/2, s=1.0, printed_middle=True)  # the middle member evaluated as the book prints it
w_quad = LL.downwash(y1=0.0, Gamma_fn=lambda y: np.sqrt(1 - 4*y**2), s=1.0)                    # direct quadrature of the downwash integral (14.13)
print(f"last member: {w_last:+.3f}   printed middle member: {w_printed:+.3f}   direct integral: {w_quad:+.3f}   [m/s]")   # +0.5, −0.5, +0.5
""", explain=r"""
**What does this show?** The corrected series and the direct integral agree at +0.5 (downward). The printed form gives −0.5:
it would make the air *rise* behind a lifting wing.""")
P("P347", "product-to-sum identity and a three-term recursion", r"""
$\cos(n+1)\gamma+\cos(n-1)\gamma=2\cos\gamma\,\cos n\gamma$. A rule that gives each term of a sequence from the two before it,
$I_{n+1}=2cI_n-I_{n-1}$, fixes the whole sequence once two starting values are known (linear recurrences: Ch. 10, P236; proof
by induction — true for the first cases, and true for the next whenever true for the previous ones: Ch. 12, P306).""", r"""
cg, I = np.cos(1.0), [0.0, np.pi]                     # cg = cos(gamma_1) with gamma_1 = 1 rad; the starting values I_0 and I_1
for n in range(1, 4): I.append(2*cg*I[n] - I[n-1])    # the recursion
print(I[3], np.pi*np.sin(3*1.0)/np.sin(1.0))          # both 0.5269
""")
D("D15", "Glauert's integral, and the downwash of the series", ref="14.19")
nb.worked_example("two sines by hand", r"""
Γ₁ = 12 and Γ₃ = 1 m²/s on a span $s$ = 12 m.

1. Mid-span, γ₁ = 90°: sin γ₁ = 1, sin 3γ₁ = −1, so Γ = 12 − 1 = 11 m²/s.
2. Downwash there from $w(y_1)=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19): (1/24)[1 × 12 × 1 + 3 × 1 × (−1)] = 9/24 = 0.375 m/s.
3. Near the port tip, γ₁ = 30° ($y_1$ = −5.2 m): sin γ₁ = 0.5, sin 3γ₁ = 1, Γ = 6 + 1 = 7 m²/s.
4. Downwash there: (1/24)[12 × 0.5/0.5 + 3 × 1 × 1/0.5] = (12 + 6)/24 = 0.75 m/s.
5. The $n$ = 1 term alone gives 0.5 m/s everywhere; the small Γ₃ lowers the downwash in the middle and raises it by half
   toward the tips.""")
nb.code(r"""
Gn2 = np.array([12.0, 0.0, 1.0])                             # the coefficients of the hand example: Γ1 = 12, Γ2 = 0, Γ3 = 1 [m²/s]
for g1 in (np.pi/2, np.pi/6):                                # mid-span (90°) and near the port tip (30°)
    print(f"γ1 = {np.degrees(g1):.0f}° ({g1:.3f} rad): Γ = {LL.circulation_series(Gn=Gn2, gamma=g1):.1f} m²/s, w = {LL.downwash_series(Gn=Gn2, gamma1=g1, s=12.0):.3f} m/s")   # 11 and 0.375; 7 and 0.75
print("Glauert's integral at γ1 = 60°, n = 0…4, by quadrature:", np.round([LL.glauert_integral(n=n, gamma1=np.pi/3, method="quad") for n in range(5)], 6))     # 0, π, π, 0, −π
print("                               the closed form:        ", np.round([LL.glauert_integral(n=n, gamma1=np.pi/3, method="closed") for n in range(5)], 6))   # the same
""", explain=r"""
1. `LL.circulation_series` is $\Gamma=\sum\Gamma_n\sin(n\gamma)$ (14.18) and `LL.downwash_series` its downwash
   $w(y_1)=\dfrac1{2s}\sum n\Gamma_n\dfrac{\sin(n\gamma_1)}{\sin\gamma_1}$ (14.19): the hand example's 11 m²/s and 0.375 m/s at mid-span, 7 m²/s and 0.75 m/s near the tip.
2. `LL.glauert_integral` evaluates the integral of the derivation numerically (`"quad"`, the singular part subtracted) and
   by the closed form $\pi\sin n\gamma_1/\sin\gamma_1$: 0, π, π, 0, −π for $n$ = 0 … 4 at γ₁ = 60°.
3. Both routes agree to the digits shown.""")
scratch(r"""
gam_of = lambda y: np.arccos(-2*y/12.0)                      # span angle of a station y on a 12 m span
Gfn = lambda y: 12*np.sin(gam_of(y)) + np.sin(3*gam_of(y))   # the two-sine loading as a plain function of y [m²/s]
y_test = np.array([-5.0, -2.0, 0.0, 3.0, 5.5])               # five stations [m]
w_quad5 = np.array([LL.downwash(y1=y1, Gamma_fn=Gfn, s=12.0) for y1 in y_test])   # direct principal-value integral, Eq. (14.13)
w_ser5 = LL.downwash_series(Gn=Gn2, gamma1=LL.span_angle(y=y_test, s=12.0), s=12.0)   # the series, Eq. (14.19)
assert np.allclose(w_quad5, w_ser5, rtol=1e-6)               # same numbers → the series is the integral, done once
print(np.round(w_quad5, 4))                                  # the downwash at the five stations [m/s]
""", r"""
The loading of the hand example is handed to the quadrature of the original integral as a plain function of $y$; the result
equals the sine-series formula at all five stations.""")
nb.figure(r"""
g_in = np.linspace(0.02, np.pi - 0.02, 400)                  # span angles, staying clear of the tips
fig, (axm, axb) = plt.subplots(1, 2, figsize=(9.6, 3.4))     # left: downwash of single modes; right: the sign check
for n_, col in zip((1, 2, 3), (COLORS["teal"], COLORS["rose"], COLORS["blue"])):   # three single-sine loadings
    Gn_one = np.zeros(n_); Gn_one[-1] = 1.0                  # only Γ_n = 1 m²/s
    axm.plot(LL.span_position(gamma=g_in, s=1.0), LL.downwash_series(Gn=Gn_one, gamma1=g_in, s=1.0), color=col, label=f"n = {n_}")   # its downwash along the span (s = 1 m)
axm.set(xlabel="spanwise position $y/s$ [–]", ylabel="downwash $w$ [m/s]", ylim=(-5, 5), title="Only the first sine has a uniform downwash")   # labels
axm.legend(fontsize=8)                                       # which curve is which
axb.bar(["last member\n(correct)", "direct\nintegral", "middle member\nas printed"], [w_last, w_quad, w_printed], color=[COLORS["teal"], COLORS["blue"], COLORS["rose"]])   # the three numbers of the slip check
axb.axhline(0, color="k", lw=0.8)                            # zero line
axb.set(ylabel="mid-span downwash [m/s]", title="slip #1: the printed sign is wrong")   # labels
plt.show()                                                   # display
""",
    see="Left: the downwash produced by the loadings sin γ, sin 2γ and sin 3γ alone. Right: three bars — the corrected "
        "series, the direct integral, and the middle member as the book prints it.",
    read="The $n$ = 1 curve is flat: an elliptic loading sinks the air evenly. The others swing between down and up and grow "
         "toward the tips (the tip values are ±$n^2/2s$). On the right, two bars point down (+0.5) and the printed form "
         "points up.",
    change="…a little sin 3γ is added to the ellipse? The flat line gains a wiggle: less downwash in the middle, more at the tips.")
whatif(r"""…we do not know the $\Gamma_n$? Each strip of the wing tells us one relation between its circulation and its downwash. C12.""")


# ---------------------------------------------------------------------------------------------------------------------
core("C12", r"The lifting-line equation: $\dfrac K2\,Uc(\alpha+\beta)=\displaystyle\sum_{n=1}^{\infty}\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ (14.20)",
     "Given a wing's shape and angle, how do you find the circulation along it?")
problem(r"""
A designer draws a wing — this span, this taper, this twist — and wants its lift and its drag. Each strip of the wing behaves
like the two-dimensional section of C06, but at the *reduced* angle the downwash leaves it. The downwash comes from the
circulation of all the strips, and the circulation of each strip depends on its downwash: a loop. Written with the sine
series the loop is a small set of linear equations, and a computer closes it in a millisecond.""")
idea(r"""
   chord c(y), angle α(y) ──►  [ each strip:  Γ = (K/2) U c (α − w/U + β) ]  ◄──►  [ w from all the Γ_n ]
                                                   │
                                                   ▼
                          N unknowns Γ_n ,  N stations γ_j   ──►  np.linalg.solve  ──►  Γ(y), w(y), L, D_i
""")
P("P345", "integral equation", r"""
An equation in which the unknown *function* appears under an integral sign. After discretising, it is a matrix equation with
the unknown on both sides: $x=r-Mx$, i.e. $(I+M)x=r$. (`np.linalg.solve(A, r)` solves $Ax=r$ — primed in Ch. 1, P57; `np.eye(n)`
is the identity matrix and `M @ x` a matrix–vector product.)""", r"""
M = np.array([[0.2, 0.1], [0.1, 0.2]]); r2 = np.array([1.0, 1.0])   # x appears on both sides of x = r2 - M @ x
x2 = np.linalg.solve(np.eye(2) + M, r2)                             # move it to one side and solve
print(x2, r2 - M @ x2)                                              # the same vector twice: [0.769 0.769]
""")
remind("C12")
D("D16", "The lifting-line equation: from the section law to a linear system", ref="14.20")
confusion(r"""(traps T12 and T3) $K$ is the section's lift-curve slope **per radian**; `K=` is a required keyword of every `LL`
function, so that 2π is never assumed silently. And β in $C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.16) is *minus the zero-lift angle* of the
section the wing has — for a Zhukhovsky section it coincides with the camber angle of C06.""")
note("N52 [B]", r"""
**The procedure, written out for N = 3.** Keep three terms and demand the equation at the three stations γ = 45°, 90°, 135°
($\gamma_j=j\pi/(N+1)$ — never a tip, where sin γ = 0; trap T10): $\mathsf A_{jn}\Gamma_n=r_j$ with
$\mathsf A_{jn}=\Big(1+\dfrac{nKc_j}{4s\sin\gamma_j}\Big)\sin n\gamma_j$ and $r_j=\dfrac K2Uc_j(\alpha_j+\beta_j)$ (ours, formalising the book's description).""")
nb.code(r"""
A3, g3 = LL.lifting_line_matrix(c=wd["c"], s=wd["s"], K=wd["K"], N=3)                           # the 3 × 3 matrix and its three stations for the demo wing (chord 1 m, span 6 m)
r3 = LL.lifting_line_rhs(c=wd["c"], s=wd["s"], U=wd["U"], alpha=wd["alpha"], K=wd["K"], N=3)    # the right-hand side (K/2) U c α at each station [m²/s]
print("stations γ_j [°]:", np.degrees(g3).round(1), "=", g3.round(3), "rad")                    # 45°, 90°, 135°
print("matrix A (rows = stations, columns = sines):\n", A3.round(3))                            # [0.969 1.740 1.493], [1.262 0 −1.785], [0.969 −1.740 1.493]
print("right-hand side r:", r3.round(4), " → Γ_n =", np.linalg.solve(A3, r3).round(4), "m²/s")  # 0.2742 three times → [0.2487, 0, 0.0222]
""", explain=r"""
**What does this show?** Three stations, three unknowns. The matrix entries are the bracket of
$\dfrac K2\,Uc(\alpha+\beta)=\sum\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ (14.20) times $\sin n\gamma_j$; the right-hand side is the same 0.2742 m²/s at
every station of an untwisted rectangle. The solution has Γ₂ = 0 (a symmetric wing) and a small positive Γ₃.""")
nb.worked_example("one term, one station, by hand", r"""
A rectangular wing with $s$ = 6 m, $c$ = 1 m (Λ = 6), $K$ = 2π, $U$ = 1 m/s, α = 5° = 0.0873 rad, β = 0. Keep only Γ₁ and use the
mid-span station γ = 90° (sin γ = 1).

1. The bracket of $\dfrac K2\,Uc(\alpha+\beta)=\sum\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ (14.20): 1 + 2π × 1/(4 × 6) = 1.262.
2. Left side: π × 1 × 1 × 0.0873 = 0.2742.
3. Γ₁ = 0.2742/1.262 = 0.2173 m²/s.
4. Lift coefficient: $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21) divided by ½ρU²sc gives πΓ₁/(2Uc) = 0.341.
5. The two-dimensional section alone would give 2π × 0.0873 = 0.548: the finite wing delivers 62 % of it in this crude
   estimate (72 % when enough terms are kept — next cell).""")
nb.code(r"""
sol = LL.solve_lifting_line(c=LL.planform(kind="rectangular", s=wd["s"], c_root=wd["c"]), s=wd["s"], U=wd["U"], alpha=wd["alpha"], K=wd["K"], N=32, rho=wd["rho"])   # the demo wing with 32 sine terms
print(f"C_L = {sol['CL']:.4f}, slope dC_L/dα = {sol['slope']:.3f} per radian = {sol['slope']*np.pi/180:.4f} per degree")   # 0.3954; 4.530 per radian
print(f"C_Di = {sol['CDi']:.5f}, span efficiency = {sol['e_span']:.3f}, condition number of the matrix = {sol['cond']:.0f}")   # 0.00869, 0.954, about 40
print("Γ_n/Γ_1 for n = 1…5:", (sol["Gn"][:5]/sol["Gn"][0]).round(3))    # 1, 0, 0.122, 0, 0.026
print(f"section alone: C_L = {wd['K']*wd['alpha']:.4f}; the finite wing delivers {100*sol['CL']/(wd['K']*wd['alpha']):.0f} % of it")   # 0.548; 72 %
""", explain=r"""
1. `LL.solve_lifting_line` builds the matrix and right-hand side for N = 32, solves for the $\Gamma_n$ and returns the loading,
   the downwash and the coefficients.
2. The wing's lift coefficient is 0.3954 at α = 5° (0.0873 rad): a slope of 4.530 per radian, 72 % of the section's 2π.
3. The induced drag coefficient is 0.00869 and the span efficiency 0.954 (both explained in C13).
4. The *condition number* (how much the solve amplifies small errors — collocation, Ch. 6 P162) is modest.
5. The even coefficients vanish: a symmetric wing has a symmetric loading. Γ₃ is 12 % of Γ₁.""")
scratch(r"""
c_, s_, K_, U_, al_ = wd["c"], wd["s"], wd["K"], wd["U"], wd["alpha"]   # chord [m], span [m], lift slope [1/rad], speed [m/s], angle of attack [rad]
N4 = 4                                                       # keep four sine terms
g4 = np.array([(j + 1)*np.pi/(N4 + 1) for j in range(N4)])   # four stations γ_j = jπ/5: 36°, 72°, 108°, 144°
A4 = np.zeros((N4, N4)); r4 = np.zeros(N4)                   # empty matrix and right-hand side
for j in range(N4):                                          # one row per station
    r4[j] = K_/2*U_*c_*al_                                   # left side of Eq. (14.20): (K/2) U c (α + β), with β = 0
    for n in range(1, N4 + 1):                               # one column per sine term
        A4[j, n-1] = (1 + n*K_*c_/(4*s_*np.sin(g4[j])))*np.sin(n*g4[j])   # the bracket of Eq. (14.20) times sin(nγ_j)
Gn4 = np.linalg.solve(A4, r4)                                # the four coefficients [m²/s]
assert np.allclose(A4, LL.lifting_line_matrix(c=c_, s=s_, K=K_, N=4)[0])                       # same matrix → the library does exactly this
assert np.allclose(Gn4, LL.solve_lifting_line(c=c_, s=s_, U=U_, alpha=al_, K=K_, N=4)["Gn"])   # and the same solution
print("Γ_n (N = 4):", Gn4.round(4), "m²/s")                  # Γ1 = 0.2474, Γ3 = 0.0297; the even ones are zero
""", r"""
1. Two nested loops fill the matrix entry by entry from $\dfrac K2\,Uc(\alpha+\beta)=\sum\Big(1+\dfrac{nKc}{4s\sin\gamma}\Big)\Gamma_n\sin(n\gamma)$ (14.20).
2. `np.linalg.solve` returns the four coefficients.
3. The asserts confirm that `LL.lifting_line_matrix` and `LL.solve_lifting_line` do the same.""")
nb.code(r"""
g_chk = np.linspace(0.05, np.pi - 0.05, 400)                 # 400 span angles between the stations, to test the equation there
rows = []                                                    # one row per truncation
for N in (1, 2, 4, 8, 16, 32, 64):                           # number of sine terms kept
    sN = LL.solve_lifting_line(c=wd["c"], s=wd["s"], U=wd["U"], alpha=wd["alpha"], K=wd["K"], N=N)      # solve with N terms
    res = LL.lifting_line_residual(Gn=sN["Gn"], gamma=g_chk, c=wd["c"], s=wd["s"], U=wd["U"], alpha=wd["alpha"], K=wd["K"])   # how badly Eq. (14.17) is violated between the stations [m²/s]
    rows.append({"N": N, "C_L": sN["CL"], "C_Di": sN["CDi"], "e_span": sN["e_span"], "cond": sN["cond"], "largest residual / Γ1": np.abs(res).max()/sN["Gn"][0]})   # collect
print(pd.DataFrame(rows).round(5).to_string(index=False))    # the numbers settle after a handful of terms
""", explain=r"""
**What does this show?** A convergence table. The lift coefficient is 0.3413 with one term (the hand example), 0.3935 with
four, and 0.39535 from sixteen on; the span efficiency settles at 0.954. The last column is the largest mismatch of
$\Gamma(y)=\dfrac K2\,U\,c(y)\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.17) *between* the stations: it shrinks as terms are added (slowly — a rectangular tip
is a sharp corner). **A few terms usually suffice** for the forces.""")
nb.plotly(r"""
def loading_taper(taper):                                    # loading and downwash of a wing of aspect ratio 6 (span 6 m, area 6 m²) for one taper ratio
    c_root = 2*6.0/(6.0*(1 + taper))                         # root chord that keeps the area at 6 m²
    st = LL.solve_lifting_line(c=LL.planform(kind="tapered", s=6.0, c_root=c_root, taper=taper), s=6.0, U=1.0, alpha=wd["alpha"], K=2*np.pi, N=24)   # solve
    G1_ = st["Gn"][0]                                        # its first coefficient: the elliptic loading with the same lift has this amplitude
    return {"Γ(y) [m²/s] (clockwise)": (st["y"], st["Gamma"]),
            "elliptic loading of the same lift": (st["y"], LL.elliptic_circulation(y=st["y"], s=6.0, Gamma1=G1_)),
            "10 × w(y)/U": (st["y"], 10*st["w"]/1.0),
            "10 × uniform downwash Γ1/(2sU)": (st["y"], 10*LL.elliptic_downwash(Gamma1=G1_, s=6.0) + 0*st["y"])}   # four curves
figF7 = slider_figure(loading_taper, "taper ratio (tip chord / root chord)", np.linspace(0.1, 1.0, 19), unit="",
                      xlabel="spanwise position y [m]", ylabel="Γ [m²/s] and 10·w/U [–]", yrange=[-0.05, 0.75],
                      title="Taper reshapes the loading; near 0.4 it is almost elliptic")   # 19 slider positions
figF7.show()                                                 # display
for tp in (1.0, 0.8, 0.6, 0.4, 0.2, 0.1):                    # six taper ratios
    st = LL.solve_lifting_line(c=LL.planform(kind="tapered", s=6.0, c_root=2*6.0/(6.0*(1 + tp)), taper=tp), s=6.0, U=1.0, alpha=wd["alpha"], K=2*np.pi, N=24)   # solve
    print(f"taper {tp:3.1f}: span efficiency = {st['e_span']:.3f}, C_L = {st['CL']:.4f}")   # best near 0.4
""", explain=r"""
1. `loading_taper` solves the lifting-line equation for one taper ratio at fixed span and area, and returns the loading, the
   elliptic loading that would give the same lift (same Γ₁), ten times the downwash ratio $w/U$ and its uniform ghost.
2. The printed lines list the span efficiency for six taper ratios.""")
nb.figure_notes(
    see="Loading (first curve) against its elliptic ghost, and downwash (third curve) against its uniform ghost, as the taper "
        "ratio goes from a nearly pointed tip (0.1) to a rectangle (1.0).",
    read="A rectangle carries too much toward the tips: its loading is fuller than the ellipse and its downwash rises "
         "outboard. A sharply tapered wing does the opposite. Near a taper ratio of 0.4 the two curves nearly coincide and "
         "the downwash is nearly flat — the printed span efficiency peaks there at 0.99.",
    change="…you twisted the tips down a little (washout)? The loading would move inboard, as taper does (the live cell and N63).")
nb.live(r"""
def show_wing(taper=0.6, tip_twist_deg=0.0, AR=6.0, N=12):   # redraw for the slider values
    s_ = np.sqrt(AR*6.0)                                     # span for this aspect ratio at a fixed area of 6 m², from Λ = s²/A [m]
    c_root = 2*6.0/(s_*(1 + taper))                          # root chord that keeps the area [m]
    twist = lambda y: np.deg2rad(tip_twist_deg)*np.abs(2*y/s_)   # extra angle growing linearly from 0 at the root to the tip value [rad]
    st = LL.solve_lifting_line(c=LL.planform(kind="tapered", s=s_, c_root=c_root, taper=taper), s=s_, U=1.0, alpha=wd["alpha"], K=2*np.pi, twist=twist, N=int(N))   # solve
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(9.0, 3.2)) # left: loading; right: coefficient bars
    axl.plot(st["y"], st["Gamma"], color=COLORS["amber"], label="Γ(y)")   # the loading
    axl.plot(st["y"], LL.elliptic_circulation(y=st["y"], s=s_, Gamma1=st["Gn"][0]), color=COLORS["muted"], ls="--", label="elliptic, same lift")   # its elliptic ghost
    axl.set(xlabel="y [m]", ylabel="Γ [m²/s] (clockwise)", title=f"C_L = {st['CL']:.3f}, C_Di = {st['CDi']:.5f}, e_span = {st['e_span']:.3f}")   # three numbers
    axl.legend(fontsize=8)                                   # which curve is which
    nn = np.arange(1, len(st["Gn"]) + 1)                     # mode numbers 1…N
    axr.bar(nn, nn*st["Gn"]**2, color=[COLORS["teal"]] + [COLORS["rose"]]*(len(nn) - 1))   # n Γ_n²: the first bar carries the lift, the others only cost drag
    axr.set(xlabel="n", ylabel="n Γ_n² [m⁴/s²]", title="what each sine adds to the induced drag")   # labels
    plt.show()                                               # display
live(show_wing, taper=(0.1, 1.0, 0.1), tip_twist_deg=(-5.0, 5.0, 0.5), AR=(3.0, 20.0, 1.0), N=(1, 24, 1))   # four sliders (kernel only)
""", explain=r"""
**What does the code above do?** Four sliders — taper ratio, twist at the tip (degrees), aspect ratio and number of terms —
re-solve the lifting-line equation and redraw the loading with its elliptic ghost and the bars $n\Gamma_n^2$. The slider figure
above carries the taper part of this on the web page.""")
note("N63 [B]", r"""
**Twist.** A *rectangular* wing can carry an elliptic loading if its angle varies along the span like an ellipse:
$\alpha(y)=\alpha_t+\alpha_r\sqrt{1-(2y/s)^2}$ with tip angle $\alpha_t=\dfrac{Kc\,\alpha_r}{4s}$ and $\Gamma_1=\tfrac12KUc\,\alpha_r$ (ours; stated, then confirmed by the
solver in the next cell). Lowering the tips' angle like this is called *washout*.""")
slip(9, "the span as $b$ in the two exercises on this", "$s$ (here $b$ is the constant of the map $z=\\zeta+b^2/\\zeta$).")
nb.code(r"""
tw = LL.elliptic_twist_rectangular(alpha_r=wd["alpha"], c=wd["c"], s=wd["s"], K=wd["K"])    # the tip angle that makes a rectangle's loading elliptic
twist_fn = lambda y: tw["alpha_t"] + wd["alpha"]*np.sqrt(1 - (2*y/wd["s"])**2)              # the whole angle distribution α(y) [rad] (a lambda is a one-line function)
st_tw = LL.solve_lifting_line(c=wd["c"], s=wd["s"], U=wd["U"], alpha=0.0, K=wd["K"], twist=twist_fn, N=8)   # solve with that twist (alpha = 0: the twist carries the whole angle)
print(f"tip angle = {tw['alpha_t']:.5f} rad = {np.degrees(tw['alpha_t']):.2f}°; root angle = {np.degrees(tw['alpha_t'] + wd['alpha']):.2f}°")   # 1.31° at the tip, 6.31° at the root
print("Γ_n:", st_tw["Gn"].round(6), "m²/s; span efficiency =", round(st_tw["e_span"], 6))    # one coefficient only: [0.2742, 0, 0, …]; 1.0
""", explain=r"""
**What does this show?** With the tip at 1.31° and the root at 6.31° the solver returns a single coefficient, Γ₁ = 0.2742
m²/s $=\tfrac12KUc\,\alpha_r$, and a span efficiency of exactly 1: an elliptic loading on a rectangular wing.""")
note("N64 [C]", r"""
A small aircraft flying along a large one's tip vortex feels a rolling moment
$M=\displaystyle\int_{-s/2}^{s/2}\rho U\,y\,\Gamma(y)\,dy$ — the reason for wake-turbulence separation at airports. (The exercise that poses it takes $w$
positive *upward* — trap T5.) One number with our own inputs:""")
nb.code(r"""
rm = ch14.wake_vortex_rolling_moment(rho=1.225, U=ac["U_cruise"], K=ac["K"], s=ac["s"], c=1.25, Gamma_v=300.0, ell=1.0)   # our light aircraft in a vortex of 300 m²/s with a 1 m core
print(f"rolling moment = {rm['M']:.0f} N m (quadrature: {rm['M_quad']:.0f}); for scale, weight × half-span = {ac['W']*ac['s']/2:.0f} N m")   # comparable to the weight acting at the wing tip
""")
note("N68 [B]", r"""
**When it fails.** Lifting-line results are approximate. They need a large aspect ratio and an unswept wing (assumptions 1 and
2 of N39), a wake that stays flat (assumption 3), and sections below stall — the section law
$C_L=K\Big[\alpha-\dfrac{w(y)}{U}+\beta(y)\Big]$ (14.16) is a straight line. Extensions cover several wings, propellers, ground effect and bird formations.""")
whatif(r"""…we only want the wing's total lift and drag, not Γ(y)? They follow from the coefficients directly: C13.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C13", r"Lift and induced drag from the coefficients: $L=\dfrac{\pi s}{4}\,\rho U\Gamma_1$ (14.21), $D_i=\dfrac\pi8\,\rho\displaystyle\sum_{n=1}^{N}n\,\Gamma_n^2$ (14.22)",
     "Which parts of the loading lift the aircraft, and which only cost drag?")
problem(r"""
Split a wing's loading into its sine shapes. The first is a single hump. The others wiggle: they add lift in one place and
remove it in another, and over the whole span their lift adds up to nothing. But every wiggle sheds extra trailing vorticity,
and trailing vorticity costs drag. So the lift is in the first shape alone — and every other shape is a pure cost.""")
nb.md(r"""
#### The idea

| Term of the loading | Net lift | Induced drag |
|---|---|---|
| $n=1$ | all of it | $\propto\Gamma_1^2$ |
| $n=2$ | 0 | $\propto2\,\Gamma_2^2$ |
| $n=3$ | 0 | $\propto3\,\Gamma_3^2$ |
| $n$ | 0 | $\propto n\,\Gamma_n^2$ |""")
P("P348", "a double sum collapsed by orthogonality", r"""
Multiplying two sums gives every pairing: $(\sum_na_n)(\sum_mb_m)=\sum_n\sum_ma_nb_m$ — use two different index letters. If the pairings
are then weighted by a factor that is 1 when $n=m$ and 0 otherwise (the Kronecker delta $\delta_{nm}$), only the diagonal survives:
$\sum_n\sum_ma_nb_m\delta_{nm}=\sum_na_nb_n$. (`np.outer(p, q)` is the table of all products $p_nq_m$; `np.eye(3)` is $\delta_{nm}$ as a table.)""", r"""
p3 = np.array([1., 2., 3.]); q3 = np.array([4., 5., 6.])   # two short lists
print(np.sum(np.outer(p3, q3)*np.eye(3)), np.sum(p3*q3))   # 32.0 both: the cross terms are gone
""")
D("D17", "Lift from the first coefficient", ref="14.21")
D("D18", "Induced drag from the coefficients", ref="14.22")
slip(7, r"""the upper limit $N$ in $D_i=\dfrac\pi8\rho\displaystyle\sum_{n=1}^{N}n\Gamma_n^2$ (14.22) but ∞ in its series for Γ and $w$""",
     r"""either, used consistently: the two force formulas are written for the truncated solution of $N$ terms.""")
note("N54 [B]", r"""
For a given lift (a given Γ₁) the induced drag is least when every higher coefficient is zero. **Not in the book — ours:** the
**span efficiency** $e_{\text{span}}\equiv\dfrac{\Gamma_1^2}{\sum n\Gamma_n^2}\le1$, with $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$.

> ⚠️ **Common confusion (trap T3):** aeronautics writes this as *e*; this chapter already uses $e$ for the thickness parameter,
> so we write $e_{\text{span}}$ (`e_span` in code).""")
nb.figure(r"""
s12 = LL.solve_lifting_line(c=wd["c"], s=wd["s"], U=wd["U"], alpha=wd["alpha"], K=wd["K"], N=12)   # the rectangular demo wing with 12 terms
nn = np.arange(1, 13)                                        # mode numbers 1…12
terms = nn*s12["Gn"]**2                                      # n Γ_n²: what each sine adds to the induced drag [m⁴/s²]
fig, ax = plt.subplots(figsize=(6.2, 3.3))                   # one panel
ax.bar(nn, 100*terms/terms.sum(), color=[COLORS["teal"]] + [COLORS["rose"]]*11)   # the first bar (teal) carries the lift; the others (rose) only cost drag
ax.set(xlabel="sine term $n$", ylabel=r"share of the induced drag, $n\Gamma_n^2/\sum n\Gamma_n^2$ [%]", yscale="log", title="Lift lives in the first sine; the others are pure drag")   # logarithmic axis
plt.show()                                                   # display
print("shares for n = 1, 3, 5 [%]:", (100*terms[[0, 2, 4]]/terms.sum()).round(1), "; e_span =", round(LL.span_efficiency(Gn=s12["Gn"]), 3))   # 95.4, 4.2, 0.3; 0.954
""",
    see="One bar per sine term of the rectangular wing's loading: its share of the induced drag, on a logarithmic axis.",
    read="The teal bar ($n$ = 1) is 95.4 % of the drag and *all* of the lift. The rose bars — 4.2 % from $n$ = 3, 0.3 % from "
         "$n$ = 5 — buy no lift at all. The even terms are zero for a symmetric wing. The span efficiency is the teal share: "
         "0.954.",
    change="…the wing were tapered to 0.4? The rose bars would shrink to under 1 % in total (the explainer below shows it).")
nb.worked_example("two sines again", r"""
Γ₁ = 12, Γ₃ = 1 m²/s, $s$ = 12 m, ρ = 1.2 kg/m³, $U$ = 50 m/s.

1. $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21) = (π × 12/4) × 1.2 × 50 × 12 = 6786 N — Γ₃ does not appear.
2. $D_i=\dfrac\pi8\rho\sum n\Gamma_n^2$ (14.22) = (π/8) × 1.2 × (1 × 144 + 3 × 1) = 0.4712 × 147 = 69.3 N.
3. Without Γ₃ it would be 0.4712 × 144 = 67.9 N.
4. Span efficiency: 144/147 = 0.980 — the small third harmonic costs 2 % more drag and buys no lift.""")
nb.code(r"""
Gn = sol["Gn"]                                               # the 32 coefficients of the rectangular demo wing from C12 [m²/s]
L_ser = LL.lift_from_series(Gn=Gn, s=wd["s"], rho=wd["rho"], U=wd["U"])   # Eq. (14.21): L = (πs/4) ρ U Γ1 [N]
Di_ser = LL.induced_drag_from_series(Gn=Gn, rho=wd["rho"])   # Eq. (14.22): D_i = (π/8) ρ Σ n Γ_n² [N]
e_sp = LL.span_efficiency(Gn=Gn)                             # Γ1² / Σ n Γ_n²
CDi = LL.induced_drag_coefficient(CL=sol["CL"], AR=6.0, e_span=e_sp)   # C_L²/(π Λ e_span)
print(f"L = {L_ser:.3f} N, D_i = {Di_ser:.4f} N (ρ = 1 kg/m³, U = 1 m/s), e_span = {e_sp:.3f}, C_Di = {CDi:.5f} (solver: {sol['CDi']:.5f})")   # 1.186, 0.0261, 0.954, 0.00869
""", explain=r"""
1. `LL.lift_from_series` is $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21): it reads only the first coefficient.
2. `LL.induced_drag_from_series` is $D_i=\dfrac\pi8\rho\sum n\Gamma_n^2$ (14.22): it reads all of them.
3. `LL.span_efficiency` and `LL.induced_drag_coefficient` turn them into $e_{\text{span}}$ = 0.954 and $C_{D_i}$ = 0.00869 — the
   same numbers the solver reported in C12.""")
scratch(r"""
rho_, U_, s_ = wd["rho"], wd["U"], wd["s"]                   # density [kg/m³], speed [m/s], span [m]
gq = np.linspace(0, np.pi, 4001)[1:-1]                       # span angles, tips excluded (sin γ = 0 there)
yq = -(s_/2)*np.cos(gq)                                      # the stations in metres
Gam_q = LL.circulation_series(Gn=Gn, gamma=gq)               # Eq. (14.18): the loading at those stations [m²/s]
w_q = LL.downwash_series(Gn=Gn, gamma1=gq, s=s_)             # Eq. (14.19): the downwash there [m/s]
L_mine = np.trapezoid(rho_*U_*Gam_q, yq)                     # Eq. (14.15), first integral, by the trapezoid rule [N]
Di_mine = np.trapezoid(rho_*w_q*Gam_q, yq)                   # Eq. (14.15), second integral [N]
assert np.allclose([L_mine, Di_mine], [L_ser, Di_ser], rtol=2e-3)   # the integrals equal the two closed forms (to the accuracy of the trapezoid rule near the tips)
print(f"L = {L_mine:.4f} N, D_i = {Di_mine:.5f} N by direct integration")   # 1.186 and 0.0261
""", r"""
The two integrals $L=\int\rho U\Gamma\,dy$, $D_i=\int\rho w\Gamma\,dy$ (14.15) are evaluated directly with the trapezoid rule and compared with the two
closed forms of this block.""")
explainer("lifting_line_solver", "What makes one wing better than another?",
          "changing taper, twist and aspect ratio and watching the coefficient bars, the loading against its elliptic ghost, "
          "the downwash and the efficiency respond together shows the method's inputs and outputs at once.",
          ["Go from the rectangle to taper 0.4: what happens to the rose bars?",
           "Press 'rectangle with elliptic twist': how many bars are left?",
           "Raise N from 1 to 12 and watch the numbers settle.",
           "Click the n = 3 bar to see its shape on the span."])
whatif(r"""…we switch every rose bar off? C14.""")

# ---------------------------------------------------------------------------------------------------------------------
core("C14", r"The elliptic loading: $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26)",
     "What is the cheapest way to carry a given weight on a given span — and why are glider wings so long?")
problem(r"""
An albatross glides for hours on wings three metres across and a hand wide. A sparrow could not. The reason is in one
formula: the drag a wing pays for its lift is inversely proportional to the *square of its span*. Spread the same lift over
twice the span and the bill drops to a quarter. And among all ways of spreading the lift along a given span, one is cheapest:
the loading shaped like half an ellipse, which makes the air sink at the same speed everywhere along the wing.""")
idea(r"""
                  _.-""""""""""""-._
   Γ(y):       .'                    '.          half an ellipse: one sine of the span angle
              /                        \
   ──────────┴──────────────────────────┴──────  span
   w(y):     ────────────────────────────         flat: w = Γ1/2s everywhere

   one sine  ⇒  uniform downwash  ⇒  no part of the wing works in worse air than another
""")
D("D19", "The elliptic optimum: loading, downwash, drag and drag coefficient", ref="14.26")
note("N58 [B]", r"""
$dL=\rho U\Gamma\,dy=C_L\,\tfrac12\rho U^2c\,dy$ (14.27): with a uniform downwash, one section and no twist, every strip has the same $C_L$, so the
chord must be proportional to the circulation: an **elliptic planform** gives the elliptic loading.""")
confusion(r"""(trap T11) the optimum is a statement about the *loading* Γ(y), not about the outline. An elliptic outline
produces it only if untwisted with one section; a rectangle can produce it with twist (N63, C12); and "optimum" means least
induced drag at a fixed lift **and span**, with a flat wake.""")
D("D20", "The elliptic planform and the lift slope of a finite wing", ref="14.27")
nb.worked_example("the price of lift by hand", r"""
**N62 [C]** (the pattern of the book's level-flight exercise, with our own aircraft of C01: $W$ = 10 000 N, $s$ = 10 m, $A$ = 20 m²,
Λ = 5, ρ = 1.25 kg/m³ — a round value near the sea-level standard atmosphere of Ch. 1 — $U$ = 40 m/s, $C_L$ = 0.5.)

1. $D_i=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25) = 2 × 10⁸/(π × 1.25 × 1600 × 100) = 318 N.
2. As a coefficient, $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26) = 0.25/15.7 = 0.0159; check: 0.0159 × 1000 Pa × 20 m² = 318 N ✓.
3. First coefficient from $L=\dfrac{\pi s}{4}\rho U\Gamma_1$ (14.21): Γ₁ = 4 × 10 000/(π × 10 × 1.25 × 40) = 25.5 m²/s.
4. Downwash $w=\Gamma_1/2s$ (14.24) = 1.27 m/s, a tilt of 1.27/40 = 0.032 rad = 1.8°; lift × tilt = 10 000 × 0.0318 = 318 N ✓.
5. Circulation half-way out along each wing ($y$ = ±2.5 m): 25.5 × √(1 − 0.25) = 22.1 m²/s.
6. Double the span at the same area (Λ = 20): 318/4 = 80 N.""")
nb.code(r"""
Di_slow = LL.elliptic_induced_drag(L=ac["W"], rho=1.225, U=ac["U_slow"], s=ac["s"])       # Eq. (14.25): induced drag of our aircraft at 30 m/s [N]
Di_fast = LL.elliptic_induced_drag(L=ac["W"], rho=1.225, U=ac["U_cruise"], s=ac["s"])     # the same at 55 m/s [N]
print(f"D_i = {Di_slow:.1f} N at 30 m/s, {Di_fast:.1f} N at 55 m/s (ratio {Di_slow/Di_fast:.2f})")   # 324.8 and 96.6: it falls as 1/U²
print(f"C_Di at C_L = 1.088, Λ = 9.6: {LL.induced_drag_coefficient(CL=1.088, AR=9.6):.4f}")           # Eq. (14.26): 0.0393
print(f"lift slope: {LL.finite_wing_lift_slope(K=2*np.pi, AR=9.6):.3f} per radian at Λ = 9.6, {LL.finite_wing_lift_slope(K=2*np.pi, AR=6.0):.3f} at Λ = 6 (section: {2*np.pi:.3f})")   # 5.200 and 4.712
ew = LL.elliptic_wing(alpha=wd["alpha"], AR=6.0, K=2*np.pi)  # an untwisted elliptic wing of aspect ratio 6 at α = 5°
print(f"elliptic wing, Λ = 6, α = 5° ({wd['alpha']:.4f} rad): C_L = {ew['CL']:.4f}, C_Di = {ew['CDi']:.5f}, w/U = {ew['w_over_U']:.4f}")   # 0.4112, 0.00897, 0.0218
""", explain=r"""
1. `LL.elliptic_induced_drag` is $D_i=\dfrac{2L^2}{\pi\rho U^2s^2}$ (14.25): 324.8 N at 30 m/s and 96.6 N at 55 m/s for our aircraft.
2. Flying 1.8 times faster cuts the induced drag by 3.4: it falls as $1/U^2$.
3. `LL.induced_drag_coefficient` is $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26).
4. `LL.finite_wing_lift_slope` is the slope $K/(1+K/\pi\Lambda)$ of the derivation: 5.200 per radian at Λ = 9.6 and 4.712 (three
   quarters of 2π) at Λ = 6. `LL.elliptic_wing` collects the elliptic wing's numbers.""")
scratch(r"""
w_five = [LL.downwash(y1=y1, Gamma_fn=ell, s=12.0) for y1 in (-5.0, -2.0, 0.0, 3.0, 5.5)]    # direct quadrature of Eq. (14.13) for the half-ellipse at five stations [m/s]
assert np.allclose(w_five, LL.elliptic_downwash(Gamma1=1.0, s=12.0))                         # all equal Γ1/2s = 0.04167: the downwash is uniform
s_el = LL.solve_lifting_line(c=LL.planform(kind="elliptic", s=6.0, c_root=4/np.pi), s=6.0, U=1.0, alpha=wd["alpha"], K=2*np.pi, N=8)   # an elliptic planform of aspect ratio 6 through the solver of C12
assert np.allclose(s_el["Gn"][1:], 0, atol=1e-12)            # every coefficient but the first is zero
assert np.isclose(s_el["CL"], ew["CL"])                      # and the lift coefficient is the closed form's
print("downwash at five stations:", np.round(w_five, 5), "m/s; Γ_n of the elliptic wing:", s_el["Gn"].round(4)[:3], "…; C_L =", round(s_el["CL"], 4))   # 0.04167 ×5; [0.2618, 0, 0]; 0.4112
""", r"""
Two independent checks. (i) The principal-value integral $w(y_1)=\dfrac{1}{4\pi}\displaystyle\int\dfrac{d\Gamma}{dy}\dfrac{dy}{(y_1-y)}$ (14.13), evaluated numerically for a half-ellipse,
returns the same downwash at five different stations. (ii) The general solver, given an elliptic planform, returns one
coefficient only and the closed-form lift coefficient.""")
nb.figure(r"""
fig = fig_planforms_loadings()                               # three planforms of equal span and area (rectangle, taper 0.4, ellipse) with their Γ(y), w(y) and span efficiency
plt.show()                                                   # display
""",
    see="Three wings of the same span and area — a rectangle, a tapered wing and an ellipse — each with the loading and the "
        "downwash the lifting-line solver gives it.",
    read="Only the ellipse has a flat downwash, and its span efficiency is 1.000. The rectangle's downwash rises toward the "
         "tips (0.954); the tapered wing is close to the ellipse (0.991).",
    change="…the rectangle were twisted as in N63? Its downwash would flatten and its efficiency would reach 1 as well.")
nb.plotly(r"""
CLs = np.linspace(0.0, 1.5, 151)                             # lift coefficients on the horizontal axis
def induced_polar(AR):                                       # induced drag coefficient against lift coefficient for one aspect ratio
    return {"C_Di = C_L²/(πΛ), elliptic loading": (CLs, LL.induced_drag_coefficient(CL=CLs, AR=AR))}   # Eq. (14.26)
figF8a = slider_figure(induced_polar, "aspect ratio Λ", np.arange(3.0, 30.1, 1.0), unit="", xlabel="lift coefficient C_L [–]",
                       ylabel="induced drag coefficient C_Di [–]", yrange=[0, 0.25], title="The price of lift falls as the wing gets slender")   # 28 slider positions
figF8a.show()                                                # display
print(pd.DataFrame({"Λ": [3, 6, 9.6, 20, 30], "lift slope [1/rad]": [LL.finite_wing_lift_slope(K=2*np.pi, AR=A_) for A_ in (3, 6, 9.6, 20, 30)],
                    "C_Di at C_L = 1": [LL.induced_drag_coefficient(CL=1.0, AR=A_) for A_ in (3, 6, 9.6, 20, 30)]}).round(4).to_string(index=False))   # two numbers per aspect ratio
""", explain=r"""
1. `induced_polar` is $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26) as a curve; the slider changes the aspect ratio.
2. The table lists the finite wing's lift slope $K/(1+K/\pi\Lambda)$ and the induced drag coefficient at $C_L$ = 1.""")
nb.figure_notes(
    see="A parabola — induced drag coefficient against lift coefficient — that flattens as the aspect-ratio slider moves up.",
    read="At a fixed lift coefficient the induced drag is inversely proportional to Λ: at $C_L$ = 1 it is 0.106 for a stubby "
         "Λ = 3 and 0.016 for a glider's Λ = 20. The same slenderness also raises the lift slope toward the section's 2π.",
    change="…we add the drag that has nothing to do with lift — friction and form drag? Then there is a best speed to fly: C15.")
whatif(r"""…we add the drag that has nothing to do with lift — friction and form drag? Then there is a best speed to fly: C15.""")


# =====================================================================================================================
# A.7 §14.7 Lift and Drag Characteristics of Airfoils — C15
# =====================================================================================================================
nb.section("14.7", "Lift and Drag Characteristics of Airfoils", intro=r"""
**What is this section about?** What measured lift and drag curves look like, where they leave the theory (stall), and how the
drag of a wing splits into a part that comes from friction and shape and a part that is the price of lift.""")
core("C15", r"The drag breakdown and the polar: $C_D=C_{D0}+\dfrac{C_L^2}{\pi\Lambda e_{\text{span}}}$ (ours, built on $C_{D_i}=\dfrac{C_L^2}{\pi\Lambda}$ (14.26))",
     "Fly slowly and one kind of drag is large; fly fast and another is. Where is the best speed?")
problem(r"""
A glider pilot knows a number called the best-glide speed. Slower than that the glider sinks faster; faster than that it also
sinks faster. Two drags are fighting. Friction and the wake of the body grow with the square of the speed. The price of lift
does the opposite: the slower you fly, the harder the wing must work to hold the same weight, and the induced drag grows as
one over the square of the speed. The total has a minimum where the two are equal.""")
idea(r"""
   total drag ─┬─ friction drag (shear stress on the skin) ───────────┐
               │                                                      ├─ profile drag  C_D0  (what the section alone has)
               └─ pressure drag ─┬─ form drag (the wake of the shape) ─┘
                                 └─ induced drag (the wing's tips)   C_L² / (π Λ e_span)
""", r"""Friction and form drag were Chapter 9 (§9.8) and recap R02 in C01; the induced part is C10–C14.""")
note("N70 [B]", r"""
**Measured section curves.** A straight range $C_L\simeq K(\alpha+\beta)$, lift at α = 0 from the camber, a maximum, a drop at stall; the
drag is small and nearly flat at first, then rises roughly as α² and jumps at stall. The figure is a **qualitative model**
whose every constant is ours.""")
nb.figure(r"""
al = np.deg2rad(np.linspace(-6.0, 20.0, 261))                # angles of attack [rad] (−6° … 20°)
stall = dict(K=5.7, beta=sec["beta"], alpha_stall=np.deg2rad(14.0), CL_post=0.9, width=np.deg2rad(2.0))   # our illustrative constants: slope 5.7/rad, stall at 14°, plateau 0.9, blended over 2°
CLm = ch14.generic_lift_curve(alpha=al, **stall)             # a generic measured lift curve
CDm = ch14.generic_drag_curve(alpha=al, CD0=0.008, AR=np.inf, e_span=1.0, CD_sep=0.15, **stall)   # the matching drag curve of a section (infinite aspect ratio: no induced drag)
fig, (axl, axd) = plt.subplots(1, 2, figsize=(9.4, 3.4))     # left: lift; right: drag
axl.plot(np.degrees(al), CLm, color=COLORS["teal"])          # lift coefficient
axl.set(xlabel="angle of attack α [°]", ylabel="$C_L$ [–]", title="Qualitative model: lift of a section")   # labels
axd.plot(np.degrees(al), CDm, color=COLORS["rose"])          # drag coefficient
axd.set(xlabel="angle of attack α [°]", ylabel="$C_D$ [–]", title="Qualitative model: drag of a section — constants ours")   # labels
plt.show()                                                   # display
""",
    see="A generic lift curve and a generic drag curve of a wing section against the angle of attack.",
    read="Lift: a straight line through zero at α = −4° (the camber), a maximum near 14°, then a drop to a plateau — stall. "
         "Drag: nearly constant while the flow is attached, then a steep rise when it separates. Lift is some hundred times "
         "the drag in the attached range; that ratio is what makes wings worth having.",
    change="…the section had a flap down? The lift curve would shift up and to the left (more camber) and the maximum rise (N73).")
note("N71 [C]", r"""
Past stall the flow leaves the upper surface near the nose and the wake is as wide as the section is tall (sketch below); the
mechanism is N13 (C07).""")
nb.code(r"""
fig = draw_stalled_section()                                 # a stalled section (a sketch, not a computation): the flow leaves the upper surface near the nose
plt.show()                                                   # display
""", tags=["figure"])
note("N72 [B]", r"""
How much lift a section reaches before it stalls depends on the chord Reynolds number $\mathrm{Re}_c=Uc/\nu$: at moderate values the
laminar layer separates before it turns turbulent (early stall); at high values it turns turbulent first, and a turbulent
layer survives a steeper pressure rise — the same physics as the drag crisis of a sphere (Ch. 9 §9.9).""")
nb.code(r"""
Re_slow = ch14.chord_reynolds(U=ac["U_slow"], c=1.25, nu=nu_air)       # chord Reynolds number of our wing (mean chord 1.25 m) at 30 m/s
Re_cruise = ch14.chord_reynolds(U=ac["U_cruise"], c=1.25, nu=nu_air)   # the same at 55 m/s
print(f"Re_c = {Re_slow:.2e} at 30 m/s, {Re_cruise:.2e} at 55 m/s")    # 2.57e6 and 4.71e6
print("regime at 30 m/s:", ch14.stall_regime(Re_c=Re_slow, Re_laminar=5e5, Re_turbulent=3e6))   # 'mixed' — thresholds ours, qualitative
print(f"a rough profile drag: C_D0 ≈ {ch14.profile_drag_estimate(Re_c=Re_slow, t_over_c=0.12, regime='turbulent', form_factor=1.3):.4f} for the wing alone")   # flat-plate friction on both faces × a form factor of 1.3 (ours)
""", explain=r"""
**What does this show?** Our wing flies at chord Reynolds numbers of a few million; with our (illustrative) thresholds the
slow case is "mixed". The last line shows where a profile drag coefficient comes from: flat-plate friction on both faces
(Ch. 9) times a form factor — an **approximate** estimate for the wing alone. A whole aircraft adds fuselage, tail and
landing gear, which is why we use $C_{D0}$ = 0.025 for ours.""")
note("N73 [B]", r"""
Camber lifts the whole curve; a flap is camber you can switch on; so the slowest speed at which the wing still carries the
weight, $U_{\min}=\sqrt{\dfrac{2W}{\rho A\,C_{L,\max}}}$ (ours; step 5 of the derivation in C01), drops when the flap goes down.""")
nb.code(r"""
print(f"stall speed: {ch14.stall_speed(W=ac['W'], rho=1.225, A=ac['A'], CL_max=1.4):.1f} m/s clean, {ch14.stall_speed(W=ac['W'], rho=1.225, A=ac['A'], CL_max=2.0):.1f} m/s with the flap down")   # 26.5 and 22.1 m/s
""")
P("P349", "reading a drag polar; the minimum of a/x + bx", r"""
A polar plots $C_L$ (up) against $C_D$ (across). A straight line from the origin to a point of the curve has slope $C_L/C_D$ =
L/D; the steepest such line — the tangent — marks the best glide. Any function $p/x+qx$ has its minimum where the two terms are
equal, at $x=\sqrt{p/q}$, value $2\sqrt{pq}$ (set the derivative $-p/x^2+q$ to zero — the quotient rule of Ch. 11, P265; a numerical
minimiser was primed in Ch. 7, P170).""", r"""
p_, q_ = 0.02, 0.0637                            # C_D0 and 1/(pi*AR) for AR = 5
xs = np.linspace(0.1, 1.5, 1401)                 # candidate lift coefficients
print(xs[np.argmin(p_/xs + q_*xs)], np.sqrt(p_/q_))   # 0.56 both: the minimum of C_D/C_L
""")
remind("C15")
D("D21", "The best glide and the U-shaped drag curve", ref="14.26")
nb.worked_example("a polar by hand", r"""
Λ = 5, $e_{\text{span}}$ = 1, so πΛ = 15.7; $C_{D0}$ = 0.02.

1. Best glide where induced = profile: $C_L^2/15.7$ = 0.02 → $C_L$ = √0.314 = 0.56.
2. There $C_D$ = 0.02 + 0.02 = 0.04 and L/D = 0.56/0.04 = 14.
3. Formula: ½√(15.7/0.02) = ½ × 28 = 14 ✓.
4. Fly faster, at half that lift coefficient (0.28): induced 0.005, total 0.025, L/D = 11.2, induced share 20 %.
5. Fly slower, at twice it (1.12): induced 0.08, total 0.10, L/D = 11.2, induced share 80 %.

Same glide ratio on both sides — for opposite reasons.""")
nb.code(r"""
pol = dict(CD0=ac["CD0"], AR=9.6, e_span=ac["e_span"])       # our aircraft's polar: profile drag 0.025, aspect ratio 9.6, span efficiency 0.9
print(f"C_D at C_L = 1.088: {ch14.drag_polar(CL=1.088, **pol):.4f}; parts:", {k: round(v, 4) for k, v in ch14.drag_breakdown(CL=1.088, **pol).items()})   # 0.0687 = 0.025 profile + 0.0436 induced
for U_ in (ac["U_slow"], ac["U_cruise"]):                    # the two speeds
    dv = ch14.drag_vs_speed(W=ac["W"], rho=1.225, A=ac["A"], AR=9.6, CD0=ac["CD0"], U=U_, e_span=ac["e_span"])   # drag in newtons in level flight
    print(f"U = {U_:.0f} m/s: profile {dv['profile']:.0f} N + induced {dv['induced']:.0f} N = {dv['total']:.0f} N; induced share {100*dv['induced_share']:.1f} %, L/D = {dv['LD']:.1f}")   # 207 + 361 = 568; 695 + 107 = 802
bg = ch14.best_glide(**pol)                                  # the best-glide point of the polar
U_bg = ch14.speed_min_drag(W=ac["W"], rho=1.225, A=ac["A"], AR=9.6, CD0=ac["CD0"], e_span=ac["e_span"])   # the speed at which the drag is least [m/s]
D_bg = ch14.drag_vs_speed(W=ac["W"], rho=1.225, A=ac["A"], AR=9.6, CD0=ac["CD0"], U=U_bg, e_span=ac["e_span"])["total"]   # the drag there [N]
print(f"best glide: L/D = {bg['LD_max']:.1f} at C_L = {bg['CL']:.3f}, U = {U_bg:.1f} m/s, drag {D_bg:.0f} N, glide angle {np.degrees(bg['glide_angle']):.2f}° ({bg['glide_angle']:.4f} rad)")   # 16.5, 0.824, 34.5 m/s, 546 N
""", explain=r"""
1. `ch14.drag_polar` is $C_D=C_{D0}+C_L^2/(\pi\Lambda e_{\text{span}})$; `ch14.drag_breakdown` returns its two parts as coefficients.
2. `ch14.drag_vs_speed` turns the polar into newtons for level flight: at 30 m/s the drag is 568 N, of which 64 % is
   induced; at 55 m/s it is 802 N, of which only 13 % is induced.
3. `ch14.best_glide` gives the maximum L/D = 16.5 at $C_L$ = 0.824.
4. `ch14.speed_min_drag` is the speed where that happens for our weight: 34.5 m/s, with a drag of 546 N.
5. Both of our flight speeds are off the optimum — one on the slow side, one on the fast side.""")
scratch(r"""
CD0, AR_, e_sp_ = ac["CD0"], 9.6, ac["e_span"]               # profile drag coefficient, aspect ratio, span efficiency
CLs = np.linspace(0.05, 1.6, 15501)                          # a fine scan of lift coefficients
LD = CLs/(CD0 + CLs**2/(np.pi*AR_*e_sp_))                    # lift over drag at each one, from the polar
i = np.argmax(LD)                                            # index of the best one
assert np.isclose(LD[i], bg["LD_max"], rtol=1e-6)            # same maximum → the closed form is right
assert np.isclose(CLs[i], bg["CL"], atol=1e-4)               # at the same lift coefficient
print(f"brute force: L/D max = {LD[i]:.3f} at C_L = {CLs[i]:.4f}")   # 16.475 at 0.8238
""", r"""
A brute-force scan over 15 501 lift coefficients finds the same maximum of $C_L/C_D$ as the closed form of the derivation.""")
nb.md(r"""**What does this show?** For our aircraft $(L/D)_{\max}$ = 16.5: from 1000 m of height it glides 16.5 km in still air — and
only if it flies at 34.5 m/s.""")
note("N75 [B]", r"""
The book says that induced drag is normally the largest part of a wing's drag.

> ⚠️ **Common confusion (trap T9):** that is a statement about a *flight condition*, not about wings. The induced share of
> the drag is exactly one half at the best glide, more when slower, less when faster.""")
nb.figure(r"""
U_stall = ch14.stall_speed(W=ac["W"], rho=1.225, A=ac["A"], CL_max=ac["CL_max"])   # slowest level speed [m/s]
Us = np.linspace(U_stall, 70.0, 300)                         # speeds from the stall speed to 70 m/s
CL_U = ac["W"]/(0.5*1.225*Us**2*ac["A"])                     # the lift coefficient level flight needs at each speed
share = ch14.drag_share(CL=CL_U, **pol)                      # induced fraction of the drag
fig, ax = plt.subplots(figsize=(6.4, 3.4))                   # one panel
ax.plot(Us, 100*share, color=COLORS["rose"])                 # induced share against speed
ax.axhline(50, color="0.6", ls="--", lw=1)                   # one half
ax.axvline(U_bg, color=COLORS["accent"], ls=":", lw=1.5)     # the best-glide speed
ax.set(xlabel="flight speed $U$ [m/s]", ylabel="induced share of the drag [%]", title="'Induced drag is the largest part' — only when flying slowly")   # labels
plt.show()                                                   # display
print("induced share at 30, 34.5 (best glide), 55 m/s [%]:", (100*ch14.drag_share(CL=ac["W"]/(0.5*1.225*np.array([30.0, U_bg, 55.0])**2*ac["A"]), **pol)).round(1))   # 63.6, 50.0, 13.4
""",
    see="The fraction of our aircraft's drag that is induced, from the stall speed to 70 m/s; the dotted line is the "
        "best-glide speed.",
    read="The curve passes through exactly 50 % at the best-glide speed, 34.5 m/s. Slower than that the induced part "
         "dominates (64 % at 30 m/s, more near stall); at cruise it is a small part (13 % at 55 m/s).",
    change="…the wing were more slender? The whole curve would drop: less induced drag at every speed (the sliders below).")
nb.plotly(r"""
CLp = np.linspace(0.0, 1.5, 151)                             # lift coefficients for the polar
def polar(AR):                                               # the polar for one aspect ratio (area, weight, C_D0 and e_span fixed)
    bgA = ch14.best_glide(CD0=ac["CD0"], AR=AR, e_span=ac["e_span"])   # its best-glide point
    return {"polar: C_L against C_D": (ch14.drag_polar(CL=CLp, CD0=ac["CD0"], AR=AR, e_span=ac["e_span"]), CLp),
            "tangent from the origin (slope = best L/D)": (np.array([0.0, 2*bgA["CD"]]), np.array([0.0, 2*bgA["CL"]]))}   # curve and tangent
figF8b1 = slider_figure(polar, "aspect ratio Λ", np.arange(4.0, 24.1, 1.0), unit="", xlabel="drag coefficient C_D [–]",
                        ylabel="lift coefficient C_L [–]", xrange=[0, 0.2], yrange=[0, 1.55], title="A slender wing flattens the polar")   # 21 slider positions
figF8b1.show()                                               # display

Uv = np.linspace(22.0, 70.0, 200)                            # flight speeds [m/s]
def drag_curve(AR):                                          # drag in newtons against speed for one aspect ratio at fixed area and weight
    dv = ch14.drag_vs_speed(W=ac["W"], rho=1.225, A=ac["A"], AR=AR, CD0=ac["CD0"], U=Uv, e_span=ac["e_span"])   # the three drags
    return {"profile drag (grows as U²)": (Uv, dv["profile"]), "induced drag (falls as 1/U²)": (Uv, dv["induced"]), "total": (Uv, dv["total"])}   # three curves
figF8b2 = slider_figure(drag_curve, "aspect ratio Λ", np.arange(4.0, 24.1, 1.0), unit="", xlabel="flight speed U [m/s]",
                        ylabel="drag [N]", yrange=[0, 1500], title="Two drags trade places; the minimum is where they are equal")   # 21 slider positions
figF8b2.show()                                               # display
print(pd.DataFrame([{"Λ": A_, "best L/D": ch14.best_glide(CD0=ac["CD0"], AR=A_, e_span=ac["e_span"])["LD_max"], "best-glide speed [m/s]": ch14.speed_min_drag(W=ac["W"], rho=1.225, A=ac["A"], AR=A_, CD0=ac["CD0"], e_span=ac["e_span"])} for A_ in (4.0, 9.6, 16.0, 24.0)]).round(2).to_string(index=False))   # the numbers behind four slider positions
""", explain=r"""
1. `polar` returns the polar and its tangent from the origin; the tangent's slope is the best L/D.
2. `drag_curve` returns profile, induced and total drag in newtons against speed.
3. The table lists best L/D and best-glide speed for four aspect ratios at the same area and weight.""")
nb.figure_notes(
    see="Two slider figures for the aspect ratio: the polar with its tangent, and the U-shaped drag curve with its two parts.",
    read="A slender wing flattens the polar (the tangent steepens: better L/D) and moves the crossing of the two drag curves "
         "— the minimum of the total — to a lower speed. Best L/D grows like √Λ.",
    change="…the weight doubled? The induced curve would quadruple; the best speed would rise by √2 and the best L/D stay the same.")
note("N76 [C]", r"""
Near the speed of sound local shock waves appear on the wing; sweeping it back reduces the Mach number of the flow component
normal to the leading edge, $\mathrm{Ma}_n=\mathrm{Ma}\cos\Lambda_{\text{sweep}}$ (ours). Shock waves are Ch. 15 §15.6 and thin sections in supersonic flow
§15.10. (The sweep angle shares the letter Λ with the aspect ratio.)""")
explainer("elliptic_loading_aspect_ratio", "Why are glider wings so long and thin?",
          "stretching the span at fixed area reshapes the planform, flattens the polar and slides the minimum of the U-shaped "
          "drag curve at once; the speed cursor swaps which bar is larger.",
          ["Drag the speed from stall to cruise and watch the two bars cross.",
           "Double Λ at fixed area: what happens to the best-glide speed and ratio?",
           "Switch the mode from glider to airliner.",
           "Find the speed at which the induced share is exactly 50 %."])
whatif(r"""…the wing itself moved across the stream? Then its lift would lean *forward*. §14.8 and §14.9.""")

# =====================================================================================================================
# A.8 §14.8 Propulsive Mechanisms of Fish and Birds
# =====================================================================================================================
nb.section("14.8", "Propulsive Mechanisms of Fish and Birds", intro=r"""
**What is this section about?** A wing gives a force across the stream. Move the wing sideways through the water and "across
the stream it feels" has a forward part. A tail fin and a flapping wing are lifting surfaces used as engines.""")
note("N77 [B]", r"""
**The fish tail.** The fish swims forward at $U$ while its tail, a symmetric section, sweeps sideways at $V$. The water meets
the tail along the diagonal of that triangle, at the angle θ with $V/U=\tan\theta$. The lift $L$ is perpendicular to the diagonal:
its forward part is the thrust $L\sin\theta$, delivered at the rate $UL\sin\theta$; its sideways part $L\cos\theta$ is what the muscles work
against, at the rate $VL\cos\theta$.""")
nb.code(r"""
fig = draw_velocity_triangle(kind="fish")                    # the fish's forward speed U, the tail's sideways speed V, the relative velocity, and the lift normal to it
plt.show()                                                   # display
""", tags=["figure"])
nb.figure_notes(
    see="The velocity triangle of a tail fin: forward speed $U$, sideways speed $V$, their resultant (the water the tail "
        "meets), and the lift at right angles to that resultant, split into thrust and side force.",
    read="The lift leans forward by the same angle θ by which the relative velocity leans sideways. Forward part: "
         "$L\\sin\\theta$ (thrust). Sideways part: $L\\cos\\theta$.",
    change="…the tail swept back the other way? $V$ reverses, the triangle flips, the side force reverses — and the thrust "
           "still points forward.")
nb.worked_example("a tail stroke by hand", r"""
$U$ = 1.0 m/s, $V$ = 0.75 m/s, $L$ = 10 N.

1. tan θ = 0.75, so sin θ = 0.6 and cos θ = 0.8 (a 3–4–5 triangle; θ = 36.9° = 0.644 rad).
2. Thrust = 10 × 0.6 = 6 N; side force = 10 × 0.8 = 8 N.
3. Useful power = $U$ × thrust = 1.0 × 6 = 6 W.
4. Power put in = $V$ × side force = 0.75 × 8 = 6 W.
5. They are equal — always, because $V=U\tan\theta$ makes $VL\cos\theta=UL\sin\theta$.
6. With a drag of 1 N along the diagonal: thrust 6 − 0.8 = 5.2 N, side force 8 + 0.6 = 8.6 N, efficiency
   (5.2 × 1.0)/(8.6 × 0.75) = 0.81.""")
nb.code(r"""
ideal = ch14.oscillating_foil(U=fish["U"], V=fish["V"], L=fish["L"])                # the tail with lift only
real = ch14.oscillating_foil(U=fish["U"], V=fish["V"], L=fish["L"], D=fish["D"])    # the same with a drag of 1 N along the relative velocity
for name, r_ in (("no drag  ", ideal), ("drag 1 N ", real)):    # print both
    print(f"{name}: θ = {np.degrees(r_['theta']):.1f}° ({r_['theta']:.3f} rad), thrust {r_['thrust']:.1f} N, side force {r_['side_force']:.1f} N, power out {r_['power_out']:.2f} W, power in {r_['power_in']:.2f} W, efficiency {r_['efficiency']:.3f}")   # 6, 8, 6, 6, 1.000 ; 5.2, 8.6, 5.2, 6.45, 0.806
""", explain=r"""
**What does this show?** `ch14.oscillating_foil` reproduces the hand example: without drag the tail delivers 6 N of thrust and
all the power put in comes out as useful power; with a drag of 1 N the thrust is 5.2 N and the efficiency 0.806.""")
confusion(r"""(trap T18) the perfect score is an identity of the triangle, not a prediction: the book's picture has no drag
and no wake, and treats each instant as steady.""")
note("N79 [B]", r"""
**Bird flight.** On the downstroke the wing moves down and forward through the air; the lift, normal to the relative wind,
points up *and* forward: weight support and thrust at once. On the upstroke the wing is feathered — nearly edge-on to its
relative wind — and does little.""")
nb.code(r"""
fig = draw_velocity_triangle(kind="bird_down")               # downstroke: the wing moves down, the relative wind comes from below-ahead, the lift leans forward
fig = draw_velocity_triangle(kind="bird_up")                 # upstroke: the wing is feathered, the lift is small
plt.show()                                                   # display both sketches
down = ch14.oscillating_foil(U=10.0, V=4.0, L=2.0)           # downstroke: flight speed 10 m/s, wing moving down at 4 m/s, lift 2 N
up = ch14.oscillating_foil(U=10.0, V=4.0, L=0.2)             # upstroke: a tenth of the lift
print(f"thrust: {down['thrust']:.3f} N on the downstroke, {up['thrust']:.3f} N on the upstroke; mean over the cycle {0.5*(down['thrust'] + up['thrust']):.2f} N")   # 0.743, 0.074, 0.41 (inputs ours, arbitrary)
""", tags=["figure"], explain=r"""
**What does this show?** The same triangle as the fish tail turned on its side: with a wing speed of 4 m/s at a flight speed
of 10 m/s (θ = 21.8°), 2 N of lift gives 0.74 N of thrust on the downstroke; the feathered upstroke adds little. The inputs
are ours and arbitrary — the point is the direction of the force.""")
note("N78 [C] · N80 [C]", r"""
- **N78 [C]** Whales and dolphins move a horizontal tail up and down — the same triangle turned on its side; bumps on a
  humpback's flipper delay separation (N13, N72); flexible fins change camber within a stroke. Biological flows are Ch. 16.
- **N80 [C]** Measured bird-wing shapes and a computed flight of a fruit fly are quoted in the book; the numbers stay there.
  An insect flies at so low a Reynolds number that "thin boundary layers" (N01) no longer holds.""")

# =====================================================================================================================
# A.9 §14.9 Sailing against the Wind — C16
# =====================================================================================================================
nb.section("14.9", "Sailing against the Wind", intro=r"""
**What is this section about?** The same velocity triangle a third time: a sail is a wing standing on end in the wind the
boat *feels*. As long as that wind comes from far enough off the bow, the sail's lift has a forward part — even when the boat
is heading into the wind.""")
core("C16", r"Thrust from lift: $\mathbf U_r=\mathbf U-\mathbf V$, drive $=L\sin\beta_a-D\cos\beta_a$ (ours)",
     "How can a boat sail against the wind that drives it?")
problem(r"""
"The wind pushes the sail" explains sailing downwind and nothing else. A sail trimmed close is a wing: air flows *along* both
sides and the force is mostly sideways to the wind the boat feels. That felt wind — the apparent wind — is the true wind minus
the boat's own velocity, so it always comes from further ahead than the true wind. Lift is perpendicular to it; and a force
perpendicular to a wind from ahead-and-to-the-side points partly forward. The sideways rest of the force would shove the boat
downwind, were it not for a second wing under the water: the keel.""")
idea(r"""
   plan view, the boat heading up the page

                true wind U  ↙  (from ahead and to the right, at θ off the bow)
        apparent wind U_r = U − V  ↓↙  (closer to the bow, at β_a < θ, and stronger)

                      ▲  V (boat velocity)
                     ╱█╲
        L  ◄────────  █      lift ⟂ U_r : points forward-and-left
     (teal)           █  ↘ D drag ∥ U_r
                     ╲█╱
        drive  =  forward part of L  −  backward part of D
""", r"""In this block θ is the **course angle** off the true wind, and $\beta_a$ the angle of the apparent wind off the bow (not
the camber angle).""")
remind("C16")
D("D22", "Sailing: the apparent wind, the drive, and the closest course", ref="14.15")
nb.worked_example("a beam reach by hand", r"""
True wind $U$ = 4 m/s from exactly abeam (course angle 90°), boat speed $V$ = 3 m/s.

1. The apparent wind has the part $U$ = 4 m/s from the side and the part $V$ = 3 m/s from ahead: speed 5 m/s (3–4–5).
2. Its angle off the bow: tan $\beta_a$ = 4/3, $\beta_a$ = 53.1°, sin 0.8, cos 0.6.
3. The sail gives $L$ = 100 N and $D$ = 25 N.
4. Drive = $L\sin\beta_a-D\cos\beta_a$ = 80 − 15 = 65 N.
5. Side force = $L\cos\beta_a+D\sin\beta_a$ = 60 + 20 = 80 N, carried by the keel.
6. The drive would vanish at tan $\beta_a$ = $D/L$ = 0.25, i.e. 14° off the bow.
7. Head up to 45° from the true wind: tan $\beta_a$ = 2.83/(2.83 + 3) = 0.485, $\beta_a$ = 25.9° — still more than 14°, so the boat
   still drives forward, *toward* the wind.""")
nb.code(r"""
aw = ch14.apparent_wind(U=boat["U_wind"], V=boat["V_boat"], theta=boat["theta"])   # true wind 6 m/s, boat 3 m/s, course 45° (0.785 rad) off the wind
sf = ch14.sail_forces(U=boat["U_wind"], V=boat["V_boat"], theta=boat["theta"], CL=boat["CL"], CD=boat["CD"], rho=1.225, A=boat["A_sail"])   # the sail's forces in that apparent wind
lim = ch14.sailboat_course_limit(LD_sail=boat["CL"]/boat["CD"], LD_hull=boat["LD_hull"])   # the closest course to the apparent wind: ε_A + ε_H [rad]
eq = ch14.sailboat_equilibrium(U=boat["U_wind"], theta=boat["theta"], CL=boat["CL"], CD=boat["CD"], rho_air=1.225, A_sail=boat["A_sail"], rho_water=boat["rho_water"], CDA_hull=boat["CDA_hull"], LD_hull=boat["LD_hull"])   # steady sailing: drive = hull drag (ours)
print(f"apparent wind: {aw['Ur']:.2f} m/s at {np.degrees(aw['beta_a']):.1f}° ({aw['beta_a']:.3f} rad) off the bow")   # 8.39 m/s at 30.4°
print(f"sail: L = {sf['L']:.0f} N, D = {sf['D']:.0f} N → drive {sf['drive']:.0f} N, side force {sf['lateral']:.0f} N")   # 1036, 173 → 375, 981
print(f"closest course to the apparent wind: {np.degrees(np.arctan(boat['CD']/boat['CL'])):.2f}° + {np.degrees(np.arctan(1/boat['LD_hull'])):.2f}° = {np.degrees(lim):.1f}° ({lim:.3f} rad)")   # 9.46° + 7.13° = 16.6°
print(f"steady sailing on this course: V = {eq['V']:.2f} m/s, apparent wind at {np.degrees(eq['beta_a']):.1f}°, drive {eq['drive']:.0f} N = hull drag, side force {eq['lateral']:.0f} N")   # 3.95 m/s, 27.4°, 390 N, 1207 N
print(f"the keel must work at L/D ≥ {eq['LD_hull_required']:.1f} (it has {boat['LD_hull']:.0f}: holds = {eq['keel_holds']}); ε_A + ε_H = {np.degrees(eq['eps_A']):.1f}° + {np.degrees(eq['eps_H']):.1f}° = β_a (closure {eq['closure']:.0e})")   # 3.1; 9.5° + 17.9°
print(f"velocity made good toward the wind at 3 m/s on this course: {ch14.velocity_made_good(V=boat['V_boat'], theta=boat['theta']):.2f} m/s")   # 2.12 m/s
""", explain=r"""
1. `ch14.apparent_wind`: sailing at 3 m/s at 45° to a 6 m/s wind, the boat feels 8.39 m/s from 30.4° off the bow — stronger
   and further ahead than the true wind.
2. `ch14.sail_forces`: a 20 m² sail with $C_L$ = 1.2 and $C_D$ = 0.2 gives 1036 N of lift and 173 N of drag; resolved along
   and across the course that is a drive of 375 N and a side force of 981 N.
3. `ch14.sailboat_course_limit`: the apparent wind cannot be closer to the bow than the sum of the two drag angles,
   9.46° + 7.13° = 16.6°, in the best case (a keel working at its best L/D of 8).
4. `ch14.sailboat_equilibrium` (ours) finds the boat speed at which the drive equals the hull's drag: 3.95 m/s on this
   course — a little faster than the 3 m/s assumed above. There the keel needs an L/D of only 3.1, and the polygon closes:
   $\beta_a=\varepsilon_A+\varepsilon_H$ = 9.5° + 17.9° = 27.4°.
5. `ch14.velocity_made_good` is the part of the boat's speed that points toward the wind: 2.12 m/s at 3 m/s on a 45° course.""")
scratch(r"""
Uw, V_, th = boat["U_wind"], boat["V_boat"], boat["theta"]   # true wind speed [m/s], boat speed [m/s], course angle off the wind [rad]
t_hat = np.array([np.cos(th), np.sin(th)])                   # unit vector along the course (axes: the wind blows toward −x)
n_hat = np.array([np.sin(th), -np.cos(th)])                  # unit vector across the course, to leeward
U_vec = np.array([-Uw, 0.0])                                 # the true wind as a vector
Ur = U_vec - V_*t_hat                                        # the apparent wind: true wind minus the boat's velocity
ur = Ur/np.hypot(*Ur)                                        # its direction: the drag points along it
lift_dir = np.array([-ur[1], ur[0]])                         # that direction turned by 90° …
if lift_dir @ t_hat < 0: lift_dir = -lift_dir                # … toward the bow (a @ b multiplies matching components and adds them: the dot product)
F_vec = sf["L"]*lift_dir + sf["D"]*ur                        # the sail's total force
drive, lateral = F_vec @ t_hat, F_vec @ n_hat                # its parts along and across the course
assert np.allclose([drive, lateral], [sf["drive"], sf["lateral"]])   # same numbers → the library does exactly this
print(f"drive = {drive:.1f} N, side force = {lateral:.1f} N; apparent wind {np.hypot(*Ur):.2f} m/s")   # 374.6, 980.9, 8.39
""", r"""
1. Everything is done with two-component vectors: the course direction, the leeward direction, the true wind and the boat's
   velocity.
2. The apparent wind is their difference; drag points along it and lift at right angles, on the side that faces the bow.
3. Dot products (`@`) with the course direction and the leeward direction give drive and side force — the same numbers as
   `ch14.sail_forces`.""")
nb.figure(r"""
fig = draw_sail_polygon(sf=sf)                               # force diagram with our numbers: apparent wind, lift, drag, their resultant split into drive and side force
plt.show()                                                   # display
""",
    see="The force diagram of our sail: the apparent wind (blue), the lift at right angles to it (teal), the drag along it "
        "(rose), and their resultant split into drive (along the course) and side force (across it).",
    read="The resultant leans forward of the beam by $\\beta_a-\\varepsilon_A$: the apparent-wind angle minus the sail's drag angle. "
         "That lean is the drive. The much larger side force is what the keel must cancel.",
    change="…the sail were baggy, with L/D = 3? Its drag angle would open from 9.5° to 18.4° and the boat could not point as high.")
note("N82 [B]", r"""
**Tacking.** No boat sails straight into the wind. It zig-zags at a fixed angle to it, and what counts is the speed made good
toward the destination, $V\cos\theta$ (ours).""")
nb.figure(r"""
th_deg = np.linspace(5.0, 175.0, 171)                        # course angles off the true wind [degrees]
drive_curve = ch14.sail_drive_curve(theta=np.deg2rad(th_deg), U=boat["U_wind"], V=boat["V_boat"], CL=boat["CL"], CD=boat["CD"], rho=1.225, A=boat["A_sail"])   # drive [N] on each course at a boat speed of 3 m/s
th_zero = np.interp(0.0, drive_curve[:40], th_deg[:40])      # the course at which the drive vanishes (the curve rises steadily there)
fig, (axp, axd) = plt.subplots(1, 2, figsize=(9.4, 3.6))     # left: a tacking path; right: drive against course
legs = np.array([[0, 0], [3, 3], [0, 6], [3, 9], [0, 12]])   # a zig-zag at 45° either side of the wind direction (the wind blows down the page)
axp.fill([0, -12*np.tan(np.deg2rad(th_zero)), 12*np.tan(np.deg2rad(th_zero))], [0, 12, 12], color="0.88")   # the no-go wedge either side of the wind direction
axp.plot(legs[:, 0], legs[:, 1], "-o", color=COLORS["accent"])   # the path
axp.annotate("wind", xy=(-4, 8.5), xytext=(-4, 11.5), arrowprops=dict(arrowstyle="->", color=COLORS["blue"]), color=COLORS["blue"], ha="center")   # the wind blows down the page
axp.set(xlim=(-7, 7), ylim=(-0.5, 12.5), xlabel="across the wind [boat lengths]", ylabel="toward the wind [boat lengths]", title="Tacking: progress upwind in a zig-zag")   # labels
axp.set_aspect("equal")                                      # true proportions
axd.plot(th_deg, drive_curve, color=COLORS["accent"])        # drive against course angle
axd.axhline(0, color="k", lw=0.8)                            # zero drive
axd.axvspan(0, th_zero, color="0.88")                        # courses on which the sail cannot drive the boat
axd.set(xlabel="course angle θ off the true wind [°]", ylabel="drive [N]", title="No drive closer than a certain angle")   # labels
plt.show()                                                   # display
print(f"the drive vanishes at θ = {th_zero:.1f}° off the true wind (at V = 3 m/s, sail L/D = 6)")   # the edge of the no-go zone for the sail alone
""",
    see="Left: a boat's zig-zag path toward the wind, with the wedge it cannot sail in shaded. Right: the sail's drive "
        "against the course angle at a boat speed of 3 m/s; the shaded band is where the drive is negative.",
    read="The drive is zero where the apparent wind is at the sail's drag angle (9.5° off the bow), which for this boat speed "
         "is the printed course angle off the *true* wind. From there it rises steeply, peaks with the wind somewhat ahead "
         "of the beam, and falls again toward a run. To go upwind the boat sails either side of the wedge and tacks.",
    change="…the boat were faster? The apparent wind would swing further ahead and the wedge would widen: fast boats tack "
           "through wider angles.")
nb.animation(r"""
n_frames = 60 if not FAST else 30                            # number of frames
fig, (axf, axs) = plt.subplots(1, 2, figsize=(9.2, 3.6))     # left: a fish tail through one stroke; right: a boat tacking upwind
axf.set(xlim=(-1.6, 1.6), ylim=(-1.3, 1.3), title="fish tail: thrust forward on both half-strokes", xlabel="swimming direction ←", ylabel="sideways")   # fixed axes
axf.set_aspect("equal")                                      # true proportions
(tail,) = axf.plot([], [], color="k", lw=3)                  # the tail fin, seen from above
qf = axf.quiver([0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], color=[COLORS["blue"], COLORS["teal"], COLORS["accent"]], angles="xy", scale_units="xy", scale=1, width=0.012)   # relative velocity, lift, thrust
axs.set(xlim=(-5, 5), ylim=(-0.5, 12.5), title="sail: drive forward on both tacks", xlabel="across the wind", ylabel="toward the wind ↑ (the wind blows ↓)")   # fixed axes
axs.set_aspect("equal")                                      # true proportions
(track,) = axs.plot([], [], color="0.6", lw=1)               # the path sailed so far
qs = axs.quiver([0, 0], [0, 0], [0, 0], [0, 0], color=[COLORS["blue"], COLORS["accent"]], angles="xy", scale_units="xy", scale=1, width=0.012)   # apparent wind and drive
path_x = lambda p: 3*np.abs(((p/3 + 1) % 2) - 1) - 1.5       # zig-zag: across-wind position for a distance p made good toward the wind
fig.canvas.draw(); fig.set_layout_engine("none")             # lay the figure out once, then freeze it
def update(i):                                               # draw frame i
    ph = 2*np.pi*i/n_frames                                  # phase of the tail stroke
    Vt = fish["V"]*np.cos(ph)                                # the tail's sideways velocity now [m/s]
    yt = 0.6*np.sin(ph)                                      # its sideways position
    of = ch14.oscillating_foil(U=fish["U"], V=abs(Vt) + 1e-9, L=fish["L"]*abs(Vt)/fish["V"])   # quasi-steady: lift ∝ the sideways speed
    sgn = np.sign(Vt) if Vt != 0 else 1.0                    # which way the tail is moving
    tail.set_data([0.0, 0.5], [yt, yt - 0.25*sgn*np.sin(of["theta"])])   # the fin, turned a little into its relative flow
    qf.set_offsets(np.c_[[0.0, 0.0, 0.0], [yt, yt, yt]])     # all three arrows start at the fin
    qf.set_UVC([fish["U"], -of["thrust"]/10, -of["thrust"]/10], [-Vt, -sgn*of["side_force"]/10, 0.0])   # relative water velocity; lift (thrust forward, side force against the motion); thrust alone
    p = 12*i/n_frames                                        # distance made good toward the wind
    pp = np.linspace(0, p, 80)                               # the path so far
    track.set_data(path_x(pp), pp)                           # draw it
    tack = 1.0 if ((p/3 + 1) % 2) > 1 else -1.0              # which tack: heading up-right (+1) or up-left (−1)
    head = np.array([tack*np.sin(boat["theta"]), np.cos(boat["theta"])])   # unit vector along the course, 45° off the wind direction
    wind_app = np.array([0.0, -boat["U_wind"]]) - boat["V_boat"]*head      # apparent wind = true wind (blowing down the page) − boat velocity
    qs.set_offsets(np.c_[[path_x(p), path_x(p)], [p, p]])    # both arrows start at the boat
    qs.set_UVC([wind_app[0]/4, head[0]*sf["drive"]/250], [wind_app[1]/4, head[1]*sf["drive"]/250])   # apparent wind (scaled) and drive (scaled)
show_animation(animate(update, frames=n_frames, fig=fig, interval=90), player="video")   # a smooth video
""", explain=r"""
1. Left: at each instant of the stroke the tail's sideways velocity sets the velocity triangle; `ch14.oscillating_foil`
   returns thrust and side force for it (quasi-steady: each instant treated as steady, with the lift taken proportional to
   the sideways speed). Arrows: relative water velocity (blue), lift (teal), thrust (purple).
2. Right: a boat zig-zags upwind at 45° either side of the wind; the arrows are the apparent wind (blue) and the drive
   (purple) from `ch14.sail_forces`, drawn along the course.
3. The arrows are scaled for the drawing; their directions are what matters.""")
nb.figure_notes(
    see="Left: a tail fin sweeping from side to side. Right: a boat tacking toward the wind.",
    read="The teal lift arrow swings from one side to the other with the stroke, but its forward part — the purple thrust — "
         "points the same way on both half-strokes. On the right the drive arrow likewise points along the course on both "
         "tacks, although the apparent wind comes alternately over the right and the left bow.",
    change="…the tail moved faster for the same swimming speed? θ grows, the lift leans further forward: more thrust per "
           "newton of lift, but more side force to pay for in a real, lossy fin.")
nb.plotly(r"""
th_c = np.linspace(2.0, 178.0, 177)                          # course angles off the true wind [degrees]
beta_c = np.array([ch14.apparent_wind(U=1.0, V=0.5, theta=np.deg2rad(t_))["beta_a"] for t_ in th_c])   # apparent-wind angle off the bow for a boat at half the wind speed [rad]
def drive_per_lift(LD):                                      # drive per unit lift against course angle for one sail lift-to-drag ratio
    d_ = np.sin(beta_c) - np.cos(beta_c)/LD                  # (L sin β_a − D cos β_a)/L with D/L = 1/LD
    return {"drive / lift": (th_c, d_), "no drive below this line": (th_c, 0*th_c)}   # the curve and the zero line
figF9 = slider_figure(drive_per_lift, "sail L/D", np.linspace(2.0, 12.0, 21), unit="", xlabel="course angle θ off the true wind [°]",
                      ylabel="drive per unit lift [–]", yrange=[-0.55, 1.05], title="A better sail narrows the wedge you cannot sail in (V/U = 0.5)")   # 21 slider positions
figF9.show()                                                 # display
print(pd.DataFrame([{"sail L/D": LD_, "drag angle ε_A [°]": np.degrees(np.arctan(1/LD_)), "drive vanishes at θ [°]": np.interp(0.0, np.sin(beta_c) - np.cos(beta_c)/LD_, th_c)} for LD_ in (2.0, 4.0, 6.0, 12.0)]).round(1).to_string(index=False))   # where the curve crosses zero
""", explain=r"""
1. `beta_c` is the apparent-wind angle on every course for a boat moving at half the wind speed (`ch14.apparent_wind`).
2. `drive_per_lift` is the drive of the derivation divided by the lift, $\sin\beta_a-\cos\beta_a/(L/D)$.
3. The table lists, for four sails, the drag angle and the course on which the drive vanishes.""")
nb.figure_notes(
    see="Drive per unit of lift against the course angle, with a slider for the sail's lift-to-drag ratio.",
    read="The curve crosses zero where the apparent wind is at the sail's drag angle. A poor sail (L/D = 2) cannot drive the "
         "boat closer than about 40° to the true wind at this speed; a good one (L/D = 12) gets within about 7°. The hull's "
         "own drag angle adds to this (step 8 of the derivation).",
    change="…the boat sailed faster relative to the wind? Every curve shifts right: the apparent wind moves ahead.")
explainer("thrust_from_lift", "How can a boat sail against the wind?",
          "turning the course toward the wind and watching the apparent wind swing forward, the lift rotate and the drive "
          "shrink to zero is the whole argument; switching to the fish tail or the bird's wing shows the same triangle with "
          "stream and motion exchanged.",
          ["Head up until the drive bar reaches zero and read $\\beta_a$.",
           "Give the sail a perfect L/D: where is the limit now?",
           "Switch to 'fish tail' and press play.",
           "In 'bird wing', compare the downstroke with the upstroke."])
whatif(r"""…the air were compressible? Then pressure itself travels at a finite speed and a wing can outrun its own disturbance:
Chapter 15.""")


# =====================================================================================================================
# End matter
# =====================================================================================================================
nb.pointer(r"""**S01** The chapter's exercises are not reproduced. The results the chapter relies on were taught above: the
pressure loop (C07), the flat plate (C02, C06), the circular arc (C04), leading-edge suction (C07), the half-line vortex (C09),
the horseshoe wing, its velocity field and ground effect (C09), the parabolic and elliptic loadings (C09, C14), the twisted
wing (C12), the finite-wing lift slope (C14), the rolling moment (C12), tip-vortex cores (C10), the edge of finite angle (C04)
and the projectile (C01). **S02** Literature: see the chapter's reference list in the book; the public sources used for
cross-checks are listed in `analysis/ch14.md` §8.""")
nb.summary(
    clicked=[
        "**C01** A wing is a span, an area and their ratio; in level flight the lift coefficient it needs is its weight over dynamic pressure times area.",
        "**C02** Ideal flow allows any circulation; the sharp trailing edge allows one.",
        "**C03** The wing gets its circulation by shedding the opposite one — the total never leaves zero.",
        "**C04** Raise the circle and the plate bends; move it back and the nose fattens; the tail is the point the map folds.",
        "**C05** Same potential at partner points; the speed is rescaled by the map, infinitely at the tail.",
        "**C06** Lift grows by 2π per radian of angle; camber shifts the line without tilting it.",
        "**C07** Lift is the area between two pressure curves, mostly suction near the nose.",
        "**C08** Wherever the circulation drops along the span, a vortex trails off.",
        "**C09** Each trailing filament is half a line vortex; together they make the air at the wing sink.",
        "**C10** Lift is perpendicular to sinking air, so it leans back: drag without friction.",
        "**C11** Measure the span by an angle and the loading by sines, and the singular integral is done once.",
        "**C12** N stations, N sines, one linear solve: any planform, any twist.",
        "**C13** Lift lives in the first sine; every other one only costs drag.",
        "**C14** One sine means uniform downwash and the least drag, and the price falls with the square of the span.",
        "**C15** Profile drag rises with speed, induced drag falls; the best glide is where they are equal.",
        "**C16** Move a wing across the stream and its lift leans forward: a fish tail, a wing beat, a sail.",
    ],
    feeds_forward=[
        r"**Ch. 15 (compressible flow)** — the lift slope 2π against the supersonic thin-airfoil value; the suction peak and the critical Mach number; finite wings with a compressibility-corrected slope $K$; shock waves and sweep.",
        "**Ch. 16 (biofluid mechanics)** — fins, flapping and low Reynolds numbers.",
        "**Wind-turbine and propeller blades** — the same section in a rotating relative wind (the velocity triangle of C16 once more).",
    ],
    left_out=[
        "thin-airfoil theory and Theodorsen's method (named in C04; standard aerodynamics texts)",
        "the Kármán–Trefftz section (named in C06)",
        "viscous prediction of stall and maximum lift (C07, C15 — qualitative only)",
        "unsteady aerodynamics beyond the lumped-vortex sketch (C03)",
        "swept and low-aspect-ratio wings (limits listed in C12)",
        "the measured data behind the generic curves (C15)",
    ],
)

# ---------------------------------------------------------------------------------------------------------------------
# Self-checks before saving
# ---------------------------------------------------------------------------------------------------------------------
_ALL = "\n".join(c.source for c in nb.cells)
_planned = curation_items("ch14")
_missing = [k for k in _planned if not re.search(r"\b" + k + r"\b", _ALL)]
print(f"curation ids: {len(_planned)}; not named in the notebook: {_missing}")
_EQ_REF = re.compile(r"(?:\bEqs?\.\s*\(?|\()((?:\d+|[A-D])\.\d+)(?=[\s,)–-])")
_bare = [f"cell {i}: {sorted(set(_EQ_REF.findall(c.source)))} :: {c.source[:70]!r}" for i, c in enumerate(nb.cells)
         if c.cell_type == "markdown" and _EQ_REF.search(c.source) and "$" not in c.source]
_ctrl = [f"cell {i}" for i, c in enumerate(nb.cells) if any(ch in c.source for ch in "\x07\x08\x0b\x0c")]
_unc = []
for _i, _c in enumerate(nb.cells):
    if _c.cell_type != "code" or "setup" in _c.metadata.get("tags", []):
        continue
    _ls = [l for l in _c.source.splitlines() if l.strip()]
    _no = [l for l in _ls if "#" not in l and l.strip() not in (")", "]", "}", "])", "})") and not l.strip().startswith(("else:", "try:"))]
    if _no:
        _unc.append(f"cell {_i}: {len(_no)}/{len(_ls)} lines without a comment, e.g. {_no[0][:80]!r}")
_prim = [p.lower() for p in nb.primers]
_led = []
for _concept, _first, _by, *_r in LEDGER:
    if "primer" in _by.lower():
        _key = re.sub(r"[`*$\\]", "", _concept).lower().split("(")[0].strip()
        if _key and not any(_key[:18] in p or p[:18] in _key for p in _prim):
            _led.append(_concept)
for _m in _unc:
    print("UNCOMMENTED", _m)
for _m in _led:
    print("LEDGER primer/reminder missing:", _m)
_kinds = {k: sum(k in c.metadata.get("tags", []) for c in nb.cells) for k in
          ("figure", "animation", "plotly", "live-only", "explainer", "primer", "derivation", "derivation-check", "from-scratch", "recap")}
if "--partial" in sys.argv:                                  # development aid: write what exists so far, unchecked
    import nbformat as _nbf
    for _m in _bare + _ctrl:
        print("CHECK", _m)
    _nb2 = _nbf.v4.new_notebook(cells=list(nb.cells))
    _nb2.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    _nb2.metadata["fluidpy"] = {"chapter": "ch14", "explainers": nb.explainers, "cores": sorted(nb.cores),
                                "recaps": sorted(nb.recaps), "primers": nb.primers, "derivations": nb.derivations}
    for _i, _c in enumerate(_nb2.cells):
        _c["id"] = f"ch14-{_i:03d}"
    _out = ROOT / "notebooks" / "ch14_aerodynamics.ipynb"
    _nbf.write(_nb2, _out)
    print("wrote (partial, unchecked)", _out.name, len(nb.cells), "cells", _kinds)
    sys.exit(0)
if _bare:
    raise SystemExit("markdown cells cite equation numbers without maths:\n  " + "\n  ".join(_bare))
if _ctrl:
    raise SystemExit("control characters (an eaten LaTeX backslash) in: " + ", ".join(_ctrl))
if _missing:
    raise SystemExit("curation ids with no home in the notebook: " + ", ".join(_missing))
if _led:
    raise SystemExit("Part E rows explained by a primer but no primer/reminder in the notebook:\n  " + "\n  ".join(_led))
out = nb.save()
print(f"wrote {out.relative_to(ROOT)} ({len(nb.cells)} cells: "
      f"{sum(c.cell_type == 'markdown' for c in nb.cells)} markdown, {sum(c.cell_type == 'code' for c in nb.cells)} code) {_kinds}")
print(f"derivations: {len(nb.derivations)} with {sum(d['steps'] for d in nb.derivations.values())} steps, "
      f"{sum(d['check'] for d in nb.derivations.values())} sympy checks; primers: {len([p for p in nb.primers if not p.endswith('(reminder)')])} "
      f"new + {len([p for p in nb.primers if p.endswith('(reminder)')])} reminders; explainers: {nb.explainers}")
