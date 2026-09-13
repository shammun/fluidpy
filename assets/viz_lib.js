/* ======================================================================================================
   viz_lib.js — the engine behind every fluidpy interactive explainer (viz/chNN/<slug>.html).

   One call builds a whole explainer that FITS THE WINDOW (no page scroll, no panel scroll):

     const app = Viz.app({
       title, subtitle,
       params:    { k: {label, min, max, step, value, unit, help, log, fmt, optional}, deep: {type:'toggle', ...},
                    mode: {type:'select', options:[['a','Label A'], ...], value, display:'chips'?}, kick: {type:'button', action} },
       readouts:  [ {id, label, unit, value: s => number, fmt?, tone?: v => 'pos'|'neg'|''} ],
       stage:     { animate, setup(g), step(g, s, dt), draw(g, s), reset(g, s), onPointer(g, ev, s),
                    views: [ {id, title: s => html, row: 0, flex: 1, hidePortrait?, draw(v, s, app)?, onPointer(v, ev, s, app)?} ],
                    rows: [1.2, 1] },                               // linked views sharing one state and one clock
       transport: { param:'t', min:0, max:6.28, rate: 1, end:'hold'|'loop'|'stop', hold: 2, step },  // ▶ ⏮ ⏭ ↺ scrubber speed
       modes:     { param:'system', options:[['wheel','Rotating wheel'], ['spring','Spring–mass']] },   // same idea, other system
       presets:   [ {label:'ω = 1', set:{w:1}, title} ],             // one-click special cases (chips above the stage)
       status:    s => ({text:'🎯 resonance', tone:'good'|'warn'|'bad'|'info'}),
       explore:   { intro, controls:['k','deep'], readouts:['c'], callouts:[{kind:'try'|'key'|'watch'|'warn', html}] },
       terms:     { title, items:[{id, label, color, value: s => number}], total:{label}, unit },  // click a term to focus it
       inspect:   (s, app) => html | null,                           // click-to-trace exact arithmetic (state set by onPointer)
       notes:     (s, app) => html,                                  // "Right now": regime-dependent interpretation
       explain:   { title, html: (s, app) => Viz.work.step(1,'…')+Viz.work.line(...)+Viz.work.interpret(...), live: (s, app) => ({name: value}) },
                                                                     // Explain tab ("Explanation & interpretation"; old name: calc)
       derivations: [ {id, title, short, ref:'Eq. (7.27)', goal, start:{tex, plain}, plan:['…'], uses:['product rule'],
                       steps:[ {tex, did:'divide by m', why:'…allowed because…', plain:'…in words…', live: s => tex, set, highlight, watch} ],
                       result:{tex, plain, set}, interpret: (s, app) => html, check:'units / limit / sympy', picture: true} ],
       tour:      [ {title, text, set, play, controls:['k'], readouts:['c'], eq:'disp', code:{id, lines:[3,4]}, derive:{id, step},
                     terms:true, inspect:true, notes:true, callout:{kind, html}, highlight:['readout:c'], enter(app)} ],
       equations: [ {id, title, ref:'Eq. (7.27)', tex, live: s => tex, note, symbols:[[tex, meaning, unit]]} ],
       code:      [ {id, title, ref, src:'python with {{live}} placeholders', live: s => ({live: value}), note} ],  // Code tab
       check:     [ {q, a, set?:{...}} ],
       panels:    [ {id, label, short, render(el, app)} ],          // optional extra tabs
       selftest:  () => [ {name, js, py?, expect?, rtol?, atol?} ], // parity with fluidpy (py) or JS-only (expect)
       onChange(s, key, app)
     });
   Helpers: Plot (g.plot / v.plot), Viz.field (grid, contours, heatmap, streamline, quiver, particles), Viz.num,
   Viz.work (head, line, result, table), Viz.live(name), Viz.three(view, opt) → three.js scene with orbit + pick.

   Layout: header (title | tabs | actions) + main (stage [top strip · views · transport] | side panel of the active tab).
   Fit algorithm (Viz.fit): choose data-layout (wide | portrait | landscape) from the box size → raise data-dense
   0..3 until nothing overflows → paginate long lists (equations, questions, controls) into pages → if anything
   still overflows, console.warn('VIZ-OVERFLOW', …) so tools/shot.py fails the build.

   Audit hooks for tools/shot.py: window.VIZ = { app, ready, audit(), selftest(), setTab(id), goStep(i), goDerive(i, s),
   tabs, steps, derivations, explainStats(), features }.
   URL hash deep links: #tab=equations&step=3&k=0.25   ·   #tab=derive&d=1&ds=4

   Inlined verbatim into each explainer by tools/viz_inline.py between the VIZ_LIB_JS markers — edit THIS file,
   then run `python tools/viz_inline.py --all`. No dependencies; KaTeX is loaded from a CDN with a text fallback.
   ====================================================================================================== */
/* VIZ_LIB_JS:BEGIN */
(function (global) {
  'use strict';
  var Viz = { version: '1.0.0' };
  var doc = global.document;

  /* ------------------------------------------------------------------------------------------------
     0. small utilities
     ------------------------------------------------------------------------------------------------ */
  function h(tag, attrs) {
    var el = doc.createElement(tag);
    if (attrs) {
      for (var k in attrs) {
        if (!Object.prototype.hasOwnProperty.call(attrs, k)) continue;
        var v = attrs[k];
        if (v === null || v === undefined || v === false) continue;
        if (k === 'class') el.className = v;
        else if (k === 'html') el.innerHTML = v;
        else if (k === 'text') el.textContent = v;
        else if (k.slice(0, 2) === 'on' && typeof v === 'function') el.addEventListener(k.slice(2), v);
        else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
        else el.setAttribute(k, v === true ? '' : String(v));
      }
    }
    for (var i = 2; i < arguments.length; i++) {
      var c = arguments[i];
      if (c === null || c === undefined || c === false) continue;
      if (Array.isArray(c)) c.forEach(function (x) { if (x) el.appendChild(typeof x === 'string' ? doc.createTextNode(x) : x); });
      else el.appendChild(typeof c === 'string' ? doc.createTextNode(c) : c);
    }
    return el;
  }
  function esc(s) { return String(s).replace(/[&<>"']/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]; }); }
  function clamp(x, a, b) { return x < a ? a : x > b ? b : x; }
  function lerp(a, b, t) { return a + (b - a) * t; }
  function isFn(f) { return typeof f === 'function'; }
  function safeStorage() { try { var s = global.localStorage; s.setItem('__viz', '1'); s.removeItem('__viz'); return s; } catch (e) { return null; } }
  var STORE = safeStorage();
  Viz.h = h; Viz.esc = esc; Viz.clamp = clamp; Viz.lerp = lerp;

  var SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
  /** Format a number for humans: 3 significant figures, ×10ⁿ outside [1e-3, 1e5), optional unit. */
  function fmt(v, opt) {
    opt = opt || {};
    var sig = opt.sig || 3, unit = opt.unit ? ' ' + opt.unit : '';
    if (v === null || v === undefined || (typeof v === 'number' && isNaN(v))) return '—';
    if (typeof v !== 'number') return String(v) + unit;
    if (!isFinite(v)) return (v > 0 ? '∞' : '−∞') + unit;
    if (v === 0 || (Math.abs(v) < 1e-12 && !opt.keepTiny)) return '0' + unit;   // round-off (sin π) is not information
    var a = Math.abs(v), s;
    if (a >= 1e-3 && a < 1e5) {
      s = Number(v.toPrecision(sig)).toString();
      if (Math.abs(Number(s)) >= Math.pow(10, sig)) s = Math.round(v).toString();
    } else {
      var e = Math.floor(Math.log10(a)), m = v / Math.pow(10, e);
      if (Math.abs(Number(m.toPrecision(sig))) >= 10) { m /= 10; e += 1; }
      s = Number(m.toPrecision(sig)).toString() + '×10' + String(e).split('').map(function (c) { return SUP[c]; }).join('');
    }
    return (opt.plus && v > 0 ? '+' : '') + s.replace('-', '−') + unit;
  }
  /** Same, as a TeX string (for live equations): 1.23\times10^{-4}. */
  function tnum(v, sig) {
    sig = sig || 3;
    if (typeof v !== 'number' || !isFinite(v)) return '\\text{' + fmt(v) + '}';
    if (v === 0 || Math.abs(v) < 1e-12) return '0';
    var a = Math.abs(v);
    if (a >= 1e-3 && a < 1e5) return Number(v.toPrecision(sig)).toString();
    var e = Math.floor(Math.log10(a)), m = v / Math.pow(10, e);
    return Number(m.toPrecision(sig)).toString() + '\\times10^{' + e + '}';
  }
  Viz.fmt = fmt; Viz.tnum = tnum;

  /* ------------------------------------------------------------------------------------------------
     1. numerics (mirror numpy/scipy names so the parity with fluidpy is easy to read)
     ------------------------------------------------------------------------------------------------ */
  var num = {};
  num.linspace = function (a, b, n) { var out = new Array(n); if (n === 1) { out[0] = a; return out; } for (var i = 0; i < n; i++) out[i] = a + (b - a) * i / (n - 1); return out; };
  num.arange = function (a, b, step) { var out = []; for (var x = a; step > 0 ? x < b - 1e-12 : x > b + 1e-12; x += step) out.push(x); return out; };
  num.sum = function (arr) { var s = 0; for (var i = 0; i < arr.length; i++) s += arr[i]; return s; };
  num.max = function (arr) { var m = -Infinity; for (var i = 0; i < arr.length; i++) if (arr[i] > m) m = arr[i]; return m; };
  num.min = function (arr) { var m = Infinity; for (var i = 0; i < arr.length; i++) if (arr[i] < m) m = arr[i]; return m; };
  num.trapz = function (y, x) { var s = 0; for (var i = 1; i < y.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); return s; };
  /** One classical RK4 step for y' = f(t, y); y is an array. */
  num.rk4Step = function (f, t, y, dt) {
    var n = y.length, k1 = f(t, y), y2 = new Array(n), y3 = new Array(n), y4 = new Array(n), i;
    for (i = 0; i < n; i++) y2[i] = y[i] + 0.5 * dt * k1[i];
    var k2 = f(t + 0.5 * dt, y2);
    for (i = 0; i < n; i++) y3[i] = y[i] + 0.5 * dt * k2[i];
    var k3 = f(t + 0.5 * dt, y3);
    for (i = 0; i < n; i++) y4[i] = y[i] + dt * k3[i];
    var k4 = f(t + dt, y4), out = new Array(n);
    for (i = 0; i < n; i++) out[i] = y[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    return out;
  };
  /** Integrate y' = f(t, y) with fixed-step RK4 over the times ts (substeps per interval). Returns array of states. */
  num.odeint = function (f, y0, ts, substeps) {
    substeps = substeps || 4;
    var ys = [y0.slice()], y = y0.slice();
    for (var i = 1; i < ts.length; i++) {
      var dt = (ts[i] - ts[i - 1]) / substeps, t = ts[i - 1];
      for (var j = 0; j < substeps; j++) { y = num.rk4Step(f, t, y, dt); t += dt; }
      ys.push(y.slice());
    }
    return ys;
  };
  /** Brent-style bracketing root finder (bisection + secant safeguard). Throws if not bracketed. */
  num.brentq = function (f, a, b, tol, maxit) {
    tol = tol || 1e-12; maxit = maxit || 200;
    var fa = f(a), fb = f(b);
    if (fa === 0) return a; if (fb === 0) return b;
    if (fa * fb > 0) throw new Error('brentq: root not bracketed on [' + a + ', ' + b + ']');
    var c = a, fc = fa, d = b - a, e = d;
    for (var it = 0; it < maxit; it++) {
      if (fb * fc > 0) { c = a; fc = fa; d = b - a; e = d; }
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      var tol1 = 2 * 2.2e-16 * Math.abs(b) + 0.5 * tol, xm = 0.5 * (c - b);
      if (Math.abs(xm) <= tol1 || fb === 0) return b;
      if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
        var s = fb / fa, p, q, r;
        if (a === c) { p = 2 * xm * s; q = 1 - s; }
        else { q = fa / fc; r = fb / fc; p = s * (2 * xm * q * (q - r) - (b - a) * (r - 1)); q = (q - 1) * (r - 1) * (s - 1); }
        if (p > 0) q = -q; else p = -p;
        if (2 * p < Math.min(3 * xm * q - Math.abs(tol1 * q), Math.abs(e * q))) { e = d; d = p / q; }
        else { d = xm; e = d; }
      } else { d = xm; e = d; }
      a = b; fa = fb;
      b += Math.abs(d) > tol1 ? d : (xm > 0 ? tol1 : -tol1);
      fb = f(b);
    }
    return b;
  };
  /** erf with |error| < 1.2e-7 (Numerical Recipes erfc Chebyshev fit). Good for display; parity tests use rtol ≥ 1e-6. */
  num.erfc = function (x) {
    var z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    var r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 +
      t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  };
  num.erf = function (x) { return 1 - num.erfc(x); };
  /** Complex numbers as [re, im] — for potential flow and conformal maps. */
  num.C = {
    add: function (a, b) { return [a[0] + b[0], a[1] + b[1]]; },
    sub: function (a, b) { return [a[0] - b[0], a[1] - b[1]]; },
    mul: function (a, b) { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; },
    div: function (a, b) { var d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
    scale: function (a, s) { return [a[0] * s, a[1] * s]; },
    abs: function (a) { return Math.hypot(a[0], a[1]); },
    arg: function (a) { return Math.atan2(a[1], a[0]); },
    exp: function (a) { var e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; },
    log: function (a) { return [Math.log(Math.hypot(a[0], a[1])), Math.atan2(a[1], a[0])]; },
    pow: function (a, p) { var r = Math.pow(Math.hypot(a[0], a[1]), p), t = Math.atan2(a[1], a[0]) * p; return [r * Math.cos(t), r * Math.sin(t)]; },
    conj: function (a) { return [a[0], -a[1]]; }
  };
  /** "Nice" axis ticks between a and b (about n ticks). */
  num.niceTicks = function (a, b, n) {
    n = n || 5;
    if (!(isFinite(a) && isFinite(b)) || a === b) return [a];
    var lo = Math.min(a, b), hi = Math.max(a, b), span = hi - lo, raw = span / n;
    var mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
    var step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag;
    var out = [], start = Math.ceil(lo / step - 1e-9) * step;
    for (var x = start; x <= hi + step * 1e-9; x += step) out.push(Math.abs(x) < step * 1e-9 ? 0 : x);
    return out;
  };
  Viz.num = num;

  /* ------------------------------------------------------------------------------------------------
     2. colours and colormaps
     ------------------------------------------------------------------------------------------------ */
  var CMAPS = {
    viridis: ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'],
    coolwarm: ['#3b4cc0', '#6282ea', '#8db0fe', '#b8d0f9', '#dddcdc', '#f5c4ac', '#f39b7a', '#dc5d4a', '#b40426'],
    blues: ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b'],
    ocean: ['#0b1d3a', '#12406b', '#176b95', '#2f97b3', '#6cc0c9', '#b4e3dc', '#effaf5']
  };
  function hex2rgb(hx) { var n = parseInt(hx.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
  var CMAP_RGB = {};
  /** Colormap lookup: returns [r,g,b] (0..255) for t in [0,1]. */
  function cmap(name, t) {
    var stops = CMAP_RGB[name] || (CMAP_RGB[name] = (CMAPS[name] || CMAPS.viridis).map(hex2rgb));
    t = clamp(isFinite(t) ? t : 0, 0, 1) * (stops.length - 1);
    var i = Math.min(Math.floor(t), stops.length - 2), f = t - i, a = stops[i], b = stops[i + 1];
    return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
  }
  Viz.cmap = cmap; Viz.CMAPS = CMAPS;
  /** Read a design token from viz_base.css, e.g. Viz.color('accent') → '#6c5ce7' (follows the theme). */
  Viz.color = function (name) {
    var v = getComputedStyle(doc.documentElement).getPropertyValue('--viz-' + name).trim();
    return v || '#6c5ce7';
  };
  Viz.alpha = function (color, a) {
    if (color.charAt(0) === '#') { var c = hex2rgb(color.length === 4 ? '#' + color[1] + color[1] + color[2] + color[2] + color[3] + color[3] : color); return 'rgba(' + c[0] + ',' + c[1] + ',' + c[2] + ',' + a + ')'; }
    return color;
  };

  /* ------------------------------------------------------------------------------------------------
     3. KaTeX (CDN, two mirrors, text fallback) and inline $…$ in HTML strings
     ------------------------------------------------------------------------------------------------ */
  var KATEX_VER = '0.16.11';
  var KATEX_SRC = [
    'https://cdn.jsdelivr.net/npm/katex@' + KATEX_VER + '/dist/',
    'https://cdnjs.cloudflare.com/ajax/libs/KaTeX/' + KATEX_VER + '/'
  ];
  var katexPromise = null;
  function loadKatex() {
    if (katexPromise) return katexPromise;
    katexPromise = new Promise(function (resolve) {
      if (global.katex) { resolve(true); return; }
      var i = 0, done = false;
      var timer = setTimeout(function () { if (!done) { done = true; resolve(false); } }, 6000);
      function attempt() {
        if (i >= KATEX_SRC.length) { if (!done) { done = true; clearTimeout(timer); resolve(false); } return; }
        var base = KATEX_SRC[i++];
        var link = h('link', { rel: 'stylesheet', href: base + 'katex.min.css' });
        doc.head.appendChild(link);
        var s = h('script', { src: base + 'katex.min.js' });
        s.onload = function () { if (!done && global.katex) { done = true; clearTimeout(timer); resolve(true); } };
        s.onerror = function () { link.remove(); s.remove(); attempt(); };
        doc.head.appendChild(s);
      }
      attempt();
    });
    return katexPromise;
  }
  function renderTex(el, tex, display) {
    if (global.katex) {
      try { global.katex.render(tex, el, { displayMode: !!display, throwOnError: false, strict: 'ignore', trust: false }); return; }
      catch (e) { /* fall through */ }
    }
    el.innerHTML = '<span class="viz-tex-fallback">' + esc(tex) + '</span>';
    el.setAttribute('data-tex-pending', '1');
    el._vizTex = [tex, display];
  }
  /** Turn "$a+b$" and "$$…$$" inside an HTML string into math spans (rendered by Viz.typeset). */
  function mathify(html) {
    if (!html) return '';
    return String(html)
      .replace(/\$\$([\s\S]+?)\$\$/g, function (_, t) { return '<span class="viz-tex" data-display="1" data-tex="' + esc(t) + '"></span>'; })
      .replace(/\$([^$\n]+?)\$/g, function (_, t) { return '<span class="viz-tex" data-tex="' + esc(t) + '"></span>'; });
  }
  function typeset(root) {
    var spans = (root || doc).querySelectorAll('.viz-tex[data-tex]');
    for (var i = 0; i < spans.length; i++) renderTex(spans[i], spans[i].getAttribute('data-tex'), spans[i].getAttribute('data-display') === '1');
  }
  function retypesetPending() {
    var els = doc.querySelectorAll('[data-tex-pending]');
    for (var i = 0; i < els.length; i++) { var p = els[i]._vizTex; els[i].removeAttribute('data-tex-pending'); if (p) renderTex(els[i], p[0], p[1]); }
  }
  Viz.loadKatex = loadKatex; Viz.renderTex = renderTex; Viz.mathify = mathify; Viz.typeset = typeset;

  /* ------------------------------------------------------------------------------------------------
     4. canvas + 2-D plotting in world coordinates
     ------------------------------------------------------------------------------------------------ */
  /** Axis-aligned plot in a pixel rectangle of the stage canvas. All drawing calls take WORLD coordinates. */
  function Plot(g, opt) {
    opt = opt || {};
    var ctx = g.ctx, pad = Object.assign({ l: 48, r: 12, t: 12, b: 38 }, opt.pad || {});
    var rect = opt.rect || { x: 0, y: 0, w: g.w, h: g.h };
    if (g.w < 420) { pad.l = Math.min(pad.l, 40); pad.b = Math.min(pad.b, 34); }
    var px0 = rect.x + pad.l, py0 = rect.y + pad.t, pw = Math.max(10, rect.w - pad.l - pad.r), ph = Math.max(10, rect.h - pad.t - pad.b);
    var xl = (opt.xlim || [0, 1]).slice(), yl = (opt.ylim || [0, 1]).slice();
    if (opt.equal) { // same scale on both axes: expand the shorter range
      var sx = pw / (xl[1] - xl[0]), sy = ph / (yl[1] - yl[0]), s = Math.min(sx, sy);
      var cx = (xl[0] + xl[1]) / 2, cy = (yl[0] + yl[1]) / 2;
      xl = [cx - pw / s / 2, cx + pw / s / 2]; yl = [cy - ph / s / 2, cy + ph / s / 2];
    }
    var P = {
      g: g, ctx: ctx, xlim: xl, ylim: yl, px: { x: px0, y: py0, w: pw, h: ph }, rect: rect,
      X: function (x) { return px0 + (x - xl[0]) / (xl[1] - xl[0]) * pw; },
      Y: function (y) { return py0 + ph - (y - yl[0]) / (yl[1] - yl[0]) * ph; },
      ix: function (px) { return xl[0] + (px - px0) / pw * (xl[1] - xl[0]); },
      iy: function (py) { return yl[0] + (py0 + ph - py) / ph * (yl[1] - yl[0]); },
      sx: function (dx) { return dx / (xl[1] - xl[0]) * pw; },
      sy: function (dy) { return dy / (yl[1] - yl[0]) * ph; },
      inside: function (px, py) { return px >= px0 && px <= px0 + pw && py >= py0 && py <= py0 + ph; }
    };
    P.clip = function (fn) { ctx.save(); ctx.beginPath(); ctx.rect(px0, py0, pw, ph); ctx.clip(); try { fn(P); } finally { ctx.restore(); } };
    P.frame = function () {
      ctx.save(); ctx.strokeStyle = Viz.color('border'); ctx.lineWidth = 1; ctx.strokeRect(px0 + 0.5, py0 + 0.5, pw - 1, ph - 1); ctx.restore();
    };
    P.axes = function (o) {
      o = Object.assign({ grid: true, xlabel: opt.xlabel, ylabel: opt.ylabel, xticks: null, yticks: null, zeroLines: true }, o || {});
      var fs = g.w < 420 ? 11 : 12;
      ctx.save();
      ctx.font = fs + 'px ' + getComputedStyle(doc.body).fontFamily;
      var xt = o.xticks || num.niceTicks(xl[0], xl[1], Math.max(2, Math.round(pw / 90)));
      var yt = o.yticks || num.niceTicks(yl[0], yl[1], Math.max(2, Math.round(ph / 60)));
      if (o.grid) {
        ctx.strokeStyle = Viz.color('grid'); ctx.lineWidth = 1; ctx.beginPath();
        xt.forEach(function (x) { var X = Math.round(P.X(x)) + 0.5; ctx.moveTo(X, py0); ctx.lineTo(X, py0 + ph); });
        yt.forEach(function (y) { var Y = Math.round(P.Y(y)) + 0.5; ctx.moveTo(px0, Y); ctx.lineTo(px0 + pw, Y); });
        ctx.stroke();
      }
      if (o.zeroLines) {
        ctx.strokeStyle = Viz.alpha(Viz.color('axis'), 0.55); ctx.beginPath();
        if (xl[0] < 0 && xl[1] > 0) { var X0 = Math.round(P.X(0)) + 0.5; ctx.moveTo(X0, py0); ctx.lineTo(X0, py0 + ph); }
        if (yl[0] < 0 && yl[1] > 0) { var Y0 = Math.round(P.Y(0)) + 0.5; ctx.moveTo(px0, Y0); ctx.lineTo(px0 + pw, Y0); }
        ctx.stroke();
      }
      P.frame();
      ctx.fillStyle = Viz.color('axis');
      ctx.textAlign = 'center'; ctx.textBaseline = 'top';
      xt.forEach(function (x) { ctx.fillText(fmt(x, { sig: 3 }), P.X(x), py0 + ph + 4); });
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      yt.forEach(function (y) { ctx.fillText(fmt(y, { sig: 3 }), px0 - 5, P.Y(y)); });
      ctx.fillStyle = Viz.color('text-soft');
      if (o.xlabel) { ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(o.xlabel, px0 + pw / 2, rect.y + rect.h - 1); }
      if (o.ylabel) { ctx.save(); ctx.translate(rect.x + 11, py0 + ph / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(o.ylabel, 0, 0); ctx.restore(); }
      ctx.restore();
      return P;
    };
    function style(o) {
      o = o || {};
      ctx.strokeStyle = o.color || Viz.color('accent'); ctx.fillStyle = o.fill || o.color || Viz.color('accent');
      ctx.lineWidth = o.width || 2; ctx.setLineDash(o.dash || []); ctx.globalAlpha = o.alpha === undefined ? 1 : o.alpha;
      ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    }
    P.line = function (xs, ys, o) {
      ctx.save(); style(o); ctx.beginPath(); var started = false;
      for (var i = 0; i < xs.length; i++) {
        var X = P.X(xs[i]), Y = P.Y(ys[i]);
        if (!isFinite(X) || !isFinite(Y)) { started = false; continue; }
        if (!started) { ctx.moveTo(X, Y); started = true; } else ctx.lineTo(X, Y);
      }
      ctx.stroke(); ctx.restore(); return P;
    };
    P.poly = function (pts, o) { return P.line(pts.map(function (p) { return p[0]; }), pts.map(function (p) { return p[1]; }), o); };
    P.fn = function (f, o) { o = o || {}; var n = o.n || Math.max(64, Math.round(pw)); var xs = num.linspace(o.x0 === undefined ? xl[0] : o.x0, o.x1 === undefined ? xl[1] : o.x1, n); return P.line(xs, xs.map(f), o); };
    P.fillBetween = function (xs, y1, y2, o) {
      ctx.save(); style(o); ctx.globalAlpha = o && o.alpha !== undefined ? o.alpha : 0.2; ctx.beginPath();
      for (var i = 0; i < xs.length; i++) { var Y1 = P.Y(Array.isArray(y1) ? y1[i] : y1); if (i === 0) ctx.moveTo(P.X(xs[i]), Y1); else ctx.lineTo(P.X(xs[i]), Y1); }
      for (var j = xs.length - 1; j >= 0; j--) ctx.lineTo(P.X(xs[j]), P.Y(Array.isArray(y2) ? y2[j] : y2));
      ctx.closePath(); ctx.fill(); ctx.restore(); return P;
    };
    P.dot = function (x, y, r, o) { ctx.save(); style(o); ctx.beginPath(); ctx.arc(P.X(x), P.Y(y), r || 4, 0, 2 * Math.PI); ctx.fill(); if (o && o.stroke) { ctx.strokeStyle = o.stroke; ctx.lineWidth = o.strokeWidth || 1.5; ctx.stroke(); } ctx.restore(); return P; };
    P.rectW = function (x, y, w, hh, o) { ctx.save(); style(o); var X = P.X(x), Y = P.Y(y + hh); if (o && o.fill) { ctx.fillRect(X, Y, P.sx(w), P.sy(hh)); } if (!o || o.stroke !== false) ctx.strokeRect(X, Y, P.sx(w), P.sy(hh)); ctx.restore(); return P; };
    /** Arrow from (x,y) by (dx,dy) in world units; head size in px. */
    P.arrow = function (x, y, dx, dy, o) {
      o = o || {}; var X0 = P.X(x), Y0 = P.Y(y), X1 = P.X(x + dx), Y1 = P.Y(y + dy);
      var L = Math.hypot(X1 - X0, Y1 - Y0); if (L < 0.5) return P;
      var head = Math.min(o.head || 8, L * 0.6), ang = Math.atan2(Y1 - Y0, X1 - X0);
      ctx.save(); style(o); ctx.beginPath(); ctx.moveTo(X0, Y0); ctx.lineTo(X1 - head * 0.7 * Math.cos(ang), Y1 - head * 0.7 * Math.sin(ang)); ctx.stroke();
      ctx.setLineDash([]); ctx.beginPath(); ctx.moveTo(X1, Y1);
      ctx.lineTo(X1 - head * Math.cos(ang - 0.42), Y1 - head * Math.sin(ang - 0.42));
      ctx.lineTo(X1 - head * Math.cos(ang + 0.42), Y1 - head * Math.sin(ang + 0.42));
      ctx.closePath(); ctx.fillStyle = o.color || Viz.color('accent'); ctx.fill(); ctx.restore(); return P;
    };
    P.text = function (x, y, str, o) {
      o = o || {}; ctx.save();
      ctx.font = (o.weight || 500) + ' ' + (o.size || 12.5) + 'px ' + getComputedStyle(doc.body).fontFamily;
      ctx.fillStyle = o.color || Viz.color('text'); ctx.textAlign = o.align || 'left'; ctx.textBaseline = o.baseline || 'middle';
      var X = P.X(x) + (o.dx || 0), Y = P.Y(y) + (o.dy || 0);
      if (o.halo !== false) { ctx.lineWidth = 3; ctx.strokeStyle = Viz.alpha(Viz.color('card'), 0.85); ctx.strokeText(str, X, Y); }
      ctx.fillText(str, X, Y); ctx.restore(); return P;
    };
    return P;
  }
  Viz.Plot = Plot;

  /* ------------------------------------------------------------------------------------------------
     5. flow-field tools
     ------------------------------------------------------------------------------------------------ */
  var field = {};
  /** Sample f(x, y) on an nx × ny node grid spanning xlim × ylim. Returns {v: Float64Array (row-major, j*nx+i), nx, ny, xs, ys}. */
  field.grid = function (f, xlim, ylim, nx, ny) {
    var xs = num.linspace(xlim[0], xlim[1], nx), ys = num.linspace(ylim[0], ylim[1], ny), v = new Float64Array(nx * ny);
    for (var j = 0; j < ny; j++) for (var i = 0; i < nx; i++) v[j * nx + i] = f(xs[i], ys[j]);
    return { v: v, nx: nx, ny: ny, xs: xs, ys: ys };
  };
  /** Marching squares: line segments [[x1,y1,x2,y2], …] of the level set G = level. */
  field.contour = function (G, level) {
    var segs = [], nx = G.nx, ny = G.ny, v = G.v, xs = G.xs, ys = G.ys;
    function interp(x1, y1, v1, x2, y2, v2) { var t = (level - v1) / (v2 - v1); return [x1 + t * (x2 - x1), y1 + t * (y2 - y1)]; }
    for (var j = 0; j < ny - 1; j++) for (var i = 0; i < nx - 1; i++) {
      var a = v[j * nx + i], b = v[j * nx + i + 1], c = v[(j + 1) * nx + i + 1], d = v[(j + 1) * nx + i];
      if (!(isFinite(a) && isFinite(b) && isFinite(c) && isFinite(d))) continue;
      var idx = (a > level ? 1 : 0) | (b > level ? 2 : 0) | (c > level ? 4 : 0) | (d > level ? 8 : 0);
      if (idx === 0 || idx === 15) continue;
      var x0 = xs[i], x1 = xs[i + 1], y0 = ys[j], y1 = ys[j + 1], e = [];
      var bottom = function () { return interp(x0, y0, a, x1, y0, b); }, right = function () { return interp(x1, y0, b, x1, y1, c); };
      var top = function () { return interp(x0, y1, d, x1, y1, c); }, left = function () { return interp(x0, y0, a, x0, y1, d); };
      switch (idx) {
        case 1: case 14: e = [bottom(), left()]; break;
        case 2: case 13: e = [bottom(), right()]; break;
        case 3: case 12: e = [left(), right()]; break;
        case 4: case 11: e = [right(), top()]; break;
        case 6: case 9: e = [bottom(), top()]; break;
        case 7: case 8: e = [left(), top()]; break;
        case 5: e = [bottom(), right(), left(), top()]; break;   // saddle: pick one resolution
        case 10: e = [bottom(), left(), right(), top()]; break;
      }
      for (var k = 0; k + 1 < e.length; k += 2) segs.push([e[k][0], e[k][1], e[k + 1][0], e[k + 1][1]]);
    }
    return segs;
  };
  /** Draw contour lines of G at the given levels. o: {color | colorFn(level), width, dash}. */
  field.drawContours = function (P, G, levels, o) {
    o = o || {}; var ctx = P.ctx;
    P.clip(function () {
      levels.forEach(function (lev, n) {
        var segs = field.contour(G, lev);
        ctx.save(); ctx.strokeStyle = o.colorFn ? o.colorFn(lev, n) : (o.color || Viz.color('accent'));
        ctx.lineWidth = o.width || 1.4; ctx.setLineDash(o.dash || []); ctx.beginPath();
        for (var s = 0; s < segs.length; s++) { ctx.moveTo(P.X(segs[s][0]), P.Y(segs[s][1])); ctx.lineTo(P.X(segs[s][2]), P.Y(segs[s][3])); }
        ctx.stroke(); ctx.restore();
      });
    });
  };
  /** Colour image of f(x,y) (or a grid G) behind a plot. o: {cmap, vmin, vmax, nx, ny, alpha}. */
  field.heatmap = function (P, f, o) {
    o = o || {}; var nx = o.nx || 120, ny = o.ny || Math.max(40, Math.round(nx * P.px.h / P.px.w));
    var G = f.v ? f : field.grid(f, P.xlim, P.ylim, nx, ny);
    var vmin = o.vmin, vmax = o.vmax;
    if (vmin === undefined || vmax === undefined) { var lo = Infinity, hi = -Infinity; for (var q = 0; q < G.v.length; q++) { if (isFinite(G.v[q])) { lo = Math.min(lo, G.v[q]); hi = Math.max(hi, G.v[q]); } } if (vmin === undefined) vmin = lo; if (vmax === undefined) vmax = hi; }
    var off = doc.createElement('canvas'); off.width = G.nx; off.height = G.ny;
    var octx = off.getContext('2d'), img = octx.createImageData(G.nx, G.ny);
    for (var j = 0; j < G.ny; j++) for (var i = 0; i < G.nx; i++) {
      var val = G.v[j * G.nx + i], p = ((G.ny - 1 - j) * G.nx + i) * 4, c = cmap(o.cmap || 'viridis', (val - vmin) / (vmax - vmin || 1));
      img.data[p] = c[0]; img.data[p + 1] = c[1]; img.data[p + 2] = c[2]; img.data[p + 3] = isFinite(val) ? 255 : 0;
    }
    octx.putImageData(img, 0, 0);
    var ctx = P.ctx; ctx.save(); ctx.globalAlpha = o.alpha === undefined ? 1 : o.alpha; ctx.imageSmoothingEnabled = true;
    ctx.drawImage(off, P.X(G.xs[0]), P.Y(G.ys[G.ny - 1]), P.X(G.xs[G.nx - 1]) - P.X(G.xs[0]), P.Y(G.ys[0]) - P.Y(G.ys[G.ny - 1]));
    ctx.restore();
    return { vmin: vmin, vmax: vmax };
  };
  /** Streamline through (x0, y0) of the steady field vel(x, y) → [u, v]; RK4 in arc length, both directions. */
  field.streamline = function (vel, x0, y0, o) {
    o = o || {}; var ds = o.ds || 0.02, n = o.n || 600, b = o.bounds, pts = [[x0, y0]];
    function unit(t, p) { var w = vel(p[0], p[1]), s = Math.hypot(w[0], w[1]); return s < 1e-12 ? [0, 0] : [w[0] / s, w[1] / s]; }
    function run(sign) {
      var p = [x0, y0], out = [];
      for (var k = 0; k < n; k++) {
        p = num.rk4Step(function (t, q) { var u = unit(t, q); return [sign * u[0], sign * u[1]]; }, 0, p, ds);
        if (!isFinite(p[0]) || !isFinite(p[1])) break;
        if (b && (p[0] < b[0] || p[0] > b[1] || p[1] < b[2] || p[1] > b[3])) break;
        if (o.stop && o.stop(p[0], p[1])) break;
        out.push(p.slice());
      }
      return out;
    }
    var fwd = run(1); if (o.both === false) return pts.concat(fwd);
    return run(-1).reverse().concat(pts, fwd);
  };
  /** Arrow grid of vel(x, y). o: {nx, ny, scale (world units per unit speed; auto), color, maxLen px}. */
  field.quiver = function (P, vel, o) {
    o = o || {}; var nx = o.nx || 14, ny = o.ny || Math.max(4, Math.round(nx * P.px.h / P.px.w));
    var xs = num.linspace(P.xlim[0], P.xlim[1], nx + 2).slice(1, -1), ys = num.linspace(P.ylim[0], P.ylim[1], ny + 2).slice(1, -1);
    var smax = 0, W = [];
    ys.forEach(function (y) { xs.forEach(function (x) { var w = vel(x, y); W.push([x, y, w[0], w[1]]); var s = Math.hypot(w[0], w[1]); if (isFinite(s)) smax = Math.max(smax, s); }); });
    var cell = Math.min((P.xlim[1] - P.xlim[0]) / (nx + 1), (P.ylim[1] - P.ylim[0]) / (ny + 1));
    var scale = o.scale || (smax > 0 ? 0.85 * cell / smax : 1);
    P.clip(function () { W.forEach(function (a) { if (isFinite(a[2]) && isFinite(a[3])) P.arrow(a[0], a[1], a[2] * scale, a[3] * scale, { color: o.color || Viz.color('muted'), width: 1.3, head: 6 }); }); });
    return scale;
  };
  /** Passive tracer particles advected by vel(x, y, t). o: {n, bounds:[x0,x1,y0,y1], spawn(i) → [x,y], life (s), trail}. */
  field.particles = function (vel, o) {
    o = o || {}; var n = o.n || 200, b = o.bounds || [0, 1, 0, 1], trail = o.trail || 0;
    var rng = Viz.rng(o.seed || 1);
    function spawn(i) { return o.spawn ? o.spawn(i, rng) : [lerp(b[0], b[1], rng()), lerp(b[2], b[3], rng())]; }
    var P = [];
    for (var i = 0; i < n; i++) { var s = spawn(i); P.push({ x: s[0], y: s[1], age: rng() * (o.life || 4), hist: [] }); }
    return {
      list: P,
      step: function (dt, t) {
        for (var i = 0; i < P.length; i++) {
          var p = P[i];
          var y = num.rk4Step(function (tt, q) { return vel(q[0], q[1], tt); }, t, [p.x, p.y], dt);
          if (trail) { p.hist.push([p.x, p.y]); if (p.hist.length > trail) p.hist.shift(); }
          p.x = y[0]; p.y = y[1]; p.age += dt;
          if (!isFinite(p.x) || !isFinite(p.y) || p.x < b[0] || p.x > b[1] || p.y < b[2] || p.y > b[3] || (o.life && p.age > o.life) || (o.kill && o.kill(p.x, p.y))) {
            var s = spawn(i); p.x = s[0]; p.y = s[1]; p.age = 0; p.hist = [];
          }
        }
      },
      draw: function (Pl, st) {
        st = st || {}; var ctx = Pl.ctx;
        Pl.clip(function () {
          ctx.save(); ctx.fillStyle = st.color || Viz.color('blue'); ctx.strokeStyle = Viz.alpha(st.color || Viz.color('blue'), 0.35); ctx.lineWidth = 1;
          for (var i = 0; i < P.length; i++) {
            var p = P[i];
            if (trail && p.hist.length > 1) { ctx.beginPath(); ctx.moveTo(Pl.X(p.hist[0][0]), Pl.Y(p.hist[0][1])); for (var k = 1; k < p.hist.length; k++) ctx.lineTo(Pl.X(p.hist[k][0]), Pl.Y(p.hist[k][1])); ctx.lineTo(Pl.X(p.x), Pl.Y(p.y)); ctx.stroke(); }
            ctx.beginPath(); ctx.arc(Pl.X(p.x), Pl.Y(p.y), st.r || 1.8, 0, 2 * Math.PI); ctx.fill();
          }
          ctx.restore();
        });
      }
    };
  };
  Viz.field = field;
  /** Seeded PRNG (mulberry32) so every explainer looks the same on every load and in screenshots. */
  Viz.rng = function (seed) { var a = seed >>> 0; return function () { a |= 0; a = a + 0x6D2B79F5 | 0; var t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; };

  /* ------------------------------------------------------------------------------------------------
     6. icons
     ------------------------------------------------------------------------------------------------ */
  var ICON = {
    help: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="9.5"/><path d="M9.3 9.2a2.8 2.8 0 0 1 5.4 1c0 1.9-2.7 2.4-2.7 4.1"/><circle cx="12" cy="17.6" r=".6" fill="currentColor"/></svg>',
    theme: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M20 14.5A8 8 0 1 1 9.5 4a6.5 6.5 0 0 0 10.5 10.5z"/></svg>',
    full: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/></svg>',
    reset: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 4v6h6"/><path d="M5.5 15a7.5 7.5 0 1 0 1.8-7.8L4 10"/></svg>',
    play: '<svg viewBox="0 0 24 24" fill="currentColor"><path d="M7 4.5v15l13-7.5z"/></svg>',
    pause: '<svg viewBox="0 0 24 24" fill="currentColor"><rect x="6" y="4.5" width="4" height="15" rx="1"/><rect x="14" y="4.5" width="4" height="15" rx="1"/></svg>'
  };
  Viz.ICON = ICON;

  /* ------------------------------------------------------------------------------------------------
     6b. Python syntax colouring for the Code tab ({{name}} placeholders become live values)
     ------------------------------------------------------------------------------------------------ */
  var PY_KW = /^(def|return|import|from|as|for|in|if|elif|else|while|and|or|not|None|True|False|lambda|class|with|yield|try|except|raise|pass|break|continue|is|global|assert)$/;
  function liveSpan(name) { return '<span class="viz-code-live" data-code-live="' + esc(name) + '">…</span>'; }
  function withLive(escaped) { return escaped.replace(/\{\{\s*([\w.]+)\s*\}\}/g, function (_, n) { return liveSpan(n); }); }
  function highlightPy(line) {
    var out = '', i = 0, n = line.length, m;
    while (i < n) {
      var ch = line.charAt(i), rest = line.slice(i);
      if (rest.slice(0, 2) === '{{' && (m = /^\{\{\s*([\w.]+)\s*\}\}/.exec(rest))) { out += liveSpan(m[1]); i += m[0].length; continue; }
      if (ch === '#') { out += '<span class="cm">' + withLive(esc(rest)) + '</span>'; break; }
      if (ch === '"' || ch === "'") {
        var q = rest.slice(0, 3) === ch + ch + ch ? ch + ch + ch : ch, j = rest.indexOf(q, q.length);
        var s = j < 0 ? rest : rest.slice(0, j + q.length);
        out += '<span class="st">' + withLive(esc(s)) + '</span>'; i += s.length; continue;
      }
      if ((m = /^\d*\.?\d+(?:[eE][-+]?\d+)?/.exec(rest)) && !/[A-Za-z_]/.test(line.charAt(i - 1) || '')) { out += '<span class="num">' + m[0] + '</span>'; i += m[0].length; continue; }
      if ((m = /^[A-Za-z_]\w*/.exec(rest))) {
        var word = m[0], after = line.slice(i + word.length).replace(/^\s+/, '');
        if (PY_KW.test(word)) out += '<span class="kw">' + word + '</span>';
        else if (after.charAt(0) === '(') out += '<span class="fn">' + word + '</span>';
        else out += word;
        i += word.length; continue;
      }
      out += esc(ch); i++;
    }
    return out || ' ';
  }
  Viz.highlightPy = highlightPy;

  /** Builders for the "Explain" tab: every line = formula = numbers substituted = result, with a short why.
      Put values that change with time in live spans (Viz.live('x')) and return them from calc.live(s). */
  Viz.work = {
    head: function (text) { return '<div class="viz-work-h">' + text + '</div>'; },
    line: function (formula, value, why) { return '<div class="viz-work-ln">' + formula + (value !== undefined && value !== '' ? ' = <b>' + value + '</b>' : '') + (why ? '<span class="viz-work-why">' + why + '</span>' : '') + '</div>'; },
    result: function (html, tone) { return '<div class="viz-work-res' + (tone ? ' ' + tone : '') + '">' + html + '</div>'; },
    table: function (headers, rows, current) {
      return '<table class="viz-work-table"><tr>' + headers.map(function (x) { return '<th>' + x + '</th>'; }).join('') + '</tr>' +
        rows.map(function (r, i) { return '<tr' + (i === current ? ' class="cur"' : '') + '>' + r.map(function (c) { return '<td>' + c + '</td>'; }).join('') + '</tr>'; }).join('') + '</table>';
    },
    note: function (html) { return '<p>' + html + '</p>'; },
    /* the "Explanation & interpretation" pattern of the reference mathlets (forced_damped_vibrations.html):
       numbered sections that compute every displayed quantity with the reader's numbers, then say what it means */
    step: function (n, title) { return '<div class="viz-work-h"><span class="viz-work-n">' + n + '</span>' + title + '</div>'; },
    say: function (html) { return '<p class="viz-work-say">' + html + '</p>'; },
    box: function (html) { return '<div class="viz-work-res teal">' + html + '</div>'; },
    interpret: function (html, label) { return '<div class="viz-work-interp"><b class="k">' + (label || 'What this means') + '</b> ' + html + '</div>'; },
    hint: function (html) { return '<p class="viz-work-hint">' + html + '</p>'; }
  };
  Viz.live = function (name) { return '<b data-live="' + esc(name) + '">…</b>'; };

  /* ------------------------------------------------------------------------------------------------
     6c. 3-D views with three.js (loaded on demand from a CDN; drag = orbit, wheel = zoom, click = pick)
     ------------------------------------------------------------------------------------------------ */
  var threePromise = null;
  Viz.loadThree = function () {
    if (threePromise) return threePromise;
    threePromise = new Promise(function (resolve) {
      if (global.THREE) { resolve(global.THREE); return; }
      var srcs = ['https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js', 'https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js'];
      var i = 0, done = false, timer = setTimeout(function () { if (!done) { done = true; resolve(null); } }, 10000);
      (function attempt() {
        if (i >= srcs.length) { if (!done) { done = true; clearTimeout(timer); resolve(null); } return; }
        var s = h('script', { src: srcs[i++] });
        s.onload = function () { if (!done) { done = true; clearTimeout(timer); resolve(global.THREE || null); } };
        s.onerror = function () { s.remove(); attempt(); };
        doc.head.appendChild(s);
      })();
    });
    return threePromise;
  };
  /** Turn a stage view into a 3-D scene. Resolves to {THREE, scene, camera, renderer, render(), pick(ev, objs), project(v3)}
      or null (offline: a friendly fallback note is shown). The view's 2-D canvas stays on top for labels. */
  Viz.three = function (view, opt) {
    opt = opt || {};
    return Viz.loadThree().then(function (THREE) {
      var holder = view.holder;
      if (!THREE) { holder.appendChild(h('div', { class: 'viz-3d-fallback', html: opt.fallback || 'The 3-D view needs an internet connection (three.js loads from a CDN). Everything else in this explainer still works.' })); return null; }
      var canvas = h('canvas', { class: 'viz-3d-canvas' });
      holder.insertBefore(canvas, holder.firstChild);
      view.canvas.style.pointerEvents = 'none';
      var renderer = new THREE.WebGLRenderer({ canvas: canvas, antialias: true, alpha: true });
      renderer.setPixelRatio(Math.min(global.devicePixelRatio || 1, 2));
      var scene = new THREE.Scene();
      var camera = new THREE.PerspectiveCamera(opt.fov || 40, 1, 0.01, 2000);
      var target = new THREE.Vector3().fromArray(opt.target || [0, 0, 0]);
      var orb = { r: opt.distance || 10, theta: opt.theta === undefined ? 0.8 : opt.theta, phi: opt.phi === undefined ? 1.1 : opt.phi };
      scene.add(new THREE.AmbientLight(0xffffff, 0.75));
      var dl = new THREE.DirectionalLight(0xffffff, 0.6); dl.position.set(5, 10, 7); scene.add(dl);
      function place() { var sp = Math.sin(orb.phi); camera.position.set(target.x + orb.r * sp * Math.cos(orb.theta), target.y + orb.r * Math.cos(orb.phi), target.z + orb.r * sp * Math.sin(orb.theta)); camera.lookAt(target); }
      function size() { var w = holder.clientWidth, hh = holder.clientHeight; if (w && hh) { renderer.setSize(w, hh, false); canvas.style.width = w + 'px'; canvas.style.height = hh + 'px'; camera.aspect = w / hh; camera.updateProjectionMatrix(); } }
      var drag = null, moved = 0, ray = new THREE.Raycaster(), mouse = new THREE.Vector2();
      var T = {
        THREE: THREE, scene: scene, camera: camera, renderer: renderer, orbit: orb, target: target,
        render: function () { size(); renderer.render(scene, camera); if (isFn(opt.onRender)) opt.onRender(T); },
        pick: function (e, objects) { var r = canvas.getBoundingClientRect(); mouse.x = ((e.clientX - r.left) / r.width) * 2 - 1; mouse.y = -((e.clientY - r.top) / r.height) * 2 + 1; ray.setFromCamera(mouse, camera); var hits = ray.intersectObjects(objects || scene.children, true); return hits[0] || null; },
        project: function (v3) { var p = v3.clone().project(camera); return { x: (p.x + 1) / 2 * holder.clientWidth, y: (1 - p.y) / 2 * holder.clientHeight }; }
      };
      canvas.addEventListener('pointerdown', function (e) { drag = { x: e.clientX, y: e.clientY, t: orb.theta, p: orb.phi }; moved = 0; try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ } });
      canvas.addEventListener('pointermove', function (e) { if (!drag) return; var dx = e.clientX - drag.x, dy = e.clientY - drag.y; moved = Math.max(moved, Math.abs(dx) + Math.abs(dy)); orb.theta = drag.t + dx * 0.008; orb.phi = clamp(drag.p - dy * 0.008, 0.08, Math.PI - 0.08); place(); T.render(); });
      canvas.addEventListener('pointerup', function (e) { var click = drag && moved < 5; drag = null; if (click && isFn(opt.onPick)) opt.onPick(T.pick(e), e); });
      canvas.addEventListener('wheel', function (e) { e.preventDefault(); orb.r = clamp(orb.r * (1 + Math.sign(e.deltaY) * 0.08), opt.minDistance || 1, opt.maxDistance || 200); place(); T.render(); }, { passive: false });
      canvas.style.touchAction = 'none';
      place();
      if ('ResizeObserver' in global) new ResizeObserver(function () { T.render(); }).observe(holder);
      T.render();
      return T;
    });
  };

  /* ------------------------------------------------------------------------------------------------
     7. the application shell
     ------------------------------------------------------------------------------------------------ */
  var DEFAULT_TABS = [
    { id: 'tour', label: 'Walkthrough', short: 'Guide' },
    { id: 'explore', label: 'Explore', short: 'Play' },
    { id: 'explain', label: 'Explain', short: 'Why' },
    { id: 'derive', label: 'Derivation', short: 'Derive' },
    { id: 'equations', label: 'Equations', short: 'Math' },
    { id: 'code', label: 'Code', short: 'Code' },
    { id: 'check', label: 'Check yourself', short: 'Quiz' }
  ];
  var CALLOUT_LABEL = { key: 'Key idea', watch: 'Watch', try: 'Try', warn: 'Careful' };

  function parseHash() {
    var out = {};
    try { (global.location.hash || '').replace(/^#/, '').split('&').forEach(function (kv) { if (!kv) return; var p = kv.split('='); out[decodeURIComponent(p[0])] = decodeURIComponent(p[1] || ''); }); } catch (e) { /* srcdoc */ }
    return out;
  }
  function writeHash(obj) {
    try {
      var s = Object.keys(obj).filter(function (k) { return obj[k] !== undefined && obj[k] !== null && obj[k] !== ''; })
        .map(function (k) { return encodeURIComponent(k) + '=' + encodeURIComponent(obj[k]); }).join('&');
      global.history.replaceState(null, '', '#' + s);
    } catch (e) { /* about:srcdoc and sandboxed frames refuse; harmless */ }
  }
  function now() { return (global.performance && global.performance.now) ? global.performance.now() : Date.now(); }

  Viz.app = function (cfg) {
    if (cfg.calc && !cfg.explain) cfg.explain = cfg.calc;          // `calc` is the old name of `explain`
    var app = { cfg: cfg, state: {}, defaults: {}, tab: null, step: 0, playing: false, holding: false, t: 0, controls: {}, _hl: [] };
    var mount = doc.getElementById(cfg.mount || 'app') || doc.body.appendChild(h('div', { id: 'app' }));
    mount.classList.add('viz-mount');
    var params = Object.assign({}, cfg.params || {});
    var tr = cfg.transport || null;
    if (tr) {
      if (!params[tr.param]) params[tr.param] = { label: tr.label || 'time $t$', min: tr.min || 0, max: tr.max, step: 'any', value: tr.value === undefined ? (tr.min || 0) : tr.value, unit: tr.unit || 's', optional: true };
      if (!params.speed) params.speed = { label: 'speed', min: 0.1, max: 3, step: 0.05, value: tr.speed || 1, optional: true, fmt: function (v) { return fmt(v, { sig: 2 }) + '×'; } };
    }
    if (cfg.modes && !params[cfg.modes.param]) params[cfg.modes.param] = { type: 'select', label: cfg.modes.label || 'System', options: cfg.modes.options, value: cfg.modes.value === undefined ? cfg.modes.options[0][0] : cfg.modes.value, optional: true };
    app.params = params;
    Object.keys(params).forEach(function (k) { var p = params[k]; app.defaults[k] = p.value; app.state[k] = p.value; });
    var animated = !!(cfg.stage && (cfg.stage.animate || tr));

    var theme = (STORE && STORE.getItem('viz-theme')) || cfg.theme || 'light';
    doc.documentElement.setAttribute('data-theme', theme);

    /* ---------- header ---------- */
    var tabsDef = DEFAULT_TABS.filter(function (t) {
      if (t.id === 'tour') return cfg.tour && cfg.tour.length;
      if (t.id === 'explore') return true;
      if (t.id === 'explain') return !!cfg.explain;
      if (t.id === 'derive') return !!(cfg.derivations && cfg.derivations.length);
      if (t.id === 'equations') return cfg.equations && cfg.equations.length;
      if (t.id === 'code') return cfg.code && cfg.code.length;
      if (t.id === 'check') return cfg.check && cfg.check.length;
      return false;
    });
    (cfg.panels || []).forEach(function (p) { tabsDef.splice(p.index === undefined ? tabsDef.length : p.index, 0, { id: p.id, label: p.label, short: p.short || p.label, custom: p }); });
    app.tabs = tabsDef.map(function (t) { return t.id; });

    var tabBtns = {};
    var nav = h('nav', { class: 'viz-tabs', role: 'tablist', 'aria-label': 'Sections' },
      tabsDef.map(function (t) {
        var b = h('button', { class: 'viz-tab', role: 'tab', type: 'button', 'aria-selected': 'false', 'data-tab': t.id, onclick: function () { app.setTab(t.id); } },
          h('span', { class: 'viz-tab-long', text: t.label }), h('span', { class: 'viz-tab-short', text: t.short || t.label }));
        tabBtns[t.id] = b; return b;
      }));
    var helpBtn = h('button', { class: 'viz-icon-btn', type: 'button', title: 'How to use this explainer', 'aria-label': 'Help', html: ICON.help });
    var themeBtn = h('button', { class: 'viz-icon-btn viz-optional', type: 'button', title: 'Light / dark', 'aria-label': 'Toggle theme', html: ICON.theme });
    var fullBtn = h('button', { class: 'viz-icon-btn', type: 'button', title: 'Full screen', 'aria-label': 'Full screen', html: ICON.full });
    var resetBtn = h('button', { class: 'viz-icon-btn', type: 'button', title: 'Reset all controls', 'aria-label': 'Reset', html: ICON.reset });
    var canFull = !!(doc.fullscreenEnabled || doc.webkitFullscreenEnabled);
    if (!canFull) fullBtn.hidden = true;
    var header = h('header', { class: 'viz-header' },
      h('div', { class: 'viz-titlebox' }, h('h1', { class: 'viz-title', html: mathify(cfg.title || doc.title) }), cfg.subtitle ? h('p', { class: 'viz-subtitle', html: mathify(cfg.subtitle) }) : null),
      nav, h('div', { class: 'viz-actions' }, resetBtn, themeBtn, fullBtn, helpBtn));

    /* ---------- controls (a param can be rendered in several places; all copies stay in sync) ---------- */
    function sliderPos(p, v) { return p.log ? Math.log10(v) : v; }
    function sliderVal(p, x) { var v = p.log ? Math.pow(10, x) : x; if (!p.log && p.step && p.step !== 'any') v = Math.round(v / p.step) * p.step; return Number(v.toPrecision(12)); }
    function register(key, el, sync) { (app.controls[key] = app.controls[key] || []).push({ el: el, sync: sync }); }
    function control(key, opt) {
      opt = opt || {};
      var p = params[key]; if (!p) throw new Error('Viz: unknown param "' + key + '"');
      var type = p.type || 'range', wrap, input, valEl;
      if (type === 'toggle') {
        input = h('input', { type: 'checkbox' });
        input.checked = !!app.state[key];
        input.addEventListener('change', function () { app.set(key, input.checked); });
        wrap = h('label', { class: 'viz-chip-toggle' + (p.optional ? ' viz-optional' : ''), 'data-viz-key': 'param:' + key, title: p.help || '' }, input, h('span', { html: mathify(p.label) }));
        register(key, wrap, function (v) { input.checked = !!v; });
        return wrap;
      }
      if (type === 'select' && (p.display === 'chips' || opt.chips)) {
        wrap = h('div', { class: 'viz-seg', role: 'radiogroup', 'aria-label': p.label.replace(/\$/g, ''), 'data-viz-key': 'param:' + key },
          (p.options || []).map(function (o) { return h('button', { class: 'viz-seg-btn', type: 'button', 'data-value': String(o[0]), html: mathify(o[1]), onclick: function () { app.set(key, o[0]); } }); }));
        var syncSeg = function (v) { Array.prototype.forEach.call(wrap.children, function (b) { b.setAttribute('aria-pressed', b.getAttribute('data-value') === String(v) ? 'true' : 'false'); }); };
        syncSeg(app.state[key]); register(key, wrap, syncSeg);
        return wrap;
      }
      if (type === 'select') {
        input = h('select', {}, (p.options || []).map(function (o) { return h('option', { value: String(o[0]), text: o[1].replace(/\$/g, '') }); }));
        input.value = String(app.state[key]);
        input.addEventListener('change', function () { var o = (p.options || []).filter(function (x) { return String(x[0]) === input.value; })[0]; app.set(key, o ? o[0] : input.value); });
        wrap = h('label', { class: 'viz-select' + (p.optional ? ' viz-optional' : ''), 'data-viz-key': 'param:' + key }, h('span', { html: mathify(p.label) }), input);
        register(key, wrap, function (v) { input.value = String(v); });
        return wrap;
      }
      if (type === 'button') {
        return h('button', { class: 'viz-btn' + (p.primary ? ' primary' : '') + (p.optional ? ' viz-optional' : ''), type: 'button', 'data-viz-key': 'param:' + key, html: mathify(p.label), onclick: function () { if (isFn(p.action)) p.action(app); } });
      }
      input = h('input', { type: 'range', min: sliderPos(p, p.min), max: sliderPos(p, p.max), step: p.log ? 'any' : (p.step || 'any'), 'aria-label': p.label.replace(/\$/g, '') });
      input.value = sliderPos(p, app.state[key]);
      valEl = h('output', { class: 'viz-control-value' });
      var show = function (v) { valEl.textContent = p.fmt ? p.fmt(v) : fmt(v, { sig: p.sig || 3, unit: p.unit }); };
      show(app.state[key]);
      input.addEventListener('input', function () { if (tr && key === tr.param) app.play(false); app.set(key, sliderVal(p, Number(input.value))); });
      wrap = h('div', { class: 'viz-control' + (p.optional && !opt.keep ? ' viz-optional' : ''), 'data-viz-key': 'param:' + key },
        h('span', { class: 'viz-control-label', html: mathify(p.label) }), valEl, input,
        p.help && !opt.compact ? h('span', { class: 'viz-control-help', html: mathify(p.help) }) : null);
      register(key, wrap, function (v) { input.value = sliderPos(p, v); show(v); });
      return wrap;
    }
    app.control = control;

    /* ---------- stage: top strip (modes · presets · status) + linked views + transport ---------- */
    var viewsCfg = (cfg.stage && cfg.stage.views && cfg.stage.views.length) ? cfg.stage.views : [{ id: 'main', bare: true }];
    var multi = viewsCfg.length > 1 || !viewsCfg[0].bare;
    var views = [], viewById = {};
    var rowsEl = h('div', { class: 'viz-views' });
    var rowIds = [];
    viewsCfg.forEach(function (vc) { var r = vc.row || 0; if (rowIds.indexOf(r) < 0) rowIds.push(r); });
    rowIds.sort(function (a, b) { return a - b; });
    rowIds.forEach(function (r, ri) {
      var weight = (cfg.stage && cfg.stage.rows && cfg.stage.rows[ri]) || 1;
      var rowEl = h('div', { class: 'viz-row', style: { flex: weight + ' 1 0' } });
      viewsCfg.filter(function (vc) { return (vc.row || 0) === r; }).forEach(function (vc) {
        var canvas = h('canvas', { class: 'viz-canvas', role: 'img', 'aria-label': (vc.label || vc.id || cfg.title || 'visualisation').replace(/<[^>]+>/g, '') });
        var holder = h('div', { class: 'viz-view-holder' }, canvas);
        var titleEl = vc.title ? h('div', { class: 'viz-view-title' }) : null;
        var el = h('div', { class: 'viz-view' + (vc.bare ? ' bare' : '') + (vc.hidePortrait ? ' viz-hide-portrait' : ''), 'data-view': vc.id, style: { flex: (vc.flex || 1) + ' 1 0' } }, titleEl, holder);
        rowEl.appendChild(el);
        var v = { id: vc.id, cfg: vc, el: el, holder: holder, canvas: canvas, ctx: canvas.getContext('2d'), titleEl: titleEl, w: 0, h: 0, dpr: 1, pointer: { x: 0, y: 0, down: false, inside: false } };
        v.plot = function (opt) { return Plot(v, opt); };
        v.clear = function (color) { v.ctx.save(); v.ctx.setTransform(v.dpr, 0, 0, v.dpr, 0, 0); if (color) { v.ctx.fillStyle = color; v.ctx.fillRect(0, 0, v.w, v.h); } else v.ctx.clearRect(0, 0, v.w, v.h); v.ctx.restore(); };
        views.push(v); viewById[vc.id] = v;
      });
      rowsEl.appendChild(rowEl);
    });
    var badge = h('div', { class: 'viz-stage-badge' });
    views[0].holder.appendChild(badge);
    var playBtn = h('button', { class: 'viz-icon-btn viz-play', type: 'button', title: 'Play / pause (space)', 'aria-label': 'Play or pause', html: ICON.play });

    var topEl = null, statusEl = null, presetsEl = null;
    if (cfg.modes || (cfg.presets && cfg.presets.length) || cfg.status) {
      topEl = h('div', { class: 'viz-stage-top' });
      if (cfg.modes) topEl.appendChild(control(cfg.modes.param, { chips: true }));
      if (cfg.presets && cfg.presets.length) {
        presetsEl = h('div', { class: 'viz-presets', 'aria-label': 'Presets' }, h('span', { class: 'viz-presets-label', text: cfg.presetsLabel || 'Try' }),
          cfg.presets.map(function (pr) { return h('button', { class: 'viz-preset', type: 'button', title: (pr.title || '').replace(/\$/g, ''), html: mathify(pr.label), onclick: function () { app.set(pr.set || {}); if (isFn(pr.action)) pr.action(app); } }); }));
        topEl.appendChild(presetsEl);
      }
      if (cfg.status) { statusEl = h('div', { class: 'viz-status', role: 'status' }); topEl.appendChild(statusEl); }
    }
    var transportEl = null;
    if (tr) {
      var restartBtn = h('button', { class: 'viz-icon-btn', type: 'button', title: 'Restart', 'aria-label': 'Restart', html: ICON.reset, onclick: function () { app.holding = false; holdLeft = 0; app.set(tr.param, params[tr.param].min); app.play(true); } });
      var stepBack = h('button', { class: 'viz-icon-btn viz-optional', type: 'button', title: 'Step back', 'aria-label': 'Step back', text: '⏮', onclick: function () { app.play(false); app.set(tr.param, app.state[tr.param] - (tr.step || (params[tr.param].max - params[tr.param].min) / 50)); } });
      var stepFwd = h('button', { class: 'viz-icon-btn viz-optional', type: 'button', title: 'Step forward', 'aria-label': 'Step forward', text: '⏭', onclick: function () { app.play(false); app.set(tr.param, app.state[tr.param] + (tr.step || (params[tr.param].max - params[tr.param].min) / 50)); } });
      var scrub = control(tr.param, { keep: true, compact: true }); scrub.classList.add('viz-scrub');
      var speedCtl = control('speed', { compact: true }); speedCtl.classList.add('viz-speed', 'viz-hide-portrait');
      transportEl = h('div', { class: 'viz-transport' }, playBtn, stepBack, stepFwd, restartBtn, scrub, speedCtl);
    }
    var stageEl = h('section', { class: 'viz-stage' + (multi ? ' multi' : ''), 'data-views': String(views.length) }, topEl, rowsEl, transportEl);
    if (animated && !tr) views[0].holder.appendChild(h('div', { class: 'viz-stage-bar' }, playBtn));
    var side = h('aside', { class: 'viz-side' });
    var main = h('main', { class: 'viz-main' }, stageEl, side);
    var help = h('div', { class: 'viz-help', role: 'dialog', 'aria-label': 'Help' });
    var root = h('div', { class: 'viz-app', 'data-layout': 'wide', 'data-dense': '0' }, header, main, help);
    mount.appendChild(root);
    app.root = root; app.stageEl = stageEl; app.canvas = views[0].canvas; app.badge = badge; app.views = views;

    /* ---------- readouts ---------- */
    var readoutEls = [];
    function readoutsGrid(ids) {
      var list = (cfg.readouts || []).filter(function (r) { return !ids || ids.indexOf(r.id) >= 0; });
      if (!list.length) return null;
      var grid = h('div', { class: 'viz-readouts' });
      list.forEach(function (r) {
        var v = h('div', { class: 'viz-readout-value' });
        var el = h('div', { class: 'viz-readout' + (r.optional ? ' viz-optional' : ''), 'data-viz-key': 'readout:' + r.id, title: r.help || '' }, h('div', { class: 'viz-readout-label', html: mathify(r.label) }), v);
        readoutEls.push({ r: r, el: el, v: v }); grid.appendChild(el);
      });
      return grid;
    }

    /* ---------- dynamic cards: terms · inspector · notes ---------- */
    var termEls = [];
    function termsBlock() {
      var tc = cfg.terms; if (!tc) return null;
      var box = h('div', { class: 'viz-terms' });
      tc.items.forEach(function (it, idx) {
        var val = h('span', { class: 'viz-term-val' }), fill = h('span', { class: 'viz-term-fill', style: { background: it.color || Viz.color('accent') } });
        var row = h('button', { class: 'viz-term', type: 'button', 'data-viz-key': 'term:' + it.id, title: (it.help || '').replace(/\$/g, ''), onclick: function () { app.set('focus', app.state.focus === it.id ? null : it.id); } },
          h('span', { class: 'viz-term-head' }, h('span', { class: 'viz-term-label', html: mathify(it.label) }), val), h('span', { class: 'viz-term-track' }, fill));
        row.style.setProperty('--term-color', it.color || Viz.color('accent'));
        termEls.push({ it: it, idx: idx, row: row, val: val, fill: fill }); box.appendChild(row);
      });
      if (tc.total) { var tv = h('span', { class: 'viz-term-val' }); box.appendChild(h('div', { class: 'viz-term-total' }, h('span', { html: mathify(tc.total.label) }), tv)); termEls.push({ total: true, val: tv }); }
      return box;
    }
    function updateTerms() {
      termEls = termEls.filter(function (o) { return o.val.isConnected; });   // drop bars of earlier walkthrough steps
      var tc = cfg.terms; if (!tc || !termEls.length) return;
      var vals = tc.items.map(function (it) { try { return it.value(app.state, app); } catch (e) { return NaN; } });
      var total = tc.total ? (isFn(tc.total.value) ? tc.total.value(app.state, app) : vals.reduce(function (a, b) { return a + (isFinite(b) ? b : 0); }, 0)) : null;
      var mx = Math.max.apply(null, vals.map(Math.abs).concat(total === null ? [] : [Math.abs(total)]).filter(isFinite).concat([1e-300]));
      termEls.forEach(function (o, i) {
        if (o.total) { o.val.textContent = fmt(total, { sig: tc.sig || 3, unit: tc.unit }); return; }
        var v = vals[o.idx];
        o.val.textContent = fmt(v, { sig: tc.sig || 3, unit: tc.unit });
        o.fill.style.width = (isFinite(v) ? Math.abs(v) / mx * 100 : 0) + '%';
        o.row.classList.toggle('neg', v < 0);
        o.row.setAttribute('aria-pressed', app.state.focus === o.it.id ? 'true' : 'false');
      });
    }
    var inspectBody = null, notesBody = null, lastNotes = null, lastInspect = null;
    function updateInspectNotes() {
      if (inspectBody && isFn(cfg.inspect)) {
        var html; try { html = cfg.inspect(app.state, app); } catch (e) { html = null; }
        html = html || '<span class="viz-muted">' + mathify(cfg.inspectHint || 'Click the picture to inspect a point: its exact arithmetic appears here.') + '</span>';
        if (html !== lastInspect) { lastInspect = html; inspectBody.innerHTML = mathify(html); typeset(inspectBody); }
      }
      if (notesBody && isFn(cfg.notes)) {
        var nh; try { nh = cfg.notes(app.state, app); } catch (e) { nh = ''; }
        if (nh !== lastNotes) { lastNotes = nh; notesBody.innerHTML = mathify(nh); typeset(notesBody); if (explorePager && app.tab === 'explore') explorePager.layout(); }
      }
    }

    /* ---------- pager: packs children into pages that fit the box (no scrollbars, ever) ---------- */
    function Pager(title, items, headEl, extraClass) {
      var body = h('div', { class: 'viz-pages' }), label = h('span', { class: 'viz-pager-label' });
      var prev = h('button', { class: 'viz-pager-btn', type: 'button', 'aria-label': 'Previous page', text: '‹' });
      var next = h('button', { class: 'viz-pager-btn', type: 'button', 'aria-label': 'Next page', text: '›' });
      var pager = h('span', { class: 'viz-pager', hidden: true }, prev, label, next);
      var card = h('div', { class: 'viz-card grow viz-paged' + (extraClass ? ' ' + extraClass : '') }, h('div', { class: 'viz-card-title' }, headEl || h('span', { html: mathify(title) }), pager), body);
      var pages = [0], page = 0;
      function show(p) {
        page = clamp(p, 0, pages.length - 1);
        items.forEach(function (it) { it.style.display = (it._vizPage === page) ? '' : 'none'; });
        label.textContent = (page + 1) + ' / ' + pages.length; prev.disabled = page === 0; next.disabled = page === pages.length - 1;
      }
      prev.onclick = function () { show(page - 1); }; next.onclick = function () { show(page + 1); };
      var P = {
        card: card, body: body,
        setItems: function (list) { items = list; body.innerHTML = ''; list.forEach(function (it) { body.appendChild(it); }); page = 0; },
        layout: function () {
          items.forEach(function (it) { it.style.display = ''; it._vizPage = 0; });
          var avail = body.clientHeight; if (avail <= 0) return;
          var top0 = body.getBoundingClientRect().top, pageStart = 0, pg = 0; pages = [0];
          items.forEach(function (it, i) {
            var r = it.getBoundingClientRect(), top = r.top - top0, bottom = r.bottom - top0;
            if (bottom - pageStart > avail + 0.5 && top > pageStart + 0.5) {
              // keep a sub-heading together with the item that follows it
              var prevIt = items[i - 1];
              if (prevIt && (prevIt.classList.contains('viz-subhead') || prevIt.classList.contains('viz-work-h')) && prevIt._vizPage === pg && prevIt.getBoundingClientRect().top - top0 > pageStart + 0.5) { pg += 1; pageStart = prevIt.getBoundingClientRect().top - top0; prevIt._vizPage = pg; }
              else { pg += 1; pageStart = top; }
              pages.push(pg);
            }
            it._vizPage = pg;
          });
          pages = []; for (var q = 0; q <= pg; q++) pages.push(q);
          pager.hidden = pages.length < 2; show(Math.min(page, pages.length - 1));
        }
      };
      items.forEach(function (it) { body.appendChild(it); });
      return P;
    }
    var pagers = [];
    function subhead(text) { return h('div', { class: 'viz-subhead', html: mathify(text) }); }

    /* ---------- panels ---------- */
    var panels = {};
    function panel(id) { var p = h('div', { class: 'viz-panel', role: 'tabpanel', 'data-panel': id }); side.appendChild(p); panels[id] = p; return p; }

    // Explore: one pager so it always fits — intro · controls · presets · live numbers · terms · inspector · right now · callouts
    var explorePager = null;
    (function () {
      var p = panel('explore'), ex = cfg.explore || {}, items = [];
      if (ex.intro) items.push(h('div', { class: 'viz-intro', html: mathify(ex.intro) }));
      var keys = (ex.controls || Object.keys(params)).filter(function (k) { return !(tr && (k === tr.param || k === 'speed')) && !(cfg.modes && k === cfg.modes.param); });
      if (keys.length) { items.push(subhead('Controls')); keys.forEach(function (k) { items.push(control(k)); }); }
      var ro = readoutsGrid(ex.readouts);
      if (ro) { items.push(subhead('Live numbers')); items.push(ro); }
      var tb = termsBlock();
      if (tb) { items.push(subhead(cfg.terms.title || 'Term by term')); items.push(tb); }
      if (isFn(cfg.inspect)) { inspectBody = h('div', { class: 'viz-inspect-body' }); items.push(subhead(cfg.inspectTitle || 'Inspector')); items.push(h('div', { class: 'viz-inspect' }, inspectBody)); }
      if (isFn(cfg.notes)) { notesBody = h('div', { class: 'viz-notes-body' }); items.push(subhead(cfg.notesTitle || 'Right now')); items.push(h('div', { class: 'viz-notes' }, notesBody)); }
      (ex.callouts || []).forEach(function (c) { items.push(h('div', { class: 'viz-callout ' + (c.kind || 'key') + (c.optional === false ? '' : ' viz-optional') }, h('b', { class: 'k', text: CALLOUT_LABEL[c.kind || 'key'] }), h('span', { html: mathify(c.html) }))); });
      explorePager = Pager(ex.title || 'Explore', items); pagers.push(explorePager); p.appendChild(explorePager.card);
    })();

    // Walkthrough
    var tourEls = {};
    if (cfg.tour && cfg.tour.length) (function () {
      var p = panel('tour');
      var chips = h('div', { class: 'viz-steps-chips', role: 'tablist', 'aria-label': 'Steps' }, cfg.tour.map(function (s, i) { return h('button', { class: 'viz-step-chip', type: 'button', text: String(i + 1), title: s.title.replace(/\$/g, ''), onclick: function () { app.goStep(i); } }); }));
      var count = h('div', { class: 'viz-step-count' }), title = h('h2', { class: 'viz-step-title' }), text = h('div', { class: 'viz-step-text' });
      var prev = h('button', { class: 'viz-btn', type: 'button', text: '← Back', onclick: function () { app.goStep(app.step - 1); } });
      var next = h('button', { class: 'viz-btn primary', type: 'button', text: 'Next →', onclick: function () { app.goStep(app.step + 1); } });
      // the step card is a pager: on small screens a step's extras (controls, equation, code…) continue on "›"
      var tourPager = Pager('', [title, text], h('div', { class: 'viz-step-head' }, count, chips), 'viz-step-card');
      pagers.push(tourPager);
      p.appendChild(tourPager.card);
      p.appendChild(h('div', { class: 'viz-step-nav' }, prev, next));
      tourEls = { chips: chips, count: count, title: title, text: text, extras: h('div'), prev: prev, next: next, card: tourPager.card, pager: tourPager };
    })();

    // Explain — "Explanation & interpretation" (the forced_damped_vibrations.html pattern): numbered sections that
    // compute every displayed quantity with the reader's own numbers, then interpret the current setting. Live.
    var calcPager = null, calcHtml = null, calcLive = [], popWin = null, popBtn = null;
    if (cfg.explain) (function () {
      var p = panel('explain');
      popBtn = h('button', { class: 'viz-pager-btn viz-pop-btn', type: 'button', title: 'Open this explanation in a new tab (it stays live)', 'aria-label': 'Open the explanation in a new tab', text: '↗', onclick: function () { openPop(); } });
      var head = h('span', { class: 'viz-card-head' }, h('span', { html: mathify(cfg.explain.title || 'Explanation & interpretation') }), popBtn);
      calcPager = Pager('', [], head); pagers.push(calcPager); p.appendChild(calcPager.card);
    })();
    function updateCalc(structural) {
      if (!calcPager) return;
      var E = cfg.explain;
      if (structural) {
        var html; try { html = (isFn(E) ? E : E.html)(app.state, app); } catch (e) { html = '<p>—</p>'; reportError(e); }
        if (html !== calcHtml) {
          calcHtml = html;
          var tmp = h('div', { html: mathify(html) });
          var items = Array.prototype.slice.call(tmp.children);
          calcPager.setItems(items); typeset(calcPager.body);
          calcLive = Array.prototype.slice.call(calcPager.body.querySelectorAll('[data-live]'));
          if (app.tab === 'explain') calcPager.layout();
        }
      }
      if (E.live && calcLive.length) {
        var vals; try { vals = E.live(app.state, app) || {}; } catch (e) { vals = {}; }
        calcLive.forEach(function (el) { var v = vals[el.getAttribute('data-live')]; if (v !== undefined) { v = typeof v === 'number' ? fmt(v, { sig: 4 }) : String(v); if (el.textContent !== v) el.textContent = v; } });
      }
      pushPop();
    }
    // the same explanation in a separate, live browser tab (for reading next to the picture on a big screen)
    function openPop() {
      if (popWin && !popWin.closed) { popWin.focus(); pushPop(); return; }
      var w = null; try { w = global.open('', '_blank'); } catch (e) { w = null; }
      if (!w) { popBtn.title = 'The browser blocked the new tab: allow pop-ups for this page'; popBtn.classList.add('blocked'); return; }
      popWin = w;
      var styles = Array.prototype.map.call(doc.querySelectorAll('style, link[rel="stylesheet"]'), function (el) { return el.outerHTML; }).join('\n');
      var name = esc(String(cfg.title || doc.title).replace(/\$/g, ''));
      w.document.write('<!DOCTYPE html><html lang="en" data-theme="' + (doc.documentElement.getAttribute('data-theme') || 'light') + '"><head><meta charset="utf-8">' +
        '<meta name="viewport" content="width=device-width, initial-scale=1"><title>Explanation · ' + name + '</title>' + styles +
        '<style>html,body{height:auto!important;overflow:auto!important}body{margin:0;background:var(--viz-bg)}' +
        '.viz-pop-hd{position:sticky;top:0;padding:8px 18px;font-weight:650;background:var(--viz-card);border-bottom:1px solid var(--viz-border)}' +
        '.viz-pop-hd span{font-weight:400;color:var(--viz-muted);font-size:12px;margin-left:8px}' +
        '.viz-pop{padding:12px 18px;column-width:340px;column-gap:28px;font-size:14px}.viz-pop>*{break-inside:avoid;margin-bottom:6px}</style></head>' +
        '<body><div class="viz-pop-hd">' + name + ' · explanation &amp; interpretation<span>live: follows the controls and the animation in the explainer</span></div>' +
        '<div class="viz-pop" id="viz-pop-body"></div></body></html>');
      w.document.close(); pushPop();
    }
    function pushPop() {
      if (!popWin) return;
      if (popWin.closed) { popWin = null; return; }
      try {
        var b = popWin.document.getElementById('viz-pop-body'); if (!b) return;
        b.innerHTML = calcPager.body.innerHTML;
        Array.prototype.forEach.call(b.children, function (el) { el.style.display = ''; });
      } catch (e) { popWin = null; }
    }

    // Derivation — one move per step: the line we had → what we did → the new line → why that is allowed → what it
    // says in words (+ the line with the reader's numbers). Pages: the goal · step 1…n · the result.
    var der = { i: 0, s: 0, els: null };
    function derBlock(cls, label, html) {
      return h('div', { class: 'viz-der-block ' + cls }, label ? h('b', { class: 'k', html: mathify(label) }) : null, h('span', { class: 'viz-der-txt', html: mathify(html) }));
    }
    function derLine(tex, cls, label) {
      var inner = h('div', { class: 'viz-eq-inner' }); renderTex(inner, tex, true);
      return h('div', { class: 'viz-der-line ' + (cls || '') }, label ? h('span', { class: 'viz-der-tag', text: label }) : null, h('div', { class: 'viz-eq-body' }, inner));
    }
    if (cfg.derivations && cfg.derivations.length) (function () {
      var p = panel('derive');
      var pick = cfg.derivations.length > 1 ? h('div', { class: 'viz-seg viz-der-pick', role: 'radiogroup', 'aria-label': 'Derivation' },
        cfg.derivations.map(function (d, i) { return h('button', { class: 'viz-seg-btn', type: 'button', html: mathify(d.short || d.title), onclick: function () { app.goDerive(i, 0); } }); })) : null;
      var count = h('div', { class: 'viz-step-count' }), chips = h('div', { class: 'viz-steps-chips', role: 'tablist', 'aria-label': 'Derivation steps' });
      var pg = Pager('', [], h('div', { class: 'viz-step-head' }, pick, count, chips), 'viz-step-card viz-der-card'); pagers.push(pg);
      var prev = h('button', { class: 'viz-btn', type: 'button', text: '← Back', onclick: function () { app.goDerive(der.i, der.s - 1); } });
      var next = h('button', { class: 'viz-btn primary', type: 'button', text: 'Next →' });
      p.appendChild(pg.card); p.appendChild(h('div', { class: 'viz-step-nav' }, prev, next));
      der.els = { pick: pick, count: count, chips: chips, pager: pg, prev: prev, next: next };
    })();
    /** Show derivation i at page s (0 = the goal, 1…n = steps, n+1 = the result). opts.quiet: render only (no state change). */
    app.goDerive = function (i, s, opts) {
      if (!der.els) return app;
      opts = opts || {};
      var D = cfg.derivations; i = clamp(i || 0, 0, D.length - 1);
      var d = D[i], n = d.steps.length, last = n + 1;
      s = clamp(s || 0, 0, last); der.i = i; der.s = s;
      var st = s >= 1 && s <= n ? d.steps[s - 1] : null;
      if (!opts.quiet) {
        if (app.tab !== 'derive') app.setTab('derive');
        if (s === 0 && d.set) app.set(d.set);
        if (st && st.set) app.set(st.set);
        if (s === last && d.result && d.result.set) app.set(d.result.set);
        if (st && st.play !== undefined) app.play(st.play);
      }
      root.setAttribute('data-der-picture', d.picture === false ? '0' : '1');
      // on phones only one view stays next to the derivation: d.view, else the first view
      var keep = viewById[d.view] ? d.view : views[0].id;
      views.forEach(function (v) { v.el.classList.toggle('viz-der-off', v.id !== keep); });
      Array.prototype.forEach.call(rowsEl.children, function (row) { row.classList.toggle('viz-der-off', !row.querySelector('.viz-view:not(.viz-der-off)')); });
      var E = der.els, items = [];
      if (E.pick) Array.prototype.forEach.call(E.pick.children, function (b, k) { b.setAttribute('aria-pressed', k === i ? 'true' : 'false'); });
      E.count.textContent = (s === 0 ? 'The goal' : s === last ? 'The result' : 'Step ' + s + ' of ' + n) + (d.ref ? ' · ' + d.ref : '');
      E.chips.innerHTML = '';
      for (var k = 0; k <= last; k++) (function (k) {
        var c = h('button', { class: 'viz-step-chip' + (k < s ? ' done' : ''), type: 'button', text: k === 0 ? 'G' : (k === last ? '✓' : String(k)),
          title: k === 0 ? 'The goal' : (k === last ? 'The result' : String(d.steps[k - 1].did || 'Step ' + k).replace(/<[^>]+>|\$/g, '')), onclick: function () { app.goDerive(i, k); } });
        if (k === s) c.setAttribute('aria-current', 'step');
        E.chips.appendChild(c);
      })(k);
      items.push(h('h2', { class: 'viz-step-title', html: mathify(d.title) }));
      var startTex = d.start ? (typeof d.start === 'string' ? d.start : d.start.tex) : null;
      if (s === 0) {
        if (d.goal) items.push(derBlock('viz-der-goal', 'What we want to show', d.goal));
        if (startTex) items.push(derLine(startTex, 'viz-der-cur', 'we start from'));
        if (d.start && d.start.plain) items.push(derBlock('viz-der-plain', 'In words', d.start.plain));
        if (d.plan && d.plan.length) items.push(derBlock('viz-der-plan', 'The plan', '<ol>' + d.plan.map(function (x) { return '<li>' + x + '</li>'; }).join('') + '</ol>'));
        if (d.uses && d.uses.length) items.push(derBlock('viz-der-uses', 'Tools we use', d.uses.join(' · ')));
      } else if (st) {
        var prevTex = s === 1 ? startTex : d.steps[s - 2].tex;
        if (prevTex) items.push(derLine(prevTex, 'viz-der-prev', 'we had'));
        if (st.did) items.push(derBlock('viz-der-did', '', '↓ ' + st.did));
        items.push(derLine(st.tex, 'viz-der-cur', 'now'));
        if (st.why) items.push(derBlock('viz-der-why', 'Why', st.why));
        if (st.plain) items.push(derBlock('viz-der-plain', 'In words', st.plain));
        if (isFn(st.live)) { var lv = h('div', { class: 'viz-eq-inner', 'data-der-live': '1' }); items.push(h('div', { class: 'viz-der-line viz-der-live' }, h('span', { class: 'viz-der-tag', text: 'your numbers' }), h('div', { class: 'viz-eq-body' }, lv))); }
        if (st.watch) items.push(h('div', { class: 'viz-callout watch' }, h('b', { class: 'k', text: 'Watch' }), h('span', { html: mathify(st.watch) })));
      } else {
        items.push(derBlock('viz-der-goal', 'The whole chain', 'Each line follows from the one above by the move written next to it.'));
        if (startTex) items.push(derLine(startTex, 'viz-der-chain'));
        d.steps.forEach(function (x, k) { items.push(derLine(x.tex, 'viz-der-chain', String(k + 1))); });
        if (d.result && d.result.tex) items.push(derLine(d.result.tex, 'viz-der-result', 'result'));
        if (d.result && d.result.plain) items.push(derBlock('viz-der-plain', 'In words', d.result.plain));
        if (isFn(d.interpret)) items.push(h('div', { class: 'viz-der-block viz-der-interp' }, h('b', { class: 'k', text: 'What it means right now' }), h('span', { class: 'viz-der-txt viz-der-interp-body' })));
        if (d.check) items.push(derBlock('viz-der-check', 'Check it', d.check));
      }
      E.pager.setItems(items); typeset(E.pager.card);
      E.prev.disabled = s === 0;
      E.next.textContent = s === 0 ? 'Start →' : (s === n ? 'The result →' : (s === last ? (i < D.length - 1 ? 'Next derivation →' : 'Explore →') : 'Next step →'));
      E.next.onclick = s === last ? (i < D.length - 1 ? function () { app.goDerive(i + 1, 0); } : function () { app.setTab('explore'); }) : function () { app.goDerive(i, s + 1); };
      if (!opts.quiet) {
        clearHighlights();
        var hl = (st && st.highlight) || (s === last && d.result && d.result.highlight);
        if (hl) [].concat(hl).forEach(function (key) { Array.prototype.forEach.call(root.querySelectorAll('[data-viz-key="' + key + '"]'), function (el) { el.classList.add('viz-hl'); app._hl.push(el); }); });
        updateDerDynamic(); app.draw(); app.fit(); saveHash();
      }
      return app;
    };
    function updateDerDynamic() {
      if (!der.els) return;
      var d = cfg.derivations[der.i], st = der.s >= 1 && der.s <= d.steps.length ? d.steps[der.s - 1] : null;
      var lv = der.els.pager.body.querySelector('[data-der-live]');
      if (lv && st && isFn(st.live)) { var t; try { t = st.live(app.state, app); } catch (e) { t = '\\text{—}'; } if (lv._t !== t) { lv._t = t; renderTex(lv, t, true); scaleEq(lv); } }
      var ib = der.els.pager.body.querySelector('.viz-der-interp-body');
      if (ib && isFn(d.interpret)) { var hh; try { hh = d.interpret(app.state, app); } catch (e) { hh = ''; } if (ib._h !== hh) { ib._h = hh; ib.innerHTML = mathify(hh); typeset(ib); } }
    }

    // Equations
    var eqEls = {};
    if (cfg.equations && cfg.equations.length) (function () {
      var p = panel('equations');
      var items = cfg.equations.map(function (e) {
        var body = h('div', { class: 'viz-eq-body' }), inner = h('div', { class: 'viz-eq-inner' }); body.appendChild(inner);
        var live = e.live ? h('div', { class: 'viz-eq-live' }, h('div', { class: 'viz-eq-inner' })) : null;
        var sym = e.symbols && e.symbols.length ? h('div', { class: 'viz-symbols viz-optional' }, [].concat.apply([], e.symbols.map(function (r) {
          return [h('span', { class: 'sym viz-tex', 'data-tex': r[0] }), h('span', { html: mathify(r[1]) }), h('span', { class: 'unit', text: r[2] || '' })];
        }))) : null;
        var el = h('div', { class: 'viz-eq', 'data-viz-key': 'eq:' + e.id },
          h('div', { class: 'viz-eq-head' }, h('span', { html: mathify(e.title) }), e.ref ? h('span', { class: 'viz-eq-ref', text: e.ref }) : null),
          body, live, e.note ? h('div', { class: 'viz-eq-note', html: mathify(e.note) }) : null, sym);
        eqEls[e.id] = { cfg: e, el: el, inner: inner, live: live && live.firstChild };
        renderTex(inner, e.tex, true);
        return el;
      });
      var pg = Pager('Equations · live with your numbers', items); pagers.push(pg); p.appendChild(pg.card); eqEls._pager = pg;
    })();

    // Code
    var codeCfg = {};
    function codeLines(c) { return String(c.src).replace(/^\n+|\s+$/g, '').split('\n'); }
    function codeBlock(c, from, to, hl) {
      var lines = codeLines(c), pre = h('div', { class: 'viz-code', 'data-code': c.id });
      from = from === undefined ? 0 : Math.max(0, from); to = to === undefined ? lines.length - 1 : Math.min(lines.length - 1, to);
      for (var i = from; i <= to; i++) {
        var on = hl && hl[0] <= i + 1 && i + 1 <= hl[1];
        pre.appendChild(h('div', { class: 'viz-code-ln' + (on ? ' cur' : ''), html: '<span class="viz-code-no">' + (i + 1) + '</span>' + highlightPy(lines[i]) }));
      }
      return pre;
    }
    if (cfg.code && cfg.code.length) (function () {
      var p = panel('code'), items = [];
      cfg.code.forEach(function (c) {
        codeCfg[c.id] = c;
        // title, then the code in chunks of 5 lines, so a long listing pages instead of overflowing
        items.push(h('div', { class: 'viz-code-title viz-subhead-like', 'data-viz-key': 'code:' + c.id }, h('span', { html: mathify(c.title || c.id) }), c.ref ? h('span', { class: 'viz-eq-ref', text: c.ref }) : null));
        var n = codeLines(c).length;
        for (var k = 0; k < n; k += 5) items.push(codeBlock(c, k, Math.min(n - 1, k + 4)));
        if (c.note) items.push(h('div', { class: 'viz-eq-note', html: mathify(c.note) }));
      });
      var pg = Pager(cfg.codeTitle || 'The Python behind the picture', items); pagers.push(pg); p.appendChild(pg.card);
    })();
    function updateCodeLive() {
      var els = root.querySelectorAll('[data-code-live]'); if (!els.length) return;
      var vals = {};
      Object.keys(codeCfg).forEach(function (id) { var c = codeCfg[id]; if (isFn(c.live)) { try { Object.assign(vals, c.live(app.state, app) || {}); } catch (e) { /* ignore */ } } });
      Array.prototype.forEach.call(els, function (el) { var v = vals[el.getAttribute('data-code-live')]; if (v !== undefined) { v = typeof v === 'number' ? fmt(v, { sig: 4 }) : String(v); if (el.textContent !== v) el.textContent = v; } });
    }

    // Check yourself
    if (cfg.check && cfg.check.length) (function () {
      var p = panel('check');
      var items = cfg.check.map(function (c, i) {
        var q = h('div', { class: 'viz-q' });
        q.appendChild(h('p', { class: 'viz-q-text', html: mathify((i + 1) + '. ' + c.q) }));
        var row = h('div', { class: 'viz-btn-row' });
        row.appendChild(h('button', { class: 'viz-btn viz-q-reveal', type: 'button', text: 'Show answer', onclick: function () { q.classList.add('revealed'); typeset(q); app.fit(); } }));
        if (c.set) row.appendChild(h('button', { class: 'viz-btn', type: 'button', text: 'Try it on the picture', onclick: function () { app.set(c.set); } }));
        q.appendChild(row);
        q.appendChild(h('div', { class: 'viz-q-answer', html: mathify(c.a) }));
        return q;
      });
      var pg = Pager('Check yourself', items); pagers.push(pg); p.appendChild(pg.card);
    })();

    (cfg.panels || []).forEach(function (pc) { var p = panel(pc.id); if (isFn(pc.render)) pc.render(p, app); });

    help.innerHTML = '<h2>How to use this explainer</h2><ul>' +
      (cfg.tour && cfg.tour.length ? '<li><b>Walkthrough</b> tells the story step by step: press <kbd>Next →</kbd> or use <kbd>←</kbd>/<kbd>→</kbd>. Each step sets up the picture for you.</li>' : '') +
      '<li><b>Explore</b> hands you the controls: drag a slider and watch the picture and the live numbers respond.</li>' +
      (cfg.presets && cfg.presets.length ? '<li>The chips above the picture jump to special cases worth seeing.</li>' : '') +
      (cfg.explain ? '<li><b>Explain</b> works out every number on the picture with your current settings, step by step, and says what it means (↗ opens it in its own tab on a big screen).</li>' : '') +
      (cfg.derivations && cfg.derivations.length ? '<li><b>Derivation</b> builds the key formula one small move at a time: the line before, what we did, the new line, why that is allowed, and what it says in words.</li>' : '') +
      (cfg.equations && cfg.equations.length ? '<li><b>Equations</b> shows the formulas behind the picture, with your numbers substituted.</li>' : '') +
      (cfg.code && cfg.code.length ? '<li><b>Code</b> shows the Python that computes the picture; live values appear in the comments.</li>' : '') +
      (cfg.check && cfg.check.length ? '<li><b>Check yourself</b> asks questions you can answer by experimenting.</li>' : '') +
      (animated ? '<li><kbd>Space</kbd> plays or pauses' + (tr ? '; drag the time slider to scrub.' : '.') + '</li>' : '') +
      '<li>Reset (↺) restores the starting values. ' + (canFull ? 'Full screen (⤢) gives the picture more room.' : '') + '</li>' +
      (cfg.help ? '<li>' + mathify(cfg.help) + '</li>' : '') + '</ul>';
    typeset(help);

    /* ---------- graphics context (g = first view; g.view(id) for the others) ---------- */
    var g = views[0];
    g.app = app; g.t = 0; g.dt = 0; g.views = views;
    g.view = function (id) { return viewById[id] || views[0]; };
    app.g = g;
    function resizeCanvas() {
      var changed = false, dpr = Math.min(global.devicePixelRatio || 1, 2.5);
      views.forEach(function (v) {
        var w = Math.max(1, Math.floor(v.holder.clientWidth)), hh = Math.max(1, Math.floor(v.holder.clientHeight));
        if (w !== v.w || hh !== v.h || dpr !== v.dpr) {
          v.w = w; v.h = hh; v.dpr = dpr; v.canvas.width = Math.round(w * dpr); v.canvas.height = Math.round(hh * dpr);
          v.canvas.style.width = w + 'px'; v.canvas.style.height = hh + 'px';
          v.ctx.setTransform(dpr, 0, 0, dpr, 0, 0); changed = true;
        }
      });
      if (changed && cfg.stage && isFn(cfg.stage.resize)) cfg.stage.resize(g, app.state);
      return changed;
    }
    var drawQueued = false;
    app.draw = function () { if (drawQueued) return; drawQueued = true; global.requestAnimationFrame(function () { drawQueued = false; drawNow(); }); };
    function drawNow() {
      if (!cfg.stage) return;
      views.forEach(function (v) { v.ctx.setTransform(v.dpr, 0, 0, v.dpr, 0, 0); });
      try {
        if (isFn(cfg.stage.draw)) cfg.stage.draw(g, app.state);
        views.forEach(function (v) { if (isFn(v.cfg.draw)) v.cfg.draw(v, app.state, app); });
      } catch (e) { reportError(e); }
    }
    views.forEach(function (v) {
      function pointer(ev) {
        var r = v.canvas.getBoundingClientRect(); v.pointer.x = ev.clientX - r.left; v.pointer.y = ev.clientY - r.top;
        if (ev.type === 'pointerdown') { v.pointer.down = true; try { v.canvas.setPointerCapture(ev.pointerId); } catch (e) { /* ignore */ } }
        if (ev.type === 'pointerup' || ev.type === 'pointercancel') v.pointer.down = false;
        v.pointer.inside = ev.type !== 'pointerleave';
        g.pointer = Object.assign({}, v.pointer, { view: v.id });
        // UIEvent.view is a read-only prototype getter (assigning it throws in strict mode): shadow it with an own property
        ev.viewId = v.id; ev.vizView = v;
        try { Object.defineProperty(ev, 'view', { value: v, configurable: true }); } catch (e) { /* ev.vizView still set */ }
        var redraw;
        if (isFn(v.cfg.onPointer)) redraw = v.cfg.onPointer(v, ev, app.state, app);
        else if (cfg.stage && isFn(cfg.stage.onPointer)) redraw = cfg.stage.onPointer(g, ev, app.state);
        else return;
        if (redraw !== false) { updateDynamic(true); app.draw(); }
      }
      ['pointerdown', 'pointermove', 'pointerup', 'pointercancel', 'pointerleave'].forEach(function (t) { v.canvas.addEventListener(t, pointer); });
      if (isFn(v.cfg.onPointer) || (cfg.stage && isFn(cfg.stage.onPointer))) v.canvas.style.touchAction = 'none';
    });

    /* ---------- animation loop + transport (pauses when hidden or scrolled out of view) ---------- */
    var visible = true, last = null, raf = null, holdLeft = 0, lastDyn = 0;
    app.play = function (on) {
      if (!animated) return;
      app.playing = on === undefined ? !app.playing : !!on;
      if (app.playing && tr && app.state[tr.param] >= params[tr.param].max - 1e-12 && tr.end !== 'loop') { app.holding = false; holdLeft = 0; setInternal(tr.param, params[tr.param].min); }
      playBtn.innerHTML = app.playing ? ICON.pause : ICON.play;
      playBtn.setAttribute('aria-pressed', app.playing ? 'true' : 'false');
      last = null; if (app.playing) loop(); else if (raf) { global.cancelAnimationFrame(raf); raf = null; }
    };
    function setInternal(key, v) { app._fromLoop = true; try { app.set(key, v); } finally { app._fromLoop = false; } }
    function loop() {
      if (raf) return;
      raf = global.requestAnimationFrame(function (ts) {
        raf = null;
        if (!app.playing) return;
        if (visible && !doc.hidden) {
          var dt = last === null ? 0 : Math.min(0.05, (ts - last) / 1000); last = ts;
          var speed = app.state.speed !== undefined ? app.state.speed : 1;
          g.dt = dt * speed; g.t += g.dt; app.t = g.t;
          if (tr) {
            var P = params[tr.param], v = app.state[tr.param];
            if (holdLeft > 0) { holdLeft -= dt; if (holdLeft <= 0) { app.holding = false; setInternal(tr.param, P.min); } }
            else {
              v += dt * speed * (tr.rate || (P.max - P.min) / 10);
              if (v >= P.max) {
                if (tr.end === 'loop') v = P.min + (v - P.max);
                else if (tr.end === 'stop') { v = P.max; setInternal(tr.param, v); app.play(false); return; }
                else { v = P.max; holdLeft = tr.hold === undefined ? 2 : tr.hold; app.holding = true; }
              }
              setInternal(tr.param, v);
            }
          }
          if (cfg.stage && isFn(cfg.stage.step)) { try { cfg.stage.step(g, app.state, g.dt); } catch (e) { reportError(e); app.play(false); return; } }
          if (ts - lastDyn > 110) { lastDyn = ts; updateDynamic(false); }
          drawNow();
        } else last = null;
        loop();
      });
    }
    playBtn.onclick = function () { app.play(); };
    if ('IntersectionObserver' in global) new IntersectionObserver(function (en) { visible = en[0].isIntersecting; }).observe(root);

    /* ---------- state changes ---------- */
    var liveQueued = false, lastLive = 0;
    app.set = function (key, value) {
      var changes = {};
      if (typeof key === 'object' && key !== null) changes = key; else changes[key] = value;
      var changed = [];
      Object.keys(changes).forEach(function (k) {
        var p = params[k], v = changes[k];
        if (p && (p.type || 'range') === 'range' && typeof v === 'number') v = clamp(v, p.min, p.max);
        if (app.state[k] !== v) { app.state[k] = v; changed.push(k); }
        if (app.controls[k]) app.controls[k] = app.controls[k].filter(function (c) { return c.el.isConnected; });
        (app.controls[k] || []).forEach(function (c) { c.sync(v); });
      });
      if (!changed.length) return app;
      if (tr && !app._fromLoop && changed.indexOf(tr.param) >= 0) { app.holding = false; holdLeft = 0; }
      changed.forEach(function (k) { if (isFn(cfg.onChange)) { try { cfg.onChange(app.state, k, app); } catch (e) { reportError(e); } } });
      var onlyTime = tr && changed.length === 1 && changed[0] === tr.param;
      updateReadouts();
      if (!app._fromLoop) updateDynamic(true, !onlyTime);
      queueLive(); if (!app._fromLoop) app.draw(); saveHash();
      return app;
    };
    app.reset = function () { app.holding = false; holdLeft = 0; app.set(Object.assign({}, app.defaults)); g.t = 0; app.t = 0; if (cfg.stage && isFn(cfg.stage.reset)) cfg.stage.reset(g, app.state); app.draw(); };
    resetBtn.onclick = function () { app.reset(); };
    function updateReadouts() {
      readoutEls.forEach(function (o) {
        var v; try { v = o.r.value(app.state, app); } catch (e) { v = NaN; }
        o.v.textContent = o.r.fmt ? o.r.fmt(v, app.state) : fmt(v, { sig: o.r.sig || 3, unit: o.r.unit });
        var tone = isFn(o.r.tone) ? o.r.tone(v, app.state) : ''; o.el.classList.toggle('pos', tone === 'pos'); o.el.classList.toggle('neg', tone === 'neg');
      });
    }
    /** Everything derived from the state that is not a canvas: status, presets, view titles, terms, inspector, notes,
        step-by-step working, code values. `structural` rebuilds HTML that depends on parameters (not on time). */
    function updateDynamic(force, structural) {
      if (statusEl) {
        var st; try { st = cfg.status(app.state, app); } catch (e) { st = null; }
        if (typeof st === 'string') st = { text: st };
        statusEl.hidden = !st || !st.text;
        if (st && st.text) { var sh = mathify(st.text); if (statusEl._h !== sh) { statusEl._h = sh; statusEl.innerHTML = sh; typeset(statusEl); } statusEl.setAttribute('data-tone', st.tone || 'info'); }
      }
      if (presetsEl) Array.prototype.forEach.call(presetsEl.querySelectorAll('.viz-preset'), function (b, i) {
        var set = cfg.presets[i].set || {}, on = Object.keys(set).length > 0 && Object.keys(set).every(function (k) { var a = app.state[k], bb = set[k]; return typeof a === 'number' && typeof bb === 'number' ? Math.abs(a - bb) < 1e-9 * Math.max(1, Math.abs(bb)) : a === bb; });
        b.setAttribute('aria-pressed', on ? 'true' : 'false');
      });
      views.forEach(function (v) {
        if (!v.titleEl) return;
        var t = isFn(v.cfg.title) ? v.cfg.title(app.state, app) : v.cfg.title;
        if (v.titleEl._h !== t) { v.titleEl._h = t; v.titleEl.innerHTML = mathify(t); typeset(v.titleEl); }
      });
      updateTerms();
      updateInspectNotes();
      updateCalc(structural !== false);
      updateCodeLive();
      if (app.tab === 'derive') updateDerDynamic();
      if (tourEls.card && force !== false) refreshStepExtras();
    }
    function refreshStepExtras() {
      var ins = tourEls.card.querySelector('.viz-step-inspect'); if (ins && isFn(cfg.inspect)) { var hh = cfg.inspect(app.state, app) || '<span class="viz-muted">' + mathify(cfg.inspectHint || 'Click the picture to inspect a point.') + '</span>'; if (ins._h !== hh) { ins._h = hh; ins.innerHTML = mathify(hh); typeset(ins); } }
      var nt = tourEls.card.querySelector('.viz-step-notes'); if (nt && isFn(cfg.notes)) { var nh = cfg.notes(app.state, app); if (nt._h !== nh) { nt._h = nh; nt.innerHTML = mathify(nh); typeset(nt); } }
    }
    function queueLive() {
      if (liveQueued) return; liveQueued = true;
      var wait = app.playing ? Math.max(0, 180 - (now() - lastLive)) : 0;
      setTimeout(function () {
        global.requestAnimationFrame(function () {
          liveQueued = false; lastLive = now();
          Object.keys(eqEls).forEach(function (id) { var e = eqEls[id]; if (!e || !e.live || !e.cfg.live) return; if (app.tab !== 'equations') return; try { renderTex(e.live, e.cfg.live(app.state, app), true); } catch (err) { e.live.textContent = '—'; } scaleEq(e.live); });
          var tourEq = tourEls.card && tourEls.card.querySelectorAll('[data-live-eq]');
          if (tourEq && app.tab === 'tour') for (var i = 0; i < tourEq.length; i++) { var ec = eqEls[tourEq[i].getAttribute('data-live-eq')]; if (ec && ec.cfg.live) { renderTex(tourEq[i], ec.cfg.live(app.state, app), true); scaleEq(tourEq[i]); } }
        });
      }, wait);
    }
    function scaleEq(inner) {
      if (!inner || !inner.parentNode) return;
      inner.style.transform = ''; var avail = inner.parentNode.clientWidth - 2, need = inner.scrollWidth;
      if (avail > 0 && need > avail) inner.style.transform = 'scale(' + Math.max(0.78, avail / need) + ')';
      inner.parentNode.classList.toggle('viz-eq-too-wide', avail > 0 && need * 0.78 > avail);
    }
    function scaleAllEq() { Array.prototype.forEach.call(root.querySelectorAll('.viz-eq-inner'), scaleEq); }

    /* ---------- tabs and tour ---------- */
    app.setTab = function (id) {
      if (app.tabs.indexOf(id) < 0) id = app.tabs[0];
      app.tab = id; root.setAttribute('data-tab', id);
      Object.keys(tabBtns).forEach(function (k) { tabBtns[k].setAttribute('aria-selected', k === id ? 'true' : 'false'); });
      Object.keys(panels).forEach(function (k) { panels[k].classList.toggle('is-active', k === id); });
      clearHighlights();
      if (id === 'tour') applyStepHighlights();
      queueLive();
      app.fit(); saveHash();
      return app;
    };
    function clearHighlights() { app._hl.forEach(function (el) { el.classList.remove('viz-hl'); }); app._hl = []; }
    function applyStepHighlights() {
      var s = cfg.tour && cfg.tour[app.step]; if (!s || !s.highlight) return;
      [].concat(s.highlight).forEach(function (key) { Array.prototype.forEach.call(root.querySelectorAll('[data-viz-key="' + key + '"]'), function (el) { el.classList.add('viz-hl'); app._hl.push(el); }); });
    }
    app.goStep = function (i) {
      if (!cfg.tour || !cfg.tour.length) return app;
      i = clamp(i, 0, cfg.tour.length - 1); app.step = i;
      var s = cfg.tour[i];
      if (app.tab !== 'tour') app.setTab('tour');
      if (s.set) app.set(s.set);
      if (s.play !== undefined) app.play(s.play);
      tourEls.count.textContent = 'Step ' + (i + 1) + ' of ' + cfg.tour.length;
      tourEls.title.innerHTML = mathify(s.title);
      tourEls.text.innerHTML = mathify(s.text);
      var X = tourEls.extras; X.innerHTML = '';
      (s.controls || []).forEach(function (k) { X.appendChild(control(k, { keep: true })); });
      if (s.readouts && s.readouts.length) { var ro = readoutsGrid(s.readouts); if (ro) X.appendChild(ro); }
      if (s.terms && cfg.terms) { var tb = termsBlock(); if (tb) X.appendChild(tb); }
      if (s.eq && eqEls[s.eq]) {
        var e = eqEls[s.eq].cfg;
        var box = h('div', { class: 'viz-eq viz-step-eq' }, h('div', { class: 'viz-eq-head' }, h('span', { html: mathify(e.title) }), e.ref ? h('span', { class: 'viz-eq-ref', text: e.ref }) : null));
        var b1 = h('div', { class: 'viz-eq-body' }, h('div', { class: 'viz-eq-inner' })); renderTex(b1.firstChild, e.tex, true); box.appendChild(b1);
        if (e.live) box.appendChild(h('div', { class: 'viz-eq-live' }, h('div', { class: 'viz-eq-inner', 'data-live-eq': e.id })));
        X.appendChild(box);
      }
      if (s.code && codeCfg[s.code.id || s.code]) {
        var sc = typeof s.code === 'string' ? { id: s.code } : s.code, cc = codeCfg[sc.id], ln = sc.lines || [1, 1], ctxN = sc.context === undefined ? 1 : sc.context;
        X.appendChild(h('div', { class: 'viz-code-card viz-step-code' }, codeBlock(cc, ln[0] - 1 - ctxN, ln[1] - 1 + ctxN, ln)));
      }
      if (s.derive && cfg.derivations) {
        // a walkthrough step can quote one step of a derivation and link to the full one
        var sd = typeof s.derive === 'string' ? { id: s.derive } : s.derive, di = -1;
        cfg.derivations.forEach(function (dd, k) { if (dd.id === sd.id) di = k; });
        if (di >= 0) {
          var dd = cfg.derivations[di], dk = sd.step ? clamp(sd.step, 1, dd.steps.length) : 0, dst = dk ? dd.steps[dk - 1] : null;
          var mini = h('div', { class: 'viz-der-mini' },
            h('div', { class: 'viz-der-mini-head' }, h('span', { html: mathify('Derivation · ' + (dst ? 'step ' + dk + (dst.did ? ': ' + dst.did : '') : dd.title)) }),
              h('button', { class: 'viz-btn', type: 'button', text: 'All steps →', onclick: function () { app.goDerive(di, dk); } })));
          var tex = dst ? dst.tex : (dd.result && dd.result.tex) || (typeof dd.start === 'string' ? dd.start : dd.start && dd.start.tex);
          if (tex) mini.appendChild(derLine(tex, 'viz-der-cur'));
          if (dst && dst.why) mini.appendChild(derBlock('viz-der-why', 'Why', dst.why));
          X.appendChild(mini);
        }
      }
      if (s.inspect && isFn(cfg.inspect)) X.appendChild(h('div', { class: 'viz-inspect viz-step-inspect' }));
      if (s.notes && isFn(cfg.notes)) X.appendChild(h('div', { class: 'viz-notes viz-step-notes' }));
      if (s.callout) X.appendChild(h('div', { class: 'viz-callout ' + (s.callout.kind || 'try') }, h('b', { class: 'k', text: CALLOUT_LABEL[s.callout.kind || 'try'] }), h('span', { html: mathify(s.callout.html) })));
      tourEls.pager.setItems([tourEls.title, tourEls.text].concat(Array.prototype.slice.call(X.children)));
      typeset(tourEls.card);
      Array.prototype.forEach.call(tourEls.chips.children, function (c, k) { c.classList.toggle('done', k < i); if (k === i) c.setAttribute('aria-current', 'step'); else c.removeAttribute('aria-current'); });
      tourEls.prev.disabled = i === 0;
      tourEls.next.textContent = i === cfg.tour.length - 1 ? 'Explore →' : 'Next →';
      tourEls.next.onclick = i === cfg.tour.length - 1 ? function () { app.setTab('explore'); } : function () { app.goStep(app.step + 1); };
      updateReadouts(); updateDynamic(true, true); queueLive();
      clearHighlights(); applyStepHighlights();
      if (isFn(s.enter)) { try { s.enter(app); } catch (err) { reportError(err); } }
      app.draw(); app.fit(); saveHash();
      return app;
    };

    function saveHash() {
      if (!app._ready || app._fromLoop) return;
      var o = { tab: app.tab };
      if (app.tab === 'tour') o.step = app.step + 1;
      if (app.tab === 'derive') { o.d = der.i + 1; o.ds = der.s; }
      writeHash(o);
    }

    /* ---------- keyboard ---------- */
    doc.addEventListener('keydown', function (ev) {
      var tag = (ev.target && ev.target.tagName) || '';
      if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') return;
      if (ev.key === 'ArrowRight' && app.tab === 'tour') { app.goStep(app.step + 1); ev.preventDefault(); }
      else if (ev.key === 'ArrowLeft' && app.tab === 'tour') { app.goStep(app.step - 1); ev.preventDefault(); }
      else if (ev.key === 'ArrowRight' && app.tab === 'derive') { app.goDerive(der.i, der.s + 1); ev.preventDefault(); }
      else if (ev.key === 'ArrowLeft' && app.tab === 'derive') { app.goDerive(der.i, der.s - 1); ev.preventDefault(); }
      else if (ev.key === ' ' && animated) { app.play(); ev.preventDefault(); }
      else if (ev.key === '?') { help.classList.toggle('open'); }
      else if (ev.key === 'Escape') { help.classList.remove('open'); }
    });
    helpBtn.onclick = function () { help.classList.toggle('open'); };
    themeBtn.onclick = function () {
      var t = doc.documentElement.getAttribute('data-theme') === 'dark' ? 'light' : 'dark';
      doc.documentElement.setAttribute('data-theme', t); if (STORE) STORE.setItem('viz-theme', t); app.draw();
    };
    fullBtn.onclick = function () {
      var el = doc.documentElement;
      if (doc.fullscreenElement || doc.webkitFullscreenElement) (doc.exitFullscreen || doc.webkitExitFullscreen).call(doc);
      else (el.requestFullscreen || el.webkitRequestFullscreen).call(el);
    };

    /* ---------- FIT: the no-scroll guarantee ---------- */
    function overflowing() {
      var bad = [];
      var de = doc.documentElement, rr = root.getBoundingClientRect();
      if (rr.bottom > global.innerHeight + 1 || rr.right > global.innerWidth + 1) bad.push({ el: 'app-exceeds-window', h: Math.round(rr.bottom), H: global.innerHeight, w: Math.round(rr.right), W: global.innerWidth });
      if (rr.height < global.innerHeight - 2) bad.push({ el: 'app-does-not-fill-window', h: Math.round(rr.height), H: global.innerHeight });
      if (de.scrollHeight > global.innerHeight + 1 || de.scrollWidth > global.innerWidth + 1) bad.push({ el: 'document', h: de.scrollHeight, H: global.innerHeight, w: de.scrollWidth, W: global.innerWidth });
      var els = root.querySelectorAll('.viz-header, .viz-panel.is-active, .viz-panel.is-active .viz-card, .viz-panel.is-active .viz-pages, .viz-stage, .viz-stage-top, .viz-transport');
      for (var i = 0; i < els.length; i++) {
        var e = els[i]; if (!e.offsetParent && e.className.indexOf('viz-stage') < 0) continue;
        if (e.scrollHeight > e.clientHeight + 1 || e.scrollWidth > e.clientWidth + 1) bad.push({ el: (e.className || e.tagName).split(' ').slice(0, 2).join('.'), h: e.scrollHeight, H: e.clientHeight, w: e.scrollWidth, W: e.clientWidth });
      }
      views.forEach(function (v) { if (v.el.offsetParent && v.holder.clientHeight < 60) bad.push({ el: 'view-too-small:' + v.id, h: v.holder.clientHeight }); });
      Array.prototype.forEach.call(root.querySelectorAll('.viz-panel.is-active .viz-eq-too-wide'), function (e) { bad.push({ el: 'equation-too-wide', w: e.firstChild ? e.firstChild.scrollWidth : 0, W: e.clientWidth }); });
      return bad;
    }
    app.fit = function () {
      var W = global.innerWidth, H = global.innerHeight;
      var layout = (W < 700 && H >= W * 0.9) || W < 520 ? 'portrait' : (H < 520 ? 'landscape' : 'wide');
      root.setAttribute('data-layout', layout);
      root.classList.remove('short-tabs', 'tabs-row');
      root.classList.toggle('compact-title', W < 420);
      // header: long tab labels inline → short labels inline → tabs on their own row (long, then short labels)
      var titleBox = header.firstChild;
      var cramped = function () { return nav.scrollWidth > nav.clientWidth + 1 || header.scrollWidth > header.clientWidth + 1 || (layout !== 'portrait' && titleBox.clientWidth < 200); };
      if (cramped()) root.classList.add('short-tabs');
      if (cramped() && layout === 'wide') { root.classList.remove('short-tabs'); root.classList.add('tabs-row'); if (cramped()) root.classList.add('short-tabs'); }
      var bad = [];
      for (var d = 0; d <= 3; d++) {
        root.setAttribute('data-dense', String(d));
        pagers.forEach(function (p) { if (p.card.offsetParent) p.layout(); });
        scaleAllEq();
        bad = overflowing();
        if (!bad.length) break;
      }
      resizeCanvas(); drawNow();
      app.lastAudit = { layout: layout, dense: Number(root.getAttribute('data-dense')), overflow: bad, W: W, H: H, tab: app.tab, step: app.step };
      if (bad.length) console.warn('VIZ-OVERFLOW', JSON.stringify(app.lastAudit));
      return app.lastAudit;
    };
    var fitQueued = false;
    function queueFit() { if (fitQueued) return; fitQueued = true; global.requestAnimationFrame(function () { fitQueued = false; app.fit(); }); }
    global.addEventListener('resize', queueFit);
    if ('ResizeObserver' in global) { var ro2 = new ResizeObserver(function () { if (resizeCanvas()) drawNow(); }); views.forEach(function (v) { ro2.observe(v.holder); }); }

    /* ---------- errors ---------- */
    function reportError(e) {
      console.error('VIZ-ERROR', e && e.stack ? e.stack : e);
      if (!root.querySelector('.viz-error')) stageEl.appendChild(h('div', { class: 'viz-error', text: 'Something went wrong drawing this picture: ' + (e && e.message ? e.message : e) }));
    }
    app.reportError = reportError;

    /* ---------- boot ---------- */
    var hash = parseHash();
    Object.keys(params).forEach(function (k) { if (hash[k] !== undefined && hash[k] !== '') { var p = params[k]; var v = (p.type || 'range') === 'range' ? Number(hash[k]) : (p.type === 'toggle' ? hash[k] === '1' || hash[k] === 'true' : hash[k]); if (!(typeof v === 'number' && isNaN(v))) app.state[k] = v; } });
    try { resizeCanvas(); if (cfg.stage && isFn(cfg.stage.setup)) cfg.stage.setup(g, app.state); } catch (e) { reportError(e); }
    typeset(root);
    updateReadouts(); updateDynamic(true, true);
    if (hash.tab === 'calc') hash.tab = 'explain';
    var firstTab = hash.tab && app.tabs.indexOf(hash.tab) >= 0 ? hash.tab : (cfg.startTab || app.tabs[0]);
    if (cfg.tour && cfg.tour.length) app.goStep(hash.step ? Number(hash.step) - 1 : 0);
    if (der.els) app.goDerive(hash.d ? Number(hash.d) - 1 : 0, hash.ds ? Number(hash.ds) : 0, { quiet: firstTab !== 'derive' });
    app.setTab(firstTab);
    queueLive();
    if (animated && cfg.autoplay !== false) app.play(true); else app.draw();

    app.ready = loadKatex().then(function (ok) {
      app.katex = ok; retypesetPending(); queueLive();
      return new Promise(function (res) { global.requestAnimationFrame(function () { global.requestAnimationFrame(function () { app._ready = true; app.fit(); res(app.lastAudit); }); }); });
    });

    /* ---------- audit + self-test API (used by tools/shot.py) ---------- */
    app.selftest = function () {
      var out = [];
      if (!isFn(cfg.selftest)) return out;
      var rows; try { rows = cfg.selftest(app); } catch (e) { return [{ name: 'selftest() threw', ok: false, error: String(e) }]; }
      (rows || []).forEach(function (r) {
        var row = { name: r.name, js: r.js, py: r.py || null, rtol: r.rtol === undefined ? 1e-6 : r.rtol, atol: r.atol === undefined ? 0 : r.atol };
        if (r.expect !== undefined) { var dd = Math.abs(r.js - r.expect); row.expect = r.expect; row.ok = dd <= row.atol + row.rtol * Math.abs(r.expect); }
        out.push(row);
      });
      return out;
    };
    global.VIZ = {
      app: app, ready: app.ready, tabs: app.tabs, steps: (cfg.tour || []).length,
      derivations: (cfg.derivations || []).map(function (d) { return { id: d.id, title: d.title, ref: d.ref || '', pages: d.steps.length + 2, steps: d.steps.length }; }),
      goDerive: function (i, s) { app.goDerive(i, s); return app.fit(); },
      explainStats: function () { var b = calcPager ? calcPager.body : null; return b ? { sections: b.querySelectorAll('.viz-work-h').length, interpret: b.querySelectorAll('.viz-work-interp').length, boxes: b.querySelectorAll('.viz-work-res').length, live: b.querySelectorAll('[data-live]').length } : null; },
      features: { views: views.length, explain: !!cfg.explain, calc: !!cfg.explain, derivations: (cfg.derivations || []).length, code: !!(cfg.code && cfg.code.length), terms: !!cfg.terms, inspect: isFn(cfg.inspect), notes: isFn(cfg.notes), presets: (cfg.presets || []).length, transport: !!tr, modes: !!cfg.modes, status: !!cfg.status },
      audit: function () { return app.fit(); },
      selftest: function () { return app.selftest(); },
      setTab: function (id) { app.setTab(id); return app.fit(); },
      goStep: function (i) { app.goStep(i); return app.fit(); },
      play: function (on) { app.play(on); return app.playing; },
      set: function (o) { app.set(o); return app.fit(); },
      meta: function () { var m = {}; Array.prototype.forEach.call(doc.querySelectorAll('meta[name^="viz:"]'), function (x) { m[x.getAttribute('name').slice(4)] = x.getAttribute('content'); }); return m; }
    };
    return app;
  };

  global.addEventListener('error', function (e) { console.error('VIZ-ERROR', e.message, e.filename + ':' + e.lineno); });
  global.Viz = Viz;
})(window);
/* VIZ_LIB_JS:END */
