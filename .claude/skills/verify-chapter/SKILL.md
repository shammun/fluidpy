---
name: verify-chapter
description: Re-run only the verification phase for a chapter (tests on the V1-V7 evidence ladder proportional to curation tiers, benchmark references, report), with the implementer fix loop, then re-check explainer parity and notebook execution because physics may have changed. Usage - /verify-chapter 7
disable-model-invocation: true
---

# /verify-chapter N
Same as phase 5 of `/do-chapter` (read that skill for the brief header): brief `math-verifier`; on FAIL loop with
`concept-implementer` (max 3; physics fixes only, public signatures stable); update `progress.json → verify`; commit
`chNN: verify — …`.
Because code may have changed, then run `.venv/Scripts/python.exe tools/shot.py --chapter chNN --quick` (selftest parity)
and `.venv/Scripts/python.exe tools/run_notebook.py chNN`; report both. If the chapter was already published, suggest
`/publish-chapter N`.
Print: the Verdict, counts per validation label, observed convergence orders, Open items.
