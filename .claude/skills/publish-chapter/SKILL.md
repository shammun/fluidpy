---
name: publish-chapter
description: Publish (or re-publish) one chapter - embed check, execute, write the executed notebook, the Colab twin and the full-window-explainer web page, regenerate the site index and explainer gallery, audit the page in a browser, public-repo check, commit, inspect origin, push and confirm the GitHub Pages URL. Usage - /publish-chapter 7   ·   /publish-chapter 7 --no-push   ·   /publish-chapter site (index + gallery only)
disable-model-invocation: true
---

# /publish-chapter N [--no-push] | /publish-chapter site
- `site`: `.venv/Scripts/python.exe tools/publish_notebook.py --site-only`, `tools/check_public.py`, commit
  `site: index and gallery`, push (unless `--no-push`).
- Chapter: brief the `site-publisher` agent (brief header from `/do-chapter`); then do the orchestrator's **Publish
  gate** from `/do-chapter` phase 10 (stage, `check_public.py`, commit `chNN: publish — …`, fetch/inspect origin, push,
  poll the page URL for 200, check one explainer URL and the Colab file URL). Set `publish: pass`.
- If there is no `origin` yet, stop after the local commit and print the one-time commands from GUIDE.md §1.6 for the
  user to run with a leading `!`.
Print: page URL, gallery URL, Colab URL, page/ipynb sizes, runtime, audit verdict.
