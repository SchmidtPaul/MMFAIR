# MMFAIR

Deprecated GitHub Pages site. New content lives at **biomathcontent.netlify.app**.

## Structure

- HTML files are in the repo **root** (no `docs/` subfolder)
- R Markdown site (`_site.yml`), rendered HTML is committed directly

## Deprecation Banner

- Every HTML file contains a deprecation banner (`<!-- DEPRECATION-BANNER-START -->` block)
- Banner `z-index: 999` — must stay **below** the Bootstrap navbar (`z-index: 1030`) so nav dropdowns remain clickable
- Do NOT set banner z-index >= 1000
