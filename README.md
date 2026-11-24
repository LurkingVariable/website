# aaronspivak.com

Personal website of Aaron Spivak, PhD, MBA — a biotechnology investor and computational biologist. The live site is served from this repository by GitHub Pages at [aaronspivak.com](https://aaronspivak.com/).

## What's here

- `index.html` — the homepage (about, professional background, speaking highlights, selected publications, and contact links). CSS and Schema.org Person/ProfilePage structured data are embedded in the page. The homepage has no build step and no external font or JavaScript dependencies.
- `scertf/` — [ScerTF](https://aaronspivak.com/scertf/), a client-side database of benchmarked position weight matrices for *Saccharomyces* transcription factors. Search, browse, details, downloads, and references all live under this path. `ScerTF/` is a small redirect stub so the capitalized URL still reaches the app.
- `speaking/` — an archive of selected conference appearances, with a stable summary page for each event.
- `favicon.svg` and `og-image.png` — site icon and social-sharing card (the card is not shown on the page itself).
- `robots.txt`, `sitemap.xml`, `sitemap-pages.xml`, and `scertf/sitemap.xml` — crawler directives. The root `sitemap.xml` is a sitemap index that points to the pages sitemap (homepage and speaking archive) and the ScerTF sitemap (the 316 URLs in the database app).
- `CNAME` — GitHub Pages custom-domain association for `aaronspivak.com`.
- `googleefd8f934feb68710.html` — Google Search Console verification file.

## Deployment

This is a static GitHub Pages site published from the `main` branch of this repository to the custom domain `aaronspivak.com`. There is no build pipeline: the homepage, supporting files, and the ScerTF client-side app are served as committed.

To preview locally, serve the repository root over HTTP and open `/` and `/scertf/`. For example:

```bash
python3 -m http.server 8000
```

Then visit `http://127.0.0.1:8000/` and `http://127.0.0.1:8000/scertf/`.
