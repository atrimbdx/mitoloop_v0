#!/usr/bin/env python3
"""
mitomap_smoke.py - load app/mitomap.html in headless Chromium, load a
mitomap.json and a score table, and check that the page draws without
errors. Used by tests/test_mitomap.sh when Playwright is installed.

Usage: python3 mitomap_smoke.py APP_HTML MITOMAP_JSON SCORES_TSV [CHROME]
Prints one JSON line with the checked quantities; exits 1 on failure.
"""

import asyncio
import json
import sys

from playwright.async_api import async_playwright


async def main(app, data, scores, chrome):
    async with async_playwright() as p:
        kw = {"executable_path": chrome} if chrome else {}
        b = await p.chromium.launch(**kw)
        pg = await b.new_page(viewport={"width": 1280, "height": 800})
        errors = []
        pg.on("pageerror", lambda e: errors.append(str(e)))
        pg.on("console", lambda m: errors.append(m.text)
              if m.type == "error" and "fonts.g" not in m.text
              and "ERR_" not in m.text else None)
        await pg.goto("file://" + app)
        await pg.wait_for_function("typeof S !== 'undefined' && S.data && S.inst.length > 0",
                                   timeout=15000)
        demo = await pg.evaluate("S.data.genes.length")
        await pg.set_input_files("#dataFile", data)
        await pg.wait_for_function("S.data && !S.data.meta.demo")
        await pg.set_input_files("#scoreFile", scores)
        await pg.wait_for_function("S.score && S.score.size > 0")
        await pg.wait_for_timeout(300)
        out = await pg.evaluate("""() => ({
            genes: S.data.genes.length, instances: S.inst.length,
            links: S.links.length, scored: S.score.size,
            scale: S.scale.type, drawn: S.drawn.length,
            labels: S.drawn.filter(d => d.label).length,
            unmatched: document.getElementById('unmatchedSum').textContent,
            outside: S.inst.filter(i => i.x < 0 || i.x > W || i.y < 0 || i.y > H).length,
            multi: S.byGene.filter(x => x.length > 1).length})""")
        await pg.fill("#searchInput", "MC5")
        await pg.click("#searchForm button")
        await pg.wait_for_timeout(200)
        out["selected"] = await pg.evaluate("S.sel >= 0 ? S.data.genes[S.sel].s : ''")
        out["info"] = (await pg.inner_text("#infoBody"))[:60]
        out["demo_genes"] = demo
        out["errors"] = errors
        await b.close()
        return out


if __name__ == "__main__":
    a = sys.argv
    res = asyncio.run(main(a[1], a[2], a[3], a[4] if len(a) > 4 else None))
    print(json.dumps(res, ensure_ascii=False))
    sys.exit(1 if res["errors"] else 0)
