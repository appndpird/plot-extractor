# Auto-detect plots — evaluation (26–27 Nov 2026)

`plot_detect.py` (tab 1a, **Auto-detect plots**) scored against three hand-checked grids.
Corners were simulated as hand clicks: the reference grid's outer corners plus 25 cm random
error. Scripts: `eval_detect.py` (scoring), `overlay.py` (pictures), `gui_test_detect.py`
(the full 1a → 1b loop in the real app). Console output: `eval_console.txt`.

| Trial | Stage | Plots found | Across rows (median / p90) | Along rows (median / p90) | Size vs reference | Median IoU | Flagged | Time |
|---|---|---|---|---|---|---|---|---|
| York 2026 calibration (12 x 21) | mid-season wheat | 252 / 252 | 2.3 / 5.7 cm | 8.6 / 15.7 cm | width −2 cm, length +8 cm | 0.92 | 0 | 13 s |
| NUE Northam 2026-07-03 (4 x 32) | early emergence, ~1 % green, cloud shadow | 128 / 128 (123 matched) | 4.4 cm | 4.2 / 8.6 cm | width −5 cm, length −15 cm | 0.91 | 54 | 12 s |
| York 2025 calibration (6 x 6) | barely emerged, stubble | 36 / 36 | 6.5 cm | 41 cm (one plot end) | length +81 cm | 0.75 | 16 | 2 s |
| York 2026, **no corners** | | 252 / 252 + 48 neighbouring | 2.2 cm | 8.8 cm | same | 0.92 | 46 | 31 s |

Notes on the references:

- **Northam**: 5 plots "unmatched" and the across-row p90 (58 cm) come from part of range 4
  where the reference (fitted by hand in Oct 2026) is ~0.6 m off the drilled furrow groups.
  The detected plots each hold one group of 7 furrows (`figures/northam_reference_off_plot4003.jpg`).
- **York 2025**: both grids start each plot at the wheel track; at the other end the
  reference stops ~1 m earlier, before a band where the drill rows run on through stubble
  (`figures/york2025_headland_band.jpg`). Across the rows the reference sits ~9 cm off the
  rows (known since Oct 2026).
- **York 2026** without corners also finds 4 columns of real plots outside the reference;
  with corners only the trial is used.

## Method (what made it work)

1. **Direction** from the crop rows (power at seeder spacings, summed over short segments).
2. **Columns** (plots across the rows): the trial-wide across profile with rows smoothed out;
   a comb at the plot pitch finds a narrow repeating feature of *either sign* — a bare dip
   (York) or a green strip (Northam).
3. **Ranges** (along the rows): per column, a repeating *plot box* (length L, pitch P),
   scored inside minus two narrow strips outside its ends; one straight alley line through
   all columns (staggered alleys allowed); the shortest box within 10 % of the best score;
   the pitch, half or a third of it (shortest within 75 % of the best). Three candidate
   signals, the clearest wins: **row amplitude** (rows parallel to the plot only — weeds,
   stubble, headland turns drop out), **greenness**, **soil darkness** (Northam: the drill
   rows run straight through the alleys; only the lighter, wheel-tracked soil shows them).
4. **Width** = crop rows × row spacing (rows counted between the gaps on the trial-wide
   profile), the plot definition of 1a / 1c.
5. **Each plot** nudged on its own, checked against its neighbours (same column, same
   range); confidence 0–1 with the reason (side / end unclear, little crop, out of line).

Approaches tried and dropped: one threshold on greenness (breaks thin plots), trial-wide
along profile (staggered alleys blur), narrow-dip alley comb (wide flat alleys; green
strips in alleys), half-height crop width (outer rows read narrow), per-range across
offsets (noise at 1 % cover).

## GUI tests

`gui_test_detect.py northam`: 21 / 21 (Auto-detect, fields, corners, Save & use, confidence
CSV, rows file, 1b basemap = 1a ortho, Next to check / Shift+N / Checked, list carried to
the edited grid and kept on reopening). `york2026`: 16 / 16.
