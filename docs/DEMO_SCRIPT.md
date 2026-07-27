# Plot Extractor — Demo Recording Script

A ready-to-record walkthrough. Below are two ways to capture it, then the click-by-click beats.
Target length **~3–4 minutes**. Each beat lists **what to click** and a **narration line**.

## How to record

**Option A — Steps Recorder (screenshot-per-click document).** Produces a `.zip` with an
`.mht` slideshow that captures a screenshot and description on every click — ideal as a printable
how-to. It is *not* an MP4 video.
- **Easiest:** double-click **`Record_Demo.bat`** (in this `docs\` folder). It opens Steps Recorder
  and the app together. In Steps Recorder click **Start Record**, do the beats below, then
  **Stop Record** and **Save**.
- **Hands-free:** double-click **`Record_Demo_Auto.bat`** to start recording + open the app, do the
  demo, then double-click **`Stop_Demo.bat`**. The result lands on your Desktop as
  `PlotExtractor_Demo.zip`.
- *Note:* Steps Recorder is deprecated on Windows 11 24H2 but `psr.exe` is still present on this
  machine, so it works.

**Option B — Xbox Game Bar (actual MP4 video).** Press `Win`+`G`, click the **Record** (●) button
(or `Win`+`Alt`+`R` to start/stop). Records a real video to `Videos\Captures`. Use this if you need
a true video file. OBS Studio works too.

> Tip: whichever you use, follow the same beats below. For Game Bar, read the **narration** aloud.

---

---

### 0. Intro (title card) — 10 s
- Show the presentation title slide (`PlotExtractor_Presentation.html`, or the `.pptx`).
- **Say:** "This is Plot Extractor — it turns a drone survey into one clean image per plot,
  ready for emergence counting."

### 1. Launch — 10 s
- Double-click **`PlotExtractor.bat`**. The window opens.
- **Say:** "Launch it from PlotExtractor.bat. It remembers your last settings."

### 2. Source & project — 20 s
- **Source** dropdown → keep **Metashape project**.
- **Metashape project (.psx)** → **Browse…** → pick the `.psx`.
- **Chunk / sub-project** → open the dropdown → choose **RootPheno_Day1**.
- **Say:** "Point it at the Metashape project, and pick the right chunk — this is how you choose
  Day 1 versus Day 2 in the same project."

### 3. Plots & imagery — 20 s
- **Plot boundaries** → **Browse…** → select the plot shapefile.
- Note **Plot ID attribute** auto-fills (`Plot_ID`); leave **CRS** on **auto**.
- **Raw images folder** → **Browse…** → the `JPEG\Day1` folder.
- **Output folder** → **Browse…** → a new folder.
- **Say:** "Load the plot boundaries, point at the raw photos, and choose an output folder."

### 4. Options — 25 s
- **Format** → **JPEG (small files)**, **Quality** 100.
- Tick **Raw image crop per plot**; **style** → **Rectified rectangle**.
- Tick **Orthomosaic crop per plot**.
- Tick **Convert raw crops to sRGB**; set **CPU cores** (e.g. 16).
- Tick **Auto-fit crop to detected rows**; **rows/plot** = 7; **mode** = **width_shift**.
- **Say:** "Choose your outputs and the crop style. The important one is Auto-fit — set it to seven
  rows, width-shift mode. That snaps every crop to the planted rows."

### 5. (Optional) Edit plots — 30 s
- Click **Edit plots…**.
- Click **Load ortho basemap…** → pick the orthomosaic.
- **Ctrl+A** (Select all) → **Shift+←/→** to slide plots onto the rows → **]** to widen slightly.
- Watch the **width × length** read-out in the status bar.
- Click **Save & use**.
- **Say:** "If the grid is off, open the editor, load the ortho as a basemap, and nudge every plot
  onto the rows — the status bar shows the live width and length. Save and use writes an edited
  shapefile the extractor picks up."

### 6. Test run — 25 s
- Set **Limit** = 10.
- Click **Run extraction**. Point at the **log** streaming each plot.
- **Say:** "Always test with a limit of ten first. The log streams each plot — coverage, frames
  used, and the fitted area."

### 7. Review — 20 s
- Click **Open output folder**. Open a couple of crops from **RawCrops\**.
- Open **manifest.csv**; point at `fit_rows`, `fit_area_m2`, `frames_used`, `note`.
- **Say:** "Every crop is one plot, snapped to seven rows. The manifest records coverage, the fitted
  area in square metres, and flags any plot rebuilt from the ortho."

### 8. Full run + close — 15 s
- Set **Limit** back to **0**, click **Run extraction** (or just say it).
- **Say:** "Set the limit back to zero for the full trial. A full raw-plus-ortho run self-corrects
  seams and edge plots automatically — no post-processing. That's Plot Extractor."

---

**Tips for a clean recording**
- Record at 1080p; zoom the app window so field labels are readable.
- Pause ~1 s after each click so viewers can follow the cursor.
- Use the presentation deck for the intro/outro cards.
