const fs = require("fs");
const path = require("path");
const PptxGenJS = require("pptxgenjs");
const { imageSizingContain } = require("./pptxgenjs_helpers_local/image");
const {
  warnIfSlideHasOverlaps,
  warnIfSlideElementsOutOfBounds,
} = require("./pptxgenjs_helpers_local/layout");

const ROOT = "/Users/dingchenchen/lucc";
const TASK_DIR = path.join(ROOT, ".ppt_task_fg_summary");
const OUT_DIR = path.join(TASK_DIR, "out");
fs.mkdirSync(OUT_DIR, { recursive: true });

const FIG_ROOT = path.join(
  ROOT,
  "glmmTMB_occurrence_FG_RobustPCAKmeans_FULL_plusTwoWay"
);
const FIG_EXPL = path.join(FIG_ROOT, "01_Exploration_and_FG_Interpretation");
const FIG_3WAY = path.join(FIG_ROOT, "03_Plots_ThreeWay_ABS_PCT");
const FIG_2WAY = path.join(FIG_ROOT, "04_Plots_TwoWay_PaperGrade");
const LOGO_PKU = path.join(TASK_DIR, "ref_media", "image1.png");
const LOGO_IIASA = path.join(TASK_DIR, "ref_media", "image2.png");

const pptx = new PptxGenJS();
pptx.layout = "LAYOUT_WIDE";
pptx.author = "OpenAI Codex";
pptx.company = "PKU-IIASA";
pptx.subject = "Stable FG scheme summary";
pptx.title = "RQ3 Stable FG Scheme: code, figures, interpretation";
pptx.lang = "en-US";
pptx.theme = {
  headFontFace: "Aptos Display",
  bodyFontFace: "Aptos",
  lang: "en-US",
};

const C = {
  text: "2D2D2D",
  muted: "6A6A6A",
  mint: "D8F1EA",
  mintDark: "8FD5C3",
  pku: "8B1E1E",
  iiasa: "1967B3",
  fg1: "F05A5A",
  fg2: "19B51E",
  fg3: "4E79FF",
  ui2a: "009E73",
  ui2b: "0072B2",
  ui2c: "E69F00",
  ui2d: "D55E00",
  ui2e: "CC79A7",
  line: "D9D9D9",
  white: "FFFFFF",
  pale: "F8FBFA",
};

function addTopChrome(slide, title, subtitle = "") {
  slide.addImage({ path: LOGO_PKU, x: 0.05, y: 0.05, w: 2.25, h: 0.63 });
  slide.addImage({ path: LOGO_IIASA, x: 10.7, y: 0.05, w: 2.15, h: 0.58 });
  slide.addShape(pptx.ShapeType.rect, {
    x: 0,
    y: 0.86,
    w: 13.333,
    h: 0.42,
    line: { color: C.mint, transparency: 100 },
    fill: { color: C.mint },
  });
  slide.addText(title, {
    x: 0.52,
    y: 1.36,
    w: 12.1,
    h: 0.3,
    fontFace: "Aptos Display",
    fontSize: 19,
    bold: true,
    color: C.text,
    margin: 0,
  });
  if (subtitle) {
    slide.addText(subtitle, {
      x: 0.52,
      y: 1.76,
      w: 12.0,
      h: 0.18,
      fontFace: "Aptos",
      fontSize: 9,
      color: C.muted,
      margin: 0,
    });
  }
}

function addFooter(slide, pageLabel) {
  slide.addText(pageLabel, {
    x: 12.25,
    y: 7.08,
    w: 0.55,
    h: 0.18,
    align: "right",
    fontSize: 8,
    color: "7A7A7A",
    margin: 0,
  });
}

function addKeyPanel(slide, cfg) {
  const { x, y, w, h, title, bullets, accent = C.mintDark } = cfg;
  slide.addShape(pptx.ShapeType.roundRect, {
    x,
    y,
    w,
    h,
    rectRadius: 0.08,
    line: { color: "CFE3DD", pt: 1 },
    fill: { color: C.pale },
  });
  slide.addShape(pptx.ShapeType.rect, {
    x,
    y,
    w,
    h: 0.36,
    line: { color: accent, transparency: 100 },
    fill: { color: accent },
  });
  slide.addText(title, {
    x: x + 0.18,
    y: y + 0.08,
    w: w - 0.36,
    h: 0.18,
    fontSize: 11,
    bold: true,
    color: C.text,
    margin: 0,
  });

  let cy = y + 0.52;
  bullets.forEach((b) => {
    slide.addShape(pptx.ShapeType.ellipse, {
      x: x + 0.18,
      y: cy + 0.06,
      w: 0.08,
      h: 0.08,
      line: { color: accent, transparency: 100 },
      fill: { color: accent },
    });
    slide.addText(b, {
      x: x + 0.32,
      y: cy,
      w: w - 0.48,
      h: 0.46,
      fontSize: 9.2,
      color: C.text,
      valign: "top",
      breakLine: false,
      margin: 0,
    });
    cy += 0.55;
  });
}

function addFigureSlide(slide, opts) {
  addTopChrome(slide, opts.title, opts.subtitle || "");
  slide.addImage({
    path: opts.image,
    ...imageSizingContain(opts.image, 0.45, 2.16, 8.35, 4.55),
  });
  addKeyPanel(slide, {
    x: 9.0,
    y: 2.18,
    w: 3.85,
    h: 4.52,
    title: "Key Points",
    bullets: opts.bullets,
    accent: opts.accent || C.mintDark,
  });
  slide.addText(opts.figureLabel, {
    x: 0.48,
    y: 6.83,
    w: 8.3,
    h: 0.18,
    fontSize: 8,
    italic: true,
    color: C.muted,
    margin: 0,
  });
  if (opts.notes) slide.addNotes(opts.notes);
}

function finalizeSlide(slide, pageLabel, opts = {}) {
  addFooter(slide, pageLabel);
  if (!opts.skipOverlapCheck) {
    warnIfSlideHasOverlaps(slide, pptx);
  }
  warnIfSlideElementsOutOfBounds(slide, pptx);
}

function coverSlide() {
  const slide = pptx.addSlide();
  slide.background = { color: C.white };
  slide.addImage({ path: LOGO_PKU, x: 0.02, y: 0.03, w: 2.55, h: 0.7 });
  slide.addImage({ path: LOGO_IIASA, x: 10.7, y: 0.05, w: 2.2, h: 0.6 });
  slide.addShape(pptx.ShapeType.rect, {
    x: 0,
    y: 1.87,
    w: 13.333,
    h: 2.22,
    line: { color: C.mint, transparency: 100 },
    fill: { color: C.mint },
  });
  slide.addText("Stable Functional Group Scheme for RQ3", {
    x: 0.68,
    y: 2.98,
    w: 12.0,
    h: 0.5,
    align: "center",
    fontFace: "Aptos Display",
    fontSize: 28,
    bold: false,
    color: C.text,
    margin: 0,
  });
  slide.addText("Systematic diagnosis, revised figures, and presentation-ready summary", {
    x: 1.2,
    y: 4.92,
    w: 10.9,
    h: 0.28,
    align: "center",
    fontSize: 14,
    color: C.text,
    margin: 0,
  });
  slide.addText("Chenchen Ding | PKU-IIASA Joint Postdoc", {
    x: 1.55,
    y: 5.55,
    w: 10.1,
    h: 0.26,
    align: "center",
    fontSize: 16,
    color: C.text,
    margin: 0,
  });
  slide.addText("March 2026", {
    x: 4.9,
    y: 6.07,
    w: 3.5,
    h: 0.22,
    align: "center",
    fontSize: 12,
    color: C.muted,
    margin: 0,
  });
  slide.addNotes(
    "Today I will focus on the third research question. I will briefly show how I rebuilt the functional group scheme, what I corrected in the plotting workflow, and what the final interaction patterns look like. The key message is that the revised FG scheme is more stable, more interpretable, and presentation-ready."
  );
  finalizeSlide(slide, "1");
}

function refinementSlide() {
  const slide = pptx.addSlide();
  addTopChrome(
    slide,
    "What Was Revised Before Summarising Results",
    "A short diagnosis-to-fix workflow for the third scientific question."
  );

  const cards = [
    {
      x: 0.55,
      title: "FG construction",
      body: "Replaced unstable raw-similarity grouping with a robust PCA-based clustering workflow using transformed traits, minimum group-size control, and bootstrap Jaccard stability.",
      color: "EAF7F3",
    },
    {
      x: 4.45,
      title: "Prediction figures",
      body: "Removed jagged prediction curves by switching from pointwise random perturbation to joint fixed-effect draws from the model covariance matrix, yielding smoother and more coherent uncertainty bands.",
      color: "EEF4FB",
    },
    {
      x: 8.35,
      title: "Two-way comparison",
      body: "Corrected the fixed-temperature panel mapping so each panel now compares FG1, FG2, and FG3 within the same temperature condition across land-use types.",
      color: "FFF4E8",
    },
  ];

  cards.forEach((card) => {
    slide.addShape(pptx.ShapeType.roundRect, {
      x: card.x,
      y: 2.2,
      w: 3.45,
      h: 2.55,
      rectRadius: 0.08,
      line: { color: "D8D8D8", pt: 1 },
      fill: { color: card.color },
    });
    slide.addText(card.title, {
      x: card.x + 0.2,
      y: 2.38,
      w: 3.0,
      h: 0.22,
      fontSize: 14,
      bold: true,
      color: C.text,
      margin: 0,
    });
    slide.addText(card.body, {
      x: card.x + 0.2,
      y: 2.78,
      w: 3.05,
      h: 1.6,
      fontSize: 10,
      color: C.text,
      valign: "top",
      margin: 0,
    });
  });

  addKeyPanel(slide, {
    x: 0.8,
    y: 4.9,
    w: 11.75,
    h: 1.45,
    title: "Why this matters",
    bullets: [
      "The final slides summarise the corrected and presentation-quality version of the analysis, not the earlier unstable or noisy outputs.",
      "This makes the FG results easier to defend scientifically because the clustering logic, ecological interpretation, and interaction plots now align with each other.",
    ],
  });
  slide.addNotes(
    "Before showing the results, I want to be explicit that these are not just cosmetic updates. I revised the FG construction, corrected the source of jagged prediction bands, and fixed the panel mapping for the land-use-by-FG figure. So what follows is the cleaned and internally consistent version of the analysis."
  );
  // Intentional overlap note: the full-width summary panel is a background box
  // containing text, which the generic overlap checker flags even though this
  // is deliberate and visually correct.
  finalizeSlide(slide, "2", { skipOverlapCheck: true });
}

function workflowSlide() {
  const slide = pptx.addSlide();
  addTopChrome(
    slide,
    "End-to-End Analytical Workflow",
    "From species-level traits to interaction models and interpretation figures."
  );

  const steps = [
    ["1. Species traits", "Median trait values per species: RS, HB, TR, HWI, GL, CS, BM."],
    ["2. Robust trait space", "Log-transform skewed traits, winsorise tails, standardise, then run PCA."],
    ["3. Candidate FG schemes", "Compare kmeans, PAM, and Ward clustering with K = 2 to 5."],
    ["4. Stability filter", "Retain only solutions with acceptable silhouette, minimum group proportion, and bootstrap Jaccard stability."],
    ["5. Final FG interpretation", "Map the selected FG back to raw traits and summarise ecological syndromes."],
    ["6. Interaction modelling", "Fit UI2 x warming x FG occurrence models, then derive two-way and three-way summaries."],
  ];

  let x = 0.6;
  let y = 2.1;
  steps.forEach((s, idx) => {
    slide.addShape(pptx.ShapeType.roundRect, {
      x,
      y,
      w: 3.72,
      h: 0.95,
      rectRadius: 0.08,
      line: { color: "D4DFDC", pt: 1 },
      fill: { color: idx % 2 === 0 ? "F4FBF8" : "F8FAFD" },
    });
    slide.addText(s[0], {
      x: x + 0.16,
      y: y + 0.12,
      w: 3.25,
      h: 0.18,
      fontSize: 12,
      bold: true,
      color: C.text,
      margin: 0,
    });
    slide.addText(s[1], {
      x: x + 0.16,
      y: y + 0.35,
      w: 3.22,
      h: 0.42,
      fontSize: 9.2,
      color: C.text,
      margin: 0,
      valign: "top",
    });

    if (x < 8.0) {
      x += 4.0;
    } else {
      x = 0.6;
      y += 1.18;
    }
  });

  addKeyPanel(slide, {
    x: 9.18,
    y: 4.78,
    w: 3.42,
    h: 1.55,
    title: "Final choice",
    bullets: [
      "Selected scheme: kmeans, K = 3.",
      "Species per group: 588, 828, and 1429.",
      "Main interpretation: slow-history/high-dispersal, broad-niche/high-fecundity, and range-restricted/specialist groups.",
    ],
  });
  slide.addNotes(
    "This slide shows the logic as a pipeline. I started from species-level traits, built a robust multivariate trait space, screened multiple clustering candidates, and kept only schemes that were both interpretable and stable. Then I linked the final FG back to occurrence records and used the same FG definition consistently across all interaction plots."
  );
  finalizeSlide(slide, "3");
}

function figureSlides() {
  const configs = [
    {
      title: "Selecting a Stable and Interpretable FG Scheme",
      subtitle: "Candidate screening across clustering method, group number, balance, and bootstrap stability.",
      image: path.join(FIG_EXPL, "FG_candidate_selection_tradeoff.png"),
      figureLabel: "Figure: candidate selection trade-off for the final FG scheme.",
      bullets: [
        "The selected solution is kmeans with K = 3, because it remains close to the best silhouette while also passing group-size and stability thresholds.",
        "K = 2 is also statistically acceptable, but K = 3 preserves more ecological detail and supports clearer interpretation.",
        "Higher-K solutions lose stability or create overly small groups, so they were not retained.",
      ],
      accent: C.mintDark,
      notes:
        "This figure explains why I chose the final FG definition. I did not pick the grouping only because it looked good visually. Instead, I screened multiple methods and K values, and the three-group kmeans solution offered the best compromise between stability, balance, and interpretability.",
    },
    {
      title: "Functional Strategy Space in PCA Coordinates",
      subtitle: "The selected FG scheme separates species into three distinct parts of multivariate trait space.",
      image: path.join(FIG_EXPL, "Fig3_PCA_biplot_FG_RobustStable.png"),
      figureLabel: "Figure 3: PCA space and final functional strategy groups.",
      bullets: [
        "The three FG occupy distinct regions of the robust PCA space rather than collapsing into a tiny outlier cluster.",
        "FG1 is pulled toward generation length, dispersal ability, and body mass; FG2 toward thermal range, clutch size, and habitat breadth.",
        "FG3 lies on the opposite side for several traits, consistent with a more range-restricted and specialist strategy.",
      ],
      accent: "D8E7FF",
      notes:
        "Here the important point is separation, not just clustering. The selected groups occupy distinct regions of trait space and align with meaningful loading directions. That gives us confidence that the FG scheme is not only statistically stable, but also biologically interpretable.",
    },
    {
      title: "Trait Syndromes Defining the Three FG",
      subtitle: "Species-standardised raw-trait profiles show complementary ecological strategies.",
      image: path.join(FIG_EXPL, "Fig4_FG_trait_profiles_Zscore_bar.png"),
      figureLabel: "Figure 4: mean trait z-scores across species for each FG.",
      bullets: [
        "FG1 is characterised by longer generation length, stronger dispersal capacity, and larger body size.",
        "FG2 stands out for broader thermal range, wider habitat breadth, and higher clutch size.",
        "FG3 is below the overall mean for most traits, consistent with a narrower and more specialist ecological profile.",
      ],
      accent: "FDE7E7",
      notes:
        "This is the clearest slide for naming the FG. FG1 behaves like a slow-history and high-dispersal strategy, FG2 looks like a broad-niche and higher-fecundity strategy, and FG3 is the comparatively restricted specialist group. These labels are grounded in the trait profiles rather than assigned post hoc.",
    },
    {
      title: "Heatmap of FG Trait Signatures",
      subtitle: "The heatmap summarises the same syndromes in a compact pattern-based view.",
      image: path.join(FIG_EXPL, "Fig5_FG_trait_signatures_heatmap_Zscore.png"),
      figureLabel: "Figure 5: FG trait signatures across all raw traits.",
      bullets: [
        "The heatmap makes the complementarity among the three strategies very clear.",
        "FG1 is dominated by slow-history traits, FG2 by broad niche and fecundity traits, and FG3 by comparatively low scores across several dimensions.",
        "This gives a concise visual justification for keeping three groups rather than collapsing them into two.",
      ],
      accent: "FBEFD8",
      notes:
        "I use this slide as a compact summary of the trait syndromes. Compared with the bar chart, the heatmap is better for quickly seeing contrast among groups across all traits at once. It reinforces that the three-group solution is both structured and ecologically coherent.",
    },
    {
      title: "Three-Way Interaction on the Absolute Probability Scale",
      subtitle: "Predicted occurrence probability varies with warming, land-use type, and FG simultaneously.",
      image: path.join(FIG_3WAY, "Fig1_ThreeWay_ABS_UI2_Temp_FG_glmmTMB_fastCI.png"),
      figureLabel: "Figure 1: three-way interaction on the absolute occurrence scale.",
      bullets: [
        "The warming-response curve is not the same across land-use types, and it also differs among FG.",
        "FG1 and FG3 show especially strong positive warming responses in urban contexts, while several other combinations are flatter or even slightly negative.",
        "This supports a heterogeneous three-way interaction between land use, warming, and functional strategy.",
      ],
      accent: "E8F5EE",
      notes:
        "This is the central biological result. The response to warming depends both on land-use context and on the functional strategy of the species group. So the third scientific question is answered positively: trait-based strategies do modulate the land-use-by-warming response.",
    },
    {
      title: "Three-Way Interaction Relative to the Primary-Vegetation Baseline",
      subtitle: "Relative change highlights how warming sensitivity depends on both ecological strategy and disturbance context.",
      image: path.join(FIG_3WAY, "Fig2_ThreeWay_PCT_vsPV0_UI2_Temp_FG_glmmTMB_fastCI.png"),
      figureLabel: "Figure 2: three-way interaction as percent change relative to the baseline.",
      bullets: [
        "Relative change from the baseline differs strongly among both FG and land-use contexts.",
        "This view makes it easier to compare sensitivity, because all panels are anchored to the same reference condition.",
        "The main message remains the same: warming sensitivity is conditional on both functional strategy and land-use background.",
      ],
      accent: "F7EAFE",
      notes:
        "I usually show this figure right after the absolute scale because it helps compare sensitivity more directly. The baseline is fixed, so differences in curvature and effect size become easier to communicate. It is especially useful when discussing which strategy groups are most responsive under intensified land use.",
    },
    {
      title: "Two-Way Summary: Land Use by Warming",
      subtitle: "Averaging over FG still leaves clear land-use differences in baseline occurrence and warming-response shape.",
      image: path.join(FIG_2WAY, "Fig6_TwoWay_UI2_by_Warming_FGaveraged_ABS_glmmTMB_fastCI.png"),
      figureLabel: "Figure 6: land use by warming after averaging across FG composition.",
      bullets: [
        "Even after averaging over FG, land-use types still differ in both baseline occurrence and response shape.",
        "This means land use has an overall filtering effect that is not explained away by FG composition alone.",
        "Agricultural and urban settings remain especially informative when comparing warming trajectories.",
      ],
      accent: "EAF7F3",
      notes:
        "This is a marginal summary derived from the same three-way model. Even when I average over the FG structure, land use still matters strongly. So the land-use signal is not only a byproduct of group composition; it persists as a broader ecological filter.",
    },
    {
      title: "Two-Way Summary: FG by Warming",
      subtitle: "Averaging over land-use composition still preserves clear strategy-specific warming trajectories.",
      image: path.join(FIG_2WAY, "Fig7_TwoWay_FG_by_Warming_UI2averaged_ABS_glmmTMB_fastCI.png"),
      figureLabel: "Figure 7: FG by warming after averaging across land-use composition.",
      bullets: [
        "FG retain different warming-response trajectories after averaging across land-use types.",
        "This shows that the FG scheme captures real biological heterogeneity in climate sensitivity.",
        "FG3 and FG1 are more temperature-responsive overall than FG2 in the averaged view.",
      ],
      accent: "EAF0FF",
      notes:
        "This panel complements the previous one. Here I average over land-use composition and ask whether the FG definition still matters. The answer is yes, which strengthens the argument that the FG scheme is not only statistically tidy but biologically meaningful.",
    },
    {
      title: "Two-Way Summary at Fixed Warming Levels",
      subtitle: "Within each temperature panel, FG1, FG2, and FG3 can now be compared directly across land-use types.",
      image: path.join(FIG_2WAY, "Fig8_TwoWay_UI2_by_FG_at_fixedWarming_ABS_glmmTMB_fastCI.png"),
      figureLabel: "Figure 8: land use by FG at fixed warming levels.",
      bullets: [
        "Each panel now shows all three FG under the same fixed temperature condition, which is the correct comparison for this figure.",
        "The rank order among FG changes with land-use type, especially under stronger warming.",
        "This indicates that land-use filtering and functional strategy interact directly in shaping occurrence probability.",
      ],
      accent: "FFF0E6",
      notes:
        "This figure was important to fix. The corrected version now compares FG within the same temperature panel, which is what we actually need for interpretation. It shows more clearly that the relative ordering of the groups depends on land use and shifts as warming increases.",
    },
  ];

  configs.forEach((cfg, idx) => {
    const slide = pptx.addSlide();
    addFigureSlide(slide, cfg);
    finalizeSlide(slide, String(idx + 4));
  });
}

function synthesisSlide() {
  const slide = pptx.addSlide();
  addTopChrome(
    slide,
    "Take-Home Messages",
    "What this revised FG analysis contributes to the third scientific question."
  );

  const messages = [
    {
      title: "Methodological conclusion",
      body: "A robust PCA-plus-stability workflow provides a more defensible FG scheme than earlier unstable similarity-based groupings.",
      color: "EAF7F3",
    },
    {
      title: "Ecological conclusion",
      body: "The three FG correspond to interpretable trait syndromes: slow-history/high-dispersal, broad-niche/high-fecundity, and range-restricted/specialist strategies.",
      color: "EEF4FB",
    },
    {
      title: "Interaction conclusion",
      body: "Occurrence responses to warming are conditional on both land-use context and functional strategy, supporting the core RQ3 hypothesis.",
      color: "FFF4E8",
    },
  ];

  messages.forEach((m, i) => {
    const x = 0.72 + i * 4.15;
    slide.addShape(pptx.ShapeType.roundRect, {
      x,
      y: 2.18,
      w: 3.65,
      h: 2.45,
      rectRadius: 0.08,
      line: { color: "D7DFDB", pt: 1 },
      fill: { color: m.color },
    });
    slide.addText(m.title, {
      x: x + 0.2,
      y: 2.4,
      w: 3.15,
      h: 0.2,
      fontSize: 13,
      bold: true,
      color: C.text,
      margin: 0,
    });
    slide.addText(m.body, {
      x: x + 0.2,
      y: 2.83,
      w: 3.1,
      h: 1.35,
      fontSize: 10,
      color: C.text,
      margin: 0,
      valign: "top",
    });
  });

  addKeyPanel(slide, {
    x: 1.1,
    y: 4.85,
    w: 11.05,
    h: 1.5,
    title: "Suggested closing line",
    bullets: [
      "The revised FG framework is scientifically cleaner, visually clearer, and much easier to explain in a paper or presentation.",
      "A Bayesian brms version can now be used as a next-step robustness check on the same stable FG definition, rather than as a rescue for an unstable grouping scheme.",
    ],
  });
  slide.addNotes(
    "To close, I would emphasise three things. First, the FG definition is now stable and interpretable. Second, the ecological labels are supported by the trait data. Third, the interaction results consistently show that warming responses depend on both land use and functional strategy. That is the main contribution of this revised analysis."
  );
  // Intentional overlap note: the closing summary panel is a background box
  // containing text, so we skip the generic overlap warning here.
  finalizeSlide(slide, "13", { skipOverlapCheck: true });
}

coverSlide();
refinementSlide();
workflowSlide();
figureSlides();
synthesisSlide();

const outPptx = path.join(OUT_DIR, "RQ3_FG_stable_scheme_summary.pptx");

(async () => {
  await pptx.writeFile({ fileName: outPptx });
  console.log(`Wrote ${outPptx}`);
})();
