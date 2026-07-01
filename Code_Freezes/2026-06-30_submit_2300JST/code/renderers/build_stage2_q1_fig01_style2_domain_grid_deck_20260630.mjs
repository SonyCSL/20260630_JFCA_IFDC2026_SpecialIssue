#!/usr/bin/env node
import fs from "node:fs/promises";
import path from "node:path";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const DOMAINS = ["Formula", "Brite", "Pathway"];
const DOMAIN_LABELS = new Map([["Brite", "BRITE"]]);
const SUBSETS = ["All", "PM", "SM"];
const DEFAULT_MODES = [
  "Combo6yr",
  "Period1_2015_2017",
  "Period2_2018_2020",
  "2015",
  "2016",
  "2017",
  "2018",
  "2019",
  "2020",
];
const SLIDE = { width: 1920, height: 1080 };

function parseArgs(argv) {
  const args = { modes: DEFAULT_MODES.join(","), qLabel: "q1", transform: "raw" };
  for (let i = 2; i < argv.length; i += 1) {
    const item = argv[i];
    if (!item.startsWith("--")) {
      throw new Error(`Unexpected argument: ${item}`);
    }
    const key = item.slice(2).replace(/-([a-z])/g, (_, ch) => ch.toUpperCase());
    const value = argv[i + 1];
    if (value === undefined || value.startsWith("--")) {
      throw new Error(`Missing value for ${item}`);
    }
    args[key] = value;
    i += 1;
  }
  for (const key of ["figRoot", "out"]) {
    if (!args[key]) {
      throw new Error(`Missing required --${key.replace(/[A-Z]/g, (ch) => `-${ch.toLowerCase()}`)}`);
    }
  }
  return args;
}

function splitCsv(text) {
  return text.split(",").map((item) => item.trim()).filter(Boolean);
}

async function readImage(pathname) {
  const bytes = await fs.readFile(pathname);
  return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
}

async function exists(pathname) {
  try {
    await fs.access(pathname);
    return true;
  } catch {
    return false;
  }
}

async function writeBlob(pathname, blob) {
  await fs.mkdir(path.dirname(pathname), { recursive: true });
  await fs.writeFile(pathname, new Uint8Array(await blob.arrayBuffer()));
}

function addText(slide, text, position, style = {}) {
  const shape = slide.shapes.add({
    geometry: "textbox",
    position,
    fill: "none",
    line: { style: "solid", fill: "none", width: 0 },
  });
  shape.text = text;
  shape.text.style = {
    fontSize: style.fontSize ?? 24,
    bold: style.bold ?? false,
    color: style.color ?? "#111827",
    alignment: style.alignment ?? "left",
  };
  return shape;
}

function addRule(slide, left, top, width) {
  slide.shapes.add({
    geometry: "rect",
    position: { left, top, width, height: 1 },
    fill: "#E5E7EB",
    line: { style: "solid", fill: "#E5E7EB", width: 0 },
  });
}

async function addImage(slide, pathname, alt, position, fit = "contain") {
  if (!(await exists(pathname))) {
    throw new Error(`Missing image: ${pathname}`);
  }
  slide.images.add({
    blob: await readImage(pathname),
    contentType: "image/png",
    alt,
    fit,
    position,
  });
}

function figurePath(figRoot, transform, domain, subset, mode, qLabel) {
  return path.join(
    figRoot,
    transform,
    "figures",
    "main_q12",
    domain,
    subset,
    mode,
    qLabel,
    "Stage2_iNEXT_TD_m_est",
    "CliffsDelta",
    "01_reference_bridge_boxplots_paired_nolegend.png",
  );
}

function modeTitle(mode) {
  if (mode === "Combo6yr") return "Combo6yr";
  if (mode === "Period1_2015_2017") return "Period1 2015-2017";
  if (mode === "Period2_2018_2020") return "Period2 2018-2020";
  return mode;
}

async function addModeSlide(presentation, args, mode, slideIndex, slideCount) {
  const slide = presentation.slides.add();
  slide.background.fill = "#FFFFFF";

  addText(slide, `${args.transform} | ${modeTitle(mode)}`, { left: 58, top: 42, width: 420, height: 42 }, { fontSize: 30, bold: true });
  addText(slide, "Formula / BRITE / Pathway", { left: 745, top: 52, width: 430, height: 28 }, { fontSize: 18, bold: true, alignment: "center" });
  addText(slide, `${args.qLabel} | Stage2 iNEXT TD.m.est | Cliff's delta`, { left: 1380, top: 54, width: 420, height: 24 }, { fontSize: 13, color: "#374151", alignment: "right" });

  const colLefts = [430, 820, 1210];
  const panel = { width: 370, height: 188 };
  for (const [idx, subset] of SUBSETS.entries()) {
    addText(slide, subset.toUpperCase(), { left: colLefts[idx], top: 118, width: panel.width, height: 24 }, { fontSize: 18, bold: true, alignment: "center" });
  }

  const rowTops = [160, 420, 680];
  for (const [rowIdx, domain] of DOMAINS.entries()) {
    const label = DOMAIN_LABELS.get(domain) ?? domain;
    addText(slide, label, { left: 110, top: rowTops[rowIdx] + 68, width: 230, height: 34 }, { fontSize: 18, bold: true, alignment: "center" });
    for (const [colIdx, subset] of SUBSETS.entries()) {
      const imagePath = figurePath(args.figRoot, args.transform, domain, subset, mode, args.qLabel);
      await addImage(
        slide,
        imagePath,
        `${label} ${subset} ${modeTitle(mode)}`,
        { left: colLefts[colIdx], top: rowTops[rowIdx], width: panel.width, height: panel.height },
      );
    }
    if (rowIdx < DOMAINS.length - 1) {
      addRule(slide, 195, rowTops[rowIdx] + 225, 1535);
    }
  }

  await addImage(
    slide,
    path.join(args.figRoot, "figure01_boxplot_legend_style2_horizontal_single_row.png"),
    "Figure 01 style2 legend",
    { left: 175, top: 925, width: 1540, height: 56 },
  );

  addText(
    slide,
    "Descriptive Stage2 q1 diagnostic. Boxes summarize Monte Carlo/Stage2 bootstrap distributions; cells are not independent hypothesis tests.",
    { left: 66, top: 1004, width: 1350, height: 24 },
    { fontSize: 12, color: "#6B7280" },
  );
  addText(slide, `${slideIndex}/${slideCount}`, { left: 1775, top: 1004, width: 90, height: 24 }, { fontSize: 12, color: "#6B7280", alignment: "right" });
}

async function main() {
  const args = parseArgs(process.argv);
  args.figRoot = path.resolve(args.figRoot);
  args.out = path.resolve(args.out);
  args.qaDir = path.resolve(args.qaDir || `${args.out}.qa`);
  const modes = splitCsv(args.modes);

  const presentation = Presentation.create({ slideSize: SLIDE });
  for (const [idx, mode] of modes.entries()) {
    await addModeSlide(presentation, args, mode, idx + 1, modes.length);
  }

  await fs.mkdir(path.dirname(args.out), { recursive: true });
  await fs.mkdir(args.qaDir, { recursive: true });
  const previewDir = path.join(args.qaDir, "preview");
  const layoutDir = path.join(args.qaDir, "layout");
  await fs.mkdir(previewDir, { recursive: true });
  await fs.mkdir(layoutDir, { recursive: true });

  for (const [index, slide] of presentation.slides.items.entries()) {
    const stem = `slide-${String(index + 1).padStart(2, "0")}`;
    await writeBlob(path.join(previewDir, `${stem}.png`), await presentation.export({ slide, format: "png", scale: 1 }));
    await fs.writeFile(path.join(layoutDir, `${stem}.layout.json`), await (await slide.export({ format: "layout" })).text());
  }
  await writeBlob(path.join(args.qaDir, "deck_montage.webp"), await presentation.export({ format: "webp", montage: true, scale: 1 }));
  const inspect = await presentation.inspect({ kind: "deck,slide,textbox,shape,image,notes,layout", maxChars: 300000 });
  await fs.writeFile(`${args.out}.inspect.ndjson`, inspect.ndjson);
  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(args.out);
  console.log(`[OK] wrote ${args.out}`);
  console.log(`[OK] wrote QA under ${args.qaDir}`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
