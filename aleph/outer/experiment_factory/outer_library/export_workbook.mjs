import fs from "node:fs/promises";
import path from "node:path";
import { createHash } from "node:crypto";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: export_workbook.mjs INPUT.xlsx OUTPUT.json");
}

const input = await FileBlob.load(inputPath);
const sourceSha256 = createHash("sha256")
  .update(await fs.readFile(inputPath))
  .digest("hex");
const workbook = await SpreadsheetFile.importXlsx(input);
const overview = await workbook.inspect({
  kind: "workbook,sheet",
  include: "id,name,index,range",
  maxChars: 100_000,
});
const names = [];
for (const line of overview.ndjson.split("\n")) {
  if (!line.trim()) continue;
  const item = JSON.parse(line);
  if (item.kind === "sheet" && item.name && !names.includes(item.name)) {
    names.push(item.name);
  }
}

const sheets = [];
for (const name of names) {
  const sheet = workbook.worksheets.getItem(name);
  const used = sheet.getUsedRange(true);
  sheets.push({
    name,
    address: used.address?.split("!").at(-1) ?? "A1",
    values: used.values,
  });
}
await fs.mkdir(path.dirname(outputPath), { recursive: true });
await fs.writeFile(
  outputPath,
  `${JSON.stringify({ input: inputPath, source_sha256: sourceSha256, sheets })}\n`,
);
console.log(JSON.stringify({ output: outputPath, sheet_count: sheets.length }));
