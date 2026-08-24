import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const outputPath = process.argv[2];
const inputPaths = process.argv.slice(3);
if (!outputPath || inputPaths.length === 0) {
  throw new Error("usage: export_workbooks.mjs OUTPUT.json INPUT.xlsx [...]");
}

const exports = [];
for (const inputPath of inputPaths) {
  const input = await FileBlob.load(inputPath);
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
    if (item.kind === "sheet" && item.name && !names.includes(item.name)) names.push(item.name);
  }
  const sheets = names.map((name) => {
    const used = workbook.worksheets.getItem(name).getUsedRange(true);
    return {
      name,
      address: used.address?.split("!").at(-1) ?? "A1",
      values: used.values,
    };
  });
  exports.push({ input: inputPath, sheets });
}
await fs.mkdir(path.dirname(outputPath), { recursive: true });
await fs.writeFile(outputPath, `${JSON.stringify({ workbooks: exports })}\n`);
console.log(JSON.stringify({ output: outputPath, workbook_count: exports.length }));
