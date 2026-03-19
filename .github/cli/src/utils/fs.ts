import fs from 'node:fs';

export function updateVersionInFile(filePath: string, pattern: RegExp, replacement: string): void {
  const content = fs.readFileSync(filePath, 'utf8');
  if (!pattern.test(content)) throw new Error(`Pattern not found in ${filePath}`);
  const updated = content.replace(pattern, replacement);
  if (updated !== content) fs.writeFileSync(filePath, updated);
}

export function updateVendorConfig(configPath: string, key: string, version: string): void {
  const config = JSON.parse(fs.readFileSync(configPath, 'utf8')) as Record<string, { version: string }>;
  config[key] = { version };
  fs.writeFileSync(configPath, JSON.stringify(config, null, 2) + '\n');
}

export function readVendorVersion(configPath: string, key: string): string {
  const config = JSON.parse(fs.readFileSync(configPath, 'utf8')) as Record<string, { version: string }>;
  return config[key].version;
}

export function updateJsonFile<T extends Record<string, unknown>>(
  filePath: string,
  updater: (data: T) => void
): void {
  const data = JSON.parse(fs.readFileSync(filePath, 'utf8')) as T;
  updater(data);
  fs.writeFileSync(filePath, JSON.stringify(data, null, 2) + '\n');
}
