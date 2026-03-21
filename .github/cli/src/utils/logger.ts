const isColorEnabled = Boolean(process.stdout.isTTY) && !process.env['NO_COLOR'];
const isSilent = process.env['NODE_ENV'] === 'test';

const green = (s: string) => isColorEnabled ? `\x1b[0;32m${s}\x1b[0m` : s;
const red = (s: string) => isColorEnabled ? `\x1b[0;31m${s}\x1b[0m` : s;
const yellow = (s: string) => isColorEnabled ? `\x1b[1;33m${s}\x1b[0m` : s;
const blue = (s: string) => isColorEnabled ? `\x1b[0;34m${s}\x1b[0m` : s;

const noop = (_msg: string): void => {
};
const wrap = (fn: (msg: string) => void): (msg: string) => void => isSilent ? noop : fn;

export const logger = {
  success: wrap((msg) => console.log(green(`✅ ${msg}`))),
  error: wrap((msg) => console.error(red(`❌ ${msg}`))),
  warn: wrap((msg) => console.warn(yellow(`⚠️  ${msg}`))),
  info: wrap((msg) => console.log(msg)),
  header: wrap((msg) => console.log(blue(`\n${'='.repeat(50)}\n${msg}\n${'='.repeat(50)}\n`))),
};
