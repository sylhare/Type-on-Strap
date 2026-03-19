const isColorEnabled = Boolean(process.stdout.isTTY) && !process.env['NO_COLOR'];
const isSilent = process.env['NODE_ENV'] === 'test';

const green  = (s: string) => isColorEnabled ? `\x1b[0;32m${s}\x1b[0m` : s;
const red    = (s: string) => isColorEnabled ? `\x1b[0;31m${s}\x1b[0m` : s;
const yellow = (s: string) => isColorEnabled ? `\x1b[1;33m${s}\x1b[0m` : s;
const blue   = (s: string) => isColorEnabled ? `\x1b[0;34m${s}\x1b[0m` : s;

const noop = (_msg: string): void => {};

export const logger = {
  success: isSilent ? noop : (msg: string) => console.log(green(`✅ ${msg}`)),
  error:   isSilent ? noop : (msg: string) => console.error(red(`❌ ${msg}`)),
  warn:    isSilent ? noop : (msg: string) => console.warn(yellow(`⚠️  ${msg}`)),
  info:    isSilent ? noop : (msg: string) => console.log(msg),
  header:  isSilent ? noop : (msg: string) => console.log(blue(`\n${'='.repeat(50)}\n${msg}\n${'='.repeat(50)}\n`)),
};
