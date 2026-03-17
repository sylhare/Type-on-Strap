jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import fs from 'node:fs';
import { sha256File, sha256Buffer } from '../../src/utils/hash';
import { fetchBuffer } from '../../src/utils/http';
import { validate } from '../../src/validators/katex';

const mockFs = fs as jest.Mocked<typeof fs>;
const mockSha256File = sha256File as jest.MockedFunction<typeof sha256File>;
const mockSha256Buffer = sha256Buffer as jest.MockedFunction<typeof sha256Buffer>;
const mockFetchBuffer = fetchBuffer as jest.MockedFunction<typeof fetchBuffer>;

const MOCK_CONFIG = JSON.stringify({ katex: { version: '0.16.38' } });
const MOCK_HASH = 'a'.repeat(64);
const MOCK_SCSS = 'some content\ncontent: "0.16.38"\nmore content';
const MOCK_CDN_CSS = Buffer.from('content: "0.16.38"');
const realReadFileSync = (jest.requireActual('node:fs') as typeof import('node:fs')).readFileSync;

describe('validators/katex', () => {
  beforeEach(() => {
    mockFs.existsSync = jest.fn().mockReturnValue(true);
    mockFs.readFileSync = jest.fn().mockImplementation((filePath: unknown, options?: any) => {
      if (String(filePath).endsWith('vendor.config.json')) return MOCK_CONFIG;
      if (String(filePath).endsWith('.scss')) return MOCK_SCSS;
      return realReadFileSync(filePath as any, options);
    }) as jest.MockedFunction<typeof fs.readFileSync>;
    mockSha256File.mockResolvedValue(MOCK_HASH);
    mockSha256Buffer.mockReturnValue(MOCK_HASH);
    mockFetchBuffer.mockResolvedValue(MOCK_CDN_CSS);
  });

  test('passes when all checks succeed', async () => {
    const result = await validate();
    expect(result.passed).toBe(true);
    expect(result.failures).toHaveLength(0);
  });

  test('fails katex.min.js when hashes do not match', async () => {
    mockSha256Buffer
      .mockReturnValueOnce('b'.repeat(64))
      .mockReturnValueOnce(MOCK_HASH)
      .mockReturnValue(MOCK_HASH);

    const result = await validate();
    expect(result.failures).toContain('katex.min.js');
  });

  test('fails scss version check when local version does not match', async () => {
    mockFs.readFileSync = jest.fn().mockImplementation((filePath: unknown, options?: any) => {
      if (String(filePath).endsWith('vendor.config.json')) return MOCK_CONFIG;
      if (String(filePath).endsWith('.scss')) return 'content: "0.16.00"';
      return realReadFileSync(filePath as any, options);
    }) as jest.MockedFunction<typeof fs.readFileSync>;

    const result = await validate();
    expect(result.failures).toContain('katex.scss version');
  });

  test('fails when local JS file does not exist', async () => {
    mockFs.existsSync = jest.fn().mockReturnValue(false);
    const result = await validate();
    expect(result.passed).toBe(false);
  });
});
