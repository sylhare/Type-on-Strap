jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import fs from 'node:fs';
import { sha256File, sha256Buffer } from '../../src/utils/hash';
import { fetchBuffer } from '../../src/utils/http';
import { validate } from '../../src/validators/font-awesome';

const mockFs = fs as jest.Mocked<typeof fs>;
const mockSha256File = sha256File as jest.MockedFunction<typeof sha256File>;
const mockSha256Buffer = sha256Buffer as jest.MockedFunction<typeof sha256Buffer>;
const mockFetchBuffer = fetchBuffer as jest.MockedFunction<typeof fetchBuffer>;

const MOCK_CONFIG = JSON.stringify({ fontAwesome: { version: '6.7.2' } });
const MOCK_HASH = 'a'.repeat(64);
const realReadFileSync = (jest.requireActual('node:fs') as typeof import('node:fs')).readFileSync;

describe('validators/font-awesome', () => {
  beforeEach(() => {
    mockFs.readFileSync = jest.fn().mockImplementation((filePath: unknown, options?: any) => {
      if (String(filePath).endsWith('vendor.config.json')) return MOCK_CONFIG;
      return realReadFileSync(filePath as any, options);
    }) as jest.MockedFunction<typeof fs.readFileSync>;
    mockSha256File.mockResolvedValue(MOCK_HASH);
    mockSha256Buffer.mockReturnValue(MOCK_HASH);
    mockFetchBuffer.mockResolvedValue(Buffer.from('remote content'));
  });

  test('passes when all 25 files match (6 fonts + 19 SCSS)', async () => {
    const result = await validate();
    expect(result.passed).toEqual(true);
    expect(result.failures).toHaveLength(0);
    expect(mockFetchBuffer).toHaveBeenCalledTimes(25);
  });

  test('fails when any file hash does not match', async () => {
    // First call fails
    mockSha256Buffer.mockReturnValueOnce('b'.repeat(64));
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures.length).toBeGreaterThan(0);
  });

  test('fails when a local file does not exist', async () => {
    mockSha256File.mockRejectedValue(new Error("ENOENT: no such file or directory, open 'font-file'"));
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures.length).toEqual(25);
  });

  test('uses cdnjs for font files', async () => {
    await validate();
    expect(mockFetchBuffer).toHaveBeenCalledWith(
      expect.stringContaining('cdnjs.cloudflare.com')
    );
  });

  test('uses GitHub raw for SCSS files', async () => {
    await validate();
    expect(mockFetchBuffer).toHaveBeenCalledWith(
      expect.stringContaining('raw.githubusercontent.com')
    );
  });
});
