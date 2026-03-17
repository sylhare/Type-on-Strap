jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import fs from 'node:fs';
import { sha256File, sha256Buffer } from '../../src/utils/hash';
import { fetchBuffer } from '../../src/utils/http';
import { validate } from '../../src/validators/mermaid';

const mockFs = fs as jest.Mocked<typeof fs>;
const mockSha256File = sha256File as jest.MockedFunction<typeof sha256File>;
const mockSha256Buffer = sha256Buffer as jest.MockedFunction<typeof sha256Buffer>;
const mockFetchBuffer = fetchBuffer as jest.MockedFunction<typeof fetchBuffer>;

const MOCK_CONFIG = JSON.stringify({ mermaid: { version: '11.13.0' } });
const MOCK_HASH = 'a'.repeat(64);
const realReadFileSync = (jest.requireActual('node:fs') as typeof import('node:fs')).readFileSync;

describe('validators/mermaid', () => {
  beforeEach(() => {
    mockFs.existsSync = jest.fn().mockReturnValue(true);
    mockFs.readFileSync = jest.fn().mockImplementation((filePath: unknown, options?: any) => {
      if (String(filePath).endsWith('vendor.config.json')) return MOCK_CONFIG;
      return realReadFileSync(filePath as any, options);
    }) as jest.MockedFunction<typeof fs.readFileSync>;
    mockSha256File.mockResolvedValue(MOCK_HASH);
    mockSha256Buffer.mockReturnValue(MOCK_HASH);
    mockFetchBuffer.mockResolvedValue(Buffer.from('remote content'));
  });

  test('passes when local and remote hashes match', async () => {
    const result = await validate();
    expect(result.passed).toBe(true);
    expect(result.failures).toHaveLength(0);
  });

  test('fails when hashes do not match', async () => {
    mockSha256Buffer.mockReturnValue('b'.repeat(64));
    const result = await validate();
    expect(result.passed).toBe(false);
    expect(result.failures).toContain('mermaid.min.js');
  });

  test('fails when local file does not exist', async () => {
    mockFs.existsSync = jest.fn().mockReturnValue(false);
    const result = await validate();
    expect(result.passed).toBe(false);
    expect(result.failures).toContain('mermaid.min.js');
  });

  test('requests correct CDN URL', async () => {
    await validate();
    expect(mockFetchBuffer).toHaveBeenCalledWith(
      expect.stringContaining('mermaid@11.13.0')
    );
  });
});
