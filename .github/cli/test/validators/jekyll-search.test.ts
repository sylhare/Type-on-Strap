jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import fs from 'node:fs';
import { sha256File, sha256Buffer } from '../../src/utils/hash';
import { fetchBuffer, fetchJson } from '../../src/utils/http';
import { validate } from '../../src/validators/jekyll-search';

const mockFs = fs as jest.Mocked<typeof fs>;
const mockSha256File = sha256File as jest.MockedFunction<typeof sha256File>;
const mockSha256Buffer = sha256Buffer as jest.MockedFunction<typeof sha256Buffer>;
const mockFetchBuffer = fetchBuffer as jest.MockedFunction<typeof fetchBuffer>;
const mockFetchJson = fetchJson as jest.MockedFunction<typeof fetchJson>;

const LOCAL_CONTENT = '// Simple-Jekyll-Search v1.10.0\n(function(){var module={};})();\n';
const MOCK_HASH = 'a'.repeat(64);
const realReadFileSync = (jest.requireActual('node:fs') as typeof import('node:fs')).readFileSync;

describe('validators/jekyll-search', () => {
  beforeEach(() => {
    mockFs.existsSync = jest.fn().mockReturnValue(true);
    mockFs.readFileSync = jest.fn().mockImplementation((filePath: unknown, options?: any) => {
      if (String(filePath).endsWith('.js')) return LOCAL_CONTENT;
      return realReadFileSync(filePath as any, options);
    }) as jest.MockedFunction<typeof fs.readFileSync>;
    mockSha256File.mockResolvedValue(MOCK_HASH);
    mockSha256Buffer.mockReturnValue(MOCK_HASH);
    mockFetchBuffer.mockResolvedValue(Buffer.from('remote content'));
    mockFetchJson.mockResolvedValue({ tag_name: 'v1.10.0' });
  });

  test('passes when version and content match', async () => {
    const result = await validate();
    expect(result.passed).toEqual(true);
    expect(result.failures).toHaveLength(0);
  });

  test('fails when version does not match latest', async () => {
    mockFetchJson.mockResolvedValue({ tag_name: 'v1.11.0' });
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures).toContain('version mismatch');
  });

  test('fails when content hash does not match', async () => {
    mockSha256Buffer.mockReturnValue('b'.repeat(64));
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures).toContain('content mismatch');
  });

  test('passes gracefully when GitHub API is unreachable', async () => {
    mockFetchJson.mockRejectedValue(new Error('Network error'));
    const result = await validate();
    expect(result.failures).toHaveLength(0);
  });

  test('fails when local file does not exist', async () => {
    mockFs.existsSync = jest.fn().mockReturnValue(false);
    const result = await validate();
    expect(result.passed).toEqual(false);
  });
});
