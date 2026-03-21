jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import fs from 'node:fs';
import { fetchJson } from '../../src/utils/http';
import { validate } from '../../src/validators/jekyll-search';
import { mockHashAndFetch, mockSha256Buffer, realReadFileSync } from '../helpers';

const mockFs = fs as jest.Mocked<typeof fs>;
const mockFetchJson = fetchJson as jest.MockedFunction<typeof fetchJson>;

const LOCAL_CONTENT = '// Simple-Jekyll-Search v1.10.0\n(function(){var module={};})();\n';

describe('validators/jekyll-search', () => {
  beforeEach(() => {
    mockFs.existsSync = jest.fn().mockReturnValue(true);
    mockFs.readFileSync = jest.fn().mockImplementation((filePath: unknown, options?: any) => {
      if (String(filePath).endsWith('.js')) return LOCAL_CONTENT;
      return realReadFileSync(filePath as any, options);
    }) as jest.MockedFunction<typeof fs.readFileSync>;
    mockHashAndFetch();
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
