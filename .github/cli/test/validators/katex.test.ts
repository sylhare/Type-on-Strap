jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import fs from 'node:fs';
import { validate } from '../../src/validators/katex';
import { mockReadFileSync, mockSha256File, mockSha256Buffer, mockFetchBuffer, MOCK_HASH } from '../helpers';

const mockFs = fs as jest.Mocked<typeof fs>;

const MOCK_CONFIG = JSON.stringify({ katex: { version: '0.16.38' } });
const MOCK_SCSS = 'some content\ncontent: "0.16.38"\nmore content';
const MOCK_CDN_CSS = Buffer.from('content: "0.16.38"');

describe('validators/katex', () => {
  beforeEach(() => {
    mockFs.existsSync = jest.fn().mockReturnValue(true);
    mockReadFileSync(MOCK_CONFIG, filePath => {
      if (filePath.endsWith('.scss')) return MOCK_SCSS;
    });
    mockSha256File.mockResolvedValue(MOCK_HASH);
    mockSha256Buffer.mockReturnValue(MOCK_HASH);
    mockFetchBuffer.mockResolvedValue(MOCK_CDN_CSS);
  });

  test('passes when all checks succeed', async () => {
    const result = await validate();
    expect(result.passed).toEqual(true);
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
    mockReadFileSync(MOCK_CONFIG, filePath => {
      if (filePath.endsWith('.scss')) return 'content: "0.16.00"';
    });

    const result = await validate();
    expect(result.failures).toContain('katex.scss version');
  });

  test('fails when local JS file does not exist', async () => {
    mockFs.existsSync = jest.fn().mockReturnValue(false);
    const result = await validate();
    expect(result.passed).toEqual(false);
  });
});
