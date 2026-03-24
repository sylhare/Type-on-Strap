jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import { validate } from '../../src/validators/jekyll-search';
import { mockFetchBuffer, mockHashAndFetch, mockReadFileSync, mockSha256Buffer, mockSha256File } from '../helpers';

const MOCK_CONFIG = JSON.stringify({ simpleJekyllSearch: { version: '2.1.2' } });

describe('validators/jekyll-search', () => {
  beforeEach(() => {
    mockReadFileSync(MOCK_CONFIG);
    mockHashAndFetch();
  });

  test('passes when local and remote hashes match', async () => {
    const result = await validate();
    expect(result.passed).toEqual(true);
    expect(result.failures).toHaveLength(0);
  });

  test('fails when hashes do not match', async () => {
    mockSha256Buffer.mockReturnValue('b'.repeat(64));
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures).toContain('simple-jekyll-search.min.js');
  });

  test('fails when local file does not exist', async () => {
    mockSha256File.mockRejectedValueOnce(new Error("ENOENT: no such file or directory, open 'simple-jekyll-search.min.js'"));
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures).toContain('simple-jekyll-search.min.js');
  });

  test('requests correct GitHub raw URL', async () => {
    await validate();
    expect(mockFetchBuffer).toHaveBeenCalledWith(
      expect.stringContaining('Simple-Jekyll-Search/v2.1.2/dest/simple-jekyll-search.min.js')
    );
  });
});
