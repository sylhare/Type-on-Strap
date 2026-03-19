jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import { validate } from '../../src/validators/masonry';
import { mockReadFileSync, mockSha256File, mockSha256Buffer, mockFetchBuffer, mockHashAndFetch } from '../helpers';

const MOCK_CONFIG = JSON.stringify({ masonry: { version: '4.2.2' } });

describe('validators/masonry', () => {
  beforeEach(() => {
    mockReadFileSync(MOCK_CONFIG);
    mockHashAndFetch();
  });

  test('passes when hashes match', async () => {
    const result = await validate();
    expect(result.passed).toEqual(true);
    expect(result.failures).toHaveLength(0);
  });

  test('fails when hashes do not match', async () => {
    mockSha256Buffer.mockReturnValue('b'.repeat(64));
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures).toContain('masonry.pkgd.min.js');
  });

  test('fails when local file does not exist', async () => {
    mockSha256File.mockRejectedValueOnce(new Error("ENOENT: no such file or directory, open 'masonry.pkgd.min.js'"));
    const result = await validate();
    expect(result.passed).toEqual(false);
  });

  test('requests correct CDN URL', async () => {
    await validate();
    expect(mockFetchBuffer).toHaveBeenCalledWith(
      expect.stringContaining('masonry-layout@4.2.2')
    );
  });
});
