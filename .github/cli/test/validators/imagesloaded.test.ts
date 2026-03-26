jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');

import { validate } from '../../src/validators/imagesloaded';
import { mockFetchBuffer, mockHashAndFetch, mockReadFileSync, mockSha256Buffer } from '../helpers';

const MOCK_CONFIG = JSON.stringify({ imagesLoaded: { version: '5.0.0' } });

describe('validators/imagesloaded', () => {
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
    expect(result.failures).toContain('imagesloaded.pkgd.min.js');
  });

  test('requests correct CDN URL', async () => {
    await validate();
    expect(mockFetchBuffer).toHaveBeenCalledWith(
      expect.stringContaining('imagesloaded@5.0.0')
    );
  });
});
