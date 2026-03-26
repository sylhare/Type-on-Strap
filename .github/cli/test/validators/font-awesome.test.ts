jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');

import { validate } from '../../src/validators/font-awesome';
import { mockFetchBuffer, mockHashAndFetch, mockReadFileSync, mockSha256Buffer, mockSha256File } from '../helpers';

const MOCK_CONFIG = JSON.stringify({ fontAwesome: { version: '6.7.2' } });

describe('validators/font-awesome', () => {
  beforeEach(() => {
    mockReadFileSync(MOCK_CONFIG);
    mockHashAndFetch();
  });

  test('passes when all 25 files match (6 fonts + 19 SCSS)', async () => {
    const result = await validate();
    expect(result.passed).toEqual(true);
    expect(result.failures).toHaveLength(0);
    expect(mockFetchBuffer).toHaveBeenCalledTimes(25);
  });

  test('fails when any file hash does not match', async () => {
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
