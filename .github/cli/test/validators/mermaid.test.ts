jest.mock('../../src/utils/hash');
jest.mock('../../src/utils/http');
jest.mock('fs');

import { validate } from '../../src/validators/mermaid';
import { mockFetchBuffer, mockHashAndFetch, mockReadFileSync, mockSha256Buffer, mockSha256File } from '../helpers';

const MOCK_CONFIG = JSON.stringify({ mermaid: { version: '11.13.0' } });

describe('validators/mermaid', () => {
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
    expect(result.failures).toContain('mermaid.min.js');
  });

  test('fails when local file does not exist', async () => {
    mockSha256File.mockRejectedValueOnce(new Error("ENOENT: no such file or directory, open 'mermaid.min.js'"));
    const result = await validate();
    expect(result.passed).toEqual(false);
    expect(result.failures).toContain('mermaid.min.js');
  });

  test('requests correct CDN URL', async () => {
    await validate();
    expect(mockFetchBuffer).toHaveBeenCalledWith(
      expect.stringContaining('mermaid@11.13.0')
    );
  });
});
