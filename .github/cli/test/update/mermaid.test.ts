jest.mock('../../src/utils/http');
jest.mock('../../src/utils/fs');
jest.mock('../../src/utils/logger');

import { fetchLatestVersion, resolveVersion, updateMermaid } from '../../src/update/mermaid';
import { downloadFile, fetchJson } from '../../src/utils/http';

const mockFetchJson = fetchJson as jest.MockedFunction<typeof fetchJson>;
const mockDownloadFile = downloadFile as jest.MockedFunction<typeof downloadFile>;

describe('update/mermaid', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    mockDownloadFile.mockResolvedValue(undefined);
  });

  describe('fetchLatestVersion()', () => {
    test('returns version from npm registry', async () => {
      mockFetchJson.mockResolvedValue({ version: '11.13.0' });
      expect(await fetchLatestVersion()).toEqual('11.13.0');
    });
  });

  describe('resolveVersion()', () => {
    test('returns the provided argument', async () => {
      expect(await resolveVersion('11.0.0')).toEqual('11.0.0');
    });

    test('fetches latest when no argument provided', async () => {
      mockFetchJson.mockResolvedValue({ version: '11.13.0' });
      expect(await resolveVersion()).toEqual('11.13.0');
    });
  });

  describe('updateMermaid()', () => {
    test('downloads the mermaid JS file from jsdelivr CDN', async () => {
      await updateMermaid('11.13.0');
      expect(mockDownloadFile).toHaveBeenCalledWith(
        expect.stringContaining('mermaid@11.13.0/dist/mermaid.min.js'),
        expect.stringContaining('mermaid.min.js')
      );
    });
  });
});
