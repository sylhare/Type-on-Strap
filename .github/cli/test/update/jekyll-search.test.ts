jest.mock('../../src/utils/http');
jest.mock('../../src/utils/fs');
jest.mock('../../src/utils/logger');

import { downloadFile, fetchJson } from '../../src/utils/http';
import { updateVendorConfig, updateVersionInFile } from '../../src/utils/fs';
import { fetchLatestVersion, resolveVersion, updateJekyllSearch } from '../../src/update/jekyll-search';

const mockFetchJson = fetchJson as jest.MockedFunction<typeof fetchJson>;
const mockDownloadFile = downloadFile as jest.MockedFunction<typeof downloadFile>;
const mockUpdateVendorConfig = updateVendorConfig as jest.MockedFunction<typeof updateVendorConfig>;
const mockUpdateVersionInFile = updateVersionInFile as jest.MockedFunction<typeof updateVersionInFile>;

describe('update/jekyll-search', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    mockFetchJson.mockResolvedValue({ tag_name: 'v2.1.1' });
    mockDownloadFile.mockResolvedValue(undefined);
    mockUpdateVendorConfig.mockImplementation(() => {});
    mockUpdateVersionInFile.mockImplementation(() => {});
  });

  describe('fetchLatestVersion', () => {
    test('fetches latest version from GitHub API', async () => {
      const version = await fetchLatestVersion();
      expect(mockFetchJson).toHaveBeenCalledWith(
        expect.stringContaining('api.github.com/repos/sylhare/Simple-Jekyll-Search/releases/latest')
      );
      expect(version).toBe('2.1.1');
    });

    test('strips v prefix from tag_name', async () => {
      mockFetchJson.mockResolvedValue({ tag_name: 'v2.1.2' });
      const version = await fetchLatestVersion();
      expect(version).toBe('2.1.2');
    });
  });

  describe('resolveVersion', () => {
    test('uses provided version without fetching', async () => {
      const version = await resolveVersion('2.1.2');
      expect(mockFetchJson).not.toHaveBeenCalled();
      expect(version).toBe('2.1.2');
    });

    test('strips v prefix from provided version', async () => {
      const version = await resolveVersion('v2.1.2');
      expect(version).toBe('2.1.2');
    });

    test('fetches latest when no version provided', async () => {
      await resolveVersion();
      expect(mockFetchJson).toHaveBeenCalled();
    });
  });

  describe('updateJekyllSearch', () => {
    test('downloads from correct GitHub raw URL', async () => {
      await updateJekyllSearch('2.1.2');
      expect(mockDownloadFile).toHaveBeenCalledWith(
        expect.stringContaining('Simple-Jekyll-Search/v2.1.2/dest/simple-jekyll-search.min.js'),
        expect.any(String)
      );
    });

    test('updates vendor.config.json with simpleJekyllSearch key', async () => {
      await updateJekyllSearch('2.1.2');
      expect(mockUpdateVendorConfig).toHaveBeenCalledWith(
        expect.stringContaining('vendor.config.json'),
        'simpleJekyllSearch',
        '2.1.2'
      );
    });

    test('updates head.liquid comment', async () => {
      await updateJekyllSearch('2.1.2');
      expect(mockUpdateVersionInFile).toHaveBeenCalledWith(
        expect.stringContaining('head.liquid'),
        expect.any(RegExp),
        '<!-- Simple Jekyll Search 2.1.2 -->'
      );
    });
  });
});
