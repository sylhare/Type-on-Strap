jest.mock('../../src/utils/http');
jest.mock('../../src/utils/fs');

import { fetchLatestVersion, resolveVersion, updateMermaid } from '../../src/update/mermaid';
import { downloadFile, fetchJson } from '../../src/utils/http';
import { updateVendorConfig, updateVersionInFile } from '../../src/utils/fs';

const mockFetchJson = fetchJson as jest.MockedFunction<typeof fetchJson>;
const mockDownloadFile = downloadFile as jest.MockedFunction<typeof downloadFile>;
const mockUpdateVendorConfig = updateVendorConfig as jest.MockedFunction<typeof updateVendorConfig>;
const mockUpdateVersionInFile = updateVersionInFile as jest.MockedFunction<typeof updateVersionInFile>;

describe('update/mermaid', () => {
  beforeEach(() => {
    mockDownloadFile.mockResolvedValue(undefined as unknown as void);
    mockUpdateVendorConfig.mockImplementation(() => {
    });
    mockUpdateVersionInFile.mockImplementation(() => {
    });
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
    test('downloads the mermaid JS file from CDN', async () => {
      await updateMermaid('11.13.0');
      expect(mockDownloadFile).toHaveBeenCalledWith(
        expect.stringContaining('mermaid@11.13.0'),
        expect.stringContaining('mermaid.min.js')
      );
    });

    test('updates vendor config with new version', async () => {
      await updateMermaid('11.13.0');
      expect(mockUpdateVendorConfig).toHaveBeenCalledWith(
        expect.any(String), 'mermaid', '11.13.0'
      );
    });

    test('updates head.liquid with new version', async () => {
      await updateMermaid('11.13.0');
      expect(mockUpdateVersionInFile).toHaveBeenCalledWith(
        expect.stringContaining('head.liquid'),
        expect.any(RegExp),
        '<!-- Mermaid 11.13.0 -->'
      );
    });
  });
});
