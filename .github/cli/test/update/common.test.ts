jest.mock('../../src/utils/http');
jest.mock('../../src/utils/fs');
jest.mock('../../src/utils/logger');

import { updateSingleVendorFile } from '../../src/update/common';
import { downloadFile } from '../../src/utils/http';
import { updateVendorConfig, updateVersionInFile } from '../../src/utils/fs';

const mockDownloadFile = downloadFile as jest.MockedFunction<typeof downloadFile>;
const mockUpdateVendorConfig = updateVendorConfig as jest.MockedFunction<typeof updateVendorConfig>;
const mockUpdateVersionInFile = updateVersionInFile as jest.MockedFunction<typeof updateVersionInFile>;

describe('update/common', () => {
  beforeEach(() => {
    jest.clearAllMocks();
    mockDownloadFile.mockResolvedValue(undefined);
    mockUpdateVendorConfig.mockImplementation(() => {});
    mockUpdateVersionInFile.mockImplementation(() => {});
  });

  describe('updateSingleVendorFile()', () => {
    test('downloads file using the url builder', async () => {
      await updateSingleVendorFile(
        'MyLib', 'myLib', '/path/to/vendor.js', 'npm run validate:mylib',
        '1.2.3',
        v => `https://cdn.example.com/mylib@${v}/mylib.min.js`
      );
      expect(mockDownloadFile).toHaveBeenCalledWith(
        'https://cdn.example.com/mylib@1.2.3/mylib.min.js',
        '/path/to/vendor.js'
      );
    });

    test('updates vendor config with correct key and version', async () => {
      await updateSingleVendorFile(
        'MyLib', 'myLib', '/path/to/vendor.js', 'npm run validate:mylib',
        '1.2.3',
        v => `https://cdn.example.com/mylib@${v}/mylib.min.js`
      );
      expect(mockUpdateVendorConfig).toHaveBeenCalledWith(
        expect.stringContaining('vendor.config.json'),
        'myLib',
        '1.2.3'
      );
    });

    test('updates head.liquid with display name comment', async () => {
      await updateSingleVendorFile(
        'MyLib', 'myLib', '/path/to/vendor.js', 'npm run validate:mylib',
        '1.2.3',
        v => `https://cdn.example.com/mylib@${v}/mylib.min.js`
      );
      expect(mockUpdateVersionInFile).toHaveBeenCalledWith(
        expect.stringContaining('head.liquid'),
        expect.any(RegExp),
        '<!-- MyLib 1.2.3 -->'
      );
    });
  });
});
