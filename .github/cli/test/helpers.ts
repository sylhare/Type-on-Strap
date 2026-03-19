import fs from 'node:fs';
import { sha256File, sha256Buffer } from '../src/utils/hash';
import { fetchBuffer } from '../src/utils/http';

export const realReadFileSync = (jest.requireActual('node:fs') as typeof import('node:fs')).readFileSync;

export const MOCK_HASH = 'a'.repeat(64);

export const mockSha256File = sha256File as jest.MockedFunction<typeof sha256File>;
export const mockSha256Buffer = sha256Buffer as jest.MockedFunction<typeof sha256Buffer>;
export const mockFetchBuffer = fetchBuffer as jest.MockedFunction<typeof fetchBuffer>;

/**
 * Mocks fs.readFileSync so vendor.config.json returns `vendorConfig`,
 * and an optional `extra` callback can intercept other paths.
 */
export function mockReadFileSync(
  vendorConfig: string,
  extra?: (filePath: string, options?: unknown) => string | Buffer | undefined
): void {
  const mockFs = fs as jest.Mocked<typeof fs>;
  mockFs.readFileSync = jest.fn().mockImplementation((filePath: unknown, options?: unknown) => {
    if (String(filePath).endsWith('vendor.config.json')) return vendorConfig;
    const result = extra?.(String(filePath), options);
    if (result !== undefined) return result;
    return realReadFileSync(filePath as any, options as any);
  }) as jest.MockedFunction<typeof fs.readFileSync>;
}

/**
 * Sets up the standard hash/fetch mocks used by single-file validators.
 * Call this in beforeEach alongside mockReadFileSync.
 */
export function mockHashAndFetch(): void {
  mockSha256File.mockResolvedValue(MOCK_HASH);
  mockSha256Buffer.mockReturnValue(MOCK_HASH);
  mockFetchBuffer.mockResolvedValue(Buffer.from('remote content'));
}
