import { downloadFile, fetchBuffer, fetchJson } from '../../src/utils/http';
import https from 'node:https';
import fs from 'node:fs';
import { EventEmitter } from 'node:events';

jest.mock('https');

const mockHttps = https as jest.Mocked<typeof https>;

function makeMockResponse(statusCode: number, body: string, headers: Record<string, string> = {}) {
  const res = new EventEmitter() as any;
  res.statusCode = statusCode;
  res.headers = headers;
  return { res, emit: (event: string, data?: any) => res.emit(event, data) };
}

describe('http utils', () => {
  describe('fetchBuffer()', () => {
    test('resolves with concatenated chunks on 200', async () => {
      const { res, emit } = makeMockResponse(200, '');
      mockHttps.get = jest.fn((url: any, _options: any, cb: any) => {
        cb(res);
        emit('data', Buffer.from('hello'));
        emit('data', Buffer.from(' world'));
        emit('end');
        return { on: jest.fn() } as any;
      });

      const buf = await fetchBuffer('https://example.com/file');
      expect(buf.toString()).toEqual('hello world');
    });

    test('rejects on non-200 status', async () => {
      const { res, emit } = makeMockResponse(404, '');
      mockHttps.get = jest.fn((url: any, _options: any, cb: any) => {
        cb(res);
        emit('end');
        return { on: jest.fn() } as any;
      });

      await expect(fetchBuffer('https://example.com/missing')).rejects.toThrow('HTTP 404');
    });

    test('follows redirects', async () => {
      const { res: redirect } = makeMockResponse(302, '', { location: 'https://example.com/final' });
      const { res: final, emit: emitFinal } = makeMockResponse(200, '');

      let callCount = 0;
      mockHttps.get = jest.fn((url: any, _options: any, cb: any) => {
        callCount++;
        if (callCount === 1) {
          cb(redirect);
        } else {
          cb(final);
          emitFinal('data', Buffer.from('redirected'));
          emitFinal('end');
        }
        return { on: jest.fn() } as any;
      });

      const buf = await fetchBuffer('https://example.com/redirect');
      expect(buf.toString()).toEqual('redirected');
      expect(callCount).toEqual(2);
    });
  });

  describe('fetchJson()', () => {
    test('parses response as JSON', async () => {
      const { res, emit } = makeMockResponse(200, '');
      mockHttps.get = jest.fn((url: any, _options: any, cb: any) => {
        cb(res);
        emit('data', Buffer.from('{"version":"1.0.0"}'));
        emit('end');
        return { on: jest.fn() } as any;
      });

      const data = await fetchJson<{ version: string }>('https://example.com/api');
      expect(data.version).toEqual('1.0.0');
    });
  });

  describe('downloadFile()', () => {
    test('writes fetched buffer to destination', async () => {
      const { res, emit } = makeMockResponse(200, '');
      mockHttps.get = jest.fn((url: any, _options: any, cb: any) => {
        cb(res);
        emit('data', Buffer.from('file contents'));
        emit('end');
        return { on: jest.fn() } as any;
      });
      const writeFileSyncSpy = jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});

      await downloadFile('https://example.com/file.js', '/dest/file.js');

      expect(writeFileSyncSpy).toHaveBeenCalledWith('/dest/file.js', expect.any(Buffer));
    });
  });
});
