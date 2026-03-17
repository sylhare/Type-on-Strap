import { sha256Buffer, sha256File } from '../../src/utils/hash';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';

describe('hash utils', () => {
  describe('sha256Buffer()', () => {
    test('returns a 64-char hex string', () => {
      const hash = sha256Buffer(Buffer.from('hello'));
      expect(hash).toMatch(/^[a-f0-9]{64}$/);
    });

    test('returns deterministic output for the same input', () => {
      const buf = Buffer.from('test data');
      expect(sha256Buffer(buf)).toEqual(sha256Buffer(buf));
    });

    test('returns different hashes for different inputs', () => {
      expect(sha256Buffer(Buffer.from('a'))).not.toEqual(sha256Buffer(Buffer.from('b')));
    });

    test('matches Node crypto SHA256 output', () => {
      const { createHash } = require('node:crypto');
      const expected = createHash('sha256').update(Buffer.from('abc')).digest('hex');
      expect(sha256Buffer(Buffer.from('abc'))).toEqual(expected);
    });
  });

  describe('sha256File()', () => {
    test('hashes a real temp file and matches sha256Buffer', async () => {
      const tmp = path.join(os.tmpdir(), `hash-test-${Date.now()}.txt`);
      const content = Buffer.from('file content for hashing');
      fs.writeFileSync(tmp, content);
      try {
        const fileHash = await sha256File(tmp);
        const bufHash = sha256Buffer(content);
        expect(fileHash).toEqual(bufHash);
      } finally {
        fs.unlinkSync(tmp);
      }
    });

    test('rejects for non-existent file', async () => {
      await expect(sha256File('/non/existent/file.txt')).rejects.toThrow();
    });
  });
});
