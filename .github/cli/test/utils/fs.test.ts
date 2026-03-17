import { updateVersionInFile, updateVendorConfig } from '../../src/utils/fs';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';

describe('fs utils', () => {
  describe('updateVersionInFile()', () => {
    test('replaces matching pattern in file', () => {
      const tmp = path.join(os.tmpdir(), `fs-test-${Date.now()}.txt`);
      fs.writeFileSync(tmp, 'VERSION="1.0.0"\nother content');
      try {
        updateVersionInFile(tmp, /VERSION="[^"]+"/, 'VERSION="2.0.0"');
        expect(fs.readFileSync(tmp, 'utf8')).toContain('VERSION="2.0.0"');
      } finally {
        fs.unlinkSync(tmp);
      }
    });

    test('throws if pattern not found', () => {
      const tmp = path.join(os.tmpdir(), `fs-test-${Date.now()}.txt`);
      fs.writeFileSync(tmp, 'no match here');
      try {
        expect(() =>
          updateVersionInFile(tmp, /VERSION="[^"]+"/, 'VERSION="2.0.0"')
        ).toThrow('Pattern not found');
      } finally {
        fs.unlinkSync(tmp);
      }
    });

    test('does not rewrite file if content is unchanged', () => {
      const tmp = path.join(os.tmpdir(), `fs-test-${Date.now()}.txt`);
      fs.writeFileSync(tmp, 'VERSION="1.0.0"');
      const spy = jest.spyOn(fs, 'writeFileSync');
      try {
        // Pattern matches but replacement produces same content
        updateVersionInFile(tmp, /VERSION="1.0.0"/, 'VERSION="1.0.0"');
        expect(spy).not.toHaveBeenCalledWith(tmp, expect.anything());
      } finally {
        fs.unlinkSync(tmp);
      }
    });
  });

  describe('updateVendorConfig()', () => {
    test('updates a version key in the JSON config', () => {
      const tmp = path.join(os.tmpdir(), `vendor-test-${Date.now()}.json`);
      fs.writeFileSync(tmp, JSON.stringify({ mermaid: { version: '11.0.0' } }, null, 2) + '\n');
      try {
        updateVendorConfig(tmp, 'mermaid', '11.13.0');
        const updated = JSON.parse(fs.readFileSync(tmp, 'utf8'));
        expect(updated.mermaid.version).toBe('11.13.0');
      } finally {
        fs.unlinkSync(tmp);
      }
    });

    test('preserves other keys', () => {
      const tmp = path.join(os.tmpdir(), `vendor-test-${Date.now()}.json`);
      fs.writeFileSync(tmp, JSON.stringify({ katex: { version: '0.16.0' }, mermaid: { version: '11.0.0' } }, null, 2) + '\n');
      try {
        updateVendorConfig(tmp, 'mermaid', '11.13.0');
        const updated = JSON.parse(fs.readFileSync(tmp, 'utf8'));
        expect(updated.katex.version).toBe('0.16.0');
        expect(updated.mermaid.version).toBe('11.13.0');
      } finally {
        fs.unlinkSync(tmp);
      }
    });
  });
});
