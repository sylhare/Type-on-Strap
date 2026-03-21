import fs from 'node:fs';
import * as path from 'node:path';
import * as os from 'node:os';
import { randomBytes } from 'node:crypto';
import sharp from 'sharp';
import { compressImage, convertToWebp, createThumbnail } from '../../src/images';

const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'tos-images-'));
const SOURCE_PNG = path.join(TMP, 'source.png');
const SOURCE_WIDTH = 400;
const SOURCE_HEIGHT = 200;

let sourceSize: number;

beforeAll(async () => {
  const pixels = randomBytes(SOURCE_WIDTH * SOURCE_HEIGHT * 3);
  await sharp(Buffer.from(pixels), { raw: { width: SOURCE_WIDTH, height: SOURCE_HEIGHT, channels: 3 } })
    .png({ compressionLevel: 0 })
    .toFile(SOURCE_PNG);
  sourceSize = fs.statSync(SOURCE_PNG).size;
});

afterAll(() => {
  fs.rmSync(TMP, { recursive: true, force: true });
});

describe('source image', () => {
  test('is a valid PNG', async () => {
    const meta = await sharp(SOURCE_PNG).metadata();
    expect(meta.format).toEqual('png');
  });

  test('has correct dimensions', async () => {
    const meta = await sharp(SOURCE_PNG).metadata();
    expect(meta.width).toEqual(SOURCE_WIDTH);
    expect(meta.height).toEqual(SOURCE_HEIGHT);
  });
});

describe('createThumbnail() integration', () => {
  const thumbOut = path.join(TMP, 'thumb.png');

  beforeAll(async () => {
    await createThumbnail(SOURCE_PNG, thumbOut);
  });

  test('creates the thumbnail file', () => {
    expect(fs.existsSync(thumbOut)).toEqual(true);
  });

  test('thumbnail width is exactly 50% of source', async () => {
    const meta = await sharp(thumbOut).metadata();
    expect(meta.width).toEqual(SOURCE_WIDTH / 2);
  });

  test('thumbnail height scales proportionally', async () => {
    const meta = await sharp(thumbOut).metadata();
    expect(meta.height).toEqual(SOURCE_HEIGHT / 2);
  });

  test('thumbnail is a valid image', async () => {
    const meta = await sharp(thumbOut).metadata();
    expect(meta.format).toEqual('png');
  });
});

describe('compressImage() integration', () => {
  const compressedOut = path.join(TMP, 'compressed.png');

  beforeAll(async () => {
    await compressImage(SOURCE_PNG, compressedOut);
  });

  test('creates the compressed output file', () => {
    expect(fs.existsSync(compressedOut)).toEqual(true);
  });

  test('output has the same dimensions as the source', async () => {
    const meta = await sharp(compressedOut).metadata();
    expect(meta.width).toEqual(SOURCE_WIDTH);
    expect(meta.height).toEqual(SOURCE_HEIGHT);
  });

  test('output is a valid PNG', async () => {
    const meta = await sharp(compressedOut).metadata();
    expect(meta.format).toEqual('png');
  });

  test('output is smaller than the source', () => {
    expect(fs.statSync(compressedOut).size).toBeLessThan(sourceSize);
  });
});

describe('convertToWebp() integration', () => {
  const webpOut = path.join(TMP, 'output.webp');

  beforeAll(async () => {
    await convertToWebp(SOURCE_PNG, webpOut);
  });

  test('creates the WebP output file', () => {
    expect(fs.existsSync(webpOut)).toEqual(true);
  });

  test('output format is webp', async () => {
    const meta = await sharp(webpOut).metadata();
    expect(meta.format).toEqual('webp');
  });

  test('output has the same dimensions as the source', async () => {
    const meta = await sharp(webpOut).metadata();
    expect(meta.width).toEqual(SOURCE_WIDTH);
    expect(meta.height).toEqual(SOURCE_HEIGHT);
  });

  test('output is smaller than the source', () => {
    expect(fs.statSync(webpOut).size).toBeLessThan(sourceSize);
  });
});
