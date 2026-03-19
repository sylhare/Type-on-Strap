import * as path from 'node:path';
import sharp from 'sharp';
import {
  compressImage,
  convertToWebp,
  getCompressionSettings,
  getOutputPath,
  getThumbnailSettings,
} from '../src/images';

jest.mock('sharp', () => jest.fn());
jest.mock('glob', () => ({ globSync: jest.fn() }));

const mockSharp = sharp as jest.Mock;

describe('images.ts', () => {
  let mockPipeline: {
    jpeg: jest.Mock;
    png: jest.Mock;
    webp: jest.Mock;
    resize: jest.Mock;
    metadata: jest.Mock;
    toFile: jest.Mock;
  };

  beforeEach(() => {
    mockPipeline = {
      jpeg: jest.fn().mockReturnThis(),
      png: jest.fn().mockReturnThis(),
      webp: jest.fn().mockReturnThis(),
      resize: jest.fn().mockReturnThis(),
      metadata: jest.fn().mockResolvedValue({ width: 1000, height: 600 }),
      toFile: jest.fn().mockResolvedValue({}),
    };
    mockSharp.mockReturnValue(mockPipeline);
  });

  describe('getOutputPath()', () => {
    test('remaps input path from inputBase to outputBase', () => {
      const result = getOutputPath(
        '/project/assets/img/feature-img/photo.jpg',
        '/project/assets/img/feature-img',
        '/project/assets/img/thumbnails/feature-img'
      );
      expect(result).toEqual('/project/assets/img/thumbnails/feature-img/photo.jpg');
    });

    test('preserves filename in output path', () => {
      const result = getOutputPath('/img/a/b.png', '/img/a', '/out');
      expect(path.basename(result)).toEqual('b.png');
    });

    test('preserves subdirectory structure', () => {
      const result = getOutputPath('/img/sub/dir/file.jpg', '/img', '/out');
      expect(result).toEqual('/out/sub/dir/file.jpg');
    });
  });

  describe('getCompressionSettings()', () => {
    test('returns compressionLevel for PNG', () => {
      const settings = getCompressionSettings('.png');
      expect(settings).toHaveProperty('compressionLevel');
      expect(typeof settings['compressionLevel']).toEqual('number');
    });

    test('returns quality for JPEG', () => {
      const settings = getCompressionSettings('.jpg');
      expect(settings).toHaveProperty('quality');
      expect(typeof settings['quality']).toEqual('number');
    });

    test('returns quality for WebP', () => {
      const settings = getCompressionSettings('.webp');
      expect(settings).toHaveProperty('quality');
    });

    test('accepts custom quality option', () => {
      const settings = getCompressionSettings('.jpg', { quality: 70 });
      expect(settings['quality']).toEqual(70);
    });

    test('accepts custom compressionLevel for PNG', () => {
      const settings = getCompressionSettings('.png', { compressionLevel: 9 });
      expect(settings['compressionLevel']).toEqual(9);
    });
  });

  describe('getThumbnailSettings()', () => {
    test('returns width as a number', () => {
      const result = getThumbnailSettings(1000);
      expect(typeof result.width).toEqual('number');
    });

    test('returns 50% of imageWidth by default', () => {
      expect(getThumbnailSettings(1000).width).toEqual(500);
    });

    test('rounds the width to an integer', () => {
      const result = getThumbnailSettings(101);
      expect(Number.isInteger(result.width)).toEqual(true);
    });

    test('respects custom ratio option', () => {
      expect(getThumbnailSettings(1000, { ratio: 0.25 }).width).toEqual(250);
    });
  });

  describe('compressImage()', () => {
    test('uses jpeg pipeline for .jpg files', async () => {
      await compressImage('/input/photo.jpg', '/output/photo.jpg');
      expect(mockPipeline.jpeg).toHaveBeenCalled();
      expect(mockPipeline.png).not.toHaveBeenCalled();
    });

    test('uses png pipeline for .png files', async () => {
      await compressImage('/input/photo.png', '/output/photo.png');
      expect(mockPipeline.png).toHaveBeenCalled();
      expect(mockPipeline.jpeg).not.toHaveBeenCalled();
    });

    test('uses webp pipeline for .webp files', async () => {
      await compressImage('/input/photo.webp', '/output/photo.webp');
      expect(mockPipeline.webp).toHaveBeenCalled();
    });
  });

  describe('convertToWebp()', () => {
    test('uses custom quality when provided', async () => {
      await convertToWebp('/input/photo.jpg', '/output/photo.webp', { quality: 90 });
      expect(mockPipeline.webp).toHaveBeenCalledWith({ quality: 90 });
    });
  });
});
