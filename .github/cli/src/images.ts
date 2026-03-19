import fs from 'node:fs';
import path from 'node:path';
import { globSync } from 'glob';
import sharp from 'sharp';
import { logger } from './utils/logger';

export function getOutputPath(inputPath: string, inputBase: string, outputBase: string): string {
  const relative = path.relative(inputBase, inputPath);
  return path.join(outputBase, relative);
}

export function getCompressionSettings(ext: string, opts: { compressionLevel?: number; quality?: number; progressive?: boolean } = {}): Record<string, unknown> {
  const lext = ext.slice(1).toLowerCase();
  if (lext === 'png') {
    return {
      compressionLevel: opts.compressionLevel ?? 6,
      quality: opts.quality ?? 85,
    };
  }
  return {
    quality: opts.quality ?? 85,
    progressive: opts.progressive ?? true,
  };
}

export function getThumbnailSettings(imageWidth: number, opts: { ratio?: number } = {}): { width: number } {
  const ratio = opts.ratio ?? 0.5;
  return { width: Math.round(imageWidth * ratio) };
}

export async function compressImage(inputPath: string, outputPath: string, opts: { compressionLevel?: number; quality?: number; progressive?: boolean } = {}): Promise<void> {
  const ext = path.extname(inputPath).toLowerCase();
  const settings = getCompressionSettings(ext, opts) as { compressionLevel?: number; quality?: number; progressive?: boolean };
  let pipeline = sharp(inputPath);
  if (ext === '.png') {
    pipeline = pipeline.png({ compressionLevel: settings.compressionLevel, quality: settings.quality });
  } else if (ext === '.webp') {
    pipeline = pipeline.webp({ quality: settings.quality });
  } else {
    pipeline = pipeline.jpeg({ quality: settings.quality, progressive: settings.progressive });
  }
  await pipeline.toFile(outputPath);
}

export async function createThumbnail(inputPath: string, outputPath: string, opts: { ratio?: number } = {}): Promise<void> {
  const instance = sharp(inputPath);
  const meta = await instance.metadata();
  const { width } = getThumbnailSettings(meta.width ?? 0, opts);
  await instance.resize({ width }).toFile(outputPath);
}

export async function convertToWebp(inputPath: string, outputPath: string, opts: { quality?: number } = {}): Promise<void> {
  await sharp(inputPath).webp({ quality: opts.quality ?? 85 }).toFile(outputPath);
}

if (require.main === module) {
  const cwd = process.cwd();
  const command = process.argv[2];

  async function runThumbnails(files: string[], imgBase: string, outputBase: string): Promise<void> {
    for (const file of files) {
      const full = path.join(cwd, file);
      const out = getOutputPath(full, imgBase, outputBase);
      fs.mkdirSync(path.dirname(out), { recursive: true });
      await createThumbnail(full, out);
      logger.info(`Thumbnail: ${file} -> ${path.relative(cwd, out)}`);
    }
  }

  async function run(): Promise<void> {
    switch (command) {
      case 'compress': {
        const files = globSync('assets/img/**/*.{png,jpg,jpeg,webp}', { cwd });
        for (const file of files) {
          const full = path.join(cwd, file);
          const tmp = full + '.tmp';
          try {
            await compressImage(full, tmp);
            fs.renameSync(tmp, full);
            logger.info(`Compressed: ${file}`);
          } catch (err) {
            if (fs.existsSync(tmp)) fs.unlinkSync(tmp);
            throw err;
          }
        }
        break;
      }
      case 'thumbnails': {
        const imgBase = path.join(cwd, 'assets/img/feature-img');
        const outputBase = path.join(cwd, 'assets/img/thumbnails/feature-img');
        const files = globSync('assets/img/feature-img/*', { cwd });
        await runThumbnails(files, imgBase, outputBase);
        break;
      }
      case 'thumbnails-all': {
        const imgBase = path.join(cwd, 'assets/img');
        const outputBase = path.join(cwd, 'assets/img/thumbnails');
        const files = [
          ...globSync('assets/img/*.{png,jpg,jpeg,webp}', { cwd }),
          ...globSync('assets/img/!(thumbnails)/**/*.{png,jpg,jpeg,webp}', { cwd }),
        ];
        await runThumbnails(files, imgBase, outputBase);
        break;
      }
      case 'webp': {
        const files = globSync('assets/img/**/*.{png,jpg,jpeg,gif,svg}', { cwd });
        for (const file of files) {
          const full = path.join(cwd, file);
          const out = full.replace(/\.[^.]+$/, '.webp');
          await convertToWebp(full, out);
          logger.info(`WebP: ${file} -> ${path.relative(cwd, out)}`);
        }
        break;
      }
      default:
        logger.error('Usage: npm run <compress|thumbnails|thumbnails-all|webp>');
        process.exit(1);
    }
  }

  run().catch(err => {
    logger.error(String(err));
    process.exit(1);
  });
}
