'use strict';
const fs = require('fs');
const path = require('path');
const { globSync } = require('glob');
const sharp = require('sharp');

function getOutputPath(inputPath, inputBase, outputBase) {
  const relative = path.relative(inputBase, inputPath);
  return path.join(outputBase, relative);
}

function getCompressionSettings(ext, opts = {}) {
  const lext = ext.toLowerCase().replace('.', '');
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

function getThumbnailSettings(imageWidth, opts = {}) {
  const ratio = opts.ratio ?? 0.5;
  return { width: Math.round(imageWidth * ratio) };
}

async function compressImage(inputPath, outputPath, opts = {}) {
  const ext = path.extname(inputPath).toLowerCase();
  const settings = getCompressionSettings(ext, opts);
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

async function createThumbnail(inputPath, outputPath, opts = {}) {
  const meta = await sharp(inputPath).metadata();
  const { width } = getThumbnailSettings(meta.width, opts);
  await sharp(inputPath).resize({ width }).toFile(outputPath);
}

async function convertToWebp(inputPath, outputPath, opts = {}) {
  await sharp(inputPath).webp({ quality: opts.quality ?? 85 }).toFile(outputPath);
}

module.exports = {
  getOutputPath,
  getCompressionSettings,
  getThumbnailSettings,
  compressImage,
  createThumbnail,
  convertToWebp,
};

if (require.main === module) {
  const cwd = process.cwd();
  const command = process.argv[2];

  async function run() {
    if (command === 'compress') {
      const files = globSync('assets/img/**/*.{png,jpg,jpeg,webp}', { cwd });
      for (const file of files) {
        const full = path.join(cwd, file);
        const tmp = full + '.tmp';
        await compressImage(full, tmp);
        fs.renameSync(tmp, full);
        console.log('Compressed:', file);
      }

    } else if (command === 'thumbnails') {
      const inputBase = path.join(cwd, 'assets/img/feature-img');
      const outputBase = path.join(cwd, 'assets/img/thumbnails/feature-img');
      const files = globSync('assets/img/feature-img/*', { cwd });
      for (const file of files) {
        const full = path.join(cwd, file);
        const out = getOutputPath(full, inputBase, outputBase);
        fs.mkdirSync(path.dirname(out), { recursive: true });
        await createThumbnail(full, out);
        console.log('Thumbnail:', file, '->', path.relative(cwd, out));
      }

    } else if (command === 'thumbnails-all') {
      const outputBase = path.join(cwd, 'assets/img/thumbnails');
      const imgBase = path.join(cwd, 'assets/img');
      const rootFiles = globSync('assets/img/*.{png,jpg,jpeg,webp}', { cwd });
      const subFiles = globSync('assets/img/!(thumbnails)/**/*.{png,jpg,jpeg,webp}', { cwd });
      for (const file of [...rootFiles, ...subFiles]) {
        const full = path.join(cwd, file);
        const out = getOutputPath(full, imgBase, outputBase);
        fs.mkdirSync(path.dirname(out), { recursive: true });
        await createThumbnail(full, out);
        console.log('Thumbnail:', file, '->', path.relative(cwd, out));
      }

    } else if (command === 'webp') {
      const files = globSync('assets/img/**/*.{png,jpg,jpeg,gif,svg}', { cwd });
      for (const file of files) {
        const full = path.join(cwd, file);
        const out = full.replace(/\.[^.]+$/, '.webp');
        await convertToWebp(full, out);
        console.log('WebP:', file, '->', path.relative(cwd, out));
      }

    } else {
      console.error('Usage: node images.js <compress|thumbnails|thumbnails-all|webp>');
      process.exit(1);
    }
  }

  run().catch(err => {
    console.error(err);
    process.exit(1);
  });
}
