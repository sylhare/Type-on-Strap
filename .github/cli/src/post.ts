import fs from 'node:fs';
import path from 'node:path';

export function formatDate(date: Date): string {
  return date.toLocaleDateString('en-CA');
}

export function formatTime(date: Date): string {
  return date.toLocaleTimeString('en-CA', { hour12: false });
}

export function createFilename(title: string, date: Date): string {
  return formatDate(date) + '-' + title.replace(/ /g, '-') + '.md';
}

export function createContent(title: string): string {
  return '---\n' +
    'layout: post\n' +
    'title: ' + title + '\n' +
    'tags: []\n' +
    '---';
}

export function createPost(title: string, postsDir: string, date: Date = new Date()): { filename: string; filepath: string; content: string } {
  const filename = createFilename(title, date);
  const filepath = path.join(postsDir, filename);
  const content = createContent(title);
  fs.writeFileSync(filepath, content);
  console.log('[' + formatTime(date) + '] File created: _posts/' + filename);
  return { filename, filepath, content };
}

if (require.main === module) {
  const title = process.argv.slice(2).join(' ');
  if (!title) {
    console.error('Usage: npm run post "Post Title"');
    process.exit(1);
  }
  const postsDir = path.join(process.cwd(), '_posts');
  createPost(title, postsDir);
}
