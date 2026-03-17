import fs from 'node:fs';
import path from 'node:path';
import { formatDate, formatTime, createFilename, createContent, createPost } from '../src/post';

describe('post.ts', () => {
  const fixedDate = new Date('2024-01-15T10:30:00');

  describe('formatDate()', () => {
    test('returns YYYY-MM-DD format', () => {
      expect(formatDate(fixedDate)).toEqual('2024-01-15');
    });

    test('pads month and day with leading zeros', () => {
      const date = new Date('2024-03-05T00:00:00');
      expect(formatDate(date)).toMatch(/^\d{4}-\d{2}-\d{2}$/);
    });
  });

  describe('formatTime()', () => {
    test('returns HH:MM:SS format', () => {
      expect(formatTime(fixedDate)).toMatch(/^\d{2}:\d{2}:\d{2}$/);
    });
  });

  describe('createFilename()', () => {
    test('slugifies title by replacing spaces with hyphens', () => {
      expect(createFilename('Hello World', fixedDate)).toEqual('2024-01-15-Hello-World.md');
    });

    test('prepends date to title', () => {
      const filename = createFilename('My Post Title', fixedDate);
      expect(filename).toEqual('2024-01-15-My-Post-Title.md');
    });

    test('has .md extension', () => {
      expect(createFilename('Test', fixedDate)).toMatch(/\.md$/);
    });
  });

  describe('createContent()', () => {
    test('contains layout: post', () => {
      expect(createContent('My Title')).toContain('layout: post');
    });

    test('contains the post title', () => {
      expect(createContent('My Title')).toContain('title: My Title');
    });

    test('contains empty tags array', () => {
      expect(createContent('My Title')).toContain('tags: []');
    });

    test('starts and ends with frontmatter delimiters', () => {
      const content = createContent('Test');
      expect(content.startsWith('---')).toEqual(true);
      expect(content.endsWith('---')).toEqual(true);
    });
  });

  describe('createPost()', () => {
    test('writes file with correct path', () => {
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      jest.spyOn(console, 'log').mockImplementation(() => {});

      const postsDir = '/tmp/posts';
      createPost('Hello World', postsDir, fixedDate);

      expect(fs.writeFileSync).toHaveBeenCalledWith(
        path.join(postsDir, '2024-01-15-Hello-World.md'),
        expect.any(String)
      );
    });

    test('writes correct frontmatter content', () => {
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      jest.spyOn(console, 'log').mockImplementation(() => {});

      const result = createPost('Hello World', '/tmp/posts', fixedDate);

      expect(result.content).toContain('layout: post');
      expect(result.content).toContain('title: Hello World');
      expect(result.content).toContain('tags: []');
    });

    test('returns filename, filepath, and content', () => {
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      jest.spyOn(console, 'log').mockImplementation(() => {});

      const postsDir = '/tmp/posts';
      const result = createPost('Hello World', postsDir, fixedDate);

      expect(result.filename).toEqual('2024-01-15-Hello-World.md');
      expect(result.filepath).toEqual(path.join(postsDir, '2024-01-15-Hello-World.md'));
      expect(typeof result.content).toEqual('string');
    });

    test('logs the created file name', () => {
      jest.spyOn(fs, 'writeFileSync').mockImplementation(() => {});
      const logSpy = jest.spyOn(console, 'log').mockImplementation(() => {});

      createPost('Hello World', '/tmp/posts', fixedDate);

      expect(logSpy).toHaveBeenCalledWith(
        expect.stringContaining('2024-01-15-Hello-World.md')
      );
    });
  });
});
