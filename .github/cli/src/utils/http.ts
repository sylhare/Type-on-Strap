import * as https from 'node:https';
import type { IncomingMessage } from 'node:http';
import * as http from 'node:http';
import * as fs from 'node:fs';

function getProtocol(url: string) {
  return url.startsWith('https') ? https : http;
}

function buildRequestOptions(url: string): https.RequestOptions {
  const headers: Record<string, string> = { 'User-Agent': 'type-on-strap-cli' };
  const token = process.env['GITHUB_TOKEN'];
  if (token && url.includes('api.github.com')) {
    headers['Authorization'] = `Bearer ${token}`;
  }
  return { headers };
}

export function fetchBuffer(url: string): Promise<Buffer> {
  return new Promise((resolve, reject) => {
    getProtocol(url).get(url, buildRequestOptions(url), (res: IncomingMessage) => {
      if (res.statusCode === 301 || res.statusCode === 302) {
        const location = res.headers.location;
        if (!location) return reject(new Error(`Redirect with no location for ${url}`));
        return fetchBuffer(location).then(resolve).catch(reject);
      }
      if (res.statusCode !== 200) {
        return reject(new Error(`HTTP ${res.statusCode} for ${url}`));
      }
      const chunks: Buffer[] = [];
      res.on('data', (chunk: Buffer) => chunks.push(chunk));
      res.on('end', () => resolve(Buffer.concat(chunks)));
      res.on('error', reject);
    }).on('error', reject);
  });
}

export async function fetchJson<T>(url: string): Promise<T> {
  const buf = await fetchBuffer(url);
  return JSON.parse(buf.toString('utf8')) as T;
}

export async function downloadFile(url: string, dest: string): Promise<void> {
  const buf = await fetchBuffer(url);
  fs.writeFileSync(dest, buf);
}
