import { randomBytes } from 'crypto';

const NONCE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
const DEFAULT_NONCE_LENGTH = 32;
const BYTE_LIMIT = 256 - (256 % NONCE_CHARS.length);

export function generateNonce(length: number = DEFAULT_NONCE_LENGTH): string {
  if (!Number.isInteger(length) || length <= 0) {
    throw new Error('Nonce length must be a positive integer');
  }

  let nonce = '';
  while (nonce.length < length) {
    const bytes = randomBytes(length);
    for (const byte of bytes) {
      if (byte >= BYTE_LIMIT) {
        continue;
      }
      nonce += NONCE_CHARS[byte % NONCE_CHARS.length];
      if (nonce.length === length) {
        break;
      }
    }
  }
  return nonce;
}
