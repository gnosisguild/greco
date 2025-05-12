import { defineConfig } from 'vite';
import fs from "vite-plugin-fs";

export default defineConfig({
  server: {
    port: 3000,
    headers: {
      'Cross-Origin-Embedder-Policy': 'require-corp',
      'Cross-Origin-Opener-Policy': 'same-origin',
    }
  },
  plugins: [fs()],
  resolve: {
    alias: {
      '@': '/'
    }
  },
  optimizeDeps: {
    exclude: ['@noir-lang/noir_js']
  }
}); 
