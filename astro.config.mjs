import { defineConfig } from 'astro/config';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import nordLight from './nord-light.json';

export default defineConfig({
  site: 'https://gchure.bio',
  markdown: {
    remarkPlugins: [remarkMath],
    rehypePlugins: [rehypeKatex],
    shikiConfig: {
      themes: {
        light: nordLight,
        dark: nordLight,
      },
    }
  }
});
