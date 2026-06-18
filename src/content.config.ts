import { defineCollection, z } from 'astro:content';
import { glob } from 'astro/loaders';

const publications = defineCollection({
  loader: glob({ pattern: '**/*.md', base: './src/content/publications' }),
  schema: z.object({
    title: z.string(),
    description: z.string().optional(),
    summary: z.string().optional(),
    abstract: z.string().optional(),
    date: z.string(),
    authors: z.array(z.object({
      name: z.string(),
      equal: z.boolean().optional(),
      corr: z.boolean().optional(),
    })),
    tags: z.array(z.string()).optional(),
    journal: z.string(),
    DOI: z.string().optional(),
    image: z.string().optional(),
    links: z.array(z.object({
      name: z.string(),
      url: z.string(),
    })).optional(),
  }),
});

const research = defineCollection({
  loader: glob({ pattern: '**/*.md', base: './src/content/research' }),
  schema: z.object({
    title: z.string(),
    description: z.string(),
    order: z.number().optional(),
    image: z.string().optional(),
    publications: z.array(z.string()).optional(), // Array of publication slugs
  }),
});

const software = defineCollection({
  loader: glob({ pattern: '**/*.md', base: './src/content/software' }),
  schema: z.object({
    title: z.string(),
    description: z.string(),
    order: z.number().optional(),
    image: z.string().optional(),
    stack: z.array(z.string()).optional(), // Tech stack (e.g., Python, NumPy, etc.)
    links: z.array(z.object({
      name: z.string(),
      url: z.string(),
    })).optional(),
    publications: z.array(z.string()).optional(), // Array of publication slugs
  }),
});

const writing = defineCollection({
  loader: glob({ pattern: '**/*.md', base: './src/content/writing' }),
  schema: z.object({
    title: z.string(),
    description: z.string(),
    date: z.string(),
    tags: z.array(z.string()),
    draft: z.boolean().optional(),
  }),
});

export const collections = {
  publications,
  research,
  software,
  writing,
};
