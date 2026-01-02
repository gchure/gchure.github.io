/**
 * Publication-related utilities
 */

import type { CollectionEntry } from 'astro:content';

type Publication = CollectionEntry<'publications'>;

/**
 * Gets publications by their slugs from a list of all publications
 */
export function getPublicationsBySlugs(
  slugs: string[] | undefined,
  allPublications: Publication[]
): Publication[] {
  if (!slugs) return [];
  return slugs
    .map(slug => allPublications.find(pub => pub.slug === slug))
    .filter((pub): pub is Publication => pub !== undefined);
}
