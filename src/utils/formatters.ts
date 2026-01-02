/**
 * Shared formatting utilities for the website
 */

export interface Author {
  name: string;
  equal?: boolean;
  corr?: boolean;
}

/**
 * Emphasizes "Griffin Chure" in an author string with bold styling
 */
export function emphasizeAuthor(authorString: string): string {
  return authorString.replace(
    /(Griffin Chure[†*]*)/g,
    '<strong class="emphasized-author">$1</strong>'
  );
}

/**
 * Formats a date string to "Mon DD, YYYY" format
 */
export function formatDate(dateString: string): string {
  const date = new Date(dateString);
  return date.toLocaleDateString('en-US', {
    year: 'numeric',
    month: 'short',
    day: 'numeric'
  });
}

/**
 * Formats an array of authors into a comma-separated string with markers
 * † = equal contribution, * = corresponding author
 */
export function formatAuthors(authors: Author[]): string {
  return authors
    .map(a => {
      let name = a.name;
      if (a.equal) name += '†';
      if (a.corr) name += '*';
      return name;
    })
    .join(', ');
}

/**
 * Formats authors with HTML markers (superscript)
 */
export function formatAuthorsHtml(authors: Author[]): string {
  return authors
    .map(a => {
      let name = a.name;
      if (a.equal) name += '<sup>†</sup>';
      if (a.corr) name += '<sup>*</sup>';
      return name;
    })
    .join(', ');
}
