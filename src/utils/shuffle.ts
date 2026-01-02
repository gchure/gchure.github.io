/**
 * Shared text shuffle/scramble animation utilities
 */

// Unified character set for shuffle animations
export const SHUFFLE_CHARS = '0123456789@#$%&*↗{}[]()!<>_+=';

/**
 * Animates text by scrambling characters then revealing them one by one
 * @param element - The element containing the text to animate
 * @param originalText - The original text to reveal
 * @param options - Animation options
 */
export function shuffleText(
  element: HTMLElement,
  originalText: string,
  options: {
    chars?: string;
    speed?: number;
    increment?: number;
  } = {}
) {
  const {
    chars = SHUFFLE_CHARS,
    speed = 30,
    increment = 0.5,
  } = options;

  let iteration = 0;

  const interval = setInterval(() => {
    element.textContent = originalText
      .split('')
      .map((char, index) => {
        if (char === ' ') return ' ';
        if (index < iteration) {
          return originalText[index];
        }
        return chars[Math.floor(Math.random() * chars.length)];
      })
      .join('');

    if (iteration >= originalText.length) {
      clearInterval(interval);
    }

    iteration += increment;
  }, speed);

  return interval;
}

/**
 * Attaches shuffle animation to elements on mouseenter
 * @param selector - CSS selector for elements to animate
 * @param textSelector - CSS selector for the text element within each item
 * @param textAttr - Attribute containing the original text (default: 'data-text')
 */
export function attachShuffleOnHover(
  selector: string,
  textSelector: string,
  textAttr: string = 'data-text'
) {
  document.querySelectorAll(selector).forEach((element) => {
    const originalText = element.getAttribute(textAttr) || '';
    const textElement = element.querySelector(textSelector) as HTMLElement;

    if (!textElement) return;

    element.addEventListener('mouseenter', () => {
      shuffleText(textElement, originalText);
    });
  });
}
