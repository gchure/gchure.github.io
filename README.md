# Griffin Chure Personal Website

A minimal, tech-forward personal website built with Astro, featuring light/dark themes and content management for publications and research projects.

## Features

- **Modern Design**: Clean, minimal aesthetic inspired by EvolutionaryScale and BioHub
- **Dark/Light Themes**: Automatic theme detection with manual toggle
- **Content Collections**: Markdown-based content management for publications and projects
- **Letter Shuffle Animation**: Eye-catching text animation on the landing page
- **PDF Embedding**: Full PDF embedding on publication pages
- **Responsive Design**: Mobile-friendly layouts throughout
- **Docker Support**: Containerized development environment

## Tech Stack

- **Framework**: Astro 4.x
- **Fonts**: Geist and Geist Mono via Google Fonts
- **Styling**: Scoped CSS with CSS custom properties for theming
- **Content**: Markdown with frontmatter
- **Containerization**: Docker and Docker Compose

## Getting Started

### Prerequisites

- Docker and Docker Compose
- Your profile image (profile.jpg)
- Your publication PDFs (or links to them)

### Installation

1. **Clone the repository** (or you're already here!)

2. **Add your profile image**:
   ```bash
   # Add your image to public/images/profile.jpg
   cp /path/to/your/image.jpg public/images/profile.jpg
   ```

3. **Install dependencies and start the dev server**:
   ```bash
   # Build and start the Docker container
   docker-compose up
   ```

4. **Visit your site**:
   Open http://localhost:4321 in your browser

### Without Docker

If you prefer to run without Docker:

```bash
# Install dependencies
npm install

# Start dev server
npm run dev
```

## Project Structure

```
/
├── public/
│   └── images/          # Static images (add profile.jpg here)
├── src/
│   ├── components/      # Reusable components
│   │   ├── ShuffleText.astro   # Letter shuffle animation
│   │   └── ThemeToggle.astro   # Dark/light mode toggle
│   ├── content/         # Markdown content
│   │   ├── config.ts    # Content collection schemas
│   │   ├── publications/  # Publication markdown files
│   │   └── projects/      # Project markdown files
│   ├── layouts/
│   │   └── BaseLayout.astro    # Main layout with nav
│   ├── pages/           # File-based routing
│   │   ├── index.astro         # Landing page
│   │   ├── about.astro         # About page
│   │   ├── projects.astro      # Projects listing
│   │   ├── teaching.astro      # Teaching page
│   │   └── publications/
│   │       ├── index.astro     # Publications listing
│   │       └── [slug].astro    # Individual publication pages
│   └── styles/
│       └── global.css   # Global styles and theme
├── Dockerfile
├── docker-compose.yml
└── astro.config.mjs
```

## Customizing Content

### Adding Publications

Create a new markdown file in `src/content/publications/`:

```markdown
---
title: "Your Paper Title"
authors: "Author 1, Author 2, Author 3"
journal: "Journal Name"
year: 2024
date: "2024-01-15"
pdf: "https://example.com/paper.pdf"  # Optional
publisher: "https://doi.org/..."       # Optional
doi: "10.1234/example"                 # Optional
---

Your abstract text goes here. This will be displayed on the publication detail page.

You can use multiple paragraphs and markdown formatting.
```

### Adding Projects

Create a new markdown file in `src/content/projects/`:

```markdown
---
title: "Project Title"
description: "Brief one-line description"
order: 1  # Optional: controls display order
---

### Overview

Project details in markdown format...

### Methods

More details...
```

### Updating the About Page

Edit `src/pages/about.astro` to include your personal information, research interests, and background.

### Updating Social Links

Edit the hero section in `src/pages/index.astro` to update your GitHub, Twitter, email, or add other social links.

## Color Scheme

The site uses CSS custom properties for theming. To customize colors, edit `src/styles/global.css`:

```css
:root {
  --accent: #1e90ff;  /* Change to your preferred accent color */
  /* Other color variables... */
}
```

## Building for Production

```bash
# With Docker
docker-compose exec astro npm run build

# Without Docker
npm run build
```

The built site will be in the `dist/` directory.

## Deployment

The built site is static and can be deployed to:
- Netlify
- Vercel
- GitHub Pages
- Cloudflare Pages
- Any static hosting service

Most platforms will auto-detect Astro and use the correct build settings.

## License

This template is provided as-is for your personal use.

## Support

For issues with:
- **Astro**: See [Astro documentation](https://docs.astro.build)
- **This template**: Open an issue or modify as needed

---

Built with [Astro](https://astro.build)
