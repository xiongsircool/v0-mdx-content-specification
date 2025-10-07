# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Next.js 15 application for the Shanghai Bioinformatics Center (SBC) student organization, built with v0.app. It's a content-focused website featuring announcements, technical blog posts, learning resources, and meeting information. The project uses a custom markdown loading system instead of Contentlayer2 for content management.

## Development Commands

```bash
# Development
pnpm dev          # Start development server
pnpm build        # Build for production
pnpm start        # Start production server
pnpm lint         # Run ESLint
```

## Technology Stack

- **Framework**: Next.js 15 with App Router
- **Content**: Custom markdown loading system with gray-matter and marked
- **Styling**: Tailwind CSS v4 with shadcn/ui components
- **Language**: TypeScript
- **Package Manager**: pnpm (use pnpm commands, not npm)
- **Deployment**: Vercel (auto-sync from v0.app)

## Content Architecture

The project uses a custom markdown loading system located in `lib/server-markdown-loader.ts` that:

- Parses markdown files using gray-matter for frontmatter
- Converts markdown to HTML using marked with highlight.js for syntax highlighting
- Provides type-safe interfaces for different content types
- Supports search functionality across all content types
- Handles file system operations only on the server side

### Content Types

1. **Announcements** (`content/announcements/**/*.md`)
   - Fields: title, publishedAt, summary, tags (活动/招新/比赛/讲座/通知), priority, expiresAt, attachments
   - Auto-expiring content with priority ordering
   - Supports meeting-specific fields: location, time, speaker, capacity, meetingUrl

2. **Posts** (`content/posts/**/*.md`)
   - Fields: title, publishedAt, updatedAt, excerpt, tags, authors, coverImage
   - Technical blog posts with multiple authors
   - URL structure: `/posts/[slug]` but displayed as `/[slug]` for clean URLs

3. **Learn Resources** (`content/learn/resources/**/*.md`)
   - Fields: title, updatedAt, type (文章/视频/项目/数据集/幻灯片), level (入门/进阶/实战), track, duration, prerequisites
   - Educational content organized by learning tracks
   - Supports external URLs for linked resources

4. **Meetings** (`content/meetings/**/*.md`)
   - Fields: title, date, time, location, eventType, status, speaker, capacity, registered, meetingUrl
   - Meeting information with registration tracking and video conferencing links

## Key Components

### Core Components

- **Navigation** (`components/navigation.tsx`): Responsive header with mobile sheet menu and SimpleGlobalSearch
- **MarkdownContent** (`components/markdown-content.tsx`): Handles markdown rendering with client-side copy button injection
- **ContentCard** (`components/content-card.tsx`): Generic card component for all content types
- **SimpleContentCard** (`components/simple-content-card.tsx`): Lightweight card for posts with clean URL generation

### Search System

- **SimpleGlobalSearch** (`components/simple-global-search.tsx`): Search across all content types
- Search functionality in `lib/server-markdown-loader.ts` provides filtering and ranking
- Real-time search with type-ahead suggestions

## API Routes

The project includes several API routes for dynamic content loading:

- `/api/posts` - Returns posts with fallback to direct markdown loading
- `/api/learn` - Learning resources data
- `/api/sponsors` - Sponsor information with caching
- `/api/literature-cache` - External literature data caching

## URL Structure

- Posts: `/[slug]` (clean URLs, but stored in `/posts/` directory)
- Announcements: `/announcements/[slug]`
- Learn Resources: `/learn/resources/[slug]`
- Meetings: `/meetings/[slug]`

## Development Notes

- All content is written in Chinese (zh-CN)
- Uses Inter font for Latin characters and JetBrains Mono for code
- Vercel Analytics integration for usage tracking
- Auto-synced with v0.app deployments
- Server-side only file system operations for markdown loading
- Client-side copy button injection for code blocks

## Build Configuration

- Next.js config ignores TypeScript and ESLint errors during builds (for v0.app compatibility)
- MDX support enabled with experimental features (mdxRs: true)
- Image optimization disabled for unoptimized images
- Custom markdown parsing with highlight.js integration

## Content Organization

- Content files should be placed in appropriate subdirectories under `content/`
- Use proper frontmatter as defined in `MarkdownMetadata` interface
- Tags for announcements are predefined enums: 活动, 招新, 比赛, 讲座, 通知
- Learning resources should be organized by track (e.g., "生物信息学基础", "Python编程")

## Important Considerations

- The project uses a custom markdown loader instead of Contentlayer2
- URL generation for posts uses clean URLs (`/[slug]`) but files are stored in `/posts/` directory
- Search functionality is implemented server-side in the markdown loader
- Code highlighting uses highlight.js with client-side copy button injection
- All file system operations are wrapped in server-side checks