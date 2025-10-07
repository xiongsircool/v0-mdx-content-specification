import { NextResponse } from 'next/server'
import type { NextRequest } from 'next/server'

export function middleware(request: NextRequest) {
  const { pathname } = request.nextUrl

  // Static routes that should not be rewritten
  const staticRoutes = ['about', 'learn', 'announcements', 'meetings', 'posts', 'search', 'admin', 'test-video', 'api', '_next']

  // Skip if it's an API route, static file, or root
  if (pathname === '/' ||
      pathname.startsWith('/api/') ||
      pathname.startsWith('/_next/') ||
      pathname.startsWith('/announcements/') ||
      pathname.startsWith('/learn/') ||
      pathname.startsWith('/meetings/') ||
      pathname.startsWith('/posts/') ||
      // Skip static files (files with extensions)
      /\.[a-zA-Z0-9]+$/.test(pathname)) {
    return NextResponse.next()
  }

  // Check if this is a clean URL for a content item
  const pathSegments = pathname.split('/').filter(Boolean)

  if (pathSegments.length === 1 && !staticRoutes.includes(pathSegments[0])) {
    const slug = pathSegments[0]

    // Parse the slug to determine content type
    if (slug.startsWith('announcements/') || slug.startsWith('posts/') || slug.startsWith('learn/resources/') || slug.startsWith('meetings/')) {
      // This is a full path slug - extract the type and rewrite accordingly
      const parts = slug.split('/')

      if (parts[0] === 'announcements') {
        return NextResponse.rewrite(
          new URL(`/announcements/${parts.slice(1).join('/')}`, request.url)
        )
      } else if (parts[0] === 'posts') {
        return NextResponse.rewrite(
          new URL(`/posts/${parts.slice(1).join('/')}`, request.url)
        )
      } else if (parts[0] === 'learn' && parts[1] === 'resources') {
        return NextResponse.rewrite(
          new URL(`/learn/resources/${parts.slice(2).join('/')}`, request.url)
        )
      } else if (parts[0] === 'meetings') {
        return NextResponse.rewrite(
          new URL(`/meetings/${parts.slice(1).join('/')}`, request.url)
        )
      }
    } else {
      // This is a simple slug - try to determine content type
      const datePattern = /^\d{4}-\d{2}-\d{2}/

      if (datePattern.test(slug)) {
        // Date-based slugs are likely posts or announcements
        // Try posts first
        return NextResponse.rewrite(
          new URL(`/posts/${slug}`, request.url)
        )
      } else {
        // Non-date-based slugs - default to posts
        return NextResponse.rewrite(
          new URL(`/posts/${slug}`, request.url)
        )
      }
    }
  }

  return NextResponse.next()
}

export const config = {
  matcher: [
    /*
     * Match all request paths except for the ones starting with:
     * - api (API routes)
     * - _next/static (static files)
     * - _next/image (image optimization files)
     * - announcements (announcement pages)
     * - learn (learning resources)
     * - meetings (meeting pages)
     * - posts (post pages)
     */
    '/((?!api|_next/static|_next/image|announcements|learn|meetings|posts).*)',
  ],
}