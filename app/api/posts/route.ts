import { NextRequest, NextResponse } from 'next/server'
import { getPosts } from '@/lib/server-markdown-loader'

export async function GET(request: NextRequest) {
  try {
    const { searchParams } = new URL(request.url)
    const query = searchParams.get('q')
    const tag = searchParams.get('tag')
    const author = searchParams.get('author')

    let markdownPosts = getPosts()

    // Apply filters if needed
    if (query) {
      markdownPosts = markdownPosts.filter(post =>
        post.title.toLowerCase().includes(query.toLowerCase()) ||
        post.excerpt?.toLowerCase().includes(query.toLowerCase()) ||
        post.authors?.some(author => author.toLowerCase().includes(query.toLowerCase())) ||
        post.tags?.some(tag => tag.toLowerCase().includes(query.toLowerCase()))
      )
    }

    if (tag) {
      markdownPosts = markdownPosts.filter(post =>
        post.tags?.includes(tag)
      )
    }

    if (author) {
      markdownPosts = markdownPosts.filter(post =>
        post.authors?.includes(author)
      )
    }

    // Transform markdown posts to match expected interface
    const transformedPosts = markdownPosts.map(post => ({
      _id: `post-${post.slug}`,
      _raw: {
        sourceFilePath: `content/posts/${post.slug}.md`,
        sourceFileName: `${post.slug}.md`,
        sourceFileDir: "posts",
        contentType: "markdown",
        flattenedPath: `posts/${post.slug}`
      },
      type: "Post" as const,
      title: post.title,
      publishedAt: post.publishedAt,
      updatedAt: post.updatedAt,
      excerpt: post.excerpt,
      tags: post.tags,
      authors: post.authors,
      coverImage: post.coverImage,
      url: `/posts/${post.slug}`,
      slug: post.slug,
      body: {
        raw: "",
        code: "rendered-markdown-content"
      }
    }))

    return NextResponse.json({ posts: transformedPosts, success: true })
  } catch (error) {
    console.error('Error fetching posts:', error)
    return NextResponse.json(
      { error: 'Failed to fetch posts', success: false },
      { status: 500 }
    )
  }
}

export async function POST(request: NextRequest) {
  try {
    const body = await request.json()
    const { title, excerpt, content, tags, authors, coverImage } = body

    if (!title || !excerpt || !content || !tags || !authors) {
      return NextResponse.json(
        { error: 'Missing required fields', success: false },
        { status: 400 }
      )
    }

    // Note: Creation functionality would need to be added to ServerPostManager
    // For now, return a placeholder response
    return NextResponse.json({
      message: 'Post creation not yet implemented in server manager',
      success: false
    })
  } catch (error) {
    console.error('Error creating post:', error)
    return NextResponse.json(
      { error: 'Failed to create post', success: false },
      { status: 500 }
    )
  }
}