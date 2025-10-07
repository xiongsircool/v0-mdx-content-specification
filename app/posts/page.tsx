"use client"

import { useState, useMemo, useEffect } from "react"
import { SimplePostCard } from "@/components/simple-content-card"
import type { ContentItem } from "@/lib/server-markdown-loader"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent } from "@/components/ui/card"
import { Search, Filter, BookOpen, TrendingUp } from "lucide-react"

export default function PostsPage() {
  const [posts, setPosts] = useState<ContentItem[]>([])
  const [loading, setLoading] = useState(true)
  const [searchQuery, setSearchQuery] = useState("")
  const [selectedTags, setSelectedTags] = useState<string[]>([])
  const [sortBy, setSortBy] = useState<"date" | "updated">("date")
  const [postStats, setPostStats] = useState<any>(null)

  // Load posts and stats on component mount
  useEffect(() => {
    const loadData = async () => {
      try {
        // Load posts from API
        const response = await fetch('/api/posts')
        const data = await response.json()
        setPosts(data.posts || [])
      } catch (error) {
        console.error("Failed to load data from API:", error)
        // Fallback: load posts directly from markdown loader
        try {
          const { getPosts } = await import('@/lib/server-markdown-loader')
          const postsData = getPosts()
          const transformedPosts = postsData.map(post => ({
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
          setPosts(transformedPosts)
        } catch (fallbackError) {
          console.error("Fallback loading also failed:", fallbackError)
        }
      } finally {
        setLoading(false)
      }
    }

    loadData()
  }, [])

  // Get all unique tags from posts
  const allTags = useMemo(() => {
    const tags = new Set<string>()
    posts.forEach((post) => {
      post.tags.forEach((tag) => tags.add(tag))
    })
    return Array.from(tags).sort()
  }, [posts])

  const filteredPosts = useMemo(() => {
    return posts
      .filter((post) => {
        // Search filter
        const matchesSearch =
          searchQuery === "" ||
          post.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
          post.excerpt.toLowerCase().includes(searchQuery.toLowerCase()) ||
          post.authors.some((author) => author.toLowerCase().includes(searchQuery.toLowerCase()))

        // Tag filter
        const matchesTags = selectedTags.length === 0 || selectedTags.some((tag) => post.tags.includes(tag))

        return matchesSearch && matchesTags
      })
      .sort((a, b) => {
        if (sortBy === "updated" && a.updatedAt && b.updatedAt) {
          return new Date(b.updatedAt).getTime() - new Date(a.updatedAt).getTime()
        }
        return new Date(b.publishedAt).getTime() - new Date(a.publishedAt).getTime()
      })
  }, [posts, searchQuery, selectedTags, sortBy])

  const toggleTag = (tag: string) => {
    setSelectedTags((prev) => (prev.includes(tag) ? prev.filter((t) => t !== tag) : [...prev, tag]))
  }

  return (
    <div className="min-h-screen bg-background">
      <Navigation />

      <div className="container mx-auto px-4 py-8">
        {/* Header */}
        <div className="text-center mb-12">
          <div className="flex justify-center mb-6">
            <div className="flex h-12 w-12 items-center justify-center rounded-xl bg-primary/10">
              <BookOpen className="h-6 w-6 text-primary" />
            </div>
          </div>
          <h1 className="text-4xl font-bold mb-4">技术推文</h1>
          <p className="text-xl text-muted-foreground max-w-2xl mx-auto">深度技术分享、生物信息学教程和科研经验总结</p>

          {/* Stats */}
          {postStats && (
            <div className="flex justify-center gap-8 mt-6">
              <div className="text-center">
                <div className="text-2xl font-bold text-primary">{postStats.total}</div>
                <div className="text-sm text-muted-foreground">总文章数</div>
              </div>
              <div className="text-center">
                <div className="text-2xl font-bold text-primary">{postStats.recent}</div>
                <div className="text-sm text-muted-foreground">最近30天</div>
              </div>
              <div className="text-center">
                <div className="text-2xl font-bold text-primary">{postStats.topTags?.[0]?.count || 0}</div>
                <div className="text-sm text-muted-foreground">热门标签</div>
              </div>
            </div>
          )}
        </div>

        {/* Filters */}
        <Card className="mb-8">
          <CardContent className="p-6">
            <div className="space-y-4">
              {/* Search */}
              <div className="relative">
                <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-muted-foreground" />
                <Input
                  placeholder="搜索文章标题、内容或作者..."
                  value={searchQuery}
                  onChange={(e) => setSearchQuery(e.target.value)}
                  className="pl-10"
                />
              </div>

              {/* Sort and Tag Filters */}
              <div className="flex flex-col lg:flex-row gap-4">
                {/* Sort Options */}
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <TrendingUp className="h-4 w-4 text-muted-foreground" />
                    <span className="text-sm font-medium">排序方式</span>
                  </div>
                  <div className="flex gap-2">
                    <Button
                      variant={sortBy === "date" ? "default" : "outline"}
                      size="sm"
                      onClick={() => setSortBy("date")}
                    >
                      发布时间
                    </Button>
                    <Button
                      variant={sortBy === "updated" ? "default" : "outline"}
                      size="sm"
                      onClick={() => setSortBy("updated")}
                    >
                      更新时间
                    </Button>
                  </div>
                </div>

                {/* Tag Filters */}
                <div className="flex-1 space-y-3">
                  <div className="flex items-center gap-2">
                    <Filter className="h-4 w-4 text-muted-foreground" />
                    <span className="text-sm font-medium">标签筛选</span>
                  </div>
                  <div className="flex flex-wrap gap-2">
                    {allTags.slice(0, 12).map((tag) => (
                      <Badge
                        key={tag}
                        variant={selectedTags.includes(tag) ? "default" : "outline"}
                        className="cursor-pointer transition-colors"
                        onClick={() => toggleTag(tag)}
                      >
                        {tag}
                      </Badge>
                    ))}
                    {allTags.length > 12 && (
                      <Badge variant="outline" className="cursor-pointer">
                        +{allTags.length - 12} 更多
                      </Badge>
                    )}
                  </div>
                </div>
              </div>

              {/* Clear Filters */}
              {selectedTags.length > 0 && (
                <div>
                  <Button variant="ghost" size="sm" onClick={() => setSelectedTags([])}>
                    清除标签筛选
                  </Button>
                </div>
              )}
            </div>
          </CardContent>
        </Card>

        {/* Results */}
        <div className="mb-6">
          <p className="text-sm text-muted-foreground">
            找到 {filteredPosts.length} 篇文章
            {selectedTags.length > 0 && <span> · 筛选标签: {selectedTags.join(", ")}</span>}
          </p>
        </div>

        {/* Loading State */}
        {loading && (
          <Card>
            <CardContent className="p-12 text-center">
              <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto mb-4"></div>
              <h3 className="text-xl font-semibold mb-2">加载中...</h3>
              <p className="text-muted-foreground">正在获取最新的技术推文</p>
            </CardContent>
          </Card>
        )}

        {/* Posts Grid */}
        {!loading && filteredPosts.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredPosts.map((post) => (
              <SimplePostCard key={post.slug} post={post} />
            ))}
          </div>
        ) : !loading && (
          <Card>
            <CardContent className="p-12 text-center">
              <BookOpen className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-xl font-semibold mb-2">没有找到相关文章</h3>
              <p className="text-muted-foreground mb-4">
                {searchQuery || selectedTags.length > 0 ? "请尝试调整搜索条件或筛选标签" : "目前没有发布任何技术推文"}
              </p>
              {(searchQuery || selectedTags.length > 0) && (
                <Button
                  variant="outline"
                  onClick={() => {
                    setSearchQuery("")
                    setSelectedTags([])
                  }}
                >
                  清除所有筛选
                </Button>
              )}
            </CardContent>
          </Card>
        )}
      </div>

      <Footer />
    </div>
  )
}
