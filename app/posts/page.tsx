"use client"

import { useState, useMemo } from "react"
import { allPosts } from "@/lib/contentlayer-mock"
import { PostCard } from "@/components/content-card"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent } from "@/components/ui/card"
import { Search, Filter, BookOpen, TrendingUp } from "lucide-react"

export default function PostsPage() {
  const [searchQuery, setSearchQuery] = useState("")
  const [selectedTags, setSelectedTags] = useState<string[]>([])
  const [sortBy, setSortBy] = useState<"date" | "updated">("date")

  // Get all unique tags from posts
  const allTags = useMemo(() => {
    const tags = new Set<string>()
    allPosts.forEach((post) => {
      post.tags.forEach((tag) => tags.add(tag))
    })
    return Array.from(tags).sort()
  }, [])

  const filteredPosts = useMemo(() => {
    return allPosts
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
  }, [searchQuery, selectedTags, sortBy])

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

        {/* Posts Grid */}
        {filteredPosts.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredPosts.map((post) => (
              <PostCard key={post.slug} post={post} />
            ))}
          </div>
        ) : (
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
