"use client"

import { useState, useMemo } from "react"
import { allAnnouncements } from "@/lib/contentlayer-mock"
import { AnnouncementCard } from "@/components/content-card"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent } from "@/components/ui/card"
import { Search, Filter, Megaphone } from "lucide-react"

const allTags = ["活动", "招新", "比赛", "讲座", "通知"]

export default function AnnouncementsPage() {
  const [searchQuery, setSearchQuery] = useState("")
  const [selectedTags, setSelectedTags] = useState<string[]>([])
  const [showExpired, setShowExpired] = useState(false)

  const filteredAnnouncements = useMemo(() => {
    return allAnnouncements
      .filter((announcement) => {
        // Search filter
        const matchesSearch =
          searchQuery === "" ||
          announcement.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
          announcement.summary.toLowerCase().includes(searchQuery.toLowerCase())

        // Tag filter
        const matchesTags = selectedTags.length === 0 || selectedTags.some((tag) => announcement.tags.includes(tag))

        // Expired filter
        const matchesExpired = showExpired || !announcement.isExpired

        return matchesSearch && matchesTags && matchesExpired
      })
      .sort((a, b) => {
        // Sort by priority first, then by date
        if (a.priority !== b.priority) {
          return b.priority - a.priority
        }
        return new Date(b.publishedAt).getTime() - new Date(a.publishedAt).getTime()
      })
  }, [searchQuery, selectedTags, showExpired])

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
              <Megaphone className="h-6 w-6 text-primary" />
            </div>
          </div>
          <h1 className="text-4xl font-bold mb-4">公告中心</h1>
          <p className="text-xl text-muted-foreground max-w-2xl mx-auto">
            了解最新的活动信息、招新通知、比赛公告和重要通知
          </p>
        </div>

        {/* Filters */}
        <Card className="mb-8">
          <CardContent className="p-6">
            <div className="space-y-4">
              {/* Search */}
              <div className="relative">
                <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-muted-foreground" />
                <Input
                  placeholder="搜索公告标题或内容..."
                  value={searchQuery}
                  onChange={(e) => setSearchQuery(e.target.value)}
                  className="pl-10"
                />
              </div>

              {/* Tag Filters */}
              <div className="space-y-3">
                <div className="flex items-center gap-2">
                  <Filter className="h-4 w-4 text-muted-foreground" />
                  <span className="text-sm font-medium">标签筛选</span>
                </div>
                <div className="flex flex-wrap gap-2">
                  {allTags.map((tag) => (
                    <Badge
                      key={tag}
                      variant={selectedTags.includes(tag) ? "default" : "outline"}
                      className="cursor-pointer transition-colors"
                      onClick={() => toggleTag(tag)}
                    >
                      {tag}
                    </Badge>
                  ))}
                </div>
              </div>

              {/* Show Expired Toggle */}
              <div className="flex items-center gap-2">
                <Button
                  variant={showExpired ? "default" : "outline"}
                  size="sm"
                  onClick={() => setShowExpired(!showExpired)}
                >
                  {showExpired ? "隐藏过期公告" : "显示过期公告"}
                </Button>
                {selectedTags.length > 0 && (
                  <Button variant="ghost" size="sm" onClick={() => setSelectedTags([])}>
                    清除筛选
                  </Button>
                )}
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Results */}
        <div className="mb-6">
          <p className="text-sm text-muted-foreground">
            找到 {filteredAnnouncements.length} 条公告
            {selectedTags.length > 0 && <span> · 筛选标签: {selectedTags.join(", ")}</span>}
          </p>
        </div>

        {/* Announcements Grid */}
        {filteredAnnouncements.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredAnnouncements.map((announcement) => (
              <AnnouncementCard key={announcement.slug} announcement={announcement} />
            ))}
          </div>
        ) : (
          <Card>
            <CardContent className="p-12 text-center">
              <Megaphone className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-xl font-semibold mb-2">没有找到相关公告</h3>
              <p className="text-muted-foreground mb-4">
                {searchQuery || selectedTags.length > 0 ? "请尝试调整搜索条件或筛选标签" : "目前没有发布任何公告"}
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
