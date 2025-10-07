import { getAnnouncements } from "@/lib/server-markdown-loader"
import { ContentCard } from "@/components/content-card"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent } from "@/components/ui/card"
import { Search, Filter, Megaphone } from "lucide-react"
import type { ContentItem } from "@/lib/server-markdown-loader"

interface AnnouncementsPageProps {
  searchParams: {
    q?: string
    tags?: string
    expired?: string
  }
}

export default function AnnouncementsPage({ searchParams }: AnnouncementsPageProps) {
  const announcements = getAnnouncements()

  // Get filter parameters
  const searchQuery = searchParams.q || ''
  const selectedTags = searchParams.tags ? searchParams.tags.split(',') : []
  const showExpired = searchParams.expired === 'true'

  // Filter announcements
  const filteredAnnouncements = announcements
    .filter((announcement) => {
      // Search filter
      const matchesSearch =
        searchQuery === "" ||
        announcement.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
        announcement.summary?.toLowerCase().includes(searchQuery.toLowerCase())

      // Tag filter
      const matchesTags = selectedTags.length === 0 ||
        selectedTags.some((tag) => announcement.tags?.includes(tag))

      // Expired filter
      const matchesExpired = showExpired || !announcement.isExpired

      return matchesSearch && matchesTags && matchesExpired
    })
    .sort((a, b) => {
      // Sort by priority first, then by date
      if (a.priority !== b.priority) {
        return (b.priority || 0) - (a.priority || 0)
      }
      return new Date(b.publishedAt || '').getTime() - new Date(a.publishedAt || '').getTime()
    })

  // Get all available tags
  const allTags = Array.from(new Set(
    announcements.flatMap(announcement => announcement.tags || [])
  )).sort()

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
          <h1 className="text-4xl font-bold mb-4">联盟活动</h1>
          <p className="text-xl text-muted-foreground max-w-2xl mx-auto">
            了解联盟最新的活动安排、招新信息、技术分享会和合作机会
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
                  placeholder="搜索活动标题或内容..."
                  defaultValue={searchQuery}
                  name="q"
                  className="pl-10"
                />
              </div>

              {/* Tag Filters */}
              <div className="space-y-3">
                <div className="flex items-center gap-2">
                  <Filter className="h-4 w-4 text-muted-foreground" />
                  <span className="text-sm font-medium">活动类型</span>
                </div>
                <div className="flex flex-wrap gap-2">
                  {allTags.map((tag) => (
                    <Badge
                      key={tag}
                      variant={selectedTags.includes(tag) ? "default" : "outline"}
                      className="cursor-pointer transition-colors"
                    >
                      {tag}
                    </Badge>
                  ))}
                </div>
              </div>

              {/* Results count */}
              <div className="text-sm text-muted-foreground">
                找到 {filteredAnnouncements.length} 个活动
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Results Grid */}
        {filteredAnnouncements.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredAnnouncements.map((announcement) => (
              <ContentCard key={announcement.slug} item={announcement} type="announcement" />
            ))}
          </div>
        ) : (
          <Card>
            <CardContent className="p-12 text-center">
              <Megaphone className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-xl font-semibold mb-2">没有找到相关活动</h3>
              <p className="text-muted-foreground mb-4">
                {searchQuery || selectedTags.length > 0
                  ? "请尝试调整搜索条件或筛选类型"
                  : "目前没有发布任何活动"
                }
              </p>
              {(searchQuery || selectedTags.length > 0) && (
                <Button asChild variant="outline">
                  <a href="/announcements">清除所有筛选</a>
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