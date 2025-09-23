"use client"

import { useState, useMemo } from "react"
import { allMeetings } from "@/lib/contentlayer-mock"
import { MeetingCard } from "@/components/content-card"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent } from "@/components/ui/card"
import { Search, Filter, Calendar, Video, MapPin, Users } from "lucide-react"

const allTags = ["工作坊", "文献讨论", "月度会议", "特邀讲座", "Python", "AI", "单细胞", "多组学"]
const allStatuses = ["即将开始", "进行中", "已结束"]
const allTypes = ["线上", "线下", "混合"]

export default function MeetingsPage() {
  const [searchQuery, setSearchQuery] = useState("")
  const [selectedTags, setSelectedTags] = useState<string[]>([])
  const [selectedStatuses, setSelectedStatuses] = useState<string[]>([])
  const [selectedTypes, setSelectedTypes] = useState<string[]>([])

  const filteredMeetings = useMemo(() => {
    return allMeetings
      .filter((meeting) => {
        // Search filter
        const matchesSearch =
          searchQuery === "" ||
          meeting.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
          meeting.summary.toLowerCase().includes(searchQuery.toLowerCase()) ||
          meeting.speaker?.toLowerCase().includes(searchQuery.toLowerCase())

        // Tag filter
        const matchesTags = selectedTags.length === 0 || selectedTags.some((tag) => meeting.tags.includes(tag))

        // Status filter
        const matchesStatus = selectedStatuses.length === 0 || selectedStatuses.includes(meeting.status)

        // Type filter
        const matchesType = selectedTypes.length === 0 || selectedTypes.includes(meeting.type)

        return matchesSearch && matchesTags && matchesStatus && matchesType
      })
      .sort((a, b) => {
        // Sort by date (upcoming first, then by date)
        if (a.status === "已结束" && b.status !== "已结束") return 1
        if (a.status !== "已结束" && b.status === "已结束") return -1
        return new Date(a.date).getTime() - new Date(b.date).getTime()
      })
  }, [searchQuery, selectedTags, selectedStatuses, selectedTypes])

  const toggleTag = (tag: string) => {
    setSelectedTags((prev) => (prev.includes(tag) ? prev.filter((t) => t !== tag) : [...prev, tag]))
  }

  const toggleStatus = (status: string) => {
    setSelectedStatuses((prev) => (prev.includes(status) ? prev.filter((s) => s !== status) : [...prev, status]))
  }

  const toggleType = (type: string) => {
    setSelectedTypes((prev) => (prev.includes(type) ? prev.filter((t) => t !== type) : [...prev, type]))
  }

  const clearAllFilters = () => {
    setSearchQuery("")
    setSelectedTags([])
    setSelectedStatuses([])
    setSelectedTypes([])
  }

  const hasActiveFilters =
    searchQuery || selectedTags.length > 0 || selectedStatuses.length > 0 || selectedTypes.length > 0

  return (
    <div className="min-h-screen bg-background">
      <Navigation />

      <div className="container mx-auto px-4 py-8">
        {/* Header */}
        <div className="text-center mb-12">
          <div className="flex justify-center mb-6">
            <div className="flex h-12 w-12 items-center justify-center rounded-xl bg-primary/10">
              <Calendar className="h-6 w-6 text-primary" />
            </div>
          </div>
          <h1 className="text-4xl font-bold mb-4">会议中心</h1>
          <p className="text-xl text-muted-foreground max-w-2xl mx-auto">
            参加我们的学术会议、工作坊和讨论会，与同行交流学习
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
                  placeholder="搜索会议标题、内容或主讲人..."
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

              {/* Status Filters */}
              <div className="space-y-3">
                <div className="flex items-center gap-2">
                  <Calendar className="h-4 w-4 text-muted-foreground" />
                  <span className="text-sm font-medium">状态筛选</span>
                </div>
                <div className="flex flex-wrap gap-2">
                  {allStatuses.map((status) => (
                    <Badge
                      key={status}
                      variant={selectedStatuses.includes(status) ? "default" : "outline"}
                      className="cursor-pointer transition-colors"
                      onClick={() => toggleStatus(status)}
                    >
                      {status}
                    </Badge>
                  ))}
                </div>
              </div>

              {/* Type Filters */}
              <div className="space-y-3">
                <div className="flex items-center gap-2">
                  <Video className="h-4 w-4 text-muted-foreground" />
                  <span className="text-sm font-medium">类型筛选</span>
                </div>
                <div className="flex flex-wrap gap-2">
                  {allTypes.map((type) => (
                    <Badge
                      key={type}
                      variant={selectedTypes.includes(type) ? "default" : "outline"}
                      className="cursor-pointer transition-colors"
                      onClick={() => toggleType(type)}
                    >
                      <span className="flex items-center gap-1">
                        {type === "线上" && <Video className="h-3 w-3" />}
                        {type === "线下" && <MapPin className="h-3 w-3" />}
                        {type === "混合" && <Users className="h-3 w-3" />}
                        {type}
                      </span>
                    </Badge>
                  ))}
                </div>
              </div>

              {/* Clear Filters */}
              {hasActiveFilters && (
                <div className="flex items-center gap-2">
                  <Button variant="ghost" size="sm" onClick={clearAllFilters}>
                    清除所有筛选
                  </Button>
                </div>
              )}
            </div>
          </CardContent>
        </Card>

        {/* Results */}
        <div className="mb-6">
          <p className="text-sm text-muted-foreground">
            找到 {filteredMeetings.length} 个会议
            {selectedTags.length > 0 && <span> · 标签: {selectedTags.join(", ")}</span>}
            {selectedStatuses.length > 0 && <span> · 状态: {selectedStatuses.join(", ")}</span>}
            {selectedTypes.length > 0 && <span> · 类型: {selectedTypes.join(", ")}</span>}
          </p>
        </div>

        {/* Meetings Grid */}
        {filteredMeetings.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredMeetings.map((meeting) => (
              <MeetingCard key={meeting.slug} meeting={meeting} />
            ))}
          </div>
        ) : (
          <Card>
            <CardContent className="p-12 text-center">
              <Calendar className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-xl font-semibold mb-2">没有找到相关会议</h3>
              <p className="text-muted-foreground mb-4">
                {hasActiveFilters ? "请尝试调整搜索条件或筛选标签" : "目前没有安排任何会议"}
              </p>
              {hasActiveFilters && (
                <Button variant="outline" onClick={clearAllFilters}>
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
