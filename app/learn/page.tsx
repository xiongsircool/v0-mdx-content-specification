"use client"

import { useState, useMemo, useEffect } from "react"
import { ContentCard } from "@/components/content-card"
import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Search, Filter, GraduationCap, BookOpen, Target, Users } from "lucide-react"
import type { ContentItem } from "@/lib/server-markdown-loader"

const resourceTypes = ["文章", "视频", "项目", "数据集", "幻灯片"]
const difficultyLevels = ["入门", "进阶", "实战"]

export default function LearnPage() {
  const [searchQuery, setSearchQuery] = useState("")
  const [selectedTypes, setSelectedTypes] = useState<string[]>([])
  const [selectedLevels, setSelectedLevels] = useState<string[]>([])
  const [selectedTrack, setSelectedTrack] = useState<string>("")
  const [resources, setResources] = useState<ContentItem[]>([])
  const [loading, setLoading] = useState(true)

  // Load resources on component mount
  useEffect(() => {
    const loadData = async () => {
      try {
        const response = await fetch('/api/learn')
        const resourcesData = await response.json()
        setResources(resourcesData)
      } catch (error) {
        console.error("Failed to load learn resources:", error)
      } finally {
        setLoading(false)
      }
    }

    loadData()
  }, [])

  // Get all unique tracks from resources
  const allTracks = useMemo(() => {
    const tracks = new Set<string>()
    resources.forEach((resource) => {
      tracks.add(resource.track)
    })
    return Array.from(tracks).sort()
  }, [resources])

  const filteredResources = useMemo(() => {
    return resources
      .filter((resource) => {
        // Search filter
        const matchesSearch =
          searchQuery === "" ||
          resource.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
          resource.track.toLowerCase().includes(searchQuery.toLowerCase())

        // Type filter
        const matchesType = selectedTypes.length === 0 || selectedTypes.includes(resource.type)

        // Level filter
        const matchesLevel = selectedLevels.length === 0 || selectedLevels.includes(resource.level)

        // Track filter
        const matchesTrack = selectedTrack === "" || resource.track === selectedTrack

        return matchesSearch && matchesType && matchesLevel && matchesTrack
      })
      .sort((a, b) => new Date(b.updatedAt).getTime() - new Date(a.updatedAt).getTime())
  }, [resources, searchQuery, selectedTypes, selectedLevels, selectedTrack])

  const toggleType = (type: string) => {
    setSelectedTypes((prev) => (prev.includes(type) ? prev.filter((t) => t !== type) : [...prev, type]))
  }

  const toggleLevel = (level: string) => {
    setSelectedLevels((prev) => (prev.includes(level) ? prev.filter((l) => l !== level) : [...prev, level]))
  }

  // Group resources by track for overview
  const resourcesByTrack = useMemo(() => {
    const grouped = resources.reduce(
      (acc, resource) => {
        if (!acc[resource.track]) {
          acc[resource.track] = []
        }
        acc[resource.track].push(resource)
        return acc
      },
      {} as Record<string, typeof resources>,
    )
    return grouped
  }, [resources])

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
          <h1 className="text-4xl font-bold mb-4">学习资料</h1>
          <p className="text-xl text-muted-foreground max-w-2xl mx-auto">
            我们整理的学习资料和培训内容，一起学习生物信息学
          </p>
        </div>

        {/* Learning Areas Overview */}
        <div className="mb-12">
          <h2 className="text-2xl font-bold mb-6">学习领域</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {Object.entries(resourcesByTrack).map(([track, resources]) => (
              <Card key={track} className="cursor-pointer transition-all hover:shadow-md">
                <CardHeader className="pb-3">
                  <div className="flex items-center gap-3">
                    <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-primary/10">
                      <Target className="h-5 w-5 text-primary" />
                    </div>
                    <CardTitle className="text-lg">{track}</CardTitle>
                  </div>
                </CardHeader>
                <CardContent>
                  <div className="space-y-3">
                    <div className="flex items-center gap-4 text-sm text-muted-foreground">
                      <div className="flex items-center gap-1">
                        <BookOpen className="h-4 w-4" />
                        <span>{resources.length} 个资料</span>
                      </div>
                      <div className="flex items-center gap-1">
                        <Users className="h-4 w-4" />
                        <span>
                          {resources.filter((r) => r.level === "入门").length} 入门 ·{" "}
                          {resources.filter((r) => r.level === "进阶").length} 进阶 ·{" "}
                          {resources.filter((r) => r.level === "实战").length} 实战
                        </span>
                      </div>
                    </div>
                    <Button
                      variant="outline"
                      size="sm"
                      className="w-full bg-transparent"
                      onClick={() => setSelectedTrack(selectedTrack === track ? "" : track)}
                    >
                      {selectedTrack === track ? "显示全部" : "查看资料"}
                    </Button>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        </div>

        {/* Filters */}
        <Card className="mb-8">
          <CardContent className="p-6">
            <div className="space-y-4">
              {/* Search */}
              <div className="relative">
                <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-muted-foreground" />
                <Input
                  placeholder="搜索学习资源..."
                  value={searchQuery}
                  onChange={(e) => setSearchQuery(e.target.value)}
                  className="pl-10"
                />
              </div>

              {/* Filters */}
              <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                {/* Resource Types */}
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <Filter className="h-4 w-4 text-muted-foreground" />
                    <span className="text-sm font-medium">资料类型</span>
                  </div>
                  <div className="flex flex-wrap gap-2">
                    {resourceTypes.map((type) => (
                      <Badge
                        key={type}
                        variant={selectedTypes.includes(type) ? "default" : "outline"}
                        className="cursor-pointer transition-colors"
                        onClick={() => toggleType(type)}
                      >
                        {type}
                      </Badge>
                    ))}
                  </div>
                </div>

                {/* Difficulty Levels */}
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <Target className="h-4 w-4 text-muted-foreground" />
                    <span className="text-sm font-medium">难度级别</span>
                  </div>
                  <div className="flex flex-wrap gap-2">
                    {difficultyLevels.map((level) => (
                      <Badge
                        key={level}
                        variant={selectedLevels.includes(level) ? "default" : "outline"}
                        className="cursor-pointer transition-colors"
                        onClick={() => toggleLevel(level)}
                      >
                        {level}
                      </Badge>
                    ))}
                  </div>
                </div>

                {/* Learning Tracks */}
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <Target className="h-4 w-4 text-muted-foreground" />
                    <span className="text-sm font-medium">学习领域</span>
                  </div>
                  <div className="flex flex-wrap gap-2">
                    {allTracks.slice(0, 4).map((track) => (
                      <Badge
                        key={track}
                        variant={selectedTrack === track ? "default" : "outline"}
                        className="cursor-pointer transition-colors"
                        onClick={() => setSelectedTrack(selectedTrack === track ? "" : track)}
                      >
                        {track}
                      </Badge>
                    ))}
                  </div>
                </div>
              </div>

              {/* Clear Filters */}
              {(selectedTypes.length > 0 || selectedLevels.length > 0 || selectedTrack) && (
                <div>
                  <Button
                    variant="ghost"
                    size="sm"
                    onClick={() => {
                      setSelectedTypes([])
                      setSelectedLevels([])
                      setSelectedTrack("")
                    }}
                  >
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
            找到 {filteredResources.length} 个学习资料
            {(selectedTypes.length > 0 || selectedLevels.length > 0 || selectedTrack) && (
              <span>
                {" "}
                · 筛选条件:{" "}
                {[...selectedTypes, ...selectedLevels, ...(selectedTrack ? [selectedTrack] : [])].join(", ")}
              </span>
            )}
          </p>
        </div>

        {/* Resources Grid */}
        {filteredResources.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredResources.map((resource) => (
              <ContentCard key={resource.slug} item={resource} type="resource" />
            ))}
          </div>
        ) : (
          <Card>
            <CardContent className="p-12 text-center">
              <BookOpen className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-xl font-semibold mb-2">没有找到相关资料</h3>
              <p className="text-muted-foreground mb-4">
                {searchQuery || selectedTypes.length > 0 || selectedLevels.length > 0 || selectedTrack
                  ? "请尝试调整搜索条件或筛选选项"
                  : "目前没有发布任何学习资料"}
              </p>
              {(searchQuery || selectedTypes.length > 0 || selectedLevels.length > 0 || selectedTrack) && (
                <Button
                  variant="outline"
                  onClick={() => {
                    setSearchQuery("")
                    setSelectedTypes([])
                    setSelectedLevels([])
                    setSelectedTrack("")
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
