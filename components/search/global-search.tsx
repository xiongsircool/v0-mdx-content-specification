"use client"

import { useState, useEffect, useRef } from "react"
import { useRouter } from "next/navigation"
import { Search, X, FileText, Megaphone, GraduationCap, Calendar } from "lucide-react"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent } from "@/components/ui/card"
import { allPosts, allAnnouncements, allLearnResources, allMeetings } from "contentlayer/generated"

interface SearchResult {
  id: string
  title: string
  excerpt: string
  type: "post" | "announcement" | "resource" | "meeting"
  url: string
  tags?: string[]
  date: string
}

const typeConfig = {
  post: { label: "技术推文", icon: FileText, color: "bg-blue-100 text-blue-800" },
  announcement: { label: "公告", icon: Megaphone, color: "bg-red-100 text-red-800" },
  resource: { label: "学习资源", icon: GraduationCap, color: "bg-green-100 text-green-800" },
  meeting: { label: "会议", icon: Calendar, color: "bg-purple-100 text-purple-800" },
}

export function GlobalSearch() {
  const [query, setQuery] = useState("")
  const [results, setResults] = useState<SearchResult[]>([])
  const [isOpen, setIsOpen] = useState(false)
  const [selectedType, setSelectedType] = useState<string>("all")
  const searchRef = useRef<HTMLDivElement>(null)
  const router = useRouter()

  // 搜索所有内容
  const searchContent = (searchQuery: string) => {
    if (!searchQuery.trim()) {
      setResults([])
      return
    }

    const query = searchQuery.toLowerCase()
    const allResults: SearchResult[] = []

    // 搜索技术推文
    allPosts.forEach((post) => {
      if (
        post.title?.toLowerCase().includes(query) ||
        post.excerpt?.toLowerCase().includes(query) ||
        post.tags?.some((tag) => tag?.toLowerCase().includes(query))
      ) {
        allResults.push({
          id: post.slug,
          title: post.title,
          excerpt: post.excerpt,
          type: "post",
          url: `/posts/${post.slug}`,
          tags: post.tags,
          date: post.publishedAt,
        })
      }
    })

    // 搜索公告
    allAnnouncements.forEach((announcement) => {
      if (
        announcement.title?.toLowerCase().includes(query) ||
        announcement.summary?.toLowerCase().includes(query) ||
        announcement.tags?.some((tag) => tag?.toLowerCase().includes(query))
      ) {
        allResults.push({
          id: announcement.slug,
          title: announcement.title,
          excerpt: announcement.summary,
          type: "announcement",
          url: `/announcements/${announcement.slug}`,
          tags: announcement.tags,
          date: announcement.publishedAt,
        })
      }
    })

    // 搜索学习资源
    allLearnResources.forEach((resource) => {
      const description = resource.body?.raw || ""
      if (
        resource.title?.toLowerCase().includes(query) ||
        description.toLowerCase().includes(query) ||
        resource.track?.toLowerCase().includes(query) ||
        resource.level?.toLowerCase().includes(query)
      ) {
        allResults.push({
          id: resource.slug,
          title: resource.title,
          excerpt: `${resource.level} | ${resource.track} | ${resource.duration || ""}`,
          type: "resource",
          url: `/learn/resources/${resource.slug}`,
          tags: [resource.level, resource.track].filter(Boolean),
          date: resource.updatedAt,
        })
      }
    })

    // 搜索会议
    allMeetings.forEach((meeting) => {
      if (
        meeting.title?.toLowerCase().includes(query) ||
        meeting.summary?.toLowerCase().includes(query) ||
        meeting.speaker?.toLowerCase().includes(query) ||
        meeting.tags?.some((tag) => tag?.toLowerCase().includes(query))
      ) {
        allResults.push({
          id: meeting.slug,
          title: meeting.title,
          excerpt: meeting.summary,
          type: "meeting",
          url: `/meetings`,
          tags: meeting.tags,
          date: meeting.date,
        })
      }
    })

    // 按类型筛选
    const filteredResults =
      selectedType === "all" ? allResults : allResults.filter((result) => result.type === selectedType)

    // 按相关性排序（标题匹配优先）
    filteredResults.sort((a, b) => {
      const aTitle = a.title?.toLowerCase().includes(query) ? 1 : 0
      const bTitle = b.title?.toLowerCase().includes(query) ? 1 : 0
      return bTitle - aTitle
    })

    setResults(filteredResults.slice(0, 8)) // 限制显示8个结果
  }

  useEffect(() => {
    searchContent(query)
  }, [query, selectedType])

  // 点击外部关闭搜索
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (searchRef.current && !searchRef.current.contains(event.target as Node)) {
        setIsOpen(false)
      }
    }

    document.addEventListener("mousedown", handleClickOutside)
    return () => document.removeEventListener("mousedown", handleClickOutside)
  }, [])

  const handleResultClick = (result: SearchResult) => {
    router.push(result.url)
    setIsOpen(false)
    setQuery("")
  }

  const handleViewAllResults = () => {
    router.push(`/search?q=${encodeURIComponent(query)}&type=${selectedType}`)
    setIsOpen(false)
  }

  return (
    <div ref={searchRef} className="relative">
      <div className="relative">
        <Search className="absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-muted-foreground" />
        <Input
          type="text"
          placeholder="搜索内容..."
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          onFocus={() => setIsOpen(true)}
          className="pl-10 pr-10 w-64"
        />
        {query && (
          <Button
            variant="ghost"
            size="sm"
            onClick={() => {
              setQuery("")
              setResults([])
            }}
            className="absolute right-1 top-1/2 h-6 w-6 -translate-y-1/2 p-0"
          >
            <X className="h-3 w-3" />
          </Button>
        )}
      </div>

      {isOpen && (query || results.length > 0) && (
        <Card className="absolute top-full mt-2 w-96 max-h-96 overflow-hidden shadow-lg z-50">
          <CardContent className="p-0">
            {/* 类型筛选 */}
            <div className="p-3 border-b bg-muted/50">
              <div className="flex gap-2 flex-wrap">
                <Button
                  variant={selectedType === "all" ? "default" : "outline"}
                  size="sm"
                  onClick={() => setSelectedType("all")}
                >
                  全部
                </Button>
                {Object.entries(typeConfig).map(([type, config]) => (
                  <Button
                    key={type}
                    variant={selectedType === type ? "default" : "outline"}
                    size="sm"
                    onClick={() => setSelectedType(type)}
                  >
                    <config.icon className="h-3 w-3 mr-1" />
                    {config.label}
                  </Button>
                ))}
              </div>
            </div>

            {/* 搜索结果 */}
            <div className="max-h-64 overflow-y-auto">
              {results.length > 0 ? (
                <>
                  {results.map((result) => {
                    const config = typeConfig[result.type]
                    const Icon = config.icon
                    return (
                      <div
                        key={`${result.type}-${result.id}`}
                        onClick={() => handleResultClick(result)}
                        className="p-3 hover:bg-muted/50 cursor-pointer border-b last:border-b-0"
                      >
                        <div className="flex items-start gap-3">
                          <Icon className="h-4 w-4 mt-1 text-muted-foreground flex-shrink-0" />
                          <div className="flex-1 min-w-0">
                            <div className="flex items-center gap-2 mb-1">
                              <h4 className="font-medium text-sm truncate">{result.title}</h4>
                              <Badge variant="secondary" className={`text-xs ${config.color}`}>
                                {config.label}
                              </Badge>
                            </div>
                            <p className="text-xs text-muted-foreground line-clamp-2">{result.excerpt}</p>
                            {result.tags && (
                              <div className="flex gap-1 mt-2">
                                {result.tags.slice(0, 3).map((tag) => (
                                  <Badge key={tag} variant="outline" className="text-xs">
                                    {tag}
                                  </Badge>
                                ))}
                              </div>
                            )}
                          </div>
                        </div>
                      </div>
                    )
                  })}
                  {results.length >= 8 && (
                    <div className="p-3 border-t bg-muted/30">
                      <Button variant="ghost" size="sm" onClick={handleViewAllResults} className="w-full">
                        查看全部结果
                      </Button>
                    </div>
                  )}
                </>
              ) : query ? (
                <div className="p-6 text-center text-muted-foreground">
                  <Search className="h-8 w-8 mx-auto mb-2 opacity-50" />
                  <p>未找到相关内容</p>
                </div>
              ) : null}
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  )
}
