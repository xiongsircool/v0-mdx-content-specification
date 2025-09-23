"use client"

import { useState, useEffect, Suspense } from "react"
import { useSearchParams } from "next/navigation"
import { Search, FileText, Megaphone, GraduationCap, Calendar, Filter } from "lucide-react"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { allPosts, allAnnouncements, allLearnResources, allMeetings } from "@/lib/contentlayer-mock"
import Link from "next/link"

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

function SearchPageContent() {
  const searchParams = useSearchParams()
  const [query, setQuery] = useState(searchParams.get("q") || "")
  const [selectedType, setSelectedType] = useState(searchParams.get("type") || "all")
  const [sortBy, setSortBy] = useState("relevance")
  const [results, setResults] = useState<SearchResult[]>([])

  const searchContent = (searchQuery: string) => {
    if (!searchQuery.trim()) {
      setResults([])
      return
    }

    const query = searchQuery.toLowerCase()
    const allResults: SearchResult[] = []

    // 搜索所有内容类型
    allPosts.forEach((post) => {
      if (
        post.title.toLowerCase().includes(query) ||
        post.excerpt.toLowerCase().includes(query) ||
        post.tags?.some((tag) => tag.toLowerCase().includes(query))
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

    allAnnouncements.forEach((announcement) => {
      if (
        announcement.title.toLowerCase().includes(query) ||
        announcement.excerpt.toLowerCase().includes(query) ||
        announcement.tags?.some((tag) => tag.toLowerCase().includes(query))
      ) {
        allResults.push({
          id: announcement.slug,
          title: announcement.title,
          excerpt: announcement.excerpt,
          type: "announcement",
          url: `/announcements/${announcement.slug}`,
          tags: announcement.tags,
          date: announcement.publishedAt,
        })
      }
    })

    allLearnResources.forEach((resource) => {
      if (
        resource.title.toLowerCase().includes(query) ||
        resource.description.toLowerCase().includes(query) ||
        resource.tags?.some((tag) => tag.toLowerCase().includes(query))
      ) {
        allResults.push({
          id: resource.slug,
          title: resource.title,
          excerpt: resource.description,
          type: "resource",
          url: `/learn/resources/${resource.slug}`,
          tags: resource.tags,
          date: resource.publishedAt,
        })
      }
    })

    allMeetings.forEach((meeting) => {
      if (
        meeting.title.toLowerCase().includes(query) ||
        meeting.description.toLowerCase().includes(query) ||
        meeting.speaker.toLowerCase().includes(query)
      ) {
        allResults.push({
          id: meeting.id,
          title: meeting.title,
          excerpt: meeting.description,
          type: "meeting",
          url: `/meetings`,
          date: meeting.date,
        })
      }
    })

    // 按类型筛选
    const filteredResults =
      selectedType === "all" ? allResults : allResults.filter((result) => result.type === selectedType)

    // 排序
    if (sortBy === "date") {
      filteredResults.sort((a, b) => new Date(b.date).getTime() - new Date(a.date).getTime())
    } else {
      // 按相关性排序（标题匹配优先）
      filteredResults.sort((a, b) => {
        const aTitle = a.title.toLowerCase().includes(query) ? 1 : 0
        const bTitle = b.title.toLowerCase().includes(query) ? 1 : 0
        return bTitle - aTitle
      })
    }

    setResults(filteredResults)
  }

  useEffect(() => {
    searchContent(query)
  }, [query, selectedType, sortBy])

  return (
    <div className="container mx-auto px-4 py-8">
      <div className="max-w-4xl mx-auto">
        {/* 搜索头部 */}
        <div className="mb-8">
          <h1 className="text-3xl font-bold mb-4">搜索结果</h1>

          {/* 搜索框 */}
          <div className="relative mb-6">
            <Search className="absolute left-3 top-1/2 h-5 w-5 -translate-y-1/2 text-muted-foreground" />
            <Input
              type="text"
              placeholder="搜索内容..."
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              className="pl-12 text-lg h-12"
            />
          </div>

          {/* 筛选器 */}
          <div className="flex flex-wrap gap-4 items-center">
            <div className="flex items-center gap-2">
              <Filter className="h-4 w-4" />
              <span className="text-sm font-medium">筛选:</span>
            </div>

            <Select value={selectedType} onValueChange={setSelectedType}>
              <SelectTrigger className="w-32">
                <SelectValue />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="all">全部类型</SelectItem>
                <SelectItem value="post">技术推文</SelectItem>
                <SelectItem value="announcement">公告</SelectItem>
                <SelectItem value="resource">学习资源</SelectItem>
                <SelectItem value="meeting">会议</SelectItem>
              </SelectContent>
            </Select>

            <Select value={sortBy} onValueChange={setSortBy}>
              <SelectTrigger className="w-32">
                <SelectValue />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="relevance">相关性</SelectItem>
                <SelectItem value="date">发布时间</SelectItem>
              </SelectContent>
            </Select>
          </div>
        </div>

        {/* 搜索结果 */}
        {query && (
          <div className="mb-4">
            <p className="text-muted-foreground">
              找到 <span className="font-medium">{results.length}</span> 个结果
              {query && (
                <span>
                  {" "}
                  关于 "<span className="font-medium">{query}</span>"
                </span>
              )}
            </p>
          </div>
        )}

        {results.length > 0 ? (
          <div className="space-y-6">
            {results.map((result) => {
              const config = typeConfig[result.type]
              const Icon = config.icon
              return (
                <Card key={`${result.type}-${result.id}`} className="hover:shadow-md transition-shadow">
                  <CardHeader className="pb-3">
                    <div className="flex items-start justify-between gap-4">
                      <div className="flex-1">
                        <div className="flex items-center gap-2 mb-2">
                          <Icon className="h-4 w-4 text-muted-foreground" />
                          <Badge variant="secondary" className={`text-xs ${config.color}`}>
                            {config.label}
                          </Badge>
                        </div>
                        <CardTitle className="text-xl">
                          <Link href={result.url} className="hover:text-primary transition-colors">
                            {result.title}
                          </Link>
                        </CardTitle>
                      </div>
                    </div>
                  </CardHeader>
                  <CardContent className="pt-0">
                    <p className="text-muted-foreground mb-4 line-clamp-3">{result.excerpt}</p>
                    <div className="flex items-center justify-between">
                      <div className="flex gap-2">
                        {result.tags?.slice(0, 4).map((tag) => (
                          <Badge key={tag} variant="outline" className="text-xs">
                            {tag}
                          </Badge>
                        ))}
                      </div>
                      <span className="text-sm text-muted-foreground">
                        {new Date(result.date).toLocaleDateString("zh-CN")}
                      </span>
                    </div>
                  </CardContent>
                </Card>
              )
            })}
          </div>
        ) : query ? (
          <div className="text-center py-12">
            <Search className="h-16 w-16 mx-auto mb-4 text-muted-foreground opacity-50" />
            <h3 className="text-xl font-medium mb-2">未找到相关内容</h3>
            <p className="text-muted-foreground">尝试使用不同的关键词或调整筛选条件</p>
          </div>
        ) : (
          <div className="text-center py-12">
            <Search className="h-16 w-16 mx-auto mb-4 text-muted-foreground opacity-50" />
            <h3 className="text-xl font-medium mb-2">开始搜索</h3>
            <p className="text-muted-foreground">输入关键词来搜索技术推文、公告、学习资源和会议信息</p>
          </div>
        )}
      </div>
    </div>
  )
}

export default function SearchPage() {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <SearchPageContent />
    </Suspense>
  )
}
