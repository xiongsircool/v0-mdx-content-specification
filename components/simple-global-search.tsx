"use client"

import { useState, useEffect, useRef } from "react"
import { useRouter } from "next/navigation"
import { Search, X, FileText } from "lucide-react"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent } from "@/components/ui/card"
import type { PostMetadata } from "@/lib/simple-content-loader"

interface SimpleSearchResult {
  id: string
  title: string
  excerpt: string
  url: string
  type: 'post'
  tags: string[]
  publishedAt: string
}

export function SimpleGlobalSearch() {
  const [open, setOpen] = useState(false)
  const [query, setQuery] = useState("")
  const [results, setResults] = useState<SimpleSearchResult[]>([])
  const [loading, setLoading] = useState(false)
  const router = useRouter()
  const searchRef = useRef<HTMLDivElement>(null)

  // Simple search function using posts data
  const performSearch = async (searchQuery: string) => {
    if (!searchQuery.trim()) {
      setResults([])
      return
    }

    setLoading(true)
    try {
      const response = await fetch('/api/search?q=' + encodeURIComponent(searchQuery))
      if (response.ok) {
        const data = await response.json()
        setResults(data.results || [])
      }
    } catch (error) {
      console.error('Search error:', error)
      setResults([])
    } finally {
      setLoading(false)
    }
  }

  useEffect(() => {
    const delayedSearch = setTimeout(() => {
      if (query.length > 0) {
        performSearch(query)
      } else {
        setResults([])
      }
    }, 300)

    return () => clearTimeout(delayedSearch)
  }, [query])

  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (searchRef.current && !searchRef.current.contains(event.target as Node)) {
        setOpen(false)
      }
    }

    document.addEventListener('mousedown', handleClickOutside)
    return () => document.removeEventListener('mousedown', handleClickOutside)
  }, [])

  const handleResultClick = (url: string) => {
    setOpen(false)
    setQuery("")
    setResults([])
    router.push(url)
  }

  return (
    <div ref={searchRef} className="relative">
      <Button
        variant="ghost"
        size="sm"
        className="relative h-9 w-9 p-0 md:h-10 md:w-auto md:px-3"
        onClick={() => setOpen(!open)}
      >
        <Search className="h-4 w-4" />
        <span className="hidden md:inline-block md:ml-2">搜索</span>
        <span className="sr-only">搜索</span>
      </Button>

      {open && (
        <div className="absolute top-full mt-2 w-full md:w-96 right-0 z-50">
          <Card className="shadow-lg">
            <CardContent className="p-0">
              <div className="flex items-center border-b px-3">
                <Search className="h-4 w-4 text-muted-foreground" />
                <Input
                  placeholder="搜索文章标题、内容或作者..."
                  value={query}
                  onChange={(e) => setQuery(e.target.value)}
                  className="border-0 bg-transparent focus-visible:ring-0 focus-visible:ring-offset-0"
                  autoFocus
                />
                {query && (
                  <Button
                    variant="ghost"
                    size="sm"
                    className="h-auto p-1"
                    onClick={() => setQuery("")}
                  >
                    <X className="h-4 w-4" />
                  </Button>
                )}
              </div>

              <div className="max-h-96 overflow-y-auto">
                {loading && (
                  <div className="p-4 text-center">
                    <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-primary mx-auto"></div>
                    <p className="text-sm text-muted-foreground mt-2">搜索中...</p>
                  </div>
                )}

                {!loading && query.length > 0 && results.length === 0 && (
                  <div className="p-4 text-center">
                    <FileText className="h-12 w-12 mx-auto mb-4 text-muted-foreground opacity-50" />
                    <p className="text-sm text-muted-foreground">没有找到相关内容</p>
                  </div>
                )}

                {results.length > 0 && (
                  <div className="p-2">
                    {results.map((result) => (
                      <div
                        key={result.id}
                        className="flex flex-col gap-2 rounded-lg p-3 hover:bg-muted cursor-pointer transition-colors"
                        onClick={() => handleResultClick(result.url)}
                      >
                        <div className="flex items-start justify-between gap-2">
                          <div className="flex-1">
                            <h3 className="text-sm font-medium text-balance line-clamp-1">
                              {result.title}
                            </h3>
                            <p className="text-xs text-muted-foreground line-clamp-2 mt-1">
                              {result.excerpt}
                            </p>
                          </div>
                          <Badge variant="secondary" className="text-xs">
                            文章
                          </Badge>
                        </div>
                        <div className="flex flex-wrap gap-1">
                          {result.tags.slice(0, 2).map((tag) => (
                            <Badge key={tag} variant="outline" className="text-xs">
                              {tag}
                            </Badge>
                          ))}
                        </div>
                      </div>
                    ))}
                  </div>
                )}

                {query.length === 0 && (
                  <div className="p-4 text-center">
                    <p className="text-sm text-muted-foreground">输入关键词开始搜索</p>
                  </div>
                )}
              </div>
            </CardContent>
          </Card>
        </div>
      )}
    </div>
  )
}