"use client"

import { useState, useEffect } from "react"
import { PostCard } from "@/components/content-card"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Badge } from "@/components/ui/badge"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Dialog, DialogContent, DialogHeader, DialogTitle, DialogTrigger } from "@/components/ui/dialog"
import { Textarea } from "@/components/ui/textarea"
import { Plus, Edit, Trash2, Search, Filter, BookOpen } from "lucide-react"
import type { Post } from "contentlayer/generated"
import { postManager } from "@/lib/post-manager"

export default function AdminPostsPage() {
  const [posts, setPosts] = useState<any[]>([])
  const [loading, setLoading] = useState(true)
  const [searchQuery, setSearchQuery] = useState("")
  const [isCreateDialogOpen, setIsCreateDialogOpen] = useState(false)
  const [editingPost, setEditingPost] = useState<any>(null)

  // New post form state
  const [newPost, setNewPost] = useState({
    title: "",
    excerpt: "",
    content: "",
    tags: "",
    authors: "",
    coverImage: ""
  })

  // Load posts
  useEffect(() => {
    loadPosts()
  }, [])

  const loadPosts = async () => {
    try {
      const response = await fetch('/api/posts')
      const data = await response.json()
      if (data.success) {
        setPosts(data.posts)
      }
    } catch (error) {
      console.error('Failed to load posts:', error)
    } finally {
      setLoading(false)
    }
  }

  const handleCreatePost = async () => {
    try {
      const postData = {
        title: newPost.title,
        excerpt: newPost.excerpt,
        content: newPost.content,
        tags: newPost.tags.split(',').map((tag: string) => tag.trim()),
        authors: newPost.authors.split(',').map((author: string) => author.trim()),
        coverImage: newPost.coverImage || undefined
      }

      const response = await fetch('/api/posts', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(postData)
      })

      const data = await response.json()
      if (data.success) {
        setIsCreateDialogOpen(false)
        setNewPost({
          title: "",
          excerpt: "",
          content: "",
          tags: "",
          authors: "",
          coverImage: ""
        })
        loadPosts()
      } else {
        alert('创建推文失败: ' + data.error)
      }
    } catch (error) {
      console.error('Failed to create post:', error)
      alert('创建推文失败')
    }
  }

  const handleDeletePost = async (slug: string) => {
    if (!confirm('确定要删除这篇推文吗？')) {
      return
    }

    try {
      const response = await fetch(`/api/posts/${slug}`, {
        method: 'DELETE'
      })

      const data = await response.json()
      if (data.success) {
        loadPosts()
      } else {
        alert('删除推文失败: ' + data.error)
      }
    } catch (error) {
      console.error('Failed to delete post:', error)
      alert('删除推文失败')
    }
  }

  const filteredPosts = posts.filter(post =>
    post.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
    post.excerpt.toLowerCase().includes(searchQuery.toLowerCase()) ||
    post.tags.some((tag: string) => tag.toLowerCase().includes(searchQuery.toLowerCase()))
  )

  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto px-4 py-8">
        {/* Header */}
        <div className="flex items-center justify-between mb-8">
          <div>
            <h1 className="text-3xl font-bold mb-2">推文管理</h1>
            <p className="text-muted-foreground">管理技术推文内容</p>
          </div>
          <Dialog open={isCreateDialogOpen} onOpenChange={setIsCreateDialogOpen}>
            <DialogTrigger asChild>
              <Button>
                <Plus className="h-4 w-4 mr-2" />
                新建推文
              </Button>
            </DialogTrigger>
            <DialogContent className="max-w-2xl max-h-[80vh] overflow-y-auto">
              <DialogHeader>
                <DialogTitle>创建新推文</DialogTitle>
              </DialogHeader>
              <div className="space-y-4">
                <div>
                  <label className="text-sm font-medium mb-2 block">标题</label>
                  <Input
                    value={newPost.title}
                    onChange={(e) => setNewPost({ ...newPost, title: e.target.value })}
                    placeholder="输入推文标题"
                  />
                </div>
                <div>
                  <label className="text-sm font-medium mb-2 block">摘要</label>
                  <Textarea
                    value={newPost.excerpt}
                    onChange={(e) => setNewPost({ ...newPost, excerpt: e.target.value })}
                    placeholder="输入推文摘要"
                    rows={3}
                  />
                </div>
                <div>
                  <label className="text-sm font-medium mb-2 block">标签（逗号分隔）</label>
                  <Input
                    value={newPost.tags}
                    onChange={(e) => setNewPost({ ...newPost, tags: e.target.value })}
                    placeholder="例如：Python, 生物信息学, 教程"
                  />
                </div>
                <div>
                  <label className="text-sm font-medium mb-2 block">作者（逗号分隔）</label>
                  <Input
                    value={newPost.authors}
                    onChange={(e) => setNewPost({ ...newPost, authors: e.target.value })}
                    placeholder="例如：张三, 李四"
                  />
                </div>
                <div>
                  <label className="text-sm font-medium mb-2 block">封面图片URL（可选）</label>
                  <Input
                    value={newPost.coverImage}
                    onChange={(e) => setNewPost({ ...newPost, coverImage: e.target.value })}
                    placeholder="https://example.com/image.jpg"
                  />
                </div>
                <div>
                  <label className="text-sm font-medium mb-2 block">内容</label>
                  <Textarea
                    value={newPost.content}
                    onChange={(e) => setNewPost({ ...newPost, content: e.target.value })}
                    placeholder="输入推文内容（支持Markdown）"
                    rows={10}
                  />
                </div>
                <div className="flex gap-2">
                  <Button onClick={handleCreatePost} className="flex-1">
                    创建推文
                  </Button>
                  <Button variant="outline" onClick={() => setIsCreateDialogOpen(false)}>
                    取消
                  </Button>
                </div>
              </div>
            </DialogContent>
          </Dialog>
        </div>

        {/* Search */}
        <Card className="mb-6">
          <CardContent className="p-4">
            <div className="relative">
              <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-muted-foreground" />
              <Input
                placeholder="搜索推文标题、摘要或标签..."
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                className="pl-10"
              />
            </div>
          </CardContent>
        </Card>

        {/* Posts List */}
        {loading ? (
          <Card>
            <CardContent className="p-12 text-center">
              <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto mb-4"></div>
              <h3 className="text-xl font-semibold mb-2">加载中...</h3>
              <p className="text-muted-foreground">正在获取推文列表</p>
            </CardContent>
          </Card>
        ) : (
          <div className="space-y-4">
            {filteredPosts.map((post) => (
              <Card key={post.slug}>
                <CardContent className="p-6">
                  <div className="flex items-start justify-between">
                    <div className="flex-1">
                      <h3 className="text-lg font-semibold mb-2">{post.title}</h3>
                      <p className="text-muted-foreground mb-3">{post.excerpt}</p>
                      <div className="flex flex-wrap gap-2 mb-3">
                        {post.tags.map((tag: string) => (
                          <Badge key={tag} variant="secondary" className="text-xs">
                            {tag}
                          </Badge>
                        ))}
                      </div>
                      <div className="text-sm text-muted-foreground">
                        <span>作者: {post.authors.join(', ')}</span>
                        <span className="mx-2">•</span>
                        <span>发布: {new Date(post.publishedAt).toLocaleDateString('zh-CN')}</span>
                        {post.updatedAt && (
                          <>
                            <span className="mx-2">•</span>
                            <span>更新: {new Date(post.updatedAt).toLocaleDateString('zh-CN')}</span>
                          </>
                        )}
                      </div>
                    </div>
                    <div className="flex gap-2 ml-4">
                      <Button variant="outline" size="sm">
                        <Edit className="h-4 w-4 mr-1" />
                        编辑
                      </Button>
                      <Button
                        variant="destructive"
                        size="sm"
                        onClick={() => handleDeletePost(post.slug)}
                      >
                        <Trash2 className="h-4 w-4 mr-1" />
                        删除
                      </Button>
                    </div>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        )}
      </div>
    </div>
  )
}