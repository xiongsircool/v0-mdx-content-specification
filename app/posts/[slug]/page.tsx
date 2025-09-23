import { notFound } from "next/navigation"
import { allPosts } from "@/lib/contentlayer-mock"
import { MDXContent } from "@/components/mdx-content"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { Card, CardContent } from "@/components/ui/card"
import { CalendarDays, Clock, User } from "lucide-react"
import { format } from "date-fns"
import { zhCN } from "date-fns/locale"
import Link from "next/link"
import Image from "next/image"
import { GiscusComments } from "@/components/giscus-comments"

interface PostPageProps {
  params: Promise<{
    slug: string
  }>
}

export async function generateStaticParams() {
  return allPosts.map((post) => ({
    slug: post.slug,
  }))
}

export default async function PostPage({ params }: PostPageProps) {
  const { slug } = await params
  const post = allPosts.find((post) => post.slug === slug)

  if (!post) {
    notFound()
  }

  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto px-4 py-8 max-w-4xl">
        {/* Back Navigation */}
        <div className="mb-8">
          <Link href="/posts">
            <Button variant="ghost" className="mb-4">
              ← 返回技术推文
            </Button>
          </Link>
        </div>

        {/* Cover Image */}
        {post.coverImage && (
          <div className="mb-8 overflow-hidden rounded-lg">
            <Image
              src={post.coverImage || "/placeholder.svg"}
              alt={post.title}
              width={800}
              height={400}
              className="w-full object-cover"
              unoptimized
            />
          </div>
        )}

        {/* Header */}
        <header className="mb-8">
          <div className="flex flex-wrap gap-2 mb-4">
            {post.tags.map((tag) => (
              <Badge key={tag} variant="secondary">
                {tag}
              </Badge>
            ))}
          </div>

          <h1 className="text-4xl font-bold text-balance mb-6">{post.title}</h1>

          <p className="text-xl text-muted-foreground text-pretty mb-6">{post.excerpt}</p>

          <div className="flex flex-wrap items-center gap-6 text-sm text-muted-foreground">
            <div className="flex items-center gap-2">
              <User className="h-4 w-4" />
              <span>{post.authors.join(", ")}</span>
            </div>
            <div className="flex items-center gap-2">
              <CalendarDays className="h-4 w-4" />
              <span>发布于 {format(new Date(post.publishedAt), "yyyy年MM月dd日", { locale: zhCN })}</span>
            </div>
            {post.updatedAt && (
              <div className="flex items-center gap-2">
                <Clock className="h-4 w-4" />
                <span>更新于 {format(new Date(post.updatedAt), "yyyy年MM月dd日", { locale: zhCN })}</span>
              </div>
            )}
          </div>
        </header>

        {/* Content */}
        <Card>
          <CardContent className="p-8">
            <MDXContent code={post.body.code} />
          </CardContent>
        </Card>

        {/* Author Info */}
        <Card className="mt-8">
          <CardContent className="p-6">
            <h3 className="text-lg font-semibold mb-3">作者信息</h3>
            <div className="flex flex-wrap gap-4">
              {post.authors.map((author) => (
                <div key={author} className="flex items-center gap-3 p-3 bg-muted rounded-lg">
                  <div className="h-10 w-10 rounded-full bg-primary/10 flex items-center justify-center">
                    <User className="h-5 w-5 text-primary" />
                  </div>
                  <div>
                    <p className="font-medium">{author}</p>
                    <p className="text-sm text-muted-foreground">SBC 成员</p>
                  </div>
                </div>
              ))}
            </div>
          </CardContent>
        </Card>

        {/* Giscus comment system */}
        <Card className="mt-8">
          <CardContent className="p-6">
            <GiscusComments
              repo="xiongsircool/sbc-website"
              repoId="R_kgDONJQqVw"
              category="General"
              categoryId="DIC_kwDONJQqV84CkQHZ"
              mapping="pathname"
              reactionsEnabled={true}
              lang="zh-CN"
            />
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
