import { notFound } from "next/navigation"
import { getFileBySlug } from "@/lib/server-markdown-loader"
import { MarkdownContent } from "@/components/markdown-content"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader } from "@/components/ui/card"
import { CalendarDays, Clock, ExternalLink, BookOpen, Users, Target } from "lucide-react"
import { format } from "date-fns"
import { zhCN } from "date-fns/locale"
import Link from "next/link"
import { GiscusComments } from "@/components/giscus-comments"

interface LearnResourcePageProps {
  params: Promise<{
    slug: string
  }>
}

export async function generateStaticParams() {
  // Static generation - return empty array for dynamic rendering
  return []
}

export default async function LearnResourcePage({ params }: LearnResourcePageProps) {
  const { slug } = await params
  const resourceFile = getFileBySlug('learn/resources', slug)

  if (!resourceFile) {
    notFound()
  }

  const resource = resourceFile.metadata

  const levelColors = {
    入门: "bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200",
    进阶: "bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200",
    实战: "bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200",
  }

  const typeIcons = {
    文章: BookOpen,
    视频: ExternalLink,
    项目: Target,
    数据集: Users,
    幻灯片: BookOpen,
  }

  const TypeIcon = typeIcons[resource.type] || BookOpen

  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto px-4 py-8 max-w-4xl">
        {/* Back Navigation */}
        <div className="mb-8">
          <Link href="/learn">
            <Button variant="ghost" className="mb-4">
              ← 返回学习资源
            </Button>
          </Link>
        </div>

        {/* Header */}
        <header className="mb-8">
          <div className="flex flex-wrap items-center gap-3 mb-4">
            <div className="flex items-center gap-2">
              <TypeIcon className="h-5 w-5 text-primary" />
              <Badge variant="outline">{resource.type}</Badge>
            </div>
            <Badge className={levelColors[resource.level]}>{resource.level}</Badge>
            <Badge variant="secondary">{resource.track}</Badge>
          </div>

          <h1 className="text-4xl font-bold text-balance mb-6">{resource.title}</h1>

          <div className="flex flex-wrap items-center gap-6 text-sm text-muted-foreground mb-6">
            <div className="flex items-center gap-2">
              <CalendarDays className="h-4 w-4" />
              <span>更新于 {format(new Date(resource.updatedAt), "yyyy年MM月dd日", { locale: zhCN })}</span>
            </div>
            {resource.duration && (
              <div className="flex items-center gap-2">
                <Clock className="h-4 w-4" />
                <span>预计时长：{resource.duration}</span>
              </div>
            )}
          </div>

          {/* Prerequisites */}
          {resource.prerequisites && resource.prerequisites.length > 0 && (
            <Card className="mb-6">
              <CardHeader className="pb-3">
                <h3 className="text-lg font-semibold">前置要求</h3>
              </CardHeader>
              <CardContent>
                <div className="flex flex-wrap gap-2">
                  {resource.prerequisites.map((prereq) => (
                    <Badge key={prereq} variant="outline">
                      {prereq}
                    </Badge>
                  ))}
                </div>
              </CardContent>
            </Card>
          )}

          {/* External Link */}
          {resource.externalUrl && (
            <Card className="mb-6 border-primary/20 bg-primary/5">
              <CardContent className="p-6">
                <div className="flex items-center justify-between">
                  <div>
                    <h3 className="text-lg font-semibold mb-2">外部资源</h3>
                    <p className="text-muted-foreground">此资源托管在外部平台，点击下方按钮访问。</p>
                  </div>
                  <Button asChild>
                    <a href={resource.externalUrl} target="_blank" rel="noopener noreferrer">
                      <ExternalLink className="h-4 w-4 mr-2" />
                      访问资源
                    </a>
                  </Button>
                </div>
              </CardContent>
            </Card>
          )}
        </header>

        {/* Content */}
        <Card className="mb-8">
          <CardContent className="p-8">
            <div className="max-w-none">
              <MarkdownContent content={resourceFile.htmlContent} />
            </div>
          </CardContent>
        </Card>

        {/* Learning Path Info */}
        <Card className="mb-8">
          <CardContent className="p-6">
            <h3 className="text-lg font-semibold mb-3">学习路径</h3>
            <div className="flex items-center gap-3 p-3 bg-muted rounded-lg">
              <div className="h-10 w-10 rounded-full bg-primary/10 flex items-center justify-center">
                <Target className="h-5 w-5 text-primary" />
              </div>
              <div>
                <p className="font-medium">{resource.track}</p>
                <p className="text-sm text-muted-foreground">
                  难度级别：{resource.level} | 资源类型：{resource.type}
                </p>
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Giscus comment system */}
        <Card>
          <CardContent className="p-6">
            <GiscusComments
              repo="xiongsircool/sbcshanghai"
              repoId="R_kgDOP8J5uA"
              category="Announcements"
              categoryId="DIC_kwDOP8J5uM4CwPdy"
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
