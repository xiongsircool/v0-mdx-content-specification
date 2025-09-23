import { notFound } from "next/navigation"
import { allAnnouncements } from "@/lib/contentlayer-mock"
import { MDXContent } from "@/components/mdx-content"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader } from "@/components/ui/card"
import { CalendarDays, Clock, Download, AlertTriangle } from "lucide-react"
import { format } from "date-fns"
import { zhCN } from "date-fns/locale"
import Link from "next/link"

interface AnnouncementPageProps {
  params: Promise<{
    slug: string
  }>
}

export async function generateStaticParams() {
  return allAnnouncements.map((announcement) => ({
    slug: announcement.slug,
  }))
}

export default async function AnnouncementPage({ params }: AnnouncementPageProps) {
  const { slug } = await params
  const announcement = allAnnouncements.find((announcement) => announcement.slug === slug)

  if (!announcement) {
    notFound()
  }

  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto px-4 py-8 max-w-4xl">
        {/* Back Navigation */}
        <div className="mb-8">
          <Link href="/announcements">
            <Button variant="ghost" className="mb-4">
              ← 返回公告列表
            </Button>
          </Link>
        </div>

        {/* Header */}
        <header className="mb-8">
          <div className="flex flex-wrap gap-2 mb-4">
            {announcement.tags.map((tag) => (
              <Badge key={tag} variant="secondary">
                {tag}
              </Badge>
            ))}
            {announcement.priority > 0 && <Badge variant="destructive">置顶</Badge>}
            {announcement.isExpired && (
              <Badge variant="outline" className="text-muted-foreground">
                已过期
              </Badge>
            )}
          </div>

          <h1 className="text-4xl font-bold text-balance mb-6">{announcement.title}</h1>

          <div className="flex flex-wrap items-center gap-6 text-sm text-muted-foreground">
            <div className="flex items-center gap-2">
              <CalendarDays className="h-4 w-4" />
              <span>发布于 {format(new Date(announcement.publishedAt), "yyyy年MM月dd日", { locale: zhCN })}</span>
            </div>
            {announcement.expiresAt && (
              <div className="flex items-center gap-2">
                <Clock className="h-4 w-4" />
                <span>
                  {announcement.isExpired ? "已于" : "截止至"}
                  {format(new Date(announcement.expiresAt), "yyyy年MM月dd日", { locale: zhCN })}
                </span>
              </div>
            )}
          </div>

          {announcement.isExpired && (
            <div className="mt-4 p-4 bg-yellow-50 dark:bg-yellow-950 border border-yellow-200 dark:border-yellow-800 rounded-lg">
              <div className="flex items-center gap-2 text-yellow-800 dark:text-yellow-200">
                <AlertTriangle className="h-5 w-5" />
                <span className="font-medium">此公告已过期</span>
              </div>
            </div>
          )}
        </header>

        {/* Content */}
        <Card className="mb-8">
          <CardContent className="p-8">
            <MDXContent code={announcement.body.code} />
          </CardContent>
        </Card>

        {/* Attachments */}
        {announcement.attachments && announcement.attachments.length > 0 && (
          <Card>
            <CardHeader>
              <h2 className="text-xl font-semibold">相关附件</h2>
            </CardHeader>
            <CardContent>
              <div className="space-y-3">
                {announcement.attachments.map((attachment, index) => (
                  <div key={index} className="flex items-center justify-between p-3 border rounded-lg">
                    <div className="flex items-center gap-3">
                      <Download className="h-5 w-5 text-muted-foreground" />
                      <span className="font-medium">{attachment.name}</span>
                    </div>
                    <Button asChild size="sm">
                      <a href={attachment.url} download>
                        下载
                      </a>
                    </Button>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
        )}
      </div>
    </div>
  )
}
